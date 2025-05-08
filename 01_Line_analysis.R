####################################################################################################
#### Line analysis
#### Author: Ben Gannon (benjamin.gannon@usda.gov)
#### Date Created: 07/22/2021
#### Last Modified: 05/30/2024
####################################################################################################
# Summary: ingests polylines or polygons, converts to polyline if necessary, generates sample points
# at fixed distance intervals along the lines, extracts raster data for sample point buffers, and 
# outputs point, polyline, and tabular results.
# Input data is controlled in the companion spreadsheet.
# - Analysis lines: can be polyline or polygon
# - Rasters: including options for extraction function and several data transformations
# - Parameters: sample point spacing, buffer distance for extraction, naming 
# Constraints:
# - Cannot accommodate multipart lines or polygons, converted to singlepart in script
####################################################################################################
#-> Set data paths
setwd('C:/Users/UserName/WorkingDirectory')
####################################################################################################

############################################START MESSAGE###########################################
cat('Line analysis\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Errors and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages (and install if necessary)
pd <- .libPaths()[1]
#pd <- 'C:/Users/UserName/R/R-4.4.0/library' # Or, specify library directory
packages <- c('terra','plyr','readxl')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
		require(package,lib.loc=pd,character.only=T)
	}
}

#-> Load in settings
settings <- data.frame(read_excel('01_Line_analysis_settings.xlsx',sheet='Settings'))
run_name <- settings[settings$Setting=='Run name','Value']
run_name <- gsub(',','',gsub(' ','_',run_name)) # Clean up for file naming
bdist <- as.numeric(settings[settings$Setting=='Buffer','Value'])
pdist <- as.numeric(settings[settings$Setting=='Point spacing','Value'])
rasters <- data.frame(read_excel('01_Line_analysis_settings.xlsx',sheet='Rasters'))
rasters <- rasters[rasters$Include==1,]
alines <- data.frame(read_excel('01_Line_analysis_settings.xlsx',sheet='Lines'))
alines <- alines[alines$Include==1,]

#-> Set projection for perimeter measurements and exports
proj <- paste('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83',
              '+units=m +no_defs')

#-> Load in analysis lines; store in list
AL.l <- list() # List to store analysis lines
ALE.l <- list() # List to store buffered extents for clipping rasters
for(i in 1:nrow(alines)){
	al <- suppressWarnings(vect(alines$Layer[i]))
	al <- disagg(al) # Multipart to singlepart
	al <- project(al,proj)
	if(alines$Type[i]=='polygon'){
		al <- buffer(al,width=0) # Mitigate any geometry errors
		al <- as.lines(al)
	}
	al$LID <- seq(1,nrow(al),1) # Add line ID field
	al$Length_mi <- perim(al)*0.000621371 # Get length of each line
	AL.l[[i]] <- al
	ALE.l[[i]] <- buffer(al,width=bdist*2) # Apply buffer
}

#-> Get extent needed for rasters
if(length(ALE.l) > 1){
	alextents <- do.call('rbind',ALE.l)
}else{
	alextents <- ALE.l[[1]]
}

#-> Load in rasters; process if needed; store in list
R.l <- list()
for(i in 1:nrow(rasters)){
	#-> Crop
	R.l[[i]] <- crop(rast(rasters$Raster[i]),alextents)
	#-> Reclassify
	if(!is.na(rasters$Reclass[i])){
		classes <- unlist(strsplit(rasters$Reclass[i],';'))
		from <- NA; to <- NA
		for(k in 1:length(classes)){
			csplit <- unlist(strsplit(classes[k],' '))
			from[k] <- as.numeric(csplit[1]); to[k] <- as.numeric(csplit[2])
		}
		R.l[[i]] <- classify(R.l[[i]],rcl=data.frame(from,to)) 
	}
	#-> Replace NA with zero
	if(!is.na(rasters$NA2Zero[i])){
		R.l[[i]][is.na(R.l[[i]])] <- 0
	}
}

#-> Function to generate regularly-spaced points along line centerlines
# inLine = single polyline with line id field (LID)
# pdist = point spacing distance in meters
# returns sample points for line tagged with LID and sequential point IDs (PID)
regPoints <- function(inLine,pdist){
	lcdf <- crds(inLine,df=T) # Get coordinates
	lcdf$pdist <- 0 # Field for point distances
	for(i in 2:nrow(lcdf)){ # Calculate point distances
		lcdf$pdist[i] <- sqrt((lcdf[i,'x'] - lcdf[(i-1),'x'])^2 + (lcdf[i,'y'] - lcdf[(i-1),'y'])^2)
	}
	lcdf$cdist <- cumsum(lcdf$pdist) # Calculate cumulative distance
	npoints <- round(max(lcdf$cdist)/pdist,0) # Calculate number of points to generate
	if(npoints == 0){ # Means line is shorter than desired point spacing
		bds <- max(lcdf$cdist)/2 # Use midpoint for break distance
	}
	if(npoints > 0){ # Means line is at least as long as desired point spacing
		rem <- (max(lcdf$cdist) - npoints*pdist)/2 # Get remainder
		bds <- seq(pdist/2 + rem,max(lcdf$cdist) - pdist/2 - rem,pdist) # Calc break distances
	}
	spdf.l <- list()
	for(i in 1:length(bds)){
		scp <- max(which(lcdf$cdist <= bds[i])) # Starting calculation point for interpolation
		if(lcdf$cdist[scp] < bds[i]){ # Interpolate
			sp <- lcdf[scp,]; ep <- lcdf[(scp+1),] # Starting point and end point
			xdiff <- ep$x - sp$x # x difference
			ydiff <- ep$y - sp$y # y difference
			if(xdiff != 0){
				m <- ydiff/xdiff # Calculate slope
				z <- bds[i] - lcdf$cdist[scp] # Calculate z (hypotenuse)
				x <- sqrt((z^2)/abs(1+m^2)) # Solve for x
				if(xdiff < 0){ # Add correct sign to x if it was negative
					x <- x*(-1)
				}
				y <- x*m # Solve for y
			}else{
				x <- 0
				y <- ifelse(ydiff > 0,bds[i] - lcdf$cdist[scp],-1*(bds[i] - lcdf$cdist[scp]))
			}
			spdf.l[[length(spdf.l)+1]] <- data.frame(x=sp$x+x,y=sp$y+y,cdist=bds[i]) # Save point
		}else{ # Use exact point if match
			spdf.l[[length(spdf.l)+1]] <- lcdf[scp,] # Use exact point
		}
	}
	spdf <- do.call('rbind',spdf.l) # Compile to data frame
	sps <- vect(as.matrix(spdf[,c('x','y')]),crs=crs(inLine),type='points') # Convert to spatial
	sps$LID <- inLine$LID # Transfer line ID
	sps$PID <- seq(1,nrow(sps),1) # Generate point ID
	return(sps) # Return points
}


#############################################END SET UP#############################################

###########################################START ANALYSIS###########################################

#-> Create output directory
dir.create(paste0('./',run_name))

###---> Process analysis lines and extract raster values

for(i in 1:nrow(alines)){
	
	#-> Convert lines to points for analysis
	# Uses regPoints function defined above
	al <- AL.l[[i]] # Subset analysis lines from list
	alp.l <- list() # List to store results
	for(j in 1:nrow(al)){ # Iterate through lines
		alp.l[[j]] <- regPoints(al[j,],pdist) # Generate regular points		
	}
	alp <- do.call('rbind',alp.l) # Compile results
	alpb <- buffer(alp,width=bdist) # Create associated buffers for summary stats
	
	#-> Extract raster data
	for(j in 1:nrow(rasters)){
		alpb <- project(alpb,crs(R.l[[j]])) # Match projections
		X <- extract(R.l[[j]],alpb,fun=rasters$Extract_function[j],na.rm=T)[,2] # Extract
		alpb$X <- round(X,rasters$Decimals[j]) # Round
		names(alpb)[ncol(alpb)] <- rasters$Name[j] # Rename field
	}
	alp <- merge(alp,data.frame(alpb),by=c('LID','PID'),all.x=T)
	
	#-> Analysis line point stats
	stats.l <- list() # List for storing results
	for(j in 1:nrow(rasters)){ # Iterate through rasters
		bins <- as.numeric(unlist(strsplit(rasters$Bins[j],','))) # Extract bins
		bin_labs <- unlist(strsplit(rasters$R_bin_labels[j],',')) # Extract bin labels
		x <- hist(data.frame(alp)[,rasters$Name[j]],breaks=bins,plot=F)$counts # Get counts
		xdf <- data.frame(Var=rasters$Name[j],Bin_label=bin_labs,Range=
		                  paste(bins[1:(length(bins)-1)],bins[2:length(bins)],sep='-'),Count=x)
		xdf$Percent <- round(100*(xdf$Count/sum(xdf$Count)),1) # Get percents
		xdf$Length_mi <- sum(al$Length_mi)*(xdf$Percent/100) # Get lengths by percent of total
		stats.l[[j]] <- xdf # Save to list
	}
	write.csv(do.call('rbind',stats.l),paste0('./',run_name,'/',run_name,'_',
	          gsub(' ','_',alines$Map_name[i]),'_stats.csv'),row.names=F)
	
	#-> Save analysis line points for mapping
	alp  <- project(alp,proj)
	alp_save <- alp
	for(j in 1:nrow(rasters)){ # Recode NA so it is properly exported to shapefile
		X <- data.frame(alp_save)[,rasters$Name[j]]
		X[is.na(X)] <- -1
		alp_save[,rasters$Name[j]] <- X
	}
	writeVector(alp_save,filename=paste0('./',run_name,'/',run_name,'_',gsub(' ','_',
	            alines$Map_name[i]),'_points.shp'),overwrite=T)
				
	#-> Attribute analysis lines with min, mean, and max of each metric
	for(j in 1:nrow(rasters)){
		alp$X <- data.frame(alp)[,paste(rasters$Name[j])]
		xdf <- ddply(data.frame(alp),.(LID),summarize,
					 min = min(X,na.rm=T),
					 mean = round(mean(X,na.rm=T),rasters$Decimals[j]+2),
					 max = max(X,na.rm=T))
		al <- merge(al,xdf,by='LID',all.x=T)
		al$min[is.na(al$mean)] <- -1 # Recode NA so it is properly exported to shapefile
		al$max[is.na(al$mean)] <- -1
		al$mean[is.na(al$mean)] <- -1
		names(al)[(ncol(al)-2):ncol(al)] <- paste(rasters$Name[j],c('min','mean','max'),sep='_')
		alp$X <- NULL
	}	
	
	#-> Save analysis lines for mapping and analysis
	al  <- project(al,proj)
	writeVector(al,filename=paste0('./',run_name,'/',run_name,'_',
	            gsub(' ','_',alines$Map_name[i]),'_lines.shp'),overwrite=T)				
					
}
	
############################################END ANALYSIS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END LOGGING#############################################

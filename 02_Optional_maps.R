####################################################################################################
#### Line analysis maps and graphics
#### Author: Ben Gannon (benjamin.gannon@usda.gov)
#### Date Created: 09/26/2021
#### Last Modified: 05/30/2024
####################################################################################################
# Makes simple maps for each line summary attribute. Run after completing the line analysis. Note
# that title, legend, and arrow placement is sensitive to the R and terra package versions - some
# adjustments may be needed.
####################################################################################################
#-> Set data paths
setwd('C:/Users/UserName/WorkingDirectory')
####################################################################################################

############################################START MESSAGE###########################################
cat('Line analysis maps and graphics\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Errors and Warnings (if they exist):\n')
####################################################################################################

############################################START SETUP#############################################

#-> Load packages
pd <- .libPaths()[1]
#pd <- 'C:/Users/UserName/R/R-4.4.0/library' # Or, specify library directory
packages <- c('terra','plyr','readxl','basemaps','magick','grid','gridExtra')
for(package in packages){
	if(suppressMessages(!require(package,lib.loc=pd,character.only=T))){
		install.packages(package,lib=pd,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,lib.loc=pd,character.only=T))
		require(package)
	}
}

#-> Load in settings
rasters <- data.frame(read_excel('01_Line_analysis_settings.xlsx',sheet='Rasters'))
rasters <- rasters[rasters$Include==1,]
settings <- data.frame(read_excel('01_Line_analysis_settings.xlsx',sheet='Settings'))
run_name <- settings[settings$Setting=='Run name','Value']
run_name_file <- gsub(',','',gsub(' ','_',run_name)) # Clean up for file naming
map_text <- settings[settings$Setting=='Map text','Value']
map_text <- unlist(strsplit(map_text,','))

#-> Load in results of spatial analysis
clps <- vect(paste0('./',run_name_file,'/',run_name_file,'_Planned_lines_points.shp'))
clsum <- read.csv(paste0('./',run_name_file,'/',run_name_file,'_Planned_lines_stats.csv'),
                  header=T)

#-> Reference layers
fire <- vect('./Analysis_lines/perimeter.shp')

#-> Set global parameters
bdist <- 800 # Buffer distance in meters around spatial data for mapping
proj <- paste('+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m',
              '+nadgrids=@null +wktext +no_defs')

#-> Define make square function to specify imagery extent
makeSquare <- function(inPoly){
	xmin <- as.numeric(ext(inPoly)[1])
	xmax <- as.numeric(ext(inPoly)[2])
	ymin <- as.numeric(ext(inPoly)[3])
	ymax <- as.numeric(ext(inPoly)[4])
	diff <- abs((xmax - xmin) - (ymax - ymin))
	if((xmax - xmin) <= (ymax - ymin)){
		xmin <- xmin - diff/2
		xmax <- xmax + diff/2
	}else{
		ymin <- ymin - diff/2
		ymax <- ymax + diff/2
	}
	bb <- as.polygons(ext(c(xmin,xmax,ymin,ymax)))
	crs(bb) <- crs(inPoly)
	return(bb)
}

#-> Define color transparency function
add.alpha <- function(COLORS,ALPHA){
   RGB <- col2rgb(COLORS,alpha=T)
   RGB[4,] <- round(RGB[4,]*ALPHA)
   return(rgb(RGB[1,],RGB[2,],RGB[3,],RGB[4,],maxColorValue=255))
}

#############################################END SETUP##############################################

##########################################START ANALYSIS############################################

###---> Organize data and download background layer

#-> Project spatial data
clps <- suppressWarnings(project(clps,proj))
fire <- suppressWarnings(project(fire,proj))

#-> Create bounding box
mext <- as.polygons(ext(clps))
crs(mext) <- proj
bb <- makeSquare(buffer(mext,width=bdist))

#-> Get base map
# Use map_res to adjust resolution as needed bigger number means higher resolution
# Values of 1 to 3 seem to work well for most fires
#base <- basemap_raster(ext=bb,map_service='esri',map_type='world_topo_map',verbose=F)
base <- rast(basemap_stars(ext=bb,map_service='esri',map_type='world_topo_map',verbose=F))

#-> Crop data to bounding box
fire <- crop(buffer(fire,width=0),bb)

for(i in 1:nrow(rasters)){
	
	#-> Organize data
	brks <- as.numeric(unlist(strsplit(rasters$Bins[i],',')))
	labs <- unlist(strsplit(rasters$R_bin_labels[i],','))
	colRP <- unlist(strsplit(rasters$R_ColorPalette[i],','))
	cols <- colorRampPalette(colRP)(length(brks)-1)
	clps$X <- data.frame(clps)[,rasters$Name[i]]
	
	#-> Make map
	tiff('map.tif',width=2000,height=2000,compression='lzw',type='windows',pointsize=44)
	plot(bb,border=NA,axes=F,mar=c(0.1,0.1,2,0.1)) # Set desired extent for plot
	plotRGB(base,add=T) # Basemap
	mtext(run_name,cex=1.75,line=2,adj=0.05,font=2)
	mtext(rasters$Map_name[i],cex=1.75,line=2,adj=0.95,font=2)
	plot(bb,add=T) # Plot boundary
	plot(fire,col=add.alpha('grey60',0.8),border=NA,add=T)
	plot(clps,pch=20,col=cols[findInterval(clps$X,brks,all.inside=T)],cex=0.5,add=T)
	#-> Scale bar
	sbx <- xmin(bb) + 0.025*(xmax(bb)-xmin(bb))
	sby <- ymin(bb) + 0.025*(ymax(bb)-ymin(bb))
	lbins <- c(0,1600,3200,8000,16000,32000,100000,1000000)
	llens <- c(400,800,1600,3200,8000,16000,32000)
	llabs <- c('1/4 mile','1/2 mile','1 mile','2 miles','5 miles','10 miles','20 miles')
	sbi <- findInterval((xmax(bb)-xmin(bb)),lbins)
	sbar(llens[sbi],xy=c(sbx,sby),type='line',divs=2,lwd=8,label=llabs[sbi])
	#-> North arrow
	y1 <- ymin(bb) + 0.075*(ymax(bb)-ymin(bb))
	y2 <- ymin(bb) + 0.125*(ymax(bb)-ymin(bb))
	arrows(sbx,y1,sbx,y2,length=0.25,lwd=8,xpd=T)
	#-> Legend
	legend('bottomright',c('Fire'),pch=c(20),fill=c('grey70'),border=c(NA),col=c(NA),lwd=c(NA),
	       bg='white',xpd=T)
	g <- dev.off()
	
	#-> Text block
	tiff('text.tif',width=614,height=440,compression='lzw',
	     type='windows',pointsize=22)
	t1 <- ttheme_minimal(base_size=40,core=list(fg_params=list(fontface=c('bold',rep('plain',
	                     length(map_text)-1)),hjust=0,x=0)))
	grid.arrange(tableGrob(map_text,rows=NULL,theme=t1))
	g <- dev.off()
	
	#-> Summary table
	# Can be expanded to communicate line length or estimate crew hours by category
	ot <- data.frame(labs)
	colnames(ot) <- c(rasters$Map_name[i])
	tiff('table.tif',width=1000,height=440,compression='lzw',
	     type='windows',pointsize=22)
	t2 <- ttheme_default(base_size=40,core=list(bg_params=list(fill=add.alpha(cols,0.6))))
	grid.arrange(tableGrob(ot,rows=NULL,theme=t2))
	g <- dev.off()
	
	#-> Assemble panels with magick
	map_p <- image_read('map.tif')
	text_p <- image_read('text.tif')
	table_p <- image_read('table.tif')
	RMA_p <- image_read('RMA_small.tif')
	lp <- image_append(c(text_p,table_p,RMA_p),stack=F)
	full <- image_append(c(map_p,lp),stack=T)
	image_write(full,paste0('./',run_name_file,'/',run_name_file,'_Map_',rasters$Name[i],'.jpg'),
	            format='jpeg')
	g <- file.remove(c('text.tif','table.tif','map.tif'))
	
}

###########################################END ANALYSIS#############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Close command window to proceed!\n',sep='')
############################################END LOGGING#############################################

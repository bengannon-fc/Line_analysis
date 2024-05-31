# Line analysis
This line analysis workflow is designed to extract raster data at fixed intervals along lines to generate line level summary metrics. It was developed for wildfire management to rate potential fire containment lines on multiple dimensions of opportunity, difficulty, and safety. The core sampling, extraction, and summary workflow may have other applications in natural resources management.

**Sampling framework**

The user provides singlepart polylines of interest for analysis. An option is also available for providing polygons that will be converted to perimeter polylines for analysis. Sample points are generated at fixed intervals along the lines. Any multipart polylines will be converted to singlepart because multipart polylines are fundamentally incompatible with the concept of fixed interval point sampling. Each polyline and associated sample points are assigned matching line identifier (LID) attributes for relating their data tables. Raster attributes of interest are then extracted within a buffer zone around each sample point. Each point is assigned a single value for each raster using the summary statistic of interest - generally, the mean makes sense for continuous data rasters and the maximum makes sense for binary presence absence rasters. Polyline summary statistics are then calculated from the points associated with each line including the minimum, mean, and maximum. The main output is shapefiles of the sample points and analysis lines with the extracted raster data in their attribute tables.

![image](https://github.com/bengannon-fc/Line_analysis/assets/81584637/f2f95ab4-6610-4f74-9378-c8f072675a85)

**Instructions for Use**

The line analysis script was developed and tested using the terra package version 1.7-78 (Hijmans 2024) in the R language for statistics and computing version 4.4.0 (R Core Team 2024). The script does not include advanced error handling. The end user is responsible for verifying that their inputs are of the correct type, format, and units. 

_Line analysis_

The line analysis workflow is spreadsheet-driven. The only line in the script that a user should need to modify is the working directory. Change the working directory path to match the location of the 01_Line_analysis.R script and 01_Line_analysis_settings.xlsx file. The vector sample points are reprojected within the script to match each raster as it is analyzed, so there is no need to reproject and/or resample data to match. The only requirement is that each input has a defined coordinate reference system.

The Settings worksheet is used to set the following parameters:
1) Run name (text field used for naming output directory, files, and optional maps)
2) Map text (optional map text for providing context separated by commas to break lines)
3) Buffer size in meters around sample points for raster data extraction - recommended value is 60 meters for typical applications with 30 meter resolution raster data
4) Point spacing in meters for sample point generation - recommended value is 100 meters for typical large fire analysis applications

The Lines worksheet is used to specify the vector lines to analyze:
1) Use a row for each input dataset
2) Specify the layer name for mapping, data type, and location as prompted by the comments in the header row 

The Rasters worksheet is used to specify the raster data to summarize:
1) Use a row for each input dataset
2) Specify the layer name for attributes, layer name for mapping, location, extraction function, any transformations to apply, and optional map symbology preferences as prompted by the comments in the header row 

_Optional maps_

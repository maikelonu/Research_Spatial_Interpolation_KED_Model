# Precipitation Interpolation Methods
# Instituto Tecnologico de Costa Rica (www.tec.ac.cr)
# Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)
# Luis Alexander Calvo-V (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)
# This script is structured in R (www.r-project.org)
# General purpose: Generate temporal series of average precipitation for a waterbasin using
# the geostatistical interpolation method of Universal Kriging (UK)
# Input files: "calibration.dat", "validation.dat", "interpolation.dat", "blank.asc", "attri.txt", "spatial.dat"
# Output files: "output_interpolation.csv"

# working directory is defined
setwd ("B:\\R_ITC\\OPTI_SCRIPTS\\INTER_UK\\INTER_UK_GITHUB_02_OCT_2015")

# source() and library() statements
require(automap)
require(doParallel)
require(doSNOW)
require(foreach)
require(gstat)
require(lattice)
require(maptools)
require(parallel)
require(pastecs)
require(plyr)
require(raster)
require(rgdal)
require(rgeos)
require(snow)
require(sp)

# Reading various spatial basin-attributes values and parameters
     attri <- read.delim("attri.txt", header = TRUE, sep = "\t")  # Atributes file
   spa.res <- as.numeric(attri[1, 2])  # Selected spatial resolution (m)
      xmin <- as.numeric(attri[2, 2])  # Minimun X extension
      xmax <- as.numeric(attri[3, 2])  # Maximum X extension
      ymin <- as.numeric(attri[4, 2])  # Minimun Y extension
      ymax <- as.numeric(attri[5, 2])  # Maximum Y extension
  goef.col <- as.numeric(attri[6, 2])  # Number of columns in the georeference
  geof.row <- as.numeric(attri[7, 2])  # Number of rows in the georeference
 threshold <- as.numeric(attri[8, 2])  # Minimum precipitation-interpolation threshold (mm)
nugget.all <- as.numeric(attri[9, 2])  # Nugget value for AutoKrige method using all stations
 range.all <- as.numeric(attri[10, 2])  # Range value for AutoKrige method using all stations
  sill.all <- as.numeric(attri[11, 2])  # Sill value for AutoKrige method using all stations
nugget.cal <- as.numeric(attri[12, 2])  # Nugget value for AutoKrige method using a calibration subset of the all stations
 range.cal <- as.numeric(attri[13, 2])  # Range value for AutoKrige method using a calibration subset of the all stations
  sill.cal <- as.numeric(attri[14, 2])  # Sill value for AutoKrige method using a calibration subset of the all stations
     cores <- as.numeric(attri[15, 2])  # Number of hardware cores for cluster

cluster <- makeCluster(cores)  # The cluster is defined with a number of cores
registerDoParallel(cluster)  # The cluster is registered

# Reading input files and creating data.frames
inter.all <- read.table("interpolation.dat", header = T, )  # Interpolation file, all stations are included
inter.spt <- read.table("spatial.dat", header = T, )  # UK complementary data; Z, ASPECT, SLOPE. Prepared externally by GIS

if (file.exists("calibration.dat") & file.exists("validation.dat")){  # checks if there is not a subset of the all stations
  inter.cal <- read.table("calibration.dat", header = T, )  # Calibration file, only contains calibration stations
  inter.val <- read.table("validation.dat", header = T, )  # Validation file, only contains validation stations
} else {
  inter.cal <- read.table("interpolation.dat", header = T, )  # Calibration file, only contains calibration stations
  inter.val <- read.table("interpolation.dat", header = T, )  # Validation file, only contains validation stations
}

# Defining spatial data.frame structure
colClasses <- c("double", "double", "double", "double", "double", "double")
 col.names <- c("X", "Y", "Z", "ASPECT", "SLOPE","DATE")
 frame.all <- read.table(text = "", colClasses = colClasses, col.names = col.names)  # Interpolation spatial data.frame
 frame.cal <- read.table(text = "", colClasses = colClasses, col.names = col.names)  # Calibration spatial data.frame
 frame.val <- read.table(text = "", colClasses = colClasses, col.names = col.names)  # Validation spatial data.frame

# Reading input data.frames and assigning counters
   n.obs.all <- nrow(inter.all)  # number of data rows for interpolation
n.points.all <- ncol(inter.all)  # number of columns for interpolation
   n.obs.cal <- nrow(inter.cal)  # number of data rows for calibration
n.points.cal <- ncol(inter.cal)  # number of columns for calibration
   n.obs.val <- nrow(inter.val)  # number of data rows for validation
n.points.val <- ncol(inter.val)  # number of columns for validation

# Completing values for interpolation spatial data.frame 
for (i in 1:(n.points.all - 1)) {
  frame.all[i, 1] <- inter.all[1, 1 + i]
  frame.all[i, 2] <- inter.all[2, 1 + i]
  frame.all[i, 3] <- inter.all[3, 1 + i]
  frame.all[i, 4] <- inter.all[4, 1 + i]
  frame.all[i, 5] <- inter.all[5, 1 + i]
  frame.all[i, 6] <- inter.all[6, 1 + i]
}

# Completing values for calibration spatial data.frame
for (i in 1:(n.points.cal - 1)) {
  frame.cal[i, 1] <- inter.cal[1, 1 + i]
  frame.cal[i, 2] <- inter.cal[2, 1 + i]
  frame.cal[i, 3] <- inter.cal[3, 1 + i]
  frame.cal[i, 4] <- inter.cal[4, 1 + i]
  frame.cal[i, 5] <- inter.cal[5, 1 + i]
  frame.cal[i, 6] <- inter.cal[6, 1 + i]
}

# Completing values for validation spatial data.frame
for (i in 1:(n.points.val - 1)) {
  frame.val[i, 1] <- inter.val [1, 1 + i]
  frame.val[i, 2] <- inter.val [2, 1 + i]
  frame.val[i, 3] <- inter.val [3, 1 + i]
  frame.val[i, 4] <- inter.val [4, 1 + i]
  frame.val[i, 5] <- inter.val [5, 1 + i]
  frame.val[i, 6] <- inter.val [6, 1 + i]
}

# Coordinates are assigned and transformed into SpatialPointsDataFrame
coordinates(frame.all) <- c("X", "Y")
coordinates(frame.cal) <- c("X", "Y")
coordinates(frame.val) <- c("X", "Y")

# Empty UK complementary data SpatialGrid is converted into SpatialPixels
# This example has a 100 x 100 m spatial resolution. It varies with the waterbasin under analysis 
gridded(inter.spt) =~ X + Y

# An empty data.frame is created with a selected resolution
# The extent of this SpatialGrid is valid only for waterbasin under analysis
georef <- expand.grid(X = seq(xmin, xmax, by = spa.res), Y = seq(ymin, ymax, by = spa.res), KEEP.OUT.ATTRS = F)

coordinates(georef) <- ~X + Y  # The empty data.frame is converted into SpatialPoints

# An empty SpatialGrid is created with a selected resolution
# The extent of this SpatialGrid is valid only for waterbasin under analysis
         basin.grid <- SpatialGrid(GridTopology(c(X = xmin, Y = ymin), c(spa.res, spa.res), c(goef.col, geof.row)))
         basin.grid <- SpatialPoints(basin.grid)  # The empty SpatialGrid is converted into SpatialPoints
gridded(basin.grid) <- T  # The empty SpatialGrid is converted into SpatialPixels

       blank <- paste("blank.asc", sep = "/")  # The *.ASC map is imported from GIS to be used as a black or mask
  blank.grid <- read.asciigrid(blank, as.image = FALSE)  # The map is transformed into SpatialGridDataFrame
blank.raster <- raster(blank.grid)  # The map is transformed into RasterLayer
image(blank.raster, main = "*.ASC Blanking Map", col = topo.colors(64)) # The map is printed as verification

# data.frame that will contain the output of the interpolation process is defined
parallel.container <- NULL  # Container of the parallel process 

# The start time of the process is saved
start.time <- Sys.time()

# Main counters is defined
counter.main <- (n.obs.all - 5)

# Main parallel-loop is initialized
parallel.container <- foreach(i = 1:counter.main, .combine = 'rbind', 
                              .packages = c("gstat", "sp", "raster", "automap"), .multicombine = TRUE) %dopar% {
  
# Completing values for interpolation spatial data.frame within the main loop
for (a in 1:(n.points.all - 1)) {
  frame.all [a,4] <- inter.all [i + 5, a + 1]
}
  
# Completing values for calibration spatial data.frame within the main loop
for (c in 1:(n.points.cal - 1)) {
  frame.cal [c,4] <- inter.cal [i + 5, c + 1]
}
  
# Completing values for validation spatial data.frame within the main loop
for (v in 1:(n.points.val - 1)) {
  frame.val [v,4] <- inter.val [i + 5, v + 1]
}
  
# The mean of the vertical data.frame for interpolation should be above the defined threshold
if ((mean(frame.all$DATE)) > threshold) {
    
# Total Interpolation 
        uk.all <- autoKrige((frame.all$DATE)~Z + ASPECT + SLOPE, frame.all, inter.spt, 
                            fix.values = c(nugget.all, range.all, sill.all)) # Nugget, Range and Sill can be fixed
    uk.out.all <- uk.all$krige_output  # Output instructions of Autokrige are defined
 uk.all.raster <- raster(uk.out.all)  # SpatialPixelsDataFrame is converted into RasterLayer
    resampling <- resample(blank.raster, uk.all.raster, resample = 'bilinear')  # The resampling is defined
   mask.uk.all <- mask(uk.all.raster, resampling)  # RasterLayers are created after cutting
   mean.uk.all <- cellStats(mask.uk.all, mean)  # The "mean" is extracted from RasterLayer 
rounded.uk.all <- round(mean.uk.all, 4)  # The mean is rounded to four significant digits
    
  if (rounded.uk.all < 0)  # Negative values are discarded
      rounded.uk.all <- 0
    
  # Interpolation Calibration
          uk.cal <- autoKrige((frame.cal$DATE)~Z + ASPECT + SLOPE, frame.cal, inter.spt, 
                              fix.values = c(nugget.all, range.all, sill.all)) # Nugget, Range and Sill can be fixed
      uk.out.cal <- uk.cal$krige_output  # Output instruction is defined 
   uk.cal.raster <- raster(uk.out.cal)  # SpatialPixelsDataFrame is converted into RasterLayer   
     mask.uk.cal <- mask(uk.cal.raster, resampling)  # RasterLayers are created after cutting
     mean.uk.cal <- cellStats(mask.uk.cal, mean)  # The "mean" is extracted from RasterLayer 
  rounded.uk.cal <- round(mean.uk.cal, 4)  # The mean is rounded to four significant digits
    
    if (rounded.uk.cal < 0)  # Negative values are discarded
        rounded.uk.cal <- 0
    
    # Point cross validation
    raster.georef <- raster(inter.spt)  # The grid is converted into RasterLayer
       ext.raster <- frame.val[4]  # Date and hour are extracted
     point.raster <- rasterize(ext.raster, raster.georef)  # Point values are rasterized
    
    # The layers are disaggregated from the point.raster
    upper.layer <- unstack(point.raster)
    lower.layer <- upper.layer[[2]]
    
    # The Mean Absolute Error is calculated
            mae.uk <- abs(lower.layer - (uk.cal.raster))
       mean.mae.uk <- cellStats(mae.uk, (mean))  # The mean value of deviations is determined 
    rounded.mae.uk <- round(mean.mae.uk, 4)  # The mean is rounded to four significant digits
    
    # The RMSE is calculated
            rmse.uk <- (lower.layer - (uk.cal.raster)) ^ 2
       mean.rmse.uk <- (cellStats(rmse.uk, (mean))) ^ 0.5  # The mean value of deviations is determined 
    rounded.rmse.uk <- round(mean.rmse.uk, 4)  # The mean is rounded to four significant digits
  }
  else {
     rounded.uk.all <- 0
     rounded.uk.cal <- 0
     rounded.mae.uk <- 0
    rounded.rmse.uk <- 0
  }
         date.charac <- as.character (inter.all [i + 5, 1])  # Date is saved as character
  parallel.container <- c(date.charac, rounded.uk.all, rounded.uk.cal, rounded.mae.uk, 
                          rounded.rmse.uk)  # Output matrix structure is defined
}
            temp.character <- parallel.container[, 1]  # A matrix containing "DATES" as character is created
              temp.numeric <- parallel.container[, 2 : 5]  # A matrix containing all other outputs as character is created
              temp.numeric <- apply(temp.numeric, c(1 : 2), as.numeric)  # Matrix is tranformed to numeric
               out.numeric <- as.data.frame(temp.numeric)  # A data.frame containing "DATES" as character is created
             out.character <- as.data.frame(temp.character)  # A matrix containing all other outputs as numeric is created
          output.dataframe <- cbind(out.character, out.numeric)  # Data.frames are binded
colnames(output.dataframe) <- c("DATE", "UK.ALL", "UK.CAL", "MAE.UK", "RMSE.UK")  # Columns of the output dataframe are renamed

write.csv(output.dataframe, file = "output_general.csv")  # Output dataframe is saved

# A descriptive statistics data.frame for output is created and saved
descri.output <- round(stat.desc(output.dataframe [, 2 : 5]), 3)
write.csv(descri.output, file = "descri_output.csv")

# End time is saved
  end.time <- Sys.time()  # End time is saved
time.taken <- end.time - start.time
stopCluster(cluster)
print(time.taken)

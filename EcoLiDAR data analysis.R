#Calebe Pereira Mendes
#2023-04-28
#calebepm3@gmail.com
#calebepm3@hotmail.com
#https://github.com/CalebePMendes/EcoLiDAR

### This script is divided in 2 segments.
  #The first  demonstrate how the raw data can be turned into a .las file.
  #The second show some analysis performed using the package lidR




###################################################################################################################################
##### Handling the EcoLiDAR raw data
###################################################################################################################################

rm(list = ls()) #clean the environment

require(tidyverse)
require(rlas)
require(lidR)
options(digits=22) #to avoid scientific notation



setwd("folder path")
list.files()

source("EcoLiDAR Functions.r") #Load some functions to help




### LOAD THE RAW SCAN

Lidar_raw = read.csv("raw_scan_1.csv") #load the raw scan file

Lidar_raw = laser_offset(Lidar_raw,2, right_side = T) #correct the 2 cm offset from the laser rangefinder and the axis of rotation of the rotating head.

Lidar_raw$distance = Lidar_raw$distance/100 #Change the units of distance from centimeter to meter

Lidar_raw = Lidar_raw[Lidar_raw$distance <= 40,] ##optional, exclude the points farther than 40 m
#the laser rangefinder maximum range specification is 40 meters. It can detect things farther but with reduced precision. For demonstration purposes, lets limit it to 40 meters.

Lidar_raw$time = Lidar_raw$time - Lidar_raw$time[1] ##Subtract the time at the scan start, so all timestamps are the same since the start of the scan.

Lidar_raw_mirrored = mirror_raw(Lidar_raw) #for demonstration only.
# The orientation if the direction of the spinning head is not clockwise it will not match the output of the python software provided. This function fix the issue by mirroring the bearing directions.
# Compare the bearings between Lidar_raw and Lidar_raw_mirrored.
rm(Lidar_raw_mirrored)




### Visualize the scan 
pic = angtopic(Lidar_raw) #Attention!! This function is atrociously slow!! it works...but it would benefit from some serious optimization.
pic






### CREATE A XYZ FILE

Lidar_xyz = angtoxyz(Lidar_raw)




### ADD THE SCAN COORDINATES (optional)

#The coordinates needs to be taken in field as part of the scan protocol using a GPS unit, as described in the manual.
#As long as the scan is initiated pointing to the north, as described in the manual, the axis work similar to the axis in the UTM coordinate system. x increases from west to east, y increases from south to north, z is elevation. 
  #The center of rotation of the 2 axis gimbals have coordinates x=0, y=0, z=0.
  #as long as the GPS coordinates are in a metric system (ex:EPSG:4087), it can be simply sum to the scan xyz axis.
  # In case that subcentimeter precision is required from the GPS unit (unlikely since the LiDAR itself don`t have subcentimiter precision), remember to consider possible offsets if the GPS unit is not located at the 2 axis gimbals center of rotation during the coordinate acquisition.


#ex: 1.349234 N, 103.679563 W, altitude 40.23 m

#The coordinates can be converted using the sf package, or our convenient function.
coordinates = Latlongto4087(1.349234,103.679563)
coordinates

#Sum the coordinates to the lidar scan
Lidar_xyz$x = Lidar_xyz$x + coordinates[1,"X"] #sum the x axis
Lidar_xyz$y = Lidar_xyz$y + coordinates[1,"Y"] #sum the y axis
Lidar_xyz$z = Lidar_xyz$z + 40.23 #sum the altitude




### CREATE A .LAS FILE

#Create the data using the rlas package
lasdata = data.frame(X = round(Lidar_xyz$x, digits = 4),
                     Y = round(Lidar_xyz$y, digits = 4),
                     Z = round(Lidar_xyz$z, digits = 4),
                     gpstime = Lidar_xyz$time,
                     ReturnNumber = rep(1L, nrow(Lidar_xyz)),
                     NumberOfReturns = rep(1L, nrow(Lidar_xyz)),
                     ScanDirectionFlag = rep(0L, nrow(Lidar_xyz)),
                     EdgeOfFlightline = rep(0L, nrow(Lidar_xyz)),
                     Classification = rep(0L, nrow(Lidar_xyz)),
                     ScanAngleRank = rep(0L, nrow(Lidar_xyz)),
                     UserData = rep(0L, nrow(Lidar_xyz)),
                     PointSourceID = rep(1L, nrow(Lidar_xyz)))

#header
lasheader = header_create(lasdata)

#check validity and compliance
check_las_validity(lasheader, lasdata)
check_las_compliance(lasheader, lasdata)

write.las("Lidar_pointcloud.las", lasheader, lasdata) #save the new .las file
Lidar_pointcloud = readLAS("Lidar_pointcloud.las") #read the .las file using the lidR package

plot(Lidar_pointcloud) #plot for visualization


rm(lasdata, lasheader)




#The .las files are designed for laser point clouds and is accepted by most LiDAR analysis software.
#For now, I`m not aware of any r package for point cloud alignment, so this will not be covered here. I suggest the software CloudCompare for basic point cloud handling and alignment.



###################################################################################################################################
##### Basic analysis of the point cloud
###################################################################################################################################



### Load 4 EcoLiDAR scans which were aligned and merged using the software CloudCompare.
Merged_scan = readLAS("merged.las") #read the .las file using the lidR package

plot(Merged_scan, bg = "black",size = 0.1)




### Ground classification
# Cloth Simulation Function

CSF <- classify_ground(Merged_scan, algorithm = csf(sloop_smooth = T))

gnd <- filter_ground(CSF) #isolate ground points
veg <- filter_poi(CSF, Classification == 0) #isolate the vegetation points


#Visualise
Composite_plot = plot(gnd, color = "Classification",size = 3, bg = "white")
plot(Merged_scan, size = 0.1, add = Composite_plot)




### Digital Terrain Model
# Use the ground point cloud (gnd) and a Triangular Irregular Network algorithm to create a Digital Terrain Model

dtm <- grid_terrain(gnd, res = 1, algorithm = tin()) #make the Digital Terrain Model
plot_dtm3d(dtm, bg = "white") #See the digital terrain model


#plot1
Composite_plot <- plot(veg, bg = "white", size = 0.1)
add_dtm3d(Composite_plot, dtm)


#plot2
Composite_plot = plot(gnd, color = "Classification",size = 3, bg = "white")
plot(Merged_scan, size = 0.1, add = Composite_plot)
add_dtm3d(Composite_plot, dtm)

#plot3
plot(dtm, col = gray.colors(50, 0, 1))




### Normalization
# Makes the soil flat for analytic purpose

normalized_Merged_scan <- normalize_height(Merged_scan, dtm, method = "bilinear", na.rm = T)
plot(normalized_Merged_scan, bg = "black", size = 0.01)


CSF_norm <- classify_ground(normalized_Merged_scan, algorithm = csf()) #Classify the normalized ground
gnd_norm <- filter_ground(CSF_norm) #normalized ground
veg_norm <- filter_poi(CSF_norm, Classification == 0) #normalized vegetation


#plot1
Composite_plot = plot(gnd_norm, color = "Classification",size = 3, bg = "white")
plot(normalized_Merged_scan, size = 0.1, add = Composite_plot)


#make the Digital Terrain Model of the normalized point cloud
dtm_norm <- grid_terrain(gnd_norm, res = 1, algorithm = tin())

#plot
Composite_plot = plot(veg_norm, size = 0.1, bg = "white")
add_dtm3d(Composite_plot, dtm_norm)





### Canopy Height model
#Maps the surface of the canopy
# Similar to TIN, but in the canopy and with some tricks to make it work better

chm <- grid_canopy(normalized_Merged_scan, res = 0.5, algorithm = dsmtin()) #make the Canopy Height model
plot(chm, col = height.colors(50))










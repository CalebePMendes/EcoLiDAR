#Calebe Pereira Mendes
#2023-04-28
#calebepm3@gmail.com
#calebepm3@hotmail.com
#https://github.com/CalebePMendes/EcoLiDAR


### Functions for the EcoLiDAR data


#Raw scan structure:
#four columns:"bearing"   "inclination"     "distance"  "time"
#bearing is the direction (north = 0/360)
#inclination is the vertical angle, in degrees (0 = horizon)
#distance to the target, in meters
#time is not used here, but kept for the other
require(tidyverse)
require(sf)
options(scipen = 999)





#laser_offset (correct the offset of the laser to the center of the vertical axis of the rotating head)

#raw = EcoLidar data scan (data.frame)
#offset = Distance from the laser rangefinder source to the vertical axis center (in cm)
#right_side = Side of the laser rangefinder offset in relation to the vertical axis (when looking toward the front of the rotating head)
  #default: right_side = TRUE

laser_offset = function(raw, offset, right_side){ 
  
  #True_dist
  raw$true_dist = sqrt(raw$distance^2 + offset^2)
  
  #Tan_radian
  raw$tan_radian = offset / raw$distance
  
  #Tan_degrees
  raw$tan_degrees = raw$tan_radian * 180/pi
  
  #True_dist (according with the side)
  if(right_side){
    raw$true_bearing = raw$bearing - raw$tan_degrees
    raw$true_bearing[raw$true_bearing < 0 & !is.infinite(raw$true_bearing)] = 360 + raw$true_bearing[raw$true_bearing < 0 & !is.infinite(raw$true_bearing)]
  }else{
    raw$true_bearing = raw$bearing + raw$tan_degrees
    raw$true_bearing[raw$true_bearing > 360 & !is.infinite(raw$true_bearing)] = raw$true_bearing[raw$true_bearing > 360 & !is.infinite(raw$true_bearing)] %% 360
  }
  
  # Zero is the default return value when the laser bean reflex is not detected (shot beyond measurement range).
  # Since the offset effect is smaller at higher ranges, it become negligible at large ranges.
  raw$true_dist[raw$distance == 0] = 0
  raw$true_bearing[raw$distance == 0] = raw$bearing[raw$distance == 0]

  #save the data
  raw$distance = raw$true_dist
  raw$bearing = raw$true_bearing
  
  #clean
  raw$true_dist = NULL
  raw$tan_radian = NULL
  raw$tan_degrees = NULL
  raw$true_bearing = NULL
  
  return(raw)
}




#Mirror_raw

#invert the direction of the bearing measurements (in case the rotating head is turning in a anti-clockwise direction. ex: 0/360 -> 270 -> 180 -> 90 -> 0/360)

mirror_raw = function(raw){
  
  raw$bearing = 360-raw$bearing
  return(raw)
  
}



#Latlongto4087 (Convert lat long or decimal degrees to EPSG:4087 (WGS 84 / World Equidistant Cylindrical)

Latlongto4087 = function(lat,long){
  
  if(is.character(lat)){ #Convert degrees-minutes-seconds to decimal degrees
    
    lat_deg = as.numeric(strsplit(lat, " ")[[1]][1:3])
    long_deg = as.numeric(strsplit(long, " ")[[1]][1:3])
    
    hemisfere = strsplit(lat, " ")[[1]][4]
    side = strsplit(long, " ")[[1]][4]
    
    
    lat = lat_deg[1] + lat_deg[2]/60 + lat_deg[3]/3600
    long = long_deg[1] + long_deg[2]/60 + long_deg[3]/3600
    
    if(hemisfere %in% c("S","s")){lat = 0-lat}
    if(side %in% c("W","w")){long = 0-long}
    
  }
  
  coordinates = data.frame(lat = lat, long = long)
  coordinates = st_as_sf(coordinates, coords = c("long", "lat"), crs = "+proj=longlat +datum=WGS84")
  coordinates = st_transform(coordinates, 4087)
  
  return(st_coordinates(coordinates))

}






### angtoxyz (Converts angular directions + distance to XYZ coordinates)

angtoxyz = function(raw){
  
  raw = raw[raw$distance != 0,]
  
  raw$v_ang_rad = (pi/180)* raw$inclination
  raw$z = sin(raw$v_ang_rad)*raw$distance # Z Values!!
  raw$yx_hipot = cos(raw$v_ang_rad)*raw$distance
  
  raw$bearing_rad = (pi/180)*raw$bearing
  raw$y =  cos(raw$bearing_rad)*raw$yx_hipot #y values is here
  raw$x = sin(raw$bearing_rad)*raw$yx_hipot #x values is here
  
  raw = raw[,c("y","x","z","time")]
  
  
  return(raw)
}









# angtopic (Converts angular directions + distance to a picture). But is horribly slow!! be warned!!

angtopic = function(raw){ 
  
  #raw$bearing = round(raw$bearing, 4)
  raw = raw[,-4]
  raw = raw[raw$bearing <= 360,]
  raw$bearing[raw$bearing == 360] = 0
  raw$distance[raw$distance == 0] = NA

  
  #turn the angular data into pixels of a raster (1/2 degree pixels)
  
  pixel_angular_resolution = 0.5
  horz_axis = seq(0, 360, by = pixel_angular_resolution)
  vert_axis = seq(min(raw$inclination), max(raw$inclination), by = pixel_angular_resolution)
  
  
  pic = expand.grid(horz_axis,vert_axis)
  colnames(pic) = c("horz_axis", "vert_axis")
  pic$distance = NA
  

  for(h in unique(pic$horz_axis)){
    for(v in unique(pic$vert_axis)){
      
      t = raw %>% 
        filter(bearing >= h, #bearing bigger than the angle of the pixel  start
               bearing < h+pixel_angular_resolution, #bearing smaller than the angle of the pixel end
               inclination >= v,  #inclination bigger than the angle of the pixel  start
               inclination < v+pixel_angular_resolution) #inclination smaller than the angle of the pixel end
      
      pic$distance[pic$horz_axis == h & pic$vert_axis == v] = mean(t$distance, na.rm = T) #save the average distance

    }
    print(paste(h, "/360", sep = ""))

  }
  

  Figure = pic %>% 
    ggplot(aes(x=horz_axis, y=vert_axis, fill=distance)) + 
    geom_raster() +
    scale_fill_viridis_c(option = "magma", direction = 1, na.value="white")

  return(Figure)
  
}
  
  
  
    
  
  
  
  
  








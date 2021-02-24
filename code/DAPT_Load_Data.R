#Merritt Harlan
#Discharge from Pressure Transducers (DAPT)

#Load Input Data and Create bamdata files
library(sp)
library(raster)

#Functions------
#extract elevation profiles
ortho_lines <- function(DEM, centerline_coords, scale_xs){
  ortho_lines = list()
  UTMcrs = crs(DEM)
  cl_coords_df = data.frame(x=c(1:length(centerline_coords[,1])), lon=centerline_coords[,1], 
                            lat=centerline_coords[,2])
  mean_lon = cl_coords_df$lon[1:(length(cl_coords_df$lon)-1)]+(1/2)*diff(cl_coords_df$lon)
  mean_lat = cl_coords_df$lat[1:(length(cl_coords_df$lat)-1)]+(1/2)*diff(cl_coords_df$lat)
  slope = (-1)/(diff(cl_coords_df$lat)/diff(cl_coords_df$lon))
  for (i in which(!is.na(slope))){
    left_lon = mean_lon-(cos(atan(slope))*scale_xs)
    left_lat = mean_lat-(sin(atan(slope))*scale_xs)
    right_lon = mean_lon+(cos(atan(slope))*scale_xs)
    right_lat = mean_lat+(sin(atan(slope))*scale_xs)
    ortho_lines[[i]] = SpatialLines(list(Lines(Line(cbind(c(left_lon[i],right_lon[i]),
                                                          c(left_lat[i], right_lat[i]))), ID=i)),
                                    proj4string = UTMcrs)
  }
  ortho_lines
}

setwd( "C:/Users/Merritt/Box Sync/Research (mharlan@umass.edu)/DAPT/input")
#Load Data
Sac_rch1_PT_coords <- read.csv("Pressure Transducer Coordinates/Sac_rch1_PT_coords.csv")
Sac_rch2_PT_coords <- read.csv("Pressure Transducer Coordinates/Sac_rch2_PT_coords.csv")
Sac_rch3_PT_coords <- read.csv("Pressure Transducer Coordinates/Sac_rch3_PT_coords.csv")
Olen_PT_coords <- read.csv("Pressure Transducer Coordinates/Olen_PT_coords.csv")
Tan_PT_coords <- read.csv("Pressure Transducer Coordinates/Tan_PT_coords.csv")
Will_PT_coords <- read.csv("Pressure Transducer Coordinates/Will_PT_coords.csv")
North_Sask_PT_coords <- read.csv("Pressure Transducer Coordinates/North_Sask_PT_coords.csv")

Sac_rch1_DEM <- raster("Digital Elevation Model (DEM)/Sac_DEM_rch1.tif")
Sac_rch2_DEM <- raster("Digital Elevation Model (DEM)/Sac_DEM_rch2.tif")
Sac_rch3_DEM <- raster("Digital Elevation Model (DEM)/Sac_DEM_rch3.tif")

Sac_rch1_cl <-readOGR(dsn = "River Centerline", layer = "Sac_rch1_UTM")
Sac_rch2_cl <-readOGR(dsn = "River Centerline", layer = "Sac_rch2_UTM")
Sac_rch3_cl <-readOGR(dsn = "River Centerline", layer = "Sac_rch3_cl_UTM")

Sac_orthos_rch1 <- ortho_lines(Sac_rch1_DEM, data.frame(Sac_rch1_cl@lines[[1]]@Lines[[1]]@coords), 1000)
Sac_orthos_rch2 <- ortho_lines(Sac_rch2_DEM, data.frame(Sac_rch2_cl@lines[[1]]@Lines[[1]]@coords), 1000)
Sac_orthos_rch3 <- ortho_lines(Sac_rch3_DEM, data.frame(Sac_rch3_cl@lines[[1]]@Lines[[1]]@coords), 1000)

#Find closest orthos to pressure transducers
Sac_rch1_PT_coords_sp <-SpatialPoints(Sac_rch1_PT_coords[,4:3], proj4string = crs(Sac_DEM_rch1))
Sac_rch1_dist_mat <- pointDistance(Sac_rch1_PT_coords_sp, data.frame(Sac_cl_rch1@lines[[1]]@Lines[[1]]@coords))
Sac_rch1_closest_ortho <- Sac_orthos_rch1[apply(Sac_rch1_dist_mat, 1, FUN = which.min)]

Sac_rch2_PT_coords_sp <-SpatialPoints(Sac_rch2_PT_coords[,4:3], proj4string = crs(Sac_DEM_rch2))
Sac_rch2_dist_mat <- pointDistance(Sac_rch2_PT_coords_sp, data.frame(Sac_cl_rch2@lines[[1]]@Lines[[1]]@coords))
Sac_rch2_closest_ortho <- Sac_orthos_rch2[apply(Sac_rch2_dist_mat, 1, FUN = which.min)]

Sac_rch3_PT_coords_sp <-SpatialPoints(Sac_rch3_PT_coords[,4:3], proj4string = crs(Sac_DEM_rch3))
Sac_rch3_dist_mat <- pointDistance(Sac_rch3_PT_coords_sp, data.frame(Sac_cl_rch3@lines[[1]]@Lines[[1]]@coords))
Sac_rch3_closest_ortho <- Sac_orthos_rch3[apply(Sac_rch3_dist_mat, 1, FUN = which.min)]

#copy cross section profile code at specific location
numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Sac_rch1_elev_prof <- foreach(i=1:length(Sac_rch1_closest_ortho), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Sac_DEM_rch1, Sac_rch1_closest_ortho[[i]], method = 'simple')
                              }
Sac_rch2_elev_prof <- foreach(i=1:length(Sac_rch2_closest_ortho), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Sac_DEM_rch2, Sac_rch2_closest_ortho[[i]], method = 'simple')
                              }
Sac_rch3_elev_prof <- foreach(i=1:length(Sac_rch3_closest_ortho), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Sac_DEM_rch3, Sac_rch3_closest_ortho[[i]], method = 'simple')
                              }
proc.time()-ptm
stopCluster(cl)

Sac_rch1_scale_x <- vector()
for(i in 1:length(Sac_rch1_elev_prof)){
  Sac_rch1_scale_x[i] = 2000/length(Sac_rch1_elev_prof[[i]])
}

Sac_rch2_scale_x <- vector()
for(i in 1:length(Sac_rch2_elev_prof)){
  Sac_rch2_scale_x[i] = 2000/length(Sac_rch2_elev_prof[[i]])
}

Sac_rch3_scale_x <- vector()
for(i in 1:length(Sac_rch3_elev_prof)){
  Sac_rch3_scale_x[i] = 2000/length(Sac_rch3_elev_prof[[i]])
}

###Olentangy-----
setwd(paste(home_dir,"/Olentangy", sep =""))
Olen_PT_coords <- read.csv("Olentangy_PT_coord.csv")

Olen_bounding_box <- c(min(Olen_PT_coords$Latitude)-2000, max(Olen_PT_coords$Latitude)+2000,
                       min(Olen_PT_coords$Longitude)-1000, max(Olen_PT_coords$Longitude)+1000)

setwd(paste(home_dir,"/Olentangy/DEM", sep =""))
Olen_DEM <- raster("Cropped_Olen_NED.tif")

setwd(paste(home_dir,"/Olentangy/GRWL", sep =""))
Olen_cl <- read.csv("Olen_GRWL.csv")
Olen_cl_sp <- SpatialPointsDataFrame(Olen_cl[,1:2], Olen_cl, proj4string = crs(Olen_DEM))
Olen_cl <- crop(Olen_cl_sp, Olen_bounding_box)

setwd(paste(home_dir,"/Olentangy", sep =""))
Olen_cl <-readOGR(dsn = ".", layer = "Olen_cl_UTM")

ptm<-proc.time()
Olen_orthos <- ortho_lines(Olen_DEM, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords), 500)
proc.time()-ptm


#Find closest orthos to pressure transducers
Olen_PT_coords_sp <-SpatialPoints(Olen_PT_coords[,3:4], proj4string = crs(Olen_DEM))
Olen_dist_mat <- pointDistance(Olen_PT_coords_sp, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords))
Olen_closest_ortho <- Olen_orthos[apply(Olen_dist_mat, 1, FUN = which.min)]

#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Olen_elev_prof_grwl_NED <- foreach(i=1:length(Olen_closest_ortho), 
                                   .combine='rbind', .packages=c('raster')) %dopar% {
                                     extract(Olen_DEM, Olen_closest_ortho[[i]], method = 'simple')
                                   }
proc.time()-ptm
stopCluster(cl)

Olen_scale_x <- vector()
for(i in 1:length(Olen_elev_prof_grwl_NED)){
  Olen_scale_x[i] = 1000/length(Olen_elev_prof_grwl_NED[[i]])
}

plot(Olen_DEM)
plot(Olen_PT_coords_sp, add = TRUE)
for(i in 1:length(Olen_closest_ortho)){
  plot(Olen_closest_ortho[[i]], add = TRUE, col = "blue")
}
###Willamette----
setwd(paste(home_dir,"/Willamette", sep =""))
Will_PT_coords <- read.csv("Willamette_PT_coords.csv")
Will_bounding_box <- c(min(Will_PT_coords$Longitude)-2000, max(Will_PT_coords$Longitude)+2000,
                       min(Will_PT_coords$Latitude)-1000, max(Will_PT_coords$Latitude)+1000)

setwd(paste(home_dir,"/Willamette/DEM", sep =""))
Will_DEM <- raster("Cropped_Will_NED.tif")

setwd(paste(home_dir,"/Willamette", sep =""))
Will_cl <-readOGR(dsn = ".", layer = "Will_cl_UTM")

ptm<-proc.time()
Will_orthos <- ortho_lines(Will_DEM, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords), 1000)
proc.time()-ptm

#Find closest orthos to pressure transducers
Will_PT_coords_sp <-SpatialPoints(Will_PT_coords[,3:4], proj4string = crs(Will_DEM))
Will_dist_mat <- pointDistance(Will_PT_coords_sp, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords))
Will_closest_ortho <- Will_orthos[apply(Will_dist_mat, 1, FUN = which.min)]

#copy cross section profile code at specific location
plot(Will_PT_coords_sp)
for(i in 1:length(Will_closest_ortho)){
  plot(Will_closest_ortho[[i]], add = TRUE, col = "blue")
}
numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Will_elev_prof_grwl_NED <- foreach(i=1:length(Will_closest_ortho), 
                                   .combine='rbind', .packages=c('raster')) %dopar% {
                                     extract(Will_DEM, Will_closest_ortho[[i]], method = 'simple')
                                   }
proc.time()-ptm
stopCluster(cl)

Will_scale_x <- vector()
for(i in 1:length(Will_elev_prof_grwl_NED)){
  Will_scale_x[i] = 2000/length(Will_elev_prof_grwl_NED[[i]])
}



###Tanana------
setwd(paste(home_dir,"/Tanana", sep =""))
Tan_PT_coords <- read.csv("Tanana_PT_coord.csv")

Tan_bounding_box <- c(min(Tan_PT_coords$Longitude)-10000, max(Tan_PT_coords$Longitude)+10000,
                      min(Tan_PT_coords$Latitude)-10000, max(Tan_PT_coords$Latitude)+10000)

ptm<-proc.time()
setwd(paste(home_dir,"/Tanana/DEM", sep =""))
Tan_DEM_raw <- raster("Tan_DEM_2m.tif")
Tan_DEM <- crop(Tan_DEM_raw, Tan_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/Tanana", sep =""))
Tan_cl <-readOGR(dsn = ".", layer = "Tan_cl")

ptm<-proc.time()
Tan_orthos <- ortho_lines(Tan_DEM, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords), 5000)
proc.time()-ptm

#Find closest orthos to pressure transducers
Tan_PT_coords_sp <-SpatialPoints(Tan_PT_coords[,3:4], proj4string = crs(Tan_DEM))
Tan_dist_mat <- pointDistance(Tan_PT_coords_sp, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords))
Tan_closest_ortho <- Tan_orthos[apply(Tan_dist_mat, 1, FUN = which.min)]

plot(Tan_PT_coords_sp)
for(i in 1:length(Tan_closest_ortho)){
  plot(Tan_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Tan_elev_prof_grwl_NED <- foreach(i=1:length(Tan_closest_ortho), 
                                  .combine='rbind', .packages=c('raster')) %dopar% {
                                    extract(Tan_DEM, Tan_closest_ortho[[i]], method = 'simple')
                                  }
proc.time()-ptm
stopCluster(cl)

Tan_scale_x <- vector()
for(i in 1:length(Tan_elev_prof_grwl_NED)){
  Tan_scale_x[i] = 10000/length(Tan_elev_prof_grwl_NED[[i]])
}

###North Sask----
setwd(paste(home_dir,"/North_Saskatchewan_2017", sep =""))
North_Sask_PT_coords <- read.csv("North_Sask_PT_coords.csv")
North_Sask_PT_coords_sp <-SpatialPoints(North_Sask_PT_coords[,4:3], proj4string = CRS("+init=epsg:4326"))
North_Sask_PT_coords_sp <- spTransform(North_Sask_PT_coords_sp, )

North_Sask_bounding_box <- c(min(Tan_PT_coords$Longitude)-0.1, max(Tan_PT_coords$Longitude)+0.1,
                             min(Tan_PT_coords$Latitude)-0.1, max(Tan_PT_coords$Latitude)+0.1)

ptm<-proc.time()
setwd(paste(home_dir,"/North_Saskatchewan_2017/DEM", sep =""))
North_Sask_DEM_raw <- raster("DEM.tif")
North_Sask_DEM <- crop(North_Sask_DEM_raw, North_Sask_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/North_Saskatchewan_2017", sep =""))
North_Sask_cl <-readOGR(dsn = ".", layer = "North_Sask_cl")

ptm<-proc.time()
North_Sask_orthos <- ortho_lines(North_Sask_DEM, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords), 1000)
proc.time()-ptm

#Find closest orthos to pressure transducers
North_Sask_PT_coords_sp <-SpatialPoints(North_Sask_PT_coords[,4:3], proj4string = crs(North_Sask_DEM))
North_Sask_dist_mat <- pointDistance(North_Sask_PT_coords_sp, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords))
North_Sask_closest_ortho <- North_Sask_orthos[apply(North_Sask_dist_mat, 1, FUN = which.min)]


#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
North_Sask_elev_prof_grwl_NED <- foreach(i=1:length(North_Sask_closest_ortho), 
                                         .combine='rbind', .packages=c('raster')) %dopar% {
                                           extract(North_Sask_DEM, North_Sask_closest_ortho[[i]], method = 'simple')
                                         }
proc.time()-ptm
stopCluster(cl)

North_Sask_scale_x <- 15

plot(North_Sask_PT_coords_sp)
for(i in 1:length(North_Sask_closest_ortho)){
  plot(North_Sask_closest_ortho[[i]], add = TRUE, col = "blue")
}

save(Sac_elev_prof_grwl_3m_rch1, "Sac_elev_prof_grwl_3m_rch1.RData")
save(Sac_elev_prof_grwl_3m_rch2, "Sac_elev_prof_grwl_3m_rch2.RData")
save(Sac_elev_prof_grwl_3m_rch3, "Sac_elev_prof_grwl_3m_rch3.RData")
save(Olen_elev_prof_grwl_NED, "Olen_elev_prof_grwl.RData")
save(Will_elev_prof_grwl_NED, "Will_elev_prof_grwl.RData")
save(Tan_elev_prof_grwl_NED, "Tan_elev_prof_grwl.RData")
save(North_Sask_elev_prof_grwl_NED, "North_Sask_elev_prof_grwl.RData")

###PAD 2018------
#Athabasca------
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
Athabasca_PT_coords <- read.csv("Athabasca_PT_coords.csv")

Athabasca_bounding_box <- c(min(Athabasca_PT_coords$Longitude)-2000, max(Athabasca_PT_coords$Longitude)+2000,
                            min(Athabasca_PT_coords$Latitude)-1000, max(Athabasca_PT_coords$Latitude)+1000)

ptm<-proc.time()
Athabasca_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Athabasca_DEM <- crop(Athabasca_DEM_raw, Athabasca_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/Athabasca", sep =""))
Athabasca_cl <-readOGR(dsn = ".", layer = "Athabasca_cl_UTM")

setwd(paste(home_dir,"/PAD/PAD_2018", sep =""))
PAD_2018_ADCP_UTM <-readOGR(dsn = ".", layer = "PAD_2018_ADCP_UTM")

Athabasca_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Athabasca_bounding_box))

ptm<-proc.time()
Athabasca_orthos <- ortho_lines(DEM = Athabasca_DEM, centerline_coords = 
                                  data.frame(Athabasca_cl@lines[[1]]@Lines[[1]]@coords),
                                scale_xs = 300)
proc.time()-ptm
Athabasca_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Athabasca_bounding_box))
Athabasca_qobs_df$Location <- unlist(lapply(Athabasca_qobs_df$Location, toString))
Athabasca_qobs_df <- Athabasca_qobs_df[startsWith(Athabasca_qobs_df$Location,
                                                  "ATH"),]
Athabasca_qobs_df <- Athabasca_qobs_df %>%
  mutate(day = as.Date(new_date, format="%Y-%m-%d"))
Athabasca_corr_Q <- Athabasca_qobs_df[Athabasca_qobs_df$Q>600,]
write.csv(Athabasca_corr_Q, "Athabasca_corr_Q.csv")
Athabasca_qobs_sp <-SpatialPoints(Athabasca_corr_Q[,8:9],proj4string = crs(Athabasca_DEM))
Athabasca_dist_mat_ADCP <- pointDistance(Athabasca_qobs_sp, data.frame(Athabasca_cl@lines[[1]]@Lines[[1]]@coords))
Athabasca_closest_ADCP <- Athabasca_orthos[apply(Athabasca_dist_mat_ADCP, 1, FUN = which.min)]
Athabasca_merged_ADCP <- do.call(rbind, Athabasca_closest_ADCP)
Athabasca_ADCP_df <- SpatialLinesDataFrame(Athabasca_merged_ADCP, Athabasca_corr_Q, match.ID = FALSE)
writeOGR(Athabasca_ADCP_df, dsn = ".", layer = "Athabasca_ADCP", driver = "ESRI Shapefile")

#Find closest orthos to pressure transducers
Athabasca_PT_coords_sp <-SpatialPoints(Athabasca_PT_coords[,4:5], proj4string = crs(Athabasca_DEM))
Athabasca_dist_mat <- pointDistance(Athabasca_PT_coords_sp, data.frame(Athabasca_cl@lines[[1]]@Lines[[1]]@coords))
Athabasca_closest_ortho <- Athabasca_orthos[apply(Athabasca_dist_mat, 1, FUN = which.min)]

plot(Athabasca_PT_coords_sp)
for(i in 1:length(Athabasca_closest_ortho)){
  plot(Athabasca_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Athabasca_elev_prof <- foreach(i=1:length(Athabasca_closest_ortho), 
                               .combine='rbind', .packages=c('raster')) %dopar% {
                                 extract(Athabasca_DEM, Athabasca_closest_ortho[[i]], method = 'simple')
                               }
proc.time()-ptm
stopCluster(cl)

Athabasca_scale_x <- vector()
for(i in 1:length(Athabasca_elev_prof)){
  Athabasca_scale_x[i] = 12000/length(Athabasca_elev_prof[[i]])
}

#Coupe-----
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
Coupe_PT_coords <- read.csv("Coupe_PT_coords.csv")

Coupe_bounding_box <- c(min(Coupe_PT_coords$Longitude)-2000, max(Coupe_PT_coords$Longitude)+2000,
                        min(Coupe_PT_coords$Latitude)-1000, max(Coupe_PT_coords$Latitude)+1000)

ptm<-proc.time()
Coupe_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Coupe_DEM <- crop(Coupe_DEM_raw, Coupe_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/Coupe", sep =""))
Coupe_cl <-readOGR(dsn = ".", layer = "Coupe_CL_UTM")

Coupe_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Coupe_bounding_box))

ptm<-proc.time()
Coupe_orthos <- ortho_lines(DEM = Coupe_DEM, centerline_coords = 
                              data.frame(Coupe_cl@lines[[1]]@Lines[[1]]@coords),
                            scale_xs = 200)
proc.time()-ptm
Coupe_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Coupe_bounding_box))
Coupe_qobs_df$Location <- unlist(lapply(Coupe_qobs_df$Location, toString))
Coupe_qobs_df <- Coupe_qobs_df[startsWith(Coupe_qobs_df$Location,
                                          "CP"),]
Coupe_qobs_df <- Coupe_qobs_df %>%
  mutate(day = as.Date(new_date, format="%Y-%m-%d"))
Coupe_corr_Q <- Coupe_qobs_df[Coupe_qobs_df$Q>50,]
write.csv(Coupe_corr_Q, "Coupe_corr_Q.csv")
Coupe_qobs_sp <-SpatialPoints(Coupe_corr_Q[,8:9],proj4string = crs(Coupe_DEM))
Coupe_dist_mat_ADCP <- pointDistance(Coupe_qobs_sp, data.frame(Coupe_cl@lines[[1]]@Lines[[1]]@coords))
Coupe_closest_ADCP <- Coupe_orthos[apply(Coupe_dist_mat_ADCP, 1, FUN = which.min)]
Coupe_merged_ADCP <- do.call(rbind, Coupe_closest_ADCP)
Coupe_ADCP_df <- SpatialLinesDataFrame(Coupe_merged_ADCP, Coupe_corr_Q, match.ID = FALSE)
writeOGR(Coupe_ADCP_df, dsn = ".", layer = "Coupe_ADCP", driver = "ESRI Shapefile")

#Find closest orthos to pressure transducers
Coupe_PT_coords_sp <-SpatialPoints(Coupe_PT_coords[,4:5], proj4string = crs(Coupe_DEM))
Coupe_dist_mat <- pointDistance(Coupe_PT_coords_sp, data.frame(Coupe_cl@lines[[1]]@Lines[[1]]@coords))
Coupe_closest_ortho <- Coupe_orthos[apply(Coupe_dist_mat, 1, FUN = which.min)]

plot(Coupe_PT_coords_sp)
for(i in 1:length(Coupe_closest_ortho)){
  plot(Coupe_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Coupe_elev_prof <- foreach(i=1:length(Coupe_closest_ortho), 
                           .combine='rbind', .packages=c('raster')) %dopar% {
                             extract(Coupe_DEM, Coupe_closest_ortho[[i]], method = 'simple')
                           }
proc.time()-ptm
stopCluster(cl)

Coupe_scale_x <- vector()
for(i in 1:length(Coupe_elev_prof)){
  Coupe_scale_x[i] = 10000/length(Coupe_elev_prof[[i]])
}

#Fletcher-----
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
Fletcher_PT_coords <- read.csv("Fletcher_PT_coords.csv")

Fletcher_bounding_box <- c(min(Fletcher_PT_coords$Longitude)-2000, max(Fletcher_PT_coords$Longitude)+2000,
                           min(Fletcher_PT_coords$Latitude)-1000, max(Fletcher_PT_coords$Latitude)+1000)

ptm<-proc.time()
Fletcher_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Fletcher_DEM <- crop(Fletcher_DEM_raw, Fletcher_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/Fletcher", sep =""))
Fletcher_cl <-readOGR(dsn = ".", layer = "Fletcher_CL_UTM")
Fletcher_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Fletcher_bounding_box))
ptm<-proc.time()
Fletcher_orthos <- ortho_lines(DEM = Fletcher_DEM, centerline_coords = 
                                 data.frame(Fletcher_cl@lines[[1]]@Lines[[1]]@coords),
                               scale_xs = 200)
proc.time()-ptm
Fletcher_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Fletcher_bounding_box))
Fletcher_qobs_df$Location <- unlist(lapply(Fletcher_qobs_df$Location, toString))
Fletcher_qobs_df <- Fletcher_qobs_df[startsWith(Fletcher_qobs_df$Location,
                                                "FL"),]
Fletcher_qobs_df <- Fletcher_qobs_df %>%
  mutate(day = as.Date(new_date, format="%Y-%m-%d"))
Fletcher_corr_Q <- Fletcher_qobs_df[Fletcher_qobs_df$Q >100,]
write.csv(Fletcher_corr_Q, "Fletcher_corr_Q.csv")
Fletcher_qobs_sp <-SpatialPoints(Fletcher_corr_Q[,8:9],proj4string = crs(Fletcher_DEM))
Fletcher_dist_mat_ADCP <- pointDistance(Fletcher_qobs_sp, data.frame(Fletcher_cl@lines[[1]]@Lines[[1]]@coords))
Fletcher_closest_ADCP <- Fletcher_orthos[apply(Fletcher_dist_mat_ADCP, 1, FUN = which.min)]
Fletcher_merged_ADCP <- do.call(rbind, Fletcher_closest_ADCP)
Fletcher_ADCP_df <- SpatialLinesDataFrame(Fletcher_merged_ADCP, Fletcher_corr_Q, match.ID = FALSE)
writeOGR(Fletcher_ADCP_df, dsn = ".", layer = "Fletcher_ADCP", driver = "ESRI Shapefile")

#Find closest orthos to pressure transducers
Fletcher_PT_coords_sp <-SpatialPoints(Fletcher_PT_coords[,4:5], proj4string = crs(Fletcher_DEM))
Fletcher_dist_mat <- pointDistance(Fletcher_PT_coords_sp, data.frame(Fletcher_cl@lines[[1]]@Lines[[1]]@coords))
Fletcher_closest_ortho <- Fletcher_orthos[apply(Fletcher_dist_mat, 1, FUN = which.min)]

plot(Fletcher_DEM)
plot(Fletcher_PT_coords_sp, add = TRUE)
for(i in 1:length(Fletcher_closest_ortho)){
  plot(Fletcher_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Fletcher_elev_prof <- foreach(i=1:length(Fletcher_closest_ortho), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Fletcher_DEM, Fletcher_closest_ortho[[i]], method = 'simple')
                              }
proc.time()-ptm
stopCluster(cl)

Fletcher_scale_x <- vector()
for(i in 1:length(Fletcher_elev_prof)){
  Fletcher_scale_x[i] = 16000/length(Fletcher_elev_prof[[i]])
}

#Mamawi-----
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
Mamawi_PT_coords <- read.csv("Mamawi_PT_coords.csv")

Mamawi_bounding_box <- c(min(Mamawi_PT_coords$Longitude)-2000, max(Mamawi_PT_coords$Longitude)+2000,
                         min(Mamawi_PT_coords$Latitude)-1000, max(Mamawi_PT_coords$Latitude)+1000)

ptm<-proc.time()
Mamawi_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Mamawi_DEM <- crop(Mamawi_DEM_raw, Mamawi_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/Mamawi", sep =""))
Mamawi_cl <-readOGR(dsn = ".", layer = "Mamawi_CL_UTM")

ptm<-proc.time()
Mamawi_orthos <- ortho_lines(DEM = Mamawi_DEM, centerline_coords = 
                               data.frame(Mamawi_cl@lines[[1]]@Lines[[1]]@coords),
                             scale_xs = 200)
proc.time()-ptm
Mamawi_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Mamawi_bounding_box))
Mamawi_qobs_df$Location <- unlist(lapply(Mamawi_qobs_df$Location, toString))
Mamawi_qobs_df <- Mamawi_qobs_df[startsWith(Mamawi_qobs_df$Location,
                                            "MC"),]
Mamawi_qobs_df <- Mamawi_qobs_df %>%
  mutate(day = as.Date(new_date, format="%Y-%m-%d"))
Mamawi_corr_Q <- Mamawi_qobs_df
write.csv(Mamawi_corr_Q, "Mamawi_corr_Q.csv")
Mamawi_qobs_sp <-SpatialPoints(Mamawi_corr_Q[,8:9],proj4string = crs(Mamawi_DEM))
Mamawi_dist_mat_ADCP <- pointDistance(Mamawi_qobs_sp, data.frame(Mamawi_cl@lines[[1]]@Lines[[1]]@coords))
Mamawi_closest_ADCP <- Mamawi_orthos[apply(Mamawi_dist_mat_ADCP, 1, FUN = which.min)]
Mamawi_merged_ADCP <- do.call(rbind, Mamawi_closest_ADCP)
Mamawi_ADCP_df <- SpatialLinesDataFrame(Mamawi_merged_ADCP, Mamawi_corr_Q, match.ID = FALSE)
writeOGR(Mamawi_ADCP_df, dsn = ".", layer = "Mamawi_ADCP", driver = "ESRI Shapefile")


#Find closest orthos to pressure transducers
Mamawi_PT_coords_sp <-SpatialPoints(Mamawi_PT_coords[,4:5], proj4string = crs(Mamawi_DEM))
Mamawi_dist_mat <- pointDistance(Mamawi_PT_coords_sp, data.frame(Mamawi_cl@lines[[1]]@Lines[[1]]@coords))
Mamawi_closest_ortho <- Mamawi_orthos[apply(Mamawi_dist_mat, 1, FUN = which.min)]

plot(Mamawi_DEM)
plot(Mamawi_PT_coords_sp, add = TRUE)
for(i in 1:length(Mamawi_closest_ortho)){
  plot(Mamawi_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Mamawi_elev_prof <- foreach(i=1:length(Mamawi_closest_ortho), 
                            .combine='rbind', .packages=c('raster')) %dopar% {
                              extract(Mamawi_DEM, Mamawi_closest_ortho[[i]], method = 'simple')
                            }
proc.time()-ptm
stopCluster(cl)

Mamawi_scale_x <- vector()
for(i in 1:length(Mamawi_elev_prof)){
  Mamawi_scale_x[i] = 2000/length(Mamawi_elev_prof[[i]])
}

#Rochers-----
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
Rochers_PT_coords <- read.csv("Rochers_PT_coords.csv")

Rochers_bounding_box <- c(min(Rochers_PT_coords$Longitude)-2000, max(Rochers_PT_coords$Longitude)+2000,
                          min(Rochers_PT_coords$Latitude)-1000, max(Rochers_PT_coords$Latitude)+1000)

ptm<-proc.time()
Rochers_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Rochers_DEM <- crop(Rochers_DEM_raw, Rochers_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/Rochers", sep =""))
Rochers_cl <-readOGR(dsn = ".", layer = "Rochers_CL_UTM")
Rochers_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, Rochers_bounding_box))
Rochers_corr_Q <- Rochers_qobs_df[Rochers_qobs_df$Q > 1400,]
write.csv(Rochers_corr_Q, "Rochers_corr_Q.csv")
ptm<-proc.time()
Rochers_orthos <- ortho_lines(DEM = Rochers_DEM, centerline_coords = 
                                data.frame(Rochers_cl@lines[[1]]@Lines[[1]]@coords),
                              scale_xs = 200)
proc.time()-ptm
Rochers_qobs_sp <-SpatialPoints(Rochers_corr_Q[,8:9],proj4string = crs(Rochers_DEM))
Rochers_dist_mat_ADCP <- pointDistance(Rochers_qobs_sp, data.frame(Rochers_cl@lines[[1]]@Lines[[1]]@coords))
Rochers_closest_ADCP <- Rochers_orthos[apply(Rochers_dist_mat_ADCP, 1, FUN = which.min)]
Rochers_merged_ADCP <- do.call(rbind, Rochers_closest_ADCP)
Rochers_ADCP_df <- SpatialLinesDataFrame(Rochers_merged_ADCP, Rochers_corr_Q, match.ID = FALSE)
writeOGR(Rochers_ADCP_df, dsn = ".", layer = "Rochers_ADCP", driver = "ESRI Shapefile")


#Find closest orthos to pressure transducers
Rochers_PT_coords_sp <-SpatialPoints(Rochers_PT_coords[,4:5], proj4string = crs(Rochers_DEM))
Rochers_dist_mat <- pointDistance(Rochers_PT_coords_sp, data.frame(Rochers_cl@lines[[1]]@Lines[[1]]@coords))
Rochers_closest_ortho <- Rochers_orthos[apply(Rochers_dist_mat, 1, FUN = which.min)]

plot(Rochers_DEM)
plot(Rochers_PT_coords_sp, add = TRUE)
for(i in 1:length(Rochers_closest_ortho)){
  plot(Rochers_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Rochers_elev_prof <- foreach(i=1:length(Rochers_closest_ortho), 
                             .combine='rbind', .packages=c('raster')) %dopar% {
                               extract(Rochers_DEM, Rochers_closest_ortho[[i]], method = 'simple')
                             }
proc.time()-ptm
stopCluster(cl)

Rochers_scale_x <- vector()
for(i in 1:length(Rochers_elev_prof)){
  Rochers_scale_x[i] = 2000/length(Rochers_elev_prof[[i]])
}

#QF-----
setwd(paste(home_dir,"/PAD/PAD_2018/PAD_2018_PT", sep =""))
QF_PT_coords <- read.csv("QF_PT_coords.csv")

QF_bounding_box <- c(min(QF_PT_coords$Longitude)-2000, max(QF_PT_coords$Longitude)+2000,
                     min(QF_PT_coords$Latitude)-1000, max(QF_PT_coords$Latitude)+1000)

ptm<-proc.time()
QF_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
QF_DEM <- crop(QF_DEM_raw, QF_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2018/QF", sep =""))
QF_cl <-readOGR(dsn = ".", layer = "QF_CL_UTM")
QF_qobs_df <- data.frame(crop(PAD_2018_ADCP_UTM, QF_bounding_box))
ptm<-proc.time()
QF_orthos <- ortho_lines(DEM = QF_DEM, centerline_coords = 
                           data.frame(QF_cl@lines[[1]]@Lines[[1]]@coords),
                         scale_xs = 200)
proc.time()-ptm
QF_qobs_sp <-SpatialPoints(QF_corr_Q[,8:9],proj4string = crs(QF_DEM))
QF_dist_mat_ADCP <- pointDistance(QF_qobs_sp, data.frame(QF_cl@lines[[1]]@Lines[[1]]@coords))
QF_closest_ADCP <- QF_orthos[apply(QF_dist_mat_ADCP, 1, FUN = which.min)]
QF_merged_ADCP <- do.call(rbind, QF_closest_ADCP)
QF_ADCP_df <- SpatialLinesDataFrame(QF_merged_ADCP, QF_corr_Q, match.ID = FALSE)
writeOGR(QF_ADCP_df, dsn = ".", layer = "QF_ADCP", driver = "ESRI Shapefile")


#Find closest orthos to pressure transducers
QF_PT_coords_sp <-SpatialPoints(QF_PT_coords[,4:5], proj4string = crs(QF_DEM))
QF_dist_mat <- pointDistance(QF_PT_coords_sp, data.frame(QF_cl@lines[[1]]@Lines[[1]]@coords))
QF_closest_ortho <- QF_orthos[apply(QF_dist_mat, 1, FUN = which.min)]

plot(QF_DEM)
plot(QF_PT_coords_sp, add = TRUE)
for(i in 1:length(QF_closest_ortho)){
  plot(QF_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
QF_elev_prof <- foreach(i=1:length(QF_closest_ortho), 
                        .combine='rbind', .packages=c('raster')) %dopar% {
                          extract(QF_DEM, QF_closest_ortho[[i]], method = 'simple')
                        }
proc.time()-ptm
stopCluster(cl)

QF_scale_x <- vector()
for(i in 1:length(QF_elev_prof)){
  QF_scale_x[i] = 2000/length(QF_elev_prof[[i]])
}

#Embarras-----
setwd(paste(home_dir,"/PAD/PAD_2019/PAD_2019_PT", sep =""))
Embarras_PT_coords <- read.csv("Embarras_PT_coords.csv")

Embarras_bounding_box <- c(min(Embarras_PT_coords$Longitude)-2000, max(Embarras_PT_coords$Longitude)+2000,
                           min(Embarras_PT_coords$Latitude)-1000, max(Embarras_PT_coords$Latitude)+1000)

ptm<-proc.time()
Embarras_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Embarras_DEM <- crop(Embarras_DEM_raw, Embarras_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2019/Embarras", sep =""))
Embarras_cl <-readOGR(dsn = ".", layer = "Embarras_CL_UTM")

setwd(paste(home_dir,"/PAD/PAD_2019", sep =""))
PAD_2019_ADCP_UTM <-readOGR(dsn = ".", layer = "PAD_2019_ADCP_UTM")

Embarras_qobs_df <- data.frame(crop(PAD_2019_ADCP_UTM, Embarras_bounding_box))
Embarras_corr_Q <- Embarras_qobs_df[Embarras_qobs_df$Q > 100 & Embarras_qobs_df$Q < 200, ]
write.csv(Embarras_corr_Q, "Embarras_corr_Q.csv")
ptm<-proc.time()
Embarras_orthos <- ortho_lines(DEM = Embarras_DEM, centerline_coords = 
                                 data.frame(Embarras_cl@lines[[1]]@Lines[[1]]@coords),
                               scale_xs = 100)
proc.time()-ptm
Embarras_qobs_sp <-SpatialPoints(Embarras_corr_Q[,8:9],proj4string = crs(Embarras_DEM))
Embarras_dist_mat_ADCP <- pointDistance(Embarras_qobs_sp, data.frame(Embarras_cl@lines[[1]]@Lines[[1]]@coords))
Embarras_closest_ADCP <- Embarras_orthos[apply(Embarras_dist_mat_ADCP, 1, FUN = which.min)]
Embarras_merged_ADCP <- do.call(rbind, Embarras_closest_ADCP)
Embarras_ADCP_df <- SpatialLinesDataFrame(Embarras_merged_ADCP, Embarras_corr_Q, match.ID = FALSE)
writeOGR(Embarras_ADCP_df, dsn = ".", layer = "Embarras_ADCP", driver = "ESRI Shapefile")


#Find closest orthos to pressure transducers
Embarras_PT_coords_sp <-SpatialPoints(Embarras_PT_coords[,4:5], proj4string = crs(Embarras_DEM))
Embarras_dist_mat <- pointDistance(Embarras_PT_coords_sp, data.frame(Embarras_cl@lines[[1]]@Lines[[1]]@coords))
Embarras_closest_ortho <- Embarras_orthos[apply(Embarras_dist_mat, 1, FUN = which.min)]

plot(Embarras_DEM)
plot(Embarras_PT_coords_sp, add = TRUE)
for(i in 1:length(Embarras_closest_ortho)){
  plot(Embarras_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Embarras_elev_prof <- foreach(i=1:length(Embarras_closest_ortho), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Embarras_DEM, Embarras_closest_ortho[[i]], method = 'simple')
                              }
proc.time()-ptm
stopCluster(cl)

Embarras_scale_x <- vector()
for(i in 1:length(Embarras_elev_prof)){
  Embarras_scale_x[i] = 2000/length(Embarras_elev_prof[[i]])
}

#Limon
setwd(paste(home_dir,"/PAD/PAD_2019/PAD_2019_PT", sep =""))
Limon_PT_coords <- read.csv("Limon_PT_coords.csv")

Limon_bounding_box <- c(min(Limon_PT_coords$Longitude)-2000, max(Limon_PT_coords$Longitude)+2000,
                        min(Limon_PT_coords$Latitude)-1000, max(Limon_PT_coords$Latitude)+1000)

ptm<-proc.time()
Limon_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Limon_DEM <- crop(Limon_DEM_raw, Limon_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2019/Limon", sep =""))
Limon_cl <-readOGR(dsn = ".", layer = "Limon_CL_UTM")

Limon_qobs_df <- data.frame(crop(PAD_2019_ADCP_UTM, Limon_bounding_box))
Limon_corr_Q <- Limon_qobs_df[Limon_qobs_df$Q > 30 & Limon_qobs_df$Q < 50,]
write.csv(Limon_corr_Q, "Limon_corr_Q.csv")

writeOGR(Limon_qobs_df, dsn = ".", layer = "Limon_qobs")
ptm<-proc.time()
Limon_orthos <- ortho_lines(DEM = Limon_DEM, centerline_coords = 
                              data.frame(Limon_cl@lines[[1]]@Lines[[1]]@coords),
                            scale_xs = 100)
proc.time()-ptm
Limon_qobs_sp <-SpatialPoints(Limon_corr_Q[,8:9],proj4string = crs(Limon_DEM))
Limon_dist_mat_ADCP <- pointDistance(Limon_qobs_sp, data.frame(Limon_cl@lines[[1]]@Lines[[1]]@coords))
Limon_closest_ADCP <- Limon_orthos[apply(Limon_dist_mat_ADCP, 1, FUN = which.min)]
Limon_merged_ADCP <- do.call(rbind, Limon_closest_ADCP)
Limon_ADCP_df <- SpatialLinesDataFrame(Limon_merged_ADCP, Limon_corr_Q, match.ID = FALSE)
writeOGR(Limon_ADCP_df, dsn = ".", layer = "Limon_ADCP", driver = "ESRI Shapefile")


#Find closest orthos to pressure transducers
Limon_PT_coords_sp <-SpatialPoints(Limon_PT_coords[,4:5], proj4string = crs(Limon_DEM))
Limon_dist_mat <- pointDistance(Limon_PT_coords_sp, data.frame(Limon_cl@lines[[1]]@Lines[[1]]@coords))
Limon_closest_ortho <- Limon_orthos[apply(Limon_dist_mat, 1, FUN = which.min)]

plot(Limon_DEM)
plot(Limon_PT_coords_sp, add = TRUE)
for(i in 1:length(Limon_closest_ortho)){
  plot(Limon_closest_ortho[[i]], add = TRUE, col = "blue")
}
#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Limon_elev_prof <- foreach(i=1:length(Limon_closest_ortho), 
                           .combine='rbind', .packages=c('raster')) %dopar% {
                             extract(Limon_DEM, Limon_closest_ortho[[i]], method = 'simple')
                           }
proc.time()-ptm
stopCluster(cl)

Limon_scale_x <- vector()
for(i in 1:length(Limon_elev_prof)){
  Limon_scale_x[i] = 2000/length(Limon_elev_prof[[i]])
}

#Richardson
setwd(paste(home_dir,"/PAD/PAD_2019/PAD_2019_PT", sep =""))
Richardson_PT_coords <- read.csv("Richardson_PT_coords.csv")

Richardson_bounding_box <- c(min(Richardson_PT_coords$Longitude)-2000, max(Richardson_PT_coords$Longitude)+2000,
                             min(Richardson_PT_coords$Latitude)-1000, max(Richardson_PT_coords$Latitude)+1000)

ptm<-proc.time()
Richardson_DEM_raw <- PAD_DEM <- raster("C:/Users/Merritt/Desktop/Research/AGU2018/PAD/Pad_DEM.tif")
Richardson_DEM <- crop(Richardson_DEM_raw, Richardson_bounding_box)
proc.time()-ptm

setwd(paste(home_dir,"/PAD/PAD_2019/Richardson", sep =""))
Richardson_cl <-readOGR(dsn = ".", layer = "Richardson_CL_UTM")

Richardson_qobs_df <- data.frame(crop(PAD_2019_ADCP_UTM, Richardson_bounding_box))
Richardson_corr_Q <- Richardson_qobs_df[Richardson_qobs_df$Q < 25,]
write.csv(Richardson_corr_Q, "Richardson_corr_Q.csv")
Richardson_qobs_shp <- crop(PAD_2019_ADCP_UTM, Richardson_bounding_box)
writeOGR(Richardson_qobs_shp, dsn = ".", layer = "Richardson_qobs", driver = "ESRI Shapefile")
ptm<-proc.time()
Richardson_orthos <- ortho_lines(DEM = Richardson_DEM, centerline_coords = 
                                   data.frame(Richardson_cl@lines[[1]]@Lines[[1]]@coords),
                                 scale_xs = 100)
proc.time()-ptm

#Find closest orthos to pressure transducers
Richardson_PT_coords_sp <-SpatialPoints(Richardson_PT_coords[,4:5], proj4string = crs(Richardson_DEM))
Richardson_dist_mat <- pointDistance(Richardson_PT_coords_sp, data.frame(Richardson_cl@lines[[1]]@Lines[[1]]@coords))
Richardson_closest_ortho <- Richardson_orthos[apply(Richardson_dist_mat, 1, FUN = which.min)]

plot(Richardson_DEM)
plot(Richardson_PT_coords_sp, add = TRUE)
for(i in 1:length(Richardson_closest_ortho)){
  plot(Richardson_closest_ortho[[i]], add = TRUE, col = "blue")
}

Richardson_qobs_sp <-SpatialPoints(Richardson_corr_Q[,8:9],proj4string = crs(Richardson_DEM))
Richardson_dist_mat_ADCP <- pointDistance(Richardson_qobs_sp, data.frame(Richardson_cl@lines[[1]]@Lines[[1]]@coords))
Richardson_closest_ADCP <- Richardson_orthos[apply(Richardson_dist_mat_ADCP, 1, FUN = which.min)]
Richardson_merged_ADCP <- do.call(rbind, Richardson_closest_ADCP)
Richardson_ADCP_df <- SpatialLinesDataFrame(Richardson_merged_ADCP, Richardson_corr_Q, match.ID = FALSE)
writeOGR(Richardson_ADCP_df, dsn = ".", layer = "Richardson_ADCP", driver = "ESRI Shapefile")



#copy cross section profile code at specific location

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Richardson_elev_prof <- foreach(i=1:length(Richardson_closest_ortho), 
                                .combine='rbind', .packages=c('raster')) %dopar% {
                                  extract(Richardson_DEM, Richardson_closest_ortho[[i]], method = 'simple')
                                }
proc.time()-ptm
stopCluster(cl)

Richardson_scale_x <- vector()
for(i in 1:length(Richardson_elev_prof)){
  Richardson_scale_x[i] = 2000/length(Richardson_elev_prof[[i]])
}

library(gtools)
library(bamr)
library(zoo)
library(hydroGOF)
library(sp)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)

###Functions-------
#extract elevation profiles
ortho_lines <- function(DEM, centerline_coords, scale_xs){
  ortho_lines = list()
  UTMcrs = crs(DEM)
  cl_coords_df = data.frame(x=c(1:length(centerline_coords[,1])), lon=centerline_coords[,1], 
                            lat=centerline_coords[,2])
  mean_lon = cl_coords_df$lon[1:(length(cl_coords_df$lon)-1)]+(1/2)*diff(cl_coords_df$lon)
  mean_lat = cl_coords_df$lat[1:(length(cl_coords_df$lat)-1)]+(1/2)*diff(cl_coords_df$lat)
  slope = (-1)/(diff(cl_coords_df$lat)/diff(cl_coords_df$lon))
  for (i in which(!is.na(slope))){
    left_lon = mean_lon-(cos(atan(slope))*scale_xs)
    left_lat = mean_lat-(sin(atan(slope))*scale_xs)
    right_lon = mean_lon+(cos(atan(slope))*scale_xs)
    right_lat = mean_lat+(sin(atan(slope))*scale_xs)
    ortho_lines[[i]] = SpatialLines(list(Lines(Line(cbind(c(left_lon[i],right_lon[i]),
                                                          c(left_lat[i], right_lat[i]))), ID=i)),
                                    proj4string = UTMcrs)
  }
  ortho_lines
}

#Plot elevation profiles with elevation points on top
plot_elev_prof <- function(elev_prof, wse_mat, rivname){
  for(i in 1:length(elev_prof)){
    plot(spline(elev_prof[[i]], n = 1000), main = paste0(rivname, "\n Cross-Section #", i),
         xlab = "Horizontal Distance", ylab = "Elevation (m)", type = "l", 
         ylim = c(min(c(wse_mat[,i], elev_prof[[i]]), na.rm = TRUE), 
                  max(c(wse_mat[,i], elev_prof[[i]]), na.rm = TRUE)))
    points(seq(1, length(elev_prof[[i]]), length.out = length(wse_mat[,i])), wse_mat[,i], 
           pch = 19, col = "blue")
  }
}

#Calculate Width
calc_width <- function(h, elev_prof, upper_l, lower_l, lower_r, upper_r, scale_prof){
  width = matrix(nrow = nrow(h), ncol = length(elev_prof))
  for(xs in 1:length(elev_prof)){
    cropped_elev_prof = elev_prof[[xs]][c(upper_l[xs]:upper_r[xs])]
    s_elev_prof = round(spline(cropped_elev_prof, n=10000)$y, 3)
    scale = (upper_r[xs] - upper_l[xs]) *scale_prof[xs]/10000
    s_left = round(s_elev_prof[1:round((lower_l[xs]- upper_l[xs])*(10000/(upper_r[xs] - upper_l[xs])),0)],1)
    right_start = round((lower_r[xs]-upper_l[xs])*(10000/(upper_r[xs] - upper_l[xs])),0)
    s_right = round(s_elev_prof[right_start:10000],1)
    wse = intersect(s_left, s_right)
    wse = wse[!is.na(wse)]
    width_vec = vector()
    for(i in 1:length(wse)){
      width_vec[i] = (which(s_right==wse[i])[1]+right_start - which(s_left ==wse[i])[1])*scale
    }
    fit = lm(width_vec~wse)
    width[,xs] = h[,xs]*fit$coefficients[2]+fit$coefficients[1]
  }
  width
}

#Vector of distances between pt coordinates
calc_xvec <- function(pt_coords_sp, cl){
  coords = cl@lines[[1]]@Lines[[1]]@coords
  csum = cumsum(pointDistance(coords[1:(nrow(coords)-1),], coords[2:nrow(coords),], lonlat = FALSE))
  dist_mat = pointDistance(pt_coords_sp, data.frame(cl@lines[[1]]@Lines[[1]]@coords))
  xvec = csum[apply(dist_mat, 1, FUN = which.min)]
  return(xvec)
}

#Calculate dA
#Calculate area
#' Calculate partial cross-section area from DAWG-formatted width and height matrices
#' @param w Matrix of widths
#' @param h Matrix of heights(FROM MARK)
calcdA_mat <- function(w, h) {
  stopifnot(all(dim(w) == dim(h)))
  dA <- w
  for (i in 1:nrow(dA)) {
    dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
  }
  dA
}


#' Calculate partial cross-section area from width and height vectors (time series)
#' @param w vector of widths
#' @param h vector of heights(FROM MARK)
calcdA_vec <- function(w, h) {
  words <- order(w)
  warr <- w[words]
  harr <- h[words]
  delh <- c(0, diff(as.numeric(harr)))
  delA <- cumsum(warr * delh)
  dA <- 1:length(w)
  dA[words] <- delA
  dA
}

#Calculate slope
#' Calculate slope based on height and along-river distance
#' @param xvec vector of x (along-stream) distances, measured from
#' upstream-most cross-section
#' @param hmat Matrix (space-down, time-across) of stream heights above
#'arbitrary datum.(FROM MARK)
calcslope <- function(xvec, hmat) {
  dH <- apply(hmat, 2, function(x) c(NA, -diff(x)))
  dX_vec <- c(NA, diff(xvec))
  dX <- matrix(dX_vec, ncol = ncol(hmat), nrow = nrow(hmat),
               byrow = FALSE)
  out <- dH / dX
}

bam_val <- function(bam_est, qobs){
  bam_stats = data.frame()
  detach("package:hydroGOF", unload = TRUE)
  pred = bam_validate(bam_est, qobs)$valdata$qpred
  library(hydroGOF)
  sd = exp(summary(bam_est)$summary[1:length(qobs), 2])
  rrmse = sqrt(mean(((pred - qobs)/qobs)^2))
  nash = NSE(sim = pred, obs = qobs)
  n_rmse = nrmse(sim = pred, obs = qobs)
  p_bias = pbias(sim = pred, obs = qobs)
  kge = KGE(sim = pred, obs = qobs)
  
  rrmse_lower = sqrt(mean((((pred-sd) - qobs)/qobs)^2))
  nash_lower = NSE(sim = (pred-sd), obs = qobs)
  n_rmse_lower = nrmse(sim = (pred-sd), obs = qobs)
  p_bias_lower = pbias(sim = (pred-sd), obs = qobs)
  kge_lower = KGE(sim = (pred-sd), obs = qobs)
  
  rrmse_upper = sqrt(mean((((pred+sd) - qobs)/qobs)^2))
  nash_upper = NSE(sim = (pred+sd), obs = qobs)
  n_rmse_upper = nrmse(sim = (pred+sd), obs = qobs)
  p_bias_upper = pbias(sim = (pred+sd), obs = qobs)
  kge_upper = KGE(sim = (pred+sd), obs = qobs)
  
  df = data.frame(nash = nash, rrmse = rrmse, n_rmse = n_rmse, p_bias = p_bias, 
                  kge = kge, nash_lower = nash_lower, rrmse_lower = rrmse_lower,
                  n_rmse_lower = n_rmse_lower, p_bias_lower = p_bias_lower, kge_lower = kge_lower,
                  nash_upper = nash_upper, rrmse_upper = rrmse_upper, 
                  n_rmse_upper = n_rmse_upper, p_bias_upper = p_bias_upper, kge_upper = kge_upper)
  bam_stats = rbind(df, bam_stats)
  return(list(pred,bam_stats))
}
###Sacramento_Reach_One
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Sacramento/Reach_1"))
load("Sac_WSE_rch1_daily.RData")
load("Sac_WSE_rch1_6hr.RData")
load("Sac_WSE_rch1_hourly.RData")
load("Sac_WSE_rch1_15min.RData")

Sac_rch1_PT_coords <- read.csv("Sac_rch1_PT_coords.csv")
Sac_rch1_PT_coords_sp <- SpatialPoints(Sac_rch1_PT_coords[,c(5,4)], proj4string =  CRS("+proj=utm +zone=10 +datum=WGS84"))
Sac_cl_rch1 <-readOGR(dsn = ".", layer = "Sac_rch1_Lid_cl")

setwd(paste0(home_dir, "/Sacramento/Reach_1/DEM"))
Sac_rch1_DEM_3m <- raster("Sac_DEM_rch1_3m.tif")
Sac_rch1_DEM_NED <- raster("Sac_DEM_rch1.tif")

setwd(paste0(home_dir, "/Sacramento/Reach_1"))
Sac_orthos_rch1_NED <- ortho_lines(Sac_rch1_DEM_NED, data.frame(Sac_cl_rch1@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch1_dist_mat_NED <- pointDistance(Sac_rch1_PT_coords_sp, data.frame(Sac_cl_rch1@lines[[1]]@Lines[[1]]@coords))
Sac_rch1_closest_ortho_NED <- Sac_orthos_rch1_NED[apply(Sac_rch1_dist_mat_NED, 1, FUN = which.min)]

Sac_orthos_rch1_3m <- ortho_lines(Sac_rch1_DEM_3m, data.frame(Sac_cl_rch1@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch1_dist_mat_3m <- pointDistance(Sac_rch1_PT_coords_sp, data.frame(Sac_cl_rch1@lines[[1]]@Lines[[1]]@coords))
Sac_rch1_closest_ortho_3m <- Sac_orthos_rch1_3m[apply(Sac_rch1_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Sac_rch1_elev_prof_NED <- foreach(i=1:length(Sac_rch1_closest_ortho_NED), 
                                  .combine='rbind', .packages=c('raster')) %dopar% {
                                    extract(Sac_rch1_DEM_NED, Sac_rch1_closest_ortho_NED[[i]], method = 'simple')
                                  }
Sac_rch1_elev_prof_3m <- foreach(i=1:length(Sac_rch1_closest_ortho_3m), 
                                 .combine='rbind', .packages=c('raster')) %dopar% {
                                   extract(Sac_rch1_DEM_3m, Sac_rch1_closest_ortho_3m[[i]], method = 'simple')
                                 }
proc.time()-ptm
stopCluster(cl)

Sac_rch1_scale_x_NED <- vector()
for(i in 1:length(Sac_rch1_elev_prof_NED)){
  Sac_rch1_scale_x_NED[i] = 1000/length(Sac_rch1_elev_prof_NED[[i]])
}

Sac_rch1_scale_x_3m <- vector()
for(i in 1:length(Sac_rch1_elev_prof_3m)){
  Sac_rch1_scale_x_3m[i] = 1000/length(Sac_rch1_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Sac_rch1_elev_prof_NED, Sac_WSE_rch1_daily_df, "Sacramento Reach One")
#locator(n = 4)
Sac_rch1_upper_l_NED <- c( 53.9, 57.8, 16.2,  29.1,  29.1,  61.7, 61.7,  61.7, 61.2, 75.3)
Sac_rch1_lower_l_NED <- c( 70.1, 66.7, 58.4,  59.7,  59.7,  71.2, 71.2,  71.2, 64.5, 77.3)
Sac_rch1_lower_r_NED <- c( 85.3, 82.9, 66.9,  80.2,  80.2,  82.1, 82.1,  82.1, 71.4, 89.0)
Sac_rch1_upper_r_NED <- c(106.3, 87.9, 76.2, 131.6, 131.6,  90.3, 90.3,  90.3, 74.7, 91.1)

Sac_rch1_elev_prof_NED_corr <- Sac_rch1_elev_prof_NED[c(1:4, 4, 8, 8, 8, 9, 10)]
Sac_rch1_elev_prof_NED_corr[[6]] <- Sac_rch1_elev_prof_NED[[6]] - 3.5
Sac_rch1_elev_prof_NED_corr[[7]] <- Sac_rch1_elev_prof_NED[[7]] - 3
Sac_rch1_elev_prof_NED_corr[[8]] <- Sac_rch1_elev_prof_NED[[8]] - 4
Sac_rch1_elev_prof_NED_corr[[9]] <- Sac_rch1_elev_prof_NED[[9]] - 4
Sac_rch1_elev_prof_NED_corr[[10]] <- Sac_rch1_elev_prof_NED[[10]] - 4
Sac_rch1_scale_x_NED <- Sac_rch1_scale_x_NED[c(1:4, 4, 8, 8, 8, 9, 10)]

Sac_rch1_width_NED_daily <- calc_width(h = t(Sac_WSE_rch1_daily_df), elev_prof = Sac_rch1_elev_prof_NED_corr, 
                                       upper_l = Sac_rch1_upper_l_NED, lower_l = Sac_rch1_lower_l_NED,
                                       lower_r = Sac_rch1_lower_r_NED, upper_r = Sac_rch1_upper_r_NED,
                                       scale_prof = Sac_rch1_scale_x_NED)
Sac_rch1_dA_NED_daily <- calcdA_mat(w = t(Sac_rch1_width_NED_daily), h = Sac_WSE_rch1_daily_df)

Sac_rch1_width_NED_6hr <- calc_width(h = t(Sac_WSE_rch1_6hr_df), elev_prof = Sac_rch1_elev_prof_NED_corr, 
                                     upper_l = Sac_rch1_upper_l_NED, lower_l = Sac_rch1_lower_l_NED,
                                     lower_r = Sac_rch1_lower_r_NED, upper_r = Sac_rch1_upper_r_NED,
                                     scale_prof = Sac_rch1_scale_x_NED)
Sac_rch1_dA_NED_6hr <- calcdA_mat(w = t(Sac_rch1_width_NED_6hr), h = Sac_WSE_rch1_6hr_df)

Sac_rch1_width_NED_hourly <- calc_width(h = t(Sac_WSE_rch1_hourly_df), elev_prof = Sac_rch1_elev_prof_NED_corr, 
                                        upper_l = Sac_rch1_upper_l_NED, lower_l = Sac_rch1_lower_l_NED,
                                        lower_r = Sac_rch1_lower_r_NED, upper_r = Sac_rch1_upper_r_NED,
                                        scale_prof = Sac_rch1_scale_x_NED)
Sac_rch1_dA_NED_hourly <- calcdA_mat(w = t(Sac_rch1_width_NED_hourly), h = Sac_WSE_rch1_hourly_df)

Sac_rch1_width_NED_15min <- calc_width(h = t(Sac_WSE_rch1_df), elev_prof = Sac_rch1_elev_prof_NED_corr, 
                                       upper_l = Sac_rch1_upper_l_NED, lower_l = Sac_rch1_lower_l_NED,
                                       lower_r = Sac_rch1_lower_r_NED, upper_r = Sac_rch1_upper_r_NED,
                                       scale_prof = Sac_rch1_scale_x_NED)
Sac_rch1_dA_NED_15min <- calcdA_mat(w = t(Sac_rch1_width_NED_15min), h = Sac_WSE_rch1_df)

#plot_elev_prof(Sac_rch1_elev_prof_3m, Sac_WSE_rch1_daily_df, "Sacramento Reach One")
#locator(n = 4)
Sac_rch1_upper_l_3m <- c(171.6, 183.1, 151.9, 175.8, 186.4, 156.4, 157.9, 195.2, 146.6, 201.6)
Sac_rch1_lower_l_3m <- c(193.6, 198.4, 164.0, 197.2, 208.4, 176.9, 172.6, 212.4, 163.9, 216.9)
Sac_rch1_lower_r_3m <- c(215.5, 227.6, 186.5, 219.3, 242.3, 203.7, 205.3, 246.8, 187.6, 250.3)
Sac_rch1_upper_r_3m <- c(230.8, 245.7, 199.2, 238.5, 265.9, 214.3, 211.4, 270.0, 202.7, 267.1)

Sac_rch1_width_3m_daily <- calc_width(h = t(Sac_WSE_rch1_daily_df), elev_prof = Sac_rch1_elev_prof_3m, 
                                      upper_l = Sac_rch1_upper_l_3m, lower_l = Sac_rch1_lower_l_3m,
                                      lower_r = Sac_rch1_lower_r_3m, upper_r = Sac_rch1_upper_r_3m,
                                      scale_prof = Sac_rch1_scale_x_3m)
Sac_rch1_dA_3m_daily <- calcdA_mat(w = t(Sac_rch1_width_3m_daily), h = Sac_WSE_rch1_daily_df)

Sac_rch1_width_3m_6hr <- calc_width(h = t(Sac_WSE_rch1_6hr_df), elev_prof = Sac_rch1_elev_prof_3m, 
                                    upper_l = Sac_rch1_upper_l_3m, lower_l = Sac_rch1_lower_l_3m,
                                    lower_r = Sac_rch1_lower_r_3m, upper_r = Sac_rch1_upper_r_3m,
                                    scale_prof = Sac_rch1_scale_x_3m)
Sac_rch1_dA_3m_6hr <- calcdA_mat(w = t(Sac_rch1_width_3m_6hr), h = Sac_WSE_rch1_6hr_df)

Sac_rch1_width_3m_hourly <- calc_width(h = t(Sac_WSE_rch1_hourly_df), elev_prof = Sac_rch1_elev_prof_3m, 
                                       upper_l = Sac_rch1_upper_l_3m, lower_l = Sac_rch1_lower_l_3m,
                                       lower_r = Sac_rch1_lower_r_3m, upper_r = Sac_rch1_upper_r_3m,
                                       scale_prof = Sac_rch1_scale_x_3m)
Sac_rch1_dA_3m_hourly <- calcdA_mat(w = t(Sac_rch1_width_3m_hourly), h = Sac_WSE_rch1_hourly_df)

Sac_rch1_width_3m_15min <- calc_width(h = t(Sac_WSE_rch1_df), elev_prof = Sac_rch1_elev_prof_3m, 
                                      upper_l = Sac_rch1_upper_l_3m, lower_l = Sac_rch1_lower_l_3m,
                                      lower_r = Sac_rch1_lower_r_3m, upper_r = Sac_rch1_upper_r_3m,
                                      scale_prof = Sac_rch1_scale_x_3m)
Sac_rch1_dA_3m_15min <- calcdA_mat(w = t(Sac_rch1_width_3m_15min), h = Sac_WSE_rch1_df)


xvec_rch1 <- calc_xvec(Sac_rch1_PT_coords_sp, Sac_cl_rch1)
Sac_rch1_xvec <- saveRDS(xvec_rch1, file = "Sac_rch1_xvec.rds")
Sac_rch1_slope_daily <- calcslope(rev(xvec_rch1), hmat = Sac_WSE_rch1_daily_df)
Sac_rch1_slope_daily[1,] <- colMeans(Sac_rch1_slope_daily, na.rm = TRUE)
Sac_rch1_slope_6hr <- calcslope(rev(xvec_rch1), hmat = Sac_WSE_rch1_6hr_df)
Sac_rch1_slope_6hr[1,] <- colMeans(Sac_rch1_slope_6hr, na.rm = TRUE)
Sac_rch1_slope_hourly <- calcslope(rev(xvec_rch1), hmat = Sac_WSE_rch1_hourly_df)
Sac_rch1_slope_hourly[1,] <- colMeans(Sac_rch1_slope_hourly, na.rm = TRUE)
Sac_rch1_slope_15min <- calcslope(rev(xvec_rch1), hmat = Sac_WSE_rch1_df)
Sac_rch1_slope_15min[1,] <- colMeans(Sac_rch1_slope_15min, na.rm = TRUE)

setwd(paste0(home_dir, "/Sacramento/OrigData"))
Sacqobs <- read.csv("SacQobs.csv")$Q*0.028316847
Sac_rch1_qobs_daily<- colMeans(matrix(Sacqobs[133:1284], 96))
Sac_rch1_qobs_6hr<- colMeans(matrix(Sacqobs[85:1308], 24))
Sac_rch1_qobs_hourly<- colMeans(matrix(Sacqobs[85: 1316], 4))
Sac_rch1_qobs_15min<- colMeans(matrix(Sacqobs[85: 1316], 1))

#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Sac_rch1_wse_plot <- matplot(t(Sac_WSE_rch1_daily_df), type = c("l"), col= colfunc,
                             main = "WSE", xlab = "Days",
                             ylab = "WSE (m)",
                             cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch1_slope_plot <- matplot(t(Sac_rch1_slope_daily), type = c("l"), col= colfunc,
                               main = "Slope", xlab = "Days",
                               ylab = "Slope",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch1_width_plot <- matplot(Sac_rch1_width_NED_daily, type = c("l"), col= colfunc,
                               main = "Width (NED)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch1_width_plot <- matplot(Sac_rch1_width_3m_daily, type = c("l"), col= colfunc,
                               main = "Width (Lidar)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch1_dA_plot <- matplot(t(Sac_rch1_dA_NED_daily), type = c("l"), col= colfunc,
                            main = "dA (NED)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch1_dA_plot <- matplot(t(Sac_rch1_dA_3m_daily), type = c("l"), col= colfunc,
                            main = "dA (Lidar)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot

save(Sac_rch1_width_NED_daily, file = "Sacramento_Reach_One_width_NED_daily.RData")
save(Sac_rch1_dA_NED_daily, file = "Sacramento_Reach_One_dA_NED_daily.RData")
save(Sac_rch1_slope_daily, file = "Sacramento_Reach_One_slope_daily.RData")
save(Sac_rch1_qobs_daily, file = "Sacramento_Reach_One_qobs_daily.RData")
save(Sac_rch1_width_NED_6hr, file = "Sacramento_Reach_One_width_NED_6hr.RData")
save(Sac_rch1_dA_NED_6hr, file = "Sacramento_Reach_One_dA_NED_6hr.RData")
save(Sac_rch1_slope_6hr, file = "Sacramento_Reach_One_slope_6hr.RData")
save(Sac_rch1_qobs_6hr, file = "Sacramento_Reach_One_qobs_6hr.RData")
save(Sac_rch1_width_NED_hourly, file = "Sacramento_Reach_One_width_NED_hourly.RData")
save(Sac_rch1_dA_NED_hourly, file = "Sacramento_Reach_One_dA_NED_hourly.RData")
save(Sac_rch1_slope_hourly, file = "Sacramento_Reach_One_slope_hourly.RData")
save(Sac_rch1_qobs_hourly, file = "Sacramento_Reach_One_qobs_hourly.RData")
save(Sac_rch1_width_NED_15min, file = "Sacramento_Reach_One_width_NED_15min.RData")
save(Sac_rch1_dA_NED_15min, file = "Sacramento_Reach_One_dA_NED_15min.RData")
save(Sac_rch1_slope_15min, file = "Sacramento_Reach_One_slope_15min.RData")
save(Sac_rch1_qobs_15min, file = "Sacramento_Reach_One_qobs_15min.RData")
#Run BAM----
#Default-----
Sac_rch1_data_NED_daily <- bam_data(w = t(Sac_rch1_width_NED), s = Sac_rch1_slope, dA = Sac_rch1_dA_NED, Qhat = Sac_rch1_qobs_daily)
Sac_rch1_NED_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "manning_amhg")
Sac_rch1_NED_daily_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "amhg")
Sac_rch1_NED_daily_default_man <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "manning")

Sac_rch1_NED_daily_default_amhg_val <- bam_val(Sac_rch1_NED_daily_default_amhg, Sac_rch1_qobs_daily)
Sac_rch1_NED_daily_default_man_val <- bam_val(Sac_rch1_NED_daily_default_man, Sac_rch1_qobs_daily)
Sac_rch1_NED_daily_default_man_amhg_val <- bam_val(Sac_rch1_NED_daily_default_man_amhg, Sac_rch1_qobs_daily)

Sac_rch1_data_3m_daily <- bam_data(w = t(Sac_rch1_width_3m), s = Sac_rch1_slope, dA = Sac_rch1_dA_3m, Qhat = Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "manning_amhg")
Sac_rch1_3m_daily_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "amhg")
Sac_rch1_3m_daily_default_man <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "manning")

Sac_rch1_3m_daily_default_amhg_val <- bam_val(Sac_rch1_3m_daily_default_amhg, Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_default_man_val <- bam_val(Sac_rch1_3m_daily_default_man, Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_default_man_amhg_val <- bam_val(Sac_rch1_3m_daily_default_man_amhg, Sac_rch1_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily Default NED",
     ylim = c(120, 275))
lines(Sac_rch1_NED_daily_default_amhg_val[[1]], col = "red")
lines(Sac_rch1_NED_daily_default_man_val[[1]], col = "blue")
lines(Sac_rch1_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch1_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch1_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch1_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch1_data_NED_daily <- bam_data(w = t(Sac_rch1_width_NED), s = Sac_rch1_slope, dA = Sac_rch1_dA_NED,
                                    Qhat = Sac_rch1_qobs_daily)
Sac_best_priors <- bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_logQ = log(0.75 * min(Sac_rch1_qobs_daily)),
                              upperbound_logQ = log(1.25 * max(Sac_rch1_qobs_daily)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch1_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch1_NED_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch1_NED_daily_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_NED_daily_best_man <- bam_estimate(bamdata = Sac_rch1_data_NED_daily,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_NED_daily_best_amhg_val <- bam_val(Sac_rch1_NED_daily_best_amhg, Sac_rch1_qobs_daily)
Sac_rch1_NED_daily_best_man_val <- bam_val(Sac_rch1_NED_daily_best_man, Sac_rch1_qobs_daily)
Sac_rch1_NED_daily_best_man_amhg_val <- bam_val(Sac_rch1_NED_daily_best_man_amhg, Sac_rch1_qobs_daily)

Sac_rch1_data_3m_daily <- bam_data(w = t(Sac_rch1_width_3m), s = Sac_rch1_slope, dA = Sac_rch1_dA_3m, Qhat = Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch1_3m_daily_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_3m_daily_best_man <- bam_estimate(bamdata = Sac_rch1_data_3m_daily,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_3m_daily_best_amhg_val <- bam_val(Sac_rch1_3m_daily_best_amhg, Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_best_man_val <- bam_val(Sac_rch1_3m_daily_best_man, Sac_rch1_qobs_daily)
Sac_rch1_3m_daily_best_man_amhg_val <- bam_val(Sac_rch1_3m_daily_best_man_amhg, Sac_rch1_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily best NED",
     ylim = c(135, 155))
lines(Sac_rch1_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch1_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch1_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch1_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch1_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch1_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution-----
#6hr-----
#Default----
Sac_rch1_data_NED_6hr <- bam_data(w = t(Sac_rch1_width_NED_6hr), s = Sac_rch1_slope_6hr, dA = Sac_rch1_dA_NED_6hr,
                                  Qhat = Sac_rch1_qobs_6hr)
Sac_rch1_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "manning_amhg")
Sac_rch1_NED_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "amhg")
Sac_rch1_NED_6hr_default_man <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "manning")

Sac_rch1_NED_6hr_default_amhg_val <- bam_val(Sac_rch1_NED_6hr_default_amhg, Sac_rch1_qobs_6hr)
Sac_rch1_NED_6hr_default_man_val <- bam_val(Sac_rch1_NED_6hr_default_man, Sac_rch1_qobs_6hr)
Sac_rch1_NED_6hr_default_man_amhg_val <- bam_val(Sac_rch1_NED_6hr_default_man_amhg, Sac_rch1_qobs_6hr)

Sac_rch1_data_3m_6hr <- bam_data(w = t(Sac_rch1_width_3m_6hr), s = Sac_rch1_slope_6hr, dA = Sac_rch1_dA_3m_6hr, 
                                 Qhat = Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "manning_amhg")
Sac_rch1_3m_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "amhg")
Sac_rch1_3m_6hr_default_man <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "manning")

Sac_rch1_3m_6hr_default_amhg_val <- bam_val(Sac_rch1_3m_6hr_default_amhg, Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_default_man_val <- bam_val(Sac_rch1_3m_6hr_default_man, Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_default_man_amhg_val <- bam_val(Sac_rch1_3m_6hr_default_man_amhg, Sac_rch1_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr Default NED",
     ylim = c(120, 275))
lines(Sac_rch1_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Sac_rch1_NED_6hr_default_man_val[[1]], col = "blue")
lines(Sac_rch1_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch1_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch1_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch1_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch1_data_NED_6hr <- bam_data(w = t(Sac_rch1_width_NED_6hr), s = Sac_rch1_slope_6hr, dA = Sac_rch1_dA_NED_6hr,
                                  Qhat = Sac_rch1_qobs_6hr)
Sac_best_priors <- bam_priors(bamdata= Sac_rch1_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Sac_rch1_qobs_6hr)),
                              upperbound_logQ = log(1.25 * max(Sac_rch1_qobs_6hr)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch1_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch1_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "manning_amhg", 
                                               bampriors = Sac_best_priors)
Sac_rch1_NED_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_NED_6hr_best_man <- bam_estimate(bamdata = Sac_rch1_data_NED_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_NED_6hr_best_amhg_val <- bam_val(Sac_rch1_NED_6hr_best_amhg, Sac_rch1_qobs_6hr)
Sac_rch1_NED_6hr_best_man_val <- bam_val(Sac_rch1_NED_6hr_best_man, Sac_rch1_qobs_6hr)
Sac_rch1_NED_6hr_best_man_amhg_val <- bam_val(Sac_rch1_NED_6hr_best_man_amhg, Sac_rch1_qobs_6hr)

Sac_rch1_data_3m_6hr <- bam_data(w = t(Sac_rch1_width_3m_6hr), s = Sac_rch1_slope_6hr, dA = Sac_rch1_dA_3m_6hr,
                                 Qhat = Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "manning_amhg",
                                              bampriors = Sac_best_priors)
Sac_rch1_3m_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_3m_6hr_best_man <- bam_estimate(bamdata = Sac_rch1_data_3m_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_3m_6hr_best_amhg_val <- bam_val(Sac_rch1_3m_6hr_best_amhg, Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_best_man_val <- bam_val(Sac_rch1_3m_6hr_best_man, Sac_rch1_qobs_6hr)
Sac_rch1_3m_6hr_best_man_amhg_val <- bam_val(Sac_rch1_3m_6hr_best_man_amhg, Sac_rch1_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr best NED",
     ylim = c(135, 155))
lines(Sac_rch1_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch1_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch1_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch1_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch1_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch1_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Sac_rch1_data_NED_hourly <- bam_data(w = t(Sac_rch1_width_NED_hourly), s = Sac_rch1_slope_hourly, dA = Sac_rch1_dA_NED_hourly,
                                     Qhat = Sac_rch1_qobs_hourly)
Sac_rch1_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "manning_amhg")
Sac_rch1_NED_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "amhg")
Sac_rch1_NED_hourly_default_man <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "manning")

Sac_rch1_NED_hourly_default_amhg_val <- bam_val(Sac_rch1_NED_hourly_default_amhg, Sac_rch1_qobs_hourly)
Sac_rch1_NED_hourly_default_man_val <- bam_val(Sac_rch1_NED_hourly_default_man, Sac_rch1_qobs_hourly)
Sac_rch1_NED_hourly_default_man_amhg_val <- bam_val(Sac_rch1_NED_hourly_default_man_amhg, Sac_rch1_qobs_hourly)

Sac_rch1_data_3m_hourly <- bam_data(w = t(Sac_rch1_width_3m_hourly), s = Sac_rch1_slope_hourly, dA = Sac_rch1_dA_3m_hourly, 
                                    Qhat = Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "manning_amhg")
Sac_rch1_3m_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "amhg")
Sac_rch1_3m_hourly_default_man <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "manning")

Sac_rch1_3m_hourly_default_amhg_val <- bam_val(Sac_rch1_3m_hourly_default_amhg, Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_default_man_val <- bam_val(Sac_rch1_3m_hourly_default_man, Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_default_man_amhg_val <- bam_val(Sac_rch1_3m_hourly_default_man_amhg, Sac_rch1_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly Default NED",
     ylim = c(120, 275))
lines(Sac_rch1_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Sac_rch1_NED_hourly_default_man_val[[1]], col = "blue")
lines(Sac_rch1_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch1_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch1_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch1_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch1_data_NED_hourly <- bam_data(w = t(Sac_rch1_width_NED_hourly), s = Sac_rch1_slope_hourly, dA = Sac_rch1_dA_NED_hourly,
                                     Qhat = Sac_rch1_qobs_hourly)
Sac_best_priors <- bam_priors(bamdata= Sac_rch1_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Sac_rch1_qobs_hourly)),
                              upperbound_logQ = log(1.25 * max(Sac_rch1_qobs_hourly)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch1_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch1_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "manning_amhg", 
                                                  bampriors = Sac_best_priors)
Sac_rch1_NED_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_NED_hourly_best_man <- bam_estimate(bamdata = Sac_rch1_data_NED_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_NED_hourly_best_amhg_val <- bam_val(Sac_rch1_NED_hourly_best_amhg, Sac_rch1_qobs_hourly)
Sac_rch1_NED_hourly_best_man_val <- bam_val(Sac_rch1_NED_hourly_best_man, Sac_rch1_qobs_hourly)
Sac_rch1_NED_hourly_best_man_amhg_val <- bam_val(Sac_rch1_NED_hourly_best_man_amhg, Sac_rch1_qobs_hourly)

Sac_rch1_data_3m_hourly <- bam_data(w = t(Sac_rch1_width_3m_hourly), s = Sac_rch1_slope_hourly, dA = Sac_rch1_dA_3m_hourly,
                                    Qhat = Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "manning_amhg",
                                                 bampriors = Sac_best_priors)
Sac_rch1_3m_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_3m_hourly_best_man <- bam_estimate(bamdata = Sac_rch1_data_3m_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_3m_hourly_best_amhg_val <- bam_val(Sac_rch1_3m_hourly_best_amhg, Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_best_man_val <- bam_val(Sac_rch1_3m_hourly_best_man, Sac_rch1_qobs_hourly)
Sac_rch1_3m_hourly_best_man_amhg_val <- bam_val(Sac_rch1_3m_hourly_best_man_amhg, Sac_rch1_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly best NED",
     ylim = c(135, 155))
lines(Sac_rch1_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch1_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch1_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch1_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch1_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch1_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#15min-----
#Default----
Sac_rch1_data_NED_15min <- bam_data(w = t(Sac_rch1_width_NED_15min), s = Sac_rch1_slope_15min, dA = Sac_rch1_dA_NED_15min,
                                    Qhat = Sac_rch1_qobs_15min)
Sac_rch1_NED_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "manning_amhg")
Sac_rch1_NED_15min_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "amhg")
Sac_rch1_NED_15min_default_man <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "manning")

library(hydroGOF)
Sac_rch1_NED_15min_default_amhg_val <- bam_val(Sac_rch1_NED_15min_default_amhg, Sac_rch1_qobs_15min)
Sac_rch1_NED_15min_default_man_val <- bam_val(Sac_rch1_NED_15min_default_man, Sac_rch1_qobs_15min)
Sac_rch1_NED_15min_default_man_amhg_val <- bam_val(Sac_rch1_NED_15min_default_man_amhg, Sac_rch1_qobs_15min)

Sac_rch1_data_3m_15min <- bam_data(w = t(Sac_rch1_width_3m_15min), s = Sac_rch1_slope_15min, dA = Sac_rch1_dA_3m_15min, 
                                   Qhat = Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "manning_amhg")
Sac_rch1_3m_15min_default_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "amhg")
Sac_rch1_3m_15min_default_man <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "manning")

library(hydroGOF)
Sac_rch1_3m_15min_default_amhg_val <- bam_val(Sac_rch1_3m_15min_default_amhg, Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_default_man_val <- bam_val(Sac_rch1_3m_15min_default_man, Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_default_man_amhg_val <- bam_val(Sac_rch1_3m_15min_default_man_amhg, Sac_rch1_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min Default NED",
     ylim = c(120, 275))
lines(Sac_rch1_NED_15min_default_amhg_val[[1]], col = "red")
lines(Sac_rch1_NED_15min_default_man_val[[1]], col = "blue")
lines(Sac_rch1_NED_15min_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch1_3m_15min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch1_3m_15min_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch1_3m_15min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch1_data_NED_15min <- bam_data(w = t(Sac_rch1_width_NED_15min), s = Sac_rch1_slope_15min, dA = Sac_rch1_dA_NED_15min,
                                    Qhat = Sac_rch1_qobs_15min)
Sac_best_priors <- bam_priors(bamdata= Sac_rch1_data_NED_15min, lowerbound_logQ = log(0.75 * min(Sac_rch1_qobs_15min)),
                              upperbound_logQ = log(1.25 * max(Sac_rch1_qobs_15min)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch1_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch1_NED_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch1_NED_15min_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_NED_15min_best_man <- bam_estimate(bamdata = Sac_rch1_data_NED_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_NED_15min_best_amhg_val <- bam_val(Sac_rch1_NED_15min_best_amhg, Sac_rch1_qobs_15min)
Sac_rch1_NED_15min_best_man_val <- bam_val(Sac_rch1_NED_15min_best_man, Sac_rch1_qobs_15min)
Sac_rch1_NED_15min_best_man_amhg_val <- bam_val(Sac_rch1_NED_15min_best_man_amhg, Sac_rch1_qobs_15min)

Sac_rch1_data_3m_15min <- bam_data(w = t(Sac_rch1_width_3m_15min), s = Sac_rch1_slope_15min, dA = Sac_rch1_dA_3m_15min,
                                   Qhat = Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch1_3m_15min_best_amhg <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch1_3m_15min_best_man <- bam_estimate(bamdata = Sac_rch1_data_3m_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch1_3m_15min_best_amhg_val <- bam_val(Sac_rch1_3m_15min_best_amhg, Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_best_man_val <- bam_val(Sac_rch1_3m_15min_best_man, Sac_rch1_qobs_15min)
Sac_rch1_3m_15min_best_man_amhg_val <- bam_val(Sac_rch1_3m_15min_best_man_amhg, Sac_rch1_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch1_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min best NED",
     ylim = c(135, 155))
lines(Sac_rch1_NED_15min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch1_NED_15min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch1_NED_15min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch1_3m_15min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch1_3m_15min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch1_3m_15min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Sac_rch1_daily_lower_logQ <- log(seq(exp(max(apply(log(Sac_rch1_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                     min(Sac_rch1_qobs_15min),length.out = 10))
Sac_rch1_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_logQ = Sac_rch1_daily_lower_logQ[i])
}
#upperbound_logQ
Sac_rch1_daily_upper_logQ <- log(seq(max(Sac_rch1_qobs_15min), exp(min(apply(log(Sac_rch1_width_NED_daily), 2, max)) + log(40) + log(5)),
                                     length.out = 10))
Sac_rch1_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_logQ = Sac_rch1_daily_upper_logQ[i])
}
#lowerbound_A0
Sac_rch1_daily_lower_A0 <- seq(0, 30, length.out = 10)
Sac_rch1_daily_lower_A0_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_A0 = Sac_rch1_daily_lower_A0[i])
}
#upperbound_A0
Sac_rch1_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Sac_rch1_daily_upper_A0_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_A0 = Sac_rch1_daily_upper_A0[i])
}
#lowerbound_logn
Sac_rch1_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Sac_rch1_daily_lower_logn_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_logn = Sac_rch1_daily_lower_logn[i])
}
#upperbound_logn
Sac_rch1_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Sac_rch1_daily_upper_logn_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_logn = Sac_rch1_daily_upper_logn[i])
}
#lowerbound_logQc
Sac_rch1_daily_lower_logQc <- seq(0, log(min(Sac_rch1_qobs_15min)),length.out = 10)
Sac_rch1_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_logQc = Sac_rch1_daily_lower_logQc[i])
}
#upperbound_logQc
Sac_rch1_daily_upper_logQc <- seq(log(max(Sac_rch1_qobs_15min)), 10, length.out = 10)
Sac_rch1_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_logQc = Sac_rch1_daily_upper_logQc[i])
}
#lowerbound_logWc
Sac_rch1_daily_lower_logWc <- seq(1, log(min(Sac_rch1_width_NED_daily)),length.out = 10)
Sac_rch1_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_logWc = Sac_rch1_daily_lower_logWc[i])
}
#upperbound_logWc
Sac_rch1_daily_upper_logWc <- seq(log(max(Sac_rch1_width_NED_daily)), 8, length.out = 10)
Sac_rch1_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_logWc = Sac_rch1_daily_upper_logWc[i])
}
#lowerbound_b
Sac_rch1_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Sac_rch1_daily_lower_b_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_lower_b_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, lowerbound_b = Sac_rch1_daily_lower_b[i])
}
#upperbound_b
Sac_rch1_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Sac_rch1_daily_upper_b_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_upper_b_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, upperbound_b = Sac_rch1_daily_upper_b[i])
}
#sigma_man
Sac_rch1_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Sac_rch1_daily_sigma_man_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, sigma_man = Sac_rch1_daily_sigma_man[i])
}
#sigma_amhg
Sac_rch1_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Sac_rch1_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, sigma_amhg = Sac_rch1_daily_sigma_amhg[i])
}
#logQc_hat
Sac_rch1_daily_logQc_hat <- log(seq(min(Sac_rch1_qobs_15min), max(Sac_rch1_qobs_15min), length.out = 10))
Sac_rch1_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logQc_hat = Sac_rch1_daily_logQc_hat[i])
}
#logWc_hat
Sac_rch1_daily_logWc_hat <- log(seq(min(Sac_rch1_width_NED_daily), max(Sac_rch1_width_NED_daily), length.out = 10))
Sac_rch1_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logWc_hat = Sac_rch1_daily_logWc_hat[i])
}
#b_hat
Sac_rch1_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Sac_rch1_daily_b_hat_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_b_hat_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, b_hat = Sac_rch1_daily_b_hat[i])
}
#logA0_hat
Sac_rch1_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Sac_rch1_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logA0_hat = Sac_rch1_daily_logA0_hat[i])
}
#logn_hat
Sac_rch1_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Sac_rch1_daily_logn_hat_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logn_hat = Sac_rch1_daily_logn_hat[i])
}
#logQ_sd
Sac_rch1_daily_logQ_sd <- seq(sd(log(Sac_rch1_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch1_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logQ_sd = Sac_rch1_daily_logQ_sd[i])
}
#logQc_sd
Sac_rch1_daily_logQc_sd <- seq(sd(log(Sac_rch1_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch1_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logQc_sd = Sac_rch1_daily_logQc_sd[i])
}
#logWc_Sd
Sac_rch1_daily_logWc_Sd <- seq(log(sd(Sac_rch1_width_NED_daily)), 4.712493, length.out = 10)
Sac_rch1_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logWc_sd = Sac_rch1_daily_logWc_Sd[i])
}
#b_sd
Sac_rch1_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Sac_rch1_daily_b_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_b_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, b_sd = Sac_rch1_daily_b_sd[i])
}
#logA0_sd
Sac_rch1_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Sac_rch1_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logA0_sd = Sac_rch1_daily_logA0_sd[i])
}
#logn_sd
Sac_rch1_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Sac_rch1_daily_logn_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, logn_sd = Sac_rch1_daily_logn_sd[i])
}
#Werr_sd
Sac_rch1_daily_Werr_sd <- seq(5, 15, length.out = 10)
Sac_rch1_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, Werr_sd = Sac_rch1_daily_Werr_sd[i])
}
#Serr_sd
Sac_rch1_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Sac_rch1_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, Serr_sd = Sac_rch1_daily_Serr_sd[i])
}
#dAerr_sd
Sac_rch1_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Sac_rch1_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Sac_rch1_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch1_data_NED_daily, dAerr_sd = Sac_rch1_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Sac_rch1_daily_manning_list <- list(Sac_rch1_daily_lower_logQ_priors,
                                    Sac_rch1_daily_upper_logQ_priors,
                                    Sac_rch1_daily_lower_A0_priors,
                                    Sac_rch1_daily_upper_A0_priors,
                                    Sac_rch1_daily_lower_logn_priors,
                                    Sac_rch1_daily_upper_logn_priors,
                                    Sac_rch1_daily_sigma_man_priors,
                                    Sac_rch1_daily_logA0_hat_priors,
                                    Sac_rch1_daily_logn_hat_priors,
                                    Sac_rch1_daily_logQ_sd_priors,
                                    Sac_rch1_daily_logA0_sd_priors,
                                    Sac_rch1_daily_logn_sd_priors,
                                    Sac_rch1_daily_Werr_sd_priors,
                                    Sac_rch1_daily_Serr_sd_priors,
                                    Sac_rch1_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Sac_rch1_data_NED_daily, variant = "manning", bampriors = Sac_rch1_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Sac_rch1_daily_amhg_list <- list(Sac_rch1_daily_lower_logQ_priors,
                                 Sac_rch1_daily_upper_logQ_priors,
                                 Sac_rch1_daily_lower_logQc_priors,
                                 Sac_rch1_daily_upper_logQc_priors,
                                 Sac_rch1_daily_lower_logWc_priors,
                                 Sac_rch1_daily_upper_logWc_priors,
                                 Sac_rch1_daily_lower_b_priors,
                                 Sac_rch1_daily_upper_b_priors,
                                 Sac_rch1_daily_sigma_amhg_priors,
                                 Sac_rch1_daily_logQc_hat_priors,
                                 Sac_rch1_daily_logWc_hat_priors,
                                 Sac_rch1_daily_b_hat_priors,
                                 Sac_rch1_daily_logQ_sd_priors,
                                 Sac_rch1_daily_logQc_sd_priors,
                                 Sac_rch1_daily_logWc_Sd_priors,
                                 Sac_rch1_daily_b_sd_priors,
                                 Sac_rch1_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Sac_rch1_data_NED_daily, variant = "amhg", bampriors = Sac_rch1_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Sac_rch1_daily_man_amhg_list <- list(Sac_rch1_daily_lower_logQ_priors,
                                     Sac_rch1_daily_upper_logQ_priors,
                                     Sac_rch1_daily_lower_A0_priors,
                                     Sac_rch1_daily_upper_A0_priors,
                                     Sac_rch1_daily_lower_logn_priors,
                                     Sac_rch1_daily_upper_logn_priors,
                                     Sac_rch1_daily_lower_logQc_priors,
                                     Sac_rch1_daily_upper_logQc_priors,
                                     Sac_rch1_daily_lower_logWc_priors,
                                     Sac_rch1_daily_upper_logWc_priors,
                                     Sac_rch1_daily_lower_b_priors,
                                     Sac_rch1_daily_upper_b_priors,
                                     Sac_rch1_daily_sigma_man_priors,
                                     Sac_rch1_daily_sigma_amhg_priors,
                                     Sac_rch1_daily_logQc_hat_priors,
                                     Sac_rch1_daily_logWc_hat_priors,
                                     Sac_rch1_daily_b_hat_priors,
                                     Sac_rch1_daily_logA0_hat_priors,
                                     Sac_rch1_daily_logn_hat_priors,
                                     Sac_rch1_daily_logQ_sd_priors,
                                     Sac_rch1_daily_logQc_sd_priors,
                                     Sac_rch1_daily_logWc_Sd_priors,
                                     Sac_rch1_daily_b_sd_priors,
                                     Sac_rch1_daily_logA0_sd_priors,
                                     Sac_rch1_daily_logn_sd_priors,
                                     Sac_rch1_daily_Werr_sd_priors,
                                     Sac_rch1_daily_Serr_sd_priors,
                                     Sac_rch1_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_1/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Sac_rch1_data_NED_daily, variant = "manning_amhg", bampriors = Sac_rch1_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------


###Sacramento_Reach_Two
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Sacramento/Reach_2"))
load("Sac_WSE_rch2_daily.RData")
load("Sac_WSE_rch2_6hr.RData")
load("Sac_WSE_rch2_hourly.RData")
load("Sac_WSE_rch2_15min.RData")

Sac_rch2_PT_coords <- read.csv("Sac_rch2_PT_coords.csv")
Sac_rch2_PT_coords_sp <- SpatialPoints(Sac_rch2_PT_coords[,c(5,4)], proj4string =  CRS("+proj=utm +zone=10 +datum=WGS84"))
Sac_cl_rch2 <-readOGR(dsn = ".", layer = "Sac_rch2_Lid_cl")

setwd(paste0(home_dir, "/Sacramento/Reach_2/DEM"))
Sac_rch2_DEM_3m <- raster("Sac_DEM_rch2_3m.tif")
Sac_rch2_DEM_NED <- raster("Sac_DEM_rch2.tif")

setwd(paste0(home_dir, "/Sacramento/Reach_2"))
Sac_orthos_rch2_NED <- ortho_lines(Sac_rch2_DEM_NED, data.frame(Sac_cl_rch2@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch2_dist_mat_NED <- pointDistance(Sac_rch2_PT_coords_sp, data.frame(Sac_cl_rch2@lines[[1]]@Lines[[1]]@coords))
Sac_rch2_closest_ortho_NED <- Sac_orthos_rch2_NED[apply(Sac_rch2_dist_mat_NED, 1, FUN = which.min)]

Sac_orthos_rch2_3m <- ortho_lines(Sac_rch2_DEM_3m, data.frame(Sac_cl_rch2@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch2_dist_mat_3m <- pointDistance(Sac_rch2_PT_coords_sp, data.frame(Sac_cl_rch2@lines[[1]]@Lines[[1]]@coords))
Sac_rch2_closest_ortho_3m <- Sac_orthos_rch2_3m[apply(Sac_rch2_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Sac_rch2_elev_prof_NED <- foreach(i=1:length(Sac_rch2_closest_ortho_NED), 
                                  .combine='rbind', .packages=c('raster')) %dopar% {
                                    extract(Sac_rch2_DEM_NED, Sac_rch2_closest_ortho_NED[[i]], method = 'simple')
                                  }
Sac_rch2_elev_prof_3m <- foreach(i=1:length(Sac_rch2_closest_ortho_3m), 
                                 .combine='rbind', .packages=c('raster')) %dopar% {
                                   extract(Sac_rch2_DEM_3m, Sac_rch2_closest_ortho_3m[[i]], method = 'simple')
                                 }
proc.time()-ptm
stopCluster(cl)

Sac_rch2_scale_x_NED <- vector()
for(i in 1:length(Sac_rch2_elev_prof_NED)){
  Sac_rch2_scale_x_NED[i] = 1000/length(Sac_rch2_elev_prof_NED[[i]])
}

Sac_rch2_scale_x_3m <- vector()
for(i in 1:length(Sac_rch2_elev_prof_3m)){
  Sac_rch2_scale_x_3m[i] = 1000/length(Sac_rch2_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Sac_rch2_elev_prof_NED, Sac_WSE_rch2_daily_df, "Sacramento Reach Two")
#locator(n = 4)
Sac_rch2_upper_l_NED <- c(75.0, 67.2, 59.7, 45.5, 64.2, 73.7)
Sac_rch2_lower_l_NED <- c(78.9, 73.3, 72.4, 49.3, 69.7, 77.6)
Sac_rch2_lower_r_NED <- c(82.8, 77.7, 82.8, 55.6, 82.5, 84.6)
Sac_rch2_upper_r_NED <- c(85.9, 91.3, 87.4, 57.8, 88.5, 87.7)

Sac_rch2_elev_prof_NED_corr <- Sac_rch2_elev_prof_NED
Sac_rch2_elev_prof_NED_corr[[1]] <- Sac_rch2_elev_prof_NED[[1]] - 3
Sac_rch2_elev_prof_NED_corr[[2]] <- Sac_rch2_elev_prof_NED[[2]] - 3
Sac_rch2_elev_prof_NED_corr[[3]] <- Sac_rch2_elev_prof_NED[[3]] - 3
Sac_rch2_elev_prof_NED_corr[[4]] <- Sac_rch2_elev_prof_NED[[4]] - 3
Sac_rch2_elev_prof_NED_corr[[5]] <- Sac_rch2_elev_prof_NED[[5]] - 3
Sac_rch2_elev_prof_NED_corr[[6]] <- Sac_rch2_elev_prof_NED[[6]] - 3

Sac_rch2_width_NED_daily <- calc_width(h = t(Sac_WSE_rch2_daily_df), elev_prof = Sac_rch2_elev_prof_NED_corr, 
                                       upper_l = Sac_rch2_upper_l_NED, lower_l = Sac_rch2_lower_l_NED,
                                       lower_r = Sac_rch2_lower_r_NED, upper_r = Sac_rch2_upper_r_NED,
                                       scale_prof = Sac_rch2_scale_x_NED)
Sac_rch2_width_NED_daily <- na.approx(Sac_rch2_width_NED_daily)
Sac_rch2_dA_NED_daily <- calcdA_mat(w = t(Sac_rch2_width_NED_daily), h = Sac_WSE_rch2_daily_df)

Sac_rch2_width_NED_6hr <- calc_width(h = t(Sac_WSE_rch2_6hr_df), elev_prof = Sac_rch2_elev_prof_NED_corr, 
                                     upper_l = Sac_rch2_upper_l_NED, lower_l = Sac_rch2_lower_l_NED,
                                     lower_r = Sac_rch2_lower_r_NED, upper_r = Sac_rch2_upper_r_NED,
                                     scale_prof = Sac_rch2_scale_x_NED)
Sac_rch2_width_NED_6hr <- na.approx(Sac_rch2_width_NED_6hr)
Sac_rch2_dA_NED_6hr <- calcdA_mat(w = t(Sac_rch2_width_NED_6hr), h = Sac_WSE_rch2_6hr_df)

Sac_rch2_width_NED_hourly <- calc_width(h = t(Sac_WSE_rch2_hourly_df), elev_prof = Sac_rch2_elev_prof_NED_corr, 
                                        upper_l = Sac_rch2_upper_l_NED, lower_l = Sac_rch2_lower_l_NED,
                                        lower_r = Sac_rch2_lower_r_NED, upper_r = Sac_rch2_upper_r_NED,
                                        scale_prof = Sac_rch2_scale_x_NED)
Sac_rch2_width_NED_hourly <- na.approx(Sac_rch2_width_NED_hourly)
Sac_rch2_dA_NED_hourly <- calcdA_mat(w = t(Sac_rch2_width_NED_hourly), h = Sac_WSE_rch2_hourly_df)

Sac_rch2_width_NED_15min <- calc_width(h = t(Sac_WSE_rch2_df), elev_prof = Sac_rch2_elev_prof_NED_corr, 
                                       upper_l = Sac_rch2_upper_l_NED, lower_l = Sac_rch2_lower_l_NED,
                                       lower_r = Sac_rch2_lower_r_NED, upper_r = Sac_rch2_upper_r_NED,
                                       scale_prof = Sac_rch2_scale_x_NED)
Sac_rch2_width_NED_15min <- na.approx(Sac_rch2_width_NED_15min)
Sac_rch2_dA_NED_15min <- calcdA_mat(w = t(Sac_rch2_width_NED_15min), h = Sac_WSE_rch2_df)

plot_elev_prof(Sac_rch2_elev_prof_3m, Sac_WSE_rch2_daily_df, "Sacramento Reach Two")
#locator(n = 4)
Sac_rch2_upper_l_3m <- c(195.7, 199.1, 199.9, 154.5, 183.6, 199.4)
Sac_rch2_lower_l_3m <- c(219.5, 211.8, 212.9, 164.1, 194.9, 209.8)
Sac_rch2_lower_r_3m <- c(249.4, 246.4, 251.9, 190.8, 236.8, 243.2)
Sac_rch2_upper_r_3m <- c(267.0, 266.6, 267.2, 201.0, 248.8, 258.0)

Sac_rch2_width_3m_daily <- calc_width(h = t(Sac_WSE_rch2_daily_df), elev_prof = Sac_rch2_elev_prof_3m, 
                                      upper_l = Sac_rch2_upper_l_3m, lower_l = Sac_rch2_lower_l_3m,
                                      lower_r = Sac_rch2_lower_r_3m, upper_r = Sac_rch2_upper_r_3m,
                                      scale_prof = Sac_rch2_scale_x_3m)
Sac_rch2_width_3m_daily <- na.approx(Sac_rch2_width_3m_daily)
Sac_rch2_dA_3m_daily <- calcdA_mat(w = t(Sac_rch2_width_3m_daily), h = Sac_WSE_rch2_daily_df)

Sac_rch2_width_3m_6hr <- calc_width(h = t(Sac_WSE_rch2_6hr_df), elev_prof = Sac_rch2_elev_prof_3m, 
                                    upper_l = Sac_rch2_upper_l_3m, lower_l = Sac_rch2_lower_l_3m,
                                    lower_r = Sac_rch2_lower_r_3m, upper_r = Sac_rch2_upper_r_3m,
                                    scale_prof = Sac_rch2_scale_x_3m)
Sac_rch2_width_3m_6hr <- na.approx(Sac_rch2_width_3m_6hr)
Sac_rch2_dA_3m_6hr <- calcdA_mat(w = t(Sac_rch2_width_3m_6hr), h = Sac_WSE_rch2_6hr_df)

Sac_rch2_width_3m_hourly <- calc_width(h = t(Sac_WSE_rch2_hourly_df), elev_prof = Sac_rch2_elev_prof_3m, 
                                       upper_l = Sac_rch2_upper_l_3m, lower_l = Sac_rch2_lower_l_3m,
                                       lower_r = Sac_rch2_lower_r_3m, upper_r = Sac_rch2_upper_r_3m,
                                       scale_prof = Sac_rch2_scale_x_3m)
Sac_rch2_width_3m_hourly <- na.approx(Sac_rch2_width_3m_hourly)
Sac_rch2_dA_3m_hourly <- calcdA_mat(w = t(Sac_rch2_width_3m_hourly), h = Sac_WSE_rch2_hourly_df)

Sac_rch2_width_3m_15min <- calc_width(h = t(Sac_WSE_rch2_df), elev_prof = Sac_rch2_elev_prof_3m, 
                                      upper_l = Sac_rch2_upper_l_3m, lower_l = Sac_rch2_lower_l_3m,
                                      lower_r = Sac_rch2_lower_r_3m, upper_r = Sac_rch2_upper_r_3m,
                                      scale_prof = Sac_rch2_scale_x_3m)
Sac_rch2_width_3m_15min <- na.approx(Sac_rch2_width_3m_15min)
Sac_rch2_dA_3m_15min <- calcdA_mat(w = t(Sac_rch2_width_3m_15min), h = Sac_WSE_rch2_df)


xvec_rch2 <- calc_xvec(Sac_rch2_PT_coords_sp, Sac_cl_rch2)
Sac_rch2_xvec <- saveRDS(xvec_rch2, file = "Sac_rch2_xvec.rds")
Sac_rch2_slope_daily <- calcslope(rev(xvec_rch2), hmat = Sac_WSE_rch2_daily_df)
Sac_rch2_slope_daily <- na.replace(Sac_rch2_slope_daily, mean(Sac_rch2_slope_daily))
Sac_rch2_slope_daily[1,] <- colMeans(Sac_rch2_slope_daily, na.rm = TRUE)
Sac_rch2_slope_6hr <- calcslope(rev(xvec_rch2), hmat = Sac_WSE_rch2_6hr_df)
Sac_rch2_slope_6hr <- na.approx(Sac_rch2_slope_6hr)
Sac_rch2_slope_6hr[1,] <- colMeans(Sac_rch2_slope_6hr, na.rm = TRUE)
Sac_rch2_slope_hourly <- calcslope(rev(xvec_rch2), hmat = Sac_WSE_rch2_hourly_df)
Sac_rch2_slope_hourly <- na.approx(Sac_rch2_slope_hourly)
Sac_rch2_slope_hourly[1,] <- colMeans(Sac_rch2_slope_hourly, na.rm = TRUE)
Sac_rch2_slope_15min <- calcslope(rev(xvec_rch2), hmat = Sac_WSE_rch2_df)
Sac_rch2_slope_15min <- na.approx(Sac_rch2_slope_15min)
Sac_rch2_slope_15min[1,] <- colMeans(Sac_rch2_slope_15min, na.rm = TRUE)

setwd(paste0(home_dir, "/Sacramento/OrigData"))
Sacqobs <- read.csv("SacQobs.csv")$Q*0.028316847
Sac_rch2_qobs_daily<- colMeans(matrix(Sacqobs[1330:6609], 96))
Sac_rch2_qobs_6hr<- colMeans(matrix(Sacqobs[1282:6609], 24))
Sac_rch2_qobs_hourly<- colMeans(matrix(Sacqobs[1262:6609], 4))
Sac_rch2_qobs_15min<- colMeans(matrix(Sacqobs[1261:6609], 1))

save(Sac_rch2_width_NED_daily, file = "Sacramento_Reach_Two_width_NED_daily.RData")
save(Sac_rch2_dA_NED_daily, file = "Sacramento_Reach_Two_dA_NED_daily.RData")
save(Sac_rch2_slope_daily, file = "Sacramento_Reach_Two_slope_daily.RData")
save(Sac_rch2_qobs_daily, file = "Sacramento_Reach_Two_qobs_daily.RData")
save(Sac_rch2_width_NED_6hr, file = "Sacramento_Reach_Two_width_NED_6hr.RData")
save(Sac_rch2_dA_NED_6hr, file = "Sacramento_Reach_Two_dA_NED_6hr.RData")
save(Sac_rch2_slope_6hr, file = "Sacramento_Reach_Two_slope_6hr.RData")
save(Sac_rch2_qobs_6hr, file = "Sacramento_Reach_Two_qobs_6hr.RData")
save(Sac_rch2_width_NED_hourly, file = "Sacramento_Reach_Two_width_NED_hourly.RData")
save(Sac_rch2_dA_NED_hourly, file = "Sacramento_Reach_Two_dA_NED_hourly.RData")
save(Sac_rch2_slope_hourly, file = "Sacramento_Reach_Two_slope_hourly.RData")
save(Sac_rch2_qobs_hourly, file = "Sacramento_Reach_Two_qobs_hourly.RData")
save(Sac_rch2_width_NED_15min, file = "Sacramento_Reach_Two_width_NED_15min.RData")
save(Sac_rch2_dA_NED_15min, file = "Sacramento_Reach_Two_dA_NED_15min.RData")
save(Sac_rch2_slope_15min, file = "Sacramento_Reach_Two_slope_15min.RData")
save(Sac_rch2_qobs_15min, file = "Sacramento_Reach_Two_qobs_15min.RData")
#Run BAM----

#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Sac_rch2_wse_plot <- matplot(t(Sac_WSE_rch2_daily_df), type = c("l"), col= colfunc,
                             main = "WSE", xlab = "Days",
                             ylab = "WSE (m)",
                             cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch2_slope_plot <- matplot(t(Sac_rch2_slope_daily), type = c("l"), col= colfunc,
                               main = "Slope", xlab = "Days",
                               ylab = "Slope",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch2_width_plot <- matplot(Sac_rch2_width_NED_daily, type = c("l"), col= colfunc,
                               main = "Width (NED)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch2_width_plot <- matplot(Sac_rch2_width_3m_daily, type = c("l"), col= colfunc,
                               main = "Width (Lidar)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch2_dA_plot <- matplot(t(Sac_rch2_dA_NED_daily), type = c("l"), col= colfunc,
                            main = "dA (NED)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch2_dA_plot <- matplot(t(Sac_rch2_dA_3m_daily), type = c("l"), col= colfunc,
                            main = "dA (Lidar)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
Sac_rch2_data_NED_daily <- bam_data(w = t(Sac_rch2_width_NED_daily), s = Sac_rch2_slope_daily, dA = Sac_rch2_dA_NED_daily,
                                    Qhat = Sac_rch2_qobs_daily)
Sac_rch2_NED_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "manning_amhg")
Sac_rch2_NED_daily_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "amhg")
Sac_rch2_NED_daily_default_man <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "manning")
library(hydroGOF)
Sac_rch2_NED_daily_default_amhg_val <- bam_val(Sac_rch2_NED_daily_default_amhg, Sac_rch2_qobs_daily)
Sac_rch2_NED_daily_default_man_val <- bam_val(Sac_rch2_NED_daily_default_man, Sac_rch2_qobs_daily)
Sac_rch2_NED_daily_default_man_amhg_val <- bam_val(Sac_rch2_NED_daily_default_man_amhg, Sac_rch2_qobs_daily)

Sac_rch2_data_3m_daily <- bam_data(w = t(Sac_rch2_width_3m_daily), s = Sac_rch2_slope_daily, dA = Sac_rch2_dA_3m_daily,
                                   Qhat = Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "manning_amhg")
Sac_rch2_3m_daily_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "amhg")
Sac_rch2_3m_daily_default_man <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "manning")

Sac_rch2_3m_daily_default_amhg_val <- bam_val(Sac_rch2_3m_daily_default_amhg, Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_default_man_val <- bam_val(Sac_rch2_3m_daily_default_man, Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_default_man_amhg_val <- bam_val(Sac_rch2_3m_daily_default_man_amhg, Sac_rch2_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily Default NED",
     ylim = c(100, 220))
lines(Sac_rch2_NED_daily_default_amhg_val[[1]], col = "red")
lines(Sac_rch2_NED_daily_default_man_val[[1]], col = "blue")
lines(Sac_rch2_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch2_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch2_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch2_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch2_data_NED_daily <- bam_data(w = t(Sac_rch2_width_NED_daily), s = Sac_rch2_slope_daily, dA = Sac_rch2_dA_NED_daily,
                                    Qhat = Sac_rch2_qobs_daily)
Sac_best_priors <- bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_logQ = log(0.5 * min(Sac_rch2_qobs_daily)),
                              upperbound_logQ = log(1.75 * max(Sac_rch2_qobs_daily)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch2_slope_daily, na.rm = TRUE),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch2_NED_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch2_NED_daily_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_NED_daily_best_man <- bam_estimate(bamdata = Sac_rch2_data_NED_daily,variant = "manning", bampriors = Sac_best_priors)
library(hydroGOF)
Sac_rch2_NED_daily_best_amhg_val <- bam_val(Sac_rch2_NED_daily_best_amhg, Sac_rch2_qobs_daily)
Sac_rch2_NED_daily_best_man_val <- bam_val(Sac_rch2_NED_daily_best_man, Sac_rch2_qobs_daily)
Sac_rch2_NED_daily_best_man_amhg_val <- bam_val(Sac_rch2_NED_daily_best_man_amhg, Sac_rch2_qobs_daily)

Sac_rch2_data_3m_daily <- bam_data(w = t(Sac_rch2_width_3m_daily), s = Sac_rch2_slope_daily, dA = Sac_rch2_dA_3m_daily,
                                   Qhat = Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch2_3m_daily_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_3m_daily_best_man <- bam_estimate(bamdata = Sac_rch2_data_3m_daily,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_3m_daily_best_amhg_val <- bam_val(Sac_rch2_3m_daily_best_amhg, Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_best_man_val <- bam_val(Sac_rch2_3m_daily_best_man, Sac_rch2_qobs_daily)
Sac_rch2_3m_daily_best_man_amhg_val <- bam_val(Sac_rch2_3m_daily_best_man_amhg, Sac_rch2_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily best NED",
     ylim = c(100, 200))
lines(Sac_rch2_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch2_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch2_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch2_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch2_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch2_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
Sac_rch2_data_NED_6hr <- bam_data(w = t(Sac_rch2_width_NED_6hr), s = Sac_rch2_slope_6hr, dA = Sac_rch2_dA_NED_6hr,
                                  Qhat = Sac_rch2_qobs_6hr)
Sac_rch2_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "manning_amhg")
Sac_rch2_NED_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "amhg")
Sac_rch2_NED_6hr_default_man <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "manning")

Sac_rch2_NED_6hr_default_amhg_val <- bam_val(Sac_rch2_NED_6hr_default_amhg, Sac_rch2_qobs_6hr)
Sac_rch2_NED_6hr_default_man_val <- bam_val(Sac_rch2_NED_6hr_default_man, Sac_rch2_qobs_6hr)
Sac_rch2_NED_6hr_default_man_amhg_val <- bam_val(Sac_rch2_NED_6hr_default_man_amhg, Sac_rch2_qobs_6hr)

Sac_rch2_data_3m_6hr <- bam_data(w = t(Sac_rch2_width_3m_6hr), s = Sac_rch2_slope_6hr, dA = Sac_rch2_dA_3m_6hr, 
                                 Qhat = Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "manning_amhg")
Sac_rch2_3m_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "amhg")
Sac_rch2_3m_6hr_default_man <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "manning")

Sac_rch2_3m_6hr_default_amhg_val <- bam_val(Sac_rch2_3m_6hr_default_amhg, Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_default_man_val <- bam_val(Sac_rch2_3m_6hr_default_man, Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_default_man_amhg_val <- bam_val(Sac_rch2_3m_6hr_default_man_amhg, Sac_rch2_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr Default NED",
     ylim = c(120, 275))
lines(Sac_rch2_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Sac_rch2_NED_6hr_default_man_val[[1]], col = "blue")
lines(Sac_rch2_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch2_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch2_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch2_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch2_data_NED_6hr <- bam_data(w = t(Sac_rch2_width_NED_6hr), s = Sac_rch2_slope_6hr, dA = Sac_rch2_dA_NED_6hr,
                                  Qhat = Sac_rch2_qobs_6hr)
Sac_best_priors <- bam_priors(bamdata= Sac_rch2_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Sac_rch2_qobs_6hr)),
                              upperbound_logQ = log(1.25 * max(Sac_rch2_qobs_6hr)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch2_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch2_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "manning_amhg", 
                                               bampriors = Sac_best_priors)
Sac_rch2_NED_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_NED_6hr_best_man <- bam_estimate(bamdata = Sac_rch2_data_NED_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_NED_6hr_best_amhg_val <- bam_val(Sac_rch2_NED_6hr_best_amhg, Sac_rch2_qobs_6hr)
Sac_rch2_NED_6hr_best_man_val <- bam_val(Sac_rch2_NED_6hr_best_man, Sac_rch2_qobs_6hr)
Sac_rch2_NED_6hr_best_man_amhg_val <- bam_val(Sac_rch2_NED_6hr_best_man_amhg, Sac_rch2_qobs_6hr)

Sac_rch2_data_3m_6hr <- bam_data(w = t(Sac_rch2_width_3m_6hr), s = Sac_rch2_slope_6hr, dA = Sac_rch2_dA_3m_6hr,
                                 Qhat = Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "manning_amhg",
                                              bampriors = Sac_best_priors)
Sac_rch2_3m_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_3m_6hr_best_man <- bam_estimate(bamdata = Sac_rch2_data_3m_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_3m_6hr_best_amhg_val <- bam_val(Sac_rch2_3m_6hr_best_amhg, Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_best_man_val <- bam_val(Sac_rch2_3m_6hr_best_man, Sac_rch2_qobs_6hr)
Sac_rch2_3m_6hr_best_man_amhg_val <- bam_val(Sac_rch2_3m_6hr_best_man_amhg, Sac_rch2_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr best NED",
     ylim = c(135, 155))
lines(Sac_rch2_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch2_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch2_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch2_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch2_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch2_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Sac_rch2_data_NED_hourly <- bam_data(w = t(Sac_rch2_width_NED_hourly), s = Sac_rch2_slope_hourly, dA = Sac_rch2_dA_NED_hourly,
                                     Qhat = Sac_rch2_qobs_hourly)
Sac_rch2_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "manning_amhg")
Sac_rch2_NED_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "amhg")
Sac_rch2_NED_hourly_default_man <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "manning")

Sac_rch2_NED_hourly_default_amhg_val <- bam_val(Sac_rch2_NED_hourly_default_amhg, Sac_rch2_qobs_hourly)
Sac_rch2_NED_hourly_default_man_val <- bam_val(Sac_rch2_NED_hourly_default_man, Sac_rch2_qobs_hourly)
Sac_rch2_NED_hourly_default_man_amhg_val <- bam_val(Sac_rch2_NED_hourly_default_man_amhg, Sac_rch2_qobs_hourly)

Sac_rch2_data_3m_hourly <- bam_data(w = t(Sac_rch2_width_3m_hourly), s = Sac_rch2_slope_hourly, dA = Sac_rch2_dA_3m_hourly, 
                                    Qhat = Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "manning_amhg")
Sac_rch2_3m_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "amhg")
Sac_rch2_3m_hourly_default_man <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "manning")

Sac_rch2_3m_hourly_default_amhg_val <- bam_val(Sac_rch2_3m_hourly_default_amhg, Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_default_man_val <- bam_val(Sac_rch2_3m_hourly_default_man, Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_default_man_amhg_val <- bam_val(Sac_rch2_3m_hourly_default_man_amhg, Sac_rch2_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly Default NED",
     ylim = c(120, 275))
lines(Sac_rch2_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Sac_rch2_NED_hourly_default_man_val[[1]], col = "blue")
lines(Sac_rch2_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch2_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch2_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch2_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch2_data_NED_hourly <- bam_data(w = t(Sac_rch2_width_NED_hourly), s = Sac_rch2_slope_hourly, dA = Sac_rch2_dA_NED_hourly,
                                     Qhat = Sac_rch2_qobs_hourly)
Sac_best_priors <- bam_priors(bamdata= Sac_rch2_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Sac_rch2_qobs_hourly)),
                              upperbound_logQ = log(1.25 * max(Sac_rch2_qobs_hourly)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch2_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch2_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "manning_amhg", 
                                                  bampriors = Sac_best_priors)
Sac_rch2_NED_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_NED_hourly_best_man <- bam_estimate(bamdata = Sac_rch2_data_NED_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_NED_hourly_best_amhg_val <- bam_val(Sac_rch2_NED_hourly_best_amhg, Sac_rch2_qobs_hourly)
Sac_rch2_NED_hourly_best_man_val <- bam_val(Sac_rch2_NED_hourly_best_man, Sac_rch2_qobs_hourly)
Sac_rch2_NED_hourly_best_man_amhg_val <- bam_val(Sac_rch2_NED_hourly_best_man_amhg, Sac_rch2_qobs_hourly)

Sac_rch2_data_3m_hourly <- bam_data(w = t(Sac_rch2_width_3m_hourly), s = Sac_rch2_slope_hourly, dA = Sac_rch2_dA_3m_hourly,
                                    Qhat = Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "manning_amhg",
                                                 bampriors = Sac_best_priors)
Sac_rch2_3m_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_3m_hourly_best_man <- bam_estimate(bamdata = Sac_rch2_data_3m_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_3m_hourly_best_amhg_val <- bam_val(Sac_rch2_3m_hourly_best_amhg, Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_best_man_val <- bam_val(Sac_rch2_3m_hourly_best_man, Sac_rch2_qobs_hourly)
Sac_rch2_3m_hourly_best_man_amhg_val <- bam_val(Sac_rch2_3m_hourly_best_man_amhg, Sac_rch2_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly best NED",
     ylim = c(135, 155))
lines(Sac_rch2_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch2_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch2_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch2_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch2_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch2_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#15min-----
#Default----
Sac_rch2_data_NED_15min <- bam_data(w = t(Sac_rch2_width_NED_15min), s = Sac_rch2_slope_15min, dA = Sac_rch2_dA_NED_15min,
                                    Qhat = Sac_rch2_qobs_15min)
Sac_rch2_NED_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "manning_amhg")
Sac_rch2_NED_15min_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "amhg")
Sac_rch2_NED_15min_default_man <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "manning")

library(hydroGOF)
Sac_rch2_NED_15min_default_amhg_val <- bam_val(Sac_rch2_NED_15min_default_amhg, Sac_rch2_qobs_15min)
Sac_rch2_NED_15min_default_man_val <- bam_val(Sac_rch2_NED_15min_default_man, Sac_rch2_qobs_15min)
Sac_rch2_NED_15min_default_man_amhg_val <- bam_val(Sac_rch2_NED_15min_default_man_amhg, Sac_rch2_qobs_15min)

Sac_rch2_data_3m_15min <- bam_data(w = t(Sac_rch2_width_3m_15min), s = Sac_rch2_slope_15min, dA = Sac_rch2_dA_3m_15min, 
                                   Qhat = Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "manning_amhg")
Sac_rch2_3m_15min_default_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "amhg")
Sac_rch2_3m_15min_default_man <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "manning")

library(hydroGOF)
Sac_rch2_3m_15min_default_amhg_val <- bam_val(Sac_rch2_3m_15min_default_amhg, Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_default_man_val <- bam_val(Sac_rch2_3m_15min_default_man, Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_default_man_amhg_val <- bam_val(Sac_rch2_3m_15min_default_man_amhg, Sac_rch2_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min Default NED",
     ylim = c(120, 275))
lines(Sac_rch2_NED_15min_default_amhg_val[[1]], col = "red")
lines(Sac_rch2_NED_15min_default_man_val[[1]], col = "blue")
lines(Sac_rch2_NED_15min_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch2_3m_15min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch2_3m_15min_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch2_3m_15min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch2_data_NED_15min <- bam_data(w = t(Sac_rch2_width_NED_15min), s = Sac_rch2_slope_15min, dA = Sac_rch2_dA_NED_15min,
                                    Qhat = Sac_rch2_qobs_15min)
Sac_best_priors <- bam_priors(bamdata= Sac_rch2_data_NED_15min, lowerbound_logQ = log(0.75 * min(Sac_rch2_qobs_15min)),
                              upperbound_logQ = log(1.25 * max(Sac_rch2_qobs_15min)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch2_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch2_NED_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch2_NED_15min_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_NED_15min_best_man <- bam_estimate(bamdata = Sac_rch2_data_NED_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_NED_15min_best_amhg_val <- bam_val(Sac_rch2_NED_15min_best_amhg, Sac_rch2_qobs_15min)
Sac_rch2_NED_15min_best_man_val <- bam_val(Sac_rch2_NED_15min_best_man, Sac_rch2_qobs_15min)
Sac_rch2_NED_15min_best_man_amhg_val <- bam_val(Sac_rch2_NED_15min_best_man_amhg, Sac_rch2_qobs_15min)

Sac_rch2_data_3m_15min <- bam_data(w = t(Sac_rch2_width_3m_15min), s = Sac_rch2_slope_15min, dA = Sac_rch2_dA_3m_15min,
                                   Qhat = Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch2_3m_15min_best_amhg <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch2_3m_15min_best_man <- bam_estimate(bamdata = Sac_rch2_data_3m_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch2_3m_15min_best_amhg_val <- bam_val(Sac_rch2_3m_15min_best_amhg, Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_best_man_val <- bam_val(Sac_rch2_3m_15min_best_man, Sac_rch2_qobs_15min)
Sac_rch2_3m_15min_best_man_amhg_val <- bam_val(Sac_rch2_3m_15min_best_man_amhg, Sac_rch2_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch2_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min best NED",
     ylim = c(135, 155))
lines(Sac_rch2_NED_15min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch2_NED_15min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch2_NED_15min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch2_3m_15min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch2_3m_15min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch2_3m_15min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Sac_rch2_daily_lower_logQ <- log(seq(exp(max(apply(log(Sac_rch2_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                     min(Sac_rch2_qobs_15min),length.out = 10))
Sac_rch2_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_logQ = Sac_rch2_daily_lower_logQ[i])
}
#upperbound_logQ
Sac_rch2_daily_upper_logQ <- log(seq(max(Sac_rch2_qobs_15min), exp(min(apply(log(Sac_rch2_width_NED_daily), 2, max)) + log(40) + log(5)),
                                     length.out = 10))
Sac_rch2_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_logQ = Sac_rch2_daily_upper_logQ[i])
}
#lowerbound_A0
Sac_rch2_daily_lower_A0 <- seq(0, 30, length.out = 10)
Sac_rch2_daily_lower_A0_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_A0 = Sac_rch2_daily_lower_A0[i])
}
#upperbound_A0
Sac_rch2_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Sac_rch2_daily_upper_A0_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_A0 = Sac_rch2_daily_upper_A0[i])
}
#lowerbound_logn
Sac_rch2_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Sac_rch2_daily_lower_logn_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_logn = Sac_rch2_daily_lower_logn[i])
}
#upperbound_logn
Sac_rch2_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Sac_rch2_daily_upper_logn_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_logn = Sac_rch2_daily_upper_logn[i])
}
#lowerbound_logQc
Sac_rch2_daily_lower_logQc <- seq(0, log(min(Sac_rch2_qobs_15min)),length.out = 10)
Sac_rch2_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_logQc = Sac_rch2_daily_lower_logQc[i])
}
#upperbound_logQc
Sac_rch2_daily_upper_logQc <- seq(log(max(Sac_rch2_qobs_15min)), 10, length.out = 10)
Sac_rch2_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_logQc = Sac_rch2_daily_upper_logQc[i])
}
#lowerbound_logWc
Sac_rch2_daily_lower_logWc <- seq(1, log(min(Sac_rch2_width_NED_daily)),length.out = 10)
Sac_rch2_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_logWc = Sac_rch2_daily_lower_logWc[i])
}
#upperbound_logWc
Sac_rch2_daily_upper_logWc <- seq(log(max(Sac_rch2_width_NED_daily)), 8, length.out = 10)
Sac_rch2_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_logWc = Sac_rch2_daily_upper_logWc[i])
}
#lowerbound_b
Sac_rch2_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Sac_rch2_daily_lower_b_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_lower_b_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, lowerbound_b = Sac_rch2_daily_lower_b[i])
}
#upperbound_b
Sac_rch2_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Sac_rch2_daily_upper_b_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_upper_b_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, upperbound_b = Sac_rch2_daily_upper_b[i])
}
#sigma_man
Sac_rch2_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Sac_rch2_daily_sigma_man_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, sigma_man = Sac_rch2_daily_sigma_man[i])
}
#sigma_amhg
Sac_rch2_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Sac_rch2_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, sigma_amhg = Sac_rch2_daily_sigma_amhg[i])
}
#logQc_hat
Sac_rch2_daily_logQc_hat <- log(seq(min(Sac_rch2_qobs_15min), max(Sac_rch2_qobs_15min), length.out = 10))
Sac_rch2_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logQc_hat = Sac_rch2_daily_logQc_hat[i])
}
#logWc_hat
Sac_rch2_daily_logWc_hat <- log(seq(min(Sac_rch2_width_NED_daily), max(Sac_rch2_width_NED_daily), length.out = 10))
Sac_rch2_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logWc_hat = Sac_rch2_daily_logWc_hat[i])
}
#b_hat
Sac_rch2_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Sac_rch2_daily_b_hat_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_b_hat_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, b_hat = Sac_rch2_daily_b_hat[i])
}
#logA0_hat
Sac_rch2_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Sac_rch2_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logA0_hat = Sac_rch2_daily_logA0_hat[i])
}
#logn_hat
Sac_rch2_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Sac_rch2_daily_logn_hat_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logn_hat = Sac_rch2_daily_logn_hat[i])
}
#logQ_sd
Sac_rch2_daily_logQ_sd <- seq(sd(log(Sac_rch2_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch2_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logQ_sd = Sac_rch2_daily_logQ_sd[i])
}
#logQc_sd
Sac_rch2_daily_logQc_sd <- seq(sd(log(Sac_rch2_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch2_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logQc_sd = Sac_rch2_daily_logQc_sd[i])
}
#logWc_Sd
Sac_rch2_daily_logWc_Sd <- seq(log(sd(Sac_rch2_width_NED_daily)), 4.712493, length.out = 10)
Sac_rch2_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logWc_sd = Sac_rch2_daily_logWc_Sd[i])
}
#b_sd
Sac_rch2_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Sac_rch2_daily_b_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_b_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, b_sd = Sac_rch2_daily_b_sd[i])
}
#logA0_sd
Sac_rch2_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Sac_rch2_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logA0_sd = Sac_rch2_daily_logA0_sd[i])
}
#logn_sd
Sac_rch2_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Sac_rch2_daily_logn_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, logn_sd = Sac_rch2_daily_logn_sd[i])
}
#Werr_sd
Sac_rch2_daily_Werr_sd <- seq(5, 15, length.out = 10)
Sac_rch2_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, Werr_sd = Sac_rch2_daily_Werr_sd[i])
}
#Serr_sd
Sac_rch2_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Sac_rch2_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, Serr_sd = Sac_rch2_daily_Serr_sd[i])
}
#dAerr_sd
Sac_rch2_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Sac_rch2_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Sac_rch2_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch2_data_NED_daily, dAerr_sd = Sac_rch2_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Sac_rch2_daily_manning_list <- list(Sac_rch2_daily_lower_logQ_priors,
                                    Sac_rch2_daily_upper_logQ_priors,
                                    Sac_rch2_daily_lower_A0_priors,
                                    Sac_rch2_daily_upper_A0_priors,
                                    Sac_rch2_daily_lower_logn_priors,
                                    Sac_rch2_daily_upper_logn_priors,
                                    Sac_rch2_daily_sigma_man_priors,
                                    Sac_rch2_daily_logA0_hat_priors,
                                    Sac_rch2_daily_logn_hat_priors,
                                    Sac_rch2_daily_logQ_sd_priors,
                                    Sac_rch2_daily_logA0_sd_priors,
                                    Sac_rch2_daily_logn_sd_priors,
                                    Sac_rch2_daily_Werr_sd_priors,
                                    Sac_rch2_daily_Serr_sd_priors,
                                    Sac_rch2_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Sac_rch2_data_NED_daily, variant = "manning", bampriors = Sac_rch2_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Sac_rch2_daily_amhg_list <- list(Sac_rch2_daily_lower_logQ_priors,
                                 Sac_rch2_daily_upper_logQ_priors,
                                 Sac_rch2_daily_lower_logQc_priors,
                                 Sac_rch2_daily_upper_logQc_priors,
                                 Sac_rch2_daily_lower_logWc_priors,
                                 Sac_rch2_daily_upper_logWc_priors,
                                 Sac_rch2_daily_lower_b_priors,
                                 Sac_rch2_daily_upper_b_priors,
                                 Sac_rch2_daily_sigma_amhg_priors,
                                 Sac_rch2_daily_logQc_hat_priors,
                                 Sac_rch2_daily_logWc_hat_priors,
                                 Sac_rch2_daily_b_hat_priors,
                                 Sac_rch2_daily_logQ_sd_priors,
                                 Sac_rch2_daily_logQc_sd_priors,
                                 Sac_rch2_daily_logWc_Sd_priors,
                                 Sac_rch2_daily_b_sd_priors,
                                 Sac_rch2_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Sac_rch2_data_NED_daily, variant = "amhg", bampriors = Sac_rch2_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Sac_rch2_daily_man_amhg_list <- list(Sac_rch2_daily_lower_logQ_priors,
                                     Sac_rch2_daily_upper_logQ_priors,
                                     Sac_rch2_daily_lower_A0_priors,
                                     Sac_rch2_daily_upper_A0_priors,
                                     Sac_rch2_daily_lower_logn_priors,
                                     Sac_rch2_daily_upper_logn_priors,
                                     Sac_rch2_daily_lower_logQc_priors,
                                     Sac_rch2_daily_upper_logQc_priors,
                                     Sac_rch2_daily_lower_logWc_priors,
                                     Sac_rch2_daily_upper_logWc_priors,
                                     Sac_rch2_daily_lower_b_priors,
                                     Sac_rch2_daily_upper_b_priors,
                                     Sac_rch2_daily_sigma_man_priors,
                                     Sac_rch2_daily_sigma_amhg_priors,
                                     Sac_rch2_daily_logQc_hat_priors,
                                     Sac_rch2_daily_logWc_hat_priors,
                                     Sac_rch2_daily_b_hat_priors,
                                     Sac_rch2_daily_logA0_hat_priors,
                                     Sac_rch2_daily_logn_hat_priors,
                                     Sac_rch2_daily_logQ_sd_priors,
                                     Sac_rch2_daily_logQc_sd_priors,
                                     Sac_rch2_daily_logWc_Sd_priors,
                                     Sac_rch2_daily_b_sd_priors,
                                     Sac_rch2_daily_logA0_sd_priors,
                                     Sac_rch2_daily_logn_sd_priors,
                                     Sac_rch2_daily_Werr_sd_priors,
                                     Sac_rch2_daily_Serr_sd_priors,
                                     Sac_rch2_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_2/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Sac_rch2_data_NED_daily, variant = "manning_amhg", bampriors = Sac_rch2_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

###Sacramento_Reach_Three
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Sacramento/Reach_3"))
load("Sac_WSE_rch3_daily.RData")
load("Sac_WSE_rch3_6hr.RData")
load("Sac_WSE_rch3_hourly.RData")
load("Sac_WSE_rch3_15min.RData")

Sac_rch3_PT_coords <- read.csv("Sac_rch3_PT_coords.csv")
Sac_rch3_PT_coords_sp <- SpatialPoints(Sac_rch3_PT_coords[,c(5,4)], proj4string =  CRS("+proj=utm +zone=10 +datum=WGS84"))
Sac_cl_rch3 <-readOGR(dsn = ".", layer = "Sac_rch3_Lid_cl")

setwd(paste0(home_dir, "/Sacramento/Reach_3/DEM"))
Sac_rch3_DEM_3m <- raster("Sac_DEM_rch3_3m.tif")
Sac_rch3_DEM_NED <- raster("Sac_DEM_rch3.tif")

setwd(paste0(home_dir, "/Sacramento/Reach_3"))
Sac_orthos_rch3_NED <- ortho_lines(Sac_rch3_DEM_NED, data.frame(Sac_cl_rch3@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch3_dist_mat_NED <- pointDistance(Sac_rch3_PT_coords_sp, data.frame(Sac_cl_rch3@lines[[1]]@Lines[[1]]@coords))
Sac_rch3_closest_ortho_NED <- Sac_orthos_rch3_NED[apply(Sac_rch3_dist_mat_NED, 1, FUN = which.min)]

Sac_orthos_rch3_3m <- ortho_lines(Sac_rch3_DEM_3m, data.frame(Sac_cl_rch3@lines[[1]]@Lines[[1]]@coords), 500)
Sac_rch3_dist_mat_3m <- pointDistance(Sac_rch3_PT_coords_sp, data.frame(Sac_cl_rch3@lines[[1]]@Lines[[1]]@coords))
Sac_rch3_closest_ortho_3m <- Sac_orthos_rch3_3m[apply(Sac_rch3_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Sac_rch3_elev_prof_NED <- foreach(i=1:length(Sac_rch3_closest_ortho_NED), 
                                  .combine='rbind', .packages=c('raster')) %dopar% {
                                    extract(Sac_rch3_DEM_NED, Sac_rch3_closest_ortho_NED[[i]], method = 'simple')
                                  }
Sac_rch3_elev_prof_3m <- foreach(i=1:length(Sac_rch3_closest_ortho_3m), 
                                 .combine='rbind', .packages=c('raster')) %dopar% {
                                   extract(Sac_rch3_DEM_3m, Sac_rch3_closest_ortho_3m[[i]], method = 'simple')
                                 }
proc.time()-ptm
stopCluster(cl)

Sac_rch3_scale_x_NED <- vector()
for(i in 1:length(Sac_rch3_elev_prof_NED)){
  Sac_rch3_scale_x_NED[i] = 1000/length(Sac_rch3_elev_prof_NED[[i]])
}

Sac_rch3_scale_x_3m <- vector()
for(i in 1:length(Sac_rch3_elev_prof_3m)){
  Sac_rch3_scale_x_3m[i] = 1000/length(Sac_rch3_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Sac_rch3_elev_prof_NED, Sac_WSE_rch3_daily_df, "Sacramento Reach Three")
locator(n = 4)
Sac_rch3_upper_l_NED <- c( 82.6, 50.2,  56.2, 67.7, 68.8, 46.8, 53.5, 64.2,  53.3, 46.9, 42.6, 49.3, 49.3, 65.7, 66.6, 66.2, 56.4, 46.7, 37.4)
Sac_rch3_lower_l_NED <- c( 88.9, 54.6,  59.8, 72.9, 76.0, 53.6, 56.4, 70.1,  57.6, 49.3, 45.3, 53.7, 53.7, 69.2, 70.5, 72.7, 59.9, 50.4, 41.0)
Sac_rch3_lower_r_NED <- c(116.4, 91.8, 100.8, 87.7, 88.2, 68.3, 85.2, 84.3, 107.9, 88.8, 74.1, 66.6, 66.6, 77.5, 86.0, 88.5, 72.5, 66.0, 59.5)
Sac_rch3_upper_r_NED <- c(125.6, 95.3, 104.0, 93.4, 91.8, 71.3, 88.3, 88.9, 111.0, 92.2, 77.4, 70.8, 70.8, 81.0, 90.0, 95.7, 77.5, 69.7, 63.3)

Sac_rch3_elev_prof_NED_corr <- Sac_rch3_elev_prof_NED
Sac_rch3_elev_prof_NED_corr[[12]] <- Sac_rch3_elev_prof_NED[[13]]
Sac_rch3_scale_x_NED[12] <- Sac_rch3_scale_x_NED[13]

Sac_rch3_width_NED_daily <- calc_width(h = t(Sac_WSE_rch3_daily_df), elev_prof = Sac_rch3_elev_prof_NED_corr, 
                                       upper_l = Sac_rch3_upper_l_NED, lower_l = Sac_rch3_lower_l_NED,
                                       lower_r = Sac_rch3_lower_r_NED, upper_r = Sac_rch3_upper_r_NED,
                                       scale_prof = Sac_rch3_scale_x_NED)
Sac_rch3_width_NED_daily <- na.approx(Sac_rch3_width_NED_daily)
Sac_rch3_dA_NED_daily <- calcdA_mat(w = t(Sac_rch3_width_NED_daily), h = Sac_WSE_rch3_daily_df)

Sac_rch3_width_NED_6hr <- calc_width(h = t(Sac_WSE_rch3_6hr_df), elev_prof = Sac_rch3_elev_prof_NED_corr, 
                                     upper_l = Sac_rch3_upper_l_NED, lower_l = Sac_rch3_lower_l_NED,
                                     lower_r = Sac_rch3_lower_r_NED, upper_r = Sac_rch3_upper_r_NED,
                                     scale_prof = Sac_rch3_scale_x_NED)
Sac_rch3_width_NED_6hr <- na.approx(Sac_rch3_width_NED_6hr)
Sac_rch3_dA_NED_6hr <- calcdA_mat(w = t(Sac_rch3_width_NED_6hr), h = Sac_WSE_rch3_6hr_df)

Sac_rch3_width_NED_hourly <- calc_width(h = t(Sac_WSE_rch3_hourly_df), elev_prof = Sac_rch3_elev_prof_NED_corr, 
                                        upper_l = Sac_rch3_upper_l_NED, lower_l = Sac_rch3_lower_l_NED,
                                        lower_r = Sac_rch3_lower_r_NED, upper_r = Sac_rch3_upper_r_NED,
                                        scale_prof = Sac_rch3_scale_x_NED)
Sac_rch3_width_NED_hourly <- na.approx(Sac_rch3_width_NED_hourly)
Sac_rch3_dA_NED_hourly <- calcdA_mat(w = t(Sac_rch3_width_NED_hourly), h = Sac_WSE_rch3_hourly_df)

Sac_rch3_width_NED_15min <- calc_width(h = t(Sac_WSE_rch3_df), elev_prof = Sac_rch3_elev_prof_NED_corr, 
                                       upper_l = Sac_rch3_upper_l_NED, lower_l = Sac_rch3_lower_l_NED,
                                       lower_r = Sac_rch3_lower_r_NED, upper_r = Sac_rch3_upper_r_NED,
                                       scale_prof = Sac_rch3_scale_x_NED)
Sac_rch3_width_NED_15min <- na.approx(Sac_rch3_width_NED_15min)
Sac_rch3_dA_NED_15min <- calcdA_mat(w = t(Sac_rch3_width_NED_15min), h = Sac_WSE_rch3_df)

plot_elev_prof(Sac_rch3_elev_prof_3m, Sac_WSE_rch3_daily_df, "Sacramento Reach Three")
locator(n = 4)
Sac_rch3_upper_l_3m <- c(157.7, 169.0, 117.5, 201.2, 143.3, 129.0, 143.9,  33.0, 132.9, 230.4, 182.2, 211.3,  88.1, 155.0,  63.5,  90.4, 197.4, 203.8, 106.4)
Sac_rch3_lower_l_3m <- c(175.8, 176.1, 123.6, 216.4, 156.9, 136.9, 150.7,  54.5, 159.9, 237.0, 191.0, 226.7,  94.7, 168.1,  72.7, 103.3, 204.3, 209.4, 130.2)
Sac_rch3_lower_r_3m <- c(209.0, 218.0, 216.5, 303.7, 353.8, 343.5, 209.0, 455.0, 390.5, 376.9, 262.4, 309.3, 178.2, 275.2, 179.1, 160.2, 253.8, 242.9, 262.4)
Sac_rch3_upper_r_3m <- c(231.2, 226.1, 225.0, 321.9, 371.9, 356.9, 225.0, 468.1, 401.4, 388.1, 283.4, 323.1, 184.8, 285.0, 207.6, 175.4, 266.2, 252.2, 275.2)

Sac_rch3_width_3m_daily <- calc_width(h = t(Sac_WSE_rch3_daily_df), elev_prof = Sac_rch3_elev_prof_3m, 
                                      upper_l = Sac_rch3_upper_l_3m, lower_l = Sac_rch3_lower_l_3m,
                                      lower_r = Sac_rch3_lower_r_3m, upper_r = Sac_rch3_upper_r_3m,
                                      scale_prof = Sac_rch3_scale_x_3m)
Sac_rch3_width_3m_daily <- na.approx(Sac_rch3_width_3m_daily)
Sac_rch3_dA_3m_daily <- calcdA_mat(w = t(Sac_rch3_width_3m_daily), h = Sac_WSE_rch3_daily_df)

Sac_rch3_width_3m_6hr <- calc_width(h = t(Sac_WSE_rch3_6hr_df), elev_prof = Sac_rch3_elev_prof_3m, 
                                    upper_l = Sac_rch3_upper_l_3m, lower_l = Sac_rch3_lower_l_3m,
                                    lower_r = Sac_rch3_lower_r_3m, upper_r = Sac_rch3_upper_r_3m,
                                    scale_prof = Sac_rch3_scale_x_3m)
Sac_rch3_width_3m_6hr <- na.approx(Sac_rch3_width_3m_6hr)
Sac_rch3_dA_3m_6hr <- calcdA_mat(w = t(Sac_rch3_width_3m_6hr), h = Sac_WSE_rch3_6hr_df)

Sac_rch3_width_3m_hourly <- calc_width(h = t(Sac_WSE_rch3_hourly_df), elev_prof = Sac_rch3_elev_prof_3m, 
                                       upper_l = Sac_rch3_upper_l_3m, lower_l = Sac_rch3_lower_l_3m,
                                       lower_r = Sac_rch3_lower_r_3m, upper_r = Sac_rch3_upper_r_3m,
                                       scale_prof = Sac_rch3_scale_x_3m)
Sac_rch3_width_3m_hourly <- na.approx(Sac_rch3_width_3m_hourly)
Sac_rch3_dA_3m_hourly <- calcdA_mat(w = t(Sac_rch3_width_3m_hourly), h = Sac_WSE_rch3_hourly_df)

Sac_rch3_width_3m_15min <- calc_width(h = t(Sac_WSE_rch3_df), elev_prof = Sac_rch3_elev_prof_3m, 
                                      upper_l = Sac_rch3_upper_l_3m, lower_l = Sac_rch3_lower_l_3m,
                                      lower_r = Sac_rch3_lower_r_3m, upper_r = Sac_rch3_upper_r_3m,
                                      scale_prof = Sac_rch3_scale_x_3m)
Sac_rch3_width_3m_15min <- na.approx(Sac_rch3_width_3m_15min)
Sac_rch3_dA_3m_15min <- calcdA_mat(w = t(Sac_rch3_width_3m_15min), h = Sac_WSE_rch3_df)


xvec_rch3 <- calc_xvec(Sac_rch3_PT_coords_sp, Sac_cl_rch3)
Sac_rch3_xvec <- saveRDS(xvec_rch3, file = "Sac_rch3_xvec.rds")
Sac_rch3_slope_daily <- calcslope(rev(xvec_rch3), hmat = Sac_WSE_rch3_daily_df)
Sac_rch3_slope_daily <- na.replace(Sac_rch3_slope_daily, mean(Sac_rch3_slope_daily))
Sac_rch3_slope_daily[1,] <- colMeans(Sac_rch3_slope_daily, na.rm = TRUE)
Sac_rch3_slope_6hr <- calcslope(rev(xvec_rch3), hmat = Sac_WSE_rch3_6hr_df)
Sac_rch3_slope_6hr <- na.approx(Sac_rch3_slope_6hr)
Sac_rch3_slope_6hr[1,] <- colMeans(Sac_rch3_slope_6hr, na.rm = TRUE)
Sac_rch3_slope_hourly <- calcslope(rev(xvec_rch3), hmat = Sac_WSE_rch3_hourly_df)
Sac_rch3_slope_hourly <- na.approx(Sac_rch3_slope_hourly)
Sac_rch3_slope_hourly[1,] <- colMeans(Sac_rch3_slope_hourly, na.rm = TRUE)
Sac_rch3_slope_15min <- calcslope(rev(xvec_rch3), hmat = Sac_WSE_rch3_df)
Sac_rch3_slope_15min <- na.approx(Sac_rch3_slope_15min)
Sac_rch3_slope_15min[1,] <- colMeans(Sac_rch3_slope_15min, na.rm = TRUE)

setwd(paste0(home_dir, "/Sacramento/OrigData"))
Sacqobs <- read.csv("SacQobs.csv")$Q*0.028316847
Sac_rch3_qobs_daily<- colMeans(matrix(Sacqobs[1825:7872], 96))
Sac_rch3_qobs_6hr<- colMeans(matrix(Sacqobs[1825:7920], 24))
Sac_rch3_qobs_hourly<- colMeans(matrix(Sacqobs[1825:7936], 4))
Sac_rch3_qobs_15min<- colMeans(matrix(Sacqobs[1261:7938], 1))

save(Sac_rch3_width_NED_daily, file = "Sacramento_Reach_Three_width_NED_daily.RData")
save(Sac_rch3_dA_NED_daily, file = "Sacramento_Reach_Three_dA_NED_daily.RData")
save(Sac_rch3_slope_daily, file = "Sacramento_Reach_Three_slope_daily.RData")
save(Sac_rch3_qobs_daily, file = "Sacramento_Reach_Three_qobs_daily.RData")
save(Sac_rch3_width_NED_6hr, file = "Sacramento_Reach_Three_width_NED_6hr.RData")
save(Sac_rch3_dA_NED_6hr, file = "Sacramento_Reach_Three_dA_NED_6hr.RData")
save(Sac_rch3_slope_6hr, file = "Sacramento_Reach_Three_slope_6hr.RData")
save(Sac_rch3_qobs_6hr, file = "Sacramento_Reach_Three_qobs_6hr.RData")
save(Sac_rch3_width_NED_hourly, file = "Sacramento_Reach_Three_width_NED_hourly.RData")
save(Sac_rch3_dA_NED_hourly, file = "Sacramento_Reach_Three_dA_NED_hourly.RData")
save(Sac_rch3_slope_hourly, file = "Sacramento_Reach_Three_slope_hourly.RData")
save(Sac_rch3_qobs_hourly, file = "Sacramento_Reach_Three_qobs_hourly.RData")
save(Sac_rch3_width_NED_15min, file = "Sacramento_Reach_Three_width_NED_15min.RData")
save(Sac_rch3_dA_NED_15min, file = "Sacramento_Reach_Three_dA_NED_15min.RData")
save(Sac_rch3_slope_15min, file = "Sacramento_Reach_Three_slope_15min.RData")
save(Sac_rch3_qobs_15min, file = "Sacramento_Reach_Three_qobs_15min.RData")




#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Sac_rch3_wse_plot <- matplot(t(Sac_WSE_rch3_daily_df), type = c("l"), col= colfunc,
                             main = "WSE", xlab = "Days",
                             ylab = "WSE (m)",
                             cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch3_slope_plot <- matplot(t(Sac_rch3_slope_daily), type = c("l"), col= colfunc,
                               main = "Slope", xlab = "Days",
                               ylab = "Slope",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch3_width_plot <- matplot(Sac_rch3_width_NED_daily, type = c("l"), col= colfunc,
                               main = "Width (NED)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch3_width_plot <- matplot(Sac_rch3_width_3m_daily, type = c("l"), col= colfunc,
                               main = "Width (Lidar)", xlab = "Days",
                               ylab = "Width (m)",
                               cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch3_dA_plot <- matplot(t(Sac_rch3_dA_NED_daily), type = c("l"), col= colfunc,
                            main = "dA (NED)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Sac_rch3_dA_plot <- matplot(t(Sac_rch3_dA_3m_daily), type = c("l"), col= colfunc,
                            main = "dA (Lidar)", xlab = "Days",
                            ylab = "dA (m^2)",
                            cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
Sac_rch3_data_NED_daily <- bam_data(w = t(Sac_rch3_width_NED_daily), s = Sac_rch3_slope_daily, dA = Sac_rch3_dA_NED_daily,
                                    Qhat = Sac_rch3_qobs_daily)
Sac_rch3_NED_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "manning_amhg")
Sac_rch3_NED_daily_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "amhg")
Sac_rch3_NED_daily_default_man <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "manning")
library(hydroGOF)
Sac_rch3_NED_daily_default_amhg_val <- bam_val(Sac_rch3_NED_daily_default_amhg, Sac_rch3_qobs_daily)
Sac_rch3_NED_daily_default_man_val <- bam_val(Sac_rch3_NED_daily_default_man, Sac_rch3_qobs_daily)
Sac_rch3_NED_daily_default_man_amhg_val <- bam_val(Sac_rch3_NED_daily_default_man_amhg, Sac_rch3_qobs_daily)

Sac_rch3_data_3m_daily <- bam_data(w = t(Sac_rch3_width_3m_daily), s = Sac_rch3_slope_daily, dA = Sac_rch3_dA_3m_daily,
                                   Qhat = Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "manning_amhg")
Sac_rch3_3m_daily_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "amhg")
Sac_rch3_3m_daily_default_man <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "manning")

Sac_rch3_3m_daily_default_amhg_val <- bam_val(Sac_rch3_3m_daily_default_amhg, Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_default_man_val <- bam_val(Sac_rch3_3m_daily_default_man, Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_default_man_amhg_val <- bam_val(Sac_rch3_3m_daily_default_man_amhg, Sac_rch3_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily Default NED",
     ylim = c(100, 220))
lines(Sac_rch3_NED_daily_default_amhg_val[[1]], col = "red")
lines(Sac_rch3_NED_daily_default_man_val[[1]], col = "blue")
lines(Sac_rch3_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch3_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch3_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch3_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch3_data_NED_daily <- bam_data(w = t(Sac_rch3_width_NED_daily), s = Sac_rch3_slope_daily, dA = Sac_rch3_dA_NED_daily,
                                    Qhat = Sac_rch3_qobs_daily)
Sac_best_priors <- bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_logQ = log(0.5 * min(Sac_rch3_qobs_daily)),
                              upperbound_logQ = log(1.1*max(Sac_rch3_qobs_daily)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch3_slope_daily, na.rm = TRUE),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch3_NED_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch3_NED_daily_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_NED_daily_best_man <- bam_estimate(bamdata = Sac_rch3_data_NED_daily,variant = "manning", bampriors = Sac_best_priors)
library(hydroGOF)
Sac_rch3_NED_daily_best_amhg_val <- bam_val(Sac_rch3_NED_daily_best_amhg, Sac_rch3_qobs_daily)
Sac_rch3_NED_daily_best_man_val <- bam_val(Sac_rch3_NED_daily_best_man, Sac_rch3_qobs_daily)
Sac_rch3_NED_daily_best_man_amhg_val <- bam_val(Sac_rch3_NED_daily_best_man_amhg, Sac_rch3_qobs_daily)

Sac_rch3_data_3m_daily <- bam_data(w = t(Sac_rch3_width_3m_daily), s = Sac_rch3_slope_daily, dA = Sac_rch3_dA_3m_daily,
                                   Qhat = Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch3_3m_daily_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_3m_daily_best_man <- bam_estimate(bamdata = Sac_rch3_data_3m_daily,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_3m_daily_best_amhg_val <- bam_val(Sac_rch3_3m_daily_best_amhg, Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_best_man_val <- bam_val(Sac_rch3_3m_daily_best_man, Sac_rch3_qobs_daily)
Sac_rch3_3m_daily_best_man_amhg_val <- bam_val(Sac_rch3_3m_daily_best_man_amhg, Sac_rch3_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento Daily best NED",
     ylim = c(100, 200))
lines(Sac_rch3_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch3_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch3_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch3_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch3_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch3_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
Sac_rch3_data_NED_6hr <- bam_data(w = t(Sac_rch3_width_NED_6hr), s = Sac_rch3_slope_6hr, dA = Sac_rch3_dA_NED_6hr,
                                  Qhat = Sac_rch3_qobs_6hr)
Sac_rch3_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "manning_amhg")
Sac_rch3_NED_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "amhg")
Sac_rch3_NED_6hr_default_man <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "manning")

Sac_rch3_NED_6hr_default_amhg_val <- bam_val(Sac_rch3_NED_6hr_default_amhg, Sac_rch3_qobs_6hr)
Sac_rch3_NED_6hr_default_man_val <- bam_val(Sac_rch3_NED_6hr_default_man, Sac_rch3_qobs_6hr)
Sac_rch3_NED_6hr_default_man_amhg_val <- bam_val(Sac_rch3_NED_6hr_default_man_amhg, Sac_rch3_qobs_6hr)

Sac_rch3_data_3m_6hr <- bam_data(w = t(Sac_rch3_width_3m_6hr), s = Sac_rch3_slope_6hr, dA = Sac_rch3_dA_3m_6hr, 
                                 Qhat = Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "manning_amhg")
Sac_rch3_3m_6hr_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "amhg")
Sac_rch3_3m_6hr_default_man <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "manning")

Sac_rch3_3m_6hr_default_amhg_val <- bam_val(Sac_rch3_3m_6hr_default_amhg, Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_default_man_val <- bam_val(Sac_rch3_3m_6hr_default_man, Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_default_man_amhg_val <- bam_val(Sac_rch3_3m_6hr_default_man_amhg, Sac_rch3_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr Default NED",
     ylim = c(120, 275))
lines(Sac_rch3_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Sac_rch3_NED_6hr_default_man_val[[1]], col = "blue")
lines(Sac_rch3_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch3_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch3_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch3_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch3_data_NED_6hr <- bam_data(w = t(Sac_rch3_width_NED_6hr), s = Sac_rch3_slope_6hr, dA = Sac_rch3_dA_NED_6hr,
                                  Qhat = Sac_rch3_qobs_6hr)
Sac_best_priors <- bam_priors(bamdata= Sac_rch3_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Sac_rch3_qobs_6hr)),
                              upperbound_logQ = log(1.25 * max(Sac_rch3_qobs_6hr)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch3_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch3_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "manning_amhg", 
                                               bampriors = Sac_best_priors)
Sac_rch3_NED_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_NED_6hr_best_man <- bam_estimate(bamdata = Sac_rch3_data_NED_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_NED_6hr_best_amhg_val <- bam_val(Sac_rch3_NED_6hr_best_amhg, Sac_rch3_qobs_6hr)
Sac_rch3_NED_6hr_best_man_val <- bam_val(Sac_rch3_NED_6hr_best_man, Sac_rch3_qobs_6hr)
Sac_rch3_NED_6hr_best_man_amhg_val <- bam_val(Sac_rch3_NED_6hr_best_man_amhg, Sac_rch3_qobs_6hr)

Sac_rch3_data_3m_6hr <- bam_data(w = t(Sac_rch3_width_3m_6hr), s = Sac_rch3_slope_6hr, dA = Sac_rch3_dA_3m_6hr,
                                 Qhat = Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "manning_amhg",
                                              bampriors = Sac_best_priors)
Sac_rch3_3m_6hr_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_3m_6hr_best_man <- bam_estimate(bamdata = Sac_rch3_data_3m_6hr,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_3m_6hr_best_amhg_val <- bam_val(Sac_rch3_3m_6hr_best_amhg, Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_best_man_val <- bam_val(Sac_rch3_3m_6hr_best_man, Sac_rch3_qobs_6hr)
Sac_rch3_3m_6hr_best_man_amhg_val <- bam_val(Sac_rch3_3m_6hr_best_man_amhg, Sac_rch3_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 6hr best NED",
     ylim = c(135, 155))
lines(Sac_rch3_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch3_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch3_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch3_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch3_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch3_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Sac_rch3_data_NED_hourly <- bam_data(w = t(Sac_rch3_width_NED_hourly), s = Sac_rch3_slope_hourly, dA = Sac_rch3_dA_NED_hourly,
                                     Qhat = Sac_rch3_qobs_hourly)
Sac_rch3_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "manning_amhg")
Sac_rch3_NED_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "amhg")
Sac_rch3_NED_hourly_default_man <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "manning")

Sac_rch3_NED_hourly_default_amhg_val <- bam_val(Sac_rch3_NED_hourly_default_amhg, Sac_rch3_qobs_hourly)
Sac_rch3_NED_hourly_default_man_val <- bam_val(Sac_rch3_NED_hourly_default_man, Sac_rch3_qobs_hourly)
Sac_rch3_NED_hourly_default_man_amhg_val <- bam_val(Sac_rch3_NED_hourly_default_man_amhg, Sac_rch3_qobs_hourly)

Sac_rch3_data_3m_hourly <- bam_data(w = t(Sac_rch3_width_3m_hourly), s = Sac_rch3_slope_hourly, dA = Sac_rch3_dA_3m_hourly, 
                                    Qhat = Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "manning_amhg")
Sac_rch3_3m_hourly_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "amhg")
Sac_rch3_3m_hourly_default_man <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "manning")

Sac_rch3_3m_hourly_default_amhg_val <- bam_val(Sac_rch3_3m_hourly_default_amhg, Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_default_man_val <- bam_val(Sac_rch3_3m_hourly_default_man, Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_default_man_amhg_val <- bam_val(Sac_rch3_3m_hourly_default_man_amhg, Sac_rch3_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly Default NED",
     ylim = c(120, 275))
lines(Sac_rch3_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Sac_rch3_NED_hourly_default_man_val[[1]], col = "blue")
lines(Sac_rch3_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch3_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch3_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch3_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch3_data_NED_hourly <- bam_data(w = t(Sac_rch3_width_NED_hourly), s = Sac_rch3_slope_hourly, dA = Sac_rch3_dA_NED_hourly,
                                     Qhat = Sac_rch3_qobs_hourly)
Sac_best_priors <- bam_priors(bamdata= Sac_rch3_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Sac_rch3_qobs_hourly)),
                              upperbound_logQ = log(1.25 * max(Sac_rch3_qobs_hourly)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch3_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch3_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "manning_amhg", 
                                                  bampriors = Sac_best_priors)
Sac_rch3_NED_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_NED_hourly_best_man <- bam_estimate(bamdata = Sac_rch3_data_NED_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_NED_hourly_best_amhg_val <- bam_val(Sac_rch3_NED_hourly_best_amhg, Sac_rch3_qobs_hourly)
Sac_rch3_NED_hourly_best_man_val <- bam_val(Sac_rch3_NED_hourly_best_man, Sac_rch3_qobs_hourly)
Sac_rch3_NED_hourly_best_man_amhg_val <- bam_val(Sac_rch3_NED_hourly_best_man_amhg, Sac_rch3_qobs_hourly)

Sac_rch3_data_3m_hourly <- bam_data(w = t(Sac_rch3_width_3m_hourly), s = Sac_rch3_slope_hourly, dA = Sac_rch3_dA_3m_hourly,
                                    Qhat = Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "manning_amhg",
                                                 bampriors = Sac_best_priors)
Sac_rch3_3m_hourly_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_3m_hourly_best_man <- bam_estimate(bamdata = Sac_rch3_data_3m_hourly,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_3m_hourly_best_amhg_val <- bam_val(Sac_rch3_3m_hourly_best_amhg, Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_best_man_val <- bam_val(Sac_rch3_3m_hourly_best_man, Sac_rch3_qobs_hourly)
Sac_rch3_3m_hourly_best_man_amhg_val <- bam_val(Sac_rch3_3m_hourly_best_man_amhg, Sac_rch3_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento hourly best NED",
     ylim = c(135, 155))
lines(Sac_rch3_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch3_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch3_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch3_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch3_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch3_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#15min-----
#Default----
Sac_rch3_data_NED_15min <- bam_data(w = t(Sac_rch3_width_NED_15min), s = Sac_rch3_slope_15min, dA = Sac_rch3_dA_NED_15min,
                                    Qhat = Sac_rch3_qobs_15min)
Sac_rch3_NED_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "manning_amhg")
Sac_rch3_NED_15min_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "amhg")
Sac_rch3_NED_15min_default_man <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "manning")

library(hydroGOF)
Sac_rch3_NED_15min_default_amhg_val <- bam_val(Sac_rch3_NED_15min_default_amhg, Sac_rch3_qobs_15min)
Sac_rch3_NED_15min_default_man_val <- bam_val(Sac_rch3_NED_15min_default_man, Sac_rch3_qobs_15min)
Sac_rch3_NED_15min_default_man_amhg_val <- bam_val(Sac_rch3_NED_15min_default_man_amhg, Sac_rch3_qobs_15min)

Sac_rch3_data_3m_15min <- bam_data(w = t(Sac_rch3_width_3m_15min), s = Sac_rch3_slope_15min, dA = Sac_rch3_dA_3m_15min, 
                                   Qhat = Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_default_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "manning_amhg")
Sac_rch3_3m_15min_default_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "amhg")
Sac_rch3_3m_15min_default_man <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "manning")

library(hydroGOF)
Sac_rch3_3m_15min_default_amhg_val <- bam_val(Sac_rch3_3m_15min_default_amhg, Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_default_man_val <- bam_val(Sac_rch3_3m_15min_default_man, Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_default_man_amhg_val <- bam_val(Sac_rch3_3m_15min_default_man_amhg, Sac_rch3_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min Default NED",
     ylim = c(120, 275))
lines(Sac_rch3_NED_15min_default_amhg_val[[1]], col = "red")
lines(Sac_rch3_NED_15min_default_man_val[[1]], col = "blue")
lines(Sac_rch3_NED_15min_default_man_amhg_val[[1]], col = "purple")
lines(Sac_rch3_3m_15min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Sac_rch3_3m_15min_default_man_val[[1]], col = "blue", lty = 2)
lines(Sac_rch3_3m_15min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Sac_rch3_data_NED_15min <- bam_data(w = t(Sac_rch3_width_NED_15min), s = Sac_rch3_slope_15min, dA = Sac_rch3_dA_NED_15min,
                                    Qhat = Sac_rch3_qobs_15min)
Sac_best_priors <- bam_priors(bamdata= Sac_rch3_data_NED_15min, lowerbound_logQ = log(0.75 * min(Sac_rch3_qobs_15min)),
                              upperbound_logQ = log(1.25 * max(Sac_rch3_qobs_15min)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Sac_rch3_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Sac_rch3_NED_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "manning_amhg", 
                                                 bampriors = Sac_best_priors)
Sac_rch3_NED_15min_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_NED_15min_best_man <- bam_estimate(bamdata = Sac_rch3_data_NED_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_NED_15min_best_amhg_val <- bam_val(Sac_rch3_NED_15min_best_amhg, Sac_rch3_qobs_15min)
Sac_rch3_NED_15min_best_man_val <- bam_val(Sac_rch3_NED_15min_best_man, Sac_rch3_qobs_15min)
Sac_rch3_NED_15min_best_man_amhg_val <- bam_val(Sac_rch3_NED_15min_best_man_amhg, Sac_rch3_qobs_15min)

Sac_rch3_data_3m_15min <- bam_data(w = t(Sac_rch3_width_3m_15min), s = Sac_rch3_slope_15min, dA = Sac_rch3_dA_3m_15min,
                                   Qhat = Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_best_man_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "manning_amhg",
                                                bampriors = Sac_best_priors)
Sac_rch3_3m_15min_best_amhg <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "amhg", bampriors = Sac_best_priors)
Sac_rch3_3m_15min_best_man <- bam_estimate(bamdata = Sac_rch3_data_3m_15min,variant = "manning", bampriors = Sac_best_priors)

Sac_rch3_3m_15min_best_amhg_val <- bam_val(Sac_rch3_3m_15min_best_amhg, Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_best_man_val <- bam_val(Sac_rch3_3m_15min_best_man, Sac_rch3_qobs_15min)
Sac_rch3_3m_15min_best_man_amhg_val <- bam_val(Sac_rch3_3m_15min_best_man_amhg, Sac_rch3_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Sac_rch3_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Sacramento 15min best NED",
     ylim = c(135, 155))
lines(Sac_rch3_NED_15min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Sac_rch3_NED_15min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Sac_rch3_NED_15min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Sac_rch3_3m_15min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Sac_rch3_3m_15min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Sac_rch3_3m_15min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Sac_rch3_daily_lower_logQ <- log(seq(exp(max(apply(log(Sac_rch3_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                     min(Sac_rch3_qobs_15min),length.out = 10))
Sac_rch3_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_logQ = Sac_rch3_daily_lower_logQ[i])
}
#upperbound_logQ
Sac_rch3_daily_upper_logQ <- log(seq(max(Sac_rch3_qobs_15min), exp(min(apply(log(Sac_rch3_width_NED_daily), 2, max)) + log(40) + log(5)),
                                     length.out = 10))
Sac_rch3_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_logQ = Sac_rch3_daily_upper_logQ[i])
}
#lowerbound_A0
Sac_rch3_daily_lower_A0 <- seq(0, 30, length.out = 10)
Sac_rch3_daily_lower_A0_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_A0 = Sac_rch3_daily_lower_A0[i])
}
#upperbound_A0
Sac_rch3_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Sac_rch3_daily_upper_A0_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_A0 = Sac_rch3_daily_upper_A0[i])
}
#lowerbound_logn
Sac_rch3_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Sac_rch3_daily_lower_logn_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_logn = Sac_rch3_daily_lower_logn[i])
}
#upperbound_logn
Sac_rch3_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Sac_rch3_daily_upper_logn_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_logn = Sac_rch3_daily_upper_logn[i])
}
#lowerbound_logQc
Sac_rch3_daily_lower_logQc <- seq(0, log(min(Sac_rch3_qobs_15min)),length.out = 10)
Sac_rch3_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_logQc = Sac_rch3_daily_lower_logQc[i])
}
#upperbound_logQc
Sac_rch3_daily_upper_logQc <- seq(log(max(Sac_rch3_qobs_15min)), 10, length.out = 10)
Sac_rch3_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_logQc = Sac_rch3_daily_upper_logQc[i])
}
#lowerbound_logWc
Sac_rch3_daily_lower_logWc <- seq(1, log(min(Sac_rch3_width_NED_daily)),length.out = 10)
Sac_rch3_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_logWc = Sac_rch3_daily_lower_logWc[i])
}
#upperbound_logWc
Sac_rch3_daily_upper_logWc <- seq(log(max(Sac_rch3_width_NED_daily)), 8, length.out = 10)
Sac_rch3_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_logWc = Sac_rch3_daily_upper_logWc[i])
}
#lowerbound_b
Sac_rch3_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Sac_rch3_daily_lower_b_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_lower_b_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, lowerbound_b = Sac_rch3_daily_lower_b[i])
}
#upperbound_b
Sac_rch3_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Sac_rch3_daily_upper_b_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_upper_b_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, upperbound_b = Sac_rch3_daily_upper_b[i])
}
#sigma_man
Sac_rch3_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Sac_rch3_daily_sigma_man_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, sigma_man = Sac_rch3_daily_sigma_man[i])
}
#sigma_amhg
Sac_rch3_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Sac_rch3_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, sigma_amhg = Sac_rch3_daily_sigma_amhg[i])
}
#logQc_hat
Sac_rch3_daily_logQc_hat <- log(seq(min(Sac_rch3_qobs_15min), max(Sac_rch3_qobs_15min), length.out = 10))
Sac_rch3_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logQc_hat = Sac_rch3_daily_logQc_hat[i])
}
#logWc_hat
Sac_rch3_daily_logWc_hat <- log(seq(min(Sac_rch3_width_NED_daily), max(Sac_rch3_width_NED_daily), length.out = 10))
Sac_rch3_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logWc_hat = Sac_rch3_daily_logWc_hat[i])
}
#b_hat
Sac_rch3_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Sac_rch3_daily_b_hat_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_b_hat_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, b_hat = Sac_rch3_daily_b_hat[i])
}
#logA0_hat
Sac_rch3_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Sac_rch3_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logA0_hat = Sac_rch3_daily_logA0_hat[i])
}
#logn_hat
Sac_rch3_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Sac_rch3_daily_logn_hat_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logn_hat = Sac_rch3_daily_logn_hat[i])
}
#logQ_sd
Sac_rch3_daily_logQ_sd <- seq(sd(log(Sac_rch3_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch3_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logQ_sd = Sac_rch3_daily_logQ_sd[i])
}
#logQc_sd
Sac_rch3_daily_logQc_sd <- seq(sd(log(Sac_rch3_qobs_15min)), 0.8325546, length.out = 10)
Sac_rch3_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logQc_sd = Sac_rch3_daily_logQc_sd[i])
}
#logWc_Sd
Sac_rch3_daily_logWc_Sd <- seq(log(sd(Sac_rch3_width_NED_daily)), 4.712493, length.out = 10)
Sac_rch3_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logWc_sd = Sac_rch3_daily_logWc_Sd[i])
}
#b_sd
Sac_rch3_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Sac_rch3_daily_b_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_b_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, b_sd = Sac_rch3_daily_b_sd[i])
}
#logA0_sd
Sac_rch3_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Sac_rch3_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logA0_sd = Sac_rch3_daily_logA0_sd[i])
}
#logn_sd
Sac_rch3_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Sac_rch3_daily_logn_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, logn_sd = Sac_rch3_daily_logn_sd[i])
}
#Werr_sd
Sac_rch3_daily_Werr_sd <- seq(5, 15, length.out = 10)
Sac_rch3_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, Werr_sd = Sac_rch3_daily_Werr_sd[i])
}
#Serr_sd
Sac_rch3_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Sac_rch3_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, Serr_sd = Sac_rch3_daily_Serr_sd[i])
}
#dAerr_sd
Sac_rch3_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Sac_rch3_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Sac_rch3_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Sac_rch3_data_NED_daily, dAerr_sd = Sac_rch3_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Sac_rch3_daily_manning_list <- list(Sac_rch3_daily_lower_logQ_priors,
                                    Sac_rch3_daily_upper_logQ_priors,
                                    Sac_rch3_daily_lower_A0_priors,
                                    Sac_rch3_daily_upper_A0_priors,
                                    Sac_rch3_daily_lower_logn_priors,
                                    Sac_rch3_daily_upper_logn_priors,
                                    Sac_rch3_daily_sigma_man_priors,
                                    Sac_rch3_daily_logA0_hat_priors,
                                    Sac_rch3_daily_logn_hat_priors,
                                    Sac_rch3_daily_logQ_sd_priors,
                                    Sac_rch3_daily_logA0_sd_priors,
                                    Sac_rch3_daily_logn_sd_priors,
                                    Sac_rch3_daily_Werr_sd_priors,
                                    Sac_rch3_daily_Serr_sd_priors,
                                    Sac_rch3_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Sac_rch3_data_NED_daily, variant = "manning", bampriors = Sac_rch3_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Sac_rch3_daily_amhg_list <- list(Sac_rch3_daily_lower_logQ_priors,
                                 Sac_rch3_daily_upper_logQ_priors,
                                 Sac_rch3_daily_lower_logQc_priors,
                                 Sac_rch3_daily_upper_logQc_priors,
                                 Sac_rch3_daily_lower_logWc_priors,
                                 Sac_rch3_daily_upper_logWc_priors,
                                 Sac_rch3_daily_lower_b_priors,
                                 Sac_rch3_daily_upper_b_priors,
                                 Sac_rch3_daily_sigma_amhg_priors,
                                 Sac_rch3_daily_logQc_hat_priors,
                                 Sac_rch3_daily_logWc_hat_priors,
                                 Sac_rch3_daily_b_hat_priors,
                                 Sac_rch3_daily_logQ_sd_priors,
                                 Sac_rch3_daily_logQc_sd_priors,
                                 Sac_rch3_daily_logWc_Sd_priors,
                                 Sac_rch3_daily_b_sd_priors,
                                 Sac_rch3_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Sac_rch3_data_NED_daily, variant = "amhg", bampriors = Sac_rch3_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Sac_rch3_daily_man_amhg_list <- list(Sac_rch3_daily_lower_logQ_priors,
                                     Sac_rch3_daily_upper_logQ_priors,
                                     Sac_rch3_daily_lower_A0_priors,
                                     Sac_rch3_daily_upper_A0_priors,
                                     Sac_rch3_daily_lower_logn_priors,
                                     Sac_rch3_daily_upper_logn_priors,
                                     Sac_rch3_daily_lower_logQc_priors,
                                     Sac_rch3_daily_upper_logQc_priors,
                                     Sac_rch3_daily_lower_logWc_priors,
                                     Sac_rch3_daily_upper_logWc_priors,
                                     Sac_rch3_daily_lower_b_priors,
                                     Sac_rch3_daily_upper_b_priors,
                                     Sac_rch3_daily_sigma_man_priors,
                                     Sac_rch3_daily_sigma_amhg_priors,
                                     Sac_rch3_daily_logQc_hat_priors,
                                     Sac_rch3_daily_logWc_hat_priors,
                                     Sac_rch3_daily_b_hat_priors,
                                     Sac_rch3_daily_logA0_hat_priors,
                                     Sac_rch3_daily_logn_hat_priors,
                                     Sac_rch3_daily_logQ_sd_priors,
                                     Sac_rch3_daily_logQc_sd_priors,
                                     Sac_rch3_daily_logWc_Sd_priors,
                                     Sac_rch3_daily_b_sd_priors,
                                     Sac_rch3_daily_logA0_sd_priors,
                                     Sac_rch3_daily_logn_sd_priors,
                                     Sac_rch3_daily_Werr_sd_priors,
                                     Sac_rch3_daily_Serr_sd_priors,
                                     Sac_rch3_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Sacramento/Reach_3/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Sac_rch3_data_NED_daily, variant = "manning_amhg", bampriors = Sac_rch3_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

###Willamette_Reach_Three
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Willamette"))
load("Will_WSE_daily.RData")
load("Will_WSE_6hr.RData")
load("Will_WSE_hourly.RData")
load("Will_WSE_15min.RData")

Will_PT_coords <- read.csv("Willamette_PT_coords.csv")
Will_PT_coords_sp <- SpatialPoints(Will_PT_coords[,c(3,4)], proj4string =  CRS("+proj=utm +zone=10 +datum=WGS84"))
Will_cl <-readOGR(dsn = ".", layer = "Will_cl_UTM")

setwd(paste0(home_dir, "/Willamette/DEM"))
Will_DEM_NED <- raster("Cropped_Will_NED.tif")
Will_DEM_3m <- raster()

setwd(paste0(home_dir, "/Willamette"))
Will_orthos_NED <- ortho_lines(Will_DEM_NED, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords), 500)
Will_dist_mat_NED <- pointDistance(Will_PT_coords_sp, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords))
Will_closest_ortho_NED <- Will_orthos_NED[apply(Will_dist_mat_NED, 1, FUN = which.min)]

Will_orthos_3m <- ortho_lines(Will_DEM_3m, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords), 500)
Will_dist_mat_3m <- pointDistance(Will_PT_coords_sp, data.frame(Will_cl@lines[[1]]@Lines[[1]]@coords))
Will_closest_ortho_3m <- Will_orthos_3m[apply(Will_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Will_elev_prof_NED <- foreach(i=1:length(Will_closest_ortho_NED), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Will_DEM_NED, Will_closest_ortho_NED[[i]], method = 'simple')
                              }
Will_elev_prof_3m <- foreach(i=1:length(Will_closest_ortho_3m), 
                             .combine='rbind', .packages=c('raster')) %dopar% {
                               extract(Will_DEM_3m, Will_closest_ortho_3m[[i]], method = 'simple')
                             }
proc.time()-ptm
stopCluster(cl)

Will_scale_x_NED <- vector()
for(i in 1:length(Will_elev_prof_NED)){
  Will_scale_x_NED[i] = 1000/length(Will_elev_prof_NED[[i]])
}

Will_scale_x_3m <- vector()
for(i in 1:length(Will_elev_prof_3m)){
  Will_scale_x_3m[i] = 1000/length(Will_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Will_elev_prof_NED, Will_WSE_daily_df, "Willamette Reach Three")
locator(n = 4)
Will_upper_l_NED <- c(57.1, 56.2, 71.9,  66.9, 70.3, 65.7,  69.2,  73.6,  68.0, 43.7, 62.8, 61.9, 35.0, 55.1, 64.2, 69.8, 73.6, 60.2, 71.8, 54.4, 90.5)
Will_lower_l_NED <- c(63.2, 60.9, 77.8,  74.1, 74.7, 70.1,  75.4,  78.2,  72.6, 46.9, 67.9, 65.5, 38.3, 59.5, 69.3, 75.5, 79.3, 64.0, 75.7, 58.6, 95.8)
Will_lower_r_NED <- c(75.9, 68.9, 88.6,  96.5, 92.4, 78.8, 111.1, 107.3,  96.9, 58.3, 80.1, 85.3, 59.3, 78.0, 93.4, 92.9, 89.1, 72.0, 87.2, 63.7, 114.9)
Will_upper_r_NED <- c(80.2, 73.0, 96.3, 101.2, 97.4, 86.1, 116.6, 113.8, 102.3, 60.3, 86.4, 88.4, 62.4, 82.8, 97.2, 97.2, 94.5, 76.8, 90.9, 68.6, 120.0)

Will_elev_prof_NED_corr <- Will_elev_prof_NED[c(1:20)]
#Will_elev_prof_NED_corr[[12]] <- Will_elev_prof_NED[[13]]
#Will_scale_x_NED[12] <- Will_scale_x_NED[13]

Will_width_NED_daily <- calc_width(h = t(Will_WSE_daily_df), elev_prof = Will_elev_prof_NED_corr, 
                                   upper_l = Will_upper_l_NED, lower_l = Will_lower_l_NED,
                                   lower_r = Will_lower_r_NED, upper_r = Will_upper_r_NED,
                                   scale_prof = Will_scale_x_NED)
Will_width_NED_daily <- na.approx(Will_width_NED_daily)
Will_dA_NED_daily <- calcdA_mat(w = t(Will_width_NED_daily), h = Will_WSE_daily_df)

Will_width_NED_6hr <- calc_width(h = t(Will_WSE_6hr_df), elev_prof = Will_elev_prof_NED_corr, 
                                 upper_l = Will_upper_l_NED, lower_l = Will_lower_l_NED,
                                 lower_r = Will_lower_r_NED, upper_r = Will_upper_r_NED,
                                 scale_prof = Will_scale_x_NED)
Will_width_NED_6hr <- na.approx(Will_width_NED_6hr)
Will_dA_NED_6hr <- calcdA_mat(w = t(Will_width_NED_6hr), h = Will_WSE_6hr_df)

Will_width_NED_hourly <- calc_width(h = t(Will_WSE_hourly_df), elev_prof = Will_elev_prof_NED_corr, 
                                    upper_l = Will_upper_l_NED, lower_l = Will_lower_l_NED,
                                    lower_r = Will_lower_r_NED, upper_r = Will_upper_r_NED,
                                    scale_prof = Will_scale_x_NED)
Will_width_NED_hourly <- na.approx(Will_width_NED_hourly)
Will_dA_NED_hourly <- calcdA_mat(w = t(Will_width_NED_hourly), h = Will_WSE_hourly_df)

Will_width_NED_15min <- calc_width(h = t(Will_WSE_df), elev_prof = Will_elev_prof_NED_corr, 
                                   upper_l = Will_upper_l_NED, lower_l = Will_lower_l_NED,
                                   lower_r = Will_lower_r_NED, upper_r = Will_upper_r_NED,
                                   scale_prof = Will_scale_x_NED)
Will_width_NED_15min <- na.approx(Will_width_NED_15min)
Will_dA_NED_15min <- calcdA_mat(w = t(Will_width_NED_15min), h = Will_WSE_df)

plot_elev_prof(Will_elev_prof_3m, Will_WSE_daily_df, "Willamette Reach Three")
locator(n = 4)
Will_upper_l_3m <- c(157.7, 169.0, 117.5, 201.2, 143.3, 129.0, 143.9,  33.0, 132.9, 230.4, 182.2, 211.3,  88.1, 155.0,  63.5,  90.4, 197.4, 203.8, 106.4)
Will_lower_l_3m <- c(175.8, 176.1, 123.6, 216.4, 156.9, 136.9, 150.7,  54.5, 159.9, 237.0, 191.0, 226.7,  94.7, 168.1,  72.7, 103.3, 204.3, 209.4, 130.2)
Will_lower_r_3m <- c(209.0, 218.0, 216.5, 303.7, 353.8, 343.5, 209.0, 455.0, 390.5, 376.9, 262.4, 309.3, 178.2, 275.2, 179.1, 160.2, 253.8, 242.9, 262.4)
Will_upper_r_3m <- c(231.2, 226.1, 225.0, 321.9, 371.9, 356.9, 225.0, 468.1, 401.4, 388.1, 283.4, 323.1, 184.8, 285.0, 207.6, 175.4, 266.2, 252.2, 275.2)

Will_width_3m_daily <- calc_width(h = t(Will_WSE_daily_df), elev_prof = Will_elev_prof_3m, 
                                  upper_l = Will_upper_l_3m, lower_l = Will_lower_l_3m,
                                  lower_r = Will_lower_r_3m, upper_r = Will_upper_r_3m,
                                  scale_prof = Will_scale_x_3m)
Will_width_3m_daily <- na.approx(Will_width_3m_daily)
Will_dA_3m_daily <- calcdA_mat(w = t(Will_width_3m_daily), h = Will_WSE_daily_df)

Will_width_3m_6hr <- calc_width(h = t(Will_WSE_6hr_df), elev_prof = Will_elev_prof_3m, 
                                upper_l = Will_upper_l_3m, lower_l = Will_lower_l_3m,
                                lower_r = Will_lower_r_3m, upper_r = Will_upper_r_3m,
                                scale_prof = Will_scale_x_3m)
Will_width_3m_6hr <- na.approx(Will_width_3m_6hr)
Will_dA_3m_6hr <- calcdA_mat(w = t(Will_width_3m_6hr), h = Will_WSE_6hr_df)

Will_width_3m_hourly <- calc_width(h = t(Will_WSE_hourly_df), elev_prof = Will_elev_prof_3m, 
                                   upper_l = Will_upper_l_3m, lower_l = Will_lower_l_3m,
                                   lower_r = Will_lower_r_3m, upper_r = Will_upper_r_3m,
                                   scale_prof = Will_scale_x_3m)
Will_width_3m_hourly <- na.approx(Will_width_3m_hourly)
Will_dA_3m_hourly <- calcdA_mat(w = t(Will_width_3m_hourly), h = Will_WSE_hourly_df)

Will_width_3m_15min <- calc_width(h = t(Will_WSE_df), elev_prof = Will_elev_prof_3m, 
                                  upper_l = Will_upper_l_3m, lower_l = Will_lower_l_3m,
                                  lower_r = Will_lower_r_3m, upper_r = Will_upper_r_3m,
                                  scale_prof = Will_scale_x_3m)
Will_width_3m_15min <- na.approx(Will_width_3m_15min)
Will_dA_3m_15min <- calcdA_mat(w = t(Will_width_3m_15min), h = Will_WSE_df)


xvec <- calc_xvec(Will_PT_coords_sp[c(1:20)], Will_cl)
Will_xvec <- saveRDS(xvec, file = "Will_xvec.rds")
Will_slope_daily <- calcslope((xvec), hmat = Will_WSE_daily_df)
Will_slope_daily <- na.replace(Will_slope_daily, mean(Will_slope_daily))
Will_slope_daily[1,] <- colMeans(Will_slope_daily, na.rm = TRUE)
Will_slope_6hr <- calcslope((xvec), hmat = Will_WSE_6hr_df)
Will_slope_6hr <- na.approx(Will_slope_6hr)
Will_slope_6hr[1,] <- colMeans(Will_slope_6hr, na.rm = TRUE)
Will_slope_hourly <- calcslope((xvec), hmat = Will_WSE_hourly_df)
Will_slope_hourly <- na.approx(Will_slope_hourly)
Will_slope_hourly[1,] <- colMeans(Will_slope_hourly, na.rm = TRUE)
Will_slope_15min <- calcslope((xvec), hmat = Will_WSE_df)
Will_slope_15min <- na.approx(Will_slope_15min)
Will_slope_15min[1,] <- colMeans(Will_slope_15min, na.rm = TRUE)

Willqobs <- read.csv("WillQobs.csv")$Q*0.028316847
Will_qobs_daily<- colMeans(matrix(Willqobs[1:4896], 96))
Will_qobs_6hr<- colMeans(matrix(Willqobs[1:4968], 24))
Will_qobs_hourly<- colMeans(matrix(Willqobs[1:4988], 4))
Will_qobs_15min<- colMeans(matrix(Willqobs[1:4988], 1))


save(Will_width_NED_daily, file = "Willamette_width_NED_daily.RData")
save(Will_dA_NED_daily, file = "Willamette_dA_NED_daily.RData")
save(Will_slope_daily, file = "Willamette_slope_daily.RData")
save(Will_qobs_daily, file = "Willamette_qobs_daily.RData")
save(Will_width_NED_6hr, file = "Willamette_width_NED_6hr.RData")
save(Will_dA_NED_6hr, file = "Willamette_dA_NED_6hr.RData")
save(Will_slope_6hr, file = "Willamette_slope_6hr.RData")
save(Will_qobs_6hr, file = "Willamette_qobs_6hr.RData")
save(Will_width_NED_hourly, file = "Willamette_width_NED_hourly.RData")
save(Will_dA_NED_hourly, file = "Willamette_dA_NED_hourly.RData")
save(Will_slope_hourly, file = "Willamette_slope_hourly.RData")
save(Will_qobs_hourly, file = "Willamette_qobs_hourly.RData")
save(Will_width_NED_15min, file = "Willamette_width_NED_15min.RData")
save(Will_dA_NED_15min, file = "Willamette_dA_NED_15min.RData")
save(Will_slope_15min, file = "Willamette_slope_15min.RData")
save(Will_qobs_15min, file = "Willamette_qobs_15min.RData")
#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Will_wse_plot <- matplot(t(Will_WSE_daily_df), type = c("l"), col= colfunc,
                         main = "WSE", xlab = "Days",
                         ylab = "WSE (m)",
                         cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Will_slope_plot <- matplot(t(Will_slope_daily), type = c("l"), col= colfunc,
                           main = "Slope", xlab = "Days",
                           ylab = "Slope",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Will_width_plot <- matplot(Will_width_NED_daily, type = c("l"), col= colfunc,
                           main = "Width (NED)", xlab = "Days",
                           ylab = "Width (m)",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Will_width_plot <- matplot(Will_width_3m_daily, type = c("l"), col= colfunc,
                           main = "Width (Lidar)", xlab = "Days",
                           ylab = "Width (m)",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Will_dA_plot <- matplot(t(Will_dA_NED_daily), type = c("l"), col= colfunc,
                        main = "dA (NED)", xlab = "Days",
                        ylab = "dA (m^2)",
                        cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Will_dA_plot <- matplot(t(Will_dA_3m_daily), type = c("l"), col= colfunc,
                        main = "dA (Lidar)", xlab = "Days",
                        ylab = "dA (m^2)",
                        cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
Will_data_NED_daily <- bam_data(w = t(Will_width_NED_daily), s = Will_slope_daily, dA = Will_dA_NED_daily,
                                Qhat = Will_qobs_daily)
Will_NED_daily_default_man_amhg <- bam_estimate(bamdata = Will_data_NED_daily,variant = "manning_amhg")
Will_NED_daily_default_amhg <- bam_estimate(bamdata = Will_data_NED_daily,variant = "amhg")
Will_NED_daily_default_man <- bam_estimate(bamdata = Will_data_NED_daily,variant = "manning")
library(hydroGOF)
Will_NED_daily_default_amhg_val <- bam_val(Will_NED_daily_default_amhg, Will_qobs_daily)
Will_NED_daily_default_man_val <- bam_val(Will_NED_daily_default_man, Will_qobs_daily)
Will_NED_daily_default_man_amhg_val <- bam_val(Will_NED_daily_default_man_amhg, Will_qobs_daily)

Will_data_3m_daily <- bam_data(w = t(Will_width_3m_daily), s = Will_slope_daily, dA = Will_dA_3m_daily,
                               Qhat = Will_qobs_daily)
Will_3m_daily_default_man_amhg <- bam_estimate(bamdata = Will_data_3m_daily,variant = "manning_amhg")
Will_3m_daily_default_amhg <- bam_estimate(bamdata = Will_data_3m_daily,variant = "amhg")
Will_3m_daily_default_man <- bam_estimate(bamdata = Will_data_3m_daily,variant = "manning")

Will_3m_daily_default_amhg_val <- bam_val(Will_3m_daily_default_amhg, Will_qobs_daily)
Will_3m_daily_default_man_val <- bam_val(Will_3m_daily_default_man, Will_qobs_daily)
Will_3m_daily_default_man_amhg_val <- bam_val(Will_3m_daily_default_man_amhg, Will_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette Daily Default NED",
     ylim = c(140, 300))
lines(Will_NED_daily_default_amhg_val[[1]], col = "red")
lines(Will_NED_daily_default_man_val[[1]], col = "blue")
lines(Will_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Will_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Will_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Will_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Will_data_NED_daily <- bam_data(w = t(Will_width_NED_daily), s = Will_slope_daily, dA = Will_dA_NED_daily,
                                Qhat = Will_qobs_daily)
Will_best_priors <- bam_priors(bamdata= Will_data_NED_daily, lowerbound_logQ = log(0.35 * min(Will_qobs_daily)),
                               upperbound_logQ = log(1.75 * max(Will_qobs_daily)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Will_slope_daily, na.rm = TRUE),
                               lowerbound_A0 = 30, 
                               upperbound_A0 = 1e+05)
Will_NED_daily_best_man_amhg <- bam_estimate(bamdata = Will_data_NED_daily,variant = "manning_amhg", 
                                             bampriors = Will_best_priors)
Will_NED_daily_best_amhg <- bam_estimate(bamdata = Will_data_NED_daily,variant = "amhg", bampriors = Will_best_priors)
Will_NED_daily_best_man <- bam_estimate(bamdata = Will_data_NED_daily,variant = "manning", bampriors = Will_best_priors)
library(hydroGOF)
Will_NED_daily_best_amhg_val <- bam_val(Will_NED_daily_best_amhg, Will_qobs_daily)
Will_NED_daily_best_man_val <- bam_val(Will_NED_daily_best_man, Will_qobs_daily)
Will_NED_daily_best_man_amhg_val <- bam_val(Will_NED_daily_best_man_amhg, Will_qobs_daily)

Will_data_3m_daily <- bam_data(w = t(Will_width_3m_daily), s = Will_slope_daily, dA = Will_dA_3m_daily,
                               Qhat = Will_qobs_daily)
Will_3m_daily_best_man_amhg <- bam_estimate(bamdata = Will_data_3m_daily,variant = "manning_amhg",
                                            bampriors = Will_best_priors)
Will_3m_daily_best_amhg <- bam_estimate(bamdata = Will_data_3m_daily,variant = "amhg", bampriors = Will_best_priors)
Will_3m_daily_best_man <- bam_estimate(bamdata = Will_data_3m_daily,variant = "manning", bampriors = Will_best_priors)

Will_3m_daily_best_amhg_val <- bam_val(Will_3m_daily_best_amhg, Will_qobs_daily)
Will_3m_daily_best_man_val <- bam_val(Will_3m_daily_best_man, Will_qobs_daily)
Will_3m_daily_best_man_amhg_val <- bam_val(Will_3m_daily_best_man_amhg, Will_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette Daily best NED",
     ylim = c(140, 300))
lines(Will_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Will_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Will_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Will_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Will_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Will_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
Will_data_NED_6hr <- bam_data(w = t(Will_width_NED_6hr), s = Will_slope_6hr, dA = Will_dA_NED_6hr,
                              Qhat = Will_qobs_6hr)
Will_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "manning_amhg")
Will_NED_6hr_default_amhg <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "amhg")
Will_NED_6hr_default_man <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "manning")

Will_NED_6hr_default_amhg_val <- bam_val(Will_NED_6hr_default_amhg, Will_qobs_6hr)
Will_NED_6hr_default_man_val <- bam_val(Will_NED_6hr_default_man, Will_qobs_6hr)
Will_NED_6hr_default_man_amhg_val <- bam_val(Will_NED_6hr_default_man_amhg, Will_qobs_6hr)

Will_data_3m_6hr <- bam_data(w = t(Will_width_3m_6hr), s = Will_slope_6hr, dA = Will_dA_3m_6hr, 
                             Qhat = Will_qobs_6hr)
Will_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "manning_amhg")
Will_3m_6hr_default_amhg <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "amhg")
Will_3m_6hr_default_man <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "manning")

Will_3m_6hr_default_amhg_val <- bam_val(Will_3m_6hr_default_amhg, Will_qobs_6hr)
Will_3m_6hr_default_man_val <- bam_val(Will_3m_6hr_default_man, Will_qobs_6hr)
Will_3m_6hr_default_man_amhg_val <- bam_val(Will_3m_6hr_default_man_amhg, Will_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette 6hr Default NED",
     ylim = c(120, 275))
lines(Will_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Will_NED_6hr_default_man_val[[1]], col = "blue")
lines(Will_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Will_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Will_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Will_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Will_data_NED_6hr <- bam_data(w = t(Will_width_NED_6hr), s = Will_slope_6hr, dA = Will_dA_NED_6hr,
                              Qhat = Will_qobs_6hr)
Will_best_priors <- bam_priors(bamdata= Will_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Will_qobs_6hr)),
                               upperbound_logQ = log(1.25 * max(Will_qobs_6hr)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Will_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Will_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "manning_amhg", 
                                           bampriors = Will_best_priors)
Will_NED_6hr_best_amhg <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "amhg", bampriors = Will_best_priors)
Will_NED_6hr_best_man <- bam_estimate(bamdata = Will_data_NED_6hr,variant = "manning", bampriors = Will_best_priors)

Will_NED_6hr_best_amhg_val <- bam_val(Will_NED_6hr_best_amhg, Will_qobs_6hr)
Will_NED_6hr_best_man_val <- bam_val(Will_NED_6hr_best_man, Will_qobs_6hr)
Will_NED_6hr_best_man_amhg_val <- bam_val(Will_NED_6hr_best_man_amhg, Will_qobs_6hr)

Will_data_3m_6hr <- bam_data(w = t(Will_width_3m_6hr), s = Will_slope_6hr, dA = Will_dA_3m_6hr,
                             Qhat = Will_qobs_6hr)
Will_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "manning_amhg",
                                          bampriors = Will_best_priors)
Will_3m_6hr_best_amhg <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "amhg", bampriors = Will_best_priors)
Will_3m_6hr_best_man <- bam_estimate(bamdata = Will_data_3m_6hr,variant = "manning", bampriors = Will_best_priors)

Will_3m_6hr_best_amhg_val <- bam_val(Will_3m_6hr_best_amhg, Will_qobs_6hr)
Will_3m_6hr_best_man_val <- bam_val(Will_3m_6hr_best_man, Will_qobs_6hr)
Will_3m_6hr_best_man_amhg_val <- bam_val(Will_3m_6hr_best_man_amhg, Will_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette 6hr best NED",
     ylim = c(135, 155))
lines(Will_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Will_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Will_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Will_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Will_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Will_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Will_data_NED_hourly <- bam_data(w = t(Will_width_NED_hourly), s = Will_slope_hourly, dA = Will_dA_NED_hourly,
                                 Qhat = Will_qobs_hourly)
Will_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "manning_amhg")
Will_NED_hourly_default_amhg <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "amhg")
Will_NED_hourly_default_man <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "manning")

Will_NED_hourly_default_amhg_val <- bam_val(Will_NED_hourly_default_amhg, Will_qobs_hourly)
Will_NED_hourly_default_man_val <- bam_val(Will_NED_hourly_default_man, Will_qobs_hourly)
Will_NED_hourly_default_man_amhg_val <- bam_val(Will_NED_hourly_default_man_amhg, Will_qobs_hourly)

Will_data_3m_hourly <- bam_data(w = t(Will_width_3m_hourly), s = Will_slope_hourly, dA = Will_dA_3m_hourly, 
                                Qhat = Will_qobs_hourly)
Will_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "manning_amhg")
Will_3m_hourly_default_amhg <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "amhg")
Will_3m_hourly_default_man <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "manning")

Will_3m_hourly_default_amhg_val <- bam_val(Will_3m_hourly_default_amhg, Will_qobs_hourly)
Will_3m_hourly_default_man_val <- bam_val(Will_3m_hourly_default_man, Will_qobs_hourly)
Will_3m_hourly_default_man_amhg_val <- bam_val(Will_3m_hourly_default_man_amhg, Will_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette hourly Default NED",
     ylim = c(120, 275))
lines(Will_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Will_NED_hourly_default_man_val[[1]], col = "blue")
lines(Will_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Will_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Will_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Will_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Will_data_NED_hourly <- bam_data(w = t(Will_width_NED_hourly), s = Will_slope_hourly, dA = Will_dA_NED_hourly,
                                 Qhat = Will_qobs_hourly)
Will_best_priors <- bam_priors(bamdata= Will_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Will_qobs_hourly)),
                               upperbound_logQ = log(1.25 * max(Will_qobs_hourly)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Will_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Will_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "manning_amhg", 
                                              bampriors = Will_best_priors)
Will_NED_hourly_best_amhg <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "amhg", bampriors = Will_best_priors)
Will_NED_hourly_best_man <- bam_estimate(bamdata = Will_data_NED_hourly,variant = "manning", bampriors = Will_best_priors)

Will_NED_hourly_best_amhg_val <- bam_val(Will_NED_hourly_best_amhg, Will_qobs_hourly)
Will_NED_hourly_best_man_val <- bam_val(Will_NED_hourly_best_man, Will_qobs_hourly)
Will_NED_hourly_best_man_amhg_val <- bam_val(Will_NED_hourly_best_man_amhg, Will_qobs_hourly)

Will_data_3m_hourly <- bam_data(w = t(Will_width_3m_hourly), s = Will_slope_hourly, dA = Will_dA_3m_hourly,
                                Qhat = Will_qobs_hourly)
Will_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "manning_amhg",
                                             bampriors = Will_best_priors)
Will_3m_hourly_best_amhg <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "amhg", bampriors = Will_best_priors)
Will_3m_hourly_best_man <- bam_estimate(bamdata = Will_data_3m_hourly,variant = "manning", bampriors = Will_best_priors)

Will_3m_hourly_best_amhg_val <- bam_val(Will_3m_hourly_best_amhg, Will_qobs_hourly)
Will_3m_hourly_best_man_val <- bam_val(Will_3m_hourly_best_man, Will_qobs_hourly)
Will_3m_hourly_best_man_amhg_val <- bam_val(Will_3m_hourly_best_man_amhg, Will_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette hourly best NED",
     ylim = c(135, 155))
lines(Will_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Will_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Will_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Will_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Will_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Will_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#15min-----
#Default----
Will_data_NED_15min <- bam_data(w = t(Will_width_NED_15min), s = Will_slope_15min, dA = Will_dA_NED_15min,
                                Qhat = Will_qobs_15min)
Will_NED_15min_default_man_amhg <- bam_estimate(bamdata = Will_data_NED_15min,variant = "manning_amhg")
Will_NED_15min_default_amhg <- bam_estimate(bamdata = Will_data_NED_15min,variant = "amhg")
Will_NED_15min_default_man <- bam_estimate(bamdata = Will_data_NED_15min,variant = "manning")

library(hydroGOF)
Will_NED_15min_default_amhg_val <- bam_val(Will_NED_15min_default_amhg, Will_qobs_15min)
Will_NED_15min_default_man_val <- bam_val(Will_NED_15min_default_man, Will_qobs_15min)
Will_NED_15min_default_man_amhg_val <- bam_val(Will_NED_15min_default_man_amhg, Will_qobs_15min)

Will_data_3m_15min <- bam_data(w = t(Will_width_3m_15min), s = Will_slope_15min, dA = Will_dA_3m_15min, 
                               Qhat = Will_qobs_15min)
Will_3m_15min_default_man_amhg <- bam_estimate(bamdata = Will_data_3m_15min,variant = "manning_amhg")
Will_3m_15min_default_amhg <- bam_estimate(bamdata = Will_data_3m_15min,variant = "amhg")
Will_3m_15min_default_man <- bam_estimate(bamdata = Will_data_3m_15min,variant = "manning")

library(hydroGOF)
Will_3m_15min_default_amhg_val <- bam_val(Will_3m_15min_default_amhg, Will_qobs_15min)
Will_3m_15min_default_man_val <- bam_val(Will_3m_15min_default_man, Will_qobs_15min)
Will_3m_15min_default_man_amhg_val <- bam_val(Will_3m_15min_default_man_amhg, Will_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette 15min Default NED",
     ylim = c(120, 275))
lines(Will_NED_15min_default_amhg_val[[1]], col = "red")
lines(Will_NED_15min_default_man_val[[1]], col = "blue")
lines(Will_NED_15min_default_man_amhg_val[[1]], col = "purple")
lines(Will_3m_15min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Will_3m_15min_default_man_val[[1]], col = "blue", lty = 2)
lines(Will_3m_15min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Will_data_NED_15min <- bam_data(w = t(Will_width_NED_15min), s = Will_slope_15min, dA = Will_dA_NED_15min,
                                Qhat = Will_qobs_15min)
Will_best_priors <- bam_priors(bamdata= Will_data_NED_15min, lowerbound_logQ = log(0.75 * min(Will_qobs_15min)),
                               upperbound_logQ = log(1.25 * max(Will_qobs_15min)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Will_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Will_NED_15min_best_man_amhg <- bam_estimate(bamdata = Will_data_NED_15min,variant = "manning_amhg", 
                                             bampriors = Will_best_priors)
Will_NED_15min_best_amhg <- bam_estimate(bamdata = Will_data_NED_15min,variant = "amhg", bampriors = Will_best_priors)
Will_NED_15min_best_man <- bam_estimate(bamdata = Will_data_NED_15min,variant = "manning", bampriors = Will_best_priors)

Will_NED_15min_best_amhg_val <- bam_val(Will_NED_15min_best_amhg, Will_qobs_15min)
Will_NED_15min_best_man_val <- bam_val(Will_NED_15min_best_man, Will_qobs_15min)
Will_NED_15min_best_man_amhg_val <- bam_val(Will_NED_15min_best_man_amhg, Will_qobs_15min)

Will_data_3m_15min <- bam_data(w = t(Will_width_3m_15min), s = Will_slope_15min, dA = Will_dA_3m_15min,
                               Qhat = Will_qobs_15min)
Will_3m_15min_best_man_amhg <- bam_estimate(bamdata = Will_data_3m_15min,variant = "manning_amhg",
                                            bampriors = Will_best_priors)
Will_3m_15min_best_amhg <- bam_estimate(bamdata = Will_data_3m_15min,variant = "amhg", bampriors = Will_best_priors)
Will_3m_15min_best_man <- bam_estimate(bamdata = Will_data_3m_15min,variant = "manning", bampriors = Will_best_priors)

Will_3m_15min_best_amhg_val <- bam_val(Will_3m_15min_best_amhg, Will_qobs_15min)
Will_3m_15min_best_man_val <- bam_val(Will_3m_15min_best_man, Will_qobs_15min)
Will_3m_15min_best_man_amhg_val <- bam_val(Will_3m_15min_best_man_amhg, Will_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Will_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Willamette 15min best NED",
     ylim = c(135, 155))
lines(Will_NED_15min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Will_NED_15min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Will_NED_15min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Will_3m_15min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Will_3m_15min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Will_3m_15min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Will_daily_lower_logQ <- log(seq(exp(max(apply(log(Will_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                 min(Will_qobs_15min),length.out = 10))
Will_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Will_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_logQ = Will_daily_lower_logQ[i])
}
#upperbound_logQ
Will_daily_upper_logQ <- log(seq(max(Will_qobs_15min), exp(min(apply(log(Will_width_NED_daily), 2, max)) + log(40) + log(5)),
                                 length.out = 10))
Will_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Will_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_logQ = Will_daily_upper_logQ[i])
}
#lowerbound_A0
Will_daily_lower_A0 <- seq(0, 30, length.out = 10)
Will_daily_lower_A0_priors <- list()
for(i in 1:10){
  Will_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_A0 = Will_daily_lower_A0[i])
}
#upperbound_A0
Will_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Will_daily_upper_A0_priors <- list()
for(i in 1:10){
  Will_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_A0 = Will_daily_upper_A0[i])
}
#lowerbound_logn
Will_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Will_daily_lower_logn_priors <- list()
for(i in 1:10){
  Will_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_logn = Will_daily_lower_logn[i])
}
#upperbound_logn
Will_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Will_daily_upper_logn_priors <- list()
for(i in 1:10){
  Will_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_logn = Will_daily_upper_logn[i])
}
#lowerbound_logQc
Will_daily_lower_logQc <- seq(0, log(min(Will_qobs_15min)),length.out = 10)
Will_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Will_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_logQc = Will_daily_lower_logQc[i])
}
#upperbound_logQc
Will_daily_upper_logQc <- seq(log(max(Will_qobs_15min)), 10, length.out = 10)
Will_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Will_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_logQc = Will_daily_upper_logQc[i])
}
#lowerbound_logWc
Will_daily_lower_logWc <- seq(1, log(min(Will_width_NED_daily)),length.out = 10)
Will_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Will_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_logWc = Will_daily_lower_logWc[i])
}
#upperbound_logWc
Will_daily_upper_logWc <- seq(log(max(Will_width_NED_daily)), 8, length.out = 10)
Will_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Will_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_logWc = Will_daily_upper_logWc[i])
}
#lowerbound_b
Will_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Will_daily_lower_b_priors <- list()
for(i in 1:10){
  Will_daily_lower_b_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, lowerbound_b = Will_daily_lower_b[i])
}
#upperbound_b
Will_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Will_daily_upper_b_priors <- list()
for(i in 1:10){
  Will_daily_upper_b_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, upperbound_b = Will_daily_upper_b[i])
}
#sigma_man
Will_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Will_daily_sigma_man_priors <- list()
for(i in 1:10){
  Will_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, sigma_man = Will_daily_sigma_man[i])
}
#sigma_amhg
Will_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Will_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Will_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, sigma_amhg = Will_daily_sigma_amhg[i])
}
#logQc_hat
Will_daily_logQc_hat <- log(seq(min(Will_qobs_15min), max(Will_qobs_15min), length.out = 10))
Will_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Will_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logQc_hat = Will_daily_logQc_hat[i])
}
#logWc_hat
Will_daily_logWc_hat <- log(seq(min(Will_width_NED_daily), max(Will_width_NED_daily), length.out = 10))
Will_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Will_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logWc_hat = Will_daily_logWc_hat[i])
}
#b_hat
Will_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Will_daily_b_hat_priors <- list()
for(i in 1:10){
  Will_daily_b_hat_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, b_hat = Will_daily_b_hat[i])
}
#logA0_hat
Will_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Will_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Will_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logA0_hat = Will_daily_logA0_hat[i])
}
#logn_hat
Will_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Will_daily_logn_hat_priors <- list()
for(i in 1:10){
  Will_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logn_hat = Will_daily_logn_hat[i])
}
#logQ_sd
Will_daily_logQ_sd <- seq(sd(log(Will_qobs_15min)), 0.8325546, length.out = 10)
Will_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Will_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logQ_sd = Will_daily_logQ_sd[i])
}
#logQc_sd
Will_daily_logQc_sd <- seq(sd(log(Will_qobs_15min)), 0.8325546, length.out = 10)
Will_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Will_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logQc_sd = Will_daily_logQc_sd[i])
}
#logWc_Sd
Will_daily_logWc_Sd <- seq(log(sd(Will_width_NED_daily)), 4.712493, length.out = 10)
Will_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Will_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logWc_sd = Will_daily_logWc_Sd[i])
}
#b_sd
Will_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Will_daily_b_sd_priors <- list()
for(i in 1:10){
  Will_daily_b_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, b_sd = Will_daily_b_sd[i])
}
#logA0_sd
Will_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Will_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Will_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logA0_sd = Will_daily_logA0_sd[i])
}
#logn_sd
Will_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Will_daily_logn_sd_priors <- list()
for(i in 1:10){
  Will_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, logn_sd = Will_daily_logn_sd[i])
}
#Werr_sd
Will_daily_Werr_sd <- seq(5, 15, length.out = 10)
Will_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Will_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, Werr_sd = Will_daily_Werr_sd[i])
}
#Serr_sd
Will_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Will_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Will_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, Serr_sd = Will_daily_Serr_sd[i])
}
#dAerr_sd
Will_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Will_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Will_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Will_data_NED_daily, dAerr_sd = Will_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Will_daily_manning_list <- list(Will_daily_lower_logQ_priors,
                                Will_daily_upper_logQ_priors,
                                Will_daily_lower_A0_priors,
                                Will_daily_upper_A0_priors,
                                Will_daily_lower_logn_priors,
                                Will_daily_upper_logn_priors,
                                Will_daily_sigma_man_priors,
                                Will_daily_logA0_hat_priors,
                                Will_daily_logn_hat_priors,
                                Will_daily_logQ_sd_priors,
                                Will_daily_logA0_sd_priors,
                                Will_daily_logn_sd_priors,
                                Will_daily_Werr_sd_priors,
                                Will_daily_Serr_sd_priors,
                                Will_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/Willamette/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Willamette/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Will_data_NED_daily, variant = "manning", bampriors = Will_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Will_daily_amhg_list <- list(Will_daily_lower_logQ_priors,
                             Will_daily_upper_logQ_priors,
                             Will_daily_lower_logQc_priors,
                             Will_daily_upper_logQc_priors,
                             Will_daily_lower_logWc_priors,
                             Will_daily_upper_logWc_priors,
                             Will_daily_lower_b_priors,
                             Will_daily_upper_b_priors,
                             Will_daily_sigma_amhg_priors,
                             Will_daily_logQc_hat_priors,
                             Will_daily_logWc_hat_priors,
                             Will_daily_b_hat_priors,
                             Will_daily_logQ_sd_priors,
                             Will_daily_logQc_sd_priors,
                             Will_daily_logWc_Sd_priors,
                             Will_daily_b_sd_priors,
                             Will_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/Willamette/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/Willamette/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Will_data_NED_daily, variant = "amhg", bampriors = Will_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Will_daily_man_amhg_list <- list(Will_daily_lower_logQ_priors,
                                 Will_daily_upper_logQ_priors,
                                 Will_daily_lower_A0_priors,
                                 Will_daily_upper_A0_priors,
                                 Will_daily_lower_logn_priors,
                                 Will_daily_upper_logn_priors,
                                 Will_daily_lower_logQc_priors,
                                 Will_daily_upper_logQc_priors,
                                 Will_daily_lower_logWc_priors,
                                 Will_daily_upper_logWc_priors,
                                 Will_daily_lower_b_priors,
                                 Will_daily_upper_b_priors,
                                 Will_daily_sigma_man_priors,
                                 Will_daily_sigma_amhg_priors,
                                 Will_daily_logQc_hat_priors,
                                 Will_daily_logWc_hat_priors,
                                 Will_daily_b_hat_priors,
                                 Will_daily_logA0_hat_priors,
                                 Will_daily_logn_hat_priors,
                                 Will_daily_logQ_sd_priors,
                                 Will_daily_logQc_sd_priors,
                                 Will_daily_logWc_Sd_priors,
                                 Will_daily_b_sd_priors,
                                 Will_daily_logA0_sd_priors,
                                 Will_daily_logn_sd_priors,
                                 Will_daily_Werr_sd_priors,
                                 Will_daily_Serr_sd_priors,
                                 Will_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/Willamette/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Willamette/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Will_data_NED_daily, variant = "manning_amhg", bampriors = Will_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

###OlenOlengy
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Olentangy"))
load("Olen_WSE_daily.RData")
load("Olen_WSE_6hr.RData")
load("Olen_WSE_hourly.RData")
load("Olen_WSE_15min.RData")

Olen_PT_coords <- read.csv("Olentangy_PT_coord.csv")
Olen_PT_coords_sp <- SpatialPoints(Olen_PT_coords[,c(3,4)], proj4string =  CRS("+proj=utm +zone=17 +datum=WGS84"))
Olen_cl <-readOGR(dsn = ".", layer = "Olen_cl_UTM")

setwd(paste0(home_dir, "/Olentangy/DEM"))
Olen_DEM_NED <- raster("Cropped_Olen_NED.tif")
Olen_DEM_3m <- raster()

setwd(paste0(home_dir, "/Olentangy"))
Olen_orthos_NED <- ortho_lines(Olen_DEM_NED, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords), 150)
Olen_dist_mat_NED <- pointDistance(Olen_PT_coords_sp, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords))
Olen_closest_ortho_NED <- Olen_orthos_NED[apply(Olen_dist_mat_NED, 1, FUN = which.min)]

Olen_orthos_3m <- ortho_lines(Olen_DEM_3m, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords), 150)
Olen_dist_mat_3m <- pointDistance(Olen_PT_coords_sp, data.frame(Olen_cl@lines[[1]]@Lines[[1]]@coords))
Olen_closest_ortho_3m <- Olen_orthos_3m[apply(Olen_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Olen_elev_prof_NED <- foreach(i=1:length(Olen_closest_ortho_NED), 
                              .combine='rbind', .packages=c('raster')) %dopar% {
                                extract(Olen_DEM_NED, Olen_closest_ortho_NED[[i]], method = 'simple')
                              }
Olen_elev_prof_3m <- foreach(i=1:length(Olen_closest_ortho_3m), 
                             .combine='rbind', .packages=c('raster')) %dopar% {
                               extract(Olen_DEM_3m, Olen_closest_ortho_3m[[i]], method = 'simple')
                             }
proc.time()-ptm
stopCluster(cl)

Olen_scale_x_NED <- vector()
for(i in 1:length(Olen_elev_prof_NED)){
  Olen_scale_x_NED[i] = 300/length(Olen_elev_prof_NED[[i]])
}

Olen_scale_x_3m <- vector()
for(i in 1:length(Olen_elev_prof_3m)){
  Olen_scale_x_3m[i] = 300/length(Olen_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Olen_elev_prof_NED, Olen_WSE_daily_df, "Olentangy")
locator(n = 4)
Olen_upper_l_NED <- c(12.0, 15.8, 17.4, 26.2, 18.7,  9.7, 17.7, 18.5, 18.7, 12.4, 11.2, 17.6, 15.5, 15.0, 18.4,  8.1, 11.2, 14.2, 17.7, 17.9)
Olen_lower_l_NED <- c(16.6, 18.2, 21.5, 29.0, 21.1, 12.9, 21.3, 21.6, 21.5, 17.2, 16.2, 20.4, 18.6, 20.0, 21.5, 16.6, 16.5, 17.5, 20.5, 20.7)
Olen_lower_r_NED <- c(21.7, 26.7, 24.7, 31.4, 29.6, 18.7, 25.6, 27.1, 26.2, 26.3, 22.2, 27.5, 27.1, 22.5, 27.6, 21.4, 22.3, 22.2, 22.5, 25.5)
Olen_upper_r_NED <- c(25.5, 29.0, 28.8, 34.6, 32.3, 22.4, 28.9, 32.6, 30.8, 28.9, 26.9, 30.0, 30.0, 25.5, 29.8, 23.9, 27.3, 27.8, 26.2, 32.5)

Olen_elev_prof_NED_corr <- Olen_elev_prof_NED
#Olen_elev_prof_NED_corr[[12]] <- Olen_elev_prof_NED[[13]]
#Olen_scale_x_NED[12] <- Olen_scale_x_NED[13]

Olen_width_NED_daily <- calc_width(h = t(Olen_WSE_daily_df), elev_prof = Olen_elev_prof_NED_corr, 
                                   upper_l = Olen_upper_l_NED, lower_l = Olen_lower_l_NED,
                                   lower_r = Olen_lower_r_NED, upper_r = Olen_upper_r_NED,
                                   scale_prof = Olen_scale_x_NED)
Olen_width_NED_daily <- na.approx(Olen_width_NED_daily)
Olen_dA_NED_daily <- calcdA_mat(w = t(Olen_width_NED_daily), h = Olen_WSE_daily_df)

Olen_width_NED_6hr <- calc_width(h = t(Olen_WSE_6hr_df), elev_prof = Olen_elev_prof_NED_corr, 
                                 upper_l = Olen_upper_l_NED, lower_l = Olen_lower_l_NED,
                                 lower_r = Olen_lower_r_NED, upper_r = Olen_upper_r_NED,
                                 scale_prof = Olen_scale_x_NED)
Olen_width_NED_6hr <- na.approx(Olen_width_NED_6hr)
Olen_dA_NED_6hr <- calcdA_mat(w = t(Olen_width_NED_6hr), h = Olen_WSE_6hr_df)

Olen_width_NED_hourly <- calc_width(h = t(Olen_WSE_hourly_df), elev_prof = Olen_elev_prof_NED_corr, 
                                    upper_l = Olen_upper_l_NED, lower_l = Olen_lower_l_NED,
                                    lower_r = Olen_lower_r_NED, upper_r = Olen_upper_r_NED,
                                    scale_prof = Olen_scale_x_NED)
Olen_width_NED_hourly <- na.approx(Olen_width_NED_hourly)
Olen_dA_NED_hourly <- calcdA_mat(w = t(Olen_width_NED_hourly), h = Olen_WSE_hourly_df)

Olen_width_NED_15min <- calc_width(h = t(Olen_WSE_df), elev_prof = Olen_elev_prof_NED_corr, 
                                   upper_l = Olen_upper_l_NED, lower_l = Olen_lower_l_NED,
                                   lower_r = Olen_lower_r_NED, upper_r = Olen_upper_r_NED,
                                   scale_prof = Olen_scale_x_NED)
Olen_width_NED_15min <- na.approx(Olen_width_NED_15min)
Olen_dA_NED_15min <- calcdA_mat(w = t(Olen_width_NED_15min), h = Olen_WSE_df)

plot_elev_prof(Olen_elev_prof_3m, Olen_WSE_daily_df, "Olentangy Reach Three")
locator(n = 4)
Olen_upper_l_3m <- c(157.7, 169.0, 117.5, 201.2, 143.3, 129.0, 143.9,  33.0, 132.9, 230.4, 182.2, 211.3,  88.1, 155.0,  63.5,  90.4, 197.4, 203.8, 106.4)
Olen_lower_l_3m <- c(175.8, 176.1, 123.6, 216.4, 156.9, 136.9, 150.7,  54.5, 159.9, 237.0, 191.0, 226.7,  94.7, 168.1,  72.7, 103.3, 204.3, 209.4, 130.2)
Olen_lower_r_3m <- c(209.0, 218.0, 216.5, 303.7, 353.8, 343.5, 209.0, 455.0, 390.5, 376.9, 262.4, 309.3, 178.2, 275.2, 179.1, 160.2, 253.8, 242.9, 262.4)
Olen_upper_r_3m <- c(231.2, 226.1, 225.0, 321.9, 371.9, 356.9, 225.0, 468.1, 401.4, 388.1, 283.4, 323.1, 184.8, 285.0, 207.6, 175.4, 266.2, 252.2, 275.2)

Olen_width_3m_daily <- calc_width(h = t(Olen_WSE_daily_df), elev_prof = Olen_elev_prof_3m, 
                                  upper_l = Olen_upper_l_3m, lower_l = Olen_lower_l_3m,
                                  lower_r = Olen_lower_r_3m, upper_r = Olen_upper_r_3m,
                                  scale_prof = Olen_scale_x_3m)
Olen_width_3m_daily <- na.approx(Olen_width_3m_daily)
Olen_dA_3m_daily <- calcdA_mat(w = t(Olen_width_3m_daily), h = Olen_WSE_daily_df)

Olen_width_3m_6hr <- calc_width(h = t(Olen_WSE_6hr_df), elev_prof = Olen_elev_prof_3m, 
                                upper_l = Olen_upper_l_3m, lower_l = Olen_lower_l_3m,
                                lower_r = Olen_lower_r_3m, upper_r = Olen_upper_r_3m,
                                scale_prof = Olen_scale_x_3m)
Olen_width_3m_6hr <- na.approx(Olen_width_3m_6hr)
Olen_dA_3m_6hr <- calcdA_mat(w = t(Olen_width_3m_6hr), h = Olen_WSE_6hr_df)

Olen_width_3m_hourly <- calc_width(h = t(Olen_WSE_hourly_df), elev_prof = Olen_elev_prof_3m, 
                                   upper_l = Olen_upper_l_3m, lower_l = Olen_lower_l_3m,
                                   lower_r = Olen_lower_r_3m, upper_r = Olen_upper_r_3m,
                                   scale_prof = Olen_scale_x_3m)
Olen_width_3m_hourly <- na.approx(Olen_width_3m_hourly)
Olen_dA_3m_hourly <- calcdA_mat(w = t(Olen_width_3m_hourly), h = Olen_WSE_hourly_df)

Olen_width_3m_15min <- calc_width(h = t(Olen_WSE_df), elev_prof = Olen_elev_prof_3m, 
                                  upper_l = Olen_upper_l_3m, lower_l = Olen_lower_l_3m,
                                  lower_r = Olen_lower_r_3m, upper_r = Olen_upper_r_3m,
                                  scale_prof = Olen_scale_x_3m)
Olen_width_3m_15min <- na.approx(Olen_width_3m_15min)
Olen_dA_3m_15min <- calcdA_mat(w = t(Olen_width_3m_15min), h = Olen_WSE_df)


xvec <- calc_xvec(Olen_PT_coords_sp[c(1:20)], Olen_cl)
Olen_xvec <- saveRDS(xvec, file = "Olen_xvec.rds")
Olen_slope_daily <- calcslope(rev(xvec), hmat = Olen_WSE_daily_df)
Olen_slope_daily <- na.replace(Olen_slope_daily, mean(Olen_slope_daily))
Olen_slope_daily[1,] <- colMeans(Olen_slope_daily, na.rm = TRUE)
Olen_slope_6hr <- calcslope(rev(xvec), hmat = Olen_WSE_6hr_df)
Olen_slope_6hr <- na.approx(Olen_slope_6hr)
Olen_slope_6hr[1,] <- colMeans(Olen_slope_6hr, na.rm = TRUE)
Olen_slope_hourly <- calcslope(rev(xvec), hmat = Olen_WSE_hourly_df)
Olen_slope_hourly <- na.approx(Olen_slope_hourly)
Olen_slope_hourly[1,] <- colMeans(Olen_slope_hourly, na.rm = TRUE)
Olen_slope_15min <- calcslope(rev(xvec), hmat = Olen_WSE_df)
Olen_slope_15min <- na.approx(Olen_slope_15min)
Olen_slope_15min[1,] <- colMeans(Olen_slope_15min, na.rm = TRUE)

setwd("C:/Users/Merritt/Desktop/Research/PT_Paper/Olentangy")
Olenqobs <- read.csv("Olen_qobs.csv", header = TRUE)[,1]
Olen_qobs_daily<- colMeans(matrix(Olenqobs[1:3648], 96))
Olen_qobs_6hr<- colMeans(matrix(Olenqobs[1:3696], 24))
Olen_qobs_hourly<- colMeans(matrix(Olenqobs[1:3700], 4))
Olen_qobs_15min<- colMeans(matrix(Olenqobs[1:3701], 1))

save(Olen_width_NED_daily, file = "Olentangy_width_NED_daily.RData")
save(Olen_dA_NED_daily, file = "Olentangy_dA_NED_daily.RData")
save(Olen_slope_daily, file = "Olentangy_slope_daily.RData")
save(Olen_qobs_daily, file = "Olentangy_qobs_daily.RData")
save(Olen_width_NED_6hr, file = "Olentangy_width_NED_6hr.RData")
save(Olen_dA_NED_6hr, file = "Olentangy_dA_NED_6hr.RData")
save(Olen_slope_6hr, file = "Olentangy_slope_6hr.RData")
save(Olen_qobs_6hr, file = "Olentangy_qobs_6hr.RData")
save(Olen_width_NED_hourly, file = "Olentangy_width_NED_hourly.RData")
save(Olen_dA_NED_hourly, file = "Olentangy_dA_NED_hourly.RData")
save(Olen_slope_hourly, file = "Olentangy_slope_hourly.RData")
save(Olen_qobs_hourly, file = "Olentangy_qobs_hourly.RData")
save(Olen_width_NED_15min, file = "Olentangy_width_NED_15min.RData")
save(Olen_dA_NED_15min, file = "Olentangy_dA_NED_15min.RData")
save(Olen_slope_15min, file = "Olentangy_slope_15min.RData")
save(Olen_qobs_15min, file = "Olentangy_qobs_15min.RData")

#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Olen_wse_plot <- matplot(t(Olen_WSE_daily_df), type = c("l"), col= colfunc,
                         main = "WSE", xlab = "Days",
                         ylab = "WSE (m)",
                         cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Olen_slope_plot <- matplot(t(Olen_slope_daily), type = c("l"), col= colfunc,
                           main = "Slope", xlab = "Days",
                           ylab = "Slope",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Olen_width_plot <- matplot(Olen_width_NED_daily, type = c("l"), col= colfunc,
                           main = "Width (NED)", xlab = "Days",
                           ylab = "Width (m)",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Olen_width_plot <- matplot(Olen_width_3m_daily, type = c("l"), col= colfunc,
                           main = "Width (Lidar)", xlab = "Days",
                           ylab = "Width (m)",
                           cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Olen_dA_plot <- matplot(t(Olen_dA_NED_daily), type = c("l"), col= colfunc,
                        main = "dA (NED)", xlab = "Days",
                        ylab = "dA (m^2)",
                        cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Olen_dA_plot <- matplot(t(Olen_dA_3m_daily), type = c("l"), col= colfunc,
                        main = "dA (Lidar)", xlab = "Days",
                        ylab = "dA (m^2)",
                        cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
Olen_data_NED_daily <- bam_data(w = t(Olen_width_NED_daily), s = Olen_slope_daily, dA = Olen_dA_NED_daily,
                                Qhat = Olen_qobs_daily)
Olen_NED_daily_default_man_amhg <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "manning_amhg")
Olen_NED_daily_default_amhg <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "amhg")
Olen_NED_daily_default_man <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "manning")
library(hydroGOF)
Olen_NED_daily_default_amhg_val <- bam_val(Olen_NED_daily_default_amhg, Olen_qobs_daily)
Olen_NED_daily_default_man_val <- bam_val(Olen_NED_daily_default_man, Olen_qobs_daily)
Olen_NED_daily_default_man_amhg_val <- bam_val(Olen_NED_daily_default_man_amhg, Olen_qobs_daily)

Olen_data_3m_daily <- bam_data(w = t(Olen_width_3m_daily), s = Olen_slope_daily, dA = Olen_dA_3m_daily,
                               Qhat = Olen_qobs_daily)
Olen_3m_daily_default_man_amhg <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "manning_amhg")
Olen_3m_daily_default_amhg <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "amhg")
Olen_3m_daily_default_man <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "manning")

Olen_3m_daily_default_amhg_val <- bam_val(Olen_3m_daily_default_amhg, Olen_qobs_daily)
Olen_3m_daily_default_man_val <- bam_val(Olen_3m_daily_default_man, Olen_qobs_daily)
Olen_3m_daily_default_man_amhg_val <- bam_val(Olen_3m_daily_default_man_amhg, Olen_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy Daily Default NED",
     ylim = c(0, 50))
lines(Olen_NED_daily_default_amhg_val[[1]], col = "red")
lines(Olen_NED_daily_default_man_val[[1]], col = "blue")
lines(Olen_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Olen_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Olen_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Olen_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Olen_data_NED_daily <- bam_data(w = t(Olen_width_NED_daily), s = Olen_slope_daily, dA = Olen_dA_NED_daily,
                                Qhat = Olen_qobs_daily)
Olen_best_priors <- bam_priors(bamdata= Olen_data_NED_daily, lowerbound_logQ = log(0.01 * min(Olen_qobs_daily)),
                               upperbound_logQ = log(1.75 * max(Olen_qobs_daily)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Olen_slope_daily, na.rm = TRUE),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+03)
Olen_NED_daily_best_man_amhg <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "manning_amhg", 
                                             bampriors = Olen_best_priors)
Olen_NED_daily_best_amhg <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "amhg", bampriors = Olen_best_priors)
Olen_NED_daily_best_man <- bam_estimate(bamdata = Olen_data_NED_daily,variant = "manning", bampriors = Olen_best_priors)
library(hydroGOF)
Olen_NED_daily_best_amhg_val <- bam_val(Olen_NED_daily_best_amhg, Olen_qobs_daily)
Olen_NED_daily_best_man_val <- bam_val(Olen_NED_daily_best_man, Olen_qobs_daily)
Olen_NED_daily_best_man_amhg_val <- bam_val(Olen_NED_daily_best_man_amhg, Olen_qobs_daily)

Olen_data_3m_daily <- bam_data(w = t(Olen_width_3m_daily), s = Olen_slope_daily, dA = Olen_dA_3m_daily,
                               Qhat = Olen_qobs_daily)
Olen_3m_daily_best_man_amhg <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "manning_amhg",
                                            bampriors = Olen_best_priors)
Olen_3m_daily_best_amhg <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "amhg", bampriors = Olen_best_priors)
Olen_3m_daily_best_man <- bam_estimate(bamdata = Olen_data_3m_daily,variant = "manning", bampriors = Olen_best_priors)

Olen_3m_daily_best_amhg_val <- bam_val(Olen_3m_daily_best_amhg, Olen_qobs_daily)
Olen_3m_daily_best_man_val <- bam_val(Olen_3m_daily_best_man, Olen_qobs_daily)
Olen_3m_daily_best_man_amhg_val <- bam_val(Olen_3m_daily_best_man_amhg, Olen_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy Daily best NED",
     ylim = c(0, 15))
lines(Olen_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Olen_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Olen_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Olen_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Olen_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Olen_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
Olen_data_NED_6hr <- bam_data(w = t(Olen_width_NED_6hr), s = Olen_slope_6hr, dA = Olen_dA_NED_6hr,
                              Qhat = Olen_qobs_6hr)
Olen_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "manning_amhg")
Olen_NED_6hr_default_amhg <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "amhg")
Olen_NED_6hr_default_man <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "manning")

Olen_NED_6hr_default_amhg_val <- bam_val(Olen_NED_6hr_default_amhg, Olen_qobs_6hr)
Olen_NED_6hr_default_man_val <- bam_val(Olen_NED_6hr_default_man, Olen_qobs_6hr)
Olen_NED_6hr_default_man_amhg_val <- bam_val(Olen_NED_6hr_default_man_amhg, Olen_qobs_6hr)

Olen_data_3m_6hr <- bam_data(w = t(Olen_width_3m_6hr), s = Olen_slope_6hr, dA = Olen_dA_3m_6hr, 
                             Qhat = Olen_qobs_6hr)
Olen_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "manning_amhg")
Olen_3m_6hr_default_amhg <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "amhg")
Olen_3m_6hr_default_man <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "manning")

Olen_3m_6hr_default_amhg_val <- bam_val(Olen_3m_6hr_default_amhg, Olen_qobs_6hr)
Olen_3m_6hr_default_man_val <- bam_val(Olen_3m_6hr_default_man, Olen_qobs_6hr)
Olen_3m_6hr_default_man_amhg_val <- bam_val(Olen_3m_6hr_default_man_amhg, Olen_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy 6hr Default NED",
     ylim = c(120, 275))
lines(Olen_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Olen_NED_6hr_default_man_val[[1]], col = "blue")
lines(Olen_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Olen_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Olen_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Olen_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Olen_data_NED_6hr <- bam_data(w = t(Olen_width_NED_6hr), s = Olen_slope_6hr, dA = Olen_dA_NED_6hr,
                              Qhat = Olen_qobs_6hr)
Olen_best_priors <- bam_priors(bamdata= Olen_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Olen_qobs_6hr)),
                               upperbound_logQ = log(1.25 * max(Olen_qobs_6hr)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Olen_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Olen_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "manning_amhg", 
                                           bampriors = Olen_best_priors)
Olen_NED_6hr_best_amhg <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "amhg", bampriors = Olen_best_priors)
Olen_NED_6hr_best_man <- bam_estimate(bamdata = Olen_data_NED_6hr,variant = "manning", bampriors = Olen_best_priors)

Olen_NED_6hr_best_amhg_val <- bam_val(Olen_NED_6hr_best_amhg, Olen_qobs_6hr)
Olen_NED_6hr_best_man_val <- bam_val(Olen_NED_6hr_best_man, Olen_qobs_6hr)
Olen_NED_6hr_best_man_amhg_val <- bam_val(Olen_NED_6hr_best_man_amhg, Olen_qobs_6hr)

Olen_data_3m_6hr <- bam_data(w = t(Olen_width_3m_6hr), s = Olen_slope_6hr, dA = Olen_dA_3m_6hr,
                             Qhat = Olen_qobs_6hr)
Olen_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "manning_amhg",
                                          bampriors = Olen_best_priors)
Olen_3m_6hr_best_amhg <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "amhg", bampriors = Olen_best_priors)
Olen_3m_6hr_best_man <- bam_estimate(bamdata = Olen_data_3m_6hr,variant = "manning", bampriors = Olen_best_priors)

Olen_3m_6hr_best_amhg_val <- bam_val(Olen_3m_6hr_best_amhg, Olen_qobs_6hr)
Olen_3m_6hr_best_man_val <- bam_val(Olen_3m_6hr_best_man, Olen_qobs_6hr)
Olen_3m_6hr_best_man_amhg_val <- bam_val(Olen_3m_6hr_best_man_amhg, Olen_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy 6hr best NED",
     ylim = c(135, 155))
lines(Olen_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Olen_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Olen_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Olen_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Olen_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Olen_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Olen_data_NED_hourly <- bam_data(w = t(Olen_width_NED_hourly), s = Olen_slope_hourly, dA = Olen_dA_NED_hourly,
                                 Qhat = Olen_qobs_hourly)
Olen_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "manning_amhg")
Olen_NED_hourly_default_amhg <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "amhg")
Olen_NED_hourly_default_man <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "manning")

Olen_NED_hourly_default_amhg_val <- bam_val(Olen_NED_hourly_default_amhg, Olen_qobs_hourly)
Olen_NED_hourly_default_man_val <- bam_val(Olen_NED_hourly_default_man, Olen_qobs_hourly)
Olen_NED_hourly_default_man_amhg_val <- bam_val(Olen_NED_hourly_default_man_amhg, Olen_qobs_hourly)

Olen_data_3m_hourly <- bam_data(w = t(Olen_width_3m_hourly), s = Olen_slope_hourly, dA = Olen_dA_3m_hourly, 
                                Qhat = Olen_qobs_hourly)
Olen_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "manning_amhg")
Olen_3m_hourly_default_amhg <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "amhg")
Olen_3m_hourly_default_man <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "manning")

Olen_3m_hourly_default_amhg_val <- bam_val(Olen_3m_hourly_default_amhg, Olen_qobs_hourly)
Olen_3m_hourly_default_man_val <- bam_val(Olen_3m_hourly_default_man, Olen_qobs_hourly)
Olen_3m_hourly_default_man_amhg_val <- bam_val(Olen_3m_hourly_default_man_amhg, Olen_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy hourly Default NED",
     ylim = c(120, 275))
lines(Olen_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Olen_NED_hourly_default_man_val[[1]], col = "blue")
lines(Olen_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Olen_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Olen_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Olen_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Olen_data_NED_hourly <- bam_data(w = t(Olen_width_NED_hourly), s = Olen_slope_hourly, dA = Olen_dA_NED_hourly,
                                 Qhat = Olen_qobs_hourly)
Olen_best_priors <- bam_priors(bamdata= Olen_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Olen_qobs_hourly)),
                               upperbound_logQ = log(1.25 * max(Olen_qobs_hourly)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Olen_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Olen_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "manning_amhg", 
                                              bampriors = Olen_best_priors)
Olen_NED_hourly_best_amhg <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "amhg", bampriors = Olen_best_priors)
Olen_NED_hourly_best_man <- bam_estimate(bamdata = Olen_data_NED_hourly,variant = "manning", bampriors = Olen_best_priors)

Olen_NED_hourly_best_amhg_val <- bam_val(Olen_NED_hourly_best_amhg, Olen_qobs_hourly)
Olen_NED_hourly_best_man_val <- bam_val(Olen_NED_hourly_best_man, Olen_qobs_hourly)
Olen_NED_hourly_best_man_amhg_val <- bam_val(Olen_NED_hourly_best_man_amhg, Olen_qobs_hourly)

Olen_data_3m_hourly <- bam_data(w = t(Olen_width_3m_hourly), s = Olen_slope_hourly, dA = Olen_dA_3m_hourly,
                                Qhat = Olen_qobs_hourly)
Olen_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "manning_amhg",
                                             bampriors = Olen_best_priors)
Olen_3m_hourly_best_amhg <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "amhg", bampriors = Olen_best_priors)
Olen_3m_hourly_best_man <- bam_estimate(bamdata = Olen_data_3m_hourly,variant = "manning", bampriors = Olen_best_priors)

Olen_3m_hourly_best_amhg_val <- bam_val(Olen_3m_hourly_best_amhg, Olen_qobs_hourly)
Olen_3m_hourly_best_man_val <- bam_val(Olen_3m_hourly_best_man, Olen_qobs_hourly)
Olen_3m_hourly_best_man_amhg_val <- bam_val(Olen_3m_hourly_best_man_amhg, Olen_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy hourly best NED",
     ylim = c(135, 155))
lines(Olen_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Olen_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Olen_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Olen_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Olen_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Olen_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#15min-----
#Default----
Olen_data_NED_15min <- bam_data(w = t(Olen_width_NED_15min), s = Olen_slope_15min, dA = Olen_dA_NED_15min,
                                Qhat = Olen_qobs_15min)
Olen_NED_15min_default_man_amhg <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "manning_amhg")
Olen_NED_15min_default_amhg <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "amhg")
Olen_NED_15min_default_man <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "manning")

library(hydroGOF)
Olen_NED_15min_default_amhg_val <- bam_val(Olen_NED_15min_default_amhg, Olen_qobs_15min)
Olen_NED_15min_default_man_val <- bam_val(Olen_NED_15min_default_man, Olen_qobs_15min)
Olen_NED_15min_default_man_amhg_val <- bam_val(Olen_NED_15min_default_man_amhg, Olen_qobs_15min)

Olen_data_3m_15min <- bam_data(w = t(Olen_width_3m_15min), s = Olen_slope_15min, dA = Olen_dA_3m_15min, 
                               Qhat = Olen_qobs_15min)
Olen_3m_15min_default_man_amhg <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "manning_amhg")
Olen_3m_15min_default_amhg <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "amhg")
Olen_3m_15min_default_man <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "manning")

library(hydroGOF)
Olen_3m_15min_default_amhg_val <- bam_val(Olen_3m_15min_default_amhg, Olen_qobs_15min)
Olen_3m_15min_default_man_val <- bam_val(Olen_3m_15min_default_man, Olen_qobs_15min)
Olen_3m_15min_default_man_amhg_val <- bam_val(Olen_3m_15min_default_man_amhg, Olen_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy 15min Default NED",
     ylim = c(120, 275))
lines(Olen_NED_15min_default_amhg_val[[1]], col = "red")
lines(Olen_NED_15min_default_man_val[[1]], col = "blue")
lines(Olen_NED_15min_default_man_amhg_val[[1]], col = "purple")
lines(Olen_3m_15min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Olen_3m_15min_default_man_val[[1]], col = "blue", lty = 2)
lines(Olen_3m_15min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Olen_data_NED_15min <- bam_data(w = t(Olen_width_NED_15min), s = Olen_slope_15min, dA = Olen_dA_NED_15min,
                                Qhat = Olen_qobs_15min)
Olen_best_priors <- bam_priors(bamdata= Olen_data_NED_15min, lowerbound_logQ = log(0.75 * min(Olen_qobs_15min)),
                               upperbound_logQ = log(1.25 * max(Olen_qobs_15min)),
                               lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Olen_slope),
                               lowerbound_A0 = 0, 
                               upperbound_A0 = 1e+05)
Olen_NED_15min_best_man_amhg <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "manning_amhg", 
                                             bampriors = Olen_best_priors)
Olen_NED_15min_best_amhg <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "amhg", bampriors = Olen_best_priors)
Olen_NED_15min_best_man <- bam_estimate(bamdata = Olen_data_NED_15min,variant = "manning", bampriors = Olen_best_priors)

Olen_NED_15min_best_amhg_val <- bam_val(Olen_NED_15min_best_amhg, Olen_qobs_15min)
Olen_NED_15min_best_man_val <- bam_val(Olen_NED_15min_best_man, Olen_qobs_15min)
Olen_NED_15min_best_man_amhg_val <- bam_val(Olen_NED_15min_best_man_amhg, Olen_qobs_15min)

Olen_data_3m_15min <- bam_data(w = t(Olen_width_3m_15min), s = Olen_slope_15min, dA = Olen_dA_3m_15min,
                               Qhat = Olen_qobs_15min)
Olen_3m_15min_best_man_amhg <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "manning_amhg",
                                            bampriors = Olen_best_priors)
Olen_3m_15min_best_amhg <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "amhg", bampriors = Olen_best_priors)
Olen_3m_15min_best_man <- bam_estimate(bamdata = Olen_data_3m_15min,variant = "manning", bampriors = Olen_best_priors)

Olen_3m_15min_best_amhg_val <- bam_val(Olen_3m_15min_best_amhg, Olen_qobs_15min)
Olen_3m_15min_best_man_val <- bam_val(Olen_3m_15min_best_man, Olen_qobs_15min)
Olen_3m_15min_best_man_amhg_val <- bam_val(Olen_3m_15min_best_man_amhg, Olen_qobs_15min)

#Plot Hydrograph
par(family = "serif")
plot(Olen_qobs_15min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "OlenOlengy 15min best NED",
     ylim = c(135, 155))
lines(Olen_NED_15min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Olen_NED_15min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Olen_NED_15min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Olen_3m_15min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Olen_3m_15min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Olen_3m_15min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Olen_daily_lower_logQ <- log(seq(exp(max(apply(log(Olen_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                 min(Olen_qobs_15min),length.out = 10))
Olen_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Olen_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_logQ = Olen_daily_lower_logQ[i])
}
#upperbound_logQ
Olen_daily_upper_logQ <- log(seq(max(Olen_qobs_15min), exp(min(apply(log(Olen_width_NED_daily), 2, max)) + log(40) + log(5)),
                                 length.out = 10))
Olen_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Olen_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_logQ = Olen_daily_upper_logQ[i])
}
#lowerbound_A0
Olen_daily_lower_A0 <- seq(0, 30, length.out = 10)
Olen_daily_lower_A0_priors <- list()
for(i in 1:10){
  Olen_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_A0 = Olen_daily_lower_A0[i])
}
#upperbound_A0
Olen_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Olen_daily_upper_A0_priors <- list()
for(i in 1:10){
  Olen_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_A0 = Olen_daily_upper_A0[i])
}
#lowerbound_logn
Olen_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Olen_daily_lower_logn_priors <- list()
for(i in 1:10){
  Olen_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_logn = Olen_daily_lower_logn[i])
}
#upperbound_logn
Olen_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Olen_daily_upper_logn_priors <- list()
for(i in 1:10){
  Olen_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_logn = Olen_daily_upper_logn[i])
}
#lowerbound_logQc
Olen_daily_lower_logQc <- seq(0, log(min(Olen_qobs_15min)),length.out = 10)
Olen_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Olen_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_logQc = Olen_daily_lower_logQc[i])
}
#upperbound_logQc
Olen_daily_upper_logQc <- seq(log(max(Olen_qobs_15min)), 10, length.out = 10)
Olen_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Olen_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_logQc = Olen_daily_upper_logQc[i])
}
#lowerbound_logWc
Olen_daily_lower_logWc <- seq(1, log(min(Olen_width_NED_daily)),length.out = 10)
Olen_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Olen_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_logWc = Olen_daily_lower_logWc[i])
}
#upperbound_logWc
Olen_daily_upper_logWc <- seq(log(max(Olen_width_NED_daily)), 8, length.out = 10)
Olen_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Olen_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_logWc = Olen_daily_upper_logWc[i])
}
#lowerbound_b
Olen_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Olen_daily_lower_b_priors <- list()
for(i in 1:10){
  Olen_daily_lower_b_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, lowerbound_b = Olen_daily_lower_b[i])
}
#upperbound_b
Olen_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Olen_daily_upper_b_priors <- list()
for(i in 1:10){
  Olen_daily_upper_b_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, upperbound_b = Olen_daily_upper_b[i])
}
#sigma_man
Olen_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Olen_daily_sigma_man_priors <- list()
for(i in 1:10){
  Olen_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, sigma_man = Olen_daily_sigma_man[i])
}
#sigma_amhg
Olen_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Olen_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Olen_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, sigma_amhg = Olen_daily_sigma_amhg[i])
}
#logQc_hat
Olen_daily_logQc_hat <- log(seq(min(Olen_qobs_15min), max(Olen_qobs_15min), length.out = 10))
Olen_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Olen_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logQc_hat = Olen_daily_logQc_hat[i])
}
#logWc_hat
Olen_daily_logWc_hat <- log(seq(min(Olen_width_NED_daily), max(Olen_width_NED_daily), length.out = 10))
Olen_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Olen_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logWc_hat = Olen_daily_logWc_hat[i])
}
#b_hat
Olen_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Olen_daily_b_hat_priors <- list()
for(i in 1:10){
  Olen_daily_b_hat_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, b_hat = Olen_daily_b_hat[i])
}
#logA0_hat
Olen_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Olen_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Olen_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logA0_hat = Olen_daily_logA0_hat[i])
}
#logn_hat
Olen_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Olen_daily_logn_hat_priors <- list()
for(i in 1:10){
  Olen_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logn_hat = Olen_daily_logn_hat[i])
}
#logQ_sd
Olen_daily_logQ_sd <- seq(sd(log(Olen_qobs_15min)), 0.8325546, length.out = 10)
Olen_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Olen_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logQ_sd = Olen_daily_logQ_sd[i])
}
#logQc_sd
Olen_daily_logQc_sd <- seq(sd(log(Olen_qobs_15min)), 0.8325546, length.out = 10)
Olen_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Olen_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logQc_sd = Olen_daily_logQc_sd[i])
}
#logWc_Sd
Olen_daily_logWc_Sd <- seq(log(sd(Olen_width_NED_daily)), 4.712493, length.out = 10)
Olen_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Olen_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logWc_sd = Olen_daily_logWc_Sd[i])
}
#b_sd
Olen_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Olen_daily_b_sd_priors <- list()
for(i in 1:10){
  Olen_daily_b_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, b_sd = Olen_daily_b_sd[i])
}
#logA0_sd
Olen_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Olen_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Olen_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logA0_sd = Olen_daily_logA0_sd[i])
}
#logn_sd
Olen_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Olen_daily_logn_sd_priors <- list()
for(i in 1:10){
  Olen_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, logn_sd = Olen_daily_logn_sd[i])
}
#Werr_sd
Olen_daily_Werr_sd <- seq(5, 15, length.out = 10)
Olen_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Olen_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, Werr_sd = Olen_daily_Werr_sd[i])
}
#Serr_sd
Olen_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Olen_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Olen_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, Serr_sd = Olen_daily_Serr_sd[i])
}
#dAerr_sd
Olen_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Olen_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Olen_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Olen_data_NED_daily, dAerr_sd = Olen_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Olen_daily_manning_list <- list(Olen_daily_lower_logQ_priors,
                                Olen_daily_upper_logQ_priors,
                                Olen_daily_lower_A0_priors,
                                Olen_daily_upper_A0_priors,
                                Olen_daily_lower_logn_priors,
                                Olen_daily_upper_logn_priors,
                                Olen_daily_sigma_man_priors,
                                Olen_daily_logA0_hat_priors,
                                Olen_daily_logn_hat_priors,
                                Olen_daily_logQ_sd_priors,
                                Olen_daily_logA0_sd_priors,
                                Olen_daily_logn_sd_priors,
                                Olen_daily_Werr_sd_priors,
                                Olen_daily_Serr_sd_priors,
                                Olen_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Olen_data_NED_daily, variant = "manning", bampriors = Olen_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Olen_daily_amhg_list <- list(Olen_daily_lower_logQ_priors,
                             Olen_daily_upper_logQ_priors,
                             Olen_daily_lower_logQc_priors,
                             Olen_daily_upper_logQc_priors,
                             Olen_daily_lower_logWc_priors,
                             Olen_daily_upper_logWc_priors,
                             Olen_daily_lower_b_priors,
                             Olen_daily_upper_b_priors,
                             Olen_daily_sigma_amhg_priors,
                             Olen_daily_logQc_hat_priors,
                             Olen_daily_logWc_hat_priors,
                             Olen_daily_b_hat_priors,
                             Olen_daily_logQ_sd_priors,
                             Olen_daily_logQc_sd_priors,
                             Olen_daily_logWc_Sd_priors,
                             Olen_daily_b_sd_priors,
                             Olen_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Olen_data_NED_daily, variant = "amhg", bampriors = Olen_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Olen_daily_man_amhg_list <- list(Olen_daily_lower_logQ_priors,
                                 Olen_daily_upper_logQ_priors,
                                 Olen_daily_lower_A0_priors,
                                 Olen_daily_upper_A0_priors,
                                 Olen_daily_lower_logn_priors,
                                 Olen_daily_upper_logn_priors,
                                 Olen_daily_lower_logQc_priors,
                                 Olen_daily_upper_logQc_priors,
                                 Olen_daily_lower_logWc_priors,
                                 Olen_daily_upper_logWc_priors,
                                 Olen_daily_lower_b_priors,
                                 Olen_daily_upper_b_priors,
                                 Olen_daily_sigma_man_priors,
                                 Olen_daily_sigma_amhg_priors,
                                 Olen_daily_logQc_hat_priors,
                                 Olen_daily_logWc_hat_priors,
                                 Olen_daily_b_hat_priors,
                                 Olen_daily_logA0_hat_priors,
                                 Olen_daily_logn_hat_priors,
                                 Olen_daily_logQ_sd_priors,
                                 Olen_daily_logQc_sd_priors,
                                 Olen_daily_logWc_Sd_priors,
                                 Olen_daily_b_sd_priors,
                                 Olen_daily_logA0_sd_priors,
                                 Olen_daily_logn_sd_priors,
                                 Olen_daily_Werr_sd_priors,
                                 Olen_daily_Serr_sd_priors,
                                 Olen_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/OlenOlengy/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Olen_data_NED_daily, variant = "manning_amhg", bampriors = Olen_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

###Tanana
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/Tanana"))
load("Tan_WSE_daily.RData")
load("Tan_WSE_6hr.RData")
load("Tan_WSE_hourly.RData")
load("Tan_WSE_2min.RData")

Tan_PT_coords <- read.csv("Tanana_PT_coord.csv")
Tan_PT_coords_sp <- SpatialPoints(Tan_PT_coords[,c(3,4)], proj4string =  CRS("+proj=utm +zone=6 +datum=WGS84"))
Tan_cl <-readOGR(dsn = ".", layer = "Tan_cl")

setwd(paste0(home_dir, "/Tanana/DEM"))
Tan_DEM_NED <- raster("Tan_DEM_2m.tif")
Tan_DEM_3m <- raster()

setwd(paste0(home_dir, "/Tanana"))
Tan_orthos_NED <- ortho_lines(Tan_DEM_NED, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords), 1000)
Tan_dist_mat_NED <- pointDistance(Tan_PT_coords_sp, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords))
Tan_closest_ortho_NED <- Tan_orthos_NED[apply(Tan_dist_mat_NED, 1, FUN = which.min)]

Tan_orthos_3m <- ortho_lines(Tan_DEM_3m, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords), 1000)
Tan_dist_mat_3m <- pointDistance(Tan_PT_coords_sp, data.frame(Tan_cl@lines[[1]]@Lines[[1]]@coords))
Tan_closest_ortho_3m <- Tan_orthos_3m[apply(Tan_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
Tan_elev_prof_NED <- foreach(i=1:length(Tan_closest_ortho_NED), 
                             .combine='rbind', .packages=c('raster')) %dopar% {
                               extract(Tan_DEM_NED, Tan_closest_ortho_NED[[i]], method = 'simple')
                             }
Tan_elev_prof_3m <- foreach(i=1:length(Tan_closest_ortho_3m), 
                            .combine='rbind', .packages=c('raster')) %dopar% {
                              extract(Tan_DEM_3m, Tan_closest_ortho_3m[[i]], method = 'simple')
                            }
proc.time()-ptm
stopCluster(cl)

Tan_scale_x_NED <- vector()
for(i in 1:length(Tan_elev_prof_NED)){
  Tan_scale_x_NED[i] = 2000/length(Tan_elev_prof_NED[[i]])
}

Tan_scale_x_3m <- vector()
for(i in 1:length(Tan_elev_prof_3m)){
  Tan_scale_x_3m[i] = 2000/length(Tan_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(Tan_elev_prof_NED, Tan_WSE_daily_df, "Tanana")
locator(n = 4)
Tan_upper_l_NED <- c( 164.8,  479.0, 600.0,  408.5,  466.5,  406.7,  362.3, 363.7, 381.7, 471.5, 437.2, 533.2)
Tan_lower_l_NED <- c( 241.9,  533.1, 715.8,  490.1,  562.3,  514.6,  417.6, 409.1, 441.7, 556.0, 483.2, 572.3)
Tan_lower_r_NED <- c(1079.9, 1047.9, 875.4, 1416.2, 1278.0,  918.5,  965.5, 802.8, 710.3, 823.2, 682.2, 759.4)
Tan_upper_r_NED <- c(1159.7, 1117.4, 947.5, 1500.4, 1361.2, 1089.1, 1046.0, 858.3, 747.8, 881.8, 717.2, 818.1)

Tan_elev_prof_NED_corr <- Tan_elev_prof_NED
#Tan_elev_prof_NED_corr[[12]] <- Tan_elev_prof_NED[[13]]
#Tan_scale_x_NED[12] <- Tan_scale_x_NED[13]

Tan_width_NED_daily <- calc_width(h = t(Tan_WSE_daily_df), elev_prof = Tan_elev_prof_NED_corr, 
                                  upper_l = Tan_upper_l_NED, lower_l = Tan_lower_l_NED,
                                  lower_r = Tan_lower_r_NED, upper_r = Tan_upper_r_NED,
                                  scale_prof = Tan_scale_x_NED)
Tan_width_NED_daily <- na.approx(Tan_width_NED_daily)
Tan_dA_NED_daily <- calcdA_mat(w = t(Tan_width_NED_daily), h = Tan_WSE_daily_df)

Tan_width_NED_6hr <- calc_width(h = t(Tan_WSE_6hr_df), elev_prof = Tan_elev_prof_NED_corr, 
                                upper_l = Tan_upper_l_NED, lower_l = Tan_lower_l_NED,
                                lower_r = Tan_lower_r_NED, upper_r = Tan_upper_r_NED,
                                scale_prof = Tan_scale_x_NED)
Tan_width_NED_6hr <- na.approx(Tan_width_NED_6hr)
Tan_dA_NED_6hr <- calcdA_mat(w = t(Tan_width_NED_6hr), h = Tan_WSE_6hr_df)

Tan_width_NED_hourly <- calc_width(h = t(Tan_WSE_hourly_df), elev_prof = Tan_elev_prof_NED_corr, 
                                   upper_l = Tan_upper_l_NED, lower_l = Tan_lower_l_NED,
                                   lower_r = Tan_lower_r_NED, upper_r = Tan_upper_r_NED,
                                   scale_prof = Tan_scale_x_NED)
Tan_width_NED_hourly <- na.approx(Tan_width_NED_hourly)
Tan_dA_NED_hourly <- calcdA_mat(w = t(Tan_width_NED_hourly), h = Tan_WSE_hourly_df)

Tan_width_NED_2min <- calc_width(h = t(Tan_WSE_df), elev_prof = Tan_elev_prof_NED_corr, 
                                 upper_l = Tan_upper_l_NED, lower_l = Tan_lower_l_NED,
                                 lower_r = Tan_lower_r_NED, upper_r = Tan_upper_r_NED,
                                 scale_prof = Tan_scale_x_NED)
Tan_width_NED_2min <- na.approx(Tan_width_NED_2min)
Tan_dA_NED_2min <- calcdA_mat(w = t(Tan_width_NED_2min), h = Tan_WSE_df)

plot_elev_prof(Tan_elev_prof_3m, Tan_WSE_daily_df, "Tanana Reach Three")
locator(n = 4)
Tan_upper_l_3m <- c(157.7, 169.0, 117.5, 201.2, 143.3, 129.0, 143.9,  33.0, 132.9, 230.4, 182.2, 211.3,  88.1, 155.0,  63.5,  90.4, 197.4, 203.8, 106.4)
Tan_lower_l_3m <- c(175.8, 176.1, 123.6, 216.4, 156.9, 136.9, 150.7,  54.5, 159.9, 237.0, 191.0, 226.7,  94.7, 168.1,  72.7, 103.3, 204.3, 209.4, 130.2)
Tan_lower_r_3m <- c(209.0, 218.0, 216.5, 303.7, 353.8, 343.5, 209.0, 455.0, 390.5, 376.9, 262.4, 309.3, 178.2, 275.2, 179.1, 160.2, 253.8, 242.9, 262.4)
Tan_upper_r_3m <- c(231.2, 226.1, 225.0, 321.9, 371.9, 356.9, 225.0, 468.1, 401.4, 388.1, 283.4, 323.1, 184.8, 285.0, 207.6, 175.4, 266.2, 252.2, 275.2)

Tan_width_3m_daily <- calc_width(h = t(Tan_WSE_daily_df), elev_prof = Tan_elev_prof_3m, 
                                 upper_l = Tan_upper_l_3m, lower_l = Tan_lower_l_3m,
                                 lower_r = Tan_lower_r_3m, upper_r = Tan_upper_r_3m,
                                 scale_prof = Tan_scale_x_3m)
Tan_width_3m_daily <- na.approx(Tan_width_3m_daily)
Tan_dA_3m_daily <- calcdA_mat(w = t(Tan_width_3m_daily), h = Tan_WSE_daily_df)

Tan_width_3m_6hr <- calc_width(h = t(Tan_WSE_6hr_df), elev_prof = Tan_elev_prof_3m, 
                               upper_l = Tan_upper_l_3m, lower_l = Tan_lower_l_3m,
                               lower_r = Tan_lower_r_3m, upper_r = Tan_upper_r_3m,
                               scale_prof = Tan_scale_x_3m)
Tan_width_3m_6hr <- na.approx(Tan_width_3m_6hr)
Tan_dA_3m_6hr <- calcdA_mat(w = t(Tan_width_3m_6hr), h = Tan_WSE_6hr_df)

Tan_width_3m_hourly <- calc_width(h = t(Tan_WSE_hourly_df), elev_prof = Tan_elev_prof_3m, 
                                  upper_l = Tan_upper_l_3m, lower_l = Tan_lower_l_3m,
                                  lower_r = Tan_lower_r_3m, upper_r = Tan_upper_r_3m,
                                  scale_prof = Tan_scale_x_3m)
Tan_width_3m_hourly <- na.approx(Tan_width_3m_hourly)
Tan_dA_3m_hourly <- calcdA_mat(w = t(Tan_width_3m_hourly), h = Tan_WSE_hourly_df)

Tan_width_3m_2min <- calc_width(h = t(Tan_WSE_df), elev_prof = Tan_elev_prof_3m, 
                                upper_l = Tan_upper_l_3m, lower_l = Tan_lower_l_3m,
                                lower_r = Tan_lower_r_3m, upper_r = Tan_upper_r_3m,
                                scale_prof = Tan_scale_x_3m)
Tan_width_3m_2min <- na.approx(Tan_width_3m_2min)
Tan_dA_3m_2min <- calcdA_mat(w = t(Tan_width_3m_2min), h = Tan_WSE_df)


xvec <- calc_xvec(Tan_PT_coords_sp, Tan_cl)
xvec = rev(xvec)
Tan_xvec <- saveRDS(xvec, file = "Tan_xvec.rds")
Tan_slope_daily <- calcslope(rev(xvec), hmat = Tan_WSE_daily_df)
Tan_slope_daily <- na.replace(Tan_slope_daily, mean(Tan_slope_daily))
Tan_slope_daily[1,] <- colMeans(Tan_slope_daily, na.rm = TRUE)
Tan_slope_6hr <- calcslope(rev(xvec), hmat = Tan_WSE_6hr_df)
Tan_slope_6hr <- na.approx(Tan_slope_6hr)
Tan_slope_6hr[1,] <- colMeans(Tan_slope_6hr, na.rm = TRUE)
Tan_slope_hourly <- calcslope(rev(xvec), hmat = Tan_WSE_hourly_df)
Tan_slope_hourly <- na.approx(Tan_slope_hourly)
Tan_slope_hourly[1,] <- colMeans(Tan_slope_hourly, na.rm = TRUE)
Tan_slope_2min <- calcslope(rev(xvec), hmat = Tan_WSE_df)
Tan_slope_2min <- na.approx(Tan_slope_2min)
Tan_slope_2min[1,] <- colMeans(Tan_slope_2min, na.rm = TRUE)

setwd("C:/Users/Merritt/Desktop/Research/PT_Paper/Tanana")
Tanqobs <- read.csv("TanQobs.csv", header = TRUE)$Q*0.028316847
Tan_qobs_daily<- colMeans(matrix(Tanqobs[1:2976], 96))
Tan_qobs_6hr<- colMeans(matrix(Tanqobs[1:3048], 24))
Tan_qobs_hourly<- colMeans(matrix(Tanqobs[1:3056], 4))

save(Tan_width_NED_daily, file = "Tanana_width_NED_daily.RData")
save(Tan_dA_NED_daily, file = "Tanana_dA_NED_daily.RData")
save(Tan_slope_daily, file = "Tanana_slope_daily.RData")
save(Tan_qobs_daily, file = "Tanana_qobs_daily.RData")
save(Tan_width_NED_6hr, file = "Tanana_width_NED_6hr.RData")
save(Tan_dA_NED_6hr, file = "Tanana_dA_NED_6hr.RData")
save(Tan_slope_6hr, file = "Tanana_slope_6hr.RData")
save(Tan_qobs_6hr, file = "Tanana_qobs_6hr.RData")
save(Tan_width_NED_hourly, file = "Tanana_width_NED_hourly.RData")
save(Tan_dA_NED_hourly, file = "Tanana_dA_NED_hourly.RData")
save(Tan_slope_hourly, file = "Tanana_slope_hourly.RData")
save(Tan_qobs_hourly, file = "Tanana_qobs_hourly.RData")
save(Tan_width_NED_15min, file = "Tanana_width_NED_15min.RData")
save(Tan_dA_NED_15min, file = "Tanana_dA_NED_15min.RData")
save(Tan_slope_15min, file = "Tanana_slope_15min.RData")
save(Tan_qobs_15min, file = "Tanana_qobs_15min.RData")
#Tan_qobs_2min<- colMeans(matrix(Tanqobs[1:3701], 1))

#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Tan_wse_plot <- matplot(t(Tan_WSE_daily_df), type = c("l"), col= colfunc,
                        main = "WSE", xlab = "Days",
                        ylab = "WSE (m)",
                        cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Tan_slope_plot <- matplot(t(Tan_slope_daily), type = c("l"), col= colfunc,
                          main = "Slope", xlab = "Days",
                          ylab = "Slope",
                          cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Tan_width_plot <- matplot(Tan_width_NED_daily, type = c("l"), col= colfunc,
                          main = "Width (NED)", xlab = "Days",
                          ylab = "Width (m)",
                          cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Tan_width_plot <- matplot(Tan_width_3m_daily, type = c("l"), col= colfunc,
                          main = "Width (Lidar)", xlab = "Days",
                          ylab = "Width (m)",
                          cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Tan_dA_plot <- matplot(t(Tan_dA_NED_daily), type = c("l"), col= colfunc,
                       main = "dA (NED)", xlab = "Days",
                       ylab = "dA (m^2)",
                       cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
Tan_dA_plot <- matplot(t(Tan_dA_3m_daily), type = c("l"), col= colfunc,
                       main = "dA (Lidar)", xlab = "Days",
                       ylab = "dA (m^2)",
                       cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
Tan_data_NED_daily <- bam_data(w = t(Tan_width_NED_daily), s = Tan_slope_daily, dA = Tan_dA_NED_daily,
                               Qhat = Tan_qobs_daily)
Tan_NED_daily_default_man_amhg <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "manning_amhg")
Tan_NED_daily_default_amhg <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "amhg")
Tan_NED_daily_default_man <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "manning")
library(hydroGOF)
Tan_NED_daily_default_amhg_val <- bam_val(Tan_NED_daily_default_amhg, Tan_qobs_daily)
Tan_NED_daily_default_man_val <- bam_val(Tan_NED_daily_default_man, Tan_qobs_daily)
Tan_NED_daily_default_man_amhg_val <- bam_val(Tan_NED_daily_default_man_amhg, Tan_qobs_daily)

Tan_data_3m_daily <- bam_data(w = t(Tan_width_3m_daily), s = Tan_slope_daily, dA = Tan_dA_3m_daily,
                              Qhat = Tan_qobs_daily)
Tan_3m_daily_default_man_amhg <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "manning_amhg")
Tan_3m_daily_default_amhg <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "amhg")
Tan_3m_daily_default_man <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "manning")

Tan_3m_daily_default_amhg_val <- bam_val(Tan_3m_daily_default_amhg, Tan_qobs_daily)
Tan_3m_daily_default_man_val <- bam_val(Tan_3m_daily_default_man, Tan_qobs_daily)
Tan_3m_daily_default_man_amhg_val <- bam_val(Tan_3m_daily_default_man_amhg, Tan_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana Daily Default NED",
     ylim = c(300, 2000))
lines(Tan_NED_daily_default_amhg_val[[1]], col = "red")
lines(Tan_NED_daily_default_man_val[[1]], col = "blue")
lines(Tan_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(Tan_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(Tan_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(Tan_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Tan_data_NED_daily <- bam_data(w = t(Tan_width_NED_daily), s = Tan_slope_daily, dA = Tan_dA_NED_daily,
                               Qhat = Tan_qobs_daily)
Tan_best_priors <- bam_priors(bamdata= Tan_data_NED_daily, lowerbound_logQ = log(0.75 * min(Tan_qobs_daily)),
                              upperbound_logQ = log(1.1 * max(Tan_qobs_daily)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Tan_slope_daily, na.rm = TRUE),
                              lowerbound_A0 = 100, 
                              upperbound_A0 = 1e+06)
Tan_NED_daily_best_man_amhg <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "manning_amhg", 
                                            bampriors = Tan_best_priors)
Tan_NED_daily_best_amhg <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "amhg", bampriors = Tan_best_priors)
Tan_NED_daily_best_man <- bam_estimate(bamdata = Tan_data_NED_daily,variant = "manning", bampriors = Tan_best_priors)
library(hydroGOF)
Tan_NED_daily_best_amhg_val <- bam_val(Tan_NED_daily_best_amhg, Tan_qobs_daily)
Tan_NED_daily_best_man_val <- bam_val(Tan_NED_daily_best_man, Tan_qobs_daily)
Tan_NED_daily_best_man_amhg_val <- bam_val(Tan_NED_daily_best_man_amhg, Tan_qobs_daily)

Tan_data_3m_daily <- bam_data(w = t(Tan_width_3m_daily), s = Tan_slope_daily, dA = Tan_dA_3m_daily,
                              Qhat = Tan_qobs_daily)
Tan_3m_daily_best_man_amhg <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "manning_amhg",
                                           bampriors = Tan_best_priors)
Tan_3m_daily_best_amhg <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "amhg", bampriors = Tan_best_priors)
Tan_3m_daily_best_man <- bam_estimate(bamdata = Tan_data_3m_daily,variant = "manning", bampriors = Tan_best_priors)

Tan_3m_daily_best_amhg_val <- bam_val(Tan_3m_daily_best_amhg, Tan_qobs_daily)
Tan_3m_daily_best_man_val <- bam_val(Tan_3m_daily_best_man, Tan_qobs_daily)
Tan_3m_daily_best_man_amhg_val <- bam_val(Tan_3m_daily_best_man_amhg, Tan_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana Daily best NED",
     ylim = c(300, 2000))
lines(Tan_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Tan_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(Tan_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Tan_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Tan_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Tan_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
Tan_data_NED_6hr <- bam_data(w = t(Tan_width_NED_6hr), s = Tan_slope_6hr, dA = Tan_dA_NED_6hr,
                             Qhat = Tan_qobs_6hr)
Tan_NED_6hr_default_man_amhg <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "manning_amhg")
Tan_NED_6hr_default_amhg <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "amhg")
Tan_NED_6hr_default_man <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "manning")

Tan_NED_6hr_default_amhg_val <- bam_val(Tan_NED_6hr_default_amhg, Tan_qobs_6hr)
Tan_NED_6hr_default_man_val <- bam_val(Tan_NED_6hr_default_man, Tan_qobs_6hr)
Tan_NED_6hr_default_man_amhg_val <- bam_val(Tan_NED_6hr_default_man_amhg, Tan_qobs_6hr)

Tan_data_3m_6hr <- bam_data(w = t(Tan_width_3m_6hr), s = Tan_slope_6hr, dA = Tan_dA_3m_6hr, 
                            Qhat = Tan_qobs_6hr)
Tan_3m_6hr_default_man_amhg <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "manning_amhg")
Tan_3m_6hr_default_amhg <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "amhg")
Tan_3m_6hr_default_man <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "manning")

Tan_3m_6hr_default_amhg_val <- bam_val(Tan_3m_6hr_default_amhg, Tan_qobs_6hr)
Tan_3m_6hr_default_man_val <- bam_val(Tan_3m_6hr_default_man, Tan_qobs_6hr)
Tan_3m_6hr_default_man_amhg_val <- bam_val(Tan_3m_6hr_default_man_amhg, Tan_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana 6hr Default NED",
     ylim = c(120, 275))
lines(Tan_NED_6hr_default_amhg_val[[1]], col = "red")
lines(Tan_NED_6hr_default_man_val[[1]], col = "blue")
lines(Tan_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(Tan_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(Tan_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(Tan_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Tan_data_NED_6hr <- bam_data(w = t(Tan_width_NED_6hr), s = Tan_slope_6hr, dA = Tan_dA_NED_6hr,
                             Qhat = Tan_qobs_6hr)
Tan_best_priors <- bam_priors(bamdata= Tan_data_NED_6hr, lowerbound_logQ = log(0.75 * min(Tan_qobs_6hr)),
                              upperbound_logQ = log(1.25 * max(Tan_qobs_6hr)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Tan_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Tan_NED_6hr_best_man_amhg <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "manning_amhg", 
                                          bampriors = Tan_best_priors)
Tan_NED_6hr_best_amhg <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "amhg", bampriors = Tan_best_priors)
Tan_NED_6hr_best_man <- bam_estimate(bamdata = Tan_data_NED_6hr,variant = "manning", bampriors = Tan_best_priors)

Tan_NED_6hr_best_amhg_val <- bam_val(Tan_NED_6hr_best_amhg, Tan_qobs_6hr)
Tan_NED_6hr_best_man_val <- bam_val(Tan_NED_6hr_best_man, Tan_qobs_6hr)
Tan_NED_6hr_best_man_amhg_val <- bam_val(Tan_NED_6hr_best_man_amhg, Tan_qobs_6hr)

Tan_data_3m_6hr <- bam_data(w = t(Tan_width_3m_6hr), s = Tan_slope_6hr, dA = Tan_dA_3m_6hr,
                            Qhat = Tan_qobs_6hr)
Tan_3m_6hr_best_man_amhg <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "manning_amhg",
                                         bampriors = Tan_best_priors)
Tan_3m_6hr_best_amhg <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "amhg", bampriors = Tan_best_priors)
Tan_3m_6hr_best_man <- bam_estimate(bamdata = Tan_data_3m_6hr,variant = "manning", bampriors = Tan_best_priors)

Tan_3m_6hr_best_amhg_val <- bam_val(Tan_3m_6hr_best_amhg, Tan_qobs_6hr)
Tan_3m_6hr_best_man_val <- bam_val(Tan_3m_6hr_best_man, Tan_qobs_6hr)
Tan_3m_6hr_best_man_amhg_val <- bam_val(Tan_3m_6hr_best_man_amhg, Tan_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana 6hr best NED",
     ylim = c(135, 155))
lines(Tan_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Tan_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(Tan_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Tan_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Tan_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Tan_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
Tan_data_NED_hourly <- bam_data(w = t(Tan_width_NED_hourly), s = Tan_slope_hourly, dA = Tan_dA_NED_hourly,
                                Qhat = Tan_qobs_hourly)
Tan_NED_hourly_default_man_amhg <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "manning_amhg")
Tan_NED_hourly_default_amhg <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "amhg")
Tan_NED_hourly_default_man <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "manning")

Tan_NED_hourly_default_amhg_val <- bam_val(Tan_NED_hourly_default_amhg, Tan_qobs_hourly)
Tan_NED_hourly_default_man_val <- bam_val(Tan_NED_hourly_default_man, Tan_qobs_hourly)
Tan_NED_hourly_default_man_amhg_val <- bam_val(Tan_NED_hourly_default_man_amhg, Tan_qobs_hourly)

Tan_data_3m_hourly <- bam_data(w = t(Tan_width_3m_hourly), s = Tan_slope_hourly, dA = Tan_dA_3m_hourly, 
                               Qhat = Tan_qobs_hourly)
Tan_3m_hourly_default_man_amhg <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "manning_amhg")
Tan_3m_hourly_default_amhg <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "amhg")
Tan_3m_hourly_default_man <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "manning")

Tan_3m_hourly_default_amhg_val <- bam_val(Tan_3m_hourly_default_amhg, Tan_qobs_hourly)
Tan_3m_hourly_default_man_val <- bam_val(Tan_3m_hourly_default_man, Tan_qobs_hourly)
Tan_3m_hourly_default_man_amhg_val <- bam_val(Tan_3m_hourly_default_man_amhg, Tan_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana hourly Default NED",
     ylim = c(120, 275))
lines(Tan_NED_hourly_default_amhg_val[[1]], col = "red")
lines(Tan_NED_hourly_default_man_val[[1]], col = "blue")
lines(Tan_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(Tan_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(Tan_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(Tan_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Tan_data_NED_hourly <- bam_data(w = t(Tan_width_NED_hourly), s = Tan_slope_hourly, dA = Tan_dA_NED_hourly,
                                Qhat = Tan_qobs_hourly)
Tan_best_priors <- bam_priors(bamdata= Tan_data_NED_hourly, lowerbound_logQ = log(0.75 * min(Tan_qobs_hourly)),
                              upperbound_logQ = log(1.25 * max(Tan_qobs_hourly)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Tan_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Tan_NED_hourly_best_man_amhg <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "manning_amhg", 
                                             bampriors = Tan_best_priors)
Tan_NED_hourly_best_amhg <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "amhg", bampriors = Tan_best_priors)
Tan_NED_hourly_best_man <- bam_estimate(bamdata = Tan_data_NED_hourly,variant = "manning", bampriors = Tan_best_priors)

Tan_NED_hourly_best_amhg_val <- bam_val(Tan_NED_hourly_best_amhg, Tan_qobs_hourly)
Tan_NED_hourly_best_man_val <- bam_val(Tan_NED_hourly_best_man, Tan_qobs_hourly)
Tan_NED_hourly_best_man_amhg_val <- bam_val(Tan_NED_hourly_best_man_amhg, Tan_qobs_hourly)

Tan_data_3m_hourly <- bam_data(w = t(Tan_width_3m_hourly), s = Tan_slope_hourly, dA = Tan_dA_3m_hourly,
                               Qhat = Tan_qobs_hourly)
Tan_3m_hourly_best_man_amhg <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "manning_amhg",
                                            bampriors = Tan_best_priors)
Tan_3m_hourly_best_amhg <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "amhg", bampriors = Tan_best_priors)
Tan_3m_hourly_best_man <- bam_estimate(bamdata = Tan_data_3m_hourly,variant = "manning", bampriors = Tan_best_priors)

Tan_3m_hourly_best_amhg_val <- bam_val(Tan_3m_hourly_best_amhg, Tan_qobs_hourly)
Tan_3m_hourly_best_man_val <- bam_val(Tan_3m_hourly_best_man, Tan_qobs_hourly)
Tan_3m_hourly_best_man_amhg_val <- bam_val(Tan_3m_hourly_best_man_amhg, Tan_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana hourly best NED",
     ylim = c(135, 155))
lines(Tan_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Tan_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(Tan_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Tan_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Tan_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Tan_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#2min-----
#Default----
Tan_data_NED_2min <- bam_data(w = t(Tan_width_NED_2min), s = Tan_slope_2min, dA = Tan_dA_NED_2min,
                              Qhat = Tan_qobs_2min)
Tan_NED_2min_default_man_amhg <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "manning_amhg")
Tan_NED_2min_default_amhg <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "amhg")
Tan_NED_2min_default_man <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "manning")

library(hydroGOF)
Tan_NED_2min_default_amhg_val <- bam_val(Tan_NED_2min_default_amhg, Tan_qobs_2min)
Tan_NED_2min_default_man_val <- bam_val(Tan_NED_2min_default_man, Tan_qobs_2min)
Tan_NED_2min_default_man_amhg_val <- bam_val(Tan_NED_2min_default_man_amhg, Tan_qobs_2min)

Tan_data_3m_2min <- bam_data(w = t(Tan_width_3m_2min), s = Tan_slope_2min, dA = Tan_dA_3m_2min, 
                             Qhat = Tan_qobs_2min)
Tan_3m_2min_default_man_amhg <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "manning_amhg")
Tan_3m_2min_default_amhg <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "amhg")
Tan_3m_2min_default_man <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "manning")

library(hydroGOF)
Tan_3m_2min_default_amhg_val <- bam_val(Tan_3m_2min_default_amhg, Tan_qobs_2min)
Tan_3m_2min_default_man_val <- bam_val(Tan_3m_2min_default_man, Tan_qobs_2min)
Tan_3m_2min_default_man_amhg_val <- bam_val(Tan_3m_2min_default_man_amhg, Tan_qobs_2min)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_2min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana 2min Default NED",
     ylim = c(120, 275))
lines(Tan_NED_2min_default_amhg_val[[1]], col = "red")
lines(Tan_NED_2min_default_man_val[[1]], col = "blue")
lines(Tan_NED_2min_default_man_amhg_val[[1]], col = "purple")
lines(Tan_3m_2min_default_amhg_val[[1]], col = "red", lty = 2)
lines(Tan_3m_2min_default_man_val[[1]], col = "blue", lty = 2)
lines(Tan_3m_2min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
Tan_data_NED_2min <- bam_data(w = t(Tan_width_NED_2min), s = Tan_slope_2min, dA = Tan_dA_NED_2min,
                              Qhat = Tan_qobs_2min)
Tan_best_priors <- bam_priors(bamdata= Tan_data_NED_2min, lowerbound_logQ = log(0.75 * min(Tan_qobs_2min)),
                              upperbound_logQ = log(1.25 * max(Tan_qobs_2min)),
                              lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(Tan_slope),
                              lowerbound_A0 = 0, 
                              upperbound_A0 = 1e+05)
Tan_NED_2min_best_man_amhg <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "manning_amhg", 
                                           bampriors = Tan_best_priors)
Tan_NED_2min_best_amhg <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "amhg", bampriors = Tan_best_priors)
Tan_NED_2min_best_man <- bam_estimate(bamdata = Tan_data_NED_2min,variant = "manning", bampriors = Tan_best_priors)

Tan_NED_2min_best_amhg_val <- bam_val(Tan_NED_2min_best_amhg, Tan_qobs_2min)
Tan_NED_2min_best_man_val <- bam_val(Tan_NED_2min_best_man, Tan_qobs_2min)
Tan_NED_2min_best_man_amhg_val <- bam_val(Tan_NED_2min_best_man_amhg, Tan_qobs_2min)

Tan_data_3m_2min <- bam_data(w = t(Tan_width_3m_2min), s = Tan_slope_2min, dA = Tan_dA_3m_2min,
                             Qhat = Tan_qobs_2min)
Tan_3m_2min_best_man_amhg <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "manning_amhg",
                                          bampriors = Tan_best_priors)
Tan_3m_2min_best_amhg <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "amhg", bampriors = Tan_best_priors)
Tan_3m_2min_best_man <- bam_estimate(bamdata = Tan_data_3m_2min,variant = "manning", bampriors = Tan_best_priors)

Tan_3m_2min_best_amhg_val <- bam_val(Tan_3m_2min_best_amhg, Tan_qobs_2min)
Tan_3m_2min_best_man_val <- bam_val(Tan_3m_2min_best_man, Tan_qobs_2min)
Tan_3m_2min_best_man_amhg_val <- bam_val(Tan_3m_2min_best_man_amhg, Tan_qobs_2min)

#Plot Hydrograph
par(family = "serif")
plot(Tan_qobs_2min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "Tanana 2min best NED",
     ylim = c(135, 155))
lines(Tan_NED_2min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(Tan_NED_2min_best_man_val[[1]], col = "blue", lwd = 2)
lines(Tan_NED_2min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(Tan_3m_2min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(Tan_3m_2min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(Tan_3m_2min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
Tan_daily_lower_logQ <- log(seq(exp(max(apply(log(Tan_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                min(Tan_qobs_2min),length.out = 10))
Tan_daily_lower_logQ_priors <- list()
for(i in 1:10){
  Tan_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_logQ = Tan_daily_lower_logQ[i])
}
#upperbound_logQ
Tan_daily_upper_logQ <- log(seq(max(Tan_qobs_2min), exp(min(apply(log(Tan_width_NED_daily), 2, max)) + log(40) + log(5)),
                                length.out = 10))
Tan_daily_upper_logQ_priors <- list()
for(i in 1:10){
  Tan_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_logQ = Tan_daily_upper_logQ[i])
}
#lowerbound_A0
Tan_daily_lower_A0 <- seq(0, 30, length.out = 10)
Tan_daily_lower_A0_priors <- list()
for(i in 1:10){
  Tan_daily_lower_A0_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_A0 = Tan_daily_lower_A0[i])
}
#upperbound_A0
Tan_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
Tan_daily_upper_A0_priors <- list()
for(i in 1:10){
  Tan_daily_upper_A0_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_A0 = Tan_daily_upper_A0[i])
}
#lowerbound_logn
Tan_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
Tan_daily_lower_logn_priors <- list()
for(i in 1:10){
  Tan_daily_lower_logn_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_logn = Tan_daily_lower_logn[i])
}
#upperbound_logn
Tan_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
Tan_daily_upper_logn_priors <- list()
for(i in 1:10){
  Tan_daily_upper_logn_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_logn = Tan_daily_upper_logn[i])
}
#lowerbound_logQc
Tan_daily_lower_logQc <- seq(0, log(min(Tan_qobs_2min)),length.out = 10)
Tan_daily_lower_logQc_priors <- list()
for(i in 1:10){
  Tan_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_logQc = Tan_daily_lower_logQc[i])
}
#upperbound_logQc
Tan_daily_upper_logQc <- seq(log(max(Tan_qobs_2min)), 10, length.out = 10)
Tan_daily_upper_logQc_priors <- list()
for(i in 1:10){
  Tan_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_logQc = Tan_daily_upper_logQc[i])
}
#lowerbound_logWc
Tan_daily_lower_logWc <- seq(1, log(min(Tan_width_NED_daily)),length.out = 10)
Tan_daily_lower_logWc_priors <- list()
for(i in 1:10){
  Tan_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_logWc = Tan_daily_lower_logWc[i])
}
#upperbound_logWc
Tan_daily_upper_logWc <- seq(log(max(Tan_width_NED_daily)), 8, length.out = 10)
Tan_daily_upper_logWc_priors <- list()
for(i in 1:10){
  Tan_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_logWc = Tan_daily_upper_logWc[i])
}
#lowerbound_b
Tan_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
Tan_daily_lower_b_priors <- list()
for(i in 1:10){
  Tan_daily_lower_b_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, lowerbound_b = Tan_daily_lower_b[i])
}
#upperbound_b
Tan_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
Tan_daily_upper_b_priors <- list()
for(i in 1:10){
  Tan_daily_upper_b_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, upperbound_b = Tan_daily_upper_b[i])
}
#sigma_man
Tan_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
Tan_daily_sigma_man_priors <- list()
for(i in 1:10){
  Tan_daily_sigma_man_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, sigma_man = Tan_daily_sigma_man[i])
}
#sigma_amhg
Tan_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
Tan_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  Tan_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, sigma_amhg = Tan_daily_sigma_amhg[i])
}
#logQc_hat
Tan_daily_logQc_hat <- log(seq(min(Tan_qobs_2min), max(Tan_qobs_2min), length.out = 10))
Tan_daily_logQc_hat_priors <- list()
for(i in 1:10){
  Tan_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logQc_hat = Tan_daily_logQc_hat[i])
}
#logWc_hat
Tan_daily_logWc_hat <- log(seq(min(Tan_width_NED_daily), max(Tan_width_NED_daily), length.out = 10))
Tan_daily_logWc_hat_priors <- list()
for(i in 1:10){
  Tan_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logWc_hat = Tan_daily_logWc_hat[i])
}
#b_hat
Tan_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
Tan_daily_b_hat_priors <- list()
for(i in 1:10){
  Tan_daily_b_hat_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, b_hat = Tan_daily_b_hat[i])
}
#logA0_hat
Tan_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
Tan_daily_logA0_hat_priors <- list()
for(i in 1:10){
  Tan_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logA0_hat = Tan_daily_logA0_hat[i])
}
#logn_hat
Tan_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
Tan_daily_logn_hat_priors <- list()
for(i in 1:10){
  Tan_daily_logn_hat_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logn_hat = Tan_daily_logn_hat[i])
}
#logQ_sd
Tan_daily_logQ_sd <- seq(sd(log(Tan_qobs_2min)), 0.8325546, length.out = 10)
Tan_daily_logQ_sd_priors <- list()
for(i in 1:10){
  Tan_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logQ_sd = Tan_daily_logQ_sd[i])
}
#logQc_sd
Tan_daily_logQc_sd <- seq(sd(log(Tan_qobs_2min)), 0.8325546, length.out = 10)
Tan_daily_logQc_sd_priors <- list()
for(i in 1:10){
  Tan_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logQc_sd = Tan_daily_logQc_sd[i])
}
#logWc_Sd
Tan_daily_logWc_Sd <- seq(log(sd(Tan_width_NED_daily)), 4.712493, length.out = 10)
Tan_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  Tan_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logWc_sd = Tan_daily_logWc_Sd[i])
}
#b_sd
Tan_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
Tan_daily_b_sd_priors <- list()
for(i in 1:10){
  Tan_daily_b_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, b_sd = Tan_daily_b_sd[i])
}
#logA0_sd
Tan_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
Tan_daily_logA0_sd_priors <- list()
for(i in 1:10){
  Tan_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logA0_sd = Tan_daily_logA0_sd[i])
}
#logn_sd
Tan_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
Tan_daily_logn_sd_priors <- list()
for(i in 1:10){
  Tan_daily_logn_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, logn_sd = Tan_daily_logn_sd[i])
}
#Werr_sd
Tan_daily_Werr_sd <- seq(5, 15, length.out = 10)
Tan_daily_Werr_sd_priors <- list()
for(i in 1:10){
  Tan_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, Werr_sd = Tan_daily_Werr_sd[i])
}
#Serr_sd
Tan_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
Tan_daily_Serr_sd_priors <- list()
for(i in 1:10){
  Tan_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, Serr_sd = Tan_daily_Serr_sd[i])
}
#dAerr_sd
Tan_daily_dAerr_sd <- seq(5, 15, length.out = 10)
Tan_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  Tan_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= Tan_data_NED_daily, dAerr_sd = Tan_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
Tan_daily_manning_list <- list(Tan_daily_lower_logQ_priors,
                               Tan_daily_upper_logQ_priors,
                               Tan_daily_lower_A0_priors,
                               Tan_daily_upper_A0_priors,
                               Tan_daily_lower_logn_priors,
                               Tan_daily_upper_logn_priors,
                               Tan_daily_sigma_man_priors,
                               Tan_daily_logA0_hat_priors,
                               Tan_daily_logn_hat_priors,
                               Tan_daily_logQ_sd_priors,
                               Tan_daily_logA0_sd_priors,
                               Tan_daily_logn_sd_priors,
                               Tan_daily_Werr_sd_priors,
                               Tan_daily_Serr_sd_priors,
                               Tan_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/Tanana/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Tanana/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = Tan_data_NED_daily, variant = "manning", bampriors = Tan_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
Tan_daily_amhg_list <- list(Tan_daily_lower_logQ_priors,
                            Tan_daily_upper_logQ_priors,
                            Tan_daily_lower_logQc_priors,
                            Tan_daily_upper_logQc_priors,
                            Tan_daily_lower_logWc_priors,
                            Tan_daily_upper_logWc_priors,
                            Tan_daily_lower_b_priors,
                            Tan_daily_upper_b_priors,
                            Tan_daily_sigma_amhg_priors,
                            Tan_daily_logQc_hat_priors,
                            Tan_daily_logWc_hat_priors,
                            Tan_daily_b_hat_priors,
                            Tan_daily_logQ_sd_priors,
                            Tan_daily_logQc_sd_priors,
                            Tan_daily_logWc_Sd_priors,
                            Tan_daily_b_sd_priors,
                            Tan_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/Tanana/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/Tanana/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = Tan_data_NED_daily, variant = "amhg", bampriors = Tan_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
Tan_daily_man_amhg_list <- list(Tan_daily_lower_logQ_priors,
                                Tan_daily_upper_logQ_priors,
                                Tan_daily_lower_A0_priors,
                                Tan_daily_upper_A0_priors,
                                Tan_daily_lower_logn_priors,
                                Tan_daily_upper_logn_priors,
                                Tan_daily_lower_logQc_priors,
                                Tan_daily_upper_logQc_priors,
                                Tan_daily_lower_logWc_priors,
                                Tan_daily_upper_logWc_priors,
                                Tan_daily_lower_b_priors,
                                Tan_daily_upper_b_priors,
                                Tan_daily_sigma_man_priors,
                                Tan_daily_sigma_amhg_priors,
                                Tan_daily_logQc_hat_priors,
                                Tan_daily_logWc_hat_priors,
                                Tan_daily_b_hat_priors,
                                Tan_daily_logA0_hat_priors,
                                Tan_daily_logn_hat_priors,
                                Tan_daily_logQ_sd_priors,
                                Tan_daily_logQc_sd_priors,
                                Tan_daily_logWc_Sd_priors,
                                Tan_daily_b_sd_priors,
                                Tan_daily_logA0_sd_priors,
                                Tan_daily_logn_sd_priors,
                                Tan_daily_Werr_sd_priors,
                                Tan_daily_Serr_sd_priors,
                                Tan_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/Tanana/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/Tanana/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = Tan_data_NED_daily, variant = "manning_amhg", bampriors = Tan_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

###North_Sask
#BAM data-----
home_dir <- "C:/Users/Merritt/Desktop/Research/PT_Paper"

setwd(paste0(home_dir, "/North_Saskatchewan_2017"))
load("Sask_WSE_daily.RData")
load("Sask_WSE_6hr.RData")
load("Sask_WSE_hourly.RData")
load("Sask_WSE_15min.RData")

North_Sask_PT_coords <- read.csv("North_Sask_PT_coords.csv")
North_Sask_PT_coords_sp <- SpatialPoints(North_Sask_PT_coords[,c(4,3)], proj4string =  CRS("+init=epsg:4326"))
North_Sask_PT_coords_sp <-spTransform(North_Sask_PT_coords_sp, CRS("+proj=utm +zone=13 +datum=WGS84"))
North_Sask_cl <-readOGR(dsn = ".", layer = "North_Sask_cl_UTM")

setwd(paste0(home_dir, "/North_Saskatchewan_2017/DEM"))
North_Sask_DEM_NED <- raster("North_Sask_DEM.tif")
North_Sask_DEM_3m <- raster()

setwd(paste0(home_dir, "/North_Saskatchewan_2017"))
North_Sask_orthos_NED <- ortho_lines(North_Sask_DEM_NED, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords), 500)
North_Sask_dist_mat_NED <- pointDistance(North_Sask_PT_coords_sp, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords))
North_Sask_closest_ortho_NED <- North_Sask_orthos_NED[apply(North_Sask_dist_mat_NED, 1, FUN = which.min)]

North_Sask_orthos_3m <- ortho_lines(North_Sask_DEM_3m, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords), 500)
North_Sask_dist_mat_3m <- pointDisNorth_Saskce(North_Sask_PT_coords_sp, data.frame(North_Sask_cl@lines[[1]]@Lines[[1]]@coords))
North_Sask_closest_ortho_3m <- North_Sask_orthos_3m[apply(North_Sask_dist_mat_3m, 1, FUN = which.min)]

numcores <- detectCores()
cl<-makeCluster(numcores-1)
registerDoParallel(cl)
ptm<-proc.time()
North_Sask_elev_prof_NED <- foreach(i=1:length(North_Sask_closest_ortho_NED), 
                                    .combine='rbind', .packages=c('raster')) %dopar% {
                                      extract(North_Sask_DEM_NED, North_Sask_closest_ortho_NED[[i]], method = 'simple')
                                    }
North_Sask_elev_prof_3m <- foreach(i=1:length(North_Sask_closest_ortho_3m), 
                                   .combine='rbind', .packages=c('raster')) %dopar% {
                                     extract(North_Sask_DEM_3m, North_Sask_closest_ortho_3m[[i]], method = 'simple')
                                   }
proc.time()-ptm
stopCluster(cl)

North_Sask_scale_x_NED <- vector()
for(i in 1:length(North_Sask_elev_prof_NED)){
  North_Sask_scale_x_NED[i] = 1000/length(North_Sask_elev_prof_NED[[i]])
}

North_Sask_scale_x_3m <- vector()
for(i in 1:length(North_Sask_elev_prof_3m)){
  North_Sask_scale_x_3m[i] = 1000/length(North_Sask_elev_prof_3m[[i]])
}

par(family = "serif")
plot_elev_prof(North_Sask_elev_prof_NED, Sask_WSE_daily_df, "North_Sask")
locator(n = 4)
North_Sask_upper_l_NED <- c(22.6,  7.2, 19.7, 29.4, 23.0, 23.5, 29.4,  1.7, 2.2)
North_Sask_lower_l_NED <- c(29.3, 10.0, 28.3, 35.3, 25.5, 28.9, 33.7, 23.4, 30.6)
North_Sask_lower_r_NED <- c(46.5, 43.6, 71.6, 50.2, 39.9, 44.4, 43.4, 55.2, 49.0)
North_Sask_upper_r_NED <- c(51.9, 49.4, 76.2, 53.7, 43.9, 53.8, 50.0, 65.6, 61.7)

North_Sask_elev_prof_NED_corr <- North_Sask_elev_prof_NED
#North_Sask_elev_prof_NED_corr[[12]] <- North_Sask_elev_prof_NED[[13]]
#North_Sask_scale_x_NED[12] <- North_Sask_scale_x_NED[13]

North_Sask_width_NED_daily <- calc_width(h = t(Sask_WSE_daily_df), elev_prof = North_Sask_elev_prof_NED_corr, 
                                         upper_l = North_Sask_upper_l_NED, lower_l = North_Sask_lower_l_NED,
                                         lower_r = North_Sask_lower_r_NED, upper_r = North_Sask_upper_r_NED,
                                         scale_prof = North_Sask_scale_x_NED)
North_Sask_width_NED_daily <- na.approx(North_Sask_width_NED_daily)
North_Sask_dA_NED_daily <- calcdA_mat(w = t(North_Sask_width_NED_daily), h = Sask_WSE_daily_df)

North_Sask_width_NED_6hr <- calc_width(h = t(Sask_WSE_6hr_df), elev_prof = North_Sask_elev_prof_NED_corr, 
                                       upper_l = North_Sask_upper_l_NED, lower_l = North_Sask_lower_l_NED,
                                       lower_r = North_Sask_lower_r_NED, upper_r = North_Sask_upper_r_NED,
                                       scale_prof = North_Sask_scale_x_NED)
North_Sask_width_NED_6hr <- na.approx(North_Sask_width_NED_6hr)
North_Sask_dA_NED_6hr <- calcdA_mat(w = t(North_Sask_width_NED_6hr), h = Sask_WSE_6hr_df)

North_Sask_width_NED_hourly <- calc_width(h = t(Sask_WSE_hourly_df), elev_prof = North_Sask_elev_prof_NED_corr, 
                                          upper_l = North_Sask_upper_l_NED, lower_l = North_Sask_lower_l_NED,
                                          lower_r = North_Sask_lower_r_NED, upper_r = North_Sask_upper_r_NED,
                                          scale_prof = North_Sask_scale_x_NED)
North_Sask_width_NED_hourly <- na.approx(North_Sask_width_NED_hourly)
North_Sask_dA_NED_hourly <- calcdA_mat(w = t(North_Sask_width_NED_hourly), h = Sask_WSE_hourly_df)

North_Sask_width_NED_2min <- calc_width(h = t(Sask_WSE_df), elev_prof = North_Sask_elev_prof_NED_corr, 
                                        upper_l = North_Sask_upper_l_NED, lower_l = North_Sask_lower_l_NED,
                                        lower_r = North_Sask_lower_r_NED, upper_r = North_Sask_upper_r_NED,
                                        scale_prof = North_Sask_scale_x_NED)
North_Sask_width_NED_2min <- na.approx(North_Sask_width_NED_2min)
North_Sask_dA_NED_2min <- calcdA_mat(w = t(North_Sask_width_NED_2min), h = Sask_WSE_df)

plot_elev_prof(North_Sask_elev_prof_3m, Sask_WSE_daily_df, "North_Sask Reach Three")
locator(n = 4)
North_Sask_upper_l_3m <- c(157.7, 169.0, 117.5, 201.2, 143.3, 129.0, 143.9,  33.0, 132.9, 230.4, 182.2, 211.3,  88.1, 155.0,  63.5,  90.4, 197.4, 203.8, 106.4)
North_Sask_lower_l_3m <- c(175.8, 176.1, 123.6, 216.4, 156.9, 136.9, 150.7,  54.5, 159.9, 237.0, 191.0, 226.7,  94.7, 168.1,  72.7, 103.3, 204.3, 209.4, 130.2)
North_Sask_lower_r_3m <- c(209.0, 218.0, 216.5, 303.7, 353.8, 343.5, 209.0, 455.0, 390.5, 376.9, 262.4, 309.3, 178.2, 275.2, 179.1, 160.2, 253.8, 242.9, 262.4)
North_Sask_upper_r_3m <- c(231.2, 226.1, 225.0, 321.9, 371.9, 356.9, 225.0, 468.1, 401.4, 388.1, 283.4, 323.1, 184.8, 285.0, 207.6, 175.4, 266.2, 252.2, 275.2)

North_Sask_width_3m_daily <- calc_width(h = t(Sask_WSE_daily_df), elev_prof = North_Sask_elev_prof_3m, 
                                        upper_l = North_Sask_upper_l_3m, lower_l = North_Sask_lower_l_3m,
                                        lower_r = North_Sask_lower_r_3m, upper_r = North_Sask_upper_r_3m,
                                        scale_prof = North_Sask_scale_x_3m)
North_Sask_width_3m_daily <- na.approx(North_Sask_width_3m_daily)
North_Sask_dA_3m_daily <- calcdA_mat(w = t(North_Sask_width_3m_daily), h = Sask_WSE_daily_df)

North_Sask_width_3m_6hr <- calc_width(h = t(Sask_WSE_6hr_df), elev_prof = North_Sask_elev_prof_3m, 
                                      upper_l = North_Sask_upper_l_3m, lower_l = North_Sask_lower_l_3m,
                                      lower_r = North_Sask_lower_r_3m, upper_r = North_Sask_upper_r_3m,
                                      scale_prof = North_Sask_scale_x_3m)
North_Sask_width_3m_6hr <- na.approx(North_Sask_width_3m_6hr)
North_Sask_dA_3m_6hr <- calcdA_mat(w = t(North_Sask_width_3m_6hr), h = Sask_WSE_6hr_df)

North_Sask_width_3m_hourly <- calc_width(h = t(Sask_WSE_hourly_df), elev_prof = North_Sask_elev_prof_3m, 
                                         upper_l = North_Sask_upper_l_3m, lower_l = North_Sask_lower_l_3m,
                                         lower_r = North_Sask_lower_r_3m, upper_r = North_Sask_upper_r_3m,
                                         scale_prof = North_Sask_scale_x_3m)
North_Sask_width_3m_hourly <- na.approx(North_Sask_width_3m_hourly)
North_Sask_dA_3m_hourly <- calcdA_mat(w = t(North_Sask_width_3m_hourly), h = Sask_WSE_hourly_df)

North_Sask_width_3m_2min <- calc_width(h = t(Sask_WSE_df), elev_prof = North_Sask_elev_prof_3m, 
                                       upper_l = North_Sask_upper_l_3m, lower_l = North_Sask_lower_l_3m,
                                       lower_r = North_Sask_lower_r_3m, upper_r = North_Sask_upper_r_3m,
                                       scale_prof = North_Sask_scale_x_3m)
North_Sask_width_3m_2min <- na.approx(North_Sask_width_3m_2min)
North_Sask_dA_3m_2min <- calcdA_mat(w = t(North_Sask_width_3m_2min), h = Sask_WSE_df)


xvec <- calc_xvec(North_Sask_PT_coords_sp, North_Sask_cl)

North_Sask_xvec <- saveRDS(xvec, file = "North_Sask_xvec.rds")
North_Sask_slope_daily <- calcslope((xvec), hmat = Sask_WSE_daily_df)
North_Sask_slope_daily <- na.replace(North_Sask_slope_daily, mean(North_Sask_slope_daily))
North_Sask_slope_daily[1,] <- colMeans(North_Sask_slope_daily, na.rm = TRUE)
North_Sask_slope_6hr <- calcslope((xvec), hmat = Sask_WSE_6hr_df)
North_Sask_slope_6hr <- na.approx(North_Sask_slope_6hr)
North_Sask_slope_6hr[1,] <- colMeans(North_Sask_slope_6hr, na.rm = TRUE)
North_Sask_slope_hourly <- calcslope((xvec), hmat = Sask_WSE_hourly_df)
North_Sask_slope_hourly <- na.approx(North_Sask_slope_hourly)
North_Sask_slope_hourly[1,] <- colMeans(North_Sask_slope_hourly, na.rm = TRUE)
North_Sask_slope_2min <- calcslope((xvec), hmat = Sask_WSE_df)
North_Sask_slope_2min <- na.approx(North_Sask_slope_2min)
North_Sask_slope_2min[1,] <- colMeans(North_Sask_slope_2min, na.rm = TRUE)

setwd("C:/Users/Merritt/Desktop/Research/PT_Paper/North_Saskatchewan_2017/Qobs")
North_Saskqobs <- read.csv("Sask_daily_Qobs.csv", header = TRUE)$Value
North_Sask_qobs_daily<- North_Saskqobs
North_Sask_qobs_6hr<- colMeans(matrix(North_Saskqobs[1:3048], 24))
North_Sask_qobs_hourly<- colMeans(matrix(North_Saskqobs[1:3056], 4))

save(North_Sask_width_NED_daily, file = "North_Sask_width_NED_daily.RData")
save(North_Sask_dA_NED_daily, file = "North_Sask_dA_NED_daily.RData")
save(North_Sask_slope_daily, file = "North_Sask_slope_daily.RData")
save(North_Sask_qobs_daily, file = "North_Sask_qobs_daily.RData")
save(North_Sask_width_NED_6hr, file = "North_Sask_width_NED_6hr.RData")
save(North_Sask_dA_NED_6hr, file = "North_Sask_dA_NED_6hr.RData")
save(North_Sask_slope_6hr, file = "North_Sask_slope_6hr.RData")
save(North_Sask_qobs_6hr, file = "North_Sask_qobs_6hr.RData")
save(North_Sask_width_NED_hourly, file = "North_Sask_width_NED_hourly.RData")
save(North_Sask_dA_NED_hourly, file = "North_Sask_dA_NED_hourly.RData")
save(North_Sask_slope_hourly, file = "North_Sask_slope_hourly.RData")
save(North_Sask_qobs_hourly, file = "North_Sask_qobs_hourly.RData")
save(North_Sask_width_NED_15min, file = "North_Sask_width_NED_15min.RData")
save(North_Sask_dA_NED_15min, file = "North_Sask_dA_NED_15min.RData")
save(North_Sask_slope_15min, file = "North_Sask_slope_15min.RData")
save(North_Sask_qobs_15min, file = "North_Sask_qobs_15min.RData")
#North_Sask_qobs_2min<- colMeans(matrix(North_Saskqobs[1:3701], 1))

#Plot Data----
par(mfrow=c(3,2), mai = c(0.6, 0.75, 0.5, 0.25))
colfunc <- colorRampPalette(c("cyan", "darkblue"))(10)
par(family = "serif")
Sask_WSE_plot <- matplot(t(Sask_WSE_daily_df), type = c("l"), col= colfunc,
                         main = "WSE", xlab = "Days",
                         ylab = "WSE (m)",
                         cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
North_Sask_slope_plot <- matplot(t(North_Sask_slope_daily), type = c("l"), col= colfunc,
                                 main = "Slope", xlab = "Days",
                                 ylab = "Slope",
                                 cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
North_Sask_width_plot <- matplot(North_Sask_width_NED_daily, type = c("l"), col= colfunc,
                                 main = "Width (NED)", xlab = "Days",
                                 ylab = "Width (m)",
                                 cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
North_Sask_width_plot <- matplot(North_Sask_width_3m_daily, type = c("l"), col= colfunc,
                                 main = "Width (Lidar)", xlab = "Days",
                                 ylab = "Width (m)",
                                 cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
North_Sask_dA_plot <- matplot(t(North_Sask_dA_NED_daily), type = c("l"), col= colfunc,
                              main = "dA (NED)", xlab = "Days",
                              ylab = "dA (m^2)",
                              cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
North_Sask_dA_plot <- matplot(t(North_Sask_dA_3m_daily), type = c("l"), col= colfunc,
                              main = "dA (Lidar)", xlab = "Days",
                              ylab = "dA (m^2)",
                              cex.lab = 1.5, cex.axis = 1.5, cex.main =2) #plot
#Run BAM----
#Default-----
North_Sask_data_NED_daily <- bam_data(w = t(North_Sask_width_NED_daily), s = North_Sask_slope_daily, dA = North_Sask_dA_NED_daily,
                                      Qhat = na.approx(North_Sask_qobs_daily[c(2:43)]))
North_Sask_NED_daily_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "manning_amhg")
North_Sask_NED_daily_default_amhg <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "amhg")
North_Sask_NED_daily_default_man <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "manning")
library(hydroGOF)
North_Sask_NED_daily_default_amhg_val <- bam_val(North_Sask_NED_daily_default_amhg, na.approx(North_Sask_qobs_daily[2:43]))
North_Sask_NED_daily_default_man_val <- bam_val(North_Sask_NED_daily_default_man, na.approx(North_Sask_qobs_daily[2:43]))
North_Sask_NED_daily_default_man_amhg_val <- bam_val(North_Sask_NED_daily_default_man_amhg, na.approx(North_Sask_qobs_daily[2:43]))

North_Sask_data_3m_daily <- bam_data(w = t(North_Sask_width_3m_daily), s = North_Sask_slope_daily, dA = North_Sask_dA_3m_daily,
                                     Qhat = North_Sask_qobs_daily)
North_Sask_3m_daily_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "manning_amhg")
North_Sask_3m_daily_default_amhg <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "amhg")
North_Sask_3m_daily_default_man <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "manning")

North_Sask_3m_daily_default_amhg_val <- bam_val(North_Sask_3m_daily_default_amhg, North_Sask_qobs_daily)
North_Sask_3m_daily_default_man_val <- bam_val(North_Sask_3m_daily_default_man, North_Sask_qobs_daily)
North_Sask_3m_daily_default_man_amhg_val <- bam_val(North_Sask_3m_daily_default_man_amhg, North_Sask_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask Daily Default NED",
     ylim = c(0, 300))
lines(North_Sask_NED_daily_default_amhg_val[[1]], col = "red")
lines(North_Sask_NED_daily_default_man_val[[1]], col = "blue")
lines(North_Sask_NED_daily_default_man_amhg_val[[1]], col = "purple")
lines(North_Sask_3m_daily_default_amhg_val[[1]], col = "red", lty = 2)
lines(North_Sask_3m_daily_default_man_val[[1]], col = "blue", lty = 2)
lines(North_Sask_3m_daily_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
North_Sask_data_NED_daily <- bam_data(w = t(North_Sask_width_NED_daily), s = North_Sask_slope_daily, dA = North_Sask_dA_NED_daily,
                                      Qhat = na.approx(North_Sask_qobs_daily[c(2:43)]))
North_Sask_best_priors <- bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_logQ = log(0.5 * min(na.approx(North_Sask_qobs_daily[c(2:43)]))),
                                     upperbound_logQ = log(1.1 * max(na.approx(North_Sask_qobs_daily[c(2:43)]))),
                                     lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(North_Sask_slope_daily, na.rm = TRUE),
                                     lowerbound_A0 = 0, 
                                     upperbound_A0 = 1e+03)
North_Sask_NED_daily_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "manning_amhg", 
                                                   bampriors = North_Sask_best_priors)
North_Sask_NED_daily_best_amhg <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_NED_daily_best_man <- bam_estimate(bamdata = North_Sask_data_NED_daily,variant = "manning", bampriors = North_Sask_best_priors)
library(hydroGOF)
North_Sask_NED_daily_best_amhg_val <- bam_val(North_Sask_NED_daily_best_amhg, na.approx(North_Sask_qobs_daily[c(2:43)]))
North_Sask_NED_daily_best_man_val <- bam_val(North_Sask_NED_daily_best_man, na.approx(North_Sask_qobs_daily[c(2:43)]))
North_Sask_NED_daily_best_man_amhg_val <- bam_val(North_Sask_NED_daily_best_man_amhg, na.approx(North_Sask_qobs_daily[c(2:43)]))

North_Sask_data_3m_daily <- bam_data(w = t(North_Sask_width_3m_daily), s = North_Sask_slope_daily, dA = North_Sask_dA_3m_daily,
                                     Qhat = North_Sask_qobs_daily)
North_Sask_3m_daily_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "manning_amhg",
                                                  bampriors = North_Sask_best_priors)
North_Sask_3m_daily_best_amhg <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_3m_daily_best_man <- bam_estimate(bamdata = North_Sask_data_3m_daily,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_3m_daily_best_amhg_val <- bam_val(North_Sask_3m_daily_best_amhg, North_Sask_qobs_daily)
North_Sask_3m_daily_best_man_val <- bam_val(North_Sask_3m_daily_best_man, North_Sask_qobs_daily)
North_Sask_3m_daily_best_man_amhg_val <- bam_val(North_Sask_3m_daily_best_man_amhg, North_Sask_qobs_daily)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_daily, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask Daily best NED",
     ylim = c(0, 20))
lines(North_Sask_NED_daily_best_amhg_val[[1]], col = "red", lwd = 2)
lines(North_Sask_NED_daily_best_man_val[[1]], col = "blue", lwd = 2)
lines(North_Sask_NED_daily_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(North_Sask_3m_daily_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(North_Sask_3m_daily_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(North_Sask_3m_daily_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#####Sensitivity Analysis----
#Temporal Resolution
#6hr
#Default----
North_Sask_data_NED_6hr <- bam_data(w = t(North_Sask_width_NED_6hr), s = North_Sask_slope_6hr, dA = North_Sask_dA_NED_6hr,
                                    Qhat = North_Sask_qobs_6hr)
North_Sask_NED_6hr_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "manning_amhg")
North_Sask_NED_6hr_default_amhg <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "amhg")
North_Sask_NED_6hr_default_man <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "manning")

North_Sask_NED_6hr_default_amhg_val <- bam_val(North_Sask_NED_6hr_default_amhg, North_Sask_qobs_6hr)
North_Sask_NED_6hr_default_man_val <- bam_val(North_Sask_NED_6hr_default_man, North_Sask_qobs_6hr)
North_Sask_NED_6hr_default_man_amhg_val <- bam_val(North_Sask_NED_6hr_default_man_amhg, North_Sask_qobs_6hr)

North_Sask_data_3m_6hr <- bam_data(w = t(North_Sask_width_3m_6hr), s = North_Sask_slope_6hr, dA = North_Sask_dA_3m_6hr, 
                                   Qhat = North_Sask_qobs_6hr)
North_Sask_3m_6hr_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "manning_amhg")
North_Sask_3m_6hr_default_amhg <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "amhg")
North_Sask_3m_6hr_default_man <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "manning")

North_Sask_3m_6hr_default_amhg_val <- bam_val(North_Sask_3m_6hr_default_amhg, North_Sask_qobs_6hr)
North_Sask_3m_6hr_default_man_val <- bam_val(North_Sask_3m_6hr_default_man, North_Sask_qobs_6hr)
North_Sask_3m_6hr_default_man_amhg_val <- bam_val(North_Sask_3m_6hr_default_man_amhg, North_Sask_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask 6hr Default NED",
     ylim = c(120, 275))
lines(North_Sask_NED_6hr_default_amhg_val[[1]], col = "red")
lines(North_Sask_NED_6hr_default_man_val[[1]], col = "blue")
lines(North_Sask_NED_6hr_default_man_amhg_val[[1]], col = "purple")
lines(North_Sask_3m_6hr_default_amhg_val[[1]], col = "red", lty = 2)
lines(North_Sask_3m_6hr_default_man_val[[1]], col = "blue", lty = 2)
lines(North_Sask_3m_6hr_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
North_Sask_data_NED_6hr <- bam_data(w = t(North_Sask_width_NED_6hr), s = North_Sask_slope_6hr, dA = North_Sask_dA_NED_6hr,
                                    Qhat = North_Sask_qobs_6hr)
North_Sask_best_priors <- bam_priors(bamdata= North_Sask_data_NED_6hr, lowerbound_logQ = log(0.75 * min(North_Sask_qobs_6hr)),
                                     upperbound_logQ = log(1.25 * max(North_Sask_qobs_6hr)),
                                     lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(North_Sask_slope),
                                     lowerbound_A0 = 0, 
                                     upperbound_A0 = 1e+05)
North_Sask_NED_6hr_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "manning_amhg", 
                                                 bampriors = North_Sask_best_priors)
North_Sask_NED_6hr_best_amhg <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_NED_6hr_best_man <- bam_estimate(bamdata = North_Sask_data_NED_6hr,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_NED_6hr_best_amhg_val <- bam_val(North_Sask_NED_6hr_best_amhg, North_Sask_qobs_6hr)
North_Sask_NED_6hr_best_man_val <- bam_val(North_Sask_NED_6hr_best_man, North_Sask_qobs_6hr)
North_Sask_NED_6hr_best_man_amhg_val <- bam_val(North_Sask_NED_6hr_best_man_amhg, North_Sask_qobs_6hr)

North_Sask_data_3m_6hr <- bam_data(w = t(North_Sask_width_3m_6hr), s = North_Sask_slope_6hr, dA = North_Sask_dA_3m_6hr,
                                   Qhat = North_Sask_qobs_6hr)
North_Sask_3m_6hr_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "manning_amhg",
                                                bampriors = North_Sask_best_priors)
North_Sask_3m_6hr_best_amhg <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_3m_6hr_best_man <- bam_estimate(bamdata = North_Sask_data_3m_6hr,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_3m_6hr_best_amhg_val <- bam_val(North_Sask_3m_6hr_best_amhg, North_Sask_qobs_6hr)
North_Sask_3m_6hr_best_man_val <- bam_val(North_Sask_3m_6hr_best_man, North_Sask_qobs_6hr)
North_Sask_3m_6hr_best_man_amhg_val <- bam_val(North_Sask_3m_6hr_best_man_amhg, North_Sask_qobs_6hr)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_6hr, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask 6hr best NED",
     ylim = c(135, 155))
lines(North_Sask_NED_6hr_best_amhg_val[[1]], col = "red", lwd = 2)
lines(North_Sask_NED_6hr_best_man_val[[1]], col = "blue", lwd = 2)
lines(North_Sask_NED_6hr_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(North_Sask_3m_6hr_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(North_Sask_3m_6hr_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(North_Sask_3m_6hr_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#hourly-----
#Default----
North_Sask_data_NED_hourly <- bam_data(w = t(North_Sask_width_NED_hourly), s = North_Sask_slope_hourly, dA = North_Sask_dA_NED_hourly,
                                       Qhat = North_Sask_qobs_hourly)
North_Sask_NED_hourly_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "manning_amhg")
North_Sask_NED_hourly_default_amhg <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "amhg")
North_Sask_NED_hourly_default_man <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "manning")

North_Sask_NED_hourly_default_amhg_val <- bam_val(North_Sask_NED_hourly_default_amhg, North_Sask_qobs_hourly)
North_Sask_NED_hourly_default_man_val <- bam_val(North_Sask_NED_hourly_default_man, North_Sask_qobs_hourly)
North_Sask_NED_hourly_default_man_amhg_val <- bam_val(North_Sask_NED_hourly_default_man_amhg, North_Sask_qobs_hourly)

North_Sask_data_3m_hourly <- bam_data(w = t(North_Sask_width_3m_hourly), s = North_Sask_slope_hourly, dA = North_Sask_dA_3m_hourly, 
                                      Qhat = North_Sask_qobs_hourly)
North_Sask_3m_hourly_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "manning_amhg")
North_Sask_3m_hourly_default_amhg <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "amhg")
North_Sask_3m_hourly_default_man <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "manning")

North_Sask_3m_hourly_default_amhg_val <- bam_val(North_Sask_3m_hourly_default_amhg, North_Sask_qobs_hourly)
North_Sask_3m_hourly_default_man_val <- bam_val(North_Sask_3m_hourly_default_man, North_Sask_qobs_hourly)
North_Sask_3m_hourly_default_man_amhg_val <- bam_val(North_Sask_3m_hourly_default_man_amhg, North_Sask_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask hourly Default NED",
     ylim = c(120, 275))
lines(North_Sask_NED_hourly_default_amhg_val[[1]], col = "red")
lines(North_Sask_NED_hourly_default_man_val[[1]], col = "blue")
lines(North_Sask_NED_hourly_default_man_amhg_val[[1]], col = "purple")
lines(North_Sask_3m_hourly_default_amhg_val[[1]], col = "red", lty = 2)
lines(North_Sask_3m_hourly_default_man_val[[1]], col = "blue", lty = 2)
lines(North_Sask_3m_hourly_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
North_Sask_data_NED_hourly <- bam_data(w = t(North_Sask_width_NED_hourly), s = North_Sask_slope_hourly, dA = North_Sask_dA_NED_hourly,
                                       Qhat = North_Sask_qobs_hourly)
North_Sask_best_priors <- bam_priors(bamdata= North_Sask_data_NED_hourly, lowerbound_logQ = log(0.75 * min(North_Sask_qobs_hourly)),
                                     upperbound_logQ = log(1.25 * max(North_Sask_qobs_hourly)),
                                     lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(North_Sask_slope),
                                     lowerbound_A0 = 0, 
                                     upperbound_A0 = 1e+05)
North_Sask_NED_hourly_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "manning_amhg", 
                                                    bampriors = North_Sask_best_priors)
North_Sask_NED_hourly_best_amhg <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_NED_hourly_best_man <- bam_estimate(bamdata = North_Sask_data_NED_hourly,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_NED_hourly_best_amhg_val <- bam_val(North_Sask_NED_hourly_best_amhg, North_Sask_qobs_hourly)
North_Sask_NED_hourly_best_man_val <- bam_val(North_Sask_NED_hourly_best_man, North_Sask_qobs_hourly)
North_Sask_NED_hourly_best_man_amhg_val <- bam_val(North_Sask_NED_hourly_best_man_amhg, North_Sask_qobs_hourly)

North_Sask_data_3m_hourly <- bam_data(w = t(North_Sask_width_3m_hourly), s = North_Sask_slope_hourly, dA = North_Sask_dA_3m_hourly,
                                      Qhat = North_Sask_qobs_hourly)
North_Sask_3m_hourly_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "manning_amhg",
                                                   bampriors = North_Sask_best_priors)
North_Sask_3m_hourly_best_amhg <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_3m_hourly_best_man <- bam_estimate(bamdata = North_Sask_data_3m_hourly,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_3m_hourly_best_amhg_val <- bam_val(North_Sask_3m_hourly_best_amhg, North_Sask_qobs_hourly)
North_Sask_3m_hourly_best_man_val <- bam_val(North_Sask_3m_hourly_best_man, North_Sask_qobs_hourly)
North_Sask_3m_hourly_best_man_amhg_val <- bam_val(North_Sask_3m_hourly_best_man_amhg, North_Sask_qobs_hourly)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_hourly, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask hourly best NED",
     ylim = c(135, 155))
lines(North_Sask_NED_hourly_best_amhg_val[[1]], col = "red", lwd = 2)
lines(North_Sask_NED_hourly_best_man_val[[1]], col = "blue", lwd = 2)
lines(North_Sask_NED_hourly_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(North_Sask_3m_hourly_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(North_Sask_3m_hourly_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(North_Sask_3m_hourly_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))
#2min-----
#Default----
North_Sask_data_NED_2min <- bam_data(w = t(North_Sask_width_NED_2min), s = North_Sask_slope_2min, dA = North_Sask_dA_NED_2min,
                                     Qhat = North_Sask_qobs_2min)
North_Sask_NED_2min_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "manning_amhg")
North_Sask_NED_2min_default_amhg <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "amhg")
North_Sask_NED_2min_default_man <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "manning")

library(hydroGOF)
North_Sask_NED_2min_default_amhg_val <- bam_val(North_Sask_NED_2min_default_amhg, North_Sask_qobs_2min)
North_Sask_NED_2min_default_man_val <- bam_val(North_Sask_NED_2min_default_man, North_Sask_qobs_2min)
North_Sask_NED_2min_default_man_amhg_val <- bam_val(North_Sask_NED_2min_default_man_amhg, North_Sask_qobs_2min)

North_Sask_data_3m_2min <- bam_data(w = t(North_Sask_width_3m_2min), s = North_Sask_slope_2min, dA = North_Sask_dA_3m_2min, 
                                    Qhat = North_Sask_qobs_2min)
North_Sask_3m_2min_default_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "manning_amhg")
North_Sask_3m_2min_default_amhg <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "amhg")
North_Sask_3m_2min_default_man <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "manning")

library(hydroGOF)
North_Sask_3m_2min_default_amhg_val <- bam_val(North_Sask_3m_2min_default_amhg, North_Sask_qobs_2min)
North_Sask_3m_2min_default_man_val <- bam_val(North_Sask_3m_2min_default_man, North_Sask_qobs_2min)
North_Sask_3m_2min_default_man_amhg_val <- bam_val(North_Sask_3m_2min_default_man_amhg, North_Sask_qobs_2min)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_2min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask 2min Default NED",
     ylim = c(120, 275))
lines(North_Sask_NED_2min_default_amhg_val[[1]], col = "red")
lines(North_Sask_NED_2min_default_man_val[[1]], col = "blue")
lines(North_Sask_NED_2min_default_man_amhg_val[[1]], col = "purple")
lines(North_Sask_3m_2min_default_amhg_val[[1]], col = "red", lty = 2)
lines(North_Sask_3m_2min_default_man_val[[1]], col = "blue", lty = 2)
lines(North_Sask_3m_2min_default_man_amhg_val[[1]], col = "purple", lty = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19))
#Best----
North_Sask_data_NED_2min <- bam_data(w = t(North_Sask_width_NED_2min), s = North_Sask_slope_2min, dA = North_Sask_dA_NED_2min,
                                     Qhat = North_Sask_qobs_2min)
North_Sask_best_priors <- bam_priors(bamdata= North_Sask_data_NED_2min, lowerbound_logQ = log(0.75 * min(North_Sask_qobs_2min)),
                                     upperbound_logQ = log(1.25 * max(North_Sask_qobs_2min)),
                                     lowerbound_logn = log(0.025), upperbound_logn = log(0.8), Serr_sd = 0.1 * mean(North_Sask_slope),
                                     lowerbound_A0 = 0, 
                                     upperbound_A0 = 1e+05)
North_Sask_NED_2min_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "manning_amhg", 
                                                  bampriors = North_Sask_best_priors)
North_Sask_NED_2min_best_amhg <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_NED_2min_best_man <- bam_estimate(bamdata = North_Sask_data_NED_2min,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_NED_2min_best_amhg_val <- bam_val(North_Sask_NED_2min_best_amhg, North_Sask_qobs_2min)
North_Sask_NED_2min_best_man_val <- bam_val(North_Sask_NED_2min_best_man, North_Sask_qobs_2min)
North_Sask_NED_2min_best_man_amhg_val <- bam_val(North_Sask_NED_2min_best_man_amhg, North_Sask_qobs_2min)

North_Sask_data_3m_2min <- bam_data(w = t(North_Sask_width_3m_2min), s = North_Sask_slope_2min, dA = North_Sask_dA_3m_2min,
                                    Qhat = North_Sask_qobs_2min)
North_Sask_3m_2min_best_man_amhg <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "manning_amhg",
                                                 bampriors = North_Sask_best_priors)
North_Sask_3m_2min_best_amhg <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "amhg", bampriors = North_Sask_best_priors)
North_Sask_3m_2min_best_man <- bam_estimate(bamdata = North_Sask_data_3m_2min,variant = "manning", bampriors = North_Sask_best_priors)

North_Sask_3m_2min_best_amhg_val <- bam_val(North_Sask_3m_2min_best_amhg, North_Sask_qobs_2min)
North_Sask_3m_2min_best_man_val <- bam_val(North_Sask_3m_2min_best_man, North_Sask_qobs_2min)
North_Sask_3m_2min_best_man_amhg_val <- bam_val(North_Sask_3m_2min_best_man_amhg, North_Sask_qobs_2min)

#Plot Hydrograph
par(family = "serif")
plot(North_Sask_qobs_2min, pch = 19, col = "black", xlab = "Days", ylab = expression(Discharge~ (m^{3}/s)), 
     main = "North_Sask 2min best NED",
     ylim = c(135, 155))
lines(North_Sask_NED_2min_best_amhg_val[[1]], col = "red", lwd = 2)
lines(North_Sask_NED_2min_best_man_val[[1]], col = "blue", lwd = 2)
lines(North_Sask_NED_2min_best_man_amhg_val[[1]], col = "purple", lwd = 2)
lines(North_Sask_3m_2min_best_amhg_val[[1]], col = "red", lty = 2, lwd = 2)
lines(North_Sask_3m_2min_best_man_val[[1]], col = "blue", lty = 2, lwd = 2)
lines(North_Sask_3m_2min_best_man_amhg_val[[1]], col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("AMHG", "Manning", "Manning-AMHG", "Observed"), col = c("red", "blue", "purple", "black"), 
       lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2,2,2, NA))

#BAM priors-------
#Generate Priors----
#lowerbound_logQ
North_Sask_daily_lower_logQ <- log(seq(exp(max(apply(log(North_Sask_width_NED_daily), 2, min)) + log(0.5) + log(0.5)), 
                                       min(North_Sask_qobs_2min),length.out = 10))
North_Sask_daily_lower_logQ_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_logQ_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_logQ = North_Sask_daily_lower_logQ[i])
}
#upperbound_logQ
North_Sask_daily_upper_logQ <- log(seq(max(North_Sask_qobs_2min), exp(min(apply(log(North_Sask_width_NED_daily), 2, max)) + log(40) + log(5)),
                                       length.out = 10))
North_Sask_daily_upper_logQ_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_logQ_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_logQ = North_Sask_daily_upper_logQ[i])
}
#lowerbound_A0
North_Sask_daily_lower_A0 <- seq(0, 30, length.out = 10)
North_Sask_daily_lower_A0_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_A0_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_A0 = North_Sask_daily_lower_A0[i])
}
#upperbound_A0
North_Sask_daily_upper_A0 <- seq(1e+03, 1e+06, length.out = 10)
North_Sask_daily_upper_A0_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_A0_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_A0 = North_Sask_daily_upper_A0[i])
}
#lowerbound_logn
North_Sask_daily_lower_logn <- seq(-4.6, log(0.025), length.out = 10)
North_Sask_daily_lower_logn_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_logn_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_logn = North_Sask_daily_lower_logn[i])
}
#upperbound_logn
North_Sask_daily_upper_logn <- seq(log(0.8), -1.5, length.out = 10)
North_Sask_daily_upper_logn_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_logn_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_logn = North_Sask_daily_upper_logn[i])
}
#lowerbound_logQc
North_Sask_daily_lower_logQc <- seq(0, log(min(North_Sask_qobs_2min)),length.out = 10)
North_Sask_daily_lower_logQc_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_logQc_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_logQc = North_Sask_daily_lower_logQc[i])
}
#upperbound_logQc
North_Sask_daily_upper_logQc <- seq(log(max(North_Sask_qobs_2min)), 10, length.out = 10)
North_Sask_daily_upper_logQc_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_logQc_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_logQc = North_Sask_daily_upper_logQc[i])
}
#lowerbound_logWc
North_Sask_daily_lower_logWc <- seq(1, log(min(North_Sask_width_NED_daily)),length.out = 10)
North_Sask_daily_lower_logWc_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_logWc_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_logWc = North_Sask_daily_lower_logWc[i])
}
#upperbound_logWc
North_Sask_daily_upper_logWc <- seq(log(max(North_Sask_width_NED_daily)), 8, length.out = 10)
North_Sask_daily_upper_logWc_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_logWc_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_logWc = North_Sask_daily_upper_logWc[i])
}
#lowerbound_b
North_Sask_daily_lower_b <- seq(0.00001, 0.01,length.out = 10)
North_Sask_daily_lower_b_priors <- list()
for(i in 1:10){
  North_Sask_daily_lower_b_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, lowerbound_b = North_Sask_daily_lower_b[i])
}
#upperbound_b
North_Sask_daily_upper_b <- seq(0.5, 0.8, length.out = 10)
North_Sask_daily_upper_b_priors <- list()
for(i in 1:10){
  North_Sask_daily_upper_b_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, upperbound_b = North_Sask_daily_upper_b[i])
}
#sigma_man
North_Sask_daily_sigma_man <- seq(0.2, 0.3, length.out = 10)
North_Sask_daily_sigma_man_priors <- list()
for(i in 1:10){
  North_Sask_daily_sigma_man_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, sigma_man = North_Sask_daily_sigma_man[i])
}
#sigma_amhg
North_Sask_daily_sigma_amhg <- seq(0.2, 0.3, length.out = 10)
North_Sask_daily_sigma_amhg_priors <- list()
for(i in 1:10){
  North_Sask_daily_sigma_amhg_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, sigma_amhg = North_Sask_daily_sigma_amhg[i])
}
#logQc_hat
North_Sask_daily_logQc_hat <- log(seq(min(North_Sask_qobs_2min), max(North_Sask_qobs_2min), length.out = 10))
North_Sask_daily_logQc_hat_priors <- list()
for(i in 1:10){
  North_Sask_daily_logQc_hat_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logQc_hat = North_Sask_daily_logQc_hat[i])
}
#logWc_hat
North_Sask_daily_logWc_hat <- log(seq(min(North_Sask_width_NED_daily), max(North_Sask_width_NED_daily), length.out = 10))
North_Sask_daily_logWc_hat_priors <- list()
for(i in 1:10){
  North_Sask_daily_logWc_hat_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logWc_hat = North_Sask_daily_logWc_hat[i])
}
#b_hat
North_Sask_daily_b_hat <- seq(0.1, 0.3, length.out = 10)
North_Sask_daily_b_hat_priors <- list()
for(i in 1:10){
  North_Sask_daily_b_hat_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, b_hat = North_Sask_daily_b_hat[i])
}
#logA0_hat
North_Sask_daily_logA0_hat <- seq(-0.3, 12, length.out = 10)
North_Sask_daily_logA0_hat_priors <- list()
for(i in 1:10){
  North_Sask_daily_logA0_hat_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logA0_hat = North_Sask_daily_logA0_hat[i])
}
#logn_hat
North_Sask_daily_logn_hat <- seq(log(0.025), log(0.8), length.out = 10)
North_Sask_daily_logn_hat_priors <- list()
for(i in 1:10){
  North_Sask_daily_logn_hat_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logn_hat = North_Sask_daily_logn_hat[i])
}
#logQ_sd
North_Sask_daily_logQ_sd <- seq(sd(log(North_Sask_qobs_2min)), 0.8325546, length.out = 10)
North_Sask_daily_logQ_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_logQ_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logQ_sd = North_Sask_daily_logQ_sd[i])
}
#logQc_sd
North_Sask_daily_logQc_sd <- seq(sd(log(North_Sask_qobs_2min)), 0.8325546, length.out = 10)
North_Sask_daily_logQc_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_logQc_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logQc_sd = North_Sask_daily_logQc_sd[i])
}
#logWc_Sd
North_Sask_daily_logWc_Sd <- seq(log(sd(North_Sask_width_NED_daily)), 4.712493, length.out = 10)
North_Sask_daily_logWc_Sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_logWc_Sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logWc_sd = North_Sask_daily_logWc_Sd[i])
}
#b_sd
North_Sask_daily_b_sd <- seq(0.02, 0.09, length.out = 10)
North_Sask_daily_b_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_b_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, b_sd = North_Sask_daily_b_sd[i])
}
#logA0_sd
North_Sask_daily_logA0_sd <- seq(0.2, 0.5, length.out = 10)
North_Sask_daily_logA0_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_logA0_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logA0_sd = North_Sask_daily_logA0_sd[i])
}
#logn_sd
North_Sask_daily_logn_sd <- seq(0.05, 0.25, length.out = 10)
North_Sask_daily_logn_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_logn_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, logn_sd = North_Sask_daily_logn_sd[i])
}
#Werr_sd
North_Sask_daily_Werr_sd <- seq(5, 15, length.out = 10)
North_Sask_daily_Werr_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_Werr_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, Werr_sd = North_Sask_daily_Werr_sd[i])
}
#Serr_sd
North_Sask_daily_Serr_sd <- seq(1e-06, 1e-04, length.out = 10)
North_Sask_daily_Serr_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_Serr_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, Serr_sd = North_Sask_daily_Serr_sd[i])
}
#dAerr_sd
North_Sask_daily_dAerr_sd <- seq(5, 15, length.out = 10)
North_Sask_daily_dAerr_sd_priors <- list()
for(i in 1:10){
  North_Sask_daily_dAerr_sd_priors[[i]] = bam_priors(bamdata= North_Sask_data_NED_daily, dAerr_sd = North_Sask_daily_dAerr_sd[i])
}
#Test Variants of BAM with Priors-----
#Manning
North_Sask_daily_manning_list <- list(North_Sask_daily_lower_logQ_priors,
                                      North_Sask_daily_upper_logQ_priors,
                                      North_Sask_daily_lower_A0_priors,
                                      North_Sask_daily_upper_A0_priors,
                                      North_Sask_daily_lower_logn_priors,
                                      North_Sask_daily_upper_logn_priors,
                                      North_Sask_daily_sigma_man_priors,
                                      North_Sask_daily_logA0_hat_priors,
                                      North_Sask_daily_logn_hat_priors,
                                      North_Sask_daily_logQ_sd_priors,
                                      North_Sask_daily_logA0_sd_priors,
                                      North_Sask_daily_logn_sd_priors,
                                      North_Sask_daily_Werr_sd_priors,
                                      North_Sask_daily_Serr_sd_priors,
                                      North_Sask_daily_dAerr_sd_priors)
manning_prior_list <- bam_settings()$paramnames[c(1:6, 13, 18, 19, 20, 24:28)]
for(p in 11:15){
  dir.create(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man = bam_estimate(bamdata = North_Sask_data_NED_daily, variant = "manning", bampriors = North_Sask_daily_manning_list[[p]][[i]])
    save(man, file = paste0(manning_prior_list[p], "_man_", i, ".RData"))
  }
}

#AMHG
North_Sask_daily_amhg_list <- list(North_Sask_daily_lower_logQ_priors,
                                   North_Sask_daily_upper_logQ_priors,
                                   North_Sask_daily_lower_logQc_priors,
                                   North_Sask_daily_upper_logQc_priors,
                                   North_Sask_daily_lower_logWc_priors,
                                   North_Sask_daily_upper_logWc_priors,
                                   North_Sask_daily_lower_b_priors,
                                   North_Sask_daily_upper_b_priors,
                                   North_Sask_daily_sigma_amhg_priors,
                                   North_Sask_daily_logQc_hat_priors,
                                   North_Sask_daily_logWc_hat_priors,
                                   North_Sask_daily_b_hat_priors,
                                   North_Sask_daily_logQ_sd_priors,
                                   North_Sask_daily_logQc_sd_priors,
                                   North_Sask_daily_logWc_Sd_priors,
                                   North_Sask_daily_b_sd_priors,
                                   North_Sask_daily_Werr_sd_priors)
amhg_prior_list <- bam_settings()$paramnames[c(7:12, 14:17, 21:24, 26)]
for(p in 1:15){
  dir.create(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", amhg_prior_list[p]))
  setwd(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", amhg_prior_list[p]))
  for(i in 1:10){
    amhg = bam_estimate(bamdata = North_Sask_data_NED_daily, variant = "amhg", bampriors = North_Sask_daily_amhg_list[[p]][[i]])
    save(amhg, file = paste0(amhg_prior_list[p], "_amhg_", i, ".RData"))
  }
}


#MAN-AMHG
North_Sask_daily_man_amhg_list <- list(North_Sask_daily_lower_logQ_priors,
                                       North_Sask_daily_upper_logQ_priors,
                                       North_Sask_daily_lower_A0_priors,
                                       North_Sask_daily_upper_A0_priors,
                                       North_Sask_daily_lower_logn_priors,
                                       North_Sask_daily_upper_logn_priors,
                                       North_Sask_daily_lower_logQc_priors,
                                       North_Sask_daily_upper_logQc_priors,
                                       North_Sask_daily_lower_logWc_priors,
                                       North_Sask_daily_upper_logWc_priors,
                                       North_Sask_daily_lower_b_priors,
                                       North_Sask_daily_upper_b_priors,
                                       North_Sask_daily_sigma_man_priors,
                                       North_Sask_daily_sigma_amhg_priors,
                                       North_Sask_daily_logQc_hat_priors,
                                       North_Sask_daily_logWc_hat_priors,
                                       North_Sask_daily_b_hat_priors,
                                       North_Sask_daily_logA0_hat_priors,
                                       North_Sask_daily_logn_hat_priors,
                                       North_Sask_daily_logQ_sd_priors,
                                       North_Sask_daily_logQc_sd_priors,
                                       North_Sask_daily_logWc_Sd_priors,
                                       North_Sask_daily_b_sd_priors,
                                       North_Sask_daily_logA0_sd_priors,
                                       North_Sask_daily_logn_sd_priors,
                                       North_Sask_daily_Werr_sd_priors,
                                       North_Sask_daily_Serr_sd_priors,
                                       North_Sask_daily_dAerr_sd_priors)
man_amhg_prior_list <- bam_settings()$paramnames
for(p in 1:15){
  dir.create(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", manning_prior_list[p]))
  setwd(paste0(home_dir, "/North_Sask/Prior_Sensitivity/", manning_prior_list[p]))
  for(i in 1:10){
    man_amhg = bam_estimate(bamdata = North_Sask_data_NED_daily, variant = "manning_amhg", bampriors = North_Sask_daily_manning_list[[p]][[i]])
    save(man_amhg, file = paste0(man_amhg_prior_list[p], "_man_amhg_", i, ".RData"))
  }
}

#Number of Transducers------

#Spacing of Transducers------

#Plot Hydrographs--------

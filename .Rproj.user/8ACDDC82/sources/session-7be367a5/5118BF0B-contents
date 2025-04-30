#### FUNCTION TO CREATE THE SHAPEFILE OF OCCURRENCES ####
# AUTHOR: ANDRES CAMILO MENDEZ ALZATE
###############################################


create_occ_shp <- function(data, file_output ,shp_output, validation, mask_pth , occName){
  
  cat("Importing data base \n")
  msk <- terra::rast(mask_pth)
  
  Occ <- data

  if("status" %in% names(Occ)){
    Occ <- Occ %>%
      dplyr::filter(., status == "G") %>% 
      dplyr::select(., "Longitude", "Latitude", one_of(c("Y", "ensemble")))  
  } else{
    Occ <- Occ  %>% 
      dplyr::select(., "Longitude", "Latitude", one_of(c("Y", "ensemble")))
  }
  
  
  
 
  Occ <- Occ[which(Occ$Y == occName),]
  
  cat("Removing coordiantes on the ocean/sea \n")
  Occ <- Occ[which(!is.na(terra::extract(x = msk, y = Occ[,c("Longitude", "Latitude")])[,2] )),]
  
  cat("Removing duplicated coordinates \n")
  
  Occ <- Occ[!duplicated(terra::cellFromXY(msk, Occ[, c("Longitude", "Latitude")])), ]
  
  

  #save occurrences in csv format
  write.csv(Occ, file_output, row.names = FALSE)
  
  
  coordinates(Occ) <- ~Longitude+Latitude
  crs(Occ)  <- crs(msk)
  #Occ@bbox <- matrix(raster::extent(msk), ncol = 2, byrow = T)
  cat("Saving occurrences \n")
  shapefile(Occ, shp_output, overwrite = TRUE)
  #save the same file but into the results folder
  
  
  
  #writeOGR(Occ, paste0(occDir,"/Occ.shp"), "Occ", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  cat(">>> Total number of occurrences:", nrow(Occ), " \n")
  

  
}







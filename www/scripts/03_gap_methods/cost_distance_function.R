
cost_dist_function <-  function(  
                                 friction, 
                                 mask, 
                                 Occ, 
                                 sdm_path,
                                 cost_out_path){
  
 

  if(!file.exists(cost_out_path)){
   
     msk <- raster::raster(mask)
     
     if("status" %in% names(Occ)){
       Occ <- Occ %>% 
         dplyr::filter(status == "G") %>% 
         dplyr::select(Longitude, Latitude)
     }
     
     coordinates(Occ) <- ~Longitude+Latitude
     raster::crs(Occ)  <- raster::crs(msk)
     

  
    cat("Calculating cost distance raster \n")
    #p <- shapefile(paste0(occDir, "/Occ.shp"))
    r <- raster(friction) %>%
      raster::crop( x =., y = raster::extent(msk)  ) %>% 
      raster::mask(., msk)
    
    t <- gdistance::transition(r, function(x) 1/mean(x), 8) 
    t <- gdistance::geoCorrection(t) 
    
    cost_dist <- gdistance::accCost(t, Occ) 
    cost_dist[which(cost_dist[] == Inf)] <- NA
    
    #Normalize cost distance rater to be in 0-1 scale
    
    #mask cost distance raster with SDM raster
    sdm <- raster(sdm_path)
    cost_dist <- raster::mask(cost_dist, sdm)
    #idenfy outliers using IQR 
    qls <- raster::quantile(cost_dist, na.rm = T)
    up_limit <- qls[4] + (1.5* (qls[4] - qls[2]))
    cost_dist <- cost_dist/up_limit
    cost_dist[which(cost_dist[] > 1 )] <- 1
    
    cost_dist[which(is.na(cost_dist[]) & !is.na(sdm[]))] <- max(cost_dist[],na.rm=T)
    
    cost_dist <- cost_dist*sdm
    rm(t)
  
    
  raster::writeRaster(cost_dist, filename = cost_out_path, overwrite= T)

    
  } 
  
  return(NA)
  cat('Done... \n')
}


# 

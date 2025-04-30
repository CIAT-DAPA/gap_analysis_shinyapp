
######********FUCNTION TO CALCULATE A GAP SCORE USING DELAUNAYS TRIANGULATION*****************#################
######SCRIPT CREADO POR ANDRES CAMILO MENDEZ
######THIS SCRIPT COMES WHIT ABSOLUTLELY NO WARRANTY


calc_delaunay_score <- function(
    coreDir, 
    ncores = NULL, 
    validation = FALSE, 
    pnt = NULL){
  
  
  
  cat(">>> Generating Delaunays triangulation \n \n")
  
  if(validation == FALSE){
    cat(">>> Initializing process to calculte a gap score using Delaunays Triangulation \n \n ")
  } else {
    cat(">>> Initializing process to Validate a gap score using Delaunays Triangulation through +++ LOGCV +++ \n \n ")
  }
  
   
 
   
  
  sdmDir <- switch(as.character(validation),
                   "FALSE" = list.files(paste0(coreDir, "/species_distribution"), full.names = T, pattern = "_sdm_median.tif$"),
                   "TRUE"  = list.files(paste0(coreDir, "/02_sdm_results/"), full.names = T, pattern = "_sdm_median.tif$"))
  
  outDir <- switch(as.character(validation),
                   "FALSE" = paste0(coreDir , "/gap_scores" ),
                   "TRUE"  = paste0( coreDir,"/03_gap_scores"))
  
  delaDir <- switch(as.character(validation),
                    "FALSE" = paste0(coreDir , "/gap_scores"),
                    "TRUE"  =  paste0( coreDir,"/03_gap_scores"))
  
  occDir <- switch(as.character(validation),
                   "FALSE" = paste0(coreDir, "/species_distribution/occurrences.shp"),
                   "TRUE"  = paste0(coreDir,"/01_occurrences/occurrences.shp"))
  
  cat(">>> Importing Delaunays triangulation \n \n")
  
  if(file.exists(paste0(delaDir, "/delaunay/raw_delaunay.shp")) == FALSE){
    
    delaunaypolygons(x = shapefile(occDir) , outdir = delaDir)
    
    cat("Delaunay was successfully created \n \n")
  } else {
    cat("Delaunay is already created \n \n")
  }
  
  cat("Importing SDM raster \n")  
  SDM <- raster(sdmDir)
  
  SDM[!is.na(SDM[])] <- 1
  
  delanuay <- shapefile(paste0(delaDir, "/delaunay/raw_delaunay.shp"))
  
  SDM <- raster::crop(SDM, raster::extent(delanuay))
  
  ext_sdm <- extent(SDM)
  res_sdm <- res(SDM)
  coord_sys <- crs(SDM)
  # Change path after
  delanuay@bbox <- as.matrix(ext_sdm)
  
  if(length(grep("area", names(delanuay@data))) == 0){
    delanuay$area <- raster::area(delanuay)/1000000
  }
  
  if(sum(which(delanuay@data$area == 0)) > 0){
    delanuay <- delanuay[-which(delanuay@data$area == 0), ]
    
  }
  
  if( length(grep("centroid", names(delanuay@data))) == 0){
    
    centroids <- sp::coordinates(delanuay)
    delanuay$centroid.x <- centroids[,1]
    delanuay$centroid.y <- centroids[,2]
    
  }
  
  areas <- lapply(1:nrow(delanuay@data), function(i){
    
    shp <- delanuay[i,]
    
    inters <- raster::mask(SDM, shp, inverse = F)
    
    ar <- raster::area(inters) %>%
      raster::mask(., inters) 
    
    sm <- sum(ar[!is.na(ar[])])
    return(sm)
  }) %>% 
    unlist()
  
  max.area <- max(areas, na.rm = T)
  
  #max.area <- max(raster::area(i_shp)/1000000) # al dividir por 1000000 la unidad es kilometro^2 
  
  cat(paste(">>> Max.area = ", max.area, "\n \n"))
  
  # delanuay  <- delanuay[delanuay$area >= 0.02,]
  
  r <- raster()
  raster::extent(r) <- ext_sdm
  raster::crs(r) <- coord_sys
  raster::res(r) <- res_sdm
  r[is.na(r[])]<- 0
  
  cat("Calculating distances inside each triangulation... \n \n")
  
  
    pb <- txtProgressBar( style = 3)
    
    delaDist_list <- lapply(1:(length(delanuay)), function(x){
      
      setTxtProgressBar(pb, x/length(delanuay))
      
      vertex_1 <- delanuay@polygons[[x]]@Polygons[[1]]@coords[1,] 
      vertex_1 <-  SpatialPoints(data.frame( x = vertex_1[1], y = vertex_1[2]), proj4string = coord_sys)
      
      vertex_2 <- delanuay@polygons[[x]]@Polygons[[1]]@coords[2,]
      vertex_2 <-  SpatialPoints(data.frame( x = vertex_2[1], y = vertex_2[2]), proj4string = coord_sys)
      
      vertex_3 <- delanuay@polygons[[x]]@Polygons[[1]]@coords[3,]
      vertex_3 <-  SpatialPoints(data.frame( x = vertex_3[1], y = vertex_3[2]), proj4string = coord_sys)
      
      centroid <- delanuay@data[x, c("centroid.x", "centroid.y")]
      
      #cat(paste("Processing feature: " , x, "\n"))
      cr <- raster::crop( x = r , y = extent(delanuay[x, ]) )
      rr <- raster::mask(cr, delanuay[x,])
      
      if(all(is.na(rr[])) == TRUE){rr[is.na(rr[])] <- 0}
      
      xy <- SpatialPoints(centroid, proj4string = coord_sys)
      
      dRas <- raster::distanceFromPoints(rr, xy)
      
      dVer1 <- raster::distanceFromPoints(rr, vertex_1)  + rr # dist to nearest vertex
      dVer2 <- raster::distanceFromPoints(rr, vertex_2) + rr # dist to nearest vertex
      dVer3 <- raster::distanceFromPoints(rr, vertex_3)+ rr # dist to nearest vertex
      
      dVer <- min(raster::stack(dVer1, dVer2, dVer3))
      dVer <- dVer/max(dVer[], na.rm = T)
      
      #cat( paste("max value distances: ",  max(values(dRas),na.rm = T)," \n \n" ))
      
      dRas <-  dRas + rr
      dRas <- dRas/max(dRas[], na.rm = T)
      
      s <- delanuay$area[x]/max.area
      #cat(paste("El area relativa para este feature es: ", s , "\n \n "))
      if(s > 1){s <- 1}
      sRas <- rr + s
      if(round(res(sRas)[1],8) !=  0.04166667){
        cat("Something went wrong \n") # CAMBIAR
      }
      final_rast <- sRas * (1 - dRas) * dVer
      
      #cat(paste("resolution:", res(sRas)[1], "****** \n \n \n "))
      return(dela_score = final_rast)
      #return(list(dist_to_centroid = dRas, dist_near_vertex = dVer, area_score = sRas))
    })
    close(pb)
    g <- gc(); rm(g); removeTmpFiles() 
  
  
  ### merge everything
  cat("Calculating gap score \n \n ")
  
  cat("+++ Merging rasters \n \n")
  
  out_dela_score <- do.call(merge, delaDist_list)
  
  # x <- out_dist_centroid
  # a <- out_dist_vertex
  # p <- out_area_relativa * (1 - x) * a
  
  # p2 <- out_area_relativa * mean((1-x),a)
  # p3 <- out_area_relativa * sqrt((1-x),a)
  
  #rm(x, a); g <- gc(); rm(g); removeTmpFiles(h = 24)
  
  # writeRaster(p, paste0(outDir, "/delaunay_aux.tif"), format = "GTiff", overwrite= T)
  # writeRaster(p2, paste0(outDir,"/delanuay_probs_mean.tif"), format = "GTiff", overwrite= T)
  # writeRaster(p3, paste0(outDir,"/delanuay_probs_sqrt.tif"), format = "GTiff", overwrite= T)
  
  SDM <- raster(sdmDir)
  # 
  # sdm_crop <- raster::crop(x = SDM, y = extent(p))
  # p_masked <- raster::mask(x = p, mask = sdm_crop)
  
  
  ###### calculating the gap score
  # gap_score <- sdm_crop * p_masked
  # extent(gap_score) <- extent(SDM)
  
  cat("Calculating scores for delaunay's outside points \n \n")
  #### scoring points out side of the delaunay triangulation
  
  delanuay <- rgeos::gBuffer(delanuay, byid=TRUE, width=0)
  dela_bound <-  rgeos::gUnaryUnion(delanuay, id = NULL)
  
  vertex <- dela_bound@polygons[[1]]@Polygons[[1]]@coords[]
  vertex <- vertex[-1, ]
  vertex <- SpatialPoints(data.frame(x = vertex[,1], y = vertex[,2]), proj4string = coord_sys)
  
  out <- raster::mask(x = SDM, mask = dela_bound, inverse = T)
  
  
  out2 <- out
  out2[!is.na(out2[])] <- 1
  
  dist_out <- raster::distanceFromPoints(out2 , vertex)
  dist_out <- dist_out * out2
  #### normalizing distances
  dist_out_norm <- dist_out/max(dist_out[], na.rm = T)
  ### calculating new score
  # dist_out_score <-  dist_out_norm * out
  
  #gap_score<- raster(paste0(gap_outDir,"/delanuay_gap_score.tif"))
  
  gap_score <- merge(dist_out_norm, out_dela_score)
  
  writeRaster(gap_score, paste0(outDir, "/delaunay_network.tif"), format = "GTiff", overwrite = T)
  
  gap_score <- raster::mask(gap_score, SDM)
  #idenfy outliers using IQR 
  qls <- raster::quantile(gap_score, na.rm = T)
  up_limit <- qls[4] + (1.5* (qls[4] - qls[2]))
  gap_score <- gap_score/up_limit
  gap_score[which(gap_score[] > 1 )] <- 1
  
  gap_score[which(is.na(gap_score[]) & !is.na(SDM[]))] <- max(gap_score[],na.rm=T)
  
  gap_score <- gap_score*SDM
  
  cat(paste("Done", "\n \n"))
  writeRaster(gap_score, paste0(outDir, "/network_score.tif"), format = "GTiff", overwrite = T)
  
  
}
# Example:
# a <- c("americas", "world")
# g <- c("mesoamerican", "andean")
# c <- "common_bean"
# lvl <- "lvl_1"
# results_dir <- "//dapadfs/Workspace_cluster_9/gap_analysis_landraces/runs" 
# delaunay_scoring(baseDir = baseDir, area = a[1], group = g[1], crop = c, lvl = "lvl_1", ncores = 10, validation = FALSE , pnt = NULL )

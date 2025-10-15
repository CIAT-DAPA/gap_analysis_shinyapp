pseudo_occ <- function(xy, clim_dir, aux_clim_dir, tms, nu){
  
  ## funcion para seleccionar pseudo-PRESENCIAS cuando el numero de occ es bajo,
  ## se seleccionar a partir de un SVM que indica pixeles con condiciones muy similares a las occ
  
  if(tms != 1){
    
    if("status" %in% names(xy)){
      xy$status <- NULL
    }
    if("predicted" %in% names(xy)){
      xy$predicted <- NULL
    }
    if("database_id" %in% names(xy)){
      xy$database_id <- NULL
    }
    if("source_db" %in% names(xy)){
      xy$source_db <- NULL
    }
    
    message("Generating pseudo-presences")
    ############ OCVMprofiling ####################
    clim_fls <- list.files(clim_dir, pattern = ".tif$", full.names = T)
    aux_fls  <- list.files(aux_clim_dir, pattern = ".tif$", full.names = T)
    
    if(length(aux_fls) == 0){
      varstack <- terra::rast(clim_fls)
    }else{
      varstack <- terra::rast(c(clim_fls, aux_fls))
    }
    
    proj <- terra::as.data.frame( varstack, xy = T ) %>% 
      drop_na() 
    
    nms <- names(proj[, -c(1,2)])
    
    mat        <- xy %>% 
      dplyr::mutate(Y = 1) %>% 
      dplyr::select(Y, everything(.), -Latitude, -Longitude, all_of(nms) )
    
    mod        <- e1071::svm(mat[, -1], y = NULL, type = "one-classification", nu = nu)
    pre        <- predict(mod, proj %>% dplyr::select(-x,-y))
    
    svm_msk <- terra::vect(data.frame(lon = proj$x, lat = proj$y, z = as.numeric(pre)), geom=c("lon", "lat"), crs = terra::crs(varstack))
    svm_msk <- terra::rasterize(svm_msk, varstack[[1]], field = "z") 
    svm_msk <- terra::ifel(svm_msk == 0, NA, svm_msk)
    rm( proj, mat, mod, pre)
    
    pseudo_pres_xy <- terra::spatSample(svm_msk, size = nrow(xy)*tms, na.rm = T, xy =T, values = F)
    
    to_ret <- data.frame(Y = unique(xy$Y)[1],  
                         Latitude = pseudo_pres_xy[,2], 
                         Longitude = pseudo_pres_xy[,1],
                         terra::extract(varstack, pseudo_pres_xy))
    
    to_ret <- dplyr::bind_rows(xy, to_ret)
  }else{
    to_ret <- xy
  }
  #varstack <- terra::mask(varstack, svm_msk)
  
  
  return(to_ret)
  
}

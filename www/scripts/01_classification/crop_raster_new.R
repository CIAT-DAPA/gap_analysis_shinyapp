get_c_mask_df <- function(w_msk_path, shp_c, out_path){
  
  w_mask <- terra::rast(w_msk_path) 
  
  if(any(c("sf", "sp", "SpatVector") %in% class(shp_c))){
    c_mask <- w_mask %>% 
      terra::crop(., terra::ext(shp_c)) %>% 
      terra::mask(., shp_c) 
  }else if("SpatRaster" %in% class(shp_c)){
    c_mask <- shp_c %>% 
      terra::resample(., w_mask %>% terra::crop(shp_c))
  }else if(any(c("RasterLayer", "RasterStack", "RasterBrick") %in% class(shp_c))){
    c_mask <- terra::rast(shp_c)%>% 
      terra::resample(.,  w_mask %>% terra::crop(shp_c))
  }else{
    stop("Invalid shapefile format")
  }
   
  terra::writeRaster(c_mask, out_path, overwrite = T)
  
  c_mask_df <- terra::as.data.frame(c_mask, xy = T, na.rm = F) %>% 
    data.table::as.data.table(.)
  
  names(c_mask_df)[3] <- "mask_world"
  
  c_mask_df$cellID <- terra::cellFromXY(w_mask, as.matrix(c_mask_df[, c("x", "y")]))
  
  data.table::setorder(c_mask_df, cellID)
  
  rm(w_mask, c_mask)
  
  c_mask_df <- na.omit(c_mask_df, cols = "mask_world" )
  c_mask_df <- c_mask_df[, .(cellID, x, y, mask_world)]
  
  
  return(list(c_mask_df = c_mask_df))
  
}


extract_values <- function(c_mask_df, tables_path, c_mask, out_path){
  
  
  
  fls <- list.files(tables_path, full.names = T)
  
  mn <-  min(c_mask_df$cellID) #3431# # 
  mx <-  max(c_mask_df$cellID)#4912940
  
  ids_range_tbl <- tibble(paths = fls) %>% 
    dplyr::mutate(min_range = stringr::str_extract(fls, pattern = "_[0-9]+") %>% stringr::str_replace(., "_", "") %>% as.numeric(),
                  max_range = stringr::str_extract(fls, pattern = "-[0-9]+") %>% stringr::str_replace(., "-", "") %>% as.numeric()) %>% 
    dplyr::arrange(min_range) %>% 
    dplyr::mutate(min_in = ifelse( mn >= min_range &  mn <= max_range, TRUE, FALSE  ),
                  max_in = ifelse( (mx >= min_range & mx <= max_range)  , TRUE, FALSE) )
  
  pos <- which(ids_range_tbl$min_in):which(ids_range_tbl$max_in)
  
  ids_range_tbl <- ids_range_tbl[pos, ]
  
  
  if(nrow(ids_range_tbl)> 1){
    
    df_extracted <- lapply(1:nrow(ids_range_tbl), function(i){
      
      rw <- ids_range_tbl[i, ]
      
      to_extract <- readRDS(rw$paths)
      to_extract <- to_extract[, -c("x", "y", "mask_world")]
      
      to_ret <- to_extract[c_mask_df[ cellID >= rw$min_range & cellID <= rw$max_range,   ] , on = .(cellID)]
      rm(to_extract)
      
      
      
      data.table::setcolorder(to_ret, neworder = c(c("cellID", "x", "y", "mask_world" ), setdiff(names(to_ret), c("cellID", "x", "y", "mask_world" ))))
      
      return(to_ret)
      
    }) %>% 
      data.table::rbindlist(.)
    
    
  }else{
    
    to_extract <- readRDS(ids_range_tbl$paths)
    to_extract <- to_extract[, -c("x", "y", "mask_world")]
    
    df_extracted <- to_extract[c_mask_df[ cellID >= ids_range_tbl$min_range & cellID <= ids_range_tbl$max_range,   ] , on = .(cellID)]
    rm(to_extract)
    
    
    
    data.table::setcolorder(df_extracted, neworder = c(c("cellID", "x", "y", "mask_world" ), setdiff(names(df_extracted), c("cellID", "x", "y", "mask_world" ))))
    
    
  }
  
  ##create cropped rasters
  
  nms <- names(df_extracted)[-c(1,2,3,4)]
  nms <- nms[!grepl("monthCountByTemp10|ethnicity", nms)]
  
  for(i in nms){
    vrs <- c("x", "y", i)
    rf <- terra::rast(df_extracted[, ..vrs], crs = terra::crs(c_mask)) %>% 
      terra::resample(., c_mask)
    
    stopifnot("different raster resolution. " = all(terra::res(rf) == terra::res(c_mask)))
    
    terra::writeRaster(rf, paste0(out_path, "/", i, ".tif"), overwrite = T)
  }
  
  
  return(df_extracted)
  
}



# c_mask_df <- get_c_mask_df(w_msk_path = "www/masks/mask_world.tif",
#                            shp_c = terra::vect("www/world_shape_simplified/KEN.shp"),
#                            out_path)
# 
# 
# 
# c_mask_extracted <- extract_values(c_mask_df = c_mask_df$c_mask_df, 
#                                    tables_path = "www/generic_raster_tables/",
#                                    c_mask = c_mask_df$c_mask,
#                                    out_path)
# Create sample (background, occurrences, swd) for SDMs
# Chrystiam Sosa, Julian Ramirez-Villegas, Harold Achicanoy, Andres camilo mendez
# CIAT, Nov 2018


# Profiling function
pseudoAbsences2 <- function(xy, 
                            varstack, 
                            nu = 0.5, 
                            exclusion.buffer = 20000, 
                            tms = 10){
  
  if(is.character(xy$Y)){
    xy$Y <- 1
  }
  
  xy_vec <- xy %>%
    dplyr::select(Longitude, Latitude, Y) %>% 
    terra::vect(., geom = c("Longitude", "Latitude"),  crs = terra::crs(varstack))
  
  
  #buffer of 20km width
  spol <- terra::buffer(xy_vec, width = exclusion.buffer) %>%
    terra::aggregate(., dissolve = T)
    #sf::st_as_sf()
  #remove cell within bufffers  
  to_remove <- terra::mask(varstack[[1]], spol, inverse = T) 
  varstack  <- terra::mask(varstack, to_remove)
  
  if('matrix' %in% class(varstack[[1]][!is.na(varstack[[1]])])){
    tot_cells <-  length(as.numeric(varstack[[1]][!is.na(varstack[[1]])]))
  }else{
    stop('Unssuported class in pseudoAbsences2 function.')
  }
  
  #test if sample sizes matches env raster non-NA cells
  if(tot_cells > tms * nrow(xy)){
    n_sample <- tms*nrow(xy)
  }else{
    n_sample <- round(tot_cells*0.95, 0)
  }
  
  PA.r <- biomod2::bm_PseudoAbsences(resp.var = xy_vec,
                                     expl.var = varstack,
                                     nb.rep = 1,
                                     nb.absences = n_sample,
                                     strategy = 'random', #"random"
                                     #seed.val = 1234,
                                     sre.quant = 0.025)
  
  pseudo <- dplyr::bind_cols(PA.r$xy, PA.r$env) %>% 
    tidyr::drop_na() %>% 
    dplyr::rename("Longitude" = x, "Latitude" = y)
   
  
  
  return(pseudo)
}


# Pseudo-absences options: pa_method
# "ntv_area_ecoreg": Native Area cropped by EcoRegion
# "ecoreg": EcoRegion
# "all_area": Full extent of mask
pseudoAbsences_generator <- function(data,
                                     climDir, 
                                     aux_dir,
                                     clsModel, 
                                     overwrite = F, 
                                     correlation = 0,
                                     occName,
                                     pa_method = "ntv_area_ecoreg",
                                     bg_out_path,
                                     smd_var_selected_path,
                                     mask_path,
                                     ecoreg_path = "www/masks/World_ELU_2015_5km.tif"){

  cat("Loading raster files","\n")
  
  fls <- list.files(climDir, pattern =  ".tif$", full.names = T)
  
  
  # Create background if it does not exist
  if (!file.exists(bg_out_path) | overwrite){
    msk <- terra::rast(mask_path)
    
    cat("Processing:", paste(occName), "\n")
    #file_path archivo completo
    spData            <- data#read.csv(file_path, header = T)
    #spData[,clsModel] <- tolower(spData[,clsModel])
    spData            <- spData[which(spData[,clsModel]== occName),]
    
    
    if("status" %in% names(spData)){
      spData$status <- NULL
    }
    if("predicted" %in% names(spData)){
      spData$predicted <- NULL
    }
    if("database_id" %in% names(spData)){
      spData$database_id <- NULL
    }
    if("source_db" %in% names(spData)){
      spData$source_db <- NULL
    }
   
    cat("Creating random Pseudo-absences points using: ", pa_method, "method \n")
    
    climLayers <-  terra::rast(fls)
   
    #Remove variables that are causing problems
    #vars_to_remove <- c("Yield", "Production", "Harvested", "drymonths_2_5_min", "monthCountByTemp10", "ethnicity")
    
    if(pa_method == "ecoreg"){
      
      elu     <- terra::rast(ecoreg_path) %>% 
        terra::crop(., terra::ext(msk)) %>% 
        terra::mask(., msk)
      #print(head(spData))
      regions <- terra::extract(x = elu, y = spData[,c("Longitude","Latitude")])[,2]
      regions <- sort(unique(regions))
      elu[!(elu %in% regions)] <- NA # Exclude ecoregions that are not in occurrence data
      #unsuit_bg <- OCSVMprofiling2(xy = spData, varstack = climLayers)
      random_bg <- pseudoAbsences2(xy = spData, 
                                   varstack = climLayers %>%  terra::mask(., elu),
                                   nu = 0.05,
                                   exclusion.buffer = 15000,
                                   tms = 10)
      
      
    }
    if(pa_method == "all_area"){
      #unsuit_bg <- OCSVMprofiling2(xy = unique(spData[,c("Longitude","Latitude")]), varstack = climLayers)
      random_bg <- pseudoAbsences2(xy = spData, 
                                   varstack = climLayers,
                                   nu = 0.05,
                                   exclusion.buffer = 15000,
                                   tms = 10)
     
    }
    
    stopifnot("climLayers names does not match random_bg names" = all(names(random_bg[, -c(1,2)]) %in% names(climLayers) ))
    
    cat( nrow(random_bg), "pseudo-absences generated for n =", nrow(spData), "presences\n")
    
    aux_fls <- list.files(aux_dir, full.names = T, pattern = ".tif$")
    if(length(aux_fls) != 0){
     fls <- c(fls, aux_fls) 
     climLayers <- terra::rast(fls)
    }
    
    
    random_bg <- random_bg %>% 
      #dplyr::bind_cols(., terra::extract(climLayers, .)) %>% 
      dplyr::mutate(Y = occName) %>% 
      dplyr::select(Y, Latitude, Longitude, everything(.)) 
    
    
    # Extract variable data
      
    
    swdSample <- dplyr::bind_rows(spData , random_bg) 
    
    rm(spData, climLayers)
    
    # Using choose variables algorithms (Correlation, VIF, or PCA + VIF)
    if(correlation == 1){
      cat("Using Pearson correlation approach\n")
      descrCor       <- cor(swdSample[,-c(1:3)])
      highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .75)
      swdSample      <- swdSample[,!colnames(swdSample) %in% (colnames(descrCor)[highlyCorDescr]) ]
    }
    
    if(correlation == 2){
      cat("Using VIF approach\n")
      descrCor       <- usdm::vifstep(swdSample[,-c(1:3)], th = 5)
      highlyCorDescr <- descrCor@excluded
      swdSample      <- swdSample[,!colnames(swdSample) %in% highlyCorDescr]
    }
    
    if(correlation == 3){
      cat("Using PCA + VIF approach","\n")
      z <- FactoMineR::PCA(X = swdSample[,-c(1:3)], ncp = 5, scale.unit = T, graph = F)
      # Selecting a number of components based on the cumulative ratio of variance which has more than 70%
      ncomp <- as.numeric(which(z$eig[,ncol(z$eig)] >= 70)[1])
      vars  <- rownames(z$var$cos2)[unlist(lapply(X = 1:nrow(z$var$cos2[,1:ncomp]), FUN = function(r){
        if(length(which(z$var$cos2[r,1:ncomp] >= .15)) > 0){ res <- T } else { res <- F }
        return(res)
      }))]
      descrCor       <- usdm::vifstep(swdSample[,vars], th = 10)
      highlyCorDescr <- descrCor@excluded
      swdSample      <- swdSample[,!colnames(swdSample) %in% highlyCorDescr]
    }
    cat("Saving csv files","\n")
    
    var_names <- swdSample %>% 
      dplyr::select(-Y, -Latitude, -Longitude) %>% 
      names(.)
    
    
    rm(swdSample)
    
    write.csv(random_bg, bg_out_path, row.names = F)
  
    write.csv(x = var_names, file =  smd_var_selected_path, row.names = F)
    
  } else {
    
    cat("PseudoAbsence file already created, importing it...\n")
    random_bg <- read.csv(bg_out_path, header = T)
    
  }
  return(random_bg)
}



msk <- terra::rast(paths$mask_path)

knl <- raster_kernel(mask          = paths$mask_path, 
                     occDir        = paths$occ_shp, 
                     out_path      = paths$knl_out_path, 
                     kernel_method = 2, 
                     scale         = T)


knl <- terra::rast(knl)
knl[knl == 1] <- NA

c1 <- round(input$n_points)


to_sample<- knl %>% 
  terra::as.data.frame(., xy = T) %>% 
  tidyr::drop_na() 

names(to_sample)[length(names(to_sample))] <- "kernel"

set.seed(1234)

if(!file.exists(paste0(paths$gap_valDir, "/buffer_coordinates_to_exclude.csv"))){
  points <- to_sample %>% 
    dplyr::filter(kernel == 3) %>% 
    dplyr::slice_sample(n = c1)
  
  write.csv(points, paste0(paths$gap_valDir, "/buffer_coordinates_to_exclude.csv"), row.names = F)
  
}else{
  points <- read.csv(paste0(paths$gap_valDir, "/buffer_coordinates_to_exclude.csv"))
}




rm(to_sample, knl)
#input$n_points
summ_all<-lapply(1:input$n_points,function(i){
  
  cat(">>>>>>processing point: ", i, "\n")
  
  
  out_pnt_path <- paste0(paths$gap_valDir, "/pnt", i)
  
  val_occ_out  <- paste0(out_pnt_path, "/01_occurrences") 
  val_sdm_out  <- paste0(out_pnt_path, "/02_sdm_results")
  val_gap_out  <- paste0(out_pnt_path, "/03_gap_scores")
  
  cv_fold_out     <- paste0(val_sdm_out, "/sdm_cv_folds") 
  sdm_occ_path    <- paste0(val_sdm_out,"/", paths$occName, "_sdm_median.tif")
  dela_out_folder <- paste0(val_gap_out, "/delaunay")
  cost_out_path   <- paste0(val_gap_out,"/cost_dist_score.tif")
  dela_out_path   <- paste0(val_gap_out,"/network_score.tif")
  envi_out_path   <- paste0(val_gap_out,"/environmental_score.tif")
  knl_out_path    <- paste0(val_gap_out,"/kernel.tif")
  summ_out_path   <- paste0(out_pnt_path,"/summary_raw.csv")
  
  if(!file.exists(out_pnt_path)){dir.create(out_pnt_path, recursive = T)}
  if(!file.exists(val_occ_out)){dir.create(val_occ_out, recursive = T)}
  if(!file.exists(val_sdm_out)){dir.create(val_sdm_out, recursive = T)}
  if(!file.exists(val_gap_out)){dir.create(val_gap_out, recursive = T)}
  
  if(!file.exists(cv_fold_out)){dir.create(cv_fold_out, recursive = T)}
  if(!file.exists(dela_out_folder)){dir.create(dela_out_folder, recursive = T)}
  
  
  write.csv(points[i , ], paste0(val_occ_out, "/buffer_centroid_coords.csv"), row.names = F)
  
  buff <- points[i , 1:2] %>% 
    terra::vect(., geom = c("x", "y"), crs = crs(msk)) %>% 
    terra::buffer(., width = input$bf_size*1000) 
  
  terra::writeVector(buff, filename = paste0(val_occ_out,"/buffer_shp.shp"), overwrite=TRUE)
  
  occ <- resources$spData %>% 
    dplyr::mutate(ID = 1:nrow(.))
  
  to_exclude <- terra::intersect(terra::vect(occ , geom = c("Longitude", "Latitude"), crs = crs(msk)),
                                 buff) %>% 
    as.data.frame()
  
  write.csv(to_exclude, paste0(val_occ_out,"/excluded_occ.csv"), row.names = F)
  
  occ <- occ[-to_exclude$ID, ]
  occ$ID <- NULL
  
  if(file.exists(paths$sdm_calibration)){
    params_tunned  <- Calibration_function(spData   = NULL,
                                           bg_data  = NULL,
                                           occName  = NULL,
                                           out_path = paths$sdm_calibration,
                                           ommit = F,
                                           use.maxnet = TRUE)
  }else{
    params_tunned <- list()
    params_tunned$beta <- 1
    params_tunned$feat <- c("linear", "quadratic", "hinge", "product")
  }
  
  if(!file.exists(sdm_occ_path)){
    sdm_maxent_approach_function(occName      = paths$occName,
                                 spData       = occ %>% 
                                   dplyr::select(-any_of("status")),
                                 bg_data      = resources$pseudo_abs,
                                 var_names    = resources$var_names,
                                 model_outDir = val_sdm_out,
                                 replic_path  = cv_fold_out,
                                 nFolds       = 3,
                                 climDir      = paths$generic_dir,
                                 clim_spReg   = paths$aux_dir,
                                 beta         = params_tunned$beta,
                                 feat         = params_tunned$features,
                                 doSDraster   = FALSE,
                                 varImp       = FALSE,
                                 validation   = TRUE)
  }
  
  cat("sdm calculated... \n")
  
  if("status" %in% names(resources$spData)) {
    occ <- occ %>% 
      dplyr::filter(status == "G") %>% 
      dplyr::select(-status)
  }
  
  terra::vect(occ , geom = c("Longitude", "Latitude"), crs = crs(msk)) %>% 
    terra::writeVector(., paste0(val_occ_out, "/occurrences.shp"), overwrite = T)
  
  
  cat("Calculating cost dist validation \n")
  
  if(!file.exists(cost_out_path)){
    cost_dist_function(
      cost_out_path = cost_out_path,
      friction      = "www/masks/friction_surface.tif",
      mask          = paths$mask_path,
      Occ           = occ,
      sdm_path      = sdm_occ_path
    )
  }
  cat("cost dist calculated \n")
  
  
  cat("Calculating delaunay score \n")
  if(!file.exists(dela_out_path)){
    calc_delaunay_score(
      coreDir = out_pnt_path, 
      ncores = NULL, 
      validation = TRUE, 
      pnt = NULL)
  }
  cat("Done \n")
  
  cat("calulating environmental score \n")
  if(!file.exists(envi_out_path)){
    calc_env_score(sdm_path    = sdm_occ_path, 
                   clus_method = "hclust_mahalanobis", 
                   gap_dir     = paths$gap_outDir,
                   out_dir     = val_gap_out,
                   occ_dir     = val_occ_out, 
                   env_dir     = paths$generic_dir, 
                   var_names   = resources$var_names,
                   n.sample    = 5000,
                   n.clust     = 6)
  }
  
  cat("Done \n")
  
  cat("Summarizing results \n")
  
  knl_val <- raster_kernel(mask      = paths$mask_path, 
                           occDir        = paste0(val_occ_out, "/occurrences.shp"), 
                           out_path      = knl_out_path, 
                           kernel_method =2, 
                           scale         = T)
  
  
  cords_dummy <- knl_val %>% 
    terra::rast() %>% 
    terra::mask(., buff, inverse = T) %>% 
    terra::as.data.frame(knl_val, xy = T) 
  
  names(cords_dummy)[length(names(cords_dummy))] <- "kernel"
  
  cords_dummy <- cords_dummy%>% 
    tidyr::drop_na() %>% 
    dplyr::filter(kernel != 1) %>% 
    dplyr::slice_sample(n = 100) %>% 
    dplyr:::select(x,y)
  
  
  cent <- points[i, 1:2]%>% 
    terra::vect(., geom = c("x", "y"), crs = crs(msk))
  
  #gp_m <- gap_score
  #gp_m2 <- gap score sin buffer
  
  gp_m <- terra::rast(c(cost_out_path, dela_out_path, envi_out_path))
  
  names(gp_m) <- terra::sources(gp_m) %>% 
    stringr::str_extract("([a-zA-Z]+_[a-zA-Z]+|[a-zA-Z]+)_score.tif",string = .) %>%
    stringr::str_replace(., ".tif", "")
  
  gp_m2  <- terra::mask(gp_m, buff, inverse = T)
  names(gp_m2) <- names(gp_m)
  
  big_rad <- input$bf_size
  
  radius <-  seq(big_rad/2, big_rad, length.out = 5) %>% round(., 0)
  
  if(!file.exists(summ_out_path)){
    sm_res <- lapply(radius, function(rd){
      
      cat("Processing radius :", rd, "\n")
      width = rd*1000
      all_buffs <- terra::buffer(terra::vect(cords_dummy, geom = c("x", "y"), crs = crs(terra::rast(msk))), width=width)
      no_gap_ls <- terra::extract(gp_m2, all_buffs) %>% 
        tidyr::drop_na() 
      
      
      scr <- terra::extract(gp_m, terra::buffer(cent, width=width )) %>% 
        tidyr::drop_na() %>% 
        dplyr::select(-ID)
      
      
      
      
      results <- lapply(unique(no_gap_ls$ID), function(i){
        
        no_gap <- no_gap_ls %>% 
          dplyr::filter(ID == i) %>% 
          dplyr::select(-ID)
        
        if(length(no_gap)!=0){
          
          ng <- data.frame( score = no_gap, observe = rep(0, nrow(no_gap) ))
          gap <- data.frame(score = scr, observe = rep(1, nrow(scr)))
          
          croc_summ <- dplyr::bind_rows(ng, gap, .id = NULL) %>% 
            tidyr::pivot_longer(., cols = -observe, names_to = "gap_score", values_to = "vals") %>% 
            dplyr::group_by(gap_score) %>% 
            dplyr::summarise(roc_summ(observe, score = vals)) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(id = i)
          
          
          
        }else{
          
          croc_summ <- data.frame(threshold = NA, auc = NA, sensi = NA, speci = NA, max.TSS = NA, id= i)
          
        }
        return(croc_summ)
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(radius = rd)
      
      
      return(results)
      
    })%>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(pnt = i) %>% 
      tidyr::drop_na()
    
    write.csv(sm_res, summ_out_path, row.names = F)
    
  }else{
    sm_res <- read.csv(summ_out_path, header = T)
  }
  
  
  
  cat("Done. \n")
  
  shinyWidgets::updateProgressBar(session = session, id = "pg", value = round(i/input$n_points,1)*100)
  
  return(sm_res)
}) %>% 
  dplyr::bind_rows()

score_mean <- function(x,li,ls, se, es){
  #x %>%  dplyr:: filter(., auc >= li & auc <= ls ) %>% dplyr::select(., score) %>% mean(., na.rm = TRUE)
  y <- max(  x$threshold[which(x$auc >= li & x$auc <= ls)], na.rm = TRUE)
  z <- mean(  x$sensi[which(x$auc >= li & x$auc <= ls)], na.rm = TRUE)
  w <- mean(  x$speci[which(x$auc >= li & x$auc <= ls)], na.rm = TRUE)
  return( c(y, z, w))
}


means_all <- summ_all %>% 
  dplyr::group_by(., pnt,gap_score ) %>% 
  dplyr::summarise(., auc.median = round(median(auc, na.rm = T),3)
                   , auc.mean = round(mean(auc, na.rm = T), 3)
                   #, skewness = round(skew(auc), 3)
                   , auc.sd = round(sd(auc, na.rm = TRUE), 3)
                   , lower.ic = round(t.test(auc, conf.int = TRUE, conf.level = 0.95, na.rm = TRUE)$conf.int[1], 3)
                   , upper.ic = round(t.test(auc, conf.int = TRUE, conf.level = 0.95, na.rm = TRUE)$conf.int[2], 3)
                   , threshold = round(score_mean(x = data.frame(threshold, auc, sensi, speci), li = lower.ic, ls = upper.ic)[1], 3)
                   , se.mean = round(score_mean(x = data.frame(threshold, auc, sensi, speci), li = lower.ic, ls = upper.ic)[2],3)
                   , es.mean = round(score_mean(x = data.frame(threshold, auc, sensi, speci), li = lower.ic, ls = upper.ic)[3] ,3)
  )  %>%
  dplyr::mutate(pnt = as.character(pnt)) %>% 
  dplyr::ungroup()

cat("Classifiying gap scores \n")
to_ret <- lapply(unique(means_all$gap_score), function(o){
  
  gp_s <- gsub("score.", "", o)
  
  tmp <- means_all %>% 
    dplyr::filter(gap_score == o) %>% 
    dplyr::mutate(gap_score = gsub("score.", "", gap_score))
  
  
  tmp <- tmp %>% 
    rbind( c(pnt = "Mean", gap_score= "Overall mean", colMeans(tmp[, sapply(tmp, is.numeric)], na.rm = T)))
  
  thr<- as.numeric(tmp$threshold[length(tmp$threshold)])
  rs <- terra::rast(paste0(paths$gap_outDir,"/", gp_s, ".tif"))
  
  print(thr > max(rs[], na.rm = T))
  
  if(thr > max(rs[], na.rm = T)){
    thr <- min(as.numeric(tmp$threshold), na.rm = T)
  }
  
  rs[rs <= thr]  <- NA
  rs[!is.na(rs)] <- 1
  
  terra::writeRaster(rs,paste0(paths$gap_outDir,"/", o, ".tif") ,overwrite = T)
  
  
  return(list (summ_mtrs = tmp, gap_map = rs))
})


names(to_ret) <- gsub("score.", "", unique(means_all$gap_score))

cat("Done \n")

cat("Creating gap map final \n")
resources$final_gap_map <- lapply(to_ret, function(rs){
  rs$gap_map
  
}) %>% 
  terra::rast(.) %>% 
  sum(., na.rm = T) 
terra::writeRaster(resources$final_gap_map, paste0(paths$sp_results, "/", paths$occName,"_final_gap_map.tif"), overwrite = T)


cat("Done")



rs_gaps <- resources$final_gap_map #raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/mesoamerican/americas/gap_models/gap_class_final_new.tif")
resources$covg_ul <- 100 - length(which(rs_gaps[]==3))/length(which(!is.na(rs_gaps[])))*100
resources$covg_ll <- 100 - length(which(rs_gaps[]>=2))/length(which(!is.na(rs_gaps[])))*100


resources$over_all_auc <- lapply(to_ret, function(ls){
  tmp <- ls$summ_mtrs
  tmp$auc.mean[length(tmp$auc.mean)]
}) %>% 
  unlist() %>%
  as.numeric() %>% 
  mean()


auc_bounds <- lapply(to_ret, function(ls){
  tmp <- ls$summ_mtrs
  lw <- tmp$lower.ic[length(tmp$lower.ic)]
  up <- tmp$upper.ic[length(tmp$upper.ic)]
  
  
  return(data.frame(up, lw))
}) %>% 
  dplyr::bind_rows() %>%
  dplyr::mutate(across(everything(.), as.numeric)) %>% 
  colMeans()

sm_mtrs <- lapply(to_ret, function(ls){
  ls$summ_mtrs
})

sm_mtrs$overall_simulation_metrics <- data.frame(over_all_auc         = resources$over_all_auc*100,
                                                 auc_upper_bound      = auc_bounds[1] *100,
                                                 auc_lower_bound      = auc_bounds[2] *100,
                                                 coverage_upper_bound = resources$covg_ul,
                                                 coverage_lower_bound = resources$covg_ll) %>% 
  dplyr::mutate(across(everything(.), round, 2)) %>%
  t(.) %>%
  as_tibble(., rownames = "Metric") %>% 
  dplyr::rename("Value" = up)



resources$overall_simulation_metrics <- sm_mtrs$overall_simulation_metrics

writexl::write_xlsx(sm_mtrs, paste0(paths$gap_valDir,"/simulation_results.xlsx"))
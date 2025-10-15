compute_se_es <- function(thr, obs, pred_bin) {
  cm <- table(factor(pred_bin, levels = c(0,1)),
              factor(obs, levels = c(0,1))) 
  # cm forma:       obs=0 obs=1
  # pred_bin=0      TN     FN
  # pred_bin=1      FP     TP
  TN <- cm[1,1]
  FN <- cm[1,2]
  FP <- cm[2,1]
  TP <- cm[2,2]
  sens <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) 0 else TN / (TN + FP)
  mcc  <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  return(data.frame(sensi = sens, speci = spec, mcc = mcc, threshold = thr))
}


get_sensi_especi_df <- function(obs, preds){
  
  TSS_inner_mat <- data.frame( sensi = numeric(), speci = numeric(), mcc = numeric(), threshold = numeric())
  threshlods <-  seq(0, 1, by = 0.01)
  
  for(i in seq_along(threshlods)){
    
    thr = threshlods[i]
    pred_bin <- ifelse(preds >= thr, 1, 0)
    TSS_inner_mat[i, ] <- compute_se_es(thr, obs, pred_bin)
    
  }
  
  return(TSS_inner_mat)
}





sdm_maxent_approach_function <- function(occName      = occName,
                                         spData       = spData,
                                         bg_data      = bg_data,
                                         var_names    = var_names,
                                         model_outDir = model_outDir,
                                         replic_path  ,
                                         climDir,
                                         clim_spReg,
                                         nFolds       = 5,
                                         beta         = params_tunned$beta,
                                         feat         = params_tunned$features,
                                         doSDraster   = TRUE,
                                         varImp       = TRUE,
                                         validation   = FALSE)
{
  
 
  #si el numero de occurrencias es bajo- entonces se hace
  
  
  
  spData <- dplyr::bind_rows( spData %>% 
                               dplyr::mutate(Y = 1), 
                             bg_data %>% 
                               dplyr::mutate(Y = 0)) %>% 
    dplyr::rename(!!occName := Y) %>% 
    dplyr::mutate(across(everything(.), as.numeric)) %>% 
    dplyr::select(all_of(occName),  "lon" = Longitude, "lat" = Latitude, all_of(var_names))
  
  ### prueba escalando variables
  
  for(nm in c("Accessibility", "Irrigation", "population_density_2015", "dist_h_set")){
    pos <- which(names(spData) == nm)
    spData[, pos] <- ( spData[, pos] - min( spData[, pos], na.rm = T))/(max( spData[, pos], na.rm = TRUE)- min( spData[, pos], na.rm = T))
    
  }
  
  ### fin prueba escalando
  
  cat("Loading environmental raster files \n")
  clim_vars     <- paste0(var_names, ".tif") %in% list.files(climDir, pattern = ".tif$") 
  generic_vars  <- paste0(var_names, ".tif") %in% list.files(clim_spReg, pattern = ".tif$")
  clim_layer    <- lapply(paste0(climDir, "/", var_names[clim_vars], ".tif"), raster)
  
  if(any(generic_vars)){
    generic_layer <- lapply(paste0(clim_spReg,"/", var_names[generic_vars],".tif"), raster)
    clim_layer    <- raster::stack(c(clim_layer, generic_layer))
  }else{
    clim_layer  <- raster::stack(clim_layer)
  }
  
  
  #prueba escalando variables
  clim_layer <- lapply(1:nlayers(clim_layer), function(i){
    
    if(names(clim_layer[[i]]) %in% c("Accessibility", "Irrigation", "population_density_2015", "dist_h_set")){
      
      r = clim_layer[[i]]
      to_ret <- (r - cellStats(r, stat = "min", na.rm = TRUE))/(cellStats(r, stat = "max", na.rm = TRUE)- raster::cellStats(r, stat = "min", na.rm = TRUE))
      
      
    }else{
      to_ret <- clim_layer[[i]]
    }
    return(to_ret)
    })
  
  clim_layer = raster::stack(clim_layer)
  to_save <<- clim_layer
  #fin prueba escalando variables
  
  cat("Initializing MAXENT model fitting throught cross validation. \n ")
  
  #split data  in nFolds for cross - validation process
  cvfolds <- modelr::crossv_kfold(spData, k= nFolds)
  cat("Number of folds:", nFolds, "\n")
  # Do all sdm process
  sdm_results <- cvfolds %>% dplyr::mutate(.
                                           #train 5 sdm models using Maxent and train data
                                           ,model_train = purrr::map2(.x = train, .y = .id, function(.x, .y){
                                             
                                             cat("Training MAXENT model for fold", .y, ", all presence points added to background \n")
                                             data_train <- as.data.frame(.x)
                                             
                                             feat2text <- function(feat){
                                               all_feats <- c('linear','quadratic','product','threshold','hinge')
                                               out <- c(paste0(feat, '=true'),paste0(base::setdiff(all_feats, feat), '=false'))
                                               return(out)
                                             }
                                             
                                             
                                             fit.maxent <- dismo::maxent(x = clim_layer,
                                                                         p = data_train[which(data_train[,1] == 1),c('lon','lat')],
                                                                         a = data_train[,c('lon','lat')],
                                                                         args = c(feat2text(feat), "togglelayertype=ethnicity"))
                                            
                                             return(fit.maxent)
                                             
                                           })
                                           
                                           #evaluate trained model
                                           , predictions_train = purrr::pmap(list(.x = model_train, .y = .id, .z = train), function(.x, .y, .z){
                                             cat("Predicting train data for fold", .y, "\n")
                                             train <- as.data.frame(.z)
                                             predictions <- dismo::predict(object = .x, train[, 4:ncol(train)], type = "cloglog")
                                             dt <-  data.frame(obs = factor(train[, 1]), pred = predictions)
                                             
                                             return(dt)
                                           })
                                           #calculate auc for trained model
                                           ,AUC_train = purrr::map2(.x = predictions_train, .y = .id, function(.x, .y){
                                             cat("Calculating AUC_train for model", .y,"\n")
                                             croc <- pROC::roc(response = .x$obs, predictor = .x$pred)
                                             
                                             return(as.numeric(croc$auc))
                                           } ) 
                                           #calculate max preformance measures (sensitivity, specificity and Treshold) using train data
                                           ,evaluation_train = purrr::map2(.x = predictions_train, .y = .id, function(.x, .y){
                                             
                                             cat("Calculating optimal threshold for model", .y, "\n")
                                             
                                             # croc <- pROC::roc(response = .x$obs, predictor = .x$pred)
                                             # croc_summ <- data.frame (sensi = croc$sensitivities, speci = croc$specificities, threshold =  croc$thresholds) %>% 
                                             #   round(., 3) %>% 
                                             #   dplyr::mutate(., max.TSS = sensi + speci - 1) %>% 
                                             #   dplyr::mutate(., minROCdist = sqrt((1- sensi)^2 + (speci -1)^2))
                                             # 
                                             
                                             croc_summ <- get_sensi_especi_df(obs = .x$obs,
                                                           preds = .x$pred) %>% 
                                               round(., 3) %>% 
                                               dplyr::mutate(., max.TSS = sensi + speci - 1) %>% 
                                               dplyr::mutate(., minROCdist = sqrt((1- sensi)^2 + (speci -1)^2))
                                             
                                             
                                             
                                             max.tss <- croc_summ %>% dplyr::filter(., max.TSS == max(max.TSS)) %>% 
                                               dplyr::mutate(., method = rep("max(TSS)", nrow(.)))
                                             
                                             minRoc <- croc_summ %>% 
                                               dplyr::filter(., minROCdist == min(minROCdist))%>% 
                                               dplyr::mutate(., method = rep("minROCdist", nrow(.)))
                                             
                                             croc_summ <- rbind(max.tss, minRoc) %>% 
                                               dplyr::filter(., speci == max(speci))  %>% 
                                               dplyr::sample_n(., 1)
                                             
                                             return(croc_summ)
                                           })
                                           #Make predictions using testing data
                                           , predictions_test = purrr::pmap(list(.x = test, .y = model_train, .z = .id), function(.x, .y, .z){
                                             
                                             cat("Using test data to predict model", .z," \n")
                                             test <- as.data.frame(.x)
                                             predictions <- dismo::predict(object = .y, test[, 4:ncol(test)],type = "cloglog")
                                             dt <-  data.frame(obs = factor(test[, 1]), pred = predictions)
                                             
                                             return(dt )
                                           })
                                           #Calculate AUC for testing 
                                           , AUC = map2(.x = predictions_test, .y = .id, function(.x, .y){
                                             cat("Calculating AUC for model", .y,"\n")
                                             croc <- pROC::roc(response = .x$obs, predictor = .x$pred)
                                             
                                             return(as.numeric(croc$auc))
                                           } ) 
                                           #calculate max preformance measures (sensitivity, specificity and Treshold) using max(TSS) criterion
                                           , evaluation_test = pmap(list(.x = evaluation_train, .y = .id, .z = predictions_test), function(.x, .y, .z){
                                             
                                             thr <- .x$threshold
                                             
                                             evaluation <- get_sensi_especi_df(obs = .z$obs,
                                                                              preds = .z$pred) %>% 
                                               round(., 3) %>% 
                                               dplyr::mutate(., max.TSS = sensi + speci - 1) %>% 
                                               dplyr::mutate(., minROCdist = sqrt((1- sensi)^2 + (speci -1)^2))
                                             
                                             
                                             
                                             # a <- .z %>% dplyr::filter(., pred >= thr & obs == 1) %>% nrow()
                                             # b <- .z %>% dplyr::filter(., pred >= thr & obs == 0) %>% nrow()
                                             # c <- .z %>% dplyr::filter(., pred < thr & obs == 1) %>% nrow()
                                             # d <- .z %>% dplyr::filter(., pred < thr & obs == 0) %>% nrow()
                                             # 
                                             # #senitivity and specificity
                                             # se <- a/(a+c)
                                             # es <- d/(b+d)
                                             # #Matthews correlation coefficient
                                             # den <- sqrt(a+b)*sqrt(a+c)*sqrt(d+b)*sqrt(d+c)
                                             # den <- ifelse(den  != 0 ,den, 1 )
                                             # mcc <- (a*d - b*c)/den
                                             # #Likelyhood Ratio +
                                             # lr_ps <- se/(1 - es)
                                             # #Likelihood ratio -
                                             # lr_ne <- (1 - se)/es
                                             # 
                                             # #calculate kappa index
                                             # pr_a <- (a+d)/(a+b+c+d)
                                             # pr_e <- (((a+b)/(a+b+c+d))* ((a+c)/(a+b+c+d))) + ( ((c+d)/(a+b+c+d) )* ((b+d)/(a+b+c+d) )) 
                                             # kappa <- (pr_a - pr_e)/(1 - pr_e) 
                                             # 
                                             
                                             #evaluation <- data.frame(threshold= thr, sensi = se, speci = es, matthews.cor = mcc, LR_pos = lr_ps, LR_neg = lr_ne, kappa_index = kappa)
                                             
                                             return(evaluation)
                                             # cat("Calculating optimal threshold for model", .y, "\n")
                                             # croc <- pROC::roc(response = .x$obs, predictor = .x$pred)
                                             # croc_summ <- data.frame (sensi = croc$sensitivities, speci = croc$specificities, threshold =  croc$thresholds) %>% 
                                             #   round(., 3) %>% 
                                             #   dplyr::mutate(., max.TSS = sensi + speci - 1) %>% 
                                             #   dplyr::mutate(., minROCdist = sqrt((1- sensi)^2 + (speci -1)^2))
                                             # 
                                             # max.tss <- croc_summ %>% dplyr::filter(., max.TSS == max(max.TSS)) %>% 
                                             #   dplyr::mutate(., method = rep("max(TSS)", nrow(.)))
                                             # 
                                             # minRoc <- croc_summ %>% 
                                             #   dplyr::filter(., minROCdist == min(minROCdist))%>% 
                                             #   dplyr::mutate(., method = rep("minROCdist", nrow(.)))
                                             # 
                                             # croc_summ <- rbind(max.tss, minRoc) %>% 
                                             #   dplyr::filter(., speci == max(speci))  %>% 
                                             #   dplyr::sample_n(., 1)
                                             # 
                                             # return(croc_summ)
                                           })
                                           #Calculate nAUC using both train and test data
                                           , nAUC = pmap(list(.x = train, .y = test, .z = .id), function(.x, .y, .z){
                                             cat("calculating AUC from NULL model", .z,"\n")
                                             train_dt <- as.data.frame(.x) %>% dplyr::select(., occName, starts_with("lon"), starts_with("lat")  )
                                             test_dt  <- as.data.frame(.y) %>% dplyr::select(., occName, starts_with("lon"), starts_with("lat")  )
                                             
                                             train_p <- train_dt[which(train_dt[, 1] == 1) , 2:3]
                                             train_a <- train_dt[which(train_dt[, 1] == 0) , 2:3]
                                             
                                             gd <- dismo::geoDist(p = train_p, a = train_a, lonlat=TRUE)
                                             pred <- predict(gd, test_dt[, 2:3])
                                             
                                             nAUC <- pROC::roc(response = test_dt[, 1], predictor = pred)
                                             return(as.numeric(nAUC$auc))
                                           }) 
                                           #Calculate cAUC using the formula cAUC = AUC + 0.5 - max( 0.5, nAUC)
                                           , cAUC = purrr::pmap(list(.x = AUC, .y = nAUC, .z = .id), function(.x, .y, .z){
                                             cat("Calculating AUC correction using NULL model", .z, " \n")
                                             cAUC = .x + 0.5 - max( 0.5, .y)
                                             return(cAUC)
                                           })
                                           #Project rasters using maxnet model for mean, median and sd
                                           , do.projections =  purrr::pmap(list(.x = model_train, .y = .id, .z = evaluation_train) ,function(.x, .y, .z){
                                             
                                             cat(">>> Proyecting MAXENT model", .y,"to a raster object \n")
                                             r <- raster::predict(clim_layer, .x, type = "cloglog", progress='text')
                                             writeRaster(r, paste0(replic_path, "/",occName,"_sdm_cvfold-", .y,".tif"), format="GTiff", overwrite = TRUE)
                                             #thresholding raster 
                                             if(!validation){
                                               r[which(r[] < .z$threshold)] <- NA 
                                             }
                                             
                                             writeRaster(r, paste(replic_path, "/",occName,"_sdm_thresh_cvfold-", .y,".tif",sep=""), format="GTiff", overwrite = TRUE)
                                             return(r)
                                           })
                                           
                                           
                                           
  )#end mutate
  
  
  
  
  #calculate  mean, median and sd raster from replicates 
  prj_stk <- sdm_results %>% dplyr::select(., do.projections) %>% unlist() %>% raster::stack() 
  cat("Calculating mean, median and sd for replicates \n")
  mean(prj_stk) %>% writeRaster(., paste0(model_outDir,"/", occName, "_sdm_mean.tif" ), overwrite = TRUE)
  cat("Mean raster calculated \n")
  raster::calc(prj_stk, fun = function(x) {median(x)}) %>% writeRaster(., paste0(model_outDir,"/", occName, "_sdm_median.tif" ), overwrite = TRUE)
  cat("Median raster calculated \n")
  if(doSDraster){
    raster::calc(prj_stk, fun = function(x) {sd(x)}) %>% writeRaster(., paste0(model_outDir,"/", occName, "_sdm_std.tif" ), overwrite = TRUE)
    cat("Sd raster calculated \n")
  }
  
  
  ######## calculate varImportance ########
  
  if(varImp){
    
    sdm_results <- sdm_results %>%
      dplyr::mutate(., varImp = purrr::pmap(list(.x = train, .y = .id, .z = model_train), function(.x, .y, .z){
      cat("Calculating var importance for model", .y, "... \n")
      
      data_train <- as.data.frame(.x) 
      
      pres_to_bg <- data_train[which(data_train[, 1] == 1), ]
      pres_to_bg[,1] <- rep(0, length(pres_to_bg [, 1]))
      
      p <- c(data_train[, 1], pres_to_bg[,1])#adding all presence points to background
      data <- rbind(data_train[, -c(1,2,3)], pres_to_bg[, -c(1,2,3)])#adding all presence points to background
      
      pred_all_vars <-  data.frame(obs = p, pred = dismo::predict(object = .z , data, type = "cloglog"))
      
      roc_all_vars <- pROC::roc(response = pred_all_vars$obs, predictor = pred_all_vars$pred)
      #lovout -> leave one variable out
      #calculate varImportance through randomization procedure as  the sdm package does (https://onlinelibrary.wiley.com/doi/epdf/10.1111/ecog.01881) page.373
      corTest <- c() 
      aucTest <- c()
      for(i in 1:length(var_names)){
        
        cors <- c()
        aucs <- c()
        for(k in 1:5){
          data_rand <- data 
          #randonmize var i 
          data_rand[, i] <- data_rand[(nrow(data_rand)), i]
          
          pred_rand <- data.frame(obs = p, pred = dismo::predict(object = .z , data_rand, type = "cloglog"))
          
          #varImportance by correlation
          cor_eval <- cor(pred_all_vars$pred, pred_rand$pred)
          if(cor_eval < 0  ) cor_eval <- 0
          cors[k] <- 1 - cor_eval
          
          #varImportance by AUC
          roc_rand <- pROC::roc(response = pred_rand$obs, predictor = pred_rand$pred)
          auc_rand <- (roc_all_vars$auc - roc_rand$auc)*2
          if (auc_rand > 1) auc_rand <- 1
          else if (auc_rand < 0) auc_rand <- 0
          aucs[k] <- auc_rand 
        }
        corTest[i] <- round(mean(cors, na.rm = TRUE), 4)
        aucTest[i] <- round(mean(aucs, na.rm = TRUE), 4)
      }
      
      varImp <- data.frame(vars = var_names, "corTest"= corTest, "aucTest" = aucTest)
      return(varImp)
    })
    )#end mutate
    
  }
  #save all results in an .rds file
  cat("Process Done... Saving results as .rds file in the path", paste0(model_outDir, "/sdm_results.rds"), " \n")
  saveRDS(sdm_results, paste0(model_outDir, "/sdm_results.rds"))
  return("OK")
}
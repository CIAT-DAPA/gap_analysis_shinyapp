###############################################################################################
####Creating arguments to run MaxEnt in the SDM Packages (Modified of Jorge Velasquez Script###
###############################################################################################

CreateMXArgs <- function(calibration, use.maxnet = TRUE){
  mxnt.args <- c("linear")
  if(use.maxnet){
    
    if(!is.null(calibration)){
      best.ind <- which.min(calibration$deltaAICc)
      args <- calibration[best.ind, ]
      features <- args$classes
      betamultiplier <- args$regMult
      
      if(grepl("q", features)){
        mxnt.args <- c(mxnt.args, "quadratic")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("h", features)){
        mxnt.args <- c(mxnt.args, "hinge")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("p", features)){
        mxnt.args <- c(mxnt.args, "product")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("t", features)){
        mxnt.args <- c(mxnt.args, "threshold")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      mxnt.args <- c(mxnt.args, paste0("betamultiplier=", betamultiplier))
    }else{
      mxnt.args <- c(mxnt.args, "quadratic", "hinge", "product", "", "betamultiplier=1.0")
    }
    
    mxnt.args <- mxnt.args[which(mxnt.args != "")]
    
    
  }else{
    if(!is.null(calibration)){
      
      best.ind <- which.min(calibration$deltaAICc)
      features <-as.data.frame(cbind(calibration$linear, calibration$quadratic, calibration$product, calibration$hinge, calibration$threshold))
      names(features) <- c("l", "q", "p", "h", "t")
      features <- features[best.ind,]
      features <- unlist(lapply(1:ncol(features),function(i){
        if(features[,i] == TRUE){
          x <- colnames(features[i])
        } else {
          x <- NULL
        }
        return(x)
      }))
      features <- paste(features, collapse = "")
      betamultiplier <- calibration$regMult[best.ind]
      
      if(grepl("q", features)){
        mxnt.args <- c(mxnt.args, "quadratic")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("h", features)){
        mxnt.args <- c(mxnt.args, "hinge")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("p", features)){
        mxnt.args <- c(mxnt.args, "product")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      if(grepl("t", features)){
        mxnt.args <- c(mxnt.args, "threshold")
      } else {
        mxnt.args <- c(mxnt.args, "")
      }
      mxnt.args <- c(mxnt.args, paste0("betamultiplier=", betamultiplier))
    } else {
      mxnt.args <- c(mxnt.args, "quadratic", "hinge", "product", "", "betamultiplier=1.0")
    }
    
    mxnt.args <- mxnt.args[which(mxnt.args != "")]
    
  }  
  return(mxnt.args)
  
}

letter_to_feat <- function(txt){
  
  to_ret <- if(txt == "l"){
    "linear"
  }else if(txt== "q"){
    "quadratic"
  }else if(txt == "p"){
    "product"
  }else if(txt == "h"){
    "hinge"
  }else if(txt == "t"){
    "threshold"
  }else{
    ""
  }
  
  
  return(to_ret)
  
  
  
}


###############################################################################################
####Calibration function using wright et al., 2014 approach###
###############################################################################################

Calibration_function <- function(spData, bg_data, occName, out_path, ommit, use.maxnet = TRUE){
  cat("Initializing calibration step \n")
  
  if(ommit == F){
    
    # Calibration using MaxEnt instead of Maxnet R package.
    
    if(use.maxnet){ 
      #use maxnet
      if(!file.exists(out_path)){
        
        
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
        
        
        
        spData <- spData[complete.cases(spData),]
        
        
        spData <- dplyr::bind_rows(spData %>% 
                                     dplyr::mutate(Y = 1), 
                                   bg_data %>% 
                                     dplyr::mutate(Y = 0)) %>% 
          dplyr::mutate(across(everything(.), as.numeric)) %>% 
          dplyr::mutate(Y = as.factor(Y)) #%>% dplyr::rename(!!occName := Y) # 
        
        
        for(nm in c("Accessibility", "Irrigation", "population_density_2015", "dist_h_set")){
          pos <- which(names(spData) == nm)
          spData[, pos] <- ( spData[, pos] - min( spData[, pos], na.rm = T))/(max( spData[, pos], na.rm = TRUE)- min( spData[, pos], na.rm = T))
          
        }
        
        #rm(bg_data)
        
        
        cat("Calculating best parameters for maxNet \n")
        cat("This process will take several minutes, please be patient. \n")
        
        
        compute_TSS <- function(thr, obs, pred_bin) {
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
          return(data.frame(prob = thr, se = sens, es = spec, TSS = sens + spec - 1))
        }
        
        get_TSS_df <- function(obs, preds){
          
          TSS_inner_mat <- data.frame(thr = numeric(), se = numeric(), es = numeric(), TSS = numeric())
          threshlods <-  seq(0, 1, by = 0.01)
          
          for(i in seq_along(threshlods)){
            
            thr = threshlods[i]
            pred_bin <- ifelse(preds >= thr, 1, 0)
            TSS_inner_mat[i, ] <- compute_TSS(thr, obs, pred_bin)
            
          }
          
          return(TSS_inner_mat)
        }
        
        custom_mtrs <- function(train) {
          #ajham <<- vars
          
          # train_pred <- rbind(data.frame(Y = 1, pred = vars$occs.train.pred),
          #                     data.frame(Y=0, pred = vars$bg.train.pred))
          
          train_pred <- train
          
          train_roc <- get_TSS_df(obs = train_pred$Y,
                                  preds = train_pred$pred)
          
          
          #test_PRACU = MLmetrics::PRAUC(y_pred = test_pred$pred, y_true = test_pred$Y)
          #test_roc <- pROC::roc(test_pred, response ="Y", predictor = "pred", quiet = T, ret = "all_coords", direction = "<", algorithm = 3)
          
          #test_roc$TSS <- test_roc$sensitivity+test_roc$specificity-1
          
          pos <- which.max(train_roc$TSS)
          out <- data.frame(thr = train_roc$thr[pos],
                            train_TSS = train_roc$TSS[pos], # train_PRAUC,
                            #val_TSS = max(test_roc$TSS), #test_PRACU,
                            row.names = NULL)
          #train_PRAUC <- MLmetrics::PRAUC(y_pred = train_pred$pred, y_true = train_pred$Y)
          
          #train_roc <- pROC::roc(train_pred, response ="Y", predictor = "pred", quiet = T, ret = "all_coords", direction = "<", algorithm = 3)
          #train_roc$TSS <- train_roc$sensitivity+train_roc$specificity-1
          
          #out_mtrs_df <<- append(out_mtrs_df, list(train_roc))
          
          return(out)
        }
        
        
        #spData = sp_list$spData
        
        #bg_data = sp_list$bg_data
        
        
        tune_rec <- recipes::recipe(spData %>% dplyr::select(-Latitude, -Longitude), Y ~ .)
        
        
        wk_models <- workflowsets::workflow_set(preproc = list(default = tune_rec),
                                                models =  list( default_maxent  = tidysdm::sdm_spec_maxent()),
                                                cross = T)
        
        spData <- sf::st_make_valid(spData  %>% 
                                          sf::st_as_sf(., 
                                                       coords = c("Longitude", "Latitude"), 
                                                       crs = "EPSG:4326"))
        
        cv_folds <- spatialsample::spatial_block_cv(spData, v = 5)
        
        
        
        fc_vals <- c("l", "lq", "lqp",  "lqph", "lqpht")
        rm_vals <- seq(1, 5, by = 1)
        
        custom_grid <- tidyr::crossing(
          feature_classes = fc_vals,
          regularization_multiplier = rm_vals
        )
        
        tunning_df <- plyr::llply(1:nrow(custom_grid), function(row){
          #cat("processing: ", row, "\n")
          comb <- custom_grid[row, ]
          id_txt <- paste0("fc.", comb[1], "_rm.", comb[2])
          
          model = tune::finalize_workflow(wk_models$info[[1]]$workflow[[1]], comb)
          
          
          splits_df = lapply(1:nrow(cv_folds),  function(i){
            
            split_i <- cv_folds$splits[[i]]
            
            fit_train <- parsnip::fit(model , 
                                      data = rsample::training( split_i) )
            
            preds_train <- predict(fit_train, rsample::training(split_i), type = "prob") %>%
              dplyr::select(pred = .pred_1) %>% 
              bind_cols(pred = predict(fit_train, rsample::training(split_i))) %>% 
              bind_cols(Y = rsample::training(split_i)$Y)
            
            preds_test <- predict(fit_train, rsample::testing(split_i), type = "prob") %>%
              dplyr::select(pred = .pred_1) %>% 
              bind_cols(pred = predict(fit_train, rsample::testing(split_i))) %>% 
              bind_cols(Y = rsample::testing(split_i)$Y)
            
            #train_TSS <- compute_TSS(thr = 0.5, obs = preds_train$Y, pred_bin = preds_train$.pred_class)$TSS
            
            
            res <- custom_mtrs(preds_train)#data.frame(train_TSS = train_TSS, val_TSS = val_TSS)#
            val_TSS <- compute_TSS(thr = res$thr, obs = preds_test$Y, pred_bin = ifelse(preds_test$pred > res$thr, 1, 0 )  )$TSS
            
            res <- cbind(res, val_TSS = val_TSS, tune.args = id_txt, split_id = i, classes = toupper(comb[1]), regMult = as.numeric(comb[2]))
            
            return(res)
            
          }) %>% dplyr::bind_rows(.)
          
          
          
          
          return(splits_df)
          
        }, .progress =  "text")
        
        calibration_txt <- dplyr::bind_rows(tunning_df) %>%
          dplyr::group_by(tune.args) %>% 
          dplyr::reframe(median_train_TSS = median(train_TSS, na.rm =T),
                         median_val_TSS = median(val_TSS, na.rm = T),
                         avg_train_TSS = mean(train_TSS, na.rm =T),
                         avg_val_TSS = mean(val_TSS, na.rm = T),
                         classes = unique(classes),
                         regMult = unique(regMult) ) %>% 
          dplyr::mutate(diff_median = abs((median_val_TSS - median_train_TSS)/median_train_TSS),
                        diff_avg    = abs((avg_val_TSS - avg_train_TSS)/avg_train_TSS)) 
        
        
        
        write.csv(calibration_txt,  out_path, quote = F, row.names = F)
        
        
      }else{
        cat("Calibration File already created, importing it \n")
        
        calibration_txt <- read.csv(out_path)
      }
      
      
      
    }
    
    args <- calibration_txt[which.min(calibration_txt$diff_median),] #se elige el modelo que menos sobre ajuste tenga
    
    feat <- sapply(tolower(unlist(base::strsplit(args$classes, ""))), letter_to_feat, simplify = T, USE.NAMES = F)
    beta <- unlist(args$regMult)
    
    # 
    # feat <- CreateMXArgs(calibration_txt, use.maxnet)
    # beta <- feat[(grepl("betamultiplier=", feat))]
    # beta <- as.numeric(gsub("betamultiplier=", "", beta))
    # feat <- feat[(!grepl("betamultiplier=", feat))]
    
    
  } else {
    cat("Ommiting calibration step\n")
    calibration_txt <- NULL
    cat("Setting default parameters","\n")
    beta <- 1
    feat <- c("linear", "quadratic", "hinge", "product")
    
  }
  
  
  return(list(features = feat, beta = beta))
  cat("Process done... \n")
}

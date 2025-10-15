pacman::p_load(biomod2, raster, terra, tidysdm, parsnip, tidyverse, sf, tune, future, furrr, ENMeval, SDMtune, MLmetrics)

library(terra)
# Load species occurrences (6 species available)
DataSpecies = read.csv("~/LGA_workshop/pumking/africa_custom/input_data/occurrences_files/pumking_original_data_in.csv")
DataSpecies = DataSpecies %>% 
  tidyr::drop_na() %>% 
  dplyr::select(Y,  X_WGS84=Latitude, Y_WGS84=Longitude) %>% 
  dplyr::mutate(Y = 1)
# Select the name of the studied species
myRespName <- 'Y'
# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])
# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data(bioclim_current)
myExpl <- terra::rast(list.files("~/LGA_workshop/black_fonio/africa/input_data/generic_rasters", pattern = ".tif$", full.names = T))

myResp.PA <-  myResp
myResp.PA.vect <- vect(cbind(myRespXY, myResp.PA) , geom = c("X_WGS84","Y_WGS84"))

PA.r <- bm_PseudoAbsences(resp.var = myResp.PA.vect,
                          expl.var = myExpl,
                          nb.rep = 1,
                          nb.absences = 1000,
                          strategy = 'random')
aa = PA.r$env
aa$pseudo = unlist(PA.r$pa.tab)
aa$source = "random"
aa$Y= PA.r$sp

plot(myExpl[[1]]);points(PA.r$xy[,1],  PA.r$xy[,2])

PA.s <- bm_PseudoAbsences(resp.var = myResp.PA.vect,
                          expl.var = myExpl,
                          nb.rep = 1,
                          nb.absences = 1000,
                          strategy = 'sre',
                          sre.quant = 0.025)
bb = PA.s$env
bb$pseudo = unlist(PA.s$pa.tab)
bb$source = "sre"
bb$Y = PA.s$sp


pseudo_max <- read.csv("~/LGA_workshop/pumking/africa_custom/results/curcubita_maxima/species_distribution/pseudo_absences/pseudo_absences_curcubita_maxima.csv") %>% 
  dplyr::mutate(pseudo = T, source = "lga_random", Y = 0, status = "pseudo") %>% 
  dplyr::select(-Latitude, -Longitude)


ff = dplyr::bind_rows(aa,bb) %>%
  dplyr::mutate(Y= ifelse(is.na(Y), 0, Y)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(status = ifelse(Y != 1 | is.na(Y), "pseudo", "presc"))

x11()
ff %>% 
  dplyr::bind_rows(., pseudo_max) %>% 
  dplyr::select(-Y ) %>% 
  tidyr::pivot_longer(cols = Accessibility:wind_speed , names_to = "labs", values_to = "vals") %>% 
  ggplot(aes(x = source, y= vals)) + 
  geom_boxplot(aes(fill = status))+
    facet_wrap(vars(labs), scales = "free")

############################################
####### tidysdm tunning and modelling #####
##########################################


swd <- read.csv("~/LGA_workshop/black_fonio/africa/results/black_fonio/species_distribution/pseudo_absences/pseudo_absences_black_fonio.csv") %>% 
  dplyr::mutate(Y = 0)
occ <- read.csv("~/LGA_workshop/black_fonio/africa/results/black_fonio/black_fonio_occurrences.csv") %>% 
  dplyr::mutate(Y = 1) %>% 
  dplyr::select(-status)


df <- bind_rows(swd, occ) %>% 
  dplyr::mutate(Y = as.factor(Y))


maxent_spec <- tidysdm::maxent(
  regularization_multiplier = tune::tune(),
  feature_classes = tune::tune()
)

cv <- rsample::vfold_cv(df %>% 
                          dplyr::select(-Latitude, -Longitude), v = 5, strata = "Y")

lapply(cv$splits, function(spl){
  train <- spl %>% rsample::analysis() %>%  pull(Y) %>% table %>% prop.table() %>% matrix(., nrow =1, ncol = 2, byrow = T) %>% as.data.frame() %>% mutate(lab = "train")
  test  <- spl %>% rsample::assessment() %>% pull(Y) %>% table %>% prop.table() %>% matrix(., nrow =1, ncol = 2, byrow = T) %>% as.data.frame() %>% mutate(lab = "test")
  
  bind_rows(train, test)
})

nrow(df)

system.time({
  #cv <- spatialsample::spatial_block_cv(sf::st_as_sf(df, coords = c("Longitude", "Latitude"), crs = st_crs(4326)), v = 2)
  maxent_tune_res <- maxent_spec %>%
    tune::tune_grid(Y ~ ., cv, grid = 15, control = control_grid(allow_par  = T,
                                                                parallel_over = "everything"))
  
})

#tune::collect_metrics(maxent_tune_res)
calibration <- tune::show_best(maxent_tune_res, metric = "roc_auc", n = 10) %>% 
  dplyr::select(regMult = regularization_multiplier,  
                classes = feature_classes,
                roc_auc = mean,
                metric  = .metric,
                std_err )

write.csv(calibration, "D:/OneDrive - CGIAR/Documents/tunned_sdm_parms_new.csv", row.names = F)


args <- calibration[which.max(calibration$roc_auc),]

feat <- unlist(base::strsplit(args$classes, ""))

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


params_tunned_new <- list()
params_tunned_new$features <- sapply(feat, letter_to_feat, simplify = T, USE.NAMES = F)
params_tunned_new$beta     <- unlist(args$regMult)

#### LGA

system.time({
params_tunned  <- Calibration_function(spData   = occ,
                                       bg_data  = swd,
                                       occName  = "Y",
                                       out_path = "D:/OneDrive - CGIAR/Documents/tunned_sdm_parms.csv",
                                       ommit = F,
                                       use.maxnet = TRUE)


})


#### sdm

var_names <- read.csv("~/LGA_workshop/black_fonio/africa/results/black_fonio/species_distribution/sdm_variables_selected.csv")%>% 
  dplyr::pull(x)


##### new
sdm_maxent_approach_function(occName      = "bfonio",
                             spData       =  occ %>% dplyr::select(-any_of("status")),
                             bg_data      = swd,
                             var_names    = var_names,
                             model_outDir = "~/sdm_tests/new_tune",
                             replic_path  = "~/sdm_tests/new_tune/sdm_cv_folds",
                             nFolds       = 5,
                             climDir      = "~/LGA_workshop/black_fonio/africa/input_data/generic_rasters",
                             clim_spReg   = "",
                             beta         = params_tunned_new$beta,
                             feat         = params_tunned_new$features,
                             doSDraster   = FALSE,
                             varImp       = FALSE,
                             validation   = F)

##### OLD
sdm_maxent_approach_function(occName      = "bfonio",
                             spData       =  occ %>% dplyr::select(-any_of("status")),
                             bg_data      = swd,
                             var_names    = var_names,
                             model_outDir = "~/sdm_tests/old_tune",
                             replic_path  = "~/sdm_tests/old_tune/sdm_cv_folds",
                             nFolds       = 5,
                             climDir      = "~/LGA_workshop/black_fonio/africa/input_data/generic_rasters",
                             clim_spReg   = "",
                             beta         = params_tunned$beta,
                             feat         = params_tunned$features,
                             doSDraster   = FALSE,
                             varImp       = FALSE,
                             validation   = F)





########################################################
######## NEW METHOD OF PSEUDO #########################
######################################################

myExpl <- terra::rast(list.files("~/LGA_workshop/black_fonio/africa/input_data/generic_rasters", pattern = ".tif$", full.names = T))

myResp.PA <-  myResp
myResp.PA.vect <- vect(cbind(myRespXY, myResp.PA) , geom = c("X_WGS84","Y_WGS84"))

DataSpecies = occ %>% dplyr::select(-any_of("status"))
DataSpecies = DataSpecies %>% 
  tidyr::drop_na() %>% 
  dplyr::select(Y,  X_WGS84=Latitude, Y_WGS84=Longitude) %>% 
  dplyr::mutate(Y = 1)
# Select the name of the studied species
myRespName <- 'Y'
# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])
# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

myResp.PA <-  myResp
myResp.PA.vect <- vect(cbind(myRespXY, myResp.PA) , geom = c("X_WGS84","Y_WGS84"))


PA.s <- bm_PseudoAbsences(resp.var = myResp.PA.vect,
                          expl.var = myExpl,
                          nb.rep = 1,
                          nb.absences = 220,
                          strategy = 'sre',
                          sre.quant = 0.025)
bb = PA.s$env
bb$Y = PA.s$sp
bb$Longitude = PA.s$xy[, 1]
bb$Latitude  = PA.s$xy[, 2]



swd_new = bb %>%
  dplyr::mutate(Y= ifelse(is.na(Y), 0, Y)) %>% 
  dplyr::filter(Y == 0) %>% 
  tidyr::drop_na() %>% 
  dplyr::select(Y, everything(.))


df <- bind_rows(swd_new, occ) %>% 
  dplyr::mutate(Y = as.factor(Y))


maxent_spec <- tidysdm::maxent(
  regularization_multiplier = tune::tune(),
  feature_classes = tune::tune()
)

cv <- rsample::vfold_cv(df %>% 
                          dplyr::select(-Latitude, -Longitude), v = 5, strata = "Y")


system.time({
  #cv <- spatialsample::spatial_block_cv(sf::st_as_sf(df, coords = c("Longitude", "Latitude"), crs = st_crs(4326)), v = 2)
  maxent_tune_res <- maxent_spec %>%
    tune::tune_grid(Y ~ ., cv, grid = 15, control = control_grid(allow_par  = T,
                                                                 parallel_over = "everything"))
  
})

#tune::collect_metrics(maxent_tune_res)
calibration <- tune::show_best(maxent_tune_res, metric = "roc_auc", n = 10) %>% 
  dplyr::select(regMult = regularization_multiplier,  
                classes = feature_classes,
                roc_auc = mean,
                metric  = .metric,
                std_err )

write.csv(calibration, "~/sdm_tests/new_tune/tunned_sdm_parms_new.csv", row.names = F)


args <- calibration[which.max(calibration$roc_auc),]

feat <- unlist(base::strsplit(args$classes, ""))

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


params_tunned_new <- list()
params_tunned_new$features <- sapply(feat, letter_to_feat, simplify = T, USE.NAMES = F)
params_tunned_new$beta     <- unlist(args$regMult)



##### new
sdm_maxent_approach_function(occName      = "bfonio",
                             spData       =  occ %>% dplyr::select(-any_of("status")),
                             bg_data      = swd_new,
                             var_names    = var_names,
                             model_outDir = "~/sdm_tests/new_tune",
                             replic_path  = "~/sdm_tests/new_tune/sdm_cv_folds",
                             nFolds       = 5,
                             climDir      = "~/LGA_workshop/black_fonio/africa/input_data/generic_rasters",
                             clim_spReg   = "",
                             beta         = params_tunned_new$beta,
                             feat         = params_tunned_new$features,
                             doSDraster   = FALSE,
                             varImp       = FALSE,
                             validation   = F)


system.time({
  params_tunned  <- Calibration_function(spData   = occ,
                                         bg_data  = swd_new,
                                         occName  = "Y",
                                         out_path = "~/sdm_tests/new_tune/tunned_sdm_parms.csv",
                                         ommit = F,
                                         use.maxnet = TRUE)
  
  
})
 ##### OLD
sdm_maxent_approach_function(occName      = "bfonio",
                             spData       =  occ %>% dplyr::select(-any_of("status")),
                             bg_data      = swd_new,
                             var_names    = var_names,
                             model_outDir = "~/sdm_tests/old_tune",
                             replic_path  = "~/sdm_tests/old_tune/sdm_cv_folds",
                             nFolds       = 5,
                             climDir      = "~/LGA_workshop/black_fonio/africa/input_data/generic_rasters",
                             clim_spReg   = "",
                             beta         = params_tunned$beta,
                             feat         = params_tunned$features,
                             doSDraster   = FALSE,
                             varImp       = FALSE,
                             validation   = F)

##################################################
############ Usando ENMeval #####################
################################################

varstack <- terra::rast(list.files("~/LGA_workshop/black_fonio_new/custom/input_data/generic_rasters", pattern = ".tif$", full.names = T))

occ <- read.csv("~/LGA_workshop/black_fonio/africa/results/black_fonio/black_fonio_occurrences.csv") %>% 
  #dplyr::mutate(Y = 1) %>% 
  dplyr::select(-status)


pseudo_occ <- function(xy, varstack, tms){
  
  ## funcion para seleccionar pseudo-PRESENCIAS cuando el numero de occ es bajo,
  ## se seleccionar a partir de un SVM que indica pixeles con condiciones muy similares a las occ
  
  
  ############ OCVMprofiling ####################
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
  
  #varstack <- terra::mask(varstack, svm_msk)
  
  
  return(to_ret)
 
}

occ <- pseudo_occ(xy = occ, varstack = varstack, tms  = 10)



pseudo <-pseudoAbsences_generator(data        = occ,
                                  climDir     = "~/LGA_workshop/black_fonio_new/custom/input_data/generic_rasters", 
                                  aux_dir     = "~/LGA_workshop/black_fonio_new/custom/input_data/auxiliary_rasters",
                                  clsModel    = "Y", 
                                  overwrite   = F, 
                                  occName     = "black_fonio",
                                  correlation = 3, 
                                  pa_method   = "ecoreg",
                                  bg_out_path = "D:/OneDrive - CGIAR/Documents/prueba_1/pseudo.csv",
                                  mask_path   = "~/LGA_workshop/black_fonio_new/custom/input_data/mask_custom.tif",
                                  smd_var_selected_path = "D:/OneDrive - CGIAR/Documents/prueba_1/varselected.csv",
                                  ecoreg_path = "www/masks/World_ELU_2015_5km.tif")


tune.args = list(fc = c("L",    "P",    "Q",    "H",   "LP",   "LQ" ,  "LH"  ,  "PQ" , 
                        "PH",  "QH", "LPQ",  "LPH",   "LQH",   "PQH", "LPQH"), rm = seq(5, 15, by =2 ))

feats <-   unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
tune.args <- list(fc =  toupper(c("l","lq", "lh", "lqp", "lqhp", "lqhpt")), rm = seq(1, 5, by = 1 ))

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

out_mtrs_df <- list()

custom_mtrs <- function(vars) {
  ajham <<- vars
  
  train_pred <- rbind(data.frame(Y = 1, pred = vars$occs.train.pred),
                      data.frame(Y=0, pred = vars$bg.train.pred))
  
  
  train_roc <- get_TSS_df(obs = train_pred$Y,
             preds = train_pred$pred)
  
  
  #train_PRAUC <- MLmetrics::PRAUC(y_pred = train_pred$pred, y_true = train_pred$Y)
  
  #train_roc <- pROC::roc(train_pred, response ="Y", predictor = "pred", quiet = T, ret = "all_coords", direction = "<", algorithm = 3)
  #train_roc$TSS <- train_roc$sensitivity+train_roc$specificity-1
  
  
  test_pred <- rbind(data.frame(Y =1, pred = vars$occs.val.pred),
                     data.frame(Y =0, pred = vars$bg.val.pred))
  
  
  #test_PRACU = MLmetrics::PRAUC(y_pred = test_pred$pred, y_true = test_pred$Y)
  #test_roc <- pROC::roc(test_pred, response ="Y", predictor = "pred", quiet = T, ret = "all_coords", direction = "<", algorithm = 3)
  
  #test_roc$TSS <- test_roc$sensitivity+test_roc$specificity-1
  
  test_roc <- get_TSS_df(obs = test_pred$Y,
                        preds = test_pred$pred)
  
  out <- data.frame(train_TSS = max(train_roc$TSS), # train_PRAUC,
                    val_TSS = max(test_roc$TSS), #test_PRACU,
                    row.names = NULL)
  
  out_mtrs_df <<- append(out_mtrs_df, list(train_roc))
  
  return(out)
}

system.time({
  eval_res <- ENMeval::ENMevaluate(occs = occ[, c("Longitude", "Latitude")], 
                                   envs = varstack, 
                                   bg = pseudo[, c("Longitude", "Latitude")],
                                   algorithm = 'maxnet', 
                                   partitions = 'block',
                                   partition.settings = list(
                                     orientation = "lat_lon", 
                                     kfolds = 5)
                                   ,tune.args = tune.args,
                                   raster.preds = F, user.eval = custom_mtrs, doClamp = T)
  
})


eval_res@results



eval_res@results.partitions %>% 
  group_by(tune.args) %>% 
  dplyr::reframe(median_val_auc = median(auc.val, na.rm = T),
                 median_train_TSS = median(train_TSS, na.rm =T),
                 median_val_TSS = median(val_TSS, na.rm = T),
                 avg_val_auc = mean(auc.val, na.rm = T),
                 avg_train_TSS = mean(train_TSS, na.rm =T),
                 avg_val_TSS = mean(val_TSS, na.rm = T)) %>% 
  as.data.frame() %>% 
  dplyr::mutate(diff_median = abs((median_val_TSS - median_train_TSS)/median_train_TSS),
                diff_avg = abs((avg_val_TSS - avg_train_TSS)/avg_train_TSS)) %>%  View
  pull(diff_median) %>% which.min()


beta = 1
feats = c("linear", "hinge")

mod.seq <- eval.models(eval_res)$fc.LQ_rm.5
mod.seq$betas
x11();plot(mod.seq, type = "cloglog")
vr <- "aridityIndexThornthwaite" 
mod.complex.mrc = maxnet::response.plot(mod.seq, 
                      v = vr, 
                      type = "cloglog",
                      plot = FALSE)
mod.complex.mrc %>% 
  dplyr::rename( var := vr) %>% 
ggplot(., aes(x = var, y = pred)) + 
  geom_line() +
  labs(title = "Response curve")+
  xlab(vr)+
  ylab("cloglog prediction") +
  theme_bw()

eval_res@results.partitions %>% 
  ggplot()+
  geom_line(aes(x = fold, y = train_TSS), colour = "blue")+
  geom_line(aes(x = fold, y = val_TSS), colour = "red")+
  facet_wrap(~tune.args)+
  theme_bw()


var_names = read.csv("~/LGA_workshop/black_fonio_new/custom/results/black_fonio/species_distribution/sdm_variables_selected.csv", stringsAsFactors = F) %>% 
  dplyr::pull(x)


occ <- dplyr::bind_cols(occ[, 1:3], apply(occ[,-c(1:3)], 2, function(col){
  (col - min(col, na.rm = T)[1])/(max(col, na.rm = T)[1]- min(col, na.rm = T)[1])
  
}) )

pseudo <- dplyr::bind_cols(pseudo[, 1:3], apply(pseudo[,-c(1:3)], 2, function(col){
  (col - min(col, na.rm = T)[1])/(max(col, na.rm = T)[1]- min(col, na.rm = T)[1])
  
}) )

sdm_maxent_approach_function(occName      = "black_fonio",
                             spData       =  occ %>% dplyr::select(-any_of("status")),
                             bg_data      = pseudo,
                             var_names    = var_names,
                             model_outDir = "~/prueba_1/sdm_tests",
                             replic_path  = "~/prueba_1/sdm_tests/sdm_cv_folds",
                             nFolds       = 5,
                             climDir      = "~/LGA_workshop/black_fonio_new/custom/input_data/generic_rasters",
                             clim_spReg   = "",
                             beta         = beta,
                             feat         = feats,
                             doSDraster   = FALSE,
                             varImp       = FALSE,
                             validation   = F)



rds <- readRDS("~/prueba_1/sdm_tests/sdm_results.rds")

to_save[[1]]

rds$test
rds$predictions_train[[1]] %>% summary()

rds$evaluation_train

sapply(rds$evaluation_test, function(df)(df[which.max(df$max.TSS), ]))







ps_cm <- read.csv("D:/OneDrive - CGIAR/Documents/prueba_1/pseudo.csv")
ps_af <- read.csv("D:/OneDrive - CGIAR/Documents/prueba_1/pseudo_africamsk.csv")


gg = ps_cm %>% 
  dplyr::mutate(source  = "cm") %>% 
  dplyr::bind_rows(., ps_af %>% dplyr::mutate(source = "af") %>% dplyr::select(-Y)) %>% 
  dplyr::bind_rows(., occ %>% dplyr::mutate(source = "occ") %>% dplyr::select(-Y)) %>% 
  dplyr::select(-Y) %>% 
  tidyr::pivot_longer(., -source, names_to = "var", values_to = "val") %>% 
  dplyr::filter(!var %in% c("Latitude", "Longitude")) %>% 
  ggplot(aes(x = source, y = val, fill = source))+
  geom_boxplot() +
  facet_wrap(~var, scales = "free")+
  theme_bw()

  ggsave(plot = gg, "D:/OneDrive - CGIAR/Documents/prueba_1/bx1.png", dpi = 300, width = 12, height = 8)

  
###########################################################
#### USANDO ENSAMBLE DE MODELOS DE LA CV #################
#########################################################
# ADAPTADO DE https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.70019
  





#####################################################
######### USANDO RANDOMFOREST ######################
###################################################
  
  
  
  occ_grid = eval_res@occs.grp
  bg_grid = eval_res@bg.grp
  
  no.forest <- 25
  no.trees <- 100
  nVars <- 8
  
  for(level in levels(occ_grid)){
    
    occ_id_test  <- which(occ_grid %in% level)
    occ_id_train <- which(!occ_grid %in% level)
    
    bg_id_test   <- which(bg_grid %in% level)
    bg_id_train  <- which(!bg_grid %in% level)
    
    
    bg_train <- pseudo %>% 
      dplyr::slice(bg_id_train) %>% 
      dplyr::mutate(Y = 0)
    
    bg_test <- pseudo %>% 
      dplyr::slice(bg_id_test) %>% 
      dplyr::mutate(Y =0)
    
    occ_train <- occ %>% 
      dplyr::slice(occ_id_train) %>% 
      dplyr::mutate(Y = 1)
    occ_test <- occ %>% 
      dplyr::slice(occ_id_test) %>% 
      dplyr::mutate(Y = 1)
    
    train_df <- dplyr::bind_rows(occ_train, bg_train) %>% 
      dplyr::sample_n(., size = 220) %>% 
      #dplyr::mutate(Y = factor(Y, levels = c(0, 1))) %>% 
      dplyr::select(-"Latitude", -"Longitude")
    
    test_df <- dplyr::bind_rows(occ_test, bg_test) %>% 
      #dplyr::mutate(Y = factor(Y, levels = c(0, 1)))%>% 
      dplyr::select(-"Latitude", -"Longitude")
    
    vrs <- names(occ %>% dplyr::select(-"Y", -"Latitude", -"Longitude"))
    
    model1 <- as.formula(paste('factor(Y) ~', paste(paste(vrs), collapse = '+', sep =' ')))
    
   
    #rfmodel <- randomForest::randomForest(model1, data = train_df, ntree = 500, na.action = na.omit, nodesize = 2)
    pesos <- c("0" = 1, "1" = 10000)
    
    rfmodel <- randomForest::tuneRF(x = train_df %>% dplyr::select(-Y),
                                   y = factor(train_df$Y),
                                   ntreeTry = 100,
                                   stepFactor = 0.5,
                                   doBest = T,
                                   classwt = pesos )
    
    pred_train <- predict(rfmodel, train_df %>% dplyr::select(-Y), type = "prob")[,2]
    pred_test  <- predict(rfmodel, test_df %>% dplyr::select(-Y), type = "prob")[, 2]
    
    
    prAUC_train  <- MLmetrics::PRAUC(y_pred = pred_train, y_true = train_df$Y)
    
    prAUC_test <- MLmetrics::PRAUC(y_pred = pred_test, y_true = test_df$Y)
    
    
    auc_train <- dplyr::bind_cols(Y = train_df$Y, pred = pred_train)
    
    
    
    res_roc <- pROC::roc(auc_train, response ="Y", predictor = "pred", quiet = T, ret = "all_coords", direction = "<", algorithm = 3)
    
    res_roc$TSS = res_roc$sensitivity + res_roc$specificity -1
    
    thr <- as.numeric(quantile(env_pred, probs = 0.7))#res_roc$threshold[which.max(res_roc$TSS)]
    
    
    hist(env_pred)
    abline(v = thr, col = "red")
    
    env_df = terra::as.data.frame(varstack, na.rm = F)
    non_na <- which(rowSums(is.na(env_df)) == 0)
    env_pred = predict(rfmodel, env_df[non_na,], type = "prob")[,2]
    
    full_pred <- rep(NA, nrow(env_df))
    full_pred[non_na] <- env_pred
    
    ref <- varstack[[1]]
    values(ref) <- full_pred
    ref_thr <- terra::ifel(ref <= thr, NA, ref)
    
    nrow(ref_thr[!is.na(ref_thr)])
    
    predicted <- as.numeric(predict(rfmodel, envtest))
    observed <- as.vector(envtest[,'pb'])
    auc[[repe]] <- auc(observed, predicted) 
    
    
    
  }
  
  
  model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  
  
  
  





















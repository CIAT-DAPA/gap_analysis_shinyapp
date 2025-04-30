pacman::p_load(biomod2, terra, tidysdm, parsnip, tidyverse, sf, tune, future, furrr)

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








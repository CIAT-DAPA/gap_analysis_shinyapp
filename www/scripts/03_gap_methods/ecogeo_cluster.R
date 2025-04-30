
#SCRIPT CREATED BY: ANDRES CAMILO MENDEZ
#THIS SCRIPT COMES WITH ABSOLUTELY NO WARRANTY
#FEEL FREE TO USE IN THE GAP ANALISYS FOR LANDRACES PROJECT

#THIS SCRIPT ALLOW YOU TO CONSTRUCT A ECO-GEOGRAPHIC CLUSTER USING MAHALANOBIS DISTANCES AND WARD METHOD TO CONSOLIDATE THE CLUSTERS


                #### FUNCTION TO CREATE ECOGEOGRAFICAL CLUSTER##############
ecogeo_clustering <- function(n.sample = 10000, 
                              var_names, 
                              k.clust = 11,
                              sdm_path,
                              climDir = climDir,
                              gap_outDir = gap_outDir){
  

  # Analysis region: "americas", "world"
  
  
  
  cat("Importing SDM raster \n   \n \n")

  SDM <- raster(sdm_path) 
  
  #selec rasters
  
  # sdm_obj <- read.sdm(paste(sp_Dir,"/sdm.sdm",sep=""))
  # var_names <- names(sdm_obj@data@features)[2:ncol(sdm_obj@data@features)]
  # 
  #var_names <- paste0(var_names,".tif")

  # vars <- c("Accessibility"        ,    "Altitude"       ,          "aridityIndexThornthwaite" ,"bio_14"  ,                
  #  "bio_15"  ,                "bio_18"    ,               "bio_19"    ,               "bio_2"   ,                
  # "climaticMoistureIndex" ,   "continentality"     ,      "dist_rivers"         ,     "embergerQ" ,              
  # "Irrigation"          ,    "minTempWarmest"        ,   "monthCountByTemp10"     ,  "PETDriestQuarter" ,       
  # "PETWarmestQuarter"      ,  "PETWettestQuarter"     ,   "Physical.area")
  
  pos <- which(list.files(path = climDir ) %in% paste0(var_names, ".tif"))
  
  path <- list.files(path = climDir, full.names=TRUE)[pos]
  cat( "Importing WorldClim and Envirem rasters \n   \n  \n"  )
  cat( "Take care of your RAM  \n"  )
  
  environ_df <- terra::rast(path) %>% 
    terra::mask(., mask = terra::rast(SDM)) %>% 
    terra::as.data.frame(., xy = T, na.rm = F) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(., ID = 1:nrow(.)) %>%
    dplyr::select(., ncol(.) , 1:(ncol(.)-1) )
  
 
  rm(pos, var_names, environ_list)


rownames(environ_df) <- environ_df$ID

zeroVar <- caret::nearZeroVar(environ_df)
if(length(zeroVar) != 0){
  environ_df <- environ_df[, -zeroVar]
}

#environ <- environ[complete.cases(environ),]

#environ_scaled <- scale(environ_df_in, center = T, scale = T)

####
cat( paste("Starting clustering process whit: ", n.sample," Sample size \n   \n \n")  )

set.seed(100)

if(n.sample >= nrow(environ_df)){
  
  df_temp <- environ_df
  rownames(df_temp) <- rownames(environ_df)
 
  
}else{
  
  muestra <- sample(1:nrow(environ_df), n.sample)
  df_temp <- environ_df[muestra,]
  rownames(df_temp) <- rownames(environ_df[muestra, ])
  
  
}

cat( "Calculating Mahalanobis distances... \n   \n \n")



mahaRed_dist <-  distances::distances(as.matrix(df_temp[, 4:ncol(df_temp)]), normalize = "mahalanobize")



cat("Hierarchical Clustering to distances using the WARD method \n   \n \n")

clust_hc <- fastcluster::hclust(distances::distance_matrix(mahaRed_dist), method = "ward.D") 

rm(mahaRed_dist); g <- gc(); rm(g)

memb <- cutree(clust_hc, k= k.clust)

rm(clust_hc); g <- gc(); rm(g)

environR_clust <- data.frame( df_temp, clust= factor(memb) )


if(n.sample < nrow(environ_df)){
  
  cat("Starting assignation of occurrences to each cluster using RF \n   \n \n")
  cat("This can take several minutes \n   \n \n")
  
  ctrol2 <-  trainControl(method = "LGOCV", p = 0.8, number = 5, savePredictions = T, verboseIter = TRUE )
  
  set.seed(825)
  cat("Fitting Random Forest model ...\n   \n \n")
  
  tunegrid <- expand.grid(mtry = 8:10)
  FDA <- train(clust ~ ., data = environR_clust[, 4:ncol(environR_clust)],  method = 'rf', ntree = 1000, tuneGrid = tunegrid, trControl = ctrol2)
  cat("finishing Random Forest ...\n   \n \n")
  
  cat("Classifying the rest of occurrences \n \n \n")
  to_assign <- environ_df[-muestra, ]
  rownames(to_assign) <- rownames(environ_df[-muestra,])
  
  to_assign$clust <-  predict(FDA, newdata = to_assign[,4:ncol(to_assign)]  )
  cat("converting predictions to a factor \n")
  to_assign$clust <- as.factor(to_assign$clust)
  
  environR_clust$clust <- as.factor(environR_clust$clust)
  
  cat("Binding both dataframes \n \n \n")
  environR_clust <- dplyr::bind_rows(environR_clust, to_assign)
  
}


cat( "Starting rasterization and saving process \n   \n  \n")

row.names(environR_clust) <- NULL
 SPF <- sp::SpatialPointsDataFrame( coords = environR_clust[,2:3], data = data.frame( cluster = environR_clust[,ncol(environR_clust)] ), proj4string = crs(SDM) )

r <- raster()
extent(r) <- extent(SDM)
crs(r) <- crs(SDM)
res(r) <- res(SDM)

to_rasterize <- raster::rasterize( x = SPF, y = r, field = as.numeric(SPF$cluster)  )
#raster mahalanobis
writeRaster(to_rasterize,filename= paste0( gap_outDir, "/ecogeo_hclust_mahalanobis.tif") , format="GTiff", overwrite = T )
 gc()
cat(paste("Process Done... Pls check the Path:", gap_outDir,"\n   \n \n"))

return(to_rasterize)


}# ENDCLUSTER FUNCTION

 



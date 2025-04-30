#####################################
# KERNEL DENSITY RASTER ESTIMATION  #
# Author: Chrystian C. Sosa         #
# Date: 2017/02/12                  #
#####################################
#https://cran.r-project.org/web/packages/adehabitatHR/adehabitatHR.pdf


###############

raster_kernel <- function(mask, occDir, out_path, kernel_method, scale){
  

if(!file.exists(out_path) ){
  
  cat("Reading mask and occurrences shapefile \n")
  msk <- raster(mask)
  
  ### Reading occurrences 
  
  occurrences <- shapefile(occDir) %>% 
    sp::SpatialPoints()
  crs(occurrences) <- crs(msk)
  
  #####
   if(kernel_method==2){
    cat("Using Adehabitat Kernel UD version","\n")
      
      df_mask <- raster::as.data.frame(msk, xy = TRUE )
      df_mask_sp <- sp::SpatialPoints(df_mask[, c(1,2)])
      crs(df_mask_sp) <- crs(msk)
      spixels <- sp::SpatialPixels(df_mask_sp)
      
      
      
    kernel <- adehabitatHR::kernelUD(xy = occurrences, h = "LSCV", grid= spixels, kern = "bivnorm")
    
    
    kernel<- raster::raster(kernel)
    kernel <- raster::mask(kernel, mask = msk)
    
    if(scale==T){
      kernel <- kernel/max(kernel[],na.rm=T)
      #kernel <- raster::scale(kernel,center=F,scale = T)
    } else {
      kernel <- kernel
      }
    
    } 
  ### Rasterizing density object
  
  crs(kernel) <- crs(msk)
 
  cat("Creating kernel classes raster.. \n")
  
  kernel[kernel[] == 0] <- NA
  kernel <- kernel * 10000
  qVal_1 <- raster::quantile(x = kernel[], probs = c(.9, 1), na.rm = T)
  knl_temp <- kernel
  knl_temp[which(knl_temp[] <= qVal_1[1])] <- NA
  qVal_2 <- raster::quantile(x = knl_temp, probs = c(.6, .95), na.rm = TRUE)
  kernel_class <- raster::reclassify(kernel, c(-Inf,qVal_1[1],1, qVal_1[1],qVal_2[2],2, qVal_2[2],Inf,3))
  
  ### Saving raster object
  cat("Saving raster objects","\n")
  raster::writeRaster(kernel_class, out_path , format = "GTiff")
  
  
}else{
  cat("kernel raster is already created... /n")
  kernel_class <- raster::raster(out_path)
 
}  
 
  
  cat("     ","\n")
  cat("DONE!","\n")
  cat("     ","\n")
  
  return(kernel_class)
  
}



#speciesList <- c("Mesoamerican","Andean")

### Calling paths to perform analysis
# root <- "U:"
# input_dir <- paste0(root,"/","Input_data")
# sdm_dir <- paste0(input_dir,"/","SDMs")
# out_dir <- input_dir
# ### Reading mask
# mask <- raster::raster(paste0(input_dir,"/","mask_wb_c_ant.tif"))
# 
# 
# ### Running for a Species list
# lapply(1:length(speciesList), function(i){
# 
# species=speciesList[[i]]
# 
# ### Reading occurrences
# occurrences <- read.csv(paste0(sdm_dir,"/","occurrences","/","occ_",species,".csv"),header=T)
# occurrences <- occurrences[,c("lon","lat")]
# 
# x<-RASTER_KERNEL(species,mask,occurrences,out_dir,spatstat=F,scale=T)
# })
# #x<-RASTER_KERNEL(species,mask,occurrences,out_dir,spatstat=F,standardize=T)

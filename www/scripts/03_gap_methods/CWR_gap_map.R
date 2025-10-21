# function to calculate gap map from SDM as for CWR methodology

CWR_gap_map <- function(sdm_path, occ_shp, out_file){
  #sdm_path <- "~/results_gap_analysis_2025/brassica_napusv2/pak_no_chn/results/brassica_napus/species_distribution/brassica_napus_sdm_median.tif"
  #occ_shp <- "~/results_gap_analysis_2025/brassica_napusv2/pak_no_chn/results/brassica_napus/species_distribution/occurrences.shp"

  sdm <- terra::rast(sdm_path)
  mtx <- matrix(c(-Inf, 0.7, 2,
                  0.7, 1, 3), ncol = 3, byrow = T)
  
  sdm <- terra::classify(sdm, mtx)
  occ <- terra::vect(occ_shp)
  
  terra::crs(occ) <- terra::crs(sdm)
  
  buff <- terra::buffer(occ, width = 50000 ) %>% 
    terra::rasterize(., sdm)
  
  buff[!is.na(buff)] <- 1
  
  gap_map <- terra::mask(sdm, buff, inverse = T) 
  gap_map <- sum(gap_map, buff, na.rm =T)
  
  gap_map[gap_map == 0 ] <- NA
  gap_map[gap_map >= 4] <- NA
  
  terra::writeRaster(gap_map, out_file , overwrite = T)
  
  return(gap_map)
    
}
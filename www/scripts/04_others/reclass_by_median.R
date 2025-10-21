
reclass_by_median <- function(rast, q = 0.5) {
  # Verifica que el objeto sea un SpatRaster
  if (!inherits(rast, "SpatRaster")) {
    stop("El argumento debe ser un objeto SpatRaster.")
  }
  
  # Calcula el cuantíl 50 (mediana)
  q50 <- terra::global(rast, fun = quantile, probs = q, na.rm = TRUE)[1, 1]
  
  # Reclasifica: 0 si <= mediana, 1 si > mediana
  reclass_matrix <- matrix(c(-Inf, q50, 0,
                             q50,  Inf, 1),
                           ncol = 3, byrow = TRUE)
  
  rast_reclass <- terra::classify(rast, rcl = reclass_matrix)
  
  return(rast_reclass)
}
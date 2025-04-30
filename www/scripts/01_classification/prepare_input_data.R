####################
# function to prepare the main input datafile and do predictions if it is needed
# authors: Andres camilo mendez, Harold achicanoy
###################


prepare_input_data <- function(data, 
    #data_path = choose.files( caption = "Select a valid .csv file"), 
                               col_number = NULL, 
                               mask = mask,
                               env_paths,
                               aux_paths){
  
  msk  <- terra::rast(mask)
  #data <- read.csv(data_path, header = T)
  
  
  #in case of col_number empty request a value from the user
  if(is.null(col_number)){
    warning("please enter the column number of response variable", immediate. = TRUE, noBreaks. = T)
    col_number <- readline(prompt="Enter column number of response variable: ") 
  }
  #check if the number enter by de user has the correct format(only a number)
  if(!grepl("^[0-9]+$", col_number)){
    warning("Only numbers are allowed", immediate. = TRUE, noBreaks. = T)
    col_number <- readline(prompt="Enter column number of response variable: ")
  }
  
  #select only the response variable and lat / long and create an ID column
  data <- data %>% 
    dplyr::select(., as.integer(col_number), 
                  matches("declat|^[L|l]atitude|lat", ignore.case = T), 
                  matches("declon|^[L|l]ongitude|lon|(L|l)ng", ignore.case = T),
                  everything(.)) 
  
  #change names
  names(data)[1:3] <- c("Y", "Latitude", "Longitude")
  cat("Cleaning missing lat/lon \n")
  data <- data %>%
    dplyr::filter(., !is.na(Latitude) | !is.na(Longitude)) 
  
  
  if("status" %in% names(data)){
    status <- data %>% dplyr::select(., matches("status", ignore.case = T)) %>% pull(1)
    data <- data %>% dplyr::select(-matches("status" , ignore.case = T))
  }else{status <- NULL}
  
  if("database_id" %in% names(data)){
    database_id <- data %>% dplyr::select(., matches("database_id", ignore.case = T))  %>% pull(1)
    data <- data %>% dplyr::select(-matches("database_id", ignore.case = T))
  }else{database_id <- NULL}
  
  if("source_db" %in% names(data)){
    source_db <- data %>% dplyr::select(., matches("source_db", ignore.case = T))  %>% pull(1)
    data <- data %>% dplyr::select(-matches("source_db", ignore.case = T))
  }else{source_db <- NULL}
   
  
  nms_to_remove <- data %>% 
    dplyr::select(-Y, -Longitude, -Latitude) %>% 
    names(.)
  
  data <- data%>% 
    dplyr::mutate(ID = 1:nrow(.))
  
  cat("Removing coordinates on Oceans/Seas \n")
  ocean_coords <- is.na(terra::extract(x = msk, y = data[,c("Longitude", "Latitude")])[,2] )
  cat("Removing duplicated coordinates \n")
  
  if(length(table(na.omit(data$Y))) == 1 & any(is.na(data$Y))){
    class_miss <- is.na(data$Y)
  }else{
    class_miss <- rep(FALSE, nrow(data))
  }
  
  #identificar duplicados por grupo y en los NA's
  if(length(table(na.omit(data$Y))) >= 2){
    
    rep <- lapply(unique(data$Y), function(gr){

      if(!is.na(gr)){
        tmp_df <- data %>% 
          dplyr::filter(Y == gr) %>% 
          dplyr::mutate(lon_lat_dups = duplicated(data.frame(Latitude, Longitude)),
                        msk_dups = duplicated(terra::cellFromXY(msk, .[, c("Longitude", "Latitude")])))  
          
        tmp_df$final_dups <- tmp_df %>% 
          dplyr::select(lon_lat_dups, msk_dups) %>% 
          apply(., MARGIN = 1, any)
        
        tmp_df <- tmp_df %>% 
          dplyr::select(-lon_lat_dups, -msk_dups)
        
      }else{
        tmp_df <- data %>% 
          dplyr::filter(is.na(Y)) %>% 
          dplyr::mutate(final_dups = duplicated(data.frame(Latitude, Longitude))) 
        
      }
      
      return(tmp_df)
      
    }) %>% 
      dplyr::bind_rows(.) %>% 
      dplyr::arrange(ID) %>% 
      dplyr::pull(final_dups)
      
  }else{
    rep <- duplicated(terra::cellFromXY(msk, data[, c("Longitude", "Latitude")]))
  }
  
  ##identificar duplicados ENTRE los grupos y los NA's
  
  if(length(table(na.omit(data$Y))) >= 2 & any(is.na(data$Y))){
    
    na_df <- data %>% 
      dplyr::filter(is.na(Y))
    
    unq_df <- data %>% 
      dplyr::filter(!is.na(Y)) %>% 
      dplyr::filter(!duplicated(data.frame(Latitude, Longitude))) %>% 
      dplyr::mutate(dm = duplicated(terra::cellFromXY(msk, .[, c("Longitude", "Latitude")]))) %>%
      dplyr::filter(!dm) %>%
      dplyr::select(-dm)
    
    tmp_df <- bind_rows(unq_df, na_df) %>% 
      dplyr::mutate(na_rep   = duplicated(data.frame(Latitude, Longitude)),
                    na_rep2  = duplicated(terra::cellFromXY(msk, .[, c("Longitude", "Latitude")])) )
    
    tmp_df$na_rep_f <- tmp_df %>% 
      dplyr::select(na_rep, na_rep2) %>% 
      apply(., MARGIN = 1, any)
    
    ids <- tmp_df %>% 
      dplyr::filter(na_rep_f) %>% 
      pull(ID)
    
    na_rep <- rep(FALSE, nrow(data))
    na_rep[ids] <- TRUE
    
    rm(tmp_df,na_df, unq_df, ids)
    
  }else{
    na_rep <- rep(FALSE, nrow(data))
  }
  
  
  cat("Extracting values from rasters \n")
  #climate extraction
  ax <- list.files(aux_paths, full.names = T)
  env <- list.files(env_paths, full.names = T)
  if(length(ax) != 0 ){
    fls <- c(env, ax)
  }else{
    fls <- env
  }
  
   current_clim_layer <- terra::rast(fls)
  
   
   
  clim_table <- terra::extract(current_clim_layer, data[, c("Longitude", "Latitude")])

  
  full_data <- data %>% 
        dplyr::bind_cols(., clim_table %>% dplyr::select(-matches("^ID", ignore.case = T))) 
  
  missings_vals <- complete.cases(clim_table %>% dplyr::select(-matches("^ID", ignore.case = T)))
  rm(clim_table)
  
  invalid <- ocean_coords | rep | !missings_vals | class_miss | na_rep

  full_data <- droplevels(full_data)

    #add status column if it exists
    
  if(!is.null(status)){
    full_data$status <- status
  }
  if(!is.null(database_id)){
    full_data$database_id <- database_id
  }
  if(!is.null(source_db)){
    full_data$source_db <- source_db
  }
 
    names(full_data)[1] <- "Y"
    full_data$ID <- NULL
    full_data$invalid <- invalid
    
    
    
   cat("cleaned done \n") 
    
    return(list(full_data = full_data, names_out  = nms_to_remove))
  
}#end fucntion

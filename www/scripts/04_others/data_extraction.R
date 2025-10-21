#conver raster to splitted tables
if(FALSE){
  msk <- msk <- terra::rast("D:/gap_analysis_shinyapps/www/masks/mask_world.tif")
#terra::rast("X:/runs/input_data/mask/mask_world.tif")
  stk <- terra::rast(purrr::map(list.files("C:/Users/acmendez/Downloads/new_rasters/", pattern = ".tif$", full.names = T), terra::rast))
  
  
  msk_df <- terra::as.data.frame(msk, xy = T, na.rm = F)
  msk_df <- msk_df %>% 
    dplyr::mutate( cellID = terra::cellFromXY(msk, as.matrix(msk_df[c("x", "y")]))) %>% 
    dplyr::select(cellID, everything(.))
  
  dfs <- terra::as.data.frame(stk, xy = T, na.rm = F)
  
  dfs_ed <-dfs %>% 
    dplyr::select(-x, -y) %>% 
    dplyr::bind_cols(msk_df) %>%
    dplyr::select(cellID, x, y, mask_world, everything(.)) %>% 
    tidyr::drop_na(mask_world) %>% 
    dplyr::arrange(cellID)
  
  nm <- round(nrow(dfs_ed)/4, 0)
  
  
  
  to_split <- c(rep("cat1", nm), rep("cat2", nm), rep("cat3", nm), rep("cat4", nm))
  dfs_ed$to_split <- to_split[1:nrow(dfs_ed)]
  
  splitted_df <- split(dfs_ed,dfs_ed$to_split)
  
  lapply(splitted_df, function(i){
    dt <- i %>% 
      data.table::as.data.table(.)
    
    dt <- dt[, -c("to_split")]
    
    dt <- dt[ , lapply(.SD, function(i){ifelse(is.nan(i), rep(NA, length(i)), i )}), .SDcols = names(dt)]
    
    
    
    saveRDS(dt, paste0("X:/others/generic_raster_tables/table_", min(dt[, cellID]),"-",max(dt[, cellID]),".rds" ))
    
    return(NULL)
    
  })
}


pacman::p_load(terra, tidyverse, data.table)

get_c_mask_df <- function(w_msk_path, shp_c, out_path){
  
  stopifnot("Load shapefile using terra package" = "SpatVector" %in% class(shp_c))
  
  w_mask <- terra::rast(w_msk_path) 
  
  crs_template <- terra::crs(w_mask)
  res_template <- terra::res(w_mask)
  
  c_mask <- w_mask %>% 
    terra::crop(., terra::ext(shp_c)) %>% 
    terra::mask(., shp_c)
  
  terra::writeRaster(c_mask, out_path, overwrite = T)
  
  c_mask_df <- terra::as.data.frame(c_mask, xy = T, na.rm = F) %>% 
    data.table::as.data.table(.)
  
  c_mask_df$cellID <- terra::cellFromXY(w_mask, as.matrix(c_mask_df[, c("x", "y")]))
  
  data.table::setorder(c_mask_df, cellID)
  
  rm(w_mask, c_mask)
  
  c_mask_df <- na.omit(c_mask_df, cols = "mask_world" )
  c_mask_df <- c_mask_df[, .(cellID, x, y, mask_world)]
  
  
  return(list(c_mask_df = c_mask_df, c_mask = c_mask))
  
}


extract_values <- function(c_mask_df, tables_path, c_mask, out_path){
  
  
  
  fls <- list.files(tables_path, full.names = T)
  
  mn <-  min(c_mask_df$cellID) #3431# # 
  mx <-  max(c_mask_df$cellID)#4912940
  
  ids_range_tbl <- tibble(paths = fls) %>% 
    dplyr::mutate(min_range = stringr::str_extract(fls, pattern = "_[0-9]+") %>% stringr::str_replace(., "_", "") %>% as.numeric(),
                  max_range = stringr::str_extract(fls, pattern = "-[0-9]+") %>% stringr::str_replace(., "-", "") %>% as.numeric()) %>% 
    dplyr::arrange(min_range) %>% 
    dplyr::mutate(min_in = ifelse( mn >= min_range &  mn <= max_range, TRUE, FALSE  ),
                  max_in = ifelse( (mx >= min_range & mx <= max_range)  , TRUE, FALSE) )
  
  pos <- which(ids_range_tbl$min_in):which(ids_range_tbl$max_in)
  
  ids_range_tbl <- ids_range_tbl[pos, ]
  
  
  if(nrow(ids_range_tbl)> 1){
    
    df_extracted <- lapply(1:nrow(ids_range_tbl), function(i){
      
      rw <- ids_range_tbl[i, ]
      
      to_extract <- readRDS(rw$paths)
      to_extract <- to_extract[, -c("x", "y", "mask_world")]
      
      to_ret <- to_extract[c_mask_df[ cellID >= rw$min_range & cellID <= rw$max_range,   ] , on = .(cellID)]
      rm(to_extract)
      
      
      
      data.table::setcolorder(to_ret, neworder = c(c("cellID", "x", "y", "mask_world" ), setdiff(names(to_ret), c("cellID", "x", "y", "mask_world" ))))
      
      return(to_ret)
      
    }) %>% 
      data.table::rbindlist(.)
    
    
  }else{
    
    to_extract <- readRDS(ids_range_tbl$paths)
    to_extract <- to_extract[, -c("x", "y", "mask_world")]
    
    df_extracted <- to_extract[c_mask_df[ cellID >= ids_range_tbl$min_range & cellID <= ids_range_tbl$max_range,   ] , on = .(cellID)]
    rm(to_extract)
    
    
    
    data.table::setcolorder(df_extracted, neworder = c(c("cellID", "x", "y", "mask_world" ), setdiff(names(df_extracted), c("cellID", "x", "y", "mask_world" ))))
    
    
  }
  
  ##create cropped rasters
  
  nms <- names(df_extracted)[-c(1,2,3,4)]
  for(i in nms){
    vrs <- c("x", "y", i)
    rf <- terra::rast(df_extracted[, ..vrs], crs = terra::crs(c_mask)) %>% 
      terra::resample(., c_mask)
    
    stopifnot("different raster resolution. " = all(terra::res(rf) == terra::res(c_mask)))
    
    terra::writeRaster(rf, paste0(out_path, "/", i, ".tif"), overwrite = T)
  }
  
  
  return(df_extracted)
  
}



c_mask_df <- get_c_mask_df(w_msk_path = "www/masks/mask_world.tif",
                           shp_c = terra::vect("www/world_shape_simplified/KEN.shp"),
                           out_path)



c_mask_extracted <- extract_values(c_mask_df = c_mask_df$c_mask_df, 
                                   tables_path = "www/generic_raster_tables/",
                                   c_mask = c_mask_df$c_mask,
                                   out_path)





pacman::p_load(RInno)

#RInno::install_inno(quick_start_pack = T, keep_install_file=TRUE)


 
 
 get_R2 <- function (app_dir = getwd(), R_version = paste0(">=", R.version$major, ".", R.version$minor)) {
   if (!dir.exists(app_dir)) 
     stop(glue::glue("{app_dir} does not exist."), call. = FALSE)
   R_version <- sanitize_R_version(R_version, clean = TRUE)
   latest_R_version <- readLines("https://cran.rstudio.com/bin/windows/base/", 
                                 warn = F) %>% stringr::str_extract("[1-4]\\.[0-9]+\\.[0-9]+") %>% 
     stats::na.omit() %>% unique()
   old_R_versions <- readLines("https://cran.rstudio.com/bin/windows/base/old/", 
                               warn = F) %>% stringr::str_extract("[1-4]\\.[0-9]+\\.[0-9]+") %>% 
     stats::na.omit()
   if (latest_R_version == R_version) {
     base_url <- glue::glue("https://cran.r-project.org/bin/windows/base/R-{R_version}-win.exe")
   }
   else {
     base_url <- glue::glue("https://cran.r-project.org/bin/windows/base/old/{R_version}/R-{R_version}-win.exe")
   }
   filename <- file.path(app_dir, glue::glue("R-{R_version}-win.exe"))
   if (file.exists(filename)) {
     cat("Using the copy of R already included:\n", 
         filename, "\n")
   }
   else {
     cat(glue::glue("Downloading R-{R_version} ...\n"))
     if (!R_version %in% c(latest_R_version, old_R_versions)) 
       stop(glue::glue("That version of R ({R_version}) is not listed on CRAN."), 
            call. = F)
     tryCatch(curl::curl_download(base_url, filename), error = function(e) {
       cat(glue::glue("\n          {base_url} is not a valid URL.\n\n          This is likely to have happened because there was a change in the URL.\n\n          This might have already been fixed in the latest version of RInno. Install it with remotes::install_github('ficonsulting/RInno').\n\n          If this doesn't help please submit an issue: {packageDescription('RInno', fields = 'BugReports')}\n\n          - Thanks!"))
     })
     if (!file.exists(filename)) 
       stop(glue::glue("{filename} failed to download."), 
            call. = FALSE)
   }
 }
 
 code_section2 <- function (iss, R_version = paste0(">=", R.version$major, 
                                   ".", R.version$minor)) 
 {
   if (length(R_version) == 0) {
     R_version = paste0(">=", R.version$major, ".", 
                        R.version$minor)
   }
   R_version <- sanitize_R_version(R_version)
   R_versions <- c(unique(stats::na.omit(stringr::str_extract(readLines("https://cran.rstudio.com/bin/windows/base/", 
                                                                        warn = F), "[1-4]\\.[0-9]+\\.[0-9]+"))), stats::na.omit(stringr::str_extract(readLines("https://cran.rstudio.com/bin/windows/base/old/", 
                                                                                                                                                               warn = F), "[1-4]\\.[0-9]+\\.[0-9]+")))
   inequality <- substr(R_version, 1, attr(regexpr("[<>=]+", 
                                                   R_version), "match.length"))
   R_version <- gsub("[<>=]", "", R_version)
   version_specs <- paste0("numeric_version('", R_versions, 
                           "')", inequality, "numeric_version('", R_version, 
                           "')")
   if (!R_version %in% R_versions && interactive()) 
     stop(glue::glue("R version - {R_version} - was not found on CRAN. Please use `R_version` to specify one that is or let us know if you think you received this message in error: \n\nhttps://github.com/ficonsulting/RInno/issues"), 
          call. = FALSE)
   results <- unlist(lapply(version_specs, function(x) eval(parse(text = x))))
   acceptable_R_versions <- paste0(glue::glue("RVersions.Add('{R_versions[results]}');"), 
                                   collapse = "\n  ")
   code_file <- paste0(readLines(system.file("installation/code.iss", 
                                             package = "RInno")), collapse = "\n")
   glue::glue("{iss}\n  {code_file}\n  // Initialize the values of supported versions\n  RVersions := TStringList.Create; // Make a new TStringList object reference\n  // Add strings to the StringList object\n  {acceptable_R_versions}\n\nend;\n\n// Procedure called by InnoSetup when it is closing\nprocedure DeinitializeSetup();\nbegin\n  RVersions.Free;\nend;\n  ")
 }
 
 create_app2 <- function (app_name = "myapp", app_dir = getwd(), dir_out = "RInno_installer", 
                          pkgs = c("jsonlite", "shiny", "magrittr"), 
                          pkgs_path = "bin", repo = "https://cran.rstudio.com", 
                          remotes = "none", locals = NULL, app_repo_url = "none", 
                          auth_user = "none", auth_pw = "none", auth_token = NULL,#github_pat(), 
                          user_browser = "electron", include_R = FALSE, include_Pandoc = FALSE, 
                          include_Chrome = FALSE, include_Rtools = FALSE, R_version = paste0(">=", 
                                                                                             R.version$major, ".", R.version$minor), Pandoc_version = rmarkdown::pandoc_version(), 
                          Rtools_version = "3.5", overwrite = TRUE, force_nativefier = TRUE, 
                          nativefier_opts = c(), ...) 
 {
   dots <- list(...)
   if (!is.null(locals)) {
     warning("locals is deprecated. Please use pkgs instead.", 
             call. = FALSE)
     pkgs <- pkgs %>% standardize_pkgs %>% add_pkgs(locals)
   }
   if (user_browser != "electron") {
     warning(glue::glue("user_browser = {glue::double_quote(user_browser)} will be deprecated in the next release. Please use user_browser = \"electron\" in the future."), 
             call. = FALSE)
   }
   if (include_Chrome) {
     warning("include_Chrome will be deprecated in the next release. Please use user_browser = \"electron\"", 
             call. = FALSE)
   }
   if (!is.null(dots$ping_site)) {
     warning("ping_site is deprecated in favor of self-contained dependency management in the .exe.", 
             call. = FALSE)
   }
   if (class(app_name) != "character") 
     stop("app_name must be a character.", call. = FALSE)
   if (class(dir_out) != "character") 
     stop("dir_out must be a character.", call. = FALSE)
   logicals <- c(include_Chrome = class(include_Chrome), include_Pandoc = class(include_Pandoc), 
                 include_R = class(include_R), include_Rtools = class(include_Rtools), 
                 overwrite = class(overwrite))
   failed_logicals <- !logicals %in% "logical"
   if (any(failed_logicals)) {
     stop(glue::glue("{names(logicals[which(failed_logicals)])} must be TRUE/FALSE."), 
          call. = F)
   }
   if (!dir.exists(app_dir)) 
     dir.create(app_dir)
   R_version <- sanitize_R_version(R_version)
   copy_installation(app_dir, overwrite)
   if (include_R) 
     get_R2(app_dir, R_version)
   if (include_Pandoc) 
     get_Pandoc(app_dir, Pandoc_version)
   if (include_Chrome) 
     get_Chrome(app_dir)
   if (include_Rtools) 
     get_Rtools(app_dir, Rtools_version, R_version)
   if (user_browser == "electron" && interactive()) {
     if (force_nativefier) {
       #dots$app_icon <- "default.ico"
       nativefy_app(app_name, app_dir, nativefier_opts, 
                    app_icon = dots$app_icon)
     }
     else {
      
       if (!dir.exists(file.path(app_dir, "nativefier-app"))) 
         nativefy_app(app_name, app_dir, nativefier_opts, 
                      app_icon = dots$app_icon)
       cat("\nUsing previously built electron app...\n")
     }
   }
   create_bat(app_name, app_dir)
   create_config(app_name, app_dir, pkgs = pkgs, pkgs_path = pkgs_path, 
                 remotes = remotes, repo = repo, error_log = dots$error_log, 
                 app_repo_url = app_repo_url, auth_user = auth_user, auth_pw = auth_pw, 
                 auth_token = auth_token, user_browser = user_browser)
   start_iss(app_name) %>% directives_section(include_R, R_version, 
                                              include_Pandoc, Pandoc_version, include_Chrome, include_Rtools, 
                                              Rtools_version, app_version = dots$app_version, publisher = dots$publisher, 
                                              main_url = dots$main_url) %>% setup_section(app_dir, 
                                                                                          dir_out, app_version = dots$app_version, default_dir = dots$default_dir, 
                                                                                          privilege = dots$privilege, info_before = dots$info_before, 
                                                                                          info_after = dots$info_after, setup_icon = dots$setup_icon, 
                                                                                          inst_pw = dots$inst_pw, license_file = dots$license_file, 
                                                                                          pub_url = dots$pub_url, sup_url = dots$sup_url, upd_url = dots$upd_url) %>% 
     languages_section %>% tasks_section(desktop_icon = dots$desktop_icon) %>% 
     icons_section(app_dir, app_desc = dots$app_desc, app_icon = dots$app_icon, 
                   prog_menu_icon = dots$prog_menu_icon, desktop_icon = dots$desktop_icon) %>% 
     files_section(app_name, app_dir, user_browser, file_list = dots$file_list) %>% 
     run_section(dots$R_flags) %>% code_section2(R_version) %>% 
     writeLines(file.path(app_dir, paste0(app_name, ".iss")))
   check_app(app_dir, pkgs_path)
 }
 
 c("shiny", "shinydashboard", "leaflet", "raster", "rgdal", "rgeos", "sp", "rsconnect", "ggplot2", 
   "shinyFiles", "shinyBS", "shinyjs", "yaml", "shinyWidgets", "rmarkdown", "bsplus", "tidyverse", "shinydashboardPlus", "DT",
   "tcltk", "adehabitatHR",   "raster", "rgdal", "sdm", "dismo",  "rgeos", "distances",   "sp", "shiny", 
   "tidyverse", "rlang", "sf", "gdistance", "caret", "earth", "fastcluster",  "FactoMineR", "deldir",
   "bindrcpp",  "pROC", "maxnet", "usdm", "mltools", "ISLR", "nnet","ranger",  "googleVis", "terra")
 
 create_app2(app_name = "myapp3", 
             app_dir = "D:/OneDrive - CGIAR/Documents/R/R-4.0.2/library/RInno/app", 
             pkgs = c("shiny", "datasets"),
             include_R = F,
             include_Pandoc = T)
 2
 compile_iss()
 
 
 RInno::compile_iss()
 
 
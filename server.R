#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

suppressMessages(if(!require(pacman)){install.packages("pacman");library(pacman)}else{library(pacman)})
pacman::p_load(tcltk, adehabitatHR,   raster, sdm, dismo, distances,   sp, shiny,
               tidyverse, rlang, sf, gdistance, earth, fastcluster, deldir,
                bindrcpp,  pROC, maxnet, mltools, ISLR, nnet, HDclassif, rJava,
               ranger, plotly, terra, ModelMetrics, e1071, caret, 
               writexl, rmarkdown, knitr, kableExtra, googleVis, tidysdm, spatialsample,
               rsample, workflowsets, tune, parsnip, biomod2)
#usdm
#rminer
#gbm
#caretEnsemble
#caret
#FactoMineR

#check if tinytex is installed
# if(!tinytex::is_tinytex()){
#   tinytex::install_tinytex()
# }
# setting global variables
g <- gc(reset = T); rm(list = ls()); options(warn = -1); options(scipen = 999)

shp <- sf::st_read("www/world_shape_simplified/all_countries_simplified_new.shp") #%>% 
#shp_new <- sf::st_read("www/world_shape_simplified/all_countries_simplified_new.shp")
  #dplyr::filter(ISO3 != "ATA")

scrDir <- "www/scripts"
coor_sys <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
options(shiny.maxRequestSize=100*1024^2)
options(scipen = 999)

source("www/scripts/04_others/helpers.R", local = TRUE)
source("www/scripts/04_others/03_generate_report.R", local = TRUE)


av_fls <- list.files("www/scripts", recursive = TRUE, pattern = ".R$", full.names = T)
av_fls <- av_fls[!grepl("04_others", av_fls)]
av_fls <- av_fls[!grepl("validation_workflow.R", av_fls)]
lapply(av_fls, source)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  
#**********************************
##### Reactive values lists #####
#*********************************
  
  resources <- reactiveValues()
  #shp_custm <- reactiveValues()
  clicklist <- reactiveValues(ids = vector(), names = vector())
  paths     <- reactiveValues()

#**********************************
### monitoring folder system ####
#**********************************

  ##monitorear los folders de la mascara si existen o no
shiny::observe({
  
  req(paths$cropDir)
  
  drs <- list.dirs(paths$cropDir, recursive = F, full.names = F)
  
  if(length(drs) != 0){
    lg <- lapply(1:length(drs), function(i){
     ret <- file.exists(paste0(paths$cropDir, "/" ,drs[i], "/input_data/mask_", drs[i], ".tif"))
     return(ret)
    }) %>% unlist
    
    resources$mask_avaliables <- drs[lg]
   
  output$avaliable_masks <- renderUI({
    shinyWidgets::pickerInput(
      inputId = "mask_avaliable",
      label   = "Select an existing mask: ", 
      choices = resources$mask_avaliables,
      options = list(
        style = "btn-primary")
    )
  })  
   
    
  }else{
    resources$mask_avaliables <- character(0)
  }
  
})
 
  #monitorear todo el sistema de archivos 
  shiny::observe({

  req(paths$results_dir,  paths$input_data_dir)
    #results per race dir
       
    if(!file.exists(paths$results_dir)){dir.create(paths$results_dir, recursive = TRUE)}
            #child folders
    if(!is.null(paths$occName)){
      paths$sp_results <- paste(paths$results_dir, paths$occName, sep = .Platform$file.sep);if(!file.exists(paths$sp_results)){dir.create(paths$sp_results, recursive = TRUE)}
      paths$gap_outDir <- paste0(paths$sp_results, "/gap_scores");if(!file.exists(paths$gap_outDir)){dir.create(paths$gap_outDir, recursive = TRUE)}
      paths$delaDir    <- paste0(paths$gap_outDir, "/delaunay");if(!file.exists(paths$delaDir)){dir.create(paths$delaDir, recursive = TRUE)}
      paths$sdm_outDir <- paste0(paths$sp_results, "/species_distribution");if(!file.exists(paths$sdm_outDir)){dir.create(paths$sdm_outDir, recursive = TRUE)}
      paths$cvfoldsDir <- paste0(paths$sdm_outDir, "/sdm_cv_folds");if(!file.exists(paths$cvfoldsDir)){dir.create(paths$cvfoldsDir, recursive = TRUE)}
      paths$pseudoDir  <- paste0(paths$sdm_outDir, "/pseudo_absences");if(!file.exists(paths$pseudoDir)){dir.create(paths$pseudoDir, recursive = TRUE)}
      paths$gap_valDir <- paste0(paths$sp_results, "/gaps_validation");if(!file.exists(paths$gap_valDir)){dir.create(paths$gap_valDir, recursive = TRUE)}
      paths$occ_group  <- paste0(paths$sp_results, "/", paths$occName,"_occurrences.csv")
      paths$occ_shp    <- paste0(paths$sdm_outDir, "/occurrences.shp")
      paths$occ_csv    <- paste0(paths$sdm_outDir, "/occurrences.csv")
      paths$pseudo_file<- paste0(paths$pseudoDir,  "/pseudo_absences_", paths$occName, ".csv")
      paths$smd_var_selected <- paste0(paths$sdm_outDir, "/sdm_variables_selected.csv")
      paths$sdm_calibration  <- paste0(paths$sdm_outDir, "/tunned_sdm_parms.csv")
      paths$sdm_occ_path     <- paste0(paths$sdm_outDir, "/", paths$occName, "_sdm_median.tif")
      paths$knl_out_path     <- paste0(paths$gap_outDir, "/kernel_density_classes.tif")
      paths$cost_out_path    <- paste0(paths$gap_outDir, "/cost_dist_score.tif")
      paths$dela_out_path    <- paste0(paths$gap_outDir, "/network_score.tif")
      paths$envi_out_path    <- paste0(paths$gap_outDir, "/environmental_score.tif")
      
      paths$cost_class_out_path    <- paste0(paths$gap_outDir, "/cost_dist_gap_map.tif")
      paths$dela_class_out_path    <- paste0(paths$gap_outDir, "/network_gap_map.tif")
      paths$envi_class_out_path    <- paste0(paths$gap_outDir, "/environmental_gap_map.tif")
      paths$final_gap_map_path     <- paste0(paths$sp_results, "/", paths$occName,"_final_gap_map.tif")
      paths$CWR_gap_map_path       <- paste0(paths$sp_results, "/", paths$occName,"_CWR_gap_map.tif")
      paths$out_report_file        <- paste0(paths$sp_results, "/results_report.Rmd")
      
      if(file.exists(paths$occ_group)){
        resources$spData <-  read.csv(paths$occ_group, header = T, stringsAsFactors = F)
      }
      if(file.exists(paths$pseudo_file)){
        resources$pseudo_abs <- read.csv(paths$pseudo_file, header = T, stringsAsFactors = F)
      }
      if(file.exists(paths$smd_var_selected)){
        resources$var_names <- read.csv(paths$smd_var_selected) %>% 
          dplyr::pull(x)
      }
      if(file.exists(paths$CWR_gap_map_path)){
        resources$check_gap_scores <- "error"
      }
      
      #save automatically inputs and results to be restored 
      saveRDS(paths, paste0(paths$results_dir, "/", paths$occName ,"_Rsession.rds") )
      print("actualice el archivo rds")
    }
             #input data dir folders creation
    
    if(!file.exists(paths$input_data_dir)){dir.create(paths$input_data_dir, recursive = TRUE)}          
            #child Folders
            paths$occDir       <- paste0(paths$input_data_dir, "/occurrences_files");if(!file.exists(paths$occDir)){dir.create(paths$occDir, recursive = TRUE)}
            paths$generic_dir  <- paste0(paths$input_data_dir, "/generic_rasters");if(!file.exists(paths$generic_dir)){dir.create(paths$generic_dir, recursive = TRUE)} 
            paths$aux_dir      <- paste0(paths$input_data_dir, "/auxiliary_rasters");if(!file.exists(paths$aux_dir)){dir.create(paths$aux_dir, recursive = TRUE)}
            paths$mask_path    <- paste(paths$input_data_dir, paste0("mask_" , paths$region_name, ".tif") , sep = .Platform$file.sep)
            paths$db_file      <- paste0(paths$occDir, "/", paths$crop_name, "_occurrence_DB.csv")
            paths$original_path<- paste0(paths$occDir, "/", paths$crop_name, "_original_data_in.csv")
            paths$db_summ      <- paste0(paths$occDir, "/", paths$crop_name, "_summary.rds")
            
            if(length(list.files(.libPaths(), pattern = "dismo", full.names = T)) == 1){
              paths$maxent_jar   <-paste0(list.files(.libPaths(), pattern = "dismo", full.names = T), "/java/maxent.jar")
            }
             
            
            if(file.exists(paths$db_file)){
              resources$cleaned_data <- read.csv(paths$db_file, header = T, stringsAsFactors = F)
              resources$group_names <- unique(na.omit(resources$cleaned_data$Y))
            
            }
            
           
 
}, priority = 0)
  


 #app info seccion 
  output$messageMenu <- renderMenu({
    
    crop    <- ifelse(is.null(paths$crop_name), "Not specified yet", paths$crop_name)
    occName <- ifelse(is.null(paths$occName), "Not specified yet", paths$occName)
    region  <- ifelse(is.null(paths$region_name), "Not specified yet", paths$region_name)
    
    msgs <- list(messageItem(from = "Crop Name" , message = crop,    icon = icon("seedling","fas fa" , verify_fa = FALSE)),
                 messageItem(from = "Group Name", message = occName, icon = icon("sitemap","fas fa" , verify_fa = FALSE)),
                 messageItem(from = "Region"    , message = region,  icon = icon("map-pin","fas fa" , verify_fa = FALSE)) )
   
   
     
     icono <- icon("circle-info", verify_fa = FALSE)
   
    # This is equivalent to calling:
    #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = "messages", .list = msgs, icon = icono)
    
  })
  
#************************* 
#### Restore session #####
#*************************
  output$restoreSession <- renderMenu({
    msgs2 <- list(notificationItem(
      text = "Restore session",
      icon = icon("rotate-left","fas fa") 
    ) %>%
      tagAppendAttributes(., id = "restore"))
    
    dropdownMenu(type = "notifications", .list = msgs2, icon = icon("user-gear","fas fa", verify_fa = FALSE) )
  })
  
  #when user click on restore then open a modal dialog to import rsession file
    shinyjs::onclick("restore", expr = function(){
      output$modal1 <- renderUI({
        showModal(modalDialog(
          title = tags$strong("Restore Session"),
          tags$div(
            tags$div(style = "float:left;", tags$img(src = 'restore-icon.jpg', eigth = "45px", width = "45px")),
            tags$div(style = "margin-left: 50px;",tags$h5("
            Did your last Session closed unexpectedly or went something wrong? You can restore the variables and results from a previous session, or start a
            new session if you need to process other groups, class or genetic group."))
          ),
          tags$hr(),
          fileInput("restorePath", "Select  Rsession.rds file:"),
          bsButton("accept_restore", size="default",label = "Restore", block = F, style="primary"),
          bsTooltip(id = "Restore", title = "Restore Previous Session", placement = "right", trigger = "hover"),
          easyClose = FALSE,
          footer = modalButton("Cancel")
        ))
    })
  })

    #save recursively response variables
    observeEvent(input$accept_restore, {

      tryCatch({
        
        rsession_path <- input$restorePath$datapath
        if(nchar(rsession_path != 0)){
          
          rsd <- readRDS(rsession_path)
          #selec main variables to restore
       
          paths$results_dir    <- rsd$results_dir
          paths$input_data_dir <- rsd$input_data_dir
          paths$occName        <- rsd$occName
          paths$crop_name      <- rsd$crop_name
          paths$region_name    <- rsd$region_name
         
          
          sendSweetAlert(
            session = session,
            title = "Session restored !!",
            text = paste("Your Session for", paths$occName," was successfully restored."),
            type = "success"
          )
          print("Reactive values were uploaded.")
          
        }else{
          safeError("Not a valid file Path detected.")
        }
        
      }, error = function(e){
        sendSweetAlert(
          session = session,
          title = "Error !!",
          text = paste("Your Session could not be restored.  \n", e),
          type = "error"
        )
      })
  
    })
#*************************
##### Root folders ###### 
#***********************
    
  #create button to select container folder to create dirs and download data
  shinyDirChoose(input , id =  "select_path_btn", updateFreq = 0, session = session,
                 defaultPath = "", roots =  c('Documents' = Sys.getenv("HOME"), 'Local Disk' = Sys.getenv("HOMEDRIVE") ))
  
  
  #boton para seleccionar el directorios
  observeEvent(input$select_path_btn, {
    text_path <- parseDirPath(roots =  c('Documents' = Sys.getenv("HOME"), 'Local Disk' = Sys.getenv("HOMEDRIVE") ), input$select_path_btn)
    
    paths$root_path <- as.character(text_path)
    
    resources$n_crops <- list.dirs(paths$root_path, full.names = F, recursive = F)
    
    output$crop_def <- renderUI({
     
        if(length(resources$n_crops) == 0){
          textInput(inputId =  "set_crop_name", label = "Write Crop name") %>%
            short_info(input = ., place = "top",title = "Name of the major crop")
        }else{
          
          shinyWidgets::pickerInput(
            inputId = "set_crop_name",
            label = "Crops found:", 
            choices = c(as.character(resources$n_crops), "add new crop"),
            options = list(
              style = "btn-primary")
          )
          
         
          
        }
        
      
      
    })
    
    
    #.GlobalEnv$baseDir <- as.character(text_path)
    updateTextInput(session, "selected.root.folder", 
                    label = "Folder selected",
                    value = as.character(text_path)
    )
    
  })
 

  #boton para crear los directorios y descargar los scripts desde github
  observeEvent(input$create_dirs,{
  
    if(!is.null(paths$root_path)){
     
      paths$crop_name <- as.character(input$set_crop_name)
      
      if( nchar(input$set_crop_name2) != 0){
        paths$crop_name <- input$set_crop_name2
      }
      #paths$occName <- input$set.level.name
     
      #print(paths$crop_name)
      if(nchar(paths$crop_name) > 0  ){
        
        ### create crop dir
        
        paths$cropDir <- paste0(paths$root_path,"/", paths$crop_name)
        if(!file.exists(paths$cropDir)){dir.create(paths$cropDir, recursive = TRUE)}
        
        sendSweetAlert(
          session = session,
          title = "Success !!",
          text = "Main directories were successfully created.",
          type = "success"
        )
        
        
        updateButton(session, "create_dirs",label = "Selected", style = "success")
      }else{   
        sendSweetAlert(
          session = session,
          title = "Error !!",
          text = "Please write a valid Crop or Race name.",
          type = "error"
        )
      
      }
        

    }else{
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, select a valid Folder path from Choose Dir buttom.",
        type = "error"
      )
    }
 
    
  })#end boton para crear directorios y descargar scripts
 
  
  output$map_selector <- renderLeaflet({
    leaflet("map_selector") %>% 
      setView(lat= 0, lng = 0, zoom = 1) %>% 
      addTiles(options = providerTileOptions(noWrap = TRUE) )%>% 
      addPolygons(data = shp,  stroke = T,color = "#4682B4", weight = 1, smoothFactor = 0.2,  #adicionar el archivo .shp al mapa y hacer que brillen cuadno son seleccionados
                  label = ~as.character(shp$ISO3),
                  layerId = as.character(shp$ISO3),
                  opacity = 0.5, 
                  fillOpacity = 0.5,
                  fillColor = "#D6DEE6",
                  highlightOptions = highlightOptions(color = "#696262", weight = 2,
                                                      bringToFront = F) )
      
  })
  
##*************************************************
####  features shape selector or MASK creator ####
#************************************************** 
observeEvent(c(input$map_selector_shape_click, input$area_selector),{
    proxy <- leafletProxy("map_selector")
    print(input$area_selector)
    cont <- as.numeric(input$area_selector)
    
    if(cont == 0 | cont == 8){
      resources$shp_custm <- shp %>% 
        dplyr::filter(REGION != cont)
      proxy %>% 
        clearShapes() %>% 
        addPolygons(data= resources$shp_custm, stroke = T,color = "#2690EF", weight = 1, smoothFactor = 0.2,  #adicionar el archivo .shp al mapa y hacer que brillen cuadno son seleccionados
                    opacity = 0.5, fillOpacity = 0.5,
                    layerId = resources$shp_custm$ISO3,
                    fillColor = "#D6DEE6",
                    highlightOptions = highlightOptions(color = "#696262", weight = 2,
                                                        bringToFront = TRUE) )
    }else{
      resources$shp_custm <- shp %>% 
        dplyr::filter(REGION == cont )
      proxy %>% 
        clearShapes() %>% 
        addPolygons(data= resources$shp_custm, stroke = T,color = "#2690EF", weight = 1, smoothFactor = 0.2,  #adicionar el archivo .shp al mapa y hacer que brillen cuadno son seleccionados
                    opacity = 0.5, fillOpacity = 0.5,
                    layerId = resources$shp_custm$ISO3,
                    fillColor = "#D6DEE6",
                    highlightOptions = highlightOptions(color = "#696262", weight = 2,
                                                        bringToFront = TRUE) )
    }
    
    
    
    if(as.numeric(input$area_selector) == 8){
      
      click <- input$map_selector_shape_click
      # Store click IDS in reactiveValue
      clicklist$ids <- c(clicklist$ids, click$id)
  
      if(any(duplicated(clicklist$ids))){
        
        val <- clicklist$ids[duplicated(clicklist$ids)]
        clicklist$ids <- clicklist$ids[!clicklist$ids %in% val]
      }
      
      resources$shp_custm <- shp %>% 
        dplyr::filter(ISO3 %in% clicklist$ids)
      
      clicklist$names <- resources$shp_custm$NAME
     
      proxy %>% 
        addPolygons( data= resources$shp_custm,
                     layerId = resources$shp_custm$ISO3,
                     color = "#444444", weight = 0.5, smoothFactor = 0.5,
                     opacity = 1.0, fillOpacity = 0.5,
                     fillColor ="yellow" )
    }
    
  }) 
 
  #actualizar los paises seleccionados
shiny::observe({ 
    if(as.numeric(input$area_selector) == 8){
      updateAwesomeCheckboxGroup(session, inputId = "chk_bx_gr", label = "Countries selected:", choices = clicklist$names, selected = clicklist$names)
   
    }
  })

 observeEvent(input$create_mask,{
   
   
   if(is.null(paths$crop_name)){
     sendSweetAlert(
       session = session,
       title = "Error !!",
       text = "Please, write a Crop name and select a root folder from your computer.",
       type = "error"
     )
     system.time(0.5)
     
     updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
     
   }else if(grepl("\\W|\\$|\\?|^[0-9]", input$mask_name)){
     sendSweetAlert(
       session = session,
       title = "Error !!",
       text = "Especial characters found in region name.",
       type = "error"
     )
     
   }else{
   
  
      if(nchar(as.character(input$mask_name)) != 0){
      
          
     paths$region_name <- as.character(input$mask_name)
     paths$baseDir     <- paste0(paths$cropDir, "/", paths$region_name)
     if(!file.exists(paths$baseDir)){dir.create(paths$baseDir, recursive = TRUE)}
    
     paths$input_data_dir     <- paste0(paths$baseDir,"/input_data")
     if(!file.exists(paths$input_data_dir)){dir.create(paths$input_data_dir, recursive = TRUE)}
     
     paths$generic_dir  <- paste0(paths$input_data_dir, "/generic_rasters")
     if(!file.exists(paths$generic_dir)){dir.create(paths$generic_dir, recursive = TRUE)}
    
     
     paths$results_dir    <- paste0(paths$baseDir, "/results/")
     if(!file.exists(paths$results_dir)){dir.create(paths$results_dir, recursive = TRUE)}
    
     
     paths$mask_path <- paste(paths$input_data_dir, paste0("mask_" , paths$region_name, ".tif") , sep = .Platform$file.sep)
    
     
     ## poner bussy indicator
     withBusyIndicatorServer("create_mask", {
       
      #print(terra::ext(resources$shp_custm))
       c_mask_df <- get_c_mask_df(w_msk_path = "www/masks/mask_world.tif",
                                  shp_c = terra::vect(resources$shp_custm),
                                  out_path = paths$mask_path)
       
       resources$c_mask_df <- c_mask_df$c_mask_df
       
       c_mask_extracted <- extract_values(c_mask_df = resources$c_mask_df, 
                                          tables_path = "www/generic_raster_tables/",
                                          c_mask = terra::rast(paths$mask_path),
                                          out_path = paths$generic_dir)
       
     # terra::rast("www/masks/mask_world.tif") %>% 
     #   terra::crop(., y = terra::ext(shp_custm)) %>%
     #   terra::ext(., shp_custm) %>%
     #   terra::writeRaster(., paths$mask_path, overwrite = T)
     
       }) 
     updateButton(session, "create_mask",label = "Mask created",style = "success")
    
     
   }else{
     sendSweetAlert(
       session = session,
       title = "Warning !!",
       text = "Please,write a valid mask name.",
       type = "warning"
     )
     
   }
 }#end else check global vars e.j paths$baseDir 
   
 })
 
 
 #####import mask in case that it are already created###
 
 observeEvent(input$import_mask,{
  
   req(input$mask_import)
   
   tryCatch(
     {
       r <- terra::rast(input$mask_import$datapath)
       rast_name <- input$mask_import$name
       
       
       if(is.null(paths$crop_name) ){
         sendSweetAlert(
           session = session,
           title = "Error !!",
           text = "Please, write a Crop name and select a root folder from your computer.",
           type = "error"
         )
         system.time(0.5)
         
         updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
         
       }else{
         
       if(all(grepl(".tif$", rast_name), grepl("^mask_", rast_name))){ 
         #define region rv when mask is uploaded
         withBusyIndicatorServer("import_mask",{  
           
         if(is.null(paths$region_name)){ 
           
           reg <- gsub("^mask_|.tif$", "", rast_name)
         print(reg)
         paths$region_name <- reg
         
         paths$baseDir     <- paste0(paths$cropDir, "/", paths$region_name)
         if(!file.exists(paths$baseDir)){dir.create(paths$baseDir, recursive = TRUE)}
         
         paths$input_data_dir     <- paste0(paths$baseDir,"/input_data")
         if(!file.exists(paths$input_data_dir)){dir.create(paths$input_data_dir, recursive = TRUE)}
         
         paths$generic_dir  <- paste0(paths$input_data_dir, "/generic_rasters")
         if(!file.exists(paths$generic_dir)){dir.create(paths$generic_dir, recursive = TRUE)}
         
        
         paths$results_dir    <- paste0(paths$baseDir, "/results/")
         if(!file.exists(paths$results_dir)){dir.create(paths$results_dir, recursive = TRUE)}
        
         paths$mask_path <- paste(paths$input_data_dir, paste0("mask_" , paths$region_name, ".tif") , sep = .Platform$file.sep)
         
         source(paste("www/scripts/01_classification/crop_raster_new.R", sep = ""), local = F)
         c_mask_df <- get_c_mask_df(w_msk_path = "www/masks/mask_world.tif",
                                    shp_c = r,
                                    out_path = paths$mask_path)
         
         resources$c_mask_df <- c_mask_df$c_mask_df
         
         
         
         c_mask_extracted <- extract_values(c_mask_df = resources$c_mask_df, 
                                            tables_path = "www/generic_raster_tables/",
                                            c_mask = terra::rast(paths$mask_path),
                                            out_path = paths$generic_dir)
         
         
         }
       })
         updateButton(session, "import_mask",label = "Imported",style = "success")
         
       }else{
         stop(safeError("Invalid file extension (only .tif files are allowed)"))
       }
         
         
       } 
     },
     error = function(e) {
       # return a safeError if a parsing error occurs
       #stop(safeError(e))
       sendSweetAlert(
         session = session,
         title = "Error !!",
         text = "Please, Select a valid raster file.",
         type = "error"
       )
     }
   )
   
 })
 
 ##### seleccionar una mascara si ya esta creada
 observeEvent(input$select_mask, {
   
   
   if(is.null(paths$crop_name)){
     sendSweetAlert(
       session = session,
       title = "Error !!",
       text = "Please, write a Crop name and select a root folder from your computer.",
       type = "error"
     )
     system.time(0.5)
     
     updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
     
   }else{
   
     withBusyIndicatorServer("select_mask",{  
       
     paths$region_name <- as.character(input$mask_avaliable)
     paths$baseDir     <- paste0(paths$cropDir, "/", paths$region_name)
     if(!file.exists(paths$baseDir)){dir.create(paths$baseDir, recursive = TRUE)}
     
     paths$input_data_dir     <- paste0(paths$baseDir,"/input_data")
     if(!file.exists(paths$input_data_dir)){dir.create(paths$input_data_dir, recursive = TRUE)}
     
     paths$generic_dir  <- paste0(paths$input_data_dir, "/generic_rasters")
     if(!file.exists(paths$generic_dir)){dir.create(paths$generic_dir, recursive = TRUE)}
     
     
     paths$results_dir    <- paste0(paths$baseDir, "/results/")
     if(!file.exists(paths$results_dir)){dir.create(paths$results_dir, recursive = TRUE)}
     
     
     paths$mask_path <- paste(paths$input_data_dir, paste0("mask_" , paths$region_name, ".tif") , sep = .Platform$file.sep)
     
     if(length(list.files(paths$generic_dir)) == 0){
       source(paste("www/scripts/01_classification/crop_raster_new.R", sep = ""), local = F)
       c_mask_df <- get_c_mask_df(w_msk_path = "www/masks/mask_world.tif",
                                  shp_c = terra::vect(resources$shp_custm),
                                  out_path = paths$mask_path)
       
       resources$c_mask_df <- c_mask_df$c_mask_df
       
       
       
       c_mask_extracted <- extract_values(c_mask_df = resources$c_mask_df, 
                                          tables_path = "www/generic_raster_tables/",
                                          c_mask = terra::rast(paths$mask_path),
                                          out_path = paths$generic_dir)
     }
     
     updateButton(session, "select_mask",label = "File selected",style = "success")
     
     })#end with bussy indicator
     
     }
   
 })
 
#*******************************************
##### CARGAR ARCHIVOS AUXILIARES      #####
#*******************************************
 
 observeEvent(input$load_rast,{
   
   req(input$aux_files)
   tryCatch({
     if(is.null(paths$crop_name)|is.null(paths$mask_path) | !file.exists(paths$mask_path) | is.null(paths$aux_dir)){
       sendSweetAlert(
         session = session,
         title = "Error !!",
         text = "Please, write a Crop name and select a root folder from your computer.",
         type = "error"
       )
       system.time(1)
       updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
     }
     
     withBusyIndicatorServer("load_rast",{ 
     
     files_in  <- input$aux_files$datapath#probar como importar varios archivos
     files_nms <- input$aux_files$name
     
     purrr::map2(.x = files_in, .y = files_nms , function(.x, .y){
       c_mask <- terra::rast(paths$mask_path)
       
       r_in <- terra::rast(.x)
       
       if(all(terra::res(r_in) == terra::res(c_mask))){
         r_in <- r_in %>% 
           terra::crop(., terra::ext(c_mask)) %>% 
           terra::mask(., c_mask)
       }else{
         r_in <- r_in %>% 
           terra::crop(., terra::ext(c_mask)) %>% 
           terra::resample(., c_mask) %>% 
           terra::mask(., c_mask)
       }
       
       terra::writeRaster(r_in, paste0(paths$aux_dir,"/", .y))
       
       
     })
     
     
     }) 
     updateButton(session, "load_rast",label = "Files loaded",style = "success")
   },
   error = function(e) {
     # return a safeError if a parsing error occurs
     #stop(safeError(e))
     sendSweetAlert(
       session = session,
       title = "Error !!",
       text = e,
       type = "error"
     )
   })
   
   
   
 })

 
 
#***********************************
#### Passport data processing ######
#**********************************
 

  # change tab when data is loaded
 bd <- reactive({
    req(input$data_in)
    tryCatch(
      {
        df.raw <- read.csv(input$data_in$datapath, header = TRUE, sep = ",", na.strings = c("","NA")) 
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        sendSweetAlert(
          session = session,
          title = "Error when loading database!!",
          text = e,
          type = "error"
        )
        return(NA)
      }
    )
    return(df.raw)
  
  })
  
  
 output$data_prev <- DT::renderDataTable({
   req(bd())
   
  DT::datatable(bd() %>%
               dplyr::slice(., 1:20), options = list(pageLength =5, scrollX = TRUE, searching = FALSE)) %>%
     DT::formatStyle(columns = input$col_number, backgroundColor  = "#B8D7F7")

 })
 
 output$groups_count <- renderValueBox({
  req(bd())
   valueBox(
     paste0(bd() %>% dplyr::pull(input$col_number) %>% unique %>% length(.)), "Groups", icon = icon("object-group","fas fa", verify_fa = FALSE),
     color = "purple"
   )
 }) 
 output$total_records <- renderValueBox({
   req(bd())
   valueBox(
     paste0(bd() %>% nrow(.)), "Rows", icon = icon("bars","fas fa", verify_fa = FALSE)
     
   )
 }) 
 output$na_percent <- renderValueBox({
   x <- (bd() %>% dplyr::pull(input$col_number) %>% is.na(.) %>% sum(.))/ nrow(bd())
   x <- round(x,3)*100
   color <- ifelse(x <= 30, "green", ifelse(x<= 75, "yellow", ifelse(x > 75, "red", "NA")))
   valueBox(
     paste0(x, "%" ), "Missing Values", icon = icon("percent","fas fa", verify_fa = FALSE),
     color = color
   )
 }) 
 

 observeEvent(input$data_in,{
   updateTabsetPanel(session, inputId = "tab_passport", selected = "Preview data" )
 })
 
 
#### boton para procesar la base de datos
  observeEvent(input$prepare_data, {
    req(input$data_in)
  
    #prevent errors when no crop name and race are defined by the user
    if(is.null(paths$crop_name) | is.null(paths$input_data_dir)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, write a Crop name and select a root folder from your computer.",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
    }else if(is.null(paths$mask_path)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, import a valid region raster mask..",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Geographic area" )
    }else{
  
  clas_res <- list()
   if(file.exists( paths$original_path )){file.remove( paths$original_path)}
   if(file.exists( paths$db_file)){file.remove(paths$db_file)}
   if(file.exists( paths$db_summ)){file.remove(paths$db_summ)}
    
   withBusyIndicatorServer("prepare_data", {
     
     
     
     source(paste("www/scripts/01_classification/classification_function_new_models.R", sep = ""), local = F)
     source(paste("www/scripts/01_classification/prepare_input_data.R", sep = ""), local = F)
     
    tryCatch({
        
      
       data_out <- prepare_input_data(data = bd(),#data_path = input$data_in$datapath,
                                      col_number = input$col_number,
                                      mask = paths$mask_path,
                                      env_paths = paths$generic_dir,
                                      aux_paths = paths$aux_dir)
       
       
      #guardar archivo de datos originale
       
       resources$n_invalid <- data_out$full_data %>% 
         dplyr::filter(invalid) %>% 
         nrow(.)
       
      resources$n_rows <- nrow(data_out$full_data)
       
      write.csv(data_out$full_data, paths$original_path, row.names = F)

     
      
      resources$cleaned_data    <- data_out$full_data %>% 
        dplyr::filter(!invalid) %>% 
        dplyr::select(Y,
                      everything(.),
                      -any_of(data_out$names_out),
                      matches("^status"), 
                      matches("^database_id"), 
                      matches("^source_db"),
                      -invalid) 
      
      
      predicted <- rep(FALSE, nrow(resources$cleaned_data))
      
      
      if(any(is.na(resources$cleaned_data$Y)) & !input$do_preds){
        shiny::safeError("Missing values found in response variable!!. Please, check do prediction parameter")
        sendSweetAlert(
          session = session,
          title = "Warning!",
          text = "Missing values found in response variable!!. Please, check do prediction parameter",
          type = "warning"
        )
        
      }else if(!any(is.na(resources$cleaned_data$Y)) & input$do_preds){
        sendSweetAlert(
          session = session,
          title = "Warning!",
          text = "No Missing values found in response variable. Can not train ML models to predict.",
          type = "warning"
        )
      }else if(any(is.na(resources$cleaned_data$Y)) & input$do_preds){
        
        
          var_to_exclude <- c("drymonths_2_5_min", "ethnicity", "monthCountByTemp10")
        
          
          
          to_model <- data_out$full_data %>%
            dplyr::filter(!is.na(Y)) %>% 
            dplyr::select(-all_of(data_out$names_out))
          
          #asignar clase con mayor freq a duplicados
          
          get_predominant_class <- function(k){
            names(which.max(table(k)))
          }
          new_labs <- to_model %>% 
            dplyr::group_by(Latitude, Longitude) %>% 
            dplyr::summarise(lab =  get_predominant_class(Y)) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(key = paste0(Latitude,";", Longitude)) %>% 
            dplyr::select(-Latitude, -Longitude)
          
          # i <- to_test[184, c("Latitude", "Longitude")]
          # i <- as.numeric(i)
          #   to_model %>% 
          #     dplyr::filter(., !is.na(Y)) %>% 
          #     dplyr::filter(., Latitude == i[1] & Longitude == i[2] ) %>% 
          #     dplyr::pull(., Y) %>%
          #     table()
           
          to_model<- to_model %>% 
            dplyr::mutate(key = paste0(Latitude,";", Longitude)) %>% 
            dplyr::left_join(., new_labs, by = c("key")) %>% 
            dplyr::mutate(Y = lab) %>% 
            dplyr::filter(!invalid) %>% 
            dplyr::select(Y,
                        -Latitude,
                        -Longitude,
                        -matches("^status"), 
                        -matches("^database_id"), 
                        -matches("^source_db"),
                        where(is.numeric),
                        -any_of(var_to_exclude))
          
          
          if(is.factor(to_model$Y)){
            to_model$Y <- droplevels(to_model$Y)
          }
          
          if(any(table(na.omit(to_model$Y)) < 20)){
            sendSweetAlert(
              session = session,
              title = "Unbalanced groups!!",
              text = "Class with less than 20 records. can not fit ML models.",
              type = "warning"
            )
          }else{
            
            to_predict <- resources$cleaned_data  %>% 
              dplyr::filter(is.na(Y)) %>% 
              dplyr::select(Y,
                            -Latitude,
                            -Longitude,
                            -matches("^status"), 
                            -matches("^database_id"), 
                            -matches("^source_db"),
                            where(is.numeric),
                            -any_of(var_to_exclude)) 
            

            clas_res <- classification_fun(df              = to_model %>%
                                           dplyr::mutate(Y = factor(Y)),
                                           standardize_all = T,
                                           sampling_mthd   = "none",
                                           omit_correlated = T,
                                           top_variables   = 5,
                                           external_df     = to_predict%>%
                                             dplyr::mutate(Y = factor(Y)))

            
            predicted[is.na(resources$cleaned_data$Y)] <- TRUE
            resources$cleaned_data$Y[is.na(resources$cleaned_data$Y)] <-  clas_res$External_data_predictions %>%
              dplyr::pull(., ensemble) %>%
              as.character()
            

          }

        
      }
      
      resources$cleaned_data  <- resources$cleaned_data %>% 
        dplyr::mutate(predicted = predicted)
      
      resources$group_names <- names(table(resources$cleaned_data$Y))
        
      write.csv(resources$cleaned_data,   paths$db_file, row.names = F)
      
      
      
      clas_res$invalid <- resources$n_invalid
      clas_res$n_rows  <-resources$n_rows
      
      saveRDS(clas_res, paths$db_summ)
      
      updateButton(session, "prepare_data",label = "Done",style = "success")
      updateTabsetPanel(session, inputId = "tab_passport", selected = "Results" )
      
      
    }, error = function(e){
      
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = paste("Cleaning database process failed \n ", e),
        type = "error"
      )
      
    })
    
    })
    
    }#end else
    
    
  })#end observeEvent clean data
  
  observeEvent(input$tab_passport, {
    req(paths$db_file)
    if(file.exists(paths$db_file) & file.exists(paths$db_summ)){
      
      clas_res <- readRDS(paths$db_summ)
      
      resources$n_invalid <- clas_res$invalid
      resources$n_rows <- clas_res$n_rows
      
      
      output$data_out <-  renderDT({
        req(paths$db_file)
        if(file.exists(paths$db_file)){
          resources$cleaned_data <- read.csv(paths$db_file, header =T, stringsAsFactors = F)
        }
        
        datatable(resources$cleaned_data, options = list(pageLength =5, scrollX = TRUE, searching = FALSE)) 
      })
      
      
      output$infbox <- renderUI({
        
        list(
          fluidRow(
            valueBox(value =  round((resources$n_invalid/resources$n_rows)*100, 1),
                     subtitle = "% Not Valid records",
                     icon     = icon("percent","fas fa", verify_fa = FALSE),
                     color    = "aqua",
                     width    = 15)
          ),
          fluidRow(
            valueBox( value    = format(resources$n_invalid, nsmall=0, big.mark=".", decimal.mark = ","),
                      subtitle = "Total Not valid records",
                      icon     = icon("square-xmark","fas fa", verify_fa = FALSE),
                      color    = "light-blue",
                      width    = 15) 
          ),
          fluidRow(
            valueBox( value    = format(resources$n_rows, nsmall=0, big.mark=".", decimal.mark = ","),
                      subtitle = "Total records",
                      icon     = icon("bars","fas fa", verify_fa = FALSE),
                      color    = "olive",
                      width    = 15) 
          )
        )
        
        
        
      })
      
      output$gchart1 <- renderGvis({
        gvisPieChart(resources$cleaned_data %>% 
                       dplyr::group_by(Y) %>% count(),
                     options = list(is3D = "true",
                                    width=400,
                                    height=300,
                                    title = "Groups counts"))
      })
      
      output$lmap1 <- renderLeaflet({
        
        dt <- resources$cleaned_data %>% 
          dplyr::select(Y, Latitude, Longitude) 
        
        pp<- dt
        coordinates(pp) <- ~Latitude+Longitude
        crs(pp) <- coor_sys
        bbx <- pp@bbox
        
        cent <- as.numeric(rowSums(bbx)/2)
        pal <- colorFactor(palette = "Set1", domain = unique(dt$Y))
        
        leaflet("lmap1") %>% 
          setView(lat= cent[1], lng = cent[2], zoom = 3) %>% 
          addTiles(options = providerTileOptions(noWrap = TRUE) )%>% 
          addCircles(data = dt, 
                     radius =  ~rep(7000, nrow(dt)),
                     color = ~pal(Y), 
                     stroke = F, 
                     fillOpacity = 0.8,
                     label = ~as.character(Y)) %>%
          addLegend(data = dt, 
                    position = "bottomright", 
                    colors = ~pal(unique(Y)), 
                    labels = ~factor(unique(Y)),
                    title = "Class names",
                    opacity = 0.8)
      })
      
      output$pca_res <- renderDT({
        
        if(!is.null(clas_res$Testing_CM$ensemble)){
          to_print <- round(clas_res$Testing_CM$ensemble[[4]], 2) %>% 
            as_tibble(., rownames = "Group")
          DT::datatable(to_print, options = list(pageLength =5, scrollX = TRUE, searching = FALSE))
          
        }else{
          safeError("No ML results file found.")
        }
        
      }, res = 100)
      
      
      
    }
   
    
  })

  
  #******************************
  ######### SDM MODELLING ########
  #******************************
  
  ## check if occ file exists
  
  output$groups_picker1 <- renderUI({
  
    if(is.null(resources$group_names)){
      chois <- "none"
    }else{
      chois <- resources$group_names
    }
    tagList(
      pickerInput(
        inputId = "select_group",
        label = tags$h4(tags$b("Select Group/Class/Race to process: ")), 
        choices = chois,
        options = list(
          style = "btn-primary")
      )
      
      
    )
    
    
  })
 
  ##***************************************##
  ##**** SDM modelling *****************##
  ##**************************************##
  
  
  
  observeEvent(input$run_sdm, {
    
    if(!file.exists(paths$maxent_jar)){
      file.copy(from =  "www/maxent.jar", 
                to = paths$maxent_jar)
      
    }
    
    
    #req(resources$pseudo_abs, resources$spData)
    
    if(is.null(paths$crop_name) | is.null(paths$input_data_dir)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, write a Crop name and select a root folder from your computer.",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
    }else if(is.null(paths$mask_path)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, import a valid region raster mask..",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Geographic area" )
    }else{
      
        source("www/scripts/01_classification/create_occ_shp.R")
        source("www/scripts/02_sdm_modeling/background_points.R")
        source("www/scripts/02_sdm_modeling/calibration_function.R")
        source("www/scripts/02_sdm_modeling/tuning_maxNet.R")
        source("www/scripts/02_sdm_modeling/sdm_maxent_java_approach_function.R")
        source("www/scripts/02_sdm_modeling/pseudo_occ.R")
      
      paths$occName <- as.character(input$select_group)
      
      
      
      withBusyIndicatorServer("run_sdm", { 
        
        if(paths$occName != "none"){
          paths$sp_results <- paste(paths$results_dir, paths$occName, sep = .Platform$file.sep);if(!file.exists(paths$sp_results)){dir.create(paths$sp_results, recursive = TRUE)}
          paths$gap_outDir <- paste0(paths$sp_results, "/gap_scores");if(!file.exists(paths$gap_outDir)){dir.create(paths$gap_outDir, recursive = TRUE)}
          paths$delaDir    <- paste0(paths$gap_outDir, "/delaunay");if(!file.exists(paths$delaDir)){dir.create(paths$delaDir, recursive = TRUE)}
          paths$sdm_outDir <- paste0(paths$sp_results, "/species_distribution");if(!file.exists(paths$sdm_outDir)){dir.create(paths$sdm_outDir, recursive = TRUE)}
          paths$pseudoDir  <- paste0(paths$sdm_outDir, "/pseudo_absences");if(!file.exists(paths$pseudoDir)){dir.create(paths$pseudoDir, recursive = TRUE)}
          paths$cvfoldsDir <- paste0(paths$sdm_outDir, "/sdm_cv_folds");if(!file.exists(paths$cvfoldsDir)){dir.create(paths$cvfoldsDir, recursive = TRUE)}
          paths$gap_valDir <- paste0(paths$sp_results, "/gaps_validation");if(!file.exists(paths$gap_valDir)){dir.create(paths$gap_valDir, recursive = TRUE)}
          paths$occ_group  <- paste0(paths$sp_results, "/", paths$occName,"_occurrences.csv")
          paths$occ_ids    <- paste0(paths$sp_results, "/", paths$occName,"_occurrences_IDs.csv")
          paths$occ_shp    <- paste0(paths$sdm_outDir, "/occurrences.shp")
          paths$occ_csv    <- paste0(paths$sdm_outDir, "/occurrences.csv")
          paths$pseudo_file<- paste0(paths$pseudoDir,  "/pseudo_absences_", paths$occName, ".csv")
          paths$smd_var_selected <- paste0(paths$sdm_outDir, "/sdm_variables_selected.csv")
          paths$sdm_calibration  <- paste0(paths$sdm_outDir, "/tunned_sdm_parms.csv")
        }
        
        
        #reso  <<- reactiveValuesToList(resources)
        #pths  <<- reactiveValuesToList(paths)
        create_occ_shp(data         = resources$cleaned_data,
                       file_output  = paths$occ_csv,
                       shp_output   = paths$occ_shp,
                       validation   = FALSE,
                       mask_pth     = paths$mask_path,
                       occName      = paths$occName)
        
        print(paths$mask_path)
        
        
        resources$spData <- resources$cleaned_data %>% 
          dplyr::filter(Y == paths$occName) %>% 
          dplyr::select(-any_of(c("predicted", "database_id", "source_db")))
        
        ##Generar pseudo-presencias cunado n es pequeño
        
        nrows <- nrow(resources$spData)
        
        tms <- dplyr::case_when(
          nrows <= 30 ~ 10,
          nrows > 30 & nrows <= 50 ~ 5,
          nrows >50 & nrows <=80 ~ 3,
          .default = 1
        )
          
      
          resources$spData <- pseudo_occ(xy = resources$spData, 
                                         clim_dir = paths$generic_dir, 
                                         aux_clim_dir = paths$aux_dir ,
                                         tms  = tms,
                                         nu = 0.05)
          
        
        resources$pseudo_abs <- pseudoAbsences_generator(data        = resources$spData,
                                                         climDir     = paths$generic_dir, 
                                                         aux_dir     = paths$aux_dir,
                                                         clsModel    = "Y", 
                                                         overwrite   = F, 
                                                         occName     = paths$occName,
                                                         correlation = 3, 
                                                         pa_method   = "ecoreg",
                                                         bg_out_path = paths$pseudo_file,
                                                         mask_path   = paths$mask_path,
                                                         smd_var_selected_path = paths$smd_var_selected,
                                                         ecoreg_path = "www/masks/World_ELU_2015_5km.tif")
        #spDatax <<- resources$cleaned_data
        #pseudo_abs <<- resources$pseudo_abs
        
        resources$var_names <- read.csv(paths$smd_var_selected, stringsAsFactors = F) %>% 
          dplyr::pull(x)
        
       
        
        write.csv(resources$spData, paths$occ_group ,row.names = F)
        
        resources$cleaned_data %>% 
          dplyr::filter(Y == paths$occName) %>% 
          write.csv(.,
                    paths$occ_ids,
                    row.names = F)
        
        spData <<- resources$spData %>% 
          dplyr::select(-any_of("status"))
                        
        bg_data <<- resources$pseudo_abs
        #if(input$calib == 1){
          ##no olvidar cambiar el numvero de betas a testear en TrainMaxnet function lo reduje para probar
          params_tunned  <- Calibration_function(spData   = resources$spData %>% 
                                                   dplyr::select(-any_of("status")),
                                                 bg_data  = resources$pseudo_abs,
                                                 occName  = paths$occName,
                                                 out_path = paths$sdm_calibration,
                                                 ommit = F,
                                                 use.maxnet = TRUE)

        #}else{
         # params_tunned <- list()
        #  params_tunned$beta <- input$betamp
         # params_tunned$features <- input$feats
          
      #  }
        
        if(file.exists(paths$smd_var_selected)){
          resources$var_names <- read.csv(paths$smd_var_selected, stringsAsFactors = F) %>% 
            dplyr::pull(x)
        }else{
          #safeError("No pseudo-abs file is avialiable.")
          sendSweetAlert(
            session = session,
            title = "Error !!",
            text = "Psuedo-absences file is missing.",
            type = "error"
          )
          system.time(1)
          updateNavbarPage(session, inputId = "nvpage_tab2", selected = "Pseudo-absences" )
        }
        
        
       
       
        sdm_maxent_approach_function(occName      = paths$occName,
                                     spData       = resources$spData%>% 
                                       dplyr::select(-any_of("status")),
                                     bg_data      = resources$pseudo_abs,
                                     var_names    = resources$var_names,
                                     model_outDir = paths$sdm_outDir,
                                     replic_path  = paths$cvfoldsDir,
                                     nFolds       = input$nfolds,
                                     climDir      = paths$generic_dir,
                                     clim_spReg   = paths$aux_dir,
                                     beta         = params_tunned$beta,
                                     feat         = params_tunned$features,
                                     doSDraster   = FALSE,
                                     varImp       = FALSE,
                                     validation   = FALSE)
        
        
        
      })
        #paths$sdm_occ_path
        
    }
    

    updateButton(session, "run_sdm",label = "Done",style = "success")
    updateTabsetPanel(session, inputId = "sdm_res", selected = "Results" )
    
  })

  
  output$map3 <- renderLeaflet({
    
    mp <- leaflet() %>% 
      addTiles()
    if(!is.null(paths$sdm_occ_path) & file.exists(paths$sdm_occ_path)){
      tmp_r <- terra::rast(paths$sdm_occ_path)
      
      
      pal <- colorNumeric("RdYlBu", terra::values(tmp_r),
                          na.color = "transparent", 
                          reverse = TRUE)
      
      mp<- mp %>% 
        addRasterImage(x =  tmp_r, colors = pal) %>%
        addLegend(pal = pal, 
                  values = terra::values(tmp_r),
                  title = "Probability map", 
                  group = "sdm", 
                  position = "bottomleft")
      
    }
    
    
  })
  
#**********************************  
##### GAP SCORES              ####
#********************************

  
  observeEvent(input$calc_all_gps,{
    
    if(is.null(paths$crop_name) | is.null(paths$input_data_dir)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, write a Crop name and select a root folder from your computer.",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Working directory" )
    }else if(is.null(paths$mask_path)){
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Please, import a valid region raster mask..",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab1", selected = "Geographic area" )
    }else if(is.null(paths$occName)){
      
      sendSweetAlert(
        session = session,
        title = "Error !!",
        text = "Group/race/class not selected.",
        type = "error"
      )
      system.time(1)
      updateNavbarPage(session, inputId = "nvpage_tab2", selected = "Pseudo-absences" )
      
    }else{
      
      withBusyIndicatorServer("calc_all_gps", { 
        
        print(">>> Calculating Accessibility geoscore: ")
        
        source("www/scripts/03_gap_methods/cost_distance_function.R")
        
        tryCatch(expr = {
          
          cost_dist_function(
            cost_out_path = paths$cost_out_path,
            friction      = "www/masks/friction_surface.tif",
            mask          = paths$mask_path,
            Occ           = resources$spData,
            sdm_path      = paths$sdm_occ_path
          )
          

          print(">>> Accessibility Done")
          
          print(">>> Calculating Connectivity geoscore: ")
          
          source("www/scripts/03_gap_methods/delaunay.R")
          source("www/scripts/03_gap_methods/delaunay_geo_score.R")
          
          
          calc_delaunay_score(
            coreDir = paths$sp_results, 
            ncores = NULL, 
            validation = FALSE, 
            pnt = NULL)
          
          print(">>> Connectivity Done")
          
          
          print(">>> Calculating Environmental geoscore: ")
          source("www/scripts/03_gap_methods/ecogeo_cluster.R")
          source("www/scripts/03_gap_methods/env_distance.R")
          
          calc_env_score(sdm_path = paths$sdm_occ_path, 
                         clus_method = "hclust_mahalanobis", 
                         gap_dir   = paths$gap_outDir, 
                         occ_dir   = paths$sdm_outDir, 
                         env_dir   = paths$generic_dir, 
                         var_names = resources$var_names,
                         n.sample  = input$nsample,
                         n.clust   = input$nclust)
          
          print(">>> Environmental Done")
          
          updateButton(session, "calc_all_gps",label = "Done",style = "success")
          updateTabsetPanel(session, inputId = "all_res", selected = "Results" )
          
          
          
        }, error = function(e){
          base::message(e)
          resources$check_gap_scores <- "error"
          
        })
        
        # in case gap score process fails THEN 
        if(!is.null(resources$check_gap_scores)){

          sendSweetAlert(
            session = session,
            title = "Error !!",
            text = "Gap scores process Failed",
            type = "error"
          )
          
          source("www/scripts/03_gap_methods/CWR_gap_map.R")
          
          resources$final_gap_map <- CWR_gap_map(sdm_path = paths$sdm_occ_path, 
                                                 occ_shp  = paths$occ_shp, 
                                                 out_file = paths$CWR_gap_map_path)
          
          resources$covg_ul <- 100 - length(which(resources$final_gap_map[]==3))/length(which(!is.na(resources$final_gap_map[])))*100
          resources$covg_ll <- 100 - length(which(resources$final_gap_map[]>=2))/length(which(!is.na(resources$final_gap_map[])))*100
          
        }
        
        
        
        
      })
      
      
      
    }
    
  })
  
  
  output$map4 <- renderLeaflet({
    req(paths$cost_out_path)

        mp <- leaflet() %>% 
      addTiles()
    
    if(!is.null(paths$cost_out_path) & file.exists(paths$cost_out_path)){
      
      tmp_r <- terra::rast(paths$cost_out_path)
      
      
      pal <- colorNumeric("magma", terra::values(tmp_r),
                          na.color = "transparent", 
                          reverse = F)
      
      mp<- mp %>% 
        addRasterImage(x =  tmp_r, colors = pal) %>%
        addLegend(pal = pal, 
                  values = terra::values(tmp_r),
                  title = "Cost distance", 
                  group = "cost", 
                  position = "bottomleft")
      
      
    }
    
  })

#################################
###### delaunay geo score ######
###############################
  
  
  output$map5 <- renderLeaflet({
    req(paths$dela_out_path)
    mp <- leaflet() %>% 
      addTiles()
    
    if(!is.null(paths$dela_out_path) & file.exists(paths$dela_out_path)){
      
      tmp_r <- terra::rast(paths$dela_out_path)
      
      
      pal <- colorNumeric("Greens", terra::values(tmp_r),
                          na.color = "transparent", 
                          reverse = F)
      
      mp<- mp %>% 
        addRasterImage(x =  tmp_r, colors = pal) %>%
        addLegend(pal = pal, 
                  values = terra::values(tmp_r),
                  title = "Delaunay geo score", 
                  group = "dela", 
                  position = "bottomleft")
      
      
    }
    
  })
  
#################################
##### environmental score ######
###############################
  
 
 output$map6 <- renderLeaflet({
   req(paths$envi_out_path)
   mp <- leaflet() %>% 
     addTiles()
   
   if(!is.null(paths$envi_out_path) & file.exists(paths$envi_out_path)){
     
     tmp_r <- terra::rast(paths$envi_out_path)
     
     
     pal <- colorNumeric("Greens", terra::values(tmp_r),
                         na.color = "transparent", 
                         reverse = F)
     
     mp<- mp %>% 
       addRasterImage(x =  tmp_r, colors = pal) %>%
       addLegend(pal = pal, 
                 values = terra::values(tmp_r),
                 title = "Environmental geo score", 
                 group = "dela", 
                 position = "bottomleft")
     
     
   }
   
 })
  

 
 ############################################
 ####### validation function ##############
 ########################################
 
 
 observeEvent(input$do_validation, {
   
   if(is.null(paths$occName)){
     
     sendSweetAlert(
       session = session,
       title = "Error !!",
       text = "Group/race/class not selected.",
       type = "error"
     )
     system.time(1)
     updateNavbarPage(session, inputId = "nvpage_tab2", selected = "Pseudo-absences" )
     
   }else{
     
     if(!is.null(resources$check_gap_scores)){
       # If gap score are not available THEN
       
       if(file.exists(paths$CWR_gap_map_path)){
         
         source("www/scripts/04_others/03_generate_report_CWR.R")
         generate_report_CWR(out_file = paths$out_report_file, 
                             paths = reactiveValuesToList(paths)  , 
                             resources = reactiveValuesToList(resources)  )
         
         updateButton(session, "do_validation",label = "Done",style = "success")
         updateTabsetPanel(session, inputId = "val_res", selected = "Results" )
         
         
       }else{
         sendSweetAlert(
           session = session,
           title = "Error !!",
           text = "Final Gap Map not Found!",
           type = "error"
         )
       }
       
       
       
     }else{
       
       withBusyIndicatorServer("do_validation", {
        
         
         tryCatch({
           
           
           #aqui hay un error que no he podido arreglar
           source("www/scripts/03_gap_methods/validation_workflow.R", local = TRUE)
           
           ####  GENERATE REPORT
           generate_report(out_file = paths$out_report_file,
                           paths = reactiveValuesToList(paths)  ,
                           resources = reactiveValuesToList(resources)  )

           
         }, error = function(e){
          
           print(paste0("We got an error in validation: ", e))
           resources$check_validation <- "error"
           
         })
         
         
         if(!is.null(resources$check_validation )){
           
           sendSweetAlert(
             session = session,
             title = "Error !!",
             text = "Cannot perfrom validation, low number of accessions!",
             type = "error"
           )
           
           cat("Classifying gap scores using 0.5 quantile \n")
           
           source("www/scripts/04_others/reclass_by_median.R")
           source("www/scripts/04_others/03_generate_report_CWR.R")

           if(file.exists(paths$cost_out_path) & file.exists(paths$dela_out_path ) & file.exists(paths$envi_out_path) ){
             

             q = 0.5
             cost <- terra::rast(paths$cost_out_path) %>% 
               reclass_by_median(., q = q)
             dela <- terra::rast(paths$dela_out_path)%>% 
               reclass_by_median(., q = q)
             envi <- terra::rast(paths$envi_out_path)%>% 
               reclass_by_median(., q = q)
             
             resources$final_gap_map <- sum(cost, dela, envi, na.rm = T)
             resources$final_gap_map[resources$final_gap_map == 0] <- NA
             #levels(resources$final_gap_map) <- data.frame(ID = 1:3, clase = c("Low", "Medium", "high"))
             
             terra::writeRaster(resources$final_gap_map, paste0(paths$sp_results, "/", paths$occName,"_final_gap_map.tif"), overwrite = T)
             
             
             rs_gaps <- resources$final_gap_map #raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/mesoamerican/americas/gap_models/gap_class_final_new.tif")
             resources$covg_ul <- 100 - length(which(rs_gaps[]==3))/length(which(!is.na(rs_gaps[])))*100
             resources$covg_ll <- 100 - length(which(rs_gaps[]>=2))/length(which(!is.na(rs_gaps[])))*100
             
             generate_report_CWR(out_file = paths$out_report_file, 
                                 paths = reactiveValuesToList(paths)  , 
                                 resources = reactiveValuesToList(resources)  )
             
             #plet(gap_map, col = c("1" = "lightgrey", "2"="blue","3"= "red"))
             
           }
           
           
         }
        
         
       })
     }
   }
   
   updateButton(session, "do_validation",label = "Done",style = "success")
   updateTabsetPanel(session, inputId = "val_res", selected = "Results" )
   
 })
 
 
 output$map7 <- renderLeaflet({
   
   mp <- leaflet() %>% 
     addTiles()
   
   if(file.exists(paths$CWR_gap_map_path)){
     
     tmp_r <- terra::rast(paths$CWR_gap_map_path)
     
     mp<- mp %>% 
       addRasterImage(x =  tmp_r, colors = c("#d9d9d9", "#80bfff", "#ff471a" ), opacity = 0.7) %>%
       addLegend(colors = c("#d9d9d9", "#80bfff", "#ff471a" ),
                 labels = c("Low pro. Gap", "Medium prob. Gap", "High prob. Gap"),
                 values = terra::values(tmp_r),
                 title = "Delaunay geo score", 
                 group = "dela", 
                 position = "bottomleft")
     
   }else{
    
     if(file.exists(paths$final_gap_map_path)){
       
       tmp_r <- terra::rast(paths$final_gap_map_path)
       
     }
     
     mp<- mp %>% 
       addRasterImage(x =  tmp_r, colors = c("#d9d9d9", "#80bfff", "#ff471a" ), opacity = 0.7) %>%
       addLegend(colors = c("#d9d9d9", "#80bfff", "#ff471a" ),
                 labels = c("Low pro. Gap", "Medium prob. Gap", "High prob. Gap"),
                 values = terra::values(tmp_r),
                 title = "Delaunay geo score", 
                 group = "dela", 
                 position = "bottomleft")
     
   }
  
   
   
   
   mp
   
 })
 
 

 df_aucs <- eventReactive(paths$gap_valDir, {
   
   if(file.exists(paste0(paths$gap_valDir,"/simulation_results.xlsx"))){
     
     to_ret <- list()
     pth  <- paste0(paths$gap_valDir,"/simulation_results.xlsx")
     shts <- readxl::excel_sheets(pth)
     
     overall_simulation_metrics <-  lapply(shts, function(sht){
       
       if(grepl("overall", sht)){
         fls <- readxl::read_excel(pth, sheet = sht) %>% as.data.frame()
         to_ret <- data.frame(gap_score = fls[1,1],
                              auc.mean = fls[1,2],
                              lower.ic = fls[2,2],
                              upper.ic = fls[3,2]
         )
       }else{
         
         fls <- readxl::read_excel(pth, sheet = sht)
         
         nm <- fls[1,2]
         to_ret <- data.frame(gap_score = nm, fls[nrow(fls), c(4,6,7)])
       }
       
       return(to_ret)
     }) %>%
       do.call(rbind , .)
     
     overall_simulation_metrics$gap_score <- factor(overall_simulation_metrics$gap_score, levels = rev(overall_simulation_metrics$gap_score))
     
     to_ret$aucs <-overall_simulation_metrics
     
     covs <- readxl::read_excel(pth, "overall_simulation_metrics" )
     covs <- data.frame(Metric = "Coverage", average_coverage = mean(covs$Value[(length(covs$Value)-1):length(covs$Value)]),
                        upper = covs$Value[(length(covs$Value)-1)],
                        lower = covs$Value[length(covs$Value)])
     
     to_ret$covs <- covs
     
    
     return(to_ret)
     
   }
   
   })
 
 output$plt1  <- renderPlot({
 
   
   ggplot(data= df_aucs()$aucs, aes(y=gap_score, x=auc.mean, xmin=lower.ic, xmax=upper.ic)) +
     geom_point() +
     geom_errorbarh(height=.1) +
     xlab("Area Under ROC Cuve")+
     ylab("Gap score")+
     geom_label(aes(y = gap_score, x = auc.mean, label = auc.mean, vjust = -0.5))+
     labs(title  = "Area Under ROC Curve (AUC)")+
     theme_bw()
   
  
   
 })
 
 
 output$plt2  <- renderPlot({ 
   
   ggplot(data= df_aucs()$covs, aes(y=Metric, x=average_coverage, xmin=lower, xmax=upper)) +
     geom_point() +
     geom_errorbarh(height=.1) +
     xlab("coverage of existing germplasm collections")+
     ylab("")+
     geom_label(aes(y = Metric, x = average_coverage, label = average_coverage, vjust = -0.5))+
     labs(title  = "Coverage (%)")+
     theme_bw()
   
  
   
   })

 

 output$renderReport <- renderUI({
   
   req(paths$out_report_file)
   
   if(file.exists(paths$out_report_file)){
     pth = gsub(".Rmd", ".html", paths$out_report_file)
     htmltools::includeHTML(pth)
     #tags$iframe(srcdoc = readLines(pth), scrolling = "yes", width = "100%", height = "600px")
     
    
   }
 })
 

 
}#end server function everything











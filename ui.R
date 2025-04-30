#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#library(pacman)
suppressMessages(if(!require(pacman)){install.packages("pacman");library(pacman)}else{library(pacman)})

pacman::p_load(shiny, shinydashboard, leaflet, raster, sp, rsconnect, ggplot2, 
               shinyFiles, shinyBS, shinyjs, yaml, shinyWidgets, rmarkdown, bsplus, tidyverse, 
               shinydashboardPlus, DT )
# Define UI for application that draws a histogram
urls <- read.csv("www/downloadable_files.csv")
source("www/helpers.R", local = TRUE)


header <- shinydashboardPlus::dashboardHeader(
  title = "Gap analysis UI",
  dropdownMenuOutput("messageMenu"),
  dropdownMenuOutput("restoreSession")
  # tags$li(a(href = 'http://shinyapps.company.com',
  #           icon("power-off"),
  #           title = "Back to Apps Home"),
  #         class = "dropdown")
)

sidebar <- shinydashboardPlus::dashboardSidebar(
              sidebarMenu(id = "menu", 
                           #menuItem(" Introduction", tabName = "intro", icon = icon("sun", "far fa", verify_fa = FALSE), selected = T),
                           menuItem(" Wizards", tabName = "tab1", icon = icon("folder-tree")),
                           menuItem(" Gap scores", tabName = "tab2", icon = icon("chart-area", verify_fa = FALSE))
                           #menuItem(" SDM", tabName = "tab3", icon = icon("earth-americas","fas fa", verify_fa = FALSE))
              ) 
                           
                           )
body <- dashboardBody(uiOutput("modal1"),
  tabItems(
    tabItem(tabName = "intro", tags$html('<p> Nothing to see here now </p>') ,
            tags$head(HTML("<style> .well{
                        border: 1px solid rgb(45, 152, 214);
                        position: relative;
                        border-radius: 3px;
                        background: #ffffff;
                        border-top: 4px solid #3F89B8;
                        margin-bottom: 20px;
                        width: 100%;
                        box-shadow: 0 1px 1px rgba(0,0,0,.1);}

            .modal-content {
                           position: relative;
                           background-color: #fff;
                           -webkit-background-clip: padding-box;
                           background-clip: padding-box;
                           border: 1px solid #999;
                           border: 1px solid rgba(0,0,0,.2);
                           border-radius: 6px;
                           outline: 0;
                           -webkit-box-shadow: 0 3px 9px rgba(0,0,0,.5);
                           box-shadow: 0 3px 9px rgba(0,0,0,.5);
                           }
                 </style>"))
            
            
    ),
    tabItem(tabName = "tab1", 
            use_bs_tooltip(),
            use_bs_popover(),
            withMathJax(),
            navbarPage(title =  tags$i(class = "fa fa-gear",  style = "color: darkgray"),
                       windowTitle = "navp1",
                       id = "nvpage_tab1",
                       tabPanel(title = "Working directory",
                                sidebarLayout(
                                  sidebarPanel(width = 4, id = "write_crop_info",
                                              h3( tags$strong("Working directory set up")),
                                              textInput(inputId =  "selected.root.folder", label = "select folder", value = "No folder path selected")%>%
                                                short_info(input = ., place = "top",title = "Define the root folder path, where all inputs and outputs will be stored. Press \"Browse\" btn and select a folder from your computer."),
                                              shinyDirButton(id = "select_path_btn", label = "Browse", title = "Select root folder", buttonType = "default",
                                                             class = NULL, icon = NULL, style = NULL),
                                              uiOutput("crop_def"),
                                              conditionalPanel("input.set_crop_name == 'add new crop'",
                                                               textInput(inputId =  "set_crop_name2", label = "Write Crop name2") %>%
                                                                 short_info(input = ., place = "top",title = "Name of the major crop")),
                                               
                                               # textInput(inputId =  "set.level.name", label = "Please write Race name")%>%
                                               #   short_info(input = ., place = "top",title = "Name of the Race, genetic group, sub-group,class, etc.. to process."),
                                               
                                               div(id ="separator", style = "width:100px; height:15px"),
                                               bsButton("create_dirs", size="default",label = "Select", block = F, style="primary"),
                                               bsTooltip(id = "create_dirs", title = "Create directories", placement = "right", trigger = "hover")
                                               
                                  ),
                                  mainPanel(
                                    tags$head(HTML(
                                      '<style> 
                             .card {
                             box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);
                             padding: 8px;
                             border: 1px solid #ccc;
                             border-radius: 5px; /* 5px rounded corners */
                             max-width: 800px;
                             margin: 0 auto;
                             background-color:white;
                             }
                             </style>'
                                    )),
                                    div(id = "text_1", class= "card",
                                        shiny::includeMarkdown("Rmarkdown_files/arrange_dir_text.Rmd")
                                    )
                                   
                                    
                                  )#end main panel
                                  
                                  
                                )
                                
                       ),
                       tabPanel("Define study area",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3("Geographic area"),
                                    radioGroupButtons("choose_1", label = "", choices  = c("Create" = 1, "Import"= 2, "Select" = 3), status = "primary", justified = T ),
                                    conditionalPanel(condition = "input.choose_1 == '1'",
                                                     
                                                     selectInput("area_selector", label = "Select one region:", choices = c("World" = 0,
                                                                                                                            "America" = 19,
                                                                                                                            "Europe" = 150,
                                                                                                                            "Asia" = 142,
                                                                                                                            "Africa" = 2,
                                                                                                                            "Oceania" = 9,
                                                                                                                            "---Custom region---" = 8), selected = 0),
                                                     #haciendo el conditional panel
                                                     conditionalPanel( condition = "input.area_selector == 8",
                                                                       awesomeCheckboxGroup("chk_bx_gr", "Countries selected:", choices = "None", selected = "None", 
                                                                                           
                                                                                           status = "primary")
                                                                       
                                                     ),
                                                     textInput("mask_name", label = "Write file name")%>%
                                                       short_info(input = ., place = "top",title = "Name of the Geographical area to analyze."),
                                                     useShinyjs(),
                                                     withBusyIndicatorUI(
                                                       button = bsButton("create_mask", size="default",label = "Create", style="primary")
                                                       
                                                     ),
                                                     bsTooltip(id = "create_mask", title = "Create Raster Mask", placement = "left", trigger = "hover")),
                                    conditionalPanel(condition = "input.choose_1 == '2'",
                                                     fileInput("mask_import", "Import raster file:",multiple = FALSE,accept = c("image/tiff"), buttonLabel = icon("magnifying-glass","fa-solid", verify_fa = FALSE))%>%
                                                       short_info(input = ., place = "top",title = "If you have already created a raster mask, please import it."),
                                                     withBusyIndicatorUI(
                                                     bsButton("import_mask", size="default",label = "Import", block = F, style="primary")
                                                     )),
                                    conditionalPanel(condition = "input.choose_1 == '3'", 
                                                     uiOutput("avaliable_masks"),
                                                     withBusyIndicatorUI(
                                                       bsButton("select_mask", size="default",label = "Select", block = F, style="primary")
                                                     ))
                                                     
                                    
                                  ),
                                  mainPanel(
                                    div(id= "text_2", class = "card",
                                        shiny::includeMarkdown("Rmarkdown_files/mask_creator_text.Rmd")
                                    ),
                                    conditionalPanel(condition = "input.choose_1 == '1'",
                                                     leafletOutput("map_selector")
                                                     )
                                  )#end main panel
                                )#end sidebar layout
                                
                                
                       ),
                       tabPanel("Covariates",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3("Auxiliary files"),
                                    pickerInput(
                                      inputId = "quest_aux",
                                      label = "Add external/auxiliary raster files: ", 
                                      choices = c("No", "Yes"),
                                      options = list(
                                        style = "btn-primary")
                                    ),
                                    conditionalPanel("input.quest_aux == 'Yes'",
                                                     fileInput(inputId = "aux_files", label = "Select Raster files", 
                                                               accept = c("image/tiff"), 
                                                               multiple = TRUE, 
                                                               buttonLabel = icon("magnifying-glass","fa-solid", verify_fa = FALSE))
                                    ),
                                    
      
                                    #div(style="float:right;width:95px;background-color:lightblue")
                                        withBusyIndicatorUI(
                                          bsButton("load_rast", size="urldefault",label = "Load", block = F, style="primary")
                                        ),bsTooltip(id = "Load rasters", title = "Import auxiliary rasters", placement = "right", trigger = "hover"),
                                  
                                    verbatimTextOutput("log_output_1")
                                    
                                  ),
                                  mainPanel(
                                    div(id= "text_dw",  class = "card",
                                        shiny::includeMarkdown("Rmarkdown_files/download_rasters_text.Rmd")
                                    )
                                    
                                  )
                                )
                       ),
                       tabPanel("Passport data",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3("Data base set up"),
                                    tags$hr(),
                                    fileInput("data_in",
                                              label = "1. Select .cvs database",
                                              multiple = FALSE,
                                              buttonLabel = icon("magnifying-glass","fa-solid", verify_fa = FALSE),
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")
                                    )%>%
                                      short_info(input = ., place = "top",title = "Browse passport data file (only .CSV files are accepted)."),
                                    numericInput("col_number", label = "2. Write response column number", value = 1, min = 1, max= 999)%>%
                                      short_info(input = ., place = "top",title = "Column number of (race name, genetic group name, etc.) variale"),
                                    awesomeCheckbox(
                                      inputId = "do_preds",
                                      label = "Predict missing groups", 
                                      value = FALSE
                                    ),
                                    conditionalPanel("input.do_preds",
                                                     pickerInput(
                                                       inputId = "samp_mth",
                                                       label = "Sampling method: ", 
                                                       choices = c("none", "down", "up"),
                                                       selected = "none",
                                                       options = list(
                                                         style = "btn-primary")
                                                     )),
                                    withBusyIndicatorUI(
                                      bsButton("prepare_data", size="default",label = "Set up database", block = F, style="primary")
                                      
                                    )
                                  ),
                                  mainPanel(
                                    tabBox(id = "tab_passport",width = 12,
                                           tabPanel("Description",
                                                    div(id= "text_3", class = "card",
                                                        shiny::includeMarkdown("Rmarkdown_files/database_creator.Rmd")
                                                        
                                                    )
                                                    
                                                    
                                           ),
                                           tabPanel("Parameters",
                                                    div(id = "param_1", class = "card",
                                                        shiny::includeMarkdown("Rmarkdown_files/database_creator_parameters.Rmd")
                                                    )
                                                    
                                           ),
                                           tabPanel("Preview data",
                                                    fluidRow(
                                                    DT::dataTableOutput("data_prev"),
                                                    tags$hr(),
                                                    valueBoxOutput("total_records"),
                                                    valueBoxOutput("groups_count"),
                                                    valueBoxOutput("na_percent")
                                                    )
                                                    
                                           ),
                                           tabPanel("Results",
                                                    tags$div(id = "tbl_output", style = "height:800px",
                                                           box( title = "Data base output",
                                                                    status = "primary",
                                                                    width = NULL,
                                                                    solidHeader = FALSE,
                                                                    collapsible = T,
                                                                    closable = FALSE,
                                                                    enable_dropdown = FALSE,
                                                                    DT::DTOutput("data_out")
                                                           ),
                                                           tags$hr(),
                                                           box(title = "Database Summary",
                                                               width = NULL,
                                                               status= "primary",
                                                               accordion(
                                                                 id="acordion1",
                                                                 accordionItem(
                                                                   title = "Counts",
                                                                   status  = "success",
                                                                   solidHeader = FALSE,
                                                                   collapsed = TRUE,
                                                                   fluidRow(
                                                                         column(8, 
                                                                                htmlOutput("gchart1")  
                                                                              ),
                                                                       column(4,
                                                                              uiOutput("infbox")
                                                                              )
                                                                   )
                                                                   
                                                                 ),
                                                                 accordionItem(
                                                                   title = "Valid records Map",
                                                                   status  = "success",
                                                                   collapsed = TRUE,
                                                                   solidHeader = FALSE,
                                                                   leafletOutput("lmap1")
                                                                 ),
                                                                 accordionItem(
                            
                                                                           title = "ML performance metrics",
                                                                   status  = "success",
                                                                   solidHeader = FALSE,
                                                                   collapsed = TRUE,
                                                                   DTOutput("pca_res")
                                                                 )
                                                                )
                                                               )
                                                           )
                                                    
                                                    )#end tabPanel
                                    )#end tabbox
                                    
                                  )#end mainpanel
                                  
                                )#end sidebar layout
                       )
                       
            )#end navarpage
    ),
    tabItem(tabName = "tab2",
            navbarPage(title = icon("earth-americas","fas fa", verify_fa = FALSE),
                       windowTitle = "navp2",
                       id = "nvpage_tab2",
                       tabPanel("SDM modelling",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3(tags$strong("Landrace spatial distribution")),
                                    tags$hr(),
                                    uiOutput("groups_picker1"),
                                    tags$hr(),
                                    h3(tags$strong("Model settings")),
                                    sliderInput("nfolds", "Number of cross validation folds:",
                                                min = 5,
                                                max = 10,
                                                value = 5),
                                    tags$hr(),
                                    # radioButtons("dosdrast",label = "Calculate SD raster:",
                                    #              choices = c("Yes" = TRUE, "No" = FALSE)),
                                    # radioButtons("varimp",label = "Calculate Var importance:",
                                    #              choices = c("Yes" = TRUE, "No" = FALSE)),
                                    withBusyIndicatorUI(
                                      bsButton("run_sdm", size="default",label = "Run model", block = F, style="primary")
                                      
                                    )
                                    
                                    
                                  ),
                                  mainPanel(
                                    tabBox(id = "sdm_res",width = 12,
                                           tabPanel("Description",
                                                    shiny::includeMarkdown("Rmarkdown_files/sdm_modelling_text.Rmd")),
                                           tabPanel("Results",
                                                    leafletOutput("map3")
                                                    )
                                    )
                                    
                                  )
                                  
                                )
                         
                       ),
                       tabPanel("Geo-Scores",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3(tags$strong("Accessibility geo score")),
                                    pickerInput(
                                      inputId = "friction",
                                      label = "Select friction surface file:", 
                                      choices = c("friction.tif"),
                                      selected = 1,
                                      options = list(
                                        style = "btn-primary")
                                    ),
                                    tags$hr(),
                                    h3(tags$strong("Connectivity geo score")),
                                    pickerInput(
                                      inputId = "Occ_sel",
                                      label = "Select occurrence shp file:", 
                                      choices = c("occurreces.shp"),
                                      selected = 1,
                                      options = list(
                                        style = "btn-primary")
                                    ),
                                    tags$hr(),
                                    h3(tags$strong("Environmental geo score")),
                                    numericInput(inputId = "nclust",
                                                 label = "Number of environmental cluster:",
                                                 value = 15,
                                                 min =2,
                                                 max = 45,
                                                 step = 1),
                                    sliderInput(inputId = "nsample",
                                                label   = "Coordinates to sample: ",
                                                min     = 2000,
                                                max     = 20000,
                                                value   = 5000
                                    ),
                                    withBusyIndicatorUI(
                                      bsButton("calc_all_gps", size="default",label = "Calculate", block = F, style="primary")
                                      
                                    )
                                    
                                  ),
                                  mainPanel(
                                    tabBox(width = 12,
                                           id = "all_res",
                                           tabPanel("Description",
                                                    div(id = "desc_1", class = "card",
                                                        shiny::includeMarkdown("Rmarkdown_files/all_gp_scores_desc.Rmd")
                                                    )
                                           ),
                                           tabPanel("Results",
                                                    tagList(
                                                      h3(tags$strong("Map: Accessibility geo score")),
                                                      leafletOutput("map4"),
                                                      tags$hr(),
                                                      h3(tags$strong("Map: Connectivity geo score")),
                                                      leafletOutput("map5"),
                                                      tags$hr(),
                                                      h3(tags$strong("Map: Environmental geo score")),
                                                      leafletOutput("map6")
                                                    )
                                                    # div(id = "params_2", class = "card",
                                                    #    leafletOutput("map4")
                                                    # )
                                                    
                                           )
                                           
                                    )
                                  )
                                )
                       ),
                       tabPanel("Geo score assessment",
                         sidebarLayout(
                           sidebarPanel(
                             h3(tags$strong("Simulation:")),
                             sliderInput(inputId = "n_points",
                                         label   = "Number of simulations:",
                                         min     = 5,
                                         max     = 10,
                                         value   = 5 ),
                             sliderInput(inputId = "bf_size",
                                         label   = "Buffer size (in Km):",
                                         min     = 50,
                                         max     = 200,
                                         value   = 100 ,
                                         post    = "Km"),
                             tags$hr(),
                             tags$h5(tags$b("Progress: ")),
                             shinyWidgets::progressBar(id = "pg", value = 0, display_pct = TRUE),
                             withBusyIndicatorUI(
                               bsButton("do_validation", size="default",label = "Calculate", block = F, style="primary")
                               
                             )
                           ),
                           mainPanel(
                             tabBox(width = 12,
                                    id = "val_res",
                                    tabPanel("Description",
                                             div(id = "desc_x1", class = "card",
                                                 shiny::includeMarkdown("Rmarkdown_files/vali_descp.Rmd")
                                                 
                                             )
                                    ),
                                    tabPanel("Results",
                                                 accordion(
                                                   id="acordion2",
                                                   accordionItem(title = "Gap Map",
                                                                 collapsed = T,
                                                                 leafletOutput("map7") ),
                                                   accordionItem(title = "Summary metrics",
                                                                 fluidRow(
                                                                   plotOutput("plt1"),
                                                                   plotOutput("plt2")
                                                                   )
                                                                 ),
                                                   accordionItem(title = "Report",
                                                                 fluidRow(
                                                                   uiOutput("renderReport")
                                                                 )
                                                     
                                                   )
                                             )
                                             
                                             
                                    )
                                    
                             )
                             
                           )
                           
                         )
                         
                       )
                       
            )
            
            
    )
    
    )#end tabitems
  
)#end body


shinydashboardPlus::dashboardPage(skin = "green",
  header,
  sidebar,#dashboardSidebar(disable = F),
  body
)


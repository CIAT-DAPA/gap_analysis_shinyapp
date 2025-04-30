# All the code in this file needs to be copied to your Shiny app, and you need
# to call `withBusyIndicatorUI()` and `withBusyIndicatorServer()` in your app.
# You can also include the `appCSS` in your UI, as the example app shows.

# =============================================

# Set up a button to have an animated loading indicator and a checkmark
# for better user experience
# Need to use with the corresponding `withBusyIndicator` server function
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    htmltools::span(
      class = "btn-loading-container",
      hidden(
        img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("circle-exclamation"),
              tags$b("Error: "),
              htmltools::span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                     time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

appCSS <- "
.btn-loading-container {
  margin-left: 10px;
  font-size: 1.2em;
}
.btn-done-indicator {
  color: green;
}
.btn-err {
  margin-top: 10px;
  color: red;
}
"

short_info <- function(input, title, place){
  input %>% shinyInput_label_embed(shiny_iconlink(name = "circle-info") %>% bs_embed_tooltip( title = title, placement = place))
}


check_folder_estructure<- function(root_path){
  
  n_crops <- list.dirs(root_path, full.names = F, recursive = F)
  
  dirs_in <- data.frame(pth = list.dirs(root_path, full.names = T, recursive = T), 
             dir_names = list.dirs(root_path, full.names = F, recursive = T)%>% 
               stringr::str_extract(., "[0-9A-Za-z_-]+$"))
  
  
  dirs_must <- c("")
}


# 
# df <- data.frame(id = 1:50, resultado = c(rep(1, 19), rep(0, 31)),  inicio_sintomas = c(sample(6:15, 19, replace = T), rep(NA, 31)) , fin_sintomas = c(sample(3:9, 19, replace = T), rep(NA, 31)) , sexo = sample(c("M", "F"), 50, replace = T), edad = rnorm(50, mean = 25, sd = 6) %>% round(., 0), estrato = c(sample(1:4, 19, replace = T, prob = c(0.4, 0.3, 0.2, 0.1) ), sample(2:5, 31, replace = T, prob = c(0.1, 0.2, 0.3, 0.4))) )
# 
# dfs <- lapply(c(15, 25, 35), function(i){
#   
#   df %>% 
#     dplyr::slice_sample(., n = i)
#   
#   
# })
# 
# names(dfs) <- c("15_indvs", "25_indvs", "35_indvs")
# dfs$`50_indvs` <- df
# 
# p_load(writexl)
# writexl::write_xlsx(dfs, path ="D:/OneDrive - CGIAR/Documents/database_incidencia.xlsx")

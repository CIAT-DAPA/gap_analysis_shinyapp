


generate_report <- function(out_file,
                            paths,
                            resources){
  

  cat("Generating HTML report... \n")
  print(resources$final_gap_map)
  # paths <- list()
  # resources <- list()
  # paths$occName <- "Black_fonio"
  # paths$original_path <- "D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/input_data/occurrences_files/black_fonio_original_data_in.csv"
  # resources$spData <- read.csv("D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/black_fonio_occurrences.csv")
  # paths$sdm_occ_path <- "D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/species_distribution/black_fonio_sdm_median.tif"
  # resources$final_gap_map <- terra::rast("D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/black_fonio_final_gap_map.tif")
  # resources$overall_simulation_metrics <- readxl::read_excel("D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/gaps_validation/simulation_results.xlsx", 
  #                                                            sheet = "overall_simulation_metrics")
  # paths$sdm_outDir <- "D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/species_distribution"
  # resources$over_all_auc <- 0.6862
  # paths$occ_csv <-  "D:/OneDrive - CGIAR/Documents/prueba_1/black_fonio/africa/results/black_fonio/species_distribution/occurrences.csv"
  # 
  occ_G <- read.csv(paths$occ_csv)
  sdm_res_pth <- paste0(paths$sdm_outDir, "/sdm_results.rds")
  sp_occ <- terra::vect(paths$occ_shp)#terra::vect(resources$cleaned_data, geom = c("Longitude", "Latitude" ), crs = "epsg:4326" )
  sdm <- terra::rast(paths$sdm_occ_path)
  sdm_res <- readRDS(sdm_res_pth)
  
  auc_m <- round(mean(unlist(sdm_res$AUC)),3)*100
  auc_txt <- dplyr::case_when(auc_m <= 0.5 ~ "lower", 
                              auc_m > 0.5 & auc_m <= 0.75 ~ "regular",
                              auc_m > 0.75 ~ "good")
  
  auc_tbl <- data.frame("CV Fold" = paste0("Fold ", 1:length(sdm_res$AUC)),
                        "AUC"  = paste0(round(unlist(sdm_res$AUC),3)*100, "%"))
  auc_tbl[nrow(auc_tbl)+1, ] <-   c("Average",
                                    paste0(auc_m, "%"))
  
  #final gap map
  auc_txt_gap <- dplyr::case_when(resources$over_all_auc <= 0.5 ~ "lower", 
                                  resources$over_all_auc > 0.5 & resources$over_all_auc <= 0.75 ~ "regular",
                                  resources$over_all_auc > 0.75 ~ "good")
  
  pal <- colorBin("RdYlBu", sdm, bins = c(0, .2, .4, .6, .8, 1), na.color = "#00000000", reverse = T)
  full_df_nrow <- length(count.fields(paths$original_path))-1
  
  up_cov <- resources$overall_simulation_metrics$Value[4]
  lw_cov <- resources$overall_simulation_metrics$Value[5]
  me_cov <- (up_cov+lw_cov)/2
  
    
    
    cat(
      "---
title: \"Gap Analysis Landrace Results\"
date: \"\`r format(Sys.Date(), '%m-%d-%Y\')\`\"
output: html_document
---

\`\`\`{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
\`\`\`

## Summary

Crop landraces have been defined as “dynamic population(s) of a 
cultivated plant that has historical origin, distinct identity and 
lacks formal crop improvement, as well as often being genetically diverse, 
locally adapted and associated with traditional farming systems” . 
A landrace can be further classified as autochthonous when grown in the 
original location where it developed its unique genetic and socioeconomic 
characteristics through grower selection and allochthonous when introduced 
from another region and then locally adapted (Ramirez-Villegas, et al., 2020).

---

## Results for Landracre Group: \`r paths$occName \`

### 1. Spatial distribution:

A data cleaning process was applied, during which accessions outside the continental mask and those with duplicated coordinates were removed. 
The cleaned passport data used for the gap analysis comprises a total of ***\`r nrow(resources$spData) \`*** accessions out of a total of ***\`r full_df_nrow \`***,
representing  the ***\`r paste0(round(nrow(resources$spData)/full_df_nrow*100, 0), '%') \` *** of the records. The next graph shows the modelled spatial
distribution for the specie. Used accession are shown as darker points.

\`\`\`{r gra1, echo = FALSE, fig.cap='Species distribution map. Dark points indicate the locations of accessions.'}

leaflet(resources$spData) %>% 
    addTiles() %>% 
    addCircleMarkers(lng = ~Longitude,
               lat = ~Latitude,
               radius = 5,
               fillColor = 'black',
               stroke = F,
               fillOpacity   = 0.8) %>% 
    addRasterImage(x = sdm,
                   colors = pal,
                   opacity = 0.8) %>% 
    addLegend(pal = pal,
              values = terra::values(sdm),
              title = \"Probability\")
  

\`\`\`


\`\`\`{r tbl1, echo = F}

knitr::kable(auc_tbl, align = 'cc', caption = 'Cross validation Area under ROC curve (AUC)')

\`\`\`


Distribution model achieved an Area under ROC curve of ***\`r paste0(auc_m, \"%\")\`***  which indicate a ***\`r auc_txt \`*** level performance,
this metric indicated that: if we randomly select one accessions from the pool of presecences and pseudo-absences, the model
will have a \`r paste0(auc_m, \"%\")\` chance of correctly classifying it as either a presence or pseudo-absence. 

---

### 2. Gap Analysis Results:

The map highlights priority areas for seed collection by genebanks, illustrating regions where genetic resources are 
underrepresented or absent. These areas are identified using three gap scores, which assess the geographic 
and environmental diversity of existing ex situ conservation collections relative to the specie modeled distribution. 

The following map highlights areas in red as potential gaps in genebank collections. 
These areas represent regions where all three gap scores overlap, indicating a high likelihood of 
finding uncollected sample materials. Dark points on map represent the location of accession that are already collected.


\`\`\`{r gpmap1, echo = F, caption = 'Result map for gap analysis.'}

    leaflet(occ_G) %>% 
     addTiles() %>% 
     addCircleMarkers(lng = ~Longitude,
               lat = ~Latitude,
               radius = 3,
               fillColor = 'black',
               stroke = F,
               fillOpacity   = 0.8) %>% 
       addRasterImage(x =  resources$final_gap_map, colors = c(\"#d9d9d9\", \"#80bfff\", \"#ff471a\" ), opacity = 0.7) %>%
       addLegend(colors = c(\"#d9d9d9\", \"#80bfff\", \"#ff471a\" ),
          labels = c(\"Low\", \"Medium\", \"High\"),
          values = terra::values(resources$final_gap_map),
          title = \"Gap probability\", 
          group = \"dela\", 
          position = \"bottomleft\")


\`\`\`

A process of simulation it's implemented to obtain some perfomance metrics to support in gap maps interpretation.
Througth the AUC we determine that the estimated gap areas accounts for an averaged precision of ***\`r paste0(round(resources$over_all_auc, 3)*100, '%') \`*** which 
indicate ***\`r auc_txt_gap\`*** performance. 

---

### 3. Coverage of existing germplasm collections.

Coverage is a metric used to assess the level of diversity of a species that is preserved in a 
GenBank collection. It is defined as the proportion of the species' distribution area 
identified as a gap divided by the total area of the species' distribution.

Coverage is calculated based on the agreement of gap scores:

- The upper limit is estimated using high-probability gaps, where all three gap scores align.
- The lower limit is determined using medium-probability gaps, where at least two gap scores align.
- Coverage value is calculated as the mid point between the upper and lower limit.


\`\`\`{r gauge1, echo = F}

plot_ly(
    domain = list(x = c(0, 1), y = c(0, 1)),
    value  = me_cov,
    title  = list(text = \"Coverage (%)\"),
    type   = \"indicator\",
    mode   = \"gauge+number\",
    delta = list(reference = lw_cov),
    gauge = list(
      bar = list(color = \"black\"),
      axis =list(range = list(NULL, 100)),
      steps = list(
        list(range = c(0, 33), color = \"#ff6443\"),
        list(range = c(33, 66), color = \"#FFC300\"),
        list(range = c(66, 100), color = \"#37e300\")) ) 
      
    ) %>% 
    plotly::layout(margin = list(l=20,r=30))

\`\`\`

The ***\`r paste0(me_cov, '%')\`*** of the specie distribution is already covered by the Genebanks collections.




---

# References

- Ramirez‐Villegas, J., Khoury, C. K., Achicanoy, H. A., Mendez, A. C., Diaz, M. V., Sosa, C. C., ... & Guarino, L. (2020). A gap analysis modelling framework to prioritize collecting for ex situ conservation of crop landraces. Diversity and Distributions, 26(6), 730-742.


"  
      ,  file = out_file)
  
  rmarkdown::render(out_file,  "html_document")

  #,this value could vary from [***\`r paste0(lw_cov, '%')\`*** - ***\`r paste0(up_cov, '%')\`***].
  
  cat("Process done. Final report generated at:" , out_file,".html  \n")
  #file.remove(out_file)

}#END FUNCTION

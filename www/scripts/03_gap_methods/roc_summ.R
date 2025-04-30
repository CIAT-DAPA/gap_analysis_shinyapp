roc_summ <- function(observe, score){
  
  if(length(levels(factor(observe))) > 1){
    #val <- pROC::roc(response = factor(bd$observe), predictor = bd$score )
    croc <- suppressMessages( pROC::roc(response =observe, predictor = score))
    croc_summ <- data.frame (sensi = croc$sensitivities, speci = croc$specificities, threshold =  croc$thresholds) %>% 
      round(., 3) %>% 
      dplyr::mutate(., max.TSS = sensi + speci - 1) %>% 
      dplyr::mutate(., minROCdist = sqrt((1- sensi)^2 + (speci -1)^2))
    
    max.tss <- croc_summ %>% dplyr::filter(., max.TSS == max(max.TSS)) %>% 
      dplyr::mutate(., method = rep("max(TSS)", nrow(.)))
    
    minRoc <- croc_summ %>% 
      dplyr::filter(., minROCdist == min(minROCdist))%>% 
      dplyr::mutate(., method = rep("minROCdist", nrow(.)))
    
    croc_summ <- rbind(max.tss, minRoc) %>% 
      dplyr::filter(., speci == max(speci))  %>% 
      dplyr::sample_n(., 1) %>%
      dplyr::mutate(auc = round(croc$auc,3)) %>%
      dplyr::select(threshold, auc, sensi, speci, max.TSS)
    
    
  }else{
    croc_summ <- data.frame(threshold = NA, auc = NA, sensi = NA, speci = NA, max.TSS = NA)
    
  }
  return(croc_summ)
}

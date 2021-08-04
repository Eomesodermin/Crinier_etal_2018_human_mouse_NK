

Predict.IFN <- function(input.seurat, 
                        gene.signature = IFNsig,
                        sig.name = "IFNsig", 
                        upregulated = "TRUE",
                        deviation.multiplyer = 2,
                        quantile.threshold = 0.85){
  
  
  #browser()
  
  # Add module scores
  input.seurat <- AddModuleScore(input.seurat, 
                                 features = list(gene.signature), 
                                 name = sig.name,
                                 ctrl = 80,
                                 seed = 42)
  
  # Seurat always appends "1" to signature ID. therefore 
  seurat.sig.name <- paste0(sig.name, "1")
  
  # Calculate 75% quantile for threshold
  stats.vals <- input.seurat@meta.data %>%
    summarize(mean = mean(get(seurat.sig.name)), 
              n = n(), 
              standard_deviation = sd(get(seurat.sig.name)),
              median = median(get(seurat.sig.name)),
              qs = quantile(get(seurat.sig.name), c(quantile.threshold)))
  
  
  threshold.val <- stats.vals$median + (stats.vals$standard_deviation * deviation.multiplyer)
  
  if(upregulated){
  output.meta <- input.seurat@meta.data %>%
    dplyr::mutate(Prediction = case_when(get(seurat.sig.name) <= threshold.val ~ "Other", 
                                         get(seurat.sig.name) > threshold.val ~ "ILC1-like"))
  
  }else{
    
    output.meta <- input.seurat@meta.data %>%
      dplyr::mutate(Prediction = case_when(get(seurat.sig.name) <= threshold.val ~ "ILC1-like", 
                                           get(seurat.sig.name) > threshold.val ~ "Other"))
    
  }
  
  
  
  
  # Add Prediction to metadata of input.seurat  
  input.seurat <- AddMetaData(input.seurat, 
                              metadata = output.meta$Prediction, 
                              col.name = paste0(sig.name, "_Prediction"))
  
  
  return(input.seurat)
  
}

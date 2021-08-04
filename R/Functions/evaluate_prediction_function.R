

evaluate.prediction <- function(input.seurat, 
                                sig.name = "IFNsig",
                                validation.column, 
                                validation.cluster){
  
  
  output.temp <- input.seurat@meta.data %>%
    mutate(response_var = case_when(get(validation.column) == get("validation.cluster") ~ 1,
                                    TRUE ~ 0), 
           prediction_var = case_when(get(paste0(sig.name, "_Prediction")) == "IFN_cell" ~ 1, 
                                      TRUE ~ 0))
  
  roc(response = output.temp$response_var,
      predictor = output.temp$prediction_var, 
      ci = TRUE, 
      ci.alpha = 0.9,
      stratified = FALSE,
      plot = TRUE, 
      auc.polygon = TRUE, 
      max.auc.polygon = TRUE, 
      grid = TRUE, 
      print.auc = TRUE, 
      show.thres = TRUE)
  
  
  
  
################################
# calculate false positives 
################################
  
  false.positives <- input.seurat@meta.data %>%
    dplyr::filter(get(validation.column) != get("validation.cluster") & get(paste0(sig.name, "_Prediction")) == "IFN_cell") %>%
    dplyr::select(get("validation.column"), paste0(sig.name, "_Prediction")) %>%
    dplyr::count(get(validation.column))
  
  ################################
  # calculate false Negatives 
  ################################
  
  false.negatives <- input.seurat@meta.data %>%
    dplyr::filter(get(validation.column) == get("validation.cluster") & get(paste0(sig.name, "_Prediction")) != "IFN_cell") %>%
    dplyr::select(get("validation.column"), paste0(sig.name, "_Prediction")) %>%
    dplyr::count(get(validation.column))
  
  
  
  ################################
  # calculate total number of cells in no IFN group 
  ################################
  
  false.positives.total <- input.seurat@meta.data %>%
    dplyr::filter(get(validation.column) != get("validation.cluster")) %>%
    dplyr::select(get("validation.column"), paste0(sig.name, "_Prediction")) %>%
    dplyr::count(get(validation.column))
  
  ##################################################
  # calculate total number of cells in ifn group 
  ##################################################
  
  
  false.negatives.total <- input.seurat@meta.data %>%
    dplyr::filter(get(validation.column) == get("validation.cluster")) %>%
    dplyr::select(get("validation.column"), paste0(sig.name, "_Prediction")) %>%
    dplyr::count(get(validation.column))
  
  
  # Plot frequency 
  
  plot.df <- matrix(nrow = 1, ncol = 2)
  
  rownames(plot.df) <- c("nval")
  colnames(plot.df) <- c("False_Pos_freq", "False_Neg_freq")
  plot.df[1] <- sum(false.positives$n)/sum(false.positives.total$n) * 100
  plot.df[2] <- false.negatives$n/false.negatives.total$n * 100
  
  
  print(
    barplot(plot.df)
  )
  
  print(plot.df)
  
}
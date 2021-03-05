# Single cell RNA-sequencing functions 

# UMAP optimisation function 
UMAP.optimise <- function(input.seurat, 
                          output.dir = NULL, 
                          init.min.dist = 0.001, 
                          init.neigh = 5,
                          min.dist.step = 2,
                          neigh.step = 20){
  # Info: 
  # n.neighbor sensible values = 5 - 50 
  # min.dist sensible values = 0.001 - 0.5
  
  # Create directory if one isnt already created
  if(is.null(output.dir)){
    if(!dir.exists("output/Optimising_UMAP")){
      dir.create("output/Optimising_UMAP", 
                 recursive = TRUE)}
    
    output.dir <- "output/Optimising_UMAP/"
  }
  
  
  
  # set starting min.dist
  min.dist.val <- init.min.dist
  
  while(min.dist.val < 1){
    
    # set starting neighbor value
    n.neigh.val <- init.neigh
    
    while(n.neigh.val <= 50){
      
      print(paste0("Calculating UMAP for min.dist = ", min.dist.val, " & n.neighbor = ", n.neigh.val))
      
      # Calculate UMAP
      temp <- RunUMAP(object = input.seurat,
                      reduction = "pca",
                      dims = 1:20,
                      umap.method = "uwot",
                      n.neighbors = n.neigh.val,
                      min.dist = min.dist.val, 
                      seed.use = 42)
      
      # Plot UMAP
      print(UMAPPlot(object = temp,
                     label = TRUE, 
                     label.size = 4) + 
              ggtitle(paste0("UMAP Min dist = ", min.dist.val, " n.val = ", n.neigh.val)) + 
              NoLegend())
      
      dev.copy(pdf, paste0(output.dir, "UMAP_Min_dist_", min.dist.val, "_neighval_", n.neigh.val, ".pdf"))
      dev.off()
      
      n.neigh.val <- n.neigh.val + neigh.step
    }
    min.dist.val <- min.dist.val * min.dist.step
  }
  
  rm(temp)
  
}







# Export function
save.data.frame.function <- function(df, path, title){
  
  x <- as.data.frame(as.matrix(df))
  
  write.table(x, 
              paste0(path, title, ".txt"),
              sep = "\t",
              quote = FALSE, 
              row.names = TRUE, 
              col.names = NA)
}





#################
# Volcano Plots
#################

# edited enhancedvolcano() function to have the default variables for aesthetics I want

clean.data <- function(input.data, group.id){
  
  output.data <- input.data %>%
    dplyr::filter(cluster == paste0(group.id))
  
  rownames(output.data) <- output.data$gene
  
  output.data <- output.data %>%
    dplyr::select(avg_log2FC, p_val_adj)
  
  colnames(output.data) <- c("logFC", "FDR")
  return(output.data)
  
}



colour.points <- function(volcano.data, 
                          increase.col = "Red",
                          decreased.col = "Blue",
                          FDR.cutoff = 0.05,
                          logFC.cutoff = 0.25){
  
  # set the base colour
  keyvals <- rep('grey50', nrow(volcano.data))
  
  # set the base name/label as 'NS'
  names(keyvals) <- rep('NS', nrow(volcano.data))
  
  # modify keyvals for vars meeting FDR and LogFC threshold
  
  # Increased
  keyvals[which(volcano.data$logFC > logFC.cutoff & volcano.data$FDR < FDR.cutoff)] <- increase.col
  names(keyvals)[which(volcano.data$logFC > logFC.cutoff & volcano.data$FDR < FDR.cutoff)] <- 'Increased'
  
  # Decreased
  keyvals[which(volcano.data$logFC < -logFC.cutoff & volcano.data$FDR < FDR.cutoff)] <- decreased.col
  names(keyvals)[which(volcano.data$logFC < -logFC.cutoff & volcano.data$FDR < FDR.cutoff)] <- 'Decreased'
  
  return(keyvals)
  
}




# transparent colour function
transparent.col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  return(t.col)
}






clean.data <- function(sample.data, metadata){
  
  # Set rownames as gene IDs
  rownames(sample.data) <- sample.data$X
  
  # remove former gene ID column
  sample.data <- sample.data[,-grep("^X$", colnames(sample.data))]
  
  # Filter counts matrix for cells found in metadata
  keep.logic <- colnames(sample.data) %in% metadata$X
  sample.data <- sample.data[ , keep.logic]
  
  return(sample.data)
}

IKAP.modifed <- function (sobj, pcs = NA, pc.range = 20, k.max = NA, r.kmax.est = 1.5, 
          out.dir = "./output/IKAP", scale.data = TRUE, confounders = c("nUMI", 
                                                                 "percent.mito"), plot.decision.tree = TRUE, random.seed = 42) 
{
  dir.create(out.dir, recursive = T)


  
  cat("Running PCA ... \n")
  if (is.na(pcs)) {
    sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
    pc.change <- which(abs(diff(sobj@reductions$pca@stdev)/sobj@reductions$pca@stdev[2:length(sobj@reductions$pca@stdev)]) > 
                         0.1)
    while (length(pc.change) > 0 && max(pc.change) + pc.range + 
           2 > length(sobj@reductions$pca@stdev)) {
      sobj <- RunPCA(sobj, npcs = max(pc.change) + 
                       pc.range + 2 + 10, verbose = F)
      pc.change <- which(abs(diff(sobj@reductions$pca@stdev)/sobj@reductions$pca@stdev[2:length(sobj@reductions$pca@stdev)]) > 
                           0.1)
    }
    pcs <- if (length(pc.change) == 0){ 
      2:(pc.range + 2)
    }else{ (max(pc.change) + 2):(max(pc.change) + pc.range + 2)
  }
  }else {
    sobj <- RunPCA(sobj, npcs = max(pcs))
  }
  
  
  if (is.na(k.max)) {
    cat("Determine k.max.\n")

    sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:min(pcs))
    
    sobj <-  FindClusters(object = sobj, resolution = r.kmax.est, verbose = FALSE)
    
    k.min.pc <- length(unique(sobj@active.ident))
    
  
    sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:max(pcs))
    
    sobj <-  FindClusters(object = sobj, resolution = r.kmax.est, verbose = FALSE)
    k.max.pc <- length(unique(sobj@active.ident))
  
    
    k.max <- as.integer((k.min.pc + k.max.pc)/2)
    cat("k.max =", k.max, "\n")
  }
  gap.gain <- data.frame(matrix(NA, ncol = k.max - 1, nrow = length(pcs)))
  colnames(gap.gain) <- as.character(2:k.max)
  rownames(gap.gain) <- paste0(pcs)
  cat("Perform clustering for every nPC:\n")
  for (npc in pcs) {
    clusterings <- BottomUpMerge.mod(sobj, k.max, npc, random.seed)
    gap.stat <- GapStatistic(sobj@reductions$pca@cell.embeddings[, 
                                                         1:npc], clusterings)
    names(clusterings) <- paste0("PC", npc, "K", 1:k.max)
    sobj@meta.data <- cbind(sobj@meta.data, as.data.frame(clusterings)[, 
                                                                       2:k.max])
    gap.gain[as.character(npc), ] <- diff(gap.stat$gap)
  }
  candidates <- SelectCandidate(gap.gain)
  cat("Compute marker gene lists ... \n")
  markers.all <- ComputeMarkers.mod(sobj, gap.gain, candidates, 
                                out.dir)
  cat("Build decision tree ... \n")
  summary.rpart <- DecisionTree(sobj, markers.all, out.dir, 
                                plot.decision.tree)
  cat("Plotting summary ... \n")
  PlotSummary(gap.gain, summary.rpart, markers.all, out.dir)
  return(sobj)
}


#Warning: The following arguments are not used: reduction.type, dims.use, print.output, save.SNN
#Suggested parameter: reduction instead of reduction.type; dims instead of dims.use; verbose instead of print.output





BottomUpMerge.mod <- function (sobj, k.max, npc, random.seed) 
{
  k.clustering <- 0
  clusterings <- list()
  clust.r <- 1
  
  sobj <- FindNeighbors(sobj, 
                        reduction = "pca", 
                        dims = 1:npc)
  
  sobj <-  FindClusters(object = sobj, resolution = clust.r, verbose = FALSE)
  
  
  cat("Iteration for nPC =", npc, ", r = 1.0")
  while (length(unique(sobj@active.ident)) < k.max) {
    clust.r <- clust.r + 0.2
    cat(",", clust.r)
    sobj <- FindClusters(object = sobj, 
                         resolution = clust.r, 
                         verbose = FALSE,
                         random.seed = random.seed)
  }
  cat("\n")
  clusterings[[(k.clustering <- length(unique(sobj@active.ident)))]] <- as.character(sobj@active.ident)
  while (k.clustering > 2) {
    merged <- NearestCluster(sobj@reductions$pca@cell.embeddings[, 
                                                         1:npc], clusterings[[k.clustering]])
    clustering.merged <- clusterings[[k.clustering]]
    clustering.merged[which(clustering.merged == merged[1])] <- merged[2]
    clusterings[[k.clustering - 1]] <- as.character(as.integer(as.factor(clustering.merged)))
    k.clustering <- k.clustering - 1
  }
  clusterings[[1]] <- rep("1", nrow(sobj@reductions$pca@cell.embeddings))
  return(clusterings[1:k.max])
}




ComputeMarkers.mod <- function (sobj, gap.gain, candidates, out.dir) 
{
  markers.all <- list()
  out.xls <- list()
  out.xls$gap.gain <- cbind(PC_K = rownames(gap.gain), gap.gain)
  for (i in 1:length(candidates$k)) {
    clustering.label <- paste0("PC", candidates$pc[i], "K", 
                               candidates$k[i])
    sobj <- SetIdent(sobj, value = clustering.label)
    sobj <- RunUMAP(sobj, dims = 1:candidates$pc[i])
    ggsave(UMAPPlot(sobj), filename = paste0(out.dir, 
                                                            "/", clustering.label, "_UMAP.pdf"))
    sobj.markers <- FindAllMarkers(object = sobj, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.25)
    sobj.markers$AUROC <- NA
    for (j in 1:nrow(sobj.markers)) {
      sobj.markers$AUROC[j] <- roc.curve(scores.class0 = sobj@assays$RNA@data[sobj.markers$gene[j], 
      ], weights.class0 = sobj@active.ident == sobj.markers$cluster[j])$auc
    }
    top.10 <- sobj.markers %>% group_by(cluster) %>% top_n(10, 
                                                           avg_log2FC)
    ggsave(DoHeatmap(object = sobj, 
                     genes = top.10$gene, 
                     #slim.col.label = TRUE, 
                     #remove.key = TRUE, 
                     size = 7), 
           filename = paste0(out.dir, "/", clustering.label, 
                             "_DE_genes_LCF.png"), units = "in", width = 12, 
           height = 8)
    out.xls[[clustering.label]] <- sobj.markers[, c("gene", 
                                                    "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", 
                                                    "cluster", "AUROC")]
    markers.all[[clustering.label]] <- sobj.markers
  }
  WriteXLS(out.xls, ExcelFileName = paste0(out.dir, "/data.xls"))
  saveRDS(markers.all, file = paste0(out.dir, "/markers.all.rds"))
  return(markers.all)
}

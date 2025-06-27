library(ggplot2)
library(cowplot)
library(ggplotify)
library(dplyr)
library(tidyr)

prepare_data <- function(immgen_so,
                         res=0.1,
                         group.by="label.main",
                         cluster_id=1) {
  bridge_violin <- VlnPlot(immgen_so, 
                           pt.size = 0, 
                           features = paste0("Clustering_",res,"_", cluster_id), 
                           group.by = group.by, 
                           ncol = 4)
  score.data <- bridge_violin$data
  colnames(score.data) <- c("ModuleScore", "CellType")
  score.data$CellType <- 
    gsub(pattern = "[^a-zA-Z0-9]", replacement = "_", x = score.data$CellType)
  return(score.data)
}

analyze_pairwise <- function(score_data, test_function=pairwise.t.test, ...) {
  pairwise_results <- test_function(score_data$ModuleScore, score_data$CellType, ...)
  pairwise_table <- data.frame(pairwise_results$p.value)
  name_1 <- colnames(pairwise_table)[1]
  name_n <- rownames(pairwise_table)[nrow(pairwise_table)]
  pairwise_table <- pairwise_table %>% tibble::add_row(!!name_1:=NA, .before=1)
  rownames(pairwise_table)[1] <- name_1
  pairwise_table <- pairwise_table %>% tibble::add_column(!!name_n:=NA)
  pairwise_table <- pairwise_table[,rownames(pairwise_table)]
  
  for(r in 1:(ncol(pairwise_table)-1)){
    for(c in r:nrow(pairwise_table)){
      pairwise_table[r,c] <- pairwise_table[c,r]
    }
  }
  
  pt <- pairwise_table %>% 
    tibble::add_column(Cell=rownames(pairwise_table), .before=1) %>% 
    tidyr::gather(key="CellType", value="PValue",-1)
  
  sign_data <- pt %>% na.omit() %>%
    tibble::add_column(tr_PValue=if_else(.$PValue==0, 300, -log10(.$PValue)))
  
  return(sign_data)
}

generate_module_score_plot <- function(score_data, custom.theme=theme_classic()) {
  plot <- ggplot(score_data, aes(x=CellType, y=ModuleScore, fill=CellType)) +
    geom_violin(scale="width") +
    geom_boxplot(width=.5, fill="white") + 
    custom.theme + 
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())+ 
    coord_flip() + 
    NoLegend()
  return(plot)
}

generate_pvalue_plot <- function(sign_data, custom.theme=theme_classic()) {
  plot <- ggplot(sign_data, aes(x=CellType, y=tr_PValue, fill=CellType)) +
    geom_violin(scale="width") +
    geom_boxplot(width=.5, fill="white") +
    custom.theme + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
    ggtitle("Distribution of pairwise t.test p-values") + 
    geom_hline(yintercept = -log10(.05))+
    coord_flip() + 
    NoLegend()
  return(plot)
}

combine_plots <- 
  function(pvalue_plot, 
           module_score_plot, 
           cluster_id, 
           rel_widths=c(2.5,2)) {
    combined_plot <- 
      cowplot::plot_grid(pvalue_plot, 
                         module_score_plot + 
                           ggtitle("Distribution of module scores"), 
                         rel_widths = rel_widths, 
                         ncol=2) %>%
      ggplotify::as.ggplot() + 
      ggtitle(paste0("Cluster ", (cluster_id - 1))) + 
      theme(plot.title=element_text(hjust=.5, face="bold", color="#000080")) 
    return(combined_plot)
  }

generate_score_and_pvalue_plots <- 
  function(immgen_so, 
           res=0.1, 
           cluster_id, 
           test=pairwise.t.test,
           group.by="label.main", 
           rel_widths=c(2.5,2), 
           custom.theme=theme_classic(),...) {
    score_data <-prepare_data(immgen_so,
                              res=res, 
                              group.by=group.by, 
                              cluster_id = cluster_id)
    sign_data <- analyze_pairwise(score_data, test_function=test, ...)
    module_score_plot <- 
      generate_module_score_plot(score_data, 
                                 custom.theme=custom.theme)
    pvalue_plot <- generate_pvalue_plot(sign_data, custom.theme=custom.theme)
    combined_plot <- 
      combine_plots(pvalue_plot, 
                    module_score_plot, cluster_id,
                    rel_widths)
    return(combined_plot)
  }

generate_all_score_and_pvalue_plots <- 
  function(immgen_so, 
           res=0.1, 
           group.by="label.main", 
           rel_widths=c(2.5,2), 
           test=pairwise.t.test, 
           custom.theme=theme_classic(), ...) {
    clusters <- colnames(immgen_so@meta.data)
    # select all columns with the pattern paste0("Clustering_",res,"_")
    clusters <- clusters[grepl(paste0("Clustering_",res,"_"), x=clusters)]
    combined_plots <- vector("list", length(clusters))
    for(cluster_id in 1:length(clusters)){
      combined_plots[[cluster_id]] <- 
        generate_score_and_pvalue_plots(immgen_so, 
                                        res=res, 
                                        cluster_id = cluster_id, 
                                        group.by=group.by, 
                                        test=test,
                                        rel_widths=rel_widths, 
                                        custom.theme=custom.theme, ...)
    }
    
    return(combined_plots)
  }

# Define the function
AddClusterAnnotations <- function(seuratObject, annotationDF, clusterColumnName, annotationColumnName) {
  # Extract cluster IDs from the Seurat object's metadata
  clusterIDs <- seuratObject@meta.data[[clusterColumnName]]
  
  # Prepare a data frame with cluster IDs and cell names
  cellClustersDF <- data.frame(seurat_clusters = clusterIDs, names = rownames(seuratObject@meta.data))
  
  
  # Merge to map annotations to each cell based on cluster IDs
  cellAnnotations <- merge(cellClustersDF, 
                           annotationDF, 
                           by = "seurat_clusters",
                           all.x = TRUE)
  rownames(cellAnnotations) <- cellAnnotations$names
  # Prepare the annotations for addition: keep only the annotation column, with row names as cell names
  cellAnnotationsForAdd <- cellAnnotations[, annotationColumnName, drop = FALSE]
  
  cellAnnotationsForAdd <- cellAnnotationsForAdd[match(rownames(seuratObject@meta.data), rownames(cellAnnotationsForAdd)), ]
  # Use AddMetaData to add the annotations to the Seurat object's metadata
  seuratObject <- AddMetaData(seuratObject, metadata = cellAnnotationsForAdd, col.name = annotationColumnName)
  
  # Return the updated Seurat object
  return(seuratObject)
}


compute_module_scores_and_pairwise_plots <- function(seurat_obj, 
                                                     top_markers_list, 
                                                     cluster_name,
                                                     prefix = "CellTypeScore", 
                                                     group.by = "label.top",
                                                     test_func = pairwise.t.test) {
  # Load necessary libraries
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(cowplot)
  library(ggplotify)
  
  # Extract markers for the specified cluster
  markers <- top_markers_list[[cluster_name]] %>% pull(gene)
  
  # Compute module scores
  so <- AddModuleScore(object = seurat_obj, features = list(markers), name = prefix)
  
  # Rename module score column
  module_col <- paste0(prefix, "1")
  colnames(so@meta.data)[which(colnames(so@meta.data) == module_col)] <- "ModuleScore"
  
  
  # Prepare data for pairwise comparisons
  score_data <- so@meta.data %>%
    dplyr::select(all_of(group.by), ModuleScore) %>%
    dplyr::rename(CellType = all_of(group.by))
  
  # Clean group names
  score_data$CellType <- gsub("[^a-zA-Z0-9]", "_", score_data$CellType)
  
  
  # Perform pairwise t-tests and structure results using analyze_pairwise
  analyze_pairwise <- function(score_data, test_function = pairwise.t.test, ...) {
    pairwise_results <- test_function(score_data$ModuleScore, score_data$CellType, ...)
    pairwise_table <- data.frame(pairwise_results$p.value)
    
    name_1 <- colnames(pairwise_table)[1]
    name_n <- rownames(pairwise_table)[nrow(pairwise_table)]
    
    pairwise_table <- pairwise_table %>% tibble::add_row(!!name_1 := NA, .before = 1)
    rownames(pairwise_table)[1] <- name_1
    pairwise_table <- pairwise_table %>% tibble::add_column(!!name_n := NA)
    print(pairwise_table)
    pairwise_table <- pairwise_table[, rownames(pairwise_table)]
    
    for (r in 1:(ncol(pairwise_table)-1)) {
      for (c in r:nrow(pairwise_table)) {
        pairwise_table[r, c] <- pairwise_table[c, r]
      }
    }
    
    pt <- pairwise_table %>% 
      tibble::add_column(Cell = rownames(pairwise_table), .before = 1) %>% 
      tidyr::gather(key = "CellType", value = "PValue", -1)
    
    sign_data <- pt %>% na.omit() %>%
      tibble::add_column(tr_PValue = if_else(.$PValue == 0, 300, -log10(.$PValue)))
    
    return(sign_data)
  }
  
  pw_table_complete <- analyze_pairwise(score_data)
  
  # Generate violin plot for module score distributions
  module_score_plot <- ggplot(score_data, aes(x = CellType, y = ModuleScore, fill = CellType)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.5) +
    coord_flip() +
    theme_classic() +
    NoLegend()
  
  
  # Generate violin plot for pairwise t-test p-values
  pvalue_plot <- ggplot(pw_table_complete, aes(x = CellType, y = tr_PValue, fill = CellType)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.5) +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05)) +
    NoLegend()
  
  
  # Combine the two plots side-by-side
  combined_plot <- cowplot::plot_grid(pvalue_plot, module_score_plot, ncol = 2, rel_widths = c(2.5, 2)) %>%
    ggplotify::as.ggplot() +
    ggtitle(paste0("Module Score and Pairwise Comparisons (", cluster_name, ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#000080"))
  
  
  # Return updated Seurat object, pairwise t-test results, and combined plot
  return(list(seurat_object = so,
              pairwise_results = pw_table_complete,
              combined_plot = combined_plot,
              top_markers = top_markers_list[[cluster_name]]))
}



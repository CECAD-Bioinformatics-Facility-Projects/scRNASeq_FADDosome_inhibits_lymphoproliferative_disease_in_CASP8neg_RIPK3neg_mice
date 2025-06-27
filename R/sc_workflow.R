setClass("SeuratWorkflow")

setClass(
  "ConcreteSeuratWorkflow",
  contains = "SeuratWorkflow"
)

setGeneric(
  name = "run_analysis",
  def = function(object, directory, id=1, dims=1:10, resolution=0.5) {
    standardGeneric("run_analysis")
  }
)


setMethod(
  f = "run_analysis",
  signature = "ConcreteSeuratWorkflow",
  definition = function(object, directory, id=1, dims=1:10, resolution=0.5) {
    
    # check if directory exists
    if(fs::dir_exists(directory)){
      reading_target <- tar_target_raw(
        paste0("tenx_Data","_",id),
        substitute(
          expr = Seurat::Read10X(directory),
          env = list(directory = directory)
        )
      )
    } else {
      if(is.data.frame(directory)){
        reading_target <- tar_target_raw(
          paste0("tenx_Data","_",id),
          substitute(
            expr = directory,
            env = list(directory = directory)
          )
        )
      } else{
        message("Directory does not exist")
        stop()
      }
    }
    
    list(
      reading_target,
      tar_target_raw(
        paste0("SO","_",id),        
        substitute(
          expr = Seurat::CreateSeuratObject(
            counts = tenx_Data,
            project = "SeuratAnalysis"
          ),
          env = list(
            tenx_Data = as.name(paste0("tenx_Data","_",id))
          )
        )
      ),
      tar_target_raw(
        paste0("SO_Normalized_Data","_",id),
        substitute(
          expr = Seurat::NormalizeData(SO),
          env = list(
            SO = as.name(paste0("SO","_",id))
          )
        )
      ),
      tar_target_raw(
        paste0("SO_Variable_Features","_",id),
        substitute(
          expr = Seurat::FindVariableFeatures(SO_Normalized_Data),
          env = list(
            SO_Normalized_Data = as.name(paste0("SO_Normalized_Data","_",id))
          )
        )
      ),
      tar_target_raw(
        # name with id
        paste0("SO_Scaled_Data","_",id),
        substitute(
          expr = Seurat::ScaleData(SO_Variable_Features),
          env = list(
            SO_Variable_Features = as.name(paste0("SO_Variable_Features","_",id))
          )
        )
      ),
      tar_target_raw(
        # name with id
        paste0("SO_PCA","_",id),
        substitute(
          expr = Seurat::RunPCA(SO_Scaled_Data),
          env = list(
            SO_Scaled_Data = as.name(paste0("SO_Scaled_Data","_",id))
          )
        )
      ),
      tar_target_raw(
        # Name with id
        paste0("SO_Neighbors","_",id),
        substitute(
          expr = Seurat::FindNeighbors(SO_PCA, dims = dims),
          env = list(
            SO_PCA = as.name(paste0("SO_PCA","_",id)),
            dims = dims
          )
        )
      ),
      tar_target_raw(
        paste0("SO_Clusters_",id),
        substitute(
          expr = Seurat::FindClusters(SO_Neighbors, resolution = resolution),
          env = list(
            SO_Neighbors = as.name(paste0("SO_Neighbors","_",id)),
            resolution = resolution
          )
        )
      ),
      tar_target_raw(
        paste0("SO_UMAP","_",id),
        substitute(
          expr = Seurat::RunUMAP(SO_Clusters, dims = dims),
          env = list(
            SO_Clusters = as.name(paste0("SO_Clusters_",id)),
            dims = dims
          )
        )
      )
    )
  }
)


find_all_markers <- function(seurat_obj, meta_column, ...){
  groups <- seurat_obj@meta.data %>% pull(meta_column) %>% unique() %>% as.character()
  print(groups)
  markers <- vector("list", length(groups))
  for(g in 1:length(groups)){
    print(paste0("Group ",groups[g]))
    markers[[g]] <- FindMarkers(object = seurat_obj, group.by=meta_column, ident.1=groups[g], ...)
    write.table(markers[[g]], file = paste0("markers_",groups[g],".tsv"), sep="\t", quote = FALSE)
    print("before")
    print(head(markers[[g]]))
    print("after")
    if(is.null(markers[[g]])){
      next
    }
    if(nrow(markers[[g]]) == 0){
      next
    }
    markers[[g]]$gene <- rownames(markers[[g]])
    markers[[g]]$population <- groups[g] 
  }
  return(do.call("rbind", markers))
}


find_all_markers_safe <- function(seurat_obj, meta_column, ...) {
  # 1. Extract cluster IDs as a sorted character vector
  clusters <- unique(seurat_obj@meta.data[[meta_column]])
  clusters <- sort(as.character(clusters))  # ensure character & sorted
  
  # 2. Optional: define a safe version of FindMarkers
  safe_find_markers <- function(object, cluster_id, group_var, ...) {
    # a. tryCatch approach:
    tryCatch(
      {
        res <- FindMarkers(
          object   = object,
          group.by = group_var,
          ident.1  = cluster_id,
          ...
        )
        # Return a data frame even if empty
        if (!is.null(res) && nrow(res) > 0) {
          res$gene       <- rownames(res)
          res$population <- cluster_id
        }
        return(res)
      },
      error = function(e) {
        message("Error in FindMarkers for cluster ", cluster_id, ": ", e$message)
        return(NULL) # return NULL on error
      }
    )
    
    # Alternatively, you can use purrr::possibly:
    # Possibly(FindMarkers, otherwise = NULL)(object, group.by = group_var, ident.1 = cluster_id, ...)
  }
  
  # 3. Iterate over each cluster ID, returning a list of data frames
  markers_list <- lapply(clusters, function(cid) {
    message("Computing markers for cluster: ", cid)
    safe_find_markers(object = seurat_obj, cluster_id = cid, group_var = meta_column, ...)
  })
  
  # 4. Filter out NULL or zero-row results, then row-bind
  markers_list <- markers_list[!vapply(markers_list, is.null, logical(1))]  # remove NULLs
  markers_list <- markers_list[vapply(markers_list, nrow, numeric(1)) > 0]  # remove 0-row
  
  # 5. Combine everything into one data frame
  if (length(markers_list) == 0) {
    message("No markers found for any cluster.")
    return(data.frame())
  }
  
  markers_df <- do.call(rbind, markers_list)
  return(markers_df)
}

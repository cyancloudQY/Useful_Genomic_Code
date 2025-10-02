## Function Helps Read single cell RNA-seq data into R


return_seurat_object <- function(files_name, files_path){
  sample_object_list <- c()
  for (i in 1:length(files_name)){
    print(files_name[i])
    print(paste0('Doing', i,' sample'))
    ### CHANGE THIS PART 
    cur_10X <- ReadMtx(
      mtx = paste0(files_path, files_name[i],"_matrix.mtx.gz"),
      cells =  paste0(files_path, files_name[i],"_barcodes.tsv.gz"),
      features =  paste0(files_path, files_name[i],"_features.tsv.gz")
    )
    ###
    cur_object <- CreateSeuratObject(cur_10X)
    cur_object@meta.data$GEO_ident <- files_name[i]
    print(paste0('reading:', files_name[i]))
    
    ############### mark doublet #############
    print('Doing Dedoublet...')
    dedoub_object <- NormalizeData(cur_object)
    dedoub_object <- FindVariableFeatures(dedoub_object, selection.method = "vst", nfeatures = 2000)
    dedoub_object <- ScaleData(dedoub_object)
    dedoub_object <- RunPCA(dedoub_object, features = VariableFeatures(object = dedoub_object))
    doublet_rate <- ncol(dedoub_object) * 8 * 1e-6
    nExp_poi <- doublet_rate * ncol(dedoub_object)
    dedoub_object <- doubletFinder(dedoub_object, PCs = 1:15, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
    
    print(dedoub_object@meta.data)
    
    colnames(dedoub_object@meta.data)[3] <- 'DF_value'
    colnames(dedoub_object@meta.data)[4] <- 'DF_result'
    
    
    cur_object@meta.data$DF_value <- dedoub_object@meta.data$DF_value
    cur_object@meta.data$DF_result <- dedoub_object@meta.data$DF_result
    ##########################################
    
    
    
    assign(files_name[i], cur_object)
    sample_object_list <- c(sample_object_list, get(files_name[i]))
    print(sample_object_list[i])
    
  }
  
  merged_object <- merge(x = sample_object_list[[1]], 
                         y = sample_object_list[-1])
  
  merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])
  return(merged_object)
}

add_meta <- function(seurat_object, meta_add, left_column, right_column){
  meta <- seurat_object@meta.data
  meta <- rownames_to_column(meta)
  meta <- left_join(meta,
                    meta_add, 
                    by = join_by(left_column == right_column))
  meta <- column_to_rownames(meta, 'rowname')
  seurat_object@meta.data <- meta
  return(seurat_object)
  
}


quality_control <- function(seurat_object, target_file){
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
  ggsave(paste0(target_file, 'qc_before.pdf'))
  
  before <- nrow(seurat_object@meta.data)
  print(before)
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 300 & 
                             percent.mt < 20
                            )
  
  VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
  ggsave(paste0(target_file, 'qc_after.pdf'))
  after <- nrow(seurat_object@meta.data)
  print(after)
  
  
  return(seurat_object)
  
}


dimension_reduction <- function(seurat_object, batch){
  seurat_object <- NormalizeData(seurat_object)
  
  ### 
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  ElbowPlot(seurat_object, ndims = 50)
  
  
  harmony_object <- RunHarmony(seurat_object, 
                               group.by.vars = batch,
                               max_iter = 50
  )
  
  harmony_object<- FindNeighbors(harmony_object, dims = 1:50, reduction = 'harmony')
  harmony_object <- FindClusters(harmony_object, resolution = 0.3)
  harmony_object<- RunUMAP(harmony_object, dims = 1:50, reduction = 'harmony')
  
  
  
  return(harmony_object)
}


consider_batch <- function(seurat_object, batch){
  harmony_object <- RunHarmony(seurat_object, 
                               group.by.vars = batch,
                               max_iter = 50
  )
  
  harmony_object<- FindNeighbors(harmony_object, dims = 1:50, reduction = 'harmony')
  harmony_object <- FindClusters(harmony_object, resolution = 0.3)
  harmony_object<- RunUMAP(harmony_object, dims = 1:50, reduction = 'harmony')
  
  return(harmony_object)
  
}

redo_reduce_dimension <- function(seurat_object, pc,res){
  seurat_object <- FindNeighbors(seurat_object, dims = 1:pc)
  seurat_object <- FindClusters(seurat_object, resolution = res)
  seurat_object <- RunUMAP(seurat_object, dims = 1:pc)
  
  return(seurat_object)
}

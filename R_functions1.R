# library(Seurat)
# library(dplyr)
# 
# 
# 
# 
# 


# function to reorder metadata according to counts column order
reorder_metadata <- function(seurat_obj){
    # get ordering of counts columns
    ordering <- data.frame(
        reordering_id=1:dim(seurat_obj@meta.data)[1],
        cell_names=colnames(LayerData(seurat_obj, layer = 'counts')))
    # reorder with a temp metadta
    temp_meta <- seurat_obj@meta.data
    temp_meta$cell_names <- rownames(temp_meta)

    # don't want duplicate columns, so delede reordering_id column if there is one
    if("reordering_id" %in% colnames(temp_meta)) {
        temp_meta <- subset(temp_meta, select=-c(reordering_id))
        }
    
    # merge the ordering in
    #print(temp_meta[1:3,])
    temp_meta <- merge(temp_meta, ordering, by.x = 'cell_names', by.y = 'cell_names')
    # sort by the ordering
    #print(temp_meta[1:3,])
    temp_meta <- temp_meta[order(temp_meta$reordering_id),]
    # reset the indexes to the cell_names, then set the metadata accordingly
    rownames(temp_meta) <- temp_meta$cell_names

    # remove reordering_id column
    temp_meta <- subset(temp_meta, select=-c(reordering_id))
    
    seurat_obj@meta.data <- temp_meta
    print('reordered!')
    return(seurat_obj)
    }

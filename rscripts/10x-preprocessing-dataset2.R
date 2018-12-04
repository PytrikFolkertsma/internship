
.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/')

library(Seurat)
library(magrittr)
library(dplyr)

n.pcs <- 21

resolutions <- c(0.5, 0.7, 1, 1.5)
#resolutions <- c(0.3, 0.5, 0.7, 1, 1.2, 1.5, 2)
#resolutions <- c(0.5)

################################################################################
# 
# print('LOADING DATA')
# 
# df.10x <- Read10X("/data/sc-10x/data-runs/171120-scheele-adipose/agg-180831-unnormalized/outs/filtered_gene_bc_matrices_mex/hg19")
# seurobj <- CreateSeuratObject(df.10x, min.cells = 3, min.genes = 200, is.expr = 0)
# 
# ################################################################################
# 
# #print('ADDING METADATA')
# 
# timepoints <- c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
# 
# #extract the indices from the cellnames
# sample_agg_idx <- as.numeric(sapply(strsplit(seurobj@cell.names, split = "-"), '[[', 2))
# 
# #assign the correct sample names
# df.metadata <- data.frame(row.names=seurobj@cell.names,
#                           timepoint=timepoints[sample_agg_idx]
#                           )
# 
# seurobj <- AddMetaData(seurobj, df.metadata)
# 
# print('nr of cells per sample')
# print(seurobj@meta.data %>% count(timepoint))
# 
# print('filter out S6')
# seurobj <- SubsetData(seurobj, cells.use=rownames(seurobj@meta.data)[seurobj@meta.data$timepoint != 'T6'])
# 
# ###ADD TIME-COMBINED METADATA
# seurobj@meta.data['time-combined'] <- unlist(lapply(seurobj@meta.data$timepoint, function(x){
#   if (x %in% c('T1', 'T2', 'T3')){
#     return('T1T2T3')
#   } else {
#     return('T4T5')
#   }
# }))
# 
# ################################################################################
# 
# print('QUALITY CONTROL')
# 
# print('calculating percent.ribo and percent.mito')
# ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(seurobj@data), value = TRUE)
# mito.genes <- grep(pattern = "^MT-", x = rownames(seurobj@data), value = TRUE, ignore.case=TRUE)
# percent.ribo <- Matrix::colSums(seurobj@raw.data[ribo.genes, ])/Matrix::colSums(seurobj@raw.data)
# percent.mito <- Matrix::colSums(seurobj@raw.data[mito.genes, ])/Matrix::colSums(seurobj@raw.data)
# seurobj <- AddMetaData(seurobj, metadata=percent.ribo, col.name="percent.ribo")
# seurobj <- AddMetaData(seurobj, metadata=percent.mito, col.name="percent.mito")
# 
# print('filter cells')
# #seurobj <- FilterCells(seurobj, subset.names=c("percent.mito"), low.thresholds = c(-Inf), high.thresholds = c(0.1))
# 
# print('normalizing data')
# seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", scale.factor = 1e4)
# 
# print('scaling data')
# #TODO: maybe regress out nUMI per timepoint
# #seurobj <- ScaleData(seurobj, vars.to.regress = c("nUMI", "percent.mito"), do.par=T, num.cores=10)
# seurobj <- ScaleData(seurobj, do.par=T, num.cores=10)
# 
# ################################################################################

seurobj <- readRDS('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/10x-180831-noreg')

print('CLUSTERING')

print('finding variable genes')
seurobj <- FindVariableGenes(seurobj, do.plot=F)

print('nr of variable genes:')
print(length(seurobj@var.genes))

print('running pca')
seurobj <- RunPCA(seurobj, pcs.compute=50, do.print=F)

print('running clustering')
for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", dims.use = 1:n.pcs, resolution = res, print.output = 0, save.SNN = TRUE, force.recalc=T)
}

################################################################################


print('CALCULATING CC SCORES')

cc.genes <- readLines('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/downloads/regev_lab_cell_cycle_genes.txt')
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurobj <- CellCycleScoring(seurobj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

################################################################################

print('RUNNING TSNE')

seurobj <- RunTSNE(seurobj, reduction.use='pca', dims.use=1:n.pcs)

################################################################################

print('SAVING DATASET')

print('saving dataset as 10x-180831-noreg')
saveRDS(seurobj, 'data/10x-180831-noreg')

###ALSO ADD: Split data into T1T2T3 - T4T5? 
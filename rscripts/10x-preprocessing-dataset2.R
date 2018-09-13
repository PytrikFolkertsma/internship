
.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/')

library(Seurat)
library(magrittr)
library(dplyr)

n.pcs <- 15

resolutions <- c(0.4, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# 
# ################################################################################
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
# #samples_info <- read.table('data/downloads/180406-Cell ID and 10x sample Index_final-extracted columns.txt', sep='\t', header=T, stringsAsFactors=F)
# #samples_info.ordered <- rbind(samples_info[13:14,], samples_info[1:12,])
# 
# #get the metadata
# #sample_names <- samples_info.ordered$Sample.name..Single.Cell.ID.
# #diff <- samples_info.ordered$Diff
# #ucp1.ctrl <- samples_info.ordered$UCP1.ctrl
# #ucp1.ne <- samples_info.ordered$UCP1.NE
# #bmi <- samples_info.ordered$BMI
# #age <- samples_info.ordered$AGE
# timepoints <- c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
# 
# #get the sample names consistent
# #sample_names[1] <- 'Supra_4'
# #sample_names[2] <- 'Subq_4'
# 
# #get second labeling (only supra, subq, peri and visce)
# #sample_names2 <- unlist(lapply(sample_names, function(x){
# #  return(substring(x, 0, nchar(x)-2))
# #}))
# 
# #extract the indices from the cellnames
# sample_agg_idx <- as.numeric(sapply(strsplit(seurobj@cell.names, split = "-"), '[[', 2))
# 
# #assign the correct sample names
# df.metadata <- data.frame(row.names=seurobj@cell.names,
#                           timepoint=timepoints[sample_agg_idx]
#                           # sample_name2=sample_names2[sample_agg_idx],
#                           # sample_name=sample_names[sample_agg_idx],
#                           # diff=diff[sample_agg_idx],
#                           # ucp1.ctrl=ucp1.ctrl[sample_agg_idx],
#                           # ucp1.ne=ucp1.ne[sample_agg_idx],
#                           # bmi=bmi[sample_agg_idx],
#                           # age=age[sample_agg_idx]
#                           )
# 
# seurobj <- AddMetaData(seurobj, df.metadata)
# 
# print('nr of cells per sample')
# print(seurobj@meta.data %>% count(timepoint))
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
# #VlnPlot(seurobj, c("nGene"), group.by='timepoint', nCol = 1, point.size.use=-1)
# #VlnPlot(seurobj, c("percent.mito"), group.by='timepoint', nCol = 1, point.size.use=-1)
# #VlnPlot(seurobj, c("nUMI"), group.by='timepoint', nCol = 1, point.size.use=-1)
# 
# print('filter out S6')
# seurobj <- SubsetData(seurobj, cells.use=rownames(seurobj@meta.data)[which(!(seurobj@meta.data$timepoint %in% 'T6'))])
# 
# print('filter cells')
# seurobj <- FilterCells(seurobj, subset.names=c("percent.mito"), low.thresholds = c(-Inf), high.thresholds = c(0.1))
# 
# print('normalizing data')
# seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", scale.factor = 1e4)
# 
# print('scaling data (regressing out nUMI and percent.mito)')
# seurobj <- ScaleData(seurobj, vars.to.regress = c("nUMI", "percent.mito"), do.par=T, num.cores=10)
# 
# ################################################################################
# 
# print('CLUSTERING')
# 
# print('finding variable genes')
# seurobj <- FindVariableGenes(seurobj, do.plot=F)
# 
# print('nr of variable genes:')
# print(length(seurobj@var.genes))
# 
# print('running pca')
# seurobj <- RunPCA(seurobj, pcs.compute=50, do.print=F)
# 
# print('saving dataset as 10x-180831')
# saveRDS(seurobj, 'data/10x-180831')

###

seurobj <- readRDS('data/10x-180831')

print('running clustering')
for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", dims.use = 1:n.pcs, resolution = res, print.output = 0, save.SNN = TRUE)
}

################################################################################


print('CALCULATING CC SCORES')

cc.genes <- readLines('data/downloads/regev_lab_cell_cycle_genes.txt')
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurobj <- CellCycleScoring(seurobj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

################################################################################

print('RUNNING TSNE')

seurobj <- RunTSNE(seurobj, reduction.use='pca', dims.use=1:n.pcs)

################################################################################

print('SAVING DATASET')

print('saving dataset as 10x-180831')
saveRDS(seurobj, 'data/10x-180831')


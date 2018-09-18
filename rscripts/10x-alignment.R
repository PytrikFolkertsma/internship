
### Runs CCA, aligns subspaces 30 CC's, runs tSNE on 15.
### Discarded cells are saved in data/alignment_discarded-cells

### 180504 should be aligned on ... cc's
### 180831 should be aligned on 12 cc's 
### 180831_T1_old should be aligned on ... cc's

n.ccs <- 30
n.ccs.use <- 12

.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/scripts-10x-analysis/rscripts/')

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run the alignment on.'),
  make_option(c('-c', '--column'), type='character', help='Column to align on (e.g. sample_name or timepoint).')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$file) || is.null(opt$column)){
  print_help(optionparser)
  quit()
}

################################################################################

library(Rcpp)
library(magrittr)
library(Seurat)

################################################################################

print('LOADING DATA')

seurobj_all10x <- readRDS(opt$file)

################################################################################

print('SPLITTING DATASET PER SAMPLE')

samples <- unlist(seurobj_all10x@meta.data[opt$column])
sample_names <- unique(samples)
print('sample names:')
print(sample_names)

objectlist <- list()
variable_genes <- c()

for (i in c(1:length(sample_names))){
  print(paste(i, '- splitting data and finding variable genes for', sample_names[i]))
  objectlist[[i]] <- SubsetData(seurobj_all10x, cells.use=seurobj_all10x@cell.names[which(samples %in% sample_names[i])])
  objectlist[[i]] <- FindVariableGenes(objectlist[[i]], do.plot = FALSE)
  variable_genes <- c(variable_genes, rownames(head(objectlist[[i]]@hvg.info, n=1000)))
}

unique.variable_genes <- unique(variable_genes)
print(paste('length variable genes:', length(unique.variable_genes)))

################################################################################

print('RUNNING MULTICCA')

data <- RunMultiCCA(object.list = objectlist, genes.use=unique.variable_genes, num.ccs=n.ccs)

################################################################################

print('ALIGNING SUBSPACES')

#Discard cells whose expression profile cannot be well-explained by low-dimensional CCA, 
#compared to low-dimensional PCA. Here cells are discarded when the ratio of the variance 
#explained by CCA is smaller than 0.5 (compared to the variance explained by PCA).

print('Discarding cells whose expression profile cannot be well explained by low-dim CCA compared to low-dim PCA.')
data <- CalcVarExpRatio(data, reduction.type = "pca", grouping.var = opt$column, dims.use = 1:10)
data.all.save <- data
data <- SubsetData(object = data, subset.name = "var.ratio.pca", accept.low = 0.5)
data.discard <- SubsetData(object = data.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

print('Saving object discarded cells')
saveRDS(data.discard, paste(opt$file, 'cca-discardedcells', sep='-'))

print('Aligning subspaces')
data.aligned <- AlignSubspace(data, reduction.type = "cca", grouping.var = opt$column, dims.align = 1:n.ccs)

################################################################################

print('RUNNING TSNE')
data.aligned <- RunTSNE(data.aligned, reduction.use='cca.aligned', dims.use=1:n.ccs.use)

print('SAVING SEURAT OBJECT')
saveRDS(data.aligned, paste(opt$file, 'aligned', sep='-'))


### Runs CCA, aligns subspaces 30 CC's, runs tSNE.
### Discarded cells are saved in data/alignment_discarded-cells

.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/scripts-10x-analysis/rscripts/')

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run the alignment on.')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$file)){
  print_help(optionparser)
  quit()
}

################################################################################

library(Rcpp)
library(magrittr)
library(Seurat)

n.ccs <- 30
#ccs.tsne <- c(10, 20, 30) #nr of cc's to use for tsne

################################################################################

print('LOADING DATA')

seurobj_all10x <- readRDS(opt$file)

################################################################################

print('SPLITTING DATASET PER SAMPLE')

samples <- seurobj_all10x@meta.data$timepoint
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
data <- CalcVarExpRatio(data, reduction.type = "pca", grouping.var = "sample_name", dims.use = 1:10)
data.all.save <- data
data <- SubsetData(object = data, subset.name = "var.ratio.pca", accept.low = 0.5)
data.discard <- SubsetData(object = data.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

print('Saving object discarded cells')
saveRDS(data.discard, paste(opt$file, 'cca-discardedcells', sep='-'))

print('Aligning subspaces')
data.aligned <- AlignSubspace(data, reduction.type = "cca", grouping.var = "sample_name", dims.align = 1:n.ccs)

################################################################################

print('RUNNING TSNE')
data.aligned <- RunTSNE(data.aligned, reduction.use='cca.aligned', dims.use=1:15)

print('SAVING SEURAT OBJECT')
saveRDS(data.aligned, paste(opt$file, 'aligned', sep='-'))

# for (i in ccs.tsne){
#   data.aligned <- RunTSNE(data.aligned, reduction.use = "cca.aligned", dims.use = 1:i)
#   print(paste('saving object: ', 'data/10x__ccregout_cca30-aligned_tsne-1-', i, '.rds', sep=''))
#   saveRDS(data.aligned, paste('data/10x_ccregout_cca30-aligned_tsne-1-', i, '.rds', sep=''))
# }


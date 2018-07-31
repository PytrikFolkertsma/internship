.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/rscripts/')

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run cell cycle regression on.')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$file)){
  print_help(optionparser)
  quit()
}

library(Seurat)

print('loading Seurat object')
data <- readRDS(opt$file)

print('scaling data')
data <- ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score", "nUMI", "percent.mito"))

print('running tsne')
data <- RunPCA(data, pcs.compute=12) #nr pc's determined from elbow plot
data <- RunTSNE(data, reduction.use='pca', dims.use=1:12)

print('save Seurat object')
saveRDS(data, paste(opt$file, 'ccregout', sep='-'))


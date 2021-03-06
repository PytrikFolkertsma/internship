################################################################################
#
# Regresses out the effects of cell cycle of a given Seurat object. 
# object. Saves the regressed out data as filename-ccregout'.
#
################################################################################

resolutions <- c(0.5, 0.7, 1, 1.5)

library(optparse)
library(Seurat)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run cell cycle regression on.'),
  make_option(c('-n', '--npcs'), type='character', help='Number of PCs to use for clustering and tSNE after cell cycle regression. Default is 12.', default=12)
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$file)){
  print_help(optionparser)
  quit()
}

print('Loading Seurat object...')
data <- readRDS(opt$file)

print('Regressing out cell cycle effects...')
data <- ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score", "nUMI", "percent.mito"))

print('Running PCA...')
data <- RunPCA(data, pcs.compute=n.pcs)

print('Running clustering...')
for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", dims.use = 1:opt$npcs, resolution = res, print.output = 0, save.SNN = TRUE)
}

print('Running tSNE...')
data <- RunTSNE(data, reduction.use='pca', dims.use=1:opt$npcs)

print(paste('Saving dataset: ../output/', opt$file, '-ccregout', sep=''))
saveRDS(data, paste('../output/', opt$file, '-ccregout', sep=''))


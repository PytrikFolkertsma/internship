
################################################################################
#
# Creates a Seurat object of the 10x 180504 samples. Executes the following steps:
# - Adds all metadata
# - Filters, normalizes and scales data (filters genes 200-9000, 
#   percent.mito < 0.08, nUMI < 110000, regresses out nUMI and percent.mito)
# - Computes clusters (using 15 PC's based on elbow plot) of several resolutions 
#   and adds them to metadata
# - Calculates cell cycle scores and adds them to metadata
# - Performs tSNE (15 PC's (determined before from elbow plot), default perplexity)
#
################################################################################

.libPaths('/home/cbmr/pytrik/libraries/')

library(Seurat)
library(magrittr)
library(dplyr)
library(optparse)

n.pcs <- 15
resolutions <- c(0.5, 0.7, 1, 1.5)

option_list <- list(
  make_option(c('-o', '--output_folder'), type='character', help='Path to the output folder.')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if(is.null(opt$output_folder)){
  print_help(optionparser)
  quit()
}

################################################################################

print('LOADING DATA')

df.10x <- Read10X("/data/sc-10x/data-runs/171120-scheele-adipose/agg-180504-unnormalized/outs/filtered_gene_bc_matrices_mex/hg19")
seurobj <- CreateSeuratObject(df.10x, min.cells = 3, min.genes = 200, is.expr = 0)

################################################################################

print('ADDING METADATA')

samples_info <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/downloads/180406-Cell ID and 10x sample Index_final-extracted columns.txt', sep='\t', header=T, stringsAsFactors=F)
samples_info.ordered <- rbind(samples_info[13:14,], samples_info[1:12,])

#get the metadata
sample_names <- samples_info.ordered$Sample.name..Single.Cell.ID.
diff <- samples_info.ordered$Diff
ucp1.ctrl <- samples_info.ordered$UCP1.ctrl
ucp1.ne <- samples_info.ordered$UCP1.NE
bmi <- samples_info.ordered$BMI
age <- samples_info.ordered$AGE

#get the sample names consistent
sample_names[1] <- 'Supra_4'
sample_names[2] <- 'Subq_4'

#get depot labeling (supra, subq, peri or visce)
sample_names2 <- unlist(lapply(sample_names, function(x){
  return(substring(x, 0, nchar(x)-2))
}))

#extract the indices from the cellnames
sample_agg_idx <- as.numeric(sapply(strsplit(seurobj@cell.names, split = "-"), '[[', 2))

#create metadata df
df.metadata <- data.frame(row.names=seurobj@cell.names,
                          sample_name2=sample_names2[sample_agg_idx],
                          sample_name=sample_names[sample_agg_idx],
                          diff=diff[sample_agg_idx],
                          ucp1.ctrl=ucp1.ctrl[sample_agg_idx],
                          ucp1.ne=ucp1.ne[sample_agg_idx],
                          bmi=bmi[sample_agg_idx],
                          age=age[sample_agg_idx])

seurobj <- AddMetaData(seurobj, df.metadata)

print('nr of cells per sample')
print(seurobj@meta.data %>% count(sample_name))

################################################################################

print('QUALITY CONTROL')

print('calculating percent.mito')
mito.genes <- grep(pattern = "^MT-", x = rownames(seurobj@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(seurobj@raw.data[mito.genes, ])/Matrix::colSums(seurobj@raw.data)
seurobj <- AddMetaData(seurobj, metadata=percent.mito, col.name="percent.mito")

print('filter cells')
seurobj <- FilterCells(seurobj, subset.names=c("nGene", "percent.mito", "nUMI"), low.thresholds = c(200, -Inf, -Inf), high.thresholds = c(9000, 0.08, 110000))

print('normalizing data')
seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", scale.factor = 1e4)

print('scaling data (regressing out nUMI and percent.mito)')
seurobj <- ScaleData(seurobj, vars.to.regress = c("nUMI", "percent.mito"), do.par=T, num.cores=10)

################################################################################

print('CLUSTERING')

print('finding variable genes')
seurobj <- FindVariableGenes(seurobj, do.plot=F)
print('running pca')
seurobj <- RunPCA(seurobj, pcs.compute=50, do.print=F)

print('running clustering')
for (res in resolutions){
  seurobj <- FindClusters(seurobj, reduction.type = "pca", dims.use = 1:n.pcs, resolution = res, print.output = 0, save.SNN = TRUE)
}

################################################################################

print('CALCULATING CC SCORES')

cc.genes <- readLines('/projects/pytrik/sc_adipose/data/downloads/regev_lab_cell_cycle_genes.txt')
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
seurobj <- CellCycleScoring(seurobj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

################################################################################

print('RUNNING TSNE')

seurobj <- RunTSNE(seurobj, reduction.use='pca', dims.use=1:n.pcs)

################################################################################

print('SAVING DATASET')

print('saving dataset as 10x-180504')
saveRDS(seurobj, paste(opt$output_folder, '/10x-180504'))


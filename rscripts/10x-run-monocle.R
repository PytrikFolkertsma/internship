.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/scripts-10x-analysis/rscripts/')

r.seed <- 66

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run Monocle on.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory.'),
  make_option(c('-d', '--depots'), action='store_true', default=F, help='OPTIONAL. If flag -d is given, Monocle is run on individual depots.'),
  make_option(c('-r', '--regressout'), type='character', help="OPTIONAL. Choose from 'pm-umi', 'pm-umi-cc' or 'cc'.")
  #add option for downsampling
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$file) || is.null(opt$outdir)){
  print_help(optionparser)
  quit()  
}

if (!is.null(opt$regressout)){
  if (!(opt$regressout == 'pm-umi' || opt$regressout == 'cc' || opt$regressout == 'pm-umi-cc')){
    print("For regressout (-r) choose from 'pm-umi', 'pm-umi-cc' or 'cc'.")
    print_help(optionparser)
    quit()
  }
}

###################################################################################

library(monocle)
library(Seurat)

print('LOADING SEURAT OBJECT')
seurat_object <- readRDS(opt$file)
seurat_object <- SubsetData(SetAllIdent(seurat_object, id='sample_name'), max.cells.per.ident=1000, random.seed=r.seed)

setwd(opt$outdir)
output_prefix <- '10x'

run_monocle_workflow <- function(data, output_name){
  
  print('Converting Seurat object to CDS...')
  cds <- importCDS(data)

  print('Estimating size factors and dispersions...')
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  ###TRAJECTORY STEP 1: Choose genes that define a cell's process
  # There are a few ways to do this:
  # - For time series data: Do a differential gene test between cells collected at
  #   the beginning of the experiment and those collected at the end.
  # - For non time series data: 
  #     - dpFeature: cluster cells and do a differential gene test between clusters 
  #       (not doing this here because cells cluster mostly to their own sample, and
  #       also because the DEtest takes an extremely long time). Maybe try out: 
  #       use subtissue as 'cluster'. 
  #     - Select genes with high dispersion across cells (used below). However, what 
  #       I think might happen now is that those are mostly cell cycle genes. We can 
  #       regress out the cell cycle effects when building the trajectory, but it will
  #       still use the same genes for ordering. 
  #     - Use known marker genes.

  print('Selecting genes for ordering with high dispersion...')
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  print(paste('Nr of genes selected:', length(ordering_genes)))

  print('Seting ordering filter...')
  cds <- setOrderingFilter(cds, ordering_genes)

  ###TRAJECTORY STEP 2: Reduce data dimensionality

  print('Reducing data to two dimensions with DDRTree...')

  if (is.null(opt$regressout)){
    cds <- reduceDimension(cds, max_components = 2, reduction.method = 'DDRTree')
  } else if (opt$regressout == 'pm-umi'){
    cds <- reduceDimension(cds, max_components = 2, reduction.method = 'DDRTree', residualModelFormulaStr='~percent.mito + nUMI')
  } else if (opt$regressout == 'cc') {
    cds <- reduceDimension(cds, max_components = 2, reduction.method = 'DDRTree', residualModelFormulaStr='~S.Score + G2M.Score')
  } else if (opt$regressout == 'pm-umi-cc'){
    cds <- reduceDimension(cds, max_components = 2, reduction.method = 'DDRTree', residualModelFormulaStr='~percent.mito + nUMI + S.Score + G2M.Score')
  }

  ###TRAJECTORY STEP 3: Order cells along the trajectory
  
  print('Ordering cells along trajectory...')
  cds <- orderCells(cds)

  if (is.null(opt$regressout)){
    print(paste("Saving CDS as ", output_name, " in ", opt$outdir, "'.", sep=''))
    saveRDS(cds, output_name)
  } else {
    print(paste('Saving CDS as ', paste(output_name, opt$regressout, sep='-'), ' in ', opt$outdir, sep=''))
    saveRDS(cds, paste(output_name, opt$regressout, sep='-'))
  }
}

if (opt$depots){
  print('RUN ON DEPOTS')
  subtissues <- seurat_object@meta.data$sample_name2
  for (subtissue in unique(subtissues)){
    print(subtissue)
    subset <- SubsetData(seurat_object, cells.use=seurat_object@cell.names[which(subtissues %in% subtissue)])
    output_name <- paste(output_prefix, 'monocle', subtissue, sep='-')
    run_monocle_workflow(subset, output_name)
  }
} else {
  print('RUN ON ALL')
  output_name <- paste(output_prefix, 'monocle', 'downsampled', r.seed, sep='-')
  run_monocle_workflow(seurat_object, output_name)
}
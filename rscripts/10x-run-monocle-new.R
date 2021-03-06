.libPaths('/home/cbmr/pytrik/libraries/')
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/scripts-10x-analysis/rscripts/')

r.seed <- 11

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run Monocle on.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory.'),
  make_option(c('-r', '--regressout'), type='character', help="OPTIONAL. Choose from 'pm-umi', 'pm-umi-cc' or 'cc'.", default='none'),
  make_option(c('-g', '--genelist'), type='character', help='Path to a list of genes to use for ordering.')
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$file) || is.null(opt$outdir) || is.null(opt$genelist)){
  print_help(optionparser)
  quit()  
}

if (!(opt$regressout == 'none')){
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

output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]

run_monocle_workflow <- function(data, output_name){
  
  print('Converting Seurat object to CDS...')
  cds <- importCDS(data)

  print('Estimating size factors and dispersions...')
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  ###TRAJECTORY STEP 1: Choose genes that define a cell's process
  # There are a few ways to do this:
  # - For time series data: Do a differential gene test between cells collected at
  #   the beginning of the experiment and those collected at the end. Did this the
  #   first time, but only the general adipogenesis genes are found then.
  # - For non time series data: 
  #     - (Used here) Use differentially exressed genes between clusters for ordering. 
  #       Monocle's dpFeature does this, but this takes an extremly long time so I 
  #       obtained the DE genes for cluster res.1.5 with Seurat. 
  #     - Select genes with high dispersion across cells. 
  #     - Use known marker genes.

  if (is.null(opt$genelist)){
    print('Selecting genes for ordering with high dispersion...')
    disp_table <- dispersionTable(cds)
    ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    print(paste('Nr of genes selected:', length(ordering_genes)))
  } else {
    ordering_genes <- read.table(opt$genelist)$V1
  }

  print(paste('Nr of genes:', length(ordering_genes)))
  
  print('Seting ordering filter...')
  cds <- setOrderingFilter(cds, ordering_genes)

  ###TRAJECTORY STEP 2: Reduce data dimensionality

  print('Reducing data to two dimensions with DDRTree...')

  if (opt$regressout == 'none'){
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
  
  ###SAVE CDS
  
  setwd(opt$outdir)
  
  if (opt$regressout == 'none'){
    print(paste("Saving CDS as ", output_name, " in ", opt$outdir, "'.", sep=''))
    saveRDS(cds, output_name)
  } else {
    print(paste('Saving CDS as ', paste(output_name, opt$regressout, sep='-'), ' in ', opt$outdir, sep=''))
    saveRDS(cds, paste(output_name, opt$regressout, sep='-'))
  }
}

#############################################

run_monocle_workflow(seurat_object, paste(output_prefix, 'monocle', sep='-'))

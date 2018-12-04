
setwd('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/scripts-10x-analysis/rscripts/')

r.seed <- 66

library(optparse)

option_list <- list(
  make_option(c('-f', '--file'), type='character', help='Path to the dataset to run Monocle on.'),
  make_option(c('-o', '--outdir'), type='character', help='Output directory.'),
  make_option(c('-d', '--depots'), action='store_true', default=F, help='OPTIONAL. If flag -d is given, Monocle is run on individual depots.'),
  make_option(c('-r', '--regressout'), type='character', help="OPTIONAL. Choose from 'pm-umi', 'pm-umi-cc' or 'cc'.", default='none'),
  
  #make_option(c('-s', '--samplingdown'), type='character', help='OPTIONAL. Integer specifying the number of cells per sample to downsample on.')
  #add option for downsampling
  #change -d to -c: column to subset on. 
)

optionparser <- OptionParser(option_list=option_list)
opt <- parse_args(optionparser)

if (is.null(opt$file) || is.null(opt$outdir)){
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

library(Seurat)
library(monocle)

print('LOADING SEURAT OBJECT')
seurat_object <- readRDS(opt$file)

if(!is.null(opt$samplingdown)){
  seurat_object <- SubsetData(SetAllIdent(seurat_object, id='sample_name'), max.cells.per.ident=1000, random.seed=r.seed)
}

setwd(opt$outdir)
output_prefix <- unlist(strsplit(opt$file, '/'))[length(unlist(strsplit(opt$file, '/')))]

run_monocle_workflow <- function(data, output_name){
  
  print('Converting Seurat object to CDS...')
  cds <- importCDS(data)
  
  #if column sample_name exists, replace with samplename, otherwise Monocle will complain about duplicate column names.
  colnames(pData(data)) <- replace(colnames(pData(data)), colnames(pData(data)) =='sample_name', 'samplename')

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
  #       use subtissue as cluster. 
  #     - Select genes with high dispersion across cells (used below). However, what 
  #       I think might happen now is that those are mostly cell cycle genes. We can 
  #       regress out the cell cycle effects when building the trajectory, but it will
  #       still use the same genes for ordering. 
  #     - Use known marker genes.

  print('Selecting genes for ordering with high dispersion...')
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  print(paste('Nr of genes selected:', length(ordering_genes)))
  
  if (opt$regressout == 'cc' || opt$regressout == 'pm-umi-cc'){
    #remove cell cycle genes from the ordering_genes list if cell cycle effects will be regressed out. 
    cc.genes <- readLines('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/downloads/regev_lab_cell_cycle_genes.txt')
    ordering_genes <- ordering_genes[! ordering_genes %in% cc.genes]
  }

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
  
  if (opt$regressout == 'none'){
    print(paste("Saving CDS as ", output_name, " in ", opt$outdir, "'.", sep=''))
    saveRDS(cds, output_name)
  } else {
    print(paste('Saving CDS as ', paste(output_name, opt$regressout, sep='-'), ' in ', opt$outdir, sep=''))
    saveRDS(cds, paste(output_name, opt$regressout, sep='-'))
  }
}

#############################################

if (opt$depots){
  print('RUN ON DEPOTS')
  depots <- seurat_object@meta.data$depot
  for (d in unique(depots)){
    print(d)
    subset <- SubsetData(seurat_object, cells.use=seurat_object@cell.names[which(subtissues == d)])
    if(is.null(opt$samplingdown)){
      output_name <- paste(output_prefix, 'monocle', d, sep='-')
    } else {
      output_name <- paste(output_prefix, 'monocle', 'downsampled', r.seed, d, sep='-')
    }
    run_monocle_workflow(subset, output_name)
  }
} else {
  print('RUN ON ALL')
  if(is.null(opt$samplingdown)){
    output_name <- paste(output_prefix, 'monocle', sep='-')
  } else {
    output_name <- paste(output_prefix, 'monocle', 'downsampled', r.seed, sep='-')
  }
  run_monocle_workflow(seurat_object, output_name)
}
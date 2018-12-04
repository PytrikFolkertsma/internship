
data <- readRDS('../output/10x-180831')

getDemuxletForSample <- function(i, outputdir){
  demuxlet <- read.table(paste('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/demuxlet/demuxlet_out/', outputdir, '/180831_10x_s', i, '.best', sep=''), header=T)

  demuxlet$correct_barcode <- paste(unlist(sapply(strsplit(as.character(demuxlet$BARCODE), '-'), '[[', 1)), '-', i, sep='')
  rownames(demuxlet) <- demuxlet$correct_barcode
  cells <- rownames(data@meta.data)[data@meta.data$timepoint == paste('T', i, sep='')]
  demuxlet_filtered <- demuxlet[demuxlet$correct_barcode %in% cells, ]
  
  num_sng_dbl <- sapply(strsplit(as.character(demuxlet_filtered$BEST), '-'), '[[', 1)

  print(paste('SAMPLE', i))
  print(table(num_sng_dbl))
  return(demuxlet_filtered)
}

getAllDemuxletResults <- function(outdir){
  demuxlet_list <- list()
  for (i in 1:5){
    demuxlet <- getDemuxletForSample(i, outdir)
    demuxlet_list[[i]] <- demuxlet
  }
  demuxlet_all <- do.call(rbind, unname(demuxlet_list))
  Viedemuxlet_all$label <- as.character(demuxlet_all$BEST)
  demuxlet_all$label[startsWith(demuxlet_all$label, 'DBL')] <- 'NA'
  demuxlet_all$label[startsWith(demuxlet_all$label, 'AMB')] <- 'NA'

  data <- AddMetaData(data, demuxlet_all['label'], col.name='label')

  print(dim(demuxlet_all))
  
  plot_grid(
    DimPlot(data, cells.highlight = rownames(data@meta.data)[grep('44B', data@meta.data$label)], reduction.use='tsne', cols.use='gray', cols.highlight = 'blue'),
    DimPlot(data, cells.highlight = rownames(data@meta.data)[grep('BAT14', data@meta.data$label)], reduction.use='tsne', cols.use='gray', cols.highlight = 'blue'),
    DimPlot(data, cells.highlight = rownames(data@meta.data)[grep('1AF', data@meta.data$label)], reduction.use='tsne', cols.use='gray', cols.highlight = 'blue'),
    DimPlot(data, cells.highlight = rownames(data@meta.data)[grep('13a', data@meta.data$label)], reduction.use='tsne', cols.use='gray', cols.highlight = 'blue'),
    labels=c('44B', 'BAT14', '1AF', '13a')
  )
}

getAllDemuxletResults('demuxlet_vcf_plink_bed-updated_sorted')
getAllDemuxletResults('demuxlet_vcf-imputed-genotypes')
getAllDemuxletResults('demuxlet_vcf-imputed-genotypes-alpha0.1-0.5')
getAllDemuxletResults('demuxlet_vcf_tutorial')


#SNG-13a_13a
#SNG-1AF_1AF
#SNG-44B_44B
#SNG-BAT14_BAT14

#SNG-H1_44B
#SNG-H2_BAT14
#SNG-H3_1AF
#SNG-H4_13a


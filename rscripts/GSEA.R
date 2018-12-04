library(gProfileR)

GSEA_BEAM <- function(x){
  
  all_gsea <- list()
  
  for (i in 1:x){
    print(i)
    gene_symbols <- read.table(paste('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/BEAM_results/', x, 'clusters/', 'cluster', i, '.txt', sep=''), header=T)[,1]
    gsea <- gprofiler(as.vector(gene_symbols), organism='hsapiens', significant=T, ordered_query = F, src_filter='GO')
    gsea$cluster <- i
    gsea <- gsea[order(gsea$p.value),]
    all_gsea[[i]] <- gsea
    write.table(data.frame(gsea$term.id, gsea$p.value), paste('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/BEAM_results/BEAM_GSEA/', x, 'clusters/revigo_list_cluster', i, '.txt', sep=''), sep='\t', row.names=F, quote=F)
  }
  
  df.all_gsea <- bind_rows(all_gsea)
  df.all_gsea <- df.all_gsea[,c(ncol(df.all_gsea),1:(ncol(df.all_gsea)-1))]
  write.table(df.all_gsea, paste('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/BEAM_results/BEAM_GSEA/', x, 'clusters/BEAM_', x, 'clusters_GSEAresults.txt', sep=''), sep='\t', quote=F)
}

GSEA_BEAM(6)
GSEA_BEAM(8)
GSEA_BEAM(10)

GSEA_BEAM('6clusters')
GSEA_BEAM('8clusters')
GSEA_BEAM('10clusters')

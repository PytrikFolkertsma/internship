library(gProfileR)
library(biomaRt)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

files <- dir(path = '/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/wgcna-all/tables/', pattern="module_genes.csv", full.names = T, recursive = F, ignore.case = T)

all_gsea <- list()

for (i in 1:length(files)){
  filename <- unlist(strsplit(files[i], '/'))[length(unlist(strsplit(files[i], '/')))]
  modulename <- unlist(strsplit(filename, '_'))[3]
  gene_ids <- read.table(files[i], header=T)[,1]
  gene_symbols <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = gene_ids, mart = mart)
  #take the first 50 genes
  gsea <- gprofiler(gene_symbols[,2][1:50], organism='hsapiens')[1:10,]
  all_gsea[[modulename]] <- gsea
}

df.all_gsea <- bind_rows(all_gsea, .id='module')



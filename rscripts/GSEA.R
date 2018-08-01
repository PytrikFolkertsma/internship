library(gProfileR)
library(biomaRt)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
seagreen_genes <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/wgcna-all/tables/10x-adipocytes_SeuratProject_seagreen_module_genes.csv', header=T)[,1]



gene_list <- getBM(filters="ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=seagreen_genes,mart=mart)
gsea <- gprofiler(gene_list[,2], organism='hsapiens')


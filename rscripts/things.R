



cells_ebf2 <- WhichCells(all10x, subset.name = "EBF2", accept.low = 0.7)

all10x <- AddMetaData(all10x, data.frame(row.names=cells_ebf2, ebf2=rep('ebf2', length(cells_ebf2))))

VlnPlot(all10x, group.by='ebf2', features.plot=c('EBF2'), point.size.use =0.1, x.lab.rot=T)

plot_grid(plot_cell_trajectory(peri, color_by='samplename', cell_size=0.5),
          plot_cell_trajectory(visce, color_by='samplename', cell_size=0.5),
          plot_cell_trajectory(subq, color_by='samplename', cell_size=0.5),
          plot_cell_trajectory(supra, color_by='samplename', cell_size=0.5),
          labels = c('Peri', 'Visce', 'Subq', 'Supra'))

plot_grid(plot_cell_trajectory(peri, color_by='Phase', cell_size=0.5),
          plot_cell_trajectory(visce, color_by='Phase', cell_size=0.5),
          plot_cell_trajectory(subq, color_by='Phase', cell_size=0.5),
          plot_cell_trajectory(supra, color_by='Phase', cell_size=0.5),
          labels = c('Peri', 'Visce', 'Subq', 'Supra'))

plot_grid(plot_cell_trajectory(peri, color_by='ebf2', cell_size=0.5),
          plot_cell_trajectory(visce, color_by='ebf2', cell_size=0.5),
          plot_cell_trajectory(subq, color_by='ebf2', cell_size=0.5),
          plot_cell_trajectory(supra, color_by='ebf2', cell_size=0.5),
          labels = c('Peri', 'Visce', 'Subq', 'Supra'))


matrix<-subset0_STIM@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["gene of interest",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})

library(data.table)
demuxlet <- fread("/projects/pytrik/demuxlet.tutorial/out.best")




#get genes for monocle ordering
markers_T1T2T3_res1 <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-noreg-T1T2T3_res.1.5_negbinom', header=T)
markers_T4T5_res1 <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-noreg-T4T5_res.1.5_negbinom', header=T)

markers_T1T2T3_res1 <- markers_T1T2T3_res1[markers_T1T2T3_res1$p_val_adj < 0.05,]
markers_T4T5_res1 <- markers_T4T5_res1[markers_T4T5_res1$p_val_adj < 0.05,]

length(markers_T1T2T3_res1$gene) #4552 genes
length(unique(markers_T1T2T3_res1$gene)) #1726 unique genes

length(markers_T4T5_res1$gene) #4906 genes
length(unique(markers_T4T5_res1$gene)) #1579 unique genes

length(intersect(markers_T4T5_res1$gene, markers_T1T2T3_res1$gene)) #843 shared genes

all_genes <- unique(c(as.vector(markers_T4T5_res1$gene), as.vector(markers_T1T2T3_res1$gene)))
length(all_genes) #2462

write.table(all_genes, col.names=F, row.names=F, file='/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-noreg-T1T2T3-T4T5_res.1.5_negbinom_filtered-pval_genes-only')
#write.table(all_genes, col.names=F, row.names=F, file='/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T1T2T3-T4T5_res.1_negbinom_filtered-pval_genes-only')
####

#Do the same, but filter the DE gene table for the top 10/20 genes per cluster.

top50 <- markers_T1T2T3_res1 %>% arrange(desc(avg_logFC)) %>% group_by(cluster) %>% top_n(50, avg_logFC) %>% arrange(cluster)
length(unique(top50$gene)) #463 genes

top50_2 <- markers_T4T5_res1 %>% arrange(desc(avg_logFC)) %>% group_by(cluster) %>% top_n(50, avg_logFC) %>% arrange(cluster)
length(unique(top50_2$gene)) #420 genes

#166 intersection genes. So adding them together should give 717 unique genes.

top50_combined <- unique(c(as.vector(top50$gene), as.vector(top50_2$gene)))

write.table(top50_combined, col.names=F, row.names=F, file='/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T1T2T3-T4T5_res.1_negbinom_filtered-pval_top50logFC-per-cluster_genes-only', quote=F)


##########################################

#DE genes T1T2T3, T4, T5 combined (res1)
markers_T1T2T3_res1 <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T1T2T3_res.1_negbinom', header=T)
T4_markers_res1 <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T4_res.1_negbinom', header=T)
T5_markers_res1 <- read.table('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T5_res.1_negbinom', header=T)

markers_T1T2T3_res1 <- markers_T1T2T3_res1[markers_T1T2T3_res1$p_val_adj < 0.05,]
T4_markers_res1 <- T4_markers_res1[T4_markers_res1$p_val_adj < 0.05,]
T5_markers_res1 <- T5_markers_res1[T5_markers_res1$p_val_adj < 0.05,]

all_genes <- unique(c(as.vector(markers_T1T2T3_res1$gene), as.vector(T4_markers_res1$gene), as.vector(T5_markers_res1$gene)))
write.table(all_genes, col.names=F, row.names=F, file='/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/markergenes/180831/markers_10x-180831-T1T2T3-T4-T5_res.1_negbinom_filtered-pval_genes-only', quote=F)

###############

t <- as.data.frame(table(discarded.cells@meta.data$sample_name))
barplot(t$Freq, names.arg=t$Var1, las=2, main='Discarded cells from alignment')



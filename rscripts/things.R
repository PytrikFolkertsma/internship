



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




load('/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/10x.RData')

n_cells <- length(rownames(all10x@meta.data))
cols <- c('sample_name', 'sample_name2', 'Phase')

#get each sample name, subtissue and cell cycle phase
cols_dummies <- c()
for (i in 1:length(cols)) {
  cols_dummies <- c(cols_dummies, as.vector(unique(unlist(all10x@meta.data[cols[i]]))))
}

#remove Supra4 and Subq4
#cols_dummies <- cols_dummies[!cols_dummies %in% c('Supra_4', 'Subq_4')]

#create dataframe, every sample name, subtissue and cell cycle phase is a column
DF <- as.data.frame(matrix(NA, ncol = length(cols_dummies), nrow = n_cells), col.names=cols_dummies, row.names=rownames(seurat_obj@meta.data))
colnames(DF) <- cols_dummies 

#fill dataframe. 
for (col in cols){
  for (v in as.vector(unique(unlist(all10x@meta.data[col])))){
    print(v)
    DF[,v] <- unlist(lapply(as.vector(unlist(all10x@meta.data[col])), function(x){
      if (x == v){
        return(1)
      } else {
        return(0)
      }
    }))
  }
}

brown <- unlist(lapply(all10x@meta.data$sample_name2, function(x){
  if (x == 'Supra' || x == 'Peri'){
    return(1)
  } else {
    return(0)
  }
}))

white <- unlist(lapply(all10x@meta.data$sample_name2, function(x){
  if (x == 'Visce' || x == 'Subq'){
    return(1)
  } else {
    return(0)
  }
}))

cluster12 <- unlist(lapply(all10x@meta.data$res.0.5, function(x){
  if (x == 12){
    return(1)
  } else {
    return(0)
  }
}))

DF['brown'] <- brown
DF['white'] <- white
DF['cluster12'] <- cluster12

all10x <- AddMetaData(all10x, DF)

save(all10x, file='/projects/pytrik/sc_adipose/analyze_10x_fluidigm/data/10x.RData')

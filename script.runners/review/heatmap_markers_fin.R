library(gplots)
library(RColorBrewer)

load("./r.data.files/second_proc/sub.list.fin_genes.rda") # limma_fin.r for processing
dropSamples(sub.list, which(sub.list$groups==0)) -> sub.list
orderBy(sub.list, "groups") -> sub.list
markers <- c("MKI67", "TNFRSF8", "TIA1", "CD4", "CXCL13", "ICOS", "IL17A",
             "CD163", "GZMB", "TBX21", "CXCR3", "GATA3", "LEF1", "CCR4",
             "CXCL12", "CD8A", "CTGF", "IL4", "CXCL5", "CPA3", "MZB1",
             "FAM46C", "KCNN3", "AICDA", "PRF1")
sub.list$exprs.mat[which(rownames(sub.list$exprs.mat) %in% markers),] -> df

# scale and center
t(df) -> df
apply(df, 2, function(y) y - mean(y)) -> df
t(df) -> df

# arrange df 
my.colors <- colorRampPalette(colorRampPalette(brewer.pal(11,"RdBu")[-c(4,5,7,8)])(50))
all.colors <- sub.list$groups
all.colors[which(sub.list$groups %in% 1)] <- "red"
all.colors[which(sub.list$groups %in% 2)] <- "cyan"
all.colors[which(sub.list$groups %in% 3)] <- "magenta"
all.colors[which(sub.list$groups %in% 4)] <- "blue"

colnames(df) <- rep("", ncol(df))
heatmap.2(df, trace="none", margins=c(10,10),
          col=my.colors(50)[50:1], yaxt="n",
          Colv=F, ColSideColors=all.colors)

df.bin <- df
apply(df.bin, 1, function(x){
  bin <- rep(0, length(x))
  bin[which(x > median(x)+sd(x))] <- 1
  return(bin)
}) -> bin
t(bin) -> df.bin

colnames(df.bin) <- rep("", ncol(df.bin))
markers <- c("MKI67", "TNFRSF8", "TIA1", "CD4", "CXCL13", "ICOS", "IL17A",
             "CD163", "GZMB", "TBX21", "CXCR3", "GATA3", "LEF1", "CCR4",
             "CXCL12", "CD8A", "CTGF", "IL4", "CXCL5", "CPA3", "MZB1",
             "FAM46C", "KCNN3", "AICDA", "PRF1")
c(2,1,3,2,4,2,2,3,3,3,3,1,1,1,1,4,4,4,4,4,4,2,2,2,3) -> marker.labs

fin.ord <- c()
for (i in 1:4){
  c(fin.ord, which(rownames(df.bin) %in% markers[which(marker.labs %in% i)])) -> fin.ord
}
df.bin[fin.ord,] -> df.bin

col.pos <- c("red", "cyan", "magenta", "blue")
row.cols <- rep(NA, nrow(df.bin))
for (i in 1:nrow(df.bin)){
  row.cols[i] <- col.pos[marker.labs[which(markers %in% rownames(df.bin)[i])]]
}
heatmap.2(df.bin, trace="none", margins=c(10,10),
          col=c("white", "darkgrey"), yaxt="n",
          Colv=F, Rowv=F, ColSideColors=all.colors,
          RowSideColors = row.cols)

# create heatmap for iqbal data
load("./subtyping/r.data.files/external/sub.list.fin_genes.rda")

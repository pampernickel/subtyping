library(ggplot2)
library(rgl)

# load("./r.data.files/second_proc/sub.list.normalized.adjusted.combat.cov_QC.noSex.mapped.symbol.entrezID.w.groups.Rdata")
load("./r.data.files/second_proc/sub.list.fin_genes.rda")
df <- sub.list$exprs.mat
apply(df, 1, var) -> vars
which(vars > 0.25) -> ind
df[ind,] -> df

pca = prcomp(df,scale=T)
x= pca$rotation[,1]
y= pca$rotation[,2]
z= pca$rotation[,3]

all.colors <- sub.list$groups
all.colors[which(sub.list$groups %in% 0)] <- "white"
all.colors[which(sub.list$groups %in% 1)] <- "red"
all.colors[which(sub.list$groups %in% 2)] <- "cyan"
all.colors[which(sub.list$groups %in% 3)] <- "magenta"
all.colors[which(sub.list$groups %in% 4)] <- "blue"
plot3d(x,y,z,
       col=all.colors,
       size=2.0, type='s', labcex=1.2, alpha=0.7)

# create plot with tfh demarcated
tfh_annot <- read.csv("./annotations/tfh_annot.csv")
tfh_annot$ID[intersect(which(tfh_annot$AITL.PTCL_NOS.PTCL.F %in% c("PTCL-NOS", "PTCL-F")), 
                       which(tfh_annot$TFH %in% 1))] -> tfh

all.colors <- sub.list$groups
#all.colors[which(sub.list$groups %in% 0)] <- "white"
#all.colors[which(sub.list$groups %in% 1)] <- "red"
#all.colors[which(sub.list$groups %in% 2)] <- "cyan"
#all.colors[which(sub.list$groups %in% 3)] <- "magenta"
#all.colors[which(sub.list$groups %in% 4)] <- "blue"

all.colors[which(sub.list$color.ind %in% "red")] <- "cornflowerblue"
all.colors[which(colnames(sub.list$exprs.mat) %in% "PTCL_NOS")] <- "red"
all.colors[which(sub.list$tenomic %in% tfh)] <- "black"
plot3d(x,y,z,
       col=all.colors,
       size=2.0, type='s', labcex=1.2, alpha=0.7)

# pca of external data set
# load("./r.data.files/external/sub.list.batch1.Rdata")
# load("./r.data.files/external/cores.Iqbal.Rdata")

# load("./r.data.files/external/results.resICL.Rdata")
# sub.list.batch1$groups <- rep(NA, ncol(sub.list.batch1$exprs))
# resICL[[2]]$itemConsensus[which(resICL[[2]]$itemConsensus$k %in% 4),] -> groups
# 
# for (i in 1:ncol(sub.list.batch1$exprs)){
#   groups[which(groups$item %in% colnames(sub.list.batch1$exprs)[i]),] -> loc
#   loc$cluster[which(loc$itemConsensus %in% max(loc$itemConsensus))] -> sub.list.batch1$groups[i]
# }

# sub.list.batch1$groups <- rep(0, ncol(sub.list.batch1$exprs))
# for (i in 1:length(cores.Iqbal)){
#  sub.list.batch1$groups[which(colnames(sub.list.batch1$exprs) %in% cores.Iqbal[[i]])] <- i
#}

load("./r.data.files/external/sub.list.fin_genes.rda") # see limma_fin.r for generation
df <- sub.list.batch1$exprs
apply(df, 1, var) -> vars
which(vars > 0.25) -> ind
df[ind,] -> df

pca = prcomp(df,scale=T)
x= pca$rotation[,1]
y= pca$rotation[,2]
z= pca$rotation[,3]
all.colors <- sub.list.batch1$groups
all.colors[which(sub.list.batch1$groups %in% 0)] <- "white"
all.colors[which(sub.list.batch1$groups %in% 1)] <- "red"
all.colors[which(sub.list.batch1$groups %in% 2)] <- "cyan"
all.colors[which(sub.list.batch1$groups %in% 3)] <- "magenta"
all.colors[which(sub.list.batch1$groups %in% 4)] <- "blue"
plot3d(x,y,z,
       col=all.colors,
       size=2.0, type='s', labcex=1.2, alpha=0.7)
# :::
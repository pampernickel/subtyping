`%ni%` = Negate(`%in%`)

# TENOMIC data
load("./r.data.files/second_proc/sub.list.normalized.adjusted.combat.cov_QC.noSex.mapped.symbol.entrezID.w.groups.Rdata")
dim(sub.list$exprs.mat)
table(colnames(sub.list$exprs.mat))
sub.list$tfh.stat[which(sub.list$color.ind %in% "red")] <- 2
orderBy(sub.list, "tfh.stat") -> sub.list
sub.list$tenomic[which(sub.list$color.ind %in% "red")] -> aitls
sub.list$tenomic[which(sub.list$color.ind %ni% "red")] -> o

# check how tfh-like are distributed across subgroups
tfh_annot <- read.csv("./annotations/tfh_annot.csv")
tfh_annot$ID[intersect(which(tfh_annot$AITL.PTCL_NOS.PTCL.F %in% c("PTCL-NOS", "PTCL-F")), 
                       which(tfh_annot$TFH %in% 1))] -> tfh
sub.list$groups[which(sub.list$tenomic %in% tfh)]
sub.list$tenomic[intersect(which(sub.list$tenomic %in% tfh), which(sub.list$groups %in% 0))]

# external data + external data annotation
load("./r.data.files/external/sub.list.batch1.Rdata")
iq <- read.csv("./annotations/filelist.ordered.dates.csv")
iq_annots <- read.csv("./annotations/gse58445.iqbal.csv")
gsub("/export/scratch/npietros//Projects/PTCL-Project/ExternalData/GSE58445/", "", iq$V1) -> iq$V1
dim(sub.list.batch1$exprs) # 54675 x 93
iq[which(iq$V1 %in% colnames(sub.list.batch1$exprs)),] -> iq
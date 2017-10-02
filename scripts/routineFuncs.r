`%ni%` = Negate(`%in%`)

dropSamples <- function(list, ind){
  # list is of the form
  # found in 20130902.full.list.genes.rda
  if (length(ind)==0){
    stop("0 length index. No changes made to your list.", call.=FALSE)
  }
  # --- for each list element, remove the indicated index
  for (i in 1:length(list)){
    if (is.matrix(list[[i]])){
      # --- drop columns
      list[[i]][,-ind] -> list[[i]]
    } else {
      list[[i]][-ind] -> list[[i]]
    }
  }
  
  return(list)
}

subSample <- function(list, ind){
  if (length(ind)==0){
    stop("0 length index. No changes made to your list.", call.=FALSE)
  }
  
  # --- for each list element, remove the indicated index
  for (i in 1:length(list)){
    if (is.matrix(list[[i]])){
      list[[i]][,ind] -> list[[i]]
    } else {
      list[[i]][ind] -> list[[i]]
    }
  }
  return(list)
}

orderBy <- function(list, by=c("exprs.mat", "tenomic", "affy.no", "sorted", "hybrid.date", "tet2", "idh", "dnmt", "tfh.stat", "perc.tumor", "cell.type", "pathos", "centre","subgroup","class","hybrid.date.corrected")){
  vect <- list[[which(names(list) %in% by)]]
  #print(vect)
  for (i in 1:length(list)){
    if (is.matrix(list[[i]])){
      list[[i]][,order(vect)] -> list[[i]]
    } else {
      list[[i]][order(vect)] -> list[[i]]
    }
  }
  return(list)
}

runLimma <- function(df, lab, eBayes=F){
  if (!is.loaded("limma")){
    library(limma)
  }
  
  # create switch for limma where unique(lab) > 2
  toptab.GRP <- NA
  if (length(unique(lab))==2){
    ff <- rep(0,ncol(df))
    for (i in 1:length(unique(lab))){
      ff[which(lab %in% unique(lab)[i])] <- i
    }
    design <- model.matrix(~ -1+factor(ff))
    colnames(design) <-c("GRP1", "GRP2") 
    fit <- lmFit(df, design) 
    contrast.matrix <- makeContrasts( 
      GRP.1.2 = GRP2 - GRP1,
      levels=design
    )
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2, trend=eBayes) 
    toptab.GRP = topTable(fit2, coef=c(1), number = dim(df)[1], adjust="BH")
  } else if (length(unique(lab)) > 2) {
    # case of pairwise comparisons across three groups
    design <- model.matrix(~ 0+factor(lab))
    colnames(design) <- paste("group", unique(lab), sep="")
    combn(paste("group", unique(lab), sep=""), 2) -> pw
    apply(pw, 2, function(x) paste(x[1],"-",x[2],sep="")) -> x
    contrast.matrix <- makeContrasts(x, 
                                     levels=paste("group", unique(lab), sep=""))
    fit <- lmFit(df, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    toptab.GRP = topTable(fit2, coef=c(1), number = dim(df)[1], adjust="BH")
  }
  return(toptab.GRP)
}

getTopTab <- function(res, p.val=.05, mode=c("%", "logFC"), param=""){
  # method for retrieving top hits from limma toptab
  res[which(res$adj.P.Val <= p.val),] -> pres
  if (mode %in% "%" && param != ""){
    # filter based on top x% adj.p.val (get separately for up and dn genes) 
    pres[which(pres$logFC < 0),] -> dn
    pres[which(pres$logFC >= 0),] -> up
    up[1:ceiling(nrow(up)*param),] -> up
    dn[1:ceiling(nrow(dn)*param),] -> dn
    rbind(up, dn) -> press
  } else if (mode %in% "logFC" && param != ""){
    pres[which(abs(pres$logFC) >= param),] -> press
  } else if (mode %ni% c("%", "logFC") ||
             param == ""){
    stop("Missing one ore more input parameters (mode/param).")
  }
  return(press)
}

setLimmaRun <- function(mat, g, e=""){
  # mat: exprs.mat
  # g: list containing columns of indices of mat; the length of g corresponds
  # with the number of subgroups; in the case of length(g) == 2, this results
  # in a pairwise comparison; in the case of length(g) > 2, the default mode
  # would be of the form x vs. rest; alternative mode to follow;
  # e: exclusion list; should be of the same length as g
  if (!is.loaded("limma")) require(limma)
  lapply(1:length(g), function(x){
    lab <- rep(0, ncol(mat))
    lab[g[[x]]] <- 1
    if (e %ni% "" && length(g) == length(e)){
      lab[e[[x]]] <- NA
    } else {
      stop("Exclusion list should have the same length as the group argument.")
    }
    mat[,which(lab %in% c(0,1))] -> mat
    lab[which(lab %in% c(0,1))] -> lab
    return(runLimma(mat, lab))
  }) -> limma.res
  return(limma.res)
}

compareLimmaRes <- function(l1, l2, nn){
  if (length(l1) != length(l2)) stop("Lists do not have equal lengths.")
  lapply(1:length(l1), function(x){
    if (nrow(l1[[x]]) != nrow(l2[[x]])) stop("Data frames do not have an equal number of rows.")
    l1[[x]][order(rownames(l1[[x]])),which(colnames(l1[[x]]) %in% nn)] -> l1x
    l2[[x]][order(rownames(l2[[x]])),which(colnames(l2[[x]]) %in% nn)] -> l2x
    cbind.data.frame(l1x,l2x,rep(paste("Subgroup" , x, "vs. REST", sep=" "), length(l1x))) -> t
    return(t)
  }) -> res
  do.call("rbind", res) -> res
  colnames(res) <- c("V1", "V2", "comparison")
  return(res)
}

plotComp <- function(res, xlab, ylab){
  # plot function for compareLimmaRes() results
  if (!is.loaded("ggplot2")) require(ggplot2)
  ggplot(res, aes(x=V1, y=V2))+
    geom_point(size=0.5)+facet_wrap(~comparison)+
    theme_bw()+xlab(xlab)+ylab(ylab)+
    theme(axis.title = element_text(size=18),
          axis.text=element_text(size=16),
          strip.text=element_text(size=17))
}

# ::: functions for RNASeq (just-in-case)
createLabs <- function(df, labs){
  colData <- matrix("", nrow=ncol(df), ncol=2)
  colnames(colData) <- c("condition", "type")
  colData[,2] <- rep("paired-end", nrow(colData))
  colData[which(labs %in% 1),1] <- "set2"
  colData[which(labs %in% 0),1] <- "set1"
  as.data.frame(colData) -> colData
  return(colData)
}

runDESeq <- function(df, colData){
  if (!is.loaded("DESeq2")){
    library(DESeq2)
  }
  dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = colData,
                                design = ~ condition)
  featureData <- data.frame(gene=rownames(df))
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds$condition <- factor(dds$condition, levels=c("set1","set2"))
  dds$condition <- relevel(dds$condition, ref="set1")
  dds$condition <- droplevels(dds$condition)
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}

# ::: other useful functions for checks
plotGenes <- function(mat, g, gene, type){
  # plots gene expression across groups
  # mat: expression matrix
  # g: groups
  # gene: list of genes of interest
  if (!is.loaded("ggplot2")) require(ggplot2)
  if (!is.loaded("reshape2")) require(reshape2)
  if (!is.loaded("gplots")) require(gplots)
  if (!is.loaded("RColorBrewer")) require(RColorBrewer)
  
  if (type == "boxplot"){
    lapply(1:length(g), function(x){
      melt(mat[which(rownames(mat) %in% gene),g[[x]]]) -> t
      cbind.data.frame(t, rep(paste("Group",x,sep="_"), nrow(t))) -> df
      colnames(df) <- c(colnames(df)[1:(ncol(df)-1)], "group")
      return(df)
    }) -> res
    do.call("rbind", res) -> res
    ggplot(res, aes(x=group, y=value))+
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15, height=0.01, alpha=0.5)+
      facet_wrap(~Var1)+
      theme_bw()
  } else if (type == heatmap){
    mat[which(rownames(mat) %in% gene),] -> t
    my.colors <- colorRampPalette(colorRampPalette(brewer.pal(11,"RdBu")[-c(4,5,7,8)])(50))
    
    # create RowSide colors based on g
    heatmap.2(t(t), scale="col", col=my.colors,
              key=TRUE, keysize=0.5, trace="none", 
              sepwidth=0.01, #colsep=1:ncol(assembled), sepcolor='white',
              margins = c(10,10),
              cexRow=1,cexCol=1)
  }
}


# ::: create binarized heatmap
binarizeMarkers <- function(df, groups, sampleNames, markers, sdFactor){
  if (!is.loaded("gplots")) require(gplots)
  sapply(markers$marker, function(x) 
    which(rownames(df) %in% x)) -> rord
  df[rord,] -> df
  
  # currently hard-coded  
  color.hash <- cbind(c(1:4), c("red", "cyan", "magenta", "blue"))
  colnames(color.hash) <- c("group", "color")
  col.colors <- rep(NA, ncol(df))
  for (i in 1:4){
    col.colors[which(groups %in% i)] <- color.hash[which(color.hash[,1] %in% i),2]
  }
  
  row.colors <- rep(NA, nrow(df))
  for (i in 1:4){
    row.colors[which(markers$group %in% i)] <- color.hash[which(color.hash[,1] %in% i),2]
  }
  
  # binarization; currently uses the median as the baseline limit
  # for binarization; the restriction can be adjusted with respect to the
  # sd in a vector
  factor <- sdFactor
  apply(df, 1, function(x)
    median(x)+factor*sd(x)) -> medx
  
  df.bin <- matrix(0, nrow=nrow(df), ncol=ncol(df))
  for (i in 1:nrow(df)){
    df.bin[i,which(df[i,] > medx[i])] <- 1
  }
  colnames(df.bin) <- sampleNames
  rownames(df.bin) <- rownames(df)
  
  # then within each group, order acc. to highest expression of markers
  fin.ord <- c()
  for (i in 1:4){
    df.bin[which(rownames(df.bin) %in% 
                   markers$marker[which(markers$group %in% i)]),
           which(groups %in% i)] -> x
    colSums(x) -> s
    colnames(x)[order(s)] -> n.ord
    c(fin.ord, n.ord) -> fin.ord
  }
  
  as.numeric(sapply(fin.ord, function(x) 
    which(colnames(df.bin) %in% x))) -> fin.ord
  df.bin[,fin.ord] -> df.bin
  
  heatmap.2(df.bin, dendrogram="none", Rowv=F, Colv=F,
            trace="none", col=c("grey95", "grey50"),
            ColSideColors=col.colors,
            RowSideColors=row.colors,
            labCol="", cexRow=1.5, margin=c(7,7),
            colsep=getChangePoint(as.numeric(as.factor(col.colors))), 
            rowsep=getChangePoint(as.numeric(as.factor(row.colors))), sepcolor = "grey12")
}

# ::: heatamp label creation
createColLabs <- function(labs, groups){
  color.hash <- cbind(c(1:4), c("red", "cyan", "magenta", "blue"))
  colnames(color.hash) <- c("group", "color")
  col.colors <- rep(NA, ncol(df))
  for (i in 1:4){
    col.colors[which(groups %in% i)] <- color.hash[which(color.hash[,1] %in% i),2]
  }
  return(col.colors)
}

getChangePoint <- function(d) {
  p <- cumsum(rle(d)$lengths)
  p[-length(p)]
}

# ::: functions for enrichment analyses
runRomer <- function(df, labs, gene.sets, rots){
  # --- params: MsigDB signature set, iset, sub.list, vector for generating design matrix
  if (!is.loaded("limma")) library(limma)
  if (!is.loaded("parallel")) library(parallel)
  ff <- rep(0,ncol(df))
  ff[which(labs %in% 1)] <- 1
  design <- cbind(Intercept=1,Group=labs)
  
  # index: list corresponding to the rows of the gene sets
  # adapt if gene.sets is a list of lists or not
  if (is.list(gene.sets[[1]])){
    lapply(gene.sets, function(x) 
      ids2indices(x, rownames(df))) -> indices
  } else if (is.character(gene.sets[[1]]) || is.numeric(gene.sets[[1]])){
    list(gene.sets) -> gene.sets
    names(gene.sets) <- "GS1"
    lapply(gene.sets, function(x) 
      ids2indices(x, rownames(df))) -> indices
  }
  
  # select signatures of interest: 
  mclapply(indices, function(x)
    romer(y=df, index=x, design=design, set.statistic = "mean", 
          nrot=rots), mc.cores=detectCores()-1) -> rr
  
  lapply(rr, function(x) 
    x[setdiff(c(which(x[,2] <= 0.05),
                which(x[,3] <= 0.05)), which(x[,4] <= 0.05)),]) -> rr.sub
  list(rr, rr.sub) -> res
  names(res) <- c("all", "sub")
  return(res)
}

createGeneSets <- function(dir){
  # specify directory from which to create gene sets; this function version has been
  # designed to handle MSigDB signatures
  list.files(dir, full.names = T)[grep(".txt", list.files(dir, full.names = T))] -> f
  lapply(f, function(x){
    readLines(x) -> ff
    lapply(ff, function(x) unlist(strsplit(x, "\t"), recursive=F)) -> ff
    lapply(ff, function(x) x[1]) -> names(ff)
    lapply(ff, function(x) x[-c(1:2)]) -> names(ff)
  }) -> gene.sets
  sapply(strsplit(f, "\\/"), function(x) x[length(x)]) -> names(gene.sets)
  return(gene.sets)
}

createGSEAset <- function(gene.list, format=c("on-the-fly", "gmt")){
  # from a gene list, create a gene set compatible with the standalone
  # GSEA program; file is a .gmt file; change extension manually
  if (format == "gmt"){
    # has some issues
    max(sapply(gene.list, function(x) length(unique(x)))) -> ml
    df <- matrix("", nrow=ml, ncol=length(gene.list))
    colnames(df) <- names(gene.list)
    for (i in 1:length(gene.list)){
      df[1:length(unique(gene.list[[i]])),i] <- unique(gene.list[[i]])
    }
  } else {
    # for pasting on the `on-the-fly' gene set option
    file.create("./signatures/gene.list.txt")
    sink("./signatures/gene.list.txt")
    for (x in 1:length(gene.list)){
      cat(paste(toupper(names(gene.list)[x]), "\n", sep=""))
      for (i in 1:length(gene.list[[x]])){
        cat(paste(gene.list[[x]][i], "\n", sep=""))
      }
    }
    sink()
  }
  
  
  print("Writing output to signatures/gene_list.csv")
  # write.csv(as.data.frame(df[,1:10]), file="./signatures/gene_list.csv")
}

barPlot <- function(df, genes, labs, mode){
  # create a bar plot based on the mean or median expression of 
  # a list of gene sets
  if (is.list(genes)){
    lapply(genes, function(x){
      lapply(unique(labs), function(y){
        df[which(rownames(df) %in% x),which(labs %in% y)] -> df.sub
        res <- NA
        if (mode %in% "mean"){
          cbind.data.frame(rowMeans(df.sub), rep(y, nrow(df.sub))) -> res
          colnames(res) <- c("m", "l")
        } else if (mode %in% "median"){
          cbind.data.frame(apply(df.sub, 1, median),rep(y, nrow(df.sub))) -> res
          colnames(res) <- c("m", "l")
        }
        return(res)
      }) -> res
      do.call("rbind.data.frame", res) -> res
      return(res)
    }) -> res
  } else if (is.character(genes)){
    lapply(unique(labs), function(y){
      df[which(rownames(df) %in% genes),which(labs %in% y)] -> df.sub
      if (mode %in% "mean"){
        cbind.data.frame(rowMeans(df.sub), rep(y, nrow(df.sub))) -> res
        colnames(res) <- c("m", "l")
      } else if (mode %in% "median"){
        cbind.data.frame(apply(df.sub, 1, median),rep(y, nrow(df.sub))) -> res
        colnames(res) <- c("m", "l")
      }
      return(res)
    }) -> res
    do.call("rbind.data.frame", res) -> res
  }
  return(res)
}

sigPlot <- function(df, genes, labs, mode=c("mean", "median"), by=c("gene", "sample")){
  if (!is.loaded("ggplot2")) require(ggplot2)
  if (by == "gene"){
    # plot the median value for each gene across samples with the same label
    barPlot(df, genes, labs, mode) -> bp
    lapply(1:length(bp), function(x){ 
      cbind.data.frame(bp[[x]], rep(names(genes)[x], nrow(bp[[x]]))) -> r
      colnames(r) <- c(colnames(bp[[x]]), "s")
      return(r)
    }) -> res
    do.call("rbind.data.frame", res) -> res
    as.factor(res$l) -> res$l
    ggplot(res, aes(x=l, y=m))+
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.05, height=0.05, alpha=0.5)+
      facet_wrap(~s)+xlab("Group")+ylab(paste(mode, "expression",sep=" "))+
      theme_bw()+
      theme(axis.title = element_text(size=18),
            axis.text=element_text(size=16),
            strip.text=element_text(size=17))
  } else if (by == "sample"){
    # plot the median value of each signature per sample, then
    # group the samples by label
    sigExpression(df, genes, labs, mode, by="sample") -> res
    ggplot(res, aes(x=lab, y=expression))+
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.05, height=0.05, alpha=0.5)+
      facet_wrap(~signature)+xlab("Group")+ylab(paste(mode, "expression",sep=" "))+
      theme_bw()+
      theme(axis.title = element_text(size=18),
            axis.text=element_text(size=16),
            strip.text=element_text(size=17))
  }
}

sigExpression <- function(df, genes, labs, mode, by){
  # returns mean or median expression of gene set per sample
  # in long format
  lapply(1:length(genes), function(z){
    df[which(rownames(df) %in% genes[[z]]),] -> sub
    if (by %in% "sample"){
      if (length(sub) > length(labs)){
        if (mode == "mean"){
          cbind(labs, rep(names(genes)[z], length(labs)), colMeans(sub)) -> res
        } else {
          cbind(labs, rep(names(genes)[z], length(labs)), 
                apply(sub, 2, function(y) median(y))) -> res
        }
        colnames(res) <- c("lab", "signature", "expression")
      } else {
        cbind(labs, rep(names(genes)[z], length(labs)), sub) -> res
        colnames(res) <- c("lab", "signature", "expression")
      }
    } else {
      # summarize results for each label into a single value
      if (length(sub) > length(labs)){
        if (mode == "mean"){
          sapply(unique(labs), function(x) 
            median(apply(sub[,which(labs %in% unique(labs)[x])], 2, function(y) median(y)))) -> smed
          cbind(unique(labs), names(genes)[z], smed) -> res
          colnames(res) <- c("lab", "signature", "expression")
        } else {
          sapply(unique(labs), function(x) 
            mean(apply(sub[,which(labs %in% unique(labs)[x])], 2, function(y) mean(y)))) -> smed
          cbind(unique(labs), names(genes)[z], smed) -> res
          colnames(res) <- c("lab", "signature", "expression")
        }
      } else {
        if (mode == "median"){
          sapply(unique(labs), function(x) median(sub[which(labs %in% unique(labs)[x])])) -> res
        } else {
          sapply(unique(labs), function(x) mean(sub[which(labs %in% unique(labs)[x])])) -> res
        }
        cbind(unique(labs), names(genes)[z], smed) -> res
        colnames(res) <- c("lab", "signature", "expression")
      }
    }
    return(res)
  }) -> res
  do.call("rbind.data.frame", res) -> res
  as.factor(res$lab) -> res$lab
  as.numeric(as.character(res$expression)) -> res$expression
  return(res)
}

boxPlot <- function(df, genes, labs){
  # create boxplot of multiple genes
  lapply(genes, function(x){
    lapply(unique(labs), function(y){
      t(rbind(df[which(rownames(df) %in% x), which(labs %in% y)], 
              rep(y, length(which(labs %in% y))))) -> temp
    }) -> temp
    do.call("rbind", temp) -> temps
    cbind(rep(x, nrow(temps)), temps) -> temp
    colnames(temp) <- c("genes", "exp.level", "group")
    return(temp)
  }) -> res
  do.call("rbind.data.frame", res) -> res
  as.numeric(as.character(res$exp.level)) -> res$exp.level
  
  # then check if expression levels are counts or not
  if (range(res$exp.level)[2]-range(res$exp.level)[1] > 100){
    log2(res$exp.level+1) -> res$exp.level
  }
  
  ggplot(res, aes(x=group, y=exp.level))+
    geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.05, height=0.05, alpha=0.5)+
    facet_wrap(~genes)+xlab("Group")+ylab("log2 expression level")+
    theme_bw()+
    theme(axis.title = element_text(size=18),
          axis.text=element_text(size=16),
          strip.text=element_text(size=17))
}

createHeatmap <- function(df, labs, gene.sets, metric=c("cor", "t")){
  # create a heatmap based on contents of a list "gene.sets"
  # where the order of the samples and genes are held constant
  # based on the labels and gene sets
  if (!is.loaded("gplots")) library(gplots)
  if (!is.loaded("RColorBrewer")) library(RColorBrewer)
  
  getLeadingEdge(df, labs, gene.sets, metric) -> le
  unique(unlist(le)) -> all.genes
  df[which(rownames(df) %in% all.genes),] -> df
  
  # color by pathway
  color.sig <- data.frame(names(gene.sets), color=brewer.pal(length(names(gene.sets)), "Spectral"))
  as.character(sapply(rownames(df), function(x)
    names(gene.sets)[which(sapply(gene.sets, function(y) ifelse(x %in% y, T, F)) %in% T)[1]])) -> sig.a
  as.character(merge(as.data.frame(sig.a), color.sig, by.x="sig.a", by.y="names.gene.sets.")[,2]) -> row.colors
  
  col.colors <- rep(NA, ncol(df))
  r.color.sig <- brewer.pal(length(unique(labs)), "Paired")
  for (i in 1:length(unique(labs))){
    col.colors[which(labs %in% unique(labs)[i])] <- r.color.sig[i]
  }
  
  my.colors <- colorRampPalette(colorRampPalette(brewer.pal(11,"RdBu")[-c(4,5,7,8)])(25))
  heatmap.2(as.matrix(df), dendrogram="none", Rowv=F, Colv=F,
            trace="none", col=my.colors(50)[50:1],
            ColSideColors=col.colors,
            RowSideColors=row.colors,cexRow=1.5, margin=c(7,7),
            scale="row")
}

getLeadingEdge <- function(df, labs, gene.sets, metric){
  # use either a simple test (cor/t/ranks) to approximate a leading edge
  # function; in the case of the rank mode, this function yields
  # both the leading edge per gene set, and implicitly, the enrichment
  # per contrast of the form x vs REST
  if (!is.loaded("reshape2")) library("reshape2")
  if (metric %in% c("t", "cor")){
    lapply(gene.sets, function(x){
      le <- list()
      df[which(rownames(df) %in% x),] -> sub
      
        apply(sub, 1, function(y){
          res <- NA
          if (var(as.numeric(y)) > 0 && metric == "t"){
            t.test(as.numeric(y[which(labs %in% 0)]), 
                   as.numeric(y[which(labs %in% 1)]))$p.value -> res
          } else if (var(as.numeric(y)) > 0 && metric == "cor"){
            cor.test(as.numeric(y), labs, method="s")$p.value -> res
          }
          return(res)
        }) -> res
        names(res)[which(res <= 0.05)]  -> le
        return(le)
    }) -> le
  } else if (metric %in% "rank"){
    # get ranks of genes in the list compared to all other
    # genes; in x vs REST mode; get those with a tendency
    # towards enrichment in the group of interest
    lapply(unique(labs), function(y){
      labs.mod <- rep(0, length(labs))
      labs.mod[which(labs %in% y)] <- 1
      runLimma(df, labs.mod) -> lr
      lr[rev(order(lr$logFC)),] -> lr
      lapply(genes, function(z){
        which(rownames(lr) %in% z) -> ranks
        if (median(ranks) < quantile(1:nrow(lr))[2]){
          # return genes which rank highest
          rownames(lr)[ranks][which(ranks <= median(ranks))] -> res
        }
      }) -> le
    }) -> le
    lapply(le, function(x) 
      x[which(sapply(x, function(y) length(y)) > 0)]) -> le
    names(le) <- unique(labs)
  }
  return(le)
}

visSigOverlaps <- function(sub, mode=c("UP", "DN")){
  # ::: Function to get pairwise overlap between lists of gene signatures; 
  # ::: group UP with UP and DN with DN; those without any up or down are grouped with UP
  c(grep("UP", names(sub)), setdiff(1:length(names(sub)), 
                                    c(grep("UP", names(sub)), grep("DN", names(sub))))) -> up
  grep("DN", names(sub)) -> dn
  sub[up] -> sub.up
  sub[dn] -> sub.dn  
  
  df <- NA
  if (mode == "UP"){
    df <- matrix(NA, nrow=length(up), ncol=length(up))
    colnames(df) <- rownames(df) <- names(sub.up)
  } else {
    df <- matrix(NA, nrow=length(dn), ncol=length(dn))
    colnames(df) <- rownames(df) <- names(sub.dn)
  }
  
  for (i in 1:nrow(df)){
    for (j in 1:ncol(df)){
      if (i != j){
        length(intersect(sub[[which(names(sub) %in% colnames(df)[i])]],
                         sub[[which(names(sub) %in% colnames(df)[j])]]))/
          max(length(sub[[which(names(sub) %in% colnames(df)[i])]]),
              length(sub[[which(names(sub) %in% colnames(df)[j])]]))-> df[i,j]
      }
    }
  }
  
  return(df)
}

# ::: signature creation (*.grp format, GSEA-compatible) from limma results (need to implement 
# alternative for results of DESeq)
subsetLimmaRes <- function(res, perc=.05){
  res[which(res$logFC < 0),] -> dn
  res[which(res$logFC > 0),] -> up
  dn[order(dn$adj.P.Val),] -> dn
  up[order(up$adj.P.Val),] -> up
  dn[1:floor(perc*nrow(dn)),] -> dn
  up[1:floor(perc*nrow(up)),] -> up
  list(up, dn) -> res
  names(res) <- c("up", "dn")
  return(res)
}

createGRP <- function(res, perc=.05, sig.name=""){
    # creates UP/DN file based on top x% of up/downregulated genes; genes ordered
    # by p-value
  subset(res, perc) -> resl
  resl$up -> up
  resl$dn -> dn
  fileConn<-file(paste(sig.name, "_up.grp", sep=""))
  writeLines(c(paste("#", sig.name, "_UP", sep=""), 
               rownames(up)), fileConn)
  close(fileConn)
    
  fileConn<-file(paste(sig.name, "_dn.grp", sep=""))
  writeLines(c(paste("#", sig.name, "_DN", sep=""), 
               rownames(dn)), fileConn)
  close(fileConn)
}

createGRPFromList <- function(genes, tag=c("up", "down"), sig.name=""){
  # creates UP/DN file based on top x% of up/downregulated genes; genes ordered
  # by p-value
  if (tag == "up"){
    fileConn<-file(paste(sig.name, "_up.grp", sep=""))
    writeLines(c(paste("#", sig.name, "_UP", sep=""), 
                 genes), fileConn)
    close(fileConn)
  }
  
  if (tag == "down"){
  fileConn<-file(paste(sig.name, "_dn.grp", sep=""))
    writeLines(c(paste("#", sig.name, "_DN", sep=""), 
                 genes), fileConn)
    close(fileConn)
  }
}

createGMTFromList <- function(geneset, path, minSize){
  # given a list of gene sets, create a .gmt file
  # get longest set in geneset
  sapply(geneset, function(x) length(x)) -> gsl
  geneset[which(gsl >= minSize)] -> geneset
  mat <- matrix("", ncol=max(gsl)+2, nrow=length(geneset))
  mat[,1] <- mat[,2] <- names(geneset)
  for (i in 1:length(geneset)){
    geneset[[i]] -> cg
    mat[i,3:(length(cg)+2)] <- cg
  }
  write.csv(mat, file=path)
}
## Functions
ReadData <- function(f) {
  read.table(gzfile(f), row.names=1, header=TRUE, sep="\t", check.names=FALSE)
} 


BuildCellMetaData <- function(x, is.expr=1){ 
  data.frame(Timepoint=sapply(strsplit(names(x), ".", fixed=TRUE), function(x) (x[1])), 
    Genes=apply(x, 2, function(x) sum(x>is.expr)), 
    Transcipts=colSums(x))
}


FilterGenes <- function(x, n.cells=5, is.expr=1) {
  x[rowSums(x>=is.expr) >= n.cells,]
}


NormalizeData <- function(x, scale.factor=10000) { 
  log((sweep(x, 2, apply(x, 2, function(x) sum(x)), "/")*scale.factor)+1)
}


ScaleData <- function(x, do.scale=TRUE, do.center=TRUE){
  bin.size = 1000
  max.bin = floor(nrow(x)/bin.size) + 1
  for(i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
    my.inds <- my.inds[my.inds <= nrow(x)]
    new.data <- t(scale(t(as.matrix(x[my.inds, ])), center=do.center, scale=do.scale))
    new.data[new.data>10] <- 10
    x[my.inds, ] <- new.data
  }
  x
}


RunTSNE <- function(x, perp=30, id=50, iter=2000, theta=0, dims=2){
  require(Rtsne)
  Rtsne(x, perplexity=perp, theta=theta, max_iter=iter, initial_dims=id, dims=2, pca=TRUE, check_duplicates=FALSE, verbose=TRUE) 
}


AinB <- function(a,b) {
  a[a%in%b]
}


LogExpMean <- function(x) {
  log(mean(exp(x)-1)+1)
}


LogVarDivMean <- function(x) {
  log(var(exp(x)-1)/mean(exp(x)-1))
} 


Pct.Exp <- function(x){
  round(apply(x, 1, function(x)return(length(x[x>1])/length(x))),3)
}


MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2>max] <- max
  data2[data2<min] <- min
  data2
}


DiffLRT <- function(x, y, df=1) {
  lrtX <- BimodLikData(x)
  lrtY <- BimodLikData(y)
  lrtZ <- BimodLikData(c(x,y))
  lrt.diff <- 2*(lrtX+lrtY-lrtZ)
  pchisq(lrt.diff, df, lower.tail=F)
}


BimodLikData <- function(x,xmin=0) {
  x1 <- x[x<=xmin]
  x2 <- x[x>xmin]
  xal <- MinMax(length(x2)/length(x), min=1e-5, max=(1-1e-5))
  likA <- length(x1) * log(1-xal)
  mysd <- sd(x2)
  if(length(x2)<2) { mysd <- 1 }
  likB <- length(x2) * log(xal) + sum(dnorm(x2, mean(x2), mysd, log=TRUE))
  likA+likB
}


NegBinomDiffExpTest <- function(data1, data2, mygenes, slots=1) {
  require(parallel)
  p.val <- unlist(mclapply(mygenes, function(x) DiffLRT(as.numeric(data1[x,]), as.numeric(data2[x,])), mc.cores=slots))
  p.val[is.na(p.val)] <- 1
  avg.diff <- unlist(mclapply(mygenes, function(x) LogExpMean(as.numeric(data1[x,])) - LogExpMean(as.numeric(data2[x,])), mc.cores=slots))
  toRet <- data.frame(row.names=mygenes, cbind(p.val, avg.diff))
  toRet[order(toRet$p.val),]
}


DiffExp <- function(x, data, out, cores=1){
  one.v.all.results <- list()
  for(j in unique(subset(x, Cluster!="Noise" & Cluster!="Small_Cluster")$Cluster)) {
    message(paste0("Processing cluster: ", j))  
    cells.in.cluster = data[,x$Cluster==j]
    cells.outside.cluster = data[,x$Cluster!=j]
    
    # do gene selection (based on percent expressed)
    min.pct = 0.25 ; min.diff.pct = -Inf ; thresh.use = 0.25
    data.cells.in.cluster = round(apply(cells.in.cluster, 1, function(x)return(length(x[x>1])/length(x))),3)
    data.cells.outside.cluster = round(apply(cells.outside.cluster, 1, function(x)return(length(x[x>1])/length(x))),3)
    data.alpha = cbind(data.cells.in.cluster, data.cells.outside.cluster); colnames(data.alpha)=c(paste0("pct.exp.",j), "pct.exp.others")
    alpha.min = apply(data.alpha, 1, max); names(alpha.min) = rownames(data.alpha)
    genes.use=names(which(alpha.min > min.pct))
    alpha.diff = alpha.min - apply(data.alpha, 1, min) 
    genes.use = names(which(alpha.min > min.pct & alpha.diff > min.diff.pct))
    
    # gene selection (based on average difference)
    data.1 = apply(cells.in.cluster, 1, LogExpMean)
    data.2 = apply(cells.outside.cluster, 1, LogExpMean)
    total.diff = (data.1 - data.2)
    genes.diff = names(which(abs(total.diff) > thresh.use))
    genes.use = AinB(genes.use, genes.diff)
    
    results = NegBinomDiffExpTest(cells.in.cluster, cells.outside.cluster, genes.use, slots=cores)
    results = cbind(results, data.alpha[rownames(results),])
    results = results[with(results, order(p.val,-avg.diff)),]
    
    other.data = data.frame(rowMeans(cells.in.cluster), rowMeans(cells.outside.cluster), data.1, data.2)
    colnames(other.data) = c(paste0("mean.exp.",j), "mean.exp.others", paste0("log.exp.mean.",j), "log.exp.mean.others")
    
    to.return = cbind(results, other.data[row.names(results),])
    id <- paste0("C",j)
    one.v.all.results[[id]] <- to.return  
    WriteData(to.return, gzfile(paste0(out, ".negbinom.cluster.v.all.",id,".tsv.gz")))
  }
  one.v.all.results
}
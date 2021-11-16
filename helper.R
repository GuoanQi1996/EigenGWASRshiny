# Source helpers ----
source('EigenGWAS_Friends.R')
FN="arab";
PC=5;
inbred=F;
#Function below controls dataLoading
#RunEigenGWAS(FN, PC, inbred, "gear.jar");

RunPlink <- function(dat, PC, inbred=T, plink2) {
  layout(matrix(1:6, 2, 3))
  #make-grm
  grmCmd=paste(plink2, "--bfile ", dat, "--make-grm-gz --out ", dat)
  system(grmCmd)
  
  gz=gzfile(paste0(dat, ".grm.gz"))
  grm=read.table(gz, as.is = T)
  Ne=-1/mean(grm[grm[,1]!=grm[,2], 4]/2)
  Me=1/var(grm[grm[,1]!=grm[,2], 4]/2)
  print(paste("Ne=", format(Ne, digits = 2), "Me=", format(Me, digits = 2)))
  ## [1] "Ne= 293 Me= 395"
  hist(grm[grm[,1]!=grm[,2],4]/2, main="Pairwise relatedness", xlab="Relatedness score", breaks = 50)
  
  #pca
  pcaCmd=paste(plink2, "--bfile ", dat, "--pca 10 --out ", dat)
  system(pcaCmd)
  barplot(main="Top 10 eigenvalue", read.table(paste0(dat, ".eigenval"), as.is = T)[,1]/2, border = F)
  abline(h=1, col="red", lty=2, lwd=3)
  
  pc=read.table(paste0(dat, ".eigenvec"), as.is = T)
  plot(pc[,3], pc[,4], xlab="Eigenvector 1", ylab="Eigenvector 2", bty="n", main="Eigenspace", bty="n", col=ifelse(pc[,3]>0, "red", "blue"), pch=16, cex=0.5)
  
  #EigenGWAS
  for(i in 1:PC) {
    outRoot=paste0(dat, ".", i)
    liCmd=paste0(plink2, " --linear --bfile ", dat, " --pheno ", dat, ".eigenvec --out ", outRoot)
    system(liCmd)
  }

  # #plot
  # EigenRes=read.table(paste0(dat, ".assoc.linear"), as.is = T, header = T)
  # EigenRes$Praw=EigenRes$P
  # gc=qchisq(median(EigenRes$P), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
  # print(paste("GC = ", format(gc, digits = 4)))
  # ## [1] "GC =  9.047"
  # EigenRes$P=pchisq(qchisq(EigenRes$Praw, 1, lower.tail = F)/gc, 1, lower.tail = F)
  # manhattan(EigenRes, title="EigenGWAS 1", pch=16, cex=0.3, bty='n')
  # 
  # #QQplot
  # chiseq=rchisq(nrow(EigenRes), 1)
  # qqplot(chiseq, qchisq(EigenRes$Praw, 1, lower.tail = F), xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ", chi[1]^2)), bty="n", col="grey", pch=16, cex=0.5)
  # points(sort(chiseq), sort(qchisq(EigenRes$P, 1, lower.tail = F)), col="black", pch=16, cex=0.5)
  # legend("topleft", legend = c("Raw", "GC correction"), pch=16, cex=0.5, col=c("grey", "black"), bty='n')
  # abline(a=0, b=1, col="red", lty=2)
  
  
}
#FileSave functions

# Prepare the compact summary of p-values for QQ-plot
qqPlotCache = function(pvalues, ntests = NULL, ismlog10 = FALSE){
  if(is.null(ntests))
    ntests = length(pvalues);
  
  if(ismlog10) {
    ypvs = pvalues;
  } else {
    ypvs = -log10(pvalues);
  }
  xpvs = -log10(seq_along(ypvs) / (ntests+1));
  
  if(is.unsorted(-ypvs))
    ypvs = sort.int(ypvs, decreasing = TRUE);
  
  if(length(ypvs)*2 > ntests) {
    lambda =
      qchisq(p = 0.1^ypvs[ntests/2], df = 1, lower.tail = FALSE) / 
      qchisq(p = 0.5, df = 1, lower.tail = FALSE);
  } else {
    lambda = NULL;
  }
  
  if(length(ypvs) > 1000) {
    # need to filter a bit, make the plotting faster
    levels = as.integer( (xpvs - xpvs[1])/(tail(xpvs,1) - xpvs[1]) * 2000);
    keep = c(TRUE, diff(levels)!=0);
    levels = as.integer( (ypvs - ypvs[1])/(tail(ypvs,1) - ypvs[1]) * 2000);
    keep = keep | c(TRUE, diff(levels)!=0);
    keep = which(keep);
    ypvs = ypvs[keep];
    xpvs = xpvs[keep];
  } else {
    keep = seq_along(ypvs)
  }
  
  qq = list(
    xpvs = xpvs, # p-values expected under null
    ypvs = ypvs, # observer p-values
    keep = keep, # indices of the preserved p-values
    ntests = ntests, # Number of tests
    lambda = lambda  # Estimate of inflation factor lambda
  );
  class(qq) = "qqPlotInfo";
  return(qq);
}

# Create QQ-plot from p-values or prepared summary
qqPlotQ = function(
  x, 
  ntests = NULL, 
  ismlog10 = FALSE, 
  ci.level = 0.05, 
  ylim = NULL, 
  newplot = TRUE, 
  col = "grey", 
  cex = 0.5, 
  yaxmax = NULL, 
  lwd = 3, 
  col.band = "red",
  makelegend = TRUE,
  title = "",
  xlab = expression(
    paste("\u2013", " log"[10]*"(", italic("P"), ") expected")),
  ylab = expression(
    paste("\u2013", " log"[10]*"(", italic("P"), "), observed"))
){
  
  # Get compact summary of p-values for QQ-plot
  if( methods::is(x, "qqPlotInfo") ){
    qq = x;
  } else {
    qq = qqPlotCache(pvalues = x, ntests = ntests, ismlog10 = ismlog10);
  }
  
  # Axis ranges
  mx = head(qq$xpvs,1) * 1.05;
  if( is.null(ylim) ) {
    my = max(mx, head(qq$ypvs,1) * 1.05) ;
    ylim = c(0, my);
  } else {
    my = ylim[2];
  }
  if(is.null(yaxmax))
    yaxmax = floor(my);
  
  if(newplot){
    plot(
      x = NA, 
      y = NA, 
      ylim = ylim, 
      xlim = c(0, mx), 
      xaxs = "i", 
      yaxs = "i", 
      xlab = xlab,
      ylab = ylab,
      axes = FALSE);
    axistep = floor(mx/5);
    if(axistep<1){
      axistep = 1
    }
    axis(1, seq(0, mx, axistep), lwd = lwd);
    axistep = floor(yaxmax/5);
    axis(2, seq(0, yaxmax, axistep), lwd = lwd);
    abline(a = 0, b = 1, col = col.band, lwd = lwd, lty =2);
    points(qq$xpvs, qq$ypvs, col = col, cex = cex, pch = 19);
    
    if( !is.null(ci.level) ){
      if( (ci.level>0)&(ci.level<1) ){
        quantiles = qbeta(
          p = rep(c(ci.level/2,1-ci.level/2), each=length(qq$xpvs)), 
          shape1 = qq$keep, 
          shape2 = qq$ntests - qq$keep + 1);
        quantiles = matrix(quantiles, ncol=2);
        
        lines( qq$xpvs, -log10(quantiles[,1]), col = col.band, lwd = lwd, lty =2);
        lines( qq$xpvs, -log10(quantiles[,2]), col = col.band, lwd = lwd, lty =2);
      }
    }
    if(makelegend){
      if( !is.null(ci.level) ){
        legend(
          "topleft", 
          legend = c(
            expression(paste(italic("Raw P"), " value")),
            expression(paste(lambda[GC],italic(" corrected P")," value")),
            sprintf("%.0f%% CI",100-ci.level*100)),
          lwd = c(0, 0, lwd), 
          pch = c(19, 19, NA_integer_), 
          lty = c(0, 0, 2), 
          col = c(col, "darkblue", col.band),
          box.col = "transparent",
          bg = "transparent");
      } else {
        legend(
          "topleft", 
          legend = c(
            expression(paste(italic("Raw P"), " value")),
            expression(paste(lambda[GC],italic(" corrected P")," value"))
          ),
          lwd = c(0, 0), 
          pch = c(19, 19), 
          lty = c(0, 0), 
          col = c(col, "darkblue"),
          box.col = "transparent",
          bg = "transparent");
      }
    }
    title(title);
  } else {
    points(qq$xpvs, qq$ypvs, col = "darkblue", cex = cex, pch = 19);
  }
  
  #if( !is.null(qq$lambda) ){
  #    lastr = sprintf("%.3f", qq$lambda);
  #    legend("bottom", legend = bquote(lambda == .(lastr)), bty = "n")
  #}
  return(invisible(qq));
}

# Prepare the compact summary of p-values for Manhattan plot
manhattanCache = function(
  pvalues,
  chr,
  pos,
  ismlog10 = FALSE,
  chrmargins = 1){
  # chr = locs[,1]; pos = locs[,2]; pvalues = mwas[,3]
  # z = getMWASandLocations(param)
  # chr = z$chr; pos = z$start; pvalues = z$`p-value`; 
  # ismlog10 = FALSE; chrmargins = 0;
  
  stopifnot( length(pvalues) == length(chr) );
  stopifnot( length(pvalues) == length(pos) );
  
  chr = factor(chr)
  
  # max of each chromosome
  poslist = split(pos, chr, drop = FALSE);
  poslist[vapply(poslist, length, 0)==0] = list(0);
  chrmax = vapply(poslist, max, 0) + 0;# + chrmargins;
  
  # chromosome starts on the plot
  names(chrmax) = NULL;
  offsets = c(0, cumsum(chrmax)) + chrmargins;
  names(offsets)[seq_along(poslist)] = levels(chr);
  
  # within plot coordinates
  x0 = offsets[unclass(chr)] + pos;
  if(ismlog10) {
    y0 = pvalues;
  } else {
    y0 = -log10(pvalues);
  }
  
  # Prune the data
  yfac = as.integer(y0*100)+1L;
  yorder = sort.list(yfac);
  levels(yfac) = as.character(seq_len(max(yfac)));
  class(yfac) = "factor";
  
  ygroup = split(seq_along(yfac), yfac);
  for( i in seq_along(ygroup)){ # i=1
    if( length(ygroup[[i]]) > 300 ){
      ygroup[[i]] = sample(ygroup[[i]], size = 300, replace = FALSE);
    }
  }
  # sum(vapply(ygroup, length, 0))
  keep = unlist(ygroup, use.names = FALSE);
  
  # Color code
  colindex = unclass(chr);
  
  # Chromosome names
  chrnames = gsub("chr", "", levels(chr));
  
  # Return minimum object;
  man = list(
    x = x0[keep],
    y = y0[keep],
    colindex = colindex[keep],
    offsets = offsets,
    chrnames = chrnames,
    chrmargins = chrmargins
  );
  class(man) = "manPlotInfo";
  return(man);
}

# Create Manhattan plot from prepared summary
manhattan = function(
  man,
  ylim = NULL,
  colorSet = c("steelblue4", "#2C82D1", "#4CB2D1"),
  yaxmax = NULL,
  significant = 1e-5,
  lwd = 3,
  cex = 1,
  title = ""){
  
  if( !methods::is(man, "manPlotInfo") )
    stop("The \"man\" parameter is not produced by manPlotPrepare().");
  
  # Axis ranges
  if(is.null(ylim)) {
    my = max(man$y) * 1.05;
    ylim = c(0,my);
  } else {
    my = ylim[2];
  }
  if(is.null(yaxmax))
    yaxmax = floor(my);
  
  # Plot frame
  plot(
    x = NA,
    y = NA, 
    xlim = c(0, tail(man$offsets,1)), 
    ylim = ylim, 
    xaxs = "i", 
    yaxs = "i",
    xlab = "Chromosome", 
    ylab = expression(
      paste("\u2013", " log"[10]*"(", italic("p"), ")")),
    axes = FALSE);
  axis(
    side = 1,
    at = man$offsets,
    labels = rep("", length(man$offsets)),
    lwd = lwd);
  axis(
    side = 1,
    at = (man$offsets[-1] + man$offsets[-length(man$offsets)])/2,
    labels = man$chrnames,
    tick = FALSE,
    lwd = lwd);
  axis(
    side = 2,
    at = seq(0, yaxmax, floor(yaxmax/5)),
    lwd = lwd);
  
  # plot points in palette color
  oldPal = palette(colorSet);
  points(
    x = man$x,
    y = man$y,
    pch = 16,
    col = ((man$colindex-1L) %% length(colorSet)) + 1L,
    cex = cex);
  abline(h=-log10(significant),lty=2,col='red',lwd=1)
  palette(oldPal);
  title(title);
}


#GRM
GRM_SaveAsPdf <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsJpeg <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsBmp <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsPng <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

#Eigenvalue
Eigenvalue_SaveAsPdf <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsJpeg <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsBmp <- function(path,arg1,arg2,arg3){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsPng <- function(path,arg1,arg2,arg3){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

#miamiPlot
miamiPlot_SaveAsPdf <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsJpeg <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsBmp <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsPng <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

#EigenGWASPlot
EigenGWASPlot_SaveAsPdf <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsJpeg <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsBmp <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsPng <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

# R script to read the GRM binary file
# GCTA
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
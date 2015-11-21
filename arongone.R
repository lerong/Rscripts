wdpath = function(){
  HOME_drive = ifelse(.Platform$OS.type=="windows","//mymdafiles/usersdqs3/lli11","/home/lli11") 
  HN_drive = ifelse(.Platform$OS.type=="windows","//q1prphn/projects/lli11","/projects/lli11")
  Bioinfo2 = ifelse(.Platform$OS.type=="windows", "//d1prpccifs/bcb/bioinfo/bioinfo2", "/data/bioinfo2")
  assign("HOME_drive",HOME_drive,envir = .GlobalEnv)
  assign("HN_drive",HN_drive,envir = .GlobalEnv)
  assign("Bioinfo2",Bioinfo2,envir = .GlobalEnv)
}

## function to calculate fold change
o2fc = function(expData,gf){
  expData = as.matrix(expData)
  gf = as.factor(gf)
  data_m = t(apply(expData, 1, function(x) {tapply(x, gf, mean)}))
  del = data_m[,levels(gf)[1]] - data_m[,levels(gf)[2]]
  fc = matrix(sign(del)*2^abs(del),ncol=1)
  colnames(fc) = paste(levels(gf)[1],levels(gf)[2],"FC",sep=".")
  colnames(data_m) = paste(colnames(data_m),"mean",sep="-")
  return(cbind(data_m,fc))
}
scatter_plot0 = function(x,y,method= method,...){
  x = as.numeric(x)
  y = as.numeric(y)
  cor_rlt = cor.test(x, y, method= method)
  est = round(cor_rlt$estimate,4)
  pvalue  = format.pval(cor_rlt$p.value,4)
  
  lm_xy = summary(lm(y~x))
  p_lm = format.pval(lm_xy$coefficients[2,4],4)
  r2_lm = round(lm_xy$adj.r.squared,4)
  
  mymain = paste(method, "correlation: rho =",est,", pvalue =",pvalue,
                 "\nLinear Model: ", "Adj-Rsquare = ",r2_lm,", pvalue =", p_lm)
  op = par(mar=c(5, 4, 4, 2) + 0.5,mfrow=c(1,1))
  plot(x,y,pch=16,col="gray",panel.first = grid(8, 8),...)
  title(mymain,font.main = 10,col.main= "black",cex.main=1)
  abline(lm(y~x),col="red",lty=3,lwd=2)
  par(op)
}

scatter_plot = function(x,y,...){
  # x and y are two numerical vector
  cor_pearson  = cor.test(x, y, method="pearson") 
  cor_spearman = cor.test(x, y, method="spearman",exact = FALSE)
  lm_xy = summary(lm(y~x))
  e_pearson = round(cor_pearson$estimate,3)
  e_spearman = round(cor_spearman$estimate,3)
  p_pearson  = format.pval(cor_pearson$p.value,3)
  p_spearman = format.pval(cor_spearman$p.value,3)
  p_lm = format.pval(lm_xy$coefficients[2,4],3)
  r2_lm = round(lm_xy$adj.r.squared,3)
  mymain = paste("Pearson Correlation: ", "r = ",e_pearson,", p =",p_pearson,
                 "\nSpearman Correlation: ", " rho = ",e_spearman,", p =",p_spearman,
                 "\nLinear Model: ", "Adj-Rsquare = ",r2_lm,", p =", p_lm )
  op = par(mar=c(5, 5, 6, 5))
  plot(x,y,pch=16,panel.first = grid(8, 8),...)
  title(mymain,font.main = 10,col.main= "black",cex.main=1)
  abline(lm(y~x),col="red",lty=3,lwd=2)
  par(op)
}

stat.bar = function(stat_rlt,xlab,ylab,title,pcut = 0.05,cex.y = 5){
  stat_rlt$Significance = pval2factor(stat_rlt$pvalue)
  stat_rlt$labs = rownames(stat_rlt)
  sid = topvalue(stat_rlt$pvalue,ind.cut = pcut)
  stat_rlt = stat_rlt[sid,]
  p = ggplot(stat_rlt,aes(x = reorder(labs,correlation),y=correlation,fill=Significance),environment = environment()) 
  p+ geom_bar(stat ="identity") + coord_flip() + xlab(xlab) + ylab(ylab) + labs(title=title)+
    #scale_x_discrete(breaks = reorder(stat_rlt$labs,stat_rlt$correlation)[sid], labels = stat_rlt$labs[sid])+
    theme_bw() + theme(axis.text.y = element_text(colour="grey20",size=cex.y,angle=0,hjust=1,vjust=0))
}

stat_bar = function(stat_rlt,xlab,ylab,title,cex.y = 5){
  stat_rlt$Significance = pval2factor(stat_rlt$pvalue)
  stat_rlt$labs = rownames(stat_rlt)
  p = ggplot(stat_rlt,aes(x = reorder(labs,correlation),y=correlation,fill=Significance),environment = environment()) 
  p+ geom_bar(stat ="identity") + coord_flip() + xlab(xlab) + ylab(ylab) + labs(title=title)+
    #scale_x_discrete(breaks = reorder(stat_rlt$labs,stat_rlt$correlation)[sid], labels = stat_rlt$labs[sid])+
    theme_bw() + theme(axis.text.y = element_text(colour="grey20",size=cex.y,angle=0,hjust=1,vjust=0))
}


box_plot = function(x,class,srt=45,...){
  aov_out = aov(x ~ class)
  aov_summary = summary(aov_out)
  mean_grp = tapply(x,class,mean)
  
  n_grp = length(unique(class))
  test_name = ifelse(n_grp>2,"ANOVA","t-test")
  test_p = format.pval(aov_summary[[1]]$'Pr(>F)'[1],4)
  
  mymain = paste(test_name,"p-value: ", test_p)
  
  op =  par(mar = c(10, 4, 5, 2) + 0.2)
  boxplot(x~class,main=mymain,xaxt="n",...)
  points(mean_grp,col="blue",pch=3)
  stripchart(x~class,vertical=TRUE,pch=16,cex=1,method='jitter', add=TRUE,col=2) 
  
  axis(1,1:nlevels(class),labels = F)
  text(x = 1:nlevels(class), par("usr")[3] - 0.5, labels = levels(class), srt = srt, adj = 1, xpd = TRUE)
  par(op)
  
  rlt = list(TukeyHSD =  TukeyHSD(aov_out))
  rlt
}




gg_qqplot <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                      labels = names(x),lsize = 3){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z,...)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point(shape=1,size=4) +  
    geom_abline(intercept = coef[1], slope = coef[2],colour="red",size=0.58,ltype="dashed") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) + 
    labs(x = "Theoretical Quantiles",y="Sample Quantiles",title=paste("Q-Q ",distribution,sep=""))
  #theme_classic()
  if(!is.null(labels)) p <- p + geom_text(aes(label = label,family="mono"), size=lsize,colour="blue")
  print(p)
  coef
}

#-------------------------------------------------------------------------------------------------------------------------
# Calculate the correlation between y and a vector of random variable column by column
# input: y--vector; xmtx--matrix
# output: matrix of correlation estimate and p-values
# updated: 06-15-2015
#-------------------------------------------------------------------------------------------------------------------------
vec_cor = function(y,xmtx,method="spearman",...){
  rlt = matrix(0, ncol(xmtx), 2, byrow=TRUE)
  rownames(rlt) = colnames(xmtx)
  #colnames(rlt) = paste(method,c("r","pvalue"),sep = "-")
  colnames(rlt) = c("correlation","pvalue")
  for(i in 1:nrow(rlt)){
    xi = as.numeric(as.vector(xmtx[,i]))  
    ######################  spearman correlaton
    cortestresult = cor.test(xi, y, method=method,... )
    rlt[i,1] = cortestresult$estimate
    rlt[i,2] = cortestresult$p.value
  }
  data.frame(rlt)
}

mtx_cor = function(mtxA,mtxB,method="pearson",...){
  rlt = list()
  cor_rlt = pval_rlt =  matrix(NA, nrow(mtxA),nrow(mtxB))
  for (i in 1:nrow(mtxA)){
    vec_rlt = vec_cor(matrix(mtxA[i,],ncol=1),t(mtxB),method=method,...)
    cor_rlt[i,] = vec_rlt[,1]
    pval_rlt[i,] = vec_rlt[,2]
  }
  rlt$corMtx= cor_rlt
  rlt$pMtx = pval_rlt
  return(rlt)
}



vec_ttest = function(xmtx,class,...){
  class = as.factor(class)
  rlt = vector()
  for (i in 1:ncol(xmtx)){
    tt = t.test(xmtx[,i]~class,...)
    exname = c("statistic","conf.int","p.value","estimate")
    ttout = unlist(tt[exname])
    # fold change
    del = unname(diff(tt$estimate))
    fc = sign(del)*2^abs(del)
    names(fc) = paste(levels(class)[2],levels(class)[1],"FC",sep="-")
    ttout = c(ttout,fc)
    rlt = rbind(rlt,ttout)
  }
  rownames(rlt) = colnames(xmtx)
  data.frame(rlt)
}


plot_ttest = function(x,class,...){
  class = as.factor(class)
  tt = t.test(x ~ class)
  p = format.pval(tt$p.value,4)
  mymain = paste("t-test p-value: ", p)
  boxplot(x ~ class,main=mymain,...)
  points(tt$estimate,col="blue",pch=3)
  stripchart(x ~ class,vertical=TRUE,pch=18,cex=1,method='jitter', add=TRUE,col=2) 
}


vec_annova = function(xmtx,class){
  class = as.factor(class)
  rlt = vector()
  an =  aov(xmtx[,3] ~ class)
  aov_summary = summary(an)
  round(aov_summary[[1]]$'Pr(>F)'[1],4)
}

#-------------------------------------------------------------------------------------------------------------------------
# Generate significant tables 
#-------------------------------------------------------------------------------------------------------------------------
sig_tab = function(pvec, sig_vec=c(0.01, 0.05, 0.1, 0.2), method="FDR",...){
  require(ClassComparison)
  #require(ggplot2)
  tmpBum = Bum(pvec)
  selectors = sapply(sig_vec, function(alpha) selectSignificant(tmpBum, alpha, by=method))
  counts = sapply(sig_vec, function(alpha) countSignificant(tmpBum, alpha, by=method))
  cuts = sapply(sig_vec, function(alpha) cutoffSignificant(tmpBum, alpha, by=method))
  temp = data.frame(sig_vec, counts, cuts)
  colnames(temp) = c(paste(method, "s",sep=''), "Num_Sig", "Corresponding p-value")
  rlt = list()
  rlt$sigTab = temp
  rlt$selectors = selectors
  rlt
}

#-------------------------------------------------------------------------------------------------------------------------
# Based on  Li Shen 's function
#' fake.data <- c(runif(700), rbeta(300, 0.3, 1))
#' gg_Bum_Tab(fake.data, main="", heights=c(0.75, 0.25))
#-------------------------------------------------------------------------------------------------------------------------
# run bum and return sid
bum_table = function(pvals, method='FDR', 
                     fdrs=c(0.01, 0.05, 0.10, 0.20), pin = c(0.001,0.01,0.05,0.1),
                     res=100, heights=c(3,1), main="Bum Plot",plot=TRUE) {
  require(gridExtra)
  require(ClassComparison)
  #par(mar=c(5, 4, 2, 2) + 0.1,mfrow=c(1,1))
  a = Bum(pvals)
  counts = sapply(fdrs, function(alpha) countSignificant(a, alpha, by=method))
  cuts = sapply(fdrs, function(alpha) cutoffSignificant(a, alpha, by=method))
  mat = data.frame(fdrs, counts, format.pval(cuts,4))
  colnames(mat) = c(paste(method, "s",sep=''), "Num_Sig", "Corresponding Pvalues")
  
  # bum plot
  if(plot == TRUE){
  fdat = data.frame(pvals=pvals)
  xvals = (0:res)/res
  fit = a@lhat + (1 - a@lhat) * dbeta(xvals, a@ahat, 1)
  betaf = data.frame(xvals=xvals, fit=fit)[-1,]
  den.plot = ggplot2::ggplot(fdat, aes(x=pvals)) + 
    geom_histogram(aes(y = ..density..), colour="black", fill = "white", binwidth = 0.01) +
    geom_line(data=betaf, aes(xvals,fit), colour="darkgreen") +
    geom_hline(yintercept = a@pihat, colour="blue") +
    labs(x="P Values", y="Density", title=main) +
    theme_classic()
  
  # output table
  ptab = data.frame(Pvalues=pin, Num_In_Set=sapply(pin,function(x)length(which(pvals<=x))))
  ptab.plot = tableGrob(ptab)
  ftab.plot = tableGrob(mat)
  grid.arrange(den.plot, grid.arrange(ftab.plot,ptab.plot,nrow=1), ncol=1, heights=heights)
  }
  
  # select gene id
  gid.adjust = sapply(fdrs,function(alpha)selectSignificant(a,alpha,by=method))
  colnames(gid.adjust) = paste("gidAdj",fdrs,sep="")
  gid.p = sapply(pin, function(x)pvals <= x)
  colnames(gid.p) = paste("gidP",pin,sep="")
  gid = data.frame(gid.adjust,gid.p)
  return(gid)
}


pval2factor = function(pvals,breaks=c(0,0.001,0.01,0.05,1)){
  cut(pvals,breaks=breaks,labels=c("p<0.001","p<0.01","p<0.05","p>0.05"),include.lowest=TRUE)
}

topvalue = function(pvals,ind.cut=0.05,fdr.cut=NA){
  gid = which(pvals<=ind.cut)
  if(!is.na(fdr.cut)){
    tmpBum = Bum(pvals)
    selectors = selectSignificant(tmpBum, fdr.cut, by="FDR")
    gid = which(selectors)
  }
  return(gid)
}

adjust.p = function(pvals,cut = 1, method = "none"){
  p.adj = p.adjust(pvals,method = method)
  gid = which(p.adj <= cut)
  return(gid)
}

#-------------------------------------------------------------------------------------------------------------------------
# calculate EMT score from Pan Tong
#--------------------------------------------------------------------------------------------------------------------------
getEMTscore = function(mat=matSel, status=status, scale=FALSE) {
  if(scale) mat <- t(scale(t(mat))) 
  if(length(status)!=nrow(mat)) stop("Gene length does not equal to E/M status!\n")
  Mmean = apply(mat[status=='M', ], 2, mean, na.rm=T) # NA removed
  Emean = apply(mat[status=='E', ], 2, mean, na.rm=T)
  EMTscore = Mmean-Emean
  names(EMTscore) = colnames(mat)
  #browser()
  #heatmap.2(mat, col=greenred(30), trace='none', symbreaks=T)
  EMTscore
}

library(RJSONIO) # function to extract mirna (mature) target genes from targetHub using any prediction algorithm
targetHub.byMethods <- function(mir.name, method) {
  tarhub.matmirna <- 'http://app1.bioinformatics.mdanderson.org/tarhub/_design/basic/_view/by_matureMIRmethod'
  method <- gsub("\\+", "%2B", method)
  data.link <- gsub("\\\"", "%22", paste(tarhub.matmirna, '?key=["', mir.name, '","', method ,'"]&include_docs=true', sep=''))
  json.data <- paste(readLines(data.link), collapse='') #get data from targetHub
  target_data <- RJSONIO::fromJSON(json.data) #convert json formatted data to list
  target_Info <- NA
  if(is.list(target_data) & (length(target_data$rows) > 0)) {
    target_data <- target_data$rows
    target_Info <- matrix(nrow = length(target_data), ncol = 3)
    colnames(target_Info) <- c("Gene_Symbol", "Gene_EntrezID", "miRNA")
    for(i in 1:length(target_data)) {
      target_Info[i,1] = target_data[[i]]$doc$Gene_Symbol
      target_Info[i,2] = target_data[[i]]$doc$Gene_EntrezID
      target_Info[i,3] = target_data[[i]]$doc$miRNA
    }    
  }
  target_Info
}


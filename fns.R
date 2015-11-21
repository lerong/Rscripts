# cut vector into equal parts in order
chunk = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# drop NA in vectors
dropNA = function(x)x[!is.na(x)]

# cut p values into sig categories
p2f = function(pvals,breaks=c(0,0.001,0.01,0.05,1)){
  cut(pvals,breaks=breaks,labels=c("p<0.001","p<0.01","p<0.05","p>0.05"),include.lowest = TRUE)
}

# pdf output
pdf_out = function(file,family="mono", pointsize=16, width=16,height=10,...){
  pdf(file = file,family=family,pointsize = pointsize,width=width,height =height,...)
}


ccorr = function(y,xmtx,method="spearman",...){
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
  rlt = data.frame(rlt)
  rlt$sig = p2f(rlt$pvalue)
  rlt
}


row_ttest = function(data,class,pair=NULL,na.rm=TRUE,var.equal=TRUE,base=2){
  require(ClassComparison)
  if(is.null(pair)){
    if(var.equal) res = MultiTtest(data=data,classes = class,na.rm = na.rm) else
      res = MultiTtestUnequal(data=data,classes = class)
  }else{# paired data
    res = MultiTtestPaired(data=data,classes = class,pairing = pair)
  }
  # calculate group mean and fc
  v = class==levels(class)[1]
  am = matrixMean(data[, v],na.rm = na.rm)
  bm = matrixMean(data[, !v],na.rm = na.rm)
  
  if(!is.null(base)){
    del = bm - am
    fc = sign(del)*base^abs(del)
  }else{# if original scale and positive, use standard foldchange
    if(all(data>0,na.rm = na.rm))fc = foldchange(bm,am) 
    else {
      warning("Please check data for FC")
      fc = NA
    }
  }
  rlt = data.frame(t.stats = res@t.statistics,p.values=res@p.values ,amean=am,bmean=bm,
                   fc = round(fc,2),Significance=p2f(res@p.values))
  rlt
}



row_wilcoxon = function(data,class){
  
}

row_lm = function(formula,class,data){

  
}



# correlation

row_cortest = function(xmtx,yvec,method="pearson",use = "pairwise"){
  require(psych)
  xmtx = t(xmtx)
  yvec = as.numeric(matrix(yvec,ncol = 1))
  res = corr.test(x=xmtx,y=yvec,use=use,method=method,adjust="none",alpha=alpha,ci=FALSE)
  rlt = data.frame(correlation = res$r,pvalue=res$p,significance = p2f(res$p))
  rlt
}

mtxRow_cortest = function(xmtx, ymtx, method = "pearson", use = "pairwise", 
                       adjust = "none", alpha = 0.05, ci = FALSE){
  require(psych)
  x = t(xmtx)
  y = t(ymtx) 
  res = corr.test(x=x,y=y,use=use,method=method,adjust=adjust,alpha=alpha,ci=ci)
  return(res)
}








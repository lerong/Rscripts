
# 1 run bum and return sid
bum_table = function(pvals, method='FDR', 
                     fdrs=c(0.01, 0.05, 0.10, 0.20), pin = c(0.001,0.01,0.05,0.1),
                     res=100, heights=c(3,1), main="Bum Plot") {
  require(gridExtra)
  require(ClassComparison)
  #par(mar=c(5, 4, 2, 2) + 0.1,mfrow=c(1,1))
  a = Bum(pvals)
  counts = sapply(fdrs, function(alpha) countSignificant(a, alpha, by=method))
  cuts = sapply(fdrs, function(alpha) cutoffSignificant(a, alpha, by=method))
  mat = data.frame(fdrs, counts, format.pval(cuts,4))
  colnames(mat) = c(paste(method, "s",sep=''), "Num_Sig", "Corresponding Pvalues")
  
  # bum plot
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
  # select gene id
  gid.adjust = sapply(fdrs,function(alpha)selectSignificant(a,alpha,by=method))
  colnames(gid.adjust) = sprintf("gidAdj_%02d",fdrs*100)
  gid.p = sapply(pin, function(x)pvals <= x)
  colnames(gid.p) = sprintf("gidP_%03d",pin*1000)
  gid = data.frame(gid.adjust,gid.p)
  return(gid)
}

# 2
rowttests = function(m,fac,na.rm=FALSE,var.equal=TRUE){
  require(ClassComparison)
  v <- fac == levels(fac)[1]
  if (na.rm) {
    an <- apply(m[, v], 1, function(y) sum(!is.na(y)))
    bn <- apply(m[, !v], 1, function(y) sum(!is.na(y)))
  }
  else {
    an <- sum(v)
    bn <- sum(!v)
  }
  am <- matrixMean(m[, v],na.rm = na.rm)
  av <- matrixVar(m[, v], am,na.rm = na.rm)
  bm <- matrixMean(m[, !v],na.rm = na.rm)
  bv <- matrixVar(m[, !v], bm,na.rm = na.rm)
  
  if(var.equal){
    tt = (am - bm)/sqrt(((an - 1) * av + (bn - 1) * bv)/(an + bn - 2))/sqrt(1/an + 1/bn)
    df <- apply(m, 1, function(x) sum(!is.na(x)) - 2)
    p.values <- 2 * pt(-abs(tt), df)
  }else{
    tt <- (am - bm)/sqrt(av/an + bv/bn)
    u <- bv/av
    #df <- trunc((1/an + u/bn)^2/(1/(an^2 * (an - 1)) + u^2/(bn^2 * (bn - 1))))
    df <- (1/an + u/bn)^2/(1/(an^2 * (an - 1)) + u^2/(bn^2 * (bn - 1)))
    p.values <- 2 * pt(-abs(tt), df)
  }
  
  rlt = data.frame(t.stats = tt,p.values = p.values, df = round(df,2), amean = am, bmean = bm )
  rownames(rlt) = rownames(m)  
  rlt
}

# 
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

## pdf_output

pdf_output = function(file,...){
  pdf(file = file,family="mono", pointsize=16, width=16, height=10,...)
}
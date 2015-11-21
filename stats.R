
rowttests = function(m,fac,na.rm=FALSE,var.equal=TRUE,fc.title=FALSE){
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
  
  # for expression data fold change
  del = bm - am
  fc = sign(del)*2^abs(del)
  rlt = data.frame(t.stats = tt,p.values = p.values, df = round(df,2), amean = am, bmean = bm, fc ,Significance = p2f(p.values) )
  if(fc.title){
    colnames(rlt) = c("t.stats","p.values","df",
                      paste("mean","(",levels(fac)[1],")",sep=""),
                      paste("mean","(",levels(fac)[2],")",sep=""),
                      paste("(",levels(fac)[2],"|",levels(fac)[1],")",": fc",sep=""),"sig")
  }
  rownames(rlt) = rownames(m)  
  rlt
}
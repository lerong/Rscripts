#library(dplyr)
library(ggplot2)
library(squash)
library(colorRamps)
library(RColorBrewer)



## qqplot using ggplot2:  from stackoverflow
ggplot2.qqplot <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,labels = names(x)){
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
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) 
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  print(p)
  coef
}


makecmap.cat0 = function (x, colFn = NULL,colvec=NA,...) {
  x = as.factor(x)
  ncat = nlevels(x)
  if (ncat < 2 )
    warning("check the number of levels of x")
  if (is.na(colvec) && length(colvec) != levels(x)){
    stop("number of input colors shoule match the level of x")
  }
  if(!is.null(colFn)){
    myColors = colFn(ncat,...)
  } else{
    myColors = colvec 
  }
  list(x = x, xcolors = myColors[x] ,breaks = levels(x),colors = myColors)
}

makecmap.cat1 = function(mat,colSet ="Set1"){
  rlt = list()
  mat = as.matrix(mat)
  nc = ncol(mat)
  nset = length(colSet)
  if( nset != nc) {
    avail.set = rownames(subset(brewer.pal.info,category=="qual"))
    colSet = c(colSet, sample(setdiff(avail.set,colSet),nc-nset))
  }
  ncol = apply(mat,2,function(x)length(unique(x)))
  for (i in 1:nc){
    x = as.factor(mat[,i])
    nl = nlevels(x)
    myColors = brewer.pal(max(nlevels(x),3),colSet[i])[1:nl]
    rlt[[i]] = list(x = x, xcolors = myColors[x] ,breaks = levels(x),colors = myColors)
  }
  names(rlt) = colnames(mat)
  matcolor = sapply(rlt,function(a)a$xcolors)
  list(colorMtx = matcolor,map = rlt)
}

vlegend = function (map, title = NA, side = 2, stretch = 1, x, y,border=NA, cex.title=1){
  opar <- par(mar = c(0, 0, 0, 0))
  #opar = par(fig=c(0, 1, 0, 1), oma=c(0.1, 0.1, 0.1, 0.1), mar=c(0, 0, 0, 0),xpd = NA)
  #opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
  on.exit(par(opar))
  n <- length(map$breaks)+1
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,from = "inches"))
  dx <- dy * aspect
  labs <- format(map$breaks)
  maxlabwidth <- max(strwidth(labs))
  if (missing(x)) {
    x <- grconvertX(1, from = "nfc") - (2 * dx)
    if (side == 4)
      x <- x - maxlabwidth - dx
  }
  else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if (missing(y))
    y <- par("usr")[3] + dy
  ybord <- y + ((1:(n)) * dy * stretch)
  rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = border)
  #box("figure",lty = 2)
  if (side == 4) {
    xtext <- x + dx
    text(x = x, y = ybord[n] + (1* dy), title, adj = c(0,0), cex=cex.title)
  }
  if (side == 2) {
    xtext <- x
    text(x = x + dx, y = ybord[n] + (1 * dy), title, adj = c(1,0), cex=cex.title)
  }
  text(x = xtext, y = ybord[-1]-diff(ybord)/2, labels = labs, pos = side, cex=cex.title)
  res = list()
  res$position = c(x,y)
  invisible(res)
}

hlegend = function (map, title = NA, side = 1, stretch = 1, x, y,border=NA, cex.title=1) {
  opar <- par(mar = c(0, 0, 0, 0))
  #opar = par(fig=c(0, 1, 0, 1), oma=c(0.1, 0.1, 0.1, 0.1), mar=c(0, 0, 0, 0),xpd = NA)
  on.exit(par(opar))
  n <- length(map$breaks)+1
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,from = "inches"))
  dx <- dy * aspect
  labs <- format(map$breaks)
  labwidth <- strwidth(labs)
  if (missing(x)) {
    x <- grconvertX(0, from = "nfc") + dx + (0.5 * strwidth(format(map$breaks[1])))
  }
  else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if (missing(y)) 
    y <- grconvertY(0, from = "nfc") + (2 * dy)
  xbord <- x + ((0:(n - 1)) * dx * stretch)
  rect(xbord[-1], y, xbord[-n], y + dy, col = map$colors, border = border)
  if (side == 1) {
    ytext <- y
    text(x = x, y = y + (1.5 * dy), title, adj = c(0, 0),cex=cex.title)
  }
  if (side == 3) {
    ytext <- y + dy
    text(x = x, y = y - (0.5 * dy), title, adj = c(0, 1),cex=cex.title)
  }
  text(x = xbord[-n] + diff(xbord)/2, y = ytext, labels = labs, pos = side,cex=cex.title)
  res = list()
  res$position = c(x,y)
  invisible(res)
}


hkeybar = function (map, title = NA, side = 1, stretch = 1, x, y, skip, wh,cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  #opar <- par(xpd = T)
  #opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
  opar <- par(mar = c(0, 0, 0, 0))
  on.exit(par(opar))
  n <- length(map$breaks)
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2, from = "inches"))
  dx <- dy * aspect
  labs <- format(map$breaks)
  labwidth <- strwidth(labs)
  if (missing(x)) {
    x <- grconvertX(0, from = "nfc") + dx + (0.5 * strwidth(format(map$breaks[1])))
  }
  else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if (missing(y)) 
    y <- grconvertY(0, from = "nfc") + (2 * dy)
  xbord <- x + ((0:(n - 1)) * dx * stretch)
  if (missing(wh)) {
    if (missing(skip)) {
      for (i in 1:min(n, 20)) {
        if ((n - 1)%%i == 0) {
          step <- (n - 1)/i
          wh.tmp <- seq(1, n, by = step)
          maxlabwidth <- max(strwidth(format(map$breaks[wh.tmp])))
          if (maxlabwidth + dx < dx * step * stretch) 
            wh <- wh.tmp
        }
      }
    }
    else {
      wh <- seq(1, n, by = skip)
    }
  }
  rect(xbord[-n], y, xbord[-1], y + dy, col = map$colors, 
       border = NA)
  if (side == 1) {
    ytext <- y
    text(x = x, y = y + (1.5 * dy), title, adj = c(0, 0),cex=cex.title)
  }
  if (side == 3) {
    ytext <- y + dy
    text(x = x, y = y - (0.5 * dy), title, adj = c(0, 1),cex=cex.title)
  }
  text(x = xbord[wh], y = ytext, labels = format(map$breaks[wh]), pos = side,cex=cex.title)
  res = list()
  res$position = c(x,y)
  invisible(res)
}

vkeybar = function (map, title = NA, side = 2, stretch = 1, x, y, skip, wh, cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  #opar <- par(xpd = NA)
  #opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
  opar <- par(mar = c(0, 0, 0, 0))
  on.exit(par(opar))
  n <- length(map$breaks)
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2, 
                                                                   from = "inches"))
  dx <- dy * aspect
  if (missing(wh)) {
    if (missing(skip)) {
      for (i in 1:min(n, 20)) {
        if ((n - 1)%%i == 0) {
          step <- (n - 1)/i
          wh.tmp <- seq(1, n, by = step)
          if (strheight("A") * 1.2 < dy * step * stretch) 
            wh <- wh.tmp
        }
      }
    }
    else {
      wh <- seq(1, n, by = skip)
    }
  }
  labs <- format(map$breaks[wh])
  maxlabwidth <- max(strwidth(labs))
  if (missing(x)) {
    x <- grconvertX(1, from = "nfc") - (2 * dx)
    if (side == 4) 
      x <- x - maxlabwidth - dx
  }
  else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if (missing(y)) 
    y <- par("usr")[3] + dy
  ybord <- y + ((0:(n - 1)) * dy * stretch)
  rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
  if (side == 4) {
    xtext <- x + dx
    text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0, 0),cex=cex.title)
  }
  if (side == 2) {
    xtext <- x
    text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1, 0),cex=cex.title)
  }
  text(x = xtext, y = ybord[wh], labels = labs, pos = side,cex=cex.title)
  res = list()
  res$position = c(x,y)
  invisible(res)
}

#-------------------------------------------------------------------------------------------------------------------------
# rfunctions for plots and visulizations
# ---showpanel.R
# ---cmap.cat.R
#--------------------------------------------------------------------------------------------------------------------------
require(ggplot2)
require(gplots)
require(squash)
require(RColorBrewer)
library(colorRamps)

#-------------------------------------------------------------------------------------------------------------------------
# given col to check
# col: vector of colors, usullay outputs of certain color pallete 
# eg: showpanel(greenred(64))
#--------------------------------------------------------------------------------------------------------------------------
showpanel = function(col,...){
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n",...)
}

#
# overlay a (empty) graph over the complete plotting area and adding the legend to that
# http://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend = function(...) {
  opar = par(fig=c(0, 1, 0, 1), oma=c(2, 2, 2, 2), mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}




## prepare the side bar of heatmap
## continouse bar and discrete bar
## continouse bar: colorRamps
## discrete: RColorBrewer
catcmap2 = function(mat,colSet ="Set1"){
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
    
#-------------------------------------------------------------------------------------------------------------------------
# Generate a color map from a factor to a set of colors. Corresponding to makecmap in library(squash)
# automatically generate a mapping that can be used in combination with cmap to represent category data with colors in a consistent way
#--------------------------------------------------------------------------------------------------------------------------
catcmap = function (x, colFn = NULL,colvec=NA,...) {
  x = as.factor(x)
  ncat = nlevels(x)
  if (ncat < 2 )
    warning("check the number of levels of x, should > 1")
  if (!is.na(colvec) && length(colvec) != nlevels(x)){
    stop("number of input colors shoule match the level of x")
  }
  if(!is.null(colFn)){
    myColors = colFn(ncat,...)
  } else{
    myColors = colvec 
  }
  list(x = x, xcolors = myColors[x] ,breaks = levels(x),colors = myColors)
}

hkeybar = function (map, title = NA, side = 1, stretch = 1.4, x, y, skip, wh,cex = 1,cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  #opar <- par(xpd = T)
  opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
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
  text(x = xbord[wh], y = ytext, labels = format(map$breaks[wh]), pos = side,cex=cex)
}

vkeybar = function (map, title = NA, side = 2, stretch = 1.4, x, y, skip, wh, cex = 1,cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  #opar <- par(xpd = NA)
  opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
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
  text(x = xtext, y = ybord[wh], labels = labs, pos = side,cex=cex)
}
vlegend = function (map, title = NA, side = 2, stretch = 1.4, x, y,border=NA, cex=1.0, cex.title=1.0){
  #opar <- par(xpd = NA)
  #opar = par(fig=c(0, 1, 0, 1), oma=c(0.1, 0.1, 0.1, 0.1), mar=c(0, 0, 0, 0),xpd = NA)
  opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
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
  if (side == 4) {
    xtext <- x + dx
    text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0,0), cex=cex.title)
  }
  if (side == 2) {
    xtext <- x
    text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1,0), cex=cex.title)
  }
  text(x = xtext, y = ybord[-1]-diff(ybord)/2, labels = labs, pos = side, cex=cex)
}

hlegend = function (map, title = NA, side = 1, stretch = 1.4, x, y,cex=1, cex.title=1.0, border=NA) {
  opar <- par(fig=c(0, 1, 0, 1),xpd = NA)
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
  text(x = xbord[-n] + diff(xbord)/2, y = ytext, labels = labs, pos = side,cex=cex)
}
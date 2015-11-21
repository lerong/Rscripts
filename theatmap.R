library(gtools)
# to do list
# 1- plot the main image use cimage so can drop the border parts
theatmap = function (x, 
                      # dendrogram control
                      Rowv = TRUE, Colv = (if (symm) "Rowv" else TRUE), 
                      
                      # distance and cluster methods
                      row.dist = "euclidean", col.dist = "euclidean",
                      row.clust = "ward.D2", col.clust = "ward.D2",
                      
                      dendrogram = c("both","row", "column", "none"), 
                      reorderfun = function(d, w) reorder(d, w), 
                      symm = FALSE,
                      
                      # data scaling
                      scale = c("none", "row", "column"), 
                      na.rm = TRUE, 
                      ### added: triming data
                      trim = FALSE,trim.size = 2,
                      
                      # image plot
                      revC = identical(Colv, "Rowv"),
                      add.expr,
                      breaks = pretty, #splitting points for binning x into colors, changed to be breaks functions or vector
                      symbreaks = any(x < 0, na.rm = TRUE) || scale !=  "none", 
                      
                      # main image colors
                      col = redblue,
                      ncolor = 20,
                      # block sepration
                      #colsep, rowsep, 
                      sepcolor = "white", 
                      sepwidth = c(0.05, 0.05), 
                      
                      # cell labeling
                      cellnote, notecex = 1, notecol = "cyan", 
                      na.color = par("bg"), 
                      
                      # level trace: dropped
                      #trace = c("column", "row", "both", "none"), 
                      #tracecol = "cyan", hline = median(breaks), 
                      #vline = median(breaks), linecol = tracecol, 
                      
                      # Row/Column Labeling
                      margins = c(10,10,4,5), 
                      ColSideColors, RowSideColors,th=0.2,
                      cexRow = 0.2 + 1/log10(nr), 
                      cexCol = 0.2 + 1/log10(nc), 
                      labRow = NULL, labCol = NULL, 
                      srtRow = NULL, srtCol = NULL, 
                      adjRow = c(0, NA), adjCol = c(NA, 0), 
                      offsetRow = 0.5, offsetCol = 0.5, 
                      colRow = NULL, colCol = NULL,
                      
                      # color key + density info
                      #density.info = c("histogram",  "density", "none"), denscol = tracecol, 
                      #symkey = any(x <  0, na.rm = TRUE) || symbreaks, 
                      #densadj = 0.25,
                      key = TRUE, keysize = 1, 
                      key.title = NULL,
                      key.stretch = 1,
                     #key.xlab = NULL, key.ylab = NULL, 
                      #key.xtickfun = NULL, key.ytickfun = NULL, 
                      #key.par = list(), 
                      
                      # plot labels 
                      main = NULL, xlab = NULL, ylab = NULL, 
                      # plot layout
                      lmat = NULL, lhei = NULL, lwid = NULL, 
                      
                      extrafun = NULL,
                      
                      ## added options
                      barwidth = 1, # the size of barside
                      gap = 0.1, # the dist between main image and side bar
                      sborder = NA, # the border in sidebar
                      ...) {
  
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  trime = function(xmtx,trim.size){
    xmtx[xmtx > trim.size] = trim.size
    xmtx[xmtx <  -trim.size] = -trim.size
    xmtx
  }
  
  retval <- list()
  scale <- if (symm && missing(scale)) "none" else match.arg(scale)
  
  dendrogram <- match.arg(dendrogram)
  
  # get the colFn
  if (length(col) == 1 && is.character(col))  col <- get(col, mode = "function")
  
  if (!missing(breaks) && (scale != "none")) 
    warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", 
            "Please consider using only one or the other.")
  
  if (is.null(Rowv) || is.na(Rowv)) Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv)) Colv <- FALSE  else if (all(Colv == "Rowv")) Colv <- Rowv
  
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) stop("`x' must be a numeric matrix")
  
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 4)  stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote)) cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  
  if (!inherits(Rowv, "dendrogram")) {
    if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) 
        dendrogram <- "column"
      else dendrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  
  if (!inherits(Colv, "dendrogram")) {
    if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) 
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd > nr)) stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr) nr <- length(rowInd)
  }
  else if (is.integer(Rowv)) {
    #distr <- row.distfun(x)
    distr = distanceMatrix(t(x), metric = row.dist)
    #hcr <- hclustfun(distr)
    hcr = hclust(distr,method = row.clust)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    #distr <- distfun(x)
    distr = distanceMatrix(t(x), metric = row.dist)
    #hcr <- hclustfun(distr)
    hcr = hclust(distr,method = row.clust)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd)) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd > nc)) stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc) nc <- length(colInd)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc) stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    #distc <- distfun(if (symm) x else t(x))
    distc = distanceMatrix(x, metric = col.dist)
    #hcc <- hclustfun(distc)
    hcc = hclust(distc,method = col.clust)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    #distc <- distfun(if (symm) x else t(x))
    distc = distanceMatrix(x, metric = col.dist)
    #hcc <- hclustfun(distc)
    hcc = hclust(distc,method = col.clust)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd)) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  
  # annotation labels
  if (is.null(labRow)) 
    labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] 
    else rownames(x) 
  else labRow <- labRow[rowInd]
  if (is.null(labCol)) 
    labCol <- if (is.null(colnames(x))) (1:nc)[colInd] 
    else colnames(x) 
  else labCol <- labCol[colInd]
  # 
  if (!is.null(colRow)) colRow <- colRow[rowInd]
  if (!is.null(colCol)) colCol <- colCol[colInd]
  
  if (scale == "row") {
    #retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    #x <- sweep(x, 1, rm)
    #retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    #x <- sweep(x, 1, sx, "/")
    x = t(scale(t(x)))
    if(trim) x = trime(x,trim.size)
  }
  else if (scale == "column") {
    #retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    #x <- sweep(x, 2, rm)
    #retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    #x <- sweep(x, 2, sx, "/")
    x = scale(x)
    if(trim) x = trime(x,trim.size)
  }
  

  if (is.function(col)){
    map.x = makecmap( x, n=ncolor,breaks = breaks, colFn = col,symm = symbreaks,col.na = na.color)
    col = map.x$colors
    breaks = map.x$breaks
    } else if(length(breaks) != length(col)+1){ 
    stop("please check the length of breaks to match col")
    map.x = list(breaks = breaks,colors=col,right=FALSE,include.lowest =FALSE)
    breaks = map.x$breaks
    col = map.x$col
    }
  
  if (missing(lhei) || is.null(lhei)) lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid)) lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (nrow(ColSideColors) != nc)         #LR
        stop("nrow of 'ColSideColors' must be  ==  ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], barwidth*ncol(ColSideColors), lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (nrow(RowSideColors) != nr) 
        stop("nrow of 'RowSideColors' must be == nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], barwidth*ncol(RowSideColors), lwid[2])
    }
    n = sum(is.na(lmat))
    #lmat[is.na(lmat)] <- max(lmat,na.rm = T)+ 1:n
    #lmat[is.na(lmat)] <- 0
    lmat[is.na(lmat)] <- max(lmat,na.rm=T)
  }
  
  if (length(lhei) != nrow(lmat)) stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat)) stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  
  # start to plot now
  # the cimage is plot of image(t(x))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, gap), las = 2)
    cimage(zcol = t(RowSideColors[rowInd, , drop = FALSE]),axes = FALSE,xlab=NA,ylab=NA,border = sborder)
    axis(1, at = 1:ncol(RowSideColors), labels = colnames(RowSideColors), tick = FALSE,line = -0.5)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(gap, 0, 0, margins[2]), las = 2)
    cimage(zcol = as.matrix(ColSideColors)[colInd, , drop = FALSE],axes = FALSE,xlab=NA,ylab=NA,border = sborder)
    axis(2, at = 1:ncol(ColSideColors), labels = colnames(ColSideColors), tick = FALSE,line = -0.5)
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  

  x <- t(x)
  xcolor = cmap(x,map.x)
  
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }else iy <- 1:nr
  
  #image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks)
  cimage(zcol = xcolor,axes = FALSE,xlab=NA,ylab=NA,border = sepcolor)
  
  retval$carpet <- x
  
  if (exists("ddr")) retval$rowDendrogram <- ddr
  if (exists("ddc")) retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$color <- col
  
  
  if (is.null(srtCol) && is.null(colCol)) 
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], padj = adjCol[2])
  else {
    if (is.null(srtCol) || is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol)) 
      adjCol = c(1, NA)
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * strheight("M"), labels = labCol, adj = adjCol, cex = cexCol, srt = srtCol, col = colCol)
      print(colCol)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow) && is.null(colRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
  }
  else {
    if (is.null(srtRow) || is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2, line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
           srt = srtRow, col = colRow)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  
  # overall label for x and y
  if (!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25)
  # some annotations
  if (!missing(add.expr)) eval(substitute(add.expr))
  
  # plot seperator of plot: can be updated!!!

  if (!missing(cellnote)) text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)
  
  # plot row ddr
  par(mar = c(margins[1], margins[3], 0, gap))
  if (dendrogram %in% c("both", "row")) {
    flag <- try(plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == 
          "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
        stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }
  else {
    par(mar=c(margins[1],0,0,0))
    plot.new()
    if(key) hkeybar0(x = 0.1, y= 0.7, map.x,stretch = key.stretch,title = key.title,cex.title = keysize,cex = keysize)
    #box()
  }
  
  par(mar = c(0, 0, (tm=if (!is.null(main)) 1.2*margins[4] else margins[4]), margins[2]))
  if (dendrogram %in% c("both", "column")) {
    flag <- try(plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == 
          "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
        stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }else {
    par(mar=c(0,0,tm,margins[2]))
    plot.new()
    if (key) hkeybar0(x = 0, y= if (!is.null(main)) 0.7 else 0.6 , map.x,stretch = key.stretch,title = key.title,cex.title = keysize,cex = keysize)
    #box()
  } 
  
  if (!is.null(main)) title(main, cex.main = 1.5 * op[["cex.main"]])
  
  if (dendrogram == "both"){
    par(mar=c(0,0,0,0))
    plot.new()
    if(key)hkeybar0(x = 0.2, y= 0.5, map.x,stretch = key.stretch, title = key.title,cex = keysize,cex.title = keysize)
    #box()
  } 
  
  #par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                  high = retval$breaks[-1], color = retval$color)
  if (!is.null(extrafun)) extrafun()
  invisible(retval)
}

hkeybar0 = function (map, title = NA, side = 1, stretch = 1.4, x, y, skip, wh,cex = 1,cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  opar <- par(xpd = TRUE)
  #opar <- par(xpd = NA,mar = c(0, 0, 0, 0))
  #opar <- par(fig=c(0, 1, 0, 1),xpd = TRUE,new = TRUE)
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
  rect(xbord[-n], y, xbord[-1], y + dy, col = map$colors, border = NA)
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

hkeybar0 = function (map, title = NA, side = 1, stretch = 1.4, x, y, skip, wh,cex = 1,cex.title = 1) {
  if (!missing(skip) && !missing(wh)) 
    stop("cannot specify both 'skip' and 'wh'")
  opar <- par(xpd = TRUE)
  #opar <- par(xpd = NA,mar = c(0, 0, 0, 0))
  #opar <- par(fig=c(0, 1, 0, 1),xpd = TRUE,new = TRUE)
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
  rect(xbord[-n], y, xbord[-1], y + dy, col = map$colors, border = NA)
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


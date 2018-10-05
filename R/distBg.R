#' @import HiTC
#' @export
# Modified getExpectedCountsMean() from 
# HiTC package: ttps://rdrr.io/bioc/HiTC/src/R/normalize_hiC.R
getExpectedCountsMeanSd <- function (x, logbin = TRUE, 
                                     step = 1.05, filter.low = 0.05) {
    stopifnot(inherits(x,"HTCexp"))

    xdata <- intdata(x)
    N <- dim(xdata)[1]
    if (logbin) {
        bins <- HiTC:::logbins(from = 1, to = N, step = step)
        bins <- as.vector(Rle(values = bins, lengths = c(diff(bins),1)))
        stopifnot(length(bins) == N)
    } else {
        bins <- 1:N
    }
    message("Estimate expected using mean contact frequency per genomic distance ...")
    xdata <- as.matrix(xdata)
    rc <- colSums(xdata, na.rm = TRUE)
    print(quantile(rc[which(rc > 0)], probs = filter.low))
    rc <- which(rc <= ceiling(quantile(rc[which(rc > 0)], probs = filter.low)))
    rr <- rowSums(xdata, na.rm = TRUE)
    rr <- which(rr <= ceiling(quantile(rr[which(rr > 0)], probs = filter.low)))
    xdata[rr, ] <- NA
    xdata[, rc] <- NA
    rows <- matrix(rep.int(bins, N), nrow = N)

    d <- matrix(bins[1 + abs(col(rows) - row(rows))], nrow = N) - 1
    d[upper.tri(d)] <- -d[upper.tri(d)]
    if (isSymmetric(xdata)) { d[upper.tri(d)] <- NA }
    
    mi <- split(xdata, d)
    milen <- lapply(mi, length)

    mimean <- lapply(mi, mean, na.rm = TRUE)
    miexp <- lapply(1:length(milen), function(i) { rep(mimean[[i]], milen[[i]]) })
    names(miexp) <- names(mi)
    expmat <- as(matrix(unsplit(miexp, d), nrow = nrow(xdata), ncol = ncol(xdata)), "Matrix")
    if (isSymmetric(xdata)) { expmat <- forceSymmetric(expmat, uplo = "L") }
    colnames(expmat) <- colnames(xdata)
    rownames(expmat) <- rownames(xdata)
    expmat[rr, ] <- NA
    expmat[, rc] <- NA
    
    misd <- lapply(mi, sd , na.rm=TRUE)
    mistdev <- lapply(1:length(milen), function(i) { rep(misd[[i]], milen[[i]]) })
    names(mistdev) <- names(mi)
    sdmat <- as(matrix(unsplit(mistdev, d), nrow = nrow(xdata), ncol = ncol(xdata)), "Matrix")
    if (isSymmetric(xdata)) { sdmat <- forceSymmetric(sdmat, uplo = "L") }
    colnames(sdmat) <- colnames(xdata)
    rownames(sdmat) <- rownames(xdata)
    sdmat[rr, ] <- NA
    sdmat[, rc] <- NA
    
    return(list(exp.interaction = expmat, stdev.estimate = sdmat))
}


#' @import HiTC
#' @export
# From HiTC package: https://rdrr.io/bioc/HiTC/src/R/normalize_hiC.R
getExpectedCountsLoess <- function(x, span=0.01, bin=0.005, stdev=FALSE, plot=FALSE){
    stopifnot(inherits(x,"HTCexp"))

    xdata <- as.matrix(intdata(x))
    rc <- which(colSums(xdata, na.rm=TRUE)==0)
    rr <- which(rowSums(xdata, na.rm=TRUE)==0)

    ## rm line with only zeros
    xdata[rr,] <- NA
    xdata[,rc] <- NA
      
    ydata <- as.vector(xdata)
    ydata[which(is.na(ydata))] <- 0
    xdata.dist <- as.vector(intervalsDist(x))
    o<- order(xdata.dist)
    xdata.dist <- xdata.dist[o]
    ydata <- ydata[o]
    
    delta <- bin*diff(range(xdata.dist))
    ######################
    ## Lowess Fit
    ######################
    message("Lowess fit ...")
    #lowess.fit <- .C("lowess", x = as.double(xdata.dist), as.double(ydata), 
    #                length(ydata), as.double(span), as.integer(3), as.double(delta), 
    #                y = double(length(ydata)), double(length(ydata)), double(length(ydata)), PACKAGE = "stats")$y

    lowess.fit <-lowess(x=xdata.dist, y=ydata, f=span, delta=delta)$y
    
    y1 <- sort(ydata)
    y1 <-  quantile(y1[which(y1>1)], probs=0.99)

    if (plot){
        par(font.lab=2, mar=c(4,4,1,1))

        ##plotIntraDist(ydata, xdata.dist, xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="Counts", main="", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, col="gray", frame=FALSE)
        
        plot(x=xdata.dist, y=ydata,  xlab="Genomic Distance (bp)",  ylim=c(0,y1), ylab="Counts", main="", cex=0.5, cex.lab=0.7, pch=20, cex.axis=0.7, col="gray", frame=FALSE)
        points(x=xdata.dist[order(lowess.fit)], y=sort(lowess.fit), type="l", col="red")
    }
    lowess.mat <- Matrix(lowess.fit[order(o)], nrow=length(y_intervals(x)), byrow=FALSE)
    rownames(lowess.mat) <- id(y_intervals(x))
    colnames(lowess.mat) <- id(x_intervals(x))

    ######################
    ## Variance estimation
    ######################
    stdev.mat <- NULL
    if (stdev){
        message("Standard deviation calculation ...")
        ##interpolation
        ind <- getDeltaRange(delta, xdata.dist)
        lx <- length(xdata.dist)
        Q <- floor(lx*span)
        stdev.delta <- unlist(mclapply(1:length(ind), function(k){
            i <- ind[k]
            x1 <- xdata.dist[i]
            
            ## Neighbors selection 2*Q
            ll <- i-Q-1
            lr <- i+Q-1
            if (ll<0) ll=0
            if (lr>lx) lr=lx
            xdata.dist.sub <- xdata.dist[ll:lr]
            ydata.sub <- ydata[ll:lr]
            ## Select the Q closest distances
            d <- abs(x1-xdata.dist.sub)
            o2 <- order(d)[1:Q]
            x2 <- xdata.dist.sub[o2]
            y2 <- ydata.sub[o2]
            ## Distance between x and other points
            dref <- d[o2] 
            drefs <- dref/max(abs(dref-x1)) ##max(dref) - NS
            ## Tricube weigths and stdev calculation
            w <- tricube(drefs)
            sqrt <- w*(y2-lowess.fit[i])^2
            
            stdev <- sqrt(sum(sqrt)/
                          (((length(sqrt)-1) * sum(w))/length(sqrt)))
        }))

        if (plot){
            points(x=xdata.dist[ind], y=lowess.fit[ind], col="black", cex=.8, pch="+")
            legend(x="topright", lty=c(1,NA), pch=c(NA,"+"), col=c("red","black"),legend=c("Lowess fit","Interpolation points"), cex=.8, bty="n")
        }
        
        ## Approximation according to delta
        stdev.estimate <- approx(x=xdata.dist[ind], y=stdev.delta, method="linear", xout=xdata.dist)$y
        stdev.mat <- matrix(stdev.estimate[order(o)], nrow=length(y_intervals(x)), byrow=FALSE)
        rownames(stdev.mat) <- id(y_intervals(x))
        colnames(stdev.mat) <- id(x_intervals(x))
    }
    
    ## Put NA at rc and cc
    lowess.mat[rr,] <- NA
    lowess.mat[,rc] <- NA
 
    return(list(exp.interaction=lowess.mat,stdev.estimate=stdev.mat))
}
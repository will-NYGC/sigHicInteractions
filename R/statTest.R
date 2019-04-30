#' @export
# Binomial p-value
# 1.  use ( expected counts in bin / marginal total ) as p, correct expected counts with bias
# 2.  use ( marginal observed counts ) as n
# 3.  use ( observed bin counts ) as x
# 4.  pbinom(x,n,p) == p.value
# 5.  convert -log10 p
# bin1: bed/vcf focal point; bin2: other interacting bin to evaluate
binomP <- function(chr,bin1,bins2,obs,exp,sample) {  
        
    stopifnot(class(obs)=="HTCexp")
    stopifnot(class(exp)=="list")
        
    o <- intdata(obs)[bin1,]
    e <- exp$exp.interaction[bin1,]

    p.values <- pbinom(o,sum(o,na.rm=TRUE),e/sum(e,na.rm=TRUE),lower.tail=FALSE)
    names(p.values) <- names(e)

    q.values <- p.adjust(p.values,method='BH')
    q.values[bin1] <- NA
    q.values <- -log10(q.values[bins2])
    q.values <- as.numeric(gsub("Inf",325,q.values))
    return(q.values)
}


#' @export
# Poisson p-value
# 1.  use ( expected counts in bin / marginal total ) as p, correct expected counts with bias
# 2.  use ( marginal observed counts ) as n
# 3.  use ( observed bin counts ) as x
# 4.  ppois(x,n) == p.value
# 5.  convert -log10 p
# bin1: bed/vcf focal point; bin2: other interacting bin to evaluate
poisP <- function(chr,bin1,bins2,obs,exp,sample) {  

    stopifnot(class(obs)=="HTCexp")
    stopifnot(class(exp)=="list")
    
    o <- intdata(obs)[bin1,]
    e <- exp$exp.interaction[bin1,]
    
    p.values <- ppois(o,e,lower.tail=FALSE)
    names(p.values) <- names(e)

    q.values <- p.adjust(p.values,method='BH')
    q.values[bin1] <- NA
    q.values <- -log10(q.values[bins2])
    q.values <- as.numeric(gsub("Inf",325,q.values))
    return(q.values)
}

#' @import fitdistrplus
#' @export
# Weibull p-value
weibullP <- function(chr,bin1,bins2,obs,exp,sample) {  # bin1: bed/vcf focal point; bin2: other interacting bin to evaluate
    
    stopifnot(class(obs)=="HTCexp")
    stopifnot(class(exp)=="list")

    o <- intdata(obs)[bin1,]

    quant95 <- quantile(as.vector(exp$exp.interaction), probs = .95, na.rm = TRUE)
    exp$exp.interaction[exp$exp.interaction > quant95] <- NA
    exp$exp.interaction[exp$exp.interaction == 0] <- NA
    
    e.mu <- exp$exp.interaction[bin1,]
    e.sd <- exp$stdev.estimate[bin1,]
    z <- (o-e.mu)/e.sd
    fit <- fitdistrplus::fitdist(e.mu[!is.na(e.mu)], 'weibull')

    # if (runif(1) > .95) {  # Plot 5% of the Weibull fits for interactions randomly selected

    #     random <- rweibull(sum(!is.na(z)),shape=fit[[1]][1],scale=fit[[1]][2])
    #     h.z <- hist(z,plot=FALSE,breaks=30)
    #     span <- h.z$breaks[2] - h.z$breaks[1]
    #     h.random <- hist(random,plot=FALSE,breaks=seq(min(h.z$breaks,min(random-span)),max(h.z$breaks,max(random+span)),span))
        
    #     png(file=paste(opt$prefix,paste0('bin',bin1),'weibullFit.png'),height=400,width=600)
    #     plot(h.z$mids,h.z$counts,col='red',type='b',main=paste("Weibull fit for",sample,paste0("bin",bin1)),
    #          xlab='Z-score bins',ylab='Counts')
    #     points(h.random$mids,h.random$counts,col='black',type='b')
    #     legend('topright',legend=c('Observed','Fit'),lty=1,col=c('red','black'))
    #     dev.off()
    # }
        
    p.values <- pweibull(o,shape=fit$estimate[1],scale=fit$estimate[2],lower.tail=FALSE)
    names(p.values) <- names(e.mu)

    q.values <- p.adjust(p.values,method='BH')
    q.values[bin1] <- NA
    q.values <- -log10(q.values[bins2])
    q.values <- as.numeric(gsub("Inf",325,q.values))
    return(q.values)
}

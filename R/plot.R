#' @import ggplot2 gridExtra RColorBrewer TxDb.Hsapiens.UCSC.hg19.knownGene org.Hs.eg.db TxDb.Mmusculus.UCSC.mm10.knownGene org.Mm.eg.db plyr data.table reshape2 csaw
#' @export
makePlot <- function(range,region,p,transcript,exon,domains,anno,loci,abs,output,cols = NULL, cols2 = NULL) {

    # Specifically ordered for RajarajanP
    if (is.null(cols)) cols <- brewer.pal(n = length(unique(p$variable)), name = 'Set1')
    if (is.null(cols2)) cols2 <- brewer.pal(n = length(unique(anno$state)), name = 'Set2')
    if (!file.exists(dirname(output))) { dir.create(dirname(output),recursive=TRUE) }
    
    rsid <- range$id
    min <- min(start(region))
    max <- max(end(region))

    window <- end(abs[1]) - start(abs[1]) + 1
    breaks <- seq(window,max(end(abs[as.vector(seqnames(abs))==as.vector(seqnames(region)),])),
                  max(window,50000))  #Min 5kb grid width

    # heights: domains / loci / p-values / anno / genes
    heights <- c(2,1,8,max(2,max(anno$height)*.5),max(2,max(transcript$height)/10))
    
    if (nrow(domains)>0) {
        plot.domains <- ggplot(data=domains,aes(xmin=start,xmax=end,ymin=height-.3,
                                   ymax=height+.3,fill=sample,alpha=alpha)) +
                                       geom_rect(colour='black',size=.1)
    } else {
        plot.domains <- ggplot(data=domains,aes(xmin=start,xmax=end)) +
            geom_rect(colour='black',size=.1)
    }
    plot.domains <- plot.domains + scale_fill_manual(values=cols,drop=FALSE) +
        scale_alpha(range=c(.75,1)) +
        scale_x_continuous(expand=c(0,0)) +
        ylab("Domains") +
        coord_cartesian(xlim=c(min,max)) +
        guides(alpha=FALSE) +
        theme(legend.position='top',
              legend.title=element_blank(),
              legend.key.size=unit(1.5,'lines'),
              legend.text=element_text(size=14),
              panel.background=element_blank(),
              axis.text=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks=element_blank())

    plot.loci <- ggplot(data=loci,aes(xmin=start,xmax=end,ymin=height+.35,ymax=height-.35,fill=focus)) +
        geom_rect(colour='black',size=.25) +
        geom_text(aes(label=id,x=end,y=height),size=5,hjust=0,colour='black') +
        scale_fill_manual(values=c('red','lightpink')) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(name="ROI",expand=c(0,0),limits=c(.5,2.5)) +
        coord_cartesian(xlim=c(min,max)) +
        theme(legend.position='none',
              panel.background=element_blank(),
              axis.text=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks=element_blank())
    
    plot.p <- ggplot(data=p,aes(x=(start+end)/2,y=value,colour=variable))

    if (max - min <= 5e6) {
        plot.p <- plot.p + geom_line(size=1.25,alpha=.65) + geom_point(size=2.5,alpha=.65) +
            scale_x_continuous(name=paste(seqnames(range),"(Mb)"),
                               breaks=breaks,labels=sprintf("%.2f",breaks/1e6),expand=c(0,0))
    } else {
        plot.p <- plot.p + geom_line(size=1,alpha=.85) +
            scale_x_continuous(name=paste(seqnames(range),"(Mb)"),
                               breaks=seq(min(breaks),max(breaks),5e6),
                               labels=sprintf("%.0f",seq(min(breaks),max(breaks),5e6)/1e6),expand=c(0,0))
    }
    y.max <- max(ggplot_build(plot.p)$layout$panel_scales_y[[1]]$range$range[2], 5)
    
    plot.p <- plot.p +
        scale_colour_manual(values=cols,drop=FALSE) +            
        scale_y_continuous(name="-log q-value", limits=c(0, y.max)) +
        coord_cartesian(xlim=c(min,max)) +
        theme(legend.position='none',
              panel.background=element_blank(),
              panel.grid.major.x=element_line(size=.25,linetype='dotted',colour='black'),
              panel.grid.major.y=element_line(size=.25,linetype='dotted',colour='black'),
              axis.text.x=element_text(size=14,colour='black',angle=45,vjust=1,hjust=1),
              axis.text.y=element_blank(),
              axis.title.x=element_text(size=24,face='bold',margin=margin(t=12,b=0)),
              axis.title.y=element_text(size=16),
              axis.ticks.x=element_line(colour='black'),
              axis.ticks.y=element_blank())

    if (end(range)-start(range) == 0) {
        plot.p <- plot.p + geom_vline(xintercept=start(range),size=.75,linetype='dashed',colour='red')
    } else {
        plot.p <- plot.p +
            annotate('rect', size=.75,linetype='dashed',
                     colour='red',fill='red',alpha=.1,
                     xmin=start(range),xmax=end(range),
                     ymin=0,ymax=y.max)
    }
    plot.p <- plot.p + annotate("text",size=5,hjust=0,vjust=-.5,
                                x=min,                                                                                           
                                y=y.max,
                                hjust = -.35,
                                colour='black',
                                fontface = 'bold',
                                label=sprintf("max = %.2f", y.max)) +
                       geom_hline(size=.25,
                                  yintercept = y.max,
                                  colour='black',
                                  linetype = 'dotted')
                        
    
    plot.anno <- ggplot(data=anno, aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax,fill=state)) +
        geom_rect(colour=NA,size=.1) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(name="Annotations",limits=c(min(anno$ymin),max(anno$ymax)+.65)) +
        scale_fill_manual(values=cols2) +
#        scale_fill_manual(values=c('lightcoral','cornflowerblue','forestgreen','gray80','darkorange','violet')) +
        coord_cartesian(xlim=c(min,max)) +
        theme(legend.position='bottom',
              legend.title=element_blank(),
              legend.key.size=unit(1,'lines'),
              legend.text=element_text(size=12),
              panel.background=element_blank(),
              axis.text=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks=element_blank())
    
    for (i in 1:length(unique(anno$sample))) {
        plot.anno <- plot.anno + annotate("text",x=(min+max)/2,y=i+.65,label=unique(anno$sample)[i],size=3.5,vjust=1)
    }

    plot.genes <- ggplot(data=transcript,aes(x=start,xend=end,y=height,yend=height)) + geom_segment(colour='blue',size=.2) +
        geom_text(size=4.5,colour='blue',vjust=-.5,#position=position_jitter(w=0,h=max(transcript$height)/3),
                  aes(x=(start+end)/2,y=height+.15,label=gene_id)) +
        geom_rect(data=exon,colour='blue',fill='blue',aes(xmin=start,xmax=end,ymin=height-.3,ymax=height+.3)) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(name="Gene models",expand=c(.15,0)) +
        coord_cartesian(xlim=c(min,max)) +
        theme(panel.background=element_blank(),
              axis.text=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks=element_blank())
    
    pdf(file=output,height=sum(heights),width=24)
    grid.arrange(plot.domains,plot.loci,plot.p,plot.anno,plot.genes,ncol=1,
                 layout_matrix=matrix(c(rep(1,heights[1]),
                     rep(2,heights[2]),
                     rep(3,heights[3]),
                     rep(4,heights[4]),
                     rep(5,heights[5]))))
    dev.off()
}

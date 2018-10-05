#' @import GenomicRanges
#' @export
getBins <- function(query ,  abs ,  midpoint = FALSE) {

    stopifnot(inherits(query,"GRanges"))
    stopifnot(inherits(abs,"GRanges"))

    if (midpoint == TRUE) {
        cat("Using midpoint of ranges to define bins.\n", file = stderr())
        start(query) <- end(query) <- (start(query)+end(query))/2
    }
    hits <- GenomicRanges::findOverlaps(query = query, subject = abs)
    return(mcols(abs)$name[subjectHits(hits)])
}

#' @import GenomicRanges data.table
#' @export
# Assumes a GTF-like object, uses gr$source to reduce, otherwise must specify
collapseRanges <- function(gr,key='ucsc_id', flank = NULL) {  
    
    stopifnot(inherits(gr,"GRanges"))

    gr$key <- paste(paste(mcols(gr)$source, mcols(gr)[[key]]))
    
    tx.gr <- if ('type' %in% names(mcols(gr)) && 'transcript' %in% gr$type) {
        tx.gr <- gr[gr$type=='transcript'] 
    } else {
        tx.gr <- gr
    }

    # Pick longest isoform
    longest.gr <- GRanges()
    for (id in unique(mcols(tx.gr)[[key]])) {
        id.gr <- tx.gr[ mcols(tx.gr)[[key]] == id ]
        longest.gr <- c(longest.gr, id.gr[order(width(id.gr)[1], decreasing= TRUE)])
    }
    
    new.gr <- GRanges()
    while (length(longest.gr)>0) {
        first.gr <- longest.gr[1]
        # first.gr <- range(longest.gr[mcols(longest.gr)[[key]]==mcols(longest.gr[1])[[key]]])
        # mcols(first.gr) <- mcols(longest.gr[1])
        first.gr$source <- sprintf("%s:%d-%d(%s)",seqnames(first.gr),start(first.gr),end(first.gr),strand(first.gr))
        
        expanded.first.gr <- first.gr
        expanded.longest.gr <- longest.gr
        if (!is.null(flank)) {
            start(expanded.first.gr) <- start(expanded.first.gr) - flank
            end(expanded.first.gr) <- end(expanded.first.gr) + flank
            start(expanded.longest.gr) <- start(expanded.longest.gr) - flank
            end(expanded.longest.gr) <- end(expanded.longest.gr) + flank
        }

        overlaps <- findOverlaps(expanded.first.gr,expanded.longest.gr,ignore.strand=TRUE)
        group.gr <- longest.gr[ subjectHits(overlaps[queryHits(overlaps)==1]) ]
        group.gr$height <- 3 * (as.numeric(factor(mcols(group.gr)[[key]])) - 1)
        new.gr <- append(new.gr,group.gr)
        longest.gr <- longest.gr[!(mcols(longest.gr)[[key]] %in% mcols(group.gr)[[key]])]
    }

    if ('type' %in% names(mcols(gr)) && 'exon' %in% gr$type) {
        exon.gr <- gr[gr$type=='exon']
        dt <- data.table(unique(data.frame(new.gr)[,c('key', key, 'height')]))
        setkeyv(dt, key)
        mcols(exon.gr)$height <- sapply(mcols(exon.gr)[[key]], function(x) dt[x]$height)
        new.gr <- c(new.gr,exon.gr)
    } 
    return(new.gr)
}
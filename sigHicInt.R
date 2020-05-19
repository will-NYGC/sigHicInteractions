#!/usr/bin/env Rscript
suppressMessages(library(argparse))


# Get command-line parameters
parser <- argparse::ArgumentParser()

parser$add_argument("-m", "--matrix", dest = "matrix",
                    type = "character", required = TRUE, nargs = "+",
                    help = "HiC-Pro 3-column raw sparse matrix [1 or more]")
parser$add_argument("-n", "--names", dest = "names",
                    type = "character", required = FALSE, nargs = "+",
                    help = "Names to associate with matrix files of same length as -m [Default: strips extension of matrix file basenames]")
parser$add_argument("-k", "--colors", dest = "colors",
                    type = "character", required = FALSE, nargs = "*", default = NULL,
                    help = "Colors to associate with Names [Default: brewer \"Set1\" palette]")
parser$add_argument("-a", "--abs", dest = "abs",
                    type = "character", required = TRUE,
                    help = "HiC-Pro \"abs.bed\" file describing bins")
parser$add_argument("-r", "--roi", dest = "roi",
                    type = "character", required = TRUE,
                    help = "Genomic loci to bin and test for significant interactions")
parser$add_argument("-l", "--roi_label", dest = "roi_label",
                    type = "character", required = FALSE, default = "ROI",
                    help = "Label for ROI's in plot [Default: \"ROI\"]")
parser$add_argument("-u", "--roi_column", dest = "roi_column",
                    type = "character", required = FALSE, default = NULL,
                    help = "GRanges mcols column to use as ID")
parser$add_argument("-o", "--output_prefix", dest = "output_prefix",
                    type = "character", required = FALSE, default = 'sigInt',
                    help = "Output prefix [Default: \'sigInt\']")
parser$add_argument("-d", "--domains", dest = "domains",
                    type = "character", required = FALSE, nargs = "*",
                    help = "Domain annotations, e.g. TADs or A compartments, using the first \"word\" in filename as label [Default: Rajarajan, P. Glia/NPC/Neuron 500kb res A-compartment calls]")
parser$add_argument("-e", "--anno", dest = "anno",
                    type = "character", required = FALSE, nargs = "*",
                    help = "Genome annotations, e.g. ChromHMM, using the first \"word\" in filename as label [Default: REMC E081/E082 Fetal brain ChromHMM]")
parser$add_argument("-i", "--anno_label", dest = "anno_label",
                    type = "character", required = FALSE, nargs = "+",
                    help = "Labels for annotation tracks")
parser$add_argument("-q", "--anno_colors", dest = "anno_colors",
                    type = "character", required = FALSE, default = NULL, nargs = "+",
                    help = "Default annotation color [Default: brewer.pal(n, \"Dark2\")")
parser$add_argument("-f", "--flank", dest = "flank",
                    type = "integer", required = FALSE, default = 1e6,
                    help = "Flanking region to evaluate for significance [Default: 1e6]")                    
parser$add_argument("-j", "--genes", dest = "genes",
                    type = "character", required = FALSE,
                    help = "Gene annotations, [Default: UCSC hg19 knownCanonical (custom GTF-like format)]")                    
parser$add_argument("-s", "--stat", dest = "stat",
                    type = "character", required = FALSE,
                    default = "binom", choices = c('binom', 'pois', 'weibull'),
                    help = "Significance test [binom (default) / pois / weibull]")
parser$add_argument("-b", "--distBg", dest = "distBg",
                    type = "character", required = FALSE,
                    default = "meanSd", choices = c('meanSd', 'loess'),
                    help = "Method to estimate background interaction distribution [meanSd (default) / loess]")
parser$add_argument("-c", "--chroms", dest = "chroms", nargs = "*",
                    type = "character", required = FALSE,                    
                    help = "Chromosomes to process")
parser$add_argument("-z", "--force", dest = "force",
                    default = FALSE, action = 'store_true',
                    help = "Overwrite *.scores.txt file if it exists")
parser$add_argument("-t", "--top", dest = "top",
                    type = "integer", required = FALSE,
                    help = "Only process <top> BED Regions of Interest file")

args <- parser$parse_args()
saveRDS(args, paste0(args$output_prefix, ".args.RDS"))

# Set-up
library(sigHicInteractions)
library(rtracklayer)
library(GenomicRanges)
library(HiTC)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(reshape2)
library(plyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

sparse.list <- args$matrix

# Get names to associate with matrices
if (!is.null(args$names)) {
    names <- args$names
} else {
    names <- gsub("\\.\\w+$", "", basename(args$matrix))
}

# Create output directory
if ( !file.exists(dirname(args$output_prefix)) ) { 
    dir.create(dirname(args$output_prefix),recursive=TRUE)
}

# Save parameters for debugging
saveRDS(args, paste0(args$output_prefix, ".args.RDS"))

# Create scores file to append to
# Fail if non-empty b/c appending scores to the file
# It should be empty or non-existent.
if (file.exists(paste0(args$output_prefix,'.scores.txt')) &
    file.info( paste0(args$output_prefix,'.scores.txt') )$size != 0 ) {

    if (args$force == TRUE) {
        file.remove(paste0(args$output_prefix, '.scores.txt'))
    } else {
        stop(sprintf("Scores file %s exist. Will not overwrite. Delete to proceed.",
                     paste0(args$output_prefix,'.scores.txt')))
    }
} else {
    file.create(paste0(args$output_prefix,'.scores.txt'))
}


# Load ROI from BED file
ranges <- import.bed(args$roi)
if (!is.null(args$roi_column)) {
    mcols(ranges)$id <- names(ranges) <- mcols(ranges)[[args$roi_column]]
} else {
    mcols(ranges)$id <- names(ranges) <- sprintf("%s:%d-%d",seqnames(ranges),start(ranges),end(ranges))
}

# args$top or args$chroms, then subselect ranges
if (!is.null(args$top)) { ranges <- ranges[seq(args$top)] }
if (!is.null(args$chroms) ) {
    chroms <- args$chroms
    message(sprintf("Processing %s...", paste(chroms, collapse=",")))
} else {
    chroms <- unique(seqnames(ranges))
}


# Get Hi-C bins from HiC-Pro abs.bed file
abs <- import.bed(args$abs)


# Get Gene models
if (is.null(args$genes)) {
    args$genes <- system.file("extdata",
                              "hg19_ucsc_knownCanonical.v2.gtf",
                              package = "sigHicInteractions")
    db <- org.Hs.eg.db
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
} else if (!is.null(args$genes) & grepl("mm10", args$genes)) {
    db <- org.Mm.eg.db
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
} else if (!is.null(args$genes) & grepl("hg38", args$genes)) {
    db <- org.Hs.eg.db
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
} else {
    db <- org.Hs.eg.db
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
}
genes <- import.gff2(con=args$genes,format='gff2')


# Get Domain ranges
if (is.null(args$domains)) {
    args$domains <- c(system.file("extdata",
                                  "Glia.bin500kb.domainA.bed",
                                  package = "sigHicInteractions"),
                      system.file("extdata",
                                  "Neuron.bin500kb.domainA.bed",
                                  package = "sigHicInteractions"),
                      system.file("extdata",
                                  "NPC.bin500kb.domainA.bed",
                                  package = "sigHicInteractions"))
}
domains <- lapply(args$domains, import.bed)
# Use first "word" of file basename as label
names(domains) <- gsub("^(\\w+)?\\..*", "\\1", basename(args$domains))
domains <- lapply(seq_along(domains),function (x) {
    y <- domains[[x]]
    mcols(y)$sample <- names[x]
    return(y)
})

# Concatenate all regions into a single GRanges
domains.gr <- Reduce(c,domains)


# Get Annotations (e.g., ChromHMM)
if (is.null(args$anno)) {
    args$anno <- c(system.file("extdata",
                               "E081_Fetal_Brain_Male.6_coreMarks_statebyline.txt.gz",
                               package = "sigHicInteractions"),
                   system.file("extdata",
                               "E082_Fetal_Brain_Female.6_coreMarks_statebyline.txt.gz",
                               package = "sigHicInteractions"))
}
anno <- lapply(args$anno,function(x) read.table(gzfile(x,open='r'),sep="\t"))
# Use first "word" of file basename as label
if (is.null(args$anno_label)) {
    names(anno) <- gsub("^(.+?)(\\.|_).*", "\\1", basename(args$anno))
} else {
    names(anno) <- args$anno_label
}
anno <- lapply(seq_along(anno),function(x) GRanges(seqnames=anno[[x]][,1],
                                                   IRanges(start=anno[[x]][,2],
                                                           end=anno[[x]][,3]),
                                                   state=anno[[x]][,4],
                                                   sample=names(anno)[x]))
# Concatenate annotations into single GRanges
anno.gr <- Reduce(c,anno)


# Get statistical test
if (args$stat == "binom") {
    test <- binomP
} else if (args$stat == "pois") {
    test <- poisP
} else if (args$stat == "weibull") {
    test <- weibullP
} else {
    stop("Must indicate whether to test with binom/pois/weibull.")
}

# Load HiC contact matrix
obs.full <- lapply(sparse.list,function(x) importC(x,xgi.bed=args$abs,rm.trans=TRUE,lazyload=TRUE))

# Loop through ROI's and get sig scores
for (chrom in chroms) {

    # Create ROI subset of full ranges by chromosome
    ranges2 <- ranges[seqnames(ranges)==chrom,]

    # Subset Hi-C matrix by chrom, convert diag --> symmetrtic matrix
    obs <- lapply(obs.full,function (x) reduce(x,chr=chrom,cis=TRUE,trans=FALSE))
    obs <- lapply(obs,function(x) HTClist(lapply(x,function(y) forceSymmetric(y))))

    # Get expected counts accordinding to args$distBg
    if (args$distBg=='meanSd') {
        # Custom version of HiTC:::getExpectedCountsMean that also returns stdev
        exp <- lapply(obs,function (x) lapply(x,getExpectedCountsMeanSd)) 
    } else if ( args$distBg=='loess') {
        exp <- lapply(obs,function (x) lapply(x,getExpectedCounts,method='loess',
                                              stdev=TRUE,asList=TRUE,plot=FALSE))
    }
    
    names(obs) <- names(exp) <- names

    # Iterate over ranges2, ROI's from chrom
    for (i in seq_along(ranges2)) {

        # GRanges object, ranges3, with a single range --> chop into bins
        ranges3 <- ranges2[i]      
        bins1 <- getBins(ranges3,abs)

        # Iterate over each bin as anchor and evaluate interaction with all target bins in region
        for (bin in seq_along(bins1)) {
            
            suffix <- paste0("bin",bin,"of",length(bins1))   
                                        # This is the anchor bin and all interactions within flanks to this will be evaluated         
            bin1 <- bins1[bin]
            range <- abs[abs$name==bin1]
            range$id <- ranges2[i]$id

            # Define span of target bins to pair with anchor (bin1)
            region <- ranges3
            start(region) <- max(1, mean(c(start(ranges3), end(ranges3))) - as.numeric(args$flank))  # Shift by left and right flank
            end(region) <- mean(c(start(ranges3), end(ranges3))) + as.numeric(args$flank)  
            c <- seqnames(region)
            
            # Establish minimum plot width of ROI loci, ie if too narrow, make it look wider
            min.width <- .0025 * (end(region) - start(region))
            
            # Extract bin IDs from regions and annotate with overlapping genes
            bins2 <- getBins(region, abs)
            bins2.genes <- csaw::detailRanges(abs[abs$name %in% bins2],orgdb=db, txdb=txdb)[[1]]
            names(bins2.genes) <- bins2
            bins2.genes <- as.list(bins2.genes)
            
            # Run test() on each sample in obs
            p.matrix <- lapply(seq_along(obs),
                               function(x) {
                                   test(c,bin1,bins2,
                                        obs[[x]][[paste0(c,c)]],
                                        exp[[x]][[paste0(c,c)]],
                                        names(obs)[1])
                               })
            names(p.matrix) <- names(obs)

            # Convert to GRanges, then data.frame and dump scores to file
            p.gr <- abs[abs$name %in% bins2, ]
            mcols(p.gr) <- cbind(p.gr$name,data.frame(p.matrix))   
            colnames(mcols(p.gr)) <- c('name', names)         
            lines <- t( apply(data.frame(p.gr),1,function(x) {
                return( c(range$id,
                          paste(bins1,collapse=","),
                          bin1,
                          sprintf("%s:%s-%s",seqnames(range),start(range),end(range)),
                          bins2.genes[[bin1]],
                          x[6],
                          sprintf("%s:%s-%s",x[1],x[2],x[3]),
                          bins2.genes[[as.character(x[6])]],
                          x[7:ncol(data.frame(p.gr))]) )
            }))
            colnames(lines) <- c('ROI','ROI.bins','anchor.bin','anchor.bin.coord','anchor.bin.genes',
                                    'target.bin','target.bin.coord','target.bin.genes',names(p.matrix))

            write.table(lines,file=paste0(args$output_prefix,'.scores.txt'),sep="\t",
                        col.names=ifelse(file.exists(paste0(args$output_prefix,'.scores.txt')) &
                            file.info( paste0(args$output_prefix,'.scores.txt') )$size != 0,
                            FALSE, TRUE),
                        row.names=FALSE,quote=FALSE,append=TRUE)                    
            
            # Create GRanges of bins in region and append p-values to mcols()
            mcols(p.gr) <- data.frame(p.matrix)
            p.df <- data.frame(p.gr)
            colnames(p.df) <- c('seqnames', 'start', 'end', 'width', 'strand', names)
  
            # Melt p-values GRanges to a data.frame for ggplot
            p.df <- melt(data=p.df,id.vars=c(1:3),measure.vars=c(6:ncol(data.frame(p.df))))

            # For RajarajanP, ordered labels
            p.df$variable <- factor(p.df$variable, levels=names(obs))
            
            # Get genes in region and convert into transcript and exon data.frames
            sub <- subsetByOverlaps(genes,region)
            if (length(sub) > 0) {
                key <- ifelse('gene_id' %in% colnames(mcols(sub)), 'gene_id', 'id')
                cr <- data.frame(collapseRanges(sub, key = key,
                                 flank = width(region)/20))
                genes.df <- cr[c(1:3,7,10:ncol(cr))]   
#                genes.df$height <- as.numeric(as.factor(genes.df[[key]]))     
            } else {
                genes.df <- data.frame(seqnames = rep(NA,2),
                                        start = rep(0,2),
                                        end = rep(0,2),
                                        type = c('transcript', 'exon'),
                                        gene_id = rep(NA, 2),
                                        ucsc_id = rep(NA, 2),
                                        ensembl_id = rep(NA, 2),
                                        height = rep(1,2))
            }
            transcript.df <- genes.df[genes.df$type=='transcript',]
            exon.df <- genes.df[genes.df$type=='exon',]
            
            # Get Domains in region and convert into data.frame
            domains.df <- data.frame(subsetByOverlaps(domains.gr,region))
            # For RajarajanP, ordered labels
            domains.df$sample <- factor(domains.df$sample,levels=names)
            domains.df$height <- sapply(domains.df$sample, function(x) which(names == x))
            domains.df <- ddply(domains.df, .(sample), function(x) cbind(x,rep(c(0,1), length.out=nrow(x))))
            colnames(domains.df)[ncol(domains.df)] <- 'alpha'
            
            # Get annotations in region and convert into data.frame
            anno.df <- data.frame(subsetByOverlaps(anno.gr,region))
            states <- levels(anno.gr$state)
            anno.df$state <- factor(anno.df$state,levels=states[order(as.numeric(gsub("(\\d+)\\..*","\\1",states)))])
            anno.df$height <- as.numeric(factor(anno.df$sample))

            # For extdata ChromHMM calls from REMC, make Enhancer annotations bigger so more visible
            anno.df$ymin <- apply(anno.df,1,function(x) {
                if (x[6] == '3. Enhancer') {
                    return(as.numeric(x[8])-.3)
                } else {
                    return(as.numeric(x[8])-.2)
                }
            })
            anno.df$ymax <- apply(anno.df,1,function(x) {
                if (x[6] == '3. Enhancer') {
                    return(as.numeric(x[8])+.3)
                } else {
                    return(as.numeric(x[8])+.2)
                }
            })
            
            # Get anchor BED ranges2, e.g. PGC haplotype and convert into data.frame
            # Widen really narrow ROI's using min.width
            loci.gr <- subsetByOverlaps(ranges2,region)
            short.idx <- which(width(loci.gr) < min.width)
            loci.gr[short.idx] <- resize(loci.gr[short.idx], width = min.width)
            
            range.gr <- loci.gr[loci.gr$id==range$id]
            loci.gr <- loci.gr[loci.gr$id!=range$id]
            
            # If focusing on 1 ROI but another ROI falls into region, color the differently
            if (length(loci.gr) > 0) {
                id <- loci.gr$id
                mcols(loci.gr) <- NULL
                if (length(loci.gr)==1) {
                    loci.gr$id <- paste0(" ",id)                
                } else {
                    loci.gr$id[length(loci.gr)] <- sprintf(" %s", args$roi_label)
                }
                loci.gr$height <- 2
            }
            
            mcols(range.gr) <- NULL
            range.gr$id <- paste0(' ',names(range.gr))
            range.gr$height <- 1
            
            loci.df <- data.frame(c(range.gr,loci.gr))
            loci.df$focus <- factor(apply(loci.df,1,function(x) {
                if (!is.na(x[["id"]])) {
                    return(ifelse(x[["id"]]!=paste0(" ",range$id),2,1))
                } else {
                    return(2)
                }
            }))

            # Make plot
            cat(sprintf("%d. Plotting: %s, %s\n",i,range$id,suffix),file=stderr())            
            output <- paste(args$output_prefix,
                            paste0(ranges2[i],"/",basename(args$output_prefix),".",ranges2[i]),
                            suffix,"pdf",sep=".")
            makePlot(range,region,p.df,transcript.df,exon.df,domains.df,anno.df,loci.df,abs,output,args$colors, args$anno_colors)
        }
    }
}

#!/usr/bin/env Rscript
args <- commandArgs()
scriptPath <- strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName <- basename(scriptPath)
args <- args[-(1:5)]
usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -g=refGenome -gtf=gtf1,gtf2 -mb=mixedBam -hb=haploBam -c=chrom -cs=chromStart -ce=chromEnd
Option:
    Commmon:
    -f|refGenome  str     The path to file contains the genes you want to enrich, or the file list separated by comma.
    -g|gtfs      str     The sample name you want to specify whose order should consistent with the gene file(s).
    -mb|mixedBam         str     The gene2go file which used server as background in the enrichment.
    -hb|haploBams         str     The output file name the plot will be stored.
    -c|chrom    str
    -cs|chromStart   int
    -ce|chromEnd    int
")
    q(save="no")
}
parseArg = function(arg, pattern, msg){
    pattern = paste0('^-', pattern, '=')
    if(grepl(pattern, arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        if(is.na(arg.split[2])){
            stop(paste0('Please specify the value of -', msg))
        }else{
            return(arg.split[2])
        }
    }
}

refGenome <- ""
gtfs <- ""
mixedBam <- ""
haploBams <- ""
chrom <- ""
chromStart <- 0
chromEnd <- 0

if(length(args) >= 1){
    for(i in 1:length(args)){
        arg <- args[i]
        if(arg == '-h' || arg == '-help') usage()
        tmp <- parseArg(arg, 'f', 'refGenome'); if(!is.null(tmp)) refGenome <- tmp
        tmp <- parseArg(arg, 'g', 'gtfs'); if(!is.null(tmp)) gtfs <- tmp
        tmp <- parseArg(arg, 'mb', 'mixedBam'); if(!is.null(tmp)) mixedBam <- tmp
        tmp <- parseArg(arg, 'hb', 'haploBams'); if(!is.null(tmp)) haploBams <- tmp
        tmp <- parseArg(arg, 'c', 'chrom'); if(!is.null(tmp)) chrom <- tmp
        tmp <- parseArg(arg, 'cs', 'chromStart'); if(!is.null(tmp)) chromStart <- as.numeric(tmp)
        tmp <- parseArg(arg, 'ce', 'chromEnd'); if(!is.null(tmp)) chromEnd <- as.numeric(tmp)
    }
}


library(rtracklayer)
library(GenomicFeatures)
library(Gviz)

plotAllelicAsStructure <- function(refGenome, gtfs, mixedBam, haploBams, chrom, chromStart, chromEnd, title){
    gtrack <- GenomeAxisTrack(cex = 1)
    sTrack <- SequenceTrack(refGenome)
    gtfs <- unlist(strsplit(gtfs, ","))
    haploIso1Track <- GeneRegionTrack(gtfs[1], name="Haplotype1", transcriptAnnotation = "transcript", stackHeight = 0.5)
    haploIso2Track <- GeneRegionTrack(gtfs[2], name="Haplotype2", transcriptAnnotation = "transcript", stackHeight = 0.5)
    covTrack <- AlignmentsTrack(mixedBam, isPaired = TRUE, min.height = 20, size=20, type="coverage")
    haploBams <- unlist(strsplit(haploBams, ","))
    haploReads1Track <- AlignmentsTrack(haploBams[1], isPaired = TRUE, min.height = 20, size = 80, type="pileup")
    haploReads2Track <- AlignmentsTrack(haploBams[2], isPaired = TRUE, min.height = 20, size = 80, type="pileup")
    plotTracks(c(covTrack, haploIso1Track, haploReads1Track, haploIso2Track, haploReads2Track, sTrack, gtrack), chromosome = chrom, from = chromStart, to = chromEnd, cex=0.5, min.height=2, fontsize=10, extend.left=0.02, extend.right=0.02, main=title, cex.main=1)
}

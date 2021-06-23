#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: plotRscriptStrs.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-05-13 21:53:28
Last modified: 2021-05-13 21:53:28
'''

plotReadsCorrectedEvalStr = '''
    library(ggplot2)
    library(scales)
    plotReadsCorrectedEval <- function(rawMappedBed, correctMappedBed, outPdf){
        rawMapped <- read.csv(rawMappedBed, sep="\t", header=FALSE)
        correctMapped <- read.csv(correctMappedBed, sep="\t", header=FALSE)
        rawMapped$accuracy <- (rawMapped$V13 * rawMapped$V14)/100
        correctMapped$accuracy <- (correctMapped$V13 * correctMapped$V14)/100
        rawMapped$type <- "Raw"
        correctMapped$type <- "Corrected"
        pdf(outPdf)
        p <- ggplot(rbind(rawMapped,correctMapped), aes(x=accuracy, fill=type)) + 
        geom_histogram(alpha=0.6, position = 'identity', binwidth=2) + 
        scale_fill_manual(values=c("#69b3a2", "#404080")) + 
        ggtitle("Reads Correction Evaluation") + xlab("Accuracy (%)") + 
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5, size=20), text = element_text(size=12), 
              panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"), strip.background = element_blank(), 
              strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8), 
              legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines")) +
        scale_y_continuous(labels = label_number(scale = 1 / 1e4))
        print(p)
        dev.off()
    }
'''

plotAllelicAsStructureStr = '''
    library(rtracklayer)
    library(GenomicFeatures)
    library(Gviz)
    
    plotAllelicAsStructure <- function(refGenome, gtfs, mixedBam, haploBams, ngsBam, chrom, chromStart, chromEnd, htStarts, htEnds, color, outName){
        options(ucscChromosomeNames=FALSE)
        pdf(paste0(outName, ".pdf"), height=15)
        gtrack <- GenomeAxisTrack(cex = 1)
        sTrack <- SequenceTrack(refGenome)
        gtfs <- unlist(strsplit(gtfs, ","))
        gtf1_db <- makeTxDbFromGFF(gtfs[1], format = "gtf")
        gtf2_db <- makeTxDbFromGFF(gtfs[2], format = "gtf")
        haploIso1Track <- GeneRegionTrack(gtf1_db, name="Haplotype1 Isoforms", transcriptAnnotation = "transcript", stackHeight = 0.5)
        haploIso2Track <- GeneRegionTrack(gtf2_db, name="Haplotype2 Isoforms", transcriptAnnotation = "transcript", stackHeight = 0.5)
        covTrack <- AlignmentsTrack(mixedBam, isPaired = FALSE, type="coverage", name="Coverage")
        haploBams <- unlist(strsplit(haploBams, ","))
        haploReads1Track <- AlignmentsTrack(haploBams[1], isPaired = FALSE, type="pileup", name="Haplotype1 Alignments", stackHeight = 0.5)
        haploReads2Track <- AlignmentsTrack(haploBams[2], isPaired = FALSE, type="pileup", name="Haplotype2 Alignments", stackHeight = 0.5)
        htStarts <- as.numeric(unlist(strsplit(htStarts, ","))) - 50
        htEnds <- as.numeric(unlist(strsplit(htEnds, ","))) + 50
        ht <- HighlightTrack(trackList=list(haploIso1Track, haploReads1Track, haploIso2Track, haploReads2Track), start=htStarts, end=htEnds, chromosome=chrom, col=c(color), alpha=0.8)
        if(ngsBam!=""){
            # ngsCovTrack <- AlignmentsTrack(ngsBam, isPaired = TRUE, min.height = 20, size = 80, type="coverage", name="Short Reads Alignments")
            # plotTracks(c(covTrack, ngsCovTrack, ht, sTrack, gtrack), chromosome = chrom, from = chromStart, to = chromEnd, cex=0.5, min.height=2, fontsize=10, extend.left=0.02, extend.right=0.02, main=outName, cex.main=1)
            plotTracks(c(covTrack, ht, sTrack, gtrack), chromosome = chrom, from = chromStart, to = chromEnd, cex=0.5, fontsize=10, extend.left=0.15, extend.right=0.02, main=outName, cex.main=1, sizes=c(0.05,0.1,0.35,0.1,0.35,0.025,0.025))
        }else{
            plotTracks(c(covTrack, ht, sTrack, gtrack), chromosome = chrom, from = chromStart, to = chromEnd, cex=0.5, fontsize=10, extend.left=0.15, extend.right=0.02, main=outName, cex.main=1, sizes=c(0.05,0.1,0.35,0.1,0.35,0.025,0.025))
        }
        dev.off()
    }
'''

plotAsCountStatisticsStr = '''
    library(ggplot2)
    plotAsCountStatistics <- function(AsAnnotationFile, outPdf){
        as_anno_data <- read.csv(AsAnnotationFile, sep="\t")
        as_anno_data$AS_type <- factor(as_anno_data$AS_type, levels=c("SE", "A5SS", "A3SS", "IR", "APA", "PA"))
        pdf(outPdf)
        p <- ggplot(as_anno_data, aes(x=Annotation, y=Count, fill=Annotation)) + 
            geom_bar(stat = 'identity') + 
            facet_grid(~ AS_type) + 
            coord_cartesian(xlim=c(0.8,2.2)) + 
            theme_bw() + ggtitle("AS Type Statistics") + xlab("AS types") + 
            scale_y_continuous(expand = c(0.02, 0), trans="log10") + 
            theme(plot.title = element_text(hjust = 0.5, size=20), panel.border = element_blank(), 
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"), 
                  strip.background = element_blank(), legend.title = element_blank(), 
                  legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"), 
                  panel.spacing = unit(0, "lines"))
        print(p)
        dev.off()
    }
'''

plotAsDinucleotideStatisticsStr = '''
    library(ggplot2)
    library(dplyr)
    plotAsDinucleotideStatistics <- function(AsDinucleotideFile, outPdf){
        as_splice_data <- read.csv(AsDinucleotideFile, sep="\t")
        as_splice_data[is.na(as_splice_data)] <- 0
        as_splice_data$AS_type <- factor(as_splice_data$AS_type, levels=c("SE", "A5SS", "A3SS", "IR"))
        as_splice_data$Category <- factor(as_splice_data$Category, levels=c("Inc", "Exc"))
        pdf(outPdf)
        p <- ggplot(as_splice_data, aes(x=Category, y=Frequency, fill=Dinucleotide)) + 
            geom_bar(stat = 'identity') + 
            facet_grid(~ AS_type) + 
            coord_cartesian(xlim=c(0.8,2.2)) + 
            theme_bw() + ggtitle("AS Dinucleotide Statistics") + xlab("AS types") + 
            scale_y_continuous(expand = c(0.02, 0)) + 
            theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                  strip.background = element_blank(), legend.title = element_blank(), 
                  legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"), 
                  panel.spacing = unit(0, "lines"))
        print(p)
        dev.off()
    }
    
'''

plotNovelHqASStr = '''
    library(ggplot2)
    plotNovelHqAS <- function(isoformScoreFile, outPdf){
        isoformScore <- read.csv(isoformScoreFile, sep="\t", header=FALSE)
        names(isoformScore) <- c("Gene", "Isos", "Count", "Total_count", "Freq", "Annotation")
        isoformScore <- isoformScore[order(-isoformScore$Freq),]
        isoformScoreNovel <- isoformScore[isoformScore$Annotation == "Novel", ]
        isoformScoreAnno <- isoformScore[isoformScore$Annotation == "Annotated", ]
        
        isoformScoreNovel$Isos <- factor(isoformScoreNovel$Isos, levels=as.vector(isoformScoreNovel$Isos))
        isoformScoreAnno$Isos <- factor(isoformScoreAnno$Isos, levels=as.vector(isoformScoreAnno$Isos))
        pdf(outPdf)
        p1 <- ggplot(data=isoformScoreNovel, aes(x=Isos, y=Freq, group=1)) + 
            geom_point() + 
            geom_line() + 
            theme_bw() + ggtitle("Annotated Isoforms Rank") + 
            theme(plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_blank(), 
                axis.title.x=element_blank(), text = element_text(size=12), panel.border = element_blank(), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), strip.background = element_blank(), 
                strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8), 
                legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines"))
        print(p1)
        p2 <- ggplot(data=isoformScoreAnno, aes(x=Isos, y=Freq, group=1)) + 
            geom_point() +
            geom_line() +
            theme_bw() + ggtitle("Novel Isoforms Rank") + 
            theme(plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_blank(), 
                axis.title.x=element_blank(), text = element_text(size=12), panel.border = element_blank(), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), strip.background = element_blank(), 
                strip.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8), 
                legend.key.size = unit(0.5, "cm"), panel.spacing = unit(0, "lines"))
        print(p2)
        dev.off()
    }
'''

plotDiffASStr = '''
    library(ggplot2)
    plotDiffAS <- function(diffAsFile, outPdf){
        diffAS <- read.csv(diffAsFile, sep="\t")
        pdf(outPdf)
        p <- ggplot(diffAS, aes(x=AS_type, y=Count, fill=AS_type)) + geom_bar(stat = 'identity') + 
            theme_bw() + ggtitle("Diff AS Type Statistics") + xlab("AS types") + 
            scale_y_continuous(expand = c(0.02, 0)) + 
            theme(plot.title = element_text(hjust = 0.5, size=20),
                panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_blank(), 
                legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"), 
                panel.spacing = unit(0, "lines"))
        print(p)
        dev.off()
    }
    
'''

plotTargetGenesGoEnrichmentStr = '''
    library(GO.db)
    library(stringr)
    library(plyr)
    library(ggplot2)
    library(AnnotationDbi)
    library(clusterProfiler)
    
    plotTargetGenesGoEnrichment <- function(targetGeneFiles, sampleNames, gene2goFile, outName){
        targetGeneFiles <- unlist(strsplit(targetGeneFiles, ","))
        sampleNames <- unlist(strsplit(sampleNames, ","))
        gene2go <- read.delim(gene2goFile, header=T, sep="\t")
        go2term <- unique(toTable(GOTERM)[,c(1,3)])
        go2gene <- gene2go[,c(2,1)]
        names(targetGeneFiles) <- sampleNames
        outPdf <- paste0(outName, ".pdf")
        if(length(sampleNames)>1){
            r_bind <- data.frame()
            for (i in seq(length(targetGeneFiles))){
                sample <- names(targetGeneFiles[i])
                genes <- read.csv(targetGeneFiles[i], header=FALSE)
                genes <- as.vector(genes[,1])
                goRes <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff=1, maxGSSize = 100000)
                goResult <- summary(goRes)
                # goResult <- goResult[order(goResult$p.adjust),]
                goResult <- goResult[which(goResult$p.adjust<=0.05),]
                goResult$Ontology <- Ontology(as.vector(goResult$ID))
                goResult <- goResult[c(1, ncol(goResult), 2:(ncol(goResult)-1))]
                goResult$Ontology <- revalue(goResult$Ontology, c("BP"="Biological process", "MF"="Molecular function", "CC"="Cellular component"))
                goResult <- goResult[order(goResult$Ontology, goResult$pvalue),]
                goResult$Description <- factor(goResult$Description, levels = as.vector(goResult$Description))
                outFile <- paste0(sample, ".goEnrichResults.txt")
                write.table(goResult, file=outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
                goResult$Cluster <- sample
                goResult$group <- sample
                goResult$level <- sample
                r_bind <- rbind(r_bind, goResult)
            }
    
            r_bind <- na.omit(r_bind)
        
            r_bind$level <- factor(r_bind$level, levels=sampleNames)
        
            r_bind_new <- new("compareClusterResult", compareClusterResult = r_bind)
            pdf(outPdf, width=10, height=8)
            p <- dotplot(r_bind_new, showCategory=30, x=~group) + scale_color_continuous(low='purple', high='green') + 
                scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + 
                scale_size(range=c(0, 5)) + ggplot2::facet_grid(~level)
            print(p)
            dev.off()
        }else{
            sample <- sampleNames
            genes <- read.csv(targetGeneFiles[[1]][1], header=FALSE)
            genes <- as.vector(genes[,1])
            goRes <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff=0.05, maxGSSize = 100000)
            goResult <- summary(goRes)
            goResult <- goResult[which(goResult$p.adjust<=0.05),]
            goResult$Ontology <- Ontology(as.vector(goResult$ID))
            goResult <- goResult[c(1, ncol(goResult), 2:(ncol(goResult)-1))]
            goResult$Ontology <- revalue(goResult$Ontology, c("BP"="Biological process", "MF"="Molecular function", "CC"="Cellular component"))
            goResult <- goResult[order(goResult$Ontology, goResult$pvalue),]
            goResult$Description <- factor(goResult$Description, levels = as.vector(goResult$Description))
            outFile <- paste0(sample, ".goEnrichResults.txt")
            write.table(goResult, file=outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
            pdf(outPdf, width=10, height=8)
            p <- dotplot(goRes, showCategory=30) + scale_color_continuous(low='purple', high='green') + 
                scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size(range=c(0, 5))
            print(p)
            dev.off()
        }
        
    }
    
'''

plotPaTailASStr = '''
library(RColorBrewer)
library(grid)
library(gridBase)
library(ggplot2)
library(valr)

convertstrandinfo <- function(strandvector)
{
  if (class((strandvector)) == "character")
  {
    b = rep(1,length((strandvector)))
    b[which((strandvector) == "-")] = -1
    strandvector = b
  }
  return (strandvector)
}

plotGenes <- function(geneinfo=NULL, chrom=NULL, chromstart=NULL,chromend=NULL,
                      col="red",bheight=0.3,lheight=0.3,bentline=TRUE,
                      packrow=TRUE,maxrows=10000, colorbygroup=NULL,
                      types="exon",wigglefactor=0.05,
                      labeltext=TRUE,labeloffset=0.4,fontsize=.7,fonttype=2,labelat="middle",...)
{

  # convert strand info
  if (ncol(geneinfo) >= 6)
  {
    geneinfo[,6] = convertstrandinfo(geneinfo[,6])
  }

  # Define a function that merges overlapping regions
  mergetypes <- function(exons)
  {
    newexons = c()
    for (subtype in names(table(exons[,8])))
    {
      # get subtype
      subexons = exons[which(exons[,8] == subtype),]
      subexons = subexons[order(subexons[,2]),]

      # skip if none
      if (nrow(subexons) == 0)
      {
        next
      }

      # check for overlap
      i = 0
      for (j in (1:nrow(subexons)))
      {
        i = i + 1
        curstart = subexons[j,2]
        curstop = subexons[j,3]
        if (i > 1)
        {

          if (subexons[j,2] <= curstop)
          {
            subexons[j,2] = curstart
            subexons[j,3] = max(subexons[j,3],curstop)
          }
        }
      }

      # make sure all starts and stops are the same
      for (startpos in names(table(subexons[,2])))
      {
        subexons[which(subexons[2] == startpos),3] = max(subexons[which(subexons[2] == startpos),3])
      }

      # remove duplicate rows
      subexons = subexons[!duplicated(subexons),]
      newexons = rbind(newexons,subexons)
    }
    return (newexons)
  }


  # Define a function that plots the gene structure given a y-value and all exons
  plottranscript <- function(exons,col,yvalue,bheight,lheight,bentline=TRUE,border="black",
                             bprange,labeltext=TRUE,labeloffset=0.4,
                             fontsize=.7,fonttype=1,labelat="middle",...)
  {
    strand = exons[1,6]

    # if label is true add the label

   labellocation = mean(c(max(exons[,c(2,3)]),min(exons[,c(2,3)])))
   if (labeltext == TRUE)
   {
       labellocation =  min(exons[,c(2,3)])

      adj = 0.0
      if (strand == 1)
      {
        adj = 0.0
      }
      text(labellocation,yvalue+labeloffset,labels=exons[1,4],adj=adj,cex=fontsize,font=fonttype)
   }

    # make sure coordinates are in the correct order
    min = apply(exons[,c(2,3)],1,min)
    max = apply(exons[,c(2,3)],1,max)
    exons[,2] = min
    exons[,3] = max

    # merge types
    exons = mergetypes(exons)

    # sort by order
    exons = exons[order(exons[,2],exons[,3]),]

    # combine all types of exons for proper line plotting
    allexons =  exons
    allexons$types = "exon"
    allexons = mergetypes(allexons)

    if (nrow(allexons) > 1)
    {
      # organize line info
      linecoords = cbind(allexons[(1:nrow(allexons)-1),3],allexons[(2:nrow(allexons)),2])
      linecoords = cbind(linecoords,apply(linecoords,1,mean))

      # plot lines
      top = yvalue
      if (bentline == TRUE)
      {
        top = yvalue + lheight
      }

      for (i in (1:nrow(linecoords)))
      {
        if ((linecoords[i,2] - linecoords[i,1]) > 0 )
        {
          segments(x0 = linecoords[i,1],
                   x1 = linecoords[i,3],
                   y0 = yvalue,
                   y1 = top,
                   col     = "black",
                   ...)
          segments(x0 = linecoords[i,3],
                   x1 = linecoords[i,2],
                   y0 = top,
                   y1 = yvalue,
                   col     = "black",
                   ...)
          if (strand == -1){
              if (i == 1){
                  arrows(x0 = linecoords[i,1],
                  x1 = linecoords[i,3],
                  y0 = yvalue,
                  y1 = top,
                  col     = "black", angle=30, length=0.1, code=1)
              }
          }else{
              if (i == nrow(linecoords)){
                  arrows(x0 = linecoords[i,3],
                  x1 = linecoords[i,2],
                  y0 = top,
                  y1 = yvalue,
                  col     = "black", angle=30, length=0.1)
              }
            }
        }
      }
    }

    # add boxes
    for (i in (1:nrow(exons)))
    {
      height = bheight
      if (exons[i,7] == "utr")
      {
        height = bheight/2
      }

      # boxes
        rect(xleft    = exons[i,2],
        ybottom  = yvalue - height,
        xright   = exons[i,3],
        ytop    = yvalue + height,
        col     = col,
        border  = col,
        ...)
    }
  }



  # Define a function th determines which row to plot a gene on
  checkrow <- function(data,alldata,maxrows,wiggle=0)
  {

    startcol = 2
    stopcol  = 3


    # strand = data[1,5]

    for (row in (1:maxrows))
    {

        thestart = as.numeric(data[startcol]) - wiggle
        thestop  = as.numeric(data[stopcol])  + wiggle
      #   if (plotgenetype == "box")
      #   {
      # }

      if (nrow(alldata[which(alldata$plotrow == row &
                               ((thestart >= alldata[,startcol] & thestart <= alldata[,stopcol]) |
                                  (thestop >= alldata[,startcol] & thestop <= alldata[,stopcol]) |
                                  (thestart <= alldata[,startcol] & thestop >= alldata[,stopcol]))),]) == 0)
      {
        return (row)
      }
    }
    return (NA)
  }

  # remove unwanted columns
    geneinfo = geneinfo[,seq(1,7)]
    colnames(geneinfo) = c("chrom", "start", "end", "transcript", "score", "strand", "group")

  # establish start and stop columns
  startcol = 2
  stopcol  = 3

  # add types
  geneinfo$types = types

  # remove lines with NA
  geneinfo = geneinfo[which(is.na(geneinfo[,2]) ==FALSE),]

  # color by group
  if (is.null(colorbygroup) == FALSE){
      geneinfo$group = as.factor(geneinfo$group)
      if (length(levels(geneinfo$group)) < 3){
          colors = brewer.pal(3, "Set1")
      }else{
          colors = brewer.pal(length(levels(geneinfo$group)), "Set1")
      }
    geneinfo$colors = colors[geneinfo$group]
  }

  if (is.null(colorbygroup) == TRUE)
  {
    geneinfo$colors = col
  }

  # set xlim
  if (is.null(chrom) == TRUE | is.null(chromstart) == TRUE | is.null(chromend) == TRUE)
  {
    chrom      = geneinfo[1,1]
    chromstart = min(geneinfo[,c(startcol,stopcol)])
    chromend   = max(geneinfo[,c(startcol,stopcol)])
    extend     = abs(chromend - chromstart) * 0.04
    chromstart = chromstart - extend
    chromend   = chromend   + extend
  }
  else{
      extend     = abs(chromend - chromstart) * 0.04
      chromstart = chromstart - extend
      chromend   = chromend + extend
  }

  # get bprange
  bprange = abs(chromend - chromstart)

  # set wiggle
  wiggle = bprange * wigglefactor

  # get number of geneinfo
  numberofgeneinfo = length(names(table(geneinfo[,4])))
  namesgeneinfo    = names(table(geneinfo[,4]))

  # sort by length
  starts = c()
  stops  = c()
  sizes  = c()
  strands = c()


  # collect the info for each transcript
  for (i in (1:numberofgeneinfo))
  {
    subgeneinfo  = geneinfo[which(geneinfo[,4] == namesgeneinfo[i]),]
    starts = c(starts,min(subgeneinfo[,2:3]))
    stops  = c(stops, max(subgeneinfo[,2:3]))
    sizes  = c(sizes, stops[i] - starts[i])
    strands = c(strands,subgeneinfo[1,6])
  }

  transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
  transcriptinfo = transcriptinfo[order(sizes,decreasing=TRUE),]



  # get row information
  if (packrow == TRUE)
  {
    transcriptinfo$plotrow = 0
    for (i in (1:nrow(transcriptinfo)))
    {
      transcriptinfo$plotrow[i] = checkrow(transcriptinfo[i,],transcriptinfo,maxrows=maxrows,wiggle=wiggle)
    }
  }
    rowNum <- nrow(transcriptinfo)
  if (packrow == FALSE)
  {
    transcriptinfo$plotrow = seq(1:rowNum)
  }
    rowHeight = 0.5
    if (is.null(rowHeight) == FALSE){
        medianRow <- which(transcriptinfo$plotrow == quantile(transcriptinfo$plotrow, .5, type = 1))
        newMaxHeight <- transcriptinfo$plotrow[medianRow] + (rowNum - medianRow) * rowHeight
        newMinHeight <- transcriptinfo$plotrow[medianRow] - (medianRow - 1) * rowHeight
        transcriptinfo$plotrow <- seq(newMinHeight, newMaxHeight, by=rowHeight)
    }



  # make the empty plot
  offsettop = 0.5

  # filter out rows above max row
  transcriptinfo = transcriptinfo[which(is.na(transcriptinfo$plotrow)==FALSE),]

  # filter out transcrits that don't overlap region
  transcriptinfo = transcriptinfo[which((transcriptinfo[,2] > chromstart & transcriptinfo[,2] < chromend)
                       | (transcriptinfo[,3] > chromstart & transcriptinfo[,3] < chromend)),]

  if (nrow(transcriptinfo) == 0)
  {
    toprow = 1
  }
  if (nrow(transcriptinfo) > 0)
  {
    toprow = max(transcriptinfo$plotrow)
  }

  plot(c(1,1),xlim=c(chromstart,chromend),ylim=c(0.5,(toprow  + offsettop)),type ='n',bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i")

  if (nrow(transcriptinfo) > 0)
  {
    for (i in (1:nrow(transcriptinfo)))
    {
      subgeneinfo  = geneinfo[which(geneinfo[,4] == transcriptinfo[i,1]),]
      plottranscript(subgeneinfo,col=subgeneinfo$colors[1],yvalue=transcriptinfo$plotrow[i],bheight=bheight,lheight=lheight,bentline=bentline,border=col,
                     bprange=bprange,
                     labeltext=labeltext,labeloffset=labeloffset,fontsize=fontsize,fonttype=fonttype,labelat=labelat)
    }
  }

}

plotPaTailAsStructure <- function(sigFile, outPdf){
    inFile <- sigFile
    data <- read.csv(inFile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    data <- data[!duplicated(data), ]
    data[, 14] <- as.factor(data[, 14])
    pdf(outPdf)
    for (i in levels(data[, 14])){
        subData <- data[data[, 14] == i, ]
        groupInfo <- subData[, 15]
        gene <- unique(subData[, 13])
        names(groupInfo) <- subData[, 4]
        bed12Data <- subData[, 1:12]
        names(bed12Data) <- c("chrom", "start", "end", "name", "score", "strand", "cds_start", "cds_end", "item_rgb", "exon_count", "exon_sizes", "exon_starts")
        bed6Data <- bed12_to_exons(bed12Data)
        bed6Data <- as.data.frame(bed6Data)
        bed6Data <- bed6Data[order(bed6Data$name), ]
        bed6Data$group <- ""
        for (j in names(groupInfo)){
            bed6Data[which(bed6Data$name==j), ]$group <- groupInfo[[j]]
        }
        paLenData <- subData[, c(4,15,16,17)]
        names(paLenData) <- c("isoform", "group", "paLen", "count")
        s <- strsplit(paLenData$paLen, split = ",")
        paLenData <- data.frame(isoform = rep(paLenData$isoform, sapply(s, length)), group = rep(paLenData$group, sapply(s, length)), paLen = as.integer(unlist(s)), count = rep(paLenData$count, sapply(s, length)), stringsAsFactors = FALSE)
    
        chrom <- bed6Data[1,1]
        chromStart <- min(bed6Data$start, bed6Data$end)
        chromEnd <- max(bed6Data$start, bed6Data$end)
    
        # sigPic <- gsub("/", "_", i)
        # sigPic <- paste0(file.path(outdir, sigPic), ".pdf")
    
        par(mfrow = c(2,1), mar=c(0.5,0.5,0,0), oma=c(0,0,0,0))
        plot.new()
    
        colorbygroup <- TRUE
        plotGenes(bed6Data, chrom, chromStart, chromEnd, maxrows=50,bentline=FALSE,col="red", labeloffset=0.2,fontsize=1, bheight = 0.08, lheight = 0.4, packrow = FALSE, labelat = "start", colorbygroup=colorbygroup)
    
        if (is.null(colorbygroup) == FALSE){
            paLenData$isoform = as.factor(paLenData$isoform)
            if (length(levels(paLenData$isoform)) < 3){
              colors = brewer.pal(3, "Set1")
            }else{
              colors = brewer.pal(length(levels(paLenData$isoform)), "Set1")
            }
        }
        labelColor <- colors[paLenData$isoform]
        paLenData$color <- as.factor(labelColor)
        # labelColor <- colors[0:length(levels(paLenData$isoform))]
        vp <- viewport(height = unit(0.5,"npc"), width=unit(0.95, "npc"),  just = c("left","top"), y = 1, x = 0)
    
        maxPaLen <- max(paLenData$paLen)
        p <- ggplot(paLenData, aes(x=group, y=paLen, label=paste0("n = ", count), fill=group)) +
            geom_violin() + geom_boxplot(width=0.1) + ggtitle(gene) + ylab("Poly(A) tail length") +
            geom_text(y=maxPaLen+1, vjust=0, size=4, fontface="plain", color=labelColor) +
            coord_cartesian(ylim=c(0,maxPaLen), clip="off") + theme_bw() +
            theme(plot.margin=unit(c(2,0,0,1), "lines"), text = element_text(size=12, face="bold"),
            plot.title = element_text(hjust = 0.5, vjust = 5), axis.title.x=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
            scale_fill_manual(values = levels(paLenData$color))
        print(p, vp = vp)
    
    }
    dev.off()
}


'''
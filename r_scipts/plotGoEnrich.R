#!/usr/bin/env Rscript

args <- commandArgs()
scriptPath <- strsplit(args[4], '=', fixed = T)[[1]][2]
scriptName <- basename(scriptPath)
args <- args[-(1:5)]
usage = function(){
    cat(paste0("Usage: ", scriptName) )
    cat(" -g=targetGeneFiles -s=sampleNames -bg=gene2go -o=outfile
Option:
    Commmon:
    -g|targetGeneFiles  str     The path to file contains the genes you want to enrich, or the file list separated by comma.
    -s|sampleNames      str     The sample name you want to specify whose order should consistent with the gene file(s).
    -bg|gene2go         str     The gene2go file which used server as background in the enrichment.
    -of|outfile         str     The output file name the plot will be stored.
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

targetGeneFiles <- ""
sampleNames <- ""
gene2goFile <- ""
outfile <- ""

if(length(args) >= 1){
    for(i in 1:length(args)){
        arg <- args[i]
        if(arg == '-h' || arg == '-help') usage()
        tmp <- parseArg(arg, 'g(enes)', 'targetGeneFiles'); if(!is.null(tmp)) targetGeneFiles <- tmp
        tmp <- parseArg(arg, 's(amples)', 'sampleNames'); if(!is.null(tmp)) sampleNames <- tmp
        tmp <- parseArg(arg, 'bg', 'gene2go'); if(!is.null(tmp)) gene2goFile <- tmp
        tmp <- parseArg(arg, 'of', 'outfile'); if(!is.null(tmp)) outfile <- tmp
    }
}

library(GO.db)
library(stringr)
library(plyr)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)

targetGeneFiles <- strsplit(targetGeneFiles, ",")
sampleNames <- unlist(strsplit(sampleNames, ","))
gene2go <- read.delim(gene2goFile, header=T, sep="\t")
go2term <- unique(toTable(GOTERM)[,c(1,3)])
go2gene <- gene2go[,c(2,1)]

names(targetGeneFiles) <- sampleNames
r_bind <- data.frame()
for (i in seq(length(targetGeneFiles))){
    sample <- names(targetGeneFiles[i])
    genes <- read.csv(targetGeneFiles[i][[1]], header=FALSE)
    goRes <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff=1, maxGSSize = 100000)
    goResult <- summary(goRes)
    # goResult <- goResult[order(goResult$p.adjust),]
    goResult <- goResult[which(goResult$p.adjust<=0.05),]
    goResult$Ontology <- Ontology(as.vector(goResult$ID))
    goResult <- goResult[c(1, ncol(goResult), 2:(ncol(goResult)-1))]
    goResult$Ontology <- revalue(goResult$Ontology, c("BP"="Biological process", "MF"="Molecular function", "CC"="Cellular component"))
    goResult <- goResult[order(goResult$Ontology, goResult$pvalue),]
    goResult$Description <- factor(goResult$Description, levels = as.vector(goResult$Description))
    goResult$Cluster <- sample
    goResult$group <- sample
    goResult$level <- sample
    r_bind <- rbind(r_bind, goResult)
}

r_bind <- na.omit(r_bind)

r_bind$level <- factor(r_bind$level, levels=sampleNames)

r_bind_new <- new("compareClusterResult", compareClusterResult = r_bind)
pdf(outfile, width=10, height=8)
p <- dotplot(r_bind_new, showCategory=30, x=~group) + scale_color_continuous(low='purple', high='green') + scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size(range=c(0, 5)) + ggplot2::facet_grid(~level)
print(p)
dev.off()
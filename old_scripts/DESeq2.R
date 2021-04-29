#!/bin/env Rscript
library(DESeq2)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

#prepare data
expFile <- args[1]
newColNames <- args[2]
selectColNames <- args[3]
conditions <- args[4]
expData <- read.table(expFile, row=1, header=TRUE, quote="\t", sep="\t", skip=1)
expData <- expData[7:ncol(expData)]

newColNames <- unlist(strsplit(newColNames, ","))
selectColNames <- unlist(strsplit(selectColNames, ","))
conditions <- unlist(strsplit(conditions, ","))
condition1 <- unique(conditions)[1]
condition2 <- unique(conditions)[2]
names(expData) <- newColNames

subexpData <- subset(expData, select=selectColNames)
subSampleCond <- data.frame(condtions = conditions, row.names = selectColNames)
row_sub <- !apply(subexpData, 1, function(row) all(row==0))
subexpData <- subexpData[row_sub, ]

CountTable <- DESeqDataSetFromMatrix(countData = subexpData, colData = subSampleCond, design = ~conditions)
dds <- DESeq(CountTable)
contrastV <- c("conditions", condition1, condition2)
res <- results(dds, contrast=contrastV)

baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == condition1]
if(is.vector(baseA)){
	baseMeanA <- as.data.frame(baseA)
}else{
	baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- condition1

baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == condition2]
if(is.vector(baseB)){
	baseMeanB <- as.data.frame(baseB)
}else{
	baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- condition2
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- cbind(gene_id=rownames(res), res)
res$padj[is.na(res$padj)] <- 1
res <- res[order(res$padj),]
write.table(res, file=paste(condition1, "vs", condition2, ".DESeq2.result.xls", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

#up-regulated and down-regulated gene output
selectCol <- c('gene_id', condition1, condition2, 'log2FoldChange','padj')
res_de <- subset(res, res$padj<0.1, select=selectCol)
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
write.table(as.data.frame(res_de_up), file=paste(condition1, 'vs', condition2, '.DESeq2.up-regulated.result.xls',sep=''), sep="\t", quote=FALSE, row.names=FALSE)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)
write.table(as.data.frame(res_de_dw), file=paste(condition1, 'vs', condition2, '.DESeq2.down-regulated.result.xls',sep=''), sep="\t", quote=FALSE, row.names=FALSE)
de_gene <- c(as.vector(res_de_up$gene_id), as.vector(res_de_dw$gene_id))
write.table(as.data.frame(de_gene), file=paste(condition1, '_', condition2, '.DESeq2.genes.txt', sep=''), sep="\t", quote=F, row.names=F, col.names=F)

#volcano plot
vp <- res[c('log2FoldChange','padj')]
vp$FDR <- -log10(vp$padj)
vp$change <- as.factor(ifelse(vp$padj<0.01 & abs(vp$log2FoldChange)>2, ifelse(vp$log2FoldChange>2, "Up", "Down"), "NoDiff"))
ggplot(vp,aes(x=log2FoldChange, y=FDR))+
  geom_point(aes(color=change))+
  geom_hline(yintercept=-log10(0.01), linetype=2)+
  geom_vline(xintercept=c(-2,2), linetype=2)+
  scale_colour_manual(values = c('red', 'green', 'grey'), limits= c("Up", "Down", "NoDiff"))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(paste(condition1, 'vs', condition2, 'volcano plot', sep=' '))+
  ylab('-log10(padj)')+
  xlab('log2(FC)')+
  theme(plot.title = element_text(hjust=0.5,vjust=2,lineheight=.8, face="bold",size=16,margin = margin(b=20,t=10)),
        axis.text = element_text(size=14,color='black'),
        axis.title=element_text(size=16),
        axis.title.x = element_text(margin = margin(t=20,b=10)),
        axis.title.y = element_text(margin = margin(r=20)),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        aspect.ratio=1)
ggsave(paste(condition1, 'vs', condition2, '.volcano.pdf', sep=''))

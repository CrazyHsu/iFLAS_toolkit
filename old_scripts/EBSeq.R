#!/bin/env Rscript
library(EBSeq)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

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
# subSampleCond <- data.frame(condtions = conditions, row.names = selectColNames)
row_sub <- !apply(subexpData, 1, function(row) all(row==0))
subexpData <- subexpData[row_sub, ]

Sizes.norep <- MedianNorm(subexpData)

#data(GeneMat)
#GeneMat.norep=GeneMat[,c(1,6)]
#Sizes.norep=MedianNorm(GeneMat.norep)
EBOut.norep=EBTest(Data=subexpData, Conditions=as.factor(rep(c("C1","C2"))), sizeFactors=Sizes.norep, maxround=5)
EBDERes.norep=GetDEResults(EBOut.norep)
write.table(EBDERes.norep$DEfound, file=paste0(condition1, "vs", condition2, "_DEG_FDR005.txt"), quote=F,sep="\t",row.names=F,col.names=F)
GeneFC.norep=PostFC(EBOut.norep)

C1_Mean <- unlist(EBOut.norep$C1Mean)
C2_Mean <- unlist(EBOut.norep$C2Mean)
FoldChange <-  GeneFC.norep$PostFC
log2FoldChange <- log(FoldChange, 2)
FDR <- EBDERes.norep$PPMat[, 1]

FDR <- FDR[which(names(FDR) %in% names(C1_Mean))]
EBSeq_stat <- data.frame(C1_Mean, C2_Mean, log2FoldChange, FDR)
write.table(EBSeq_stat, file=paste0(condition1, "vs", condition2, "_results.txt"), sep = '\t', col.names = NA, quote = FALSE)
EBSeq_stat$change <- as.factor(ifelse(EBSeq_stat$FDR<0.05 & abs(EBSeq_stat$log2FoldChange)>2, ifelse(EBSeq_stat$log2FoldChange>2, "Up", "Down"), "NoDiff"))

EBSeq_stat$FDR <- -log10(EBSeq_stat$FDR)
ggplot(EBSeq_stat,aes(x=log2FoldChange, y=FDR))+
  geom_point(aes(color=change))+
  geom_hline(yintercept=-log10(0.05), linetype=2)+
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

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: isoViewerFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 19:49:40
Last modified: 2019-12-08 19:49:41
'''
import os, subprocess
from commonFuncs import *
from commonObjs import *
import cPickle as pickle
import shutil

def makeIsoPreCorr(testGene, targetGeneBed, projectName=None, sampleName=None, dirSpec=None):
    targetGenePreCorrBedFile = "preCorr.%s.bed" % (testGene)
    targetGenePreCorrGpeFile = "preCorr.%s.gpe" % (testGene)
    targetGenePreCorrGtfFile = "preCorr.%s.gtf" % (testGene)
    targetGenePreCorrGffFile = "preCorr.%s.gff" % (testGene)
    flncAlnSortedBed = os.path.join(dirSpec.out_dir, projectName, sampleName, "isoseq3", "flnc.mm2.sorted.bed12")
    cmd = '''echo -e "%s" |bedtools intersect -a %s -b - -wa -s | awk -v targetGene=%s '{print $0"\t"targetGene}' > %s''' % (
    targetGeneBed, flncAlnSortedBed, testGene, targetGenePreCorrBedFile)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "bed2gpe.pl -g 13 %s > %s" % (targetGenePreCorrBedFile, targetGenePreCorrGpeFile)
    subprocess.call(cmd, shell=True)
    # targetGeneMinpos, targetGeneMaxpos = int(targetGeneBed.split("\t")[1]), int(targetGeneBed.split("\t")[2])
    # targetGpeObj = GenePredObj(targetGenePreCorrGpeFile)
    # redefineBorderGPE = targetGpeObj.redefineBorder(targetGeneMinpos, targetGeneMaxpos, outFile="tmp.gpe")
    # targetGenePreCorrGpeFile = redefineBorderGPE
    itemCounts = getFileRowCounts(targetGenePreCorrGpeFile)
    relativeSize = resizeTrackRatio(itemCounts)
    cmd = "genePredToGtf file %s %s -source=iFLAS" % (targetGenePreCorrGpeFile, targetGenePreCorrGtfFile)
    subprocess.call(cmd, shell=True)
    cmd = "gene_model_to_splicegraph.py -m %s -o %s" % (targetGenePreCorrGtfFile, targetGenePreCorrGffFile)
    subprocess.call(cmd, shell=True)
    return targetGenePreCorrGffFile, itemCounts, relativeSize

def makeNGSIso(testGene, targetGeneBed, projectName=None, sampleName=None, dirSpec=None):
    targetGeneNGSGtfFile = "ngsAssemble.%s.gtf" % (testGene)
    targetGeneNGSGffFile = "ngsAssemble.%s.gff" % (testGene)
    targetGeneNGSGpeFile = "ngsAssemble.%s.gpe" % (testGene)
    stringtieMergeGTF = os.path.join(dirSpec.out_dir, projectName, sampleName, "RNA-seq", "reassembly", "stringtie_merged.gtf")
    cmd = "gtfToGenePred {} stringtie_merged.gp -genePredExt".format(stringtieMergeGTF)
    subprocess.call(cmd, shell=True)
    stringtieMergeBED = "stringtie_merged.bed"
    GenePredObj("stringtie_merged.gp", bincolumn=False).toBed(outFile=stringtieMergeBED)
    # cmd = "genePredToBed stringtie_merged.gp {}".format(stringtieMergeBED)
    # subprocess.call(cmd, shell=True)
    cmd = '''echo -e "%s" |bedtools intersect -a %s -b - -wa -s | awk -v targetGene=%s '{print $0"\t"targetGene}' > tmp.bed ''' % (targetGeneBed, stringtieMergeBED, testGene)
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "bed2gpe.pl -g 13 tmp.bed > {}".format(targetGeneNGSGpeFile)
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file %s %s -source=iFLAS" % (targetGeneNGSGpeFile, targetGeneNGSGtfFile)
    itemCounts = getFileRowCounts(targetGeneNGSGpeFile)
    relativeSize = resizeTrackRatio(itemCounts)
    subprocess.call(cmd, shell=True)
    cmd = "gene_model_to_splicegraph.py -m %s -o %s -a" % (targetGeneNGSGtfFile, targetGeneNGSGffFile)
    subprocess.call(cmd, shell=True)
    # os.remove("tmp.gpe")
    # os.remove("tmp.gtf")
    return targetGeneNGSGffFile, itemCounts, relativeSize

def parallelPlotterAnno(gene, gpeTargetGenePickle, sampleTargetGenePickle, tgsSample, dirSpec):
    isoViewDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "isoViewer")
    resolveDir(os.path.join(isoViewDir, gene))

    sampleTargetGeneObj = pickle.loads(sampleTargetGenePickle)
    targetGeneChrom = sampleTargetGeneObj.chrom
    targetGeneMinpos = sampleTargetGeneObj.minpos
    targetGeneMaxpos = sampleTargetGeneObj.maxpos
    targetGeneStrand = sampleTargetGeneObj.strand

    # Gene model section
    geneModelGPE = "%s.gpe" % (gene)
    geneModelGTF = "%s.gtf" % (gene)
    geneModelGFF = "%s.gff" % (gene)

    gpeTargetGene = pickle.loads(gpeTargetGenePickle)
    geneModelGPEOut = open(geneModelGPE, "w")
    gpeMinposList, gpeMaxposList = [targetGeneMinpos], [targetGeneMaxpos]
    for i in gpeTargetGene:
        gpeMinposList.append(i.txStart)
        gpeMaxposList.append(i.txEnd)
        print >> geneModelGPEOut, i
    geneModelGPEOut.close()
    plotMinpos = min(gpeMinposList)
    plotMaxpos = max(gpeMaxposList)
    targetGeneBed = "{}\t{}\t{}\t{}\t.\t{}".format(targetGeneChrom, plotMinpos, plotMaxpos, gene,
                                                   sampleTargetGeneObj.strand)
    targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
    os.system("genePredToGtf file %s %s -source=iFLAS" % (geneModelGPE, geneModelGTF))
    os.system("gene_model_to_splicegraph.py -m %s -o %s" % (geneModelGTF, geneModelGFF))
    geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                    gene_name=gene, relative_size=5.0,
                                    title_string="Gene Model for %gene", hide=True)
    geneModelItemCount = getFileRowCounts(geneModelGPE)
    geneModelRelativeSize = resizeTrackRatio(geneModelItemCount)
    geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                       source_file=geneModelGFF, relative_size=geneModelRelativeSize,
                                       title_string="Gene Model for %gene [{}({})]".format(
                                           targetGeneRegion,
                                           targetGeneStrand))

    # Corrected tgs reads from the pipeline
    postCorrIsoforms = sampleTargetGeneObj.reads
    postCorrIsoGPE = "%s.iFLAS.gpe" % gene
    postCorrIsoGTF = "%s.iFLAS.gtf" % gene
    postCorrIsoGFF = "%s.iFLAS.gff" % gene
    postCorrIsoGPEOut = open(postCorrIsoGPE, "w")
    for read in postCorrIsoforms:
        print >> postCorrIsoGPEOut, postCorrIsoforms[read].go_to_gpe()
    postCorrIsoGPEOut.close()
    os.system("genePredToGtf file %s %s -source=iFLAS" % (postCorrIsoGPE, postCorrIsoGTF))
    postCorrIsoItemCounts = len(postCorrIsoforms)
    postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
    os.system("gene_model_to_splicegraph.py -m %s -o %s -a" % (postCorrIsoGTF, postCorrIsoGFF))
    postCorrIsoPlotType = "splice_graph" if postCorrIsoItemCounts < 30 else "isoforms"
    postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse]", plot_type=postCorrIsoPlotType,
                                  source_file=postCorrIsoGFF, relative_size=postCorrIsoRelativeSize,
                                  title_string="Corrected isoforms and AS events in %s from TGS data" % gene)

    # tgs flnc reads before process of the pipeline
    print targetGeneBed
    preCorrIsoGFF, preCorrIsoItemCounts, preCorrIsoRelativeSize = makeIsoPreCorr(gene,
                                                                                 targetGeneBed,
                                                                                 projectName=tgsSample.projectName,
                                                                                 sampleName=tgsSample.sampleName,
                                                                                 dirSpec=dirSpec)
    preCorrIsoPlot = PlotSection(section_name="[PreIsoform]", plot_type="isoforms",
                                 relative_size=preCorrIsoRelativeSize,
                                 source_file=preCorrIsoGFF,
                                 title_string="All TGS isoforms before correction")

    # Isoforms assembled from hisat2 + stringtie
    ngsAssembleGFF, ngsAssembleItemCounts, ngsAssembleRelativeSize = makeNGSIso(gene, targetGeneBed,
                                                                                projectName=tgsSample.projectName,
                                                                                sampleName=tgsSample.sampleName,
                                                                                dirSpec=dirSpec)
    ngsAssemblePlot = PlotSection(section_name="[ngsAssembleIsoform]", plot_type="isoforms",
                                  relative_size=ngsAssembleRelativeSize,
                                  source_file=ngsAssembleGFF, title_string="NGS assembled isoforms")

    figOut = gene + ".pdf"
    cfgOut = open(gene + ".cfg", "w")
    majorItemCount = geneModelItemCount + postCorrIsoItemCounts + preCorrIsoItemCounts + ngsAssembleItemCounts
    figHeight = 10 if majorItemCount <= 50 else 15 if majorItemCount <= 150 else 20
    mainSec = MainSection(fout=figOut, height=figHeight)
    print >> cfgOut, mainSec.printStr()
    print >> cfgOut, geneModelHidePlot.printStr()
    print >> cfgOut, geneModelVisiblePlot.printStr()
    print >> cfgOut, postCorrIsoPlot.printStr()
    print >> cfgOut, preCorrIsoPlot.printStr()
    print >> cfgOut, ngsAssemblePlot.printStr()

    # Abundance from sam file in NGS pipeline
    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        ngsSams = []
        for n in range(len(tgsSample.ngsLeftReads.split(";"))):
            repeatName = "repeat" + str(n)
            bamFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq",
                                   "alignment/{}/{}.sorted.bam".format(repeatName, repeatName))
            targetSam = "{}.{}.sam".format(repeatName, gene)
            cmd = "samtools view %s %s > %s" % (bamFile, targetGeneRegion, targetSam)
            subprocess.call(cmd, shell=True)
            ngsSams.append(targetSam)
            readsPlot = PlotSection(section_name="[Reads_%s]" % (repeatName), plot_type="read_depth",
                                    source_file=targetSam, relative_size=5.0,
                                    title_string="%s Read Coverage in sample %s" % (gene, repeatName))
            print >> cfgOut, readsPlot.printStr()

    cfgOut.close()
    os.system("plotter.py %s.cfg" % gene)
    # removeFiles(ngsSams)

def parallelPlotterNovel(gene, sampleTargetGenePickle, tgsSample, dirSpec):
    isoViewDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "isoViewer")
    resolveDir(os.path.join(isoViewDir, gene))

    sampleTargetGeneObj = pickle.loads(sampleTargetGenePickle)
    targetGeneChrom = sampleTargetGeneObj.chrom
    targetGeneMinpos = sampleTargetGeneObj.minpos
    targetGeneMaxpos = sampleTargetGeneObj.maxpos
    targetGeneStrand = sampleTargetGeneObj.strand

    # Gene model section
    geneModelGPE = "%s.gpe" % (gene)
    geneModelGTF = "%s.gtf" % (gene)
    geneModelGFF = "%s.gff" % (gene)

    geneModelGPEOut = open(geneModelGPE, "w")
    geneMinposList, geneMaxposList = [targetGeneMinpos], [targetGeneMaxpos]
    for read in sampleTargetGeneObj.reads:
        geneMinposList.append(sampleTargetGeneObj.reads[read].start)
        geneMaxposList.append(sampleTargetGeneObj.reads[read].end)
        print >> geneModelGPEOut, sampleTargetGeneObj.reads[read].go_to_gpe()
    geneModelGPEOut.close()
    plotMinpos = min(geneMinposList)
    plotMaxpos = max(geneMaxposList)
    targetGeneBed = "{}\t{}\t{}\t{}\t.\t{}".format(targetGeneChrom, plotMinpos, plotMaxpos, gene,
                                                   sampleTargetGeneObj.strand)
    targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
    os.system("genePredToGtf file %s %s -source=iFLAS" % (geneModelGPE, geneModelGTF))
    os.system("gene_model_to_splicegraph.py -m %s -o %s" % (geneModelGTF, geneModelGFF))
    geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                    gene_name=gene, relative_size=5.0,
                                    title_string="Gene Model for %gene", hide=True)

    # Corrected tgs reads from the pipeline
    postCorrIsoGPE = "%s.iFLAS.gpe" % gene
    postCorrIsoGTF = "%s.iFLAS.gtf" % gene
    postCorrIsoGFF = "%s.iFLAS.gff" % gene
    shutil.copy(geneModelGPE, postCorrIsoGPE)
    shutil.copy(geneModelGTF, postCorrIsoGTF)

    postCorrIsoItemCounts = len(sampleTargetGeneObj.reads)
    postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
    os.system("gene_model_to_splicegraph.py -m %s -o %s -a" % (postCorrIsoGTF, postCorrIsoGFF))
    postCorrIsoPlotType = "splice_graph" if postCorrIsoItemCounts < 30 else "isoforms"
    postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse]", plot_type=postCorrIsoPlotType,
                                  source_file=postCorrIsoGFF, relative_size=postCorrIsoRelativeSize,
                                  title_string="Corrected isoforms and AS events in {} [{}({})] from TGS data".format(
                                      gene, targetGeneRegion, targetGeneStrand))

    # Uncorrected tgs reads before process of the pipeline
    preCorrIsoGFF, preCorrIsoItemCounts, preCorrIsoRelativeSize = makeIsoPreCorr(gene,
                                                                                 targetGeneBed,
                                                                                 projectName=tgsSample.projectName,
                                                                                 sampleName=tgsSample.sampleName,
                                                                                 dirSpec=dirSpec)
    preCorrIsoPlot = PlotSection(section_name="[PreIsoform]", plot_type="isoforms",
                                 relative_size=preCorrIsoRelativeSize,
                                 source_file=preCorrIsoGFF,
                                 title_string="All TGS isoforms before correction")

    # Isoforms assembled from hisat2 + stringtie
    ngsAssembleGFF, ngsAssembleItemCounts, ngsAssembleRelativeSize = makeNGSIso(gene, targetGeneBed,
                                                                                projectName=tgsSample.projectName,
                                                                                sampleName=tgsSample.sampleName,
                                                                                dirSpec=dirSpec)
    ngsAssemblePlot = PlotSection(section_name="[ngsAssembleIsoform]", plot_type="isoforms",
                                  relative_size=ngsAssembleRelativeSize,
                                  source_file=ngsAssembleGFF, title_string="NGS assembled isoforms")

    figOut = gene + ".pdf"
    cfgOut = open(gene + ".cfg", "w")
    majorItemCount = postCorrIsoItemCounts + preCorrIsoItemCounts + ngsAssembleItemCounts
    figHeight = 10 if majorItemCount <= 50 else 15 if majorItemCount <= 150 else 20
    mainSec = MainSection(fout=figOut, height=figHeight)
    print >> cfgOut, mainSec.printStr()
    print >> cfgOut, geneModelHidePlot.printStr()
    # print >> cfgOut, geneModelVisiblePlot.printStr()
    print >> cfgOut, postCorrIsoPlot.printStr()
    print >> cfgOut, preCorrIsoPlot.printStr()
    print >> cfgOut, ngsAssemblePlot.printStr()

    # Abundance from sam file in NGS pipeline
    if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
        ngsSams = []
        for n in range(len(tgsSample.ngsLeftReads.split(";"))):
            repeatName = "repeat" + str(n)
            bamFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "RNA-seq",
                                   "alignment/{}/{}.sorted.bam".format(repeatName, repeatName))
            targetSam = "{}.{}.sam".format(repeatName, gene)
            cmd = "samtools view %s %s > %s" % (bamFile, targetGeneRegion, targetSam)
            subprocess.call(cmd, shell=True)
            ngsSams.append(targetSam)
            readsPlot = PlotSection(section_name="[Reads_%s]" % (repeatName), plot_type="read_depth",
                                    source_file=targetSam, relative_size=5.0,
                                    title_string="%s Read Coverage in sample %s" % (gene, repeatName))
            print >> cfgOut, readsPlot.printStr()

    cfgOut.close()
    os.system("plotter.py %s.cfg" % gene)
    # removeFiles(ngsSams)
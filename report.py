#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: report.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:53
Last modified: 2021-04-29 16:20:53
'''
from commonFuncs import *
import pandas as pd
import matplotlib.pyplot as plt
import PyPDF2
import seaborn as sns

def readsCorrectedEval(dataObj=None):
    filtrationDir = os.path.join(dataObj.out_dir, dataObj.project_name, dataObj.sample_name, "filtration")
    rawMappedBed = os.path.join(filtrationDir, "raw.mapped.addCVandID.bed12+")
    correctMappedBed = os.path.join(filtrationDir, "mapped.addCVandID.bed12+")
    rawMapped = pd.read_csv(rawMappedBed, sep="\t", header=None)
    correctMapped = pd.read_csv(correctMappedBed, sep="\t", header=None)
    rawMapped["accuracy"] = rawMapped.iloc[:, 12] * rawMapped.iloc[:, 13]
    correctMapped["accuracy"] = correctMapped.iloc[:, 12] * correctMapped.iloc[:, 13]
    for x in [rawMapped, correctMapped]:
        sns.distplot(x.accuracy, kde=False, bins=range(0, 100, 10))
    plt.savefig("readsCorrectResult.pdf")
    # g.savefig("readsCorrectResult.pdf")

def gcAcrossRead(fxFile, outFile, interval=20):
    import math
    from Bio.SeqUtils import GC
    out = open(outFile, "w")
    fileType = validateFaAndFqFile(fxFile)
    for seq in SeqIO.parse(fxFile, fileType):
        if len(seq) < interval:
            chunkSize = 1
        else:
            chunkSize = int(math.ceil(len(seq)/interval))
        gcList = [seq.name]
        for i in xrange(0, len(seq), chunkSize):
            gcList.append(GC(seq[i:i+chunkSize].seq))
        print >>out, "\t".join(map(str, gcList))
    out.close()

def readsContentEval(dataObj=None, refParams=None):
    flncFx = dataObj.data_processed_location

    cmd = "seqkit fx2tab -n -g {} | cut -f 1,2 > GC_of_raw_flnc_reads.log".format(flncFx)
    subprocess.call(cmd, shell=True)
    cmd = '''cut -f2 GC_of_raw_flnc_reads.log | distrCurve.R -d -m='GC Content of raw flnc Reads' -x='Binned GC%' -y='Fraction of Reads' -v=50 -p=GC_of_raw_flnc_reads.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    gcAcrossRead(flncFx, "GC_across_raw_flnc_read.log")
    cmd = '''cut -f2- GC_across_raw_flnc_read.log | box.R -stack -nJ -ho=50 -m='GC Content across raw flnc Reads' -x=Interval -y=GC% -oS=0.5 -w=11 -p=GC_across_raw_flnc_read.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

    cmd = "seqkit fx2tab -n -l {} | cut -f 2 > readsLength.lst".format(flncFx)
    subprocess.call(cmd, shell=True)
    cmd = "gpe2bed {} | bedLength.pl | cut -f 13 > referenceGeneLength.lst".format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True)
    cmd = '''distrCurves.R -x1=0 -x2=10000 -d -x='Binned Length (limited in 0-10000)' -w=15 *.lst -b=150 -p=LengthDistribution.curve.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)
    cmd = '''boxes.R -ng -no *.lst -p=LengthDistribution.box.pdf 2>/dev/null'''
    subprocess.call(cmd, shell=True)

def report(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Generate plot report for project {} sample {}...".format(projectName, sampleName)
    if dataObj.use_fmlrc2:
        readsCorrectedEval(dataObj=dataObj)
    readsContentEval(dataObj=dataObj)
    # preprocess (readsCorrectResult.pdf)
    # mapping
    # filtration
    # collapse
    # identify_as (asType_anno_novel, splice_pattern)
    # visual_as (单独页面, all target gene visualization.pdf)
    # rank_as (novel_as_rank.pdf)
    # allelic_as (Gviz.pdf, 单独页面, 所有的as都存在一个pdf中)
    # palen_as (polyaTailLength.pdf, asTypeCountDistribute.pdf)
    # diff_as (as_distribute.pdf)
    # go (go_enrichment.pdf)
    print getCurrentTime() + " Generate plot report for project {} sample {} done!".format(projectName, sampleName)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: relationshipBetweenPAlenAndAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-14 18:07:51
Last modified: 2020-05-14 18:07:51
'''
import subprocess, argparse, itertools
import scipy.stats as stats
from commonObjs import *

def isoformWithDiffPalen(collapsedGroupFile, palenFile, filterByCount=10):
    currentFileDir = os.path.dirname(os.path.realpath(__file__))
    flncReads2Palen = getDictFromFile(palenFile, sep="\t", valueCol=10)
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    gene2readPalen = {}
    for isoform in collapsedTrans2reads:
        geneId = ".".join(isoform.split(".")[0:2])
        reads = collapsedTrans2reads[isoform]
        reads2palen = dict(zip(reads, [int(flncReads2Palen[i]) for i in reads]))
        if geneId not in gene2readPalen:
            gene2readPalen[geneId] = {isoform: reads2palen}
        else:
            gene2readPalen[geneId][isoform] = reads2palen

    for g in gene2readPalen:
        allReadsInGene = sum([len(gene2readPalen[g][i]) for i in gene2readPalen[g]])
        for a, b in itertools.combinations(gene2readPalen[g].keys(), 2):
            if len(gene2readPalen[g][a].keys()) >= filterByCount and len(gene2readPalen[g][b].keys()) >= filterByCount:
                aPalen = map(int, gene2readPalen[g][a].values())
                bPalen = map(int, gene2readPalen[g][b].values())
                stat_val1, p_val1 = stats.ttest_ind(aPalen, bPalen)
                stat_val2, p_val2 = stats.kruskal(aPalen, bPalen)
                if float(p_val1) <= 0.05 and float(p_val2) <= 0.05:
                    # print "\t".join(map(str, [allReadsInGene, g, a, len(aPalen), np.mean(aPalen), np.median(aPalen), b, len(bPalen), np.mean(bPalen), np.median(bPalen)]))
                    fileOut = "_".join([g, a, b]) + ".txt"
                    out = open("_".join([g, a, b]) + ".txt", "w")
                    for i in gene2readPalen[g][a]:
                        print >>out, "\t".join(map(str, [g, a, i, gene2readPalen[g][a][i]]))
                    for j in gene2readPalen[g][b]:
                        print >>out, "\t".join(map(str, [g, b, j, gene2readPalen[g][b][j]]))
                    out.close()
                    # plotGeneStructurePath = os.path.join(currentFileDir, "violin_paired.R")
                    # cmd = "Rscript {} {}".format(plotGeneStructurePath, fileOut)
                    cmd = "Rscript violin_paired.R {}".format(fileOut)
                    subprocess.call(cmd, shell=True)

def relationshipBetweenPAlenAndAS(asFile, collapsedGroupFile, flncFile, isoformBed, asType="SE", filterByCount=5, sigFile=None, sigPicOutDir=None):
    currentFileDir = os.path.dirname(os.path.realpath(__file__))
    flncReads2Palen = getDictFromFile(flncFile, sep=",", valueCol=5)
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    isoformDict = getDictFromFile(isoformBed, keyCol=4, sep="\t")

    asType = asType.upper()
    out = open(sigFile, "w")
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            elif asType == "PA":
                inclusionIsos = []
                exclusionIsos = []
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            for item in itertools.product(inclusionIsos, exclusionIsos):
                inclusionReads = collapsedTrans2reads[item[0]]
                exclusionReads = collapsedTrans2reads[item[1]]
                if filterByCount:
                    if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
                        continue
                inclusionReads2palen = [int(flncReads2Palen[i]) for i in inclusionReads]
                exclusionReads2palen = [int(flncReads2Palen[i]) for i in exclusionReads]
                stat_val, p_val = stats.ttest_ind(inclusionReads2palen, exclusionReads2palen)
                if float(p_val) <= 0.05:
                    print >> out, "\t".join(isoformDict[item[0]][0:12]) + "\t" + "\t".join(
                        ["inclusion", ",".join(map(str, inclusionReads2palen)), str(len(inclusionReads)), "_".join(item)])
                    print >> out, "\t".join(isoformDict[item[1]][0:12]) + "\t" + "\t".join(
                        ["exclusion", ",".join(map(str, exclusionReads2palen)), str(len(exclusionReads)), "_".join(item)])
    out.close()

    if os.stat(sigFile).st_size != 0:
        plotGeneStructurePath = os.path.join(currentFileDir, "plotGeneStructure.R")
        cmd = "Rscript {} -sig={} -od={}".format(plotGeneStructurePath, sigFile, sigPicOutDir)
        subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-i", type=str,
                        help="All isoforms from iFLAS pipeline")
    parser.add_argument("-f", type=str,
                        help="The flnc reads in bed format from iFLAS pipeline")
    parser.add_argument("-a", type=str,
                        help="AS file from iFLAS pipeline")
    parser.add_argument("-c", type=str,
                        help="The collapsed isoform file from cDNA_Cupcake")
    parser.add_argument("-o", type=str, default="palenAndAS.sig.bed12+",
                        help="Output file")
    parser.add_argument("-oPic", type=str, default="sig_pics",
                        help="Output directory of pictures")
    parser.add_argument("-asType", type=str, default="SE",
                        help="The alternative splicing type")
    parser.add_argument("-filterByCount", type=int, default=0,
                        help="Filter the collapsed isoforms by the count of supporting flnc reads")
    args = parser.parse_args()
    relationshipBetweenPAlenAndAS(args.a, args.c, args.f, args.i, asType=args.asType, sigFile=args.o, sigPicOutDir=args.oPic)

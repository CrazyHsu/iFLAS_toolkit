#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: paAnalysisFuncs.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-02-27 18:15:30
Last modified: 2020-02-27 18:15:30
'''

import shutil, subprocess, os, pybedtools
import pandas as pd
from multiprocessing import Pool
from collections import Counter
from commonObjs import *

def paCluster(readsBedList, distance=20, windowSize=3, manner="mode", outHandle=None):
    chroms = [i.chrom for i in readsBedList]
    strands = [i.strand for i in readsBedList]
    if len(set(strands)) != 1: return
    chrom = list(set(chroms))[0]
    strand = list(set(strands))[0]
    if strand == "+":
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromEnd)
        paRangeEnd = sortedPAs[0].chromEnd
        paRangeStart = paRangeEnd - 1
    else:
        sortedPAs = sorted(readsBedList, key=lambda x: x.chromStart)
        paRangeStart = sortedPAs[0].chromStart
        paRangeEnd = paRangeStart + 1
    readNames = [sortedPAs[0].name]
    relPaSites = [1]
    for i in range(1, len(sortedPAs)):
        tmpRead = sortedPAs[i]
        currentPa = tmpRead.chromEnd if strand == "+" else tmpRead.chromStart + 1
        if currentPa - paRangeEnd <= distance:
            paRangeEnd = currentPa
            readNames.append(tmpRead.name)
            relPaSites.append(currentPa - paRangeStart)
        else:
            if strand == "-":
                relPaSites = [paRangeEnd - paRangeStart - i + 1 for i in relPaSites][::-1]
            paCounter = Counter(relPaSites)
            currentCount = reduce(lambda x, y: x + y, [paCounter[z+1] for z in range(windowSize) if z + 1 in paCounter])
            maxCount = currentCount
            for j in range(2, max(paCounter.keys()) - windowSize + 1):
                currentCount = 0
                for z in range(j, j + windowSize):
                    if z + 1 in paCounter:
                        currentCount += paCounter[z + 1]
                if currentCount > maxCount: maxCount = currentCount
            if manner == "mode":
                relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: (x[1], x[0]),
                                       reverse=True)
            elif manner == "downstream":
                relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: x[0], reverse=True)
            else:
                relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: x[0])

            mainPaSite = relSite2Count[0][0]
            if strand == "+":
                paSite = paRangeStart + mainPaSite
            else:
                readNames = readNames[::-1]
                paSite = paRangeEnd - mainPaSite + 1
            readCount = len(readNames)
            freq = []
            for tmp in range(max(relPaSites)):
                if tmp + 1 in paCounter:
                    freq.append(round(paCounter[tmp + 1]/float(len(relPaSites)), 2))
                else:
                    freq.append(round(0, 2))
            print >> outHandle, "\t".join(map(str, [chrom, paRangeStart, paRangeEnd, ",".join(readNames), readCount, strand,
                                               paSite - 1, paSite, ",".join(map(str, relPaSites)), len(relSite2Count),
                                               round(maxCount / float(readCount), 2), "\t".join(map(str, freq))]))
            paRangeStart = currentPa - 1
            paRangeEnd = currentPa
            readNames = [tmpRead.name]
            relPaSites = [1]

    if strand == "-":
        relPaSites = [paRangeEnd - paRangeStart - i + 1 for i in relPaSites][::-1]
    paCounter = Counter(relPaSites)
    currentCount = reduce(lambda x, y: x + y, [paCounter[z + 1] for z in range(windowSize) if z + 1 in paCounter])
    maxCount = currentCount
    for j in range(2, max(paCounter.keys()) - windowSize + 1):
        currentCount = 0
        for z in range(j, j + windowSize):
            if z + 1 in paCounter:
                currentCount += paCounter[z + 1]
        if currentCount > maxCount: maxCount = currentCount
    if manner == "mode":
        relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: (x[1], x[0]),
                               reverse=True)
    elif manner == "downstream":
        relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: x[0], reverse=True)
    else:
        relSite2Count = sorted([(k, v) for k, v in paCounter.iteritems()], key=lambda x: x[0])

    mainPaSite = relSite2Count[0][0]
    if strand == "+":
        paSite = paRangeStart + mainPaSite
    else:
        readNames = readNames[::-1]
        paSite = paRangeEnd - mainPaSite + 1
    readCount = len(readNames)
    freq = []
    for tmp in range(max(relPaSites)):
        if tmp + 1 in paCounter:
            freq.append(round(paCounter[tmp + 1] / float(len(relPaSites)), 2))
        else:
            freq.append(round(0, 2))
    print >> outHandle, "\t".join(map(str, [chrom, paRangeStart, paRangeEnd, ",".join(readNames), readCount, strand,
                                            paSite - 1, paSite, ",".join(map(str, relPaSites)), len(relSite2Count),
                                            round(maxCount / float(readCount), 2), "\t".join(map(str, freq))]))

def getPaCluster(readsBed=None, tofuGroup=None, filterByCount=0, threads=None, paClusterOut=None):
    gene2reads = {}
    readsDict = BedFile(readsBed, type="bed12+").reads
    with open(tofuGroup) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            gene = ".".join(infoList[0].split(".")[:2])
            if gene not in gene2reads:
                gene2reads[gene] = {"allReads": []}
            gene2reads[gene]["allReads"].extend(infoList[1].split(","))
    outHandle = open(paClusterOut, "w")
    for z in gene2reads:
        if filterByCount:
            if len(gene2reads[z]["allReads"]) < filterByCount:
                continue
        readsList = gene2reads[z]["allReads"]
        readsBedList = [readsDict[r] for r in readsList if r in readsDict]
        if len(readsBedList) == 0: continue
        chroms = [r.chrom for r in readsBedList]
        mostChrom = Counter(chroms).most_common(1)[0][0]
        readsBedList = [r for r in readsBedList if r.chrom == mostChrom]
        paCluster(readsBedList, outHandle=outHandle)
    outHandle.close()

def pa(polish_flnc_cluster=None, bed12=None, tofu_group=None, filterByCount=None, threads=None, out=None):
    polish_flnc_cluster_dict = {}
    # trans2reads = {}
    gene2reads = {}
    bed12readsDict = BedFile(bed12, type="bed12+").reads
    with open(polish_flnc_cluster) as f:
        for i in f.readlines()[1:]:
            infoList = i.strip("\n").split(",")
            if infoList[0] not in polish_flnc_cluster_dict:
                polish_flnc_cluster_dict[infoList[0]] = []
            polish_flnc_cluster_dict[infoList[0]].append(infoList[1])
    with open(tofu_group) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            gene = ".".join(infoList[0].split(".")[:2])
            if gene not in gene2reads:
                gene2reads[gene] = {"allReads": []}
            readsList = []
            for j in infoList[1].split(","):
                if j in polish_flnc_cluster_dict:
                    gene2reads[gene].update({j: polish_flnc_cluster_dict[j]})
                    readsList.extend(polish_flnc_cluster_dict[j])
            gene2reads[gene]["allReads"].extend(readsList)
    outHandle = open(out, "w")
    for z in gene2reads:
        if filterByCount:
            if len(gene2reads[z]["allReads"]) < filterByCount:
                continue
        readsList = gene2reads[z]["allReads"]
        readsBedList = [bed12readsDict[r] for r in readsList if r in bed12readsDict]
        if len(readsBedList) == 0: continue
        chroms = [r.chrom for r in readsBedList]
        mostChrom = Counter(chroms).most_common(1)[0][0]
        readsBedList = [r for r in readsBedList if r.chrom == mostChrom]
        paCluster(readsBedList, outHandle=outHandle)
    outHandle.close()

def motifAroundPA(bed6plus=None, up1=100, down1=100, up2=100, down2=100, refFasta=None, chrLenFile=None):
    singleNucleotideMotif = ["A", "T", "C", "G"]
    sixNucleotideMotif = ["AATAAA", "AAATAA", "ATAAAA", "ATTAAA", "ATAAAT", "TAATAA",
                          "ATAAAG", "AAAATA", "CAATAA", "ATAAAC", "AAAAAA", "AAAAAG"]
    chrLenDict = {}
    with open(chrLenFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            chrLenDict[infoList[0]] = int(infoList[1])
    with open(bed6plus) as f:
        singleNucleotideUpBedList = []
        singleNucleotideDownBedList = []
        sixNucleotideUpBedList = []
        sixNucleotideDownBedList = []
        for i in f.readlines():
            bedObj = Bed6Plus(i)
            chrom, start, end, strand = bedObj.chrom, bedObj.chromStart, bedObj.chromEnd, bedObj.strand
            name = bedObj.name
            if strand == "+":
                if start > up1 and end + down1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up1, start, name, ".", strand])))
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down1, name, ".", strand])))
                if start > up2 and end + down2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, start - up2, start + 5, name, ".", strand])))
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, end, end + down2 + 5, name, ".", strand])))
            else:
                if start > down1 and end + up1 < chrLenDict[chrom]:
                    singleNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up1, name, ".", strand])))
                    singleNucleotideDownBedList.append(" ".join(map(str, [chrom, start - down1 - 1, start-1, name, ".", strand])))
                if start > down2 and end + up2 + 5 < chrLenDict[chrom]:
                    sixNucleotideUpBedList.append(" ".join(map(str, [chrom, end, end + up2 + 5, name, ".", strand])))
                    sixNucleotideDownBedList.append(" ".join(map(str, [chrom, start - down2 - 5, start, name, ".", strand])))
        singleNucleotideUpBedObj = pybedtools.BedTool("\n".join(singleNucleotideUpBedList), from_string=True)
        singleNucleotideDownBedObj = pybedtools.BedTool("\n".join(singleNucleotideDownBedList), from_string=True)
        sixNucleotideUpBedObj = pybedtools.BedTool("\n".join(sixNucleotideUpBedList), from_string=True)
        sixNucleotideDownBedObj = pybedtools.BedTool("\n".join(sixNucleotideDownBedList), from_string=True)

        singleNucleotideUpBedFastaRes = singleNucleotideUpBedObj.sequence(refFasta, name=True, tab=True, s=True)
        singleNucleotideDownBedFastaRes = singleNucleotideDownBedObj.sequence(refFasta, name=True, tab=True, s=True)
        sixNucleotideUpBedFastaRes = sixNucleotideUpBedObj.sequence(refFasta, name=True, tab=True, s=True)
        sixNucleotideDownBedFastaRes = sixNucleotideDownBedObj.sequence(refFasta, name=True, tab=True, s=True)

        singleNucleotideUpBedFastaDict = {}
        singleNucleotideDownBedFastaDict = {}
        sixNucleotideUpBedFastaDict = {}
        sixNucleotideDownBedFastaDict = {}
        for res in str(open(singleNucleotideUpBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            singleNucleotideUpBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i+1] for i in range(0, len(infoList[1]), 1)]
        for res in str(open(singleNucleotideDownBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            singleNucleotideDownBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 1] for i in range(0, len(infoList[1]), 1)]
        for res in str(open(sixNucleotideUpBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            sixNucleotideUpBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 6] for i in range(0, len(infoList[1]) - 5, 1)]
        for res in str(open(sixNucleotideDownBedFastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = res.split("\t")
            sixNucleotideDownBedFastaDict[infoList[0]] = [infoList[1].upper()[i:i + 6] for i in range(0, len(infoList[1]) - 5, 1)]
        singleNucleotideUpDf = pd.DataFrame.from_dict(singleNucleotideUpBedFastaDict, orient="index")
        singleNucleotideDownDf = pd.DataFrame.from_dict(singleNucleotideDownBedFastaDict, orient="index")
        sixNucleotideUpDf = pd.DataFrame.from_dict(sixNucleotideUpBedFastaDict, orient="index")
        sixNucleotideDownDf = pd.DataFrame.from_dict(sixNucleotideDownBedFastaDict, orient="index")

        for i in singleNucleotideMotif:
            singleNucleotideOut = open(i+".nucleotide", "w")
            singleNucleotideUpDf1 = singleNucleotideUpDf == i
            upPercent = singleNucleotideUpDf1.sum() / float(len(singleNucleotideUpDf1))
            for j in range(len(upPercent)):
                print >>singleNucleotideOut, "\t".join(map(str, [j - up1, upPercent[j]]))

            singleNucleotideDownDf1 = singleNucleotideDownDf == i
            downPercent = singleNucleotideDownDf1.sum() / float(len(singleNucleotideDownDf1))
            for j in range(len(downPercent)):
                print >>singleNucleotideOut, "\t".join(map(str, [j + 1, downPercent[j]]))
            singleNucleotideOut.close()

        for i in sixNucleotideMotif:
            sixNucleotideOut = open(i + ".PAS", "w")
            sixNucleotideUpDf1 = sixNucleotideUpDf == i
            upPercent = sixNucleotideUpDf1.sum() / float(len(sixNucleotideUpDf1))
            for j in range(len(upPercent)):
                print >>sixNucleotideOut, "\t".join(map(str, [j - up2, upPercent[j]]))

            sixNucleotideDownDf1 = sixNucleotideDownDf == i
            downPercent = sixNucleotideDownDf1.sum() / float(len(sixNucleotideDownDf1))
            for j in range(len(downPercent)):
                print >>sixNucleotideOut, "\t".join(map(str, [j + 1, downPercent[j]]))
            sixNucleotideOut.close()


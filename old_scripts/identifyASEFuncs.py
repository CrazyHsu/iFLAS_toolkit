#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identifyASEFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 19:35:06
Last modified: 2019-12-08 19:35:07
'''

from commonObjs import *
from collections import Counter
import pybedtools
# from findAS import *
import pandas as pd

def intronBoundaryMatch(myIntronChain, compIntronChain, skippedIntron=False, offset=0):
    returnFlag = False
    compIntronBoundList = listStr2Int(re.split('[-|;]', compIntronChain))
    if skippedIntron:
        left, right = listStr2Int(re.split('[-|;]', myIntronChain)[1:-1])
        for i in range(1, len(compIntronBoundList)-1, 2):
            itemLeft = compIntronBoundList[i]
            itemRight = compIntronBoundList[i+1]
            if abs(itemLeft-left) <= offset and abs(itemRight - right) <= offset:
                returnFlag = True
    else:
        keepBoundList = listStr2Int(re.split('[-|;]', myIntronChain)[1:-1])
        for i in range(1, len(compIntronBoundList)-1, len(keepBoundList)):
            keepFlag = True
            for j in range(len(keepBoundList)):
                if abs(keepBoundList[i] - compIntronBoundList[i+j]) > offset:
                    keepFlag = False
                    break
            returnFlag = keepFlag
    return returnFlag

def findAS(gene2ReadsDict, outASType=None, anno=True, out=None, isoform2reads=None):
    findASType = {"IR": findIR, "SE": findSE, "A3SS": findA3SS, "A5SS": findA5SS}
    # if outASType and outASType.upper() not in findASType:
    #     raise Exception("You should input the right AS type!")
    if anno:
        if outASType:
            if outASType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[outASType](gene2ReadsDict, outASType, out, isoform2reads=isoform2reads)
        else:
            for asType in findASType:
                findASType[asType](gene2ReadsDict, isoform2reads=isoform2reads)
    else:
        newGene2ReadsDict = {}
        for myChr, myValue in gene2ReadsDict.iteritems():
            for myStrand in myValue:
                sortedStrandV = sorted(myValue[myStrand], key=lambda x: x.chromStart)
                count = 1
                fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                gene2reads = Gene2Reads(fakeGeneName)
                gene2reads.update(sortedStrandV[0])
                clusterEnd = sortedStrandV[0].chromEnd
                for i in range(1, len(sortedStrandV)):
                    myStart, myEnd = sortedStrandV[i].chromStart, sortedStrandV[i].chromEnd
                    if myStart < clusterEnd:
                        gene2reads.update(sortedStrandV[i])
                        if myEnd > clusterEnd:
                            clusterEnd = myEnd
                    else:
                        newGene2ReadsDict[fakeGeneName] = gene2reads
                        count += 1
                        fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                        gene2reads = Gene2Reads(fakeGeneName)
                        gene2reads.update(sortedStrandV[i])
                        clusterEnd = sortedStrandV[i].chromEnd
                newGene2ReadsDict[fakeGeneName] = gene2reads
        if outASType:
            if outASType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[outASType](newGene2ReadsDict, outASType, out, isoform2reads=isoform2reads, anno=False)
        else:
            for asType in findASType:
                findASType[asType](newGene2ReadsDict, isoform2reads=isoform2reads)

def findIR(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        irDict = {}
        for rName1 in rNames:
            exonStarts, exonEnds = readsDict[rName1].exonStarts, readsDict[rName1].exonEnds
            for exon1 in range(len(exonStarts) - 1):
                IR = "%d-%d" % (exonEnds[exon1], exonStarts[exon1 + 1])
                if IR in irDict: continue
                irDict[IR] = {}
                if "spliced" not in irDict[IR]:
                    irDict[IR].__setitem__("spliced", [readsDict[rName1]])
                else:
                    irDict[IR]["spliced"].append(readsDict[rName1])
                for rName2 in rNames:
                    if rName1 == rName2: continue
                    exons2 = readsDict[rName2].exonChain.split(";")
                    for exon2 in range(len(exons2)):
                        exon2Start, exon2End = [int(x) for x in exons2[exon2].split("-")]
                        # if exon2Start < exonEnds[exon1] and exonStarts[exon1 + 1] < exon2End and abs(exonStarts - exonEnds[exon1]) > offset and abs(exonStarts[exon1 + 1] - exon2End) > offset:
                        if exon2Start < exonEnds[exon1] and exonStarts[exon1 + 1] < exon2End \
                                and abs(exon2Start - exonEnds[exon1]) > offset \
                                and abs(exonStarts[exon1 + 1] - exon2End) > offset:
                            if "retention" not in irDict[IR]:
                                irDict[IR].__setitem__("retention", [readsDict[rName2]])
                                break
                            else:
                                irDict[IR]["retention"].append(readsDict[rName2])
                        if exon2End == exonEnds[exon1]:
                            if exon2 + 1 >= len(exons2): break
                            exon2NStart, exon2NEnd = [int(x) for x in exons2[exon2 + 1].split("-")]
                            if exon2NStart == exonStarts[exon1 + 1]:
                                irDict[IR]["spliced"].append(readsDict[rName2])
                        if exon2Start > exonEnds[exon1 + 1]: break
        newIrDict = {}
        for ir in irDict:
            if "retention" in irDict[ir]:
                newIrDict[ir] = irDict[ir]
                if outAS:
                    irStart, irEnd = [int(x) for x in ir.split("-")]
                    if anno:
                        transDict = gene2ReadsDict[gName].trans
                        overlapWithGene = 0
                        for tran in transDict:
                            if irStart < transDict[tran].chromEnd and irEnd > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            retentionReads = [x.name for x in irDict[ir]["retention"]]
                            splicedReads = [x.name for x in irDict[ir]["spliced"]]
                            psi = float(len(retentionReads))/(len(retentionReads)+len(splicedReads))*1000
                            if isoform2reads:
                                rawRetentionReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in retentionReads]))
                                rawSplicedReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in splicedReads]))
                                psi = float(len(rawRetentionReads))/(len(rawRetentionReads)+len(rawSplicedReads))*1000
                            print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName+":"+ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
                    else:
                        retentionReads = [x.name for x in irDict[ir]["retention"]]
                        splicedReads = [x.name for x in irDict[ir]["spliced"]]
                        psi = float(len(retentionReads)) / (len(retentionReads) + len(splicedReads)) * 1000
                        if isoform2reads:
                            rawRetentionReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in retentionReads]))
                            rawSplicedReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in splicedReads]))
                            psi = float(len(rawRetentionReads))/(len(rawRetentionReads)+len(rawSplicedReads))*1000
                        print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName + ":" + ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
        gene2ReadsDict[gName].asDict.update({"IR": newIrDict})


def findSE(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        exonDict = {}
        for rName1 in rNames:
            exons = readsDict[rName1].exonChain.split(";")
            for l in range(len(exons) - 2):
                for r in range(l + 2, len(exons)):
                    lExonEnd = int(exons[l].split("-")[1])
                    rExonStart = int(exons[r].split("-")[0])
                    SEs = ";".join(exons[l + 1:r])
                    seChain = "%d@%s@%d" % (lExonEnd, SEs, rExonStart)
                    if seChain in exonDict: continue
                    exonDict[seChain] = {}
                    if "keep" not in exonDict[seChain]:
                        exonDict[seChain].__setitem__("keep", [readsDict[rName1]])
                    else:
                        exonDict[seChain]["keep"].append(readsDict[rName1])
                    for rName2 in rNames:
                        if rName2 == rName1: continue
                        skipPat = "-%d;%d-" % (lExonEnd, rExonStart)
                        keepPat = "-%d;%s;%d-" % (lExonEnd, SEs, rExonStart)
                        # if intronBoundaryMatch(skipPat, readsDict[rName2].exonChain, skippedIntron=True, offset=offset):
                        #     if "skip" not in exonDict[seChain]:
                        #         exonDict[seChain].__setitem__("skip", [readsDict[rName2]])
                        #     else:
                        #         exonDict[seChain]["skip"].append(readsDict[rName2])
                        # elif intronBoundaryMatch(keepPat, readsDict[rName2].exonChain, skippedIntron=False, offset=offset):
                        #     exonDict[seChain]["keep"].append(readsDict[rName2])
                        if re.search(skipPat, readsDict[rName2].exonChain):
                            if "skip" not in exonDict[seChain]:
                                exonDict[seChain].__setitem__("skip", [readsDict[rName2]])
                            else:
                                exonDict[seChain]["skip"].append(readsDict[rName2])
                        elif re.search(keepPat, readsDict[rName2].exonChain):
                            exonDict[seChain]["keep"].append(readsDict[rName2])
        newExonDict = {}
        for seChain1 in exonDict.keys():
            if "skip" in exonDict[seChain1]:
                newExonDict[seChain1] = exonDict[seChain1]
                if outAS:
                    keepReads = [x.name for x in exonDict[seChain1]["keep"]]
                    skipReads = [x.name for x in exonDict[seChain1]["skip"]]
                    keepReadsCount = len(keepReads)
                    skipReadsCount = len(skipReads)
                    lSpliceSite, SEs1, rSpliceSite = seChain1.split("@")
                    SElist = SEs1.split(";")
                    starts, ends = [], []
                    for exon in SElist:
                        exonStart, exonEnd = listStr2Int(exon.split("-"))
                        starts.append(exonStart)
                        ends.append(exonEnd)
                    if anno:
                        overlapWithGene = 0
                        transDict = gene2ReadsDict[gName].trans
                        for tran in transDict:
                            if starts[0] < transDict[tran].chromEnd and ends[-1] > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                            blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                            psi = float(keepReadsCount)/(keepReadsCount + skipReadsCount) * 1000
                            if isoform2reads:
                                rawKeepReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in keepReads]))
                                rawSkipReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in skipReads]))
                                psi = float(len(rawKeepReads)) / (len(rawKeepReads) + len(rawSkipReads)) * 1000
                            print >> out, "\t".join(map(str, [chrom, str(starts[0]), str(ends[-1]),
                                                              "%s:%s" % (gName, seChain1), psi, strand, starts[0],
                                                              ends[-1], "0,0,0", len(starts), blockSizes,
                                                              blockRelStarts, gName, lSpliceSite, rSpliceSite,
                                                              ",".join(keepReads), len(keepReads), ",".join(skipReads),
                                                              len(skipReads)]))
                    else:
                        blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                        blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                        psi = float(keepReadsCount) / (keepReadsCount + skipReadsCount) * 1000
                        if isoform2reads:
                            rawKeepReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in keepReads]))
                            rawSkipReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in skipReads]))
                            psi = float(len(rawKeepReads)) / (len(rawKeepReads) + len(rawSkipReads)) * 1000
                        print >> out, "\t".join(map(str,
                                                    [chrom, str(starts[0]), str(ends[-1]), "%s:%s" % (gName, seChain1),
                                                     psi, strand, starts[0], ends[-1], "0,0,0", len(starts), blockSizes,
                                                     blockRelStarts, gName, lSpliceSite, rSpliceSite,
                                                     ",".join(keepReads), len(keepReads), ",".join(skipReads),
                                                     len(skipReads)]))
        gene2ReadsDict[gName].asDict.update({"SE": newExonDict})


def findA3SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 1 if strand == "+" else 0
        for r1 in range(len(rNames) - 1):
            exonsChain1 = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain1.split(";")
            indexEnd1 = len(exons1) if strand == "+" else len(exons1) - 1
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].exonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) if strand == "+" else len(exons2) - 1
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "-":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon2End - int(exons1[i + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon1End - int(exons2[j + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon1Start - int(exons2[j - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon2Start - int(exons1[i - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            if outAS:
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
                incJunc = altDict[alt]["incJunc"]
                excJunc = altDict[alt]["excJunc"]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                if isoform2reads:
                    rawIncReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in incReads]))
                    rawExcReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in excReads]))
                    psi = float(len(rawIncReads)) / (len(rawIncReads) + len(rawExcReads)) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount), \
                                          "inc:{}-{}".format(incJunc[0], incJunc[1]), \
                                          "exc:{}-{}".format(excJunc[0], excJunc[1])]))

        gene2ReadsDict[gName].asDict.update({"A3SS": altDict})


def findA5SS(gene2ReadsDict, outAS=None, out=None, anno=True, offset=0, isoform2reads=None):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 0 if strand == "+" else 1
        for r1 in range(len(rNames) - 1):
            exonsChain1 = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain1.split(";")
            indexEnd1 = len(exons1) - 1 if strand == "+" else len(exons1)
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].exonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) - 1 if strand == "+" else len(exons2)
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "+":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon2End - int(exons1[i + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]) \
                                    and abs(exon1End - exon2End) > offset \
                                    and abs(exon1End - int(exons2[j + 1].split("-")[0])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (exon2End, int(exons2[j + 1].split("-")[0]))})
                                    altDict[exonBoundChain].update({"incJunc": (exon1End, int(exons1[i + 1].split("-")[0]))})
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon1Start - int(exons2[j - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]) \
                                    and abs(exon1Start - exon2Start) > offset \
                                    and abs(exon2Start - int(exons1[i - 1].split("-")[1])) > offset:
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    altDict[exonBoundChain].update({"excJunc": (int(exons1[i - 1].split("-")[1]), exon1Start)})
                                    altDict[exonBoundChain].update({"incJunc": (int(exons2[j - 1].split("-")[1]), exon2Start)})
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            incJunc = altDict[alt]["incJunc"]
            excJunc = altDict[alt]["excJunc"]
            if outAS:
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                if isoform2reads:
                    rawIncReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in incReads]))
                    rawExcReads = list(itertools.chain.from_iterable([isoform2reads[x] for x in excReads]))
                    psi = float(len(rawIncReads)) / (len(rawIncReads) + len(rawExcReads)) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount), \
                                          "inc:{}-{}".format(incJunc[0], incJunc[1]), \
                                          "exc:{}-{}".format(excJunc[0], excJunc[1])]))
        gene2ReadsDict[gName].asDict.update({"A5SS": altDict})

def getAS_bak(gpeDict, geneObj, annoDict, novelDict, asType, offset=0):
    # annoDict, novelDict = {}, {}
    chrom, strand, reads = geneObj.chrom, geneObj.strand, geneObj.reads.values()
    geneObjName, geneObjMinpos, geneObjMaxpos = geneObj.geneName, geneObj.minpos, geneObj.maxpos
    if chrom in gpeDict:
        isNovel = 1
        if strand in gpeDict[chrom]:
            for trans in gpeDict[chrom][strand]:
                if geneObjMinpos < trans[1] and geneObjMaxpos > trans[0]:
                    transStarts, transEnds = trans[2].exonStarts, trans[2].exonEnds
                    for read in reads:
                        readStarts, readEnds, blockCount = read.blockStarts, read.blockEnds, read.blockNum
                        if blockCount < 2: continue
                        consensusIntronN = getConsensusIntronN(readStarts, readEnds, transStarts, transEnds, asType, offset)
                        if consensusIntronN >= 1:
                            transName, geneName = trans[2].transName, trans[2].geneName
                            if geneName not in annoDict:
                                gene2Reads = Gene2Reads(geneObjName)
                                gene2Reads.update(read)
                                gene2Reads.trans = {transName: trans}
                                annoDict[geneName] = gene2Reads
                            else:
                                gene2Reads = annoDict[geneName]
                                if read.readName not in [i.readName for i in gene2Reads.reads.values()]:
                                    gene2Reads.update(read)
                                gene2Reads.trans.update({transName: trans})
                                annoDict[geneName] = gene2Reads
                            isNovel = 0
                if trans[0] > geneObjMaxpos: break
        if isNovel == 1:
            if chrom not in novelDict:
                novelDict[chrom] = {strand: [geneObj]}
            elif strand not in novelDict[chrom]:
                novelDict[chrom][strand] = [geneObj]
            else:
                novelDict[chrom][strand].append(geneObj)
    else:
        '''
        The situation that the chrom is not in annotation file(gpe file)
        '''
        pass

def getAS(annoBedRes, novelBedRes, offset=0, reference=0):
    annoDict, novelDict = {}, {}
    isNovelDict = {}
    for line in annoBedRes:
        record = str(line).strip("\n").split("\t")
        read = Bed12("\t".join(record[:12]))
        trans = ReadLineStruc("\t".join(record[12:]))
        if reference and record[12] != record[-1]: continue
        if read.name not in isNovelDict:
            isNovelDict[read.name] = {"read": read, "isNovel": 1}
        transStarts, transEnds = trans.exonStarts, trans.exonEnds
        readStarts, readEnds, blockCount = read.exonStarts, read.exonEnds, read.blockCount
        # if blockCount < 2: continue
        consensusIntronN = getConsensusIntronNfromAS(readStarts, readEnds, transStarts, transEnds, offset)
        if consensusIntronN >= 1:
            transName, geneName = trans.name, trans.geneName
            if geneName not in annoDict:
                gene2reads = Gene2Reads(geneName)
                gene2reads.update(read)
                gene2reads.trans = {transName: trans}
                annoDict[geneName] = gene2reads
            else:
                gene2reads = annoDict[geneName]
                # if read.name not in gene2reads.readNames:
                gene2reads.update(read)
                gene2reads.trans.update({transName: trans})
                annoDict[geneName] = gene2reads
            isNovelDict[read.name]["isNovel"] = 0
        else:
            if blockCount == 1:
                transName, geneName = trans.name, trans.geneName
                if isOverlap((read.chromStart, read.chromEnd), (trans.chromStart, trans.chromEnd)):
                    overlap = getOverlapOfTuple([(read.chromStart, read.chromEnd)], [(trans.chromStart, trans.chromEnd)])
                    sumOverlap = sum([x[1] - x[0]for x in overlap])
                    if sumOverlap/float(read.chromEnd-read.chromStart) >= 0.5:
                        if geneName not in annoDict:
                            gene2reads = Gene2Reads(geneName)
                            gene2reads.update(read)
                            gene2reads.trans = {transName: trans}
                            annoDict[geneName] = gene2reads
                        else:
                            gene2reads = annoDict[geneName]
                            gene2reads.update(read)
                            gene2reads.trans.update({transName: trans})
                            annoDict[geneName] = gene2reads
                        isNovelDict[read.name]["isNovel"] = 0

    for line in novelBedRes:
        read = Bed12(str(line))
        if read.chrom not in novelDict:
            novelDict[read.chrom] = {read.strand: [read]}
        elif read.strand not in novelDict[read.chrom]:
            novelDict[read.chrom][read.strand] = [read]
        else:
            novelDict[read.chrom][read.strand].append(read)

    for i in isNovelDict:
        if isNovelDict[i]["isNovel"]:
            read = isNovelDict[i]["read"]
            if read.chrom not in novelDict:
                novelDict[read.chrom] = {read.strand: [read]}
            elif read.strand not in novelDict[read.chrom]:
                novelDict[read.chrom][read.strand] = [read]
            else:
                novelDict[read.chrom][read.strand].append(read)
    return annoDict, novelDict

def findASmain_bak(asType="IR", annoDict=None, novelDict=None, outFile=None, isoform2reads=None):
    out = open(outFile, "w")
    findAS(annoDict, outASType=asType, anno=True, out=out, isoform2reads=isoform2reads)
    findAS(novelDict, outASType=asType, anno=False, out=out, isoform2reads=isoform2reads)
    out.close()

def findASmain(asType="IR", gpeFile=None, sourceFile=None, outFile=None):
    transBedList = GenePredObj(gpeFile, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    with open(sourceFile) as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            readsBedList.append("\t".join(readStruc.record[:12]))
    readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)

    annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
    novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
    annoDict, novelDict = getAS(annoBedRes, novelBedRes, asType, offset=0)
    out = open(outFile, "w")
    findAS(annoDict, outASType=asType, anno=True, out=out)
    findAS(novelDict, outASType=asType, anno=False, out=out)
    out.close()

# def findASmain_bak(asType="IR", gpeFile=None, sourceFile=None, outFile=None):
#     gpeObj = GenePredObj(gpeFile, False)
#     gpeDict = gpeObj.genePredDict
#     sampleObj = Sample()
#     sample = "temp"
#     with open(sourceFile) as f:
#         for line in f:
#             readStruc = ReadLineStruc(line)
#             sampleObj.update(sample, sourceFile, readStruc)
#     geneDict, novelDict = {}, {}
#     for g in sampleObj.sample2Gene[sample]:
#         getAS(gpeDict, sampleObj.sample2Gene[sample][g], geneDict, novelDict, offset=0)
#     out = open(outFile, "w")
#     findAS(geneDict, outASType=asType, anno=True, out=out)
#     findAS(novelDict, outASType=asType, anno=False, out=out)
#     out.close()

def getAnnoASList(inFile, outFile, PA=False, append=False, uniq=False):
    with open(inFile) as f:
        if not append:
            out = open(outFile, "w")
        else:
            out = open(outFile, "w+")
        asLst = []
        if not PA:
            for line in f:
                if "Novel" in line: continue
                lineInfo = line.strip("\n").split("\t")
                asEvent = "{}:{}:{}".format(lineInfo[0], lineInfo[3], lineInfo[5])
                if uniq:
                    if asEvent not in asLst:
                        asLst.append(asEvent)
                    else:
                        continue
                print >>out, asEvent
        else:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                paPos = lineInfo[3].split(",")
                if len(paPos) > 1:
                    asEvent = ":".join(lineInfo[0:4])
                    if uniq:
                        if asEvent not in asLst:
                            asLst.append(asEvent)
                        else:
                            continue
                    print >>out, asEvent
        out.close()

def getASstatistics(asType="IR", asFile=None, annoFile=None, novelFile=None, outFile=None):
    out = open(outFile, "w")
    if asType == "PA":
        print >> out, "#Chr\tStrand\tKnown or Novel\tGene\tPA Sites"
        with open(annoFile) as f1:
            for line in f1:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Known", lineInfo[2], lineInfo[3]])
        with open(novelFile) as f2:
            for line in f2:
                lineInfo = line.strip("\n").split("\t")
                if len(lineInfo[3].split(",")) > 1:
                    print >>out, "\t".join([lineInfo[0], lineInfo[1], "Novel", lineInfo[2], lineInfo[3]])
    else:
        asDict = {}
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                key = ":".join([lineInfo[0], lineInfo[3], lineInfo[5]])
                asDict[key] = "\t".join([key, line.strip("\n")])
        annoList = filter(originFile=annoFile, targetFile=asDict, returnFlag=True)
        novelList = filter(originFile=novelFile, targetFile=asDict, returnFlag=True)
        os.remove("filtered.txt")
        geneCol = 7
        if asType == "IR":
            print >> out, "##AS ID is composed of Gene:Retained intron start-Retained intron end"
        elif asType == "SE":
            geneCol = 13
            print >> out, "##AS ID is composed of Gene:Left flanking constitutive exon end@Alternative exons locus@Right flanking constitutive exon start"
            print >> out, "##Alternative exons locus is composed of Alternative exon1 start-Alternative exon1 end[;Alternative exon2 start-Alternative exon2 end[;Alternative exon3 start-Alternative exon3 end...]"
        elif asType == "A3SS":
            print >> out, "##AS ID is composed of Gene:Alternative 5' splicing region start-Alternative 5' splicing region end"
        elif asType == "A5SS":
            print >> out, "##AS ID is composed of Gene:Alternative 3' splicing region start-Alternative 3' splicing region end"
        print >> out, "#Chr\tStrand\tKnown or Novel\tAS ID\tGene"

        for i in annoList:
            tmpList = asDict[i].strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Known", tmpList[4], tmpList[geneCol]])
        for i in novelList:
            tmpList = asDict[i].strip().split("\t")
            print >> out, "\t".join([tmpList[1], tmpList[6], "Novel", tmpList[4], tmpList[geneCol]])

    out.close()

def getPSICorrBetweenPBandNGS(file1, file2, outPdf=None, outCor=None):
    f1df = pd.read_csv(file1, sep="\t", header=None)
    f2df = pd.read_csv(file2, sep="\t", header=None)
    f1sub = f1df.iloc[:, [3,4]]
    f2sub = f2df.iloc[:, [3,4]]
    f1sub.rename({3: 'asEvent', 4: 'psi'}, axis="columns", inplace=True)
    f2sub.rename({3: 'asEvent', 4: 'psi'}, axis="columns", inplace=True)
    f1sub.loc[:, "psi"] = f1sub.loc[:, "psi"].apply(lambda x: float(x)/1000)
    f2sub.loc[:, "psi"] = f2sub.loc[:, "psi"].apply(lambda x: 1 - float(x) / 1000)
    mergedDf = pd.merge(f1sub, f2sub, on=["asEvent"])
    mergedDf.iloc[:, [1, 2]].to_csv("tmpPSI.txt", sep="\t", header=None, index=None)
    cmd = "correlation.R -x=PB -y=NGS -p={} tmpPSI.txt >{} 2>/dev/null".format(outPdf, outCor)
    subprocess.call(cmd, shell=True)

def drawSSmotif(asMotif=None, outPrefix=None):
    with open(asMotif) as f:
        lineList = f.readlines()
        mySum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList])
        tmp = open("tmp.txt", "w")
        if len(lineList) < 4:
            for i in lineList:
                lineInfo = i.strip("\n").split("\t")
                print >> tmp, "\t".join([lineInfo[0], str(float(lineInfo[1])/mySum)])
        else:
            for i in lineList[0:3]:
                lineInfo = i.strip("\n").split("\t")
                print >> tmp, "\t".join([lineInfo[0], str(float(lineInfo[1]) / mySum)])
            otherSum = sum([int(i.strip("\n").split("\t")[1]) for i in lineList[3:]])
            print >> tmp, "\t".join(["Other", str(float(otherSum)/mySum)])
        tmp.close()
        cmd = "cat tmp.txt | bar.R -fillV=V1 -fp -lgPos=top -w=12 -p={}.ssMotif.pdf 2>/dev/null".format(outPrefix)
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        os.remove("tmp.txt")
        os.remove(asMotif)

def irConfirmByJunc_1(irBed, juncBed, offset=0, juncBedAssigned=False, ngsOutFile=None, pbOutFile=None):
    junc = BedFile(juncBed, type="bed12+")
    ir = BedFile(irBed, type="bed6")
    ngsOut = open(ngsOutFile, "w")
    pbOut = open(pbOutFile, "w")
    if juncBedAssigned:
        colIsIr = -3
    else:
        colIsIr = -1
    juncDict = {}
    for i in junc.reads:
        juncRecord = junc.reads[i]
        juncTuple = (juncRecord.exonEnds[0], juncRecord.exonStarts[1])
        if juncRecord.chrom not in juncDict:
            juncDict[juncRecord.chrom] = {juncTuple: juncRecord.record}
        else:
            juncDict[juncRecord.chrom].update({juncTuple: juncRecord.record})
    for j in ir.reads:
        irRecord = ir.reads[j]
        irTuple = (irRecord.chromStart, irRecord.chromEnd)
        if not irRecord.chrom in juncDict: continue
        if irTuple in juncDict[irRecord.chrom]:
            if juncDict[irRecord.chrom][irTuple][colIsIr] == "yes":
                print >>ngsOut, "\t".join(juncDict[irRecord.chrom][irTuple])
                print >>pbOut, "\t".join(irRecord.record)
        else:
            juncTuples = juncDict[irRecord.chrom]
            overlappedJunc = []
            for t in juncTuples:
                if irTuple[0] > t[1]:break
                if isOverlap(irTuple, t) and abs(irTuple[0] - t[0]) <= offset and abs(irTuple[1] - t[1]):
                    overlappedJunc.append(t)
            if overlappedJunc:
                for z in overlappedJunc:
                    print >>ngsOut, "\t".join(juncDict[irRecord.chrom][z])
                print >>pbOut, "\t".join(irRecord.record)
    ngsOut.close()
    pbOut.close()

def irConfirmByJunc(irBed, juncBed, ngsOutFile=None, pbOutFile=None):
    junc = BedFile(juncBed, type="bed12")
    ir = BedFile(irBed, type="bed6")
    ngsOut = open(ngsOutFile, "w")
    pbOut = open(pbOutFile, "w")
    juncDict = {}
    for i in junc.reads:
        juncRecord = junc.reads[i]
        juncTuple = (juncRecord.exonEnds[0], juncRecord.exonStarts[1])
        juncDict[juncTuple] = juncRecord.record
    for j in ir.reads:
        irRecord = ir.reads[j]
        irTuple = (irRecord.chromStart, irRecord.chromEnd)
        if irTuple in juncDict:
            print >>ngsOut, "\t".join(juncDict[irTuple])
            print >>pbOut, "\t".join(irRecord.record)
    ngsOut.close()
    pbOut.close()

def AnSSconfirmByJunc():
    pass

def splicesite2seq(refFasta, spliceSite, noFiveEnd=False, noThreeEnd=False, outFile=None):
    donorList, acceptorList = [], []
    exonDist, intronDist = 0, 2
    with open(spliceSite) as f:
        count = 0
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            chrom, leftSite, rightSite, strand = infoList[0], infoList[1], infoList[2], infoList[3]
            count += 1
            juncName = "junction_" + str(count)
            if infoList[3] == "+":
                if not noFiveEnd:
                    leftSite = int(leftSite)
                    donorBed = "\t".join(map(str, [chrom, leftSite - exonDist, leftSite + intronDist, juncName, ".", strand]))
                    donorList.append(donorBed)
                if not noThreeEnd:
                    rightSite = int(rightSite)
                    acceptorBed = "\t".join(map(str, [chrom, rightSite - intronDist, rightSite + exonDist, juncName, ".", strand]))
                    acceptorList.append(acceptorBed)
            else:
                if not noFiveEnd:
                    rightSite = int(rightSite)
                    donorBed = "\t".join(map(str, [chrom, rightSite - intronDist, rightSite + exonDist, juncName, ".", strand]))
                    donorList.append(donorBed)
                if not noThreeEnd:
                    leftSite = int(leftSite)
                    acceptorBed = "\t".join(map(str, [chrom, leftSite - exonDist, leftSite + intronDist, juncName, ".", strand]))
                    acceptorList.append(acceptorBed)
    donorBedObj = pybedtools.BedTool("\n".join(donorList), from_string=True)
    acceptorBedObj = pybedtools.BedTool("\n".join(acceptorList), from_string=True)
    junc2seq = {}
    if not noFiveEnd and not noThreeEnd:
        donorBedGetfastaRes = donorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        acceptorBedGetfastaRes = acceptorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(donorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
        for i in str(open(acceptorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = junc2seq[junc] + "-" + baseSeq
    elif not noFiveEnd and noThreeEnd:
        donorBedGetfastaRes = donorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(donorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
    elif noFiveEnd and not noThreeEnd:
        acceptorBedGetfastaRes = acceptorBedObj.sequence(refFasta, name=True, tab=True, s=True)
        for i in str(open(acceptorBedGetfastaRes.seqfn).read()).split("\n")[:-1]:
            infoList = str(i).strip("\n").split("\t")
            junc, baseSeq = infoList[0].split(":")[0], infoList[1]
            junc2seq[junc] = baseSeq
    d_tmp = [(v, k) for k, v in Counter(junc2seq.values()).iteritems()]
    d_tmp.sort(reverse=True)
    out = open(outFile, "w")
    for v, k in d_tmp:
        print >> out, "\t".join(map(str, [k, v]))
    out.close()

def splicesite2figure(refParams=None):
    refGenome = refParams.ref_genome
    splicesite2seq(refGenome, "IR.splicesite", outFile="IR.tmp")
    drawSSmotif("IR.tmp", "IR")
    splicesite2seq(refGenome, "IR.anno.splicesite", outFile="IR.anno.tmp")
    drawSSmotif("IR.anno.tmp", "IR.anno")
    splicesite2seq(refGenome, "IR.novel.splicesite", outFile="IR.novel.tmp")
    drawSSmotif("IR.novel.tmp", "IR.novel")
    splicesite2seq(refGenome, "SE.inc.splicesite", outFile="SE.inc.tmp")
    drawSSmotif("SE.inc.tmp", "SE.inc")
    splicesite2seq(refGenome, "SE.exc.splicesite", outFile="SE.exc.tmp")
    drawSSmotif("SE.exc.tmp", "SE.exc")
    splicesite2seq(refGenome, "A5SS.inc.splicesite", noThreeEnd=True, outFile="A5SS.inc.tmp")
    drawSSmotif("A5SS.inc.tmp", "A5SS.inc")
    splicesite2seq(refGenome, "A5SS.exc.splicesite", noThreeEnd=True, outFile="A5SS.exc.tmp")
    drawSSmotif("A5SS.exc.tmp", "A5SS.exc")
    splicesite2seq(refGenome, "A3SS.inc.splicesite", noFiveEnd=True, outFile="A3SS.inc.tmp")
    drawSSmotif("A3SS.inc.tmp", "A3SS.inc")
    splicesite2seq(refGenome, "A3SS.exc.splicesite", noFiveEnd=True, outFile="A3SS.exc.tmp")
    drawSSmotif("A3SS.exc.tmp", "A3SS.exc")

def getSpliceSite(asType=None, asFile=None, outFile=None):
    if asType == "IR":
        out = open(outFile, "w")
        with open(asFile) as f:
            for line in f:
                lineInfo = line.strip("\n").split(":")
                posList = lineInfo[2].split("-")
                print >>out, "\t".join([lineInfo[0], posList[0], posList[1], lineInfo[3]])
        out.close()
    elif asType == "SE":
        cmd = "seDecompose.pl        confident.SE.lst   >SE.inc.splicesite   2>SE.exc.splicesite"
        subprocess.call(cmd, shell=True)
    elif asType == "A3SS":
        cmd = "anssDecompose.pl -n 5 confident.A5SS.lst >A5SS.inc.splicesite 2>A5SS.exc.splicesite"
        subprocess.call(cmd, shell=True)
    elif asType == "A5SS":
        cmd = "anssDecompose.pl -n 3 confident.A3SS.lst >A3SS.inc.splicesite 2>A3SS.exc.splicesite"
        subprocess.call(cmd, shell=True)

def getDist2TTS(refParams=None, paGroup=None):
    with open(paGroup) as f:
        out = open("pbPA.bed6", "w")
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if lineInfo[5] == "+":
                print >>out, "\t".join(map(str, [lineInfo[0], int(lineInfo[2])-1, lineInfo[2]] + lineInfo[3:]))
            else:
                print >>out, "\t".join(map(str, [lineInfo[0], lineInfo[1], int(lineInfo[1])+1] + lineInfo[3:]))
        out.close()
    cmd = '''
        bedtools closest -a <(sort -k1,1 -k2,2n pbPA.bed6) -b <(gpeFeature.pl --tts {}|
        sort -k1,1 -k2,2n) -s -D a | select.pl -i 13,4 | sort -u | tee pbPA2TTS.tsv | 
        cut -f1 | box.R -ng -nJ -no -y='Distance to TTS' -p=pbPA2TTS.pdf
    '''.format(refParams.ref_gpe)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

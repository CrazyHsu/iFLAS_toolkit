#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: findAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-08-13 19:23:01
Last modified: 2019-08-13 19:23:02
'''
from commonFuncs import *
from commonObjs import Gene2Reads

def findAS(gene2ReadsDict, outASType=None, anno=True, out=None):
    findASType = {"IR": findIR, "SE": findSE, "A3SS": findA3SS, "A5SS": findA5SS}
    # if outASType and outASType.upper() not in findASType:
    #     raise Exception("You should input the right AS type!")
    if anno:
        if outASType:
            if outASType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[outASType](gene2ReadsDict, outASType, out)
        else:
            for asType in findASType:
                findASType[asType](gene2ReadsDict)
    else:
        newGene2ReadsDict = {}
        for myChr, myValue in gene2ReadsDict.iteritems():
            for myStrand in myValue:
                sortedStrandV = sorted(myValue[myStrand], key=lambda x: x.minpos)
                count = 1
                fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                gene2reads = Gene2Reads(fakeGeneName)
                gene2reads.updateFromGene2Reads(sortedStrandV[0])
                clusterEnd = sortedStrandV[0].maxpos
                for i in range(1, len(sortedStrandV)):
                    myStart, myEnd = sortedStrandV[i].minpos, sortedStrandV[i].maxpos
                    if myStart < clusterEnd:
                        gene2reads.updateFromGene2Reads(sortedStrandV[i])
                        if myEnd > clusterEnd: clusterEnd = myEnd
                    else:
                        newGene2ReadsDict[fakeGeneName] = gene2reads
                        count += 1
                        fakeGeneName = "%s:%s:%s_%d" % (myChr, myStrand, "Novel", count)
                        gene2reads = Gene2Reads(fakeGeneName)
                        gene2reads.updateFromGene2Reads(sortedStrandV[i])
                newGene2ReadsDict[fakeGeneName] = gene2reads
        if outASType:
            if outASType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                findASType[outASType](newGene2ReadsDict, outASType, out, anno=False)
        else:
            for asType in findASType:
                findASType[asType](newGene2ReadsDict)
        # for asType in findASType:
        #     findASType[asType](newGene2ReadsDict)
        gene2ReadsDict = newGene2ReadsDict
    for gene in gene2ReadsDict:
        gene2ReadsDict[gene].getCombinationEvents()
        gene2ReadsDict[gene].construcFakeRef()
    return gene2ReadsDict


def findIR(gene2ReadsDict, outAS=None, out=None, anno=True):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        irDict = {}
        for rName1 in rNames:
            exonStarts, exonEnds = readsDict[rName1].blockStarts, readsDict[rName1].blockEnds
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
                    exons2 = readsDict[rName2].readExonChain.split(";")
                    for exon2 in range(len(exons2)):
                        exon2Start, exon2End = [int(x) for x in exons2[exon2].split("-")]
                        if exon2Start < exonEnds[exon1] and exonStarts[exon1 + 1] < exon2End:
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
                            if irStart < transDict[tran][1] and irEnd > transDict[tran][0]:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            retentionReads = [x.readName for x in irDict[ir]["retention"]]
                            splicedReads = [x.readName for x in irDict[ir]["spliced"]]
                            psi = float(len(retentionReads))/(len(retentionReads)+len(splicedReads))*1000
                            print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName+":"+ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
                    else:
                        retentionReads = [x.readName for x in irDict[ir]["retention"]]
                        splicedReads = [x.readName for x in irDict[ir]["spliced"]]
                        psi = float(len(retentionReads)) / (len(retentionReads) + len(splicedReads)) * 1000
                        print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName + ":" + ir, psi, strand, gName, ",".join(retentionReads), len(retentionReads), ",".join(splicedReads), len(splicedReads)]))
        gene2ReadsDict[gName].asDict.update({"IR": newIrDict})


def findSE(gene2ReadsDict, outAS=None, out=None, anno=True):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        exonDict = {}
        for rName1 in rNames:
            exons = readsDict[rName1].readExonChain.split(";")
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
                        if re.search(skipPat, readsDict[rName2].readExonChain):
                            if "skip" not in exonDict[seChain]:
                                exonDict[seChain].__setitem__("skip", [readsDict[rName2]])
                            else:
                                exonDict[seChain]["skip"].append(readsDict[rName2])
                        elif re.search(keepPat, readsDict[rName2].readExonChain):
                            exonDict[seChain]["keep"].append(readsDict[rName2])
        newExonDict = {}
        for seChain1 in exonDict.keys():
            if "skip" in exonDict[seChain1]:
                newExonDict[seChain1] = exonDict[seChain1]
                if outAS:
                    keepReads = [x.readName for x in exonDict[seChain1]["keep"]]
                    skipReads = [x.readName for x in exonDict[seChain1]["skip"]]
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
                            if starts[0] < transDict[tran][1] and ends[-1] > transDict[tran][0]:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                            blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                            psi = float(keepReadsCount)/(keepReadsCount + skipReadsCount) * 1000
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
                        print >> out, "\t".join(map(str,
                                                    [chrom, str(starts[0]), str(ends[-1]), "%s:%s" % (gName, seChain1),
                                                     psi, strand, starts[0], ends[-1], "0,0,0", len(starts), blockSizes,
                                                     blockRelStarts, gName, lSpliceSite, rSpliceSite,
                                                     ",".join(keepReads), len(keepReads), ",".join(skipReads),
                                                     len(skipReads)]))
        gene2ReadsDict[gName].asDict.update({"SE": newExonDict})


def findA3SS(gene2ReadsDict, outAS=None, out=None, anno=True):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 1 if strand == "+" else 0
        for r1 in range(len(rNames) - 1):
            exonsChain = gene2ReadsDict[gName].reads[rNames[r1]].readExonChain
            exons1 = exonsChain.split(";")
            indexEnd1 = len(exons1) if strand == "+" else len(exons1) - 1
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].readExonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) if strand == "+" else len(exons2) - 1
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "-":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]):
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]):
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]):
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]):
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            if outAS:
                incReads = [x.readName for x in altDict[alt]["inc"]]
                excReads = [x.readName for x in altDict[alt]["exc"]]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount)]))
        gene2ReadsDict[gName].asDict.update({"A3SS": altDict})


def findA5SS(gene2ReadsDict, outAS=None, out=None, anno=True):
    for gName in gene2ReadsDict:
        readsDict = gene2ReadsDict[gName].reads
        chrom = gene2ReadsDict[gName].chrom
        strand = gene2ReadsDict[gName].strand
        rNames = readsDict.keys()
        if len(rNames) < 2: continue
        altDict = {}
        indexStart = 0 if strand == "+" else 1
        for r1 in range(len(rNames) - 1):
            exonsChain = gene2ReadsDict[gName].reads[rNames[r1]].readExonChain
            exons1 = exonsChain.split(";")
            indexEnd1 = len(exons1) - 1 if strand == "+" else len(exons1)
            for r2 in range(r1 + 1, len(rNames)):
                exonsChain2 = gene2ReadsDict[gName].reads[rNames[r2]].readExonChain
                exons2 = exonsChain2.split(";")
                indexEnd2 = len(exons2) - 1 if strand == "+" else len(exons2)
                for i in range(indexStart, indexEnd1):
                    exon1Start, exon1End = listStr2Int(exons1[i].split("-"))
                    for j in range(indexStart, indexEnd2):
                        exon2Start, exon2End = listStr2Int(exons2[j].split("-"))
                        if exon1Start < exon2End and exon1End > exon2Start:
                            if strand == "+":
                                if exon1End < exon2End and exon2End < int(exons1[i + 1].split("-")[0]):
                                    exonBoundChain = "%d-%d" % (exon1End, exon2End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                elif exon1End > exon2End and exon1End < int(exons2[j + 1].split("-")[0]):
                                    exonBoundChain = "%d-%d" % (exon2End, exon1End)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                            else:
                                if exon1Start < exon2Start and exon1Start > int(exons2[j - 1].split("-")[1]):
                                    exonBoundChain = "%d-%d" % (exon1Start, exon2Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r2]])
                                elif exon1Start > exon2Start and exon2Start > int(exons1[i - 1].split("-")[1]):
                                    exonBoundChain = "%d-%d" % (exon2Start, exon1Start)
                                    if exonBoundChain not in altDict:
                                        altDict[exonBoundChain] = {"inc": set(), "exc": set()}
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                                    else:
                                        altDict[exonBoundChain]["exc"].add(readsDict[rNames[r1]])
                                        altDict[exonBoundChain]["inc"].add(readsDict[rNames[r2]])
                        if exon2Start > exon1End: break
        for alt in altDict:
            altStart, altEnd = alt.split("-")
            altDict[alt]["inc"] = list(altDict[alt]["inc"])
            altDict[alt]["exc"] = list(altDict[alt]["exc"])
            if outAS:
                incReads = [x.readName for x in altDict[alt]["inc"]]
                excReads = [x.readName for x in altDict[alt]["exc"]]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount)]))
        gene2ReadsDict[gName].asDict.update({"A5SS": altDict})

def getAS(gpeDict, geneObj, geneDict, novelDict, offset=0):
    # geneDict, novelDict = {}, {}
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
                        consensusIntronN = getConsensusIntronNfromAS(readStarts, readEnds, transStarts, transEnds, offset)
                        if consensusIntronN >= 1:
                            transName, geneName = trans[2].transName, trans[2].geneName
                            if geneName not in geneDict:
                                gene2Reads = Gene2Reads(geneObjName)
                                gene2Reads.update(read)
                                gene2Reads.trans = {transName: trans}
                                geneDict[geneName] = gene2Reads
                            else:
                                gene2Reads = geneDict[geneName]
                                if read.readName not in [i.readName for i in gene2Reads.reads.values()]:
                                    gene2Reads.update(read)
                                gene2Reads.trans.update({transName: trans})
                                geneDict[geneName] = gene2Reads
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
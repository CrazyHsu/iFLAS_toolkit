#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: paCluster.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-02-24 16:53:00
Last modified: 2021-02-24 16:53:00
'''

from collections import Counter
import subprocess

class Bed12(object):
    "BED12 format gene structure."

    class BedError(Exception):
        "Error in manipulating Bed12 structures"

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    def __init__(self, line=None):
        if line:
            self.record = line.strip().split("\t")
            (self.chrom, self.chromStart, self.chromEnd, self.name,
             self.score, self.strand, self.thickStart, self.thickEnd,
             self.itemRgb, self.blockCount, self.blockSizes,
             self.blockStarts) = self.record[:12]

            self.chromStart = int(self.chromStart)
            self.chromEnd = int(self.chromEnd)
            self.score = int(float(self.score))
            self.thickStart = int(self.thickStart)
            self.thickEnd = int(self.thickEnd)
            self.blockCount = int(self.blockCount)

            self.blockSizes = self.blockSizes.strip(",").split(",")
            self.blockStarts = self.blockStarts.strip(",").split(",")

            assert len(self.blockStarts) == len(self.blockSizes)
            assert len(self.blockStarts) == self.blockCount

            for i in range(self.blockCount):
                self.blockSizes[i] = int(self.blockSizes[i])
                self.blockStarts[i] = int(self.blockStarts[i])

            self.exonStarts = []
            self.exonEnds = []
            for i in range(self.blockCount):
                self.exonStarts.append(self.chromStart + self.blockStarts[i])
                self.exonEnds.append(self.exonStarts[i] + self.blockSizes[i])

            self.exons = self.parse_exon()
            self.introns = self.parse_intron()
            self.exonChain = ';'.join(map(lambda x: str(self.exons[x][0])+"-"+str(self.exons[x][1]), range(len(self.exons))))
            self.juncChain = ";".join(map(lambda x: str(self.introns[x][0])+"-"+str(self.introns[x][1]), range(len(self.introns))))
        else:
            self.empty()

    def empty(self):
        "return an empty class with all values None, '', or []."
        self.chrom = ""
        self.chromStart = self.chromEnd = 0
        self.name = ""
        self.score = 0
        self.strand = ""
        self.thickStart = self.thickEnd = 0
        self.itemRgb = "0,0,0"
        self.blockCount = 0
        self.blockSizes = []
        self.blockStarts = []

    def __len__(self):
        "the length of transcript"
        return sum(self.blockSizes)

    def cds_len(self):
        "return cds/thick length"
        return sum([min(ed, self.thickEnd) - max(st, self.thickStart)
                    for (st, ed) in zip(self.exonStarts, self.exonEnds)
                    if (ed - self.thickStart) * (st - self.thickEnd) < 0])

    def parse_exon(self):
        "return a list of exon pos [(st, ed), (st, ed) , ... ]"
        exons = []
        for i in range(self.blockCount):
            st = self.exonStarts[i]
            ed = self.exonEnds[i]
            exons.append((st, ed))
        return exons

    def parse_intron(self):
        "return a list of intron pos [(st, ed], (st, ed], ... ]"
        introns = []
        for i in range(self.blockCount - 1):
            st = self.exonEnds[i]
            ed = self.exonStarts[i + 1]
            introns.append((st, ed))
        return introns

    def utr_5_len(self):
        "return the length of 5'UTR"
        if self.strand == "+":
            utr_5_len = sum([min(ed, self.thickStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.thickStart])
        else:
            utr_5_len = sum([ed - max(st, self.thickEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.thickEnd])
        return utr_5_len

    def utr_3_len(self):
        "return the length of 3'UTR"
        if self.strand == "-":
            utr_3_len = sum([min(ed, self.thickStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.thickStart])
        else:
            utr_3_len = sum([ed - max(st, self.thickEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.thickEnd])
        return utr_3_len

    def has_intron(self, intron):
        if intron.chrom == self.chrom and intron.strand == self.strand:
            for i in range(self.blockCount - 1):
                if (self.exonEnds[i] == intron.start and
                        self.exonStarts[i + 1] == intron.end):
                    return True
        return False

    def __repr__(self):
        "return a line of bed12 format, without newline ending"
        fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
                  self.name, str(self.score), self.strand, str(self.thickStart),
                  str(self.thickEnd), self.itemRgb, str(self.blockCount)]

        blockSizesline = ','.join(repr(i) for i in self.blockSizes)
        blockStartsline = ','.join(repr(i) for i in self.blockStarts)
        fields.extend((blockSizesline, blockStartsline))
        return "\t".join(fields)

    def go_to_gpe(self):
        "return a gpe with same structure"
        gpe = GenePredExtLine()
        gpe.bin = None
        gpe.transName = self.name
        gpe.chrom = self.chrom
        gpe.strand = self.strand
        gpe.txStart = self.chromStart
        gpe.txEnd = self.chromEnd
        gpe.cdsStart = self.thickStart
        gpe.cdsEnd = self.thickEnd
        gpe.exonCount = self.blockCount
        gpe.exonStarts = self.exonStarts
        gpe.exonEnds = self.exonEnds
        gpe.score = self.score
        gpe.geneName = "."
        gpe.cdsStartStat = "."
        gpe.cdsEndStat = "."
        gpe.exonFrames = ["."]
        return gpe

    def cds_bed(self):
        "return a new Bed12 object of cds region of self"
        assert self.thickStart < self.thickEnd
        newbed = Bed12(self.__repr__())
        newBlockSizes = []
        newBlockStarts = []
        newBlockCount = newbed.blockCount
        for i in range(newbed.blockCount):
            if newbed.thickEnd < newbed.exonStarts[i]:
                newBlockCount -= 1
            elif (newbed.thickStart <= newbed.exonStarts[i] and
                  newbed.thickEnd <= newbed.exonEnds[i]):
                newBlockSizes.append(newbed.thickEnd - newbed.exonStarts[i])
                newBlockStarts.append(newbed.exonStarts[i] - newbed.thickStart)
            elif (newbed.thickStart >= newbed.exonStarts[i] and
                  newbed.thickEnd <= newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.exonStarts[i])
                newBlockStarts.append(0)
            elif (newbed.thickStart >= newbed.exonStarts[i] and
                  newbed.thickEnd > newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.thickStart)
                newBlockStarts.append(0)
            elif (newbed.thickStart < newbed.exonStarts[i] and
                  newbed.thickEnd > newbed.exonEnds[i]):
                newBlockSizes.append(newbed.exonEnds[i] - newbed.exonStarts[i])
                newBlockStarts.append(newbed.exonStarts[i] - newbed.thickStart)
            elif newbed.thickStart > newbed.exonEnds[i]:
                newBlockCount -= 1
            else:
                raise self.BedError("Un-expected transcript structure.")
        assert len(newBlockSizes) == len(newBlockStarts) == newBlockCount
        newbed.blockSizes = newBlockSizes
        newbed.blockStarts = newBlockStarts
        newbed.blockCount = newBlockCount
        newbed.chromStart = newbed.thickStart
        newbed.chromEnd = newbed.thickEnd
        # renew exonStarts and exonEnds, in case further use
        return Bed12(newbed.__repr__())

# class Bed6Plus(Bed6):
#     def __init__(self, line=None):
#         if line:
#             Bed6.__init__(self, line)
#             self.otherList = self.record[6:]
#
#     def __str__(self):
#         return Bed6.__str__(self) + "\t" + "\t".join(self.otherList)

class Bed12Plus(Bed12):
    def __init__(self, line=None):
        if line:
            Bed12.__init__(self, line)
            self.otherList = self.record[12:]

    def __str__(self):
        return Bed12.__str__(self) + "\t" + "\t".join(self.otherList)

class BedFile(object):
    def __init__(self, bedFile, type=None):
        self.bedFile = bedFile
        self.reads = self.getReadsInfo(type)

    def getReadsInfo(self, type):
        readsDict = {}
        with open(self.bedFile) as f:
            for line in f:
                if type == "bed12":
                    b = Bed12(line)
                    readsDict.__setitem__(b.name, b)
                elif type == "bed12+":
                    b = Bed12Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict

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

getPaCluster(readsBed="Zm00001d027231.reads.assigned.unambi.bed12+", tofuGroup="tofu.collapsed.group.txt", threads=1, paClusterOut="test.bed8+")
paBed6 = open("test.PA.bed6+", "w")
with open("test.bed8+") as f:
    for line in f.readlines():
        lineInfo = line.strip("\n").split("\t")
        print >> paBed6, "\t".join(map(str, [lineInfo[0], lineInfo[6], lineInfo[7]] + lineInfo[3:6] + lineInfo[8:]))
paBed6.close()
# cmd = "paGroup.pl deSingleExonIsoform.bed12+ >isoform.paGrouped.tsv 2>isoform.paGrouped.bed6"
# subprocess.call(cmd, shell=True)
cmd = '''3endRevise.pl -p PA.bed6+ <(awk '{if($10>1){print}}' Zm00001d027231.reads.assigned.unambi.bed12+ | cut -f 1-12,15) | paGroup.pl >test.reads.paGrouped.tsv 2>test.reads.paGrouped.bed6'''
subprocess.call(cmd, shell=True, executable="/bin/bash")
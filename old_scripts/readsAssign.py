#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: readsAssign.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-04-29 21:07:14
Last modified: 2020-04-29 21:07:14
'''
import pybedtools, argparse, re
# from commonFuncs import *
# from commonObjs import *
def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]

def getBlockLength(blockList):
    return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))

def getOverlapOfTuple(tupleListA, tupleListB):
    res = []
    i = j = 0
    while i < len(tupleListA) and j < len(tupleListB):
        lo = max(tupleListA[i][0], tupleListB[j][0])
        hi = min(tupleListA[i][1], tupleListB[j][1])
        if lo <= hi:
            res.append([lo, hi])
        if tupleListA[i][1] < tupleListB[j][1]:
            i += 1
        else:
            j += 1
    return res

class GenePredExtLine(object):
    ' a line of gpe file'

    def __init__(self, line="", bincolumn=True):
        'initialize each field; attribute blockSizes and blockStarts as BED12.'
        if line:
            self.record = line.strip().split("\t")
            self.bin = None
            if bincolumn == True:
                self.bin = self.record.pop(0)
            self.transName = self.record[0]
            self.chrom = self.record[1]
            self.strand = self.record[2]
            self.txStart = int(self.record[3])
            self.txEnd = int(self.record[4])
            self.cdsStart = int(self.record[5])
            self.cdsEnd = int(self.record[6])
            self.exonCount = int(self.record[7])
            self.exonStarts = [int(i) for i in self.record[8].strip(',').split(',')]
            self.exonEnds = [int(i) for i in self.record[9].strip(',').split(',')]
            self.score = int(float(self.record[10]))
            self.geneName = self.record[11]
            self.cdsStartStat = self.record[12]
            self.cdsEndStat = self.record[13]
            self.exonFrames = [int(c) for c in
                               self.record[14].strip(",").split(",")]
            self.blockSizes = []
            self.blockStarts = []
            for i in range(self.exonCount):
                self.blockSizes.append(self.exonEnds[i] - self.exonStarts[i])
                self.blockStarts.append(self.exonStarts[i] - self.txStart)

            self.exons = self.parse_exon()
            self.introns = self.parse_intron()
        else:
            self.empty()

    def empty(self):
        "construct an empty gpe instance with all fields None, [], or 0"
        self.bin = None
        self.transName = ""
        self.chrom = ""
        self.strand = ""
        self.txStart = 0
        self.txEnd = 0
        self.cdsStart = 0
        self.cdsEnd = 0
        self.exonCount = 0
        self.exonStarts = []
        self.exonEnds = []
        self.score = 0
        self.geneName = ""
        self.cdsStartStat = ""
        self.cdsEndStat = ""
        self.exonFrames = []

    def __len__(self):
        "return total length of transcript"
        return sum([ed - st for st, ed in zip(self.exonStarts, self.exonEnds)])

    def copy(self):
        "return a new object of self"
        if self.bin:
            has_bin = True
        else:
            has_bin = False
        return GenePredExtLine(repr(self), bincolumn=has_bin)

    def cds_len(self):
        "return cds length"
        return sum([min(ed, self.cdsEnd) - max(st, self.cdsStart)
                    for (st, ed) in zip(self.exonStarts, self.exonEnds)
                    if (ed - self.cdsStart) * (st - self.cdsEnd) < 0])

    def utr_5_len(self):
        "return the length of 5'UTR"
        if self.strand == "+":
            utr_5_len = sum([min(ed, self.cdsStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.cdsStart])
        else:
            utr_5_len = sum([ed - max(st, self.cdsEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.cdsEnd])
        return utr_5_len

    def utr_3_len(self):
        "return the length of 3'UTR"
        assert self.is_standard()
        if self.strand == "-":
            utr_3_len = sum([min(ed, self.cdsStart) - st
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if st < self.cdsStart])
        else:
            utr_3_len = sum([ed - max(st, self.cdsEnd)
                             for st, ed in zip(self.exonStarts, self.exonEnds)
                             if ed > self.cdsEnd])
        return utr_3_len

    def is_standard(self):
        "check if all fields in gpe are standard (i.e. no abnormal positions)"
        "the standards might be modified to accommodate specific filters"
        if (self.txStart < self.txEnd and
                self.cdsStart < self.cdsEnd and
                self.exonCount > 0):
            return True
        else:
            return False

    def is_complete(self):
        "return true if cdsStartStat and cdsEndStat are cmpl, else False"
        if self.cdsStartStat == self.cdsEndStat == "cmpl":
            return True
        else:
            return False

    def parse_exon(self):
        "return a list of exon pos [(st, ed), (st, ed) , ... ]"
        exons = []
        for i in range(self.exonCount):
            st = self.exonStarts[i]
            ed = self.exonEnds[i]
            exons.append((st, ed))
        return exons

    def parse_intron(self):
        "return a list of intron pos [(st, ed], (st, ed], ... ]"
        introns = []
        for i in range(self.exonCount - 1):
            st = self.exonEnds[i]
            ed = self.exonStarts[i + 1]
            introns.append((st, ed))
        return introns

    def has_intron(self, intron):
        "determine if the argument is one of self's introns"
        if self.chrom == intron.chrom and self.strand == intron.strand:
            for pos in self.introns:
                if pos[0] == intron.st and pos[1] == intron.ed:
                    return True
        return False

    def go_to_bed(self, plus=False, gene=False):
        "return a Bed12 object of same gene structure."
        if plus == True or gene == True:
            bed = Bed12Plus()
        else:
            bed = Bed12()
        bed.chrom = self.chrom
        bed.chromStart = self.txStart
        bed.chromEnd = self.txEnd
        bed.name = self.transName
        bed.score = self.score
        bed.itemRgb = "0,0,0"
        bed.strand = self.strand
        bed.thickStart = self.cdsStart
        bed.thickEnd = self.cdsEnd
        bed.blockCount = self.exonCount
        bed.blockSizes = self.blockSizes
        bed.blockStarts = self.blockStarts
        bed.exonStarts = self.exonStarts
        bed.exonEnds = self.exonEnds
        if gene == True:
            bed.otherList = [self.geneName]
        return bed

    def __repr__(self):
        "return the line generating this gpe object without the last newline."
        outputlist = [self.transName, self.chrom, self.strand, repr(self.txStart),
                      repr(self.txEnd), repr(self.cdsStart), repr(self.cdsEnd),
                      repr(self.exonCount)]
        exonStarts_seq = ",".join([repr(i) for i in self.exonStarts])
        exonEnds_seq = ",".join([repr(i) for i in self.exonEnds])
        outputlist.extend((exonStarts_seq, exonEnds_seq))
        exonFrames_seq = ",".join([repr(i) for i in self.exonFrames])
        outputlist.extend((repr(self.score), self.geneName, self.cdsStartStat,
                           self.cdsEndStat, exonFrames_seq))
        if self.bin:
            return self.bin + "\t" + "\t".join(outputlist)
        else:
            return "\t".join(outputlist)

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

class Bed12Plus(Bed12):
    def __init__(self, line=None):
        if line:
            Bed12.__init__(self, line)
            self.otherList = self.record[12:]

    def __str__(self):
        return Bed12.__str__(self) + "\t" + "\t".join(self.otherList)

def getConsensusIntronN1(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    j, consensusN = 0, 0
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronStarts1[i]: break
            if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
                    intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                consensusN += 1
                j += 1
    return consensusN

def readsAssign(transBedFile, readsBedFile, offset=10, minConsenesusIntronN=1, minCoverageOnRead=0.9, singleLine=True, transColNum=13, readsColNum=13, outPrefix="readsAssign", group=True):
    transBedObj = pybedtools.BedTool(transBedFile)
    readsBedObj = pybedtools.BedTool(readsBedFile)
    intersectRes = readsBedObj.intersect(transBedObj, wa=True, wb=True)
    notintersectRes = readsBedObj.intersect(transBedObj, v=True)
    matchedDict = {}
    outList = []
    for i in intersectRes:
        infoList = str(i).strip("\n").split("\t")
        if readsColNum > 12:
            readsBed = Bed12Plus("\t".join(infoList[:readsColNum]))
        else:
            readsBed = Bed12("\t".join(infoList[:readsColNum]))
        if transColNum < 13:
            raise Exception("The colmun count of your transcript bed should large than 12")
        transBed = Bed12Plus("\t".join(infoList[readsColNum:readsColNum+transColNum]))
        readStart, readEnd, readsExons, transExons = readsBed.chromStart, readsBed.chromEnd, readsBed.exons, transBed.exons
        overlapLength = getBlockLength(getOverlapOfTuple(readsExons, transExons))
        coverageOnRead = overlapLength/float(getBlockLength(readsExons))
        coverageOnTrans = overlapLength/float(getBlockLength(transExons))
        transGene = transBed.record[12]
        if len(readsBed.exons) == 1:
            if transBed.exonStarts[0] < readEnd and readEnd <= transBed.exonEnds[0] or \
                    transBed.exonStarts[-1] <= readStart and readStart < transBed.exonEnds[-1]:
                inEndExonAndExtension = 1
            else:
                inEndExonAndExtension = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "matchedTransExonCount": [len(transBed.exons)],
                                              "inEndExonAndExtension": [inEndExonAndExtension], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedTransExonCount"].append(len(transBed.exons))
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["inEndExonAndExtension"].append(inEndExonAndExtension)
        else:
            consensusIntronN = getConsensusIntronN1(readsBed.exonStarts, readsBed.exonEnds, transBed.exonStarts, transBed.exonEnds, offset=offset)
            readsJuncChain, transJuncChain = readsBed.juncChain, transBed.juncChain
            if re.search(readsJuncChain + "$", transJuncChain) or re.search("^" + readsJuncChain, transJuncChain):
                juncChainFlank = 1
            else:
                juncChainFlank = 0

            if readsBed.name not in matchedDict:
                matchedDict[readsBed.name] = {"matchedTrans": [transBed.name], "matchedCovOnReads": [coverageOnRead],
                                              "matchedCovOnTrans": [coverageOnTrans], "matchedGenes": [transGene],
                                              "consensusIntronN": [consensusIntronN],
                                              "juncChainFlank": [juncChainFlank], "readsBed": readsBed}
            else:
                matchedDict[readsBed.name]["matchedTrans"].append(transBed.name)
                matchedDict[readsBed.name]["matchedCovOnReads"].append(coverageOnRead)
                matchedDict[readsBed.name]["matchedCovOnTrans"].append(coverageOnTrans)
                matchedDict[readsBed.name]["matchedGenes"].append(transGene)
                matchedDict[readsBed.name]["consensusIntronN"].append(consensusIntronN)
                matchedDict[readsBed.name]["juncChainFlank"].append(juncChainFlank)

    for readName in matchedDict:
        if singleLine:
            if "juncChainFlank" in matchedDict[readName]:
                readType = "I"
                readIntronN = len(matchedDict[readName]["readsBed"].introns)
                realMatchedTrans, realConsensusIntronN, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                        (readIntronN < minConsenesusIntronN and readIntronN == matchedDict[readName]["consensusIntronN"][j]) or \
                        (readIntronN >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN) or \
                        (matchedDict[readName]["matchedCovOnReads"][j] >= minCoverageOnRead):
                        readType = "E"
                        realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                        realConsensusIntronN.append(matchedDict[readName]["consensusIntronN"][j])
                        realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                        realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                        realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(realMatchedTrans), ",".join(realMatchedGene),
                                    ",".join(map(str, realConsensusIntronN)), ",".join(map(str, realMatchedCovOnReads)),
                                    ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]),
                                    ",".join(map(str, matchedDict[readName]["consensusIntronN"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
            else:
                readType = "I"
                realMatchedTrans, realMatchedCovOnReads, realMatchedCovOnTrans, realMatchedGene = [], [], [], []
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        continue
                    readType = "E"
                    realMatchedTrans.append(matchedDict[readName]["matchedTrans"][j])
                    realMatchedCovOnReads.append(matchedDict[readName]["matchedCovOnReads"][j])
                    realMatchedCovOnTrans.append(matchedDict[readName]["matchedCovOnTrans"][j])
                    realMatchedGene.append(matchedDict[readName]["matchedGenes"][j])
                if readType == "E":
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(map(str, realMatchedTrans)), ",".join(map(str, realMatchedGene)), "NA",
                                    ",".join(map(str, realMatchedCovOnReads)), ",".join(map(str, realMatchedCovOnTrans))])
                else:
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    ",".join(matchedDict[readName]["matchedTrans"]),
                                    ",".join(matchedDict[readName]["matchedGenes"]), "NA",
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnReads"])),
                                    ",".join(map(str, matchedDict[readName]["matchedCovOnTrans"]))])
        else:
            if "juncChainFlank" in matchedDict[readName]:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    readType = "I"
                    exonNum = len(matchedDict[readName]["readsBed"].exons)
                    if matchedDict[readName]["juncChainFlank"][j] == 1 or \
                            (exonNum < minConsenesusIntronN and exonNum - 1 == matchedDict[readName]["consensusIntronN"][j]) or \
                            (exonNum >= minConsenesusIntronN and matchedDict[readName]["consensusIntronN"][j] >= minConsenesusIntronN):
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j], matchedDict[readName]["matchedGenes"][j],
                                    matchedDict[readName]["consensusIntronN"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])
            else:
                for j in range(len(matchedDict[readName]["matchedTrans"])):
                    if matchedDict[readName]["matchedTransExonCount"][j] > 1 and \
                        matchedDict[readName]["inEndExonAndExtension"][j] == 0 and \
                        matchedDict[readName]["matchedCovOnReads"][j] < minCoverageOnRead:
                        readType = "I"
                    else:
                        readType = "E"
                    outList.append(["\t".join(matchedDict[readName]["readsBed"].record), readType,
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnReads"][j], "NA",
                                    matchedDict[readName]["matchedTrans"][j],
                                    matchedDict[readName]["matchedCovOnTrans"][j]])

    for i in notintersectRes:
        outList.append([str(i).strip("\n"), "IG", "NA", "NA", "NA", "NA", "NA"])
        # print str(i).strip("\n") + "\t" + "\t".join(["IG", "NA", "NA", "NA", "NA", "NA"])

    assignOut = open(outPrefix + ".bed12+", "w")
    sortedOutList = sorted(outList, key=lambda x: (x[0].split("\t")[0], int(x[0].split("\t")[1])))
    for item in sortedOutList:
        print >> assignOut, "\t".join(map(str, item))
    assignOut.close()

    # if group:
    #     unambiOut = open(outPrefix + ".unambi.bed12+", "w")
    #     ambiOut = open(outPrefix + ".ambiguous.bed12+", "w")
    #     novelList = []
    #     for item in sortedOutList:
    #         if item[3] == "NA" or item[1] != "E":
    #             novelList.append(item[0].split("\t"))
    #         else:
    #             uniqGenes = list(set(item[3].split(",")))
    #             if len(uniqGenes) == 1:
    #                 print >> unambiOut, item[0] + "\t" + uniqGenes[0]
    #             else:
    #                 print >> ambiOut, item[0] + "\t" + ",".join(uniqGenes)
    #     sortedNovelList = sorted(novelList, key=lambda x: int(x[1]))
    #     inc, clusterEnd = 0, 0
    #     for novelItem in sortedNovelList:
    #         myStart, myEnd = int(novelItem[1]), int(novelItem[2])
    #         if myStart < clusterEnd:
    #             if myEnd > clusterEnd:
    #                 clusterEnd = myEnd
    #         else:
    #             inc += 1
    #             clusterEnd = myEnd
    #         print >> unambiOut, "\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc]))
    #     unambiOut.close()
    #     ambiOut.close()

    if group:
        unambiOut = open(outPrefix + ".unambi.bed12+", "w")
        # ambiOut = open(outPrefix + ".ambiguous.bed12+", "w")
        novelList = []
        for item in sortedOutList:
            if item[3] == "NA" or item[1] != "E":
                novelList.append(item[0].split("\t"))
            else:
                uniqGenes = list(set(item[3].split(",")))
                if len(uniqGenes) == 1:
                    print >> unambiOut, item[0] + "\t" + uniqGenes[0]
                else:
                    novelList.append(item[0].split("\t"))
                    # print >> ambiOut, item[0] + "\t" + ",".join(uniqGenes)
        sortedNovelList = sorted(novelList, key=lambda x: int(x[1]))
        inc, clusterEnd = 0, 0
        for novelItem in sortedNovelList:
            myStart, myEnd = int(novelItem[1]), int(novelItem[2])
            if myStart < clusterEnd:
                if myEnd > clusterEnd:
                    clusterEnd = myEnd
            else:
                inc += 1
                clusterEnd = myEnd
            print >> unambiOut, "\t".join(novelItem) + "\t" + ":".join(map(str, [novelItem[0], novelItem[5], inc]))
        unambiOut.close()
        # ambiOut.close()

ref_bed = "/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/Zea_mays.B73_RefGen_v4.50.bed"
isoBed = "B73_merged.processed.ignore_id_removed.bed12+"
readsAssign(ref_bed, isoBed, readsColNum=13, outPrefix="ont_merged.reads.assigned", group=True)

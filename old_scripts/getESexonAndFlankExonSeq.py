#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: getESexonAndFlankExonSeq.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-04-08 00:46:51
Last modified: 2020-04-08 00:46:51
'''
import sys, pybedtools


def mergeIntervals(arr):
    # Sorting based on the increasing order
    # of the start intervals
    arr.sort(key=lambda x: x[0])

    m = []
    s = -10000
    max = -100000
    for i in range(len(arr)):
        a = arr[i]
        if a[0] > max:
            if i != 0:
                m.append([s, max])
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]

    if max != -100000 and [s, max] not in m:
        m.append([s, max])
    return m

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
            self.exonChain = ';'.join(map(lambda x: str(self.exonStarts[x])+"-"+str(self.exonEnds[x]), range(self.blockCount)))
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
                elif type == "bed6":
                    b = Bed6(line)
                    readsDict.__setitem__(b.name, b)
                elif type == "bed6+":
                    b = Bed6Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict


esFile = sys.argv[1]
readsFile = sys.argv[2]
refFasta = sys.argv[3]

esDict = {}
with open(esFile) as f:
    for i in f.readlines():
        infoList = i.strip("\n").split("\t")
        chrom = infoList[0]
        esName = infoList[3]
        esList = esName.split("@")[1].split(";")
        incReads = infoList[15].split(",")
        esDict[esName] = {"esList": esList, "incReads": incReads, "chrom": chrom}

reads = BedFile(readsFile, type="bed12+").reads
# Bed12().

esSurroundingInfo = {}
for i in esDict:
    chrom = esDict[i]["chrom"]
    readsObj = [reads[z] for z in esDict[i]["incReads"]]
    firstEsExon =  (int(esDict[i]["esList"][0].split("-")[0]), int(esDict[i]["esList"][0].split("-")[1]))
    lastEsExon = (int(esDict[i]["esList"][-1].split("-")[0]), int(esDict[i]["esList"][-1].split("-")[1]))
    leftExon = (0, 0)
    rightExon = (0, 0)
    tmpExons = readsObj[0].exons
    for x in readsObj[1:]:
        tmpExons = mergeIntervals(tmpExons + x.exons)

    esList = [[int(j.split("-")[0]), int(j.split("-")[1])] for j in esDict[i]["esList"]]
    bedList = []

    outFile = "{}.bed".format(i)
    out = open(outFile, "w")
    for j in tmpExons:
        if j in  esList:
            print >> out, "\t".join(map(str, [chrom, j[0], j[1], "esExon"]))
        else:
            print >> out, "\t".join(map(str, [chrom, j[0], j[1], "normal"]))
    out.close()
    bedObj = pybedtools.BedTool(outFile)
    bedRes = bedObj.sequence(refFasta, name=True)
    bedFastaOut = open("{}.bed.fasta".format(i), "w")
    for res in str(open(bedRes.seqfn).read()).split("\n")[:-1]:
        print >> bedFastaOut, res
    bedFastaOut.close()
    #     sortedReadsExons = sorted(x.exons)
    #     indexFirstEsExon = sortedReadsExons.index(firstEsExon)
    #     indexLastEsExon = sortedReadsExons.index(lastEsExon)
    #     if sortedReadsExons[indexFirstEsExon-1][1] - sortedReadsExons[indexFirstEsExon-1][0] > leftExon[1] - leftExon[0]:
    #         leftExon = sortedReadsExons[indexFirstEsExon-1]
    #     if sortedReadsExons[indexLastEsExon+1][1] - sortedReadsExons[indexLastEsExon+1][0] > rightExon[1] - rightExon[0]:
    #         rightExon = sortedReadsExons[indexLastEsExon+1]
    # # for j in esDict[i]["esList"]:
    # #     pass
    # # print i
    # outFile = "{}.bed".format(i)
    # out = open(outFile, "w")
    # print >> out, "\t".join(map(str, [chrom, leftExon[0], leftExon[1]]))
    # for j in esDict[i]["esList"]:
    #     print >> out, "\t".join(map(str, [chrom, j.split("-")[0], j.split("-")[1]]))
    # print >> out, "\t".join(map(str, [chrom, rightExon[0], rightExon[1]]))
    # out.close()
    # bedObj = pybedtools.BedTool(outFile)
    #
    # bedRes = bedObj.sequence(refFasta)
    # bedFastaOut = open("{}.bed.fasta".format(i), "w")
    # for res in str(open(bedRes.seqfn).read()).split("\n")[:-1]:
    #     print >> bedFastaOut, res
    # bedFastaOut.close()

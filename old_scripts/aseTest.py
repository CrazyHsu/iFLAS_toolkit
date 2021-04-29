#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: aseTest.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-04-02 20:51:23
Last modified: 2020-04-02 20:51:23
'''
from sys import maxint as MAXINT
import re, copy
import pybedtools
AStypes = ["SE", "IR", "A3SS", "A5SS"]

# from commonObjs import *
# from findAS import *
def getCnsSite(intronStarts, transIntronStarts, offset):
    '''
    get consensus intron boundary site between reference and reads info
    :return: count of consensus site num
    '''
    s2, cnsN = 0, 0
    for intronStart1 in intronStarts:
        for i in range(s2, len(transIntronStarts)):
            intronStart2 = transIntronStarts[i]
            if intronStart2 > intronStart1 + offset: break
            if intronStart1 - offset <= intronStart2 and intronStart2 <= intronStart1 + offset:
                cnsN += 1
                s2 += 1
    return cnsN

class Bed6(object):
    "BED6 format gene structure"

    class BedError(Exception):
        "Error in manipulating bed12 structures"

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    def __init__(self, line=None):
        if line:
            self.record = line.strip().split("\t")
            (self.chrom, self.chromStart, self.chromEnd, self.name,
             self.score, self.strand) = self.record[:6]
            self.chromStart = int(self.chromStart)
            self.chromEnd = int(self.chromEnd)
            try:
                self.score = int(float(self.score))
            except ValueError:
                pass
            self.pos = (str(self.chrom) + ":" + str(self.chromStart) +
                        "-" + str(self.chromEnd))
        else:
            self.empty()

    def empty(self):
        (self.chrom, self.chromStart, self.chromEnd, self.name,
         self.score, self.strand) = ("", 0, 0, "", 0, "")

    def toBed12(self):
        line = "\t".join([repr(self), repr(self.chromStart),
                          repr(self.chromEnd), "255,0,0", "1",
                          repr(self.chromEnd - self.chromStart), "0"])
        return Bed12(line)

    def __repr__(self):
        "return a line of bed6 format, without newline ending"
        fields = [self.chrom, str(self.chromStart), str(self.chromEnd),
                  self.name, str(self.score), self.strand]
        return "\t".join(fields)


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

class GenePredObj(object):
    """
    build gene prediction object from gpe file
    """
    def __init__(self, gpeFile, bincolumn=False):
        self.gpeFile = gpeFile
        self.bincolumn = bincolumn
        self.geneName2gpeObj = {}
        self.genePredDict = self.buildObj()

    def buildObj(self, geneId=None):
        tmpDict = {}
        with open(self.gpeFile) as f:
            for i in f:
                gpeObj = GenePredExtLine(i, bincolumn=self.bincolumn)
                chrom, strand, start, end = gpeObj.chrom, gpeObj.strand, gpeObj.txStart, gpeObj.txEnd
                if geneId != None:
                    gpeObj.transName = gpeObj.transName.replace(gpeObj.geneName, geneId)
                    gpeObj.geneName = geneId
                # if geneId != None:
                #     gpeObj.geneName = geneId
                geneName = gpeObj.geneName

                if geneName not in self.geneName2gpeObj:
                    self.geneName2gpeObj[geneName] = [gpeObj]
                else:
                    self.geneName2gpeObj[geneName].append(gpeObj)

                if chrom not in tmpDict:
                    tmpDict[chrom] = {strand: [[start, end, gpeObj]]}
                elif strand not in tmpDict[chrom]:
                    tmpDict[chrom][strand] = [[start, end, gpeObj]]
                else:
                    tmpDict[chrom][strand].append([start, end, gpeObj])
            for k, v in tmpDict.iteritems():
                for s in v:
                    v[s] = sorted(v[s], key=lambda x: (x[0], x[1]))
        return tmpDict

    # def buildGenePredObj(self):
    #     return self.buildObj()

    def changeGeneId(self, geneId):
        self.geneName2gpeObj = {}
        self.genePredDict = self.buildObj(geneId=geneId)

    def gpeObj2file(self, fout, geneName):
        f = open(fout, "w")
        gps = self.geneName2gpeObj[geneName]
        for i in gps:
            print >>f, i
        f.close()

    def getGeneExonLength(self):
        exonLength = 0
        for g in self.geneName2gpeObj:
            exonList = []
            for l in self.geneName2gpeObj[g]:
                exonList.extend(l.exons)
            exonList = set(exonList)
            exonLength += getBlockLength(exonList)
        return exonLength

    def getGeneIntronLength(self):
        intronLength = 0
        for g in self.geneName2gpeObj:
            intronList = []
            for l in self.geneName2gpeObj[g]:
                intronList.extend(l.introns)
            intronList = set(intronList)
            intronLength += getBlockLength(intronList)
        return intronLength

    def redefineBorder(self, specificMinpos, specificMaxpos, outFile="tmp.gpe"):
        out = open(outFile, "w")
        for i in self.geneName2gpeObj:
            for j in self.geneName2gpeObj[i]:
                records = copy.deepcopy(j.record)
                if j.txStart < specificMinpos:
                    for z in range(len(j.exonEnds)):
                        if j.exonStarts[z] >= specificMinpos and j.exonEnds[z-1] <= specificMaxpos:
                            records[3] = records[5] = str(j.exonStarts[z])
                            records[7] = str(j.exonCount - z)
                            records[8] = ",".join(map(str, j.exonStarts[z:]))
                            records[9] = ",".join(map(str, j.exonEnds[z:]))
                        elif j.exonStarts[z] <= specificMinpos and specificMinpos < j.exonEnds[z]:
                            records[3] = records[5] = str(specificMinpos)
                            records[7] = str(j.exonCount - z)
                            records[8] = ",".join(map(str, [specificMinpos] + j.exonStarts[z+1:]))
                            records[9] = ",".join(map(str, j.exonEnds[z:]))
                if j.txEnd > specificMaxpos:
                    for z in range(len(j.exonEnds)):
                        if j.exonEnds[z-1] <= specificMaxpos and j.exonStarts[z] >= specificMaxpos:
                            records[4] = records[6] = str(j.exonEnds[z])
                            records[7] = str(j.exonCount - len(j.exonStarts) + z + 1)
                            records[8] = ",".join(map(str, j.exonStarts[:z]))
                            records[9] = ",".join(map(str, j.exonEnds[:z]))
                        elif j.exonStarts[z] <= specificMaxpos and specificMaxpos < j.exonEnds[z]:
                            records[4] = records[6] = str(specificMaxpos)
                            records[7] = str(j.exonCount - len(j.exonStarts) + z + 1)
                            records[8] = ",".join(map(str, j.exonStarts[:z+1]))
                            records[9] = ",".join(map(str, j.exonEnds[:z] + [specificMaxpos]))
                recordStr = "\t".join(records)
                gpeLine = GenePredExtLine(recordStr, bincolumn=False)
                print >>out, str(gpeLine)
        out.close()
        return outFile

    def toBed(self, plus=False, gene=False, outFile=None):
        if outFile:
            out = open(outFile, "w")
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bed = gp[2].go_to_bed()
                        print >>out, bed
            out.close()
        else:
            bedList = []
            for chrom in self.genePredDict:
                for strand in self.genePredDict[chrom]:
                    for gp in self.genePredDict[chrom][strand]:
                        bedList.append(str(gp[2].go_to_bed(gene=gene)))
            return bedList

def getSizes(starts, ends):
    return map(lambda x: int(ends[x]) - int(starts[x]), range(len(starts)))

def getRelStarts(blockStarts):
    return map(lambda x: int(blockStarts[x]) - int(blockStarts[0]), range(len(blockStarts)))

def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]

def getBlockLength(blockList):
    return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))

def getFileRowCounts(inFile):
    return len(open(inFile).readlines())

def getColFromFile(myFile, col, sep=None):
    with open(myFile) as f:
        myDict = {}
        for line in f:
            lineInfo = line.strip("\n").split("\t")
            if sep:
                myDict.fromkeys(lineInfo[col - 1].split(sep), "")
            else:
                myDict[lineInfo[col - 1]] = ""
        return myDict

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

class Bed6Plus(Bed6):
    def __init__(self, line=None):
        if line:
            Bed6.__init__(self, line)
            self.otherList = self.record[6:]

    def __str__(self):
        return Bed6.__str__(self) + "\t" + "\t".join(self.otherList)

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

def getExonFramesForGPE(gpeObj):
    exonFrame, exonFrames = 0, []
    if gpeObj.strand == "+":
        for i in range(gpeObj.exonCount):
            if gpeObj.exonStarts[i] < gpeObj.cdsStart:
                exonFrame = -1 if gpeObj.exonEnds[i] < gpeObj.cdsStart else 0
            elif gpeObj.exonStarts[i] == gpeObj.cdsStart:
                exonFrame = 0
            else:
                if gpeObj.exonStarts[i] < gpeObj.cdsEnd:
                    if gpeObj.exonStarts[i - 1] < gpeObj.cdsStart:
                        exonFrame = (gpeObj.exonEnds[i - 1] - gpeObj.cdsStart) % 3
                    else:
                        exonFrame = (gpeObj.exonEnds[i - 1] - gpeObj.exonStarts[i - 1] + exonFrames[i - 1] - 3) % 3
                else:
                    exonFrame = -1
            exonFrames.append(exonFrame)
    else:
        for i in range(gpeObj.exonCount - 1, -1, -1):
            if gpeObj.exonEnds[i] > gpeObj.cdsEnd:
                exonFrame = -1 if gpeObj.exonStarts[i] > gpeObj.cdsEnd else 0
            elif gpeObj.exonEnds[i] == gpeObj.cdsEnd:
                exonFrame = -1 if gpeObj.cdsStart == gpeObj.cdsEnd else 0
            else:
                if gpeObj.exonEnds[i] > gpeObj.cdsStart:
                    if gpeObj.exonEnds[i + 1] > gpeObj.cdsEnd:
                        exonFrame = (gpeObj.cdsEnd - gpeObj.exonStarts[i + 1]) % 3
                    else:
                        exonFrame = (gpeObj.exonEnds[i + 1] - gpeObj.exonStarts[i + 1] + exonFrames[
                            gpeObj.exonCount - i - 2] - 3) % 3
                else:
                    exonFrame = -1
            exonFrames.append(exonFrame)
    return exonFrames

def getAbsStarts(start, blockRelStarts):
    return map(lambda x: start + x, blockRelStarts)

def sortTupleList(tupleList):
    return sorted(tupleList)

def mergeTupleList(tupleList):
    sortedTupleList = sortTupleList(tupleList)
    newTupleList = sortedTupleList[0:1]
    for myTuple in sortedTupleList[1:]:
        lastTuple = newTupleList.pop()
        if lastTuple[1] >= myTuple[0]:
            newTuple = (lastTuple[0], max(myTuple[1], lastTuple[1]))
            newTupleList.append(newTuple)
        else:
            newTupleList.extend([lastTuple, myTuple])
    return newTupleList

def getAbsLoc(start, blockSizes, blockRelStarts):
    '''
    get absolute location of the exons in transcripts through: start(first base location), blockSizes(length of each \
    exon), blockRelStarts(the relative position related to the first base, i.e. start)
    :returns: blockStarts and blockEnds, the list of starts of each starts and ends
    '''
    blockStarts = getAbsStarts(start, blockRelStarts)
    blockEnds = map(lambda x: blockStarts[x] + blockSizes[x], range(len(blockSizes)))
    return blockStarts, blockEnds

class ReadLineStruc(Bed12):
    """
    get Bed12+ read structure
    """
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.record = line.strip().split("\t")
        # self.chrom = self.record[0]
        # self.start, self.end = int(self.record[1]), int(self.record[2])
        # self.readName = self.record[3]
        # self.score = float(self.record[4])
        # self.strand = self.record[5]
        # self.blockNum = int(self.record[9])
        # self.blockSizes = [int(i) for i in self.record[10].strip(',').split(',')]
        # self.blockRelStarts = [int(i) for i in self.record[11].strip(',').split(',')]
        # self.blockStarts, self.blockEnds = getAbsLoc(self.start, self.blockSizes, self.blockRelStarts)
        self.exonLen = sum(self.blockSizes)
        self.readLen = abs(self.chromStart - self.chromEnd)
        # self.readExonList = map(lambda x: (self.blockStarts[x], self.blockEnds[x]), range(self.blockNum))
        self.geneName = self.record[12]

    def go_to_gpe(self):
        gpe = GenePredExtLine()
        gpe.bin = None
        gpe.transName = self.name
        gpe.chrom = self.chrom
        gpe.strand = self.strand
        gpe.txStart = self.thickStart
        gpe.txEnd = self.thickEnd
        gpe.cdsStart = self.chromStart
        gpe.cdsEnd = self.chromEnd
        gpe.exonCount = self.blockCount
        gpe.exonStarts = self.exonStarts
        gpe.exonEnds = self.exonEnds
        gpe.score = 0
        gpe.geneName = self.geneName
        gpe.cdsStartStat = "unk"
        gpe.cdsEndStat = "unk"
        gpe.exonFrames = getExonFramesForGPE(gpe)
        return gpe

    def __str__(self):
        return "%s:%s:%d-%d (%s)" % (self.name, self.chrom, self.chromStart, self.chromEnd, self.geneName)

    __repr__ = __str__

class Gene2Reads(object):
    def __init__(self, geneName):
        self.minpos = MAXINT
        self.maxpos = 0
        self.geneName = geneName
        self.readNames = []
        self.trans = {}
        self.reads = {}
        self.asDict = dict.fromkeys(AStypes, {})

    def update(self, readStruc):
        self.chrom = readStruc.chrom
        self.strand = readStruc.strand
        self.minpos = min(self.minpos, readStruc.chromStart)
        self.maxpos = max(self.maxpos, readStruc.chromEnd)
        self.readNames.append(readStruc.name)
        self.reads.update({readStruc.name: readStruc})

    def updateFromGene2Reads(self, gene2ReadsObj):
        self.chrom = gene2ReadsObj.chrom
        self.strand = gene2ReadsObj.strand
        self.minpos = min(self.minpos, gene2ReadsObj.minpos)
        self.maxpos = max(self.maxpos, gene2ReadsObj.maxpos)
        self.trans.update(gene2ReadsObj.trans)
        self.reads.update(gene2ReadsObj.reads)
        self.readNames.append(gene2ReadsObj.readNames)

    def getCombinationEvents(self):
        asType2asKey = {"A5SS": "exc", "A3SS": "exc", "SE": "skip", "IR": "spliced"}
        combineEvent = "combineEvent"
        read2AS = {}
        self.asDict[combineEvent] = {}
        for asType in self.asDict:
            for event in self.asDict[asType]:
                if asType2asKey[asType] in self.asDict[asType][event]:
                    for read in self.asDict[asType][event][asType2asKey[asType]]:
                        if read not in read2AS:
                            read2AS[read] = {asType: [event]}
                        elif asType not in read2AS[read]:
                            read2AS[read].update({asType: [event]})
                        else:
                            read2AS[read][asType].append(event)
        for r in read2AS:
            if len(read2AS[r]) == 1: continue
            for t, l in read2AS[r].items():
                for e in l:
                    self.asDict[t][e][asType2asKey[t]].remove(r)
                    if r not in self.asDict[combineEvent]:
                        self.asDict[combineEvent][r] = {t: [e]}
                    elif t not in self.asDict[combineEvent][r]:
                        self.asDict[combineEvent][r].update({t: [e]})
                    else:
                        self.asDict[combineEvent][r][t].append(e)
                    # self.asDict[combineEvent].add(r)
                    if len(self.asDict[t][e][asType2asKey[t]]) == 0:
                        del self.asDict[t][e]

    def construcFakeRef(self):
        asType2asKey = {"A5SS": "exc", "A3SS": "exc", "SE": "skip", "IR": "spliced"}
        for asType in asType2asKey:
            for event in self.asDict[asType]:
                reads = self.asDict[asType][event][asType2asKey[asType]]
                maxLen = reads[0].exonLen
                maxLenRead = reads[0]
                for i in range(1, len(reads)):
                    if reads[i].exonLen > maxLen:
                        maxLen = reads[i].exonLen
                        maxLenRead = reads[i]
                records = copy.deepcopy(maxLenRead.record)
                if asType == "SE":
                    excExonList = map(lambda x: (int(x.split("-")[0]), int(x.split("-")[1])), event.split("@")[1].split(";"))
                    excExonStarts = map(lambda x: x[0], excExonList)
                    excExonEnds = map(lambda x: x[1], excExonList)
                    newExonStarts = sorted(maxLenRead.blockStarts + excExonStarts)
                    newExonEnds = sorted(maxLenRead.blockEnds + excExonEnds)
                    # newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
                    # newExonRelStarts = map(lambda x: newExonEnds[x] - maxLenRead.start, range(len(newExonStarts)))
                    # records[3] = records[3]+"_"+asType+"_fake"
                    # records[9] = len(newExonSizes)
                    # records[10] = newExonRelStarts
                    # recordLine = "\t".join(records)
                    # fakeRead = ReadLineStruc(recordLine)
                else:
                    # records = maxLenRead.records
                    excStart, excEnd = map(int, re.split("-", event))
                    mergedExonList = mergeTupleList(maxLenRead.readExonList + [(excStart, excEnd)])
                    newExonStarts = map(lambda x: x[0], mergedExonList)
                    newExonEnds = map(lambda x: x[1], mergedExonList)
                newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
                newExonRelStarts = map(lambda x: newExonStarts[x] - maxLenRead.start, range(len(newExonStarts)))
                records[3] = records[3] + "_" + asType + "_fake"
                records[9] = len(newExonSizes)
                records[10] = ",".join(map(str, newExonSizes))
                records[11] = ",".join(map(str, newExonRelStarts))
                recordLine = "\t".join(map(str, records))
                fakeRead = ReadLineStruc(recordLine)
                self.asDict[asType][event][asType2asKey[asType]].append(fakeRead)
        fakeReadDict = {}
        for r in self.asDict["combineEvent"]:
            newExonStarts = r.blockStarts
            newExonEnds = r.blockStarts
            records = copy.deepcopy(r.record)
            newExonList = []
            for t in self.asDict["combineEvent"][r]:
                for e in self.asDict["combineEvent"][r][t]:
                    if t == "SE":
                        excExonList = map(lambda x: (int(x.split("-")[0]), int(x.split("-")[1])),
                                          e.split("@")[1].split(";"))
                        newExonList.extend(excExonList)
                        # excExonStarts = map(lambda x: int(x.split("-")[0]), e.split("@")[1].split(";"))
                        # excExonEnds = map(lambda x: int(x.split("-")[1]), e.split("@")[1].split(";"))
                        # newExonStarts = sorted(newExonStarts + excExonStarts)
                        # newExonEnds = sorted(newExonEnds + excExonEnds)
                        # newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
                        # newExonRelStarts = map(lambda x: newExonEnds[x] - maxLenRead.start, range(len(newExonStarts)))
                    else:
                        excStart, excEnd = map(int, re.split("-", e))
                        newExonList.append((excStart, excEnd))
            mergedExonList = mergeTupleList(r.readExonList + newExonList)
            newExonStarts = map(lambda x: x[0], mergedExonList)
            newExonEnds = map(lambda x: x[1], mergedExonList)
            newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
            newExonRelStarts = map(lambda x: newExonStarts[x] - r.start, range(len(newExonStarts)))
            records[3] = records[3] + "_combineEvent_fake11"
            records[9] = len(newExonSizes)
            records[10] = ",".join(map(str, newExonSizes))
            records[11] = ",".join(map(str, newExonRelStarts))
            recordLine = "\t".join(map(str, records))
            fakeRead = ReadLineStruc(recordLine)
            fakeReadDict[fakeRead] = {}
        self.asDict["combineEvent"].update(fakeReadDict)

    def __str__(self):
        return "%s:%s-%s (%s) (%s) (%d reads)" % \
               (self.chrom, self.minpos, self.maxpos, self.geneName, self.strand, len(self.reads))

    __repr__ = __str__

class Position(object):
    "chr:st-ed, zero based"

    def __init__(self, line):
        self.chrom, st_ed = line.strip().split(":")
        self.chromStart, self.chromEnd = (int(c) for c in st_ed.split("-"))
        assert self.chromStart < self.chromEnd

    def __repr__(self):
        return "%s:%d-%d" % (self.chrom, self.chromStart, self.chromEnd)

def findAS_bak(gene2readsDict, outASType=None, anno=True, out=None):
    findASType = {"IR": findIR, "SE": findSE, "A3SS": findA3SS, "A5SS": findA5SS}
    if anno:
        pass
    else:
        pass
    pass


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
        # outHand.close()
        if outASType:
            if outASType.upper() not in findASType:
                raise Exception("You should input the right AS type!")
            else:
                print "prepare to novel"
                # outHand = open("tmp.irtest2", "w")
                # for gg in newGene2ReadsDict:
                #     for rr in newGene2ReadsDict[gg].reads:
                #         print >>outHand, gg,"\t", rr, "\t","\t".join(newGene2ReadsDict[gg].reads[rr].record)
                # outHand.close()
                findASType[outASType](newGene2ReadsDict, outASType, out, anno=False)
        else:
            for asType in findASType:
                findASType[asType](newGene2ReadsDict)
        # for asType in findASType:
        #     findASType[asType](newGene2ReadsDict)
    #     gene2ReadsDict = newGene2ReadsDict
    # for gene in gene2ReadsDict:
    #     gene2ReadsDict[gene].getCombinationEvents()
    #     gene2ReadsDict[gene].construcFakeRef()
    # return gene2ReadsDict


def findIR(gene2ReadsDict, outAS=None, out=None, anno=True):
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
                            if irStart < transDict[tran].chromEnd and irEnd > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            retentionReads = [x.name for x in irDict[ir]["retention"]]
                            splicedReads = [x.name for x in irDict[ir]["spliced"]]
                            psi = float(len(retentionReads))/(len(retentionReads)+len(splicedReads))*1000
                            print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName + ":" + ir, psi, strand, gName,
                                                              ",".join(retentionReads), len(retentionReads),
                                                              ",".join(splicedReads), len(splicedReads)]))
                    else:
                        retentionReads = [x.name for x in irDict[ir]["retention"]]
                        splicedReads = [x.name for x in irDict[ir]["spliced"]]
                        psi = float(len(retentionReads))/(len(retentionReads)+len(splicedReads))*1000
                        print >> out, "\t".join(map(str, [chrom, irStart, irEnd, gName + ":" + ir, psi, strand, gName,
                                                          ",".join(retentionReads), len(retentionReads),
                                                          ",".join(splicedReads), len(splicedReads)]))
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
                        transDict = gene2ReadsDict[gName].trans
                        overlapWithGene = 0
                        for tran in transDict:
                            if starts[0] < transDict[tran].chromEnd and ends[-1] > transDict[tran].chromStart:
                                overlapWithGene = 1
                                break
                        if overlapWithGene == 1:
                            blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                            blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                            psi = float(keepReadsCount)/(keepReadsCount + skipReadsCount) * 1000
                            print >>out, "\t".join(map(str, [chrom, str(starts[0]), str(ends[-1]), \
                                                             "%s:%s" % (gName, seChain1), psi, strand, \
                                                             starts[0], ends[-1], "0,0,0", len(starts), blockSizes, blockRelStarts, \
                                                             gName, lSpliceSite, rSpliceSite, ",".join(keepReads), len(keepReads), \
                                                             ",".join(skipReads), len(skipReads)]))
                    else:
                        blockSizes = ",".join([str(x) for x in getSizes(starts, ends)])
                        blockRelStarts = ",".join([str(x) for x in getRelStarts(starts)])
                        psi = float(keepReadsCount) / (keepReadsCount + skipReadsCount) * 1000
                        print >> out, "\t".join(map(str, [chrom, str(starts[0]), str(ends[-1]), \
                                                          "%s:%s" % (gName, seChain1), psi, strand, \
                                                          starts[0], ends[-1], "0,0,0", len(starts), blockSizes,
                                                          blockRelStarts, \
                                                          gName, lSpliceSite, rSpliceSite, ",".join(keepReads),
                                                          len(keepReads), \
                                                          ",".join(skipReads), len(skipReads)]))
        gene2ReadsDict[gName].asDict.update({"SE": newExonDict})

def listStr2Int(myList):
    return [int(i) for i in myList]

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
            exonsChain = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain.split(";")
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
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
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
            exonsChain = gene2ReadsDict[gName].reads[rNames[r1]].exonChain
            exons1 = exonsChain.split(";")
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
                incReads = [x.name for x in altDict[alt]["inc"]]
                excReads = [x.name for x in altDict[alt]["exc"]]
                incReadsCount, excReadsCount = len(incReads), len(excReads)
                psi = float(incReadsCount)/(incReadsCount + excReadsCount) * 1000
                print >>out, "\t".join(map(str, [chrom, str(altStart), str(altEnd), \
                                          "%s:%s" % (gName, alt), str(psi), strand, gName, ",".join(incReads), \
                                          str(incReadsCount), ",".join(excReads), str(excReadsCount)]))
        gene2ReadsDict[gName].asDict.update({"A5SS": altDict})

def getConsensusIntronN1(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    totalCnsN = 0
    for asType in AStypes:
        if asType in ["IR", "SE"]:
            j, consensusN = 0, 0
            for i in range(len(intronStarts1)):
                for k in range(j, len(intronStarts2)):
                    if intronStarts2[k] > intronStarts1[i]: break
                    if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
                            intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                        consensusN += 1
                        j += 1
            totalCnsN += consensusN
        elif asType in ["A3SS", "A5SS"]:
            cnsIntronStartN = getCnsSite(intronStarts1, intronStarts2, offset)
            cnsIntronEndN = getCnsSite(intronEnds1, intronEnds2, offset)
            totalCnsN += cnsIntronStartN + cnsIntronEndN
    return totalCnsN

def getConsensusIntronN(exonStarts1, exonEnds1, exonStarts2, exonEnds2, asType, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    if asType in ["IR", "SE"]:
        j, consensusN = 0, 0
        for i in range(len(intronStarts1)):
            for k in range(j, len(intronStarts2)):
                if intronStarts2[k] > intronStarts1[i]: break
                if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
                        intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                    consensusN += 1
                    j += 1
        return consensusN
    elif asType in ["A3SS", "A5SS"]:
        cnsIntronStartN = getCnsSite(intronStarts1, intronStarts2, offset)
        cnsIntronEndN = getCnsSite(intronEnds1, intronEnds2, offset)
        return cnsIntronStartN + cnsIntronEndN


def getAS_bak(annoBedRes, novelBedRes, asType, offset=0):
    annoDict, novelDict = {}, {}
    isNovelDict = {}
    for line in annoBedRes:
        record = str(line).strip("\n").split("\t")
        read = Bed12("\t".join(record[:12]))
        trans = ReadLineStruc("\t".join(record[12:]))
        if read.name not in isNovelDict:
            isNovelDict[read.name] = {"read": read, "isNovel": 1}
        transStarts, transEnds = trans.exonStarts, trans.exonEnds
        readStarts, readEnds, blockCount = read.exonStarts, read.exonEnds, read.blockCount
        # if blockCount < 2: continue
        consensusIntronN = getConsensusIntronN(readStarts, readEnds, transStarts, transEnds, asType, offset)
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
                        consensusIntronN = getConsensusIntronN1(readStarts, readEnds, transStarts, transEnds, offset)
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

# def getAS_anno(gpeFile=None, sourceFile=None):
#     transBedList = GenePredObj(gpeFile, bincolumn=False).toBed(gene=True)
#     transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
#     readsBedList = []
#     with open(sourceFile) as f:
#         for line in f.readlines():
#             readStruc = Bed12(line)
#             readsBedList.append("\t".join(readStruc.record[:12]))
#     readsBedObj = pybedtools.BedTool("\n".join(readsBedList), from_string=True)
#     annoBedRes = readsBedObj.intersect(transBedObj, wa=True, wb=True, s=True)
#     novelBedRes = readsBedObj.intersect(transBedObj, v=True, s=True)
#     annoDict, novelDict = getAS_bak(annoBedRes, novelBedRes, asType, offset=0)
#     return annoDict, novelDict

def findASmain_bak(asType="IR", gpeFile=None, sourceFile=None, outFile=None):
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
    annoDict, novelDict = getAS_bak(annoBedRes, novelBedRes, asType, offset=0)
    out = open(outFile, "w")
    findAS(annoDict, outASType=asType, anno=True, out=out)
    findAS(novelDict, outASType=asType, anno=False, out=out)
    out.close()

# annoDict, novelDict = getAS_anno(gpeFile="/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.gpe", sourceFile="/data/CrazyHsu_data/Projects/isoseq/wangLab_data_test/wangLab_data_test_20200326/Zheng58/8/ASE/deSingleExonRead.bed12+")
# findASmain_bak(asType="IR", annoDict=annoDict, novelDict=novelDict, outFile="ir.test")

findASmain_bak(asType="IR", gpeFile="/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.gpe", sourceFile="/data/CrazyHsu_data/Projects/isoseq/wangLab_data_test/wangLab_data_test_20200326/Zheng58/8/ASE/deSingleExonRead.bed12+", outFile="ir.test")
# findASmain_bak(asType="SE", gpeFile="/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.gpe", sourceFile="/data/CrazyHsu_data/Projects/isoseq/wangLab_data_test/wangLab_data_test_20200326/Zheng58/8/ASE/deSingleExonRead.bed12+", outFile="se.test")
# findASmain_bak(asType="A5SS", gpeFile="/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.gpe", sourceFile="/data/CrazyHsu_data/Projects/isoseq/wangLab_data_test/wangLab_data_test_20200326/Zheng58/8/ASE/deSingleExonRead.bed12+", outFile="a5ss.test")
# findASmain_bak(asType="A3SS", gpeFile="/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.gpe", sourceFile="/data/CrazyHsu_data/Projects/isoseq/wangLab_data_test/wangLab_data_test_20200326/Zheng58/8/ASE/deSingleExonRead.bed12+", outFile="a3ss.test")


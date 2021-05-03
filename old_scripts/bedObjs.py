#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: commonObjs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-07-27 14:45:39
Last modified: 2019-07-27 14:45:42
'''
from geneFeatures import *

class Bed6(object):
    def __init__(self, line=None):
        if line:
            self.record = line.strip("\n").split("\t")
            self.chrom, self.start, self.end, self.name, self.score, self.strand = self.record[:6]
            self.start = int(self.start)
            self.end = int(self.end)
            self.score = int(self.score) if self.score != '.' else self.score
        else:
            self.empty()

    def empty(self):
        (self.chrom, self.start, self.end, self.name, self.score, self.strand) = ('', 0, 0, '', 0, '')

    def __str__(self):
        return "%s:%d-%d:%s,%s" % (self.chrom, self.start, self.end, self.strand, str(self.score))

    __repr__ = __str__

class Bed12(object):
    def __init__(self, line=None):
        if line:
            self.record = line.strip("\n").split("\t")
            (self.chrom, self.start, self.end, self.name,
             self.score, self.strand, self.thickStart, self.thickEnd,
             self.rgb, self.blockCount, self.blockSizes,
             self.blockStarts) = self.record[:12]
            self.start = int(self.start)

            self.end = int(self.end)
            self.score = int(self.score)
            self.thickStart = int(self.thickStart)
            self.thickEnd = int(self.thickEnd)
            self.blockCount = int(self.blockCount)

            self.blockSizes = map(int, self.blockSizes.strip(",").split(","))
            self.blockStarts = map(int, self.blockStarts.strip(",").split(","))

            self.exonStarts = map(lambda x: x+self.start, self.blockStarts)
            self.exonEnds = map(lambda x: self.exonStarts[x]+self.blockSizes[x], range(self.blockCount))
            self.exonList = map(lambda x: [self.exonStarts[x], self.exonEnds[x]], range(self.blockCount))
            self.minpos = min(self.exonStarts)
            self.maxpos = max(self.exonEnds)
        else:
            self.empty()

    def empty(self):
        self.chrom = ""
        self.start, self.end = 0, 0
        self.name = ""
        self.score = 0
        self.strand = ""
        self.thickStart = self.thickEnd = 0
        self.rgb = "0,0,0"
        self.blockCount = 0
        self.blockSizes, self.blockStarts = [], []
        self.exonStarts, self.exonEnds = [], []
        self.exonList = []
        self.minpos, self.maxpos = 0, 0

    def __len__(self):
        return sum(self.blockSizes)

    def __repr__(self):
        fields = [self.chrom, str(self.start), str(self.end),
                  self.name, str(self.score), self.strand, str(self.thickStart),
                  str(self.thickEnd), self.rgb, str(self.blockCount)]
        blockSizesLine = ','.join(repr(i) for i in self.blockSizes)
        blockStartsLine = ','.join(repr(i) for i in self.blockStarts)
        fields.extend([blockSizesLine, blockStartsLine])
        return "\t".join(fields)

class BedLoader(object):
    def __init__(self, bedFile, sourceType):
        self.sourceType = sourceType
        self.genes = []
        self.geneNames = []
        self.loadBedFile(bedFile)

    def loadBedFile(self, bedFile):
        with open(bedFile) as f:
            for line in f.readlines():
                lineObj = RefBedLine(line) if self.sourceType == "ref" else AltBedLine(line)
                fakeTransId = lineObj.geneName + "_" + lineObj.asType
                if lineObj.geneName not in self.geneNames:
                    gene = GeneInfo(lineObj)
                    gene.exonDict = {lineObj.asType: ExonInfo(lineObj.exonList, transId=fakeTransId, strand=lineObj.strand, splicePos=lineObj.splicePos)}
                    # gene.exonDict = {lineObj.asType: ExonInfo(lineObj.exonList)}
                    gene.minposList = [lineObj.minpos]
                    gene.maxposList = [lineObj.maxpos]
                    self.geneNames.append(lineObj.geneName)
                    self.genes.append(gene)
                else:
                    gene = self.genes[-1]
                    gene.exonDict.update({lineObj.asType: ExonInfo(lineObj.exonList, transId=fakeTransId, strand=lineObj.strand, splicePos=lineObj.splicePos)})
                    # gene.exonDict.update({lineObj.asType: ExonInfo(lineObj.exonList)})
                    gene.minposList.append(lineObj.minpos)
                    gene.maxposList.append(lineObj.maxpos)

class AltBedLine(Bed12):
    def __init__(self, line):
        Bed12.__init__(self, line)
        record = line.strip("\n").split("\t")
        self.geneName = record[12]
        self.transList = record[14].split(",")
        self.asType = record[13]
        self.splicePos = record[15]
        self.juncPos = record[16]

class RefBedLine(Bed12):
    def __init__(self, line):
        Bed12.__init__(self, line)
        record = line.strip("\n").split("\t")
        self.geneName = record[12]
        self.transList = record[13].split(",")
        self.asType = record[14]
        self.splicePos = None
        self.juncPos = None


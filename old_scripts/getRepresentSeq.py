#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: getRepresentSeq.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-02-27 11:55:49
Last modified: 2021-02-27 11:55:49
'''

import pybedtools, subprocess

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
        return readsDict


readEnumerate = "reads.enumerate.txt"
isoBed = "tofu.collapsed.bed12+"
genome = "Zea_mays.AGPv4.36.dna.genome.fa"
isoGroupedFile = "tofu.collapsed.group.txt"

isoGroupDict = {}
with open(isoGroupedFile) as f:
    for i in f.readlines():
        infoList = i.strip("\n").split("\t")
        isoGroupDict[infoList[0]] = infoList[1].split(",")

allIsos = BedFile(isoBed, type="bed12+")

allBedOut = open("representSeqIso.all.bed12+", "w")
embryoBedOut = open("representSeqIso.embryo.bed12+", "w")
endoBedOut = open("representSeqIso.endo.bed12+", "w")
rootBedOut = open("representSeqIso.root.bed12+", "w")
embryoRootBedOut = open("representSeqIso.embryo-root.bed12+", "w")
embryoEndoBedOut = open("representSeqIso.embryo-endo.bed12+", "w")
endoRootBedOut = open("representSeqIso.endo-root.bed12+", "w")

newBedInfoList = []
with open(readEnumerate) as f:
    for i in f.readlines():
        infoList = i.strip("\n").split("\t")
        isos = [allIsos.reads[x] for x in infoList[3].split(",")]
        sortedIsos = sorted(isos, key=lambda x: len(isoGroupDict[x]), reverse=True)
        # sortedIsos = sorted(sorted(isos, key=lambda x: x.chromStart), key=lambda x: x.chromEnd, reverse=True)
        repIso = sortedIsos[0]
        if infoList[-1] == "all":
            print >> allBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "embryo":
            print >> embryoBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "endo":
            print >> endoBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "root":
            print >> rootBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "embryo-root":
            print >> embryoRootBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "embryo-endo":
            print >> embryoEndoBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
        if infoList[-1] == "endo-root":
            print >> endoRootBedOut, "\t".join(repIso.record[0:12]) + "\t" + infoList[3] + ":" + infoList[-1]
allBedOut.close()
embryoBedOut.close()
endoBedOut.close()
rootBedOut.close()
embryoRootBedOut.close()
embryoEndoBedOut.close()
endoRootBedOut.close()

cmd = "bed2gpe.pl -g 13 representSeqIso.all.bed12+ > representSeqIso.all.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.embryo.bed12+ > representSeqIso.embryo.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.endo.bed12+ > representSeqIso.endo.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.root.bed12+ > representSeqIso.root.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.embryo-root.bed12+ > representSeqIso.embryo-root.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.embryo-endo.bed12+ > representSeqIso.embryo-endo.gpe"
cmd = cmd + "\n" + "bed2gpe.pl -g 13 representSeqIso.endo-root.bed12+ > representSeqIso.endo-root.gpe"
subprocess.call(cmd, shell=True)

cmd = "genePredToGtf file representSeqIso.all.gpe representSeqIso.all.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.embryo.gpe representSeqIso.embryo.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.endo.gpe representSeqIso.endo.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.root.gpe representSeqIso.root.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.embryo-root.gpe representSeqIso.embryo-root.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.embryo-endo.gpe representSeqIso.embryo-endo.gtf"
cmd = cmd + "\n" + "genePredToGtf file representSeqIso.endo-root.gpe representSeqIso.endo-root.gtf"
subprocess.call(cmd, shell=True)

cmd = "gffread -w representSeqIso.all.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.all.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.embryo.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.embryo.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.endo.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.endo.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.root.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.root.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.embryo-root.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.embryo-root.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.embryo-endo.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.embryo-endo.gtf"
cmd = cmd + "\n" + "gffread -w representSeqIso.endo-root.exon.fa -g Zea_mays.AGPv4.36.dna.genome.fa representSeqIso.endo-root.gtf"
subprocess.call(cmd, shell=True)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: asEnumerate4_for_reference.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-02-27 22:12:58
Last modified: 2021-02-27 22:12:58
'''
import itertools

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

class BedFile(object):
    def __init__(self, bedFile, type=None):
        self.bedFile = bedFile
        self.reads = self.getReadsInfo(type)

    def getReadsInfo(self, type):
        readsDict = {}
        with open(self.bedFile) as f:
            for line in f:
                if type == "bed12+":
                    b = Bed12Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict

def isOverlap(listA, listB):
    sortedList = sorted([listA, listB])
    return True if sortedList[0][1] >= sortedList[1][0] else False

def asEnumerate1(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile):
    # isoform2reads = {}
    # with open(isoform2readsFile) as f:
    #     for line in f.readlines():
    #         isoform, reads = line.strip("\n").split("\t")
    #         isoform2reads[isoform] = reads.split(",")
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    gene2isoDict = {}
    for isoName in isoformBed:
        geneName = isoformBed[isoName].otherList[0]
        if geneName not in gene2isoDict:
            gene2isoDict[geneName] = [isoformBed[isoName]]
        else:
            gene2isoDict[geneName].append(isoformBed[isoName])
    juncCombDict = {}
    with open(irFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
            retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
            geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((retentionStart, retentionEnd))

    with open(seFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
            juncStart = int(name.split(":")[-1].split("@")[0])
            juncEnd = int(name.split(":")[-1].split("@")[-1])
            geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((juncStart, juncEnd))

    with open(a3ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
            excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
            excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
            geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))

    with open(a5ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
            excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
            excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
            geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))

    paDict = {}
    with open(paFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, strand, geneName, paSites = infoList[0:4]
            paSites = [int(x) for x in paSites.split(",")]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            if geneName in gene2isoDict:
                if strand == "+":
                    isoEnds = [x.chromEnd for x in gene2isoDict[geneName]]
                    pa = (min(isoEnds), max(isoEnds))
                else:
                    isoStarts = [x.chromStart for x in gene2isoDict[geneName]]
                    pa = (min(isoStarts), max(isoStarts))
            else:
                if strand == "+":
                    pa = (paSites[0] - 1 - 20, paSites[-1] + 20)
                else:
                    pa = (paSites[0] - 20, paSites[-1] + 1 + 20)
            juncCombDict[uniqueName].append(pa)
            paDict[geneName] = [pa, paSites]

    for uniqueName in juncCombDict:
        chrom, strand, geneName = uniqueName.split("@")
        if geneName not in gene2isoDict: continue
        isosInGene = gene2isoDict[geneName]
        juncs = sorted(juncCombDict[uniqueName])
        tmpJuncStart, tmpJuncEnd = juncs[0]
        newJuncs = []
        for junc in juncs[1:]:
            if isOverlap((tmpJuncStart, tmpJuncEnd), junc):
                tmpJuncStart = tmpJuncStart if tmpJuncStart < junc[0] else junc[0]
                tmpJuncEnd = tmpJuncEnd if tmpJuncEnd > junc[1] else junc[1]
            else:
                newJuncs.append((tmpJuncStart, tmpJuncEnd))
                tmpJuncStart = junc[0]
                tmpJuncEnd = junc[1]
        newJuncs.append((tmpJuncStart, tmpJuncEnd))
        newJuncs = sorted(newJuncs)
        # print newJuncs
        # print [x.name for x in isosInGene]
        eventCombination = {}
        isoformAssign = {}
        for junc in newJuncs:
            eventName = "event_{}-{}".format(junc[0], junc[1])
            eventCombination[eventName] = {}
            for iso in isosInGene:
                if not isOverlap(junc, (iso.chromStart, iso.chromEnd)):
                    continue
                juncsInIso = iso.introns
                if junc in juncsInIso:
                    choiceName = "choice_{}".format("-".join(map(str, junc)))
                    if choiceName not in eventCombination[eventName]:
                        eventCombination[eventName][choiceName] = [iso]
                    else:
                        eventCombination[eventName][choiceName].append(iso)
                else:
                    overlapedJuncs = []
                    for tmpJunc in juncsInIso:
                        if isOverlap(tmpJunc, junc):
                            overlapedJuncs.append(tmpJunc)
                    if overlapedJuncs:
                        choiceName = "choice_{}".format(";".join(["-".join(map(str, x)) for x in overlapedJuncs]))
                    else:
                        if geneName in paDict and isOverlap(junc, paDict[geneName][0]):
                            if strand == "+":
                                dist2pa = [abs(x - iso.chromEnd) for x in paDict[geneName][1]]
                            else:
                                dist2pa = [abs(x - iso.chromStart) for x in paDict[geneName][1]]
                            index = dist2pa.index(min(dist2pa))
                            choiceName = "choice_pa{}".format(index)
                        else:
                            choiceName = "choice_ir"
                    if choiceName not in eventCombination[eventName]:
                        eventCombination[eventName][choiceName] = [iso]
                    else:
                        eventCombination[eventName][choiceName].append(iso)
                if iso.name not in isoformAssign:
                    isoformAssign[iso.name] = {eventName: choiceName}
                else:
                    isoformAssign[iso.name].update({eventName: choiceName})
        # allEventChoice = []
        combCount = 1
        for eventName in eventCombination:
            allChoices = eventCombination[eventName]
            tmpList = []
            for choiceName in allChoices:
                tmpList.append((eventName, choiceName))
            # allEventChoice.append(tmpList)
            combCount *= len(tmpList)
        # combinations = list(itertools.product(*allEventChoice))
        print "\t".join(map(str, [chrom, strand, geneName, len(eventCombination), ";".join(eventCombination.keys()), "simCount", combCount]))
        print "\t".join(map(str, [chrom, strand, geneName, len(eventCombination), ";".join(eventCombination.keys()), "realCount", len(isosInGene)]))


irFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/characterization/IR.reference.bed6+"
seFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/characterization/SE.reference.bed6+"
a3ssFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/characterization/A3SS.reference.bed6+"
a5ssFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/characterization/A5SS.reference.bed6+"
paFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/characterization/paGroup.reference.tsv"
isoformFile = "/data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.36/Zea_mays.AGPv4.36.bed"
# isoform2readsFile = "../filtration/collapse/tofu.collapsed.group.txt"
asEnumerate1(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile)

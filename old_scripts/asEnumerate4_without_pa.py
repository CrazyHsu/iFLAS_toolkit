#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: asEnumerate4_without_pa.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-02-27 21:59:14
Last modified: 2021-02-27 21:59:14
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

def asEnumerate1(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
    isoform2reads = {}
    with open(isoform2readsFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
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
        # print isoformAssign
        # for eventName in eventCombination:
        #     for choiceName in eventCombination[eventName]:
        #         print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]]
        resultDict = {}
        allEventChoice = []
        choiceFreqList = []
        allReads = []
        for eventName in eventCombination:
            if eventName not in resultDict:
                resultDict[eventName] = {}
            allChoices = eventCombination[eventName]
            allIsosInEvent = list(set(itertools.chain.from_iterable(allChoices.values())))
            allReadsInEvent = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allIsosInEvent])))
            allReads.extend(allReadsInEvent)
            tmpList = []
            for choiceName in allChoices:
                tmpList.append((eventName, choiceName))
                readsInChoice = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allChoices[choiceName]])))
                choiceFreq = len(readsInChoice)/float(len(allReadsInEvent))
                resultDict[eventName].update({choiceName: {"readsInChoice": readsInChoice, "choiceFreq": choiceFreq}})
                # print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]], len(readsInChoice), len(allReadsInEvent)
            allEventChoice.append(tmpList)
            choiceFreqList.append(",".join(map(str, [resultDict[eventName][x]["choiceFreq"] for x in resultDict[eventName]])))
        combinations = list(itertools.product(*allEventChoice))
        # print combinations, len(combinations)
        # print choiceFreqList
        combination2iso = {}
        identifiedComb = []
        for comb in combinations:
            combination2iso[comb] = {"iso": [], "freq": "NA", "freqList": [], "reads": []}
            for iso in isoformAssign:
                if sorted(isoformAssign[iso].items()) == sorted(comb):
                    identifiedComb.append(comb)
                    combination2iso[comb]["iso"].append(iso)
                    freqList = []
                    for i in comb:
                        freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
                    combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
                    combination2iso[comb].update({"freqList": freqList})
                else:
                    freqList = []
                    for i in comb:
                        freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
                    combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
                    combination2iso[comb].update({"freqList": freqList})
            readsInComb = set(itertools.chain.from_iterable([isoform2reads[x] for x in combination2iso[comb]["iso"]]))
            counts = len(readsInComb)
            combination2iso[comb].update({"counts": counts})
            combination2iso[comb].update({"reads": readsInComb})
        # print combination2iso.values()
        identifiedComb = list(set(identifiedComb))
        combCounts = sorted([(combination2iso[x]["iso"], combination2iso[x]["counts"]) for x in identifiedComb], key=lambda x: x[1], reverse=True)
        majorCombCount = 0
        tmpSum = 0
        isoList = []
        for x in combCounts:
            majorCombCount += 1
            tmpSum += x
            isoList.append(x[0])
            if tmpSum/float(sum([y[1] for y in combCounts])) >= 0.75:
                break
        print "\t".join(map(str, [chrom, strand, geneName, len(eventCombination), ";".join(eventCombination.keys()), \
            ";".join(choiceFreqList), len(combinations), \
            ",".join(map(str, [combination2iso[x]["counts"] for x in identifiedComb])), \
            ",".join(map(str, [combination2iso[x]["freq"] for x in identifiedComb])), \
            len(list(set(identifiedComb))), len(list(set(allReads))), majorCombCount, ";".join([",".join(x) for x in isoList])]))



def asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
    isoform2reads = {}
    with open(isoform2readsFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    # juncsDict = getJuncsFromIsoformFile(isoformFile, isoform2reads)
    juncCombDict = {}
    with open(irFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
            retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            excJuncComb = "{}-{}".format(retentionStart, retentionEnd)
            incJuncComb = "ir_inc:{}-{}".format(retentionStart, retentionEnd)
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: excIsos})
                juncCombDict[uniqueName].update({incJuncComb: incIsos})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: excIsos})
                else:
                    juncCombDict[uniqueName][excJuncComb].extend(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: incIsos})
                else:
                    juncCombDict[uniqueName][incJuncComb].extend(incIsos)
    with open(seFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
            # blockSizes, relStarts = infoList[10], infoList[11]
            juncStart = int(name.split(":")[1].split("@")[0])
            juncEnd = int(name.split(":")[1].split("@")[-1])
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            incJuncComb = name.split(":")[1]
            excJuncComb = "{}-{}".format(juncStart, juncEnd)
            # skipedExonStart, skipedExonEnd, score = int(skipedExonStart), int(skipedExonEnd), float(score)
            # blockSizes = [int(x) for x in blockSizes.split(",")]
            # relStarts = [int(x) for x in relStarts.split(",")]
            # skipedExonStarts, skipedExonEnds = getAbsLoc(skipedExonStart, blockSizes, relStarts)

            incIsos = infoList[15].split(",")
            excIsos = infoList[17].split(",")
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: excIsos})
                juncCombDict[uniqueName].update({incJuncComb: incIsos})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: excIsos})
                else:
                    juncCombDict[uniqueName][excJuncComb].extend(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: incIsos})
                else:
                    juncCombDict[uniqueName][incJuncComb].extend(incIsos)
    with open(a3ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
            incJuncComb = infoList[-2].split(":")[-1]
            excJuncComb = infoList[-1].split(":")[-1]
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: excIsos})
                juncCombDict[uniqueName].update({incJuncComb: incIsos})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: excIsos})
                else:
                    juncCombDict[uniqueName][excJuncComb].extend(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: incIsos})
                else:
                    juncCombDict[uniqueName][incJuncComb].extend(incIsos)
    with open(a5ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
            incJuncComb = infoList[-2].split(":")[-1]
            excJuncComb = infoList[-1].split(":")[-1]
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: excIsos})
                juncCombDict[uniqueName].update({incJuncComb: incIsos})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: excIsos})
                else:
                    juncCombDict[uniqueName][excJuncComb].extend(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: incIsos})
                else:
                    juncCombDict[uniqueName][incJuncComb].extend(incIsos)
    paDict = {}
    with open(paFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, strand, geneName, paSites = infoList[0:4]
            freqs = infoList[6]
            paSites = [int(x) for x in paSites.split(",")]
            freqs = [float(x) for x in freqs.split(",")]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if strand == "+":
                pa = "pa:{}-{}".format(paSites[0] - 1, paSites[-1])
            else:
                pa = "pa:{}-{}".format(paSites[0], paSites[-1] + 1)
            paDict[uniqueName] = {pa: [freqs, paSites]}

    for key in juncCombDict:
        juncCombs = juncCombDict[key].keys()
        eventsInJuncsDict = {}
        for i in range(len(juncCombs)):
            if "ir_inc" in juncCombs[i] or "pa" in juncCombs[i]:
                juncLeftI = int(juncCombs[i].split(":")[1].split("-")[0])
                juncRightI = int(juncCombs[i].split(":")[1].split("-")[-1])
            elif "@" in juncCombs[i]:
                juncLeftI = int(juncCombs[i].split("@")[0])
                juncRightI = int(juncCombs[i].split("@")[-1])
            else:
                juncLeftI = int(juncCombs[i].split("-")[0])
                juncRightI = int(juncCombs[i].split("-")[1])
            tmpJuncLeft, tmpJuncRight = juncRightI, juncRightI
            clusteredJunc = []

    for key in juncCombDict:
        chrom, strand, geneName = key.split(":")
        juncCombs = juncCombDict[key].keys()
        eventsInJuncsDict = {}
        for i in range(len(juncCombs)):
            # irInc = []
            # if "ir_inc" in juncCombs[i]:
            #     irInc.append(juncCombs[i])
            if "ir_inc" in juncCombs[i]:
                juncLeftI = int(juncCombs[i].split(":")[1].split("-")[0])
                juncRightI = int(juncCombs[i].split(":")[1].split("-")[-1])
            elif "@" in juncCombs[i]:
                juncLeftI = int(juncCombs[i].split("@")[0])
                juncRightI = int(juncCombs[i].split("@")[-1])
            else:
                juncLeftI = int(juncCombs[i].split("-")[0])
                juncRightI = int(juncCombs[i].split("-")[1])
            tmpJuncLeft, tmpJuncRight = juncLeftI, juncRightI
            clusteredJunc = [juncCombs[i]]
            for j in range(i+1, len(juncCombs)):
                if "ir_inc" in juncCombs[j]:
                    juncLeftJ = int(juncCombs[j].split(":")[1].split("-")[0])
                    juncRightJ = int(juncCombs[j].split(":")[1].split("-")[-1])
                elif "@" in juncCombs[j]:
                    juncLeftJ = int(juncCombs[j].split("@")[0])
                    juncRightJ = int(juncCombs[j].split("@")[-1])
                else:
                    juncLeftJ = int(juncCombs[j].split("-")[0])
                    juncRightJ = int(juncCombs[j].split("-")[1])
                # if "ir_inc" in juncCombs[j]:
                #     irInc.append(juncCombs[j])
                if isOverlap((tmpJuncLeft, tmpJuncRight), (juncLeftJ, juncRightJ)):
                    tmpJuncLeft = tmpJuncLeft if tmpJuncLeft < juncLeftJ else juncLeftJ
                    tmpJuncRight = tmpJuncRight if tmpJuncRight > juncRightJ else juncRightJ
                    clusteredJunc.append(juncCombs[j])
            if bool(eventsInJuncsDict) == False:
                eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
                continue
            flag = 0
            for junc in eventsInJuncsDict.keys():
                if isOverlap(junc, (tmpJuncLeft, tmpJuncRight)):
                    if (tmpJuncLeft, tmpJuncRight) not in eventsInJuncsDict[junc]:
                        eventsInJuncsDict[junc].extend(clusteredJunc)
                    flag = 1
            if flag == 0:
                eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
            # if bool(eventsInJuncsDict) == False:
            #     eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
            # else:
            #     # eventsInJuncsDict.update({(tmpJuncLeft, tmpJuncRight): clusteredJunc})
            #     flag = 0
            #     for junc in eventsInJuncsDict.keys():
            #         if isOverlap(junc, (tmpJuncLeft, tmpJuncRight)):
            #             if (tmpJuncLeft, tmpJuncRight) not in eventsInJuncsDict[junc]:
            #                 eventsInJuncsDict[junc].append(juncCombs[j])
            #             flag = 1
            #     if flag == 0:
            #         eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
        eventCount = 0
        combinationDict = {}
        isoformAssign = {}
        print
        print eventsInJuncsDict
        for junc in eventsInJuncsDict:
            eventsInJuncsDict[junc] = list(set(eventsInJuncsDict[junc]))
        print
        print eventsInJuncsDict
        for juncComb in eventsInJuncsDict:
            # irTmpJuncLeft, irTmpJuncRight = 0, 0
            irDict = {}
            # irClusteredList = []
            irCount = 0
            eventName = "event_{}".format(eventCount)
            eventCount += 1
            choiceCount = 0
            if eventName not in combinationDict:
                combinationDict[eventName] = {}
            for tmpJunc in eventsInJuncsDict[juncComb]:
                if "ir_inc" in tmpJunc:
                    tmpJuncLeft = int(tmpJunc.split(":")[1].split("-")[0])
                    tmpJuncRight = int(tmpJunc.split(":")[1].split("-")[1])
                    if irCount not in irDict:
                        irDict[irCount] = {"events": [tmpJunc], "junc": [tmpJuncLeft, tmpJuncRight]}
                    else:
                        if isOverlap((tmpJuncLeft, tmpJuncRight), irDict[irCount]["junc"]):
                            irTmpJuncLeft = tmpJuncLeft if tmpJuncLeft < irDict[irCount]["junc"][0] else irDict[irCount]["junc"][0]
                            irTmpJuncRight = tmpJuncRight if tmpJuncRight > irDict[irCount]["junc"][1] else irDict[irCount]["junc"][1]
                            irDict[irCount]["events"].append(tmpJunc)
                            irDict[irCount]["junc"] = [irTmpJuncLeft, irTmpJuncRight]
                        else:
                            irCount += 1
                            irDict[irCount] = {"events": [tmpJunc], "junc": [tmpJuncLeft, tmpJuncRight]}

                    # if irTmpJuncLeft == 0 and irTmpJuncRight == 0:
                    #     irTmpJuncLeft = tmpJunc.split(":")[1].split("-")[0]
                    #     irTmpJuncRight = tmpJunc.split(":")[1].split("-")[1]
                    # else:
                    #     irTmpJuncLeft = irTmpJuncLeft if irTmpJuncLeft < tmpJuncLeft else tmpJuncLeft
                    #     irTmpJuncRight = irTmpJuncRight if irTmpJuncRight > tmpJuncRight else tmpJuncRight
                    # irClusteredList.append(tmpJunc)
                else:
                    choiceName = "choice_{}".format(choiceCount)
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})

                    if "@" in tmpJunc:
                        leftJunc, rightJunc = int(tmpJunc.split("@")[0]), int(tmpJunc.split("@")[-1])
                    else:
                        leftJunc, rightJunc = int(tmpJunc.split("-")[0]), int(tmpJunc.split("-")[1])
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in juncCombDict[key][tmpJunc]])))
                    # tmpIsos.extend(juncCombDict[key][tmpJunc])
                    tmpIsos = list(set(juncCombDict[key][tmpJunc]))
                    combinationDict[eventName][choiceName].update({"isoforms": tmpIsos})
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    combinationDict[eventName][choiceName].update({"junc": tmpJunc})
                    for iso in tmpIsos:
                        if iso not in isoformAssign:
                            isoformAssign[iso] = {"event": [eventName], "choice": [choiceName]}
                        else:
                            isoformAssign[iso]["event"].append(eventName)
                            isoformAssign[iso]["choice"].append(choiceName)
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    choiceCount += 1

            for ir in irDict:
                leftJunc, rightJunc = irDict[ir]["junc"]
                choiceName = "choice_{}_ir".format(choiceCount)
                tmpIsos = list(set(itertools.chain.from_iterable([juncCombDict[key][x] for x in irDict[ir]["events"]])))
                tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in tmpIsos])))
                combinationDict[eventName][choiceName] = {"reads": tmpReads, "junc": "{}-{}".format(leftJunc, rightJunc), "isoforms": tmpIsos}
                for iso in tmpIsos:
                    if iso not in isoformAssign:
                        isoformAssign[iso] = {"event": [eventName], "choice": [choiceName]}
                    else:
                        isoformAssign[iso]["event"].append(eventName)
                        isoformAssign[iso]["choice"].append(choiceName)
                tmpComb = "{}_{}".format(eventName, choiceName)
                choiceCount += 1

        paName = "event_pa"
        paCount = 0
        for paSite in paDict[key][-1]:
            paChoice = "choice_pa{}".format(paCount)
            tmpIsos = []
            for iso in isoformAssign:
                if isoformBed[iso].chromEnd == paSite:
                    isoformAssign[iso]["event"].append(paName)
                    isoformAssign[iso]["choice"].append(paChoice)
                    tmpIsos.append(iso)
                    continue
            paCount += 1
            tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in tmpIsos])))
            if paName not in combinationDict:
                combinationDict[paName] = {paChoice: {"reads": tmpReads, "junc": "{}-{}".format(paSite - 1, paSite), "isoforms": tmpIsos}}
            else:
                combinationDict[paName].update({paChoice: {"reads": tmpReads, "junc": "{}-{}".format(paSite-1, paSite), "isoforms": tmpIsos}})

        for eventName in combinationDict:
            for choiceName in combinationDict[eventName]:
                print eventName, choiceName, combinationDict[eventName][choiceName]["junc"]

# irFile = "Zm00001d050245.IR.bed6+"
# seFile = "Zm00001d050245.SE.PB.bed12+"
# a3ssFile = "Zm00001d050245.A3SS.PB.bed6+"
# a5ssFile = "Zm00001d050245.A5SS.PB.bed6+"
# paFile = "Zm00001d050245.PA.tsv"
# isoformFile = "Zm00001d050245.deSingleExonIsoform.bed12+"

irFile = "../ASE/PB/IR.bed6+"
seFile = "../ASE/PB/SE.bed12+"
a3ssFile = "../ASE/PB/A3SS.bed6+"
a5ssFile = "../ASE/PB/A5SS.bed6+"
paFile = "test.pa"
isoformFile = "../ASE/deSingleExonIsoform.bed12+"
isoform2readsFile = "../filtration/collapse/tofu.collapsed.group.txt"
asEnumerate1(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile)

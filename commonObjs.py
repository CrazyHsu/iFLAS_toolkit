import copy, re
from sys import maxint as MAXINT

AStypes = ["IR", "SE", "A3SS", "A5SS"]

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

    # def getCombinationEvents(self):
    #     asType2asKey = {"A5SS": "exc", "A3SS": "exc", "SE": "skip", "IR": "spliced"}
    #     combineEvent = "combineEvent"
    #     read2AS = {}
    #     self.asDict[combineEvent] = {}
    #     for asType in self.asDict:
    #         for event in self.asDict[asType]:
    #             if asType2asKey[asType] in self.asDict[asType][event]:
    #                 for read in self.asDict[asType][event][asType2asKey[asType]]:
    #                     if read not in read2AS:
    #                         read2AS[read] = {asType: [event]}
    #                     elif asType not in read2AS[read]:
    #                         read2AS[read].update({asType: [event]})
    #                     else:
    #                         read2AS[read][asType].append(event)
    #     for r in read2AS:
    #         if len(read2AS[r]) == 1: continue
    #         for t, l in read2AS[r].items():
    #             for e in l:
    #                 self.asDict[t][e][asType2asKey[t]].remove(r)
    #                 if r not in self.asDict[combineEvent]:
    #                     self.asDict[combineEvent][r] = {t: [e]}
    #                 elif t not in self.asDict[combineEvent][r]:
    #                     self.asDict[combineEvent][r].update({t: [e]})
    #                 else:
    #                     self.asDict[combineEvent][r][t].append(e)
    #                 # self.asDict[combineEvent].add(r)
    #                 if len(self.asDict[t][e][asType2asKey[t]]) == 0:
    #                     del self.asDict[t][e]
    #
    # def construcFakeRef(self):
    #     asType2asKey = {"A5SS": "exc", "A3SS": "exc", "SE": "skip", "IR": "spliced"}
    #     for asType in asType2asKey:
    #         for event in self.asDict[asType]:
    #             reads = self.asDict[asType][event][asType2asKey[asType]]
    #             maxLen = reads[0].exonLen
    #             maxLenRead = reads[0]
    #             for i in range(1, len(reads)):
    #                 if reads[i].exonLen > maxLen:
    #                     maxLen = reads[i].exonLen
    #                     maxLenRead = reads[i]
    #             records = copy.deepcopy(maxLenRead.record)
    #             if asType == "SE":
    #                 excExonList = map(lambda x: (int(x.split("-")[0]), int(x.split("-")[1])), event.split("@")[1].split(";"))
    #                 excExonStarts = map(lambda x: x[0], excExonList)
    #                 excExonEnds = map(lambda x: x[1], excExonList)
    #                 newExonStarts = sorted(maxLenRead.blockStarts + excExonStarts)
    #                 newExonEnds = sorted(maxLenRead.blockEnds + excExonEnds)
    #                 # newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
    #                 # newExonRelStarts = map(lambda x: newExonEnds[x] - maxLenRead.start, range(len(newExonStarts)))
    #                 # records[3] = records[3]+"_"+asType+"_fake"
    #                 # records[9] = len(newExonSizes)
    #                 # records[10] = newExonRelStarts
    #                 # recordLine = "\t".join(records)
    #                 # fakeRead = ReadLineStruc(recordLine)
    #             else:
    #                 # records = maxLenRead.records
    #                 excStart, excEnd = map(int, re.split("-", event))
    #                 mergedExonList = mergeTupleList(maxLenRead.readExonList + [(excStart, excEnd)])
    #                 newExonStarts = map(lambda x: x[0], mergedExonList)
    #                 newExonEnds = map(lambda x: x[1], mergedExonList)
    #             newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
    #             newExonRelStarts = map(lambda x: newExonStarts[x] - maxLenRead.start, range(len(newExonStarts)))
    #             records[3] = records[3] + "_" + asType + "_fake"
    #             records[9] = len(newExonSizes)
    #             records[10] = ",".join(map(str, newExonSizes))
    #             records[11] = ",".join(map(str, newExonRelStarts))
    #             recordLine = "\t".join(map(str, records))
    #             fakeRead = ReadLineStruc(recordLine)
    #             self.asDict[asType][event][asType2asKey[asType]].append(fakeRead)
    #     fakeReadDict = {}
    #     for r in self.asDict["combineEvent"]:
    #         newExonStarts = r.blockStarts
    #         newExonEnds = r.blockStarts
    #         records = copy.deepcopy(r.record)
    #         newExonList = []
    #         for t in self.asDict["combineEvent"][r]:
    #             for e in self.asDict["combineEvent"][r][t]:
    #                 if t == "SE":
    #                     excExonList = map(lambda x: (int(x.split("-")[0]), int(x.split("-")[1])),
    #                                       e.split("@")[1].split(";"))
    #                     newExonList.extend(excExonList)
    #                     # excExonStarts = map(lambda x: int(x.split("-")[0]), e.split("@")[1].split(";"))
    #                     # excExonEnds = map(lambda x: int(x.split("-")[1]), e.split("@")[1].split(";"))
    #                     # newExonStarts = sorted(newExonStarts + excExonStarts)
    #                     # newExonEnds = sorted(newExonEnds + excExonEnds)
    #                     # newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
    #                     # newExonRelStarts = map(lambda x: newExonEnds[x] - maxLenRead.start, range(len(newExonStarts)))
    #                 else:
    #                     excStart, excEnd = map(int, re.split("-", e))
    #                     newExonList.append((excStart, excEnd))
    #         mergedExonList = mergeTupleList(r.readExonList + newExonList)
    #         newExonStarts = map(lambda x: x[0], mergedExonList)
    #         newExonEnds = map(lambda x: x[1], mergedExonList)
    #         newExonSizes = map(lambda x: newExonEnds[x] - newExonStarts[x], range(len(newExonStarts)))
    #         newExonRelStarts = map(lambda x: newExonStarts[x] - r.start, range(len(newExonStarts)))
    #         records[3] = records[3] + "_combineEvent_fake11"
    #         records[9] = len(newExonSizes)
    #         records[10] = ",".join(map(str, newExonSizes))
    #         records[11] = ",".join(map(str, newExonRelStarts))
    #         recordLine = "\t".join(map(str, records))
    #         fakeRead = ReadLineStruc(recordLine)
    #         fakeReadDict[fakeRead] = {}
    #     self.asDict["combineEvent"].update(fakeReadDict)

    def __str__(self):
        return "%s:%s-%s (%s) (%s) (%d reads)" % \
               (self.chrom, self.minpos, self.maxpos, self.geneName, self.strand, len(self.reads))

    __repr__ = __str__


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

    def changeGeneId(self, geneId):
        self.geneName2gpeObj = {}
        self.genePredDict = self.buildObj(geneId=geneId)

    def gpeObj2file(self, fout, geneName):
        f = open(fout, "w")
        gps = self.geneName2gpeObj[geneName]
        for i in gps:
            print >>f, i
        f.close()

    def getBlockLength(self, blockList):
        return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))

    def getGeneExonLength(self):
        exonLength = 0
        for g in self.geneName2gpeObj:
            exonList = []
            for l in self.geneName2gpeObj[g]:
                exonList.extend(l.exons)
            exonList = set(exonList)
            exonLength += self.getBlockLength(exonList)
        return exonLength

    def getGeneIntronLength(self):
        intronLength = 0
        for g in self.geneName2gpeObj:
            intronList = []
            for l in self.geneName2gpeObj[g]:
                intronList.extend(l.introns)
            intronList = set(intronList)
            intronLength += self.getBlockLength(intronList)
        return intronLength

    # def redefineBorder(self, specificMinpos, specificMaxpos, outFile="tmp.gpe"):
    #     out = open(outFile, "w")
    #     for i in self.geneName2gpeObj:
    #         for j in self.geneName2gpeObj[i]:
    #             records = copy.deepcopy(j.record)
    #             if j.txStart < specificMinpos:
    #                 for z in range(len(j.exonEnds)):
    #                     if j.exonStarts[z] >= specificMinpos and j.exonEnds[z-1] <= specificMaxpos:
    #                         records[3] = records[5] = str(j.exonStarts[z])
    #                         records[7] = str(j.exonCount - z)
    #                         records[8] = ",".join(map(str, j.exonStarts[z:]))
    #                         records[9] = ",".join(map(str, j.exonEnds[z:]))
    #                     elif j.exonStarts[z] <= specificMinpos and specificMinpos < j.exonEnds[z]:
    #                         records[3] = records[5] = str(specificMinpos)
    #                         records[7] = str(j.exonCount - z)
    #                         records[8] = ",".join(map(str, [specificMinpos] + j.exonStarts[z+1:]))
    #                         records[9] = ",".join(map(str, j.exonEnds[z:]))
    #             if j.txEnd > specificMaxpos:
    #                 for z in range(len(j.exonEnds)):
    #                     if j.exonEnds[z-1] <= specificMaxpos and j.exonStarts[z] >= specificMaxpos:
    #                         records[4] = records[6] = str(j.exonEnds[z])
    #                         records[7] = str(j.exonCount - len(j.exonStarts) + z + 1)
    #                         records[8] = ",".join(map(str, j.exonStarts[:z]))
    #                         records[9] = ",".join(map(str, j.exonEnds[:z]))
    #                     elif j.exonStarts[z] <= specificMaxpos and specificMaxpos < j.exonEnds[z]:
    #                         records[4] = records[6] = str(specificMaxpos)
    #                         records[7] = str(j.exonCount - len(j.exonStarts) + z + 1)
    #                         records[8] = ",".join(map(str, j.exonStarts[:z+1]))
    #                         records[9] = ",".join(map(str, j.exonEnds[:z] + [specificMaxpos]))
    #             recordStr = "\t".join(records)
    #             gpeLine = GenePredExtLine(recordStr, bincolumn=False)
    #             print >>out, str(gpeLine)
    #     out.close()
    #     return outFile

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

class MergedTgsSample(object):
    # def __init__(self):
    #     self.projectName = None
    #     self.sampleName = None
    #     self.refStrain = None
    #
    #     self.tgsPlat = None
    #     self.tgsStrategy = None
    #     self.tgsProcessedData = None
    #     self.tgsPrimer = None
    #     self.dataLocation = None
    #
    #     self.ngsLeftReads = None
    #     self.ngsRightReads = None
    #     self.ngsReadPair = "paired"
    #     self.ngsReadsLength = 150
    #     self.ngsJunctions = None
    #     self.ngsCondition = None
    #
    #     self.runThreads = None
    #     self.runMemory = None

    def __init__(self):
        self.project_name = None
        self.sample_name = None
        self.ref_strain = None
        self.strain = None
        # self.condition = None

        self.tgs_plat = None
        self.strategy = None
        self.data_location = None
        self.data_processed_location = None
        # self.primer = None
        self.polya_location = None

        self.ngs_left_reads = None
        self.ngs_right_reads = None
        self.ngs_reads_paired = "paired"
        self.ngs_reads_length = None
        self.ngs_junctions = None

        self.use_fmlrc2 = True
        self.single_run_threads = 1

    # def getMergedName(self):
    #     self.mergedName = "{}_all".format(self.strain)

    def __str__(self):
        return "%s:%s_all:%s, %s:%s:%s, %s:%s:%s:%s:%s, %d" % \
               (self.project_name, self.sample_name, self.ref_strain, self.tgs_plat, self.strategy,
                self.data_processed_location, self.ngs_left_reads, self.ngs_right_reads, self.ngs_reads_paired,
                self.ngs_reads_length, self.ngs_junctions, self.single_run_threads)

    __repr__ = __str__

class ReadLineStruc(Bed12):
    """
    get Bed12+ read structure
    """
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.record = line.strip().split("\t")
        self.exonLen = sum(self.blockSizes)
        self.readLen = abs(self.chromStart - self.chromEnd)
        self.geneName = self.record[12]

    def getExonFramesForGPE(self, gpeObj):
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

    def go_to_gpe(self):
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
        gpe.score = 0
        gpe.geneName = self.geneName
        gpe.cdsStartStat = "unk"
        gpe.cdsEndStat = "unk"
        gpe.exonFrames = self.getExonFramesForGPE(gpe)
        return gpe

    def __str__(self):
        return "%s:%s:%d-%d (%s)" % (self.name, self.chrom, self.chromStart, self.chromEnd, self.geneName)

    __repr__ = __str__


#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: alleleSpecificReadsAssign.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-03-15 20:46:29
Last modified: 2021-03-15 20:46:29
'''
import sys, os, itertools, glob, copy
import pandas as pd

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
            self.exonChain = self.chrom + ":" + ';'.join(map(lambda x: str(self.exons[x][0])+"-"+str(self.exons[x][1]), range(len(self.exons))))
            self.juncChain = self.chrom + ":" + ";".join(map(lambda x: str(self.introns[x][0])+"-"+str(self.introns[x][1]), range(len(self.introns))))
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

def getIso2reads(myFile):
    iso2reads = {}
    with open(myFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            if infoList[0] not in iso2reads:
                iso2reads[infoList[0]] = infoList[1].split(",")
    return iso2reads

def parseHaplotype(cleanHapInfo, cleanHapHumanReadable, vcfHapDict, iso2reads):
    partialDF = pd.read_csv(cleanHapHumanReadable, skiprows=[0], index_col=0, sep="\t")
    index = list(partialDF.index)
    columns = list(partialDF.columns)
    returnDict = {}
    hap2reads = {}
    with open(cleanHapInfo) as f:
        for i in f.readlines()[1:]:
            infoList = i.strip("\n").split(",")
            if infoList[2] not in hap2reads:
                hap2reads[infoList[2]] = []
            hap2reads[infoList[2]].append(infoList[0])
    for iso in columns:
        for i in index:
            parent = vcfHapDict[i]
            if iso not in returnDict:
                returnDict[iso] = {parent: list(set([x+"/ccs" for x in hap2reads[i]]) & set(iso2reads[iso]))}
            else:
                returnDict[iso].update({parent: list(set([x+"/ccs" for x in hap2reads[i]]) & set(iso2reads[iso]))})

    return returnDict

def parseVCF(cleanVCF):
    refStr, altStr = "", ""
    with open(cleanVCF) as f:
        for i in f.readlines():
            if i.startswith("#"): continue
            infoList = i.strip("\n").split("\t")
            refStr += infoList[3]
            altStr += infoList[4]
    return {refStr: "B73", altStr: "Ki11"}

def getAlleleSpecificInfo(targetDir, iso2reads):
    cleanVCF = os.path.join(targetDir, "phased.partial.cleaned.vcf")
    cleanHapInfo = os.path.join(targetDir, "phased.partial.cleaned.hap_info.txt")
    cleanHapHumanReadable = os.path.join(targetDir, "phased.partial.cleaned.human_readable.txt")
    vcfHapDict = parseVCF(cleanVCF)
    try:
        cleanHapDict = parseHaplotype(cleanHapInfo, cleanHapHumanReadable, vcfHapDict, iso2reads)
        return cleanHapDict
    except:
        return []

def readsAssign(myDir, myName=None, isHybrid=False):
    isoBedFile = os.path.join(myDir, "filtration/collapse/tofu.collapsed.assigned.unambi.bed12+")
    iso2readsFile = os.path.join(myDir, "filtration/collapse/tofu.collapsed.group.txt")
    iso2reads = getIso2reads(iso2readsFile)
    isoBedObj = BedFile(isoBedFile, type="bed12+")
    gene2iso = {}
    for iso in isoBedObj.reads:
        if isoBedObj.reads[iso].otherList[0] not in gene2iso:
            gene2iso[isoBedObj.reads[iso].otherList[0]] = []
        gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])
    if not isHybrid:
        mergedIsoDict = {myName: {}}
        for gene in gene2iso:
            tmpDict = {}
            refGene = ""
            for iso in gene2iso[gene]:
                refGene = iso.otherList[-1]
                if iso.juncChain not in tmpDict:
                    tmpDict[iso.juncChain] = [iso]
                else:
                    tmpDict[iso.juncChain].append(iso)

            for tmp in tmpDict:
                isos = tmpDict[tmp]
                mergedIsoName = "+".join([x.name for x in isos])
                reads = list(itertools.chain.from_iterable([iso2reads[x.name] for x in isos]))
                sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                repIso = copy.copy(sortedIsos[0])
                repIso.name = mergedIsoName
                juncChain = repIso.juncChain
                mergedIsoDict[myName][juncChain] = [isos, repIso, reads, refGene, "parent", {myName: reads}]
        return mergedIsoDict
    else:
        mergedIsoDict = {myName: {}}
        for gene in gene2iso:
            tmpDict = {}
            refGene = ""
            for iso in gene2iso[gene]:
                refGene = iso.otherList[-1]
                if iso.juncChain not in tmpDict:
                    tmpDict[iso.juncChain] = [iso]
                else:
                    tmpDict[iso.juncChain].append(iso)

            lociDir = glob.glob(os.path.join(myDir, "alleleSpecificRelatedAS/by_loci/{}_size*".format(gene)))
            flag = True
            if len(lociDir) == 0:
                flag = False
                cleanHap = ""
            else:
                cleanHap = os.path.join(lociDir[0], "phased.partial.cleaned.human_readable.txt")
            if os.path.exists(cleanHap) and flag:
                cleanHapDict = getAlleleSpecificInfo(lociDir[0], iso2reads)
                if len(cleanHapDict) == 0:
                    flag = False
                if flag:
                    for tmp in tmpDict:
                        isos = tmpDict[tmp]
                        reads = list(itertools.chain.from_iterable([iso2reads[x.name] for x in isos]))
                        # mergedIsoName = "+".join([x.name for x in isos])
                        # sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                        # repIso = copy.copy(sortedIsos[0])
                        # repIso.name = mergedIsoName
                        # juncChain = repIso.juncChain
                        # B73_reads = list(itertools.chain.from_iterable([cleanHapDict[x.name]["B73"] for x in isos if x.name in cleanHapDict and "B73" in cleanHapDict[x.name]]))
                        # Ki11_reads = list(itertools.chain.from_iterable([cleanHapDict[x.name]["Ki11"] for x in isos if x.name in cleanHapDict and "Ki11" in cleanHapDict[x.name]]))
                        # mergedIsoDict[myName][juncChain] = [isos, repIso, reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": Ki11_reads}]

                        B73_isos_names = []
                        Ki11_isos_names = []
                        filteredIsos = [x for x in isos if x.name in cleanHapDict]
                        for i in cleanHapDict:
                            if "B73" in cleanHapDict[i]:
                                B73_isos_names.append(i)
                            if "Ki11" in cleanHapDict[i]:
                                Ki11_isos_names.append(i)
                        B73_reads = list(itertools.chain.from_iterable([cleanHapDict[x.name]["B73"] for x in filteredIsos if x.name in B73_isos_names]))
                        Ki11_reads = list(itertools.chain.from_iterable([cleanHapDict[x.name]["Ki11"] for x in filteredIsos if x.name in Ki11_isos_names]))

                        B73_isos = sorted([i for i in filteredIsos if i.name in B73_isos_names], key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                        Ki11_isos = sorted([i for i in filteredIsos if i.name in Ki11_isos_names], key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                        if len(B73_isos) > 0 and len(Ki11_isos) > 0:
                            B73_repIso = copy.copy(B73_isos[0])
                            B73_repIso.name = "+".join(B73_isos_names)
                            B73JuncChain = B73_repIso.juncChain
                            Ki11_repIso = copy.copy(Ki11_isos[0])
                            Ki11_repIso.name = "+".join(Ki11_isos_names)
                            Ki11JuncChain = Ki11_repIso.juncChain
                            if B73JuncChain == Ki11JuncChain:
                                tmpIso = sorted([B73_repIso, Ki11_repIso], key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                                mergedIsoDict[myName][B73JuncChain] = [filteredIsos, tmpIso[0], tmpIso[0], reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": Ki11_reads}]
                            else:
                                mergedIsoDict[myName][B73JuncChain] = [filteredIsos, B73_repIso, "NA", reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": "NA"}]
                                mergedIsoDict[myName][Ki11JuncChain] = [filteredIsos, "NA", Ki11_repIso, reads, refGene, "hybrid", {"B73": "NA", "Ki11": Ki11_reads}]
                        elif len(B73_isos) == 0 and len(Ki11_isos) > 0:
                            Ki11_repIso = copy.copy(Ki11_isos[0])
                            Ki11_repIso.name = "+".join(Ki11_isos_names)
                            Ki11JuncChain = Ki11_repIso.juncChain
                            mergedIsoDict[myName][Ki11JuncChain] = [filteredIsos, "NA", Ki11_repIso, reads, refGene, "hybrid", {"B73": "NA", "Ki11": Ki11_reads}]
                        elif len(B73_isos) > 0 and len(Ki11_isos) == 0:
                            B73_repIso = copy.copy(B73_isos[0])
                            B73_repIso.name = "+".join(B73_isos_names)
                            B73JuncChain = B73_repIso.juncChain
                            mergedIsoDict[myName][B73JuncChain] = [filteredIsos, B73_repIso, "NA", reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": "NA"}]
                        # if len(Ki11_isos):
                        #     Ki11_repIso = copy.copy(Ki11_isos[0])
                        #     Ki11_repIso.name = "+".join(Ki11_isos_names)
                        #     Ki11JuncChain = Ki11_repIso.juncChain
                        # B73_repIso = copy.copy(B73_isos[0])
                        # Ki11_repIso = copy.copy(Ki11_isos[0])
                        # B73_repIso.name = "+".join(B73_isos_names)
                        # Ki11_repIso.name = "+".join(Ki11_isos_names)
                        # B73JuncChain = B73_repIso.juncChain
                        # Ki11JuncChain = Ki11_repIso.juncChain
                        # if B73JuncChain == Ki11JuncChain:
                        #     tmpIso = sorted([B73_repIso, Ki11_repIso], key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                        #     mergedIsoDict[myName][B73JuncChain] = [filteredIsos, tmpIso[0], tmpIso[0], reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": Ki11_reads}]
                        # else:
                        #     mergedIsoDict[myName][B73JuncChain] = [filteredIsos, B73_repIso, "NA", reads, refGene, "hybrid", {"B73": B73_reads, "Ki11": "NA"}]
                        #     mergedIsoDict[myName][Ki11JuncChain] = [filteredIsos, "NA", Ki11_repIso, reads, refGene, "hybrid", {"B73": "NA", "Ki11": Ki11_reads}]
            else:
                flag = False
                    # print [x.name for x in isos]
                    # print cleanHapDict
                    # mergedIsoDict[myName][juncChain] = [isos, repIso, reads, refGene, "hybrid", {"B73": cleanHapDict[mergedIsoName]["B73"], "Ki11": cleanHapDict[mergedIsoName]["Ki11"]}]
            if not flag:
                for tmp in tmpDict:
                    isos = tmpDict[tmp]
                    mergedIsoName = "+".join([x.name for x in isos])
                    reads = list(itertools.chain.from_iterable([iso2reads[x.name] for x in isos]))
                    sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                    repIso = copy.copy(sortedIsos[0])
                    repIso.name = mergedIsoName
                    juncChain = repIso.juncChain
                    mergedIsoDict[myName][juncChain] = [isos, repIso, "NA", reads, refGene, "hybrid", {"B73": reads, "Ki11": "NA"}]
        return mergedIsoDict

def getColValue(myFile, column, sep = "\t"):
    return list(pd.read_csv(myFile, sep=sep, header=None).iloc[:, column])

def readsTissueAssign(allReads, strain = "B73", sample2reads=None):
    tissues = ["embryo", "endo", "root"]
    readsTissueAssignDict = {}
    for i in tissues:
        sample = strain + "_" + i
        intersectList = list(sample2reads[sample] & set(allReads))
        readsTissueAssignDict[i] = len(intersectList)
    return readsTissueAssignDict

def getSampleReads():
    tissues = ["embryo", "endo", "root"]
    strains = ["B73", "Ki11", "B73_Ki11", "Ki11_B73"]
    baseDir = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test_v4.50/test1/"
    sample2reads = {}
    for i in tissues:
        for j in strains:
            sample = j + "_" + i
            flncBedFile = os.path.join(baseDir, sample, "rawCorrection/flnc.mm2.sorted.bed12")
            mySet = set(getColValue(flncBedFile, 3, sep="\t"))
            sample2reads.update({sample: mySet})
    return sample2reads

def print_test(myDict, myName=""):
    for i in myDict[myName]:
        print i, myDict[myName][i]

def tgsAssign():
    B73_dir = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test_v4.50/test1/B73_all"
    Ki11_dir = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test_v4.50/test1/Ki11_all"
    B73_Ki11_dir = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test_v4.50/test1/B73_Ki11_all"
    Ki11_B73_dir = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test_v4.50/test1/Ki11_B73_all"
    B73_dict = readsAssign(B73_dir, myName="B73", isHybrid=False)
    Ki11_dict = readsAssign(Ki11_dir, myName="Ki11", isHybrid=False)
    B73_Ki11_dict = readsAssign(B73_Ki11_dir, myName="B73_Ki11", isHybrid=True)
    Ki11_B73_dict = readsAssign(Ki11_B73_dir, myName="Ki11_B73", isHybrid=True)
    sample2reads = getSampleReads()
    print_test(B73_dict, myName="B73")
    print_test(Ki11_dict, myName="Ki11")
    print_test(B73_Ki11_dict, myName="B73_Ki11")
    print_test(Ki11_B73_dict, myName="Ki11_B73")
    print "before set"
    allJuncChain = list(set(B73_dict["B73"]) | set(Ki11_dict["Ki11"]) | set(B73_Ki11_dict["B73_Ki11"]) | set(Ki11_B73_dict["Ki11_B73"]))

    no_tissue_specific_out = open("no_tissue_specific.iso_count.txt", "w")
    tissue_specific_out = open("tissue_specific.iso_count.txt", "w")
    merged_iso_out = open("all_sample_merged_iso.bed", "w")
    name_mapping_out = open("name_mapping.txt", "w")
    print "enter for loop"
    count = 0
    merged_iso_out_B73_Ki11_B73 = open("all_sample_merged_iso.B73_Ki11_B73.bed", "w")
    merged_iso_out_B73_Ki11_Ki11 = open("all_sample_merged_iso.B73_Ki11_Ki11.bed", "w")
    merged_iso_out_Ki11_B73_B73 = open("all_sample_merged_iso.Ki11_B73_B73.bed", "w")
    merged_iso_out_Ki11_B73_Ki11 = open("all_sample_merged_iso.Ki11_B73_Ki11.bed", "w")

    for juncChain in allJuncChain:
        if juncChain in B73_dict["B73"]:
            B73_reads, B73_refGene, B73_repIso = B73_dict["B73"][juncChain][5]["B73"], B73_dict["B73"][juncChain][3], B73_dict["B73"][juncChain][1]
        else:
            B73_reads, B73_refGene, B73_repIso = "NA", "NA", "NA"
        if juncChain in Ki11_dict["Ki11"]:
            Ki11_reads, Ki11_refGene, Ki11_repIso = Ki11_dict["Ki11"][juncChain][5]["Ki11"], Ki11_dict["Ki11"][juncChain][3], Ki11_dict["Ki11"][juncChain][1]
        else:
            Ki11_reads, Ki11_refGene, Ki11_repIso = "NA", "NA", "NA"



        if juncChain in B73_Ki11_dict["B73_Ki11"]:
            B73_Ki11_reads_B73, B73_Ki11_repIso_B73 = B73_Ki11_dict["B73_Ki11"][juncChain][6]["B73"], B73_Ki11_dict["B73_Ki11"][juncChain][1]
            if B73_refGene == "NA" and Ki11_refGene == "NA":
                B73_refGene = B73_Ki11_dict["B73_Ki11"][juncChain][4]
        else:
            B73_Ki11_reads_B73, B73_Ki11_repIso_B73 = "NA", "NA"

        if juncChain in B73_Ki11_dict["B73_Ki11"]:
            B73_Ki11_reads_Ki11, B73_Ki11_repIso_Ki11 = B73_Ki11_dict["B73_Ki11"][juncChain][6]["Ki11"], B73_Ki11_dict["B73_Ki11"][juncChain][2]
            if B73_refGene == "NA" and Ki11_refGene == "NA":
                B73_refGene = B73_Ki11_dict["B73_Ki11"][juncChain][4]
        else:
            B73_Ki11_reads_Ki11, B73_Ki11_repIso_Ki11 = "NA", "NA"



        if juncChain in Ki11_B73_dict["Ki11_B73"]:
            Ki11_B73_reads_B73, Ki11_B73_repIso_B73 = Ki11_B73_dict["Ki11_B73"][juncChain][6]["B73"], Ki11_B73_dict["Ki11_B73"][juncChain][1]
            if B73_refGene == "NA" and Ki11_refGene == "NA":
                Ki11_refGene = Ki11_B73_dict["Ki11_B73"][juncChain][4]
        else:
            Ki11_B73_reads_B73, Ki11_B73_repIso_B73 = "NA", "NA"

        if juncChain in Ki11_B73_dict["Ki11_B73"]:
            Ki11_B73_reads_Ki11, Ki11_B73_repIso_Ki11 = Ki11_B73_dict["Ki11_B73"][juncChain][6]["Ki11"], Ki11_B73_dict["Ki11_B73"][juncChain][2]
            if B73_refGene == "NA" and Ki11_refGene == "NA":
                Ki11_refGene = Ki11_B73_dict["Ki11_B73"][juncChain][4]
        else:
            Ki11_B73_reads_Ki11, Ki11_B73_repIso_Ki11 = "NA", "NA"

        B73_count = len(B73_reads) if B73_reads != "NA" else 0
        Ki11_count = len(Ki11_reads) if Ki11_reads != "NA" else 0
        B73_Ki11_count_B73 = len(B73_Ki11_reads_B73) if B73_Ki11_reads_B73 != "NA" else 0
        B73_Ki11_count_Ki11 = len(B73_Ki11_reads_Ki11) if B73_Ki11_reads_Ki11 != "NA" else 0
        Ki11_B73_count_B73 = len(Ki11_B73_reads_B73) if Ki11_B73_reads_B73 != "NA" else 0
        Ki11_B73_count_Ki11 = len(Ki11_B73_reads_Ki11) if Ki11_B73_reads_Ki11 != "NA" else 0

        if B73_reads != "NA":
            B73_readsTissueAssigned = readsTissueAssign(B73_reads, strain="B73", sample2reads=sample2reads)
        else:
            B73_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}
        if Ki11_reads != "NA":
            Ki11_readsTissueAssigned = readsTissueAssign(Ki11_reads, strain="Ki11", sample2reads=sample2reads)
        else:
            Ki11_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}
        if B73_Ki11_reads_B73 != "NA":
            B73_Ki11_reads_B73_readsTissueAssigned = readsTissueAssign(B73_Ki11_reads_B73, strain="B73_Ki11", sample2reads=sample2reads)
        else:
            B73_Ki11_reads_B73_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}
        if B73_Ki11_reads_Ki11 != "NA":
            B73_Ki11_reads_Ki11_readsTissueAssigned = readsTissueAssign(B73_Ki11_reads_Ki11, strain="B73_Ki11", sample2reads=sample2reads)
        else:
            B73_Ki11_reads_Ki11_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}
        if Ki11_B73_reads_B73 != "NA":
            Ki11_B73_reads_B73_readsTissueAssigned = readsTissueAssign(Ki11_B73_reads_B73, strain="Ki11_B73", sample2reads=sample2reads)
        else:
            Ki11_B73_reads_B73_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}
        if Ki11_B73_reads_Ki11 != "NA":
            Ki11_B73_reads_Ki11_readsTissueAssigned = readsTissueAssign(Ki11_B73_reads_Ki11, strain="Ki11_B73", sample2reads=sample2reads)
        else:
            Ki11_B73_reads_Ki11_readsTissueAssigned = {"embryo": 0, "endo": 0, "root": 0}

        tissues = ["embryo", "endo", "root"]
        B73_readsTissueAssigned_count = [B73_readsTissueAssigned[i] for i in tissues]
        Ki11_readsTissueAssigned_count = [Ki11_readsTissueAssigned[i] for i in tissues]
        B73_Ki11_reads_B73_readsTissueAssigned_count = [B73_Ki11_reads_B73_readsTissueAssigned[i] for i in tissues]
        B73_Ki11_reads_Ki11_readsTissueAssigned_count = [B73_Ki11_reads_Ki11_readsTissueAssigned[i] for i in tissues]
        Ki11_B73_reads_B73_readsTissueAssigned_count = [Ki11_B73_reads_B73_readsTissueAssigned[i] for i in tissues]
        Ki11_B73_reads_Ki11_readsTissueAssigned_count = [Ki11_B73_reads_Ki11_readsTissueAssigned[i] for i in tissues]

        all_sample_repIso_list = [B73_repIso, Ki11_repIso, B73_Ki11_repIso_B73, B73_Ki11_repIso_Ki11,
                                  Ki11_B73_repIso_B73, Ki11_B73_repIso_Ki11]
        all_sample_repIso_sorted = sorted([i for i in all_sample_repIso_list if i != "NA"], key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)

        B73_repIso_name = B73_repIso.name if B73_repIso != "NA" else "NA"
        Ki11_repIso_name = Ki11_repIso.name if Ki11_repIso != "NA" else "NA"
        B73_Ki11_repIso_B73_name = B73_Ki11_repIso_B73.name if B73_Ki11_repIso_B73 != "NA" else "NA"
        B73_Ki11_repIso_Ki11_name = B73_Ki11_repIso_Ki11.name if B73_Ki11_repIso_Ki11 != "NA" else "NA"
        Ki11_B73_repIso_B73_name = Ki11_B73_repIso_B73.name if Ki11_B73_repIso_B73 != "NA" else "NA"
        Ki11_B73_repIso_Ki11_name = Ki11_B73_repIso_Ki11.name if Ki11_B73_repIso_Ki11 != "NA" else "NA"
        merged_repIso_name = "B73_{}:Ki11_{}:B73_Ki11_B73_{}:B73_Ki11_Ki11_{}:Ki11_B73_B73_{}:Ki11_B73_Ki11_{}".format(
            B73_repIso_name, Ki11_repIso_name, B73_Ki11_repIso_B73_name, B73_Ki11_repIso_Ki11_name,
            Ki11_B73_repIso_B73_name, Ki11_B73_repIso_Ki11_name)
        all_sample_repIso = all_sample_repIso_sorted[0]
        count += 1
        simple_name = "iflas_merged_isoform.{}".format(count)
        all_sample_repIso.name = simple_name
        print >>name_mapping_out, "{}\t{}".format(simple_name, merged_repIso_name)
        if B73_refGene.startswith("Zm") or Ki11_refGene.startswith("Zm"):
            if B73_refGene.startswith("Zm") and not Ki11_refGene.startswith("Zm"):
                ref_gene = B73_refGene
            elif not B73_refGene.startswith("Zm") and Ki11_refGene.startswith("Zm"):
                ref_gene = Ki11_refGene
            elif B73_refGene == Ki11_refGene:
                ref_gene = B73_refGene
            else:
                ref_gene = "ambiguous1"
            print >>no_tissue_specific_out, "\t".join(map(str, [juncChain, ref_gene, simple_name, B73_count, Ki11_count, B73_Ki11_count_B73,
                                                                B73_Ki11_count_Ki11, Ki11_B73_count_B73, Ki11_B73_count_Ki11]))
            tissueCountsList = [juncChain, ref_gene, simple_name]
            tissueCountsList += B73_readsTissueAssigned_count
            tissueCountsList += Ki11_readsTissueAssigned_count
            tissueCountsList += B73_Ki11_reads_B73_readsTissueAssigned_count
            tissueCountsList += B73_Ki11_reads_Ki11_readsTissueAssigned_count
            tissueCountsList += Ki11_B73_reads_B73_readsTissueAssigned_count
            tissueCountsList += Ki11_B73_reads_Ki11_readsTissueAssigned_count
            print >>tissue_specific_out, "\t".join(map(str, tissueCountsList))
        else:
            print >>no_tissue_specific_out, "\t".join(map(str, [juncChain, "ambiguous2", simple_name, B73_count, Ki11_count, B73_Ki11_count_B73,
                                                                B73_Ki11_count_Ki11, Ki11_B73_count_B73, Ki11_B73_count_Ki11]))
            tissueCountsList = [juncChain, "ambiguous2", simple_name]
            tissueCountsList += B73_readsTissueAssigned_count
            tissueCountsList += Ki11_readsTissueAssigned_count
            tissueCountsList += B73_Ki11_reads_B73_readsTissueAssigned_count
            tissueCountsList += B73_Ki11_reads_Ki11_readsTissueAssigned_count
            tissueCountsList += Ki11_B73_reads_B73_readsTissueAssigned_count
            tissueCountsList += Ki11_B73_reads_Ki11_readsTissueAssigned_count
            print >>tissue_specific_out, "\t".join(map(str, tissueCountsList))
            all_sample_repIso.otherList[0] = "ambiguous2"

        if sum(B73_Ki11_reads_B73_readsTissueAssigned_count) > 0:
            print >> merged_iso_out_B73_Ki11_B73, str(all_sample_repIso)
        if sum(B73_Ki11_reads_Ki11_readsTissueAssigned_count) > 0:
            print >> merged_iso_out_B73_Ki11_Ki11, str(all_sample_repIso)
        if sum(Ki11_B73_reads_B73_readsTissueAssigned_count) > 0:
            print >> merged_iso_out_Ki11_B73_B73, str(all_sample_repIso)
        if sum(Ki11_B73_reads_Ki11_readsTissueAssigned_count) > 0:
            print >> merged_iso_out_Ki11_B73_Ki11, str(all_sample_repIso)
        print >> merged_iso_out, str(all_sample_repIso)

    no_tissue_specific_out.close()
    tissue_specific_out.close()
    merged_iso_out.close()
    name_mapping_out.close()
    merged_iso_out_B73_Ki11_B73.close()
    merged_iso_out_B73_Ki11_Ki11.close()
    merged_iso_out_Ki11_B73_B73.close()
    merged_iso_out_Ki11_B73_Ki11.close()
    print "function end"


def ngsAssign():
    pass

def main():
    tgsAssign()
    # ngsAssign()

if __name__ == '__main__':
    main()
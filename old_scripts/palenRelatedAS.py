#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: palenRelatedAS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-02 16:01:42
Last modified: 2020-05-02 16:01:42
'''
import subprocess, argparse, itertools
import scipy.stats as stats
from commonObjs import *
from Config import *
# from relationshipBetweenPAlenAndAS import relationshipBetweenPAlenAndAS, isoformWithDiffPalen

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
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
                if type == "bed12+":
                    b = Bed12Plus(line)
                    readsDict.__setitem__(b.name, b)
        return readsDict

def relationshipBetweenPAlenAndAS(asFile, collapsedGroupFile, flncReads2Palen, asType="SE", filterByCount=5, sigFile=None, sigPicOutDir=None):
    currentFileDir = os.path.dirname(os.path.realpath(__file__))
    # flncReads2Palen = getDictFromFile(flncFile, sep=",", valueCol=5)
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    # isoformDict = getDictFromFile(isoformBed, keyCol=4, sep="\t")

    asType = asType.upper()
    out = open(sigFile, "w")
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            elif asType == "PA":
                inclusionIsos = []
                exclusionIsos = []
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            for item in itertools.product(inclusionIsos, exclusionIsos):
                inclusionReads = collapsedTrans2reads[item[0]]
                exclusionReads = collapsedTrans2reads[item[1]]
                inclusionReads2palen = [int(flncReads2Palen[i]) for i in inclusionReads if i in flncReads2Palen]
                exclusionReads2palen = [int(flncReads2Palen[i]) for i in exclusionReads if i in flncReads2Palen]
                if filterByCount:
                    if len(inclusionReads2palen) < filterByCount or len(exclusionReads2palen) < filterByCount:
                        continue
                stat_val, p_val = stats.ttest_ind(inclusionReads2palen, exclusionReads2palen)
                if float(p_val) <= 0.001:
                    print >> out, "\t".join(records) + "\t" + "\t".join(
                        ["inclusion", ",".join(map(str, inclusionReads2palen)), str(len(inclusionReads2palen)), "_".join(item)])
                    print >> out, "\t".join(records) + "\t" + "\t".join(
                        ["exclusion", ",".join(map(str, exclusionReads2palen)), str(len(exclusionReads2palen)), "_".join(item)])
    out.close()

    if os.stat(sigFile).st_size != 0:
        plotGeneStructurePath = os.path.join(currentFileDir, "plotGeneStructure.R")
        cmd = "Rscript {} -sig={} -od={}".format(plotGeneStructurePath, sigFile, sigPicOutDir)
        # subprocess.call(cmd, shell=True)

def palenRelatedAS(tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Analysis the relationship between poly-A length and AS events for project {} entry {}...".format(
        tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    groupDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
    resolveDir(os.path.join(groupDir, "palenRelatedAS"))

    collapsedGroupFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    flncFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "pacbioDemultiplexed", tgsSample.uniqName+".flnc.report.csv")
    palenFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "nanoporeRetrieve", tgsSample.uniqName+".palen.tsv")
    isoformBed = os.path.join(groupDir, "ASE", "deSingleExonIsoform.AS.confident.bed12+")

    isoformWithDiffPalen(collapsedGroupFile, palenFile)

    seFile = os.path.join(groupDir, "ASE", "PB", "SE.confident.bed12+")
    irFile = os.path.join(groupDir, "ASE", "PB", "IR.confident.bed6+")
    a5ssFile = os.path.join(groupDir, "ASE", "PB", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(groupDir, "ASE", "PB", "A3SS.confident.bed6+")

    relationshipBetweenPAlenAndAS(seFile, collapsedGroupFile, flncFile, isoformBed, asType="SE", filterByCount=5,
                                  sigFile="SE.palenAndAS.sig.bed12+", sigPicOutDir="SE.sig_pics")
    relationshipBetweenPAlenAndAS(irFile, collapsedGroupFile, flncFile, isoformBed, asType="IR", filterByCount=5,
                                  sigFile="IR.palenAndAS.sig.bed12+", sigPicOutDir="IR.sig_pics")
    relationshipBetweenPAlenAndAS(a5ssFile, collapsedGroupFile, flncFile, isoformBed, asType="A5SS", filterByCount=5,
                                  sigFile="A5SS.palenAndAS.sig.bed12+", sigPicOutDir="A5SS.sig_pics")
    relationshipBetweenPAlenAndAS(a3ssFile, collapsedGroupFile, flncFile, isoformBed, asType="A3SS", filterByCount=5,
                                  sigFile="A3SS.palenAndAS.sig.bed12+", sigPicOutDir="A3SS.sig_pics")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Analysis the relationship between poly-A length and AS events for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)

def getPalenAS(flncReads2Palen, tgsSample=None, refParams=None, dirSpec=None, collapsedTrans2reads=None, mergedIsoDict=None, filterByCount=10, merged=False):
    gene2readPalen = {}
    prevDir = os.getcwd()
    if not merged:
        resolveDir("no_merge")
    else:
        resolveDir("merged")

    if not merged:
        for isoform in collapsedTrans2reads:
            geneId = ".".join(isoform.split(".")[0:2])
            reads = collapsedTrans2reads[isoform]
            readsWithPolyaLen = [i for i in reads if i in flncReads2Palen]
            reads2palen = dict(zip(readsWithPolyaLen, [int(flncReads2Palen[i]) for i in readsWithPolyaLen]))
            if geneId not in gene2readPalen:
                gene2readPalen[geneId] = {isoform: reads2palen}
            else:
                gene2readPalen[geneId].update({isoform: reads2palen})
    else:
        for isoform in mergedIsoDict:
            geneId = mergedIsoDict[isoform][2]
            reads = mergedIsoDict[isoform][0]
            readsWithPolyaLen = [i for i in reads if i in flncReads2Palen]
            reads2palen = dict(zip(readsWithPolyaLen, [int(flncReads2Palen[i]) for i in readsWithPolyaLen]))
            if geneId not in gene2readPalen:
                gene2readPalen[geneId] = {isoform: reads2palen}
            else:
                gene2readPalen[geneId].update({isoform: reads2palen})

    count = 0
    for g in gene2readPalen:
        # allReadsInGene = sum([len(gene2readPalen[g][i]) for i in gene2readPalen[g]])
        for a, b in itertools.combinations(gene2readPalen[g].keys(), 2):
            if len(gene2readPalen[g][a].keys()) >= filterByCount and len(
                    gene2readPalen[g][b].keys()) >= filterByCount:
                aPalen = map(int, gene2readPalen[g][a].values())
                bPalen = map(int, gene2readPalen[g][b].values())
                stat_val1, p_val1 = stats.ttest_ind(aPalen, bPalen)
                stat_val2, p_val2 = stats.kruskal(aPalen, bPalen)
                if float(p_val1) <= 0.001 and float(p_val2) <= 0.001:
                    # print "\t".join(map(str, [allReadsInGene, g, a, len(aPalen), np.mean(aPalen), np.median(aPalen), b, len(bPalen), np.mean(bPalen), np.median(bPalen)]))
                    if merged:
                        count += 1
                        fileOut = "{}_{}.txt".format(g, count)
                    else:
                        fileOut = "_".join([g, a, b]) + ".txt"
                    out = open(fileOut, "w")
                    for i in gene2readPalen[g][a]:
                        print >> out, "\t".join(map(str, [g, a, i, gene2readPalen[g][a][i]]))
                    for j in gene2readPalen[g][b]:
                        print >> out, "\t".join(map(str, [g, b, j, gene2readPalen[g][b][j]]))
                    out.close()
                    # plotGeneStructurePath = os.path.join(currentFileDir, "violin_paired.R")
                    # cmd = "Rscript {} {}".format(plotGeneStructurePath, fileOut)
                    cmd = "Rscript violin_paired.R {}".format(fileOut)
                    # subprocess.call(cmd, shell=True)

    if merged:
        cmd = "getSingleIr.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.IR.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getSE.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.SE.bed12+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getA5SS.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.A5SS.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getA3SS.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.A3SS.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        irFile = os.path.join(os.getcwd(), "mergedIso.IR.bed6+")
        seFile = os.path.join(os.getcwd(), "mergedIso.SE.bed12+")
        a5ssFile = os.path.join(os.getcwd(), "mergedIso.A5SS.bed6+")
        a3ssFile = os.path.join(os.getcwd(), "mergedIso.A3SS.bed6+")
        collapsedGroupFile = os.path.join(os.getcwd(), "../mergedIso2Reads.group.txt")
    else:
        groupDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
        seFile = os.path.join(groupDir, "ASE", "PB", "SE.confident.bed12+")
        irFile = os.path.join(groupDir, "ASE", "PB", "IR.confident.bed6+")
        a5ssFile = os.path.join(groupDir, "ASE", "PB", "A5SS.confident.bed6+")
        a3ssFile = os.path.join(groupDir, "ASE", "PB", "A3SS.confident.bed6+")
        collapsedGroupFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    relationshipBetweenPAlenAndAS(irFile, collapsedGroupFile, flncReads2Palen, asType="IR", filterByCount=5,
                                  sigFile="IR.palenAndAS.sig.bed12+", sigPicOutDir="IR.sig_pics")
    relationshipBetweenPAlenAndAS(seFile, collapsedGroupFile, flncReads2Palen, asType="SE", filterByCount=5,
                                  sigFile="SE.palenAndAS.sig.bed12+", sigPicOutDir="SE.sig_pics")
    relationshipBetweenPAlenAndAS(a5ssFile, collapsedGroupFile, flncReads2Palen, asType="A5SS", filterByCount=5,
                                  sigFile="A5SS.palenAndAS.sig.bed12+", sigPicOutDir="A5SS.sig_pics")
    relationshipBetweenPAlenAndAS(a3ssFile, collapsedGroupFile, flncReads2Palen, asType="A3SS", filterByCount=5,
                                  sigFile="A3SS.palenAndAS.sig.bed12+", sigPicOutDir="A3SS.sig_pics")
    os.chdir(prevDir)

def mergeIso(isoBedFile, collapsedTrans2reads):
    isoBedObj = BedFile(isoBedFile, type="bed12+")
    gene2iso = {}
    for iso in isoBedObj.reads:
        if isoBedObj.reads[iso].otherList[0] not in gene2iso:
            gene2iso[isoBedObj.reads[iso].otherList[0]] = []
        gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])
    mergedIsoDict = {}
    mergedIso2ReadsGroupOut = open("mergedIso2Reads.group.txt", "w")
    mergedIso2ReadsBed = open("mergedIso2Reads.bed12", "w")
    for gene in gene2iso:
        tmpDict = {}
        for iso in gene2iso[gene]:
            if iso.juncChain not in tmpDict:
                tmpDict[iso.juncChain] = [iso]
            else:
                tmpDict[iso.juncChain].append(iso)

        for tmp in tmpDict:
            isos = tmpDict[tmp]
            mergedIsoName = "+".join([x.name for x in isos])
            reads = list(itertools.chain.from_iterable([collapsedTrans2reads[x.name] for x in isos]))
            sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
            repIso = sortedIsos[0]
            repIso.name = mergedIsoName
            mergedIsoDict[mergedIsoName] = [reads, repIso, gene]
            print >> mergedIso2ReadsGroupOut, "{}\t{}".format(mergedIsoName, ",".join(reads))
            print >> mergedIso2ReadsBed, str(repIso)
    mergedIso2ReadsBed.close()
    mergedIso2ReadsGroupOut.close()
    return mergedIsoDict

def palenRelatedAS1(tgsSample=None, refParams=None, dirSpec=None, filterByCount=10):
    print str(datetime.datetime.now()) + " Analysis the relationship between poly-A length and AS events for project {} entry {}...".format(
        tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    groupDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
    resolveDir(os.path.join(groupDir, "palenRelatedAS"))

    collapsedGroupFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    isoBedFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    mergedIsoDict = mergeIso(isoBedFile, collapsedTrans2reads)

    palenFile = os.path.join(dirSpec.out_dir, tgsSample.projectName, "retrieveData", "nanoporeRetrieve", tgsSample.sampleName, "polya_results.tsv")
    palen = pd.read_csv(palenFile, sep="\t")
    palen_pass = palen.loc[palen.qc_tag == "PASS", ]
    flncReads2Palen = dict(zip(palen_pass.readname, palen_pass.polya_length))

    getPalenAS(flncReads2Palen, tgsSample=tgsSample, refParams=refParams, dirSpec=dirSpec, collapsedTrans2reads=collapsedTrans2reads, filterByCount=filterByCount, merged=False)
    getPalenAS(flncReads2Palen, tgsSample=tgsSample, refParams=refParams, dirSpec=dirSpec, mergedIsoDict=mergedIsoDict, filterByCount=filterByCount, merged=True)

    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Analysis the relationship between poly-A length and AS events for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refInfoParams = defaultCfg.refInfoParams
    mecatParams = defaultCfg.mecatParams
    ccsParams = defaultCfg.ccsParams
    optionTools = defaultCfg.optionTools
    lncRNAcpatParams = defaultCfg.lncRNAcpatParams
    dirSpec = defaultCfg.dirParams

    projectCfg1 = ProjectCfg1(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    tgsDict = projectCfg1.tgsDict
    # tgsList = projectCfg1.tgsList
    # tgsStrain2sample = projectCfg1.tgsStrain2sample
    tgsPlat2sample = projectCfg1.tgsPlat2sample
    # tgsRefStrain2sample = projectCfg1.tgsRefStrain2sample
    tgsGroupedSample = projectCfg1.tgsGroupedSample

    ngsDict = projectCfg1.ngsDict

    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
        # initSysSetting(refParams=refParams, optionTools=optionTools, dirSpec=dirSpec)
    initResourceSetting(optionTools=optionTools)

    # retrieveTGSdata(tgsPlat2sample, dirSpec, refInfoParams, ccsParams, optionTools)

    sampleDict = {}
    myproject2refStrain = {}
    for projectName in tgsDict:
        sampleDict[projectName] = []
        myproject2refStrain[projectName] = {}
        for tgsPlat in tgsPlat2sample[projectName]:
            if tgsPlat == "pacbio":
                demultiplexDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "pacbioDemultiplexed")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        myproject2refStrain[projectName][refStrain] = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            for i in groupedSampleName.split("+"):
                                if i not in ngsDict[projectName][strainName]: continue
                                sampleName = strainName + "_" + i
                                ngsSample = ngsDict[projectName][strainName][i]

                                mySample = MergedTgsSample()
                                mySample.projectName = projectName
                                mySample.sampleName = sampleName
                                mySample.refStrain = refStrain
                                mySample.getUniqName()
                                mySample.tgsPlat = tgsPlat
                                mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                                tgsUniqName = "{}_{}_{}".format(projectName, strainName, groupedSampleName)
                                mySample.tgsProcessedData = os.path.join(demultiplexDir,
                                                                         "{}.{}.flnc.bam".format(tgsUniqName, i))
                                if "endo" in i:
                                    cmd = "samtools cat {}/{}.{}.flnc.bam {}/test1_All_endo_B73_endo+Ki11_endo+B73_Ki11_endo+Ki11_B73_endo.{}_endo.flnc.bam > {}/{}.{}.merged.flnc.bam".format(
                                        demultiplexDir, tgsUniqName, i, demultiplexDir, strainName, demultiplexDir,
                                        tgsUniqName, i)
                                    # subprocess.call(cmd, shell=True)
                                    mySample.tgsProcessedData = os.path.join(demultiplexDir,
                                                                             "{}.{}.merged.flnc.bam".format(tgsUniqName,
                                                                                                            i))
                                mySample.tgsPrimer = os.path.join(demultiplexDir,
                                                                  "{}.{}.primers.fa".format(tgsUniqName, i))
                                mySample.ngsLeftReads = ngsSample.leftReads
                                mySample.ngsRightReads = ngsSample.rightReads
                                mySample.ngsPaired = ngsSample.ngsReadPair
                                mySample.ngsReadsLength = ngsSample.readLength

                                sampleDict[projectName].append(mySample)
                                myproject2refStrain[projectName][refStrain].append(mySample)
            else:
                nanoporeDir = os.path.join(dirSpec.out_dir, projectName, "retrieveData", "nanoporeRetrieve")
                for refStrain in tgsDict[projectName]:
                    for strainName in tgsDict[projectName][refStrain]:
                        myproject2refStrain[projectName][refStrain] = []
                        for groupedSampleName in tgsGroupedSample[projectName]:
                            sampleName = strainName + "_" + groupedSampleName
                            ngsSample = ngsDict[projectName][strainName][groupedSampleName]

                            mySample = MergedTgsSample()
                            mySample.projectName = projectName
                            mySample.sampleName = sampleName
                            mySample.refStrain = refStrain
                            mySample.getUniqName()
                            mySample.tgsPlat = tgsPlat
                            mySample.dataLocation = tgsGroupedSample[projectName][groupedSampleName].dataLocation
                            mySample.tgsStrategy = tgsGroupedSample[projectName][groupedSampleName].tgsStrategy
                            mySample.tgsProcessedData = os.path.join(nanoporeDir, sampleName, "nanoporeRawSeq.fq")
                            mySample.ngsLeftReads = ngsSample.leftReads
                            mySample.ngsRightReads = ngsSample.rightReads
                            mySample.ngsPaired = ngsSample.ngsReadPair
                            mySample.ngsReadsLength = ngsSample.readLength

                            sampleDict[projectName].append(mySample)
                            myproject2refStrain[projectName][refStrain].append(mySample)

    for projectName in sampleDict:
        for tmpRun in sampleDict[projectName]:
            refParams = refInfoParams[tmpRun.refStrain]
            if tmpRun.ngsLeftReads or tmpRun.ngsRightReads:
                pass
                # makeHisat2Index(refParams=refParams, dirSpec=dirSpec, optionTools=optionTools)

    pybedtools.set_tempdir(dirSpec.tmp_dir)
    projectNum = sum([len(sampleDict[i]) for i in sampleDict])
    # pool = MyPool(processes=projectNum)
    for projectName in sampleDict:
        for tmpRun in sampleDict[projectName]:
            refParams = refInfoParams[tmpRun.refStrain]
            tmpRun.threads = int(optionTools.threads / float(projectNum))
            tmpRun.memory = str(int(optionTools.memory[0:-1]) / float(projectNum)) + "M"
            palenRelatedAS1(tmpRun, refParams, dirSpec, 10)
            # pool.apply_async(palenRelatedAS1, (tmpRun, refParams, dirSpec, 10))
    # pool.close()
    # pool.join()

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    # integratedRun(defaultCfg=defaultCfg)
    #    integratedRun1(defaultCfg=defaultCfg)
    main(defaultCfg=defaultCfg)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: getJuncCountFromNGSassembly.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-01-13 20:24:07
Last modified: 2020-01-13 20:24:08
'''

import subprocess, pysam, re, itertools, os, psutil, pybedtools
from multiprocessing import Pool
import numpy as np
from sklearn.cluster import KMeans
from scipy import stats

######################## classes ##############################
class TopHatJunction(object):
    def __init__(self):
        self.chrom = None
        self.start = None
        self.end = None
        self.name = None
        self.juncSupport = None
        self.strand = None
        self.juncLeft = None
        self.juncRight = None

    def getBlockSize(self):
        return (self.juncLeft - self.start, self.end - self.juncRight)

    def getRelStarts(self):
        return (0, self.juncRight - self.start)

    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.start, self.end, self.name, self.juncSupport, self.strand,
                                   self.start, self.end, "255,0,0", 2,
                                   ",".join(map(str, self.getBlockSize())), ",".join(map(str, self.getRelStarts()))]))

    __repr__ = __str__

######################## functions ##############################
def listInt2Str(myList):
    return [str(i) for i in myList]

def getCountsFromHisat2MappingLog(logFile):
    conOneTimePat = re.compile("(.*) (.*) aligned concordantly exactly 1 time")
    conMoreOneTimePat = re.compile("(.*) (.*) aligned concordantly >1 times")
    disconOneTimePat = re.compile("(.*) (.*) aligned discordantly 1 time")
    disconExactOneTimePat = re.compile("(.*) (.*) aligned exactly 1 time")
    disconMoreOneTimePat = re.compile("(.*) (.*) aligned >1 times")
    count = 0
    with open(logFile) as f:
        for line in f.readlines():
            if conOneTimePat.search(line):
                count += int(conOneTimePat.search(line).group(1).strip())
                continue
            if conMoreOneTimePat.search(line):
                count += int(conMoreOneTimePat.search(line).group(1).strip())
                continue
            if disconOneTimePat.search(line):
                count += int(disconOneTimePat.search(line).group(1).strip())
                continue
            if disconExactOneTimePat.search(line):
                count += int(disconExactOneTimePat.search(line).group(1).strip())
                continue
            if disconMoreOneTimePat.search(line):
                count += int(disconMoreOneTimePat.search(line).group(1).strip())
                continue
    return count

def parallelJuncCount2(bedFile, bamFile, outFile, allReadCounts=5000000, filterCount=None):
    out = open(outFile, "w")
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    allReadsStream = pybedtools.BedTool(bamFile).bam_to_bed(split=True)
    with open(bedFile) as f:
        count = 1
        readsBamList = []
        juncFlankBedList = []
        juncBedList = []
        juncNameList = []
        juncInfoDict = {}
        for line in f:
            lineList = line.strip("\n").split("\t")
            chrom, start, end, strand = lineList[0], int(lineList[1]), int(lineList[2]), lineList[3]
            reads = bam.fetch(chrom, start, end)
            # allReadsBam = "{}_{}_{}.bam".format(chrom, start, end)
            # if allReadsBam in readsBamList: continue
            # readsBamList.append(allReadsBam)
            # allReadsOut = pysam.AlignmentFile(allReadsBam, "wb", template=bam)
            readListWithJunc = []
            # print >> debug, chrom, start, end
            for read in reads:
                if read.cigarstring and "N" in read.cigarstring:
                    # print >> debug, read
                    if re.match("^\d+S", read.cigarstring):
                        intronList = [
                            (read, read.get_blocks()[(j - 2) / 2][1] - 1, read.get_blocks()[(j - 2) / 2 + 1][0]) for j
                            in range(len(read.cigartuples)) if read.cigartuples[j][0] == 3]
                    else:
                        intronList = [
                            (read, read.get_blocks()[(j - 1) / 2][1] - 1, read.get_blocks()[(j - 1) / 2 + 1][0]) for j
                            in range(len(read.cigartuples)) if read.cigartuples[j][0] == 3]
                    for intron in intronList:
                        if start == intron[1] and end == intron[2]:
                            readListWithJunc.append(intron[0])
            #     allReadsOut.write(read)
            # allReadsOut.close()
            if len(readListWithJunc) == 0: continue
            if filterCount and len(readListWithJunc) < filterCount: continue
            readsMinPos = min([i.reference_start for i in readListWithJunc])
            readsMaxPos = max([i.reference_end for i in readListWithJunc])
            junc = TopHatJunction()
            junc.chrom, junc.start, junc.end, junc.strand, junc.juncSupport = chrom, readsMinPos, readsMaxPos, strand, len(
                readListWithJunc)
            if "chr" in chrom.lower():
                junc.name = "JUNC_{}_".format(chrom.upper()) + "{0:08d}".format(count)
            else:
                junc.name = "JUNC_{}_".format(chrom) + "{0:08d}".format(count)
            junc.juncLeft, junc.juncRight = start+1, end

            juncLeftExonToIntron = "{}\t{}\t{}\t{}\t{}\t{}".format(chrom, junc.juncLeft - 3, junc.juncLeft + 2, junc.name, 0, strand)
            juncRightIntronToExon = "{}\t{}\t{}\t{}\t{}\t{}".format(chrom, junc.juncRight - 2, junc.juncRight + 3, junc.name, 0, strand)
            juncBed = "{}\t{}\t{}\t{}\t{}\t{}".format(chrom, junc.juncLeft, junc.juncRight, junc.name, 0, strand)
            juncFlankBedList.append(juncLeftExonToIntron)
            juncFlankBedList.append(juncRightIntronToExon)
            juncBedList.append(juncBed)
            juncNameList.append(junc.name)
            juncInfoDict.update({junc.name: junc})
            count += 1

        juncFlankBedStr = "\n".join(juncFlankBedList)
        juncFlankBedObj = pybedtools.BedTool(juncFlankBedStr, from_string=True)
        juncFlankCoverAllDf = juncFlankBedObj.coverage(allReadsStream).to_dataframe()
        juncFlankCoverF1Df = juncFlankBedObj.coverage(allReadsStream, f=1).to_dataframe()

        juncBedStr = "\n".join(juncBedList)
        juncBedObj = pybedtools.BedTool(juncBedStr, from_string=True)
        juncCoverDf = juncBedObj.coverage(allReadsStream, d=True).to_dataframe()
        juncCoverDfGrouped = juncCoverDf.groupby(juncCoverDf.name)

        juncFlankCoverF1Df = juncFlankCoverF1Df.loc[juncFlankCoverF1Df.blockCount == 1, ]
        juncFlankFilterF1Df = juncFlankCoverF1Df.groupby("name").size().astype(int).reset_index(name='Count')
        juncNameWithFlankF1 = juncFlankFilterF1Df.loc[juncFlankFilterF1Df.Count == 2,].name
        juncFlankCoverF1Df = juncFlankCoverF1Df.loc[juncFlankCoverF1Df.name.isin(juncNameWithFlankF1), ]
        juncFlankCoverF1DfGrouped = juncFlankCoverF1Df.groupby(juncFlankCoverF1Df.name)
        juncFlankCoverAllDfGrouped = juncFlankCoverAllDf.groupby(juncFlankCoverAllDf.name)

        for juncName in juncNameList:
            juncIntronWithReadsSupp, juncIntronEvenlyDist, juncIntronLowCoverCount = "no", "no", 0
            juncCv, juncFPKM, juncIntronRetion = 0, 0, "no"
            if juncName in juncNameWithFlankF1.values:
                leftJunc = juncFlankCoverF1DfGrouped.get_group(juncName).iloc[0,]
                rightJunc = juncFlankCoverF1DfGrouped.get_group(juncName).iloc[1,]
                leftJuncSupp = int(leftJunc.thickStart)
                rightJuncSupp = int(rightJunc.thickStart)
                leftJuncInfo = "{}-{}:{}".format(leftJunc.start, leftJunc.end, leftJuncSupp)
                rightJuncInfo = "{}-{}:{}".format(rightJunc.start, rightJunc.end, rightJuncSupp)

                tmpJuncCoverDf = juncCoverDfGrouped.get_group(juncName)
                scoreResetIndex = tmpJuncCoverDf.loc[:, "thickEnd"].reset_index(name="count")
                if len(scoreResetIndex.query("count>=3"))/float(len(tmpJuncCoverDf)) > 0.95 and leftJuncSupp >= 3 and rightJuncSupp >= 3:
                    juncIntronWithReadsSupp = "yes"
                    if np.std(tmpJuncCoverDf["thickEnd"]):
                        tmpJuncCoverDf.loc[:, "zscore"] = stats.zscore(tmpJuncCoverDf.loc[:, "thickEnd"])
                        tmpJuncCoverDf = tmpJuncCoverDf[abs(tmpJuncCoverDf["zscore"]) <= 3]
                        juncBaseCover = tmpJuncCoverDf.loc[:, ["thickEnd"]]
                        juncCv = round(abs(np.std(tmpJuncCoverDf.loc[:, "thickEnd"]) / np.mean(tmpJuncCoverDf.loc[:, "thickEnd"])), 3)
                        if juncCv > 0.2:
                            estimator = KMeans(n_clusters=2).fit(juncBaseCover)
                            juncIntronLowCoverCount = int(np.ceil(np.min(estimator.cluster_centers_)))
                        else:
                            juncIntronEvenlyDist = "yes"
                            juncIntronLowCoverCount = int(np.median(tmpJuncCoverDf.loc[:, "thickEnd"]))
                    else:
                        juncIntronEvenlyDist = "yes"
                        juncIntronLowCoverCount = int(np.median(tmpJuncCoverDf.loc[:, "thickEnd"]))

                    juncFPKM = round(1e9*juncIntronLowCoverCount/(len(tmpJuncCoverDf)*allReadCounts), 3)
                    if juncIntronLowCoverCount <= min([leftJuncSupp, rightJuncSupp]) * 1.5:
                        if juncFPKM >= 0.1 or juncIntronLowCoverCount >= 10:
                            juncIntronRetion = "yes"
                    else:
                        juncIntronRetion = "no"
                print >> out, "\t".join(listInt2Str(
                    [str(juncInfoDict[juncName]), leftJuncInfo, rightJuncInfo, juncIntronWithReadsSupp, juncIntronLowCoverCount,
                     juncIntronEvenlyDist, juncCv, juncFPKM, juncIntronRetion]))
            else:
                leftJunc = juncFlankCoverAllDfGrouped.get_group(juncName).iloc[0,]
                rightJunc = juncFlankCoverAllDfGrouped.get_group(juncName).iloc[1,]
                leftJuncSupp = int(leftJunc.thickStart)
                rightJuncSupp = int(rightJunc.thickStart)
                leftJuncInfo = "{}-{}:{}".format(leftJunc.start, leftJunc.end, leftJuncSupp)
                rightJuncInfo = "{}-{}:{}".format(rightJunc.start, rightJunc.end, rightJuncSupp)
                print >> out, "\t".join(listInt2Str(
                    [str(juncInfoDict[juncName]), leftJuncInfo, rightJuncInfo, juncIntronWithReadsSupp, juncIntronLowCoverCount,
                     juncIntronEvenlyDist, juncCv, juncFPKM, juncIntronRetion]))
    out.close()
    bam.close()

def removeFiles(myDir, fileList):
    for f in fileList:
        os.remove(os.path.join(myDir, f.strip("\n")))

def mergeFiles(fileList, refGtf=None, outFile="mergedDefault.txt", mode="intersect", sort=True, sortByCol=0):
    mergedLineList = []
    if refGtf:
        cmd = "hisat2_extract_splice_sites.py {} > tmp.ss".format(refGtf)
        subprocess.call(cmd, shell=True)
        fileList.append("tmp.ss")
    if mode == "intersect":
        for line in itertools.chain.from_iterable(itertools.imap(open, fileList)):
            mergedLineList.append(line.strip("\n"))
        mergedLineList = list(set(mergedLineList))
    elif mode == "union":
        mergedLineSet = set(open(fileList[0]).readlines())
        for i in fileList[1:]:
            mergedLineSet = mergedLineSet & set(open(i).readlines())
        mergedLineList = list(mergedLineSet)
    if sort:
        mergedLineList = sorted(mergedLineList, key=lambda line: int(line.strip("\n").split("\t")[sortByCol]))
    out = open(outFile, "w")
    for line in mergedLineList:
        print >> out, line
    out.close()

def getJuncCountFromNGSassembly(assembledGTF, bamFile, allReadCounts=5000000, outFile=None, threads=None, filterCount=None, generateTrack=False, juncList=None):
    if juncList:
        mergeFiles(juncList, refGtf=assembledGTF, outFile="tmp.ss", sortByCol=1, mode="intersect")
    else:
        cmd = "hisat2_extract_splice_sites.py {} > tmp.ss".format(assembledGTF)
        subprocess.call(cmd, shell=True)
    cmd = '''awk '{close(f);f=$1}{print > "chr"f".bed"}' tmp.ss'''
    subprocess.call(cmd, shell=True)
    processNum = threads if threads else psutil.cpu_count() * 3 / 4
    pool = Pool(processes=processNum)

    # bedFiles = os.popen("ls chr*")
    bedFiles = [i for i in os.listdir(".") if i.startswith("chr")]
    # random.shuffle(bedFiles)
    juncBedList = []
    chrBedList = []
    bamList = []
    baiList = []

    for i in bedFiles:
        chrId = os.path.splitext(i)[0].split("chr")[1]
        singleBam = "{}.bam".format(chrId)
        cmd = "samtools view -b {} {} -o {} && samtools index {}".format(bamFile, chrId, singleBam, singleBam)
        subprocess.call(cmd, shell=True)
        bamList.append(singleBam)
        baiList.append(singleBam + ".bai")

    for singleBed in bedFiles:
        juncBed = os.path.splitext(singleBed)[0] + ".junctions.bed"
        chrId = os.path.splitext(singleBed)[0].split("chr")[1]
        singleBam = "{}.bam".format(chrId)
        # pool.apply_async(parallelJuncCount, (singleBed, bamFile, juncBed, allReadCounts, filterCount))
        pool.apply_async(parallelJuncCount2, (singleBed, singleBam, juncBed, allReadCounts, filterCount))
        juncBedList.append(juncBed)
        chrBedList.append(singleBed.strip("\n"))
    pool.close()
    pool.join()

    with open(outFile, "w") as out:
        if generateTrack:
            out.write('#track name=junctions description="iFLAS TopHat2-like junctions"\n')
        for line in itertools.chain.from_iterable(itertools.imap(open, juncBedList)):
            out.write(line)
    removeFiles(os.getcwd(), juncBedList)
    removeFiles(os.getcwd(), chrBedList)
    removeFiles(os.getcwd(), bamList)
    removeFiles(os.getcwd(), baiList)

def main():
    hisat2mappingLogList = ["../../log/B73_Ki11_endo.repeat0.hisat2.log", "../../log/B73_Ki11_endo.repeat1.hisat2.log"]
    juncList = ["../alignment/repeat0/repeat0.ss", "../alignment/repeat1/repeat1.ss"]
    allMappedReadCounts = sum([getCountsFromHisat2MappingLog(i) for i in hisat2mappingLogList])
    getJuncCountFromNGSassembly("stringtie_merged.gtf", "tmp.bam", allReadCounts=allMappedReadCounts,
                                outFile="junctions.bed12", threads=10, filterCount=3, juncList=juncList)

if __name__ == '__main__':
    main()

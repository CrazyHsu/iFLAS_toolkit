#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: commonFuncs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-07-19 16:10:36
Last modified: 2019-07-19 16:10:37
'''
import re, os, subprocess, psutil, datetime
import venn, itertools
import numpy as np
import pandas as pd
import pysam, array, random, gzip, pybedtools
from sklearn.cluster import KMeans
from scipy import stats
import multiprocessing
from multiprocessing import Pool
from Bio import SeqIO, Seq, SeqRecord
from collections import Counter
# from commonObjs import Gene2Reads

#=============== Global variables =================
AStypes = ["SE", "IR", "A3SS", "A5SS"]

class TopHatJunction(object):
    def __init__(self, line=""):
        if line:
            self.record = line.strip("\n").split("\t")
            self.chrom = self.record[0]
            self.start = int(self.record[1])
            self.end = int(self.record[2])
            self.name = self.record[3]
            self.juncSupport = int(self.record[4])
            self.strand = self.record[5]
            self.blockSizes = listStr2Int(self.record[10].split(","))
            self.relStarts = listStr2Int(self.record[11].split(","))
            self.juncLeft = self.start + self.blockSizes[0]
            self.juncRight = self.start + self.relStarts[1]
        else:
            self.empty()

    def empty(self):
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

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

#=============== Functions=========================
def bam2fasta(bamFile, fastaFile):
    with open(fastaFile, "w") as outHandle:
        SeqIO.write(bam2rec(bamFile), outHandle, "fasta")

def bam2rec(inFile):
    bamFile = pysam.Samfile(inFile, "rb", check_sq=False)
    for read in bamFile:
        seq = Seq.Seq(str(read.seq))
        rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
        yield rec

def samProcess(samFile, isBam=False, outPrefix="out", name2pass=False, toFq=False, toFa=False):
    if isBam:
        bamHandle = pysam.AlignmentFile(samFile, "rb", check_sq=False)
    else:
        bamHandle = pysam.AlignmentFile(samFile, "r", check_sq=False)
    name2passDict = {}
    if toFq and toFa:
        outFq = open(outPrefix + ".fq", "w")
        outFa = open(outPrefix + ".fa", "w")
        for line in bamHandle:
            if name2pass:
                name2passDict[line.query_name] = line.get_tag("np")
            print >> outFq, "@{}\n{}\n+\n{}".format(line.query_name, line.query_sequence, pysam.qualities_to_qualitystring(line.query_qualities))
            print >> outFa, ">{}\n{}".format(line.query_name, line.query_sequence)
        outFa.close()
        outFq.close()
    elif toFa:
        outFa = open(outPrefix + ".fa", "w")
        for line in bamHandle:
            if name2pass:
                name2passDict[line.query_name] = line.get_tag("np")
            print >> outFa, ">{}\n{}".format(line.query_name, line.query_sequence)
            # print >> outFa, "{}".format(line.query_sequence)
        outFa.close()
    elif toFq:
        outFq = open(outPrefix + ".fq", "w")
        for line in bamHandle:
            if name2pass:
                name2passDict[line.query_name] = line.get_tag("np")
            print >> outFq, "@{}\n{}\n+\n{}".format(line.query_name, line.query_sequence,
                                                    pysam.qualities_to_qualitystring(line.query_qualities))
        outFq.close()

    if name2pass:
        name2passOut = open(outPrefix + ".name2pass.tsv", "w")
        for i in name2passDict:
            print >>name2passOut, "{}\t{}".format(i, name2passDict[i])
        name2passOut.close()
    bamHandle.close()

def batchMakeDir(dirList):
    for i in dirList:
        if not os.path.exists(i):
            os.makedirs(i)

def changeBamByFq(bamFile, fqFile, outBam):
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    fqRec = {}
    out = pysam.AlignmentFile(outBam, "wb", template=bam)
    with open(fqFile) as f:
        for rec in SeqIO.parse(f, "fastq"):
            fqRec[rec.id] = rec.seq
    for line in bam:
        line.query_sequence = str(fqRec[line.query_name])
        line.query_qualities = array.array('B', [20]*len(line.query_sequence))
        line.set_tag("qe", len(line.query_sequence), value_type="i")
        out.write(line)
    bam.close()
    out.close()

def changeUnclusterFastaName(clusterFile, flncFa, outFile, outClusterReport):
    "full_length_coverage=num;length=lec"
    cluster = pd.read_csv(clusterFile)
    out = open(outFile, "w")
    fastaOut = SeqIO.FastaIO.FastaWriter(out, wrap=None)
    outCluster = open(outClusterReport, "w")
    flncRecords = SeqIO.index(flncFa, "fasta")
    unclusteredList = list(set(flncRecords.keys()).difference(set(cluster.iloc[:,1].values)))
    transNum = len(cluster.iloc[:, 0].drop_duplicates())
    count = 0
    tmpRecordList = []
    for i in unclusteredList:
        tmpRecord = flncRecords[i]
        print >>outCluster, ",".join(["transcript/"+str(transNum+count), tmpRecord.id, "FL"])
        tmpRecord.id = "transcript/" + str(transNum + count) + " full_length_coverage=1;length=" + str(len(tmpRecord.seq))
        tmpRecord.description = ""
        count += 1
        tmpRecordList.append(tmpRecord)
    fastaOut.write_file(tmpRecordList)
    out.close()
    outCluster.close()

def checkArgs(gpeFile=None, readsFile=None, AStype=None, offset=None):
    if gpeFile: validateFile(gpeFile)
    if readsFile: validateFile(readsFile)
    assert AStype in AStypes
    assert str(offset).isdigit()

def checkFast5Files(fast5Dir):
    if os.path.isfile(fast5Dir):
        if fast5Dir.endswith("fq") or fast5Dir.endswith("fastq"):
            return "fastq"
    else:
        fileCount = len([x for x in os.listdir(fast5Dir) if not x.endswith(".txt")])
        fast5FileCount = 0
        fqFileCount = 0
        for i in os.listdir(fast5Dir):
            if i.endswith("fast5"):
                fast5FileCount += 1
            if i.endswith("fq") or i.endswith("fastq"):
                fqFileCount += 1
        if fast5FileCount == fileCount:
            return "fast5"
        if fqFileCount == fileCount:
            return "fastq"
        if fast5FileCount != fileCount and fqFileCount != fileCount:
            raise Exception("Please check the {} directory which is not a pure fast5 or fastq-containing directory".format(fast5Dir))

def checkHisat2IndexExist(dir):
    if len(os.listdir(dir)) == 0:
        return False
    else:
        ht2Flag = 0
        for i in os.listdir(dir):
            if i.endswith("exon") or i.endswith("ss") or i.endswith("log"):
                continue
            elif i.endswith("ht2"):
                ht2Flag = 1
                continue
            else:
                return False
        if ht2Flag:
            return True
        else:
            return False

def clusteredBamAddSingleSeq(clusterFile, flncFa, clusteredBam, mergedPolishBam):
    cmd = "sed '1d' {} | cut -f 1 -d ' ' | seqkit grep -n -v -f - -w 0 {} > unclustered.fa".format(clusterFile, flncFa)
    subprocess.call(cmd, shell=True)
    outSam = "uncluster.sam"
    clusterBamRead = pysam.AlignmentFile(clusteredBam, "rb", check_sq=False)
    firstLine = clusterBamRead.next()
    clusterCount = int(pysam.idxstats(clusteredBam).strip("\n").split("\t")[-1])
    out = open(outSam, "w")
    clusterBamRead.close()
    with open("unclustered.fa") as f:
        for rec in SeqIO.parse(f, "fasta"):
            query_name = "transcript/"+str(clusterCount)
            flag = 4
            rname = "*"
            pos = 0
            mapq = 255
            cigar = "*"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = str(rec.seq)
            qual = "*"
            tags = "\t".join([("ib:Z:0,1,1"), ("im:Z:"+str(rec.id)), ("is:i:1"), ("zm:i:"+str(clusterCount)),
                              ("RG:Z:"+firstLine.get_tag("RG"))])
            line = "\t".join(map(str, [query_name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags]))
            print >>out, line
            clusterCount += 1
    out.close()
    cmd = "samtools view -h {}| cat - uncluster.sam | samtools view -bS > {}".format(clusteredBam, mergedPolishBam)
    subprocess.call(cmd, shell=True)

def filter(originFile=None, targetFile=None, originField=1, targetField=1, mode="i", outFile=None, returnFlag=False):
    originList, targetList = [], []
    if isinstance(originFile, list):
        originList = originFile
    elif os.path.isfile(originFile):
        with open(originFile) as f:
            originList = [line.strip("\n").split("\t")[originField - 1] for line in f.readlines()]

    returnList = []
    outFile = outFile if outFile else "filtered.txt"
    out = open(outFile, "w")
    if isinstance(targetFile, list):
        targetList = targetFile
        if mode == "i":
            returnList = list(set(originList) & set(targetList))
        else:
            returnList = list(set(targetList) - set(originList))

        for i in returnList:
            print >>out, i
    elif isinstance(targetFile, dict):
        if mode == "i":
            returnList = list(set(originList) & set(targetFile.keys()))
        else:
            returnList = list(set(targetFile.keys()) - set(originList))

        for i in returnList:
            print >>out, i
    elif os.path.isfile(targetFile):
        with open(targetFile) as f:
            targetDict = {}
            for line in f:
                lineInfo = line.strip("\n").split("\t")
                targetDict[lineInfo[targetField-1]] = line.strip("\n")
            if mode == "i":
                returnList = list(set(originList) & set(targetDict.keys()))
                # for line in f:
                #     lineInfo = line.strip("\n").split("\t")
                #     if lineInfo[targetField-1] in originList:
                #         print >>out, line.strip("\n")
                #         returnList.append(line.strip("\n"))
            else:
                returnList = list(set(targetDict.keys()) - set(originList))
                # for line in f:
                #     lineInfo = line.strip("\n").split("\t")
                #     if lineInfo[targetField-1] not in originList:
                #         print >>out, line.strip("\n")
                #         returnList.append(line.strip("\n"))
            for i in returnList:
                print >> out, targetDict[i]
    out.close()
    if returnFlag:
        return returnList

def fusionTransAnalysis():
    print str(datetime.datetime.now()) + " Find fusion transcripts..."
    prevDir = os.getcwd()
    resolveDir("fusionTrans")
    inputFa = "../aln.sorted.fa"
    inputSam = "../aln.sorted.sam"
    cmd = "fusion_finder.py --input {} -s {} -o aln.fusionTrans.txt".format(inputFa, inputSam)
    # subprocess.call(cmd, shell=True)
    print str(datetime.datetime.now()) + " Find fusion transcripts done!"
    os.chdir(prevDir)

def generateRMTAScompFile(targetDir, sampleNames, outFile):
    bamStr = ','.join([os.path.join(targetDir, i, i+".sorted.bam") for i in sampleNames])
    out = open(outFile, "w")
    print >>out, bamStr
    out.close()

def getSampleCfg(sampleCfgFile):
    sampleDict = {}
    with open(sampleCfgFile) as f:
        for line in f:
            if line.startswith("#"): continue
            info = line.strip().split("\t")
            sampleDict[info[0]] = info[1]
    return sampleDict

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

def getAttribute(key, default, **args):
    return default if key not in args else args[key]

def getAbsStarts(start, blockRelStarts):
    return map(lambda x: start + x, blockRelStarts)

def getAbsLoc(start, blockSizes, blockRelStarts):
    '''
    get absolute location of the exons in transcripts through: start(first base location), blockSizes(length of each \
    exon), blockRelStarts(the relative position related to the first base, i.e. start)
    :returns: blockStarts and blockEnds, the list of starts of each starts and ends
    '''
    blockStarts = getAbsStarts(start, blockRelStarts)
    blockEnds = map(lambda x: blockStarts[x] + blockSizes[x], range(len(blockSizes)))
    return blockStarts, blockEnds


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

def getConsensusIntronNfromAS(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
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
                    if intronStarts2[k] > intronEnds1[i]: break
                    if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronStarts1[i] + offset and \
                            intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                        consensusN += 1
                        j += 1
            totalCnsN += consensusN
        elif asType in ["A3SS", "A5SS"]:
            cnsIntronStartN = getCnsSite(intronStarts1, intronStarts2, offset)
            cnsIntronEndN = getCnsSite(intronEnds1, intronEnds2, offset)
            totalCnsN += cnsIntronStartN + cnsIntronEndN
    return totalCnsN

def getConsensusIntronN1(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    j, consensusN = 0, 0
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronStarts1[i]: break
            if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
                    intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                consensusN += 1
                j += 1
    return consensusN
    # for asType in AStypes:
    #     if asType in ["IR", "SE"]:
    #         j, consensusN = 0, 0
    #         for i in range(len(intronStarts1)):
    #             for k in range(j, len(intronStarts2)):
    #                 if intronStarts2[k] > intronStarts1[i]: break
    #                 if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronEnds1[i] + offset and \
    #                         intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
    #                     consensusN += 1
    #                     j += 1
    #         totalCnsN += consensusN
    #     elif asType in ["A3SS", "A5SS"]:
    #         cnsIntronStartN = getCnsSite(intronStarts1, intronStarts2, offset)
    #         cnsIntronEndN = getCnsSite(intronEnds1, intronEnds2, offset)
    #         totalCnsN += cnsIntronStartN + cnsIntronEndN
    # return totalCnsN

def getDictFromFile(myFile, sep="\t", inlineSep=None, keyCol=1, valueCol=None):
    with open(myFile) as f:
        myDict = {}
        for line in f.readlines():
            infoList = line.strip("\n").split(sep)
            key = infoList[keyCol-1]
            if valueCol:
                if inlineSep:
                    value = infoList[valueCol-1].split(inlineSep)
                else:
                    value = infoList[valueCol-1]
            else:
                value = infoList
            myDict[key] = value
        return myDict

def getFxSequenceId(fxFile, isFa=False, isFq=False, outFile=None):
    recordId = []
    if isFa:
        for rec in SeqIO.parse(fxFile, "fasta"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId
    if isFq:
        for rec in SeqIO.parse(fxFile, "fastq"):
            recordId.append(rec.id)
        if outFile:
            with open(outFile, 'w') as f:
                for item in recordId:
                    print >> f, item
        else:
            return recordId

def getJuncCountFromNGSassembly(assembledGTF, bamFile, allReadCounts=5000000, outFile=None, threads=None, filterCount=None, generateTrack=False, juncList=None):
    if juncList:
        mergeFiles(juncList, refGtf=assembledGTF, outFile="tmp.ss", sortByCol=1, mode="intersect")
    else:
        cmd = "hisat2_extract_splice_sites.py {} > tmp.ss".format(assembledGTF)
        subprocess.call(cmd, shell=True)
    cmd = '''awk '{close(f);f=$1}{print > "chr_"f".bed"}' tmp.ss'''
    subprocess.call(cmd, shell=True)
    processNum = threads if threads else psutil.cpu_count() * 3 / 4
    pool = Pool(processes=processNum)

    # bedFiles = os.popen("ls chr*")
    bedFiles = [i for i in os.listdir(".") if i.startswith("chr_")]
    # random.shuffle(bedFiles)
    juncBedList = []
    chrBedList = []
    bamList = []
    baiList = []

    for i in bedFiles:
        chrId = os.path.splitext(i)[0].split("chr_")[1]
        singleBam = "{}.bam".format(chrId)
        cmd = "samtools view -b {} {} -o {} && samtools index {}".format(bamFile, chrId, singleBam, singleBam)
        subprocess.call(cmd, shell=True)
        bamList.append(singleBam)
        baiList.append(singleBam + ".bai")

    for singleBed in bedFiles:
        juncBed = os.path.splitext(singleBed)[0] + ".junctions.bed"
        chrId = os.path.splitext(singleBed)[0].split("chr_")[1]
        singleBam = "{}.bam".format(chrId)
        # pool.apply_async(parallelJuncCount, (singleBed, bamFile, juncBed, allReadCounts, filterCount))
        pool.apply_async(parallelJuncCount2, (singleBed, singleBam, juncBed, allReadCounts, filterCount))
        juncBedList.append(juncBed)
        chrBedList.append(singleBed.strip("\n"))
    pool.close()
    pool.join()

    # pool = Pool(processes=processNum)
    #
    # # bedFiles = os.popen("ls chr*")
    # bedFiles = [i for i in os.listdir(".") if i.startswith("chr")]
    # # random.shuffle(bedFiles)
    # juncBedList = []
    # chrBedList = []
    # for i in bedFiles:
    #     chrId = os.path.splitext(i)[0].split("chr")[1]
    #     singleBam = "{}.bam".format(chrId)
    #     cmd = "samtools view -b {} {} -o {} && samtools index {}".format(bamFile, chrId, singleBam, singleBam)
    #     subprocess.call(cmd, shell=True)
    #
    # for singleBed in ["chr1.bed"]:
    #     juncBed = os.path.splitext(singleBed)[0] + ".junctions.bed"
    #     chrId = os.path.splitext(singleBed)[0].split("chr")[1]
    #     singleBam = "{}.bam".format(chrId)
    #     # pool.apply_async(parallelJuncCount, (singleBed, bamFile, juncBed, allReadCounts, filterCount))
    #     pool.apply_async(parallelJuncCount, (singleBed, singleBam, juncBed, allReadCounts, filterCount))
    #     # pool.apply_async(parallelJuncCount1, (singleBed, singleBam, juncBed, allReadCounts, filterCount))
    #     # parallelJuncCount1(singleBed, singleBam, juncBed, allReadCounts, filterCount)
    #     juncBedList.append(juncBed)
    #     chrBedList.append(singleBed.strip("\n"))
    # pool.close()
    # pool.join()

    # juncSuppOut = open("juncSupp.txt", "w")
    # for j in multiResults:
    #     leftJuncSupp, rightJuncSupp = j.get()
    #     for i in range(len(leftJuncSupp)):
    #         leftInfo = leftJuncSupp[i].split("\t")
    #         rightInfo = rightJuncSupp[i].split("\t")
    #         print >> juncSuppOut, "".join([leftInfo[3], leftInfo[6], rightInfo[6]])
    #         # print >> juncSuppOut, "\t".join(leftInfo) + "\t" + "\t".join(rightInfo)
    # juncSuppOut.close()
    with open(outFile, "w") as out:
        if generateTrack:
            out.write('#track name=junctions description="iFLAS TopHat2-like junctions"\n')
        for line in itertools.chain.from_iterable(itertools.imap(open, juncBedList)):
            out.write(line)
    removeFiles(os.getcwd(), juncBedList)
    removeFiles(os.getcwd(), chrBedList)
    removeFiles(os.getcwd(), bamList)
    removeFiles(os.getcwd(), baiList)
    # removeFiles(os.getcwd(), ["tmp.ss"])

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

def getSubSamByName(samFile, nameList, isBam=False, nameListIsFile=False, outPrefix=None, sort=True, threads=4):
    # n = get_names(options.names)
    n = nameList
    if nameListIsFile:
        with open(nameList, 'r') as infile:
            n = infile.read().splitlines()
    if isBam:
        samHandle = pysam.AlignmentFile(samFile, 'rb', check_sq=False)
        suffix = ".bam"
    else:
        tmpBam = "tmp.bam"
        pysam.view("-o", tmpBam, "--output-fmt", "BAM", samFile, catch_stdout=False)
        samHandle = pysam.AlignmentFile(tmpBam, 'rb', check_sq=False)
        suffix = ".sam"
    name_indexed = pysam.IndexedReads(samHandle)
    name_indexed.build()
    header = samHandle.header.copy()
    if isBam:
        out = pysam.Samfile(outPrefix + suffix, 'wb', header=header)
    else:
        out = pysam.Samfile(outPrefix + suffix, 'w', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)
    out.close()
    if not isBam:
        os.remove("tmp.bam")
    if sort:
        if isBam:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), outPrefix + suffix, catch_stdout=False)
        else:
            pysam.sort("-o", outPrefix + ".sorted" + suffix, "-@", str(threads), "--output-fmt", "SAM", outPrefix + suffix, catch_stdout=False)
        os.remove(outPrefix + suffix)
        return outPrefix + ".sorted" + suffix
    else:
        return outPrefix + suffix

def isOverlap(listA, listB):
    sortedList = sorted([listA, listB])
    return True if sortedList[0][1] >= sortedList[1][0] else False

def initRefSetting(refParams=None, dirSpec=None):
    originDir = os.getcwd()
    if refParams.ref_genome == None:
        raise Exception("You must specify a genome file!")
    refParams.ref_genome = os.path.abspath(os.path.join(originDir, refParams.ref_genome))
    validateFile(refParams.ref_genome)

    if refParams.ref_gtf == None and refParams.ref_gff == None:
        raise Exception("You must specify a GTF or GFF3 file!")
    elif refParams.ref_gtf == None and refParams.ref_gff != None:
        refParams.ref_gff = os.path.abspath(os.path.join(originDir, refParams.ref_gff))
        validateFile(refParams.ref_gff)
        gffDir = os.path.dirname(refParams.ref_gff)
        ref_gtf = os.path.abspath(os.path.join(gffDir, "iFLAS.gtf"))
        cmd = "gffread -T -o {} {}".format(ref_gtf, refParams.ref_gff)
        subprocess.call(cmd, shell=True)
        refParams.ref_gtf = ref_gtf

    gtfDir = os.path.dirname(refParams.ref_gtf)
    if refParams.ref_size == None:
        refParams.ref_size = os.path.join(gtfDir, "iFLAS.len.txt")
        sizeOut = open(refParams.ref_size)
        for seqRecord in SeqIO.parse(refParams.ref_genome, "fasta"):
            print >> sizeOut, "\t".join([seqRecord.id, seqRecord.seq])
        sizeOut.close()
    else:
        refParams.ref_size = os.path.join(originDir, refParams.ref_size)

    if refParams.ref_gpe == None:
        refParams.ref_gpe = os.path.join(gtfDir, "iFLAS.gpe")
        cmd = "gtfToGenePred {} {} -genePredExt".format(refParams.ref_gtf, refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_gpe = os.path.abspath(os.path.join(originDir, refParams.ref_gpe))

    if refParams.ref_bed == None:
        refParams.ref_bed = os.path.join(gtfDir, "iFLAS.annotation.bed")
        cmd = "gpe2bed.pl -g {} > {}".format(refParams.ref_gpe, refParams.ref_bed)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_bed = os.path.abspath(os.path.join(originDir, refParams.ref_bed))
    refParams.ref_mm2_index = os.path.abspath(os.path.join(originDir, refParams.ref_mm2_index))

    if dirSpec.tmp_dir == None:
        dirSpec.tmp_dir = "/tmp"
    else:
        if not os.path.exists(dirSpec.tmp_dir):
            batchMakeDir([dirSpec.tmp_dir])
    dirSpec.out_dir = os.path.abspath(os.path.join(originDir, dirSpec.out_dir))
    if not os.path.exists(dirSpec.out_dir):
        batchMakeDir([dirSpec.out_dir])

def initResourceSetting(optionTools=None):
    maxMem = psutil.virtual_memory().total >> 20
    if str(optionTools.memory).endswith("M"):
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif str(optionTools.memory).endswith("G"):
        optionTools.memory = str(int(optionTools.memory[0:-1]) * 1024) + "M"
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif optionTools.memory.isdigit():
        if int(optionTools.memory) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    else:
        raise Exception("Please check the format of memory you specified in optionTool section!")

    maxThreads = psutil.cpu_count()
    if str(optionTools.threads).isdigit():
        if int(optionTools.threads) > maxThreads * 0.75:
            optionTools.threads = int(maxThreads * 0.75)
        else:
            optionTools.threads = int(optionTools.threads)
    else:
        raise Exception("The number of threads should be digital, please check it!")

def initSysSetting(refParams=None, optionTools=None, dirSpec=None):
    # global originDir
    originDir = os.getcwd()
    if refParams.ref_genome == None:
        raise Exception("You must specify a genome file!")
    refParams.ref_genome = os.path.abspath(os.path.join(originDir, refParams.ref_genome))
    validateFile(refParams.ref_genome)

    if refParams.ref_gtf == None and refParams.ref_gff == None:
        raise Exception("You must specify a GTF or GFF3 file!")
    elif refParams.ref_gtf == None and refParams.ref_gff != None:
        refParams.ref_gff = os.path.abspath(os.path.join(originDir, refParams.ref_gff))
        validateFile(refParams.ref_gff)
        gffDir = os.path.dirname(refParams.ref_gff)
        ref_gtf = os.path.abspath(os.path.join(gffDir, "iFLAS.gtf"))
        cmd = "gffread -T -o {} {}".format(ref_gtf, refParams.ref_gff)
        subprocess.call(cmd, shell=True)
        refParams.ref_gtf = ref_gtf

    gtfDir = os.path.dirname(refParams.ref_gtf)
    if refParams.ref_size == None:
        refParams.ref_size = os.path.join(gtfDir, "iFLAS.len.txt")
        sizeOut = open(refParams.ref_size)
        for seqRecord in SeqIO.parse(refParams.ref_genome, "fasta"):
            print >>sizeOut, "\t".join([seqRecord.id, seqRecord.seq])
        sizeOut.close()
    else:
        refParams.ref_size = os.path.join(originDir, refParams.ref_size)

    if refParams.ref_gpe == None:
        refParams.ref_gpe = os.path.join(gtfDir, "iFLAS.gpe")
        cmd = "gtfToGenePred {} {} -genePredExt".format(refParams.ref_gtf, refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_gpe = os.path.abspath(os.path.join(originDir, refParams.ref_gpe))

    if refParams.ref_bed == None:
        refParams.ref_bed = os.path.join(gtfDir, "iFLAS.annotation.bed")
        cmd = "gpe2bed.pl -g {} > {}".format(refParams.ref_gpe, refParams.ref_bed)
        subprocess.call(cmd, shell=True)
    else:
        refParams.ref_bed = os.path.abspath(os.path.join(originDir, refParams.ref_bed))
    refParams.ref_mm2_index = os.path.abspath(os.path.join(originDir, refParams.ref_mm2_index))

    if dirSpec.tmp_dir == None:
        dirSpec.tmp_dir = "/tmp"
    else:
        if not os.path.exists(dirSpec.tmp_dir):
            batchMakeDir([dirSpec.tmp_dir])
    dirSpec.out_dir = os.path.abspath(os.path.join(originDir, dirSpec.out_dir))
    if not os.path.exists(dirSpec.tmp_dir):
        batchMakeDir([dirSpec.tmp_dir])

    maxMem = psutil.virtual_memory().total >> 20
    if str(optionTools.memory).endswith("M"):
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif str(optionTools.memory).endswith("G"):
        optionTools.memory = str(int(optionTools.memory[0:-1]) * 1024) + "M"
        if int(optionTools.memory[0:-1]) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    elif optionTools.memory.isdigit():
        if int(optionTools.memory) > maxMem * 0.75:
            optionTools.memory = str(maxMem * 0.75) + "M"
    else:
        raise Exception("Please check the format of memory you specified in optionTool section!")

    maxThreads = psutil.cpu_count()
    if str(optionTools.threads).isdigit():
        if int(optionTools.threads) > maxThreads * 0.75:
            optionTools.threads = int(maxThreads * 0.75)
        else:
            optionTools.threads = int(optionTools.threads)
    else:
        raise Exception("The number of threads should be digital, please check it!")

def listStr2Int(myList):
    return [int(i) for i in myList]

def listInt2Str(myList):
    return [str(i) for i in myList]

def listDigital2Str(myList):
    return [str(i) for i in myList]

def makeLink(sourcePath, targetPath):
    if os.path.islink(targetPath) or os.path.exists(targetPath):
        os.remove(targetPath)
    os.symlink(sourcePath, targetPath)

def makeHisat2Index(refParams=None, dirSpec=None, optionTools=None):
    hisat2indexDir = os.path.join(dirSpec.out_dir, "hisat2indexDir")
    if not os.path.exists(hisat2indexDir):
        os.makedirs(hisat2indexDir)

    if not checkHisat2IndexExist(hisat2indexDir):
        print str(datetime.datetime.now()) + " Start Hisat2 indexing..."
        batchMakeDir([hisat2indexDir])
        gtfPrefix = os.path.splitext(os.path.basename(refParams.ref_gtf))[0]
        refAnnoGTF = refParams.ref_gtf
        refAnnoSS = "{}/{}.ss".format(hisat2indexDir, gtfPrefix)
        cmd = "hisat2_extract_splice_sites.py {} >{}".format(refAnnoGTF, refAnnoSS)
        subprocess.call(cmd, shell=True)
        refAnnoExon = "{}/{}.exon".format(hisat2indexDir, gtfPrefix)
        cmd = "hisat2_extract_exons.py {} >{}".format(refAnnoGTF, refAnnoExon)
        subprocess.call(cmd, shell=True)
        cmd = "hisat2-build -p {} --ss {} --exon {} {} {}/{} 1>{}/hisat2build.log 2>&1".format(optionTools.threads,
                                                                                               refAnnoSS, refAnnoExon,
                                                                                               refParams.ref_genome,
                                                                                               hisat2indexDir,
                                                                                               gtfPrefix, hisat2indexDir)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " End Hisat2 indexing..."

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

def mergeTupleList1(tupleList):
    newList = []
    minPos = min([i[0] for i in tupleList])
    maxPos = max([i[1] for i in tupleList])
    posDict = dict.fromkeys(range(minPos, maxPos+1), 0)
    for t in tupleList:
        for p in range(t[0], t[1]+1):
            posDict[p] = 1
    flag = 0
    tmpMin = maxPos
    tmpMax = minPos
    for j in range(minPos, maxPos+1):
        if posDict[j] == 1:
            tmpMin = j if j < tmpMin else tmpMin
            tmpMax = j if j > tmpMax else tmpMax
            flag = 1
        else:
            if flag == 1:
                newList.append((tmpMin, tmpMax))
                tmpMin = maxPos
                tmpMax = minPos
                flag = 0
    newList.append((tmpMin, tmpMax))
    return newList

def merge2exonList1(listA, listB):
    return mergeTupleList1(listA + listB)

def merge2exonList(listA, listB):
    k = 0
    newList = []
    for i in range(len(listA)):
        for j in range(k, len(listB)):
            if listA[i][0] > listB[j][1]:
                continue
            if listA[i][1] < listB[j][0]:
                break
            exonBound = ()
            if listA[i][1] >= listB[j][0] and listA[i][1] <= listB[j][1]:
                exonBound = (min(listA[i][0], listB[j][0]), listB[j][1])
            elif listA[i][0] <= listB[j][1] and listA[i][0] >= listB[j][0]:
                exonBound = (listB[j][0], max(listA[i][1], listB[j][1]))
            elif listA[i][0] < listB[j][0] and listA[i][1] > listB[j][1]:
                exonBound = (listA[i][0], listA[i][1])
            elif listA[i][0] >= listB[j][0] and listA[i][1] <= listB[j][1]:
                exonBound = (listB[j][0], listB[j][1])
            if exonBound:
                newList.append(exonBound)
                print exonBound
        else:
            newList.append(listA[i])
    return mergeTupleList(newList)

def parallelJuncCount(bedFile, bamFile, outFile, allReadCounts=5000000, filterCount=None):
    out = open(outFile, "w")
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    with open(bedFile) as f:
        count = 1
        readsBamList = []
        for line in f:
            lineList = line.strip("\n").split("\t")
            chrom, start, end, strand = lineList[0], int(lineList[1]), int(lineList[2]), lineList[3]
            reads = bam.fetch(chrom, start, end)
            allReadsBam = "{}_{}_{}.bam".format(chrom, start, end)
            if allReadsBam in readsBamList: continue
            readsBamList.append(allReadsBam)
            allReadsOut = pysam.AlignmentFile(allReadsBam, "wb", template=bam)
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
                allReadsOut.write(read)
            allReadsOut.close()
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

            juncLeftExonToIntron = "{}\t{}\t{}\t{}".format(chrom, junc.juncLeft - 3, junc.juncLeft + 2, strand)
            juncRightIntronToExon = "{}\t{}\t{}\t{}".format(chrom, junc.juncRight - 2, junc.juncRight + 3, strand)
            allReadsStream = pybedtools.BedTool(allReadsBam).bam_to_bed(split=True)

            leftBedObj = pybedtools.BedTool(juncLeftExonToIntron, from_string=True)
            rightBedObj = pybedtools.BedTool(juncRightIntronToExon, from_string=True)
            leftCover = leftBedObj.coverage(allReadsStream, f=1)
            rightCover = rightBedObj.coverage(allReadsStream, f=1)
            leftJuncSupp = int(str(leftCover[0]).strip("\n").split("\t")[4])
            rightJuncSupp = int(str(rightCover[0]).strip("\n").split("\t")[4])
            leftJuncInfo = "{}-{}:{}".format(junc.juncLeft - 3, junc.juncLeft + 2, leftJuncSupp)
            rightJuncInfo = "{}-{}:{}".format(junc.juncRight - 2, junc.juncRight + 3, rightJuncSupp)

            juncBed = "{}\t{}\t{}".format(chrom, junc.juncLeft, junc.juncRight)
            juncBedObj = pybedtools.BedTool(juncBed, from_string=True)
            juncCoverDf = juncBedObj.coverage(allReadsStream, d=True).to_dataframe()

            scoreResetIndex = juncCoverDf.loc[:, "score"].reset_index(name="count")
            juncWithReadsSupp, juncEvenlyDist, juncLowCoverCount = "no", "no", 0
            juncCv, juncFPKM, juncIntronRetion = 0, 0, "no"
            if len(scoreResetIndex.query("count>=3"))/float(len(juncCoverDf)) > 0.95 and leftJuncSupp >= 3 and rightJuncSupp >= 3:
                juncWithReadsSupp = "yes"
                if np.std(juncCoverDf["score"]):
                    juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
                    juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
                    juncBaseCover = juncCoverDf.loc[:, ["score"]]
                    juncCv = round(abs(np.std(juncCoverDf.loc[:, "score"]) / np.mean(juncCoverDf.loc[:, "score"])), 3)
                    if juncCv > 0.2:
                        estimator = KMeans(n_clusters=2).fit(juncBaseCover)
                        juncLowCoverCount = int(np.ceil(np.min(estimator.cluster_centers_)))
                    else:
                        juncEvenlyDist = "yes"
                        juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))

                    juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
                    if juncFPKM >= 1:
                        juncIntronRetion = "yes"
                else:
                    juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))
                    juncEvenlyDist = "yes"
                    juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
                    juncIntronRetion = "yes"

            # juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
            # juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
            #
            # juncBaseCover = juncCoverDf.loc[:, ["score"]]
            # juncCv = abs(np.std(juncBaseCover)/np.mean(juncBaseCover))
            # if juncCv > 0.25:
            #     estimator = KMeans(n_clusters=2).fit(juncBaseCover)

            print >> out, "\t".join(listInt2Str([str(junc), leftJuncInfo, rightJuncInfo, juncWithReadsSupp, juncLowCoverCount, juncEvenlyDist, juncCv, juncFPKM, juncIntronRetion]))
            count += 1
            # removeFiles(os.getcwd(), [leftBam, rightBam])
        removeFiles(os.getcwd(), readsBamList)
    out.close()
    bam.close()

def parallelJuncCount1(bedFile, bamFile, outFile, allReadCounts=5000000, filterCount=None):
    out = open(outFile, "w")
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    allReadsStream = pybedtools.BedTool(bamFile).bam_to_bed(split=True)
    with open(bedFile) as f:
        count = 1
        readsBamList = []
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

            juncLeftExonToIntron = "{}\t{}\t{}\t{}".format(chrom, junc.juncLeft - 3, junc.juncLeft + 2, strand)
            juncRightIntronToExon = "{}\t{}\t{}\t{}".format(chrom, junc.juncRight - 2, junc.juncRight + 3, strand)


            leftBedObj = pybedtools.BedTool(juncLeftExonToIntron, from_string=True)
            rightBedObj = pybedtools.BedTool(juncRightIntronToExon, from_string=True)
            leftCover = leftBedObj.coverage(allReadsStream, f=1)
            rightCover = rightBedObj.coverage(allReadsStream, f=1)
            leftJuncSupp = int(str(leftCover[0]).strip("\n").split("\t")[4])
            rightJuncSupp = int(str(rightCover[0]).strip("\n").split("\t")[4])
            leftJuncInfo = "{}-{}:{}".format(junc.juncLeft - 3, junc.juncLeft + 2, leftJuncSupp)
            rightJuncInfo = "{}-{}:{}".format(junc.juncRight - 2, junc.juncRight + 3, rightJuncSupp)

            juncBed = "{}\t{}\t{}".format(chrom, junc.juncLeft, junc.juncRight)
            juncBedObj = pybedtools.BedTool(juncBed, from_string=True)
            juncCoverDf = juncBedObj.coverage(allReadsStream, d=True).to_dataframe()

            scoreResetIndex = juncCoverDf.loc[:, "score"].reset_index(name="count")
            juncWithReadsSupp, juncEvenlyDist, juncLowCoverCount = "no", "no", 0
            juncCv, juncFPKM, juncIntronRetion = 0, 0, "no"
            if len(scoreResetIndex.query("count>=3"))/float(len(juncCoverDf)) > 0.95 and leftJuncSupp >= 3 and rightJuncSupp >= 3:
                juncWithReadsSupp = "yes"
                if np.std(juncCoverDf["score"]):
                    juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
                    juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
                    juncBaseCover = juncCoverDf.loc[:, ["score"]]
                    juncCv = round(abs(np.std(juncCoverDf.loc[:, "score"]) / np.mean(juncCoverDf.loc[:, "score"])), 3)
                    if juncCv > 0.2:
                        estimator = KMeans(n_clusters=2).fit(juncBaseCover)
                        juncLowCoverCount = int(np.ceil(np.min(estimator.cluster_centers_)))
                    else:
                        juncEvenlyDist = "yes"
                        juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))

                    juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
                    if juncFPKM >= 1:
                        juncIntronRetion = "yes"
                else:
                    juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))
                    juncEvenlyDist = "yes"
                    juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
                    juncIntronRetion = "yes"

            # juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
            # juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
            #
            # juncBaseCover = juncCoverDf.loc[:, ["score"]]
            # juncCv = abs(np.std(juncBaseCover)/np.mean(juncBaseCover))
            # if juncCv > 0.25:
            #     estimator = KMeans(n_clusters=2).fit(juncBaseCover)

            print >> out, "\t".join(listInt2Str([str(junc), leftJuncInfo, rightJuncInfo, juncWithReadsSupp, juncLowCoverCount, juncEvenlyDist, juncCv, juncFPKM, juncIntronRetion]))
            count += 1
            # removeFiles(os.getcwd(), [leftBam, rightBam])
        removeFiles(os.getcwd(), readsBamList)
    out.close()
    bam.close()

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


    #         leftBedObj = pybedtools.BedTool(juncLeftExonToIntron, from_string=True)
    #         rightBedObj = pybedtools.BedTool(juncRightIntronToExon, from_string=True)
    #         leftCover = leftBedObj.coverage(allReadsStream, f=1)
    #         rightCover = rightBedObj.coverage(allReadsStream, f=1)
    #         leftJuncSupp = int(str(leftCover[0]).strip("\n").split("\t")[4])
    #         rightJuncSupp = int(str(rightCover[0]).strip("\n").split("\t")[4])
    #         leftJuncInfo = "{}-{}:{}".format(junc.juncLeft - 3, junc.juncLeft + 2, leftJuncSupp)
    #         rightJuncInfo = "{}-{}:{}".format(junc.juncRight - 2, junc.juncRight + 3, rightJuncSupp)
    #
    #         juncBed = "{}\t{}\t{}".format(chrom, junc.juncLeft, junc.juncRight)
    #         juncBedObj = pybedtools.BedTool(juncBed, from_string=True)
    #         juncCoverDf = juncBedObj.coverage(allReadsStream, d=True).to_dataframe()
    #
    #         scoreResetIndex = juncCoverDf.loc[:, "score"].reset_index(name="count")
    #         juncWithReadsSupp, juncEvenlyDist, juncLowCoverCount = "no", "no", 0
    #         juncCv, juncFPKM, juncIntronRetion = 0, 0, "no"
    #         if len(scoreResetIndex.query("count>=3"))/float(len(juncCoverDf)) > 0.95 and leftJuncSupp >= 3 and rightJuncSupp >= 3:
    #             juncWithReadsSupp = "yes"
    #             if np.std(juncCoverDf["score"]):
    #                 juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
    #                 juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
    #                 juncBaseCover = juncCoverDf.loc[:, ["score"]]
    #                 juncCv = round(abs(np.std(juncCoverDf.loc[:, "score"]) / np.mean(juncCoverDf.loc[:, "score"])), 3)
    #                 if juncCv > 0.2:
    #                     estimator = KMeans(n_clusters=2).fit(juncBaseCover)
    #                     juncLowCoverCount = int(np.ceil(np.min(estimator.cluster_centers_)))
    #                 else:
    #                     juncEvenlyDist = "yes"
    #                     juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))
    #
    #                 juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
    #                 if juncFPKM >= 1:
    #                     juncIntronRetion = "yes"
    #             else:
    #                 juncLowCoverCount = int(np.median(juncCoverDf.loc[:, "score"]))
    #                 juncEvenlyDist = "yes"
    #                 juncFPKM = round(1e9*juncLowCoverCount/(len(juncCoverDf)*allReadCounts), 3)
    #                 juncIntronRetion = "yes"
    #
    #         # juncCoverDf["zscore"] = stats.zscore(juncCoverDf.loc[:, "score"])
    #         # juncCoverDf = juncCoverDf[juncCoverDf["zscore"] <= 3]
    #         #
    #         # juncBaseCover = juncCoverDf.loc[:, ["score"]]
    #         # juncCv = abs(np.std(juncBaseCover)/np.mean(juncBaseCover))
    #         # if juncCv > 0.25:
    #         #     estimator = KMeans(n_clusters=2).fit(juncBaseCover)
    #
    #         print >> out, "\t".join(listInt2Str([str(junc), leftJuncInfo, rightJuncInfo, juncWithReadsSupp, juncLowCoverCount, juncEvenlyDist, juncCv, juncFPKM, juncIntronRetion]))
    #         count += 1
    #         # removeFiles(os.getcwd(), [leftBam, rightBam])
    #     removeFiles(os.getcwd(), readsBamList)
    # out.close()
    # bam.close()


def parseBlockChain(exonChain):
    '''
    parse block chain info, convert the exon chain, in the format of "122-144;344-455;567-778", to blockStarts and blockEnds
    :return:
    '''
    blockStarts, blockEnds = [], []
    for block in exonChain.split(";"):
        blockStart, blockEnd = block.split("-")
        blockStarts.append(int(blockStart))
        blockEnds.append(int(blockEnd))
    return blockStarts, blockEnds

def renameNGSdata2fastp(tgsSample=None):
    if tgsSample.ngsPaired == "paired":
        leftReadsRepeats = [i.strip() for i in tgsSample.ngsLeftReads.split(";")]
        rightReadsRepeats = [i.strip() for i in tgsSample.ngsRightReads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            newLeftReadsRepeats = []
            newRightReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                if len(leftReads) != len(rightReads):
                    raise Exception("please input the right paired NGS reads for the analysis")
                else:
                    newLeftReads = []
                    newRightReads = []
                    for j in range(len(leftReads)):
                        leftReadsDir = os.path.dirname(leftReads[j])
                        rightReadsDir = os.path.dirname(rightReads[j])
                        leftReadsBase = os.path.basename(leftReads[j])
                        rightReadsBase = os.path.basename(rightReads[j])
                        newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                        newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                        newLeftReads.append(newLeft)
                        newRightReads.append(newRight)
                    newLeftReadsRepeats.append(",".join(newLeftReads))
                    newRightReadsRepeats.append(",".join(newRightReads))
            tgsSample.ngsLeftReads = ";".join(newLeftReadsRepeats)
            tgsSample.ngsRightReads = ";".join(newRightReadsRepeats)
    else:
        if tgsSample.ngsLeftReads and tgsSample.ngsRightReads == None:
            leftReadsRepeats = [i.strip() for i in tgsSample.ngsLeftReads.split(";")]
            newLeftReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                newLeftReads = []
                for j in range(len(leftReads)):
                    leftReadsDir = os.path.dirname(leftReads[j])
                    leftReadsBase = os.path.basename(leftReads[j])
                    newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                    newLeftReads.append(newLeft)
                newLeftReadsRepeats.append(",".join(newLeftReads))
            tgsSample.ngsLeftReads = ";".join(newLeftReadsRepeats)
        elif tgsSample.ngsRightReads and tgsSample.ngsLeftReads == None:
            rightReadsRepeats = [i.strip() for i in tgsSample.ngsRightReads.split(";")]
            newRightReadsRepeats = []
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                newRightReads = []
                for j in range(len(rightReads)):
                    rightReadsDir = os.path.dirname(rightReads[j])
                    rightReadsBase = os.path.basename(rightReads[j])
                    newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                    newRightReads.append(newRight)
                newRightReadsRepeats.append(",".join(newRightReads))
            tgsSample.ngsRightReads = ";".join(newRightReadsRepeats)
        else:
            raise Exception("The NGS data seem not to be single, please check it")

def preprocessNGSdata(tgsSample=None):
    if tgsSample.ngsPaired == "paired":
        leftReadsRepeats = [i.strip() for i in tgsSample.ngsLeftReads.split(";")]
        rightReadsRepeats = [i.strip() for i in tgsSample.ngsRightReads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            newLeftReadsRepeats = []
            newRightReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                if len(leftReads) != len(rightReads):
                    raise Exception("please input the right paired NGS reads for the analysis")
                else:
                    newLeftReads = []
                    newRightReads = []
                    for j in range(len(leftReads)):
                        leftReadsDir = os.path.dirname(leftReads[j])
                        rightReadsDir = os.path.dirname(rightReads[j])
                        leftReadsBase = os.path.basename(leftReads[j])
                        rightReadsBase = os.path.basename(rightReads[j])
                        newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                        newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                        cmd = "fastp -i {} -I {} -o {} -O {} -w {} -q 20 -l {} 2>/dev/null"
                        cmd = cmd.format(leftReads[j], rightReads[j], newLeft, newRight, tgsSample.threads, int(tgsSample.ngsReadsLength)-30)
                        subprocess.call(cmd, shell=True)
                        newLeftReads.append(newLeft)
                        newRightReads.append(newRight)
                    newLeftReadsRepeats.append(",".join(newLeftReads))
                    newRightReadsRepeats.append(",".join(newRightReads))
            tgsSample.ngsLeftReads = ";".join(newLeftReadsRepeats)
            tgsSample.ngsRightReads = ";".join(newRightReadsRepeats)
    else:
        if tgsSample.ngsLeftReads and tgsSample.ngsRightReads == None:
            leftReadsRepeats = [i.strip() for i in tgsSample.ngsLeftReads.split(";")]
            newLeftReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                newLeftReads = []
                for j in range(len(leftReads)):
                    leftReadsDir = os.path.dirname(leftReads[j])
                    leftReadsBase = os.path.basename(leftReads[j])
                    newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                    cmd = "fastp -i {} -o {} -w {} -q 20 -l {} 2>/dev/null"
                    cmd.format(leftReads[j], newLeft, tgsSample.threads, int(tgsSample.ngsReadsLength)-30)
                    subprocess.call(cmd, shell=True)
                    newLeftReads.append(newLeft)
                newLeftReadsRepeats.append(",".join(newLeftReads))
            tgsSample.ngsLeftReads = ";".join(newLeftReadsRepeats)
        elif tgsSample.ngsRightReads and tgsSample.ngsLeftReads == None:
            rightReadsRepeats = [i.strip() for i in tgsSample.ngsRightReads.split(";")]
            newRightReadsRepeats = []
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                newRightReads = []
                for j in range(len(rightReads)):
                    rightReadsDir = os.path.dirname(rightReads[j])
                    rightReadsBase = os.path.basename(rightReads[j])
                    newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                    cmd = "fastp -i {} -o {} -w {} -q 20 -l {} 2>/dev/null"
                    cmd.format(rightReads[j], newRight, tgsSample.threads, int(tgsSample.ngsReadsLength)-30)
                    subprocess.call(cmd, shell=True)
                    newRightReads.append(newRight)
                newRightReadsRepeats.append(",".join(newRightReads))
            tgsSample.ngsRightReads = ";".join(newRightReadsRepeats)
        else:
            raise Exception("The NGS data seem not to be single, please check it")

def removeFiles(myDir, fileList):
    # print fileList
    for f in fileList:
        # print f, myDir
        # print os.path.join(myDir, f.strip("\n"))
        os.remove(os.path.join(myDir, f.strip("\n")))

def resizeTrackRatio(itemCounts):
    originRelativeSize = 10.0
    if itemCounts <= 10:
        return originRelativeSize
    elif 10 <= itemCounts and itemCounts <= 100:
        return 2 * originRelativeSize * np.log10(itemCounts) - originRelativeSize
    else:
        return originRelativeSize * 3

def resolveDir(dirName, chdir=True):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if chdir:
        os.chdir(dirName)

def sortTupleList(tupleList):
    return sorted(tupleList)

def validateFile(myFile):
    if not os.path.exists(myFile):
        raise Exception("File '%s' not found! Please input again!" % myFile)

    if not os.path.isfile(myFile):
        raise Exception("File '%s' is not a file! Please input again!" % myFile)

    return True

def validateDir(myDir):
    if not os.path.exists(myDir):
        raise Exception("Dir '%s' not found! Please input again!" % myDir)

    if not os.path.isdir(myDir):
        raise Exception("Dir '%s' is not a directory! Please input again!" % myDir)

    return True

def vennPlot(outFile, **args):
    if len(args) == 2:
        file1, file2 = args
        list1 = [i.strip("\n") for i in open(file1).readlines()]
        list2 = [i.strip("\n") for i in open(file2).readlines()]
        labels = venn.get_labels([list1, list2], fill=["number"])
        fig, ax = venn.venn2(labels, names=[file1, file2])
        fig.savefig(outFile)
        # fig.show()
    if len(args) == 3:
        file1, file2, file3 = args
        list1 = [i.strip("\n") for i in open(file1).readlines()]
        list2 = [i.strip("\n") for i in open(file2).readlines()]
        list3 = [i.strip("\n") for i in open(file3).readlines()]
        labels = venn.get_labels([list1, list2, list3], fill=['number'])
        fig, ax = venn.venn3(labels, names=[file1, file2, file3])
        fig.savefig(outFile)
        # fig.show()
    if len(args) == 4:
        file1, file2, file3, file4 = args
        list1 = [i.strip("\n") for i in open(file1).readlines()]
        list2 = [i.strip("\n") for i in open(file2).readlines()]
        list3 = [i.strip("\n") for i in open(file3).readlines()]
        list4 = [i.strip("\n") for i in open(file4).readlines()]
        labels = venn.get_labels([list1, list2, list3, list4], fill=['number'])
        fig, ax = venn.venn4(labels, names=[file1, file2, file3, file4])
        fig.savefig(outFile)
        # fig.show()
    if len(args) == 5:
        file1, file2, file3, file4, file5 = args
        list1 = [i.strip("\n") for i in open(file1).readlines()]
        list2 = [i.strip("\n") for i in open(file2).readlines()]
        list3 = [i.strip("\n") for i in open(file3).readlines()]
        list4 = [i.strip("\n") for i in open(file4).readlines()]
        list5 = [i.strip("\n") for i in open(file5).readlines()]
        labels = venn.get_labels([list1, list2, list3, list4, list5], fill=['number'])
        fig, ax = venn.venn5(labels, names=[file1, file2, file3, file4, file5])
        fig.savefig(outFile)
        # fig.show()
    if len(args) == 6:
        file1, file2, file3, file4, file5, file6 = args
        list1 = [i.strip("\n") for i in open(file1).readlines()]
        list2 = [i.strip("\n") for i in open(file2).readlines()]
        list3 = [i.strip("\n") for i in open(file3).readlines()]
        list4 = [i.strip("\n") for i in open(file4).readlines()]
        list5 = [i.strip("\n") for i in open(file5).readlines()]
        list6 = [i.strip("\n") for i in open(file6).readlines()]
        labels = venn.get_labels([list1, list2, list3, list4, list5, list6], fill=['number'])
        fig, ax = venn.venn6(labels, names=[file1, file2, file3, file4, file5, file6])
        fig.savefig(outFile)
        # fig.show()
    if len(args) > 6:
        print "The number of files you input larger than 6, I can't figure out the venn plot!"
#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: isoformIdentityByJunc.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-01-15 21:33:29
Last modified: 2021-01-15 21:33:29
'''

# from commonObjs import *
import argparse, sys
from itertools import islice
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool
import random
import cPickle as pickle

parser = argparse.ArgumentParser()
parser.add_argument("-r", dest="refBed", type=str, default=False, required=True,
                    help="The reference bed used to be compared")
parser.add_argument("-q", dest="queryBed", type=str, default=False, required=True,
                    help="The query bed used to compare")
parser.add_argument("-w", dest="offset", type=int, default=3,
                    help="The offset used to tolerant the bias in splice site")
parser.add_argument("-t", dest="threads", type=int, default=10,
                    help="The threads number to reduce running time")
args = parser.parse_args()

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

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
            self.score = str(self.score)
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
        return readsDict

def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]

def getConsensusIntronN(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
    '''
    get consensus intron number, i.e. the intron info of reads is identical to the reference
    :return: the count of the consensus introns
    '''
    intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
    intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
    j, consensusN = 0, 0
    junc1, junc2 = [], []
    for i in range(len(intronStarts1)):
        for k in range(j, len(intronStarts2)):
            if intronStarts2[k] > intronStarts1[i]: break
            if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronStarts1[i] + offset and \
                    intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
                consensusN += 1
                j += 1
                junc1.append((intronStarts1[i], intronEnds1[i]))
                junc2.append((intronStarts2[k], intronEnds2[k]))
    return consensusN, junc1, junc2

def getBlockLength(blockList):
    return sum(map(lambda x: int(x[1]) - int(x[0]), blockList))

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

def juncMatchByFiveEnd(refJunc, queryJunc, queryBedExons, strand, overlapRatio2query):
    matchedIndex = []
    if strand == "+":
        for i in range(len(refJunc)):
            for j in range(len(queryJunc)):
                if refJunc[i] == queryJunc[j]:
                    matchedIndex.append(i)
    else:
        for i in reversed(range(len(refJunc))):
            for j in reversed(range(len(queryJunc))):
                if refJunc[i] == queryJunc[j]:
                    matchedIndex.append(i)
    if overlapRatio2query >= 0.8:
        if max(matchedIndex) - min(matchedIndex) + 1 != len(queryJunc):
            return "exon_alt_in_middle"
        else:
            refJuncOverlapQueryExon1 = getBlockLength(getOverlapOfTuple(refJunc[0: min(matchedIndex)], queryBedExons))
            refJuncOverlapQueryExon2 = getBlockLength(getOverlapOfTuple(refJunc[max(matchedIndex) + 1: len(refJunc)], queryBedExons))
            if refJuncOverlapQueryExon1 > 10 or refJuncOverlapQueryExon2 > 10:
                return "incorrect"
            else:
                return "partial_truncated"
    else:
        return "incorrect"

def chunksDict(data, chunkNum=10):
    dataLen = len(data)
    chunkSize = dataLen/chunkNum
    it = iter(data)
    for i in xrange(0, len(data), chunkSize):
        yield {k: data[k] for k in islice(it, chunkSize)}

def chunksList(lst, chunkNum=10):
    """Yield successive n-sized chunks from lst."""
    dataLen = len(lst)
    chunkSize = dataLen/chunkNum
    for i in range(0, dataLen, chunkSize):
        yield lst[i:i + chunkSize]

def process1(queryPickle, refPickle, offset=10):
    resultDict = {}
    queryBedDict = pickle.loads(queryPickle)
    refBedDict = pickle.loads(refPickle)
    for chrom in queryBedDict:
        for n in queryBedDict[chrom]:
            queryBed = queryBedDict[chrom][n]
            queryBedExons = queryBed.parse_exon()
            queryBedIntrons = queryBed.parse_intron()
            if chrom not in refBedDict: continue
            for m in refBedDict[chrom]:
                refBed = refBedDict[chrom][m]
                refBedExons = refBed.parse_exon()
                if not getOverlapOfTuple(queryBedExons, refBedExons): continue
                overlapLength = getBlockLength(getOverlapOfTuple(queryBedExons, refBedExons))
                refExonLength = getBlockLength(refBedExons)
                queryExonLength = getBlockLength(queryBedExons)
                overlapRatio = round(float(overlapLength) / refExonLength, 3)
                overlapRatio2query = round(float(overlapLength) / queryExonLength, 3)
                if overlapRatio < 0.5 and overlapRatio2query < 0.5:
                    continue
                if len(queryBedExons) == 1:
                    if len(refBedExons) != 1:
                        if queryBed.name not in resultDict:
                            resultDict[queryBed.name] = ["monoExon_truncated", 0, refBed.name, overlapRatio,
                                                         overlapRatio2query]
                        else:
                            if resultDict[queryBed.name][0] == "monoExon":
                                continue
                            if overlapRatio * overlapRatio2query * overlapRatio2query > resultDict[queryBed.name][3] * resultDict[queryBed.name][4] * resultDict[queryBed.name][4]:
                                resultDict[queryBed.name] = ["monoExon_truncated", 1, refBed.name, overlapRatio,
                                                             overlapRatio2query]
                    else:
                        if queryBed.name not in resultDict:
                            resultDict[queryBed.name] = ["monoExon", 0, refBed.name, overlapRatio, overlapRatio2query]
                        else:
                            if overlapRatio * overlapRatio2query * overlapRatio2query > resultDict[queryBed.name][3] * resultDict[queryBed.name][4] * resultDict[queryBed.name][4]:
                                resultDict[queryBed.name] = ["monoExon", 1, refBed.name, overlapRatio,
                                                             overlapRatio2query]
                    continue

                consensusIntronN, consensusIntron1, consensusIntron2 = getConsensusIntronN(refBed.exonStarts,
                                                                                           refBed.exonEnds,
                                                                                           queryBed.exonStarts,
                                                                                           queryBed.exonEnds,
                                                                                           offset=offset)
                refBedIntrons = refBed.parse_intron()
                if consensusIntronN == len(refBedIntrons) and consensusIntronN == len(queryBedIntrons):

                    if overlapRatio2query >= 0.7 or overlapRatio >= 0.7:
                        if queryBed.name not in resultDict:
                            resultDict[queryBed.name] = ["complete_identity", consensusIntronN, refBed.name,
                                                         overlapRatio, overlapRatio2query]
                        else:
                            if overlapRatio > resultDict[queryBed.name][3] or overlapRatio2query > \
                                    resultDict[queryBed.name][4] or resultDict[queryBed.name][0] != "complete_identity":
                                resultDict[queryBed.name] = ["complete_identity", consensusIntronN, refBed.name,
                                                             overlapRatio, overlapRatio2query]
                    else:
                        resultDict[queryBed.name] = ["incorrect", consensusIntronN, "none5;" + refBed.name,
                                                     overlapRatio, overlapRatio2query]
                else:
                    if queryBed.name in resultDict and resultDict[queryBed.name][0] == "complete_identity":
                        continue
                    if consensusIntronN == len(queryBedIntrons):
                        annotation = juncMatchByFiveEnd(refBedIntrons, consensusIntron1, queryBedExons, queryBed.strand, overlapRatio2query)
                        if queryBed.name not in resultDict:
                            resultDict[queryBed.name] = [annotation, consensusIntronN, refBed.name, overlapRatio,
                                                         overlapRatio2query]
                        else:
                            if resultDict[queryBed.name][0] == "incorrect":
                                if consensusIntronN > resultDict[queryBed.name][1]:
                                    resultDict[queryBed.name] = [annotation, consensusIntronN, refBed.name, overlapRatio,
                                                                 overlapRatio2query]
                                elif overlapRatio * overlapRatio2query * overlapRatio2query > resultDict[queryBed.name][3] * resultDict[queryBed.name][4] * resultDict[queryBed.name][4]:
                                    resultDict[queryBed.name] = [annotation, consensusIntronN, refBed.name,
                                                                 overlapRatio, overlapRatio2query]
                                # resultDict[queryBed.name] = ["partial_identity", consensusIntronN, refBed.name]
                            else:
                                if annotation != "incorrect":
                                    if consensusIntronN > resultDict[queryBed.name][1]:
                                        resultDict[queryBed.name] = [annotation, consensusIntronN, refBed.name,
                                                                     overlapRatio, overlapRatio2query]
                                    elif overlapRatio * overlapRatio2query * overlapRatio2query > resultDict[queryBed.name][3] * resultDict[queryBed.name][4] * resultDict[queryBed.name][4]:
                                        resultDict[queryBed.name] = [annotation, consensusIntronN, refBed.name,
                                                                     overlapRatio, overlapRatio2query]
                    else:
                        if queryBed.name not in resultDict:
                            resultDict[queryBed.name] = ["incorrect", consensusIntronN, "none1;" + refBed.name,
                                                         overlapRatio, overlapRatio2query]
                        else:
                            if resultDict[queryBed.name][0] != "incorrect":
                                continue
                            if consensusIntronN > resultDict[queryBed.name][1]:
                                resultDict[queryBed.name] = ["incorrect", consensusIntronN, "none2;" + refBed.name,
                                                             overlapRatio, overlapRatio2query]
                            elif overlapRatio * overlapRatio2query * overlapRatio2query > resultDict[queryBed.name][3] * resultDict[queryBed.name][4] * resultDict[queryBed.name][4]:
                                resultDict[queryBed.name] = ["incorrect", consensusIntronN, "none3;" + refBed.name,
                                                             overlapRatio, overlapRatio2query]
            if len(queryBedExons) == 1:
                if queryBed.name not in resultDict:
                    resultDict[queryBed.name] = ["incorrect", 0, "none3", 0, 0]
            else:
                if queryBed.name not in resultDict:
                    resultDict[queryBed.name] = ["incorrect", len(queryBedIntrons), "none4", 0, 0]
    return resultDict

def restoreByChrom(myBed):
    myDict = {}
    myList = []
    for n in myBed.reads:
        if myBed.reads[n].chrom not in myDict:
            myDict[myBed.reads[n].chrom] = {n: myBed.reads[n]}
        else:
            myDict[myBed.reads[n].chrom].update({n: myBed.reads[n]})
        myList.append([myBed.reads[n].chrom, {n: myBed.reads[n]}])
    return myDict, myList

def main():
    refBed = BedFile(args.refBed, type="bed12")
    queryBed = BedFile(args.queryBed, type="bed12")
    refBedDict, refBedList = restoreByChrom(refBed)
    queryBedDict, queryBedList = restoreByChrom(queryBed)
    random.shuffle(queryBedList)
    queryList = []
    for q in chunksList(queryBedList, chunkNum=args.threads):
        tmpDict = {}
        for i in q:
            if i[0] not in tmpDict:
                tmpDict.update({i[0]: i[1]})
            else:
                tmpDict[i[0]].update(i[1])
        queryList.append(tmpDict)

    pool = Pool(processes=args.threads)
    resultList = []
    refPickle = pickle.dumps(refBedDict)
    for tmpDict in queryList:
        queryPickle = pickle.dumps(tmpDict)
        tmpRes = pool.apply_async(process1, (queryPickle, refPickle, args.offset,))
        resultList.append(tmpRes)

    for i in resultList:
        i.wait()
    finalRes = {}
    for i in resultList:
        if i.ready():
            if i.successful():
                tmpRes = i.get()
                for j in tmpRes:
                    print j + "\t" + "\t".join(map(str, tmpRes[j]))
                finalRes.update(tmpRes)
    for x in ["complete_identity", "monoExon", "partial_truncated", "monoExon_truncated", "incorrect", "exon_alt_in_middle"]:
        count = 0
        for j in finalRes:
            if finalRes[j][0] == x:
                count += 1
        print >> sys.stderr, "\t".join(map(str, [x, count, "collapsed"]))
    supportList = set()
    for j in finalRes:
        if finalRes[j][0] in ["complete_identity", "monoExon", "partial_truncated", "monoExon_truncated", "exon_alt_in_middle"]:
            supportList.add(finalRes[j][2])
    print >> sys.stderr, "\t".join(map(str, ["covered", len(supportList), "covered"]))

if __name__ == '__main__':
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: juncOverlapBetweenFiles.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-01-27 13:01:32
Last modified: 2021-01-27 13:01:32
'''

import sys
import numpy as np
import pandas as pd
from copy import copy

def getIntrons(exonStarts, exonEnds):
    return exonEnds[:-1], exonStarts[1:]

# def getConsensusIntronN(exonStarts1, exonEnds1, exonStarts2, exonEnds2, offset):
#     '''
#     get consensus intron number, i.e. the intron info of reads is identical to the reference
#     :return: the count of the consensus introns
#     '''
#     intronStarts1, intronEnds1 = getIntrons(exonStarts1, exonEnds1)
#     intronStarts2, intronEnds2 = getIntrons(exonStarts2, exonEnds2)
#     j, consensusN = 0, 0
#     junc1, junc2 = [], []
#     for i in range(len(intronStarts1)):
#         for k in range(j, len(intronStarts2)):
#             if intronStarts2[k] > intronStarts1[i]: break
#             if intronStarts1[i] - offset <= intronStarts2[k] and intronStarts2[k] <= intronStarts1[i] + offset and \
#                     intronEnds1[i] - offset <= intronEnds2[k] and intronEnds2[k] <= intronEnds1[i] + offset:
#                 consensusN += 1
#                 j += 1
#                 junc1.append((intronStarts1[i], intronEnds1[i]))
#                 junc2.append((intronStarts2[k], intronEnds2[k]))
#     return consensusN, junc1, junc2

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

def compareJuncIdentity(juncList, junc2, offset):
    for junc1 in juncList:
        consensusN = 0
        if junc1[0][0] > junc2[-1][-1]: break
        if len(junc1) != len(junc2): continue
        if not getOverlapOfTuple(junc1, junc2): continue
        for i in range(len(junc1)):
            if junc1[i][0] - offset <= junc2[i][0] and junc2[i][0] <= junc1[i][0] + offset and \
                    junc1[i][1] - offset <= junc2[i][1] and junc2[i][1] <= junc1[i][1] + offset:
                consensusN += 1
        if consensusN == len(junc1) and consensusN == len(junc2):
            return True
    else:
        return False


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
        return tuple(exons)

    def parse_intron(self):
        "return a list of intron pos [(st, ed], (st, ed], ... ]"
        introns = []
        for i in range(self.blockCount - 1):
            st = self.exonEnds[i]
            ed = self.exonStarts[i + 1]
            introns.append((st, ed))
        return tuple(introns)

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

def getJuncList(bedFile):
    reads = BedFile(bedFile, type="bed12").reads
    juncDict = {}
    junc2id = {}
    for i in reads:
        junc = reads[i].parse_intron()
        chrom = reads[i].chrom
        if chrom not in juncDict:
            juncDict[chrom] = []
            junc2id[chrom] = {}
        if len(junc) >= 1:
            juncDict[chrom].append(junc)
            if junc not in junc2id[chrom]:
                junc2id[chrom][junc] = [i]
            else:
                junc2id[chrom][junc].append(i)
    for c in juncDict:
        # juncDict[c] = sorted(list(np.unique(juncDict[c])))
        juncDict[c] = tuple(set(juncDict[c]))
    return juncDict, junc2id

def getAllJuncList(bed2junc, offset):
    allJuncDict = {}
    for i in bed2junc:
        for c in bed2junc[i]:
            if c not in allJuncDict:
                allJuncList = []
                allJuncDict[c] = {}

                for j in bed2junc[i][c]:
                    if len(allJuncList) == 0:
                        allJuncList.append(j)
                    else:
                        if not compareJuncIdentity(allJuncList, j, offset=offset):
                            allJuncList.append(j)
                # allJuncList = sorted(allJuncList)
                for k in range(len(allJuncList)):
                    allJuncDict[c][c + "_junc_"+str(k)] = allJuncList[k]
            else:
                allJuncList = sorted(allJuncDict[c].values())
                newAllJuncList = copy(allJuncList)
                tmpList = []
                for j in bed2junc[i][c]:
                    if len(newAllJuncList) == 0:
                        newAllJuncList.append(j)
                    else:
                        if not compareJuncIdentity(newAllJuncList, j, offset=offset):
                            newAllJuncList.append(j)
                            tmpList.append(j)
                            newAllJuncList = sorted(list(set(newAllJuncList)))
                for k in range(len(tmpList)):
                    allJuncDict[c][c + "_junc_"+str(k+len(allJuncList))] = tmpList[k]
    return allJuncDict

def compare(allJuncDict, bed2junc, junc2idDict, offset):
    juncDistribution = {}

    for c in allJuncDict:
        for i in bed2junc:
            if c not in bed2junc[i]:
                continue
            for k in bed2junc[i][c]:
                for j in allJuncDict[c]:
                    if compareJuncIdentity([allJuncDict[c][j]], k, offset=offset):
                        if j not in juncDistribution:
                            juncDistribution[j] = dict.fromkeys(bed2junc.keys(), 0)
                            juncDistribution[j][i] += 1
                        else:
                            juncDistribution[j][i] += 1
                        print >>sys.stderr, "\t".join([j, i]) + "\t" + "\t".join(junc2idDict[i][c][k])
    return juncDistribution

def main():
    bedFiles = sys.argv[1:]
    offset = 10
    bed2junc = {}
    junc2idDict = {}
    for i in bedFiles:
        juncDict, junc2id = getJuncList(i)
        bed2junc[i] = juncDict
        junc2idDict[i] = junc2id
    allJuncDict = getAllJuncList(bed2junc, offset=offset)
    juncDistribution = compare(allJuncDict, bed2junc, junc2idDict, offset=offset)
    print "junc_id" + "\t" + "\t".join(bedFiles)
    for j in juncDistribution:
        print j + "\t" + "\t".join(map(str, [juncDistribution[j][k] for k in bedFiles]))


if __name__ == '__main__':
    main()
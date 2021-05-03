#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: filterSpecificExonInEs.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-03-31 17:05:24
Last modified: 2020-03-31 17:05:24
'''

import argparse

def processEsFile(esFile):
    myDict = {}
    recordToExon = {}
    with open(esFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            chrom = records[0]
            esInfoList = records[3].split("@")
            esExons = esInfoList[1].split(";")
            recordToExon["\t".join(records)] = []
            for i in esExons:
                exon = "{}:{}".format(chrom, i)
                if exon not in myDict:
                    myDict[exon] = [records]
                else:
                    myDict[exon].append(records)
                recordToExon["\t".join(records)].append(exon)
    return myDict, recordToExon

def processBedFile(bedFile):
    myDict = {}
    with open(bedFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            chrom, chromStart = records[0], int(records[1])
            blockSizes = [int(i) for i in records[10].split(",")]
            exonStarts = [chromStart + int(i) for i in records[11].split(",")]
            exonEnds = [exonStarts[i] + blockSizes[i] for i in range(len(blockSizes))]
            for i in range(len(exonStarts)):
                exon = "{}:{}-{}".format(chrom, exonStarts[i], exonEnds[i])
                if exon not in myDict:
                    myDict[exon] = ""
    return myDict

def filter(fref=None, source=None, target=None, source_all=None, target_all=None, mode=None, outFile=None):
    out = open(outFile, "w")
    refExonDict = {}
    with open(fref) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            transName, chrom = records[0], records[1]
            exonStarts = [int(i) for i in records[8].strip(',').split(',')]
            exonEnds = [int(i) for i in records[9].strip(',').split(',')]
            for i in range(len(exonStarts)):
                exon = "{}:{}-{}".format(chrom, exonStarts[i], exonEnds[i])
                if exon not in refExonDict:
                    refExonDict[exon] = [transName]
                else:
                    refExonDict[exon].append(transName)
    sourceDict, sourceRecordToExonDict = processEsFile(source)
    targetDict, targetRecordToExonDict = processEsFile(target)
    sourceAllDict = processBedFile(source_all)
    targetAllDict = processBedFile(target_all)
    refExons = set(refExonDict.keys())
    sourceExons = set(sourceDict.keys())
    sourceAllExons = set(sourceAllDict.keys())
    targetExons = set(targetDict.keys())
    targetAllExons = set(targetAllDict.keys())
    if mode == "e1":
        filterExons = targetExons - (sourceAllExons | refExons)
        for i in filterExons:
            for j in targetDict[i]:
                if len(set(targetRecordToExonDict["\t".join(j)]) & (sourceAllExons | refExons)) == 0:
                    print >> out, "\t".join(j)
    elif mode == "e2":
        filterExons = sourceExons - targetAllExons
        for i in filterExons:
            if i in sourceDict:
                for j in sourceDict[i]:
                    if len(set(sourceRecordToExonDict["\t".join(j)]) & targetAllExons) == 0:
                        print >> out, "\t".join(j)
    else:
        filterExons = targetExons & (refExons | sourceExons)
        for i in filterExons:
            for j in targetDict[i]:
                print >>out, "\t".join(j)
    out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-ref", type=str,
                        help="Reference gene structure in GPE format")
    parser.add_argument("-s", type=str,
                        help="source file")
    parser.add_argument("-sa", type=str,
                        help="Target file")
    parser.add_argument("-t", type=str,
                        help="Target file")
    parser.add_argument("-ta", type=str,
                        help="Target file")
    parser.add_argument("-m", type=str,
                        help="Mode of filtering")
    parser.add_argument("-o", type=str,
                        help="Output file")
    args = parser.parse_args()
    filter(fref=args.ref, source=args.s, target=args.t, source_all=args.sa, target_all=args.ta, mode=args.m, outFile=args.o)
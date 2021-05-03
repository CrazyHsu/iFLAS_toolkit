#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: PAandAS.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-02-28 17:58:59
Last modified: 2020-02-28 17:58:59
'''
import subprocess, argparse
from commonObjs import *

def relationshipBetweenPAandAS():
    cmd = '''sed '1d' ../isoseq3/flnc.report.csv | tr ',' '\t' | cut -f 1,5 > flnc.pa.len.txt'''
    subprocess.call(cmd, shell=True)

def clusterByFiveEnd(reads, strand):
    if strand == "+":
        firstFiveEndExon = list(reads[0].exons[0])
        for i in range(1, len(reads)):
            if isOverlap(firstFiveEndExon, reads[i].exons[0]):
                firstFiveEndExon[0] = min(firstFiveEndExon[0], reads[i].exons[0][0])
                firstFiveEndExon[1] = max(firstFiveEndExon[1], reads[i].exons[0][1])
            else:
                return False
        return firstFiveEndExon
    elif strand == "-":
        firstFiveEndExon = list(reads[0].exons[-1])
        for i in range(1, len(reads)):
            if isOverlap(firstFiveEndExon, reads[i].exons[-1]):
                firstFiveEndExon[0] = min(firstFiveEndExon[0], reads[i].exons[-1][0])
                firstFiveEndExon[1] = max(firstFiveEndExon[1], reads[i].exons[-1][1])
            else:
                return False
        return firstFiveEndExon

def asProcess(asFile=None, allReadsFile=None, asType=None, out=False):
    allReads = BedFile(allReadsFile, type="bed12+").reads
    outHandle = open(out, "w")

    with open(asFile) as f:
        for line in f.readlines():
            if asType == "SE":
                records = line.strip("\n").split("\t")
                strand = records[5]
                inclusionReadStr, exclusionReadStr = records[15], records[17]
                inclusionReadList = inclusionReadStr.split(",")
                exclusionReadList = exclusionReadStr.split(",")
                inclusionReads = [allReads[i] for i in inclusionReadList]
                exclusionReads = [allReads[i] for i in exclusionReadList]
                if len(inclusionReads) == 1 and len(exclusionReads) == 1: continue
                inclusionFirstFiveEndExon = clusterByFiveEnd(inclusionReads, strand)
                exclusionFirstFiveEndExon = clusterByFiveEnd(exclusionReads, strand)
                if inclusionFirstFiveEndExon and exclusionFirstFiveEndExon and not isOverlap(inclusionFirstFiveEndExon, exclusionFirstFiveEndExon):
                    print >>outHandle, "\t".join(records)
            elif asType == "PA":
                pass
            else:
                pass
    outHandle.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-all", type=str,
                        help="All transcript from iFLAS pipeline")
    parser.add_argument("-a", type=str,
                        help="AS file from iFLAS pipeline")
    parser.add_argument("-t", type=str,
                        help="AS type for the AS file")
    parser.add_argument("-o", type=str,
                        help="Output file")
    args = parser.parse_args()
    asProcess(asFile=args.a, allReadsFile=args.all, asType=args.t, out=args.o)
    # asProcess(asFile=args.as, all=args.all, asType=args.t, out=args.o)
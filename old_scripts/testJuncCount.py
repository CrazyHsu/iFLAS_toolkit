#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: testJuncCount.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-02-25 17:40:31
Last modified: 2020-02-25 17:40:31
'''

import subprocess, psutil, os, itertools, pysam, re
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

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

def removeFiles(myDir, fileList):
    for f in fileList:
        os.remove(os.path.join(myDir, f.strip("\n")))

def parallelJuncCount(bedFile, bamFile, outFile, filterCount=None):
    out = open(outFile, "w")
    bam = pysam.AlignmentFile(bamFile, "rb", check_sq=False)
    with open(bedFile) as f:
        count = 1
        for line in f:
            lineList = line.strip("\n").split("\t")
            chrom, start, end, strand = lineList[0], int(lineList[1]), int(lineList[2]), lineList[3]
            reads = bam.fetch(chrom, start, end)
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
            if len(readListWithJunc) == 0: continue
            readsMinPos = min([i.reference_start for i in readListWithJunc])
            readsMaxPos = max([i.reference_end for i in readListWithJunc])
            junc = TopHatJunction()
            junc.chrom, junc.start, junc.end, junc.strand, junc.juncSupport = chrom, readsMinPos, readsMaxPos, strand, len(
                readListWithJunc)
            if "chr" in chrom.lower():
                junc.name = "JUNC_{}_".format(chrom.upper()) + "{0:08d}".format(count)
            else:
                junc.name = "JUNC_{}_".format(chrom) + "{0:08d}".format(count)
            junc.juncLeft, junc.juncRight = start, end
            print >> out, junc
            count += 1
    out.close()
    bam.close()

def getJuncCountFromNGSassembly(group, assembledGTF, bamFile, outFile=None, threads=None, filterCount=None, generateTrack=False):
    prevDir = os.getcwd()
    os.chdir("/data/CrazyHsu_data/Projects/isoseq/FLAS_alpha/test_AGPv4.36_20200221/project3/"+group+"/RNA-seq/reassembly")
    cmd = "hisat2_extract_splice_sites.py {} > tmp.ss".format(assembledGTF)
    subprocess.call(cmd, shell=True)
    cmd = '''
            rm chr* && awk '{close(f);f=$1}{print > "chr"f".bed"}' tmp.ss
        '''
    subprocess.call(cmd, shell=True)
    processNum = threads if threads else psutil.cpu_count() * 3 / 4
    try:
        pool = Pool(processes=processNum)
        print "into parallel"
        multiResults = []
        bedFiles = os.popen("ls chr*")
        juncBedList = []
        chrBedList = []
        for singleBed in bedFiles:
            juncBed = singleBed.strip("\n").split(".")[0] + ".junctions.bed"
            singleRunRes = pool.apply_async(parallelJuncCount, (singleBed.strip("\n"), bamFile, juncBed, filterCount))
            multiResults.append(singleRunRes)
            juncBedList.append(juncBed)
            chrBedList.append(singleBed.strip("\n"))
        for j in multiResults:
            j.wait()
        with open(outFile, "w") as out:
            if generateTrack:
                out.write('track name=junctions description="iFLAS TopHat2-like junctions"\n')
            for line in itertools.chain.from_iterable(itertools.imap(open, juncBedList)):
                out.write(line)
        removeFiles(os.getcwd(), juncBedList)
        removeFiles(os.getcwd(), chrBedList)
        removeFiles(os.getcwd(), ["tmp.bam", "tmp.bam.bai", "tmp.ss"])
        os.chdir(prevDir)
    except Exception as e:
        print e

pool = MyPool(2)
# pool = MyPool(processes=2)
multiRes = []
for i in ["2", "4"]:
    # prevDir = os.getcwd()
    # os.chdir("/data/CrazyHsu_data/Projects/isoseq/FLAS_alpha/test_AGPv4.36_20200221/project3/"+i+"/RNA-seq/reassembly")
    singleRunRes = pool.apply_async(getJuncCountFromNGSassembly, (i, "stringtie_merged.gtf", "tmp.bam", "junctions.bed12", 48))
    # getJuncCountFromNGSassembly("stringtie_merged.gtf", "tmp.bam", outFile="junctions.bed12", threads=96)
    multiRes.append(singleRunRes)
    # os.chdir(prevDir)

for j in multiRes:
    j.wait()

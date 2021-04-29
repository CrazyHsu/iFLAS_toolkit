#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: merge_splited_collapse.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-18 14:42:26
Last modified: 2021-04-18 14:42:26
'''

'''
1       PacBio  transcript      50871   52827   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    50871   51217   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    51371   51436   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    51886   52004   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    52137   52294   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    52391   52546   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
1       PacBio  exon    52668   52827   .       -       .       gene_id "PB.2"; transcript_id "PB.2.5";
'''

import os, sys, shutil
from Bio import SeqIO

class GTFLine(object):
    """
    Gene Transfer Format.
    Attributes: seqname (chromsome number or scaffold), source, feature, start,
                end, score, strand, frame, attributes (key-value pairs).
    Methods:
    """

    def __init__(self, line=""):
        if line:
            self.record = line.strip().split("\t")
            self.seqname, self.source, self.feature = self.record[:3]
            self.start = int(self.record[3])
            self.end = int(self.record[4])
            self.score = int(float(self.record[5])) if self.record[5] != "." else None
            self.strand = self.record[6] if self.record[6] != "." else None
            assert self.strand in ("+", "-", None)
            self.frame = int(self.record[7]) if self.record[7] != "." else None
            assert self.frame in (0, 1, 2, None)
            attri = self.record[8].strip().strip(";").split("; ")
            self.attributes = {}
            for pair in attri:
                key, value = pair.split(" ", 1)
                value = value.strip("\"")
                self.attributes[key] = value
        else:
            self.empty()

    def empty(self):
        self.record = []
        self.seqname, self.source, self.feature = "", "", ""
        self.start = self.end = 0
        self.score = 0
        self.strand = ''
        self.frame = '.'
        self.attributes = {}

    def __repr__(self):
        attriLst = []
        for key, value in self.attributes.items():
            attriLst.append("%s \"%s\"" % (str(key), str(value)))
        attriStr = "; ".join(attriLst)
        return "\t".join([self.seqname, self.source, self.feature,
                          repr(self.start), repr(self.end),
                          repr(self.score) if self.score else ".",
                          self.strand if self.strand else ".",
                          repr(self.frame) if self.frame else ".",
                          attriStr])

def changeGffName(inputGff, tmpFile, addNum):
    tmpOut = open(tmpFile, "w")
    with open(inputGff) as f:
        for line in f.readlines():
            gtfLine = GTFLine(line)
            geneInfoList = gtfLine.attributes["gene_id"].split(".")
            transInfoList = gtfLine.attributes["transcript_id"].split(".")
            gtfLine.attributes["gene_id"] = "PB.{}".format(int(geneInfoList[1]) + addNum)
            gtfLine.attributes["transcript_id"] = "PB.{}.{}".format(int(transInfoList[1]) + addNum, transInfoList[2])
            print >>tmpOut, gtfLine
    tmpOut.close()

def changeGroupName(inputGroupFile, tmpFile, addNum):
    tmpOut = open(tmpFile, "w")
    with open(inputGroupFile) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            transInfoList = infoList[0].split(".")
            newTrans = "PB.{}.{}".format(int(transInfoList[1]) + addNum, transInfoList[2])
            print >>tmpOut, "{}\t{}".format(newTrans, infoList[1])
    tmpOut.close()

def changeFaName(inputFaFile, tmpFile, addNum):
    tmpOut = open(tmpFile, "w")
    for record in SeqIO.parse(inputFaFile, "fasta"):
        recordInfo = record.id.split("|")
        transInfoList = recordInfo[0].split(".")
        newTrans = "PB.{}.{}".format(int(transInfoList[1]) + addNum, transInfoList[2])
        newName = "{}|{}".format(newTrans, "|".join(recordInfo[1:]))
        record.id = newName
        SeqIO.write(record, tmpOut, "fasta")
    tmpOut.close()


targetDir = sys.argv[1]
filePrefix = sys.argv[2]
num = int(sys.argv[3])

def main():
    goodGff = os.path.join(targetDir, filePrefix + ".collapsed.good.gff")
    goodFuzzyGff = os.path.join(targetDir, filePrefix + ".collapsed.good.gff.unfuzzy")
    badGff = os.path.join(targetDir, filePrefix + ".collapsed.bad.gff")
    groupFile = os.path.join(targetDir, filePrefix + ".collapsed.group.txt")
    groupFuzzyFile = os.path.join(targetDir, filePrefix + ".collapsed.group.txt.unfuzzy")
    faFile = os.path.join(targetDir, filePrefix + ".collapsed.rep.fa")

    tmpGoodGff = os.path.join(targetDir, "tmp.good.gff")
    tmpGoodFuzzyGff = os.path.join(targetDir, "tmp.good.fuzzy.gff")
    tmpBadGff = os.path.join(targetDir, "tmp.bad.gff")
    tmpGroupFile = os.path.join(targetDir, "tmp.group.txt")
    tmpGroupFuzzyGff = os.path.join(targetDir, "tmp.group.fuzzy.gff")
    tmpFaFile = os.path.join(targetDir, "tmp.fa")

    changeGffName(goodGff, tmpGoodGff, num)
    changeGffName(goodFuzzyGff, tmpGoodFuzzyGff, num)
    changeGffName(badGff, tmpBadGff, num)
    changeGroupName(groupFile, tmpGroupFile, num)
    changeGroupName(groupFuzzyFile, tmpGroupFuzzyGff, num)
    changeFaName(faFile, tmpFaFile, num)

    # shutil.move(tmpGoodGff, goodGff)
    # shutil.move(tmpGoodFuzzyGff, goodFuzzyGff)
    # shutil.move(tmpBadGff, badGff)
    # shutil.move(tmpGroupFile, groupFile)
    # shutil.move(tmpGroupFuzzyGff, groupFuzzyFile)
    # shutil.move(tmpFaFile, faFile)

if __name__ == '__main__':
    main()

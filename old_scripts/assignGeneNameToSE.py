#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: assignGeneNameToSE.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-01-03 16:07:53
Last modified: 2020-01-03 16:07:54
'''

import sys
from commonObjs import Position, GenePredExtLine, Bed12


class SE(Bed12):
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.geneName, self.se_info = self.name.split(":")
        self.supCount = int(self.record[12])
        self.totalCount = int(self.record[13])
        self.juncPoses = [Position(j) for j in self.record[14].split(",")]

    def __repr__(self):
        return "\t".join([Bed12.__repr__(self), repr(self.supCount),
                          repr(self.totalCount),
                          ",".join([repr(j) for j in self.juncPoses])])


def find_refs(se, refLst):
    refs = []
    for ref in refLst:
        if ref.strand != se.strand:
            continue
        else:
            for junc in se.juncPoses:
                if (junc.chromStart in ref.exonEnds or
                junc.chromEnd in ref.exonStarts):
                    refs.append(ref)
                    break
    return refs

def change_name(se, refDct, errorOut):
    "change name if necessary"
    if se.chrom not in refDct:
        return 1
    refs = find_refs(se, refDct[se.chrom])
    if refs:
        geneNames = set([ref.geneName for ref in refs])
        if len(geneNames) > 1:
            return 2
        else:
            gName = geneNames.pop()
            print >> errorOut, "%s\t%s" % (se.geneName, gName)
            se.geneName = gName
            se.name = ":".join([se.geneName, se.se_info])
            return 0
    else:
        return 1

def assignGeneNameToSE(fref, has_bin, fes, outFile, errorOutFile):
    # load gene structures
    refDct = {}
    with open(fref) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ref = GenePredExtLine(line, bincolumn=has_bin)
            if ref.chrom in refDct:
                refDct[ref.chrom].append(ref)
            else:
                refDct[ref.chrom] = [ref]
    for chrom in refDct:
        refDct[chrom].sort(key=lambda r: r.txStart)

    out = open(outFile, "w")
    errorOut = open(errorOutFile, "w")
    with open(fes) as f:
        for line in f:
            if line.startswith("#"):
                print >>out, line,
                continue
            se = SE(line)
            if se.geneName.startswith("Novel"):
                change_name(se, refDct, errorOut)
                print >>out, se
            else:
                print >>out, line,
    out.close()
    errorOut.close()

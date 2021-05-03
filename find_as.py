#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: find_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:16:33
Last modified: 2021-04-29 16:16:33
'''
from commonFuncs import *
from commonObjs import *
import pybedtools

def find_perl():
    pass

def find_python():
    pass

def find_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {}...".format(
        projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dataObj.out_dir, projectName, sampleName)
    logDir = os.path.join(baseDir, "log")

    workDir = os.path.join(baseDir, "as_events")
    resolveDir(workDir)
    isoformBed = os.path.join(baseDir, "collapse", "isoformGrouped.bed12+")
    tofuGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    if dataObj.ngsJunctions == None:
        dataObj.ngsJunctions = os.path.join(baseDir, "rna-seq", "reassembly", "junctions.bed")

    if not os.path.exists("PB"):
        os.makedirs("PB")
    transBedList = GenePredObj(refParams.ref_gpe, bincolumn=False).toBed(gene=True)
    transBedObj = pybedtools.BedTool("\n".join(transBedList), from_string=True)
    readsBedList = []
    isoform2reads = {}
    with open("tofu.collapsed.group.txt") as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    with open("isoformGrouped.bed12+") as f:
        for line in f.readlines():
            readStruc = Bed12(line)
            if len(isoform2reads[readStruc.name]) < 2: continue
            readsBedList.append("\t".join(readStruc.record[:12]))


    print getCurrentTime() + " Alternative splicing events identifying for project {} sample {}...".format(
        projectName, sampleName)

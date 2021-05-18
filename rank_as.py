#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: rank_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:19:25
Last modified: 2021-04-29 16:19:25
'''
from commonFuncs import *
from rank_as_functions import *

def rank_as(dataObj=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Start analysis independence between alternative splicing(AS) isoforms for project {} sample {}...".format(projectName, sampleName)
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    prevDir = os.getcwd()
    indepCombDir = os.path.join(baseDir, "isoIdentity")
    resolveDir(indepCombDir)
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "IR.confident.bed6+")
    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "SE.confident.bed12+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A3SS.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A5SS.confident.bed6+")
    paFile = None
    isoformFile = os.path.join(baseDir, "collapse", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    # asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile)
    novelIsoformFile = os.path.join(os.path.join(baseDir, "collapse", "isoformGrouped.novel.bed12+"))
    scoreAsIsoform(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, collapsedTrans2reads, novelIsoformFile)
    os.chdir(prevDir)
    print getCurrentTime() + " End analysis independence between alternative splicing(AS) isoforms for project {} sample {}!".format(projectName, sampleName)
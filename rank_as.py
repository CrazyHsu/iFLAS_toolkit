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
    indepCombDir = os.path.join(baseDir, "indepComb")
    resolveDir(indepCombDir)
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "IR.bed6+")
    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "SE.bed12+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "A3SS.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "A5SS.bed6+")
    paFile = None
    isoformFile = os.path.join(baseDir, "as_events", "isoformGrouped.bed12+")
    isoform2readsFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile)
    os.chdir(prevDir)
    print getCurrentTime() + " End analysis independence between alternative splicing(AS) isoforms for project {} sample {}!".format(projectName, sampleName)
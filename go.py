#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: go.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:45
Last modified: 2021-04-29 16:20:45
'''
from commonFuncs import *
from rpy2 import robjects

def go(args):
    print getCurrentTime() + " Perform GO enrichment for target genes..."
    targetGeneFile = args.targetGeneFile
    gene2goFile = args.gene2goFile
    sampleName = args.sampleName
    # scriptDir = os.path.dirname(os.path.realpath(__file__))
    if not validateFile(targetGeneFile):
        validateIndex = []
        if "," in targetGeneFile:
            fileList = targetGeneFile.split(",")
            sampleList = sampleName.split(",")
            for x in range(len(fileList)):
                validateIndex.append(x)
            if len(validateIndex) != 0:
                validateTargetGeneFiles = [fileList[x] for x in validateIndex]
                validateSampleNames = [sampleList[x] for x in validateIndex]
            else:
                raise Exception("The path of gene file you input seems not correct, please check it!")
        else:
            raise Exception("Please input the correct gene file you want to perform GO enrichment!")
    else:
        validateTargetGeneFiles = targetGeneFile
        validateSampleNames = sampleName
    from plotRscriptStrs import plotTargetGenesGoEnrichmentStr
    robjects.r(plotTargetGenesGoEnrichmentStr)
    robjects.r.plotTargetGenesGoEnrichment(",".join(validateTargetGeneFiles), ",".join(validateSampleNames), gene2goFile, args.out)
    # cmd = "Rscript {}/plotGoEnrich.R -g={} -s={} -bg={} -o={}".format(scriptDir, ",".join(validateTargetGeneFiles), ",".join(validateSampleNames), gene2goFile, args.out)
    # subprocess.call(cmd, shell=True)
    print getCurrentTime() + " Perform GO enrichment for target genes done!"

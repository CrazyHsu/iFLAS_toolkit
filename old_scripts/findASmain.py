#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: findASmain.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-09-08 13:55:17
Last modified: 2019-09-08 13:55:18
'''
import argparse
from commonObjs import Sample, ReadLineStruc, GenePredObj
from findAS import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-g", type=str, default="raw.gpe",
                    help="Reference gene structure in GPE format")
parser.add_argument("-a", type=str, default="IR",
                    help="Ths type of AS")
parser.add_argument("-s", type=str, default="deSingleExonRead.bed12+",
                    help="The config file which contain the sample name and the responding path of deSingleExonRead.bed12+ ")
parser.add_argument("-o", type=str, default=False,
                    help="The output file of the result. Default: print to the stdout")
args = parser.parse_args()

asType, gpeFile, sourceFile, outFile = args.a, args.g, args.s, args.o
def findASmain(asType="IR", gpeFile=None, sourceFile=None, outFile=False):
    gpeObj = GenePredObj(gpeFile, False)
    gpeDict = gpeObj.genePredDict
    sampleObj = Sample()
    sample = "temp"
    with open(sourceFile) as f:
        for line in f:
            readStruc = ReadLineStruc(line)
            sampleObj.update(sample, sourceFile, readStruc)
    geneDict, novelDict = {}, {}
    for g in sampleObj.sample2Gene[sample]:
        getAS(gpeDict, sampleObj.sample2Gene[sample][g], geneDict, novelDict, offset=0)
        # if sampleObj.sample2Gene[sample][g].geneName == testGene:
        #     pass
    out = open(outFile, "w")
    findAS(geneDict, outASType=asType, anno=True, out=out)
    findAS(novelDict, outASType=asType, anno=False, out=out)
    out.close()

def main():
    findASmain(asType=asType, gpeFile=gpeFile, sourceFile=sourceFile, outFile=outFile)

if __name__ == '__main__':
    main()

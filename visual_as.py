#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: visual_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:16:55
Last modified: 2021-04-29 16:16:55
'''
# from commonFuncs import *
from commonObjs import *
from multiprocessing import Pool
from visual_as_functions import *

def visual_as(dataObj=None, targetGenes=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in project {} sample {}...".format(projectName, sampleName)
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in project {} sample {}...".format(projectName, sampleName)
        targetGenes = processTargetGenes(targetGenes)

    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    gpeObj = GenePredObj(refParams.ref_gpe, False)
    resolveDir(os.path.join(baseDir, "isoViewer"))
    tgsIsoFile = os.path.join(baseDir, "as_events", "ordinary_as", "isoformGrouped.AS.confident.bed12+")
    readsDict = {}
    selectedGenes = {}
    with open(tgsIsoFile) as f:
        for line in f:
            readStruc = ReadLineStruc(line)
            if readStruc.chrom not in readsDict:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom] = {readStruc.strand: {readStruc.geneName: gene2reads}}
            elif readStruc.strand not in readsDict[readStruc.chrom]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand] = {readStruc.geneName: gene2reads}
            elif readStruc.geneName not in readsDict[readStruc.chrom][readStruc.strand]:
                gene2reads = Gene2Reads(readStruc.geneName)
                gene2reads.update(readStruc)
                readsDict[readStruc.chrom][readStruc.strand].update({readStruc.geneName: gene2reads})
            else:
                readsDict[readStruc.chrom][readStruc.strand][readStruc.geneName].update(readStruc)
            if targetGenes == None:
                selectedGenes.update({readStruc.geneName: ""})
            else:
                if readStruc.geneName in targetGenes:
                    selectedGenes.update({readStruc.geneName: ""})

    try:
        pool = Pool(processes=dataObj.single_run_threads)
        for chrom in readsDict:
            for strand in readsDict[chrom]:
                for geneName in readsDict[chrom][strand]:
                    if geneName in selectedGenes:
                        gene2readsObj = readsDict[chrom][strand][geneName]
                        sampleTargetGenePickle = pickle.dumps(gene2readsObj)
                        gpeTargetGenePickle = pickle.dumps(gpeObj.geneName2gpeObj[gene2readsObj.geneName])
                        pool.apply_async(parallelPlotterAnno, (gene2readsObj.geneName, gpeTargetGenePickle,
                                                               sampleTargetGenePickle, dataObj, dirSpec))
                        # if gene2readsObj.geneName in gpeObj.geneName2gpeObj:
                        #     gpeTargetGenePickle = pickle.dumps(gpeObj.geneName2gpeObj[gene2readsObj.geneName])
                        #     pool.apply_async(parallelPlotterAnno, (gene2readsObj.geneName, gpeTargetGenePickle,
                        #                                            sampleTargetGenePickle, tgsSample, dirSpec))
                        # else:
                        #     pool.apply_async(parallelPlotterNovel, (gene2readsObj.geneName, sampleTargetGenePickle,
                        #                                             tgsSample, dirSpec))
        pool.close()
        pool.join()
    except Exception as e:
        print e
    os.chdir(prevDir)
    print getCurrentTime() + " Visualize the gene structure compared to the reference genome for project {} sample {} done!".format(projectName, sampleName)

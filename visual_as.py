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
import glob

def visual_as_merge(dataToProcess=None, targetGenes=None, refParams=None, dirSpec=None):
    isoViewerDir = os.path.join(dirSpec.out_dir, "isoViewer_sample_merged")
    resolveDir(isoViewerDir)
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in all samples..."
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in all samples..."
        targetGenes = processTargetGenes(targetGenes)
    refGpeObj = GenePredObj(refParams.ref_gpe, False)
    for geneName in targetGenes:
        mergedData = {}
        for dataObj in dataToProcess:
            tmpDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "isoViewer", geneName)
            sampleName = "{}_{}".format(dataObj.project_name, dataObj.sample_name)
            if validateDir(tmpDir):
                tmpDict = {"tgs": [], "ngs": {}}
                isoformGff = os.path.abspath(os.path.join(tmpDir, "{}.iFLAS.gff".format(geneName)))
                isoformGpe = os.path.abspath(os.path.join(tmpDir, "{}.iFLAS.gpe".format(geneName)))
                ngsSamList = glob.glob(os.path.join(tmpDir, "repeat*.{}.sam".format(geneName)))
                if validateFile(isoformGff):
                    tmpDict["tgs"] = [isoformGff, isoformGpe]
                for i in ngsSamList:
                    if validateFile(os.path.abspath(i)):
                        repeat = os.path.basename(i).split(".")[0]
                        tmpDict["ngs"].update({repeat: os.path.abspath(i)})
                mergedData.update({sampleName: tmpDict})
        if mergedData:
            resolveDir(os.path.join(isoViewerDir, geneName))

            isoformGpes = [mergedData[x]["tgs"][1] for x in mergedData]
            minPos = []
            maxPos = []
            for x in isoformGpes:
                gpeObjs = GenePredObj(x, False).geneName2gpeObj[geneName]
                minPos = [gpeObj.txStart for gpeObj in gpeObjs]
                maxPos = [gpeObj.txEnd for gpeObj in gpeObjs]

            refGeneObj = refGpeObj.geneName2gpeObj[geneName]

            targetGeneChrom = refGeneObj[0].chrom
            targetGeneStrand = refGeneObj[0].strand
            plotMinpos = min([x.txStart for x in refGeneObj] + minPos)
            plotMaxpos = max([x.txEnd for x in refGeneObj] + maxPos)

            # Gene model section
            geneModelGPE = "{}.gpe".format(geneName)
            geneModelGTF = "{}.gtf".format(geneName)
            geneModelGFF = "{}.gff".format(geneName)

            geneModelGPEOut = open(geneModelGPE, "w")
            for i in refGeneObj:
                print >> geneModelGPEOut, i
            geneModelGPEOut.close()
            targetGeneRegion = "{}:{}-{}".format(targetGeneChrom, plotMinpos, plotMaxpos)
            os.system("genePredToGtf file {} {} -source=iFLAS".format(geneModelGPE, geneModelGTF))
            os.system("gene_model_to_splicegraph.py -m {} -o {} 2>/dev/null".format(geneModelGTF, geneModelGFF))
            sectionToPlot = []
            geneModelHidePlot = PlotSection(section_name="[GeneModelGraph]", source_file=geneModelGTF,
                                            gene_name=geneName, relative_size=5.0,
                                            title_string="Gene Model for %gene", hide=True)
            sectionToPlot.append(geneModelHidePlot)
            geneModelItemCount = getFileRowCounts(geneModelGPE)
            geneModelRelativeSize = resizeTrackRatio(geneModelItemCount)
            geneModelVisiblePlot = PlotSection(section_name="[GeneModelIsoformsGraph]", plot_type="isoforms",
                                               source_file=geneModelGFF, relative_size=geneModelRelativeSize,
                                               title_string="Gene Model for %gene [{}({})]".format(
                                                   targetGeneRegion, targetGeneStrand))
            sectionToPlot.append(geneModelVisiblePlot)
            isoItemCounts = 0
            for sampleName in mergedData:
                tmpSample = mergedData[sampleName]
                postCorrIsoItemCounts = getFileRowCounts(tmpSample["tgs"][1])
                isoItemCounts += postCorrIsoItemCounts
                postCorrIsoRelativeSize = resizeTrackRatio(postCorrIsoItemCounts)
                postCorrIsoPlotType = "isoforms"
                postCorrIsoPlot = PlotSection(section_name="[AllReadsCollapse_{}]".format(sampleName), plot_type=postCorrIsoPlotType,
                                              source_file=tmpSample["tgs"][0], relative_size=postCorrIsoRelativeSize,
                                              title_string="Corrected isoforms and AS events in {} from {} data".format(geneName, sampleName))
                sectionToPlot.append(postCorrIsoPlot)
                for repeatName in tmpSample["ngs"]:
                    ngsSam = tmpSample["ngs"][repeatName]
                    ngsPlot = PlotSection(section_name="[Reads_{}_{}]".format(repeatName, sampleName), plot_type="read_depth",
                                            source_file=ngsSam, relative_size=5.0,
                                            title_string="{} Read Coverage in {} {}".format(geneName, sampleName, repeatName))
                    sectionToPlot.append(ngsPlot)
            figOut = geneName + ".pdf"
            cfgOut = open(geneName + ".cfg", "w")
            majorItemCount = geneModelItemCount + isoItemCounts
            figHeight = 20 if majorItemCount <= 50 else 30 if majorItemCount <= 150 else 40
            mainSec = MainSection(fout=figOut, height=figHeight)
            print >> cfgOut, mainSec.printStr()
            for sec in sectionToPlot:
                print >> cfgOut, sec.printStr()
            cfgOut.close()
            os.system("plotter.py {}.cfg 2>/dev/null".format(geneName))
    if targetGenes == None:
        print getCurrentTime() + " Visualize all gene structure compared to the reference genome in all samples done!"
    else:
        print getCurrentTime() + " Visualize selected gene structure compared to the reference genome in all samples done!"

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
                        pool.apply_async(parallelPlotter, (gene2readsObj.geneName, gpeTargetGenePickle,
                                                               sampleTargetGenePickle, dataObj, dirSpec))
        pool.close()
        pool.join()
    except Exception as e:
        print e
    os.chdir(prevDir)
    print getCurrentTime() + " Visualize the gene structure compared to the reference genome for project {} sample {} done!".format(projectName, sampleName)

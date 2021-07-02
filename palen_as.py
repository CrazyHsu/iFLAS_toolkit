#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: palen_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:25
Last modified: 2021-04-29 16:20:25
'''
from commonFuncs import *
from commonObjs import *
import pandas as pd
import scipy.stats as stats
import itertools

def getIsosFromAsFile(asFile, asType="IR"):
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")
        return inclusionIsos, exclusionIsos

def polyaLenCalling(dataObj=None, refParams=None, dirSpec=None):
    baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
    fast5Dir = dataObj.data_location
    nanoporeFileType = checkFast5Files(fast5Dir)
    if nanoporeFileType == "fast5":
        basecallingOut = os.path.join(baseDir, "preprocess", "nanopore", "guppy_basecalling")
        nanoporeRawFq = os.path.join(baseDir, "preprocess", "nanopore", "rawFlnc.fq")

        cmd = "nanopolish index --directory={} --sequencing-summary={}/sequencing_summary.txt {}".format(fast5Dir, basecallingOut, nanoporeRawFq)
        subprocess.call(cmd, shell=True)
        mappedBam = os.path.join(baseDir, "mapping", "flnc.mm2.sorted.bam")

        cmd = "nanopolish polya --threads={} --reads={} --bam={} --genome={} > polya_results.tsv".format(dataObj.single_run_threads, nanoporeRawFq, mappedBam, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
        dataObj.polya_location = os.path.join(os.getcwd(), "polya_results.tsv")
        cmd = '''grep 'PASS' polya_results.tsv | cut -f 9 | hist.R -x='poly(A) tail Length' -y=Density -b=1 -d -x1=0 -x2=100 -p=polyaTailLength.pdf 2>/dev/null'''
        subprocess.call(cmd, shell=True)
    else:
        raise Exception("Please specify the correct directory contain fast5 files!")

def getASpairedIsoforms(asFile, collapsedGroupFile, isoformFile, asType="SE", filterByCount=0, mergeByJunc=True):
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    isoBedObj = BedFile(isoformFile, type="bed12+")
    asPairs = {}
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            if mergeByJunc:
                incJuncCombDict = {}
                excJuncCombDict = {}
                for incIso in inclusionIsos:
                    incIsoObj = isoBedObj.reads[incIso]
                    if len(incIsoObj.introns) > 1:
                        if incIsoObj.juncChain not in incJuncCombDict:
                            incJuncCombDict[incIsoObj.juncChain] = [incIso]
                        else:
                            incJuncCombDict[incIsoObj.juncChain].append(incIso)
                    else:
                        if "monoExon" not in incJuncCombDict:
                            incJuncCombDict["monoExon"] = [incIso]
                        else:
                            incJuncCombDict["monoExon"].append(incIso)
                for excIso in exclusionIsos:
                    excIsoObj = isoBedObj.reads[excIso]
                    if len(excIsoObj.introns) > 1:
                        if excIsoObj.juncChain not in excJuncCombDict:
                            excJuncCombDict[excIsoObj.juncChain] = [excIso]
                        else:
                            excJuncCombDict[excIsoObj.juncChain].append(excIso)
                    else:
                        if "monoExon" not in excJuncCombDict:
                            excJuncCombDict["monoExon"] = [excIso]
                        else:
                            excJuncCombDict["monoExon"].append(excIso)

                for item in itertools.product(incJuncCombDict.values(), excJuncCombDict.values()):
                    newInclusionIsos = [x for x in item[0] if collapsedTrans2reads[x] >= filterByCount]
                    newExclusionIsos = [x for x in item[1] if collapsedTrans2reads[x] >= filterByCount]
                    if len(newInclusionIsos) == 0 or len(newExclusionIsos) == 0:
                        continue
                    '''the isoform PB.x.x split by "." to determine gene PB.x'''
                    inclusionGene = [".".join(x.split(".")[:2]) for x in newInclusionIsos]
                    exclusionGene = [".".join(x.split(".")[:2]) for x in newExclusionIsos]
                    uniqGene = list(set(inclusionGene+exclusionGene))
                    if len(uniqGene) == 1:
                        if uniqGene[0] not in asPairs:
                            asPairs[uniqGene[0]] = [[newInclusionIsos, newExclusionIsos, records[3]]]
                        else:
                            asPairs[uniqGene[0]].append([newInclusionIsos, newExclusionIsos, records[3]])
            else:
                newInclusionIsos = [x for x in inclusionIsos if len(collapsedTrans2reads[x]) >= filterByCount]
                newExclusionIsos = [x for x in exclusionIsos if len(collapsedTrans2reads[x]) >= filterByCount]
                if len(newInclusionIsos) == 0 or len(newExclusionIsos) == 0:
                    continue
                inclusionGene = [".".join(x.split(".")[:2]) for x in newInclusionIsos]
                exclusionGene = [".".join(x.split(".")[:2]) for x in newExclusionIsos]
                uniqGene = list(set(inclusionGene + exclusionGene))
                if len(uniqGene) == 1:
                    if uniqGene[0] not in asPairs:
                        asPairs[uniqGene[0]] = [[newInclusionIsos, newExclusionIsos, records[3]]]
                    else:
                        asPairs[uniqGene[0]].append([newInclusionIsos, newExclusionIsos, records[3]])
                # for item in itertools.product(inclusionIsos, exclusionIsos):
                #     inclusionReads = collapsedTrans2reads[item[0]]
                #     exclusionReads = collapsedTrans2reads[item[1]]
                #     if filterByCount:
                #         if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
                #             continue
                #     '''the isoform PB.x.x split by "." to determine gene PB.x'''
                #     inclusionGene = ".".join(item[0].split(".")[:2])
                #     exclusionGene = ".".join(item[1].split(".")[:2])
                #     if inclusionGene == exclusionGene:
                #         if inclusionGene not in asPairs:
                #             asPairs[inclusionGene] = [[[item[0]], [item[1]], records[3]]]
                #         else:
                #             asPairs[inclusionGene].append([[item[0]], [item[1]], records[3]])
        return asPairs

def getPalenAS(flncReads2Palen, isoformFile, readsFile, collapsedTrans2reads=None, asPairs=None, filterByCount=10, mergeByJunc=False):
    prevDir = os.getcwd()
    if not mergeByJunc:
        resolveDir("no_mergeByJunc")
    else:
        resolveDir("mergeByJunc")

    isoBedObj = BedFile(isoformFile, type="bed12+")
    readBedObj = BedFile(readsFile, type="bed12+")
    for asType in asPairs:
        sigFile = "{}.palenAndAS.sig.bed12+".format(asType)
        sigOut = open(sigFile, "w")
        count = 0
        for gene in asPairs[asType]:
            for item in asPairs[asType][gene]:
                inclusionIsos = item[0]
                exclusionIsos = item[1]
                inclusionReads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in inclusionIsos])
                exclusionReads = itertools.chain.from_iterable([collapsedTrans2reads[x] for x in exclusionIsos])
                if asType == "IR":
                    tmpReads = []
                    retention = [map(int, re.split("[:|-]", item[2])[-2:])]
                    for x in inclusionReads:
                        overlap = getOverlapOfTuple(readBedObj.reads[x].exons, retention)
                        if getBlockLength(overlap) == getBlockLength(retention):
                            tmpReads.append(x)
                    inclusionReads = tmpReads
                inclusionReads2palen = [[x, flncReads2Palen[x]] for x in inclusionReads if x in flncReads2Palen]
                exclusionReads2palen = [[x, flncReads2Palen[x]] for x in exclusionReads if x in flncReads2Palen]
                if len(inclusionReads2palen) < filterByCount or len(exclusionReads2palen) < filterByCount: continue
                aPalen = map(float, [x[1] for x in inclusionReads2palen])
                bPalen = map(float, [x[1] for x in exclusionReads2palen])
                stat_val1, p_val1 = stats.ttest_ind(aPalen, bPalen)
                stat_val2, p_val2 = stats.kruskal(aPalen, bPalen)
                if float(p_val1) <= 0.001 and float(p_val2) <= 0.001:
                    count += 1
                    fileOut = "{}_{}_{}.txt".format(asType, gene, count)
                    out = open(fileOut, "w")
                    for i in inclusionReads2palen:
                        print >>out, "\t".join(map(str, [gene, asType, "_".join(inclusionIsos), i[0], i[1]]))
                    for j in exclusionReads2palen:
                        print >>out, "\t".join(map(str, [gene, asType, "_".join(exclusionIsos), j[0], j[1]]))
                    out.close()
                    # cmd = "Rscript violin_paired.R {}".format(fileOut)
                    # subprocess.call(cmd, shell=True)

                    # isoBedObj = BedFile(isoformFile, type="bed12+")
                    inclusionIsosObj = [isoBedObj.reads[x] for x in inclusionIsos]
                    exclusionIsosObj = [isoBedObj.reads[x] for x in exclusionIsos]
                    incSortedIsos = sorted(inclusionIsosObj, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                    excSortedIsos = sorted(exclusionIsosObj, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
                    incRepIso = copy.copy(incSortedIsos[0])
                    excRepIso = copy.copy(excSortedIsos[0])
                    incRepIso.name = "_".join(inclusionIsos)
                    excRepIso.name = "_".join(exclusionIsos)

                    print >> sigOut, str(incRepIso) + "\t" + "\t".join(["_".join(inclusionIsos)+"+"+"_".join(exclusionIsos), "inclusion", ",".join(map(str, aPalen)), str(len(aPalen))])
                    print >> sigOut, str(excRepIso) + "\t" + "\t".join(["_".join(inclusionIsos)+"+"+"_".join(exclusionIsos), "exclusion", ",".join(map(str, bPalen)), str(len(bPalen))])
        sigOut.close()
    os.chdir(prevDir)

def getPalenAPA(paClusterBed, flncReads2Palen, filterByCount, palenAPA):
    gene2paReads = {}
    with open(paClusterBed) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            reads = infoList[3].split(",")
            if len(reads) < filterByCount: continue
            gene = infoList[8]
            paSite = "{}_{}".format(infoList[0], infoList[3])
            if gene not in gene2paReads:
                gene2paReads[gene] = {paSite: reads}
            else:
                gene2paReads[gene].update({paSite: reads})
    out = open(palenAPA, "w")
    for gene in gene2paReads:
        if len(gene2paReads[gene].keys()) < 2: continue
        validPaSites = {}
        for item in itertools.combinations(gene2paReads[gene].keys(), 2):
            aPalen = [float(flncReads2Palen[x]) for x in gene2paReads[gene][item[0]] if x in flncReads2Palen]
            bPalen = [float(flncReads2Palen[x]) for x in gene2paReads[gene][item[1]] if x in flncReads2Palen]
            if len(aPalen) < filterByCount or len(bPalen) < filterByCount: continue
            stat_val1, p_val1 = stats.ttest_ind(aPalen, bPalen)
            stat_val2, p_val2 = stats.kruskal(aPalen, bPalen)
            if float(p_val1) <= 0.001 and float(p_val2) <= 0.001:
                if item[0] not in validPaSites:
                    validPaSites[item[0]] = aPalen
                if item[1] not in validPaSites:
                    validPaSites[item[1]] = bPalen

        if len(validPaSites) < 2: continue
        paNumStr = ";".join(map(str, [len(validPaSites[x]) for x in validPaSites.keys()]))
        paLenStr = ";".join([",".join(map(str, validPaSites[x])) for x in validPaSites.keys()])
        print >> out, "\t".join([gene, ";".join(validPaSites.keys()), paNumStr, paLenStr])
    out.close()

def palen_as(dataObj=None, refParams=None, dirSpec=None, filterByCount=10, dataToProcess=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify functional poly(A) tail length related to AS for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "palenAS"))

    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    readsFile = os.path.join(baseDir, "mapping", "flnc.addCVandID.bed12+")
    isoformFile = os.path.join(baseDir, "refine", "isoformGrouped.bed12+")
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)

    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "IR.confident.bed6+")
    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "SE.confident.bed12+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A3SS.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A5SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, isoformFile, asType="SE", mergeByJunc=True)
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, isoformFile, asType="IR", mergeByJunc=True)
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, isoformFile, asType="A5SS", mergeByJunc=True)
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, isoformFile, asType="A3SS", mergeByJunc=True)
    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    palenFile = dataObj.polya_location
    palenFileList = []
    if palenFile == None:
        for tmpObj in dataToProcess:
            polyaLenCalling(dataObj=tmpObj, refParams=refParams, dirSpec=dirSpec)
            palenFileList.append(tmpObj.polya_location)
    elif isinstance(palenFile, dict):
        palenFileList.extend(palenFile.values())
        for tmpObj in dataToProcess:
            if tmpObj.sample_name not in palenFile:
                polyaLenCalling(dataObj=tmpObj, refParams=refParams, dirSpec=dirSpec)
                palenFileList.append(tmpObj.polya_location)
    else:
        palenFileList.append(palenFile)

    palen = pd.concat([pd.read_csv(x, sep="\t", dtype={"contig": "string"}) for x in palenFileList])
    palen_pass = palen.loc[palen.qc_tag == "PASS", ]
    flncReads2Palen = dict(zip(palen_pass.readname, palen_pass.polya_length))

    getPalenAS(flncReads2Palen, isoformFile, readsFile, collapsedTrans2reads=collapsedTrans2reads, asPairs=asPairs, filterByCount=filterByCount, mergeByJunc=True)
    # getPalenAS(flncReads2Palen, dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, collapsedTrans2reads=collapsedTrans2reads, filterByCount=filterByCount, merged=False)
    # getPalenAS(flncReads2Palen, dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, mergedIsoDict=mergedIsoDict, filterByCount=filterByCount, merged=True)

    os.chdir(prevDir)
    print getCurrentTime() + " Identify functional poly(A) tail length related to AS for project {} entry {} done!".format(projectName, sampleName)

    print getCurrentTime() + " Identify functional poly(A) tail length related to APA for project {} sample {}...".format(projectName, sampleName)
    resolveDir(os.path.join(baseDir, "palenAS", "palenAPA"))
    paClusterBed = os.path.join(baseDir, "as_events", "pa", "PA.bed6+")
    palenAPA = "apaRelatedPalen.txt"
    getPalenAPA(paClusterBed, flncReads2Palen, filterByCount, palenAPA)
    os.chdir(prevDir)

    # palen_apa(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, filterByCount=filterByCount, dataToProcess=dataToProcess)
    print getCurrentTime() + " Identify functional poly(A) tail length related to APA for project {} sample {} done".format(projectName, sampleName)

def palen_apa(dataObj=None, refParams=None, dirSpec=None, filterByCount=10, dataToProcess=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify functional poly(A) tail length related to APA for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "palenAPA"))

    palenFile = dataObj.polya_location
    palenFileList = []
    if palenFile == None:
        for tmpObj in dataToProcess:
            polyaLenCalling(dataObj=tmpObj, refParams=refParams, dirSpec=dirSpec)
            palenFileList.append(tmpObj.polya_location)
    elif isinstance(palenFile, dict):
        palenFileList.extend(palenFile.values())
        for tmpObj in dataToProcess:
            if tmpObj.sample_name not in palenFile:
                polyaLenCalling(dataObj=tmpObj, refParams=refParams, dirSpec=dirSpec)
                palenFileList.append(tmpObj.polya_location)
    else:
        palenFileList.append(palenFile)

    palen = pd.concat([pd.read_csv(x, sep="\t", dtype={"contig": "string"}) for x in palenFileList])
    palen_pass = palen.loc[palen.qc_tag == "PASS", ]
    flncReads2Palen = dict(zip(palen_pass.readname, palen_pass.polya_length))
    paClusterBed = os.path.join(baseDir, "as_events", "pa", "PA.bed6+")
    palenAPA = "apaRelatedPalen.txt"
    getPalenAPA(paClusterBed, flncReads2Palen, filterByCount, palenAPA)

    os.chdir(prevDir)
    print getCurrentTime() + " Identify functional poly(A) tail length related to APA for project {} sample {}...".format(projectName, sampleName)

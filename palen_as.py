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

def mergeIso(isoBedFile, collapsedTrans2reads):
    isoBedObj = BedFile(isoBedFile, type="bed12+")
    gene2iso = {}
    for iso in isoBedObj.reads:
        if isoBedObj.reads[iso].otherList[0] not in gene2iso:
            gene2iso[isoBedObj.reads[iso].otherList[0]] = []
        gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])
    mergedIsoDict = {}
    mergedIso2ReadsGroupOut = open("mergedIso2Reads.group.txt", "w")
    mergedIso2ReadsBed = open("mergedIso2Reads.bed12", "w")
    for gene in gene2iso:
        tmpDict = {}
        for iso in gene2iso[gene]:
            if iso.juncChain not in tmpDict:
                tmpDict[iso.juncChain] = [iso]
            else:
                tmpDict[iso.juncChain].append(iso)

        for tmp in tmpDict:
            isos = tmpDict[tmp]
            mergedIsoName = "+".join([x.name for x in isos])
            reads = list(itertools.chain.from_iterable([collapsedTrans2reads[x.name] for x in isos]))
            sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
            repIso = sortedIsos[0]
            repIso.name = mergedIsoName
            mergedIsoDict[mergedIsoName] = [reads, repIso, gene]
            print >> mergedIso2ReadsGroupOut, "{}\t{}".format(mergedIsoName, ",".join(reads))
            print >> mergedIso2ReadsBed, str(repIso)
    mergedIso2ReadsBed.close()
    mergedIso2ReadsGroupOut.close()
    return mergedIsoDict

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

def relationshipBetweenPAlenAndAS(asFile, collapsedGroupFile, flncReads2Palen, asType="SE", filterByCount=5, sigFile=None, sigPicOutDir=None):
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    # flncReads2Palen = getDictFromFile(flncFile, sep=",", valueCol=5)
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    # isoformDict = getDictFromFile(isoformBed, keyCol=4, sep="\t")

    asType = asType.upper()
    out = open(sigFile, "w")
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            elif asType == "PA":
                inclusionIsos = []
                exclusionIsos = []
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            for item in itertools.product(inclusionIsos, exclusionIsos):
                inclusionReads = collapsedTrans2reads[item[0]]
                exclusionReads = collapsedTrans2reads[item[1]]
                inclusionReads2palen = [int(flncReads2Palen[i]) for i in inclusionReads if i in flncReads2Palen]
                exclusionReads2palen = [int(flncReads2Palen[i]) for i in exclusionReads if i in flncReads2Palen]
                if filterByCount:
                    if len(inclusionReads2palen) < filterByCount or len(exclusionReads2palen) < filterByCount:
                        continue
                stat_val, p_val = stats.ttest_ind(inclusionReads2palen, exclusionReads2palen)
                if float(p_val) <= 0.001:
                    print >> out, "\t".join(records) + "\t" + "\t".join(
                        ["inclusion", ",".join(map(str, inclusionReads2palen)), str(len(inclusionReads2palen)), "_".join(item)])
                    print >> out, "\t".join(records) + "\t" + "\t".join(
                        ["exclusion", ",".join(map(str, exclusionReads2palen)), str(len(exclusionReads2palen)), "_".join(item)])
    out.close()

    if os.stat(sigFile).st_size != 0:
        plotGeneStructurePath = os.path.join(scriptDir, "plotGeneStructure.R")
        cmd = "Rscript {} -sig={} -od={}".format(plotGeneStructurePath, sigFile, sigPicOutDir)
        # subprocess.call(cmd, shell=True)


def getPalenAS(flncReads2Palen, dataObj=None, refParams=None, dirSpec=None, collapsedTrans2reads=None, mergedIsoDict=None, filterByCount=10, merged=False):
    baseDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name)
    gene2readPalen = {}
    prevDir = os.getcwd()
    if not merged:
        resolveDir("no_merge")
    else:
        resolveDir("merged")

    if not merged:
        for isoform in collapsedTrans2reads:
            geneId = ".".join(isoform.split(".")[0:2])
            reads = collapsedTrans2reads[isoform]
            readsWithPolyaLen = [i for i in reads if i in flncReads2Palen]
            reads2palen = dict(zip(readsWithPolyaLen, [int(flncReads2Palen[i]) for i in readsWithPolyaLen]))
            if geneId not in gene2readPalen:
                gene2readPalen[geneId] = {isoform: reads2palen}
            else:
                gene2readPalen[geneId].update({isoform: reads2palen})
    else:
        for isoform in mergedIsoDict:
            geneId = mergedIsoDict[isoform][2]
            reads = mergedIsoDict[isoform][0]
            readsWithPolyaLen = [i for i in reads if i in flncReads2Palen]
            reads2palen = dict(zip(readsWithPolyaLen, [int(flncReads2Palen[i]) for i in readsWithPolyaLen]))
            if geneId not in gene2readPalen:
                gene2readPalen[geneId] = {isoform: reads2palen}
            else:
                gene2readPalen[geneId].update({isoform: reads2palen})

    count = 0
    for g in gene2readPalen:
        # allReadsInGene = sum([len(gene2readPalen[g][i]) for i in gene2readPalen[g]])
        for a, b in itertools.combinations(gene2readPalen[g].keys(), 2):
            if len(gene2readPalen[g][a].keys()) >= filterByCount and len(
                    gene2readPalen[g][b].keys()) >= filterByCount:
                aPalen = map(int, gene2readPalen[g][a].values())
                bPalen = map(int, gene2readPalen[g][b].values())
                stat_val1, p_val1 = stats.ttest_ind(aPalen, bPalen)
                stat_val2, p_val2 = stats.kruskal(aPalen, bPalen)
                if float(p_val1) <= 0.001 and float(p_val2) <= 0.001:
                    # print "\t".join(map(str, [allReadsInGene, g, a, len(aPalen), np.mean(aPalen), np.median(aPalen), b, len(bPalen), np.mean(bPalen), np.median(bPalen)]))
                    if merged:
                        count += 1
                        fileOut = "{}_{}.txt".format(g, count)
                    else:
                        fileOut = "_".join([g, a, b]) + ".txt"
                    out = open(fileOut, "w")
                    for i in gene2readPalen[g][a]:
                        print >> out, "\t".join(map(str, [g, a, i, gene2readPalen[g][a][i]]))
                    for j in gene2readPalen[g][b]:
                        print >> out, "\t".join(map(str, [g, b, j, gene2readPalen[g][b][j]]))
                    out.close()
                    # plotGeneStructurePath = os.path.join(currentFileDir, "violin_paired.R")
                    # cmd = "Rscript {} {}".format(plotGeneStructurePath, fileOut)
                    cmd = "Rscript violin_paired.R {}".format(fileOut)
                    # subprocess.call(cmd, shell=True)

    if merged:
        cmd = "getSingleIr.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.IR.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getSE.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.SE.bed12+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getA5SS.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.A5SS.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        cmd = "getA3SS.pl -g {} ../mergedIso2Reads.bed12 > mergedIso.A3SS.bed6+".format(refParams.ref_gpe)
        subprocess.call(cmd, shell=True)
        irFile = os.path.join(os.getcwd(), "mergedIso.IR.bed6+")
        seFile = os.path.join(os.getcwd(), "mergedIso.SE.bed12+")
        a5ssFile = os.path.join(os.getcwd(), "mergedIso.A5SS.bed6+")
        a3ssFile = os.path.join(os.getcwd(), "mergedIso.A3SS.bed6+")
        collapsedGroupFile = os.path.join(os.getcwd(), "../mergedIso2Reads.group.txt")
    else:
        seFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "SE.confident.bed12+")
        irFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "IR.confident.bed6+")
        a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "A5SS.confident.bed6+")
        a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "A3SS.confident.bed6+")
        collapsedGroupFile = os.path.join(baseDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    relationshipBetweenPAlenAndAS(irFile, collapsedGroupFile, flncReads2Palen, asType="IR", filterByCount=5,
                                  sigFile="IR.palenAndAS.sig.bed12+", sigPicOutDir="IR.sig_pics")
    relationshipBetweenPAlenAndAS(seFile, collapsedGroupFile, flncReads2Palen, asType="SE", filterByCount=5,
                                  sigFile="SE.palenAndAS.sig.bed12+", sigPicOutDir="SE.sig_pics")
    relationshipBetweenPAlenAndAS(a5ssFile, collapsedGroupFile, flncReads2Palen, asType="A5SS", filterByCount=5,
                                  sigFile="A5SS.palenAndAS.sig.bed12+", sigPicOutDir="A5SS.sig_pics")
    relationshipBetweenPAlenAndAS(a3ssFile, collapsedGroupFile, flncReads2Palen, asType="A3SS", filterByCount=5,
                                  sigFile="A3SS.palenAndAS.sig.bed12+", sigPicOutDir="A3SS.sig_pics")
    os.chdir(prevDir)


def palen_as(dataObj=None, refParams=None, dirSpec=None, filterByCount=10, sampleMerged=False):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify functional poly(A) tail length related to AS for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "palenRelatedAS"))

    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    isoBedFile = os.path.join(baseDir, "collapse", "tofu.collapsed.assigned.unambi.bed12+")
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=2)
    mergedIsoDict = mergeIso(isoBedFile, collapsedTrans2reads)

    palenFile = dataObj.polya_location
    if palenFile == None or not validateFile(palenFile):
        polyaLenCalling(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)

    palen = pd.read_csv(palenFile, sep="\t")
    palen_pass = palen.loc[palen.qc_tag == "PASS", ]
    flncReads2Palen = dict(zip(palen_pass.readname, palen_pass.polya_length))

    getPalenAS(flncReads2Palen, dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, collapsedTrans2reads=collapsedTrans2reads, filterByCount=filterByCount, merged=False)
    getPalenAS(flncReads2Palen, dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, mergedIsoDict=mergedIsoDict, filterByCount=filterByCount, merged=True)

    os.chdir(prevDir)
    print getCurrentTime() + " Identify functional poly(A) tail length related to AS for project {} entry {} done!".format(projectName, sampleName)

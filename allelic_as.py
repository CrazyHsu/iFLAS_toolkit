#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: allelic_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:04
Last modified: 2021-04-29 16:20:04
'''
from commonFuncs import *
from ConfigParser import RawConfigParser
from multiprocessing import Pool
from scipy.stats import chi2_contingency
import StringIO, itertools, pysam, glob
import pandas as pd

min_bq = 13
ploidy = 2

def getASpairedIsoforms(asFile, collapsedGroupFile, asType="SE", filterByCount=0):
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=1)
    asPairs = {}
    with open(asFile) as f:
        for line in f.readlines():
            records = line.strip("\n").split("\t")
            if asType == "SE":
                inclusionIsos = records[15].split(",")
                exclusionIsos = records[17].split(",")
            # elif asType == "PA":
            #     inclusionIsos = []
            #     exclusionIsos = []
            else:
                inclusionIsos = records[7].split(",")
                exclusionIsos = records[9].split(",")

            for item in itertools.product(inclusionIsos, exclusionIsos):
                inclusionReads = collapsedTrans2reads[item[0]]
                exclusionReads = collapsedTrans2reads[item[1]]
                if filterByCount:
                    if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
                        continue
                '''the isoform PB.x.x split by "." to determine gene PB.x'''
                inclusionGene = ".".join(item[0].split(".")[:2])
                exclusionGene = ".".join(item[1].split(".")[:2])
                isoPair = [item[0], item[1]]
                if inclusionGene == exclusionGene:
                    if inclusionGene not in asPairs:
                        asPairs[inclusionGene] = [isoPair]
                    else:
                        asPairs[inclusionGene].append(isoPair)
        return asPairs

def makeAbundanceFile(groupFile, outFile=None):
    with open(groupFile) as f:
        cidInfo = {}
        for line in f.readlines():
            pbid, members = line.strip().split('\t')
            for cid in members.split(','):
                cidInfo[cid] = pbid
        if outFile:
            out = open(outFile, "w")
            print >> out, "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print >> out, "\t".join([i, "NA", "Y", "unique", cidInfo[i]])
            out.close()
        else:
            print "id\tlength\tis_fl\tstat\tpbid"
            for i in cidInfo:
                print "\t".join([i, "NA", "Y", "unique", cidInfo[i]])

def runPhaser(targetDir):
    curDir = os.getcwd()
    os.chdir(targetDir)
    configFile = "config"
    with open(configFile) as f:
        config_string = StringIO.StringIO("[dummy_section]\n" + f.read())
    config_parser = RawConfigParser()
    config_parser.readfp(config_string)
    strand = config_parser.get("dummy_section", "ref_strand")

    cmd = "minimap2 -ax splice fake.fasta ccs.fastq >ccs.sam 2>/dev/null"
    subprocess.call(cmd, shell=True)

    pysam.sort("-o", "ccs.sorted.bam", "ccs.sam", catch_stdout=False)
    pileupRes = pysam.mpileup("--min-BQ", str(min_bq), "-f", "fake.fasta", "-s", "ccs.sorted.bam")
    pileupResOut = open("ccs.mpileup", "w")
    pileupResOut.write(pileupRes)
    pileupResOut.close()

    # cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand {} -o phased.nopartial -n {} 1>nopartial.stdout.txt".format(
    #     strand, ploidy)
    # subprocess.call(cmd, shell=True)
    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} -o phased.partial -n {} 1>partial.stdout.txt".format(
        strand, ploidy)
    subprocess.call(cmd, shell=True)

    os.chdir(curDir)

def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    resultDict = {"partial": {}, "nopartial": {}}
    for isoIndex in range(len(isoPairs)):
        # partial
        if partial:
            partialHaplotype = "phased.partial.cleaned.human_readable.txt"
            if os.path.exists(partialHaplotype):
                partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
                indexPairs = list(itertools.combinations(partialDF.index, 2))
                if len(set(isoPairs[isoIndex]) & set(partialDF.columns)) == 2:
                    for pairIndex in range(len(indexPairs)):
                        subPartialDF = partialDF.reindex(index=indexPairs[pairIndex], columns=isoPairs[isoIndex])
                        if (pd.DataFrame.max(subPartialDF) >= 5).all() and (pd.DataFrame.min(subPartialDF) <= 3).all():
                            partialChi2Results = chi2_contingency(subPartialDF)
                            partialChi2Pvalue = partialChi2Results[1]
                            if partialChi2Pvalue <= 0.001:
                                combination = "_".join(map(str, ["partial", isoIndex, pairIndex]))
                                rows = subPartialDF.index
                                columns = subPartialDF.columns
                                resultDict["partial"].update({combination: dict(zip(rows, columns))})


        # no-partial
        if nopartial:
            nopartialHaplotype = "phased.nopartial.cleaned.human_readable.txt"
            if os.path.exists(nopartialHaplotype):
                nopartialDF = pd.read_csv(nopartialHaplotype, skiprows=[0], index_col=0, sep="\t")
                indexPairs = list(itertools.combinations(nopartialDF.index, 2))
                if len(set(isoPairs[isoIndex]) & set(nopartialDF.columns)) == 2:
                    for pairIndex in range(len(indexPairs)):
                        subNopartialDF = nopartialDF.reindex(index=indexPairs[pairIndex], columns=isoPairs[isoIndex])
                        nopartialChi2Results = chi2_contingency(subNopartialDF)
                        nopartialChi2Pvalue = nopartialChi2Results[1]
                        if nopartialChi2Pvalue <= 0.05:
                            combination = "_".join(map(str, ["nopartial", isoIndex, pairIndex]))
                            rows = subNopartialDF.index
                            columns = subNopartialDF.columns
                            resultDict["nopartial"].update({combination: dict(zip(rows, columns))})

        if not partial and not nopartial:
            raise Exception("You must specify either partial or nopartial file")
    return resultDict

def allelic_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify allelic-specific AS events for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "allelicSpecificRelatedAS"))

    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    # select_loci_to_phase = os.path.join(cDNA_Cupcake_dir, "phasing", "utils", "select_loci_to_phase.py")
    processedFa = os.path.join(baseDir, "filtration", "processed.fa")
    collapsedGff = os.path.join(baseDir, "filtration", "collapse", "tofu.collapsed.good.gff")

    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "SE.confident.bed12+")
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "IR.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A3SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, asType="SE")
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, asType="IR")
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, asType="A5SS")
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, asType="A3SS")

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    # generate "by_loci" directory
    # removeFiles()
    cmd = "select_loci_to_phase.py {} {} {} {} -c 10 1>/dev/null 2>&1".format(refParams.ref_genome, processedFa,
                                                                              collapsedGff, readStatFile)
    subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    multiResults = []
    pool = Pool(processes=dataObj.single_run_threads)
    # for i in lociDir:
    #     runPhaser(i, asPairs)
    for i in lociDir:
        singleRunRes = pool.apply_async(runPhaser, (i))
        multiResults.append(singleRunRes)
    pool.close()
    pool.join()

    resultDict = {}
    for i in lociDir:
        prevDir = os.getcwd()
        os.chdir(i)
        fakeGene = os.path.basename(i).split("_")[0]
        resultDict[fakeGene] = {}
        for asType in asPairs:
            if fakeGene not in asPairs[asType]:
                continue
            isoPairs = asPairs[asType][fakeGene]
            sigRelatedDict = relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=False, isoPairs=isoPairs)
            resultDict[fakeGene].update({asType: sigRelatedDict})
        os.chdir(prevDir)

    '''the structure of resultDict is resultDict = {fakeGene: {asType: {partialCategory: {combinationName: {haplotype1: PB.X.1, haplotype2: PB.X.2, ...}}}}}'''
    partialOut = open("partialAsRelatedHaplotype.txt", "w")
    # nopartialOut = open("nopartialAsRelatedHaplotype.txt", "w")
    for g in resultDict:
        for a in resultDict[g]:
            for p in resultDict[g][a]:
                for c in resultDict[g][a][p]:
                    for h in resultDict[g][a][p][c]:
                        if p == "partial":
                            print >> partialOut, "\t".join([g, a, h, resultDict[g][a][p][c][h], c])
    partialOut.close()
    os.chdir(prevDir)
    print getCurrentTime() + " Identify allelic-specific AS events for project {} sample {} done!".format(projectName, sampleName)

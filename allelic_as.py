#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: allelic_as.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-04-29 16:20:04
Last modified: 2021-04-29 16:20:04
'''
from commonFuncs import *
from commonObjs import *
from ConfigParser import RawConfigParser
from multiprocessing import Pool
from scipy.stats import chi2_contingency
import StringIO, itertools, glob
import pandas as pd

min_bq = 13
ploidy = 2

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
                    newInclusionIsos = [x for x in item[0] if len(collapsedTrans2reads[x]) >= filterByCount]
                    newExclusionIsos = [x for x in item[1] if len(collapsedTrans2reads[x]) >= filterByCount]
                    if len(newInclusionIsos) == 0 or len(newExclusionIsos) == 0:
                        continue
                    '''the isoform PB.x.x split by "." to determine gene PB.x'''
                    inclusionGene = [".".join(x.split(".")[:2]) for x in newInclusionIsos]
                    exclusionGene = [".".join(x.split(".")[:2]) for x in newExclusionIsos]
                    uniqGene = list(set(inclusionGene+exclusionGene))
                    if len(uniqGene) == 1:
                        if uniqGene[0] not in asPairs:
                            asPairs[uniqGene[0]] = [[newInclusionIsos, newExclusionIsos]]
                        else:
                            asPairs[uniqGene[0]].append([newInclusionIsos, newExclusionIsos])
            else:
                for item in itertools.product(inclusionIsos, exclusionIsos):
                    inclusionReads = collapsedTrans2reads[item[0]]
                    exclusionReads = collapsedTrans2reads[item[1]]
                    if filterByCount:
                        if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
                            continue
                    '''the isoform PB.x.x split by "." to determine gene PB.x'''
                    inclusionGene = ".".join(item[0].split(".")[:2])
                    exclusionGene = ".".join(item[1].split(".")[:2])
                    if inclusionGene == exclusionGene:
                        if inclusionGene not in asPairs:
                            asPairs[inclusionGene] = [[[item[0]], [item[1]]]]
                        else:
                            asPairs[inclusionGene].append([[item[0]], [item[1]]])

            # for item in itertools.product(inclusionIsos, exclusionIsos):
            #     inclusionReads = collapsedTrans2reads[item[0]]
            #     exclusionReads = collapsedTrans2reads[item[1]]
            #     if filterByCount:
            #         if len(inclusionReads) < filterByCount or len(exclusionReads) < filterByCount:
            #             continue
            #     '''the isoform PB.x.x split by "." to determine gene PB.x'''
            #     inclusionGene = ".".join(item[0].split(".")[:2])
            #     exclusionGene = ".".join(item[1].split(".")[:2])
            #     isoPair = [item[0], item[1]]
            #     if inclusionGene == exclusionGene:
            #         if inclusionGene not in asPairs:
            #             asPairs[inclusionGene] = [isoPair]
            #         else:
            #             asPairs[inclusionGene].append(isoPair)
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

    cmd = "samtools sort ccs.sam > ccs.sorted.bam 2>/dev/null; samtools mpileup --min-BQ {} -f fake.fasta -s ccs.sorted.bam > ccs.mpileup 2>/dev/null".format(min_bq)
    subprocess.call(cmd, shell=True)
    # pysam.sort("-o", "ccs.sorted.bam", "ccs.sam", catch_stdout=False)
    # pileupRes = pysam.mpileup("--min-BQ", str(min_bq), "-f", "fake.fasta", "-s", "ccs.sorted.bam")
    # pileupResOut = open("ccs.mpileup", "w")
    # pileupResOut.write(pileupRes)
    # pileupResOut.close()

    # cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand {} -o phased.nopartial -n {} 1>nopartial.stdout.txt".format(
    #     strand, ploidy)
    # subprocess.call(cmd, shell=True)
    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} " \
          "-o phased.partial -n {} --bhFDR 0.01 -p 0.01 1>partial.stdout.txt".format(strand, ploidy)
    subprocess.call(cmd, shell=True)

    os.chdir(curDir)

def relationshipBetweenAlleleSpecifiAndAS(partial=False, nopartial=False, isoPairs=None):
    resultDict = {"partial": {}, "nopartial": {}}
    for isoPair in isoPairs:
        if partial:
            partialHaplotype = "phased.partial.cleaned.human_readable.txt"
            if os.path.exists(partialHaplotype):
                partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
                if not set(isoPair[0] + isoPair[1]).issubset(set(partialDF.columns)): continue
                incIsoDf = partialDF.loc[:, isoPair[0]].sum(axis=1)
                excIsoDf = partialDF.loc[:, isoPair[1]].sum(axis=1)
                mergedDf = pd.concat([incIsoDf, excIsoDf], axis=1)
                if (pd.DataFrame.max(mergedDf) >= 5).all() and (pd.DataFrame.min(mergedDf) <= 3).all():
                    if len(mergedDf[(mergedDf.T != 0).any()]) == 1: continue
                    chi2Result = chi2_contingency(mergedDf)
                    chi2Pvalue = chi2Result[1]
                    if chi2Pvalue <= 0.001:
                        combination = "+".join(map(str, ["partial", "_".join(isoPair[0]), "_".join(isoPair[1])]))
                        haploInc = incIsoDf.idxmax(axis=0)
                        haploExc = excIsoDf.idxmax(axis=0)
                        resultDict["partial"].update({combination: dict(zip([haploInc, haploExc], ["_".join(isoPair[0]), "_".join(isoPair[1])]))})
        
        if nopartial:
            nopartialHaplotype = "phased.nopartial.cleaned.human_readable.txt"
            if os.path.exists(nopartialHaplotype):
                nopartialDF = pd.read_csv(nopartialHaplotype, skiprows=[0], index_col=0, sep="\t")
                if not set(isoPair[0] + isoPair[1]).issubset(set(nopartialDF.columns)): continue
                incIsoDf = nopartialDF.loc[:, isoPair[0]].sum(axis=1)
                excIsoDf = nopartialDF.loc[:, isoPair[1]].sum(axis=1)
                mergedDf = pd.concat([incIsoDf, excIsoDf], axis=1)
                if (pd.DataFrame.max(mergedDf) >= 5).all() and (pd.DataFrame.min(mergedDf) <= 3).all():
                    if len(mergedDf[(mergedDf.T != 0).any()]) == 1: continue
                    chi2Result = chi2_contingency(mergedDf)
                    chi2Pvalue = chi2Result[1]
                    if chi2Pvalue <= 0.001:
                        combination = "+".join(map(str, ["nopartial", "_".join(isoPair[0]), "_".join(isoPair[1])]))
                        haploInc = incIsoDf.idxmax(axis=0)
                        haploExc = excIsoDf.idxmax(axis=0)
                        resultDict["nopartial"].update({combination: dict(zip([haploInc, haploExc], [["inc", "_".join(isoPair[0])], ["exc", "_".join(isoPair[1])]]))})

        if not partial and not nopartial:
            raise Exception("You must specify either partial or nopartial file")
    return resultDict

def allelic_as(dataObj=None, refParams=None, dirSpec=None):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Identify allelic-specific AS events for project {} sample {}...".format(projectName, sampleName)
    prevDir = os.getcwd()
    baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
    resolveDir(os.path.join(baseDir, "allelicAS"))

    collapsedGroupFile = os.path.join(baseDir, "collapse", "tofu.collapsed.group.txt")
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    processedFq = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", dataObj.tgs_plat.lower(), "rawFlnc.fq")
    collapsedGff = os.path.join(baseDir, "collapse", "tofu.collapsed.good.gff")
    isoformFile = os.path.join(baseDir, "refine", "isoformGrouped.bed12+")

    seFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "SE.confident.bed12+")
    irFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "IR.confident.bed6+")
    a5ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(baseDir, "as_events", "ordinary_as", "PB", "A3SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, isoformFile, asType="SE")
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, isoformFile, asType="IR")
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, isoformFile, asType="A5SS")
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, isoformFile, asType="A3SS")

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    cmd = "rm -R by_loci"
    subprocess.call(cmd, shell=True)
    cmd = "select_loci_to_phase.py {} {} {} {} -c 10 1>/dev/null 2>&1"
    cmd = cmd.format(refParams.ref_genome, processedFq, collapsedGff, readStatFile)
    subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    # multiResults = []
    pool = Pool(processes=dataObj.single_run_threads)
    for i in lociDir:
        pool.apply_async(runPhaser, (i))
        # multiResults.append(singleRunRes)
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
    for gene in resultDict:
        for asType in resultDict[gene]:
            for partial in resultDict[gene][asType]:
                for comb in resultDict[gene][asType][partial]:
                    haplos = resultDict[gene][asType][partial][comb]
                    haplotype1, haplotype2 = haplos.keys()
                    print >>partialOut, "\t".join([gene, asType, haplotype1, haplos[haplotype1], haplotype2, haplos[haplotype2]])
                    # for haplo in resultDict[gene][asType][partial][comb]:
                    #     haplo2isos = resultDict[gene][asType][partial][comb][haplo]
                    #     print >> partialOut, "\t".join([gene, asType, haplo, haplo2isos])
    partialOut.close()
    os.chdir(prevDir)
    print getCurrentTime() + " Identify allelic-specific AS events for project {} sample {} done!".format(projectName, sampleName)

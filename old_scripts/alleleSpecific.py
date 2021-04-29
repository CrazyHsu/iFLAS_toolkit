#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: alleleSpecific.py.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-13 20:46:45
Last modified: 2020-05-13 20:46:45
'''
import argparse, glob, os, pysam
import StringIO
from Config import *
from commonObjs import *
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from relationshipBetweenAlleleSpecificAndAS import relationshipBetweenAlleleSpecifiAndAS

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

cDNA_Cupcake_dir = "/data/CrazyHsu_data/software/cDNA_Cupcake"
min_bq = 13
ploidy = 2

##############################################
def getASpairedIsoforms(asFile, collapsedGroupFile, asType="SE", filterByCount=0):
    collapsedTrans2reads = getDictFromFile(collapsedGroupFile, sep="\t", inlineSep=",", valueCol=1)
    asPairs = {}
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

def runPhaser(targetDir, asPairs):
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

    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand {} -o phased.nopartial -n {} 1>nopartial.stdout.txt".format(
        strand, ploidy)
    subprocess.call(cmd, shell=True)
    cmd = "run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand {} -o phased.partial -n {} 1>partial.stdout.txt".format(
        strand, ploidy)
    subprocess.call(cmd, shell=True)

    os.chdir(curDir)

def msa(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
        elif "?" == seq1[i] or "?" == seq2[i]:
            score += 0
        else:
            score -= 1
    return score

def getHapScore(seq):
    score = len(seq)
    if "?" in Counter(seq):
        score = score - Counter(seq)["?"]
    return score

def ASAP1():
    import pandas as pd
    from scipy.stats import chi2_contingency
    import itertools, glob, os
    lociDir = glob.glob("../../B73_Ki11_embryo/alleleSpecificRelatedAS/by_loci/*size*")
    for j in lociDir:
        partialHaplotype = os.path.join(j, "phased.partial.cleaned.human_readable.txt")
        if not os.path.exists(partialHaplotype): continue
        partialDF = pd.read_csv(partialHaplotype, skiprows=[0], index_col=0, sep="\t")
        partialDF = partialDF.loc[pd.DataFrame.max(partialDF, axis=1) >= 5, pd.DataFrame.max(partialDF) >= 5]
        if partialDF.empty: continue
        indexPairs = list(itertools.combinations(partialDF.index, 2))
        isoPairs = list(itertools.combinations(partialDF.columns, 2))
        for pairIndex in range(len(indexPairs)):
            for isoIndex in range(len(isoPairs)):
                subPartialDF = partialDF.loc[indexPairs[pairIndex], isoPairs[isoIndex]]
                try:
                    if (pd.DataFrame.max(subPartialDF) >= 5).all() and (pd.DataFrame.min(subPartialDF) <= 3).all():
                        partialChi2Results = chi2_contingency(subPartialDF)
                        partialChi2Pvalue = partialChi2Results[1]
                        if partialChi2Pvalue <= 0.001:
                            print isoPairs[isoIndex]
                except:
                    continue
                    # combination = "_".join(map(str, ["partial", isoIndex, pairIndex]))
                    # rows = subPartialDF.index
                    # columns = subPartialDF.columns[pd.DataFrame(subPartialDF).idxmax()]
                    # resultDict["partial"].update({combination: dict(zip(rows, columns))})

def ASAP():
    import operator
    lociDir = glob.glob("by_loci/*size*")
    for j in lociDir:
        partial = os.path.join(j, "phased.partial.human_readable.txt")
        if not os.path.exists(partial): continue
        partialDF = pd.read_csv(partial, skiprows=[0], index_col=0, sep="\t")
        partialDF = partialDF.loc[pd.DataFrame.max(partialDF, axis=1) >= 5, pd.DataFrame.max(partialDF) >= 5]
        if partialDF.empty: continue
        catergoryDict = {"cat1": {partialDF.index[0]: getHapScore(partialDF.index[0])}}
        for i in partialDF.index[1:]:
            if "cat2" not in catergoryDict:
                tmpSeq = max(catergoryDict["cat1"].iteritems(), key=operator.itemgetter(1))[0]
                if msa(i, tmpSeq) < 0:
                    catergoryDict.update({"cat2": {i: getHapScore(i)}})
                else:
                    catergoryDict["cat1"].update({i: getHapScore(i)})
            else:
                tmpSeq1 = max(catergoryDict["cat1"].iteritems(), key=operator.itemgetter(1))[0]
                tmpSeq2 = max(catergoryDict["cat2"].iteritems(), key=operator.itemgetter(1))[0]
                if msa(tmpSeq1, i) < 0 and msa(tmpSeq2, i) > 0:
                    catergoryDict["cat2"].update({i: getHapScore(i)})
                elif msa(tmpSeq1, i) > 0 and msa(tmpSeq2, i) < 0:
                    catergoryDict["cat1"].update({i: getHapScore(i)})
                else:
                    catergoryDict.update({"nocat": {i: 0}})

        if len(catergoryDict.keys()) == 1: continue
        for item1 in itertools.product(*[catergoryDict["cat1"], catergoryDict["cat2"]]):
            for item2 in list(itertools.combinations(partialDF.columns, 2)):
                subPartialDF = partialDF.reindex(index=item1, columns=item2)
                try:
                    if (pd.DataFrame.max(subPartialDF) >= 5).all() and (pd.DataFrame.min(subPartialDF) <= 3).all():
                        partialChi2Results = chi2_contingency(subPartialDF)
                        partialChi2Pvalue = partialChi2Results[1]
                        if partialChi2Pvalue <= 0.0001:
                            print item2
                except:
                    continue

        # isoPairs = list(itertools.combinations(partialDF.columns, 2))

        asIso = []
        for j in partialDF.columns:
            tmpCol = partialDF.loc[:, j]
            tmpCol = tmpCol[tmpCol >= 5]
            if list(tmpCol.index) and (set(tmpCol.index).issubset(set(catergoryDict["cat1"].keys())) or set(tmpCol.index).issubset(set(catergoryDict["cat2"].keys()))):
                asIso.append(j)

        if len(asIso) < 1: continue
        if len(asIso) == 1:
            if (partialDF.loc[:, asIso[0]].max() - partialDF.loc[:, asIso[0]].min())/float(partialDF.loc[:, asIso[0]].max()) >= 0.75:
                print asIso[0]
        else:
            for item1 in itertools.product(*[catergoryDict["cat1"], catergoryDict["cat2"]]):
                for item2 in list(itertools.combinations(partialDF.columns, 2)):
                    subPartialDF = partialDF.reindex(index=item1, columns=item2)
                    try:
                        if (pd.DataFrame.max(subPartialDF) >= 5).all() and (pd.DataFrame.min(subPartialDF) <= 3).all():
                            partialChi2Results = chi2_contingency(subPartialDF)
                            partialChi2Pvalue = partialChi2Results[1]
                            if partialChi2Pvalue <= 0.0001:
                                print item2
                    except:
                        continue
        # for

def alleleSpecific(refParams=None, tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Analysis the relationship between allele specific and AS events for project {} entry {}...".format(
        tgsSample.projectName, tgsSample.sampleName)
    # prevDir = os.getcwd()
    groupDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName)
    resolveDir(os.path.join(groupDir, "alleleSpecificRelatedAS"))

    collapsedGroupFile = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    readStatFile = "tofu.collapsed.read_stat.txt"
    makeAbundanceFile(collapsedGroupFile, outFile=readStatFile)

    # select_loci_to_phase = os.path.join(cDNA_Cupcake_dir, "phasing", "utils", "select_loci_to_phase.py")
    flncFq = os.path.join(groupDir, "isoseq3", "flnc.fq")
    collapsedGff = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.good.gff")

    seFile = os.path.join(groupDir, "ASE", "PB", "SE.confident.bed12+")
    irFile = os.path.join(groupDir, "ASE", "PB", "IR.confident.bed6+")
    a5ssFile = os.path.join(groupDir, "ASE", "PB", "A5SS.confident.bed6+")
    a3ssFile = os.path.join(groupDir, "ASE", "PB", "A3SS.confident.bed6+")

    seAsPairs = getASpairedIsoforms(seFile, collapsedGroupFile, asType="SE")
    irAsPairs = getASpairedIsoforms(irFile, collapsedGroupFile, asType="IR")
    a5ssAsPairs = getASpairedIsoforms(a5ssFile, collapsedGroupFile, asType="A5SS")
    a3ssAsPairs = getASpairedIsoforms(a3ssFile, collapsedGroupFile, asType="A3SS")

    asPairs = {"SE": seAsPairs, "IR": irAsPairs, "A5SS": a5ssAsPairs, "A3SS": a3ssAsPairs}

    # generate "by_loci" directory
    # removeFiles()
    cmd = "select_loci_to_phase.py {} {} {} {} -c 40 1>/dev/null 2>&1".format(refParams.ref_genome, flncFq, collapsedGff, readStatFile)
    subprocess.call(cmd, shell=True)

    lociDir = glob.glob("by_loci/*size*")
    multiResults = []
    pool = Pool(processes=tgsSample.threads)
    for i in lociDir:
        runPhaser(i, asPairs)
    # for i in lociDir:
    #     singleRunRes = pool.apply_async(runPhaser, (i, refParams, asPairs))
    #     multiResults.append(singleRunRes)
    # for j in multiResults:
    #     j.wait()

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
            sigRelatedDict = relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, isoPairs=isoPairs)
            resultDict[fakeGene].update({asType: sigRelatedDict})
        os.chdir(prevDir)

    '''the structure of resultDict is resultDict = {fakeGene: {asType: {partialCategory: {combinationName: {haplotype1: PB.X.1, haplotype2: PB.X.2, ...}}}}}'''
    partialOut = open("partialAsRelatedHaplotype.txt", "w")
    nopartialOut = open("nopartialAsRelatedHaplotype.txt", "w")
    for g in resultDict:
        for a in resultDict[g]:
            for p in resultDict[g][a]:
                for c in resultDict[g][a][p]:
                    for h in resultDict[g][a][p][c]:
                        if p == "partial":
                            print >> partialOut, "\t".join([g, a, h, resultDict[g][a][p][c][h], c])
                        else:
                            print >> nopartialOut, "\t".join([g, a, h, resultDict[g][a][p][c][h], c])
    partialOut.close()
    nopartialOut.close()
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=seAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=irAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=a5ssAsPairs)
    # relationshipBetweenAlleleSpecifiAndAS(partial=True, nopartial=True, asPairs=a3ssAsPairs)


    print str(datetime.datetime.now()) + " Analysis the relationship between allele specific and AS events for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    poolNum = sum([len(projectsDict[i].keys()) for i in projectsDict])
    # try:
    #     pool = MyPool(processes=poolNum)
    #     multiResults = []
    #     for project in projects:
    #         if project.tgsDataDir not in tgsAnalysisList:
    #             tgsAnalysisList.append(project.tgsDataDir)
    #             singleRunRes = pool.apply_async(alleleSpecific, (refParams, project))
    #             multiResults.append(singleRunRes)
    #     for j in multiResults:
    #         j.wait()
    # except Exception as e:
    #     print e

    project = projects[0]
    alleleSpecific(refParams, project)

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
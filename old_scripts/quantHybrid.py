#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: quantHybrid.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-04 17:33:29
Last modified: 2019-12-04 17:33:30
'''

import datetime, os, argparse
from multiprocessing import Pool

from commonFuncs import *
from commonObjs import *
from Config import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

###############################################
def hisat2Mapping(indexDir, gtfPrefix, project=None, dirSpec=None):
    prevDir = os.getcwd()
    resolveDir(os.path.join(dirSpec.out_dir, project.projectName, "hybrid", "alignment", project.sampleName))
    logDir = os.path.join(dirSpec.out_dir, project.projectName, "hybrid", "log")
    if project.ngsReadPair == "paired":
        leftReadsRepeats = project.ngsLeftReads.split(";")
        rightReadsRepeats = project.ngsRightReads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = ",".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = ",".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -1 {} -2 {} --dta -p {} --max-intronlen 20000 --novel-splicesite-outfile {}.ss " \
              "--un-conc-gz {}.unmapped.fastq.gz 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam ".format(
            indexDir, gtfPrefix, leftReads, rightReads, project.threads, repeatName, repeatName, logDir,
            project.sampleName, repeatName, project.threads, repeatName)
            subprocess.call(cmd, shell=True)

            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, project.threads)
            subprocess.call(cmd, shell=True)
            # pysam.index("{}.sorted.bam".format(n.sampleName), catch_stdout=False)
            cmd = "samtools index -@ {} {}.sorted.bam".format(project.threads, repeatName)
            subprocess.call(cmd, shell=True)
            # bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            # mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            # juncList.append("{}/{}.ss".format(os.getcwd(), repeatName))
            # hisat2mappingLogList.append("{}/{}.{}.hisat2.log".format(logDir, tgsSample.sampleName, repeatName))
            os.chdir("../")
    else:
        if project.leftReads and project.rightReads == None:
            singleReadsRepeats = project.leftReads.split(";")
        elif project.leftReads == None and project.rightReads:
            singleReadsRepeats = project.rightReads.split(";")
        else:
            raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            repeatName = "repeat" + str(i)
            resolveDir(repeatName)
            cmd = "hisat2 -x {}/{} -U {} --dta -p {} --max-intronlen 20000 --novel-splicesite-outfile {}.ss " \
                  "--un-conc {}.unmmaped.fastq 2>{}/{}.{}.hisat2.log | samtools sort -@ {} -o {}.sorted.bam".format(
                indexDir, gtfPrefix, singleReads, project.threads, repeatName, repeatName, logDir,
                project.sampleName, repeatName, project.threads, repeatName
            )
            subprocess.call(cmd, shell=True)
            cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(repeatName, repeatName, project.threads)
            subprocess.call(cmd, shell=True)
            cmd = "samtools index -@ {} {}.sorted.bam".format(project.threads, repeatName)
            subprocess.call(cmd, shell=True)
            # bamList.append("{}/{}.sorted.bam".format(os.getcwd(), repeatName))
            # mergeList.append("{}/{}.gtf".format(os.getcwd(), repeatName))
            # juncList.append("{}/{}.ss".format(os.getcwd(), repeatName))
            # hisat2mappingLogList.append("{}/{}.{}.hisat2.log".format(logDir, tgsSample.sampleName, repeatName))
            os.chdir("../")
    # cmd = "hisat2 -x {}/{} -1 {} -2 {} --dta -p {} --max-intronlen 20000 | samtools sort -@ {} -o {}.sorted.bam".format(
    #     indexDir, gtfPrefix, project.ngsLeftReads, project.ngsRightReads, project.threads, project.threads,
    #     project.sampleName)
    # subprocess.call(cmd, shell=True)
    # cmd = "stringtie {}.sorted.bam -o {}.gtf -p {}".format(project.sampleName, project.sampleName, project.threads)
    # subprocess.call(cmd, shell=True)
    # cmd = "samtools index -@ {} {}.sorted.bam".format(project.threads, project.sampleName)
    # subprocess.call(cmd, shell=True)
    os.chdir(prevDir)

def mergeIsoforms(refParams=None, mergeProjects=None, projects=None, dirSpec=None, projectName=None, refStrain=None):
    tmpDict = {}
    mergedIso2ReadsBed = open("all_sample_merged_iso.bed", "w")
    if not mergeProjects:
        mergeProjects = projects
    for i in mergeProjects:
        unambiBed12 = os.path.join(dirSpec.out_dir, i.projectName, i.sampleName, "filteration", "collapse",
                           "tofu.collapsed.assigned.unambi.bed12+")
        isoBedObj = BedFile(unambiBed12, type="bed12+")
        gene2iso = {}

        for iso in isoBedObj.reads:
            iso.name = "{}_{}".format(i.sampleName, iso.name)
            if isoBedObj.reads[iso].otherList[0] not in gene2iso:
                gene2iso[isoBedObj.reads[iso].otherList[0]] = []
            gene2iso[isoBedObj.reads[iso].otherList[0]].append(isoBedObj.reads[iso])

        for gene in gene2iso:
            for iso in gene2iso[gene]:
                if iso.juncChain not in tmpDict:
                    tmpDict[iso.juncChain] = [iso]
                else:
                    tmpDict[iso.juncChain].append(iso)

    for tmp in tmpDict:
        isos = tmpDict[tmp]
        mergedIsoName = "+".join([x.name for x in isos])
        sortedIsos = sorted(isos, key=lambda x: abs(x.chromStart - x.chromEnd), reverse=True)
        repIso = copy.copy(sortedIsos[0])
        repIso.name = mergedIsoName
        print >> mergedIso2ReadsBed, str(repIso)
    mergedIso2ReadsBed.close()
    return os.path.join(os.getcwd(), "all_sample_merged_iso.bed")

def quantifyHybrid1(refParams=None, mergeProjects=None, projects=None, dirSpec=None, projectName=None, refStrain=None):
    print str(datetime.datetime.now()) + " Quantify merged TGS long reads for project {} with the refStrain is {} with short reads...".format(projectName, refStrain)
    prevDir = os.getcwd()
    hybridDir = os.path.join(dirSpec.out_dir, projectName, "hybrid")
    resolveDir(hybridDir)
    currentThreads = sum([i.threads for i in projects]) if sum([i.threads for i in projects]) < 64 else 64
    mergedIsoBed = mergeIsoforms(projects=projects, dirSpec=dirSpec)
    cmd = "bed2gpe.pl -g 13 {} > all_sample_merged_iso.gpe".format(mergedIsoBed)
    subprocess.call(cmd, shell=True)
    cmd = "genePredToGtf file all_sample_merged_iso.gpe all_sample_merged_iso.gtf"
    subprocess.call(cmd, shell=True)
    cmd = "gffread -g {} -w all_sample_merged_iso.fa all_sample_merged_iso.gtf".format(refParams.ref_genome)
    subprocess.call(cmd, shell=True)
    cmd = "salmon index -t all_sample_merged_iso.fa -i all_sample_merged_iso.salmon_index -k 31 -p {}".format(currentThreads)
    subprocess.call(cmd, shell=True)



def quantifyHybrid(refParams=None, mergeProjects=None, projects=None, dirSpec=None, projectName=None, refStrain=None):
    print str(datetime.datetime.now()) + " Quantify merged TGS long reads for project {} with the refStrain is {} with short reads...".format(projectName, refStrain)
    # remapping to new gtf annotation
    prevDir = os.getcwd()

    hybridDir = os.path.join(dirSpec.out_dir, projectName, "hybrid")
    resolveDir(hybridDir)
    gtfList = []
    if mergeProjects:
        for i in mergeProjects:
            gtf = os.path.join(dirSpec.out_dir, i.projectName, i.sampleName, "filteration", "collapse",
                               "tofu.collapsed.good.gff")
            gtfList.append(gtf)
    else:
        for i in projects:
            gtf = os.path.join(dirSpec.out_dir, i.projectName, i.sampleName, "filteration", "collapse",
                               "tofu.collapsed.good.gff")
            gtfList.append(gtf)
    gtfStr = " ".join(gtfList)
    stringtieGTF = os.path.join(hybridDir, "stringtie.merged.gtf")
    gtfPrefix = "stringtie.merged"
    currentThreads = sum([i.threads for i in projects]) if sum([i.threads for i in projects]) < 64 else 64
    # cmd = "stringtie --merge -G {} -p {} -F 0 -T 0 -f 0 -o {} {}".format(refParams.ref_gtf, currentThreads,
    #                                                                      stringtieGTF, gtfStr)
    cmd = "stringtie --merge -p {} -F 0 -T 0 -f 0 -o {} {}".format(currentThreads, stringtieGTF, gtfStr)
    subprocess.call(cmd, shell=True)
    indexDir = os.path.join(hybridDir, "indexDir")
    resolveDir(indexDir)
    if not checkHisat2IndexExist(indexDir):
        print str(datetime.datetime.now()) + " start hisat2 indexing for new gtf annotation"
        cmd = "hisat2_extract_splice_sites.py {} >{}.ss".format(stringtieGTF, gtfPrefix)
        subprocess.call(cmd, shell=True)
        cmd = "hisat2_extract_exons.py {} >{}.exon".format(stringtieGTF, gtfPrefix)
        subprocess.call(cmd, shell=True)
        cmd = "hisat2-build -p {} --ss {}.ss --exon {}.exon {} {}".format(currentThreads, gtfPrefix, gtfPrefix,
                                                                          refParams.ref_genome, gtfPrefix)
        subprocess.call(cmd, shell=True)
        print str(datetime.datetime.now()) + " end hisat2 indexing for new gtf annotation"
    os.chdir("../")

    currDir = os.getcwd()
    alignDir = os.path.join(hybridDir, "alignment")
    resolveDir(alignDir)
    pool = Pool(processes=len(projects))
    for project in projects:
        pool.apply_async(hisat2Mapping, (indexDir, gtfPrefix, project, dirSpec))
    pool.close()
    pool.join()

    sample2bamDict = {"bamList": [], "bam2sample": {}, "sample2bam": {}}
    # bamList = []
    # bam2sample = {}
    for project in projects:
        projectDir = os.path.join(dirSpec.out_dir, project.projectName, "hybrid", "alignment", project.sampleName)
        if project.ngsReadPair == "paired":
            leftReadsRepeats = project.ngsLeftReads.split(";")
            for i in range(len(leftReadsRepeats)):
                repeatName = "repeat" + str(i)
                bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
                sample2bamDict["bamList"].append(bamFile)
                sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
                sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]
        else:
            if project.leftReads and project.rightReads == None:
                singleReadsRepeats = project.leftReads.split(";")
            elif project.leftReads == None and project.rightReads:
                singleReadsRepeats = project.rightReads.split(";")
            else:
                raise Exception(
                    "The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

            for i in range(len(singleReadsRepeats)):
                repeatName = "repeat" + str(i)
                bamFile = os.path.join(projectDir, repeatName, "{}.sorted.bam".format(repeatName))
                sample2bamDict["bamList"].append(bamFile)
                sample2bamDict["bam2sample"][bamFile] = [project.sampleName + "_" + repeatName]
                sample2bamDict["sample2bam"][project.sampleName + "_" + repeatName] = [bamFile]

    os.chdir(currDir)
    bamStr = " ".join(sample2bamDict["bamList"])

    quantDir = os.path.join(hybridDir, "quant")
    resolveDir(quantDir)
    # quantify gene
    cmd = "featureCounts -a {} -G {} -o featureCounts.gene.txt {} -J -T {}".format(stringtieGTF, refParams.ref_genome,
                                                                                   bamStr, currentThreads)
    subprocess.call(cmd, shell=True)
    # quantify transcript
    cmd = "featureCounts -a {} -G {} -o featureCounts.transcript.txt {} -J -g transcript_id -T {}".format(stringtieGTF,
                                                                                                          refParams.ref_genome,
                                                                                                          bamStr,
                                                                                                          currentThreads)
    subprocess.call(cmd, shell=True)
    geneFeatureCountFile = os.path.join(os.getcwd(), "featureCounts.gene.txt")
    transFeatureCountFile = os.path.join(os.getcwd(), "featureCounts.trans.txt")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Quantify merged TGS long reads with short reads done!"
    return geneFeatureCountFile, transFeatureCountFile, sample2bamDict
    # return expGeneFileDict, expTransFileDict, projectSampleDict

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    gtfDict = {}
    tgsAnalysisList = []
    for project in projects:
        if project.tgsDataDir not in tgsAnalysisList:
            tgsAnalysisList.append(project.tgsDataDir)
            gtf = os.path.join(refParams.out_dir, project.projectName, project.group, "filteration", "collapse", "unambi.gtf")
            if project.projectName not in gtfDict:
                gtfDict[project.projectName] = [gtf]
            else:
                gtfDict[project.projectName].append(gtf)
    quantifyHybrid(gtfDict, refParams=refParams, projects=projects)

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

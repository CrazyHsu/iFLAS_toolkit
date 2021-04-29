#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: processRawReadsAndCorr.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 15:41:52
Last modified: 2019-12-08 15:41:53
'''

import datetime, argparse, shutil, sys
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool

from Config import *
from commonFuncs import *
from commonObjs import NoDaemonProcess, MyPool
from preprocessAndCorrFuncs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################

# def processRawReadsAndCorrection(tgsSample=None, refParams=None, ccsParams=None, mecatParams=None,
#                                  useProovread=None, projectsDict=None, threads=None):
#     tgsProject, ngsProject = tgsSample.tgsProject, tgsSample.ngsProject
#     refParams.threads = threads if threads else refParams.threads
#     prevDir = os.getcwd()
#     workDir = os.path.join(refParams.out_dir, tgsSample.projectName, tgsSample.group)
#     resolveDir(workDir)
#     print str(datetime.datetime.now()) + " Start raw TGS/NGS reads processing for project {} group {}...".format(
#         tgsSample.projectName, tgsSample.group)
#     if tgsProject != None:
#         tgsPlat = tgsProject.tgsPlat
#         if tgsPlat.lower() == "pacbio":
#             inputBam = processReadsFromPacbio(tgsSample=tgsSample, refParams=refParams, ccsParams=ccsParams)
#             # inputBam, subreads = "CCS.bam", "input.subreadset.xml"
#             correctWithDifferentTools(inputBam=inputBam, tgsSample=tgsSample, refParams=refParams,
#                                       useProovread=useProovread, projectsDict=projectsDict, tgsPlat="pacbio")
#             # pass
#         elif tgsPlat.lower() == "nanopore":
#             seqNeedCorrect = processReadsFromNanopore(tgsSample=tgsSample, refParams=refParams)
#             # seqNeedCorrect = os.path.join(tgsSample.tgsDataDir, "nanopore.fq")
#             correctWithDifferentTools(seqNeedCorrect, tgsSample=tgsSample, refParams=refParams, mecatParams=mecatParams,
#                                       useProovread=useProovread, projectsDict=projectsDict, tgsPlat="nanopore")
#         else:
#             raise Exception("Please input the right TGS sequencing strategy, PacBio or Nanopore!")
#
#     if ngsProject != None:
#         # pass
#         processReadsFromNGS(tgsSample=tgsSample, refParams=refParams, projectsDict=projectsDict)
#
#     if tgsProject == None and ngsProject == None:
#         raise Exception("Please input the data you want to be analyzed. "
#                         "By setting the value of 'tgs_dir' or 'ngs_cfg_file' in research sections")
#     print str(datetime.datetime.now()) + " Raw TGS/NGS reads processing for project {} group {} done!".format(
#         tgsSample.projectName, tgsSample.group)
#     os.chdir(prevDir)

def correctionAndAlignToReference(tgsSample=None, refParams=None, mecatParams=None, useProovread=None, dirSpec=None):
    prevDir = os.getcwd()
    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    workDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "rawCorrection")
    resolveDir(workDir)
    resolveDir(logDir, chdir=False)
    print str(datetime.datetime.now()) + " Start raw TGS reads processing for project {} entry {}...".format(
        tgsSample.projectName, tgsSample.sampleName)

    if not os.path.exists(refParams.ref_mm2_index):
        cmd = "minimap2 -d {} -t {} {}".format(refParams.ref_mm2_index, tgsSample.threads, refParams.ref_genome)
        subprocess.call(cmd, shell=True)
    if tgsSample.tgsProcessedData != None:
        tgsPlat = tgsSample.tgsPlat
        if tgsPlat.lower() == "pacbio":
            if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
                samProcess(tgsSample.tgsProcessedData, isBam=True, outPrefix="flnc", toFq=True, toFa=True)
                if useProovread:
                    # hybrid strategy will use proovread to correct
                    flncFq = "flnc.fq"
                    proovreadCorr(seqNeedCorrect=flncFq, tgsSample=tgsSample, dirSpec=dirSpec)
                    proovreadFq = os.path.join(os.getcwd(), "proovread", "proovread.fq")
                    processedFlncFq = proovreadFq

                    cmd = "minimap2 -ax splice:hq -G 100k -uf {} -t {} --secondary=no --MD {} >aln.proovread.sam 2>{}/{}.proovreadFlnc.mm2.log".format(
                        refParams.ref_mm2_index, tgsSample.threads, processedFlncFq, logDir, tgsSample.sampleName)
                    subprocess.call(cmd, shell=True)
                    cmd = "minimap2 -ax splice:hq -G 100k -uf {} -t {} --secondary=no --MD flnc.fq >aln.rawFlnc.sam 2>{}/{}.rawFlnc.mm2.log".format(
                        refParams.ref_mm2_index, tgsSample.threads, logDir, tgsSample.sampleName)
                    subprocess.call(cmd, shell=True)
                    pysam.sort("-o", "aln.proovread.sorted.sam", "--output-fmt", "SAM", "-@", str(refParams.threads),
                               "aln.proovread.sam", catch_stdout=False)
                    alignSam = "aln.proovread.sam"
                else:
                    processedFlncFq = "flnc.fq"
                    cmd = "minimap2 -ax splice:hq -G 100k -uf {} -t {} --secondary=no --MD {} >aln.rawFlnc.sam 2>{}/{}.rawFlnc.mm2.log".format(
                        refParams.ref_mm2_index, tgsSample.threads, processedFlncFq, logDir, tgsSample.sampleName)
                    subprocess.call(cmd, shell=True)
                    pysam.sort("-o", "aln.rawFlnc.sorted.sam", "--output-fmt", "SAM", "-@", str(tgsSample.threads),
                               "aln.rawFlnc.sam", catch_stdout=False)
                    alignSam = "aln.rawFlnc.sam"
                # cmd = "(samtools view -H {};samtools view -f 16 -F 4079 {};samtools view -f 0 -F 4095 {}) | samtools sort -@ {} > flnc.mm2.sorted.bam".format(
                #     flncSam, flncSam, flncSam, refParams.threads)
                # subprocess.call(cmd, shell=True)
                # pybedtools.BedTool("flnc.mm2.sorted.bam").bam_to_bed(bed12=True).saveas("flnc.mm2.sorted.bed12")
            else:
                processedFlncFq = "flnc.fq"
                cmd = "minimap2 -ax splice:hq -G 100k -uf {} -t {} --secondary=no --MD {} >aln.rawFlnc.sam 2>{}/{}.rawFlnc.mm2.log".format(
                    refParams.ref_mm2_index, tgsSample.threads, processedFlncFq, logDir, tgsSample.sampleName)
                subprocess.call(cmd, shell=True)
                pysam.sort("-o", "aln.rawFlnc.sorted.sam", "--output-fmt", "SAM", "-@", str(tgsSample.threads), "aln.rawFlnc.sam",
                           catch_stdout=False)
                alignSam = "aln.rawFlnc.sam"

            # correctWithDifferentTools(inputBam=tgsSample.tgsProcessedData, tgsSample=tgsSample, refParams=refParams,
            #                           useProovread=useProovread, tgsPlat="pacbio")
        elif tgsPlat.lower() == "nanopore":
            cmd = "minimap2 -ax splice -G 100k -k14 -uf {} --secondary=no --MD {} -t {} >aln.rawFlnc.sam 2>{}/{}.rawFlnc.mm2.log".format(
                refParams.ref_mm2_index, tgsSample.tgsProcessedData, tgsSample.threads, logDir, tgsSample.sampleName)
            subprocess.call(cmd, shell=True)
            # pysam.sort("-o", "aln.rawFlnc.sorted.sam", "--output-fmt", "SAM", "-@", str(tgsSample.threads),
            #            "aln.rawFlnc.sam", catch_stdout=False)
            alignSam = "aln.rawFlnc.sam"
            makeLink(tgsSample.tgsProcessedData, "flnc.fa")
            if tgsSample.ngsLeftReads or tgsSample.ngsRightReads:
                fmlrc2Corr(tgsSample=tgsSample, dirSpec=dirSpec)
                tgsSample.tgsProcessedData = os.path.join(os.getcwd(), "fmlrc", "nanopore.fmlrc_corrected.fasta")
                cmd = "minimap2 -ax splice -G 100k -k14 -uf {} --secondary=no --MD {} -t {} >aln.fmlrc_corrected.sam 2>{}/{}.fmlrc_corrected.mm2.log".format(
                    refParams.ref_mm2_index, tgsSample.tgsProcessedData, tgsSample.threads, logDir, tgsSample.sampleName)
                subprocess.call(cmd, shell=True)
                # pysam.sort("-o", "aln.fmlrc_corrected.sorted.sam", "--output-fmt", "SAM", "-@", str(tgsSample.threads),
                #            "aln.fmlrc_corrected.sam", catch_stdout=False)
                cmd = "bamToBed -i <(samtools view -bS aln.fmlrc_corrected.sam) -bed12 > tmp.bed12"
                subprocess.call(cmd, shell=True, executable="/bin/bash")
                cmd = '''filter.pl -o <(awk 'OFS="\t"{if($3-$2>50000){print}}' tmp.bed12) aln.fmlrc_corrected.sam -1 4 > tmp.sam'''
                subprocess.call(cmd, shell=True, executable="/bin/bash")
                cmd = "mv tmp.sam aln.fmlrc_corrected.sam"
                subprocess.call(cmd, shell=True)
                alignSam = "aln.fmlrc_corrected.sam"
                makeLink("fmlrc/nanopore.fmlrc_corrected.fasta", "flnc.fa")
                # mecatCorr(tgsSample.tgsProcessedData, tgsSample=tgsSample, refParams=refParams, mecatParams=mecatParams,
                #           useProovread=useProovread, dirSpec=dirSpec)
            # else:
            #     mecatCorr(tgsSample.tgsProcessedData, tgsSample=tgsSample, refParams=refParams, mecatParams=mecatParams,
            #               dirSpec=dirSpec)
            else:
                cmd = "bamToBed -i <(samtools view -bS aln.rawFlnc.sam) -bed12 > tmp.bed12"
                subprocess.call(cmd, shell=True, executable="/bin/bash")
                cmd = '''filter.pl -o <(awk 'OFS="\t"{if($3-$2>50000){print}}' tmp.bed12) aln.rawFlnc.sam -1 4 > tmp.sam'''
                subprocess.call(cmd, shell=True, executable="/bin/bash")
                cmd = "mv tmp.sam aln.rawFlnc.sam"
                subprocess.call(cmd, shell=True)
        else:
            raise ValueError("The platform you input can't be identified, please input the proper one")

        cmd = "(samtools view -H {};samtools view -f 16 -F 4079 {};samtools view -f 0 -F 4095 {}) | samtools sort -@ {} > flnc.mm2.sorted.bam".format(
            alignSam, alignSam, alignSam, tgsSample.threads)
        subprocess.call(cmd, shell=True)
        pybedtools.BedTool("flnc.mm2.sorted.bam").bam_to_bed(bed12=True).saveas("flnc.mm2.sorted.bed12")
        samProcess("flnc.mm2.sorted.bam", isBam=True, outPrefix="flnc.mm2.sorted", toFa=True)
    
    if tgsSample.tgsProcessedData == None and tgsSample.ngsLeftReads == None or tgsSample.ngsRightReads == None:
        raise Exception("Please input the data you want to be analyzed. "
                        "By setting the value of 'tgs_dir' or 'ngs_cfg_file' in cfg file which the section_type should be 'researchCfg'")
    print str(datetime.datetime.now()) + " Raw TGS reads processing for project {} entry {} done!".format(
        tgsSample.projectName, tgsSample.sampleName)
    os.chdir(prevDir)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg1(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    refParams = defaultCfg.refParams
    mecatParams = defaultCfg.mecatParams
    ccsParams = defaultCfg.ccsParams
    optionTools = defaultCfg.optionTools
    tgsAnalysisList = []
    pybedtools.set_tempdir(refParams.tmp_dir)
    try:
        for project in projects:
            if project.ngsProject:
                makeHisat2Index(refParams=refParams)
                break
        pool = MyPool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(processRawReadsAndCorrection, (
                                                project, refParams, ccsParams, mecatParams, optionTools.use_proovread, projectsDict))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)


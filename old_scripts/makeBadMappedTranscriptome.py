#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: makeBadMappedTranscriptome.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-05-19 16:33:25
Last modified: 2020-05-19 16:33:25
'''

import argparse, datetime
from multiprocessing import Pool
from Bio import SeqIO

from commonFuncs import *
from commonObjs import *
from Config import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

# input  = sys.argv[1] # ex: hq.fasta
# prefix = sys.argv[2] # ex: hq.no5merge
#
# group_file = prefix + '.collapsed.group.txt'
# out_file   = prefix + '.badmapped.fasta'
#
# if not os.path.exists(group_file):
#     print >> sys.stderr, "Collapsed group file {0} does not exist! Abort!".format(group_file)
#
# good = set()
# for line in open(group_file):
#     a,b=line.strip().split()
#     for x in b.split(','): good.add(x)
#
# fa_or_fq = 'fasta' if input.upper().split('.')[-1] in ('FA', 'FASTA') else 'fastq'
# f = open(out_file, 'w')
# print >> sys.stderr, "Reading transcript file {0}....".format(input)
# for r in SeqIO.parse(open(input),fa_or_fq):
#     if r.id not in good:
#         f.write(">{0}\n{1}\n".format(r.id, r.seq))
# f.close()
#
# print >> sys.stderr, "Unmapped or badly mapped transcripts written to: {0}".format(out_file)

def getBadOrUnmappedSeq(inSeq, collapsedGroup, outFile="cogent.badmapped.fasta"):
    good = set()
    with open(collapsedGroup) as f:
        for line in f.readlines():
            a, b = line.strip("\n").split("\t")
            for x in b.split(","): good.add(x)
    faOrFq = "fasta" if inSeq.upper().split('.')[-1] in ("FA", "FASTA", "FSA") else "fastq"
    f = open(outFile, "w")
    for r in SeqIO.parse(open(inSeq), faOrFq):
        if r.id not in good:
            print >> f, ">{}\n{}\n".format(r.id, r.seq)
    f.close()
    # return outFile

def makeBadMappedTranscriptome(refParams=None, paramObj=None):
    print str(datetime.datetime.now()) + " Make unmapped or badly-mapped transcript for project {} group {} done!".format(
        paramObj.projectName, paramObj.group)
    groupDir = os.path.join(refParams.out_dir, paramObj.projectName, paramObj.group)
    cogentDir = os.path.join(groupDir, "cogent")
    resolveDir(cogentDir)

    inSeq = os.path.join(groupDir, "isoseq3", "flnc.fa")
    collapsedGroup = os.path.join(groupDir, "filtration", "collapse", "tofu.collapsed.group.txt")
    cogentInputPrefix = "cogentInput"
    cogentInputFa = cogentInputPrefix + ".fa"
    getBadOrUnmappedSeq(inSeq, collapsedGroup, outFile=cogentInputFa)
    inputNum = len(getFxSequenceId(cogentInputFa))
    k = 30
    if inputNum <= 20000:
        cmd = "run_mash.py -k {} --cpus={} {}".format(k, refParams.threads, cogentInputFa)
        subprocess.call(cmd, shell=True)
        cmd = "process_kmer_to_graph.py -c {}.weights {}.fa {}.fa.s1000k{}.dist cogent_out cogentOut".format(cogentInputPrefix, cogentInputPrefix, cogentInputPrefix, k)
        subprocess.call(cmd, shell=True)
        cmd = "generate_batch_cmd_for_Cogent_reconstruction.py -G {} -S cogentOut cogent_out/ > reconstruction.cmd".format(refParams.ref_mm2)
        subprocess.call(cmd, shell=True)
        cmd = "parallel -j {} < reconstruction.cmd 2>/dev/null".format(refParams.threads)
        subprocess.call(cmd, shell=True)
    else:
        makeLink(cogentInputFa, "isoseq_flnc.fasta")
        cmd = "run_preCluster.py --cpus={}".format(refParams.threads)
        subprocess.call(cmd, shell=True)
        chunkSize = 10
        chunkThreads = refParams.threads / chunkSize
        cmd = "generate_batch_cmd_for_Cogent_family_finding.py --cpus={} --cmd_filename=batchCmd preCluster.cluster_info.csv preCluster_out cogent_out".format(chunkThreads)
        subprocess.call(cmd, shell=True)
        num_lines = sum(1 for line in open("batchCmd"))
        chunkLines = num_lines / chunkSize
        cmd = "split --lines={} batchCmd batchCmd_split -da 4".format(chunkLines)
        subprocess.call(cmd, shell=True)
        cmd = '''ls batchCmd_split* | parallel -j %d "bash {}"''' % chunkThreads
        subprocess.call(cmd, shell=True, executable="/bin/bash")
        cmd = "generate_batch_cmd_for_Cogent_reconstruction.py -G {} -S cogentOut cogent_out/ > reconstruction.cmd".format(refParams.ref_mm2)
        subprocess.call(cmd, shell=True)
    cmd = "cat cogent_out/*/cogent2.renamed.fasta > cogent.fake_genome.fasta"
    subprocess.call(cmd, shell=True)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    refParams = defaultCfg.refParams
    # mecatParams = defaultCfg.mecatParams
    # ccsParams = defaultCfg.ccsParams
    # optionTools = defaultCfg.optionTools
    tgsAnalysisList = []
    try:
        # for project in projects:
        #     if project.ngsProject:
        #         makeHisat2Index(refParams=refParams)
        #         break
        pool = MyPool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(makeBadMappedTranscriptome, (refParams, project))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

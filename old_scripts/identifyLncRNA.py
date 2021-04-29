#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: identifyLncRNA.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-05 17:05:31
Last modified: 2019-12-05 17:05:34
'''

import datetime, argparse
from multiprocessing import Pool
from Config import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################

def lncrnaCPAT(refParams=None, tgsSample=None, lncRNAcpat=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Identify lncRNAs for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    lncRNADir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "lncRNA")
    resolveDir(lncRNADir)
    seqToBeIdentify = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "filtration", "sqanti", "sqanti_artifact_removed.fa")
    if lncRNAcpat[tgsSample.refStrain].model == None or lncRNAcpat[tgsSample.refStrain].hexamer == None:
        coding_seq = lncRNAcpat[tgsSample.refStrain].coding_seq
        noncoding_seq = lncRNAcpat[tgsSample.refStrain].noncoding_seq
        cmd = "make_hexamer_tab.py -c {} -n {} > hexamer.tsv 2>/dev/null".format(coding_seq, noncoding_seq)
        subprocess.call(cmd, shell=True)
        cmd = "make_logitModel.py -c {} -n {} -x hexamer.tsv -o {} 2>/dev/null".format(coding_seq, noncoding_seq, tgsSample.sampleName)
        subprocess.call(cmd, shell=True)
        cmd = "cpat.py -r {} -g {} -d {}.logit.RData -x hexamer.tsv -o lncRNA_CPAT 2>/dev/null".format(refParams.ref_genome, seqToBeIdentify, tgsSample.sampleName)
        subprocess.call(cmd, shell=True)
    else:
        model = lncRNAcpat[tgsSample.refStrain].model
        hexamer = lncRNAcpat[tgsSample.refStrain].hexamer
        cmd = "cpat.py -r {} -g {} -d {} -x {} -o lncRNA_CPAT 2>/dev/null".format(refParams.ref_genome, seqToBeIdentify, model, hexamer)
        subprocess.call(cmd, shell=True)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Identify lncRNAs for project {} group {} done".format(tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    lncRNAcpatParams = defaultCfg.lncRNAcpatParams
    tgsAnalysisList = []
    try:
        pool = Pool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(lncrnaCPAT, (refParams, project, lncRNAcpatParams))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

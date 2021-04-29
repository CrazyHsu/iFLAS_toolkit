#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: fusionTransAnalysis.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-01-02 21:56:54
Last modified: 2020-01-02 21:56:58
'''
from commonFuncs import *
from multiprocessing import Pool

from Config import *
import argparse

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

def fusionTransAnalysis(tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Find fusion transcripts for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    fusionTransDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "fusionTrans")
    logDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "log")
    resolveDir(fusionTransDir)
    inputFa = "../rawCorrection/flnc.fa"
    inputSam = "../rawCorrection/aln.rawFlnc.sorted.sam"
    cmd = "fusion_finder.py --input {} -s {} -o aln.fusionTrans 1>{}/fusion.log 2>&1".format(inputFa, inputSam, logDir)
    subprocess.call(cmd, shell=True)
    print str(datetime.datetime.now()) + " Find fusion transcripts for project {} group {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def main():
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    try:
        pool = Pool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(fusionTransAnalysis, (refParams, project))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main()

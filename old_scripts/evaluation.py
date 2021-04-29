#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: evaluation.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 16:42:10
Last modified: 2019-12-08 16:42:11
'''

import argparse, datetime, random, pybedtools
from multiprocessing import Pool

from Config import *
from commonObjs import *
from evaluationFuncs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################
def evaluation(refParams=None, tgsSample=None, dirSpec=None):
    preCorrEval(tgsSample=tgsSample, dirSpec=dirSpec)
    postCorrEval(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    pybedtools.set_tempdir(refParams.tmp_dir)
    try:
        pool = MyPool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(evaluation, (refParams, project))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

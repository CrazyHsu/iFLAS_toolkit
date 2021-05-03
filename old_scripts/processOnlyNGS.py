#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: processOnlyNGS.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-02-01 19:53:39
Last modified: 2020-02-01 19:53:39
'''

import argparse
import multiprocessing
import multiprocessing.pool
from multiprocessing import Pool

from Config import *
from preprocessAndCorrFuncs import *
from commonObjs import NoDaemonProcess, MyPool

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    refParams = defaultCfg.refParams
    tgsAnalysisList = []
    try:
        pool = MyPool(processes=len(projects))
        multiResults = []
        makeHisat2Index(refParams=refParams)
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                singleRunRes = pool.apply_async(processReadsFromNGS, (project, refParams, projectsDict))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

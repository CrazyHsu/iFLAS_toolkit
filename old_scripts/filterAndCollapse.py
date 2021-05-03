#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: filterAndCollapse.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2019-12-08 16:18:45
Last modified: 2019-12-08 16:18:48
'''

import argparse, pybedtools
from multiprocessing import Pool

from Config import *
from filterAndCollapseFuncs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

##############################################

def filterAndCollpase(refParams=None, tgsSample=None, optionTools=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Filter and collapsing for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    filtrationDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "filtration")
    resolveDir(filtrationDir)
    filterByReferenceInfo(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    # if optionTools.use_sqanti:
    #     sqantiRemoveArtifacts(refParams=refParams, pbSample=pbSample)
    #     collapse_bak(pbSample=pbSample, refParams=refParams, useSqanti=True)
    # else:
    #     collapse_bak(pbSample=pbSample, refParams=refParams, useSqanti=False)
    if optionTools.use_sqanti:
        sqantiRemoveArtifacts(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
        collapse(tgsPlat=tgsSample.tgsPlat, refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    else:
        collapse_no_sqanti(refParams=refParams, tgsSample=tgsSample, dirSpec=dirSpec)
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Filter and collapsing for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
    researchCfg = defaultCfg.researchCfg
    refParams = defaultCfg.refParams
    optionTools = defaultCfg.optionTools
    projectCfg = ProjectCfg(researchCfg.tgs_cfg_file, researchCfg.ngs_cfg_file)
    projects = projectCfg.projects
    projectsDict = projectCfg.projectsDict
    initSysSetting(refParams=refParams, projectsDict=projectsDict)
    tgsAnalysisList = []
    pybedtools.set_tempdir(refParams.tmp_dir)
    try:
        pool = Pool(processes=len(projects))
        multiResults = []
        for project in projects:
            if project.tgsDataDir not in tgsAnalysisList:
                tgsAnalysisList.append(project.tgsDataDir)
                workDir = os.path.join(refParams.out_dir, project.projectName, project.group)
                if project.ngsProject:
                    project.ngsJunctions = os.path.join(workDir, "RNA-seq", "reassembly", "junctions.bed12")
                singleRunRes = pool.apply_async(filterAndCollpase, (refParams, project, optionTools))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)
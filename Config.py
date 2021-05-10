#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: Config.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2019-08-29 19:27:02
Last modified: 2019-08-29 19:27:02
'''
from ConfigParser import ConfigParser, RawConfigParser
from commonFuncs import *

REF_SECTION = "refSection"
# MECAT_PARAMS_SECTION = "mecatParams"
CCS_PARAMS_SECTION = "ccsCalling"
MINIMAP2_SECTION = "minimap2"
COLLAPSE_SECTION = "collapse"
# LNCRNA_CPAT_SECTION = "lncRNAcpat"
DATA_CFG_SECTION = "dataCfg"
OPTION_SECTION = "optionTools"
BUILDIN_SECTION = "buildInTools"
DIR_SECTION = "dirs"

SECTION_TYPE_LIST = [i.lower() for i in REF_SECTION, CCS_PARAMS_SECTION, DATA_CFG_SECTION, OPTION_SECTION,
                                        BUILDIN_SECTION, DIR_SECTION]

SECTIONTYPE = "section_type"
REFSTRAIN = "ref_strain"

REF_TAGS = [REFGENOME, REFSIZE, REFANNOGFF, REFANNOGTF, REFANNOGPE, REFBED, REFMM2INDEX] \
    = ["ref_genome", "ref_size", "ref_gff", "ref_gtf", "ref_gpe", "ref_bed", "ref_mm2_index"]
CCS_TAGS = [CCSMINREADLENGTH, CCSMINREADSCORE, CCSMINSUBREADLENGTH, CCSMINCCSLENGTH, CCSMINPREDICTEDACCURACY,
            CCSMINPASS] \
    = ["min_read_length", "min_read_score", "min_subread_length", "min_ccs_length", "min_predicted_accuracy",
       "min_pass"]
MINIMAP2_TAGS = [MM2INDEX, MAXINTRONLENGTH] \
    = ["mm2_index", "max_intron_length"]
DATA_TAGS = [PROJECTNAME, SAMPLENAME, STRAIN, CONDITION, TGSPLAT, STRATEGY, DATALOTAION, PRIMER, DATAPROCESSEDLOCATION,
             POLYALOCATION, NGSLEFTREADS, NGSRIGHTREADS, NGSREADSPAIRED, NGSREADSLENGTH, NGSJUNCTIONS, USEFMLRC2,
             SINGLERUNTHREADS] \
    = ["project_name", "sample_name", "strain", "condition", "tgs_plat", "strategy", "data_location",
       "data_processed_location", "ngs_left_reads", "ngs_right_reads", "ngs_reads_paired", "ngs_reads_length",
       "ngs_junctions", "use_fmlrc2", "single_run_threads"]
COLLAPSE_TAGS = [MAXIDENTITY, MAXCOVERAGE, FLCOVERAGE, MAXFUZZYJUNCTION, MAX5DIFF, MAX3DIFF, DUNMERGE5SHORTER] \
    = ["max_indentity", "max_coverage", "fl_coverage", "max_fuzzy_junction", "max_5_diff", "max_3_diff", "dun_merge_5_shorter"]
OPTION_TAGS = [THREADS, MEMORY, MERGEDATAFROMSAMESTRAIN, GENE2GO] \
    = ["threads", "memory", "merge_data_from_same_strain", "gene2go"]
BUILDIN_TAGS = [PYTHON, R, BASH] \
    = ["python_location", "r_location", "bash_location"]
DIR_TAGS = [TMPDIR, OUTDIR] \
    = ["tmp_dir", "out_dir"]

REF_TAGS = [SECTIONTYPE, REFSTRAIN] + REF_TAGS
DATA_TAGS = [SECTIONTYPE, REFSTRAIN] + DATA_TAGS
CCS_TAGS = [SECTIONTYPE] + CCS_TAGS
MINIMAP2_TAGS = [SECTIONTYPE] + MINIMAP2_TAGS
COLLAPSE_TAGS = [SECTIONTYPE] + COLLAPSE_TAGS
OPTION_TAGS = [SECTIONTYPE] + OPTION_TAGS
BUILDIN_TAGS = [SECTIONTYPE] + BUILDIN_TAGS
DIR_TAGS = [SECTIONTYPE] + DIR_TAGS
VALID_TAGS = [SECTIONTYPE] + REF_TAGS + CCS_TAGS + DATA_TAGS + OPTION_TAGS + BUILDIN_TAGS + DIR_TAGS

BOOLEAN_TAGS = [USEFMLRC2, DUNMERGE5SHORTER, MERGEDATAFROMSAMESTRAIN]
FLOAT_TAGS = [CCSMINREADSCORE, CCSMINPREDICTEDACCURACY, MAXIDENTITY, MAXCOVERAGE]
INTEGER_TAGS = [CCSMINREADLENGTH, CCSMINSUBREADLENGTH, CCSMINCCSLENGTH, CCSMINPASS, FLCOVERAGE, MAXFUZZYJUNCTION,
                MAX5DIFF, MAX3DIFF, MAXINTRONLENGTH, THREADS, SINGLERUNTHREADS, NGSREADSLENGTH]
STRING_TAGS = [SECTIONTYPE, REFSTRAIN, STRAIN, CONDITION, REFGENOME, REFSIZE, REFANNOGFF, REFANNOGTF, REFBED, REFMM2INDEX, OUTDIR, MEMORY]

FILE_TAGS = [REFGENOME, REFSIZE, REFANNOGFF, REFANNOGTF, REFBED]

TAG_TYPE = {}
for t in BOOLEAN_TAGS: TAG_TYPE[t] = 'boolean'
for t in FLOAT_TAGS: TAG_TYPE[t] = 'float'
for t in INTEGER_TAGS: TAG_TYPE[t] = 'int'
for t in STRING_TAGS: TAG_TYPE[t] = 'string'


# ===================== Functions =====================
# def getSampleCond(sampleCondFile):
#     with open(sampleCondFile) as f:
#         myDict = {}
#         for i in f.readlines():
#             if i.startswith("#"): continue
#             projectName, refStrain, sampleName, condition = i.strip("\n").split("\t")
#             if projectName not in myDict:
#                 myDict[projectName] = {refStrain: {condition: [sampleName]}}
#             elif refStrain not in myDict[projectName]:
#                 myDict[projectName][refStrain] = {condition: [sampleName]}
#             elif condition not in myDict[projectName][refStrain]:
#                 myDict[projectName][refStrain][condition] = [sampleName]
#             else:
#                 myDict[projectName][refStrain][condition].append(sampleName)
#         return myDict


# def getCfg(cfgFile, seqStrategy="tgs"):
#     cfgDict = {}
#     with open(cfgFile) as f:
#         for i in f.readlines():
#             if i.startswith("#"): continue
#             n = TGSsample(i) if seqStrategy == "tgs" else NGSsample(i)
#             if n.projectName not in cfgDict:
#                 cfgDict[n.projectName] = {n.groupName: {n.sampleName: n}}
#             elif n.groupName not in cfgDict[n.projectName]:
#                 cfgDict[n.projectName][n.groupName] = {n.sampleName: n}
#             else:
#                 cfgDict[n.projectName][n.groupName][n.sampleName] = n
#     return cfgDict


# ===================== Classes =====================
class CcsCalling(object):
    def __init__(self, type):
        self.section_type = type
        self.min_read_length = 50
        self.min_read_score = 0.75
        self.min_subread_length = 50
        self.min_ccs_length = 50
        self.min_predicted_accuracy = 0.9
        self.min_pass = 2

    def __setattr__(self, key, value):
        if key != "section_type" and key not in CCS_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, CCS_PARAMS_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])

class Minimap2Section(object):
    def __init__(self, type):
        self.section_type = type
        self.mm2_index = None
        self.max_reads_length = 50000

    def __setattr__(self, key, value):
        if key != "section_type" and key not in MINIMAP2_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, MINIMAP2_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])

class CollapseSection(object):
    def __init__(self, type):
        self.section_type = type
        self.max_identity = None
        self.max_coverage = 50000
        self.max_fuzzy_junction = 5
        self.max_5_diff = 1000
        self.max_3_diff = 500
        self.fl_coverage = 2
        self.dun_merge_5_shorter = True

    def __setattr__(self, key, value):
        if key != "section_type" and key not in COLLAPSE_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, COLLAPSE_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])

class OptionSection(object):
    def __init__(self, type):
        self.section_type = type
        self.threads = 4
        self.memory = "4000M"
        self.merge_data_from_same_strain = True

    def __setattr__(self, key, value):
        if key != "section_type" and key not in OPTION_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, OPTION_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])


class RefSection(object):
    """Encapsulates the meta-information provided for a set of plots."""

    def __init__(self, type):
        self.section_type = type
        self.ref_strain = None
        self.ref_genome = None
        self.ref_size = None
        self.ref_gpe = None
        self.ref_gff = None
        self.ref_gtf = None
        self.ref_bed = None
        self.ref_mm2_index = "ref.mm2"

    def __setattr__(self, key, value):
        if key != "section_type" and key not in REF_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, REF_SECTION))
        self.__dict__[key] = value

    def __eq__(self, o):
        return self.__str__() == o.__str__()

    def __hash__(self):
        return self.__str__().__hash__()

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])


class DataCfg(object):
    def __init__(self, type):
        self.section_type = type
        self.project_name = None
        self.sample_name = None
        self.ref_strain = None
        self.strain = None
        self.condition = None

        self.tgs_plat = None
        self.strategy = None
        self.data_location = None
        self.data_processed_location = None
        self.primer = None
        self.polya_location = None

        self.ngs_left_reads = None
        self.ngs_right_reads = None
        self.ngs_reads_paired = "paired"
        self.ngs_reads_length = None
        self.ngs_junctions = None

        self.use_fmlrc2 = True
        self.single_run_threads = 1

    def __setattr__(self, key, value):
        if key != 'section_type' and key not in DATA_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, DATA_CFG_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])


class BuildInTools(object):
    def __init__(self, type):
        self.section_type = type
        # self.python_location = None
        # self.r_location = None
        self.bash_location = "/bin/bash"

    def __setattr__(self, key, value):
        if key != 'section_type' and key not in BUILDIN_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, BUILDIN_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])


class Dirs(object):
    def __init__(self, type):
        self.section_type = type
        self.out_dir = "outDir"
        self.tmp_dir = None

    def __setattr__(self, key, value):
        if key != 'section_type' and key not in DIR_TAGS:
            raise ValueError("Unrecognized %s tag in section %s" % (key, DIR_SECTION))
        self.__dict__[key] = value

    def __str__(self):
        keys = sorted(self.__dict__.keys())
        return ','.join(['%s=%s' % (x, str(self.__dict__[x])) for x in keys])


class Config(object):
    refInfoParams = {}
    ccsParams = None
    minimap2Params = None
    collapseParams = None
    dataToProcess = []
    optionTools = None
    buildInTools = None
    dirParams = None

    def __init__(self, cfgFile=None, **args):
        self.config = ConfigParser()
        if cfgFile:
            self.config.read(cfgFile)
            self.validate()
            self.instantiate()

    def getIds(self):
        """Returns a list of all plot sections found in the config file."""
        return [s for s in self.config.sections()]

    def getValue(self, section, name, default=None):
        """Returns a value from the configuration."""
        tag = name.lower()

        try:
            value = self.config.get(section, name)
        except Exception:
            return default

        try:
            if tag in FLOAT_TAGS:
                return float(value)
            elif tag in INTEGER_TAGS:
                return int(value)
            elif tag in BOOLEAN_TAGS:
                return value.lower() in ['t', 'true', '0']
            else:  # STRING_TAGS
                return value
        except ValueError as e:
            raise ValueError('Invalid assignment in config file: %s is %s, should be %s\n' % (tag, value, TAG_TYPE[tag]))

    def instantiate(self):
        """Instantiates all objects related to the configuration."""
        for sec in self.getIds():
            sectionType = self.getValue(sec, "section_type")
            if sectionType.lower() not in SECTION_TYPE_LIST:
                raise ValueError(
                    "The section type %s in %s isn't in %s" % (sectionType, sec, ",".join(SECTION_TYPE_LIST)))

            if sectionType.lower() == REF_SECTION.lower():
                s = RefSection(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.refInfoParams[s.ref_strain] = s

            if sectionType.lower() == DATA_CFG_SECTION.lower():
                s = DataCfg(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.dataToProcess.append(s)

            if sectionType.lower() == CCS_PARAMS_SECTION.lower():
                s = CcsCalling(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.ccsParams = s

            if sectionType.lower() == MINIMAP2_SECTION.lower():
                s = Minimap2Section(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.minimap2Params = s

            if sectionType.lower() == COLLAPSE_SECTION.lower():
                s = CollapseSection(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.collapseParams = s

            if sectionType.lower() == OPTION_SECTION.lower():
                s = OptionSection(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.optionTools = s

            if sectionType.lower() == DIR_SECTION.lower():
                s = Dirs(sectionType)
                for o in self.config.options(sec):
                    s.__dict__[o] = self.getValue(sec, o)
                self.dirParams = s

    def validate(self):

        validSet = set(VALID_TAGS)
        for secId in self.getIds():
            configSet = set(self.config.options(secId))
            if "section_type" not in configSet:
                raise ValueError("There must be a 'section_type' option in %s section" % secId)
            badOptions = configSet - validSet
            if badOptions:
                raise ValueError('Unrecognized options found in %s section: %s\nValid options are: %s' % (
                    secId, ', '.join(badOptions), ', '.join(validSet)))

            for o in self.config.options(secId):
                value = self.getValue(secId, o)
                if o in FILE_TAGS:
                    validateFile(value)
                # elif o in DIR_TAGS:
                #     validateDir(value)


# class NGSsample(object):
#     def __init__(self, cfgLine):
#         self.records = cfgLine.strip().split("\t")
#         self.projectName = self.records[0]
#         self.groupName = self.records[1]
#         self.sampleName = self.records[2]
#         self.leftReads = self.records[3]
#         self.rightReads = self.records[4]
#         try:
#             self.ngsReadPair = self.records[5].lower()
#         except IndexError as e:
#             self.ngsReadPair = "paired"
#         self.readLength = self.records[6]
#
#     def __str__(self):
#         return "%s:%s:%s:%s-%s; %s-end sequencing" % (
#             self.projectName, self.groupName, self.sampleName, self.leftReads, self.rightReads, self.ngsReadPair)
#
#     __repr__ = __str__
#
#
# class NGSsample1():
#     def __init__(self, cfgLine):
#         self.records = cfgLine.strip().split("\t")
#         self.projectName = self.records[0]
#         self.strainName = self.records[1]
#         self.sampleName = self.records[2]
#         self.leftReads = self.records[3]
#         self.rightReads = self.records[4]
#         try:
#             self.ngsReadPair = self.records[5].lower()
#         except IndexError as e:
#             self.ngsReadPair = "paired"
#         self.readLength = self.records[6]
#
#     def __str__(self):
#         return "%s:%s:%s:%s-%s; %s-end sequencing" % (
#             self.projectName, self.strainName, self.sampleName, self.leftReads, self.rightReads, self.ngsReadPair)
#
#     __repr__ = __str__
#
#
# class TGSsample(object):
#     def __init__(self, cfgLine):
#         self.records = cfgLine.strip().split("\t")
#         if len(self.records) < 4:
#             raise Exception("You should provide valid TGS config file with at least 4 columns, the format like: "
#                             "project\tgroup\tsampleName\tdataDir. And if you dont't provide the value in the "
#                             "field of 5th, 6th, 7th, we will use the default setting which can broke the "
#                             "pipeline for some unknown reasons.")
#         self.projectName = self.records[0]
#         self.groupName = self.records[1]
#         self.sampleName = self.records[2]
#         self.dataDir = self.records[3]
#         self.primer = self.records[4]
#         self.tgsPlat = self.records[5]
#         self.strategy = self.records[6]
#
#     def __str__(self):
#         return "%s:%s:%s:%s:%s:%s:%s" % (
#             self.projectName, self.groupName, self.sampleName, self.dataDir, self.tgsPlat, self.strategy, self.primer)
#
#     __repr__ = __str__
#
#
# class TGSsample1():
#     def __init__(self, cfgLine):
#         self.records = cfgLine.strip().split("\t")
#         if len(self.records) < 4:
#             raise Exception("You should provide valid TGS config file with at least 4 columns, the format like: "
#                             "project\tgroup\tsampleName\tdataDir. And if you dont't provide the value in the "
#                             "field of 5th, 6th, 7th, we will use the default setting which can broke the "
#                             "pipeline for some unknown reasons.")
#         self.projectName = self.records[0]
#         self.refStrain = self.records[1]
#         self.strainName = self.records[2]
#         self.groupedSampleName = self.records[3]
#         self.dataLocation = self.records[4]
#         self.tgsPlat = self.records[5].lower()
#         self.tgsStrategy = self.records[6]
#         self.primers = self.records[7]
#         if len(self.records) == 9:
#             self.subPrimers = self.records[8]
#         else:
#             self.subPrimers = 1
#         self.primersUsed = None
#
#     def getUniqName(self):
#         self.uniqName = "{}_{}_{}".format(self.projectName, self.strainName, self.groupedSampleName)
#
#     def __str__(self):
#         return "%s:%s:%s:%s:%s:%s:%s:%s" % (
#             self.projectName, self.refStrain, self.strainName, self.groupedSampleName, self.dataLocation, self.tgsPlat,
#             self.tgsStrategy, self.primers)
#
#     __repr__ = __str__
#
#
# class Project(object):
#     def __init__(self):
#         self.projectName = None
#         self.group = None
#         self.sampleName = None
#
#         self.tgsProject = None
#         self.tgsDataDir = None
#         self.tgsPrimer = None
#         self.tgsPlat = None
#         self.tgsStrategy = None
#
#         self.ngsProject = None
#         self.ngsLeftReads = None
#         self.ngsRightReads = None
#         self.ngsSingleReads = None
#         self.ngsReadPair = None
#         self.ngsJunctions = None
#
#
# class ProjectCfg(object):
#     def __init__(self, tgsCfgFile, ngsCfgFile):
#         self.tgsCfgFile = tgsCfgFile
#         self.ngsCfgFile = ngsCfgFile
#         self.projects = self.mergeCfg()
#         self.projectsDict = self.getProjectBound()
#
#     def mergeCfg(self):
#         projects = []
#         if self.tgsCfgFile:
#             validateFile(self.tgsCfgFile)
#             tgsCfgDict = getCfg(self.tgsCfgFile, seqStrategy="tgs")
#             if self.ngsCfgFile:
#                 validateFile(self.ngsCfgFile)
#                 ngsCfgDict = getCfg(self.ngsCfgFile, seqStrategy="ngs")
#             else:
#                 ngsCfgDict = {}
#             for projectName in tgsCfgDict:
#                 for group in tgsCfgDict[projectName]:
#                     for sampleName in tgsCfgDict[projectName][group]:
#                         tgsProject = tgsCfgDict[projectName][group][sampleName]
#                         myProject = Project()
#                         myProject.projectName = projectName
#                         myProject.group = group
#                         myProject.sampleName = sampleName
#
#                         myProject.tgsProject = tgsProject
#                         myProject.tgsDataDir = tgsProject.dataDir
#                         myProject.tgsPrimer = tgsProject.primer
#                         myProject.tgsPlat = tgsProject.tgsPlat
#                         myProject.tgsStrategy = tgsProject.strategy
#                         if projectName in ngsCfgDict and group in ngsCfgDict[projectName] and sampleName in \
#                                 ngsCfgDict[projectName][group]:
#                             ngsProject = ngsCfgDict[projectName][group][sampleName]
#                             myProject.ngsProject = ngsProject
#                             myProject.ngsLeftReads = ngsProject.leftReads
#                             myProject.ngsRightReads = ngsProject.rightReads
#                             myProject.ngsReadPair = ngsProject.ngsReadPair
#                         projects.append(myProject)
#         else:
#             if self.ngsCfgFile:
#                 validateFile(self.ngsCfgFile)
#                 ngsCfgDict = getCfg(self.ngsCfgFile, seqStrategy="ngs")
#             else:
#                 raise Exception("You should write either valid tgs or ngs file path in the config file")
#             for projectName in ngsCfgDict:
#                 for group in ngsCfgDict[projectName]:
#                     for sampleName in ngsCfgDict[projectName][group]:
#                         ngsProject = ngsCfgDict[projectName][group][sampleName]
#                         myProject = Project()
#                         myProject.projectName = projectName
#                         myProject.group = group
#                         myProject.sampleName = sampleName
#
#                         myProject.ngsProject = ngsProject
#                         myProject.ngsLeftReads = ngsProject.leftReads
#                         myProject.ngsRightReads = ngsProject.rightReads
#                         myProject.ngsReadPair = ngsProject.ngsReadPair
#                         projects.append(myProject)
#         return projects
#
#     def getProjectBound(self):
#         projectsDict = {}
#         for project in self.projects:
#             if project.projectName not in projectsDict:
#                 projectsDict[project.projectName] = {
#                     project.group: {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}}}
#             elif project.group not in projectsDict[project.projectName]:
#                 projectsDict[project.projectName].update(
#                     {project.group: {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}}})
#             else:
#                 projectsDict[project.projectName][project.group].update(
#                     {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}})
#         return projectsDict
#
#
# class ProjectCfg1(object):
#     def __init__(self, tgsCfgFile, ngsCfgFile):
#         self.tgsCfgFile = tgsCfgFile
#         self.ngsCfgFile = ngsCfgFile
#         self.tgsDict, self.tgsList, self.tgsStrain2sample, self.tgsPlat2sample, self.tgsRefStrain2sample, self.tgsGroupedSample = self.getTgsCfg()
#         self.ngsDict, self.ngsList = self.getNgsCfg()
#         # self.projects = self.mergeCfg()
#         # self.projectsDict = self.getProjectBound()
#
#     def getTgsCfg(self):
#         validateFile(self.tgsCfgFile)
#         tgsDict = {}
#         tgsList = []
#         tgsStrain2sample = {}
#         tgsPlat2sample = {}
#         tgsGroupedSample = {}
#         tgsRefStrain2sample = {}
#         with open(self.tgsCfgFile) as f:
#             for i in f.readlines():
#                 if i.startswith("#"): continue
#                 t = TGSsample1(i)
#                 if t.tgsPlat == "pacbio":
#                     self.getGroupedPrimersFile(t)
#                 t.getUniqName()
#
#                 if t.projectName not in tgsDict:
#                     tgsDict[t.projectName] = {t.refStrain: {t.strainName: {t.groupedSampleName: t}}}
#                 elif t.refStrain not in tgsDict[t.projectName]:
#                     tgsDict[t.projectName][t.refStrain] = {t.strainName: {t.groupedSampleName: t}}
#                 elif t.strainName not in tgsDict[t.projectName][t.refStrain]:
#                     tgsDict[t.projectName][t.refStrain][t.strainName] = {t.groupedSampleName: t}
#                 else:
#                     tgsDict[t.projectName][t.refStrain][t.strainName][t.groupedSampleName] = t
#
#                 tgsList.append(t)
#
#                 if t.projectName not in tgsPlat2sample:
#                     tgsPlat2sample[t.projectName] = {t.tgsPlat: [t]}
#                 elif t.tgsPlat not in tgsPlat2sample[t.projectName]:
#                     tgsPlat2sample[t.projectName][t.tgsPlat] = [t]
#                 else:
#                     tgsPlat2sample[t.projectName][t.tgsPlat].append(t)
#
#                 if t.projectName not in tgsStrain2sample:
#                     tgsStrain2sample[t.projectName] = {t.strainName: [t]}
#                 elif t.strainName not in tgsStrain2sample[t.projectName]:
#                     tgsStrain2sample[t.projectName][t.strainName] = [t]
#                 else:
#                     tgsStrain2sample[t.projectName][t.strainName].append(t)
#
#                 if t.projectName not in tgsRefStrain2sample:
#                     tgsRefStrain2sample[t.projectName] = {t.refStrain: [t]}
#                 elif t.refStrain not in tgsRefStrain2sample[t.projectName]:
#                     tgsRefStrain2sample[t.projectName][t.refStrain] = [t]
#                 else:
#                     tgsRefStrain2sample[t.projectName][t.refStrain].append(t)
#
#                 if t.projectName not in tgsGroupedSample:
#                     tgsGroupedSample[t.projectName] = {t.groupedSampleName: t}
#                 else:
#                     tgsGroupedSample[t.projectName][t.groupedSampleName] = t
#         return tgsDict, tgsList, tgsStrain2sample, tgsPlat2sample, tgsRefStrain2sample, tgsGroupedSample
#
#     def getNgsCfg(self):
#         validateFile(self.ngsCfgFile)
#         ngsDict = {}
#         ngsList = []
#         with open(self.ngsCfgFile) as f:
#             for i in f.readlines():
#                 if i.startswith("#"): continue
#                 n = NGSsample1(i)
#                 if n.projectName not in ngsDict:
#                     ngsDict[n.projectName] = {n.strainName: {n.sampleName: n}}
#                 elif n.strainName not in ngsDict[n.projectName]:
#                     ngsDict[n.projectName][n.strainName] = {n.sampleName: n}
#                 else:
#                     ngsDict[n.projectName][n.strainName][n.sampleName] = n
#                 ngsList.append(n)
#         return ngsDict, ngsList
#
#     def getGroupedPrimersFile(self, tmpRun):
#         from Bio import SeqIO
#         records = SeqIO.parse(tmpRun.primers, "fasta")
#         primerDir = os.path.dirname(os.path.abspath(tmpRun.primers))
#         if tmpRun.subPrimers == 1:
#             primersUsed = os.path.join(primerDir, tmpRun.projectName + tmpRun.groupedSampleName + ".primers.fa")
#         else:
#             primersUsed = tmpRun.primers
#
#         if not os.path.exists(primersUsed):
#             primersOut = open(primersUsed, "w")
#             recordsDict = SeqIO.to_dict(records)
#             ids = ["5p"]
#             for i in tmpRun.groupedSampleName.split("_"):
#                 ids.append("_".join([tmpRun.strainName, i, "3p"]))
#             resultRecords = [recordsDict[id_] for id_ in ids]
#             SeqIO.write(resultRecords, primersOut, "fasta")
#             primersOut.close()
#         tmpRun.primersUsed = primersUsed
#
#     def mergeCfg(self):
#         projects = []
#         if self.tgsCfgFile:
#             validateFile(self.tgsCfgFile)
#             tgsCfgDict = getCfg(self.tgsCfgFile, seqStrategy="tgs")
#             if self.ngsCfgFile:
#                 validateFile(self.ngsCfgFile)
#                 ngsCfgDict = getCfg(self.ngsCfgFile, seqStrategy="ngs")
#             else:
#                 ngsCfgDict = {}
#             for projectName in tgsCfgDict:
#                 for group in tgsCfgDict[projectName]:
#                     for sampleName in tgsCfgDict[projectName][group]:
#                         tgsProject = tgsCfgDict[projectName][group][sampleName]
#                         myProject = Project()
#                         myProject.projectName = projectName
#                         myProject.group = group
#                         myProject.sampleName = sampleName
#
#                         myProject.tgsProject = tgsProject
#                         myProject.tgsDataDir = tgsProject.dataDir
#                         myProject.tgsPrimer = tgsProject.primer
#                         myProject.tgsPlat = tgsProject.tgsPlat
#                         myProject.tgsStrategy = tgsProject.strategy
#                         if projectName in ngsCfgDict and group in ngsCfgDict[projectName] and sampleName in \
#                                 ngsCfgDict[projectName][group]:
#                             ngsProject = ngsCfgDict[projectName][group][sampleName]
#                             myProject.ngsProject = ngsProject
#                             myProject.ngsLeftReads = ngsProject.leftReads
#                             myProject.ngsRightReads = ngsProject.rightReads
#                             myProject.ngsReadPair = ngsProject.ngsReadPair
#                         projects.append(myProject)
#         else:
#             if self.ngsCfgFile:
#                 validateFile(self.ngsCfgFile)
#                 ngsCfgDict = getCfg(self.ngsCfgFile, seqStrategy="ngs")
#             else:
#                 raise Exception("You should write either valid tgs or ngs file path in the config file")
#             for projectName in ngsCfgDict:
#                 for group in ngsCfgDict[projectName]:
#                     for sampleName in ngsCfgDict[projectName][group]:
#                         ngsProject = ngsCfgDict[projectName][group][sampleName]
#                         myProject = Project()
#                         myProject.projectName = projectName
#                         myProject.group = group
#                         myProject.sampleName = sampleName
#
#                         myProject.ngsProject = ngsProject
#                         myProject.ngsLeftReads = ngsProject.leftReads
#                         myProject.ngsRightReads = ngsProject.rightReads
#                         myProject.ngsReadPair = ngsProject.ngsReadPair
#                         projects.append(myProject)
#         return projects
#
#     def getProjectBound(self):
#         projectsDict = {}
#         for project in self.projects:
#             if project.projectName not in projectsDict:
#                 projectsDict[project.projectName] = {
#                     project.group: {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}}}
#             elif project.group not in projectsDict[project.projectName]:
#                 projectsDict[project.projectName].update(
#                     {project.group: {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}}})
#             else:
#                 projectsDict[project.projectName][project.group].update(
#                     {project.sampleName: {"tgs": project.tgsProject, "ngs": project.ngsProject}})
#         return projectsDict
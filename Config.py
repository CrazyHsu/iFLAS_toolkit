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
                                        BUILDIN_SECTION, DIR_SECTION, MINIMAP2_SECTION, COLLAPSE_SECTION]

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
    = ["project_name", "sample_name", "strain", "condition", "tgs_plat", "strategy", "data_location", "primer",
       "data_processed_location", "polya_location", "ngs_left_reads", "ngs_right_reads", "ngs_reads_paired",
       "ngs_reads_length", "ngs_junctions", "use_fmlrc2", "single_run_threads"]
COLLAPSE_TAGS = [MINIDENTITY, MINCOVERAGE, FLCOVERAGE, MAXFUZZYJUNCTION, MAX5DIFF, MAX3DIFF, DUNMERGE5SHORTER] \
    = ["min_identity", "min_coverage", "fl_coverage", "max_fuzzy_junction", "max_5_diff", "max_3_diff", "dun_merge_5_shorter"]
OPTION_TAGS = [THREADS, MEMORY, GENE2GO, MERGEDATAFROMSAMESTRAIN] \
    = ["threads", "memory", "gene2go", "merge_data_from_same_strain"]
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
VALID_TAGS = [SECTIONTYPE] + REF_TAGS + CCS_TAGS + DATA_TAGS + MINIMAP2_TAGS + COLLAPSE_TAGS + OPTION_TAGS + BUILDIN_TAGS + DIR_TAGS

BOOLEAN_TAGS = [USEFMLRC2, DUNMERGE5SHORTER, MERGEDATAFROMSAMESTRAIN, DUNMERGE5SHORTER]
FLOAT_TAGS = [CCSMINREADSCORE, CCSMINPREDICTEDACCURACY, MINIDENTITY, MINCOVERAGE]
INTEGER_TAGS = [CCSMINREADLENGTH, CCSMINSUBREADLENGTH, CCSMINCCSLENGTH, CCSMINPASS, FLCOVERAGE, MAXFUZZYJUNCTION,
                MAX5DIFF, MAX3DIFF, MAXINTRONLENGTH, THREADS, SINGLERUNTHREADS, NGSREADSLENGTH]
STRING_TAGS = [SECTIONTYPE, REFSTRAIN, STRAIN, CONDITION, REFGENOME, REFSIZE, REFANNOGFF, REFANNOGTF, REFBED, REFMM2INDEX, OUTDIR, MEMORY]

FILE_TAGS = [REFGENOME, REFSIZE, REFANNOGFF, REFANNOGTF, REFBED]
NONE_TAG = "None"

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
        self.max_intron_length = 50000

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
        self.min_identity = 0.9
        self.min_coverage = 0.9
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
                if value.lower() in ['t', 'true']:
                    return True
                elif value.lower() in ['0', 'f', 'false']:
                    return False
            elif tag == NONE_TAG:
                return None
            else:  # STRING_TAGS
                return value
        except ValueError as e:
            raise ValueError('Invalid assignment in config file: %s is %s, should be %s\n' % (tag, value, TAG_TYPE[tag]))

    def instantiate(self):
        """Instantiates all objects related to the configuration."""
        for sec in self.getIds():
            sectionType = self.getValue(sec, "section_type")
            if sectionType.lower() not in SECTION_TYPE_LIST:
                raise ValueError("The section type %s in %s isn't in %s" % (sectionType, sec, ",".join(SECTION_TYPE_LIST)))

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

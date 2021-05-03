#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: readsEnumerate.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-02-24 19:55:42
Last modified: 2021-02-24 19:55:42
'''
from collections import Counter
import itertools
import scipy.stats as stats

readsGroup = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/filtration/collapse/tofu.collapsed.group.txt"
a3ssFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/PB/A3SS.bed6+"
a5ssFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/PB/A5SS.bed6+"
irFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/PB/IR.bed6+"
seFile = "/data/CrazyHsu_data/Projects/isoseq/iflas_20210109_wangbo/iflas_test/test1/B73_all/ASE/PB/SE.bed12+"

embryo_id ="embryo.id.txt"
endo_id = "endo.id.txt"
root_id = "root.id.txt"

eventEnumerate = "all_event_enu.txt_without_pa"

def readsAssign(myFile):
    myDict = {}
    with open(myFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            myDict[infoList[0]] = infoList[1]
    return myDict

def getIsoPairs(asFile, asType="IR"):
    myDict = {}
    with open(asFile) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            if asType not in myDict:
                myDict[asType] = {}
            ase = ":".join([infoList[0], infoList[3]])
            if asType != "SE":
                isos1 = infoList[7].split(",")
                isos2 = infoList[9].split(",")
                myDict[asType][ase] = list(itertools.product(*[isos1, isos2]))
            else:
                isos1 = infoList[15].split(",")
                isos2 = infoList[17].split(",")
                myDict[asType][ase] = list(itertools.product(*[isos1, isos2]))
    return myDict

def calPvalue(listA, listB):
    try:
        g, p, dof, expctd = stats.chi2_contingency([listA, listB])
    except:
        return 1
    return p

def readsEnumerate():
    iso2reads = {}
    with open(readsGroup) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            iso = infoList[0]
            reads = infoList[1].split(",")
            iso2reads[iso] = reads
    embryoDict = readsAssign(embryo_id)
    endoDict = readsAssign(endo_id)
    rootDict = readsAssign(root_id)

    mergeDict = {}
    mergeDict.update(embryoDict)
    mergeDict.update(endoDict)
    mergeDict.update(rootDict)

    with open(eventEnumerate) as f:
        for i in f.readlines():
            infoList = i.strip("\n").split("\t")
            if int(infoList[3]) < 2: continue
            isos = infoList[-1].split(";")
            if isos == [""]: continue
            chrom = infoList[0]
            gene = infoList[2]
            event = infoList[4]
            for j in isos:
                reads = list(itertools.chain.from_iterable([iso2reads[x] for x in j.split(",")]))
                if len(reads) < 5: continue
                readsAssigned = [(i, mergeDict[i]) for i in reads]
                assignedCount = Counter([x[1] for x in readsAssigned])
                embryoCount = assignedCount["embryo"] if "embryo" in assignedCount else 0
                endoCount = assignedCount["endo"] if "endo" in assignedCount else 0
                rootCount = assignedCount["root"] if "root" in assignedCount else 0
                meanCount = sum([embryoCount, endoCount, rootCount])/float(3)
                g, p, dof, expctd = stats.chi2_contingency([[embryoCount, endoCount, rootCount], [meanCount, meanCount, meanCount]])
                category = "all"
                if p <= 0.05:
                    countList = [embryoCount, endoCount, rootCount]
                    countListReverse = sorted(countList, reverse=True)
                    if countListReverse[0] >= countListReverse[1] * 1.5:
                        if countList.index(countListReverse[0]) == 0:
                            category = "embryo"
                        elif countList.index(countListReverse[0]) == 1:
                            category = "endo"
                        else:
                            category = "root"
                    else:
                        if countList.index(countListReverse[0]) == 0 and countList.index(countListReverse[1]) == 1 or \
                            countList.index(countListReverse[0]) == 1 and countList.index(countListReverse[1]) == 0:
                            category = "embryo-endo"
                        elif countList.index(countListReverse[0]) == 0 and countList.index(countListReverse[1]) == 2 or \
                            countList.index(countListReverse[0]) == 2 and countList.index(countListReverse[1]) == 0:
                            category = "embryo-root"
                        else:
                            category = "endo-root"
                print "\t".join(map(str, [gene, chrom, event, j, embryoCount, endoCount, rootCount, p, category]))

readsEnumerate()


#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: indepComb.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-01-03 14:36:49
Last modified: 2020-01-03 14:36:49
'''
import subprocess, argparse
from multiprocessing import Pool
from collections import Counter
from Config import *
from commonObjs import *

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("-c", dest="default_cfg", type=str, default="default.cfg",
                    help="The config file which contain the necessary parameters for reference information "
                         "and mapping settings")
args = parser.parse_args()

def isoCountSimParallele(infoList):
    number = 10000
    eventProbs, isoNum, readCounts = infoList[5], infoList[9], infoList[10]
    events = eventProbs.split(";")
    isoCount = []
    for n in xrange(number):
        combinationDict = {}
        for r in xrange(int(readCounts)):
            eventNames = []
            for e in xrange(len(events)):
                choices = ["choice_{}".format(t) for t in xrange(len(events[e].split(",")))]
                currentChoice = np.random.choice(choices, p=[float(x) for x in events[e].split(",")])
                currentEventName = "event_{}_{}".format(e, currentChoice)
                eventNames.append(currentEventName)
            combination = "-".join(eventNames)
            combinationDict[combination] = ''
        isoCount.append(len(combinationDict))
    isoCounter = Counter(isoCount)
    simCount = isoCounter.most_common(1)[0][0]
    sigSum = 0
    for i in isoCounter:
        if i <= int(isoNum):
            sigSum += isoCounter[i]
    pValue = sigSum / float(number)
    return "\t".join(infoList) + "\t" + str(simCount) + "\t" + str(pValue)

def isoCountSim(fileIn):
    pool = Pool(processes=20)
    with open(fileIn) as f:
        lineList = f.readlines()
        resultList = []
        for line in lineList:
            infoList = line.strip("\n").split("\t")
            if int(infoList[3]) < 3: continue
            tmpRes = pool.apply_async(isoCountSimParallele, (infoList,))
            resultList.append(tmpRes)
        for i in resultList:
            i.wait()
        for i in resultList:
            if i.ready():
                if i.successful():
                    tmpRes = i.get()
                    print tmpRes

def getJuncsFromIsoformFile(isoformFile, isoform2reads):
    juncsDict = {}
    isoformBed = BedFile(isoformFile, type="bed12+")
    allIso = isoformBed.reads
    for r in allIso:
        geneName = allIso[r].otherList[0]
        juncs = allIso[r].parse_intron()
        juncs = sorted(juncs)
        if geneName not in juncsDict:
            juncsDict[geneName] = {}
        for j in juncs:
            if j not in juncsDict[geneName]:
                juncsDict[geneName][j] = {allIso[r].name: isoform2reads[allIso[r].name]}
            else:
                juncsDict[geneName][j].update({allIso[r].name: isoform2reads[allIso[r].name]})
            if "juncsList" not in juncsDict[geneName]:
                juncsDict[geneName] = {"juncsList": set()}
            juncsDict[geneName]["juncsList"].add(j)
    return juncsDict

def juncOneSideMatch(myJunc, allSortedEvents):
    flag = False
    isoList = []
    for ase in allSortedEvents:
        juncLeft, juncRight = ase[0], ase[1]
        if myJunc[0] == juncLeft or myJunc[1] == juncRight:
            isoList.append(ase[-1])
            flag = True
    return flag, list(itertools.chain.from_iterable(isoList))

def getIsosInEvent(junc, allSortedEvents):
    flag = False
    isoList = []
    for ase in allSortedEvents:
        juncLeft, juncRight = ase[0], ase[1]
        if isOverlap((juncLeft, juncRight), junc):
            isoList.append(ase[-1])
            flag = True
    return flag, list(set(itertools.chain.from_iterable(isoList)))

# def juncOneSideMatch1(ase, junc):
#     flag = False
#     isoList = []
#     for ase in allSortedEvents:
#         juncLeft, juncRight = ase[0], ase[1]
#         if myJunc[0] == juncLeft or myJunc[1] == juncRight:
#             isoList.append(ase[-1])
#             flag = True
#     return flag, list(itertools.chain.from_iterable(isoList))

def asEnumerate1(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
    isoform2reads = {}
    with open(isoform2readsFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    # juncsDict = getJuncsFromIsoformFile(isoformFile, isoform2reads)
    juncCombDict = {}
    with open(irFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
            retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            excJuncComb = "{}-{}".format(retentionStart, retentionEnd)
            incJuncComb = "ir_inc:{}-{}".format(retentionStart, retentionEnd)
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                else:
                    juncCombDict[uniqueName][excJuncComb].append(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
                else:
                    juncCombDict[uniqueName][incJuncComb].append(incIsos)
    with open(seFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
            # blockSizes, relStarts = infoList[10], infoList[11]
            juncStart = int(name.split(":")[1].split("@")[0])
            juncEnd = int(name.split(":")[1].split("@")[-1])
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            incJuncComb = name.split(":")[1]
            excJuncComb = "{}-{}".format(juncStart, juncEnd)
            # skipedExonStart, skipedExonEnd, score = int(skipedExonStart), int(skipedExonEnd), float(score)
            # blockSizes = [int(x) for x in blockSizes.split(",")]
            # relStarts = [int(x) for x in relStarts.split(",")]
            # skipedExonStarts, skipedExonEnds = getAbsLoc(skipedExonStart, blockSizes, relStarts)

            incIsos = infoList[15].split(",")
            excIsos = infoList[17].split(",")
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                else:
                    juncCombDict[uniqueName][excJuncComb].append(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
                else:
                    juncCombDict[uniqueName][incJuncComb].append(incIsos)
    with open(a3ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
            incJuncComb = infoList[-2].split(":")[-1]
            excJuncComb = infoList[-1].split(":")[-1]
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                else:
                    juncCombDict[uniqueName][excJuncComb].append(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
                else:
                    juncCombDict[uniqueName][incJuncComb].append(incIsos)
    with open(a5ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
            incJuncComb = infoList[-2].split(":")[-1]
            excJuncComb = infoList[-1].split(":")[-1]
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = {}
                juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
            else:
                if excJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({excJuncComb: [excIsos]})
                else:
                    juncCombDict[uniqueName][excJuncComb].append(excIsos)
                if incJuncComb not in juncCombDict[uniqueName]:
                    juncCombDict[uniqueName].update({incJuncComb: [incIsos]})
                else:
                    juncCombDict[uniqueName][incJuncComb].append(incIsos)
    paDict = {}
    with open(paFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, strand, geneName, paSites = infoList[0:4]
            freqs = infoList[6]
            paSites = [int(x) for x in paSites.split(",")]
            freqs = [float(x) for x in freqs.split(",")]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if strand == "+":
                paDict[uniqueName] = [paSites[0] - 1, paSites[-1], "PA", freqs, paSites]
            else:
                paDict[uniqueName] = [[paSites[0], paSites[-1] + 1, "PA", freqs, paSites]]

    for key in juncCombDict:
        chrom, strand, geneName = key.split(":")
        juncCombs = juncCombDict[key].keys()
        eventsInJuncsDict = {}
        for i in range(len(juncCombs)):
            if "@" in juncCombs[i]:
                juncLeftI = int(juncCombs[i].split(":")[1].split("-")[0])
                juncRightI = int(juncCombs[i].split(":")[1].split("-")[-1])
            else:
                juncLeftI = int(juncCombs[i].split("-")[0])
                juncRightI = int(juncCombs[i].split("-")[1])
            irInc = []
            if "ir_inc" in juncCombs[i]:
                irInc.append(juncCombs[i])
            tmpJuncLeft, tmpJuncRight = juncLeftI, juncRightI
            clusteredJunc = [i]
            for j in range(i+1, len(juncCombs)):
                if "@" in juncCombs[j]:
                    juncLeftJ = int(juncCombs[j].split(":")[1].split("-")[0])
                    juncRightJ = int(juncCombs[j].split(":")[1].split("-")[-1])
                else:
                    juncLeftJ = int(juncCombs[j].split("-")[0])
                    juncRightJ = int(juncCombs[j].split("-")[1])
                if "ir_inc" in juncCombs[j]:
                    irInc.append(juncCombs[j])
                if isOverlap((tmpJuncLeft, tmpJuncRight), (juncLeftJ, juncRightJ)):
                    tmpJuncLeft = tmpJuncLeft if tmpJuncLeft < juncLeftJ else juncLeftJ
                    tmpJuncRight = tmpJuncRight if tmpJuncRight > juncRightJ else juncRightJ
                    clusteredJunc.append(juncCombs[j])
            if bool(eventsInJuncsDict) == False:
                eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
                continue
            # flag = 0
            # for junc in eventsInJuncsDict.keys():
            #     if isOverlap(junc, (tmpJuncLeft, tmpJuncRight)):
            #         if (tmpJuncLeft, tmpJuncRight) not in eventsInJuncsDict[junc]:
            #             eventsInJuncsDict[junc].append(juncCombs[j])
            #         flag = 1
            # if flag == 0:
            #     eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
        eventCount = 0
        combinationDict = {}
        isoformAssign = {}
        for juncComb in eventsInJuncsDict:
            # irTmpJuncLeft, irTmpJuncRight = 0, 0
            irDict = {}
            # irClusteredList = []
            irCount = 0
            eventName = "event_{}".format(eventCount)
            eventCount += 1
            choiceCount = 0
            for tmpJunc in eventsInJuncsDict[juncComb]:
                if "ir_inc" in tmpJunc:
                    tmpJuncLeft = tmpJunc.split(":")[1].split("-")[0]
                    tmpJuncRight = tmpJunc.split(":")[1].split("-")[1]
                    if irCount not in irDict:
                        irDict[irCount] = {"events": [tmpJunc], "junc": [tmpJuncLeft, tmpJuncRight]}
                    else:
                        if tmpJuncLeft < irDict[irCount]["junc"][0] or tmpJuncRight > irDict[irCount]["junc"][1]:
                            irTmpJuncLeft = tmpJuncLeft if tmpJuncLeft < irDict[irCount]["junc"][0] else irDict[irCount]["junc"][0]
                            irTmpJuncRight = tmpJuncRight if tmpJuncRight > irDict[irCount]["junc"][1] else irDict[irCount]["junc"][1]
                            irDict[irCount]["events"].append(tmpJunc)
                            irDict[irCount]["junc"] = [irTmpJuncLeft, irTmpJuncRight]
                        else:
                            irCount += 1
                            irDict[irCount] = {"events": [tmpJunc], "junc": [tmpJuncLeft, tmpJuncRight]}

                    # if irTmpJuncLeft == 0 and irTmpJuncRight == 0:
                    #     irTmpJuncLeft = tmpJunc.split(":")[1].split("-")[0]
                    #     irTmpJuncRight = tmpJunc.split(":")[1].split("-")[1]
                    # else:
                    #     irTmpJuncLeft = irTmpJuncLeft if irTmpJuncLeft < tmpJuncLeft else tmpJuncLeft
                    #     irTmpJuncRight = irTmpJuncRight if irTmpJuncRight > tmpJuncRight else tmpJuncRight
                    # irClusteredList.append(tmpJunc)
                else:
                    choiceName = "choice_{}".format(choiceCount)
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})

                    leftJunc, rightJunc = int(tmpJunc.split("-")[0]), int(tmpJunc.split("-")[1])
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in juncCombDict[key][tmpJunc]])))
                    # tmpIsos.extend(juncCombDict[key][tmpJunc])
                    tmpIsos = list(set(itertools.chain.from_iterable([x for x in juncCombDict[key][tmpJunc]])))
                    combinationDict[eventName][choiceName].update({"isoforms": tmpIsos})
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    combinationDict[eventName][choiceName].update({"junc": (leftJunc, rightJunc)})
                    for iso in tmpIsos:
                        if iso not in isoformAssign:
                            isoformAssign[iso] = {"event": [eventName], "choice": [choiceName]}
                        else:
                            isoformAssign[iso]["event"].append(eventName)
                            isoformAssign[iso]["choice"].append(choiceName)
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    choiceCount += 1

            for ir in irDict:
                leftJunc, rightJunc = irDict[ir]["junc"]
                choiceName = "choice_{}_ir".format(choiceCount)
                tmpIsos = list(set(itertools.chain.from_iterable([juncCombDict[key][x] for x in irDict[ir]["events"]])))
                tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in tmpIsos])))
                combinationDict[eventName][choiceName] = {"reads": tmpReads, "junc": (leftJunc, rightJunc), "isoforms": tmpIsos}
                for iso in tmpIsos:
                    if iso not in isoformAssign:
                        isoformAssign[iso] = {"event": [eventName], "choice": [choiceName]}
                    else:
                        isoformAssign[iso]["event"].append(eventName)
                        isoformAssign[iso]["choice"].append(choiceName)
                tmpComb = "{}_{}".format(eventName, choiceName)
                choiceCount += 1

        paName = "event_pa"
        paCount = 0
        for paSite in paDict[key][-1]:
            paChoice = "choice_pa{}".format(paCount)
            tmpIsos = []
            for iso in isoformAssign:
                if isoformBed[iso].chromEnd == paSite:
                    isoformAssign[iso]["event"].append(paName)
                    isoformAssign[iso]["choice"].append(paChoice)
                    tmpIsos.append(iso)
                    continue
            paCount += 1
            tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in tmpIsos])))
            combinationDict[paName][paChoice] = {"reads": tmpReads, "junc": (paSite-1, paSite), "isoforms": tmpIsos}

        for eventName in combinationDict:
            for choiceName in combinationDict[eventName]:
                print eventName, choiceName

            # irTmpJuncLeft, irTmpJuncRight = irInc[0].split(":")[1].split("-")[0], irInc[0].split(":")[1].split("-")[1]
            # irInJuncsDict = {}
            # for ir1 in range(len(irInc)):
            #     ir1JuncLeft = irInc[ir1].split(":")[1].split("-")[0]
            #     ir1JuncRight = irInc[ir1].split(":")[1].split("-")[1]
            #     clusteredIrJunc = [irInc[ir1]]
            #     for ir2 in range(ir1+1, len(irInc)):
            #         ir2JuncLeft = irInc[ir2].split(":")[1].split("-")[0]
            #         ir2JuncRight = irInc[ir2].split(":")[1].split("-")[1]
            #         if isOverlap((ir1JuncLeft, ir1JuncRight), (ir2JuncLeft, ir2JuncRight)):
            #             irTmpJuncLeft = irTmpJuncLeft if irTmpJuncLeft < ir2JuncLeft else ir2JuncLeft
            #             irTmpJuncRight = irTmpJuncRight if irTmpJuncRight > ir2JuncRight else ir2JuncRight
            #             clusteredIrJunc.append(irInc[ir2])
            #         if bool(irInJuncsDict) == False:
            #             irInJuncsDict[(irTmpJuncLeft, irTmpJuncRight)] = clusteredIrJunc
            #         flag = 0
            #         for junc in irInJuncsDict.keys():
            #             if isOverlap(junc, (irTmpJuncLeft, irTmpJuncRight)):
            #                 if (irTmpJuncLeft, irTmpJuncRight) not in irInJuncsDict[junc]:
            #                     irInJuncsDict[junc].append(irInc[ir2])
            #                 flag = 1
            #         if flag == 0:
            #             irInJuncsDict[(irTmpJuncLeft, irTmpJuncRight)] = clusteredIrJunc


        # juncsInGene = sorted(list(juncsDict[geneName]["juncsList"]))
        # # allSortedEvents = sorted(asDict[key], key=lambda x: (x[0], x[1]))
        # eventsInJuncsDict = {}
        # # tmpJuncLeft, tmpJuncRight = juncsInGene[0]
        # # pas = []
        # for i in range(len(juncsInGene)):
        #     tmpJuncLeft, tmpJuncRight = juncsInGene[i]
        #     clusteredJunc = [(tmpJuncLeft, tmpJuncRight)]
        #     for j in range(i + 1, len(juncsInGene)):
        #         if isOverlap((tmpJuncLeft, tmpJuncRight), juncsInGene[j]):
        #             tmpJuncLeft = tmpJuncLeft if tmpJuncLeft < juncsInGene[j][0] else juncsInGene[j][0]
        #             tmpJuncRight = tmpJuncRight if tmpJuncRight > juncsInGene[j][1] else juncsInGene[j][1]
        #             clusteredJunc.append(juncsInGene[j])
        #     if bool(eventsInJuncsDict) == False:
        #         eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
        #         continue
        #     flag = 0
        #     for junc in eventsInJuncsDict.keys():
        #         if isOverlap(junc, (tmpJuncLeft, tmpJuncRight)):
        #             if (tmpJuncLeft, tmpJuncRight) not in eventsInJuncsDict[junc]:
        #                 eventsInJuncsDict[junc].append((tmpJuncLeft, tmpJuncRight))
        #             flag = 1
        #     if flag == 0:
        #         eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc

        # events = juncCombDict[key].keys()
        # for i in events:
        #     pass




def asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
    isoform2reads = {}
    with open(isoform2readsFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
    asDict = {}
    with open(irFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
            retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            allIsos = list(set(itertools.chain(incIsos, excIsos)))
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in asDict:
                asDict[uniqueName] = [[retentionStart, retentionEnd, "IR", score/1000, allIsos]]
            else:
                asDict[uniqueName].append([retentionStart, retentionEnd, "IR", score/1000, allIsos])
    with open(seFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
            blockSizes, relStarts = infoList[10], infoList[11]
            skipedExonStart, skipedExonEnd, score = int(skipedExonStart), int(skipedExonEnd), float(score)
            blockSizes = [int(x) for x in blockSizes.split(",")]
            relStarts = [int(x) for x in relStarts.split(",")]
            skipedExonStarts, skipedExonEnds = getAbsLoc(skipedExonStart, blockSizes, relStarts)
            juncStart = int(name.split(":")[1].split("@")[0])
            juncEnd = int(name.split(":")[1].split("@")[-1])
            incIsos = infoList[15].split(",")
            excIsos = infoList[17].split(",")
            allIsos = list(set(itertools.chain(incIsos, excIsos)))
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in asDict:
                asDict[uniqueName] = [[juncStart, juncEnd, "SE", (1000-score)/1000, skipedExonStarts, skipedExonEnds, allIsos]]
            else:
                asDict[uniqueName].append([juncStart, juncEnd, "SE", (1000 - score) / 1000, skipedExonStarts, skipedExonEnds, allIsos])
    with open(a3ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
            a3ssStart, a3ssEnd, score = int(a3ssStart), int(a3ssEnd), float(score)
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            allIsos = list(set(itertools.chain(incIsos, excIsos)))
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in asDict:
                asDict[uniqueName] = [[a3ssStart, a3ssEnd, "A3SS", score / 1000, allIsos]]
            else:
                asDict[uniqueName].append([a3ssStart, a3ssEnd, "A3SS", score / 1000, allIsos])
    with open(a5ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
            a5ssStart, a5ssEnd, score = int(a5ssStart), int(a5ssEnd), float(score)
            incIsos = infoList[7].split(",")
            excIsos = infoList[9].split(",")
            allIsos = list(set(itertools.chain(incIsos, excIsos)))
            geneName = name.split(":")[0]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if uniqueName not in asDict:
                asDict[uniqueName] = [[a5ssStart, a5ssEnd, "A5SS", score / 1000, allIsos]]
            else:
                asDict[uniqueName].append([a5ssStart, a5ssEnd, "A5SS", score / 1000, allIsos])
    with open(paFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, strand, geneName, paSites = infoList[0:4]
            freqs = infoList[6]
            paSites = [int(x) for x in paSites.split(",")]
            freqs = [float(x) for x in freqs.split(",")]
            uniqueName = "{}:{}:{}".format(chrom, strand, geneName)
            if strand == "+":
                if uniqueName not in asDict:
                    asDict[uniqueName] = [[paSites[0]-1, paSites[-1], "PA", freqs, paSites]]
                else:
                    asDict[uniqueName].append([paSites[0] - 1, paSites[-1], "PA", freqs, paSites])
            else:
                if uniqueName not in asDict:
                    asDict[uniqueName] = [[paSites[0], paSites[-1]+1, "PA", freqs, paSites]]
                else:
                    asDict[uniqueName].append([paSites[0], paSites[-1]+1, "PA", freqs, paSites])

    juncsDict = getJuncsFromIsoformFile(isoformFile, isoform2reads)

    for key in asDict:
        chrom, strand, geneName = key.split(":")
        juncsInGene = sorted(list(juncsDict[geneName]["juncsList"]))
        allSortedEvents = sorted(asDict[key], key=lambda x: (x[0], x[1]))
        eventsInJuncsDict = {}
        # tmpJuncLeft, tmpJuncRight = juncsInGene[0]
        # pas = []
        for i in range(len(juncsInGene)):
            tmpJuncLeft, tmpJuncRight = juncsInGene[i]
            clusteredJunc = [(tmpJuncLeft, tmpJuncRight)]
            for j in range(i+1, len(juncsInGene)):
                if isOverlap((tmpJuncLeft, tmpJuncRight), juncsInGene[j]):
                    tmpJuncLeft = tmpJuncLeft if tmpJuncLeft < juncsInGene[j][0] else juncsInGene[j][0]
                    tmpJuncRight = tmpJuncRight if tmpJuncRight > juncsInGene[j][1] else juncsInGene[j][1]
                    clusteredJunc.append(juncsInGene[j])
            if bool(eventsInJuncsDict) == False:
                eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc
                continue
            flag = 0
            for junc in eventsInJuncsDict.keys():
                if isOverlap(junc, (tmpJuncLeft, tmpJuncRight)):
                    if (tmpJuncLeft, tmpJuncRight) not in eventsInJuncsDict[junc]:
                        eventsInJuncsDict[junc].append((tmpJuncLeft, tmpJuncRight))
                    flag = 1
            if flag == 0:
                eventsInJuncsDict[(tmpJuncLeft, tmpJuncRight)] = clusteredJunc

        combinationDict = {}
        combinations = []
        eventBoundDict = {}
        isoList = []
        # eventList = []
        # eventFreqsList = []
        # allReadsInJuncList = []
        sortedEvents = sorted(eventsInJuncsDict.keys())
        for i in range(len(sortedEvents)):
            # tmpAllReadsInJunc = []
            eventName = "event_{}".format(i)

            junc = sortedEvents[i]
            overlapFlag, allIsosInEvent = getIsosInEvent(junc, allSortedEvents)
            if overlapFlag:
                count = 0
                tmpIsos = []
                for tmpJunc in eventsInJuncsDict[junc]:
                    count += 1
                    choiceName = "choice_{}".format(count)
                    if eventName not in combinationDict:
                        combinationDict[eventName] = {}
                    eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})
                    isosInJunc = juncsDict[geneName][tmpJunc].keys()
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in isosInJunc])))
                    tmpIsos.extend(isosInJunc)
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    combinations.append(tmpComb)
                leftIsos = list(set(allIsosInEvent) - set(tmpIsos))
                if leftIsos:
                    count += 1
                    choiceName = "choice_{}".format(count)
                    if eventName not in combinationDict:
                        combinationDict[eventName] = {}
                    eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in leftIsos])))
                    tmpIsos.extend(leftIsos)
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    # isoList.extend(ase[-1])
                    combinations.append(tmpComb)

                for ase in allSortedEvents:
                    if ase[2] == "PA": continue
                    if isOverlap((ase[0], ase[1]), junc):
                        count += 1
                        choiceName = "choice_{}".format(count)
                        if eventName not in combinationDict:
                            combinationDict[eventName] = {}
                        eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                        if choiceName not in combinationDict[eventName]:
                            combinationDict[eventName].update({choiceName: {}})
                        tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in ase[-1]])))
                        tmpIsos.extend(ase[-1])
                        combinationDict[eventName][choiceName].update({"reads": tmpReads})
                        tmpComb = "{}_{}".format(eventName, choiceName)
                        # isoList.extend(ase[-1])
                        combinations.append(tmpComb)
                leftIsos = list(set(allIsosInEvent) - set(tmpIsos))
                if leftIsos:
                    count += 1
                    choiceName = "choice_{}".format(count)
                    if eventName not in combinationDict:
                        combinationDict[eventName] = {}
                    eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in leftIsos])))
                    tmpIsos.extend(leftIsos)
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    # isoList.extend(ase[-1])
                    combinations.append(tmpComb)


            count = 0
            for ase in allSortedEvents:
                if ase[2] == "PA": continue
                if isOverlap((ase[0], ase[1]), junc):
                    count += 1
                    choiceName = "choice_{}".format(count)
                    if eventName not in combinationDict:
                        combinationDict[eventName] = {}
                    eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                    if choiceName not in combinationDict[eventName]:
                        combinationDict[eventName].update({choiceName: {}})
                    tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in ase[-1]])))
                    combinationDict[eventName][choiceName].update({"reads": tmpReads})
                    tmpComb = "{}_{}".format(eventName, choiceName)
                    isoList.extend(ase[-1])
                    combinations.append(tmpComb)


            tmpJuncList = eventsInJuncsDict[sortedEvents[i]]
            for r in range(len(tmpJuncList)):
                junc = tmpJuncList[r]
                matchFlag, allIsosInEvent = juncOneSideMatch(junc, allSortedEvents)
                if not matchFlag:
                    continue
                if eventName not in combinationDict:
                    combinationDict[eventName] = {}
                eventBoundDict[eventName] = "-".join(map(str, sortedEvents[i]))
                # tmpJuncList.append(junc)
                isoNames = juncsDict[geneName][junc].keys()
                isoList.extend(isoNames)
                choiceName = "choice_{}".format(r)
                if choiceName not in combinationDict[eventName]:
                    combinationDict[eventName].update({choiceName: {}})
                tmpReads = list(set(itertools.chain.from_iterable([isoform2reads[x] for x in allIsosInEvent])))
                combinationDict[eventName][choiceName].update({"reads": tmpReads})
                tmpComb = "{}_{}".format(eventName, choiceName)
                combinations.append(tmpComb)
                # for r in range(len(isoNames)):
                #     choiceName = "choice_{}".format(r)
                #     if choiceName not in combinationDict[eventName]:
                #         combinationDict[eventName].update({choiceName: {}})
                #     combinationDict[eventName][choiceName].update({"reads": juncsDict[geneName][junc][isoNames[r]]})
                #     tmpComb = "{}_{}".format(eventName, choiceName)
                #     combinations.append(tmpComb)

        isoList = list(set(isoList))
        for eventName in eventBoundDict:
            allReadsInEvent = list(set(itertools.chain.from_iterable([combinationDict[eventName][choice]["reads"] for choice in combinationDict[eventName]])))
            choices = combinationDict[eventName]
            for choice in choices:
                readsCount = len(list(combinationDict[eventName][choice]["reads"]))
                choiceFreq = readsCount / float(len(list(allReadsInEvent)))
                combinationDict[eventName][choice].update({"readsCount": readsCount})
                combinationDict[eventName][choice].update({"choiceFreq": choiceFreq})
            choiceCountStr = ",".join(map(str, [combinationDict[eventName][i]["readsCount"] for i in combinationDict[eventName]]))
            choiceFreqStr = ",".join(map(str, [combinationDict[eventName][i]["choiceFreq"] for i in combinationDict[eventName]]))
            choiceReadsComb = "|".join([",".join(combinationDict[eventName][choice]["reads"]) for choice in combinationDict[eventName]])
            combinationDict[eventName].update({"choiceCountStr": choiceCountStr, "choiceFreqStr": choiceFreqStr, "choiceReadsComb": choiceReadsComb, "allReadsInEvent": len(list(allReadsInEvent))})

        print chrom, strand, geneName, len(combinationDict), ";".join([eventBoundDict[x] for x in combinationDict]), \
            ";".join([combinationDict[x]["choiceFreq"] for x in combinationDict]), str(len(combinations)), \
            ";".join([combinationDict[x]["readsCount"] for x in combinationDict]), \
            ";".join([combinationDict[x]["choiceReadsComb"] for x in combinationDict]), str(len(isoList)), \
            str(len(list(set(itertools.chain([combinationDict[x]["choiceReadsComb"] for x in combinationDict])))))

        # print chrom, strand, geneName, len(combinationDict), \
        #     ";".join([eventBoundDict[x] for x in combinationDict]), \
        #         eventBoundDict[x] for x in combinationDict
            # choiceFreq = [len(combinationDict[eventName][choice]["reads"])/float(len(list(allReadsInEvent))) for choice in combinationDict[eventName]]
            # combinationDict[eventName].update({"choiceFreq": choiceFreq})

        #             tmpReads = juncsDict[geneName][junc][isoNames[r]]
        #             tmpAllReadsInJunc.append(tmpReads)
        #
        #         allReadsInJunc = list(itertools.chain([juncsDict[geneName][junc][r] for r in juncsDict[geneName][junc]]))
        #         eventFreqs = [len(juncsDict[geneName][junc][r])/float(len(allReadsInJunc)) for r in juncsDict[geneName][junc]]
        #         eventFreqsList.append(eventFreqs)
        #         for r in juncsDict[geneName][junc]:
        #             tmpAllReadsInJunc.append(juncsDict[geneName][junc][r])
        #         allReadsInJuncList.append(tmpAllReadsInJunc)
        #     eventList.append(tmpJuncList)
        #
        # paFreqs = [i[3] for i in pas]
        # eventFreqsList.append(paFreqs)

        # print chrom, strand, geneName, len(eventList), ",".join(eventList), ";".join([",".join(i) for i in eventFreqsList]),



def indepComb(tgsSample=None, dirSpec=None):
    print str(datetime.datetime.now()) + " Start analysis independence between alternative splicing isoforms for project {} entry {}...".format(tgsSample.projectName, tgsSample.sampleName)
    prevDir = os.getcwd()
    indepCombDir = os.path.join(dirSpec.out_dir, tgsSample.projectName, tgsSample.sampleName, "indepComb")
    resolveDir(indepCombDir)
    cmd = '''
          asEnumerate.pl -s <(grep -v Novel ../ASE/NGS/SE.PB.bed12+) \
                         -5 <(grep -v Novel ../ASE/NGS/A5SS.PB.bed6+) \
                         -3 <(grep -v Novel ../ASE/NGS/A3SS.PB.bed6+) \
                         -i <(grep -v Novel ../ASE/PB/IR.bed6+) \
                         -p <(awk '{len=split($4,array,",");if(len>1){print}}' ../PA/reads.paGrouped.tsv) \
                         ../PA/deSingleExonRead.3endRevised.bed12+ >enu.tsv
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          asSummary.pl enu.tsv >asSummary.tsv 2>geneAsSummary.tsv
          cut -f1-6 geneAsSummary.tsv | transpose.R | pheatmap.R -header -rower -noClusterR -w=15 -colorL=lightblue -colorH=darkblue -colorN=9 -p=geneAsSummary.pdf 2>/dev/null
          select.pl --index 1,6,3,2,4,5,7- geneAsSummary.tsv >geneAsSummary.reArrange.tsv
    '''
    subprocess.call(cmd, shell=True)
    cmd = '''
          ( head -n1 geneAsSummary.reArrange.tsv
            tail -n+2 geneAsSummary.reArrange.tsv | sort -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n
          ) | cut -f1-6 | transpose.R | pheatmap.R -header -rower -noClusterR -noClusterC -w=15 -colorL=lightblue -colorH=darkblue -colorN=9 -p=geneAsSummary2.pdf 2>/dev/null
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          cut -f7 geneAsSummary.tsv | tail -n +2 | sort | uniq -c | tr -s ' ' | tr ' ' '\t' | select.pl --index 3,2 | sort -k1,1n >geneCountInAsNumber.tsv
          (awk '$1<=5' geneCountInAsNumber.tsv; awk 'BEGIN{sum=0}$1>5{sum+=$2}END{print ">5\t"sum}' geneCountInAsNumber.tsv) | pie.R -p=geneCountInAsNumber.pdf
          cut -f8 geneAsSummary.tsv | tail -n +2 | sort | uniq -c | tr -s ' ' | tr ' ' '\t' | select.pl --index 3,2 | sort -k1,1n >geneCountInIsoformNumber.tsv
          (awk '0<$1&&$1<=5' geneCountInIsoformNumber.tsv; awk 'BEGIN{sum=0}$1>5{sum+=$2}END{print ">5\t"sum}' geneCountInIsoformNumber.tsv) | pie.R -p=geneCountInIsoformNumber.pdf
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = "awk '$4>1 && $12>0 && $13>=5' enu.tsv | tee enu.subset.tsv | isoCountSim.pl -n 10000 -p ge | sum.R -s=1,2,3,4,5,6,7,8 -b=1 -q=0.01,0.5,0.99 >isoCountSim.tsv"
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          sort -k7,7n -k10,10n -k9,9n isoCountSim.tsv | tee bin2d.tsv | \
          perl -e '
              my $i = 1;
              while(<>){
                my @fields = split;
                for my $simuCount(@fields[11..$#fields]){ print "$i\t$simuCount\n"; }
                $i++;
              }' >simuIsoformCount.tsv
          awk '{print NR"\t"$7}' bin2d.tsv >realIsoformCount.tsv
          isoCountSim_bin2d.R <simuIsoformCount.tsv -r=realIsoformCount.tsv -w=15 -y='Isoform Number' -x1=0 -x2=$(wc -l realIsoformCount.tsv) -p=isoCountSim_bin2d.pdf
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          cut -f9,10 enu.subset.tsv | multiNomialTest.R -m -n=1000000 >EMT.log 2>EMT.tsv
          paste <(select.pl --index 7,9-11,8 isoCountSim.tsv) <(grep -v There EMT.tsv) | sort -k1,1n -k3,3n -k2,2n >range.tsv
          isoCountSim.R <range.tsv -w=50 -size=0.5 -barWidth=1 -x1=0 -x2=$(wc -l range.tsv) -x=Genes -y='Isoform Number' -p=isoCountSim.pdf
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          sort -k5,5n range.tsv | awk 'BEGIN{OFS="\t"}{print NR,$5,"MC p-value"; print NR,$6,"Distance"; print NR,$7,"Spearman R"; print NR,$8,"Spearman p-value"; print NR,$9,"EMT p-value"; print NR,$10,"EMT corrected p-value"}' | point.R -colorV=V3 -alpha=0.8 -w=15 -ho=0.05 -x='Gene Index' -y=Value -p=freq.point.pdf 2>/dev/null
          awk 'BEGIN{OFS="\t"}{print "MC p-value",$5; print "Distance",$6; print "Spearman R",$7; print "Spearman p-value",$8; print "EMT p-value",$9; print "EMT corrected p-value",$10}' range.tsv | box.R -fillV=V1 -nJ -w=12 -ho=0.05 -noGuide -x=Metrics -y=Value -p=freq.box.pdf 2>/dev/null
          awk 'BEGIN{OFS="\t"}{print $5,"MC p-value"; print $6,"Distance"; print $7, "Spearman R"; print $8,"Spearman p-value"; print $9,"EMT p-value";print $10,"EMT corrected p-value"}' range.tsv | distrCurve.R -colorV=V2 -w=15 -v=0.05 -x=Value -y=Density -d -p=freq.distriCurve.pdf 2>/dev/null
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    cmd = '''
          (
            echo -e "#Chr\tStrand\tGene\tASE Count\tASE ID\tASE (Inclusion) Ratio\tIsoform Frequency (RNA-seq)\tIsoform Coverage (Iso-seq)\tReal Isoform Count\tTotal Iso-seq Coverage\tMonte Carlo p-value\tSpearman R\tSpearman p-value\tEMT p-value\tEMT corrected p-value\tASE Association"
            paste <(cut -f1-6,9,10,12,13 enu.subset.tsv) <(cut -f8 isoCountSim.tsv) <(grep -v There EMT.tsv|cut -f2-) | awk '{if($11<0.05&&$15<0.05){type="Dependent"}else{type="Independent"};print $0"\t"type}'
          ) >statistic.indepComb.tsv
    '''
    subprocess.call(cmd, shell=True, executable="/bin/bash")
    os.chdir(prevDir)
    print str(datetime.datetime.now()) + " Analysis independence between alternative splicing isoforms for project {} entry {} done!".format(tgsSample.projectName, tgsSample.sampleName)

def main(defaultCfg=None):
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
                singleRunRes = pool.apply_async(indepComb, (refParams, project))
                multiResults.append(singleRunRes)
        for j in multiResults:
            j.wait()
    except Exception as e:
        print e

if __name__ == '__main__':
    defaultCfg = Config1(args.default_cfg)
    main(defaultCfg=defaultCfg)

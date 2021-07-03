from commonFuncs import *
from commonObjs import *
from collections import Counter
import itertools

# def asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
#     isoform2reads = {}
#     with open(isoform2readsFile) as f:
#         for line in f.readlines():
#             isoform, reads = line.strip("\n").split("\t")
#             isoform2reads[isoform] = reads.split(",")
#     isoformBed = BedFile(isoformFile, type="bed12+").reads
#     gene2isoDict = {}
#     for isoName in isoformBed:
#         geneName = isoformBed[isoName].otherList[0]
#         if geneName not in gene2isoDict:
#             gene2isoDict[geneName] = [isoformBed[isoName]]
#         else:
#             gene2isoDict[geneName].append(isoformBed[isoName])
#     juncCombDict = {}
#     with open(irFile) as f:
#         for line in f.readlines():
#             if line.startswith("#"): continue
#             infoList = line.strip("\n").split("\t")
#             chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
#             retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
#             geneName = name.split(":")[0]
#             uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
#             if uniqueName not in juncCombDict:
#                 juncCombDict[uniqueName] = []
#             juncCombDict[uniqueName].append((retentionStart, retentionEnd))
#
#     with open(seFile) as f:
#         for line in f.readlines():
#             if line.startswith("#"): continue
#             infoList = line.strip("\n").split("\t")
#             chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
#             juncStart = int(name.split(":")[-1].split("@")[0])
#             juncEnd = int(name.split(":")[-1].split("@")[-1])
#             geneName = name.split(":")[0]
#             uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
#             if uniqueName not in juncCombDict:
#                 juncCombDict[uniqueName] = []
#             juncCombDict[uniqueName].append((juncStart, juncEnd))
#
#     with open(a3ssFile) as f:
#         for line in f.readlines():
#             if line.startswith("#"): continue
#             infoList = line.strip("\n").split("\t")
#             chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
#             excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
#             excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
#             geneName = name.split(":")[0]
#             uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
#             if uniqueName not in juncCombDict:
#                 juncCombDict[uniqueName] = []
#             juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))
#
#     with open(a5ssFile) as f:
#         for line in f.readlines():
#             if line.startswith("#"): continue
#             infoList = line.strip("\n").split("\t")
#             chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
#             excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
#             excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
#             geneName = name.split(":")[0]
#             uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
#             if uniqueName not in juncCombDict:
#                 juncCombDict[uniqueName] = []
#             juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))
#
#     paDict = {}
#     if paFile:
#         with open(paFile) as f:
#             for line in f.readlines():
#                 if line.startswith("#"): continue
#                 infoList = line.strip("\n").split("\t")
#                 chrom, strand, geneName, paSites = infoList[0:4]
#                 paSites = [int(x) for x in paSites.split(",")]
#                 uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
#                 if uniqueName not in juncCombDict:
#                     juncCombDict[uniqueName] = []
#                 if geneName in gene2isoDict:
#                     if strand == "+":
#                         isoEnds = [x.chromEnd for x in gene2isoDict[geneName]]
#                         pa = (min(isoEnds), max(isoEnds))
#                     else:
#                         isoStarts = [x.chromStart for x in gene2isoDict[geneName]]
#                         pa = (min(isoStarts), max(isoStarts))
#                 else:
#                     if strand == "+":
#                         pa = (paSites[0] - 1 - 20, paSites[-1] + 20)
#                     else:
#                         pa = (paSites[0] - 20, paSites[-1] + 1 + 20)
#                 juncCombDict[uniqueName].append(pa)
#                 paDict[geneName] = [pa, paSites]
#
#     for uniqueName in juncCombDict:
#         chrom, strand, geneName = uniqueName.split("@")
#         if geneName not in gene2isoDict: continue
#         isosInGene = gene2isoDict[geneName]
#         juncs = sorted(juncCombDict[uniqueName])
#         tmpJuncStart, tmpJuncEnd = juncs[0]
#         newJuncs = []
#         for junc in juncs[1:]:
#             if isOverlap((tmpJuncStart, tmpJuncEnd), junc):
#                 tmpJuncStart = tmpJuncStart if tmpJuncStart < junc[0] else junc[0]
#                 tmpJuncEnd = tmpJuncEnd if tmpJuncEnd > junc[1] else junc[1]
#             else:
#                 newJuncs.append((tmpJuncStart, tmpJuncEnd))
#                 tmpJuncStart = junc[0]
#                 tmpJuncEnd = junc[1]
#         newJuncs.append((tmpJuncStart, tmpJuncEnd))
#         newJuncs = sorted(newJuncs)
#         # print newJuncs
#         # print [x.name for x in isosInGene]
#         eventCombination = {}
#         isoformAssign = {}
#         for junc in newJuncs:
#             eventName = "event_{}-{}".format(junc[0], junc[1])
#             eventCombination[eventName] = {}
#             for iso in isosInGene:
#                 if not isOverlap(junc, (iso.chromStart, iso.chromEnd)):
#                     continue
#                 juncsInIso = iso.introns
#                 if junc in juncsInIso:
#                     choiceName = "choice_{}".format("-".join(map(str, junc)))
#                     if choiceName not in eventCombination[eventName]:
#                         eventCombination[eventName][choiceName] = [iso]
#                     else:
#                         eventCombination[eventName][choiceName].append(iso)
#                 else:
#                     overlapedJuncs = []
#                     for tmpJunc in juncsInIso:
#                         if isOverlap(tmpJunc, junc):
#                             overlapedJuncs.append(tmpJunc)
#                     if overlapedJuncs:
#                         choiceName = "choice_{}".format(";".join(["-".join(map(str, x)) for x in overlapedJuncs]))
#                     else:
#                         if geneName in paDict and isOverlap(junc, paDict[geneName][0]):
#                             if strand == "+":
#                                 dist2pa = [abs(x - iso.chromEnd) for x in paDict[geneName][1]]
#                             else:
#                                 dist2pa = [abs(x - iso.chromStart) for x in paDict[geneName][1]]
#                             index = dist2pa.index(min(dist2pa))
#                             choiceName = "choice_pa{}".format(index)
#                         else:
#                             choiceName = "choice_ir"
#                     if choiceName not in eventCombination[eventName]:
#                         eventCombination[eventName][choiceName] = [iso]
#                     else:
#                         eventCombination[eventName][choiceName].append(iso)
#                 if iso.name not in isoformAssign:
#                     isoformAssign[iso.name] = {eventName: choiceName}
#                 else:
#                     isoformAssign[iso.name].update({eventName: choiceName})
#         # print isoformAssign
#         # for eventName in eventCombination:
#         #     for choiceName in eventCombination[eventName]:
#         #         print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]]
#         resultDict = {}
#         allEventChoice = []
#         choiceFreqList = []
#         allReads = []
#         for eventName in eventCombination:
#             if eventName not in resultDict:
#                 resultDict[eventName] = {}
#             allChoices = eventCombination[eventName]
#             allIsosInEvent = list(set(itertools.chain.from_iterable(allChoices.values())))
#             allReadsInEvent = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allIsosInEvent])))
#             allReads.extend(allReadsInEvent)
#             tmpList = []
#             for choiceName in allChoices:
#                 tmpList.append((eventName, choiceName))
#                 readsInChoice = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allChoices[choiceName]])))
#                 choiceFreq = len(readsInChoice)/float(len(allReadsInEvent))
#                 resultDict[eventName].update({choiceName: {"readsInChoice": readsInChoice, "choiceFreq": choiceFreq}})
#                 # print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]], len(readsInChoice), len(allReadsInEvent)
#             allEventChoice.append(tmpList)
#             choiceFreqList.append(",".join(map(str, [resultDict[eventName][x]["choiceFreq"] for x in resultDict[eventName]])))
#         combinations = list(itertools.product(*allEventChoice))
#         # print combinations, len(combinations)
#         # print choiceFreqList
#         combination2iso = {}
#         identifiedComb = []
#         for comb in combinations:
#             combination2iso[comb] = {"iso": [], "freq": "NA", "freqList": [], "reads": []}
#             for iso in isoformAssign:
#                 if sorted(isoformAssign[iso].items()) == sorted(comb):
#                     identifiedComb.append(comb)
#                     combination2iso[comb]["iso"].append(iso)
#                     freqList = []
#                     for i in comb:
#                         freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
#                     combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
#                     combination2iso[comb].update({"freqList": freqList})
#                 else:
#                     freqList = []
#                     for i in comb:
#                         freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
#                     combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
#                     combination2iso[comb].update({"freqList": freqList})
#             readsInComb = set(itertools.chain.from_iterable([isoform2reads[x] for x in combination2iso[comb]["iso"]]))
#             counts = len(readsInComb)
#             combination2iso[comb].update({"counts": counts})
#             combination2iso[comb].update({"reads": readsInComb})
#         # print combination2iso.values()
#         identifiedComb = list(set(identifiedComb))
#         combCounts = sorted([(combination2iso[x]["iso"], combination2iso[x]["counts"]) for x in identifiedComb], key=lambda x: x[1], reverse=True)
#         majorCombCount = 0
#         tmpSum = 0
#         isoList = []
#         for x in combCounts:
#             majorCombCount += 1
#             tmpSum += x
#             isoList.append(x[0])
#             if tmpSum/float(sum([y[1] for y in combCounts])) >= 0.75:
#                 break
#         print "\t".join(map(str, [chrom, strand, geneName, len(eventCombination), ";".join(eventCombination.keys()), \
#             ";".join(choiceFreqList), len(combinations), \
#             ",".join(map(str, [combination2iso[x]["counts"] for x in identifiedComb])), \
#             ",".join(map(str, [combination2iso[x]["freq"] for x in identifiedComb])), \
#             len(list(set(identifiedComb))), len(list(set(allReads))), majorCombCount, ";".join([",".join(x) for x in isoList])]))

def scoreAsIsoform(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, collapsedTrans2reads, annoIsoformFile):
    # isoform2reads = {}
    isoform2reads = getDictFromFile(collapsedTrans2reads, sep="\t", inlineSep=",", valueCol=2)
    annoIsoformDict = getDictFromFile(annoIsoformFile, sep="\t", keyCol=4)
    # with open(isoform2readsFile) as f:
    #     for line in f.readlines():
    #         isoform, reads = line.strip("\n").split("\t")
    #         isoform2reads[isoform] = reads.split(",")
    isoformBed = BedFile(isoformFile, type="bed12+").reads
    gene2isoDict = {}
    for isoName in isoformBed:
        geneName = isoformBed[isoName].otherList[0]
        if geneName not in gene2isoDict:
            gene2isoDict[geneName] = [isoformBed[isoName]]
        else:
            gene2isoDict[geneName].append(isoformBed[isoName])
    juncCombDict = {}
    with open(irFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, retentionStart, retentionEnd, name, score, strand = infoList[0:6]
            retentionStart, retentionEnd, score = int(retentionStart), int(retentionEnd), float(score)
            geneNameList = [".".join(x.split(".")[0:2]) for x in infoList[7].split(",") + infoList[9].split(",")]
            geneName = Counter(geneNameList).most_common()[0][0]
            # geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((retentionStart, retentionEnd))

    with open(seFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, skipedExonStart, skipedExonEnd, name, score, strand = infoList[0:6]
            juncStart = int(name.split(":")[-1].split("@")[0])
            juncEnd = int(name.split(":")[-1].split("@")[-1])
            geneNameList = [".".join(x.split(".")[0:2]) for x in infoList[16].split(",") + infoList[18].split(",")]
            geneName = Counter(geneNameList).most_common()[0][0]
            # geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((juncStart, juncEnd))

    with open(a3ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a3ssStart, a3ssEnd, name, score, strand = infoList[0:6]
            excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
            excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
            geneNameList = [".".join(x.split(".")[0:2]) for x in infoList[7].split(",") + infoList[9].split(",")]
            geneName = Counter(geneNameList).most_common()[0][0]
            # geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))

    with open(a5ssFile) as f:
        for line in f.readlines():
            if line.startswith("#"): continue
            infoList = line.strip("\n").split("\t")
            chrom, a5ssStart, a5ssEnd, name, score, strand = infoList[0:6]
            excJuncStart = int(infoList[-1].split(":")[-1].split("-")[0])
            excJuncEnd = int(infoList[-1].split(":")[-1].split("-")[1])
            geneNameList = [".".join(x.split(".")[0:2]) for x in infoList[7].split(",") + infoList[9].split(",")]
            geneName = Counter(geneNameList).most_common()[0][0]
            # geneName = name.split(":")[0]
            uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
            if uniqueName not in juncCombDict:
                juncCombDict[uniqueName] = []
            juncCombDict[uniqueName].append((excJuncStart, excJuncEnd))

    paDict = {}
    if paFile:
        with open(paFile) as f:
            for line in f.readlines():
                if line.startswith("#"): continue
                infoList = line.strip("\n").split("\t")
                chrom, strand, geneName, paSites = infoList[0:4]
                paSites = [int(x) for x in paSites.split(",")]
                uniqueName = "{}@{}@{}".format(chrom, strand, geneName)
                if uniqueName not in juncCombDict:
                    juncCombDict[uniqueName] = []
                if geneName in gene2isoDict:
                    if strand == "+":
                        isoEnds = [x.chromEnd for x in gene2isoDict[geneName]]
                        pa = (min(isoEnds), max(isoEnds))
                    else:
                        isoStarts = [x.chromStart for x in gene2isoDict[geneName]]
                        pa = (min(isoStarts), max(isoStarts))
                else:
                    if strand == "+":
                        pa = (paSites[0] - 1 - 20, paSites[-1] + 20)
                    else:
                        pa = (paSites[0] - 20, paSites[-1] + 1 + 20)
                juncCombDict[uniqueName].append(pa)
                paDict[geneName] = [pa, paSites]

    asEnumerateOut= open("asEnumerate.txt", "w")
    isoformScoreOut= open("isoformScore.txt", "w")
    for uniqueName in juncCombDict:
        chrom, strand, geneName = uniqueName.split("@")
        if geneName not in gene2isoDict: continue
        isosInGene = gene2isoDict[geneName]
        juncs = sorted(juncCombDict[uniqueName])
        tmpJuncStart, tmpJuncEnd = juncs[0]
        newJuncs = []
        for junc in juncs[1:]:
            if isOverlap((tmpJuncStart, tmpJuncEnd), junc):
                tmpJuncStart = tmpJuncStart if tmpJuncStart < junc[0] else junc[0]
                tmpJuncEnd = tmpJuncEnd if tmpJuncEnd > junc[1] else junc[1]
            else:
                newJuncs.append((tmpJuncStart, tmpJuncEnd))
                tmpJuncStart = junc[0]
                tmpJuncEnd = junc[1]
        newJuncs.append((tmpJuncStart, tmpJuncEnd))
        newJuncs = sorted(newJuncs)
        eventCombination = {}
        isoformAssign = {}
        for junc in newJuncs:
            eventName = "event_{}-{}".format(junc[0], junc[1])
            eventCombination[eventName] = {}
            for iso in isosInGene:
                if not isOverlap(junc, (iso.chromStart, iso.chromEnd)):
                    continue
                juncsInIso = iso.introns
                if junc in juncsInIso:
                    choiceName = "choice_{}".format("-".join(map(str, junc)))
                    if choiceName not in eventCombination[eventName]:
                        eventCombination[eventName][choiceName] = [iso]
                    else:
                        eventCombination[eventName][choiceName].append(iso)
                else:
                    overlapedJuncs = []
                    for tmpJunc in juncsInIso:
                        if isOverlap(tmpJunc, junc):
                            overlapedJuncs.append(tmpJunc)
                    if overlapedJuncs:
                        choiceName = "choice_{}".format(";".join(["-".join(map(str, x)) for x in overlapedJuncs]))
                    else:
                        if geneName in paDict and isOverlap(junc, paDict[geneName][0]):
                            if strand == "+":
                                dist2pa = [abs(x - iso.chromEnd) for x in paDict[geneName][1]]
                            else:
                                dist2pa = [abs(x - iso.chromStart) for x in paDict[geneName][1]]
                            index = dist2pa.index(min(dist2pa))
                            choiceName = "choice_pa{}".format(index)
                        else:
                            if not getBlockLength(getOverlapOfTuple(isoformBed[iso.name].exons, [junc])) == junc[1] - junc[0]: continue
                            choiceName = "choice_ir"
                    if choiceName not in eventCombination[eventName]:
                        eventCombination[eventName][choiceName] = [iso]
                    else:
                        eventCombination[eventName][choiceName].append(iso)
                if iso.name not in isoformAssign:
                    isoformAssign[iso.name] = {eventName: choiceName}
                else:
                    isoformAssign[iso.name].update({eventName: choiceName})
        resultDict = {}
        allEventChoice = []
        choiceFreqList = []
        allReads = []
        for eventName in eventCombination:
            if eventName not in resultDict:
                resultDict[eventName] = {}
            allChoices = eventCombination[eventName]
            allIsosInEvent = list(set(itertools.chain.from_iterable(allChoices.values())))
            allReadsInEvent = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allIsosInEvent])))
            allReads.extend(allReadsInEvent)
            tmpList = []
            for choiceName in allChoices:
                tmpList.append((eventName, choiceName))
                readsInChoice = list(
                    set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allChoices[choiceName]])))
                choiceFreq = len(readsInChoice) / float(len(allReadsInEvent))
                resultDict[eventName].update({choiceName: {"readsInChoice": readsInChoice, "choiceFreq": choiceFreq}})
            allEventChoice.append(tmpList)
            choiceFreqList.append(",".join(map(str, [resultDict[eventName][x]["choiceFreq"] for x in resultDict[eventName]])))

        for n in range(1, len(allEventChoice) + 1):
            tmpComb = list(itertools.combinations(allEventChoice, n))
            for m in range(len(tmpComb)):
                combinations = list(itertools.product(*(tmpComb[m])))
                combination2iso = {}
                identifiedComb = []
                for comb in combinations:
                    combination2iso[comb] = {"iso": [], "freq": "NA", "freqList": [], "reads": []}
                    for iso in isoformAssign:
                        if sorted(isoformAssign[iso].items()) == sorted(comb):
                            identifiedComb.append(comb)
                            combination2iso[comb]["iso"].append(iso)
                            freqList = []
                            for i in comb:
                                freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
                            combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
                            combination2iso[comb].update({"freqList": freqList})
                        else:
                            freqList = []
                            for i in comb:
                                freqList.append(resultDict[i[0]][i[1]]["choiceFreq"])
                            combination2iso[comb].update({"freq": round(reduce(lambda x, y: x * y, freqList), 4)})
                            combination2iso[comb].update({"freqList": freqList})
                    readsInComb = set(itertools.chain.from_iterable([isoform2reads[x] for x in combination2iso[comb]["iso"]]))
                    counts = len(readsInComb)
                    combination2iso[comb].update({"counts": counts})
                    combination2iso[comb].update({"reads": readsInComb})
                identifiedComb = list(set(identifiedComb))
                combCounts = sorted([(combination2iso[x]["iso"], combination2iso[x]["counts"]) for x in identifiedComb],
                                    key=lambda x: x[1], reverse=True)
                # majorCombCount = 0
                # tmpSum = 0
                # isoList = []
                # for x in combCounts:
                #     majorCombCount += 1
                #     tmpSum += x[1]
                #     isoList.append(x[0])
                #     if tmpSum / float(sum([y[1] for y in combCounts])) >= 0.75:
                #         break
                if len(identifiedComb) == 0: continue
                choiceFreqList = []
                events = list(set([x[0] for x in list(itertools.chain.from_iterable(identifiedComb))]))
                for e in events:
                    choiceFreqList.append(",".join(map(str, [resultDict[e][x]["choiceFreq"] for x in resultDict[e]])))
                print >>asEnumerateOut, "\t".join(map(str, [chrom, strand, geneName, len(events),
                                          ";".join(events), \
                                          ";".join(choiceFreqList), len(combinations), \
                                          ",".join(map(str, [combination2iso[x]["counts"] for x in identifiedComb])), \
                                          ",".join(map(str, [combination2iso[x]["freq"] for x in identifiedComb])), \
                                          len(list(set(identifiedComb))), len(set(allReads)),
                                          ";".join([",".join(x[0]) for x in combCounts])]))
                for item in combCounts:
                    simFreq = 0
                    for x in combination2iso:
                        if combination2iso[x]["iso"] == item[0]:
                            simFreq = combination2iso[x]["freq"]
                    annoFlag = 0
                    tmpJuncDict = {}
                    for i in item[0]:
                        if isoformBed[i].juncChain not in tmpJuncDict:
                            tmpJuncDict[isoformBed[i].juncChain] = [i]
                        else:
                            tmpJuncDict[isoformBed[i].juncChain].append(i)
                    for j in tmpJuncDict:
                        for x in tmpJuncDict[j]:
                            if x in annoIsoformDict:
                                annoFlag = 1
                                break
                        if annoFlag == 1:
                            annotation = "Annotated"
                        else:
                            annotation = "Novel"
                        allReadCount = len(set(allReads))
                        tmpReads = set(itertools.chain.from_iterable([isoform2reads[x] for x in tmpJuncDict[j]]))
                        print >> isoformScoreOut, "\t".join(map(str, [geneName, ",".join(tmpJuncDict[j]), len(tmpReads),
                                                                      allReadCount, float(len(tmpReads))/allReadCount,
                                                                      simFreq, annotation]))
    asEnumerateOut.close()
    isoformScoreOut.close()
    return "asEnumerate.txt", "isoformScore.txt"

def getIsoTPM(quant_sf):
    iso2tpmDict = {}
    with open(quant_sf) as f:
        for line in f.readlines()[1:]:
            iso = line.strip("\n").split("\t")[0]
            tpm = line.strip("\n").split("\t")[3]
            iso2tpmDict.update({iso: float(tpm)})
    return iso2tpmDict

def quantIsoformWithSalmon(isoformScoreFile, isoformFile, collapsedTrans2reads, dataObj, refParams, dirSpec):
    isoform2reads = getDictFromFile(collapsedTrans2reads, sep="\t", inlineSep=",", valueCol=2)
    representIsoOut = open("representIso.bed", "w")
    with open(isoformScoreFile) as f:
        isoformBed = BedFile(isoformFile, type="bed12+").reads
        for line in f.readlines():
            isoforms = line.strip("\n").split("\t")[1]
            repIsoName = ""
            maxCount = 0
            for x in isoforms.split(","):
                if len(isoform2reads[x]) > maxCount:
                    maxCount = len(isoform2reads[x])
                    repIsoName = x
            repIso = copy.copy(isoformBed[x])
            repIso.name = isoforms
            print >>representIsoOut, str(repIso) + "\t" + repIsoName
    representIsoOut.close()

    cmd = '''cut -f 1-12 representIso.bed | bedtools getfasta -fi {} -bed - -name -split | seqkit replace -w 0 -p "(.*?):(.*)" -r '$1' > representIso.fa'''.format(refParams.ref_genome)
    subprocess.call(cmd, shell=True, executable="/bin/bash")

    logDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "log")
    cmd = "salmon index -t representIso.fa -i representIso_index --keepDuplicates -p {} 1>{}/salmon.index.log 2>&1".format(dataObj.single_run_threads, logDir)
    subprocess.call(cmd, shell=True)

    from preprocess import renameNGSdata2fastp, processRnaseq
    processRnaseq(dataObj=dataObj, threads=dataObj.single_run_threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
    renameNGSdata2fastp(dataObj=dataObj)
    quant_sf_list = []
    iso2tpm = {}
    if dataObj.ngs_reads_paired == "paired":
        leftReadsRepeats = dataObj.ngs_left_reads.split(";")
        rightReadsRepeats = dataObj.ngs_right_reads.split(";")
        for i in range(len(leftReadsRepeats)):
            leftReads = " ".join([r.strip() for r in leftReadsRepeats[i].split(",")])
            rightReads = " ".join([r.strip() for r in rightReadsRepeats[i].split(",")])
            salmonOut = "repeat{}.salmon_quant".format(i)
            cmd = "salmon quant -l A -i representIso_index -1 {} -2 {} -p {} -o {} --consensusSlack 0.5 --preMergeChainSubThresh 0.9 1>{}/{}_log 2>&1"
            cmd = cmd.format(leftReads, rightReads, dataObj.single_run_threads, salmonOut, logDir, salmonOut)
            subprocess.call(cmd, shell=True)
            quant_sf_list.append(os.path.abspath("{}/quant.sf".format(salmonOut)))
            iso2tpm.update({salmonOut: getIsoTPM("{}/quant.sf".format(salmonOut))})
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            singleReadsRepeats = dataObj.ngs_left_reads.split(";")
        elif dataObj.ngs_left_reads == None and dataObj.ngs_right_reads:
            singleReadsRepeats = dataObj.ngs_right_reads.split(";")
        else:
            raise Exception("The NGS reads type you input is 'single', but the actually it seems like 'paired' one, please check it!")

        for i in range(len(singleReadsRepeats)):
            singleReads = ",".join([i.strip() for i in singleReadsRepeats[i].split(",")])
            salmonOut = "repeat{}.salmon_quant".format(i)
            cmd = "salmon quant -l A -i representIso_index -r {} -p {} -o {} --consensusSlack 0.5 --preMergeChainSubThresh 0.9 1>{}/{}_log 2>&1"
            cmd = cmd.format(singleReads, dataObj.single_run_threads, salmonOut, logDir, salmonOut)
            subprocess.call(cmd, shell=True)
            quant_sf_list.append(os.path.abspath("{}/quant.sf".format(salmonOut)))
            iso2tpm.update({salmonOut: getIsoTPM("{}/quant.sf".format(salmonOut))})

    newScoreOut = open("isoformScore.tpm.txt", "w")
    print >> newScoreOut, "\t".join(["Gene", "Isoform", "PB_count", "All_PB_count", "Freq_real", "Freq_comb", "TPM_min", "TPM_max", "TPM_mean", "Annotation"])
    with open(isoformScoreFile) as f:
        for line in f.readlines():
            infoList = line.strip("\n").split("\t")
            iso = infoList[1]
            tpmList = [iso2tpm[x][iso] for x in iso2tpm]
            minTPM = min(tpmList)
            maxTPM = max(tpmList)
            meanTPM = sum(tpmList)/float(len(tpmList))
            print >> newScoreOut, "\t".join(infoList[0:-1]) + "\t" + "\t".join(map(str, [minTPM, maxTPM, meanTPM])) + "\t" + infoList[-1]
    newScoreOut.close()

def getHqIsoCombs():
    isoformScoreFile = os.path.join(os.getcwd(), "isoformScore.tpm.txt")
    isoformFile = os.path.join("../refine/", "tofu.collapsed.assigned.bed12+")
    isoformScoreDict = {}
    with open(isoformScoreFile) as f:
        for line in f.readlines()[1:]:
            infoList = line.strip("\n").split("\t")
            if int(infoList[3]) < 10: continue
            if infoList[0] not in isoformScoreDict:
                # isoformScoreDict[infoList[0]] = {}
                isoformScoreDict[infoList[0]] = {infoList[-1]: [infoList]}
            elif infoList[-1] not in isoformScoreDict[infoList[0]]:
                isoformScoreDict[infoList[0]][infoList[-1]] = [infoList]
            else:
                isoformScoreDict[infoList[0]][infoList[-1]].append(infoList)

    isoformBed = BedFile(isoformFile, type="bed12+").reads
    out = open("hqIsoCombs.txt", "w")
    for g in isoformScoreDict:
        if len(isoformScoreDict[g]) < 2: continue
        annoIsos = isoformScoreDict[g]["Annotated"]
        novelIsos = isoformScoreDict[g]["Novel"]
        annoIsos = sorted(annoIsos, key=lambda x: int(x[2]), reverse=True)
        flag = 0
        for iso in novelIsos:
            if int(iso[2]) > int(annoIsos[0][2]) * 1.5:
                flag = 1
                if len(isoformBed[iso[1].split(",")[0]].introns) > 1:
                    print >> out, "\t".join(iso) + "\t" + "{}:{}".format(isoformBed[iso[1].split(",")[0]].chrom, isoformBed[iso[1].split(",")[0]].juncChain)
                else:
                    refGene = list(set(isoformBed[iso[1].split(",")[0]].otherList[3].split(",")))[0]
                    print >> out, "\t".join(iso) + "\t" + "{}:{}".format(isoformBed[iso[1].split(",")[0]].chrom, "Novel_"+refGene)
            elif int(annoIsos[0][2]) < int(iso[2]) and int(iso[2]) < 1.5 * int(annoIsos[0][2]):
                if float(iso[-2]) > float(annoIsos[0][-2]):
                    flag = 1
                    if len(isoformBed[iso[1].split(",")[0]].introns) > 1:
                        print >> out, "\t".join(iso) + "\t" + "{}:{}".format(isoformBed[iso[1].split(",")[0]].chrom, isoformBed[iso[1].split(",")[0]].juncChain)
                    else:
                        refGene = list(set(isoformBed[iso[1].split(",")[0]].otherList[3].split(",")))[0]
                        print >> out, "\t".join(iso) + "\t" + "{}:{}".format(isoformBed[iso[1].split(",")[0]].chrom, "Novel_" + refGene)
        if flag:
            if len(isoformBed[annoIsos[0][1].split(",")[0]].introns) > 1:
                print >> out, "\t".join(annoIsos[0]) + "\t" + "{}:{}".format(isoformBed[annoIsos[0][1].split(",")[0]].chrom, isoformBed[annoIsos[0][1].split(",")[0]].juncChain)
            else:
                refGene = list(set(isoformBed[annoIsos[0][1].split(",")[0]].otherList[3].split(",")))[0]
                print >> out, "\t".join(annoIsos[0]) + "\t" + "{}:{}".format(isoformBed[annoIsos[0][1].split(",")[0]].chrom, "Known_" + refGene)
            # print >> out, "\t".join(annoIsos[0]) + "\t" + "{}:{}".format(isoformBed[annoIsos[0][1].split(",")[0]].chrom, isoformBed[annoIsos[0][1].split(",")[0]].juncChain)
    out.close()
# irFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050245.IR.PB.bed6+"
# seFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050245.SE.PB.bed12+"
# a3ssFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050245.A3SS.PB.bed6+"
# a5ssFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050245.A5SS.PB.bed6+"
# paFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/test.py.pa"
# isoformFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050245.deSingleExonIsoform.bed12+"
# isoform2readsFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/tofu.collapsed.group.txt"
#
irFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/Zm00001d050929.IR.bed6+"
seFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/test.py.pa"
a3ssFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/test.py.pa"
a5ssFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/test.py.pa"
paFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/test.py.pa"
isoformFile = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/PB.7806.bed12+"
collapsedTrans2reads = "/home/xufeng/xufeng/iso-seq/iFLAS_toolkit/test_data/PB.7806.collapsed.group.txt"
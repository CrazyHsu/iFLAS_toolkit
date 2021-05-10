from commonFuncs import *
from commonObjs import *
import itertools

def asEnumerate(irFile, seFile, a3ssFile, a5ssFile, paFile, isoformFile, isoform2readsFile=None):
    isoform2reads = {}
    with open(isoform2readsFile) as f:
        for line in f.readlines():
            isoform, reads = line.strip("\n").split("\t")
            isoform2reads[isoform] = reads.split(",")
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
            geneName = name.split(":")[0]
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
            geneName = name.split(":")[0]
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
            geneName = name.split(":")[0]
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
            geneName = name.split(":")[0]
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
        # print newJuncs
        # print [x.name for x in isosInGene]
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
                            choiceName = "choice_ir"
                    if choiceName not in eventCombination[eventName]:
                        eventCombination[eventName][choiceName] = [iso]
                    else:
                        eventCombination[eventName][choiceName].append(iso)
                if iso.name not in isoformAssign:
                    isoformAssign[iso.name] = {eventName: choiceName}
                else:
                    isoformAssign[iso.name].update({eventName: choiceName})
        # print isoformAssign
        # for eventName in eventCombination:
        #     for choiceName in eventCombination[eventName]:
        #         print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]]
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
                readsInChoice = list(set(itertools.chain.from_iterable([isoform2reads[i.name] for i in allChoices[choiceName]])))
                choiceFreq = len(readsInChoice)/float(len(allReadsInEvent))
                resultDict[eventName].update({choiceName: {"readsInChoice": readsInChoice, "choiceFreq": choiceFreq}})
                # print eventName, choiceName, [x.name for x in eventCombination[eventName][choiceName]], len(readsInChoice), len(allReadsInEvent)
            allEventChoice.append(tmpList)
            choiceFreqList.append(",".join(map(str, [resultDict[eventName][x]["choiceFreq"] for x in resultDict[eventName]])))
        combinations = list(itertools.product(*allEventChoice))
        # print combinations, len(combinations)
        # print choiceFreqList
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
        # print combination2iso.values()
        identifiedComb = list(set(identifiedComb))
        combCounts = sorted([(combination2iso[x]["iso"], combination2iso[x]["counts"]) for x in identifiedComb], key=lambda x: x[1], reverse=True)
        majorCombCount = 0
        tmpSum = 0
        isoList = []
        for x in combCounts:
            majorCombCount += 1
            tmpSum += x
            isoList.append(x[0])
            if tmpSum/float(sum([y[1] for y in combCounts])) >= 0.75:
                break
        print "\t".join(map(str, [chrom, strand, geneName, len(eventCombination), ";".join(eventCombination.keys()), \
            ";".join(choiceFreqList), len(combinations), \
            ",".join(map(str, [combination2iso[x]["counts"] for x in identifiedComb])), \
            ",".join(map(str, [combination2iso[x]["freq"] for x in identifiedComb])), \
            len(list(set(identifiedComb))), len(list(set(allReads))), majorCombCount, ";".join([",".join(x) for x in isoList])]))

from collections import defaultdict
import time
import matplotlib.pyplot as plt


def createCharDictionary():
    charDictionary = {}
    with open('blosum_62.txt', 'r', encoding='utf8') as f:
        tableLine = 0
        for line in f:
            if line[0] != '#':
                words = [word.format() for word in line.split()]
                num = len(words)
                i = 0
                for word in line.split():
                    charDictionary[word] = i
                    i += 1
                return charDictionary


def createCostTable():
    with open('blosum_62.txt', 'r', encoding='utf8') as f:
        firstTableLine, tableLine = True, 0
        for line in f:
            if line[0] != '#':
                if firstTableLine:
                    size = len(charDictionary)
                    costTable = [[0 for x in range(size)] for y in range(size)]
                    firstTableLine = False
                else:
                    words = [word.format() for word in line.split()]
                    key = words[0]
                    val = charDictionary[key]
                    for i in range(1, len(words)):
                        costTable[tableLine][i - 1] = words[i]
                    tableLine += 1
    return costTable


def createString(fileName):
    string = ''
    with open(fileName, 'r', encoding='utf8') as f:
        for line in f:
            if line[0] != '>':
                string += line[:-1]
    return string


def databasePreProcess(database, w):
    wmers = defaultdict(list)
    for i in range(0, len(database) - w + 1):
        key = database[i:i + w]
        wmers[key].append(i)
    wmers = dict((k, tuple(v)) for k, v in wmers.items())
    return wmers


def calculateScore(str1, str2):
    score = 0;
    for i in range(0, len(str1)):
        key1, key2 = str1[i], str2[i]
        val1, val2 = charDictionary[key1], charDictionary[key2]
        score += int(costTable[val1][val2])
    return score


def findNextNeighbors(neighbors, key, potentialN, T, w, j):  ### Bonus ###
    for k, v in charDictionary.items():
        potentialNeighbor = potentialN[:j] + k + potentialN[j + 1:]
        score = calculateScore(key, potentialNeighbor)
        if score >= T:
            if not potentialNeighbor in neighbors[key]:
                neighbors[key].append(potentialNeighbor)
            if j + 1 < w:
                neighbors = findNextNeighbors(neighbors, key, potentialNeighbor, T, w, j + 1)
    return neighbors


def findingNeighbors(query, w, T):
    neighbors = defaultdict(list)
    for i in range(0, len(query) - w + 1):
        key = query[i:i + w]
        for j in range(0, w):
            neighbors = findNextNeighbors(neighbors, key, key, T, w, j)
    neighbors = dict((k, tuple(v)) for k, v in neighbors.items())
    return neighbors


def findingHits(wmers, neighbors):
    hits = defaultdict(list)
    for k_wmerQ, v_wmerQ in neighbors.items():
        for neighborString in v_wmerQ:
            if neighborString in wmers:
                value = wmers[neighborString]
                for v in value:
                    hits[k_wmerQ].append(v)
    hits = dict((k, tuple(v)) for k, v in hits.items())
    return hits


def checkExistance(MSPs, msp):
    # msp = [startIndexDBmax, endIndexDBmax, startIndexQmax, endIndexQmax, maxScore]
    if msp in MSPs:
        return MSPs
    toAdd = True
    for m in MSPs:
        # msp is in m and msp's score is higher - delete m and add msp
        if msp[0] >= m[0] and msp[1] <= m[1] and msp[2] >= m[2] and msp[3] <= m[3] and msp[4] > m[4]:
            MSPs.remove(m)
            MSPs.append(msp)
            toAdd = False
        # msp is in m and m's score is higher - no actions
        elif msp[0] >= m[0] and msp[1] <= m[1] and msp[2] >= m[2] and msp[3] <= m[3] and msp[4] < m[4]:
            toAdd = False
        # m is in msp and msp's score is higher - delete m and add msp
        elif msp[0] <= m[0] and msp[1] >= m[1] and msp[2] <= m[2] and msp[3] >= m[3] and msp[4] > m[4]:
            MSPs.remove(m)
            MSPs.append(msp)
            toAdd = False
        # m is in msp and m's score is higher - no actions
        elif msp[0] <= m[0] and msp[1] >= m[1] and msp[2] <= m[2] and msp[3] >= m[3] and msp[4] < m[4]:
            toAdd = False
    if toAdd:
        MSPs.append(msp)
    return MSPs


def ExtentionOfHSPs(query, database, hits, X, wmersQ, w):
    # MSP is a list [startDBindex,endDBindex,startQindex,endQindex,score]
    MSPs = []
    for k_hits, v_hits in hits.items():
        for k_wmersQ, v_wmersQ in wmersQ.items():
            if k_hits == k_wmersQ:
                for v_wmer in v_wmersQ:
                    for v_hit in v_hits:
                        startIndexQmax = v_wmer
                        endIndexQmax = v_wmer + w
                        startIndexDBmax = v_hit
                        endIndexDBmax = v_hit + w
                        maxScore = calculateScore(query[startIndexQmax:endIndexQmax],
                                                  database[startIndexDBmax:endIndexDBmax])
                        startIndexQ = startIndexQmax
                        endIndexQ = endIndexQmax
                        startIndexDB = startIndexDBmax
                        endIndexDB = endIndexDBmax
                        score = maxScore
                        while startIndexQ >= 1 and startIndexDB >= 1:
                            startIndexQ -= 1
                            startIndexDB -= 1
                            pairScore = calculateScore(query[startIndexQ], database[startIndexDB])
                            newScore = score + pairScore
                            if maxScore - X > newScore:
                                break
                            if newScore >= maxScore:
                                maxScore = newScore
                                startIndexQmax = startIndexQ
                                startIndexDBmax = startIndexDB
                            score = newScore
                        score = maxScore
                        while endIndexQ <= (len(query) - 1) and endIndexDB <= (len(database) - 1):
                            pairScore = calculateScore(query[endIndexQ], database[endIndexDB])
                            newScore = score + pairScore
                            if maxScore - X > newScore:
                                break
                            if newScore >= maxScore:
                                maxScore = newScore
                                endIndexQmax = endIndexQ + 1
                                endIndexDBmax = endIndexDB + 1
                            score = newScore
                            endIndexQ += 1
                            endIndexDB += 1
                        msp = [startIndexDBmax, endIndexDBmax, startIndexQmax, endIndexQmax, maxScore]
                        MSPs = checkExistance(MSPs, msp)
    return MSPs


def topNMSP(MSPs, N):
    MSPs.sort(key=lambda x: x[-1], reverse=True)
    if N == -1:
        return MSPs
    return MSPs[:N]


# finding 3 top scoring MSPs
charDictionary = createCharDictionary()
costTable = createCostTable()
query = createString('query1.fasta')
database = createString('db1.fasta')

# parameters - chosen by the user
w = 3
T = 10
X = 5
N = 6

startTime = time.time()

wmers = databasePreProcess(database,w)
wmersQ = databasePreProcess(query,w)
neighbors = findingNeighbors(query,w,T)
hits = findingHits(wmers,neighbors)
MSPs = ExtentionOfHSPs(query,database,hits,X,wmersQ,w)
topNMSPs = topNMSP(MSPs,N)

endTime = time.time()
runtime = endTime - startTime

print('top 5 MSPs: ',topNMSPs)
print('format: [DB start, DB end, query start, query end, score]')
print('runtime (seconds): ',runtime)
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


def createCalculationTable(query, database):
    calculationTable = [[0 for x in range(len(database) + 1)] for y in range(len(query) + 1)]
    maxScore, final_i, final_j = float('-inf'), 0, 0
    for i in range(1, len(query) + 1):
        for j in range(1, len(database) + 1):
            keyQ, keyD = query[i - 1], database[j - 1]
            valQ, valD = charDictionary[keyQ], charDictionary[keyD]
            goDownCost = costTable[valQ][len(costTable) - 1]
            goRightCost = costTable[valD][len(costTable) - 1]
            goVerticalCost = costTable[valQ][valD]
            calculationTable[i][j] = max(0, calculationTable[i][j - 1] + int(goRightCost),
                                         calculationTable[i - 1][j] + int(goDownCost),
                                         calculationTable[i - 1][j - 1] + int(goVerticalCost))
            if calculationTable[i][j] > maxScore:
                maxScore = calculationTable[i][j]
                final_i, final_j = i, j
    return maxScore, final_i, final_j, calculationTable


def findStrings(query, database, maxScore, final_i, final_j, calculationTable):
    strQ, strD, i, j = '', '', final_i, final_j
    while calculationTable[i][j] > 0:
        keyQ, keyD = query[i - 1], database[j - 1]
        valQ, valD = charDictionary[keyQ], charDictionary[keyD]
        goDownCost = costTable[valQ][len(costTable) - 1]
        goRightCost = costTable[valD][len(costTable) - 1]
        goVerticalCost = costTable[valQ][valD]
        if calculationTable[i][j] == calculationTable[i][j - 1] + int(goRightCost):
            strD = database[j - 1] + strD
            strQ = '_' + strQ
            j -= 1
        elif calculationTable[i][j] == calculationTable[i - 1][j] + int(goDownCost):
            strD = '_' + strD
            strQ = query[i - 1] + strQ
            i -= 1
        else:
            strD = database[j - 1] + strD
            strQ = query[i - 1] + strQ
            j -= 1
            i -= 1
    return i, j, strQ, strD


startTime = time.time()
charDictionary = createCharDictionary()
costTable = createCostTable()
query = createString('query2.fasta')
database = createString('db1.fasta')
maxScore, final_i, final_j, calculationTable = createCalculationTable(query, database)
first_i, first_j, strQ, strD = findStrings(query, database, maxScore, final_i, final_j, calculationTable)
endTime = time.time()
runtime = endTime - startTime

print(strQ)
print(strD)
print('score = ',maxScore)
print('s1 (',first_i,',',final_i,')')
print('s2 (',first_j,',',final_j,')')
print('runtime is: ',runtime/60,' minutes')
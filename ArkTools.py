#!/usr/bin/env python3
from collections import defaultdict
import re
import statistics
import numpy as np
import os
import random


def AAI(seq1, seq2):
    counter = 0
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                counter += 1
    return (counter/len(seq1))*100


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count/len(seq)


def ID(seq1, seq2):
    counter = 0
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                counter += 1
    return (counter/len(seq1))*100


def ISsearch(ls):
    IS = "unknown-class"
    for i in ls:
        if re.findall(r'IS', i):
            IS = i
            break

    if re.findall(r'like', IS):
        IS = IS.split("-")[0]
    return IS


def LCMSMS(trypsinFragsList):
    count = 0
    for i in trypsinFragsList:
        if len(i) >= 7:
            count += 1
    return count


def Num(ls):
    outputNum = "0"
    for i in ls:
        try:
            int(i)
            outputNum = i
        except ValueError:
            pass
    return outputNum


def Parrot(seq):
    ls = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        ls.append(codon)
    return ls


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def RemoveLeadingSpaces(line):
    counter = 0
    newLine = ''
    for i in line:
        if i == " ":
            if counter == 0:
                pass
            else:
                newLine += i
        else:
            counter += 1
            newLine += i
    return newLine


def SUM(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def SeqCoord(seq, start, end):
    return seq[start-1:end]


def Strip(ls):
    outList = []
    for i in ls:
        gene = i.split("|")[0]
        outList.append(gene)
    return outList


def Sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def Unique(ls):
    unqList = []
    for i in ls:
        if i not in unqList:
            unqList.append(i)
    return unqList


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def alnCheckC(seq1, seq2, slack):
    count = 0
    countGaps = 0
    for i in range(len(seq1)-1, len(seq1)-slack, -1):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                count += 1
        else:
            countGaps += 1
    return count/slack, countGaps/slack


def alnCheckN(seq1, seq2, slack):
    count = 0
    countGaps = 0
    for i in range(0, slack):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                count += 1
        else:
            countGaps += 1
    return count/slack, countGaps/slack


def anot(ls):
    outList = []
    for i in ls[1:]:
        if "[" not in list(i):
            outList.append(i + " ")
        else:
            break
    outStr = "".join(outList)
    return outStr[0:len(outStr)-1]


def asNumeric(list):
    newList = []
    for i in list:
        newList.append(float(i))
    return newList


def ave(ls):
    count = 0
    for i in ls:
        try:
            count += float(i)
        except ValueError:
            pass
    return count/len(ls)


def backTrim(seq1, seq2):
    '''
        seq1 needs to be the reference
        '''
    newSeq1 = seq1
    newSeq2 = seq2
    for i in range(len(seq1), 0, -3):
        if lastItem(newSeq2) == "-":
            newSeq1 = seq1[0:i]
            newSeq2 = seq2[0:i]
        else:
            break

    return newSeq1, newSeq2


def breakdown(ls):
    Dict = defaultdict(list)
    for i in ls:
        Dict[i].append("i")
    Dict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for j in Dict.keys():
        Dict2[j] = len(Dict[j])/len(ls)
    return Dict2


def capitalizeCodon(codon):
    codonOut = ''
    for i in codon:
        if i == "a":
            codonOut += "A"
        elif i == "g":
            codonOut += "G"
        elif i == "c":
            codonOut += "C"
        elif i == "t":
            codonOut += "T"
        elif i == "u":
            codonOut += "U"
        else:
            codonOut += i
    return codonOut


def checkDFE1(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"]:
                count += 1
    return count


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def codonTable(seq):
    Dict = defaultdict(lambda: defaultdict(list))
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            Dict[CodonTable[codon]][codon].append(codon)
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return Dict


def collectGaps(seq):
    ls = []
    word = ''
    for i in seq:
        if i != "-":
            if len(word) > 0:
                ls.append(word)
            word = ''
        else:
            word += "-"
    return ls


def combineLists(ls1, ls2):
    combinedList = []
    for i in ls1:
        if i not in combinedList:
            combinedList.append(i)
    for i in ls2:
        if i not in combinedList:
            combinedList.append(i)
    return combinedList


def compareLists(ls1, ls2):
    sharedList = []
    for i in ls2:
        if i in ls1:
            sharedList.append(i)
    return sharedList


def core(listOfLists):
    Dict = defaultdict(list)
    for i in listOfLists:
        ls = i
        for j in ls:
            Dict[j].append("1")
    coreList = []
    totalList = []
    for i in Dict.keys():
        totalList.append(i)
        if len(Dict[i]) == len(listOfLists):
            coreList.append(i)
    return coreList, totalList


def customfilter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def cutter(string):
    outString = ""
    for i in string:
        if i in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            break
        else:
            outString += i
    return outString


def deAligner(inList):
    outList = []
    for i in inList:
        if i != "-":
            outList.append(i)
    outStr = "".join(outList)
    return outStr


def deExt(string):
    ls = string.split(".")
    ls2 = ls[0:len(ls) - 1]
    outstring = "".join(ls2)
    return outstring


def deString(string):
    newString = ''
    for i in string:
        try:
            newString += str(int(i))
        except ValueError:
            break
    return newString


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def diff(num1, num2):
    lower = sorted([num1, num2])[0]
    higher = sorted([num1, num2])[1]
    return higher - lower


def digitize(string):
    outStr = ''
    for i in string:
        try:
            int(i)
            outStr += str(i)
        except ValueError:
            pass
    return (int(outStr))


def extender(listOfcoords):
    counter = 0
    max = 0
    outlist = []
    count = 0
    for i in listOfcoords:
        if counter == 0:
            counter += 1
            start = int(i[0])
            end = int(i[1])
        else:
            if int(i[0]) < end:
                if int(i[1]) > end:
                    end = int(i[1])
                    count = 0
            else:
                count += 1
                outlist.append([start, end])
                start = int(i[0])
                end = int(i[1])
    outlist.append([start, end])
    return outlist


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fastaRename(fasta_file):
    counter = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                counter += 1
                header = header + "_" + str(counter)
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                counter += 1
                header = header + "_" + str(counter)
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def findDictKey(regex, dict):
    count = 0
    for i in dict.keys():
        if re.findall(regex, dict[i]):
            count += 1
    if count > 0:
        return True


def firstNonspace(ls):
    for i in ls:
        if i != "":
            break
    return i


def firstNum(string):
    outputNum = []
    for i in string:
        try:
            int(i)
            outputNum.append(i)
        except ValueError:
            break
    Num = "".join(outputNum)
    return Num


def gc(seq):
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            gc += 1
    return gc / len(seq)


def howMany(list, item):
    counter = 0
    for i in list:
        if i == item:
            counter += 1
    return counter


def howManyNot(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def listFetcher(ls, term):
    for i in ls:
        # print(i)
        for j in i.split("; "):
            # print(j)
            for k in j.split(", "):
                # print(k)
                if re.findall(term, k):
                    return k


def listFinder(ls, term):
    counter = 0
    for i in ls:
        for j in i.split("; "):
            if re.findall(term, j):
                counter += 1
        for j in i.split(", "):
            if re.findall(term, j):
                counter += 1
    if counter > 0:
        return True


def listSearch(item, ls):
    count = 0
    for i in ls:
        if re.findall(i, item):
            count += 1
    if count > 0:
        return True


def localize(item, ls):
    count = 0
    for i in ls:
        if i == item:
            break
        else:
            count += 1
    return count


def mass(seq):
    dal = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    dal["A"] = 89
    dal["R"] = 174
    dal["N"] = 132
    dal["D"] = 133
    dal["C"] = 121
    dal["Q"] = 146
    dal["E"] = 147
    dal["G"] = 75
    dal["H"] = 155
    dal["I"] = 131
    dal["L"] = 131
    dal["K"] = 146
    dal["M"] = 149
    dal["F"] = 165
    dal["P"] = 115
    dal["S"] = 105
    dal["T"] = 119
    dal["W"] = 204
    dal["Y"] = 181
    dal["V"] = 117
    dal["*"] = 0
    mw = 0
    for i in seq:
        mw += int(dal[i])
    return mw/1000


def mean(ls):
    count = 0
    for i in ls:
        count += int(i)
    ave = count/len(ls)
    return ave


def mostCommon(Dict):
    for i in sorted(Dict.keys()):
        highest = 0
        for j in Dict:
            aa = j
            count = Dict[j]
            if count > highest:
                highest = count
                consensus = aa
            if count == highest:
                if aa != "-":
                    consensus = aa
        return consensus


def numerateList(ls):
    total = 0
    for i in ls:
        total += len(i)
    return total


def parse(ls, id):
    str = ""
    for i in ls:
        if re.findall(id, i):
            str = i
    return str


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def readCompile(listOfcoords, range):
    outList = []
    start = range[0]
    end = range[1]
    count = 0
    for i in listOfcoords:
        if i[0] >= start and i[1] <= end:
            outList.append(count)
            count += 1
        else:
            count += 1
    return outList


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def reject_outliers(data, m):
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - m * s < e < u + m * s)]
    return filtered


def remove(stringOrlist, item):
    emptyList = []
    for i in stringOrlist:
        if i != item:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeString(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeList(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    # outString = "".join(emptyList)
    return emptyList


def removeInt(ls):
    newList = []
    for i in ls:
        try:
            number = float(i)
        except ValueError:
            newList.append(i)
    return newList


def removeRegex(ls, regex):
    newList = []
    for i in ls:
        if re.findall(regex, i):
            pass
        else:
            newList.append(i)
    return newList


def replaceString(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def replaceList(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    return emptyList


def reverse(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "C"
        elif nucleotide == "G":
            nucleotide = "G"
        elif nucleotide == "T":
            nucleotide = "T"
        elif nucleotide == "A":
            nucleotide = "A"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"

        elif nucleotide == "g":
            nucleotide = "c"
        elif nucleotide == "T":
            nucleotide = "a"
        elif nucleotide == "a":
            nucleotide = "t"
        elif nucleotide == "c":
            nucleotide = "g"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Reverse(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    Dict = defaultdict(lambda: defaultdict(list))
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def ribosome2(seq):
    NTs = ['t', 'c', 'a', 'g']
    stopCodons = ['taa', 'tag', 'tga']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def search(string, listofitems):
    outList = []
    for i in string:
        if i in listofitems:
            outList.append(i)
    return outList


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def startFinder(seq):
    start = ''
    for i in range(len(seq)):
        codon = (seq[i:i+3])
        if not re.findall(r'-', codon):
            start = i
            break
    return start


def stopfinder(seq):
    newseq = ''
    gaps = ''
    for i in seq:
        if i != '-':
            newseq += gaps
            newseq += i
            gaps = ''
        else:
            gaps += i
    return newseq, len(gaps)


def sum(list):
    sum = 0
    for i in list:
        i = int(i)
        sum += i
    return sum


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def tet(seq):
    Dict = defaultdict(list)
    NTs = ['T', 'C', 'A', 'G']
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    if NTs[i] in ["A", "G", "C", "T"] and NTs[j] in ["A", "G", "C", "T"] and NTs[k] in ["A", "G", "C", "T"] and NTs[l] in ["A", "G", "C", "T"]:
                        tet = NTs[i] + NTs[j] + NTs[k] + NTs[l]
                        Dict[tet] = []

    totalKmers = 0
    for m in range(len(seq)):
        TET = (seq[m:m+4])
        if len(TET) == 4:
            if TET[0] in ["A", "G", "C", "T"] and TET[1] in ["A", "G", "C", "T"] and TET[2] in ["A", "G", "C", "T"] and TET[3] in ["A", "G", "C", "T"]:
                Dict[TET].append(TET)
                totalKmers += 1

    return Dict, totalKmers


def tetfreq(seq):
    tet_Dict = defaultdict(list)
    Nucs = ["T","A","G","C"]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if Nucs[a] in ["A", "G", "C", "T"] and Nucs[b] in ["A", "G", "C", "T"] and Nucs[c] in ["A", "G", "C", "T"] and Nucs[d] in ["A", "G", "C", "T"]:
                        tet = Nucs[a] + Nucs[b] + Nucs[c] + Nucs[d]
                        tet_Dict[tet] = []

    total = 0
    for e in range(len(seq)):
        f = (seq[e:e+4])
        if len(f) == 4:
            if f[0] in ["A", "G", "C", "T"] and f[1] in ["A", "G", "C", "T"] and f[2] in ["A", "G", "C", "T"] and f[3] in ["A", "G", "C", "T"]:
                tet_Dict[f].append(f)
                total += 1
    totalkmers = total
    return tet_Dict, total


def thirdToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-2]:
        x = i
    return x


def topRange(listOfRanges):
    max = 0
    for i in listOfRanges:
        length = i[1] - i[0]
        if length > max:
            max = length
            bestPair = [i[0], i[1]]
    return [max, bestPair]



def trim(seq):
    newSeq = ''
    counter = 0
    for i in seq:
        if i == "-":
            if counter == 0:
                pass
            else:
                newSeq += i

        else:
            counter += 1
            newSeq += i
    return newSeq


def trimEnd(seq):
    for i in range(1, len(seq)):
        if seq[len(seq)-i] == "-":
            pass
        else:
            break
    index = i-1
    return seq[0:len(seq)-index]


def trimmer(seq1, seq2):
    '''
    seq1 needs to be the reference
    '''
    counter = 0
    newSeq1 = ""
    newSeq2 = ""
    for i in range(0, len(seq1), 3):
        if re.findall(r'-', seq1[i:i+3]) or re.findall(r'-', seq2[i:i+3]):
            if counter == 0:
                pass
            else:
                newSeq1 += seq1[i:i + 3]
                newSeq2 += seq2[i:i + 3]
        else:
            counter += 1
            newSeq1 += seq1[i:i+3]
            newSeq2 += seq2[i:i + 3]
    return newSeq1, newSeq2


def trypsin(seq):
    peptide = ''
    peptideLS = []
    for j in seq:
        if j in ["R", "L"]:
            peptide += j
            peptideLS.append(peptide)
            peptide = ""
        else:
            peptide += j
    peptideLS.append(peptide)
    return peptideLS


def unq(ls):
    new = []
    for i in ls:
        if i not in new:
            new.append(i)
    return new



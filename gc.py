#!/usr/bin/env python3
import statistics
from collections import defaultdict
import argparse
import textwrap
import re
import sys


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def deString(string):
    newString = ''
    for i in string:
        try:
            newString += str(int(i))
        except ValueError:
            break
    return newString


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


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
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


parser = argparse.ArgumentParser(
    prog="snp_table.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber;
    Middle Author Bioinformatics
    Please send comments and inquiries to ark@midauthorbio.com

    *******************************************************
    '''))

parser.add_argument('-fna', type=str, help="")

parser.add_argument('-faa', type=str, help="")

parser.add_argument('-out', type=str, help="")

args = parser.parse_args()


# out = open(args.out, "w")
try:
    contigs = open(args.fna)
    contigs = fasta(contigs)
    faa = open(args.faa)
    faa = fasta(faa)
except FileNotFoundError:
    sys.exit()

GClist = []
totalLength = 0
for i in contigs.keys():
    seq = contigs[i]
    GC = (seq.count("G") + seq.count("C")) / len(seq)
    GClist.append(GC)
    totalLength += len(seq)
header = i

try:
    GCmean = round(statistics.mean(GClist), 4)
    GCvar = round(statistics.stdev(GClist), 4)
except statistics.StatisticsError:
    GCmean = round(GClist[0], 4)
    GCvar = 0

Dict = defaultdict(list)
SEQ = ""
for i in faa.keys():
    seq = faa[i]
    SEQ += seq

A = round(SEQ.count("A") / len(SEQ), 4)
G = round(SEQ.count("G") / len(SEQ), 4)
V = round(SEQ.count("V") / len(SEQ), 4)
I = round(SEQ.count("I") / len(SEQ), 4)
L = round(SEQ.count("L") / len(SEQ), 4)
M = round(SEQ.count("M") / len(SEQ), 4)
F = round(SEQ.count("F") / len(SEQ), 4)
Y = round(SEQ.count("Y") / len(SEQ), 4)
W = round(SEQ.count("W") / len(SEQ), 4)
S = round(SEQ.count("S") / len(SEQ), 4)
T = round(SEQ.count("T") / len(SEQ), 4)
N = round(SEQ.count("N") / len(SEQ), 4)
Q = round(SEQ.count("Q") / len(SEQ), 4)
C = round(SEQ.count("C") / len(SEQ), 4)
P = round(SEQ.count("P") / len(SEQ), 4)
R = round(SEQ.count("R") / len(SEQ), 4)
H = round(SEQ.count("H") / len(SEQ), 4)
K = round(SEQ.count("K") / len(SEQ), 4)
D = round(SEQ.count("D") / len(SEQ), 4)
E = round(SEQ.count("E") / len(SEQ), 4)

accession = args.fna
accession = accession.split(".")[0] + "." + accession.split(".")[1]
print(str(accession) + "\t" + str(GCmean) + "\t" + str(totalLength) + "\t" + str(A) + "\t" + str(G) + "\t" + str(V) + "\t" +
      str(I) + "\t" + str(L) + "\t" + str(M) + "\t" + str(F) + "\t" + str(Y) + "\t" + str(W) + "\t" + str(S) + "\t" +
      str(T) + "\t" + str(N) + "\t" + str(Q) + "\t" + str(C) + "\t" + str(P) + "\t" + str(R) + "\t" +
      str(H) + "\t" + str(K) + "\t" + str(D) + "\t" + str(E) + "\t" + header)


















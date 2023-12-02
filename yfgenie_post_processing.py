#!/usr/bin/env python3
from collections import defaultdict
import re
import statistics
import numpy as np
import os
from ArkTools import *


# function to get unique values
def Unique(list1):
    # initialize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list


def AminoacylatioCheck(list):
    counter = 0
    if "TIGR00440" in list and "CytCoxidase_aa3_coxB" in list:
        counter = 1
    if "Cyt_aa3_quinol_oxidase_QoxA" in list and "Cyt_aa3_quinol_oxidase_QoxB" in list:
        counter = 1
    if "Cyt_bd_oxidase_CydB.hmm" in list:
        counter = 1
    if "CytCoxidase_cbb3_ccoN" in list and "CytCoxidase_cbb3_ccoO" in list:
        counter = 1
    return counter


def OxCheck(list):
    counter = 0
    if "CytCoxidase_aa3_coxA" in list and "CytCoxidase_aa3_coxB" in list:
        counter = 1
    if "Cyt_aa3_quinol_oxidase_QoxA" in list and "Cyt_aa3_quinol_oxidase_QoxB" in list:
        counter = 1
    if "Cyt_bd_oxidase_CydB.hmm" in list:
        counter = 1
    if "CytCoxidase_cbb3_ccoN" in list and "CytCoxidase_cbb3_ccoO" in list:
        counter = 1
    return counter


def NdhCheck(list):
    counter = 0
    if len(list) > 1:
        counter = 1
    return counter

def Sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count

def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def ribosome(seq):
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
            prot.append("")
    protein = ("".join(prot))
    return protein, CodonTable


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

#########################################################################################
### CHAPTER 1
catDict = defaultdict(lambda: defaultdict(list))
keyDict = defaultdict(lambda: defaultdict(lambda: '-'))
key = open("HMM_key.csv")
for i in key:
    ls = i.rstrip().split(",")
    if ls[2] != "category":
        keyDict[ls[3].split(".")[0]]["gene"] = ls[0]
        keyDict[ls[3].split(".")[0]]["cat"] = ls[1]
        keyDict[ls[3].split(".")[0]]["section"] = ls[2]

        catDict[ls[2]][ls[1]].append(ls[3].split(".")[0])

ncbiDict = defaultdict(lambda: 'EMPTY')
ncbi = open("ncbi_assembly_info.1.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.2.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.3.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.4.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

catDict3 = defaultdict(list)
catDict2 = defaultdict(lambda: '-')
cats = open("CATEGORIES2.txt")
for i in cats:
    if re.findall(r':', i):
        category = (i.rstrip().split(":")[0])
    else:
        ls = (i.rstrip().split(" "))
        for j in ls:
            if re.findall(r'hmm', j) or re.findall(r'HMM', j):
                catDict2[j.split(".")[0]] = category
                catDict3[category].append(j.split(".")[0])

print("---")
allDict = defaultdict(lambda: defaultdict(lambda: '0'))
mainDir = os.listdir("HMM_results")
for i in mainDir:
    if i != ".DS_Store":
        catDir = os.listdir("HMM_results/%s" % i)
        for j in catDir:
            if re.findall(r'csv', j):
                accession = j.split(".")[0]
                csv = open("HMM_results/%s/%s" % (i, j))
                for k in csv:
                    if not re.match(r'#', k):
                        ls = k.rstrip().split(",")
                        if ls[4] != "HMM_accession":
                            if ls[4] == "-":
                                KEY = ls[5].split(".")[0]
                                if keyDict[KEY]["cat"] != "-":
                                    allDict[accession][KEY] = "1"
                            else:
                                KEY = ls[4]
                                if keyDict[KEY]["cat"] != "-":
                                    allDict[accession][KEY] = "1"

print("---")
ncDict = defaultdict(lambda: defaultdict(lambda: '0'))
mainDir = os.listdir("HMM_results")
for i in mainDir:
    if i != ".DS_Store":
        catDir = os.listdir("HMM_results/%s" % i)
        for j in catDir:
            if re.findall(r'csv', j):
                accession = j.split(".")[0]
                csv = open("HMM_results/%s/%s" % (i, j))
                for k in csv:
                    if not re.match(r'#', k):
                        ls = k.rstrip().split(",")
                        if ls[4] != "HMM_accession":
                            if ls[4] == "-":
                                KEY = ls[5].split(".")[0]
                                if keyDict[KEY]["cat"] != "-":
                                    ncDict[accession][KEY] = "1"
                            else:
                                KEY = ls[4]
                                if keyDict[KEY]["cat"] != "-":
                                    ncDict[accession][KEY] = "1"

print("---")
tcDict = defaultdict(lambda: defaultdict(lambda: '0'))
mainDir = os.listdir("HMM_results/nc")
for i in mainDir:
    if i != ".DS_Store":
        catDir = os.listdir("HMM_results/%s" % i)
        for j in catDir:
            if re.findall(r'csv', j):
                accession = j.split(".")[0]
                csv = open("HMM_results/%s/%s" % (i, j))
                for k in csv:
                    if not re.match(r'#', k):
                        ls = k.rstrip().split(",")
                        if ls[4] != "HMM_accession":
                            if ls[4] == "-":
                                KEY = ls[5].split(".")[0]
                                if keyDict[KEY]["cat"] != "-":
                                    tcDict[accession][KEY] = "1"
                            else:
                                KEY = ls[4]
                                if keyDict[KEY]["cat"] != "-":
                                    tcDict[accession][KEY] = "1"


catDict4 = defaultdict(lambda: defaultdict(list))
for i in catDict3.keys():
    for j in catDict3[i]:
        if keyDict[j]["cat"] != "-":
            catDict4[i][keyDict[j]["cat"]].append(j)

gff = open("GFF_results/tmrna.gff")
accession = ""
for i in gff:
    if re.match(r'GC', i):
        accession = (i.rstrip().split(".")[0])
    else:
        ls = i.rstrip().split("\t")
        if len(ls) > 2:
            tcDict[accession]["ssrA"] = "1"
            ncDict[accession]["ssrA"] = "1"
            allDict[accession]["ssrA"] = "1"

gff = open("GFF_results/rrna.gff")
accession = ""
accessions = []
for i in gff:
    if re.match(r'GC', i):
        accession = (i.rstrip().split(".")[0])
        if accession not in accessions:
            accessions.append(accession)
    else:
        ls = i.rstrip().split("\t")
        if len(ls) > 2:
            if re.findall(r'rrs', ls[8]) or re.findall(r'16S', ls[8]):
                tcDict[accession]["16S"] = "1"
                ncDict[accession]["16S"] = "1"
                allDict[accession]["16S"] = "1"

            if re.findall(r'rrl', ls[8]) or re.findall(r'23S', ls[8]):
                tcDict[accession]["23S"] = "1"
                ncDict[accession]["23S"] = "1"
                allDict[accession]["23S"] = "1"

            if re.findall(r'rrf', ls[8]) or re.findall(r'5S', ls[8]):
                tcDict[accession]["5S"] = "1"
                ncDict[accession]["5S"] = "1"
                allDict[accession]["5S"] = "1"

catDict5 = defaultdict(list)
catDict5["envelope"] = ["fatty_acid_biosynthesis", "phospholipid_biosynthesis", "PG_CellWall_division", "bam_complex", "sec_translocon"]
catDict5["PMF"] = ["ndh_dehydrogenase", "terminal_ox", "atp_synthase"]
catDict5["replication"] = ["holoenzyme", "replication"]
catDict5["transcription"] = ["RNApol", "sigma", "termination"]
catDict5["translation"] = ["tRNA_modification", "Aminoacylation", "Elongation_factor", "Release_factor", "Ribosome", "tmRNA"]
replication = [""]

out = open("final_summary.csv", "w")
out.write("accession,sp,section,category,hmm,gene,tc,nc,all,tcCol,ncCol,allCol\n")
for genome in accessions:
    for i in catDict5.keys():
        for j in catDict5[i]:

            if i == "PMF" and j == "terminal_ox":
                section = "PMF"
                category = "terminal_ox"

                tcList = []
                for k in catDict4[i][j]:
                    if tcDict[genome][k] == "1":
                        tcList.append(k)

                ncList = []
                for k in catDict4[i][j]:
                    if ncDict[genome][k] == "1":
                        ncList.append(k)

                allList = []
                for k in catDict4[i][j]:
                    if allDict[genome][k] == "1":
                        allList.append(k)

                tc = 0
                nc = 0
                all = 0
                tcCol = "missing"
                ncCol = "missing"
                allCol = "missing"
                if OxCheck(tcList) == 1:
                    tc = 1
                    tcCol = "terminal_ox"
                if OxCheck(ncList) == 1:
                    nc = 1
                    ncCol = "terminal_ox"
                if OxCheck(allList) == 1:
                    all = 1
                    allCol = "terminal_ox"

                out.write(
                    genome + "," + str(ncbiDict[genome]) + "," + section + "," + category + "," + "oxidase" + "," + "oxidase" +
                    "," + str(tc) + "," + str(nc) + "," + str(all) + "," + str(tcCol) + "," + str(ncCol) + "," + str(allCol) + "\n")

            elif i == "PMF" and j == "ndh_dehydrogenase":
                section = "PMF"
                category = "ndh_dehydrogenase"

                tcList = []
                for k in catDict4[i][j]:
                    if tcDict[genome][k] == "1":
                        tcList.append(k)

                ncList = []
                for k in catDict4[i][j]:
                    if ncDict[genome][k] == "1":
                        ncList.append(k)

                allList = []
                for k in catDict4[i][j]:
                    if allDict[genome][k] == "1":
                        allList.append(k)

                tc = 0
                nc = 0
                all = 0
                tcCol = "missing"
                ncCol = "missing"
                allCol = "missing"
                if NdhCheck(tcList) == 1:
                    tc = 1
                    tcCol = "ndh_dehydrogenase"
                if NdhCheck(ncList) == 1:
                    nc = 1
                    ncCol = "ndh_dehydrogenase"
                if NdhCheck(allList) == 1:
                    all = 1
                    allCol = "ndh_dehydrogenase"

                out.write(
                    genome + "," + str(ncbiDict[genome]) + "," + section + "," + category + "," + "dehydrogenase" + "," + "dehydrogenase" +
                    "," + str(tc) + "," + str(nc) + "," + str(all) + "," + str(tcCol) + "," + str(ncCol) + "," + str(allCol) + "\n")

            else:
                for k in catDict4[i][j]:
                    tc = tcDict[genome][k]
                    nc = ncDict[genome][k]
                    all = allDict[genome][k]

                    if tc == "0":
                        tcCol = "missing"
                    else:
                        tcCol = j

                    if nc == "0":
                        ncCol = "missing"
                    else:
                        ncCol = j

                    if all == "0":
                        allCol = "missing"
                    else:
                        allCol = j
                    if k.split(".")[0] not in ["TIGR03168", "TIGR01990", "TIGR00749", "TIGR01391",
                                   "TIGR02071", "TIGR02614", "TIGR03703", "TIGR01120",
                                   "TIGR01307", "TIGR00016", "TIGR00651", "TIGR00759",
                                   "TIGR01348", "PF01206", "TIGR03010", "TIGR03011", "TIGR03012",
                                   "TIGR03342", "TIGR00568", "TIGR00571", "TIGR00618", "TIGR00618",
                                   "TIGR00634", "TIGR01074", "TIGR01075", "TIGR01389", "TIGR02248",
                                   "PF00817", "TIGR00664", "TIGR00678", "TIGR03420", "TIGR04418", "TIGR00690",
                                   "TIGR02394", "TIGR02395", "TIGR02479", "TIGR02939", "TIGR01955", "TIGR00094", "TIGR00113", "TIGR00452",
                                   "TIGR00740", "TIGR00742", "TIGR00731", "TIGR01029", "TIGR00452",
                                   "TIGR00448", "TIGR00276", "TIGR00063", "TIGR00057", "TIGR00197",
                                   "TIGR00462", "TIGR03821", "TIGR00575", "TIGR03197", "TIGR03150", "TIGR02074"]:

                        if k in ["TIGR02013", "TIGR02027"]:
                            out.write(genome + "," + str(ncbiDict[genome]) + "," + i + "," + j + "," + k + "," + str(keyDict[k]["gene"]) +
                                  "," + str(tc) + "," + str("1") + "," + str(all) + "," + str(tcCol) + "," + str("RNApol") + "," + str(allCol) + "\n")
                        else:
                            out.write(genome + "," + str(ncbiDict[genome]) + "," + i + "," + j + "," + k + "," + str(
                                keyDict[k]["gene"]) +
                                      "," + str(tc) + "," + str(nc) + "," + str(all) + "," + str(tcCol) + "," + str(
                                ncCol) + "," + str(allCol) + "\n")
out.close()


#########################################################################################
### CHAPTER 2
ncbiDict = defaultdict(lambda: 'EMPTY')
ncbi = open("ncbi_assembly_info.1.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.2.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.3.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

ncbi = open("ncbi_assembly_info.4.tsv")
for i in ncbi:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ncbiDict[ls[0].split(".")[0]] = ls[7]

# TBLOUT IS THE DIRECTORY CONTAINING RAW OUTPUT FROM HMMER
tblDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '-')))
tblout = os.listdir("tblout")
for i in tblout:
    if re.findall(r'hmm', i):
        accession = i.split("hmmhits_")[1].split(".")[0]
        file = open(f"tblout/{i}/PF01653.tblout")
        match = 0
        for j in file:
            if not re.match(r'#', j):
                ls = delim(j.rstrip())
                locus = ls[0]
                hmm = ls[3]
                e = float(ls[4])
                bit = float(ls[5])
                if e < 0.0000000001:
                    match += 1

        if match > 0:
            tblDict[accession]["ligA"]["presence"] = "1"
            tblDict[accession]["ligA"]["color"] = "replication"
        else:
            tblDict[accession]["ligA"]["presence"] = "0"
            tblDict[accession]["ligA"]["color"] = "missing"


shareDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '-')))
shared = open("presence-absence-matrix-2.csv")
for i in shared:
    ls = i.rstrip().split(",")
    shareDict[ls[0].split(".")[0]][ls[4]]["presence"] = ls[6]
    shareDict[ls[0].split(".")[0]][ls[4]]["col"] = ls[7]

count = 0
geneDict = defaultdict(lambda: defaultdict(lambda: '-'))
acylDict = defaultdict(lambda: defaultdict(lambda: '-'))
out = open("final_summary.aminoacyl_updated.csv", "w")
final = open("final_summary.csv")
for i in final:
    ls = i.rstrip().split(",")
    if ls[1] == "sp":
        out.write(i.rstrip() + "\n")
    else:
        if ls[3] in ["Aminoacylation"]:
            if ls[4] in shareDict[ls[0]]:
                presence = (shareDict[ls[0]][ls[4]]["presence"])
                col = (shareDict[ls[0]][ls[4]]["col"])
                out.write(",".join(ls[0:6]) + "," + presence + "," + col + "\n")
            else:
                acylDict[ls[0]][ls[5]] = ls[7]
                out.write(",".join(ls[0:6]) + "," + ls[7] + "," + ls[10] + "\n")

        elif ls[3] in ["Release_factor"]:
            if ls[4] in shareDict[ls[0]]:
                presence = (shareDict[ls[0]][ls[4]]["presence"])
                col = (shareDict[ls[0]][ls[4]]["col"])
                out.write(",".join(ls[0:6]) + "," + presence + "," + col + "\n")
            else:
                out.write(",".join(ls[0:6]) + "," + ls[7] + "," + ls[10] + "\n")

        elif ls[3] in ["sigma"]:
            if ls[4] in shareDict[ls[0]]:
                presence = (shareDict[ls[0]][ls[4]]["presence"])
                col = (shareDict[ls[0]][ls[4]]["col"])
                out.write(",".join(ls[0:6]) + "," + presence + "," + col + "\n")
            else:
                out.write(",".join(ls[0:6]) + "," + ls[7] + "," + ls[10] + "\n")

        elif ls[2] in ["envelope"]:
            if ls[5] in ["secF", "secD"]:
                geneDict[ls[0]][ls[5]] = ls[7]
                out.write(",".join(ls[0:6]) + "," + ls[8] + "," + ls[11] + "\n")
            elif ls[5] in ["pgpA", "pgpB"]:
                geneDict[ls[0]][ls[5]] = ls[7]
                out.write(",".join(ls[0:6]) + "," + ls[8] + "," + ls[11] + "\n")
            else:
                if ls[4] in shareDict[ls[0]]:
                    presence = (shareDict[ls[0]][ls[4]]["presence"])
                    col = (shareDict[ls[0]][ls[4]]["col"])
                    out.write(",".join(ls[0:6]) + "," + presence + "," + col + "\n")
                else:
                    out.write(",".join(ls[0:6]) + "," + ls[8] + "," + ls[11] + "\n")

        elif ls[2] in ["transcription"]:
            if ls[4] in shareDict[ls[0]]:
                presence = (shareDict[ls[0]][ls[4]]["presence"])
                col = (shareDict[ls[0]][ls[4]]["col"])
                out.write(",".join(ls[0:6]) + "," + presence + "," + col + "\n")
            else:
                out.write(",".join(ls[0:6]) + "," + ls[8] + "," + ls[11] + "\n")

        else:
            out.write(",".join(ls[0:6]) + "," + ls[8] + "," + ls[11] + "\n")
out.close()

acylDict2 = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '0')))
for i in acylDict.keys():
    if "1" in [acylDict[i]["proS2"], acylDict[i]["proS"]]:
        acylDict2[i]["proS"]["presence"] = "1"
        acylDict2[i]["proS"]["color"] = "Aminoacylation"
    else:
        acylDict2[i]["proS"]["presence"] = "0"
        acylDict2[i]["proS"]["color"] = "missing"

    if "1" in [acylDict[i]["asnS"], acylDict[i]["glnS"]]:
        acylDict2[i]["glnS-asnS"]["presence"] = "1"
        acylDict2[i]["glnS-asnS"]["color"] = "Aminoacylation"
    else:
        acylDict2[i]["glnS-asnS"]["presence"] = "0"
        acylDict2[i]["glnS-asnS"]["color"] = "missing"

geneDict2 = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '-')))
for i in geneDict.keys():
    if "1" in [geneDict[i]["secD"], geneDict[i]["secF"]]:
        geneDict2[i]["secDF"]["presence"] = "1"
        geneDict2[i]["secDF"]["color"] = "sec_translocon"
    else:
        geneDict2[i]["secDF"]["presence"] = "0"
        geneDict2[i]["secDF"]["color"] = "missing"

    if "1" in [geneDict[i]["pgpA"], geneDict[i]["pgpB"]]:
        geneDict2[i]["pgpAB"]["presence"] = "1"
        geneDict2[i]["pgpAB"]["color"] = "phospholipid_biosynthesis"
    else:
        geneDict2[i]["pgpAB"]["presence"] = "0"
        geneDict2[i]["pgpAB"]["color"] = "missing"

dupDict = defaultdict(lambda: defaultdict(lambda: '-'))
out = open("final_summary.aminoacyl_updated2.csv", "w")
final = open("final_summary.aminoacyl_updated.csv")
for i in final:
    ls = i.rstrip().split(",")
    if ls[1] == "sp":
        out.write(i.rstrip() + "\n")
    else:

        if ls[5] in ["proS2", "proS"]:
            if "proS" not in dupDict[ls[0]]:
                out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + "proS" + "," + "proS" +"," +
                          str(acylDict2[ls[0]]["proS"]["presence"]) + "," + str(acylDict2[ls[0]]["proS"]["color"]) + "\n")
                dupDict[ls[0]]["proS"] = "done"

        elif ls[5] in ["secF", "secD"]:
            if "secDF" not in dupDict[ls[0]]:
                out.write(
                    ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + "secDF" + "," + "secDF" + "," + str(
                            geneDict2[ls[0]]["secDF"]["presence"]) + "," + str(geneDict2[ls[0]]["secDF"]["color"]) + "\n")
                dupDict[ls[0]]["secDF"] = "done"

        elif ls[5] in ["pgpA", "pgpB"]:
            if "pgpAB" not in dupDict[ls[0]]:
                out.write(
                    ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + "pgpAB" + "," + "pgpAB" + "," + str(
                        geneDict2[ls[0]]["pgpAB"]["presence"]) + "," + str(geneDict2[ls[0]]["pgpAB"]["color"]) + "\n")
                dupDict[ls[0]]["pgpAB"] = "done"

        elif ls[5] == "ligA":
            presence = tblDict[ls[0]][ls[5]]["presence"]
            color = tblDict[ls[0]][ls[5]]["color"]

            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + "," + "DNA_ligase_aden" + "," +
                      "DNA_ligase_aden" + "," + str(presence) + "," + str(color) + "\n")
        else:
            out.write(i.rstrip() + "\n")
out.close()

riboDict = defaultdict(lambda: defaultdict(lambda: '0'))
ribosomes = open("ribosomal_proteins.txt")
for i in ribosomes:
    ls = i.rstrip().split(" ")
    if len(ls) > 1:
        riboDict[ls[0].split(".")[0]][ls[1]] = "1"

trnaDict = defaultdict(lambda: defaultdict(lambda: "0"))
indexDict = defaultdict(lambda: "0")
accession = ""
trna = open("tRNAs_replaced.txt")
for i in trna:
    ls = i.rstrip().split("\t")
    accession = ls[0].split(".")[0]
    if ls[1] == "OrganismName":
        count = 0
        for j in ls[2:]:
            indexDict[count] = j
            count += 1
    else:
        count = 0
        for j in ls[2:]:
            aa = indexDict[count]
            # print(str(count) + "\t\t" + str(aa) + "\t\t" + j)
            trnaDict[accession][aa] = j
            count += 1
        # print("\n\n")


out = open("final_summary.aminoacyl_updated2.ribo_updated.csv", "w")
final = open("final_summary.aminoacyl_updated2.csv")
for i in final:
    ls = i.rstrip().split(",")
    if ls[3] == "Ribosome" and ls[4] not in ["23S", "16S", "5S"]:
        if ls[5] == "L7":
            RP = "L7/L12"
        else:
            RP = ls[5]

        if ls[6] != "1":
            if riboDict[ls[0]][RP] == "1":
                presence = "1"
                color = "Ribosome"
            else:
                presence = "0"
                color = "missing"

            out.write(",".join(ls[0:6]) + "," + presence + "," + color + "\n")
        else:
            out.write(i.rstrip() + "\n")

    elif ls[4] in ["23S", "16S", "5S"]:
        out.write(",".join(ls[0:6]) + ",1,rRNA\n")

    elif ls[5] in ["ssrA", "smpB"]:
        out.write(",".join(ls[0:6]) + ",1,tmRNA\n")

    else:
        out.write(i.rstrip() + "\n")

for i in trnaDict.keys():
    for j in trnaDict[i]:
        if trnaDict[i][j] == "1":
            presence = "1"
            color = "tRNA"
        else:
            presence = "0"
            color = "missing"
        out.write(i + "," + str(ncbiDict[i]) + ",translation,Ribosome," + j + "," + j + "," + presence + "," + color + "\n")

out.close()







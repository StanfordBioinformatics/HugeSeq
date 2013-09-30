#!/bin/env python
#usage: $0 <sampleID, e.g. LP6005243-DNA_A08>

import sys

# Gets header from the original VCF file
def writeHeader(fileIn, fileOut):
    line=fileIn.readline()
    while line and line.startswith("#"):
        fileOut.write(line)
        if line.startswith("#CHROM"): return None
        else: line=fileIn.readline()

# Gets next valid line from file and splits by tab.
def getNextLineValues(file):
    line = file.readline()
    while line and line.startswith("#"):
        line = file.readline()
    if line and len(line.strip()) > 0: return line.split("\t")
    return None

# Converts chr value to an interger
# Example input: "chr1" or  "chrX"
# Example output: 1 or 23
def getChrNumber(chrStr):
    val = chrStr[3:]
    if (val.isdigit()):
        return int(val)
    elif val.lower() == "x":
        return 23
    elif val.lower() == "y":
        return 24
    elif val.lower() == "m":
        return 25
    else:
        # something wrong
        return 26

# Generates a comparable integer key based on chr value and position value.
# Example input: ["chr1","10150", ..,]
# Example output: 1000010150
def getKey(values):
    chr = getChrNumber(values[0])
    pos = int(values[1])
    return chr * 1000000000 + pos

# Returns remaining lines of a file
def getRemainingLines(file, output_file):
    line = file.readline()
    while line:
        output_file.write(line)
        line = file.readline()
    return None

# Main
def main():

    print(sys.argv[1])
    
    fileA = open(sys.argv[2], "r")
    fileB = open(sys.argv[3], "r")
    fileCombine = open(sys.argv[1], "w")

    writeHeader(fileA, fileCombine)

    valuesA = getNextLineValues(fileA)
    valuesB = getNextLineValues(fileB)

    while valuesA and valuesB:
        keyA = getKey(valuesA)
        keyB = getKey(valuesB)
        if keyA < keyB:
            newline = "\t".join(valuesA)
            fileCombine.write(newline)
            valuesA = getNextLineValues(fileA)
        elif keyA > keyB:
            newline = "\t".join(valuesB)
            fileCombine.write(newline)
            valuesB = getNextLineValues(fileB)
        elif keyA == keyB:
            newlineA = "\t".join(valuesA)
            fileCombine.write(newlineA)
            newlineB = "\t".join(valuesB)
            fileCombine.write(newlineB)
            valuesA = getNextLineValues(fileA)
            valuesB = getNextLineValues(fileB)

    if valuesA:
        newline = "\t".join(valuesA)
        fileCombine.write(newline)
        getRemainingLines(fileA, fileCombine)
    elif valuesB:
        newline = "\t".join(valuesB)
        fileCombine.write(newline)
        getRemainingLines(fileB, fileCombine)


    for i in [fileA, fileB, fileCombine]:
        i.flush()
        i.close()

main()

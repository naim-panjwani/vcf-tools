#!/usr/local/anaconda3/bin/python3

import pandas as pd
import numpy as np
import gzip
import sys
from tqdm import tqdm
from datetime import datetime
import subprocess
import argparse

###### Helper functions ######
def getNumHeaderLines(vcf_filename, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    with gzip.open(vcf_filename, 'rb') as f:
        nextline = f.readline().decode('utf-8')
        num_lines_checked += 1
        while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
            num_header_lines += 1
            nextline = f.readline().decode('utf-8')
            num_lines_checked += 1
    return num_header_lines

def readHeader(infile):
    header = ""
    numHeaderLines = getNumHeaderLines(infile)
    with gzip.open(infile, 'rb') as f:
        nextline = f.readline().decode('utf-8')
        while nextline[0:1] == "#":
            header += nextline
            nextline = f.readline().decode('utf-8')
    return header

def readCHROMLine(infile):
    header = readHeader(infile)
    chromLine = ""
    for line in header.split('\n'):
        firstWord = line.split('\t')[0]
        if firstWord == "#CHROM":
            chromLine = line
    return chromLine

def updateID(infile, IDupdateFilename):
    IDupdateFile = pd.read_csv(IDupdateFilename, sep="\t")
    chromLine = readCHROMLine(infile)
    fullheader = readHeader(vcfFilename)
    fullheaderlist = [x for x in fullheader.split('\n') if x!='']    
    fullheaderlist.remove(chromLine)
    newIDs = ''
    for ID in chromLine.split('\t')[9:]:
        oldFID = ID.split('_')[0]
        oldIID = ID.split('_')[1]
        newID = '_'.join(IDupdateFile.loc[ (IDupdateFile['oldFID'] == oldFID) & (IDupdateFile['oldIID'] == oldIID) ][[ 'newFID', 'newIID' ]].values.tolist()[0])
        newIDs += '\t' + newID
    newChromLine = '\t'.join(chromLine.split('\t')[0:9]) + newIDs
    fullheaderlist.append(newChromLine)
    return '\n'.join(fullheaderlist) + '\n'

def writeHeader(headerText, outfilename):
    with open(outfilename, 'w') as f_out:
        f_out.write(headerText)

def appendVCFbody(infile, outfilename):
    numHeaderLines = getNumHeaderLines(infile)
    vcfchunks = pd.read_csv(infile, compression='gzip', sep='\t', skiprows=numHeaderLines, header=0, chunksize=500000)
    for vcfchunk in vcfchunks:
        vcfchunk.to_csv(outfilename, sep='\t', index=False, mode='a', encoding='utf-8', header=False)


#######################
###### MAIN
######################
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Update header with new IDs')
    parser.add_argument('vcfFilename', help="VCF filename")
    parser.add_argument('IDupdateFilename', help="A tsv file with columns oldFID, oldIID, newFID, newIID")
    parser.add_argument('outfilename', help="Output filename")
    args = parser.parse_args()

    #IDupdateFilename = "03-ID_update_to_BIOJUME_ID.txt"
    vcfFilename = args.vcfFilename
    IDupdateFilename = args.IDupdateFilename
    outfilename = args.outfilename
#    tmp = datetime.now().strftime("%H%M%S%f-")

#    outfile = tmp + vcfFilename.replace('.gz','')
    outfile = outfilename.replace('.gz','')
    newHeader = updateID(vcfFilename, IDupdateFilename)
    writeHeader(newHeader, outfile)
    appendVCFbody(vcfFilename, outfile)
    subprocess.run(args=['bgzip', outfile])
    subprocess.run(args=['tabix', '-p', 'vcf', outfile+'.gz'])
#    newVCFname = '03-' + vcfFilename.replace('02-','').replace('.vcf.gz','') + '_id_update.vcf.gz'
#    subprocess.run(args=['mv', f'{outfile}.gz', newVCFname])
#    subprocess.run(args=['mv', f'{outfile}.gz.tbi', f'{newVCFname}.tbi'])



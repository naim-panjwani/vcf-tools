#!/usr/bin/python3.6

import numpy as np
import pandas as pd
import gzip
import subprocess, os, sys
import argparse
from datetime import datetime
from tqdm import tqdm


###########################################################################
###### Helper functions #######
###########################################################################
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

def writeHeader(headerText, outfilename):
    with open(outfilename, 'w') as f_out:
        f_out.write(headerText)
###########################################################################


#######################
###### MAIN
######################
#vcf_filename = "03-JME_Round1_and_2_chrMT_dosage.vcf.gz"
#samples_to_extract_filename = "04-samples_to_extract_from_dosage.txt"
#outfilename = "05-JME_Round1_and_2_chrMT_dosage_BIS_only"
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Subset VCF file to the desired samples given')
    parser.add_argument('vcf_filename', help='VCF filename')
    parser.add_argument('samples_to_extract_filename', help='File with the sample names to extract')
    parser.add_argument('outfilename', help="Desired output filename (without the .vcf.gz extension)")
    args = parser.parse_args()

    vcf_filename = args.vcf_filename
    samples_to_extract_filename = args.samples_to_extract_filename
    outfilename = args.outfilename.replace('.vcf.gz','')
    logfilename = outfilename+'.log'

    old_stdout = sys.stdout
    log_file = open(logfilename, "w")
    sys.stdout = log_file

    print(datetime.now().strftime('%c'))
    print('vcf_filename: '+vcf_filename)
    print('samples_to_extract_filename: '+samples_to_extract_filename)
    print('outfilename: '+outfilename)
    print('logfilename: '+logfilename)

    print("Determining columns to subset from VCF file")
    numHeaderLines = getNumHeaderLines(vcf_filename)
    samples_to_subset = pd.read_csv(samples_to_extract_filename, sep="\t", header=None)
    samples_to_subset = list(samples_to_subset.iloc[:,0])
    chromLine = readCHROMLine(vcf_filename).split('\t')
    columns_to_subset = list(np.arange(9))
    columns_to_subset.extend([i for i,e in enumerate(chromLine) if e in samples_to_subset])
    print("Columns to subset:")
    print(columns_to_subset)

    print('Writing new VCF header')
    newHeader = '\t'.join(list(pd.Series(chromLine).iloc[columns_to_subset]))
    newHeader += '\n'
    writeHeader(newHeader, outfilename+'.vcf')
    print('Loading VCF file')
    vcfchunks = pd.read_csv(vcf_filename, sep="\t", compression="gzip", skiprows=numHeaderLines, header=0, chunksize=500000)
    for vcfchunk in tqdm(vcfchunks):
        subsetted_vcf = vcfchunk.iloc[:, columns_to_subset] 
        subsetted_vcf.to_csv(outfilename+'.vcf', sep="\t", index=False, mode='a', encoding='utf-8', header=False)
    print('Compressing subsetted VCF')
    subprocess.run(args=['bgzip', outfilename+'.vcf'])
    print('Tabix indexing subsetted VCF')
    subprocess.run(args=['tabix', '-p', 'vcf', outfilename+'.vcf.gz'])
    num_samples = subsetted_vcf.shape[1] - 9
    total_samples = len(chromLine) - 9
    print(f'Subsetting done. {num_samples}/{total_samples} samples were extracted')
    print(datetime.now().strftime('%c'))

    sys.stdout = old_stdout
    log_file.close()



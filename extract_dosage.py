#!/usr/local/anaconda3/bin/python3

import pandas as pd
import numpy as np
import argparse
import subprocess, os, sys
import gzip
from tqdm import tqdm
from datetime import datetime

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
    if num_lines_checked >= num_lines_to_check:
        raise Exception(f'Number of header lines seems to exceed {num_lines_to_check}. Could this be right?')
    return num_header_lines

def writeHeader(infile, outfilename):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfilename, 'wb') as f_out:
            nextline = f_in.readline().decode('utf-8')
            while nextline[0:1] == "#":
                f_out.write(bytes(nextline, 'UTF-8'))
                nextline = f_in.readline().decode('utf-8')

def writeCHROMline(infile, outfilename):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfilename, 'wb') as f_out:
            nextline = f_in.readline().decode('utf-8')
            while nextline.split('\t')[0] != "#CHROM" or nextline.split('\t')[0] != "CHROM":
                nextline = f_in.readline().decode('utf-8')
            f_out.write(bytes(nextline, 'UTF-8'))


if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Extract dosage field from BEAGLE 5 imputed file')
  parser.add_argument('inputvcf', help='The imputed VCF file to extract the DS field from')
  parser.add_argument('outvcf', help='The output VCF file name')
  args = parser.parse_args()

  print(datetime.now().strftime('%c'))
  print('args.inputvcf: ' + args.inputvcf.replace('.gz','').replace('.vcf','')+'.vcf.gz')
  print('args.outvcf: ' + args.outvcf.replace('.gz','').replace('.vcf','')+'.vcf')
  inputvcf = args.inputvcf
  outvcf = args.outvcf.replace('.gz','').replace('.vcf','')+'.vcf'
  logfilename = outvcf.replace('.gz','').replace('.vcf','')+'.log'

  old_stdout = sys.stdout
  log_file = open(logfilename, "w")
  sys.stdout = log_file

  print('Reading file')
  numHeaderLines = getNumHeaderLines(inputvcf)
  vcf = pd.read_csv(inputvcf, skiprows=numHeaderLines, compression='gzip', encoding='utf-8', sep="\t", chunksize=100000)
  print('Extracting dosage')
  writeHeader(inputvcf, outvcf)
  for vcfchunk in tqdm(vcf):
    vcfformatcol = vcfchunk.iloc[:,8].replace('GT:DS:GP','DS')
    vcfgenos = vcfchunk.iloc[:,9:]
    for i in np.arange(vcfgenos.shape[1]):
      vcfgenos.iloc[:,i] = [x[1] for x in list(vcfgenos.iloc[:,i].str.split(':'))]
    newvcf = pd.concat([ vcfchunk.iloc[:,0:8], vcfformatcol, vcfgenos ], axis=1)
    #print('Writing to file')
    # newvcf.to_csv(outvcf, index=False, sep="\t", compression='gzip', encoding='utf-8', mode='a', header=False) # Doesn't work on CentOS 7
    newvcf.to_csv(outvcf, index=False, sep="\t", encoding='utf-8', mode='a', header=False)
  print('Compressing')
  subprocess.run(args=['bgzip', outvcf])
  subprocess.run(args=['tabix', '-p', 'vcf', outvcf+'.gz'])
  print('Done')
  print(datetime.now().strftime('%c'))
  sys.stdout = old_stdout
  log_file.close()



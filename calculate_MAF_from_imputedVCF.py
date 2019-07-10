#!/usr/local/bin/python3.3
# Usage: python calculate_MAF_from_imputedVCF.py -i <inputfile> -o <outputfile>
# PRE: inputfile must be gzipped (or bgzipped) and must be a vcf with hard called genotypes colon and then dosage such as 0|1:0.1:*
#      That is, GT:DS:*
#      If multi-allelic (separated by commas on the same row), it will calculate all the allele frequencies and return the highest one (ie. second major allele)
# POST: A file with the first 9 columns of the VCF plus the MAF from the hard called genotypes and the MAF from dosages

import sys, getopt, gzip

def convertToMAF(AF):
  if(AF>0.5):
    return 1-AF
  else:
    return AF

def getDosage(astring, allele):
  astring2 = astring.split(':')[1]
  if(',' in astring2):
    dosages = astring2.split(',')
    return(float(dosages[allele]))
  else:
    return(float(astring2))

def getGenos(astring):
  geno1 = int(astring.split(':')[0][0])
  geno2 = int(astring.split(':')[0][2])
  return geno1, geno2

def isMultiallelic(individuals_cols):
  astring = individuals_cols[0]
  astring2 = astring.split(':')[1]
  if(',' in astring2):
    return(True)
  return(False)

def getMaxAlleles(individuals_cols):
  astring = individuals_cols[0]
  astring2 = astring.split(':')[1]
  if(',' in astring2):
    return(len(astring2.split(',')) + 1)
  else:
    return(2)

def getAlleleCount(individuals_cols, allele):
  numAlleles = 0
  for individual in individuals_cols:
    geno1, geno2 = getGenos(individual)
    if(geno1 == allele): numAlleles += 1
    if(geno2 == allele): numAlleles += 1
  return(numAlleles)

def getDosageSums(individuals_cols, allele):
  dosageSum = 0
  for individual in individuals_cols:
    dosageSum += getDosage(individual, allele)
  return(dosageSum)

def main():
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError as err:
      print(str(err))
      print('calculate_MAF_from_imputedVCF.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('calculate_MAF_from_imputedVCF.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      else:
        assert False, "unhandled option"
   print('Input file is ', inputfile)
   print('Output file is ', outputfile)
   with gzip.open(inputfile,'rb') as f_in:
     with open(outputfile, 'wb') as f_out:
       for line in f_in:
         line = line.decode('utf-8')
         if(line[0:2] != "##"):
           columns = line.split('\t')
           columns = columns[0:8]
           if(line[0] == "#"):
             columns.append("MAF_dosage")
             columns.append("MAF_geno")
           else:
             individuals_columns = line.split('\t')[9:]
             dosageSum = 0
             genosSum = 0
             numIndividuals = len(individuals_columns)
             if(isMultiallelic(individuals_columns)):
               numAlleles = getMaxAlleles(individuals_columns)
               alleleCounts = list()
               for i in range(0,numAlleles):
                 alleleCounts.append(getAlleleCount(individuals_columns, i))
               del alleleCounts[alleleCounts.index(max(alleleCounts))] # removing the major allele count
               secondMajorAlleleIndex = alleleCounts.index(max(alleleCounts))
               dosageSums = getDosageSums(individuals_columns, secondMajorAlleleIndex)
               MAF_geno = convertToMAF(alleleCounts[secondMajorAlleleIndex] / (2*numIndividuals))
               MAF_dosage = convertToMAF(dosageSums / (2*numIndividuals))
               columns.append(str(MAF_dosage))
               columns.append(str(MAF_geno))
             else:
               for individual in individuals_columns:
                 dosage = getDosage(individual, 0)
                 geno1, geno2 = getGenos(individual)
                 dosageSum += dosage
                 genosSum = genosSum + geno1 + geno2
               MAF_dosage = convertToMAF(dosageSum / (numIndividuals*2))
               MAF_geno = convertToMAF(genosSum / (numIndividuals*2))
               columns.append(str(MAF_dosage))
               columns.append(str(MAF_geno))
           newline = '\t'.join(columns)
           newline += '\n'
           f_out.write(bytes(newline, 'UTF-8'))


if __name__ == '__main__':
  main()


#!/usr/local/bin/python3.3

import sys, getopt, gzip

def convertToMAF(AF):
  if(AF>0.5):
    return 1-AF
  else:
    return AF

def main():
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError as err:
      print(str(err))
      print('select_most_common_multiallelic.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('select_most_common_multiallelic.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      else:
        assert False, "unhandled option"
   print('Input file is ', inputfile)
   print('Output file is ', outputfile)
   with open(inputfile,'rb') as f_in:
     with open(outputfile, 'wb') as f_out:
       for line in f_in:
         line = line.decode('utf-8')
         columns = line.split('\t')
         if("," in columns[4]):
           alleles = columns[4].split(',')
           AFs = columns[6].split(',')
           maf_curr = 0
           best_allele = ''
           for i in range(0,len(AFs)):
             if(convertToMAF(float(AFs[i])) > maf_curr):
               maf_curr = convertToMAF(float(AFs[i]))
               best_allele = alleles[i]
           columns[4] = best_allele
           columns[6] = str(maf_curr)
           columns[7] = str(maf_curr)
           modLine = '\t'.join(columns)
           modLine += '\n'
           f_out.write(bytes(modLine, 'UTF-8'))
         else:
           f_out.write(bytes(line,'UTF-8'))

if __name__ == '__main__':
  main()


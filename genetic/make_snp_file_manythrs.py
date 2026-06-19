
#python3 extract_interesting_rs_v3.py |grep rs| perl -pe 's/\e\[?.*?[\@-~]//g' |grep -v "\[" > myrsnums.txt
#cat myrsnums.txt|while read F; do printf "'$F',";done
import os
import sys
import scipy.io
import time

genelist = 'genes_all'
thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10]
if len(sys.argv) > 1:
  genelist = sys.argv[1]

#Load the list of all rs numbers within the set of genes
A=scipy.io.loadmat('filtered_ripke2020_v3_'+genelist+'.mat')
isignificant = [i for i in range(0,len(A['pvals'])) if float(A['pvals'][i]) < max(thrs)]
print("len(isignificant)="+str(len(isignificant))+", len(A['pvals'])="+str(len(A['pvals'])))

#Save the names of genes, and gene index, p-value and rs number for each given rs number that has p-value smaller than largest threshold in thrs
geneNamesAll = A['geneNamesAll']
igenesFound = [A['igenes'][0][i] for i in isignificant]
pvalsFound = [A['pvals'][i] for i in isignificant]
rssFound = [A['rss'][i] for i in isignificant]
rsnums = [x[0:x.find(' ')] if x.find(' ') > -1 else x for x in rssFound]

#Read the variants names to get the SNP identificator
input_file = open('../hrc.grch37.variants.csv','r')
line = input_file.readline()

iline = 0
lines_found = []
pvals_found = []
irsnums_found = []
while len(line) > 0:
  if iline%1000000 == 999999:
    print(str(iline+1)+'/40000000 done, time = '+str(time.time()))
  line = input_file.readline()
  fields = line.split('\t')
  if len(line) == 0:
    break
  if fields[2] in rsnums:
    lines_found.append(line)
    irsnum = [i for i in range(0,len(rsnums)) if rsnums[i] == fields[2]][0]
    irsnums_found.append(irsnum)
    pvals_found.append(pvalsFound[irsnum])
    print(line)
  iline = iline + 1

input_file.close()

for ithr in range(0,len(thrs)):
  mythr = thrs[ithr]
  output_file = open('snps_of_interest_'+genelist+'_'+str(mythr)+'.txt','w')
  output_file2 = open('rs_of_interest_'+genelist+'_'+str(mythr)+'.txt','w')
  #output_file3 = open('geneticBD/rsonly_of_interest_'+genelist+'_'+str(mythr)+'.txt','w')

  for iline in range(0,len(lines_found)):
    if float(pvals_found[iline]) > mythr:
      continue
    line = lines_found[iline]
    if line.find('rs') > -1:
      fields = line.split('\t')
      mystr = fields[0]+':'+fields[1]+":"+fields[3]+":"+(fields[4] if fields[4][-1] != '\n' else fields[4][0:-1])
      print(mystr)
      output_file.write(mystr+'\n')
      output_file2.write(fields[2]+' '+mystr+' '+geneNamesAll[igenesFound[irsnums_found[iline]]]+'\n')
      #output_file3.write(fields[2]+'\n')

  output_file.close()
  output_file2.close()
  #output_file3.close()

qwe

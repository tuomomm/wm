#module load SciPy-bundle/2024.05-gfbf-2024a
import scipy.io
import numpy
import sys
from os.path import exists
import os
from datetime import datetime

genelist = 'genes_all'
#imeas_str = '0' #0=dsp_r,1=dsp_r_1,2=lns_r,3=MX_lns_r,4=MX_wms_r
#cmeas_str = '1'

thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 1e-10]
if len(sys.argv) > 1:
  genelist = sys.argv[1]
#if len(sys.argv) > 3:
#  imeas_str = sys.argv[3]
#if len(sys.argv) > 4:
#  cmeas_str = sys.argv[4]

#imeas = [int(float(x)) for x in imeas_str.split(',')]
#cmeas = [float(x) for x in imeas_str.split(',')]

PRS_vecs_all = []
IDs_all = []
for ithr in range(0,len(thrs)):
  mythr = thrs[ithr]
  if not exists('geneticBD/'+genelist+'_'+str(mythr)+'_PRS.best'):
    PRS_vecs_all.append([])
    continue
  #qwe
  input_file = open('geneticBD/'+genelist+'_'+str(mythr)+'_PRS.best','r')
  line = input_file.readline()
  IDs = []
  PRSs = []
  while len(line) > 0:
    line = input_file.readline()
    if len(line) < 2:
      break
    fields = line.split(' ')
    #print(str(fields))
    IDs.append(fields[0])
    PRSs.append(float(fields[3][0:-1])) #last character removed, it's '\n'
  input_file.close()

  IDs_all.append(IDs[:])
  PRS_vecs_all.append(PRSs[:])
PRS_vecs_all = numpy.array(PRS_vecs_all)
  
#Check that the IDs are the same in all threshold p-value files:
for ithr in range(0,min(len(thrs),len(IDs_all))):
  if any([IDs_all[ithr][i] != IDs[i] for i in range(0,len(IDs))]):
    print('ID mismatch ithr = '+str(ithr))
    qwe

#Find the el-phys IDs for each genotype ID
input_file = open('subjs_WMdata_genotypeID.txt','r')
line = 'noobnoob\n'
IDdict = {}
IDdict2 = {}
while len(line) > 0:
  line = input_file.readline().split('\n')[0]
  if len(line) < 3:
    break
  fields = line.split(' ')
  IDdict[fields[1]] = fields[0]
  IDdict2[fields[0]] = fields[1]
input_file.close()
#print(str(IDdict))

data_allsubjects = []
print('Loading TOPstudy_baseline_cogvars_4Tuomo_211024_ztransformed.csv')
input_file = open('TOPstudy_baseline_cogvars_4Tuomo_211024_ztransformed.csv')
line = input_file.readline()
while len(line) > 0:
  line = input_file.readline().split('\n')[0]
  if len(line) == 0:
    continue
  data_allsubjects.append(line.split(','))

ID_allsubjects = [data_allsubjects[i][0] for i in range(0,len(data_allsubjects))]
WMmeasures_allsubjects = [ [(float(data_allsubjects[i][8+j]) if data_allsubjects[i][8+j] not in ['NA','N/A'] else numpy.nan) for i in range(0,len(data_allsubjects))] for j in range(2,4)]
ages_allsubjects = [data_allsubjects[i][2] for i in range(0,len(data_allsubjects))]
sexes_allsubjects = [data_allsubjects[i][3] for i in range(0,len(data_allsubjects))]
diags_allsubjects = [data_allsubjects[i][4] for i in range(0,len(data_allsubjects))]
days_since_2000_allsubjects = [(datetime.strptime(data_allsubjects[i][5], "%Y-%m-%d")-datetime(2000, 1, 1)).days if 'NA' not in data_allsubjects[i][5] else numpy.nan for i in range(0,len(data_allsubjects))]
IQ_allsubjects = [data_allsubjects[i][6] for i in range(0,len(data_allsubjects))]
educat_allsubjects = [data_allsubjects[i][7] for i in range(0,len(data_allsubjects))]

IDsMissing = [] #this will remain empty
WMmeasures_vec = []
age_vec = []
sex_vec = []
group_vec = []
subject_vec = []
days_since_2000_vec = []
IQ_vec = []
educat_vec = []
for iID in range(0,len(IDs)):
  isubjs = [i for i in range(0,len(ID_allsubjects)) if ID_allsubjects[i] == IDdict[IDs[iID]]]
  if len(isubjs) == 0:
    print('No subj found')
    continue
  if len(isubjs) > 1:
    print('Many subjs found')
    continue
  isubj = isubjs[0]
  if ages_allsubjects[isubj] in ['NA','N/A'] or sexes_allsubjects[isubj] in ['NA','N/A'] or diags_allsubjects[isubj] in ['NA','N/A']:
    print('Demographics data missing for ID='+ID_allsubjects[isubj])
    continue
  WMmeasures_vec.append([WMmeasures_allsubjects[i][isubj] for i in range(0,2)])
  age_vec.append(float(ages_allsubjects[isubj]))
  sex_vec.append([1 if sexes_allsubjects[isubj] == 'male' else (0 if sexes_allsubjects[isubj] == 'female' else -1)])
  group_vec.append([1 if diags_allsubjects[isubj] == 'SZ' else (0 if diags_allsubjects[isubj] == 'HC' else -1)])
  #subject_vec.append(ID_allsubjects[isubj])
  subject_vec.append(IDs[iID])
  days_since_2000_vec.append(float(days_since_2000_allsubjects[isubj]))
  IQ_vec.append(float(IQ_allsubjects[isubj]))
  educat_vec.append(float(educat_allsubjects[isubj]) if 'NA' not in educat_allsubjects[isubj] else numpy.nan)

#IDs = [IDs[iID] for iID in range(0,len(IDs)) if iID not in IDsMissing]
#PRS_vecs_all = [PRS_vecs_all[:,iID] for iID in range(0,len(PRS_vecs_all[0,:])) if iID not in IDsMissing]

PRS_vecs_all = numpy.array(PRS_vecs_all).T
#PRS_vecs_all = [ [PRS_vecs_all[i][j] for j in range(0,len(PRS_vecs_all[0])) ] for i in range(0,len(PRS_vecs_all)) ]
print("PRS_vecs_all shape = "+str(numpy.array(PRS_vecs_all).shape))
    
scipy.io.savemat('assoc/BDPRSs_othermeasures_ztransformed_LNSonly_diags_'+genelist+'.mat', {'WMmeasures_vec': WMmeasures_vec, 'age_vec': age_vec, 'group_vec': group_vec, 'sex_vec': sex_vec, 'subjects': subject_vec, 'PRS_vecs': PRS_vecs_all, 'thrs': thrs, 'IQ_vec': IQ_vec, 'days_since_2000_vec': days_since_2000_vec, 'educat_vec': educat_vec})

qwe

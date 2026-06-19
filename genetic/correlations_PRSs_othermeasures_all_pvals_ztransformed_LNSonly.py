#cp lin_regressions3.py correlations.py : just calculate correlations between SNP data and EEG phenotypes, compare with random.
#cp lin_regressions2.py lin_regressions3.py : use only healthy controls
#cp lin_regressions_randomgenes_noremoval_compare_sklearn.py lin_regressions2.py
import sys
import numpy
import scipy.io
import time
import scipy.stats
import pandas
from sklearn import linear_model
from sklearn_extension2 import LinearRegression
from os.path import exists

imethod = 0 #0: lin regr., 1: Ridge regr., 2: Lasso regr.
group = 'HC'

if len(sys.argv) > 1:
  group = sys.argv[1]
if len(sys.argv) > 2:
  imethod = int(float(sys.argv[2]))

attrs = ['genes_all_modeled_synaptic', 'genes_full', 'genes_all_modeled_ionchannels2', 'genes_all_modeled_ionchannels2_CACNA1C', 'genes_all_modeled_ionchannels2_CACNA1D', 'genes_all_modeled_ionchannels2_CACNA1I', 'genes_all_modeled_ionchannels2_KCNB1', 'genes_all_modeled_ionchannels2_KCNJ6', 'genes_all_modeled_ionchannels2_KCND3', 'genes_all_modeled_ionchannels2_HCN1', 'genes_all_modeled_ionchannels2_ATP2A2', 'genes_new', 'genes_new_ionchannels', 'genes_new_synaptic', 'genes_new_ionchannels_kcn', 'genes_new_ionchannels_cacn', 'genes_new_ionchannels_scn_and_hcn', 'genes_new_synaptic_others', 'genes_new_synaptic_PKA', 'genes_new_synaptic_PKC', 'genes_new_synaptic_PKAPKC', 'genes_new_synaptic_PKAPKC_noPP']

corrs = []
betas_pvals = []
Ns = []
thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 1e-10]

WMmeasures = ['WAIS','MCCB']
WMmeasures_all = ['WAIS','MCCB','Combined']
attrs_saved_all = []
for iWMgroup in range(0,len(WMmeasures_all)):
  corrs_thisWMgroup = []
  PRS_pvals_thisWMgroup = []
  betas_pvals_thisWMgroup = []
  Ns_thisWMgroup = []
  attrs_saved = []
  for iattr in range(0,len(attrs)):
    myattr = attrs[iattr]
    corrs_thisattr = []
    PRSs_pvals_thisattr = []
    betas_pvals_thisattr = []
    Ns_thisattr = []
    if not exists('assoc/PRSs_othermeasures_ztransformed_LNSonly_diags_'+myattr+'.mat'):
      print('assoc/PRSs_othermeasures_ztransformed_LNSonly_diags_'+myattr+'.mat does not exist')
      continue
    MAT_WM = scipy.io.loadmat('assoc/PRSs_othermeasures_ztransformed_LNSonly_diags_'+myattr+'.mat')

    Yall = numpy.array(MAT_WM['WMmeasures_vec'])
    Yavg = numpy.array([[Yall[i,1] if numpy.isnan(Yall[i,0]) else (Yall[i,0] if numpy.isnan(Yall[i,1]) else 0.5*(Yall[i,0]+Yall[i,1])) for i in range(0,Yall.shape[0])]]).T

    Yall = numpy.c_[Yall,Yavg] #columns 0-4 separate exps, columns 5-6 similar experiments combined, columns 7-9 two experiments of different modality combined, column 10 trhee experiments combined

    if len(MAT_WM['PRS_vecs']) == 1:
      MAT_WM['PRS_vecs'] = MAT_WM['PRS_vecs'][0]
    if len(MAT_WM['group_vec']) == 1:
      MAT_WM['group_vec'] = MAT_WM['group_vec'][0]

    subjs_reject = ['FAM001_UMRK59311'] # From unrel.king.cutoff.out.id : too close relative to a psych patient
    isubj_reject = [i for i in range(0,len(MAT_WM['subjects'])) if (MAT_WM['subjects'][i] if MAT_WM['subjects'][i].find(' ') == -1 else MAT_WM['subjects'][i][0:MAT_WM['subjects'][i].find(' ')]) in subjs_reject]
    isubj_reject = isubj_reject + [i for i in range(0,len(MAT_WM['subjects'])) if numpy.isnan(MAT_WM['age_vec'][0][i])] #116 participants had undeteremined  age, sex and group, Remove these.
    isubj_reject = isubj_reject + [i for i in range(0,len(MAT_WM['subjects'])) if numpy.isnan(Yall[i,iWMgroup])] # Remove participants that don't have the requested WM index measured. This varies from 0 to hundreds of subjects

    if group == 'HC':
      isubjs = [i for i in range(0,len(MAT_WM['group_vec'])) if MAT_WM['group_vec'][i] == 0 and i not in isubj_reject]
    elif group == 'SZ':
      isubjs = [i for i in range(0,len(MAT_WM['group_vec'])) if MAT_WM['group_vec'][i] == 1 and i not in isubj_reject]
    elif group == 'both':
      isubjs = [i for i in range(0,len(MAT_WM['group_vec'])) if MAT_WM['group_vec'][i] in [0,1] and i not in isubj_reject]
    else:
      print('Group unknown!')
    print("len(isubj_reject) = "+str(len(isubj_reject))+" / "+str(len(MAT_WM['group_vec'][0])))

    if len(MAT_WM['age_vec']) == len(MAT_WM['group_vec']) and type(MAT_WM['age_vec'][0]) is numpy.ndarray:
      MAT_WM['age_vec'] = MAT_WM['age_vec'].T[0]
    ages = [float(MAT_WM['age_vec'][i] if len(MAT_WM['age_vec']) > 1 else MAT_WM['age_vec'][0][i]) for i in isubjs]
    if len(MAT_WM['sex_vec']) == len(MAT_WM['group_vec']) and type(MAT_WM['sex_vec'][0]) is numpy.ndarray:
      MAT_WM['sex_vec'] = MAT_WM['sex_vec'].T[0]
    sexes = [MAT_WM['sex_vec'][i] if len(MAT_WM['sex_vec']) > 1 else MAT_WM['sex_vec'][0][i] for i in isubjs]
    IQs = [MAT_WM['IQ_vec'][i] if len(MAT_WM['IQ_vec']) > 1 else MAT_WM['IQ_vec'][0][i] for i in isubjs]
    days = [MAT_WM['days_since_2000_vec'][i] if len(MAT_WM['days_since_2000_vec']) > 1 else MAT_WM['days_since_2000_vec'][0][i] for i in isubjs]
    educat = [MAT_WM['educat_vec'][i] if len(MAT_WM['educat_vec']) > 1 else MAT_WM['educat_vec'][0][i] for i in isubjs]
    diags = [MAT_WM['group_vec'][i] for i in isubjs]
    diags = numpy.array(diags).reshape(len(diags))
    
    Y = [Yall[i,iWMgroup] for i in isubjs]
    for ithr in range(0,len(thrs)):
      mythr = thrs[ithr]
      MAT_WM_PRS = MAT_WM['PRS_vecs']
      if len(MAT_WM_PRS) == 1 and len(MAT_WM_PRS[0]) > 1:
        MAT_WM_PRS = MAT_WM_PRS[0]
      if len(MAT_WM_PRS) != len(thrs):
        MAT_WM_PRS = MAT_WM_PRS.T
      MAT_WM_PRS = MAT_WM_PRS[ithr]
      if len(MAT_WM_PRS) == 0:
        continue
      if len(MAT_WM_PRS) == 1:
        try:
          MAT_WM_PRS = MAT_WM_PRS[0]
        except:
          print('MAT_WM_PRS = '+str(MAT_WM_PRS))
          pass
      if len(MAT_WM_PRS) < len(isubjs):
        #print("i/len(thrs)="+str(ithr)+"/"+str(len(thrs))+" continue, len(MAT_WM_PRS) = "+str(len(MAT_WM_PRS)))
        try:
          print(myattr+", len(MAT_WM_PRS[0]) = "+str(len(MAT_WM_PRS[0]))+", len(isubjs) = "+str(len(isubjs)))
        except:
          print("Problem with "+myattr+", ithr="+str(ithr))
          pass
        #continue
      #else:
      #  print(myattr+", len(MAT_WM_PRS) = "+str(len(MAT_WM_PRS))+", len(isubjs) = "+str(len(isubjs)))
      PRSs = numpy.array([MAT_WM_PRS[i] for i in isubjs])

      thiscorr = numpy.corrcoef(PRSs,Y)[0,1]

      #if ielphys == 3 or ielphys == 5 or ielphys == 9 or ielphys == 13:
      #print("Correlations "+elphystarget+" ("+group+"), "+myattr+" SNPs, thr="+str(mythr)+", corr="+'{:.4f}'.format(thiscorr))

      Xdf = pandas.DataFrame(numpy.array([PRSs,ages,sexes]).T,columns = ['PRS','Age','Sex'])
      Ydf = pandas.DataFrame(Y,columns = ["WM group "+str(iWMgroup)])
      Xdf.insert(0,'const',[1 for i in range(0,len(Y))])

      model = LinearRegression()
      model.fit(Xdf, Ydf)
      model = LinearRegression().fit(Xdf, Ydf)
      r_sq = model.score(Xdf, Ydf)
      #print(model.coef_)
      #print(model.p)
      
      corrs_thisattr.append(thiscorr)
      PRSs_pvals_thisattr.append(model.p[0][1])
      betas_pvals_thisattr.append([model.coef_[0][:],model.p[0][:]])
      Ns_thisattr.append(len(Y))
    for ithr in range(1,len(Ns_thisattr)):
      if Ns_thisattr[ithr] != Ns_thisattr[0]:
        print("Error: N mismatch!")
    if len(corrs_thisattr) == 0:
      print("Problem with "+myattr+", corrs_thisattr="+str(corrs_thisattr))
      continue
    print(  "Correlations "+"WM group "+str(iWMgroup)+" ("+group+"), "+myattr+(" "*(len('genes_lips_synaptic_ABCDEFGHIJKLMNOPQR')-len(myattr)))+\
            " SNPs, min pval="+'{:.4f}'.format(numpy.min(PRSs_pvals_thisattr))+" ("+''.join(['§' if x < 0.00357 else ('*' if x < 0.05 else '\'') for x in PRSs_pvals_thisattr])+" "*(13-len(PRSs_pvals_thisattr))+"), av corr="+'{:.4f}'.format(numpy.mean(corrs_thisattr))+" +- "+'{:.4f}'.format(numpy.std(corrs_thisattr))+" (min "+'{:.4f}'.format(numpy.min(corrs_thisattr))+" max "+'{:.4f}'.format(numpy.max(corrs_thisattr))+"; N="+str(Ns_thisattr[0])+'; corr(IQ,WM)='+'{:.2f}'.format(numpy.corrcoef(IQs,Y)[0,1])+'; corr(days,WM)='+'{:.2f}'.format(numpy.corrcoef(days,Y)[0,1])+'; corr(sex,WM)='+'{:.2f}'.format(numpy.corrcoef(sexes,Y)[0,1]))
    corrs_thisWMgroup.append(corrs_thisattr[:])
    betas_pvals_thisWMgroup.append(betas_pvals_thisattr[:])
    Ns_thisWMgroup.append(Ns_thisattr[:])
    attrs_saved.append(myattr)
    print("len(ages)="+str(len(ages))+",len(Y)="+str(len(Y))+",len(sexes)="+str(len(sexes))+",ages="+str(numpy.mean(ages))+"+-"+str(numpy.std(ages))+",Nmale="+str(sum([1 for x in sexes if x==1]))+",Nfemale="+str(sum([1 for x in sexes if x==0]))+" (#SZ = "+str(sum([1 for x in diags if x==1]))+", #HC = "+str(sum([1 for x in diags if x==0]))+")")
  print("")
  r_days, p_value_days = scipy.stats.pearsonr(days, Y)
  r_ages, p_value_ages = scipy.stats.pearsonr(ages, Y)
  r_IQs, p_value_IQs = scipy.stats.pearsonr(IQs, Y)
  r_sexes, p_value_sexes = scipy.stats.pearsonr(sexes, Y)
  r_diags, p_value_diags = scipy.stats.pearsonr(diags, Y)
  print("iWMgroup = "+str(iWMgroup)+" Correlation between WM and days with significance: "+str(numpy.corrcoef(days,Y)[0,1])+"="+str(r_days)+", p="+str(p_value_days))
  print("iWMgroup = "+str(iWMgroup)+" Correlation between WM and age with significance: "+str(numpy.corrcoef(ages,Y)[0,1])+"="+str(r_ages)+", p="+str(p_value_ages))
  print("iWMgroup = "+str(iWMgroup)+" Correlation between WM and IQ with significance: "+str(numpy.corrcoef(IQs,Y)[0,1])+"="+str(r_IQs)+", p="+str(p_value_IQs))
  print("iWMgroup = "+str(iWMgroup)+" Correlation between WM and sex with significance: "+str(numpy.corrcoef(sexes,Y)[0,1])+"="+str(r_sexes)+", p="+str(p_value_sexes))
  print("iWMgroup = "+str(iWMgroup)+" Correlation between WM and diagnosis with significance: "+str(numpy.corrcoef(diags,Y)[0,1])+"="+str(r_diags)+", p="+str(p_value_diags))
  print("iWMgroup = "+str(iWMgroup)+" p-value between WM in SCZ versus HC: "+str(scipy.stats.ranksums([Y[i] for i in range(0,len(Y)) if diags[i]==0],[Y[i] for i in range(0,len(Y)) if diags[i]==1])[1])+", N="+str(len([Y[i] for i in range(0,len(Y)) if diags[i]==0]))+","+str(len([Y[i] for i in range(0,len(Y)) if diags[i]==1])))
  print("iWMgroup = "+str(iWMgroup)+" p-value between WM in male versus female: "+str(scipy.stats.ranksums([Y[i] for i in range(0,len(Y)) if sexes[i]==0],[Y[i] for i in range(0,len(Y)) if sexes[i]==1])[1])+", N="+str(len([Y[i] for i in range(0,len(Y)) if sexes[i]==0]))+","+str(len([Y[i] for i in range(0,len(Y)) if sexes[i]==1])))
  corrs.append(corrs_thisWMgroup[:])
  betas_pvals.append(betas_pvals_thisWMgroup[:])
  Ns.append(Ns_thisWMgroup[:])
  attrs_saved_all.append(attrs_saved[:])
scipy.io.savemat('corrs_ztransformed_LNSonly_againstAgeSex_'+group+'.mat',{'corrs': corrs,'betas_pvals': betas_pvals, 'Ns': Ns, 'attrs': attrs_saved_all})
qwe

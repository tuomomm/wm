#cp drawfigcorrs.py drawfig5.py
from pylab import *
import scipy.io
import mytools

thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 1e-10]

diaggroup = 'both'

f,axs = subplots(2,1)
for iax in range(0,2):
  axs[iax].set_position([0.08+0.49*iax,0.56*1,0.42,0.37])
  axs[iax].tick_params(axis='both', which='major', labelsize=4, length=2)
  for axis in ['top','bottom','left','right']:
    axs[iax].spines[axis].set_linewidth(0.0)
  for axis in ['left']:
    axs[iax].spines[axis].set_linewidth(0.3)

axarr = axs.reshape(prod(axs.shape),).tolist()
axarr[0].set_position([0.06,       0.12+0.44+0.24,0.25,0.15])
axarr[1].set_position([0.06+0.29,  0.12+0.44+0.24,0.25,0.15])

for icohort in range(0,2):
  A = scipy.io.loadmat('../genetic/corrs_ztransformed_LNSonly_againstAgeSex_'+diaggroup+'.mat')

  PRSnames = A['attrs'][0]
  for iname in range(0,len(PRSnames)):
    ifirstspace = PRSnames[iname].find(' ')
    if ifirstspace > -1:
      PRSnames[iname] = PRSnames[iname][0:ifirstspace]

  axlabels = ['WAIS-III','MCCB']
  ithr = 9

  igroup = icohort

  PRSnames_to_consider = ['genes_full','genes_new_ionchannels','genes_new_synaptic','genes_new','genes_new_ionchannels_scn_and_hcn','genes_new_ionchannels_kcn','genes_new_ionchannels_cacn','genes_all_modeled_ionchannels2','genes_new_synaptic_PKA','genes_new_synaptic_PKC','genes_new_synaptic_others','genes_new_synaptic_PKAPKC','genes_new_synaptic_PKAPKC_noPP','genes_all_modeled_synaptic']

  cols = ['#88AA88','#6666FF','#FF2222','#FF66FF','#6666FF','#6666FF','#6666FF','#6666FF','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222']
  PRSlabel_names = ['Full genome', 'Ion channels', 'Plasticity', 'Plasticity+\nion channels', 'V.G. Na+ and\nHCN channels','V.G. K+ channels','V.G. Ca2+\nchannels (CACN*)','Only modelled\nion channels','Plasticity via PKA','Plasticity via PKC','Plasticity, others','Plasticity PKA+PKC','Plasticity PKA+PKC no PP','Plasticity, only modelled']

  edgecolors = ['#888888','#9999FF','#FF9999','#FF99FF','#9999FF','#9999FF','#9999FF','#9999FF','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999']
  hatches = [None,None,None,None,'|||','----','////','oooo','||||','----','////','++++','xxxx','oooo']


  corrs = A['corrs'][igroup]
  betas_pvals = A['betas_pvals'][igroup]
  mincorr = inf
  maxcorr = -inf
  for iPRS in range(0,len(PRSnames_to_consider)):
    iPRSname = -1
    for iPRS0 in range(0,len(PRSnames)):
      if PRSnames[iPRS0] == PRSnames_to_consider[iPRS]:
        iPRSname = iPRS0
        break
    if iPRSname == -1:
      print("Error: PRS not found, "+PRSnames_to_consider[iPRS])
      exit()
    axarr[icohort].bar(iPRS,corrs[iPRSname][0][ithr],facecolor=cols[iPRS],hatch=hatches[iPRS],edgecolor=edgecolors[iPRS],linewidth=0,label=PRSnames[iPRS])
    printZero = 0 if (corrs[iPRSname][0][ithr] < 0 and abs(corrs[iPRSname][0][ithr]) > 0.07) or corrs[iPRSname][0][ithr] > 0 else 1
    fontweight = 'bold' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else 'normal'
    print('igroup = '+str(igroup)+', '+PRSnames[iPRSname]+', ithr = '+str(ithr)+': pval = '+str(betas_pvals[iPRSname][ithr][1][1]))
    axarr[icohort].text(iPRS,min(0,corrs[iPRSname][0][ithr])-0.001,'p='+'{:.3f}'.format(betas_pvals[iPRSname][ithr][1][1])+('$\\ast\!\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05/14 else ('$\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else '')), fontsize=6,rotation=90,ha='center',va='top',fontweight=fontweight)
    maxcorr = max(maxcorr, corrs[iPRSname][0][ithr])
    mincorr = min(mincorr, corrs[iPRSname][0][ithr])
    

  axarr[icohort].set_xlim([-0.6,13.6])
  axarr[icohort].plot([-1,14],[0,0],'k-',lw=0.3)
  axarr[icohort].set_xticks([])
  axarr[icohort].set_ylim([-0.15,0.04])

  axarr[icohort].set_title(axlabels[icohort],fontsize=7)
  print(str([mincorr,maxcorr]))

for iax in range(0,2):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.01, pos.y1 +0.01, chr(ord('A')+5+iax), fontsize=11)
axarr[0].set_ylabel('Corr. coeff.',fontsize=6,labelpad=1)
axarr[1].set_yticklabels([])

PRSlabel_names = ['Full genome', 'Ion channels', 'Plasticity', 'Plasticity+\nion channels', 'V.G. Na+ and\nHCN channels','V.G. K+ channels','V.G. Ca2+\nchannels (CACN*)','Only modelled\nion channels','Plasticity via PKA','Plasticity via PKC','Plasticity, others','Plasticity PKA+PKC','Plasticity PKA+PKC no PP','Plasticity, only modelled']
labelsToPlot = ['   '+x.replace('\n','\n   ') for x in PRSlabel_names]
edgecolorsToPlot = ['#888888','#9999FF','#FF9999','#FF99FF','#9999FF','#9999FF','#9999FF','#9999FF','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999']
widthsToPlot = [0.0]*len(PRSlabel_names)
hatchesToPlot = [None,None,None,None,'|||','----','////','oooo','||||','----','////','++++','xxxx','oooo']
colsToPlot = ['#88AA88','#6666FF','#FF2222','#FF66FF','#6666FF','#6666FF','#6666FF','#6666FF','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222']

f.savefig("figS9F-G.pdf")
      

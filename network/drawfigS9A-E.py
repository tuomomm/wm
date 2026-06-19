from pylab import *
import scipy.io
import mytools

thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 1e-10]

diaggroup = 'both'
covariates = ['againstAgeSex','againstAgeSex','againstAgeSexIq','againstAge','nocovariates']
ithrs = [12, 5, 9, 9, 9]

axlabels = ['Against age\nand sex,\nSNP thresh. 1e-10','Against age\nand sex,\nSNP thresh. 1e-6','Against age,\nsex and IQ,\nSNP thresh. 1e-8','Against only age,\nSNP thresh. 1e-8','No covariates,\nSNP thresh. 1e-8']

f,axs = subplots(2,2)
for iax in range(0,2):
  for iay in range(0,2):
    axs[iay,iax].set_position([0.08+0.49*iax,0.56*(1-iay),0.42,0.37])
    axs[iay,iax].tick_params(axis='both', which='major', labelsize=4, length=2)
    for axis in ['top','bottom','left','right']:
      axs[iay,iax].spines[axis].set_linewidth(0.0)
    for axis in ['left']:
      axs[iay,iax].spines[axis].set_linewidth(0.3)

axarr = axs.reshape(prod(axs.shape),).tolist()
axarr[0].set_position([0.06,  0.12+0.44+0.23,0.14,0.11])
axarr[1].set_position([0.06+0.16,  0.12+0.44-0.1,0.77,0.15+0.24])
axarr[2].set_visible(False)
axarr[3].set_visible(False)

for icov in range(0,len(covariates)):
  A = scipy.io.loadmat('../genetic/corrs_ztransformed_LNSonly_'+covariates[icov]+'_'+diaggroup+'.mat')
  ithr = ithrs[icov]

  PRSnames = A['attrs'][0]
  for iname in range(0,len(PRSnames)):
    ifirstspace = PRSnames[iname].find(' ')
    if ifirstspace > -1:
      PRSnames[iname] = PRSnames[iname][0:ifirstspace]

  igroup = 2
  groupNames = ['WAIS','MCCB','Combined']


  PRSnames_to_consider = ['genes_full','genes_new_ionchannels','genes_new_synaptic','genes_new','genes_new_ionchannels_scn_and_hcn','genes_new_ionchannels_kcn','genes_new_ionchannels_cacn','genes_all_modeled_ionchannels2','genes_new_synaptic_PKA','genes_new_synaptic_PKC','genes_new_synaptic_others','genes_new_synaptic_PKAPKC','genes_new_synaptic_PKAPKC_noPP','genes_all_modeled_synaptic']

  cols = ['#88AA88','#6666FF','#FF2222','#FF66FF','#6666FF','#6666FF','#6666FF','#6666FF','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222','#FF2222']
  PRSlabel_names = ['Full genome', 'Ion channels', 'Plasticity', 'Plasticity+\nion channels', 'V.G. Na+ and\nHCN channels','V.G. K+ channels','V.G. Ca2+\nchannels (CACN*)','Only modelled\nion channels','Plasticity via PKA','Plasticity via PKC','Plasticity, others','Plasticity PKA+PKC','Plasticity PKA+PKC no PP','Plasticity, only modelled']
  PRSlabel_names = ['Full genome', 'Ion channels', 'Plasticity', 'Plasticity and ion channels', 'V.G. Na+ and HCN channels','V.G. K+ channels','V.G. Ca2+ channels (CACN*)','Only modelled ion channels','Plasticity via PKA','Plasticity via PKC','Plasticity, others','Plasticity PKA+PKC','Plasticity PKA+PKC no PP','Plasticity, only modelled']

  edgecolors = ['#888888','#9999FF','#FF9999','#FF99FF','#9999FF','#9999FF','#9999FF','#9999FF','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999','#FF9999']
  hatches = [None,None,None,None,'|||','----','////','oooo','||||','----','////','++++','xxxx','oooo']


  corrs = A['corrs'][igroup]
  betas_pvals = A['betas_pvals'][igroup]

  for iPRS in range(0,len(PRSnames_to_consider)):
    iPRSname = -1
    for iPRS0 in range(0,len(PRSnames)):
      if PRSnames[iPRS0] == PRSnames_to_consider[iPRS]:
        iPRSname = iPRS0
        break
    if iPRSname == -1:
      print("Error: PRS not found, "+PRSnames_to_consider[iPRS])
      exit()
    if icov==0:
      axarr[0].bar(iPRS,corrs[iPRSname][0][ithr],facecolor=cols[iPRS],hatch=hatches[iPRS],edgecolor=edgecolors[iPRS],linewidth=0,label=PRSnames[iPRS])
      axarr[1].bar(0.3,0.6,bottom=14-iPRS,width=0.2,facecolor=cols[iPRS],hatch=hatches[iPRS],edgecolor=edgecolors[iPRS],linewidth=0,label=PRSnames[iPRS])
      axarr[1].text(0,14-iPRS+0.3,PRSlabel_names[iPRS], fontsize=5.3,ha='right',va='center')

    printZero = 0 if (corrs[iPRSname][0][ithr] < 0 and abs(corrs[iPRSname][0][ithr]) > 0.07) or corrs[iPRSname][0][ithr] > 0 else 1
    fontweight = 'bold' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else 'normal'
    print('igroup = '+str(igroup)+', '+PRSnames[iPRSname]+', ithr = '+str(ithr)+': pval = '+str(betas_pvals[iPRSname][ithr][1][1]))
    axarr[1].text(1+2.5*icov,14-iPRS+0.3,'p='+'{:.3f}'.format(betas_pvals[iPRSname][ithr][1][1])+('$\\ast\!\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05/14 else ('$\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else '')), fontsize=5.3,ha='left',va='center',fontweight=fontweight)


  axarr[0].set_xlim([-0.6,13.6])
  axarr[0].plot([-1,14],[0,0],'k-',lw=0.3)
  axarr[0].set_xticks([])
  axarr[0].set_ylim([-0.15,0.01])
  axarr[1].text(0.8+2.5*icov,17,axlabels[icov],ha='left',va='top',fontsize=5.4)
  axarr[1].set_xlim([0.2,12.6])

pos = axarr[1].get_position()
for iax in range(0,5):
  f.text(pos.x0 +0.002 + 0.157*iax, pos.y1+0.024, chr(ord('A')+iax), fontsize=11)
axarr[0].set_ylabel('Corr. coeff.',fontsize=6,labelpad=1)
axarr[1].set_yticklabels([])
axarr[2].set_yticklabels([])
for q in ['top','bottom','left','right']:
  axarr[1].spines[q].set_linewidth(0.0)
axarr[1].set_xticks([])
axarr[1].set_yticks([])
f.savefig("figS9A-E.pdf")
      

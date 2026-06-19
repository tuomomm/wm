#cp drawfigcorrs.py drawfig5.py
from pylab import *
import scipy.io
import mytools

thrs = [5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 1e-10]

diaggroup = 'both'
A = scipy.io.loadmat('../genetic/corrs_ztransformed_LNSonly_againstAgeSex_'+diaggroup+'.mat')

PRSnames = A['attrs'][0]
for iname in range(0,len(PRSnames)):
  ifirstspace = PRSnames[iname].find(' ')
  if ifirstspace > -1:
    PRSnames[iname] = PRSnames[iname][0:ifirstspace]

groupNames = ['dsp_r','dsp_r_1','lns_r','MX_lns_r','MX_wms_r','dsp_avg','lns_avg','avg(dsp,lns)','avg(dsp,wms)','avg(lns,wms)','avg(dsp,lns,wms)']

PRSnames_to_consider = ['genes_all_modeled_ionchannels2','genes_all_modeled_ionchannels2_CACNA1C','genes_all_modeled_ionchannels2_CACNA1D','genes_all_modeled_ionchannels2_CACNA1I',
                        'genes_all_modeled_ionchannels2_ATP2A2','genes_all_modeled_ionchannels2_HCN1','genes_all_modeled_ionchannels2_KCNB1','genes_all_modeled_ionchannels2_KCND3','genes_all_modeled_ionchannels2_KCNJ6',
                        'genes_all_modeled_ionchannels2_KCNQ3','genes_all_modeled_ionchannels2_KCNMA1','genes_all_modeled_ionchannels2_SCN1B']

cols = ['#6666FF']*len(PRSnames_to_consider)
PRSlabel_names = ['All modelled\nion channels','CACNA1C','CACNA1D','CACNA1I','ATP2A2','HCN1','KCNB1','KCND3','KCNJ6','KCNQ3','KCNMA1','SCN1B']

edgecolors = ['#9999FF','#FF9999','#FF9999','#FF9999','#888888','#99FF99','#FFFF99','#FFFF99','#FFFF99','#FFFF99','#FFFF99','#FFAFFF']
hatches = ['oooo','++++','////','xxxx','++++','++++','++++','----','||||','////','xxxx','++++']

ithrs = [9,12,5]
igroups = [2,2,2]
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
axarr[0].set_position([0.08,0.56+0.1,0.24,0.3])
axarr[1].set_position([0.08+0.32*1,0.56+0.20,0.20,0.15])
axarr[2].set_position([0.08+0.32*1,0.56-0.06,0.20,0.15])
axarr[3].set_visible(False)

for iithr in range(0,len(igroups)):
  ithr = ithrs[iithr]
  igroup = igroups[iithr]
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
      axarr[iithr].text(iPRS,0.001*(1-2*0),PRSlabel_names[iPRS]+': N/A', fontsize=5.5,rotation=90,ha='center',va='bottom' if not 0 else 'top',fontweight=fontweight,color='#888888')
      continue
    if ithr >= len(corrs[iPRSname][0]):
      print("No data for "+PRSnames_to_consider[iPRS]+" ithr="+str(ithr))
      axarr[iithr].text(iPRS,0.001*(1-2*0),PRSlabel_names[iPRS]+': N/A', fontsize=5.5,rotation=90,ha='center',va='bottom' if not 0 else 'top',fontweight=fontweight,color='#888888')
      continue
    axarr[iithr].bar(iPRS,corrs[iPRSname][0][ithr],facecolor=cols[iPRS],hatch=hatches[iPRS],edgecolor=edgecolors[iPRS],linewidth=0,label=PRSnames[iPRS])
    printBelow = 1 if (corrs[iPRSname][0][ithr] < 0 and abs(corrs[iPRSname][0][ithr]) > 0.07) or corrs[iPRSname][0][ithr] > 0 else 0
    if iithr > 0:
      printBelow = 1 if (corrs[iPRSname][0][ithr] < 0 and abs(corrs[iPRSname][0][ithr]) > 0.15) or corrs[iPRSname][0][ithr] > 0 else 0
      
    fontweight = 'bold' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else 'normal'
    axarr[iithr].text(iPRS,0.001*(1-2*printBelow),'p='+'{:.3f}'.format(betas_pvals[iPRSname][ithr][1][1])+('$\\ast\!\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05/14 else ('$\\ast$' if betas_pvals[iPRSname][ithr][1][1] < 0.05 else '')), fontsize=6,rotation=90,ha='center',va='bottom' if not printBelow else 'top',fontweight=fontweight)

    print('igroup = '+str(igroup)+', '+PRSnames[iPRSname]+', ithr = '+str(ithr)+': pval = '+str(betas_pvals[iPRSname][ithr][1][1]))
  if iithr == 0:
    axarr[iithr].set_ylim([-0.1,0.07])
    axarr[iithr].text(1.0,0.075,'SNP P-val threshold '+str(thrs[ithr]),ha='left',va='top',fontsize=7)
  else:
    axarr[iithr].set_ylim([-0.12,0.07])
    axarr[iithr].text(-0.4-0.6*(iithr==2),0.17,'SNP P-val threshold '+str(thrs[ithr]),ha='left',va='top',fontsize=7)
  axarr[iithr].set_xlim([-0.6,11.6])
  axarr[iithr].plot([-1,14],[0,0],'k-',lw=0.3)
  axarr[iithr].set_xticks([])

labelsToPlot = ['   '+x.replace('\n','\n   ') for x in PRSlabel_names]
colsToPlot = cols[:]
hatchesToPlot = hatches[:]
edgecolorsToPlot = edgecolors[:]
widthsToPlot = [0.0]*len(PRSlabel_names)

myleg = mytools.mylegendbars_distry(f,[0.02,0.5,0.33,0.13],labelsToPlot,[4,4,4],4.2,0.1,0.5,colsToPlot,hatchesToPlot,edgecolorsToPlot,widthsToPlot,myfontsize=5)
for q in ['top','bottom','left','right']:
  myleg.spines[q].set_linewidth(0.0)
myleg.set_xlim([-0.3,12.0])
myleg.set_facecolor('none')

for iax in range(0,3):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.05, pos.y1 - 0.01+0.02*(iax==0), chr(ord('A')+iax), fontsize=11)
  axarr[iax].set_ylabel('Corr. coeff.',fontsize=6)
  
f.savefig("fig6A-C.pdf")
      

from pylab import *
import scipy.io
import mytools


genes = ['ATP2A2','CACNA1C','CACNA1D','CACNA1I','HCN1','KCNB1','KCND3','KCNJ6','KCNQ3','KCNMA1','SCN1B']
cols = ['#6666FF']*len(genes)
edgecolors = ['#888888','#FF9999','#FF9999','#FF9999','#99FF99','#FFFF99','#FFFF99','#FFFF99','#FFFF99','#FFFF99','#FFAFFF']
hatches = ['++++','++++','////','xxxx','++++','++++','----','||||','////','xxxx','++++']


input_file = open('../genetic/SCZ_MetaBrain_MR_Results.txt','r')
firstline = input_file.readline()
headers = firstline.split('\t')
line = input_file.readline()
MRdata = []
while len(line) > 0:
  data = line.split('\t')
  MRdata.append(data[:])
  line = input_file.readline()
input_file.close()

f,axs = subplots(1,2)
axarr = axs.reshape(prod(axs.shape),).tolist()
for iax in [0,1]:
  axarr[iax].set_position([0.1+0.29*iax,0.6,0.2,0.29])
  axarr[iax].plot([-1,11],[0,0],'k-',lw=0.6)
  axs[iax].tick_params(axis='both', which='major', labelsize=4, length=2)
  for axis in ['top','bottom','left','right']:
    axs[iax].spines[axis].set_linewidth(0.0)
  for axis in ['left']:
    axs[iax].spines[axis].set_linewidth(0.3)

nPlotted = 0
genesPlotted = []
for iline in range(0,len(MRdata)):
  #if MRdata[iline][0] != 'Whole_Brain':
  #  continue
  for igene in range(0,len(genes)):
    if MRdata[iline][3] == genes[igene] and MRdata[iline][2] == 'SCZ': #MR_outcome = SCZ, MR_exposure = gene expression
      fontweight = 'bold' if float(MRdata[iline][8]) < 0.05 else 'normal'
      printBelow = float(MRdata[iline][8]) < 0.75
      axarr[0].bar(nPlotted,-log(float(MRdata[iline][8]))/log(10),facecolor=cols[igene],hatch=hatches[igene],edgecolor=edgecolors[igene],linewidth=0)
      axarr[0].text(nPlotted,0.2*(not printBelow)-0.2*printBelow,genes[igene]+('$\\ast\!\\ast$' if float(MRdata[iline][8]) < 0.05/len(genes) else ('$\\ast$' if float(MRdata[iline][8]) < 0.05 else '')),fontsize=6,rotation=90,ha='center',va='bottom' if not printBelow else 'top',fontweight=fontweight)
      nPlotted = nPlotted + 1
      genesPlotted.append(genes[igene])
for igene in range(0,len(genes)):
  if genes[igene] not in genesPlotted:
    axarr[0].text(nPlotted,0.2,genes[igene]+': N/A',fontsize=5.5,rotation=90,ha='center',va='bottom' if not printBelow else 'top',fontweight=fontweight,color='#888888')
    nPlotted = nPlotted + 1
  
nPlotted = 0
genesPlotted = []
for iline in range(0,len(MRdata)):
  #if MRdata[iline][0] != 'Whole_Brain':
  #  continue
  for igene in range(0,len(genes)):
    if MRdata[iline][2] == genes[igene] and MRdata[iline][3] == 'SCZ': #MR_outcome = gene expression, MR_exposure = SCZ
      fontweight = 'bold' if float(MRdata[iline][8]) < 0.05 else 'normal'
      print(str(MRdata[iline]))
      printBelow = float(MRdata[iline][8]) < 0.75
      axarr[1].bar(nPlotted,-log(float(MRdata[iline][8]))/log(10),facecolor=cols[igene],hatch=hatches[igene],edgecolor=edgecolors[igene],linewidth=0)
      axarr[1].text(nPlotted,0.2*(not printBelow)-0.2*printBelow,genes[igene]+('$\\ast\!\\ast$' if float(MRdata[iline][8]) < 0.05/len(genes) else ('$\\ast$' if float(MRdata[iline][8]) < 0.05 else '')),fontsize=6,rotation=90,ha='center',va='bottom' if not printBelow else 'top',fontweight=fontweight)
      nPlotted = nPlotted + 1
      genesPlotted.append(genes[igene])
for igene in range(0,len(genes)):
  if genes[igene] not in genesPlotted:
    axarr[1].text(nPlotted,0.2,genes[igene]+': N/A',fontsize=5.5,rotation=90,ha='center',va='bottom' if not printBelow else 'top',fontweight=fontweight,color='#888888')
    nPlotted = nPlotted + 1

for iax in [0,1]:
  axarr[iax].set_xlim([-1,11])
  axarr[iax].set_xticks([])
  axarr[iax].set_ylabel('-log$_{10}$ $p$',fontsize=7)

axarr[0].set_title('xQTL -> Trait',fontsize=7)
axarr[1].set_title('Trait -> xQTL',fontsize=7)
for iax in range(0,2):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.05, pos.y1 + 0.01, chr(ord('D')+iax), fontsize=11)

f.savefig("fig6D-E.pdf")      

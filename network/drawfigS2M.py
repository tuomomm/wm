from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

#error("Run this on SAGA, the files take more than 5GB!")
def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0,whitecol='#FFFFFF',blackcol='#000000'):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'-',linewidth=lw2,color=blackcol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'-',linewidth=lw2,color=blackcol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'-',color=whitecol,zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

addition = ''

cols = ['#FF3333','#3333FF']
dimcols = ['#FFBBBB','#BBBBFF']
dimgray = '#CCCCCC'
cols_compartments = ['#000000','#888888','#FF4444','#4444FF','#33FF33']

ipop = 3 #Stimulated pyramidal cells

Nsum = 100
Ncellsperpop = [582,97,97,291,291]
FRThrCoeff = 2.0

areas = ['PFC','ACC']
iareas = [1,0] # Plot first ACC, then PFC

channels = ['CTRL', 'nax_PConly', 'iar_PConly', 'kap', 'cat', 'kdr_PConly', 'serca', 'cagk_ikc', 'km', 'cal_PConly', 'GABAB', 'PConly']
channel_titles = ['CTRL','Fast Na$^+$', 'HCN', 'A-type K$^+$', 'T-type Ca$^{2+}$', 'DR K$^+$', 'SERCA', 'BK', 'M-type K$^+$', 'L-type Ca$^{2+}$', 'GABAB', 'Combination']
channel_genes = ['','\n(SCN1B)','\n(HCN1)','\n(KCND3)','\n(CACNA1I)','\n(KCNB1)','\n(ATP2A2)','\n(KCNMA1)','\n(KCNQ3)','\n(CACNA1C,\nCACNA1D)','\n(KCNJ6)','']

f,axs = subplots(4,3)
axarr = axs.reshape(prod(axs.shape),).tolist()
for iay in range(0,3):
  for iax in range(0,4):
    axarr[iay*4+iax].set_position([0.1+0.074*(iax+4*iay),0.75,0.06,0.15])
for iax in range(0,len(axarr)):
  axarr[iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)


ipops = [3]
for iipop in [0]:
  ipop = ipops[iipop]
  
  groupmeans_all = []
  for iiarea in range(0,len(iareas)):
    iarea = iareas[iiarea]
    #Non-subject-wise:
    thisarea_templ = areas[iarea]
    Cavecs_all = []
    CaapicMeans_all = []
    CaapicBaselines_all = []
    for ichannel in range(0,len(channels)):
      print('ichannel = '+str(ichannel))
      thischannel = channels[ichannel]
      thisarea = thisarea_templ+'mild' if thischannel in ['PConly','nax_PConly','kap'] else thisarea_templ
      if 'PFC' in thisarea and thischannel in ['kap','cagk_ikc'] or 'ACC' in thisarea and thischannel in ['iar_PConly','nax_PConly']:
        print('Continue, ichannel='+str(ichannel))
        continue
      bins_all = []
      Cavecs = []
      CaapicMeans = []
      CaapicBaselines = []
      for iseed in range(1,11):
      
        filename = 'FR_sim_net_long_withSK3_'+thisarea+'_'+thischannel+addition+'_altgain_seed'+str(iseed)+'_withCaApic_higherres.mat'
        if thischannel == '':
          filename = 'FR_sim_net_long_withSK3_'+thisarea+addition+'_altgain_seed'+str(iseed)+'_withCaApic_higherres.mat'
        elif thischannel == 'CTRL':
          filename = 'FR_sim_net_long_withSK3_CTRL_altgain_seed'+str(iseed)+'_withCaApic_higherres.mat'

        binwidth = 10
        if not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')

        if not exists(filename):
          print(filename+' does not exist')
          continue
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        Caapics = A['Caapics']
        meanCaTimeCourse = mean([Caapics[i,:] for i in range(0,Caapics.shape[0]) if Caapics[i,0] != 0],axis=0)
        Cavecs.append(meanCaTimeCourse)
        CaapicMeans.append(mean(meanCaTimeCourse[1000:]))
        CaapicBaselines.append(mean(meanCaTimeCourse[7000:10000]))

      Camean = [mean([Cavecs[isamp][it] for isamp in range(0,len(Cavecs))]) for it in range(0,len(Caapics[0]))]

      Cavecs_all.append(Cavecs[:])
      CaapicMeans_all.append(CaapicMeans[:])
      CaapicBaselines_all.append(CaapicBaselines[:])
      print(' appending Cavecs_all.append')

      ichannelToPlot = len(Cavecs_all)-1

      Cabaseline = CaapicBaselines_all[ichannelToPlot]
      if ichannel == 0:
        Cabaseline_CTRL = CaapicBaselines_all[0]

      axarr[ichannel].bar(0,mean(Cabaseline_CTRL)*1e6,facecolor='#888888')
      axarr[ichannel].bar(1+iiarea,mean(Cabaseline)*1e6,facecolor=cols[iarea])
      axarr[ichannel].set_title(channel_titles[ichannel],fontsize=6)
      if 'Comb' in channel_titles[ichannel]:
        axarr[ichannel].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5)
      else:
        axarr[ichannel].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5,pad=2.5)
      axarr[ichannel].set_xlim([-0.8,2.8])
      axarr[ichannel].set_ylim([229,249])
      drawdiscontinuity(axarr[ichannel],232,1,-0.8,0.2,1.0,0.8,'#FFFFFF')
      drawdiscontinuity(axarr[ichannel],232,1.4,0,0.28,1.0,0.8,'#FFFFFF','#888888')
      drawdiscontinuity(axarr[ichannel],232,1.4,1+iiarea,0.28,1.0,0.8,'#FFFFFF',cols[iarea])
      axarr[ichannel].set_yticks([235,240,245])
      axarr[ichannel].set_yticklabels([])
      axarr[ichannel].set_xticks([])
for iax in range(0,1):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.04, pos.y1 - 0.0, chr(ord('A')+iax), fontsize=9)

axarr[0].set_yticklabels(['235','240','245'])
  
for iax in [0]:
  axarr[iax].set_ylabel('Apical\n[Ca$^{2+}$] (nM)',fontsize=6)

  
f.savefig("figS2M.pdf")

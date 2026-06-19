from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

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
channel_genes = ['','\n(SCN1B)','\n(HCN1)','\n(KCND3)','\n(CACNA1I)','\n(KCNB1)','\n(ATP2A2)','\n(KCNMA1)','\n(KCNQ3)','\n(CACNA1C,CACNA1D)','\n(KCNJ6)','']
f,axs = subplots(4,3)
axarr = axs.reshape(prod(axs.shape),).tolist()
axnew = []
for iay in range(0,3):
  for iax in range(0,4):
    axarr[iay*4+iax].set_position([0.1+0.23*iax,0.78-0.235*iay,0.18,0.15])
    axnew.append(f.add_axes([0.1+0.23*iax+0.15,0.78-0.235*iay+0.07,0.025,0.07]))
for iax in range(0,len(axarr)):
  axarr[iax].tick_params(axis='both', which='major', labelsize=4)
  axnew[iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
    axnew[iax].spines[axis].set_linewidth(0.3)
axnew[0].set_visible(False)


ipops = [3]
for iipop in [0]:
  ipop = ipops[iipop]
  
  groupmeans_all = []
  for iiarea in range(0,len(iareas)):
    iarea = iareas[iiarea]
    #Non-subject-wise:
    thisarea_templ = areas[iarea]
    bins_all_all = []
    FRvecs_all = []
    for ichannel in range(0,len(channels)):
      thischannel = channels[ichannel]
      thisarea = thisarea_templ+'mild' if thischannel in ['PConly','nax_PConly','kap'] else thisarea_templ
      if 'PFC' in thisarea and thischannel in ['kap','cagk_ikc'] or 'ACC' in thisarea and thischannel in ['iar_PConly','nax_PConly']:
        continue
      bins_all = []
      FRvecs = []
      for iseed in range(1,7):
      
        filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+'_'+thischannel+addition+'_altgain_seed'+str(iseed)+'.mat'
        if thischannel == '':
          filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+addition+'_altgain_seed'+str(iseed)+'.mat'
        elif thischannel == 'CTRL':
          filename = 'FRandMUA_sim_net_long_withSK3_CTRL'+addition+'_altgain_seed'+str(iseed)+'.mat'

        binwidth = 1000
        if not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')

        if not exists(filename):
          print(filename+' does not exist')
          continue
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        FRs = zeros(int(len(A['FRvecs'][ipop])/Nsum)+1,)
        for ibin in range(0,len(FRs)):
          FRs[ibin] = sum(A['FRvecs'][ipop][Nsum*ibin:Nsum*(ibin+1)])/(Ncellsperpop[ipop]*Nsum/100) #Divided by the number of cells in the population, corrected by the bin size (if bin size is 1 s, i.e. Nsum=100, then we have the frequency in Hz without normalization, but if 10 s, i.e. Nsum=1000, then we have to also divide by 10 to get spikes/s)
          baselineFR = mean(FRs[:])
        FRvecs.append(FRs[:])
      FRvecs_all.append(FRvecs[:])

      FRmean = [mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))]
      if thischannel == 'CTRL':
        FRbaseline = mean(FRmean[3:10])
      #Determine the time above 2.0*baseline_HC:                                                                                                                        
      FRThr = FRThrCoeff*FRbaseline
      istarts = [i for i in range(0,len(FRmean)-1) if FRmean[i] < FRThr and FRmean[i+1] >= FRThr and i >= 10]
      iends = [i for i in range(0,len(FRmean)-1) if FRmean[i] >= FRThr and FRmean[i+1] < FRThr and i >= 10]
      if len(istarts) > 0 and len(iends) > 0:
        startt = 10*Nsum*(istarts[0]+(FRThr-FRmean[istarts[0]])/(FRmean[istarts[0]+1]-FRmean[istarts[0]]))
        endt = 10*Nsum*(iends[-1]+(FRmean[iends[-1]]-FRThr)/(FRmean[iends[-1]]-FRmean[iends[-1]+1]))
        FRbeyondThr = (endt-startt)/1000 #change from milliseconds to second                                                                         
      elif len(istarts) > 0:
        FRbeyondThr = inf
      else:
        FRbeyondThr = 0

      if ichannel == 0:
        FRbeyondThr_CTRL = FRbeyondThr
      else:  
        axnew[ichannel].bar(0,FRbeyondThr_CTRL,facecolor='#888888')
        axnew[ichannel].bar(1+iiarea,FRbeyondThr,facecolor=cols[iarea])
        axnew[ichannel].set_ylabel('time above\n'+'{:.2f}'.format(FRThr)+' Hz (s)',fontsize=5)
        print('time above\n'+'{:.2f}'.format(FRThr)+' Hz (s)')
        
      axarr[ichannel].plot(10*Nsum*array(range(0,len(FRs)))/1000-10, [mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))],label='SCZ, '+areas[iarea],color=cols[iarea],lw=0.4)

      FRvecs = FRvecs_all[0]
      if iarea == 1:
        axarr[ichannel].plot(10*Nsum*array(range(0,len(FRs)))/1000-10, [mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))],'k-',lw=0.4)
      else:
        axarr[ichannel].plot(10*Nsum*array(range(0,len(FRs)))/1000-10, [mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))],'k-',label='CTRL',lw=0.4)
      if 'Comb' in channel_titles[ichannel]:
        axarr[ichannel].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5)
      else:
        axarr[ichannel].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5,pad=2.5)

for iax in range(0,len(axarr)):
  axarr[iax].set_xlim([-9,89])
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.04, pos.y1 - 0.0, chr(ord('A')+iax), fontsize=9)
  axnew[iax].set_xlim([-0.5,2.5])
  axnew[iax].set_xticks([])

axarr[0].legend(fontsize=5)
for iax in [0,4,8]:
  axarr[iax].set_ylabel('FR (spikes/s)',fontsize=6)
for iax in [8,9,10,11]:
  axarr[iax].set_xlabel('time (s)',fontsize=6)
  
f.savefig("figS2A-L.pdf")

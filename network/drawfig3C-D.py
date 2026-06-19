#cp drawfig3_supp0b_FRvec.py drawfig3_supp0b.py #13.1.2026
from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

addition = ''

def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

f,axarr = subplots(3,3)
for iax in range(0,3):
  for iay in range(0,2):
    axarr[iay,iax].set_position([0.06+0.158*iax,0.737-0.31*iay,0.11,0.23])
    boxoff(axarr[iay,iax])
for iax in range(0,3):
  boxoff(axarr[2,iax])
axarr[2,0].set_position([0.16,0.22,0.326,0.11])
axarr[2,1].set_position([0.16,0.12,0.326,0.08])
axarr[2,2].set_position([0.16,0.01,0.326,0.08])

def mybar(ax,x,y,facecolor=[],linewidth=0.3,w=0.4):
  qs = quantile(y, [0,0.25,0.5,0.75,1])
  polygon = Polygon(array([[x-w,x+w,x+w,x-w],[qs[1],qs[1],qs[3],qs[3]]]).T)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  if type(facecolor) is not list or len(facecolor) > 0:
    p.set_facecolor(facecolor)
  p.set_edgecolor('#000000')
  p.set_linewidth(0.3)
  ax.add_collection(p)
  if isinf(qs[4]) or isnan(qs[4]):
    a2 = []
    a2.append(ax.plot([x-w,x+w,x,x,x-w,x+w,x,x],[qs[0],qs[0],qs[0],qs[2],qs[2],qs[2],qs[2],qs[3]+(qs[3]-qs[1])*0.5],'k-',lw=linewidth))
    a2.append(ax.plot([x,x],[qs[3]+(qs[3]-qs[1])*0.5,qs[3]+(qs[3]-qs[1])*1.0],'k--',lw=linewidth,dashes=(1.2,1.2)))
  else:
    a2 = ax.plot([x-w,x+w,x,x,x-w,x+w,x,x,x-w,x+w],[qs[0],qs[0],qs[0],qs[2],qs[2],qs[2],qs[2],qs[4],qs[4],qs[4]],'k-',lw=linewidth)
    print(str(qs[4]))
  return [p,a2]

cols = ['#FF3333','#3333FF']
dimcols = ['#FFBBBB','#BBBBFF']
dimgray = '#CCCCCC'
cols_compartments = ['#000000','#888888','#FF4444','#4444FF','#33FF33']

bands = [[0.5,4.0],[4.0,8.0],[8.0,13.0],[13.0,30],[30,70],[70,150]]
bandNames = ['Delta','Theta','Alpha','Beta','Low\ngamma','High\ngamma']

Nsum = 100
Ncellsperpop = [582,97,97,291,291]

areas = ['PFC','ACC']
iareas = [0,1] # Plot first ACC, then PFC

channels = ['CTRL', 'PConly']
channel_titles = ['CTRL','Combination']

binWidth = 100
binWidthF = 20

ipops = [3,1,2]
if True:
  
  groupmeans_all = []
  
  for iiarea in range(0,len(iareas)):
    iarea = iareas[iiarea]
    #Non-subject-wise:
    thisarea = areas[iarea]+'mild'
    bins_all_all = []
    Cavecs_all = []
    FRs_all_all = []
    FRs_PSDs_all_all = []
    bandPSDs_all_all = []
    for ichannel in range(0,len(channels)):
      thischannel = channels[ichannel]
      if thisarea == 'PFCmild' and thischannel in ['kap','cagk_ikc'] or thisarea == 'ACCmild' and thischannel in ['iar_PConly','nax_PConly']:
        continue
      bins_all = []
      Cavecs = []
      FRs_all = []
      FRs_PSDs_all = []
      bandPSDs_all = []
      for iseed in range(1,11):
      
        filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+'_'+thischannel+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'
        if thischannel == '':
          filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'
        elif thischannel == 'CTRL':
          filename = 'FRandMUA_sim_net_long_withSK3_CTRL'+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'

        if not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')

        if not exists(filename):
          print(filename+' or '+filename.replace('FRandMUA','FR')+' does not exist')
          continue
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        FRs = A['FRvecs']
        FRs_PSDs = [abs(fft(FRs[ipop,:]/Ncellsperpop[ipop]*100))**2 for ipop in range(0,5)]
          
        FRs_all.append(FRs[:])
        FRs_PSDs_all.append(FRs_PSDs[:])
            
      dt = 5.0 #the  dt of FRvec is 5.0 ms
      dt_binned = dt*binWidth
      df = 1/len(FRs[0])/dt*1000 #the df of FFT signal in Hz
      df_binned = df*binWidthF
      bandPSDs_all = []
      for iipop in range(0,len(ipops)):
        ipop = ipops[iipop]
        FRs_this = mean(array(FRs_all),axis=0)[ipop]
        if binWidth == 1:
          FRs_binned = FRs_this/Ncellsperpop[ipop]*100
        else:
          FRs_binned = mean(reshape(FRs_this[0:binWidth*int(len(FRs_this)/binWidth)],[int(len(FRs_this)/binWidth),binWidth]),axis=1)/Ncellsperpop[ipop]*100 #Multiply by 100 to get Hz because dt is 0.01 s
          
        if iarea == 0 or ichannel == 1:
          axarr[0,iipop].plot(dt_binned*array(range(0,len(FRs_binned)))/1000-10,FRs_binned,label='SCZ, '+areas[iarea] if ichannel == 1 else 'CTRL',color=cols[iarea] if ichannel == 1 else '#000000',lw=0.4)

        FRs_PSD_this = mean(array(FRs_PSDs_all),axis=0)[ipop]
        FRs_PSD_this_nonmean = array(FRs_PSDs_all)[:,ipop,:]
        if binWidthF == 1:
          FRs_binned_PSDs_binned = FRs_PSD_this
        else:
          FRs_binned_PSDs_binned = mean(reshape(FRs_PSD_this[0:binWidthF*int(len(FRs_PSD_this)/binWidthF)],[int(len(FRs_PSD_this)/binWidthF),binWidthF]),axis=1)

        if iarea == 0 or ichannel == 1:
          axarr[1,iipop].semilogy(df_binned*array(range(0,len(FRs_binned_PSDs_binned))),FRs_binned_PSDs_binned,label='SCZ, '+areas[iarea] if ichannel == 1 else 'CTRL',color=cols[iarea] if ichannel == 1 else '#000000',lw=0.4)
        axarr[1,iipop].set_xlim([0,100])

        bandPSDs = []
        for iband in range(0,len(bands)):
          ifs = [i for i in range(0,FRs_PSD_this.shape[0]) if bands[iband][0] <= df*i < bands[iband][1]]
          ydata = mean([FRs_PSD_this_nonmean[:,i] for i in ifs],axis=0)
          if iband == 0 and (iarea == 0 or ichannel == 1):
            axarr[2,iipop].bar(5*iband+iiarea if ichannel==1 else 5*iband+2,mean(ydata),facecolor=cols[iarea] if ichannel == 1 else '#000000',label='SCZ, '+areas[iarea] if ichannel == 1 else 'CTRL')
            axarr[2,iipop].plot([5*iband+iiarea if ichannel==1 else 5*iband+2]*2,[mean(ydata)-std(ydata),mean(ydata)+std(ydata)],'k-',lw=0.3)
          else:
            axarr[2,iipop].bar(5*iband+iiarea if ichannel==1 else 5*iband+2,mean(ydata),facecolor=cols[iarea] if ichannel == 1 else '#000000')
            axarr[2,iipop].plot([5*iband+iiarea if ichannel==1 else 5*iband+2]*2,[mean(ydata)-std(ydata),mean(ydata)+std(ydata)],'k-',lw=0.3)
          bandPSDs.append(mean(array([FRs_PSD_this_nonmean[:,i] for i in ifs]),axis=0))
        bandPSDs_all.append(bandPSDs[:])
      FRs_all_all.append(FRs_all[:])
      FRs_PSDs_all_all.append(FRs_PSDs_all[:])
      bandPSDs_all_all.append(bandPSDs_all[:])
        
    for iipop in range(0,3):
      for iband in range(0,6):
        pval = scipy.stats.ranksums(bandPSDs_all_all[0][iipop][iband],bandPSDs_all_all[1][iipop][iband])[1]
        print('FR PSD '+thisarea+' iipop='+str(iipop)+' '+bandNames[iband].replace('\n','')+' pval='+str(pval))
        if pval < 0.05/6:
          if iipop==0:
            ylevel = 1.3e5 if iband != 3 else 4e5
            ydelta = 0.1e5
          elif iipop==1:
            ylevel = 4e5 if iband in [0,4,5] else (6.5e5 if iband != 3 else 1e6)
            ydelta = 0.4e5
          elif iipop==2:
            ylevel = 1.0e6 if iband != 3 else 4.9e6
            ydelta = 0.2e6
          axarr[2,iipop].plot([5*iband+1-iiarea,5*iband+1-iiarea,5*iband+2,5*iband+2],[ylevel+ydelta*5*iiarea,ylevel+ydelta*5*iiarea+ydelta,ylevel+ydelta*5*iiarea+ydelta,ylevel+ydelta*5*iiarea],'k-',lw=0.3)
          axarr[2,iipop].text(mean([5*iband+(1-iiarea),5*iband+2]),ylevel+ydelta*5*iiarea+ydelta,'*',fontsize=5,ha='center',va='center')

for iax in range(0,3):
  for iay in [0,1,2]:
    axarr[iay,iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[0,iax].spines[axis].set_linewidth(0.3)
    axarr[1,iax].spines[axis].set_linewidth(0.3)
    axarr[2,iax].spines[axis].set_linewidth(0.3)
  axarr[0,iax].set_xlim([-9,90])
  axarr[0,iax].set_xlabel('$t$ (s)',fontsize=6)
  axarr[1,iax].set_xlabel('f (Hz)',fontsize=6)
  axarr[1,iax].set_ylim([1e4,1e8])
axarr[0,0].set_ylabel('Firing rate (Hz)',fontsize=6)
axarr[1,0].set_ylabel('PSD (A.U.)',fontsize=6)

axarr[1,0].set_title('Stimulated PCs',fontsize=6)
axarr[1,1].set_title('LTS cells',fontsize=6)
axarr[1,2].set_title('PV+ cells',fontsize=6)

axarr[0,0].legend(fontsize=4)
axarr[1,0].legend(fontsize=4)
axarr[2,2].legend(fontsize=4,ncol=2)

for iband in range(0,len(bands)):
  axarr[2,0].text(5*iband+1,6.4e5,bandNames[iband],fontsize=6,ha='center',va='top')

axarr[2,0].set_ylim([0,6e5])
axarr[2,0].set_yticks([0,2e5,4e5])
axarr[2,0].set_yticklabels(['0','2$\\times10^5$','4$\\times10^5$'])

axarr[2,1].set_ylim([0,13e5])
axarr[2,1].set_yticks([0,5e5,10e5])
axarr[2,1].set_yticklabels(['0','5$\\times10^5$','1$\\times10^6$'])

axarr[2,2].set_ylim([0,6.5e6])
axarr[2,2].set_yticks([0,3e6,6e6])
axarr[2,2].set_yticklabels(['0','3$\\times10^6$','6$\\times10^6$'])

for iax in range(0,3):
  axarr[2,iax].set_xticks([])
  axarr[2,iax].set_xlim([-0.5,27.5])

axarr[2,0].set_ylabel('Band-wise\nmean PSD\n(Stimulated\nPCs)',fontsize=6,rotation=0,loc='center',va='center',labelpad=20)
axarr[2,1].set_ylabel('Band-wise\nmean PSD\n(LTS cells)',fontsize=6,rotation=0,loc='center',va='center',labelpad=20)
axarr[2,2].set_ylabel('Band-wise\nmean PSD\n(PV+ cells)',fontsize=6,rotation=0,loc='center',va='center',labelpad=20)

toplotchr = [axarr[0,0],axarr[1,0],axarr[2,0]]
letters = []
for ichr in range(1,len(toplotchr)):
  pos = toplotchr[ichr].get_position()
  letters.append(f.text(pos.x0 - 0.055 - 0.1*(ichr==2), pos.y1 - 0.03*(ichr==2), chr(ord('B')+ichr), fontsize=12))

axarr[0,0].set_visible(False)
axarr[0,1].set_visible(False)
axarr[0,2].set_visible(False)

f.savefig("fig3C-D.pdf")

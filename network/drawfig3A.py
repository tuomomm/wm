#cp drawfig3_Ca_withSK3_altgain.py drawfig3_supp0.py
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

f,axarr = subplots(1,2)
boxoff(axarr[0])

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


axarr[0].set_position([0.06,0.66,0.426,0.15])
axnew = f.add_axes([0.4,0.77,0.07,0.05])
boxoff(axnew)

    
cols = ['#FF3333','#3333FF']
dimcols = ['#FFBBBB','#BBBBFF']
dimgray = '#CCCCCC'
cols_compartments = ['#000000','#888888','#FF4444','#4444FF','#33FF33']

ipop = 3 #Stimulated pyramidal cells

Nsum = 100
Ncellsperpop = [582,97,97,291,291]

areas = ['PFC','ACC']
iareas = [1,0] # Plot first ACC, then PFC

channels = ['CTRL', 'PConly']
channel_titles = ['CTRL','Combination']

FRThrCoeff = 2.0

if True:  
  groupmeans_all = []
  for iiarea in [1,0]:
    iarea = iareas[iiarea]
    #Non-subject-wise:
    thisarea = areas[iarea]+'mild'
    bins_all_all = []
    Cavecs_all = []
    CasomaMeans_all = []
    for ichannel in range(0,len(channels)):
      thischannel = channels[ichannel]
      if thisarea == 'PFCmild' and thischannel in ['kap','cagk_ikc'] or thisarea == 'ACCmild' and thischannel in ['iar_PConly','nax_PConly']:
        continue
      bins_all = []
      Cavecs = []
      CasomaMeans = []
      for iseed in range(1,11):
      
        filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+'_'+thischannel+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'
        if thischannel == '':
          filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'
        elif thischannel == 'CTRL':
          filename = 'FRandMUA_sim_net_long_withSK3_CTRL'+addition+'_altgain_seed'+str(iseed)+'_withCa_higherres.mat'

        binwidth = 1000
        if not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')

        if not exists(filename):
          print(filename+' does not exist')
          continue
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        Casomas = A['Casomas']
        meanCaTimeCourse = mean([Casomas[i,:] for i in range(0,Casomas.shape[0]) if Casomas[i,0] != 0],axis=0)
        Cavecs.append(meanCaTimeCourse)
        CasomaMeans.append(mean(meanCaTimeCourse[1000:]))

      Camean = [mean([Cavecs[isamp][it] for isamp in range(0,len(Cavecs))]) for it in range(0,len(Casomas[0]))]

      Cavecs_all.append(Cavecs[:])
      CasomaMeans_all.append(CasomaMeans[:])

    Camean_CTRL = [mean([Cavecs_all[0][isamp][it] for isamp in range(0,len(Cavecs_all[0]))]) for it in range(0,len(Casomas[0]))]
    Camean = [mean([Cavecs_all[1][isamp][it] for isamp in range(0,len(Cavecs_all[1]))]) for it in range(0,len(Casomas[0]))]
    binWidth = 10
    dt_binned = 0.5*binWidth
    Camean_CTRL_binned = mean(reshape(Camean_CTRL[0:binWidth*int(len(Camean_CTRL)/binWidth)],[int(len(Camean_CTRL)/binWidth),binWidth]),axis=1)
    Camean_binned = mean(reshape(Camean[0:binWidth*int(len(Camean)/binWidth)],[int(len(Camean)/binWidth),binWidth]),axis=1)
    if iarea == 1:
      axarr[0].plot(10*binWidth*array(range(0,len(Camean_CTRL_binned)))/1000-10, 1e6*Camean_CTRL_binned,'k-',label='CTRL',lw=0.4)
    axarr[0].plot(10*binWidth*array(range(0,len(Camean_binned)))/1000-10, 1e6*Camean_binned,label='SCZ, '+areas[iarea],color=cols[iarea],lw=0.4)
    pval = scipy.stats.ranksums(CasomaMeans_all[0],CasomaMeans_all[1])[1]
    print(thisarea+' pval='+str(pval))
    axnew.bar(iiarea,mean(CasomaMeans_all[1])*1e6,facecolor=cols[iarea],label='SCZ, '+areas[iarea])
    axnew.plot([iiarea,iiarea],[mean(CasomaMeans_all[1])*1e6-std(CasomaMeans_all[1])*1e6,mean(CasomaMeans_all[1])*1e6+std(CasomaMeans_all[1])*1e6],'k-',lw=0.3)
    if iarea == 1:
      axnew.bar(2,mean(CasomaMeans_all[0])*1e6,facecolor='#000000',label='CTRL')
      axnew.plot([2,2],[mean(CasomaMeans_all[0])*1e6-std(CasomaMeans_all[0])*1e6,mean(CasomaMeans_all[0])*1e6+std(CasomaMeans_all[0])*1e6],'k-',lw=0.3)

    if pval < 0.05:
      axnew.plot([iiarea,iiarea,2,2],[x+3.5*(iiarea==0) for x in [13,14,14,13]],'k-',lw=0.2)
      axnew.text(mean([iarea,2]),14+3.5*iiarea,'*',fontsize=5,ha='center',va='center')

    
axarr[0].set_ylim([4,18])
axarr[0].set_xlim([-9,90])
handles, labels = axarr[0].get_legend_handles_labels()
l=axarr[0].legend([handles[idx] for idx in [1,0,2]],[labels[idx] for idx in [1,0,2]],fontsize=4,ncol=3,loc=8)

axarr[1].set_visible(False)

axarr[0].tick_params(axis='both', which='major', labelsize=4)
axnew.tick_params(axis='both', which='major', labelsize=4)
for axis in ['top','bottom','left','right']:
  axarr[0].spines[axis].set_linewidth(0.3)
  axnew.spines[axis].set_linewidth(0.3)

axarr[0].set_ylabel('[Ca$^{2+}$] (nM)',fontsize=6)
axarr[0].set_xlabel('$t$ (s)',fontsize=6)
axnew.set_ylabel('[Ca$^{2+}$]\n(nM)',fontsize=5)

pos = axarr[0].get_position()
f.text(pos.x0 - 0.06, pos.y1 - 0.01, chr(ord('A')+0), fontsize=12)


f.savefig("fig3A.pdf")

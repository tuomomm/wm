#cp drawfig3_Ca_withSK3_altgain.py drawfig3_supp0.py
from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

#error("Run this on SAGA, the files take more than 5GB!")

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

f,axarr = subplots(1,3)
boxoff(axarr[0])
boxoff(axarr[1])
boxoff(axarr[2])

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

axarr[0].set_position([0.076,0.62,0.11,0.26])
axarr[1].set_position([0.226,0.62,0.11,0.26])
axarr[2].set_position([0.376,0.62,0.11,0.26])
axnew = []
for iax in range(0,3):
  axnew.append(f.add_axes([0.076+0.075+0.15*iax - 0.02*(iax>0),0.64+0.14*(iax>0),0.02,0.09]))
  boxoff(axnew[iax],'bottom')
  axnew[iax].set_xticks([])
    
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
  print('Loading data/sim_net_long_withSK3_PFCmild_altgain_pc_*.npz to get the lctyID vector')
  ld = {}
  for pcid in range(0,8):
    ld[pcid]=load('data/sim_net_long_withSK3_PFCmild_altgain_pc_'+str(pcid)+'.npz')
  lctyID = ld[0]['lctyID']

  ipops = [3,1,2]
  inds = [[],[],[],[],[]]
  inds8 = [i for i in range(0,len(lctyID)) if lctyID[i] == 8]
  inds9 = [i for i in range(0,len(lctyID)) if lctyID[i] == 9]
  inds11 = [i for i in range(0,len(lctyID)) if lctyID[i] == 11]
  inds12 = [i for i in range(0,len(lctyID)) if lctyID[i] == 12]
  inds13 = [i for i in range(0,len(lctyID)) if lctyID[i] == 13]
  inds14 = [i for i in range(0,len(lctyID)) if lctyID[i] == 14]
  inds15 = [i for i in range(0,len(lctyID)) if lctyID[i] == 15]
  inds19 =[i for i in range(0,len(lctyID)) if lctyID[i] == 19]
  inds21 =[i for i in range(0,len(lctyID)) if lctyID[i] == 21]
  inds24 =[i for i in range(0,len(lctyID)) if lctyID[i] == 24]
  for ilctyID in range(0,len(lctyID)):
    if ilctyID in inds8+inds12+inds13+inds19: #'r', i.e., FRvecs[0]
      inds[0].append(ilctyID)
    if ilctyID in inds9+inds14+inds21: #'g', i.e., FRvecs[2]
      inds[1].append(ilctyID)
    if ilctyID in inds11+inds15+inds24: #'b', i.e., FRvecs[1]
      inds[2].append(ilctyID)
    if inds8[0] <= ilctyID < inds8[int(len(inds8)/2)] or inds12[0] <= ilctyID < inds12[int(len(inds12)/2)] or inds13[0] <= ilctyID < inds13[int(len(inds13)/2)] or inds19[0] <= ilctyID < inds19[int(len(inds19)/2)]: #half of 'r', i.e., FRvecs[3]
      inds[3].append(ilctyID)
    if inds8[int(len(inds8)/2)] <= ilctyID < inds8[-1] or inds12[int(len(inds12)/2)] <= ilctyID < inds12[-1] or inds13[int(len(inds13)/2)] <= ilctyID < inds13[-1] or inds19[int(len(inds19)/2)] <= ilctyID < inds19[-1]: #half of 'r', i.e., FRvecs[4]
      inds[4].append(ilctyID)
    
  groupmeans_all = []
  for iiarea in [1,0]:
    iarea = iareas[iiarea]
    #Non-subject-wise:
    thisarea = areas[iarea]+'mild'
    bins_all_all = []
    Cavecs_all = []
    NMDAMeans_all = []
    for ichannel in range(0,len(channels)):
      thischannel = channels[ichannel]
      if thisarea == 'PFCmild' and thischannel in ['kap','cagk_ikc'] or thisarea == 'ACCmild' and thischannel in ['iar_PConly','nax_PConly']:
        continue
      bins_all = []
      Cavecs = []
      NMDAMeans = []
      for iseed in range(1,11):
      
        filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+'_'+thischannel+addition+'_altgain_withallNMDA_seed'+str(iseed)+'.mat' #'_withCa_higherres.mat'
        if thischannel == '':
          filename = 'FRandMUA_sim_net_long_withSK3_'+thisarea+addition+'_withallNMDA_altgain_seed'+str(iseed)+'.mat' #'_withCa_higherres.mat'
        elif thischannel == 'CTRL':
          filename = 'FRandMUA_sim_net_long_withSK3_CTRL'+addition+'_altgain_withallNMDA_seed'+str(iseed)+'.mat' #'_withCa_higherres.mat'

        binwidth = 1000
        if not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')

        if not exists(filename):
          print(filename+' does not exist')
          continue
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        NMDAs = A['NMDAcurrs'][0]+A['NMDAcurrs'][1]+A['NMDAcurrs'][2]+A['NMDAcurrs'][3]+A['NMDAcurrs'][4] #Sum up the NMDA currents. They are all in nA.
        Cavecs_thisseed = []
        NMDAMeans_thisseed = []
        for iipop in range(0,3):
          ipop = ipops[iipop]
          meanCaTimeCourse = mean([NMDAs[i,:] for i in range(0,NMDAs.shape[0]) if i in inds[ipop]],axis=0)
          Cavecs_thisseed.append(meanCaTimeCourse) 
          NMDAMeans_thisseed.append(mean(meanCaTimeCourse[1000:]))
        Cavecs.append(Cavecs_thisseed[:])
        NMDAMeans.append(NMDAMeans_thisseed[:])
      Camean = [[mean([Cavecs[isamp][iipop][it] for isamp in range(0,len(Cavecs))]) for it in range(0,len(NMDAs[0]))] for iipop in range(0,len(ipops))]

      Cavecs_all.append(Cavecs[:])
      NMDAMeans_all.append(NMDAMeans[:])

    for iipop in range(0,len(ipops)):
      Camean_CTRL = [mean([Cavecs_all[0][isamp][iipop][it] for isamp in range(0,len(Cavecs_all[0]))]) for it in range(0,len(NMDAs[0]))]
      Camean = [mean([Cavecs_all[1][isamp][iipop][it] for isamp in range(0,len(Cavecs_all[1]))]) for it in range(0,len(NMDAs[0]))]
      binWidth = 10
      dt_binned = 0.5*binWidth
      Camean_CTRL_binned = mean(reshape(Camean_CTRL[0:binWidth*int(len(Camean_CTRL)/binWidth)],[int(len(Camean_CTRL)/binWidth),binWidth]),axis=1)
      Camean_binned = mean(reshape(Camean[0:binWidth*int(len(Camean)/binWidth)],[int(len(Camean)/binWidth),binWidth]),axis=1)
      if iarea == 1:
        axarr[iipop].plot(10*binWidth*array(range(0,len(Camean_CTRL_binned)))/1000-10, 1e3*Camean_CTRL_binned,'k-',label='CTRL',lw=0.4)
      axarr[iipop].plot(10*binWidth*array(range(0,len(Camean_binned)))/1000-10, 1e3*Camean_binned,label='SCZ, '+areas[iarea],color=cols[iarea],lw=0.4)
      NMDAMeans0 = [NMDAMeans_all[0][j][iipop] for j in range(0,len(NMDAMeans_all[0]))]
      NMDAMeans1 = [NMDAMeans_all[1][j][iipop] for j in range(0,len(NMDAMeans_all[1]))]
      pval = scipy.stats.ranksums(NMDAMeans0,NMDAMeans1)[1]
      print(thisarea+' pval='+str(pval))
      print(str(NMDAMeans0)+' '+str(NMDAMeans1))
      axnew[iipop].bar(iiarea,mean(NMDAMeans1)*1e3,facecolor=cols[iarea],label='SCZ, '+areas[iarea])
      axnew[iipop].plot([iiarea,iiarea],[mean(NMDAMeans1)*1e3-std(NMDAMeans1)*1e3,mean(NMDAMeans1)*1e3+std(NMDAMeans1)*1e3],'k-',lw=0.3)
      if iarea == 1:
        axnew[iipop].bar(2,mean(NMDAMeans0)*1e3,facecolor='#000000',label='CTRL')
        axnew[iipop].plot([2,2],[mean(NMDAMeans0)*1e3-std(NMDAMeans0)*1e3,mean(NMDAMeans0)*1e3+std(NMDAMeans0)*1e3],'k-',lw=0.3)

      if pval < 0.05:
        if iipop == 0:
          axnew[iipop].plot([iiarea,iiarea,2,2],[x-0.350*(iiarea==0) for x in [-1.600,-1.700,-1.700,-1.600]],'k-',lw=0.2)
          axnew[iipop].text(mean([iarea,2]),-1.960-0.350*iiarea,'*',fontsize=5,ha='center',va='center')
        elif iipop == 1:
          axnew[iipop].plot([iiarea,iiarea,2,2],[x-8*(iiarea==0) for x in [-32,-34,-34,-32]],'k-',lw=0.2)
          axnew[iipop].text(mean([iarea,2]),-39.2-8*iiarea,'*',fontsize=5,ha='center',va='center')
        elif iipop == 2:
          axnew[iipop].plot([iiarea,iiarea,2,2],[x-2*(iiarea==0) for x in [-8,-8.5,-8.5,-8]],'k-',lw=0.2)
          axnew[iipop].text(mean([iarea,2]),-9.8-2*iiarea,'*',fontsize=5,ha='center',va='center')

for iipop in range(0,3):
  axarr[iipop].set_xlim([-9,90])

  axarr[iipop].tick_params(axis='both', which='major', labelsize=4)
  axnew[iipop].tick_params(axis='both', which='major', labelsize=4, pad=0.5, length=1)
  for axis in ['top','bottom','left','right']:
    axarr[iipop].spines[axis].set_linewidth(0.3)
    axnew[iipop].spines[axis].set_linewidth(0.3)

  axarr[iipop].set_xlabel('$t$ (sec)',fontsize=6)

  pos = axarr[iipop].get_position()
  if iipop == 0:
    f.text(pos.x0 - 0.06, pos.y1 - 0.03, chr(ord('B')+iipop), fontsize=12)

  axnew[iipop].xaxis.tick_top()
  axnew[iipop].xaxis.set_label_position('top')
  axnew[iipop].set_ylabel('Avg.\n$I_{\\mathrm{NMDA}}$ (nA)',fontsize=5,labelpad=1.0)

axarr[0].set_title('Stimulated PCs\n',fontsize=6)
axarr[1].set_title('LTS cells\n',fontsize=6)
axarr[2].set_title('PV+ cells\n',fontsize=6)

handles, labels = axarr[2].get_legend_handles_labels()
l=axarr[0].legend([handles[idx] for idx in [1,0,2]],[labels[idx] for idx in [1,0,2]],fontsize=4,loc='upper center', bbox_to_anchor=(0.5, 1.1)) #loc='lower left')

axarr[0].set_ylabel('$I_{\\mathrm{NMDA}}$ (nA)',fontsize=6)


f.savefig("fig3B.pdf")

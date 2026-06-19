from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

f,axarr = subplots(6,1)

axarr[0].set_position([0.1,0.78,0.28,0.14])
axarr[1].set_position([0.48,0.78,0.3,0.14])
axarr[2].set_position([0.1,0.511,0.68,0.2])
axarr[3].set_position([0.1,0.356,0.68,0.12])
axarr[4].set_position([0.1,0.211,0.68,0.12])
axarr[5].set_position([0.1,0.066,0.68,0.12])

for iax in range(0,len(axarr)):
  axarr[iax].tick_params(axis='both', which='major', labelsize=4,direction='out',width=0.4,length=2)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  boxoff(axarr[iax])

axnew = f.add_axes([0.3,0.6,0.3,0.1])
axnew.tick_params(axis='both', which='major', labelsize=4,direction='out',width=0.4,length=2)
for axis in ['top','bottom','left','right']:
  axnew.spines[axis].set_linewidth(0.3)
boxoff(axnew)
    
cols = mytools.colorsredtolila(4,0.8)
cols_compartments = ['#000000','#888888','#AA0000','#FF00FF','#4444FF']
filename_seed1 = 'sim_onepyr_withSK3_seed1_withIh'


#Panel D: Comparison of predictions with Guan et al. data.
guandata = [[0, -0.5176876617774013],
            [100.39705048213273, 17.911993097497856],
            [199.6596710153148, 27.748058671268353],
            [300.6239364719229, 35.409836065573785],
            [399.8865570051049, 45.24590163934427],
            [502.5524673851389, 53.63244176013805],
            [600.1134429948952, 56.42795513373598],
            [701.6449234259785, 62.01898188093185],
            [800.3403289846852, 67.19585849870579],
            [900.737379466818, 75.89301121656601],
            [1001.1344299489507, 80.96635030198448]]
guandata_apamin = [[1.701644923425988, -0.31061259706642375],
                   [99.82983550765744, 22.05349439171701],
                   [198.52524106636417, 35.20276100086281],
                   [299.4895065229723, 48.766177739430546],
                   [399.3193420306296, 59.016393442622956],
                   [497.4475326148612, 64.81449525452977],
                   [599.5462280204198, 75.37532355478862],
                   [699.9432785025525, 81.06988783433995],
                   [800.3403289846852, 89.87057808455566],
                   [898.4685195689166, 97.84296807592753],
                   [1002.2688598979013, 103.12338222605695]]

givennames = ['withSK3_95apamin','withSK3_fullapamin','withSK3_apamin','withSK3']
amps = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
fs_all = []
for iSK in [0,1,2,3]:
  fs = []
  amps_plotted = []
  givenname = givennames[iSK]
  for iamp in range(0,len(amps)):
    if not exists('15dec31_1_'+givenname+'_longDC_amp'+str(amps[iamp])+'.mat'):
      print('15dec31_1_'+givenname+'_longDC_amp'+str(amps[iamp])+'.mat does not exist')
      continue
    print('Loading 15dec31_1_'+givenname+'_longDC_amp'+str(amps[iamp])+'.mat')
    A = scipy.io.loadmat('15dec31_1_'+givenname+'_longDC_amp'+str(amps[iamp])+'.mat')
    if len(A['spikes']) == 0:
      A['spikes'] = [A['spikes']]
    nSp = len([x for x in A['spikes'][0] if 1000 <= x < 4000])
    fs.append(nSp/3.0)
    amps_plotted.append(amps[iamp])
    print('nSp = '+str(nSp))
  fs_all.append(fs[:])
  if iSK == 3:
    axarr[1].plot(amps_plotted,fs,'r-',lw=0.5,color='#888888')
  elif iSK == 2:
    axarr[1].plot(amps_plotted,fs,'r-',lw=0.5,color='#FF9999')
  elif iSK == 1:
    polygon = Polygon(array([amps_plotted+amps_plotted[::-1], fs_all[0]+fs_all[1][::-1]]).T)
    p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
    p.set_facecolor('#FFEEEE')
    p.set_edgecolor(None)
    axarr[1].add_collection(p)
    

axarr[1].plot([guandata[i][0]/1000 for i in range(0,len(guandata))], [guandata[i][1] for i in range(0,len(guandata))], 'k--.', color='#000000',lw=0.4,ms=1,mew=1)
axarr[1].plot([guandata_apamin[i][0]/1000 for i in range(0,len(guandata_apamin))], [guandata_apamin[i][1] for i in range(0,len(guandata_apamin))], 'k--.', color='#FF0000',lw=0.4,ms=1,mew=1)

l = mytools.mylegend(f,[0.485,0.87,0.165,0.0375],['k--.','r--.'],['exp. data, apamin','exp. data, CTRL'],nx=1,dx=1,yplus=0.5,yplustext=0.35,colors=['#FF0000','#000000'],linewidths=[0.4,0.4],myfontsize=5.5)
try:
  l.get_children()[0].set_mew(1)
  l.get_children()[0].set_ms(1)
  l.get_children()[2].set_mew(1)
  l.get_children()[2].set_ms(1)
except:
  print('matplotlib issue, no mew/ms set')
l.set_xlim([0.1,3.4])
for axis in ['top','bottom','left','right']:
  l.spines[axis].set_visible(False)

l2 = mytools.mylegend(f,[0.61,0.785,0.165,0.0375],['k-','r-'],['sim. data, 98% block of SK','sim. data, CTRL'],nx=1,dx=1,yplus=0.5,yplustext=0.35,colors=['#FF9999','#888888'],linewidths=[0.4,0.4],myfontsize=5.5)
l2.set_xlim([0.1,3.4])
for axis in ['top','bottom','left','right']:
  l2.spines[axis].set_visible(False)


#Panel C: Comparison of predictions with Almog & Korngreen 2009 data.
print('Loading 15dec31_1_withSK3_VClamp_r0.02.mat')
A = scipy.io.loadmat('15dec31_1_withSK3_VClamp_r0.02.mat')
spikes = A['spikes']
times = A['times']
somaical = A['somaical']
somaican = A['somaican']
somaicat = A['somaicat']

recvdt = 0.1
axarr[0].plot([-10+recvdt*i for i in range(0,len(array(somaical)[0]))],1e3*array(somaical)[0],lw=0.5)
axarr[0].plot([-10+recvdt*i for i in range(0,len(array(somaical)[0]))],1e3*array(somaican)[0],lw=0.5)
axarr[0].plot([-10+recvdt*i for i in range(0,len(array(somaical)[0]))],1e3*array(somaicat)[0],lw=0.5)

cal_meancurr = mean(somaical[0][100:600])
can_meancurr = mean(somaican[0][100:600])
cat_meancurr = mean(somaicat[0][100:600])

axarr[0].text(-3,-18,'L-type AUC: '+'{:.1f}'.format(100*cal_meancurr/(cal_meancurr+can_meancurr+cat_meancurr))+'% (sim.), 47% (exp.)',color='#1f77b4',fontsize=5.5)
axarr[0].text(-3,-23,'N-type AUC: '+'{:.1f}'.format(100*can_meancurr/(cal_meancurr+can_meancurr+cat_meancurr))+'% (sim.), 27% (exp.)',color='#ff7f0e',fontsize=5.5)
axarr[0].text(-3,-28,'R-/T-type AUC: '+'{:.1f}'.format(100*cat_meancurr/(cal_meancurr+can_meancurr+cat_meancurr))+'% (sim.), 26% (exp.)',color='#2ca02c',fontsize=5.5)

#Panels E-H: 
binwidth = 1000
binwidth_cai = 10
bins_all = []
binsHCNm_all = []
binslcai_all = []
binscai_all = []
for iseed in range(1,11):
  filename = filename_seed1.replace('seed1','seed'+str(iseed))
  print('Loading '+filename)
  A = scipy.io.loadmat(filename)
  if iseed < 4 and len(A['spikes']) > 0:
    axarr[4].plot((A['spikes'][0]-10000)/1000,[-170*iseed for x in A['spikes'][0]],'r.',mew=0.4,lw=0.4,ms=0.4,color=cols[iseed])
    axarr[4].plot((A['times'][0]-10000)/1000,A['Vsoma'][0]-170*iseed,'b-',lw=0.1,color=cols[iseed],label='seed '+str(iseed))
    print(str(len(A['spikes'][0]))+' spikes')
  bins = zeros(int(100000/binwidth)+1,)
  for ispike in range(0,len(A['spikes'][0])):
    bins[int(A['spikes'][0][ispike]/binwidth)] = bins[int(A['spikes'][0][ispike]/binwidth)]+1
  binsHCNm_comps = []
  binslcai_comps = []
  binscai_comps = []
  for icomp in range(0,5): #0-4, soma, basal or apical dendrite
    binsHCNm = zeros(int(100000/A['times'][0][1]/binwidth)+1,)
    binslcai = zeros(int(100000/A['times'][0][1]/binwidth)+1,)
    binscai = zeros(int(1000/A['times'][0][1]/binwidth_cai)+1,)
    for it in range(0,int(len(A['times'][0])/binwidth)):
      binsHCNm[it] = mean(A['HCNm'][icomp][it*binwidth:(it+1)*binwidth])
      binslcai[it] = mean(A['lcai'][icomp][it*binwidth:(it+1)*binwidth])
      binscai[it] = mean(A['cai'][icomp][it*binwidth_cai:(it+1)*binwidth_cai])
    binsHCNm_comps.append(binsHCNm[:])
    binslcai_comps.append(binslcai[:])
    binscai_comps.append(binscai[:])
  bins_all.append(bins[:])
  binsHCNm_all.append(binsHCNm_comps[:])
  binslcai_all.append(binslcai_comps[:])
  binscai_all.append(binscai_comps[:])

polygon = Polygon(array([(binwidth*array(list(range(0,int(100000/binwidth)+1))+list(range(int(100000/binwidth),-1,-1)))-10000)/1000,
                         r_[mean(array(bins_all),axis=0)/(binwidth/1000)-std(array(bins_all),axis=0)/(binwidth/1000),mean(array(bins_all),axis=0)[::-1]/(binwidth/1000)+std(array(bins_all),axis=0)[::-1]/(binwidth/1000)]]).T)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
p.set_facecolor('#CCCCCC')
p.set_edgecolor(None)
axarr[5].add_collection(p)

axarr[5].plot((binwidth*array(range(0,int(100000/binwidth)+1))-10000)/1000,mean(array(bins_all),axis=0)/(binwidth/1000),lw=0.4,color='#000000',label='Average firing rate')

for icomp in range(0,5):
  compname = 'Soma' if icomp == 0 else ('Basal' if icomp == 1 else 'Apic '+str(icomp-1))
  axarr[2].plot((0.1*binwidth*array(range(0,len(binscai_all[0][0])))-10000)/1000,[mean([binscai_all[iseed][icomp][it] for iseed in range(0,len(binscai_all))]) for it in range(0,len(binscai_all[0][icomp]))],lw=0.4,color=cols_compartments[icomp],label=compname)
  axnew.plot((0.1*binwidth*array(range(0,len(binscai_all[0][0])))-10000)/1000,[mean([binscai_all[iseed][icomp][it] for iseed in range(0,len(binscai_all))]) for it in range(0,len(binscai_all[0][icomp]))],lw=0.4,color=cols_compartments[icomp],label=compname)
  axarr[3].plot((0.1*binwidth*array(range(0,len(binsHCNm_all[0][0])))-10000)/1000,[mean([binsHCNm_all[iseed][icomp][it] for iseed in range(0,len(binsHCNm_all))]) for it in range(0,len(binsHCNm_all[0][icomp]))],lw=0.4,color=cols_compartments[icomp],label=compname)

axarr[4].set_ylim([-670,90])
axarr[4].plot([4.500,4.500],[-30,70],'k-',lw=0.5)
axarr[4].text(5.000,20,'100 mV',va='center',ha='left',fontsize=5.5)
for axis in ['top','left','right']:
  axarr[4].spines[axis].set_visible(False)
axarr[4].set_yticks([])

axarr[1].set_xlabel('$I$ (nA)',fontsize=5.5)
axarr[0].set_xlabel('time (ms)',fontsize=5.5)
axarr[2].legend(fontsize=5.5,frameon=False)
axarr[3].legend(fontsize=5.5,ncol=3,frameon=False)
axarr[4].legend(fontsize=5.5,ncol=3,frameon=False)
axarr[5].set_xlabel('time (s)',fontsize=5.5)

axarr[1].set_ylabel('Firing rate\n(spikes/s)',fontsize=5.5)
axarr[0].set_ylabel('Ca$^{2+}$ current\n($\mu$A/cm$^2$)',fontsize=5.5)
axarr[2].set_ylabel('Cytosolic [Ca2+] (mM)',fontsize=5.5)
axarr[3].set_ylabel('HCN activation\nvariable (A.U.)',fontsize=5.5)
axarr[4].set_ylabel('Somatic memb.\npot.',fontsize=5.5)
axarr[5].set_ylabel('Firing rate\n(spikes/s)',fontsize=5.5)

axnew.set_xlim([-0.3,3.0])
axnew.set_ylabel('Cytosolic\n[Ca2+] (mM)',fontsize=5.5)
axnew.set_xlabel('time (s)',fontsize=5.5)

for iax in range(0,6):
  axarr[iax].set_xlim([-9.000,89.000])
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.09, pos.y1 - 0.01, chr(ord('C')+iax), fontsize=12)
axarr[1].set_xlim([0,1])
axarr[0].set_xlim([-5,55])

for iax in [2,3,4]:
  axarr[iax].set_xticklabels([])

for iax in range(0,5):
  if iax < 4:
    thisax = axarr[2+iax]
  else:
    thisax = axnew
  qs = thisax.get_ylim()
  polygon = Polygon(array([[0,0,0+0.1,0+0.1],[qs[0],qs[1],qs[1],qs[0]] if iax != 2 else [qs[0],-50,-50,qs[0]]]).T)
  print("iax = "+str(iax)+" qs = "+str(qs))
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
  p.set_facecolor('#CCCCCC')
  p.set_edgecolor(None)
  thisax.add_collection(p)
  
f.savefig("fig8C-H.pdf")

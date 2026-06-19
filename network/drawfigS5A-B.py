from pylab import *
import scipy.io
import mytools
from os.path import exists
import time
f,axarr = subplots(1,3)
Nsum = 100
Ncellsperpop = [582,97,97,291,291]

rateCoeffs = [0, 0.05, 0.1, 0.15, 0.2, 0.3]

axarr[0].set_position([0.09,0.8,0.25,0.15])
axarr[1].set_position([0.41,0.8,0.35,0.15])
axarr[2].set_position([0.6,0.915,0.14,0.03])

for iax in range(0,3):
  axarr[iax].tick_params(axis='both', which='major',labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)

FRmeans_all = []
errors = []
for irateCoeff in range(-1,len(rateCoeffs)):
  rateCoeff = rateCoeffs[irateCoeff]
  FRs_thisCoeff = []
  for myseed in range(1,11):
    if irateCoeff >= 0:
      filename = 'FRandMUA_sim_net_ExpWMnoNM_long_withSK3_CTRL_altgain_k'+str(rateCoeff).replace('.','_')+'_seed'+str(myseed)+'.mat'
    else:
      filename = 'FRandMUA_sim_net_long_withSK3_CTRL_altgain_seed'+str(myseed)+'.mat'

    if not exists(filename):
      print(filename)
      continue
    A = scipy.io.loadmat(filename)
    print('Loading '+filename)
    ipop = 3
    Nsum = 100

    FRs = zeros(int(len(A['FRvecs'][ipop])/Nsum)+1,)
    for ibin in range(0,len(FRs)):
      FRs[ibin] = sum(A['FRvecs'][ipop][Nsum*ibin:Nsum*(ibin+1)])/(Ncellsperpop[ipop]*Nsum/100) #Divided by the number of cells in the population, corrected by the bin size (if bin size is 1 s, i.e. Nsum=100, then we have the frequency in Hz without normalization, but if 10 s, i.e. Nsum=1000, then we have to also divide by 10 to get spikes/s)
    #  axarr[ipop].plot(10*Nsum*array(range(0,len(FRs))),FRs,lw=0.1,label=filename)
    FRs_thisCoeff.append(FRs[:])
  FRmeans_all.append(mean(array(FRs_thisCoeff),axis=0))
  if irateCoeff >= 0:
    errors.append(mean(abs(mean(array(FRs_thisCoeff),axis=0)[0:-2] - FRmeans_all[0][2:])))

axarr[0].plot(range(0,len(FRmeans_all[0])),FRmeans_all[0],'k-',lw=0.5)

yReal = FRmeans_all[0]
baseline = mean(yReal[-20:])
decayTimeConsts = [0.5*x for x in range(1,60)]
ks = [1+0.1*x for x in range(0,40)]
errors_decay_all = []
timenow = time.time()
for idecayTimeConst in range(0,len(decayTimeConsts)):
  if time.time() - timenow > 3:
    timenow = time.time()
    print('idecayTimeConst = '+str(idecayTimeConst)+'/'+str(len(decayTimeConsts)))
  decayTimeConst = decayTimeConsts[idecayTimeConst]
  errors_decay = []
  for ik in range(0,len(ks)):
    k = ks[ik]
    yProp = [baseline if t < 2 else baseline+baseline*k*exp(-t/decayTimeConst) for i,t in enumerate(10*Nsum*array(range(0,len(FRs)))/1000-10)]
    errors_decay.append(mean([abs(yProp[i]-yReal[i]) for i in range(0,len(yReal))]))
  errors_decay_all.append(errors_decay[:])
minErrors_Decay = [min(x) for x in errors_decay_all]
idecayTimeConst = argmin(minErrors_Decay)
decayTimeConst = decayTimeConsts[idecayTimeConst]
ik = argmin(errors_decay_all[idecayTimeConst])
k = ks[ik]
print("decayTimeConst = "+str(decayTimeConst)+" seconds, k = "+str(k))
axarr[0].plot(range(0,len(FRmeans_all[0])),[baseline if t<2 else baseline+baseline*k*exp(-t/decayTimeConst) for i,t in enumerate(10*Nsum*array(range(0,len(FRs)))/1000-10)],'r--',lw=0.5,dashes=(2,2))



axarr[1].plot(range(-2,-2+len(FRmeans_all[0])),FRmeans_all[0],'k-',lw=0.5)
cols = mytools.colorsredtolila(9,0.8)
axarr[2].plot(rateCoeffs,errors,'k-',lw=0.8)
for irateCoeff in range(0,len(rateCoeffs)):
  if rateCoeffs[irateCoeff] != 0.1:
    axarr[1].plot(range(0,len(FRmeans_all[0])),FRmeans_all[irateCoeff+1],'--',lw=0.5,color=cols[3+irateCoeff],dashes=(2,1))
  else:
    axarr[1].plot(range(0,len(FRmeans_all[0])),FRmeans_all[irateCoeff+1],'-',lw=0.5,color=cols[3+irateCoeff])
  axarr[2].plot(rateCoeffs[irateCoeff],errors[irateCoeff],'.',lw=1.6,mew=1.6,ms=1.6,color=cols[3+irateCoeff])

axarr[0].set_xlim([0,95])
axarr[1].set_xlim([0,95])
axarr[0].set_xlabel('time (s)',fontsize=6)
axarr[1].set_xlabel('time (s)',fontsize=6)
axarr[0].set_ylabel('FR (spikes/s)     ',fontsize=6)
axarr[1].set_ylabel('FR (spikes/s)     ',fontsize=6)
axarr[2].set_xlabel('rate coeff. (A.U.)',fontsize=5)
axarr[2].set_ylabel('avg. error\n(spikes/s)',fontsize=5,rotation=0,ha='right',va='center')
axarr[2].set_ylim([0.2,1.0])
axarr[2].set_yticks([0.3,0.6,0.9])

for iax in range(0,2):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.04, pos.y1 - 0.015, chr(ord('A')+iax), fontsize=9)

f.savefig("figS5A-B.pdf")


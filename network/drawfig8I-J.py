from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection

A = scipy.io.loadmat('data/sim_net_long_withSK_altgain_savespikes_seed1_ld.mat')

f,axarr = subplots(2,1)
axarr[0].set_position([0.1,0.5,0.8,0.4])
axarr[1].set_position([0.1,0.25,0.8,0.2])

#From simdat.py
from neuron import h
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, comment")
h("xwindows=1.0")
h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
h.load_file("labels.hoc") # has variables needed by network
from labels import *
def getcolors (lspks,lids,lctyID):
  lclrs = []
  for i in range(len(lids)):
    ty = lctyID[int(lids[i])]
    if IsLTS(ty): lclrs.append('b')
    elif ice(ty): lclrs.append('g')
    else: lclrs.append('r')
  return lclrs

Nskipped = 5 #TODO: implement this in the coloured version

lclrs = getcolors(A['ld']['lspks'][0][0][0],A['ld']['lids'][0][0][0],A['ld']['lctyID'][0][0][0])

for ipatch in range(0,3):
  if ipatch == 0:
    polygon = Polygon(array([[-3,-3,103,103],[0,258,258,0]]).T)
  elif ipatch == 1:
    polygon = Polygon(array([[-3,-3,103,103],[258,578,578,258]]).T)
  else:
    polygon = Polygon(array([[-3,-3,103,103],[578,778,778,578]]).T)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
  p.set_facecolor('#CCCCCC' if ipatch == 0 else ('#DDDDDD' if ipatch == 1 else '#EEEEEE'))
  p.set_edgecolor(None)
  axarr[0].add_collection(p)

for ipatch in range(0,4):
  if ipatch == 0:
    polygon = Polygon(array([[-3,-3,103,103],[0,96,96,0]]).T)
  elif ipatch == 1:
    polygon = Polygon(array([[-3,-3,103,103],[256,328,328,256]]).T)
  elif ipatch == 2:
    polygon = Polygon(array([[-3,-3,103,103],[401,448,448,401]]).T)
  else:
    polygon = Polygon(array([[-3,-3,103,103],[576,651,651,576]]).T)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
  p.set_facecolor('#CC8888' if ipatch == 0 else ('#DDAAAA' if ipatch <=2 else '#EECCCC'))
  p.set_edgecolor(None)
  axarr[0].add_collection(p)


sz = 1.0
#From simdat.py
lspks,lids,lctyID=A['ld']['lspks'][0][0][0][::Nskipped]/1e3,A['ld']['lids'][0][0][0][::Nskipped],A['ld']['lctyID'][0][0][0][::Nskipped]
#lclrs = getcolors(lspks,lids,lctyID)
axarr[0].scatter(lspks,lids,s=0.1*sz**2,c=lclrs[::Nskipped],marker='s',lw=0.4)
axarr[0].text(-0.8,130,'LAYER 6',ha='center',va='center',fontsize=6,rotation=90)
axarr[0].text(-0.8,400,'LAYER 5',ha='center',va='center',fontsize=6,rotation=90)
axarr[0].text(-0.8,670,'LAYER 2/3',ha='center',va='center',fontsize=6,rotation=90)

for iax in range(0,2):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.02, pos.y1 - 0.01, chr(ord('I')+iax), fontsize=12)
  axarr[iax].tick_params(axis='both', which='major', labelsize=5)
for axis in ['bottom']:
  axarr[0].spines[axis].set_linewidth(0.3)
for axis in ['left','right','top']:
  axarr[0].spines[axis].set_visible(False)
for axis in ['bottom','left','right','top']:
  axarr[1].spines[axis].set_linewidth(0.3)

axarr[0].get_xaxis().tick_bottom()

axarr[0].set_yticks([])
axarr[0].set_ylim([-1,778])
axarr[0].set_xlim([-1.6,99])
axarr[0].set_xlabel('time (s)',fontsize=6)

tstop = 100000
FRvecs = [zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,)]
for ispike in range(0,len(A['ld']['lspks'][0][0][0])):
  ivec = 0 if lclrs[ispike] == 'r' else (2 if lclrs[ispike] == 'g' else 1)
  FRvecs[ivec][int(A['ld']['lspks'][0][0][0][ispike]/10)] = FRvecs[ivec][int(A['ld']['lspks'][0][0][0][ispike]/10)] + 1 #Count number of spikes in the neuron class in each 10ms bin
  if ivec == 0:
    if A['ld']['lids'][0][0][0][ispike] < 96 or (A['ld']['lids'][0][0][0][ispike] >= 256 and A['ld']['lids'][0][0][0][ispike] < 328) or (A['ld']['lids'][0][0][0][ispike] >= 401 and A['ld']['lids'][0][0][0][ispike] < 448) or (A['ld']['lids'][0][0][0][ispike] >= 576 and A['ld']['lids'][0][0][0][ispike] < 651): #Stimulated half
      FRvecs[3][int(A['ld']['lspks'][0][0][0][ispike]/10)] = FRvecs[3][int(A['ld']['lspks'][0][0][0][ispike]/10)] + 1
    else: #Non-stimulated half
      FRvecs[4][int(A['ld']['lspks'][0][0][0][ispike]/10)] = FRvecs[4][int(A['ld']['lspks'][0][0][0][ispike]/10)] + 1

Nsum = 100 #Use bin size 1 s
Ncellsperpop = [582,97,97,290,292]
popstyles = ['r-','b-','g-','r-','m-']  #Changed after old4 - PV is green, LTS blue
poplabels = ['Pyr','LTS','PV+','Pyr, stimulated','Pyr, non-stimulated']
for ipop in range(1,5):
  FRs = zeros(int(len(FRvecs[ipop])/Nsum)+1,)
  for ibin in range(0,len(FRs)):
    FRs[ibin] = sum(FRvecs[ipop][Nsum*ibin:Nsum*(ibin+1)])/Ncellsperpop[ipop] #Number of spikes in each bin of 1 second divided by the number of cells gives us the average firing rate of a single neuron of the given pop
  axarr[1].plot(10*Nsum*array(range(0,len(FRs)))/1000,FRs,popstyles[ipop],lw=0.3,label=poplabels[ipop])

axarr[1].set_ylabel('Firing rate (spikes/s)',fontsize=6)
axarr[1].set_xlabel('time (s)',fontsize=6)
axarr[1].set_xlim([-1.6,99])
axarr[1].legend(fontsize=6)

f.set_size_inches(10,6)
print('#SPIKES PLOTTED = '+str(len(lspks)))
f.savefig("fig8I-J.pdf")


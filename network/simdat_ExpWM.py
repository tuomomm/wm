from pylab import *
import sys
import os
import string
from neuron import *
from datetime import datetime
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, comment")
h.simname=simname = "mtlhpc"
h.allfiles=allfiles = "geom.hoc pyinit.py geom_ExpWM.py mpisim_ExpWM.py"
h.simfiles=simfiles = "pyinit.py geom_ExpWM.py mpisim_ExpWM.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "10dec15"
h.output_file=output_file = "data/sim_net_ExpWM"
h.uname=uname = "x86_64"
h("templates_loaded=0")
templates_loaded=0
h("xwindows=1.0")
h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
CTYPi = 50.0; STYPi = 12.0
from pyinit import *
from neuron import h, gui
from nqs import *
from labels import *
import random
from pylab import *
import matplotlib.gridspec as gridspec
import configparser

Vector=h.Vector
newfigs=False
dodrawbounds=True
defbinsz=10
stimdel = 500 # delay from stim to measured firing rate peak response
myhost = os.uname()[1]
defnCPU=8 # this should match -np 8 , where 8 is number of cores used in mpiexec
print('Using ' , defnCPU, ' cpus by default.')

config = configparser.ConfigParser()

# determine config file name
def setfcfg ():
  fcfg = "netcfg_ExpWM.cfg" # default config file name
  for i in range(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  print("config file is " , fcfg)
  return fcfg

fcfg=setfcfg() # config file name
config.read(fcfg)

def conffloat (base,var): return float(config.get(base,var))
def confint (base,var): return int(config.get(base,var))
def confstr (base,var): return config.get(base,var)

tstart = 0 
tstop = conffloat('run','tstop') + tstart
binsz = conffloat("run","binsz") # bin size (in milliseconds) for MUA
recdt = conffloat('run','recdt')
recvdt = conffloat('run','recvdt')
vtslow=Vector(); vtslow.indgen(tstart,tstop-recdt,recdt); vtslow=numpy.array(vtslow.to_python())
vtfast=Vector(); vtfast.indgen(tstart,tstop-recvdt,recvdt); vtfast=numpy.array(vtfast.to_python())
useWM = confint('mem','useWM')
simstr = config.get('run','simstr')
baselinest = float(config.get("mem","baset")) #Baseline Start
nMem = confint('mem','nMem')
startt = conffloat('mem','startt')

def makefname (simstr,pcid): return 'data/' + simstr + '_pc_' + str(pcid) + '.npz'

ix,ixe=[1e9 for i in range(CTYPi)],[-1e9 for i in range(CTYPi)]
numc=[0 for i in range(CTYPi)]
allcells,icells,ecells=0,0,0

# reload the variables from config file
def reloadvars ():
  global tstart,loadstate,tstop,binsz,recdt,recvdt,vtslow,vtfast,useWM,simstr,baselinest,numc,allcells,ecells,icells,startt
  tstart = conffloat('run','loadtstop')
  loadstate = confint('run','loadstate')
  if loadstate == 0: tstart = 0 # only use previous end time if loading state
  tstop = conffloat('run','tstop') + tstart
  binsz = conffloat("run","binsz") # bin size (in milliseconds) for MUA
  recdt = conffloat('run','recdt')
  recvdt = conffloat('run','recvdt')
  vtslow=Vector(); vtslow.indgen(tstart,tstop-recdt,recdt); vtslow=numpy.array(vtslow.to_python())
  vtfast=Vector(); vtfast.indgen(tstart,tstop-recvdt,recvdt); vtfast=numpy.array(vtfast.to_python())
  useWM = confint('mem','useWM')
  simstr = config.get('run','simstr')
  baselinest = float(config.get("mem","baset")) #Baseline Start
  numc=[0 for i in range(CTYPi)]
  allcells,ecells,icells=0,0,0
  startt = conffloat('mem','startt')

# setup start/end indices for cell types
def makeix (lctyID):
  global ix,ixe,allcells,ecells,icells
  allcells,ecells,icells=0,0,0
  for i in range(CTYPi):
    ix[i]=1e9
    ixe[i]=-1e9
    numc[i]=0
  for i in range(len(lctyID)):
    ty = lctyID[i]
    numc[ty]+=1
    allcells+=1
    ix[ty] = min(ix[ty],i)
    ixe[ty] = max(ixe[ty],i)
    if h.ice(ty): icells+=1
    else: ecells+=1

lfastvar = ['volt', 'iAM', 'iNM', 'iGB', 'iGA','ina','ik','ica','ih']

# reads data from files saved by mpisim_ExpWM.py - must get nhost correct
def readdat (simstr,nhost):
  ld = {} # dict for data from a single host
  for pcid in range(nhost):
    fn = makefname(simstr,pcid)
    ld[pcid]=numpy.load(fn)
  return ld
  
# load/concat data from mpisim_ExpWM.py (which saves data from each host separately) - must get nhost correct
def loaddat (simstr,nhost,quiet=False):
  ld = readdat(simstr,nhost) # dict for data from a single host
  lids,lspks=array([]),array([]) # stitch together spike times
  for pcid in range(nhost):
    lids=concatenate((lids,ld[pcid]['lids']))
    lspks=concatenate((lspks,ld[pcid]['lspks']))
  ldout={} # concatenated output from diff hosts
  ldout['lids']=lids; ldout['lspks']=lspks;
  for k in ['lctyID','vt','lX', 'lY', 'lZ']: ldout[k]=ld[0][k][:]
  makeix(ldout['lctyID'])
  ncells = len(ld[0]['lctyID'])
  vt = ld[0]['vt'][:]
  locvarlist = []
  for loc in ['soma','Adend3']:
    for var in ['cai','Ihm','Ihp1','volt','caer','caCB','ip3cyt', 'iAM', 'iNM', 'iGB', 'iGA','ina','ik','ica','ih']:
      locvarlist.append(loc + '_' + var)
  locvarlist.append('Bdend_volt')
  for k in locvarlist:
    if not k in ld[0]: continue # skip?
    if not quiet: print(k)
    if k.endswith('volt') or lfastvar.count(k.split('_')[1])>0: sz = len(vtfast)
    else: sz = len(vtslow)
    ldout[k] = zeros((ncells,sz))
    ldk = ldout[k]
    for pcid in range(nhost):
      d = ld[pcid][k]
      lE = ld[pcid]['lE']; 
      for row,cdx in enumerate(lE): ldk[cdx,:] = d[row,0:sz]
      if k == 'soma_volt': # special case for soma voltage of I cells
        lI=ld[pcid]['lI']; dI=ld[pcid]['soma_voltI']; row=0
        for row,cdx in enumerate(lI): ldk[cdx,:] = dI[row,0:sz]
  for pcid in range(nhost):
    del ld[pcid].f
    ld[pcid].close()
    del ld[pcid]
  return ldout

# colors for spikes for raster -- based on cell type
def getcolors (lspks,lids,lctyID):
  lclrs = []
  for i in range(len(lids)):
    ty = lctyID[int(lids[i])]
    if IsLTS(ty): lclrs.append('b')
    elif ice(ty): lclrs.append('g')
    else: lclrs.append('r')
  return lclrs

# draw vertical lines indicating start/stop of a stimulus (units of seconds)
def drawBounds (ymin,ymax,row=0):
  if not dodrawbounds: return
  if not useWM: return
  nqm.tog("DB")
  MemONT = nqm.getcol("startt",row).x[row]/1e3
  MemOFFT = nqm.getcol("endt",row).x[row]/1e3
  axvline(MemONT,ymin,ymax,color='k',linewidth=3)
  axvline(MemOFFT,ymin,ymax,color='k',linewidth=3)

# get complement of cells in vid (all cells of any type in lty that's not in vid)
def getcomp (vid,ix,ixe,lty=[E2,E5R,E5B,E6]):
  vcomp = Vector()
  for ct in lty:
    for i in range(ix[ct],ixe[ct]+1):
      if i not in vid: vcomp.append(i)
  return vcomp

# get cells in vid if they're of specified type (in lty)
def getact (vid,ix,ixe,lty=[E2,E5R,E5B,E6]):
  vact = Vector()
  for ct in lty:
    for i in range(ix[ct],ixe[ct]+1):
      if i in vid: vact.append(i)
  return vact

# return a list of which cells get signal. each item in output corresponds to one row
# of nqm. first element in the item is list of cells that gets input and second element
# is its complement.
def getMemList (nq,ix,ixe):
  nq.tog("DB")
  memlist = [] #list of lists of which pyr cells to increase the input signal to.
  for i in range(int(nq.v[0].size())):
    vid = Vector()
    vid.copy(nq.get("vid",i).o[0])
    memlist.append([])
    memlist[-1].append(vid)
    memlist[-1].append(getcomp(vid,ix,ixe,lty=[E2,E5R,E5B,E6]))
  return memlist

# sets up NQS with memory stimuli
def setupnqm (fn=config.get("mem","nqm")):
  if os.path.exists(fn): # user specified a file that exists? then load it
    print('reading nqm from fn = ', fn)
    nqm = h.NQS(fn)
  elif os.path.exists('meminput/' + simstr + '_nqm.nqs'): # was the nqs saved from mpisim_ExpWM.py ? load it
    print('reading nqm from ' , 'meminput/' + simstr + '_nqm.nqs')
    nqm = h.NQS('meminput/' + simstr + '_nqm.nqs')
  else: # last resort - reconstruct from params in config file
    print('reconstructing nqm...')
    nqm = h.NQS('vid','startt','endt','rate','w'); nqm.odec('vid')
    lvid = []; nMem = confint('mem','nMem')
    memstartt = conffloat('mem','startt')
    memintert = conffloat('mem','intert')
    memdurt = conffloat('mem','durt')
    try:
      memfrac = conffloat('mem','frac')
    except:
      memfrac = 0.5
    memrate = conffloat('mem','rate')
    memW = conffloat('mem','weight')
    startt = memstartt; endt = startt+memdurt
    for i in range(nMem):
      lvid.append(h.Vector()); vtmp=Vector()
      for ct in [E2,E5R,E5B,E6]:
        if i > 0: vtmp.indgen(int(ix[ct]+numc[ct]*memfrac),ixe[ct],1)
        else: vtmp.indgen(ix[ct],int(ix[ct]+numc[ct]*memfrac-1),1)
        lvid[-1].append(vtmp)
      nqm.append(lvid[-1],startt,endt,memrate,memW)
      startt += (memintert+memdurt);
      endt = startt + memdurt
  memlist = getMemList(nqm,ix,ixe) # setup memlist
  nMem = len(memlist)
  return nqm,memlist,nMem

# return x with only two decimal places after x
def twop (x): return round(x,2)


# get an NQS with spikes (snq)
def getsnq (ld):
  lspks,lids,lctyID=ld['lspks'],ld['lids'],ld['lctyID']
  snq = NQS('id','t','ty','ice')
  snq.v[0].from_python(lids)
  snq.v[1].from_python(lspks)
  for i in range(len(lids)):
    snq.v[2].append(lctyID[int(lids[i])])
    snq.v[3].append(h.ice(lctyID[int(lids[i])]))
  ld['snq']=snq
  return snq





# make a multiunit activity vector from specified cell IDs
def getMUA (ld,IDs,tstart,tstop,binsz=5):
  #print 'getMUA: binsz=',binsz,',tstart=',tstart,',tstop=',tstop
  snq=None
  if 'snq' not in ld: getsnq(ld)
  snq=ld['snq']; snq.verbose=0; MUA=Vector()
  if type(IDs) == int:
    if IDs >= 0:
      for idx in range(ix[IDs],ixe[IDs]+1):
        if snq.select("id",idx) > 0:
          vt = snq.getcol("t")
          vh = vt.histogram(tstart,tstop,binsz)
          if MUA.size() < vh.size():
            MUA.copy(vh)
          else:
            MUA.add(vh)
    else: # all inhib cells, indicated with IDs < 0
      if snq.select("ice",1) > 0:
        vt = snq.getcol("t")
        vh = vt.histogram(tstart,tstop,binsz)
        if MUA.size() < vh.size():
          MUA.copy(vh)
        else:
          MUA.add(vh)
  elif type(IDs) == tuple:
    for ct in IDs:
      for idx in range(ix[ct],ixe[ct]+1):
        if snq.select("id",idx) > 0:
          vt = snq.getcol("t")
          vh = vt.histogram(tstart,tstop,binsz)
          if MUA.size() < vh.size(): MUA.copy(vh)
          else: MUA.add(vh)
  else:
    for idx in IDs: # list of cells to use
      if snq.select("id",idx) > 0:
        vt = snq.getcol("t")
        vh = vt.histogram(tstart,tstop,binsz)
        if MUA.size() < vh.size():
          MUA.copy(vh)
        else:
          MUA.add(vh)
  snq.tog("DB"); snq.verbose=1
  return MUA

def nzratVec (v1,v2):
  vout = Vector(v1.size())
  for i in range(int(v1.size())):
    if v2.x[i] != 0.0: vout.x[i] = v1.x[i] / v2.x[i]
  return vout

#
def ToHz (vmua,ncells,binsz):
  vmua.div(ncells)
  vmua.mul(1e3/binsz)

# make MUA and store in network as NQS and also as a dictionary
#  binsz is not global - passed in as an arg!
def makeMUA (ld,memlist,nMem,row,tstart,tstop,binsz,norm=True,lpops=[E5R,I5,I5L]):
  if 'nqMUA' in ld: nqsdel(ld['nqMUA'])
  vact,vna = memlist[row][0],memlist[row][1] # IDs of activated and non-activated pyramidal cells
  nq=NQS("row","tys","ty","vMUA"); nq.odec("vMUA"); nq.strdec("tys")
  ld['nqMUA']=nq; dMUA=[]
  diffMUA,ratMUA = Vector(),Vector()
  dMUA.append({})
  for i in lpops:
    vmua = getMUA(ld,i,tstart,tstop,binsz)
    if type(i) == tuple: 
      ncells = sum(numc[j] for j in i)
      nq.append(row,CTYP[i[0]],sum(j for j in i),vmua)
    else: 
      ncells = numc[i]
      nq.append(row,CTYP[i],i,vmua)
    if norm: ToHz(vmua,ncells,binsz)
    dMUA[-1][i] = Vector()
    dMUA[-1][i].copy(vmua)
    if type(i) != tuple: dMUA[-1][i].label("M"+str(row)+"_"+CTYP[i])
  vmua = getMUA(ld,-1,tstart,tstop,binsz) # a MUA for all I cells (-1 to getMUA indicates all I cells)
  if norm: ToHz(vmua,icells,binsz)
  dMUA[-1][-666] = Vector() # inhib cells indicated by -666
  dMUA[-1][-666].copy(vmua)
  dMUA[-1][-666].label("M"+str(row)+"_"+"I")
  nq.append(row,"I",-666,vmua) # inhib cells indicated by "I" and -666
  v1,v2=Vector(),Vector()
  vmua = getMUA(ld,vact,tstart,tstop,binsz) # activated E cell MUA
  if norm: ToHz(vmua,vact.size(),binsz) # rates
  diffMUA.copy(vmua); v1.copy(vmua)
  dMUA[-1][-1] = Vector()
  dMUA[-1][-1].copy(vmua)
  nq.append(row,"act",-1,vmua) # activated pyramidal cell MUA
  vmua = getMUA(ld,vna,tstart,tstop,binsz) # nonactivated E cell MUA
  if norm: ToHz(vmua,vna.size(),binsz) # rates
  diffMUA.sub(vmua); v2.copy(vmua)
  ratMUA0,ratMUA1 = nzratVec(v1,v2),nzratVec(v2,v1)
  dMUA[-1][-2] = Vector()
  dMUA[-1][-2].copy(vmua)
  nq.append(row,"nonact",-2,vmua) # non-activated pyramidal cell MUA
  dMUA[-1][-1].label("M"+str(row)+"_PyrMUA: Activated")
  dMUA[-1][-2].label("M"+str(row)+"_PyrMUA: NOT_Activated")
  ld['dMUA'+str(binsz)]=dMUA; ld['diffMUA'+str(binsz)]=diffMUA;

ld,nqm,memlist,nMem=[None for i in range(4)]

def drawraster (ld,sz=2):
  lspks,lids,lctyID=ld['lspks']/1e3,ld['lids'],ld['lctyID']
  try:
    lclrs = ld['lclrs']
  except:    
    lclrs = getcolors(lspks,lids,lctyID)
    ld['lclrs'] = lclrs

  import scipy.io
  FRvecs = [zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,),zeros(int(tstop/10)+1,)]
  #limits = [0,cpernet[E6]/2],[cpernet[E6]+cpernet[I6]+cpernet[I6L],cpernet[E6]+cpernet[I6]+cpernet[I6L]+cpernet[E5B]/2],[cpernet[E6]+cpernet[I6]+cpernet[I6L]+cpernet[E5B],cpernet[E6]+cpernet[I6]+cpernet[I6L]+cpernet[E5B]+cpernet[E5R]/2],[cpernet[E6]+cpernet[I6]+cpernet[I6L]+cpernet[E5B]+cpernet[E5R]+cpernet[I5]+cpernet[I5L],cpernet[E6]+cpernet[I6]+cpernet[I6L]+cpernet[E5B]+cpernet[E5R]+cpernet[I5]+cpernet[I5L]+cpernet[E2]/2]
  for ispike in range(0,len(lspks)):
    ivec = 0 if lclrs[ispike] == 'r' else (1 if lclrs[ispike] == 'b' else 2)
    FRvecs[ivec][int(lspks[ispike]*100)] = FRvecs[ivec][int(lspks[ispike]*100)] + 1
    if ivec == 0:
      if lids[ispike] < 96 or (lids[ispike] >= 256 and lids[ispike] < 328) or (lids[ispike] >= 401 and lids[ispike] < 448) or (lids[ispike] >= 576 and lids[ispike] < 651): #Stimulated half
        FRvecs[3][int(lspks[ispike]*100)] = FRvecs[3][int(lspks[ispike]*100)] + 1
      else: #Non-stimulated half
        FRvecs[4][int(lspks[ispike]*100)] = FRvecs[4][int(lspks[ispike]*100)] + 1
  row = 0
  makeMUA(ld,memlist,nMem,row,tstart,tstop,10.0)
  MUAkeys = list(ld['dMUA10.0'][0].keys())
  MUAs = [array(ld['dMUA10.0'][0][MUAkeys[i]]) for i in range(0,len(MUAkeys))]
  scipy.io.savemat('FRandMUA_'+simstr+'.mat', {'FRvecs': FRvecs, 'MUAkeys': MUAkeys, 'MUAs': MUAs})


# draw sim data (designed for baseline sim but works with others too if their data is loaded)
def baseDraw (drawty='raster',BINSZ=defbinsz,mint=5e3/1e3,maxt=20e3/1e3,setyl=True,skipdiff=True,lbinsz=[10,100],nrow=-1,noff=0,lyr=-1):
  global binsz; idx=1; locs=['soma', 'Adend3']; tx,ty=-.15,1.06; naxbins=5; fsz=18
  drawraster(ld,sz=2)

def mydraw ():
  baseDraw('rastmua')

ld=loaddat(simstr,defnCPU); nMem=0;  nqm,memlist,nMem=setupnqm();
mydraw()



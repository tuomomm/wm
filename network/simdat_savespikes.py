import sys
import os
import string
from neuron import *
from datetime import datetime
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, comment")
h.simname=simname = "mtlhpc"
h.allfiles=allfiles = "geom.hoc pyinit.py geom.py mpisim.py"
h.simfiles=simfiles = "pyinit.py geom.py mpisim.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "10dec15"
h.output_file=output_file = "data/15dec31.1"
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
import scipy.io

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
  fcfg = "netcfg.cfg" # default config file name
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

# reads data from files saved by mpisim.py - must get nhost correct
def readdat (simstr,nhost):
  ld = {} # dict for data from a single host
  for pcid in range(nhost):
    fn = makefname(simstr,pcid)
    ld[pcid]=numpy.load(fn)
  return ld
  
# load/concat data from mpisim.py (which saves data from each host separately) - must get nhost correct
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
		
#
def drawraster (ld,sz=2):
  if newfigs: figure();
  lspks,lids,lctyID=ld['lspks']/1e3,ld['lids'],ld['lctyID']
  try:
    lclrs = ld['lclrs']
  except:    
    lclrs = getcolors(lspks,lids,lctyID)
    ld['lclrs'] = lclrs
  scatter(lspks,lids,s=0.1*sz**2,c=lclrs,marker='.',lw=0.5)
  allcells = len(lctyID)
  for row in range(nMem): drawBounds(0,allcells,row)
  xlim((tstart/1e3,tstop/1e3)); ylim((0,allcells)); tight_layout(); xlabel('Time (s)',fontsize=18); title('raster');

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
  elif os.path.exists('meminput/' + simstr + '_nqm.nqs'): # was the nqs saved from mpisim.py ? load it
    print('reading nqm from ' , 'meminput/' + simstr + '_nqm.nqs')
    nqm = h.NQS('meminput/' + simstr + '_nqm.nqs')
  else: # last resort - reconstruct from params in config file
    print('reconstructing nqm...')
    nqm = h.NQS('vid','startt','endt','rate','w'); nqm.odec('vid')
    lvid = []; nMem = confint('mem','nMem')
    memstartt = conffloat('mem','startt')
    memintert = conffloat('mem','intert')
    memdurt = conffloat('mem','durt')
    memfrac = conffloat('mem','frac')
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

# get an array of times of interest
def gettimes (nqm,rampt=3e3):
  nqm.tog("DB")
  times = [] # periods of interest
  vstartt,vendt = Vector(),Vector()
  vstartt.copy(nqm.getcol("startt"))
  vendt.copy(nqm.getcol("endt"))
  if nMem > 0:
    times = [(baselinest,vstartt[0])]
    times.append((vstartt[0],vendt[0]))
    times.append((vendt[0],vendt[0]+rampt))
    for row in range(1,nMem):
      times.append((vendt[row-1]+rampt,vstartt[row]))
      times.append((vstartt[row],vendt[row]))
      times.append((vendt[row],vendt[row]+rampt))
    times.append((vendt[-1]+rampt,tstop))
  return times

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

# print spike rates in each population in the periods of interest (baseline, signal, recall)
def prspkcount (snq,nqm,rampt=3e3,times=None,quiet=False):
  snq.verbose=0
  if times is None:
    times = gettimes(nqm,rampt) # periods of interest
  lfa,lfn=[],[]
  for row in range(nMem):
    lfa.append([]); lfn.append([]);
    vact,vna = memlist[row][0],memlist[row][1]
    for (startt,endt) in times: 
      fa = round(1e3*snq.select("id","EQW",vact,"t","[]",startt,endt) / ((endt-startt)*vact.size()),2)
      lfa[-1].append(fa)
      try:
        fn = round(1e3*snq.select("id","EQW",vna,"t","[]",startt,endt) / ((endt-startt)*vna.size()),2)
      except:
        print(" cant assign fn")
        fn = 0.0
      lfn[-1].append(fn)
      if not quiet:
        stro = "t=" + str(startt) + "-" + str(endt) + ". ActE:" + str(fa) + " Hz. NonActE:" + str(fn) + " Hz."
        for ct in [I2, I2L, I5, I5L, I6, I6L]:
          fb = round(1e3*snq.select("id","[]",ix[ct],ixe[ct],"t","[]",startt,endt) / ((endt-startt)*numc[ct]),2)
          stro += ' ' + CTYP[ct] + ':' + str(fb)
        fi = round(1e3*snq.select("ice",1,"t","[]",startt,endt) / ((endt-startt)*icells),2)
        stro += ' I:' + str(fi)
        for ct in [E2, E5R, E5B, E6]:
          fb = round(1e3*snq.select("id","[]",ix[ct],ixe[ct],"t","[]",startt,endt) / ((endt-startt)*numc[ct]),2)
          stro += ' ' + CTYP[ct] + ':' + str(fb)
        fe = round(1e3*snq.select("ice",0,"t","[]",startt,endt) / ((endt-startt)*ecells),2)
        stro += ' E:' + str(fe)
        print(stro)
  snq.verbose=1
  return lfa,lfn

# print spike counts & draw activated and non-activated rates
def drawrat (ld,incsz=500.0,winsz=1e3,quiet=False):
  snq=getsnq(ld);
  times=[]
  for i in range(2*int(tstop/winsz)):times.append( (i*incsz,(i*incsz)+winsz) )
  lfa,lfn = prspkcount(snq,nqm,rampt=1e3,times=times,quiet=quiet)
  rat = []; tt=[];
  for i in range(len(lfa[0])): rat.append(lfa[0][i]/lfn[0][i])
  for tup in times: tt.append(0.5*(tup[0]+tup[1])) # time midpoints
  plot(tt,rat,'b',linewidth=2);plot(tt,rat,'b.',markersize=10);tight_layout();show() #

# getfnq - make an NQS with ids, firing rates, types
def getfnq(snq,lctyID,skipms=200):
  snq.verbose=0; snq.tog("DB"); 
  fnq = h.NQS("id","freq","ty")
  tf = tstop - skipms # duration we're considering for frequency calc
  for i in range(allcells):
    n = float( snq.select("t",">",skipms,"id",i) ) # number of spikes
    fnq.append(i, n*1e3/tf, lctyID[i])
  snq.verbose=1
  return fnq

# pravgrates - print average firing rates using fnq
def pravgrates(fnq,skipms=200):
  fnq.verbose=0
  ty=0; tf=float( tstop - skipms )
  for ty in range(CTYPi):
    if numc[ty] < 1: continue
    fnq.select("ty",ty)
    vf = fnq.getcol("freq")
    if vf.size() > 1: print(CTYP[ty], " avg rate = ", vf.mean(), "+/-", vf.stderr(), " Hz")
    else: print(CTYP[ty], " avg rate = ", vf.mean(), "+/-", 0.0 , " Hz")
    ty += 1
  fnq.verbose=1

#
def drawcellcal (ld,cdx):
  lk = ['Adend3_cai', 'Adend3_ip3cyt', 'Adend3_caer', 'soma_cai', 'soma_ip3cyt', 'soma_caer']
  for i,k in enumerate(lk):
    plot(vtslow/1e3,ld[k][cdx,:],linewidth=2)
  legend(lk,loc='best')
  xlabel('Time (s)'); ylabel('i');

# draws synaptic currents
def drawcellcurr (ld,cdx):
  #clrs = ['k','b','g','r']
  subplot(1,3,1); # synaptic currents
  lk = ['Adend3_iAM', 'Adend3_iNM', 'Adend3_iGA', 'Adend3_iGB', 'soma_iGA']
  for i,k in enumerate(lk): plot(vtfast/1e3,ld[k][cdx,:],linewidth=2)
  legend(lk,loc='best'); xlabel('Time (s)'); ylabel('i');
  for g,sec in zip([2,3],['Adend3','soma']):
    subplot(1,3,g)
    li = [sec + s for s in ['_ik','_ina','_ica','_ih']]
    for k,clr in zip(li,['b','g','r', 'c']): plot(vtfast/1e3,ld[k][cdx,:],clr,linewidth=2)
    legend(['ik','ina','ica','ih'],loc='best'); xlabel('Time (ms)'); ylabel(sec + '_i')

#
def drawcellvolt (ld,cdx,clr='g',ln=1,sec='soma'):
  vt = vtfast; sv=ld[sec+'_volt'][cdx,:]
  plot(vt,sv,clr,linewidth=ln)
  xlim((tstart,tstop))

# gets min,max value from the cell on the key. ld is from loaddat
def getminmaxy (ld,cellidx,loc,key,sidx,eidx):
  da=ld[loc+'_'+key]
  miny=min(da[cellidx,sidx:eidx])
  maxy=max(da[cellidx,sidx:eidx])
  return miny,maxy

# make integrated Ih p1 for the populations of pyramidal cells (by section)
# lyr==-1 means all layers; lyr==2,5,6 means only E cells from that layer
def makeIntegMem (ld,var,memlist,nMem,lyr=-1):
  d = []
  for row in range(nMem):
    vact,vna = memlist[row][0],memlist[row][1] # IDs of activated and non-activated pyramidal cells
    if lyr == 2:
      vact = getact(vact,ix,ixe,lty=[E2])
      vna = getact(vna,ix,ixe,lty=[E2])
    elif lyr == 5:
      vact = getact(vact,ix,ixe,lty=[E5R,E5B])
      vna = getact(vna,ix,ixe,lty=[E5R,E5B])
    elif lyr == 6:
      vact = getact(vact,ix,ixe,lty=[E6])
      vna = getact(vna,ix,ixe,lty=[E6])
    d.append({})
    for s in ["soma", "Adend3"]:
      d[-1][s+"_act_"+var] = getInteg(ld,var,s,vact)
      d[-1][s+"_nonact_"+var] = getInteg(ld,var,s,vna)
  ld['d'+var]=d
  return d

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

#
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
  ld['ratMUA0'+str(binsz)],ld['ratMUA1'+str(binsz)]=ratMUA0,ratMUA1

# draws the MUA on top of the raster plot -- binsz is in milliseconds, and
#  should be same value as used when calling makeMUA
def drawMUA (ld,memlist,nMem,row,binsz,ipop=[],skipdiff=False,plotrows=1,plotoff=0):
  try: 
    dMUA=ld['dMUA'+str(binsz)]; diffMUA=ld['diffMUA'+str(binsz)]
  except:
    makeMUA(ld,memlist,nMem,row,tstart,tstop,binsz)
    dMUA,diffMUA=ld['dMUA'+str(binsz)],ld['diffMUA'+str(binsz)];
  if newfigs: figure(); 
  lwidth=1;
  if binsz >= 100: lwidth=4;
  title('MUA'); tt=linspace(tstart/1e3,tstop/1e3,int(dMUA[0][-1].size())); subplot(plotrows,1,1+plotoff)
  for i in ipop: plot(tt,dMUA[row][i].as_numpy(),'g',linewidth=lwidth)
  plot(tt,dMUA[row][-1].as_numpy(),'c',linewidth=lwidth);
  plot(tt,dMUA[row][-2].as_numpy(),'m',linewidth=lwidth)
  if nMem > 1 and not skipdiff: plot(tt,diffMUA.as_numpy(),'y',linewidth=lwidth)
  xlim((tstart/1e3,tstop/1e3)); grid(True);
  for r in range(nMem): drawBounds(-10,10,r)
  plot(tt,dMUA[row][-666].as_numpy(),color=(0.5,0.5,0.5),linewidth=lwidth) # all I cells
  xlabel('Time (s)',fontsize=18);
  xlim((tstart/1e3,tstop/1e3));
  print(tstart/1e3,tstop/1e3)
  for r in range(nMem): drawBounds(-10,10,r)
  tight_layout(); 
###########################################################3

# reload variables, reset config file name
def myreload (Simstr,ncpu=defnCPU):
  global fcfg,simstr
  print('simstr is ' , simstr)
  simstr = Simstr
  print('simstr is ' , simstr)
  fcfg = 'backupcfg/' + simstr + '.cfg'
  if len(config.read(fcfg))==0:
    print('myreload ERRA: could not read ', fcfg, '!!!')
    return None
  reloadvars()
  print(fcfg)
  try:
    ld=loaddat(simstr,ncpu);
    nqm,memlist,nMem=setupnqm();
    return ld,nqm,memlist,nMem
  except:
    print('could not reload data from' , Simstr)
    return None

ld,nqm,memlist,nMem=[None for i in range(4)]

# draw sim data (designed for baseline sim but works with others too if their data is loaded)
def baseDraw (drawty='raster',BINSZ=defbinsz,mint=5e3/1e3,maxt=20e3/1e3,setyl=True,skipdiff=True,lbinsz=[10,100],nrow=-1,noff=0,lyr=-1):
  global binsz; idx=1; locs=['soma', 'Adend3']; tx,ty=-.15,1.06; naxbins=5; fsz=18
  if drawty == 'rastmua':
    ax=subplot(2,1,1);ax.locator_params(nbins=naxbins);
    dodrawbounds = True
    drawraster(ld,sz=2); ylabel('Cell ID'); title('');xlabel('');
    xlim((mint,maxt))
    dodrawbounds = False
    for BINSZ in lbinsz:
      drawMUA(ld,memlist,nMem,0,BINSZ,ipop=[],skipdiff=skipdiff,plotrows=2,plotoff=1); 
      ax=subplot(2,1,2); ax.locator_params(nbins=naxbins); ylabel('MUA (Hz)',fontsize=fsz); title('');
      text(tx,ty,'b',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
      xlim((mint,maxt));grid(True);
      #ylim((0,50));
  elif drawty == 'mua':
    print('mua' , ' binsz : ' , binsz , ' BINSZ: ' , BINSZ)
    try: 
      if binsz != BINSZ: del ld['dMUA'],ld['diffMUA'],ld['nqMUA'];
    except: pass
    binsz=BINSZ;
    plotrows=2
    if nrow > 0: plotrows=nrow
    drawMUA(ld,memlist,nMem,0,BINSZ,plotrows=plotrows,ipop=[],skipdiff=skipdiff,plotoff=noff); 
    txstr=['a','b']; msty = ['E', 'I']
    for i in range(2):
      ax=subplot(plotrows,1,i+1); ax.locator_params(nbins=naxbins); ylabel(msty[i]+' MUA (Hz)'); title('');
      text(tx,ty,txstr[i],fontsize=14,fontweight='bold',transform=ax.transAxes);
      xlim((mint,maxt));grid(True);

def mydraw ():
  baseDraw('rastmua')
  for i in range(2,3,1): subplot(2,1,i); #ylim((0,45))  
  subplot(2,1,1);title('');
  tight_layout()
  ltxt=['a','b']; fsz=20
  # need to draw letters on after tight_layout or they won't fit/display properly
  for i in range(2):
    ax=subplot(2,1,i+1);
    tx,ty=-0.0225,0.95
    text(tx,ty,ltxt[i],fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  #get rid of vertical ticks and ylabels for raster
  ax=subplot(2,1,1); ax.grid(True); ax.set_yticklabels([]); ax.set_yticks([]); ax.set_ylabel('')
  fsz=18
  # and now draw text indicating layer
  tx,ty=-0.04,0.825
  text(tx,ty,'L2/3',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  tx,ty=-0.04,0.45
  text(tx,ty,'L5',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  tx,ty=-0.04,0.12
  text(tx,ty,'L6',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  fsz=20
  subplot(2,1,2); ylabel('MUA (Hz)',fontsize=fsz,fontweight='bold'); xlabel('Time (s)',fontsize=fsz,fontweight='bold')
  for i in range(1,3,1): subplot(2,1,i); xlim((8,14))
  ax.figure.savefig(output_file+'.eps')

ld=loaddat(simstr,defnCPU); nMem=0;  nqm,memlist,nMem=setupnqm();
scipy.io.savemat(output_file+'_ld.mat',{'ld': ld})
ion()
mydraw()



from neuron import *
import sys
import os
import string
h("strdef simname, allfiles, simfiles, output_file, datestr, uname, comment")
h.simname=simname = "onepyr"
h.allfiles=allfiles = "pyinit.py geom_withSK3.py onepyr.py"
h.simfiles=simfiles = "pyinit.py geom_withSK3.py onepyr.py"
h("runnum=1")
h.datestr=datestr = "10dec15"
h.output_file=output_file = "data/15dec31.2"
h.uname=uname = "x86_64"
h("templates_loaded=0")
h("xwindows=1.0")
h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
h("proc setMemb () { }") # so e_pas will not get modified
from neuron import h
from pyinit import *
from labels import *
from geom_withSK3 import *
from nqs import *
import random
from pylab import *
from datetime import datetime
from neuron import rxd
import shutil
import scipy.io
import mytools
Vector = h.Vector
# nsubseg == number of subsegments per segment, subsegum == subseg size in microns
# nsubseg > 0 takes precedence over subsegum
nsubseg,subsegum = dconf['nsubseg'],dconf['subsegum'] 

if len(sys.argv) > 1:
  for j in range(1,len(sys.argv)):
    print("""dconf['"""+sys.argv[j].split('=')[0]+"""'] = """+sys.argv[j].split('=')[1])
    dconf[sys.argv[j].split('=')[0]] = sys.argv[j].split('=')[1]

for k in ['tstop','iseed','pseed','wseed']:
  if type(dconf[k]) == str:
    dconf[k] = float(dconf[k])
for k in ['simstr', 'saveout', 'saveoutSpikes', 'recdt', 'recvdt', 'verbose', 'indir', 'outdir']:
  if type(dconf[k]) == str: exec(k + ' = \'' + str(dconf[k]) + '\'')
  else: exec(k + ' = ' + str(dconf[k]))  

# make dir, catch exceptions
def safemkdir (dn):
  try:
    os.mkdir(dn)
    return True
  except OSError:
    if not os.path.exists(dn):
      print('could not create', dn)
      return False
    else:
      return True

safemkdir('data')
safemkdir('meminput')

# backup the config file
def backupcfg (simstr):
  safemkdir('backupcfg')
  fout = 'backupcfg/' + simstr + '.cfg'
  if os.path.exists(fout): os.remove(fout)  
  try:
    shutil.copy2(fcfg, fout) # fcfg created in geom.py via conf.py
  except:
    pass

backupcfg(simstr) # backup the config file

# return True if s exists in file fname
def grepstr (fname, s):
  found = False
  try:
    fp = open(fname,'r')
    lns = fp.readlines()
    for ln in lns:
      if ln.count(s) > 0:
        found = True
        break
    fp.close()
  except:
    pass
  return found

h.tstop = tstop = dconf['tstop']
tstart = 0 # only use previous end time if loading state
h.dt = dconf['dt']
h.steps_per_ms = 1/h.dt
h.celsius = 37
h.fracca_MyExp2SynNMDABB = dconf['nmfracca'] # fraction of NMDA current that is from calcium
wmatex = numpy.zeros( ( STYPi,) ) # external weights
ratex = numpy.zeros( (STYPi,) )  # external rates
mGLURR = dconf['mGLURR'] # ratio of mGLUR weights to AM2 weights
mGLURRWM = dconf['mGLURRWM']; 
NMAMRWM = dconf['NMAMRWM']; AMRWM = dconf['AMRWM'] # WM related weight ratios
GARWM = dconf['GARWM']
EXGain = dconf['EXGain'] # gain for noise inputs , default 1.0
sgrhzEE = dconf['sgrhzEE'] # external E inputs to E cells; 1000 is default
sgrhzEI = dconf['sgrhzEI'] # external E inputs to I cells
sgrhzIE = dconf['sgrhzIE'] # external I inputs to E cells
sgrhzII = dconf['sgrhzII'] # external I inputs to I cells
sgrhzNME = dconf['sgrhzNME'] # external NM inputs to E cells; 10 is default
sgrhzMGLURE = dconf['sgrhzMGLURE'] # external mGLUR inputs to E cells
sgrhzGB2 = dconf['sgrhzGB2'] # external inputs onto E cell GB2 synapses

ratex[AM2]=sgrhzEE
ratex[NM2]=sgrhzNME
ratex[GA]=sgrhzIE
ratex[GA2]=sgrhzIE
ratex[GB2]=sgrhzGB2
wmatex[AM2] = 0.02e-3 
wmatex[NM2] = 0.02e-3 
wmatex[GA] =  0.2e-3 
wmatex[GA2] = 0.2e-3 
wmatex[GB2] = 5e-3
for sy in range(STYPi): wmatex[sy] *= EXGain # apply gain control

cell = PyrAdr(0,0,0,0,E5R)
ce = [cell]
rxdsec = [s for s in cell.all_sec] # Section list for use with rxd

#
def setnsubseg ():
  global mGLURRWM
  if nsubseg > 0:
    print('setting rxd nsubseg to ', dconf['nsubseg'])
    rxd.set_solve_type(nsubseg=dconf['nsubseg'])
    mGLURRWM /= float(nsubseg)
    print('adjusted mGLURRWM to ' , mGLURRWM)
  elif subsegum > 0.0:
    for sec in cell.all_sec:
      ns = int(sec.L / subsegum)
      if ns > 1:
        if ns % 2 == 0: ns += 1
        print('L = ', sec.L, ' subseg = ', ns)
        rxd.set_solve_type(sec, nsubseg=ns)
        if sec == cell.Adend3:
          mGLURRWM /= float(ns)
          print('adjusted mGLURRWM to ' , mGLURRWM)
      else:
        print('L = ', sec.L, ' not using subseg')
  rxd.options.subseg_interpolation = 1
  rxd.options.subseg_averaging = 1

if nsubseg > 0 or subsegum > 0.0: setnsubseg()

# get a recording vector for the synaptic current
def recveccurr (sy):
  try:
    return saferecord(sy.syn._ref_i, recvdt)
  except:
    return saferecord(sy.syn._ref_iNMDA, recvdt)

cell.dsy = {}; dsstr = ['AMPA', 'NMDA', 'GABAS', 'mGLUR', 'GABAB']
cell.vsy = {}
#
def addsyns ():
  ds = {}; ds[cell.Adend3]='Adend3'; ds[cell.Adend2]='Adend2'; ds[cell.Adend1]='Adend1'; ds[cell.Bdend]='Bdend'
  for sec in [cell.Adend3, cell.Adend2, cell.Adend1, cell.Bdend]:
    llsy = [];
    nloc = sec.nseg
    llvsy = []; # for recording currents
    for i,seg in enumerate(sec):
      if seg.x == 0.0 or seg.x == 1.0: continue # skip endpoints
      lsy = []; loc = seg.x; lvsy = [] #AMPA, NMDA, GABAA_slow, GABAB
      lsy.append(Synapse(sect=sec,loc=loc,tau1=0.05,tau2=5.3,e=0)); lvsy.append(recveccurr(lsy[-1]))#AMPA
      lsy.append(SynapseNMDA(sect=sec,loc=loc,tau1NMDA=15,tau2NMDA=150,r=1,e=0)) # NMDA
      lvsy.append(recveccurr(lsy[-1]))
      lsy.append(Synapse(sect=sec,loc=loc,tau1=0.2,tau2=20,e=-80)) # GABAA_slow
      lvsy.append(recveccurr(lsy[-1]))
      lsy.append(SynapsemGLUR(sect=sec,loc=loc)) # mGLUR
      for node in ip3.nodes(seg): node.include_flux(lsy[-1].syn._ref_rip3 ) # all the sub-segments get flux
      lsy.append(SynapseGABAB(sect=sec,loc=loc)) # GABAB
      lvsy.append(recveccurr(lsy[-1]))
      llsy.append(lsy); llvsy.append(lvsy)
    cell.dsy[sec] = llsy; cell.vsy[sec] = llvsy
  sec = cell.soma; llsy = []; nloc = sec.nseg; llvsy = []
  for i,seg in enumerate(sec):
    if seg.x == 0.0 or seg.x == 1.0: continue # skip endpoints
    lsy = []; loc = seg.x; lvsy = []
    lsy.append(Synapse(sect=sec,loc=loc,tau1=0.07,tau2=9.1,e=-80)) # GABAA_fast
    lvsy.append(recveccurr(lsy[-1]))
    lsy.append(Synapse(sect=sec,loc=loc,tau1=0.05,tau2=5.3,e=0) ) # AMPA
    lvsy.append(recveccurr(lsy[-1]))
    lsy.append(SynapseNMDA(sect=sec,loc=loc,tau1NMDA=15,tau2NMDA=150,r=1,e=0)) # NMDA
    lvsy.append(recveccurr(lsy[-1]))
    llsy.append(lsy); llvsy.append(lvsy);
  cell.dsy[sec] = llsy; cell.vsy[sec] = llvsy;

[cyt,er,cyt_er_membrane,ca,caextrude,serca,leak,CB,caCB,buffering]=[None for i in range(10)]

vt=Vector(); vt.record(h._ref_t);
vdt=Vector(); vdt.record(h._ref_dt);

lncrec,ltimevec,lidvec=[],[],[] # spike recorders
lX,lY,lZ,lctyID=[0],[0],[0],[E2] # E2 pyramidal cell - nothing particular about the E2
timevec,idvec = h.Vector(),h.Vector()
ncrec = h.NetCon(ce[-1].soma(0.5)._ref_v, None, sec=ce[-1].soma)
ncrec.record(timevec,idvec,0)
ncrec.threshold = 0 # 10 mV is default, lower it
ltimevec.append(timevec); lidvec.append(idvec); lncrec.append(ncrec) # record the spikes at soma

rxd.options.use_reaction_contribution_to_jacobian = False # faster (checked a few days before 10/16/13)
fc, fe = 0.83, 0.17 # cytoplasmic, er volume fractions
cyt = rxd.Region(rxdsec, nrn_region='i', geometry=rxd.FractionalVolume(fc, surface_fraction=1))
er  = rxd.Region(rxdsec, geometry=rxd.FractionalVolume(fe))
cyt_er_membrane = rxd.Region(rxdsec, geometry=rxd.ScalableBorder(1))
caDiff = 0.233
ca = rxd.Species([cyt, er], d=caDiff, name='ca', charge=2, initial=dconf['cacytinit'])
caexinit = dconf['caexinit']
caextrude = rxd.Rate(ca, (caexinit-ca[cyt])/taurcada, regions=cyt, membrane_flux=False)

ip3 = rxd.Species(cyt, d=0.283, name='ip3', initial=0.0)
# action of IP3 receptor
Kip3=0.13; Kact=0.4
minf = ip3[cyt] * 1000. * ca[cyt] / (ip3[cyt] + Kip3) / (1000. * ca[cyt] + Kact)
ip3r_gate_state = rxd.State(cyt_er_membrane, initial=0.8)
h_gate = ip3r_gate_state[cyt_er_membrane]
k = dconf['gip3'] * (minf * h_gate) ** 3 
ip3r = rxd.MultiCompartmentReaction(ca[er]!=ca[cyt], k, k, membrane=cyt_er_membrane)    
# IP3 receptor gating
ip3rg = rxd.Rate(h_gate, (1. / (1 + 1000. * ca[cyt] / (0.4)) - h_gate) / 400.0)
# IP3 degradation - moves towards baseline level (ip3_init)
ip3degTau = 1000 # 1000 ms
ip3deg = rxd.Rate(ip3, (0.0-ip3[cyt])/ip3degTau, regions=cyt, membrane_flux=False)

def setmGLURflux (): # mGLUR synapses generate ip3 that is fed into rxd ip3
  for c in ce:
    if ice(c.ty): continue
    for syn,seg in zip([c.Adend3mGLUR.syn,c.Adend2mGLUR.syn,c.Adend1mGLUR.syn],[c.Adend3(0.5), c.Adend2(0.5), c.Adend1(0.5)]):
      for node in ip3.nodes(seg): 
        node.include_flux(syn._ref_rip3)

def setrecip3 ():
  for c in ce:
    if ice(c.ty): continue
    c.soma_ip3cyt = Vector(tstop/h.dt)
    c.soma_ip3cyt.record( ip3[cyt].nodes(c.soma)(0.5)[0]._ref_concentration)
    c.Adend3_ip3cyt = Vector(tstop/h.dt)
    c.Adend3_ip3cyt.record( ip3[cyt].nodes(c.Adend3)(0.5)[0]._ref_concentration)

# SERCA pump: pumps ca from cyt -> ER
Kserca = 0.1 # Michaelis constant for SERCA pump
gserca = dconf['gserca']
#Tuomo: Didn't find a way to change the reaction rate after creation of this object, so I change the gserca before the creation of this object although the rest of the parameter changes are done later on.
dconfkeys = list(dconf.keys())
for ikey in range(0,len(dconfkeys)):
  if dconfkeys[ikey] == 'gCoeff_serca':
    gserca = gserca * float(dconf[dconfkeys[ikey]])
    print("gserca multiplied by "+str(dconf[dconfkeys[ikey]]))

serca = rxd.MultiCompartmentReaction(ca[cyt]>ca[er],gserca*(1e3*ca[cyt])**2/(Kserca**2+(1e3*ca[cyt])**2),membrane=cyt_er_membrane,custom_dynamics=True)
gleak = dconf['gleak']   # leak channel: bidirectional ca flow btwn cyt <> ER
leak = rxd.MultiCompartmentReaction(ca[er]!=ca[cyt], gleak, gleak, membrane=cyt_er_membrane)

def setreccaer (): # setup recording of ca[er] for each pyramidal cell in Adend3,soma center
  for c in ce:
    if ice(c.ty): continue
    c.soma_caer = Vector(tstop/h.dt)
    c.soma_caer.record( ca[er].nodes(c.soma)(0.5)[0]._ref_concentration )
    c.Adend3_caer = Vector(tstop/h.dt)
    c.Adend3_caer.record( ca[er].nodes(c.Adend3)(0.5)[0]._ref_concentration)

CB_init = dconf["CB_init"]
CB_frate = dconf["CB_frate"]
CB_brate = dconf["CB_brate"]
CBDiff = 0.043   # um^2 / msec
CB = rxd.Species(cyt,d=CBDiff,name='CB',charge=0,initial=CB_init) # CalBindin (Anwar)
caCB = rxd.Species(cyt,d=CBDiff,name='caCB',charge=0,initial=0.0) # Calcium-CB complex
kCB = [CB_frate, CB_brate] # forward,backward buffering rates
buffering = rxd.Reaction(ca+CB != caCB, kCB[0], kCB[1], regions=cyt)

def setreccacb (): # setup recording of caCB for each pyramidal cell in Adend3,soma center
  for c in ce:
    if ice(c.ty): continue
    c.soma_caCB = Vector(tstop/h.dt)
    c.soma_caCB.record( caCB.nodes(c.soma)(0.5)[0]._ref_concentration)
    c.Adend3_caCB = Vector(tstop/h.dt)
    c.Adend3_caCB.record( caCB.nodes(c.Adend3)(0.5)[0]._ref_concentration)

setreccaer() # NB: only record from RXD variables after ALL species setup!
setreccacb() # otherwise, the pointers get messed up.
setrecip3()
setmGLURflux()


#Tuomo: find the name of soma and set the hot zone calcium channels
handle = 0
somasec = 0
hotzonetype = 1 #0:hot zone in second to last. 1:hot zone in last. 2: hot zone in two last ones
hotzonecoeff = 100.0
for sec in h.allsec():
  for seg in sec:
    for pp in seg.point_processes():
      if 'IClamp' in pp.hname():
        handle = pp
        somasec = sec

handle.delay = 0
handle.dur = 0
handle.amp = 0

midpointsec = 0
h("objref apicstim")
for my1stsec in somasec.children():
  for my2ndsec in my1stsec.children():
    if hotzonetype in [0,2]:
      #Create the hot zone                                                                                                                                                                                                                  
      print('Setting hot zone in '+my2ndsec.name())
      h(my2ndsec.name()+' gcalbar_cal = gcalbar_cal * '+str(hotzonecoeff))
      h(my2ndsec.name()+' gcatbar_cat = gcatbar_cat * '+str(hotzonecoeff))
    for my3rdsec in my2ndsec.children():
      if hotzonetype in [1,2]:
        print('Setting hot zone in '+my2ndsec.name())
        #Create the hot zone                                                                                                                                                                                     
        h(my3rdsec.name()+' gcalbar_cal = gcalbar_cal * '+str(hotzonecoeff))
        h(my3rdsec.name()+' gcatbar_cat = gcatbar_cat * '+str(hotzonecoeff))

addsyns() # add syns here since mGLUR depends on ip3 creation

# stimuli

#creates NetStims (and associated NetCon,Random) - provide 'noise' inputs
#returns next useable value of sead
def makeNoiseNetStim (cel,nsl,ncl,nrl,nrlsead,syn,w,ISI,time_limit,sead):
  ns = h.NetStim()
  ns.interval = ISI
  ns.noise = 1			
  ns.number = 2 * time_limit / ISI  # create enough spikes for extra time, in case goes over limit
  ns.start = 0 # tstart
  if type(syn) == str: nc = h.NetCon(ns,cel.__dict__[syn].syn)
  else: nc = h.NetCon(ns,syn)
  nc.delay = h.dt * 2 # 0
  nc.weight[0] = w
  rds = h.Random()
  rds.Random123(sead, 0, 0)
  ns.noiseFromRandom(rds)
  ns.start = 0.
  nsl.append(ns)
  ncl.append(nc)
  nrl.append(rds)
  nrlsead.append(sead)
  cel.infod[syn] = [ns,nc] #store pointers to NetStim,NetCon for easier manipulation

#
def makeNoiseNetStims (simdur,rdmseed):
  nsl = [] #NetStim List
  ncl = [] #NetCon List
  nrl = [] #Random List for NetStims
  nrlsead = [] #List of seeds for NetStim randoms
  for cell in ce: # go through cell types, check weights,rates of inputs
    ct = cell.ty # get cell type code  
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
      llsy = cell.dsy[sec]
      for lsy in llsy:
        for i,sy in enumerate([AM2,NM2,GA2,AM,GB2]):
          if ratex[sy] > 0. and wmatex[sy] > 0. and sy != AM: # AM now for soma AMPA
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[sy],(1e3/ratex[sy]),simdur,rdmseed); rdmseed+=1
    sec = cell.Bdend; llsy = cell.dsy[sec];
    for lsy in llsy:
      for i,sy in enumerate([AM2,NM2,GA2]):
        if ratex[sy] > 0. and wmatex[sy] > 0. and sy != AM: # AM now for soma AMPA
          makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[sy],(1e3/ratex[sy]),simdur,rdmseed); rdmseed+=1
    sec = cell.soma; llsy = cell.dsy[sec];
    for i,sy in enumerate([GA,AM,NM]):
      if ratex[sy] > 0. and wmatex[sy] > 0.:
        for lsy in llsy:
          makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[sy],(1e3/ratex[sy]),simdur,rdmseed); rdmseed+=1
  return nsl,ncl,nrl,nrlsead

nsl,ncl,nrl,nrlsead = makeNoiseNetStims(tstart+tstop,dconf['iseed'])

# setup wm-related stimuli
useWM = dconf['useWM']
rnnsl = [] # based on row number -- should always have the same number of entries as rows in the nqm 
rnncl = [] # based on row number times number of cells in the row
nqm=None; rnrds = []; rnseed = []; fihSIGns=None; 
def setupWMInputs ():
  global nqm,rnnsl,rnncl,rnrds,rnseed,fihSIGns
  if os.path.exists(dconf['nqm']): # read nqm if it exists
    nqm = h.NQS(dconf['nqm'])
  else: # setup nqm here if nqm file path non-existant
    nqm = h.NQS('vid','startt','endt','rate','w'); nqm.odec('vid')
    lvid = []; nMem = dconf['nMem']; 
    memstartt = dconf['memstartt']; memintert = dconf['memintert']; memdurt = dconf['memdurt']
    memstartt += tstart # this only has an effect when loadstate != 0
    startt = memstartt; endt = startt+memdurt
    if verbose: print('startt:',startt,',memstartt:',memstartt,',tstart:',tstart,',loadstate:',loadstate)
    for i in range(nMem): # number of stimuli
      lvid.append(h.Vector()); vtmp=Vector()
      vtmp.append(0) # only one cell
      lvid[-1].append(vtmp)
      nqm.append(lvid[-1],startt,endt,dconf['memrate'],dconf['memW'])
      startt += (memdurt+memintert);
      if i == nMem-2 and dconf['lastmemdurt'] != memdurt:
        endt = startt + dconf['lastmemdurt']
      else:
        endt = startt + memdurt;
  # backup the nqm file for later retrieval
  fnqwm = 'meminput/' + simstr + '_nqm.nqs'
  if os.path.exists(fnqwm):
    os.remove(fnqwm)
  nqm.sv(fnqwm)
  if verbose:
    print('this is nqm:')
    nqm.pr()

  def createNetStims (vid,rate,w,startt,endt,seed=1234,syn='Adend3AMPA'):
    global rnnsl,rnncl,rnrds,rnseed # based on row number in nqm
    rnnsl.append( [] ) # a list for each row of nqm
    rnncl.append( [] )
    rnrds.append( [] )
    rnseed.append( [] )
    sidx = -1
    slist=[cell.Adend3, cell.Adend2, cell.Adend1]
    if syn.count('AMPA'): sidx=0
    elif syn.count('NMDA'): sidx=1
    elif syn.count('GABAss'): sidx=2
    elif syn.count('mGLUR'): sidx=3
    if syn.count('0'): slist=[cell.Adend3, cell.Adend2, cell.Adend1]
    elif syn.count('1'): slist=[cell.Adend1]
    elif syn.count('2'): slist=[cell.Adend2]
    elif syn.count('3'): slist=[cell.Adend3]
    if sidx == -1:
      ns = h.NetStim()
      ns.start = startt
      ns.number = (endt-startt) * rate / 1e3
      if verbose: print('createNetStims:',startt,endt,ns.number,rate,w)
      ns.interval = 1e3 / rate
      ns.noise = 1
      rnnsl[-1].append(ns)
      nc = h.NetCon(ns,ce[0].__dict__[syn].syn)
      nc.delay = h.dt * 2
      nc.weight[0] = w
      rnncl[-1].append(nc)
      rds = h.Random()
      rds.Random123(seed, 0, 0)
      ns.noiseFromRandom(rds)
      rnrds[-1].append(rds)
      rnseed[-1].append(seed)
    else:
      tmpseed = seed     
      for sec in slist:
        llsy = cell.dsy[sec]
        for lsy in llsy:
          ns = h.NetStim()
          ns.start = startt
          ns.number = (endt-startt) * rate / 1e3
          if verbose: print('createNetStims:',startt,endt,ns.number,rate,w)
          ns.interval = 1e3 / rate
          ns.noise = 1
          rnnsl[-1].append(ns)
          nc = h.NetCon(ns,lsy[sidx].syn)
          nc.delay = h.dt * 2
          nc.weight[0] = w
          rnncl[-1].append(nc)
          rds = h.Random()
          rds.Random123(tmpseed, 0, 0)
          ns.noiseFromRandom(rds)
          rnrds[-1].append(rds)
          rnseed[-1].append(tmpseed)
          tmpseed += 1

  def createNetStimsFromNQ (nqm,row,seed=1234,syn='Adend3AMPA',wfctr=1.0):
    nqm.tog("DB")
    vid = h.Vector()
    rate = nqm.getcol("rate").x[row]
    w = nqm.getcol("w").x[row] * wfctr # wfctr is weight scaling factor
    vid.copy(nqm.get("vid",row).o[0])
    startt = float(nqm.getcol("startt").x[row]) # + tstart
    endt = float(nqm.getcol("endt").x[row]) # + tstart
    createNetStims(vid,rate,w,startt,endt,seed,syn)

  def setStims ():
    global nqm
    lapicIDX = []
    try:
      lapicIDX = [int(dconf['apicIDX'])]
    except:
      for i in dconf['apicIDX'].split(','): lapicIDX.append(int(i))     
    nqm.tog("DB"); sz = int(nqm.v[0].size()); lastmGLURON = dconf['lastmGLURON']; nMemInhib = dconf['nMemInhib']
    if lastmGLURON:
      for i in range(sz):
        if i == sz-1:
          for j in lapicIDX: createNetStimsFromNQ(nqm,i,seed=(i+1)*12345,syn='Adend'+str(j)+'mGLUR',wfctr=mGLURRWM)
        else:
          for j in lapicIDX:
            createNetStimsFromNQ(nqm,i,seed=(i+1)*12345,syn='Adend'+str(j)+'AMPA',wfctr=AMRWM)
            createNetStimsFromNQ(nqm,i,seed=(i+1)*12345,syn='Adend'+str(j)+'NMDA',wfctr=NMAMRWM)
    elif nMemInhib > 0:
      for i in range(0,sz-nMemInhib,1):
        for j in lapicIDX:
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'AMPA',wfctr=AMRWM)
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'NMDA',wfctr=NMAMRWM)
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'mGLUR',wfctr=mGLURRWM)
        for i in range(sz-nMemInhib,sz,1):
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'GABAss',wfctr=GARWM)
    else:
      for i in range(sz):
        for j in lapicIDX:
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'AMPA',wfctr=AMRWM)
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'NMDA',wfctr=NMAMRWM)
          createNetStimsFromNQ(nqm,i,seed=j*(i+1)*12345,syn='Adend'+str(j)+'mGLUR',wfctr=mGLURRWM)

  def initSigNetStims ():
    for i in range(len(rnnsl)):
      for j in range(len(rnnsl[i])):
        rds = rnrds[i][j]
        sead = rnseed[i][j]
        rds.MCellRan4(sead,sead)
        rds.negexp(1)			

  fihSIGns = h.FInitializeHandler(0, initSigNetStims)
  setStims() # create the inputs based on contents of nqm

#this should be called @ beginning of each sim - done in an FInitializeHandler
def init_NetStims ():
  for i in range(len(nrl)):
    rds = nrl[i]
    sead = nrlsead[i]
    rds.MCellRan4(sead,sead)
    rds.negexp(1)			

fihns = h.FInitializeHandler(0, init_NetStims)

#
def printt ():
  sys.stdout.write('\rt: {0} ...'.format(h.t))
  sys.stdout.flush()

# handler for printing out time during simulation run
def fi():
  for i in range(int(tstart),int(tstart+tstop),100): h.cvode.event(i, printt)

fih = h.FInitializeHandler(1, fi)

##### RECORDINGS

# get x locations in a line, ordered from Bdend to Adend3
def getlx ():
  lx = []; startx = 0
  for sdx in [1,0,2,3,4]:
    sec = cell.all_sec[sdx]
    LX = [node.x for node in ca[cyt].nodes(sec)]
    if sdx == 1:
      LX.reverse()
      for x in LX: lx.append(-x*sec.L + startx)
    else:
      for x in LX: lx.append(x*sec.L + startx)
    if sdx != 1: startx += sec.L
  return lx

# cell.snames # ['soma', 'Bdend', 'Adend1', 'Adend2', 'Adend3']
#
lcai = []; lcaer = []; lv = []; lip3 = []
for sdx in [1,0,2,3,4]:
  sec = cell.all_sec[sdx]
  lx = [node.x for node in ca[cyt].nodes(sec)]
  if sdx == 1: lx.reverse() # since Bdend connected at other end
  for x in lx:
    if len(ca[er].nodes(sec)(x)) > 0:
      vc = Vector(); vc.record(ca[cyt].nodes(sec)(x)[0]._ref_concentration)
      lcai.append(vc)
      vc = Vector(); vc.record(ca[er].nodes(sec)(x)[0]._ref_concentration)
      lcaer.append(vc)
      vc = Vector(); vc.record(sec(x)._ref_v)
      lv.append(vc)
      vc = Vector(); vc.record(ip3.nodes(sec)(x)[0]._ref_concentration)
      lip3.append(vc)

# plot currents
def plotcurr (dstr='Adend3',dsec=cell.Adend3):
  figure(); subplot(1,3,1)
  tv = arange(0,tstop,recvdt)
  llvsy = cell.vsy[dsec]
  idx = int(len(llvsy) / 2)
  clr = ['r','g','b','k']
  for i in range(4): plot(tv,llvsy[idx][i],clr[i])
  legend(['AM','NM','GA','GB'],loc='best')
  xlabel('Time (ms)'); ylabel('i')  
  subplot(1,3,2)
  li = [cell.__dict__[dstr + s] for s in ['_ik','_ina','_ica','_ih']]
  for vec,clr in zip(li,['b','g','r', 'c']): plot(tv,vec.as_numpy(),clr)
  legend(['ik','ina','ica','ih'],loc='best'); xlabel('Time (ms)'); ylabel(dstr + ' i')
  subplot(1,3,3)
  for vec,clr in zip([cell.soma_ik,cell.soma_ina,cell.soma_ica,cell.soma_ih],['b','g','r', 'c']): plot(tv,vec.as_numpy(),clr)
  legend(['ik','ina','ica','ih'],loc='best'); xlabel('Time (ms)'); ylabel('soma i')

# get sag - for hyperpolarizing pulses - difference between minimum and steadystate
def getsag (dat,sampr,stimon=0.5,verbose=False):
  vpeak = amin(dat[int(stimon*sampr):int((stimon+0.1)*sampr)])
  vsteady = mean(dat[int((stimon+0.89)*sampr):int((stimon+0.99)*sampr)])
  sag = (vpeak-vsteady)/vsteady
  if verbose: print('sag:',sag,'peak:',vpeak, 'steady:',vsteady)
  return sag

#
def getrin (dat,dt,ic):
  sidx,eidx = int(ic.delay/dt),int((ic.delay+ic.dur)/dt)
  if ic.amp < 0: dv = dat[sidx] - amin(dat[sidx+1:eidx])
  else: dv = amax(dat[sidx+1:eidx]) - dat[sidx]
  di = abs(ic.amp)
  print('dv:',dv,'di:',di,'rin:',dv/di)
  return dv/di

#
def plotout (dstr='Adend3',xl=((0,tstop/1e3))):
  lx = getlx(); cellL = [lx[0],lx[-1]]
  tv = arange(0,tstop,h.dt)
  figure(); npa=numpy.array(lcai);
  subplot(2,2,1); imshow(npa,origin='lower',aspect='auto',extent=(0,tstop/1e3,lx[0],lx[-1]))
  title(r'$Ca^{2+}_{cyt}$'); colorbar(); xlim(xl)
  subplot(2,2,2); npb=numpy.array(lcaer);
  imshow(npb,origin='lower',aspect='auto',extent=(0,tstop/1e3,lx[0],lx[-1]))
  title(r'$Ca^{2+}_{ER}$'); colorbar(); xlim(xl)
  subplot(2,2,3); npv=numpy.array(lv);
  imshow(npv,origin='lower',aspect='auto',extent=(0,tstop/1e3,lx[0],lx[-1]))
  title('v (mV)'); colorbar(); xlim(xl)
  subplot(2,2,4); npi=numpy.array(lip3);
  imshow(npi,origin='lower',aspect='auto',extent=(0,tstop/1e3,lx[0],lx[-1]))
  title(r'$IP_3 (mM)$'); colorbar(); xlim(xl)
  tvec = vt.as_numpy()
  if len(tv) == cell.soma_volt.size(): tvec = tv
  figure(); 
  subplot(1,4,1); plot(tvec/1e3,cell.soma_cai,'b'); plot(tvec/1e3,cell.__dict__[dstr+'_cai'],'r');title(r'$Ca^{2+}_{cyt}$');
  legend(['soma','dend'],loc='best'); xlim(xl)
  subplot(1,4,2);plot(tvec/1e3,cell.soma_Ihp1,'b');plot(tvec/1e3,cell.__dict__[dstr+'_Ihp1'],'r');title(r'$cAMP$');legend(['soma','dend'],loc='best'); xlim(xl)
  subplot(1,4,3); plot(tvec/1e3,cell.soma_Ihm,'b'); plot(tvec/1e3,cell.__dict__[dstr+'_Ihm'],'r');title(r'$g_h$');legend(['soma','dend'],loc='best'); xlim(xl)
  subplot(1,4,4);
  plot(tvec/1e3,cell.soma_volt,'b'); plot(tvec/1e3,cell.__dict__[dstr+'_volt'],'r'); title('v'); legend(['soma','dend'],loc='best');  xlim(xl)

def plotfixeddt (vec,clr,lw): plot(linspace(0,tstop/1e3,vec.size()),vec,clr,linewidth=lw)

#
def plotnew (dstr='Adend3',dsec=cell.Adend3,xl=((0,tstop/1e3)),useyl=True):
  naxbins=5; tx,ty=-.07,1.04; fsz=20; lw=2
  global vt
  vtime = vt.as_numpy() / 1e3
  figure(); 
  # synaptic current i
  ax=subplot(2,4,1); ax.locator_params(nbins=naxbins); text(tx,ty,'a',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  llvsy = cell.vsy[dsec]
  idx = int(len(llvsy) / 2)
  clr = ['r','g','b']
  if not dconf['cvodeactive']:
    for i in range(3): plotfixeddt(llvsy[idx][i],clr[i],lw)
  else:
    for i in range(3): plot(vtime,llvsy[idx][i],clr[i],linewidth=lw)
  #legend([r'$i_{AMPA}$',r'$i_{NMDA}$',r'$i_{GABAA}$'],loc='best')
  #xlabel('Time (s)'); 
  title(r'$syn_i$ (nA)',fontsize=fsz); xlim(xl); 
  if useyl: ylim((-1.22,.06)); #ylim((-.64,.06))
  # cai
  ax=subplot(2,4,2); ax.locator_params(nbins=naxbins);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_cai,'b',lw)
    plotfixeddt(cell.__dict__[dstr+'_cai'],'r',lw);
  else:
    plot(vtime,cell.soma_cai,'b',linewidth=lw); 
    plot(vtime,cell.__dict__[dstr+'_cai'],'r',linewidth=lw);
  title(r'$Ca^{2+}_{cyt}$ (mM)',fontsize=fsz); xlim(xl); 
  if useyl: ylim((0,0.0094))  
  #legend(['soma','dend'],loc='best'); 
  text(tx,ty,'b',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  # ip3
  ax=subplot(2,4,3); ax.locator_params(nbins=naxbins);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_ip3cyt,'b',lw);
    plotfixeddt(cell.__dict__[dstr+'_ip3cyt'],'r',lw);
  else:
    plot(vtime,cell.soma_ip3cyt,'b',linewidth=lw);
    plot(vtime,cell.__dict__[dstr+'_ip3cyt'],'r',linewidth=lw);
  title(r'$IP_3$ (mM)',fontsize=fsz); xlim(xl); 
  if useyl: ylim((0,0.214))  
  text(tx,ty,'c',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  # caer
  ax=subplot(2,4,4); ax.locator_params(nbins=naxbins);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_caer,'b',lw); 
    plotfixeddt(cell.__dict__[dstr+'_caer'],'r',lw);
  else:
    plot(vtime,cell.soma_caer,'b',linewidth=lw); 
    plot(vtime,cell.__dict__[dstr+'_caer'],'r',linewidth=lw);
  title(r'$Ca^{2+}_{ER}$ (mM)',fontsize=fsz); xlim(xl); 
  if useyl: ylim((0,1.3))
  text(tx,ty,'d',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  # cAMP
  ax=subplot(2,4,5); ax.locator_params(nbins=naxbins);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_Ihp1,'b',lw); 
    plotfixeddt(cell.__dict__[dstr+'_Ihp1'],'r',lw);
  else:
    plot(vtime,cell.soma_Ihp1,'b',linewidth=lw); 
    plot(vtime,cell.__dict__[dstr+'_Ihp1'],'r',linewidth=lw);
  title(r'$cAMP$',fontsize=fsz); xlim(xl); 
  if useyl: ylim((0,0.025))
  text(tx,ty,'e',fontsize=fsz,fontweight='bold',transform=ax.transAxes); xlabel('Time (s)',fontsize=fsz);
  # gh
  ax=subplot(2,4,6); ax.locator_params(nbins=naxbins);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_Ihm,'b',lw); 
    plotfixeddt(cell.__dict__[dstr+'_Ihm'],'r',lw);
  else:
    plot(vtime,cell.soma_Ihm,'b',linewidth=lw); 
    plot(vtime,cell.__dict__[dstr+'_Ihm'],'r',linewidth=lw);
  title(r'$g_h$',fontsize=fsz); xlim(xl); 
  if useyl: ylim((0.011,0.168))
  text(tx,ty,'f',fontsize=fsz,fontweight='bold',transform=ax.transAxes); xlabel('Time (s)',fontsize=fsz);
  # dend ions
  ax=subplot(2,4,7); ax.locator_params(nbins=naxbins); 
  text(tx,ty,'g',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  #li = [cell.__dict__[dstr + s] for s in ['_ik','_ina','_ica','_ih']]
  #for vec,clr in zip(li,['b','g','r','c']): plot(tv/1e3,vec.as_numpy(),clr)
  li = [cell.__dict__[dstr + s] for s in ['_ih']]
  if not dconf['cvodeactive']:
    for vec,clr in zip(li,['r']): plotfixeddt(vec,clr,lw)
  else:
    for vec,clr in zip(li,['r']): plot(vtime,vec,clr,linewidth=lw)
  #legend([r'$i_K$',r'$i_{Na}$',r'$i_{Ca}$',r'$i_h$'],loc='best'); 
  xlabel('Time (s)',fontsize=fsz); title(r'$i_h$ (nA)',fontsize=fsz); xlim(xl); 
  if useyl: ylim((-.102,-.022))
  if False:
    # soma ions
    ax=subplot(2,4,7); ax.locator_params(nbins=naxbins); 
    text(tx,ty,'h',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
    for vec,clr in zip([cell.soma_ik,cell.soma_ina],['b','r']): plot(vtime,vec,clr,linewidth=lw)
    #legend([r'$i_K$',r'$i_{Na}$'],loc='best'); 
    xlabel('Time (s)',fontsize=fsz); title(r'$soma_i$ (nA)',fontsize=fsz); xlim(xl); 
    if useyl: ylim((-15.5,11.5))
  # voltage
  ax=subplot(2,4,8); ax.locator_params(nbins=naxbins); text(tx,ty,'h',fontsize=fsz,fontweight='bold',transform=ax.transAxes);
  if not dconf['cvodeactive']:
    plotfixeddt(cell.soma_volt,'b',lw); 
    plotfixeddt(cell.__dict__[dstr+'_volt'],'r',lw); 
  else:
    plot(vtime,cell.soma_volt,'b',linewidth=lw); 
    plot(vtime,cell.__dict__[dstr+'_volt'],'r',linewidth=lw); 
  title('Vm',fontsize=fsz); xlim(xl)
  xlabel('Time (s)',fontsize=fsz); 
  if useyl: ylim((-85,40))

if dconf['cvodeactive']: h.cvode.active(1) # turn on variable time-step integration

lspks,lids=array([]),array([])

# interpolate voltage recorded in simulation to a fixed grid (dt millisecond spacing)
# output time,voltage returned
def interpvolt (tsrc,vsrc,dt):
  tdest = h.Vector(); tdest.indgen(0,tstop,dt)
  vdest = h.Vector(); vdest.interpolate(tdest,tsrc,vsrc)
  return tdest,vdest

mechanisms = ['pas','iar','nax','kdr','kap','cal','can','cat','ikc','km','cagk','kdmc','im','calts','ihlts','kcalts','Kdrbwb','Nafbwb']
for ikey in range(0,len(dconfkeys)):
  if dconfkeys[ikey][0:7] == 'gCoeff_' and dconfkeys[ikey] not in ['gCoeff_serca', 'gCoeff_gbar_GABAB']:
    imech_chosen = -1
    for imech in range(0,len(mechanisms)):
      if mechanisms[imech] in dconfkeys[ikey]:
        if imech_chosen > -1 and mechanisms[imech] in mechanisms[imech_chosen]: #Don't use a shorter mechanism if a longer one has been assigned already
          continue
        imech_chosen = imech
    if imech_chosen == -1:
      print('Mechanism for '+dconfkeys[ikey]+' not found!')
      quit()
    print('forall if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][7:]+' = '+dconf[dconfkeys[ikey]]+'*'+dconfkeys[ikey][7:])
    h('forall if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][7:]+' = '+dconf[dconfkeys[ikey]]+'*'+dconfkeys[ikey][7:])
  elif dconfkeys[ikey][0:11] == 'gCoeffApic_' and dconfkeys[ikey] not in ['gCoeffApic_serca', 'gCoeffApic_gbar_GABAB']:
    imech_chosen = -1
    for imech in range(0,len(mechanisms)):
      if mechanisms[imech] in dconfkeys[ikey]:
        if imech_chosen > -1 and mechanisms[imech] in mechanisms[imech_chosen]: #Don't use a shorter mechanism if a longer one has been assigned already
          continue
        imech_chosen = imech
    if imech_chosen == -1:
      print('Mechanism for '+dconfkeys[ikey]+' not found!')
      quit()
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
      print(sec.name()+' if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][11:]+' = '+dconf[dconfkeys[ikey]]+'*'+dconfkeys[ikey][11:])
      h(sec.name()+' if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][11:]+' = '+dconf[dconfkeys[ikey]]+'*'+dconfkeys[ikey][11:])
  elif dconfkeys[ikey][0:7] == 'paraCh_':
    imech_chosen = -1
    for imech in range(0,len(mechanisms)):
      if mechanisms[imech] in dconfkeys[ikey]:
        if imech_chosen > -1 and mechanisms[imech] in mechanisms[imech_chosen]: #Don't use a shorter mechanism if a longer one has been assigned already
          continue
        imech_chosen = imech
    if imech_chosen == -1:
      print('Mechanism for '+dconfkeys[ikey]+' not found!')
      quit()
    print('forall if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][7:]+' = '+dconf[dconfkeys[ikey]])
    h('forall if(ismembrane("'+mechanisms[imech_chosen]+'")) '+dconfkeys[ikey][7:]+' = '+dconf[dconfkeys[ikey]])
  elif dconfkeys[ikey] == 'gCoeff_gbar_GABAB':
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1, cell.Bdend]:
      mydsy = cell.dsy[sec]
      for isyn1 in range(0,len(mydsy)):
         for isyn2 in range(0,len(mydsy[isyn1])):
           if 'GABAB' in mydsy[isyn1][isyn2].syn.hname():
             print(mydsy[isyn1][isyn2].syn.hname() + " gbar multiplied by "+str(dconf[dconfkeys[ikey]]))
             mydsy[isyn1][isyn2].syn.gbar = mydsy[isyn1][isyn2].syn.gbar * float(dconf[dconfkeys[ikey]])
             print(mydsy[isyn1][isyn2].syn.hname() + " gbar now "+str(mydsy[isyn1][isyn2].syn.gbar))
  elif dconfkeys[ikey] == 'Cm_apic':
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
      sec.cm = float(dconf[dconfkeys[ikey]])
print("dconfkeys done")    

#
def myrun ():
  global fit # pointer to FinitializeHandler needs to be global so doesn't disappear
  dastart = datetime.now()
  print('started at:',dastart)
  if dconf['useWM']: setupWMInputs()
  h.stdinit()
  ca[er].concentration = dconf['caerinit'] # 100e-6
  ca[cyt].concentration = dconf['cacytinit'] # 100e-6
  if dconf['cvodeactive']: h.cvode.re_init()
  h.continuerun(h.t+tstop) # run for tstop
  daend = datetime.now()
  print('finished ',tstop,' ms sim at:',daend)
  dadiff = daend - dastart;
  print('runtime:',dadiff, '(',tstop/1e3,' s)')
  # concatenate the results so can view/save all at once
  global lspks,lids
  lspks,lids=array([]),array([])
  for i in range(len(ltimevec)):
    lspks=concatenate((lspks,ltimevec[i]))
    lids=concatenate((lids,lidvec[i]))    
  ion()

myrun()

# save data - output path based on simstr 
def savedata (simstr):
  safemkdir(outdir) # make sure the output directory exists
  fn = outdir + '/' + simstr + '.npz'
  print(' saving to ' , fn)
  ne,ni,szslow,szfast = 0,0,0,0
  lE,lI=[],[]
  for c in ce:
    if ice(c.ty):
      lI.append(c.ID)
      ni += 1
    else:
      szslow = int(c.soma_cai.size()) # only E cells
      lE.append(c.ID)
      ne += 1
    szfast = int(c.soma_volt.size())
  lE=array(lE) # lE is list of E cell IDs from this host
  lI=array(lI) # Li is list of I cell IDs from this host
  soma_volt = zeros((ne,szfast)); Adend3_volt = zeros((ne,szfast)); Bdend_volt=zeros((ne,szfast));
  soma_cai = zeros((ne,szslow)); Adend3_cai = zeros((ne,szslow))
  soma_Ihm = zeros((ne,szslow)); Adend3_Ihm = zeros((ne,szslow))
  soma_Ihp1 = zeros((ne,szslow)); Adend3_Ihp1 = zeros((ne,szslow))
  soma_voltI = zeros((ni,szfast));
  soma_caer = zeros((ne,szslow));  Adend3_caer = zeros((ne,szslow));
  soma_caCB=zeros((ne,szslow)); Adend3_caCB=zeros((ne,szslow));
  soma_ip3cyt=zeros((ne,szslow)); Adend3_ip3cyt=zeros((ne,szslow));
  saveExtra = dconf['saveExtra']
  if saveExtra:
    Adend3_iAM = zeros((ne,szfast)); Adend3_iNM = zeros((ne,szfast))
    Adend3_iGB = zeros((ne,szfast)); Adend3_iGA = zeros((ne,szfast))
    soma_iGA = zeros((ne,szfast))
    [Adend3_ina, Adend3_ik, Adend3_ica, Adend3_ih] = [zeros((ne,szfast)) for i in range(4)]
    [soma_ina, soma_ik, soma_ica, soma_ih] = [zeros((ne,szfast)) for i in range(4)]
  cdx = 0; idx = 0;
  for c in ce:
    if ice(c.ty):
      soma_voltI[idx,:] = c.soma_volt.to_python()
      idx += 1
      continue
    soma_volt[cdx,:] = c.soma_volt.to_python()
    Adend3_volt[cdx,:] = c.Adend3_volt.to_python()
    Bdend_volt[cdx,:] = c.Bdend_volt.to_python()
    soma_cai[cdx,:] = c.soma_cai.to_python()
    Adend3_cai[cdx,:] = c.Adend3_cai.to_python()
    soma_Ihm[cdx,:] = c.soma_Ihm.to_python()
    Adend3_Ihm[cdx,:] = c.Adend3_Ihm.to_python()
    soma_Ihp1[cdx,:] = c.soma_Ihp1.to_python()
    Adend3_Ihp1[cdx,:] = c.Adend3_Ihp1.to_python()
    soma_caer[cdx,:] = c.soma_caer.to_python()
    Adend3_caer[cdx,:] = c.Adend3_caer.to_python()
    soma_caCB[cdx,:] = c.soma_caCB.to_python()
    Adend3_caCB[cdx,:] = c.Adend3_caCB.to_python()
    soma_ip3cyt[cdx,:] = c.soma_ip3cyt.to_python()
    Adend3_ip3cyt[cdx,:] = c.Adend3_ip3cyt.to_python()        
    if saveExtra:
      Adend3_iAM[cdx,:] = c.Adend3_iAM.to_python()
      Adend3_iNM[cdx,:] = c.Adend3_iNM.to_python()
      Adend3_iGB[cdx,:] = c.Adend3_iGB.to_python()
      Adend3_iGA[cdx,:] = c.Adend3_iGA.to_python()
      soma_iGA[cdx,:] = c.soma_iGA.to_python()
      Adend3_ina[cdx,:] = c.Adend3_ina.to_python()
      Adend3_ik[cdx,:] = c.Adend3_ik.to_python()
      Adend3_ica[cdx,:] = c.Adend3_ica.to_python()
      Adend3_ih[cdx,:] = c.Adend3_ih.to_python()
      soma_ina[cdx,:] = c.soma_ina.to_python()
      soma_ik[cdx,:] = c.soma_ik.to_python()
      soma_ica[cdx,:] = c.soma_ica.to_python()
      soma_ih[cdx,:] = c.soma_ih.to_python()    
    cdx += 1
  if saveExtra: # with somaVolt, synaptic currents
    numpy.savez(fn,lctyID=array(lctyID),lX=array(lX),lY=array(lY),lZ=array(lZ),vt=numpy.array(vt.to_python()),lspks=lspks,lids=lids,lE=lE,lI=lI,soma_volt=soma_volt,soma_voltI=soma_voltI,Adend3_volt=Adend3_volt,Bdend_volt=Bdend_volt,soma_cai=soma_cai,Adend3_cai=Adend3_cai,soma_Ihm=soma_Ihm,Adend3_Ihm=Adend3_Ihm,soma_Ihp1=soma_Ihp1,Adend3_Ihp1=Adend3_Ihp1,soma_caer=soma_caer,Adend3_caer=Adend3_caer,soma_caCB=soma_caCB,Adend3_caCB=Adend3_caCB,Adend3_ip3cyt=Adend3_ip3cyt,soma_ip3cyt=soma_ip3cyt,Adend3_iAM=Adend3_iAM,Adend3_iNM=Adend3_iNM,Adend3_iGB=Adend3_iGB,Adend3_iGA=Adend3_iGA,soma_iGA=soma_iGA,Adend3_ina=Adend3_ina,Adend3_ik=Adend3_ik,Adend3_ica=Adend3_ica,Adend3_ih=Adend3_ih,soma_ina=soma_ina,soma_ik=soma_ik,soma_ica=soma_ica,soma_ih=soma_ih)
  else: # no somaVolt, synaptic currents
    numpy.savez(fn,lctyID=array(lctyID),lX=array(lX),lY=array(lY),lZ=array(lZ),vt=numpy.array(vt.to_python()),lspks=lspks,lids=lids,lE=lE,lI=lI,Adend3_volt=Adend3_volt,Bdend_volt=Bdend_volt,soma_cai=soma_cai,Adend3_cai=Adend3_cai,soma_Ihm=soma_Ihm,Adend3_Ihm=Adend3_Ihm,soma_Ihp1=soma_Ihp1,Adend3_Ihp1=Adend3_Ihp1,soma_caer=soma_caer,Adend3_caer=Adend3_caer,soma_caCB=soma_caCB,Adend3_caCB=Adend3_caCB,Adend3_ip3cyt=Adend3_ip3cyt,soma_ip3cyt=soma_ip3cyt)

if saveout:
  print('saving data')
  savedata(simstr)
if saveoutSpikes:
  spikes = mytools.spike_times([recvdt*i for i in range(0,len(array(cell.soma_volt)))],array(cell.soma_volt))
  scipy.io.savemat(simstr+'.mat',{'spikes': spikes, 'Vsoma': array(cell.soma_volt), 'times': [recvdt*i for i in range(0,len(array(cell.soma_volt)))]})

# get an array of times of interest
def gettimes (nqm,rampt=3e3):
  nqm.tog("DB")
  nMem = int(nqm.size())
  baselinest = dconf['baset']
  times = [] # periods of interest
  vstartt,vendt = Vector(),Vector()
  vstartt.copy(nqm.getcol("startt"))
  vendt.copy(nqm.getcol("endt"))
  #vstartt.add(tstart); vendt.add(tstart)
  if nMem > 0:
    times = [(baselinest,vstartt[0])]
    times.append((vstartt[0],vendt[0]))
    times.append((vendt[0],vendt[0]+rampt))
    for row in range(1,nMem):
      times.append((vendt[row-1]+rampt,vstartt[row]))
      times.append((vstartt[row],vendt[row]))
      times.append((vendt[row],vendt[row]+rampt))
    lastt = vendt[-1]+rampt
    inct = 1e3; wint = 1e3
    while lastt < tstop:
      times.append((lastt,lastt+wint))
      lastt += inct
  return times

# get an NQS with spikes (snq)
def getsnq ():
  snq = NQS('id','t','ty','ice')
  snq.v[0].from_python(lids)
  snq.v[1].from_python(lspks)
  for i in range(len(lids)):
    snq.v[2].append(lctyID[int(lids[i])])
    snq.v[3].append(h.ice(lctyID[int(lids[i])]))
  return snq

# return x with only two decimal places after x
def twop (x): return round(x,2)

# print spike rates in each population in the periods of interest (baseline, signal, recall)
def prspkcount (snq,nqm,rampt=1e3,times=None,quiet=False):
  snq.verbose=0
  if times is None: times = gettimes(nqm,rampt) # periods of interest
  lfa=[]
  nqm.tog("DB")
  nMem = int(nqm.size())
  for row in range(nMem):
    lfa.append([]); 
    vact = Vector(1); vact.x[0] = 0
    for (startt,endt) in times: 
      fa = twop(1e3*snq.select("id","EQW",vact,"t","[]",startt,endt) / ((endt-startt)*vact.size()))
      lfa[-1].append(fa)
      if not quiet:
        print("t=",startt,"-",endt,". ActE:",fa,"Hz.")
  snq.verbose=1
  return lfa

# print spike rates in each population in the periods of interest (baseline, signal, recall)
def getspkcount (snq,binsz,quiet=False):
  snq.verbose=0
  times = arange(0,h.tstop,binsz)
  lfa=[]
  for startt in times:
    lfa.append([]); 
    vact = Vector(1); vact.x[0] = 0
    fa = twop(1e3*snq.select("id","EQW",vact,"t","[]",startt,startt+binsz) / (binsz*vact.size()))
    lfa[-1].append(fa)
    if not quiet: print("t=",startt,"-",startt+binsz,". ActE:",fa,"Hz.")
  snq.verbose=1
  return times,lfa

#
def getonecellfrd (lfa,binsz,t1,t2):
  a1 = mean(lfa[int(t1[0]/binsz):int(t1[1]/binsz)])
  a2 = mean(lfa[int(t2[0]/binsz):int(t2[1]/binsz)])
  if a1 > 0.: frd = (a2 - a1) / a1
  else: frd = 0.
  return a1,a2,frd

# print spike counts & draw activated and non-activated rates
def drawrat (ld,incsz=500.0,winsz=1e3,quiet=False):
  snq=getsnq(ld);
  times=[]
  for i in range(2*int(tstop/winsz)):times.append( (i*incsz,(i*incsz)+winsz) )
  lfa,lfn = prspkcount(snq,nqm,rampt=1e3,times=times,quiet=quiet)
  rat = []; tt=[];
  for i in range(len(lfa[0])): rat.append(lfa[0][i]/lfn[0][i])
  for tup in times: tt.append(0.5*(tup[0]+tup[1])) # time midpoints
  plot(tt,rat,'b',linewidth=2);plot(tt,rat,'bo',markersize=10);tight_layout();show() #

h("forall psection()")
if type(dconf['doquit']) is str and int(dconf['doquit']):
  print('Quitting, doquit = '+str(dconf['doquit']))
  quit()

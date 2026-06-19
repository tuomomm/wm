from mpi4py import MPI
from neuron import *
import sys
import os
import string
import time

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

## Set up MPI: see Hines and Carnevale 2008
# MPI
pc = h.ParallelContext() # MPI: Initialize the ParallelContext class
nhosts = int(pc.nhost()) # Find number of hosts
if pc.id()==0: 
  print('\nSetting up network...')
  safemkdir('data')
  safemkdir('meminput')
# MPI

h("strdef simname, allfiles, simfiles, output_file, datestr, uname, comment")
h.simname=simname = "CaHDemo"
h.allfiles=allfiles = "geom.hoc pyinit.py mpisim.py"
h.simfiles=simfiles = "pyinit.py geom_withSK3_saveExtra.py mpisim.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "15dec31"
h.output_file=output_file = "data/15dec31.1"
h.uname=uname = "x86_64"
h("templates_loaded=0")
templates_loaded=0
h("xwindows=1.0")

h.xopen("nrnoc.hoc")
h.xopen("init.hoc")
h("proc setMemb () { }") # so e_pas will not get modified

CTYPi = 50.0
STYPi = 12.0
from pyinit import *
from labels import *

delm = numpy.zeros( (CTYPi, CTYPi) )
deld = numpy.zeros( (CTYPi, CTYPi) )
pmat = numpy.zeros( (CTYPi, CTYPi) )
synloc = numpy.zeros( (CTYPi, CTYPi) )

from geom_withSK3_saveExtra import *
from nqs import *
import random
from pylab import *
from datetime import datetime

iargv_start = 1
if sum([sys.argv[i] == '-mpi' for i in range(0,len(sys.argv))]) > 0:
  iargv_start = [i for i in range(0,len(sys.argv)) if sys.argv[i] == '-mpi'][0]+2
if len(sys.argv) > iargv_start:
  for j in range(iargv_start,len(sys.argv)):
    if sys.argv[j].endswith(".cfg"):
      continue
    print("""dconf['"""+sys.argv[j].split('=')[0]+"""'] = """+sys.argv[j].split('=')[1])
    dconf[sys.argv[j].split('=')[0]] = sys.argv[j].split('=')[1]

for k in ['tstop','iseed','pseed','wseed','dt']:
  if type(dconf[k]) == str:
    dconf[k] = float(dconf[k])
    
#########################################################################
# global params
verbose = dconf['verbose']
ISEED = dconf['iseed']
WSEED = dconf['wseed']
PSEED = dconf['pseed']
scale = dconf['scale']
gGID = 0 # global ID for cells
pmatscale = 1.0 # 1.0 / scale
spiketh = 0 # spike threshold, 10 mV is NetCon default, lower it for FS cells
simstr = dconf['simstr']
saveout = dconf['saveout']
recdt = dconf['recdt']
recvdt = dconf['recvdt']
indir = dconf['indir']
outdir = dconf['outdir']
dconf['saveExtra'] = 1

# backup the config file
def backupcfg (simstr):
  safemkdir('backupcfg')
  fout = 'backupcfg/' + simstr + '.cfg'
  if os.path.exists(fout): os.system('rm ' + fout)  
  try:
    os.system('cp ' + fcfg + ' ' + fout) # fcfg created in geom.py via conf.py
  except:
    pass
  
if pc.id()==0: backupcfg(simstr) # backup the config file

h.tstop = tstop = dconf['tstop']
tstart = 0 # only use previous end time if loading state
h.dt = dconf['dt']
h.steps_per_ms = 1/h.dt
h.v_init = -65
h.celsius = 37
h.fracca_MyExp2SynNMDABB = dconf['nmfracca'] # fraction of NMDA current that is from calcium
rdmsec = dconf['rdmsec']

EEGain = dconf['EEGain']
EIGainFS = dconf['EIGainFS']
EIGainLTS = dconf['EIGainLTS']
IEGain = dconf['IEGain']
IIGain = dconf['IIGain']
IIGainLTSFS =  IIGain / 10.0
IIGainFSLTS =  IIGain 
IIGainLTSLTS = IIGain 
IIGainFSFS =   IIGain * 50
GB2R = dconf['GB2R']; 
NMAMREE = dconf['NMAMREE'] # default of 0.1
NMAMREI = dconf['NMAMREI'] # default of 0.1
mGLURR = dconf['mGLURR'] # ratio of mGLUR weights to AM2 weights
mGLURRWM = dconf['mGLURRWM']; NMAMRWM = dconf['NMAMRWM']; AMRWM = dconf['AMRWM'] # WM related weight ratios
if dconf['nsubseg'] > 0: mGLURRWM /= float(dconf['nsubseg'])
cpernet = []  # cells of a given type for network
wmat = numpy.zeros( (CTYPi, CTYPi, STYPi) ) # internal weights
wmatex = numpy.zeros( (CTYPi, STYPi) ) # external weights
ratex = numpy.zeros( (CTYPi, STYPi) )  # external rates
EXGain = dconf['EXGain']
sgrhzEE = dconf['sgrhzEE'] # external E inputs to E cells
sgrhzEI = dconf['sgrhzEI'] # external E inputs to I cells
sgrhzIE = dconf['sgrhzIE'] # external I inputs to E cells
sgrhzII = dconf['sgrhzII'] # external I inputs to I cells
sgrhzNME = dconf['sgrhzNME'] # external NM inputs to E cells
sgrhzNMI = dconf['sgrhzNMI'] # external NM inputs to I cells
sgrhzMGLURE = dconf['sgrhzMGLURE'] # external mGLUR inputs to E cells
sgrhzGB2 = dconf['sgrhzGB2'] # external inputs onto E cell GB2 synapses

SWIRE = 1; NOWIRE = 2; wirety = dconf['wirety']
if wirety == "swire": wirety = SWIRE    # spatial wiring
else: wirety = NOWIRE # wiring off

# params for swire
colside = 45 # *sqrt(scale)
slambda = 15
axonalvelocity = 10000 # axonal velocity in um/ms -- this is 10 mm/s
#########################################################################

# setwmatex - set weights of external inputs to cells
def setwmatex ():
  for ct in range(CTYPi):
    for sy in range(STYPi):
      ratex[ct][sy]=0
      wmatex[ct][sy]=0
  for ct in range(CTYPi):
    if cpernet[ct] <= 0: continue
    if IsLTS(ct): # dendrite-targeting interneurons (LTS cells)
      ratex[ct][AM2]=sgrhzEI
      ratex[ct][NM2] = sgrhzNMI
      ratex[ct][GA]=sgrhzII
      ratex[ct][GA2]=sgrhzII
      wmatex[ct][AM2] = 0.02e-3 
      wmatex[ct][NM2] = 0.02e-3
      wmatex[ct][GA]=  0
      wmatex[ct][GA2]= 0.2e-3 
    elif ice(ct): # soma-targeting interneurons (basket/FS cells)
      ratex[ct][AM2]=sgrhzEI
      ratex[ct][NM2] = sgrhzNMI
      ratex[ct][GA]=sgrhzII
      ratex[ct][GA2]=sgrhzII
      wmatex[ct][AM2] = 0.02e-3 * 5.0
      wmatex[ct][NM2] = 0.02e-3 * 5.0
      wmatex[ct][GA]= 0
      wmatex[ct][GA2]= 0.2e-3 
    else: # E cells
      ratex[ct][AM]=sgrhzMGLURE 
      ratex[ct][AM2]=sgrhzEE
      ratex[ct][NM2]=sgrhzNME
      ratex[ct][GA]=sgrhzIE
      ratex[ct][GA2]=sgrhzIE
      ratex[ct][GB2]=sgrhzGB2
      wmatex[ct][AM2] = 0.02e-3
      wmatex[ct][NM2] = 0.02e-3
      wmatex[ct][GA] = 0.2e-3
      wmatex[ct][GA2] = 0.2e-3 
      wmatex[ct][GB2] = 5e-3
    for sy in range(STYPi): wmatex[ct][sy] *= EXGain # apply gain control

# set number of cells of a type in the network at scale==1
def setcpernet ():
  global cpernet
  cpernet = []
  for i in range(CTYPi): cpernet.append(0)
  cpernet[E2]  = 150
  cpernet[I2]  =  25
  cpernet[I2L] =  25
  cpernet[E5R] =  95
  cpernet[E5B] =  145
  cpernet[I5]  =  40
  cpernet[I5L] =  40
  cpernet[E6] =   192
  cpernet[I6] =   32
  cpernet[I6L] =  32

# synapse locations DEND SOMA AXON
def setsynloc ():
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      if ice(ty1):
        if IsLTS(ty1):
          synloc[ty1][ty2]=DEND # distal [GA2] - from LTS
        else:
          synloc[ty1][ty2]=SOMA # proximal [GA] - from FS
      else:
        synloc[ty1][ty2]=DEND # E always distal. use AM2,NM2

# setdelmats -- setup delm,deld
def setdelmats ():
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      if synloc[ty1][ty2]==DEND and ice(ty2):
        # longer delays at dendrites of interneurons since they are currently single compartment
        delm[ty1][ty2]=4
        deld[ty1][ty2]=1
      else:
        delm[ty1][ty2]=2.0
        deld[ty1][ty2]=0.2

# weight params
def setwmat ():
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      for sy in range(STYPi):
        wmat[ty1][ty2][sy]=0
  wmat[E2][E2][AM2]=0.66
  wmat[E2][E5B][AM2]=0.36
  wmat[E2][E5R][AM2]=0.93
  wmat[E2][I5L][AM2]=0.36 
  wmat[E2][E6][AM2]=0
  wmat[E2][I2L][AM2]=0.23
  wmat[E2][I2][AM2] = 0.23
  wmat[E5B][E2][AM2]=0.26
  wmat[E5B][E5B][AM2]=0.66
  wmat[E5B][E5R][AM2]=0   # pulled from Fig. 6D, Table 1 of (Kiritani et al., 2012)
  wmat[E5B][E6][AM2]=0.66
  wmat[E5B][I5L][AM2]=0   # ruled out by (Apicella et al., 2012) Fig. 7
  wmat[E5B][I5][AM2]=0.23  # (Apicella et al., 2012) Fig. 7F (weight = 1/2 x weight for E5R->I5)
  wmat[E5R][E2][AM2]=0.66
  wmat[E5R][E5B][AM2]=0.66
  wmat[E5R][E5R][AM2]=0.66
  wmat[E5R][E6][AM2]=0.66
  wmat[E5R][I5L][AM2]=0    # ruled out by (Apicella et al., 2012) Fig. 7
  wmat[E5R][I5][AM2]=0.46  # (Apicella et al., 2012) Fig. 7E (weight = 2 x weight for E5B->I5)
  wmat[E6][E2][AM2]=0
  wmat[E6][E5B][AM2]=0.66
  wmat[E6][E5R][AM2]=0.66
  wmat[E6][E6][AM2]=0.66
  wmat[E6][I6L][AM2]=0.23
  wmat[E6][I6][AM2]=0.23
  wmat[I2L][E2][GA2]=0.83
  wmat[I2L][E2][GB2]=0.83 * GB2R
  wmat[I2L][I2L][GA2]=1.5
  wmat[I2L][I2][GA2]=1.5
  wmat[I2][E2][GA]=1.5
  wmat[I2][I2L][GA]=1.5
  wmat[I2][I2][GA]=1.5
  wmat[I5L][E5B][GA2]=0.83
  wmat[I5L][E5B][GB2]=0.83 * GB2R
  wmat[I5L][E5R][GA2]=0.83
  wmat[I5L][E5R][GB2]=0.83 * GB2R  
  wmat[I5L][I5L][GA2]=1.5
  wmat[I5L][I5][GA2]=1.5
  wmat[I5][E5B][GA]=1.5
  wmat[I5][E5R][GA]=1.5
  wmat[I5][I5L][GA]=1.5
  wmat[I5][I5][GA]=1.5
  wmat[I6L][E6][GA2]=0.83
  wmat[I6L][E6][GB2]=0.83 * GB2R
  wmat[I6L][I6L][GA2]=1.5
  wmat[I6L][I6][GA2]=1.5
  wmat[I6][E6][GA]=1.5
  wmat[I6][I6L][GA]=1.5
  wmat[I6][I6][GA]=1.5
  # gain control
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      for sy in range(STYPi):
        if wmat[ty1][ty2][sy] > 0:
          if ice(ty1): # I -> X
            if ice(ty2):
              if IsLTS(ty1): # LTS -> I
                if IsLTS(ty2): # LTS -> LTS
                  gn = IIGainLTSLTS
                else: # LTS -> FS
                  gn = IIGainLTSFS
              else: # FS -> I
                if IsLTS(ty2): # FS -> LTS
                  gn = IIGainFSLTS
                else: # FS -> FS
                  gn = IIGainFSFS
            else: # I -> E
              gn = IEGain
          else: # E -> X
            if ice(ty2): # E -> I
              if IsLTS(ty2): # E -> LTS
                gn = EIGainLTS
              else: # E -> FS
                gn = EIGainFS
            else: # E -> E
              gn = EEGain
              if sy==AM2: 
                wmat[ty1][ty2][AM] = wmat[ty1][ty2][AM2] * mGLURR
                if verbose: print('AM2:',wmat[ty1][ty2][AM2],'mGLURR:',wmat[ty1][ty2][AM])
            if sy==AM2:
              if ice(ty2): # E -> I
                wmat[ty1][ty2][NM2] = wmat[ty1][ty2][AM2] * NMAMREI
              else: # E -> E
                wmat[ty1][ty2][NM2] = wmat[ty1][ty2][AM2] * NMAMREE
          wmat[ty1][ty2][sy] *= gn 

# print weight matrix
def prwmat ():
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      for sy in range(STYPi):
        if wmat[ty1][ty2][sy] > 0:
          print("wmat[",CTYP[ty1],"][",CTYP[ty2],"][",STYP[sy],"]=",wmat[ty1][ty2][sy])

# print connection probability matrix
def prpmat ():
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      if pmat[ty1][ty2] > 0:
        print("pmat[",CTYP[ty1],"][",CTYP[ty2],"]=",pmat[ty1][ty2])

def setpmat ():
  for ii in range(CTYPi):
    for jj in range(CTYPi): pmat[ii][jj]=0
  pmat[E2][E2]=0.2     # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E2][E5B]=0.8    # strong by wiring matrix in (Weiler et al., 2008)
  pmat[E2][E5R]=0.8    # strong by wiring matrix in (Weiler et al., 2008)
  pmat[E2][I5L]=0.51   # L2/3 E -> L5 LTS I (justified by (Apicella et al., 2012)
  pmat[E2][E6]=0       # none by wiring matrix in (Weiler et al., 2008)
  pmat[E2][I2L]=0.51
  pmat[E2][I2]=0.43
  pmat[E5B][E2]=0          # none by wiring matrix in (Weiler et al., 2008)
  pmat[E5B][E5B]=0.04 * 4  # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 5 
  pmat[E5B][E5R]=0         # set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 5 
  pmat[E5B][E6]=0        #  none by suggestion of Ben and Gordon over phone
  pmat[E5B][I5L]=0         # ruled out by (Apicella et al., 2012) Fig. 7
  pmat[E5B][I5]=0.43 
  pmat[E5R][E2]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E5R][E5B]=0.21 * 4  # need to set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 5
  pmat[E5R][E5R]=0.11 * 4  # need to set using (Kiritani et al., 2012) Fig. 6D, Table 1, value x 5 
  pmat[E5R][E6]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E5R][I5L]=0         # ruled out by (Apicella et al., 2012) Fig. 7
  pmat[E5R][I5]=0.43
  pmat[E6][E2]=0           # none by wiring matrix in (Weiler et al., 2008)
  pmat[E6][E5B]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E6][E5R]=0.2        # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E6][E6]=0.2         # weak by wiring matrix in (Weiler et al., 2008)
  pmat[E6][I6L]=0.51
  pmat[E6][I6]=0.43
  pmat[I2L][E2]=0.35
  pmat[I2L][I2L]=0.09
  pmat[I2L][I2]=0.53
  pmat[I2][E2]=0.44
  pmat[I2][I2L]=0.34
  pmat[I2][I2]=0.62
  pmat[I5L][E5B]=0.35
  pmat[I5L][E5R]=0.35
  pmat[I5L][I5L]=0.09
  pmat[I5L][I5]=0.53
  pmat[I5][E5B]=0.44
  pmat[I5][E5R]=0.44
  pmat[I5][I5L]=0.34
  pmat[I5][I5]=0.62
  pmat[I6L][E6]=0.35
  pmat[I6L][I6L]=0.09
  pmat[I6L][I6]=0.53
  pmat[I6][E6]=0.44
  pmat[I6][I6L]=0.34
  pmat[I6][I6]=0.62
  for ii in range(CTYPi):
    for jj in range(CTYPi): pmat[ii][jj]*=pmatscale

numc = [0 for i in range(CTYPi)]; # number of cells of a type
ix = [0 for i in range(CTYPi)]; #starting index of a cell type (into self.ce list)
ixe = [0 for i in range(CTYPi)]; #ending index of a cell type
allcells,ecells,icells = 0,0,0
div = zeros( (CTYPi, CTYPi) )
conv = zeros( (CTYPi, CTYPi) )
syty1 = zeros( (CTYPi, CTYPi) ) # stores synapse codes (from labels.py)
syty2 = zeros( (CTYPi, CTYPi) ) # stores synapse code (from labels.py)
syty3 = zeros( (CTYPi, CTYPi) ) # stores synapse code (from labels.py)
sytys1 = {} # dictionary of synapse names
sytys2 = {} # dictionary of synapse names
sytys3 = {} # dictionary of synapse names
SOMA = 0; BDEND = 1; ADEND1 = 2; ADEND2 = 3; ADEND3 = 4;
dsecnames = ['soma','Bdend','Adend1','Adend2','Adend3']

def setdivmat ():
  import math
  for ty1 in range(CTYPi):
    for ty2 in range(CTYPi):
      if pmat[ty1][ty2] > 0.0: 
        div[ty1][ty2] =  math.ceil(pmat[ty1][ty2]*numc[ty2])
        conv[ty1][ty2] = int(0.5 + pmat[ty1][ty2]*numc[ty1])

# setup cell-type-to-cell-type synapse-type information
def setsyty ():
  for ty1 in range(CTYPi): # go thru presynaptic types
    for ty2 in range(CTYPi): # go thru postsynaptic types
      syty1[ty1][ty2] = syty2[ty1][ty2] = syty3[ty1][ty2] = -1 # initialize to invalid
      if numc[ty1] <= 0 or numc[ty2] <= 0: continue
      if ice(ty1): # is presynaptic type inhibitory?
        if IsLTS(ty1): # LTS -> X
          syty1[ty1][ty2] = GA2 # code for dendritic gabaa synapse
          if ice(ty2): # LTS -> Io
            sytys1[(ty1,ty2)] = "GABAss"
          else: # LTS -> E
            syty2[ty1][ty2] = GB2 # code for denritic gabab synapse
            sytys1[(ty1,ty2)] = "GABAs"
            sytys2[(ty1,ty2)] = "GABAB"
        else: # BAS -> X
          syty1[ty1][ty2] = GA # code for somatic gabaa synapse
          sytys1[(ty1,ty2)] = "GABAf"
      else: # E -> X
        syty1[ty1][ty2] = AM2 # code for dendritic ampa synapse
        syty2[ty1][ty2] = NM2 # code for dendritic nmda synapse
        if ice(ty2): # E -> I
          sytys1[(ty1,ty2)] = "AMPA"
          sytys2[(ty1,ty2)] = "NMDA"
        else: # E -> E
          sytys1[(ty1,ty2)] = "AMPA"
          sytys2[(ty1,ty2)] = "NMDA"
          sytys3[(ty1,ty2)] = "mGLUR"
          syty3[ty1][ty2] = AM # use AM for now -- really for mGLUR

lctyID,lctyClass = [],[]

# setup some convenient data structures
def setix (scale):
  import math
  global allcells,ecells,icells
  for i in range(CTYPi):
    numc[i] = int(math.ceil(cpernet[i]*scale))
    if numc[i] > 0:
      ty = PyrAdr
      if ice(i):
        if IsLTS(i): ty = LTS
        else: ty = FS
      for j in range(numc[i]):
        lctyClass.append(ty)
        lctyID.append(i)
      allcells += numc[i]
      if ice(i): icells += numc[i]
      else: ecells += numc[i]
  sidx = 0
  for i in range(CTYPi):
    if numc[i] > 0:
      ix[i] = sidx
      ixe[i] = ix[i] + numc[i] - 1
      sidx = ixe[i] + 1
  setdivmat()
  setsyty()

# setcellpos([pseed,network diameter in microns])
def setcellpos (pseed=4321,cside=colside):
  rdm = Random(); rdm.Random123(pseed, 0, 0)
  cellsnq = NQS("id","ty","ice","xloc","yloc","zloc")
  cellsnq.clear(allcells) # alloc space
  lX,lY,lZ=[],[],[]
  for i in range(allcells):    
    ctyp = lctyID[i]
    if ctyp == E2 or ctyp == I2 or ctyp == I2L: # If L2/3 cell...
      zmin = 143.0  # L2/3 upper bound (microns)
      zmax = 451.8  # L2/3 lower bound (microns)                
    elif ctyp == E5R: # Else, if L5 corticostriatal cell (goes in 5A or 5B)...
      zmin = 451.8   # L5A upper bound (microns)
      zmax = 1059.0  # L5B lower bound (microns)                
    elif ctyp == E5B: # Else, if L5 corticospinal cell (goes in 5B)...
      zmin = 663.6   # L5B upper bound (microns)
      zmax = 1059.0  # L5B lower bound (microns)                
    elif ctyp == I5 or ctyp == I5L: # Else, if L5 interneuron (goes in 5A or 5B)...
      zmin = 451.8   # L5A upper bound (microns)
      zmax = 1059.0  # L5B lower bound (microns)	        
    elif ctyp == E6 or ctyp == I6 or ctyp == I6L: # If L6 cell...
      zmin = 1059.0  # L6 upper bound (microns) 
      zmax = 1412.0  # L6 lower bound, white-matter (microns)
    [x,y,z]=[rdm.uniform(0,cside),rdm.uniform(0,cside),rdm.uniform(zmin,zmax)]
    cellsnq.append(i,ctyp,ice(ctyp),x,y,z); lX.append(x); lY.append(y); lZ.append(z);
  return cellsnq,lX,lY,lZ

setcpernet() # setup number of cells per network
setwmatex() # setup matrices of external inputs
setsynloc() # setup synapse location matrices
setdelmats() # setup delay matrices
setwmat() # setup weight matrix
setpmat() # setup connectivity matrix
setix(scale)
cellsnq,lX,lY,lZ=setcellpos()
ce = [] # cells on the host
gidvec = [] # gids of cells on the host
lncrec,ltimevec,lidvec=[],[],[] # spike recorders
dlids = {} # map from gid back to ce index

# create the cells
pcID = int(pc.id()); 
maxcells=0
cperhost = int(allcells/nhosts)
maxcells = cperhost
extra = allcells - cperhost*nhosts
if extra > 0: # check if any remainder cells
  if pcID < extra: # first hosts get extra cell
    maxcells += 1 # assign an extra cell if any remainders
    gid = pcID * (cperhost + 1)
  else: # rest of hosts do not
    gid = extra*(cperhost+1) + (pcID-extra) * cperhost
else: # even division? all hosts get equal cells
  gid = pcID * cperhost
for i in range(maxcells):
  ct = lctyID[gid]
  cell = lctyClass[gid](0+i*50,0,0,gid,ct)
  cell.x,cell.y,cell.z = lX[gid],lY[gid],lZ[gid]
  dlids[gid] = len(ce) # map from gid back to ce index
  ce.append(cell)
  gidvec.append(gid)
  pc.set_gid2node(gid,pcID)
  timevec,idvec = h.Vector(),h.Vector()
  ncrec = h.NetCon(ce[-1].soma(0.5)._ref_v, None, sec=ce[-1].soma)
  ncrec.record(timevec,idvec,gid)
  ncrec.threshold = spiketh # 10 mV is default, lower it for FS cells
  ltimevec.append(timevec); lidvec.append(idvec); lncrec.append(ncrec)
  pc.cell(gid,lncrec[-1],1) # 1 as 3rd arg means this cell can be source for events
  gid += 1
  
print(('  Number of cells on node %i: %i' % (pcID,len(ce))))
print('pcid:',pcID,'maxcells:',maxcells)
pc.barrier()

# next - do the wiring

# wire the network using NQS table
nccl = []
def wirenq (cnq):
  global nccl
  nccl = [] # NetCon list for connections between cells 
  cnq.tog("DB")
  vid1,vid2,vwt1,vwt2,vdel,vsec=cnq.getcol("id1"),cnq.getcol("id2"),cnq.getcol("wt1"),cnq.getcol("wt2"),cnq.getcol("del"),cnq.getcol("sec")
  vwt3 = cnq.getcol("wt3")
  for i in range(int(cnq.v[0].size())):
    prid = int(vid1[i])
    poid = int(vid2[i])
    if not pc.gid_exists(poid): continue # only make the connection on a node that has the target
    ty1 = lctyID[prid] # ce[prid].ty
    ty2 = lctyID[poid] # ce[poid].ty
    sname = dsecnames[int(vsec[i])] # which section is the synapse on?
    syn = sname + sytys1[(ty1,ty2)]
    wt1 = vwt1[i]
    delay = vdel[i]
    targ = ce[dlids[poid]]
    nc1 = pc.gid_connect(prid, targ.__dict__[syn].syn)
    nc1.delay = delay; nc1.weight[0] = wt1; nc1.threshold = spiketh; nccl.append(nc1)
    wt2 = vwt2[i]
    if wt2 > 0: # two synapses? (i.e., AMPA and NMDA)
      syn = sname + sytys2[(ty1,ty2)]
      nc2 = pc.gid_connect(prid, targ.__dict__[syn].syn)
      nc2.delay = delay; nc2.weight[0] = wt2; nc2.threshold = spiketh; nccl.append(nc2)
    wt3 = vwt3[i]
    if wt3 > 0: # three synapses? (i.e., AMPA and NMDA and mGLUR)
      if verbose: print('mGLUR synapse wt3 > 0:',wt3)
      syn = sname + sytys3[(ty1,ty2)]
      nc3 = pc.gid_connect(prid, targ.__dict__[syn].syn)
      nc3.delay = delay; nc3.weight[0] = wt3; nc3.threshold = spiketh; nccl.append(nc3)

#
def picksec (prty, poty, rdm):
  if ice(poty): return SOMA
  if ice(prty): # I -> E
    if IsLTS(prty): # LTS -> E
      if rdmsec: return rdm.discunif(ADEND1,ADEND3)
      else: return ADEND3
    else:
      return SOMA
  else: # E -> E
    if rdmsec: return rdm.discunif(BDEND,ADEND3)
    else: return ADEND3

# swire - spatial wiring: wires the network using pmat and cell positions
#                    (wiring probability affected by distance btwn cells)
#  slambda (global) specifies length-constant for spatially-dependent fall-off in wiring probability
#  at one slambda away, probability of connect is value in pmat
def swire (wseed):
  global slambda
  from math import sqrt,exp
  [vidx,vdel,vtmp,vwt1,vwt2,vwt3,vprob] = [Vector() for x in range(7)]
  z = 0
  if slambda <= 0:
    print("swire WARN: invalid slambda=", slambda, "setting slambda to ", colside/3)
    slambda=colside/3
  slambdasq = slambda**2 # using squared distance
  h.vrsz(1e4,vidx,vdel,vtmp)
  rdm=Random(); rdm.Random123(wseed, 0, 0) #initialize random # generator
  rdm.uniform(0,1)
  vprob.resize(allcells**2); vprob.setrand(rdm)
  pdx=0 # index into vprob
  connsnq=NQS("id1","id2","del","wt1","wt2","wt3","sec")
  connsnq.clear(1e3*allcells)
  for prid in range(allcells): 
    vrsz(0,vidx,vdel,vwt1,vwt2,vwt3)
    prty=lctyID[prid]
    ic1=ice(prty)
    for poty in range(0,CTYPi):
      if numc[poty] > 0 and pmat[prty][poty]>0:
        pbase = pmat[prty][poty]
        for poid in range(ix[poty],ixe[poty]+1): # go thru postsynaptic cells
          if prid==poid: continue # no self-connects
          ic2=ice(lctyID[poid])
          dx = lX[prid] - lX[poid]
          dy = lY[prid] - lY[poid]
          dz = lZ[prid] - lZ[poid]
          dsq = dx**2 + dy**2 # Connectivity fall-off only depends in intra-layer distance
          ds = sqrt(dsq + dz**2) # Axonal delay depends on all quantities
          prob = pbase * exp(-sqrt(dsq)/slambda) # probability of connect
          if prob >= vprob[pdx]: # rdm.uniform(0,1)
            mindelay = delm[prty][poty]-deld[prty][poty]
            maxdelay = delm[prty][poty]+deld[prty][poty]
            delay=rdm.uniform(mindelay,maxdelay) # synaptic delay
            delay += ds/axonalvelocity # add axonal delay 
            vidx.append(poid); vdel.append(delay)
            if syty1[prty][poty]>=0: vwt1.append(wmat[prty][poty][int(syty1[prty][poty])])
            else: vwt1.append(0)
            if syty2[prty][poty]>=0: vwt2.append(wmat[prty][poty][int(syty2[prty][poty])])
            else: vwt2.append(0)
            if syty3[prty][poty]>=0: vwt3.append(wmat[prty][poty][int(syty3[prty][poty])])
            else: vwt3.append(0)
          pdx += 1
    for ii in range(int(vidx.size())): connsnq.append(prid,vidx[ii],vdel[ii],vwt1[ii],vwt2[ii],vwt3[ii],picksec(prty , lctyID[int(vidx[ii])], rdm))
  wirenq(connsnq) # do the actual wiring based on self.connsnq
  return connsnq

if 'wnq' in dconf: # use a pre-specified wiring?
  connsnq = NQS(dconf['wnq'])
  wirenq(connsnq) # wire cells together with NQS
else:
  if wirety == SWIRE: connsnq=swire(WSEED)
  else: connsnq = None

pc.barrier() # wait for wiring to get completed

# setup rxd if it's ON
# get list of all Sections associated with an excitatory cell
def getesec ():
  esec = []
  for cell in ce:
    if ice(cell.ty): continue
    for s in cell.all_sec: esec.append(s)
  return esec

# get list of sections to use for rxd
def getrxdsec ():
  rxdsec = getesec() # E cell sections
  return rxdsec
  
def pcidpr (s): 
  global pcID
  print('host',pcID,':',s)

### RXD ###
rxdsec = [] # Section list for use with rxd
[cyt,er,cyt_er_membrane,ca,caextrude,serca,leak,CB,caCB,buffering]=[None for i in range(10)]
rxdsec=getrxdsec() # using rxd ?
pc.barrier()
if len(rxdsec) > 0: # only use rxd if there are viable Sections
  pcidpr('setting up rxd...')
  from neuron import rxd
  pcidpr('imported rxd...')
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
  ip3degTau = 1000 # ms for ip3 degrade
  ip3deg = rxd.Rate(ip3, (0.0-ip3[cyt])/ip3degTau, regions=cyt, membrane_flux=False) 

  def setmGLURflux (): # mGLUR synapses generate ip3 that is fed into rxd ip3
    thistime = time.time()
    for c in ce:
      if ice(c.ty): continue
      for syn,seg in zip([c.Adend3mGLUR.syn,c.Adend2mGLUR.syn,c.Adend1mGLUR.syn],[c.Adend3(0.5), c.Adend2(0.5), c.Adend1(0.5)]):
        for node in ip3.nodes(seg):
          node.include_flux(syn._ref_rip3)
          if time.time() - thistime > 60:
            thistime = time.time()
            print(str(c)+" "+str(syn)+" "+str(seg)+" "+str(node))

  def setrecip3 ():
    for c in ce:
      if ice(c.ty): continue
      c.soma_ip3cyt = Vector(tstop/h.dt)
      c.soma_ip3cyt.record( ip3[cyt].nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_ip3cyt = Vector(tstop/h.dt)
      c.Adend3_ip3cyt.record( ip3[cyt].nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )
    
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
  print("serca done")
  gleak = dconf['gleak']   # leak channel: bidirectional ca flow btwn cyt <> ER
  leak = rxd.MultiCompartmentReaction(ca[er]!=ca[cyt], gleak, gleak, membrane=cyt_er_membrane)
  print("leak done")
  
  def setreccaer (): # setup recording of ca[er] for each pyramidal cell in Adend3,soma center
    for c in ce:
      if ice(c.ty): continue
      c.soma_caer = Vector(tstop/h.dt)
      c.soma_caer.record( ca[er].nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_caer = Vector(tstop/h.dt)
      c.Adend3_caer.record( ca[er].nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )

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
      c.soma_caCB.record( caCB.nodes(c.soma)(0.5)[0]._ref_concentration, recdt )
      c.Adend3_caCB = Vector(tstop/h.dt)
      c.Adend3_caCB.record( caCB.nodes(c.Adend3)(0.5)[0]._ref_concentration, recdt )

  setreccaer() # NB: only record from RXD variables after ALL species setup!
  print("setreccaer done")
  setreccacb() # otherwise, the pointers get messed up.
  print("setreccacb done")
  setrecip3()
  print("setrecip3 done")
  setmGLURflux()
  print("setmGLURflux done")

#Tuomo: find the name of soma and set the hot zone calcium channels
somasec = 0
hotzonetype = 1 #0:hot zone in second to last. 1:hot zone in last. 2: hot zone in two last ones
hotzonecoeff = 100.0

for i in range(0,len(lctyID)):
  if i==0:
    print("dlids.keys() = "+str(list(dlids.keys())))
  #print("rank "+str(MPI.COMM_WORLD.Get_rank())+", pcID "+str(pcID)+" lctyID["+str(i)+"] = "+str(lctyID[i])+", len(lctyID) = "+str(len(lctyID))+", len(dlids) = "+str(len(list(dlids.keys())))+", dlids.keys()[0]="+str(list(dlids.keys())[0]))
  if i in list(dlids.keys()):
    ce_ind = dlids[i]
    print("i = "+str(i)+", ce_ind = "+str(ce_ind)+", rank "+str(MPI.COMM_WORLD.Get_rank())+", pcID "+str(pcID))
    if lctyID[i] in [E5R,E5B]:
      print("Setting hot zone in ce["+str(i)+"]")
      if hotzonetype in [0,2]:
        h(ce[ce_ind].Adend2.name()+' gcalbar_cal = gcalbar_cal * '+str(hotzonecoeff))
        h(ce[ce_ind].Adend2.name()+' gcatbar_cat = gcatbar_cat * '+str(hotzonecoeff))
      if hotzonetype in [1,2]:
        h(ce[ce_ind].Adend3.name()+' gcalbar_cal = gcalbar_cal * '+str(hotzonecoeff))
        h(ce[ce_ind].Adend3.name()+' gcatbar_cat = gcatbar_cat * '+str(hotzonecoeff))
  
# setup inputs - first noise inputs
def getsyns ():
  syns = {} # mapping of synapse names, first index is ice, second is synapse code
  syns[ (0,AM) ] = ["Adend3mGLUR","Adend2mGLUR","Adend1mGLUR"]
  syns[ (0,AM2) ] = ["Adend3AMPA","Adend2AMPA","Adend1AMPA","BdendAMPA"]
  syns[ (1,AM2) ] = "somaAMPA"
  syns[ (0,NM2) ] = ["Adend3NMDA","Adend2NMDA","Adend1NMDA","BdendNMDA"]
  syns[ (1,NM2) ] = "somaNMDA"
  syns[ (0,GB2) ] = ["Adend3GABAB","Adend2GABAB","Adend1GABAB"]
  syns[ (0,GA2) ] = ["Adend3GABAs","Adend2GABAs","Adend1GABAs"]
  syns[ (1,GA2) ] = "somaGABAss"
  syns[ (0,GA) ] = "somaGABAf"
  syns[ (1,GA) ] = "somaGABAf"
  return syns

dsstr = ['AMPA', 'NMDA', 'GABAS', 'mGLUR', 'GABAB']

# get a recording vector for the synaptic current
def recveccurr (sy):
  if not dconf['saveExtra']: return h.Vector()
  vec = h.Vector()
  try:
    vec.record(sy.syn._ref_i,recvdt)
  except:
    vec.record(sy.syn._ref_iNMDA,recvdt)
  return vec

# adds synapses across dendritic fields for the E cells
def addsyns ():
  for cell in ce:
    cell.dsy = {}; cell.vsy = {}
    if ice(cell.ty): continue
    ds = {}; ds[cell.Adend3]='Adend3'; ds[cell.Adend2]='Adend2'; ds[cell.Adend1]='Adend1'; ds[cell.Bdend]='Bdend'
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1, cell.Bdend]:
      llsy = [];
      nloc = sec.nseg
      llvsy = []; # for recording currents
      for i,seg in enumerate(sec):
        if seg.x == 0.0 or seg.x == 1.0: continue # skip endpoints
        lsy = []; loc = seg.x; lvsy = [] #AMPA, NMDA, GABAA_slow, GABAB
        #print 'loc:',loc
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

addsyns()
print("addsyns done")
#creates NetStims (and associated NetCon,Random) - provide 'noise' inputs
#returns next useable value of sead
def makeNoiseNetStim (cel,nsl,ncl,nrl,nrlsead,syn,w,ISI,time_limit,sead):
  ns = h.NetStim()
  ns.interval = ISI
  ns.noise = 1			
  ns.number = 2 * time_limit / ISI  # create enough spikes for extra time, in case goes over limit
  ns.start = tstart
  if type(syn) == str: nc = h.NetCon(ns,cel.__dict__[syn].syn)
  else: nc = h.NetCon(ns,syn)
  nc.delay = h.dt * 2 # 0
  nc.weight[0] = w
  rds = h.Random()
  rds.Random123(sead, 0, 0)
  rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
  ns.noiseFromRandom(rds)
  ns.start = 0.
  nsl.append(ns)
  ncl.append(nc)
  nrl.append(rds)
  nrlsead.append(sead)
  cel.infod[syn] = [ns,nc] #store pointers to NetStim,NetCon for easier manipulation

def makeNoiseNetStims (simdur,rdmseed):
  print("Making NetStims")
  nsl = [] #NetStim List
  ncl = [] #NetCon List
  nrl = [] #Random List for NetStims
  nrlsead = [] #List of seeds for NetStim randoms
  syns = getsyns() ; 
  for cell in ce: # go through cell types, check weights,rates of inputs
    ct = cell.ty # get cell type code
    if ice(ct): # only has 1 compartment
      for sy in range(STYPi):
        if wmatex[ct][sy] <= 0.0 or ratex[ct][sy] <= 0: continue
        syn = syns[(ice(ct),sy)]
        if type(syn) == list:
          for idx,SYN in enumerate(syn):
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,SYN,wmatex[ct][sy],1e3/ratex[ct][sy],simdur,rdmseed*(cell.ID+1)*(idx+1))
        else:
          makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,syn,wmatex[ct][sy],1e3/ratex[ct][sy],simdur,rdmseed*(cell.ID+1))
    else: # E cells - need to distribute noise over all sections
      for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
        llsy = cell.dsy[sec]
        for lsy in llsy:
          for i,sy in enumerate([AM2,NM2,GA2,AM,GB2]):
            if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0. and sy != AM: # AM now for soma AMPA
              makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+1));
      sec = cell.Bdend; llsy = cell.dsy[sec];
      for lsy in llsy:
        for i,sy in enumerate([AM2,NM2,GA2]):
          if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0. and sy != AM: # AM now for soma AMPA
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+4)); 
      sec = cell.soma; llsy = cell.dsy[sec];
      for i,sy in enumerate([GA,AM,NM]):
        if ratex[ct][sy] > 0. and wmatex[ct][sy] > 0.:
          for lsy in llsy:
            makeNoiseNetStim(cell,nsl,ncl,nrl,nrlsead,lsy[i].syn,wmatex[ct][sy],(1e3/ratex[ct][sy]),simdur,rdmseed*(cell.ID+1)*(i+7)); rdmseed+=1
  return nsl,ncl,nrl,nrlsead

nsl,ncl,nrl,nrlsead = makeNoiseNetStims(tstart+tstop,ISEED)

pc.barrier() # wait for completion of NetStim creation
print("barrier done")

# setup wm-related stimuli
useWM = dconf['useWM']
rnnsl = [] # based on row number -- should always have the same number of entries as rows in the nqm 
rnncl = [] # based on row number times number of cells in the row
nqm=None; rnrds = []; rnseed = []; fihSIGns=None; 
def setupWMInputs ():
  global nqm,rnnsl,rnncl,rnrds,rnseed,fihSIGns
  if pcID == 0: print('setting up WM inputs')
  if os.path.exists(dconf['nqm']): # read nqm if it exists
    nqm = h.NQS(dconf['nqm'])
  else: # setup nqm here if nqm file path non-existant
    nqm = h.NQS('vid','startt','endt','rate','w'); nqm.odec('vid')
    lvid = []; nMem = dconf['nMem']; 
    memstartt = dconf['memstartt']; memintert = dconf['memintert']; memdurt = dconf['memdurt']
    memstartt += tstart # this only has an effect when loadstate != 0
    startt = memstartt; endt = startt+memdurt
    if verbose and pcID==0: print('startt:',startt,',memstartt:',memstartt,',tstart:',tstart,',loadstate:',loadstate)
    for i in range(nMem): # number of stimuli
      lvid.append(h.Vector()); vtmp=Vector()
      for ctdx,ct in enumerate([CTYP.index(x) for x in dconf['pops'].split(',')]): # default of E2,E5R,E5B,E6
        memfrac = float(dconf['memfrac'].split(',')[ctdx])
        if dconf['memSame']: # same population stim'ed each time? 
          vtmp.indgen(ix[ct],int(ix[ct]+numc[ct]*memfrac-1),1)
        else:
          if i % 2 == 1: # odd population
            vtmp.indgen(int(ix[ct]+numc[ct]*memfrac),ixe[ct],1)
          else: # even population
            vtmp.indgen(ix[ct],int(ix[ct]+numc[ct]*memfrac-1),1)
        lvid[-1].append(vtmp)
      nqm.append(lvid[-1],startt,endt,dconf['memrate'],dconf['memW'])
      startt += (memdurt+memintert);
      if i == nMem-2 and dconf['lastmemdurt'] != memdurt:
        endt = startt + dconf['lastmemdurt']
      else:
        endt = startt + memdurt;

  if pcID == 0: # backup the nqm file for later retrieval
    fnqwm = 'meminput/' + simstr + '_nqm.nqs'
    if os.path.exists(fnqwm):
      print('removing prior nqm file')
      cmd = 'rm ' + fnqwm
      os.system(cmd)
    print('backing up nqm to ' , fnqwm)
    nqm.sv(fnqwm)
    if verbose:
      print('this is nqm:')
      nqm.pr()

  def getslist (cell,syn):
    sidx = -1
    if ice(cell.ty): return [cell.soma], sidx
    slist=[cell.Adend3, cell.Adend2, cell.Adend1]
    if syn.count('AMPA'): sidx=0
    elif syn.count('NMDA'): sidx=1
    elif syn.count('GABAss'): sidx=2
    elif syn.count('mGLUR'): sidx=3
    if syn.count('0'): slist=[cell.Adend3, cell.Adend2, cell.Adend1]
    elif syn.count('1'): slist=[cell.Adend1]
    elif syn.count('2'): slist=[cell.Adend2]
    elif syn.count('3'): slist=[cell.Adend3]
    return slist, sidx

  def createNetStims (vid,rate,w,startt,endt,seed=1234,syn='Adend3AMPA'):
    global rnnsl,rnncl,rnrds,rnseed # based on row number in nqm
    rnnsl.append( [] ) # a list for each row of nqm
    rnncl.append( [] )
    rnrds.append( [] )
    rnseed.append( [] )
    for idx in vid:
      if not pc.gid_exists(idx): continue # only make the connection on a node that has the target
      if verbose: print('idx is ' , idx, ' dlid = ' , dlids[idx], ' type is ' , CTYP[lctyID[int(idx)]])
      cell = ce[dlids[idx]]
      if ice(cell.ty) and syn.count('mGLUR') > 0: continue
      slist, sidx = getslist(cell, syn)
      if sidx == -1:
        ns = h.NetStim()
        ns.start = startt
        ns.number = (endt-startt) * rate / 1e3
        if verbose: print('createNetStims:',startt,endt,ns.number,rate,w)
        ns.interval = 1e3 / rate
        ns.noise = 1
        rnnsl[-1].append(ns)
        if ice(cell.ty):
          if syn.count('AMPA'):
            nc = h.NetCon(ns,cell.__dict__['somaAMPA'].syn)
          elif syn.count('NMDA'):
            nc = h.NetCon(ns,cell.__dict__['somaNMDA'].syn)
          if IsLTS(cell.ty):
            nc.weight[0] = w / 4.0
          else:
            nc.weight[0] = w / 4.0
        else:
          nc = h.NetCon(ns,cell.__dict__[syn].syn)
          nc.weight[0] = w
        nc.delay = h.dt * 2
        rnncl[-1].append(nc)

        rds = h.Random()
        rds.Random123(seed, 1, idx)
        rds.negexp(1)
        ns.noiseFromRandom(rds)
        rnrds[-1].append(rds)
        rnseed[-1].append(seed*idx)
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
            rds.Random123(tmpseed, 1, idx)
            rds.negexp(1)
            ns.noiseFromRandom(rds)
            rnrds[-1].append(rds)
            rnseed[-1].append(tmpseed)
            tmpseed += 1
    if verbose: vid.printf()

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
    nqm.tog("DB")
    sz = int(nqm.v[0].size()) # number of WM stims
    lastmGLURON = dconf['lastmGLURON']
    lapicIDX = []
    try:
      lapicIDX = [int(dconf['apicIDX'])]
      if lapicIDX[0] == 0: lapicIDX = [1,2,3] # 0 means all apical dends
    except:
      for i in dconf['apicIDX'].split(','): lapicIDX.append(int(i))     
    if lastmGLURON:
      print('lastmGLURON!')
      for i in range(sz):
        if i == sz-1:
          for j in lapicIDX:
            createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'mGLUR',wfctr=mGLURRWM)
        else:
          for j in lapicIDX:
            createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'AMPA',wfctr=AMRWM)
            createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'NMDA',wfctr=NMAMRWM)
    else:
      for i in range(sz):
        for j in lapicIDX:
          createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'AMPA',wfctr=AMRWM)
          createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'NMDA',wfctr=NMAMRWM)
          createNetStimsFromNQ(nqm,i,seed=(j+1)*(i+1)*12345,syn='Adend'+str(j)+'mGLUR',wfctr=mGLURRWM)

  def initSigNetStims ():
    if verbose: print('in initSigNetStims')
    for i in range(len(rnnsl)):
      for j in range(len(rnnsl[i])):
        rds = rnrds[i][j]
        sead = rnseed[i][j]
        rds.Random123(sead, 2, i*100000+j)
        rds.negexp(1)
  fihSIGns = h.FInitializeHandler(0, initSigNetStims)
  setStims() # create the inputs based on contents of nqm
  print("setStims done")

#this should be called @ beginning of each sim - done in an FInitializeHandler
def init_NetStims ():
  print('node ' , pc.id() , ' in init_NetStims')
  for i in range(len(nrl)):
    rds = nrl[i]
    sead = nrlsead[i]
    rds.Random123(sead, 1, i)
    rds.negexp(1)
fihns = h.FInitializeHandler(0, init_NetStims)
print("h.FInitializeHandler done")

# handler for printing out time during simulation run
def fi():
  for i in range(int(tstart),int(tstart+tstop),100): h.cvode.event(i, "print " + str(i))

if pc.id() == 0: fih = h.FInitializeHandler(1, fi)

vt=Vector(); vt.record(h._ref_t);

pc.barrier() # wait for NetStims to get setup 
print("barrier done")


if dconf['cvodeactive']:
  print('cvode on')
  h.cvode.active(1) # turn on variable time-step integration
  h.cvode.atol(1e-8)
  h.cvode.rtol(1e-8)

mechanisms = ['pas','iar','nax','kdr','kap','cal','can','cat','ikc','km','cagk','kdmc','im','calts','ihlts','icalts','kcalts','Kdrbwb','Nafbwb','cal','SK_E2']
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
  elif dconfkeys[ikey] == 'gCoeff_gbar_GABAB': #GABAB currents are only in the pyramidal cells
    for isyn in range(0,len(h.GABAB)):
      h.GABAB[isyn].gbar = h.GABAB[isyn].gbar * float(dconf[dconfkeys[ikey]])
  elif dconfkeys[ikey] == 'Cm_apic':
    for sec in [cell.Adend3, cell.Adend2, cell.Adend1]:
      sec.cm = float(dconf[dconfkeys[ikey]])
print("dconfkeys done")
  
####################################################################################
### simulation run here 
def myrun ():
  pc.set_maxstep(10)
  dastart,daend=None,None
  if pc.id()==0:
    dastart = datetime.now()
    print('started at:',dastart)
  if useWM: setupWMInputs()
  h.stdinit()
  if len(rxdsec)>0: # any sections with rxd?
    ca[er].concentration = dconf['caerinit'] # 100e-6
    ca[cyt].concentration = dconf['cacytinit'] # 100e-6
  for i in range(0,len(h.GABAB)):
    print("h.GABAB["+str(i)+"].gbar = "+str(h.GABAB[i].gbar))
  if pcID==0: print('current time is ' , h.t)
  pc.psolve(h.t+tstop) # run for tstop
  pc.barrier() # Wait for all hosts to get to this point
  if pc.id()==0:
    daend = datetime.now()
    print('finished ',tstop,' ms sim at:',daend)
    dadiff = daend - dastart;
    print('runtime:',dadiff, '(',tstop/1e3,' s)')

if dconf['dorun']: myrun()
print("myrun done")

# concatenate the results so can view/save all at once
lspks,lids=array([]),array([])
for host in range(nhosts): # is this loop required? can't just post messages from given host?
  if host == pc.id():
    for i in range(len(ltimevec)):
      try:
        lspks=concatenate((lspks,array(ltimevec[i])))
        lids=concatenate((lids,array(lidvec[i])))
      except:
        print("couldn't concatenate, host="+str(host)+", i="+str(i))

# save data - output path based on simstr and pcid
def savedata (simstr,pcid):
  safemkdir(outdir)
  fn = outdir + '/' + simstr + '_pc_' + str(pcid) + '.npz'
  print('host ' , pcid, ' saving to ' , fn)
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
    Adend2_iNM = zeros((ne,szfast)); Adend1_iNM = zeros((ne,szfast))
    Bdend_iNM = zeros((ne,szfast)); soma_iNM = zeros((ne,szfast))
    soma_iGA = zeros((ne,szfast))
    soma_iNM_i = zeros((ni,szfast))
    [Adend3_ina, Adend3_ik, Adend3_ica, Adend3_ih] = [zeros((ne,szfast)) for i in range(4)]
    [soma_ina, soma_ik, soma_ica, soma_ih] = [zeros((ne,szfast)) for i in range(4)]
  cdx = 0; idx = 0;
  for c in ce:
    if ice(c.ty):
      soma_voltI[idx,:] = c.soma_volt.to_python()
      if saveExtra:
        soma_iNM_i[idx,:] = c.soma_iNM.to_python()
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
      Adend2_iNM[cdx,:] = c.Adend2_iNM.to_python()
      Adend1_iNM[cdx,:] = c.Adend1_iNM.to_python()
      Bdend_iNM[cdx,:] = c.Bdend_iNM.to_python()
      soma_iNM[cdx,:] = c.soma_iNM.to_python()
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
    numpy.savez(fn,lctyID=array(lctyID),lX=array(lX),lY=array(lY),lZ=array(lZ),vt=vt.as_numpy(),lspks=lspks,\
                  lids=lids,lE=lE,lI=lI,soma_volt=soma_volt,soma_voltI=soma_voltI,Adend3_volt=Adend3_volt,Bdend_volt=Bdend_volt,\
                  soma_cai=soma_cai,Adend3_cai=Adend3_cai,soma_Ihm=soma_Ihm,Adend3_Ihm=Adend3_Ihm,soma_Ihp1=soma_Ihp1,\
                  Adend3_Ihp1=Adend3_Ihp1,soma_caer=soma_caer,Adend3_caer=Adend3_caer,soma_caCB=soma_caCB,Adend3_caCB=Adend3_caCB,\
                  Adend3_ip3cyt=Adend3_ip3cyt,soma_ip3cyt=soma_ip3cyt,Adend3_iAM=Adend3_iAM,Adend3_iNM=Adend3_iNM,Adend3_iGB=Adend3_iGB,\
                  Adend3_iGA=Adend3_iGA,soma_iGA=soma_iGA,Adend3_ina=Adend3_ina,Adend3_ik=Adend3_ik,Adend3_ica=Adend3_ica,Adend3_ih=Adend3_ih,\
                  soma_ina=soma_ina,soma_ik=soma_ik,soma_ica=soma_ica,soma_ih=soma_ih,soma_iNM=soma_iNM,Bdend_iNM=Bdend_iNM,Adend1_iNM=Adend1_iNM,Adend2_iNM=Adend2_iNM,soma_iNM_i=soma_iNM_i)
  else: # no somaVolt, synaptic currents
    numpy.savez(fn,lctyID=array(lctyID),lX=array(lX),lY=array(lY),lZ=array(lZ),vt=vt.as_numpy(),lspks=lspks,lids=lids,lE=lE,lI=lI,Adend3_volt=Adend3_volt,Bdend_volt=Bdend_volt,soma_cai=soma_cai,Adend3_cai=Adend3_cai,soma_Ihm=soma_Ihm,Adend3_Ihm=Adend3_Ihm,soma_Ihp1=soma_Ihp1,Adend3_Ihp1=Adend3_Ihp1,soma_caer=soma_caer,Adend3_caer=Adend3_caer,soma_caCB=soma_caCB,Adend3_caCB=Adend3_caCB,Adend3_ip3cyt=Adend3_ip3cyt,soma_ip3cyt=soma_ip3cyt)

pc.barrier()
####################################################################################
if saveout: # save the sim data
  if pcID == 0: print('saving data')
  savedata(simstr,pcID)

pc.runworker()
pc.done()

if nhosts > 1: h.quit() # this means was likely running in batch mode

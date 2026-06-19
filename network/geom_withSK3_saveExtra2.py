#cp geom_withSK_COEFFL1_2_N0_3_T0_84_SOMA0_25_DEND0_25.py geom_withSK3.py
from math import exp
import sys
from pyinit import *
from labels import *
h.celsius = 37
h('load_file("pywrap.hoc")')
from conf import *
# determine config file name
def setfcfg ():
  fcfg = "netcfg.cfg" # default config file name
  for i in range(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  print("config file is " , fcfg)
  return fcfg

fcfg=setfcfg() # config file name
dconf = readconf(fcfg)
print("Overwriting cagk_gbar by 0.0012*36 - geom_calhayact_withSK.py")
dconf['cagk_gbar'] = 0.0012*36
taurcada = dconf['taurcada']
h.cac_iar = 0.006
h.k4_iar = dconf['iark4'] # original: 0.008 
ihginc = dconf['ihginc']; h.ginc_iar = ihginc
ihscale = dconf['ihscale']
ilscale = dconf['ilscale']
ikmscale = dconf['ikmscale']
recdt = dconf['recdt']
recvdt = dconf['recvdt']
saveExtra = 1 #dconf['saveExtra'] # whether to save extra variables (synaptic currents, soma voltage, etc.)
erevh = dconf['erevh']
spaceum = dconf['spaceum']
h_lambda = dconf['h_lambda']
h_gbar = dconf['h_gbar'] # orig: 0.00002
cagk_gbar = dconf['cagk_gbar'] # orig: 0.009
ikc_gkbar = dconf['ikc_gkbar'] # orig: 0.003
expk = dconf['expk']
cabar = dconf['cabar'] # 0.001
tau1NMDAEE=15; tau2NMDAEE=150;
tau1NMDAEI=15; tau2NMDAEI=150;
##

# if rdt > 0 use fixed interval for recording, else let cvode determine it
def saferecord (var, rdt):
  if rdt > 0.0:
    vrec = h.Vector(h.tstop/rdt + 1)
    vrec.record(var,rdt)
  else:
    vrec = h.Vector()
    vrec.record(var)
  return vrec

# metabotropic glutamate receptor
class SynapsemGLUR:
  def __init__(self,sect,loc):
    self.syn = h.mGLUR(loc, sec=sect)

class Synapse:
  def __init__(self, sect, loc, tau1, tau2, e):
    self.syn		= h.MyExp2SynBB(loc, sec=sect)
    self.syn.tau1	= tau1
    self.syn.tau2	= tau2
    self.syn.e		= e 
		
class SynapseNMDA:
  def __init__(self, sect, loc, tau1NMDA, tau2NMDA, r, e):
    self.syn			= h.MyExp2SynNMDABB(loc, sec=sect)
    self.syn.tau1NMDA	= tau1NMDA
    self.syn.tau2NMDA	= tau2NMDA 
    self.syn.r			= r
    self.syn.e			= e 

# gabab based on 1995 PNAS paper by Destexhe
class SynapseGABAB:
  def __init__(self, sect, loc):
    self.syn = h.GABAB(loc, sec=sect)
		
###############################################################################
# General Cell
###############################################################################
class Cell:
  "General cell"
  def __init__ (self,x,y,z,ID,ty):
    self.x=x
    self.y=y
    self.z=z
    self.ID=ID
    self.ty = ty
    self.snames = [] # list of section names
    self.all_sec = []
    self.add_comp('soma',True)
    self.set_morphology()
    self.set_conductances()
    self.set_synapses()
    self.set_inj()
    self.calc_area()
    self.infod = {} # dictionary for storing indices into nsl,ncl
    self.poID = [] # list of postsynaptic IDs (indices into Network's ce)
    self.poNC = [] # list of pointers to postsynaptic NetCons
    self.poSY = [] # synapse type (code from labels.py)

  # saves information on a synapse to another cell
  #  poid is postsynaptic id, nc is NetCon, syty is synapse code (from labels.py)
  def savesyinfo (self,poid,nc,syty):
    self.poID.append(poid)
    self.poNC.append(nc)
    self.poSY.append(syty)
    
  # get number of outgoing connections
  def getdvi (self): return len(self.poID)
  def set_morphology (self): pass			
  def set_conductances (self): pass		
  def set_synapses (self): pass		
  def set_inj (self): self.somaInj = h.IClamp(0.5, sec=self.soma)	
		
  def add_comp (self, name, rec):
    self.snames.append( name )
    self.__dict__[name] = h.Section()
    self.all_sec.append(self.__dict__[name])
    #self.all_sec_ref.append(h.SectionRef(sec=self.all_sec[-1]))
    # Record voltage
    if rec:
      self.__dict__[name+"_volt"] = saferecord(self.__dict__[name](0.5)._ref_v, recvdt)
      self.__dict__[name+"_volt"].label(name+"_volt")
      
  def plot_volt (self, name, fig=1):
    figure(fig)
    volt = self.__dict__[name+"_volt"].to_python()
    plot(arange(len(volt))*h.dt, volt)
		
  def calc_area (self):
    self.total_area = 0
    self.n = 0
    for sect in self.all_sec:
      self.total_area += h.area(0.5,sec=sect)
      self.n+=1
      
###############################################################################
# Soma-targeting interneuron (fast-spiking Basket Cell -- Bas)
###############################################################################
class Bas (Cell):
  "Basket cell"	
  def set_morphology(self):
    total_area = 10000 # um2
    self.soma.nseg  = 1
    self.soma.cm    = 1      # uF/cm2
    diam = sqrt(total_area) # um
    L    = diam/pi  # um			
    h.pt3dclear(sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
			
  def set_conductances(self):
    self.soma.insert('pas')
    self.soma.e_pas = -65     # mV
    self.soma.g_pas = 0.1e-3  # S/cm2 
    self.soma.insert('Nafbwb')
    self.soma.insert('Kdrbwb')
	   
  def set_synapses(self):
    self.somaGABAf=Synapse(sect=self.soma,loc=0.5,tau1=0.07,tau2=9.1,e=-80);
    self.somaGABAss=Synapse(sect=self.soma,loc=0.5,tau1=20,tau2=40,e=-80);
    self.somaAMPA=Synapse(sect=self.soma,loc=0.5,tau1=0.05,tau2=5.3,e=0);
    self.somaNMDA=SynapseNMDA(sect=self.soma,loc=0.5, tau1NMDA=tau1NMDAEI,tau2NMDA=tau2NMDAEI,r=1,e=0);
    if saveExtra: self.reccurr() #This should be in init function, but this is easier

  def reccurr (self):
    self.soma_iNM = saferecord(self.somaNMDA.syn._ref_iNMDA, recvdt)

###############################################################################
# Dendrite-targeting interneuron (LTS Cell)
###############################################################################
class Lts (Cell):
  "LTS cell"   
  def set_morphology(self):
    total_area = 10000 # um2
    self.soma.nseg  = 1
    self.soma.cm    = 1      # uF/cm2
    diam = sqrt(total_area) # um
    L    = diam/pi  # um
    h.pt3dclear(sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
    h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
	
  def set_conductances(self):
    self.soma.insert('pas')
    self.soma.e_pas = -65     # mV
    self.soma.g_pas = 0.1e-3  # S/cm2 
    self.soma.insert('Nafbwb')
    self.soma.insert('Kdrbwb')
    self.soma.insert('icalts')
    self.soma.insert('kcalts')
    self.soma.insert('ihlts')
    self.soma.insert('calts') # calcium extrusion
    
  def set_synapses(self):
    self.somaGABAf 	= Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
    self.somaGABAss	= Synapse(    sect=self.soma, loc=0.5, tau1=20,	  tau2=40, e=-80)
    self.somaAMPA 	= Synapse(    sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
    self.somaNMDA 	= SynapseNMDA(sect=self.soma, loc=0.5, tau1NMDA=tau1NMDAEI, tau2NMDA=tau2NMDAEI, r=1, e=0)
    if saveExtra: self.reccurr() #This should be in init function, but this is easier

  def reccurr (self):
    self.soma_iNM = saferecord(self.somaNMDA.syn._ref_iNMDA, recvdt)
    
LTS = Lts
FS = Bas
		
###############################################################################
# Pyramidal Cell
###############################################################################
class PyrAdr (Cell):
  "Pyramidal cell"
  def __init__(self,x,y,z,ID,ty):
    Cell.__init__(self,x,y,z,ID,ty)
    self.set_props()
    lrec = ['soma','Adend3']
    if saveExtra:
      lrec.append('Adend2'); lrec.append('Adend1'); lrec.append('Bdend')
    for sec in lrec:
      self.reccai(sec); self.recIhm(sec); self.recIhp1(sec)
    if saveExtra: self.reccurr()

  # turn on recording of Ih m in section with given name
  def recIhm (self, name):
    sn = name + "_Ihm"
    self.__dict__[sn] = saferecord(self.__dict__[name](0.5).iar._ref_m, recdt)
    self.__dict__[sn].label(sn)

  # turn on recording of Ih p1 in section with given name
  def recIhp1 (self, name):
    sn = name + "_Ihp1"
    self.__dict__[sn] = saferecord(self.__dict__[name](0.5).iar._ref_p1, recdt)
    self.__dict__[sn].label(sn)

  # turn on recording of calcium concentration in section with given name
  def reccai (self, name):
    sn = name + "_cai"
    self.__dict__[sn] = saferecord(self.__dict__[name](0.5)._ref_cai, recdt)
    self.__dict__[sn].label(sn)

  def set_morphology(self):
    self.add_comp('Bdend',True)
    self.add_comp('Adend1',saveExtra)
    self.add_comp('Adend2',saveExtra)
    self.add_comp('Adend3',True)
    self.apic = [self.Adend1, self.Adend2, self.Adend3]
    self.basal = [self.Bdend]
    sec = self.soma; sec.L = 20.0; sec.diam = 20.0
    if self.ty == E5R or self.ty == E5B: apicL = 300.0
    else: apicL = 150.0
    for sec in self.apic:
      sec.L = apicL; sec.diam = 2.0
    self.Bdend.L = 200.0; self.Bdend.diam = 2.0

    self.Bdend.connect(self.soma,    0, 0)
    self.Adend1.connect(self.soma,   1, 0)
    self.Adend2.connect(self.Adend1, 1, 0)
    self.Adend3.connect(self.Adend2, 1, 0)

    if spaceum > 0.0:
      for sec in self.all_sec:
        ns = int(sec.L / spaceum)
        if ns % 2 == 0: ns += 1
        sec.nseg = ns

  def set_props (self):
    Vrest       = -79.8 
    h.v_init = -79.8 
    cap         = 1.0
    rall        = 150.0
    rm          = 10e3 
    p_ek          = -85.0 
    p_ena        = 55.0 
    gbar_h      = h_gbar * ihscale 
    gbar_kdmc   = 0.00085
    kdmc_gbar_somam = 20
    sh_nax = 0.0
    gbar_nax    = 0.027 * 3 
    nax_gbar_somam = 5
    gbar_kdr    = 0.007 * 3
    kdr_gbar_somam = 5
    h.a0n_kdr     = 0.0075 
    h.nmax_kdr    = 20.0 
    sh_kap = 0.0
    gbar_kap = 0.1 * 3
    kap_gbar_somam = 5
    h.vhalfn_kap  = 35.0 
    h.nmin_kap    = 0.4 
    h.lmin_kap    = 5.0 
    h.tq_kap      = -45.0 
    km_gmax = 0.1
    cal_gcalbar = cabar*1.2 
    can_gcanbar = cabar*0.3 
    cat_gcatbar = cabar*0.84 
    cal_gbar_somam = can_gbar_somam = cat_gbar_somam = 0.25
    cal_gbar_bdendm = can_gbar_bdendm = cat_gbar_bdendm = 0.25
    ikc_gbar_dendm = 0.25
    for sec in self.all_sec:
      # erev
      sec.ek = p_ek # K+ current reversal potential (mV)
      sec.ena = p_ena # Na+ current reversal potential (mV)
      sec.g_pas = 1.0/rm
      sec.Ra = rall
      sec.cm = cap
      sec.e_pas = Vrest
      # Ih
      sec.eh = erevh
      for seg in sec:
        seg.iar.k2 = 1e-4 # 1e-5 # 1e-4;
        seg.iar.ghbar = gbar_h
      sec.gbar_nax = gbar_nax
      sec.sh_nax = sh_nax
      sec.gbar_kdr = gbar_kdr
      sec.gbar_kap = gbar_kap
      sec.sh_kap = sh_kap
      sec.gSK_E2bar_SK_E2 = cagk_gbar
    soma = self.soma
    soma.gbar_kdmc  = gbar_kdmc * kdmc_gbar_somam
    soma.gbar_nax = gbar_nax * nax_gbar_somam
    soma.gbar_kdr = gbar_kdr * kdr_gbar_somam
    soma.gbar_kap = gbar_kap * kap_gbar_somam
    soma.gkbar_ikc = ikc_gkbar
    soma.gcalbar_cal = cal_gcalbar * cal_gbar_somam
    soma.gcanbar_can = can_gcanbar * can_gbar_somam
    soma.gcatbar_cat = cat_gcatbar * cat_gbar_somam
    soma.gSK_E2bar_SK_E2 = cagk_gbar
    h.distance(0,0.5,sec=self.soma) # middle of soma is origin for distance
    for sec in self.apic:
      sec.gcalbar_cal = cal_gcalbar
      sec.gcanbar_can = can_gcanbar
      sec.gcatbar_cat = cat_gcatbar
      sec.gkbar_ikc = ikc_gkbar * ikc_gbar_dendm
      #sec.gbar_cagk = cagk_gbar
      #sec.gSK_E2bar_SK_E2 = cagk_gbar
      for seg in sec:
        d = h.distance(seg.x,sec=sec)
        seg.iar.ghbar = gbar_h * exp(d/h_lambda)
        if expk:
          seg.gmax_km = km_gmax * exp(d/h_lambda)
          seg.gbar_kap = soma.gbar_kap * exp(d/h_lambda)
          seg.gbar_kdr = soma.gbar_kdr * exp(d/h_lambda)
        else:
          sec.gmax_km = km_gmax # slow voltage-dependent non-inactivating K+
          sec.gbar_kap = soma.gbar_kap;
          sec.gbar_kdr = soma.gbar_kdr
    self.apic[2].cm = 2.0
    Bdend = self.Bdend
    Bdend.gcalbar_cal = cal_gcalbar * cal_gbar_bdendm
    Bdend.gcanbar_can = can_gcanbar * can_gbar_bdendm
    Bdend.gcatbar_cat = cat_gcatbar * cat_gbar_bdendm
    Bdend.gkbar_ikc = ikc_gkbar * ikc_gbar_dendm
    #Bdend.gbar_cagk = cagk_gbar
    #Bdend.gSK_E2bar_SK_E2 = cagk_gbar
    Bdend.gbar_kap = soma.gbar_kap; Bdend.gbar_kdr = soma.gbar_kdr
    Bdend.gmax_km = km_gmax

  def set_conductances (self): # insert the conductances
    for sec in self.all_sec:
      sec.insert('k_ion')
      sec.insert('na_ion')
      sec.insert('ca_ion')
      sec.insert('pas') # passive
      sec.insert('iar') # H channel in Ih.mod
      sec.insert('nax') # Na current
      sec.insert('kdr') # K delayed rectifier current
      sec.insert('kap') # K-A current
      # calcium-related channels
      sec.insert('cal') # cal_mig.mod
      sec.insert('can') # can_mig.mod
      sec.insert('cat') # cat_mig.mod
      sec.insert('ikc') # IC.mod - ca and v dependent k channel
      sec.insert('SK_E2') # SK_E2.mod                           <       
    soma = self.soma; self.soma.insert('kdmc') # K-D current in soma
    for sec in self.apic:
      sec.insert('km') # km.mod
      #sec.insert('cagk') # cagk.mod 
      #sec.insert('SK_E2') # SK_E2.mod 
    self.Bdend.insert('km') # km.mod
    #self.Bdend.insert('cagk') # cagk.mod 
    #self.Bdend.insert('SK_E2') # SK_E2.mod 
		
  def set_synapses(self):
    erevgaba = -80
    self.somaGABAf = Synapse(sect=self.soma,loc=0.5,tau1=0.07,tau2=9.1,e=erevgaba)
    self.somaAMPA = Synapse(sect=self.soma,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.somaNMDA = SynapseNMDA(sect=self.soma,loc=0.5, tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    bdsyloc = 0.5 
    self.BdendAMPA = Synapse(sect=self.Bdend,loc=bdsyloc,tau1=0.05, tau2=5.3,e=0)    
    self.BdendNMDA = SynapseNMDA(sect=self.Bdend,loc=bdsyloc,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend1GABAs = Synapse(sect=self.Adend1,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend2GABAs = Synapse(sect=self.Adend2,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend3GABAs = Synapse(sect=self.Adend3,loc=0.5,tau1=0.2,tau2=20,e=erevgaba)
    self.Adend3GABAf = Synapse(sect=self.Adend3,loc=0.5,tau1=0.07,tau2=9.1,e=erevgaba)
    self.Adend3AMPA = Synapse(sect=self.Adend3,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend3NMDA = SynapseNMDA(sect=self.Adend3,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend2AMPA = Synapse(sect=self.Adend2,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend2NMDA = SynapseNMDA(sect=self.Adend2,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend1AMPA = Synapse(sect=self.Adend1,loc=0.5,tau1=0.05,tau2=5.3,e=0)
    self.Adend1NMDA = SynapseNMDA(sect=self.Adend1,loc=0.5,tau1NMDA=tau1NMDAEE,tau2NMDA=tau2NMDAEE,r=1,e=0)
    self.Adend3mGLUR = SynapsemGLUR(sect=self.Adend3,loc=0.5)
    self.Adend3GABAB = SynapseGABAB(sect=self.Adend3,loc=0.5)
    self.Adend2mGLUR = SynapsemGLUR(sect=self.Adend2,loc=0.5)
    self.Adend2GABAB = SynapseGABAB(sect=self.Adend2,loc=0.5)
    self.Adend1mGLUR = SynapsemGLUR(sect=self.Adend1,loc=0.5)
    self.Adend1GABAB = SynapseGABAB(sect=self.Adend1,loc=0.5)

  # record some of the synaptic currents
  def reccurr (self):
    self.Adend3_iAM = saferecord(self.Adend3AMPA.syn._ref_i, recvdt)
    self.Adend3_iNM = saferecord(self.Adend3NMDA.syn._ref_iNMDA, recvdt)
    self.Adend3_iGB = saferecord(self.Adend3GABAB.syn._ref_i, recvdt)
    self.Adend3_iGA = saferecord(self.Adend3GABAs.syn._ref_i, recvdt)
    self.soma_iGA = saferecord(self.somaGABAf.syn._ref_i, recvdt)
    self.Adend3_ina = saferecord(self.Adend3(0.5)._ref_ina, recvdt)
    self.Adend3_ik = saferecord(self.Adend3(0.5)._ref_ik, recvdt)
    self.Adend3_ica = saferecord(self.Adend3(0.5)._ref_ica, recvdt)
    self.Adend3_ih = saferecord(self.Adend3(0.5)._ref_ih, recvdt)
    self.Adend2_iNM = saferecord(self.Adend2NMDA.syn._ref_iNMDA, recvdt) #Tuomo: added for _saveExtra
    self.Adend2_ina = saferecord(self.Adend2(0.5)._ref_ina, recvdt)
    self.Adend2_ik = saferecord(self.Adend2(0.5)._ref_ik, recvdt)
    self.Adend2_ica = saferecord(self.Adend2(0.5)._ref_ica, recvdt)
    self.Adend2_ih = saferecord(self.Adend2(0.5)._ref_ih, recvdt)
    self.Adend1_iNM = saferecord(self.Adend1NMDA.syn._ref_iNMDA, recvdt) #Tuomo: added for _saveExtra
    self.Adend1_ina = saferecord(self.Adend1(0.5)._ref_ina, recvdt)
    self.Adend1_ik = saferecord(self.Adend1(0.5)._ref_ik, recvdt)
    self.Adend1_ica = saferecord(self.Adend1(0.5)._ref_ica, recvdt)
    self.Adend1_ih = saferecord(self.Adend1(0.5)._ref_ih, recvdt)
    self.Bdend_iNM = saferecord(self.BdendNMDA.syn._ref_iNMDA, recvdt) #Tuomo: added for _saveExtra
    self.Bdend_ina = saferecord(self.Bdend(0.5)._ref_ina, recvdt)
    self.Bdend_ik = saferecord(self.Bdend(0.5)._ref_ik, recvdt)
    self.Bdend_ica = saferecord(self.Bdend(0.5)._ref_ica, recvdt)
    self.Bdend_ih = saferecord(self.Bdend(0.5)._ref_ih, recvdt)
    self.soma_iNM = saferecord(self.somaNMDA.syn._ref_iNMDA, recvdt) #Tuomo: added for _saveExtra
    self.soma_ina = saferecord(self.soma(0.5)._ref_ina, recvdt)
    self.soma_ik = saferecord(self.soma(0.5)._ref_ik, recvdt)
    self.soma_ica = saferecord(self.soma(0.5)._ref_ica, recvdt)
    self.soma_ih = saferecord(self.soma(0.5)._ref_ih, recvdt)

#######################################
#      some utils to avoid the h.     #
vlk = h.vlk
Vector = h.Vector
NQS = h.NQS
gg = h.gg
ge = h.ge
Random = h.Random
List = h.List
Matrix = h.Matrix
nqsdel = h.nqsdel
Graph = h.Graph
vrsz = h.vrsz
allocvecs = h.allocvecs
NetCon = h.NetCon
NetStim = h.NetStim
#######################################


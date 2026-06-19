import configparser
import io

# default config as string
def_config = """
[seed]
iseed = 1234
wseed = 4321
pseed = 4321
[netsyn]
NMAMREE = 0.1
NMAMREI = 0.1
mGLURR = 0.0
GB2R = 7.5
rdmsec = 1
nmfracca = 0.13
[chan]
ihginc = 2.0
ihscale = 1.0
ilscale = 1.0
ikmscale = 1.0
iark2fctr = 1.0
iark4 = 0.008
erevh = -30.0
h_lambda = 325.0
h_gbar = 0.0025
cagk_gbar = 0.0001
ikc_gkbar = 0.003
expk = 1
cabar = 0.005
[cada]
taur = 5
[run]
indir = data
outdir = data
tstop = 20000.0
dt = 0.1
saveout = 0
saveoutSpikes = 1
simstr = 134mar10_mpi_rxd_test_
dorun = 1
doquit = 0
dodraw = 0
verbose = 0
saveExtra = 0
cvodeactive=0
recdt = 10.0
recvdt = 0.1
binsz = 5
loadtstop = 0
saveconns = 0
[rxd]
CB_frate=5.5
CB_brate=0.0026
CB_init=0.2
gip3 = 120400.0
gserca = 4.0
gleak = 3.0
cacytinit = 100e-6
caerinit = 1.25
caexinit = 0.0
spaceum = 0.0
nsubseg = 0
subsegum = 0.0
[net]
scale=1.0
wirety = swire
IIGain = 0.75e-4
IEGain = 0.1
EIGainFS = 0.002
EIGainLTS = 0.002
EEGain = 0.00025
[mem]
baset = 1e3
nqm = None
useWM = 1
startt = 5e3
nMem = 1
nMemInhib = 0
intert = 4e3
durt = 100
rate = 500
weight = 0.0036
frac = 0.5,0.5,0.5,0.5
mGLURRWM = 1000.0
NMAMRWM = 1.0
AMRWM = 1.0
GARWM = 10.0
same = 0
lastmGLURON = 0
lastdurt = 1e3
apicIDX = 0
pops = E2,E5R,E5B,E6
[stim]
EXGain = 1.0
noise = 1
sgrhzNMI = 300.0
sgrhzNME = 300.0
sgrhzEE = 800.0
sgrhzEI = 800.0
sgrhzIE = 150.0
sgrhzII = 150.0
sgrhzMGLURE = 0.0
sgrhzGB2 = 0.0
"""

# write config file starting with defaults and new entries
# specified in section (sec) , option (opt), and value (val)
# saves to output filepath fn
def writeconf (fn,sec,opt,val):
  conf = configparser.ConfigParser()
  conf.optionxform = str
  conf.readfp(io.BytesIO(def_config)) # start with defaults
  # then change entries by user-specs
  for i in range(len(sec)): conf.set(sec[i],opt[i],val[i])
  # write config file
  with open(fn, 'wb') as cfile: conf.write(cfile)

# read config file
def readconf (fn="netcfg.cfg"):
  config = configparser.ConfigParser()
  config.optionxform = str
  config.read(fn)
  def conffloat (base,var,defa): # defa is default value
    val = defa
    try: val=config.getfloat(base,var)
    except: pass
    return val
  def confint (base,var,defa):
    val = defa
    try: val=config.getint(base,var)
    except: pass
    return val
  def confstr (base,var,defa):
    val = defa
    try: val = config.get(base,var)
    except: pass
    return val
  d = {}
  d['iseed'] = confint("seed","iseed",1234)
  d['wseed'] = confint("seed","wseed",4321)
  d['pseed'] = confint("seed","pseed",4321)
  d['NMAMREE'] = conffloat("netsyn","NMAMREE",0.1)
  d['NMAMREI'] = conffloat("netsyn","NMAMREI",0.1)
  d['mGLURR'] = conffloat("netsyn","mGLURR",0.0)
  d['GB2R'] = conffloat("netsyn","GB2R",7.5)
  d['nmfracca'] = conffloat("netsyn","nmfracca", 0.13)
  d['rdmsec'] = confint("netsyn","rdmsec", 1)
  d['erevh'] = conffloat("chan","erevh",-30.0)
  d['h_lambda'] = conffloat("chan","h_lambda",325.0)
  d['h_gbar'] = conffloat("chan","h_gbar",0.006)
  d['cagk_gbar'] = conffloat("chan","cagk_gbar",0.0001)
  d['ikc_gkbar'] = conffloat("chan","ikc_gkbar",0.003)
  d['ihscale'] = conffloat("chan","ihscale", 1.0)
  d['ihginc'] = conffloat("chan","ihginc", 2.0)
  d['ilscale'] = conffloat("chan","ilscale", 1.0)
  d['ikmscale'] = conffloat("chan", "ikmscale", 1.0)
  d['iark2fctr'] = conffloat("chan", "iark2fctr",1.0)
  d['iark4'] = conffloat("chan", "iark4",0.008)
  d['expk'] = confint("chan","expk",1)
  d['cabar'] = conffloat("chan","cabar",0.005)
  d['taurcada'] = conffloat("cada", "taur", 5.0)
  d['outdir'] = confstr("run","outdir", "data")
  d['indir'] = confstr("run","indir", "data")
  d['tstop'] = conffloat("run","tstop", 20000.0)
  d['dt'] = conffloat("run","dt",0.1)
  d['saveout'] = conffloat("run","saveout",1)
  d['saveoutSpikes'] = conffloat("run","saveoutSpikes",1)
  d['simstr'] = confstr("run","simstr","13nov5_mpi_rxd_test_")
  d['dorun'] = confint("run","dorun",1)
  d['recdt'] = conffloat("run","recdt",10.0)
  d['recvdt'] = conffloat("run","recvdt",1.0)
  d['binsz'] = conffloat("run","binsz",5)
  for k in ['doquit','verbose','saveExtra','cvodeactive','dodraw']:d[k] = confint("run",k,0)
  d['CB_frate'] = conffloat("rxd","CB_frate", 5.5)
  d['CB_brate'] = conffloat("rxd","CB_brate", 0.0026)
  d['CB_init'] = conffloat("rxd","CB_init", 0.2)
  d['gip3'] = conffloat("rxd","gip3",120400.0)
  d['gserca'] = conffloat("rxd","gserca",4.0)
  d['gleak'] = conffloat("rxd","gleak",3.0)
  d['caerinit'] = conffloat("rxd","caerinit",1.25)
  d['cacytinit'] = conffloat("rxd","cacytinit",100e-6)
  d['caexinit'] = conffloat("rxd","caexinit",0.0)
  d['spaceum'] = conffloat("rxd","spaceum",0.0)
  d['nsubseg'] = confint("rxd","nsubseg",0)
  d['subsegum'] = conffloat("rxd","subsegum",0.0)
  d['scale'] = conffloat("net","scale",1.0)
  d['wirety'] = confstr("net","wirety","swire")
  d['IIGain'] = conffloat("net","IIGain",0.75e-4)
  d['IEGain'] = conffloat("net","IEGain",0.015)
  d['EIGainFS'] = conffloat("net","EIGainFS",0.002)
  d['EIGainLTS'] = conffloat("net","EIGainLTS",0.002)
  d['EEGain'] = conffloat("net","EEGain",0.00025)
  d['baset'] = conffloat("mem","baset",1e3)
  d['nqm'] = confstr("mem","nqm","None")
  d['useWM'] = confint("mem","useWM",1)
  d['nMem'] = confint("mem","nMem",1)
  d['nMemInhib'] = confint("mem","nMemInhib",0)
  d['memintert'] = conffloat("mem","intert",4e3)
  d['memdurt'] = conffloat("mem","durt",100)
  d['memrate'] = conffloat("mem","rate",500)
  d['memW'] = conffloat("mem","weight",0.0036)
  d['memstartt'] = conffloat("mem","startt",5e3)
  d['memfrac'] = confstr("mem","frac","0.5,0.5,0.5,0.5")
  d['mGLURRWM'] = conffloat("mem","mGLURRWM",1000.0)
  d['NMAMRWM'] = conffloat("mem","NMAMRWM",1.0)
  d['AMRWM'] = conffloat("mem","AMRWM",1.0)
  d['GARWM'] = conffloat("mem","GARWM",10.0)
  d['memSame'] = confint("mem","same",0)
  d['lastmGLURON'] = confint("mem","lastmGLURON",0) # whether last stim is only for mGLUR (rest of mGLUR off)
  d['lastmemdurt'] = conffloat("mem","lastdurt",100)
  d['apicIDX'] = confstr("mem","apicIDX","0")
  d['pops'] = confstr("mem","pops","E2,E5R,E5B,E6")
  d['EXGain'] = conffloat("stim","EXGain",1.0)
  d['noise'] = conffloat("stim","noise",1.0)
  d['sgrhzNME'] = conffloat("stim","sgrhzNME",300.0)
  d['sgrhzNMI'] = conffloat("stim","sgrhzNMI",300.0)
  d['sgrhzEE'] = conffloat("stim","sgrhzEE",800.0)
  d['sgrhzIE'] = conffloat("stim","sgrhzIE",150.0)
  d['sgrhzEI'] = conffloat("stim","sgrhzEI",800.0)
  d['sgrhzII'] = conffloat("stim","sgrhzII",150.0)
  d['sgrhzMGLURE'] = conffloat("stim","sgrhzMGLURE",0.0)
  d['sgrhzGB2'] = conffloat("stim","sgrhzGB2",0.0)
  if config.has_option("net","wnq"): d['wnq'] = config.get("net","wnq")
  return d

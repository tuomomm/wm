
IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

#Load the data for CACNA1I expression in ACC
source ./CM_wmNew_samps_IC_ACC_PConly_isSign0_isImputed1_Both_CACNA1I.sh

for isim in `seq 0 480` 
do
 for myseed in 1 
 do
  EXT=_long_withSK3_ACCmildsubj${isim}_altgain_PConly_CACNA1Ionly
  if [ -f subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat ]
  then
   echo "FRandMUA_sim_net${EXT}_seed${myseed}.mat exists"
   continue
  fi
  echo "FRandMUA_sim_net${EXT}_seed${myseed}.mat does not exist"
  cp netcfg.cfg netcfg${EXT}.cfg

  cat netcfg.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/"| sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/"> netcfg${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" > geom${EXT}.py
  echo """import numpy""" > mpisim${EXT}.py
  cat mpisim_withSK3_l5hzonly.py | sed "s/from geom_withSK3/from geom${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" >> mpisim${EXT}.py
  echo """from pylab import *""" > simdat${EXT}.py
  cat simdat_saveFRandMUAonly_allowNonArray.py | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" >> simdat${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${BLOCKEDSUBJS[isim]} simstr=sim_net${EXT}"
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${BLOCKEDSUBJS[isim]} simstr=sim_net${EXT}

  echo "python3 simdat${EXT}.py"
  python3 simdat${EXT}.py

  if [ -f FRandMUA_sim_net${EXT}.mat ]
  then
    echo "mv FRandMUA_sim_net${EXT}.mat subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat"
    mv FRandMUA_sim_net${EXT}.mat subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat
  elif [ -f FR_sim_net${EXT}.mat ]
  then
    echo "mv FR_sim_net${EXT}.mat subjs/FR_sim_net${EXT}_seed${myseed}.mat"
    mv FR_sim_net${EXT}.mat subjs/FR_sim_net${EXT}_seed${myseed}.mat
  else
    echo "FR files FRandMUA_sim_net${EXT}.mat or FR_sim_net${EXT}.mat not found"
  fi

  echo "rm data/sim_net${EXT}_pc_?.npz" 
  rm data/sim_net${EXT}_pc_?.npz #Delete raw data files to save space

  echo "rm netcfg${EXT}.cfg geom${EXT}.py mpisim${EXT}.py simdat${EXT}.py"
  rm netcfg${EXT}.cfg geom${EXT}.py mpisim${EXT}.py simdat${EXT}.py #Delete additional scripts that were made here
 done 
done

#Load the data for CACNA1I expression in PFC
source ./CM_wmNew_samps_IC_PFC_PConly_isSign0_isImputed1_Both_CACNA1I.sh

for isim in `seq 0 430` 
do
 for myseed in 1 
 do
  EXT=_long_withSK3_PFCmildsubj${isim}_altgain_PConly_CACNA1Ionly
  if [ -f subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat ]
  then
   echo "FRandMUA_sim_net${EXT}_seed${myseed}.mat exists"
   continue
  fi
  echo "FRandMUA_sim_net${EXT}_seed${myseed}.mat does not exist"
  cp netcfg.cfg netcfg${EXT}.cfg

  cat netcfg.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/"| sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/"> netcfg${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" > geom${EXT}.py
  echo """import numpy""" > mpisim${EXT}.py
  cat mpisim_withSK3_l5hzonly.py | sed "s/from geom_withSK3/from geom${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" >> mpisim${EXT}.py
  echo """from pylab import *""" > simdat${EXT}.py
  cat simdat_saveFRandMUAonly_allowNonArray.py | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" >> simdat${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${BLOCKEDSUBJS[isim]} simstr=sim_net${EXT}"
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${BLOCKEDSUBJS[isim]} simstr=sim_net${EXT}

  echo "python3 simdat${EXT}.py"
  python3 simdat${EXT}.py

  if [ -f FRandMUA_sim_net${EXT}.mat ]
  then
    echo "mv FRandMUA_sim_net${EXT}.mat subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat"
    mv FRandMUA_sim_net${EXT}.mat subjs/FRandMUA_sim_net${EXT}_seed${myseed}.mat
  elif [ -f FR_sim_net${EXT}.mat ]
  then
    echo "mv FR_sim_net${EXT}.mat subjs/FR_sim_net${EXT}_seed${myseed}.mat"
    mv FR_sim_net${EXT}.mat subjs/FR_sim_net${EXT}_seed${myseed}.mat
  else
    echo "FR files FRandMUA_sim_net${EXT}.mat or FR_sim_net${EXT}.mat not found"
  fi

  echo "rm data/sim_net${EXT}_pc_?.npz" 
  rm data/sim_net${EXT}_pc_?.npz #Delete raw data files to save space

  echo "rm netcfg${EXT}.cfg geom${EXT}.py mpisim${EXT}.py simdat${EXT}.py"
  rm netcfg${EXT}.cfg geom${EXT}.py mpisim${EXT}.py simdat${EXT}.py #Delete additional scripts that were made here
 done 
done

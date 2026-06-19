
IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

#Each simulation takes 1 minute or so
for myseed in 1 2 3 4 5 6 7 8 9 10 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
do
    echo "python3 onepyr_noplot_withSK3_withIh.py tstop=100000 iseed=${myseed} wseed=${myseed} pseed=${myseed} simstr=sim_onepyr_withSK3_seed${myseed} IEGain=${IEGain} EIGainFS=${EIGainFS} EIGainLTS=${EIGainLTS} EEGain=${EEGain}"
    python3 onepyr_noplot_withSK3_withIh.py tstop=100000 iseed=${myseed} wseed=${myseed} pseed=${myseed} simstr=sim_onepyr_withSK3_seed${myseed} IEGain=${IEGain} EIGainFS=${EIGainFS} EIGainLTS=${EIGainLTS} EEGain=${EEGain}
done


#Each DC simulation takes a few seconds
for EXT in withSK3_longDC_apamin withSK3_longDC_fullapamin withSK3_longDC_95apamin withSK3_longDC
do    
  for amp in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
  do
    echo "python3 onepyr_noplot_${EXT}.py $amp"
    python3 onepyr_noplot_${EXT}.py $amp
  done

done
echo "python3 onepyr_noplot_withSK3_VClamp.py" #Run the VClamp experiments
python3 onepyr_noplot_withSK3_VClamp.py

#Network simulations: each simulation takes a couple of hours when parallellized to 8 cores
for myseed in 1
do
    
  EXT=_long_withSK_altgain_savespikes_seed${myseed}
  cp netcfg.cfg netcfg${EXT}.cfg

  cat netcfg.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/" | sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/" > netcfg${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" > geom${EXT}.py
  echo """import numpy""" > mpisim${EXT}.py
  cat mpisim_withSK3_l5hzonly.py | sed "s/from geom_withSK3/from geom${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" >> mpisim${EXT}.py
  echo """from pylab import *""" > simdat${EXT}.py
  cat simdat_savespikes.py | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" >> simdat${EXT}.py

  echo "mpiexec -np 8 nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 simstr=sim_net${EXT}"
  mpiexec -np 8 nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 simstr=sim_net${EXT}

  echo "python3 simdat${EXT}.py"
  python3 simdat${EXT}.py

  #echo "rm data/sim_net${EXT}_pc_?.npz"
  #rm data/sim_net${EXT}_pc_?.npz #Don't delete these files because they are, exceptionally, used for plotting
  
done

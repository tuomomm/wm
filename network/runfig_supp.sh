#Simulations of a single PC where no hot zone is implemented. Gives results similar to those with the hot zone.
for amp in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
    echo "python3 onepyr_noplot_withSK3_longDC_95apamin_nohotzone.py $amp" 
    python3 onepyr_noplot_withSK3_longDC_95apamin_nohotzone.py $amp 

    echo "python3 onepyr_noplot_withSK3_longDC_apamin_nohotzone.py $amp" 
    python3 onepyr_noplot_withSK3_longDC_apamin_nohotzone.py $amp 

    echo "python3 onepyr_noplot_withSK3_longDC_fullapamin_nohotzone.py $amp" 
    python3 onepyr_noplot_withSK3_longDC_fullapamin_nohotzone.py $amp 
done
echo "python3 onepyr_noplot_withSK3_nohotzone_VClamp.py"
python3 onepyr_noplot_withSK3_nohotzone_VClamp.py


#Simulations where the SCZ-ACC and SCZ-PFC parameter changes apply to interneurons as well
EXTS=("_PFCmild" "_ACCmild")

PARCHANGES=("gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gkdr_Kdrbwb=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gna_Nafbwb=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gca_icalts=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gh_ihlts=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gkdr_Kdrbwb=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gca_icalts=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912")

IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

for myseed in 1 2 3 4 5 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
do
    
 for isim in 0 1
 do
  EXT=_long_withSK3${EXTS[isim]}_altgain
  cp netcfg.cfg netcfg${EXT}.cfg

  cat netcfg.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/" | sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/" > netcfg${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" > geom${EXT}.py
echo """import numpy""" > mpisim${EXT}.py
  cat mpisim_withSK3_l5hzonly.py | sed "s/from geom_withSK3/from geom${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" >> mpisim${EXT}.py
  echo """from pylab import *""" > simdat${EXT}.py
  cat simdat_saveFRandMUAonly.py | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" >> simdat${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${PARCHANGES[isim]} simstr=sim_net${EXT}"
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${PARCHANGES[isim]} simstr=sim_net${EXT}

  echo "python3 simdat${EXT}.py"
  python3 simdat${EXT}.py

  echo "mv FRandMUA_sim_net${EXT}.mat FRandMUA_sim_net${EXT}_seed${myseed}.mat"
  mv FRandMUA_sim_net${EXT}.mat FRandMUA_sim_net${EXT}_seed${myseed}.mat
  
  echo "rm data/sim_net${EXT}_pc_?.npz"
  rm data/sim_net${EXT}_pc_?.npz
  
 done 
done


#Simulations where the HCN-mediated response to WM-inducing stimulus is replaced by external AMPAR-mediated increased synaptic inputs. Only CTRL with various kCoeffs.
EXTS=("_CTRL")

PARCHANGES=("")

IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

isim=0
for kCoeff in 0 0.05 0.1 0.15 0.2 0.3
do
    
 kCoeffStr="${kCoeff/./_}"

 for myseed in 1 2 3 4 5 6 7 8 9 10 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
 do
  EXT=_long_withSK3${EXTS[isim]}_altgain_k${kCoeffStr}
    
  cat netcfg_ExpWM.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net_ExpWMnoNM${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/"| sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/" | sed "s/kExpWM = 3.3/kExpWM = $kCoeff/" > netcfg_ExpWMnoNM${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg_ExpWMnoNM${EXT}.cfg/" > geom_ExpWMnoNM${EXT}.py

  echo """import numpy""" > mpisim_ExpWMnoNM${EXT}.py
  cat mpisim_long_withSK3_CTRL_altgain_ExpWMnoNM_verbose.py | sed "s/from geom_long_withSK3_CTRL_altgain/from geom_ExpWMnoNM${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim_long_withSK3_CTRL_altgain.py/mpisim_ExpWMnoNM${EXT}.py/" | sed "s/geom.py/geom_ExpWMnoNM${EXT}.py/" | sed "s/15dec31.1/sim_net_ExpWMnoNM${EXT}/" >> mpisim_ExpWMnoNM${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim_ExpWMnoNM${EXT}.py netcfg_ExpWMnoNM${EXT}.cfg tstop=100000 kExpWM=$kCoeff ${PARCHANGES[isim]} simstr=sim_net_ExpWMnoNM${EXT}" #The parameter change kExpWM=$kCoeff has to be given at the command
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim_ExpWMnoNM${EXT}.py netcfg_ExpWMnoNM${EXT}.cfg tstop=100000 kExpWM=$kCoeff ${PARCHANGES[isim]} simstr=sim_net_ExpWMnoNM${EXT}        #line, conf.py doesn't extract it from the cfg file

  cat simdat_ExpWM.py | sed "s/mpisim_ExpWM.py/mpisim_ExpWMnoNM${EXT}.py/" | sed "s/sim_net_ExpWM/sim_net_ExpWMnoNM${EXT}/" | sed "s/netcfg_ExpWM.cfg/netcfg_ExpWMnoNM${EXT}.cfg/" | sed "s/geom_ExpWM.py/geom_ExpWMnoNM${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" > simdat_ExpWMnoNM${EXT}.py

  echo "python3 simdat_ExpWMnoNM${EXT}.py"
  python3 simdat_ExpWMnoNM${EXT}.py

  echo "mv FRandMUA_sim_net_ExpWMnoNM${EXT}.mat FRandMUA_sim_net_ExpWMnoNM${EXT}_seed${myseed}.mat"
  mv FRandMUA_sim_net_ExpWMnoNM${EXT}.mat FRandMUA_sim_net_ExpWMnoNM${EXT}_seed${myseed}.mat
 done
done


#Simulations where the HCN-mediated response to WM-inducing stimulus is replaced by external AMPAR-mediated increased synaptic inputs. Both CTRL and SCZ-ACC and SCZ-PFC with 0.1 kCoeff.
EXTS=("_CTRL" "_PFCmild" "_ACCmild")

PARCHANGES=("" "gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gkdr_Kdrbwb=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gna_Nafbwb=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gca_icalts=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gh_ihlts=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gkdr_Kdrbwb=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gca_icalts=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912")

IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

for isim in 0 1 2
do
    
 kCoeff=0.1
 kCoeffStr="${kCoeff/./_}"

 for myseed in 1 2 3 4 5 6 7 8 9 10 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
 do
  EXT=_long_withSK3${EXTS[isim]}_altgain_k${kCoeffStr}
    
  cat netcfg_ExpWM.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net_ExpWMnoNM${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/"| sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/" | sed "s/kExpWM = 3.3/kExpWM = $kCoeff/" > netcfg_ExpWMnoNM${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg_ExpWMnoNM${EXT}.cfg/" > geom_ExpWMnoNM${EXT}.py

  echo """import numpy""" > mpisim_ExpWMnoNM${EXT}.py
  cat mpisim_long_withSK3_CTRL_altgain_ExpWMnoNM_verbose.py | sed "s/from geom_long_withSK3_CTRL_altgain/from geom_ExpWMnoNM${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim_long_withSK3_CTRL_altgain.py/mpisim_ExpWMnoNM${EXT}.py/" | sed "s/geom.py/geom_ExpWMnoNM${EXT}.py/" | sed "s/15dec31.1/sim_net_ExpWMnoNM${EXT}/" >> mpisim_ExpWMnoNM${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim_ExpWMnoNM${EXT}.py netcfg_ExpWMnoNM${EXT}.cfg tstop=100000 kExpWM=$kCoeff ${PARCHANGES[isim]} simstr=sim_net_ExpWMnoNM${EXT}" #The parameter change kExpWM=$kCoeff has to be given at the command
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim_ExpWMnoNM${EXT}.py netcfg_ExpWMnoNM${EXT}.cfg tstop=100000 kExpWM=$kCoeff ${PARCHANGES[isim]} simstr=sim_net_ExpWMnoNM${EXT}        #line, conf.py doesn't extract it from the cfg file

  cat simdat_ExpWM.py | sed "s/mpisim_ExpWM.py/mpisim_ExpWMnoNM${EXT}.py/" | sed "s/sim_net_ExpWM/sim_net_ExpWMnoNM${EXT}/" | sed "s/netcfg_ExpWM.cfg/netcfg_ExpWMnoNM${EXT}.cfg/" | sed "s/geom_ExpWM.py/geom_ExpWMnoNM${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" > simdat_ExpWMnoNM${EXT}.py

  echo "python3 simdat_ExpWMnoNM${EXT}.py"
  python3 simdat_ExpWMnoNM${EXT}.py

  echo "mv FRandMUA_sim_net_ExpWMnoNM${EXT}.mat FRandMUA_sim_net_ExpWMnoNM${EXT}_seed${myseed}.mat"
  mv FRandMUA_sim_net_ExpWMnoNM${EXT}.mat FRandMUA_sim_net_ExpWMnoNM${EXT}_seed${myseed}.mat

 done
done

#Simulations where EE gain is adjusted according to the baseline synaptic conductance data from the simulations of CTRL, SCZ-ACC and SCZ-PFC (expression of plasticity-regulating genes)

EXTS=("_PFCmild_adjEEB" "_ACCmild_adjEEB"
      "_PFCmild_PConly_adjEEB" "_ACCmild_PConly_adjEEB"
      "_PFCmild_nax_PConly_adjEEB" "_ACCmild_kap_adjEEB"
      "_PFCmild_nax_all_adjEEB" )

PARCHANGES=("gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gkdr_Kdrbwb=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gna_Nafbwb=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gca_icalts=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gh_ihlts=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gkdr_Kdrbwb=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gca_icalts=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912"
          "gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912"
            "gCoeff_gbar_nax=0.9658" "gCoeff_gbar_kap=1.0236"
            "gCoeff_gbar_nax=0.9658 gCoeff_gna_Nafbwb=0.9658" )

EEGAINS=(0.00050548 0.00052416 0.00050548 0.00052416 0.00050548 0.00052416 0.00050548)


for myseed in 1 2 3 4 5 6 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
do
    
 for isim in 0 #This will be replaced in _submitnow.sh
 do
  EXT=_long_withSK3${EXTS[isim]}_altgain
  cp netcfg.cfg netcfg${EXT}.cfg

  IEGain=0.05
  EIGainFS=0.008
  EIGainLTS=0.004
  EEGain=${EEGAINS[isim]}

  cat netcfg.cfg |  sed "s/saveout = 0/saveout = 1/" | sed "s/simstr = 15dec31_1/simstr = sim_net${EXT}/" |  sed "s/recvdt = 0.1/recvdt = 10/" | sed "s/tstop = 14000/tstop = 100000/" | sed "s/iseed = 1234/iseed = ${myseed}/" | sed "s/seed = 4321/seed = ${myseed}/" | sed "s/IEGain = 0.1/IEGain = ${IEGain}/" | sed "s/EIGainFS = 0.002/EIGainFS = ${EIGainFS}/" | sed "s/EIGainLTS = 0.002/EIGainLTS = ${EIGainLTS}/" | sed "s/EEGain = 0.00025/EEGain = ${EEGain}/" > netcfg${EXT}.cfg

  cat geom_withSK3.py | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" > geom${EXT}.py
echo """import numpy""" > mpisim${EXT}.py
  cat mpisim_withSK3_l5hzonly.py | sed "s/from geom_withSK3/from geom${EXT}/" | sed "s/if os.path.exists(dconf/if False: #os.path.exists(dconf/" | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" >> mpisim${EXT}.py
  echo """from pylab import *""" > simdat${EXT}.py
  cat simdat_saveFRandMUAonly.py | sed "s/mpisim.py/mpisim${EXT}.py/" | sed "s/15dec31.1/sim_net${EXT}/" | sed "s/netcfg.cfg/netcfg${EXT}.cfg/" | sed "s/geom.py/geom${EXT}.py/" | sed "s/8,14/8,100/" | sed "s/lw=0.5/lw=0.25/" | sed "s/lwidth=1/lwidth=0.25/" | sed "s/lwidth=4/lwidth=1/" >> simdat${EXT}.py

  echo "mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${PARCHANGES[isim]} simstr=sim_net${EXT}"
  mpiexec -np 8 --bind-to core:overload-allowed nrniv -python -mpi mpisim${EXT}.py netcfg${EXT}.cfg tstop=100000 ${PARCHANGES[isim]} simstr=sim_net${EXT}

  echo "python3 simdat${EXT}.py"
  python3 simdat${EXT}.py

  echo "mv FRandMUA_sim_net${EXT}.mat FRandMUA_sim_net${EXT}_seed${myseed}.mat"
  mv FRandMUA_sim_net${EXT}.mat FRandMUA_sim_net${EXT}_seed${myseed}.mat
  
  echo "rm data/sim_net${EXT}_pc_?.npz"
  rm data/sim_net${EXT}_pc_?.npz
  
 done 
done

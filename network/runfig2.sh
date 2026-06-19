
EXTS=("_PFC_GABAB" "_PFC_kdr_PConly" "_PFCmild_nax_PConly" "_PFC_cal_PConly" "_PFC_cat" "_PFC_iar_PConly" "_PFC_km" "_PFC_serca"
      "_ACC_GABAB" "_ACC_cagk_ikc" "_ACCmild_kap" "_ACC_kdr_PConly" "_ACC_cal_PConly" "_ACC_cat" "_ACC_km" "_ACC_serca"
      "_CTRL"
      "_PFCmild_PConly" "_ACCmild_PConly")


PARCHANGES=("gCoeff_gbar_GABAB=1.1656" "gCoeff_gbar_kdr=1.0756" "gCoeff_gbar_nax=0.9658" "gCoeff_gcalbar_cal=1.0269" "gCoeff_gcatbar_cat=1.0286" "gCoeff_ghbar_iar=1.0607" "gCoeff_gmax_km=1.0368" "gCoeff_serca=1.0790"
            "gCoeff_gbar_GABAB=1.1498" "gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713" "gCoeff_gbar_kap=1.0236" "gCoeff_gbar_kdr=1.1020" "gCoeff_gcalbar_cal=1.1158" "gCoeff_gcatbar_cat=1.1120" "gCoeff_gmax_km=1.0389" "gCoeff_serca=1.0912"
            ""
            "gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912")


IEGain=0.05
EIGainFS=0.008
EIGainLTS=0.004
EEGain=0.0005

for myseed in 1 #2 3 4 5 6 7 8 9 10 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
do
    
 for isim in `seq 0 18` 
 do
     EXT=_long_withSK3${EXTS[isim]}_altgain
     if [ -f FRandMUA_sim_net${EXT}_seed${myseed}.mat ]
     then
	 echo FRandMUA_sim_net${EXT}_seed${myseed}.mat exists
     else
	 echo FRandMUA_sim_net${EXT}_seed${myseed}.mat does not exist
     
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
     fi
 done 
done


#drawfig1.py: channels = ['CTRL', 'nax_PConly', 'iar_PConly', 'kap', 'cat', 'kdr_PConly', 'serca', 'cagk_ikc', 'km', 'cal_PConly', 'GABAB', 'PConly']

EXTS=("_PFCmild_PConly" "_ACCmild_PConly" "_CTRL"
      "_PFC_cal_PConly" "_PFC_iar_PConly" "_PFC_cat" "_PFC_kdr_PConly" "_PFC_km" "_PFC_serca" "_PFCmild_nax_PConly" "_PFC_GABAB" "_ACC_cal_PConly" "_ACC_cat" "_ACC_kdr_PConly" "_ACC_km" "_ACC_serca" "_ACCmild_kap" "_ACC_cagk_ikc" "_ACC_GABAB"
     )

PARCHANGES=("gCoeff_gbar_GABAB=1.1656 gCoeff_gbar_kdr=1.0756 gCoeff_gbar_nax=0.9658 gCoeff_gcalbar_cal=1.0269 gCoeff_gcatbar_cat=1.0286 gCoeff_ghbar_iar=1.0607 gCoeff_gmax_km=1.0368 gCoeff_serca=1.0790" "gCoeff_gbar_GABAB=1.1498 gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713 gCoeff_gbar_kap=1.0236 gCoeff_gbar_kdr=1.1020 gCoeff_gcalbar_cal=1.1158 gCoeff_gcatbar_cat=1.1120 gCoeff_gmax_km=1.0389 gCoeff_serca=1.0912" ""
            "gCoeff_gcalbar_cal=1.0269" "gCoeff_ghbar_iar=1.0607" "gCoeff_gcatbar_cat=1.0286" "gCoeff_gbar_kdr=1.0756" "gCoeff_gmax_km=1.0368" "gCoeff_serca=1.0790" "gCoeff_gbar_nax=0.9658" "gCoeff_gbar_GABAB=1.1656" "gCoeff_gcalbar_cal=1.1158" "gCoeff_gcatbar_cat=1.1120" "gCoeff_gbar_kdr=1.1020" "gCoeff_gmax_km=1.0389" "gCoeff_serca=1.0912" "gCoeff_gbar_kap=1.0236" "gCoeff_gbar_cagk=1.0713 gCoeff_gkbar_ikc=1.0713" "gCoeff_gbar_GABAB=1.1498")


for myseed in 1 2 3 4 5 #The seeds have to be done in the same thread, otherwise they will confuse each others' I/O
do
    
 for isim in `seq 0 18`
 do
  EXT=_long_withSK3${EXTS[isim]}
  if [ -f sim_onepyr${EXT}_altgain_seed1.mat ]
  then
      echo "sim_onepyr${EXT}_altgain_seed1.mat exists"
  else
      echo "sim_onepyr${EXT}_altgain_seed1.mat does not exist"
      echo "python3 onepyr_noplot_withSK3.py tstop=100000 iseed=${myseed} wseed=${myseed} pseed=${myseed} ${PARCHANGES[isim]} simstr=sim_onepyr${EXT}_altgain_seed${myseed}"
      python3 onepyr_noplot_withSK3.py tstop=100000 iseed=${myseed} wseed=${myseed} pseed=${myseed} ${PARCHANGES[isim]} simstr=sim_onepyr${EXT}_altgain_seed${myseed}
  fi
 done
done


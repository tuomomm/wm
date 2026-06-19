TASK_ID=0

TSHORT=27000000
ONSET=24040000

FREQS=(0.001 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0)
NSTIMS=(0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400)

CAFLUX=120.0 #100.0 #125.0 #150.0
LFLUX=5.0
GLUFLUX=10.0
ACHFLUX=10.0
initfile=None

source ./CMcomb_samps_nonIC_ACC.sh #This file contains the relative gene-expression values for each subject. Use the CommonMind data to generate this sensitive file.

GLUR=GluR1,GluR1_memb,GluR2,GluR2_memb
GLURCOEFF=0.5,0.5,1.5,1.5
ifreq=16
FREQ=${FREQS[ifreq]}
NSTIM=${NSTIMS[ifreq]}
THISBLOCKED=${GLUR},${BLOCKED}

for isamp in `seq 0 480`
do
  THISBLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS[isamp]}
  if [ -f nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_ACCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat ]
  then
  echo "nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_ACCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat exists"
  else
  echo "nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_ACCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat does not exist"
  echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $THISBLOCKED $THISBLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_ACCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}"
  python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $THISBLOCKED $THISBLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_ACCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}
  fi

done

for isamp in `seq 0 429`
do
  THISBLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS[isamp]}
  if [ -f nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_PFCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat ]
  then
  echo "nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_PFCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat exists"
  else
  echo "nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_PFCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat does not exist"
  echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $THISBLOCKED $THISBLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_PFCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}"
  python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $THISBLOCKED $THISBLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_PFCx${GLURCOEFF},${isamp}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}
  fi

done

#The list of varied parameters and their coefficients. In the end the combination variants.
BLOCKEDSEPS=('CK' 'Calbin,CalbinC' 'DAGK' 'Gi' 'Gqabg' 'Ng' 'PDE4' 'PKA' 'PKC' 'DAGK' 'Gi' 'NCX' 'Ng' 'PDE4' 'PKA' 'PKC' 'PLA2' 'PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,Ng,PDE4,PKA,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1')
BLOCKEDSEPCOEFFS=(0.959 0.985,0.985 0.998 0.971 1.063 0.962 1.012 1.048 0.996 1.039 0.984 0.999 0.944 0.977 1.039 1.037 0.933 0.966 0.959,0.985,0.985,0.998,0.971,1.063,0.962,1.012,1.048,0.996 1.039,0.984,0.999,0.944,0.977,1.039,1.037,0.933,0.966)

TSHORT=27000000
ONSET=24040000

FREQS=(0.001 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0)
NSTIMS=(0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000)

CAFLUX=120.0
LFLUX=5.0
GLUFLUX=10.0
ACHFLUX=10.0
initfile=None

GLUR=GluR1,GluR1_memb,GluR2,GluR2_memb
GLURCOEFF=0.5,0.5,1.5,1.5

for iblock in `seq 0 19`
do
    
  BLOCKED=${GLUR},${BLOCKEDSEPS[iblock]}
  BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDSEPCOEFFS[iblock]}

  ifreq=16
  FREQ=${FREQS[ifreq]}
  NSTIM=${NSTIMS[ifreq]}
  if [ -f nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat ]
  then
  echo "nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat exists"
  else
  echo "nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat does not exist"
  echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}.mat"
  python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}
  fi
done

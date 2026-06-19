#This script is fully dependent on the file structure at the TSD. Use it only to see how PRSice was used

#module load plink/1.90b6.2
module load PLINK/1.90b_6.2-x86_64
setNames=(genes_all_modeled_synaptic genes_full genes_all_modeled_ionchannels2 genes_all_modeled_ionchannels2_CACNA1C genes_all_modeled_ionchannels2_CACNA1D genes_all_modeled_ionchannels2_CACNA1I genes_all_modeled_ionchannels2_KCNB1 genes_all_modeled_ionchannels2_KCNJ6 genes_all_modeled_ionchannels2_KCND3 genes_all_modeled_ionchannels2_HCN1 genes_all_modeled_ionchannels2_ATP2A2 genes_new genes_new_ionchannels genes_new_synaptic genes_new_ionchannels_kcn genes_new_ionchannels_cacn genes_new_ionchannels_scn_and_hcn genes_new_synaptic_others genes_new_synaptic_PKA genes_new_synaptic_PKC genes_new_synaptic_PKAPKC genes_new_synaptic_PKAPKC_noPP)

pwd
for myset in full ${setNames[@]}
do
for THR in 5e-10 0.0005 5e-05 5e-08 5e-09 1e-10 1e-09 1e-08 1e-07 5e-07 1e-06 5e-06 1e-05 0.0001 
do
  if [ -f genetic/${myset}_${THR}_PRS.best ]
  then
    echo "genetic/${myset}_${THR}_PRS.best already exists"
    continue
  fi
    
  cd i/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out genetic/${myset}_${THR}_i --keep subj_IDs_of_interest.txt --recode
  cd ../../ii/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out genetic/${myset}_${THR}_ii --keep subj_IDs_of_interest.txt --recode
  cd ../../iii/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out genetic/${myset}_${THR}_iii --keep subj_IDs_of_interest.txt --recode

  cd genetic
  for q in ${myset}_${THR}_i ${myset}_${THR}_ii ${myset}_${THR}_iii
  do
      echo "plink --file ${q} --make-bed --out ${q}_new"
      plink --file ${q} --make-bed --out ${q}_new
  done
  
  
  for q in ${myset}_${THR}
  do
    ls ${q}_i*_new*bed |sed 's/.bed//' > mergelist${myset}.txt
    echo "plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged"
    plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged
  done

  PRSice_linux --base PGC_SCZ_0518_EUR_noTOP.sumstats --beta --target ${myset}_${THR}_merged --pheno DUMMY_PHENO.txt --stat BETA --pvalue PVAL --print-snp --no-full --snp VARIANT_ID --bar-levels "1" --fastscore --out ${myset}_${THR}_PRS --no-clump --maf 0.01 --thread 8
  PRSice_linux --base PGC_SCZ_0518_EUR_noTOP.sumstats --beta --target ${myset}_${THR}_merged --pheno DUMMY_PHENO.txt --stat BETA --pvalue PVAL --print-snp --no-full --snp VARIANT_ID --bar-levels "1" --fastscore --out ${myset}_${THR}_PRS --no-clump --maf 0.01 --thread 8 --extract ${myset}_${THR}_PRS.valid
  echo "rm ${myset}_${THR}_PRS.valid"
  rm ${myset}_${THR}_PRS.valid

done

done
setNames=(genes_all_modeled_synaptic genes_full genes_all_modeled_ionchannels2 genes_all_modeled_ionchannels2_CACNA1C genes_all_modeled_ionchannels2_CACNA1D genes_all_modeled_ionchannels2_CACNA1I genes_all_modeled_ionchannels2_KCNB1 genes_all_modeled_ionchannels2_KCNJ6 genes_all_modeled_ionchannels2_KCND3 genes_all_modeled_ionchannels2_HCN1 genes_all_modeled_ionchannels2_ATP2A2 genes_new genes_new_ionchannels genes_new_synaptic genes_new_ionchannels_kcn genes_new_ionchannels_cacn genes_new_ionchannels_scn_and_hcn genes_new_synaptic_others genes_new_synaptic_PKA genes_new_synaptic_PKC genes_new_synaptic_PKAPKC genes_new_synaptic_PKAPKC_noPP)
module load PLINK/1.90b_6.2-x86_64

pwd
for myset in full ${setNames[@]}
do
for THR in 5e-10 0.0005 5e-05 5e-08 5e-09 1e-10 1e-09 1e-08 1e-07 5e-07 1e-06 5e-06 1e-05 0.0001 
do
  if [ -f geneticBD/${myset}_${THR}_PRS.best ]
  then
    echo "${myset}_${THR}_PRS.best already exists"
    continue
  fi
    
  cd i/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out geneticBD/${myset}_${THR}_i --keep subj_IDs_of_interest.txt --recode
  cd ../../ii/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out geneticBD/${myset}_${THR}_ii --keep subj_IDs_of_interest.txt --recode
  cd ../../iii/bed
  plink --bfile all --extract snps_of_interest_${myset}_${THR}.txt --out geneticBD/${myset}_${THR}_iii --keep subj_IDs_of_interest.txt --recode

  cd geneticBD
  for q in ${myset}_${THR}_i ${myset}_${THR}_ii ${myset}_${THR}_iii
  do
      echo "plink --file ${q} --make-bed --out ${q}_new"
      plink --file ${q} --make-bed --out ${q}_new
  done
  
  
  for q in ${myset}_${THR}
  do
    ls ${q}_i*_new*bed |sed 's/.bed//' > mergelist${myset}.txt
    echo "plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged"
    plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged
  done

  PRSice_linux --base daner_pgc4_bd_eur_no23andMe_noTOP_HRCfrq_with_var --chr CHR --bp BP --A1 A1 --A2 A2 --target ${myset}_${THR}_merged --pheno DUMMY_PHENO.txt --stat OR --or --pvalue P --print-snp --no-full --snp SNP --bar-levels "1" --fastscore --out ${myset}_${THR}_PRS --no-clump --maf 0.01 --thread 8 --no-clump --fastscore --no-full


done

done

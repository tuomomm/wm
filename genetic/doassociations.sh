#SCZ PRSs
for q in genes_all_modeled_synaptic genes_full genes_all_modeled_ionchannels2 genes_all_modeled_ionchannels2_CACNA1C genes_all_modeled_ionchannels2_CACNA1D genes_all_modeled_ionchannels2_CACNA1I genes_all_modeled_ionchannels2_KCNB1 genes_all_modeled_ionchannels2_KCNJ6 genes_all_modeled_ionchannels2_KCND3 genes_all_modeled_ionchannels2_HCN1 genes_all_modeled_ionchannels2_ATP2A2 genes_new genes_new_ionchannels genes_new_synaptic genes_new_ionchannels_kcn genes_new_ionchannels_cacn genes_new_ionchannels_scn_and_hcn genes_new_synaptic_others genes_new_synaptic_PKA genes_new_synaptic_PKC genes_new_synaptic_PKAPKC genes_new_synaptic_PKAPKC_noPP
do
    echo python3 extract_othermeasures_PRS_transformed_LNSonly.py $q
    python3 extract_othermeasures_PRS_transformed_LNSonly.py $q
done

for group in both SZ HC
do
    python3 correlations_PRSs_othermeasures_all_pvals_ztransformed_LNSonly.py $group                    #The main results: determine beta coefficients and p-values for models where the LNS measures are fitted against the PRS, age and sex.
    python3 correlations_PRSs_othermeasures_notagainstsex_all_pvals_ztransformed_LNSonly.py $group      #The same but against PRS and age only
    python3 correlations_PRSs_othermeasures_notagainstsexorage_all_pvals_ztransformed_LNSonly.py $group #The same but against PRS (and intercept) only
    python3 correlations_PRSs_othermeasures_againstIQtoo_all_pvals_ztransformed_LNSonly.py $group       #Supplementary results: determine beta coefficients and p-values for models where the LNS measures are fitted against the PRS, age, sex and IQ.
done

#BD PRSs
for q in genes_all_modeled_synaptic genes_full genes_all_modeled_ionchannels2 genes_all_modeled_ionchannels2_CACNA1C genes_all_modeled_ionchannels2_CACNA1D genes_all_modeled_ionchannels2_CACNA1I genes_all_modeled_ionchannels2_KCNB1 genes_all_modeled_ionchannels2_KCNJ6 genes_all_modeled_ionchannels2_KCND3 genes_all_modeled_ionchannels2_HCN1 genes_all_modeled_ionchannels2_ATP2A2 genes_new genes_new_ionchannels genes_new_synaptic genes_new_ionchannels_kcn genes_new_ionchannels_cacn genes_new_ionchannels_scn_and_hcn genes_new_synaptic_others genes_new_synaptic_PKA genes_new_synaptic_PKC genes_new_synaptic_PKAPKC genes_new_synaptic_PKAPKC_noPP
do
    echo python3 extract_othermeasures_BDPRS_transformed_LNSonly.py $q
    python3 extract_othermeasures_BDPRS_transformed_LNSonly.py $q
done

for group in both SZ HC
do
    python3 correlations_BDPRSs_othermeasures_all_pvals_ztransformed_LNSonly.py $group #Determine beta coefficients and p-values for models where the LNS measures are fitted against the PRS, age and sex.
done


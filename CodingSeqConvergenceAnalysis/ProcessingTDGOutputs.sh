#!/bin/bash/

# Script to process the output of convergence analysis via TDG09. 

cat viviparity_convergence_hyp_TDG.sh.o* > ALL_CONVERGENT_AA_RESULTS.txt
sed -n '/LrtResults/,/FullResults/{/FullResults/!p}' ALL_CONVERGENT_AA_RESULTS.txt > SUMMARY_AA_Results.txt
awk '/^Lr/{print "LrtResults_" ++i; next}{print}' < SUMMARY_AA_Results.txt > SUMMARY_AA_Results_headerchange.txt
ALL_CONVERGENT_AA_RESULTS.txt | grep "AlignmentFile" | cut -d'/' -f 11 > AA_names_inorder.txt
less SUMMARY_AA_Results_headerchange.txt | grep "Lrt" > Lrt_names.txt
paste Lrt_names.txt AA_names_inorder.txt > Lrt_names_AA_names.txt
awk 'FNR==NR {dict[$1]=$2; next} {$1=($1 in dict) ? dict[$1] : $1}1' Lrt_names_AA_names.txt SUMMARY_AA_Results_headerchange.txt > Convergence_SUMMARY_fixed_headers.txt
sed -i 's/-//g' Convergence_SUMMARY_fixed_headers.txt
less Convergence_SUMMARY_fixed_headers.txt | awk '{print $6}' > FDR_site_results.txt
sed -i 's/LRT,//g' FDR_site_results.txt
sed -i 's/>//g' FDR_site_results.txt
awk '!NF{$0=">"}1' FDR_site_results.txt > FDR_site_results_1.txt
mv FDR_site_results_1.txt FDR_site_results.txt
sed -i 's/>//g' FDR_site_results.txt | head
sed -i 's/>//g' FDR_site_results.txt
sed -i '/^$/d' FDR_site_results.txt

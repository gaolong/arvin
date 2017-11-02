

arvin_software_dir=./
arvin_annotation_data_dir=$arvin_software_dir/arvin_annotation_data/
arvin_scripts_dir=$arvin_software_dir/process_features_scripts/

input_snps_file=$1
csi_ann_output_file=$2
im_pet_output_file=$3
gwava_features_file=$4
funseq_features_file=$5
output_dir=$6

mkdir -p $output_dir

##################Merge CSI-ANN and IM-PET output##################

echo "Processing ARVIN and IM-PET outputs"
echo "..."
gene_symbol_mapping_file=$arvin_annotation_data_dir/gene_id_mapping.txt
combined_score_file=$output_dir/enhancer_and_ep_scores.txt
perl $arvin_scripts_dir/combine_enhancer_and_ep_scores.pl $csi_ann_output_file $im_pet_output_file $gene_symbol_mapping_file $combined_score_file 
echo "ARVIN and IM-PET outputs processed"


#Overlap input snps with enhancers
echo "Overlapping input SNPs with enhancers."
echo "..."
snp_gene_file=$output_dir/snp_target_gene.txt
perl $arvin_scripts_dir/overlap_snps_with_ep.pl $input_snps_file $combined_score_file $snp_gene_file
echo "SNPs are overlapped with enhancers."

##Process GWAVA and FunSeq features######

echo "Processing GWAVA and FunSeq features."
echo "..."
perl $arvin_scripts_dir/process_funseq1_output.pl $input_snps_file $funseq_features_file $output_dir/features_funseq.csv
perl $arvin_scripts_dir/merge_gwava_and_funseq_features.pl $gwava_features_file $output_dir/features_funseq.csv $output_dir/features_gwava_funseq.csv
echo "GWAVA and FunSeq features processed."


###Compute  disruption##################
echo "Computing TF binding disruption."
echo "..."
echo "Reading sequence..."
genome_dir=$arvin_annotation_data_dir/hg19/
disruption_dir=$output_dir/disruption/
mkdir -p  $disruption_dir
snp_seq_file=$disruption_dir/sequence.txt
perl $arvin_scripts_dir/disruption_snp_get_surrounding_regions.pl $input_snps_file $genome_dir $snp_seq_file
echo "Calculating log odds scores..."
pwm_dir=$arvin_annotation_data_dir/pwm/
TFM_cutoff_file=$arvin_annotation_data_dir/TFM_p_value_4_e_7_cutoff_scores.txt
TFBS_score_dir=$disruption_dir/TFBS_scores/
TFBS_score_changes_file=$disruption_dir/TFBS_score_changes.txt
mkdir -p  $TFBS_score_dir; 
perl $arvin_scripts_dir/disruption_calculate_log_odds_ratio.pl $TFM_cutoff_file $pwm_dir $snp_seq_file $TFBS_score_dir $TFBS_score_changes_file
echo "Calculating TFBS disruption p values..."
Rscript $arvin_scripts_dir/calculate_empirical_p_values.R $TFBS_score_changes_file $arvin_annotation_data_dir/score_changes_sampling/ $disruption_dir/disruption_p.txt
cp $disruption_dir/disruption_p.txt $output_dir/

echo "TF binding disruption computed."

###########################################
echo "Deleting temporary files."
rm -f $output_dir/enhancer_and_ep_scores.txt
rm -f $output_dir/features_funseq.csv
rm -f -r $output_dir/disruption
echo "Temporary files deleted."

echo ""
echo "*******************************************************************"
echo "ARVIN preprocessing step is completed."
echo "Check the following output files:"
echo "TFBS disruption: $output_dir/disruption_p.txt "
echo "GWAVA and Funseq features: $output_dir/features_gwava_funseq.csv"
echo "SNP gene interactions: $output_dir/snp_target_gene.txt"
echo "*******************************************************************"











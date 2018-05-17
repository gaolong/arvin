
#This is the ARVIN package directory
arvin_software_dir=/home/programs/arvin/

#This is the directory where ARVIN perl scripts are stored
scripts_dir=$arvin_software_dir/scripts/

#This is the ARVIN annotation data directory
arvin_annotation_data_dir=$arvin_software_dir/arvin_annotation_data/


input_snps_file=$1 #Input SNPs in bed format
ep_interaction_file=$2 #IM-PET enhancer-promoter interactions 
gwava_features_file=$3 #GWAVA SNP annotations 
funseq_features_file=$4 #FunSeq SNP annotations 
de_file=$5 #Differential expression results for genes 
intermediate_dir=$6 #Output directory
output_dir=$7 #Output directory

#Make the output diectory if it does not exist
mkdir -p $intermediate_dir
mkdir -p $output_dir

#Merge enhancer and enhancer promoter interactions into a matrix#
gene_symbol_mapping_file=$arvin_annotation_data_dir/gene_id_symbol_entrez_mapping.txt
epg_score_file=$intermediate_dir/epg_scores.txt
perl $scripts_dir/process_ep_scores.pl $ep_interaction_file $gene_symbol_mapping_file $epg_score_file 
echo "EP interactions processed"


#Overlap input snps with enhancers
echo "Overlapping input SNPs with enhancers."
echo "..."
snp_gene_file=$intermediate_dir/snp_target_gene.bed.txt
perl $scripts_dir/overlap_snps_with_ep.pl $input_snps_file $epg_score_file $snp_gene_file
echo "SNPs are overlapped with enhancers."
awk 'NR>1 {OFS="\t"; print $6, $8, $9}' $intermediate_dir/snp_target_gene.bed.txt > $intermediate_dir/snp_gene_edge_weights.txt 
cp $intermediate_dir/snp_target_gene.bed.txt $output_dir/


#Process differential expression p-values to use as node weights

#User wants to use precomputed differential expression results from the paper (Gao and Uzun et al, Nature Communications, 2018)
if [ "$de_file" = "CRH" ] || [ "$de_file" = "MS" ] || [ "$de_file" = "PSO" ] || [ "$de_file" = "RA" ] || [ "$de_file" = "SLE" ] || [ "$de_file" = "T1D" ] || [ "$de_file" = "ULC" ] ; then
  precomputed_de=$de_file
  echo "Using precomputed in differential expression data for $precomputed_de"
  cp $arvin_annotation_data_dir/built_in_de_data/gene_node_weights.de_log_p.normalized.${precomputed_de}.txt $intermediate_dir/gene_node_weights.de_log_p.normalized.txt
else # User desires to use his/her own differential expression p-values
  echo "Using file: $de_file for differential expression p-values."
  Rscript $scripts_dir/normalize_dge.R $de_file $arvin_annotation_data_dir/entrez_ids.txt $intermediate_dir/gene_node_weights.de_log_p.normalized.txt
fi





#Process FunSeq features and merge them with GWAVA#
echo "Processing GWAVA and FunSeq features."
echo "..."
perl $scripts_dir/process_funseq_output.pl $input_snps_file $funseq_features_file $intermediate_dir/features_funseq.csv
#perl $scripts_dir/merge_gwava_and_funseq_features.pl $gwava_features_file $intermediate_dir/features_funseq.csv $intermediate_dir/features_gwava_funseq.csv
gwava_features_cmd_file=$intermediate_dir/gwava_features_file.web.csv
perl $scripts_dir/convert_gwava_web_to_gwava_cmd.pl $gwava_features_file $intermediate_dir/gwava_features_file.cmd.csv
perl $scripts_dir/merge_gwava_and_funseq_features.cmd.pl $intermediate_dir/gwava_features_file.cmd.csv $intermediate_dir/features_funseq.csv $intermediate_dir/features_gwava_funseq.csv
sed s/","/"\t"/g  $intermediate_dir/features_gwava_funseq.csv | sed 's/\<NA\>/0.0/g' >  $intermediate_dir/features_gwava_funseq.txt
echo "GWAVA and FunSeq features processed."


#Compute  disruption#
echo "Computing TF binding disruption."
echo "..."
echo "Reading sequence..."
genome_dir=$arvin_annotation_data_dir/hg19/
disruption_dir=$intermediate_dir/disruption/
mkdir -p  $disruption_dir
snp_seq_file=$disruption_dir/sequence.txt
#Get the genomic sequence around the SNPs to test whether there is a TF binding site overlapping#
perl $scripts_dir/disruption_snp_get_surrounding_regions.pl $intermediate_dir/snp_target_gene.bed.txt $genome_dir $snp_seq_file
echo "Calculating log odds scores..."
pwm_dir=$arvin_annotation_data_dir/pwm/
TFM_cutoff_file=$arvin_annotation_data_dir/TFM_p_value_4_e_7_cutoff_scores.txt
TFBS_score_dir=$disruption_dir/TFBS_scores/
TFBS_score_changes_file=$disruption_dir/TFBS_score_changes.txt
mkdir -p  $TFBS_score_dir; 
#Compute the log odds scores for reference and alternate alleles of the SNPs which overlap with a significant TF binding site,
#and measure the difference between the scores of reference and alternate alleles 
perl $scripts_dir/disruption_calculate_log_odds_ratio.pl $TFM_cutoff_file $pwm_dir $snp_seq_file $TFBS_score_dir $TFBS_score_changes_file
echo "Calculating TFBS disruption p values..."
#Calculate empirical p-values for the differences between log odds scores by comparing the calcualted difference with the background
Rscript $scripts_dir/disruption_calculate_empirical_p_values.R $TFBS_score_changes_file $arvin_annotation_data_dir/score_changes_sampling/ $intermediate_dir/disruption_p.txt
awk '{OFS="\t"; print $1, $2, $3, $4, $5, $8, $9}' $intermediate_dir/disruption_p.txt > $output_dir/disruption_p.txt

perl $scripts_dir/disruption_add_missing_snps.pl $intermediate_dir/snp_target_gene.bed.txt $intermediate_dir/disruption_p.txt $intermediate_dir/disruption_p.all.txt 

awk 'NR>1 {OFS="\t"; print $1, $8}'  $intermediate_dir/disruption_p.all.txt > $intermediate_dir/snp_node_weights.disruption_log_p.txt 



echo "TF binding disruption computed."



cd $intermediate_dir/
python $scripts_dir/generate_network_input.py \
   $intermediate_dir/snp_node_weights.disruption_log_p.txt \
   $intermediate_dir/snp_gene_edge_weights.txt \
   $intermediate_dir/gene_node_weights.de_log_p.normalized.txt \
   $arvin_annotation_data_dir/humannet.txt \
   $intermediate_dir/Edge.txt \
   $intermediate_dir/Node.txt 



Rscript $scripts_dir/run_scoring.R \
        $arvin_annotation_data_dir/Aug22_OPT_classification.RData \
        $arvin_annotation_data_dir/Aug22_CV_classification.RData \
        $intermediate_dir/snp_node_weights.disruption_log_p.txt \
        $intermediate_dir/features_gwava_funseq.txt \
        $intermediate_dir/Edge.txt \
        $intermediate_dir/Node.txt \
        $output_dir/arvin_scores.txt 



#Delete the temporary files
#echo "Deleting temporary files."
rm -f -r $intermediate_dir/disruption
#echo "Temporary files deleted."

#Report for the output files
echo ""
echo "*******************************************************************"
echo "ARVIN completed running."
echo "Check the following output files."
echo "ARVIN scores: $output_dir/arvin_scores.txt "
echo "SNP gene interactions: $output_dir/snp_target_gene.txt"
echo "TFBS disruption: $output_dir/disruption_p.txt "
echo "*******************************************************************"








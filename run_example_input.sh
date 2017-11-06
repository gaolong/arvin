mkdir -p arvin_output
output_dir=$PWD/arvin_output/

cd /home/uzuny/arvin/revision/software/v0.5/

./process_features.sh \
   /mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_2/snps.bed \
   /mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_2/Enhancers.Treg.txt \
   /mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_2/EP_Interactions.Treg.bed \
   /mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_2/gwava_output.csv \
   /mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_2/funseq_output.bed \
   $output_dir

#input_snps_file=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_data/gold_esnps.bed
#csi_ann_output_file=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_data/CSI-ANN.txt
#im_pet_output_file=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_data/IM-PET.txt
#gwava_features_file=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_data/features_gwava.csv
#funseq_features_file=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/example_input_data/features_funseq.bed
#output_dir=/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.5/output_data/gold_esnps/




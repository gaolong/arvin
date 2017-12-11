output_dir=$PWD/arvin_output/
rm -rf $output_dir
mkdir -p $output_dir


./process_features.sh \
   example_input/snps.bed \
   example_input/Enhancers.Treg.txt \
   example_input/EP_Interactions.Treg.bed \
   example_input/gwava_output.csv \
   example_input/funseq_output.bed \
   $output_dir





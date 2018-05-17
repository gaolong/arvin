cd /home/programs/arvin/example_2/

input_dir=/home/programs/arvin/example_2/input/
intermediate_dir=/home/programs/arvin/example_2/intermediate/
output_dir=/home/programs/arvin/example_2/output/

/home/programs/arvin/scripts/arvin_wrapper.sh \
    $input_dir/input_snps_file.bed \
    /home/programs/arvin/arvin_annotation_data/built_in_ep_data/EP_AUTOIMMUNE.bed \
    $input_dir/gwava_features_file.tsv \
    $input_dir/funseq_features_file.bed \
    T1D \
    $intermediate_dir/ \
    $output_dir/



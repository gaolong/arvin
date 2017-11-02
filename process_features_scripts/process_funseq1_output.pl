use strict;

#my $funseq_output_file="/home/uzuny/arvin/revision/funseq/output_psnps_hgmd_0.5/Output.bed";
#my $corrected_output_file="/home/uzuny/arvin/revision/funseq/output_psnps_hgmd_0.5/Output.corrected.bed";

my $input_snps_file=$ARGV[0];
my $funseq_output_file=$ARGV[1];
my $corrected_output_file=$ARGV[2];

open(SNPS, $input_snps_file) or die("Cannot read $input_snps_file \n ");
open(funseq, $funseq_output_file) or die("Cannot read $funseq_output_file \n ");
open(OUT, ">$corrected_output_file") or die("Cannot WRITE $corrected_output_file \n ");

my %funseq_results = ();
my @feat_indeces = (9,11,12,10,13,7);

while(<funseq>)
{
   chomp;
   my @array = split /[;\t]/;

   my $key = $array[0]."\t".$array[2];
  
   my @feature_array = ();

   for(my $i=0; $i<scalar(@feat_indeces); $i++)
   {
      my $feature = $array[$feat_indeces[$i]];
      if($feature eq "." || $feature eq "0")
      {
			$feature_array[$i] = "0";
      } 
      else
      {
			$feature_array[$i] = "1";
      }
   }
   $funseq_results{$key} = join(",", @feature_array);
#   print OUT join("\t", @array);

}

#print OUT "SNP_id,is_annotated_in_encode,is_sensitive,is_ultrasensitive,is_motif_breaking,target_gene_is_hub,gerp,is_ultra_conserved,is_hot_region,target_gene_known,is_recurrent,is_recurrent_in_cancer_dataset,noncoding_score\n";
print OUT "SNP_id,is_annotated_in_encode,is_sensitive,is_ultrasensitive,is_motif_breaking,target_gene_known,target_gene_is_hub\n";

while(<SNPS>)
{
   chomp;
   my @array = split /\t/;
   my $snp_id = $array[5];
   my $key = $array[0]."\t".$array[2];
   print OUT $snp_id;
   
   if(exists($funseq_results{$key}) )
   {
      print OUT ",".$funseq_results{$key}."\n";
   } 
   else
   {
      for(my $i=0; $i<scalar(@feat_indeces); $i++)
      {
         print OUT ",NA";
      }#for
      print OUT "\n";
   }#else

}#while


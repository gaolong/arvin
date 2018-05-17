use strict;

#This script reads the FunSeq output, which is in bed format and the features are separated with semicolons
#and writes the features into output file in csv format

#1st argument is the input snps file in bed format
my $input_snps_file=$ARGV[0];

#2nd argument is the FunSeq output format in bed format
my $funseq_output_file=$ARGV[1];

#3rd argument is the output file containing FunSeq features in csv format
my $corrected_output_file=$ARGV[2];

#Open the input files for reading
open(SNPS, $input_snps_file) or die("Cannot read $input_snps_file \n ");
open(funseq, $funseq_output_file) or die("Cannot read $funseq_output_file \n ");

#Open the output file for writing
open(OUT, ">$corrected_output_file") or die("Cannot WRITE $corrected_output_file \n ");

#Hash array storing the FunSeq results
my %funseq_results = ();

#Feature indeces 
my @feat_indeces = (9,11,12,10,13,7);

while(<funseq>)
{
   chomp;
	#Split the fields separated by tabs and semicolons into array
   my @array = split /[;\t]/;

   my $key = $array[0]."\t".$array[2];
  
   my @feature_array = ();

   for(my $i=0; $i<scalar(@feat_indeces); $i++)
   {
      my $feature = $array[$feat_indeces[$i]];
		#If the field is "." turn into 0
      if($feature eq "." || $feature eq "0")
      {
			$feature_array[$i] = "0";
      }#Otherwise, it is 1
      else
      {
			$feature_array[$i] = "1";
      }
   }
   $funseq_results{$key} = join(",", @feature_array);

}

#Print the header to output file
print OUT "SNP_id,is_annotated_in_encode,is_sensitive,is_ultrasensitive,is_motif_breaking,target_gene_known,taget_gene_is_hub\n";

#Read the input SNPs in bed format
#If there is a matching FunSeq SNP, print the features in csv format
#If not, print NAs as features
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




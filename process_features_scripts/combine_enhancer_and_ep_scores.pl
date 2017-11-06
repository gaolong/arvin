# November 2017 by YU
# This script processes enhancer prediction results and enhancer-promoter interactions 
# and combines enhancer and ep-prediction scores into one file

use strict;

#1st argument: Enhancer prediction results file in tsv format: 
#chrom enhancer_center enhancer_score
my $csi_ann_output_file = $ARGV[0]; 

#2nd argument: Enhancer-promoter interaction results file in tsv format: 
#chrom enhancer_start enhancer_end transcript_id ep_score
my $im_pet_output_file = $ARGV[1];

#3rd argument: File for mapping transcript ids to gene symbols and gene ids in tsv format:
# transcript_id gene_symbol gene_id
my $gene_symbol_mapping_file = $ARGV[2];

#4th arguemnt is the output file containing enhancer and ep prediction scores into one file in tsv format:
#chrom enhancer_start enhancer_end enhancer_score transcript_id gene_symbol gene_id ep_score enhancer_ep_score
my $combined_score_file = $ARGV[3];

#Open the input files for reading
open(CSI, $csi_ann_output_file) or die("Error: Cannot read CSI-ANN output file: $csi_ann_output_file \n");
open(IM, $im_pet_output_file) or die("Error: Cannot read IM-PET output file: $im_pet_output_file \n");
open(GENE, $gene_symbol_mapping_file) or die("Error: Cannot read gene symbol mapping file: $gene_symbol_mapping_file \n");

#Open the output file for writing
open(COMBINED, ">$combined_score_file") or die("Error: Cannot write output file: $combined_score_file \n");

#Hash array to store the enhancer scores
my %enh_scores = ();

#Read the enhancer records and store the scores into hash array
while(<CSI>)
{
   chomp;   
   my @array = split /\t/;
   my $chrom = $array[0];
   my $center = $array[1];
   my $score = $array[2];
   	
   $enh_scores{$chrom."\t".$center} = $score;	
} 

#Hash array to store the enhancer scores
my %genes = ();

#Read the gene mapping and store them into hash array
while(<GENE>)
{
   chomp;   
   my @array = split /\t/;
   my $transcript_id = $array[0];
   my $gene_symbol = $array[1];
   my $gene_id = $array[2];
   	
   $genes{$transcript_id} = $gene_symbol."\t".$gene_id;	
} 

#Print the output file header
print COMBINED "chrom"."\t"."start"."\t"."end"."\t"."enh_score"."\t"."transcript_id"."\t"."gene_symbol"."\t"."gene_id"."\t"."ep_score"."\t"."enh_ep_score"."\n";

#Read the EP-interactions, map the transcripts to genes and write the records to output file
<IM>;
while(<IM>)
{
   chomp;   
   my @array = split /\t/;
   my $chrom = $array[0];
   my $start = $array[1];
   my $end = $array[2];
   my $transcript_id = $array[3];
   my $ep_score = $array[4];	
   #Find the enhancer center	   	
	my $center = ($start + $end) / 2;	
		
	my $enh_score;
	my $gene;

	if(exists($enh_scores{$chrom."\t".$center}) )
   {
   	$enh_score = $enh_scores{$chrom."\t".$center};	
	}else #Enhancer score is not found for an EP interaction prediction
   {
		die("Error: Enhancer prediction score not found for the enhancer with coordinates: $chrom:$start-$end \n");
   }

	if(exists($genes{$transcript_id}) )
   {
   	$gene = $genes{$transcript_id};	
	}else
	{
		next;
	}	

	my $enh_ep_score = $enh_score * $ep_score;

	
	print COMBINED $chrom."\t".$start."\t".$end."\t".$enh_score."\t".$transcript_id."\t".$gene."\t".$ep_score."\t".$enh_ep_score."\n"; 

} 





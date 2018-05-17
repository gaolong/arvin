# November 2017 by YU
# This script processes enhancer prediction results and enhancer-promoter interactions 
# and combines enhancer and ep-prediction scores into one file

use strict;



#1st argument: Enhancer-promoter interaction results file in tsv format: 
#chrom enhancer_start enhancer_end transcript_id ep_score
my $im_pet_combined_file = $ARGV[0];

#2nd argument: File for mapping transcript ids to gene symbols and gene ids in tsv format:
# transcript_id gene_symbol gene_id
my $gene_symbol_mapping_file = $ARGV[1];

#3rd arguemnt is the output file containing enhancer and ep prediction scores into one file in tsv format:
#chrom enhancer_start enhancer_end enhancer_score transcript_id gene_symbol gene_id ep_score enhancer_ep_score
my $combined_file = $ARGV[2];

#Open the input files for reading

open(IM, $im_pet_combined_file) or die("Error: Cannot read IM-PET output file: $im_pet_combined_file \n");
open(GENE, $gene_symbol_mapping_file) or die("Error: Cannot read gene symbol mapping file: $gene_symbol_mapping_file \n");

#Open the output file for writing
open(COMBINED, ">$combined_file") or die("Error: Cannot write output file: $combined_file \n");
 


#Hash array to store the gene_symbols
my %gene_symbols = ();
my %entrez_ids = ();

#Read the gene mapping and store them into hash array
print("Reading gene_symbols \n");
while(<GENE>)
{
   chomp;   
   my @array = split /\t/;
   my $transcript_id = $array[0];
   my $gene_symbol = $array[1];
   my $entrez_id = $array[2];
   	
   $gene_symbols{$transcript_id} = $gene_symbol;	
   $entrez_ids{$transcript_id} = $entrez_id;	
} 

#Print the output file header
print COMBINED "chrom"."\t"."start"."\t"."end"."\t"."transcript_id"."\t"."gene_symbol"."\t"."entrez_id"."\t"."ep_score"."\n";

print("Reading EP interactions \n");
#Read the EP-interactions, map the transcripts to gene_symbols and write the records to output file
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
	my $gene_symbol = "XXX";
	my $entrez_id = "YYY";

	if(exists($gene_symbols{$transcript_id}) )
	{
		$gene_symbol = $gene_symbols{$transcript_id};	
      $entrez_id = $entrez_ids{$transcript_id};	
	}else
	{
		next;
	}	

	print COMBINED $chrom."\t".$start."\t".$end."\t".$transcript_id."\t".$gene_symbol."\t".$entrez_id."\t".$ep_score."\n"; 

} 





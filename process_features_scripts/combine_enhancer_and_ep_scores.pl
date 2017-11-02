use strict;

my $csi_ann_output_file = $ARGV[0];
my $im_pet_output_file = $ARGV[1];
my $gene_symbol_mapping_file = $ARGV[2];
my $combined_score_file = $ARGV[3];

open(CSI, $csi_ann_output_file) or die("Error: Cannot read CSI-ANN output file: $csi_ann_output_file \n");
open(IM, $im_pet_output_file) or die("Error: Cannot read IM-PET output file: $im_pet_output_file \n");
open(GENE, $gene_symbol_mapping_file) or die("Error: Cannot read gene symbol mapping file: $gene_symbol_mapping_file \n");
open(COMBINED, ">$combined_score_file") or die("Error: Cannot write output file: $combined_score_file \n");

my %enh_scores = ();

while(<CSI>)
{
   chomp;   
   my @array = split /\t/;
   my $chrom = $array[0];
   my $center = $array[1];
   my $score = $array[2];
   	
   $enh_scores{$chrom."\t".$center} = $score;	
} 

my %genes = ();

while(<GENE>)
{
   chomp;   
   my @array = split /\t/;
   my $transcript_id = $array[0];
   my $gene_symbol = $array[1];
   my $gene_id = $array[2];
   	
   $genes{$transcript_id} = $gene_symbol."\t".$gene_id;	
} 

print COMBINED "chrom"."\t"."start"."\t"."end"."\t"."enh_score"."\t"."transcript_id"."\t"."gene_symbol"."\t"."gene_id"."\t"."ep_score"."\t"."enh_ep_score"."\n";

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
	my $center = ($start + $end) / 2;	
		
	my $enh_score;
	my $gene;

	if(exists($enh_scores{$chrom."\t".$center}) )
   {
   	$enh_score = $enh_scores{$chrom."\t".$center};	
	}else
   {
		die("Error: Prediction score not found for the enhancer with coordinates: $chrom:$start-$end \n");
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










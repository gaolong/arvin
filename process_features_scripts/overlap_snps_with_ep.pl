use strict;

my $snp_file = $ARGV[0];
my $combined_score_file = $ARGV[1];
my $snp_gene_score_file = $ARGV[2];


open(SNP_POSITIONS, $snp_file) or die("Error: Cannot read CSI-ANN output file: $snp_file \n");
open(COMBINED_SCORE, $combined_score_file) or die("Error: Cannot open enhancer-ep score file: $combined_score_file \n");
open(SNP_SCORE, ">$snp_gene_score_file") or die("Error: Cannot write snp score file: $snp_gene_score_file \n");

my %combined_scores = ();

while(<COMBINED_SCORE>)
{
   chomp;   
   my @array = split /\t/;
   my $chrom = $array[0];
   my $start = $array[1];
   my $end = $array[2];
   my $gene_symbol = $array[5];
   my $gene_id = $array[6];

   my $combined_score = $array[scalar(@array)-1];
	my $key = $chrom."\t".$start."\t".$end."\t".$gene_symbol."\t".$gene_id;	   	

	if(exists($combined_scores{$key}))
   {
      if($combined_score > $combined_scores{$key})
      {
         $combined_scores{$key} = $combined_score;
      }
   }
	else
   {
		 $combined_scores{$key} = $combined_score;	
	} 
	
  
} 

my %snp_gene_scores = ();

while(<SNP_POSITIONS>)
{
   chomp;   
   my @array = split /\t/;
   my $snp_chrom = $array[0];
   my $snp_start = $array[1];
   my $snp_end = $array[2];
   my $snp_ref = $array[3];
   my $snp_alt = $array[4];
   my $snp_id = $array[5];

	foreach my $key (keys %combined_scores)
	{
		my @ep = split /\t/, $key;
		my $enh_chrom = $ep[0];
		my $enh_start = $ep[1];
		my $enh_end = $ep[2];
		my $gene_symbol = $ep[3];
		my $gene_id = $ep[4];
		my $combined_score = $combined_scores{$key};
		
		if($snp_chrom eq $enh_chrom && $snp_start > $enh_start && $snp_end < $enh_end)
		{
			#print "SNP overlapping with enhancer: ".$_."\t".$key."\n";	
			
			if( exists($snp_gene_scores{$snp_id."\t".$gene_symbol} ) )
			{
				if($combined_score > $snp_gene_scores{$snp_id."\t".$gene_symbol} )
				{
					$snp_gene_scores{$snp_id."\t".$gene_symbol."\t".$gene_id} = $combined_score;
				}

			}#if( exists
			else
			{
				$snp_gene_scores{$snp_id."\t".$gene_symbol."\t".$gene_id} = $combined_score;
			}

		}#if($snp_chrom
		
	}#foreach my $key (keys %combined_scores)	
   
}#while(<SNP_POSITIONS>)

print SNP_SCORE "snp_id"."\t"."gene_symbol"."\t"."entrez_gene_id"."\t"."interaction_score"."\n";
foreach my $key (keys %snp_gene_scores)
{
	print SNP_SCORE $key."\t".$snp_gene_scores{$key}."\n";
}





use strict;

my $snp_file = $ARGV[0];
my $disruption_file = $ARGV[1];
my $output_file = $ARGV[2];

open(SNP, $snp_file) or die("Cannot read $snp_file \n");
open(DISRUPTION, $disruption_file) or die("Cannot read $disruption_file \n");
open(OUT, ">$output_file") or die("Cannot write $output_file \n");

my %snps = ();
while(<DISRUPTION>)
{
   print OUT;
   my @array = split /\t/;
   my $snp_id = $array[0];
   $snps{$snp_id} = 1;
}


while(<SNP>)
{
   chomp;
   my @array = split /\t/;
   my $snp_id = $array[5];
   if(exists($snps{$snp_id}) )
   {
      next
   }
   
   my $ref = $array[3];
   my $alt = $array[4];
   
   print OUT $snp_id."\t".$ref."\t".$alt."\t"."NA"."\t"."NA"."\t"."0.9999"."\t"."0.9999"."\t"."0.00001"."\t"."0"."\t"."0"."\n";  

}#while




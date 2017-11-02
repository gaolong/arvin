use strict;

use Env;
use List::MoreUtils qw(first_index);

#$arvin_dir = $ENV{'ARVIN_DIR'};
#print $arvin_dir."\n";

my $gwava_file = $ARGV[0];
my $funseq_file = $ARGV[1];
my $outfile = $ARGV[2];

#$gwava_file="$arvin_dir/features/features_gwava/gwava_features_psnps_hgmd.bed";
#$funseq_file="$arvin_dir/features/features_funseq/funseq_features_psnps_hgmd.bed";
#$outfile = "$arvin_dir/features/features_gwava_and_funseq_psnps_hgmd.bed";


unless(open (GWAVA_FILE, $gwava_file) ) { die ("cannot open $gwava_file \n"); }
unless(open (FUNSEQ_FILE, $funseq_file) ) { die ("cannot open $funseq_file \n"); }
unless(open (OUTFILE, ">$outfile") ) { die ("cannot open $outfile \n"); }

my $gwava_header = <GWAVA_FILE>;

my $gwava_header="snp_id".$gwava_header;
#print $gwava_header;
chomp $gwava_header;
my @array = split(/,/, $gwava_header);
splice @array, scalar(@array) - 1, 1;
$gwava_header = join(",",@array);

my %features_gwava = ();

while(my $line = <GWAVA_FILE>)
{
  chomp $line;
  my @array = split(/,/, $line);
  my $snp_id = $array[0];
  splice @array, scalar(@array) -1, 1;
  $features_gwava{$snp_id} = join(",",@array);
  #print $key."\n";
}

my $funseq_header = <FUNSEQ_FILE>;

#print $funseq_header;
my @array = split(/,/, $funseq_header);

splice @array, 0, 1;
my $funseq_header = join(",",@array);
#print "$funseq_header \n";

my $header = $gwava_header.",".$funseq_header ;
#print $header."\n";
print OUTFILE $header;

while(my $line = <FUNSEQ_FILE>)
{
	my @array = split(/,/, $line);
   
   my $snp_id = $array[0];

   if(exists $features_gwava{$snp_id})
   {
		 splice @array, 0, 1;
		 my $features = join(",",@array);
		 print OUTFILE $features_gwava{$snp_id}.",".$features;
   }
   
}#while



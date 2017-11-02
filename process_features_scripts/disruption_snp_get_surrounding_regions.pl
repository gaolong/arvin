#this code is used to prepate +- 25 bp region boundaries 

use strict;

my $snp_file = $ARGV[0];
my $genome_dir = $ARGV[1];
my $out_file = $ARGV[2];


unless (open (SNP_FILE, $snp_file)){ die ("Cannot read input SNP file $snp_file\n"); }  
unless (open (OUT_FILE, ">$out_file")){ die ("Cannot write sequence file $out_file\n"); }  

#chr16	rs104895488	50744886	50744887	G	A	IBD
#chr1	rs4651138	183001311	183001312	C	A	+
#chr1	rs7418179	858801	A	G
while(my $line=<SNP_FILE>)
{
   chomp($line);
   if(!($line =~ /chr23/))
   {
		my @array = split(/\t/,$line);
		my $chrom = lc $array[0];
		my $start = $array[2]-25;
		my $end = $array[2]+25;
		my $ref = $array[3];
		my $alt = $array[4];
		
		my $snp_id = $array[5];

		if($chrom eq "chrx"){
		    $chrom = "chrX";
		}
		if($chrom eq "chry"){
		    $chrom = "chrY";
		}

      my $chrom_file = $genome_dir."/".$chrom;
     
      open(CHROM_FILE, $chrom_file) or die("Cannot read $chrom_file \n");
      my $seq = "";
      read CHROM_FILE, $seq, $start;
      read CHROM_FILE, $seq, 50;
      close(CHROM_FILE);
      #print ">".$chrom."\t".$start."\t".$end."\t".$snp_id."\t".$ref."\t".$alt."\n".uc($seq)."\n";
      print OUT_FILE $chrom.":".$start."-".$end."\t".$snp_id."\t".$chrom."\t".($start+24)."\t".$ref."\t".$alt."\t".uc($seq)."\n";
   }

}

close(OUT_FILE);





#This script is used to read +- 25 bp surroundng region of SNPs
#The obtained sequence is used to scan TF binding specificity in subsequent steps 

use strict;

#1st argument is the list of SNPs in bed format as follows:
#chr9	22124476	22124477	A	G	rs10757278
my $snp_file = $ARGV[0];

#2nd argument is the directory containing full genome sequence 
#Separated into files, each file containing one chromosome sequence
#It is in the annotation data directory
my $genome_dir = $ARGV[1];

#3rd argument is the output file that containes the +- 25 bp surrounding region
#genomic sequence together with the reference and alternate alleles
my $out_file = $ARGV[2];


unless (open (SNP_FILE, $snp_file)){ die ("Cannot read input SNP file $snp_file\n"); }  
unless (open (OUT_FILE, ">$out_file")){ die ("Cannot write sequence file $out_file\n"); }  

#Read the SNP file line by line
while(my $line=<SNP_FILE>)
{
   chomp($line);
   if(!($line =~ /chr23/ or $line =~ /chrom/ ) )
   {
		#Solit the line into array
		my @array = split(/\t/,$line);
		#Obtain the +- 25 bp surrounding coordinates
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

		#Open the chromosome file
      my $chrom_file = $genome_dir."/".$chrom;
      open(CHROM_FILE, $chrom_file) or die("Cannot read $chrom_file \n");
      my $seq = "";
		#Read until the start of +- 25 bp region 
      read CHROM_FILE, $seq, $start;
      read CHROM_FILE, $seq, 50;
      close(CHROM_FILE);

		#Print the sequence into the output file
      print OUT_FILE $chrom.":".$start."-".$end."\t".$snp_id."\t".$chrom."\t".($start+24)."\t".$ref."\t".$alt."\t".uc($seq)."\n";
   }

}

#Close the output file
close(OUT_FILE);





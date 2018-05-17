use strict;

#my $gwava_web_file = '/home/uzuny/arvin/oklahoma/data_01/data_from_julio/gwava_outputSLE.csv';
#my $gwava_cmd_file = '/home/uzuny/arvin/oklahoma/data_01/intermediate/gwava_annot.cmd.csv';

#my $gwava_web_file = '/home/uzuny/arvin/oklahoma/data_03/intermediate/temp.csv';
#my $gwava_cmd_file = '/home/uzuny/arvin/oklahoma/data_01/intermediate/temp.cmd.csv';

my $gwava_web_file = $ARGV[0];
my $gwava_cmd_file = $ARGV[1];

open(WEB, $gwava_web_file) or die("Cannot read $gwava_web_file \n");
open(CMD, ">$gwava_cmd_file") or die("Cannot write $gwava_cmd_file \n");


my $header = <WEB>;
print $header;
print "*************\n";
print "*************\n";
chomp $header;


$header =~ s/Chromosome/chr/; 
$header =~ s/Position/end/;
$header =~ s/TSS distance/tss_dist/;
$header =~ s/SS distance/ss_dist/;
$header =~ s/\%GC/GC/;
#$header =~ s/in_cpg/cpg_island/;
$header =~ s/Average GERP/avg_gerp/;
$header =~ s/GERP/gerp/;
$header =~ s/Average het/avg_het/;
$header =~ s/Average DAF/avg_daf/;
#$header =~ s///;
#$header =~ s///;
#$header =~ s///;
#$header =~ s///;


print $header ."\n";

$header =~ s/\"//g;

print "-----\n";
print $header."\n";
print "-----\n";

my @array = split /\t/, $header;
my $id = $array[1];
my $chrom = $array[2];
my $end = $array[3];
my $start = "start";

print "*************\n";
print "*************\n";

print CMD ",$chrom,$end,$start,cls";

for(my $i=7; $i<scalar(@array); $i++)
{
   print CMD ",".$array[$i];  
}
print CMD "\n";

while(<WEB>)
{
   chomp;
      
   #print ;
   #print "\n";
	my @array = split /\t/;
	my $id = $array[1];
	my $chrom = $array[2];
	my $end = $array[3];
   $end =~ s/,//g;
	my $start = $end-1;
   $chrom = "chr".$chrom;
   
	print CMD "$id,$chrom,$end,$start,1";

	for(my $i=7; $i<scalar(@array); $i++)
	{
         my $field = $array[$i];
			if($field eq "-")
			{
				$field = "0.0";
			}
			if(index($field, "%") != -1)
			{
				$field =~ s/\%//;
				$field = $field/100;
			}
		print CMD ",".$field; 
       
	}
	print CMD "\n";  

}#while









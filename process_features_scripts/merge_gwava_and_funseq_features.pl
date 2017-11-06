use strict;

#This script reads the outputs of GWAVA and FunSeq and combines them into a single file in matrix format, as output

#1st argument is the GWAVA annotation output, containing the sequence features in csv format
my $gwava_file = $ARGV[0];

#2nd argument is the post-processed FunSeq output, in csv format
my $funseq_file = $ARGV[1];

#3rd argument is the output file in csv format, which contains the GWAVA output with FunSeq
my $outfile = $ARGV[2];

#Open the input files for reading
unless(open (GWAVA_FILE, $gwava_file) ) { die ("cannot open $gwava_file \n"); }
unless(open (FUNSEQ_FILE, $funseq_file) ) { die ("cannot open $funseq_file \n"); }

#Open the output file for writing
unless(open (OUTFILE, ">$outfile") ) { die ("cannot open $outfile \n"); }

#Read GWAVA header line
my $gwava_header = <GWAVA_FILE>;
chomp $gwava_header;

my $gwava_header="snp_id".$gwava_header;

#Omit the class field in the end
my @array = split(/,/, $gwava_header);
splice @array, scalar(@array) - 1, 1;
$gwava_header = join(",",@array);

#Hash array storing the GWAVA features
my %features_gwava = ();

#Read the GWAVA file line by line and store it into array, omiting the last field, which is class
while(my $line = <GWAVA_FILE>)
{
  chomp $line;
  my @array = split(/,/, $line);
  my $snp_id = $array[0];
  splice @array, scalar(@array) -1, 1;
  $features_gwava{$snp_id} = join(",",@array);
}

#Hash array storing the FunSeq features
my $funseq_header = <FUNSEQ_FILE>;
my @array = split(/,/, $funseq_header);

#Omit the snp id field from the header
splice @array, 0, 1;
my $funseq_header = join(",",@array);

#Merge the headers and prin it into output file
my $header = $gwava_header.",".$funseq_header ;
print OUTFILE $header;

#Read the FunSeq output, merge with GWAVA and print it into output file
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



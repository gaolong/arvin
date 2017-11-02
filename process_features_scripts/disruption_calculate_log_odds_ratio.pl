use List::MoreUtils qw(uniq);
use Math::Round;


#this code is used to obtain the snp log ratio score

#using the window slid to obtain the snp string
#for each snp we obtain two scores 
my $tf_file = $ARGV[0];

#step 1: first obtain the strings slided windows 
$pwm_dir = $ARGV[1];
opendir(DIR,  $pwm_dir ) || die("Cannot open input dir: $pwm_dir\n");  

#my $seq_file = '/Shared/Tan/yuzun/diabetes/disruption/data/TF_1773/background/snp_1000Genome_1M_sequence';
my $seq_file = $ARGV[2];
unless (open (SEQ_FILE,$seq_file)){ die ("Cannot open input sequence file $seq_file\n"); }  

#chr1:150994910-150994959	rs187922202	chr1	150994934	T(4)	G(5)
#TAAGAGACGTAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCACTGGCGCGA(6) [point 25]
#generating the 


#my $matrix_folder = '/Shared/Tan/yuzun/diabetes/disruption/data/TF_1773/pwm/';

my $output_folder = $ARGV[3];

my $score_changes_file = $ARGV[4];

my %cutoffs = ();
my %tf_len = ();
my %tf_pwm = ();

unless (open (TF_FILE,$tf_file)){ die ("Cannot open input file $tf_file\n"); }  

unless (open (SCORE_CHANGES_FILE, ">$score_changes_file")){ die ("Cannot write output file $score_changes_file\n"); }  

#while (my $file_name = readdir(DIR)) 
#{ next if ($file_name =~ m/^\./);

while(my $line=<TF_FILE>)
{
	chomp($line);
	my @array_TFM_cutoff = split(/\s+/,$line);
	my $file_name = $array_TFM_cutoff[0];
	$file_name =~ s/meme/pwm/;
	$file_name = $file_name.'.pwm';
	#print $file_name."\n";


	chomp($file_name);

	my $file_path = $pwm_dir.$file_name;

	my $motif_name = $file_name;
	$motif_name =~ s/.pwm//;

	#my $file_name = $matrix_folder."$line";

	unless (open (MOTIF_FILE,$file_path)){ die ("Cannot open input pwm file: $file_path\n"); }  
	my $pwm_matrix;

	while(my $line1=<MOTIF_FILE>)
	{
		chomp($line1);
		my @array = split(/\s/,$line1);
		#print $line1."\n";
		if($array[0] =~ /A/)
		{ 
			$pwm_matrix = $array[1];
			for(my $i=2;$i<scalar(@array);$i++)
			{
				$pwm_matrix[0][$i-1] = $array[$i];
				#print $pwm_matrix[0][$i-1]."\n";
				$pwm_matrix = $pwm_matrix."\t".$array[$i];
			}
				                        
		}elsif($array[0] =~ /C/)
		{ 
			for(my $i=1;$i<scalar(@array);$i++)
			{
				$pwm_matrix = $pwm_matrix."\t".$array[$i];
			}
				                        
		}elsif($array[0] =~ /G/)
		{ 
			for(my $i=1;$i<scalar(@array);$i++)
			{
				$pwm_matrix = $pwm_matrix."\t".$array[$i];
			}
				                        
			}else
			{ 
				for(my $i=1;$i<scalar(@array);$i++)
				{
					$pwm_matrix = $pwm_matrix."\t".$array[$i];
				}
				                        
			} 
		$tf_len{$motif_name} = scalar(@array) -1;

	}

	$tf_pwm{$motif_name} = $pwm_matrix;
	$cutoffs{$motif_name} = $array_TFM_cutoff[2];

	close(MOTIF_FILE);

}

$motif_counter = 0;
my $start = time;

foreach my $key (keys %tf_pwm)
{
		$motif_counter++;
		#print $motif_counter." : ".$key."\n";

		my $output_file = $output_folder.$key;

		open (OUTFILE, ">$output_file") || die("Cannot open output file: $output_file\n"); 

		my @pwm_matrix = $tf_pwm{$key};
		my $motif_len = $tf_len{$key};
		my $seq_file = $seq_file;

		
		#print $motif_len."\n";
		#print $tf_pwm{$key}."\n";

		my $value = $tf_pwm{$key};
		my @array = split(/\t/,$value); 
		#print"********************************\n";
		#print $value."\n";
		#print scalar(@array)."\n";
		my @pwm_matrix;
		my $index_1 = 0;
		my $index_2 = 0;
		for(my $i=0;$i<scalar(@array);$i++)
		{
			$pwm_matrix[$index_1][$index_2] = $array[$i];
			$index_2 = $index_2 +1;

			if($index_2%$motif_len eq 0)
			{
			$index_1 = $index_1+1;
			$index_2 = $index_2%$motif_len;
			}  

		}


	## print scalar(@{$pwm_matrix[0]})."\n";
	#   for(my $i=0;$i<4;$i++)
	#    {
	#      print $i."\n";
	#      for(my $j=0;$j<$motif_len;$j++)
	#      {
	#        print $pwm_matrix[$i][$j]." ";
	#       }
	#       print "\n";
	#       }
	# print"********************************\n";

	#*************************** here we do the window slid ****************************
	#chr1:150994910-150994959	rs187922202	chr1	150994934	T	G	
	#TAAGAGACGTAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCACTGGCGCGA

	#chr1:94240734-94240783	rs12062272	chr1	94240758	G	C	CTAGTGCCTGTGAACTTACCTAACGACCAGAGGCAGGTGAATCTCATGGT

	unless (open (SEQ_FILE,$seq_file)){ die ("cannot open input file $seq_file\n"); } 
	  
	while(my $line1=<SEQ_FILE>)
	{
		chomp($line1);
		my @array = split(/\t+/,$line1);
		my @alts = split(/,/,$array[5]);
      my $snp_id = $array[1];		
      my $ref = $array[4];
      my $alt = $array[5];

		for(my $k=0;$k<scalar(@alts);$k++)
		{
		   
			my $score = $array[0]."\t".$array[1]."\t".$array[2]."\t".$array[3]."\t".$array[4]."\t".$alts[$k]."\t".$motif_len;

			for(my $i=0;$i<$motif_len;$i++)
			{

            
				my $start = 25-$motif_len+$i;
				#my $end = $start +$motif_len-1;
				#print $motif_len."\t".$i."\n";
				#print $start."\t".$end."\n";
				my $ref_seq = '';
				my $alt_seq = '';
            
				my $temp_seq = substr $array[6], $start, $motif_len;
            #print "tempseq: $temp_seq\n";
				my $temp_seq_start = substr $temp_seq, 0, $motif_len -$i-1;
				my $temp_seq_end = substr $temp_seq, $motif_len-$i, $i;
            #print "\n";
				#print "1:".$temp_seq_start."\n";
				#print "2:".$temp_seq_end."\n";
				$ref_seq = $temp_seq_start.$array[4].$temp_seq_end;
				$alt_seq = $temp_seq_start.$alts[$k].$temp_seq_end;
				#print "\n";
				#print $ref_seq."\n";
				#print $alt_seq."\n";

				# here we have the window slid, then calculate the scores for each 
				my $score_ref = 0;
				my $score_alt = 0;
				
            #print $pwm_matrix[0][$j]." -- ".$pwm_matrix[3][$j]." \n ";
         
				for($j=0;$j<length($ref_seq);$j++)
				{

					my $temp_dna = substr $ref_seq, $j, 1;
               #print $j.":".$temp_dna."\n";
					if($temp_dna =~ /A/)
					{
					$score_ref = $score_ref +$pwm_matrix[0][$j]; 
					}elsif($temp_dna =~ /C/)
					{
					$score_ref = $score_ref +$pwm_matrix[1][$j];
					}elsif($temp_dna =~ /G/)
					{
					$score_ref = $score_ref +$pwm_matrix[2][$j];
					}else
					{
					$score_ref = $score_ref +$pwm_matrix[3][$j];
					}
               #print $score_ref."\n";
				}

				#print $score_ref."\n";

				for($j=0;$j<length($alt_seq);$j++)
				{
					my $temp_dna = substr $alt_seq, $j, 1;
					#print $j."\t".$temp_dna."\n";
					if($temp_dna =~ /A/)
					{
					$score_alt = $score_alt +$pwm_matrix[0][$j]; 
					}elsif($temp_dna =~ /C/)
					{
					$score_alt = $score_alt +$pwm_matrix[1][$j];
					}elsif($temp_dna =~ /G/)
					{
					$score_alt = $score_alt +$pwm_matrix[2][$j];
					}else
					{
					$score_alt = $score_alt +$pwm_matrix[3][$j];
					}
				}
       
				my $alt = $alts[$k];
	         my $relative_position = $motif_len - $i;
				my $motif = $key; 
            my $threshold = $cutoffs{$motif};
           #print $score_alt."\n";
				$score = $score."\t".nearest(0.0001,$score_ref)."/".nearest(0.0001,$score_alt);

				if( ( ( $score_ref > $threshold) | ($score_alt > $threshold) ) && ( $score_ref > 0 ) && ( $score_alt > 0 ) )
				{


					my $difference_real = nearest(0.0001,($score_alt - $score_ref) );
					my $difference_abs = abs($difference_real);
			

					print SCORE_CHANGES_FILE $motif."\t".$snp_id."\t".$ref."\t".$alt."\t".$difference_abs."\t".$difference_real."\t".$relative_position."\n";
				}
			#*********************** score *******************************

			}#for(my $i=0;$i<$motif_len;$i++)


			print OUTFILE $score."\n";

		}

	}#while(my $line1=<SEQ_FILE>)

	#$elapse = time - $start;
	#$avg = $elapse / $motif_counter; 
	#print "\t Elapse time: $elapse \t Avg time: $avg\n";

}


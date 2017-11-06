use List::MoreUtils qw(uniq);
use Math::Round;


#This code is used to obtain the SNP log ratio score

#Using a sliding window along the surrounding sequence to obtain the SNP log ratio score  
#For each SNP, calculate the log odds ratio score for each window overlapping with the SNP

#The file containing TFM-P value cutoff for the motifs corresponding to p=4e-7
my $TFM_cutoff_file = $ARGV[0];

#1st argument is the directory for position-weight matrix files
my $pwm_dir = $ARGV[1];

#2nd argument is the file containing the genomic sequence around the SNPs
my $seq_file = $ARGV[2];

#3rd argument is the output folder for the results
my $output_dir = $ARGV[3];

#4th argument is the output file for storing the results
my $score_changes_file = $ARGV[4];


open(SEQ_FILE,$seq_file) ||  die ("Cannot open input sequence file $seq_file\n"); 
open(TF_FILE,$TFM_cutoff_file) || die ("Cannot open input file $TFM_cutoff_file\n");
opendir(DIR,  $pwm_dir ) || die("Cannot open input dir: $pwm_dir\n");   
open(SCORE_CHANGES_FILE, ">$score_changes_file") || die ("Cannot write output file $score_changes_file\n"); 


#step 1: first obtain the sequence around the SNPs

#Hash array TFM-P value cutoffs
my %cutoffs = ();

#Hash array TF lengths
my %tf_len = ();

#Hash array for position weight matrices
my %tf_pwm = ();


#Read the TFM-Pvalue cutoff scores for p=4e-7
#Read the position weight matrices consequtively into hash array
while(my $line=<TF_FILE>)
{
	chomp($line);
	my @array_TFM_cutoff = split(/\s+/,$line);
	my $motif_name = $array_TFM_cutoff[0];
	my $file_path = $pwm_dir.$motif_name.'.pwm';

	#Open position matrix file to read the matrix into hash array
	unless (open (MOTIF_FILE, $file_path)){ die ("Cannot open input pwm file: $file_path\n"); }  
	my $pwm_matrix;

	while(my $line1=<MOTIF_FILE>)
	{
		chomp($line1);
		my @array = split(/\s/,$line1);
		if($array[0] =~ /A/)
		{ 
			$pwm_matrix = $array[1];
			for(my $i=2;$i<scalar(@array);$i++)
			{
				$pwm_matrix[0][$i-1] = $array[$i];
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

#For each motif position weight matrix, calculate the log odd score 
#using a sliding window that shifts one bp at each iteration.
#For the first iteration, SNP is at the end of the motif
#For the last iteration, SNP is at the beginning
#Calculate log odds score for each position
#using both reference and alternate allele for each SNP
$motif_counter = 0;
my $start = time;
foreach my $key (keys %tf_pwm)
{
	$motif_counter++;

	my $output_file = $output_dir.$key;

	open (OUTFILE, ">$output_file") || die("Cannot open output file: $output_file\n"); 

	my @pwm_matrix = $tf_pwm{$key};
	my $motif_len = $tf_len{$key};
	my $seq_file = $seq_file;

	my $value = $tf_pwm{$key};
	my @array = split(/\t/,$value); 

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
	

	#Open the file containing genomic sequences around the SNPs
	unless (open (SEQ_FILE, $seq_file)){ die ("Cannot open input sequence file $seq_file\n"); } 
	  
	while(my $line1=<SEQ_FILE>)
	{
		chomp($line1);
		my @array = split(/\t+/, $line1);
		my @alts = split(/,/, $array[5]);
      my $snp_id = $array[1];		
      my $ref = $array[4];
      my $alt = $array[5];

		#Calculate the TF binding score for reference allele and 
		#for each alternate allele (may be more than one, comma separated)
		for(my $k=0;$k<scalar(@alts);$k++)
		{
		   
			my $score = $array[0]."\t".$array[1]."\t".$array[2]."\t".$array[3]."\t".$array[4]."\t".$alts[$k]."\t".$motif_len;

			for(my $i=0;$i<$motif_len;$i++)
			{  
				my $start = 25-$motif_len+$i;

				my $ref_seq = '';
				my $alt_seq = '';
            
				my $temp_seq = substr $array[6], $start, $motif_len;

				my $temp_seq_start = substr $temp_seq, 0, $motif_len -$i-1;
				my $temp_seq_end = substr $temp_seq, $motif_len-$i, $i;

				$ref_seq = $temp_seq_start.$array[4].$temp_seq_end;
				$alt_seq = $temp_seq_start.$alts[$k].$temp_seq_end;


				# Here we have the window slid, then calculate the scores for each 
				my $score_ref = 0;
				my $score_alt = 0;
				
         	
				#Calculate the TF binding score for the reference sequence
				for($j=0;$j<length($ref_seq);$j++)
				{

					my $temp_dna = substr $ref_seq, $j, 1;

					if($temp_dna =~ /A/)
					{
						$score_ref = $score_ref + $pwm_matrix[0][$j]; 
					}elsif($temp_dna =~ /C/)
					{
						$score_ref = $score_ref + $pwm_matrix[1][$j];
					}elsif($temp_dna =~ /G/)
					{
						$score_ref = $score_ref + $pwm_matrix[2][$j];
					}else
					{
						$score_ref = $score_ref + $pwm_matrix[3][$j];
					}

				}

				#Calculate the TF binding score for the alternate sequence
				for($j=0;$j<length($alt_seq);$j++)
				{
					my $temp_dna = substr $alt_seq, $j, 1;
					if($temp_dna =~ /A/)
					{
						$score_alt = $score_alt + $pwm_matrix[0][$j]; 
					}elsif($temp_dna =~ /C/)
					{
						$score_alt = $score_alt + $pwm_matrix[1][$j];
					}elsif($temp_dna =~ /G/)
					{
						$score_alt = $score_alt + $pwm_matrix[2][$j];
					}else
					{
						$score_alt = $score_alt + $pwm_matrix[3][$j];
					}
				}
       
				my $alt = $alts[$k];
	         my $relative_position = $motif_len - $i;
				my $motif = $key; 
            my $threshold = $cutoffs{$motif};

				$score = $score."\t".nearest(0.0001,$score_ref)."/".nearest(0.0001,$score_alt);

				if( ( ( $score_ref > $threshold) | ($score_alt > $threshold) ) && ( $score_ref > 0 ) && ( $score_alt > 0 ) )
				{
					my $difference_real = nearest(0.0001,($score_alt - $score_ref) );
					my $difference_abs = abs($difference_real);	

					#Write the result into output file
					print SCORE_CHANGES_FILE $motif."\t".$snp_id."\t".$ref."\t".$alt."\t".$difference_abs."\t".$difference_real."\t".$relative_position."\n";
				}


			}#for(my $i=0;$i<$motif_len;$i++)


			print OUTFILE $score."\n";

		}

	}#while(my $line1=<SEQ_FILE>)



}


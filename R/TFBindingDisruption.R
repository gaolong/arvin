#'Extract TF name from the file name which is separated from motif name with @: AHR@M6139_1.02
#
#'@param filename file name containing the score changes between reference and alternate sequences around the input SNPs
#'@export
getTFname <- function(filename){
  input <- as.character(filename)
  tfname <- strsplit(input, "@")[[1]][1]
  return(tfname)
}


#'Calculate transcription factor disruption p values for the SNPs
#
#'@param score_changes_file file name containing the score changes between reference and alternate sequences around the input SNPs
#'@param bg_score_changes_dir directory cotaining the background score changes between the reference and alternate sequences aroung the 1000 Genome SNPs
#'@param output_filename the path to the file name to write the output
#'@export
get_disruption_p_value <- function(score_changes_file, bg_score_changes_dir, output_filename)
{
  #Read the input SNPs file into a data frame
  input_snps = read.table(file=score_changes_file,head=F,sep="\t");
  #Get the motifs
  motifs = unique(input_snps[,1])
  
  output <- data.frame()
  
  #Read the motifs one by one
  for(i in 1:length(motifs))
  {
    #Get the motif
    motif = as.character(motifs[i])
    
    #Get the SNPs overlapping with TFBS for the motif
    snp_subset = input_snps[input_snps$V1==motif,]
    
    #Get the background scores for the motif
    bg_score_changes_file = paste0(bg_score_changes_dir, motif);
    if(! file.exists(bg_score_changes_file) )
    {
      next
    }
  
    #Read the background SNPs table
    bg_snps = read.table(file=bg_score_changes_file,head=F,sep="\t",fill = T);
    bg_snps = bg_snps[complete.cases(bg_snps),]
    bg_snps_count <- nrow(bg_snps)
    
    n <- nrow(snp_subset)
    snp_subset$p_values <- c()
    
    #For each SNP with a TFBS matching the motif, calculate the p value
    for(j in 1:n)
    {
      motif <- as.character(snp_subset[j,1])
      snp <- snp_subset[j,2]
      score_change_abs <- as.numeric(snp_subset[j,5])
      score_change_signed <- as.numeric(snp_subset[j,6])
      position <- as.numeric(snp_subset[j,7])
      
      count_larger <- sum(bg_snps$V3>score_change_abs)  
      
      p_value <- count_larger / bg_snps_count
      
      snp_subset[j,8] <- p_value
      
    }	 
    
    #Add a pseudonumner for 0 p-values to avoid infinite log values
    snp_subset$V8[snp_subset$V8 == 0] <- 0.0000000000001
    q_values <- p.adjust(snp_subset$V8, method = 'BH')
    q_values[q_values==0] <- 0.000000001
    snp_subset$V9 <- q_values
    output <- rbind(output, snp_subset)
    
  }
  
  output$V10 = -log10(output$V9)
  
  tf_motif <- data.frame(do.call('rbind', strsplit(as.character(output$V1),'@',fixed=TRUE)))
  output2 <- cbind(tf_motif, output)
  output2$V1 <- NULL
  
  #Prepare the output
  colnames(output2) <- c('TF','motif', 'snp', 'ref', 'alt', 'score_change', 'score_change_abs', 'position', 'p', 'q', 'logq')
  output3 <- output2[FALSE,]
  snps <- unique(as.character(output2$snp))
  
  #If the SNP disrupts multiple motifs, select the disruption with minimum p-value 
  for(i in 1:length(snps))
  {
    snp = snps[i]
    temp = output2[output2$snp == snp, ]
    min_p_index = which.min(temp$p)
    output3 = rbind(output3, temp[min_p_index,])
  }
  
  #Select the output fields
  output4 <- output3[, c('snp', 'ref', 'alt', 'TF', 'p', 'q', 'logq')]
  colnames(output4) <- c('snp_id', 'ref', 'alt', 'TF', 'disrupion_p', 'disruption_q', 'log_disruption_q')
  
  #Write the output file
  write.table(output4, file=output_filename, row.names = F, col.names = T, quote=F, sep = "\t")

}



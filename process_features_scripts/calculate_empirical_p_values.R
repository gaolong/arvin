#count_tfs <- 1961

args<-commandArgs(TRUE);
score_changes_file <- args[1] ;
bg_score_changes_dir  <- args[2] ; 
output_filename <- args[3] ; 


#score_changes_file <- "/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.2/meta_data/gold_esnps//disruption//TFBS_score_changes.txt" ;
#bg_score_changes_dir  <-  "/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/source_data//score_changes_sampling/" ; 
#output_filename <-  "/mnt/isilon/cbmi/tan_lab/uzuny/arvin/revision/software/v0.2/meta_data/gold_esnps//disruption//disruption_p.txt" ; 


getTFname <- function(filename){
  input <- as.character(filename)
  tfname <- strsplit(input, "@")[[1]][1]
  return(tfname)
}

input_snps = read.table(file=score_changes_file,head=F,sep="\t");
#print(input_snps[1:5,])
#print(dim(input_snps))
motifs = unique(input_snps[,1])

#bg_snp_scores <- bg_snps$V4

#filecon = file(pvalue_filename,"w")

output <- data.frame()

for(i in 1:length(motifs))
  #for(i in 1:dim(input_snps)[1])
{
  motif = as.character(motifs[i])
  snp_subset = input_snps[input_snps$V1==motif,]
  bg_score_changes_file = paste0(bg_score_changes_dir, motif);
  if(! file.exists(bg_score_changes_file) )
  {
    next
  }
  #cat(paste(as.character(i), motif) ) 
  #cat('\n')
  bg_snps = read.table(file=bg_score_changes_file,head=F,sep="\t",fill = T);
  bg_snps = bg_snps[complete.cases(bg_snps),]
  #print(dim(bg_snps))
  #print(bg_snps[1:5,])
  bg_snps_count <- nrow(bg_snps)
  
  n <- nrow(snp_subset)
  snp_subset$p_values <- c()
  
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

colnames(output2) <- c('TF','motif', 'snp', 'ref', 'alt', 'score_change', 'score_change_abs', 'position', 'p', 'q', 'logq')
output3 <- output2[FALSE,]
snps <- unique(as.character(output2$snp))

for(i in 1:length(snps))
{
  snp = snps[i]
  temp = output2[output2$snp == snp, ]
  min_p_index = which.min(temp$p)
  output3 = rbind(output3, temp[min_p_index,])
}

output4 <- output3[, c('snp', 'ref', 'alt', 'TF', 'p', 'q', 'logq')]
colnames(output4) <- c('snp_id', 'ref', 'alt', 'TF', 'disrupion_p', 'disruption_q', 'log_disruption_q')
write.table(output4, file=output_filename, row.names = F, col.names = T, quote=F, sep = "\t")

#save.image("/home/uzuny/arvin/revision/software/v0.2/meta_data/gold_esnps/disruption/disruption_p.RData")



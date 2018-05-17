
args <- commandArgs(trailingOnly = TRUE)

#Input DGE file
input_file  <- args[1]
mapping_file <- args[2]
output_file <- args[3]

#input_file  = "/home/uzuny/t1d/apply_arvin/data/TR/de.txt"
#mapping_file = "/home/uzuny/arvin/arvin_package/v1.0/arvin_annotation_data/entrez_ids.txt"
#output_file = "/home/uzuny/t1d/apply_arvin/data/TH/de.normalized.entrez_id.txt"

input = read.table(input_file, sep = "\t")
head(input)
zeros = which(input[,2] == 0)
ones = which(input[,2] == 1)

head(zeros)
head(ones)

df1 = input
df1[zeros, 2] = 0.000000001
df1[ones, 2]  = 0.999999999
head(input)
head(df1)
colnames(df1) = c('gene_symbol', 'p_value')

gene_map = read.table(file = mapping_file, head=T);
dim(gene_map)
head(gene_map)

df2 = merge(df1, gene_map, by.x = "gene_symbol", by.y = "gene_symbol")
head(df2)

log_p = -log10(df2$p_value)
max_log_p = max(log_p)
min_log_p = min(log_p)

normalized_log_p = (log_p - min_log_p) / (max_log_p - min_log_p)
max(normalized_log_p)
min(normalized_log_p)

df3 = data.frame(df2$entrez_id, normalized_log_p)
head(df3)
dim(df3)

colnames(df3) = c('entrez_id',  'normalized_log_de_p_value')
write.table(df3, file=output_file, sep = '\t', quote = F, col.names = T, row.names = F)



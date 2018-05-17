library("igraph")
library("ARVIN")
library("stringr")

args <- commandArgs(T)

opt_classification_file = args[1]
cv_classification_file = args[2]
disruption_file = args[3]
features_file = args[4]
edge_file = args[5]
node_file = args[6]
output_score_file = args[7]


print("Loading the model...")
load(opt_classification_file)
load(cv_classification_file)
###################Treg

set.seed(123)

print("Building the network...")
Net <- makeNet(edge_file, node_file)

print("Computing the network features...")
topoFeature <-NetFeature(Net, node_file, edge_file, disruption_file)
colnames(topoFeature)[6] <- "DE"
head(topoFeature)

print("Adding the sequence features...")
gf <- read.table(features_file, header=TRUE)
row.names(gf) <- gf$snp_id
features <- cbind(topoFeature, gf[row.names(topoFeature), 2:dim(gf)[2]])
head(features)
features <- subset(features, !is.na(rowMeans(features)))
head(features)

print("Computing ARVIN scores...")
prob <- predMod(features, rfFit_Com_TSS_opt)
output <- prob[,1]
names(output) <- row.names(prob)

print("Writing the output...")
write.table(output, output_score_file, quote=F, col.names=F, sep="\t")
print("Done ...")


print("ARVIN prediction scores are in here:")
print(output_score_file)

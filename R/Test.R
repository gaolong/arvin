#library(arvin)
#source("R/NetConfig.R")
#source("R/ModSearch.R")
#source("R/NetworkFeatureCalculation.R")
edgeFile <- "example_input/EdgeFile.txt"
nodeFile <- "example_input/NodeFile.txt"

NetComp <- function(edgeFile, nodeFile){
Net <- makeNet(edgeFile, nodeFile)
node_data <- read.table(nodeFile)
node_data <- subset(node_data, node_data$V3 == "eSNP")
eSNP_seeds <- as.character(node_data$V1)#use all SNPs as seeds for search 
node_data <- NULL#free memory

V_weight <- V(Net)$weight#put the node weight into an array to save time
names(V_weight) <- V(Net)$name#set the names for the array
E_adj <- as_adj(Net,attr="weight")#use adj matrix to save time for getting edge weight
colnames(E_adj) <- names(V_weight)
row.names(E_adj) <- names(V_weight)

Adj_List <- c()
batch <- 1000#split the matrix into a few parts in case of a memory error
index <- floor(dim(E_adj)[1]/batch)
for(i in 1:index){
  ifelse(i != index,
         end <- i * batch,
         end <- dim(E_adj)[1]
  )
  start <- (i-1) * batch + 1
  cur_list <- apply(E_adj[start:end,], 1, function(x) x[x!=0])
  Adj_List <- append(Adj_List, cur_list)
  
}

mod_sets <- list()
for(i in 1:length(eSNP_seeds)){
  temp_seed <- eSNP_seeds[i]
  temp_mod <- module_search_ind(V_weight, Adj_List, temp_seed, E_adj)
  mod_sets <- append(mod_sets, list(temp_mod))
}

#####################################merge identified modules
merged_modules <- merge_modules(mod_sets, 0.1,  V_weight, Adj_List)
}
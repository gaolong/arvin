#'Compute network-based features needed for risk variant prediction
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes the path for the node file 
#'@param edgeFile the path for the network file
#'@return a matrix containing the different type network features 
#'@export
NetFeature <- function(Net, nodeFile, edgeFile, snpFile){
  #edgeFile <- "example_input/EdgeFile.txt"
  #nodeFile <- "example_input/NodeFile.txt"
  
  ########################################preprocess the data
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
  #############################################################
  
  Net <- makeNet(edgeFile, nodeFile)
  edge_data <- read.table(edgeFile)
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  snp_match <- Nodes
  names(snp_match) <- snps
  bet_vals <- BetFeature(Net, edge_data)
  close_vals <- CloseFeature(Net, edge_data)
  page_vals <- PageFeature(Net, edge_data)
  wd_vals <- WDFeature(Adj_List, edge_data)
  snpDist <- snpFeature(snpFile, eSNP_seeds)
  mod_vals <- ModuleFeature(Adj_List, E_adj, eSNP_seeds, V_weight, Nodes)
  
  FeatureMatrix <- data.frame(mod_vals[eSNP_seeds], bet_vals[eSNP_seeds], 
                              close_vals[eSNP_seeds], page_vals[eSNP_seeds],
                              wd_vals[eSNP_seeds], snpDist[eSNP_seeds])
  row.names(FeatureMatrix) <- eSNP_seeds
  colnames(FeatureMatrix) <- c("Module_Score", "Betweenness", "Closeness","Pagerank", "Weighted_Degree","SNP_Disrutption")
  return(FeatureMatrix)
}

#'Compute betweenness centrality
#'
#'This function computes the betweenness centrality for a list of candidate snps
#'
#'@param Net a graph object representing the input network
#'@param edge_data a dataframe specifying the edge information
#'@return a vector of betweenness centrality values
#'@export
BetFeature <- function(Net, edge_data){
  nm <- length(V(Net))
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  bet_vals <- estimate_betweenness(Net, Nodes, directed = FALSE, cutoff=5)
  unhit <- setdiff(Nodes, names(bet_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_bet_vals <- c(bet_vals, unhit_val)
  names(c_bet_vals) <- edge_data[,1]
  return(c_bet_vals/nm)
}

#'Compute closeness centrality
#'
#'This function computes the closeness centrality for a list of candidate snps
#'
#'@param Net a graph object representing the input network
#'@param edge_data a dataframe specifying the edge information
#'@return a vector of closeness centrality values
#'@export
CloseFeature <- function(Net, edge_data){
  nm <- length(V(Net))
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  close_vals <- estimate_closeness(Net, Nodes, cutoff=5)
  unhit <- setdiff(Nodes, names(close_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_close_vals <- c(close_vals, unhit_val)
  names(c_close_vals) <- edge_data[,1]
  return(c_close_vals*nm)
}

#'Compute pagerank centrality 
#'
#'This function computes the pagerank centrality for a list of candidate snps
#'
#'@param Net a graph object representing the input network
#'@param edge_data a dataframe specifying the edge information
#'@return a vector of pagerank centrality values
#'@export
PageFeature <- function(Net, edge_data){
  nm <- length(V(Net))
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  page_vals <- page_rank(Net, vids=Nodes)$vector
  unhit <- setdiff(Nodes, names(page_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_page_vals <- c(page_vals, unhit_val)
  names(c_page_vals) <- edge_data[,1]
  return(c_page_vals*nm)
}

#'Compute TF binding disruption score
#'
#'This function computes the TF binding disruption score for a list of candidate snps
#'
#'@param SNP_file the path of the file with SNP TF disruption score information
#'@param eSNP_seeds a vector of candidate eSNPs
#'@return a vector of TF binding disruption scores
#'@export
snpFeature <- function(SNP_file, eSNP_seeds){
  snp_tmp <- read.table(SNP_file)
  snp_vals <- snp_tmp[,2]
  names(snp_vals) <- snp_tmp[,1]
  return(snp_vals[eSNP_seeds])
}

#'Compute weighted degree
#'
#'This function computes the weighted degree for a list of candidate snps
#'
#'@param Net a graph object representing the input network
#'@param edge_data a dataframe specifying the edge information
#'@return a vector of weighted degree values
#'@export
WDFeature <- function(Adj_List, edge_data){
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  wd_vals <- unlist(lapply(Adj_List[Nodes], sum))
  unhit <- setdiff(Nodes, names(wd_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_wd_vals <- c(wd_vals, unhit_val)
  names(c_wd_vals) <- edge_data[,1]
  return(c_wd_vals)
}

#'Compute module scores
#'
#'This function computes the module scores for a list of candidate snps
#'
#'@param Net a graph object representing the input network
#'@param edge_data a dataframe specifying the edge information
#'@return a vector of module scores
#'@export
ModuleFeature <- function(Adj_List, E_adj, eSNP_seeds, V_weight, Nodes){
  mod_sets <- list()
  for(i in 1:length(eSNP_seeds)){
    temp_seed <- eSNP_seeds[i]
    temp_mod <- module_search_ind(V_weight, Adj_List, temp_seed, E_adj)
    mod_sets <- append(mod_sets, list(temp_mod))
  }
  
  #####################################merge identified modules
  merged_modules <- merge_modules(mod_sets, 0.1,  V_weight, Adj_List)
  snp_vec <- c()
  score_vec <- c()
  for(i in 1:length(merged_modules)){
    for(j in 1:length(merged_modules[[i]]$mem)){
      cur_mem <- merged_modules[[i]]$mem[j]
      if(cur_mem %in% eSNP_seeds & !cur_mem %in% snp_vec){
        snp_vec <- c(snp_vec, cur_mem)
        score_vec <- c(score_vec, merged_modules[[i]]$avg_score)
      }
    }
  }
  mod_val <- score_vec
  names(mod_val) <- snp_vec
  unhit <- setdiff(eSNP_seeds, snp_vec)#for snps which are not in any modules assign module score 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  mod_val <- c(mod_val, unhit_val)
  return(mod_val)
}
#'Compute network-based features needed for risk variant prediction
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 
NetFeature <- function(Net, nodeFile, edgeFile){
  edgeFile <- "example_input/EdgeFile.txt"
  nodeFile <- "example_input/NodeFile.txt"
  
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
  snpDist <- snpFeature(Net, eSNP_seeds)
  mod_vals <- ModuleFeature(Adj_List, E_adj, eSNP_seeds, V_weight, Nodes)
  
  #names(bet_vals) <- names(snp_match)
  #names(close_vals) <- names(snp_match)
  #names(page_vals) <- names(snp_match)
  #names(wd_vals) <- names(snp_match)
  
  FeatureMatrix <- data.frame(mod_vals[eSNP_seeds], bet_vals[eSNP_seeds], 
                              close_vals[eSNP_seeds], page_vals[eSNP_seeds],
                              wd_vals[eSNP_seeds], snpDist[eSNP_seeds])
  row.names(FeatureMatrix) <- eSNP_seeds
  colnames(FeatureMatrix) <- c("Module_Score", "Betweenness", "Closeness","Pagerank", "Weighted_Degree","SNP_Disrutption")
  return(FeatureMatrix)
}

#'Compute betweenness centrality for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 
BetFeature <- function(Net, edge_data){
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  bet_vals <- estimate_betweenness(Net, Nodes, directed = FALSE, cutoff=5)
  unhit <- setdiff(Nodes, names(bet_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_bet_vals <- c(bet_vals, unhit_val)
  names(c_bet_vals) <- edge_data[,1]
  return(c_bet_vals)
}

#'Compute closeness centrality for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 
CloseFeature <- function(Net, edge_data){
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  close_vals <- estimate_closeness(Net, Nodes, cutoff=5)
  unhit <- setdiff(Nodes, names(close_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_close_vals <- c(close_vals, unhit_val)
  names(c_close_vals) <- edge_data[,1]
  return(c_close_vals)
}

#'Compute pagerank centrality for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 
PageFeature <- function(Net, edge_data){
  edge_data <- subset(edge_data, edge_data$V4 == "EP")
  snps <- as.character(edge_data[,1])
  Nodes <- as.character(edge_data[,2])
  page_vals <- page_rank(Net, vids=Nodes)$vector
  unhit <- setdiff(Nodes, names(page_vals))#for snps which are not in any modules assign 0 to them
  unhit_val <- rep(0, length(unhit))
  names(unhit_val) <- unhit
  c_page_vals <- c(page_vals, unhit_val)
  names(c_page_vals) <- edge_data[,1]
  return(c_page_vals)
}

#'Compute module score for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 
#'

#'Compute pagerank centrality for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features
snpFeature <- function(Net, eSNP_seeds){
  set.seed(1)
  snp_vals <- runif(length(eSNP_seeds), 0, 1)
  names(snp_vals) <- eSNP_seeds
  return(snp_vals)
}

#'Compute weighted degree for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 

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

#'Compute module score for a list of nodes
#'
#'This function ...
#'
#'@param Net a graph object representing the input network
#'@param Nodes a list of node names 
#'@return a list object containing the different type network features 

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
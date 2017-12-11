library("igraph")#igraph contains basic graph operations and some simple graph algorithms
#library("Rcpp")

##################################################################################
#given a list of modules
#this function returns a matrix for pairwise jaccard index
#cppFunction('NumericMatrix pair_jac(List mod_list){
#            int num = mod_list.size();
#            NumericMatrix mat(num,num);
#            for(int i = 0; i < num - 1; ++i){
#            CharacterVector v1 = mod_list[i];
#            for(int j = i + 1; j < num; ++j){
#            CharacterVector v2 = mod_list[j];
#            IntegerVector v3 = match(v1,v2);
#            int n = v3.size();
#            int nx = v1.size();
#            int ny = v2.size();
#            int count = 0;
#            for (int k = 0; k < n; ++k){
#            count += 1 - IntegerVector::is_na(v3[k]);
#            }
#            float jac = (float)(count)/(float)(nx + ny - count);
#            mat(i,j) = jac;
#            mat(j,i) = jac;
#            }
#            }
#            return mat;
#            }')

#'Search modules from a list of seeds/nodes
#'
#'This function searches attached modules based on the list of seeds
#'
#'@param Net a graph object representing the input network
#'@param V_weight a numeric vector of node weight
#'@param Adj_List a list object representing the adjacency list of the input network
#'@param E_adj the adjacency matrix
#'@return a list object containing modules and related statistics
#'@export
module_search_ind <- function(V_weight, Adj_List, ind_seed, E_adj){
  ptm_1 <- proc.time()
  cur_seed <- ind_seed#get the seed name
  #print("current seed:")
  #print(cur_seed)
  cur_mod <- list(mem=c(cur_seed), seed=cur_seed, score=V_weight[cur_seed], last_mem=cur_seed, neighbor=c(), score_table=c(), size=1, pval=-1, fdr=-1, state="NULL", com_size=0)#initialize a module object in list format
  #last_mem is the last gene which was added into current module
  #neighbor is the list of all possible neighbor genes of members in current module
  #score table records neighbors that will be potentially added into current module and the second column is the gain of score if add this node
  while(cur_mod$state == "NULL"){#run thAdj_List_randis loop until there's no node can be added into current module
    #print(cur_mod$mem)
    cand_neighbor <- c()#get all neighbors of genes in current module
    #iteratively obtain neighbors for all genes
    
    cur_mem <- cur_mod$last_mem#obtain a given node
    cur_neighbor <- names(Adj_List[[cur_mem]])
    cand_neighbor <- union(cur_mod$neighbor, cur_neighbor)#add newly identified neighbors to the entire set	
    cand_neighbor <- setdiff(cand_neighbor, cur_mod$mem)#get rid of nodes that are already in this module
    cand_neighbor <- setdiff(cand_neighbor, names(cur_mod$score_table))#also get rid of potential neighbors that have been considered
    cur_mod$neighbor <- cand_neighbor
    
    if(length(cand_neighbor) > 5){
      #if the potential set is not empty, at most pick top 5
      NW_cand <- V_weight[cand_neighbor]#get the node weight of candidate nodes
      NW_cand <- sort(NW_cand, decreasing=T)#rank the array by node weight
      cand_neighbor <- names(NW_cand[1:5])#only use the at most top 5 nodes
    }
    cur_mod <- module_update(V_weight, Adj_List, cur_mod, cand_neighbor, E_adj)
    cur_mod$com_size <- get_com_size(cur_mod$mem, V_weight, Adj_List)
    if(cur_mod$size > 100){
      #if(cur_mod$com_size > 1000){
      break
    }
  }
  
  #print(cur_mod$mem)
  #print("Time to find this module:")
  #print(proc.time() - ptm_1)
  return(cur_mod)
  
}

##################################################################################
#This function was designed to do module membership checking
#Net, this is an igraph obeject with node and edge weight specified
module_update <- function(V_weight, Adj_List, cur_mod, can_neighbor, E_adj){
  new_mod <- cur_mod#updated module, but now it is just a copy of current module
  #each iteration the following loop tries to update new_mod
  neighbor_neighbor <- c()#neighbors of a given neighbor
  if(length(can_neighbor) != 0){#if this set of neighbors is not empty
    #score list for candidate genes only
    add_score_list <- rep(-1, length(can_neighbor))
    for(i in 1:length(can_neighbor)){
      cur_nb <- can_neighbor[i]#get a given neighbor
      neighbor_neighbor <- names(Adj_List[[cur_nb]])
      #print(proc.time() - ptm_1)
      con_mem <- intersect(cur_mod$mem, neighbor_neighbor)#to get overlap with module genes so we can know which ones are connected to this neighbor		
      gain <- V_weight[cur_nb] + sum(Adj_List[[cur_nb]][con_mem])
      add_score_list[i] <- gain
    }
    names(add_score_list) <- can_neighbor
    cur_mod$score_table <- c(cur_mod$score_table, add_score_list)
  }
  #print("length")
  #print(length(cur_mod$score_table))
  if(length(cur_mod$score_table) > 0){
    max_gain <- max(cur_mod$score_table)
    max_index <- which.max(cur_mod$score_table)#pick the node which can introduce largest score gain
    new_node <- names(cur_mod$score_table[max_index])#get the name of the node
    #if(max_gain  >= 0.01 * cur_mod$score){
    #if(max_gain  > 0){
    #if(T){
    hit_GO <- FALSE
    if(cur_mod$size < 20){
      hit_GO <- TRUE
    }
    else{
      cur_overlap <- intersect(names(Adj_List[[new_node]]), cur_mod$mem)
      if((cur_mod$score + max_gain)/(cur_mod$com_size + length(cur_overlap) + 1)  >=   cur_mod$score/cur_mod$com_size | cur_mod$size > 50){
        hit_GO <- TRUE
      }
    }
    if(hit_GO){
      #max_index <- which.max(cur_mod$score_table)#pick the node which can introduce largest score gain
      #new_node <- names(cur_mod$score_table[max_index])#get the name of the node
      max_gain <- cur_mod$score_table[max_index]#get the gain of score
      cur_mod$mem <- c(cur_mod$mem,new_node)#update the membership
      cur_mod$size <- cur_mod$size + 1#update the size
      cur_mod$score <- cur_mod$score + max_gain#update the score
      cur_mod$last_mem <- new_node#indicate which has been added
      cur_mod$score_table <- cur_mod$score_table[-max_index]#remove that added node
      #update the score table with regard to the new node since the topology of the module has changed
      neighbor_neighbor <- names(Adj_List[[new_node]])
      poten_overlap <- intersect(names(cur_mod$score_table), neighbor_neighbor)#check how many potential nodes that are connected to this newly added node, and then update their corresponding score
      if(length(poten_overlap) > 0){
        cur_mod$score_table[poten_overlap] <-  cur_mod$score_table[poten_overlap] + E_adj[new_node, poten_overlap]
      }
    }#if this module can be still updated
    else{#there is no neighbor can satistied the criterion
      cur_mod$state <- "Stop"#set the state flag to "stop"
    }
    if(length(can_neighbor) == 0 & length(cur_mod$score_table) == 0 ){
      #there is no additional nodes that are connected to current module
      cur_mod$state <- "Stop"#set the state flag to "stop"
    }
  }#check the score table length
  else{#if the score table is already empty
    cur_mod$state <- "Stop"
  }
  #print("Time for update all!")
  #print(proc.time() - ptm_1)
  return(cur_mod)
}





###################R code for computing pairwise overlap between modules
pair_jac <- function(mod_mem){
  names(mod_mem) <- as.character(1:length(mod_mem))
  nms <- combn(names(mod_mem) , 2 , FUN = paste0 , collapse = "" , simplify = FALSE )
  # Make the combinations of list elements
  ll <-  combn( mod_mem , 2 , simplify = FALSE )
  jac <- unlist(lapply( ll , function(x){
    l1 <- length(x[[1]])
    l2 <- length(x[[2]])
    overlap <- sum(x[[1]] %in% x[[2]])
    overlap/(l1 + l2 - overlap)
    }))
  len <- length(mod_mem)
  jac_matrix <- matrix(rep(1, len * len), nrow=len)
  count <- 1
  for(i in 1:(len-1)){
    for(j in (i + 1):len){
      cur_val <- jac[count]
      count <- count + 1
      jac_matrix[i,j] <- cur_val
      jac_matrix[j,i] <- cur_val
    }
  }
  return(jac_matrix)
}

##################################################################################
#cluster modules based on pairwise jaccard index
merge_modules <- function(mod_sets, tree_height=0.5, V_weight, Adj_List){
  mod_mem <- lapply(mod_sets, function(x) x[["mem"]])#get the list of member names only
  #print("Time for running pairwise jac index: ")
  system.time(pj <- pair_jac(mod_mem))#compute pairwise jaccard index
  #set.seed(1)
  #m_num <- length(mod_mem)
  #pj <- matrix(runif(m_num * m_num, 0, 1), nrow=m_num)
  #print(dim(pj))
  #print(mod_sets)
  #print("Time for clustering modules: ")
  system.time(hc <- hclust(as.dist(1-pj)))#cluster all modules based on pairwise overlap
  mem_clr <- cutree(hc, h = tree_height)#cut the tree by the height parameter
  #mem_clr <- cutree(hc, k = length(mod_mem))#cut the tree by the height parameter
  names(mem_clr) <- 1:length(mem_clr)
  num <- max(mem_clr)#get how many modules that are still remaining
  ptm_1 <- proc.time() 
  up_mod <- list()#list of modules after merging
  for(i in 1:num){
    cur_mem_list <- mod_mem[which(mem_clr==i)]
    cur_mem <- unique(unlist(cur_mem_list))#merge all modules with relatively large overlap
    #initialize a module
    cur_mod <- list(mem=cur_mem, score=0,avg_score=0, last_mem="Merged", neighbor=c(), score_table=c(), size=length(cur_mem), pval=-1, fdr=-1, state="NULL", gene_num=0, com_size=0, com_score=0, gscore=0, gcomscore=0)
    cur_mod$score <- compute_mod_score(cur_mem, V_weight, Adj_List)#compute the score for a given module
    cur_mod$avg_score <- cur_mod$score/cur_mod$size
    cur_mod$gene_num <- length(cur_mem)
    num_snp <- cur_mod$size - cur_mod$gene_num
    if(cur_mod$gene_num > 5){
      up_mod <- append(up_mod, list(cur_mod))#append a new module to the updated module list
    }	
  }#for
  
  #print("Finish module merging!")
  return(up_mod)
}

##################################################################################
#compute the score for a module given the set of member names, node weight and adj list
compute_mod_score <- function(mem, V_weight, Adj_List){
  score <- 0
  score <- score + sum(V_weight[mem])#get all node_score
  Sub_adj_list <- Adj_List[mem]#get subset of the adj list only for the module members only
  edge_score <- 0
  for(i in 1:length(mem)){
    this_mem <- mem[i]
    common_mem <- intersect(mem, names(Sub_adj_list[[this_mem]]))
    edge_score <- edge_score + sum(Sub_adj_list[[this_mem]][common_mem])
  }
  score <- score + edge_score/2#each edge has been counted twice
  return(score)
}

##################################################################################
#get combine size which is the number of nodes plus number of edges
get_com_size <- function(mem, V_weight, Adj_List){
  com_size <- 0
  com_size <- length(mem)#get the size of nodes
  edge_score <- 0
  Sub_adj_list <- Adj_List[mem]#get subset of the adj list only for the module members only
  num_edge <- 0
  for(i in 1:length(mem)){
    this_mem <- mem[i]
    common_mem <- intersect(mem, names(Sub_adj_list[[this_mem]]))
    num_edge <- num_edge + length(common_mem)
  }
  com_size <- com_size + num_edge/2#each edge has been counted twice and it is shared by two nodes
  return(com_size)
}



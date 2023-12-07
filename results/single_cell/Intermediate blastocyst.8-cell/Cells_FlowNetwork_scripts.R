library(igraph); 
subnet <- read.table('/nfs/team205/jb62/other/Stip1-in-early-development/results/single_cell//Intermediate blastocyst.8-cell/Cells_FlowNetwork_all_paths_subnet.txt', sep='	', header = TRUE);
totalFlow <- read.table('/nfs/team205/jb62/other/Stip1-in-early-development/results/single_cell//Intermediate blastocyst.8-cell/Cells_FlowNetwork_all_paths_totalFlow.txt', sep='	', header = TRUE);
colnames(subnet) <- c('from', 'to', 'weight');
iG <- graph.data.frame(subnet, directed=FALSE);
V(iG)$label <- V(iG)$name;
V(iG)[as.vector(totalFlow$cell)]$totalFlows <- totalFlow$flows;
V(iG)$degree <- igraph::degree(iG);
neighS <- V(induced.subgraph(graph=iG,vids=unlist(neighborhood(graph=iG,order=1,nodes='s'))))$name;
neighT <- V(induced.subgraph(graph=iG,vids=unlist(neighborhood(graph=iG,order=1,nodes='t'))))$name;
V(iG)$props = rep("INTER", length(V(iG)$name));
V(iG)[match(neighS, V(iG)$name)]$props = "SOURCE";
V(iG)[match(neighT, V(iG)$name)]$props = "TARGET";
iG2 <- delete.vertices(iG, c('s','t'));
write.graph(iG2,'/nfs/team205/jb62/other/Stip1-in-early-development/results/single_cell//Intermediate blastocyst.8-cell/Cells_FlowNetwork_all_paths_subnet.gml', format = 'gml');
#--------------------------------------------------------------------------

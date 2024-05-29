# Install and load the CodeDepends package
library(CodeDepends)
library(igraph)
library(graph)
library(Rgraphviz)


# Example script content
sc <- readScript("~/test_script.R")

#plot graph
g <- makeVariableGraph(info = getInputs(sc))
callg <- makeTaskGraph(info = getInputs(sc))


# Extract nodes
nodes <- nodes(g)

# Manually extract edges from the edgeL slot
edges_list <- lapply(nodes, function(node) {
  node_edges <- edgeL(g)[[node]]
  if (!is.null(node_edges)) {
    to_nodes <- node_edges$edges
    if (length(to_nodes) > 0) {
      return(data.frame(from = rep(node, length(to_nodes)), to = nodes[to_nodes]))
    }
  }
})

# Clean the list to remove NULL entries
edges_list <- edges_list[!sapply(edges_list, is.null)]

# Combine all edges into a single data frame
edges_df <- do.call(rbind, edges_list)

# Create an igraph graph from these nodes and edges
ig <- graph_from_data_frame(edges_df, directed = TRUE, vertices = data.frame(name = nodes))

# Now you can use igraph's functions to plot and analyze the graph
plot(ig)


#alternate graph

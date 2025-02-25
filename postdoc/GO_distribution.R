# Load necessary libraries
library(data.table)

# Path to the downloaded GAF file
gaf_file <- 
gaf <- fread("/Users/nikgvetr/Downloads/goa_human.gaf", skip = "!")[, .(gene_id = V2, go_id = V5)]
gaf <- gaf[!duplicated(gaf),]
gaf_gene <- sapply(split(gaf$go_id, gaf$gene_id), length)
gaf_go <- sapply(split(gaf$gene_id, gaf$go_id), length)
hist(gaf_go)
hist(gaf_gene, breaks = 100)

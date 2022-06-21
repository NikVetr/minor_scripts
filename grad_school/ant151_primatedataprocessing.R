d <- d_orig <- readLines("data/10kTrees_finalFile_version3.nex")
splits <- c(7,which(d == ""))
d <- lapply(1:(length(splits)-1), function(s) d[(splits[s]+1):(splits[s+1]-1)])
taxa <- c("Homo_sapiens", "Homo_sapiens_neanderthalensis", "Pan_troglodytes_troglodytes", 
          "Gorilla_gorilla_gorilla", "Papio_hamadryas", "Cebus_albifrons", "Tarsius_syrichta", "Lemur_catta", 
          "Nycticebus_bengalensis", "Galeopterus_variegatus", "Mandrillus_sphinx", "Pan_paniscus", "Pongo_pygmaeus")
d <- lapply(1:17, function(gene) d[[gene]][is.element(sapply(1:302, function(taxon) (strsplit(d[[gene]][taxon], " ")[[1]][1])), taxa)])
ntaxa <- length(taxa)
has_seqs <- sapply(1:17, function(gene) (setdiff(1:ntaxa, 
                                                 grep(d[[gene]], pattern = paste0(rep("N", 500), collapse = "")))))
sapply(1:length(has_seqs), function(gene) length(has_seqs[[gene]]))
sapply(1:length(has_seqs), function(gene) setdiff(1:ntaxa,has_seqs[[gene]]))

#genes 2 and 7 cover all taxa
genes <- c(2,7)
tot_nchar <- nchar(strsplit(d[[genes[1]]][1], " ")[[1]][2]) + nchar(strsplit(d[[genes[2]]][1], " ")[[1]][2])
d_orig[4] <- gsub(x = d_orig[4], pattern = "17972", replacement = tot_nchar)
newLines <- list(d_orig[1:7], 
                 lapply(1:length(genes), function(gene) d_orig[(splits[genes[gene]]+1):(splits[genes[gene]+1]-1)]  ), 
                 d_orig[length(d_orig)], "\nEND;")
filename <- paste0("data/10kTrees_finalFile_version3_genes", paste0(genes, collapse = "+"), ".nex")
writeLines(text = unlist(newLines), con = filename)
write.nexus.data(read.nexus.data(filename), interleaved = F, file = filename)

plot(phangorn::maxCladeCred(read.tree("~/data/JC.trees")))

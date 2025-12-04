#set top level parameters
n_total_pathways <- 2E3
n_total_genes <- 2E4
n_genes_in_pathways <- sample(5:50, n_total_pathways, replace = T)
n_DE_genes <- 50

#simulate independent DE genes and pathway genes
genes_in_pathways <- lapply(n_genes_in_pathways, function(n_genes){
  factor(sample(1:n_total_genes, n_genes, replace = T),
         levels = 1:n_total_genes)
})
DE_genes <- factor(x = sample(n_total_genes, size = n_DE_genes, replace = F),
                   levels = 1:n_total_genes)

#run hypergeometric tests
results <- do.call(rbind, lapply(1:n_total_pathways, function(pathway_i){
  ct_table <- matrix(c(sum(DE_genes %in% genes_in_pathways[[pathway_i]]), 
                       sum(!DE_genes %in% genes_in_pathways[[pathway_i]]), 
                       sum(!(1:n_total_genes %in% DE_genes) & (1:n_total_genes %in% genes_in_pathways[[pathway_i]])), 
                       sum(!(1:n_total_genes %in% DE_genes) & !(1:n_total_genes %in% genes_in_pathways[[pathway_i]]))), 
                     nrow=2, 
                     dimnames=list(Pathway=c("in","out"),DE=c("DE","not_DE")))
  fisher_out <- fisher.test(ct_table, alternative = "greater")
  return(data.frame(intersect_size = ct_table[1,1],
                    hyper_pval = fisher_out$p.value,
                    OR = fisher_out$estimate))
}))

#subset to non-empty intersects
results_sub <- results[results$intersect_size > 0,]

#post-hoc multiplicity adjustment
results$adj_pval <- p.adjust(results$hyper_pval, method = "BH")
results_sub$adj_pval <- p.adjust(results_sub$hyper_pval, method = "BH")

#visualize
p_breaks <- 0:100/100
hist(results$hyper_pval, breaks = p_breaks)
hist(results$adj_pval, breaks = p_breaks)

#subset to non-empty intersects
hist(results_sub$hyper_pval, breaks = p_breaks)
hist(results_sub$adj_pval, breaks = p_breaks)

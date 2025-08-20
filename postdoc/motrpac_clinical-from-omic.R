library(tidyverse)
library(tidymodels)
library(glmnet)
library(MotrpacHumanPreSuspensionData)
library(data.table)

# load molecular data
data("cln_raw_cardiopulmonary_ex_test")
molecular_data <- load_qc()
str_lol <- function(x, level = 0, n_names = Inf) {
  indent <- function(lvl) paste(rep("  ", lvl), collapse = "")
  
  format_names <- function(nms, max_n, indent_str) {
    if (is.null(nms)) return("unnamed")
    n_total <- length(nms)
    if (n_total <= max_n) {
      return(paste(nms, collapse = ", "))
    } else {
      shown <- nms[1:max_n]
      omitted <- n_total - max_n
      return(paste(paste(shown, collapse = ", "), sprintf("... [+ %d names]", omitted)))
    }
  }
  
  recurse <- function(x, level, path = list()) {
    if (!is.list(x)) {
      if (length(path) < 2) {
        stop("List is not deep enough to go up two levels.")
      }
      
      lvl_up <- level - 2
      container <- path[[length(path) - 1]]
      nms <- names(container)
      
      cat(indent(lvl_up), "== Summary two levels up ==\n")
      if (is.null(nms)) {
        cat(indent(lvl_up), sprintf("List at level %d: unnamed list of length %d\n", lvl_up + 1, length(container)))
      } else {
        name_str <- format_names(nms, n_names, indent(lvl_up))
        cat(indent(lvl_up), sprintf("List at level %d: named list of length %d (from element [[1]]) with names: %s\n",
                                    lvl_up + 1, length(nms), name_str))
      }
      
      types <- sapply(container, function(e) class(e)[1])
      cat(indent(lvl_up), "  Classes of elements:", paste(unique(types), collapse = ", "), "\n")
      
      dims <- sapply(container, function(e) {
        d <- dim(e)
        if (is.null(d)) "" else paste(d, collapse = "x")
      })
      dims <- dims[nzchar(dims)]
      if (length(dims)) {
        cat(indent(lvl_up), "  Dimensions of elements:", paste(unique(dims), collapse = ", "), "\n")
      }
      
      cat(indent(lvl_up), "  Example str() of element [[1]]:\n")
      str_out <- capture.output(str(container[[1]]))
      cat(paste0(indent(lvl_up), "  ", str_out, collapse = "\n"), "\n")
      return(invisible(NULL))
    }
    
    nms <- names(x)
    if (is.null(nms)) {
      cat(indent(level), sprintf("Level %d: unnamed list of length %d\n", level + 1, length(x)))
    } else {
      name_str <- format_names(nms, n_names, indent(level))
      cat(indent(level), sprintf("Level %d: named list of length %d (from element [[1]]) with names: %s\n",
                                 level + 1, length(nms), name_str))
    }
    
    if (length(x) == 0) {
      cat(indent(level), "  (empty list)\n")
      return(invisible(NULL))
    }
    
    recurse(x[[1]], level + 1, append(path, list(x)))
  }
  
  recurse(x, level, list())
}
str_lol(molecular_data, n_names = 20)

#subset molecular data
tissue <- "blood"
ome <- "transcript-rna-seq"
md <- molecular_data[[tissue]][[ome]]$qc_norm

md$ENSG <- gsub("\\..*", "", rownames(md))

#load phenotypic data
pd <- fread("~/data/smontgom/motrpac_precovid_sedentary_clinical_obs.txt")

#merge the two
data("pheno")
ph <- as.data.table(pheno$data)
map <- ph[, .(labelID, pid, visit_code, Timepoint, study)]

## restrict omics to Adult‑Sedentary baseline blood samples
sel <- map[study == "Adult Sedentary" & visit_code == "ADU_SCP"]
keep <- intersect(sel$labelID, colnames(md))   # columns present in matrix

# -------------------------- #
# 1. subset & rename columns #
# -------------------------- #
md_sub <- md[ , keep ]                     # still a data.frame / matrix
colnames(md_sub) <- sel[match(keep, labelID), pid]

# md_sub now has participants as columns, row = gene
# ---------------------------------------- #
# 2. merge with your clin_wide data.table   #
# ---------------------------------------- #
expr_long <- as.data.table(t(md_sub), keep.rownames = "pid")
# t() → pid as first column, genes as variables (wide)

clin_full <- merge(clin_wide, expr_long, by = "pid", all.x = TRUE)

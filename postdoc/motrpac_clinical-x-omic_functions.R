mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}
build_sparse_mats <- function(x, row_names = NULL, verbose = F) {
  # x: list of length R; each element is a named list of numeric named vectors
  # returns: named list of Matrix::dgCMatrix, rows = length(x), cols = union of feature names per category
  
  if (!is.list(x) || length(x) == 0) stop("'x' must be a non-empty list.")
  R <- length(x)
  if (is.null(row_names)) {
    row_names <- if (!is.null(names(x)) && any(nzchar(names(x)))) names(x) else paste0("row", seq_len(R))
  }
  if (length(row_names) != R) stop("'row_names' must have same length as 'x'.")
  
  # categories are the names of each lower-level list (e.g., "phenotypes", ...)
  categories <- unique(unlist(lapply(x, function(li) {
    if (is.null(li)) character(0) else names(li)
  })))
  categories <- categories[nzchar(categories)]
  if (length(categories) == 0) stop("No lower-level categories (names) found.")
  
  if (verbose) message("found categories: ", paste(categories, collapse = ", "))
  
  # precompute per-category union of column names in order of first appearance
  colnames_by_cat <- setNames(vector("list", length(categories)), categories)
  for (cat in categories) {
    seen <- character(0)
    for (r in seq_len(R)) {
      li <- x[[r]]
      if (!is.null(li) && !is.null(li[[cat]])) {
        nm <- names(li[[cat]])
        nm <- nm[nzchar(nm)]
        if (length(nm)) {
          # append in first-appearance order
          add <- nm[!nm %in% seen]
          if (length(add)) seen <- c(seen, add)
        }
      }
    }
    colnames_by_cat[[cat]] <- seen
    if (verbose) message(sprintf("category '%s': %d unique columns", cat, length(seen)))
  }
  
  # build each sparse matrix via triplet form
  out <- setNames(vector("list", length(categories)), categories)
  for (cat in categories) {
    cols <- colnames_by_cat[[cat]]
    C <- length(cols)
    if (C == 0L) {
      # degenerate case: category present by name but no named entries; make a 0-column sparse matrix
      out[[cat]] <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                         dims = c(R, 0), dimnames = list(row_names, character(0)))
      next
    }
    
    # build index+value vectors
    ii <- integer(0)
    jj <- integer(0)
    xx <- numeric(0)
    
    col_index <- setNames(seq_len(C), cols)
    
    for (r in seq_len(R)) {
      li <- x[[r]]
      if (!is.null(li) && !is.null(li[[cat]])) {
        v <- li[[cat]]
        if (!is.numeric(v)) stop("all leaf entries must be numeric named vectors")
        nm <- names(v)
        if (is.null(nm) || any(!nzchar(nm))) stop("leaf vectors must have non-empty names")
        # keep only names that are in union (defensive)
        keep <- nm %in% cols
        if (any(keep)) {
          nm <- nm[keep]
          v  <- v[keep]
          ii <- c(ii, rep.int(r, length(v)))
          jj <- c(jj, unname(col_index[nm]))
          xx <- c(xx, unname(v))
        }
      }
    }
    
    if (verbose) {
      message(sprintf("assembling category '%s': %d nonzeros", cat, length(xx)))
    }
    
    mat <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                dims = c(R, C),
                                dimnames = list(row_names, cols))
    out[[cat]] <- mat
  }
  
  out
}
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
load_qc_local <- function(data_dir, tissue, ome, version = "v1.2") {
  # maps are straight from the directory/file naming you showed
  ome <- tolower(ome)
  tissue <- tolower(tissue)
  
  # 1) where to look + filename tokens
  folder_key <- switch(ome,
                       "transcript-rna-seq" = list(folder = "transcriptomics", key = "transcript-rna-seq", suffix = "log-cpm"),
                       "prot-ol"            = list(folder = "proteomics",     key = "prot-ol",            suffix = "log2"),
                       "metab"              = list(folder = "metabolomics-targeted", key = "metab-",      suffix = "log2"),
                       "metabolomics-targeted" = list(folder = "metabolomics-targeted", key = "metab-",   suffix = "log2"),
                       "metabolomics-untargeted" = list(folder = "metabolomics-untargeted", key = "metab-u-", suffix = "log2"),
                       stop("unsupported ome: ", ome)
  )
  
  # 2) tissue token depends on folder (their naming differs by assay)
  tissue_tok <- switch(folder_key$folder,
                       "transcriptomics" = switch(tissue,
                                                  "blood" = "t04-blood-rna",
                                                  "muscle" = "t06-muscle",
                                                  "adipose" = "t11-adipose",
                                                  stop("unsupported tissue for transcriptomics: ", tissue)
                       ),
                       "proteomics" = switch(tissue,
                                             "plasma" = "t02-plasma",
                                             "adipose" = "t07-adipose",
                                             "muscle"  = "t10-muscle",
                                             stop("unsupported tissue for proteomics: ", tissue)
                       ),
                       "metabolomics-targeted" = switch(tissue,
                                                        "plasma"  = "t02-plasma",
                                                        "adipose" = "t11-adipose",
                                                        "muscle"  = "t10-muscle", # if present
                                                        stop("unsupported tissue for metabolomics-targeted: ", tissue)
                       ),
                       "metabolomics-untargeted" = switch(tissue,
                                                          "plasma"  = "t02-plasma",
                                                          "adipose" = "t11-adipose",
                                                          stop("unsupported tissue for metabolomics-untargeted: ", tissue)
                       )
  )
  
  # 3) construct exact filenames (no search, no regex)
  qc_fp   <- file.path(
    data_dir, folder_key$folder, "qc-norm",
    paste0("human-precovid-sed-adu_", tissue_tok, "_",
           folder_key$key, "_qc-norm_", folder_key$suffix, "_", version, ".txt")
  )
  meta_fp <- file.path(
    data_dir, folder_key$folder, "metadata",
    paste0("human-precovid-sed-adu_", tissue_tok, "_",
           folder_key$key, "_metadata_samples_", version, ".txt")
  )
  
  if (!file.exists(qc_fp))   stop("qc file not found: ", qc_fp)
  if (!file.exists(meta_fp)) stop("metadata file not found: ", meta_fp)
  
  # 4) read (base R or data.table if you already use it)
  fread_exists <- requireNamespace("data.table", quietly = TRUE)
  read_tab <- function(fp) {
    if (fread_exists) data.table::fread(fp, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
    else read.table(fp, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  qc   <- read_tab(qc_fp)
  meta <- read_tab(meta_fp)
  
  # 5) make features rownames if first column is the feature id (common pattern)
  if (ncol(qc) > 1 && !anyDuplicated(qc[[1]]) && !(colnames(qc)[1] %in% colnames(meta))) {
    rn <- as.character(qc[[1]])
    qc <- qc[, -1, drop = FALSE]
    rownames(qc) <- rn
  }
  
  # 6) return in the same shape you used before
  out <- list()
  out[[tissue]] <- list()
  out[[tissue]][[ome]] <- list(
    qc_norm = qc,
    sample_metadata = meta
  )
  out
}

rep_by_group <- function(w, p) unlist(mapply(rep, w, p, SIMPLIFY = TRUE), use.names = FALSE)

zs <- function(M, mu, sd) {
  M <- t(t(M) - mu); M <- t(t(M) / sd); M
}

# fast univariate OLS loo predictions matrix and test preds
build_unilasso_stage2 <- function(X_tr, X_te, y_tr) {
  n <- nrow(X_tr); p <- ncol(X_tr)
  
  xty  <- as.numeric(t(X_tr) %*% y_tr)      # p
  xtx  <- colSums(X_tr^2)                   # p
  beta1 <- xty / pmax(xtx, .Machine$double.eps)
  sgn   <- sign(beta1); sgn[sgn == 0] <- 1
  
  Yhat  <- sweep(X_tr, 2, beta1, `*`)       # n x p
  Hdiag <- sweep(X_tr^2, 2, xtx, `/`)       # n x p
  
  ## FIX HERE: residuals must be (y - yhat), no extra sweep subtraction
  # OLD (buggy): E <- sweep(matrix(y_tr, n, p), 2, rep(1, p), `-`) - Yhat
  E <- matrix(y_tr, n, p) - Yhat            # n x p
  
  denom <- pmax(1 - Hdiag, 1e-8)
  Yhat_loo <- Yhat - E / denom              # PRESS identity
  
  Z_tr <- sweep(Yhat_loo, 2, sgn, `*`)
  Z_te <- sweep(X_te,     2, beta1, `*`)
  Z_te <- sweep(Z_te,     2, sgn,   `*`)
  
  list(Z_tr = Z_tr, Z_te = Z_te, sgn = sgn, beta1 = beta1)
}

residualize_feats <- function(output, input) {
  input <- as.matrix(cbind(Intercept = 1, input))
  XtX <- crossprod(as.matrix(input))
  B   <- chol2inv(chol(XtX)) %*% crossprod(input, output)  # q x p
  R   <- output - input %*% B                              # n x p
  list(R = R, B = B)
}

calibrate_linear <- function(y_train_z, yhat_train_z, yhat_test_z) {
  v <- var(yhat_train_z)
  if (v == 0) { b <- 1; a <- mean(y_train_z) - mean(yhat_train_z) } else {
    b <- cov(y_train_z, yhat_train_z) / v
    a <- mean(y_train_z) - b * mean(yhat_train_z)
  }
  list(a=a, b=b,
       yhat_train_cal_z = a + b * yhat_train_z,
       yhat_test_cal_z = a + b * yhat_test_z)
}

ifelse2 <- function(test, yes, no){if(test){return(yes)}else{return(no)}}

# compact duplicated column names by taking the single non-NA per row
compact_by_name <- function(df, verbose = TRUE) {
  # ensure column names exist
  if (is.null(colnames(df))) stop("data frame must have column names")
  
  # split column indices by name
  groups <- split(seq_along(df), colnames(df))
  n <- nrow(df)
  
  # helper to merge one group (vectorized; no apply over rows)
  merge_group <- function(idx, name) {
    # fast path: single column
    if (length(idx) == 1L) return(df[[idx]])
    
    # check type consistency to avoid character coercion surprises
    cls <- vapply(idx, function(j) class(df[[j]])[1L], character(1))
    if (length(unique(cls)) > 1L && verbose) {
      message("warning: mixed types in columns named '", name,
              "': ", paste(unique(cls), collapse = ", "),
              " -> coercing via as.matrix()")
    }
    
    m <- as.matrix(df[idx])  # becomes numeric/logical if consistent; character otherwise
    nz <- !is.na(m)
    
    # detect violation: >1 non-NA in a row
    bad <- rowSums(nz) > 1L
    if (any(bad) && verbose) {
      message("violation in '", name, "': ", sum(bad), " rows have multiple non-NA values")
    }
    
    # pick the first non-NA per row (arbitrary if >1; already warned)
    # max.col returns 1L for all-zero rows, so guard those
    has_any <- rowSums(nz) > 0L
    pick <- integer(n); pick[has_any] <- max.col(nz[has_any, , drop = FALSE], ties.method = "first")
    
    out <- rep(NA, n)
    if (any(has_any)) {
      out[has_any] <- m[cbind(which(has_any), pick[has_any])]
    }
    
    # attempt to cast back to a more specific type when possible
    if (is.numeric(type.convert(out, as.is = TRUE))) {
      suppressWarnings(out2 <- as.numeric(out))
      return(out2)
    } else if (all(out %in% c("TRUE","FALSE", NA))) {
      return(out == "TRUE")
    } else {
      return(out)
    }
  }
  
  # build the compacted data frame
  out_list <- mapply(merge_group, groups, names(groups), SIMPLIFY = FALSE, USE.NAMES = TRUE)
  out <- as.data.frame(out_list, check.names = FALSE, optional = TRUE)
  rownames(out) <- rownames(df)
  out
}

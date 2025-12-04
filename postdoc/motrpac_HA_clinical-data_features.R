library(data.table)
library(ape)

#load sed data first
sd_ph <- as.data.frame(fread("~/data/motrpac_precovid_sedentary_clinical_obs.txt"))

#print out column names for language model to parse
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
ha_paths <- list.files(data_dir)
for(i in seq_along(ha_paths)){
  ha_path <- paste0(data_dir, ha_paths[i])
  ha_ph <- as.data.frame(fread(ha_path))
  cat(paste0("\n\n", ha_path, ":\n"))
  cat(paste0("\"", colnames(ha_ph), "\","))
}

cat(paste0("\"", colnames(sd_ph), "\","))

data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_sets/"
ha_paths <- list.files(data_dir)
colnms <- list()
for(i in seq_along(ha_paths)){
  ha_path <- paste0(data_dir, ha_paths[i])
  ha_ph <- as.data.frame(fread(ha_path))
  colnms[[ha_path]] <- colnames(ha_ph)
}
table(table(unlist(colnms)))
head(sort(table(unlist(colnms)), T), 12)

#can try combining like strings?
all_strings <- unique(unlist(colnms))
str(all_strings)
all_strings[grepl("hrexer", all_strings)] #as an exanple

#### LLM approaches ####

# Build a string-similarity tree and plot to a multi-page PDF
# - Distance: Levenshtein via base::adist (or stringdist if you flip the switch)
# - Tree: hclust (average/complete/ward.D2) or NJ (unrooted), optional midpoint root
# - Readability: larger device, label cex/offset, ladderize, and cluster-based pagination
#
# Example quick run:
#   p <- string_similarity_tree_pdf(
#          all_strings,
#          outfile = "string_tree.pdf",
#          tree_method = "nj",
#          midpoint = TRUE,
#          tips_per_page = 250,
#          label_cex_full = 0.28,
#          label_cex_zoom = 0.85,
#          debug = TRUE
#        )

string_similarity_tree_pdf <- function(strings,
                                       outfile = "string_tree.pdf",
                                       # distance
                                       use_stringdist = FALSE,
                                       # tree
                                       tree_method = c("hclust", "nj"),
                                       hclust_method = "average",
                                       midpoint = FALSE,
                                       ladderize_tree = TRUE,
                                       # plotting controls
                                       label_cex_full = NULL,   # cex on the full tree page
                                       label_cex_zoom = NULL,   # cex on paginated zoom pages
                                       label_offset = NULL,     # tip label offset in branch-length units
                                       height_per_tip_full = 0.14, # inches per tip on full page
                                       height_per_tip_zoom = 0.20, # inches per tip on zoom pages
                                       pdf_width_in = 12,
                                       tips_per_page = 300,     # set to NULL/NA to disable pagination
                                       plot_types = c("phylogram", "fan", "unrooted"),
                                       debug = TRUE) {
  stopifnot(is.character(strings))
  strings <- strings[!is.na(strings)]
  n <- length(strings)
  if (n < 2L) stop("need at least 2 strings")
  
  # Ensure unique tip labels (ape requires uniqueness); preserve originals
  if (any(duplicated(strings))) {
    if (debug) message("duplicated labels found; making unique with '·' suffix")
    tip_labels <- make.unique(strings, sep = "·")
  } else {
    tip_labels <- strings
  }
  
  # ---------- distance ----------
  if (debug) message("distance: Levenshtein (", if (use_stringdist) "stringdist::stringdistmatrix" else "base::adist", ") on ", n, " strings")
  if (use_stringdist) {
    if (!requireNamespace("stringdist", quietly = TRUE)) stop("set use_stringdist=FALSE or install.packages('stringdist')")
    dm_full <- stringdist::stringdistmatrix(strings, strings, method = "lv")
  } else {
    dm_full <- utils::adist(strings, partial = FALSE, ignore.case = FALSE)
  }
  storage.mode(dm_full) <- "double"
  d <- stats::as.dist(dm_full)
  
  # ---------- tree ----------
  tree_method <- match.arg(tree_method)
  if (tree_method == "hclust") {
    if (debug) message("tree: hclust(method='", hclust_method, "') → as.phylo()")
    hc <- stats::hclust(d, method = hclust_method)
    hc$labels <- tip_labels
    ph <- ape::as.phylo(hc)
  } else {
    if (debug) message("tree: neighbor-joining (unrooted, non-ultrametric)")
    ph <- ape::nj(d)
    ph$tip.label <- tip_labels
  }
  
  # optional midpoint rooting
  if (isTRUE(midpoint)) {
    rooted <- FALSE
    if (requireNamespace("phangorn", quietly = TRUE)) {
      if (debug) message("midpoint rooting via phangorn::midpoint()")
      ph <- phangorn::midpoint(ph)
      rooted <- TRUE
    } else if (requireNamespace("phytools", quietly = TRUE)) {
      if (debug) message("midpoint rooting via phytools::midpoint.root()")
      ph <- phytools::midpoint.root(ph)
      rooted <- TRUE
    } else if (debug) {
      message("midpoint requested but 'phangorn'/'phytools' not installed; skipping")
    }
    if (ladderize_tree && rooted) ph <- ape::ladderize(ph, right = TRUE)
  } else if (ladderize_tree) {
    # ladderize even if not midpointed (helpful for hclust)
    ph <- ape::ladderize(ph, right = TRUE)
  }
  
  # ---------- sizing heuristics ----------
  # full-page text size: smaller for large n unless user overrides
  if (is.null(label_cex_full)) {
    label_cex_full <- max(0.22, min(1.0, 160 / n))
  }
  if (is.null(label_cex_zoom)) {
    label_cex_zoom <- max(0.40, min(1.0, 200 / min(n, tips_per_page %||% n)))
  }
  `%||%` <- function(a,b) if (!is.null(a) && !is.na(a)) a else b
  
  # compute device heights
  height_full <- max(8.5, min(100, n * height_per_tip_full))
  # zoom height chosen by cluster size later
  
  # label offset in branch-length units
  # default: small fraction of tree depth
  if (is.null(label_offset)) {
    # estimate depth as max root-to-tip distance
    depth <- tryCatch(max(ape::node.depth.edgelength(ph)[1:n]), error = function(e) 1)
    label_offset <- 0.03 * depth
  }
  
  # ---------- PDF ----------
  if (debug) message("writing PDF: ", outfile)
  grDevices::pdf(outfile, width = pdf_width_in, height = height_full, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  # helper to set xlim with space for labels
  plot_with_room <- function(tree, type = "phylogram", cex = 0.8, main = "", zoom = FALSE) {
    depth <- tryCatch(max(ape::node.depth.edgelength(tree)[1:length(tree$tip.label)]), error = function(e) 1)
    # leave extra space for labels; fan/unrooted ignore xlim
    extra <- if (type == "phylogram") max(0.15 * depth, label_offset * 6) else 0
    par(mar = c(1, 1, 3, 1))
    ape::plot.phylo(tree,
                    type = type,
                    cex = cex,
                    no.margin = TRUE,
                    label.offset = label_offset,
                    x.lim = if (type == "phylogram") c(0, depth + extra) else NULL,
                    main = main)
    if (type == "phylogram") ape::axisPhylo()
  }
  
  # ---------- Page 1: full tree views ----------
  if ("phylogram" %in% plot_types) {
    plot_with_room(ph, type = "phylogram", cex = label_cex_full,
                   main = sprintf("String similarity tree (%s; %d tips)", tree_method, n))
  }
  if ("fan" %in% plot_types) {
    plot_with_room(ph, type = "fan", cex = label_cex_full,
                   main = "Fan view")
  }
  if ("unrooted" %in% plot_types) {
    plot_with_room(ph, type = "unrooted", cex = label_cex_full,
                   main = "Unrooted view")
  }
  
  # ---------- Pagination: cluster zoom pages ----------
  if (!is.null(tips_per_page) && !is.na(tips_per_page) && tips_per_page < n) {
    if (debug) message("paginating into zoom pages; tips_per_page = ", tips_per_page)
    # get clusters by cutting an hclust (consistent with your distance)
    # if tree_method == nj, we still make an hclust to cut
    hc_for_cut <- stats::hclust(d, method = hclust_method)
    k <- ceiling(n / tips_per_page)
    grp <- stats::cutree(hc_for_cut, k = k)
    
    for (g in sort(unique(grp))) {
      tips <- tip_labels[grp == g]
      sub <- ape::keep.tip(ph, tips)
      if (ladderize_tree) sub <- ape::ladderize(sub, right = TRUE)
      h_zoom <- max(7.5, min(80, length(tips) * height_per_tip_zoom))
      grDevices::dev.control("enable")
      grDevices::dev.set(grDevices::dev.cur())
      grDevices::dev.new(width = pdf_width_in, height = h_zoom) # new page
      
      main_txt <- sprintf("Cluster %d/%d  (%d tips)", g, k, length(tips))
      if ("phylogram" %in% plot_types) {
        plot_with_room(sub, type = "phylogram", cex = label_cex_zoom, main = paste0(main_txt, " — phylogram"))
      }
      if ("fan" %in% plot_types) {
        plot_with_room(sub, type = "fan", cex = label_cex_zoom, main = paste0(main_txt, " — fan"))
      }
      if ("unrooted" %in% plot_types) {
        plot_with_room(sub, type = "unrooted", cex = label_cex_zoom, main = paste0(main_txt, " — unrooted"))
      }
    }
  } else if (debug) {
    message("pagination disabled (tips_per_page >= n or NULL)")
  }
  
  if (debug) {
    message("done. Saved to: ", normalizePath(outfile, mustWork = FALSE))
  }
  invisible(list(phy = ph, dist = d))
}

# build string-distance trees and plot to a PDF with readable labels
# default: ultrametric (UPGMA/hclust -> as.phylo)
# optional: unrooted NJ + midpoint root
string_tree_pdf2 <- function(strings,
                             outfile = "string_tree.pdf",
                             tree_type = c("ultrametric", "midpoint"),  # "midpoint" = NJ + midpoint root
                             hclust_method = "average",
                             width_in = 14,                 # page width in inches
                             height_per_tip_in = 0.045,     # inches per tip for page height scaling
                             min_height_in = 9,             # minimum page height
                             cex = NULL,                    # label size; if NULL, auto based on n
                             min_cex = 0.55,                # floor for auto cex
                             max_cex = 1.10,                # cap for auto cex
                             label_offset_in = 0.06,        # extra space between tip and tip label (inches)
                             ladderize_tree = TRUE,         # tidy left-right orientation
                             add_fan_page = TRUE,           # add a circular "fan" page
                             debug = TRUE) {
  if (!is.character(strings) || length(strings) < 2L) stop("strings must be character, length >= 2")
  keep <- !is.na(strings)
  if (any(!keep)) {
    if (debug) message("dropping ", sum(!keep), " NA entries")
    strings <- strings[keep]
  }
  
  # unique tip labels (phylo requires uniqueness)
  lbl <- make.unique(strings, sep = "·")
  
  # --- distance (Levenshtein) ---
  if (debug) message("computing Levenshtein distances (adist) for ", length(strings), " strings ...")
  dm <- utils::adist(strings, partial = FALSE, ignore.case = FALSE)
  storage.mode(dm) <- "double"
  rownames(dm) <- lbl
  colnames(dm) <- lbl
  d <- stats::as.dist(dm)
  attr(d, "Labels") <- lbl
  
  # --- pick tree type ---
  tree_type <- match.arg(tree_type)
  if (tree_type == "ultrametric") {
    if (debug) message("hclust(method = '", hclust_method, "') → as.phylo (ultrametric)")
    hc <- stats::hclust(d, method = hclust_method)
    hc$labels <- lbl
    ph <- ape::as.phylo(hc)
  } else {
    if (debug) message("building unrooted NJ tree ...")
    ph <- ape::nj(d)
    if (ladderize_tree) ph <- ape::ladderize(ph)
    # midpoint rooting (phangorn → phytools → message)
    if (requireNamespace("phangorn", quietly = TRUE)) {
      if (debug) message("midpoint rooting via phangorn::midpoint()")
      ph <- phangorn::midpoint(ph)
    } else if (requireNamespace("phytools", quietly = TRUE)) {
      if (debug) message("midpoint rooting via phytools::midpoint.root()")
      ph <- phytools::midpoint.root(ph)
    } else {
      warning("No midpoint-rooting package found (install 'phangorn' or 'phytools'); keeping NJ unrooted.")
    }
  }
  
  # optional ladderize for prettier layout
  if (ladderize_tree) ph <- ape::ladderize(ph)
  
  # --- sizing / margins ---
  n <- length(ph$tip.label)
  if (is.null(cex)) {
    # heuristic: scale cex gently with n; keep within [min_cex, max_cex]
    cex <- max(min_cex, min(max_cex, 0.9 * sqrt(200 / n)))
  }
  height_in <- max(min_height_in, n * height_per_tip_in)
  
  # compute right margin (inches) so long labels don’t get clipped
  # strwidth() needs a device; open a tiny off-screen device briefly
  grDevices::pdf(NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  max_chars <- max(nchar(ph$tip.label), na.rm = TRUE)
  char_w_in <- strwidth("W", units = "inches", cex = cex)
  # allow a bit of breathing room
  right_margin_in <- max(1, min(6, char_w_in * max_chars * 1.05))
  # label offset is additional space between the terminal node and the label
  label_offset <- label_offset_in
  
  # open the real PDF
  if (debug) message(
    sprintf("PDF size: %.1f × %.1f in  |  cex=%.2f  |  right margin=%.1f in  |  label.offset=%.2f in",
            width_in, height_in, cex, right_margin_in, label_offset)
  )
  grDevices::pdf(outfile, width = width_in, height = width_in, onefile = TRUE)
  
  # page 1: base R dendrogram, horizontal (only for ultrametric case using hclust heights)
  # if (tree_type == "ultrametric") {
  #   hc <- stats::hclust(d, method = hclust_method)
  #   hc$labels <- lbl
  #   op <- par(no.readonly = TRUE)
  #   # set margins in inches so labels have space on the right
  #   par(mai = c(0.6, 0.6, 0.6, right_margin_in))
  #   plot(stats::as.dendrogram(hc),
  #        horiz = TRUE,
  #        main = "String similarity tree (Levenshtein + UPGMA)",
  #        xlab = "edit distance",
  #        ylab = "",
  #        cex = cex)
  #   par(op)
  # }
  # 
  # # page 2: ape phylogram (good for reading labels)
  # op <- par(no.readonly = TRUE)
  # par(mai = c(0.4, 0.4, 0.8, right_margin_in))
  # ttl <- if (tree_type == "ultrametric") {
  #   "Phylogram (UPGMA; ultrametric)"
  # } else {
  #   "Phylogram (NJ; midpoint-rooted if available)"
  # }
  # ape::plot.phylo(ph,
  #                 type = "phylogram",
  #                 direction = "rightwards",
  #                 cex = cex,
  #                 label.offset = label_offset,
  #                 no.margin = TRUE,
  #                 main = ttl)
  # ape::axisPhylo()
  # par(op)
  
  # page 3 (optional): circular fan (often easier to eyeball large clusters)
  if (add_fan_page) {
    op <- par(no.readonly = TRUE)
    par(mai = c(0.4, 0.4, 0.8, 0.4))
    ape::plot.phylo(ph,
                    type = "fan",
                    cex = cex,
                    no.margin = TRUE,
                    main = "Fan view")
    par(op)
  }
  
  grDevices::dev.off()
  
  if (debug) {
    message("wrote: ", normalizePath(outfile, mustWork = FALSE))
  }
  invisible(list(dist = d, phylo = ph))
}

# ultrametric with bigger labels and generous spacing:
res_full <- string_tree_pdf2(all_strings,
                             outfile = "string_tree_readable.pdf",
                             tree_type = "ultrametric", 
                             width_in = 90,
                             height_per_tip_in = 0.25,  # try 0.05–0.06 if still tight
                             cex = 0.5,                   # bump to 0.8–1.0 if you want larger text
                             label_offset_in = 0.2,
                             debug = TRUE)

# midpoint-rooted NJ (non-ultrametric):
res_mid <- string_tree_pdf2(all_strings,
                            outfile = "string_tree_midpoint.pdf",
                            tree_type = "midpoint",
                            width_in = 80,
                            height_per_tip_in = 0.25,  # try 0.05–0.06 if still tight
                            cex = 0.5,                   # bump to 0.8–1.0 if you want larger text
                            label_offset_in = 0.2,
                            debug = TRUE)


#### process tree for text output ####

# longest common prefix / suffix (char-wise)
lcp <- function(x) {
  if (!length(x)) return("")
  p <- x[1]
  for (s in x[-1]) {
    m <- min(nchar(p), nchar(s)); i <- 0L
    while (i < m && substr(p, i+1L, i+1L) == substr(s, i+1L, i+1L)) i <- i + 1L
    p <- substr(p, 1L, i); if (p == "") break
  }
  p
}
lcsuf <- function(x) {
  if (!length(x)) return("")
  revs <- vapply(x, function(s) paste(rev(strsplit(s, "", fixed=TRUE)[[1]]), collapse=""), "")
  pr <- lcp(revs)
  paste(rev(strsplit(pr, "", fixed=TRUE)[[1]]), collapse="")
}

# numeric ranges like {0..5,7,9..11}
nums_to_brace <- function(nums) {
  nums <- sort(unique(as.integer(nums)))
  if (!length(nums)) return("{}")
  runs <- list(); start <- nums[1]; prev <- nums[1]
  for (k in nums[-1]) {
    if (k == prev + 1L) prev <- k else { runs[[length(runs)+1L]] <- c(start, prev); start <- k; prev <- k }
  }
  runs[[length(runs)+1L]] <- c(start, prev)
  parts <- vapply(runs, function(r) if (r[1]==r[2]) as.character(r[1]) else paste0(r[1],"..",r[2]), "")
  paste0("{", paste(parts, collapse=","), "}")
}

# top anchor substrings: frequent k-grams (k = 3..5) that occur in >= min_anchor_prop of strings
find_anchors <- function(strings, k_min=3L, k_max=5L, min_anchor_prop=0.6, max_anchors=3L) {
  n <- length(strings); if (!n) return(character(0))
  k_min <- max(1L, k_min); k_max <- max(k_min, k_max)
  # build doc-frequency map without exploding memory
  df <- new.env(parent = emptyenv())
  add_df <- function(key) { df[[key]] <- if (is.null(df[[key]])) 1L else df[[key]] + 1L }
  
  # per string, add unique substrings
  for (s in strings) {
    L <- nchar(s); if (L == 0) next
    seen <- new.env(parent = emptyenv())
    for (k in k_min:k_max) {
      if (k > L) break
      for (i in 1:(L - k + 1L)) {
        sub <- substr(s, i, i + k - 1L)
        if (is.null(seen[[sub]])) { seen[[sub]] <- TRUE; add_df(sub) }
      }
    }
  }
  # collect eligible anchors
  keys <- ls(df)
  if (!length(keys)) return(character(0))
  freq <- vapply(keys, function(k) df[[k]], integer(1))
  keep <- keys[freq >= ceiling(min_anchor_prop * n)]
  if (!length(keep)) return(character(0))
  
  # prefer longer & more frequent; remove anchors that are substrings of stronger ones
  ord <- order(nchar(keep), freq[match(keep, keys)], decreasing = TRUE)
  keep <- keep[ord]
  chosen <- character(0)
  for (cand in keep) {
    if (!any(vapply(chosen, function(z) grepl(z, cand, fixed=TRUE), logical(1)))) {
      chosen <- c(chosen, cand)
      if (length(chosen) >= max_anchors) break
    }
  }
  if (!length(chosen)) return(character(0))
  
  # order by median position across strings where they appear
  pos_median <- function(a) {
    pp <- suppressWarnings(vapply(strings, function(s) {
      p <- regexpr(a, s, fixed=TRUE); if (p[1] == -1L) NA_integer_ else as.integer(p[1])
    }, integer(1)))
    stats::median(pp[is.finite(pp)], na.rm = TRUE)
  }
  ord2 <- order(vapply(chosen, pos_median, numeric(1)))
  chosen[ord2]
}

# summarize a clade's strings into a compact pattern
summarize_clade <- function(strings,
                            max_enum = 10L,
                            min_anchor_prop = 0.6,
                            max_anchors = 3L) {
  if (!length(strings)) return(list(pattern = "", interesting = FALSE, kind = "empty"))
  
  pre <- lcp(strings)
  suf <- lcsuf(strings)
  core <- substring(strings, nchar(pre)+1L, nchar(strings) - nchar(suf))
  
  # strip obviously useless pre/suf (single punctuation chars)
  if (nchar(pre) <= 1L && grepl("^[^A-Za-z0-9]$", pre)) pre <- ""
  if (nchar(suf) <= 1L && grepl("^[^A-Za-z0-9]$", suf)) suf <- ""
  core <- substring(strings, nchar(pre)+1L, nchar(strings) - nchar(suf))
  
  # 1) pure numeric core → range
  if (all(grepl("^[0-9]+$", core))) {
    return(list(pattern = paste0(pre, nums_to_brace(as.integer(core)), suf),
                interesting = TRUE, kind = "numeric"))
  }
  # 2) small set of discrete cores → brace set
  ucore <- unique(core)
  if (length(ucore) <= max_enum && all(nchar(ucore) <= 50)) {
    return(list(pattern = paste0(pre, "{", paste(sort(ucore), collapse=","), "}", suf),
                interesting = TRUE, kind = "enum"))
  }
  
  # 3) anchors: pick frequent k-grams; stitch as pre*anchor1*anchor2*...*suf
  anchors <- find_anchors(strings, k_min=3L, k_max=5L,
                          min_anchor_prop = min_anchor_prop, max_anchors = max_anchors)
  if (length(anchors)) {
    body <- paste(anchors, collapse="*")
    pat <- paste0(if (nchar(pre)) paste0(pre, "*") else "",
                  body,
                  if (nchar(suf)) paste0("*", suf) else "")
    # collapse any accidental "**"
    pat <- gsub("\\*\\*+", "*", pat)
    # deem “interesting” if letters/digits cover at least 4 chars and '*' ratio < 0.6
    letters_digits <- gsub("[^A-Za-z0-9]+", "", pat)
    star_ratio <- nchar(gsub("[^*]", "", pat)) / max(1L, nchar(pat))
    ok <- nchar(letters_digits) >= 4L && star_ratio < 0.6
    return(list(pattern = pat, interesting = ok, kind = "anchors"))
  }
  
  # fallback: prefix/suffix only, if meaningful
  if (nchar(pre) + nchar(suf) >= 4L) {
    return(list(pattern = paste0(pre, "*", suf), interesting = TRUE, kind = "affix"))
  }
  
  list(pattern = "<diverse>", interesting = FALSE, kind = "none")
}

# descendants & depths
desc_tips <- function(phy, node) {
  N <- length(phy$tip.label)
  kids <- phy$edge[phy$edge[,1]==node, 2]
  out <- integer(0)
  for (k in kids) {
    if (k <= N) out <- c(out, k) else out <- c(out, desc_tips(phy, k))
  }
  out
}
node_depths <- function(phy) {
  N <- length(phy$tip.label); root <- N + 1L
  adj <- split(phy$edge[,2], phy$edge[,1])
  depth <- integer(max(phy$edge))
  q <- root; depth[root] <- 0L
  while (length(q)) {
    n <- q[1]; q <- q[-1]
    kids <- adj[[as.character(n)]]
    if (length(kids)) {
      for (c in kids) { depth[c] <- depth[n] + 1L; q <- c(q, c) }
    }
  }
  depth
}

# main: produce pruned, indented summaries
summarize_tree_text2 <- function(phy,
                                 strings_by_label,   # named by tip labels
                                 min_clade_size = 12L,
                                 max_enum = 12L,
                                 min_anchor_prop = 0.65,
                                 max_anchors = 3L,
                                 max_depth = Inf,
                                 show_counts = TRUE) {
  if (is.null(names(strings_by_label)) ||
      !all(phy$tip.label %in% names(strings_by_label))) {
    stop("strings_by_label must be a named character vector keyed by phy$tip.label")
  }
  N <- length(phy$tip.label); root <- N + 1L
  depth <- node_depths(phy)
  out <- character(0)
  
  emit <- function(node, parent_pat = NULL) {
    if (depth[node] > max_depth) return(invisible(NULL))
    tips <- desc_tips(phy, node)
    if (length(tips) < min_clade_size) return(invisible(NULL))
    
    labs <- phy$tip.label[tips]
    strs <- unname(strings_by_label[labs])
    summ <- summarize_clade(strs, max_enum=max_enum,
                            min_anchor_prop=min_anchor_prop,
                            max_anchors=max_anchors)
    
    # decide if worth printing
    pat <- summ$pattern
    same_as_parent <- !is.null(parent_pat) && identical(parent_pat, pat)
    dull <- !summ$interesting || pat %in% c("<diverse>", "", "*", "_", ".*")
    if (!same_as_parent && !dull) {
      indent <- paste0(rep("  ", depth[node]), collapse = "")
      hdr <- paste0(indent, "- ", pat, if (show_counts) paste0("  [n=", length(tips), "]") else "")
      out <<- c(out, hdr)
      parent_pat <- pat
    }
    # recurse regardless (children may be interesting even if parent wasn't)
    kids <- phy$edge[phy$edge[,1]==node, 2]
    for (k in kids) if (k > N) emit(k, parent_pat)
  }
  
  emit(root, parent_pat = NULL)
  out
}


phy <- res_full$phylo  # or wherever your phylo lives

# map tip labels to original strings (often identical)
strings_by_label <- structure(all_strings, names = all_strings)

lines2 <- summarize_tree_text2(
  phy,
  strings_by_label,
  min_clade_size = 15,     # up to you
  max_enum = 15,
  min_anchor_prop = 0.7,   # raise for tighter anchors; lower to be more permissive
  max_anchors = 3,         # 2–4 tends to read best
  max_depth = Inf,
  show_counts = TRUE
)

cat(paste(lines2, collapse = "\n"))
library(extrafont)
library(stringr)
library(psych)
library(corrplot)


#### functions ####

source("~/scripts/minor_scripts/postdoc/1Drepel_v2.R")
cell_text <- function(node, ns){
  str_squish(paste(xml_text(xml_find_all(node, ".//w:t", ns)), collapse = " "))
}

find_row_value <- function(doc_xml, pattern, ns) {
  # iterate over every table row
  rows <- xml_find_all(doc_xml, ".//w:tbl//w:tr", ns)
  for (row in rows) {
    cells  <- xml_find_all(row, ".//w:tc", ns)
    texts  <- map_chr(cells, cell_text, ns = ns)
    hit    <- which(str_detect(tolower(texts), pattern))
    if (length(hit)) {
      # take value from *next* cell if possible, else from the same cell
      idx <- if (hit[1] < length(cells)) hit[1] + 1 else hit[1]
      return(texts[idx])
    }
  }
  NA_character_
}

get_text <- function(node, ns) {
  # Collect every paragraph (<w:p>) inside this node …
  pars <- xml_find_all(node, ".//w:p", ns)
  
  # … collapse each paragraph’s runs (<w:t>) and explicit line-breaks (<w:br/>)
  para_vec <- map_chr(pars, function(p) {
    parts <- xml_find_all(p, ".//w:t | .//w:br", ns)
    paste(
      vapply(
        parts,
        function(x) if (xml_name(x) == "br") "\n" else xml_text(x),
        character(1)
      ),
      collapse = ""
    )
  })
  
  # join paragraphs with a single \n
  str_trim(paste(para_vec, collapse = "\n"))
}

# A robust function to map a score to the barplot's x-axis coordinates
map_score_to_x <- function(score) {
  bar_centers <- bp[-1] # Centers for scores 1, 2, 3, 4, 5
  # Linearly interpolate between the centers of the bars
  return(approx(x = 1:likert_score_max, y = bar_centers, xout = score)$y)
}

device_usr_coords <- function() {
  usr <- par("usr")     # current plot's user coords: c(x1, x2, y1, y2)
  plt <- par("plt")     # current plot region in normalized device coords: c(xleft, xright, ybottom, ytop)
  
  # compute scale: user units per NDC unit
  x_per_ndc <- (usr[2] - usr[1]) / (plt[2] - plt[1])
  y_per_ndc <- (usr[4] - usr[3]) / (plt[4] - plt[3])
  
  # extrapolate to full device range in user units
  x_min <- usr[1] - plt[1] * x_per_ndc
  x_max <- usr[2] + (1 - plt[2]) * x_per_ndc
  y_min <- usr[3] - plt[3] * y_per_ndc
  y_max <- usr[4] + (1 - plt[4]) * y_per_ndc
  
  # return named vector
  c(x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
}

# generic "label → value" hunter
grab_value <- function(doc_xml, label_re, ns) {
  #scan every table cell
  cells <- xml_find_all(doc_xml, ".//w:tc", ns)
  for (cell in cells) {
    text <- get_text(cell, ns)
    if (str_detect(str_to_lower(text), label_re)) {
      # same-cell form: "Label: value"
      after <- str_trim(str_remove(text,
                                   regex(label_re, ignore_case = TRUE))) |>
        str_remove("^[:\\-–]+\\s*")
      if (after != "") return(after)  # got value
      
      # next-cell form
      sib <- xml_find_first(cell, "following-sibling::w:tc", ns)
      if (!is.na(sib)) return(get_text(sib, ns))
    }
  }
  
  #fall back: label in a paragraph
  pars <- xml_find_all(doc_xml, ".//w:p", ns)
  for (p in pars) {
    txt <- get_text(p, ns)
    if (str_detect(str_to_lower(txt), label_re)) {
      after <- str_trim(str_remove(txt,
                                   regex(label_re, ignore_case = TRUE))) |>
        str_remove("^[:\\-–]+\\s*")
      if (after != "") return(after)
      
      nxt <- xml_find_first(p, "following-sibling::w:p", ns)
      if (!is.na(nxt)) return(get_text(nxt, ns))
    }
  }
  
  NA_character_
}

# txtwrap: wrap a string for base-graphics plotting using user-coordinate width
txtwrap <- function(string, cex = 1, width, break_long_words = TRUE, debug = FALSE) {
  
  if(length(string) > 1){
    return(sapply(string, function(si) txtwrap(si, cex = cex, width = width, 
                                               break_long_words = break_long_words,
                                               debug = debug)))
  }
  
  # check that a device with user coords is active
  ok <- TRUE
  tryCatch({ par("usr") }, error = function(e) { ok <<- FALSE })
  if (!ok) {
    stop("txtwrap requires an active plot with user coordinates. ",
         "Call plot.new(); plot.window(xlim=..., ylim=...) or draw a plot first.")
  }
  if (!is.character(string) || length(string) != 1L) {
    stop("string must be a single character string")
  }
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) {
    stop("width must be a single positive number in user units")
  }
  if (!is.numeric(cex) || length(cex) != 1L || !is.finite(cex) || cex <= 0) {
    stop("cex must be a single positive number")
  }
  
  # helper: measure a candidate line's width in user units
  line_width_user <- function(s) {
    # empty strings have ~0 width
    if (nchar(s) == 0L) return(0)
    strwidth(s, units = "user", cex = cex)
  }
  
  # wrap one paragraph (no internal '\n')
  wrap_paragraph <- function(p) {
    words <- strsplit(p, "\\s+", perl = TRUE)[[1L]]
    if (length(words) == 0L) return("")  # blank paragraph
    
    lines <- character(0L)
    current <- ""
    
    i <- 1L
    while (i <= length(words)) {
      w <- words[i]
      
      # if current line is empty, try to place the word (or break it)
      if (nchar(current) == 0L) {
        if (line_width_user(w) <= width) {
          current <- w
          if (debug) cat("[wrap] start line with:", shQuote(w),
                         "->", sprintf("%.4f <= %.4f\n", line_width_user(current), width))
          i <- i + 1L
        } else if (break_long_words) {
          # break the single word greedily to fit width
          # grow a chunk char-by-char until it would exceed width
          chars <- strsplit(w, "", fixed = TRUE)[[1L]]
          start <- 1L
          # ensure at least one character per line
          while (start <= length(chars)) {
            chunk <- chars[start]
            j <- start + 1L
            while (j <= length(chars)) {
              candidate <- paste0(paste0(chunk, collapse = ""), chars[j])
              if (line_width_user(candidate) <= width) {
                chunk <- c(chunk, chars[j])
                j <- j + 1L
              } else {
                break
              }
            }
            piece <- paste0(chunk, collapse = "")
            if (debug) cat("[wrap-long]", shQuote(piece),
                           "->", sprintf("%.4f <= %.4f\n", line_width_user(piece), width))
            lines <- c(lines, piece)
            start <- j
          }
          current <- ""
          i <- i + 1L
        } else {
          # place the long word on its own line, allow overflow
          if (debug) cat("[wrap] long word placed without breaking:", shQuote(w),
                         "->", sprintf("%.4f > %.4f\n", line_width_user(w), width))
          lines <- c(lines, w)
          current <- ""
          i <- i + 1L
        }
      } else {
        # try appending with a space
        candidate <- paste(current, w)
        if (line_width_user(candidate) <= width) {
          current <- candidate
          if (debug) cat("[wrap] append:", shQuote(w),
                         "->", sprintf("%.4f <= %.4f\n", line_width_user(current), width))
          i <- i + 1L
        } else {
          # commit current line and continue
          lines <- c(lines, current)
          current <- ""
          if (debug) cat("[wrap] line committed:", shQuote(lines[length(lines)]), "\n")
        }
      }
    }
    
    if (nchar(current) > 0L) {
      lines <- c(lines, current)
      if (debug) cat("[wrap] final line:", shQuote(current), "\n")
    }
    
    paste(lines, collapse = "\n")
  }
  
  # split by existing newlines and wrap each paragraph independently
  paragraphs <- strsplit(string, "\n", fixed = TRUE)[[1L]]
  wrapped <- vapply(paragraphs, wrap_paragraph, character(1L))
  paste(wrapped, collapse = "\n")
}

# Install and load necessary packages
library(extrafont)

#### read in data ####
randomize_data <- F
base_w_pix <- 500
nchar_pl <- 115
parmar <- c(12,6,6,6)
base_dir <- "~/Documents/Documents - nikolai/cam_review/"
results_dir <- paste0(base_dir, "results/")
input_dir <- paste0(base_dir, "input/")
# input_name <- c("2023 WAI Executive Director Evaluation (Responses) - Form Responses 1.csv")[1]
year <- 2024
ED_included <- F

if(ED_included){
  file_path_ED_self <- paste0(input_dir, year, " ED Self-Assessment.txt")
  d_all <- as.data.frame(data.table::fread(file = file_path_ED_self))
  ED_index <- nrow(d_all)
  d <- d_all
} else {
  input_name <- paste0(year, " WAI Executive Director Evaluation (Responses).txt")
  file_path <- paste0(input_dir, input_name)
  d <- as.data.frame(data.table::fread(file_path))
}

likert <- c("Strongly Disagree", "Disagree", "Neither Agree nor Disagree" ,"Agree" ,"Strongly Agree")
d[d == 0] <- NA
d[d == "N/A"] <- NA
d[d == ""] <- NA
for(ri in likert){
  d[d == ri] <- which(likert == ri)
}

# replace faulty encoding of dash and apostrophe :/
colnames(d) <- vapply(colnames(d), function(s) {
  r <- charToRaw(enc2native(s))
  out <- raw()
  for (b in r) {
    if (b == as.raw(0xD5)) {
      out <- c(out, charToRaw("’"))  # or charToRaw("'") for ASCII apostrophe
    } else if (b == as.raw(0xD0)) {
      out <- c(out, charToRaw("–"))
    } else {
      out <- c(out, b)
    }
  }
  rawToChar(out)
}, character(1), USE.NAMES = FALSE)
colnames(d) <- as.character(sapply(colnames(d), function(x){
  if(substr(x, nchar(x), nchar(x)) == "’"){
    return(substr(x, 1, nchar(x)-1))
  } else {
    return(x)
  }
}))
colnames(d) <- as.character(sapply(colnames(d), function(x){
  if(substr(x, 1, 1) == "\n"){
    return(substr(x, 2, nchar(x)))
  } else {
    return(x)
  }
}))

#did this manually before
# col.names <- c("Timestamp",
#                "The ED clearly communicates relevant expectations and goals to me and my team.",
#                "The ED responds to requests or communications in a timely manner.",
#                "The ED understands and helps me to understand how my work contributes to organizational strategy.",
#                "The Executive Director clearly communicates motivations and justifications for their decisions. I understand why they make the choices they do, even when I don’t agree with them.",
#                "I have a clear understanding of how the ED perceives my work performance at WAI.",
#                "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
#                "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
#                "In 2023, the ED made excellent progress toward the ORCAs for which they are primarily responsible.",
#                "The ED does not trust my expertise in my specific work area, and tends to micromanage projects that would better be left to me.",
#                "I am inspired by the ED’s efforts to boost team morale and spirit. Their enthusiasm seems sincere and resonates with me.",
#                "I often disagree with decisions made by the ED in my work area.",
#                "I would feel comfortable approaching the ED directly with suggestions on how they could improve in their role at WAI.",
#                "I would feel comfortable approaching the ED directly about something impeding my success and happiness at WAI that they are not directly responsible for.",
#                "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
#                "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
#                "The ED is sufficiently familiar with the technical aspects of wild animal welfare science to be able to provide meaningful input and direction for research at WAI.",
#                "The ED is effective at time management – they appropriately prioritize important tasks and are able judge less important ones to ignore or delegate.",
#                "I am satisfied with the amount of time I interact with the ED.",
#                "The ED provides regular feedback in ways that effectively facilitate my own improvement, both in the form of constructive criticism and in positive affirmation of jobs well done.",
#                "The ED has been effective at obtaining funding for WAI’s continued operation and growth.",
#                "I am proud of how the ED represents WAI when engaging with external stakeholders.",
#                "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
#                "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
#                "The ED fosters an inclusive work environment, both in supporting staff and in intervening against intolerant work dynamics.",
#                "I have been made uncomfortable by the ED in the last 12 months.",
#                "The ED is proactive in identifying ways in which the WAI work environment can be improved.",
#                "Disagreements with the ED about approach, mission, and tactics are handled appropriately and respectfully.",
#                "The ED successfully moderates disagreements among staff.",
#                "The ED treats all employees equitably and without favoritism.",
#                "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
#                "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
#                "In what important areas has the ED been most successful, relative to your expectations for their role? In what ways would they likely outperform a different ED?",
#                "In what important areas has the ED been least successful, relative to your expectations for their role? In what ways would a different ED likely outperform them?",
#                "Imagine you have the opportunity to direct how the ED spends 40 hours of their focused time over the next year, without affecting their other responsibilities. How would you allocate these hours to maximize the ED's effectiveness in their role at WAI? Responses can include, but are not limited to, professional development courses, wellness activities, skill-building workshops, or any other accessible experiences that you believe would enhance their performance and contribution to WAI. Please be as specific as possible in your suggestions.",
#                "As above, but corresponding to only 1 hour of time.",
#                "Please put any comments on the survey as a whole here! Responses here will not be shared with Cam directly and will only be visible to the board. You may also use this space to privately respond to the above questions.",
#                "Score")
# colnames(d) <- col.names
inverted_qs <- colnames(d)[c(10, 12, 26)]

d <- d[,-1]
d <- d[-1,]
splits <- c(-1, 
            which(grepl(colnames(d), pattern = "Please use this space for comments on the above statements or your selections that will be shared with Cam")),
            35)
ds <- lapply(1:(length(splits)-1), function(i){
  d[,(splits[i] + 2):(splits[i+1])]
})
names(ds) <- c("Communication", 
               "Leadership", 
               "Operational_Effectiveness", 
               "Workplace_Culture_and_Environment", 
               "Free_Response")
cat_key <- do.call(rbind, sapply(1:length(ds), function(i) cbind(colnames(ds[[i]]), names(ds)[i])))
cat_key <- setNames(cat_key[,2], cat_key[,1])
cat_cols <- setNames(RColorBrewer::brewer.pal(5, "Dark2"), names(ds))

#### some quick analyses of these data ####

#construct table
dn <- apply(d, 2, as.numeric)
dn <- dn[,!apply(apply(dn, 2, is.na), 2, all)]
dn[,inverted_qs] <- 6 - dn[,inverted_qs]

if(randomize_data){
  for(i in 1:ncol(dn)){
    dn[,i] <- sample(x = c(NA, 1:5), size = nrow(dn), replace = T, prob = c(exp(rnorm(6))))
  }
}

#remove Cam from data
if(ED_included){
  d_ED <- d[ED_index-1,]
  d <- d[-(ED_index-1),]
  dn_ED <- dn[ED_index-1,, drop = F]
  dn <- dn[-(ED_index-1),]
  ED_scores <- apply(dn_ED, 2, mean, na.rm = TRUE)
}

#compute metrics
item_means <- apply(dn, 2, mean, na.rm = TRUE)
item_beta_means <- (apply(dn, 2, sum, na.rm = TRUE) + 1) / 
  (apply(!apply(dn, 2, is.na), 2, sum) + 1)
item_beta_shapes <- cbind(1 + apply(dn, 2, sum, na.rm = TRUE),
                     1 + (apply(dn, 2, function(x) sum(!is.na(x)))) * 5)

item_90CI_beta <- do.call(rbind, lapply(setNames(1:nrow(item_beta_shapes), nm = rownames(item_beta_shapes)), function(i){
  qbeta(p = c(0.05, 0.95), 
        shape1 = item_beta_shapes[i,1], 
        shape2 = item_beta_shapes[i,2] - item_beta_shapes[i,1])
})) * 5

item_50CI_beta <- do.call(rbind, lapply(setNames(1:nrow(item_beta_shapes), nm = rownames(item_beta_shapes)), function(i){
  qbeta(p = c(0.25, 0.75), 
        shape1 = item_beta_shapes[i,1], 
        shape2 = item_beta_shapes[i,2] - item_beta_shapes[i,1])
})) * 5

item_beta_median <- do.call(rbind, lapply(setNames(1:nrow(item_beta_shapes), nm = rownames(item_beta_shapes)), function(i){
  qbeta(p = c(0.5), 
        shape1 = item_beta_shapes[i,1], 
        shape2 = item_beta_shapes[i,2] - item_beta_shapes[i,1])
})) * 5

item_90CI_betabinom <- do.call(rbind, lapply(setNames(1:nrow(item_beta_shapes), nm = rownames(item_beta_shapes)), function(i){
  TailRank::qbb(p = c(0.05, 0.95), N = 5,
        u = item_beta_shapes[i,1], 
        v = item_beta_shapes[i,2] - item_beta_shapes[i,1])
}))

item_beta_sds <- apply(item_beta_shapes, 1, function(x){
  sqrt(x[1] * x[2] / (x[1]+x[2])^2 / (x[1] + x[2] + 1)) * 5
})


item_sds <- apply(dn, 2, sd, na.rm = TRUE)
cronbachs_alpha <- alpha(dn, check.keys = T)
item_mean_correlations <- setNames(sapply(1:ncol(dn), function(i){
  cor(dn[,i], apply(dn[,-i], 1, mean, na.rm = T), use = "complete.obs")
}), colnames(dn))
item_corrmat <- cor(dn, use = "pair")
item_order <- hclust(as.dist(1-item_corrmat))$order
item_corrmat <- item_corrmat[item_order,item_order]
# rownames(item_corrmat) <- colnames(item_corrmat) <- match(colnames(item_corrmat), colnames(dn))
# corrplot(as.matrix(item_corrmat), method = "circle")

person_corrmat <- cor(t(dn), use = "pair") #represents agreement between parties
person_order <- hclust(as.dist(1-person_corrmat))$order
person_corrmat <- person_corrmat[person_order,person_order]
corrplot(as.matrix(person_corrmat), method = "square")

#### prep fig metadata ####
col_dat <- data.frame(comms_h = numeric(4), 
                       base_h = numeric(4),
                       comms = character(4))

for(si in 1:4){
    
  subds <- ds[[si]][,1:(length(ds[[si]])-1)]
  
  base_h <- 275 * ncol(subds) / 72
  cairo_pdf(tempfile(),
            width = base_w_pix / 72, height = base_h, family = "Crimson Text")
  
  comms <- ds[[si]][,length(ds[[si]])]
  
  #mark Cam's response here
  if(ED_included){
    if(!is.na(comms[ED_index-1])){
      comms[ED_index-1] <- paste0("Cam's Reflection:\n\n", comms[ED_index-1])
    }
  }
  #fix encoding issue
  comms <- vapply(comms, function(s) {
    r <- charToRaw(enc2native(s))
    out <- raw()
    for (b in r) {
      if (b == as.raw(0xD5)) {
        out <- c(out, charToRaw("’"))  # or charToRaw("'")
      } else if (b == as.raw(0xD0)) {
        out <- c(out, charToRaw("-"))
      } else {
        out <- c(out, b)
      }
    }
    rawToChar(out)
  }, character(1))
  
  if(randomize_data){
    for(i in 1:length(comms)){
      if(!is.na(comms[i])){
        comms[i] <- paste0(sample(strsplit(comms[i], "")[[1]]), collapse = "")
      }
      
    }
  }
  comms <- comms[!is.na(comms) & comms != "NA"]
  comms <- paste0(comms, collapse = "~~\n")
  
  #redact identifying info
  comms <- gsub(" Cat.", " REDACTED.", comms)
  
  comms_parts <- str_split(comms, "\n")[[1]]
  comms_wrapped <- sapply(comms_parts, str_wrap, width = nchar_pl)
  comms <- paste(comms_wrapped, collapse = "\n")
  
  comms <- gsub("~~", paste0("\n\n", paste0(rep("-", nchar_pl + 10), collapse = ""),"\n\n"), comms)
  comms <- gsub("\n\n\n", "\n\n", comms)
  comms_h <- strheight(comms, units = "inches")
  
  dev.off()
  
  col_dat$comms_h[si] <- comms_h
  col_dat$comms[si] <- comms
  col_dat$base_h[si] <- base_h

}

#write sumstats to disk
num_responses_per_category <- do.call(rbind, lapply(apply(dn, 2, factor, levels = 1:5, simplify = F), table))
colnames(num_responses_per_category) <- paste0("# respondents answering \"", likert, "\"")
reversed_inds <- which(rownames(num_responses_per_category) %in% inverted_qs)
for(ri in reversed_inds){
  num_responses_per_category[ri,] <- rev(num_responses_per_category[ri,])
}

item_data <- data.frame(
  label = names(item_beta_means),
  beta_mean = item_beta_means,
  beta_sds = item_beta_sds, 
  beta_90CI_lb = item_90CI_beta[,1],
  beta_90CI_ub = item_90CI_beta[,2],
  beta_50CI_lb = item_50CI_beta[,1],
  beta_50CI_ub = item_50CI_beta[,2],
  beta_shape_1 = item_beta_shapes[,1],
  beta_shape_2 = item_beta_shapes[,2],
  beta_median = item_beta_median,
  "sorted_scores (1=worst, 5=best)" = sapply(apply(dn, 2, sort), paste0, collapse = "; "),
  num_responses_per_category)
data.table::fwrite(item_data, file = paste0(results_dir, year, "_item_summary_data.csv"))

#### one big figure ####
if(randomize_data){
  cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/example_output.pdf"),
            width = base_w_pix / 72 * 6, height = max(col_dat$base_h + col_dat$comms_h), family = "Crimson Text")
} else {
  cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", year, "_ED-Eval_Results_WAI.pdf"),
            width = base_w_pix / 72 * 6, height = max(col_dat$base_h + col_dat$comms_h), family = "Crimson Text")
}
par(xpd = NA, c(12,5,5,2))

txt_cex <- 1.3
row_heights <- c(0.175, rep(1, 6))
row_heights <- c(row_heights, sum(row_heights) * max((col_dat$comms_h / col_dat$base_h)[sapply(ds, ncol) - 1 == 6]))
nmc <- sapply(ds, ncol) - 1
nrows <- nmc + 2
layout_mat <- do.call(cbind, lapply(nrows, function(x) c(1:x, rep(x, 8 - x))))[,-5]
layout_mat <- layout_mat + matrix(rep(cumsum(c(0, apply(layout_mat, 2, max)[
  -ncol(layout_mat)])), each = nrow(layout_mat)), ncol = 4)
layout_mat <- cbind(layout_mat, c(max(layout_mat) + 1,
                                  rep(max(layout_mat)+2, 3),
                                  rep(max(layout_mat)+3, 3),
                                  max(layout_mat)+4
                                  ))
col_widths <- c(rep(1, ncol(layout_mat)-1), 2)
layout(layout_mat, heights = row_heights, widths = col_widths)

for(si in 1:4){
  
subds <- ds[[si]][,1:(length(ds[[si]])-1)]
subds_reorder <- order(item_means[colnames(subds)])
subds <- subds[,subds_reorder]
  
col_dat$comms_h[si] -> comms_h
col_dat$comms[si] -> comms
col_dat$base_h[si] -> base_h

# cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", names(ds)[si], ".pdf"),
#           width = 500 / 72, height = base_h + comms_h, family = "Crimson Text")

par(xpd = NA, mar = c(0,0,0,0))

plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
sec_name <- names(ds)[si]
sec_name_cex <- c(5, 5, 4, 4)
text(gsub("_", " ", sec_name), y = 0.4, x = 0, cex = sec_name_cex[si], col = cat_cols[sec_name])

par(mar = parmar)
likert_names <- c("N/A", "Strongly Disagree", "Disagree", "Neither Agree nor Disagree", "Agree", "Strongly Agree")
for(i in 1:ncol(subds)){

  qtext <- colnames(subds)[i]
  vals <- as.numeric(subds[,i])
  n_na <- sum(is.na(vals))
  if(qtext %in% inverted_qs){
    vals <- 6 - vals
    likert_names <- c("N/A", "Strongly Disagree", "Disagree", "Neither Agree nor Disagree", "Agree", "Strongly Agree")
    likert_names[2:6] <- rev(likert_names[2:6])
  } else {
    likert_names <- c("N/A", "Strongly Disagree", "Disagree", "Neither Agree nor Disagree", "Agree", "Strongly Agree")
  }
  counts <- table(factor(vals, levels = 1:5))
  counts <- c(n_na, counts)
  
  
  # plot it
  bp <- barplot(counts, names.arg = likert_names, xlab = "", ylab = "Count", 
                yaxt = "n", xaxt = "n", col = c("grey90", rep("#00c1b2", 5)), cex.lab = 1.75)
  qlab <- str_wrap(paste0((1:length(subds_reorder))[subds_reorder][i], ". ", qtext), width = 90)
  text(label = qlab, 
       x = mean(par("usr")[1:2]), 
       y = par("usr")[4] + diff(par("usr")[3:4])/4 - strheight(qlab),
       pos = 3, cex = 1.5, col = "darkred")
  axis(2, at = seq(0, max(counts), by = 1), 
       las = 2, labels = seq(0, max(counts), by = 1), cex.axis = 1.5)
  text(bp, par("usr")[3] - diff(par("usr")[3:4])/5, srt = 30, adj = 1, 
       labels = likert_names, xpd = TRUE, xpd = NA, cex = 1.5)
  segments(x0 = bp, x1 = bp, y0 = par("usr")[3], y1 = par("usr")[3] - diff(par("usr")[3:4])/6,
           lty = 3, lwd = 2, col = "grey")
  
  #label mean, 
  mean_val <- item_beta_means[qtext]
  unit_intervals <- c(bp) + diff(bp[1:2]) / 2
  mean_loc <- unit_intervals[1] + mean_val * diff(unit_intervals[1:2])
  item_90CI_beta_locs <- unit_intervals[1] + item_90CI_beta[qtext,] * diff(unit_intervals[1:2])
  item_90CI_betabinom_locs <- unit_intervals[1] + item_90CI_betabinom[qtext,] * diff(unit_intervals[1:2])
  
  ybox_locs <- c(par("usr")[3] - diff(par("usr")[3:4])/10, 
                 par("usr")[3] - diff(par("usr")[3:4])/30)
  segments(x0 = item_90CI_betabinom_locs[1], x1 = item_90CI_betabinom_locs[2], 
           y0 = mean(ybox_locs), y1 = mean(ybox_locs),
           lwd = 2)
  rect(xleft = item_90CI_beta_locs[1], xright = item_90CI_beta_locs[2], 
       ybottom = ybox_locs[1], ytop = ybox_locs[2],
       lwd = 2, col = "white")
  segments(x0 = mean_loc, x1 = mean_loc, 
           y0 = ybox_locs[1], y1 = ybox_locs[2],
       lwd = 2)
  # points(x=mean_loc, y=par("usr")[3] - diff(par("usr")[3:4])/10, pch = 23, cex = 3, col = 1, bg = "#c1000f")
  
  #plot Cam's self-assessment
  if(ED_included){
    ED_score <- ED_scores[qtext]
    ED_loc <- mean(unit_intervals[ED_score + c(0,1)])
    points(x = ED_loc, y = counts[ED_score+1], pch = 18, col = 1, cex = 5.5, xpd = NA)
    points(x = ED_loc, y = counts[ED_score+1], pch = 18, col = "#b306c2", cex = 5, xpd = NA)
  }
  
  #annotation text
  cronbachs_alpha_diff <- cronbachs_alpha$total$raw_alpha - cronbachs_alpha$alpha.drop[qtext,"raw_alpha"]
  
  annot <- latex2exp::TeX(paste0("$\\hat{\\mu}$ = ", round(item_means[qtext], 2), 
                                 " (#", rank(-item_means)[qtext], "/", ncol(dn), ")"))
  annot_y <- max(counts) - strheight(annot, cex = 1.25) / 2
  text(x = par("usr")[1] + ifelse(which.max(counts) == 1, bp[2] - diff(bp[1:2]) / 2, 0), 
       y = annot_y, labels = annot, pos = 4, cex = 1.25)
  # points(x= strwidth(annot, cex = 1.25) - 0.1275 + ifelse(which.max(counts) == 1, bp[2] - diff(bp[1:2]) / 2, 0),
  #        y=annot_y + diff(par("usr")[3:4])/60, pch = 23, cex = 1.5, col = 1, bg = "#c1000f")
  
  annot_y <- annot_y - strheight(annot, cex = 1.25) * 1.1
  annot <- latex2exp::TeX(paste0("$\\hat{\\sigma}$ = ", round(item_sds[qtext], 2), 
                                 " (#", rank(-item_sds)[qtext], "/", ncol(dn), ")"))
  text(x = par("usr")[1] + ifelse(which.max(counts) == 1, bp[2] - diff(bp[1:2]) / 2, 0), 
       y = annot_y, labels = annot, pos = 4, cex = 1.25)
  
  annot_y <- annot_y - strheight(annot, cex = 1.25) * 1.1
  annot <- latex2exp::TeX(paste0("$\\Delta$$\\alpha$ = ", #for diff in cronbach's alpha
                                 round(cronbachs_alpha$total$raw_alpha, 2), " - ", 
                                 round(cronbachs_alpha$total$raw_alpha - cronbachs_alpha_diff, 2), 
                                 "$\\approx$ ", round(cronbachs_alpha_diff, 2), 
                                 " (#", rank(cronbachs_alpha$alpha.drop[,"raw_alpha"])[which(qtext == rownames(cronbachs_alpha$alpha.drop))], "/", ncol(dn), ")"))
  text(x = par("usr")[1] + ifelse(which.max(counts) == 1, bp[2] - diff(bp[1:2]) / 2, 0), 
       y = annot_y, labels = annot, pos = 4, cex = 1.25)
  
  annot_y <- annot_y - strheight(annot, cex = 1.25) * 1.1
  annot <- latex2exp::TeX(paste0("item-mean \\textit{r} = ", round(item_mean_correlations[qtext], 2), 
                                 " (#", rank(-item_mean_correlations)[qtext], "/", ncol(dn), ")"))
  text(x = par("usr")[1] + ifelse(which.max(counts) == 1, bp[2] - diff(bp[1:2]) / 2, 0), 
       y = annot_y, labels = annot, pos = 4, cex = 1.25)
  
}


adj_y <- c(0,0,0,0)
plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
text("Please use this space for comments on the above statements or your selections that will be shared with Cam:", 
     y = 1.075 + adj_y[si], x = -0.1, cex = txt_cex, font = 2, xpd = NA, col = 1)
text(comms, y = 1 + adj_y[si] - 
       strheight(comms, units = "user", family = "Crimson Text") * txt_cex * 1.2 / 2 -
       strheight(" ", units = "user", family = "Crimson Text") * txt_cex * 3, 
     x = -1.3, cex = txt_cex, pos = c(4), xpd = NA, col = "grey20")


}

par(xpd = NA, mar = c(0,0,0,0))

plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
# text(gsub("_", " ", names(ds)[5]), y = 0.4, x = 0, cex = 5, col = "#00c1b2")

# plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
# yloc <- par("usr")[4]
# for(i in 1:ncol(ds[[5]])){
# # for(i in 1){
#   for(j in 1:2){
#     
#     if(j == 1){
#       comms <- colnames(ds[[5]])[i]  
#     } else {
#       comms <- ds[[5]][,i]
#     }
#     
#     comms <- comms[!is.na(comms)]
#     comms <- paste0(comms, collapse = "~~\n")
#     
#     #redact identifying info
#     comms <- gsub(" Cat.", " REDACTED.", comms)
#     
#     nchar_pl <- 240
#     comms_parts <- str_split(comms, "\n")[[1]]
#     comms_wrapped <- sapply(comms_parts, str_wrap, width = nchar_pl)
#     comms <- paste(comms_wrapped, collapse = "\n")
#     
#     comms <- gsub("~~", paste0("\n\n", paste0(rep("-", nchar_pl + 10), collapse = ""),"\n\n"), comms)
#     comms <- gsub("\n\n\n", "\n\n", comms)
#     comms_h <- strheight(comms, units = "inches")
#     
#     
#     yadj <- c(1,1,1,1) / 10 / 1.15
#     text(comms, pos = 4,
#          y = yloc - 
#            strheight(comms, units = "user", family = "Crimson Text", cex = txt_cex) / 2 - 
#            ifelse(j == 2, yadj[i], 0) -
#            strheight(" ", units = "user", family = "Crimson Text", cex = txt_cex) * 3, 
#          x = -1.1, cex = txt_cex, font = ifelse(j==1, 2, 1), xpd = NA)
#     
#     
#     yloc <- yloc - 
#       strheight(comms, units = "user", family = "Crimson Text", cex = txt_cex) -
#       strheight(" ", units = "user", family = "Crimson Text", cex = txt_cex)  * 3 - 
#       ifelse(j == 2, yadj[i], 0)
#   }
# }

# corrplot(as.matrix(item_corrmat), method = "circle")
# Create a custom color palette
custom_palette <- colorRampPalette(c("#B200C1", "white", "#00C1B2"))

# Generate the corrplot with the custom color palette
corrplot(as.matrix(item_corrmat), method = "square", 
         col = custom_palette(200),  # Use 200 colors for a smooth gradient
         tl.cex = 1.6,              
         tl.col = "black",          
         tl.pos = "n",              
         addgrid.col = NA,          
         cl.cex = 3,                
         mar = c(0, 0, 15, 15), cl.length = 5, cl.align.text = "l", diag = F)
text(latex2exp::TeX("Pearson's \\textit{r} Between Individual Items"), 
     y = mean(par("usr")[3:4]), srt = 270, cex = 5, x = par("usr")[2] - diff(par("usr")[1:2])/25, 
     col = "#00c1b2")


truncate_label <- function(label, max_length = 65) {
  if (nchar(label) > max_length) {
    words <- unlist(strsplit(label, " "))
    truncated_label <- ""
    for (word in words) {
      if (nchar(paste0(truncated_label, word, sep = " ")) > max_length) break
      truncated_label <- paste0(truncated_label, word, " ")
    }
    return(paste0(trimws(truncated_label), "..."))
  } else {
    return(label)
  }
}

text(1:ncol(item_corrmat), ncol(item_corrmat) + 1, srt = 30, pos = 4,
     labels = sapply(colnames(item_corrmat), truncate_label), xpd = TRUE, xpd = NA, cex = 1.5,
     col = cat_cols[cat_key[colnames(item_corrmat)]])
legend(par("usr")[2] - diff(par("usr")[1:2])/15, par("usr")[4] + diff(par("usr")[3:4])/15, box.lwd = 0, col = cat_cols[-5], 
       legend = gsub("_|and_Environment", " ", names(cat_cols))[-5], pch = 15, pt.cex = 3, cex = 2)
rect(xleft = 1:ncol(item_corrmat) - 1/2, 
     ybottom = ncol(item_corrmat):1 - 1/2, 
     xright = 1:ncol(item_corrmat) + 1/2, 
     ytop = ncol(item_corrmat):1 + 1/2,
     col = cat_cols[cat_key[colnames(item_corrmat)]], border = "white")
text(1:ncol(item_corrmat), ncol(item_corrmat):1, 
     labels = sapply(colnames(item_corrmat), function(x) match(x, colnames(ds[[cat_key[x]]]))), 
     col = "white", cex = 3)

#histograms for average person-rating, all ratings, item-mean correlations


# Generate the corrplot with the custom color palette
corrplot(as.matrix(person_corrmat), method = "square", 
         col = custom_palette(200),  # Use 200 colors for a smooth gradient
         tl.cex = 1.6,              
         tl.col = "black",          
         tl.pos = "n",              
         addgrid.col = NA,          
         cl.cex = 3,                
         mar = c(10, 0, 5, 15), cl.length = 5, cl.align.text = "l", diag = F)
rect(xleft = 1:ncol(person_corrmat) - 1/2, 
     ybottom = ncol(person_corrmat):1 - 1/2, 
     xright = 1:ncol(person_corrmat) + 1/2, 
     ytop = ncol(person_corrmat):1 + 1/2,
     col = "grey", border = "white")
text(x = 1:ncol(person_corrmat), 
     y = ncol(person_corrmat):1, 
     labels = 1:ncol(person_corrmat),
     col = "white", cex = 4)
text(latex2exp::TeX("Pearson's \\textit{r} Between Individual Respondents"), 
     y = mean(par("usr")[3:4]), srt = 270, cex = 5, x = par("usr")[2] - diff(par("usr")[1:2])/25, 
     col = "#00c1b2")

#add in a Legend
yloc <- -1.5
txt_cex <- 5
text(latex2exp::TeX("Legend"), pos = 4,
     y = yloc, srt = 0, cex = txt_cex, x = 1, 
     col = "#00c1b2", font = 2)
yloc <- yloc - strheight(latex2exp::TeX("Legend"), cex = txt_cex) * 1.1

t2w <- list(latex2exp::TeX(paste0("$\\hat{\\mu}$ = Sample Mean")),
            latex2exp::TeX(paste0("item-mean \\textit{r} is Pearson's \\textit{r} between item score and overall score")),
            latex2exp::TeX(paste0("#s in (parentheses) indicate ranking")),
            latex2exp::TeX(paste0("Cronbach's $\\alpha$ measures internal consistency")),
            latex2exp::TeX(paste0("$\\hat{\\sigma}$ = Sample Standard Deviation")),
            latex2exp::TeX(paste0(ifelse(ED_included, "◆ = Cam's Self-Assessment", "")))
            )
            
txt_cex <- 2
for(i in 1:length(t2w)){
  text(t2w[[i]], pos = 4,
       y = yloc, srt = 0, cex = txt_cex, x = 1, 
       col = 1, font = 2)
  yloc <- yloc - strheight(latex2exp::TeX("Legend"), cex = txt_cex) * 1.2
  
}


mean_loc <- 5
item_90CI_beta_locs <- c(4,6)
item_90CI_betabinom_locs <- c(2,8)

ybox_locs <- c(par("usr")[3] - diff(par("usr")[3:4])/2.625, 
               par("usr")[3] - diff(par("usr")[3:4])/2.875)
segments(x0 = item_90CI_betabinom_locs[1], x1 = item_90CI_betabinom_locs[2], 
         y0 = mean(ybox_locs), y1 = mean(ybox_locs),
         lwd = 2)
rect(xleft = item_90CI_beta_locs[1], xright = item_90CI_beta_locs[2], 
     ybottom = ybox_locs[1], ytop = ybox_locs[2],
     lwd = 2, col = "white")
segments(x0 = mean_loc, x1 = mean_loc, 
         y0 = ybox_locs[1], y1 = ybox_locs[2],
         lwd = 2)

text(x = mean_loc, y = ybox_locs[1] - 0.5, labels = "posterior mean of score", pos = 1)
segments(x0 = mean_loc, x1 = mean_loc, y0 = ybox_locs[1]-0.05, 
         y1 = ybox_locs[1]-0.6, lty = 3, col = "grey", lwd = 2)
text(x = item_90CI_beta_locs[1], y = ybox_locs[1], 
     labels = latex2exp::TeX(paste0("$5^{th}$ posterior quantile of mean score")), pos = 2,
     srt = 90)
text(x = item_90CI_beta_locs[2] + 0.25, y = ybox_locs[1], 
     labels = latex2exp::TeX(paste0("$95^{th}$ posterior quantile of mean score")), pos = 2,
     srt = 90)
text(x = item_90CI_betabinom_locs[1], y = mean(ybox_locs), 
     labels = latex2exp::TeX(paste0("$5^{th}$ posterior predictive quantile of score")), pos = 2,
     srt = 90)
text(x = item_90CI_betabinom_locs[2] + 0.25, y = mean(ybox_locs), 
     labels = latex2exp::TeX(paste0("$95^{th}$ posterior predictive quantile of score")), pos = 2,
     srt = 90)

dev.off()


fr <- d[,(ncol(d)-5):(ncol(d)-2)]

for(i in 1:length(fr)){
  cat(paste0(colnames(fr)[i], "\n\n\n\n"))
  resps <- fr[,i]
  resps <- resps[!is.na(resps)]
  cat(paste0(resps, collapse = "\n\n"))
  cat(paste0("\n\n----\n\n"))
  
  
}

#### compare years ####

# with just two years, we can do a scatterplot
years <- 2023:2024
dats <- lapply(setNames(years, years), function(yi){
  as.data.frame(data.table::fread(paste0(results_dir, yi, "_item_summary_data.csv")))
})
nitems <- nrow(dats[[1]])

#check to make sure we are selecting the same items
equal_labels_across_years <- dats[[1]]$label == dats[[2]]$label
cbind(dats[[1]]$label, dats[[2]]$label)[!equal_labels_across_years,]

#preprocess
use_median <- T
use_50CI <- T
if(use_median){
  centers <- data.frame(x = dats[[1]]$beta_median, y = dats[[2]]$beta_median)
} else {
  centers <- data.frame(x = dats[[1]]$beta_mean, y = dats[[2]]$beta_mean)
}
if(use_50CI){
  intervals_x <- data.frame(x0 = dats[[1]]$beta_50CI_lb, x1 = dats[[1]]$beta_50CI_ub)
  intervals_y <- data.frame(y0 = dats[[2]]$beta_50CI_lb, y1 = dats[[2]]$beta_50CI_ub)
} else {
  intervals_x <- data.frame(x0 = dats[[1]]$beta_90CI_lb, x1 = dats[[1]]$beta_90CI_ub)
  intervals_y <- data.frame(y0 = dats[[2]]$beta_90CI_lb, y1 = dats[[2]]$beta_90CI_ub)
}
delta <- centers$y - centers$x
order_items <- order(delta)
rank_items <- rank(delta)

#plot main scatterplot
if(randomize_data){
  cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/example_output_comparison.pdf"),
            width = 500 / 72, height = 1200 / 72, family = "Crimson Text")
} else {
  cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", paste0(years, collapse = "-"), "_ED-Eval_Results_WAI.pdf"),
            width = 520 / 72, height = 1200 / 72, family = "Crimson Text")
}

layout(cbind(c(1,2)), heights = c(1,1.5))
plot(centers[,1], centers[,2],
     xlab = paste0(years[1], " Results"), ylab = paste0(years[2], " Results"), 
     xlim = range(c(intervals_x, centers$x)),
     ylim = range(c(intervals_y, centers$y)), 
     col = cat_cols[cat_key[dats[[as.character(year)]]$label]],
     cex = 2, cex.lab = 1.25
)

#get CIs in there
segments(x0 = intervals_x$x0, x1 = intervals_x$x1,
         y0 = centers$y, y1 = centers$y, col = cat_cols[cat_key[dats[[as.character(year)]]$label]])
segments(y0 = intervals_y$y0, y1 = intervals_y$y1,
         x0 = centers$x, x1 = centers$x, col = cat_cols[cat_key[dats[[as.character(year)]]$label]])
abline(0,1,lty=2,lwd=2,col="darkgrey")

#plot numbers for labels
text(centers[,1], centers[,2], labels = rank_items, col = "white", cex = 0.7, font = 2)

#color plot for good or bad result
usr <- par("usr")
x1 <- usr[1];  x2 <- usr[2]
y1 <- usr[3];  y2 <- usr[4]
polygon(
  x = c(x1, x1, x2, x2),
  y = c(x1, y2, y2, x2),
  col=adjustcolor("green", 0.05), border=NA
)

polygon(
  x = c(x1, x1, x2, x2),
  y = c(x1, y1, y1, x2),
  col=adjustcolor("red", 0.05), border=NA
)

title(main = paste0(years[1], " to ", years[2], " Comparison"), cex.main = 2, line = 2)

#legend
legend(
  "topleft",
  legend = c("posterior median", "posterior 50% CI", "1-to-1 line",
             "increased", "decreased"),
  col = c("grey30", "grey30", "grey30",
          adjustcolor("green", 0.1), adjustcolor("red", 0.1)),
  pch = c(19, NA, NA, 15, 15),
  pt.cex = c(2, NA, NA, 2, 2),
  lty = c(NA, 1, 2, NA, NA),   # lines for CI and 1-to-1
  lwd = c(NA, 1, 2, NA, NA),   # thickness for CI and 1-to-1
)
box(lwd = 2)

legend(par("usr")[1], par("usr")[4] + max(strheight(trimws(gsub("_|and_Environment", " ", names(cat_cols))[-5]))) * 3, 
       legend = trimws(gsub("_|and_Environment", " ", names(cat_cols))[-5]), 
       col = cat_cols[-5], pch = 15, pt.cex = 2, cex = 1, horiz = TRUE, xpd = TRUE, box.lwd = 0, bty = "n", 
       text.width = strwidth(trimws(gsub("_|and_Environment", " ", names(cat_cols))[-5])))

#get the prob of each delta
beta_diff_prob_pos <- function(s1, s2, n = 1E4){
  b1 <- rbeta(n, s1[1], s1[2])
  b2 <- rbeta(n, s2[1], s2[2])
  mean((b2 - b1) > 0)
}

prob_improvement <- sapply(1:nitems, function(i){
  beta_diff_prob_pos(
    s1 = c(dats[[1]]$beta_shape_1[i], dats[[1]]$beta_shape_2[i]),
    s2 = c(dats[[2]]$beta_shape_1[i], dats[[2]]$beta_shape_2[i])
  )
})
  
#now plot the item label
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
pr_delta <- round(prob_improvement[order_items], 2)
pr_delta[delta[order_items] < 0] <- 1 - pr_delta[delta[order_items] < 0]
pr_delta <- as.character(pr_delta)
pr_delta[nchar(pr_delta) == 3] <- paste0(pr_delta[nchar(pr_delta) == 3], "0")
text_to_plot <- data.frame(rank = rank_items[order_items], 
                           label = dats[[1]]$label[order_items],
                           col = cat_cols[cat_key[dats[[as.character(year)]]$label]], 
                           delta = delta[order_items],
                           pr_delta = pr_delta)

nchar_pl <- 90
curr_h <- 1.14
leading <- 0.02
for(i in 1:nrow(text_to_plot)){
  points(-0.1, curr_h, 
         col = cat_cols[cat_key[text_to_plot$label[i]]], 
         cex = 2, pch = 19, xpd = NA)
  text(-0.1, curr_h, labels = text_to_plot$rank[i], col = "white", xpd = NA, cex = 0.75)
  
  #calculate text to plot
  comms <- paste0("(", round(text_to_plot$delta[i], 2), ") ", text_to_plot$label[i],
                  " (Pr(Δ) ≈ ", text_to_plot$pr_delta[i],")")
  comms_parts <- str_split(comms, "\n")[[1]]
  comms_wrapped <- sapply(comms_parts, str_wrap, width = nchar_pl)
  comms <- paste(comms_wrapped, collapse = "\n")
  comms_h <- strheight(comms, units = "user", cex = 1)
  
  #plot text
  text(-0.0875, curr_h - comms_h/2, labels = comms, col = 1, xpd = NA, cex = 1, pos = 4)
  
  curr_h <- curr_h - comms_h - leading
}


dev.off()

### DOT PLOT ####

years_to_plot <- 2023:2024
dats <- lapply(setNames(years_to_plot, years_to_plot), function(yi){
  as.data.frame(data.table::fread(paste0(results_dir, yi, "_item_summary_data.csv")))
})

# small helper for device→user coords (for the context label)
device_usr_coords <- function() {
  usr <- par("usr"); plt <- par("plt")
  x_per_ndc <- (usr[2] - usr[1]) / (plt[2] - plt[1])
  y_per_ndc <- (usr[4] - usr[3]) / (plt[4] - plt[3])
  c(x_min = usr[1] - plt[1] * x_per_ndc,
    x_max = usr[2] + (1 - plt[2]) * x_per_ndc,
    y_min = usr[3] - plt[3] * y_per_ndc,
    y_max = usr[4] + (1 - plt[4]) * y_per_ndc)
}

# expected Likert labels (as used in your script)
likert_levels <- c("Strongly Disagree",
                   "Disagree",
                   "Neither Agree nor Disagree",
                   "Agree",
                   "Strongly Agree")

# try to grab 5 count columns by fuzzy-matching these labels; otherwise build counts from sorted_scores
extract_count_matrix <- function(D) {
  # 1) try header-match for each Likert level
  col_idx <- integer(0)
  for (lv in likert_levels) {
    hits <- grep(lv, names(D), ignore.case = TRUE, fixed = FALSE)
    if (length(hits) == 0) { col_idx <- c(col_idx, NA_integer_); next }
    # prefer columns that look numeric across rows
    if (length(hits) > 1) {
      num_like <- sapply(hits, function(h) mean(grepl("^[0-9]+$", as.character(D[[h]])))) # crude score
      hits <- hits[which.max(num_like)]
    }
    col_idx <- c(col_idx, hits[1])
  }
  
  if (!any(is.na(col_idx))) {
    counts <- as.matrix(D[, col_idx, drop = FALSE])
    storage.mode(counts) <- "integer"
    counts[is.na(counts)] <- 0L
    return(list(counts = counts, source = "headers"))
  }
  
  # 2) fallback: rebuild counts from the "sorted_scores..." column
  sc_col <- grep("^sorted_scores", names(D), ignore.case = TRUE, value = TRUE)
  if (!length(sc_col)) stop("Could not find count columns or a 'sorted_scores' column to reconstruct counts.")
  sc_col <- sc_col[1]
  
  # parse each row’s semicolon-separated scores and count 1..5
  counts <- t(sapply(D[[sc_col]], function(s) {
    if (is.na(s) || !nzchar(s)) return(rep.int(0L, 5))
    toks <- unlist(strsplit(s, ";"))
    toks <- trimws(toks)
    toks <- toks[nzchar(toks)]
    toks <- toks[toks %in% c("1","2","3","4","5")]
    tab  <- tabulate(as.integer(toks), nbins = 5)
    as.integer(tab)
  }))
  colnames(counts) <- likert_levels
  return(list(counts = counts, source = "sorted_scores"))
}

# colors
lhs_col <- "#BE4D25"; rhs_col <- "#2596BE"

dot_outfile <- paste0(results_dir, min(years_to_plot), "-", max(years_to_plot), "_dot-plot_multi.pdf")
cairo_pdf(dot_outfile, width = 25, height = 6 * length(years_to_plot), family = "Crimson Text")
on.exit(try(dev.off(), silent = TRUE), add = TRUE)
par(mar = c(1, 75, 4, 2), xpd = TRUE, mfrow = c(length(years_to_plot), 1))
cat("[dot plot] writing -> ", dot_outfile, "\n")

for (yi in years_to_plot) {
  D <- dats[[as.character(yi)]]
  
  # counts (matrix n_items x 5), robustly obtained
  ext <- extract_count_matrix(D)
  counts <- ext$counts
  cat("[dot plot]", yi, "counts source:", ext$source, "\n")
  
  # 50% CI + median (fallbacks if needed)
  if ("beta_median" %in% names(D)) {
    pmeans <- as.numeric(D$beta_median)
  } else if ("beta_mean" %in% names(D)) {
    pmeans <- as.numeric(D$beta_mean)
  } else {
    pmeans <- rowMeans(D[, c("beta_50CI_lb","beta_50CI_ub")], na.rm = TRUE)
  }
  if (!all(c("beta_50CI_lb","beta_50CI_ub") %in% names(D))) {
    CIs <- cbind(pmeans, pmeans)  # degenerate if missing
  } else {
    CIs <- as.matrix(D[, c("beta_50CI_lb","beta_50CI_ub")])
  }
  
  # order by posterior median (desc)
  ord <- order(pmeans, decreasing = TRUE)
  labs   <- as.character(D$label)[ord]
  pmeans <- pmeans[ord]
  CIs    <- CIs[ord, , drop = FALSE]
  counts <- counts[ord, , drop = FALSE]
  
  n_items <- nrow(counts)
  n_cat   <- ncol(counts)  # should be 5
  
  # display labels (wrap only the middle category for readability)
  disp_labs <- colnames(counts)
  disp_labs[disp_labs == "Neither Agree nor Disagree"] <- "Neither Agree\nnor Disagree"
  
  # category colors
  dp_cols <- colorRampPalette(c(lhs_col, rhs_col))(n_cat)
  
  # geometry (width driven by max dots in any item/category)
  resp_per_block <- max(counts, na.rm = TRUE); if (!is.finite(resp_per_block) || resp_per_block < 1) resp_per_block <- 1
  block_w  <- resp_per_block + 1
  left_pad <- 2
  row_ht   <- 0.8
  
  plot(NA, NA,
       xlim = c(0, left_pad + block_w * n_cat),
       ylim = c(0.5 * row_ht, (n_items + 0.5) * row_ht),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  pusr <- par("usr")
  
  # vertical boundaries + headers
  vlines  <- left_pad + seq(0, block_w * n_cat, by = block_w)
  centers <- left_pad + (0:(n_cat - 1)) * block_w + block_w / 2
  text(x = centers, y = pusr[4] + diff(pusr[3:4]) / 20, labels = disp_labs, cex = 0.9)
  text(x = pusr[1], y = pusr[4] + diff(pusr[3:4]) / 30, labels = yi, cex = 2, pos = 2)
  
  # left axis (item text)
  y_pos <- seq_len(n_items) * row_ht
  lax_labs <- labs
  inv_ind <- labs %in% inverted_qs
  lax_labs[inv_ind] <- paste0("Invert(", lax_labs[inv_ind], ")")
  y_ind <- grepl("2024", x = names(cat_key))
  year_entry <- cat_key[y_ind]
  names(year_entry) <- gsub("2024|2024,", paste0(yi, ","), names(year_entry))
  cat_key_ext <- c(cat_key, year_entry)
  entry_cols <- cat_cols[cat_key_ext[labs]]
  text(x = vlines[1], y = y_pos, labels = lax_labs, pos = 2, cex = 1,
       col = entry_cols)
  
  #mark the likert scale inversionas
  points(x = rep(tail(vlines, 1) + diff(vlines)[1]/10, sum(inv_ind)), y = y_pos[inv_ind], pch = "◆", cex = 2)
  
  # guides
  abline(v = vlines, col = rgb(0.6, 0.6, 0.6, 0.5), lwd = 3, xpd = FALSE)
  segments(x0 = vlines[1], x1 = tail(vlines, 1), y0 = y_pos, y1 = y_pos,
           col = adjustcolor(entry_cols, 0.5), lwd = 2, lty = 2, xpd = FALSE)
  
  if(yi == years_to_plot[1]){
    legend(par("usr")[1] - diff(par("usr")[1:2]) * 1.15, par("usr")[4] + diff(par("usr")[3:4])/6, 
           box.lwd = 0, col = cat_cols[-5], legend = gsub("_|and_Environment", " ", names(cat_cols))[-5], 
           pch = 15, pt.cex = 2, cex = 1.5, bty = "n", horiz = F, text.width = 20, ncol = 2)
    
  }
  
  # map score in [1..5] to x via centers (for CI/median overlay)
  map_score_to_x <- function(score_vec) approx(x = 1:n_cat, y = centers, xout = score_vec)$y
  
  # dots (counts) + 50% CI line + median dot
  dot_cex <- 1.5
  for (i in seq_len(n_items)) {
    for (cat_i in seq_len(n_cat)) {
      n <- counts[i, cat_i]
      if (!is.finite(n) || n <= 0) next
      x <- left_pad + (cat_i - 1) * block_w + (block_w - n) / 2 + seq_len(n) - 0.5
      y <- rep(y_pos[i], n)
      points(x, y, pch = 19, cex = dot_cex, col = dp_cols[cat_i])
    }
    # x_ci <- map_score_to_x(CIs[i, ])
    # x_mu <- map_score_to_x(pmeans[i])
    # segments(x0 = x_ci[1], x1 = x_ci[2], y0 = y_pos[i], y1 = y_pos[i], col = "grey30", lwd = 2)
    # points(x_mu, y_pos[i], pch = 19, cex = 0.9, col = "grey10")
  }
  
  # small context note on first panel
  if (yi == years_to_plot[1]) {
    dpusr <- device_usr_coords()
    text("     = Likert Scale Inverted",
         x = dpusr["x_min"], y = dpusr["y_max"] - diff(dpusr[c("y_min","y_max")]) / 20,
         pos = 4, cex = 1.5)
    points(x = dpusr["x_min"] + diff(dpusr[c("x_min", "x_max")]) / 80, 
           y = dpusr["y_max"] - diff(dpusr[c("y_min","y_max")]) / 21, pch = "◆", cex = 2.5)
    
  }
}
dev.off()

### RANKED ITEM SLOPEGRAPH ####
base_col <- "#323840"
add_average_scores <- T
width_lab <- 1.35
heigh_lab_ws <- 0.35

# we already built 'dats' above when you made the two-year scatter; re-use it
# dats[["2023"]] and dats[["2024"]] each have columns: label, beta_median, etc.
stopifnot(all(c("2023", "2024") %in% names(dats)))

m23 <- setNames(dats[["2023"]][,"beta_median"], dats[["2023"]][,"label"])
m24 <- setNames(dats[["2024"]][,"beta_median"], dats[["2024"]][,"label"])
items <- intersect(names(m23), names(m24))
stopifnot(length(items) >= 2)  # unit test: need at least 2 common items

m23 <- m23[items]; m24 <- m24[items]
items23 <- names(sort(m23, decreasing = TRUE))
items24 <- names(sort(m24, decreasing = TRUE))
y23 <- setNames(seq_along(items23), items23)
y24 <- setNames(seq_along(items24), items24)

d_rank <- y23[items] - y24[items]
Dmax   <- max(1L, max(abs(d_rank)))

# gradient: blue (worse rank) → grey (no change) → orange (improved)
col_fun <- function(d) {
  if (!length(d)) return(character(0))
  t <- (pmax(-Dmax, pmin(Dmax, d)) + Dmax) / (2 * Dmax)
  rgb_mat <- colorRamp(c("#2E718F", "grey50", "#D95D26"))(t)
  rgb(rgb_mat[,1]/255, rgb_mat[,2]/255, rgb_mat[,3]/255)
}

rank_outfile <- paste0(results_dir, "2023-2024_rank-slope.pdf")
cairo_pdf(rank_outfile, width = 28, height = max(6, length(items) * 0.4), family = "Crimson Text")
on.exit(try(dev.off(), silent = TRUE), add = TRUE)

par(mar = c(2, 24, 5, 34), xpd = NA)
xL <- 0; xR <- 1
plot(NA, NA, xlim = c(-0.5, 2), ylim = c(0.5, length(items) + 0.5),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

# left labels (2023): "rank. item"
labs_map <- items23
inv_ind <- labs_map %in% inverted_qs
labs_map[inv_ind] <- paste0("Invert(", labs_map[inv_ind], ")")
labs_map <- setNames(labs_map, items23)

labs23 <- paste0(seq_along(items23), ". ", gsub(":", ".", labs_map[items23], fixed = TRUE))
labs23 <- txtwrap(labs23, width = width_lab, cex = 0.9)
h23    <- strheight(labs23, cex = 0.9) + 0.35
yloc23 <- setNames(length(items) - cumsum(c(0, (h23[-length(h23)] + h23[-1]) / 2)), items23)

# right labels (2024): "(Δ) rank. item", colored by Δ
dr_vec <- as.numeric(y23[items24] - y24[items24])
sgn_txt <- ifelse(dr_vec > 0, paste0("+", dr_vec), ifelse(dr_vec < 0, as.character(dr_vec), "0"))
labs24  <- paste0("(", sgn_txt, ") ", seq_along(items24), ". ", gsub(":", ".", labs_map[items24], fixed = TRUE))
labs24  <- txtwrap(labs24, width = width_lab, cex = 0.9)
h24     <- strheight(labs24, cex = 0.9) + 0.35
yloc24  <- setNames(length(items) - cumsum(c(0, (h24[-length(h24)] + h24[-1]) / 2)), items24)
yloc24  <- yloc24 - (h24[1] - min(h24)) / 2
col24   <- col_fun(dr_vec)

text(xL - 0.01, yloc23, labels = labs23, adj = 1, cex = 0.9, col = "#323840")
text(xR + 0.01, yloc24, labels = labs24, adj = 0, cex = 0.9, col = col24)

# connectors
for (nm in items) {
  bezier_curve(x0 = xL, y0 = yloc23[[nm]],
               x1 = xR, y1 = yloc24[[nm]],
               col = col_fun(d_rank[[nm]]), lwd = 2.5, lty = 1)
}

usr <- par("usr")
text(xL - 0.75, usr[4] + 0.7, labels = "2023", cex = 2.0, font = 2, col = "#323840")
text(xR + 0.75, usr[4] + 0.7, labels = "2024", cex = 2.0, font = 2, col = "#323840")

# Δrank legend bar
lg_x0 <- 0.35; lg_x1 <- 0.65; lg_y0 <- usr[4] + 0.25; lg_y1 <- usr[4] + 0.75
xs <- seq(lg_x0, lg_x1, length.out = 2 * Dmax + 2)
cols_leg <- col_fun(seq(-Dmax, Dmax, length.out = length(xs)))
for (i in 1:(length(xs) - 1)) rect(xs[i], lg_y0, xs[i + 1], lg_y1, col = cols_leg[i], border = NA)
rect(lg_x0, lg_y0, lg_x1, lg_y1, border = "grey30")
text(c(lg_x0, (lg_x0 + lg_x1) / 2, lg_x1), lg_y0 - 0.35, labels = c(paste0("-", Dmax), "0", paste0("+", Dmax)), cex = 0.9)
text((lg_x0 + lg_x1) / 2, lg_y1 + 0.4, expression(Delta * " rank"), cex = 1.0)

#add average scores
#add average scores
if(add_average_scores){
  ## short colored stubs extending from 2024 labels
  w24 <- strwidth(labs24, cex = 0.9) + 0.025
  x_right <- max(xR + 0.01 + w24)
  segments(y0 = yloc24, y1 = yloc24,
           x0 = xR + 0.01 + w24,
           x1 = x_right,
           col = col24, lwd = 2.5)
  
  ## position the vertical "average score (2024)" axis
  usr <- par("usr")
  axis_x <- x_right + 0.4
  axis_y <- range(unname(yloc24))
  
  ## draw the axis
  segments(x0 = axis_x, x1 = axis_x,
           y0 = axis_y[1], y1 = axis_y[2],
           col = base_col, lwd = 2.5)
  
  ## define 1–5 Likert mapping (labels on left of axis)
  tick_l <- 0.03
  yvals  <- 1:5
  likert_labels <- c("Strongly Disagree",
                     "Disagree",
                     "Neither Agree nor Disagree",
                     "Agree",
                     "Strongly Agree")
  
  ylocfunc <- function(score){
    # clamp and map [1,5] → [axis_y[1], axis_y[2]]
    s <- pmin(5, pmax(1, score))
    ( (s - min(yvals)) / (max(yvals) - min(yvals)) ) * diff(axis_y) + axis_y[1]
  }
  
  # major ticks + labels
  labx <- axis_x
  laby <- ylocfunc(yvals)
  segments(x0 = axis_x + tick_l, x1 = axis_x, y0 = laby, y1 = laby,
           col = base_col, lwd = 2)
  text(x = labx + tick_l, y = laby, labels = likert_labels, pos = 4, cex = 0.9, col = base_col)
  
  # minor ticks at half-steps
  minor <- setdiff(seq(1, 5, by = 0.5), yvals)
  if(length(minor)){
    segments(x0 = axis_x + tick_l/2, x1 = axis_x,
             y0 = ylocfunc(minor), y1 = ylocfunc(minor),
             col = adjustcolor(base_col, 0.8), lwd = 1)
  }
  
  ## connect each item’s 2024 label to its average score on the axis
  for (nm in items) {
    s <- m24[[nm]]                 # posterior median on 1–5 scale from CSV
    if (is.na(s)) next
    y_target <- ylocfunc(s)
    bezier_curve(x0 = x_right, y0 = yloc24[[nm]],
                 x1 = axis_x - 0.025, y1 = y_target,
                 col = col_fun(d_rank[[nm]]), lwd = 2.5,
                 p = 0.5, n = 128, k = 1, ends = "flat")
  }
}


# context note
x_dev_left   <- grconvertX(0, from = "ndc", to = "user")
y_dev_bottom <- grconvertY(0, from = "ndc", to = "user")
dx <- grconvertX(0.01, from = "ndc", to = "user") - grconvertX(0, from = "ndc", to = "user")
dy <- grconvertY(0.015, from = "ndc", to = "user") - grconvertY(0, from = "ndc", to = "user")
text(x = x_dev_left + dx, y = y_dev_bottom + dy,
     labels = "(ranked by posterior median)", adj = c(0, 0), cex = 0.8, col = adjustcolor(1, 0.65))

dev.off()

#### END ####

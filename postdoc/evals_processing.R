library(extrafont)
library(stringr)
library(psych)
library(corrplot)


# Install and load necessary packages
library(extrafont)

#### read in data ####
randomize_data <- F
base_dir <- "~/Documents/Documents - nikolai/cam_review/"
results_dir <- paste0(base_dir, "results/")
input_dir <- paste0(base_dir, "input/")
# input_name <- c("2023 WAI Executive Director Evaluation (Responses) - Form Responses 1.csv")[1]
year <- 2023
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
            width = 500 / 72, height = base_h, family = "Crimson Text")
  
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
  
  nchar_pl <- 115
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
            width = 500 / 72 * 6, height = max(col_dat$base_h + col_dat$comms_h), family = "Crimson Text")
} else {
  cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", year, "_ED-Eval_Results_WAI.pdf"),
            width = 500 / 72 * 6, height = max(col_dat$base_h + col_dat$comms_h), family = "Crimson Text")
}
par(xpd = NA, mar = c(0,0,0,0))

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

par(mar = c(12,5,5,2))
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

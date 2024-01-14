library(extrafont)
library(stringr)
library(psych)
library(corrplot)


# Install and load necessary packages
library(extrafont)

#### functions ####
text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}
#### read in data ####

d <- read.csv("~/Documents/Documents - nikolai/cam_review/2023 WAI Executive Director Evaluation (Responses) - Form Responses 1.csv")
likert <- c("Strongly Disagree", "Disagree", "Neither Agree nor Disagree" ,"Agree" ,"Strongly Agree")
d[d == "N/A"] <- NA
d[d == ""] <- NA
for(ri in likert){
  d[d == ri] <- which(likert == ri)
}
col.names <- c("Timestamp",
               "The ED clearly communicates relevant expectations and goals to me and my team.",
               "The ED responds to requests or communications in a timely manner.",
               "The ED understands and helps me to understand how my work contributes to organizational strategy.",
               "The Executive Director clearly communicates motivations and justifications for their decisions. I understand why they make the choices they do, even when I don’t agree with them.",
               "I have a clear understanding of how the ED perceives my work performance at WAI.",
               "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
               "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
               "In 2023, the ED made excellent progress toward the ORCAs for which they are primarily responsible.",
               "The ED does not trust my expertise in my specific work area, and tends to micromanage projects that would better be left to me.",
               "I am inspired by the ED’s efforts to boost team morale and spirit. Their enthusiasm seems sincere and resonates with me.",
               "I often disagree with decisions made by the ED in my work area.",
               "I would feel comfortable approaching the ED directly with suggestions on how they could improve in their role at WAI.",
               "I would feel comfortable approaching the ED directly about something impeding my success and happiness at WAI that they are not directly responsible for.",
               "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
               "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
               "The ED is sufficiently familiar with the technical aspects of wild animal welfare science to be able to provide meaningful input and direction for research at WAI.",
               "The ED is effective at time management – they appropriately prioritize important tasks and are able judge less important ones to ignore or delegate.",
               "I am satisfied with the amount of time I interact with the ED.",
               "The ED provides regular feedback in ways that effectively facilitate my own improvement, both in the form of constructive criticism and in positive affirmation of jobs well done.",
               "The ED has been effective at obtaining funding for WAI’s continued operation and growth.",
               "I am proud of how the ED represents WAI when engaging with external stakeholders.",
               "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
               "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
               "The ED fosters an inclusive work environment, both in supporting staff and in intervening against intolerant work dynamics.",
               "I have been made uncomfortable by the ED in the last 12 months.",
               "The ED is proactive in identifying ways in which the WAI work environment can be improved.",
               "Disagreements with the ED about approach, mission, and tactics are handled appropriately and respectfully.",
               "The ED successfully moderates disagreements among staff.",
               "The ED treats all employees equitably and without favoritism.",
               "Please use this space for comments on the above statements or your selections that will be shared with Cam directly:",
               "Please use this space for comments on the above statements or your selections that will be privately visible only to the board:",
               "In what important areas has the ED been most successful, relative to your expectations for their role? In what ways would they likely outperform a different ED?",
               "In what important areas has the ED been least successful, relative to your expectations for their role? In what ways would a different ED likely outperform them?",
               "Imagine you have the opportunity to direct how the ED spends 40 hours of their focused time over the next year, without affecting their other responsibilities. How would you allocate these hours to maximize the ED's effectiveness in their role at WAI? Responses can include, but are not limited to, professional development courses, wellness activities, skill-building workshops, or any other accessible experiences that you believe would enhance their performance and contribution to WAI. Please be as specific as possible in your suggestions.",
               "As above, but corresponding to only 1 hour of time.",
               "Please put any comments on the survey as a whole here! Responses here will not be shared with Cam directly and will only be visible to the board. You may also use this space to privately respond to the above questions.",
               "Score")
inverted_qs <- col.names[c(10, 12, 26)]

colnames(d) <- col.names
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


#### some quick analyses of these data ####

#construct table
dn <- apply(d, 2, as.numeric)
dn <- dn[,!apply(apply(dn, 2, is.na), 2, all)]
dn[,inverted_qs] <- 6 - dn[,inverted_qs]

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
comms <- comms[!is.na(comms)]
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

#### one big figure ####

cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/2023-12_ED-Eval_Results_WAI.pdf"),
          width = 500 / 72 * 6, height = max(col_dat$base_h + col_dat$comms_h), family = "Crimson Text")
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
text(gsub("_", " ", sec_name), y = 0.4, x = 0, cex = sec_name_cex[si], col = "#00c1b2")

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
       pos = 3, cex = 1.5)
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
  annot <- latex2exp::TeX(paste0("$\\Delta$ Cronbach's $\\alpha$ = ", 
                                 round(cronbachs_alpha$total$raw_alpha, 2), " - ", 
                                 round(cronbachs_alpha$total$raw_alpha - cronbachs_alpha_diff, 2), 
                                 " $\\approx$ ", round(cronbachs_alpha_diff, 2), 
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
     y = 1.075 + adj_y[si], x = -0.1, cex = txt_cex, font = 2, xpd = NA)
text(comms, y = 1 + adj_y[si] - 
       strheight(comms, units = "user", family = "Crimson Text") * txt_cex * 1.2 / 2 -
       strheight(" ", units = "user", family = "Crimson Text") * txt_cex * 3, 
     x = -1.3, cex = txt_cex, pos = c(4), xpd = NA)


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
         mar = c(0, 0, 15, 15), cl.length = 5, cl.align.text = "l")
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
     labels = sapply(colnames(item_corrmat), truncate_label), xpd = TRUE, xpd = NA, cex = 1.5)


# Generate the corrplot with the custom color palette
corrplot(as.matrix(person_corrmat), method = "square", 
         col = custom_palette(200),  # Use 200 colors for a smooth gradient
         tl.cex = 1.6,              
         tl.col = "black",          
         tl.pos = "n",              
         addgrid.col = NA,          
         cl.cex = 3,                
         mar = c(10, 0, 5, 15), cl.length = 5, cl.align.text = "l")
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
            latex2exp::TeX(paste0("$\\hat{\\sigma}$ = Sample Standard Deviation"))
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

ybox_locs <- c(par("usr")[3] - diff(par("usr")[3:4])/3, 
               par("usr")[3] - diff(par("usr")[3:4])/3.25)
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

#####


fr <- d[,(ncol(d)-5):(ncol(d)-2)]

for(i in 1:length(fr)){
  cat(paste0(colnames(fr)[i], "\n\n\n\n"))
  resps <- fr[,i]
  resps <- resps[!is.na(resps)]
  cat(paste0(resps, collapse = "\n\n"))
  cat(paste0("\n\n----\n\n"))
  
  
}

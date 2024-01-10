library(extrafont)
library(stringr)

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


#### print ####

for(si in 1:4){
  
subds <- ds[[si]][,1:(length(ds[[si]])-1)]

base_h <- 250 * ncol(subds) / 72
cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", names(ds)[si], ".pdf"),
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

cairo_pdf(paste0("~/Documents/Documents - nikolai/cam_review/results/", names(ds)[si], ".pdf"),
          width = 500 / 72, height = base_h + comms_h, family = "Crimson Text")


par(xpd = NA, mar = c(0,0,0,0))
row_heights <- c(0.175, rep(1, ncol(subds)))
row_heights <- c(row_heights, sum(row_heights) * comms_h / base_h)
layout(matrix(1:(length(subds) + 2), ncol = 1), heights = row_heights)

plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
sec_name <- names(ds)[si]
sec_name_cex <- c(5, 5, 4, 4)
text(gsub("_", " ", sec_name), y = 0.4, x = 0, cex = sec_name_cex[si], col = "#00c1b2")

par(mar = c(9,5,5,2))
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
  qlab <- str_wrap(paste0(i, ". ", qtext), width = 90)
  text(label = qlab, 
       x = mean(par("usr")[1:2]), 
       y = par("usr")[4] + diff(par("usr")[3:4])/6 - strheight(qlab),
       pos = 3, cex = 1.5)
  axis(2, at = seq(0, max(counts), by = 1), las = 2, labels = seq(0, max(counts), by = 1), cex.axis = 1.5)
  text(bp, par("usr")[3] - 0.5, srt = 30, adj = 1, labels = likert_names, xpd = TRUE, xpd = NA, cex = 1.5)
  
}


adj_y <- c(0,0,0,0.2)
plot(NA,NA, xlim = c(-1,1), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", xpd = NA, yaxt = "n", frame = F)
text("Please use this space for comments on the above statements or your selections that will be shared with Cam:", 
     y = 1.075 + adj_y[si], x = -0.1, cex = 1.3, font = 2, xpd = NA)
txt_cex <- 1.3
text(comms, y = 1 + adj_y[si] - 
       strheight(comms, units = "user", family = "Crimson Text") * txt_cex * 1.2 / 2 -
       strheight(" ", units = "user", family = "Crimson Text") * txt_cex * 3, 
     x = -1.3, cex = txt_cex, pos = c(4), xpd = NA)


dev.off()

}


fr <- d[,(ncol(d)-5):(ncol(d)-2)]

for(i in 1:length(fr)){
  cat(paste0(colnames(fr)[i], "\n\n\n\n"))
  resps <- fr[,i]
  resps <- resps[!is.na(resps)]
  cat(paste0(resps, collapse = "\n\n"))
  cat(paste0("\n\n----\n\n"))
  
  
}

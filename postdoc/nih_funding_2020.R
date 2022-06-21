library(XML)

#load data
if(!exists("d_nih")){
  d_nih <- read.csv("~/data/nih_awards_2020.csv")
  d_nih$FUNDING <- gsub("\\$", "", d_nih$FUNDING)
  d_nih$FUNDING <- gsub(",", "", d_nih$FUNDING)
  d_nih$FUNDING <- gsub(" ", "", d_nih$FUNDING)
  d_nih$log10FUNDING <- log10(as.numeric(d_nih$FUNDING))
  d_nih <- d_nih[!is.na(d_nih$log10FUNDING),]
  d_nih <- d_nih[d_nih$log10FUNDING != 0,]
}

if(!exists("d_nsf")){
  dir_nsf <- "~/data/nsf_2020/"
  files <- list.files(dir_nsf)
  d_nsf <- data.frame(FUNDING = rep(NA, length(files)))
  for(i in 1:length(files)){
    if(i %% 100 == 0){cat(paste0(i, " "))}
    d_nsf$FUNDING[i] <- as.numeric(xmlToDataFrame(xmlRoot(xmlParse(paste0(dir_nsf, files[i]))))$AwardAmount)
  }
  d_nsf$log10FUNDING <- log10(d_nsf$FUNDING)
  d_nsf <- d_nsf[!is.na(d_nsf$log10FUNDING),]
  d_nsf <- d_nsf[d_nsf$log10FUNDING != 0,]
  d_nsf <- d_nsf[d_nsf$log10FUNDING != -Inf,]
}

par(mfrow = c(2,1))

#nih
hist(d_nih$log10FUNDING, breaks = seq(0, ceiling(max(d_nih$log10FUNDING)), length.out = 100), xaxt = 'n', cex.main = 1.5,
     ylim = c(0, 15000), main = paste0("NIH 2020 Distribution of Funding Amounts\n(", nrow(d_nih), " awards, $", 
                                       round(sum(as.numeric(d_nih$FUNDING)) / 1E9, 2) , " Billion USD in total)"),
     xlab = latex2exp::TeX("Award Amount ($log_{10}(USD)$)"))
axis(side=1, at=seq(0,ceiling(max(d_nih$log10FUNDING)), by = 1), labels=10^seq(0,ceiling(max(d_nih$log10FUNDING)), by = 1))
segments(x1 = log10(250000), x0 = log10(250000), y1 = 14000, y0 = 0, col = 2, lwd = 2)

prop_less <- round(mean(d_nih$log10FUNDING < log10(250000)) * 100)
text(log10(250000), y = 13500, labels = paste0(prop_less, "%"), srt = 90, pos = 2, font = 2)
text(log10(250000), y = 13500, labels = paste0(100-prop_less, "%"), srt = 270, pos = 4, font = 2)
text(log10(250000), y = 14000, labels = "$250,000", srt = 0, pos = 3, font = 3)

#nsf
hist(d_nsf$log10FUNDING, breaks = seq(0, ceiling(max(d_nsf$log10FUNDING)), length.out = 100), xaxt = 'n', cex.main = 1.5,
     ylim = c(0, 2000), main = paste0("NSF 2020 Distribution of Funding Amounts\n(", nrow(d_nsf), " awards, $", 
                                      round(sum(d_nsf$FUNDING) / 1E9, 2) , " Billion USD in total)"),
     xlab = latex2exp::TeX("Award Amount ($log_{10}(USD)$)"))
axis(side=1, at=seq(0,ceiling(max(d_nsf$log10FUNDING)), by = 1), labels=10^seq(0,ceiling(max(d_nsf$log10FUNDING)), by = 1))
segments(x1 = log10(250000), x0 = log10(250000), y1 = 1900, y0 = 0, col = 2, lwd = 2)

prop_less <- round(mean(d_nsf$log10FUNDING < log10(250000)) * 100)
text(log10(250000), y = 1800, labels = paste0(prop_less, "%"), srt = 90, pos = 2, font = 2)
text(log10(250000), y = 1800, labels = paste0(100-prop_less, "%"), srt = 270, pos = 4, font = 2)
text(log10(250000), y = 1900, labels = "$250,000", srt = 0, pos = 3, font = 3)

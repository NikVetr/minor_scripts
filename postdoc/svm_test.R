n <- 3000
p <- 100
train_test_rat <- 0.8
out_b <- 0
dat <- data.frame(preds = matrix(rnorm(n*p, sd = 0.1), n, p), 
                  out = as.factor(sample(rep(c(0,1), n/2))))
dat <- cbind(dat[,1:p] + matrix(as.numeric(dat$out)-1, n, p, byrow = F) * out_b, out = dat$out)

out_inds <- lapply(unique(dat$out), function(bi) which(dat$out == bi))
train_dat <- dat[unlist(lapply(out_inds, function(bis) bis[1:(round(train_test_rat * length(bis)))])), ]
test_dat <- dat[unlist(lapply(out_inds, function(bis) bis[(round(train_test_rat * length(bis) + 1):length(bis))])), ]

fit <- e1071::svm(out ~ ., dat = train_dat, kernel = "radial")

#out of sample
results <- caret::confusionMatrix(data = predict(fit, test_dat[,!"out" == colnames(test_dat)]), reference = test_dat$out)
round(results$overall[c("Accuracy", "AccuracyPValue", "AccuracyLower", "AccuracyUpper")] * 100, 2)

#in sample
results <- caret::confusionMatrix(data = predict(fit, train_dat[,!"out" == colnames(train_dat)]), reference = train_dat$out)
round(results$overall[c("Accuracy", "AccuracyPValue", "AccuracyLower", "AccuracyUpper")] * 100, 2)


#also try LDA
LDA.1 <- suppressWarnings(MASS::lda(train_dat$out ~ ., dat = train_dat))
LDA_pred.1 <- as.matrix(test_dat[,!"out" == colnames(test_dat)]) %*% LDA.1$scaling
LDA_pred.1.g0 <- LDA_pred.1[test_dat$out == 0]
LDA_pred.1.g1 <- LDA_pred.1[test_dat$out == 1]
LDA.1.breaks <- seq2(range(c(LDA_pred.1.g0, LDA_pred.1.g1)) + c(-1E-6,1E-6), length.out = 50)
LDA_pred.1.freq_max <- max(c(hist(LDA_pred.1.g0, breaks = LDA.1.breaks, plot = F)$counts, 
                             hist(LDA_pred.1.g1, breaks = LDA.1.breaks, plot = F)$counts)) * 1.2
hist(LDA_pred.1.g0, breaks = LDA.1.breaks, col = adjustcolor("blue", 0.5),
     ylim = c(0,LDA_pred.1.freq_max), xlab = "Score on LD Axis 1", main = "out \"Uncorrected\" LD Scores\n")
hist(LDA_pred.1.g1, breaks = LDA.1.breaks, col = adjustcolor("orange", 0.5), add = T)
legend(x = "topright", legend = c("Group 1", "Group 2"),
       pch = 15, col = adjustcolor(c("blue","orange"), 0.5), pt.bg = 1)
title(main = paste0("\n(difference in means = ", round(abs(mean(LDA_pred.1.g0) - mean(LDA_pred.1.g1)), 3),")"), col.main = 2)

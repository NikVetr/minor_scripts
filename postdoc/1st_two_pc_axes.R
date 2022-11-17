# set.seed(1)

#functions
seq2 <- function(r, length.out, ...) seq(r[1], r[2], length.out = length.out)

#high level params
n <- 2E2 #number of individuals (eg humans) in the sample
m <- 5E2 #number of analytes / outcomes
p <- 5E1 #number of variables that systematically structure outcomes
pervasive_interaction <- F
batch_heteroskedasticity_ratio <- 1
batch_extension_ratio <- 0
batch_var <- 0
train_test_rat <- 0.8 #80/20 train/test split
b <- matrix(rnorm(m*p), p, m) #true coefficients, effect of variables on each outcome
x <- matrix(rnorm(n*p), n, p) #values for each variable for each individual

#big ol batch effect
x[,1] <- sample(rep(c(0,1), each=n/2))
b[1,] <- rnorm(m, 0, sqrt(batch_var))

#interaction effects with batch
if(pervasive_interaction){
  b <- rbind(b, matrix(rnorm(m*(p-1)), p-1, m))
  x <- cbind(x, apply(x[,-1], 2, function(i) i * (x[,1]-0.5))) #subtract 0.5 to ensure equal variances across batch in generative model
}

#now simulate analytes
y <- x %*% b #y observed without iid noise
e <- matrix(rnorm(n*m, sd = sqrt(p/2)) * 1, n, m) #residual error for y, with ~1/3rd variance indep
e <- e + diag(x[,1]) %*% e * batch_heteroskedasticity_ratio #heteroskedasticity
ye <- y + e #the y we observe, with residual noise

#batch effect depends on analyte value
ye <- ye + diag(x[,1]) %*% ye * batch_extension_ratio

#standardize ye for later analysis
ye_std <- ye %*% diag(1/apply(ye, 2, sd))
ye_std <- ye_std - matrix(apply(ye_std, 2, mean), n, m, byrow = T)

#pre-specify train and test sets
batch_inds <- lapply(unique(x[,1]), function(bi) which(x[,1] == bi))
train_inds <- unlist(lapply(batch_inds, function(bis) bis[1:(round(train_test_rat * length(bis)))]))
test_inds <- setdiff(1:n, train_inds)

#get batch residuals
ye.resids <- do.call(cbind, lapply(1:m, function(i){
  fit <- lm(ye[,i] ~ x[,1])
  fit$residuals
}))
colnames(ye) <- colnames(y) <- colnames(ye.resids) <- colnames(ye_std) <- paste0("y.", 1:m)

#perform PCA
pca.1 <- prcomp(ye)
scaled_PC_scores.1 <- pca.1$x %*% diag(1/pca.1$sdev) #get everything into unit variance
scaled_PC_scores.1 <- scaled_PC_scores.1[,1:(p-1)] #maximum p axes of systematic (non-noise related) variation

pca.2 <- prcomp(ye.resids)
scaled_PC_scores.2 <- pca.2$x %*% diag(1/pca.2$sdev)
scaled_PC_scores.2 <- scaled_PC_scores.2[,1:(p-1)]

#now do the plotting
layout(mat = rbind(c(1,2,5,7,9), c(3,4,6,8,10)))
plot(scaled_PC_scores.1[,1:2], col = adjustcolor(c("blue","orange"), 0.5)[x[,1] + 1],
     xlab = "PC Axis 1", ylab = "PC Axis 2", main = "Batch \"Uncorrected\" PCA", pch = 19)
legend(x = "topleft", legend = c("Group 1", "Group 2"),
       pch = 19, col = adjustcolor(c("blue","orange"), 0.5))
PCA.1.diffs <- sapply(1:ncol(scaled_PC_scores.1), function(i) mean(scaled_PC_scores.1[x[,1] == 0,i]) - mean(scaled_PC_scores.1[x[,1] == 1,i]))
hist.1 <- hist(PCA.1.diffs, breaks = 20, plot = F)
plot(hist.1, xlab = "Difference in Means between Batches on PC Axes", 
     main = "Batch \"Uncorrected\" PCA", col = "grey")
PCA.1.1_i <- which.min(abs((PCA.1.diffs[1] - hist.1$mids)))
points(hist.1$mids[PCA.1.1_i], y = hist.1$density[PCA.1.1_i] + diff(par("usr")[3:4])/ 15, 
       pch = 25, col = 2, bg = 2)
segments(x0 = hist.1$mids[PCA.1.1_i], x1 = hist.1$mids[PCA.1.1_i], 
         y0 = hist.1$density[PCA.1.1_i] + diff(par("usr")[3:4])/ 15, 
         y1 = hist.1$density[PCA.1.1_i] + diff(par("usr")[3:4])/ 5, 
         pch = 25, col = 2, bg = 2, lwd = 2)
lab <- "Difference on first PC"
# text(x = hist.1$mids[PCA.1.1_i] - diff(par("usr")[1:2])/ 30, 
#      y = hist.1$density[PCA.1.1_i] + diff(par("usr")[3:4])/ 4.5,
#      labels = lab, srt = 90, col = 2, xpd = NA, pos = 4)
text(x = hist.1$mids[PCA.1.1_i] - diff(par("usr")[1:2])/ 30, 
     y = hist.1$density[PCA.1.1_i] + diff(par("usr")[3:4])/ 4,
     labels = lab, srt = 45, col = 2, xpd = NA, pos = 4)


plot(scaled_PC_scores.2[,1:2], col = adjustcolor(c("blue","orange"), 0.5)[x[,1] + 1],
     xlab = "PC Axis 1", ylab = "PC Axis 2", main = "Batch \"Corrected\" PCA", pch = 19)
legend(x = "topleft", legend = c("Group 1", "Group 2"),
       pch = 19, col = adjustcolor(c("blue","orange"), 0.5))
PCA.2.diffs <- sapply(1:ncol(scaled_PC_scores.1), function(i) mean(scaled_PC_scores.2[x[,1] == 0,i]) - mean(scaled_PC_scores.2[x[,1] == 1,i]))
hist(PCA.2.diffs, breaks = 20, xlab = "Difference in Means between Batches on PC Axes", 
     main = "Batch \"Corrected\" PCA\n")
title(main = "\n(note scale on horizontal axis)", col.main = 2)

#perform PCA in the train set and project test set onto axes
pca.train.1 <- prcomp(ye[train_inds,])
scaled_PC_scores.train.1 <- pca.train.1$x %*% diag(1/pca.train.1$sdev) #get everything into unit variance
scaled_PC_scores.train.1 <- scaled_PC_scores.train.1[,1:(p-1)] #maximum p axes of systematic (non-noise related) variation

scaled_PC_scores.test.1 <- ye[train_inds,] %*% pca.train.1$rotation
scaled_PC_scores.test.1 <- scaled_PC_scores.test.1 %*% diag(1/pca.train.1$sdev) #get everything into unit variance
scaled_PC_scores.test.1 <- scaled_PC_scores.train.1[,1:(p-1)] #maximum p axes of systematic (non-noise related) variation

pca.train.2 <- prcomp(ye.resids[train_inds,])
scaled_PC_scores.train.2 <- pca.train.2$x %*% diag(1/pca.train.2$sdev) #get everything into unit variance
scaled_PC_scores.train.2 <- scaled_PC_scores.train.2[,1:(p-1)] #maximum p axes of systematic (non-noise related) variation

scaled_PC_scores.test.2 <- ye[train_inds,] %*% pca.train.2$rotation
scaled_PC_scores.test.2 <- scaled_PC_scores.test.2 %*% diag(1/pca.train.2$sdev) #get everything into unit variance
scaled_PC_scores.test.2 <- scaled_PC_scores.train.2[,1:(p-1)] #maximum p axes of systematic (non-noise related) variation

#now residualize the train and test sets seperately too
ye.resids.train.test <- do.call(cbind, lapply(1:m, function(i){
  fit <- lm(ye[train_inds,i] ~ x[train_inds,1])
  c(train_resids = fit$residuals, 
    test_resids = ye[test_inds,i] - (coef(fit)[1] + coef(fit)[2] * x[test_inds,1]))
}))
ye.resids.train <- ye.resids.train.test[grep("train", rownames(ye.resids.train.test)),]
ye.resids.test <- ye.resids.train.test[grep("test", rownames(ye.resids.train.test)),]
ye.resids.train.std <- sapply(1:m, function(ci) (ye.resids.train[,ci] - mean(ye.resids.train[,ci])) / sd(ye.resids.train[,ci]))
ye.resids.test.std <- sapply(1:m, function(ci) (ye.resids.test[,ci] - mean(ye.resids.train[,ci])) / sd(ye.resids.train[,ci]))

#also standardize original ye
ye.train <- ye[train_inds,]
ye.test <- ye[test_inds,]
ye.train.std <- sapply(1:m, function(ci) (ye.train[,ci] - mean(ye.train[,ci])) / sd(ye.train[,ci]))
ye.test.std <- sapply(1:m, function(ci) (ye.test[,ci] - mean(ye.train[,ci])) / sd(ye.train[,ci]))

#confirm LDA does not work
LDA.1 <- suppressWarnings(MASS::lda(as.factor(x[train_inds,1]) ~ ye.train.std))
LDA_pred.1 <- ye.test.std %*% LDA.1$scaling
LDA_pred.1.g0 <- LDA_pred.1[x[test_inds,1] == 0]
LDA_pred.1.g1 <- LDA_pred.1[x[test_inds,1] == 1]
LDA.1.breaks <- seq2(range(c(LDA_pred.1.g0, LDA_pred.1.g1)) + c(-1E-6,1E-6), length.out = 50)
LDA_pred.1.freq_max <- max(c(hist(LDA_pred.1.g0, breaks = LDA.1.breaks, plot = F)$counts, 
                             hist(LDA_pred.1.g1, breaks = LDA.1.breaks, plot = F)$counts)) * 1.2
hist(LDA_pred.1.g0, breaks = LDA.1.breaks, col = adjustcolor("blue", 0.5),
     ylim = c(0,LDA_pred.1.freq_max), xlab = "Score on LD Axis 1", main = "Batch \"Uncorrected\" LD Scores\n")
hist(LDA_pred.1.g1, breaks = LDA.1.breaks, col = adjustcolor("orange", 0.5), add = T)
legend(x = "topright", legend = c("Group 1", "Group 2"),
       pch = 15, col = adjustcolor(c("blue","orange"), 0.5), pt.bg = 1)
title(main = paste0("\n(diff in means = ", round(abs(mean(LDA_pred.1.g0) - mean(LDA_pred.1.g1)), 1),", log10p-val = ",
                    round(log10(t.test(LDA_pred.1.g0, LDA_pred.1.g1)$p.value), 1), ")"), col.main = 2)

#now the residualized LDA
LDA.2 <- suppressWarnings(MASS::lda(as.factor(x[train_inds,1]) ~ ye.resids.train.std))
LDA_pred.2 <- ye.resids.test.std %*% LDA.2$scaling
LDA_pred.2.g0 <- LDA_pred.2[x[test_inds,1] == 0]
LDA_pred.2.g1 <- LDA_pred.2[x[test_inds,1] == 1]

LDA.2.breaks <- seq2(range(c(LDA_pred.2.g0, LDA_pred.2.g1)) + c(-1E-6,1E-6), length.out = 50)
LDA_pred.2.freq_max <- max(c(hist(LDA_pred.2.g0, breaks = LDA.2.breaks, plot = F)$counts, 
                             hist(LDA_pred.2.g1, breaks = LDA.2.breaks, plot = F)$counts)) * 1.2
hist(LDA_pred.2.g0, breaks = LDA.2.breaks, col = adjustcolor("blue", 0.5), 
     ylim = c(0,LDA_pred.2.freq_max), xlab = "Score on LD Axis 1", main = "Batch \"Corrected\" LD Scores\n")
legend(x = "topright", legend = c("Group 1", "Group 2"),
       pch = 15, col = adjustcolor(c("blue","orange"), 0.5), pt.bg = 1)
title(main = paste0("\n(diff in means = ", round(abs(mean(LDA_pred.2.g0) - mean(LDA_pred.2.g1)), 1),", log10p-val = ",
                    round(log10(t.test(LDA_pred.2.g0, LDA_pred.2.g1)$p.value), 1), ")"), col.main = 2)
hist(LDA_pred.2.g1, breaks = LDA.2.breaks, col = adjustcolor("orange", 0.5), add = T)

#now let's try SVMs
dat.1 <- list(train = data.frame(batch = as.factor(x[train_inds,1]), analytes = ye[train_inds,]),
              test = data.frame(batch = as.factor(x[test_inds,1]), analytes =  ye[test_inds,]))
dat.2 <- list(train = data.frame(batch = as.factor(x[train_inds,1]), analytes = ye.resids.train),
              test = data.frame(batch = as.factor(x[test_inds,1]), analytes =  ye.resids.test))
random_batch <- sample(x[,1])
dat.3 <- list(train = data.frame(batch = as.factor(random_batch[train_inds]), analytes = ye[train_inds,]),
              test = data.frame(batch = as.factor(random_batch[test_inds]), analytes =  ye[test_inds,]))
dat.4 <- list(train = data.frame(batch = as.factor(random_batch[train_inds]), analytes = ye.resids.train),
              test = data.frame(batch = as.factor(random_batch[test_inds]), analytes =  ye.resids.test))
dats <- list(uncorrected = dat.1, 
             corrected = dat.2, 
             uncorrected.random_batch_label = dat.3, 
             corrected.random_batch_label = dat.4)

tune_cost <- F
printp <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
svm_accuracy <- parallel::mclapply(seq_along(dats), function(dat_i){
  data.frame(t(sapply(c(linear = "linear", radial = "radial"), function(kern){
    printp(paste0(names(dats)[[dat_i]], ", ", kern, "\n"))
    train_dat <- dats[[dat_i]]$train
    test_dat <- dats[[dat_i]]$test
    if(tune_cost){
      tune.res <- tune(svm, batch ~ ., data = dat ,kernel = kern, 
                       ranges = list(cost=10^(-3:3)))
      best.cost <- as.numeric(tune.res$best.parameters)
    } else {
      best.cost <- 1
    }
    fit <- e1071::svm(batch ~ ., data = train_dat, kernel = kern, cost = best.cost)
    results <- caret::confusionMatrix(reference = test_dat$batch, 
                                      data = predict(fit, test_dat[,!"batch" == colnames(test_dat)]))
    results$overall[c("Accuracy", "AccuracyPValue", "AccuracyLower", "AccuracyUpper")]
  })))
}, mc.cores = 4)
names(svm_accuracy) <- names(dats)


for(i in 1:length(svm_accuracy)){
  plot(x = 1:2, y = svm_accuracy[[i]]$Accuracy, ylim = c(0,1), xlab = "Kernel", 
       xaxt = "n", xlim = c(0.5, 2.5), ylab = "Accuracy", pch = 19,
       main = paste0("SVM, ", tools::toTitleCase(gsub("\\.", "\n", gsub("_", " ", names(svm_accuracy)[i])))))
  segments(x0 = 1:2, x1 = 1:2, y0 = svm_accuracy[[i]]$AccuracyLower, y1 = svm_accuracy[[i]]$AccuracyUpper, xpd = NA, lwd = 2)
  abline(h = 0.5, col = adjustcolor(1, 0.5), lty = 2)
  segments(x0 = 1:2, x1 = 1:2, y0 = c(-0.04,-0.04), y1 = c(-0.08,-0.08), xpd = NA)
  text(x = 1:2, y = c(-0.08,-0.08), labels = c("linear", "radial"), xpd = NA, pos = 1)
  text(x = 1:2, y = svm_accuracy[[i]]$Accuracy, pos = 4,
       labels = round(log10(svm_accuracy[[i]]$AccuracyPValue), 1), 
       xpd = NA, col = 2, font = 2)
  legend(x = "bottomright", legend = c("Accuracy", "95% CI"),
         pch = c(19, NA), lwd = c(NA, 2))
  text(x = mean(par("usr")[1]), y = par("usr")[3] + diff(par("usr")[3:4])/30, 
       latex2exp::TeX("\\textbf{# is log$_{10}$(p-value)}"), pos = 4, xpd = NA, col = 2)
}


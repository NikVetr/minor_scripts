
#### load data ####
if(!exists("x")){
  x <- read.csv("~/data/transplant_data/proteinGroupsMeta - Corrected_removed 6AT2_changed5NT2toT3.csv")
  x <- x[x$Group != "",]  
}

table(x$TimePoint)
prots <- unique(x$proteinNames)
genes <- unique(x$gene)
groups <- unique(x$Group)
prot2gene <- as.data.frame(do.call(rbind, strsplit(unique(paste0(x$proteinNames, "~", x$gene)), "~")))
colnames(prot2gene) <- c("prot", "gene")
unique(x$majorityProteinNames)
plot(quantile(log10(x$intensity), 0:10000/10000, na.rm = T), type = "l")
mean(x$intensity < 1E2, na.rm = T)
head(sort(x$intensity.total), n = 1000)
head(sort(x$intensity), n = 1000)
rle(sort(x$intensity))

table(x$TimePoint)
table(x$Patient)
table(x$BioMaterial)
summary(x)


hist(x$intensity)
mean(x$intensity == 0, na.rm = T)


### test out this VGAM package ####
# library(VGAM)
# 
# N <- 10
# sex <- sample(0:1, N, T)
# cdata <- data.frame(y = rpois(N, exp(1.25 + sex * 0.25)), sex = sex)
# L <- round(quantile(cdata$y, 0.25))
# cdata <- transform(cdata,
#                    cY = pmax(L, y),
#                    lcensored = y <  L)  # Note y < L, not cY == L or y <= L
# cdata <- transform(cdata, status = ifelse(lcensored, 0, 1))
# with(cdata, table(cY))
# with(cdata, table(lcensored))
# with(cdata, table(print(SurvS4(cY, status, type = "left"))))  # Check
# cdata <- cbind(cdata)
# fit <- vglm(SurvS4(cY, status, type = "left") ~ 1 + sex, cens.poisson,
#             data = cdata, trace = TRUE)
# coef(fit, matrix = TRUE)
# fit@misc$p
# fits <- summary(fit)
# fits
# fits@coef3[,4]
# (1 - pnorm(fits@coef3[,3] * sign(fits@coef3[,3]))) * 2


#### munge into desired data frame ####
L <- 44780
U <- Inf
d <- data.frame(count = x$intensity)
d$cens <- 1
d$cens[d$count == 0] <- 0
d$count[d$count == 0] <- L
d$prot <- match(x$proteinNames, prots)
d$time <- x$TimePoint
d$group <- match(x$Group, groups)
# d$time1 <- as.integer(d$time >= 1)
# d$time2 <- as.integer(d$time >= 2)
# d$time3 <- as.integer(d$time >= 3)
d$group1 <- as.integer(d$group == 1)
d$group2 <- as.integer(d$group == 2)
d$group3 <- as.integer(d$group == 3)
d$time1 <- as.integer(d$time == 1)
d$time2 <- as.integer(d$time == 2)
d$time3 <- as.integer(d$time == 3)
# d$group1 <- as.integer(d$group == 1)
# d$group2 <- as.integer(d$group == 2)
# d$group3 <- as.integer(d$group == 3)

#### apply VGAM to renal transplant data ####
# results <- lapply(1:length(groups), function(group){
#   prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
#     if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
#     dsub <- d[d$prot == prot & d$group == group,]
#     fit <- try(vglm(SurvS4(count, cens, type = "left") ~ 1 + time1 + time2 + time3, cens.poisson,
#                     data = dsub, trace = F, link = "loglink", maxit = 50), silent = T)
#     if(class(fit) == "try-error") return(NA)
#     fits <- summary(fit)
#     pvals <- (1 - pnorm(fits@coef3[,3] * sign(fits@coef3[,3]))) * 2
#   }, mc.cores = 8))
# })
# hist(do.call(rbind, results))



d$log_count <- log(d$count)
# results <- lapply(1:length(groups), function(group){
#   prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
#     if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
#     dsub <- d[d$prot == prot & d$group == group,]
#     # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
#     fit <- try(censReg::censReg(log_count ~ 1 + time1 + time2 + time3, left = log(L), right = Inf, data = dsub), silent = T)
#     if(class(fit) == "try-error") return(NA)
#     fits <- summary(fit)
#     # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
#     pvals <- coef(fits)[,4]
#   }, mc.cores = 8))
# })
# hist(do.call(rbind, results)[,2:4], breaks = 1000)
# ps <- do.call(rbind, results)[,2:4]
# min(p.adjust(ps, method = "BH"), na.rm = T)


# results <- lapply(1:length(groups), function(group){
#   prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
#     if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
#     time_out <- sapply(1:3, function(time_i){
#       dsub <- d[d$prot == prot & d$group == group & (d$time == time_i | d$time == 0),]
#       # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
#       fit <- try(censReg::censReg(log_count ~ 1 + time1, left = log(L), right = Inf, data = dsub), silent = T)
#       if(length(class(fit)) == 1 && class(fit) == "try-error") return(NA)
#       fits <- summary(fit)
#       # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
#       pvals <- coef(fits)[2,4]
#     })
#   }, mc.cores = 8))
# })
# hist(do.call(rbind, results), breaks = 1000, xlab = "p.values", main = "")
# ps <- do.call(rbind, results)
# hist(ps[,3])
# min(p.adjust(ps, method = "BH"), na.rm = T)
# min(p.adjust(ps[,2], method = "BH"), na.rm = T)
# 
# 
# results <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
#   if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
#   dsub <- d[d$prot == prot,]
#   # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
#   fit <- try(censReg::censReg(log_count ~ 1 + group2 + group3 + time1 + time2 + time3 + time1*group2 + time2*group2 + time3*group2 + time1*group3 + time2*group3 + time3*group3 , 
#                               left = log(L), right = Inf, data = dsub), silent = T)
#   if(length(class(fit)) == 1 && class(fit) == "try-error") return(NA)
#   fits <- summary(fit)
#   # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
#   pvals <- coef(fits)[grep(":", rownames(coef(fits))),4]
# }, mc.cores = 8))
# 
# hist(results, breaks = 1000, xlab = "p.values", main = "")
# ps <- do.call(rbind, results)
# hist(ps[,3])
# min(p.adjust(ps, method = "BH"), na.rm = T)
# min(p.adjust(ps[,2], method = "BH"), na.rm = T)
# 
# #### try to invert the problem -- elastic net to predict rejection? ####
# library(caret)
# library(glmnet)
# 
# #convert from long to wide
# d$group_name <- as.factor(groups[d$group])
# d$group_name <- as.factor(groups[d$group])
# d$prot_name <- as.factor(prots[d$prot])
# d$prot_id <- as.factor(d$prot)
# d$group_patient <- as.factor(paste0(x$Group, ".", x$Patient))
# dsub <- d[d$time==0,c("prot_id", "group_name", "count", "group_patient")]
# dw <- reshape(dsub, timevar = "prot_id", idvar = "group_patient", direction= "wide", v.names = "count")
# dw <- dw[,-match("group_patient", colnames(dw))]
# y <- dw$group_name
# x <- as.matrix(dw[,-1])
# 
# traincv <- train(group_name ~ ., data = dsub, method = "multinom", trControl = trainControl(method = "cv", number = 5), trace = T)
# fit <- glmnet(x = x, dw$group_name, family = "multinomial", type.multinomial = "ungrouped",nfolds =3)
# coef(traincv)

#### ok, what about a Bayesian DE model ####
dsub <- d[d$time %in% c(0,2),]
# dsub <- d
dat <- list(n = nrow(dsub),
          n_prot = length(unique(dsub$prot)),
          n_group = length(unique(dsub$group)),
          n_time = length(unique(dsub$time)),
          censor_threshold = L,
          log_count = log(dsub$count),
          uncens = dsub$cens,
          group = dsub$group,
          time = match(dsub$time, sort(unique(dsub$time))),
          prot = dsub$prot
)

base = "basic_2timepoint_censored_poisson"
# STAN model
#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

stan_program <- '
data {
    int<lower=1> n;
    int<lower=1> n_prot;
    int<lower=1> n_group;
    int<lower=1> n_time;
    int<lower=0> censor_threshold;
    //int<lower=0> count[n];
    real log_count[n];
    int<lower=0, upper=1> uncens[n];
    int<lower=1, upper=n_group> group[n];
    int<lower=1, upper=n_time> time[n];
    int<lower=1, upper=n_prot> prot[n];
}
transformed data {
    //real log_count[n] = log(count);
    real log_censor_threshold = log(censor_threshold);
}
parameters {
    real B_raw[n_time, n_group, n_prot];
    real<lower=0> B_sd[n_time, n_group];
    real B_mean[n_time, n_group];
    real<lower=0> sigma;
}
transformed parameters {
    real B[n];
    for(i in 1:n){
      real B_temp = 0;
      for(j in 1:time[i]){
        if(group[i] == 1){
          B_temp = B_temp + B_mean[time[i], group[i]] + B_raw[time[i], group[i], prot[i]] * B_sd[time[i], group[i]];
        } else {
          B_temp = B_temp + B_mean[time[i], 1] + B_mean[time[i], group[i]] + 
                   B_raw[time[i], 1, prot[i]] * B_sd[time[i], 1] +
                   B_raw[time[i], group[i], prot[i]] * B_sd[time[i], group[i]];
        }
      }
      B[i] = B_temp;
    }
}
model {
    //priors
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,5);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,5);
    
    //likelihood
    log_count ~ normal(B, sigma);
    
    //normal model for noncensored data
    //count_0 ~ poisson(exp(t0_mean + t0_prot[prot_comp_0] + t0_indiv[indiv_comp_0]));
    //count_1 ~ poisson(exp(t1_mean + t1_prot[prot_comp_1] + t1_indiv[indiv_comp_1]));
    
    //censored obs
    //target += poisson_lcdf(censor_threshold | exp(t0_mean + t0_prot[prot_cens_0] + t0_indiv[indiv_cens_0]));
    //target += poisson_lcdf(censor_threshold | exp(t1_mean + t1_prot[prot_cens_1] + t1_indiv[indiv_cens_1]));
}
'

if(!exists("curr_stan_program") || stan_program != curr_stan_program){
  curr_stan_program <- stan_program
  f <- write_stan_file(stan_program)
}
mod <- cmdstan_model(f)

#fit model
write_stan_file(stan_program, dir = "~/Desktop/", basename = paste0(base))
write_stan_json(dat, paste0("~/Desktop/", paste0(base, ".json")))
fit_model <- T
if(fit_model){
  out <- mod$sample(chains = 4, iter_sampling = 1E2, iter_warmup = 1E2, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 1, init = 0.1, max_treedepth = 20, thin = 2)
  summ <- out$summary()
  print(summ[order(summ$ess_bulk),])
  print(summ[order(summ$rhat, decreasing = T),])
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
}

samps <- data.frame(as_draws_df(out$draws()))

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

#combine extracts
xsub <- x[,c("proteinNames", "gene", "intensity", "Patient", "Group", "TimePoint", "Extract", "experiment")]
xsplit <- split(xsub, interaction(x$proteinNames, x$Patient, x$TimePoint, x$Group))
nxsplit <- sapply(xsplit, nrow)
table(nxsplit)
repeat_experiments <- as.integer(which(nxsplit == 3))
xsplit[[repeat_experiments[1]]]

xcomb <- parallel::mclapply(xsplit, function(xs){
  if(nrow(xs) == 0){
    return(integer(0))
  } else if(nrow(xs) == 3){
    nxs <- xs[!(xs$experiment %in% gsub("r", "", xs$experiment[grepl("r", xs$experiment)])),]
  } else{
    nxs <- xs
  }
  return(cbind(nxs[1,c("proteinNames", "gene", "Patient", "Group", "TimePoint")], intensity = sum(nxs$intensity)))
}, mc.cores = 12)
xcomb <- do.call(rbind, xcomb)
x <- xcomb



### test out this VGAM package ####
library(VGAM)

N <- 10
sex <- sample(0:1, N, T)
cdata <- data.frame(y = rpois(N, exp(1.25 + sex * 0.25)), sex = sex)
L <- round(quantile(cdata$y, 0.25))
cdata <- transform(cdata,
                   cY = pmax(L, y),
                   lcensored = y <  L)  # Note y < L, not cY == L or y <= L
cdata <- transform(cdata, status = ifelse(lcensored, 0, 1))
with(cdata, table(cY))
with(cdata, table(lcensored))
with(cdata, table(print(SurvS4(cY, status, type = "left"))))  # Check
cdata <- cbind(cdata)
fit <- vglm(SurvS4(cY, status, type = "left") ~ 1 + sex, cens.poisson,
            data = cdata, trace = TRUE)
coef(fit, matrix = TRUE)
fit@misc$p
fits <- summary(fit)
fits
fits@coef3[,4]
(1 - pnorm(fits@coef3[,3] * sign(fits@coef3[,3]))) * 2

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
d$days <- x$delta
d$patient <- x$Patient
extracts <- unique(x$Extract)
d$extract <- match(x$Extract, extracts)
# d$group1 <- as.integer(d$group == 1)
# d$group2 <- as.integer(d$group == 2)
# d$group3 <- as.integer(d$group == 3)

#### apply VGAM to renal transplant data ####
try_freq_tests <- T
if(try_freq_tests){
  results_1 <- lapply(1:length(groups), function(group){
    prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
      if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
      dsub <- d[d$prot == prot & d$group == group,]
      fit <- try(vglm(SurvS4(count, cens, type = "left") ~ 1 + time1 + time2 + time3, cens.poisson,
                      data = dsub, trace = F, link = "loglink", maxit = 50), silent = T)
      if(class(fit) == "try-error") return(NA)
      fits <- summary(fit)
      pvals <- (1 - pnorm(fits@coef3[,3] * sign(fits@coef3[,3]))) * 2
    }, mc.cores = 8))
  })
  ps_1 <- do.call(rbind, results_1)
  hist(ps_1)
  min(p.adjust(ps_1, method = "BH"), na.rm = T)
  
  d$log_count <- log(d$count)
  results_2 <- lapply(1:length(groups), function(group){
    prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
      if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
      dsub <- d[d$prot == prot & d$group == group,]
      # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
      fit <- try(censReg::censReg(log_count ~ 1 + time1 + time2 + time3, left = log(L), right = Inf, data = dsub), silent = T)
      if(class(fit) == "try-error") return(NA)
      fits <- summary(fit)
      # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
      pvals <- coef(fits)[,4]
    }, mc.cores = 8))
  })
  hist(do.call(rbind, results_2)[,2:4], breaks = 1000)
  ps_2 <- do.call(rbind, results_2)[,2:4]
  min(p.adjust(ps_2, method = "BH"), na.rm = T)
  
  results_3 <- lapply(1:length(groups), function(group){
    prot_out <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
      if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
      time_out <- sapply(1:3, function(time_i){
        dsub <- d[d$prot == prot & d$group == group & (d$time == time_i | d$time == 0),]
        # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
        fit <- try(censReg::censReg(log_count ~ 1 + time1, left = log(L), right = Inf, data = dsub), silent = T)
        if(length(class(fit)) == 1 && class(fit) == "try-error") return(NA)
        fits <- summary(fit)
        # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
        pvals <- coef(fits)[2,4]
      })
    }, mc.cores = 8))
  })
  hist(do.call(rbind, results_3), breaks = 1000, xlab = "p.values", main = "")
  ps_3 <- do.call(rbind, results_3)
  hist(unlist(ps_3[,3]))
  min(p.adjust(ps_3, method = "BH"), na.rm = T)
  min(p.adjust(unlist(ps_3[,1]), method = "BH"), na.rm = T)
  
  results_4 <- do.call(rbind, parallel::mclapply(1:length(prots), function(prot){
    if(prot %% 100 == 0){system(sprintf('echo "%s "', paste0(prot, collapse="")))}
    dsub <- d[d$prot == prot,]
    # fit <- vglm(log_count ~ 1 + time1 + time2 + time3, tobit(Lower = log(L)), data = dsub)
    fit <- try(censReg::censReg(log_count ~ 1 + group2 + group3 + time1 + time2 + time3 + time1*group2 + time2*group2 + time3*group2 + time1*group3 + time2*group3 + time3*group3 ,
                                left = log(L), right = Inf, data = dsub), silent = T)
    if(length(class(fit)) == 1 && class(fit) == "try-error") return(NA)
    fits <- summary(fit)
    # pvals <- (1 - pnorm(@coef3[,3] * sign(fits@coef3[,3]))) * 2
    pvals <- coef(fits)[grep(":", rownames(coef(fits))),4]
  }, mc.cores = 8))
  
  hist(results_4, breaks = 1000, xlab = "p.values", main = "")
  ps_4 <- results_4
  hist(ps_4[,3])
  min(p.adjust(ps_4, method = "BH"), na.rm = T)
  min(p.adjust(ps_4[,2], method = "BH"), na.rm = T)
}

# #### try to invert the problem -- elastic net to predict rejection? ####
try_inversion <- F
if(try_inversion){
  library(caret)
  library(glmnet)
  
  #convert from long to wide
  d$group_name <- as.factor(groups[d$group])
  d$group_name <- as.factor(groups[d$group])
  d$prot_name <- as.factor(prots[d$prot])
  d$prot_id <- as.factor(d$prot)
  d$group_patient <- as.factor(paste0(x$Group, ".", x$Patient))
  dsub <- d[d$time==0,c("prot_id", "group_name", "count", "group_patient")]
  dw <- reshape(dsub, timevar = "prot_id", idvar = "group_patient", direction= "wide", v.names = "count")
  dw <- dw[,-match("group_patient", colnames(dw))]
  y <- dw$group_name
  x <- as.matrix(dw[,-1])
  
  traincv <- train(group_name ~ ., data = dsub, method = "multinom", trControl = trainControl(method = "cv", number = 5), trace = T)
  fit <- glmnet(x = x, dw$group_name, family = "multinomial", type.multinomial = "ungrouped",nfolds =3)
  coef(traincv)
  coef(fit)
}

#do basic EDA to find var within prot & timepoint between subj vs between prot within subj & timepoint

sapply(sort(unique(d$time)), function(ti){
  mean(sapply(sort(unique(d$prot)), function(pi){
    ds <- d$log_count[d$time == ti & d$prot == pi]
    var(ds)
  }))
})

group_indices <- setNames(sort(unique(d$group)), groups)
time_indices <- setNames(sort(unique(d$time)), paste0("t", sort(unique(d$time))))
protein_indices <- setNames(sort(unique(d$prot)), prots)
patient_indices <- setNames(sort(unique(d$patient)), paste0("patient_", sort(unique(d$patient))))
extract_indices <- setNames(sort(unique(d$extract)), extracts)

#sample sds
ei <- extract_indices[1]

#within proteins
sapply(group_indices, function(gi){
  sapply(time_indices, function(ti){
    mean(sapply(protein_indices, function(pi){
      ds <- d$log_count[d$time == ti & d$prot == pi & d$group == gi & (is.null(d$extract) || d$extract == ei)]
      sd(ds)
    }), na.rm = T)
  })
})

#within patients
sapply(group_indices, function(gi){
  sapply(time_indices, function(ti){
    mean(sapply(patient_indices, function(pi){
      ds <- d$log_count[d$time == ti & d$patient == pi & d$group == gi & (is.null(d$extract) || d$extract == ei)]
      var(ds)
    }), na.rm = T)
  })
})

#within just timepoint and group
sapply(group_indices, function(gi){
  sapply(time_indices, function(ti){
    ds <- d$log_count[d$time == ti & d$group == gi & (is.null(d$extract) || d$extract == ei)]
    sd(ds)
  })
})

#sample sds vs t0
sapply(group_indices, function(gi){
  sapply(time_indices[-1], function(ti){
    mean(sapply(protein_indices, function(pi){
      d0 <- d[d$time == 0 & d$prot == pi & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      ds <- d[d$time == ti & d$prot == pi & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      sd(ds$log_count - d0$log_count)
    }), na.rm = T)
  })
})

sapply(group_indices, function(gi){
  sapply(time_indices[-1], function(ti){
    mean(sapply(patient_indices, function(pi){
      d0 <- d[d$time == 0 & d$patient == pi & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      ds <- d[d$time == ti & d$patient == pi & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      d0 <- d0[match(ds$prot, d0$prot),]
      var(ds$log_count - d0$log_count)
    }), na.rm = T)
  })
})

sapply(group_indices, function(gi){
  sapply(time_indices[-1], function(ti){
      d0 <- d[d$time == 0 & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      ds <- d[d$time == ti & d$group == gi & (is.null(d$extract) || d$extract == ei),]
      sd(ds$log_count - d0$log_count)
  })
})

#look at extract differences?
sapply(group_indices, function(gi){
  sapply(time_indices, function(ti){
    var(c(unlist(sapply(patient_indices, function(pi){
        ds_e1 <- d[d$time == ti & d$patient == pi & d$group == gi & d$extract == 1,]
        ds_e2 <- d[d$time == ti & d$patient == pi & d$group == gi & d$extract == 2,]
        ds_e1 <- ds_e1[match(ds_e2$prot, ds_e1$prot),]
        out <- ds_e1$log_count - ds_e2$log_count
        out
    }))))
  })
})


#### ok, what about a Bayesian DE model ####
# dsub <- d[d$time %in% c(0,2),]
dsub <- d
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

# STAN model
#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

base = "basic_stairway_censored_normal_approx"
stan_program <- '
functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}
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
    real log_censor_threshold = log(censor_threshold);
    int cens_idx[count_elem(uncens, 0)] = which_elem(uncens, 0);
    real cens_log_count[count_elem(uncens, 0)] = log_count[cens_idx];
    int uncens_idx[count_elem(uncens, 1)] = which_elem(uncens, 1);
    real uncens_log_count[count_elem(uncens, 1)] = log_count[uncens_idx];
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
      int t = time[i];  int g = group[i];  int p = prot[i];
      for(j in 1:t){
        if(group[i] == 1){
          B_temp += B_mean[t, g] + B_raw[t, g, p] * B_sd[t, g];
        } else {
          B_temp += B_mean[j, 1] + B_mean[j, g] + 
                   B_raw[j, 1, p] * B_sd[j, 1] +
                   B_raw[j, g, p] * B_sd[j, g];
        }
      }
      B[i] = B_temp;
    }
}
model {
    //priors
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,2);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,2);
    
    //uncensored obs likelihood
    uncens_log_count ~ normal(B[uncens_idx], sigma);
    
    //censored obs
    target += normal_lcdf(log_censor_threshold | B[cens_idx], sigma);
    
}
'

base = "basic_stairway_censored_normal_approx_spikeNslab"
stan_program <- '
functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}
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
    real log_censor_threshold = log(censor_threshold);
    
    int<lower=0> n_cens = count_elem(uncens, 0);
    int cens_idx[n_cens] = which_elem(uncens, 0);
    real cens_log_count[n_cens] = log_count[cens_idx];
    
    int<lower=0> n_uncens = count_elem(uncens, 1);
    int uncens_idx[n_uncens] = which_elem(uncens, 1);
    real uncens_log_count[n_uncens] = log_count[uncens_idx];
}
parameters {
    real B_raw[n_time, n_group, n_prot];
    real<lower=0> B_sd[n_time, n_group];
    real B_mean[n_time, n_group];
    real<lower=0> sigma;
    
    real<lower=0, upper=1> gamma[n_group, n_prot];
    real<lower=0, upper=1> gamma_mean[n_group];
    real<lower=0.1> gamma_conc[n_group];
}
transformed parameters {
    real B_1[n];
    for(i in 1:n){
      real B_temp_1 = 0;
      int t = time[i];  int g = group[i];  int p = prot[i];
      for(j in 1:t){
        if(group[i] == 1){
          B_temp_1 += B_mean[t, g] + B_raw[t, g, p] * B_sd[t, g];
        } else {
          B_temp_1 += B_mean[j, 1] + B_mean[j, g] + 
                   B_raw[j, 1, p] * B_sd[j, 1] +
                   B_raw[j, g, p] * B_sd[j, g];
        }
      }
      B_1[i] = B_temp_1;
    }
    
    real B_2[n];
    for(i in 1:n){
      real B_temp_2 = 0;
      int t = time[i];  int g = group[i];  int p = prot[i];
      for(j in 1:t){
        B_temp_2 += B_mean[j, 1] + B_raw[j, 1, p] * B_sd[j, 1];
      }
      B_2[i] = B_temp_2;
    }
}
model {
    //(hyper)priors for coef
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,2);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,2);
    
    //(hyper)priors for mixture param
    gamma_mean ~ beta(1, 1);
    gamma_conc ~ pareto(0.1, 1.5);
    for(j in 1:n_group){
      gamma[j,] ~ beta(gamma_mean[j] * gamma_conc[j], (1 - gamma_mean[j]) * gamma_conc[j]);
    }
    
    //uncensored obs likelihood -- need to marginalize over probability each protein has 0 deviation
    for (i in 1:n_uncens) {
         target += log_mix(gamma[group[uncens_idx[i]], prot[uncens_idx[i]]],
                      normal_lpdf(uncens_log_count[i] | B_1[uncens_idx[i]], sigma), 
                      normal_lpdf(uncens_log_count[i] | B_2[uncens_idx[i]], sigma));
    }
    
    //censored obs
    for (i in cens_idx) {
         target += log_mix(gamma[group[i], prot[i]],
                      normal_lcdf(log_censor_threshold | B_1[i], sigma), 
                      normal_lcdf(log_censor_threshold | B_2[i], sigma));
    }
}
'

dsub <- d[d$group == 3 & d$time %in% c(0,2),]
dat <- list(n = nrow(dsub),
            n_prot = length(unique(dsub$prot)),
            n_group = length(unique(dsub$group)),
            n_time = length(unique(dsub$time)),
            censor_threshold = L,
            log_count = log(dsub$count),
            uncens = dsub$cens,
            group = rep(1, nrow(dsub)),
            time = match(dsub$time, sort(unique(dsub$time))),
            prot = dsub$prot
)
base = "basic_stairway_censored_normal_approx_spikeNslab_1group"
stan_program <- '
functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}
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
    real log_censor_threshold = log(censor_threshold);
    
    int<lower=0> n_cens = count_elem(uncens, 0);
    int cens_idx[n_cens] = which_elem(uncens, 0);
    real cens_log_count[n_cens] = log_count[cens_idx];
    
    int<lower=0> n_uncens = count_elem(uncens, 1);
    int uncens_idx[n_uncens] = which_elem(uncens, 1);
    real uncens_log_count[n_uncens] = log_count[uncens_idx];
}
parameters {
    real B_raw[n_time, n_group, n_prot];
    real<lower=0> B_sd[n_time, n_group];
    real B_mean[n_time, n_group];
    real<lower=0> sigma;
    
    real<lower=0, upper=1> gamma[n_group, n_prot];
    real<lower=0, upper=1> gamma_mean[n_group];
    real<lower=0.1> gamma_conc[n_group];
}
transformed parameters {

}
model {
    //(hyper)priors for coef
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,2);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,2);
    
    //(hyper)priors for mixture param
    gamma_mean ~ beta(1, 1);
    gamma_conc ~ pareto(0.1, 1.5);
    for(j in 1:n_group){
      gamma[j,] ~ beta(gamma_mean[j] * gamma_conc[j], (1 - gamma_mean[j]) * gamma_conc[j]);
    }
    
    //combine the effects
    real B_1[n];
    for(i in 1:n){
      real B_temp_1 = 0;
      int t = time[i];  int g = group[i];  int p = prot[i];
      for(j in 1:t){
        B_temp_1 += B_mean[t, g] + B_raw[t, g, p] * B_sd[t, g];
      }
      B_1[i] = B_temp_1;
    }
    
    real B_2[n];
    for(i in 1:n){
      int t = time[i];  int g = group[i];  int p = prot[i];
      B_2[i] = B_mean[1, g] + B_raw[1, g, p] * B_sd[1, g];
    }
    
    //uncensored obs likelihood
    for (i in uncens_idx) {
         target += log_mix(gamma[group[i], prot[i]],
                      normal_lpdf(log_count[i] | B_1[i], sigma), 
                      normal_lpdf(log_count[i] | B_2[i], sigma));
    }
    
    //censored obs
    for (i in cens_idx) {
         target += log_mix(gamma[group[i], prot[i]],
                      normal_lcdf(log_censor_threshold | B_1[i], sigma), 
                      normal_lcdf(log_censor_threshold | B_2[i], sigma));
    }
}
'

#ok, let's ignore the spike and slab for now
focal_group <- 1
timepoints <- c(0,2)
dsub <- d[d$group == focal_group & d$time %in% timepoints,]
dat <- list(n = nrow(dsub),
            n_prot = length(unique(dsub$prot)),
            n_group = length(unique(dsub$group)),
            n_time = length(unique(dsub$time)),
            censor_threshold = L,
            log_count = log(dsub$count),
            uncens = dsub$cens,
            group = rep(1, nrow(dsub)),
            time = match(dsub$time, sort(unique(dsub$time))),
            prot = dsub$prot
)
base <- paste0("stairway_censored_normal_approx_group-", focal_group, "_timepoints-", paste0(timepoints, collapse = ","))
stan_program <- '
functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}
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
    real log_censor_threshold = log(censor_threshold);
    
    int<lower=0> n_cens = count_elem(uncens, 0);
    int cens_idx[n_cens] = which_elem(uncens, 0);
    real cens_log_count[n_cens] = log_count[cens_idx];
    
    int<lower=0> n_uncens = count_elem(uncens, 1);
    int uncens_idx[n_uncens] = which_elem(uncens, 1);
    real uncens_log_count[n_uncens] = log_count[uncens_idx];
}
parameters {
    real B_raw[n_time, n_group, n_prot];
    real<lower=0> B_sd[n_time, n_group];
    real B_mean[n_time, n_group];
    real<lower=0> sigma[n_time];
}
model {
    //(hyper)priors for coef
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,2);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    sigma ~ normal(0,2);
    
    //combine the effects
    real B[n];
    for(i in 1:n){
      real B_temp = 0;
      int t = time[i];  int g = group[i];  int p = prot[i];
      for(j in 1:t){
        B_temp += B_mean[t, g] + B_raw[t, g, p] * B_sd[t, g];
      }
      B[i] = B_temp;
    }
    
    //uncensored obs likelihood
    uncens_log_count ~ normal(B[uncens_idx], sigma[time[uncens_idx]]);
    
    //censored obs
    target += normal_lcdf(log_censor_threshold | B[cens_idx], sigma[time[cens_idx]]);
    
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
  out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.95, refresh = 10, init = 0.1, max_treedepth = 15, thin = 2)
  summ <- out$summary()
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
  save(summ, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.diag")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.diag")))
}
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])

#examine posterior output
samps <- data.frame(as_draws_df(out$draws()))
# pairs(samps[,c("B_sd.2.1.", "gamma_conc.1.", "gamma_mean.1.")], pch = 19, col = c(rep(adjustcolor(2, 0.5), 1000),rep(adjustcolor(1, 0.5), 3000)))
# samps <- samps[1001:4000,]
subset_samps <- function(include = "", exclude = NULL, samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  if(length(exclude) != 0){
    excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
    incl_inds <- setdiff(incl_inds, excl_inds)  
  } 
  return(samps[,incl_inds])
}
prop_greater_than_0 <- function(x) mean(x>0)


#### exploring output from simplest 2 timepoint model ####

sapply(1:3, function(fg) length(unique(d[d$group == fg,"prot"]))) #same prots in each group
all_samps <- lapply(1:3, function(fg_i){
  timepoints <- c(0,2)
  base <- paste0("stairway_censored_normal_approx_group-", fg_i, "_timepoints-", paste0(timepoints, collapse = ","))
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
  samps <- data.frame(as_draws_df(out$draws()))
  samps
})
names(all_samps) <- groups

B_samps <- lapply(all_samps, function(samps){
  B_raw <- subset_samps("B_raw", samps = samps)
  B_sd <- subset_samps("B_sd", samps = samps)
  B_mean <- subset_samps("B_mean", samps = samps)
  
  B_raw_foc <- B_raw[,grep("B_raw.2.", colnames(B_raw))]
  B_nat <- crossprod(diag(B_sd$B_sd.2.1.), as.matrix(B_raw_foc))
  B <- B_nat + B_mean$B_mean.2.1.
  B
})

hist(apply(B_samps$ACR, 2, prop_greater_than_0))
hist(apply(B_samps$Normal, 2, prop_greater_than_0))
hist(apply(B_samps$CAMBR, 2, prop_greater_than_0))

hist(all_samps$Normal$B_mean.2.1.)
hist(all_samps$ACR$B_mean.2.1.)
hist(all_samps$CAMBR$B_mean.2.1.)

#### exploring output from spike and slab model ####

B_raw <- subset_samps("B_raw", samps = samps)
B_sd <- subset_samps("B_sd", samps = samps)
B_mean <- subset_samps("B_mean", samps = samps)
# gamma <- subset_samps("gamma", exclude = c("mean", "conc"), samps = samps)

hist(B_sd[,2])

vec2array <- function(samples, sep = "\\."){
  inds <- apply(do.call(rbind, strsplit(names(samples), sep))[,-1], 2, as.numeric)
  samples_array <- array(NA, c(sapply(apply(inds, 2, unique), length), nrow(samples)))
  for(i in 1:nrow(inds)){
    if(i %% 100 == 0){cat(paste0(i, " "))}
    expr <- paste0(
      "samples_array[", 
      paste0("inds[i,", 1:ncol(inds), "]", collapse = ", "),
      ",] <- samples[,i]"
    )
    eval(parse(text = expr))
  }
  samples_array
}

B_mean_array <- vec2array(B_mean)
B_raw_array <- vec2array(B_raw)
B_sd_array <- vec2array(B_sd)
# gamma_array <- vec2array(gamma)
# hist(gamma_array[1,,]) #uniform, as expected?
# hist(apply(gamma_array[1,,], 1, mean))
# curve(dbeta(x = x, 
#             shape1 = mean(samps$gamma_conc.1. * samps$gamma_mean.1.), 
#             shape2 = mean(samps$gamma_conc.1. * (1-samps$gamma_mean.1.))), 
#       0, 1)
# 
# hist(gamma_array[2,,]) 
# hist(gamma_array[3,,]) 

B_relative_array <- B_raw_array
B_array <- B_raw_array
dims <- dim(B_array) #time, group, prot
B_above_0 <- array(dim = dims[1:3])
B_above_0_relative <- array(dim = dims[1:3])
for(t in 1:dims[1]){
  for(g in 1:dims[2]){
    for(p in 1:dims[3]){
      B_temp <- rep(0, dims[length(dims)])
      B_temp_relative <- rep(0, dims[length(dims)])
      for(tc in 1:t){
         B_temp <- B_temp + B_raw_array[tc,g,p,] * B_sd_array[tc,g,] + B_mean_array[tc,g,]
         B_temp_relative <- B_temp_relative + B_raw_array[tc,g,p,] * B_sd_array[tc,g,]
      }
      B_array[t,g,p,] <- B_temp
      B_relative_array[t,g,p,] <- B_temp_relative
      B_above_0[t,g,p] <- mean(B_temp > 0)
      B_above_0_relative[t,g,p] <- mean(B_temp_relative > 0)
    }
  }
}


#### exploring output from original model ####

groups
group_cols <- c("black", "blue", "green")

par(mfrow=c(2,2))
hist(B_above_0_relative[1,-1,], breaks = 0:20/20, main = "time-1 protein-specific deviations from\naverage deviation from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0_relative[1,2,], breaks = 0:20/20, add = T, col = group_cols[2])
hist(B_above_0_relative[2,-1,], breaks = 0:20/20, main = "time-2 protein-specific deviations from\naverage deviation from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0_relative[2,2,], breaks = 0:20/20, add = T, col = group_cols[2])
hist(B_above_0_relative[3,-1,], breaks = 0:20/20, main = "time-3 protein-specific deviations from\naverage deviation from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0_relative[3,2,], breaks = 0:20/20, add = T, col = group_cols[2])
hist(B_above_0_relative[4,-1,], breaks = 0:20/20, main = "time-4 protein-specific deviations from\naverage deviation from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0_relative[4,2,], breaks = 0:20/20, add = T, col = group_cols[2])

# plot(B_above_0_relative[1,-1,], B_above_0_relative[4,-1,])

par(mfrow=c(2,2))
hist(B_above_0[1,-1,], breaks = 0:20/20, main = "time-1 protein-specific deviations from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0[1,2,], breaks = 0:20/20, add = T, col = group_cols[2])
hist(B_above_0[2,-1,], breaks = 0:20/20, main = "time-2 protein-specific deviations from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0[2,2,], breaks = 0:20/20, add = T, col = group_cols[2])
legend("topright", col = c(1,"blue","green"), legend = groups, pch = 15)
hist(B_above_0[3,-1,], breaks = 0:20/20, main = "time-3 protein-specific deviations from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0[3,2,], breaks = 0:20/20, add = T, col = group_cols[2])
hist(B_above_0[4,-1,], breaks = 0:20/20, main = "time-4 protein-specific deviations from \'Normal\' trajectory",
     xlab = "probability effect is above 0", ylab = "# proteins", col = group_cols[3])
hist(B_above_0[4,2,], breaks = 0:20/20, add = T, col = group_cols[2])


par(mfrow=c(1,1))
yvals <- do.call(rbind, lapply(1:3, function(g) cumsum(sapply(1:4, function(t) mean(B_mean_array[t,g,])))))
yvals[2,] <- yvals[2,] + yvals[1,]
yvals[3,] <- yvals[3,] + yvals[1,]
plot(1:4, yvals[1,], type = "l", ylim = range(yvals), xaxt = "n",
     xlab = "timepoint", ylab = latex2exp::TeX("$log_e(posterior\\_mean)"), lwd = 2)
axis(1, 1:4, 1:4)
lines(yvals[2,], col = "blue", lwd = 2)
lines(yvals[3,], col = "green", lwd = 2)
legend("topright", col = c(1,"blue","green"), legend = groups, lwd = 2)

groups
plot(sapply(0:3, function(t) mean(d$log_count[d$group == 3 & d$time == t & d$cens == 1])), type = "l")


#### let's try a version that incorporates patient effects ####
# dsub <- d[d$time %in% c(0,2),]
dsub <- d
dat <- list(n = nrow(dsub),
            n_prot = length(unique(dsub$prot)),
            n_group = length(unique(dsub$group)),
            n_time = length(unique(dsub$time)),
            n_patient = length(unique(dsub$patient)),
            censor_threshold = L,
            log_count = log(dsub$count),
            uncens = dsub$cens,
            group = dsub$group,
            time = match(dsub$time, sort(unique(dsub$time))),
            prot = dsub$prot,
            patient = dsub$patient
)

base = "basic_stairway_censored_normal_approx_patientRE"
# STAN model
#load libraries
library(MASS)
library(mvtnorm)
library(cmdstanr)
library(posterior)

stan_program <- '
functions {
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }
  
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
}
data {
    int<lower=1> n;
    int<lower=1> n_prot;
    int<lower=1> n_group;
    int<lower=1> n_time;
    int<lower=1> n_patient;
    int<lower=0> censor_threshold;
    //int<lower=0> count[n];
    real log_count[n];
    int<lower=0, upper=1> uncens[n];
    int<lower=1, upper=n_group> group[n];
    int<lower=1, upper=n_time> time[n];
    int<lower=1, upper=n_prot> prot[n];
    int<lower=1, upper=n_patient> patient[n];
}
transformed data {
    real log_censor_threshold = log(censor_threshold);
    int cens_idx[count_elem(uncens, 0)] = which_elem(uncens, 0);
    real cens_log_count[count_elem(uncens, 0)] = log_count[cens_idx];
    int uncens_idx[count_elem(uncens, 1)] = which_elem(uncens, 1);
    real uncens_log_count[count_elem(uncens, 1)] = log_count[uncens_idx];
}
parameters {
    real B_raw[n_time, n_group, n_prot];
    real<lower=0> B_sd[n_time, n_group];
    real B_mean[n_time, n_group];

    real B_pat_raw[n_patient, n_prot];
    real<lower=0> B_pat_sd[n_patient];
    real B_pat_mean[n_patient];
    
    real<lower=0> sigma;
}
transformed parameters {
    real B[n];
    for(i in 1:n){
      real B_temp = 0;
      int t = time[i];  int g = group[i];  int p = prot[i]; int pat = patient[i];
      B_temp += B_pat_mean[pat] + B_pat_raw[pat, p] * B_pat_sd[pat];
      for(j in 1:t){
        if(group[i] == 1){
          B_temp += B_mean[t, g] + B_raw[t, g, p] * B_sd[t, g];
        } else {
          B_temp += B_mean[j, 1] + B_mean[j, g] + 
                   B_raw[j, 1, p] * B_sd[j, 1] +
                   B_raw[j, g, p] * B_sd[j, g];
        }
      }
      B[i] = B_temp;
    }
}
model {
    //priors
    for(i in 1:n_time){
      B_mean[i,] ~ normal(0,10);
      B_sd[i,] ~ normal(0,2);
      for(j in 1:n_group){
        B_raw[i,j,] ~ std_normal();  
      }
    }
    
    for(pat in 1:n_patient){
      B_pat_raw[pat,] ~ std_normal();
    }
    B_pat_sd ~ normal(0,2);
    B_pat_mean ~ normal(0,5);
    
    sigma ~ normal(0,2);
    
    //uncensored obs likelihood
    uncens_log_count ~ normal(B[uncens_idx], sigma);
    
    //censored obs
    target += normal_lcdf(log_censor_threshold | B[cens_idx], sigma);
    
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
  out <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 1E3, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.85, refresh = 10, init = 0.1, max_treedepth = 15, thin = 1)
  summ <- out$summary()
  print(summ[order(summ$ess_bulk),])
  print(summ[order(summ$rhat, decreasing = T),])
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
}




#### starting over with the model, building to a spike and slab with within-patient effects, combined extracts ####
dsub <- d[d$group == 3 & d$time %in% c(0,2),]

#populate missing obs as censored? or remove patients without both timepoints?
missing_data_patients <- as.integer(do.call(rbind, strsplit(names(which(sapply(
  split(dsub$count, interaction(dsub$patient, dsub$time)), length) == 0)), "\\."))[,1])
dsub <- dsub[!(dsub$patient %in% missing_data_patients),]

dat <- list(n = nrow(dsub),
            n_prot = length(unique(dsub$prot)),
            n_group = length(unique(dsub$group)),
            n_time = length(unique(dsub$time)),
            n_pat = length(unique(dsub$patient)),
            censor_threshold = L,
            log_count = log(dsub$count),
            uncens = dsub$cens,
            group = rep(1, nrow(dsub)),
            time = match(dsub$time, sort(unique(dsub$time))),
            prot = dsub$prot,
            pat = match(dsub$patient, sort(unique(dsub$patient)))
)


base <- paste0("stairway_do-over_", focal_group, "_timepoints-", paste0(timepoints, collapse = ","))
stan_model <- "~/scripts/minor_scripts/postdoc/renal_transplant/simple_stairstep_difference_model.stan"
mod <- cmdstan_model(stan_model)

#fit model
write_stan_json(dat, paste0("~/Desktop/", paste0(base, ".json")))
fit_model <- T
if(fit_model){
  out <- mod$sample(chains = 4, iter_sampling = 2E3, iter_warmup = 2E3, data = dat, parallel_chains = 4, 
                    adapt_delta = 0.95, refresh = 10, init = 0.1, max_treedepth = 15, thin = 2)
  summ <- out$summary()
  save(out, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.fit")))
  save(summ, file = paste0("~/Desktop/", paste0(base, ".cmdStanR.diag")))
} else {
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.fit")))
  load(paste0("~/Desktop/", paste0(base,".cmdStanR.diag")))
}
print(summ[order(summ$ess_bulk),])
print(summ[order(summ$rhat, decreasing = T),])


#####

#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0259331
n = 1E4
Illness <- rnorm(n)
Numeracy <- rnorm(n)
Illness_Income <- -1
Numeracy_Income <- 2
Income <- Illness * Illness_Income + Numeracy * Numeracy_Income + rnorm(n)
Illness_Satisfaction <- -2
Numeracy_Satisfaction <- 0
Income_Satisfaction <- 2
Income_x_Numeracy_Satisfaction <- 0
Satisfaction <- Illness * Illness_Satisfaction + 
  Numeracy * Numeracy_Satisfaction + 
  Income * Income_Satisfaction + 
  Numeracy * Income * Income_x_Numeracy_Satisfaction + 
  rnorm(n)
summary(lm(Satisfaction ~ 0 + Income + Income:Numeracy))

#####

library(dagitty)
library(ggplot2)
library(ggdag)
dag <- dagitty(x = "dag{Satisfaction <- Illness -> Income
                      Numeracy -> Income
                      Income -> Satisfaction}")
adjustmentSets(dag, exposure = "USD", outcome = "JOY")

tidy_dag <- tidy_dagitty(dag)
is_Illnesser(tidy_dag, "ILL", "USD", "JOY")

ggdag(dag, layout = "circle", 
      text_size = 6, 
      node_size = 20, 
      text_col = 1,
      node = F) + theme_dag(plot.margin = unit(c(1,1,1,1), "cm")) +
      coord_cartesian(clip = "off")
      
par("usr")
text(x = par("usr")[1], y = par("usr")[3], labels = "+", col = "red", cex = 2)
text(x = 100, y = 20, labels = "+")



#####
par(mar = c(4,4,4,4)*1.2)
intercept = 5.27
numeracy_satisfaction <- 0.02
log10income_satisfaction <- 2.67
numeracyXlog10income_satisfaction <- 0.33
numeracy <- 0:8
income <- 1:1500*100
satisfaction <- sapply(numeracy, function(numeracy_score) intercept +
  numeracy_score * numeracy_satisfaction + 
  log10(income) * log10income_satisfaction + 
  numeracy_score * log10(income) * numeracyXlog10income_satisfaction)
plot(income / 1E3, satisfaction[,1], type = "l", ylim = range(c(satisfaction)), ylab = "Income Satisfaction", xlab = "Income ($k)")
for(i in 2:length(numeracy)){lines(income / 1E3, satisfaction[,i])}

#####

# n <- 5000
# x <- rnorm(n)
# z <- rnorm(n) + x
# y <- x * 0.5 + 7 + rnorm(n) + x^2 * 0.7 - exp(x) * 0.2 - z^3 * 0.1
# plot(x,y)
# mean(x)
# mean(y)
# predict(lm(y~x), list(x = mean(x)))


intercept <- 5.27
numeracy_satisfaction <- 0.02
log10income_satisfaction <- 2.67
numeracyXlog10income_satisfaction <- 0.33
numeracy <- 0:8
income <- 1:1500*100
mean_income <- 62400
mean_satisfaction <- 5.6
mean_numeracy <- 3.6
satisfaction <- sapply(numeracy, function(numeracy_score) 
  numeracy_score * numeracy_satisfaction + 
    log10(income) * log10income_satisfaction + 
    numeracy_score * log10(income) * numeracyXlog10income_satisfaction)
satisfaction <- sapply(1:ncol(satisfaction), function(si) satisfaction[,si] - satisfaction[which.min(abs(income - mean_income)), si] + intercept)
plot(income / 1E3, satisfaction[,1], type = "l", ylim = c(0,10), xlim = c(0,150), ylab = "Income Satisfaction", xlab = "Income ($k)")
for(i in 2:length(numeracy)){lines(income / 1E3, satisfaction[,i])}

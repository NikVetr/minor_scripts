library(dagitty)
library(ggdag)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridGraphics)



#functions
condMVN <- function(means, cov, obs, inds){
  A <- cov[-inds, -inds]
  C <- cov[inds, inds]
  B <- cov[-inds, inds]
  condCov <- A - B %*% solve(C) %*% t(B)
  condMean <- as.vector(means[-inds] + B %*% solve(C) %*% (obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}


dag_txt <- 'dag{
Time [pos="0,2.5"]
Age [pos="0,2"]
U_1 [pos="-1,1.5"]
U_2 [pos="1,1.5"]
Dis_1 [pos="-1,1"]
Dis_2 [pos="1,1"]
X_1 [pos="-2,0"]
X_2 [pos="-1,0"]
X_3 [pos="0,0"]
X_4 [pos="1,0"]
X_5 [pos="2,0"]
Time -> Age
Age -> Dis_1
Age -> Dis_2
U_1 -> Dis_1
U_2 -> Dis_2
Dis_1 -> X_1
Dis_1 -> X_2
Dis_1 -> X_3
Dis_1 -> X_4
Dis_2 -> X_2
Dis_2 -> X_3
Dis_2 -> X_4
Dis_2 -> X_5
}'


dag <- dagitty::dagitty(dag_txt)
tidy_dag <- tidy_dagitty(dag)
cols <- tail(nationalparkcolors::park_palette("GeneralGrant"), 4)
p1 <- ggdag(tidy_dag) + theme_dag_blank() + 
  geom_dag_point(col = c(cols[1], rep(cols[2], 2), cols[1], rep(cols[3], 2), rep(cols[4], 5))) + 
  geom_dag_text(col = "white", size = 3.5)
p2 <- ~{
  par(mar = c(0,0,0,0)); 
  plot(NULL, xaxt = "n", xlab = "", yaxt = "n", ylab = "",  bty="n", xlim = c(0,1), ylim = c(0,1))
  points(x = rep(-1, 4), y = 6:9/20 + 0.5, xpd = NA, pch = 19, cex = 4, col = cols[c(4,2,3,1)])
  text(x = -0.9, y = 6:9/20 + 0.5, labels = rev(c("Temporal", "Unobserved Causes", "Diseases", "Analytes")), pos = 4, xpd = NA)
  }
plot_grid(p1, p2, label_size = 0, rel_widths = c(5,1))

ggdag_dseparated(tidy_dag, from = "Age", to = "Time") + theme_dag_blank()
ggdag_dconnected(tidy_dag, from = "Age", to = "Time", controlling_for = paste0("X_", 1:5)) + theme_dag_blank()


R <- round(impliedCovarianceMatrix(dag, b.default = 0.5), 4)
cov2cor(R)
R["Age", "Dis_2"] / R["Age", "Age"]
conditional_dist <- condMVN(rep(0, dim(R)[1]), R, 
                            obs = rep(0, sum(grepl("X", rownames(R)))), 
                            grep("X", rownames(R)))
cov2cor(conditional_dist$cov)
conditional_dist$means
conditional_dist$cov["Age", "Dis_2"] / conditional_dist$cov["Age", "Age"]


#alternate model

dag_txt <- 'dag{
Age [pos="-1,1"]
Org [pos="0,0"]
Dis [pos="1,1"]
U [pos="2,2"]
X_1 [pos="-1,-1"]
X_2 [pos="0,-1"]
X_3 [pos="1,-1"]

Age -> Dis [beta=0.5]
Age -> Org [beta=-0.25]
U -> Dis [beta=0.5]
Dis -> Org [beta=-0.5]
Org -> X_1
Org -> X_2
Org -> X_3
}'

labels = c("Age" = "Chronological\nAge", 
           "Org" = "Organ\nHealth",
           "Dis" = "Disease",
           "U" = "Unobserved\nCause",
           "X_1" = "Protein 1",
           "X_2" = "Protein 2",
           "X_3" = "Protein 2")

dag <- dagitty::dagitty(dag_txt)
tidy_dag <- tidy_dagitty(dag)

cols <- tail(nationalparkcolors::park_palette("GeneralGrant"), 5)
ggdag(tidy_dag) + theme_dag_blank() + 
  geom_dag_point(col = c(cols[1], rep(cols[2], 1), rep(cols[3], 1), cols[4], rep(cols[5], 3))) + 
  geom_dag_text(col = "white", size = 3.5)
ggdag_drelationship(tidy_dag, from = "Org", to = "Age") + theme_dag_blank()
ggdag_drelationship(tidy_dag, from = "Org", to = "Age", controlling_for = paste0("X_", 1:3)) + theme_dag_blank()
impliedConditionalIndependencies(dag)

R <- round(impliedCovarianceMatrix(dag, b.lower = -0.5, b.upper = 0.5), 4)
conditional_dist <- condMVN(rep(0, dim(R)[1]), R, 
                            obs = rep(0, sum(grepl("X", rownames(R)))), 
                            grep("X", rownames(R)))


cov2cor(R)["Age", "U"]
cov2cor(conditional_dist$cov)["Age", "U"]

R["Age", "Dis"] / R["Age", "Age"]
conditional_dist$cov["Age", "Dis"] / conditional_dist$cov["Age", "Age"]

d <- simulate_data(tidy_dag, N = 1E3)
X <- as.matrix(d[,grep("X_", colnames(d))])
cor(d$U, lm(d$Age ~ X)$resid)
cor(d$Dis, lm(d$Age ~ X)$resid)

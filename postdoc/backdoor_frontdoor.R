library(dagitty)
library(ggdag)
library(ggplot2)

example <- dagify(x ~ u,
                  m ~ x,
                  y ~ u + m,
                  labels = c("x" = "Smoking",
                             "y" = "Cancer",
                             "m" = "Tar",
                             "u" = "Genotype"),
                  latent = "u",
                  exposure = "x",
                  outcome = "y")
ggdag(example)


dag <- downloadGraph("dagitty.net/m331")
ggdag(dag) +
  labs(title = "We only need to measure W to estimate the effect of X on Y")


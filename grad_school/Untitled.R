d <- read.csv("data/animalsuffering_lives_perproduct.csv")
d <- d[,-c(3,4)]
d$log.lives <- log10(d$direct.lives.per.serving)
d$log.days <- log10(d$total.suffering.per.serving)
plot(d$log.days, d$log.lives)
sort((d$log.days))

#some quick code exploring Bob's question here: https://statmodeling.stat.columbia.edu/2022/10/07/mira-magic-a-probabilistic-card-trick/#comments
# Kids divide up a standard deck of 52 playing cards and lay them out end to end in any order. The cards are big enough for kids to stand on, and laying 52 of them out creates a winding path through the exhibition hall. The exhibition is just full of these nice touches-everything is highly interactive and scaled for not just one child, but groups of children.
# Next, each kid chooses a starting point from among the first five cards and stands on the card they have chosen. At this point, the science educator leading the tour tells them where they will all end up. Each child proceeds along the path of cards the path (down the deck), each time advancing by the value of the card on which they were standing (cards 1-9 have their natural values, whereas ace, 10, jack, queen, king all count as 1).
# How does the educator know where they'll all end up and what's the probability that the trick works?

out <- data.frame(t(replicate(1E2, {
  x <- c(rep(1:9, 4), rep(rep(1, 5), 4))
  y <- sample(x)
  i0 <- sample(1:5, 1)
  i <- i0
  keep_going <- T
  while(keep_going){
    if((i + y[i]) > length(x)){
      keep_going <- F
    } else {
      i <- i + y[i]
    }
  }
  c(i_0 = i0, i_n = i)
})))

par(mfrow = c(5,1))
lapply(split(out$i_n, out$i_0), hist, breaks = min(out$i_n):max(out$i_n))


n <- 1E2
x <- c(rep(1:9, 4), rep(rep(1, 4), 4))
k <- length(x)
out <- data.frame(t(replicate(n, {
  y <- sample(x)
  endvals <- sapply(1:k, function(i0){
    i <- i0
    keep_going <- T
    while(keep_going){
      if(((i + y[i]) - length(x)) > 1E-3){
        keep_going <- F
      } else {
        i <- i + y[i]
      }
    }
    i
  })
  c(starting_value_ = y[1:k], ending_index_ = endvals)
})))
ends <- out[,grep("ending_index", colnames(out))]
mean(sapply(apply(ends, 1, unique), length) == 1)

all_together <- c(1, sapply(2:length(x), function(i) mean(sapply(apply(ends[,1:i], 1, unique), length) == 1)))
plot(1:length(x), all_together, type = "l", lwd = 2, 
     xlab = "# of kids at the start", 
     ylab = paste0("frequency they're on the same card at the end over ", n, " runs"))

plot(1:(length(x)-1), diff(all_together), type = "l", lwd = 2)
plot(1:(length(x)-2), diff(diff(all_together)), type = "l", lwd = 2)

just_2_in_5 <- mean(apply(ends[,1:5], 1, function(x) length(unique(sample(x, 2))) == 1))
just_2_in_5

hist(sapply(apply(ends, 1, unique), length), breaks = 1:9, probability = T, xlab = "number of unique ending positions", main = "")

firstn <- 10
plot(
  mean(apply(ends[,1:firstn], 1, function(x) length(unique(sample(x, firstn))) == 1))^(2:firstn/firstn),
  sapply(2:firstn, function(nchild) mean(apply(ends[,1:firstn], 1, function(x) length(unique(sample(x, nchild))) == 1))),
  pch = 19, xlab = "raising to fractional power", ylab = "sampling fresh"
)
abline(0,1,col=2,lty=2)


#look at the distribution of coupling times for pairs starting in first 5 (or other params)
n_rep <- 1E5
deck <- c(rep(1:9, 4), rep(rep(1, 4), 4))
n_cards <- length(deck)
n_players <- 2
max_start_index <- 5
max_len <- sum(cumsum(sort(deck)) <= 52)
catchup <- F
stay_on_board <- F
coupling_times <- replicate(n_rep, {
  shuffled_deck <- sample(deck)
  i0s <- sample(1:max_start_index, n_players)
  ivs <- lapply(1:n_players, function(indiv){
    iv <- matrix(NA, max_len)
    round <- 1
    i <- i0s[indiv]
    keep_going <- T
    while(keep_going){
      iv[round] <- i
      prop_i <- i + shuffled_deck[i]
      if((prop_i - n_cards) > 1E-3){
        keep_going <- F
      } else {
        round <- round + 1
        i <- prop_i
      }
    }
    if(stay_on_board){
      iv[is.na(iv)] <- max(iv[!is.na(iv)])
    }
    t(iv)
  })
  if(catchup){
    together <- Reduce(intersect, ivs)
  } else {
    together <- ivs[[1]][which(apply(do.call(rbind, ivs), 2, function(x) length(unique(x))) < 1.1)]
  }
  together <- together[!is.na(together)]
  ifelse(length(together) > 0.01, min(together), n_cards+1)
})

if(catchup){
  laymat <- matrix(c(1,1,1,2,2,1,2,2,1), 3, 3)
  layout(laymat)
} else {
  laymat <- matrix(1, 5, 5)
  laymat[1:3, 2:4] <- 2
  layout(laymat, widths = c(1,rep(2,ncol(laymat)-1)))
}
par(cex = 1.1)
hist(coupling_times, breaks = 1:(n_cards+1), probability = T, xlim = c(1, n_cards), 
     col = c(rep(adjustcolor(3, 0.5), n_cards-1), adjustcolor(2, 0.5)),
     xlab = "coupling time", ylab = "mass", xaxt = "n",
     main = paste0("PMF & CDF of coupling card indices for ", n_players, " players starting on first ", max_start_index, " cards")
     )
axis(1, at=0:n_cards, labels=rep("", n_cards+1))
xlab_vals1 <- seq(1, n_cards, by = 2)
text(y = par("usr")[3] - diff(par("usr")[3:4])/20, x = xlab_vals1 - 0.4, labels= xlab_vals1, xpd = NA)


hist.out <- hist(coupling_times, breaks = 1:(n_cards+1), plot = F)
hist.out$density <- cumsum(hist.out$density)
plot(hist.out, freq = F,
     col = c(rep(adjustcolor(3, 0.5), n_cards-1), adjustcolor(2, 0.5)),
     xlab = "coupling time", ylab = "cumulative mass", xaxt = "n", main = "")
axis(1, at=0:n_cards, labels=rep("", n_cards+1))
xlab_vals2 <- seq(1, n_cards, by = 4)
text(y = par("usr")[3] - diff(par("usr")[3:4])/10, x = xlab_vals2 - 0.4, labels= xlab_vals2, xpd = NA)

legend("topleft", legend = c("coupled before game ended", "game ended before coupling"),
       pt.bg = adjustcolor(c(3, 2), 0.5), pch = 22, pt.cex = 1.5, cex = 0.75)

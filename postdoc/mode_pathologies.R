x <- 0:1000/100
dens <- dexp(x, rate = 50)


x <- c(-rev(x), x)
dens <- c(rev(dens), dens)


dens <- dens + dnorm(x = x, mean = 6, sd = 1) * 90


dens <- dens[x > -1]
x <- x[x > -1]


par(mar = c(1,3,1,1))
plot(x,dens,type = "l", xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.lab = 2)
mtext(text = "Density", side = 2, cex = 2, line = 1)
polygon(x, dens, col = adjustcolor("darkblue", 0.2))

#####

x <- 0:1000/100
dens <- dexp(x, rate = 0)


x <- c(-rev(x), x)
dens <- c(rev(dens), dens)


dens <- dens + dnorm(x = x, mean = 6, sd = 1) * 40
dens <- dens + dnorm(x = x, mean = -6, sd = 1) * 40


dens <- dens[x > -10]
x <- x[x > -10]


par(mar = c(1,3,1,1))
plot(x,dens,type = "l", xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.lab = 2, ylim = range(dens)*1.3)
mtext(text = "Density", side = 2, cex = 2, line = 1)
polygon(x, dens, col = adjustcolor("darkgreen", 0.2))

segments(x0 = 0, x1 = 0, y0 = 0, y1 = max(dens) * 1.1, lwd = 3, lty = 2)
text(0, y = max(dens)*1.1, pos =3, "0", cex = 3)


#####

dir.create(path = "~/Documents/Documents - nikolai/bathtub_animation")
file.remove(paste0("~/Documents/Documents - nikolai/bathtub_animation/frames/", 
                   list.files("~/Documents/Documents - nikolai/bathtub_animation/frames", pattern = "*.png")))
dir.create(path = "~/Documents/Documents - nikolai/bathtub_animation/frames")

x <- -5000:3000/100
# dens <- apply(t(sapply(1:9, function(i) dnorm(x = x, mean = -10 + i*2, sd = 0.4) * i^0.2)), 2, sum)
dens <- dnorm(x = x, mean = -10, sd = 10) * 0.9 + dnorm(x = x, mean = 10, sd = 1) * 0.1
dens <- dens / sum(diff(c(0, x)) * dens)

dens_levels <- seq(max(dens) / 0.85, 0, length.out = 300)

for(ni in 1:length(dens_levels)){
# for(ni in 1:1){
  
  cat(paste0(ni, " "))
  
  dens_level <- dens_levels[ni]

  png(filename = paste0("~/Documents/Documents - nikolai/bathtub_animation/frames/", 
                        paste0(rep(0, 5-nchar(((ni - 1)) + 1)), collapse = ""), ((ni - 1)) + 1,".png"),
      width = 1400, height = 800, res = 200)
  
  dens_above_water <- which(dens > dens_level)
  mass_above_water <- sum(dens[dens_above_water] * diff(c(0,x))[dens_above_water])
  intervals_contained <- cbind(dens_above_water[which(diff(c(0,dens_above_water)) != 1)],
                               c(dens_above_water[which(diff(c(dens_above_water)) != 1)], 
                                 dens_above_water[length(dens_above_water)]))
  
  par(mar = c(1,3,1,1))
  plot(x,dens,type = "l", xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.lab = 2, ylim = range(dens)*1.3)
  mtext(text = "Density", side = 2, cex = 2, line = 1)
  polygon(c(x, rev(x)), c(dens, rep(0, length(x))), col = adjustcolor("darkgreen", 0.2))
  if(nrow(intervals_contained) != 0){
    for(i in 1:nrow(intervals_contained)){
      inds <- intervals_contained[i,1]:intervals_contained[i,2]
      polygon(c(x[inds], rev(x[inds])), c(dens[inds], rep(0, length(inds))), col = adjustcolor("red", 0.5))
    }
  }
  polygon(x = c(x, rev(x)), y = c(rep(dens_level, length(x)), rep(0, length(x))), col = adjustcolor("darkblue", 0.3))
  
  perc_contained <- round(mass_above_water * 100)
  # perc_contained <- paste0(paste0(rep(0, 3-nchar(perc_contained)), collapse = ""), perc_contained)
  text(paste0(perc_contained, "% covered"), col = "darkred", font = 1,
       x = par("usr")[1] + strwidth(paste0(rep(0, 3-nchar(perc_contained)), collapse = ""), units = "user"), y = par("usr")[4] - diff(par("usr")[3:4]) / 20, pos = 4)

  dev.off()
  
}

file.remove(paste0("~/Documents/Documents - nikolai/bathtub_animation/", 
                   list.files("~/Documents/Documents - nikolai/bathtub_animation", pattern = "*.mp4")))
system(paste0("cd \'Documents/Documents - nikolai/bathtub_animation\'; ffmpeg -r ", 
              60," -f image2 -s 1500x900 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p forw_bathtub_animation.mp4"))

#reverse and append and loop
system(paste0("cd \'Documents/Documents - nikolai/bathtub_animation\'; ", 
              "ffmpeg -i forw_bathtub_animation.mp4 -vf reverse rev_bathtub_animation.mp4; ",
              "touch input.txt;",
              "echo \"file forw_bathtub_animation.mp4\nfile rev_bathtub_animation.mp4\" > input.txt;",
              "ffmpeg -f concat -i input.txt -codec copy once_bathtub_animation.mp4; ",
              "ffmpeg -stream_loop 1 -i once_bathtub_animation.mp4 -c copy bathtub_animation.mp4"))

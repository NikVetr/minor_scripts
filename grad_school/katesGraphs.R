df <- read.csv(file = "C:\\Users\\Nikolai\\Downloads\\countF.csv", header = TRUE)

d <- as.matrix(df)
for (i in 1:length(d[1,])){
  for(j in 1:length(d[,1])){
    if(d[j,i] > 1 & !is.na(d[j,i])){
      d[j,i] <- d[j,i] - 2
    }
  }
}

par(mfrow = c(3,4))
for (i in 1:length(d[1,])){
  plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", xlab = "value", main = colnames(df[i]), axes = F)
  lines(density(d[!is.na(d[,i]),i]), col = 3, lwd = 1); 
  Axis(side=2, labels=T)
  
  #  abline(v = 1, col = 2, lwd = 2)
  #  text(labels = paste0(round(sum(d[!is.na(d[,i]),i] > 1)/length(d[!is.na(d[,i]),i]) * 100, 2), "%"), x = 2, y = 0.1)
  #  text(labels = paste0(round(sum(d[!is.na(d[,i]),i] < 1)/length(d[!is.na(d[,i]),i]) * 100, 2), "%"), x = 0, y = 0.1)
}

par(mfrow = c(2,3))

# par(mfrow = c(1,1))
i <- 1
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[i]), ", ", colnames(df[2]), ", ", colnames(df[3])), xlim = c(0,170), 
     ylim = c(0, 0.2), axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,170,10), at = seq(0,170,10))
axis(2)
par(new=TRUE)
i<-2
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = c(0,170), ylim = c(0, 0.2), axes = FALSE, ylab = ""); par(new=TRUE)
i<-3
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 0, 1, 0.5), xlim = c(0,170), ylim = c(0, 0.2), axes = F, ylab = ""); 
legend(legend = c(colnames(df[i]), colnames(df[2]), colnames(df[3])), x = 132, y = 0.2, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), lty=1, cex=0.8)


# par(mfrow = c(1,1))
i <- 5
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[5]), ", ", colnames(df[6])), xlim = c(0,120), 
     ylim = c(0, 0.35), axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,120,10), at = seq(0,120,10))
axis(2)
par(new=TRUE)
i<-6
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = c(0,120), ylim = c(0, 0.35), axes = FALSE, ylab = ""); par(new=TRUE)
legend(legend = c(colnames(df[5]), colnames(df[6])), x = 93, y = 0.35, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5)), lty=1, cex=0.8)


# par(mfrow = c(1,1))
i <- 8
xlim <- c(0,100)
ylim <- c(0,0.4)
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[i]), ", ", colnames(df[9]), ", ", colnames(df[10])), xlim = xlim, 
     ylim = ylim, axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,100,10), at = seq(0,100,10))
axis(2)
par(new=TRUE)
i<-9
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = xlim, ylim = ylim, axes = FALSE, ylab = ""); par(new=TRUE)
i<-10
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 0, 1, 0.5), xlim = xlim, ylim = ylim, axes = F, ylab = ""); 
legend(legend = c(colnames(df[8]), colnames(df[9]), colnames(df[10])), x = 80, y = 0.4, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), lty=1, cex=0.8)


# par(mfrow = c(1,1))
i <- 12
xlim <- c(0,150)
ylim <- c(0,0.25)
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[12]), ", ", colnames(df[13]), ", ", colnames(df[14])), xlim = xlim, 
     ylim = ylim, axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,xlim[2],10), at = seq(0,xlim[2],10))
axis(2)
par(new=TRUE)
i<-13
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = xlim, ylim = ylim, axes = FALSE, ylab = ""); par(new=TRUE)
i<-14
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 0, 1, 0.5), xlim = xlim, ylim = ylim, axes = F, ylab = ""); 
legend(legend = c(colnames(df[12]), colnames(df[13]), colnames(df[14])), x = 125, y = 0.25, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), lty=1, cex=0.8)


# par(mfrow = c(1,1))
i <- 16
xlim <- c(0,100)
ylim <- c(0,0.4)
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[16]), ", ", colnames(df[17]), ", ", colnames(df[18])), xlim = xlim, 
     ylim = ylim, axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,xlim[2],10), at = seq(0,xlim[2],10))
axis(2)
par(new=TRUE)
i<-17
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = xlim, ylim = ylim, axes = FALSE, ylab = ""); par(new=TRUE)
i<-18
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 0, 1, 0.5), xlim = xlim, ylim = ylim, axes = F, ylab = ""); 
legend(legend = c(colnames(df[16]), colnames(df[17]), colnames(df[18])), x = 80, y = 0.4, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), lty=1, cex=0.8)


# par(mfrow = c(1,1))
i <- 20
xlim <- c(0,80)
ylim <- c(0,0.3)
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), ylab = "frequency", 
     xlab = "value", main = paste0(colnames(df[20]), ", ", colnames(df[21]), ", ", colnames(df[22])), xlim = xlim, 
     ylim = ylim, axes = F, col = rgb(1, 0, 0, 0.5)); par(new=TRUE)
axis(1, labels = seq(0,xlim[2],10), at = seq(0,xlim[2],10))
axis(2)
par(new=TRUE)
i<-21
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 1, 0, 0.5), xlim = xlim, ylim = ylim, axes = FALSE, ylab = ""); par(new=TRUE)
i<-22
plot(table(d[!is.na(d[,i]),i])/length(d[!is.na(d[,i]),i]), col=rgb(0, 0, 1, 0.5), xlim = xlim, ylim = ylim, axes = F, ylab = ""); 
legend(legend = c(colnames(df[16]), colnames(df[17]), colnames(df[18])), x = 66, y = 0.3, 
       col = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), lty=1, cex=0.8)



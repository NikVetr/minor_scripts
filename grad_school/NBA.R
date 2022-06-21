#load in data
nbaPts <- read.csv("Downloads/Seasons_Stats.csv")
nbaPls <- read.csv("Downloads/player_data.csv")
nbaMon <- read.csv("Downloads/NBA_season1718_salary.csv")

#get heights and names
nba <- nbaPls[c("name", "height")]
nba$height <- sapply(1:length(nba$height), function(x) as.numeric(strsplit(as.character(nba$height[x]), "-")[[1]])[1]*12 + as.numeric(strsplit(as.character(nba$height[x]), "-")[[1]])[2])

#process modern data
nbaMon$height <- sapply(1:length(nbaMon$season17_18), function(x) nba$height[as.character(nba$name) == as.character(nbaMon$Player[x])])
nbaMon <- nbaMon[sapply(1:length(nbaMon$height), function(x) length(nbaMon$height[[x]])) == 1,]
nbaMon$height <- as.numeric(unlist(nbaMon$height))

#process historical data
nba$games <- sapply(1:length(nba$name), function(x) sum(nbaPts$G[as.character(nbaPts$Player) == as.character(nba$name[x])]))
nba$pts <- sapply(1:length(nba$name), function(x) sum(nbaPts$PTS[as.character(nbaPts$Player) == as.character(nba$name[x])]))
nba$min <- sapply(1:length(nba$name), function(x) sum(nbaPts$MP[as.character(nbaPts$Player) == as.character(nba$name[x])],na.rm = T))
nba$eff <- sapply(1:length(nba$name), function(x) mean(nbaPts$PER[as.character(nbaPts$Player) == as.character(nba$name[x])],na.rm = T))
nba$PointsPerGame <- nba$pts / nba$games
nba$PointsPerMinute <- nba$pts / nba$min
nba <- nba[nba$games != 0 & nba$min != 0,]

#create figure
dev.off()
par(mfrow = c(1,4))
plot(nba$height, nba$PointsPerMinute, col = rgb(0,0,0,0.2), ylim = c(0,1.6), xlab = "Player Height (inches)", ylab = "Points Scored Per Minute Played")
text(paste0("r = ", round(cor(nba$height, nba$PointsPerMinute, use = "com"), 4)), x = par("usr")[2]*0.975, y = par("usr")[4]*0.98)
plot(nba$height, nba$PointsPerGame, col = rgb(0,0,0,0.2), xlab = "Player Height (inches)", ylab = "Points Scored Per Game Played")
text(paste0("r = ", round(cor(nba$height, nba$PointsPerGame, use = "com"), 4)), x = par("usr")[2]*0.975, y = par("usr")[4]*0.98)
mtext("Comparing Player Performance Metrics to Player Height (after the Selection Filter of Entry into the NBA)", outer = T, cex = 1, side = 3, line = -3)
plot(nba$height, nba$eff, col = rgb(0,0,0,0.2), xlab = "Player Height (inches)", ylab = "Mean Player Efficiency Rating (whatever that is)")
text(paste0("r = ", round(cor(nba$height, nba$eff, use = "com"), 4)), x = par("usr")[2]*0.977, y = par("usr")[4]*0.975)
plot(nbaMon$height, nbaMon$season17_18/1E6, col = rgb(0,0,0,0.5), xlab = "Player Height (inches)", ylab = "Player Salary (2017-2018; Millions USD)")
text(paste0("r = ", round(cor(nbaMon$height, nbaMon$season17_18, use = "com"), 4)), x = par("usr")[2]*0.985, y = par("usr")[4]*0.98)

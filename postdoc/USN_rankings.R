d <- readLines("~/data/top500_USN.txt")
e <- grep("ENROLLMENT", d)
s <- c(-1, e + 1)
d <- lapply(1:(length(s)-1), function(i) d[(s[i]+2):s[i+1]])

sd(sapply(d, function(x) length(x))) == 0

d <- do.call(rbind,d)
d <- as.data.frame(d[,c(1,2,4,10,12)])
d[,c("country", "city")] <- do.call(rbind, strsplit(d$V2, "\\|"))
d$rank <- as.integer(gsub("#", "", d$V3))
d$score <- as.numeric(d$V4)
d$pop <- as.numeric(gsub(",", "", d$V5))
d$name <- d[,1]
d <- d[,-(1:5)]
d$country <- gsub("Russia", "Russian Federation", d$country)
d$country <- gsub("Hong Kong", "Hong Kong SAR, China", d$country)
d$country <- gsub("Egypt", "Egypt, Arab Rep.", d$country)
d$country <- gsub("Iran", "Iran, Islamic Rep.", d$country)
d$country <- gsub("Taiwan", "China", d$country)
d$country <- gsub("South Korea", "Korea, Rep.", d$country)
d$country <- gsub("Macau", "Macao SAR, China", d$country)
d$country <- gsub("Slovakia", "Slovak Republic", d$country)

countries <- unique(d$country)

mean_attendance <- sapply(setNames(countries, countries), function(x) mean(d$pop[d$country == x], na.rm = T))
d$pop[is.na(d$pop)] <- mean_attendance[match(d$country[is.na(d$pop)], names(mean_attendance))]
d <- d[!is.na(d$pop),]

gdp <- read.csv("~/data/API_NY.GDP.MKTP.CD_DS2_en_csv_v2_3011433/API_NY.GDP.MKTP.CD_DS2_en_csv_v2_3011433.csv")
gdp_pc <- read.csv("~/data/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_3011517/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_3011517.csv")
gdp$X2020[is.na(gdp$X2020)] <- gdp$X2019[is.na(gdp$X2020)]
gdp$geom_mean <- sapply(1:nrow(gdp), function(ri) exp(mean(as.numeric(log(gdp[ri,union(grep("X1", colnames(gdp)), grep("X2", colnames(gdp)))])), na.rm = T)))
gdp_pc$geom_mean <- sapply(1:nrow(gdp), function(ri) exp(mean(as.numeric(log(gdp_pc[ri,union(grep("X1", colnames(gdp_pc)), grep("X2", colnames(gdp_pc)))])), na.rm = T)))

d$gdp <- gdp$geom_mean[match(d$country, gdp$Country.Name)]
d$gdp_pc <- gdp_pc$geom_mean[match(d$country, gdp_pc$Country.Name)]
# sapply(d$country[is.na(d$gdp)], function(x) gdp$Country.Name[grep(x, gdp$Country.Name)])


if(all(colnames(gdp) == colnames(gdp_pc))){
  popsize <- gdp[,union(grep("X1", colnames(gdp_pc)), grep("X2", colnames(gdp_pc)))] / gdp_pc[,union(grep("X1", colnames(gdp_pc)), grep("X2", colnames(gdp_pc)))]
  popsize <- data.frame(country = gdp$Country.Name,
                   popsize = round(as.numeric(unlist(sapply(1:nrow(popsize), function(ri) ifelse(!all(is.na(popsize[ri,])), 
                                                                      popsize[ri,][max(which(!is.na(popsize[ri,])))], NA))))))
}

d$scoreXpop <- d$score * d$pop
percap_score <- sapply(setNames(countries, countries), function(x) sum(d$scoreXpop[d$country == x]) / sum(d$pop[d$country == x]) )
#alternatively
# percap_score <- sapply(setNames(countries, countries), function(x) sum(d$scoreXpop[d$country == x]) / popsize$popsize[popsize$country == x] )
percap_geom_gdp <- gdp_pc$geom_mean[match(names(percap_score), gdp_pc$Country.Name)]

library(ggrepel)
dat <- data.frame(per_capita_score = percap_score, per_capita_gdp = percap_geom_gdp)
dat$country <- rownames(dat)
dat <- dat[complete.cases(dat),]

p <- ggplot(dat, aes(per_capita_gdp, per_capita_score, label = country)) +
  geom_point(color = "red") + theme_bw() + labs(y = "Average Global Score per Student", x = "Geometric Mean GDP per Capita (1960 - 2020)")
p2 <- p + geom_text_repel() + 
  labs(title = paste0("Comparing Average School Quality Across Countries (r = ", round(cor(dat[,1:2])[1,2], 2), ")"))  + geom_smooth(method = "lm", se = TRUE)
p2

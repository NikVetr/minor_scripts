act_full <- read.csv("C:/Users/nikol/Downloads/atusact-0319/atusact_0319.dat")
act <- act_full[,c("TUCASEID", "TUACTIVITY_N", "TUACTDUR24")]
who <- read.csv("C:/Users/nikol/Downloads/atuswho-0319/atuswho_0319.dat")
who_codes <- read.csv("C:/Users/nikol/Downloads/TUWHO_CODES.csv")
colnames(who_codes) <- c("code", "value")
who_codes <- rbind(who_codes, (c(-1, "NA")))
roster <- read.csv("C:/Users/nikol/Downloads/atusrost-0319/atusrost_0319.dat")
roster <- roster[roster$TERRP %in% c(18,19),]

#some quick sanity checks
act <- act[order(act$TUCASEID),]
who <- who[order(who$TUCASEID),]
setdiff(act$TUCASEID, who$TUCASEID)

inds_act <- c(0, cumsum(table(act$TUCASEID))) #household inds
inds_who <- c(0, cumsum(table(who$TUCASEID))) #household inds
length(inds_act) == length(inds_who)
sapply(1:(length(inds_act) - 1), function(x) sum(act$TUACTDUR24[(inds_act[x]+1):(inds_act[x+1])])) #so are households individuals? yes I think so

#get respondent age and sex in there
act$age <- roster$TEAGE[match(act$TUCASEID, roster$TUCASEID)]
act$sex <- roster$TESEX[match(act$TUCASEID, roster$TUCASEID)]

#get who they were with in there
act <- cbind(act, matrix(0, nrow = nrow(act), ncol = length(who_codes$code)))
colnames(act) <- c(colnames(act)[1:3], c("age", "sex"), who_codes$code)
milestones <- round(seq(1, length(inds_act), length.out = 100))
for(i in 1:(length(inds_act) - 1)){
  if(i %in% milestones){
    print(round(i / length(inds_act) * 100))
  }
  subact <- act[(inds_act[i]+1):(inds_act[i+1]),]
  subwho <- who[(inds_who[i]+1):(inds_who[i+1]),]
  for(a in sort(unique(subact$TUACTIVITY_N))){
    people <- table(as.character(subwho$TUWHO_CODE[subwho$TUACTIVITY_N == a]))
    subact[subact$TUACTIVITY_N == a, names(people)] <- as.numeric(people)
  }
  subact -> act[(inds_act[i]+1):(inds_act[i+1]),]
  
}

write.table(act, file = "H:/act.txt")
write.table(who, file = "H:/who.txt")


#visualize network of who Americans tend to hang out with in terms of age
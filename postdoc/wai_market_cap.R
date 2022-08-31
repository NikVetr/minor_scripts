
mcap <- read.csv("~/data/companiesmarketcap.com - Companies ranked by Market Cap - CompaniesMarketCap.com.csv")
exploit <- read.csv("~/data/crueltyfreeinvesting_companies_that_exploit.csv")
colnames(exploit) <- exploit[1,]
exploit <- exploit[-1,]
nexploit <- read.csv("~/data/crueltyfreeinvesting_companies_that_do_not_exploit.csv")
colnames(nexploit) <- nexploit[1,]
nexploit <- nexploit[-1,]

expl_inds <- match(exploit$`Stock Symbol`, mcap$Symbol)
expl_inds <- expl_inds[!is.na(expl_inds)]
mcap_exploit <- mcap[expl_inds,]

testing_inds <- match(exploit$`Stock Symbol`[exploit$`Animal Usage` == "Animal Testing"], mcap$Symbol)
testing_inds <- testing_inds[!is.na(testing_inds)]
mcap_testing <- mcap[testing_inds,]

nexpl_inds <- match(nexploit$`Stock Symbol`, mcap$Symbol)
nexpl_inds <- nexpl_inds[!is.na(nexpl_inds)]
mcap_nexploit <- mcap[nexpl_inds,]

expl_total <- sum(mcap_exploit$marketcap)
nexpl_total <- sum(mcap_nexploit$marketcap)

nexpl_total / (nexpl_total + expl_total)

nexpl_total / sum(mcap$marketcap)

sum(mcap$marketcap) - expl_total - nexpl_total


d <- as.data.frame(readxl::read_excel("G:\\Downloads\\WIID_19Dec2018.xlsx"))
d <- d[d$year == "2020" & d$currency == "Euro", c("country", "population", "mean")]

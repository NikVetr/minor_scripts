#input produced using opendetex (https://github.com/pkubowicz/opendetex)
# sed 's/\\autocite/\\cite/g' motrpac_companion.tex > motrpac_companion_simplecite.tex 
# detex -n /Users/nikgvetr/repos/opendetex/motrpac_companion_simplecite.tex > motrpac_companion_plain.txt

library(quanteda)
library(quanteda.textstats)

txt <- readLines("~/repos/opendetex/motrpac_companion_plain.txt")
txt <- txt[nchar(txt) > 15 & grepl(" ", txt) & grepl("[.?!…]$", txt)]
corp <- corpus(paste0(txt, collapse = " "))
sentences <- tokens(corp, what = "sentence")
strings <- as.list(sentences)[[1]]
strings <- strings[stringr::str_count(strings, " ") > 3]
flesch <- data.frame(txt = strings, score = textstat_readability(corpus(strings), measure = "Flesch")$Flesch)

flesch <- flesch[order(flesch$score),]
hard <- flesch[flesch$score < 10,]
write.csv(hard, "/Volumes/2TB_External/low_readability_sentences.csv")

#edits for chatGPT suggestions
new_txt <- read.csv("/Volumes/2TB_External/low_readability_sentences.csv")
textstat_readability(corpus(new_txt$new_txt), measure = "Flesch")$Flesch
cat(paste0("\n\nold: ", new_txt$txt, "\nnew: ", new_txt$new_txt))

library(rvest)
library(tidytext)
library(tidyverse)
library(tibble)

#functions

get_thesaurus_synonyms <- function(word, levels = c(1,2,3)){
  
  word <- gsub(" ", "%20", word)

  url <- paste0("https://www.thesaurus.com/browse/", word)
  html <- read_html(url)

  synonyms <- trimws(c(
    if(1 %in% levels){html_text(html_elements(html, xpath = "//a[@class='css-1kg1yv8 eh475bn0']"))}else{character()},
    if(2 %in% levels){html_text(html_elements(html, xpath = "//a[@class='css-1gyuw4i eh475bn0']"))}else{character()},
    if(3 %in% levels){html_text(html_elements(html, xpath = "//a[@class='css-1n6g4vv eh475bn0']"))}else{character()}
  ))
  

  
  # html_text(html_elements(html, xpath = "//a[@data-linkid='nn1ov4']"))
  
  return(synonyms)

}

get_thesaurus_synonyms_recursive <- function(word, n = 2, levels = c(1,2,3)){
  
  synonyms <- get_thesaurus_synonyms(word, levels)
  
  if(n > 1){
    synonyms <- c(synonyms, sapply(synonyms, function(syn_word) get_thesaurus_synonyms_recursive(syn_word, n-1, levels)))
  }
  
  return(as.character(unlist(synonyms)))
  
}

keywords <- c("slice", "punch", "kill", "stab")
physical_violence_words <- unlist(lapply(keywords, function(word) get_thesaurus_synonyms_recursive(word, 2, 2)))
word_counts <- table(physical_violence_words)
physical_violence_words <- names(word_counts)[word_counts >= quantile(word_counts, probs = 0.9)]

dictionary <- rbind(get_sentiments("nrc"), tibble(word = physical_violence_words, sentiment = "violence"))
dictionary <- dictionary[dictionary$sentiment %in% c("fear", "anger", "violence"),]

count_sentiments <- function(tokens, dictionary){
  sentiment_options <- unique(dictionary$sentiment)
  sentiments <- tokens %>% inner_join(dictionary, by = "word") %>% count(sentiment) %>% spread(sentiment, n, fill = 0)
  sentiments_sparse <- data.frame(t(rep(0, length(sentiment_options))))
  colnames(sentiments_sparse) <- sentiment_options
  if(ncol(sentiments) > 0){
    sentiments_sparse[1,colnames(sentiments)] <- sentiments
  }
  return(sentiments_sparse)
}


setwd("~")
setwd("WtC")
story <-  data.frame(value = paste(readLines("WtC.txt", warn = F, encoding = "UTF-8")))
story_sentence_tokens <- unnest_tokens(story, sentence, value, token = "sentences")
sentence_tokens <- lapply(story_sentence_tokens$sentence, function(sentence) unnest_tokens(data.frame(value = sentence), word, value, token = "words"))

n_sentences = 1E4
sentence_sentiments_orig <-  do.call(rbind, 
                               lapply(sentence_tokens[1:n_sentences], function(sentence_tokens_indiv) unlist(count_sentiments(sentence_tokens_indiv, dictionary))))
sentence_sentiments <- as.data.frame(sentence_sentiments_orig)

#reweigh by sentence length to get sentence density
# n_words_per_sentence <- sapply(sentence_tokens[1:n_sentences], function(x) nrow(x))
# density_weigher <- diag(1/n_words_per_sentence)^0.1
# sentence_sentiments <- density_weigher %*% as.matrix(sentence_sentiments)

double_exponential_spreader <- function(x, rate = 1){
  weight_matrix <- dexp(as.matrix(dist(1:length(x))), rate = rate)
  weight_matrix <- weight_matrix * 1 / weight_matrix[1,1]
  weight_matrix %*% x 
}

sentence_sentiments <- apply(sentence_sentiments, 2, double_exponential_spreader)

sentence_sentiments_order <- apply(sentence_sentiments, 2, rank)
div_by_max <- function(x){x/max(x)}
sentence_sentiments_quantile <- apply(sentence_sentiments_order, 2, div_by_max)
sentence_sentiments_prod <- apply(sentence_sentiments_quantile[,c("fear","anger")], 1, max) * sentence_sentiments_quantile[,c("violence")]
story_sentence_tokens$sentence[1:n_sentences][order(sentence_sentiments_prod, decreasing = T)[1:20]]

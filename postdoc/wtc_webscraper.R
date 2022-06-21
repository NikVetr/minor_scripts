#dir.create("WtC")
setwd("~")
setwd("WtC")
#dir.create("comments")
library(RedditExtractoR)
library(data.table)
library(tidytext)
library(tidyverse)
source("~/double_exponential_smoother.R")


get_chapters_and_comments <- F
if(get_chapters_and_comments){
  urls <- reddit_urls(search_terms = "worth the candle, ch", subreddit = "rational", sort_by = "new", page_threshold = 4)
  urls$title <- tolower(urls$title)
  urls$chapters <- as.character(gsub(",", "", sapply(sapply(urls$title, function(title) strsplit(title, " ch ")[[1]][2]), function(title) strsplit(title, " ")[[1]][1])))
  urls <- urls[!is.na(urls$chapters),]
  urls$chapters <- gsub(pattern = "77\02378", replacement = "77", urls$chapters)
  urls$chapters <- gsub(pattern = "/", replacement = "-", urls$chapters)
  
  fwrite(urls, file = "url_info.txt")
  for(url_i in 1:nrow(urls)){
    print(urls$chapters[url_i])
    fwrite(data.table(reddit_content(urls$URL[url_i])), file = paste0("comments/WtC_", urls$chapters[url_i], ".txt"))
  }
  
  for(url_i in 1:nrow(urls)){
    print(urls$chapters[url_i])
    fread(file = paste0("comments/WtC_", urls$chapters[url_i], ".txt"))
  }
  
  comments <- rbindlist(lapply(paste0("comments/WtC_", rev(urls$chapters), ".txt"), fread), use.names = TRUE)
  comments$title <- tolower(comments$title)
  comments$chapters <- as.character(gsub(",", "", sapply(sapply(comments$title, function(title) strsplit(title, " ch ")[[1]][2]), function(title) strsplit(title, " ")[[1]][1])))
  comments$chapters <- gsub(pattern = "77\02378", replacement = "77", comments$chapters)
  comments$chapters <- gsub(pattern = "/", replacement = "-", comments$chapters)
  fwrite(comments, file = "comments.txt")
  
  # comments[grep(x = comments$user, "phylogenik"),][4,]
  story <- readLines("WtC.txt", warn = F)
  lines_with_chapter <- grep("Chapter ", story)
  lines_with_chapter <- lines_with_chapter[sapply(paste0("Chapter ", 1:246, ":"), function(chap_num) grep(chap_num, story[lines_with_chapter]))]
  # plot(diff(lines_with_chapter), type = "l")
  for(chi in 1:(length(lines_with_chapter)-1)){
    print(chi)
    writeLines(story[lines_with_chapter[chi]:(lines_with_chapter[chi+1]-1)],con = paste0("chapters/WtC_ch-", chi, ".txt"))
  }
  writeLines(story[lines_with_chapter[length(lines_with_chapter)]:length(story)],con = paste0("chapters/WtC_ch-", length(lines_with_chapter), ".txt"))
  
}

urls <- fread("url_info.txt")
comments <- fread("comments.txt")
comments$chapters <- as.numeric(sapply(strsplit(comments$chapters, "-"), function(x) x[1]))
comments <- comments[order(comments$chapters),]
comments <- comments[-grep(pattern = "not everything", comments$title),]
comments <- comments[-grep(pattern = "unicorn in the room", comments$title),]
comments <- comments[-grep(pattern = "666 - hells", comments$title),]

starting_chapters <- sort(unique(comments$chapters))
comment_text <- sapply(starting_chapters, function(chpts) paste0(comments$comment[comments$chapters == chpts], collapse = " "))

#find word uniqueness info
story <-  data.frame(value = paste(readLines("WtC.txt", warn = F)))
story_tokens <- story %>% unnest_tokens(word, value)
word_freq <- read.csv("C:\\Users\\nikol\\Downloads\\unigram_freq.csv\\unigram_freq.csv")
word_freq$freq <- word_freq$count / sum(word_freq$count)
word_tibble <- tibble(word = word_freq$word, sentiment = word_freq$word)
word_sentiments <- story_tokens %>% inner_join(word_tibble) %>% count(sentiment) %>% spread(sentiment, n, fill = 0)
word_sentiments <- word_sentiments / nrow(story_tokens)
word_sentiments <- word_sentiments / word_freq$freq[match(names(word_sentiments), word_freq$word)]
# word_sentiments <- word_sentiments / min(word_sentiments)
word_sentiments <- word_sentiments[order(word_sentiments, decreasing = T)]
word_sentiments <- data.frame(word = names(word_sentiments), weighted_freq = as.numeric(word_sentiments))
characters <- c("Juniper", "Joon", "Amaryllis", "Mary", "Fenn", "Bethel", "Solace", "Raven", "Valencia", "Val", 
                "Grakhuil", "Grak", "Arthur", "Uther", "Reimer", "Tiff", "Maddie", "Tom", "Craig", "Ropey", "Doona", "Pallida", "Lisi", "Aumann")
characters <- tolower(characters)
word_sentiments <- word_sentiments[!c(word_sentiments$word %in% characters),]


png(filename = "wtc_word_frequency.png", width = 1000, height = 12000)
ggplot(word_sentiments[1:750,], aes(x=weighted_freq, y=reorder(word, weighted_freq))) + 
  geom_bar(stat = "identity", width=0.2) + theme_bw(base_size=24) + scale_x_log10() +
  xlab("Word Frequency Relative to Web Corpus") + ylab("Word") + ggtitle("Worth the Candle Word Frequencies")
dev.off()


not_appearing <- story_tokens$word[!(story_tokens$word %in% word_freq$word)]
sort(table(not_appearing), decreasing = T)



#process comment-based data
comment_summary <- data.frame()
characters <- c("Juniper", "Joon", "Amaryllis", "Mary", "Fenn", "Bethel", "Solace", "Raven", "Valencia", "Val", 
                "Grakhuil", "Grak", "Arthur", "Uther", "Reimer", "Tiff", "Maddie", "Tom", "Craig")
character_tibble <- tibble(word = tolower(characters), sentiment = characters)
sentiment_options <- unique(get_sentiments("nrc")$sentiment)
sentiment_options_bing <- unique(get_sentiments("bing")$sentiment)
for(chi in 1:length(comment_text)){
  if(chi != length(comment_text)){chapters_included <- (starting_chapters[chi]:(starting_chapters[chi+1]-1))} else {chapters_included <-starting_chapters[chi]}
  print(chi)
  chapters <- data.frame(value = comment_text[chi])
  chapters_tokens <- chapters %>% unnest_tokens(word, value)
  # chapters_tokens %>% count(word, sort = T) %>% head
 
  #get chapter text
  actual_chapter <- data.frame(value = readLines(con = paste0("chapters/WtC_ch-", chapters_included[1], ".txt")))
  actual_chapter_tokens <- actual_chapter %>% unnest_tokens(word, value)
  if(length(chapters_included) > 1){
    for(i in 2:length(chapters_included)){
      actual_chapter <- data.frame(value = readLines(con = paste0("chapters/WtC_ch-", chapters_included[i], ".txt")))
      actual_chapter_tokens_extra <- actual_chapter %>% unnest_tokens(word, value)
      actual_chapter_tokens <- rbind(actual_chapter_tokens, actual_chapter_tokens_extra)
    }
  }
  
  #first do the chapter text
  char_sentiments <- actual_chapter_tokens %>% inner_join(character_tibble) %>% count(sentiment) %>% spread(sentiment, n, fill = 0)
  char_sentiments_sparse <- data.frame(t(rep(0, length(characters))))
  colnames(char_sentiments_sparse) <- characters
  if(ncol(char_sentiments_sparse) > 0){
    char_sentiments_sparse[1,colnames(char_sentiments)] <- char_sentiments
  }
  
  sentiments <- actual_chapter_tokens %>% inner_join(get_sentiments("nrc")) %>% count(sentiment) %>% spread(sentiment, n, fill = 0) %>% mutate(sentiment = positive - negative) 
  sentiments_sparse <- data.frame(t(rep(0, length(sentiment_options))))
  colnames(sentiments_sparse) <- sentiment_options
  if(ncol(sentiments_sparse) > 0){
    sentiments_sparse[1,colnames(sentiments)] <- sentiments
  }
  
  bing_sentiments <- actual_chapter_tokens %>% inner_join(get_sentiments("bing")) %>% count(sentiment) %>% spread(sentiment, n, fill = 0) %>% mutate(sentiment = positive - negative) 
  bing_sentiments_sparse <- data.frame(t(rep(0, length(sentiment_options_bing))))
  colnames(bing_sentiments_sparse) <- sentiment_options_bing
  if(ncol(bing_sentiments_sparse) > 0){
    bing_sentiments_sparse[1,colnames(bing_sentiments)] <- bing_sentiments
  }
  names(bing_sentiments_sparse) <- paste0("bing_", names(bing_sentiments_sparse))
  
  actual_chapter_output <- cbind(words = nrow(actual_chapter_tokens), sentiments_sparse, bing_sentiments_sparse, char_sentiments_sparse)
  names(actual_chapter_output) <- paste0("actual-text_", names(actual_chapter_output))
  
  #now do the comments
  
  char_sentiments <- chapters_tokens %>% inner_join(character_tibble) %>% count(sentiment) %>% spread(sentiment, n, fill = 0)
  char_sentiments_sparse <- data.frame(t(rep(0, length(characters))))
  colnames(char_sentiments_sparse) <- characters
  if(ncol(char_sentiments_sparse) > 0){
    char_sentiments_sparse[1,colnames(char_sentiments)] <- char_sentiments
  }
  
  sentiments <- chapters_tokens %>% inner_join(get_sentiments("nrc")) %>% count(sentiment) %>% spread(sentiment, n, fill = 0) %>% mutate(sentiment = positive - negative) 
  sentiments_sparse <- data.frame(t(rep(0, length(sentiment_options))))
  colnames(sentiments_sparse) <- sentiment_options
  if(ncol(sentiments_sparse) > 0){
    sentiments_sparse[1,colnames(sentiments)] <- sentiments
  }
  
  bing_sentiments <- chapters_tokens %>% inner_join(get_sentiments("bing")) %>% count(sentiment) %>% spread(sentiment, n, fill = 0) %>% mutate(sentiment = positive - negative) 
  bing_sentiments_sparse <- data.frame(t(rep(0, length(sentiment_options_bing))))
  colnames(bing_sentiments_sparse) <- sentiment_options_bing
  if(ncol(bing_sentiments_sparse) > 0){
    bing_sentiments_sparse[1,colnames(bing_sentiments)] <- bing_sentiments
  }
  names(bing_sentiments_sparse) <- paste0("bing_", names(bing_sentiments_sparse))
  
  comment_output <- cbind(words = nrow(chapters_tokens), sentiments_sparse, bing_sentiments_sparse, char_sentiments_sparse)
  names(comment_output) <- paste0("comments_", names(comment_output))
  
  comment_summary <- rbind(comment_summary, cbind(chapters = chi, comment_output, actual_chapter_output))
}


comment_summary[,"comments_Juniper-Joon"] <- comment_summary[,"comments_Juniper"] + comment_summary[,"comments_Joon"]
comment_summary[,"comments_Amaryllis-Mary"] <- comment_summary[,"comments_Amaryllis"] + comment_summary[,"comments_Mary"]
comment_summary[,"comments_Arthur-Uther"] <- comment_summary[,"comments_Arthur"] + comment_summary[,"comments_Uther"]
comment_summary[,"comments_Grakhuil-Grak"] <- comment_summary[,"comments_Grakhuil"] + comment_summary[,"comments_Grak"]
comment_summary[,"comments_Valencia-Val"] <- comment_summary[,"comments_Valencia"] + comment_summary[,"comments_Val"]
comment_summary <- comment_summary[,-match(c("comments_Juniper","comments_Joon","comments_Amaryllis","comments_Mary","comments_Valencia","comments_Val","comments_Grakhuil","comments_Grak","comments_Arthur","comments_Uther"), colnames(comment_summary))]

comment_summary[,"actual-text_Juniper-Joon"] <- comment_summary[,"actual-text_Juniper"] + comment_summary[,"actual-text_Joon"]
comment_summary[,"actual-text_Amaryllis-Mary"] <- comment_summary[,"actual-text_Amaryllis"] + comment_summary[,"actual-text_Mary"]
comment_summary[,"actual-text_Arthur-Uther"] <- comment_summary[,"actual-text_Arthur"] + comment_summary[,"actual-text_Uther"]
comment_summary[,"actual-text_Grakhuil-Grak"] <- comment_summary[,"actual-text_Grakhuil"] + comment_summary[,"actual-text_Grak"]
comment_summary[,"actual-text_Valencia-Val"] <- comment_summary[,"actual-text_Valencia"] + comment_summary[,"actual-text_Val"]
comment_summary <- comment_summary[,-match(c("actual-text_Juniper","actual-text_Joon","actual-text_Amaryllis","actual-text_Mary","actual-text_Valencia","actual-text_Val","actual-text_Grakhuil","actual-text_Grak","actual-text_Arthur","actual-text_Uther"), colnames(comment_summary))]

comment_inds <- setdiff(grep(colnames(comment_summary), pattern = "comments_"), grep(colnames(comment_summary), pattern = "words"))
text_inds <- setdiff(grep(colnames(comment_summary), pattern = "actual-text_"), grep(colnames(comment_summary), pattern = "words"))

comment_summary[,comment_inds] <- diag(1/comment_summary$comments_words) %*% as.matrix(comment_summary[,comment_inds])
comment_summary[,text_inds] <- diag(1/comment_summary$`actual-text_words`) %*% as.matrix(comment_summary[,text_inds])
normalize <- function(x){(x - mean(x))/sd(x)}
comment_summary[,-1] <- apply(comment_summary[,-1],2,normalize)

par(mfrow = c(5,6), mar = c(3,4,4,2))
for(ci in 2:ncol(comment_summary)){
  plot(comment_summary$chapters, comment_summary[,ci], type = "l", xlab = "Chapter", main = toupper(colnames(comment_summary)[ci]),
       ylab = "Normalized Sentiment Density")
  lines(comment_summary$chapters,
        dexp_smooth(comment_summary$chapters, comment_summary[,ci], r = 0.5),
        col = 2, lwd = 3)
}

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reordered_cormat <- reorder_cormat(cor(apply(comment_summary[,-1], 2, diff), method = "spearman"))
colours <- rep("purple", nrow(reordered_cormat))
colours[grep("actual-text", rownames(reordered_cormat))] <- "darkgreen"
melted_cormat <- reshape2::melt(reordered_cormat)
melted_cormat$value <- round(melted_cormat$value, 2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1, colour = colours),
        axis.text.y = element_text(colour = colours))+
  coord_fixed()



#process chapter-specfic sentiments
sentiment_chapter_summary <- data.frame()
for(chi in 1:246){
  print(chi)
  chapter <- data.frame(value = readLines(con = paste0("chapters/WtC_ch-", chi, ".txt")))
  chapter_tokens <- chapter %>% unnest_tokens(word, value)
  # chapter_tokens %>% count(word, sort = T) %>% head
  sentiments <- chapter_tokens %>%
    inner_join(get_sentiments("nrc")) %>% # pull out only sentiment words
    count(sentiment) %>% # count the # of positive & negative words
    spread(sentiment, n, fill = 0) %>% # made data wide rather than narrow
    mutate(sentiment = positive - negative) # # of positive words - # of negative owrds
  bing_sentiments <- chapter_tokens %>% inner_join(get_sentiments("bing")) %>% count(sentiment) %>% spread(sentiment, n, fill = 0) %>% mutate(sentiment = positive - negative) 
  names(bing_sentiments) <- paste0("bing_", names(bing_sentiments))
  sentiment_chapter_summary <- rbind(sentiment_chapter_summary, cbind(chapter = chi, words = nrow(chapter_tokens), sentiments, bing_sentiments))
}

sentiment_chapter_summary[,-(1:2)] <- diag(1/sentiment_chapter_summary$words) %*% as.matrix(sentiment_chapter_summary[,-(1:2)])

normalize <- function(x){(x - mean(x))/sd(x)}
sentiment_chapter_summary[,-(1:2)] <- apply(sentiment_chapter_summary[,-(1:2)],2,normalize)

par(mfrow = c(4,4), mar = c(3,4,4,2))
for(ci in 2:ncol(sentiment_chapter_summary)){
  plot(sentiment_chapter_summary$chapter, sentiment_chapter_summary[,ci], type = "l", xlab = "Chapter", main = toupper(colnames(sentiment_chapter_summary)[ci]),
       ylab = "Normalized Sentiment Density")
  lines(sentiment_chapter_summary$chapter,
        dexp_smooth(sentiment_chapter_summary$chapter, sentiment_chapter_summary[,ci], r = 0.1),
        col = 2, lwd = 3)
}


# plot(sentiment_chapter_summary$bing_positive, sentiment_chapter_summary$positive)
# plot(sentiment_chapter_summary$anger, sentiment_chapter_summary$disgust)

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

melted_cormat <- reshape2::melt(reorder_cormat(cor(apply(sentiment_chapter_summary[,-1], 2, diff))))
melted_cormat$value <- round(melted_cormat$value, 2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#process chapter-specfic character frequencies
characters <- c("Juniper", "Joon", "Amaryllis", "Mary", "Fenn", "Bethel", "Solace", "Raven", "Valencia", "Val", 
                "Grakhuil", "Grak", "Arthur", "Uther", "Reimer", "Tiff", "Maddie", "Tom", "Craig")
character_tibble <- tibble(word = tolower(characters), sentiment = characters)
character_chapter_summary <- data.frame()
for(chi in 1:246){
  print(chi)
  chapter <- data.frame(value = readLines(con = paste0("chapters/WtC_ch-", chi, ".txt")))
  chapter_tokens <- chapter %>% unnest_tokens(word, value)
  # chapter_tokens %>% count(word, sort = T) %>% head
  char_sentiments <- chapter_tokens %>% inner_join(character_tibble) %>% count(sentiment) %>% spread(sentiment, n, fill = 0)
  char_sentiments_sparse <- data.frame(t(rep(0, length(characters))))
  colnames(char_sentiments_sparse) <- characters
  if(ncol(char_sentiments_sparse) > 0){
    char_sentiments_sparse[1,colnames(char_sentiments)] <- char_sentiments
  }
  character_chapter_summary <- rbind(character_chapter_summary, cbind(chapter = chi, words = nrow(chapter_tokens), char_sentiments_sparse))
}
character_chapter_summary[,"Juniper-Joon"] <- character_chapter_summary[,"Juniper"] + character_chapter_summary[,"Joon"]
character_chapter_summary[,"Amaryllis-Mary"] <- character_chapter_summary[,"Amaryllis"] + character_chapter_summary[,"Mary"]
character_chapter_summary[,"Arthur-Uther"] <- character_chapter_summary[,"Arthur"] + character_chapter_summary[,"Uther"]
character_chapter_summary[,"Grakhuil-Grak"] <- character_chapter_summary[,"Grakhuil"] + character_chapter_summary[,"Grak"]
character_chapter_summary[,"Valencia-Val"] <- character_chapter_summary[,"Valencia"] + character_chapter_summary[,"Val"]
character_chapter_summary <- character_chapter_summary[,-match(c("Juniper", "Joon", "Amaryllis", "Mary", "Valencia", "Val", "Grakhuil", "Grak", "Arthur", "Uther"), colnames(character_chapter_summary))]


character_chapter_summary[,-(1:2)] <- diag(1/character_chapter_summary$words) %*% as.matrix(character_chapter_summary[,-(1:2)])
character_chapter_summary[,-(1:2)] <- apply(character_chapter_summary[,-(1:2)],2,normalize)
par(mfrow = c(4,4), mar = c(3,4,4,2))
for(ci in 2:ncol(character_chapter_summary)){
  plot(character_chapter_summary$chapter, character_chapter_summary[,ci], type = "l", xlab = "Chapter", main = toupper(colnames(character_chapter_summary)[ci]),
       ylab = "Normalized Character Density")
  lines(character_chapter_summary$chapter,
        dexp_smooth(character_chapter_summary$chapter, character_chapter_summary[,ci], r = 0.1),
        col = 2, lwd = 3)
}

melted_cormat <- reshape2::melt(reorder_cormat(cor(apply(cbind(sentiment_chapter_summary[,-1], character_chapter_summary[,-c(1:2)]), 2, diff), method = "spearman")))
melted_cormat$value <- round(melted_cormat$value, 2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

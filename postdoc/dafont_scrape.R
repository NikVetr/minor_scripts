# Load the necessary libraries
library(rvest)
library(dplyr)
library(purrr)

# The URL you want to scrape
url_parts <- c("https://www.dafont.com/top.php?page=", "&fpp=200")


# Use CSS selectors to scrape the font block. The CSS selector '.lv1left' is used to find the font block.

# Extract font names and urls
df <- lapply(1:5, function(i){
  
  print(i)
  
  url <- paste0(url_parts[1], i, url_parts[2])
  
  webpage <- read_html(url)

  font_block <- webpage %>%
    html_nodes(".lv1left") 
  
  df <- data.frame(
    font_names = font_block %>%
      html_nodes("a:nth-child(2)") %>% # The :nth-child(2) selector gets the second <a> tag within the div
      html_text(trim = TRUE),
    
    font_url = font_block %>%
      html_nodes("a:nth-child(2)") %>%
      html_attr("href") %>%
      map_chr(~ paste0("https://www.dafont.com/", .))
  )
  df
})

d <- do.call(rbind, df)


d$font_id <- gsub(".+.com/|.font.+", "", d$font_url)
d$dload_link <- paste0("https://dl.dafont.com/dl/?f=", d$font_id)
d$idn <- table(d$font_id)[d$font_id]
d$rank <- 1:nrow(d)
d <- d[order(d$font_id),]
d$occurrence <- unlist(lapply(rle(d$font_id)$lengths, function(i) 1:i))
d$occurrence[d$idn == 1] <- ""
d <- d[order(d$rank),]

#find font name embeddings
library(text)
font_name_embeddings <- parallel::mclapply(1:nrow(d), function(i){
    system(sprintf('echo "%s "', paste0(i, collapse="")))
    textEmbed(d$font_names[i])$x
  }, mc.cores = 12)
font_name_embeddings <- do.call(rbind, font_name_embeddings)
font_name_embeddings <- as.data.frame(font_name_embeddings)
rownames(font_name_embeddings) <- paste0(d$font_id, c("", ".")[(d$occurrence != "") + 1], d$occurrence)
colnames(font_name_embeddings) <- paste0("Dim", 1:ncol(font_name_embeddings))
PCA <- prcomp(font_name_embeddings)

write.csv(font_name_embeddings, "~/dafont_top1000_name-embeddings.txt", row.names = T)
data.table::fwrite(d, "~/dafont_top1000.txt")

d <- data.table::fread("~/dafont_top1000.txt")
font_name_embeddings <- read.csv("~/dafont_top1000_name-embeddings.txt", 
                                 row.names = 1)

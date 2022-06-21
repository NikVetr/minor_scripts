library(RSelenium)
library(pdftools)
library(tesseract)
library(data.table)
library(tabulizer)

zipcodeR::zipcodes <- search_radius(lat = 0, lng = 0, radius = 100000)
zipcodes <- zipcodes$zipcode

driver <- rsDriver(browser=c("chrome"))
remote_driver <- driver[["client"]]
remote_driver$open()

#try it out for real
# link <- 'https://990finder.foundationcenter.org/990results.aspx?990_type=A&fn=Animal&st=&zp=&ei=&fy=&action=Search'

year <- 2019
link_start <- "https://990finder.foundationcenter.org/990results.aspx?990_type=&fn=&st=&zp="
link_middle <- "&ei=&fy="
link_end <- "&action=Search"

for(zipi in 1:length(zipcodes)){
  
  link <- paste0(link_start, zipcodes[zipi], link_middle, year, link_end)
  
  remote_driver$navigate(link)
  html <- remote_driver$getPageSource()[[1]]
  html <- strsplit(html, "\n")[[1]]
  n_records <- html[grep(html, pattern = "Total results") + 1]
  n_records <- as.numeric(strsplit(strsplit(x = n_records, ">")[[1]][2], "<")[[1]][1])
  n_pages <- ceiling(n_records / 100)
  
  if(n_records == 0){next()}
  
  pdf_data <- data.frame()
  for(page_i in 1:n_pages){
    
    cat(paste0(zipi, " / ", length(zipcodes), "; ", page_i, " / ", n_pages, "\n"))
    
    #snag page info and process to retrive URLs
    html <- remote_driver$getPageSource()[[1]]
    html <- strsplit(html, "\n")[[1]]
    html <- html[grep("pdf", html)]
    pdf_urls <- as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = "href=\"//")[[1]][2], split = "\\.pdf")[[1]][1], ".pdf")))
    pdf_name <- as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = "pdf\">")[[1]][2], split = "</a>")[[1]][1])))
    pdf_money <- as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = "align=\"right\">")[[1]][3], split = "</td><td")[[1]][1])))
    pdf_state <- as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = "</a></td><td>")[[1]][2], split = paste0("</td><td>", year))[[1]][1])))
    pdf_zipcode <- rep(zipcodes[zipi], length(pdf_urls))
    pdf_year <- rep(year, length(pdf_urls))
    pdf_form <- trimws(as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = paste0(year, "</td><td>"))[[1]][2], split = "</td><td align=\"right\"")[[1]][1]))))
    pdf_EIN <- trimws(as.character(sapply(html, function(x) paste0(strsplit(strsplit(x, split = "style=\"white-space:nowrap;\">")[[1]][2], split = "</td>")[[1]][1]))))
    
    pdf_data_add <- data.frame(url = pdf_urls, name = pdf_name, total_assets = pdf_money, state = pdf_state, zipcode = pdf_zipcode, year = pdf_year, form = pdf_form, EIN = pdf_EIN)
    pdf_data <- rbind(pdf_data, pdf_data_add)
    
    #navigate to next page
    if(page_i < n_pages){
      webElem <- remote_driver$findElements(using = "link text", as.character(ifelse((page_i + 1) %% 10 == 1, "...", (page_i + 1))))
      webElem <- webElem[[ifelse((page_i + 1) %% 10 == 1 & (page_i + 1) > 20, 2, 1)]]
      webElem$clickElement()
    }  
    
    #wait for page to load
    Sys.sleep(1)
    
  }
  
  fwrite(pdf_data, file = paste0("G:\\990_pdfs\\", year, "_zipcodes\\", zipcodes[zipi], ".txt"))

}

pdf_dir <- "G:\\990_pdfs\\"
if(!dir.exists(pdf_dir)){dir.create(pdf_dir)}

year <- 2019
all_links <- data.table()
for(zipi in 1:length(zipcodes)){
  if(zipi %% 100 == 0){cat(paste0(zipi, " "))}
  if(file.exists(paste0(pdf_dir, year, "_zipcodes\\", zipcodes[zipi], ".txt"))){
    all_links <- rbind(all_links, fread(paste0(pdf_dir, year, "_zipcodes\\", zipcodes[zipi], ".txt")))
  }
}
fwrite(x = all_links, paste0(pdf_dir, year, "_links.txt"))


all_links <- fread(paste0(pdf_dir, year, "_links.txt"))
for(i in 1:nrow(all_links)){
  if(i %% 100 == 0){cat(paste0(i, " "))}
  if(!dir.exists(paste0(pdf_dir, year, "\\", all_links$zipcode[i]))){dir.create(paste0(pdf_dir, year, "\\", all_links$zipcode[i]))}
  try(download.file(paste0("https://", all_links$url[i]), 
                    destfile = paste0(pdf_dir, year, "\\", all_links$zipcode[i], "\\", 
                                      gsub(' ', "_", all_links$name[i]), "_EIN_", all_links$EIN[i], "_", all_links$form[i], ".pdf"), mode = "wb"))
}

# pdf_urls <- unique(pdf_urls)
# for(i in 1:length(pdf_urls)){
#   print(i)
#   try(download.file(paste0("https://", pdf_urls[i]), destfile = paste0(pdf_dir, i, ".pdf"), mode = "wb"))
# }

mgrep <- function(patterns, text){
  sapply(patterns, function(pattern) grep(pattern = pattern, x = text))
}

all_nonzero_length <- function(x){
  all(sapply(x, function(xi) length(xi) > 0))
}

#read pdfs back in and extract information from them
all_links <- fread(paste0(pdf_dir, year, "_links.txt"))
all_links <- all_links[all_links$form == "990",]
year <- 2019
i = sample(1:9000, 1)

file = paste0(pdf_dir, year, "\\", all_links$zipcode[i], "\\", 
       gsub(' ', "_", all_links$name[i]), "_EIN_", all_links$EIN[i], "_", all_links$form[i], ".pdf")


pdf_info <- pdf_info(file)
target_pages <- which(unlist(sapply(1:pdf_info$pages, function(page_i) all_nonzero_length(mgrep(apdf_machine[[page_i]]$text, patterns = c("Compensation", "Officer", "VII", "Contractors"))))))


#pdf_data solution
apdf_machine <- pdf_data(file)
a = apdf_machine[target_pages][[1]]
tl <- as.numeric(a[min(which(a$text == "Form")), c("x", "y")])
br <- as.numeric(a[max(which(a$text == "Form")), c("x", "y")])


#pdf_text attempt at solution
pdf_text_strings <- pdf_text(file)
pdf_text_strings_split <- strsplit((pdf_text_strings[[7]]), "\n")[[1]]
lines_with_dots <- grep("\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.", pdf_text_strings_split)
if(length(lines_with_dots) == 0){lines_with_dots <- grep("-------", pdf_text_strings_split)}
pdfd <- data.frame(
  name = pdf_text_strings_split[lines_with_dots-1],
  compensation = pdf_text_strings_split[lines_with_dots],
  title = pdf_text_strings_split[lines_with_dots+1]
)

trimws_internal <- function(x){
  x_trim <- lapply(x, function(xi) strsplit(xi, " ")[[1]])
  orig_ws <- lapply(x_trim, function(xi) which(sapply(xi, function(xii) xii == "")))

  #get rid of runs of length > 5 and sub in a space at that location
  orig_ws <- lapply(orig_ws, function(xi) xi[-(cumsum(rle(diff(xi))$lengths)[rle(diff(xi))$values == 1 & rle(diff(xi))$lengths > 5] - 1)]) 
  no_ws <- lapply(1:length(x_trim), function(xi) x_trim[[xi]][-orig_ws[[xi]]])
  
  t(sapply(no_ws, function(xi) c(paste0(xi[1:(which(xi == "")-1)], collapse = "_"), 
    paste0(xi[(which(xi == ""):length(xi))], collapse = ""))))
}
  
trimws_internal(pdfd$name)

# apdf_ocr <- pdf_ocr_data(paste0(pdf_dir, i, ".pdf"))
# out <- extract_tables(file, pages = target_pages, columns = list(6,6), guess = T, method = "decide")
# out

nameframe <- apdf_machine[[1]][abs(apdf_machine[[1]]$y - 156) < 3 & (apdf_machine[[1]]$x > 140 & apdf_machine[[1]]$x < 400 ),]
nameframe <- nameframe[order(nameframe$x),]
paste0(nameframe$text, collapse = " ")

apdf <- unlist(strsplit(apdf, "\n"))
intersect(apdf[grep(pattern = "Name", apdf)], apdf[grep(pattern = ":", apdf)])

findKeywords <- function(dir, patt){
  files <- list.files(dir)
  for(i in length(files)){
    doc <- readLines(paste0(dir, files[i]))
    if(any(grepl(pattern = patt, doc))){
      print(files[i])
      lines <- which(grepl(pattern = patt, readLines(paste0(dir, files[i]))))
      for(j in 1:length(lines)){
        print(paste0("\t\t", lines[j]))
      }
    }
  }
}
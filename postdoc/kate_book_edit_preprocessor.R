
#specify file paths
# source_dir <- "/Users/nikgvetr/Documents/Documents - nikolai/all_creatures_great_and_small/"
# source_dir <- "/Users/nikgvetr/Documents/Documents - nikolai/nothing_much_happens/"
source_dir <- "/Users/nikgvetr/Documents/Documents - nikolai/all_things_bright_and_beautiful/"
# source_dir <- "/Users/nikgvetr/Documents/Documents - nikolai/all_things_wise_and_wonderful/"
source_txt_files <- list.files(source_dir, full.names = T)[endsWith(list.files(source_dir), "txt")]
contents_file <- source_txt_files[endsWith(source_txt_files, "contents.txt")]
if(any(endsWith(source_txt_files, "processed.txt"))){
  processed_file <- source_txt_files[endsWith(source_txt_files, "processed.txt")]
  file.remove(processed_file)
  source_txt_files <- source_txt_files[!endsWith(source_txt_files, "processed.txt")]
}
input_file <- source_txt_files[!endsWith(source_txt_files, "contents.txt")]
book_title <- gsub("\\.txt", "", tail(strsplit(input_file, "/")[[1]], 1))
output_file <- gsub("\\.txt", "_processed\\.txt", input_file)

# Read the entire file as a single string
text <- readLines(input_file, warn = FALSE)
text <- trimws(text)
text_bl <- text == ""
text <- text[!text_bl]

# recover section headers
sections <- readLines(contents_file, warn = FALSE)
sections <- trimws(sections)
sections <- sections[sections != ""]

#if sections have some other format from the contents file...
if(mean(tolower(sections) %in% tolower(text)) < 0.5){
  sections <- JATSdecoder::text2num(sections)
}

# Split the text into sections, including the section headers
section_locs <- setNames(match(tolower(sections), tolower(text)), sections)
section_df <- data.frame(section = sections, 
                         start = as.numeric(section_locs),
                         end = c(as.numeric(section_locs)[-1] - 1, length(text))
)

title_info <- ""
if(section_df$start[1] > 1){
  title_info <- paste0(text[1:(section_df$start[1]-1)], collapse = "\n")
}

# Now merge them back as section-header + section-body pairs
insert_chapter <- F
insert_chapter_number <- F
keep_section_header <- F
change_allcaps_to_sentence_case <- T
fixed_sections <- c()
for (i in seq_along(sections)) {
  # Get the section header
  section_header <- sections[i]
  
  # Get the corresponding section body (even index in sections_split, after the first element)
  section_body <- text[section_df$start[i]:section_df$end[i]]
  
  #get appropriate case for the start of each section?
  if(change_allcaps_to_sentence_case){
    section_body[2] <- stringr::str_to_sentence(section_body[2])
  }
  
  # Fix line breaks within the section body
  # fixed_section_body <- gsub("([^.?!\"'”’\\)])[ \t]*\\n[ \t]*([^\n])", "\\1 \\2", section_body, perl = TRUE)
  fixed_section_body <- paste0(section_body[-1], collapse = "\n")
  
  # Add the fixed section (header + body) to the list
  section_header_name <- paste0(ifelse(insert_chapter, "Chapter ", ""), 
                                ifelse(insert_chapter_number, paste0(i, ": "), ""),
                                ifelse(keep_section_header, section_header, ""))
  
  fixed_sections <- c(fixed_sections, 
                      paste0(section_header_name, 
                             ifelse(section_header_name != "", "\n\n", ""), 
                             fixed_section_body))
}

# Recombine the fixed sections
fixed_sections[[1]] <- paste(c(title_info, fixed_sections[[1]]), collapse = "\n\n")
final_text <- paste(fixed_sections, collapse = "\n\n")

# Write the final, fixed text to the output file
writeLines(final_text, output_file)

#write the individual sections
section_names <- paste0(book_title, "-", 
                        gsub(" ", "-", tolower(gsub("\n| ", "_", sections))))
sections_dir <- paste0(source_dir, "sections/")
subsections_dir <- paste0(source_dir, "sections_sub/")
if(!dir.exists(sections_dir)){
  dir.create(sections_dir)
  dir.create(subsections_dir)
}

section_paths <- paste0(sections_dir, seq_along(section_names), "-", section_names,".txt")
subsection_path <- paste0(subsections_dir, 1, "-", section_names[1],".txt")
writeLines(fixed_sections[1], subsection_path)
for (i in seq_along(section_paths)) {
  writeLines(fixed_sections[i], section_paths[i])
}

nchars_per_section <- setNames(sapply(fixed_sections, nchar), section_names)
print(nchars_per_section)
setNames(paste0("$", round(nchars_per_section / 1194 * 0.015, 2)), section_names)

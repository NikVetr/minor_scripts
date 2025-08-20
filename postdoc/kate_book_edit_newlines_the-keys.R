
for(i in 1:2){
  
  # Define file paths
  if(i == 1){
    input_file <- "/Users/nikgvetr/Documents/Documents - nikolai/the_keys/the_keys_kate-gates.txt"    
  } else {
    input_file <- "/Users/nikgvetr/Documents/Documents - nikolai/the_keys/the_keys_kate-gates_fixed.txt"
  }
  
  output_file <- "/Users/nikgvetr/Documents/Documents - nikolai/the_keys/the_keys_kate-gates_fixed.txt"
  
  # Read the entire file as a single string
  text <- readLines(input_file, warn = FALSE)
  
  # Combine the text with explicit line breaks for easier processing
  text_combined <- paste(text, collapse = "\n")
  
  # Define the list of section headers, including the split lines
  sections <- c(
    "Prologue",
    "Part 1\nGathering",
    "Part 2\nCapture and Escape",
    "Part 3\nWhen Vardelien Dies",
    "Part 4\nCreatures of the Mines",
    "Part 5\nBattle of Darien… Yet again",
    "Part 6\nDissembled Feelings",
    "Epilogue"
  )
  
  # Create a single regex pattern to match any section header (preserves line breaks)
  section_pattern <- paste0("(", paste(sections, collapse = "|"), ")")
  
  # Split the text into sections, including the section headers
  sections_split <- strsplit(text_combined, section_pattern, perl = TRUE)[[1]]
  
  # Find the actual section headers (they will be the odd-indexed elements after the first)
  section_headers <- unlist(regmatches(text_combined, gregexpr(section_pattern, text_combined, perl = TRUE)))
  
  # Now merge them back as section-header + section-body pairs
  fixed_sections <- c()
  for (i in seq_along(section_headers)) {
    # Get the section header
    section_header <- section_headers[i]
    
    # Get the corresponding section body (even index in sections_split, after the first element)
    section_body <- sections_split[i + 1]
    
    # Fix line breaks within the section body
    fixed_section_body <- gsub("([^.?!\"'”’\\)])[ \t]*\\n[ \t]*([^\n])", "\\1 \\2", section_body, perl = TRUE)
    
    # Add the fixed section (header + body) to the list
    fixed_sections <- c(fixed_sections, paste0(section_header, fixed_section_body, collapse = "\n\n"))
  }
  
  # Recombine the fixed sections
  final_text <- paste(c("The Keys, by Katherine Gates.\n\n", fixed_sections), collapse = "\n\n")
  
  # Write the final, fixed text to the output file
  writeLines(final_text, output_file)
  
}

#write the individual sections
section_names <- paste0("the-keys_", gsub(" ", "-", tolower(gsub("\n", "_", section_headers))))
section_paths <- paste0("/Users/nikgvetr/Documents/Documents - nikolai/the_keys/sections/", section_names,".txt")
for (i in seq_along(section_headers)) {
  section_path <- section_paths[i]
  writeLines(fixed_sections[i], section_path)
}

nchars_per_section <- setNames(sapply(fixed_sections, nchar), section_names)
print(nchars_per_section)
setNames(paste0("$", round(nchars_per_section / 1194 * 0.015, 2)), section_names)

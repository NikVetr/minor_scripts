# library(xlsx)
library(openxlsx)

#remove all currently loaded objects
curr_objects <- ls(); curr_objects <- ls()
rm(list = curr_objects)

# Define the directory containing the RData files
file_dir <- "~/motrpac_companion_r3/source_data/"

# Get the paths of the RData files
source_data_paths <- list.files(path = file_dir, full.names = TRUE)
if(any(endsWith(source_data_paths, "xlsx"))){
  file.remove(source_data_paths[endsWith(source_data_paths, "xlsx")])
  source_data_paths <- source_data_paths[!endsWith(source_data_paths, "xlsx")]
}

#do some quick processing of write.xlsx
# write.xlsx.2 <- function(x, file, sheetName, append = T){
#   has_rownames <- !is.null(rownames(x))
#   has_colnames <- !is.null(colnames(x))
#   write.xlsx(x = x, file = file, 
#              sheetName = sheetName, append = append, 
#              col.names = has_colnames, row.names = has_rownames, showNA = F)
# }

write.xlsx.2 <- function(x, file, sheetName, obj_class, obj_name) {
  # Check if the file already exists
  if (file.exists(file)) {
    # Load existing workbook if appending
    wb <- loadWorkbook(file)
  } else {
    # Create a new workbook if not appending
    wb <- createWorkbook()
  }
  
  #add obj_class and object name to sheet
  x <- as.data.frame(x)
  x <- cbind(obj_name = c(obj_name, rep(NA, nrow(x)-1)),
             obj_class = c(obj_class, rep(NA, nrow(x)-length(obj_class))),
             NA,
             x
  )
  colnames(x)[3] = "."
  
  # Add a new worksheet with the specified sheet name
  if(nchar(sheetName) > 30){sheetName <- substr(sheetName, 1, 30)}
  addWorksheet(wb, sheetName)
  
  # Write data to the worksheet
  writeData(wb, sheetName, x, colNames = !is.null(colnames(x)), rowNames = !is.null(rownames(x)), )
  
  # Save the workbook
  saveWorkbook(wb, file, overwrite = TRUE)
}


#shorten sheet substrings (31 char limit total)
shorten_string <- function(string){
  string <- ifelse(nchar(string) < 8, 
                   string, 
                   paste0(substr(strsplit(string, "_|-")[[1]], 1, 3), collapse = "_"))
  if(nchar(string) > 8){
    string <- substr(string, 1, 9)
  }
  return(string)
}

ragged_to_df <- function(list_of_vectors) {
  # Find the maximum length of the vectors in the list
  max_length <- max(sapply(list_of_vectors, length))
  
  # Function to pad a vector with NAs to match the max length
  pad_vector <- function(vec, max_len) {
    length(vec) <- max_len
    vec
  }
  
  # Apply padding to each vector in the list
  padded_list <- lapply(list_of_vectors, pad_vector, max_length)
  
  # Convert the padded list to a data frame
  data.frame(padded_list)
}

# process lists function
write_list <- function(lst, file_name, base_sheet_name, obj_class, obj_name) {
  
  for (elem_name in names(lst)) {
    
    #retrieve element
    elem <- lst[[elem_name]]
    
    #construct the sheet name
    sheet_name <- paste0(base_sheet_name, ".", shorten_string(elem_name))
    
    #check if list
    if(any(class(elem) %in% "histogram")){
      write.xlsx.2(x = ragged_to_df(elem), file = file_name, sheetName = sheet_name, obj_class = obj_class, obj_name = obj_name)
    } else if(any(class(elem) %in% c("list"))){
      write_list(lst = elem, 
                 file_name = file_name, 
                 base_sheet_name = sheet_name, obj_class = obj_class, obj_name = obj_name)
    #check if array
    } else if (is.array(elem) && length(dim(elem)) > 2) {
      write_slices(elem, file_name, sheet_name, obj_class = obj_class, obj_name = obj_name)
    #check if vector or table and convert + write if that is the case
    } else if (is.vector(elem) | any(class(elem) %in% c("table"))) {
      elem <- as.data.frame(elem)
      write.xlsx.2(x = elem, file = file_name, sheetName = sheet_name, obj_class = obj_class, obj_name = obj_name)
    #check if matrix or df and write
    } else if (is.matrix(elem) || is.data.frame(elem)) {
      write.xlsx.2(x = elem, file = file_name, sheetName = sheet_name, obj_class = obj_class, obj_name = obj_name)
    } else {
    #throw error for debugging
      stop(paste0("ERROR: ", sheet_name, " cannot be written in ", file_name))
    }
    
  }
}

#process arrays
write_slices <- function(obj, file_name, base_sheet_name, obj_class, obj_name) {
  obj_dims <- dim(obj)
  
  if (is.list(obj)) {
    # Handle lists separately
    write_list(obj, file_name, base_sheet_name, obj_class = obj_class, obj_name = obj_name)
  } else if (length(obj_dims) > 2) {
    array_list <- apply(obj, 1, function(x) list(x))
    for (i in seq_along(array_list)) {
      new_obj <- array_list[[i]][[1]]
      new_base_sheet_name <- paste0(base_sheet_name, ".", i)
      write_slices(new_obj, file_name, new_base_sheet_name, obj_class = obj_class, obj_name = obj_name)
    }
  } else {
    write.xlsx.2(as.matrix(obj), file = file_name, sheetName = base_sheet_name, obj_class = obj_class, obj_name = obj_name)
  }
}

base_objects <- ls(); base_objects <- ls() #call twice to get self in there
#### Iterate through each RData file ####
for (path in source_data_paths) {
  cat(paste0("\n\n", path, ": "))
  
  #load the RData file
  load(path)
  
  #get all object names in the environment
  all_objects <- ls()
  all_objects <- setdiff(all_objects, base_objects)
  
  #filter out functions
  compatible_objects <- all_objects[!sapply(all_objects, function(x) is.function(get(x)))]
  
  #filter out indices
  compatible_objects <- setdiff(compatible_objects, c("i", "j", "k"))
  obj_lengths <- sapply(compatible_objects, function(x) length(get(x)))
  obj_name_lengths <- nchar(compatible_objects)
  has_i_in_name <- grepl(x = compatible_objects, "i")
  ends_in_underscore_ind <- substr(compatible_objects, obj_name_lengths - 1, obj_name_lengths) %in% c("_i", "_j", "_k")
  object_sizes <- sapply(compatible_objects, function(x) object.size(get(x)))
  compatible_objects <- compatible_objects[!((obj_lengths == 1) &
                                               ( ((obj_name_lengths == 2) & has_i_in_name) | ends_in_underscore_ind ) & 
                                               (object_sizes < 1000))]
  
  
  #sort by object size
  object_sizes <- sapply(compatible_objects, function(x) object.size(get(x)))
  compatible_objects <- compatible_objects[order(object_sizes, decreasing = T)]
  
  # Define the Excel file name based on the RData file
  excel_file_name <- paste0(file_dir, tools::file_path_sans_ext(basename(path)), ".xlsx")
  if(file.exists(excel_file_name)){
    file.remove(excel_file_name)
  }
  
  # Process each object
  for (obj_name in compatible_objects) {
    cat(paste0(obj_name, ", "))
    
    obj <- get(obj_name)
    obj_class <- class(obj)
    sheet_name <- shorten_string(obj_name)
    
    if(any(obj_class %in% "list")){
      write_list(lst = obj, file_name = excel_file_name, base_sheet_name = sheet_name, obj_class = obj_class, obj_name = obj_name)
    } else if (is.array(obj) && length(dim(obj)) > 2) {
      write_slices(obj, excel_file_name, sheet_name, obj_class = obj_class, obj_name = obj_name)
    } else {
      # Write 2D objects directly to the Excel file
      if (is.vector(obj)) {
        obj <- matrix(obj, ncol = 1)  # Convert to a single-column matrix
      }
      if (is.array(obj)) {
        obj <- as.matrix(obj)  # Convert to a matrix, bc 2D arrays play poorly with write.xlsx
      }
      
      #thin mcmc samples for excel to behave
      if(nrow(obj) > 1E3){
        obj <- obj[seq(1, nrow(obj), length.out = 1E3),]
      }
      write.xlsx.2(x = obj, file = excel_file_name, sheetName = obj_name, obj_class = obj_class, obj_name = obj_name)
    }
  }
  
  # Clear the environment
  objects_to_remove <- setdiff(ls(), base_objects)
  rm(list = objects_to_remove)
  gc()  # Garbage collection to free up memory
}

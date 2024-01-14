library(xlsx)

# Define the directory containing the RData files
file_dir <- "~/motrpac_companion_r3/source_data/"

# Get the paths of the RData files
source_data_paths <- list.files(path = file_dir, full.names = TRUE)

# Function to process multi-dimensional arrays
write_slices <- function(obj, file_name, base_sheet_name) {
  obj_dims <- dim(obj)
  
  if (length(obj_dims) > 2) {
    # Iterate through the first dimension
    for (i in seq(obj_dims[1])) {
      # Slice the object for the current dimension
      new_obj <- obj[i, , , drop = FALSE]
      
      # Construct new base sheet name for the slice
      new_base_sheet_name <- paste0(base_sheet_name, ".", i)
      
      # Recursive call with the sliced object
      write_slices(new_obj, file_name, new_base_sheet_name)
    }
  } else {
    # Write the 2D object to an Excel sheet
    write.xlsx(obj, file = file_name, sheetName = base_sheet_name, row.names = FALSE, append = TRUE)
  }
}

base_objects <- ls()
base_objects <- ls()

# Iterate through each RData file
for (path in source_data_paths) {
  # Load the RData file
  load(path)
  
  # Get all object names in the environment
  all_objects <- ls()
  all_objects <- setdiff(all_objects, base_objects)
  
  # Filter out non-compatible objects
  compatible_objects <- all_objects[!sapply(all_objects, function(x) is.function(get(x)))]
  
  # Define the Excel file name based on the RData file
  excel_file_name <- paste0(file_dir, tools::file_path_sans_ext(basename(path)), ".xlsx")
  
  # Process each object
  for (obj_name in compatible_objects) {
    obj <- get(obj_name)
    
    # Check if the object is a multi-dimensional array
    if (is.array(obj) && length(dim(obj)) > 2) {
      write_slices(obj, excel_file_name, obj_name)
    } else {
      # Write 2D objects directly to the Excel file
      if (is.vector(obj) && is.character(obj)) {
        obj <- matrix(obj, ncol = 1)  # Convert to a single-column matrix
      }
      write.xlsx(obj, file = excel_file_name, sheetName = obj_name, row.names = FALSE, append = TRUE)
    }
  }
  
  # Clear the environment
  rm(list = compatible_objects)
  gc()  # Garbage collection to free up memory
}

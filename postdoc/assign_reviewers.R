nms <- "Christine	Nik
Anam	Nik
Christine	Mark
Nik	Anam
Anam	Christine
Nik	Mark
Anam	Mark
Christine	Nik
Nik	Anam
Anam	Christine
Nik	Anam
Nik	Anam"

nms <- matrix(strsplit(nms, "\n|\t")[[1]], ncol = 2, byrow = T)
reviewers <- c("Nik", "Anam", "Christine", "Ignacio")


# Build list of valid reviewers per row
valid_reviewers <- apply(nms, 1, function(row) setdiff(reviewers, row))

# Try multiple assignments until a balanced one is found
assign_reviewer <- function(valid_reviewers, reviewers, n_per_reviewer = 3, max_tries = 10000) {
  for (try in 1:max_tries) {
    # Create pool: 3 of each reviewer
    pool <- rep(reviewers, each = n_per_reviewer)
    assignment <- rep(NA, length(valid_reviewers))
    
    # Shuffle the order in which we try to assign
    row_order <- sample(seq_along(valid_reviewers))
    
    success <- TRUE
    for (i in row_order) {
      options <- intersect(valid_reviewers[[i]], pool)
      if (length(options) == 0) {
        success <- FALSE
        break
      }
      chosen <- sample(options, 1)
      assignment[i] <- chosen
      pool <- pool[pool != chosen | duplicated(pool)]  # remove one instance
    }
    
    if (success) return(assignment)
  }
  stop("Failed to find a valid assignment after many tries.")
}

assigned <- assign_reviewer(valid_reviewers, reviewers)

# Combine result
result <- cbind(nms, Reviewer = assigned)
print(result)
cat(paste0(assigned, sep = "\n"))

#check constraints
table(result[,3])
all(!apply(result, 1, function(x) x[3] %in% x[1:2]))

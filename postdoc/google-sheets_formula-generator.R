generate_formula <- function(col, prev_col, swap_row, start_row, end_row) {
  paste0(
    "=IF(OR(ISBLANK(M$", swap_row, "), ISBLANK(N$", swap_row, ")), ", prev_col, start_row, 
    ", IF(", prev_col, start_row, "=M$", swap_row, 
    ", INDEX(", prev_col, "$", start_row, ":", prev_col, "$", end_row, 
    ", MATCH(N$", swap_row, ", ", prev_col, "$", start_row, ":", prev_col, "$", end_row, ", 0)), IF(", 
    prev_col, start_row, "=N$", swap_row, 
    ", INDEX(", prev_col, "$", start_row, ":", prev_col, "$", end_row, 
    ", MATCH(M$", swap_row, ", ", prev_col, "$", start_row, ":", prev_col, "$", end_row, ", 0)), ", 
    prev_col, start_row, ")))"
  )
}

columns <- LETTERS[which(LETTERS=="M"):which(LETTERS=="Z")]
swap_rows <- 22:34
start_row <- 76
end_row <- 103

formulas <- sapply(1:14, function(i) {
  generate_formula(columns[i + 1], columns[i], swap_rows[i], start_row, end_row)
})

cat(paste0(formulas, collapse = "\n"))

=IF(OR(ISBLANK(M$22), ISBLANK(N$22)), 
    M76, 
    IF(M76=M$22, 
       INDEX(M$76:M$103, MATCH(N$22, M$76:M$103, 0)), 
       IF(M76=N$22, 
          INDEX(M$76:M$103, MATCH(M$22, M$76:M$103, 0)), 
          M76)))

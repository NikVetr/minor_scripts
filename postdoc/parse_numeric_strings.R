n <- 1E3
nums <- sapply(abs(rcauchy(n, 0, 5)), function(x) paste0(round(x, sample(0:3, 1, T)), sample(c("K", "M", "B", "T"), 1)))

multi_gsub <- function(pattern, replacement, x){
  if(length(pattern) != 0){
    multi_gsub(pattern[-1], replacement, gsub(pattern[1], replacement, x))
  } else {
    return(x)
  }
}

conv_str_to_num <- function(nums){
  nums_num <- multi_gsub(c(LETTERS, letters), replacement = "", x = nums)
  nums_char <- substr(nums, nchar(nums_num)+1, nchar(nums))
  numbers <- as.numeric(nums_num) * 10^(1:4*3)[match(nums_char, c("K", "M", "B", "T"))]
  return(numbers)
}
  
nums
conv_str_to_num(nums)

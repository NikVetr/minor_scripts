#make it replicable
set.seed(1)
match_index <- 4

#load libraries
library(rJava)
library(mailR)
library(rvest)
library(dplyr)
library(data.table)
library(purrr)


#snag staff & board bios + contact info
staff_url <- "https://www.wildanimalinitiative.org/staff"
staff_html <- read_html(staff_url)

extract_individual_info <- function(node) {
  bio <- html_nodes(node, "div.accordion-item__description p") %>%
    html_text(trim = TRUE) %>%
    paste(collapse = " ") # Collapses the multiple p elements into a single string
  
  # Extract email addresses from the href attribute and remove the 'mailto:' prefix
  email <- html_attr(html_node(node, "a[href^='mailto:']"), "href") %>% 
    sub("mailto:", "", .)
  
  tibble(bio = bio, email = email)
}

# Extract the staff information
staff_df <- staff_html %>%
  html_nodes("li.accordion-item") %>%    
  map_df(extract_individual_info)
staff_df$names <- (staff_html %>% html_nodes("div.sqs-html-content h4") %>% html_text(trim = TRUE))[1:nrow(staff_df)]
staff_df$bio <- gsub(" Contact.*", "",  staff_df$bio)
staff_df$bio <- sapply(1:nrow(staff_df), function(i) trimws(gsub(staff_df$email[i], "", staff_df$bio[i])))
staff_df$bio <- gsub("’", "\'",  staff_df$bio)


#remove Cam from consideration
staff_df <- staff_df[staff_df$names != "Cameron Meyer Shorb",]

#now do board
board_url <- "https://www.wildanimalinitiative.org/board"
board_html <- read_html(board_url)

board_df <- board_html %>%
  html_nodes("li.accordion-item") %>%    
  map_df(extract_individual_info)
board_df$names <- (board_html %>% html_nodes("div.sqs-html-content h4") %>% html_text(trim = TRUE))[1:nrow(board_df)]
board_df$email <- c("nikgvetr@stanford.edu",
                    "joshyou12@gmail.com",
                    "anami.nguyen@protonmail.com",
                    "itmoore@vt.edu",
                    "christinedperry337@gmail.com")
board_df$bio <- gsub(" Learn more →", "", board_df$bio)
board_df$bio <- gsub("’", "\'",  board_df$bio)


staff <- setNames(object = staff_df$email,
                  nm = sapply(strsplit(staff_df$names, " "), `[`, 1))
board <- setNames(object = board_df$email,
                  nm = sapply(strsplit(board_df$names, " "), `[`, 1))
both <- c(staff, board)
both_inv <- setNames(names(both), both)

#read in old meetings
old_matches <- do.call(rbind, lapply(1:(match_index-1), function(i) data.table::fread(paste0("~/Documents/Documents - nikolai/WAI_Chat_", i,".txt"))))

#match em together
possible_matches <- data.frame(board = rep(board, times = length(staff)), 
                               staff = rep(staff, each = length(board)))
possible_matches <- possible_matches[!apply(possible_matches, 1, paste0, collapse = "~") %in% apply(old_matches, 1, paste0, collapse = "~"),]
possible_matches <- possible_matches[sample(1:nrow(possible_matches)),]
matches <- data.frame(board = NA, staff = staff)
matches$board <- NA
remaining_matches <- possible_matches
for(i in 1:nrow(matches)){
  si <- matches$staff[i]
  opt <- sample(which(remaining_matches$staff == si))
  pair <- remaining_matches$board[opt[1]]
  matches$board[i] <- pair
  remaining_matches <- remaining_matches[!((remaining_matches$staff == si) | (remaining_matches$board == pair)),]
  if(i %% length(board) == 0){
    remaining_matches <- possible_matches
  }
}

#confirm there are no duplicates
sapply(1:nrow(matches), function(i){
  any((matches$board[i] == old_matches$board) & (matches$staff[i] == old_matches$staff))
})

#write current matches to file
data.table::fwrite(matches, paste0("~/Documents/Documents - nikolai/WAI_Chat_", match_index,".txt"))

#compile email data
email_df <- data.table(greetings = sapply(1:nrow(matches), function(i) paste0("Hi ", both_inv[matches$board[i]], " and ", both_inv[matches$staff[i]],"!\n\n")),
                       email1 = matches$board[i], email2 = matches$staff[i])
new_here <- setdiff(both, unique(unlist(old_matches)))
new_here <- setdiff(new_here, "mal.graham@wildanimalinitiative.org")

#write the emails
old_ways <- F
if(old_ways){
  message_body <- "So as Michelle mentioned, one of our goals this year is to better facilitate communication between WAI staff and board.\n
  To that end, we (the board) were thinking to try pairing members of each up for two short, 15-30 minute chats. I went ahead and matched the six board members to the six staff members at uniform, and am sending this email out to each pair. Whichever member of your pair sees this email first, please share your availability to chat with your partner sometime during this next month and schedule a meeting, preferably via Zoom or other video conferencing app.\n
  Before meeting, staff members can consult information provided in the attached document about board members, and board members can look here (https://www.wildanimalinitiative.org/about-us) for more information about staff. Please come prepared with two questions about your chat partner's specific role at the WAI.\n
  Happy chatting! 
  -Nik\n"
  
  message_body <- "As mentioned in the last email, we were hoping to orchestrate two of these paired chat sessions, and so here is the long-awaited sequel!\n
  Like before, these should take between 15-30 minutes. Staff can consult the attached document for info about your randomly chosen board member, and board members can look at information contained at this link: https://www.wildanimalinitiative.org/about-us.\n
  If you're reading this, please reach out to your chat partner with a few good times. Happy chatting!\n
  Best,\nNik"
  
  message_body <- sapply(1:nrow(matches), function(i) {
    new_peeps <- unlist(intersect(new_here, matches[i,]))
    paste0("We're on to the next Volume ", utils::as.roman(match_index), " of our bi-annual WAI Staff + Board Meet & Greet!\n\n",
           ifelse(length(new_peeps) == 0,"As a quick recap: the goal here is to", 
           paste0(ifelse(length(new_peeps) == 1, both_inv[new_peeps[1]], paste0(both_inv[new_peeps[1]], " and ", both_inv[new_peeps[2]])), " -- it looks like", 
                  ifelse(length(new_peeps) == 0, " neither of", "")," you've never done one of these before.",
                  " In 2021, we thought that to better facilitate board-staff communication, we could do a little icebreaker where we would")), 
           " randomly pair staff members up with board members for informal, 15-20 minute chats.\n\n",
           "What will you talk about? Whatever you want! Your interests, hobbies, goals? What brought you to WAI? Your favorite animals? Sky's the limit.\n\n",
           "Staff can consult the attached document for info about your randomly chosen board member, and board members can look at staff information contained at this link: https://www.wildanimalinitiative.org/team\n\n",
           "Whoever's reading this first, please reach out to your chat partner with your best available times & your preferred communication medium. Happy chatting!\n\nBest,\nNik")
    })
  
  #combine emails
  emails <- lapply(1:nrow(matches), function(i) paste0(email_df$greetings[i], message_body[i]))

} else { #made with GPT assistance
  
  # Function to generate message body for each pair
  generate_message <- function(newcomers, match_index, pair_names, staff_df, board_df, email_addresses) {
    
    message <- paste0("Hi ", pair_names[1], " and ", pair_names[2],"!\n\n")
    
    if (length(newcomers) == 2) {
      message <- paste0(message, 
        "Welcome to your first WAI Staff + Board Meet & Greet! My records indicate that both ", 
        pair_names[1], " and ", pair_names[2], " joined since the last time I sent out these emails, which is to say they've never done one of these before. Very exciting! \n\nThe "
      )
    } else if (length(newcomers) == 1) {
      message <- paste0(message, 
        "Welcome to Volume ", utils::as.roman(match_index), " of our bi-annual WAI Staff + Board Meet & Greet!\n\n",
        
        "We're thrilled to have ", pair_names[newcomers[1]], " joining us this round. For our returning participant, ", 
        pair_names[setdiff(unlist(email_addresses), newcomers)], ", we're glad to have you back!\n\nIf you recall, ", 
        pair_names[setdiff(unlist(email_addresses), newcomers)], ", the "
      )
    } else {
      message <- paste0(message, 
        "Welcome back to Volume ", utils::as.roman(match_index), " of our bi-annual WAI Staff + Board Meet & Greet!\n\n",
        
        "It's great to have both ", pair_names[1], " and ", pair_names[2], " back for another round of insightful chats.\n\nAs a reminder, the "
      )
    }
    
    message <- paste0(message, 
    
    "aim of these 15-20 minute chats is to foster open dialogue and a strong relationship between our staff and board members. You can talk about whatever you want - your interests, hobbies, goals, what brought you to WAI, your favorite animals, or anything else that comes to mind.\n\n",
    
    "I've taken the liberty of snagging both of your current bios from https://www.wildanimalinitiative.org/. Here they are reproduced below:\n\n\t",
    
    pair_names[email_addresses$staff], ": ", staff_df$bio[staff_df$email == email_addresses$staff], "\n\n\t",
    
    pair_names[email_addresses$board], ": ", board_df$bio[board_df$email == email_addresses$board], "\n\n",
    
    "I'm also attaching a short document that describes the general purpose of WAI's board.\n\n",
    
    "Whichever of you reads this email first, please reach out to your chat partner with your best available times and your preferred communication medium.\n\n",
    
    "Happy chatting!\nNik")
    
    
    return(message)
  }
  
  # Generate the emails
  email_df <- data.table(email1 = matches$board, email2 = matches$staff)
  
  emails <- sapply(1:nrow(matches), function(i) {
    pair_names <- c(both_inv[matches$board[i]], both_inv[matches$staff[i]])
    newcomers <- intersect(new_here, unlist(matches[i,]))
    generate_message(newcomers, match_index, pair_names, staff_df, board_df, matches[i,])
  })
  
}

cat(emails[[1]])

#snag password
pass <- readLines("~/data/test_letters.txt") #generate from https://myaccount.google.com/u/1/apppasswords

for(i in 1:nrow(matches)){
# for(i in 7:8){
  print(i)
  Sys.sleep(1)
  send.mail(from = "nikolai.vetr@gmail.com",
            to = unlist(matches[i,]),
            # to = c("nlashinsky@ucdavis.edu", "nikolaiveter@gmail.com"),
            subject = paste0("WAI Staff + Board Paired Meet & Greet, Vol. ", utils::as.roman(match_index)),
            body = emails[[i]],
            smtp = list(host.name = "smtp.gmail.com", port = 587,
                        user.name = "nikolai.vetr@gmail.com",
                        passwd = pass, ssl = TRUE),
            authenticate = TRUE,
            send = TRUE,
            attach.files = c("~/data/WAI_Board_Profile_and_Purpose_Summary.pdf", "~/data/wai_logo.jpg"))
}

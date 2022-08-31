#make it replicable
set.seed(1)
match_index <- 3

#load libraries
library(rJava)
library(mailR)
library(data.table)

#specify staff & board contact info
staff <- setNames(object = c("luke.hecht@wildanimalinitiative.org", 
           "simon.liedholm@wildanimalinitiative.org", 
           "vittoria.elliott@wildanimalinitiative.org",
           "tanaiia.hall@wildanimalinitiative.org",
           "cat.kerr@wildanimalinitiative.org",
           "amy.klarup@wildanimalinitiative.org",
           "jason.orlando@wildanimalinitiative.org",
           "emily.sharp@wildanimalinitiative.org",
           "mark.onley@wildanimalinitiative.org",
           "suzanne.vanarsdale@wildanimalinitiative.org"),
           nm = c("Luke", "Simon", "Vittoria", "Tanaiia", "Cat", "Amy", "Jason", "Emily", "Mark", "Suzanne"))
board <- setNames(object = c("itmoore@vt.edu",
            "nikgvetr@stanford.edu",
            "joshyou12@gmail.com",
            "stien.vanderploeg@wildanimalinitiative.org",
            "christinedperry337@gmail.com"),
            nm = c("Ignacio", "Nik", "Josh", "Stien", "Christine"))
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
  opt <- which(remaining_matches$staff == si)
  pair <- remaining_matches$board[opt[1]]
  matches$board[i] <- pair
  remaining_matches <- remaining_matches[!((remaining_matches$staff == si) | (remaining_matches$board == pair)),]
  if(i %% length(board) == 0){
    remaining_matches <- possible_matches
  }
}

#write current matches to file
data.table::fwrite(matches, paste0("~/Documents/Documents - nikolai/WAI_Chat_", match_index,".txt"))

#compile email data
email_df <- data.table(greetings = sapply(1:nrow(matches), function(i) paste0("Hi ", both_inv[matches$board[i]], " and ", both_inv[matches$staff[i]],"!\n\n")),
                       email1 = matches$board[i], email2 = matches$staff[i])

message_body <- "So as Michelle mentioned, one of our goals this year is to better facilitate communication between WAI staff and board.\n
To that end, we (the board) were thinking to try pairing members of each up for two short, 15-30 minute chats. I went ahead and matched the six board members to the six staff members at uniform, and am sending this email out to each pair. Whichever member of your pair sees this email first, please share your availability to chat with your partner sometime during this next month and schedule a meeting, preferably via Zoom or other video conferencing app.\n
Before meeting, staff members can consult information provided in the attached document about board members, and board members can look here (https://www.wildanimalinitiative.org/about-us) for more information about staff. Please come prepared with two questions about your chat partner's specific role at the WAI.\n
Happy chatting! 
-Nik\n"

message_body <- "As mentioned in the last email, we were hoping to orchestrate two of these paired chat sessions, and so here is the long-awaited sequel!\n
Like before, these should take between 15-30 minutes. Staff can consult the attached document for info about your randomly chosen board member, and board members can look at information contained at this link: https://www.wildanimalinitiative.org/about-us.\n
If you're reading this, please reach out to your chat partner with a few good times. Happy chatting!\n
Best,\nNik"

new_here <- setdiff(both, unique(unlist(old_matches)))
message_body <- sapply(1:nrow(matches), function(i) {
  new_peeps <- unlist(intersect(new_here, matches[i,]))
  paste0("So we're up to Volume ", utils::as.roman(match_index), " of our bi-annual WAI Staff + Board Meet & Greet!\n\n",
         ifelse(length(new_peeps) == 0,"As a quick recap: the goal here is to", 
         paste0(ifelse(length(new_peeps) == 1, both_inv[new_peeps[1]], paste0(both_inv[new_peeps[1]], " and ", both_inv[new_peeps[2]])), " -- it looks like", 
                ifelse(length(new_peeps) == 0, " neither of", "")," you've never done one of these before.",
                " In 2021, we thought that to better facilitate board-staff communication, we could do a little icebreaker where we would")), 
         " randomly pair staff members up with board members for informal, 20-30 minute chats.\n\n",
         "What will you talk about? Whatever you want! Your interests, hobbies, goals? What brought you to WAI? Your favorite animals? Sky's the limit.\n\n",
         "Staff can consult the attached document for info about your randomly chosen board member, and board members can look at staff information contained at this link: https://www.wildanimalinitiative.org/team\n\n",
         "Whoever's reading this first, please reach out to your chat partner with your best available times & your preferred communication medium. Happy chatting!\n\nBest,\nNik")
  })

#combine emails
emails <- lapply(1:nrow(matches), function(i) paste0(email_df$greetings[i], message_body[i]))

#snag password
pass <- readLines("~/data/test_letters.txt") #generate from https://myaccount.google.com/u/1/apppasswords

for(i in 1:nrow(matches)){
# for(i in 3:4){
  print(i)
  Sys.sleep(1)
  send.mail(from = "nikolai.vetr@gmail.com",
            to = c(unlist(matches[i,]), "nikolai.vetr@gmail.com"),
            subject = paste0("WAI Staff + Board Paired Meet & Greet, Vol. ", utils::as.roman(match_index)),
            body = emails[[i]],
            smtp = list(host.name = "smtp.gmail.com", port = 587,
                        user.name = "nikolai.vetr@gmail.com",
                        passwd = pass, ssl = TRUE),
            authenticate = TRUE,
            send = TRUE,
            attach.files = "~/data/WAI_Board_Profile_and_Purpose_Summary.pdf")
}

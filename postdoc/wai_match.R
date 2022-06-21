#make it replicable
set.seed(1)

#load libraries
library(gmailr) #I used this but see https://mailtrap.io/blog/r-send-email/ for more options, or https://gmailr.r-lib.org/ for usage guide
library(data.table)

#specify staff & board contact info
staff <- cbind(emails = c("hollis.howe@wildanimalinitiative.org", 
           "luke.hecht@wildanimalinitiative.org", 
           "simon.liedholm@wildanimalinitiative.org", 
           "cameronms@wildanimalinitiative.org",
           "vittoria.elliott@gmail.com",
           "tanaiia.hall@wildanimalinitiative.org"),
           names = c("Hollis", "Luke", "Simon", "Cam", "Vittoria", "Tanaiia"))
board <- cbind(emails = c("itmoore@vt.edu",
            "nikgvetr@stanford.edu",
            "joshyou12@gmail.com",
            "stien.vanderploeg@wildanimalinitiative.org",
            "christinedperry337@gmail.com",
            "ehatch475@gmail.com"),
            names = c("Ignacio", "Nik", "Josh", "Stien", "Christine", "Emily"))

#match em together
matches <- cbind(sample(board[,"names"]), sample(staff[,"names"]))
matches <- cbind(matches, matches[c(length(staff[,"names"]), 1:(length(staff[,"names"]) - 1)),2]) #get second round
matches

#compile email data
email_df <- data.table(greetings = sapply(1:nrow(matches), function(i) paste0("Hi ", matches[i,1], " and ", matches[i,2],"!\n\n")),
                       email1 = sapply(1:nrow(matches), function(i) board[match(matches[i,1], board[,2]),1]),
                       email2 = sapply(1:nrow(matches), function(i) staff[match(matches[i,2], staff[,2]),1]))

message_body <- "So as Michelle mentioned, one of our goals this year is to better facilitate communication between WAI staff and board.\n
To that end, we (the board) were thinking to try pairing members of each up for two short, 15-30 minute chats. I went ahead and matched the six board members to the six staff members at uniform, and am sending this email out to each pair. Whichever member of your pair sees this email first, please share your availability to chat with your partner sometime during this next month and schedule a meeting, preferably via Zoom or other video conferencing app.\n
Before meeting, staff members can consult information provided in the attached document about board members, and board members can look here (https://www.wildanimalinitiative.org/about-us) for more information about staff. Please come prepared with two questions about your chat partner's specific role at the WAI.\n
Happy chatting! 
-Nik\n"

#authenticate access, get credentials (json file or key / secret pair) from https://console.cloud.google.com/apis/credentials/consent
gm_auth_configure(path = "~/nikolai.vetr@gmail.com.json")
gm_auth(email = "nikolai.vetr@gmail.com")

#construct email objects
emails <- lapply(1:nrow(email_df), function(i){
  gm_mime() %>%
  gm_to(c("nikgvetr@stanford.edu", "nlashinsky@ucdavis.edu")) %>%
  # gm_to(c(email_df$email1[i], email_df$email2[i])) %>%
  gm_from("nikolai.vetr@gmail.com") %>%
  gm_bcc("nikolai.vetr@gmail.com") %>%
  gm_subject("WAI Staff + Board Paired Meet & Greet") %>%
  gm_text_body(paste0(email_df$greetings[i], message_body)) %>%
  gm_attach_file(filename = "~/data/WAI_Board_Profile_and_Purpose_Summary.pdf")
})

#send emails
for(i in 1:length(emails)){
  cat(paste0(i, " "))
  gm_send_message(emails[[i]])
}

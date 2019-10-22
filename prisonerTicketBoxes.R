survivals <- vector(length = 1e4)
for(i in 1:length(survivals)){
  if(i %% 1e3 == 0){print(i)}
  prisoners <- 1:5
  ticketBoxes <- rbind(1:100, sample(1:100, 100, F)); rownames(ticketBoxes) <- c("Box", "Ticket")
  for(prisoner in prisoners){
    survive <- T
    success <- F
    searchBox <- prisoner
    ticket <- ticketBoxes[2,searchBox]
    if(ticket == prisoner){
      success <- T
    } else {
      for(attempt in 2:80){
        searchBox <- ticket
        ticket <- ticketBoxes[2,searchBox]
        if(ticket == prisoner){
          success <- T
          break
        }
      }
    }
    if (success == F){
      survive <- F
      break
    }
  }
  survivals[i] <- survive
}
sum(survivals)/length(survivals)
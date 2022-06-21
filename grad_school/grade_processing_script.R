g <- read.csv("~/Downloads/2020-06-16T1556_Grades-ANT_151_A01_SQ_2020.csv")
maxPts <- g[1,]
g <- g[-c(1,length(g[,1])),]
mean(g$Final.Exam..495790.[!is.na(g$Final.Exam..495790.)]) / maxPts$Final.Exam..495790.
mean(g$Midterm.II.Curved..494678.[!is.na(g$Final.Exam..495790.)]) / maxPts$Midterm.II.Curved..494678.
mean(g$Midterm.1..476580.[!is.na(g$Final.Exam..495790.)]) / maxPts$Midterm.1..476580.

g$Final.Exam..495790.[is.na(g$Final.Exam..495790.)] <- 0
g$Midterm.II.Curved..494678.[is.na(g$Midterm.II.Curved..494678.) ] <- 0
g$Midterm.1..476580.[is.na(g$Midterm.1..476580.)] <- 0

g[is.na(g)] <- 0

exams <- cbind(
  g$Final.Exam..495790. / maxPts$Final.Exam..495790.,
  g$Midterm.II.Curved..494678. / maxPts$Midterm.II.Curved..494678.,
  g$Midterm.1..476580. / maxPts$Midterm.1..476580.
)

to_drop <- (sapply(1:length(g[,1]), function(s) which.min(exams[s,])))

ex <- cbind(
  g$Final.Exam..495790. / maxPts$Final.Exam..495790. * 200,
  g$Midterm.II.Curved..494678. / maxPts$Midterm.II.Curved..494678. * 150,
  g$Midterm.1..476580. / maxPts$Midterm.1..476580 * 150.
)

ex[to_drop == 1,1] <- 0
ex[to_drop == 2,2] <- 0
ex[to_drop == 3,3] <- 0
exam_tot <- apply(ex, 1, sum)
total_pts <- rep(800, length(exam_tot)); total_pts[to_drop != 1] <- 850

raw_scores <- exam_tot +
  g$Syllabus.Worksheet..463691. + 
  g$Lab.1.Quiz..468664. +
  g$Lab.2.Quiz..472023. +
    g$Lab.3..Part.1..474951. + g$Lab.3..Part.2..474952. +
  g$Lab.4.Quiz..477801. +
  g$Lab.5.Quiz..481359. +
  g$Lab.6.Quiz..483453. +
  g$Lab.7.Quiz..483460. +
  g$Lab.8.Quiz..490236. +
  g$Lab.Final..493186. / maxPts$Lab.Final..493186. * 100 +
  g$Term.Paper.Outline..482813. / maxPts$Term.Paper.Outline..482813. * 100 +
  g$Term.Paper..493498. / maxPts$Term.Paper..493498. * 200

raw_percents <- (raw_scores / total_pts)
ec_incl <- raw_percents + g$Meme.Extra.Credit..487620. / 1000 + g$Outreach.Extra.Credit..487621. / 1000
mean(raw_percents[raw_percents > 0.35])

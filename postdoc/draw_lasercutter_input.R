library(extrafont)
library(extrafontdb)

#import fonts
# font_import()
# loadfonts()
names(pdfFonts())[grepl("oboto", names(pdfFonts()))]

#read in and process names
labnames <- trimws(readLines("~/montgomery_lab_names.txt", warn = F))
labnames <- gsub(" ", "\n", labnames)

#specify dimensions of plotting (in inches)
dim_grid <- c(8,8)
dim_unit <- c(2,1)
dim_spacing <- c(1,1)
start_loc <- c(1,1)

#construct matrix of centers
centers <- as.matrix(expand.grid(1:dim_grid[1], 1:dim_grid[2]))
centers <- centers %*% diag(dim_unit + dim_spacing)
centers <- t(t(centers) + start_loc)

#specify font information
font_name <- "Source Serif Pro"
cex <- 0.75

#now do the plotting
svg(filename = "~/lab_names_grid.svg", width = 12, height = 12)
plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = F,
     xlim = c(0, max(centers[,1]) + dim_unit[1] + dim_spacing[1]),
     ylim = c(0, max(centers[,2]) + dim_unit[2] + dim_spacing[2])
)

for(i in 1:length(names)){
  text(centers[i,1], centers[i,2], labels = labnames[i], family = font_name, cex = cex, xpd = NA)  
}

dev.off()

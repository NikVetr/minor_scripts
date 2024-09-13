library(sf)

#functions
xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

#simple demonstration of drawing a perpendicular line
plot(0,0, ylim = c(-100,100))
m <- 20
abline(0,m)
xyr <- xyrat()
mv <- m*xyr^2
mvi <- -1/mv
abline(0,mvi, col = 2)

#now to build up the polylines() function

# create a horizontal variable
npts <- 200
x <- seq(0, 2 * pi, length.out = npts)

# and a vertical variable
y <- sin(x)

# Plot the base plot without lines
plot(x,y, type = "l", 
     xlim = range(x) + diff(range(x))/10 * c(-1,1), 
     ylim = range(y) + diff(range(y))/2 * c(-1,1))

# Calculate the aspect ratio
xyr <- xyrat()

# calculate slopes of lines
m <- diff(y) / diff(x)
m <- c(m, m[length(m)])

#adjust for visual distortion
mv <- m * xyr^2

#calculate perpendicular slope
mvi <- -1/mv

#plot tangent lines at different points
npts <- length(x)
pts <- sample(1:length(x), npts)

points(x[pts], y[pts], col = 2)
for(i in pts){
  abline(y[i] - m[i] * x[i], m[i], col = 2)
}

#plot perpendicular lines at different points
for(i in pts){
  abline(y[i] - -1/m[i] * x[i], -1/m[i], col = 3)
}

#plot visually perpendicular lines at different points
for(i in pts){
  abline(y[i] - mvi[i] * x[i], mvi[i], col = 4)
}

#plot visually perpendicular line *segments* at different points
seglen <- 1
for(i in pts){
  segments(x0 = x[i] - seglen / 2,
           x1 = x[i] + seglen / 2,
           y0 = y[i] - seglen / 2 * mvi[i],
           y1 = y[i] + seglen / 2 * mvi[i],
           col = 2, lwd = 2)
}

#plot visually perpendicular line *segments* with same visual lengths (on horiz axis) at different points
target_seglen <- 0.5
for(i in pts){
  endpts <- c(x0 = x[i] - target_seglen / 2, 
              x1 = x[i] + target_seglen / 2, 
              y0 = y[i] - target_seglen / 2 * mvi[i], 
              y1 = y[i] + target_seglen / 2 * mvi[i])
  actual_len <- sqrt((endpts["y1"] - endpts["y0"])^2 * xyr^2 + 
    (endpts["x1"] - endpts["x0"])^2)
  seglen <- target_seglen^2 / actual_len
  segments(x0 = x[i] - seglen / 2, 
           x1 = x[i] + seglen / 2, 
           y0 = y[i] - seglen / 2 * mvi[i], 
           y1 = y[i] + seglen / 2 * mvi[i],
           col = 3, lwd = 3)
}

#create polygon for these endpts
target_seglen <- 0.5
endpts <- cbind(x0 = x - target_seglen / 2, 
            x1 = x + target_seglen / 2, 
            y0 = y - target_seglen / 2 * mvi, 
            y1 = y + target_seglen / 2 * mvi)
actual_lens <- sqrt((endpts[,"y1"] - endpts[,"y0"])^2 * xyr^2 + 
                     (endpts[,"x1"] - endpts[,"x0"])^2)
seglens <- target_seglen^2 / actual_lens
#multiply by sign to always go in same dir from pt
poly_pts <- data.frame(x0 = x - seglens / 2 * sign(mvi),
                  x1 = x + seglens / 2 * sign(mvi), 
                  y0 = y - seglens / 2 * mvi * sign(mvi), 
                  y1 = y + seglens / 2 * mvi * sign(mvi))

polygon(x = c(poly_pts$x0, rev(poly_pts$x1)),
        y = c(poly_pts$y0, rev(poly_pts$y1)))

#create polygon with variable width for these endpts
target_seglen <- 0.1 + abs(y) / 4
endpts <- cbind(x0 = x - target_seglen / 2, 
                x1 = x + target_seglen / 2, 
                y0 = y - target_seglen / 2 * mvi, 
                y1 = y + target_seglen / 2 * mvi)
actual_lens <- sqrt((endpts[,"y1"] - endpts[,"y0"])^2 * xyr^2 + 
                      (endpts[,"x1"] - endpts[,"x0"])^2)
seglens <- target_seglen^2 / actual_lens
#multiply by sign to always go in same dir from pt
poly_pts <- data.frame(x0 = x - seglens / 2 * sign(mvi),
                       x1 = x + seglens / 2 * sign(mvi), 
                       y0 = y - seglens / 2 * mvi * sign(mvi), 
                       y1 = y + seglens / 2 * mvi * sign(mvi))

polygon(x = c(poly_pts$x0, rev(poly_pts$x1)),
        y = c(poly_pts$y0, rev(poly_pts$y1)), 
        border = 1, col = adjustcolor(1, 0.2))

#create polygon with variable width for these endpts
#but fix the problem with internal loops / bowties
target_seglen <- 0.1 + abs(y)^2*2.5
endpts <- cbind(x0 = x - target_seglen / 2, 
                x1 = x + target_seglen / 2, 
                y0 = y - target_seglen / 2 * mvi, 
                y1 = y + target_seglen / 2 * mvi)
actual_lens <- sqrt((endpts[,"y1"] - endpts[,"y0"])^2 * xyr^2 + 
                      (endpts[,"x1"] - endpts[,"x0"])^2)
seglens <- target_seglen^2 / actual_lens
#multiply by sign to always go in same dir from pt
poly_pts <- data.frame(x0 = x - seglens / 2 * sign(mvi),
                       x1 = x + seglens / 2 * sign(mvi), 
                       y0 = y - seglens / 2 * mvi * sign(mvi), 
                       y1 = y + seglens / 2 * mvi * sign(mvi))
poly_coords <- data.frame(x = c(poly_pts$x0, rev(poly_pts$x1), poly_pts$x0[1]),
                          y = c(poly_pts$y0, rev(poly_pts$y1), poly_pts$y0[1]))

#this can be a complex polygon, so we should make it simple

#can use a concave hull?
poly_coords_simple <- as.data.frame(concaveman::concaveman(as.matrix(poly_coords)))
colnames(poly_coords_simple) <- c("x", "y")

#or find the union of the component quadrilaterals 
quads_sf <- lapply(1:(npts-1), function(quad) {
  quad_pts <- poly_pts[quad:(quad+1),]
  quad_pts_big <- cbind(c(quad_pts$x0, rev(quad_pts$x1), quad_pts$x0[1]), 
                        c(quad_pts$y0, rev(quad_pts$y1), quad_pts$y0[1]))
  # if(quad %% 5 == 0){polygon(quad_pts_big)}
  
  #convert to the st format
  return(st_polygon(list(quad_pts_big))  
)
})

multipolygon <- st_sfc(quads_sf)
multipolygon <- st_make_valid(multipolygon)
multipolygon <- st_buffer(multipolygon, dist = 1E-6) #some artefacting from floating point arithmetic
simple_polygon <- st_union(multipolygon)
poly_coords_simple <- as.data.frame(as.matrix(simple_polygon[[1]]))
colnames(poly_coords_simple) <- c("x", "y")

# polygon(poly_coords$x, poly_coords$y,
#         border = 5, col = adjustcolor(1, 0))
polygon(poly_coords_simple$x, poly_coords_simple$y, 
        border = 1, col = adjustcolor(1, 0.2))


#### now put it all in a function ####
polylines <- function(x, y = NULL, lwd, col = NA, border = NULL, complex = F, ...){
  if(is.null(y)){
    y <- x[,2]
    x <- x[,1]
  }
  
  # Calculate the current aspect ratio
  xyr <- xyrat()
  
  # calculate slopes of lines
  m <- diff(y) / diff(x)
  m <- c(m, m[length(m)])
  
  #adjust for visual distortion
  mv <- m * xyr^2
  
  #calculate perpendicular slope
  mvi <- -1/mv
  
  #find length of vectors
  npts <- length(x)
  pts <- 1:length(x)
  
  #calculate endpoints of vertices
  endpts <- cbind(x0 = x - lwd / 2, 
                  x1 = x + lwd / 2, 
                  y0 = y - lwd / 2 * mvi, 
                  y1 = y + lwd / 2 * mvi)
  actual_lens <- sqrt((endpts[,"y1"] - endpts[,"y0"])^2 * xyr^2 + 
                        (endpts[,"x1"] - endpts[,"x0"])^2)
  lwds <- lwd^2 / actual_lens
  
  #multiply by sign of slope to always go in same dir from pt
  poly_pts <- data.frame(x0 = x - lwds / 2 * sign(mvi),
                         x1 = x + lwds / 2 * sign(mvi), 
                         y0 = y - lwds / 2 * mvi * sign(mvi), 
                         y1 = y + lwds / 2 * mvi * sign(mvi))
  
  #this can be a complex polygon, so we should make it simple
  #by finding the union of the component quadrilaterals 
  if(complex){
    poly_coords <- data.frame(x = c(poly_pts$x0, rev(poly_pts$x1), poly_pts$x0[1]),
                              y = c(poly_pts$y0, rev(poly_pts$y1), poly_pts$y0[1]))
    polygon(poly_coords$x, poly_coords$y, 
            border = border, col = col, ...)
    
  } else {
    
    quads_sf <- lapply(1:(npts-1), function(quad) {
      quad_pts <- poly_pts[quad:(quad+1),]
      quad_pts_big <- cbind(c(quad_pts$x0, rev(quad_pts$x1), quad_pts$x0[1]), 
                            c(quad_pts$y0, rev(quad_pts$y1), quad_pts$y0[1]))
      
      #convert to the st format
      return(st_polygon(list(quad_pts_big))  
      )
    })
    
    #unite and buffer the resulting polygons
    multipolygon <- st_sfc(quads_sf)
    multipolygon <- st_make_valid(multipolygon)
    multipolygon <- st_buffer(multipolygon, dist = 1E-6) #some artefacting from floating point arithmetic
    simple_polygon <- st_union(multipolygon)
    poly_coords_simple <- as.data.frame(as.matrix(simple_polygon[[1]]))
    colnames(poly_coords_simple) <- c("x", "y")
    
    #plot
    polygon(poly_coords_simple$x, poly_coords_simple$y, 
            border = border, col = col, ...)
    
  }
}

#### demonstration! ####

# create a horizontal variable
npts <- 200
x <- seq(0, 2 * pi, length.out = npts)

# and a vertical variable
y <- sin(x)

# Define variable line widths
lwds <- 1 + 1.75 * abs(sin(x))^2

# Plot the base plot without lines
plot(x,y, type = "l", 
     xlim = range(x) + diff(range(x))/10 * c(-1,1), 
     ylim = range(y) + diff(range(y))/2 * c(-1,1))


polylines(x, y, lwd = lwds, col = adjustcolor(1, 0.2), complex = F)

#### ggplot comparison ####
# Load required libraries
# Load required libraries
library(ggplot2)
library(dplyr)
library(svglite)

# Create a sample dataset
set.seed(123)
data <- data.frame(
  x = seq(1, 10, by = 0.1),
  y = sin(seq(1, 10, by = 0.1))
)

# Calculate thickness as a function of distance from zero
data <- data %>%
  mutate(thickness = abs(y) + 0.5) # Example function of distance from zero

# Define the color with transparency
line_color <- adjustcolor("black", alpha.f = 0.2) # Adjusts the color 'black' to be semi-transparent

# Plot the data with variable thickness
p <- ggplot(data, aes(x = x, y = y)) +
  geom_path(aes(linewidth = thickness), lineend = "round", color = line_color) +
  scale_linewidth_continuous(range = c(0.5, 3)) + # Adjust the linewidth range if needed
  theme_minimal() +
  labs(title = "Line with Variable Thickness (ggplot2)",
       x = "X-axis",
       y = "Y-axis",
       linewidth = "Thickness")

# Save the plot as an SVG file
ggsave("~/line_with_variable_thickness.svg", plot = p, device = "svg")

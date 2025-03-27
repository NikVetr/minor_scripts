library(sf)

# Outer circle (radius = 1)
n <- 100
theta <- seq(0, 2*pi, length.out=n)
outer_circle <- cbind(cos(theta), sin(theta))

# Inner circle (radius = 0.5, hole)
inner_circle <- cbind(0.5 * cos(theta), 0.5 * sin(theta))[n:1,]

#combine circles
both_circles <- rbind(outer_circle, inner_circle)

#close circles
outer_circle <- rbind(outer_circle, outer_circle[1,])
inner_circle <- rbind(inner_circle, inner_circle[1,])
both_circles <- rbind(both_circles, both_circles[1,])

# Note: reversing inner circle points ensures proper hole orientation

# Create the polygon with a hole
circle_with_hole <- st_polygon(list(outer_circle, inner_circle))
circle_with_hole <- st_polygon(list(both_circles))

# Validate and convert to sf object
circle_sf <- st_sfc(circle_with_hole) |> st_make_valid()

# Plot
plot(circle_sf, col="black", border=NA, asp=1)


#### try with two inner circles ####

library(sf)

# Outer ellipsoid
n <- 100
theta <- seq(0, 2*pi, length.out=n)
theta_rev <- rev(theta)
outer_circle <- cbind(cos(theta) * 3, sin(theta))

# inner circles, meant to represent holes in the outer ellipsoid
inner_circle_1 <- cbind(-1 + 0.5 * cos(theta), 0.5 * sin(theta))[n:1,]
inner_circle_2 <- cbind(1 + 0.5 * cos(theta), 0.5 * sin(theta))[n:1,]
inner_circle_3 <- cbind(0 + 0.25 * cos(theta), 0.5 + 0.25 * sin(theta))[n:1,]

#combine circles
all_circles <- rbind(outer_circle, inner_circle_1, inner_circle_2, inner_circle_3)

#close circles
all_circles <- rbind(all_circles, all_circles[1,])
outer_circle <- rbind(outer_circle, outer_circle[1,])
inner_circle_1 <- rbind(inner_circle_1, inner_circle_1[1,])
inner_circle_2 <- rbind(inner_circle_2, inner_circle_2[1,])
inner_circle_3 <- rbind(inner_circle_3, inner_circle_3[1,])

# Create the polygon with a hole
circle_with_holes <- st_polygon(list(all_circles))

# Validate and convert to sf object
circle_sf <- st_sfc(circle_with_holes) |> st_make_valid()

# Plot
plot(circle_sf, border=NA, asp=1, col = adjustcolor(1, 0.5))

#### alternatively  ####

# Create a polygon with holes
circle_with_holes <- st_polygon(list(outer_circle, inner_circle_1, inner_circle_2, inner_circle_3))

# Convert to sf object
circle_sf <- st_sfc(circle_with_holes) |> st_make_valid()

# Plot
plot(circle_sf, border=NA, asp=1, col=adjustcolor(1, 0.5))

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

df <- do.call(rbind, nm$bpoly_info[[1]]$bpoly)
df <- rbind(df, df[1,])

# Example of converting dataframe points into exterior + hole explicitly
exterior_coords <- as.matrix(df[1:n1, ])  # first ring (outer polygon)
hole_coords     <- as.matrix(df[(n1+1):(n1+n2), ])  # second ring (hole)

# Create an sf polygon explicitly
glyph_poly <- st_polygon(list(exterior_coords, hole_coords))
glyph_poly <- st_polygon(list(as.matrix(df)))

# Check validity and repair if needed
glyph_poly <- st_make_valid(glyph_poly)

# Make an sf object for convenience
glyph_sf <- st_sfc(glyph_poly)

plot(glyph_sf, col = "black", border = NA)

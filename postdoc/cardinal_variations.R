library(plotwidgets)
cr <- rgb2col(t(t(c(140,21,21))))
plot(1,1,col=cr,cex = 10, pch = 19)
cr_hsl <- rgb2hsl(t(t(c(140,21,21))))
cr <- hsl2col(cr_hsl)
plot(1,1,col=cr,cex = 10, pch = 19)

cg <- hsl2col(cr_hsl + t(t(c(120,0,0))))
col2rgb(cg)
cb <- hsl2col(cr_hsl + t(t(c(120,0,0))) * 2)
col2rgb(cb)
crgb <- c(cr, cg, cb)
ca <- adjustcolor(c(red = cr, blue = cb, green = cg), 0.25)


plot(1,1,col=ca[1],cex = 10, pch = 19)

blend_with_color <- function(base_color, blend_color = "white", alpha) {
  # Split the base color into red, green, and blue components
  base_red <- col2rgb(base_color)["red", ] / 255
  base_green <- col2rgb(base_color)["green", ] / 255
  base_blue <- col2rgb(base_color)["blue", ] / 255
  
  # Split the blend color into red, green, and blue components
  blend_red <- col2rgb(blend_color)["red", ] / 255
  blend_green <- col2rgb(blend_color)["green", ] / 255
  blend_blue <- col2rgb(blend_color)["blue", ] / 255
  
  # Blend each component with the corresponding component of the blend color
  blended_red <- (alpha * base_red) + ((1 - alpha) * blend_red)
  blended_green <- (alpha * base_green) + ((1 - alpha) * blend_green)
  blended_blue <- (alpha * base_blue) + ((1 - alpha) * blend_blue)
  
  # Convert the blended components back to a hex color
  blended_color <- rgb(blended_red, blended_green, blended_blue)
  return(list(rgb = rbind(blended_red, blended_green, blended_blue), hex = blended_color))
}


plot(1:3,1:3,col=blend_with_color(crgb, alpha = 0.1)$hex,cex = 10, pch = 19)
plot(1:3,1:3,col=blend_with_color(crgb, "black", alpha = 0.4)$hex,cex = 10, pch = 19)
plot(1:3,1:3,col=crgb,cex = 10, pch = 19)

round(blend_with_color(crgb, "white", alpha = 0.1)$rgb * 255)
round(blend_with_color(crgb, "black", alpha = 0.4)$rgb * 255)
blend_with_color(crgb, "black", alpha = 0.4)

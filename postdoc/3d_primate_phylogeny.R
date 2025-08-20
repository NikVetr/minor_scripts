############################################################################
## 0.  Packages & data  ────────────────────────────────────────────────────
library(ape)
library(phytools)
library(rgl)

data(primate.tree)
tr <- primate.tree
mrca_node  <- getMRCA(tr, c("Gorilla_gorilla", "Macaca_mulatta"))
tr <- extract.clade(tr, node = mrca_node)        # keep only catarrhines

ntip  <- Ntip(tr)
N     <- ntip + tr$Nnode                         # total vertices

############################################################################
## 1.  Rectangular coordinates (depth & tip order)                         #
############################################################################
depth <- numeric(N)
for (e in seq_len(nrow(tr$edge))) {              # cladewise ⇒ parent first
  p <- tr$edge[e, 1];  c <- tr$edge[e, 2]
  depth[c] <- depth[p] + tr$edge.length[e]
}

xpos <- numeric(N)
xpos[seq_len(ntip)] <- seq_len(ntip)             # even spacing for tips
for (node in rev(unique(tr$edge[, 1]))) {
  kids <- tr$edge[tr$edge[, 1] == node, 2]
  xpos[node] <- mean(xpos[kids])
}
theta <- function(x) 2 * pi * (x - 1) / ntip

############################################################################
## 2.  Build 3-D segment tables                                            #
############################################################################
## (a) parent-child radial lines
radial_segments <- lapply(seq_len(nrow(tr$edge)), function(e) {
  p <- tr$edge[e, 1];  c <- tr$edge[e, 2]
  th <- theta(xpos[c])
  matrix(c(depth[p]*cos(th), depth[p]*sin(th), depth[p],
           depth[c]*cos(th), depth[c]*sin(th), depth[c]),
         ncol = 3, byrow = TRUE)
})

## (b) constant-radius arcs joining sisters
arc_segments <- lapply(unique(tr$edge[, 1]), function(node) {
  kids <- tr$edge[tr$edge[, 1] == node, 2]
  if (length(kids) < 2) return(NULL)
  r   <- depth[node]
  th0 <- theta(min(xpos[kids]));  th1 <- theta(max(xpos[kids]))
  if (th1 < th0) th1 <- th1 + 2*pi            # wrap through 0°
  ts  <- seq(th0, th1, length.out = 80)
  cbind(r*cos(ts), r*sin(ts), r)
})
arc_segments <- arc_segments[!sapply(arc_segments, is.null)]

############################################################################
## 3.  Render with rgl                                                     #
############################################################################

open3d()
bg3d("white"); aspect3d(1, 1, 1)

## draw radial lines
for (seg in radial_segments)
  segments3d(seg, color = "black", lwd = 2)      # :contentReference[oaicite:0]{index=0}

## draw arcs
for (arc in arc_segments)
  lines3d(arc, color = "black", lwd = 2)         # :contentReference[oaicite:1]{index=1}


############################################################################
## 0.  Packages + (tree-building code stays exactly as in your script)    ##
############################################################################

# 1. Set output dirs
vid_dir   <- "~/tree_viz"             # change as needed
frame_dir <- file.path(vid_dir, "frames")
output_mp4 <- file.path(vid_dir, "catarrhine_spin.mp4")
fps <- 60

# 2. Create dir
dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Open 3D window and draw tree (assume radial_segments and arc_segments exist)
close3d()
open3d(windowRect = c(0, 0, 1000, 1000))
bg3d("white"); aspect3d(1, 1, 1)
for (seg in radial_segments) segments3d(seg, color = "black", lwd = 2)
for (arc in arc_segments)   lines3d(arc,  color = "black", lwd = 2)
view3d(theta = 0, phi = -80, fov = 0)

# 4. Save images
spin_fun <- spin3d(axis = c(0.1, 0, 1), rpm = 12)
movie3d(spin_fun,
        duration = 10,
        fps = fps,
        dir = frame_dir,
        movie = "",             # disables built-in FFmpeg
        clean = FALSE,          # keep PNGs
        startTime = 0)

# 5. Run ffmpeg manually (frame naming: 00001.png, 00002.png, ...)
system(paste0(
  "cd ", normalizePath(frame_dir), " && ",
  "ffmpeg -y -framerate ", fps,
  " -i %03d.png -vcodec libx264 -crf 20 -pix_fmt yuv420p ",
  shQuote(normalizePath(output_mp4))
))

library(parallel)
library(rgl)

# Create output folder
dir.create("~/rgl_snapshots", showWarnings = FALSE)

# Function to create a 3D plot and save snapshot
render_sphere <- function(i) {
  library(rgl)  # must re-load inside worker
  open3d(useNULL = F)  # render off-screen
  text3d(0,0,0,i)
  snapshot_filename <- sprintf("~/rgl_snapshots/sphere_%02d.png", i)
  rgl.snapshot(snapshot_filename)
  rgl.close()
  return(snapshot_filename)
}

# Parallel execution
cl <- makeCluster(4)
clusterExport(cl, varlist = c("render_sphere"))
clusterEvalQ(cl, library(rgl))
clusterExport(cl, "render_sphere")


image_files <- parLapply(cl, 1:12, render_sphere)
stopCluster(cl)

# Now load images in main session
print(image_files)

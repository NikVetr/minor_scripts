# motrpac_letters_curved_no_bite_caps.R
# run per-letter with rstudio "run code section" (cmd-shift-t)
# set write_to_svg <- TRUE and run the export section at the bottom

write_to_svg <- FALSE
out_dir <- "letters_svg"
sheet_file <- file.path(out_dir, "alphabet_all.svg")

# ---------- core primitives ----------

.lolli_r <- 0.135
.lolli_w <- 0.105

.arc_n_head <- 160
.arc_n_bite <- 110
.arc_n_stem <- 180

pt <- function(x, y) c(x, y)

.plot_canvas <- function() {
  op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
}

.arc_pts <- function(cx, cy, r, a0, a1, n) {
  aa <- seq(a0, a1, length.out = n)
  list(x = cx + r * cos(aa), y = cy + r * sin(aa), a = aa)
}

# circumcenter of 3 points
.circumcenter <- function(a, b, p) {
  ax <- a[1]; ay <- a[2]
  bx <- b[1]; by <- b[2]
  px <- p[1]; py <- p[2]
  d <- 2 * (ax * (by - py) + bx * (py - ay) + px * (ay - by))
  if (abs(d) < 1e-10) return(c(NA_real_, NA_real_))
  ax2 <- ax*ax + ay*ay
  bx2 <- bx*bx + by*by
  px2 <- px*px + py*py
  ux <- (ax2 * (by - py) + bx2 * (py - ay) + px2 * (ay - by)) / d
  uy <- (ax2 * (px - bx) + bx2 * (ax - px) + px2 * (bx - ax)) / d
  c(ux, uy)
}

.norm_ang <- function(a) (a %% (2 * pi))

.on_ccw_arc <- function(th0, th1, th_test) {
  th0 <- .norm_ang(th0); th1 <- .norm_ang(th1); th_test <- .norm_ang(th_test)
  if (th1 < th0) th1 <- th1 + 2 * pi
  if (th_test < th0) th_test <- th_test + 2 * pi
  th_test >= th0 && th_test <= th1
}

.on_directed_arc <- function(a0, a1, t) {
  if (a1 >= a0) {
    t2 <- t
    while (t2 < a0) t2 <- t2 + 2 * pi
    while (t2 > a1) t2 <- t2 - 2 * pi
    return(t2 >= a0 && t2 <= a1)
  }
  t2 <- t
  while (t2 > a0) t2 <- t2 - 2 * pi
  while (t2 < a1) t2 <- t2 + 2 * pi
  t2 <= a0 && t2 >= a1
}

.pick_head_arc <- function(aL, aR, stem_dir_ang) {
  d1 <- (aR - aL) %% (2 * pi)
  c1 <- c(aL, aL + d1)
  c2 <- c(aL, aL - (2 * pi - d1))
  len1 <- abs(c1[2] - c1[1])
  len2 <- abs(c2[2] - c2[1])
  c1_has <- .on_directed_arc(c1[1], c1[2], stem_dir_ang)
  c2_has <- .on_directed_arc(c2[1], c2[2], stem_dir_ang)
  
  if (!c1_has && len1 >= pi) return(c1)
  if (!c2_has && len2 >= pi) return(c2)
  if (len1 >= len2) return(c1)
  c2
}

.pick_bite_arc <- function(aR, aL, toward_candy_ang) {
  d <- (aL - aR) %% (2 * pi)
  c1 <- c(aR, aR + d)
  c2 <- c(aR, aR - (2 * pi - d))
  len1 <- abs(c1[2] - c1[1])
  len2 <- abs(c2[2] - c2[1])
  c1_has <- .on_directed_arc(c1[1], c1[2], toward_candy_ang)
  c2_has <- .on_directed_arc(c2[1], c2[2], toward_candy_ang)
  
  if (len1 <= pi && c1_has) return(c1)
  if (len2 <= pi && c2_has) return(c2)
  if (len1 <= len2) return(c1)
  c2
}

plot_circle <- function(center, r = .lolli_r, col = "black") {
  symbols(center[1], center[2], circles = r, inches = FALSE, add = TRUE, bg = col, fg = NA)
  invisible(NULL)
}

# plot_lollipop: candy at 'candy' end, bite at 'bite' end; optional curvature via bulge
# bulge = signed sagitta in plot coordinates: 0 straight, >0 bends left, <0 bends right
plot_lollipop <- function(candy, bite, col = "black", bulge = 0) {
  r <- .lolli_r
  w <- .lolli_w
  w2 <- w / 2
  tt <- sqrt(max(r * r - w2 * w2, 0))
  
  dx <- bite[1] - candy[1]
  dy <- bite[2] - candy[2]
  L <- sqrt(dx * dx + dy * dy)
  if (L < 1e-6) return(invisible(NULL))
  
  # straight stem
  if (abs(bulge) < 1e-8) {
    ux <- dx / L; uy <- dy / L
    vx <- -uy; vy <- ux
    
    pL <- c(candy[1] + ux * tt + vx * w2, candy[2] + uy * tt + vy * w2)
    pR <- c(candy[1] + ux * tt - vx * w2, candy[2] + uy * tt - vy * w2)
    qL <- c(bite[1] - ux * tt + vx * w2, bite[2] - uy * tt + vy * w2)
    qR <- c(bite[1] - ux * tt - vx * w2, bite[2] - uy * tt - vy * w2)
    
    a_pL <- atan2(pL[2] - candy[2], pL[1] - candy[1])
    a_pR <- atan2(pR[2] - candy[2], pR[1] - candy[1])
    stem_dir_ang <- atan2(uy, ux)
    head_ab <- .pick_head_arc(a_pL, a_pR, stem_dir_ang)
    head_pts <- .arc_pts(candy[1], candy[2], r, head_ab[1], head_ab[2], .arc_n_head)
    
    a_qR <- atan2(qR[2] - bite[2], qR[1] - bite[1])
    a_qL <- atan2(qL[2] - bite[2], qL[1] - bite[1])
    toward_candy_ang <- atan2(-uy, -ux)
    bite_ab <- .pick_bite_arc(a_qR, a_qL, toward_candy_ang)
    bite_pts <- .arc_pts(bite[1], bite[2], r, bite_ab[1], bite_ab[2], .arc_n_bite)
    
    x <- c(head_pts$x, qR[1], bite_pts$x, pL[1])
    y <- c(head_pts$y, qR[2], bite_pts$y, pL[2])
    polygon(x, y, col = col, border = NA)
    return(invisible(NULL))
  }
  
  # curved stem via arc circle through candy, bite, and a bulge point
  ux_ch <- dx / L; uy_ch <- dy / L
  nx <- -uy_ch; ny <- ux_ch
  mid <- c((candy[1] + bite[1]) / 2, (candy[2] + bite[2]) / 2)
  pbul <- c(mid[1] + bulge * nx, mid[2] + bulge * ny)
  
  cc <- .circumcenter(candy, bite, pbul)
  if (any(!is.finite(cc))) return(plot_lollipop(candy, bite, col = col, bulge = 0))
  Rc <- sqrt((candy[1] - cc[1])^2 + (candy[2] - cc[2])^2)
  if (!is.finite(Rc) || Rc < 1e-6) return(plot_lollipop(candy, bite, col = col, bulge = 0))
  
  th0 <- atan2(candy[2] - cc[2], candy[1] - cc[1])
  th1 <- atan2(bite[2] - cc[2], bite[1] - cc[1])
  thP <- atan2(pbul[2] - cc[2], pbul[1] - cc[1])
  
  dir <- if (.on_ccw_arc(th0, th1, thP)) +1 else -1
  
  th0u <- th0
  th1u <- th1
  if (dir == +1) while (th1u < th0u) th1u <- th1u + 2 * pi
  if (dir == -1) while (th1u > th0u) th1u <- th1u - 2 * pi
  
  dth <- tt / Rc
  th_start <- th0u + dir * dth
  th_end <- th1u - dir * dth
  if ((dir == +1 && th_end <= th_start) || (dir == -1 && th_end >= th_start)) {
    return(plot_lollipop(candy, bite, col = col, bulge = 0))
  }
  
  aa <- seq(th_start, th_end, length.out = .arc_n_stem)
  cx <- cc[1] + Rc * cos(aa)
  cy <- cc[2] + Rc * sin(aa)
  
  if (dir == +1) {
    tx <- -sin(aa); ty <- cos(aa)
  } else {
    tx <- sin(aa); ty <- -cos(aa)
  }
  tnorm <- sqrt(tx*tx + ty*ty)
  tx <- tx / tnorm; ty <- ty / tnorm
  vx <- -ty; vy <- tx
  
  left_x <- cx + vx * w2
  left_y <- cy + vy * w2
  right_x <- cx - vx * w2
  right_y <- cy - vy * w2
  
  uAx <- tx[1]; uAy <- ty[1]
  vAx <- vx[1]; vAy <- vy[1]
  uBx <- tx[length(tx)]; uBy <- ty[length(ty)]
  vBx <- vx[length(vx)]; vBy <- vy[length(vy)]
  
  pL <- c(candy[1] + uAx * tt + vAx * w2, candy[2] + uAy * tt + vAy * w2)
  pR <- c(candy[1] + uAx * tt - vAx * w2, candy[2] + uAy * tt - vAy * w2)
  qL <- c(bite[1] - uBx * tt + vBx * w2, bite[2] - uBy * tt + vBy * w2)
  qR <- c(bite[1] - uBx * tt - vBx * w2, bite[2] - uBy * tt - vBy * w2)
  
  a_pL <- atan2(pL[2] - candy[2], pL[1] - candy[1])
  a_pR <- atan2(pR[2] - candy[2], pR[1] - candy[1])
  stem_dir_ang <- atan2(uAy, uAx)
  head_ab <- .pick_head_arc(a_pL, a_pR, stem_dir_ang)
  head_pts <- .arc_pts(candy[1], candy[2], r, head_ab[1], head_ab[2], .arc_n_head)
  
  a_qR <- atan2(qR[2] - bite[2], qR[1] - bite[1])
  a_qL <- atan2(qL[2] - bite[2], qL[1] - bite[1])
  toward_candy_ang <- atan2(-uBy, -uBx)
  bite_ab <- .pick_bite_arc(a_qR, a_qL, toward_candy_ang)
  bite_pts <- .arc_pts(bite[1], bite[2], r, bite_ab[1], bite_ab[2], .arc_n_bite)
  
  right_boundary_x <- c(pR[1], right_x, qR[1])
  right_boundary_y <- c(pR[2], right_y, qR[2])
  left_boundary_x <- c(qL[1], rev(left_x), pL[1])
  left_boundary_y <- c(qL[2], rev(left_y), pL[2])
  
  x <- c(head_pts$x, right_boundary_x, bite_pts$x, left_boundary_x)
  y <- c(head_pts$y, right_boundary_y, bite_pts$y, left_boundary_y)
  polygon(x, y, col = col, border = NA)
  invisible(NULL)
}

# "bar" made of two lollipops meeting bite-to-bite at midpoint
# ensures each lollipop has only one candy end
plot_bar <- function(a, b, col = "black", bulge = 0) {
  m <- c((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)
  plot_lollipop(a, m, col = col, bulge = bulge)
  plot_lollipop(b, m, col = col, bulge = -bulge)
  invisible(NULL)
}

# ---------- anchors (starting guess) ----------
xL <- 0.22; xC <- 0.50; xR <- 0.78
yB <- 0.22; yM <- 0.50; yT <- 0.78
xL2 <- 0.34; xR2 <- 0.66
yB2 <- 0.34; yT2 <- 0.66

LETTER_EXPR <- list()

.register_letter <- function(letter, expr) {
  LETTER_EXPR[[letter]] <<- expr
  if (!isTRUE(write_to_svg)) {
    cat("#### ", letter, " ####\n", sep = "")
    eval(expr, envir = parent.frame())
  }
  invisible(NULL)
}

# ---- A ----
.register_letter("A", quote({
  .plot_canvas()
  top <- pt(xC, yT)
  bl <- pt(xL, yB)
  br <- pt(xR, yB)
  plot_bar(top, bl)
  plot_bar(top, br)
  plot_bar(pt(xC - 0.10, yM), pt(xC + 0.10, yM))
}))

# ---- B ----
.register_letter("B", quote({
  .plot_canvas()
  spine_t <- pt(xL, yT)
  spine_m <- pt(xL, yM)
  spine_b <- pt(xL, yB)
  plot_bar(spine_t, spine_m)
  plot_bar(spine_m, spine_b)
  plot_bar(spine_t, pt(xR2 + 0.12, yT2), bulge = 0.10)
  plot_bar(spine_m, pt(xR2 + 0.12, yT2), bulge = -0.08)
  plot_bar(spine_m, pt(xR2 + 0.12, yB2), bulge = 0.08)
  plot_bar(spine_b, pt(xR2 + 0.12, yB2), bulge = -0.10)
}))

# ---- C ----
.register_letter("C", quote({
  .plot_canvas()
  top <- pt(xR2 + 0.08, yT2)
  bot <- pt(xR2 + 0.08, yB2)
  plot_bar(top, bot, bulge = 0.22)
}))

# ---- D ----
.register_letter("D", quote({
  .plot_canvas()
  lt <- pt(xL2, yT2)
  lb <- pt(xL2, yB2)
  rt <- pt(xR2 + 0.12, yT2)
  rb <- pt(xR2 + 0.12, yB2)
  plot_bar(lt, lb)
  plot_bar(rt, rb, bulge = -0.18)
  plot_bar(lt, rt, bulge = -0.06)
  plot_bar(lb, rb, bulge = 0.06)
}))

# ---- E ----
.register_letter("E", quote({
  .plot_canvas()
  t <- pt(xL2, yT2); m <- pt(xL2, yM); b <- pt(xL2, yB2)
  plot_bar(t, m); plot_bar(m, b)
  plot_bar(t, pt(xR2 + 0.10, yT2))
  plot_bar(m, pt(xR2 + 0.10, yM))
  plot_bar(b, pt(xR2 + 0.10, yB2))
}))

# ---- F ----
.register_letter("F", quote({
  .plot_canvas()
  t <- pt(xL2, yT2); m <- pt(xL2, yM); b <- pt(xL2, yB2)
  plot_bar(t, m); plot_bar(m, b)
  plot_bar(t, pt(xR2 + 0.10, yT2))
  plot_bar(m, pt(xR2 + 0.10, yM))
}))

# ---- G ----
.register_letter("G", quote({
  .plot_canvas()
  top <- pt(xR2 + 0.08, yT2)
  bot <- pt(xR2 + 0.08, yB2)
  plot_bar(top, bot, bulge = 0.22)
  plot_bar(pt(xC + 0.06, yM), pt(xC + 0.14, yM))
}))

# ---- H ----
.register_letter("H", quote({
  .plot_canvas()
  lt <- pt(xL2, yT2); lm <- pt(xL2, yM); lb <- pt(xL2, yB2)
  rt <- pt(xR2, yT2); rm <- pt(xR2, yM); rb <- pt(xR2, yB2)
  plot_bar(lt, lb)
  plot_bar(rt, rb)
  plot_bar(lm, rm)
}))

# ---- I ----
.register_letter("I", quote({
  .plot_canvas()
  tl <- pt(xL2, yT2); tr <- pt(xR2, yT2)
  bl <- pt(xL2, yB2); br <- pt(xR2, yB2)
  tc <- pt(xC, yT2); bc <- pt(xC, yB2)
  plot_bar(tl, tr)
  plot_bar(bl, br)
  plot_bar(tc, bc)
}))

# ---- J ----
.register_letter("J", quote({
  .plot_canvas()
  tl <- pt(xC - 0.10, yT2); tr <- pt(xC + 0.10, yT2)
  mid <- pt(xC + 0.10, yM)
  bot <- pt(xC + 0.10, yB2)
  hook <- pt(xC - 0.06, yB2 - 0.06)
  plot_bar(tl, tr)
  plot_bar(tr, mid)
  plot_bar(mid, bot)
  plot_bar(bot, hook, bulge = 0.10)
}))

# ---- K ----
.register_letter("K", quote({
  .plot_canvas()
  t <- pt(xL2, yT2); m <- pt(xL2, yM); b <- pt(xL2, yB2)
  plot_bar(t, b)
  plot_bar(m, pt(xR2 + 0.08, yT2))
  plot_bar(m, pt(xR2 + 0.08, yB2))
}))

# ---- L ----
.register_letter("L", quote({
  .plot_canvas()
  t <- pt(xC - 0.08, yT2)
  b <- pt(xC - 0.08, yB2)
  r <- pt(xR2 + 0.10, yB2)
  plot_bar(t, b)
  plot_bar(b, r)
}))

# ---- M ----
.register_letter("M", quote({
  .plot_canvas()
  lb <- pt(xL, yB); lt <- pt(xL, yT)
  mid <- pt(xC, yM)
  rt <- pt(xR, yT); rb <- pt(xR, yB)
  plot_bar(lb, lt)
  plot_bar(lt, mid)
  plot_bar(mid, rt)
  plot_bar(rt, rb)
  plot_circle(pt(xR - 0.06, yB - 0.06), r = .lolli_r * 0.30)
}))

# ---- N ----
.register_letter("N", quote({
  .plot_canvas()
  lt <- pt(xL2, yT2); lb <- pt(xL2, yB2)
  rt <- pt(xR2, yT2); rb <- pt(xR2, yB2)
  plot_bar(lt, lb)
  plot_bar(rt, rb)
  plot_bar(lt, rb)
}))

# ---- O ----
.register_letter("O", quote({
  .plot_canvas()
  a <- pt(xC, yT2)
  b <- pt(xR2 + 0.10, yM)
  c <- pt(xC, yB2)
  d <- pt(xL2 - 0.10, yM)
  plot_bar(a, b, bulge = -0.14)
  plot_bar(b, c, bulge = -0.14)
  plot_bar(c, d, bulge = -0.14)
  plot_bar(d, a, bulge = -0.14)
}))

# ---- P ----
.register_letter("P", quote({
  .plot_canvas()
  stem_t <- pt(xL2, yT2)
  stem_b <- pt(xL2, yB2)
  plot_bar(stem_t, stem_b)
  
  rt <- pt(xR2 + 0.12, yT2)
  rm <- pt(xR2 + 0.12, yM)
  plot_bar(rt, rm, bulge = -0.16)
  plot_bar(stem_t, rt, bulge = -0.05)
  plot_bar(pt(xL2, yM), rm, bulge = 0.05)
}))

# ---- Q ----
.register_letter("Q", quote({
  .plot_canvas()
  a <- pt(xC, yT2)
  b <- pt(xR2 + 0.10, yM)
  c <- pt(xC, yB2)
  d <- pt(xL2 - 0.10, yM)
  plot_bar(a, b, bulge = -0.14)
  plot_bar(b, c, bulge = -0.14)
  plot_bar(c, d, bulge = -0.14)
  plot_bar(d, a, bulge = -0.14)
  plot_bar(pt(xC + 0.12, yB2 + 0.02), pt(xC + 0.22, yB2 - 0.10), bulge = -0.05)
}))

# ---- R ----
.register_letter("R", quote({
  .plot_canvas()
  stem_t <- pt(xL2, yT2)
  stem_b <- pt(xL2, yB2)
  plot_bar(stem_t, stem_b)
  
  rt <- pt(xR2 + 0.12, yT2)
  rm <- pt(xR2 + 0.12, yM)
  plot_bar(rt, rm, bulge = -0.16)
  plot_bar(stem_t, rt, bulge = -0.05)
  plot_bar(pt(xL2, yM), rm, bulge = 0.05)
  
  plot_bar(pt(xL2, yM), pt(xR2 + 0.10, yB2 - 0.02))
}))

# ---- S ----
.register_letter("S", quote({
  .plot_canvas()
  top <- pt(xC + 0.12, yT2)
  mid <- pt(xC - 0.05, yM)
  bot <- pt(xC - 0.12, yB2)
  plot_bar(top, mid, bulge = 0.14)
  plot_bar(mid, bot, bulge = -0.14)
}))

# ---- T ----
.register_letter("T", quote({
  .plot_canvas()
  tl <- pt(xL2, yT2); tr <- pt(xR2, yT2)
  plot_bar(tl, tr)
  plot_bar(pt(xC, yT2), pt(xC, yB2))
}))

# ---- U ----
.register_letter("U", quote({
  .plot_canvas()
  lt <- pt(xL2, yT2); lb <- pt(xL2, yB2)
  rt <- pt(xR2, yT2); rb <- pt(xR2, yB2)
  plot_bar(lt, lb)
  plot_bar(rt, rb)
  plot_bar(lb, rb, bulge = -0.10)
}))

# ---- V ----
.register_letter("V", quote({
  .plot_canvas()
  tl <- pt(xL2, yT2)
  tr <- pt(xR2, yT2)
  b <- pt(xC, yB2)
  plot_bar(b, tl)
  plot_bar(b, tr)
}))

# ---- W ----
.register_letter("W", quote({
  .plot_canvas()
  lt <- pt(xL, yT); lb <- pt(xL, yB)
  mid <- pt(xC, yM)
  rt <- pt(xR, yT); rb <- pt(xR, yB)
  plot_bar(lt, lb)
  plot_bar(lb, mid)
  plot_bar(mid, rb)
  plot_bar(rb, rt)
}))

# ---- X ----
.register_letter("X", quote({
  .plot_canvas()
  c <- pt(xC, yM)
  plot_bar(c, pt(xL2, yT2))
  plot_bar(c, pt(xR2, yT2))
  plot_bar(c, pt(xL2, yB2))
  plot_bar(c, pt(xR2, yB2))
}))

# ---- Y ----
.register_letter("Y", quote({
  .plot_canvas()
  c <- pt(xC, yM)
  plot_bar(c, pt(xL2, yT2))
  plot_bar(c, pt(xR2, yT2))
  plot_bar(c, pt(xC, yB2))
}))

# ---- Z ----
.register_letter("Z", quote({
  .plot_canvas()
  tl <- pt(xL2, yT2); tr <- pt(xR2, yT2)
  bl <- pt(xL2, yB2); br <- pt(xR2, yB2)
  plot_bar(tl, tr)
  plot_bar(tr, bl)
  plot_bar(bl, br)
}))

# ---- EXPORT ALL (run this section only when write_to_svg <- TRUE) ----
if (isTRUE(write_to_svg)) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (ch in names(LETTER_EXPR)) {
    f <- file.path(out_dir, paste0(ch, ".svg"))
    grDevices::svg(filename = f, width = 4.0, height = 4.0, bg = "transparent")
    .plot_canvas()
    eval(LETTER_EXPR[[ch]])
    grDevices::dev.off()
  }
  
  grDevices::svg(filename = sheet_file, width = 13 * 2.0, height = 2 * 2.0, bg = "transparent")
  op <- par(mfrow = c(2, 13), mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  for (ch in names(LETTER_EXPR)) {
    .plot_canvas()
    eval(LETTER_EXPR[[ch]])
  }
  grDevices::dev.off()
  
  message("wrote svgs to: ", normalizePath(out_dir))
  message("wrote sheet:   ", normalizePath(sheet_file))
}

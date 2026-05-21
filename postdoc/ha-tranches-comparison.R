d <- readxl::read_excel("~/Data/rprt_TrancheNumbers_ParticipantsByTranche_2025-08-25_15-27-45.xls", sheet = 1)
library(ggrepel)
library(ggtext)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(treemapify)
library(treemap)
library(forcats)
library(patchwork)
library(cowplot)
library(purrr)
library(ggfittext)
library(ggh4x)

#lookup for MoTrPAC human cohorts
rg_lookup <- tibble::tribble(
  ~randomizationGroup, ~population,  ~activity, ~modality,
  "ADUControl",        "Adult",      "SED",     "Control",
  "ADUEndur",          "Adult",      "SED",     "EE",
  "ADUResist",         "Adult",      "SED",     "RE",
  "ATHEndur",          "Adult",      "HA",      "EE",
  "ATHResist",         "Adult",      "HA",      "RE",
  "PED Control",       "Pediatric",  "SED",     "Control",
  "PED EE",            "Pediatric",  "SED",     "EE",
  "PED HA",            "Pediatric",  "HA",      "EE",
  "PED LA",            "Pediatric",  "SED",     "EE"  # LA grouped with SED
)

# clean names (Excel can add spaces)
names(d) <- stringr::str_trim(names(d))

d <- d %>%
  select(-any_of(c("population","activity","modality"))) %>%  # drop derived if present
  left_join(rg_lookup, by = "randomizationGroup")

#check tranche info
dt <- d[,grepl("tranche_[0-9]", colnames(d))]
rowSums(dt) - d$ct_tranche_all
d$ct_ppts_left

#recode age to be consistent
d$age_grp[d$age_grp == "14+"] <- "14-17"
d$ct_tranche_999 <- d$ct_ppts_left

#pivot to long
dl <- d %>%
  pivot_longer(
    cols = matches("^ct_tranche_\\d+$"),
    names_to = "tranche",
    names_pattern = "^ct_tranche_(\\d+)$",
    values_to = "n_participants"
  ) %>%
  mutate(
    tranche = as.integer(tranche)
  ) %>%
  select(
    ord, sex_grp, age_grp,
    population, activity, modality,
    tranche, n_participants
  )
colnames(dl) <- gsub("_grp", "", colnames(dl))

# aggregate by key variables
age_levels <- c("10-13","14-17","18-39","40-59","60+")
dl_plot <- dl %>%
  mutate(
    # normalize text + NAs
    age  = str_replace_all(age, "[\u2013\u2014]", "-"),  # en/em dash -> hyphen
    age  = if_else(is.na(age) | age == "", "Unknown", age),
    sex  = if_else(is.na(sex) | sex == "", "Unknown", sex),
    population = if_else(is.na(population) | population == "", "Unknown", population),
    activity   = if_else(is.na(activity)   | activity == "",   "Unknown", activity),
    modality   = if_else(is.na(modality)   | modality == "",   "Unknown", modality),
    # set ordered factor with all bins (prevents "unknown level" warnings)
    age = factor(age, levels = age_levels, ordered = TRUE),
    
    tranche_group = case_when(
      tranche == 5     ~ "Tranche 5",
      tranche %in% 0:4 ~ "Tranches 0–4",
      tranche == 999   ~ "Left to process",
      TRUE             ~ "Other"
    ),
    tranche_group = factor(tranche_group, levels = c("Tranches 0–4","Tranche 5"))
  ) %>%
  filter(tranche_group %in% c("Tranches 0–4","Tranche 5")) %>%
  group_by(population, activity, modality, age, tranche_group) %>%
  summarise(n = sum(n_participants, na.rm = TRUE), .groups = "drop")

#and plot
#### same size facets ####
png("~/Pictures/MoTrPAC_HA/tranche_sample-size_comparison.png",
    width = 2400, height = 2000, res = 250, pointsize = 10)

ggplot(
  dl_plot,
  aes(area = n, fill = tranche_group, alpha = age, label = paste0(age, "\n(N=", n, ")"))
) +
  geom_treemap() +
  geom_treemap_text(
    colour = "white",
    place = "centre",
    grow = FALSE,        # fixed size across all tiles
    reflow = TRUE,
    min.size = 0         # raise/lower this to taste
  ) +
  # facet_grid(population + activity ~ modality) +
  ggh4x::facet_nested(
    rows  = vars(population, activity),
    cols  = vars(modality),
    switch = "y"      # row strips on the right
    # scales, space etc. can stay default
  ) +
  scale_fill_manual(values = c("Tranches 0–4" = "#5DA5DA", "Tranche 5" = "#F15854")) +
  scale_alpha_manual(
    values = c("10-13" = 0.6, "14-17" = 0.7, "18-39" = 0.8, "40-59" = 0.9, "60+" = 1),
    drop = FALSE
  ) +
  guides(alpha = guide_legend(order = 2), fill = guide_legend(order = 1)) +
  labs(
    title = "Stratifying N across Tranche 5 vs Tranches 0–4",
    subtitle = "Area = participants; color = tranche; shade = age",
    fill = "Tranche group", alpha = "Age group",
    caption = "(areas comparable within-panel only, not between panels)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

dev.off()

#### diff size facets ####

# clean dl_plot2 for this branch
age_levels <- c("10-13","14-17","18-39","40-59","60+")
dl_plot2 <- dl %>%
  mutate(
    age = str_trim(str_replace_all(age, "[\u2013\u2014]", "-")),
    tranche_group = case_when(
      tranche == 5     ~ "Tranche 5",
      tranche %in% 0:4 ~ "Tranches 0–4",
      tranche == 999   ~ "Left to process",
      TRUE             ~ "Other"
    )
  ) %>%
  filter(tranche_group %in% c("Tranches 0–4","Tranche 5"),
         age %in% age_levels) %>%
  mutate(
    age           = factor(age, levels = age_levels, ordered = TRUE),
    tranche_group = factor(tranche_group, levels = c("Tranches 0–4","Tranche 5")),
    population    = factor(population, levels = c("Adult","Pediatric")),
    activity      = factor(activity,   levels = c("SED","HA")),
    modality      = factor(modality,   levels = c("Control","EE","RE"))
  ) %>%
  group_by(population, activity, modality, age, tranche_group) %>%
  summarise(n = sum(n_participants, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    facet_key = paste(population, activity, modality, sep = "|"),
    tile_id   = paste(age, tranche_group, sep = "|") # unique within facet
  )

# facet totals and scale factors
facet_totals <- dl_plot2 %>%
  group_by(facet_key) %>%
  summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop")
max_total <- max(facet_totals$total_n, na.rm = TRUE)
s_map <- setNames(sqrt(facet_totals$total_n / max_total), facet_totals$facet_key)

# layout per facet using treemap; returns 0..1 coords inside the facet
layout_one <- function(df) {
  if (nrow(df) == 0) return(NULL)
  # treemap wants a plain data.frame
  tmp <- as.data.frame(df[, c("tile_id","n")])
  tm  <- treemap::treemap(tmp, index = "tile_id", vSize = "n",
                          algorithm = "squarified", draw = FALSE)
  # tm$tm has x0,y0,w,h plus original columns
  coords <- tm$tm[, c("tile_id","x0","y0","w","h")]
  # join back to keep metadata
  out <- dplyr::left_join(df, coords, by = "tile_id") %>%
    mutate(
      xmin = x0, xmax = x0 + w,
      ymin = y0, ymax = y0 + h
    )
  out
}

# build rectangles for all facets and apply the scale factor s
lst <- split(dl_plot2, dl_plot2$facet_key)

dl_rects <- dplyr::bind_rows(lapply(names(lst), function(key) {
  df <- lst[[key]]
  s  <- s_map[[key]]
  lay <- layout_one(df)
  if (is.null(lay)) return(NULL)
  lay %>%
    mutate(
      xmin_s = 0.5 + (xmin - 0.5) * s,
      xmax_s = 0.5 + (xmax - 0.5) * s,
      ymin_s = 0.5 + (ymin - 0.5) * s,
      ymax_s = 0.5 + (ymax - 0.5) * s,
      lab    = paste0(as.character(age), "\n(N=", n, ")")
    ) %>%
    select(population, activity, modality, age, tranche_group, n,
           xmin_s, xmax_s, ymin_s, ymax_s, lab)
}))

# positions for "Total N" label at top of scaled content in each facet
content_tops <- dl_rects %>%
  group_by(population, activity, modality) %>%
  summarise(
    x_mid   = (min(xmin_s, na.rm = TRUE) + max(xmax_s, na.rm = TRUE)) / 2,
    y_top   = max(ymax_s, na.rm = TRUE),
    total_n = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_lab = y_top + 0.09)

# 6) plot
ragg::agg_png("~/Pictures/MoTrPAC_HA/tranche_sample-size_comparison_scaled.png",
              width = 3200, height = 3000, res = 300, pointsize = 11)

ggplot(dl_rects) +
  geom_rect(
    aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s,
        fill = tranche_group, alpha = age),
    colour = "white", size = 0.3
  ) +
  ggfittext::geom_fit_text(
    aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s, label = lab),
    reflow   = TRUE,
    grow     = TRUE,
    min.size = 0,
    contrast = TRUE
  ) +
  geom_text(
    data = content_tops,
    aes(x = x_mid, y = y_lab, label = paste0("Total N = ", total_n)),
    inherit.aes = FALSE,
    vjust = 1,              # anchor to top edge
    fontface = "bold",
    size = 4
  ) +
  # facet_grid(population + activity ~ modality) +
  ggh4x::facet_nested(
    rows  = vars(population, activity),
    cols  = vars(modality),
    switch = "y"      # row strips on the right
    # scales, space etc. can stay default
  ) +
  scale_fill_manual(values = c("Tranches 0–4" = "#5DA5DA", "Tranche 5" = "#F15854"), drop = FALSE) +
  scale_alpha_manual(values = c("10-13"=0.6,"14-17"=0.7,"18-39"=0.8,"40-59"=0.9,"60+"=1), drop = FALSE) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  guides(
    fill  = guide_legend(order = 1, nrow = 1, byrow = TRUE, title.position = "top", label.position = "right"),
    alpha = guide_legend(order = 2, nrow = 1, byrow = TRUE, title.position = "top", label.position = "right")
  ) +
  theme_void(base_size = 16) +  # nukes axes; bumps base size a bit
  theme(
    # facet strip labels (columns: Control/EE/RE; rows: Adult/Pediatric and SED/HA)
    strip.text.x = element_text(face = "bold", size = 18, margin = margin(t = 6, b = 6)),
    strip.text.y = element_text(face = "bold", size = 18, margin = margin(l = 6, r = 6)),
    
    # title/subtitle spacing
    plot.title = element_text(size = 22, face = "bold", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 16, margin = margin(b = 10)),
    
    # legend sizing/spacing (less cramped)
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(16, "pt"),
    legend.spacing.x = unit(60, "pt"),
    legend.box.spacing = unit(8, "pt"),
    
    # breathing room around the whole figure and between panels
    plot.margin = margin(12, 18, 12, 18),
    panel.spacing = unit(8, "pt")
  )+
  labs(
    title = "Stratifying N across Tranche 5 vs Tranches 0–4",
    subtitle = "Area = participants; color = tranche; shade = age",
    fill = "Tranche group", alpha = "Age group"
  )

dev.off()

#### split by sex ####

age_levels <- c("10-13","14-17","18-39","40-59","60+")
sex_levels <- c("Female","Male")  # add "Unknown" if you actually have it

# counts by facet × sex × age × tranche
dl_plot3 <- dl %>%
  mutate(
    age = str_trim(str_replace_all(age, "[\u2013\u2014]", "-")),
    tranche_group = case_when(
      tranche == 5     ~ "Tranche 5",
      tranche %in% 0:4 ~ "Tranches 0–4",
      tranche == 999   ~ "Left to process",
      TRUE             ~ "Other"
    ),
    sex = factor(sex, levels = sex_levels)
  ) %>%
  filter(tranche_group %in% c("Tranches 0–4","Tranche 5"),
         age %in% age_levels,
         !is.na(sex)) %>%
  mutate(
    age           = factor(age, levels = age_levels, ordered = TRUE),
    tranche_group = factor(tranche_group, levels = c("Tranches 0–4","Tranche 5")),
    population    = factor(population, levels = c("Adult","Pediatric")),
    activity      = factor(activity,   levels = c("SED","HA")),
    modality      = factor(modality,   levels = c("Control","EE","RE")),
    facet_key     = paste(population, activity, modality, sep = "|")
  ) %>%
  group_by(population, activity, modality, sex, age, tranche_group, facet_key) %>%
  summarise(n = sum(n_participants, na.rm = TRUE), .groups = "drop")

# facet totals + scale (same as before)
facet_totals <- dl_plot3 %>%
  group_by(facet_key) %>%
  summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop")
max_total <- max(facet_totals$total_n, na.rm = TRUE)
s_map <- setNames(sqrt(facet_totals$total_n / max_total), facet_totals$facet_key)

# layout: hier index = sex -> age -> tranche_group
layout_one <- function(df, prefer_horizontal = TRUE) {
  # guard
  if (nrow(df) == 0) return(list(leaf = NULL, sex = NULL))

  # build minimal input for treemap + a deterministic ordering key
  tmp <- as.data.frame(df[, c("sex","age","tranche_group","n")])

  # make sure factors are ordered as you want to see them
  tmp$age <- factor(tmp$age, levels = age_levels, ordered = TRUE)
  tmp$tranche_group <- factor(tmp$tranche_group, levels = c("Tranches 0–4","Tranche 5"))

  # numeric order key: tranche first (to cluster colors), then age within tranche
  tmp$ord <- as.integer(tmp$tranche_group) * 10L + as.integer(tmp$age)

  # optional: a stable row order helps remove run-to-run jitter when sizes tie
  tmp <- tmp[order(tmp$sex, tmp$tranche_group, tmp$age, -tmp$n), ]

  # prefer 'pivotSize' (ordered treemap) and make the canvas wider via aspRatio
  tm <- treemap::treemap(
    tmp,
    index     = c("sex","tranche_group","age"),  # <-- tranche above age = color blocks
    vSize     = "n",
    algorithm = "pivotSize",                      # ordered + stable
    sortID    = "ord",                            # honors our tranche->age order
    aspRatio  = if (prefer_horizontal) 16/9 else NA, # bias toward horizontal rows
    draw      = FALSE
  )

  # rectangles by level (sex=1, tranche_group=2, age=3)
  rects <- tm$tm

  # leaves (age-level tiles)
  leaf <- rects[rects$level == 3, c("sex","age","tranche_group","x0","y0","w","h")]
  leaf <- df %>%
    dplyr::left_join(leaf, by = c("sex","age","tranche_group")) %>%
    dplyr::mutate(xmin = x0, xmax = x0 + w, ymin = y0, ymax = y0 + h)

  # sex blocks for outlines
  sgrp <- rects[rects$level == 1, c("sex","x0","y0","w","h")] %>%
    dplyr::mutate(xmin = x0, xmax = x0 + w, ymin = y0, ymax = y0 + h) %>%
    dplyr::select(sex, xmin, xmax, ymin, ymax)

  # basic debugging prints so you can verify behavior quickly
  message(sprintf(
    "layout_one(): alg=pivotSize, aspRatio=%s; n_leaf=%d",
    if (prefer_horizontal) "16/9" else "device", nrow(leaf)
  ))

  list(leaf = leaf, sex = sgrp)
}


# build + scale per facet
lst <- split(dl_plot3, dl_plot3$facet_key)

built <- lapply(names(lst), function(key) {
  df <- lst[[key]]
  s  <- s_map[[key]]
  lay <- layout_one(df)
  if (is.null(lay$leaf)) return(NULL)
  
  leaf <- lay$leaf %>%
    mutate(
      xmin_s = 0.5 + (xmin - 0.5) * s,
      xmax_s = 0.5 + (xmax - 0.5) * s,
      ymin_s = 0.5 + (ymin - 0.5) * s,
      ymax_s = 0.5 + (ymax - 0.5) * s,
      lab    = paste0(as.character(age), "\n(N=", n, ")"),
      facet_key = key
    ) %>%
    select(population, activity, modality, sex, age, tranche_group, n,
           xmin_s, xmax_s, ymin_s, ymax_s, lab)
  
  sgrp <- lay$sex %>%
    mutate(
      xmin_s = 0.5 + (xmin - 0.5) * s,
      xmax_s = 0.5 + (xmax - 0.5) * s,
      ymin_s = 0.5 + (ymin - 0.5) * s,
      ymax_s = 0.5 + (ymax - 0.5) * s,
      facet_key = key
    ) %>%
    separate_wider_delim(facet_key, delim = "|",
                         names = c("population","activity","modality"),
                         too_few = "align_start")  # reconstruct facet cols
  
  list(leaf = leaf, sex = sgrp)
})

# combine
dl_leaf <- bind_rows(lapply(built, `[[`, "leaf"))
dl_sgrp <- bind_rows(lapply(built, `[[`, "sex"))

gap <- 0.03  # width of the gutter (panel is 0..1 wide)
# determine left vs right sex block per facet, and compute new x-bounds
s_bounds <- dl_sgrp %>%
  dplyr::group_by(population, activity, modality) %>%
  dplyr::arrange(xmin_s, .by_group = TRUE) %>%
  dplyr::mutate(side = if_else(dplyr::row_number() == 1, "left", "right")) %>%
  dplyr::mutate(
    s_xmin_old = xmin_s,
    s_xmax_old = xmax_s,
    s_ymin     = ymin_s,
    s_ymax     = ymax_s,
    s_xmin_new = if_else(side == "left",  s_xmin_old, s_xmin_old + gap/2),
    s_xmax_new = if_else(side == "left",  s_xmax_old - gap/2, s_xmax_old)
  ) %>%
  dplyr::select(population, activity, modality, sex,
                s_xmin_old, s_xmax_old, s_ymin, s_ymax, s_xmin_new, s_xmax_new)

# nudged sex block outlines (use these for the bold borders)
dl_sgrp_n <- s_bounds %>%
  dplyr::transmute(population, activity, modality, sex,
                   xmin_s = s_xmin_new, xmax_s = s_xmax_new,
                   ymin_s = s_ymin,     ymax_s = s_ymax)

# remap leaf rectangles within each sex block to the nudged bounds
dl_leaf_n <- dl_leaf %>%
  dplyr::left_join(s_bounds,
                   by = c("population","activity","modality","sex")) %>%
  dplyr::mutate(
    w_old = pmax(1e-9, s_xmax_old - s_xmin_old),
    scale = (s_xmax_new - s_xmin_new) / w_old,
    xmin_s = s_xmin_new + (xmin_s - s_xmin_old) * scale,
    xmax_s = s_xmin_new + (xmax_s - s_xmin_old) * scale
  ) %>%
  dplyr::select(-w_old, -scale)

# mark tiny tiles and precompute centers
dl_leaf_n <- dl_leaf_n %>%
  mutate(
    w    = pmax(1e-9, xmax_s - xmin_s),
    h    = pmax(1e-9, ymax_s - ymin_s),
    area = w * h,
    xc   = (xmin_s + xmax_s)/2,
    yc   = (ymin_s + ymax_s)/2,
    # heuristic thresholds; tweak to taste
    tiny = area < 0.015 | w < 0.07 | h < 0.07,
    lab_out = paste0(as.character(age), " (N=", n, ")")
  )

## side labels computed from *nudged* blocks
pad_mf  <- 0.04   # how far labels stick out (panel is 0..1 wide)
xlim_pad <- pad_mf + 0.01  # extra canvas so labels aren't clipped

dl_sgrp_lab <- dl_sgrp_n %>%
  dplyr::group_by(population, activity, modality) %>%
  dplyr::arrange(xmin_s, .by_group = TRUE) %>%
  dplyr::mutate(
    side = c("left", "right"),
    x_lab = if_else(side == "left",  xmin_s - pad_mf, xmax_s + pad_mf),
    y_lab = (ymin_s + ymax_s) / 2,
    ang   = if_else(side == "left",  90, 270)
  ) %>%
  dplyr::ungroup()

# colors to match your plot
col_t04 <- "#5DA5DA"        # Tranches 0–4
col_t5  <- "#F15854"        # Tranche 5
col_sex <- c(Female = "forestgreen", Male = "mediumpurple3")

# totals per facet × sex for each tranche group
sex_totals <- dl_plot3 %>%
  dplyr::group_by(population, activity, modality, sex) %>%
  dplyr::summarise(
    n_0_4 = sum(dplyr::if_else(tranche_group == "Tranches 0–4", n, 0), na.rm = TRUE),
    n_5   = sum(dplyr::if_else(tranche_group == "Tranche 5",     n, 0), na.rm = TRUE),
    .groups = "drop"
  )

# positions + pretty HTML label for side text
sex_labels <- dl_sgrp_lab %>%
  dplyr::left_join(sex_totals,
                   by = c("population","activity","modality","sex")) %>%
  dplyr::mutate(
    sex_col = unname(col_sex[as.character(sex)]),
    label_html = paste0(
      "<span style='color:", sex_col, "; font-weight:700'>", sex, "</span> ",
      "<span style='color:black'>(</span>",
      "<span style='color:", col_t04, "'>", n_0_4, "</span>",
      "<span style='color:black'> + </span>",
      "<span style='color:", col_t5,  "'>", n_5,   "</span>",
      "<span style='color:black'>)</span>"
    )
  )

# top-of-content label per facet (from leaves)
# --- breakdown by tranche within each facet (0–4 vs 5) ---

facet_breakdown <- dl_plot3 %>%
  dplyr::group_by(population, activity, modality, tranche_group) %>%
  dplyr::summarise(n_tranche = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = tranche_group,
    values_from = n_tranche,
    values_fill = 0
  ) %>%
  dplyr::rename(n_04 = `Tranches 0–4`, n_5 = `Tranche 5`) %>%
  dplyr::mutate(total_n = n_04 + n_5)

# --- label anchor positions from your scaled rectangles ---
content_pos <- dl_leaf %>%
  dplyr::group_by(population, activity, modality) %>%
  dplyr::summarise(
    x_mid = (min(xmin_s, na.rm = TRUE) + max(xmax_s, na.rm = TRUE)) / 2,
    y_top = max(ymax_s, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(y_lab = pmin(1.18, y_top + 0.10))  # small headroom

# --- join & build rich-text label with tranche-colored numbers ---
content_tops_colored <- content_pos %>%
  dplyr::left_join(facet_breakdown, by = c("population","activity","modality")) %>%
  dplyr::mutate(
    label_md = paste0(
      "Total N = ",
      "<span style='color:#5DA5DA;'>", n_04, "</span>",
      " + ",
      "<span style='color:#F15854;'>", n_5,  "</span>",
      " = ", total_n
    )
  )


ragg::agg_png("~/Pictures/MoTrPAC_HA/tranche_sample-size_comparison_scaled_xSex.png",
              width = 3200, height = 3000, res = 300, pointsize = 11)

ggplot() +
  # leaf rectangles (nudged)
  geom_rect(
    data = dl_leaf_n,
    aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s,
        fill = tranche_group, alpha = age),
    colour = "white", size = 0.25
  ) +
  ggfittext::geom_fit_text(
    data = dl_leaf_n,
    aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s, label = lab),
    reflow = TRUE, grow = TRUE, min.size = 0, contrast = TRUE
  ) +
  # bold sex-block borders (nudged) — color by sex if you want
  geom_rect(
    data = dl_sgrp_n,
    aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s, colour = sex),
    inherit.aes = FALSE, fill = NA, linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c(Female = "forestgreen", Male = "mediumpurple3"),   # pick your own hexes
    guide  = "none"                                      # or name = "Sex (border)" to show a legend
  ) +
  # optional sex labels at centers
  ggtext::geom_richtext(
    data = sex_labels,
    aes(x = x_lab, y = y_lab, label = label_html, angle = ang),
    inherit.aes   = FALSE,
    hjust = 0.5, vjust = 0.5,
    size  = 2.5,
    fontface = "bold",
    # make it look like plain text (no box)
    label.padding = grid::unit(0.5, "pt"),
    label.r       = grid::unit(0, "pt"),
    label.size    = 0,
    fill          = NA
  ) +
  # facetting with merged population label on the right
  ggh4x::facet_nested(
    rows  = vars(population, activity),
    cols  = vars(modality),
    switch = "y"
  ) +
  scale_fill_manual(values = c("Tranches 0–4" = "#5DA5DA", "Tranche 5" = "#F15854"),
                    drop = FALSE, name = "Tranche group") +
  scale_alpha_manual(values = c("10-13"=0.6,"14-17"=0.7,"18-39"=0.8,"40-59"=0.9,"60+"=1),
                     drop = FALSE, name = "Age group") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-0.05,1.05), ylim = c(0,1.2), clip = "off") +
  theme_void(base_size = 16) +
  theme(
    legend.position = c(0.8, 0.05),
    legend.justification = c(0, 0),
    legend.spacing.x = unit(8, "pt"),
    legend.box.spacing = unit(12, "pt"),
    strip.text.x      = element_text(face = "bold", size = 16, margin = margin(t = 6, b = 6)),
    strip.text.y.left = element_text(face = "bold", size = 16, margin = margin(l = 6, r = 6)),
    strip.placement = "outside",
    plot.title.position = "plot",  # anchor title to the plot area, not the panel
    plot.title = element_text(size = 22, face = "bold", margin = margin(b = 6), hjust = 0),
    plot.subtitle = element_text(size = 16, margin = margin(b = 10), hjust = 0),
    plot.margin = margin(12, 18, 12, 18),
    panel.spacing = unit(8, "pt")
  ) +
  labs(
    title = "Stratifying N across Tranche 5 vs Tranches 0–4",
    subtitle = "Area = # participants; color = tranche; shade = age",
    fill = "Tranche group", alpha = "Age group"
  ) +
  # facet-level totals anchored to the top of content
  ggtext::geom_richtext(
    data = content_tops_colored,
    aes(x = x_mid, y = y_lab, label = label_md),
    inherit.aes   = FALSE,
    vjust         = 1,
    size          = 4,
    fontface      = "bold",
    fill          = NA,          # no white box
    label.color   = NA,          # no border
    label.padding = grid::unit(rep(0, 4), "pt")
  )

dev.off()


# make two additional scaled xSex figures with the same styling/annotations as your existing output

plot_scaled_xSex_compare <- function(dl,
                                     tranche_a, tranche_b,
                                     label_a, label_b,
                                     out_png, title_txt) {
  age_levels <- c("10-13","14-17","18-39","40-59","60+")
  sex_levels <- c("Female","Male")
  
  # colors: keep your existing palette mapping (blue vs red)
  col_a <- "#5DA5DA"
  col_b <- "#F15854"
  col_sex <- c(Female = "forestgreen", Male = "mediumpurple3")
  
  # counts by facet × sex × age × tranche_group
  dl_plot3 <- dl %>%
    dplyr::mutate(
      age = stringr::str_trim(stringr::str_replace_all(age, "[\u2013\u2014]", "-")),
      tranche_group = dplyr::case_when(
        tranche %in% tranche_a ~ label_a,
        tranche %in% tranche_b ~ label_b,
        TRUE ~ NA_character_
      ),
      sex = factor(sex, levels = sex_levels)
    ) %>%
    dplyr::filter(
      !is.na(tranche_group),
      age %in% age_levels,
      !is.na(sex)
    ) %>%
    dplyr::mutate(
      age           = factor(age, levels = age_levels, ordered = TRUE),
      tranche_group = factor(tranche_group, levels = c(label_a, label_b)),
      population    = factor(population, levels = c("Adult","Pediatric")),
      activity      = factor(activity,   levels = c("SED","HA")),
      modality      = factor(modality,   levels = c("Control","EE","RE")),
      facet_key     = paste(population, activity, modality, sep = "|")
    ) %>%
    dplyr::group_by(population, activity, modality, sex, age, tranche_group, facet_key) %>%
    dplyr::summarise(n = sum(n_participants, na.rm = TRUE), .groups = "drop")
  
  # facet totals + scale
  facet_totals <- dl_plot3 %>%
    dplyr::group_by(facet_key) %>%
    dplyr::summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop")
  max_total <- max(facet_totals$total_n, na.rm = TRUE)
  s_map <- stats::setNames(sqrt(facet_totals$total_n / max_total), facet_totals$facet_key)
  
  # layout: sex -> tranche_group -> age (same as your script)
  layout_one <- function(df, prefer_horizontal = TRUE) {
    if (nrow(df) == 0) return(list(leaf = NULL, sex = NULL))
    
    tmp <- as.data.frame(df[, c("sex","age","tranche_group","n")])
    tmp$age <- factor(tmp$age, levels = age_levels, ordered = TRUE)
    tmp$tranche_group <- factor(tmp$tranche_group, levels = c(label_a, label_b))
    
    tmp$ord <- as.integer(tmp$tranche_group) * 10L + as.integer(tmp$age)
    tmp <- tmp[order(tmp$sex, tmp$tranche_group, tmp$age, -tmp$n), ]
    
    tm <- treemap::treemap(
      tmp,
      index     = c("sex","tranche_group","age"),
      vSize     = "n",
      algorithm = "pivotSize",
      sortID    = "ord",
      aspRatio  = if (prefer_horizontal) 16/9 else NA,
      draw      = FALSE
    )
    
    rects <- tm$tm
    
    leaf <- rects[rects$level == 3, c("sex","age","tranche_group","x0","y0","w","h")]
    leaf <- df %>%
      dplyr::left_join(leaf, by = c("sex","age","tranche_group")) %>%
      dplyr::mutate(xmin = x0, xmax = x0 + w, ymin = y0, ymax = y0 + h)
    
    sgrp <- rects[rects$level == 1, c("sex","x0","y0","w","h")] %>%
      dplyr::mutate(xmin = x0, xmax = x0 + w, ymin = y0, ymax = y0 + h) %>%
      dplyr::select(sex, xmin, xmax, ymin, ymax)
    
    list(leaf = leaf, sex = sgrp)
  }
  
  # build + scale per facet
  lst <- split(dl_plot3, dl_plot3$facet_key)
  built <- lapply(names(lst), function(key) {
    df <- lst[[key]]
    s  <- s_map[[key]]
    lay <- layout_one(df)
    if (is.null(lay$leaf)) return(NULL)
    
    leaf <- lay$leaf %>%
      dplyr::mutate(
        xmin_s = 0.5 + (xmin - 0.5) * s,
        xmax_s = 0.5 + (xmax - 0.5) * s,
        ymin_s = 0.5 + (ymin - 0.5) * s,
        ymax_s = 0.5 + (ymax - 0.5) * s,
        lab    = paste0(as.character(age), "\n(N=", n, ")"),
        facet_key = key
      ) %>%
      dplyr::select(population, activity, modality, sex, age, tranche_group, n,
                    xmin_s, xmax_s, ymin_s, ymax_s, lab)
    
    sgrp <- lay$sex %>%
      dplyr::mutate(
        xmin_s = 0.5 + (xmin - 0.5) * s,
        xmax_s = 0.5 + (xmax - 0.5) * s,
        ymin_s = 0.5 + (ymin - 0.5) * s,
        ymax_s = 0.5 + (ymax - 0.5) * s,
        facet_key = key
      ) %>%
      tidyr::separate_wider_delim(
        facet_key, delim = "|",
        names = c("population","activity","modality"),
        too_few = "align_start"
      )
    
    list(leaf = leaf, sex = sgrp)
  })
  
  dl_leaf <- dplyr::bind_rows(lapply(built, `[[`, "leaf"))
  dl_sgrp <- dplyr::bind_rows(lapply(built, `[[`, "sex"))
  
  # gutter between sex blocks + remap leaves into nudged blocks (same as your script)
  gap <- 0.03
  s_bounds <- dl_sgrp %>%
    dplyr::group_by(population, activity, modality) %>%
    dplyr::arrange(xmin_s, .by_group = TRUE) %>%
    dplyr::mutate(side = dplyr::if_else(dplyr::row_number() == 1, "left", "right")) %>%
    dplyr::mutate(
      s_xmin_old = xmin_s,
      s_xmax_old = xmax_s,
      s_ymin     = ymin_s,
      s_ymax     = ymax_s,
      s_xmin_new = dplyr::if_else(side == "left",  s_xmin_old, s_xmin_old + gap/2),
      s_xmax_new = dplyr::if_else(side == "left",  s_xmax_old - gap/2, s_xmax_old)
    ) %>%
    dplyr::select(population, activity, modality, sex,
                  s_xmin_old, s_xmax_old, s_ymin, s_ymax, s_xmin_new, s_xmax_new)
  
  dl_sgrp_n <- s_bounds %>%
    dplyr::transmute(population, activity, modality, sex,
                     xmin_s = s_xmin_new, xmax_s = s_xmax_new,
                     ymin_s = s_ymin,     ymax_s = s_ymax)
  
  dl_leaf_n <- dl_leaf %>%
    dplyr::left_join(s_bounds, by = c("population","activity","modality","sex")) %>%
    dplyr::mutate(
      w_old = pmax(1e-9, s_xmax_old - s_xmin_old),
      scale = (s_xmax_new - s_xmin_new) / w_old,
      xmin_s = s_xmin_new + (xmin_s - s_xmin_old) * scale,
      xmax_s = s_xmin_new + (xmax_s - s_xmin_old) * scale
    ) %>%
    dplyr::select(-w_old, -scale)
  
  # side label anchors (female left / male right)
  pad_mf  <- 0.04
  xlim_pad <- pad_mf + 0.01
  
  dl_sgrp_lab <- dl_sgrp_n %>%
    dplyr::group_by(population, activity, modality) %>%
    dplyr::arrange(xmin_s, .by_group = TRUE) %>%
    dplyr::mutate(
      side = c("left", "right"),
      x_lab = dplyr::if_else(side == "left",  xmin_s - pad_mf, xmax_s + pad_mf),
      y_lab = (ymin_s + ymax_s) / 2,
      ang   = dplyr::if_else(side == "left",  90, 270)
    ) %>%
    dplyr::ungroup()
  
  # totals per facet × sex (a vs b) for side labels
  sex_totals <- dl_plot3 %>%
    dplyr::group_by(population, activity, modality, sex) %>%
    dplyr::summarise(
      n_a = sum(dplyr::if_else(tranche_group == label_a, n, 0), na.rm = TRUE),
      n_b = sum(dplyr::if_else(tranche_group == label_b, n, 0), na.rm = TRUE),
      .groups = "drop"
    )
  
  sex_labels <- dl_sgrp_lab %>%
    dplyr::left_join(sex_totals, by = c("population","activity","modality","sex")) %>%
    dplyr::mutate(
      sex_col = unname(col_sex[as.character(sex)]),
      label_html = paste0(
        "<span style='color:", sex_col, "; font-weight:700'>", sex, "</span> ",
        "<span style='color:black'>(</span>",
        "<span style='color:", col_a, "'>", n_a, "</span>",
        "<span style='color:black'> + </span>",
        "<span style='color:", col_b, "'>", n_b, "</span>",
        "<span style='color:black'>)</span>"
      )
    )
  
  # top-of-content label per facet: Total N = a + b = total (colored)
  facet_breakdown <- dl_plot3 %>%
    dplyr::group_by(population, activity, modality, tranche_group) %>%
    dplyr::summarise(n_tranche = sum(n), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = tranche_group,
      values_from = n_tranche,
      values_fill = 0
    )
  
  facet_breakdown <- facet_breakdown %>%
    dplyr::mutate(
      n_a = .data[[label_a]],
      n_b = .data[[label_b]],
      total_n = n_a + n_b
    )
  
  content_pos <- dl_leaf %>%
    dplyr::group_by(population, activity, modality) %>%
    dplyr::summarise(
      x_mid = (min(xmin_s, na.rm = TRUE) + max(xmax_s, na.rm = TRUE)) / 2,
      y_top = max(ymax_s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(y_lab = pmin(1.18, y_top + 0.10))
  
  content_tops_colored <- content_pos %>%
    dplyr::left_join(facet_breakdown, by = c("population","activity","modality")) %>%
    dplyr::mutate(
      label_md = paste0(
        "Total N = ",
        "<span style='color:", col_a, ";'>", n_a, "</span>",
        " + ",
        "<span style='color:", col_b, ";'>", n_b, "</span>",
        " = ", total_n
      )
    )
  
  ragg::agg_png(out_png, width = 3200, height = 3000, res = 300, pointsize = 11)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = dl_leaf_n,
      ggplot2::aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s,
                   fill = tranche_group, alpha = age),
      colour = "white", linewidth = 0.25
    ) +
    ggfittext::geom_fit_text(
      data = dl_leaf_n,
      ggplot2::aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s, label = lab),
      reflow = TRUE, grow = TRUE, min.size = 0, contrast = TRUE
    ) +
    ggplot2::geom_rect(
      data = dl_sgrp_n,
      ggplot2::aes(xmin = xmin_s, xmax = xmax_s, ymin = ymin_s, ymax = ymax_s, colour = sex),
      inherit.aes = FALSE, fill = NA, linewidth = 1.2
    ) +
    ggplot2::scale_colour_manual(values = col_sex, guide = "none") +
    ggtext::geom_richtext(
      data = sex_labels,
      ggplot2::aes(x = x_lab, y = y_lab, label = label_html, angle = ang),
      inherit.aes = FALSE,
      fill = NA, label.color = NA, size = 2
    ) +
    ggtext::geom_richtext(
      data = content_tops_colored,
      ggplot2::aes(x = x_mid, y = y_lab, label = label_md),
      inherit.aes = FALSE,
      fill = NA, label.color = NA, fontface = "bold", size = 4
    ) +
    ggh4x::facet_nested(
      rows = ggplot2::vars(population, activity),
      cols = ggplot2::vars(modality),
      switch = "y"
    ) +
    ggplot2::scale_fill_manual(values = c(stats::setNames(col_a, label_a),
                                          stats::setNames(col_b, label_b)),
                               drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c("10-13"=0.6,"14-17"=0.7,"18-39"=0.8,"40-59"=0.9,"60+"=1),
                                drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(-xlim_pad, 1 + xlim_pad), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.2), expand = c(0,0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE, title.position = "top"),
      alpha = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE, title.position = "top")
    ) +
    ggplot2::theme_void(base_size = 16) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold", size = 18, margin = ggplot2::margin(t = 6, b = 6)),
      strip.text.y = ggplot2::element_text(face = "bold", size = 18, margin = ggplot2::margin(l = 6, r = 6)),
      plot.title = ggplot2::element_text(size = 22, face = "bold", margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(size = 16, margin = ggplot2::margin(b = 10)),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 13, face = "bold"),
      legend.text  = ggplot2::element_text(size = 12),
      legend.key.size = grid::unit(16, "pt"),
      legend.spacing.x = grid::unit(60, "pt"),
      legend.box.spacing = grid::unit(8, "pt"),
      plot.margin = ggplot2::margin(12, 18, 12, 18),
      panel.spacing = grid::unit(8, "pt")
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = "Area = participants; color = tranche; shade = age",
      fill = "Tranche group", alpha = "Age group"
    )
  
  print(p)
  dev.off()
  message(sprintf("wrote: %s", out_png))
}

# 1) Tranche 0 vs Tranches 1–4
plot_scaled_xSex_compare(
  dl = dl,
  tranche_a = 0, tranche_b = 1:4,
  label_a = "Tranche 0", label_b = "Tranches 1–4",
  out_png = "~/Pictures/MoTrPAC_HA/tranche_sample-size_comparison_scaled_xSex_0_vs_1-4.png",
  title_txt = "Stratifying N across Tranche 0 vs Tranches 1–4"
)

# 2) Tranche 0 vs Tranches 1–5
plot_scaled_xSex_compare(
  dl = dl,
  tranche_a = 0, tranche_b = 1:5,
  label_a = "Tranche 0", label_b = "Tranches 1–5",
  out_png = "~/Pictures/MoTrPAC_HA/tranche_sample-size_comparison_scaled_xSex_0_vs_1-5.png",
  title_txt = "Stratifying N across Tranche 0 vs Tranches 1–5"
)
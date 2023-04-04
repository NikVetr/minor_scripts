#functions
logit <- function(p) log(p / (1-p))
invlogit <- function(x){ out <- exp(x) / (1 + exp(x)); ifelse(is.na(out), 1, out)}
softmax <- function(x){ ex <- exp(x); ex / sum(ex)}
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}
polar2cart <- function(t, r){
  if(length(t) == 1 & length(r) == 1){
    return(c(r*cos(t), r * sin(t)))  
  } else if(length(t) == 1){
    return(do.call(rbind, lapply(r, function(ri) polar2cart(t,ri))))
  } else if(length(r) == 1){
    return(do.call(rbind, lapply(t, function(ti) polar2cart(ti,r))))
  } else {
    if(length(t) > length(r)){
      rep(r, times = ceiling(length(t) / length(r)))[1:length(t)]
    } else {
      rep(t, times = ceiling(length(r) / length(t)))[1:length(r)]
    }
    return(do.call(rbind, lapply(1:length(t), function(i) polar2cart(t[i],r[i]))))
  }
}
xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

penalize_matches <- function(x, match_penalty_matrix){
  sum(match_penalty_matrix[x,x][upper.tri(match_penalty_matrix[x,x])])
}

#can also write function to penalize geography, but this needs all the tables and needs to be written efficiently
#can also incorporate "waypoints" through table indices, eg to get to bathrooms etc. and not just raw dists
penalize_geography <- function(prop_tables_inds, prop_indiv_inds, table_distance_cost, focal_tables, table_closenesses, tables){
    focal_indivs <- t(apply(prop_indiv_inds, 1, function(x) c(focal_tables[[1]][x[1]], focal_tables[[2]][x[2]])))
    compl_tables <- tables[-prop_tables_inds] #cos of symmetry of distance function
    distance_scores <- t(apply(prop_indiv_inds, 1, function(x) sapply(seq_along(x), function(indiv_i){
      penalize_geography_indiv(indiv = x[indiv_i],
                               indiv_table = names(tables[prop_tables_inds[indiv_i]]),
                               compl_tables = compl_tables,
                               table_closenesses = table_closenesses,
                               table_distance_cost = table_distance_cost)
      
    })))
  return(apply(distance_scores, 2, sum))
}

penalize_geography_indiv <- function(indiv, indiv_table, compl_tables, table_closenesses, table_distance_cost){
  sum(sapply(seq_along(compl_tables), function(table_j){
    sum(sapply(compl_tables[[table_j]], function(tj_indiv){
      pair_cost <- table_distance_cost[indiv, tj_indiv]
      dist_tables <- table_closenesses[indiv_table, names(compl_tables[table_j])]
      -pair_cost * dist_tables
    }))
  }))
}

#different objective functions

#want to maximize uniformity
eval_tables <- function(x, vals){
  tabx <- table(vals[x])
  sum(1-1/(2^tabx)) #geometric series 1/2 + 1/4 + 1/8 etc.
}

#want to maximize spread
eval_tables <- function(x, vals){
  var(vals[x]) 
}

#want to maximize spread *&* uniformity
eval_tables <- function(x, vals){
  -sum(dist(vals[x])^0.5) 
}

#### parameters ####

#seating parameters
n <- 100
ncat <- 4
nseats <- 6
consider_geography <- T
ntables <- max(1, floor(n / nseats))
cat_probs <- gtools::rdirichlet(1, ncat:1^1.5 * 5)
indivs <- data.frame(t(rmultinom(n, 1, cat_probs)))
colnames(indivs) <- paste0("c.", 1:ncat)
indivs$category <- apply(indivs, 1, which.max)

#table locations
sqrt_tabs <- ceiling(sqrt(ntables))
rad <- 0.5

#tables in a grid
xlocs <- floor((1:ntables + 1E-6) %% sqrt_tabs)
ylocs <- sqrt_tabs - ceiling((1:ntables + 1E-6) / sqrt_tabs)

#tables in a spiral
n_circles <- 10
n_gran <- 2000
xylocs <- polar2cart(t = seq(0, n_circles * 2 * pi, length.out = n_gran), 
                     seq(rad, rad * n_circles * 2.25, length.out = n_gran))
ti <- rep(NA, ntables)
ti[1] <- 1
for(i in 2:ntables){
  curr_loc <- xylocs[ti[i-1],]
  new_ti_disp <- min(which(sqrt(apply(t(t(xylocs[(ti[i-1]+1):n_gran,]) - curr_loc)^2, 1, sum)) > (rad*2.25)))
  new_ti <- ti[i-1] + new_ti_disp
  ti[i] <- new_ti
}
xylocs <- xylocs[ti,]
xlocs <- xylocs[,1] + sqrt_tabs / 2
ylocs <- xylocs[,2] + sqrt_tabs / 2

#get distance matrix between tables
table_distances <- as.matrix(dist(xylocs))
# table_distances <- log10(table_distances)
table_distances <- table_distances / max(table_distances)
table_closenesses <- 1 - table_distances
table_closenesses <- table_closenesses^2 #sound travels as a shell propto square of distance

#sampling parameters
logit_scale <- 1
niter <- 0.5E2 * n
n_total_rec <- 1000

#### start & run chain ####
init_match <- ceiling(1:n/nseats)
if(n%%nseats != 0){
  init_match[(ntables*nseats+1):n] <- 1:(n%%nseats)
}
tables <- init_tables <- split(seq_along(init_match), init_match)

#suppose we hate this initial location
match_penalty_matrix <- Reduce("+", lapply(tables, function(x){
  xmat <- diag(n) * -10
  xmat[x,x] <- -10
  xmat
}))

#or maybe we want to cluster similar tables together geographically
table_distance_cost <- (length(unique(indivs$category)) - as.matrix(dist(indivs$category)))
# table_distance_cost <- diag(n)
# table_distance_cost[1,10] <- table_distance_cost[10,1] <- 100
#PROBLEM -- currently, geography draws individuals together instead of pushing them apart

#initialize container variables
recorded_scores <- rep(NA, n_total_rec)
curr_table_score <- sapply(tables, eval_tables, vals = indivs$category) +
  sapply(tables, penalize_matches, match_penalty_matrix = match_penalty_matrix)
if(consider_geography){
  curr_table_score <- curr_table_score + sapply(1:ntables, function(table_i){
    sum(sapply(seq_along(tables[[table_i]]), function(indiv_i){
      sapply((1:ntables)[-table_i], function(table_j){
        penalize_geography_indiv(indiv = best_tables[[table_i]][indiv_i], 
                                 compl_tables = best_tables[table_j], 
                                 indiv_table = names(best_tables[table_i]), 
                                 table_closenesses = table_closenesses,
                                 table_distance_cost = table_distance_cost)
      })
    }))
  })
}
curr_total_score <- sum(curr_table_score)
best_total_score <- curr_total_score
best_tables <- tables

#run chain
for(i in 1:niter){
  if(i %% 100 == 0){print(i)}
  
  #pick tables to swap in proportion to their current score
  # prop_tables_inds <- sample(1:ntables, 2, replace = F, prob = rep(1/ntables, ntables))
  curr_table_score_std <- max(curr_table_score) - curr_table_score
  curr_table_score_std <- curr_table_score_std / max(curr_table_score_std)
  prop_tables_inds <- sample(1:ntables, 2, replace = F, prob = softmax(curr_table_score_std))
  
  max_possible_swap <- min(sapply(tables[prop_tables_inds], length)) - 1
  n_indiv_to_swap <- sample(1:max_possible_swap, size = 1, prob = dexp(1:max_possible_swap, rate = 0.4))
  
  prop_indiv_inds <- cbind(sample(1:length(tables[[prop_tables_inds[1]]]), n_indiv_to_swap), 
                       sample(1:length(tables[[prop_tables_inds[2]]]), n_indiv_to_swap))
  
  curr_tables <- prop_tables <- list(
                      tables[[prop_tables_inds[1]]],
                      tables[[prop_tables_inds[2]]]
                      )
  
  for(swap_i in 1:n_indiv_to_swap){
    prop_tables[[1]][prop_indiv_inds[swap_i, 1]] <- curr_tables[[2]][prop_indiv_inds[swap_i, 2]]
    prop_tables[[2]][prop_indiv_inds[swap_i, 2]] <- curr_tables[[1]][prop_indiv_inds[swap_i, 1]]  
  }
  
  # if(consider_geography){
  #   curr_subtable_score <- sapply(curr_tables, eval_tables, vals = indivs$category) +
  #     sapply(curr_tables, penalize_matches, match_penalty_matrix = match_penalty_matrix) +
  #     penalize_geography(prop_tables_inds = prop_tables_inds, prop_indiv_inds = prop_indiv_inds,
  #                        table_distance_cost = table_distance_cost, focal_tables = curr_tables,
  #                        table_closenesses = table_closenesses, tables = tables)    
  # } else {
  #  curr_subtable_score <- curr_table_score[prop_tables_inds]  
  # }
  
  curr_subtable_score <- curr_table_score[prop_tables_inds] 

  prop_subtable_score <- sapply(prop_tables, eval_tables, vals = indivs$category) +
    sapply(prop_tables, penalize_matches, match_penalty_matrix = match_penalty_matrix) +
    ifelse2(consider_geography,
      penalize_geography(prop_tables_inds = prop_tables_inds, prop_indiv_inds = prop_indiv_inds,
                         table_distance_cost = table_distance_cost, focal_tables = prop_tables,
                         table_closenesses = table_closenesses, tables = tables),
    0)
  
  diff_subtable_score <- prop_subtable_score - curr_subtable_score
  diff_score <- sum(diff_subtable_score)
  
  prob_swap <- invlogit(diff_score * logit_scale)
  swap <- sample(c(T,F), 1, prob = c(prob_swap, 1-prob_swap))
  if(swap){
    tables[[prop_tables_inds[1]]] <- prop_tables[[1]]
    tables[[prop_tables_inds[2]]] <- prop_tables[[2]]
    curr_total_score <- curr_total_score + diff_score
    curr_table_score[prop_tables_inds] <- curr_table_score[prop_tables_inds] + diff_subtable_score
    if(curr_total_score > best_total_score){
      best_total_score <- curr_total_score
      best_tables <- tables
    }
  }
  
  if(i %% floor(niter / n_total_rec) == 1){
    rec_ind <- floor(i / floor(niter / n_total_rec)) + 1
    recorded_scores[rec_ind] <- curr_total_score
  }
  
}

#quick check of success
par(mar = c(5,4,2,2))
plot(tail(recorded_scores, ceiling(n_total_rec * 0.8)), type = "l")
mean(sapply(init_tables, function(x) length(table(indivs$category[x]))))
mean(sapply(best_tables, function(x) length(table(indivs$category[x]))))

#### now optimize seating *within* each table ####

#### visualize ####

#plot params
person_cols <- viridisLite::inferno(ncat, alpha = 0.5)
table_cols <- colorRampPalette(c("red", "green"))(100)
plot_table_edges <- T
n_pts_per_wedge <- 5
cex_person <- 125 / n
cex_table <- 25 / ntables 

#table goodness
intratable_compatibilities <- sapply(best_tables, eval_tables, vals = indivs$category) +
  sapply(tables, penalize_matches, match_penalty_matrix = match_penalty_matrix)
intratable_compatibilities_std <- intratable_compatibilities - min(intratable_compatibilities)
intratable_compatibilities_std <- ceiling(((intratable_compatibilities_std / max(intratable_compatibilities_std)) * 0.999 + 0.001) * 100)

#table loc goodness
intertable_compatibilities <- do.call(rbind, lapply(1:(ntables-1), function(table_i){
  c(rep(0, table_i), sapply((table_i+1):ntables, function(table_j){
    mean(sapply(seq_along(best_tables[[table_i]]), function(indiv_i){
      penalize_geography_indiv(indiv = best_tables[[table_i]][indiv_i], 
                               compl_tables = best_tables[table_j], 
                               indiv_table = names(best_tables[table_i]), 
                               table_closenesses = table_closenesses,
                               table_distance_cost = table_distance_cost)
    }))
  }))
}))
intertable_compatibilities <- rbind(intertable_compatibilities, 0) + t(rbind(intertable_compatibilities, 0))
intertable_compatibilities_std <- intertable_compatibilities - min(intertable_compatibilities)
intertable_compatibilities_std <- ceiling(((intertable_compatibilities_std / max(intertable_compatibilities_std)) * 0.999 + 0.001) * 100)
intertable_compatibilities_std <- 101 - intertable_compatibilities_std

edge_thresh <- quantile(intertable_compatibilities[upper.tri(intertable_compatibilities)], 
                        probs = c(0.001, 0.999))
edges_to_plot <- which(((intertable_compatibilities < edge_thresh[1] | 
                         intertable_compatibilities > edge_thresh[2])) &
                         table_distances < quantile(table_distances[upper.tri(table_distances)], probs = 0.1), 
                       arr.ind = T)
edges_to_plot <- edges_to_plot[edges_to_plot[,1] > edges_to_plot[,2],]
if(is.null(nrow(edges_to_plot))){edges_to_plot <- cbind(edges_to_plot[1], edges_to_plot[2])}



# a = sapply(seq_along(best_tables), function(table_i){
#   sapply(seq_along(best_tables), function(table_j){
#     mean(sapply(seq_along(best_tables[[table_i]]), function(indiv_i){
#       penalize_geography_indiv(indiv = best_tables[[table_i]][indiv_i], 
#                                compl_tables = best_tables[table_j], 
#                                indiv_table = names(best_tables[table_i]), 
#                                table_closenesses = table_closenesses,
#                                table_distance_cost = indiv_penalty_matrix)
#     }))
#   })
# })


#now do the actual plotting
par(mar = rep(1,4))
plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="",
     xlim = range(xlocs) + c(-rad/2, rad/2), ylim = range(ylocs) + c(-rad/2, rad/2))
# lines(xylocs + sqrt_tabs/2, xpd = NA)

#initial edges
if(length(edges_to_plot) != 0){
  for(k in 1:nrow(edges_to_plot)){
    i <- edges_to_plot[k,1]
    j <- edges_to_plot[k,2]
    segments(x0 = xlocs[i], y0 = ylocs[i], x1 = xlocs[j], y1 = ylocs[j],
             col = table_cols[intertable_compatibilities_std[i, j]], lwd = 2)
  }
}


#now the nodes
for(i in 1:ntables){
  xloc <- xlocs[i]
  yloc <- ylocs[i]
  x <- best_tables[[i]]
  wedge_angles <- seq(0, 2*pi, length.out = length(x) + 1)
  table_degrees <- polar2cart(seq(0, 2*pi, length.out = n_pts_per_wedge * (length(x) + 1)), rad)
  # polygon(x = c(xloc, xloc + wedge_degrees[,1] * xyrat(), xloc),
  #         y = c(yloc, yloc + wedge_degrees[,2], yloc), xpd = NA,
  #         col = person_cols[indivs$category[x[j]]]
  # )
  
  for(j in 1:length(x)){
    wedge_degrees <- polar2cart(seq(wedge_angles[j], wedge_angles[j+1], length.out = n_pts_per_wedge), rad)
    polygon(x = c(xloc, xloc + wedge_degrees[,1] * xyrat(), xloc),
            y = c(yloc, yloc + wedge_degrees[,2], yloc), xpd = NA,
            col = "white"
    ) #<--- inefficient, move to commented code above
    polygon(x = c(xloc, xloc + wedge_degrees[,1] * xyrat(), xloc),
            y = c(yloc, yloc + wedge_degrees[,2], yloc), xpd = NA,
            col = person_cols[indivs$category[x[j]]]
    )
    text(labels = x[j], 
         x = xloc + polar2cart(mean(c(wedge_angles[j], wedge_angles[j+1])), rad / 1.5)[1] * xyrat(),
         y = yloc + polar2cart(mean(c(wedge_angles[j], wedge_angles[j+1])), rad / 1.5)[2],
         cex = cex_person, xpd = NA)
  }
}

#plot the edges
if(length(edges_to_plot) != 0){
  for(k in 1:nrow(edges_to_plot)){
    i <- edges_to_plot[k,1]
    j <- edges_to_plot[k,2]
    segments(x0 = xlocs[i], y0 = ylocs[i], x1 = xlocs[j], y1 = ylocs[j],
             col = adjustcolor(table_cols[intertable_compatibilities_std[i, j]], 0.3),
             lty = 2, lwd = 2)
  }
}

#plot the table labels
for(i in 1:ntables){
  xloc <- xlocs[i]
  yloc <- ylocs[i]
  x <- best_tables[[i]]
  lab_wedge_degrees <- polar2cart(seq(0, 2*pi, length.out = 10), rad / 3)
  polygon(x = xloc + lab_wedge_degrees[,1] * xyrat(),
          y = yloc + lab_wedge_degrees[,2], xpd = NA,
          col = "white")
  polygon(x = xloc + lab_wedge_degrees[,1] * xyrat(),
          y = yloc + lab_wedge_degrees[,2], xpd = NA,
          col = adjustcolor(table_cols[intratable_compatibilities_std[i]], 0.2))
  text(xloc, yloc, i, cex = cex_table, font = 2, xpd = NA, col = table_cols[intratable_compatibilities_std[i]])
}
  
legend(x = "topright",
       pt.bg = person_cols, legend = 1:ncat, pch = 24, pt.cex = 1.5, box.lwd = 0, title = "Category")



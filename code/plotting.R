choose_phenotype_colors <- function(membership, nspp){
  # Choose colors
  # max 6 species, 9 phenotypes each
  mycolors <- character(0)
  sppcolors <- character(0)
  sppdarkcolors <- character(0)
  palettenames <- c("Blues", "Reds", "Purples", "Oranges", "Greens",  "Reds", "Greys")
  for (i in 1:nspp){
    npheno <- sum(membership == i)
    tmpcol <- brewer.pal(9, palettenames[i])[8-(0:(npheno - 1))]
    # make color transparent
    sppdarkcolors <- c(sppdarkcolors, tmpcol[floor(npheno / 2 + 1)])
    sc <- col2rgb(tmpcol[floor(npheno / 2 + 1)])/255
    sppcolors <- c(sppcolors, rgb(sc[1],sc[2],sc[3], 0.5))
    
    mycolors <- c(mycolors, tmpcol)
  }
  return(list(mycolors = mycolors, sppcolors = sppcolors, sppdarkcolors = sppdarkcolors))
}

plot_matrix_H <- function(m, H){
  n <- length(m) # number of species
  nm <- sum(m) # total number of phenotypes
  # create an expanded vector of memberships
  membership <- numeric(0)
  labels <- character(0)
  for (i in 1:n) {
    membership <- c(membership, rep(i, m[i]))
    labels <- c(labels, paste0(letters[i], " ", 1:m[i]))
  }
  tmp <- choose_phenotype_colors(membership, n)
  mycolors <- tmp$mycolors
  sppcolors <- tmp$sppcolors
  colnames(H) <- labels
  rownames(H) <- labels
  mH <- reshape::melt(H)
  mH$X1 <- factor(mH$X1, labels)
  mH$X2 <- factor(mH$X2, rev(labels))
  pl_H <- ggplot(mH) + 
    aes(x = X1, y = X2, fill = value) + 
    geom_tile(colour = "grey") + 
    #scale_y_reverse() + 
    scale_fill_gradientn("probability of winning",
      colours = c(muted("red"), "white", muted("blue")), 
      values = rescale(c(-1, 0, 1))) + 
    #theme_void() + 
    theme(legend.position = "bottom") + 
    coord_equal() + 
    scale_x_discrete("", position = "top") + 
    scale_y_discrete("") + 
    theme(axis.text.x = element_text(colour = mycolors, face = "bold"),
          axis.text.y = element_text(colour = rev(mycolors), face = "bold"))
  return(pl_H)
}

plot_matrix_Q <- function(m, Q){
  n <- length(m) # number of species
  nm <- sum(m) # total number of phenotypes
  # create an expanded vector of memberships
  membership <- numeric(0)
  labels <- character(0)
  for (i in 1:n) {
    membership <- c(membership, rep(i, m[i]))
    labels <- c(labels, paste0(letters[i], " ", 1:m[i]))
  }
  tmp <- choose_phenotype_colors(membership, n)
  mycolors <- tmp$mycolors
  sppcolors <- tmp$sppcolors
  colnames(Q) <- labels
  rownames(Q) <- labels
  mQ <- reshape::melt(Q)
  mQ$X1 <- factor(mQ$X1, labels)
  mQ$X2 <- factor(mQ$X2, rev(labels))
  pl_Q <- ggplot(mQ) + 
    aes(x = X1, y = X2, fill = value) + 
    geom_tile(colour = "grey") + 
    #scale_y_reverse() + 
    scale_fill_gradient("probability", low = "white", high = "black") +
    #theme_void() + 
    theme(legend.position = "bottom") + 
    coord_equal() + 
    scale_x_discrete("", position = "top") + 
    scale_y_discrete("") + 
    theme(axis.text.x = element_text(colour = mycolors, face = "bold"),
          axis.text.y = element_text(colour = rev(mycolors), face = "bold"))
  return(pl_Q)
}

plot_competition_graph_phenotypes <- function(
  m, # vector where membership[i] is the species phenotype i belongs to 
  H # matrix where Hij = prob i beats j
  ){
  # get colors
  n <- length(m) # number of species
  nm<-sum(m) # total number of phenotypes
  # create an expanded vector of memberships
  membership <- numeric(0)
  for (i in 1:n) {
    membership <- c(membership, rep(i, m[i]))
  }
  tmp <- choose_phenotype_colors(membership, n)
  mycolors <- tmp$mycolors
  sppcolors <- tmp$sppcolors
  # Build graph
  Pay <- H - t(H)
  Pay[Pay <= 0] <- 0
  g <- graph_from_adjacency_matrix(Pay, weighted = TRUE, mode = "directed")
  V(g)$color <- mycolors
  coords <- layout_in_circle(g, order = V(g))
  return(
  plot(g, mark.groups = split(1:vcount(g), membership), layout = coords, 
       mark.col = sppcolors, mark.border = sppcolors,
       edge.width= 3 * E(g)$weight, edge.color = "black",vertex.label=NA))
}

plot_dynamics <- function(comdyn, 
                          m, 
                          normalize = TRUE, 
                          trans = "identity"){
  n <- length(m)
  nm <- sum(m)
  membership <- numeric(0)
  for (i in 1:n) {
    membership <- c(membership, rep(i, m[i]))
  } 
  tmp <- choose_phenotype_colors(membership, n)
  # rescale to within 0,1 to account for round-off error?
  if(normalize){
    rs <- apply(comdyn[,-1], 1, sum)
    comdyn[,-1] <- sweep(comdyn[,-1], 1, rs, FUN="/")
  }
  abund_mat<- data.frame(comdyn) %>% melt(id = "time") %>% dplyr::rename(phenotype = variable)
  return(ggplot(abund_mat, aes(x = as.numeric(time), y = value + 1e-16, colour = phenotype)) + 
           geom_line(alpha = 0.7, size = 1) + 
           theme(legend.position = "none") + 
           scale_y_continuous(trans = trans, name = "Time") + 
           xlab("Relative abundance") + 
           scale_colour_manual(values = tmp$mycolors) + 
           geom_hline(yintercept = 0))
}



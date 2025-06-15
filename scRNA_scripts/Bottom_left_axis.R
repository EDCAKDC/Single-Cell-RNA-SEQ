left_axes <- function(scRNA) {
  # Extract UMAP coordinates
  umap_coords <- Embeddings(object = scRNA, reduction = 'umap') %>%
    data.frame()
  
  # Preview first few coordinates
  head(umap_coords, 3)
  
  # Calculate the lower bound for axes
  lower_bound <- floor(min(min(umap_coords$UMAP_1), min(umap_coords$UMAP_2))) - 2
  
  # Calculate relative line length
  line_length <- abs(0.3 * lower_bound) + lower_bound
  
  # Calculate midpoint for axis labels
  midpoint <- abs(0.3 * lower_bound) / 2 + lower_bound
  
  # Construct axis lines data
  axes <- data.frame(
    x = c(lower_bound, lower_bound, lower_bound, line_length),
    y = c(lower_bound, line_length, lower_bound, lower_bound),
    group = c(1, 1, 2, 2),
    label = rep(c('UMAP_2', 'UMAP_1'), each = 2)
  )
  
  # Construct axis labels
  label <- data.frame(
    lab = c('UMAP_2', 'UMAP_1'),
    angle = c(90, 0),
    x = c(lower_bound - 3, midpoint),
    y = c(midpoint, lower_bound - 2.5)
  )
  
  return(list(axes = axes, label = label))
}

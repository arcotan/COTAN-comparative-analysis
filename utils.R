library(dplyr)

write_clustering = function(outdir, label_df, cell_col, cluster_col) {
  to_write = label_df[c(cell_col, cluster_col)]
  colnames(to_write)[colnames(to_write) == cell_col] = "cell"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/clustering_labels.csv", sep=""), row.names = FALSE)
}

write_markers = function(outdir, marker_df, gene_col, cluster_col, order_by_col, is_higher_better, top_k) {
  order_coef = -1
  if (is_higher_better) {
    order_coef = 1
  }
  marker_df = marker_df %>%
    group_by(get(cluster_col)) %>%
    slice_max(n = top_k, order_by = (order_coef*get(order_by_col))) %>% 
    mutate(rank = row_number(), ties.method = "first")
  marker_df = data.frame(marker_df)
  to_write = marker_df[c(gene_col, cluster_col, "rank")]
  colnames(to_write)[colnames(to_write) == gene_col] = "gene"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/markers.csv", sep=""), row.names = FALSE)
}
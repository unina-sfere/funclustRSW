test_data <- read_csv("test_data.csv")
test_mat <- as.matrix(test_data)

matplot((data$X), col = cl_5_hc, type = "l")

cl_5_hc_new <- cl_5_hc
cl_5_hc_new[cl_5_hc == 1] <- 3
cl_5_hc_new[cl_5_hc == 2] <- 1
cl_5_hc_new[cl_5_hc == 3] <- 2


delta <- diff(data$grid)[1]

cluster_assignment <- sapply(1:3, function(ii) {
  colSums(((test_mat -
              (get_centroids_df(cl_5_hc_new, data, "filtering B-spline hc") %>% 
                 filter(cluster == ii) %>% 
                 pull(Resistance))) ^ 2) * delta)
}) %>% 
  apply(1, which.min)

distances <- lapply(1:3, function(ii) {
  colSums(((data$X[, cl_5_hc_new == ii] -
              (get_centroids_df(cl_5_hc_new, data, "filtering B-spline hc") %>% 
                 filter(cluster == ii) %>% 
                 pull(Resistance))) ^ 2) * delta)
}) %>% 
  setNames(1:3)


distance_test <- sapply(1:3, function(ii) {
  colSums(((test_mat -
              (get_centroids_df(cl_5_hc_new, data, "filtering B-spline hc") %>% 
                 filter(cluster == ii) %>% 
                 pull(Resistance))) ^ 2) * delta)
})

limits <- sapply(distances, function(x) quantile(x, .99))

control_chart_df <-  lapply(1:3, function(ii) 
  data.frame(distance = distance_test[cluster_assignment == ii, ii],
             cluster = ii,
             limit = limits[ii])
) %>% 
  bind_rows %>% 
  mutate(obs = rownames(.),
         cluster = factor(cluster)) %>% 
  arrange(obs) %>% 
  mutate(obs = 1:n()) 

control_chart <- control_chart_df %>% 
  ggplot +
  geom_line(aes(obs, distance), col = "black") +
  # geom_point(aes(obs, limit, col = cluster), shape = "-", size = 10) +
  geom_line(aes(obs, limit), col = "black", lty = 2) +
  geom_point(aes(obs, distance, col = cluster)) +
  theme_bw() +
  scale_x_continuous(breaks = 1:100) +
  scale_color_brewer(palette = "Set1")
control_chart
# ggsave("control_chart.png", control_chart, width = 5, height = 3)
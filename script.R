# libraries and functions ----------------------------------------------------
if (!("curvclust" %in% installed.packages())) {
  install.packages("./curvclust_0.0.1.tar.gz", repos = NULL, type="source")
}
library(curvclust)
library(funHDDC)
library(clValid)
library(mclust)
library(NbClust)
library(MASS)
library(fda)
library(parallel)
library(pbmcapply)
library(matrixcalc)
library(pbapply)
library(fda.usc)
library(readxl)
library(tidyverse)
library(tidyfun)
library(mgcv)
library(pbmcapply)
library(factoextra)
library(gridExtra)
library(cowplot)

source("functions.R")
source("fclust.R")

# load data ---------------------------------------------------------------
df <- read.csv("data.csv")
names(df) <- gsub("X", "", names(df))
data <- list(X = as.matrix(df[, - 1]), grid = df$x)

basis<- create.bspline.basis(c(0,1), nbasis=200)
loglam         = seq(-10, -4, 0.25)
Gcvsave        = numeric()
for(i in 1:length(loglam)){
  fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
  Sm.i    = smooth.basis(data$grid, data$X, fdPari)
  Gcvsave[i] = sum(Sm.i$gcv)
  
}
lambda_s=10^loglam[which.min(Gcvsave)]
fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
data$X_fd<-smooth.basis(data$grid, data$X, fdPari)$fd
# plot(data$X_fd)


# plots -------------------------------------------------------------------
plot_analysis <- data$X %>%
  as.data.frame %>%
  mutate(x = 1:nrow(data$X)) %>%
  pivot_longer(-x) %>%
  ggplot +
  geom_line(
    aes(x, value, group = name),
    alpha = .75,
    size = .3,
    col = "darkgrey"
  ) +
  theme_bw() +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  xlab("Time [ms]") +
  ylim(c(260, 320)) +
  theme(aspect.ratio=1)
# plot_analysis
# ggsave("plots/plot_analysis.png", plot_analysis, width = 5, height = 4)

df_ref <- read.csv("reference_profile.csv") %>% 
  mutate(time = 1:n())
which_min <- df_ref %>% 
  filter(time < 100) %>% 
  pull(Resistance) %>% 
  which.min()
this_min <- df_ref$Resistance[which_min]
which_max <- df_ref %>% 
  mutate(time = 1:n()) %>%
  filter(time > 50) %>% 
  filter(Resistance == max(Resistance)) %>% 
  pull(time)
this_max <- df_ref$Resistance[which_max]
which_end <- 238
this_end <- last(df_ref$Resistance)
which_mp <- which.max(diff(df_ref$Resistance))
this_mp <- df_ref$Resistance[which_mp]
plot_ref_data <- df_ref %>% 
  ggplot +
  theme_bw() +
  ylab(expression("Resistance " ~ paste("[", mu, Omega, "]"))) +
  xlab("Time [ms]") +
  # ylim(c(260, 320)) +
  theme(aspect.ratio=1) +
  geom_segment(aes(x = which_max, xend = which_max, y = this_min, yend = this_max), 
               lty = 2, 
               size = .25,
               colour = "darkgrey",
               arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  geom_text(aes(x = which_max + 9, y = (this_min + this_max) / 2, label = "(I)"),
            colour = "darkgrey") +
  geom_segment(aes(x = which_min, xend = which_max, y = this_min, yend = this_min), 
               lty = 2, 
               size = .25,
               colour = "darkgrey",
               arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  geom_segment(aes(x = which_max, xend = 238, y = this_min, yend = this_min), 
               lty = 2, 
               size = .25,
               colour = "darkgrey",
               arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  geom_segment(aes(x = 238, xend = 238, y = this_end - 3, yend = this_min + 3), 
               lty = 2, 
               size = .25,
               colour = "darkgrey") +
  geom_text(aes(x = (which_min + which_max) / 2, y = this_min - 1, label = "(II)"), 
            colour = "darkgrey") +
  # geom_text(aes(x = which_end + 12, y = this_end, label = "(III)")) +
  geom_text(aes(x = (238 - which_min) / 2 + which_min, y = this_min - 1, label = "(IV)"), 
            colour = "darkgrey") +
  geom_point(aes(which_mp, this_mp), colour = "darkgrey") +
  geom_point(aes(which_max, this_max), colour = "darkgrey") +
  geom_line(aes(time, Resistance)) +
  geom_text(aes(x = which_mp - 12, y = this_mp, label = "MP"), colour = "darkgrey") +
  geom_text(aes(x = which_max, y = this_max + 1.35, label = "NF"), colour = "darkgrey") +
  scale_y_continuous(breaks = c(270, this_end, 280, 290),
                     labels = c("270", "(III)", "280", "290")) +
  geom_hline(aes(yintercept = this_end), lty = 2, colour = "darkgrey") +
  theme_bw() +
  theme(panel.grid = element_blank())
plot_ref_data
# ggsave("plots/plot_ref_data_new.pdf", plot_ref_data, width = 5, height = 4)



num_cluster_seq <- 2:10

###############################
# Clustering methods ------------------------------------------------------
###############################
#############
# Adaptive ------------------------------------------
#############
# fclust James Sugar (2003) JASA ------------------------------------------
nrows <- nrow(data$X)
seq_vec <- seq(1,nrows,by=13)
dat <- data
dat$X <- data$X[seq_vec, ]
dat$grid <- data$grid[seq_vec]


ncores <- 1
tic <- proc.time()
mod1<-fit_fclust_ms(data=dat, num_cluster_seq = num_cluster_seq, dim_seq = c(5, 10, 15),
                    ncores = ncores)
mod1$toc <- proc.time() - tic

# Giacofci (2012) Biometrics curvclust -----------------------------------
parameters <- expand.grid(
  structures = c("constant", "group", "scale.location", "group.scale.location", "none"), 
  mixed = c(TRUE, FALSE), reduction = c(TRUE, FALSE)) %>% 
  filter(!(structures == "none" & mixed),
         !(structures != "none" & !mixed)) %>% 
  mutate(structures = as.character(structures))

mod2 <- vector(mode = "list", length = nrow(parameters))
names(mod2) <- 
  sapply(1:nrow(parameters), function(ii) parameters[ii, ] %>%
           mutate(
             structures = paste0("structure_", structures),
             mixed = paste0("mixed_", mixed),
             reduction = paste0("reduction_", reduction)
           ) %>% 
           paste0(collapse = " "))

for (ii in 1:nrow(parameters)) {
  tic <- proc.time()
  mod2[[ii]] <- curvclust_ms(data=data,
                             num_cluster_seq = c(1, num_cluster_seq),
                             structure = parameters$structures[ii],
                             mixed = parameters$mixed[ii], 
                             reduction = parameters$reduction[ii])
  mod2[[ii]]$toc <- proc.time() - tic
}

bic <- mod2 %>% sapply(function(x) x$BIC)

## n. clusters
bic %>% 
  as.data.frame %>% 
  mutate(K = c(1, num_cluster_seq)) %>% 
  pivot_longer(- K, values_to = "BIC", names_to = "model") %>% 
  ggplot +
  geom_line(aes(K, BIC, col = model))

arrind <- which(bic == max(bic), arr.ind = TRUE)
modopt <- mod2[[arrind[2]]]
Kopt <- which.max(modopt$BIC)[1] + 1

# Bouveyron (2011) ADAC -------------------------------
## FPCA scores funHDDC
tic <- proc.time()
mod3<-fit_funHDDC_ms(data=data,num_cluster_seq = num_cluster_seq,
                     model = c('AkjBkQkDk','AkjBQkDk', 'AkBkQkDk', 
                               'ABkQkDk', 'AkBQkDk', 'ABQkDk'),
                     threshold_seq = c(.2, .5, .9),
                     nb.rep = 20)
mod3$toc <- proc.time() - tic


#############
# Distance based ------------------------------------------
#############
# Distance based (fda.usc) ------------------------------------------------
dat <- data
colnames(dat$X) <- 1:ncol(dat$X)
set.seed(0)
tic <- proc.time()
mod4<-distance_ms(data = dat, num_cluster_seq = num_cluster_seq,
                  met="other", ncores = 1)
mod4$toc <- proc.time() - tic


#############
# Filtering ------------------------------------------
#############
# Filtering B-spline ------------------------------------------------------
tic <- proc.time()
## 12 basis functions elbow Gcv
mod5 <- fil_bspline_ms_nclust(data = data, num_cluster_seq = 2:10, nbasis = 12)
mod5$toc <- proc.time() - tic



# Filtering PCA -----------------------------------------------------------
tic <- proc.time()
mod6 <- fil_fpca_ss_nbclust(data = dat, num_cluster_seq = num_cluster_seq, per_comp = 0.99)
mod6$toc <- proc.time() - tic


#############
# Raw data ------------------------------------------
#############
# Raw data ------------------------------------------------
tic <- proc.time()
dat <- data
dat$X <- dat$X[seq(1, 238, by = 13), ]
mod7<-raw_ms_nbclust(data = dat, num_cluster_seq = num_cluster_seq)
mod7$toc <- proc.time() - tic



#############################
# Results -----------------------------------------------------------------
#############################
## mod1 fclust
cl_1 <- mod1$class_opt_BIC

## mod2 curvclust
bic <- mod2 %>% sapply(function(x) x$BIC)
arrind <- which(bic == max(bic), arr.ind = TRUE)
modopt <- mod2[[arrind[2]]]
Kopt <- which.max(modopt$BIC)[1] + 1
cl_2 <- modopt$class[[which.max(modopt$BIC)[1]]]

## mod3 adaptive funHDDC
cl_3 <- mod3$mod$class

## mod4 distance-based
cl_4 <- mod4$mod_opt_sil$clus

## mod5 filtering B-spline
cl_5_km <- mod5$mod_opt$ind_km$cluster
cl_5_hc <- mod5$mod_opt$ind_hc$cluster
cl_5_mb <- mod5$mod_opt$mod_opt$classification

## mod6 filtering PCA
cl_6_km <- mod6$mod_opt$ind_km$cluster
cl_6_hc <- mod6$mod_opt$ind_hc$cluster
cl_6_mb <- mod6$mod_opt$mod_opt$classification

## mod7 raw
cl_7_km <- mod7$mod_opt$ind_km$cluster
cl_7_hc <- mod7$mod_opt$ind_hc$cluster
cl_7_mb <- mod7$mod_opt$mod_opt$classification

get_centroids_df <- function(cl, data, method) {
  df <- data$X %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(class = cl) %>% 
    group_by(class) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-class) %>% 
    as.matrix() %>% 
    t() %>%
    as.data.frame %>% 
    setNames(names(table(cl))) %>% 
    mutate(time = 1:n()) %>% 
    pivot_longer(- time, names_to = "cluster", values_to = "Resistance") %>% 
    mutate(mod = method) %>% 
    mutate(cluster = factor(cluster))
  levels(df$cluster) <- data.frame(ranked_clusters = df %>% 
                                     filter(time == 21) %>% 
                                     arrange(desc(Resistance)) %>% 
                                     pull(cluster)) %>% 
    mutate(clusters = 1:n()) %>% 
    arrange(ranked_clusters) %>% 
    pull(clusters)
  df
  
}


get_clustered_funs_df <- function(cl, data, method) {
  
  df_centroids <- data$X %>% 
    t %>% as.data.frame %>% 
    mutate(class = cl) %>% 
    group_by(class) %>% 
    summarise_all(mean) %>% 
    dplyr::select(- class) %>% as.matrix %>% t %>%
    as.data.frame %>% 
    setNames(names(table(cl))) %>% 
    mutate(time = 1:n()) %>% 
    pivot_longer(- time, names_to = "cluster", values_to = "Resistance") %>% 
    mutate(mod = method) %>% 
    mutate(cluster = factor(cluster))
  
  df <- data$X %>% 
    as.data.frame %>% 
    mutate(time = 1:n()) %>% 
    pivot_longer(- time, names_to = "obs", values_to = "Resistance") %>% 
    inner_join(data.frame(obs = colnames(data$X),
                          cluster = cl), by = "obs") %>% 
    mutate(cluster = as.character(cluster),
           mod = method) %>% 
    mutate(cluster = factor(cluster))
  
  levels(df$cluster) <- data.frame(ranked_clusters = df_centroids %>% 
                                     filter(time == 21) %>% 
                                     arrange(desc(Resistance)) %>% 
                                     pull(cluster)) %>% 
    mutate(clusters = 1:n()) %>% 
    arrange(ranked_clusters) %>% 
    pull(clusters)
  
  df %>% 
    mutate(cluster = cluster %>% as.character %>% factor)
}

lay <- rbind(c(rep(1, 496), rep(NA, 104)),
             c(rep(2, 600)),
             c(rep(3, 496), rep(NA, 104)))
p_cen <- grid.arrange(
  bind_rows(
    get_centroids_df(cl_3, data, "adaptive\nfunHDDC"),
    get_centroids_df(cl_4, data, "distance-based"),
    get_centroids_df(cl_7_hc, data, "raw\nhierarchical"),
    get_centroids_df(cl_2, data, "adaptive\ncurvclust")) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot() +
    geom_line(aes(time, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  bind_rows(
    get_centroids_df(cl_5_hc, data, "filtering B-spline\nhierarchical"),
    get_centroids_df(cl_5_km, data, "filtering B-spline\nk-means"),
    get_centroids_df(cl_6_hc, data, "filtering FPCA\nhierarchical"),
    get_centroids_df(cl_6_km, data, "filtering FPCA\nk-means"),
    get_centroids_df(cl_7_km, data, "raw\nk-means"),
  ) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot +
    geom_line(aes(time, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  bind_rows(
    get_centroids_df(cl_1, data, "adaptive\nfclust"),
    get_centroids_df(cl_5_mb, data, "filtering B-spline\nmodel-based"),
    get_centroids_df(cl_6_mb, data, "filtering FPCA\nmodel-based"),
    get_centroids_df(cl_7_mb, data, "raw\nmodel-based")
  ) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot +
    geom_line(aes(time, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  layout_matrix = lay
)
# ggsave("plots/centroids.pdf", p_cen, width = 10, height = 8)

lay <- rbind(c(rep(1, 496), rep(NA, 104)),
             c(rep(2, 600)),
             c(rep(3, 496), rep(NA, 104)))
p_funs <- grid.arrange(
  
  bind_rows(
    get_clustered_funs_df(cl_3, data, "adaptive\nfunHDDC"),
    get_clustered_funs_df(cl_4, data, "distance-based"),
    get_clustered_funs_df(cl_7_hc, data, "raw\nhierarchical"),
    get_clustered_funs_df(cl_2, data, "adaptive\ncurvclust"),
  ) %>% 
    ggplot +
    geom_line(aes(time, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  
  bind_rows(
    get_clustered_funs_df(cl_5_hc, data, "filtering B-spline\nhierarchical"),
    get_clustered_funs_df(cl_5_km, data, "filtering B-spline\nk-means"),
    get_clustered_funs_df(cl_6_hc, data, "filtering FPCA\nhierarchical"),
    get_clustered_funs_df(cl_6_km, data, "filtering FPCA\nk-means"),
    get_clustered_funs_df(cl_7_km, data, "raw\nk-means"),
  ) %>% 
    ggplot +
    geom_line(aes(time, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  
  bind_rows(
    get_clustered_funs_df(cl_1, data, "adaptive\nfclust"),
    get_clustered_funs_df(cl_5_mb, data, "filtering B-spline\nmodel-based"),
    get_clustered_funs_df(cl_6_mb, data, "filtering FPCA\nmodel-based"),
    get_clustered_funs_df(cl_7_mb, data, "raw\nmodel-based")
  ) %>% 
    ggplot +
    geom_line(aes(time, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("Time [ms]") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
    theme(aspect.ratio=1),
  layout_matrix = lay
)
# ggsave("plots/col_functions.png", p_funs, width = 10, height = 8)


# pad ---------------------------------------------------------------------
pad_profiles <- read.csv("pad_profiles.csv")

df_centr <- data$X %>% 
  t %>% as.data.frame %>% 
  mutate(class = cl_5_hc) %>% 
  group_by(class) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-class) %>% as.matrix %>% t %>%
  as.data.frame %>% 
  setNames(names(table(cl_5_hc))) %>% 
  mutate(x = 1:n()) %>% 
  pivot_longer(-x, names_to = "cluster", values_to = "Resistance") %>% 
  mutate(mod = "filtering B-spline hc") %>% 
  mutate(cluster = factor(cluster)) %>% 
  mutate(obs = paste("centroid", cluster))
ranked_clusters <- data.frame(ranked_clusters = df_centr %>% 
                                filter(x == 21) %>% 
                                arrange(desc(Resistance)) %>% 
                                pull(cluster)) %>% 
  mutate(clusters = 1:n()) %>% 
  arrange(ranked_clusters)
levels(df_centr$cluster) <- ranked_clusters$clusters

df_obs <- data$X[, pad_profiles$pad_profiles] %>% 
  as.data.frame %>% 
  mutate(x = 1:n()) %>% 
  pivot_longer(- x, names_to = "obs", values_to = "Resistance") %>% 
  inner_join(data.frame(cluster = cl_5_hc, obs = as.character(seq_along(cl_5_hc))),
             by = "obs") %>%
  inner_join(data.frame(`wear level` = pad_profiles$wear.level,
             obs = as.character(pad_profiles$pad_profiles)),
             by = "obs") %>% 
  mutate(cluster = cluster %>% as.character %>% factor,
         type = "observation")
levels(df_obs$cluster) <- ranked_clusters$clusters

plot_pad0 <- df_obs %>% mutate(type = "wear level") %>% 
  mutate(cluster = ifelse(cluster == 1, "just renewed",
                          ifelse(cluster == 2, "medium wear", "severe wear"))) %>% 
  mutate(cluster = factor(cluster, levels = c("just renewed", "medium wear", "severe wear"))) %>% 
  rename("wear level" = cluster) %>%
  ggplot +
  geom_line(aes(x, Resistance, 
                group = obs, 
                lty = `wear level`,
                lwd = `wear level`
                )) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  xlab("Time [ms]") +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  theme(legend.position = "none") +
  scale_linetype_manual(values = c(1,2,3)) + 
  scale_size_manual(values = c(.3,.5,.6)) +
  ylim(c(260, 320)) +
  theme(aspect.ratio = 1)
plot_pad0
# ggsave("plots/plot_pad0.pdf", plot_pad0, width = 5, height = 4)

plot_pad <- df_obs %>% 
  mutate(cluster = factor(as.character(cluster))) %>% 
  left_join(pad_profiles %>% 
              rename(obs = pad_profiles) %>% 
              mutate(obs = as.character(obs)), by = c("obs", "wear.level")) %>% 
  rename("wear level" = "wear.level") %>% 
  ggplot +
  geom_line(aes(x, Resistance, col = cluster), lwd = 1,
            data = df_centr %>% 
              mutate(cluster = factor(as.character(cluster)))) +
  geom_line(aes(x, Resistance, group = obs, 
                col = cluster, 
                lty = `wear level`,
                lwd = `wear level`
  )) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  xlab("Time [ms]") +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
  scale_linetype_manual(values = 1:3) +
  scale_size_manual(values = c(.3,.5,.6)) +
  ylim(c(260, 320)) +
  theme(aspect.ratio = 1)
plot_pad
# ggsave("plots/plot_pad0_col.pdf", plot_pad, width = 6.5, height = 4)

plot_pad <- bind_rows(
  df_centr %>% mutate(type = "filtering B-spline hierarchical\ncentroids"),
  df_obs %>% mutate(type = "observations with\nwear information")
) %>% 
  mutate(cluster = ifelse(cluster == 1, "1/just renewed",
                          ifelse(cluster == 2, "2/medium wear", 
                                 "3/severe wear"))) %>%
  rename("cluster/wear level" = cluster) %>%
  ggplot +
  geom_line(aes(x, Resistance, group = obs, 
                col = `cluster/wear level`, 
                lty = type,
                lwd = type
                )) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  xlab("Time [ms]") +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))
         # lty = FALSE,
         # lwd = FALSE
         ) +
  scale_size_manual(values = c("filtering B-spline hierarchical\ncentroids" = 1,
                               "observations with\nwear information" = .6)) +
  ylim(c(260, 320))
  # facet_wrap(~ type)
plot_pad
# ggsave("plots/plot_pad_wrap.pdf", plot_pad, width = 9, height = 4)
# ggsave("plots/plot_pad_all.pdf", plot_pad, width = 6, height = 4)

#############################
# CFR multivariate --------------------------------------------------------
#############################
df_rev <- bind_cols(
  data$X %>%
    as.data.frame %>%
    mutate(x = 1:238) %>%
    pivot_longer(-x) %>% 
    filter(x < 50) %>% 
    group_by(name) %>% 
    summarise(min = min(value),
              which_min = x[which.min(value)]),
  
  data$X %>%
    as.data.frame %>%
    mutate(x = 1:238) %>%
    pivot_longer(-x) %>% 
    filter(x > 25) %>% 
    group_by(name) %>% 
    summarise(max = max(value),
              which_max = x[which.max(value)]) %>% 
    select(-name),
  
  data$X %>%
    as.data.frame %>%
    mutate(x = 1:238) %>%
    pivot_longer(-x) %>% 
    group_by(name) %>% 
    summarise(end = tail(value, 1)) %>% 
    select(-name)
) %>% 
  mutate(ampl_diff = max - min,
         phase_diff = which_max - which_min) %>% 
  mutate(name = factor(name, levels = colnames(data$X))) %>% 
  arrange(name)


# plot_analysis + 
#   geom_point(aes(which_min, min), data = df_rev, size = .25) +
#   geom_point(aes(which_max, max), data = df_rev, size = .25, col = "red") +
#   geom_point(aes(238, end), data = df_rev, size = .25, col = "darkgreen")

B <- df_rev %>% 
  select(ampl_diff, phase_diff, end)

mod_km <- NbClust(data = B,  
                  method = "kmeans",
                  min.nc = 2, 
                  max.nc = 10)
## Best number of clusters is 3
mod_hc <- NbClust(data = B,  
                  method = "ward.D2",min.nc = 2, 
                  max.nc = 10)
## Best number of clusters is 3

model_based <- mclustBIC(B,G=c(1,num_cluster_seq))
plot(model_based)
## Best number of clusters is 4

mod_opt_hc <- hcut(B, k = 3)
mod_opt_km <- kmeans(B, centers = 3)
mod_opt_mod <- Mclust(B, G = 4)

cl_hc <- factor(mod_opt_hc$cluster) %>% 
  fct_recode("1" = "2",
             "2" = "1") %>% 
  fct_relevel("1", "2", "3")
cl_hc
cl_km <- factor(mod_opt_km$cluster)
cl_mod <- factor(mod_opt_mod$classification) %>% 
  fct_recode("1" = "4",
             "2" = "1",
             "4" = "2") %>% 
  fct_relevel("1", "2", "3", "4")

df_features <- data$X %>%
  as.data.frame %>%
  mutate(x = 1:238) %>%
  pivot_longer(-x) %>% 
  inner_join(
    df_rev %>% 
      mutate(`features hierarchical` = cl_hc, `features k-means` = cl_km, `features model-based` = cl_mod),
    by = "name")

p <- df_features %>%  
  pivot_longer(c(`features hierarchical`, `features k-means`, `features model-based`), 
               names_to = "method", 
               values_to = "cluster") %>% 
  ggplot() +
  geom_line(aes(x, value, group = name, col = cluster), alpha = .25) +
  theme_bw() +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  xlab("Time [ms]") +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~method) +
  guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1)))
p
# ggsave("plots/features.pdf", p, width = 9, height = 3.5)  

# silhouette index --------------------------------------------------------
dist_l2 <- semimetric.basis(data$X_fd)

basis <- create.bspline.basis(c(0, 1), nbasis = 12)
loglam <- seq(-10, -4, 0.25)
Gcvsave <- numeric()
for (i in 1:length(loglam)) {
  fdPari     = fdPar(basis, Lfdobj = 2, 10^loglam[i])
  Sm.i       = smooth.basis(data$grid, 
                            data$X, 
                            fdPari)
  Gcvsave[i] = sum(Sm.i$gcv)
}
lambda_s=10^loglam[which.min(Gcvsave)]
fdPari  = fdPar(basis, Lfdobj=2,0)
X_fd<-smooth.basis(data$grid, 
                   data$X, 
                   fdPari)$fd

B_spline<-t(X_fd$coefs)
dimnames(B_spline) <- NULL

set.seed(0)
bootsrap_sil <- mclapply(1:25, function(ii) {
  rows <- sample(1:538, replace = TRUE)
  B <- df_rev %>% 
    slice(rows) %>% 
    select(ampl_diff, phase_diff, end)
  
  mod_opt_hc <- hcut(B, k = 3)
  mod_opt_km <- kmeans(B, centers = 3)
  mod_opt_mod <- Mclust(B, G = 4)
  
  cl_hc <- factor(mod_opt_hc$cluster)
  cl_km <- factor(mod_opt_km$cluster)
  cl_mod <- factor(mod_opt_mod$classification)
  
  mod_opt_hc_bspline <- hcut(B_spline[rows, ], k = 3)
  mod_opt_km_bspline <- kmeans(B_spline[rows, ], centers = 3)
  
  data.frame(
    functional_bspline_hc   = mean(silhouette(mod_opt_hc_bspline$cluster, dist_l2[rows, rows])[, 3]),
    functional_bspline_km   = mean(silhouette(mod_opt_km_bspline$cluster, dist_l2[rows, rows])[, 3]),
    multivariate_hc         = mean(silhouette(as.numeric(cl_hc), dist_l2[rows, rows])[, 3]),
    multivariate_km         = mean(silhouette(as.numeric(cl_km), dist_l2[rows, rows])[, 3]),
    multivariate_mb         = mean(silhouette(as.numeric(cl_mod), dist_l2[rows, rows])[, 3])
  )
}, mc.cores = detectCores()) %>% 
  bind_rows()

p <- bootsrap_sil %>% 
  rename("multivariate model-based" = "multivariate_mb",
         "multivariate hierarchical" = "multivariate_hc",
         "multivariate k-means" = "multivariate_km",
         "functional filtering B-spline k-means" = "functional_bspline_km",
         "functional filtering B-spline hierarchical" = "functional_bspline_hc",
  ) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "silhouette index") %>% 
  ggplot() +
  geom_boxplot(aes(x = `silhouette index`, y = method))
p
# ggsave("plots/multiv_vs_functional.pdf", p, width = 6, height = 3)

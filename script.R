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
  ylim(c(260, 320))
# plot_analysis
# ggsave("plots/plot_analysis.png", plot_analysis, width = 5, height = 4)

plot_ref_data <- read.csv("reference_profile.csv") %>% 
  mutate(time = 1:n()) %>%
  ggplot +
  geom_line(aes(time, Resistance)) +
  theme_bw() +
  ylab(expression("Resistance " ~ paste("[", mu, Omega, "]"))) +
  xlab("Time [ms]") +
  ylim(c(260, 320))
# plot_ref_data
# ggsave("plots/plot_ref_data.pdf", plot_ref_data, width = 5, height = 4)



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

get_centroids_df <- function(cl, method) {
  df <- data$X %>% 
    t %>% as.data.frame %>% 
    mutate(class = cl) %>% 
    group_by(class) %>% 
    summarise_all(mean) %>% 
    dplyr::select(-class) %>% as.matrix %>% t %>%
    as.data.frame %>% 
    setNames(names(table(cl))) %>% 
    mutate(x = seq(0, 1, length = n())) %>% 
    pivot_longer(-x, names_to = "cluster", values_to = "Resistance") %>% 
    mutate(mod = method) %>% 
    mutate(cluster = factor(cluster))
  levels(df$cluster) <- data.frame(ranked_clusters = df %>% 
                                     filter(x > .0885, x < .089) %>% 
                                     arrange(desc(Resistance)) %>% 
                                     pull(cluster)) %>% 
    mutate(clusters = 1:n()) %>% 
    arrange(ranked_clusters) %>% 
    pull(clusters)
  df
  
}
lay <- rbind(c(rep(1, 495), rep(NA, 105)),
             c(rep(2, 600)),
             c(rep(3, 495), rep(NA, 105)))
p_cen <- grid.arrange(
  bind_rows(
    get_centroids_df(cl_3, "adaptive\nfunHDDC"),
    get_centroids_df(cl_4, "distance-based"),
    get_centroids_df(cl_7_hc, "raw\nhierarchical"),
    get_centroids_df(cl_2, "adaptive\ncurvclust"),
  ) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot +
    geom_line(aes(x, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  bind_rows(
    get_centroids_df(cl_5_hc, "filtering B-spline\nhierarchical"),
    get_centroids_df(cl_5_km, "filtering B-spline\nk-means"),
    get_centroids_df(cl_6_hc, "filtering FPCA\nhierarchical"),
    get_centroids_df(cl_6_km, "filtering FPCA\nk-means"),
    get_centroids_df(cl_7_km, "raw\nk-means"),
  ) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot +
    geom_line(aes(x, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  bind_rows(
    get_centroids_df(cl_1, "adaptive\nfunclust"),
    get_centroids_df(cl_5_mb, "filtering B-spline\nmodel-based"),
    get_centroids_df(cl_6_mb, "filtering FPCA\nmodel-based"),
    get_centroids_df(cl_7_mb, "raw\nmodel-based")
  ) %>% 
    mutate(cluster = cluster %>% as.character %>% factor) %>% 
    ggplot +
    geom_line(aes(x, Resistance, col = cluster)) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  layout_matrix = lay
)
# ggsave("plots/centroids.pdf", p_cen, width = 10, height = 10)

get_clustered_funs_df <- function(cl, method) {
  
  df_centroids <- data$X %>% 
    t %>% as.data.frame %>% 
    mutate(class = cl) %>% 
    group_by(class) %>% 
    summarise_all(mean) %>% 
    dplyr::select(- class) %>% as.matrix %>% t %>%
    as.data.frame %>% 
    setNames(names(table(cl))) %>% 
    mutate(x = seq(0, 1, length = n())) %>% 
    pivot_longer(- x, names_to = "cluster", values_to = "Resistance") %>% 
    mutate(mod = method) %>% 
    mutate(cluster = factor(cluster))
  
  df <- data$X %>% 
    as.data.frame %>% 
    mutate(x = seq(0, 1, length = n())) %>% 
    pivot_longer(- x, names_to = "obs", values_to = "Resistance") %>% 
    inner_join(data.frame(obs = colnames(data$X),
                          cluster = cl), by = "obs") %>% 
    mutate(cluster = as.character(cluster),
           mod = method) %>% 
    mutate(cluster = factor(cluster))
  
  levels(df$cluster) <- data.frame(ranked_clusters = df_centroids %>% 
                                     filter(x > .0885, x < .089) %>% 
                                     arrange(desc(Resistance)) %>% 
                                     pull(cluster)) %>% 
    mutate(clusters = 1:n()) %>% 
    arrange(ranked_clusters) %>% 
    pull(clusters)
  
  df %>% 
    mutate(cluster = cluster %>% as.character %>% factor)
}

lay <- rbind(c(rep(1, 502), rep(NA, 98)),
             c(rep(2, 600)),
             c(rep(3, 502), rep(NA, 98)))
p_funs <- grid.arrange(
  
  bind_rows(
    get_clustered_funs_df(cl_3, "adaptive\nfunHDDC"),
    get_clustered_funs_df(cl_4, "distance-based"),
    get_clustered_funs_df(cl_7_hc, "raw\nhierarchical"),
    get_clustered_funs_df(cl_2, "adaptive\ncurvclust"),
  ) %>% 
    ggplot +
    geom_line(aes(x, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  
  bind_rows(
    get_clustered_funs_df(cl_5_hc, "filtering B-spline\nhierarchical"),
    get_clustered_funs_df(cl_5_km, "filtering B-spline\nk-means"),
    get_clustered_funs_df(cl_6_hc, "filtering FPCA\nhierarchical"),
    get_clustered_funs_df(cl_6_km, "filtering FPCA\nk-means"),
    get_clustered_funs_df(cl_7_km, "raw\nk-means"),
  ) %>% 
    ggplot +
    geom_line(aes(x, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  
  bind_rows(
    get_clustered_funs_df(cl_1, "adaptive\nfunclust"),
    get_clustered_funs_df(cl_5_mb, "filtering B-spline\nmodel-based"),
    get_clustered_funs_df(cl_6_mb, "filtering FPCA\nmodel-based"),
    get_clustered_funs_df(cl_7_mb, "raw\nmodel-based")
  ) %>% 
    ggplot +
    geom_line(aes(x, Resistance, group = obs, col = cluster), alpha = .25) +
    theme_bw() +
    facet_wrap(~ mod, nrow = 1) + 
    scale_color_brewer(palette = "Set1") +
    xlab("time") +
    ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
    guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))),
  layout_matrix = lay
)
# ggsave("plots/col_functions.png", p_funs, width = 10, height = 10)


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
  mutate(x = seq(0, 1, length = n())) %>% 
  pivot_longer(-x, names_to = "cluster", values_to = "Resistance") %>% 
  mutate(mod = "filtering B-spline hc") %>% 
  mutate(cluster = factor(cluster)) %>% 
  mutate(obs = paste("centroid", cluster))
ranked_clusters <- data.frame(ranked_clusters = df_centr %>% 
                                filter(x > .0885, x < .089) %>% 
                                arrange(desc(Resistance)) %>% 
                                pull(cluster)) %>% 
  mutate(clusters = 1:n()) %>% 
  arrange(ranked_clusters)
levels(df_centr$cluster) <- ranked_clusters$clusters

df_obs <- data$X[, pad_profiles$pad_profiles] %>% 
  as.data.frame %>% 
  mutate(x = seq(0, 1, length = n())) %>% 
  pivot_longer(- x, names_to = "obs", values_to = "Resistance") %>% 
  inner_join(data.frame(cluster = cl_5_hc, obs = as.character(seq_along(cl_5_hc))),
             by = "obs") %>%
  inner_join(data.frame(`wear level` = pad_profiles$wear.level,
             obs = as.character(pad_profiles$pad_profiles)),
             by = "obs") %>% 
  mutate(cluster = cluster %>% as.character %>% factor,
         type = "observation")
levels(df_obs$cluster) <- ranked_clusters$clusters

plot_pad <- df_obs %>% mutate(type = "wear level") %>% 
  mutate(cluster = ifelse(cluster == 1, "renewed",
                          ifelse(cluster == 2, "medium wear", "severe wear"))) %>% 
  mutate(cluster = factor(cluster, levels = c("renewed", "medium wear", "severe wear"))) %>% 
  rename("wear level" = cluster) %>%
  ggplot +
  geom_line(aes(x, Resistance, 
                group = obs, 
                lty = `wear level`
                )) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  xlab("time") +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  theme(legend.position = "none") +
  scale_linetype_manual(values = c(1,2,3)) +
  ylim(c(260, 320))
plot_pad
# ggsave("plots/plot_pad0.pdf", plot_pad, width = 5, height = 4)


plot_pad <- bind_rows(
  df_centr %>% mutate(type = "Filtering B-spline hierarchical\ncentroid"),
  df_obs %>% mutate(type = "observation")
) %>% 
  mutate(cluster = ifelse(cluster == 1, "1/renewed",
                          ifelse(cluster == 2, "2/medium wear", "3/severe wear"))) %>%
  rename("cluster/wear level" = cluster) %>%
  ggplot +
  geom_line(aes(x, Resistance, 
                group = obs, 
                col = `cluster/wear level`,
                lty = type, lwd = type
  )) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  xlab("time") +
  ylab(expression("Resistance "~paste("[",mu,Omega,"]"))) +
  guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
  scale_size_manual(values = c("Filtering B-spline hierarchical\ncentroid" = .5, observation = 1)) +
  ylim(c(260, 320)) 
  # facet_wrap(~ type)
plot_pad
# ggsave("plots/plot_pad_wrap.pdf", plot_pad, width = 9, height = 4)
# ggsave("plots/plot_pad_all.pdf", plot_pad, width = 6, height = 4)

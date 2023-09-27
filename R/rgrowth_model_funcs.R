# This script is a series of functions required to generate species parameters (e.g., population growth rate, carrying capacity, dispersal, etc.)
# The script can be run to generate a single species by commenting-in the lines indicated by "comment in"
# This script will not sample or degrade species abundances
# Leave untouched to source the functions in the rgrowth_model.R file.

library(tidyverse) # For manipulating tibbles and dataframes and plots
library(ggpubr) #for plotting (ggarrange)
library(truncnorm) #MSM and TruncatedNormal packages also available
library(data.table) #for rbindlist
library(pracma) # For moving average
library(tidyselect) # ...
library(Bchron) # Age-depth Modelling

### comment in ###
# source('force_funcs_1.1.r')
# walk_force <- gen.force.walk(time = 5000)
# lin_force <- gen.force.lin(time = 5000, start.val = 98, interval = 0.003)
# forces_ex <- list(walk_force, lin_force)
###

### Create a vector of disturbance occurrence for each time-step
disturbance_chance_func <- function(forces_df, disturbance_chance = 10){
  dist_chance <- ifelse(runif(lengths(forces_df[[1]])) < 1 / disturbance_chance, 1, 0)
}

### comment in ###
#disturbance_chance <- disturbance_chance_func(forces_ex)
###

### For each environmental driver, the parameters controlling the tolerance to each driver is generated
static_vals_func <- function(proxy_mean = 100, proxy_sd = 0, proxy_tol = 5, #parameters for distribution from which the mean is drawn
                             proxy_range_mean = 5, proxy_range_sd = 1, proxy_range_tol = 5, #parameters for distribution from which the sd is drawn
                             dist_tol = 5,
                             growth_rate_min = 1.01, growth_rate_max = 1.1){

  static_vals <- tibble(proxy_mean_a = ifelse(proxy_sd !=0, proxy_mean - (proxy_sd * proxy_tol), -Inf), #using above parameters to create lower bound truncation of mean, must be less than mean so mu is positive
                        proxy_mean_b = ifelse(proxy_sd !=0, proxy_mean + (proxy_sd * proxy_tol), Inf), #using above parameters to create upper bound truncation of mean, must be less than mean
                        #proxy_range_a = ifelse(proxy_range_sd !=0, proxy_range_mean - (proxy_range_sd * proxy_range_tol), -Inf), #using above parameters to create lower bound truncation of sd
                        proxy_range_a = 0, #Must be less than SD, SD must be positive.
                        proxy_range_b = ifelse(proxy_range_sd !=0, proxy_range_mean + (proxy_range_sd * proxy_range_tol), Inf), #using above parameters to create upper bound truncation of sd
                        #lambda_growth_rate = runif(1, min = growth_rate_min, max = growth_rate_max), #radomises per forcing if included in this func. Randomised per spp using offset_func
                        #dens = 300,
                        mu = rtruncnorm(1, a = proxy_mean_a, b = proxy_mean_b, mean = proxy_mean, sd = proxy_sd),
                        stdev = rtruncnorm(1, a = 0, b = proxy_range_b, mean = proxy_range_mean, sd = proxy_range_sd), #CAREFUL with stdev, MUST not be negative otherwis NAs produced.
                        dist_a = (mu - (stdev * dist_tol)),
                        dist_b = (mu + (stdev * dist_tol)),
                        mx_dnorm = dnorm(mu, mu, stdev))
}
### comment in ###
# static_vals <- apply(as.data.frame(forces_ex), MARGIN = 2, function(x){
#   static_vals_func(proxy_mean = 100, proxy_sd = 5, proxy_tol = 5, #parameters for distribution from which the mean is drawn
#                    proxy_range_mean = 10, proxy_range_sd = 2, proxy_range_tol = 5, #parameters for distribution from which the sd is drawn
#                    dist_tol = 5)})
# static_vals
###

### For each environmental driver, a maximum growth rate is generated.
growth_rate_func <- function(static_vals_lst, growth_rate_min = 1.05, growth_rate_max = 1.1){
  lambda_growth_rate <- runif(1:length(static_vals_lst), growth_rate_min, growth_rate_max)
  print(lambda_growth_rate)
  max_lambda_growth_rate <- prod(lambda_growth_rate) # maximum possible spp growth rate if all environmental drivers are at the optimum of their respective tolerances
  lambda_growth_rate <- sort(lambda_growth_rate, decreasing = TRUE)
  print(lambda_growth_rate)
  inds <- order(sapply(static_vals_lst, "[[", "stdev"))
  print(inds)
  static_vals_lst <- Map(cbind, static_vals_lst, lambda_growth_rate = lambda_growth_rate[order(inds)])
  ### slightly shorter script below achieves same. Maybe less readable?
  # static_vals_lst <- Map(cbind, static_vals_lst[order(unlist(lapply(static_vals_lst, `[[`, "stdev")))],
  #                        lambda_growth_rate = lambda_growth_rate[order(-lambda_growth_rate)])[names(static_vals_lst)]
  static_vals_lst <- lapply(static_vals_lst, mutate, max_growth_rate = max_lambda_growth_rate)
  static_vals_lst
}
### comment in ###
# static_vals <- growth_rate_func(static_vals)
# static_vals
###

########## Create single K per spp (not one per forcing)
k_func <- function(static_vals_lst, k_mean = 1000, k_sd = 100, k_tol = 5){
  k_df <- tibble(
    k_a = ifelse(k_sd !=0, k_mean - (k_sd * k_tol), -Inf),
    k_b = ifelse(k_sd !=0, k_mean + (k_sd * k_tol), Inf),
    k = rtruncnorm(1, a = k_a, b = k_b, mean = k_mean, sd = k_sd))
  k_df

  static_vals_lst <- lapply(static_vals_lst, function(x){
    x <- x %>%
      mutate(k = k_df$k)})
  static_vals_lst
}
### comment in ###
# static_vals <- k_func(static_vals)
# static_vals
###

########## calculate dnorm for each forcing per timestep
d_niche_func <- function(dynamic_df, static_df){
  dynamic_df[["d_niche"]] = dtruncnorm(dynamic_df[[1]], a = static_df[["dist_a"]], b = static_df[["dist_b"]],
                                 mean = static_df[["mu"]], sd = static_df[["stdev"]])
  dynamic_df
}
### comment in ###
# forces_ex <- Map(d_niche_func, forces_ex, static_vals)
# lapply(forces_ex, head)
###

### dnorm plot using randomised offset.
dnorm_plot_func <- function(static_vals_lst){
dnorm_dfs_static_vals_lst <- lapply(static_vals_lst, function(x){
  dnorm_df <- dnorm(seq(x$dist_a, x$dist_b, 0.1), x$mu, x$stdev) * (x$lambda_growth_rate / x$mx_dnorm )-1
  dnorm_df <- as.data.frame(dnorm_df) %>%
    mutate(niche_range = seq(x$dist_a, x$dist_b, 0.1),
           index = 1:n())

  positive_range <- range(dnorm_df[,2][dnorm_df[,1] > 0], na.rm = TRUE)
  x$lower_positive_limit <- positive_range[1]
  x$upper_positive_limit <- positive_range[2]
  list(dnorm_df, x)
  #print(positive_range)

})
#dnorm_df <- data.frame(lapply(dnorm_df, "length<-", max(lengths(dnorm_df))))  #this code binds a list of vectors of unequal lengths to a df and fills excess rows with NA
#list(dnorm_lst, static_vals_lst)
dnorm_dfs_static_vals_lst
}
### comment in ###
# dnorm_dfs_static_vals_lst <- dnorm_plot_func(static_vals)
# static_vals <- lapply(dnorm_dfs_static_vals_lst, function(x) {x <- x[[2]]})
# dnorm_df_lst <- lapply(dnorm_dfs_static_vals_lst, function(x) {x <- x[[1]]})
###


########## calculate environmental effect for each forcing per timestep
lambda_envif_func <- function(dynamic_df, static_df){
  dynamic_df[["lambda_envif"]] = (dynamic_df[["d_niche"]] * (static_df[["lambda_growth_rate"]] / static_df[["mx_dnorm"]]))
  dynamic_df
}
### comment in ###
# forces_ex <- Map(lambda_envif_func, forces_ex, static_vals) # used to have argument for min_dist_height
# lapply(forces_ex, head)
###

########## weight each forcing with vector wi. wi must = number of forces
wi_vals_func <- function(dynamic_df, wi = NA){

  if (any(is.na(wi))) {
    wi <- rep(1, length(dynamic_df))
  }

  weights <- (wi/(sum(wi))) * length(wi)
  if(length(weights) != length(dynamic_df)){stop(paste0("weights length not equal to forces"))}
  #weights <- as.list(weights)
  print("multiplier (sums to number of forces)")
  print(weights)
  # proportion_weights <- map_dbl(wi, function(x){(x / sum(wi))*100})
  # print("Percent of total weight")
  # print(proportion_weights)
  wi_envif <- map2(dynamic_df, weights, ~ ((.x$lambda_envif-1) * .y)+1)#Does lambda need to convert to r before weighting? YES
  dynamic_df <- map2(dynamic_df, wi_envif, ~ cbind(.x, .y))
  dynamic_df <- lapply(dynamic_df, function(x) {
    names(x)[names(x) == ".y"] <- "wi_envif"
    x
  })
  dynamic_df
}
### comment in ###
# forces_ex <- wi_vals_func(forces_ex)
# lapply(forces_ex, head)
###

########## calculate total environmental effect as the product of effect per forcing
tot_r_envif_func <-  function(lst_of_df) {
  lst_of_df <- lapply(lst_of_df, function(x) {
    x <- x %>%
      select(wi_envif)
  })
  df <- as.data.frame(lst_of_df)
  df <- df %>%
    mutate(tot_r_envif = apply(df, 1, prod)-1)
  #can use select here and only retain tot_envif - currently retaining all for checking calculations
}
### comment in ###
# tot_env_df <- tot_r_envif_func(forces_ex)
# head(tot_env_df)
##########

########## dispersal
# disp_func <- function(df, disp_chance, disp_mean, disp_sd){
#   df <- df %>%
#     mutate(disp_chance = sample(1:disp_chance, size = length(df[,1]), replace = TRUE),
#            disp_val = case_when(disp_chance == 1 ~ rnorm(length(df[,1]), disp_mean, disp_sd), TRUE ~ as.numeric(0)))
# }
#tot_env_df <- disp_func(df = tot_env_df, disp_chance = 50, disp_mean = 100, disp_sd = 10)

disp_func <- function(df, disp_chance, disp_min, disp_max){
  df <- df %>%
    mutate(disp_val = ifelse(runif(length(df[,1])) < 1 / disp_chance, runif(length(df[,1]), disp_min, disp_max), 0))
}
### comment in ###
# tot_env_df <- disp_func(df = tot_env_df, disp_chance = 1000, disp_min = 10, disp_max = 100)
###


##########
kapacity_func <- function(dynamic_df, static_vals_lst, k_offtake_min, k_offtake_max){
  dynamic_df <- dynamic_df %>%
    mutate(kapacity = static_vals_lst[[1]][["k"]] * runif(n(), k_offtake_min, k_offtake_max))
  dynamic_df
}
### comment in ###
# tot_env_df <- kapacity_func(tot_env_df, static_vals, k_offtake_min = 0.8, k_offtake_max = 1.2)
# head(tot_env_df)
###

########## generate a vector of disturbances (future version will use beta distribution)
disturbance_func <- function(dynamic_df, disturbance_occurance){

  magnitude <- 1 - rgamma(length(dynamic_df[, 1]), shape = 1, rate = 8)
  for (i in 1:length(magnitude)) {
    while (magnitude[i] <= 0 | magnitude[i] >= 1) {
      magnitude[i] <-  1 - rgamma(1, shape = 1, rate = 8)
    }
    magnitude
  }
  #magnitude

  dynamic_df <- dynamic_df %>%
    mutate(disturbance_occ = disturbance_occurance,
           disturbance = ifelse(disturbance_occ == 1, magnitude, 1))
  dynamic_df

}
### comment in ###
# tot_env_df <- disturbance_func(tot_env_df, disturbance_chance)
# head(tot_env_df, 20)
###

########## Generate a single species
##### r can be scaled with time, Lambda cannot Gotelli, primer in ecology pg 13

N_func <- function(df, N, static_vals_lst){
  df <- df %>%
    mutate(N = static_vals_lst[[1]][["k"]]-1)
  for (i in 2:length(df[,1])) {
    if(df$N[i-1] >= df$kapacity[i-1] & df$tot_r_envif[i-1] <= 0 & df$tot_r_envif[i-1] > -0.5){
      modifier = 1
    } else{
      if(df$N[i-1] >= df$kapacity[i-1] & df$tot_r_envif[i-1] <= -0.5){
        modifier = abs(df$tot_r_envif[i-1])*2
      } else {
        modifier = 0
      }
    }
    df$N[i] <- ((df$N[i-1] + (((modifier + df$tot_r_envif[i-1]) * df$N[i-1]) * (1-(df$N[i-1]/df$kapacity[i-1])))) * df$disturbance[i-1]) + df$disp_val[i-1]
  }
  df$N <- ifelse(df$N < 1, 0, df$N)
  if(any(df$N == df$kapacity)){stop("N == K. N for eternity!")}
  return(df)
}

##### comment in #####
# tot_env_df <- N_func(tot_env_df, 100, static_vals)
# head(tot_env_df, 20)
# tot_env_df_plot <- tot_env_df %>%
#   mutate(time = 1:length(tot_env_df[,1]))
# tot_env_df[is.na(tot_env_df$N)] <- 0
#
# N_plot <- ggplot(tot_env_df_plot, aes(x = time, y = N))+
#   geom_area(fill = "grey")+ #consider fill by relative abundance
#   coord_flip()+
#   #  facet_grid(.~proxies, scales = "free")+
#   #facet_grid(.~species)+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   theme(panel.grid.major = element_line(colour = "black"),
#         panel.grid.minor = element_blank()) +
#   xlab("Time")+
#   ylab("Abundance")
# N_plot
##########
#######################################

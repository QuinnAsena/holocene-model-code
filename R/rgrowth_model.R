# This script is the primary run function that calls functions from the other scripts to generate S species
# the S species are then degraded and sampled
# Example parameters have been provided for the script to run.
# The script (if untouched) will take between 13-15 minutes to run a single replicate
# Approximately 7 GB of data will be used in R's memory
# approximately 350 MB of data will be saved

list.of.packages <- c("testthis", "devtools", "tidyverse", "ggpubr", "truncnorm", "data.table", "pracma", "tidyselect", "Bchron", "rlist", "renv", "fst", "config")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(ggpubr) #for plotting (ggarrange)
library(truncnorm) #MSM and TruncatedNormal packages also available
library(data.table) #for rbindlist
library(pracma) # for moving average
library(tidyselect)
library(Bchron) #For age-depth
library(rlist)
library(renv)
library(fst)
library(config)

start_time <- Sys.time()

source('force_funcs_1.1.r')
source('force_scaling_funcs_1.1.R')
source('rgrowth_model_funcs.R')
source('degradation_funcs_1.1.R')


# Example input parameters may be called
other_params <- list(
primary_directory = "model_dir",
seed = 44000,
no_model_runs = 1,
no_species = 200,
chance_of_disturbance = 500,
chr_vec_of_forcings = c("walk", "sstep", "ppower", "linear", "llogistic", "ssin",
             "eexp", "ppulse", "nnoisenorm", "nnoiseuni", "uunimodal")
)

force_params_tib <- tibble(
time.step.param = 5000,
step.start.val.param = 98, step.flip.time.param = 3000, step.end.val.param = 102, step.stdev.param = 0,
walk.mu.param = 100, walk.stdev.param = 0.1,
power.onset.param = 1000, power.start.val.param = 1, power.interval.param = 0.001, power.gradient.param = 1, power.power.param = 2,
lin.start.val.param = 98, lin.interval.param = 0.0015, lin.end.val.param = time.step.param, lin.stdev.param = 0,
logistic.onset.param = 1000, logistic.r.param = 0.003, logistic.K.param = 100, logistic.N.param = 1, logistic.stdev.param = 0, logistic.high.param = 102, logistic.low.param = 98,
sin.onset.param = 1, sin.start.val.param = 2, sin.interval.param = 0.02, sin.amp.param = 0.05, sin.stdev.param = 0.005,
exp.start.val.param = 1, exp.onset.param = 750, exp.interval.param = 0.1, exp.stdev.param = 0,
pulse.start.val.param = 1, pulse.flip.time.param = 20, pulse.period.val.param = 5, pulse.flip.back.param = 40, pulse.end.val.param = 1, pulse.stdev.param = 0,
noisenorm.start.val.param = 100, noisenorm.interval.param = 0,  noisenorm.stdev.param = 0,
noiseuni.mn.param = 95, noiseuni.mx.param = 105,
unimodal.start.val.param = -10, unimodal.end.val.param = 10, unimodal.amp.param = 0.5, unimodal.stdev.param = 0.5,
force_weight = list(c(0.15,0.85))
)
# species params
species_params_tib <- tibble(
proxy_mu_param = 100, proxy_stdev_param = 20, proxy_tolerance_param = 5,
proxy_range_mu_param = 20, proxy_range_stdev_param = 1, proxy_range_tolerance_param = 5,
dist_tolerance_param = 5,
growth_rate_minimum_param = 1.01, growth_rate_maximum_param = 1.1,
k_mu_param = 4000, k_stdev_param = 500, k_tolerance_param = 5,
dispersal_chance_param = 100, dispersal_max_param = 10, dispersal_min_param = 1,
k_offtake_min_param = 0.8, k_offtake_max_param = 1)

#degradation params
degredation_params_list <- list(
mix_wind = c(2:10),
# core chr params
depth_power = 1, #5000^1.27 = 49,853 close to radiocarbon dating limit
# age-depth params
dating_res = c(500),

var_acc_mean = 0,
var_acc_sd = 1,
var_acc_low = 0.2,
var_acc_high = 2.5,
lin_acc_cm_yr = 0.2,
lin_acc_acc_rate = 0.00,
lin_acc_gradient = 1.0001,
sample_res = c(1:10),
thickness = 1,

count_res = c(seq(from = 100, to = 1000, by = 100))
)
#
model_settings <- list(other_params, force_params_tib, species_params_tib, degredation_params_list, other_params$seed)

primary_directory <- other_params$primary_directory
if (dir.exists(primary_directory) == FALSE) {
dir.create(paste0(primary_directory))
}

directory_func <- function(directory = "model_multirun_") {
  path <- paste0(primary_directory, "/", directory, format(Sys.time(), '%Y-%m-%d_%H-%M-%S'))
  dir.create(path)
  path
}
saved_location <- directory_func()
print(saved_location)

if (dir.exists(saved_location) == FALSE) {
  stop("No directory to save output")
} else {
  #dir.create(multirun_model_directory)
  model_runs <- vector("list", other_params$no_model_runs)
  names(model_runs) <- paste0("model_", 1:length(model_runs))
  #model_runs_list <- lapply(model_runs, function(output){
}

for (i in 1:length(model_runs)) {
  dir.create(paste0(saved_location, "/", "model_run_",  i, collapse = ""))

# Generate the driving environment according to parameter list
force.function <- function(forcings = other_params$chr_vec_of_forcings) {
  forces_result <- c()
  for (force_name in forcings) {
    force <- switch(
      force_name,
      step = gen.force.step(time = force_params_tib$time.step.param, start.val = force_params_tib$step.start.val.param, flip.time = force_params_tib$step.flip.time.param, end.val = force_params_tib$step.end.val.param, stdev = force_params_tib$step.stdev.param),
      walk = gen.force.walk(time = force_params_tib$time.step.param, mu = force_params_tib$walk.mu.param, stdev = force_params_tib$walk.stdev.param),
      power = gen.force.power(time = force_params_tib$time.step.param, onset = force_params_tib$power.onset.param, start.val = force_params_tib$power.start.val.param, interval = force_params_tib$power.interval.param, gradient = force_params_tib$power.gradient.param, power = force_params_tib$power.power.param),
      linear = gen.force.lin(time = force_params_tib$time.step.param, start.val = force_params_tib$lin.start.val.param, interval = force_params_tib$lin.interval.param, end.val = force_params_tib$lin.end.val.param, stdev = force_params_tib$lin.stdev.param),
      logistic = gen.force.logistic(time = force_params_tib$time.step.param, onset = force_params_tib$logistic.onset.param, r = force_params_tib$logistic.r.param, K = force_params_tib$logistic.K.param, N = force_params_tib$logistic.N.param, stdev = force_params_tib$logistic.stdev.param,  high = force_params_tib$logistic.high.param, low = force_params_tib$logistic.low.param),
      sin = gen.force.sin(time = force_params_tib$time.step.param, onset = force_params_tib$sin.onset.param, start.val = force_params_tib$sin.start.val.param, interval = force_params_tib$sin.interval.param, amp = force_params_tib$sin.amp.param, stdev = force_params_tib$sin.stdev.param),
      exp = gen.force.exponential(time = force_params_tib$time.step.param, start.val = force_params_tib$exp.start.val.param, onset = force_params_tib$exp.onset.param, interval = force_params_tib$exp.interval.param, stdev = force_params_tib$exp.stdev.param),
      pulse = gen.force.pulse(time = force_params_tib$time.step.param, start.val = force_params_tib$pulse.start.val.param, flip.time = force_params_tib$pulse.flip.time.param, period.val = force_params_tib$pulse.period.val.param, flip.back = force_params_tib$pulse.flip.back.param, end.val = force_params_tib$pulse.end.val.param, stdev = force_params_tib$pulse.stdev.param),
      noisenorm = gen.force.noisenorm(time = force_params_tib$time.step.param, start.val = force_params_tib$noisenorm.start.val.param, interval = force_params_tib$noisenorm.interval.param, stdev = force_params_tib$noisenorm.stdev.param),
      noiseuni = gen.force.noiseuni(time = force_params_tib$time.step.param, mn = force_params_tib$noiseuni.mn.param, mx = force_params_tib$noiseuni.mx.param),
      unimodal = gen.force.uni(time = force_params_tib$time.step.param, start.val = force_params_tib$unimodal.start.val.param, end.val = force_params_tib$unimodal.end.val.param, amp = force_params_tib$unimodal.amp.param, stdev = force_params_tib$unimodal.stdev.param)
    )

    if (class(force) == "data.frame") {
      forces_result = append(forces_result, list(force))
    } else {
      warning(paste0(force_name, " is not defined"))
    }
  }
  forces_result
}
forces_ex <- force.function(other_params$chr_vec_of_forcings)
class(forces_ex)
lapply(forces_ex, head)
lapply(forces_ex, tail)
matplot(as.data.frame(forces_ex), type = 'l', main = "unscaled forces")

# Generate disturbances
disturbance_chance <- disturbance_chance_func(forces_ex, disturbance_chance = other_params$chance_of_disturbance)

# Generate species
spp_func <- function(lst_of_forces) # N_func CURRENTLY USING N = K. N_1 NOT USED AS ARGUMENT.
{
  static_vals <- apply(as.data.frame(lst_of_forces), MARGIN = 2, function(x){
    static_vals_func(proxy_mean = species_params_tib$proxy_mu_param, proxy_sd = species_params_tib$proxy_stdev_param, proxy_tol = species_params_tib$proxy_tolerance_param, #parameters for distribution from which the mean is drawn
                     proxy_range_mean = species_params_tib$proxy_range_mu_param, proxy_range_sd = species_params_tib$proxy_range_stdev_param, proxy_range_tol = species_params_tib$proxy_range_tolerance_param, #parameters for distribution from which the sd is drawn
                     dist_tol = species_params_tib$dist_tolerance_param,
                     growth_rate_min = species_params_tib$growth_rate_minimum_param, growth_rate_max = species_params_tib$growth_rate_maximum_param
                     )})
  static_vals <- growth_rate_func(static_vals, growth_rate_min = species_params_tib$growth_rate_minimum_param, growth_rate_max = species_params_tib$growth_rate_maximum_param)
  static_vals <- k_func(static_vals, k_mean = species_params_tib$k_mu_param, k_sd = species_params_tib$k_stdev_param, k_tol = species_params_tib$k_tolerance_param)
  dnorm_dfs_static_vals_lst <- dnorm_plot_func(static_vals)
  static_vals <- lapply(dnorm_dfs_static_vals_lst, function(x) {x <- x[[2]]})
  dnorm_df_lst_plot <- lapply(dnorm_dfs_static_vals_lst, function(x) {x <- x[[1]]})

  lst_of_forces <- Map(d_niche_func, lst_of_forces, static_vals)
  lst_of_forces <- Map(lambda_envif_func, lst_of_forces, static_vals)
  lst_of_forces <- wi_vals_func(lst_of_forces, wi = force_params_tib$force_weight[[1]])
  print(lapply(lst_of_forces, head))
  ##########
  tot_env_df <- tot_r_envif_func(lst_of_forces)
  tot_env_df <- disp_func(df = tot_env_df, disp_chance = species_params_tib$dispersal_chance_param, disp_max = species_params_tib$dispersal_max_param, disp_min = species_params_tib$dispersal_min_param)
  tot_env_df <- kapacity_func(tot_env_df, static_vals, species_params_tib$k_offtake_min_param, species_params_tib$k_offtake_max_param)
  tot_env_df <- disturbance_func(tot_env_df, disturbance_chance)
  tot_env_df <- N_func(tot_env_df, static_vals_lst = static_vals)
  print(head(tot_env_df))
  list(tot_env_df, static_vals, dnorm_df_lst_plot)
}
### For EACH proxy, model outputs 1 time-series df of species population, 1 df of static values, and 1 df of niche distribution per forcing
n_spp_list <- vector("list", other_params$no_species)
names(n_spp_list) <- paste0("spp_", 1:length(n_spp_list))
n_spp_list <- lapply(n_spp_list, function(x){spp_func(forces_ex)})

### saves a list of all species dfs only for the rest of the model
species_list <- lapply(n_spp_list, function(x){
  spp_lst <- x[[1]]
})

# Burn-in now done at analysis stage
# species_list <- burn_in_func(species_list, other_params$burn_in)

### saves a list of all species static values dfs only for diagnosis
species_static_list <- lapply(n_spp_list, function(x){
  species_static_list <- x[[2]]
})
### saves a list of all species niche dfs only for plotting
dnorm_list <- lapply(n_spp_list, function(x){
  spp_lst <- x[[3]]
})

dnorm_df_lst <- lapply(dnorm_list, function(x){
  x <- rbindlist(x, idcol = "forcing")
})
lapply(dnorm_df_lst, head)
lapply(dnorm_df_lst, tail)


N_list <- lapply(species_list, function(x) {
  x <- x %>%
    select("N")
})
N_df <- unlist(N_list, recursive = F)
N_df <- bind_cols(N_df)
# Burn-in now done at analysis stage
#N_df <- burn_in_func(N_df, other_params$burn_in)

# Generate core
var_acc_vec <- variable_acc_func(N_df, mu = degredation_params_list$var_acc_mean, stdev = degredation_params_list$var_acc_sd, low = degredation_params_list$var_acc_low, high = degredation_params_list$var_acc_high)
lin_acc_vec <- lin_acc_func(N_df, cm_yr = degredation_params_list$lin_acc_cm_yr, acc_rate = degredation_params_list$lin_acc_acc_rate, gradient = degredation_params_list$lin_acc_gradient)
lin_var_acc_vec <- rowMeans(cbind(lin_acc_vec, var_acc_vec))
core_charactaristics <- core_func(spp_df = N_df, acc_type = "pre_def", pre_acc_vec = lin_var_acc_vec)
head(core_charactaristics)
tail(core_charactaristics)

##### Sampling and degradation of species
N_df_mix <- mix_run_func(spp_df = N_df, wind = degredation_params_list$mix_wind)
N_df_mix_multi_sam_multi <- lapply(N_df_mix, sample_depth_func_run_func, core_chr_df = core_charactaristics, thick = degredation_params_list$thickness, freq = degredation_params_list$sample_res)
N_df_mix_multi_sam_multi <- unlist(N_df_mix_multi_sam_multi, recursive = F)
names(N_df_mix_multi_sam_multi) <- gsub(pattern = "\\.", replacement = "", x = names(N_df_mix_multi_sam_multi))

N_df_mix_multi_count_multi <- lapply(N_df_mix, count_run_func, cnt = degredation_params_list$count_res)
N_df_mix_multi_count_multi <- unlist(N_df_mix_multi_count_multi, recursive = F)

N_df_mix_multi_count_multi <- lapply(N_df_mix_multi_count_multi, function(x){
  cbind(sim_time = 1:nrow(x), x)
})

names(N_df_mix_multi_count_multi) <- gsub(pattern = "\\.", replacement = "", x = names(N_df_mix_multi_count_multi))

N_df_mix_multi_sam_multi_count_multi <- lapply(N_df_mix_multi_sam_multi, function(x) {
  x <- count_run_func(x[,3:ncol(x)], cnt = degredation_params_list$count_res)
})
N_df_mix_multi_sam_multi_count_multi <- unlist(N_df_mix_multi_sam_multi_count_multi, recursive = F)
names(N_df_mix_multi_sam_multi_count_multi) <- gsub(pattern = "\\.", replacement = "", x = names(N_df_mix_multi_sam_multi_count_multi))

N_df_mix_multi_sam_multi_count_multi <- lapply(N_df_mix_multi_sam_multi_count_multi, function(x){
  cbind(N_df_mix_multi_sam_multi[[1]][,1:2], x)
})

# Age-depth modelling
age_depth <- age_depth_func(core_chr_df = core_charactaristics, res = degredation_params_list$dating_res)
end_time <- Sys.time()
end_time - start_time


age_depth_cut <- lapply(seq_along(1:3), function(i){
  depth_ind <- rev(N_df_mix_multi_sam_multi[[1]]$depth_groups_sam_1)
  age_depth[[i]] <- age_depth[[i]][depth_ind,]
})


# Various checksums from model testing
check_sum_mix <- sapply(N_df_mix, sum, na.rm = T)
print(sprintf("Check sums for raw and mixed dataframes: %.12g", check_sum_mix))


if (length(N_df_mix_multi_sam_multi) < 10) {
  check_sum_sam_dfs <- 1:length(N_df_mix_multi_sam_multi)
} else {
  check_sum_sam_dfs <- seq(1, length(N_df_mix_multi_sam_multi), length.out = 10)
}
check_sum_mix_sam <- sapply(N_df_mix_multi_sam_multi[check_sum_sam_dfs], sum, na.rm = T)
print(sprintf("Check sums for mixed and sampled dataframes: %.12g", check_sum_mix_sam))



if (length(N_df_mix_multi_count_multi) < 10) {
  check_sum_count_dfs <- 1:length(N_df_mix_multi_count_multi)
} else {
  check_sum_count_dfs <- seq(1, length(N_df_mix_multi_count_multi), length.out = 10)
}
check_sum_mix_count <- sapply(N_df_mix_multi_count_multi[check_sum_count_dfs], sum, na.rm = T)
print(sprintf("Check sums for mixed and counted dataframes: %.12g", check_sum_mix_count))



if (length(N_df_mix_multi_sam_multi_count_multi) < 10) {
  check_sum_mix_sam_count_dfs <- 1:length(N_df_mix_multi_sam_multi_count_multi)
} else {
  check_sum_mix_sam_count_dfs <- seq(1, length(N_df_mix_multi_sam_multi_count_multi), length.out = 10)
}
check_sum_mix_sam_count_dfs <- sapply(N_df_mix_multi_sam_multi_count_multi[check_sum_mix_sam_count_dfs], sum, na.rm = T)
print(sprintf("Check sums for mixed, sampled and counted dataframes: %.12g", check_sum_mix_sam_count_dfs))


absolute_check_sum <- sum(check_sum_mix, check_sum_mix_sam, check_sum_mix_count, check_sum_mix_sam_count_dfs)
print(sprintf("Absolute check sum (sum of all check sums) %.12g", absolute_check_sum))


end_time <- Sys.time()
time_taken <- end_time - start_time
print(sprintf("Took %.3f secs or mins to run.", time_taken))

#model_runs[[i]] <- list(species_list, species_static_list, dnorm_df_lst, species_list_mix, core_charactaristics, age_depth, species_list_mix_sam, species_df_mix_sam_cnt, forces_ex, saved_location, model_settings)
model_runs[[i]] <-
  list(
    species_list = species_list,
    species_static_list = species_static_list,
    dnorm_df_lst = dnorm_df_lst,
    N_df = N_df,
    N_df_mix = N_df_mix,
    N_df_mix_multi_sam_multi = N_df_mix_multi_sam_multi,
    N_df_mix_multi_count_multi = N_df_mix_multi_count_multi,
    N_df_mix_multi_sam_multi_count_multi = N_df_mix_multi_sam_multi_count_multi,
    core_charactaristics = core_charactaristics,
    age_depth = age_depth,
    age_depth_cut = age_depth_cut,
    forces_ex = forces_ex,
    saved_location = saved_location,
    model_settings = model_settings,
    time_taken = time_taken
  )



saveRDS(model_runs, file = paste0(saved_location, "/", "model_run_",  i, "/", "all_output.rds", collapse = ""))

write_fst(N_df, path = paste0(saved_location, "/", "model_run_",  i, "/", "N_df.fst", collapse = ""))
write_fst(core_charactaristics, path = paste0(saved_location, "/", "model_run_",  i, "/", "core_charactaristics.fst", collapse = ""))
write_fst(age_depth_cut[[1]], path = paste0(saved_location, "/", "model_run_",  i, "/", "age_depth_cut_agesonly.fst", collapse = ""))

saveRDS(
  list(
    N_df = N_df,
    N_df_mix = N_df_mix,
    N_df_mix_multi_sam_multi = N_df_mix_multi_sam_multi,
    N_df_mix_multi_count_multi = N_df_mix_multi_count_multi,
    N_df_mix_multi_sam_multi_count_multi = N_df_mix_multi_sam_multi_count_multi
  ),
  file = paste0(saved_location, "/", "species_output.rds", collapse = "")
)


# if (dir.exists(degredation_params_list$core_directory)) {
# file.copy(from = degredation_params_list$core_directory, to =  paste0(saved_location, "/", "model_run_",  i, collapse = ""), recursive = TRUE)
# }
  model_runs
}



# This script is a set of functions to sample and degrade the simulated data
# functions include simulated mixing, sub-sampling, proxy counting and age-depth modelling

list.of.packages <- c("tidyverse", "Bchron", "mgcv", "data.table", "DescTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(Bchron)
library(mgcv)
library(data.table)
library(DescTools)

##### mixing function
## Function to mix data using a rolling average over a given window width
mix_func <- function(spp_df, wind = c(2)){
  mixed_df <- apply(spp_df, MARGIN = 2, function(x){
    movavg(x, n = wind, type = "w")
  })
  mixed_df <- as_tibble(mixed_df)
  colnames(mixed_df) <- paste0(colnames(mixed_df), "_mix_", wind)
  mixed_df
}

### mixing wrapper function
## Function to apply mixing function at different window widths
mix_run_func <- function(wind = c(2, 3), spp_df = N_df) {
  results_list <- mapply(
    mix_func,
    wind = wind,
    MoreArgs = list(spp_df = spp_df),
    SIMPLIFY = F
  )
  names(results_list) <- paste0("N_df_mix_", wind)
  results_list <- append(list(spp_df), results_list)
  names(results_list)[1] <- paste0("N_df")
  results_list
}

### Function to simulate counting individuals across specie abundances.
count_func <-  function(spp_df, cnt) {
  spp_names <- names(spp_df)
  spp_mat <- as.matrix(spp_df)
  nc <- ncol(spp_mat)
  spp_cols <- 1:nc
  for (i in 1:nrow(spp_df)) {
    spp_row <- spp_mat[i,]
    if(all(is.na(spp_row)) == FALSE && sum(spp_row) > cnt) {
      # randomly sample columns
      r <- sample(spp_cols, cnt, replace=T, prob=spp_row)
      # count the number of occurrences of a column in the
      # sampling list. Fill in with zero if need be.
      spp_mat[i,] <- tabulate(r, nc)
    }
  }
  spp_df = as.data.frame(spp_mat)
  names(spp_df) <- paste0(spp_names, "_cnt_", cnt)
  return(spp_df)
}

## Applies count function over multiple resolutions
count_run_func <- function(cnt = c(300, 400), spp_df = N_df) {
  results_list <- mapply(count_func,
                         cnt = cnt,
                         MoreArgs = list(
                           spp_df = spp_df
                         ),
                         SIMPLIFY = F

  )
  names(results_list) <- paste0("_cnt_", cnt)
  results_list
}

# Converts data to relative abundances
relative_abund_func <- function(spp_df) {
  spp_df[is.na(spp_df)] <- 0 # Required so that if row contains a single NA value. Row does not equal NA
  #print(spp_df)
  spp_df <- (spp_df / rowSums(spp_df))*100 # will return NaN for rows that equal 0, NA, or NaN
  spp_df[is.na(spp_df)] <- 0
  spp_df
}

## Accumulation and age-depth model functions
# cm_yr should be 0 for constant sed rate and < 0 for
# If a constant acc rate is desired then acc_rate = 0 and gradient = 1
# If a change by a constant proportion each time-step then acc_rate = 0 and gradient = the desired rate of increase/decrese. Be careful as a negative gradient of say 0.5 will halve the starting input
# If a change by a fixed number per tine-step is desired then acc_rate = eg 0.1.

### Starting from bottom of core:
## Linear accumuation rate
lin_acc_func <- function(spp_df, cm_yr = 0.2, acc_rate = 0.00, gradient = 1.0001) {
  acc <- seq(from = cm_yr, by = acc_rate, length.out = nrow(spp_df))
  if(gradient != 1) {
    acc[2] <- acc[2] * gradient
    for (i in 3:length(acc)) {
      acc[i] <- acc[i-1] * gradient
    }
  }
  acc[length(acc)] <- 0
  if(any(acc < 0)){
    stop("zero or negative accumulation rate")
  }
  return(acc)
}

# Variable accumulation rate using a smoothed random walk
variable_acc_func <- function(spp_df, mu = 0, stdev = 1, low = 0.2, high = 2.5) {
  walk_acc <- cumsum(rnorm(nrow(spp_df), mu, stdev))
  # smooth
  time <-  1:length(walk_acc)
  acc_smooth <- gam(walk_acc ~ s(time))
  walk_acc <- acc_smooth$fitted.values
  # re-scale
  walk_acc <- ((high - low) * (walk_acc - min(walk_acc)) / (max(walk_acc) - min(walk_acc))) + low
  walk_acc[length(walk_acc)] <- 0
  return(walk_acc)
}

# To avoid doing decimal to a power (0.2^2) the sequence starts with 1 and is re-scaled
# Exponential accumulation rate function
exp_acc_func <- function(spp_df, cm_yr = 1, acc_rate = 0.002, power = 2, low = 0.2, high = 1) {
  acc <- (seq(from = cm_yr, by = acc_rate, length.out = nrow(spp_df)))^power
  acc <- ((high - low) * (acc - min(acc)) / (max(acc) - min(acc))) + low
  acc[length(acc)] <- 0
  return(acc)
}

#####
# function to combine everything into a dataframe of accumulation rate and simulation time
core_func <- function(spp_df, acc_type, pre_acc_vec){
  if(acc_type == "linear") {
    acc_vec <- lin_acc_func(spp_df)
  } else {
    if(acc_type == "variable"){
      acc_vec <- variable_acc_func(spp_df)
    } else {
      if(acc_type == "lin_var"){
        acc_vec1 <- lin_acc_func(spp_df)
        acc_vec2 <- variable_acc_func(spp_df)
        acc_vec <- rowMeans(cbind(acc_vec1, acc_vec2))
        # acc_vec <- (acc_vec1+acc_vec2) /2
      } else {
        if(acc_type == "exponential") {
         acc_vec <- exp_acc_func(spp_df)
        } else{
          if(acc_type == "exp_var") {
            acc_vec1 <- exp_acc_func(spp_df)
            acc_vec2 <- variable_acc_func(spp_df)
            acc_vec <- rowMeans(cbind(acc_vec1, acc_vec2))
          } else {
            if(acc_type == "pre_def") {acc_vec <- pre_acc_vec}
          }
        }
      }
    }
  }
  df <- tibble(
    time_age = nrow(spp_df):1,
    lin_acc_depth = nrow(spp_df):1,
    error = round(rnorm(nrow(spp_df), 61, 8.67)),
    labID = paste("id_", nrow(spp_df):1),
    sim_time = 1:nrow(spp_df),
    sim_acc = acc_vec,
    cum_depth = rev(cumsum(rev(sim_acc))), #Reversing for depth from top of core
    cum_depth_round = round(cum_depth),
    sim_acc_diff = 0
  )
  df$sim_acc_diff[(nrow(df)-1):1] <- diff(df$cum_depth[(nrow(df)):1]) / diff((df$sim_time[(nrow(df)):1]))
  df$error[length(df$error)] <- 1
  df$labID[length(df$labID)] <- "top_1"
  return(df)
}


# Function to sample core at given depth intervals (depth intervals are even)
sample_depth_func <- function(core_chr_df, spp_df, freq = 5, thick = 1) {

  core_chr_df <- core_chr_df %>%
    mutate(depth_groups = cut(core_chr_df$cum_depth, breaks = seq(from = 0, to = ceiling(core_chr_df$cum_depth[1]), by=thick), include.lowest = T, labels = F))
  print(head(core_chr_df, 10))
  print(tail(core_chr_df, 10))

  unik <- !duplicated(rev(core_chr_df$depth_groups))

  time_sliced <- core_chr_df$time_age[seq_along(core_chr_df$time_age)[rev(unik)]]
  print(head(time_sliced))
  print(tail(time_sliced))

  spp_df <- cbind(depth_groups = core_chr_df$depth_groups, spp_df)
  print(head(spp_df[,1:10]))
  print(tail(spp_df[,1:10]))
  print(dim(spp_df))
  spp_df <- aggregate(.~ depth_groups, spp_df, sum, simplify = T) ### CARFUL, AGGREGATE RE-ORDERS AS ASCENDING
  spp_df <- spp_df[order(spp_df$depth_groups, decreasing = T), ]
  colnames(spp_df) <- paste0(colnames(spp_df), "_sam_", freq)


  freq_vec <- seq(from = 1, to = nrow(spp_df), by = freq)
  print(head(freq_vec))
  #spp_df <- spp_df[freq_vec, ]
  spp_df[-freq_vec,] <- NA # Fill everything else with NA so that dim matches age-depth output dim
  spp_df <- cbind(time_sliced, spp_df)

  print(dim(spp_df))

  print(head(spp_df[,1:10]))
  print(tail(spp_df[,1:10]))
  print(dim(spp_df))

  spp_df

}

##### Function to map over multiple chosen sampling intervals
sample_depth_func_run_func <- function(core_chr_df = core_chr, spp_df = N_df, freq = c(2, 3), thick = 1) {
  results_list <- mapply(
    sample_depth_func,
    freq = freq,
    MoreArgs = list(spp_df = spp_df,
                    core_chr_df = core_chr_df,
                    thick = 1),
    SIMPLIFY = F
  )
  names(results_list) <- paste0("_sam_", freq)
  results_list
}


### Age-depth model function to sample simulated core at gen depths and generate estimated ages per simulation time-step
age_depth_func <- function(core_chr_df, res = 200) {
  if (res == 0) {
    return(core_chr_df)
  }
  ind_func <- function(from = 0, to = core_chr_df$cum_depth[1], by = res)
  {
    vec <- do.call(what = seq, args = list(from, to, by))
    if (tail(vec, 1) != to) {
      return(c(vec, to))
    } else {
      return(vec)
    }
  }
  ind <- ind_func()
  print(ind)

  ind <- sapply(ind, function(x) {
    Closest(core_chr_df$cum_depth, x, which = T)
  })
  print(ind)

  chron_df <- core_chr_df[ind,]
  print(chron_df)
  print(core_chr_df)
  chron_df$calCurves <-  c("normal", rep("shcal13", nrow(chron_df)-1) )
  chron_df <- chron_df[order(chron_df$cum_depth), ]
  print(chron_df)

  chronOut = with(chron_df,
                  Bchronology(ages=time_age,
                              ageSds=error,
                              calCurves=calCurves,
                              positions=cum_depth,
                              positionThicknesses = rep(0, nrow(chron_df)),
                              ids=labID,
                              predictPositions = seq(0, core_chr_df$cum_depth[1], 1)))

  ages <- summary(chronOut)
  colnames(ages) <- paste0("ages_", colnames(ages))

  acc <- summary(chronOut, type = "acc_rate")
  colnames(acc) <- paste0("acc_", colnames(acc))

  sed <- summary(chronOut, type = "sed_rate", useExisting = F)
  colnames(acc) <- paste0("sed_", colnames(acc))
  output <- list(ages = ages, acc_rate = acc, sed_rate = sed, chron_out = chronOut)

}

# Mapply agee depth function over multile resolutions
age_depth_run_func <- function(core_chr_df, res = c(100, 500)) {
  results_list <- mapply(age_depth_func,
                         res = res,
                         MoreArgs = list(
                           core_chr_df = core_chr_df
                         ),
                         SIMPLIFY = F
  )
}

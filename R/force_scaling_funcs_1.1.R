## This script is a simple scaling function
## these functions are used to scale the environmental drivers to a desired mean/range
## Drivers may need to be scaled to ensure their values remain within a reasonable range with respect to species optima/tolerance

###### As two functions
gen_midrange <- function(mu_mean = 100, s_dev = 1, tol = 5) {
  new_midrange <- rnorm(1, mu_mean, s_dev)
  while (new_midrange < (mu_mean - (s_dev * tol)) |
         new_midrange > (mu_mean + (s_dev * tol))) {
    new_midrange <- rnorm(1, mu_mean, s_dev)
  }
#  cat(paste(new_midrange, '\n'))
  new_midrange
}

rescale_force <- function(force_mat, mu = 100, stdev = 1, tolerance = 5){
    {
    scaling_factor <- apply(force_mat, MARGIN = 2, function(x){gen_midrange(mu_mean = mu, s_dev =  stdev, tol = tolerance)/ ((max(x)+min(x))/2)})
    scaling_factor
    scaled_force <- data.frame(mapply('*', as.data.frame(force_mat), scaling_factor))
    scaled_force
  }
}

##########################################
#### As single function, can input df or mat, will output df
rescale_force1 <- function(force_mat, mu_mean = 100, s_dev = 20, tol = 5){
  gen_midrange <- function(){
    new_midrange <- rnorm(1, mu_mean, s_dev)
    while (new_midrange < (mu_mean - (s_dev*tol)) | new_midrange > (mu_mean + (s_dev*tol))) {
      new_midrange <- rnorm(1, mu_mean, s_dev)
    }
    cat(paste(new_midrange, '\n'))
    new_midrange
  }
  {
    scaling_factor <- apply(force_mat, MARGIN = 2, function(x){gen_midrange()/ ((max(x)+min(x))/2)})
    scaling_factor
    scaled_force <- data.frame(mapply('*', as.data.frame(force_mat), scaling_factor))
    scaled_force
  }
}

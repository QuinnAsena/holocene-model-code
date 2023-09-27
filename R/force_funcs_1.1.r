# These functions generate patterns that can be used for environmental drivers in the primary rgrowth_model script
# Function inputs need to be tuned (or scaled using force_scaling_funcs_1.1.R) with respect to species optima/tolerance

library("tidyverse")

# Generate an abrupt shift at a given time-step
gen.force.step <- function(df = "YES", n = 1, type = "step", rand = "norm", start.val = 60, end.val = 62, flip.time = 50, time = 100, stdev = 0)
{
  step.series <- vector(mode = "numeric", length = time) ##length should equal n or time? Both run but may have differences for n > 1
  if (type == "step") {
    step.series <- vapply(1:n, function(x) {
      step.series[1:flip.time] <- start.val
      step.series[(flip.time + 1):time] <- end.val
      step.series
    }, FUN.VALUE = numeric(length(step.series)))
    colnames(step.series) <- paste0("step_", 1:n)
    step.series # have to call step.series here
    print(class(step.series))

    if (rand == "norm"){
      step.series <- apply(step.series, MARGIN = 2, function(x){
        x <- rnorm(x, mean = x, stdev)
      })
      print(class(step.series))
    }

      if(df == "YES"){
        print(class(step.series))
        step.series <- as.data.frame(step.series)
        print(class(step.series))
      }
  }
    step.series
}
#gen.force.step()
# step_force <- gen.force.step()
# step_force
# matplot(step_force, type = "l")
#matplot(gen.force.step(), type = 'l')
########################################## PULSE ##########################################
# Generate a 'pulse' in driver effect
gen.force.pulse <- function(df = "YES", n = 1, type = "pulse", rand = "norm", start.val = 1, flip.time = 20, period.val = 5, flip.back = 40, end.val = 1, time = 100, stdev = 0)
{
  pulse.series <- vector(mode = "numeric", length = time)
  if (type == "pulse")
  {
    pulse.series <- vapply(1:n, function(x){
    pulse.series[1:flip.time] <- start.val
    pulse.series[(flip.time+1):flip.back] <- period.val
    pulse.series[flip.back:time] <- end.val
    print(pulse.series)
    }, FUN.VALUE = numeric(length(pulse.series)))
    colnames(pulse.series) <- paste0("pulse_", 1:n)
    pulse.series # have to call pulse.series here

    if (rand == "norm"){
      pulse.series <- apply(pulse.series, MARGIN = 2, function(x){
        x <- rnorm(x, mean = x, stdev)
      })
    }
      if(df == "YES"){
        print(class(pulse.series))
        pulse.series <- as.data.frame(pulse.series)
        print(class(pulse.series))
      }
    }
    pulse.series
}



########################################## SIN ##########################################
# generate a wave pattern with a given period and frequency
gen.force.sin <- function(df = "YES", n = 1, type = "sin", rand = "norm", start.val = 2, time = 100, interval = 0.1, amp = 5, stdev = 0, onset = 15) #What eq gives control of frequeny?
{
 sin.series <- vector(mode = "numeric", length = time)
  if (type == "sin")
  {
    sin.series <- vapply(1:n, function(x){
    #sin.series <- seq(from = start.val, length = time, by = interval)
    sin.series[1:onset] <- start.val
    sin.series[(onset+1):time] <- seq(from = start.val, length = (time-onset), by = interval)
    #print(sin.series)
    #print(length(sin.series[(onset+1):time]))
    sin.series <- amp * sin(sin.series)                     #### This give sin of start.val:onset
    sin.series <- (sin.series + abs(min(sin.series))+1) ### shift to positive values required for rescale function
    }, FUN.VALUE = numeric(length(sin.series)))
    colnames(sin.series) <- paste0("sin_", 1:n)
    sin.series # have to call sin.series here

    if (rand == "norm"){
      sin.series <- apply(sin.series, MARGIN = 2, function(x){
        x <- rnorm(x, mean = x, stdev)
      })
    }
      if(df == "YES"){
        print(class(sin.series))
        sin.series <- as.data.frame(sin.series)
        print(class(sin.series))
      }
    }
    sin.series
}

########################################## LINEAR ##########################################
# an increasing or decreasing presure at a constant rate
gen.force.lin <- function(df = "YES", n = 1, type = "linear", rand = "norm", start.val = 0, interval = 1,
                          end.val = time, time = 100, stdev = 0)
{
  lin.series <- vector(mode = "numeric", length = time)
  if (type == "linear")
  {
    lin.series <- vapply(1:n, function(x){
    lin.series <- seq(from = start.val, by = interval, length = time)
#    lin.series <- gradient*lin.series
    #lin.series[lin.series >= end.val] <- end.val
    lin.series[end.val:time] <- lin.series[end.val]
    lin.series <- lin.series
    }, FUN.VALUE = numeric(length(lin.series)))
colnames(lin.series) <- paste0("lin_", 1:n)
lin.series # have to call lin.series here

if (rand == "norm"){
  lin.series <- apply(lin.series, MARGIN = 2, function(x){
    x <- rnorm(x, mean = x, stdev)
  })
}
  if(df == "YES"){
    print(class(lin.series))
    lin.series <- as.data.frame(lin.series)
    print(class(lin.series))
  }
}
lin.series
}

# linear_force <- gen.force.lin(start.val = 98.09800, interval = 0.003003, time = 5000)
# matplot(linear_force, type = 'l')
#matplot(gen.force.lin(), type = 'l')
#plot(l, type="l", xlab = "", ylab = "", xaxt="n", yaxt="n", tck=0, lwd=4, box(lwd=4))

###

### A sigmoidal shaped shift in conditions with a given rate and onset
gen.force.logistic <- function(df = "YES", n = 1, type = "logistic", rand = "norm", onset = 30, r = 0.2, K = 100, N = 0.1, time = 100, stdev = 0, rescale = "TRUE", high = 104, low = 96)
{
  logistic.series <- vector(length = time, mode = "numeric")
  if (type == "logistic")
  {
    logistic.series <- vapply(1:n, function(x)
    {
      logistic.series[1:onset] <- N
#      logistic.series[onset+1] <- N+r*N*(1-N/K)
      for (i in (onset+1):time)
      {
        logistic.series[i] <- logistic.series[i-1] + (r * logistic.series[i-1]*(1-logistic.series[i-1]/K))
      }
      logistic.series
    }, FUN.VALUE = numeric(length(logistic.series)))
    colnames(logistic.series) <- paste0("logistic_", 1:n)
  }

  if (rand == "norm"){
    logistic.series <- apply(logistic.series, MARGIN = 2, function(x){
      x <- rnorm(x, mean = x, stdev)
    })
  }
  if(rescale){
    logistic.series <- ((high - low) * (logistic.series - min(logistic.series)) / (max(logistic.series) - min(logistic.series))) + low
  }

  if(df == "YES"){
    print(class(logistic.series))
    logistic.series <- as.data.frame(logistic.series)
    print(class(logistic.series))
    print(logistic.series)
  }

  logistic.series
}
# logistic_force <- gen.force.logistic()
# matplot(logistic_force, type = "l")
#matplot(gen.force.logistic(), type = 'l')


########################################## Exponential ##########################################
# Exponentially increasing pressure. Tends to be too powerful!
gen.force.exponential <- function(df = "YES", n = 1, rand = "norm", type = "exponential", onset = 60, start.val = 1, interval = 0.1, time = 100, stdev = 0)
  {
    exp.series <- vector(mode = "numeric", length = time)
    if (type == "exponential")
    {
      exp.series <- vapply(1:n, function(x) {
        exp.series[1:onset] <- start.val
        exp.series[(onset+1):time] <- seq(from = start.val, by = interval, length = (time-onset))
        exp.series <- exp(exp.series)
      }, FUN.VALUE = numeric(length(exp.series)))
      colnames(exp.series) <- paste0("exp_", 1:n)
      exp.series
      print(exp.series)
      if (rand == "norm") {
        exp.series <- apply(exp.series, MARGIN = 2, function(x){
          x <- rnorm(x, mean = x, stdev)
        })
      }
        if(df == "YES"){
          print(class(exp.series))
          exp.series <- as.data.frame(exp.series)
          print(class(exp.series))
        }
      }
      exp.series
    }
# exp_force <- gen.force.exponential()
# matplot(exp_force, type = 'l')
#matplot(gen.force.exponential(), type = 'l')

################################################
## Ignore this function. Exponential driver pressure way to much!
# gen_force_super_duper_mega_exponential <- function(type = "exponential", start.val = 1, time = 100, mu = 0, sd = 0)
# {
#   exp_series <- vector(length = time, mode = "numeric")
#   if (type == "exponential")
#
#     for (i in (start.val+1):time) {
#       exp_series[i] <- exp(exp_series[i-1])
#     }
#   exp_series
# }
# super_duper_mega_exponential_force <- gen_force_super_duper_mega_exponential()
# plot(super_duper_mega_exponential_force, type = 'l')


########################################## Power ##########################################
###### increasing at a given power: y = ax^b
gen.force.power <- function(df = "YES", n = 1, rand = "norm", type = "power", start.val = 1, interval = 0.1, onset = 10, gradient = 1, time = 100, power = 2, stdev = 0)
{
  power.series <- vector(mode = "numeric", length = time)
  if (type == "power")
  {
    power.series <- vapply(1:n, function(x)
    {
      power.series[1:onset] <- start.val
      power.series[(onset+1):time] <- seq(from = start.val, by = interval, length = (time-onset))
      power.series[(onset+1):time] <- (power.series[(onset+1):time]*gradient)^power

      power.series
    }, FUN.VALUE = numeric(length(power.series)))
    colnames(power.series) <- paste0("power_", 1:n)
    power.series

    if (rand == "norm") {
      power.series <- apply(power.series, MARGIN = 2, function(x){
        x <- rnorm(x, mean = x, stdev)
      })
    }
      if(df == "YES"){
        print(class(power.series))
        power.series <- as.data.frame(power.series)
        print(class(power.series))
      }
    }
    power.series
}
# power_forc <- gen.force.power()
# matplot(power_forc, type="l")
#matplot(gen.force.power(), type = 'l')

######  Experimenting with different equations 1^b, 2^b, 3^b...
gen.force.power1 <- function(df = "YES", n = 1, rand = "norm", type = "power1", onset = 10, start.val = 1, interval = 0.1, time = 100, power = 5, stdev = 0)
{
  power.series <- vector(mode = "numeric", length = time)
  if (type == "power1")
  {
    power.series <- vapply(1:n, function(x)
    {
      power.series[1:onset] <- start.val
      power.series[(onset+1):time] <- seq(from = start.val, by = interval, length = (time-onset))
      print(power.series)
      power.series[(onset+1):time] <- power.series[(onset+1):time]^power
      print(power.series)
    }, FUN.VALUE = numeric(length(power.series)))
    colnames(power.series) <- paste0("power_", 1:n)
    power.series

    if (rand == "norm") {
      power.series <- apply(power.series, MARGIN = 2, function(x){
        x <- rnorm(x, mean = x, stdev)
      })
    }
      if(df == "YES"){
        print(class(power.series))
        power.series <- as.data.frame(power.series)
        print(class(power.series))
      }
    }
    power.series
}

# power_forc <- gen.force.power1()
# matplot(power_forc, type="l")
#matplot(gen.force.power1(), type = 'l')


########################################## NOISEUNI ##########################################
# Uniform (white) noise driver
gen.force.noiseuni <- function(df = "YES", n = 1, type = "noiseuni", mn = 0.7, mx = 1, time = 100)
{
  noiseuni.series <- vector(mode = "numeric", length = time)
  if (type == "noiseuni")
  {
    noiseuni.series <- vapply(1:n, function(x)
      {
      #noiseuni.series <- seq(from = start.val, by = interval, length = time)
      noiseuni.series <- runif(noiseuni.series, min = mn, max = mx)
      #runif(noiseuni.series, min = mn, max = mx)
    }, FUN.VALUE = numeric(length(noiseuni.series)))
    colnames(noiseuni.series) <- paste0("noiseuni_", 1:n)
    noiseuni.series
  }
  if(df == "YES"){
    print(class(noiseuni.series))
    noiseuni.series <- as.data.frame(noiseuni.series)
    print(class(noiseuni.series))
  }
  noiseuni.series
}

# uniform_force <- gen.force.noiseuni()
# matplot(uniform_force, type="l")
#matplot(gen.force.noiseuni(), type = 'l')



########################################## NOISENORM ##########################################
# Gaussian noise driver
gen.force.noisenorm <- function(df = "YES", n = 1, start.val = 1, interval = 0, type = "noisenorm", stdev = 0.08, time = 100)
{
  noisenorm.series <- vector(mode = "numeric", length = time)
  if (type == "noisenorm")
  {
    noisenorm.series <- vapply(1:n, function(x)
    {
      noisenorm.series <- seq(from = start.val, by = interval, length = time)
      #noisenorm.series <- rnorm(noisenorm.series, mean = noisenorm.series, sd = stdev)
      rnorm(noisenorm.series, mean = noisenorm.series, sd = stdev) #stdev dows not scale with increasing series
    }, FUN.VALUE = numeric(length(noisenorm.series)))
    colnames(noisenorm.series) <- paste0("runf_", 1:n)
    noisenorm.series
  }
  if(df == "YES"){
    print(class(noisenorm.series))
    noisenorm.series <- as.data.frame(noisenorm.series)
    print(class(noisenorm.series))
  }
  noisenorm.series
}

# noisenormform_force <- gen.force.noisenorm()
# matplot(noisenormform_force, type = "l")
#matplot(gen.force.noisenorm(), type = 'l')


########################################## WALK  ##########################################
# Random walk (red) noise driver
gen.force.walk <- function(df = "YES", n = 1, type = "walk", time = 100, mu = 100, stdev = 0.15)
{
  walk.series <- vector(length = time, mode = "numeric")
  if (type == "walk")
  {
    walk.series <- vapply(1:n, function(x)
    {
      walk.series[1] <- rnorm(1, mu, stdev)
      for (i in 2:time)
      {
        walk.series[i] <- rnorm(1, walk.series[i-1], stdev)
      }
      walk.series
    }, FUN.VALUE = numeric(length(walk.series)))
    colnames(walk.series) <- paste0("walk_", 1:n)
  }
  print(class(walk.series))
  if(df == "YES"){
    print(class(walk.series))
    walk.series <- as.data.frame(walk.series)
    print(class(walk.series))
  }
  walk.series
}
# walk_force <- gen.force.walk()
# matplot(walk_force, type = "l")
#matplot(gen.force.walk(), type = 'l')


########################################## unimodal  ##########################################
# A 'U' shaped driver
gen.force.uni <- function(df = "YES", n = 1, rand = "norm", type = "unimodal", start.val = -10, end.val = 10, amp = 0.5, time = 100, stdev = 0)
{
  uni.series <- vector(length = time, mode = "numeric")
  if (type == "unimodal")
  {
  uni.series <- vapply(1:n, function(x)
    {
    uni.series <- seq(from = start.val, to = end.val, length = time)
    uni.series <- -uni.series^2
    uni.series <- amp * uni.series
    uni.series <- (uni.series + abs(min(uni.series))+1)
    }, FUN.VALUE = numeric(length(uni.series)))
  colnames(uni.series) <- paste0("uni_", 1:n)
  uni.series

  if (rand == "norm") {
    uni.series <- apply(uni.series, MARGIN = 2, function(x){
      x <- rnorm(x, mean = x, stdev)
    })
    if(df == "YES"){
      print(class(uni.series))
      uni.series <- as.data.frame(uni.series)
      print(class(uni.series))
    }
  }
  uni.series

  }
  uni.series
}
# uni_force <- gen.force.uni()
# matplot(uni_force, type = "l")
#matplot(gen.force.uni(), type = 'l')



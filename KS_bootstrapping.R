set.seed(1)
num_obs = 100
obs = rnorm(num_obs) # eg. data at each reef 
mean(obs)  # e.g. probabilities at each reef in sector

est = mean(sample(obs, num_obs)) # sample one lot of data and gets one estimate

#repeats 1000 times, so you get 1000 estimates
num_resamples = 1000
vec = c()
for(i in 1:num_resamples){
  vec = c(vec, mean(sample(obs, num_obs, replace = TRUE)))
}

quantile(vec, c(0.025, 0.975))

#rough rough example

years = 2001:2021
distyears <- c(2003, 2004, 2011)
sample(years, length(distyears))
# estimate your probability from these observations
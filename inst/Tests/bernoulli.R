bernoulli.obs <- simulate_loggrowth(growth = 0.8, carry.cap = 500, movement = 0.25, sigma = 0.1,
                                    initial.pop = 200, initial.range = 0.3, initial.sigma=0.01, 
                                    timesteps = 4,sample.type = "Bernoulli", npoints = 100, obs.prob = 0.01,
                                    boundaries = c(0,1), debug = T)
ggplot()+
  gg(bernoulli.obs$animal_obs, aes(colour = obs))+
  facet_wrap(~time)

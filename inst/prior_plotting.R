library(tidyverse)
k <- seq(0,2000, length.out = 1001)
k.df <- data.frame(k = rep(k,times = 5), shape = rep(c(0.5,1,1.5,2,2.5), each = length(k)))

k.df %>% 
  mutate(rate = shape/1000) %>%
  mutate(value = dgamma(k, shape, rate))%>%
  ggplot()+
  geom_line(aes(x = k, y = value, colour = as.factor(shape)))



k.ln <- data.frame(k = rep(k,times = 5), sd = rep(c(0.2,0.4,0.6,0.8,1), each = length(k)))


k.ln %>% mutate(mu = log(1000)) %>%
  mutate(value = dnorm(log(k), mu, sd))%>%
  ggplot()+
  geom_line(aes(x = k, y = value, colour = as.factor(sd)))


ggplot()+
  geom_line(aes(x = k, y = dnorm(log(k), log(1000)-0.5*0.25*0.25, 0.25)))
r <- seq(0,1.5, by = 0.01)
r.ln <- data.frame(r = rep(r, times = 5), sd = rep(c(0.1,0.2,0.25,0.5,0.75), each = length(r)))
r.ln %>% mutate(mu = log(0.8)-0.5*sd*sd) %>%
  mutate(value = dnorm(log(r), mu, sd))%>%
  ggplot()+
  geom_line(aes(x = r, y = value, colour = as.factor(sd)))

ggplot()+
  geom_line(aes(x = k, y = dnorm(log(k), log(477)-0.5*0.4**2, 0.4)))+
  ylab("Density")+xlab("Carrying Capacity")


ggplot()+
  geom_line(aes(x = k, y = dgamma(k, 1, 1/488)))+
  xlab("Carrying Capacity")+ylab("Density")
ggsave("gammaprior.png")


k.prior <- data.frame(value = 700:1300,parameter = "k")
k.prior <- mutate(k.prior, density = dnorm(log(value), log(1000)-0.5*0.1*0.1, 0.1))
r.prior <- data.frame(value = seq(0.4,1.2,by = 0.01),parameter = "growth")
r.prior <- mutate(r.prior, density = dnorm(log(value), log(0.8) - 0.5*0.1*0.1, 0.1))
move.prior <- data.frame(value = seq(-1,1,by = 0.01),parameter = "move")
move.prior <- mutate(move.prior, density = dnorm(value, 0.1))
sigma.prior <- data.frame(value = seq(0,0.25,by = 0.01), parameter = "sigma")
sigma.prior <- mutate(sigma.prior, density = dnorm(log(value), log(0.1) - 0.5*0.1*0.1, 0.1))


priors.df <- rbind(k.prior, r.prior)#, move.prior, sigma.prior)
ggplot(priors.df)+
  geom_line(aes(x = value, y = density))+
  facet_wrap(~parameter, scales = "free")
ggsave("priors.png")
# SPDE --------------------------------------------------------------------
dens_prior_range = function(rho_0, p_alpha, upper)
{
  # compute the density of the PC prior for the
  # range rho of the Matern field
  # rho_0 and p_alpha are defined such that
  # P(rho<rho_0) = p_alpha
  rho = seq(0, upper, length.out =100)
  alpha1_tilde = -log(p_alpha) * rho_0
  dens_rho =  alpha1_tilde / rho^2 * exp(-alpha1_tilde / rho)
  return(data.frame(x = rho, y = dens_rho))
}

dens_prior_sd = function(sigma_0, p_sigma)
{
  # compute the density of the PC prior for the
  # sd sigma of the Matern field
  # sigma_0 and p_sigma are defined such that
  # P(sigma>sigma_0) = p_sigma
  sigma = seq(0, 0.1, length.out =100)
  alpha2_tilde = -log(p_sigma)/sigma_0
  dens_sigma = alpha2_tilde* exp(-alpha2_tilde * sigma) 
  return(data.frame(x = sigma, y = dens_sigma))
}

ggplot()+
  geom_line(data = dens_prior_range(0.05,0.05,0.6), aes(x,y), colour = "red")+
  geom_line(data = dens_prior_range(0.05,0.01,0.6), aes(x,y), colour = "green")+
  geom_line(data = dens_prior_range(0.1,0.05,0.6), aes(x,y), colour = "blue")

ggplot()+
  geom_line(data = dens_prior_sd(0.25,0.1), aes(x,y))

ggplot()+
  geom_line(data = dens_prior_sd(0.1,0.1), aes(x,y, colour = "0.1"))+
  geom_line(data = dens_prior_sd(0.11,0.1), aes(x,y, colour = "0.125"))+
  geom_line(data = dens_prior_sd(0.1,0.15), aes(x,y, colour = "0.15"))+
  geom_vline(xintercept = 0.015)
  #geom_line(data = dens_prior_sd(0.15,0.1), aes(x,y, colour = "0.15"))

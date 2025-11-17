library(tidyverse)
k <- seq(0,1000)
k.df <- data.frame(k = rep(k,times = 5), shape = rep(c(0.5,1,1.5,2,2.5), each = length(k)))

k.df %>% 
  mutate(rate = shape/500) %>%
  mutate(value = dgamma(k, shape, rate))%>%
  ggplot()+
  geom_line(aes(x = k, y = value, colour = as.factor(shape)))


lognorm <- function(x,mu,sd){
  return(exp(-((log(x)-mu)**2)/(2*sd**2))/(x*sd*sqrt(2*pi)))
}

k.ln <- data.frame(k = rep(k,times = 10), sd = rep(seq(0.1,1,by = 0.1), each = length(k)))
k.ln %>% mutate(mu = log(500)-0.5*sd*sd) %>%
  mutate(value = lognorm(k, mu, sd))%>%
  ggplot()+
  geom_line(aes(x = k, y = value, colour = as.factor(sd)))

k.ln %>% mutate(mu = log(500)-0.5*sd*sd) %>%
  mutate(value = dnorm(log(k), mu, sd))%>%
  ggplot()+
  geom_line(aes(x = k, y = value, colour = as.factor(sd)))


r <- seq(0,1.5, by = 0.01)
r.ln <- data.frame(r = rep(r, times = 5), sd = rep(c(0.1,0.25,0.5,1,2), each = length(r)))
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
  sigma = seq(0, sigma_0*10, length.out =100)
  alpha2_tilde = -log(p_sigma)/sigma_0
  dens_sigma = alpha2_tilde* exp(-alpha2_tilde * sigma) 
  return(data.frame(x = sigma, y = dens_sigma))
}

ggplot()+
  geom_line(data = dens_prior_range(0.025,0.01,0.6), aes(x,y))+
  geom_line(data = dens_prior_range(0.05,0.05,0.6), aes(x,y), colour = "red")+
  geom_line(data = dens_prior_range(0.05,0.01,0.6), aes(x,y), colour = "green")

ggplot()+
  geom_line(data = dens_prior_sd(0.25,0.1), aes(x,y))

ggplot()+
  geom_line(data = dens_prior_sd(1,0.1), aes(x,y, colour = "1"))+
  geom_line(data = dens_prior_sd(0.5,0.1), aes(x,y, colour = "0.5"))

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
  geom_line(aes(x = k, y = dnorm(log(k), log(493)-0.5**3, 0.5)))+
  ylab("Density")+xlab("Carrying Capacity")

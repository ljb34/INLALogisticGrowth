remotes::install_github("ljb34/INLAlogisticgrowth", force = F)
library(INLA)
library(inlabru)
library(INLAloggrowth)
library(sf)
library(dplyr)
#Simulate data
out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 1000, movement = 0.15, sigma = 0.1,
                               initial.pop = 500, initial.range = 0.2, initial.sigma=0.01,
                               timesteps = 4,sample.type = "LGCP", boundaries = c(0,1), debug = T,
                               max.edge = 0.1)
#Create meshes
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                       max.edge = c(0.1,0.25), offset = c(-0.1, 0.75))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)

#Fit initial model to get starting values
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.25, 0.1),
                              prior.range = c(0.05, 0.025))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear", mean.linear = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))-1

fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "cp",samplers = bnd)

#Find fitted values on mesh nodes
index <- min(which(stringr::str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))

#Get initial precision
n.nodes <- fit0$misc$configs$nconfig
nodes <- data.frame(log.prob=rep(NA,n.nodes))
mat_list <- list()
for(i in 1:n.nodes){
  nodes[i,]<- fit0$misc$configs$config[[i]]$log.posterior
  mat_list[[i]] <- fit0$misc$configs$config[[i]]$Q[1:mesh_obs$n,1:mesh_obs$n]
}
nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
  dplyr::mutate(weight.prob = weight/sum(weight))
Q <- Reduce("+", Map(function(m,w) w*m, mat_list,nodes$weight.prob))

initial.precision <-Diagonal(mesh_obs$n,
                             1/((fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))%*%Q

#Estimate starting linearisation point
initial.linpoint <- log(logit.nest(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]),
                                   0.8, 1000,4)$x)

priors <- list(cc = c(log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,]))-0.5*0.15*0.15,0.15),
               growth = c(log(0.8)-0.5*0.15*0.15,0.15),move = c(0.15,0.15),sigma = c(log(0.1),0.15))

##Warning! Takes about an hour to run
iterated.fit.lgcp <- iterate.cgeneric.fit.lgcp(data = out.lgcp$animal_obs, 
                                               smesh = mesh_obs, tmesh = mesh_time,
                                               samplers = bnd,
                                               prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                               prior.precision = initial.precision, 
                                               priors = priors,
                                               max.iter = 10, gamma = 0.4,
                                               stop.crit = 0.01,
                                               initial.linpoint = initial.linpoint, 
                                               initial.growth = log(0.8),
                                               initial.log.sigma = log(0.1),
                                               initial.move.const = 0.15,
                                               initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),
                                               saveall = T, options = list(verbose = T))
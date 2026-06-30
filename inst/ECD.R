library(data.table)
library(sf)
library(tidyverse)
library(stringr)
library(INLA)
library(inlabru)
library(fmesher)
library(INLAloggrowth)

bbs_dir <- "Eurasian Collared Dove/States/States"
files <- list.files(
  bbs_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
bbs <- rbindlist(
  lapply(files, fread),
  fill = TRUE
)

ecd <- bbs %>%
  filter(
    Year >= 1989,
    Year <= 2003,
    AOU == 22860
  )

routes <- read.csv("Eurasian Collared Dove/Routes.csv")

ecd_sf <- ecd %>%
  left_join(
    routes,
    by = c(
      "StateNum" = "StateNum",
      "Route" = "Route"
    )
  ) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)%>%
  st_transform(crs = 5070)
routes_with_ecd <- ecd_sf %>%
  st_drop_geometry() %>%
  distinct(StateNum, Route)

routes_subset <- routes %>% st_as_sf(coords = c("Longitude","Latitude"), crs = 4326)%>%
  st_transform(crs = 5070) %>%
  semi_join(routes_with_ecd, by = c("StateNum", "Route"))
years <- 1989:2003

route_years <- routes_subset %>%
  distinct(StateNum, Route, geometry) %>%
  tidyr::crossing(Year = years)
ecd_sf_full <- route_years %>%
  left_join(
    st_drop_geometry(ecd_sf),
    by = c("StateNum", "Route", "Year")
  ) %>%
  mutate(
    SpeciesTotal = replace_na(SpeciesTotal, 0)
  ) %>%
  st_as_sf()%>%
  st_transform(crs = 5070)
ecd_sf_1997 <- ecd_sf_full %>% mutate(time = Year -1997) %>% filter(time >= 1)

hull <- fm_nonconvex_hull(ecd_sf_1997, convex = -0.05)
mesh_obs <- fm_mesh_2d(loc = ecd_sf_1997, 
                       boundary = hull,
                       max.edge = c(100e3,1000e3 ), cutoff = 50e3, crs = 5070)
ggplot()+gg(mesh_obs)+gg(ecd_sf_1997)
mesh_obs$n
matern <- inla.spde2.pcmatern(
  mesh_obs,
  prior.range = c(150e3, 0.01),  
  prior.sigma = c(10, 0.05)
)

cmp <-  SpeciesTotal ~ smooth(geometry, model = matern) + Intercept(1)
fit1 <-  bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == 1,], family = "poisson",
             domain = list(geometry = fm_subdivide(mesh_obs, 2)), options = list(control.compute = list(dic = TRUE),
                                                                                 control.inla=list(control.vb=list(emergency=30))))
fit2 <- bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == 1,], family = "zeroinflatedpoisson1",
            domain = list(geometry = fm_subdivide(mesh_obs, 2)), options = list(control.compute = list(dic = TRUE),
                                                                                control.inla=list(control.vb=list(emergency=30))))
fit3 <- bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == 1,], family = "zeroinflatedpoisson0",
            domain = list(geometry = fm_subdivide(mesh_obs, 2)), options = list(control.compute = list(dic = TRUE),
                                                                                control.inla=list(control.vb=list(emergency=30))))
fit4 <- bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == 1,], family = "zeroinflatednbinomial1",
            domain = list(geometry = fm_subdivide(mesh_obs, 2)), options = list(control.compute = list(dic = TRUE),
                                                                                control.inla=list(control.vb=list(emergency=30))))
fit5 <- bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == 1,], family = "zeroinflatednbinomial0",
            domain = list(geometry = fm_subdivide(mesh_obs, 2)), options = list(control.compute = list(dic = TRUE),
                                                                                control.inla=list(control.vb=list(emergency=30))))

deltaIC(fit1, fit2, fit3, fit4, fit5) #zero inflated poisson is best family

fit <- fit2
mesh_locs <- st_as_sf(data.frame(x = mesh_obs$loc[,1], y = mesh_obs$loc[,2]), coords = c("x", "y"), crs = 5070)
pred0 <- predict(fit, mesh_locs, ~data.frame(matern = smooth + Intercept))

index <- min(which(stringr::str_sub(rownames(fit$summary.fitted.values),8,8)!= "A"))
n.nodes <- fit$misc$configs$nconfig
nodes <- data.frame(log.prob=rep(NA,n.nodes))
mat_list <- list()
#print(summary(fit0))
for(i in 1:n.nodes){
  nodes[i,]<- fit$misc$configs$config[[i]]$log.posterior
  Q <- fit$misc$configs$config[[i]]$Q[1:mesh_obs$n, 1:mesh_obs$n]
  dQ <- diag(Q)
  Q <- Q + Matrix::t(Q)
  diag(Q) <- dQ
  mat_list[[i]] <- Q
}

nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
  dplyr::mutate(weight.prob = weight/sum(weight))
Q <- Reduce("+", Map(function(m,w) w*m, mat_list,nodes$weight.prob))

initial.precision <-Diagonal(mesh_obs$n,
                             1/((fit$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))%*%Q
initial.linpoint <- pred0$mean
for(i in 2:6){
  fiti <- bru(cmp, data = ecd_sf_1997[ecd_sf_1997$time == i,], family = "zeroinflatedpoisson1",
              domain = list(geometry = fm_subdivide(mesh_obs, 2)))
  predi <- predict(fiti, mesh_locs, ~smooth + Intercept)
  initial.linpoint <- c(initial.linpoint, predi$mean)
  print(summary(fiti))
}

precomputes <- list(data = ecd_sf_1997, initial.precision = initial.precision, initial.linpoint = initial.linpoint,
             prior.mean = pred0$mean, mesh = mesh_obs, bnd = hull)

mesh_time <- fm_mesh_1d(1:6)

priors <- list(cc = c(log(10),5),
               growth = c(log(1.5),5),move = c(log(5e4),5),sigma = c(log(2),5)) 

cmp <- SpeciesTotal ~ -1
fit <- iterate.fit.custom(formula = cmp, family = "zeroinflatedpoisson1", data = precomputes$data, smesh=precomputes$mesh, tmesh=mesh_time, samplers = precomputes$bnd,
                          prior.mean = precomputes$prior.mean, gamma = 0.5, stop.crit = 0.1,
                          max.iter = 5,
                          prior.precision = precomputes$initial.precision,
                          priors = priors, initial.linpoint = precomputes$initial.linpoint, 
                          initial.growth=log(1.64), domain = list(geometry = fm_subdivide(precomputes$mesh,1), time = 1:6),
                          initial.carry.cap=log(10), initial.move.const = log(5e4),
                          initial.log.sigma = log(2), options = list(verbose = T, control.inla = list(strategy = "gaussian", int.strategy = "eb")))
#saveRDS(fit, "Eurasian Collared Dove/ecd_gamma0.5_moreint.RData")
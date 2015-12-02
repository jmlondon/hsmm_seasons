library(mhsmm)
load("predData.rda")
id = levels(factor(predData$deployid))
predData$ind_dry = ifelse(predData$percent_dry>100/3, 1, 0)
predData$y_disp_km = predData$y_disp/1000
predData$x_disp_km = predData$x_disp/1000

tmp1 = predData#[predData$deployid==id[1],]

### MHSMM functions
source("mhsmm_functions.R")

# Number of states
J = 4

# Initial vals
init0 = rep(1/J, J)
B0 = list(
  p = rep(0.1, J),
  lambda = rep(20, J), 
  mu = rep(list(c(0,0)), J),
  sigma = rep(list(50*diag(2)), J)
)
P0 = exp(-abs(row(diag(J)) - col(diag(J))))
diag(P0) = 0
P0 = sweep(P0, 1, rowSums(P0), "/")
S0 = list(lambda=rep(100, J), shift=rep(0,J), type="poisson")

start_val = hsmmspec(
  init = init0, 
  transition = P0, 
  parms.emission = B0, 
  sojourn = S0,
  dens.emission = dtelem.hsmm, 
  mstep = mstep.telem
)

# Fit model
data = list(x=with(tmp1, as.matrix(tmp1[,c("ind_dry","num_dives","x_disp_km","y_disp_km")])), N=table(tmp1$deployid))
fit = hsmmfit(data, start_val, mstep = mstep.telem)
summary(fit)


# # Number of states
# J = 3
# 
# # Initial vals
# init0 = rep(1/J, J)
# B0 = list(
#   p = rep(0.5, J),
#   lambda = rep(20, J), 
#   mu = rep(-10, J),
#   sigma = rep(10, J)
# )
# P0 = exp(-abs(row(diag(J)) - col(diag(J))))
# diag(P0) = 0
# P0 = sweep(P0, 1, rowSums(P0), "/")
# S0 = list(lambda=rep(100, J), shift=rep(100,J), type="poisson")
# 
# start_val = hsmmspec(
#   init = init0, 
#   transition = P0, 
#   parms.emission = B0, 
#   sojourn = S0,
#   dens.emission = dtelem.hsmm, 
#   mstep = mstep.telem
#   )
# 
# # Fit model
# data = list(x=with(tmp1, as.matrix(tmp1[,c("ind_dry","num_dives","y_disp_km")])), N=table(tmp1$deployid))
# fit = hsmmfit(data, start_val, mstep = mstep.telem)
# summary(fit)

# Predict states
predData$state = predict(fit, data)$s

library(ggplot2)
ggplot(data=predData) + geom_path(aes(y=state, x=datadatetime)) + facet_wrap(~deployid, ncol=1)

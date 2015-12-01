mstep.bern = function(x, wt){
  k = ncol(wt)
  p = numeric(k)
  for (i in 1:k) p[i] = weighted.mean(x, wt[, i], na.rm=TRUE)
  list(p = p)
}

mstep.pois2 = function (x, wt) 
{
  k = ncol(wt)
  lambda = numeric(k)
  for (i in 1:k) lambda[i] = weighted.mean(x, wt[, i], na.rm=TRUE)
  list(lambda = lambda)
}

mstep.telem = function(x, wt){
  k = ncol(wt)
  emission = list()
  emission["p"] = mstep.bern(x[,1], wt)
  emission["lambda"] = mstep.pois2(x[,2], wt)
  emission = append(emission, mstep.norm(x[,3], wt))
  emission$sigma = sqrt(emission$sigma)
  emission
}

dtelem.hsmm = function(x, j, model){
  d1 = dbinom(x[,1], 1, model$parms.emission$p[j], TRUE)
  d1[is.na(d1)] = 0
  d2 = dpois(x[,2], model$parms.emission$lambda[j], TRUE)
  d2[is.na(d2)] = 0
  d3 = dnorm(x[,3], model$parms.emission$mu[j], model$parms.emission$sigma[j], TRUE)
  d3[is.na(d3)] = 0
  exp(d1+d2+d3)
}

rtelem.hsmm = function(j, model){
  out = vector(length=3)
  out[1] = rbinom(1, 1, model$parms.emission$p[[j]])
  out[2] = rpois(1, model$parms.emission$lambda[[j]])
  out[3] = rnorm(1, model$parms.emission$mu[[j]], model$parms.emission$sigma[[j]])
  out
}



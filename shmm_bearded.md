# Estimating Seasonal Behavior States from Biologging Sensor Data
`r format(Sys.time(), '%d %B, %Y')`  




## Introduction

Since 2006, over 128 ribbon and spotted seals have been captured and outfitted with bio-logging devices that measure hourly wet/dry proportions. In addition, seven adult bearded seals were released with bio-logging devices. The critical life history events of these seals are all associated with increased haul-out behavior and changes in dive behavior. Additionally, migratory-like movements are also key indications of seasonal states. Hidden Markov models require a geometrically distributed sojourn time in a given state. Hidden semi-Markov models allow an arbitrary sojourn distribution --- the duration an animal spends in a state can depend on the time it has already spent in that state. O’Connell et al (2010) applied Hidden semi-Markov models to the estrus detection in dairy cows and developed the `mhsmm` library for R to support similar analyses.


## Methods and Analysis

We will start by loading the bearded seal data from the `kotzeb0912` package and the ribbon and spotted seal data from the `pepteldata` package. We will also rely on several other packages for the data munging and analysis. Data within these packages --- `kotzeb0912` and `pepteldata`) --- have been processed through the Wildlife Computers Data Portal and additionally munged for more efficient incorporation into analyses like this.



### Model Seal Movement

The first step will be to load in the location data and run these data through the `crawl` package to model seal movement and estimate locations every 6 hours. The `kotzeb0912` and `pepteldata` have slightly different organizational structures and column headings (will be fixed in the near future) so we need to do some work before merging them together.


```r
data("kotzeb0912_locs")
data("kotzeb0912_depths")
data("kotzeb0912_timelines")
data("peptel_locations")
data("peptel_drytimes")
data("peptel_divedepth")

kotzeb0912_locs <- kotzeb0912_locs %>% 
  dplyr::select(deployid,date_time=unique_posix,instr,
                quality,latitude,longitude,error_radius,
                error_semimajor_axis,
                error_semiminor_axis,
                error_ellipse_orientation) %>% 
  dplyr::filter(instr == "Mk10")  

kotzeb0912_depths <- kotzeb0912_depths %>% 
  dplyr::filter(deployid %in% kotzeb0912_locs$deployid) %>% 
  dplyr::rename(date_time=datadatetime)

kotzeb0912_timelines <- kotzeb0912_timelines %>% 
  dplyr::filter(deployid %in% kotzeb0912_locs$deployid) %>% 
  dplyr::rename(date_time=datadatetime)

peptel_locations <- peptel_locations %>% 
  dplyr::select(deployid,date_time=datadatetime,instr=tagtype,
                quality,latitude,longitude,error_radius,
                error_semimajor_axis = error_semi_major_axis,
                error_semiminor_axis = error_semi_minor_axis,
                error_ellipse_orientation) %>% 
  dplyr::mutate(date_time = lubridate::with_tz(date_time,"UTC")) %>% 
  dplyr::filter(instr %in% c("SPLASH","MK10A","MK10","ST16"),
                !deployid %in% kotzeb0912_locs$deployid)

peptel_ids <- peptel_locations %>% group_by(deployid) %>% 
  filter(n()>30) %>% select(deployid)

peptel_locations <- peptel_locations %>% 
  filter(deployid %in% peptel_ids$deployid)

library(xts)
date_unique <-peptel_locations %>% 
  group_by(deployid) %>%
  do(unique_date = xts::make.time.unique(.$date_time,eps=1)) %>%
  tidyr::unnest(unique_date) %>%
  mutate(unique_posix = as.POSIXct(.$unique_date,origin='1970-01-01 00:00:00',tz='UTC')) %>%
  dplyr::arrange(deployid,unique_posix) %>% 
  dplyr::select(unique_posix) %>% 
  ungroup()

peptel_locations <- peptel_locations %>% arrange(deployid,date_time) %>%
  bind_cols(date_unique) %>% select(-date_time) %>% 
  rename(date_time = unique_posix)

peptel_drytimes <- peptel_drytimes %>% 
  dplyr::filter(deployid %in% peptel_locations$deployid) %>% 
  dplyr::mutate(datadatetime = with_tz(datadatetime,"UTC")) %>% 
  dplyr::select(deployid,date_time=datadatetime,percent_dry=dry)

peptel_divedepth <- peptel_divedepth %>% 
  tidyr::gather(bin,num_dives, starts_with('bin')) %>%
  dplyr::mutate(date_time = with_tz(datadatetime, "UTC")) %>% 
  dplyr::select(deployid,
                instr=tagtype,
                date_time,
                -datadatetime,
                histtype,bin,num_dives) %>%
  dplyr::select(deployid,date_time,num_dives,bin) %>% 
  dplyr::filter(deployid %in% peptel_locations$deployid) %>% 
  dplyr::arrange(deployid,date_time,bin) 

locs <- bind_rows(peptel_locations,kotzeb0912_locs)
depths <- bind_rows(peptel_divedepth,kotzeb0912_depths)
drytimes <- bind_rows(peptel_drytimes,kotzeb0912_timelines)

locs <- locs %>% 
  mutate(species = ifelse(grepl("^EB.",deployid),'Bearded seal',NA),
         species = ifelse(grepl("^PL.",deployid),'Spotted seal',species),
         species = ifelse(grepl("^HF.",deployid),'Ribbon seal',species)) %>% 
  dplyr::filter(!is.na(species)) %>% arrange(deployid,date_time)

depths <- depths %>% 
  mutate(species = ifelse(grepl("^EB.",deployid),'Bearded seal',NA),
         species = ifelse(grepl("^PL.",deployid),'Spotted seal',species),
         species = ifelse(grepl("^HF.",deployid),'Ribbon seal',species))%>% 
  dplyr::filter(!is.na(species)) %>% arrange(deployid,date_time)

drytimes <- drytimes %>% 
  mutate(species = ifelse(grepl("^EB.",deployid),'Bearded seal',NA),
         species = ifelse(grepl("^PL.",deployid),'Spotted seal',species),
         species = ifelse(grepl("^HF.",deployid),'Ribbon seal',species))%>% 
  dplyr::filter(!is.na(species)) %>% arrange(deployid,date_time)
```

Prior to modeling the movement with `crawl`, we will use the `argosfilter` package to pass the Argos locations through a course speed filter (vmax=3.5 m/s).


```r
# speedfilter using the paralell package for multi-core speed
cfilter<-mclapply(split(locs,locs$deployid),function(x) sdafilter(
  lat=x$latitude, lon=x$longitude, dtime=x$date_time,
  lc=x$quality, ang=-1,vmax=3.5),mc.preschedule=F,mc.cores=3)
cfilter<-do.call("c",cfilter)
cfilter<-as.vector(cfilter)
locs$filtered <- cfilter

locs_data <- filter(locs,filtered=="not") %>% arrange(deployid,date_time) %>% 
  as.data.frame(.)
```

A few of the bearded seals were deployed with GPS capable tags and so we need to merge the GPS data with our Argos data.


```r
data("kotzeb0912_gps")

kotzeb0912_gps <- kotzeb0912_gps %>% 
  mutate(species = "Bearded seal") %>% 
  select(deployid,species,date_time,latitude,longitude) %>% 
  filter(!is.na(latitude)) %>% 
  mutate(error_semimajor_axis = 50,error_semiminor_axis = 50, error_ellipse_orientation = 1) 

locs_data <- locs_data %>% bind_rows(kotzeb0912_gps) %>% 
  filter(!grepl("^EB2004.",deployid)) %>% 
  filter(!grepl("^EB2005.",deployid)) %>% 
  filter(!grepl("^EB2006.",deployid)) %>% 
  filter(!grepl("^PL2005.",deployid)) %>% 
  filter(!grepl("^PL2006.",deployid)) %>% 
  arrange(deployid,date_time)

locs_data_ellipse<- filter(locs_data,!is.na(error_semimajor_axis)) %>% 
  arrange(deployid,date_time) %>% 
  as.data.frame(.)

locs_data_argos <- filter(locs_data,is.na(error_semimajor_axis),
                          !deployid %in% locs_data_ellipse$deployid) %>% 
  arrange(deployid,date_time) %>% 
  as.data.frame(.)
```

The specifications for the movement model are fairly straightforward. Note we are using the development version of `crawl` and incorporating the error ellipse structure. The GPS locations were all given a fixed error of 50m.


```r
data <- as.data.frame(locs_data_ellipse)
coordinates(data) = ~longitude+latitude
proj4string(data) = CRS("+proj=longlat +datum=WGS84")

data=spTransform(data, CRS("+init=epsg:3571"))

## Loop over PTT and fit CTCRW models
ids = unique(data@data$deployid)

library(doParallel)
n.cores <- detectCores()
registerDoParallel(cores=3)

model_fits_ellipse <-
  foreach(i = 1:length(ids)) %dopar% {
  id_data = subset(data,deployid==ids[i])
  diag_data = model.matrix(
      ~ error_semimajor_axis + error_semiminor_axis + error_ellipse_orientation,
      id_data@data
    )[,-1]
    
    id_data@data = cbind(id_data@data, 
                         argosDiag2Cov(
                           diag_data[,1], 
                           diag_data[,2], 
                           diag_data[,3]))
    
    init = list(a = c(coordinates(id_data)[1,1],0,
                      coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                           5000 ^ 2, 10 * 3600 ^ 2)))
    
  fit <- crwMLE(
      mov.model =  ~ 1,
      err.model = list(
        x =  ~ ln.sd.x - 1, 
        y =  ~ ln.sd.y - 1, 
        rho =  ~ error.corr
      ),
      data = id_data,
      Time.name = "date_time",
      initial.state = init,
      fixPar = c(1,1,NA,NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }

names(model_fits_ellipse) <- ids

data <- as.data.frame(locs_data_argos)
data$quality <- factor(data$quality,
                       levels=c("3", "2", "1","0","A","B"))

coordinates(data) = ~longitude+latitude
proj4string(data) = CRS("+proj=longlat +datum=WGS84")

data=spTransform(data, CRS("+init=epsg:3571"))

## Loop over PTT and fit CTCRW models
ids = unique(data@data$deployid)

library(doParallel)
n.cores <- detectCores()
registerDoParallel(cores=3)

model_fits_argos <-
  foreach(i = 1:length(ids)) %dopar% {
  id_data = subset(data,deployid==ids[i])
    
    init = list(a = c(coordinates(id_data)[1,1],0,
                      coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                           5000 ^ 2, 10 * 3600 ^ 2)))
    
    fixPar = c(log(250), log(500), log(1500), rep(NA,4), NA)
    
    constr=list(lower=c(rep(log(1500),2), rep(-Inf,3)),
            upper=rep(Inf,5))
    
    ln.prior = function(theta){-abs(theta[4]-3)/0.5}
    
  fit <- crwMLE(
      mov.model =  ~ 1,
      err.model = list(
        x =  ~ quality - 1
      ),
      data = id_data,
      Time.name = "date_time",
      initial.state = init,
      fixPar = fixPar,
      constr = constr,
      prior = ln.prior,
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }

names(model_fits_argos) <- ids


# print(model_fits_ellipse)
# print(model_fits_argos)
```

### Predict locations every 6 hours

This function predicts the regular-timed -- in this case, 6-hourly -- locations along the movement path using the posterior mean and variance of the track. 


```r
predData_ellipse <- foreach(i = 1:length(model_fits_ellipse)) %dopar% {
  if(!inherits(model_fits_ellipse[[i]],"try-error")) {
  model_fits_ellipse[[i]]$data$date_time <- lubridate::with_tz(
    model_fits_ellipse[[i]]$data$date_time,"UTC")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits_ellipse[[i]]$data$date_time),"day"),
    lubridate::floor_date(max(model_fits_ellipse[[i]]$data$date_time),"day"),
    "6 hours")
  tmp = crwPredict(model_fits_ellipse[[i]], predTime=predTimes)
  } else return(NA)
}

predData_ellipse <- Filter(function(x) inherits(x,"crwPredict"),predData_ellipse)

predData_ellipse <- dplyr::bind_rows(predData_ellipse) %>% 
  filter(locType=="p")

predData_argos <- foreach(i = 1:length(model_fits_argos)) %dopar% {
  if(!inherits(model_fits_argos[[i]],"try-error")) {
  model_fits_argos[[i]]$data$date_time <- lubridate::with_tz(
    model_fits_argos[[i]]$data$date_time,"UTC")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits_argos[[i]]$data$date_time),"day"),
    lubridate::floor_date(max(model_fits_argos[[i]]$data$date_time),"day"),
    "6 hours")
  tmp = crwPredict(model_fits_argos[[i]], predTime=predTimes)
  } else return(NA)
}

predData_argos <- Filter(function(x) inherits(x,"crwPredict"),predData_argos)

predData_argos <- dplyr::bind_rows(predData_argos) %>% 
  filter(locType=="p")

predData <- dplyr::bind_rows(predData_ellipse,predData_argos)

predData$predTimes <- intToPOSIX(predData$TimeNum)
```


### Merge predicted locations with haul-out data

The haul-out behavior timelines are provided as hour percent-dry values. We want to group these values into 6-hour blocks that are centered on our 6-hourly predictions. To accomplish this, we will setup a grouping column that assigns a unique integer to each 6 hour period.







### Merge predicted locations with dive behavior

Dive behavior data are transmitted as dive histograms that represent the distribution of dives across predetermined depth bins for a given 6-hour period. To simplify things, we will sum the number of dives across all bins less than 10m. This should represent the most likely number of foraging dives over a 6 hour period. Future analysis may consider using time-at-depth instead of number of dives.


```r
diveData <- depths %>% 
  filter(bin != "bin1") %>% 
  mutate(date_time = date_time + hours(3)) %>% 
  group_by(deployid,date_time) %>% 
  summarise(num_dives=sum(num_dives,na.rm=TRUE)) %>% 
  ungroup()

predData <- dplyr::left_join(predData,diveData) %>% 
  arrange(deployid,date_time) 
```

### Calculate x and y displacement for each step

We'll use the dplyr::lag() function to calculate the difference in x and y when compared with the previous x and y values. This will give us a measure of northing and easting displacement at each time step. These values will form the basis for our xy-displacement multi-variate normal parameter in the model.


```r
predData <- predData %>% 
  group_by(deployid) %>% 
  mutate(x_disp = order_by(date_time, mu.x - lag(mu.x)), 
         y_disp = order_by(date_time, mu.y - lag(mu.y))
  ) %>% 
  filter(!is.na(x_disp) | !is.na(y_disp)) %>% 
  ungroup()
```

We can use a custom function to calculate the compass bearing between those two points. Compass bearing provides a more user friendly description of movement compared to x/y displacement.


```r
anglefun <- function(xx,yy,bearing=TRUE,as.deg=FALSE){
  ## calculates the compass bearing of the line between two points
  ## xx and yy are the differences in x and y coordinates between two points
  ## Options:
  ## bearing = FALSE returns +/- pi instead of 0:2*pi
  ## as.deg = TRUE returns degrees instead of radians
  c = 1
  if (as.deg){
    c = 180/pi
  }
  
  b<-sign(xx)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(yy<0)*pi+atan(xx/yy)
  if(bearing){
    #return a compass bearing 0 to 2pi
    #if bearing==FALSE then a heading (+/- pi) is returned
    tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  }
  return(tempangle*c)
}

predData <- predData %>% rowwise() %>% 
  mutate(bearing = anglefun(x_disp,y_disp,as.deg=TRUE)) %>% 
  ungroup()
```

### Add deploy_day for aligning deployments

Since these deployments occurred over multiple years, we will want a convenient way to align the deployments. The simplest way to do this is to use day-of-year integers and add 365 to days in January-May. At this point, we will also do some additional filtering to make sure we don't have any corrupt or other outlier parameters from the tag data.


```r
id_apr1 <- predData %>% 
  group_by(deployid) %>% summarise(start_yr = year(min(date_time))) %>% 
  mutate(apr1 = lubridate::ymd(paste(start_yr,"04-01",sep="-"))) %>% 
  select(-start_yr) %>% ungroup()

predData <- predData %>% left_join(id_apr1) %>% 
  mutate(days_since_apr1 = as.numeric(difftime(date_time,apr1,units="days")))


predData <- predData %>% 
  mutate(ho_binary=ifelse(percent_dry>=33.3,1,-1),
         y_disp = ifelse(abs(y_disp)>100000,NA,y_disp),
         x_disp = ifelse(abs(x_disp)>100000,NA,x_disp),
         num_dives = ifelse(num_dives > 150,NA,num_dives)) 

predData <- predData %>% 
  filter(days_since_apr1 <=400)

select_ids <- predData %>% group_by(deployid) %>% 
  filter(n()>720) %>% select(deployid)

remove_ids <- c("PL2009_2005_06L0126","HF2005_5898_05L0001","HF2009_1019_09L0040",
                "HF2005_5880_05L0003","HF2005_5882_05L0008","HF2005_5884_05L0009",
                "HF2005_5890_05L0011","HF2005_5894_05L0013","HF2005_5896_05L0015")

predData <- predData %>% 
  filter(deployid %in% select_ids$deployid,
         !deployid %in% remove_ids)

load('pep_age.rda')

predData <- predData %>% left_join(pep_age)

save(predData,file='predData.rda')
```

We now present a series of four plots showing the raw, 'observed' parameter values that will go into the model: haul-out status, number of dives, y-displacement, and x-displacement.

![](shmm_bearded_files/figure-html/plot-1-1.png) ![](shmm_bearded_files/figure-html/plot-1-2.png) ![](shmm_bearded_files/figure-html/plot-1-3.png) 

![](shmm_bearded_files/figure-html/plot-2-1.png) ![](shmm_bearded_files/figure-html/plot-2-2.png) ![](shmm_bearded_files/figure-html/plot-2-3.png) 

![](shmm_bearded_files/figure-html/plot-3-1.png) 

![](shmm_bearded_files/figure-html/plot-4-1.png) 

### multivariate semi-hidden markov

We are going to use the mshmm package to run our multi-variate semi-hidden markov model. For the model, we will use 6-hour time steps and the following parameters:

#. haul-out status (Bernouli; cuttoff at 33.3% dry)
#. number of dives below 10m (Poisson)
#. x-y displacement (multi-variate normal)

### bearded seal fit


```r
library(mhsmm)

dd <- predData_eb

id = levels(factor(dd$deployid))
dd$ind_dry = ifelse(dd$percent_dry>100/3, 1, 0)
dd$y_disp_km = dd$y_disp/1000
dd$x_disp_km = dd$x_disp/1000

tmp1 = dd#[predData$deployid==id[1],]

### MHSMM functions
source("mhsmm_functions.R")
```

After loading the package and doing a little bit of data cleaning, we need to setup the initial values and parameters. Initial investigation and model comparison with AIC suggests that 4 states provides the best model.


```r
# Number of states
J = 4

# Initial vals
init0 = rep(1/J, J)
B0 = list(
  p = c(0.05,0.05,0.07,0.1),
  #p = c(0.1, 0.1, 0.1),
  lambda = rep(30, J), 
  mu = rep(list(c(0,0)), J),
  sigma = rep(list(20*diag(2)), J)
)
P0 = exp(-abs(row(diag(J)) - col(diag(J))))
diag(P0) = 0
P0 = sweep(P0, 1, rowSums(P0), "/")
S0 = list(lambda=rep(150, J), shift=rep(200,J), type="poisson")

start_val = hsmmspec(
  init = init0, 
  transition = P0, 
  parms.emission = B0, 
  sojourn = S0,
  dens.emission = dtelem.hsmm.mv, 
  mstep = mstep.telem.mv
)
```

Now, we can fit the model and predict


```r
# Fit model
data = list(x=with(tmp1, as.matrix(tmp1[,c("ind_dry","num_dives","x_disp_km","y_disp_km")])), N=table(tmp1$deployid))
#M <- 600
fit = hsmmfit(data, start_val, mstep = mstep.telem.mv)

dd$state = predict(fit, data)$s
dd$state = as.factor(dd$state)

fit_bearded <- fit
dd_bearded <- dd
```


### ribbon seal fit


```r
library(mhsmm)

dd <- predData_hf %>% filter(age != "YOY")

id = levels(factor(dd$deployid))
dd$ind_dry = ifelse(dd$percent_dry>100/3, 1, 0)
dd$y_disp_km = dd$y_disp/1000
dd$x_disp_km = dd$x_disp/1000

tmp1 = dd#[predData$deployid==id[1],]

### MHSMM functions
source("mhsmm_functions.R")
```

After loading the package and doing a little bit of data cleaning, we need to setup the initial values and parameters. Initial investigation and model comparison with AIC suggests that 4 states provides the best model.


```r
# Number of states
J = 3

# Initial vals
init0 = rep(1/J, J)
B0 = list(
  #p = c(0.05,0.05,0.07,0.1),
  p = c(0.1, 0.1, 0.1),
  lambda = rep(30, J), 
  mu = rep(list(c(0,0)), J),
  sigma = rep(list(20*diag(2)), J)
)
P0 = exp(-abs(row(diag(J)) - col(diag(J))))
diag(P0) = 0
P0 = sweep(P0, 1, rowSums(P0), "/")
S0 = list(lambda=rep(150, J), shift=rep(240,J), type="poisson")

start_val = hsmmspec(
  init = init0, 
  transition = P0, 
  parms.emission = B0, 
  sojourn = S0,
  dens.emission = dtelem.hsmm.mv, 
  mstep = mstep.telem.mv
)
```

Now, we can fit the model and predict


```r
# Fit model
data = list(x=with(tmp1, as.matrix(tmp1[,c("ind_dry","num_dives","x_disp_km","y_disp_km")])), N=table(tmp1$deployid))
#M <- 600
fit = hsmmfit(data, start_val, mstep = mstep.telem.mv)

dd$state = predict(fit, data)$s
dd$state = as.factor(dd$state)

fit_ribbon <- fit
dd_ribbon <- dd
```

### YOY spotted seal fit


```r
library(mhsmm)

dd <- predData_pl %>% filter(age == "YOY")

id = levels(factor(dd$deployid))
dd$ind_dry = ifelse(dd$percent_dry>100/3, 1, 0)
dd$y_disp_km = dd$y_disp/1000
dd$x_disp_km = dd$x_disp/1000

tmp1 = dd#[predData$deployid==id[1],]

### MHSMM functions
source("mhsmm_functions.R")
```

After loading the package and doing a little bit of data cleaning, we need to setup the initial values and parameters. Initial investigation and model comparison with AIC suggests that 4 states provides the best model.


```r
# Number of states
J = 3

# Initial vals
init0 = rep(1/J, J)
B0 = list(
  #p = c(0.05,0.05,0.07,0.1),
  p = c(0.1, 0.1, 0.1),
  lambda = rep(50, J), 
  mu = rep(list(c(0,0)), J),
  sigma = rep(list(20*diag(2)), J)
)
P0 = exp(-abs(row(diag(J)) - col(diag(J))))
diag(P0) = 0
P0 = sweep(P0, 1, rowSums(P0), "/")
S0 = list(lambda=rep(100, J), shift=rep(240,J), type="poisson")

start_val = hsmmspec(
  init = init0, 
  transition = P0, 
  parms.emission = B0, 
  sojourn = S0,
  dens.emission = dtelem.hsmm.mv, 
  mstep = mstep.telem.mv
)
```

Now, we can fit the model and predict


```r
# Fit model
data = list(x=with(tmp1, as.matrix(tmp1[,c("ind_dry","num_dives","x_disp_km","y_disp_km")])), N=table(tmp1$deployid))
#M <- 600
fit = hsmmfit(data, start_val, mstep = mstep.telem.mv)

dd$state = predict(fit, data)$s
dd$state = as.factor(dd$state)

fit_spotted_yoy <- fit
dd_spotted_yoy <- dd
```

## Results

### Adult Bearded Seals

The results of the model fit are provided


```r
summary(fit_bearded)
```

```
## 
## Starting distribution = 
## [1] 8.6e-01 2.5e-16 1.4e-01 2.1e-17
## 
## Transition matrix = 
##      [,1] [,2] [,3] [,4]
## [1,]    0 0.33 0.67 0.00
## [2,]    0 0.00 0.50 0.50
## [3,]    0 0.29 0.00 0.71
## [4,]    0 1.00 0.00 0.00
## 
## Sojourn distribution parameters = 
## $lambda
## [1] 589.77507 286.63513 155.12801  13.15015
## 
## $shift
## [1]   1   1   1 374
## 
## $type
## [1] "poisson"
## 
## 
## Emission distribution parameters = 
## $p
## [1] 0.04148159 0.06957714 0.10952382 0.09925803
## 
## $lambda
## [1] 42.55516 37.51165 32.50710 28.38349
## 
## $mu
## $mu[[1]]
##   x_disp_km   y_disp_km 
## -0.41117401  0.04501103 
## 
## $mu[[2]]
##  x_disp_km  y_disp_km 
##  0.6945921 -0.5436585 
## 
## $mu[[3]]
##  x_disp_km  y_disp_km 
## -0.3194037 -2.3123997 
## 
## $mu[[4]]
##    x_disp_km    y_disp_km 
## -0.001394739 -0.123325748 
## 
## 
## $sigma
## $sigma[[1]]
##           x_disp_km y_disp_km
## x_disp_km 43.578295  1.825728
## y_disp_km  1.825728 47.899414
## 
## $sigma[[2]]
##           x_disp_km y_disp_km
## x_disp_km 19.883499  1.084859
## y_disp_km  1.084859 22.391858
## 
## $sigma[[3]]
##           x_disp_km y_disp_km
## x_disp_km 58.627959  6.269766
## y_disp_km  6.269766 76.753714
## 
## $sigma[[4]]
##            x_disp_km  y_disp_km
## x_disp_km 12.3047922 -0.2391395
## y_disp_km -0.2391395 16.1575045
```

And, we can plot the state assignments for each bearded seal

![](shmm_bearded_files/figure-html/plot-5-1.png) 

We can also combine the state assignments across our tagged seals by taking the majority state for each time step. The color transparency is set to the proportion of animals represented by that majority state. Note that, some states may not be present in the combined graph.

![](shmm_bearded_files/figure-html/plot-6-1.png) 

Now that we have some seasonal states assigned, we can examine the distribution of various behaviors across those different states.

![](shmm_bearded_files/figure-html/plot-7-1.png) 

![](shmm_bearded_files/figure-html/plot-8-1.png) 

![](shmm_bearded_files/figure-html/plot-9-1.png) 


### Adult and Sub-Adult Ribbon Seals

The results of the model fit are provided


```r
summary(fit_ribbon)
```

```
## 
## Starting distribution = 
## [1] 0.38 0.24 0.38
## 
## Transition matrix = 
##      [,1] [,2] [,3]
## [1,] 0.00 0.29 0.71
## [2,] 0.50 0.00 0.50
## [3,] 0.71 0.29 0.00
## 
## Sojourn distribution parameters = 
## $lambda
## [1] 76.01070 45.35333 76.01070
## 
## $shift
## [1] 249 615 249
## 
## $type
## [1] "poisson"
## 
## 
## Emission distribution parameters = 
## $p
## [1] 0.1119314 0.1289413 0.1119314
## 
## $lambda
## [1] 40.55069 14.90648 40.55069
## 
## $mu
## $mu[[1]]
##   x_disp_km   y_disp_km 
## 0.066615608 0.006877309 
## 
## $mu[[2]]
##   x_disp_km   y_disp_km 
## -0.40881341 -0.08987032 
## 
## $mu[[3]]
##   x_disp_km   y_disp_km 
## 0.066615608 0.006877309 
## 
## 
## $sigma
## $sigma[[1]]
##            x_disp_km  y_disp_km
## x_disp_km 55.8775464 -0.9962993
## y_disp_km -0.9962993 73.3356691
## 
## $sigma[[2]]
##           x_disp_km y_disp_km
## x_disp_km  21.12194  -2.38025
## y_disp_km  -2.38025  31.18892
## 
## $sigma[[3]]
##            x_disp_km  y_disp_km
## x_disp_km 55.8775464 -0.9962993
## y_disp_km -0.9962993 73.3356691
```

And, we can plot the state assignments for each ribbon seal

![](shmm_bearded_files/figure-html/plot-10-1.png) 

We can also combine the state assignments across our tagged seals by taking the majority state for each time step. The color transparency is set to the proportion of animals represented by that majority state. Note that, some states may not be present in the combined graph.

![](shmm_bearded_files/figure-html/plot-11-1.png) 

Now that we have some seasonal states assigned, we can examine the distribution of various behaviors across those different states.

![](shmm_bearded_files/figure-html/plot-12-1.png) 

![](shmm_bearded_files/figure-html/plot-13-1.png) 

![](shmm_bearded_files/figure-html/plot-14-1.png) 


### Young-of-Year Spotted Seals

The results of the model fit are provided


```r
summary(fit_spotted_yoy)
```

```
## 
## Starting distribution = 
## [1] 0.38 0.25 0.38
## 
## Transition matrix = 
##      [,1] [,2] [,3]
## [1,]  0.0  0.2  0.8
## [2,]  0.5  0.0  0.5
## [3,]  0.8  0.2  0.0
## 
## Sojourn distribution parameters = 
## $lambda
## [1]   4.884429 327.639782   4.884429
## 
## $shift
## [1] 1 1 1
## 
## $type
## [1] "poisson"
## 
## 
## Emission distribution parameters = 
## $p
## [1] 0.38223824 0.03015299 0.38223824
## 
## $lambda
## [1] 32.88793 88.66581 32.88793
## 
## $mu
## $mu[[1]]
##   x_disp_km   y_disp_km 
##  0.02508787 -0.13498423 
## 
## $mu[[2]]
##  x_disp_km  y_disp_km 
## 0.01759796 0.04012512 
## 
## $mu[[3]]
##   x_disp_km   y_disp_km 
##  0.02508787 -0.13498423 
## 
## 
## $sigma
## $sigma[[1]]
##           x_disp_km y_disp_km
## x_disp_km 45.240413  1.244627
## y_disp_km  1.244627 35.962200
## 
## $sigma[[2]]
##           x_disp_km y_disp_km
## x_disp_km  83.83158  -0.47421
## y_disp_km  -0.47421  66.70490
## 
## $sigma[[3]]
##           x_disp_km y_disp_km
## x_disp_km 45.240413  1.244627
## y_disp_km  1.244627 35.962200
```

And, we can plot the state assignments for each young spotted seal

![](shmm_bearded_files/figure-html/plot-15-1.png) 

We can also combine the state assignments across our tagged seals by taking the majority state for each time step. The color transparency is set to the proportion of animals represented by that majority state. Note that, some states may not be present in the combined graph.

![](shmm_bearded_files/figure-html/plot-16-1.png) 

Now that we have some seasonal states assigned, we can examine the distribution of various behaviors across those different states.

![](shmm_bearded_files/figure-html/plot-17-1.png) 

![](shmm_bearded_files/figure-html/plot-18-1.png) 

![](shmm_bearded_files/figure-html/plot-19-1.png) 



## Discussion



## References

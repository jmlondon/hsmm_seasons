# Estimating Seasonal Behavior States from Biologging Sensor Data
`r format(Sys.time(), '%d %B, %Y')`  




## Introduction

Since 2006, over 128 ribbon and spotted seals have been captured and outfitted with bio-logging devices that measure hourly wet/dry proportions. In addition, seven adult bearded seals were released with bio-logging devices. The critical life history events are all associated with increased haul-out behavior. Hidden Markov models require a geometrically distributed sojourn time in a given state. Hidden semi-Markov models allow an arbitrary sojourn distribution. In other words, the Hidden semi-Markov approach allows for estimation of states that persist longer than the observation interval. O’Connell et al (2010) applied Hidden semi-Markov models to the estrus detection in dairy cows and developed the mhsmm library for R to support similar analyses.

We expect ribbon seals to have at least 3 distinct annual states: an open water or pelagic period characterized by limited to no haul-out behavior, pupping/breeding/molting period characterized by increased haul-out frequency and duration, a transition period leading up to the pupping/breeding/molting period with increasing haul-out behavior as sea ice forms and stabilizes within the Bering Sea. These periods are generally known to correspond to July-November, April-June and December-March, but the specifics are unknown. Spotted seals should mirror some of these patterns although the open water period will be likely replaced with a near-shore period. BeardedWhile, initially, the focus will be on a single observation channel (proportion dry per hour), a multivariate approach would allow inclusion of dive behavior. Additionally, we would expect variability between age and sex classes and for key covariates (e.g. sea-ice concentration or distance to sea-ice edge) to influence state assignment.


## Methods and Analysis

We will start by loading the bearded seal data from the `kotzeb0912` package along with other packages we'll need for the data munging and analysis.


```r
library(kotzeb0912)
library(sp)
library(rgdal)
library(trip)
library(crawl)
library(argosfilter)
library(parallel)
library(dplyr)
library(lubridate)
```

### Model Seal Movement

The first step will be to load in the location data and run these data through the `crawl` package to model seal movement and estimate locations every 6 hours.

Prior to modeling the movement with `crawl`, we will use the `argosfilter` package to pass the Argos locations through a course speed filter (vmax=3.5 m/s).


```r
data(kotzeb0912_locs)

locs <- dplyr::filter(kotzeb0912_locs,instr=="Mk10") %>% 
  dplyr::arrange(deployid,unique_posix)
# speedfilter using the paralell package for multi-core speed
cfilter<-mclapply(split(locs,locs$deployid),function(x) sdafilter(
  lat=x$latitude, lon=x$longitude, dtime=x$unique_posix,
  lc=x$quality, ang=-1,vmax=3.5),mc.preschedule=F,mc.cores=3)
cfilter<-do.call("c",cfilter)
cfilter<-as.vector(cfilter)
locs$filtered <- cfilter

data <- filter(locs,filtered=="not", !is.na(error_semimajor_axis)) %>% arrange(deployid,unique_posix) %>% as.data.frame(.)
```

A few of the bearded seals were deployed with GPS capable tags and so we need to merge the GPS data with our Argos data.


```r
data(kotzeb0912_gps)

kotzeb0912_gps <- kotzeb0912_gps %>% 
  select(deployid,date_time,latitude,longitude) %>% 
  filter(!is.na(latitude)) %>% 
  rename(unique_posix=date_time) %>% 
  mutate(error_semimajor_axis = 50,error_semiminor_axis = 50, error_ellipse_orientation = 1) 

data <- data %>% bind_rows(kotzeb0912_gps) %>% 
  arrange(deployid,unique_posix)
```



```r
data <- as.data.frame(data)
coordinates(data) = ~longitude+latitude
proj4string(data) = CRS("+proj=longlat +datum=WGS84")

data=spTransform(data, CRS("+init=epsg:3571"))

## Loop over PTT and fit CTCRW models
ids = unique(data@data$deployid)

library(doParallel)
n.cores <- detectCores()
registerDoParallel(cores=3)

model_fits <-
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
      Time.name = "unique_posix",
      initial.state = init,
      fixPar = c(1,1,NA,NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }

names(model_fits) <- ids

print(model_fits)
```

### Predict locations every 6 hours

This function predicts the regular-timed -- in this case, 6-hourly -- locations along the movement path using the posterior mean and variance of the track. 


```r
predData <- foreach(i = 1:length(model_fits)) %dopar% {

  model_fits[[i]]$data$unique_posix <- lubridate::with_tz(
    model_fits[[i]]$data$unique_posix,"GMT")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix),"day"),
    lubridate::floor_date(max(model_fits[[i]]$data$unique_posix),"day"),
    "6 hours")
  tmp = crwPredict(model_fits[[i]], predTime=predTimes)
}


predData <- dplyr::bind_rows(predData) %>% 
  filter(locType=="p")

predData$predTimes <- intToPOSIX(predData$TimeNum)
```


### Merge predicted locations with haul-out data

The haul-out behavior timelines are provided as hour percent-dry values. We want to group these values into 6-hour blocks that are centered on our 6-hourly predictions. To accomplish this, we will setup a grouping column that assigns a unique integer to each 6 hour period.


```r
data(kotzeb0912_timelines)

kotzeb0912_timelines <- kotzeb0912_timelines %>% 
  filter(deployid %in% unique(predData$deployid))
```


```r
ids <- unique(kotzeb0912_timelines$deployid)
timelineData <- foreach(i = 1:length(ids)) %dopar% {
  t_dat <- kotzeb0912_timelines %>% 
    filter(deployid == ids[i]) %>% 
    arrange(deployid,datadatetime)
  t_seq <- data.frame(datadatetime = seq(
    lubridate::ceiling_date(min(t_dat$datadatetime),"day"),
    lubridate::floor_date(max(t_dat$datadatetime),"day"),
    "1 hour"),deployid=ids[i])
  t_dat <- merge(t_dat,t_seq,all=TRUE)
  rep_len <- ceiling((nrow(t_dat)-3)/6)
  t_group <- c(rep(1,3),rep(2:rep_len,each=6))
  t_group <- data.frame(group=t_group[1:nrow(t_dat)])
  t_dat <-cbind(t_dat,t_group) %>% filter(!is.na(group))
}

timelineData <- dplyr::bind_rows(timelineData) %>% 
  group_by(deployid,group) %>% 
  summarise(percent_dry=mean(percent_dry,na.rm=TRUE),
            datadatetime=ceiling_date(median(datadatetime),"hour")) %>% 
  select(deployid,datadatetime,percent_dry)
```


```r
predData <- dplyr::inner_join(predData,
                              timelineData,
                              by = c("predTimes" = "datadatetime", 
                                     "deployid")) %>%
  tbl_df() %>%
  select(deployid,predTimes,mu.x,mu.y,percent_dry) %>% 
  arrange(deployid, predTimes) %>% 
  rename(datadatetime=predTimes)
```

### Merge predicted locations with dive behavior

Dive behavior data are transmitted as dive histograms that represent the distribution of dives across predetermined depth bins for a given 6-hour period. To simplify things, we will sum the number of dives across all bins less than 10m. This should represent the most likely number of foraging dives over a 6 hour period.


```r
data(kotzeb0912_depths)

diveData <- kotzeb0912_depths %>% 
  filter(limits != "10.000000") %>% 
  filter(num_dives < 80) %>% 
  mutate(datadatetime = datadatetime + lubridate::hours(3)) %>% 
  group_by(deployid,datadatetime) %>% 
  summarise(num_dives=sum(num_dives,na.rm=TRUE))

predData <- dplyr::left_join(predData,
                             diveData) %>% 
  tbl_df() %>% 
  arrange(deployid,datadatetime)
```

### Calculate x and y displacement for each step

We'll use the dplyr::lag() function to calculate the difference in x and y when compared with the previous x and y values. This will give us a measure of northing and easting movement.


```r
predData <- predData %>% 
  group_by(deployid) %>% 
  mutate(x_disp = order_by(datadatetime, mu.x - lag(mu.x)), 
         y_disp = order_by(datadatetime, mu.y - lag(mu.y))
  ) %>% 
  filter(!is.na(x_disp) | !is.na(y_disp))
```

and the we can use a custom function to calculate the compass bearing between those two points


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

predData <- predData %>% 
  mutate(bearing = anglefun(x_disp,y_disp,as.deg=TRUE))
```

### Add Day-of-Year and Deployment-Day for aligning deployments

Since these deployments occurred over multiple years, we will want a convenient way to align the deployments. The simplest way to do this is to use day-of-year integers and pick start/end values that provide a sensible amount of overlap between the 7 deployments.

First, we'll look at the start dates for each deployment


```r
predData %>% group_by(deployid) %>% 
  filter(row_number(datadatetime) == 1)
```

```
## Source: local data frame [7 x 9]
## Groups: deployid [7]
## 
##              deployid        datadatetime     mu.x     mu.y percent_dry
##                 (chr)              (time)    (dbl)    (dbl)       (dbl)
## 1 EB2009_3000_06A1346 2009-06-24 12:00:00 775682.6 -2490010    3.000000
## 2 EB2009_3001_06A1332 2009-06-26 12:00:00 762830.2 -2478453    3.000000
## 3 EB2009_3002_06A1357 2009-06-27 12:00:00 766298.2 -2478069    3.000000
## 4 EB2011_3000_10A0219 2011-06-17 12:00:00 754848.2 -2472746    4.166667
## 5 EB2011_3001_10A0552 2011-06-18 12:00:00 771161.0 -2485771   67.666667
## 6 EB2011_3002_10A0200 2011-06-19 12:00:00 758008.8 -2476040    3.000000
## 7 EB2012_3003_09A0888 2012-07-05 12:00:00 725878.4 -2486480    3.000000
## Variables not shown: num_dives (dbl), x_disp (dbl), y_disp (dbl), bearing (dbl)
```

And, now, we'll look at the end dates for each deployment


```r
predData %>% group_by(deployid) %>% 
  filter(row_number(datadatetime) == max(row_number(datadatetime)))
```

```
## Source: local data frame [7 x 9]
## Groups: deployid [7]
## 
##              deployid        datadatetime     mu.x     mu.y percent_dry
##                 (chr)              (time)    (dbl)    (dbl)       (dbl)
## 1 EB2009_3000_06A1346 2010-04-19 00:00:00 883464.7 -2699166       100.0
## 2 EB2009_3001_06A1332 2010-03-08 18:00:00 888421.5 -2691736         1.0
## 3 EB2009_3002_06A1357 2010-02-25 18:00:00 399555.0 -3353737         0.5
## 4 EB2011_3000_10A0219 2012-03-31 18:00:00 710867.6 -3048724         3.0
## 5 EB2011_3001_10A0552 2012-01-31 18:00:00 417752.5 -2773844         3.0
## 6 EB2011_3002_10A0200 2012-03-17 18:00:00 448075.6 -3192078         2.0
## 7 EB2012_3003_09A0888 2013-03-09 18:00:00 688416.9 -2794321         3.0
## Variables not shown: num_dives (dbl), x_disp (dbl), y_disp (dbl), bearing (dbl)
```

From this information is seems sensible to have all of the deployments start on 05 July (186th day of the year). 


```r
predData <- predData %>% 
  mutate(deploy_day = ifelse(lubridate::yday(datadatetime)<150, 
                      lubridate::yday(datadatetime)+365 +
                        lubridate::hour(datadatetime)/24,
                      lubridate::yday(datadatetime) +
                        lubridate::hour(datadatetime)/24)) %>% 
  mutate(corrupt = ifelse(deployid == "EB2009_3001_06A1332" & 
                            datadatetime > lubridate::ymd_hms("2010-01-28 03:00:00"),
                          TRUE,FALSE
                          )) %>% 
  filter(abs(y_disp) < 100000, abs(x_disp) < 100000) %>% 
  filter(!corrupt)

save(predData,file='predData.rda')
```


```r
library(ggplot2)
library(uswebr)

p1 <- ggplot(predData, aes(x=deploy_day, y=percent_dry)) + geom_point(alpha=0.2) + facet_grid(deployid ~ .) + 
  usweb_theme()
p1
```

![](shmm_bearded_files/figure-html/plot-1-1.png) 


```r
p2 <- ggplot(predData, aes(x=deploy_day, y=num_dives)) + 
  geom_point(alpha=0.2) + facet_grid(deployid ~ .) + 
  usweb_theme()
p2
```

![](shmm_bearded_files/figure-html/plot-2-1.png) 


```r
predData <- predData %>% mutate(direction = ifelse(
  y_disp >= 0, "north","south"))

ggplot(predData,aes(x=deploy_day, y=y_disp/1000)) +
  geom_bar(aes(color=factor(direction)),alpha=0.2,stat="identity") +
  facet_grid(deployid ~ .) + usweb_theme()
```

![](shmm_bearded_files/figure-html/plot-3-1.png) 

## Results



## Discussion



## References

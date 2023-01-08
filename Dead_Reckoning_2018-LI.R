#### R code to calculate the Dead Reckoning of narwhals tracks in 2018 (and 2017) ####
## Used as the dataset to do further analysis in at least 2 papers mentioned in README.md
## It is used to estimate the effects of seismic activities on the natural behavior of narwhals in East Greenland
##
## We compute parameters ODBA/VeDBA, distance to shore and ship, and Area
## (Stroke rate (per minute) would be done by Outi Tervo).
## 
## Don't use MCMC (Stan) in this version, 
## since it's very slow and doesn't really improve the precision.
## 
## Use LINEAR INTERPOLATION instead of flownoise

rm(list = ls(all.names = TRUE)); closeAllConnections()
Sys.setenv(TZ='UTC'); options(digits = 12)
options(scipen=999)
options(warn=2)

library(data.table); library(lubridate)

#
#### 1a) Convert from Lat/Lon data (WGS84) into UTM function, and vice versa ####
# Zone 24N seems not be the good choice for Scorebysund
# http://spatialreference.org/ref/epsg/32624/
# but zone 26N or even zone 27N
# That's why Rene et al. in Lauge Koch ship use 27N1
# http://spatialreference.org/ref/epsg/32626/
# But note that whole Greenland uses 24N, because it's near the middle
# 26N is the better choice for small region like Scorebysund,
# to avoid distortion
# 
# Cureently have to use 24N to fit coast line, 
# but need to be change in the future

## Convert from Lat/Lon data (WGS84) into UTM function ##
# Input:
# - x: longitude time series
# - y: latitude time series
# - zoneID: UTM zone ID
# Output: 
# - a dataframe of UTM coordinates
LonLat_to_UTM = function(x,y,zoneID=24) {
  
  library(sp); library(data.table)
  dt0 = data.table(Lon=x,Lat=y)
  coordinates(dt0) = ~Lon+Lat
  # Setting default projection [https://gis.stackexchange.com/a/3342]
  proj4string(dt0) = CRS('+init=epsg:4326')
  # Projection argument
  projection_arg = paste0('+proj=utm +zone=',zoneID,
                          ' +datum=WGS84 +units=m +no_defs')
  data.table(data.frame(spTransform(dt0, CRS(projection_arg))))
}

## Convert UTM into Lat/Lon data (WGS84) function ##
# Input:
# - x: longitude time series
# - y: latitude time series
# - zoneID: UTM zone ID
# Output: 
# - a dataframe of longitude/latitude
UTM_to_LonLat = function(x,y,zoneID=24) {
  
  library(sp); library(data.table)
  dt0 = data.table(Lon=x,Lat=y)
  coordinates(dt0) = ~Lon+Lat
  # Setting default projection [https://gis.stackexchange.com/a/3342]
  # Projection argument
  projection_arg = paste0('+proj=utm +zone=',zoneID,' +datum=WGS84 +units=m +no_defs')
  proj4string(dt0) = CRS(projection_arg)
  data.table(data.frame(spTransform(dt0, CRS('+init=epsg:4326'))))
}

#### 1b) Read MT files (auxiliary) & join them together ####

library(mt)

# Extract processed accelerometer/magnetometer data pre-filtered by MATLAB
# Input:
# - inp: path of input file
# - s: indentified string to select the data inside the whole data, 
#      among accelerometer data ('HDI','HDJ','HDK'),
#      or magnetometer data ('HDX','HDY','HDZ')
# Output:
# - the dataframe of 
MAT_to_data_table = function (inp,s) {

  library(stringr)
  sdt = paste0('^',s,'.*data$')
  x = inp[str_detect( names(inp),sdt )]
  # TEST if HDIs fit together
  sdt = paste0('^',s,'.*start$')
  t = ymd_hms(inp[str_detect( names(inp),sdt )])
  if ( all(diff(as.numeric(t)) == 
           unlist(lapply(lapply(x, as.vector),length))[1:(length(t)-1)] ) ) {
    y = as.numeric(do.call(c,x))
    dt = data.table( DateTime = seq(t[1],t[1]+length(y)-1,by=1), Val=y )
  } else {
    stop( 'Something wrong with decimate, check it !' )
    rm(list = ls(all.names = TRUE)); closeAllConnections()
  }
  
  switch(s,
     HDI = { setnames( dt, c('DateTime','accelX') ) },
     HDJ = { setnames( dt, c('DateTime','accelY') )
     },
     HDK = { setnames( dt, c('DateTime','accelZ') )
     },
     HDX = { setnames( dt, c('DateTime','MagX') ) },
     HDY = { setnames( dt, c('DateTime','MagY') )
     },
     HDZ = { setnames( dt, c('DateTime','MagZ') )
     },
     { print('File is not indentified!') }
  )

  dt
}

# Select folder by browsing
folder_aux = dirname(file.choose())
# The files are sorted in alphabetical order, on the full path if full.names = TRUE.
MTfiles = list.files(folder_aux,pattern="*.MT$",full.names=T,recursive=F)
MT_name = list.files(folder_aux,pattern="*.MT$",full.names=F,recursive=F)

# Read filtered data from MATLAB
library(R.matlab)
Acc_1Hz_MATLAB = readMat( paste0(folder_aux,'/HD.mat') )

# Join HDI, HDJ, HDK together
# BE CAREFUL, need to check order AND number of rows !!
# First, join HDI00000, HDI00001; and then Y and Z
hdi = MAT_to_data_table(Acc_1Hz_MATLAB,'HDI')
hdj = MAT_to_data_table(Acc_1Hz_MATLAB,'HDJ')
hdk = MAT_to_data_table(Acc_1Hz_MATLAB,'HDK')

# Since AccelometerX (hdi), AccelometerY (hdj), AccelometerZ (hdk) may start and end at different times
# we need to fit them together at same start & end times
time_start = max(hdi$DateTime[1],hdj$DateTime[1],hdk$DateTime[1])
time_end   = min(hdi[.N]$DateTime,hdj[.N]$DateTime,hdk[.N]$DateTime)
hdi = hdi[DateTime >= time_start & DateTime <= time_end] 
hdj = hdj[DateTime >= time_start & DateTime <= time_end]
hdk = hdk[DateTime >= time_start & DateTime <= time_end]
hdi[,accelY:=hdj$accelY]; hdi[,accelZ:=hdk$accelZ]
setcolorder(hdi,c('DateTime','accelX','accelY','accelZ'))
rm(hdj); rm(hdk)

#### 1c) Join Accelometer (hdi), Magnetometer (hdx), TDR data (hdp) together ####
#### All is 1 Hz, need to change if prefer different resolution

# FIR: https://www.minidsp.com/applications/dsp-basics/fir-vs-iir-filtering
# Decimation: https://dspguru.com/dsp/faqs/multirate/decimation/
dat = hdi
setnames(dat, c('DateTime','AccX','AccY','AccZ') )
# Convert unit from mg to g (standard gravity g = 9.80665 m/s^2)
dat[,AccX := AccX/1000]; dat[,AccY := AccY/1000]; dat[,AccZ := AccZ/1000]

# Remove static component of accelerometer data using moving average
time_mv = 2 # Moving average time of 2 seconds
# See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050556#s2
dat[,Ax := AccX-frollmean(AccX,1+time_mv)]
dat[,Ay := AccY-frollmean(AccY,1+time_mv)]
dat[,Az := AccZ-frollmean(AccZ,1+time_mv)]

# ODBA and VeDBA calculated as in the article
dat[,ODBA  := abs(Ax) + abs(Ay) + abs(Az)]
dat[,VeDBA := sqrt(Ax^2 + Ay^2 + Az^2)]

# Remove NAs created by moving average in ODBA/VeDBA
dat = na.omit(dat)
# Only keep relevant data 
dat[,c(2:7):=NULL]
rm(hdi)

#### 2a) Detect "wrong" whale's positions ####
## Some whale positions is on land, so we find these positions

# Read positions from a data of positions of all whales
whale_LonLat = fread('Data/Fieldwork 2018/nar_pos_all_2018.txt')
whale_LonLat[,Datetime:=NULL]
setnames(whale_LonLat,c('Lat','Lon','Ind','DateTime'))
setcolorder(whale_LonLat, c('Ind','Lon','Lat','DateTime'))

# This coast line is fixed using QGIS (not all), following advice of Karl Z.
coast = fread('coast/Land_EG-nodes.csv')
coast$x = as.numeric(coast$x); coast$y = as.numeric(coast$y)
att = fread('coast/Land_EG-att.csv')
setnames(att, c('id','Zone_ID','TYPE','KODE'))
coast[att, Zone_ID:=Zone_ID, on=('shapeid==id') ]
coast[,shapeid:=NULL]
setnames(coast, c('X','Y','Zone_ID') )
# Convert to UTM coordinates to calculate the distances more easily
whale_UTM = LonLat_to_UTM( whale_LonLat$Lon,whale_LonLat$Lat,24 )
whale_UTM[,DateTime := whale_LonLat$DateTime]
whale_UTM[,No := seq_len(.N)]
setnames( whale_UTM, c( 'X','Y','DateTime','No' ) )

# By looking on QGIS, we identify "visually" what zones are close to whales
# Do it to make it faster
# TODO: make procedure automatically
coast_ID = unique(coast$Zone_ID)
nearest_ID = setdiff( coast_ID, c(coast_ID[coast_ID > 45300], 45158,45159,45160,45235) )

Shores = vector("list", length=length(nearest_ID))
names(Shores) = nearest_ID
for (i in seq_len(length(Shores)) ) {
  Shores[[i]] = coast[Zone_ID==nearest_ID[i] ]
}

# this allows to detect outliers, due to GPS errors or coast errors
# then put them inside dataset of whale positions
library(Rcpp); library(RcppArmadillo); sourceCpp('func.cpp')
x = find_wrong_position(whale_UTM[,c('X','Y')], Shores, nearest_ID)
whale_UTM[,wrong_pos:=x$d]
whale_UTM[,wrong_id:= ifelse(x$ID!=0, x$ID, 0) ]
which(x$d < 0)


#### 2b) Dead reckoning based on interpolation ####

# extract_whale_name: Extract the name of whale from directory path
# Input:
# - the path to the directory path having the dataset
# Output:
# - the whale name in the path
extract_whale_name = function(path) {
  library(stringr)
  w = str_locate_all(path,'/')[[1]][,1]; x = length(w)
  w = str_sub(path,w[x-1]+1,w[x]-1)
  if (str_detect(w,' ')) # if found space symbol
    return( str_sub(w,6,str_length(w)) ) # assumme that there id has 4 character and 1 space
  w
}

# Extract the name of whale from directory path
# Based on the structure of data, it lies between the last splash and the one before
whale_name = extract_whale_name(folder_aux)
# Only extract the position of the wanted whale
whale_LonLat = whale_LonLat[Ind == whale_name]
whale_LonLat[,Ind := NULL]

# Limit to the range of dataset 'dat' above
whale_LonLat[,DateTime := ymd_hms( whale_LonLat$DateTime )]
whale_LonLat = na.omit(whale_LonLat) # nrow(whale_LonLat) : 2033 rows

## Limit the dataset to pre-defined onset times belows 
#         Ind                  V1
# 1:   Asgeir 2018-08-30 04:24:16
# 2: Frederik 2018-08-29 17:01:00
# 3:    Helge 2017-08-19 11:48:36
# 4:    Kyrri 2018-09-01 00:12:29
# 5:     Nemo 2018-08-29 12:25:12
# 6:    Siggi 2018-08-29 04:32:06
switch(whale_name,
       Asgeir = {start_time = '2018-08-30 04:24:16'},
       Frederik = {start_time = '2018-08-29 17:01:00'},
       Helge18 = {start_time = '2018-09-01 00:34:30'},
       Kyrri = {start_time = '2018-09-01 00:12:29'},
       Nemo = {start_time = '2018-08-29 12:25:12'},
       Siggi = {start_time = '2018-08-29 04:32:06'},
)
start_time = ymd_hms(start_time)
st_time = whale_LonLat$DateTime[ max(which(whale_LonLat$DateTime >= start_time)[1] - 1, 0) ]
end_time = whale_LonLat$DateTime[which(whale_LonLat$DateTime > ymd_hms('2018-09-02 00:00:00'))[1]]
whale_LonLat = whale_LonLat[ DateTime >= st_time & DateTime <=  end_time ]

# Convert from Lon/Lat into UTM
# Here use zone 24N
whale_UTM = LonLat_to_UTM( whale_LonLat$Lon,whale_LonLat$Lat )
whale_UTM[,DateTime := whale_LonLat$DateTime]
whale_UTM[,No := seq_len(.N)]

# Extract first longitude and latitude
whale_1st_lon = whale_UTM$Lon[1]; whale_1st_lat = whale_UTM$Lat[1]

# Function to interpolate whale's dead reckoning at every second
# compute how far the whale move 
# in each direction X and Y horizontally every second
# 
# Input (Default projection is zone 24N): 
# - whale_LonLat: data.frame contains original Longtitude & Latitude of FGPS
# - whale_1st_lon: first longitude position in UTM
# - whale_1st_lat: first latitude position in UTM
#
# Ouput:
# - a dataframe consisting interpolated positions and speed at every second
track_simulation = function(whale_pos,LON0,LAT0)  {
  
  posFGPS = copy(whale_UTM)
  # Movement distances (X,Y axis) between time intervals showed in whale_UTM
  n = nrow(posFGPS)
  disX = diff( posFGPS$Lon ); disY = diff( posFGPS$Lat )
  delta = diff( unclass(whale_UTM$DateTime) ) # time intervals between FGPSs (secs)
  # Average speed
  V_cor.X = disX/delta # m/s
  V_cor.Y = disY/delta # m/s
  whale_UTM[,diffTime := round(c(delta/60,0),digits=2) ] # conver to Minutes
  
  px = py = 0
  for ( j in (1:(n-1)) ) {
    x1 = rep( V_cor.X[j], delta[j] )
    y1 = rep( V_cor.Y[j], delta[j] )
    px = c(px,x1); py = c(py,y1) 
  }
  dt0 = data.table(DateTime = seq(posFGPS[1]$DateTime,posFGPS[.N]$DateTime, by='secs'),
                   X = cumsum(px), Y = cumsum(py) )
  dt0[,`:=`(V_cor.X = px,V_cor.Y = py) ]  # Meter
  dt0[,`:=`(X=X+LON0,Y=Y+LAT0)]
  whale_simulated_path = dt0[,c('DateTime','X','Y','V_cor.X','V_cor.Y')]  
  whale_simulated_path
}

# interpolate the new data of 1Hz resolution
whale_LonLat[,DateTime:=ymd_hms(DateTime)]
x = seq(whale_LonLat[1]$DateTime,whale_LonLat[.N]$DateTime,by='secs')
track = track_simulation(whale_UTM,whale_1st_lon,whale_1st_lat)
# Convert back to Longtitude/Latitude coordinates to put on the original data
L = UTM_to_LonLat(track$X,track$Y)
track[,`:=`(Lon=L$Lon,Lat=L$Lat)]

# TEST: if there is any NA 
options(show.error.locations = TRUE)
if (anyNA(track)) {
  stop( "NAs exist in 'track', fix it before continuing !!" ) }
whale_UTM[track,`:=`(X0=X,Y0=Y), on=.(DateTime)]

library(geosphere)
whale_LonLat[, diffDist:=distGeo( whale_LonLat[,c('Lon','Lat')],
                                  whale_LonLat[c(2:.N,.N),c('Lon','Lat')] )]
whale_LonLat[, diffDist_sec:=diffDist/c(diff(as.numeric(whale_LonLat$DateTime)),1)]


#### 3a) Calculate distance from whale to Scorebysund coast ####

## Write a small sample of each 10 seconds to 
## estimate what regions (by their IDs) are close the whale 
# fwrite(track[seq(1,.N,by=10),c('X','Y') ],
#        paste0('Data/Fieldwork 2018/',whale_name,'_10sec.csv'))
##
## Then map it into QGIS to choose 'visually' such regions (i.e. polygons)
## This semi-automatic way is nessesary for now 
## because don't have time now to make it faster automatically,
## it's too slow to calculate distances to all polygons
## https://stackoverflow.com/q/52423060
## https://stackoverflow.com/q/19892564

# OLD METHOD: 
# Read full Greenland and extract only Scorebysund
#
# ID of regions of Scoreby Sund
# ID_Scorebysund = c('45151','45152','45153','45154','45155','45156','45157','45158','45159','45160','45168','45175','45176','45177','45178','45179','45180','45181','45182','45183','45184','45185','45186','45187','45188','45189','45190','45231','45232','45233','45234','45235','45236','45237','45238','45239','45240','45241','45242','45243','45244','45245','45246','45247','45248','45249','45311','45312','45313','45314','45315','45316','45317','45318','45319','45320','45321','45322','45323','45324','45325','45326')
# coast = fread('coast/Land-nodes.csv')
# coast$x = as.numeric(coast$x); coast$y = as.numeric(coast$y)
# att = fread('coast/Land-att.csv')
# setnames(att, c('id','Zone_ID','TYPE','KODE'))
# coast[att, Zone_ID:=Zone_ID, on=('shapeid==id') ]; coast[,shapeid:=NULL]
# setnames(coast, c('X','Y','Zone_ID') )
# coast = coast[Zone_ID %in% ID_Scorebysund]


# Make a dataframe of nearest zones of each whale
# https://stackoverflow.com/a/26803445
nearest_zone = fread('Data/Fieldwork 2018/nearest_zone.csv')
nearest_zone_ID = as.numeric( unlist(strsplit(nearest_zone[Whale==whale_name]$Zone_ID, '_')) ) # Different for each whale
shore = coast[Zone_ID %in% nearest_zone_ID]
shore[,Zone_ID:=as.factor(Zone_ID)]


# Sampling each 'sampling_time' seconds because too slow for each second
# Generate almost 'regular' sequences, i.e
# if 'to'-'from' divides by 'by', then it returns regular sequence
# otherwise, return regular sequence but include the last number 'to'
seq_include_last = function(from,to,by){
  if (to < from + by) 
    return( c(from,to) )
  s = seq(from,to,by)
  return( c(s[-length(s)],to)) 
}
sampling_time = 1 # second
track[, DateTime:=as.numeric(DateTime) ]
whale_partial_path = track[seq_include_last(1,.N,by=sampling_time)]
coast_distance_file = paste0(whale_name,'_distance_to_coast_',
                             sampling_time,'_sec.rds')
Shores = vector("list", length=length(nearest_zone_ID))
names(Shores) = as.character(nearest_zone_ID)
for (i in seq_len(length(Shores)) ) {
  Shores[[ as.character(nearest_zone_ID[i]) ]] =
    shore[ Zone_ID==nearest_zone_ID[i] ]
}


## Now compute to the distances of whales positions to closest coasts
## Check if file exist: if yes, read it; if not, create it
if (!file.exists(paste0('Data/Fieldwork 2018/',coast_distance_file))) {
  
  ptm = proc.time()
  dt = path_to_coast(whale_partial_path, Shores,
                     as.character(nearest_zone_ID) ) ; setDT(dt)
  print(proc.time() - ptm)  
  # Write output
  saveRDS(dt, paste0('Data/Fieldwork 2018/',coast_distance_file))
} else {
  dt = readRDS(paste0('Data/Fieldwork 2018/',coast_distance_file))
}


## TEST: to see unrealistic results
which(diff(dt$Distance)>6*sampling_time)
range(diff(dt$Distance))

# add distances to track dataset
dt[,DateTime:=as_datetime(DateTime)]
dt[,Zone_ID:=as.character(Zone_ID)]
track[,DateTime:=as_datetime(DateTime)]
track[dt, DistanceShore:=Distance, on=.(DateTime) ]
track = track[ DateTime >= start_time & DateTime < ymd_hms('2018-09-02 00:00:00'),]
# track[dat, `:=`(ODBA=ODBA,VeDBA=VeDBA), on=.(DateTime) ]
# Write the ODBA/VeDBA to a sepate file
fwrite(dat, file=paste0('Data/Fieldwork 2018/',whale_name,'_DBA.csv.gz'))

#### 4. Distance to ship in LINE of SIGHT (LoS) ####
## Calculate the distances between whales and ships
## LoS: the duration where the whales can "see" the ship, 
## i.e. no lands between them

## Read ship positions
ship = fread('Data/Fieldwork 2018/all_LAKO_pos.csv')
ship[, c('N','E'):=NULL ]
ship[, Datetime:=dmy_hms(Datetime) ]
x = LonLat_to_UTM(ship$Long,ship$Lat)
ship[, X_24N:=x$Lon]; ship[, Y_24N:=x$Lat]

## Read trials & LoS info
library(readxl)
Concat_date_time = function(d,t) {
  ymd_hms(paste(as.character(d), hour(t), minute(t), second(t))) }
trial = read_excel('Data/Fieldwork 2018/Overview_Acou_HTR_GPS_Seismik_LOS.xlsx',
                   sheet = 'LOS + TRIALS'); setDT(trial)
trial = na.omit(trial,cols=c('LOS')) # remove if not Light Of Sight
trial = trial[ Ind == unique(trial$Ind)[sapply(unique(trial$Ind), 
                                      function(x) str_detect(whale_name,x))] ]

options(warn=0) # turn off error mode, due to concat NAs
trial[,StartSei:=Concat_date_time(Start_sei,StarTime_sei)]
trial[,EndSei:=Concat_date_time(End_sei,EndTime_sei)]
trial[,c("Start_sei","StarTime_sei","End_sei","EndTime_sei"):=NULL]

## Mark when LoS & Seismic start and end in 'track' data
track[, LoS:=0]; track[, Sei:=0]
track[, DateTime:=as_datetime(DateTime) ]
track[, c('V_cor.X','V_cor.Y'):=NULL ]
for (i in seq_len(nrow(trial))) {
  track[DateTime >= trial$Start[i] & DateTime <= trial$End[i], LoS:=1 ]
  track[DateTime >= trial$StartSei[i] & DateTime <= trial$EndSei[i], Sei:=1]
}
track[ship, ship_X:=X_24N, on=('DateTime==Datetime')  ]
track[ship, ship_Y:=Y_24N, on=('DateTime==Datetime')  ]
track[ LoS==1, Dist_Ship:=sqrt((X-ship_X)^2 + (Y-ship_Y)^2) ]
track[ship, Lon_ship:=Long, on=('DateTime==Datetime')  ]
track[ship, Lat_ship:=ship$Lat, on=('DateTime==Datetime')  ]
# Then calculate the distances between ship and whales using Longtitude/Latitude coordinates
track[, Dist_Ship_LonLat:=distGeo(track[,c('Lon','Lat')],
                                  track[,c('Lon_ship','Lat_ship')]) ]
track[ LoS!=1, Dist_Ship_LonLat:=NA ]
track[LoS==1,diff_dist:= Dist_Ship_LonLat-Dist_Ship]

# if there is any missing positions of ship, interpolate by spline
if ( anyNA(track[LoS==1]$Dist_Ship) ) {
  library(zoo)
  track[LoS==1, Dist_Ship:=na.spline(track[LoS==1]$Dist_Ship)]
  track[LoS==1, Dist_Ship_LonLat:=na.spline(track[LoS==1]$Dist_Ship_LonLat)]
}
track[, `:=`(ship_X=NULL,ship_Y=NULL)  ] 

# Finally write the final dataset
fwrite(track, file=paste0('Data/Fieldwork 2018/',whale_name,'_dataset.csv.gz'))

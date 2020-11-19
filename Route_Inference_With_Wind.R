#' This script can be used to infer migration routes of birds tagged with geolocators.
#' Before this code is used, make sure the database contains the wind data of 
#' 1. the time span you are analyzing
#' 2. the geographical region you are analyzing
#' Authors: Mike Werfeli, Peter Ranacher, Felix Liechti and Martins Briedis




#Load required packeages
library(GeoLocTools)
setupGeolocation()
library(GeoLight)
library(SGAT)
library(lubridate)
library(maptools); data(wrld_simpl)
library(rworldmap); newmap <-getMap(resolution = "low")
library(RColorBrewer)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
require(maps)
library(sp)
require(viridis)
library(RPostgreSQL)
library(PAMLr)
library(TwGeos)
library(SGAT)
library(zoo)
library(geodist)
library(psych)

# Se up functions needed for Inference (Priors)

##' Function to generate a land mask in raster form. This will be used to push bird positions on land rather than on sea.
##' @param xlim,ylim parameters which define map extent
##' @return raster indicating probabilities. (sea very low land 1)
earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE, index) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1
  rs = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
             rasterize(wrld_simpl, r, 1, silent = TRUE), 
             rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  # make the movement raster the same resolution as the stationary raster, but allow the bird to go anywhere by giving all cells a value of 1
  rm = rs; rm[] = 1
  
  # stack the movement and stationary rasters on top of each other
  mask = stack(rs, rm)
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  mask = as.array(mask)[nrow(mask):1,,sort(unique(index)),drop=FALSE]
  
  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin), index)]
}

################## Sensor Data  Analysis ###########################


##' Function to read and plot sensor data (light measurements)
##' @param date dates on which sensor was operating
##' @return plot of light data
sensorIMG_2 <- function (date, sensor_data, tz = "UTC", plotx = TRUE, ploty = TRUE, 
                         labelx = TRUE, labely = TRUE, offset = 0, dt = NA, xlab = "Hour", 
                         ylab = "Date", cex = 2, col = c("black", viridis::viridis(90)), 
                         ...) {
  dts <- c(5, 10, 15, 20, 30, 60, 90, 120, 180, 240, 300, 360, 
           400, 480, 540, 600, 720, 900, 960, 1200)
  if (is.na(dt)) {
    dt <- mean(diff(as.numeric(date)))
    dt <- dts[which.min(abs(dt - dts))]
  }
  tmin <- .POSIXct(as.POSIXct(as.Date(date[1])) + offset * 
                     60 * 60, tz)
  if (as.numeric(tmin) > as.numeric(date[1])) 
    tmin <- tmin - 24 * 60 * 60
  tmax <- .POSIXct(as.POSIXct(as.Date(date[length(date)])) + 
                     offset * 60 * 60, tz)
  if (as.numeric(tmax) < as.numeric(date[length(date)])) 
    tmax <- tmax + 24 * 60 * 60
  sensor_data1 <- approx(as.numeric(date), sensor_data, seq(as.numeric(tmin) + 
                                                              dt/2, as.numeric(tmax) - dt/2, dt))$y
  
  
  
  m <- 24 * 60 * 60/dt
  n <- length(sensor_data1)/m
  sensor_data2 <- matrix(sensor_data1, m, n)
  
  sensor_data3 <- cbind(sensor_data2[-nrow(sensor_data2),],sensor_data2[-1,])
  sensor_data3 <- rbind(sensor_data2[,-ncol(sensor_data2)],sensor_data2[,-1])
  
  hour <- seq(offset, offset + 48, length = m*2 )
  day <- seq(tmin, tmax, length = n-1 )
  rotate <- function(x) t(apply(x, 2, rev))
  
  
  image(hour, as.numeric(day), rotate(t(sensor_data3)), col = col, 
        axes = F, xlab = "", ylab = "", ...)
  
  
  if (plotx) {
    axis(1, at = seq(0, 48, by = 4), labels = seq(0, 48, 
                                                  by = 4)%%24, cex.axis = cex)
  }
  if (labelx == TRUE) 
    mtext(xlab, side = 1, line = 3, cex = cex)
  if (ploty) {
    axis.POSIXct(4, at = seq(tmin, tmax, length = n/60), 
                 labels = as.Date(seq(tmax, tmin, length = n/60)), 
                 las = 1, cex.axis = cex)
  }
  if (labely == TRUE) 
    mtext(ylab, side = 2, line = 1.2, cex = cex)
  box()
}

########### set up parameters needed #######

#Where is data stored
folderName<- "Folder where bird data is located"
birdName<- "name of the bird"
Data <- paste0("Enter Directory here",folderName)
# list all light files
ID.list<-list.files(paste0("Enter Directory here", folderName),pattern=".glf",recursive = T) # Change to gle or glf - the one you use for analyses
ID <- ID.list[1]
ID

glf_file<- paste0("Enter Directory here",folderName, "/", folderName, ".glf")
twl_CSVfile<- paste0("Enter Directory here",folderName,"/", folderName,"_twl.csv")
PAM_filepath<- paste0("Enter Directory here",folderName,"/")
Activity_path<- "Enter Directory here"

#parameter for PAM analysis

# Start date of observation period --> used to crop pam data
# make sure the cropping period is in the correct date format
start = as.POSIXct("2015-07-15","%Y-%m-%d", tz="UTC")

# End date of observation period --> used to crop pam data
end = as.POSIXct("2016-04-12","%Y-%m-%d", tz="UTC") 

# Parameter for Database connection
user<- "username"
password<- "Password"
dbname<- "db name"
host<- "db host"

#Parameter for result storage
Route_name<- 'WM'
Storage_path_route<-paste0('Enter Directory here','.png')
LatAndLonPlot_Path<- paste0('Enter Directory here','.png')
Summaryfile_Path<- paste0("Enter Directory here",".csv")


WindSpeed_Path<- paste0('Enter Directory here','.png')
WindSupport_Path<- paste0('','.png')

WindData_Path<- paste0("",".csv")
WindSpeedCSV_Path<- paste0("Enter Directory here","_WindSpeeds.csv") 
WindSuportCSV_Path<- paste0("Enter Directory here","_WindSuports.csv")

################ Data preparation ###################


#set up lat und lon start pts
lon.calib <- 7.368
lat.calib <- 46.234

# Light data analysis and calibration

offset <- 12 

# Load raw data
raw <- glfTrans(glf_file)
names(raw) <- c("Date", "Light")
raw$Light[raw$Light<0] <- 0 
raw$Light  <- log(raw$Light+0.0001) + abs(min(log(raw$Light+0.0001)))


twl <- read.csv(twl_CSVfile)
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "UTC")
twl <- twl[!twl$Deleted,]
head(twl)


# clipping raw data to relevant extent  
raw <- subset(raw, Date>=min(twl$Twilight) & Date<=max(twl$Twilight)) 
head(raw)

# SGAT Group model----

# 1. Calibration----
lightImage( tagdata = raw,
            offset = offset,     
            zlim = c(0, max(raw$Light)))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")

tm.calib <- as.POSIXct(c("2015-07-15 00:00", "2015-09-10 00:00", "2016-03-12 00:00", "2016-04-12 00:00"), tz = "UTC",format="%Y-%m-%d %H:%M") # Selecting calibration period(s) 
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange") # chech if them make sense on the light graph

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2] | Twilight>=tm.calib[3] & Twilight<=tm.calib[4])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

zenith  <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Alternatively use Hill-Ekstrom calibration from the wintering sites
startDate <- "2015-09-10"
endDate   <- "2016-03-12"

start = min(which(as.Date(twl$Twilight) == startDate))
end = max(which(as.Date(twl$Twilight) == endDate))

(zenith_sd <- findHEZenith(twl, tol=0.05, range=c(start,end)))

lightImage( tagdata = raw,
            offset = offset,     
            zlim = c(0,5))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, zenith = zenith_sd,lwd = 2, col = "orange") # Visually check if new calibration-light curves fit with the light-dark data

# 2. Define groups of sunevents = stationary positions
# Defining movement and stationary periods----

library(PAMLr)
library(dplyr)

# PAM Data import
ID.list = list.files(paste0("Enter Directory here"),include.dirs=T)
ID.list
ID = ID.list[3]
ID

pathname = paste0("Enter Directory here/")
measurements = c(".pressure", 
                 ".glf",
                 ".acceleration", 
                 ".temperature")

PAM_data = importPAM(pathname, measurements)

# Cropping the data

# make sure the cropping period is in the correct date format
start = as.POSIXct("2015-07-15","%Y-%m-%d", tz="UTC")
# end = max(PAM_data$light$date)
end = as.POSIXct("2016-04-12","%Y-%m-%d", tz="UTC") # 20-June for full data# Crop the data
PAM_data= cutPAM(PAM_data,start,end)

length(PAM_data$pressure$date)
altitude = altitudeCALC(P = PAM_data$pressure$obs)
plot(PAM_data$pressure$date[2:17000], PAM_data$pressure$obs[2:17000], type="o",pch=16, xlab="Date", ylab="Air Pressure (hPa)", main = "Air Pressure en route")
alt_dta<- cbind(PAM_data$pressure$date, PAM_data$pressure$obs)



# Classify based on pressure
twl2 = GeoLight::twilightCalc(PAM_data$light$date, PAM_data$light$obs,
                              LightThreshold = 2, 
                              ask = FALSE)
# TOclassify = reducePAM(PAM_data , "pressure", interp = TRUE, summary="median")
success <- FALSE

#This will repeat until the assigned groups (clusters) are in correct order with lowest number for lowest-level activity, for example, 1-restig, 2-active, 3-migrating, 
while(!success){
  
  TOclassify = pamPREP(dta = PAM_data,
                       method= "pressure",
                       Pdiff_thld = 1.4,
                       # light_thld = 2, 
                       twl = twl2,
                       interp = F)
  
  # Get rid of any extra NAs 
  TOclassify = TOclassify[complete.cases(TOclassify),]
  
  
  # Select the columns to use as predictors in the model
  predictor = TOclassify[, c(#"duration",
    "night_P_diff",
    # "median_pitch",
    "total_daily_P_change",
    "sum_activity",
    "cum_pressure_change",
    "cum_altitude_change"
    # "sd_temp"
  )]
  
  
  # Perform the classification
  classification = classifyPAM(predictor,
                               states=2, 
                               method = "hmm")
  
  # Convert the events to data
  pressure_classification = classification2PAM(from = TOclassify$start,
                                               to = TOclassify$end,
                                               classification = classification$cluster,
                                               addTO = PAM_data$pressure)
  pressure_classification = data.frame(cluster=pressure_classification)
  # Find the pressure difference for each state
  P_state1 = median(TOclassify$total_daily_P_change[classification$cluster == 1])
  P_state2 = median(TOclassify$total_daily_P_change[classification$cluster == 2])
  P_states = c(P_state1, P_state2)
  # Allocate the state with the highest pressure difference to migration
  Mig_state = which(P_states == max(P_states))
  # Now add this information to the classification
  pressure_classification$states = pressure_classification$cluster
  pressure_classification$states[pressure_classification$cluster == Mig_state] = "Migration"
  pressure_classification$states[pressure_classification$cluster != Mig_state] = "Active"
  pressure_classification$states[is.na(pressure_classification$cluster)] = "Unclassified"
  # Associate each state with a color
  pressure_classification$colour = pressure_classification$cluster
  pressure_classification$colour[pressure_classification$cluster == Mig_state] = "orange"
  pressure_classification$colour[pressure_classification$cluster != Mig_state] = "royalblue3"
  pressure_classification$colour[is.na(pressure_classification$cluster)] = "black"
  # Store for later
  endurance_classification = pressure_classification
  
  pressure_classification$cluster = pressure_classification$cluster+1
  pressure_classification$cluster[is.na(pressure_classification$cluster)]=1
  
  success <- unique(pressure_classification$cluster)[2] - unique(pressure_classification$cluster)[1] == 1
  print(success)
}

axis_scale<-format(as.Date(ifelse(as.numeric(format(as.Date(seq(PAM_data$acceleration$date[1],max(PAM_data$acceleration$date),by=60*60*24)), format = "%d"))==1,
                                  as.Date(seq(PAM_data$acceleration$date[1],max(PAM_data$acceleration$date),by=60*60*24)),NA),origin = "1970-01-01"),"%e-%b")
axis_scale2<-as.Date(ifelse(as.numeric(format(as.Date(seq(PAM_data$acceleration$date[1],max(PAM_data$acceleration$date),by=60*60*24)), format = "%d"))==1,
                            as.Date(seq(max(PAM_data$acceleration$date),PAM_data$acceleration$date[1],by=-60*60*24)),NA),origin = "1970-01-01")
geo_twl <- export2GeoLight(twl)
t.lub <- ymd_hms(geo_twl$tFirst[geo_twl$type==1])
rise <- hour(t.lub) + minute(t.lub)/60
t.lub2 <- ymd_hms(geo_twl$tFirst[geo_twl$type==2])
set<-hour(t.lub2) + minute(t.lub2)/60

par( mfrow= c(1,1), oma=c(0,0,0,5))
col=c("black","royalblue4","gold")
sensorIMG_2(PAM_data$pressure$date, ploty = F,
            pressure_classification$cluster, main = paste0("Activity ",substr(ID,1,4)),
            col=col, cex=1.2, cex.main = 2)
axis.POSIXct(4, at = rev(axis_scale2),
             labels = rev(axis_scale[!is.na(axis_scale)]), 
             las = 1, 
             cex.axis = 1)
lines(rev(geo_twl$tFirst[geo_twl$type==1])~rise,col=0)
lines(rev(geo_twl$tFirst[geo_twl$type==2])~set,col=0)
lines((rise+24),rev(geo_twl$tFirst[geo_twl$type==1]),col=0)
lines((set+24),rev(geo_twl$tFirst[geo_twl$type==2]),col=0)

pressure_classification$date <- PAM_data$pressure$date
migration <- pressure_classification[pressure_classification$states=="Migration",]
migration$day <- as.Date(migration$date)


for(i in 2:nrow(twl)){
  twl$group2[i] <- !any(pressure_classification$cluster[which(PAM_data$pressure$date %within% interval(twl$Twilight[i-1], twl$Twilight[i]))]==3)
}
twl$group2[1]=1


# Create a vector which indicates which numbers sites as 111123444444567888889
g <- c(1)
for (i in 2:length(twl$group2)){
  if(twl$group2[i]==1){g<-c(g,g[length(g)])}else{
    g<-c(g,g[length(g)]+1)}
}

gr <- g

# Check how noon times of stationary periods look like - good way to visually validate if stationary-movment designation makes sense
sites<-vector()
for (i in 1:max(gr)){
  a <- length(which(gr==i))
  if (a == 1){sites<-c(sites,0)}else{sites<-c(sites,rep(i,a))}
}
for (j in 1:length(unique(sites[sites>0]))){
  a<-unique(sites[sites>0])
  sites[which(sites==a[j])]=j
}
sites<-sites[1:(length(sites)-1)]

behav   <- ifelse(gr%in%c(1:max(gr))[as.data.frame(table(gr))[,2]>1], TRUE, FALSE)
sitenum <- ifelse(behav, gr, 0)
sitenum <- sitenum[sitenum==0 | !duplicated(sitenum)]
sitenum[sitenum>0]=1:length(sitenum[sitenum>0])

stationary <- sitenum>0
twl$group<-gr

geo_twl <- export2GeoLight(twl)
for(i in 1:nrow(geo_twl)){
  geo_twl$noon[i]<-median(c(geo_twl$tFirst[i],geo_twl$tSecond[i]))
}
geo_twl$noon <- as.POSIXct(geo_twl$noon, tz = "GMT",origin = "1970-01-01")
t.lub <- ymd_hms(geo_twl$noon[geo_twl$type==1])
h.lub <- hour(t.lub) + minute(t.lub)/60
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2) )
jd = insol::JD(geo_twl$noon[geo_twl$type==1])
plot((h.lub+insol::eqtime(jd)/60)~geo_twl$noon[geo_twl$type==1],
     col=ifelse(sites[seq(2, length(sites), 2)]==0,1,sites[seq(2, length(sites), 2)]),   # as a rule of thumb - if a streight horizontal line can be drawn 
     pch=ifelse(sites[seq(2, length(sites), 2)]==0,21,16),                               # through each of the differently coloured points, your stationary periods are ok
     xlab="time",main="Midday by stationary periods",ylab="time of day")                 # if there is no difference along the y-axis between consecutive colours - those stationary periods can be suspecious (may belong to one stationary site) and maybe should have been merged


# simle checkup of daylength around equinox and winter solstice
geo_twl$daylength <- difftime(geo_twl$tSecond,geo_twl$tFirst)
mean(geo_twl$daylength[geo_twl$type==1 & geo_twl$tFirst>as.POSIXct("2015-09-15") & geo_twl$tFirst<as.POSIXct("2015-09-25")])
mean(geo_twl$daylength[geo_twl$type==1 & geo_twl$tFirst>as.POSIXct("2015-12-15") & geo_twl$tFirst<as.POSIXct("2015-12-25")])



# This is an important next step which is new in comparsion to the earlier methods of geolocation... 
# get flight duration
ID.list2 = list.files(Activity_path,pattern="_act",include.dirs=T) # this is the PAMLr output file of flight classification
ID.list2
ID2 = ID.list2[2]
ID2

dta<-read.csv(paste0(Activity_path,ID2))
dta$date <- as.Date(dta$start,tz="UTC" ) #"%m/%d/%Y"
dta$start <- as.POSIXct(dta$start, tz="UTC") #"%m/%d/%Y %H:%M",
dta$end <- as.POSIXct(dta$end,tz="UTC")     #"%m/%d/%Y %H:%M",
dta$logger_ID <- substring(ID,1,4)

# Subset only the relevant variables
alt <- dta

t_idx<- numeric(0)
twl$flight_duration <- 0
for( i in which(twl$group2==0)){
  t_idx<- c(t_idx, i)
  a <- which(alt$start %within% interval(twl$Twilight[i-1], twl$Twilight[i]))
  b <- which(alt$end %within% interval(twl$Twilight[i-1], twl$Twilight[i]))
  
  if (length(a)==1 & length(b)==1 & max(a)==max(b)) {twl$flight_duration[i] <- difftime(alt$end[b],alt$start[a],units = "hours")}
  if (length(a)>1 & length(b)>1 & max(a)==max(b)) {twl$flight_duration[i] <- sum(difftime(alt$end[b],alt$start[a],units = "hours"))}
  if (length(a)==1 & length(b)==0) {twl$flight_duration[i] <- difftime(twl$Twilight[i],alt$start[a],units = "hours")}
  if (length(a)==0 & length(b)==1) {twl$flight_duration[i] <- difftime(alt$end[b],twl$Twilight[i-1],units = "hours")}
  if (length(a)>length(b)) {twl$flight_duration[i] <- sum(difftime(alt$end[b],alt$start[b],units = "hours"),difftime(twl$Twilight[i],alt$start[max(a)],units = "hours"))}
  if (length(a)<length(b)) {twl$flight_duration[i] <- sum(difftime(alt$end[min(b)],twl$Twilight[i-1],units = "hours"),difftime(alt$end[a],alt$start[a],units = "hours"))}
  if (length(a)==1 & length(b)==1 & max(a)!=max(b)) {twl$flight_duration[i] <- sum(difftime(alt$end[b],twl$Twilight[i-1],units = "hours"),difftime(twl$Twilight[i],alt$start[a],units = "hours"))}
  if (length(a)>1 & length(a)==length(b) & max(a)>max(b)) {twl$flight_duration[i] <- sum(difftime(alt$end[min(b)],twl$Twilight[i-1],units = "hours"),
                                                                                         difftime(alt$end[a[which(a%in%b)]],alt$start[a[which(a%in%b)]],units = "hours"),
                                                                                         difftime(twl$Twilight[i],alt$start[max(a)],units = "hours"))}
  
}
sum(twl$flight_duration)
sum(alt$Duration..h.)

#Get timestamps for route inference
wind_t<- vector("double", length(t_idx))
stopover_t<- vector("double", length(t_idx))
start_mig<- vector("double", length(t_idx))
end_mig<- vector("double", length(t_idx))
for (i in (1:length(t_idx))){
  t_starts<- twl$Twilight[t_idx[i]]
  start_mig[i]<- t_starts
  t_ends<- twl$Twilight[t_idx[i+1]]
  end_mig[i]<- t_ends
  t_end_stopover<- twl$Twilight[t_idx[i]-1]
  t_stamp_mig<- as.POSIXct(as.numeric((t_starts)+(as.numeric(t_ends)-as.numeric(t_starts))/2), origin="1970-01-01",tz="UTC")
  t_stamp_stopover<- as.POSIXct(as.numeric((t_starts)+(as.numeric(t_end_stopover)-as.numeric(t_starts))/2), origin="1970-01-01",tz="UTC")
  
  wind_t[i]<- t_stamp_mig
  stopover_t[i]<- t_stamp_stopover
}

wind_t[length(t_idx)]<- twl$Twilight[t_idx[length(t_idx)]+2]
stopover_t[length(t_idx)]<- twl$Twilight[t_idx[length(t_idx)]+2]
start_mig[length(t_idx)]<- twl$Twilight[t_idx[length(t_idx)]+2]
end_mig[length(t_idx)]<- twl$Twilight[t_idx[length(t_idx)]+3]
wind_t[length(t_idx)+1]<- twl$Twilight[t_idx[length(t_idx)]+6]
start_mig[length(t_idx)+1]<- twl$Twilight[t_idx[length(t_idx)]+6]
end_mig[length(t_idx)+1]<- twl$Twilight[t_idx[length(t_idx)]+6]
stopover_t[length(t_idx)]<- twl$Twilight[t_idx[length(t_idx)]+6]

#test plot
plot(PAM_data$pressure$date[2:17000], altitude[2:17000], type="o",pch=16, xlab="Date", ylab="Altitude (m)", main='Altitude en route')
alt_dta<- cbind(PAM_data$pressure$date, PAM_data$pressure$obs)



plot(PAM_data$pressure$date[2:17000], PAM_data$pressure$obs[2:17000], type="o",pch=16, xlab="Date", ylab="Altitude (m)")

#get flight time on each migratory bout
flt_t = twl$flight_duration[twl$flight_duration>0]

# Here starts the normal SGAT analyses
tol=0.08# adjust to a larger value to extrapolate more during the equinoq times

# This will draw your initial track - play around with the tol value to find something that looks ok-ish.
path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol=tol) # Here I use the zenith value from in-habitas calibration, alternatively replace zenith with zenith_sd for Hill-ekstrom calibration
x0 <- path$x
z0 <- trackMidpts(x0)
plot(x0, type = "n", xlab = "", ylab = "")
plot(wrld_simpl, col = "grey90", add = T)
points(x0, pch=19, col="cornflowerblue", type = "o")
points(lon.calib, lat.calib, pch = 16, cex = 2.5, col = "firebrick")
box()


x0 <- cbind(tapply(path$x[,1],twl$group,median), 
            tapply(path$x[,2],twl$group,median))

points(x0, pch=19, col="firebrick", type = "o")
points(x0, col = ifelse(stationary,rainbow(max(sitenum)),"firebrick"), pch = 16, cex = ifelse(stationary, 2, 0.1)) # larger dots = medians of stationary sites

#Set up connection to wind Data Base
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user=user, password=password, 
                 dbname = dbname, host = host)

# Define known locations and set them as fixed locations. Ficed locations are not changed in BI
fixedx <- rep_len(FALSE, length.out = nrow(x0))
fixedx[1] <- TRUE # first stationary site - fixed
fixedx[nrow(x0)] <- T # last stationary site - fixed
x0[fixedx,1] <- lon.calib
x0[fixedx,2] <- lat.calib
z0 <- trackMidpts(x0)

# Define movment model
# variable beta conatins alpha and beta parameter which are used to describe the gamma distribution of the birds air speed
# The mean air speed is computed as: alpha/beta
# In this case mean air speed is 50 km/h following Bruderer et. al., 2001; Bolch et.al., 1982; Liechti et. al. 1995
beta  <- c(30,.6)

matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h",
        main="movment model")

# Define landmask - stationary sites only on land, not allowed in the ocean
xlim <- range(x0[,1]+c(-5,5))
ylim <- range(0,55)

index = ifelse(stationary, 1, 2)

mask <- earthseaMask(xlim, ylim, n = 2,index=index)


## Define the log prior for x and z
logp <- function(p) {
  f <- mask(p)
  ifelse(!f | is.na(f), -1000, log(1))
}

# Start the Estelle Model

model <- groupedThresholdModel(twl$Twilight,
                               twl$Rise,
                               group = twl$group, #This is the group vector for each time the bird was at a point
                               twilight.model = "ModifiedGamma",
                               alpha = alpha,
                               beta =  beta,
                               x0 = x0, # meadian point for each greoup (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zenith0+.5, # Here I use the zenith0 value from in-habitat calibration. Using zenith or zenith_sd is not recommended here and will give you wrong results as they are medina SEAs
                               logp.x = logp, # land sea mask
                               dt = flt_t, # This is the new part here - movment is allowed only for the number of hours that PAM analyses said the bird was flying
                               fixedx = fixedx)
# define the error shape
x.proposal <- mvnorm(S = diag(c(0.025, 0.025)), n = nrow(x0))
z.proposal <- mvnorm(S = diag(c(0.025, 0.025)), n = nrow(z0))

##' Read wind data from data base
#Unique time stams
time_list <- dbGetQuery(con, paste0("SELECT DISTINCT tms FROM", 
                                    "\"WindTable\"", ";"))

#Find closest time stamp available in db for data retrieval
#compute step size in t stamps stored in db
step_size <- time_list$tms[2] - time_list$tms[1]
start_t_idx <- time_list$tms[1]/step_size

#find closest time stamp to start times of each migratory bout
start_mig_t_idx <- round(start_mig/step_size)
start_list_idx <- start_mig_t_idx - start_t_idx
tm_start <- time_list$tms[start_list_idx]

#find closest time stamp to end times of each migratory bout
end_mig_t_idx <- round(end_mig/step_size)
end_list_idx <- end_mig_t_idx - start_t_idx
tm_end <- time_list$tms[end_list_idx]

#Find all time stamps between each start and end time
alt_t_idx<- seq(1,length(tm_start))
tm<- lapply(alt_t_idx, FUN= function(x){
  time_list[t_span<-which(time_list>= tm_start[x] & time_list<= tm_end[x]),1]
})


#set up rounding function to be called in alt finder
mround <- function(x, base) {
  base * round(x/base)
}


#altitude values flown while migration
altFinder<-function(tm){
  alt<- lapply(tm, FUN=function(x){unlist(lapply(x, FUN = function(y){altitude<-alt_dta[which(abs(as.numeric(y) - alt_dta[,1]) == min(abs(as.numeric(y) - alt_dta[,1]))),2]}),recursive = T, use.names = T)
    })
}

alt<-altFinder(tm)

alt_idx<- seq(1,length(alt))
altcalc<- function(alt_idx){
  #alt<- unlist(lapply(alt_idx, FUN = function(x){mean(alt[x])}),recursive = T, use.names = T)
  alt<-  unlist(lapply(alt_idx, FUN = function(x){median(alt[[x]])}),recursive = T, use.names = T)
  alt <- mround(alt, 25)
  alt <- ifelse(alt == 625, 650, alt)
  alt <- ifelse(alt == 675, 650, alt)
  alt <- ifelse(alt == 725, 700, alt)
  alt <- ifelse(alt > 1000, 1000, alt)  
}

alt<- altcalc(alt_idx)
alt

#compute bird flying directions on current chain
directions <- function(x, z) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad * x[, 2L])
  sinx2 <- sin(rad * x[, 2L])
  cosz2 <- cos(rad * z[, 2L])
  sinz2 <- sin(rad * z[, 2L])
  dLon1 <- rad * (x[-n, 1L] - z[, 1L])
  dLon2 <- rad * (x[-1, 1L] - z[, 1L])
  b1 <- atan2(sin(dLon1) * cosx2[-n], cosz2 * 
                sinx2[-n] - sinz2 * cosx2[-n] * cos(dLon1))
  b2 <- atan2(sin(dLon2) * cosx2[-1L], cosz2 * 
                sinx2[-1L] - sinz2 * cosx2[-1L] * cos(dLon2))
  bird_vec_angle <- ifelse(((b1 + b2)/2) > 0, 
                           (b1 + b2)/2, ((((b1 + b2)/2) * (180/pi) + 
                                            360) * (pi/180)))
}

a<-unique(sort(alt))
str_1 <-paste(paste("u_",a, ",v_",a, sep="") ,",",collapse = "") 
db_heights<- substr(str_1, 1,nchar(str_1)-1)

#compute one time stamp per migratory bout 
tmfk<- function(time_list){
  #'Function to compute the closest available time stamp (related to a wind measurement) in the wind DB
  #'@param time_list list of unique time stamps in DB
  #'@return list of timestamp at which the wind speed is queried from the DB
#Find closest time stamp available in db for data retrieval
step_size <- time_list$tms[2] - time_list$tms[1]
start_t_idx <- time_list$tms[1]/step_size
t_stamp_idx <- round(wind_t/step_size)
list_idx <- t_stamp_idx - start_t_idx
tm_stamp <- time_list$tms[list_idx]
}
tm<-tmfk(time_list)

#create query text for times
str_1<- paste("tms=", tm," OR ", sep="",collapse = "")
tm_str<- substr(str_1,1,nchar(str_1)-3)

# Retrieve all distinct points from db
query_winds <- paste0("SELECT ",db_heights,", tms, lats, lons
                      FROM ", "\"WindTable\""," AS WT
                      WHERE ",
                      tm_str,"
                      ORDER BY tms, 
                      lats, 
                      lons
                      ;")

query_pts <- paste0("SELECT DISTINCT lats, lons
                    FROM ", "\"WindTable\"","
                    WHERE tms = 1463011200
                    ORDER BY lats, lons;")
test <- dbGetQuery(con, query_pts)
test_wind<- dbGetQuery(con, query_winds)

#bind wind positions
wind_pos_distinct <- test[, c("lons","lats")]

#parameters to extract lat and lon columns
latcol<- length(test_wind[1,])-1L
loncol<- length(test_wind[1,])
#retrieve wind positions retruned from DB to find the correct wind at a given position
latWind<- test_wind[,latcol][which(test_wind$tms == tm[1])]
lonWind<-test_wind[,loncol][which(test_wind$tms == tm[1])]

wind_positions_DB<- cbind(lonWind,latWind)


#define a index list indicating the position number
position_idx<- seq(1:(length(z0[,1])))

#threshold value to set the length of the list segments to be concatenated
th_listconcat<- length(wind_positions_DB[,1])

#generate subset of the windtable, to speed up wind data retrieval
#list of times en route in lenght of available positions
tms<- rep(tm[1], length.out = th_listconcat)
for (i in position_idx[-1]){
  tms<-append(tms, rep(tm[i],length.out=th_listconcat), after = length(tms))
}


#list of winds at required heights for each position in u direction
u_list<- test_wind[,paste("u_",alt[1], collapse="",sep="")][which(test_wind$tms == tm[1])]

for (i in position_idx[-1]){
  u_list<-append(u_list, test_wind[,paste("u_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])], after = length(u_list))
}

#list of winds at requiredheights for each position in v direction
v_list<- test_wind[,paste("v_",alt[1], collapse="",sep="")][which(test_wind$tms == tm[1])]

for (i in position_idx[-1]){
  v_list<-append(v_list, test_wind[,paste("u_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])], after = length(v_list))
}

#create a geoindex
idx<- rep(1:th_listconcat)
geoIndex<- lapply(idx,FUN=function(x){
  which(test_wind[,latcol] == test[x,1] & test_wind[,loncol] == test[x,2])
})

#bind all lists to one df
prep_data<- cbind(tms,u_list,v_list)


wind_model <- function(x,z, bird_vec_angle, dt) {
  #' Function which computes the probability of each position along the migration track. THe probability describes, how likely a bird flys to this point whilst migrating
  #' from the breeding site to the wintering site and back, given the data recorded while bird migration.
  #' @param x chain of x positions along the track. These points are positions estimated by light location
  #' @param z chain of z positions intermediate between the x positions computed applying the trackMidPts function in SGAT
  #' @param bird_vec_angle list of direction angles in radian, indicating  the bird flight direction
  #' @param dt flt time given by bird activity which was categorized as migration flight
  #' @return prob_list list of  probability values indicating how probable a point is flown to given the data available
  #' 
  #Compute speed necessary to fly between the proposed points in the given flt time
  spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
  
  #create subset for distance matrix
  #all indexes of the wind positions + and - 3 degree around the bird position are returned. These positions are then used to compute
  #wind speed and direction 
  #define search radius around bird position in degrees
  sRad<- 30
  pos_subset<- lapply(position_idx, FUN = function(x){which(wind_positions_DB[,1] <= z[x,] + sRad & wind_positions_DB[,2] <= z[x,] + sRad & wind_positions_DB[,1] >= z[x,]-sRad & wind_positions_DB[,2] >= z[x,]-sRad)})
  
  # Compute distance from each Point in x to each wind measurement position  found within +/- 3 degrees in wind_position_DB data table.
  dist <- lapply(position_idx, FUN=function(x){
    #start points of distance matrix
    from<-wind_positions_DB[pos_subset[[x]],]
    to<-z[x,]
    suppressMessages(geodist(from,to,measure='cheap'))
  })
  
  k_nn_data <- lapply(position_idx,function(x){
    k <- 9  
    # Find the value of the k th smallest distance 
    # (partial = k implies  that the sorting stops once the kth nearest neighbor has been found)
    k_th <- sort(dist[[x]], partial = k)[k]
    # Find the indices of all elements smaller than k_nn (not ordered!!)
    k_nn_idx <- which(dist[[x]] <= k_th)
    #Filter the nine closest distances
    k_nn_dist<- dist[[x]][k_nn_idx]
    return_data<- cbind(k_nn_idx, k_nn_dist)
    return(return_data)
  } )
  
  #filter the wind speeds at the correct height at each position in x and y direction
  wind_speeed_data<-lapply(position_idx,FUN= function(x){
    nn_positions<- pos_subset[[x]][k_nn_data[[x]][,1]]
    wind_speed_x<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],2]})
    wind_speed_y<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],3]})
    #bind the x and y column to one data frame for simple return and preocessing
    wind_data<-cbind(unlist(wind_speed_x), unlist(wind_speed_y))
    return(wind_data)
  })
  
  #Weights computed for weighted mean. The closest wind measurement gets the highest weight, while the furthest wind measurement gets the smalest weight
  dist_min <- lapply(position_idx, FUN=function(x){
    min(k_nn_data[[x]][,2])})
  
  dist_max <- lapply(position_idx, FUN=function(x){
    max(k_nn_data[[x]][,2])})
  
  weights<- lapply(position_idx, FUN=function(x){(k_nn_data[[x]][,2]- dist_max[[x]]) / (dist_min[[x]] - dist_max[[x]])})
  
  
  weig_wind_speed<- lapply(position_idx,FUN=function(x){
    #compute the weighted speeds on each position along the chain
    weig_speeds<- wind_speeed_data[[x]] * weights[[x]]
    return(weig_speeds)
    
  })
  
  #sum the weights at each position along the track
  sum_weights<- lapply(weights,sum)
  
  #Compute the weighted mean of each wind in x and y direction
  wind_speeds_xy<- lapply(position_idx,FUN=function(x){
    #sum the weighted speeds and devide by the sum of weights -> weighted mean                  
    swspeeds<-apply(weig_wind_speed[[x]], 2,sum) / sum_weights[[x]]})
  
  #parameter to transform m/s to km/h
  ms2kmh<- 3.6
  
  #compute wind speed and wind direction at each position
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2)/2 * ms2kmh}),recursive = T, use.names = T)
  
  wind_directions<- unlist(lapply(position_idx,FUN = function(x){atan2(wind_speeds_xy[[x]][[2]],wind_speeds_xy[[x]][[1]])}), recursive = T, use.names = T)
  
  
  ##' compute the wind speed sigma
  ##' sigma is used to compute the alpha and beta parameter to characterize the gamma distribution --> log probability
  
  #define percentage of data at mu and at two sigma from mu
  nmu<-50
  twosgm<- 2.5
  #Parameter deciding if a distribution is onesided or two sided
  onesided<- 2
  #wind speed sigma in m/s set to 3.7 m/s (Alessandrini et. al. 2012)
  #if Weibul distributed insert following line:
  #((wind_speeds/50)*2.5)/2
  sigma_wind_spd <- 3.7
  
  #compute the angel between the bird flying vector and the wind direction
  alpha <- ifelse(bird_vec_angle > wind_directions, 
                  bird_vec_angle - wind_directions, wind_directions - 
                    bird_vec_angle)
  alpha <- ifelse(alpha <= (pi/2), alpha, ifelse(alpha<=pi, pi - alpha, alpha-pi))
  
  #retrieve the direction of the vector describing the wind compenstaion a bird has to fly in order to fly to its desired destination
  dir_comp <- ifelse(bird_vec_angle > wind_directions, bird_vec_angle - (pi/2), bird_vec_angle - ((3/2) * pi))
  
  dir_comp <- ifelse(dir_comp > 0, dir_comp, dir_comp + 2 * pi)
  
  #compute magnitude of wind compensation vector
  mag_comp <- cos(pi - (pi/2) - alpha[3]) * wind_speeds
  
  #After computing the compensation vector, the wind add parallel to the bird flying vector is computed
  wind_add <- ifelse(alpha<=pi | alpha>= ((2/3)*pi), sqrt(wind_speeds^2 - mag_comp^2),-sqrt(wind_speeds^2 - mag_comp^2))
  
  #Sum of wind partition parallel to bird flying vector and bird speed
  new_b_speed <- (cos(alpha)*wind_speeds) + sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2)
  
  #sigma of bird speed according to Bruderer er al, Bolch & Brudrer and Liechti & Bruderer
  bird_speed_sigma <- 3.5 
  
  #standard deviation of brid ground speed (Bird airspeed + wind speed)
  std <- bird_speed_sigma + sigma_wind_spd
  
  #Beta applied to describe gamma distribution
  beta_dist <- new_b_speed/(std^2)
  
  
  #Alpha applied to describe gamma distribution
  alpha_dist <- new_b_speed * beta_dist
  
  #compute probabilities of each point along the migration route
  prob_list <- dgamma(spd, alpha_dist, beta_dist, 
                      log = T)
  
  return(prob_list)
}

#change the SpeedGammaModel in SGAT to make sure our Movement model defined here is called when probabilities are computed
#remember to copy paste the code from SpeedGammaModel.R and paste it between the {}
trace(speedGammaModel,edit = T)

#After tracing the function remember to rerun the model (groupedThresholdModel)
model <- groupedThresholdModel(twl$Twilight,
                               twl$Rise,
                               group = twl$group, #This is the group vector for each time the bird was at a point
                               twilight.model = "ModifiedGamma",
                               alpha = alpha,
                               beta =  beta,
                               x0 = x0, # meadian point for each greoup (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zenith0+.5, # Here I use the zenith0 value from in-habitat calibration. Using zenith or zenith_sd is not recommended here and will give you wrong results as they are medina SEAs
                               logp.x = logp, # land sea mask
                               dt = flt_t, # This is the new part here - movment is allowed only for the number of hours that PAM analyses said the bird was flying
                               fixedx = fixedx)


######## Route inference section ##################
#Start the BI and infer the routes!!!

t1<- Sys.time()
t1

# Fit the model
fit<- estelleMetropolis(model, x.proposal, z.proposal, iters = 1000, thin = 20)
# plot outcome

par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
plot(newmap, xlim=xlim, ylim=ylim)
xm <- locationMean(fit$x)
lines(xm, col = "blue")
points(xm, col = ifelse(stationary,rainbow(max(sitenum)),1), pch = 16, cex = ifelse(stationary, 2, 0.5))
points(lon.calib, lat.calib, pch = 13, col = 2,cex=2.5,lwd=3)
box()
t2<-Sys.time()

t<-t2-t1
t

# Tune the model
# use output from last run
x0 <- chainLast(fit$x)[[1]]
z0 <- chainLast(fit$z)[[1]]




# NB! If this next model doesn't run, there is something wrong with either calibration or stationary-movement period designation. 
# A quick way to see if calibration may be wrong is to increase the zenith0 value and see if that helps the model to run - see the "+ 0.5" I did in this case
# if that does not help then some of your stationary periods likely include movement and you have to go back and redefine them 

model <- groupedThresholdModel(twl$Twilight,
                               twl$Rise,
                               group = twl$group, #This is the group vector for each time the bird was at a point
                               twilight.model = "Gamma", # this bit gere is what differs from the previously set model. "ModifiedGamma" will work always, but "Gamma" will fail when designation between stationary-movement phases are likely worong
                               alpha = alpha,
                               beta =  beta,
                               x0 = x0, # meadian point for each greoup (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zenith0+.5, # I often see that zenith0 for wintering period is typically a bit higher that zenith0 for breeding site, thus I manually add anywhere between 0.5 to 2 zenith0 - your model will also fail if the original zenith0 value is too low for the wintering period - just test adding a bit of the model fails
                               logp.x = logp, # land sea mask
                               dt = twl$flight_duration[twl$flight_duration>0],
                               fixedx = fixedx)

t1<-Sys.time()
t1

for (k in 1:5) {
  x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
  z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)
  fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x)[[1]],
                           z0 = chainLast(fit$z)[[1]], iters = 300, thin = 20)
}
t2<- Sys.time()
t<-t2-t1
t


## Check if MCMC chains mix - you generally want a messy plot with no obvious white areas somewhere in the middle
opar<- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!fixedx, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!fixedx, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)
# plot outcome
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
plot(newmap, xlim=xlim, ylim=ylim)
xm <- locationMean(fit$x)
lines(xm, col = "blue")
points(xm, col = ifelse(stationary,rainbow(max(sitenum)),1), pch = 16, cex = ifelse(stationary, 2, 0.5))


points(lon.calib, lat.calib, pch = 13, col = 2,cex=2.5,lwd=3)
box()

# Final run
x.proposal <- mvnorm(chainCov(fit$x), s = 0.3)
z.proposal <- mvnorm(chainCov(fit$z), s = 0.3)

t1<-Sys.time()
t1

fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 2000, thin = 20, chain = 1)
t2<-Sys.time()
t<-t2-t1
t

# Summarize the results----
sm <- locationSummary(fit$x, time=fit$model$time)


# Plot the results----
# pdf(paste0("PAM analyses/Results/Tracking/", substring(ID,1,5), "_SGAT_GroupSummary1.pdf"))

png(file=Storage_path_route)
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) )
colours <- c("black",colorRampPalette(c("blue","yellow","red"))(max(sitenum)))

r <- raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,# make empty raster of the extent
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = proj4string(wrld_simpl))

s <- slices(type = "intermediate", mcmc = fit, grid = r)
# sk <- slice(s, sliceIndices(s))
sk <- log(SGAT::slice(s, sliceIndices(s))+1)

plot(wrld_simpl, xlim=xlim, ylim=ylim)
plot(sk, useRaster = F,col = c("transparent", rev(viridis::viridis(50,alpha=.8))),add=T)

with(sm[sitenum>0,], arrows(`Lon.50%`, `Lat.50%`+`Lat.sd`, `Lon.50%`, `Lat.50%`-`Lat.sd`, length = 0, lwd = 2.5, col = "firebrick"))
with(sm[sitenum>0,], arrows(`Lon.50%`+`Lon.sd`, `Lat.50%`, `Lon.50%`-`Lon.sd`, `Lat.50%`, length = 0, lwd = 2.5, col = "firebrick"))

sz <- locationSummary(fit$z, time=fit$model$time)
track<-matrix(ncol=2)
for (i in 1:length(sitenum)){
  if (sitenum[i]>0){track<-rbind(track,as.numeric(sm[i,c("Lon.50%","Lat.50%")]))}
  else{track<-rbind(track,as.numeric(sz[i,c("Lon.50%","Lat.50%")]))}
}
track<-track[2:nrow(track),]
lines(track[,1], track[,2], col = "darkorchid4", lwd = 2)

points(sm[,"Lon.50%"][sitenum>0], sm[,"Lat.50%"][sitenum>0],pch=21,
       bg=rainbow(max(sitenum)), cex = 2.5, col = "firebrick", lwd = 2.5)

points(sm[,"Lon.50%"][sitenum>0], sm[,"Lat.50%"][sitenum>0],pch=97:(97+max(sitenum)-2),
       cex = 1)



dev.off()


track2<-matrix(ncol=2)
for (i in 1:length(sitenum)){
  if (sitenum[i]>0){track2<-rbind(track2,matrix(rep(as.numeric(sm[i,c("Lon.50%","Lat.50%")]),length(gr[gr==unique(gr)[i]])),ncol=2,byrow = T))}
  else{track2<-rbind(track2,as.numeric(sz[i,c("Lon.50%","Lat.50%")]))}
}
track2<-track2[2:nrow(track2),]

#Look at latitude and longitude over time
temp<-data.frame()
for (n in 1:length(twl$group)){
  t<-rbind(sm[twl$group[n],])
  temp<-rbind(temp,t)
}
temp<-cbind(twl[,1:3],twl$group,temp)
colnames(temp)[4] <- "movement_vector"
temp$`Lon.50%`<-track2[,1]
temp$`Lat.50%`<-track2[,2]



png(file=LatAndLonPlot_Path)
par(mfrow=c(2,1),mar = c(4, 4, 2, 2))
plot(temp$`Lat.50%`~twl$Twilight, pch=ifelse(sites==0,1,16), ylim=ylim,
     col=ifelse(sites==0,1,rainbow(max(sitenum))[replace(sites,sites==0,NA)]),
     type="o",ylab="latitude",xlab=NA)
plot(temp$`Lon.50%`~twl$Twilight, pch=ifelse(sites==0,1,16),ylim=xlim,
     col=ifelse(sites==0,1,rainbow(max(sitenum))[replace(sites,sites==0,NA)]),
     type="o",ylab="longitude",xlab=NA)

dev.off()

#get harmonic mean of coordinates of all iterations
idx_hm<- seq(1,length(x0[,1]))
lat_hm<- unlist(lapply(idx_hm, FUN=function(x){harmonic.mean(fit$x[[1]][x,2,])}))
lon_hm<- unlist(lapply(idx_hm, FUN=function(x){harmonic.mean(fit$x[[1]][x,1,])}))

x_hm<- cbind(lon_hm,lat_hm)

#compute middle pts amd flit directionof harmonic mean route
z_hm<- trackMidpts(x_hm)

bird_vec_angle<- directions(x_hm,z_hm)

#place harmonic mean of lat and lon
sm<-cbind(sm, lat_hm)
sm<-cbind(sm,lon_hm)

## Visualization of results and data
x0 <- track
z0 <- trackMidpts(x0)

spd <- pmax.int(trackDist2(x0, z0), 1e-06)/flt_t
bird_vec_angle<- directions(x0,z0)

prob_list_log<- wind_model(x0,z0, bird_vec_angle, flt_t)
prob_list_log[length(prob_list_log)+1]<- prob_list_log[1]

sm<- cbind(sm, prob_list_log)

wind_model <- function(x,z, bird_vec_angle, dt) {
  #' Function which computes the probability of each position along the migration track. THe probability describes, how likely a bird flys to this point whilst migrating
  #' from the breeding site to the wintering site and back, given the data recorded while bird migration.
  #' @param x chain of x positions along the track. These points are positions estimated by light location
  #' @param z chain of z positions intermediate between the x positions computed applying the trackMidPts function in SGAT
  #' @param bird_vec_angle list of direction angles in radian, indicating  the bird flight direction
  #' @param dt flt time given by bird activity which was categorized as migration flight
  #' @return prob_list list of  probability values indicating how probable a point is flown to given the data available
  #' 
  #Compute speed necessary to fly between the proposed points in the given flt time
  spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
  
  #create subset for distance matrix
  #all indexes of the wind positions + and - 3 degree around the bird position are returned. These positions are then used to compute
  #wind speed and direction 
  #define search radius around bird position in degrees
  sRad<- 25
  pos_subset<- lapply(position_idx, FUN = function(x){which(wind_positions_DB[,1] <= z[x,] + sRad & wind_positions_DB[,2] <= z[x,] + sRad & wind_positions_DB[,1] >= z[x,]-sRad & wind_positions_DB[,2] >= z[x,]-sRad)})
  
  # Compute distance from each Point in x to each wind measurement position  found within +/- 3 degrees in wind_position_DB data table.
  dist <- lapply(position_idx, FUN=function(x){
    #start points of distance matrix
    from<-wind_positions_DB[pos_subset[[x]],]
    to<-z[x,]
    suppressMessages(geodist(from,to,measure='cheap'))
  })
  
  k_nn_data <- lapply(position_idx,function(x){
    k <- 9  
    # Find the value of the k th smallest distance 
    # (partial = k implies  that the sorting stops once the kth nearest neighbor has been found)
    k_th <- sort(dist[[x]], partial = k)[k]
    # Find the indices of all elements smaller than k_nn (not ordered!!)
    k_nn_idx <- which(dist[[x]] <= k_th)
    #Filter the nine closest distances
    k_nn_dist<- dist[[x]][k_nn_idx]
    return_data<- cbind(k_nn_idx, k_nn_dist)
    return(return_data)
  } )
  
  #filter the wind speeds at the correct height at each position in x and y direction
  wind_speeed_data<-lapply(position_idx,FUN= function(x){
    nn_positions<- pos_subset[[x]][k_nn_data[[x]][,1]]
    wind_speed_x<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],2]})
    wind_speed_y<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],3]})
    #bind the x and y column to one data frame for simple return and preocessing
    wind_data<-cbind(unlist(wind_speed_x), unlist(wind_speed_y))
    return(wind_data)
  })
  
  #Weights computed for weighted mean. The closest wind measurement gets the highest weight, while the furthest wind measurement gets the smalest weight
  dist_min <- lapply(position_idx, FUN=function(x){
    min(k_nn_data[[x]][,2])})
  
  dist_max <- lapply(position_idx, FUN=function(x){
    max(k_nn_data[[x]][,2])})
  
  weights<- lapply(position_idx, FUN=function(x){(k_nn_data[[x]][,2]- dist_max[[x]]) / (dist_min[[x]] - dist_max[[x]])})
  
  
  weig_wind_speed<- lapply(position_idx,FUN=function(x){
    #compute the weighted speeds on each position along the chain
    weig_speeds<- wind_speeed_data[[x]] * weights[[x]]
    return(weig_speeds)
    
  })
  
  #sum the weights at each position along the track
  sum_weights<- lapply(weights,sum)
  
  #Compute the weighted mean of each wind in x and y direction
  wind_speeds_xy<- lapply(position_idx,FUN=function(x){
    #sum the weighted speeds and devide by the sum of weights -> weighted mean                  
    swspeeds<-apply(weig_wind_speed[[x]], 2,sum) / sum_weights[[x]]})
  
  #parameter to transform m/s to km/h
  ms2kmh<- 3.6
  
  #compute wind speed and wind direction at each position
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2)/2 * ms2kmh}),recursive = T, use.names = T)
  
  wind_directions<- unlist(lapply(position_idx,FUN = function(x){atan2(wind_speeds_xy[[x]][[2]],wind_speeds_xy[[x]][[1]])}), recursive = T, use.names = T)
  
  
  ##' compute the wind speed sigma
  ##' sigma is used to compute the alpha and beta parameter to characterize the gamma distribution --> log probability
  
  #Weibul distributed wind (is approximated with a normal distribution which is defined between 0 and 2*mu) 
  nmu<-50
  twosgm<- 2.5
  #wind speed sigma in m/s set to 3.7 m/s (Alessandrini et. al. 2012)
  #if Weibul distributed insert following line:
  #((wind_speeds/50)*2.5)/2
  
  sigma_wind_spd <- 3.7
  
  #compute the angel between the bird flying vector and the wind direction
  alpha <- ifelse(bird_vec_angle > wind_directions, 
                  bird_vec_angle - wind_directions, wind_directions - 
                    bird_vec_angle)
  alpha <- ifelse(alpha <= (pi/2), alpha, ifelse(alpha<=pi, pi - alpha, alpha-pi))
  
  #retrieve the direction of the vector describing the wind compenstaion a bird has to fly in order to fly to its desired destination
  dir_comp <- ifelse(bird_vec_angle > wind_directions, bird_vec_angle - (pi/2), bird_vec_angle - ((3/2) * pi))
  
  dir_comp <- ifelse(dir_comp > 0, dir_comp, dir_comp + 2 * pi)
  
  #compute magnitude of wind compensation vector
  #mag_comp <- cos(pi - (pi/2) - alpha[3]) * wind_speeds
  mag_comp <- sin(alpha)*wind_speeds
  
  #After computing the compensation vector, the ground speed vector is computed
  new_b_speed <- (cos(alpha)*wind_speeds) + sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2)
  
  #sigma of bird speed according to Bruderer er al, Bolch & Brudrer and Liechti & Bruderer
  bird_speed_sigma <- 3.5 
  
  #standard deviation of brid ground speed (Bird airspeed + wind speed)
  std <- bird_speed_sigma + sigma_wind_spd
  
  #Beta applied to describe gamma distribution
  beta_dist <- new_b_speed/(std^2)
  
  
  #Alpha applied to describe gamma distribution
  alpha_dist <- new_b_speed * beta_dist
  
  #compute probabilities of each point along the migration route
  prob_list <- dgamma(spd, alpha_dist, beta_dist, 
                      log = F)
  
  return(prob_list)
}

prob_list<- wind_model(x0,z0, bird_vec_angle, flt_t)
prob_list[length(prob_list)+1]<- prob_list[1]

sm<- cbind(sm, prob_list)


write.csv(sm, Summaryfile_Path, row.names = F)


plot(prob_list,type="l", main="Probability of trajectories along the migration path", ylab="Probability", xlab="Position index")

wind_model <- function(x,z, bird_vec_angle, dt) {
  #' Function which computes the probability of each position along the migration track. THe probability describes, how likely a bird flys to this point whilst migrating
  #' from the breeding site to the wintering site and back, given the data recorded while bird migration.
  #' @param x chain of x positions along the track. These points are positions estimated by light location
  #' @param z chain of z positions intermediate between the x positions computed applying the trackMidPts function in SGAT
  #' @param bird_vec_angle list of direction angles in radian, indicating  the bird flight direction
  #' @param dt flt time given by bird activity which was categorized as migration flight
  #' @return prob_list list of  probability values indicating how probable a point is flown to given the data available
  #' 
  #Compute speed necessary to fly between the proposed points in the given flt time
  spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
  
  #create subset for distance matrix
  #all indexes of the wind positions + and - 3 degree around the bird position are returned. These positions are then used to compute
  #wind speed and direction 
  #define search radius around bird position in degrees
  sRad<- 25
  pos_subset<- lapply(position_idx, FUN = function(x){which(wind_positions_DB[,1] <= z[x,] + sRad & wind_positions_DB[,2] <= z[x,] + sRad & wind_positions_DB[,1] >= z[x,]-sRad & wind_positions_DB[,2] >= z[x,]-sRad)})
  
  # Compute distance from each Point in x to each wind measurement position  found within +/- 3 degrees in wind_position_DB data table.
  dist <- lapply(position_idx, FUN=function(x){
    #start points of distance matrix
    from<-wind_positions_DB[pos_subset[[x]],]
    to<-z[x,]
    suppressMessages(geodist(from,to,measure='cheap'))
  })
  
  k_nn_data <- lapply(position_idx,function(x){
    k <- 9  
    # Find the value of the k th smallest distance 
    # (partial = k implies  that the sorting stops once the kth nearest neighbor has been found)
    k_th <- sort(dist[[x]], partial = k)[k]
    # Find the indices of all elements smaller than k_nn (not ordered!!)
    k_nn_idx <- which(dist[[x]] <= k_th)
    #Filter the nine closest distances
    k_nn_dist<- dist[[x]][k_nn_idx]
    return_data<- cbind(k_nn_idx, k_nn_dist)
    return(return_data)
  } )
  
  #filter the wind speeds at the correct height at each position in x and y direction
  wind_speeed_data<-lapply(position_idx,FUN= function(x){
    nn_positions<- pos_subset[[x]][k_nn_data[[x]][,1]]
    wind_speed_x<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],2]})
    wind_speed_y<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],3]})
    #bind the x and y column to one data frame for simple return and preocessing
    wind_data<-cbind(unlist(wind_speed_x), unlist(wind_speed_y))
    return(wind_data)
  })
  
  #Weights computed for weighted mean. The closest wind measurement gets the highest weight, while the furthest wind measurement gets the smalest weight
  dist_min <- lapply(position_idx, FUN=function(x){
    min(k_nn_data[[x]][,2])})
  
  dist_max <- lapply(position_idx, FUN=function(x){
    max(k_nn_data[[x]][,2])})
  
  weights<- lapply(position_idx, FUN=function(x){(k_nn_data[[x]][,2]- dist_max[[x]]) / (dist_min[[x]] - dist_max[[x]])})
  
  
  weig_wind_speed<- lapply(position_idx,FUN=function(x){
    #compute the weighted speeds on each position along the chain
    weig_speeds<- wind_speeed_data[[x]] * weights[[x]]
    return(weig_speeds)
    
  })
  
  #sum the weights at each position along the track
  sum_weights<- lapply(weights,sum)
  
  #Compute the weighted mean of each wind in x and y direction
  wind_speeds_xy<- lapply(position_idx,FUN=function(x){
    #sum the weighted speeds and devide by the sum of weights -> weighted mean                  
    swspeeds<-apply(weig_wind_speed[[x]], 2,sum) / sum_weights[[x]]})
  
  #compute wind speed and wind direction at each position
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2)/2 }),recursive = T, use.names = T)
  
}

wind_speeds<- wind_model(x0,z0,bird_vec_angle,flt_t)
hist(wind_speeds, main="Histogram of wind speeds on the migrationroute", xlab="Wind speeds [km/h]")
png(file=WindSpeed_Path)
barplot(wind_speeds,names.arg = c(1:length(wind_speeds)),
        main="Wind speeds encountered en route",
        xlab="Position index",
        ylab="Winds speed [m/s]")
dev.off()
mw<- mean(wind_speeds)
sdw<- sd(wind_speeds)
minw<-min(wind_speeds)
maxw<-max( wind_speeds)

wind_data<- cbind(mw, sdw, minw, maxw)

wind_speed_data<- cbind(c(1:length(wind_speeds)), wind_speeds)
windspeedfile<- write.csv(wind_speed_data, WindSpeedCSV_Path, row.names = F)

wind_model <- function(x,z, bird_vec_angle, dt) {
  #' Function which computes the probability of each position along the migration track. THe probability describes, how likely a bird flys to this point whilst migrating
  #' from the breeding site to the wintering site and back, given the data recorded while bird migration.
  #' @param x chain of x positions along the track. These points are positions estimated by light location
  #' @param z chain of z positions intermediate between the x positions computed applying the trackMidPts function in SGAT
  #' @param bird_vec_angle list of direction angles in radian, indicating  the bird flight direction
  #' @param dt flt time given by bird activity which was categorized as migration flight
  #' @return prob_list list of  probability values indicating how probable a point is flown to given the data available
  #' 
  #Compute speed necessary to fly between the proposed points in the given flt time
  spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
  
  #create subset for distance matrix
  #all indexes of the wind positions + and - 3 degree around the bird position are returned. These positions are then used to compute
  #wind speed and direction 
  #define search radius around bird position in degrees
  sRad<- 25
  pos_subset<- lapply(position_idx, FUN = function(x){which(wind_positions_DB[,1] <= z[x,] + sRad & wind_positions_DB[,2] <= z[x,] + sRad & wind_positions_DB[,1] >= z[x,]-sRad & wind_positions_DB[,2] >= z[x,]-sRad)})
  
  # Compute distance from each Point in x to each wind measurement position  found within +/- 3 degrees in wind_position_DB data table.
  dist <- lapply(position_idx, FUN=function(x){
    #start points of distance matrix
    from<-wind_positions_DB[pos_subset[[x]],]
    to<-z[x,]
    suppressMessages(geodist(from,to,measure='cheap'))
  })
  
  k_nn_data <- lapply(position_idx,function(x){
    k <- 9  
    # Find the value of the k th smallest distance 
    # (partial = k implies  that the sorting stops once the kth nearest neighbor has been found)
    k_th <- sort(dist[[x]], partial = k)[k]
    # Find the indices of all elements smaller than k_nn (not ordered!!)
    k_nn_idx <- which(dist[[x]] <= k_th)
    #Filter the nine closest distances
    k_nn_dist<- dist[[x]][k_nn_idx]
    return_data<- cbind(k_nn_idx, k_nn_dist)
    return(return_data)
  } )
  
  #filter the wind speeds at the correct height at each position in x and y direction
  wind_speeed_data<-lapply(position_idx,FUN= function(x){
    nn_positions<- pos_subset[[x]][k_nn_data[[x]][,1]]
    wind_speed_x<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],2]})
    wind_speed_y<- lapply(nn_positions,FUN=function(p){test_wind[geoIndex[[p]][x],3]})
    #bind the x and y column to one data frame for simple return and preocessing
    wind_data<-cbind(unlist(wind_speed_x), unlist(wind_speed_y))
    return(wind_data)
  })
  
  #Weights computed for weighted mean. The closest wind measurement gets the highest weight, while the furthest wind measurement gets the smalest weight
  dist_min <- lapply(position_idx, FUN=function(x){
    min(k_nn_data[[x]][,2])})
  
  dist_max <- lapply(position_idx, FUN=function(x){
    max(k_nn_data[[x]][,2])})
  
  weights<- lapply(position_idx, FUN=function(x){(k_nn_data[[x]][,2]- dist_max[[x]]) / (dist_min[[x]] - dist_max[[x]])})
  
  
  weig_wind_speed<- lapply(position_idx,FUN=function(x){
    #compute the weighted speeds on each position along the chain
    weig_speeds<- wind_speeed_data[[x]] * weights[[x]]
    return(weig_speeds)
    
  })
  
  #sum the weights at each position along the track
  sum_weights<- lapply(weights,sum)
  
  #Compute the weighted mean of each wind in x and y direction
  wind_speeds_xy<- lapply(position_idx,FUN=function(x){
    #sum the weighted speeds and devide by the sum of weights -> weighted mean                  
    swspeeds<-apply(weig_wind_speed[[x]], 2,sum) / sum_weights[[x]]})
  
  #parameter to transform m/s to km/h
  ms2kmh<- 3.6
  
  #compute wind speed and wind direction at each position
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2)}),recursive = T, use.names = T)
  
  wind_directions<- unlist(lapply(position_idx,FUN = function(x){atan2(wind_speeds_xy[[x]][[2]],wind_speeds_xy[[x]][[1]])}), recursive = T, use.names = T)
  
  
  ##' compute the wind speed sigma
  ##' sigma is used to compute the alpha and beta parameter to characterize the gamma distribution --> log probability
  
  #define percentage of data at mu and at two sigma from mu
  nmu<-50
  twosgm<- 30
  #Parameter deciding if a distribution is onesided or two sided
  onesided<- 2
  sigma_wind_spd <- (wind_speeds - ((wind_speeds/nmu) * twosgm)/onesided)
  
  #compute the angel between the bird flying vector and the wind direction
  alpha <- ifelse(bird_vec_angle > wind_directions, 
                  bird_vec_angle - wind_directions, wind_directions - 
                    bird_vec_angle)
  alpha <- ifelse(alpha <= (pi/2), alpha, ifelse(alpha<=pi, pi - alpha, alpha-pi))
  
  #retrieve the direction of the vector describing the wind compenstaion a bird has to fly in order to fly to its desired destination
  dir_comp <- ifelse(bird_vec_angle > wind_directions, bird_vec_angle - (pi/2), bird_vec_angle - ((3/2) * pi))
  
  dir_comp <- ifelse(dir_comp > 0, dir_comp, dir_comp + 2 * pi)
  
  #compute magnitude of wind compensation vector
  mag_comp <- cos(pi - (pi/2) - alpha[3]) * wind_speeds
  
  #After computing the compensation vector, the wind add parallel to the bird flying vector is computed
  wind_add <- ifelse(alpha<=pi, sqrt(wind_speeds^2 - mag_comp^2),-sqrt(wind_speeds^2 - mag_comp^2))
  
  #Sum of wind partition parallel to bird flying vector and bird speed
}

wind_add<-wind_model(x0,z0,bird_vec_angle,flt_t)

hist(wind_add, main="Histogram of wind partition parallel to bird flight vector", xlab="Speed [km/h]")
png(file=WindSupport_Path)
barplot(wind_add,names.arg = c(1:length(wind_add)),
        main="Wind speeds parallel to flight vector encountered en route",
        xlab="Position index",
        ylab="Winds speed [m/s]")
dev.off()
mwa<- mean(wind_add)
sdwa<- sd(wind_add)
minwa<- min(wind_add)
maxwa<- max(wind_add)

wind_suport_data<- cbind(c(1:length(wind_speeds)), wind_add)
windsuportfile<- write.csv(wind_suport_data, WindSuportCSV_Path, row.names = F)

wind_add_data<- cbind(mwa, sdwa, minwa, maxwa)
wind_data_csv<-rbind(wind_data, wind_add_data)
write.csv(wind_data_csv, WindData_Path, row.names = F)


# Save Results - important to save the final outcome = fit !!!
# save(fit,
#      file=paste0("PAM analyses/Results/Tracking/", substring(ID,1,5),"_SGAT_Groupfit.RData"),
#      compress = T)
# 
# # Save a simplified CSV file with main results
# write.csv(temp,
#           paste0("PAM analyses/Results/Tracking/", substring(ID,1,5), "_SGAT_GroupSummary.csv"),
#           row.names = F)
# 
# 


# END


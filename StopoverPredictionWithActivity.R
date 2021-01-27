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
folderName<- "16BG_20170612"
birdName<- "16BG"
Data <- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/data/",folderName)
# list all light files
ID.list<-list.files(paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/data/", folderName),pattern=".glf",recursive = T) # Change to gle or glf - the one you use for analyses
ID <- ID.list[1]
ID

glf_file<- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/data/",folderName, "/", folderName, ".glf")
twl_CSVfile<- paste0("~/MSC_Thesis/R_stuff/Msc_Thesis/data/",folderName,"/", folderName,"_twl.csv")
PAM_filepath<- paste0("~/MSC_Thesis/R_stuff/Msc_Thesis/data/",folderName,"/")
Activity_path<- "P:/home/Documents/MSC_Thesis/Route_Inferrence/Activity/"



#parameter for PAM analysis

# Start date of observation period --> used to crop pam data
# make sure the cropping period is in the correct date format
start = as.POSIXct("2015-07-15","%Y-%m-%d", tz="UTC")

# End date of observation period --> used to crop pam data
end = as.POSIXct("2016-06-30","%Y-%m-%d", tz="UTC") 

# Parameter for Database connection
user<- "mwerfeli"
password<- "M*4XgiFWLU4E+62r9GtM"
dbname<- "sc24"
host<- "sc24.geo.uzh.ch"

#Parameter for result storage
Route_name<- 'NW'
Storage_path_route<-paste0('P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/',birdName,'/VoWa/Route_', Route_name,'.png')
LatAndLonPlot_Path<- paste0('P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/',birdName,'/VoWa/LatLon_', Route_name,'.png')
Summaryfile_Path<- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/",birdName,"/VoWa/",birdName,"_sumary_", Route_name,".csv")


WindSpeed_Path<- paste0('P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/',birdName,'/VoWa/WindSpeed_', Route_name,'.png')
WindSupport_Path<- paste0('P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/',birdName,'/VoWa/WindSupport_', Route_name,'.png')


WindData_Path<- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/",birdName,"/VoWa/",birdName, "_",Route_name,".csv")
WindSpeedCSV_Path<- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/",birdName,"/VoWa/",birdName, "_",Route_name,"_WindSpeeds.csv") 
WindSuportCSV_Path<- paste0("P:/home/Documents/MSC_Thesis/Route_Inferrence/Results/",birdName,"/VoWa/",birdName, "_",Route_name,"_WindSuports.csv")

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

start_twl<- as.character(twl$Twilight[1])
end_twl<- as.character(tail(twl$Twilight,1))


# SGAT Group model----

# 1. Calibration----
lightImage( tagdata = raw,
            offset = offset,     
            zlim = c(0, max(raw$Light)))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset, lwd = 2, col = "orange")


tm.calib <- as.POSIXct(c("2016-07-20 00:00", "2016-08-12 00:00", "2017-04-15 00:00", "2017-04-28 00:00"), tz = "UTC",format="%Y-%m-%d %H:%M") # Selecting calibration period(s) 
abline(v = tm.calib, lwd = 2, lty = 2, col = "red") # chech if them make sense on the light graph

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2] | Twilight>=tm.calib[3] & Twilight<=tm.calib[4])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise, lon.calib, lat.calib, method = "gamma")

zenith  <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Alternatively use Hill-Ekstrom calibration from the wintering sites
startDate <-tm.calib[2]

endDate   <-tm.calib[3]


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
ID.list = list.files(paste0("~/MSC_Thesis/R_stuff/Msc_Thesis/data/",folderName,"/"),include.dirs=T)
ID.list
#Set idx at position where GLF file is found in the id list --> you want to read the glf file
ID = ID.list[3]
ID

pathname = paste0("~/MSC_Thesis/R_stuff/Msc_Thesis/data/",folderName,"/")
#Make sure all file names are correct
measurements = c(".pressure", 
                 ".glf",
                 ".acceleration", 
                 ".temperature")

PAM_data = importPAM(pathname, measurements)

# Cropping the data

# make sure the cropping period is in the correct date format
start = tm.calib[1]
# end = max(PAM_data$light$date)
end = tm.calib[4]
PAM_data= cutPAM(PAM_data,start,end)

length(PAM_data$pressure$date)
altitude = altitudeCALC(P = PAM_data$pressure$obs)
plot(PAM_data$pressure$date[2:80000], PAM_data$pressure$obs[2:80000], type="o",pch=16, xlab="Date", ylab="Air Pressure (hPa)", main = "Air Pressure en route")
alt_dta<- cbind(PAM_data$pressure$date, PAM_data$pressure$obs)

tz<- "UTC"
# Because we know stopover and migration periods we can add this to inform calibration period
# get flight duration
ID.list2 = list.files(Activity_path,pattern="_act",include.dirs=T) # this is the PAMLr output file of flight classification
ID.list2
#Select the corect twl file
ID2 = ID.list2[9]
ID2

timetable <- read.csv(paste0(Activity_path,ID2))

MigTable<- timetable

timetable<- timetable[timetable$Duration..h.>2,]

StopTable<- data.frame("start" = timetable$end[-length(timetable$end)],
                       "end" = timetable$start[-1])

StopTable$duration =  difftime(format(StopTable$end, format = "%Y-%m-%d %H:%M:%S", tz = tz),
                               format(StopTable$start, format = "%Y-%m-%d %H:%M:%S", tz = tz), 
                               tz="GMT", units = "hours")

# only count it as a stopover if the bird was there for more than 2 days
# stopover_timetable <- stopover_timetable[stopover_timetable$duration > 72,]
days <- 2
timetable <-  StopTable[StopTable$duration > days*24,] #  stopover_timetable #stopover_timetable[stopover_timetable$duration > days*24,]

if(MigTable$end[length(MigTable$end)]>end_twl){
  MigTable<- MigTable[seq(start_twl<MigTable$start[1],min(which(MigTable$end>end_twl))-1),]
}

############### Add flight times from activity timetable to twl file ######

twl$dt <- 0
twl$startMig<-0
twl$endMig<-0

origin = "1970-01-01"
for (i in  1:(nrow(MigTable))){
  # print(i)
  for (k in 2:nrow(twl)){
    if (MigTable$start[i] >= twl$Twilight[k] &  MigTable$start[i] <= twl$Twilight[k+1] & MigTable$end[i] <= twl$Twilight[k+1]){
      twl$dt[k] <-  twl$dt[k]+difftime(MigTable$end[i],MigTable$start[i], unit="hours")
      twl$startMig[k]<-as.numeric(strptime(MigTable$start[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
      twl$endMig[k]<- as.numeric(strptime(MigTable$end[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
    }
    if (MigTable$start[i] >= twl$Twilight[k] &  MigTable$start[i] <= twl$Twilight[k+1]  & MigTable$end[i] >= twl$Twilight[k+1]  & MigTable$Duration..h.[i]-(difftime(twl$Twilight[k+1],MigTable$start[i], unit="hours")) <= difftime(twl$Twilight[k+2],twl$Twilight[k+1], unit="hours")){
      
      twl$dt[k] <- twl$dt[k]+difftime(twl$Twilight[k+1],MigTable$start[i], unit="hours")
      twl$dt[k+1] <- twl$dt[k+1]+difftime(MigTable$end[i],twl$Twilight[k+1], unit="hours")
      twl$startMig[k]<- as.numeric(strptime(MigTable$start[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
      twl$endMig[k]<- as.numeric(twl$Twilight[k+1])
      twl$startMig[k+1]<- as.numeric(twl$Twilight[k+1])
      twl$endMig[k+1]<- as.numeric(strptime(MigTable$end[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
    }
    if (MigTable$start[i] >= twl$Twilight[k] &  MigTable$start[i] <= twl$Twilight[k+1] & MigTable$end[i] >= twl$Twilight[k+1]  & MigTable$Duration..h.[i]-(difftime(twl$Twilight[k+1],MigTable$start[i], unit="hours")) >= difftime(twl$Twilight[k+2],twl$Twilight[k+1], unit="hours")){
      twl$dt[k] <- twl$dt[k]+difftime(twl$Twilight[k+1],MigTable$start[i], unit="hours")
      twl$dt[k+1] <- twl$dt[k+1]+difftime(twl$Twilight[k+2],twl$Twilight[k+1], unit="hours")
      twl$dt[k+2] <- twl$dt[k+2]+difftime(MigTable$end[i],twl$Twilight[k+2], unit="hours")
      twl$startMig[k]<- as.numeric(strptime(MigTable$start[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
      twl$endMig[k]<- twl$Twilight[k+1]
      twl$startMig[k+1]<- as.numeric(twl$Twilight[k+1])
      twl$endMig[k+1]<- as.numeric(twl$Twilight[k+2])
      twl$startMig[k+2]<- as.numeric(twl$Twilight[k+2])
      twl$endMig[k+2]<- as.numeric(strptime(MigTable$end[i],format = "%Y-%m-%d %H:%M:%S", tz=tz))
    }
  }
}



twl$dt[twl$dt==0] <- 0.5 

# Get rid of any extra lines all over the place
twl <- subset(twl, (Twilight >= as.POSIXct(timetable$start[1],tz=tz) & 
                      Twilight <= as.POSIXct(timetable$end[nrow(timetable)],tz=tz) & Deleted == F))

twl$idx<- c(1:length(twl$dt))

index <- which(timetable$duration == max(timetable$duration)) # 5# 

#find start of stopover
start_index <- min(which(twl$Twilight >= as.POSIXct(timetable$start[index], tz=tz)))
end_index <- min(which(twl$Twilight >= as.POSIXct(timetable$end[index], tz=tz)))

# do the fireplot
(zenith_sd <- findHEZenith(twl, tol = 0.01, range=c(start_index,end_index-55)) )

#Set up connection to wind Data Base
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user=user, password=password, 
                 dbname = dbname, host = host)


Site <- rep(0, len=length(twl$Twilight))
stp <- timetable 
for (i in 1:nrow(stp)){
  Site <- ifelse(twl$Twilight > as.POSIXct(stp$start[i],tz=tz) & twl$Twilight < as.POSIXct(stp$end[i],tz=tz), 
                 i,Site)
}
grouped <- Site>0

colours <- c("black",colorRampPalette(c("blue","yellow","red"))(max(Site)))

# Create a vector which indicates which numbers sites as 111123444444567888889
twl$group <- makeGroups(Site)
g <- twl$group 
gr<- g

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


data <- data.frame(unique(cbind(from=twl$group[1:(length(twl$group)-1)], 
                                to=twl$group[2:length(twl$group)],
                                duration=twl$dt[1:(length(twl$group)-1)],
                                stopover=Site[1:(length(twl$group)-1)],
                                MigStart=twl$startMig[1:(length(twl$group)-1)],
                                MigEnd=twl$endMig[1:(length(twl$group)-1)],
                                idx=twl$idx[1:(length(twl$group)-1)])))

data <- data[data$from != data$to,]
head(data)

dt<- data.frame(cbind(data$from,data$to,data$duration,data$stopover))

names(dt)[1] <- "from"
names(dt)[2] <- "to"
names(dt)[3] <- "duration"
names(dt)[4] <- "idx"

Mig_start<-data$MigStart
Mig_end<- data$MigEnd
twl_idx<- data$idx

idx<-which(Mig_start ==0)
twl_idx<- twl_idx[idx]

Mig_start[idx]<- as.numeric(as.POSIXct(twl$Twilight[twl_idx],format = "%Y-%m-%d %H:%M:%S", tz=tz))
Mig_end[idx]<- as.numeric(as.POSIXct(twl$Twilight[twl_idx+1],format = "%Y-%m-%d %H:%M:%S", tz=tz))



wind_t<- (Mig_start+Mig_end)/2

behaviour <- c()
for (i in 1:max(g)){
  behaviour<- c(behaviour, which(g==i)[1])
}

stationary <- grouped[behaviour]
sitenum <- cumsum(stationary==T)
sitenum[stationary==F] <- 0


# Here starts the normal SGAT analyses
tol=0.12 #adjust to a larger value to extrapolate more during the equinoq times

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

beta  <- c(4, 0.1) # c(2.2, 0.03)
matplot(0:150, dgamma(0:150, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")
# Define known locations and set them as fixed locations. Ficed locations are not changed in BI
fixedx <- rep_len(FALSE, length.out = nrow(x0))
fixedx[1] <- TRUE # first stationary site - fixed
fixedx[nrow(x0)] <- T # last stationary site - fixed
x0[fixedx,1] <- lon.calib
x0[fixedx,2] <- lat.calib
z0 <- trackMidpts(x0)

# plot it
xlim <- range(x0[,1]+c(-5,5))
ylim <- range(0,55)
plot(x0, type = "n", xlab = "", ylab = "" ,xlim=xlim, ylim=ylim)
plot(wrld_simpl, col = "grey95", add = T)
points(x0, pch=19, col="cornflowerblue", type = "o")
points(x0, pch=21, bg=colours[sitenum+1], 
       col="black",cex = ifelse(sitenum>0, 3, 0))
text(x0[sitenum>0,1], x0[sitenum>0,2],  sitenum[sitenum>0])
box()


index = ifelse(stationary, 1, 2)

mask <- earthseaMask(xlim, ylim, n = 2,index=index)

logp <- function(p) {
  f <- mask(p)
  ifelse(!f | is.na(f), -1000, f)
}

# define the error shape
x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(x0))
z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(z0))

MigTimeTable<- cbind(Mig_start, Mig_end)




##' Read wind data from data base
#Unique time stams
time_list <- dbGetQuery(con, paste0("SELECT DISTINCT tms FROM", 
                                    "\"WindTable\"", ";"))

#Find closest time stamp available in db for data retrieval
#compute step size in t stamps stored in db
step_size <- time_list$tms[2] - time_list$tms[1]
start_t_idx <- time_list$tms[1]/step_size

#find closest time stamp to start times of each migratory bout
start_mig_t_idx <- round(Mig_start/step_size)
start_list_idx <- start_mig_t_idx - start_t_idx
tm_start <- time_list$tms[start_list_idx]

#find closest time stamp to end times of each migratory bout
end_mig_t_idx <- round(Mig_end/step_size)
end_list_idx <- end_mig_t_idx - start_t_idx
tm_end <- time_list$tms[end_list_idx]

#Find all time stamps between each start and end time
alt_t_idx<- seq(1,length(tm_start))
tm<- lapply(alt_t_idx, FUN= function(x){
  ifelse(tm_start[x]<tm_end[x], time_list[t_span<-which(time_list>= tm_start[x] & time_list<= tm_end[x]),1], time_list[t_span<-which(time_list>= tm_end[x] & time_list<= tm_start[x]),1])
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
  alt <- ifelse(alt <= 625, 650, alt)
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
tm_str<- substr(str_1,1,nchar(str_1)-4)

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

#Sort problems out with duplicate time stamps
if(length(anyDuplicated(tm)) == 1){
  duplicateIDX<- which(duplicated(tm) == T)
  for(i in duplicateIDX) { 
    newrow<- test_wind[(((i*th_listconcat)-th_listconcat)+1):(i*th_listconcat),]
    row<-(((i*th_listconcat)-th_listconcat)+1)
    test_wind[seq(row+length(newrow[,1]), nrow(test_wind)+length(newrow[,1])),]<- test_wind[seq(row, nrow(test_wind)),]
    test_wind[seq(row, (row+length(newrow[,1])-1)),]<- newrow
  }
}

#generate subset of the windtable, to speed up wind data retrieval
#list of times en route in lenght of available positions
tms<- rep(tm[1], length.out = th_listconcat)
for (i in position_idx[-1]){
  tms<-append(tms, rep(tm[i],length.out=th_listconcat), after = length(tms))
}


#list of winds at required heights for each position in u direction
u_list<- test_wind[,paste("u_",alt[1], collapse="",sep="")][which(test_wind$tms == tm[1])]

for (i in c(2: length(position_idx))){
  if(any(i == duplicateIDX)==T ){
    u_list<-append(u_list, test_wind[,paste("u_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])[1:length(newrow[,1])]], after = length(u_list))
  }
  else if(any(i == (duplicateIDX+1))==T){ 
    u_list<-append(u_list, test_wind[,paste("u_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])[1:length(newrow[,1])]], after = length(u_list))
  }
  else{
    u_list<-append(u_list, test_wind[,paste("u_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])], after = length(u_list))
  }
  
}

#list of winds at required heights for each position in v direction
v_list<- test_wind[,paste("v_",alt[1], collapse="",sep="")][which(test_wind$tms == tm[1])]

for (i in c(2: length(position_idx))){
  if(any(i == duplicateIDX)==T){
    v_list<-append(v_list, test_wind[,paste("v_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])[1:length(newrow[,1])]], after = length(v_list))
  }
  else if(any(i == (duplicateIDX+1))==T){ 
    v_list<-append(v_list, test_wind[,paste("v_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])[1:length(newrow[,1])]], after = length(v_list))
  }
  else{
    v_list<-append(v_list, test_wind[,paste("v_",alt[i], collapse="",sep="")][which(test_wind$tms == tm[i])], after = length(v_list))
  }
  
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
  sRad<- 200
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
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2) * ms2kmh}),recursive = T, use.names = T)
  
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
  #alpha <- ifelse(alpha <= (pi/2), alpha, ifelse(alpha<=pi, pi - alpha, alpha-pi))
  
  #retrieve the direction of the vector describing the wind compenstaion a bird has to fly in order to fly to its desired destination
  dir_comp <- ifelse(bird_vec_angle > wind_directions, bird_vec_angle - (pi/2), bird_vec_angle - ((3/2) * pi))
  
  dir_comp <- ifelse(dir_comp > 0, dir_comp, dir_comp + 2 * pi)
  
  #compute magnitude of wind compensation vector
  mag_comp <- cos(pi - (pi/2) - alpha[3]) * wind_speeds
  
  #After computing the compensation vector, the wind add parallel to the bird flying vector is computed
  wind_add <- ifelse(alpha<=pi | alpha>= ((2/3)*pi), sqrt(wind_speeds^2 - mag_comp^2),-sqrt(wind_speeds^2 - mag_comp^2))
  
  #Sum of wind partition parallel to bird flying vector and bird speed
  new_b_speed <- (cos(alpha)*wind_speeds) + ifelse(is.nan(sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2)), 0,sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2))
  
  #sigma of bird speed according to Bruderer er al, Bolch & Brudrer and Liechti & Bruderer
  bird_speed_sigma <- 3.5 
  
  #standard deviation of brid ground speed (Bird airspeed + wind speed)
  std <- bird_speed_sigma + sigma_wind_spd
  
  #Beta applied to describe gamma distribution
  beta_dist <- new_b_speed/(std^2)
  
  
  #Alpha applied to describe gamma distribution
  alpha_dist <- new_b_speed * beta_dist
  spd<- ifelse(is.nan(spd),0, spd)
  
  #compute probabilities of each point along the migration route
  prob_list <- dgamma(spd, alpha_dist, beta_dist, 
                      log = T)
  prob_list<- ifelse(is.na(prob_list),-1000,prob_list)
  
  return(prob_list)
}


#change the SpeedGammaModel in SGAT to make sure our Movement model defined here is called when probabilities are computed
#remember to copy paste the code from SpeedGammaModel.R and paste it between the {}
trace(speedGammaModel,edit = T)
#untrace(speedGammaModel)
trace(estelleMetropolis,edit=T)
untrace(estelleMetropolis)

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
                               dt = dt$duration, # This is the new part here - movment is allowed only for the number of hours that PAM analyses said the bird was flying
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
                               twilight.model = "ModifiedGamma", # this bit gere is what differs from the previously set model. "ModifiedGamma" will work always, but "Gamma" will fail when designation between stationary-movement phases are likely worong
                               alpha = alpha,
                               beta =  beta,
                               x0 = x0, # meadian point for each greoup (defined by twl$group)
                               z0 = z0, # middle points between the x0 points
                               zenith = zenith0+.5, # I often see that zenith0 for wintering period is typically a bit higher that zenith0 for breeding site, thus I manually add anywhere between 0.5 to 2 zenith0 - your model will also fail if the original zenith0 value is too low for the wintering period - just test adding a bit of the model fails
                               logp.x = logp, # land sea mask
                               dt = dt$duration,
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
sm_z <- locationSummary(fit$z,time = fit$model$time)
z<- cbind(sm_z$Lon.mean,sm_z$Lat.mean)
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

plot(wrld_simpl, xlim=xlim, ylim=c(-5,50))
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

#get harmonic mean of solar positioning
idx<-c(1:length(x[,1]))

hm_list<- unlist(lapply(idx, FUN=function(p){
  harmonic.mean(logDF[p,],na.rm = T,zero = F)
}))

sm<- cbind(sm,hm_list)

## Visualization of results and data
x0 <- track
z0 <- trackMidpts(x0)
x<- cbind(sm$Lon.mean,sm$Lat.mean)
spd_z <- pmax.int(trackDist2(x, z), 1e-06)/dt$duration
spd_z0<- pmax.int(trackDist2(x0, z0), 1e-06)/dt$duration

bird_vec_angle<- directions(x0,z0)

prob_list_log<- wind_model(x0,z0, bird_vec_angle, dt$duration)
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
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2) * ms2kmh}),recursive = T, use.names = T)
  
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
  new_b_speed <- (cos(alpha)*wind_speeds) + ifelse(is.nan(sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2)), 0,sqrt((beta[1] * (1/beta[2]))^2 - mag_comp^2))
  
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

prob_list<- wind_model(x0,z0, bird_vec_angle, dt$duration)
prob_list[length(prob_list)+1]<- prob_list[1]
spd[length(spd)+1]<- 0
spd<- spd/3.6
sm<- cbind(sm, prob_list)
sm<- cbind(sm, spd)

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
  sRad<- 5000
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
  wind_speeds<- unlist(lapply(position_idx,FUN=function(x){sqrt(wind_speeds_xy[[x]][[1]]^2 + wind_speeds_xy[[x]][[2]]^2)}),recursive = T, use.names = T)
  
}

wind_speeds<- wind_model(x0,z0,bird_vec_angle,dt$duration)
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

wind_add<-wind_model(x0,z0,bird_vec_angle,dt$duration)

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
#


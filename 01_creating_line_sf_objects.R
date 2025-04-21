#load libraries
library(sf)
library(dplyr)
library(geodata)
library(tidyverse)
library(rnaturalearth)
library(gganimate)
library(raster)
library(ggmap)
library(smoothr)
library(lubridate)
library(strex)
library(ggnewscale)

#read in data
telem<-read.csv("./data/telem.total.csv")

telem$Date<-as.Date(telem$Date,format="%Y-%m-%d")
telem$Year<-year(telem$Date)

#add a season column for both the telemetry and removal dataframes
#define a function to identify the breeding season
season<-function(x){ifelse(
  x[['Date']] >= ymd(paste0(x[['Year']],"-01-01")) & 
    x[['Date']] <= ymd(paste0(x[['Year']],"-10-14")),
  paste(as.numeric(x[['Year']])-1,"-",x[['Year']],sep=""), 
  paste(x[['Year']],"-",as.numeric(x[['Year']])+1,sep=""))
}

#apply the function
telem$season<-apply(telem, 1, season)

#add a unique id-season identifier
telem<-telem %>% mutate(id_season = paste(telem$Snake_ID, "_",telem$season, sep=""))

#remove id-seasons with less than 1 observation
telem<-telem %>% group_by(id_season) %>% mutate(n=n()) %>% filter(n>1) %>% filter(id_season!="PYBI_110_2024-2025")

telem.ind<-telem %>% filter(Snake_ID=="PYBI_82")
# x<- c(telem.ind$lon)
# y<- c(telem.ind$lat)
# coord2<- cbind(x,y)
telem<-telem %>% arrange(id_season)
id.season<-unique(telem$id_season)

sf_lines<-list()
buff_lines<-list()
buff_area<-list()
for(i in 1:length(id.season)){
  telem.ind<-telem %>% filter(id_season==id.season[i]) #isolate the telemetry points for a single id-season
  sf_object <- st_as_sf(telem.ind, coords = c("Long", "Lat"),crs=4269) #create an sf line from that telemetry data
  sf_object<-st_transform(sf_object, 26917) #reproject the data
  L1 <- st_combine(sf_object) %>% st_cast("LINESTRING") #convert the sf object to a linestring
  L2<-smoothr::smooth(L1, method = "ksmooth",smoothness=2) #smooth the line
  line <- st_as_sf(L2) #reconvert to sf object
  buff_line <- st_buffer(line, dist = 50,endCapStyle="ROUND")
  line<-st_transform(line, crs=4326) #reproject line
  buff_line<-st_transform(buff_line, crs=4326) #reproject line
  sf_lines[i]<-line
  buff_lines[i]<-buff_line
  buff_area[i]<-(st_area(buff_line))*3.86102e-7 #area in sq.mi
}  
names(sf_lines)<-id.season
names(buff_lines)<-id.season
names(buff_area)<-id.season

#combining sf objects by id
# create a list of the same sf objects but with ids as names for grouping purposes
sf_lines2<-sf_lines
ids<-unique(telem$Snake_ID)
names(sf_lines2)<-ids

#convert the list of sf objects to a dataframe for grouping purposes
sf.df<-data.frame(value = sf_lines2, stringsAsFactors = T)
sf.df<-t(sf.df) #transpose df
sf.df<-as.data.frame(sf.df) #reconvert back to df
sf.df<-tibble::rownames_to_column(sf.df,"id_season") #turn row names into column
colnames(sf.df)[2]<-"geometry" #assign name to geometry column
sf.df$id_season<-id.season #assign id.season column
sf.df$id<-str_before_nth(id.season, "_", 2) #isolate just the id part of the id.season identifier
sf.df<-sf.df %>% group_by(id) #group by ids

#combining buffer objects by id
# create a list of the same sf objects but with ids as names for grouping purposes
buff_lines2<-buff_lines
names(buff_lines2)<-ids

#convert the list of buff objects to a dataframe for grouping purposes
buff.df<-data.frame(value = buff_lines2, stringsAsFactors = T)
buff.df<-t(buff.df) #transpose df
buff.df<-as.data.frame(buff.df) #reconvert back to df
buff.df<-tibble::rownames_to_column(buff.df,"id_season") #turn row names into column
colnames(buff.df)[2]<-"geometry" #assign name to geometry column
buff.df$id_season<-id.season #assign id.season column
buff.df$id<-str_before_nth(id.season, "_", 2) #isolate just the id part of the id.season identifier
buff.df<-buff.df %>% group_by(id) #group by ids

#loop to create a list of combined sf objects and buffers for each snake id of all sf objects
sf.combined<-list()
buff.combined<-list()
for(i in 1:length(ids)){
  sfs<-sf.df %>% filter(id==ids[i]) #create a dataframe of the data filtered by snake id
  buffs<-buff.df %>% filter(id==ids[i]) #create a dataframe of the data filtered by snake id
  sf.combined[[i]]<-sfs %>% st_sf(sf_column_name = 'geometry') #combine those dataframes into a single sf object 
  buff.combined[[i]]<-buffs %>% st_sf(sf_column_name = 'geometry') #combine those dataframes into a single sf object 
}

names(sf.combined)<-ids #assign names
names(buff.combined)<-ids #assign names

#combining areas by id
# create a list of the same sf objects but with ids as names for grouping purposes
buff_area2<-buff_area
names(buff_area2)<-ids

#convert the list of buff objects to a dataframe for grouping purposes
area.df<-data.frame(value = buff_area2, stringsAsFactors = T)
area.df<-t(area.df) #transpose df
area.df<-as.data.frame(area.df) #reconvert back to df
area.df<-tibble::rownames_to_column(area.df,"id_season") #turn row names into column
colnames(area.df)[2]<-"area.sqmi" #assign name to geometry column
area.df$id_season<-id.season #assign id.season column
area.df$id<-str_before_nth(id.season, "_", 2) #isolate just the id part of the id.season identifier
area.df<-area.df %>% group_by(id) #group by ids
area.summ<-area.df %>% group_by(id) %>% summarise(total.area=sum(area.sqmi))


#creating overall track lines for an animal's entire track
cum.lines<-list()
for(i in 1:length(ids)){
  telem.ind<-telem %>% filter(Snake_ID==ids[i]) #isolate the telemetry points for a single id-season
  sf_object <- st_as_sf(telem.ind, coords = c("Long", "Lat"),crs=4269) #create an sf line from that telemetry data
  sf_object<-st_transform(sf_object, 26917) #reproject the data
  L1 <- st_combine(sf_object) %>% st_cast("LINESTRING") #convert the sf object to a linestring
  L2<-smoothr::smooth(L1, method = "ksmooth",smoothness=2) #smooth the line
  line <- st_as_sf(L2) #reconvert to sf object
  line<-st_transform(line, crs=4326) #reproject line
  cum.lines[i]<-line
}  
names(cum.lines)<-ids


#create a custom color gradient to highlight tracks over time (length.out=maximum number of seasons a snake was tracked for)
cc <- scales::seq_gradient_pal("red", "white", "Lab")(seq(0,1,length.out=6)) 
cc2 <- scales::seq_gradient_pal("blue", "white", "Lab")(seq(0,1,length.out=6)) 

#set API
ggmap::register_google(key = "AIzaSyBeZDy6_OZBaN3ewQZxDmKgh8Xcq519IKQ")

#create basemap
base<-get_map(c(left=min(st_coordinates(cum.lines[names(cum.lines)=="PYBI_82"][[1]])[, 1])-.03,
                bottom=min(st_coordinates(cum.lines[names(cum.lines)=="PYBI_82"][[1]])[, 2]),
                right=max(st_coordinates(cum.lines[names(cum.lines)=="PYBI_82"][[1]])[, 1]),
                top=max(st_coordinates(cum.lines[names(cum.lines)=="PYBI_82"][[1]])[, 2])),
                source="google",maptype="satellite",zoom=12)

#create ggmap object from base map
base<-ggmap(base)
base

#plot data
#sf.combined[names(sf.combined)=="PYBI_83"][[1]][[2]][1] OPTION FOR PLOTTING A SINGLE SEASON'S DATA FOR A GIVEN INDIVIDUAL
base+
 #geom_sf(data = buff.combined[names(buff.combined)=="PYBI_83"][[1]],inherit.aes=FALSE,show.legend = FALSE, aes(color=id_season,geometry = geometry), linewidth=1)+ #option for showing buffer lines
  geom_sf(data = cum.lines[names(cum.lines)=="PYBI_83"][[1]],inherit.aes=FALSE,show.legend = FALSE, aes(color="blue",geometry = geometry), linewidth=1)+
  scale_colour_manual(values="blue")+
  # new_scale_color() + 
  # geom_sf(data = cum.lines[names(cum.lines)=="PYBI_83"][[1]],inherit.aes=FALSE,show.legend = FALSE, aes(color=id_season,geometry = geometry), linewidth=1)+
  # scale_colour_manual(values=cc2)+
  xlim(min(st_coordinates(cum.lines[names(cum.lines)=="PYBI_83"][[1]])[, 1])-.05,max(st_coordinates(cum.lines[names(cum.lines)=="PYBI_83"][[1]])[, 1])+.05)+
  ylim(min(st_coordinates(cum.lines[names(cum.lines)=="PYBI_83"][[1]])[, 2]-.02),max(st_coordinates(cum.lines[names(cum.lines)=="PYBI_83"][[1]])[, 2]))+
  labs(x="Longitude",y="Latitude",title="Denis Total Movement")+
  theme(plot.title = element_text(hjust = 0.5))
#annotate(geom="text", x=-81.55, y=26.04, label="PYBI_83 2021-2024",
#        color="white")

ggsave(filename="./images/lines_denis_total.png")

  
  
####################################################BUFFER CODE####################################################################

#create a line sf object from the simulated data
sf_object <- st_as_sf(telem.ind, coords = c("Long", "Lat"),crs=4269)

#convert sf object to utms
sf_object<-st_transform(sf_object, 26917)

#create a linestring object of the telemetry data
L1 <- st_combine(sf_object) %>% st_cast("LINESTRING")
L2<-smoothr::smooth(L1, method = "ksmooth",smoothness=2)
plot(L2, axes = TRUE)

#convert linestring to sf object
line <- st_as_sf(L2)
plot(line, axes = TRUE)
crs(line)

#buffer around the line
buff_line <- st_buffer(line, dist = 50,endCapStyle="ROUND")

#transform crs to utm
line<-st_transform(line, crs=4326)
buff_line<-st_transform(buff_line,crs=4326)
 
#check buffer crs
crs(buff_line)
crs(line)
plot(buff_line)

#calculate area of buffer in sq km
area<-(st_area(buff_line))/1000000

##plotting the data

#plot base map
base

#plot data
cc <- scales::seq_gradient_pal("blue", "white", "Lab")(seq(0,1,length.out=10)) #create a custom color gradient first
base+
  #geom_sf(data = sf.combined,inherit.aes=FALSE,show.legend = FALSE, aes(color=id_season,geometry = geometry), linewidth=1)+
  #scale_colour_manual(values=cc)+
geom_sf(data = buff_line,inherit.aes=FALSE,show.legend = FALSE,fill="red")+
  #geom_sf(data = line,inherit.aes=FALSE,show.legend = FALSE,color="red",linewidth=2)+
  # geom_sf(data = sf_lines[[1]],inherit.aes=FALSE,show.legend = FALSE,color="darkred",linewidth=1)+
  # geom_sf(data = sf_lines[[2]],inherit.aes=FALSE,show.legend = FALSE,color="red",linewidth=1)+
  # geom_sf(data = sf_lines[[3]],inherit.aes=FALSE,show.legend = FALSE,color="orangered",linewidth=1)+
  # geom_sf(data = sf_lines[[4]],inherit.aes=FALSE,show.legend = FALSE,color="orange",linewidth=1)+
  # geom_sf(data = sf_lines[[5]],inherit.aes=FALSE,show.legend = FALSE,color="gold",linewidth=1)+
  # geom_sf(data = sf_lines[[6]],inherit.aes=FALSE,show.legend = FALSE,color="yellow",linewidth=1)+
  # geom_sf(data = sf_lines[[7]],inherit.aes=FALSE,show.legend = FALSE,color="yellowgreen",linewidth=1)+
  # geom_sf(data = sf_lines[[8]],inherit.aes=FALSE,show.legend = FALSE,color="green",linewidth=1)+
  # geom_sf(data = sf_lines[[9]],inherit.aes=FALSE,show.legend = FALSE,color="darkgreen",linewidth=1)+
  # geom_sf(data = sf_lines[[10]],inherit.aes=FALSE,show.legend = FALSE,color="blue",linewidth=1)+
  xlim(-81.72,-81.62)+
  ylim(25.97,26.05)
  #annotate(geom="text", x=-81.55, y=26.04, label="PYBI_82 2021-2024",
   #        color="white")
p

ggsave(p,filename="./images/line_test_Dion.png")

plot(sf_lines[[1]])

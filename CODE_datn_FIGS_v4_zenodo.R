####Install packages####
#load relevant packages
install.packages("ggplot2")
install.packages("plyr")
install.packages("dplyr")
install.packages("reshape2")
install.packages("tidyr")
install.packages("zyp")
install.packages("dismo")
install.packages("gbm")
install.packages("ggmap")
install.packages("maptools")
install.packages("MASS")
install.packages("data.table")
install.packages("scales")
install.packages("lubridate")
install.packages("ggridges")
install.packages("viridis")
install.packages("RColorBrewer")
install.packages("gridExtra")
install.packages("spatstat")
install.packages("ggridges")
install.packages("viridis")
install.packages("cowplot")
install.packages("sf")
install.packages("ggmap")
install.packages("raster")
install.packages("maptools")
install.packages("maps")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("rgeos")
install.packages("bestNormalize")
install.packages("mgcv")

#load relevant packages####
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(zyp)
library(dismo)
library(gbm)
library(ggmap)
library(maptools)
library(MASS)
library(data.table)
library(scales)
library(lubridate)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(spatstat)
library(ggridges)
library(viridis)
library(cowplot)
library(sf)
library(ggmap)
library(raster)
library(maptools)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(bestNormalize)
library(mgcv)
library(ggnewscale)
library(fields)
library(wesanderson)
library(inlmisc)
library(forcats)
library(ggpubr)
library(ggbeeswarm)

#Fig 1: Ridgeplots####
Vanern.data<-fread("D:/DATN/RefitData_v6/Vanern.csv",stringsAsFactors = T)
Vanern.data$date<-as_date(Vanern.data$date)
Vanern.data<-Vanern.data[year!=1973]
Giles.data<-fread("D:/DATN/RefitData_v6/Giles.csv",stringsAsFactors = T)
Giles.data$date<-as_date(Giles.data$date)
Giles.data<-Giles.data[year!=1988&
                         year!=1989&
                         year!=2014]
Trout.data<-fread("D:/DATN/RefitData_v6/TroutLake.csv",stringsAsFactors = T)
Trout.data$date<-as_date(Trout.data$date)
Trout.data<-Trout.data[year!=1981]
Kinneret.data<-fread("D:/DATN/RefitData_v6/Kinneret.csv",stringsAsFactors = T)
Kinneret.data$date<-as_date(Kinneret.data$date)
Kinneret.data<-Kinneret.data[,c(1:7,9:27)]
Kinneret.data<-Kinneret.data[year!=2014]
ridge.data<-rbindlist(list(Vanern.data,Trout.data,Giles.data,Kinneret.data),fill=TRUE)
ridge.data<-ridge.data[is.finite(fit)==TRUE&
                         is.finite(voldays)==TRUE]


hyps<-fread("D:/DATN/DATA_datn_hypsModel_v6.csv",header=TRUE,stringsAsFactors=TRUE,sep=",")


#set density difference data
rm(densitydif.all)
#lakename<-droplevels(unique(ridge.data$lake)[1])

for (lakename in unique(ridge.data$lake)){
  fitmin<-min(ridge.data[lake==lakename]$fit,na.rm=TRUE)
  fitmax<-max(ridge.data[lake==lakename]$fit,na.rm=TRUE)
  rowlake<-which(hyps$lake==paste0(lakename,sep=""))
  maxdepth<-hyps[rowlake]$MaxDepth
  sa<-hyps[rowlake]$SA
  volume<-hyps[rowlake]$MeanDepth*sa
  
  ridge.data.lake<-ridge.data[lake==lakename]
ridge.data.lake<-ridge.data.lake[sample(nrow(ridge.data.lake),nrow(ridge.data.lake),prob=ridge.data.lake$voldays,replace=TRUE)]
  
  ridge.data.lake$voldayssum<-0
  medianyear<-median(ridge.data.lake$year)
  ridge.data.lake<-ridge.data.lake[is.finite(fit)==TRUE&
                                     is.finite(voldays)==TRUE]
  ridge.data.lake[year<medianyear]$voldayssum<-sum(ridge.data.lake[year<medianyear]$voldays)
  ridge.data.lake[year>medianyear]$voldayssum<-sum(ridge.data.lake[year>medianyear]$voldays)
  ridge.data.lake<-ridge.data.lake[,voldays2:=(voldays/voldayssum)]
  
  densityX<-density(ridge.data.lake$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8)$x
  densityA<-density(ridge.data.lake[year<medianyear]$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8#,
                    #weights=ridge.data.lake[year<medianyear]$voldays2
                    )$y
  densityB<-density(ridge.data.lake[year>medianyear]$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8#,
                    #weights=ridge.data.lake[year>medianyear]$voldays2
                    )$y
  densitydif<-data.table(densityX,densityA,densityB)
  densitydif$lake<-lakename
  densitydif<-densitydif[,dif:=100*(2*(densityB-densityA))/sum(densitydif$densityA,densitydif$densityB)]
  #densitydif<-densitydif[,dif:=1000000000*dif*volume/(max(ridge.data.lake$doy)-min(ridge.data.lake$doy))]
  densitydif$difsign<-(densitydif$dif)<0
  densitydif<-densitydif[order(densityX)]
  densitydif$difsigngroup<-rleid(densitydif$difsign)
  
  # if the merged dataset does exist, append to it
  if (exists("densitydif.all")){densitydif.all<-rbindlist(list(densitydif.all,densitydif),fill=TRUE)}
  # if the merged dataset doesn't exist, create it
  if (!exists("densitydif.all")){densitydif.all<-densitydif}
}


#densitydif.all<-rbind
densitydif.all$lake<-as.factor(densitydif.all$lake)
densitydif.all<-densitydif.all[order(lake,densityX)]
head(densitydif.all)

densitydif.all[,.(turnover=sum(abs(dif))/2),.(lake)]


#set ridge data
ridge.data<-ridge.data[order(lake,year)]
ridge.data$year<-as.factor(ridge.data$year)
ridge.data[is.integer(as.numeric(year)/5)==FALSE]$year<-droplevels(ridge.data[is.integer(as.numeric(year)/5)==FALSE]$year)

#find non-overlap values for each lake
DATN_predict.m.mean<-fread("D:/DATN/DATA_predict.m.mean_v51.csv")
head(DATN_predict.m.mean)
mean(DATN_predict.m.mean[lake=="Vanern"&smd==1&zmd==1]$nov.std.predict.resid)
mean(DATN_predict.m.mean[lake=="Kinneret"&smd==1&zmd==1]$nov.std.predict.resid)
mean(DATN_predict.m.mean[lake=="Giles"&smd==1&zmd==1]$nov.std.predict.resid)
mean(DATN_predict.m.mean[lake=="TroutLake"&smd==1&zmd==1]$nov.std.predict.resid)

a<-ggplot(ridge.data[lake=="Vanern"], aes(x=fit,y=year,fill=..x..,
))+ 
  scale_fill_viridis(option = "B",breaks=c(.3,.4,.5,.6,.7),limits=c(-4,31.5))+
  geom_density_ridges_gradient(scale=3)+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Vanern"]$fit),
                                            max(ridge.data[lake=="Vanern"]$fit)))+
  scale_y_discrete(breaks=c(1960,1971,1980,1990,2000,2010))+
  theme_bw()+
  facet_wrap(~lake,nrow=1,scales="free")+
  theme(legend.position = "none",
        text = element_text(size=12),
        axis.text=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y="Year")

b<-ggplot()+
  geom_ribbon(data=densitydif.all[lake=="Vanern"],
              aes(ymin=0,
                  ymax=dif,
                  x=densityX,
                  fill=(dif<0),
                  group=difsigngroup))+
  geom_point(data=densitydif.all[lake=="Vanern"],
             aes(y=dif,
                 x=densityX),
             size=0.3)+
  scale_fill_manual(values=c("#43a2ca","#e34a33"))+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Vanern"]$fit),
                                            max(ridge.data[lake=="Vanern"]$fit)))+
  scale_y_continuous(limits=c(min(densitydif.all$dif),max(densitydif.all$dif)))+
  theme_bw()+
  annotate("text",x=16.5,y=-0.3,label="Thermal non-overlap = 6.9%",size=3.3)+
  labs(x="Temperature (°C)",
       #=expression(atop("Change in average","daily volume ("~m^3*")")))+
       y="Change in volume 
(% of total lake volume)")+
  guides(fill=FALSE)

c<-ggplot(ridge.data[lake=="Kinneret"], aes(x=fit,y=year,fill=..x..))+ 
  scale_fill_viridis(option = "B",breaks=c(.3,.4,.5,.6,.7),limits=c(-4,31.5))+
  geom_density_ridges_gradient(scale=3)+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Kinneret"]$fit),
                                            max(ridge.data[lake=="Kinneret"]$fit)))+
  scale_y_discrete(breaks=c(1960,1971,1980,1990,2000,2010))+
  theme_bw()+
  facet_wrap(~lake,nrow=1,scales="free")+
  theme(legend.position = "none",
        text = element_text(size=12),
        axis.text=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())+
  labs(y="Year")

d<-ggplot()+
  geom_ribbon(data=densitydif.all[lake=="Kinneret"],
              aes(ymin=0,
                  ymax=dif,
                  x=densityX,
                  fill=(dif<0),
                  group=difsigngroup))+
  geom_point(data=densitydif.all[lake=="Kinneret"],
             aes(y=dif,
                 x=densityX),
             size=0.3)+
  scale_fill_manual(values=c("#43a2ca","#e34a33"))+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Kinneret"]$fit),
                                            max(ridge.data[lake=="Kinneret"]$fit)))+
  scale_y_continuous(limits=c(min(densitydif.all$dif),max(densitydif.all$dif)))+
  theme_bw()+
  annotate("text",x=26,y=-0.3,label="Thermal non-overlap = 12.1%",size=3.3)+
  theme(axis.title.y=element_blank())+
  labs(x="Temperature (°C)",
       y="Change")+
  guides(fill=FALSE)

e<-ggplot(ridge.data[lake=="Giles"], aes(x=fit,y=year,fill=..x..))+ 
  scale_fill_viridis(option = "B",breaks=c(.3,.4,.5,.6,.7),limits=c(-4,31.5))+
  geom_density_ridges_gradient(scale=1.53)+#,quantile_lines=TRUE,quantiles=2)+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Giles"]$fit),
                                            max(ridge.data[lake=="Giles"]$fit)))+
  scale_y_discrete(breaks=c(1960,1971,1980,1990,2000,2010))+
  theme_bw()+
  facet_wrap(~lake,nrow=1,scales="free")+
  theme(legend.position = "none",
        text = element_text(size=12),
        axis.text=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y="Year")

f<-ggplot()+
  geom_ribbon(data=densitydif.all[lake=="Giles"],
              aes(ymin=0,
                  ymax=dif,
                  x=densityX,
                  fill=(dif<0),
                  group=difsigngroup))+
  geom_point(data=densitydif.all[lake=="Giles"],
             aes(y=dif,
                 x=densityX),
             size=0.3)+
  scale_fill_manual(values=c("#43a2ca","#e34a33"))+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="Giles"]$fit),
                                            max(ridge.data[lake=="Giles"]$fit)))+
  scale_y_continuous(limits=c(min(densitydif.all$dif),max(densitydif.all$dif)))+
  theme_bw()+
  annotate("text",x=20,y=-0.3,label="Thermal non-overlap = 8.1%",size=3.3)+
  labs(x="Temperature (°C)",
       #=expression(atop("Change in average","daily volume ("~m^3*")")))+
       y="Change in volume 
(% of total lake volume)")+
  guides(fill=FALSE)

g<-ggplot(ridge.data[lake=="TroutLake"], aes(x=fit,y=year,fill=..x..))+ 
  scale_fill_viridis(option = "B",breaks=c(.3,.4,.5,.6,.7),limits=c(-4,31.5))+
  geom_density_ridges_gradient(scale=3)+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="TroutLake"]$fit),
                                            max(ridge.data[lake=="TroutLake"]$fit)))+
  scale_y_discrete(breaks=c(1960,1971,1980,1990,2000,2010))+
  theme_bw()+
  facet_wrap(~lake,nrow=1,scales="free")+
  theme(legend.position = "none",
        text = element_text(size=12),
        axis.text=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())+
  labs(y="Year")

h<-ggplot()+
  geom_ribbon(data=densitydif.all[lake=="TroutLake"],
              aes(ymin=0,
                  ymax=dif,
                  x=densityX,
                  fill=(dif<0),
                  group=difsigngroup))+
  geom_point(data=densitydif.all[lake=="TroutLake"],
             aes(y=dif,
                 x=densityX),
             size=0.3)+
  scale_fill_manual(values=c("#43a2ca","#e34a33"))+
  scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data[lake=="TroutLake"]$fit),
                                            max(ridge.data[lake=="TroutLake"]$fit)))+
  scale_y_continuous(limits=c(min(densitydif.all$dif),max(densitydif.all$dif)))+
  theme_bw()+
  annotate("text",x=19,y=-0.3,label="Thermal non-overlap = 8.9%",size=3.3)+
  theme(axis.title.y=element_blank())+
  labs(x="Temperature (°C)",
       y="Change")+
  guides(fill=FALSE)

plot_grid(a,c,b,d,e,g,f,h,nrow=4,ncol=2,rel_heights = c(8,3,6,3),align="v")

###Fig 1 SI: Make a ridge plot for each lake####
lakenames_lil<-fread("D:/DATN/LakeNames_lil_v1.csv")
lakenames_big<-fread("D:/DATN/LakeNames_big_v1.csv")
lakenames<-c(lakenames_lil$LakeNames_lil,lakenames_big$LakeNames_big)
lakenames<-lakenames[lakenames!="Chignik"]
lakenames<-lakenames[lakenames!="Pear"]
lakenames<-lakenames[lakenames!="Lynx"]
lakenames<-lakenames[lakenames!="Hidden"]
lakenames<-lakenames[lakenames!="Hawley"]
lakenames<-lakenames[order(lakenames)]

hyps<-fread("D:/DATN/DATA_datn_hypsModel_v6.csv",header=TRUE,stringsAsFactors=TRUE,sep=",")
DATN_predict.m.mean<-fread("D:/DATN/DATA_predict.m.mean_v51.csv")
head(DATN_predict.m.mean)
lakenames<-unique(DATN_predict.m.mean$lake)

#calculate density ridges and compile plots
for (lakename in unique(lakenames)){
  ridge.data<-fread(paste0("D:/DATN/RefitData_v6/",lakename,".csv"),stringsAsFactors = T)
  ridge.data<-ridge.data[is.finite(fit)==TRUE&
                           is.finite(voldays)==TRUE]
  
  rowlake<-which(hyps$lake==paste0(lakename,sep=""))
  maxdepth<-hyps[rowlake]$MaxDepth
  sa<-hyps[rowlake]$SA
  volume<-hyps[rowlake]$MeanDepth*sa
  ed<-round(mean(DATN_predict.m.mean[lake==lakename&smd==1&zmd==1]$nov.std.predict.resid),2)
  
  fitmin<-floor(min(ridge.data[lake==lakename]$fit,na.rm=TRUE)*10)/10
  fitmax<-ceiling(max(ridge.data[lake==lakename]$fit,na.rm=TRUE)*10)/10
  ridge.data$voldayssum<-0
  medianyear<-median(ridge.data$year)
  ridge.data<-ridge.data[is.finite(fit)==TRUE&
                                     is.finite(voldays)==TRUE]
  ridge.data[year<medianyear]$voldayssum<-sum(ridge.data[year<medianyear]$voldays)
  ridge.data[year>medianyear]$voldayssum<-sum(ridge.data[year>medianyear]$voldays)
  ridge.data<-ridge.data[,voldays2:=(voldays/voldayssum)]
  
  densityX<-density(ridge.data$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8)$x
  densityA<-density(ridge.data[year<medianyear]$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8,
                    weights=ridge.data[year<medianyear]$voldays2)$y
  densityB<-density(ridge.data[year>medianyear]$fit,
                    from=fitmin,
                    to=fitmax,
                    n=((fitmax-fitmin)*10+1),
                    adjust=1.8,
                    weights=ridge.data[year>medianyear]$voldays2)$y
  densitydif<-data.table(densityX,densityA,densityB)
  densitydif$lake<-lakename
  densitydif<-densitydif[,dif:=100*(2*(densityB-densityA))/sum(densitydif$densityA,densitydif$densityB)]
  #densitydif<-densitydif[,dif:=1000000000*dif*volume/(max(ridge.data$doy)-min(ridge.data$doy))]
  densitydif$difsign<-(densitydif$dif)<0
  densitydif<-densitydif[order(densityX)]
  densitydif$difsigngroup<-rleid(densitydif$difsign)
  
  #make year into a factor for ridgeplot
  ridge.data<-ridge.data[order(lake,year)]
  ridge.data$year<-as.factor(ridge.data$year)
  
  a<-ggplot(data=ridge.data,aes(x=fit,y=year,fill=..x..))+ 
    scale_fill_viridis(option = "B",breaks=c(.3,.4,.5,.6,.7),limits=c(fitmin,fitmax))+
    geom_density_ridges_gradient(scale=3)+
    scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data$fit),
                                              max(ridge.data$fit)))+
    scale_y_discrete()+
    theme_bw()+
    facet_wrap(~lake,nrow=1,scales="free")+
    theme(legend.position = "none",
          text = element_text(size=12),
          axis.text=element_text(size=8),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())+
    labs(y="Year")
  
  b<-ggplot()+
    geom_ribbon(data=densitydif,
                aes(ymin=0,
                    ymax=dif,
                    x=densityX,
                    fill=(dif<0),
                    group=difsigngroup))+
    geom_point(data=densitydif,
                aes(y=dif,
                    x=densityX,
                    fill=(dif<0),
                    group=difsigngroup),
               size=0.8)+
    scale_fill_manual(values=c("#43a2ca","#e34a33"))+
    scale_x_continuous(expand=c(0,0),limits=c(min(ridge.data$fit),
                                              max(ridge.data$fit)))+
    theme_bw()+
    labs(x="Temperature (°C)",
         y="Change in volume 
(% of total lake volume)")+
    ggtitle(paste0("Thermal non-overlap = ",ed,"%"))+
    guides(fill=FALSE)
  
  jpeg(paste0("D:/DATN/Supplement/S3/S3_Ridgeplot_",lakename,"_v12.jpg"),width=400,height=100+length(unique(ridge.data$year))*12,units="px")
  print(plot_grid(a,b,nrow=2,ncol=1,rel_heights = c(2,1,1),align="v"))
  dev.off()
  
  }

#Fig 3####
DATN_predict.m.mean<-fread("D:/DATN/DATA_predict.m.mean_v51.csv")

DATN_predict.m.mean.limit<-DATN_predict.m.mean[smd==20&zmd==20,.(nov=mean(nov),
                                                                 nov.std=mean(nov.std.predict.resid),
                                                                 nov.randnull=mean(nov.randnull),
                                                                 nov.std.raw=mean(nov-nov.randnull),
                                                                 lat=mean(lat),
                                                                 lon=mean(lon)),.(lake)]
DATN_predict.m.mean.limit$nov<-DATN_predict.m.mean.limit$nov/100
DATN_predict.m.mean.limit$nov.std<-DATN_predict.m.mean.limit$nov.std/100
DATN_predict.m.mean.limit$nov.std.raw<-DATN_predict.m.mean.limit$nov.std.raw/100
DATN_predict.m.mean.limit[nov.std.raw<0]$nov.std.raw<-0
DATN_predict.m.mean.limit$nov.randnull<-DATN_predict.m.mean.limit$nov.randnull/100

summary(DATN_predict.m.mean.limit)


DATN_predict.m.mean<-DATN_predict.m.mean[smd==1&zmd==1,.(nov=mean(nov),
                                                         nov.std=mean(nov.std.predict.resid),
                                                         nov.randnull=mean(nov.randnull),
                                                         nov.std.raw=mean(nov-nov.randnull),
                                                         lat=mean(lat),
                                                         lon=mean(lon)),.(lake)]
#EDdata.lakemeans<-fread("D:/DATN/DATN_EDdata.lakemeans_v48.csv",stringsAsFactors = T)
DATN_predict.m.mean$nov<-DATN_predict.m.mean$nov/100
DATN_predict.m.mean$nov.std<-DATN_predict.m.mean$nov.std/100
DATN_predict.m.mean$nov.std.raw<-DATN_predict.m.mean$nov.std.raw/100
DATN_predict.m.mean[nov.std.raw<0]$nov.std.raw<-0
DATN_predict.m.mean$nov.randnull<-DATN_predict.m.mean$nov.randnull/100
summary(DATN_predict.m.mean)

fwrite(DATN_predict.m.mean,"D:/DATN/DATA_DATN_predict.m.mean_v100.csv")

DATN_predict.m.mean.l<-data.table(value=c(DATN_predict.m.mean$nov,
                                          DATN_predict.m.mean$nov.randnull,
                                          DATN_predict.m.mean$nov.std.raw,
                                          DATN_predict.m.mean$nov.std),
                                  type=rep(c("Raw change over time","Difference among randomly 
selected years","nov.std.raw","nov.std"),each=139))
summary(DATN_predict.m.mean.l)
fwrite(DATN_predict.m.mean.l,"D:/DATN/DATA_DATN_predict.m.mean.l_v1.csv")



lakeinfo<-fread("D:/DATN/DATA_datn_HydroLakesInfo_v2.csv",header=TRUE,stringsAsFactors = TRUE) #needs editing

OVdata.spatial<-fread("D:/DATN/DATN_OVdata.spatial_v52.csv",stringsAsFactors = T)
OVdata.spatial<-OVdata.spatial[,OVmean.std:=ov+(1-ov.null)]
OVdata.spatial<-OVdata.spatial[,NOVmean.std:=1-OVmean.std]
OVdata.spatial[NOVmean.std<0]$NOVmean.std<-0
hist(OVdata.spatial$NOVmean.std)


OVdata.spatial$NOVmean.std.rank<-rank(OVdata.spatial$NOVmean.std)

#test to see if there is a simpler way
OVdata.spatiotempA<-data.table(OVdata.spatial$NOVmean.std)
OVdata.spatiotempA<-OVdata.spatiotempA[,OVtype:="Differences among lakes (n = 9,591)"]
summary(OVdata.spatiotempA)
OVdata.spatiotempC<-data.table(DATN_predict.m.mean$nov.randnull)
OVdata.spatiotempC<-OVdata.spatiotempC[,OVtype:="Comparing between two 
randomly selected year groups 
(n = 139)"]
OVdata.spatiotemp.simp<-rbind(OVdata.spatiotempA,
                         OVdata.spatiotempC,
                         use.names=FALSE)
OVdata.spatiotemp.simp$OVtype<-as.factor(OVdata.spatiotemp.simp$OVtype)

OVdata.spatiotempB<-data.table(DATN_predict.m.mean$nov.std)
OVdata.spatiotempB<-OVdata.spatiotempB[,OVtype:="Comparing between first and 
second half of time series (n = 139)"]

#find the segment data for the colored distribution for thermal nvelty
OVdata.spatiotemp.ribbon<-data.table(x=density((OVdata.spatiotempB$V1),n=1024,adjust=2)$x,
                                     y=density((OVdata.spatiotempB$V1),n=1024,adjust=2)$y)
OVdata.spatiotemp.ribbon.limit<-data.table(x=density((DATN_predict.m.mean.limit$nov.std),n=1024,adjust=2)$x,
                                     y=density((DATN_predict.m.mean.limit$nov.std),n=1024,adjust=2)$y)


world <- ne_countries(scale = "medium", returnclass = "sf")
world<-st_transform(world,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

DATN_predict.m.mean<-st_as_sf(DATN_predict.m.mean,coords = c("lon", "lat"),crs=4326)
DATN_predict.m.mean<-st_transform(DATN_predict.m.mean,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

scaleFUN0.prop <- function(x) sprintf("%.0f", x)
scaleFUN0 <- function(x) sprintf("%.0f", x*100)
scaleFUN1 <- function(x) sprintf("%.1f", x*100)
scaleFUN2 <- function(x) sprintf("%.2f", x)
scaleFUN3 <- function(x) sprintf("%.3f", x)

library(png)
library(grid)
img <- readPNG("D:/DATN/viridis.png")
g <- rasterGrob(img, interpolate=TRUE)

summary(OVdata.spatiotemp.ribbon)
summary(DATN_predict.m.mean)

unname(quantile(x=DATN_predict.m.mean$nov,probs=0))
unname(quantile(x=OVdata.spatiotemp.ribbon$x,probs=0))

#Extended data figure
ggplot()+
  geom_density(data=DATN_predict.m.mean.l[type=="Raw change over time"|type=="Difference among randomly 
selected years"],
               aes(x=value,fill=type),colour=NA,alpha=0.5)+
  scale_fill_manual(values=c("grey","black"))+
  scale_x_continuous(expand=c(0,0),
                     labels = scaleFUN0)+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Thermal habitat (% non-overlap)",
       y="Kernel Density Estimation",
       fill=NULL)+
  theme_bw()+
  theme(legend.position = c(0.6, 0.60))


a<-ggplot(data=world)+
  geom_sf(fill="darkgrey",colour="grey")+
  geom_sf(data=DATN_predict.m.mean,
          aes(fill=rank(nov.std)),
          size=2.6,
          shape = 23)+
  scale_fill_viridis(breaks=c(1,35,70,104,138),labels=c("2.8 (Min)","4.5 (Q1)","5.2 (Median)","6.9 (Q3)", "20.0 (Max)"))+
  theme_bw()+
  theme(text = element_text(size=10),
        panel.grid.major = element_line(colour="lightgrey"),
        #panel.grid.minor = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        axis.title.y =element_blank(),
        # panel.background = element_rect(fill = "lightblue",
        #                                 colour = "lightblue",
        #                                 size = 0.5, linetype = "solid"),
        legend.position = c(0.15, 0.48),
        legend.key.size = unit(0.3, "cm")
  )+
  labs(fill="Thermal habitat
change over time 
(% non-overlap)")

b<-ggplot()+
  geom_density(data=OVdata.spatiotempA,aes(x=V1,fill=OVtype),colour=NA,alpha=0.6)+
  geom_segment(data=OVdata.spatiotemp.ribbon.limit[x>0],aes(x=x,xend=x,y=0,yend=y/1.6,colour=(x)),size=.2,alpha=.05)+
  scale_colour_gradientn(colours=c(viridis_pal()(72)),
                         values=c(0,gam(sort(DATN_predict.m.mean.limit$nov.std)~s(c(1:length(DATN_predict.m.mean.limit$nov.std))),sp=.1)$fitted.values[c(TRUE,FALSE)]/max(OVdata.spatiotemp.ribbon.limit$x)-0.2,1),
                         guide=FALSE)+
  # geom_smooth(data=OVdata.spatiotemp.ribbon.limit,aes(x=x,y=y,colour=..x..),size=2,method="loess",span=.1,se=FALSE)+
  # scale_colour_gradientn(colours=c(viridis_pal()(10)),
  #                        values=c(0,gam(sort(DATN_predict.m.mean.limit$nov.std)~s(c(1:length(DATN_predict.m.mean.limit$nov.std))),sp=.1)$fitted.values[c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)]/max(OVdata.spatiotemp.ribbon.limit$x)-0.1,1),
  #                        guide=FALSE)+
  new_scale_colour()+
  scale_fill_manual(values=c("darkgrey"))+
  geom_segment(data=OVdata.spatiotemp.ribbon[x>0],aes(x=x,xend=x,y=0,yend=y/1.5,colour=(x)),size=.2,alpha=.16)+
  scale_colour_gradientn(colours=c(viridis_pal()(72)),
                         values=c(0,gam(sort(DATN_predict.m.mean$nov.std)~s(c(1:length(DATN_predict.m.mean$nov.std))),sp=.1)$fitted.values[c(TRUE,FALSE)]/max(OVdata.spatiotemp.ribbon$x),1),
                         guide=FALSE)+
  scale_x_continuous(#trans="log10",
    limits=c(0,1),
    expand=c(0,0),
    breaks=c(0,0.2,0.4,0.6,0.8,1),
    labels = scaleFUN0)+
  annotation_custom(g, xmin=.17, xmax=.23, ymin=11, ymax=13)+
  annotation_custom(g, xmin=.17, xmax=.23, ymin=9.2, ymax=11.2)+
  geom_point(aes(x=.2,y=10.2), fill="white", colour="white",shape=15,size=10,alpha=0.6)+
  annotate(geom="text",x=.42,y=11,label="    Change over time (n = 139)
           
                                 Change over time with habitat limits (n = 139)",colour="black",size=3)+
  scale_y_continuous(breaks=c(0,4,8,12),expand=c(0,0),limits=c(0,14))+
  theme_bw()+ 
  theme(text = element_text(size=10),
        legend.position = c(0.456, 0.60),
        legend.margin=margin(c(0,0,0,0)))+
  labs(x="Thermal habitat (% non-overlap)",
       y="Kernel Density Estimation",
       fill=NULL)+
  coord_fixed(ratio = .05)

ab<-plot_grid(a,b,ncol=1,rel_heights = c(2,1.5))
ab

#Fig 4: Mean ED across different sensitivities####
OVdata.lakesplits.interp<-fread("D:/DATN/DATN_OVdata.lakesplits.interp.std.time_v51.csv")
OVdata.splits.interp<-fread("D:/DATN/DATN_OVdata.splits.interp.std.time_v51.csv")

#add separate panels for individual lakes
OVdata.lakesplits.interp[lake=="Pesiojarvi"]$lake<-"Pesiöjärvi"
OVdata.lakesplits.interp_sub<-(OVdata.lakesplits.interp[lake=="Annie"|
                                                         lake=="Muggelsee"|
                                                         lake=="Pesiöjärvi"|#Pesiojarvi
                                                         lake=="Zurich"|
                                                         lake=="Tanganyika"|
                                                         lake=="Tahoe"])
OVdata.lakesplits.interp_sub$lake<-factor(OVdata.lakesplits.interp_sub$lake,levels=c("Pesiöjärvi","Muggelsee","Annie","Tahoe","Zurich","Tanganyika"))

a<-ggplot()+
  geom_raster(data=OVdata.splits.interp,
              aes(x=smd,
                  y=zmd,
                  fill=nov),
              interpolate=T)+
  scale_fill_viridis(trans="log",breaks=c(6,8,10,12,14,16,18)
    )+
  geom_contour(data=OVdata.splits.interp,
               aes(x=smd,
                   y=zmd,
                   z=nov),
               colour="white",
               breaks=c(6,7,8,9,10,11,12,13,14,15,16,17,18))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  labs(x="Seasonal habitat restriction (proportion of maximum)",
       y="Depth habitat restriction (proportion of maximum)",
       fill="Thermal habitat change
(% non-overlap)")

b<-ggplot()+
  geom_raster(data=OVdata.lakesplits.interp_sub,
              aes(x=smd,
                  y=zmd,
                  fill=nov),
              interpolate=T)+
  scale_fill_viridis(trans="log",breaks=seq(0,30,5))+
  geom_contour(data=OVdata.lakesplits.interp_sub,
               aes(x=smd,
                   y=zmd,
                   z=nov),
               colour="white",
               breaks=seq(0,30,1))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  facet_wrap(.~lake)+
#   geom_point(data = OVdata.lakesplits.interp_sub[lake=="Zurich"],aes(x=0.58,y=0.89),colour="black",size=2)+
#   geom_text(data = OVdata.lakesplits.interp_sub[lake=="Zurich"],aes(x=0.68,y=0.69,label="Planktothrix
# rubescens"),colour="black",size=4,fontface="italic")+
  # annotate(geom="text")+
  theme_bw()+
  theme(panel.spacing.x = unit(7, "mm"))+
  labs(x="Seasonal habitat restriction (proportion of maximum)",
       y="Depth habitat restriction (proportion of maximum)",
       fill="Thermal habitat change
(% non-overlap)")

plot_grid(a,b,ncol=1,rel_heights = c(2,1.7))

#Same plot but now for every lake compiled into a pdf.
OVdata.lakesplits.interp<-fread("D:/DATN/DATN_OVdata.lakesplits.interp.std.time_v51.csv",stringsAsFactors = T)
lakenames<-unique(OVdata.lakesplits.interp$lake)
lakenames<-lakenames[order(lakenames)]

pdf("D:/DATN/Supplement/S4/Shiftplots_AllLakes_v12.pdf", onefile = TRUE)
for (lakename in lakenames){
print(ggplot()+
  geom_raster(data=OVdata.lakesplits.interp[lake==lakename],
              aes(x=abs(smd),
                  y=abs(zmd),
                  fill=nov),
              interpolate=T)+
  scale_fill_viridis_c(#limits=c(0,45)#trans="log",breaks=c(.001,.002,.004,.008,.016,.032,.064,.128,.256)
    )+
  geom_contour(data=OVdata.lakesplits.interp[lake==lakename],
               aes(x=abs(smd),
                   y=abs(zmd),
                   z=nov),
               colour="white")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  labs(x="Seasonal habitat restriction (proportion of maximum)",
       y="Depth habitat restriction (proportion of maximum)",
       fill="Thermal habitat change 
(% non-overlap)")+
    ggtitle(lakename))
}
dev.off()

##Now printed to multiple jpegs
for (lakename in lakenames){
  jpeg(paste0("D:/DATN/Supplement/S4/S4_Shiftplot_",lakename,"_v12.jpg"))
  print(ggplot()+
          geom_raster(data=OVdata.lakesplits.interp[lake==lakename],
                      aes(x=abs(smd),
                          y=abs(zmd),
                          fill=nov),
                      interpolate=T)+
          scale_fill_viridis_c(#limits=c(0,45)#trans="log",breaks=c(.001,.002,.004,.008,.016,.032,.064,.128,.256)
          )+
          geom_contour(data=OVdata.lakesplits.interp[lake==lakename],
                       aes(x=abs(smd),
                           y=abs(zmd),
                           z=nov),
                       colour="white")+
          scale_y_continuous(expand=c(0,0))+
          scale_x_continuous(expand=c(0,0))+
          theme_bw()+
          labs(x="Seasonal habitat restriction (proportion of maximum)",
               y="Depth habitat restriction (proportion of maximum)",
               fill="Thermal habitat change 
(% non-overlap)")+
          ggtitle(lakename))
  dev.off()
}



####FIG 5: Effects on thermal nonoverlap from lake characteristics across gradients of depth and doy limitations####
#DATN_predict_new<-fread("D:/DATN/DATN_predict_new.means_v47.csv",header=TRUE,stringsAsFactors = TRUE)
OVdata.lakesplits.interp<-fread("D:/DATN/DATN_OVdata.lakesplits.interp.std.time_v51.csv")
DATN_predict.fit<-fread("D:/DATN/DATN_predict.fit_v43.csv",header=TRUE,stringsAsFactors = TRUE)
names(DATN_predict.fit)
lakeinfo<-unique(DATN_predict.fit[,c(1,53,57)])
OVdata.lakesplits.interp.m<-merge(OVdata.lakesplits.interp,lakeinfo)
summary(OVdata.lakesplits.interp.m)
head(OVdata.lakesplits.interp.m)
OVdata.lakesplits.interp.m<-OVdata.lakesplits.interp.m[(((smd*100)/5)%%1==0)==TRUE]
OVdata.lakesplits.interp.m<-OVdata.lakesplits.interp.m[(((zmd*100)/5)%%1==0)==TRUE]


DATN_predict_std.time<-fread("D:/DATN/DATA_predict.m.mean_v51.csv",stringsAsFactors = T) 
summary(DATN_predict_std.time)
DATN_predict_std.time.means<-DATN_predict_std.time[,.(nov=mean(nov),
                                                      nov.randnull=mean(nov.randnull)),
                                                   .(lake,zmd,smd,abslat,logMeanDepth)]

DATN_predict_new.means<-fread("D:/DATN/DATA_DATN_predict_new.means_allreps_v51.csv",stringsAsFactors = T)
DATN_predict_new.means$smd<-1-(1/(DATN_predict_new.means$smd))
DATN_predict_new.means$zmd<-1-(1/(DATN_predict_new.means$zmd))
summary(DATN_predict_new.means)
DATN_predict<-fread("D:/DATN/DATA_predict_v50.csv",header=TRUE,stringsAsFactors = TRUE)
lakeinfo<-unique(DATN_predict[,c(1,23,51)])

a<-ggplot()+
  geom_density(data=(lakeinfo),aes(logMeanDepth),colour="darkgrey",fill="darkgrey",adjust=.5)+
  scale_x_continuous(trans="log",breaks=c(1,5,25,100,500),expand=c(0,0))+
  theme_void()

b<-ggplot()+
  geom_density(data=(lakeinfo),aes(abslat),colour="darkgrey",fill="darkgrey",adjust=.5)+
  #scale_x_continuous(trans="log",breaks=c(5,10,15,20,25),expand=c(0,0))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60),expand=c(0,0))+
  theme_void()

c<-ggplot()+
  geom_smooth(data=OVdata.lakesplits.interp.m[zmd==0],
              aes(y=nov,
                  x=exp(logMeanDepth),
                  colour=smd,
                  group=smd),se=FALSE,method="loess",span=1,size=.5)+
  #scale_colour_gradientn(colors=rev(c("#005a32","#238443","#41ab5d","#78c679","#addd8e","#d9f0a3")))+
  labs(x="Mean depth (m)",
       y="Thermal habitat change (% non-overlap)",
       colour="Seasonal habitat restriction 
(proportion of maximum)")+
  scale_x_continuous(trans="log",breaks=c(1,5,25,100,500),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0,23))+
  scale_colour_gradientn(colours=rev(c("#005a32","#238443","#41ab5d","#78c679","#addd8e","#d9f0a3")))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  guides(colour=FALSE)

d<-ggplot()+
  geom_smooth(data=OVdata.lakesplits.interp.m[zmd==0],
              aes(y=nov,
                  x=abslat,
                  colour=smd,
                  group=smd),se=FALSE,method="loess",span=1,size=.5)+
  labs(x="Latitude (°N/S)",
       y="Thermal habitat change (% non-overlap)",
       colour="Seasonal habitat restriction 
(proportion of maximum)")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0,23))+
  scale_colour_gradientn(colours=rev(c("#005a32","#238443","#41ab5d","#78c679","#addd8e","#d9f0a3")),
                         breaks=c(seq(0,95,5)/100),
                         guide = "colourbar")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")+
  guides(colour = guide_legend(reverse = TRUE))

e<-ggplot()+
  geom_smooth(data=OVdata.lakesplits.interp.m[smd==0],
              aes(y=nov,
                  x=exp(logMeanDepth),
                  colour=zmd,
                  group=zmd),se=FALSE,method="loess",span=1,size=.5)+
  labs(x="Mean depth (m)",
       y="Thermal habitat change (% non-overlap)",
       colour="Depth habitat restriction 
(proportion of maximum)")+
  scale_x_continuous(trans="log",breaks=c(1,5,25,100,500),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0,23))+
  scale_colour_gradientn(colors=rev(c("#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0")))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  guides(colour=FALSE)

f<-ggplot()+
  geom_smooth(data=OVdata.lakesplits.interp.m[smd==0],
              aes(y=nov,
                  x=abslat,
                  colour=zmd,
                  group=zmd),se=FALSE,method="loess",span=1,size=.5)+
  labs(x="Latitude (°N/S)",
       y="Thermal habitat change (% non-overlap)",
       colour="Depth habitat restriction 
(proportion of maximum)")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0,23))+
  scale_colour_gradientn(colours=rev(c("#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0")),
                         breaks=c(seq(0,95,5)/100))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")+
  guides(colour = guide_legend(reverse = TRUE))

h.data<-data.table(id=1,level=c(1:length(seq(0,95,5)))/100,smd=seq(0,95,5))
h<-as_ggplot(get_legend(ggplot()+
  geom_raster(data=h.data,aes(x=id,y=level,fill=smd/100))+
  scale_fill_gradientn(colours=rev(c("#005a32","#238443","#41ab5d","#78c679","#addd8e","#d9f0a3")),
                         breaks=c(0,0.25,0.5,0.75))+
    labs(fill="Seasonal habitat restriction 
(proportion of maximum)")))
i<-as_ggplot(get_legend(ggplot()+
                          geom_raster(data=h.data,aes(x=id,y=level,fill=smd/100))+
                          scale_fill_gradientn(colours=rev(c("#7a0177","#ae017e","#dd3497","#f768a1","#fa9fb5","#fcc5c0")),
                                               breaks=c(0,0.25,0.5,0.75))+
                          labs(fill="Depth habitat restriction 
(proportion of maximum)")))
  
ace<-plot_grid(a,c,e,ncol=1,rel_heights = c(1,3,3),align="v", axis = "lr")
bdf<-plot_grid(b,d,f,ncol=1,rel_heights = c(1,3,3),align="v", axis = "lr")
hi<-plot_grid(NULL,h,i,ncol=1,rel_heights = c(1,3,3))
plot_grid(ace,bdf,hi,ncol=3,rel_widths = c(1,1,.5))


####FIGURE S#: lake warming rate/non-overlap comparison.
DATN_predict.m.mean<-fread("D:/DATN/DATA_predict.m.mean_v51.csv")
DATN_predict.m.mean<-DATN_predict.m.mean[smd==1&zmd==1,.(nov=mean(nov),
                                                         nov.std=mean(nov.std.predict.resid),
                                                         nov.randnull=mean(nov.randnull),
                                                         nov.std.raw=mean(nov-nov.randnull),
                                                         lat=mean(lat),
                                                         lon=mean(lon)),.(lake)]
#EDdata.lakemeans<-fread("D:/DATN/DATN_EDdata.lakemeans_v48.csv",stringsAsFactors = T)
DATN_predict.m.mean$nov<-DATN_predict.m.mean$nov/100
DATN_predict.m.mean$nov.std<-DATN_predict.m.mean$nov.std/100
DATN_predict.m.mean$nov.std.raw<-DATN_predict.m.mean$nov.std.raw/100
DATN_predict.m.mean[nov.std.raw<0]$nov.std.raw<-0
DATN_predict.m.mean$nov.randnull<-DATN_predict.m.mean$nov.randnull/100
summary(DATN_predict.m.mean)

#EDdata.lakemeans<-fread("D:/DATN/DATN_EDdata.lakemeans_v48.csv")
data<-fread("D:/DATN/DATA_OVtemp_v1.csv")
data<-merge(DATN_predict.m.mean,data)
head(data)
summary(data)

cor(data$nov.std,data$tempmb_slope,method="kendall")
cor.test(data$nov.std,data$tempmb_slope,method="kendall")

cor(data$nov.std,abs(data$tempmb_slope),method="kendall")
cor.test(data$nov.std,abs(data$tempmb_slope),method="kendall")

a<-ggplot()+
  geom_density(data=data,aes(x=tempmb_slope),fill="darkgrey",alpha=0.6)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line())+
  geom_vline(xintercept=0,linetype="dashed")+
  #geom_hline(yintercept=0)+
  labs(x=expression(paste("Warming rate (°C ",year^-1,")")),
       y="Kernel density estimate")
  
b<-ggplot()+
  geom_point(data=data,aes(y=nov.std,x=tempmb_slope))+
  #scale_y_continuous(trans="log",breaks=c(0.01,0.02,0.04,0.08,0.16,0.32,0.64),limits=c(0.005,1))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line())+
  geom_vline(xintercept=0,linetype="dashed")+
  #geom_hline(yintercept=0)+
  labs(#x="Warming rate (°C year-1)",
       x=expression(paste("Warming rate (°C ",year^-1,")")),
       y="Thermal non-overlap")

plot_grid(a,b,ncol=1,rel_heights = c(1,2),align="v", axis = "lr")


expression(paste("Warming rate (°C ",year^-1))


####Fig 1####
doydepth.grid.comb<-fread("D:/DATN/DATA_doydepth.grid.comb_Fig1_v1.csv",stringsAsFactors = T)

# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")
pal2<-GetColors(n = 10, scheme = "smooth rainbow", alpha = NULL, 
             start = .15, end = 1, bias = 1.6, reverse = FALSE, blind = NULL, 
             gray = FALSE)
pal4<-rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFF8A9", "#E6F598",
        "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"))
pal5<-rev(c("#630000","#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFF8A9", "#E6F598",
            "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2","#644089","#441450"))

#establish the min and max of scale 
grandmin <- floor(min(doydepth.grid.comb$temp))
grandmax <- ceiling(max(doydepth.grid.comb$temp))

#define the number of breaks.  In this case 8 +1 
mybreaks <- seq(grandmin, grandmax, length.out = grandmax-grandmin+1)

#Function to return the dersired number of colors
mycolors<- function(x) {
  colors<-colorRampPalette(rev(c("#630000","#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFF8A9", "#E6F598",
                                 "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2","#644089","#441450")))( grandmax-grandmin+1 )
  colors[1:x]
}

#Function to create labels for legend
breaklabel <- function(x){
  labels<- paste0(mybreaks[1:16], "-", mybreaks[2:17])
  labels[1:x]
}


#overall
1-overlap(list(doydepth.grid.comb[date=="Past"]$temp,
                     doydepth.grid.comb[date=="Present"]$temp),nbins=100,plot=T)$OV

#sp1
1-overlap(list(doydepth.grid.comb[date=="Past"&depth>=25&depth<=40]$temp,
                     doydepth.grid.comb[date=="Present"&depth>=25&depth<=40]$temp),nbins=100,plot=T)$OV

#sp2
1-overlap(list(doydepth.grid.comb[date=="Past"&doy>=130&doy<=203]$temp,
             doydepth.grid.comb[date=="Present"&doy>=130&doy<=203]$temp),nbins=100,plot=T)$OV

#sp3
1-overlap(list(doydepth.grid.comb[date=="Past"&doy>=220&doy<=293&depth>=35&depth<=50]$temp,
             doydepth.grid.comb[date=="Present"&doy>=220&doy<=293&depth>=35&depth<=50]$temp),nbins=100,plot=T)$OV


a<-ggplot()+
  #geom_raster(data=doydepth.grid,aes(y=depth,x=doy,fill=temp),interpolate=TRUE)+
  geom_contour_filled(data=doydepth.grid.comb[date=="Past"],aes(y=depth,x=doy,z=temp),breaks= mybreaks, show.legend = TRUE)+
  scale_fill_manual(palette=mycolors, values=breaklabel(grandmax-grandmin), name="Value", drop=FALSE) +
  #geom_contour(data=doydepth.grid,aes(y=depth,x=doy,z=temp))+
  #scale_fill_viridis(option = "inferno",limits=c(4,22))+
  # scale_fill_gradientn(colours = pal5,limits=c(min(doydepth.grid.comb$temp),
  #                                              max(doydepth.grid.comb$temp))) +
  # scale_colour_discrete(colours = pal5,limits=c(min(doydepth.grid.comb$temp),
  #                                              max(doydepth.grid.comb$temp))) +
  scale_y_reverse(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  ggtitle("Past")+
  annotate("rect",xmin = 1, xmax = 365, ymin = 25, ymax = 40,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("rect",xmin = 130, xmax = 203, ymin = 0, ymax = 60,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("rect",xmin = 220, xmax = 293, ymin = 35, ymax = 50,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("text",x=40,y=27,label="Species 1",colour="white")+
  annotate("text",x=170,y=58,label="Species 2",colour="white")+
  annotate("text",x=256,y=48,label="Species 3",colour="white")+
  guides(fill=FALSE)+
  labs(x="Day of the year",y="Depth (m)")
  

b<-ggplot()+
  #geom_raster(data=doydepth.grid.present,aes(y=depth,x=doy,fill=temp),interpolate=TRUE)+
  geom_contour_filled(data=doydepth.grid.comb[date=="Present"],aes(y=depth,x=doy,z=temp),breaks= mybreaks, show.legend = TRUE)+
  scale_fill_manual(palette=mycolors, values=breaklabel(grandmax-grandmin), name="Temperature 
(°C)", drop=FALSE) +
  #geom_contour(data=doydepth.grid,aes(y=depth,x=doy,z=temp))+
  #scale_fill_viridis(option = "inferno",limits=c(4,22))+
  # scale_fill_gradientn(colours = pal5,limits=c(min(doydepth.grid.comb$temp),
  #                                              max(doydepth.grid.comb$temp))) +
  scale_colour_discrete() +
  scale_y_reverse(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  ggtitle("Recent")+
  annotate("rect",xmin = 1, xmax = 365, ymin = 25, ymax = 40,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("rect",xmin = 130, xmax = 203, ymin = 0, ymax = 60,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("rect",xmin = 220, xmax = 293, ymin = 35, ymax = 50,fill=NA,colour="white",size=1,linetype="dashed")+
  annotate("text",x=40,y=27,label="Species 1",colour="white")+
  annotate("text",x=170,y=58,label="Species 2",colour="white")+
  annotate("text",x=256,y=48,label="Species 3",colour="white")+
  labs(fill="Temperature (°C)",x="Day of the year",y="Depth (m)")

c<-ggplot()+
  geom_point(data=doydepth.grid.comb[,.(temp=density(temp,
                                                     n=((grandmax-grandmin)*10+1),
                                                     from=grandmin,
                                                     to=grandmax,
                                                     adjust=1.8)$x,
                                        temp.vol=density(temp,
                                                         n=((grandmax-grandmin)*10+1),
                                                         from=grandmin,
                                                         to=grandmax,
                                                         adjust=1.8)$y),.(date)],
             aes(x=temp,y=(60*temp.vol*(sum(density(doydepth.grid.comb$temp,
                                                    n=((grandmax-grandmin)*10+1),
                                                    from=grandmin,
                                                    to=grandmax,
                                                    adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                 colour=date,fill=date),
             size=1.2,
             alpha=.5)+
  geom_ribbon(data=doydepth.grid.comb[,.(temp=density(temp,
                                                      n=((grandmax-grandmin)*10+1),
                                                      from=grandmin,
                                                      to=grandmax,
                                                      adjust=1.8)$x,
                                         temp.vol=density(temp,
                                                          n=((grandmax-grandmin)*10+1),
                                                          from=grandmin,
                                                          to=grandmax,
                                                          adjust=1.8)$y),.(date)],
              aes(x=temp,ymin=0,ymax=(60*temp.vol*(sum(density(doydepth.grid.comb$temp,
                                                           n=((grandmax-grandmin)*10+1),
                                                           from=grandmin,
                                                           to=grandmax,
                                                           adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                  colour=date,fill=date),
               alpha=.5)+
  ggtitle("Overall")+
  scale_x_continuous(limits=c(grandmin,grandmax),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"),guide = guide_legend(
    direction = "horizontal",
    title.position = "top"
  ))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  guides(colour=FALSE)+
  theme(legend.position = c(0.62,0.35))+
  labs(x="Temperature (°C)",y=expression("Average daily volume ("~m^3*")"),fill="Time period")+
  annotate("text",x=15,y=(.65),label="Thermal non-overlap = 34%
Warming = +1.14 °C")

d<-ggplot()+
  geom_point(data=doydepth.grid.comb[depth>=25&depth<=40,.(temp=density(temp,
                                                     n=((grandmax-grandmin)*10+1),
                                                     from=grandmin,
                                                     to=grandmax,
                                                     adjust=1.8)$x,
                                        temp.vol=density(temp,
                                                         n=((grandmax-grandmin)*10+1),
                                                         from=grandmin,
                                                         to=grandmax,
                                                         adjust=1.8)$y),.(date)],
             aes(x=temp,y=(15*temp.vol*(sum(density(doydepth.grid.comb[depth>=25&depth<=40]$temp,
                                                    n=((grandmax-grandmin)*10+1),
                                                    from=grandmin,
                                                    to=grandmax,
                                                    adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                 colour=date,fill=date),
             size=1.2,
             alpha=.5)+
  geom_ribbon(data=doydepth.grid.comb[depth>=25&depth<=40,.(temp=density(temp,
                                                      n=((grandmax-grandmin)*10+1),
                                                      from=grandmin,
                                                      to=grandmax,
                                                      adjust=1.8)$x,
                                         temp.vol=density(temp,
                                                          n=((grandmax-grandmin)*10+1),
                                                          from=grandmin,
                                                          to=grandmax,
                                                          adjust=1.8)$y),.(date)],
              aes(x=temp,ymin=0,ymax=(15*temp.vol*(sum(density(doydepth.grid.comb[depth>=25&depth<=40]$temp,
                                                               n=((grandmax-grandmin)*10+1),
                                                               from=grandmin,
                                                               to=grandmax,
                                                               adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                  colour=date,fill=date),
              alpha=.5)+
  ggtitle("Species 1")+
  scale_x_continuous(#limits=c(grandmin,grandmax),
    limits=c(grandmin,17),
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  guides(colour=FALSE,fill=FALSE)+
  labs(x="Temperature (°C)",y=expression("Average daily volume ("~m^3*")"))+
  annotate("text",x=14,y=(.15),label="Thermal non-overlap = 35%
Warming = +1.06 °C")

e<-ggplot()+
  geom_point(data=doydepth.grid.comb[doy>=130&doy<=203,.(temp=density(temp,
                                                                        n=((grandmax-grandmin)*10+1),
                                                                        from=grandmin,
                                                                        to=grandmax,
                                                                        adjust=1.8)$x,
                                                           temp.vol=density(temp,
                                                                            n=((grandmax-grandmin)*10+1),
                                                                            from=grandmin,
                                                                            to=grandmax,
                                                                            adjust=1.8)$y),.(date)],
             aes(x=temp,y=(15*temp.vol*(sum(density(doydepth.grid.comb[doy>=130&doy<=203]$temp,
                                                    n=((grandmax-grandmin)*10+1),
                                                    from=grandmin,
                                                    to=grandmax,
                                                    adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                 colour=date,fill=date),
             size=1.2,
             alpha=.5)+
  geom_ribbon(data=doydepth.grid.comb[doy>=130&doy<=203,.(temp=density(temp,
                                                                         n=((grandmax-grandmin)*10+1),
                                                                         from=grandmin,
                                                                         to=grandmax,
                                                                         adjust=1.8)$x,
                                                            temp.vol=density(temp,
                                                                             n=((grandmax-grandmin)*10+1),
                                                                             from=grandmin,
                                                                             to=grandmax,
                                                                             adjust=1.8)$y),.(date)],
              aes(x=temp,ymin=0,ymax=(15*temp.vol*(sum(density(doydepth.grid.comb[doy>=130&doy<=203]$temp,
                                                               n=((grandmax-grandmin)*10+1),
                                                               from=grandmin,
                                                               to=grandmax,
                                                               adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                  colour=date,fill=date),
              alpha=.5)+
  ggtitle("Species 2")+
  scale_x_continuous(#limits=c(grandmin,grandmax),
    expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  guides(colour=FALSE,fill=FALSE)+
  labs(x="Temperature (°C)",y=expression("Average daily volume ("~m^3*")"))+
  annotate("text",x=15,y=(.15),label="Thermal non-overlap = 47%
Warming = +2.13 °C")

f<-ggplot()+
  geom_point(data=doydepth.grid.comb[doy>=220&doy<=293&depth>=35&depth<=50,.(temp=density(temp,
                                                                      n=((grandmax-grandmin)*10+1),
                                                                      from=grandmin,
                                                                      to=grandmax,
                                                                      adjust=1.8)$x,
                                                         temp.vol=density(temp,
                                                                          n=((grandmax-grandmin)*10+1),
                                                                          from=grandmin,
                                                                          to=grandmax,
                                                                          adjust=1.8)$y),.(date)],
             aes(x=temp,y=(2.4*temp.vol*(sum(density(doydepth.grid.comb[doy>=220&doy<=293&depth>=35&depth<=50]$temp,
                                                    n=((grandmax-grandmin)*10+1),
                                                    from=grandmin,
                                                    to=grandmax,
                                                    adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                 colour=date,fill=date),
             size=1.2,
             alpha=.5)+
  geom_ribbon(data=doydepth.grid.comb[doy>=220&doy<=293&depth>=35&depth<=50,.(temp=density(temp,
                                                                       n=((grandmax-grandmin)*10+1),
                                                                       from=grandmin,
                                                                       to=grandmax,
                                                                       adjust=1.8)$x,
                                                          temp.vol=density(temp,
                                                                           n=((grandmax-grandmin)*10+1),
                                                                           from=grandmin,
                                                                           to=grandmax,
                                                                           adjust=1.8)$y),.(date)],
              aes(x=temp,ymin=0,ymax=(2.4*temp.vol*(sum(density(doydepth.grid.comb[doy>=220&doy<=293&depth>=35&depth<=50]$temp,
                                                               n=((grandmax-grandmin)*10+1),
                                                               from=grandmin,
                                                               to=grandmax,
                                                               adjust=1.8)$y)/(((grandmax-grandmin)*10)+2))),
                  colour=date,fill=date),
              alpha=.5)+
  ggtitle("Species 3")+
  scale_x_continuous(#limits=c(grandmin,grandmax),
    limits=c(6,13),
    expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  guides(colour=FALSE,fill=FALSE)+
  labs(x="Temperature (°C)",y=expression("Average daily volume ("~m^3*")"))+
  annotate("text",x=11.4,y=(.05),label="Thermal non-overlap = 65%
Warming = +1.14 °C")



ab<-plot_grid(a,b,ncol=2,rel_widths = c(1,1.19))
cdef<-plot_grid(c,d,e,f,ncol=2)
abcdef<-plot_grid(ab,cdef,ncol=1)
  abcdef


####Partial dependency plots for full BRT####
BRT.partdep<-fread("D:/DATN/DATA_DATN_BRT.partdep_v51.csv",header=T,stringsAsFactors = T)
BRT.contrib<-fread("D:/DATN/DATA_DATN_BRT.contrib_v51.csv",header=T,stringsAsFactors = T)


BRT.contrib$rep<-rep(c(1:13),each=100)
BRT.contrib[,':='(rel.inf.rs=100*rel.inf/(sum(BRT.contrib$rel.inf)/100)),.(rep)]

ggplot()+
  geom_violin(data=BRT.contrib,
              aes(y=fct_reorder(var,rel.inf.rs,.fun=mean),
                  x=rel.inf.rs,colour=var),
              alpha=0.4)+
  geom_jitter(data=BRT.contrib,
              aes(y=fct_reorder(var,rel.inf.rs,.fun=mean),
                  x=rel.inf.rs,
                  colour=var),
              height = 0.3, 
              alpha=0.2,
              #colour=NA,
              shape=16,)+
  geom_boxplot(data=BRT.contrib,
                aes(y=fct_reorder(var,rel.inf.rs,.fun=mean),
                    x=rel.inf.rs),
                colour="black",
                fill=NA,
                width = 0.15,
               size=0.5,
               outlier.shape=NA,
               coef=0)+
  scale_y_discrete(labels=rev(c("Minimum decimal year",
                                "Maximum decimal year",
                                "Latitude (° N/S)",
                                "Log(mean depth (m))",
                                "Log(Mean samples per year 
(days/year))",
                                "Trend in length of seasonal 
coverage (days/year)",
                                "Depth shift limit 
(proportion of maximum)",
                                "Maximum day of the year",
                                "Minimum day of the year",
                                "Timeseries delineation",
                                "Seasonal shift limit 
(proportion of maximum)",
                                "Trend in mean samples 
per year (days/year)")))+
  theme_bw()+
  labs(x="Relative influence",
       y="Predictor variable")+
  guides(colour=FALSE)

BRT.partdep<-merge(BRT.partdep,BRT.contrib)
BRT.partdep<-droplevels(BRT.partdep)
BRT.partdep[var=="smd"]$var.value<-abs(1-(1/BRT.partdep[var=="smd"]$var.value))
BRT.partdep[var=="zmd"]$var.value<-abs(1-(1/BRT.partdep[var=="zmd"]$var.value))

levels(BRT.partdep$var)<-c("Latitude (° N/S)",
                           "Maximum year",
                           "Minimum year",
                           "Maximum day of the year",
                           "Minimum day of the year",
                           "Trend in length of seasonal 
coverage (days/year)",
                           "Trend in mean samples 
per year (days/year)",
                           "Log(mean depth (m))",
                           "Log(Mean samples per year 
(days/year))",
                           "Seasonal shift limit 
(proportion of maximum)",
                           "Timeseries delineation",
                           "Depth shift limit 
(proportion of maximum)")

ggplot()+
  geom_line(data=BRT.partdep,aes(y=nov,x=var.value,group=j,colour=var),alpha=0.3)+
  #geom_smooth(data=BRT.partdep,aes(y=nov,x=var.value),colour="black",method="gam")+
  geom_smooth(data=BRT.partdep,aes(y=nov,x=var.value),colour="black",method="loess")+
  facet_wrap(~rev(fct_reorder(var,rel.inf.rs,.fun=mean)),scales="free",)+
  facet_wrap(~fct_reorder(var,rel.inf.rs,.fun=mean,.desc=TRUE),scales="free")+
  #facet_wrap(~var,scales="free_x")+
  theme_bw()+
  labs(y="Thermal non-overlap",
       x="Predictor variable value")+
  guides(colour=FALSE)



####Fig 6#### Example for planktothrix rubescense in zurich####
lakedata.full<-fread("D:/DATN/RefitData_v6/Zurich.csv")
lakedata.full<-lakedata.full[is.na(fit)==FALSE]
summary(lakedata.full)


yearnumber<-length(unique(lakedata.full$year))
years<-(unique(lakedata.full$year))
years<-years[order(years)]
years.i.1<-sort(sample(years,round(length(years)/2,0),replace=F))
years.i.2<-sort(years[years%in%years.i.1==FALSE])

i<-round(quantile(c(1:length(years)),probs=c(0.5)))

lakedata.full$Period<-factor()
lakedata.full[decyear<years[i]]$Period<-"1936-1975"
lakedata.full[decyear>=years[i]]$Period<-"1976-2015"

####For full lake zurich
(1-overlap(list(lakedata.full[decyear<years[i]]$fit,
               lakedata.full[decyear>=years[i]]$fit),nbins=100,plot=T)$OV)-
  (1-overlap(list(lakedata.full[year%in%years.i.1]$fit,
               lakedata.full[year%in%years.i.2]$fit),nbins=100,plot=T)$OV)

0.1966404-0.04051618
#15.6%

####For p rub habitat in lake zurich
lakedata<-lakedata.full[depth>5&
                          depth<20&
                          doy<yday("2000-09-30")&
                          doy>yday("2000-07-01")&
                          is.na(fit)==FALSE]


(1-overlap(list(lakedata[decyear<years[i]]$fit,
               lakedata[decyear>=years[i]]$fit),nbins=100)$OV)-
  (1-overlap(list(lakedata[year%in%years.i.1]$fit,
               lakedata[year%in%years.i.2]$fit),nbins=100)$OV)

lakedata$Period<-factor()
lakedata[decyear<years[i]]$Period<-"1936-1975"
lakedata[decyear>=years[i]]$Period<-"1976-2015"
lakedata[,.(fit=mean(fit)),.(Period)]
lakedata.full[,.(fit=mean(fit)),.(Period)]

a<-ggplot()+
  geom_density(data=lakedata.full,aes(x=fit,fill=Period,colour=Period),adjust=2,alpha=0.3)+
  scale_x_continuous(limits=c(2,25),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  annotate("text",x=17,y=(.18),label="  Thermal habitat change = 
16% non-overlap

Warming = +0.51 °C")+
  theme(legend.position = c(0.65,0.35))+
  labs(x="Temperature (°C)",y="Kernel density estimation",fill="Time period")+
  guides(colour=FALSE)

b<-ggplot()+
  geom_density(data=lakedata,aes(x=fit,fill=Period,colour=Period),adjust=2,alpha=0.3)+
  scale_x_continuous(limits=c(2,25),expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.14),expand=c(0,0))+
  scale_fill_manual(values=c("darkgrey","black"))+
  scale_colour_manual(values=c("darkgrey","black"))+
  theme_bw()+
  annotate("text",x=14,y=(.12),label="  Thermal habitat change = 
22% non-overlap
Warming = +1.02 °C")+
  theme(legend.position = c(0.85,0.5))+
  labs(x="Temperature (°C)",y="Kernel density estimation",fill="Time period")+
  guides(colour=FALSE,
         fill=FALSE)

plot_grid(a,b,ncol=2)

base<-c(1:100)/100
halved<-(base*0.5)/(1-base+(base*0.5))

#Simulations for Extended data Fig 1
# Create different ordered samples of the population to produce random pairings
ranges1 <- expand.grid(sd = seq(from=-3, to=3, length.out=100),
                      mean = seq(from=-3, to=3, length.out=100))
ranges2 <- expand.grid(sd = seq(from=-3, to=3, length.out=100),
                       mean = seq(from=-3, to=3, length.out=100))
ranges3 <- expand.grid(sd = seq(from=-3, to=3, length.out=100),
                       mean = seq(from=-3, to=3, length.out=100))

#library(overlapping)
# define output vectors
thermal_overlap1 <- c()
thermal_overlap2 <- c()
thermal_overlap3 <- c()
# loop over all potential values of mean and sd of temperature variation. Using a sample size
#of 10,000 for speed.
#i<-1
for(i in 1:nrow(ranges1)){
  if(i==1|i==1000|i==2000|i==3000|i==4000|i==5000|i==6000|i==7000|i==8000|i==9000){print(i)
    Sys.time()}
  #Baseline 1
  # focal time frame
  x_var <- rnorm(100000, mean = ranges1[i,'mean']+15, sd = ranges1[i,'sd']+4)
  # baseline time frame
  y_var <- rnorm(100000, mean = 15, sd = 4)
  # baseline time frame2
  y_var2 <- rnorm(100000, mean = 15, sd = 4)
  
  # estimate the overlap of the two distributions
  thermal_overlap1[i] <- 1-(overlap(list(x_var, y_var))$OV+(1-(overlap(list(y_var2, y_var))$OV)))
  
  #Baseline 2
  # focal time frame
  x_var <- rnorm(100000, mean = ranges2[i,'mean']+15, sd = ranges2[i,'sd']+6)
  # baseline time frame
  y_var <- rnorm(100000, mean = 15, sd = 6)
  # baseline time frame2
  y_var2 <- rnorm(100000, mean = 15, sd = 6)
  
  # estimate the overlap of the two distributions
  thermal_overlap2[i] <- 1-(overlap(list(x_var, y_var))$OV+(1-(overlap(list(y_var2, y_var))$OV)))
  
  #Baseline 3
  # focal time frame
  x_var <- rnorm(100000, mean = ranges3[i,'mean']+15, sd = ranges3[i,'sd']+8)
  # baseline time frame
  y_var <- rnorm(100000, mean = 15, sd = 8)
  # baseline time frame2
  y_var2 <- rnorm(100000, mean = 15, sd = 8)
  
  # estimate the overlap of the two distributions
  thermal_overlap3[i] <- 1-(overlap(list(x_var, y_var))$OV+(1-(overlap(list(y_var2, y_var))$OV)))
}


ranges1$ov<-thermal_overlap1
ranges1$type<-"Baseline 1 
(mean = 15 °C, sd = 4 °C)"
ranges1$temp<-rnorm(10000, mean = 15, sd = 4)
ranges2$ov<-thermal_overlap2
ranges2$type<-"Baseline 2 
(mean = 15 °C, sd = 6 °C)"
ranges2$temp<-rnorm(10000, mean = 15, sd = 6)
ranges3$ov<-thermal_overlap3
ranges3$type<-"Baseline 3 
(mean = 15 °C, sd = 8 °C)"
ranges3$temp<-rnorm(10000, mean = 15, sd = 8)

ranges<-data.table(rbind(ranges1,ranges2,ranges3))
head(ranges)

ranges[,ov_smooth:=gam(ov~s(sd,mean),method="REML")$fitted.values,.(type)]
a<-ggplot()+
  geom_density(data=ranges,aes(x=temp,colour=(type),fill=(type)),alpha=0.2,adjust=2,size=1)+
  scale_colour_manual(values=c("#fbb4b9","#f768a1","#ae017e"))+
  scale_fill_manual(values=c("#fbb4b9","#f768a1","#ae017e"),guide=NULL)+
  # scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  # scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  labs(x="Temperature (°C)",
       y="Kernel density estimation",
       fill=NULL,
       colour=NULL)

b<-ggplot()+
  geom_raster(data=ranges,aes(x=sd, y=mean, fill=ov_smooth),interpolate=TRUE)+
  geom_contour(data=ranges,aes(x=sd, y=mean, z=ov_smooth),colour="white",breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  scale_fill_viridis_c(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  facet_wrap(~type,ncol=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_hline(yintercept = 0,colour="white",linetype="dashed")+
  geom_vline(xintercept = 0,colour="white",linetype="dashed")+
  theme_bw()+
  labs(x="Change in standard deviation (°C)",
       y="Change in mean temperature (°C)",
       fill="Thermal 
non-overlap")

plot_grid(a,b,ncol=1,rel_heights = c(1,3))
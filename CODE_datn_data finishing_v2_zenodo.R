#load relevant packages
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(zyp)
library(dismo)
library(gbm)
library(MASS)
library(data.table)
library(scales)
library(lubridate)
library(Kendall)
library(gridExtra)
library(boot)
library(spatstat)
library(rowr)
library(mvmeta)
library(zoo)
library(imputeTS)
library(mgcv)
library(overlapping)
library(reshape)
library(infer)


#set working directory
setwd("D:/DATN/data/Raw Temp Data")
setwd("D:/DATN/data")
setwd("D:/DATN")

memory.limit(size=500000000000000)

#functions
senP<-function(x,y,reps){
  slopes<-zyp.sen(y~x)$slopes
  intercepts<-zyp.sen(y~x)$intercepts
  slope_ps<-vector('numeric')
  intercept_ps<-vector('numeric')
  n<-length(x)
  
  for (i in 1:reps){
  slope_ps[i]<-wilcox.test(sample(slopes,size=n,replace=TRUE))$p.value
  intercept_ps[i]<-wilcox.test(sample(intercepts,size=n,replace=TRUE))$p.value
  }
  
  slope_p<-median(slope_ps)
  intercept_p<-median(intercept_ps)
  result <- list(slope_p=slope_p,intercept_p=intercept_p)
  return(result)
}


rm(list=setdiff(ls(),"senP"))
#lake<-"CrystalBog"


####make OV calcualtions
lakenames_lil<-fread("D:/DATN/LakeNames_lil_v1.csv")
lakenames_big<-fread("D:/DATN/LakeNames_big_v1.csv")
lakenames<-c(lakenames_lil$LakeNames_lil,lakenames_big$LakeNames_big)
lakenames<-lakenames[lakenames!="Chignik"]
#lakenames<-c("Lake114","Lake227","Lake302","RedChalkMain")
depthseq<-fread("D:/DATN/DATA_depthseq_v3.csv")


for (lake in lakenames){
  tryCatch({
    
#get the lake data
    print(lake)
    print(Sys.time())
lakedata_refit<-fread(paste0("D:/DATN/RefitData_v6/",lake,".csv",sep=""))
lakedata_refit<-unique(lakedata_refit[,c(1,2,3,4,5,6,7,8,9,10,11,12,20,21,25,26,27)])
lakedata_refit<-lakedata_refit[is.finite(fit)==TRUE]

lakedata_refit$date<-as_date(lakedata_refit$date)
lakedata_refit$doy<-yday(lakedata_refit$date)
lakedata_refit$decyear<-decimal_date(lakedata_refit$date)
lakedata_refit$year<-year(lakedata_refit$date)
lakedata_refit$lake<-lake
colnames(lakedata_refit)[colnames(lakedata_refit)=="fit"]<-"predict_refit"

if (lake=="Erie"|
    lake=="Ontario"|
    lake=="Michigan"|
    lake=="Huron"|
    lake=="Victoria"|
    lake=="Superior") {
  lakedata_refit[,depthdays:=NULL]
  lakedata_refit<-merge(lakedata_refit,depthseq)
  lakedata_refit$voldays<-lakedata_refit$depthdays*lakedata_refit$mf_sa #in units of km^3*days
}


lakedata_refit[,sampleprob:=voldays]
lakedata_refit[is.na(sampleprob)==TRUE]$sampleprob<-0

test<-lakedata_refit[,.N,.(year)]

test95<-test[N>mean(test$N)*.95]
lakedata_refit95<-lakedata_refit[year%in%test95$year]
lakedata_refit95[,count:=.N,.(doy)]
lakedata_refit95<-lakedata_refit95[count>0.95*max(lakedata_refit95$count)]

test90<-test[N>mean(test$N)*.90]
lakedata_refit90<-lakedata_refit[year%in%test90$year]
lakedata_refit90[,count:=.N,.(doy)]
lakedata_refit90<-lakedata_refit90[count>0.95*max(lakedata_refit90$count)]

test85<-test[N>mean(test$N)*.85]
lakedata_refit85<-lakedata_refit[year%in%test85$year]
lakedata_refit85[,count:=.N,.(doy)]
lakedata_refit85<-lakedata_refit85[count>0.95*max(lakedata_refit85$count)]

lakedata_refit<-lakedata_refit95
if(nrow(lakedata_refit90)>nrow(lakedata_refit)){lakedata_refit<-lakedata_refit90}
if(nrow(lakedata_refit85)>nrow(lakedata_refit)){lakedata_refit<-lakedata_refit85}

yearnumber<-length(unique(lakedata_refit$year))
years<-(unique(lakedata_refit$year))
years<-years[order(years)]

#i<-10
# 
#prep to calcualte ED
for (i in c(round(quantile(c(1:length(years)),probs=c(0.3,0.5,0.7))))){
  tryCatch({

    print(years[i])

    for(zmd in c(1,2,3,6,9,13,20)){
      for(smd in c(1,2,3,6,9,13,20)){
        #for(smd in c(1.1)){

        #smd<-1
        #zmd<-1

        for(rep.i in c(1:10)){

          if(lake=="Erie"|
             lake=="Ontario"|
             lake=="Michigan"|
             lake=="Huron"|
             lake=="Victoria"|
             lake=="Superior"){lakedata_refit<-lakedata_refit[,c(1:20)]} else
               {lakedata_refit<-lakedata_refit[,c(1:19)]}



        voldays.lake<-unique(lakedata_refit[,.(depth,voldays)])
        voldays.lake$voldays.surfsum<-cumsum(voldays.lake$voldays)
        vol.max<-sum(voldays.lake$voldays)
        voldays.lake$zmd.i<-ceiling(voldays.lake$voldays.surfsum/(vol.max/zmd))
        lakedata_refit<-merge(lakedata_refit,voldays.lake,by=c("depth","voldays"))

        doydays.lake<-unique(lakedata_refit[,.(doy)])
        doy.max<-max(doydays.lake$doy)-min(doydays.lake$doy)+1
        doydays.lake$smd.i<-round(doydays.lake$doy/(doy.max/smd))
        doydays.lake$smd.i<-doydays.lake$smd.i-min(doydays.lake$smd.i)+1
        doydays.lake[smd.i>smd]$smd.i<-smd
        lakedata_refit<-merge(lakedata_refit,doydays.lake,by="doy")

        lakedata_yearA<-lakedata_refit[decyear<years[i]][,.SD[sample(1:.N, 10000,replace=TRUE,prob=sampleprob)],.(zmd.i,smd.i)]
        lakedata_yearB<-lakedata_refit[decyear>=years[i]][,.SD[sample(1:.N, 10000,replace=TRUE,prob=sampleprob)],.(zmd.i,smd.i)]
        colnames(lakedata_yearB)<-paste0(colnames(lakedata_yearB),"X")
        lakedata_exp<-as.data.table(cbind(data.frame(lakedata_yearA),
                                          #data.frame(lakedata_yearA2),
                                          data.frame(lakedata_yearB)))
        rm(lakedata_yearA)
        rm(lakedata_yearB)
        lakedata_exp<-lakedata_exp[is.finite(predict_refit)==TRUE&is.finite(predict_refitX)==TRUE]

        lakedata.ov<-lakedata_exp[,.(ov=overlap(list(predict_refit,predict_refitX),nbins=100)$OV),.(zmd.i,smd.i)]
        
        rm(lakedata_exp)

        lakedata_yearA<-lakedata_refit[decyear<years[i]][,.SD[sample(1:.N, 10000,replace=TRUE,prob=sampleprob)],.(zmd.i,smd.i)]
        lakedata_yearA2<-lakedata_refit[decyear<years[i]][,.SD[sample(1:.N, 10000,replace=TRUE,prob=sampleprob)],.(zmd.i,smd.i)]
        colnames(lakedata_yearA2)<-paste0(colnames(lakedata_yearA2),"2")
        lakedata_exp.null<-as.data.table(cbind(data.frame(lakedata_yearA),
                                          #data.frame(lakedata_yearA2),
                                          data.frame(lakedata_yearA2)))
        rm(lakedata_yearA)
        rm(lakedata_yearA2)
        lakedata_exp.null<-lakedata_exp.null[is.finite(predict_refit)==TRUE&is.finite(predict_refit2)==TRUE]

        lakedata.ov.null<-lakedata_exp.null[,.(ov.null=overlap(list(predict_refit,predict_refit2),nbins=100)$OV),.(zmd.i,smd.i)]
        
        rm(lakedata_exp.null)

        lakedata.ov.m<-lakedata.ov
        lakedata.ov.m$zmd<-zmd
        lakedata.ov.m$smd<-smd
        lakedata.ov.m$lake<-lake
        lakedata.ov.m$split<-i
        lakedata.ov.m$split.year<-years[i]
        lakedata.ov.m$rep<-rep.i

        # if the merged dataset does exist, append to it
        if (exists("lakedata_ovreps")){lakedata_ovreps<-rbindlist(list(lakedata_ovreps,lakedata.ov.m),fill=TRUE)}
        # if the merged dataset doesn't exist, create it
        if (!exists("lakedata_ovreps")){lakedata_ovreps<-lakedata.ov.m}

        }

        
        }
      }

    fwrite(lakedata_ovreps,paste0("D:/DATN/OV.smd.zmd/DATA_",lake,"_split",i,"_v1b.csv",sep=""))

    rm(lakedata_ovreps)
    gc()

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
       

  lakedata_warm<-lakedata_refit[,
                                .(meantemp=weighted.mean(predict_refit,voldays,na.rm=TRUE)),by=year]
  lakedata_warming<-lakedata_warm[,
                                  .(tempmb_slope=zyp.sen(meantemp~year)$coefficients[2],
                                    tempmb_p=senP(y=meantemp,x=year,reps=25)$slope_p)]
  lakedata_warming$lake<-lake  
  
  # Merge the data from years into a single file
  # if the merged dataset doesn't exist, create it
  if (!exists("lakedata_temptrends")){lakedata_temptrends<-lakedata_warming}
  
  # if the merged dataset does exist, append to it
  if (exists("lakedata_temptrends")){lakedata_temptrends<-rbindlist(list(lakedata_temptrends,lakedata_warming),fill=TRUE)}
  
  
  lakedata_timeseries<-lakedata_refit[is.finite(predict_refit)==TRUE,.(temp_mean=mean(predict_refit,na.rm=TRUE),#add a check for is.finite
                                            temp_min=min(predict_refit,na.rm=TRUE),
                                            temp_max=max(predict_refit,na.rm=TRUE),
                                            decyear_mean=mean(decyear,na.rm=TRUE),
                                            decyear_min=min(decyear,na.rm=TRUE),
                                            decyear_max=max(decyear,na.rm=TRUE),
                                            doy_mean=mean(doy,na.rm=TRUE),
                                            doy_min=min(doy,na.rm=TRUE),
                                            doy_max=max(doy,na.rm=TRUE),
                                            interp.dist_mean=mean(interp.dist,na.rm=TRUE))]
  lakedata_timeseries$lake<-lake  
  
  # Merge the data from years into a single file
  # if the merged dataset doesn't exist, create it
  if (!exists("lakedata_time")){lakedata_time<-lakedata_timeseries}
  
  # if the merged dataset does exist, append to it
  if (exists("lakedata_time")){lakedata_time<-rbindlist(list(lakedata_time,lakedata_timeseries),fill=TRUE)}
  

fwrite(lakedata_temptrends,"D:/DATN/DATN_temptrends_v42b.csv")
#fwrite(lakedata_exp_data,"D:/DATN/DATN_EDdata_v40b_crystalbog_delete.test.csv")
fwrite(lakedata_time,"D:/DATN/DATN_EDtime_v42b.csv")

gc()

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


rm(list=setdiff(ls(),"senP"))


files<-list.files(path = "D:/DATN/OV.smd.zmd", pattern=".csv",full.names = TRUE)
files<-files[grep("DATA_",files)]
files<-files[grep("_randnull_",files,invert=TRUE)]

OVdata<- do.call("rbindfill", lapply(files, FUN = function(file) {
  read.table(file, header=TRUE, sep="\t")
}))

OVdata<-do.call("rbind.fill",lapply(files,FUN=function(files){read.csv(files)}))
OVdata<-unique(OVdata)
summary(OVdata)

files<-list.files(path = "D:/DATN/OV.smd.zmd", pattern=".csv",full.names = TRUE)
files<-files[grep("DATA_",files)]

OVdata.randnull<-do.call("rbind",lapply(files,FUN=function(files){read.csv(files)}))
OVdata.randnull<-unique(OVdata.randnull)
summary(OVdata.randnull)
length(unique(OVdata.randnull$lake))
lakenames<-setdiff(unique(OVdata.randnull$lake), unique(OVdata$lake))

fwrite(OVdata,"D:/DATN/DATA_OVdata_bylake.zmd.smd.zmdi.smdi.rep.split_v1.csv")
OVdata<-fread("D:/DATN/DATA_OVdata_bylake.zmd.smd.zmdi.smdi.rep.split_v1.csv")
length(unique(OVdata$lake))

fwrite(OVdata.randnull,"D:/DATN/DATA_OVdata.randnull_bylakesplit_v1.csv")
OVdata.randnull<-fread("D:/DATN/DATA_OVdata.randnull_bylakesplit_v1.csv")

names(OVdata.randnull)[3]<-"ov.randnull"
OVdata.m<-merge(OVdata,OVdata.randnull,all.x=TRUE,by=c("zmd.i","smd.i","smd","zmd","lake","rep"))

fwrite(OVdata.m,"D:/DATN/DATA_OVdata_bylake.zmd.smd.zmdi.smdi.rep.sig.split_v1.csv")
OVdata<-fread("D:/DATN/DATA_OVdata_bylake.zmd.smd.zmdi.smdi.rep.sig.split_v1.csv")


summary(OVdata)
OVdata.means<-OVdata[,.(ov=mean(ov),
                  #ov.null=mean(ov.null,na.rm=T),
                  ov.randnull=mean(ov.randnull),
                  ov.std=median((1-ov.randnull)+ov),
                  ov.stdq=median(ov/ov.randnull)),
               .(zmd,smd,lake,split,split.year)]
summary(OVdata.means)
hist(OVdata.means$ov.std)

OVdata.medians<-OVdata[,.(ov=median(ov),
                        #ov.null=median(ov.null,na.rm=T),
                        ov.randnull=median(ov.randnull),
                        ov.std=median((1-ov.randnull)+ov),
                        ov.stdq=median(ov/ov.randnull)),
                     .(zmd,smd,lake,split,split.year)]
summary(OVdata.medians)
hist(OVdata.medians$ov.std)

OVdata.sig.test<-OVdata[smd==1&zmd==1,.(ov_p=wilcox.test(ov,ov.randnull,paired=TRUE)$p.value,
                           ov_95.l=wilcox.test(ov,ov.randnull,paired=TRUE,conf.int=T)$conf.int[1],
                           ov_95.u=wilcox.test(ov,ov.randnull,paired=TRUE,conf.int=T)$conf.int[2]),.(lake)]


fwrite(OVdata.sig.test,"D:/DATN/DATA_OVdata_sig.test.mean_v1.csv")
OVdata.sig.test<-fread("D:/DATN/DATA_OVdata_sig.test.mean_v1.csv",stringsAsFactors = T)

#Merge the overlap values with the signficiance values.
OVdata.m<-merge(OVdata.means,OVdata.sig.test,all.x=TRUE,by=c("smd","zmd","lake"))
summary(OVdata.m)
#ovdata[,count:=.N,.(lake)]
fwrite(OVdata.m,"D:/DATN/DATA_OVdata_bylake.zmd.smd.sig.split_v2.csv")
ovdata<-fread("D:/DATN/DATA_OVdata_bylake.zmd.smd.sig.split_v2.csv")

OVtime<-fread("D:/DATN/DATN_EDtime_v42.csv")
fwrite(OVtime,"D:/DATN/DATA_OVtime_v1.csv")

OVtemp<-fread("D:/DATN/DATN_temptrends_v42.csv")
fwrite(OVtemp,"D:/DATN/DATA_OVtemp_v1.csv")



####Look for important predicting lake caracteristics from hydrolakes (temperature, warming, and morphometry)####
#DATN_models<-fread("DATN_models_v12.csv",header=TRUE)
DATN_OV<-fread("D:/DATN/DATA_OVdata_bylake.zmd.smd.sig.split_v2.csv",header=TRUE,stringsAsFactors = TRUE)
#DATN_OV[lake=="Zurich"&doymaxdist==1.1&depthmaxdist==1.1,.(OVmean=mean(OVmean)-mean(OVmean.null))]

DATN_lakeinfo<-fread("D:/DATN/DATA_datn_HydroLakesInfo_v2.csv",header=TRUE,stringsAsFactors = TRUE) #neOVs OViting
DATN_time<-fread("D:/DATN/DATA_OVtime_v1.csv",header=TRUE,stringsAsFactors = TRUE)
DATN_warm<-fread("D:/DATN/DATA_OVtemp_v1.csv",header=TRUE,stringsAsFactors = TRUE)

DATN_predict<-merge(DATN_OV,unique(DATN_time),all.x=TRUE)
DATN_predict<-merge(DATN_predict,DATN_lakeinfo,all.x=TRUE)
DATN_predict<-unique(DATN_predict)
DATN_warm<-unique(DATN_warm)
DATN_warm<-DATN_warm[,lapply(.SD, mean),.(lake)]
DATN_predict<-merge(DATN_predict,DATN_warm,all.x=TRUE,allow.cartesian=TRUE)
DATN_predict<-unique(DATN_predict)

DATN_predict$abslat<-abs(DATN_predict$lat)
DATN_predict$loghl.Lake_area<-log(DATN_predict$hl.Lake_area)
DATN_predict$loghl.Elevation<-log(DATN_predict$hl.Elevation-min(DATN_predict$hl.Elevation,na.rm=TRUE)+100)
DATN_predict$loghl.Depth_avg<-log(DATN_predict$hl.Depth_avg)
DATN_predict$logMeanDepth<-log(DATN_predict$MeanDepth*1000)
DATN_predict$abstempmb_slope<-abs(DATN_predict$tempmb_slope)
DATN_predict$loghl.Res_time<-log(DATN_predict$hl.Res_time)
DATN_predict$loghl.Shore_dev<-log(DATN_predict$hl.Shore_dev)


fwrite(DATN_predict,"D:/DATN/DATA_predict_v50.csv")

DATN_predict<-fread("D:/DATN/DATA_predict_v50.csv",header=TRUE,stringsAsFactors = TRUE)
DATN_predict<-DATN_predict[order(lake,split.year)]
DATN_predict<-DATN_predict[,split.year.rank:=rleid(split.year),.(lake)]

DATN_ts.trends<-fread("D:/DATN/DATA_ansum.trends_v1.csv",stringsAsFactors = T)
DATN_predict<-merge(DATN_predict,DATN_ts.trends)

DATN_predict$logmeandoys<-log(DATN_predict$meandoys)
DATN_predict$nov<-100*(1-DATN_predict$ov)
DATN_predict$nov.randnull<-100*(1-DATN_predict$ov.randnull)


LearningRate<-c(0.8,
  0.4096,
  0.2048,
  0.1024,0.0512,
  0.0256,
  0.0128,0.0064,0.0032,0.0016,0.0008,0.0004,0.0002,0.0001,0.00005,0.00001,0.000005,0.000001,0.0000005,0.0000001,0.00000005,0.00000001)
bf<-0.1
lakes<-data.table(unique(DATN_predict$lake))

for(i in 1:10){
  
  lakes.i.samp<-data.table(sample(lakes$V1,97,replace=F))
  
  # lakes.i.samp<-lakes.i[,..i]
  names(lakes.i.samp)<-"lakes.samp"
  DATN_predict.i<-DATN_predict[lake%in%lakes.i.samp$lakes.samp]

#brt for time series characteristics, shift limits
    for(lr in LearningRate){
  tryCatch({
    
    BRT<-gbm.step(
      data=DATN_predict.i[is.na(ov)==FALSE], #maybe include cutoff for seasonal coverage
      gbm.x = c(
        #1,#lake
        #34,#continent
        #70,#spaceclustr
        3,#depth limit
        2,#doy limit
        ####5,#yearsplit
        ####10,#count
        #11,#temp mean
        15,#decyear_min
        16,#decyear_max
        18,#doy_min
        19,#doy_max
        ####24,#interpolation distance mean
        #58,#log shorline development index
        #57,#log residence time
        51,#abslat
        #52,#log surface area
        #53,#log elevation
        55,#,#log mean depth
        59,#,#yearsplit.rank
        60,#doys.trend
        61,#doys.coverage.trend
        63#,#log doy means
        #65#nov.randnull
        ),
      gbm.y = 64,#nov
      family = "gaussian", 
      tree.complexity = 18,
      learning.rate = lr,
      max.trees=10000,
      bag.fraction = bf,
      site.weights = (1-DATN_predict.i$ov_p))
    
    if(is.null(BRT)==FALSE){
      if(BRT$n.trees>1000) break}
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveRDS(BRT,paste0("D:/DATN/MODEL_OVbrt_bf.",bf,"_rep.",i,"_allPreds.mt_v3.rds"))


}

####Standardizing the nov values across time series characteristics
#Prep the data table
DATN_predict_100<-DATN_predict[rep(DATN_predict[,.I],100)]
DATN_predict_100$j<-rep(c(1:100),each=nrow(DATN_predict))
DATN_predict_100$nov.std.predict<-numeric()

DATN_predict_std<-DATN_predict[rep(DATN_predict[,.I],100)]
DATN_predict_std$j<-rep(c(1:100),each=nrow(DATN_predict))
DATN_predict_std$nov.std.predict<-numeric()

#Standardize the predictors
DATN_predict_std$doy_min<-mean(DATN_predict$doy_min)
DATN_predict_std$doy_max<-mean(DATN_predict$doy_max)
DATN_predict_std$decyear_min<-mean(DATN_predict$decyear_min)
DATN_predict_std$decyear_max<-mean(DATN_predict$decyear_max)
DATN_predict_std$yearsplit.rank<-2
DATN_predict_std$doys.trend<-0
DATN_predict_std$doys.coverage.trend<-0
DATN_predict_std$logmeandoys<-median(DATN_predict$logmeandoys,na.rm=TRUE)

#run the standardization for each BRT
for(i in c(1:100)){
  print(i)
  BRT.i<-readRDS(paste0("D:/DATN/MODEL_OVbrt_bf.0.1_rep.",i,"_allPreds.mt_v3.rds"))
  
  DATN_predict_std[j==i]$nov.std.predict<-(predict.gbm(BRT.i,
                                                       newdata=DATN_predict_std[j==i],
                                                       n.trees=BRT.i$gbm.call$best.trees, 
                                                       type="response"))
  DATN_predict_100[j==i]$nov.std.predict<-(predict.gbm(BRT.i,
                                                       newdata=DATN_predict_100[j==i],
                                                       n.trees=BRT.i$gbm.call$best.trees, 
                                                       type="response"))
  
}

fwrite(DATN_predict_std,"D:/DATN/DATA_DATN_predict_std_allreps_v51.csv")
fwrite(DATN_predict_100,"D:/DATN/DATA_DATN_predict_100_allreps_v51.csv")
DATN_predict_std<-fread("D:/DATN/DATA_DATN_predict_std_allreps_v51.csv",stringsAsFactors = T)
DATN_predict_100<-fread("D:/DATN/DATA_DATN_predict_100_allreps_v51.csv",stringsAsFactors = T)

#Merge the two datasets
hist(DATN_predict_std$nov.std.predict)
summary(DATN_predict_std)
summary(DATN_predict_100)
head(DATN_predict_100)
names(DATN_predict_std)
names(DATN_predict_100)
length(unique(DATN_predict_100$lake))
length(unique(DATN_predict_std$lake))

####Calculate values for the supplementary table

DATN_predict_100[,':='(meanresid.rank=rank(meanresid)),.(j)]

DATN_predict_100.test<-DATN_predict_100[,.(meanresid=mean(abs(nov.std.predict-nov))),.(j,lake)]
DATN_predict_100.test[,':='(meanresid.rank=rank(meanresid)),.(j)]
DATN_predict_100.test.good<-DATN_predict_100.test[meanresid.rank<98]
DATN_predict_100$TestTrain<-factor("Test")

for(i in c(1:100)){
  print(i)
  DATN_predict_100.test.good.i<-DATN_predict_100.test.good[j==i]
  DATN_predict_100[j==i&lake%in%DATN_predict_100.test.good.i$lake]$TestTrain<-"Train"
}

fwrite(DATN_predict_100,"D:/DATN/DATA_DATN_predict_100_allreps.forSItable_v51.csv")

DATN_predict_100.mean<-DATN_predict_100[,lapply(.SD, mean),
                                      by=c(colnames(DATN_predict_100)[c(1:5,11:63,68)])]

#Training data correlation between predicted and observed values
DATN_predict_100.samp<-DATN_predict_100[sample(.N,10000)]
cor(DATN_predict_100.samp[TestTrain=="Train"]$nov,
    DATN_predict_100.samp[TestTrain=="Train"]$nov.std.predict,
    method="kendall")

#Training data correlation between predicted and observed values
cor(DATN_predict_100.samp[TestTrain=="Test"]$nov,
    DATN_predict_100.samp[TestTrain=="Test"]$nov.std.predict,
    method="kendall")

#Correlation between predicted and observed values where predictions were averaged across all BRTs
cor(DATN_predict_100.mean$nov,
    DATN_predict_100.mean$nov.std.predict,
    method="kendall")

#Standard deviation (SD) of thermal non-overlap values
sd(DATN_predict_100.mean$nov)

#SD of model residuals (RMSE) where predictions were averaged across all BRTs
sd(DATN_predict_100.mean$nov-DATN_predict_100.mean$nov.std.predict)

#SD of cross validation model residuals (PRESS) where SD was averaged across all BRTs
sd(DATN_predict_100.samp[TestTrain=="Test"]$nov-DATN_predict_100.samp[TestTrain=="Test"]$nov.std.predict)

#Median absolute deviation (MAD) of thermal non-overlap values
mad(DATN_predict_100.mean$nov,constant=1)

#MAD of model residuals (RMSE) where predictions were averaged across all BRTs
mad(DATN_predict_100.mean$nov-DATN_predict_100.mean$nov.std.predict,constant=1)

#MAD of cross validation model residuals (PRESS) where SD was averaged across all BRTs
mad(DATN_predict_100.samp[TestTrain=="Test"]$nov-DATN_predict_100.samp[TestTrain=="Test"]$nov.std.predict,constant=1)

##average predictions across models and test for performance with training data
DATN_predict_100.cvTrain<-DATN_predict_100[TestTrain==Train,.(cor)]

#Average predictions across models for Test dataset and test for performance with test data
DATN_predict_100.cvTest<-DATN_predict_100


DATN_predict_100.cv<-DATN_predict_100[,,.(TestTrain)]
unique(DATN_predict_100[meanresid.rank<98000&j==1]$lake)

summary(DATN_predict_100)
head(DATN_predict_100)

goodlakes.sub<-goodlakes

#plot(BRT.i,"logMeanDepth")
summary(BRT.i$var.levels)[2]
hist(DATN_predict_100$nov.std.predict-DATN_predict_100$nov)

names(DATN_predict_100)[67]<-"nov.100.predict"
DATN_predict.m<-merge(DATN_predict_100,DATN_predict_std[,c(1,2,3,4,66,67)],ALL=T)
DATN_predict.m[,resid:=(nov-nov.100.predict)]
DATN_predict.m[,nov.std.predict.resid:=nov.std.predict+resid-median(DATN_predict$nov.randnull)]
DATN_predict.m[nov.std.predict.resid<0]$nov.std.predict.resid<-0
DATN_predict.m[nov.std.predict<0]$nov.std.predict<-0
DATN_predict.m[nov.100.predict<0]$nov.100.predict<-0

summary(DATN_predict.m)
names(DATN_predict.m)

fwrite(DATN_predict.m,"D:/DATN/DATA_predict.m_v51.csv")

DATN_predict.m.mean<-(DATN_predict.m)[,lapply(.SD, mean),
                                      by=c(colnames(DATN_predict.m)[c(1:4,6:66)])]

nrow(DATN_predict.m.mean[nov.100.predict==0])
nrow(DATN_predict.m.mean[nov.std.predict==0])
nrow(DATN_predict.m.mean[nov.std.predict.resid==0])

hist(DATN_predict.m.mean$nov.100.predict,100)
hist(DATN_predict.m.mean$nov.std.predict,100)
hist(DATN_predict.m.mean$nov.std.predict.resid,100)
hist(DATN_predict.m.mean$nov,100)

summary(DATN_predict.m.mean)
names(DATN_predict.m.mean)

fwrite(DATN_predict.m.mean,"D:/DATN/DATA_predict.m.mean_v51.csv")

####Compile the model contribution and part dep data from all of the BRTs
#Prep the data for contribs
BRT.contrib<-data.table(j=rep(c(1:100),each=13),
                        var=rep(vars,times=100),
                        rel.inf=numeric())
head(BRT.contrib)

#prep the data for the part deps
vars<-c(colnames(DATN_predict)[c(3,2,15,16,18,19,51,55,59,60,61,63,65)])
BRT.partdep<-data.table(j=rep(c(1:100),each=130),
                        var=factor(rep(vars,times=10000)),
                        var.value=numeric(),
                        nov=numeric())
summary(BRT.partdep)
head(BRT.partdep,100)

#i<-46
#var.i<-vars[11]

for(i in c(1:100)){
  print(i)
  BRT.i<-readRDS(paste0("D:/DATN/MODEL_OVbrt_bf.0.1_rep.",i,"_allPreds.mt_v3.rds"))
  
  BRT.contrib[j==i]$var<-BRT.i$contributions$var
  BRT.contrib[j==i]$rel.inf<-BRT.i$contributions$rel.inf
    
  for(var.i in vars){
    BRT.partdep[j==i&var==var.i]$var<-var.i
    BRT.partdep[j==i&var==var.i]$var.value<-plot(BRT.i,var.i,return.grid=T)[,1]
    BRT.partdep[j==i&var==var.i]$nov<-plot(BRT.i,var.i,return.grid=T)$y
  }
}

BRT.partdep<-na.omit(BRT.partdep)

fwrite(BRT.contrib,"D:/DATN/DATA_DATN_BRT.contrib_v51.csv")
fwrite(BRT.partdep,"D:/DATN/DATA_DATN_BRT.partdep_v51.csv")

####Find the interactions between shift limits and lake characteristics
#prep the data
DATN_predict_std.time<-fread("D:/DATN/DATA_predict.m.mean_v51.csv",stringsAsFactors = T)

meandepths<-seq(from = min(DATN_predict_std.time[is.na(logMeanDepth)==FALSE]$logMeanDepth), 
                to = max(DATN_predict_std.time[is.na(logMeanDepth)==FALSE]$logMeanDepth),
                by = (max(DATN_predict_std.time[is.na(logMeanDepth)==FALSE]$logMeanDepth)-
                        min(DATN_predict_std.time[is.na(logMeanDepth)==FALSE]$logMeanDepth))/1000)
lats<-seq(from = min(DATN_predict_std.time$abslat), 
          to = max(DATN_predict_std.time$abslat),
          by = (max(DATN_predict_std.time$abslat)-
                  min(DATN_predict_std.time$abslat))/1000)

smds<-unique(DATN_predict_std.time$smd)
zmds<-unique(DATN_predict_std.time$zmd)

DATN_predict_new1<-expand.grid(logMeanDepth=meandepths,smd=smds)
DATN_predict_new1$abslat<-mean(DATN_predict$abslat)
DATN_predict_new1$zmd<-1
DATN_predict_new1$type<-1

DATN_predict_new2<-expand.grid(abslat=lats,smd=smds)
DATN_predict_new2$logMeanDepth<-mean(DATN_predict$logMeanDepth)
DATN_predict_new2$zmd<-1
DATN_predict_new2$type<-2

DATN_predict_new3<-expand.grid(logMeanDepth=meandepths,zmd=zmds)
DATN_predict_new3$abslat<-mean(DATN_predict$abslat)
DATN_predict_new3$smd<-1
DATN_predict_new3$type<-3

DATN_predict_new4<-expand.grid(abslat=lats,zmd=zmds)
DATN_predict_new4$logMeanDepth<-mean(DATN_predict$logMeanDepth)
DATN_predict_new4$smd<-1
DATN_predict_new4$type<-4


DATN_predict_new<-as.data.table(rbind.fill(DATN_predict_new1,DATN_predict_new2,DATN_predict_new3,DATN_predict_new4))
DATN_predict_new<-DATN_predict_new[rep(DATN_predict_new[,.I],100)]

DATN_predict_new$doy_min<-mean(DATN_predict$doy_min)
DATN_predict_new$doy_max<-mean(DATN_predict$doy_max)
DATN_predict_new$decyear_min<-mean(DATN_predict$decyear_min)
DATN_predict_new$decyear_max<-mean(DATN_predict$decyear_max)
DATN_predict_new$split.year.rank<-2
DATN_predict_new$doys.trend<-0
DATN_predict_new$doys.coverage.trend<-0
DATN_predict_new$logmeandoys<-median(DATN_predict$logmeandoys,na.rm=TRUE)
#DATN_predict_new$nov.randnull<-0
DATN_predict_new$nov.randnull<-sample(DATN_predict$nov.randnull,nrow(DATN_predict_new),replace=TRUE)

DATN_predict_new[type==1]$abslat<-mean(DATN_predict$abslat)
DATN_predict_new[type==2]$logMeanDepth<-mean(DATN_predict$logMeanDepth)
DATN_predict_new[type==3]$abslat<-mean(DATN_predict$abslat)
DATN_predict_new[type==4]$logMeanDepth<-mean(DATN_predict$logMeanDepth)
DATN_predict_new$rep<-rep(c(1:100),each=nrow(DATN_predict_new)/100)
DATN_predict_new$nov.predict.int<-numeric()
head(DATN_predict_new)

i<-1

for(i in c(1:100)){
  print(i)
  BRT.i<-readRDS(paste0("D:/DATN/MODEL_OVbrt_bf.0.1_rep.",i,"_allPreds.mt_v3.rds"))
  DATN_predict_new[rep==i]$nov.predict.int<-predict.gbm(BRT.i,
                                                        newdata=DATN_predict_new[rep==i],
                                                        n.trees=BRT.i$gbm.call$best.trees, 
                                                        type="response")
}

DATN_predict_new$nov.predict.int<-DATN_predict_new$nov.predict.int-
  median(DATN_predict$nov.randnull)

fwrite(DATN_predict_new,"D:/DATN/DATA_DATN_predict_new_allreps_v51.csv")
DATN_predict_new<-fread("D:/DATN/DATA_DATN_predict_new_allreps_v51.csv",stringsAsFactors = T)

DATN_predict_new.means<-DATN_predict_new[,lapply(.SD, mean),
                                           by=c(colnames(DATN_predict_new)[c(1:13)])]
DATN_predict_new.means$nov.predict.int<-(DATN_predict_new.means$nov.predict.int-min(DATN_predict_new.means$nov.predict.int))*
  (1+(min(DATN_predict_new.means$nov.predict.int)/max(DATN_predict_new.means$nov.predict.int)))

fwrite(DATN_predict_new.means,"D:/DATN/DATA_DATN_predict_new.means_allreps_v51.csv")
DATN_predict_new.means<-fread("D:/DATN/DATA_DATN_predict_new.means_allreps_v51.csv",stringsAsFactors = T)

hist(DATN_predict_new.means$nov.predict.int)

####Interpolate the shift limits
DATN_predict_std.time<-fread("D:/DATN/DATA_predict.m.mean_v51.csv",stringsAsFactors = T)  
 
summary(DATN_predict_std.time) 
unique(DATN_predict_std.time[nov.std.predict.resid==0]$lake)

####Calculate the OV summary by splits separately for each lake, repeatOV with standardizOV data####
OVdata.lakesplits<-DATN_predict_std.time[,.(nov=mean(nov.std.predict.resid)),.(smd,zmd,lake)]
summary(OVdata.lakesplits)
fwrite(OVdata.lakesplits,"D:/DATN/DATN_OVdata.lakesplits.std.time_v51.csv")

#smooth the data separately for each lake, repeatOV with standardizOV data
OVdata.lakesplits<-fread("D:/DATN/DATN_OVdata.lakesplits.std.time_v51.csv",stringsAsFactors = T)
#OVdata.lakesplits<-fread("D:/DATN/DATN_OVdata.lakesplits_v51.csv",stringsAsFactors = T)

#OVdata.lakesplits[NOVmean<0]$NOVmean<-0
OVdata.lakesplits$smd<-1-(1/OVdata.lakesplits$smd)
OVdata.lakesplits$zmd<-1-(1/OVdata.lakesplits$zmd)

for(lakename in droplevels(unique(OVdata.lakesplits$lake))){
  smd<-seq(from = 0, to = 0.95, by = 0.01)
  zmd<-seq(from = 0, to = 0.95, by = 0.01)
  OVdata.splits.interp<-expand.grid(smd,zmd)
  names(OVdata.splits.interp)<-c("smd","zmd")
  fit<-gam(data=OVdata.lakesplits[lake==lakename],(nov)~s((smd),(zmd)),method="REML")
  OVdata.splits.interp$nov<-(predict(fit,OVdata.splits.interp))
  OVdata.splits.interp$lake<-lakename
  # if the mergOV dataset does exist, append to it
  if (exists("OVdata.lakesplits.interp")){OVdata.lakesplits.interp<-rbindlist(list(OVdata.lakesplits.interp,OVdata.splits.interp),fill=TRUE)}
  # if the mergOV dataset doesn't exist, create it
  if (!exists("OVdata.lakesplits.interp")){OVdata.lakesplits.interp<-OVdata.splits.interp}
}

ggplot()+
  geom_raster(data=OVdata.splits.interp,aes(x=smd,y=zmd,fill=nov))+
  scale_fill_viridis_c()+
  geom_contour(data=OVdata.splits.interp,aes(x=smd,y=zmd,z=nov))

summary(OVdata.lakesplits.interp)

fwrite(OVdata.lakesplits.interp,"D:/DATN/DATN_OVdata.lakesplits.interp.std.time_v51.csv")

####Lakewide average values####
OVdata.lakesplits.interp<-fread("D:/DATN/DATN_OVdata.lakesplits.interp.stdtime_v51.csv")
summary(OVdata.lakesplits.interp)

OVdata.lakemedians<-OVdata.lakesplits.interp[,.(nov=median(nov)),.(lake)]
summary(OVdata.lakemedians)
fwrite(OVdata.lakemedians,"D:/DATN/DATN_OVdata.lakemedians_v51.csv")

OVdata.lakemeans<-OVdata.lakesplits.interp[,.(nov=mean(nov)),.(lake)]
summary(OVdata.lakemeans)
fwrite(OVdata.lakemeans,"D:/DATN/DATN_OVdata.lakemeans_v51.csv")


####Calculate the OV summary by splits averaging over lakes, repeatOV with standardizOV data####
OVdata.lakesplits.interp<-fread("D:/DATN/DATN_OVdata.lakesplits.interp.std.time_v51.csv")
OVdata.lakesplits.interp<-OVdata.lakesplits.interp[,.(nov=median(nov)),.(smd,zmd,lake)]
names(OVdata.lakesplits.interp)
summary(OVdata.lakesplits.interp)
fit<-gam(data=OVdata.lakesplits.interp,nov~s(smd,zmd),method="REML")

smd<-seq(from = 0, to = .95, by = 0.01)
zmd<-seq(from = 0, to = .95, by = 0.01)
OVdata.splits.interp<-expand.grid(smd,zmd)
names(OVdata.splits.interp)<-c("smd","zmd")
summary(OVdata.splits.interp)
OVdata.splits.interp$nov<-predict(fit,newdata=OVdata.splits.interp)


DATN_predict_std.time<-fread("D:/DATN/DATA_predict.m.mean_v51.csv",stringsAsFactors = T)  
DATN_predict_std.time$smd<-1-(1/DATN_predict_std.time$smd)
DATN_predict_std.time$zmd<-1-(1/DATN_predict_std.time$zmd)
DATN_predict_std.time<-DATN_predict_std.time[,.(nov.std.predict.resid=median(nov.std.predict.resid)),
                                             .(smd,zmd,lake)]
names(DATN_predict_std.time)
summary(DATN_predict_std.time)

#smooth the data, repeatOV with standardizOV data
smd<-seq(from = 0, to = .95, by = 0.01)
zmd<-seq(from = 0, to = .95, by = 0.01)
OVdata.splits.interp<-expand.grid(smd,zmd)
names(OVdata.splits.interp)<-c("smd","zmd")


####Comparing shift limits and no shift limits to the overall means
#the percent change when not being able to shift across depth compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
  mean(DATN_predict.m.mean[zmd==20,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
  mean(DATN_predict.m.mean$nov.std.predict.resid))

#the percent change when not being able to shift across season compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==20,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean$nov.std.predict.resid))

#the percent change when not being able to shift across season and depth compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==20&zmd==20,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean$nov.std.predict.resid))

#the percent change when being able to shift across depth compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[zmd==1,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean$nov.std.predict.resid))

#the percent change when being able to shift across season compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==1,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean$nov.std.predict.resid))

#the percent change when being able to shift across season and depth compared to overall mean
-100*((mean(DATN_predict.m.mean$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==1&zmd==1,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean$nov.std.predict.resid))


####comparing shift limits to no shift limits
#the percent change when not being able to shift across season and depth compared to no shift limits
-100*((mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==20&zmd==1,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid))

-100*((mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==1&zmd==20,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid))

-100*((mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid)-
         mean(DATN_predict.m.mean[smd==20&zmd==20,.(nov.std.predict.resid=mean(nov.std.predict.resid)),.(lake)]$nov.std.predict.resid))/
        mean(DATN_predict.m.mean[smd==1&zmd==1]$nov.std.predict.resid))


#the volume comprised by the lakes in this study
lakeinfo<-fread("DATA_datn_lakeinfo_v3.csv")
summary(lakeinfo)
sum(lakeinfo$Volume)

#total lake surface freshwater volume
181.9-(181.9*.44)

#proportion of freshwater lake volume in this study
(sum(lakeinfo$Volume)/1000)/(181.9-(181.9*.44))

#supplementary information
data<-fread("D:/DATN/DATA_DATN_predict_100_allreps.forSItable_v51.csv")
OVdata.lakemeans<-fread("D:/DATN/DATN_OVdata.lakemeans_v51.csv",header=TRUE,stringsAsFactors = TRUE)
DATN_lakeinfo<-fread("D:/DATN/DATA_datn_HydroLakesInfo_v2.csv",header=TRUE,stringsAsFactors = TRUE) #neOVs OViting
DATN_time<-fread("D:/DATN/DATA_OVtime_v1.csv",header=TRUE,stringsAsFactors = TRUE)
DATN_warm<-fread("D:/DATN/DATA_OVtemp_v1.csv",header=TRUE,stringsAsFactors = TRUE)
OVdata.lakesig<-fread("D:/DATN/DATA_OVdata_sig.test.mean_v1.csv",stringsAsFactors = T)

DATN_SI<-merge(OVdata.lakemeans,DATN_time,all.x=TRUE)
DATN_SI<-merge(DATN_SI,DATN_lakeinfo,all.x=TRUE)
DATN_SI<-unique(DATN_SI)
DATN_warm<-unique(DATN_warm)
DATN_warm<-DATN_warm[,lapply(.SD, mean),.(lake)]
DATN_SI<-merge(DATN_SI,DATN_warm,all.x=TRUE,allow.cartesian=TRUE)
DATN_SI<-merge(DATN_SI,OVdata.lakesig,all.x=TRUE)

fwrite(DATN_SI,"D:/DATN/DATN_OVdata.lakes_SI_v6.csv")


####Find the beginning day of the year and ending day of the year for each year with data for each lake####
DATN_predict<-fread("D:/DATN/DATA_predict_v46.csv",header=TRUE,stringsAsFactors = TRUE)
names(DATN_predict)
lakenames<-as.character(unique(DATN_predict$lake))
lakename<-"Bubble"
rm(lakOVata.ansum)
rm(lakOVata.ansumALL)

for (lakename in lakenames){
  tryCatch({
  print(lakename)
  lakOVata<-fread(paste0("D:/DATN/DATN_RawBylake/DATN_",lakename,".csv"),stringsAsFactors = T)
  lakOVata_refit<-fread(paste0("D:/DATN/RefitData_v6/",lakename,".csv"),stringsAsFactors = T)
  years<-unique(lakOVata_refit$year)
  lakOVata$doy<-yday(lakOVata$date)

  lakOVata.ansum<-lakOVata[,.(doy_min.an=min(doy),
                              doy_max.an=max(doy),
                              doys=length(unique(doy))
                              ),.(year,lake)]
  lakOVata.ansum<-lakOVata.ansum[year%in%years]
  
  # if the interpolatOV dataset exhist then add to it
  if (exists("lakOVata.ansumALL")){lakOVata.ansumALL<-rbindlist(list(lakOVata.ansumALL,lakOVata.ansum),fill=TRUE)}
  # if the interpolatOV dataset doesnt exhist then create it
  if (!exists("lakOVata.ansumALL")){lakOVata.ansumALL<-lakOVata.ansum}
  
  fwrite(lakOVata.ansumALL,"D:/DATN/DATA_annualdoysummaries_v5.csv")
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

lakOVata.ansumALL<-fread("D:/DATN/DATA_annualdoysummaries_v5.csv")

DATN_doys<-unique(DATN_predict[,c(1,22,23)])
lakOVata.ansumALL<-merge(lakOVata.ansumALL,DATN_doys,by="lake",all=TRUE)
fwrite(lakOVata.ansumALL,"D:/DATN/DATA_annualdoysummaries_v6.csv")
lakOVata.ansumALL<-fread("D:/DATN/DATA_annualdoysummaries_v6.csv",stringsAsFactors = T)


summary(lakOVata.ansumALL[lake=="Bubble"])
summary(lakOVata.ansumALL)
lakOVata.ansumALL[,doy_mindif:=doy_min-doy_min.an]
lakOVata.ansumALL[,doy_maxdif:=doy_max-doy_max.an]
lakOVata.ansumALL<-na.omit(lakOVata.ansumALL)
lakOVata.ansumALL[doy_mindif>0,doy_mindif:=0]
lakOVata.ansumALL[doy_maxdif<0,doy_maxdif:=0]
lakOVata.ansumALL[,doy_maxdif:=abs(doy_maxdif)]
lakOVata.ansumALL[,doy_mindif:=abs(doy_mindif)]
lakOVata.ansumALL[,doy_range:=abs(doy_max.an-doy_min.an+1)]

fwrite(lakOVata.ansumALL,"D:/DATN/DATA_annualdoysummaries_v7.csv")
lakOVata.ansumALL<-fread("D:/DATN/DATA_annualdoysummaries_v7.csv",stringsAsFactors = T)

ansum.trends<-lakOVata.ansumALL[,.(doys.trend=zyp.sen(doys~year)$coefficients[2],
                                   doys.coverage.trend=zyp.sen(doy_range~year)$coefficients[2],
                                   meandoys=mean(doys)),.(lake)]

fwrite(ansum.trends,"D:/DATN/DATA_ansum.trends_v1.csv")
ansum.trends<-fread("D:/DATN/DATA_ansum.trends_v1.csv",stringsAsFactors = T)



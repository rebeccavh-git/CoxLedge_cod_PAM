library(plotrix)
library(suncalc)
library(lunar)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(circular)
library(pscl)
library(maptools)
library(rgdal)
library(rgeos)
library(glmmTMB)
library(emmeans)
library(doBy)
library(MASS)
library(plyr)
library(ggspatial)
library(rasterVis)
library(ggsn)
library(grid)
library(itsadug)
library(car)
setwd("C:\\Users\\mdean\\Google Drive\\Cod\\Passive Acoustics\\Mass Bay - Long Term\\")
tz_shift<-0 #LOCAL TIME ACCORDING TO P. CAIGER 12/2/18 
crp_bathy<-colorRampPalette(c("black","white"))

#############
# FUNCTIONS
source("C:\\Users\\mdean\\Google Drive\\R\\micahs_functions.r")
if(TRUE){
  norm<-function(x,type="max1"){
    if(type=="max1")return(x/max(x,na.rm=TRUE))
    if(type=="sum1")return(x/sum(x,na.rm=TRUE))
  }
  
  get_period<-function(dt){
    sunrise<-sunriset(matrix(c(-70.66,42.4),ncol=2),dt,direction="sunrise",POSIXct.out=TRUE)$time-tz_shift
    sunset<-sunriset(matrix(c(-70.66,42.4),ncol=2),dt,direction="sunset",POSIXct.out=TRUE)$time-tz_shift
    dawn<-sunrise-0.5*3600
    morn_end<-sunrise+3*3600
    eve_beg<-sunset-3*3600
    dusk<-sunset+0.5*3600
    period<-rep(NA,length(dt))
    period[dt<sunrise]<-"NIGHT1"
    period[dt>sunset]<-"NIGHT2"
    period[dt>=sunrise&dt<=sunset]<-"DAY"
    return(period)
  }
  
  get_dayhrs<-function(dt){
    sunrise<-sunriset(matrix(c(-70.66,42.68),ncol=2),dt,direction="sunrise",POSIXct.out=TRUE)$day_frac
    sunset<-sunriset(matrix(c(-70.66,42.68),ncol=2),dt,direction="sunset",POSIXct.out=TRUE)$day_frac
    return(24*(sunset-sunrise))
  }
  
  study_year<-function(d){
    d<-as.Date(d)
    yr<-as.numeric(format(d,'%Y'))-2012
    yr<-ifelse(as.numeric(format(d,'%m'))<6,yr-1,yr)
    return(paste("Y",yr,sep=""))
  }
  
  nd_update<-function(nd,var,vals){
    out<-NULL
    n<-length(vals)
    for(i in 1:n)out<-rbind(out,nd)
    out[var]<-vals
    out$Hsin<-sin(2*pi*(out$H/24))
    out$Hcos<-cos(2*pi*(out$H/24))
    out$Jsin<-sin(2*pi*(out$J/365))
    out$Jcos<-cos(2*pi*(out$J/365))
    out$Msin<-sin(out$MOON)
    out$Mcos<-cos(out$MOON)
    out$MOON2<-out$MOON
    out$MOON2[out$MOON2>pi]<-out$MOON2[out$MOON2>pi]-pi
    out$MOON2<-2*out$MOON2
    out$D<-ifelse(out$H<6|out$H>16,"night","day")
    out$Lsin<-sin(out$MOON2) #SEMI-LUNAR
    out$Lcos<-cos(out$MOON2) #SEMI-LUNAR
    return(out)
  }  
  
  make_circ<-function(d){
    out<-d
    namez<-names(d)
    if("H"%in%namez){
      out$Hsin<-sin(2*pi*(out$H/24))
      out$Hcos<-cos(2*pi*(out$H/24))
    }
    if("J"%in%namez){
      out$Jsin<-sin(2*pi*(out$J/365))
      out$Jcos<-cos(2*pi*(out$J/365))
    }
    if("MOON"%in%namez){
      out$Msin<-sin(out$MOON)
      out$Mcos<-cos(out$MOON)
      out$MOON2<-out$MOON
      out$MOON2[out$MOON2>pi]<-out$MOON2[out$MOON2>pi]-pi
      out$MOON2<-2*out$MOON2
      out$Lsin<-sin(out$MOON2) #SEMI-LUNAR
      out$Lcos<-cos(out$MOON2) #SEMI-LUNAR
    }
    return(out)
  }
  
  measure<-function(nclick=2){
    xy<-locator(nclick)
    return(sqrt(diff(xy$x)^2+diff(xy$y)^2))
  }
  
  dfrac<-function(x){
    out<-rep(NA,length(x))
    for(i in 1:length(x)){
      hms<-as.numeric(strsplit(as.character(x[i]),split=":")[[1]])
      out[i]<-hms[1]/24+hms[2]/(24*60)+hms[3]/(24*60*60)
    }
    return(out)
  }
      
  gimme_date <- function(x){
    format(as.Date(x, origin = '1970-01-01'),'%b-%d')
  }
}

################
# GIS 
if(TRUE){
  ll<-CRS("+init=epsg:4326") #NAD83/WGS84  (i.e., LAT & LON)
  masp<-CRS("+init=epsg:26986")  #MASS STATE PLANE - METER
  dem8<-raster("C:\\Users\\mdean\\Google Drive\\cod\\winter cod\\passive acoustics\\for micah\\dem_8m.tif",proj4string=masp)
  gispath<-"C:\\Users\\mdean\\Google Drive\\cod\\ibs2\\analysis\\optimal fishing\\ibs\\shapefiles"
  land2<-readOGR(gispath,"land2")
  bathy<-readOGR("C:\\Users\\mdean\\Google Drive\\GIS\\data","BATHYMGM_ARC")
  land<-readOGR("C:\\Users\\mdean\\Google Drive\\GIS\\data","newyork2mainePoly")
  eez<-readOGR("C:\\Users\\mdean\\Google Drive\\gis\\data","EEZ")
  sbnms<-readOGR("C:\\Users\\mdean\\Google Drive\\gis\\data","sbnms_py")
  xl<-c(-71,-70.25); yl<-c(42,42.75)
  dem<-raster("C:\\Users\\mdean\\Google Drive\\GIS\\data\\crm_all_mass.tif")
  #dem<-raster("W:\\GIS\\gisdata\\images\\Bathymetry\\crm_gom.tif")
  dem_crop<-crop(dem,extent(c(xl,yl)))
  dem_crop<-dem_crop*(dem_crop<0)
}

###############
# PAM DATA
if(TRUE){
  pam_raw<-read.csv("mass_bay_long_term_pam_data.csv")
  pam_raw$Cod<-tolower(pam_raw$Cod)
  pam_raw$Date<-as.Date(pam_raw$Date,tz="UTC")
  pam_raw$Year<-as.numeric(format(pam_raw$Date,'%Y'))
  pam<-subset(pam_raw,Cod=='y') #JUST VERIFIED COD 
  pam$Date<-as.Date(pam$Date,tz="UTC")
  pam$Year<-as.numeric(format(pam$Date,'%Y'))
  pam$Month<-as.numeric(format(pam$Date,'%m'))
  pam$Time1<-pam$Time_Begin/(60*60*24)
  pam$Time2<-pam$Time_End/(60*60*24)
  pam$Clock<-dfrac(pam$Clock_Begin)
  pam$Time1[pam$Year==2016]<-pam$Clock[pam$Year==2016]
  pam$Datetime<-as.POSIXct(pam$Date+pam$Time1)#UTC-5hrs
  pam$H<-floor(pam$Time1*24)
  pam$J<-as.numeric(format(pam$Date,'%j'))
  pam$W<-as.numeric(format(pam$Date,'%W'))
  pam$PERIOD<-get_period(pam$Datetime)
  pam$MOON<-lunar.phase(as.Date(pam$Date),name=FALSE)
  pam$MOON4<-lunar.phase(as.Date(pam$Date),name=4)
  pam$MOON8<-lunar.phase(as.Date(pam$Date),name=8)
  #pam<-subset(pam,Year!=2016)
  
  #SAMPLING EFFORT 
  pam_eff<-read.csv("mass_bay_pam_effort.csv")
  pam_eff$Date<-as.Date(pam_eff$Date)
  hdf<-data.frame(H=0:23)
  mdates<-merge(pam_eff,hdf)
  
  #STATIONS
  maru_g<-read.csv("maru_sites_gateway.csv"); maru_g$Type<-"Gateway"
  maru_w<-read.csv("maru_sites_winter_cod.csv"); maru_w$Type<-"Winter"
  xy<-rbind(maru_g[c("LAT","LON")],maru_w[c("LAT","LON")])
  marus<-rbind(merge(data.frame(Year=2007:2012),maru_g)[names(maru_w)],maru_w)
  marus$DEM8_DEPTH<-raster::extract(dem8,spTransform(SpatialPoints(marus[c("LON","LAT")],proj4string=ll),masp))
  marus$DEPTH<-round(-marus$DEM8_DEPTH)
  maru_ptz<-SpatialPointsDataFrame(marus[c("LON","LAT")],proj4string=ll,data=marus)
}

################
# EXPLORE DATA 
if(FALSE){
  #PLOT OF STATIONS
  xl<-c(-71,-70); yl<-c(42,42.75)
  xt<-extent(spTransform(SpatialPoints(cbind(xl,yl),proj4string=ll),masp))
  xt<-extent(c(xl,yl))
  cxt<-0.75; poz<-3
  crx<-ll
  crp_bathy<-colorRampPalette(c("black","white"))
  par(mar=c(1,1,1,1))
  plot(spTransform(land2,ll),col="transparent",border="transparent",xlim=xt[1:2],ylim=xt[3:4])
  plot(dem_crop,col=crp_bathy(256),add=TRUE,legend=FALSE)
  plot(spTransform(land,ll),col="darkgray",border="transparent",xlim=xt[1:2],ylim=xt[3:4],add=TRUE)
  lines(spTransform(subset(bathy,CONTOUR%in%c(-50)),crx),col="blue")
  plot(spTransform(maru_ptz,crx),add=TRUE,cex=0.5,pch=19,col=ifelse(maru_ptz$Type=="Winter","red","black"))
  text(spTransform(maru_ptz,crx),labels=maru_ptz$Site,pos=poz,cex=cxt,col=ifelse(maru_ptz$Type=="Winter","red","black"))
  par(mar=c(5,4,4,1))
  
  #WHAT IS THE AVG DIST MOVED FOR SITES THAT WERE AGGREGATED?
  x<-spTransform(maru_ptz,masp)
  x$X<-coordinates(x)[,1]
  x$Y<-coordinates(x)[,2]
  ag<-unique(x@data[c("Site","X","Y")])
  ag<-ag[order(ag$Site),]
  movers<-as.numeric(names(table(ag$Site)[table(ag$Site)>1]))
  for(i in movers){
    sub<-subset(ag,Site==i)
    xd<-diff(sub$X)
    yd<-diff(sub$Y)
    dist<-sqrt((xd^2)+(yd^2))
    df<-data.frame(Site=i,Dist=round(dist))
    if(i==movers[1])out<-df
    if(i!=movers[1])out<-rbind(out,df)
  }
  mean(out$Dist)
  
  #WHAT IS THE AVG BOTTEMP WITHIN A MONTH ACROSS THE ARRAY
  x2<-ptzfun(summaryBy(cbind(LAT,LON)~Site,data=maru_ptz@data,FUN=mean,keep.names=TRUE))
  plot(fvbt$Dec)
  plot(x2,add=TRUE,pch=1)
  x2$temp<-raster::extract(fvbt$Jan,x2)
  range(x2$temp)
  mean(x2$temp)
  sd(x2$temp)

  #data summary tables
  table(pam[c("Year","Channel")])
  table(marus[c("Year","Channel")])
  table(pam[c("H","Year")])
  table(pam[c("MOON8","Year")])
  barplot(table(pam$MOON8)) #WICKED LUNAR PATTERN! 
  barplot(table(pam$H)) #
}

##########################
# CIRCULAR SUMMARY PLOTS
if(TRUE){
  
  #######################
  #CIRCULAR BARPLOTS - GENERAL
  
    #CIRCULAR HOUR
    if(TRUE){
      dat<-aggregate(Cod~H,data=pam,FUN=length)
      dat$period<-ifelse(dat$H<6|dat$H>16,"night","day")
      dat$period[dat$H%in%c(6,16)]<-"sunrise/set"
      colz_diel<-c("lightgray","black","darkgray")
      p_diel<-ggplot(dat, aes(x=as.factor(H), y=Cod,fill=period)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        geom_bar(stat="identity", alpha=0.85,show.legend=FALSE) +
        scale_fill_manual(values=colz_diel)+
        # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
        ylim(-max(dat$Cod)*0.5,max(dat$Cod)*1.1)+
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust=0.5)
        ) +
        coord_polar(start = 0)+ # This makes the coordinate polar instead of cartesian.
        ggtitle("Diel")
    }
  
    #CIRCULAR - MOON
    nb<-12
    if(TRUE){
      mbinz<-seq(0,2*pi,length=nb)
      pam$MBIN<-cut(pam$MOON,mbinz)
      dat<-aggregate(Cod~MBIN,data=pam,FUN=length)
      dat$X<-dat$MBIN; dat$VALUE<-dat$Cod
      dat$LABEL<-""
      dat$LABEL[1]<-"New"
      dat$LABEL[nb/2]<-"Full"
      label_data<-dat
      angle=90-360*(as.numeric(label_data$MBIN)-1)/nrow(label_data)     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
      label_data$id<-as.numeric(label_data$MBIN)
      label_data$hjust<-ifelse( angle < -90, 1, 0)
      label_data$angle<-ifelse(angle < -90, angle+180, angle)
      label_data<-subset(label_data,LABEL!="")
      crp<-colorRampPalette(c("black","lightgray","black"))
      colz_lunar<-crp(nb-1)
      p_lunar<-ggplot(dat, aes(x=as.factor(X), y=VALUE,fill=MBIN)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        # This add the bars with a blue color
        geom_bar(stat="identity", alpha=1, show.legend=FALSE) +
        scale_fill_manual(values=colz_lunar)+
        # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
        ylim(-max(dat$VALUE)*0.5,max(dat$VALUE)*1.1)+
        # Custom the theme: no axis title and no cartesian grid
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust=0.5)
        ) +
        #scale_x_discrete(breaks=dat$X[labs$breaks],labels=labs$labels) +
        # This makes the coordinate polar instead of cartesian.
        coord_polar(start = 0)+
        ggtitle("Lunar")+
        geom_text(data=label_data, aes(x=id, y=VALUE*1.2, label=LABEL, hjust=0.5), color="black", fontface="bold", size=5, inherit.aes = FALSE ) 
  
    }
    
    #CIRCULAR ANNUAL
    if(TRUE){
      pam$JBIN<-as.numeric(format(pam$Date,'%W')); nbinz<-52
      #pam$JBIN<-as.numeric(pam$Month); nbinz<-12
      dat<-aggregate(Cod~JBIN,data=pam,FUN=length)
      blanx<-data.frame(JBIN=1:nbinz)
      dat<-merge(dat,blanx,all=TRUE)
      dat<-dat[order(dat$JBIN),]
      dat$Cod[is.na(dat$Cod)]<-0
      dat$X<-dat$JBIN; dat$VALUE<-dat$Cod
      dat$LABEL<-""
      label_data<-dat
      angle=90-360*(as.numeric(label_data$JBIN)-1)/nrow(label_data)     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
      label_data$id<-as.numeric(label_data$JBIN)
      label_data$hjust<-ifelse( angle < -90, 1, 0)
      label_data$angle<-ifelse(angle < -90, angle+180, angle)
      label_data<-subset(label_data,LABEL!="")
      crp<-colorRampPalette(c("black","black","black"))
      colz<-crp(nrow(dat))
      p_annual<-ggplot(dat, aes(x=as.factor(X), y=VALUE)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        # This add the bars with a blue color
        geom_bar(stat="identity", alpha=1) +
        scale_fill_manual(values=colz)+
        # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
        ylim(-max(dat$VALUE)*0.5,max(dat$VALUE)*1.1)+
        # Custom the theme: no axis title and no cartesian grid
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust=0.5)
        ) +
        #scale_x_discrete(breaks=dat$X[labs$breaks],labels=labs$labels) +
        # This makes the coordinate polar instead of cartesian.
        coord_polar(start = 0)+
        ggtitle("Annual")
      #geom_text(data=label_data, aes(x=id, y=VALUE+5000, label=LABEL, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, inherit.aes = FALSE ) 
  
    }
  
  grid.arrange(p_diel,p_lunar,p_annual,ncol=3)
  
  
  
  #######################
  #CIRCULAR HISTOGRAMS!!!
  dat<-subset(pam,Year!=2015)
  dat<-pam
  dat$Y<-as.factor(dat$Year)
  
  #######################
  #CIRCULAR HISTOGRAMS!!! - BY YEAR
  yrz<-sort(unique(dat$Year))
  
  #DIEL
  bw<-7
  ymax_xf<-1.1
  if(TRUE){
    pz_diel<-list()
    for(y in 1:length(yrz)){
      sub<-subset(dat,Year==yrz[y])
      brx<-seq(0,24,by=1)-0.5
      h<-hist(sub$H,breaks=brx,plot=FALSE)
      ymax_d<-max(h$counts)*ymax_xf; ymin_d<-ymax_d/-3; yspan_d<-ymax_d-ymin_d
      xmin_d<-(-0.5); xmax_d<-23.5
      monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      labz_period<-data.frame(H=c(22,10),label=c("Night","Day"))
      #xbrx<-data.frame(H=labz$MOON+0.125)
      labz<-data.frame(H=h$mids,label=h$mids,y=h$counts+yspan_d*0.1)
      labz_period<-merge(labz_period,labz[c("H","y")])
      core<-data.frame(x=c(xmin_d,xmax_d),ymin=c(ymin_d,ymin_d),ymax=c(0,0))
      colz_diel<-c(rep("black",6),"darkgray",rep("lightgray",9),"darkgray",rep("black",7))
      p<-ggplot(sub,aes(x=H))+
        #geom_vline(data=xbrx,aes(xintercept=MOON),col="gray")+
        geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill=colz_diel)+
        coord_polar()+xlim(xmin_d,xmax_d)+ylim(ymin_d,ymax_d)+
        theme_minimal()+
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm"),
          #plot.title = element_text(hjust=0.5,vjust=0.5)
          plot.title = element_blank()
        )+
        #geom_text(data=labz,aes(x=H,y=y,label=label),size=2.5,alpha=0.5)+
        #geom_text(data=labz_period,aes(x=H,y=y+!!yspan_d*0.2,label=label))+
        #annotate("text", x = 0, y = 0, label = "My Ring plot !") +
        #geom_hline(aes(yintercept=0),col="gray")+
        geom_ribbon(data=core,aes(ymin=ymin_d,ymax=0,x=x),col="white",fill="white")+
        geom_text(label=yrz[y], x=0, y=ymin_d,cex=2.5)+
        ggtitle(yrz[y])
      pz_diel[[y]]<-p
    }
    names(pz_diel)<-yrz
  }
  nr<-2; nc<-5
  lmat<-matrix(seq_len(nr*nc),nrow=nr,ncol=nc,byrow=TRUE)
  #ggsave("diel_by_year.png",marrangeGrob(pz_diel,layout_matrix=lmat,top=""),width=10,height=6)
  
  #LUNAR
  nb<-16
  ymax_xf<-1.1
  if(TRUE){
    pz_lunar<-list()
    for(y in 1:length(yrz)){
      sub<-subset(dat,Year==yrz[y])
      brx<-seq(0,2*pi,length=nb+1)
      h<-hist(sub$MOON,breaks=brx,plot=FALSE)
      ymax_l<-max(h$counts)*ymax_xf; ymin_l<-ymax_l/-3; yspan_l<-ymax_l-ymin_l
      xmin_l<-0; xmax_l<-2*pi
      monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      labz<-data.frame(MOON=c(0,0.25,0.5,0.75)*2*pi,label=c("New","Waxing","Full","Waning"))
      xbrx<-data.frame(MOON=labz$MOON+0.125*2*pi)
      core<-data.frame(x=c(xmin_l,xmax_l),ymin=c(ymin_l,ymin_l),ymax=c(0,0))
      crp<-colorRampPalette(c("black","lightgray","black"))
      colz_lunar<-crp(nb)
      p<-ggplot(sub,aes(x=MOON))+
        geom_vline(data=xbrx,aes(xintercept=MOON),col="gray")+
        geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill=colz_lunar)+
        coord_polar()+xlim(xmin_l,xmax_l)+ylim(ymin_l,ymax_l)+
        theme_minimal()+
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm")
          #plot.title = element_text(hjust=0.5,vjust=0.5)
        )+
        #geom_text(data=labz,aes(x=MOON,y=!!ymax_l*0.85,label=label))+
        geom_hline(aes(yintercept=0),col="gray")+
        geom_ribbon(data=core,aes(ymin=!!ymin_l,ymax=0,x=x),col="white",fill="white")+
        geom_text(label=yrz[y], x=0, y=ymin_l,cex=2.5)+
        ggtitle(yrz[y])
      pz_lunar[[y]]<-p
    }
    names(pz_lunar)<-yrz
  }
  nr<-2; nc<-5
  lmat<-matrix(seq_len(nr*nc),nrow=nr,ncol=nc,byrow=TRUE)
  #ggsave("lunar_by_year.png",marrangeGrob(pz_lunar,layout_matrix=lmat,top=""),width=10,height=6)
  
  #ANUUAL
  if(TRUE){
    pz_annual<-list()
    for(y in 1:length(yrz)){
      sub<-subset(dat,Year==yrz[y])
      brx<-seq(0,365,by=bw)
      brx<-unique(c(0,brx,366))
      h<-hist(sub$J,breaks=brx,plot=FALSE)
      ymax_a<-max(h$counts)*1.5; ymin_a<-ymax_a/-3
      xmax_a<-365;xmin_a<-0
      monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      labz<-data.frame(J=seq(15,350,length=12),mon=monz)
      mbrx<-data.frame(J=as.numeric(format(seq(as.Date('2010-01-01'),as.Date('2010-12-01'),by="month"),'%j'))-1)
      core<-data.frame(x=c(xmin_a,xmax_a),ymin=c(ymin_a,ymin_a),ymax=c(0,0))
      nodata<-data.frame(x=c(0,min(sub$J)),ymin=c(0,0),ymax=c(ymax_a,ymax_a))
      p<-ggplot(sub,aes(x=J))+
        geom_ribbon(data=nodata,aes(ymin=!!ymin_a,ymax=!!ymax_a,x=x),fill="gray90")+
        geom_vline(data=mbrx,aes(xintercept=J),col="gray")+
        geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill="black")+
        coord_polar()+xlim(xmin_a,xmax_a)+ylim(ymin_a,ymax_a)+
        theme_minimal()+
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust=0.5)
        )+
        geom_text(data=labz,aes(x=J,y=!!ymax_a,label=mon))+
        geom_ribbon(data=core,aes(ymin=!!ymin_a,ymax=0,x=x),col="white",fill="white")+
        geom_hline(aes(yintercept=0),col="gray")+
        ggtitle(yrz[y])
      pz_annual[[y]]<-p
    }
    names(pz_annual)<-yrz
  }
  nr<-2; nc<-5
  lmat<-matrix(seq_len(nr*nc),nrow=nr,ncol=nc,byrow=TRUE)
  #ggsave("annual_by_year.png",marrangeGrob(pz_annual,layout_matrix=lmat,top=""),width=10,height=6)

  ##################################
  #CIRCULAR HISTOGRAMS!!! - OVERALL
  
  #ANNUAL
  bw<-7
  if(TRUE){
    brx<-seq(0,365,by=bw)
    brx<-unique(c(0,brx,366))
    h<-hist(dat$J,breaks=brx,plot=FALSE)
    ymax_a<-max(h$counts)*1.5; ymin_a<-ymax_a/-3
    xmax_a<-365;xmin_a<-0
    monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    labz<-data.frame(J=seq(15,350,length=12),mon=monz)
    mbrx<-data.frame(J=as.numeric(format(seq(as.Date('2010-01-01'),as.Date('2010-12-01'),by="month"),'%j'))-1)
    core<-data.frame(x=c(xmin_a,xmax_a),ymin=c(ymin_a,ymin_a),ymax=c(0,0))
    nodata<-data.frame(x=c(0,min(dat$J)),ymin=c(0,0),ymax=c(ymax_a,ymax_a))
    p_annual<-ggplot(dat,aes(x=J))+
      geom_ribbon(data=nodata,aes(ymin=ymin_a,ymax=ymax_a,x=x),fill="gray90")+
      geom_vline(data=mbrx,aes(xintercept=J),col="gray")+
      geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill="black")+
      coord_polar()+xlim(xmin_a,xmax_a)+ylim(ymin_a,ymax_a)+
      theme_minimal()+
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"),
        plot.title = element_text(hjust=0.5)
      )+
      geom_text(data=labz,aes(x=J,y=ymax_a,label=mon))+
      geom_ribbon(data=core,aes(ymin=ymin_a,ymax=0,x=x),col="white",fill="white")+
      geom_hline(aes(yintercept=0),col="gray")+
      ggtitle("Annual")
  }
  #p_annual
  
  #LUNAR
  nb<-16
  if(TRUE){
    brx<-seq(0,2*pi,length=nb+1)
    h<-hist(dat$MOON,breaks=brx,plot=FALSE)
    ymax_l<-max(h$counts)*1.5; ymin_l<-ymax_l/-3
    xmin_l<-0; xmax_l<-2*pi
    monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    labz<-data.frame(MOON=c(0,0.25,0.5,0.75)*2*pi,label=c("New","Waxing","Full","Waning"))
    xbrx<-data.frame(MOON=labz$MOON+0.125*2*pi)
    core<-data.frame(x=c(xmin_l,xmax_l),ymin=c(ymin_l,ymin_l),ymax=c(0,0))
    crp<-colorRampPalette(c("black","lightgray","black"))
    colz_lunar<-crp(nb)
    p_lunar<-ggplot(dat,aes(x=MOON))+
      geom_vline(data=xbrx,aes(xintercept=MOON),col="gray")+
      geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill=colz_lunar)+
      coord_polar()+xlim(xmin_l,xmax_l)+ylim(ymin_l,ymax_l)+
      theme_minimal()+
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"),
        plot.title = element_text(hjust=0)
      )+
      geom_text(data=labz,aes(x=MOON,y=ymax_l*0.9,label=label))+
      geom_hline(aes(yintercept=0),col="gray")+
      geom_ribbon(data=core,aes(ymin=ymin_l,ymax=0,x=x),col="white",fill="white")+
      geom_text(label="All\nYears", x=0, y=ymin_l,cex=3)+
      ggtitle("b) Lunar")
  }
  #p_lunar
  
  #DIEL
  if(TRUE){
    brx<-seq(0,24,by=1)-0.5
    h<-hist(dat$H,breaks=brx,plot=FALSE)
    ymax_d<-max(h$counts)*1.5; ymin_d<-ymax_d/-3; yspan_d<-ymax_d-ymin_d
    xmin_d<-(-0.5); xmax_d<-23.5
    monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    labz_period<-data.frame(H=c(22,10),label=c("Night","Day"))
    #xbrx<-data.frame(H=labz$MOON+0.125)
    labz<-data.frame(H=h$mids,label=h$mids,y=h$counts+yspan_d*0.1)
    labz_period<-merge(labz_period,labz[c("H","y")])
    labz_period$y<-c(ymax_d*0.75,ymax_d)
    core<-data.frame(x=c(xmin_d,xmax_d),ymin=c(ymin_d,ymin_d),ymax=c(0,0))
    colz_diel<-c(rep("black",6),"darkgray",rep("lightgray",9),"darkgray",rep("black",7))
    p_diel<-ggplot(dat,aes(x=H))+
      #geom_vline(data=xbrx,aes(xintercept=MOON),col="gray")+
      geom_histogram(breaks=brx,show.legend=FALSE,col="white",fill=colz_diel)+
      coord_polar()+xlim(xmin_d,xmax_d)+ylim(ymin_d,ymax_d)+
      theme_minimal()+
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0)
      )+
      geom_text(data=labz,aes(x=H,y=y,label=label),size=3,alpha=0.75,col="darkgray")+
      geom_text(data=labz_period,aes(x=H,y=y,label=label))+
      geom_hline(aes(yintercept=0),col="gray")+
      geom_ribbon(data=core,aes(ymin=ymin_d,ymax=0,x=x),col="white",fill="white")+
      geom_text(label="All\nYears", x=0, y=ymin_d,cex=3)+
      ggtitle("a) Diel")
  }
  #p_diel
  
  #ggsave("cycles_full_not2015.png",plot=grid.arrange(p_diel,p_lunar,p_annual,ncol=3),width=10,height=4)
  
}  
  
#FIGURE FOR PUBLICATION 
#FULL RAW OBS DIEL LUNAR PLOT BY YEAR & WITH OVERALL PATTERNS AT LEFT 
lmat<-matrix(c(21,21,1:5,21,21,6:10,22,22,11:15,22,22,16:20),nr=4,byrow=TRUE)
pz_all<-c(pz_diel,pz_lunar)
pz_all$p_diel<-p_diel
pz_all$p_lunar<-p_lunar
#ggsave("diel_annual_patterns.jpg",marrangeGrob(pz_all,layout_matrix=lmat,top=""),width=9,height=6,dpi=600)
ggsave("diel_annual_patterns.pdf",marrangeGrob(pz_all,layout_matrix=lmat,top=""),width=9,height=6,dpi=600,colormodel="cmyk")


######################
#MAKE F-ING FACIST FACETS WORK
if(FALSE){
  brx<-seq(0,24,by=1)-0.5
  h<-hist(dat$H,breaks=brx,plot=FALSE)
  ymax_d<-max(h$counts)*1.5; ymin_d<-ymax_d/-3; yspan_d<-ymax_d-ymin_d
  xmin_d<-(-0.5); xmax_d<-23.5
  monz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  labz_period<-data.frame(H=c(22,10),label=c("Night","Day"))
  labz<-data.frame(H=h$mids,label=h$mids,y=h$counts+yspan_d*0.1)
  labz_period<-merge(labz_period,labz[c("H","y")])
  core<-data.frame(x=c(xmin_d,xmax_d),ymin=c(ymin_d,ymin_d),ymax=c(0,0))
  colz_diel<-c(rep("black",6),"darkgray",rep("lightgray",9),"darkgray",rep("black",7))
  p<-ggplot(dat,aes(x=H))+
    geom_histogram(breaks=brx,show.legend=FALSE,col="white")
    
  p+facet_grid(Year~.,scales="free")
}

##############################
# GLMs 
if(TRUE){
  
  #SUMMARIZE GRUNT-RATE (GPH) BY DATE & PERIOD
  if(TRUE){
    agd<-aggregate(View~Year+Month+Date+J+H+Channel,data=pam,FUN=length)
    agd$GPH<-agd$View
    
    #GRUNT PRESENCE
    gruntz<-merge(mdates,agd,all.x=TRUE)
    gruntz$GRUNTS<-!is.na(gruntz$GPH)
    gruntz$GPH[is.na(gruntz$GPH)]<-0
    gruntz$Year<-as.numeric(format(gruntz$Date,'%Y'))
    gruntz$MONTH<-as.numeric(format(gruntz$Date,'%m'))
    gruntz$J<-as.numeric(format(gruntz$Date,'%j'))
    gruntz$PERIOD<-ifelse(gruntz$H<6|gruntz$H>16,"night","day")
    gruntz$period[gruntz$H%in%c(6,16)]<-"sunrise/set"
    gruntz$MOON<-lunar.phase(as.Date(gruntz$Date),name=FALSE)
    gruntz$MOON4<-lunar.phase(as.Date(gruntz$Date),name=4)
    gruntz$MOON8<-lunar.phase(as.Date(gruntz$Date),name=8)
    gruntz$W<-as.numeric(format(gruntz$Date,'%W')) 
    gruntz<-merge(gruntz,marus)
    
    ag<-aggregate(GPH~DEPTH,data=gruntz,FUN=mean)
    plot(GPH~DEPTH,data=gruntz)
    plot(GPH~DEPTH,data=ag,type='o')
    
    gsub<-gruntz
    ag<-aggregate(GPH~J+Year,data=gsub,FUN=mean)
    ag$Dequiv<-ag$J+as.Date('2019-12-31')
    plot(GPH~Dequiv,data=ag,type='o')
    
  }
  
  #SET UP DATA
  if(TRUE){
    gruntz$Hsin<-sin(2*pi*(gruntz$H/24))
    gruntz$Hcos<-cos(2*pi*(gruntz$H/24))
    gruntz$Jsin<-sin(2*pi*(gruntz$J/365))
    gruntz$Jcos<-cos(2*pi*(gruntz$J/365))
    gruntz$D<-ifelse(gruntz$PERIOD=="day","day","night") #SUNRISE/SET = NIGHT 
    gruntz$Msin<-sin(gruntz$MOON)
    gruntz$Mcos<-cos(gruntz$MOON)
    gruntz$MOON2<-gruntz$MOON
    gruntz$MOON2[gruntz$MOON2>pi]<-gruntz$MOON2[gruntz$MOON2>pi]-pi
    gruntz$MOON2<-2*gruntz$MOON2
    gruntz$Lsin<-sin(gruntz$MOON2) #SEMI-LUNAR
    gruntz$Lcos<-cos(gruntz$MOON2) #SEMI-LUNAR
    gruntz$Site<-as.factor(gruntz$Site)
    gruntz$Y<-as.factor(gruntz$Year)
  }
  
  #FUNTIONS TO ALLOW FOR SHORTHAND REFERENCES TO CIRCULAR VARIABLES (Handy!!)
  #CREATE A LIST OF FORMULAE FOR ALL POSSIBLE VARIABLE COMBINATIONS
  if(TRUE){
    fun<-function(x){
      form<-paste(unlist(x),collapse="+")
      form<-gsub("S","Site",form)
      form<-gsub("H","Hsin+Hcos",form)
      form<-gsub("L1","Msin+Mcos",form)
      form<-gsub("L2","Lsin+Lcos",form)
      form<-gsub("J","Jsin+Jcos",form)
      form<-gsub("X1","Jsin:Site+Jcos:Site",form)
      form<-gsub("X2","Jsin:Y+Jcos:Y",form)
      form<-gsub("X3","D:Msin+D:Mcos",form)
      form<-gsub("X4","Y:Site",form)
      form<-gsub("rW","(1|WK)",form)
      form<-gsub("rD","(1|DAY)",form)
      form<-as.formula(paste("~",form,sep=""))
      return(form)
    }
    unfun<-function(form){
      y<-as.character(form)
      x<-y[length(y)]
      x<-gsub("Site","S",x)
      x<-gsub("Hsin","H",x)
      x<-gsub("Hcos","H",x)
      x<-gsub("Msin","L1",x)
      x<-gsub("Mcos","L1",x)
      x<-gsub("Jsin","J",x)
      x<-gsub("Jcos","J",x)
      x<-gsub("Lsin","L2",x)
      x<-gsub("Lcos","L2",x)
      x<-gsub("\\(1 \\| WK\\)","rW",x)
      x<-gsub("\\(1 \\| DAY\\)","rD",x)
      vs<-c("H","L1","J","L2","S:J","Y:J","D:M")
      for(v in vs){
        vv<-paste(v,"+",v)
        x<-gsub(vv,v,x,fixed=TRUE)  
      }
      return(x)
    }
    fun_short<-function(x){
      form<-paste(unlist(x),collapse="+")
      form<-paste("~",form,sep="")
      return(form)
    }
    vgroups<-c("Y","S","H","M","L","J","X1","X3")
    flist<-list()
    for(i in 4:length(vgroups)){
      flist<-c(flist,combn(vgroups,i,simplify=FALSE))
    }
    forms<-lapply(flist,fun)
    forms_short<-lapply(flist,fun_short)
  }
  
  #ESTABLISH DATASET FOR MODELLING
  if(TRUE){
    agn<-aggregate(GPH~Site+Y,data=gruntz,FUN=sum);names(agn)[ncol(agn)]<-"ngrunts"
    agd<-aggregate(J~Site+Y,data=gruntz,FUN=function(x){length(unique(x))})
    agt<-aggregate(H~Site+Y,data=gruntz,FUN=length);names(agt)[ncol(agt)]<-"Htot"
    agp<-aggregate(H~Site+Y,data=gruntz,subset=GPH>0,FUN=length);names(agp)[ncol(agp)]<-"Hpos"
    ag<-merge(merge(merge(agt,agp,all=TRUE),agn,all=TRUE),agd,all=TRUE)
    ag[is.na(ag)]<-0
    ag$pp<-ag$Hpos/ag$Htot
    ag$SiteY<-paste(ag$Site,ag$Y,sep="_")
    gruntz$SiteY<-paste(gruntz$Site,gruntz$Y,sep="_")
    ag2<-summaryBy(cbind(Htot,Hpos,ngrunts,J)~Site,data=ag,FUN=sum,keep.names=TRUE)
    ag2$pp<-ag2$Hpos/ag2$Htot
    ag2[order(ag2$pp),]
    
    dropz1<-subset(ag,J<10)$SiteY
    dropz2<-subset(ag2,pp<0.02)$Site
    gsub<-subset(gruntz,!Site%in%dropz2&!SiteY%in%dropz1)
    gsub$Site<-droplevels(gsub$Site)
    table(gsub[c("Year","Site")])/24

    sitez<-summaryBy(cbind(LAT,LON,DEPTH)~Site,data=marus,FUN=mean,keep.names=TRUE)
    sitez_ptz<-SpatialPointsDataFrame(sitez[c("LON","LAT")],proj4string=ll,data=sitez)
    ag_site<-merge(ag,sitez)
    agj<-aggregate(J~Site+LAT+LON,data=ag_site,FUN=max)
    agpp<-aggregate(pp~Site+LAT+LON,data=ag_site,FUN=max)
    agjp<-merge(agj,agpp)
    agjp$INMOD<-agjp$Site%in%gsub$Site
    ptz_nocod<-subset(SpatialPointsDataFrame(agjp[c("LON","LAT")],proj4string=ll,data=agjp),!INMOD)
    
    gsub$WK<-(format(gsub$Date,'%Y_%W'))
    gsub$DAY<-(format(gsub$Date,'%Y_%j'))
    gsub$WK<-factor(gsub$WK)
    gsub$DAY<-factor(gsub$DAY)
    
    #A MORE RESTRICTIVE DATASET?
    dropz2<-subset(ag,pp<0.05)$SiteY
    gsub2<-subset(gsub,!SiteY%in%dropz2)
    
  }
  gsub$Site<-droplevels(gsub$Site)
  save(gsub,file="C:/Users/mdean/Google Drive/Cod/Passive Acoustics/Mass Bay - Long Term/gsub.rdat")
  load("C:/Users/mdean/Google Drive/Cod/Passive Acoustics/Mass Bay - Long Term/gsub.rdat")  
  
  #OBSELETE - EARLIER GLMs
  if(FALSE){
    #FIND BEST GLM (PRESENCE)
    vgroups<-c("Y","S","H","L1","L2","J","X1","X4")
    bform<-fun(combn(vgroups,length(vgroups),simplify=FALSE))
    b0<-glm(GRUNTS~1,data=gsub,family=binomial)
    bfull<-update(b0,bform)
    #b_best<-step(object=bfull,direction="backward",k=2)
    #AIC(bfull,b_best)
    #b_best<-bfull
    1-logLik(bfull)/logLik(b0) #R2
    
    #FIND BEST ZI MODEL FOR GRUNT RATE
    if(FALSE){
      dist<-"negbin"
      gsub2<-subset(gsub,Site%in%c(13,3,2))
      vgroups<-c("Y","S","H","L1","L2","J","X1","X4")
      form0<-as.formula("GPH~1")
      form<-fun(combn(vgroups,length(vgroups),simplify=FALSE))
      zfull<-zeroinfl(update(form0,form),data=gsub,dist=dist)
    }
    
    #GLMMTMB VERSION
    glmmTMBControl(optCtrl=list(iter.max=10e3,eval.max=10e3))
    vgroups<-c("Y","S","H","L1","L2","J","X1")
    form0<-as.formula("GPH~1")
    form<-fun(combn(vgroups,length(vgroups),simplify=FALSE))
    zgroups<-c("Y","S","H")
    zform<-fun(combn(zgroups,length(zgroups),simplify=FALSE))
    zfull<-glmmTMB(update(form0,form),ziformula=zform,data=gsub,family=nbinom2())
    AIC(zfull); BIC(zfull)
    #summary(zfull)#
  }
  
  ###############
  #REVIEWER SUGGESTIONS - INCORPORATE RANDOM EFFECT FOR SOME TIME UNIT TO REDUCE SERIAL AUTOCORRELATION

  #CANDIDATE MODEL FORMS? - PROBABLY NOT USEFUL B/C RanFX TAKE TOO LONG
  if(FALSE){
    xY<-c("","+Y")
    xS<-c("+S")
    xH<-c("+H")
    xL<-c("","+L1","+L2","+L1+L2")
    xJ<-c("","+J")
    xX<-c("","+X1")
    xR<-c("","+rW")
    smoosh<-function(x){
      xx<-paste(x,collapse="")
      substr(xx,2,nchar(xx)) #GET RID OF THE LEADING "+"
    }
    xg<-apply(expand.grid(xY,xS,xH,xL,xJ,xX,xR),1,FUN=smoosh)
  }

  #################
  #GLMM - PRESENCE
  bform<-update(as.formula("GRUNTS~1"),fun(c("Y","S","H","L1","L2","J","X1")))
  b0<-glmmTMB(bform,data=gsub,family=binomial,control=gtmbcon)
  b1<-glmmTMB(update(bform,~.-Lsin-Lcos),data=gsub,family=binomial,control=gtmbcon)
  b2<-glmmTMB(update(bform,~.+Y:Site),data=gsub,family=binomial,control=gtmbcon)
  b3<-glmmTMB(update(bform,~.+(1|WK)),data=gsub,family=binomial,control=gtmbcon)
  b4<-glmmTMB(update(bform,~.+(1|WK)-Msin-Mcos),data=gsub,family=binomial,control=gtmbcon)
  b5<-glmmTMB(update(bform,~.+(1|WK)-Lsin-Lcos),data=gsub,family=binomial,control=gtmbcon)
  b6<-glmmTMB(update(bform,~.+(1|WK)-Msin-Mcos-Lsin-Lcos),data=gsub,family=binomial,control=gtmbcon)
  b7<-glmmTMB(update(bform,~.+(1|WK)+Y:Site),data=gsub,family=binomial,control=gtmbcon)
  b8<-glmmTMB(update(bform,~.+(1|WK)-Site:Jsin-Site:Jcos+DEPTH:Jsin+DEPTH:Jcos),data=gsub,family=binomial,control=gtmbcon)
  b9<-glmmTMB(update(bform,~.+(1|WK)-Y),data=gsub,family=binomial,control=gtmbcon)
  b10<-glmmTMB(update(bform,~.+(1|DAY)),data=gsub,family=binomial,control=gtmbcon)
  b11<-glmmTMB(update(bform,~.+(1|DAY)-Lsin-Lcos),data=gsub,family=binomial,control=gtmbcon)
  b12<-glmmTMB(update(bform,~.+(1|DAY)-Msin-Mcos),data=gsub,family=binomial,control=gtmbcon)
  b13<-glmmTMB(update(bform,~.+(1|DAY)-Msin-Mcos-Lsin-Lcos),data=gsub,family=binomial,control=gtmbcon)
  b14<-update(bmodz$b6,formula=~.+DEPTH)
  b15<-update(bmodz$b15,formula=~.-Msin-Mcos)
  b16<-update(bmodz$b15,formula=~.-Msin-Mcos+Lsin+Lcos)
  b17<-update(bmodz$b15,formula=~.+Lsin+Lcos)

  bg1<-gam(GRUNTS~Y+Site+te(DEPTH,J,bs=c("cr","cc"))+s(H,bs="cc")+s(MOON,bs="cc")+s(MOON2,bs="cc")+s(WK,bs="re"),data=gsub,family="binomial")
  bg2<-gam(GRUNTS~Y+Site+s(DEPTH)+te(LAT,LON,J,bs=c("cr","cr","cc"))+s(H,bs="cc")+s(MOON,bs="cc")+s(MOON2,bs="cc")+s(WK,bs="re"),data=gsub,family="binomial")
  bg3<-gam(GRUNTS~s(Y)+Site+s(DEPTH)+te(LAT,LON,J,bs=c("cr","cr","cc"))+s(H,bs="cc")+s(MOON,bs="cc")+s(MOON2,bs="cc")+s(WK,bs="re"),data=gsub,family="binomial")

  #############
  #GLMM - RATE
  gtmbcon<-glmmTMBControl(optCtrl=list(iter.max=10e3,eval.max=10e3))
  gtmbcon2<-glmmTMBControl(optCtrl=list(iter.max=10e3,eval.max=10e3,abs.tol=1e-20))
  form<-update(as.formula("GPH~1"),fun(c("Y","S","H","L1","L2","J","X1")))
  form<-update(form,~.+(1|WK))
  zform<-fun(c("Y","S","H"))
  z_0<-glmmTMB(form,ziformula=zform,data=gsub,family=nbinom2(),control=gtmbcon)
  z_w<-update(z_0,formula=~.+(1|WK))
  z_w1<-update(z_w,formula=~.-Msin-Mcos)
  z_w2<-update(z_w,formula=~.-Lsin-Lcos)
  z_w3<-update(z_w,formula=~.-Msin-Mcos-Lsin-Lcos)
  z_w_z1<-update(z_w,ziformula=update(zform,~.+Jsin+Jcos))
  z_w_z2<-update(z_w,ziformula=update(zform,~.+Jsin:Site+Jcos:Site))
  system.time(z_w_z3<-update(z_w,ziformula=update(zform,~.+Jsin+Jcos+Msin+Mcos)))
  system.time(z_w_z4<-update(z_w,ziformula=update(zform,~.+Jsin+Jcos+(1|WK))))
  system.time(z_w_z5<-update(z_w,ziformula=update(zform,~.+(1|WK))))
  system.time(z_w_z6<-update(z_w1,ziformula=update(zform,~.+Jsin+Jcos)))
  system.time(z_w_z7<-update(z_w1,ziformula=update(zform,~.+Jsin+Jcos+Lsin+Lcos)))
  system.time(z_w_z8<-update(z_w1,ziformula=update(zform,~.+Jsin+Jcos+Lsin+Lcos+(1|WK))))  
  system.time(z_w_z9<-update(z_w1,ziformula=update(zform,~.+Jsin+Jcos+Lsin+Lcos+(1|WK)+Site:Jsin+Site:Jcos)))
  system.time(z_w_z10<-update(z_w1,ziformula=update(zform,~.+Site:Y)))
  system.time(z_w_z11<-update(z_w_z9,formula=~.+Y:Site,ziformula=~.+Y:Site))
  system.time(z_w_z12<-update(z_w_z9,formula=~.+Y:Site+Msin+Mcos,ziformula=~.+Y:Site+Msin+Mcos))
  system.time(z_w_z13<-update(z_w_z9,formula=~.-Site:Jsin-Site:Jcos+DEPTH:Jsin+DEPTH:Jcos,ziformula=~.-Site:Jsin-Site:Jcos+DEPTH:Jsin+DEPTH:Jcos))
  system.time(z_w_z14<-update(z_w,formula=~.-Site:Jsin-Site:Jcos+DEPTH:Jsin+DEPTH:Jcos,ziformula=~.+DEPTH:Jsin+DEPTH:Jcos))
  print(Sys.time());system.time(z_w_z18<-glmmTMB(update(form,~.+Y:Site),ziformula=update(zform,~.+DEPTH:Jsin+DEPTH:Jcos+Y:Site),data=gsub,control=gtmbcon))
  print(Sys.time());system.time(z_w_z16<-glmmTMB(form,ziformula=update(zform,~.+DEPTH:Jsin+DEPTH:Jcos+Y:Site),data=gsub,control=gtmbcon))
  print(Sys.time());system.time(z_w_z17<-glmmTMB(form,ziformula=update(zform,~.+Jsin+Jcos+Lsin+Lcos+(1|WK)+Site:Jsin+Site:Jcos+Y:Site),data=gsub,control=gtmbcon))
  print(Sys.time());system.time(zx<-update(zmodz$z12,ziformula=~.-Lsin-Lcos))
  print(Sys.time());system.time(zx2<-update(zmodz$z12,formula=~.-Lsin-Lcos))
  print(Sys.time());system.time(zx3<-update(zmodz$z16,formula=~.+Msin+Mcos))
  print(Sys.time());system.time(zx4<-update(zmodz$z18,ziformula=~.+Msin+Mcos))
  print(Sys.time());system.time(zx5<-update(zmodz$z18,ziformula=~.+DEPTH))
  print(Sys.time());system.time(zx6<-update(zmodz$z18,formula=~.+DEPTH))
  print(Sys.time());system.time(zx7<-update(zmodz$z18,formula=~.+DEPTH,ziformula=~.))
  print(Sys.time());system.time(zx8<-update(zmodz$z22,formula=~.,ziformula=~.+Lsin+Lcos))
  print(Sys.time());system.time(zx9<-update(zmodz$z22,formula=~.-Msin-Mcos,ziformula=~.))
  print(Sys.time());system.time(zx6<-update(zmodz$z22,ziformula=~.-DEPTH))
  print(Sys.time());system.time(zxx<-update(zmodz$z22,ziformula=~.-Site:Jsin-Site:Jcos))
  print(Sys.time());system.time(zxx<-update(zmodz$z22,ziformula=~.-Y-Lsin-Lcos-DEPTH))
  print(Sys.time());system.time(zxx2<-update(zmodz$z22,ziformula=~.-Y-Lsin-Lcos-DEPTH-H))
  print(Sys.time());system.time(zxx3<-update(zmodz$z22,ziformula=~.-Y-Lsin-Lcos))
  print(Sys.time());system.time(zxx4<-update(zmodz$z22,ziformula=~.-Y-DEPTH))
  print(Sys.time());system.time(zxx5<-update(zmodz$z22,ziformula=~.-DEPTH-Lsin-Lcos))
  
  print(Sys.time());system.time(zd0<-update(zmodz$z22,data=gsub2))
  print(Sys.time());system.time(zd1<-update(zmodz$z18,data=gsub2))
  print(Sys.time());system.time(zd2<-update(zmodz$z20,data=gsub2))
  print(Sys.time());system.time(zd3<-update(zmodz$z16,data=gsub2))
  print(Sys.time());system.time(zd4<-update(zmodz$z12,data=gsub2))
  print(Sys.time());system.time(zd5<-update(zmodz$z17,data=gsub2))
  
  
  
  #z_d<-update(z_0,formula=~.+(1|DAY))
  #z_d1<-update(z_d,formula=~.-Msin-Mcos)
  #z_d2<-update(z_d,formula=~.-Lsin-Lcos)
  #z_d3<-update(z_d,formula=~.-Msin-Mcos-Lsin-Lcos)
  #z_d_z1<-update(z_d,ziformula=~.+Jsin+Jcos)
  #z_d_z2<-update(z_d,ziformula=~.+Site:Jsin+Site:Jcos)
  

  ###################
  #STORE & SAVE MODELS
  
  fn<-"C:\\Users\\mdean\\Google Drive\\Cod\\Passive Acoustics\\Mass Bay - Long Term\\pam_mixed_modz.rdat"
  save(gsub,bmodz,zmodz,file=fn)
  load(fn)
  
  
  #bmodz<-list(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13)
  #names(bmodz)<-lapply(bmodz,FUN=function(x){unfun(formula(x))})
  #bmodz[[length(bmodz)+1]]<-b17
  names(bmodz)<-paste("b",1:length(bmodz),sep="")
  btab<-ictab(bmodz)
  btab$mod<-names(bmodz)
  btab$form<-unlist(lapply(bmodz,function(x){unfun(x[["modelInfo"]]$allForm$formula)}))
  btab$loglik<-unlist(lapply(bmodz,logLik))
  btab$dup<-duplicated(btab[,-1])
  btab$CAND<-grepl("rW",btab$form)&!is.na(btab$AIC)&!btab$dup #A CANDIDATE MODEL FORM?
  btab$dAIC[btab$CAND]<-btab$AIC[btab$CAND]-min(btab$AIC[btab$CAND],na.rm=TRUE)
  btab$w<-round(exp(-0.5*btab$dAIC)/sum(exp(-0.5*btab$dAIC),na.rm=TRUE),4)
  btab$PRED<-btab$w>0.1&!is.na(btab$w)
  btab<-btab[order(-btab$CAND,btab$AIC),][c("mod","form","loglik","AIC","dAIC","w","df","PRED")]
  btab
  write.csv(btab,file="bmod_aic.csv")
  
  #zmodz<-list(z_0,z_w,z_w1,z_w2,z_w3,z_w_z1,z_w_z2,z_w_z3,z_w_z4,z_w_z5,z_w_z6,z_w_z9,z_w_z10,z_w_z11,z_w_z12,z_w_z13,z_w_z14)
  #names(zmodz)<-lapply(zmodz,FUN=function(x){unfun(formula(x))})
  #zmodz[[length(zmodz)+1]]<-zxx5
  names(zmodz)<-paste("z",1:length(zmodz),sep="")
  ztab<-ictab(zmodz)
  ztab$mod<-names(zmodz)
  ztab$form<-unlist(lapply(zmodz,function(x){unfun(x[["modelInfo"]]$allForm$formula)}))
  ztab$ziform<-unlist(lapply(zmodz,function(x){unfun(x[["modelInfo"]]$allForm$ziformula)}))
  ztab$loglik<-unlist(lapply(zmodz,logLik))
  ztab$dup<-duplicated(ztab[,-1])
  ztab$CAND<-grepl("rW",ztab$form)&!grepl("DEPTH",ztab$form)&!is.na(ztab$AIC)&!ztab$dup
  ztab$dAIC[ztab$CAND]<-ztab$AIC[ztab$CAND]-min(ztab$AIC[ztab$CAND],na.rm=TRUE)
  ztab$w<-round(exp(-0.5*ztab$dAIC)/sum(exp(-0.5*ztab$dAIC),na.rm=TRUE),4)
  ztab$PRED<-ztab$w>0.1&!is.na(ztab$w)
  ztab<-ztab[order(-ztab$CAND,ztab$AIC),][c("mod","form","ziform","loglik","AIC","dAIC","w","df","PRED")]
  ztab
  write.csv(ztab,file="zmod_aic.csv")
  


  #SELECT BEST MODELS FOR PLOTTING 
  bbest<-bmodz$b15
  zbest<-zmodz$z18

  #EXAMINE SERIAL AUTO-CORRELATION FOR A GIVEN MODEL
  if(FALSE){
    mod<-zbest
    acf(resid(mod)) #RAW TOTAL SERIAL AUTOCORRELATION
    mod$frame$resid<-resid(mod)
    ysx<-unique(mod$frame[c("Y","Site")])
    for(i in 1:nrow(ysx)){ 
      print(ysx[i,])
      sub<-subset(mod$frame,Y==ysx$Y[i]&Site==ysx$Site[i])
      acx<-acf(sub$resid,plot=FALSE)
      df<-data.frame(Y=ysx$Y[i],Site=ysx$Site[i],lag=acx$lag,acf=acx$acf)
      if(i==1)out<-df
      if(i>1)out<-rbind(out,df)
    }
    ag<-aggregate(acf~lag,data=out,FUN=mean)
    plot(acf~lag,data=ag,type='h',ylim=c(-0.1,1)) #AVG SERIAL AUTOCORRELATION ACROSS Y-SITEs
    abline(h=0)
    
    #LOOK AT SERIAL CORRELATION OF RANDOM WEEK FX 
    wkz<-unique(bbest$frame$WK)
    yz<-as.numeric(substr(as.character(wkz),1,4))
    nd<-expand.grid(Site=as.character(unique(gsub$Site)),J=335,MOON=pi,H=16,WK=wkz,DAY=NA,DEPTH=50)
    nd$Y<-as.numeric(substr(as.character(nd$WK),1,4))
    nd<-make_circ(nd)
    nd$pG<-predict(bbest,newdata=nd,type='response')
    nd$rG<-predict(zbest,newdata=nd,type='response')
    emm_w<-summaryBy(cbind(pG,rG)~WK,data=nd,FUN=mean,keep.names=TRUE)
    emm_w$Y<-as.numeric(substr(as.character(emm_w$WK),1,4))
    emm_w$W<-as.numeric(substr(as.character(emm_w$WK),6,8))
    plot(NA,ylim=c(0,1),xlim=c(40,53))
    yrz<-unique(emm_w$Y)
    for(y in 1:length(yrz)){
      sub<-subset(emm_w,Y==yrz[y])
      lines(pG~W,data=sub,col=y)
    }
    legend("topright",legend=yrz,lty=1,col=1:length(yrz))
    
    sub<-subset(emm_w,Y==2011)
    acf(sub$pG)
    
    
    #LOOK AT SERIAL CORRELATION OF RANDOM DAY FX 
    dz<-unique(b_d$frame$DAY)
    yz<-as.numeric(substr(as.character(dz),1,4))
    nd<-expand.grid(Site=as.character(unique(gsub$Site)),J=335,MOON=pi,H=16,DAY=dz)
    nd$Y<-as.numeric(substr(as.character(nd$DAY),1,4))
    nd<-make_circ(nd)
    nd$pG<-predict(b_d,newdata=nd,type='response')
    emm_d<-summaryBy(pG~DAY,data=nd,FUN=mean,keep.names=TRUE)
    emm_d$Y<-as.numeric(substr(as.character(emm_d$DAY),1,4))
    emm_d$D<-as.numeric(substr(as.character(emm_d$DAY),6,8))
    plot(NA,ylim=c(0,1),xlim=c(287,365))
    yrz<-unique(emm_d$Y)
    for(y in 1:length(yrz)){
      sub<-subset(emm_d,Y==yrz[y])
      lines(pG~D,data=sub,col=y)
    }
    legend("topright",legend=yrz,lty=1,col=1:length(yrz))
  }
  
}

write.csv(summary(bbest)$coefficients$cond,file="est_pars_bbest.csv")
write.csv(summary(zbest)$coefficients$cond,file="est_pars_zbest.csv")
write.csv(summary(zbest)$coefficients$zi,file="est_pars_zbest_zi.csv")

#####################
#Observed v Predicted 
bbest_mm<-bmodz[btab$mod[btab$PRED]]
zbest_mm<-zmodz[ztab$mod[ztab$PRED]]
bbest_wts<-btab$w[match(names(bbest_mm),btab$mod)]
zbest_wts<-ztab$w[match(names(zbest_mm),ztab$mod)]

#GRUNT PRESENCE
cx<-0.75
cxa<-1
lwx<-2
lwx_obs<-1.25
#jpeg(file="pred_v_obs_PRESENCE.jpg",width=6,height=4,units="in",res=600)
pdf(file="pred_v_obs_PRESENCE.pdf",width=6,height=4,colormodel="cmyk")
if(TRUE){
  clx_fit<-"black"
  #gsub$pred_GRUNTS<-predict(bbest,type='response')
  gsub$pred_GRUNTS<-predict_multimod(bbest_mm,weights=bbest_wts,type='response')
  msub<-unique(gsub[c("Site","Year")])
  sitez<-unique(subset(marus[order(-marus$LAT),],Site%in%as.character(unique(gsub$Site)))$Site)
  #sitez<-as.character(unique(gsub$Site))
  yrz<-sort(as.character(unique(gsub$Year)))
  xl<-c(282,365)
  d1<-as.Date('2015-01-01')
  d2<-as.Date('2015-12-01')
  dseq<-seq(d1,d2,by="month")
  jseq<-j2seq<-as.numeric(format(dseq,'%j'))
  j2seq[j2seq<100]<-j2seq[j2seq<100]+365
  mseq<-format(dseq,'%b')
  OND<-c(10:12)
  par(mfrow=c(length(sitez),length(yrz)),oma=c(3,3.5,3,4.5),mar=c(0,0,0,0),lend=2)
  for(s in 1:length(sitez)){
    sub<-subset(gsub,Site==sitez[s])
    sub$obs_GRUNTS<-as.numeric(sub$GRUNTS)
    ag_obs<-aggregate(obs_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
    ag_pred<-aggregate(pred_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
    max_obs<-max(ag_obs$obs_GRUNTS)
    max_pred<-max(ag_pred$pred_GRUNTS)
    ymax<-max(c(max_obs,max_pred))
    yl<-c(0,ymax)
    for(y in 1:length(yrz)){
      pred<-subset(gsub,Site==sitez[s]&Year==yrz[y])
      if(nrow(pred)==0){
        plot(1,type='n',axes=FALSE,ylim=yl,xlim=xl)
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
        box()
      }
      if(nrow(pred)>1){
        pred$obs_GRUNTS<-as.numeric(pred$GRUNTS)
        ag_obs<-aggregate(obs_GRUNTS~J+MONTH,data=pred,FUN=mean)
        ag_pred<-aggregate(pred_GRUNTS~J+MONTH,data=pred,FUN=mean)
        plot(obs_GRUNTS~J,data=ag_obs[order(ag_obs$J),],type='h',lwd=lwx_obs,col="gray",xlim=xl,yaxt='n',xaxt='n',ylim=yl)
        lines(pred_GRUNTS~J,data=ag_pred[order(ag_pred$J),],col=alpha(clx_fit,0.5),type='l',lwd=lwx)
      }
      if(s==length(sitez))axis(1,at=j2seq,labels=rep("",length(j2seq)),cex.axis=cx)
      if(s==length(sitez))mtext(1,at=j2seq[OND]+15,text=mseq[OND],cex=cx,las=2,adj=1,line=0.5)
      if(s==1)mtext(side=3,yrz[y],line=0.5,cex=cx)
      if(y==1)mtext(side=2,sitez[s],line=0.5,las=2,cex=cx)
      if(y==length(yrz))axis(4,at=c(round(ymax*0.75,2)),las=2,cex.axis=cxa)
    }
  }
  par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(0,0,0,0))
  form<-formula(bbest)
  tx<-paste(unfun(form),"   AIC =",round(AIC(bbest),1),"   BIC =",round(BIC(bbest),1))
  #mtext(side=3,tx,line=3)
  mtext(side=2,"Site Number",line=3,cex=cx)
  mtext(side=4,"Probability of Grunt Occurrence",line=3,cex=cx)
}
mtext(side=3,"a)",adj=0,cex=cx*1.5,line=3)
dev.off()

#GRUNT RATE
if(FALSE){
  cx<-0.75
  cxa<-1
  lwx<-2
  sigdigs<-1
  #jpeg(file="pred_v_obs_RATE_wk.png",width=6,height=4,units="in",res=600)
  pdf(file="pred_v_obs_RATE_wk.pdf",width=6,height=4)
  if(TRUE){
    clx_fit<-"black"
    #gsub$pred_GRUNTS<-predict(zbest,type='response')
    gsub$pred_GRUNTS<-predict_multimod(zmodz[ztab$mod[ztab$PRED]],type='response')
    msub<-unique(gsub[c("Site","Year")])
    sitez<-unique(subset(marus[order(-marus$LAT),],Site%in%as.character(unique(gsub$Site)))$Site)
    #sitez<-as.character(unique(gsub$Site))
    yrz<-sort(as.character(unique(gsub$Year)))
    xl<-c(282,365)
    d1<-as.Date('2015-01-01')
    d2<-as.Date('2015-12-01')
    dseq<-seq(d1,d2,by="month")
    jseq<-j2seq<-as.numeric(format(dseq,'%j'))
    j2seq[j2seq<100]<-j2seq[j2seq<100]+365
    OND<-c(10:12)
    mseq<-format(dseq,'%b')
    par(mfrow=c(length(sitez),length(yrz)),oma=c(3,3.5,3,4.5),mar=c(0,0,0,0),lend=2)
    for(s in 1:length(sitez)){
      sub<-subset(gsub,Site==sitez[s])
      sub$obs_GRUNTS<-sub$GPH
      ag_obs<-aggregate(obs_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
      ag_pred<-aggregate(pred_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
      max_obs<-max(ag_obs$obs_GRUNTS)
      max_pred<-max(ag_pred$pred_GRUNTS)
      ymax<-(max(c(max_obs,max_pred)))
      #ymax<-max_obs
      yl<-c(0,ymax)
      for(y in 1:length(yrz)){
        pred<-subset(gsub,Site==sitez[s]&Year==yrz[y])
        if(nrow(pred)==0){
          plot(1,type='n',axes=FALSE,ylim=yl,xlim=xl)
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
          box()
        }
        if(nrow(pred)>1){
          pred$obs_GRUNTS<-as.numeric(pred$GPH)
          ag_obs<-aggregate(obs_GRUNTS~J+MONTH,data=pred,FUN=mean)
          ag_pred<-aggregate(pred_GRUNTS~J+MONTH,data=pred,FUN=mean)
          plot(obs_GRUNTS~J,data=ag_obs[order(ag_obs$J),],type='h',lwd=3,col="gray",xlim=xl,yaxt='n',xaxt='n',ylim=yl)
          lines(pred_GRUNTS~J,data=ag_pred[order(ag_pred$J),],col=alpha(clx_fit,0.5),type='l',lwd=lwx)
        }
        #if(s==length(sitez))axis(1,at=j2seq,labels=mseq)
        #if(s==1)mtext(side=3,yrz[y],line=0.5,cex=cx)
        #if(y==1)mtext(side=2,sitez[s],line=0.5,las=2,cex=cx)
        #if(y==length(yrz))axis(4,at=c(round(ymax*0.75)),las=2)
        if(s==length(sitez))axis(1,at=j2seq,labels=rep("",length(j2seq)),cex.axis=cx)
        if(s==length(sitez))mtext(1,at=j2seq[OND]+15,text=mseq[OND],cex=cx,las=2,adj=1,line=0.5)
        if(s==1)mtext(side=3,yrz[y],line=0.5,cex=cx)
        if(y==1)mtext(side=2,sitez[s],line=0.5,las=2,cex=cx)
        if(y==length(yrz))axis(4,at=c(signif(ymax*0.75,sigdigs)),las=2,cex.axis=cxa)
      }
    }
    par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(0,0,0,0))
    form<-formula(zfull)
    tx<-paste(unfun(form),"   ",unfun(zform),"   AIC =",round(AIC(zfull),1),"   BIC =",round(BIC(zfull),1))
    #mtext(side=3,tx,line=3,cex=0.85)
    mtext(side=2,"Site Number",line=3,cex=cx)
    mtext(side=4,"Grunts per Hour",line=3,cex=cx)
  }
  mtext(side=3,"b)",adj=0,cex=cx*1.5,line=3)
  dev.off()
}

cx<-0.75
cxa<-1
lwx<-2
lwx_obs<-1.25
sigdigs<-1
#jpeg(file="pred_v_obs_LOGRATE.png",width=6,height=4,units="in",res=600)
pdf(file="pred_v_obs_LOGRATE.pdf",width=6,height=4,colormodel="cmyk")
if(TRUE){
  clx_fit<-"black"
  #gsub$pred_GRUNTS<-predict(zbest,type='response')
  gsub$pred_GRUNTS<-predict_multimod(zbest_mm,weights=zbest_wts,type='response')
  msub<-unique(gsub[c("Site","Year")])
  sitez<-unique(subset(marus[order(-marus$LAT),],Site%in%as.character(unique(gsub$Site)))$Site)
  #sitez<-as.character(unique(gsub$Site))
  yrz<-sort(as.character(unique(gsub$Year)))
  xl<-c(282,365)
  d1<-as.Date('2015-01-01')
  d2<-as.Date('2015-12-01')
  dseq<-seq(d1,d2,by="month")
  jseq<-j2seq<-as.numeric(format(dseq,'%j'))
  j2seq[j2seq<100]<-j2seq[j2seq<100]+365
  OND<-c(10:12)
  mseq<-format(dseq,'%b')
  par(mfrow=c(length(sitez),length(yrz)),oma=c(3,3.5,3,4.5),mar=c(0,0,0,0),lend=2)
  for(s in 1:length(sitez)){
    sub<-subset(gsub,Site==sitez[s])
    sub$obs_GRUNTS<-sub$GPH
    ag_obs<-aggregate(obs_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
    ag_pred<-aggregate(pred_GRUNTS~J+MONTH+Year,data=sub,FUN=mean)
    max_obs<-max(ag_obs$obs_GRUNTS)
    max_pred<-max(ag_pred$pred_GRUNTS)
    ymax<-log(max(c(max_obs,max_pred)+1))
    #ymax<-max_obs
    yl<-c(0,ymax)
    for(y in 1:length(yrz)){
      pred<-subset(gsub,Site==sitez[s]&Year==yrz[y])
      if(nrow(pred)==0){
        plot(1,type='n',axes=FALSE,ylim=yl,xlim=xl)
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightgray")
        box()
      }
      if(nrow(pred)>1){
        pred$obs_GRUNTS<-as.numeric(pred$GPH)
        ag_obs<-aggregate(obs_GRUNTS~J+MONTH,data=pred,FUN=mean)
        ag_pred<-aggregate(pred_GRUNTS~J+MONTH,data=pred,FUN=mean)
        plot(log(obs_GRUNTS+1)~J,data=ag_obs[order(ag_obs$J),],type='h',lwd=lwx_obs,col="gray",xlim=xl,yaxt='n',xaxt='n',ylim=yl)
        lines(log(pred_GRUNTS+1)~J,data=ag_pred[order(ag_pred$J),],col=alpha(clx_fit,0.5),type='l',lwd=lwx)
      }
      #if(s==length(sitez))axis(1,at=j2seq,labels=mseq)
      #if(s==1)mtext(side=3,yrz[y],line=0.5,cex=cx)
      #if(y==1)mtext(side=2,sitez[s],line=0.5,las=2,cex=cx)
      #if(y==length(yrz))axis(4,at=c(round(ymax*0.75)),las=2)
      if(s==length(sitez))axis(1,at=j2seq,labels=rep("",length(j2seq)),cex.axis=cx)
      if(s==length(sitez))mtext(1,at=j2seq[OND]+15,text=mseq[OND],cex=cx,las=2,adj=1,line=0.5)
      if(s==1)mtext(side=3,yrz[y],line=0.5,cex=cx)
      if(y==1)mtext(side=2,sitez[s],line=0.5,las=2,cex=cx)
      if(y==length(yrz))axis(4,at=c(signif(ymax*0.75,sigdigs)),las=2,cex.axis=cxa)
    }
  }
  par(mfrow=c(1,1),mar=c(4,4,4,4),oma=c(0,0,0,0))
  form<-formula(zbest)
  tx<-paste(unfun(form),"   ",unfun(zform),"   AIC =",round(AIC(zbest),1),"   BIC =",round(BIC(zbest),1))
  #mtext(side=3,tx,line=3,cex=0.85)
  mtext(side=2,"Site Number",line=3,cex=cx)
  mtext(side=4,"log(Grunts-per-Hour + 1)",line=3,cex=cx)
}
mtext(side=3,"b)",adj=0,cex=cx*1.5,line=3)
dev.off()


###################
#CIRCULAR VARIABLES FX

#EST MARGINAL MEANS
#MEAN ACROSS ALL YEARS & SITES, BUT AT MIDNIGHT, FULL MOON, NOV23
if(TRUE){

  #DIEL FX
  print("Diel FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=335,MOON=pi,H=0:24,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_h<-summaryBy(cbind(pG,rG)~H,data=nd,FUN=mean,keep.names=TRUE)
  
  #LUNAR FX
  print("Lunar FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=335,MOON=seq(0,2*pi,length=16),H=0,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_l<-summaryBy(cbind(pG,rG)~MOON,data=nd,FUN=mean,keep.names=TRUE) 
  
  #SEASONAL FX - OVERALL
  print("Seasonal FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=275:400,MOON=pi,H=0,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_j<-summaryBy(cbind(pG,rG)~J,data=nd,FUN=mean,keep.names=TRUE)
  emm_j$edate<-as.Date('2010-12-31')+emm_j$J
  
  #SEASONAL x SITE FX
  print("Seasonal X Site FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=282:365,MOON=pi,H=0,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_js<-summaryBy(cbind(pG,rG)~J+Site,data=nd,FUN=mean,keep.names=TRUE)
  emm_js$edate<-as.Date('2010-12-31')+emm_js$J
  sitz<-unique(emm_js$Site)
  peakz<-NULL
  for(s in 1:length(sitz)){
    sub<-subset(emm_js,Site==sitz[s])
    pk_pG<-sub$edate[sub$pG==max(sub$pG)]
    pk_rG<-sub$edate[sub$rG==max(sub$rG)]
    df<-data.frame(Site=sitz[s],peak_pG=pk_pG,peak_rg=pk_rG)
    peakz<-rbind(peakz,df)
  }
  
  #SITE FX
  print("Site FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=282:365,MOON=pi,H=0,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_s<-summaryBy(cbind(pG,rG)~Site,data=nd,FUN=mean,keep.names=TRUE)
  
  #ANNUAL
  print("Annual FX")
  nd<-expand.grid(Y=as.character(unique(gsub$Y)),Site=as.character(unique(gsub$Site)),J=282:365,MOON=pi,H=0,WK=NA,DAY=NA,DEPTH=50)
  nd<-make_circ(nd)
  #nd$pG<-predict(bbest,newdata=nd,type='response')
  #nd$rG<-predict(zbest,newdata=nd,type='response')
  nd$pG<-predict_multimod(bbest_mm,weights=bbest_wts,newdata=nd,type='response')
  nd$rG<-predict_multimod(zbest_mm,weights=zbest_wts,newdata=nd,type='response')
  emm_y<-summaryBy(cbind(pG,rG)~Y,data=nd,FUN=mean,keep.names=TRUE)
  #emm_y<-summaryBy(pG~Y,data=nd,FUN=mean,keep.names=TRUE)
  emm_y$Y<-as.numeric(as.character(emm_y$Y))
  emm_y<-emm_y[order(emm_y$Y),]
  plot(pG~Y,data=emm_y,type='l')
  
  #emm_ys<-summaryBy(cbind(pG,rG)~Y+Site,data=nd,FUN=mean,keep.names=TRUE)
  emm_ys<-summaryBy(pG~Y+Site,data=nd,FUN=mean,keep.names=TRUE)
  emm_ys$Y<-as.numeric(as.character(emm_ys$Y))
  emm_ys<-emm_ys[order(emm_ys$Y),]
  
  sitz<-unique(emm_ys$Site)
  emm_ys$SiteY<-paste(emm_ys$Site,emm_ys$Y,sep="_")
  plot(pG~Y,data=emm_ys,col=Site,type='n')
  for(s in 1:length(sitz)){
    ltx<-ceiling(s/8)
    lines(pG~Y,data=emm_ys,subset=Site==sitz[s],col=alpha(s,0.33),lty=ltx)
    esub<-subset(emm_ys,SiteY%in%gsub$SiteY)
    lines(pG~Y,data=esub,subset=Site==sitz[s],col=s,lwd=2,lty=ltx)
  }
  legend("topleft",legend=sitz,col=1:length(sitz),cex=0.85,lty=ceiling((1:length(sitz))/8),lwd=2)
  
}

#GGPLOT FIGURE of DIEL, LUNAR FX
if(TRUE){
  alf<-0.5
  yfrac<-0.95
  ############

  ###########
  #PRESENCE
  ymn<-min(c(emm_h$pG,emm_l$pG))
  ymx<-max(c(emm_h$pG,emm_l$pG))
  yrange<-ymx-ymn
  ymn<-ymn-(yrange*0.25)
  ymx<-ymx+(yrange*0.25)

  #DIEL PRESENCE
  xmx<-24
  labz_h<-data.frame(H=c(23,5.5,11,16.5),lab=c("Night","","Day",""))
  colz_h<-rep("black",24); colz_h[7:16]<-"white"; colz_h[c(6,17)]<-"darkgray"; 
  pg_circ_h<-ggplot()+
      geom_rect(aes(xmin=0:23,xmax=1:24,ymin=!!ymn,ymax=!!ymx*yfrac),fill=alpha(colz_h,alf))+
      coord_polar()+geom_line(data=emm_h,aes(x=H,y=pG),lwd=2)+
      ylim(ymn,ymx)+theme_minimal()+
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust=0.5,vjust=1)
      ) +
      ggtitle("Grunt Presence")+
      geom_text(data=labz_h, aes(x=H, y=!!ymx, label=lab, hjust=0.5), color="black", size=4 ) 
    
  
  #LUNAR PRESENCE
  seq_l<-seq(0,2*pi,length=32)
  crp_l<-colorRampPalette(c("black","white","black"))
  colz_l<-crp_l(length(seq_l)-1)
  labz_l<-data.frame(MOON=c(0,pi/2,pi,pi*1.5),lab=c("New Moon","","Full Moon",""))
  pg_circ_l<-ggplot()+geom_rect(aes(xmin=seq_l[-length(seq_l)],xmax=seq_l[-1],ymin=!!ymn,ymax=!!ymx*yfrac),fill=alpha(colz_l,alf))+
    coord_polar()+geom_line(data=emm_l,aes(x=MOON,y=pG),lwd=2)+ylim(ymn,ymx)+
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust=0.5,vjust=1)
    ) +
    xlab("Lunar")+
    geom_text(data=labz_l, aes(x=MOON, y=!!ymx, label=lab, hjust=0.5), color="black", size=4 ) 
  
  
  ############
  #RATE
  ymn<-min(c(emm_h$rG,emm_l$rG))
  ymx<-max(c(emm_h$rG,emm_l$rG))
  yrange<-ymx-ymn
  ymn<-ymn-(yrange*0.25)
  ymx<-ymx+(yrange*0.25)
  
  #DIEL RATE
  labz_h<-data.frame(H=c(23,5.5,11,16.5),lab=c("Night","","Day",""))
  colz_h<-rep("black",24); colz_h[7:16]<-"white"; colz_h[c(6,17)]<-"darkgray"; 
  rg_circ_h<-ggplot()+
    geom_rect(aes(xmin=0:23,xmax=1:24,ymin=!!ymn,ymax=!!ymx*yfrac),fill=alpha(colz_h,alf))+
    coord_polar()+geom_line(data=emm_h,aes(H,rG),lwd=2)+ylim(ymn,ymx)+theme_minimal()+
    theme(
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust=0.5,vjust=1)
    ) +
    ggtitle("Grunt Rate")+
    geom_text(data=labz_h, aes(x=H, y=!!ymx, label=lab, hjust=0.5), color="black", size=4 ) 
  
  #LUNAR RATE
  seq_l<-seq(0,2*pi,length=32)
  crp_l<-colorRampPalette(c("black","white","black"))
  colz_l<-crp_l(length(seq_l)-1)
  labz_l<-data.frame(MOON=c(0,pi/2,pi,pi*1.5),lab=c("New Moon","","Full Moon",""))
  rg_circ_l<-ggplot()+geom_rect(aes(xmin=seq_l[-length(seq_l)],xmax=seq_l[-1],ymin=!!ymn,ymax=!!ymx*yfrac),fill=alpha(colz_l,alf))+
    coord_polar()+geom_line(data=emm_l,aes(x=MOON,y=rG),lwd=2)+ylim(ymn,ymx)+
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust=0.5,vjust=1)
    ) +
    geom_text(data=labz_l, aes(x=MOON, y=!!ymx, label=lab, hjust=0.5), color="black", size=4 ) 
  
  grid.arrange(pg_circ_h,rg_circ_h,pg_circ_l,rg_circ_l,nrow=2)
  
  ht<-5 #HEIGHT OF PLOT IN INCHES
  nr<-2; nc<-2; lmat<-matrix(seq_len(nr*nc),nrow=nr,ncol=nc,byrow=TRUE)
  #ggsave("diel_lunar_fx.jpg",marrangeGrob(list(pg_circ_h,rg_circ_h,pg_circ_l,rg_circ_l),layout_matrix=lmat,top="",padding=unit(0,"line")),width=ht*1,height=ht,dpi=600)
  ggsave("diel_lunar_fx.pdf",marrangeGrob(list(pg_circ_h,rg_circ_h,pg_circ_l,rg_circ_l),layout_matrix=lmat,top="",padding=unit(0,"line")),width=ht*1,height=ht,dpi=600)
  
}


#OTHER VARIABLE FX

##########################
# SEASONAL + YEAR FX PLOT FOR MANUSCRIPT 
ymx_xf<-1.25
lwx<-2
ltx<-c(1,4)
clx<-c("black","black")
cxab<-1.75
#jpeg("seasonal_year_fx.jpg",width=5.25,height=5.25*1.4,units="in",res=600)
pdf("seasonal_year_fx.pdf",width=5.25,height=5.25*1.4)
par(bg="white")
par_panel(2,1,mar=c(2.5,4.5,1.5,4.5),oma=c(0,0,0,0))
if(TRUE){
  ############
  #SEASONAL FX
  
  #CALC EMMEANS by J 
  J<-seq(273,365,by=3)
  Jsin<-sin(2*pi*(J/365))
  Jcos<-cos(2*pi*(J/365))
  edate<-as.Date('2010-12-31')+J
  df_j<-data.frame(J,Jsin,Jcos,edate)
  emm_j_pg<-as.data.frame(emmeans(bbest,~Jsin+Jcos,at=list(J=J,Jsin=Jsin,Jcos=Jcos)),type='response')
  emm_j_rg<-as.data.frame(emmeans(zbest,~Jsin+Jcos,at=list(J=J,Jsin=Jsin,Jcos=Jcos)),type='response')
  emm_j_pg<-merge(df_j,emm_j_pg)
  emm_j_rg<-merge(df_j,emm_j_rg)
  
  xl<-as.Date(c('2011-10-01','2012-01-01'))
  sub<-subset(emm_j_pg,edate>=xl[1]&edate<=xl[2])
  plot(prob~edate,data=sub,type='l',col=clx[1],yaxt='n',ylab="",xlab="",lwd=lwx,lty=ltx[1],xlim=xl,ylim=c(0,max(sub$prob)*ymx_xf))
  axis(2,las=2)
  mtext(side=2,"Grunt Presence",line=3.25)
  peak_pG<-sub[sub$prob==max(sub$prob),]
  #axis.Date(3,at=peak_pG$edate)
  par(new=TRUE)
  sub<-subset(emm_j_rg,edate>=xl[1]&edate<=xl[2])
  plot(rate~edate,data=sub,type='l',col=clx[2],yaxt='n',ylab="",xlab="",lwd=lwx,lty=ltx[2],xlim=xl,xaxt='n',ylim=c(0,max(sub$rate)*ymx_xf))
  axis(4,col=clx[2],col.axis=clx[2],las=2)
  mtext(side=4,"Grunt Rate (# / Hr)",line=3,col=clx[2])
  peak_rG<-sub[sub$rate==max(sub$rate),]
  #axis.Date(3,at=peak_rG$edate,col="blue",col.axis="blue")
  legend("topright",c("Grunt Presence","Grunt Rate"),lty=ltx,col=clx,lwd=lwx)
  mtext(side=3,"a)",adj=0.05,line=-2,cex=cxab)
  
  ##########
  #YEAR  FX 
  
  #CALC EMMEANS by Y
  emm_y_pg<-as.data.frame(emmeans(bbest,~Y,type='response'))
  emm_y_rg<-as.data.frame(emmeans(zbest,~Y,type='response'))
  emm_y_pg$Y<-as.numeric(as.character(emm_y_pg$Y))
  emm_y_rg$Y<-as.numeric(as.character(emm_y_rg$Y))
  
  plot(prob~Y,data=emm_y_pg,type='l',col=alpha(clx[1],0.25),yaxt='n',ylab="",xlab="",lwd=lwx,lty=ltx[1],ylim=c(0,max(emm_y_pg$prob)*ymx_xf))
  lines(prob~Y,data=emm_y_pg,type='l',col=clx[1],lwd=lwx,subset=Y<=2012)
  axis(2,las=2)
  mtext(side=2,"Grunt Presence",line=3.25)
  
  par(new=TRUE)
  plot(rate~Y,data=emm_y_rg,type='l',col=alpha(clx[2],0.25),yaxt='n',ylab="",xlab="",lwd=lwx,lty=ltx[2],ylim=c(0,max(emm_y_rg$rate)*ymx_xf),xaxt='n')
  lines(rate~Y,data=emm_y_rg,type='l',col=clx[1],lwd=lwx,lty=ltx[2],subset=Y<=2012)
  axis(4,col=clx[2],col.axis=clx[2],las=2)
  mtext(side=4,"Grunt Rate  (# / Hr)",line=3,col=clx[2])
  par(lend=2)
  legend("topright",c("Grunt Presence","Grunt Rate","Assessment SSB"),lty=c(ltx,1),col=c(clx,alpha("black",0.25)),lwd=c(lwx,lwx,15))
  mtext(side=3,"b)",adj=0.05,line=-2,cex=cxab)
  
  ass<-read.csv("C:\\Users\\mdean\\Google Drive\\Cod\\BREP\\AFS Talk - Reno 2019\\cod_had_assessment_ssb.csv")
  ass$cod_lo<-ass$ssb_cod-2*ass$sd_cod
  ass$cod_hi<-ass$ssb_cod+2*ass$sd_cod
  par(new=TRUE)
  sub_ass<-subset(ass,year%in%(2007:2016))
  plot(ssb_cod~year,col="transparent",data=sub_ass,type='l',lwd=2,ylab="",ylim=c(0,max(sub_ass$ssb_cod)*ymx_xf),axes=FALSE)
  polygon(c(sub_ass$year,rev(sub_ass$year)),c(sub_ass$cod_lo,rev(sub_ass$cod_hi)),col=alpha("black",0.25),border="transparent")

  
  df_cor<-data.frame(year=sub_ass$year,ssb=sub_ass$ssb_cod,pG=emm_y_pg$prob,rG=emm_y_rg$rate)
  df_cor2<-subset(df_cor,year<2013)
  cor.test(df_cor$ssb,df_cor$pG)
  cor.test(df_cor2$ssb,df_cor2$pG)
  
  cor.test(df_cor$ssb,df_cor$rG)
  cor.test(df_cor2$ssb,df_cor2$rG)
  
  
}
par_unpanel()
dev.off()


#YEAR EFFECT FROM EMMEANS w/CIs
pg_trend<-as.data.frame(emmeans(bbest,~Y,type='response'))
rg_trend<-as.data.frame(emmeans(zbest,~Y,type='response'))
x<-pg_trend
plot(prob~Y,data=x)
x$Y<-as.numeric(as.character(x$Y))
sub<-subset(x,Y<=2012)
plot(prob~Y,data=x,type='l',lty=2,xlim=c(2007,2016),ylim=c(0,0.2),col="darkgray")
lines(prob~Y,data=sub,type='l')
polygon(c(sub$Y,rev(sub$Y)),c(sub$lower.CL,rev(sub$upper.CL)),border='transparent',col=alpha("black",0.25))


############
#SITE FX
ptz<-merge(sitez_ptz[c("Site","LAT","LON","DEPTH")],peakz)
ptz<-merge(ptz,emm_s)
ptz$peak_pG<-as.Date(ptz$peak_pG,origin='1970-01-01')
if(TRUE){
  poz<-2
  xt<-extent(-71,-70.25,42,42.5)
  par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(0,0,0,0))
      
  plot(spTransform(land2,ll),col="transparent",border="transparent",xlim=xt[1:2],ylim=xt[3:4])
  plot(dem_crop,col=crp_bathy(256),add=TRUE,legend=FALSE)
  plot(spTransform(land,ll),col="darkgray",border="transparent",xlim=xt[1:2],ylim=xt[3:4],add=TRUE)
  plot(ptz,cex=bubbles(ptz$pG),pch=19,col=alpha("blue",0.5),add=TRUE)
  plot(ptz_nocod,pch=19,cex=0.75,col="red",add=TRUE)
  text(spTransform(ptz,crx),labels=ptz$Site,pos=poz,cex=cxt)
  text(spTransform(ptz_nocod,crx),labels=ptz_nocod$Site,pos=poz,cex=cxt)
  
  plot(spTransform(land2,ll),col="transparent",border="transparent",xlim=xt[1:2],ylim=xt[3:4])
  plot(dem_crop,col=crp_bathy(256),add=TRUE,legend=FALSE)
  plot(spTransform(land,ll),col="darkgray",border="transparent",xlim=xt[1:2],ylim=xt[3:4],add=TRUE)
  plot(ptz,cex=bubbles(ptz$rG),pch=19,col=alpha("blue",0.5),add=TRUE)
  text(spTransform(ptz,crx),labels=ptz$Site,pos=poz,cex=cxt)
  plot(ptz_nocod,pch=19,cex=0.75,col="red",add=TRUE)
  text(spTransform(ptz_nocod,crx),labels=ptz_nocod$Site,pos=poz,cex=cxt)
  
  par(mfrow=c(1,1),mar=c(1,1,1,1),oma=c(0,0,0,0))
}

#OBSELETE
if(FALSE){
  
  #GGPLOT MAP?
  xt<-extent(c(-71,-70.25,42,42.75))
  shore<-crop(spTransform(land,ll),xt)
  shore@data$id = rownames(shore@data)
  shore.points = fortify(shore, region="id")
  shore.df = join(shore.points, shore@data, by="id")
  ptz@data$peak_pG<-as.numeric(ptz@data$peak_pG)
  ptz@data$peak_rG<-as.numeric(ptz@data$peak_rg)
  map_pg<-ggplot() + 
    geom_point(data=ptz@data,aes(x=LON,y=LAT,size=pG,col=peak_pG))+
    scale_color_gradientn(colors=c("red","blue","cyan"),labels=gimme_date)+
    scale_size(range=c(1,10))+
    xlim(xt[1],xt[2])+ylim(xt[3],xt[4])+
    geom_polygon(data=shore.df,aes(long,lat,group=group),fill="darkgray")+
    geom_text(data=ptz@data,aes(x=LON,y=LAT,label=Site))+
    coord_fixed(1.3)+
    theme_minimal()
  map_pg
  
  map_rg<-ggplot() + 
    geom_point(data=ptz@data,aes(x=LON,y=LAT,size=rG,col=peak_rG))+
    scale_color_gradientn(colors=c("red","blue","cyan"),labels=gimme_date)+
    scale_size(range=c(1,10))+
    xlim(xt[1],xt[2])+ylim(xt[3],xt[4])+
    geom_polygon(data=shore.df,aes(long,lat,group=group),fill="darkgray")+
    geom_text(data=ptz@data,aes(x=LON,y=LAT,label=Site))+
    coord_fixed(1.3)+
    theme_minimal()
  map_rg
  
  #####################
  # MAP FOR PAM PAPER 
  xl<-c(-71,-70.25); yl<-c(42,42.75)
  bathy_crop<-crop(spTransform(bathy,ll),extent(c(xl,yl)))
  bathy_df<-fortify(bathy_crop)
  bathy50_df<-fortify(subset(bathy_crop,CONTOUR%in%c(-50)))
  bathy60_df<-fortify(subset(bathy_crop,CONTOUR%in%c(-60)))
  
  xt<-extent(c(-71,-69.5,42,42.75))
  shore<-crop(spTransform(land,ll),xt)
  shore@data$id = rownames(shore@data)
  shore.points = fortify(shore, region="id")
  shore.df = join(shore.points, shore@data, by="id")
  

}

##################################################
#SPATIAL PATTERN to pG, rG, & SEASON - MANUSCRIPT
if(TRUE){
  sitez<-summaryBy(cbind(LAT,LON,DEPTH)~Site,data=marus,FUN=mean,keep.names=TRUE)
  sitez_ptz<-SpatialPointsDataFrame(sitez[c("LON","LAT")],proj4string=ll,data=sitez)
  
  ptz<-merge(sitez_ptz[c("Site","LAT","LON","DEPTH")],peakz)
  ptz<-merge(ptz,emm_s)
  ptz$peak_pG<-as.Date(ptz$peak_pG,origin='1970-01-01')
  
  asprat<-1.2
  cxt<-4
  ptz<-subset(ptz,!Site%in%c("3W","7W","18"))
  ptz_nocod<-subset(ptz_nocod,Site!=18)
  ptz$peak_pG<-as.numeric(ptz$peak_pG)
  ptz$peak_rg<-as.numeric(ptz$peak_rg)
  dbrx<-as.numeric(seq.Date(as.Date("2011-10-1"),as.Date("2012-1-1"),by="month"))
  dlim<-range(dbrx,na.rm=TRUE)+c(-3,3)
  
  xl<-c(-71,-70.25); yl<-c(42,42.75)
  bathy_crop<-crop(spTransform(bathy,ll),extent(c(xl,yl)))
  bathy_crop<-subset(bathy_crop,CONTOUR%in%seq(-10,-200,by=-10))
  bathy_df<-fortify(bathy_crop)
  bathy50_df<-fortify(subset(bathy_crop,CONTOUR%in%c(-50)))
  bathy60_df<-fortify(subset(bathy_crop,CONTOUR%in%c(-60)))
  
  xt<-extent(c(-71,-69.5,42,42.75))
  shore<-crop(spTransform(land,ll),xt)
  shore@data$id = rownames(shore@data)
  shore.points = fortify(shore, region="id")
  shore.df = join(shore.points, shore@data, by="id")
  
  clx_land<-"gray70"
  clx_land_border<-"gray30"
  clx_txt<-"black"
}

if(TRUE){
  pam_pg<-ggplot()+
    geom_path(data=bathy_df,aes(x=long,y=lat,group=group),color=alpha('gray',0.5))+
    geom_path(data=bathy50_df,aes(x=long,y=lat,group=group),color='gray',size=1)+
    geom_polygon(data=shore.df,aes(long,lat,group=group),fill=clx_land,color=clx_land_border)+
    geom_point(data=ptz@data,aes(x=LON,y=LAT,size=pG,col=peak_pG))+
    labs(size="Occurrence\nProbability",color="Peak Date")+
    guides(size=guide_legend(order=1),color=guide_colorbar(order=2))+
    geom_text(data=ptz@data,aes(x=LON,y=LAT,label=Site),nudge_x=-0.03,size=cxt)+
    scale_color_gradientn(colors=c("red","blue","cyan"),labels=gimme_date,breaks=dbrx,limits=dlim)+
    scale_size_continuous(range=c(1,10),breaks=c(0.05,0.1,0.2,0.4,0.8),limits=c(0,0.8),labels=function(x){paste(x*100,"%",sep="")})+
    #xlim(xt[1],xt[2])+ylim(xt[3],xt[4])+
    geom_point(data=ptz_nocod@data,aes(x=LON,y=LAT),size=1.5,shape=4,col="black")+
    geom_text(data=ptz_nocod@data,aes(x=LON,y=LAT,label=Site),nudge_x=-0.03,size=cxt)+
    theme_minimal()+
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = c(0.85,0.6),
      plot.margin=unit(c(0,0,0,0),"cm"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust=0,color="pink",size=24)
    )+
    coord_fixed(asprat)+
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
    annotation_custom(grobTree(textGrob("a)",x=0.05,y=0.95,hjust=0,gp=gpar(col=clx_txt,fontsize=20))))+
    scalebar(bathy_df,anchor=c(x=-70.95,y=42.05), dist=5, st.dist=0.03, st.size=4, dist_unit="km",transform=TRUE, model="WGS84",location="bottomleft",box.fill="transparent",box.color=clx_txt,st.color=clx_txt)
  
  size_brx<-1*2^(0:4)
  pam_rg<-ggplot()+
    geom_path(data=bathy_df,aes(x=long,y=lat,group=group),color=alpha('gray',0.5))+
    geom_path(data=bathy50_df,aes(x=long,y=lat,group=group),color='gray',size=1)+
    geom_polygon(data=shore.df,aes(long,lat,group=group),fill=clx_land,color=clx_land_border)+
    geom_point(data=ptz@data,aes(x=LON,y=LAT,size=rG,col=peak_rg))+
    labs(size="Grunts / Hr",color="Peak Date")+
    guides(size=guide_legend(order=1),color=guide_colorbar(order=2))+
    geom_text(data=ptz@data,aes(x=LON,y=LAT,label=Site),nudge_x=-0.03,size=cxt)+
    scale_color_gradientn(colors=c("red","blue","cyan"),labels=gimme_date,breaks=dbrx,limits=dlim)+
    scale_size_continuous(range=c(1,10),breaks=size_brx,limits=c(0,max(size_brx)))+
    geom_point(data=ptz_nocod@data,aes(x=LON,y=LAT),size=1.5,shape=4,col="black")+
    geom_text(data=ptz_nocod@data,aes(x=LON,y=LAT,label=Site),nudge_x=-0.03,size=cxt)+
    theme_minimal()+
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = c(0.85,0.6),
      plot.margin=unit(c(0,0,0,0),"cm"),
      legend.background = element_rect(fill = "white")
    )+
    coord_fixed(asprat)+
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
    annotation_custom(grobTree(textGrob("b)",x=0.05,y=0.95,hjust=0,gp=gpar(col=clx_txt,fontsize=20))))
  
  #OUTPUT COMBINED FINAL MAPS
  ht<-5 #HEIGHT OF PLOT IN INCHES
  nr<-1; nc<-2; lmat<-matrix(seq_len(nr*nc),nrow=nr,ncol=nc,byrow=TRUE)
  #ggsave("pammap.jpg",marrangeGrob(list(pam_pg,pam_rg),layout_matrix=lmat,top="",padding=unit(0,"line")),width=ht*2.2,height=ht,dpi=600)
  ggsave("pammap.pdf",marrangeGrob(list(pam_pg,pam_rg),layout_matrix=lmat,top="",padding=unit(0,"line")),width=ht*2.2,height=ht,dpi=600)
  ggsave("pammap2.pdf",marrangeGrob(list(pam_pg,pam_rg),layout_matrix=lmat,top="",padding=unit(0,"line")),width=ht*2.2,height=ht,dpi=600,colormodel="cmyk")
}

emm_s<-emmeans(zbest,~Site,type='response')

####################################
#STUDY AREA FIGURE FOR MANUSCRIPT
jpeg("study_area.jpg")
cxpt<-1.2
crp<-colorRampPalette(c("black","white"))
plot(bathy_crop,col="transparent")
plot(dem_crop,add=TRUE,col=crp(100),legend=FALSE)
plot(bathy_crop,col=alpha("black",0.2),add=TRUE)
plot(shore,add=TRUE,col="black")
plot(sitez_ptz,cex=cxpt,pch=19,add=TRUE)
text(coordinates(sitez_ptz),labels=sitez_ptz$Site,pos=2)
dev.off()


#FIG 2 - RECREATE FROM PAULS VERSION IN PDF 
ag<-aggregate(GPH~Date,data=subset(gsub,Site=='1'&Y%in%c(2009)),FUN=sum)
plot(GPH~Date,data=ag,type='h',lwd=2)



#PEAK vs DEPTH
plot(-DEPTH~as.Date(peak_pG,origin='1970-01-01'),data=ptz@data,cex=bubbles(pG),col=alpha("black",0.5),pch=19)
plot(-DEPTH~as.Date(peak_rg,origin='1970-01-01'),data=ptz@data,cex=bubbles(rG),col=alpha("black",0.5),pch=19)
sub<-subset(ptz@data,!is.na(peak_pG))
cor.test(sub$DEPTH,as.numeric(sub$peak_pG))
cor.test(sub$DEPTH,as.numeric(sub$peak_rg))

plot(LAT~as.Date(peak_pG,origin='1970-01-01'),data=ptz@data,cex=bubbles(pG),col=alpha("black",0.5),pch=19)



#ANNUAL by SITE
if(TRUE){
  nclx<-8
  lwx<-2
  sitez_ptz$peak<-NA
  xl<-range(emm_js$edate)
  par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(4,5,1,6))
  plot(pG~edate,data=emm_js,xlim=xl,ylim=c(0,1.05),type='n',xaxt='n',yaxt='n',ylab="",xlab="")
  mtext(side=2,c("Min","Max"),at=c(0,1),las=2,line=1,xlab="")
  site_list<-as.character(unique(gsub$Site))
  for(s in 1:length(site_list)){
    site_this<-sitz[s]
    ltx<-ceiling(s/8)
    sub<-subset(emm_js,Site==site_this)
    sub$z<-norm01(sub$pG)
    lines(z~edate,data=sub,type='l',col=s,lty=ltx,lwd=lwx)
    peak<-mean(sub$J[sub$z==max(sub$z)])
    sitez_ptz$peak[s]<-peak
  }
  #legend("bottomright",sitez_ptz$Site,col=1:nclx,lty=c(rep(1,8),rep(2,8)),lwd=lwx)
  mtext(side=2,"Grunt Presence",line=3,cex=1.2)
  
  plot(rG~edate,data=emm_js,xlim=xl,ylim=c(0,1.05),type='n',yaxt='n',ylab="",xlab="")
  mtext(side=2,c("Min","Max"),at=c(0,1),las=2,line=1,xlab="")
  for(s in 1:length(site_list)){
    site_this<-sitz[s]
    ltx<-ceiling(s/8)
    sub<-subset(emm_js,Site==site_this)
    sub$z<-norm01(sub$rG)
    lines(z~edate,data=sub,type='l',col=s,lty=ltx,lwd=lwx)
    peak<-mean(sub$J[sub$z==max(sub$z)])
    sitez_ptz$peak[s]<-peak
  }
  #legend("bottomright",sitez_ptz$Site,col=1:nclx,lty=c(rep(1,8),rep(2,8)),lwd=lwx)
  mtext(side=2,"Grunt Rate",line=3,cex=1.2)
  par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,1,1))
  legend("right",site_list,col=1:nclx,lty=c(rep(1,8),rep(2,8)),lwd=lwx,cex=0.8,inset=0.03)
}
aggregate(rG~Site,data=emm_js,FUN=mean)
#


ag_eff<-aggregate(GPH~DEPTH,data=gruntz,FUN=length)
ag_pos<-aggregate(GPH~DEPTH,data=gruntz,FUN=sum)
ag_eff$CSUM<-normy(cumsum(ag_eff$GPH))
ag_pos$CSUM<-normy(cumsum(ag_pos$GPH))

plot(GPH~DEPTH,data=ag_eff)
plot(log(GPH)~DEPTH,data=ag_pos)


aggregate(GPH~Year,data=gruntz,FUN=sum,subset=Site==1)

cur<-read.csv("a_buoy_current_velocity.csv")
cur$date<-as.Date(cur$date,'%m/%d/%Y')
cur$MOON<-lunar.phase(cur$date,name=FALSE)
cur$MOON4<-lunar.phase(cur$date,name=4)
cur$MOON8<-lunar.phase(cur$date,name=8)

boxplot(m10~MOON8,data=cur,las=2)
boxplot(m50~MOON8,data=cur,las=2)

tab<-table(gruntz[c("Site","Year")])
rownames(tab)
ord<-c(1,12,15:22,2:11,13:14)
write.csv(as.matrix(tab[ord,]),"pam_effort.csv")

tab2<-as.data.frame.matrix(tab[ord,])
tab2$Site<-rownames(tab2)
maru_coords<-summaryBy(cbind(LAT,LON,DEPTH)~Site,data=marus,FUN=mean,keep.names=TRUE)
maru_coords[c("LAT","LON")]<-round(maru_coords[c("LAT","LON")],4)
maru_coords$DEPTH<-round(maru_coords$DEPTH)
ppos<-function(x){
  ntot<-length(x)
  npos<-length(x[x>0])
  return(npos/ntot)
}
ag_ppos<-aggregate(GPH~Site,data=gruntz,FUN=ppos)
names(ag_ppos)[ncol(ag_ppos)]<-"PPOS"
ag_ppos$PPOS<-round(100*ag_ppos$PPOS,1)

tab3<-merge(merge(tab2,maru_coords),ag_ppos)
write.csv(tab3[ord,],"pam_effort2.csv")

#########################################################
#RESPONSE TO REVIEWERS?

ag<-aggregate(GPH~H,data=gsub,FUN=mean)
barplot(ag$GPH,names.arg=ag$H)
gsub$DAYTIME<-gsub$H%in%c(16:24,0:5)
ag<-aggregate(GPH~DAYTIME,data=gsub,FUN=mean)
barplot(ag$GPH,names.arg=ag$DAYTIME)
tt<-t.test(GPH~DAYTIME,data=gsub)
tt$estimate[2]/tt$estimate[1]

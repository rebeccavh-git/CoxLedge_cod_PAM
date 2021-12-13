####################################
# STANDARD SET OF USEFUL FUNCTIONS
# MICAH DEAN 
####################################
options(scipen=10)
mlabz<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")

####################
#REQUIRED LIBRARIES
library(rgdal)
library(rgeos)
library(maptools)
library(MuMIn)

###############################
#CALCULATION & DATA MANAGEMENT
if(TRUE){

  roundy<-function(x,digits=0){
    posneg = sign(x)
    z = abs(x)*10^digits
    z = z + 0.5
    z = trunc(z)
    z = z/10^digits
    return(z*posneg)
  }
  
  #NORMALIZE A VECTOR TO EITHER SUM TO ONE OR MAX OUT AT ONE
  normy<-function(x,type='0to1'){
    if(type=='0to1'){
      x<-x-min(x,na.rm=TRUE)
      return(x/max(x,na.rm=TRUE))
    }
    if(type=='max')return(x/max(x,na.rm=TRUE))
    if(type=='sum')return(x/sum(x,na.rm=TRUE))
  }
  
  rescale<-function(x,minval=0,maxval=1){
    xmin<-min(x,na.rm=TRUE)
    xmax<-max(x-xmin,na.rm=TRUE)
    newx<-((x-xmin)/xmax*(maxval-minval))+minval
    return(newx)
  }
  
  #LOGIT TRANSFORMATION
  logit<-function(x){log(x/(1-x))}
  
  #INVERSE LOGIT
  invlogit<-function(x){exp(x)/(1+exp(x))}
  
  #CALCULATE PERCENT ZEROS of a VECTOR OF VALUES
  pz<-function(x){
    return(length(x[x==0])/length(x))
  }
  
  #DROP COLUMNS FROM A DATAFRAME
  drop_colz<-function(df,colnames){
    return(df[names(df)[!names(df)%in%colnames]])
  }
  
  #LIST ALL THE OBJECTS in MEMORY, SORTED BY SIZE (in MB)
  mem_hogz<-function(obz,top=10){
    df<-as.data.frame(sort(sapply(obz,function(x){object.size(get(x))/1e6}),decreasing=TRUE))
    names(df)<-"MB"
    return(head(df,top))
  }
  
  #Z-TRANSFORMATION
  ztran<-function(x){
    return((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  }
  
  #EXPAND A "TALLY" DATAFRAME INTO AN "INDIVIDUAL" DATAFRAME
  expand_tally<-function(x,tallycol){
    colz<-names(x)
    out<-data.frame()
    for(this_col in colz){
      expanded<-data.frame(rep(x[[this_col]],x[[tallycol]]))
      names(expanded)<-this_col
      if(this_col==colz[1])out<-expanded
      if(this_col!=colz[1])out<-cbind(out,expanded)
    }
    out[[paste(tallycol,"_OG",sep="")]]<-out[[tallycol]]
    out[[tallycol]]<-1
    return(out)
  }
  
}

#############
#BIOLOGICAL

# VON BERT GROWTH FUNCTION
vb<-function(length=NULL,age=NULL,linf=NULL,k=NULL,t0=NULL){
  length<-ifelse(length>linf,linf,length)
  if(!is.null(length))out<-t0+(log(1-(length/linf))/-k)
  if(!is.null(age))out<-linf*(1-exp(-k*(age-t0)))
  return(out)
}

# LENGTH-WEIGHT FUNCTIONS
lw<-function(length=NULL,a=NULL,b=NULL){
  w<-a*length^b
  return(w)
}

###########################
# MODEL PREDICTION, EVAL, & DIAGNOSTICS 
if(TRUE){
  
  #OVERDISPERSION FOR A MODEL - WORKS FOR GAM, GLM, OTHERS?
  disp = function(mod = NA){
    E = residuals(mod, type='pearson')
    d = sum(E^2)/mod$df.resid
    return(d)
  }
  
  #CALC AIC or BIC FOR A LIST OF MODELS
  LIC<-function(x,crit="AIC"){
    if(crit=="AIC"){IC<-function(x,...)AIC(x,...)}
    if(crit=="BIC"){IC<-function(x,...)BIC(x,...)}
    out<-NULL
    for(i in 1:length(x)){
      ic<-try(IC(x[[i]]))
      this<-NA
      if(!'try-error'%in%class(ic))this<-ic
      out<-c(out,this)
    }
    return(out)
  }
  
  #CALC EQUIVALENT DF FOR A GAM, GLM, etc
  edf<-function(x){
    isbad<-'try-error'%in%class(try(nobs(x),silent=TRUE))
    if(!isbad)return(nobs(x)-df.residual(x))
    if(isbad)return(NA)
  }
  
  #CALC n FOR A GAM, GLM, etc
  nobz<-function(x){
    return(nobs(x))
  }
  
  #DEVIANCE EXPLAINED for a GAM
  devexp<-function(x){
    out<-NA
    if("gam"%in%class(x))out<-summary(x)$dev.expl
    return(out)
  }
  
  #CREATE AN AIC/BIC/df TABLE FOR LIST OF MODEL OBJECTS
  ictab<-function(x){
    modz<-names(x)
    aicz<-LIC(x,crit="AIC")
    bicz<-LIC(x,crit="BIC")
    dfz<-unlist(lapply(x,FUN=edf))
    dfz<-round(unlist(lapply(x,FUN=edf)),3)
    daicz<-round(aicz-min(aicz,na.rm=TRUE),1)
    dbicz<-round(bicz-min(bicz,na.rm=TRUE),1)
    wa<-aicz; wa[is.na(wa)]<-Inf;   wa<-round(Weights(wa),3)
    wb<-bicz; wb[is.na(wb)]<-Inf;   wb<-round(Weights(wb),3)
    #de<-round(unlist(lapply(x,FUN=devexp)),1)
    df<-data.frame(mod=modz,wA=wa,wB=wb,dAIC=daicz,dBIC=dbicz,AIC=aicz,BIC=bicz,edf=dfz)
    rownames(df)<-NULL
    return(df)
  }
  
  #FUNCTION TO MAKE PREDICTIONS USING MODEL AVERAGING
  predict_multimod<-function(modz,weights=rep(1,length(modz)),newdata=NULL,type='response',se.fit=FALSE){
    out<-NULL
    for(m in modz){
      this<-predict(m,newdata=newdata,type=type,se.fit=se.fit)
      out<-cbind(out,this)
    }
    weights<-weights/sum(weights) #NORMALIZE TO SUM TO 1
    return(rowSums(t(t(out)*weights)))
  }
  
  
}

####################
#GENERAL PLOTTING 
if(TRUE){
  
  #A "MOUNTAIN" - filled polygon with bottom at zero and top a Y~X
  mountain<-function(form,data,col="gray",border="transparent",add=FALSE,...){
    xcol<-as.character(form[3])
    ycol<-as.character(form[2])
    dat<-data[order(data[xcol]),]
    x<-dat[[xcol]]
    y<-dat[[ycol]]
    if(!add)plot(x,y,col='transparent',xlab=xcol,ylab=ycol,...)
    polygon(c(x,rev(x)),c(rep(0,length(y)),rev(y)),col=col,border=border)
    #box()
  }
  
  #BUBBLE CEX VALUE FOR BUBBLE PLOTS
  bubbles<-function(x,cxmin=1,cxmax=5,type='radius',maxval=max(x,na.rm=TRUE)){
    if(type=='area'){x<-x^2; maxval<-maxval^2}
    minval<-min(x,na.rm=TRUE)
    maxval<-maxval-minval
    xx<-normy(x,type='0to1')
    #out<-cxmin+((xx/maxval)*(cxmax-cxmin))
    out<-cxmin+(xx*(cxmax-cxmin))
    if(length(xx)==1)out<-mean(c(cxmin,cxmax))
    return(out)
  }
    
  #CREATE MULTIPANEL PLOT
  par_panel<-function(nrow,ncol,byrow=TRUE,oma=c(5,4,4,1),mar=c(0,0,0,0)){
    if(byrow)par(mfrow=c(nrow,ncol),oma=oma,mar=mar)
    if(!byrow)par(mfcol=c(nrow,ncol),oma=oma,mar=mar)
  }
  
  #RETURN TO PAR DEVAULTS
  par_unpanel<-function(){
    par(mfrow=c(1,1),mar=c(5,4,4,1),oma=c(0,0,0,0))
  }
  
  #PLOTS A RASTER OBJECT, BUT CAPS VALUES AT A MAX PERCENTILE 
  zplot<-function(r,p=0.99,samescale=TRUE,...){
    q<-quantile(raster::values(r),p,na.rm=TRUE)
    r[r>q]<-q
    zl<-range(raster::values(r),na.rm=TRUE)
    if(!samescale)plot(r,...)
    if(samescale)plot(r,zlim=zl,...)
  }
  
  #SAME AS ZPLOT, ONLY USING image() WHICH WORKS BETTER FOR TILED PLOTS
  zimage<-function(r,p=0.99,samescale=TRUE,...){
    q<-quantile(values(r),p,na.rm=TRUE)
    r[r>q]<-q
    zl<-range(values(r),na.rm=TRUE)
    if(!samescale)image(r,...)
    if(samescale)image(r,zlim=zl,...)
  }
  
  #MEASURE EUCLIDAN DISTANCE BETWEEN POINTS IN MAP PROJECTION UNITS 
  measure<-function(nclick=2){
    xy<-locator(nclick)
    return(sqrt(diff(xy$x)^2+diff(xy$y)^2))
  }
  
  #CAPITALIZE 1st LETTER OF EACH
  initcap <- function(x,all=TRUE) {
    s<-x
    if(all)s <- strsplit(as.character(x), " ")[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "", collapse = " ")
  }
  
}

##############
#GIS FUNCTIONS
if(TRUE){
  
  #GIS PROJECTIONS 
  ll<-CRS("+init=epsg:4326") #NAD83/WGS84  (i.e., LAT & LON)
  masp<-CRS("+init=epsg:26986")  #MASS STATE PLANE - METER
  
  #CONVERT NEFSC LAT/LON FORMAT TO DECIMAL DEGREES 
  getll<-function(x){
    n<-nchar(as.character(x)) #LENGTH OF STRING
    if(substr(x,1,1)==0){
      x<-substr(x,2,n) #DROP LEADING ZERO FOR LONGITUDE, IF PRESENT
      n<-n-1
    }
    islat<-substr(x[1],1,1)<6 #ASSUME THIS IS LATITUDE IF FIRST DIGIT IS < 6
    islon<-substr(x[1],1,1)>=6 #ASSUME THIS IS LONGITUDE IF FIRST DIGIT IS >= 6
    x<-substr(x,1,ifelse(substr(x,n,n)%in%c("N","W"),n-1,n)) #DROP TRAILING LETTER, IF PRESENT
    d<-as.numeric(substr(x,1,2)) #DEGREES
    m<-as.numeric(substr(x,3,8)) #DECIMAL MINUTES
    if(islat)return(d+m/60)
    if(islon)return(-d-m/60)
  }
  
  #CREATE SPATIAL POINTS DATAFRAME IN MASP FROM DATAFRAME, USING LAT/LON
  make_ptz<-ptzfun<-function(x,latcol="LAT",loncol="LON"){
    col_idx<-which(toupper(names(x))%in%toupper(c(latcol,loncol)))
    names(x)[col_idx]<-toupper(names(x)[col_idx])
    spTransform(SpatialPointsDataFrame(x[toupper(c(loncol,latcol))],proj4string=ll,data=x),masp)
  }
  
  #CREATE A SPATIAL LINES OBJECT FROM A DATAFRAME (that has LAT1, LAT2, LON1, LON2)
  make_linez<-function(dat){
    col_idx<-which(toupper(names(dat))%in%c("LAT1","LAT2","LON1","LON2"))
    names(dat)[col_idx]<-toupper(names(dat)[col_idx])
    dat<-subset(dat,!is.na(LAT1)&!is.na(LON1)&!is.na(LAT2)&!is.na(LON2))
    l<-vector("list", nrow(dat))
    for(i in 1:length(l)){
      l[[i]]<-Lines(list(Line(rbind(c(dat$LON1[i], dat$LAT1[i]),c(dat$LON2[i], dat$LAT2[i])))), as.character(i))
    }
    sl<-SpatialLines(l,proj4string=ll)
    linez<-spTransform(SpatialLinesDataFrame(sl,dat,match.ID=FALSE),masp)
    linez$KM<-gLength(linez,byid=TRUE)
    return(linez)
  }
  
  #CREATE SPATIAL LINES OBJECT (FROM A LAT/LON SEQUENCE)
  make_line<-function(dat){
    col_idx<-which(toupper(names(dat))%in%c("LAT","LON"))
    names(dat)[col_idx]<-toupper(names(dat)[col_idx])
    dat<-subset(dat,!is.na(LAT)&!is.na(LON))
    l<-list()
    l[[1]]<-Lines(list(Line(dat[c("LON","LAT")])), as.character(i))
    sl<-spTransform(SpatialLines(l,proj4string=ll),masp)
    sl$LENGTH_M<-gLength(sl,byid=TRUE)
    return(sl)
  }
  
  
  #CUSTOM SCALEBAR FUNCTION
  myscalebar<-function(x=0.05,y=0.06,len=50,barunits="km",prj=masp,ntics=4,ticdir="up",ticlen=0.01,cex=0.7,col="black"){
    masp<-CRS("+init=epsg:26986") #MA STATE PLANE METERS
    bl<-par("usr")[c(1,3)] #COORDS FOR BOTTOM LEFT OF PLOT
    tr<-par("usr")[c(2,4)] #COORDS FOR TOP RIGHT OF PLOT
    wh<-tr-bl #WIDTH AND HEIGHT OF PLOT
    conv<-switch(barunits,m=1,km=1000,nm=1852,mi=1609.34) #ASSUMES MAPUNITS = m
    a<-bl+c(wh[1]*x,wh[2]*y) #LEFT END OF BAR (in native projection)
    a.m<-spTransform(SpatialPoints(t(a),proj4string=prj),masp) #LEFT END in MASP
    b.m<-SpatialPoints(coordinates(spTransform(a.m,masp))+t(c(len*conv,0)),proj4string=masp) #RIGHT END in MASP
    b<-coordinates(spTransform(b.m,prj)) #RIGHT END OF BAR (in native projection)
    b[2]<-a[2] #PREVENTS CROOKED BARS
    segments(a[1],a[2],b[1],b[2],col=col) #DRAW THE BAR
    tix<-cbind(seq(from=a[1],to=b[1],length=ntics),rep(a[2],ntics)) #TIC LOCATIONS
    dir<-ifelse(ticdir=="up",1,-1)
    segments(tix[,1],tix[,2],tix[,1],tix[,2]+wh[2]*ticlen*dir,col=col) #DRAW UP TICS
    text(tix[c(1,ntics),1],tix[c(1,ntics),2]+wh[2]*3*ticlen*dir #LABEL FIRST AND LAST TICS
         ,c("0",paste(len,barunits)),adj=0.2,cex=cex,col=col)
  }
  
  #USE LAT/LON TO MAKE EXTENT IN A DIFFERENT CRS
  make_xt<-function(x,crs=ll){
    return(extent(spTransform(ptzfun(data.frame(LON=x[1:2],LAT=x[3:4])),crs)))
  }
  
  zoom_xt<-function(x,f=3){
    xt<-extent(spTransform(x,ll))
    xd<-abs(diff(xt[1:2]))
    yd<-abs(diff(xt[3:4]))
    out<-ptzfun(data.frame(LON=xt[1:2]+(xd*f*c(-1,1)),LAT=xt[3:4]+(yd*f*c(-1,1))))
    return(spTransform(out,crs(x)))
  }
  
  refgrid<-function(xt,res=1000){
    xt<-extent(spTransform(xt,masp))
    xseq<-seq(xt[1],xt[2],by=res)
    yseq<-seq(xt[3],xt[4],by=res)
    xy <- expand.grid(x=xseq,y=yseq)
    xy<-SpatialPoints(xy,masp)
    return(xy)
  }
  
  #MAKE A GRID OF POLYGONS
  poly_grid<-function(x1,x2,y1,y2,width,height,crx=ll){
    xs<-seq(x1,x2-width,width)+width/2
    ys<-seq(y1,y2-height,height)+height/2
    xy<-expand.grid(x=xs,y=ys)
    xy<-xy[order(-xy$y,xy$x),] #ORDER SO NUMBERS ARE LIKE READING A PAGE
    pts_grid<-SpatialPointsDataFrame(xy,data=xy,proj4string=crx)
    gridded(pts_grid)<-TRUE
    grid_out<-as(pts_grid, "SpatialPolygons")
    grid_out$ID<-1:length(grid_out)
    return(grid_out)
  }
  
  #DISSOLVE A spatialPolygonsDataFrame OBJECT BY A COLUMN, KEEPING ATTRIBUTES (alt to gUnaryUnion)
  dissolve <- function(SPDF,attribute){
    rownames(SPDF@data) <- sapply(slot(SPDF, "polygons"), function(x) slot(x, "ID"))   
    Temp <- gUnaryUnion(SPDF,attribute)
    IDlist <- data.frame(ID=sapply(slot(Temp, "polygons"), function(x) slot(x, "ID")))
    rownames(IDlist)  <- IDlist$ID
    SpatialPolygonsDataFrame(Temp,IDlist)
  }
  
  #FUNCTION TO CREATE DISSOLVED GROUNDFISH CLOSURE POLYGON FOR ANY MONTH/YEAR
  clhist<-read.csv("C:\\Users\\mdean\\Google Drive\\Cod\\closure_history.csv") #TABLE OF POLYGON NAMES AND WHETHER CLOSED BY YEAR, MONTH
  if(TRUE){
    crs_use<-masp
    gispath<-"C:\\Users\\mdean\\Google Drive\\GIS\\data"
    sq30<-spTransform(readOGR(gispath,"thirty_min_squares"),crs_use)
    sq30$sq<-sq30$block
    wgom<-spTransform(readOGR(gispath,"wgom"),crs_use)
    wgom2<-spTransform(readOGR(gispath,"wgom2"),crs_use)
    cashes<-spTransform(readOGR(gispath,"cashes"),crs_use)
    wccz<-spTransform(readOGR(gispath,"wccz"),crs_use)
    sccz<-spTransform(readOGR(gispath,"sccz_current"),crs_use)
    wb<-spTransform(readOGR(gispath,"gom_cspa"),crs_use)
    top125<-spTransform(readOGR(gispath,"top125"),crs_use)
    top124<-spTransform(readOGR(gispath,"top124"),crs_use)
    bump124<-spTransform(readOGR(gispath,"bump124"),crs_use)
    bump124<-spTransform(buffer(spTransform(bump124,masp),0),crs_use)
  }
  closure<-function(month,year=format(Sys.time(),'%Y'),crs=masp){
    clyr<-subset(clhist,Y1<=year&ifelse(is.na(Y2),year,Y2)>=year)
    xym<-data.frame(x=rep(-70,5),y=rep(42,5)) #A 'DUMMY' POLYGON AT A SINGLE TO CREATE SP OBJECT
    cl<-spTransform(SpatialPolygons(list(Polygons(list(Polygon(xym)),1)),proj4string=ll),crs)
    #30-MIN SQUARES
    blox<-clyr$NAME[clyr$TYPE=='sq30'&clyr[paste("M",month,sep="")]]
    if(length(blox)>0){
      sub<-subset(sq30,block%in%blox)
      cl<-spTransform(unionSpatialPolygons(sub,rep(1,length(sub))),crs)
    }
    #OTHER CLOSURES
    shps<-as.character(clyr$NAME[clyr$TYPE=='shp'&clyr[paste("M",month,sep="")]])
    if(length(shps)>0){
      for(i in 1:length(shps)){
        #cl<-gUnion(cl,spTransform(get(shps[i]),crs),drop_lower_td=TRUE)
        cl<-union(cl,spTransform(get(shps[i]),crs))
      }
    }
    #CROP TO IBS2 STUDY AREA
    cl<-gBuffer(cl,width=0)
    #cl<-crop(cl,spTransform(ibs2area,crs))
    return(cl)
  }
  
  #ADD GRATICULE TICS TO A MAP
  maptix<-function(xtix,ytix,xlabs=xtix,ylabs=ytix,crx=ll,sides=c(1,2),tcx=c(-0.25,-0.5),cx=1,mgx=0.5){
    xdir<-"W"
    ydir<-"N"
    xax_lab<-paste(-xlabs,'°',sep="")
    xax_lab[length(xlabs)]<-paste(xax_lab[length(xlabs)],xdir)
    yax_lab<-paste(ylabs,'°',sep="")
    yax_lab[length(ylabs)]<-paste(yax_lab[length(ylabs)],ydir)
    xlab_loc<-coordinates(spTransform(SpatialPoints(data.frame(xlabs,rep(mean(ylabs),length(xlabs))),proj4string=ll),crx))[,1]
    ylab_loc<-coordinates(spTransform(SpatialPoints(data.frame(rep(mean(xlabs),length(ylabs)),ylabs),proj4string=ll),crx))[,2]
    xtix_loc<-coordinates(spTransform(SpatialPoints(data.frame(xtix,rep(mean(ytix),length(xtix))),proj4string=ll),crx))[,1]
    ytix_loc<-coordinates(spTransform(SpatialPoints(data.frame(rep(mean(xtix),length(ytix)),ytix),proj4string=ll),crx))[,2]
    for(i in sides){
      if(i %in% c(1,3)){tix_loc<-xtix_loc; labs<-xax_lab; lablocs<-xlab_loc; lasx<-1}
      if(i %in% c(2,4)){tix_loc<-ytix_loc; labs<-yax_lab; lablocs<-ylab_loc; lasx<-2}
      axis(i,at=tix_loc,labels=FALSE,tcl=tcx[1],mgp=c(3,mgx,0))
      axis(i,at=lablocs,labels=labs,cex.axis=cx,las=lasx,tcl=tcx[2],mgp=c(3,mgx,0))
    }
  }

  # produces text with a shadowed border  
  shadowtext = function(x, y = NULL, labels, col = "white", bg = "grey25", theta = (1:8/4)*pi, r1 = 0.06, r2 = 0.04, ...){
    xy = xy.coords(x,y)
    xo = r1*strwidth("N")
    yo = r2*strheight("N")
    for (i in theta){
      text(xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ...)
    }
    text(xy$x, xy$y, labels, col=col, ...)
    text(xy$x, xy$y, labels, col=col, ...)
  }
  
  #SIMPLE NORTH ARROW USING LOCATOR
  north<-function(col="black"){
    xy<-locator(1)
    points(xy,pch=17,col=col)
    text(xy,"N",pos=1,col=col)
  }
  
  
}
  
##################################
# SELECTIVITY ESTIMATION 
if(TRUE){
  
  #NLL FUNCTIONS
  if(TRUE){
    
    #NLL FUN for linear
    nll_lin<-function(len,n_exp,n_tot,a,b){
      p<-a+b*len
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for FIXED 
    nll_fixed<-function(len,n_exp,n_tot,a){
      p<-1/(1+exp(-a))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for 2-PAR LOGISITIC
    nll2<-function(len,n_exp,n_tot,a,b){
      p<-1/(1+exp(-a-b*len))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for 3-PAR LOGISITIC
    nll3<-function(len,n_exp,n_tot,a,b,c){
      p<-inv.logit(c)/(1+exp(-a-b*len))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for 3-PAR LOGISITIC
    #nll3<-function(len,n_exp,n_tot,a,b,c){
      #p<-c/(1+exp(-a-b*len))
      #y<-n_exp/n_tot
      #nll<-log((p^y)*(1-p)^(1-y))*-1
      #nll<-log(dbinom(y,rep(1,length(y)),p))*-1
      #nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      #sum(nll)
    #}
    
    
    #NLL FUN FOR RICHARDS CURVE (4p)
    nllrich<-function(len,n_exp,n_tot,a,b,c,d){
      p<-inv.logit(c)*((1/(1+exp(-a-b*len)))^(1/d))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN FOR GAMMA DOME 
    nllgamma<-function(len,n_exp,n_tot,a,b,c){
      h<-0.5*(sqrt(a^2+4*b^2)-a)
      p<-c*((len/a)^(a/h))*exp((len-a)/h)
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for 4-PAR LOGISITIC -NOT WORKING!
    nll4<-function(len,n_exp,n_tot,a,b,c,d){
      p<-c*(1/(1-d))*(((1-d)/d)^d)*(exp(-d*(a+b*len))/(1+exp(-a-b*len)))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #NLL FUN for 5-PAR DOUBLE LOGISITIC
    nll5<-function(len,n_exp,n_tot,a,b,c,d,e){
      p<-inv.logit(e)/((1+exp(-a*(len-b)))*(1+exp(c*(len-d))))
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #2PAR LOG-LOGISTIC
    nll_ll2<-function(len,n_exp,n_tot,a,b){
      p<-1/(1+(len/a)^-b)
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
    
    #3PAR LOG-LOGISTIC
    nll_ll3<-function(len,n_exp,n_tot,a,b,c){
      p<-c/(1+(len/a)^-b)
      nll<-(-dbinom(n_exp,prob=p,size=n_tot,log=TRUE))
      sum(nll)
    }
  }

  #MLE FUNCTION - EXPECTS DAT be a list with these elements "len, n_exp, n_tot"
  mymle<-function(dat,newlen=dat$len,type,start=NA){
    dat$eff<-dat$sel
    #FIXED EFFICIENCY ACROSS ALL LENGTHS
    if(type=="fixed"){
      fit<-mle2(nll_fixed,start=list(a=0.5),data=dat,method='L-BFGS-B',lower=list(a=-10),upper=list(a=10))
      a<-coef(fit)[1]
      pred<-rep(1/(1+exp(-a)),length(newlen))
    }
    #2PAR LOGISTIC
    if(type=="l2"){
      fit<-mle2(nll2,start=list(a=-5,b=0.2),data=dat,method='BFGS')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      pred<-1/(1+exp(-a-b*newlen))
    }
    #3PAR LOGISTIC
    if(type=="l3"){
      #fit<-mle2(nll3,data=dat,start=list(a=-10,b=0.1,c=0),method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01),upper=list(a=Inf,b=Inf,c=0.99))
      fit<-mle2(nll3,data=dat,start=list(a=-5,b=0.2,c=3),method='BFGS')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      pred<-inv.logit(c)*(1/(1+exp(-a-b*newlen)))
    }
    #4PAR RICHARDS
    if(type=="rich"){
      fit<-mle2(nllrich,data=dat,start=list(a=-10,b=0.5,c=1,d=1),method='L-BFGS-B',lower=list(a=-Inf,b=-1e10,c=0.01,d=0.01),upper=list(a=Inf,b=1e10,c=0.99,d=Inf))
      #fit<-mle2(nllrich,data=dat,start=list(a=-10,b=0.5,c=1,d=1),method='Nelder-Mead')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      pred<-inv.logit(c)*((1/(1+exp(-a-b*newlen)))^(1/d))
    }
    
    #GAMMA DOME
    if(type=="gamma"){
      fit<-mle2(nllgamma,data=dat,start=list(a=-10,b=0.5),method='L-BFGS-B',lower=list(a=-Inf,b=-1e10,c=0.01,d=0.01),upper=list(a=Inf,b=1e10,c=0.99,d=Inf))
      #fit<-mle2(nllrich,data=dat,start=list(a=-10,b=0.5,c=1,d=1),method='Nelder-Mead')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      pred<-c*((1/(1+exp(-a-b*newlen)))^(1/d))
    }
    #4PAR LOGISTIC
    if(type=="l4"){
      fit<-mle2(nll4,start=list(a=1,b=0.5,c=max(dat$sel),d=0.5),data=dat,method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01,d=0.01),upper=list(a=Inf,b=Inf,c=1,d=0.99999))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      pred<-c*(1/(1-d))*(((1-d)/d)^d)*(exp(-d*(a+b*newlen))/(1+exp(-a-b*newlen)))
    }
    #5PAR DOUBLE LOGISTIC
    if(type=="l5"){
      #fit<-mle2(nll5,start=list(a=0.1,b=mean(dat$len)*0.66,c=0.1,d=max(dat$len),e=0.8),data=dat,method='L-BFGS-B',lower=list(a=0,b=0,c=0,d=0,e=0.1),upper=list(a=2,b=max(dat$len),c=2,d=max(dat$len)*3,e=0.99))
      fit<-mle2(nll5,start=list(a=0.1,b=60,c=0.1,d=120,e=0.8),data=dat,method='BFGS',control=list(maxit=1000))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      e<-coef(fit)[5]
      pred<-inv.logit(e)/((1+exp(-a*(newlen-b)))*(1+exp(c*(newlen-d))))
    }
    #2PAR LOG-LOGISITIC
    if(type=="ll2"){
      fit<-mle2(nll_ll2,start=list(a=25,b=0),data=dat,method='BFGS')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      pred<-1/(1+(newlen/a)^-b)
    }
    #3PAR LOG-LOGISTIC
    if(type=="ll3"){
      fit<-mle2(nll_ll3,start=list(a=25,b=0,c=mean(dat$sel)),data=dat,method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01),upper=list(a=Inf,b=Inf,c=0.99))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      pred<-c/(1+(newlen/a)^-b)
    }
    nll<-dbinom(dat$n_exp,prob=pred,size=dat$n_tot,log=TRUE)
    return(list(mod=type,fit=fit,pred=data.frame(len=newlen,p=pred,nll)))
  }
  
  #MLE USING 'formula' METHOD - ALLOWS FOR PREDICTIONS & CONF INTERVALS
  mymle2<-function(dat,newlen=dat$len,type,start=NA){
    dat$eff<-dat$sel
    #FIXED EFFICIENCY ACROSS ALL LENGTHS
    if(type=="fixed"){
      fit<-mle2(nll_fixed,start=list(a=0.5),data=dat,method='L-BFGS-B',lower=list(a=-10),upper=list(a=10))
      a<-coef(fit)[1]
      pred<-rep(1/(1+exp(-a)),length(newlen))
    }
    #2PAR LOGISTIC
    if(type=="l2"){
      fit<-mle2(n_exp~dbinom(size=n_tot,prob=1/(1+exp(-b*(len-a)))),data=dat
                ,start=list(a=50,b=0.25),method='Nelder-Mead',control=list(maxit=1000))
      
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      pred<-1/(1+exp(-b*(newlen-a)))
    }
    #3PAR LOGISTIC
    if(type=="l3"){
      #fit<-mle2(nll3,data=dat,start=list(a=-10,b=0.1,c=0),method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01),upper=list(a=Inf,b=Inf,c=0.99))
      fit<-mle2(n_exp~dbinom(size=n_tot,prob=inv.logit(c)/(1+exp(-b*(len-a)))),data=dat
                             ,start=list(a=50,b=0.25,c=2),method='Nelder-Mead',control=list(maxit=1000))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      pred<-inv.logit(c)/(1+exp(-b*(newlen-a)))
    }
    #4PAR RICHARDS
    if(type=="rich"){
      #fit<-mle2(n_exp~dbinom(size=n_tot,prob=c*((1/(1+exp(-a-b*len)))^(1/d))),data=dat,start=list(a=-10,b=0.5,c=1,d=1),method='L-BFGS-B',lower=list(a=-Inf,b=-1e10,c=0.01,d=0.01),upper=list(a=Inf,b=1e10,c=0.99,d=Inf))
      fit<-mle2(n_exp~dbinom(size=n_tot,prob=inv.logit(c)*((1/(1+exp(-b*(len-a)))^(1/d)))),data=dat
                             ,start=list(a=50,b=0.25,c=2,d=1),method='Nelder-Mead',control=list(maxit=1000))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      pred<-inv.logit(c)*((1/(1+exp(-b*(newlen-a)))^(1/d)))
    }
    
    #5PAR DOUBLE LOGISTIC
    if(type=="l5"){

      fit<-mle2(n_exp~dbinom(size=n_tot,prob=inv.logit(c)/((1+exp(-b*(len-a)))*(1+exp(-e*(len-d)))))
                ,start=list(a=50,b=0.25,c=1,d=80,e=-0.1),data=dat,method='Nelder-Mead',control=list(maxit=1000))
      
      #fit<-mle2(n_exp~dbinom(size=n_tot,prob=inv.logit(c)/((1+exp(-b*(len-a)))*(1+exp(-e*(len-d)))))
      #          ,start=list(a=50,b=0.25,c=1,d=80,e=-0.1)
      #          ,lower=list(a=0,b=0.01,c=0.1,d=50,e=-1)
      #          ,upper=list(a=80,b=1,c=1000,d=250,e=-0.01)
       #         ,data=dat,method='L-BFGS-B',control=list(maxit=1000))
      
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      e<-coef(fit)[5]
      pred<-inv.logit(c)/((1+exp(-b*(newlen-a)))*(1+exp(-e*(newlen-d))))
    }
    
    
    
    #GAMMA DOME
    if(type=="gamma"){
      fit<-mle2(n_exp~dbinom(size=n_tot,prob=((len/a)^(a/(0.5*(sqrt((a^2)+4*(b^2))-a))))*exp((a-len)/(0.5*(sqrt((a^2)+4*(b^2))-a)))),data=dat
                ,start=list(a=80,b=25),method='Nelder-Mead')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      pred<-((len/a)^(a/(0.5*(sqrt((a^2)+4*(b^2))-a))))*exp((a-newlen)/(0.5*(sqrt((a^2)+4*(b^2))-a)))
    }
    #4PAR LOGISTIC
    if(type=="l4"){
      fit<-mle2(nll4,start=list(a=1,b=0.5,c=max(dat$sel),d=0.5),data=dat,method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01,d=0.01),upper=list(a=Inf,b=Inf,c=1,d=0.99999))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      d<-coef(fit)[4]
      pred<-c*(1/(1-d))*(((1-d)/d)^d)*(exp(-d*(a+b*newlen))/(1+exp(-a-b*newlen)))
    }
    #2PAR LOG-LOGISITIC
    if(type=="ll2"){
      fit<-mle2(nll_ll2,start=list(a=25,b=0),data=dat,method='BFGS')
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      pred<-1/(1+(newlen/a)^-b)
    }
    #3PAR LOG-LOGISTIC
    if(type=="ll3"){
      fit<-mle2(nll_ll3,start=list(a=25,b=0,c=mean(dat$sel)),data=dat,method='L-BFGS-B',lower=list(a=-Inf,b=-Inf,c=0.01),upper=list(a=Inf,b=Inf,c=0.99))
      a<-coef(fit)[1]
      b<-coef(fit)[2]
      c<-coef(fit)[3]
      pred<-c/(1+(newlen/a)^-b)
    }
    nll<-dbinom(dat$n_exp,prob=pred,size=dat$n_tot,log=TRUE)
    return(list(mod=type,fit=fit,pred=data.frame(len=newlen,p=pred,nll)))
  }


  #PREDICT SELECTIVITY 
  predict_sel<-function(mod,newlen){
    type<-mod$mod
    coefs<-coef(mod$fit)
    for(i in 1:length(coefs)){assign(x=names(coefs[i]),value=coefs[i])}
    if(mod$mod=="fixed"){pred<-rep(1/(1+exp(-a)),length(newlen))}
    if(mod$mod=="l2"){pred<-1/(1+exp(-a-b*newlen))}
    if(mod$mod=="l3"){pred<-inv.logit(c)/(1+exp(-b*(newlen-a)))}
    if(mod$mod=="rich"){pred<-inv.logit(c)*((1/(1+exp(-b*(newlen-a)))^(1/d)))}
    if(mod$mod=="l4"){pred<-c*(1/(1-d))*(((1-d)/d)^d)*(exp(-d*(a+b*newlen))/(1+exp(-a-b*newlen)))}
    if(mod$mod=="l5"){pred<-inv.logit(c)/((1+exp(-b*(newlen-a)))*(1+exp(-e*(newlen-d))))}
    if(mod$mod=="ll2"){pred<-1/(1+(newlen/a)^-b)}
    if(mod$mod=="ll3"){pred<-c/(1+(newlen/a)^-b)}
    return(pred)
  }

}




  
  
  
  

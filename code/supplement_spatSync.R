### Analysis of distance decay/spatial structure in synchrony

rm(list = ls())

library(RColorBrewer)
library(ncf)
library(wsyn)
library(rgdal)

setwd(dirname(dirname(rstudioapi:::getSourceEditorContext()$path)))
##### Preliminary Data Processing/Management #####
ADChem_raw<-read.csv("./data/altm_data.csv",
                     stringsAsFactors = F)

ADChem_Dat<-ADChem_raw
ADChem_Dat<-ADChem_Dat[order(ADChem_Dat$SiteID, ADChem_Dat$ConfDate),]
ADChem_Dat$month<-as.character(substr(ADChem_Dat$ConfDate,6,7))
ADChem_Dat$month_abs<-as.character(substr(ADChem_Dat$ConfDate,1,7))

## Do some easy subsetting first to reduce number of computations needed -------
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$ConfDate),] #remove any rows without a date
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$LabID),] #remove any rows without a Lab ID
ADChem_BLANK<-ADChem_Dat[ADChem_Dat$SiteID=="BLANK",] #save the blanks separately
ADChem_Dat<-ADChem_Dat[ADChem_Dat$SiteID!="BLANK",]
ADChem_Dat<-ADChem_Dat[ADChem_Dat$Station==1,] #Remove unwanted stations
ADChem_Dat<-ADChem_Dat[!is.na(ADChem_Dat$depth),]
ADChem_Dat<-ADChem_Dat[ADChem_Dat$depth<=0.5,] #Remove unwanted depths

any(duplicated(ADChem_BLANK$ConfDate))

ADSites<-sort(unique(ADChem_Dat$SiteID))
ADMonths<-sort(unique(ADChem_Dat$month))
ADMonths_abs<-sort(unique(ADChem_Dat$month_abs))
ADCats<-rownames(ADChem_Dat)

## Do blank correction --------------------------------------------------------
ADChem_Cor<-ADChem_Dat
ADChem_Cor$cor<-rep(0, nrow(ADChem_Cor))
cor.cols<-8:24 #save indices of columns to correct
for(ii in 1:nrow(ADChem_Cor)){
  same.date<-ADChem_BLANK$ConfDate==ADChem_Cor$ConfDate[ii]
  if(sum(same.date)>0){
    ADChem_Cor[ii, cor.cols]<- ADChem_Cor[ii, cor.cols] - colMeans(ADChem_BLANK[same.date, cor.cols], na.rm=T)
    ADChem_Cor$cor[ii]<-1
  }
}

## Adjust site IDs ------------------------------------------------------------
ADChem_Cor$SiteID[ADChem_Cor$SiteID=="040747A"]<-"40747"
ADChem_Cor$SiteID[ADChem_Cor$SiteID=="040750A"]<-"40750"

# unique(ADChem_Cor$SiteID)
# Reformatting for coordination's sake NEEDS TO BE DONE BEFORE TESTING CLUSTERING FUNCTION ####
ADSites_Corr<-sort(unique(as.character(ADChem_Cor$SiteID)))

locstrue<-read.csv("./data/lake_coords.csv",stringsAsFactors = F)
locstrue<-locstrue[order(locstrue$SiteID),]

limelocs<-c("30172", "40905", "60182", "40576")

coordsSet<-data.frame(NA, nrow = 53, ncol = 3)
colnames(coordsSet)<-c("Site","X","Y")

for(i in 1:53){
  for(j in 1:3){
    coordsSet[i,j]<-locstrue[i,j]
  }
}

coords<-coordsSet[,c(2,3)]
coords<-coords[-which(coordsSet[,1] %in% limelocs),]
#print(coords)

ADSites_Corr<-coordsSet[-which(coordsSet[,1] %in% limelocs),1]

## write function to create matrices and replace NAs with monthly medians

makeDataMatrix<-function(df, variable){
  
  subdat<-df[,c("ConfDate","SiteID","month","month_abs",variable)]
  
  outmat<-matrix(NA, nrow=length(unique(df$SiteID)), ncol=length(unique(df$ConfDate)))
  sites<-sort(unique(df$SiteID))
  times<-sort(unique(df$ConfDate))
  #organize into matrix
  for(site in 1:length(sites)){
    for(tt in 1:length(times)){
      outmat[site,tt]<-mean(subdat[,variable][subdat$SiteID==sites[site] & subdat$ConfDate==times[tt]], na.rm=T)
    }
  }
  rownames(outmat)<-sites
  colnames(outmat)<-times
  
  #replace NAs with monthly median
  months<-substr(colnames(outmat),6,7)
  for(site in 1:length(sites)){
    for(tt in 1:length(times)){
      if(is.na(outmat[site,tt])){
        outmat[site,tt]<-median(outmat[site, months==months[tt]], na.rm=T)
      }
    }
  }
  return(outmat)
}
# Editing ADChem_Cor to drop limed sites
ADChem_Cor<-ADChem_Cor[-which(ADChem_Cor$SiteID %in% limelocs),]
ADChem_Cor$DOC[ADChem_Cor$DOC < 0.5] <- NA
ADChem_Cor$DOC[ADChem_Cor$DOC > 25] <-NA
# Convention altered to include Mat so Fluorine could continue to use F
DOCMat<-makeDataMatrix(ADChem_Cor, variable="DOC")
clDOC<-cleandat(DOCMat,1:ncol(DOCMat),clev=5)$cdat

#write.csv(coords, "./Reformatted data/lake_coords_used.csv", row.names=FALSE)

## get distance decay

correlog <- Sncf(x=coords$X, y=coords$Y, z=clDOC, latlon=TRUE)
#plot(correlog)

tiff("./Graphics/FigS1_Correlogram.tif", units="in", res=300, width=3.5, height=3)

par(mar=c(3,3,1.1,0.5), mgp=c(1.8,0.5,0), tcl=-0.3)

plot(correlog$real$predicted$x, correlog$real$predicted$y, type="l", lwd=2, ylim=c(-1,1),
     xlab="Distance (km)", ylab="Synchrony (Pearson correlation)")
lines(correlog$boot$boot.summary$predicted$x, correlog$boot$boot.summary$predicted$y[2,], lty=2)
lines(correlog$boot$boot.summary$predicted$x, correlog$boot$boot.summary$predicted$y[10,], lty=2)
abline(h=0, lty=3, col="grey")

dev.off()




## make cluster maps

plotClusterMap <- function(inclust,spltlvl=length(inclust$clusters),basemap=NULL,nodewgt=NULL,
                           nodesize=c(1,3),nodecol=brewer.pal(9,"Set1"), edgewgt=NULL, edgesize=c(0.1,1),
                           addlegend=TRUE, filename=NA)
{
  
  #some checking of validity of inputs
  if(!is.null(nodecol)){
    if(length(unique(nodecol)) < max(unlist(inclust$clusters[spltlvl]))){
      stop("more clusters than colors in nodecol")
    }
  }
  
  #convert inclust$coords to common format
  if(all(c("X","Y") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$X,
                                                                    Y=inclust$coords$Y)}
  if(all(c("lat","lon") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$lon,
                                                                        Y=inclust$coords$lat)}
  if(all(c("latitude","longitude") %in% names(inclust$coords))){coords<-data.frame(X=inclust$coords$longitude,
                                                                                   Y=inclust$coords$latitude)}
  
  
  
  #expand right side margin if a legend is being added
  if(addlegend){
    par.mar<-par("mar")
    mar.new<-par.mar
    mar.new[4]<-6.1
    par(mar=mar.new,xpd=T)
  }
  
  if(!is.null(basemap)){
    #give a warning if basemap extent does not contain all points
    ext<-bbox(basemap)
    if(any(c(min(coords$X)<ext[1,1],
             max(coords$X)>ext[1,2],
             min(coords$Y)<ext[2,1],
             max(coords$Y)>ext[2,2]))
    ){
      warning("basemap extent does not contain all nodes; possible coordinate system mismatch?")
    }
    plot(basemap, axes=T)
  }
  
  if(is.null(basemap)){
    plot(coords$X,coords$Y,pch=16,col=NA,xlab="X",ylab="Y")
  }
  
  x0<-rep(coords$X, each=nrow(coords))
  y0<-rep(coords$Y, each=nrow(coords))
  x1<-rep(coords$X, times=nrow(coords))
  y1<-rep(coords$Y, times=nrow(coords))
  
  if(!is.null(edgewgt)){
    if(length(edgewgt)!=1 & length(edgewgt)!=length(inclust$adj)){
      stop("invalid 'edgewgt' argument")
    }
    if(all(dim(edgewgt)==dim(inclust$adjmat))){
      edgewgts<-edgewgt
    }
    if(edgewgt=="adj"){
      ewgt<-inclust$adj-min(inclust$adj,na.rm=T)
      ewgt<-ewgt/max(ewgt,na.rm=T)
      edgewgts<-ewgt*(edgesize[2]-edgesize[1]) + edgesize[1]
    }
    if(length(edgewgt)==1 & is.numeric(edgewgt)){
      wgtthresh<-quantile(inclust$adj,edgewgt,na.rm=T)
      edgewgts<-rep(0,prod(dim(inclust$adj)))
      edgewgts[c(inclust$adj)>wgtthresh]<-1
    }
    
    segments(x0,y0,x1,y1,lwd=edgewgts)
  }
  
  if(is.null(nodewgt)){
    points(coords[,1], coords[,2], pch=21, cex=nodesize[1], bg=nodecol[unlist(inclust$clusters[spltlvl])])
  }
  else{ 
    if(nodewgt=="mod.decomp"){
      membwgt<-inclust$modres[[spltlvl]]$nodeQ-min(inclust$modres[[spltlvl]]$nodeQ)
      membwgt<-membwgt/max(membwgt)
      nodecex<-membwgt*(nodesize[2]-nodesize[1]) + nodesize[1]
    }
    if(length(nodewgt)==nrow(coords)){
      membwgt<-nodewgt-min(nodewgt)
      membwgt<-membwgt/max(membwgt)
      nodecex<-membwgt*(nodesize[2]-nodesize[1]) + nodesize[1]
    }
    points(coords[,1], coords[,2], pch=16, cex=nodecex, col=nodecol[unlist(inclust$clusters[spltlvl])])
  }
  
  if(addlegend){
    legx<-par('usr')[2] + 0.01*abs(diff(par('usr')[1:2]))
    legy1<-par('usr')[4]
    leg1<-legend(legx,legy1,legend=paste0("Module ",1:max(unlist(inclust$clusters[spltlvl]))), pch=16, col=
                   nodecol[1:max(unlist(inclust$clusters[spltlvl]))],title="Membership",bty="n")
    if(!is.null(nodewgt)){
      legy2<-legy1 - leg1$rect$h
      labs=round(c(min(inclust$modres[[spltlvl]]$nodeQ,na.rm=T),
                   mean(inclust$modres[[spltlvl]]$nodeQ,na.rm=T),
                   max(inclust$modres[[spltlvl]]$nodeQ,na.rm=T)),digits=3)
      sizes=c(min(nodecex),mean(nodecex),max(nodecex))
      leg2<-legend(legx,legy2,legend=labs,pt.cex=sizes,pch=1,title="Node weight",bty="n")
    }
    if(!is.null(edgewgt)){
      if(edgewgt=="adj" | is.matrix(edgewgt)){
        if(is.null(nodewgt)){legy3<-legy1-leg1$rect$h}
        else{legy3<-legy1-leg1$rect$h-leg2$rect$h}
        labs2=round(c(min(inclust$adj,na.rm=T),mean(inclust$adj,na.rm=T),
                      max(inclust$adj,na.rm=T)),digits=3)
        wdths=c(min(edgewgts,na.rm=T),mean(edgewgts,na.rm=T),max(edgewgts,na.rm=T))
        legend(legx,legy3,legend=labs2,lty=1,lwd=wdths,title="Edge weight",bty="n")
      }
    }
    par(mar=par.mar) #reset 'mar' graphics parameter
  }
  
  
}

spatclust <- clust(clDOC, times=1:ncol(DOCMat), coords=coords, method="pearson")

sampleDate <- as.Date(paste(substr(ADMonths_abs,1,4), substr(ADMonths_abs,6,7),"01", sep="-"))

shp <- readOGR("/Users/jonathanwalter/Library/CloudStorage/GoogleDrive-jaw3es@virginia.edu/.shortcut-targets-by-id/1fFka4nHJxFqxG2JGNeDJSZA8I6brn51a/ADK DIC (Shared)/Original Data/GIS data/AdirondackParkBoundary2017/AdirondackParkBoundary2017.shp")
shp <- spTransform(shp, CRS("+init=EPSG:4326"))


cols<-brewer.pal(9,"Set1")

tiff("./Graphics/FigS2_ClusterMap.tif", units="in", res=300, width=6.5, height=4.5)

layout(matrix(c(1,1,2,3), 2,2), widths=c(0.6,0.4))
par(mar=c(2,2,1.1,0.5), mgp=c(1.8,0.5,0), tcl=-0.3)

plotClusterMap(spatclust, basemap=shp, addlegend=FALSE)
mtext("a)",at=-75.56)

par(mar=c(3.1,3.1,1.1,1.1), mgp=c(1.8,0.5,0), tcl=-0.3)
plot(sampleDate, spatclust$mns[[1]][1,], type="l", xlab="", ylab="Module 1", col=cols[1])
mtext("b)",at=as.numeric(min(sampleDate))-800)
plot(sampleDate, spatclust$mns[[2]][2,], type="l", xlab="", ylab="Module 2", col=cols[2])
mtext("c)",at=as.numeric(min(sampleDate))-800)

dev.off()


write.csv(data.frame(siteid=ADSites_Corr, cluster=spatclust$clusters[[2]]), "clusters.csv", row.names = FALSE)


#' Assign ID numbers to samples
#'
#' \code{assignIDs} attempts to assign ID numbers to contiguous regions.
#'
#' @param data The dataset containing x- and y-coordinates.
#' @param xCoord A character string specifying the name of the x-coordinate
#' variable.
#' @param yCoord A character string specifying the name of the y-coordinate
#' variable.
#' @param neighborhood A character string specifying the neighborhood
#' structure to be considered.  The options are "primary", "secondary", and
#' "tertiary".
#'
#' @return The dataset appended with a variable called "assignID".  This
#' variable provides the assigned IDs as a numeric variable.
#'
#' @examples
#' data("TAMdata")
#' TAMdata <- assignIDs(TAMdata)

assignIDs<-function(data, xCoord='x', yCoord='y', neighborhood='tertiary'){
  data$xnew<-data[[xCoord]]-min(data[[xCoord]])
  data$xnew<-round(data$xnew)

  data$ynew<-data[[yCoord]]-min(data[[yCoord]])
  data$ynew<-round(data$ynew)

  mindist<-min(dist(cbind(data$xnew,data$ynew)))

  data$xnewer<-(data$xnew)/mindist
  data$ynewer<-(data$ynew)/mindist

  # create neighborhood structure
  data$indexForRast<-seq(1:nrow(data))

  if(neighborhood=='primary'){
    neighbors<-list()
    for(i in 1:nrow(data)){
      neighbors[[i]]<-data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+mindist]
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]+mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]+mindist])
    }
  }


  if(neighborhood=='secondary'){
    neighbors<-list()
    for(i in 1:nrow(data)){
      # primary
      neighbors[[i]]<-data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+mindist]
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]+mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]+mindist])

      # secondary
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]+(2*mindist)])

    }
  }

  if(neighborhood=='tertiary'){
    neighbors<-list()
    for(i in 1:nrow(data)){
      # primary
      neighbors[[i]]<-data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+mindist]
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]+mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]-mindist])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-mindist & data$ynew==data$ynew[i]+mindist])

      # secondary
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+mindist & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]+(2*mindist)])

      # tertiary
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(mindist) & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(3*mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(2*mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]+(mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i] & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]-(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]-(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]-(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]+(mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]+(2*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(3*mindist) & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(2*mindist) & data$ynew==data$ynew[i]+(3*mindist)])
      neighbors[[i]]<-c(neighbors[[i]], data$indexForRast[data$xnew==data$xnew[i]-(mindist) & data$ynew==data$ynew[i]+(3*mindist)])

    }
  }



  tempdat<-data
  neighbors.index<-tempdat$indexForRast
  iter<-1
  sampids<-list()
  flag.outer<-1
  while(flag.outer){
    nrast.old<-length(neighbors[[neighbors.index[1]]])
    rast.curr<-neighbors[[neighbors.index[1]]]
    flag<-1
    while(flag){
      for(i in 1:length(rast.curr)){
        rast.curr<-c(rast.curr,neighbors[[rast.curr[i]]])
      }

      rast.curr<-unique(rast.curr)
      nrast.curr<-length(rast.curr)

      flag<-ifelse(nrast.curr==nrast.old,0,1)
      rast.old<-rast.curr
      nrast.old<-length(rast.old)
    }

    tempdat<-tempdat[!(tempdat$indexForRast %in% rast.curr),]
    neighbors.index<-neighbors.index[!(neighbors.index %in% rast.curr)]

    sampids[[iter]]<-rast.curr

    iter<-iter+1

    flag.outer<-ifelse(nrow(tempdat)==0,0,1)
  }

  data$assignID<-NA
  for(i in 1:length(sampids)){
    data$assignID[data$indexForRast %in% sampids[[i]]]<-i
  }

  data<-subset(data, select=-c(xnew,ynew,xnewer,ynewer,indexForRast))




}










#' Plot the assigned ID numbers
#'
#' \code{plotIDs} plots the samples and overlays the sample IDs.
#'
#' @param data The dataset containing x- and y-coordinates and a sample
#' identifier.
#' @param xCoord A character string specifying the name of the x-coordinate
#' variable.
#' @param yCoord A character string specifying the name of the y-coordinate
#' variable.
#' @param IDVar A character string speciffying the sample ID variable.
#' @param textSize The size of the sample IDs on the plot.
#' @param TMA A TRUE/FALSE variable specifying whether or not the data came
#' from a tissue microarray (TMA).  This is only used to determine how to
#' color the samples.  If \code{TMA=TRUE}, then the samples are colored using
#' a continuous color scale.  If \code{TMA=FALSE}, then the samples are
#' colored using a discrete color scale.
#'
#' @return A plot of the tissue samples colored according to the sample ID and
#' overlaid with the ID.
#'
#' @examples
#' data("TAMdata")
#' TAMdata <- assignIDs(TAMdata)
#' plotIDs(TAMdata)

plotIDs<-function(data, xCoord='x', yCoord='y', IDVar='assignID', textSize=4, TMA=FALSE){
  annoplot<-matrix(nrow=length(unique(data[[IDVar]])), ncol = 3)
  annoplot<-as.data.frame(annoplot)
  colnames(annoplot)<-c('sampid','meanx','meany')
  annoplot$sampid<-unique(data[[IDVar]])
  for(i in 1:nrow(annoplot)){
    annoplot$meanx[i]<-mean(data[[xCoord]][data[[IDVar]]==annoplot$sampid[i]])
    annoplot$meany[i]<-mean(data[[yCoord]][data[[IDVar]]==annoplot$sampid[i]])
  }

  if(TMA==TRUE){
    ggplot(data, aes(x=x, y=y)) +
      geom_point(aes_string(color=IDVar)) +
      scale_color_continuous(low='lightblue', high='navy') +
      annotate('text',x=annoplot$meanx, y=annoplot$meany, label=annoplot$sampid, size=textSize)
  }

  if(TMA==FALSE){
    data[[IDVar]]<-as.factor(data[[IDVar]])
    ggplot(data, aes(x=x, y=y)) +
      geom_point(aes_string(color=IDVar)) +
      annotate('text',x=annoplot$meanx, y=annoplot$meany, label=annoplot$sampid, size=textSize)
  }

}



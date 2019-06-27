#' Plot a range object.
#'
#' \code{plotIDs} the semivariance estimates, covariance model, and number of
#' pairs of observations used to estimate semivariance.
#'
#' @param rangeObj An object of class \code{rangeList}.  This object is a list
#' that contains the estimated range parameter(s).  This is obtained by a call
#' to the \code{\link{estRange}} function.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' plotRangeObj(rangs)
#'
#' @return A set of plots displaying the semivariance estimates, fitted
#' covariance model, and number of pairs of observations used to estimate
#' semivariance.  If no spatial variable was provided to the
#' \code{\link{estRange}} function, then there will be two plots.  The first,
#' on top, will show the semivariance estimates and the fitted covariance
#' model.  The second, on the bottom, will show the number of data pairs used
#' to the estimate the semivariance at the corresponding distance.  If a
#' spatial variable is supplied to the \code{\link{estRange}} function, then
#' there will be two plots for every level of the spatial variable.



plotRangeObj<-function(rangeObj){

  if(is.null(rangeObj[['spatialVar']])==TRUE){
    semiline<-data.frame(rangeObj$semivarFit$u,rangeObj$semivarFit$v,rangeObj$semivarFit$n)
    colnames(semiline)<-c('dist','semiv','num')

    covfit<-matrix(nrow=1000, ncol=2)
    covfit<-as.data.frame(covfit)
    colnames(covfit)<-c('dist','fitval')
    covfit[,1]<-seq(0,max(rangeObj$semivarFit$u), length.out = nrow(covfit))
    for(i in 1:nrow(covfit)){
      covfit[i,2]<-rangeObj$covModelFit$nugget + ((rangeObj$estSig2*(1-exp(-(covfit$dist[i]/rangeObj$estRange)**2))))
    }

    p1<-ggplot(semiline, aes(x=dist, y=semiv)) +
      geom_point() +
      geom_line(data=covfit, aes(x=dist, y=fitval), color='red') +
      ylim(0,max(semiline$semiv)) +
      labs(x='distance',y='semivariance') +
      theme_bw()


    p2<-ggplot(data=semiline, aes(x=dist,y=num)) +
      geom_area(data=semiline, aes(x=dist,y=num), fill='darkblue') +
      xlim(0,max(semiline$dist)) +
      labs(x='distance',y='number of data pairs') +
      theme_bw()

    grid.arrange(p1,p2, ncol=1)
  }

  if(is.null(rangeObj[['spatialVar']])==FALSE){

    semiline<-list()
    for(i in 1:length(rangeObj$estRange)){
      semiline[[i]]<-data.frame(rangeObj$semivarFit[[i]]$u,rangeObj$semivarFit[[i]]$v,rangeObj$semivarFit[[i]]$n)
      colnames(semiline[[i]])<-c('dist','semiv','num')
    }

    covfit<-list()
    for(i in 1:length(rangeObj$estRange)){
      covfit[[i]]<-matrix(nrow=1000, ncol=2)
      covfit[[i]]<-as.data.frame(covfit)
      colnames(covfit[[i]])<-c('dist','fitval')
      covfit[[i]][,1]<-seq(0,max(rangeObj$semivarFit[[i]]$u), length.out = nrow(covfit[[i]]))

      for(j in 1:nrow(covfit[[i]])){
        covfit[[i]][j,2]<-rangeObj$covModelFit[[i]]$nugget + ((rangeObj$estSig2[i]*(1-exp(-(covfit[[i]]$dist[j]/rangeObj$estRange[i])**2))))
      }
    }

    spatvarlevels<-sort(unique(rangeObj$data[[rangeObj$spatialVar]]))

    for(i in 1:length(rangeObj$estRange)){
      eval(parse(text = c(paste0('p',i,'<-','ggplot(semiline[[',i,']], aes(x=dist, y=semiv)) +'),
                          'geom_point() +',
                          paste0('geom_line(data=covfit[[',i,']], aes(x=dist, y=fitval), color="red") +'),
                          paste0('ylim(0,max(semiline[[',i,']]$semiv)) +'),
                          paste0('labs(x="distance",y="semivariance", title="',rangeObj$spatialVar,'=',spatvarlevels[i],'") +'),
                          'theme_bw()')))

      eval(parse(text = c(paste0('p',(length(rangeObj$estRange)+i),'<-','ggplot(semiline[[',i,']], aes(x=dist, y=num)) +'),
                          paste0('geom_area(data=semiline[[',i,']], aes(x=dist, y=num), fill="darkblue") +'),
                          paste0('xlim(0,max(semiline[[',i,']]$dist)) +'),
                          'labs(x="distance",y="number of data pairs") +',
                          'theme_bw()')))
    }

    eval(parse(text = paste0('grid.arrange(',paste('p',seq(1:(length(rangeObj$estRange)*2)), sep='', collapse=', '),', nrow=2, heights=c(1.25,1))')))

  }

}

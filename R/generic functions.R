

#' Summarize rangeList object
#'
#' \code{summary.rangeList} provides basic results for an object of class
#' \code{rangeList}.
#'
#' @param object An object of class \code{rangeList}.  This is the result of a
#' call to the \code{estRange} function.
#'
#' @return The \code{summary.rangeList} prints the estimated range for
#' every level of the spatial variable provided to the \code{estRange}
#' function.  If no spatial varaible was provided, then only one estimated
#' range is printed.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' summary(rangs)

summary.rangeList<-function(object){
  spatialVar<-object[['spatialVar']]
  if(is.null(object[['spatialVar']])==TRUE){
    cat('The estimated range is ',object[['estRange']],'.', sep='')
  }

  if(is.null(object[['spatialVar']])==FALSE){
    singlerange<-rep(NA,length(object[['estRange']]))
    for(i in 1:length(singlerange)){
      singlerange[i]<-paste0('The estimated range for the "',names(object[['estRange']][i]),'" level of ',object$spatialVar, ' is ',round(object[['estRange']][i], digits=5),'.')
    }
    cat(singlerange, sep='\n')
  }
}



#' Summarize structureList object
#'
#' \code{summary.structureList} provides basic results for an object of class
#' \code{structureList}.
#'
#' @param object An object of class \code{structureList}.  This is the result
#' of a call to the \code{chooseStructures} function.
#'
#' @return For every level a spatial variable (provided in the \code{estRange}
#' function), the \code{summary.rangeList} prints the chosen support structure,
#' the total number of support points across all samples after removing
#' support points, and the number of samples.  If no spatial variable was
#' provided, only one line of information is printed.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' structs <- chooseStructures(rangs)
#' summary(structs)

summary.structureList<-function(object){
  spatialVar<-object[['spatialVar']]
  if(is.null(spatialVar)==TRUE){
    cat('The ',as.character(object[[1]]$recommendStruct), ' support structure was chosen.\n  After removing support sites, there were ',
        object[[1]]$nSupportReduced, ' support sites across\n  ', length(unique(object[[1]]$coordsU$sample)), ' samples.', sep='')
  }

  if(is.null(spatialVar)==FALSE){
    data <- object[['data']][order(object[['data']][[spatialVar]]),]
    spatvarlevels <- unique(data[[spatialVar]])

    singlestruct<-rep(NA,length(unique(object[['data']][[spatialVar]])))
    for(i in 1:length(singlestruct)){
      singlestruct[i]<-paste0(spatialVar, ' level ', '"', spatvarlevels[i], '": The ',as.character(object[[i]]$recommendStruct), ' support structure was chosen.\n  After removing support sites, there were ',
                              object[[i]]$nSupportReduced, ' support sites across\n  ', length(unique(object[[i]]$coordsU$sample)), ' samples.')
    }
    cat(singlestruct, sep='\n')
  }

}



#' Summarize PCResults object
#'
#' \code{summary.PCResults} provides basic results for an object of class
#' "PCResults".
#'
#' @param object An object of class "PCResults".  This is the result of a
#' call to the \code{runPCModel} function.
#'
#' @return The \code{summary.PCResults} function prints the matrix of
#' summary measures labeled "results" in the object.  This matrix contains
#' summary information for the model coefficients.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' structs <- chooseStructures(rangs)
#' PCdat <- createPCData(structs, trimData = FALSE,
#'                       covariates = c("secondary", "TAM", "secTAM"),
#'                       covariateTypes = c("binary", "binary", "binary"),
#'                       covariateLevels = c("sample", "raster", "raster"))
#' PCmod <- writePCModel(PCdat, multiSampsPerSubj = TRUE, typeOfZero = "censored")
#' PCresults <- runPCModel(modelObj = PCmod, PCDataObj = PCdat, slideVar='slide',
#'                         monitorCoefOnly = FALSE,
#'                         nBurnin = 15000, nIter = 40000, nThin = 25)
#' summary(PCresults)

summary.PCResults<-function(object){
  print(object[['results']])
}



#' Plot the semivariance estimates overlaid with the covariance model
#'
#' \code{plot.rangeList} uses the S3 plot methods in the \code{geoR} package
#' to plot the semivariance estimates and then overlay the fitted covariance
#' model.
#'
#' @param object An object of class "rangeList".
#'
#' @return A plot or plots or the semivariance estimates overlaid with the
#' fitted covariance model.  If there was no spatial variable provided in the
#' \code{estRange} function then only a single plot is produced.  This is
#' equivalent to the plots produced by \code{geoR}.  However, if a spatial
#' variable was provided, then \code{plot.rangeList} produces a plot for
#' every level of that variable.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' plot(rangs)
#'
#' @references Cressie, N and Hawkins, DM. 1980. Robust estimation of the
#' variogram: I. \emph{Journal of the International Association for
#' Mathematical Geology}, 12(2):115-125.

plot.rangeList<-function(object){
  if(is.null(object[['spatialVar']])==TRUE){
    plot(object[['semivarFit']])
    lines(object[['covModelFit']], col='red')
  }

  if(is.null(object[['spatialVar']])==FALSE){
    spatialVar<-object[['spatialVar']]
    spatvarlevels <- unique(object[['data']][[spatialVar]])
    nvarlevels<-length(unique(spatvarlevels))
    sampleVar<-object[['sampleVar']]
    object[['data']][[spatialVar]]<-as.factor(object[['data']][[spatialVar]])

    dimrow<-ceiling(sqrt(nvarlevels))

    flag<-1
    iter<-0
    while(flag){
      iter<-iter+1
      flag<-ifelse((dimrow*iter)>=nvarlevels,0,1)
    }

    par(mfrow=c(dimrow,iter),
        oma = c(2,4,1,1) + 0.1,
        mar = c(3,1,1,1) + 0.1)

    for(i in 1:nvarlevels){
      plot(object[['semivarFit']][[i]], main=paste0(spatialVar,'=',spatvarlevels[i]))
      mtext('semivariance',2, outer=TRUE, padj=-3)
      mtext('distance',1, outer=TRUE, padj=0)
      lines(object[['covModelFit']][[i]], col='red')
    }

  }
}



#' Overlay the tissue samples with corresponding support structures
#'
#' \code{plotPCStructure} plots the tissue sample rasters and overlays them
#' with the support structure(s) used in downstream modeling.
#'
#' @param structureObj An object of class \code{structureList}.
#' @param removeBackground A TRUE/FALSE argument specifying if the background,
#' including the x- and y-axes, should be removed.
#' @param marginProp An argument specifying the size of the margins.  It is
#' a proportion of the largest distance from the origin to the outermost
#' support sites.  This is used to help separate plots.
#' @param titleSize The font size of each title specifying the sample number.
#' @param supSiteSize The size of the support point dots.
#' @param supSiteLegendSize The size of the support sites in the legend.
#' @param rastLegendSize The size of the rasters in the legend.
#' @param rasterColors The color of the rasters.  This is a concatenated
#' character string that must have at least two elements.  Each element
#' specifies an accepted color in R.  The \code{plotPCStructure} function uses
#' the \code{colorRampPalette} function, which generates smooth transitions
#' between the colors provided.  If the number of colors provided is equal to
#' the number of levels of the spatial variable, then the plotted colors will
#' be the same as those provided.  However, if the number of colors provided
#' does not equal the number of levels of the spatial variable, then the
#' plotted colors will be interpolated from those provided.
#' @param supSiteColors The color of the support sites.  This is a
#' concatenated character string that must have at least two elements.  The
#' colors are generated in the same way as the \code{rasterColors}.
#'
#' @return A gridded plot of all tissue samples.  If a spatial variable was
#' provided in the \code{estRange} function, then each sample will be colored
#' by the levels of that variable.  Using a separate set of colors, the
#' support sites will also be colored according to the levels of the spatial
#' variable.  If no spatial variable was provided, then both the sample
#' rasters and the support sites will each be represented by a single color.
#' The last cell in the grid of plots is the legend, showing which colors are
#' used to represent the rasters and support sites.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' structs <- chooseStructures(rangs)
#' plotPCStructure(structs, supSiteSize = 1.5, marginProp = 0.25)


plotPCStructure<-function(structureObj, removeBackground=TRUE, marginProp=0.1, titleSize=11,
                          supSiteSize=3, supSiteLegendSize=5, rastLegendSize=30,
                          rasterColors=c('gray75','gray20'), supSiteColors=c('royalblue','red'),
                          plotCols=NULL){
  if (class(structureObj) != "structureList")
    stop("structureObj must of class structureList.  This object must be created using the chooseStructures function.")

  spatialVar<-structureObj[['spatialVar']]
  sampleVar<-structureObj[['sampleVar']]
  # no spatialVar
  if(is.null(structureObj[['spatialVar']])==TRUE){
    if(structureObj[[1]]$recommendStruct=='noStructure')
      stop('There are no support sites to plot.')
    samplelabs<-unique(supmat$sample)

    minxsup<-min(structureObj[[1]][['coordsU']]$x.omega)
    maxxsup<-max(structureObj[[1]][['coordsU']]$x.omega)
    minysup<-min(structureObj[[1]][['coordsU']]$y.omega)
    maxysup<-max(structureObj[[1]][['coordsU']]$y.omega)
    xvals<-c(0,0,(minxsup+(minxsup*marginProp)), (maxxsup+(maxxsup*marginProp)))
    yvals<-c((minysup+(minysup*marginProp)), (maxysup+(maxysup*marginProp)),0,0)
    dummy<-data.frame(xvals,yvals)

    for(i in 1:length(samplelabs)){
      if(removeBackground==TRUE){
        eval(parse(text = c(paste0('ggplot(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes(x=xPC, y=yPC)) +'),
                            paste0('geom_raster(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],]) +'),
                            paste0('geom_point(data=structureObj[[1]][["coordsU"]][structureObj[[1]][["coordsU"]]$sample==samplelabs[',i,'],], aes(x=x.omega, y=y.omega), color=supSiteColors[1], size=supSiteSize) +'),
                            'scale_fill_manual(values = rasterColors[1], guide=FALSE) +',
                            paste0('labs(title=paste0("Sample ",samplelabs[',i,']))'),
                            'geom_blank(data=dummy, aes(x=xvals, y=yvals)) +',
                            'theme_void()')))

      }
      if(removeBackground==FALSE){
        eval(parse(text = c(paste0('ggplot(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes(x=xPC, y=yPC)) +'),
                            paste0('geom_raster(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],]) +'),
                            paste0('geom_point(data=structureObj[[1]][["coordsU"]][structureObj[[1]][["coordsU"]]$sample==samplelabs[',i,'],], aes(x=x.omega, y=y.omega), color=supSiteColors[1], size=supSiteSize) +'),
                            'scale_fill_manual(values = rasterColors[1], guide=FALSE) +',
                            paste0('labs(title=paste0("Sample ",samplelabs[',i,']))'),
                            'geom_blank(data=dummy, aes(x=xvals, y=yvals)) +')))
      }
    }
    plotvec<-paste0('p', seq(1:length(samplelabs)))

    nplots<-length(samplelabs)

    if(is.null(plotCols)==TRUE){dimcol<-ceiling(sqrt(nplots))}
    if(is.null(plotCols)==FALSE){dimcol<-plotCols}


    eval(parse(text = paste0('grid.arrange(',paste(plotvec, collapse = ', '),', ncol=dimcol)')))
  }

  # spatialVar
  if(is.null(structureObj[['spatialVar']])==FALSE){
    spatvarlevels <- unique(structureObj[['data']][[spatialVar]])
    nvarlevels<-length(unique(spatvarlevels))
    recstructures<-rep(NA,nvarlevels)
    for(i in 1:nvarlevels){
      recstructures[i]<-as.character(structureObj[[i]]$recommendStruct)
    }
    if(all(recstructures=='noStructure'))
      stop('There are no support sites to plot.')
    supsitestoplot<-which(recstructures!='noStructure')

    sampleVar<-structureObj[['sampleVar']]
    structureObj[['data']][[spatialVar]]<-as.factor(structureObj[['data']][[spatialVar]])

    supl<-list()
    for(i in supsitestoplot){
      supl[[i]]<-structureObj[[i]][['coordsU']]
      supl[[i]][[spatialVar]]<-spatvarlevels[i]
    }
    supmat<-rbind.fill(supl)
    supmat[[spatialVar]]<-as.factor(supmat[[spatialVar]])


    # need color lists to make sure all plot colors match
    samplelabs<-unique(supmat$sample)
    samplelabs<-unique(structureObj[['data']][[sampleVar]])
    rastcol<-colorRampPalette(rasterColors)
    rastcolvec<-rastcol(nvarlevels)
    supcol<-colorRampPalette(supSiteColors)
    supcolvec<-supcol(nvarlevels)
    supcolvec2<-supcolvec[supsitestoplot]
    colorls<-list()
    rastercolorl<-list()
    supportcolorl<-list()
    for(i in 1:length(unique(samplelabs))){
      colorls[[i]]<-which(spatvarlevels %in% supmat[[spatialVar]][supmat$sample==samplelabs[i]])
      supportcolorl[[i]]<-supcolvec[colorls[[i]]]
    }

    colorlr<-list()
    colorrastl<-list()
    for(i in 1:length(unique(samplelabs))){
      colorlr[[i]]<-which(spatvarlevels %in% structureObj[['data']][[spatialVar]][structureObj[['data']][[sampleVar]]==samplelabs[i]])
      rastercolorl[[i]]<-rastcolvec[colorlr[[i]]]
    }

    minxsup<-min(supmat$x.omega)
    maxxsup<-max(supmat$x.omega)
    minysup<-min(supmat$y.omega)
    maxysup<-max(supmat$y.omega)
    xvals<-c(0,0,(minxsup+(minxsup*marginProp)), (maxxsup+(maxxsup*marginProp)))
    yvals<-c((minysup+(minysup*marginProp)), (maxysup+(maxysup*marginProp)),0,0)
    dummy<-data.frame(xvals,yvals)

    for(i in 1:length(samplelabs)){
            if(removeBackground==TRUE){
      eval(parse(text=c(paste0('p',i,'<-ggplot(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes(x=xPC, y=yPC)) +'),
                        paste0('geom_raster(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes_string(fill=spatialVar)) +'),
                        paste0('geom_point(data=supmat[supmat$sample==samplelabs[',i,'],], aes_string(x="x.omega", y="y.omega", color=spatialVar), size=supSiteSize) +'),
                        paste0('scale_color_manual(values = supportcolorl[[',i,']], guide=FALSE) +'),
                        paste0('scale_fill_manual(values = rastercolorl[[',i,']], guide=FALSE) +'),
                        paste0('labs(title=paste0("Sample ",samplelabs[',i,'])) +'),
                        paste0('geom_blank(data=dummy, aes(x=xvals, y=yvals)) +'),
                        'theme_void() +',
                        'theme(plot.title = element_text(size=titleSize))')))
      }
      if(removeBackground==FALSE){
        eval(parse(text=c(paste0('p',i,'<-ggplot(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes(x=xPC, y=yPC)) +'),
                          paste0('geom_raster(data=structureObj[["data"]][structureObj[["data"]][[sampleVar]]==samplelabs[',i,'],], aes_string(fill=spatialVar)) +'),
                          paste0('geom_point(data=supmat[supmat$sample==samplelabs[',i,'],], aes_string(x="x.omega", y="y.omega", color=spatialVar), size=supSiteSize) +'),
                          paste0('scale_color_manual(values = supportcolorl[[',i,']], guide=FALSE) +'),
                          paste0('scale_fill_manual(values = rastercolorl[[',i,']], guide=FALSE) +'),
                          paste0('labs(title=paste0("Sample ",samplelabs[',i,'])) +'),
                          paste0('geom_blank(data=dummy, aes(x=xvals, y=yvals)) +'),
                          'theme(plot.title = element_text(size=titleSize))')))
      }
    }
    plotvec<-paste0('p', seq(1:length(samplelabs)))

    nplots<-length(samplelabs)+1

    if(is.null(plotCols)==TRUE){dimcol<-ceiling(sqrt(nplots))}
    if(is.null(plotCols)==FALSE){dimcol<-plotCols}


    pblank<-ggplot(data=structureObj[['data']], aes(x=x, y=y)) +
      geom_raster(data=structureObj[['data']], aes_string(fill=spatialVar)) +
      #facet_wrap(as.formula(paste("~", sampleVar))) +
      geom_point(data=supmat, aes_string(x='x.omega', y='y.omega', color=spatialVar), size=3) +
      scale_color_manual(values = supcolvec2) +
      scale_fill_manual(values = rastcolvec) +
      guides(colour = guide_legend(override.aes = list(size=supSiteLegendSize))) +
      theme_void() +
      theme(legend.key.size = unit(rastLegendSize, 'pt'))

    justleg<-cowplot::get_legend(pblank)

    eval(parse(text = paste0('grid.arrange(',paste(plotvec, collapse = ', '),',justleg, ncol=dimcol)')))
  }


}






#' Plot the MCMC iterations for each model coefficient.
#'
#' \code{plot.PCResults} generates the MCMC iteration plots for model
#' coefficients after running a PC model using \code{runPCModel}.
#'
#' @param object An object of class \code{PCResults}.  This is the result of
#' running a model using \code{runPCModel}.
#'
#' @return A plot of the MCMC chains for each model coefficient.
#'
#' @examples
#' data("TAMdata")
#' # The dataset is trimmed only for the speed of the example
#' TAMdata <- TAMdata[TAMdata$subject < 3, ]
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)
#' structs <- chooseStructures(rangs)
#' PCdat <- createPCData(structs, trimData = FALSE,
#'                       covariates = c("secondary", "TAM", "secTAM"),
#'                       covariateTypes = c("binary", "binary", "binary"),
#'                       covariateLevels = c("sample", "raster", "raster"))
#' PCmod <- writePCModel(PCdat, multiSampsPerSubj = TRUE, typeOfZero = "censored")
#' PCresults <- runPCModel(modelObj = PCmod, PCDataObj = PCdat, slideVar='slide',
#'                         monitorCoefOnly = FALSE,
#'                         nBurnin = 15000, nIter = 40000, nThin = 25)
#' plot(PCresults)

plot.PCResults<-function(object){
  if (class(object) != "PCResults")
    stop("object must of class \"PCResults\".  This object must be created using the runPCModel function.")

  if(is.null(object[['sample']])==TRUE)
    stop('You must monitor the coefficients and keep the sample in order to plot the results.  Rerun the runPCModel function with returnSample=TRUE.')


  coefs<-object[['coefNameMap']]$coefficient
  vars<-object[['coefNameMap']]$variableName
  mcmcsamp<-object[['sample']]

  dimcol<-ceiling(sqrt(length(vars)))

  flag<-1
  iter<-0
  while(flag){
    iter<-iter+1
    flag<-ifelse((dimcol*iter)>=length(vars),0,1)
  }

  par(mfrow=c(iter, dimcol))
  for(i in 1:length(vars)){
    plot(x=seq(1:nrow(mcmcsamp[[1]])), y=mcmcsamp[[1]][,coefs[i]], type='l', ylim=c(min(c(mcmcsamp[[1]][,coefs[i]],mcmcsamp[[2]][,coefs[i]])),max(c(mcmcsamp[[1]][,coefs[i]],mcmcsamp[[2]][,coefs[i]]))), xlab='iteration', ylab='', main=vars[i])
    points(x=seq(1:nrow(mcmcsamp[[1]])), y=mcmcsamp[[2]][,coefs[i]], type='l', col='red')
  }


}


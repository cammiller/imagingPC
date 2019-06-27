

#' Estimate the range parameter
#'
#' The \code{estRange} function is a wrapper for the
#' \code{\link[geoR]{variog}} and \code{\link[geoR]{variofit}} functions in
#' the \code{geoR} package.  \code{estRange} calculates the semivariance at
#' every observed distance within samples in the rescaled data and then fits
#' those estimates to a covariance model.
#'
#' @param rScaleObj An object of class \code{rScaleList}.  This object is the
#' result of using the \code{\link{rScale}} function.
#' @param outcome A character string specifying the outcome of interest.  This
#' is the variable that will later be modeled.
#' @param spatialVar An optional character string specifying a binary or
#' categorical variable.  If a variable is input, the range will be estimated
#' for all levels of that variable.
#' @param semivEst The form of the semivariance estimator.  The options are
#' 'classical' and 'modulus'.  The classical estimator is the method of
#' moments etimator, and the modulus estimator is robust estimator from
#' Hawkins and Cressie.
#' @param logTransform A TRUE/FALSE variable indicating whether or not the
#' outcome should be log-transformed.  Imaging mass spectrometry data are
#' typically lognormally distributed, and so the default is TRUE.
#' @param covarianceModel A character string specifying the form of the
#' covariance model.  The only current option is 'gaussian'.
#'
#' @return A list of class \code{rangeList} containing the estimated range
#' and other information supplied to the \code{estRange} function.
#'
#' \describe{
#'   \item{\code{data}}{A data frame containing the data.}
#'   \item{\code{subjectVar}}{A character string denoting the subject
#'   variable.}
#'   \item{\code{sampleVar}}{A character string denoting the sample variable.}
#'   \item{\code{spatialVar}}{A character string denoting the spatial
#'   variable.}
#'   \item{\code{outcome}}{A character string denoting the outcome of
#'   interest.}
#'   \item{\code{estRange}}{A single value or vector of values representing
#'   the estimated range parameter.  If no spatial variable is given then
#'   this will be a single value.  If a spatial variable is given then this
#'   will be a vector of values, one for each level of the spatial variable.
#'   Each estimated range will be named for its corresponding spatial
#'   variable level.}
#'   \item{\code{estSig2}}{A single value or vector of values representing
#'   the estimated variance parameter \eqn{\sigma^2}.  The parameter
#'   \eqn{\sigma^2} is used to calculate the covariance function.  For more
#'   information, see the \code{\link[geoR]{cov.spatial}} function in the
#'   \code{geoR} package.  If no spatial variable is given then this will be
#'   a single value.  If a spatial variable is given then this will be a
#'   vector of values, one for each level of the spatial variable.  Each
#'   estimated range will be named for its corresponding spatial variable
#'   level.}
#'   \item{\code{semivarFit}}{The empirical variogram.  This is a result of
#'   a call to the \code{\link[geoR]{variog}} function in the \code{geoR}
#'   package.  If no spatial variable is provided, then this is a single
#'   object.  If a spatial variable is provided, then this is a list of
#'   objects, one object for every level of the spatial variable.}
#'   \item{\code{covModelFit}}{The fitted covariance model.  This is a result
#'   of a call to the \code{\link[geoR]{variofit}} function in the \code{geoR}
#'   package.  If no spatial variable is provided, then this is a single
#'   object.  If a spatial variable is provided, then this is a list of
#'   objects, one object for every level of the spatial variable.}
#' }
#'
#' @references Ribeiro, Jr., PJ and Diggle, PJ. 2018. geoR: Analysis of
#' Geostatistical Data. R package version 1.7-2.1.
#' @references Cressie, N and Hawkins, DM. 1980. Robust estimation of the
#' variogram: I. \emph{Journal of the International Association for
#' Mathematical Geology}, 12(2):115-125.
#'
#' @examples
#' data("TAMdata")
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')
#' rangs <- estRange(TAMdata, outcome = 'X1282.auc', spatialVar = 'TAM',
#'                   semivEst = 'modulus', logTransform = TRUE)



estRange <- function(rScaleObj, outcome, spatialVar = NULL, semivEst = "modulus", logTransform = TRUE, covarianceModel = "gaussian") {
    if ((semivEst %in% c("modulus", "classical")) == FALSE)
        stop("The only options for semivest are 'classical' and 'modulus'.")
    #if ((covarianceModel %in% c("gaussian")) == FALSE)
     #   stop("The only covariance model that can currently be fit is the Gaussian covariance model.")
    if ((outcome %in% colnames(rScaleObj[["data"]])) == FALSE)
        stop("The outcome must be a variable name in the dataset.  Make sure you have typed the variable name correctly and check that the variable is in the dataset.")
    if (is.null(spatialVar) == FALSE){
      if((spatialVar %in% colnames(rScaleObj[["data"]])) == FALSE)
        stop("spatialVar must be a variable name in the dataset.  Make sure you have typed the variable name correctly and check that the variable is in the dataset.")
    }
    if (class(rScaleObj) != "rScaleList")
        stop("rScaleObj must of class rScaleList.  This object must be created using the rScale function.")


    data <- rScaleObj[["data"]]
    sampleVar <- rScaleObj[["sampleVar"]]
    subjectVar <- rScaleObj[["subjectVar"]]
    xCoord <- "xPC"
    yCoord <- "yPC"

    if (is.null(spatialVar) == TRUE)
        {

            # data frame of data coordinates ---------------------------------------------
            coords.data <- data.frame(data[[sampleVar]], data[[xCoord]], data[[yCoord]], data[[outcome]])
            colnames(coords.data) <- c("sample", "xCoord", "yCoord", "outcm")

            samplelabs <- unique(coords.data$sample)
            nsamps <- length(samplelabs)
            xdist.test <- rep(NA, nsamps)
            ydist.test <- rep(NA, nsamps)
            coords.datal <- list()
            for (i in 1:nsamps) {
                # center at 0,0 to simplify computation --------------------------------------
                coords.data$xnew[coords.data$sample == samplelabs[i]] <- coords.data$xCoord[coords.data$sample == samplelabs[i]] - (max(coords.data$xCoord[coords.data$sample ==
                  samplelabs[i]]) - min(coords.data$xCoord[coords.data$sample == samplelabs[i]]))/2 - min(coords.data$xCoord[coords.data$sample ==
                  samplelabs[i]])

                coords.data$ynew[coords.data$sample == samplelabs[i]] <- coords.data$yCoord[coords.data$sample == samplelabs[i]] - (max(coords.data$yCoord[coords.data$sample ==
                  samplelabs[i]]) - min(coords.data$yCoord[coords.data$sample == samplelabs[i]]))/2 - min(coords.data$yCoord[coords.data$sample ==
                  samplelabs[i]])

                xdist.test[i] <- abs(range(coords.data$xnew[coords.data$sample == samplelabs[i]])[2])
                ydist.test[i] <- abs(range(coords.data$ynew[coords.data$sample == samplelabs[i]])[2])

                coords.datal[[samplelabs[i]]] <- coords.data[coords.data$sample == samplelabs[i], ]
            }

            xdist <- max(c(xdist.test, ydist.test))
            ydist <- max(c(xdist.test, ydist.test))

            g <- 2 * xdist + 1


            # need to determine the sd.support
            datanozero <- coords.data[coords.data$outcm != 0, ]
            datanozero <- datanozero[order(datanozero$sample), ]
            ROInums <- unique(datanozero$sample)


            q <- ceiling((g - 1) * sqrt(2)) + 1
            datanozero$xnewer <- datanozero$xnew + (datanozero$sample - 1) * (q + (g))

            ## need to determine break points
            datnewl <- list()
            undistl <- list()
            for (j in 1:length(ROInums)) {
                datnewl[[j]] <- data.frame(datanozero$xnew[datanozero$sample == ROInums[j]], datanozero$ynew[datanozero$sample == ROInums[j]])
                colnames(datnewl[[j]]) <- c("xnew", "ynew")

                undistl[[j]] <- unique(dist(datnewl[[j]]))
            }

            undist <- unlist(undistl)
            undist <- sort(undist)
            undist <- unique(undist)
            undist <- round(undist, digits = 5)
            undist <- unique(undist)
            mndist <- min(dist(undist))

            distdf <- data.frame(undist)
            colnames(distdf) <- c("undist")

            distl <- matrix(NA, ncol = 2, nrow = nrow(distdf))
            for (j in 1:nrow(distdf)) {
                distl[j, ] <- c((distdf[j, 1] - mndist/3), (distdf[j, 1] + mndist/3))
            }

            dbreaks <- rep(NA, nrow(distdf * 2))
            for (j in 1:nrow(distdf)) {
                dbreaks[((j - 1) * 2 + 1):(j * 2)] <- distl[j, ]
            }

            if (logTransform == TRUE) {
                vardat <- data.frame(datanozero$xnewer, datanozero$ynew, log(datanozero$outcm))
            }
            if (logTransform != TRUE) {
                vardat <- data.frame(datanozero$xnewer, datanozero$ynew, datanozero$outcm)
            }

            colnames(vardat) <- c("xnew", "ynew", "outcm")

            gdat <- as.geodata(vardat, coords.col = 1:2, data.col = 3)
            capture.output({
                varfit <- tryCatch(variog(geodata = gdat, estimator.type = semivEst, breaks = dbreaks, max.dist = (q - 0.5)), error = function(e) {
                  print(NA)
                })
            })

            capture.output({
                vparam.nofix <- tryCatch(suppressWarnings(variofit(varfit, cov.model = covarianceModel, fix.nugget = FALSE)), error = function(e) {
                  print(NA)
                })
            })

            outstat.range <- vparam.nofix$cov.pars[2]
            outstat.sig2 <- vparam.nofix$cov.pars[1]
            semivarFit <- varfit
            covModelFit <- vparam.nofix

            outlist <- list(data = data, subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, outcome = outcome, estRange = outstat.range,
                            estSig2 = outstat.sig2, semivarFit = semivarFit, covModelFit = covModelFit)
            class(outlist) <- "rangeList"
            return(outlist)

        }  #if spatialVar==NULL



    #### length(spatialVar)==1 ####
    if (is.null(spatialVar) == FALSE) {
        # data frame of data coordinates ---------------------------------------------
        coords.data <- data.frame(data[[sampleVar]], data[[xCoord]], data[[yCoord]], data[[outcome]], data[[spatialVar]])
        colnames(coords.data) <- c("sample", "xCoord", "yCoord", "outcm", "spatialVar")

        samplelabs <- unique(coords.data$sample)
        nsamps <- length(samplelabs)
        xdist.test <- rep(NA, nsamps)
        ydist.test <- rep(NA, nsamps)
        coords.datal <- list()
        for (i in 1:nsamps) {
            # center at 0,0 to simplify computation --------------------------------------
            coords.data$xnew[coords.data$sample == samplelabs[i]] <- coords.data$xCoord[coords.data$sample == samplelabs[i]] - (max(coords.data$xCoord[coords.data$sample ==
                samplelabs[i]]) - min(coords.data$xCoord[coords.data$sample == samplelabs[i]]))/2 - min(coords.data$xCoord[coords.data$sample ==
                samplelabs[i]])

            coords.data$ynew[coords.data$sample == samplelabs[i]] <- coords.data$yCoord[coords.data$sample == samplelabs[i]] - (max(coords.data$yCoord[coords.data$sample ==
                samplelabs[i]]) - min(coords.data$yCoord[coords.data$sample == samplelabs[i]]))/2 - min(coords.data$yCoord[coords.data$sample ==
                samplelabs[i]])

            xdist.test[i] <- abs(range(coords.data$xnew[coords.data$sample == samplelabs[i]])[2])
            ydist.test[i] <- abs(range(coords.data$ynew[coords.data$sample == samplelabs[i]])[2])

            coords.datal[[samplelabs[i]]] <- coords.data[coords.data$sample == samplelabs[i], ]
        }

        xdist <- max(c(xdist.test, ydist.test))
        ydist <- max(c(xdist.test, ydist.test))

        g <- 2 * xdist + 1


        coords.data <- coords.data[order(coords.data$spatialVar), ]
        rpvals <- sort(unique(data[[spatialVar]]))

        coords.data <- coords.data[order(coords.data$sample, coords.data$spatialVar), ]

        outstat.range <- rep(NA, length(rpvals))
        outstat.sig2 <- rep(NA, length(rpvals))
        semivarFit <- list()
        covModelFit <- list()

        for (i in 1:length(rpvals)) {
            # need to determine the sd.support need to determine the sd.support
            datanozero <- coords.data[coords.data$outcm != 0 & coords.data$spatialVar == rpvals[i], ]
            datanozero <- datanozero[order(datanozero$sample), ]
            ROInums <- unique(datanozero$sample)

            q <- ceiling((g - 1) * sqrt(2)) + 1
            datanozero$xnewer <- datanozero$xnew + (datanozero$sample - 1) * (q + (g))

            ## need to determine break points
            datnewl <- list()
            undistl <- list()
            for (j in 1:length(ROInums)) {
                datnewl[[j]] <- data.frame(datanozero$xnew[datanozero$sample == ROInums[j]], datanozero$ynew[datanozero$sample == ROInums[j]])
                colnames(datnewl[[j]]) <- c("xnew", "ynew")

                undistl[[j]] <- unique(dist(datnewl[[j]]))
            }

            undist <- unlist(undistl)
            undist <- sort(undist)
            undist <- unique(undist)
            undist <- round(undist, digits = 5)
            undist <- unique(undist)
            mndist <- min(dist(undist))

            distdf <- data.frame(undist)
            colnames(distdf) <- c("undist")

            distl <- matrix(NA, ncol = 2, nrow = nrow(distdf))
            for (j in 1:nrow(distdf)) {
                distl[j, ] <- c((distdf[j, 1] - mndist/3), (distdf[j, 1] + mndist/3))
            }

            dbreaks <- rep(NA, nrow(distdf * 2))
            for (j in 1:nrow(distdf)) {
                dbreaks[((j - 1) * 2 + 1):(j * 2)] <- distl[j, ]
            }

            if (logTransform == TRUE) {
                vardat <- data.frame(datanozero$xnewer, datanozero$ynew, log(datanozero$outcm))
            }
            if (logTransform != TRUE) {
                vardat <- data.frame(datanozero$xnewer, datanozero$ynew, datanozero$outcm)
            }

            colnames(vardat) <- c("xnew", "ynew", "outcm")

            gdat <- as.geodata(vardat, coords.col = 1:2, data.col = 3)
            capture.output({
                varfit <- tryCatch(variog(geodata = gdat, estimator.type = semivEst, breaks = dbreaks, max.dist = (q - 0.5)), error = function(e) {
                  print(NA)
                })
            })

            capture.output({
                vparam.nofix <- tryCatch(suppressWarnings(variofit(varfit, cov.model = covarianceModel, fix.nugget = FALSE)), error = function(e) {
                  print(NA)
                })
            })

            outstat.range[i] <- vparam.nofix$cov.pars[2]
            outstat.sig2[i] <- vparam.nofix$cov.pars[1]
            semivarFit[[i]] <- varfit
            covModelFit[[i]] <- vparam.nofix

        }  #i in 1:length(rpvals)

        names(outstat.range) <- rpvals
        names(outstat.sig2) <- rpvals
        outlist <- list(data = data, subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, outcome = outcome, estRange = outstat.range,
            estSig2 = outstat.sig2, semivarFit = semivarFit, covModelFit = covModelFit)
        class(outlist) <- "rangeList"
        return(outlist)

    }
}

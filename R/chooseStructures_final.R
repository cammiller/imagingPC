
#' Choose support structure(s)
#'
#' The \code{chooseStructures} function selects support structures for every
#' level of the spatial variable (\code{spatialVar} argument) specified in the
#' \code{\link{estRange}} function.  If no spatial variable is provided, then
#' only a single set of support structures is chosen, one structure for each
#' sample.  See the \strong{Details} section for more information on how the
#' structures are selected.
#'
#' @param rangeObj An object of class \code{rangeList}.  This object is a list
#' that contains the estimated range parameter(s).  This is obtained by a call
#' to the \code{\link{estRange}} function.
#' @param cDist A distance expressing how far beyond the limits of the data
#' that the support sites should be placed.  In the \code{rScale} function,
#' the data are rescaled such that data locations are a distance of 1 apart.
#' Therefore, a cDist of 0.5 is half the distance between data points.
#' @param sdWithin The maximum allowed distance between any data point and
#' its nearest support site (in terms of the number of standard deviations
#' of the smoothing kernel).  The default is 1, meaning that data must be
#' within one standard deviation of the smoothing kernel of a support site.
#' Smaller values of sdWithin lead to demser support structures.  This value
#' should never be set higher than 1.
#' @param defaultStructure The default structure if the number of support
#' points exceeds the proportion threshold of data set by
#' \code{thresholdNumSup}.  The options are '5x' or 'nextHighest'.
#' @param sdDefault The default standard deviation of the smoothing kernel.
#' The default standard deviation is used if (1) the number of support sites
#' exceeds that allowed by \code{thresholdNumSup} or (2) the recommended
#' standard deviation of the smoothing kernel exceeds the largest distance
#' between data points.  \code{sdDefault} is a value mutliplied by the maximum
#' observed distance between any data point and the nearest support site.  For
#' an \code{sdDefault} of 1, the maximum distance is mutliplied by 1, and that
#' value is used for the standard deviation of the smoothing kernel when a
#' default support structure is used.
#' @param thresholdNumSup A threshold that limits the number of support sites
#' allowed for computational efficiency.  The threshold is a proportion of the
#' observations.  A threshold of 0.5 limits the number of allowed support sites
#' to 50\% of the number of observations.
#' @param sdElim A threshold for eliminating support sites (in terms of the
#' standard deviation of the smoothing kernel).  Support sites are removed if
#' they are not within the specified number of standard deviations of a data
#' point.  This is used to remove support sites for irregularly shaped samples.
#' @param noZeroRangeStructs A TRUE/FALSE argument specifying how to handle
#' instances in which the estimated range parameter is 0.  If
#' \code{noZeroRangeStructs=TRUE} and the estimated range parameter is 0, then
#' no support structure will be used for the corresponding data.  If a spatial
#' variable was provided, then no support structure will be used for the data
#' in the corresponding level of that variable.  If no spatial variable was
#' provided, then a support structure will not be used at all.  Instead, a
#' raster-level intercept will be incorporated to account for extra noise.
#' If \code{noZeroRanges=FALSE} and the estimated range parameter is 0, then a
#' default support structure will be used.
#'
#' @return An object of class \code{structureList} that includes the
#' information below.
#' \describe{
#'   \item{\code{structure}}{A list of information about the support
#'   structures.  If a spatial variable is provided, then there will be a list
#'   for each level of that variable.  The structure list includes
#'   \code{recommendStruct} (the recommended support structure),
#'   \code{recommendSd} (the recommended standard deviation of the smoothing
#'   kernel), \code{nSupportReduced} (the number of support sites after
#'   removing those more than \code{sdElim} standard deviations away from a
#'   data point), \code{coordsU} (a data frame of the coordinates of the
#'   support sites, by sample), \code{coordsData} (a data frame of the
#'   coordinates of the data), \code{defaultStatus} (a character string
#'   informing the user if a default support structure was chosen), and
#'   \code{nRowSupport} (the number of rows of support sites for alternating
#'   support structures).}
#'   \item{\code{data}}{The dataset, appended with new x- and y-coordinates.}
#'   \item{\code{subjectVar}}{A character string denoting the subject
#'   variable.}
#'   \item{\code{sampleVar}}{A character string denoting the sample variable.}
#'   \item{\code{spatialVar}}{A character string denoting the spatial
#'   variable.}
#'   \item{\code{outcome}}{A character string denoting the outcome of
#'   interest.}
#' }
#'
#' @details The process convolution approach uses points called support
#' sites that help to account for underlying spatial structure in data.
#' Collections of support sites are called support structures.
#' How many and where to place support sites were questions that drove the
#' research leading to the creation of this function.
#'
#' The \code{chooseStructures} function uses the procedure detailed below.
#' If a spatial variable was provided in the \code{\link{estRange}} function,
#' then the procedure is performed for every level of the spaital variable.
#'
#' The procedure starts by building all of the fixed support structures.  The
#' fixed support structures are a group of five support structures ranging
#' from five to twelve support sites.  The minimum distance between each data
#' point and the closest support site is measured for all observations.
#' Support sites that are more than (\code{sdElim})*(recommended standard
#' deviation of smoothing kernel) are removed.  The total number of support
#' sites is counted across all samples, and each structure is checked to see
#' if all data are within (\code{sdWithin})*(recommended standard deviation
#' of smoothing kernel).
#'
#' Once the fixed support structures are created, an iterative loop is used to
#' generate the alternating support structures.  These are support structures
#' in which the rows of support sites have alternating numbers of support
#' sites.  For each support, the support structure is built and support sites
#' more than (\code{sdElim})*(recommended standard deviation of smoothing
#' kernel) are removed.  The support structure is then checked to see if (1)
#' all the data are within (\code{sdWithin})*(recommended standard deviation
#' of smoothing kernel) and (2) the number of support sites across all samples
#' has exceeded (\code{thresholdNumSup})*(total number of observations).  If
#' either criteria is met, the iteratie loop stops.  The alternating
#' structures are compared to the fixed structures, and the support structure
#' with the fewest support sites in which all of the data are within
#' (sdWithin)*(recommended standard deviation of smoothing kernel) is chosen
#' as the support structure.
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


chooseStructures <- function(rangeObj, cDist = 0.1, sdWithin = 1, defaultStructure = "nextHighest", sdDefault = 1, thresholdNumSup = 0.5, sdElim = 2, noZeroRangeStructs=TRUE) {
    if (class(rangeObj) != "rangeList")
        stop("rangeObj must of class rangeList  This object must be created using the estRange function.")
    if (is.numeric(cDist) == FALSE)
        stop("cDist must be numeric.")
    if (is.numeric(cDist) == TRUE & cDist < 0)
        stop("cDist cannot be negative.")
    if (is.numeric(sdWithin) == FALSE)
        stop("sdWithin must be numeric.")
    if (is.numeric(sdWithin) == TRUE & sdWithin < 0)
        stop("sdWithin cannot be negative.")
    if ((defaultStructure %in% c("5x", "nextHighest")) == FALSE)
        stop("The only options for defaultStructure are \"5x\" and \"nextHighest\".")
    if (is.numeric(sdDefault) == FALSE)
        stop("sdDefault must be numeric.")
    if (is.numeric(sdDefault) == TRUE & sdDefault < 0)
        stop("sdDefault cannot be negative.")
    if (is.numeric(thresholdNumSup) == FALSE)
        stop("thresholdNumSup must be numeric.")
    if (is.numeric(thresholdNumSup) == TRUE & thresholdNumSup < 0)
        stop("thresholdNumSup cannot be negative.")
    if (is.numeric(sdElim) == FALSE)
        stop("sdElim must be numeric.")
    if (is.numeric(sdElim) == TRUE & sdElim < 0)
        stop("sdElim cannot be negative.")

    data <- rangeObj[["data"]]
    sampleVar <- rangeObj[["sampleVar"]]
    subjectVar <- rangeObj[["subjectVar"]]
    spatialVar <- rangeObj[["spatialVar"]]
    outcome <- rangeObj[["outcome"]]
    estrange <- rangeObj[["estRange"]]


    if (is.null(spatialVar) == TRUE) {
      if(estrange==0 & noZeroRangeStructs==TRUE){
        coords.data <- data.frame(data[[sampleVar]], data$xPC, data$yPC)
        colnames(coords.data) <- c("sample", "xPC", "yPC")

        outp <- list()
        outp[["structure"]] <- list(recommendStruct='noStructure', coordsData = coords.data)
        outp[["data"]] <- data
        outp[["subjectVar"]] <- subjectVar
        outp[["sampleVar"]] <- sampleVar
        outp[["spatialVar"]] <- spatialVar
        outp[["outcome"]] <- outcome
      } else {
        nobstot <- nrow(data)
        tempestrange <- 0.5 * sqrt(2)

        # calculate the threshold for the number of support points allowed -----------
        thresh.nsup <- thresholdNumSup * nobstot

        # distance function ----------------------------------------------------------
        distance <- function(x1, x2, y1, y2) {
          sqrt((x1 - x2)^2 + (y1 - y2)^2)
        }

        # data frame of data coordinates ---------------------------------------------
        coords.data <- data.frame(data[[sampleVar]], data$xPC, data$yPC)
        colnames(coords.data) <- c("sample", "xPC", "yPC")

        samplelabs <- unique(coords.data$sample)
        nsamps <- length(samplelabs)
        xdist.test <- rep(NA, nsamps)
        ydist.test <- rep(NA, nsamps)
        nobs <- matrix(nrow = nsamps, ncol = 2)
        nobs <- as.data.frame(nobs)
        colnames(nobs) <- c("sample", "nobs")
        nobs$sample <- samplelabs
        coords.datal <- list()
        for (i in 1:nsamps) {

          xdist.test[i] <- abs(range(coords.data$xPC[coords.data$sample == samplelabs[i]])[2])
          ydist.test[i] <- abs(range(coords.data$yPC[coords.data$sample == samplelabs[i]])[2])

          nobs[i, 2] <- nrow(coords.data[coords.data$sample == samplelabs[i], ])

          coords.datal[[samplelabs[i]]] <- coords.data[coords.data$sample == samplelabs[i], ]
        }

        xdist <- max(c(xdist.test, ydist.test))
        ydist <- max(c(xdist.test, ydist.test))

        # cycle through the support structures and for each determine the maximum distance between the data and nearest support point

        ##### Fixed support structures #####
        maxdistm1 <- matrix(NA, nrow = 5 * nsamps, ncol = 9)
        maxdistm1 <- as.data.frame(maxdistm1)
        colnames(maxdistm1) <- c("sample", "mstructure", "nrow.sup", "npoints", "npoints.reduced", "maxdistm", "meansupdist", "recsd", "valid")
        maxdistm1$sample <- rep(samplelabs, each = 5)
        maxdistm1$mstructure <- rep(c("5x", "6x", "9x", "10x", "12x"), nsamps)

        distthresh <- max(dist(cbind(coords.data$xPC, coords.data$yPC)))

        coords.u <- list()

        for (i in 1:nrow(maxdistm1)) {
          # x and y coordinates for support
          if (maxdistm1$mstructure[i] == "4x") {
            m <- 4

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- xdist + cDist/sqrt(2)
            m2y <- ydist + cDist/sqrt(2)
            m3x <- -(xdist + cDist/sqrt(2))
            m3y <- -(ydist + cDist/sqrt(2))
            m4x <- xdist + cDist/sqrt(2)
            m4y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x)
            yw <- c(m1y, m2y, m3y, m4y)
          }

          if (maxdistm1$mstructure[i] == "5x") {
            m <- 5

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- xdist + cDist/sqrt(2)
            m2y <- ydist + cDist/sqrt(2)
            m3x <- 0
            m3y <- 0
            m4x <- -(xdist + cDist/sqrt(2))
            m4y <- -(ydist + cDist/sqrt(2))
            m5x <- xdist + cDist/sqrt(2)
            m5y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x, m5x)
            yw <- c(m1y, m2y, m3y, m4y, m5y)

          }

          if (maxdistm1$mstructure[i] == "6x") {
            m <- 6

            quada <- 3/4
            quadb <- (ydist + (cDist/sqrt(2)))
            quadc <- -(xdist^2 + ((2 * xdist * cDist)/sqrt(2)) + (cDist^2) + (ydist^2) + ((2 * ydist * cDist)/(sqrt(2))))

            distx <- (-quadb + sqrt(quadb^2 - 4 * quada * quadc))/(2 * quada)

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- xdist + cDist/sqrt(2)
            m2y <- ydist + cDist/sqrt(2)
            m3x <- 0
            m3y <- distx/2
            m4x <- 0
            m4y <- -distx/2
            m5x <- -(xdist + cDist/sqrt(2))
            m5y <- -(ydist + cDist/sqrt(2))
            m6x <- xdist + cDist/sqrt(2)
            m6y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x, m5x, m6x)
            yw <- c(m1y, m2y, m3y, m4y, m5y, m6y)
          }

          if (maxdistm1$mstructure[i] == "9x") {
            m <- 9

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- xdist + cDist/sqrt(2)
            m2y <- ydist + cDist/sqrt(2)
            m3x <- 0
            m3y <- (ydist/1.5)
            m4x <- -(xdist/1.5)
            m4y <- 0
            m5x <- 0
            m5y <- 0
            m6x <- (xdist/1.5)
            m6y <- 0
            m7x <- 0
            m7y <- -(ydist/1.5)
            m8x <- -(xdist + cDist/sqrt(2))
            m8y <- -(ydist + cDist/sqrt(2))
            m9x <- xdist + cDist/sqrt(2)
            m9y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x)
            yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y)
          }

          if (maxdistm1$mstructure[i] == "10x") {
            m <- 10

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- xdist + cDist/sqrt(2)
            m2y <- ydist + cDist/sqrt(2)
            m3x <- 0
            m3y <- ydist/1.25
            m4x <- -(xdist/4)
            m4y <- (ydist/4)
            m5x <- -(xdist/1.25)
            m5y <- 0
            m6x <- (xdist/1.25)
            m6y <- 0
            m7x <- (xdist/4)
            m7y <- -(ydist/4)
            m8x <- 0
            m8y <- -(ydist/1.25)
            m9x <- -(xdist + cDist/sqrt(2))
            m9y <- -(ydist + cDist/sqrt(2))
            m10x <- xdist + cDist/sqrt(2)
            m10y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x)
            yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y, m10y)
          }

          if (maxdistm1$mstructure[i] == "12x") {
            m <- 12

            m1x <- -(xdist + cDist/sqrt(2))
            m1y <- ydist + cDist/sqrt(2)
            m2x <- 0
            m2y <- ydist
            m3x <- xdist + cDist/sqrt(2)
            m3y <- ydist + cDist/sqrt(2)
            m4x <- -(xdist/3)
            m4y <- (ydist/3)
            m5x <- (xdist/3)
            m5y <- (ydist/3)
            m6x <- -(xdist)
            m6y <- 0
            m7x <- (xdist)
            m7y <- 0
            m8x <- -(xdist/3)
            m8y <- -(ydist/3)
            m9x <- (xdist/3)
            m9y <- -(ydist/3)
            m10x <- -(xdist + cDist/sqrt(2))
            m10y <- -(ydist + cDist/sqrt(2))
            m11x <- 0
            m11y <- -(ydist)
            m12x <- xdist + cDist/sqrt(2)
            m12y <- -(ydist + cDist/sqrt(2))

            xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x, m11x, m12x)
            yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y, m10y, m11y, m12y)
          }

          # Matrix of coordinates --------------------------------------------------
          coords.w <- data.frame(xw, yw)
          colnames(coords.w) <- c("xw", "yw")

          # calculate distance between data points and each support point ----------
          distm <- matrix(NA, nrow = nobs[nobs$sample == maxdistm1$sample[i], 2], ncol = m)
          # mindist <- matrix(NA, nrow = nobs[nobs$sample==maxdistm1$sample[i],2], ncol = 1)
          for (j in 1:nrow(distm)) {
            for (k in 1:m) {
              distm[j, k] <- distance(coords.datal[[maxdistm1$sample[i]]]$xPC[j], coords.w$xw[k], coords.datal[[maxdistm1$sample[i]]]$yPC[j],
                                      coords.w$yw[k])
            }
            # mindist[j, 1] <- min(distm[j, ])
          }

          mindist <- apply(distm, 1, min)
          mindistsup <- apply(distm, 2, min)
          if (estrange < tempestrange) {
            sup.reduced <- which(mindistsup <= ((tempestrange/sqrt(2)) * sdElim))
          }
          if (estrange >= tempestrange) {
            sup.reduced <- which(mindistsup <= ((estrange/sqrt(2)) * sdElim))
          }
          nsup.reduced <- length(sup.reduced)

          coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]] <- coords.w[sup.reduced, ]
          coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]]$sample <- maxdistm1$sample[i]
          colnames(coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]]) <- c("x.omega", "y.omega",
                                                                                                                             "sample")
          coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]] <- coords.u[[paste(maxdistm1$mstructure[i],
                                                                                                                                 "_sample", maxdistm1$sample[i], sep = "", collapse = "")]][, c(3, 1, 2)]


          maxdistm1[i, 6] <- max(mindist)


          # need to calculate the average distance from each support point to the -- nearest support point
          # --------------------------------------------------
          supdistmat <- matrix(NA, nrow = m, ncol = m)
          for (j in 1:m) {
            for (k in 1:m) {
              supdistmat[j, k] <- distance(xw[j], xw[k], yw[j], yw[k])
            }
          }
          supdistmat <- t(apply(supdistmat, 1, sort))
          supdistmat <- supdistmat[, -1]
          maxdistm1$meansupdist[i] <- mean(supdistmat[, 1])
          maxdistm1$npoints[i] <- as.numeric(gsub(pattern = "x", replacement = "", x = maxdistm1$mstructure[i]))
          maxdistm1$npoints.reduced[i] <- nsup.reduced

        }


        ##### Alternating Support Points #####

        maxdistm2 <- matrix(NA, nrow = 1, ncol = 9)
        maxdistm2 <- as.data.frame(maxdistm2)
        colnames(maxdistm2) <- c("sample", "mstructure", "nrow.sup", "npoints", "npoints.reduced", "maxdistm", "meansupdist", "recsd", "valid")

        flag <- 1
        nrow.sup <- 3
        validstruct <- rep(NA, nsamps)
        while (flag) {

          nsup.reduced <- rep(NA, nsamps)
          for (z in 1:nsamps) {
            if (((nrow.sup/2)%%1 == 0) == FALSE) {
              # odd
              nrow.large <- (nrow.sup + 1)/2
              nrow.small <- (nrow.sup - 1)/2
              npoints.large <- nrow.sup
              npoints.small <- nrow.sup - 1

              nsup.tot <- nrow.large * npoints.large + nrow.small * npoints.small
            }
            if (((nrow.sup/2)%%1 == 0) == TRUE) {
              # even
              nrow.large <- (nrow.sup)/2
              nrow.small <- nrow.large
              npoints.large <- nrow.sup
              npoints.small <- nrow.sup - 1

              nsup.tot <- nrow.large * npoints.large + nrow.small * npoints.small
            }

            xmaxd <- xdist + cDist/sqrt(2)
            ymaxd <- ydist + cDist/sqrt(2)

            xdisttot <- xmaxd * 2
            ydisttot <- ymaxd * 2

            distbetweenx <- xdisttot/(nrow.sup - 1)
            distbetweeny <- ydisttot/(nrow.sup - 1)


            # create a vector of number of points for each row -----------------------
            npoints <- rep(NA, nrow.sup)
            for (i in 1:length(npoints)) {
              if (((i/2)%%1 == 0) == FALSE) {
                npoints[i] <- npoints.large
              }
              if (((i/2)%%1 == 0) == TRUE) {
                npoints[i] <- npoints.small
              }
            }
            csnpoints <- c(0, cumsum(npoints))

            # create vector of y coordinates
            yw <- rep(NA, nsup.tot)
            for (i in 1:nrow.sup) {
              yw[(csnpoints[i] + 1):csnpoints[(i + 1)]] <- ymaxd - (i - 1) * distbetweeny
            }

            # create vector of x coordinates alternating rows are identical, --------- so just need to calculate the two rows
            # ---------------------------------
            xw.large <- rep(NA, npoints.large)
            for (i in 1:length(xw.large)) {
              xw.large[i] <- -xmaxd + distbetweenx * (i - 1)
            }

            start.small <- -xmaxd + 0.5 * distbetweenx
            xw.small <- rep(NA, npoints.small)
            for (i in 1:length(xw.small)) {
              xw.small[i] <- start.small + distbetweenx * (i - 1)
            }

            xw <- rep(NA, 1)
            for (i in 1:length(npoints)) {
              if (((i/2)%%1 == 0) == FALSE) {
                xw <- append(xw, xw.large)
              }
              if (((i/2)%%1 == 0) == TRUE) {
                xw <- append(xw, xw.small)
              }
            }
            xw <- xw[-1]

            coords.w <- cbind(xw, yw)
            coords.w <- as.data.frame(coords.w)

            # calculate distance between data points and each support point ----------
            distm <- matrix(NA, nrow = nobs[nobs$sample == samplelabs[z], 2], ncol = nsup.tot)
            for (j in 1:nobs[nobs$sample == samplelabs[z], 2]) {
              for (k in 1:nsup.tot) {
                distm[j, k] <- distance(coords.datal[[samplelabs[z]]]$xPC[j], coords.w$xw[k], coords.datal[[samplelabs[z]]]$yPC[j], coords.w$yw[k])
              }
            }

            mindist <- apply(distm, 1, min)
            mindistsup <- apply(distm, 2, min)
            if (estrange < tempestrange) {
              sup.reduced <- which(mindistsup <= ((tempestrange/sqrt(2)) * sdElim))
            }
            if (estrange >= tempestrange) {
              sup.reduced <- which(mindistsup <= ((estrange/sqrt(2)) * sdElim))
            }
            nsup.reduced[z] <- length(sup.reduced)

            coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]] <- coords.w[sup.reduced, ]
            coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]]$sample <- samplelabs[z]
            colnames(coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]]) <- c("x.omega", "y.omega", "sample")
            coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]] <- coords.u[[paste(nsup.tot, "alternate_sample",
                                                                                                                       samplelabs[z], sep = "", collapse = "")]][, c(3, 1, 2)]

            # need to calculate the average distance from each support point to the -- nearest support point
            # --------------------------------------------------
            supdistmat <- matrix(NA, nrow = nsup.tot, ncol = nsup.tot)
            for (j in 1:nsup.tot) {
              for (k in 1:nsup.tot) {
                supdistmat[j, k] <- distance(xw[j], xw[k], yw[j], yw[k])
              }
            }
            supdistmat <- t(apply(supdistmat, 1, sort))
            supdistmat <- supdistmat[, -1]

            maxdistm2vec <- rep(NA, 9)
            maxdistm2vec[1] <- samplelabs[z]
            maxdistm2vec[2] <- paste(nsup.tot, "alternate", sep = "", collapse = "")
            maxdistm2vec[3] <- nrow.sup
            maxdistm2vec[4] <- nsup.tot
            maxdistm2vec[5] <- nsup.reduced[z]
            maxdistm2vec[6] <- max(mindist)
            maxdistm2vec[7] <- mean(supdistmat[, 1])


            maxdistm2 <- rbind(maxdistm2, maxdistm2vec)

            validstruct[z] <- ifelse(as.numeric(maxdistm2vec[6]) <= (estrange/sqrt(2)), 1, 0)
          }  #z (sample)
          totnsup.reduced <- sum(nsup.reduced)


          flag <- ifelse((totnsup.reduced > thresh.nsup) | (sum(validstruct) == length(validstruct)), 0, 1)
          nrow.sup <- nrow.sup + 1
        }  # while
        maxdistm2 <- maxdistm2[-1, ]
        maxdistm2$nrow.sup <- as.numeric(maxdistm2$nrow.sup)
        maxdistm2$npoints <- as.numeric(maxdistm2$npoints)
        maxdistm2$npoints.reduced <- as.numeric(maxdistm2$npoints.reduced)
        maxdistm2$maxdistm <- as.numeric(maxdistm2$maxdistm)
        maxdistm2$meansupdist <- as.numeric(maxdistm2$meansupdist)


        maxdistm <- rbind(maxdistm1, maxdistm2)
        maxdistm$nrow.sup <- as.numeric(maxdistm$nrow.sup)
        maxdistm$npoints <- as.numeric(maxdistm$npoints)
        maxdistm$npoints.reduced <- as.numeric(maxdistm$npoints.reduced)
        maxdistm$maxdistm <- as.numeric(maxdistm$maxdistm)
        maxdistm$meansupdist <- as.numeric(maxdistm$meansupdist)

        mstructs <- matrix(nrow = length(unique(maxdistm$mstructure)), ncol = 2)
        mstructs <- as.data.frame(mstructs)
        colnames(mstructs) <- c("mstruct", "totnsup.reduced")
        mstructs$mstruct <- unique(maxdistm$mstructure)
        for (i in 1:nrow(mstructs)) {
          mstructs[i, 2] <- sum(maxdistm$npoints.reduced[maxdistm$mstructure == mstructs[i, 1]])
        }
        mstructs <- mstructs[mstructs$totnsup.reduced <= thresh.nsup, ]

        maxdistm <- maxdistm[maxdistm$mstructure %in% mstructs$mstruct, ]
        # order by sum of support points
        mstructs <- mstructs[order(mstructs$totnsup.reduced), ]
        mstructs.ord <- mstructs$mstruct
        maxdistm$mstructure <- factor(maxdistm$mstructure, levels = mstructs.ord)

        maxdistm <- maxdistm[order(maxdistm$mstructure), ]

        # recommended SD
        maxdistm$recsd <- estrange/sqrt(2)
        for (i in 1:nrow(maxdistm)) {
          maxdistm$valid[i] <- ifelse(maxdistm$maxdistm[i] <= maxdistm$recsd[i] * sdWithin, 1, 0)
        }

        mstructs$nvalid <- NA
        for (i in 1:nrow(mstructs)) {
          mstructs$nvalid[i] <- nrow(maxdistm[maxdistm$mstructure == mstructs$mstruct[i] & maxdistm$valid == 1, ])
        }
        validms <- mstructs$mstruct[mstructs$nvalid == nsamps]

        coordnames <- names(coords.u)


        if (estrange > distthresh) {
          suggest <- maxdistm[maxdistm$mstructure == "5x", ]
          recommendSd <- max(suggest$maxdistm) * sdDefault
          recommendStruct <- "5x"
          defaultStatus <- "found valid structure, SD exceeded threshold"
          nRowSupport <- NA
          coordsU <- do.call(rbind.data.frame, coords.u[which(grepl("5x_sample", coordnames))])
        }

        if (estrange <= distthresh) {
          if (length(validms) == 0) {
            defaultStatus <- "no valid structure"
            nRowSupport <- NA

            if (defaultStructure == "5x") {
              suggest <- maxdistm[maxdistm$mstructure == "5x", ]
              recommendSd <- max(suggest$maxdistm) * sdDefault
              recommendStruct <- "5x"
              coordsU <- do.call(rbind.data.frame, coords.u[which(grepl("5x_sample", coordnames))])
            }
            if (defaultStructure == "nextHighest") {
              suggest <- maxdistm[maxdistm$mstructure == mstructs.ord[length(mstructs.ord)], ]
              # suggest <- suggest[1, ]
              if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == TRUE) {
                nRowSupport <- suggest$nrow.sup[1]
              }
              recommendSd <- max(suggest$maxdistm) * sdDefault
              recommendStruct <- mstructs.ord[length(mstructs.ord)]
              coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(mstructs.ord[length(mstructs.ord)], coordnames))])
            }
          }


          if (length(validms) > 0) {
            suggest <- maxdistm[maxdistm$mstructure == validms[1], ]
            if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == TRUE) {
              nRowSupport <- suggest$nrow.sup[1]
              recommendSd <- suggest$recsd[1]
              recommendStruct <- suggest$mstructure[1]
              defaultStatus <- "found valid structure"
              coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(suggest$mstructure[1], coordnames))])
            }
            if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == FALSE) {
              nRowSupport <- NA
              recommendSd <- suggest$recsd[1]
              recommendStruct <- suggest$mstructure[1]
              defaultStatus <- "found valid structure"
              coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(suggest$mstructure[1], coordnames))])
            }
          }
        }
        row.names(coordsU) <- seq(1:nrow(coordsU))

        outp <- list()
        outp[["structure"]] <- list(recommendStruct = recommendStruct, recommendSd = recommendSd, nSupportReduced = nrow(coordsU), coordsU = coordsU,
                                    coordsData = coords.data, defaultStatus = defaultStatus, nRowSupport = nRowSupport)
        outp[["data"]] <- data
        outp[["subjectVar"]] <- subjectVar
        outp[["sampleVar"]] <- sampleVar
        outp[["spatialVar"]] <- spatialVar
        outp[["outcome"]] <- outcome

      }
    }



    # if separate spatial processes are considered
    if (is.null(spatialVar) == FALSE){
            data <- data[order(data[[spatialVar]]), ]
            spatvarlevels <- unique(data[[spatialVar]])
            nvarlevels <- length(spatvarlevels)
            data <- data[order(data[[sampleVar]], data[[spatialVar]]), ]


            outp <- list()
            for (w in 1:nvarlevels) {
              if(estrange[w]==0 & noZeroRangeStructs==TRUE){
                dataslim <- data[data[[spatialVar]] == spatvarlevels[w], ]

                # data frame of data coordinates ---------------------------------------------
                coords.data <- data.frame(dataslim[[sampleVar]], dataslim$xPC, dataslim$yPC)
                colnames(coords.data) <- c("sample", "xPC", "yPC")

                outp[[paste0("structure_", w)]] <- list(recommendStruct='noStructure', coordsData = coords.data)


              } else {
                dataslim <- data[data[[spatialVar]] == spatvarlevels[w], ]

                nobstot <- nrow(dataslim)
                tempestrange <- 0.5 * sqrt(2)

                # calculate the threshold for the number of support points allowed -----------
                thresh.nsup <- thresholdNumSup * nobstot

                # distance function ----------------------------------------------------------
                distance <- function(x1, x2, y1, y2) {
                  sqrt((x1 - x2)^2 + (y1 - y2)^2)
                }

                # data frame of data coordinates ---------------------------------------------
                coords.data <- data.frame(dataslim[[sampleVar]], dataslim$xPC, dataslim$yPC)
                colnames(coords.data) <- c("sample", "xPC", "yPC")

                samplelabs <- unique(coords.data$sample)
                nsamps <- length(samplelabs)
                xdist.test <- rep(NA, nsamps)
                ydist.test <- rep(NA, nsamps)
                nobs <- matrix(nrow = nsamps, ncol = 2)
                nobs <- as.data.frame(nobs)
                colnames(nobs) <- c("sample", "nobs")
                nobs$sample <- samplelabs
                coords.datal <- list()
                for (i in 1:nsamps) {

                  xdist.test[i] <- abs(range(coords.data$xPC[coords.data$sample == samplelabs[i]])[2])
                  ydist.test[i] <- abs(range(coords.data$yPC[coords.data$sample == samplelabs[i]])[2])

                  nobs[i, 2] <- nrow(coords.data[coords.data$sample == samplelabs[i], ])

                  coords.datal[[samplelabs[i]]] <- coords.data[coords.data$sample == samplelabs[i], ]
                }

                xdist <- max(c(xdist.test, ydist.test))
                ydist <- max(c(xdist.test, ydist.test))

                # cycle through the support structures and for each determine the maximum ---- distance between the data and nearest support point
                # ------------------------

                ##### Fixed support structures #####
                maxdistm1 <- matrix(NA, nrow = 5 * nsamps, ncol = 9)
                maxdistm1 <- as.data.frame(maxdistm1)
                colnames(maxdistm1) <- c("sample", "mstructure", "nrow.sup", "npoints", "npoints.reduced", "maxdistm", "meansupdist", "recsd", "valid")
                maxdistm1$sample <- rep(samplelabs, each = 5)
                maxdistm1$mstructure <- rep(c("5x", "6x", "9x", "10x", "12x"), nsamps)

                distthresh <- max(dist(cbind(coords.data$xPC, coords.data$yPC)))

                coords.u <- list()

                for (i in 1:nrow(maxdistm1)) {
                  # x and y coordinates for support
                  if (maxdistm1$mstructure[i] == "4x") {
                    m <- 4

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- xdist + cDist/sqrt(2)
                    m2y <- ydist + cDist/sqrt(2)
                    m3x <- -(xdist + cDist/sqrt(2))
                    m3y <- -(ydist + cDist/sqrt(2))
                    m4x <- xdist + cDist/sqrt(2)
                    m4y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x)
                    yw <- c(m1y, m2y, m3y, m4y)
                  }

                  if (maxdistm1$mstructure[i] == "5x") {
                    m <- 5

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- xdist + cDist/sqrt(2)
                    m2y <- ydist + cDist/sqrt(2)
                    m3x <- 0
                    m3y <- 0
                    m4x <- -(xdist + cDist/sqrt(2))
                    m4y <- -(ydist + cDist/sqrt(2))
                    m5x <- xdist + cDist/sqrt(2)
                    m5y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x, m5x)
                    yw <- c(m1y, m2y, m3y, m4y, m5y)

                  }

                  if (maxdistm1$mstructure[i] == "6x") {
                    m <- 6

                    quada <- 3/4
                    quadb <- (ydist + (cDist/sqrt(2)))
                    quadc <- -(xdist^2 + ((2 * xdist * cDist)/sqrt(2)) + (cDist^2) + (ydist^2) + ((2 * ydist * cDist)/(sqrt(2))))

                    distx <- (-quadb + sqrt(quadb^2 - 4 * quada * quadc))/(2 * quada)

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- xdist + cDist/sqrt(2)
                    m2y <- ydist + cDist/sqrt(2)
                    m3x <- 0
                    m3y <- distx/2
                    m4x <- 0
                    m4y <- -distx/2
                    m5x <- -(xdist + cDist/sqrt(2))
                    m5y <- -(ydist + cDist/sqrt(2))
                    m6x <- xdist + cDist/sqrt(2)
                    m6y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x, m5x, m6x)
                    yw <- c(m1y, m2y, m3y, m4y, m5y, m6y)
                  }

                  if (maxdistm1$mstructure[i] == "9x") {
                    m <- 9

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- xdist + cDist/sqrt(2)
                    m2y <- ydist + cDist/sqrt(2)
                    m3x <- 0
                    m3y <- (ydist/1.5)
                    m4x <- -(xdist/1.5)
                    m4y <- 0
                    m5x <- 0
                    m5y <- 0
                    m6x <- (xdist/1.5)
                    m6y <- 0
                    m7x <- 0
                    m7y <- -(ydist/1.5)
                    m8x <- -(xdist + cDist/sqrt(2))
                    m8y <- -(ydist + cDist/sqrt(2))
                    m9x <- xdist + cDist/sqrt(2)
                    m9y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x)
                    yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y)
                  }

                  if (maxdistm1$mstructure[i] == "10x") {
                    m <- 10

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- xdist + cDist/sqrt(2)
                    m2y <- ydist + cDist/sqrt(2)
                    m3x <- 0
                    m3y <- ydist/1.25
                    m4x <- -(xdist/4)
                    m4y <- (ydist/4)
                    m5x <- -(xdist/1.25)
                    m5y <- 0
                    m6x <- (xdist/1.25)
                    m6y <- 0
                    m7x <- (xdist/4)
                    m7y <- -(ydist/4)
                    m8x <- 0
                    m8y <- -(ydist/1.25)
                    m9x <- -(xdist + cDist/sqrt(2))
                    m9y <- -(ydist + cDist/sqrt(2))
                    m10x <- xdist + cDist/sqrt(2)
                    m10y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x)
                    yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y, m10y)
                  }

                  if (maxdistm1$mstructure[i] == "12x") {
                    m <- 12

                    m1x <- -(xdist + cDist/sqrt(2))
                    m1y <- ydist + cDist/sqrt(2)
                    m2x <- 0
                    m2y <- ydist
                    m3x <- xdist + cDist/sqrt(2)
                    m3y <- ydist + cDist/sqrt(2)
                    m4x <- -(xdist/3)
                    m4y <- (ydist/3)
                    m5x <- (xdist/3)
                    m5y <- (ydist/3)
                    m6x <- -(xdist)
                    m6y <- 0
                    m7x <- (xdist)
                    m7y <- 0
                    m8x <- -(xdist/3)
                    m8y <- -(ydist/3)
                    m9x <- (xdist/3)
                    m9y <- -(ydist/3)
                    m10x <- -(xdist + cDist/sqrt(2))
                    m10y <- -(ydist + cDist/sqrt(2))
                    m11x <- 0
                    m11y <- -(ydist)
                    m12x <- xdist + cDist/sqrt(2)
                    m12y <- -(ydist + cDist/sqrt(2))

                    xw <- c(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x, m11x, m12x)
                    yw <- c(m1y, m2y, m3y, m4y, m5y, m6y, m7y, m8y, m9y, m10y, m11y, m12y)
                  }

                  # Matrix of coordinates --------------------------------------------------
                  coords.w <- data.frame(xw, yw)
                  colnames(coords.w) <- c("xw", "yw")

                  # calculate distance between data points and each support point ----------
                  distm <- matrix(NA, nrow = nobs[nobs$sample == maxdistm1$sample[i], 2], ncol = m)
                  # mindist <- matrix(NA, nrow = nobs[nobs$sample==maxdistm1$sample[i],2], ncol = 1)
                  for (j in 1:nrow(distm)) {
                    for (k in 1:m) {
                      distm[j, k] <- distance(coords.datal[[maxdistm1$sample[i]]]$xPC[j], coords.w$xw[k], coords.datal[[maxdistm1$sample[i]]]$yPC[j],
                                              coords.w$yw[k])
                    }
                    # mindist[j, 1] <- min(distm[j, ])
                  }

                  mindist <- apply(distm, 1, min)
                  mindistsup <- apply(distm, 2, min)
                  if (estrange[w] < tempestrange) {
                    sup.reduced <- which(mindistsup <= ((tempestrange/sqrt(2)) * sdElim))
                  }
                  if (estrange[w] >= tempestrange) {
                    sup.reduced <- which(mindistsup <= ((estrange[w]/sqrt(2)) * sdElim))
                  }
                  nsup.reduced <- length(sup.reduced)

                  if (nsup.reduced == 0) {
                    coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]] <- coords.w[sup.reduced,
                                                                                                                                    ]
                  }

                  if (nsup.reduced > 0) {
                    coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]] <- coords.w[sup.reduced,
                                                                                                                                    ]
                    coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]]$sample <- maxdistm1$sample[i]
                    colnames(coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]]) <- c("x.omega",
                                                                                                                                       "y.omega", "sample")
                    coords.u[[paste(maxdistm1$mstructure[i], "_sample", maxdistm1$sample[i], sep = "", collapse = "")]] <- coords.u[[paste(maxdistm1$mstructure[i],
                                                                                                                                           "_sample", maxdistm1$sample[i], sep = "", collapse = "")]][, c(3, 1, 2)]
                  }

                  maxdistm1[i, 6] <- max(mindist)


                  # need to calculate the average distance from each support point to the -- nearest support point
                  # --------------------------------------------------
                  supdistmat <- matrix(NA, nrow = m, ncol = m)
                  for (j in 1:m) {
                    for (k in 1:m) {
                      supdistmat[j, k] <- distance(xw[j], xw[k], yw[j], yw[k])
                    }
                  }
                  supdistmat <- t(apply(supdistmat, 1, sort))
                  supdistmat <- supdistmat[, -1]
                  maxdistm1$meansupdist[i] <- mean(supdistmat[, 1])
                  maxdistm1$npoints[i] <- as.numeric(gsub(pattern = "x", replacement = "", x = maxdistm1$mstructure[i]))
                  maxdistm1$npoints.reduced[i] <- nsup.reduced

                }


                ##### Alternating Support Points #####

                maxdistm2 <- matrix(NA, nrow = 1, ncol = 9)
                maxdistm2 <- as.data.frame(maxdistm2)
                colnames(maxdistm2) <- c("sample", "mstructure", "nrow.sup", "npoints", "npoints.reduced", "maxdistm", "meansupdist", "recsd", "valid")

                flag <- 1
                nrow.sup <- 3
                validstruct <- rep(NA, nsamps)
                while (flag) {

                  nsup.reduced <- rep(NA, nsamps)
                  for (z in 1:nsamps) {
                    if (((nrow.sup/2)%%1 == 0) == FALSE) {
                      # odd
                      nrow.large <- (nrow.sup + 1)/2
                      nrow.small <- (nrow.sup - 1)/2
                      npoints.large <- nrow.sup
                      npoints.small <- nrow.sup - 1

                      nsup.tot <- nrow.large * npoints.large + nrow.small * npoints.small
                    }
                    if (((nrow.sup/2)%%1 == 0) == TRUE) {
                      # even
                      nrow.large <- (nrow.sup)/2
                      nrow.small <- nrow.large
                      npoints.large <- nrow.sup
                      npoints.small <- nrow.sup - 1

                      nsup.tot <- nrow.large * npoints.large + nrow.small * npoints.small
                    }

                    xmaxd <- xdist + cDist/sqrt(2)
                    ymaxd <- ydist + cDist/sqrt(2)

                    xdisttot <- xmaxd * 2
                    ydisttot <- ymaxd * 2

                    distbetweenx <- xdisttot/(nrow.sup - 1)
                    distbetweeny <- ydisttot/(nrow.sup - 1)


                    # create a vector of number of points for each row -----------------------
                    npoints <- rep(NA, nrow.sup)
                    for (i in 1:length(npoints)) {
                      if (((i/2)%%1 == 0) == FALSE) {
                        npoints[i] <- npoints.large
                      }
                      if (((i/2)%%1 == 0) == TRUE) {
                        npoints[i] <- npoints.small
                      }
                    }
                    csnpoints <- c(0, cumsum(npoints))

                    # create vector of y coordinates
                    yw <- rep(NA, nsup.tot)
                    for (i in 1:nrow.sup) {
                      yw[(csnpoints[i] + 1):csnpoints[(i + 1)]] <- ymaxd - (i - 1) * distbetweeny
                    }

                    # create vector of x coordinates alternating rows are identical, --------- so just need to calculate the two rows
                    # ---------------------------------
                    xw.large <- rep(NA, npoints.large)
                    for (i in 1:length(xw.large)) {
                      xw.large[i] <- -xmaxd + distbetweenx * (i - 1)
                    }

                    start.small <- -xmaxd + 0.5 * distbetweenx
                    xw.small <- rep(NA, npoints.small)
                    for (i in 1:length(xw.small)) {
                      xw.small[i] <- start.small + distbetweenx * (i - 1)
                    }

                    xw <- rep(NA, 1)
                    for (i in 1:length(npoints)) {
                      if (((i/2)%%1 == 0) == FALSE) {
                        xw <- append(xw, xw.large)
                      }
                      if (((i/2)%%1 == 0) == TRUE) {
                        xw <- append(xw, xw.small)
                      }
                    }
                    xw <- xw[-1]

                    coords.w <- cbind(xw, yw)
                    coords.w <- as.data.frame(coords.w)

                    # calculate distance between data points and each support point ----------
                    distm <- matrix(NA, nrow = nobs[nobs$sample == samplelabs[z], 2], ncol = nsup.tot)
                    for (j in 1:nobs[nobs$sample == samplelabs[z], 2]) {
                      for (k in 1:nsup.tot) {
                        distm[j, k] <- distance(coords.datal[[samplelabs[z]]]$xPC[j], coords.w$xw[k], coords.datal[[samplelabs[z]]]$yPC[j], coords.w$yw[k])
                      }
                    }

                    mindist <- apply(distm, 1, min)
                    mindistsup <- apply(distm, 2, min)
                    if (estrange[w] < tempestrange) {
                      sup.reduced <- which(mindistsup <= ((tempestrange/sqrt(2)) * sdElim))
                    }
                    if (estrange[w] >= tempestrange) {
                      sup.reduced <- which(mindistsup <= ((estrange[w]/sqrt(2)) * sdElim))
                    }
                    nsup.reduced[z] <- length(sup.reduced)

                    if (nsup.reduced[z] == 0) {
                      coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]] <- coords.w[sup.reduced, ]
                    }

                    if (nsup.reduced[z] > 0) {
                      coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]] <- coords.w[sup.reduced, ]
                      coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]]$sample <- samplelabs[z]
                      colnames(coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]]) <- c("x.omega", "y.omega",
                                                                                                                             "sample")
                      coords.u[[paste(nsup.tot, "alternate_sample", samplelabs[z], sep = "", collapse = "")]] <- coords.u[[paste(nsup.tot, "alternate_sample",
                                                                                                                                 samplelabs[z], sep = "", collapse = "")]][, c(3, 1, 2)]
                    }


                    # need to calculate the average distance from each support point to the -- nearest support point
                    # --------------------------------------------------
                    supdistmat <- matrix(NA, nrow = nsup.tot, ncol = nsup.tot)
                    for (j in 1:nsup.tot) {
                      for (k in 1:nsup.tot) {
                        supdistmat[j, k] <- distance(xw[j], xw[k], yw[j], yw[k])
                      }
                    }
                    supdistmat <- t(apply(supdistmat, 1, sort))
                    supdistmat <- supdistmat[, -1]

                    maxdistm2vec <- rep(NA, 9)
                    maxdistm2vec[1] <- samplelabs[z]
                    maxdistm2vec[2] <- paste(nsup.tot, "alternate", sep = "", collapse = "")
                    maxdistm2vec[3] <- nrow.sup
                    maxdistm2vec[4] <- nsup.tot
                    maxdistm2vec[5] <- nsup.reduced[z]
                    maxdistm2vec[6] <- max(mindist)
                    maxdistm2vec[7] <- mean(supdistmat[, 1])


                    maxdistm2 <- rbind(maxdistm2, maxdistm2vec)

                    validstruct[z] <- ifelse(as.numeric(maxdistm2vec[6]) <= (estrange[w]/sqrt(2)), 1, 0)
                  }
                  totnsup.reduced <- sum(nsup.reduced)


                  flag <- ifelse((totnsup.reduced > thresh.nsup) | (sum(validstruct) == length(validstruct)), 0, 1)
                  nrow.sup <- nrow.sup + 1
                }
                maxdistm2 <- maxdistm2[-1, ]
                maxdistm2$nrow.sup <- as.numeric(maxdistm2$nrow.sup)
                maxdistm2$npoints <- as.numeric(maxdistm2$npoints)
                maxdistm2$npoints.reduced <- as.numeric(maxdistm2$npoints.reduced)
                maxdistm2$maxdistm <- as.numeric(maxdistm2$maxdistm)
                maxdistm2$meansupdist <- as.numeric(maxdistm2$meansupdist)


                maxdistm <- rbind(maxdistm1, maxdistm2)
                maxdistm$nrow.sup <- as.numeric(maxdistm$nrow.sup)
                maxdistm$npoints <- as.numeric(maxdistm$npoints)
                maxdistm$npoints.reduced <- as.numeric(maxdistm$npoints.reduced)
                maxdistm$maxdistm <- as.numeric(maxdistm$maxdistm)
                maxdistm$meansupdist <- as.numeric(maxdistm$meansupdist)

                mstructs <- matrix(nrow = length(unique(maxdistm$mstructure)), ncol = 2)
                mstructs <- as.data.frame(mstructs)
                colnames(mstructs) <- c("mstruct", "totnsup.reduced")
                mstructs$mstruct <- unique(maxdistm$mstructure)
                for (i in 1:nrow(mstructs)) {
                  mstructs[i, 2] <- sum(maxdistm$npoints.reduced[maxdistm$mstructure == mstructs[i, 1]])
                }
                mstructs <- mstructs[mstructs$totnsup.reduced <= thresh.nsup, ]

                maxdistm <- maxdistm[maxdistm$mstructure %in% mstructs$mstruct, ]
                # order by sum of support points
                mstructs <- mstructs[order(mstructs$totnsup.reduced), ]
                mstructs.ord <- mstructs$mstruct
                maxdistm$mstructure <- factor(maxdistm$mstructure, levels = mstructs.ord)

                maxdistm <- maxdistm[order(maxdistm$mstructure), ]

                # recommended SD
                maxdistm$recsd <- estrange[w]/sqrt(2)
                for (i in 1:nrow(maxdistm)) {
                  maxdistm$valid[i] <- ifelse(maxdistm$maxdistm[i] <= maxdistm$recsd[i] * sdWithin, 1, 0)
                }

                mstructs$nvalid <- NA
                for (i in 1:nrow(mstructs)) {
                  mstructs$nvalid[i] <- nrow(maxdistm[maxdistm$mstructure == mstructs$mstruct[i] & maxdistm$valid == 1, ])
                }
                validms <- mstructs$mstruct[mstructs$nvalid == nsamps]

                coordnames <- names(coords.u)


                if (estrange[w] > distthresh) {
                  suggest <- maxdistm[maxdistm$mstructure == "5x", ]
                  recommendSd <- max(suggest$maxdistm) * sdDefault
                  recommendStruct <- "5x"
                  defaultStatus <- "found valid structure, SD exceeded threshold"
                  nRowSupport <- NA
                  coordsU <- do.call(rbind.data.frame, coords.u[which(grepl("5x_sample", coordnames))])
                }

                if (estrange[w] <= distthresh) {
                  if (length(validms) == 0) {
                    defaultStatus <- "no valid structure"
                    nRowSupport <- NA

                    if (defaultStructure == "5x") {
                      suggest <- maxdistm[maxdistm$mstructure == "5x", ]
                      recommendSd <- max(suggest$maxdistm) * sdDefault
                      recommendStruct <- "5x"
                      coordsU <- do.call(rbind.data.frame, coords.u[which(grepl("5x_sample", coordnames))])
                    }
                    if (defaultStructure == "nextHighest") {
                      suggest <- maxdistm[maxdistm$mstructure == mstructs.ord[length(mstructs.ord)], ]
                      # suggest <- suggest[1, ]
                      if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == TRUE) {
                        nRowSupport <- suggest$nrow.sup[1]
                      }
                      recommendSd <- max(suggest$maxdistm) * sdDefault
                      recommendStruct <- mstructs.ord[length(mstructs.ord)]
                      coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(mstructs.ord[length(mstructs.ord)], coordnames))])
                    }
                  }


                  if (length(validms) > 0) {
                    suggest <- maxdistm[maxdistm$mstructure == validms[1], ]
                    if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == TRUE) {
                      nRowSupport <- suggest$nrow.sup[1]
                      recommendSd <- suggest$recsd[1]
                      recommendStruct <- suggest$mstructure[1]
                      defaultStatus <- "found valid structure"
                      coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(suggest$mstructure[1], coordnames))])
                    }
                    if (grepl(pattern = "alternate", x = suggest$mstructure[1]) == FALSE) {
                      nRowSupport <- NA
                      recommendSd <- suggest$recsd[1]
                      recommendStruct <- suggest$mstructure[1]
                      defaultStatus <- "found valid structure"
                      coordsU <- do.call(rbind.data.frame, coords.u[which(grepl(suggest$mstructure[1], coordnames))])
                    }
                  }
                }
                row.names(coordsU) <- seq(1:nrow(coordsU))

                outp[[paste0("structure_", w)]] <- (list(recommendStruct = recommendStruct, recommendSd = recommendSd, nSupportReduced = nrow(coordsU),
                                                         coordsU = coordsU, coordsData = coords.data, defaultStatus = defaultStatus, nRowSupport = nRowSupport))
              }  # spatialVar for loop
            }


            outp[["data"]] <- data
            outp[["subjectVar"]] <- subjectVar
            outp[["sampleVar"]] <- sampleVar
            outp[["spatialVar"]] <- spatialVar
            outp[["outcome"]] <- outcome
        }  # if(is.null(spatialVar)==FALSE)

    class(outp) <- "structureList"
    return(outp)

}


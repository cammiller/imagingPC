

#' Calculate the area under the curve (AUC) across m/z values on the raster
#' level.
#'
#' \code{createAUCData} calculates the AUC for sets of m/z values.
#'
#' @param data A dataset containing columns of m/z measurements.
#'
#' @param mzThreshold An m/z threshold that determines which m/z values
#' are combined when calculating the AUC.  See \code{Details} for more
#' information on how this works.
#'
#' @return A new dataset with the AUC calculations.  Any columns with column
#' names that do not start with an "X" and then a nubmer will be kept.  All
#' original m/z columns will be removed.
#'
#' @details This function is designed to work with data exported from
#' SCiLS Lab software.  When using SCiLS make sure to export the m/z values
#' as columns, meaning every column should represent a separate m/z value,
#' and every column name is the m/z value When the data are imported into R,
#' R will add an 'X' in front of each number.  Do not remove the 'X' in front
#' as this function will automatically do that.
#'
#' The \code{createAUCData} first searches through the m/z values and
#' determines which to combine based on the \code{mzThreshold}.  The function
#' goes one-by-one through the m/z values and determines which other
#' fragments are within the \code{mzThreshold}  Those within the mzThreshold
#' are considered to represent the same fragment and are used to
#' calculate the AUC.
#'
#' The \code{mzThreshold} must be smaller than the smallest m/z distance
#' between m/z values representing different fragments.  For example,
#' if the highest m/z representing fragment 1 is 2466.005859, and the
#' smallest m/z representing fragment 2 is 2466.878418, then the threshold
#' must be less than 0.872559 (2466.878418-2466.005859=0.872559) to
#' correctly distinguish between the two fragments.
#'
#' As a result of the way m/z values are combined, this function
#' should be used for data in which the m/z distances defining each
#' fragment are less than the m/z distances between fragments.  This is often
#' the case when the data that have gone through the peak picking
#' process.  If entire spectra for spot-level data are submitted to this
#' function, all m/z values will be combined.

createAUCData <- function(data, mzThreshold = 0.3, usedMergeIMSFiles=TRUE) {

  if(usedMergeIMSFiles==TRUE){
    nanames<-colnames(data)[1:3]
    numcnames<-as.numeric(colnames(data)[4:ncol(data)])
    for(i in 1:length(numcnames)){
      colnames(data)[i+3]<-paste0('X',numcnames[i])
    }
  }
  if(usedMergeIMSFiles==FALSE){
    # get colunm names for masses ---------------------------------------------
    cnames <- colnames(data)
    nanames <- rep(NA, length(cnames))
    for (i in 1:length(numcnames)) {
      nanames[i] <- as.numeric(strsplit(cnames[i], "X")[[1]][2])
    }
    nanames <- which(is.na(nanames))
    # data2<-data[,-nanames]

    cnames <- colnames(data[, -nanames])
    numcnames <- rep(NA, length(cnames))
    for (i in 1:length(numcnames)) {
      numcnames[i] <- as.numeric(strsplit(cnames[i], "X")[[1]][2])
    }
  }


    allc <- numcnames
    allc <- sort(allc)
    nmass <- length(allc)

    res <- matrix(NA, nrow = nmass, ncol = nmass)
    row.names(res) <- allc

    # determine which centroids are within threshold of one another
    for (i in 1:length(allc)) {
        for (j in 1:nmass) {
            # for each i, is mass j within thresh?
            res[i, j] <- ifelse(allc[j] <= (allc[i] + mzThreshold) & allc[j] >= (allc[i] - mzThreshold), 1, 0)
        }
    }
    # find the stop points
    stopp <- matrix(0, nrow = nmass, ncol = 2)
    stopp[, 1] <- seq(1, nmass)
    for (i in 1:(nmass - 1)) {
        # specify the stop points
        if (res[i, (i + 1)] == 0 & res[(i + 1), i] == 0) {
            stopp[i, 2] <- 1
        }
    }

    stopp <- stopp[stopp[, 2] == 1, ]
    if (is.vector(stopp) == TRUE) {
        stopvec <- c(0, stopp[1], nmass)
    }
    if (is.vector(stopp) == FALSE) {
        stopvec <- as.vector(c(0, stopp[, 1], nmass))
    }

    minires <- list()
    for (i in 1:(length(stopvec) - 1)) {
        minires[[i]] <- row.names(res[((stopvec[i] + 1):stopvec[i + 1]), ])
    }

    for (i in 1:length(minires)) {
        names(minires)[i] <- mean(as.numeric(minires[[i]]))
    }

    aucdat <- matrix(nrow = nrow(data), ncol = length(names(minires)))

    # calculate AUC -----------------------------------------------------------
    for (j in 1:length(minires)) {
        peaks <- paste("X", minires[[j]], sep = "")

        datatemp <- data[, peaks]

        npeaks <- length(peaks)

        aucs <- matrix(nrow = nrow(datatemp), ncol = (npeaks - 1))
          for (i in 1:ncol(aucs)) {
            aucs[, i] <- (datatemp[, (i)] + datatemp[, (i + 1)]) * 0.5 * (as.numeric(gsub("X", "", colnames(datatemp)[i + 1])) - as.numeric(gsub("X",
                                                                                                                                                 "", colnames(datatemp)[i])))
          }

        aucdat[, j] <- rowSums(aucs)

    }
    colnames(aucdat) <- paste("X", round(as.numeric(names(minires)), digits = 4), sep = "")

    aucdat <- cbind(data[, nanames], aucdat)

    return(aucdat)



}

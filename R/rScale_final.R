
#' Rescale the data
#'
#' \code{rScale} rescales the coordinates of the data so that the distance
#' between observations is 1.  It also centers the set of coordinates for each
#' sample at the origin.
#'
#' @param data The dataset to be rescaled.  The dataset must be a data frame.
#' @param subjectVar A character string specifying the subject variable, if it
#' exists.  This should be specified for paired data only, meaning data in which
#' there is more than one sample for at least one subject.
#' @param sampleVar A character string specifying the sample variable.  This
#' should be a variable with unique values for each sample.
#' @param xCoord A character string specifying the name of the x-coordinate
#' variable.
#' @param yCoord A character string specifying the name of the y-coordinate
#' variable.

#' @return A list including the rescaled dataset and the names of the subject
#' and sample variables.  The subject and sample variables are carried forward
#' to subsequent functions so that they only have to be input once.
#'
#' \describe{
#'   \item{\code{data}}{The dataset, appended with new x- and y-coordinates
#'   to be used in further calculations.}
#'   \item{\code{subjectVar}}{A character string denoting the subject
#'   variable.}
#'   \item{\code{sampleVar}}{A character string denoting the sample variable.}
#' }
#'
#' @examples
#' data("TAMdata")
#' TAMdata <- rScale(TAMdata, subjectVar = 'subject', sampleVar = 'ROI',
#'                   xCoord = 'x', yCoord = 'y')



rScale <- function(data, subjectVar = NULL, sampleVar, xCoord, yCoord) {
    if (is.data.frame(data) == FALSE)
        stop("Data must be a data frame.")
    if (is.character(sampleVar) == FALSE)
        stop("The variable name for sampleVar must be a character string.")
    if (is.character(xCoord) == FALSE)
        stop("The variable name for xCoord must be a character string.")
    if (is.character(yCoord) == FALSE)
        stop("The variable name for yCoord must be a character string.")
    if (is.numeric(data[[xCoord]]) == "FALSE")
        stop("The variable xCoord must be numeric")
    if (is.numeric(data[[yCoord]]) == "FALSE")
        stop("The variable yCoord must be numeric")
    if (is.null(data) == TRUE)
        stop("You must provide a dataset")
    if (is.null(sampleVar) == TRUE)
        stop("You must provide a sample variable")


    # check if sample labels are unique across subjects if not, reassign with a new sample variable
    if (is.null(subjectVar) == FALSE) {
        unsamplab <- list()
        sublabs <- unique(data[[subjectVar]])
        nsubs <- length(sublabs)
        for (i in 1:nsubs) {
            unsamplab[[i]] <- unique(data[[sampleVar]][data[[subjectVar]] == sublabs[i]])
        }
        allsamplab <- paste(unlist(unsamplab))
        unsamponly <- unique(allsamplab)

        if (length(allsamplab) == length(unsamponly)) {
        }
        if (length(allsamplab) != length(unsamponly)) {
            data$PCsampleVar <- NA
            data$PCsampleVar[1] <- 1
            for (i in 2:nrow(data)) {
                data$PCsampleVar[i] <- ifelse(data[[sampleVar]][i] == data[[sampleVar]][i - 1], data$PCsampleVar[i - 1], (data$PCsampleVar[i - 1] +
                  1))
            }
            sampleVar <- "PCsampleVar"
        }
    }



    # check if there is more than one sample per subject, given a subjectVar is provided
    if (is.null(subjectVar) == FALSE) {
        nsampmat <- matrix(NA, nrow = length(unique(data[[subjectVar]])), ncol = 2)
        nsampmat[, 1] <- unique(data[[subjectVar]])
        for (i in 1:nrow(nsampmat)) {
            nsampmat[i, 2] <- length(unique(data[[sampleVar]][data[[subjectVar]] == nsampmat[i, 1]]))
        }
        if (any(nsampmat[, 2] > 1) == FALSE) {
            subjectVar <- NULL
        }
    }


    if (is.null(subjectVar) == TRUE) {
        data <- data[order(data[[sampleVar]]), ]

        if (is.numeric(data[[sampleVar]]) == FALSE) {
            data$PCsampleVar <- NA
            data$PCsampleVar[1] <- 1
            for (i in 2:nrow(data)) {
                data$PCsampleVar[i] <- ifelse(data[[sampleVar]][i] == data[[sampleVar]][i - 1], data$PCsampleVar[i - 1], (data$PCsampleVar[i - 1] +
                  1))
            }
            sampleVar <- "PCsampleVar"
        }
    }

    if (is.null(subjectVar) == FALSE) {
        data <- data[order(data[[subjectVar]], data[[sampleVar]]), ]

        if (is.numeric(data[[subjectVar]]) == FALSE) {
            data$PCsubjectVar <- NA
            data$PCsubjectVar[1] <- 1
            for (i in 2:nrow(data)) {
                data$PCsubjectVar[i] <- ifelse(data[[subjectVar]][i] == data[[subjectVar]][i - 1], data$PCsubjectVar[i - 1], (data$PCsubjectVar[i -
                  1] + 1))
            }
            subjectVar <- "PCsubjectVar"
        }

        if (is.numeric(data[[sampleVar]]) == FALSE) {
            data$PCsampleVar <- NA
            data$PCsampleVar[1] <- 1
            for (i in 2:nrow(data)) {
                data$PCsampleVar[i] <- ifelse(data[[sampleVar]][i] == data[[sampleVar]][i - 1], data$PCsampleVar[i - 1], (data$PCsampleVar[i - 1] +
                  1))
            }
            sampleVar <- "PCsampleVar"
        }
    }

    # data$spotind<-seq(1:nrow(data))

    minx <- min(data[[xCoord]])
    miny <- min(data[[yCoord]])

    data$xnew <- data[[xCoord]] - minx
    data$ynew <- data[[yCoord]] - miny

    data$xnew <- round(data$xnew)
    data$ynew <- round(data$ynew)

    samplelabs <- unique(data[[sampleVar]])
    nsamps <- length(samplelabs)
    mindist <- rep(NA, nsamps)
    for (i in 1:nsamps) {
        mindist[i] <- min(dist(cbind(data$xnew[data[[sampleVar]] == samplelabs[i]], data$ynew[data[[sampleVar]] == samplelabs[i]])))
    }
    mindist <- min(mindist)

    data$xnewer <- data$xnew/mindist
    data$ynewer <- data$ynew/mindist


    # center each sample at origin

    data$xPC <- NA
    data$yPC <- NA
    for (i in 1:nsamps) {
        data$xPC[data[[sampleVar]] == samplelabs[i]] <- data$xnewer[data[[sampleVar]] == samplelabs[i]] - mean(range(data$xnewer[data[[sampleVar]] ==
            samplelabs[i]]))
        data$yPC[data[[sampleVar]] == samplelabs[i]] <- data$ynewer[data[[sampleVar]] == samplelabs[i]] - mean(range(data$ynewer[data[[sampleVar]] ==
            samplelabs[i]]))
    }

    data <- subset(data, select = -c(xnew, xnewer, ynew, ynewer))

    outlist <- list(data = data, subjectVar = subjectVar, sampleVar = sampleVar)

    class(outlist) <- "rScaleList"

    return(outlist)

}

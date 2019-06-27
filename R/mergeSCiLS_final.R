

#' Merge spot-level data and region spots from SCiLS
#'
#' \code{mergeIMSFiles} merges the two files needed to create spot-level data with
#' coordinates.  The first is the spot-level data without coordinates, and the
#' second is the region spots with coordinates.  We will use the terms spots and
#' rasters interchangeably.  Though data are collected over circular areas
#' dependent on the laser diameter, the rasterization process makes data
#' representation easier.
#'
#' @param spotData A dataset containing the spot-level data in which the
#' columns represent m/z values, and the rows represent spots.
#'
#' @param regionSpots A dataset containing the region spots.  This should be
#' an exported file that only contains the spot identifiers and
#' the x- and y-coordinates.
#'
#' @param trimTops A TRUE/FALSE variable indicating whether or not the top
#' rows of each dataset should be trimmed to remove metadata before merging.
#' If metadata exists in the top rows then trimTops should be set to TRUE
#' to remove it before proceeding to the next step.  Otherwise the metadata
#' will ineterfere with downstream steps.

#' @return A merged data frame with spot-level data and corresponding
#' coordinates.
#'


mergeIMSFiles <- function(spotData, regionSpots, trimTops = TRUE) {

    if (trimTops == TRUE) {
        # find row at which metadata stops ------------------------------------ region spots
        flag <- 1
        iter <- 0
        while (flag) {
            iter <- iter + 1
            flag <- ifelse(substr(regionSpots[iter, 1], 1, 1) == "#", 1, 0)
        }
        regionSpots <- regionSpots[-c(1:(iter - 1)), ]

        flag <- 1
        iter <- 0
        while (flag) {
            iter <- iter + 1
            flag <- ifelse(as.character(regionSpots[iter]) == " y", 0, 1)
        }

        estcol <- iter
        regionSpots <- matrix(as.vector(as.character(regionSpots)), nrow = length(regionSpots)/estcol, ncol = estcol, byrow = TRUE)
        nspots_plusone <- nrow(regionSpots)

        cnames <- regionSpots[1, ]
        regionSpots <- as.data.frame(regionSpots)
        colnames(regionSpots) <- cnames
        regionSpots <- regionSpots[-1, ]
        for (i in 1:estcol) {
            regionSpots[, i] <- as.numeric(as.character(regionSpots[, i]))
        }



        # spot-level data find row at which metadata stops ------------------------------------
        flag <- 1
        iter <- 0
        while (flag) {
            iter <- iter + 1
            flag <- ifelse(substr(spotData[iter, 1], 1, 1) == "#", 1, 0)
        }
        spotData <- spotData[-c(1:(iter - 1)), ]

        spotData <- matrix(as.vector(as.character(spotData)), nrow = nspots_plusone, ncol = (length(spotData)/nspots_plusone), byrow = TRUE)

        cnames <- spotData[1, ]
        spotData <- as.data.frame(spotData)
        colnames(spotData) <- cnames
        spotData <- spotData[-1, ]
        for (i in 2:ncol(spotData)) {
            spotData[, i] <- as.numeric(as.character(spotData[, i]))
        }
    }

    colnames(spotData)[1] <- "Spot.index"
    colnames(regionSpots)[1] <- "Spot.index"

    spotData$spot2 <- NA
    for (i in 1:nrow(spotData)) {
        spotData$spot2[i] <- as.numeric(gsub("Spot ", "", spotData[i, 1]))
    }


    data <- merge(x = regionSpots, y = spotData, by.x = "Spot.index", by.y = "spot2")
    data <- subset(data, select = -c(Spot.index.y))
    colnames(data)[2:3]<-c('x','y')
    return(data)
}

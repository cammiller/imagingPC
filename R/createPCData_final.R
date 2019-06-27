

#' Create the information for running models
#'
#' The \code{createPCData} function generates all the information needed to
#' fit a process convolution model in a Bayesian setting.
#'
#' @param structureObj An object of class \code{structureList}.  This object is the
#' result of using the \code{\link{chooseStructures}} function.
#' @param covariates A concatenated character string specifying the covariates for
#' which modeling information will be generated.  Not all covariates have to be
#' included in the model that will be fit, so it is a good idea to include all
#' covariates that may be of interest.
#' @param covariateTypes A concatenated character string specifying the type of
#' variable for each covariate given in the covariates argument.  The three options
#' are "binary", "continuous", and "categorical".  Every variable in the covariates
#' argument must have a corresponding covariate type, where the nth covariate type
#' corresponds to the nth covariate in the covariates argument.
#' @param covariateLevels A concatenated character string specifying the level at
#' which each covariates exists.  The options are "subject", "sample", and
#' "raster".  As an example, for a study examining differences between tumor and
#' non-tumor samples, the level would be 'sample' since the covariate changes
#' between samples.  A level must be provided for every covariate given in the
#' argument covariates.
#' @param trimData A TRUE/FALSE argument specifying whether or not the data
#' should be trimmed.  If trimData=TRUE, the data will be subset to only include
#' the variables required to fit a model.
#'
#' @return A list containing the model and other information supplied to the
#' \code{createPCData} function.
#'
#' \describe{
#'   \item{\code{data}}{A data frame containing the data}
#'   \item{\code{nSubjs}}{The number of subjects}
#'   \item{\code{nSamps}}{The number of samples}
#'   \item{\code{cNSampsPerSubj}}{A cumulative vector of the number of samples
#'   per subject.  If no subject variable was provided to the rScale function
#'   then this will be NULL.}
#'   \item{\code{cNRastPerSamp}}{A cumulative vector of the number of rasters
#'   per sample.}
#'   \item{\code{totalRasters}}{A numeric value for the total number of
#'   rasters.}
#'   \item{\code{covs}}{A list of lists of covariate information.  For each
#'   variable in the covariates argument, a list of covariate information is
#'   created that includes elements \code{covariate} (the name of covariate),
#'   \code{type} (the type of covariate given in the covariateTypes argument),
#'   \code{level} (the level of the covariate given in the covariateLevels
#'   argument), \code{info} (the data corresponding to the covariate, given as
#'   data frames for subject-level and sample-level covariates and vectors for
#'   raster-level covariates), and \code{mapping} (a data frame mapping the
#'   given covariate values to new ones that are used in modeling).}
#'   \item{\code{KMat}}{A matrix of density values generated from the
#'   smoothing kernel function.  The matrix has rows equal to the number of
#'   rows in the dataset and columns equal to the maximum number of support
#'   sites for any of the samples.  Missing cells indicate that the support
#'   site corresponding to that column was removed for the sample
#'   corresponding to that data row.}
#'   \item{\code{rastersPerVar}}{A data frame showing the number of rasters,
#'   per level of the spatial variable, for each sample.  The data frame also
#'   shows the cumulative number of rasters.  If no spatial variable is given
#'   then this is NULL.}
#'   \item{\code{nSupportSites}}{The number of support sites per sample.  If
#'   no spatial variable was given then this is a vector.  If a spatial
#'   variable is given then this is a data frame showing the number of support
#'   sites per sample, for each level of the spatial variable.}
#'   \item{\code{nObs}}{A data frame giving the number of rasters per sample.}
#'   \item{\code{gT0SupportSites}}{If no spatial variable is provided, then
#'   this is NULL.  If a spatial variable is given the this is a list of
#'   vectors, one for each level of the spatial variable, where each vector
#'   gives the samples numbers with more than one support site for the
#'   corresponding level of the spatial variable.  If no spatial variable
#'   is given then this is NULL.}
#'   \item{\code{nVarLevels}}{The number of levels of the spatial variable.
#'   If no spatial variable is given then this is NULL.}
#'   \item{\code{subjectVar}}{A character string specifying the subject
#'   variable.}
#'   \item{\code{sampleVar}}{A character string specifying the sample
#'   variable.}
#'   \item{\code{spatialVar}}{A character string specifying the spatial
#'   variable.}
#'   \item{\code{covariates}}{A concatenated character string of the
#'   covariate names.}
#'   \item{\code{covariateTypes}}{A concatenated character string of the
#'   covariate types (binary, categorical, continuous) corresponding to the
#'   covariates.}
#'   \item{\code{covariateLevels}}{A concatenated character string of the
#'   covariate levels (subject, sample, raster) corresponding to the
#'   covariates.}
#'   \item{\code{outcome}}{A character string specifying the name of the
#'   variable to be modeled.}
#'   \item{\code{recStructures}}{An indicator where 0 means no structure
#'   was chosen (if estimmated range=0) and 1 means a structure was chosen.}
#' }
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


createPCData <- function(structureObj, covariates, covariateTypes, covariateLevels, trimData = FALSE) {

    if (class(structureObj) != "structureList")
        stop("structureObj must of class structureList.  This object must be created using the chooseStructures function.")
    if (is.character(covariates) == FALSE)
        stop("covariates must be a character string.")
    if (is.character(covariateTypes) == FALSE)
        stop("covariateTypes must be a character string.")
    if (is.character(covariateLevels) == FALSE)
        stop("covariateLevels must be a character string.")
    if (length(covariates) != length(covariateLevels))
        stop("The lengths of covariates and covariateLevels must be equal")
    if (any((covariateLevels %in% c("subject", "sample", "raster")) == FALSE))
        stop("The only options for covariateLevels are \"subject\", \"sample\", and \"raster\".")
    if (any((covariateTypes %in% c("continuous", "binary", "categorical")) == FALSE))
        stop("The only options for covariateTypes are \"continuous\", \"binary\", and \"categorical\".")

    if (all(covariates %in% colnames(structureObj[["data"]])) == FALSE) {
        covcheck <- rep(NA, length(covariates))
        for (i in 1:length(covariates)) {
            covcheck[i] <- ifelse(covariates[i] %in% colnames(structureObj[["data"]]), 1, 0)
        }
        notindat <- which(covcheck == 0)
        stop(paste0("The following covariates are not in the dataset: ", paste(covariates[notindat], sep = "", collapse = ", "), "."))
    }

    data <- structureObj[["data"]]
    subjectVar <- structureObj[["subjectVar"]]
    sampleVar <- structureObj[["sampleVar"]]
    spatialVar <- structureObj[["spatialVar"]]
    outcome <- structureObj[["outcome"]]

    distance <- function(x1, x2, y1, y2) {
        sqrt((x1 - x2)^2 + (y1 - y2)^2)
    }


    ## subject variable and spatial variable ##
    if (is.null(subjectVar) == FALSE & is.null(spatialVar) == FALSE) {
        if (trimData == TRUE) {
          if(is.null(slideVar)==TRUE){
            data <- data[, which(colnames(data) %in% unique(c(subjectVar, sampleVar, spatialVar, covariates, outcome, "xPC", "yPC")))]
          }
          if(is.null(slideVar)==FALSE){
            data <- data[, which(colnames(data) %in% unique(c(subjectVar, sampleVar, spatialVar, covariates, outcome, slideVar, "xPC", "yPC")))]
          }
        }

        data <- data[order(data[[spatialVar]]), ]
        spatvarlevels <- unique(data[[spatialVar]])
        nvarlevels <- length(spatvarlevels)
        data <- data[order(data[[subjectVar]], data[[sampleVar]], data[[spatialVar]]), ]
        samplelabs <- unique(data[[sampleVar]])
        nsamps <- length(samplelabs)
        nobs_df <- matrix(nrow = nsamps, ncol = 3)
        nobs_df <- as.data.frame(nobs_df)
        colnames(nobs_df) <- c("subject", "sample", "numrasters")
        nobs_df$sample <- samplelabs
        for (i in 1:nsamps) {
            nobs_df$subject[i] <- unique(data[[subjectVar]][data[[sampleVar]] == nobs_df[i, 2]])
            nobs_df$numrasters[i] <- nrow(data[data[[sampleVar]] == nobs_df[i, 2], ])
        }

        recstructure<-rep(NA,nvarlevels)
        for(i in 1:nvarlevels){
          recstructure[i]<-ifelse(structureObj[[i]]$recommendStruct=='noStructure',0,1)
        }

        ### cumulative sample vector
        subjlabs <- unique(data[[subjectVar]])
        nsubjs <- length(unique(data[[subjectVar]]))
        nsampssub <- rep(NA, nsubjs)
        for (i in 1:nsubjs) {
            nsampssub[i] <- length(unique(data[[sampleVar]][data[[subjectVar]] == samplelabs[i]]))
        }
        cnsampspersub <- c(0, cumsum(nsampssub))

        ### cumulative raster vector
        cnrastpersamp <- c(0, cumsum(nobs_df$numrasters))

        ### covariates need separate matrix for each covariate secondary is on the sample level
        covs <- list()
        covtemp <- list()
        for (i in 1:length(covariates)) {
            ## warnings
            if (covariateTypes[i] == "continuous") {
                if (is.numeric(data[[covariates[i]]]) == FALSE)
                  stop(paste0("Covariate ", covariates[i], " is not continuous."))
            }
            if (covariateTypes[i] == "binary") {
                if (length(unique(data[[covariates[i]]])) != 2)
                  stop(paste0("Covariate ", covariates[i], " is not binary."))
            }

            if (covariateTypes[i] == "continuous") {
                covtemp[["covariate"]] <- covariates[i]
                covtemp[["level"]] <- covariateLevels[i]

                if (covariateLevels[i] == "subject") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 2)
                  covtemp[["info"]][, 1] <- subjlabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "sample") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                  covtemp[["info"]][, 1] <- samplelabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "raster") {
                  covtemp[["info"]] <- data[[covariates[i]]]
                }
                covs[[i]] <- covtemp

            }



            if (covariateTypes[i] %in% c("binary", "categorical")) {
                if (covariateTypes[i] == "binary") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "subject") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- subjlabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])

                  }

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("sample", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- rep(NA, nrow(data))
                    for (j in 1:nrow(data)) {
                      covtemp[["info"]][j] <- ifelse(data[[covariates[i]]][j] == covvalues[1], 0, 1)
                    }
                  }
                }



                if (covariateTypes[i] == "categorical") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "subject") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- subjlabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- data[[covariates[i]]]

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(data), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(data)) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                }

            }
            covs[[paste0("covariate", i, "_", covariates[i])]] <- list(covariate = covtemp[["covariate"]], type = covariateTypes[i], level = covtemp[["level"]],
                info = covtemp[["info"]], mapping = covtemp[["mapping"]])


        }  #i loop for each covariate

        ### total number of rasters
        totrast <- max(cnrastpersamp)

        ### number of support sites per sample (from choose.structure)
        nsups <- matrix(nrow = nsamps, ncol = nvarlevels)
        nsups <- as.data.frame(nsups)
        for (i in 1:nvarlevels) {
            colnames(nsups)[i] <- paste0("nsups.", i)
        }

        for (j in 1:nvarlevels) {
          if(structureObj[[j]]$recommendStruct!="noStructure"){
            for (i in 1:nsamps) {
              nsups[i, j] <- nrow(structureObj[[j]]$coordsU[structureObj[[j]]$coordsU$sample == nobs_df$sample[i], ])
            }
          }
          if(structureObj[[j]]$recommendStruct=="noStructure"){
            for(i in 1:nsamps){
              nsups[i,j]<-ifelse(nrow(data[data[[spatialVar]]==spatvarlevels[j] & data[[sampleVar]]==nobs_df$sample[i],])==0,0,1)
            }
          }
        }
        nobs_df <- cbind(nobs_df, nsups)

        # need to know for which samples there are no observations for a covariate
        nosup_var <- list()
        for (i in 1:nvarlevels) {
            nosup_var[[i]] <- nobs_df$sample[which(nsups[, i] == 0)]
        }
        gt0sup_var <- list()
        for (i in 1:nvarlevels) {
            gt0sup_var[[i]] <- nobs_df$sample[which(nsups[, i] > 0)]
        }

        ### cumulative raster vector per covariate level
        rastpervar <- matrix(nrow = nsamps * nvarlevels, ncol = 5)
        rastpervar <- as.data.frame(rastpervar)
        colnames(rastpervar) <- c("sample", "cov.level", "numraster", "start", "stop")
        rastpervar$sample <- rep(nobs_df$sample, each = nvarlevels)
        rastpervar$cov.level <- rep(spatvarlevels)
        for (i in 1:nrow(rastpervar)) {
            rastpervar$numraster[i] <- nrow(data[data[[sampleVar]] == rastpervar$sample[i] & data[[spatialVar]] == rastpervar$cov.level[i], ])
        }
        if (rastpervar$numraster[1] == 0) {
            rastpervar$start[1] <- 0
        }
        if (rastpervar$numraster[1] > 0) {
            rastpervar$start[1] <- 1
        }
        rastpervar$stop[1] <- rastpervar$numraster[1]
        for (i in 2:nrow(rastpervar)) {
            if (rastpervar$numraster[i] > 0) {
                rastpervar$start[i] <- rastpervar$stop[i - 1] + 1
                rastpervar$stop[i] <- rastpervar$start[i] + (rastpervar$numraster[i] - 1)
            }
            if (rastpervar$numraster[i] == 0) {
                rastpervar$start[i] <- rastpervar$stop[i - 1]
                rastpervar$stop[i] <- rastpervar$stop[i - 1]
            }
        }

        zeroranges<-rep(NA,nvarlevels)
        for(i in 1:nvarlevels){
          zeroranges[i]<-ifelse(structureObj[[i]]$recommendStruct=="noStructure",1,0)
        }

        if(sum(zeroranges)==length(zeroranges)){
          Kmat<-NULL
        }

        if(sum(zeroranges)<length(zeroranges)){
          Kmat.all <- list()
          for (z in 1:length(spatvarlevels)) {
            if(structureObj[[z]]$recommendStruct!="noStructure"){
              for (i in 1:length(gt0sup_var[[z]])) {
                ucoords <- structureObj[[z]]$coordsU[structureObj[[z]]$coordsU$sample == gt0sup_var[[z]][i], ]
                datacoords <- structureObj[[z]]$coordsData[structureObj[[z]]$coordsData$sample == gt0sup_var[[z]][i], ]

                dist.mat <- matrix(NA, nrow = nrow(datacoords), ncol = nrow(ucoords))
                Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]] <- matrix(0, nrow = nrow(datacoords), ncol = nrow(ucoords))
                for (j in 1:ncol(dist.mat)) {
                  for (k in 1:nrow(dist.mat)) {
                    dist.mat[k, j] <- distance(ucoords$x.omega[j], datacoords$xPC[k], ucoords$y.omega[j], datacoords$yPC[k])
                    Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]][k, j] <- dnorm(0, dist.mat[k, j], sd = structureObj[[z]]$recommendSd)
                  }
                }
              }
            }

            if(structureObj[[z]]$recommendStruct=="noStructure"){
              for (i in 1:length(gt0sup_var[[z]])){
                Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]]<-matrix(NA,nrow=nrow(structureObj[[z]]$coordsData[structureObj[[z]]$coordsData$sample == gt0sup_var[[z]][i], ]), ncol=max(nsups))
              }

            }
          }

          rastpervar <- rastpervar[rastpervar$numraster > 0, ]

          Kmat <- as.data.frame(Kmat.all[[paste0("spatcovlevel", which(spatvarlevels == rastpervar$cov.level[1]), "_sample", rastpervar$sample[1])]])
          for (i in 2:nrow((rastpervar))) {
            Kmat <- rbind.fill(Kmat, as.data.frame(Kmat.all[[paste0("spatcovlevel", which(spatvarlevels == rastpervar$cov.level[i]), "_sample",
                                                                    rastpervar$sample[i])]]))
          }
        }

        outlist <- list(data = data, nSubjs = nsubjs, nSamps = nsamps, cNSampsPerSubj = cnsampspersub, cNRastPerSamp = cnrastpersamp, totalRasters = totrast,
            covs = covs, KMat = Kmat, rastersPerVar = rastpervar, nSupportSites = nsups, nObs = nobs_df, gT0SupportSites = gt0sup_var, nVarLevels = nvarlevels,
            subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, covariates = covariates, covariateTypes = covariateTypes, covariateLevels = covariateLevels,
            outcome = outcome, recStructures=recstructure)
        class(outlist) <- "PCDataList"
        return(outlist)

    }





    ## subject variable only ##
    if (is.null(subjectVar) == FALSE & is.null(spatialVar) == TRUE) {
        if (trimData == TRUE) {
            data <- data[, which(colnames(data) %in% unique(c(subjectVar, sampleVar, covariates, outcome, "xPC", "yPC")))]
        }

        data <- data[order(data[[subjectVar]], data[[sampleVar]]), ]
        subjlabs <- unique(data[[subjectVar]])
        nsubjs <- length(subjlabs)
        samplelabs <- unique(data[[sampleVar]])
        nsamps <- length(samplelabs)
        nobs_df <- matrix(nrow = nsamps, ncol = 3)
        nobs_df <- as.data.frame(nobs_df)
        colnames(nobs_df) <- c("subject", "sample", "numrasters")
        nobs_df$sample <- samplelabs
        for (i in 1:nsamps) {
            nobs_df$subject[i] <- unique(data[[subjectVar]][data[[sampleVar]] == nobs_df[i, 2]])
            nobs_df$numrasters[i] <- nrow(data[data[[sampleVar]] == nobs_df[i, 2], ])
        }

        nsampssub <- rep(NA, nsubjs)
        for (i in 1:nsubjs) {
            nsampssub[i] <- length(unique(data[[sampleVar]][data[[subjectVar]] == samplelabs[i]]))
        }

        ### cumulative sample vector
        nsampspersub <- rep(NA, nsubjs)
        for (i in 1:nsubjs) {
            nsampspersub[i] <- nrow(nobs_df[nobs_df$subject == subjlabs[i], ])
        }
        cnsampspersub <- c(0, cumsum(nsampspersub))

        ### cumulative raster vector
        cnrastpersamp <- c(0, cumsum(nobs_df$numrasters))

        ### covariates need separate matrix for each covariate secondary is on the sample level
        covs <- list()
        covtemp <- list()
        for (i in 1:length(covariates)) {
            ## warnings
            if (covariateTypes[i] == "continuous") {
                if (is.numeric(data[[covariates[i]]]) == FALSE)
                  stop(paste0("Covariate ", covariates[i], " is not continuous."))
            }
            if (covariateTypes[i] == "binary") {
                if (length(unique(data[[covariates[i]]])) != 2)
                  stop(paste0("Covariate ", covariates[i], " is not binary."))
            }

            if (covariateTypes[i] == "continuous") {
                covtemp[["covariate"]] <- covariates[i]
                covtemp[["level"]] <- covariateLevels[i]

                if (covariateLevels[i] == "subject") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 2)
                  covtemp[["info"]][, 1] <- subjlabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "sample") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                  covtemp[["info"]][, 1] <- samplelabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "raster") {
                  covtemp[["info"]] <- data[[covariates[i]]]
                }
                covs[[i]] <- covtemp

            }



            if (covariateTypes[i] %in% c("binary", "categorical")) {
                if (covariateTypes[i] == "binary") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "subject") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- subjlabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])

                  }

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("sample", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- rep(NA, nrow(data))
                    for (j in 1:nrow(data)) {
                      covtemp[["info"]][j] <- ifelse(data[[covariates[i]]][j] == covvalues[1], 0, 1)
                    }
                  }
                }



                if (covariateTypes[i] == "categorical") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "subject") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsubjs, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- subjlabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[subjectVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- data[[covariates[i]]]

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(data), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(data)) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                }

            }
            covs[[paste0("covariate", i, "_", covariates[i])]] <- list(covariate = covtemp[["covariate"]], type = covariateTypes[i], level = covtemp[["level"]],
                info = covtemp[["info"]], mapping = covtemp[["mapping"]])

        }  #i loop for each covariate

        ### total number of rasters
        totrast <- max(cnrastpersamp)

        if(structureObj[[1]]$recommendStruct=='noStructure'){
          nsups <- rep(0, nsamps)
          Kmat<-NULL
        }

        if(structureObj[[1]]$recommendStruct!='noStructure'){
          ### number of support sites per sample (from choose.structure)
          nsups <- rep(NA, nsamps)
          for (i in 1:nsamps) {
            nsups[i] <- nrow(structureObj[[1]]$coordsU[structureObj[[1]]$coordsU$sample == samplelabs[i], ])
          }
          nobs_df <- cbind(nobs_df, nsups)

          # gt0sup_var<-list() gt0sup_var[[1]]<-nobs_df$sample[which(nsups[i]>0)]

          Kmat.all <- list()
          for (i in 1:nsamps) {
            ucoords <- structureObj[[1]]$coordsU[structureObj[[1]]$coordsU$sample == samplelabs[i], ]
            datacoords <- structureObj[[1]]$coordsData[structureObj[[1]]$coordsData$sample == samplelabs[i], ]

            dist.mat <- matrix(NA, nrow = nrow(datacoords), ncol = nrow(ucoords))
            Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]] <- matrix(0, nrow = nrow(datacoords), ncol = nrow(ucoords))
            for (j in 1:ncol(dist.mat)) {
              for (k in 1:nrow(dist.mat)) {
                dist.mat[k, j] <- distance(ucoords$x.omega[j], datacoords$xPC[k], ucoords$y.omega[j], datacoords$yPC[k])
                Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]][k, j] <- dnorm(0, dist.mat[k, j], sd = structureObj[[1]]$recommendSd)
              }
            }
          }

          Kmat <- as.data.frame(Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[1])]])
          for (i in 2:nsamps) {
            Kmat <- rbind.fill(Kmat, as.data.frame(Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]]))
          }
        }

        recstructure<-ifelse(structureObj[[1]]$recommendStruct=='noStructure',0,1)

        rastpervar <- NULL
        nvarlevels <- NULL
        gt0sup_var <- NULL
        outlist <- list(data = data, nSubjs = nsubjs, nSamps = nsamps, cNSampsPerSubj = cnsampspersub, cNRastPerSamp = cnrastpersamp, totalRasters = totrast,
            covs = covs, KMat = Kmat, rastersPerVar = rastpervar, nSupportSites = nsups, nObs = nobs_df, gT0SupportSites = gt0sup_var, nVarLevels = nvarlevels,
            subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, covariates = covariates, covariateTypes = covariateTypes, covariateLevels = covariateLevels,
            outcome = outcome, recStructures=recstructure)
        class(outlist) <- "PCDataList"
        return(outlist)

    }





    ## spatial variable only ##
    if (is.null(subjectVar) == TRUE & is.null(spatialVar) == FALSE) {
        if (trimData == TRUE) {
            data <- data[, which(colnames(data) %in% unique(c(sampleVar, spatialVar, covariates, outcome, "xPC", "yPC")))]
        }

        data <- data[order(data[[spatialVar]]), ]
        spatvarlevels <- unique(data[[spatialVar]])
        nvarlevels <- length(spatvarlevels)
        data <- data[order(data[[sampleVar]], data[[spatialVar]]), ]
        samplelabs <- unique(data[[sampleVar]])
        nsamps <- length(samplelabs)
        nobs_df <- matrix(nrow = nsamps, ncol = 2)
        nobs_df <- as.data.frame(nobs_df)
        colnames(nobs_df) <- c("sample", "numrasters")
        nobs_df$sample <- samplelabs
        for (i in 1:nsamps) {
            nobs_df$numrasters[i] <- nrow(data[data[[sampleVar]] == nobs_df[i, 1], ])
        }

        recstructure<-rep(NA,nvarlevels)
        for(i in 1:nvarlevels){
          recstructure[i]<-ifelse(structureObj[[i]]$recommendStruct=='noStructure',0,1)
        }

        ### cumulative raster vector
        cnrastpersamp <- c(0, cumsum(nobs_df$numrasters))

        ### covariates need separate matrix for each covariate secondary is on the sample level
        covs <- list()
        covtemp <- list()
        for (i in 1:length(covariates)) {
            ## warnings
            if (covariateTypes[i] == "continuous") {
                if (is.numeric(data[[covariates[i]]]) == FALSE)
                  stop(paste0("Covariate ", covariates[i], " is not continuous."))
            }
            if (covariateTypes[i] == "binary") {
                if (length(unique(data[[covariates[i]]])) != 2)
                  stop(paste0("Covariate ", covariates[i], " is not binary."))
            }

            if (covariateTypes[i] == "continuous") {
                covtemp[["covariate"]] <- covariates[i]
                covtemp[["level"]] <- covariateLevels[i]

                if (covariateLevels[i] == "sample") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                  covtemp[["info"]][, 1] <- samplelabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "raster") {
                  covtemp[["info"]] <- data[[covariates[i]]]
                }
                covs[[i]] <- covtemp

            }



            if (covariateTypes[i] %in% c("binary", "categorical")) {
                if (covariateTypes[i] == "binary") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("sample", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- rep(NA, nrow(data))
                    for (j in 1:nrow(data)) {
                      covtemp[["info"]][j] <- ifelse(data[[covariates[i]]][j] == covvalues[1], 0, 1)
                    }
                  }
                }



                if (covariateTypes[i] == "categorical") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- data[[covariates[i]]]

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(data), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(data)) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                }

            }
            covs[[paste0("covariate", i, "_", covariates[i])]] <- list(covariate = covtemp[["covariate"]], type = covariateTypes[i], level = covtemp[["level"]],
                info = covtemp[["info"]], mapping = covtemp[["mapping"]])

        }  #i loop for each covariate

        ### total number of rasters
        totrast <- max(cnrastpersamp)

        ### number of support sites per sample (from choose.structure)
        nsups <- matrix(nrow = nsamps, ncol = nvarlevels)
        nsups <- as.data.frame(nsups)
        for (i in 1:nvarlevels) {
          colnames(nsups)[i] <- paste0("nsups.", i)
        }

        for (j in 1:nvarlevels) {
          if(structureObj[[j]]$recommendStruct!="noStructure"){
            for (i in 1:nsamps) {
              nsups[i, j] <- nrow(structureObj[[j]]$coordsU[structureObj[[j]]$coordsU$sample == nobs_df$sample[i], ])
            }
          }
          if(structureObj[[j]]$recommendStruct=="noStructure"){
            for(i in 1:nsamps){
              nsups[i,j]<-ifelse(nrow(data[data[[spatialVar]]==spatvarlevels[j] & data[[sampleVar]]==nobs_df$sample[i],])==0,0,1)
            }
          }
        }
        nobs_df <- cbind(nobs_df, nsups)

        # need to know for which samples there are no observations for a covariate
        nosup_var <- list()
        for (i in 1:nvarlevels) {
            nosup_var[[i]] <- nobs_df$sample[which(nsups[, i] == 0)]
        }
        gt0sup_var <- list()
        for (i in 1:nvarlevels) {
            gt0sup_var[[i]] <- nobs_df$sample[which(nsups[, i] > 0)]
        }

        ### cumulative raster vector per covariate level
        rastpervar <- matrix(nrow = nsamps * nvarlevels, ncol = 5)
        rastpervar <- as.data.frame(rastpervar)
        colnames(rastpervar) <- c("sample", "cov.level", "numraster", "start", "stop")
        rastpervar$sample <- rep(nobs_df$sample, each = nvarlevels)
        rastpervar$cov.level <- rep(spatvarlevels)
        for (i in 1:nrow(rastpervar)) {
            rastpervar$numraster[i] <- nrow(data[data[[sampleVar]] == rastpervar$sample[i] & data[[spatialVar]] == rastpervar$cov.level[i], ])
        }
        if (rastpervar$numraster[1] == 0) {
            rastpervar$start[1] <- 0
        }
        if (rastpervar$numraster[1] > 0) {
            rastpervar$start[1] <- 1
        }
        rastpervar$stop[1] <- rastpervar$numraster[1]
        for (i in 2:nrow(rastpervar)) {
            if (rastpervar$numraster[i] > 0) {
                rastpervar$start[i] <- rastpervar$stop[i - 1] + 1
                rastpervar$stop[i] <- rastpervar$start[i] + (rastpervar$numraster[i] - 1)
            }
            if (rastpervar$numraster[i] == 0) {
                rastpervar$start[i] <- rastpervar$stop[i - 1]
                rastpervar$stop[i] <- rastpervar$stop[i - 1]
            }
        }

        zeroranges<-rep(NA,nvarlevels)
        for(i in 1:nvarlevels){
          zeroranges[i]<-ifelse(structureObj[[i]]$recommendStruct=="noStructure",1,0)
        }

        if(sum(zeroranges)==length(zeroranges)){
          Kmat<-NULL
        }

        if(sum(zeroranges)<length(zeroranges)){
          Kmat.all <- list()
          for (z in 1:length(spatvarlevels)) {
            if(structureObj[[z]]$recommendStruct!="noStructure"){
              for (i in 1:length(gt0sup_var[[z]])) {
                ucoords <- structureObj[[z]]$coordsU[structureObj[[z]]$coordsU$sample == gt0sup_var[[z]][i], ]
                datacoords <- structureObj[[z]]$coordsData[structureObj[[z]]$coordsData$sample == gt0sup_var[[z]][i], ]

                dist.mat <- matrix(NA, nrow = nrow(datacoords), ncol = nrow(ucoords))
                Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]] <- matrix(0, nrow = nrow(datacoords), ncol = nrow(ucoords))
                for (j in 1:ncol(dist.mat)) {
                  for (k in 1:nrow(dist.mat)) {
                    dist.mat[k, j] <- distance(ucoords$x.omega[j], datacoords$xPC[k], ucoords$y.omega[j], datacoords$yPC[k])
                    Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]][k, j] <- dnorm(0, dist.mat[k, j], sd = structureObj[[z]]$recommendSd)
                  }
                }
              }
            }

            if(structureObj[[z]]$recommendStruct=="noStructure"){
              for (i in 1:length(gt0sup_var[[z]])){
                Kmat.all[[paste0("spatcovlevel", z, "_sample", gt0sup_var[[z]][i])]]<-matrix(NA,nrow=nrow(structureObj[[z]]$coordsData[structureObj[[z]]$coordsData$sample == gt0sup_var[[z]][i], ]), ncol=max(nsups))
              }

            }
          }

          rastpervar <- rastpervar[rastpervar$numraster > 0, ]

          Kmat <- as.data.frame(Kmat.all[[paste0("spatcovlevel", which(spatvarlevels == rastpervar$cov.level[1]), "_sample", rastpervar$sample[1])]])
          for (i in 2:nrow((rastpervar))) {
            Kmat <- rbind.fill(Kmat, as.data.frame(Kmat.all[[paste0("spatcovlevel", which(spatvarlevels == rastpervar$cov.level[i]), "_sample",
                                                                    rastpervar$sample[i])]]))
          }
        }

        nsubjs <- NULL
        cnsampspersub <- NULL
        outlist <- list(data = data, nSubjs = nsubjs, nSamps = nsamps, cNSampsPerSubj = cnsampspersub, cNRastPerSamp = cnrastpersamp, totalRasters = totrast,
            covs = covs, KMat = Kmat, rastersPerVar = rastpervar, nSupportSites = nsups, nObs = nobs_df, gT0SupportSites = gt0sup_var, nVarLevels = nvarlevels,
            subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, covariates = covariates, covariateTypes = covariateTypes, covariateLevels = covariateLevels,
            outcome = outcome, recStructures=recstructure)
        class(outlist) <- "PCDataList"
        return(outlist)

    }





    ## neither subject nor spatial variable ##
    if (is.null(subjectVar) == TRUE & is.null(spatialVar) == TRUE) {
        if (trimData == TRUE) {
            data <- data[, which(colnames(data) %in% unique(c(sampleVar, covariates, outcome, "xPC", "yPC")))]
        }

        data <- data[order(data[[sampleVar]]), ]
        samplelabs <- unique(data[[sampleVar]])
        nsamps <- length(samplelabs)
        nobs_df <- matrix(nrow = nsamps, ncol = 2)
        nobs_df <- as.data.frame(nobs_df)
        colnames(nobs_df) <- c("sample", "numrasters")
        nobs_df$sample <- samplelabs
        for (i in 1:nsamps) {
            nobs_df$numrasters[i] <- nrow(data[data[[sampleVar]] == nobs_df[i, 1], ])
        }

        ### cumulative raster vector
        cnrastpersamp <- c(0, cumsum(nobs_df$numrasters))

        ### covariates need separate matrix for each covariate secondary is on the sample level
        covs <- list()
        covtemp <- list()
        for (i in 1:length(covariates)) {
            ## warnings
            if (covariateTypes[i] == "continuous") {
                if (is.numeric(data[[covariates[i]]]) == FALSE)
                  stop(paste0("Covariate ", covariates[i], " is not continuous."))
            }
            if (covariateTypes[i] == "binary") {
                if (length(unique(data[[covariates[i]]])) != 2)
                  stop(paste0("Covariate ", covariates[i], " is not binary."))
            }

            if (covariateTypes[i] == "continuous") {
                covtemp[["covariate"]] <- covariates[i]
                covtemp[["level"]] <- covariateLevels[i]

                if (covariateLevels[i] == "sample") {
                  covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                  covtemp[["info"]][, 1] <- samplelabs
                  for (j in 1:nrow(covtemp[["info"]])) {
                    covtemp[["info"]][j, 2] <- unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]])
                  }
                }
                if (covariateLevels[i] == "raster") {
                  covtemp[["info"]] <- data[[covariates[i]]]
                }
                covs[[i]] <- covtemp

            }



            if (covariateTypes[i] %in% c("binary", "categorical")) {
                if (covariateTypes[i] == "binary") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 3)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("sample", "oldvalue", "newvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                      covtemp[["info"]][j, 3] <- ifelse(covtemp[["info"]][j, 2] == covvalues[1], 0, 1)
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 3)])
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["mapping"]] <- matrix(NA, nrow = length(covvalues), ncol = 2)
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- c("oldvalue", "newvalue")
                    covtemp[["mapping"]][, 1] <- covvalues
                    for (j in 1:nrow(covtemp[["mapping"]])) {
                      covtemp[["mapping"]][j, 2] <- ifelse(covtemp[["mapping"]][j, 1] == covvalues[1], 0, 1)
                    }

                    covtemp[["info"]] <- rep(NA, nrow(data))
                    for (j in 1:nrow(data)) {
                      covtemp[["info"]][j] <- ifelse(data[[covariates[i]]][j] == covvalues[1], 0, 1)
                    }
                  }
                }



                if (covariateTypes[i] == "categorical") {
                  covtemp[["covariate"]] <- covariates[i]
                  covtemp[["level"]] <- covariateLevels[i]

                  if (covariateLevels[i] == "sample") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- matrix(NA, nrow = nsamps, ncol = 2)
                    covtemp[["info"]] <- as.data.frame(covtemp[["info"]])
                    colnames(covtemp[["info"]]) <- c("subject", "oldvalue")
                    covtemp[["info"]][, 1] <- samplelabs
                    for (j in 1:nrow(covtemp[["info"]])) {
                      covtemp[["info"]][j, 2] <- as.character(unique(data[[covariates[i]]][data[[sampleVar]] == covtemp[["info"]][j, 1]]))
                    }
                    covtemp[["info"]] <- as.matrix(covtemp[["info"]][, c(1, 2)])

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(covtemp[["info"]]), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(covtemp[["info"]])) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j, 2] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                  if (covariateLevels[i] == "raster") {
                    covvalues <- sort(unique(data[[covariates[i]]]))

                    covtemp[["info"]] <- data[[covariates[i]]]

                    covtemp[["mapping"]] <- matrix(0, nrow = nrow(data), ncol = (length(covvalues) - 1))
                    covtemp[["mapping"]] <- as.data.frame(covtemp[["mapping"]])
                    colnames(covtemp[["mapping"]]) <- covvalues[2:length(covvalues)]
                    for (j in 1:nrow(data)) {
                      for (k in 1:ncol(covtemp[["mapping"]])) {
                        if (covtemp[["info"]][j] == covvalues[k + 1]) {
                          covtemp[["mapping"]][j, k] <- 1
                        }
                      }
                    }
                  }

                }


            }
            covs[[paste0("covariate", i, "_", covariates[i])]] <- list(covariate = covtemp[["covariate"]], type = covariateTypes[i], level = covtemp[["level"]],
                info = covtemp[["info"]], mapping = covtemp[["mapping"]])

        }  #i loop for each covariate

        ### total number of rasters
        totrast <- max(cnrastpersamp)

        if(structureObj[[1]]$recommendStruct=='noStructure'){
          nsups <- rep(0, nsamps)
          Kmat<-NULL
        }

        if(structureObj[[1]]$recommendStruct!='noStructure'){
          ### number of support sites per sample (from choose.structure)
          nsups <- rep(NA, nsamps)
          for (i in 1:nsamps) {
            nsups[i] <- nrow(structureObj[[1]]$coordsU[structureObj[[1]]$coordsU$sample == samplelabs[i], ])
          }
          nobs_df <- cbind(nobs_df, nsups)

          Kmat.all <- list()
          for (i in 1:nsamps) {
            ucoords <- structureObj[[1]]$coordsU[structureObj[[1]]$coordsU$sample == samplelabs[i], ]
            datacoords <- structureObj[[1]]$coordsData[structureObj[[1]]$coordsData$sample == samplelabs[i], ]

            dist.mat <- matrix(NA, nrow = nrow(datacoords), ncol = nrow(ucoords))
            Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]] <- matrix(0, nrow = nrow(datacoords), ncol = nrow(ucoords))
            for (j in 1:ncol(dist.mat)) {
              for (k in 1:nrow(dist.mat)) {
                dist.mat[k, j] <- distance(ucoords$x.omega[j], datacoords$xPC[k], ucoords$y.omega[j], datacoords$yPC[k])
                Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]][k, j] <- dnorm(0, dist.mat[k, j], sd = structureObj[[1]]$recommendSd)
              }
            }
          }


          Kmat <- as.data.frame(Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[1])]])
          for (i in 2:nsamps) {
            Kmat <- rbind.fill(Kmat, as.data.frame(Kmat.all[[paste0("spatcovlevel1_sample", samplelabs[i])]]))
          }
        }

        recstructure<-ifelse(structureObj[[1]]$recommendStruct=='noStructure',0,1)

        nsubjs <- NULL
        cnsampspersub <- NULL
        rastpervar <- NULL
        nvarlevels <- NULL
        gt0sup_var <- NULL
        outlist <- list(data = data, nSubjs = nsubjs, nSamps = nsamps, cNSampsPerSubj = cnsampspersub, cNRastPerSamp = cnrastpersamp, totalRasters = totrast,
            covs = covs, KMat = Kmat, rastersPerVar = rastpervar, nSupportSites = nsups, nObs = nobs_df, gT0SupportSites = gt0sup_var, nVarLevels = nvarlevels,
            subjectVar = subjectVar, sampleVar = sampleVar, spatialVar = spatialVar, covariates = covariates, covariateTypes = covariateTypes, covariateLevels = covariateLevels,
            outcome = outcome, recStructures=recstructure)
        class(outlist) <- "PCDataList"
        return(outlist)

    }

}










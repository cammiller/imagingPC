
#' Write the model
#'
#' The \code{writePCModel} function writes a NIMBLE model that incorporates the
#' process convolution approach.
#'
#' @param PCDataObj An object of class \code{PCDataList}.  This object is the
#' result of using the \code{\link{createPCData}} function.
#' @param typeOfZero A character string specifying the type of zero that
#' will be assumed.  The two options are "censored" (default) and "true".
#' If there are no zeros for the outcome of interest then this argument
#' is superfluous.
#' @param covariates A concatenated character string specifying the covariates
#' that will be included in the model.  When there are no zeros in the data,
#' or when zeros are assumed to be censored, there is only a single set of
#' covariates.  When the zeros are assumed to be true zeros, there are two
#' model parts, a binary and a marginalized.  The covariates specified in this
#' argument will be used in the marginalized part of the model.
#' @param covariatesForBinary A concatenated character string specifying
#' covariates that will be included in the binary part of the marginalized
#' two part model.  This argument will only be used if there are zeros in
#' the data and \code{typeOfZero='true'}.  If \code{typeOfZero='true'}, and
#' this argument is left empty, then the covariates for the binary part of
#' the model will be assumed to be the same as those for the marginalized part
#' of the model.
#' @param coefPrior A character string specifying the prior that should be
#' assumed for model coefficients.  The options are "sdunif" (default),
#' "dnorm", and "Cauchy".  See details for more information.
#' @param multiSampsPerSubj A TRUE/FALSE variable indicating whether or not
#' there are mutliple samples per subject.
#' @param errorVarianceLevel A character string specifying the level at which
#' the error variance should be estimated.  The options are "sample" (default),
#' "subject", and "overall".  As an example, if
#' \code{errorVarianceLevel="sample"} then a separate error variance is
#' estimated for each sample.
#' @param latentVarianceLevel A character string specifying the level at which
#' the latent variance should be estimated.  The options are "sample" (default),
#' "subject", and "overall".
#'
#' @return A list containing the model and other information
#' supplied to the \code{\link{createPCData}} function.
#'
#'
#'
#' \describe{
#'   \item{\code{model}}{A model object created using the
#'   \code{link[nimble]{nimbleCode}} function.}
#'   \item{\code{covariates}}{A concatenated character string of the
#'   covariate names.}
#'   \item{\code{covariatesForBinary}}{A concatenated character string of the
#'   covariate names for the binary part of the marginalized two-part model.}
#'   \item{\code{coefPrior}}{A character string specifying the assumed prior
#'   for the model coefficients.}
#'   \item{\code{multiSampsPerSubj}}{A TRUE/FALSE variable indicating whether
#'   there are multiple samples per subject.}
#'   \item{\code{errorVarianceLevel}}{A character string specifying the level
#'   (overall, subject, or sample) at which the error variance will be
#'   estimated.}
#'   \item{\code{latentVarianceLevel}}{A character string specifying the level
#'   (overall, subject, or sample) at which the latent variance will be
#'   estimated.}
#'   \item{\code{typeOfZero}}{A character string specifying the assumption
#'   made about the type of zeros in the data.}
#' }
#'
#'
#' @details For the model coefficients, three possible priors are allowed.
#'
#' \strong{sdunif}
#' For all model coefficients, including the intercept(s), the
#' prior is given by \deqn{\beta \sim \mbox{N}(0, \sigma^{2}_{\beta});
#' \sigma_{\beta} \sim \mbox{U}(0,10).}
#'
#' \strong{dnorm}
#' For all model coefficients, including the intercept(s), the
#' prior is given by \deqn{\beta \sim \mbox{N}(0, 0.00001).}
#'
#' \strong{Cauchy}
#' The Cauchy priors are specified according to the recommendations of Gelman
#' et al. (input citation).  While Gelman used the Cauchy priors in the
#' context of logistic regression, imaging mass spectrometry data is often
#' lognormally distributed, and on the log scale parameter spaces may need
#' to be constrained in a similar way.  For intercepts, the prior is
#' \deqn{\beta_{0} \sim \mbox{Cauchy}(0, \mbox{scale}=10).}
#' For other model coefficients the prior is
#' \deqn{\beta \sim \mbox{Cauchy}(0, \mbox{scale}=2.5).}
#'
#' @references de Valpine, P., D. Turek, C.J. Paciorek, C. Anderson-Bergman,
#' D. Temple Lang, and R. Bodik. 2017. Programming with models: writing
#' statistical algorithms for general model structures with NIMBLE.
#' \emph{Journal of Computational and Graphical Statistics}. 26: 403-413.
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



writePCModel <- function(PCDataObj, typeOfZero = "censored", covariates = NULL, covariatesForBinary = NULL, coefPrior = "sdunif", multiSampsPerSubj = NULL,
    errorVarianceLevel = "sample", latentVarianceLevel = "sample", addAssignNode = NULL, assignNodeLevel = NULL) {

    if (class(PCDataObj) != "PCDataList")
        stop("PCDataObj must of class PCDataList.  This object must be created using the createPCData function.")

    if ((coefPrior %in% c("sdunif", "dnorm", "dt")) == FALSE)
        stop("The only options for coefPrior are \"sdunif\", \"dnorm\", and \"Cauchy\".")
    if (any((errorVarianceLevel %in% c("subject", "sample", "overall")) == FALSE))
        stop("The only options for errorVarianceLevel are \"overall\", \"subject\", and \"sample\".")
    if (any((latentVarianceLevel %in% c("subject", "sample", "overall")) == FALSE))
        stop("The only options for latentVarianceLevel are \"overall\", \"subject\", and \"sample\".")
    if (any((typeOfZero %in% c("censored", "true")) == FALSE))
        stop("The only options for latentVarianceLevel are \"censored\" and \"true\".")
    if (multiSampsPerSubj == TRUE & is.null(PCDataObj[["nSubjs"]]) == TRUE)
        stop("There is no subject variable provided.  Make sure you specified the subject variable when you rescaled the data.")

    covs <- PCDataObj[["covs"]]

    # make sure data has been created
    if (is.null(covariates) == FALSE) {
        covsnames <- rep(NA, length(covs))
        for (i in 1:length(covs)) {
            covsnames[i] <- covs[[i]]$covariate
        }

        if (all(covariates %in% covsnames) == FALSE) {
            incovs <- rep(NA, length(covariates))
            for (i in 1:length(covariates)) {
                incovs[i] <- ifelse(covariates[i] %in% covsnames, 1, 0)
            }
            notincovs <- which(incovs == 0)
            stop(paste0("Data vectors needed for modeling have not been created for the following covariate(s): ", paste(covariates[notincovs],
                sep = "", collapse = ", "), ".  Use the createPCData function to make the necessary data vectors."))
        }
    }

    if (is.null(covariates) == FALSE) {
        covsnames <- rep(NA, length(covs))
        for (i in 1:length(covs)) {
            covsnames[i] <- covs[[i]]$covariate
        }

        if (all(covariates %in% covsnames) == FALSE) {
            incovs <- rep(NA, length(covariatesForBinary))
            for (i in 1:length(covariatesForBinary)) {
                incovs[i] <- ifelse(covariatesForBinary[i] %in% covsnames, 1, 0)
            }
            notincovs <- which(incovs == 0)
            stop(paste0("Data vectors needed for modeling have not been created for the following covariate(s) in covariatesForBinary: ", paste(covariatesForBinary[notincovs],
                sep = "", collapse = ", "), ".  Use the createPCData function to make the necessary data vectors."))
        }
    }

    data <- PCDataObj[["data"]]
    subjectVar <- PCDataObj[["subjectVar"]]
    sampleVar <- PCDataObj[["sampleVar"]]
    spatialVar <- PCDataObj[["spatialVar"]]
    rastpervar <- PCDataObj[["rastersPerVar"]]
    cNRastPerSamp <- PCDataObj[["cNRastPerSamp"]]
    recStructures<-PCDataObj[['recStructures']]
    covariate.types<-PCDataObj[['covariateTypes']]
    covariate.levels<-PCDataObj[['covariateLevels']]

    if (is.null(covariates) == TRUE) {
        covariates <- PCDataObj[["covariates"]]
    }

    covstokeep <- rep(NA, length(covs))
    for (i in 1:length(covs)) {
        covstokeep[i] <- ifelse(covs[[i]]$covariate %in% covariates, 1, 0)
    }
    covstokeep <- which(covstokeep == 1)
    covsb <- covs[covstokeep]

    ordcovnames <- rep(NA, length(covs))
    for (i in 1:length(covs)) {
        ordcovnames[i] <- covs[[i]]$covariate
    }
    covariates <- factor(covariates, levels = ordcovnames)
    covariates <- sort(covariates)
    covariates <- as.character(covariates)
    covariatesb <- covariates[covstokeep]
    covariate.typesb <- covariate.types[covstokeep]
    covariate.levelsb <- covariate.levels[covstokeep]

    outcome <- PCDataObj[["outcome"]]

    if (is.null(multiSampsPerSubj) == TRUE) {
        if (is.null(subjectVar) == TRUE) {
            multiSampsPerSubj <- FALSE
        }
        if (is.null(subjectVar) == FALSE) {
            multiSampsPerSubj <- TRUE
        }
    }

    data <- data[order(data[[sampleVar]]), ]
    samplelabs <- unique(data[[sampleVar]])
    nsamps <- length(samplelabs)
    if (is.null(spatialVar) == FALSE) {
        data <- data[order(data[[spatialVar]]), ]
        spatvarlevels <- unique(data[[spatialVar]])
        nvarlevels <- length(spatvarlevels)
    }

    if(is.null(spatialVar)==FALSE & any(recStructures==0)){
      for(i in 1:nvarlevels){
        assign(paste0('rastpervar',i), rastpervar[rastpervar$cov.level==spatvarlevels[i],])
      }
    }


    if (errorVarianceLevel == "sample") {
        errorvar.index <- "[j]"
    }
    if (errorVarianceLevel == "subject") {
        errorvar.index <- "[i]"
    }
    if (errorVarianceLevel == "overall") {
        errorvar.index <- ""
    }

    coef.indices <- matrix(nrow = 3, ncol = 2)
    coef.indices <- as.data.frame(coef.indices)
    colnames(coef.indices) <- c("level", "index")
    coef.indices$level <- c("subject", "sample", "raster")
    coef.indices$index <- c("i", "j", "k")


    nvars <- length(covsb)
    nbetasper <- rep(NA, nvars)
    for (i in 1:nvars) {
        if (covariate.typesb[i] %in% c("binary", "continuous")) {
            nbetasper[i] <- 1
        }
        if (covariate.typesb[i] == "categorical") {
            nbetasper[i] <- ncol(covsb[[i]]$mapping)
        }
    }
    nbetas <- sum(nbetasper)




    if (nrow(data[data[[outcome]] == 0, ]) == 0) {
        betanames <- list()
        for (i in 1:nvars) {
            if (covariate.typesb[i] %in% c("binary", "continuous")) {
                betanames[[i]] <- paste0("beta", i)
            }
            if (covariate.typesb[i] == "categorical") {
                betanames[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
            }
        }
        betanames <- paste(unlist(betanames))

        covarnames <- list()
        for (i in 1:nvars) {
            if (covariate.typesb[i] %in% c("binary", "continuous")) {
                covarnames[[i]] <- paste0("*", covariatesb[i], "[", coef.indices$index[coef.indices$level == covariate.levelsb[i]], "]")
            }
            if (covariate.typesb[i] == "categorical") {
                covarnames[[i]] <- paste0("*", covariatesb[i], "_", colnames(covsb[[i]]$mapping), "[", coef.indices$index[coef.indices$level ==
                  covariate.levelsb[i]], "]")
            }
        }
        covarnames <- paste(unlist(covarnames))

        modterms <- rep(NA, length(covarnames))
        for (i in 1:length(covarnames)) {
            modterms[i] <- paste0(betanames[i], covarnames[i])
        }
        modterms <- c("beta0", modterms)
        linpred <- paste(modterms, sep = "", collapse = " + ")


        if (multiSampsPerSubj == TRUE) {
            modbod <- rep(NA, 8)

            modbod[1] <- "for(i in 1:Subj){"
            modbod[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
            modbod[3] <- "for(k in (Raster[j]+1):Raster[j+1]){"
            modbod[4] <- paste("y[k]~dlnorm(mu[k], tau", errorvar.index, ")", sep = "", collapse = "")
            if(is.null(spatialVar)==TRUE & sum(recStructures)==0){
              modbod[5] <- paste("mu[k]<-", linpred, " + bi[i]", sep = "", collapse = "")
            } else {
              modbod[5] <- paste("mu[k]<-", linpred, " + bi[i] + PC[k]", sep = "", collapse = "")
            }
            modbod[6] <- "} #i (raster loop)"
            modbod[7] <- "} #j (sample loop)"
            modbod[8] <- "} #i (subject loop)"
        }

        if (multiSampsPerSubj == FALSE) {
            modbod <- rep(NA, 6)

            modbod[1] <- "for(j in 1:Sample){"
            modbod[2] <- "for(k in (Raster[j]+1):Raster[j+1]){"
            modbod[3] <- paste("y[k]~dlnorm(mu[k], tau", errorvar.index, ")", sep = "", collapse = "")
            if(is.null(spatialVar)==TRUE & recStructures==0){
              modbod[4] <- paste("mu[k]<-", linpred, sep = "", collapse = "")
            } else {
              modbod[4] <- paste("mu[k]<-", linpred, " + PC[k]", sep = "", collapse = "")
            }
            modbod[5] <- "} #i (raster loop)"
            modbod[6] <- "} #j (sample loop)"
        }


        # PC component for each sample and spatial covariate level
        if (is.null(rastpervar) == FALSE) {
            modsampPC <- rep(NA, nrow(rastpervar) * 5)
            for (i in 1:nrow(rastpervar)) {
                modsampPC[(i - 1) * 5 + 1] <- ""
                modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", rastpervar$sample[i], ", Spatial level=", rastpervar$cov.level[i], sep = "", collapse = "")
                modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", rastpervar$start[i], ":", rastpervar$stop[i], "){", sep = "", collapse = "")
                if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==0){
                  modsampPC[(i - 1) * 5 + 4]<-paste('PC[n]<-rk',which(spatvarlevels == rastpervar$cov.level[i]),'[n]', sep = "", collapse = "")
                }
                if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==1){
                  modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i],
                                                      ",1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]],Kmat[n,1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]])", sep = "", collapse = "")
                }
                modsampPC[(i - 1) * 5 + 5] <- "}"
            }
        }

        if (is.null(rastpervar) == TRUE) {
          if(recStructures==0){
            modsampPC<-''
          } else {
            modsampPC <- rep(NA, nsamps * 5)
            for (i in 1:nsamps) {
              modsampPC[(i - 1) * 5 + 1] <- ""
              modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", samplelabs[i], sep = "", collapse = "")
              modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", (cNRastPerSamp[i] + 1), ":", cNRastPerSamp[i + 1], "){", sep = "", collapse = "")
              modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", "[", samplelabs[i], ",1:nsup", "[", samplelabs[i], "]],Kmat[n,1:nsup", "[",
                                                  samplelabs[i], "]])", sep = "", collapse = "")
              modsampPC[(i - 1) * 5 + 5] <- "}"
            }
          }
        }


        # PC priors
        if (multiSampsPerSubj == TRUE & is.null(spatialVar) == TRUE & latentVarianceLevel == "sample") {
          if(recStructures==0){
            PCpriors<-''
          }
          if(recStructures==1){
            if(latentVarianceLevel=='sample'){
              PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)",
                            "sdxc[j]~dunif(0,10)", "}#j","}#i")
            }
            if(latentVarianceLevel=='subject'){
              PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[i])", "}#l",
                            "}#j","tauxc[i]<-pow(sdxc[i],-2)","sdxc[i]~dunif(0,10)", "}#i")
            }
            if(latentVarianceLevel=='overall'){
              PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l", "}#j","}#i","tauxc<-pow(sdxc,-2)",
                            "sdxc~dunif(0,10)")
            }

          }
        }
        if (multiSampsPerSubj == FALSE & is.null(spatialVar) == TRUE) {
          if(recStructures==0){
            PCpriors<-''
          }
          if(recStructures==1){
            if(latentVarianceLevel=='sample'){
              PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)", "sdxc[j]~dunif(0,10)",
                            "}#j")
            }
            if(latentVarianceLevel=='overall'){
              PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l","}#j", "tauxc<-pow(sdxc,-2)", "sdxc~dunif(0,10)","}#j")
            }

          }
        }

        if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "sample") {
          if(sum(recStructures)==length(recStructures)){
            PCpriors <- rep(NA, (2 + nvarlevels * 5))
            PCpriors[1] <- "for(i in 1:Subj){"
            PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

            for (i in 1:nvarlevels) {
              PCpriors[(i - 1) * 5 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 5] <- paste("}#l")
              PCpriors[(i - 1) * 5 + 6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
            }
            PCpriors <- c(PCpriors, "}#j", "}#i")
          }
          if(any(recStructures==0)){
            PCpriors<-''
            for(i in 1:nvarlevels){
              if(recStructures[i]==0){
                PCpriorstemp <- rep(NA, 6)
                PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                PCpriorstemp[4]<-'}#k'
                PCpriorstemp[5]<-'}#j'
                PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
              }
              if(recStructures[i]==1){
                PCpriorstemp <- rep(NA, 9)
                PCpriorstemp[1] <- "for(i in 1:Subj){"
                PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                 latentVarianceLevel], "]){", sep = "", collapse = "")
                PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                         rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                PCpriorstemp[5] <- paste("}#l")
                PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                      rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                PCpriorstemp[8]<-"}#j"
                PCpriorstemp[9]<-"}#i"
              }
              PCpriors<-c(PCpriors,PCpriorstemp)
            }
            PCpriors<-PCpriors[-1]
          }
        }

        if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "subject") {
          if(sum(recStructures)==length(recStructures)){
            PCpriors <- rep(NA, (2 + nvarlevels * 3))
            PCpriors[1] <- "for(i in 1:Subj){"
            PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

            for (i in 1:nvarlevels) {
              PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
            }
            PCpriors <- c(PCpriors, "}#j")

            PCpriors2 <- rep(NA, 2 * nvarlevels)
            for (i in 1:nvarlevels) {
              PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                     rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
              PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
            }
            PCpriors <- c(PCpriors, PCpriors2, "}#i")
          }
          if(any(recStructures==0)){
            PCpriors<-''
            for(i in 1:nvarlevels){
              if(recStructures[i]==0){
                PCpriorstemp <- rep(NA, 6)
                PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                PCpriorstemp[4]<-'}#k'
                PCpriorstemp[5]<-'}#j'
                PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
              }
              if(recStructures[i]==1){
                PCpriorstemp <- rep(NA, 9)
                PCpriorstemp[1] <- "for(i in 1:Subj){"
                PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                 latentVarianceLevel], "]){", sep = "", collapse = "")
                PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                         rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
                PCpriorstemp[5] <- "}#l"
                PCpriorstemp[6] <- "}#j"
                PCpriorstemp[7] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
                PCpriorstemp[8] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
                PCpriorstemp[9] <- "}#i"
              }
              PCpriors<-c(PCpriors,PCpriorstemp)
            }
            PCpriors<-PCpriors[-1]
          }
        }

        if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "overall") {
          if(sum(recStructures)==length(recStructures)){
            PCpriors <- rep(NA, (2 + nvarlevels * 3))
            PCpriors[1] <- "for(i in 1:Subj){"
            PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

            for (i in 1:nvarlevels) {
              PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), ")", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
            }
            PCpriors <- c(PCpriors, "}#j",'}#i')

            PCpriors2 <- rep(NA, 2 * nvarlevels)
            for (i in 1:nvarlevels) {
              PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                     rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
              PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
            }
            PCpriors <- c(PCpriors, PCpriors2)
          }
          if(any(recStructures==0)){
            PCpriors<-''
            for(i in 1:nvarlevels){
              if(recStructures[i]==0){
                PCpriorstemp <- rep(NA, 6)
                PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                PCpriorstemp[4]<-'}#k'
                PCpriorstemp[5]<-'}#j'
                PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
              }
              if(recStructures[i]==1){
                PCpriorstemp <- rep(NA, 9)
                PCpriorstemp[1] <- "for(i in 1:Subj){"
                PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                PCpriorstemp[5] <- "}#l"
                PCpriorstemp[6] <- "}#j"
                PCpriorstemp[7] <- "}#i"
                PCpriorstemp[8] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                            rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                PCpriorstemp[9] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
              }
              PCpriors<-c(PCpriors,PCpriorstemp)
            }
            PCpriors<-PCpriors[-1]
          }
        }

        if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='sample') {
          if(sum(recStructures)==length(recStructures)){
            PCpriors <- rep(NA, (1 + nvarlevels * 5))
            PCpriors[1] <- "for(j in 1:Sample){"
            for (i in 1:nvarlevels) {
              PCpriors[(i - 1) * 5 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 4] <- paste("}#l")
              PCpriors[(i - 1) * 5 + 5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
              PCpriors[(i - 1) * 5 + 6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
            }
            PCpriors <- c(PCpriors, "}#j")
          }
          if(any(recStructures==0)){
            PCpriors<-''
            for(i in 1:nvarlevels){
              if(recStructures[i]==0){
                PCpriorstemp <- rep(NA, 6)
                PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                PCpriorstemp[4]<-'}#k'
                PCpriorstemp[5]<-'}#j'
                PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
              }
              if(recStructures[i]==1){
                PCpriorstemp <- rep(NA, 7)
                PCpriorstemp[1] <- "for(j in 1:Sample){"
                PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                               rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                PCpriorstemp[4] <- paste("}#l")
                PCpriorstemp[5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                            rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                PCpriorstemp[6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                PCpriorstemp[7]<-"}#j"
              }
              PCpriors<-c(PCpriors,PCpriorstemp)
            }
            PCpriors<-PCpriors[-1]
          }
        }

        if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='overall') {
          if(sum(recStructures)==length(recStructures)){
            PCpriors <- rep(NA, (1 + nvarlevels * 5))
            PCpriors[1] <- "for(j in 1:Sample){"
            for (i in 1:nvarlevels) {
              PCpriors[(i - 1) * 3 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                          rastpervar$cov.level[i]), ")", sep = "", collapse = "")
              PCpriors[(i - 1) * 3 + 4] <- paste("}#l")
            }
            PCpriors<-c(PCpriors,"}#j")
            PCpriors2<-rep(NA, nvarlevels*2)
            PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
            PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
            PCpriors<-c(PCpriors,PCpriors2)
          }
          if(any(recStructures==0)){
            PCpriors<-''
            for(i in 1:nvarlevels){
              if(recStructures[i]==0){
                PCpriorstemp <- rep(NA, 6)
                PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                PCpriorstemp[4]<-'}#k'
                PCpriorstemp[5]<-'}#j'
                PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
              }
              if(recStructures[i]==1){
                PCpriorstemp <- rep(NA, 7)
                PCpriorstemp[1] <- "for(j in 1:Sample){"
                PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                PCpriorstemp[4] <- "}#l"
                PCpriorstemp[5]<- "}#j"
                PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                            rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                PCpriorstemp[7]<-"}#j"
              }
              PCpriors<-c(PCpriors,PCpriorstemp)
            }
            PCpriors<-PCpriors[-1]
          }
        }


        # error variance prior
        if (errorVarianceLevel == "sample" & multiSampsPerSubj == TRUE) {
            errlines <- c("for(i in 1:Subj){", "for(j in (Sample[i]+1):Sample[i+1]){", "tau[j]<-pow(sigma[j],-2)", "sigma[j]~dunif(0,10)", "}#j",
                "}#i")
        }
        if (errorVarianceLevel == "subject" & multiSampsPerSubj == TRUE) {
            errlines <- c("for(i in 1:Subj){", "tau[i]<-pow(sigma[i],-2)", "sigma[i]~dunif(0,10)", "}#i")
        }
        if (errorVarianceLevel == "sample" & multiSampsPerSubj == FALSE) {
            errlines <- c("for(j in 1:Sample){", "tau[j]<-pow(sigma[j],-2)", "sigma[j]~dunif(0,10)", "}#j")
        }
        if (errorVarianceLevel == "overall") {
            errlines <- c("sigma~dunif(0,10)")
        }
        if (errorVarianceLevel == "subject" & multiSampsPerSubj == FALSE)
            stop("If multiSampsPerSubj==FALSE then you cannot estimate the error variance on the subject level.
             If multiSampsPerSubj==FALSE is correct, then change the error variance to sample or overall.")



        # model coefficients
        if (coefPrior == "sdunif") {
            modcoef <- list()
            for (i in 1:nvars) {
                if (covariate.typesb[i] %in% c("binary", "continuous")) {
                  modcoef[[i]] <- rep(NA, 3)
                  modcoef[[i]][1] <- paste("beta", i, "~dnorm(0,taub", i, ")", sep = "", collapse = "")
                  modcoef[[i]][2] <- paste("taub", i, "<-pow(sdb", i, ",-2)", sep = "", collapse = "")
                  modcoef[[i]][3] <- paste("sdb", i, "~dunif(0,10)", sep = "", collapse = "")

                }
                if (covariate.typesb[i] == "categorical") {
                  modcoef[[i]] <- rep(NA, 3 * ncol(covsb[[i]]$mapping))
                  for (j in 1:ncol(covsb[[i]]$mapping)) {
                    modcoef[[i]][(j - 1) * 3 + 1] <- paste("beta", i, "_", j, "~dnorm(0,taub", i, "_", j, ")", sep = "", collapse = "")
                    modcoef[[i]][(j - 1) * 3 + 2] <- paste("taub", i, "_", j, "<-pow(sdb", i, "_", j, ",-2)", sep = "", collapse = "")
                    modcoef[[i]][(j - 1) * 3 + 3] <- paste("sdb", i, "_", j, "~dunif(0,10)", sep = "", collapse = "")
                  }

                  # unused code for same variance for each beta modcoef[[i]]<-rep(NA,(ncol(covsb[[i]]$mapping)+2)) for(j in 1:ncol(covsb[[i]]$mapping)){
                  # modcoef[[i]][j]<-paste('beta',i,'_',j,'~dnorm(0,taub',i,')', sep='', collapse='') }
                  # modcoef[[i]][ncol(covsb[[i]]$mapping)+1]<-paste('taub',i,'<-pow(sdb',i,',-2)', sep='', collapse = '')
                  # modcoef[[i]][ncol(covsb[[i]]$mapping)+2]<-paste('sdb',i,'~dunif(0,10)', sep='', collapse = '')
                }
            }
            modcoef <- c("beta0~dnorm(0,taub0)", "taub0<-pow(sdb0,-2)", "sdb0~dunif(0,10)", paste0(unlist(modcoef)))
        }

        if (coefPrior == "dnorm") {
            modcoef <- rep(NA, length(betanames))
            for (i in 1:length(betanames)) {
                modcoef[i] <- paste(betanames[i], "~dnorm(0,0.00001)", sep = "", collapse = "")
            }
            modcoef <- c("beta0~dnorm(0,0.00001)", modcoef)
        }

        if (coefPrior == "Cauchy") {
            modcoef <- rep(NA, length(betanames))
            for (i in 1:length(betanames)) {
                modcoef[i] <- paste(betanames[i], "~dt(0,0.16,1)", sep = "", collapse = "")
            }
            modcoef <- c("beta0~dt(0,0.01,1)", modcoef)
        }

        subjint <- c("for(i in 1:Subj){", "bi[i]~dnorm(0,taubi)", "}#i", "taubi<-pow(sdbi,-2)", "sdbi~dunif(0,10)")


        # additional nodes to compare categorical coefficients
        if ("categorical" %in% covariate.typesb) {
            addbetanodes <- list()
            nodetemp <- list()
            for (i in 1:nvars) {
                if (covariate.typesb[i] == "categorical") {
                  # to provide names, need to get mapping matrix, column of betas and column of corresponding name
                  betanames2 <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                  covarnames2 <- paste0(covariatesb[i], "_", colnames(covsb[[i]]$mapping))
                  namemap <- data.frame(betanames2, covarnames2)

                  nodetemp[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                  nodemat <- expand.grid(paste(nodetemp[[i]]), paste(nodetemp[[i]]))
                  nodemat <- nodemat[nodemat[, 1] != nodemat[, 2], ]

                  nodemat[, 3] <- 0
                  for (j in 1:nrow(nodemat)) {
                    if (as.numeric(strsplit(as.character(nodemat[j, 1]), "_")[[1]][2]) > as.numeric(strsplit(as.character(nodemat[j, 2]), "_")[[1]][2])) {
                      nodemat[j, 3] <- 1
                    }
                  }
                  nodemat <- nodemat[nodemat[, 3] == 1, ]
                  nodemat <- nodemat[order(nodemat[, 1]), ]
                  nodemat <- nodemat[, c(1, 2)]

                  nodemat <- as.data.frame(nodemat)
                  colnames(nodemat) <- c("beta1", "beta2")
                  nodemat$cov1 <- NA
                  nodemat$cov2 <- NA
                  for (j in 1:nrow(nodemat)) {
                    nodemat$cov1[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta1[j]])
                    nodemat$cov2[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta2[j]])
                  }

                  # need to subtract columns
                  addbetanodes[[i]] <- rep(NA, nrow(nodemat))
                  for (j in 1:nrow(nodemat)) {
                    addbetanodes[[i]][j] <- paste0(nodemat$cov1[j], "_minus_", nodemat$cov2[j], "<-", nodemat$beta1[j], "-", nodemat$beta2[j])
                  }
                }
            }
            addbetanodes <- paste(unlist(addbetanodes))
        }



        if (multiSampsPerSubj == TRUE) {
            if (any(covariate.typesb == "categorical") == TRUE) {
                assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, addbetanodes, subjint, "})"))))
            }
            if (any(covariate.typesb == "categorical") == FALSE) {
                assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, subjint, "})"))))
            }
        }
        if (multiSampsPerSubj == FALSE) {
            if (any(covariate.typesb == "categorical") == TRUE) {
                assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, addbetanodes, "})"))))
            }
            if (any(covariate.typesb == "categorical") == FALSE) {
                assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, "})"))))
            }
        }

    }




    if (nrow(data[data[[outcome]] == 0, ]) > 0) {
        if (typeOfZero == "censored")
            {
                betanames <- list()
                for (i in 1:nvars) {
                  if (covariate.typesb[i] %in% c("binary", "continuous")) {
                    betanames[[i]] <- paste0("beta", i)
                  }
                  if (covariate.typesb[i] == "categorical") {
                    betanames[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                  }
                }
                betanames <- paste(unlist(betanames))

                covarnames <- list()
                for (i in 1:nvars) {
                  if (covariate.typesb[i] %in% c("binary", "continuous")) {
                    covarnames[[i]] <- paste0("*", covariatesb[i], "[", coef.indices$index[coef.indices$level == covariate.levelsb[i]], "]")
                  }
                  if (covariate.typesb[i] == "categorical") {
                    covarnames[[i]] <- paste0("*", covariatesb[i], "_", colnames(covsb[[i]]$mapping), "[", coef.indices$index[coef.indices$level ==
                      covariate.levelsb[i]], "]")
                  }
                }
                covarnames <- paste(unlist(covarnames))

                modterms <- rep(NA, length(covarnames))
                for (i in 1:length(covarnames)) {
                  modterms[i] <- paste0(betanames[i], covarnames[i])
                }
                modterms <- c("beta0", modterms)
                linpred <- paste(modterms, sep = "", collapse = " + ")

                # body of the model
                if (multiSampsPerSubj == TRUE) {
                  modbod <- rep(NA, 10)

                  modbod[1] <- "for(i in 1:Subj){"
                  modbod[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                  modbod[3] <- "for(k in (Raster[j]+1):Raster[j+1]){"
                  modbod[4] <- "zeros[k]~dpois(phi[k])"
                  modbod[5] <- "phi[k]<- -loglik[k] + C "
                  modbod[6] <- paste("mu[k]<-", linpred, " + bi[i] + PC[k]", sep = "", collapse = "")
                  modbod[7] <- paste("loglik[k]<-llcal.censored(y=y[k], obs=obs[k], sigma=sigma", errorvar.index, ", mu=mu[k])", sep = "", collapse = "")
                  modbod[8] <- "} #i (raster loop)"
                  modbod[9] <- "} #j (sample loop)"
                  modbod[10] <- "} #i (subject loop)"
                }

                if (multiSampsPerSubj == FALSE) {
                  modbod <- rep(NA, 8)

                  modbod[1] <- "for(j in 1:Sample){"
                  modbod[2] <- "for(k in (Raster[j]+1):Raster[j+1]){"
                  modbod[3] <- "zeros[k]~dpois(phi[k])"
                  modbod[4] <- "phi[k]<- -loglik[k] + C "
                  modbod[5] <- paste("mu[k]<-", linpred, " + PC[k]", sep = "", collapse = "")
                  modbod[6] <- paste("loglik[k]<-llcal.censored(y=y[k], obs=obs[k], sigma=sigma", errorvar.index, ", mu=mu[k])", sep = "", collapse = "")
                  modbod[7] <- "} #k (raster loop)"
                  modbod[8] <- "} #j (sample loop)"
                }


                # PC component for each sample and spatial covariate level
                if (is.null(rastpervar) == FALSE) {
                  modsampPC <- rep(NA, nrow(rastpervar) * 5)
                  for (i in 1:nrow(rastpervar)) {
                    modsampPC[(i - 1) * 5 + 1] <- ""
                    modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", rastpervar$sample[i], ", Spatial level=", rastpervar$cov.level[i], sep = "", collapse = "")
                    modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", rastpervar$start[i], ":", rastpervar$stop[i], "){", sep = "", collapse = "")
                    if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==0){
                      modsampPC[(i - 1) * 5 + 4]<-paste('PC[n]<-rk',which(spatvarlevels == rastpervar$cov.level[i]),'[n]', sep = "", collapse = "")
                    }
                    if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==1){
                      modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i],
                                                          ",1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]],Kmat[n,1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]])", sep = "", collapse = "")
                    }
                    modsampPC[(i - 1) * 5 + 5] <- "}"
                  }
                }

                if (is.null(rastpervar) == TRUE) {
                  if(recStructures==0){
                    modsampPC<-''
                  } else {
                    modsampPC <- rep(NA, nsamps * 5)
                    for (i in 1:nsamps) {
                      modsampPC[(i - 1) * 5 + 1] <- ""
                      modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", samplelabs[i], sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", (cNRastPerSamp[i] + 1), ":", cNRastPerSamp[i + 1], "){", sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", "[", samplelabs[i], ",1:nsup", "[", samplelabs[i], "]],Kmat[n,1:nsup", "[",
                                                          samplelabs[i], "]])", sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 5] <- "}"
                    }
                  }
                }



                # PC priors
                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == TRUE & latentVarianceLevel == "sample") {
                  if(recStructures==0){
                    PCpriors<-''
                  }
                  if(recStructures==1){
                    if(latentVarianceLevel=='sample'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)",
                                    "sdxc[j]~dunif(0,10)", "}#j","}#i")
                    }
                    if(latentVarianceLevel=='subject'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[i])", "}#l",
                                    "}#j","tauxc[i]<-pow(sdxc[i],-2)","sdxc[i]~dunif(0,10)", "}#i")
                    }
                    if(latentVarianceLevel=='overall'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l", "}#j","}#i","tauxc<-pow(sdxc,-2)",
                                    "sdxc~dunif(0,10)")
                    }

                  }
                }
                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == TRUE) {
                  if(recStructures==0){
                    PCpriors<-''
                  }
                  if(recStructures==1){
                    if(latentVarianceLevel=='sample'){
                      PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)", "sdxc[j]~dunif(0,10)",
                                    "}#j")
                    }
                    if(latentVarianceLevel=='overall'){
                      PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l","}#j", "tauxc<-pow(sdxc,-2)", "sdxc~dunif(0,10)","}#j")
                    }

                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "sample") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 5))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 5 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 5] <- paste("}#l")
                      PCpriors[(i - 1) * 5 + 6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                            rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, "}#j", "}#i")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                        PCpriorstemp[5] <- paste("}#l")
                        PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                        PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[8]<-"}#j"
                        PCpriorstemp[9]<-"}#i"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "subject") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 3))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
                    }
                    PCpriors <- c(PCpriors, "}#j")

                    PCpriors2 <- rep(NA, 2 * nvarlevels)
                    for (i in 1:nvarlevels) {
                      PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                             rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
                      PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, PCpriors2, "}#i")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
                        PCpriorstemp[5] <- "}#l"
                        PCpriorstemp[6] <- "}#j"
                        PCpriorstemp[7] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
                        PCpriorstemp[8] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[9] <- "}#i"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "overall") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 3))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
                    }
                    PCpriors <- c(PCpriors, "}#j",'}#i')

                    PCpriors2 <- rep(NA, 2 * nvarlevels)
                    for (i in 1:nvarlevels) {
                      PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                          rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                      PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, PCpriors2)
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                        PCpriorstemp[5] <- "}#l"
                        PCpriorstemp[6] <- "}#j"
                        PCpriorstemp[7] <- "}#i"
                        PCpriorstemp[8] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                 rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                        PCpriorstemp[9] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='sample') {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (1 + nvarlevels * 5))
                    PCpriors[1] <- "for(j in 1:Sample){"
                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 5 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 4] <- paste("}#l")
                      PCpriors[(i - 1) * 5 + 5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                            rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, "}#j")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 7)
                        PCpriorstemp[1] <- "for(j in 1:Sample){"
                        PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("}#l")
                        PCpriorstemp[5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                        PCpriorstemp[6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[7]<-"}#j"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='overall') {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (1 + nvarlevels * 5))
                    PCpriors[1] <- "for(j in 1:Sample){"
                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("}#l")
                    }
                    PCpriors<-c(PCpriors,"}#j")
                    PCpriors2<-rep(NA, nvarlevels*2)
                    PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                        rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                    PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                    PCpriors<-c(PCpriors,PCpriors2)
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 7)
                        PCpriorstemp[1] <- "for(j in 1:Sample){"
                        PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                        PCpriorstemp[4] <- "}#l"
                        PCpriorstemp[5]<- "}#j"
                        PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                 rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                        PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[7]<-"}#j"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }


                # error variance prior
                if (errorVarianceLevel == "sample" & multiSampsPerSubj == TRUE) {
                  errlines <- c("for(i in 1:Subj){", "for(j in (Sample[i]+1):Sample[i+1]){", "sigma[j]~dunif(0,10)", "}#j", "}#i")
                }
                if (errorVarianceLevel == "sample" & multiSampsPerSubj == FALSE) {
                  errlines <- c("for(j in 1:Sample){", "sigma[j]~dunif(0,10)", "}#j")
                }
                if (errorVarianceLevel == "subject" & multiSampsPerSubj == TRUE) {
                  errlines <- c("for(i in 1:Subj){", "sigma[i]~dunif(0,10)", "}#i")
                }
                if (errorVarianceLevel == "overall") {
                  errlines <- c("sigma~dunif(0,10)")
                }
                if (errorVarianceLevel == "subject" & multiSampsPerSubj == FALSE)
                  stop("If multiSampsPerSubj==FALSE then you cannot estimate the error variance on the subject level.
             If multiSampsPerSubj==FALSE is correct, then change the error variance to sample or overall.")



                # model coefficients
                if (coefPrior == "sdunif") {
                  modcoef <- list()
                  for (i in 1:nvars) {
                    if (covariate.typesb[i] %in% c("binary", "continuous")) {
                      modcoef[[i]] <- rep(NA, 3)
                      modcoef[[i]][1] <- paste("beta", i, "~dnorm(0,taub", i, ")", sep = "", collapse = "")
                      modcoef[[i]][2] <- paste("taub", i, "<-pow(sdb", i, ",-2)", sep = "", collapse = "")
                      modcoef[[i]][3] <- paste("sdb", i, "~dunif(0,10)", sep = "", collapse = "")

                    }
                    if (covariate.typesb[i] == "categorical") {
                      modcoef[[i]] <- rep(NA, 3 * ncol(covsb[[i]]$mapping))
                      for (j in 1:ncol(covsb[[i]]$mapping)) {
                        modcoef[[i]][(j - 1) * 3 + 1] <- paste("beta", i, "_", j, "~dnorm(0,taub", i, "_", j, ")", sep = "", collapse = "")
                        modcoef[[i]][(j - 1) * 3 + 2] <- paste("taub", i, "_", j, "<-pow(sdb", i, "_", j, ",-2)", sep = "", collapse = "")
                        modcoef[[i]][(j - 1) * 3 + 3] <- paste("sdb", i, "_", j, "~dunif(0,10)", sep = "", collapse = "")
                      }

                      # unused code for same variance for each beta modcoef[[i]]<-rep(NA,(ncol(covsb[[i]]$mapping)+2)) for(j in 1:ncol(covsb[[i]]$mapping)){
                      # modcoef[[i]][j]<-paste('beta',i,'_',j,'~dnorm(0,taub',i,')', sep='', collapse='') }
                      # modcoef[[i]][ncol(covsb[[i]]$mapping)+1]<-paste('taub',i,'<-pow(sdb',i,',-2)', sep='', collapse = '')
                      # modcoef[[i]][ncol(covsb[[i]]$mapping)+2]<-paste('sdb',i,'~dunif(0,10)', sep='', collapse = '')
                    }
                  }
                  modcoef <- c("beta0~dnorm(0,taub0)", "taub0<-pow(sdb0,-2)", "sdb0~dunif(0,10)", paste0(unlist(modcoef)))
                }

                if (coefPrior == "dnorm") {
                  modcoef <- rep(NA, length(betanames))
                  for (i in 1:length(betanames)) {
                    modcoef[i] <- paste(betanames[i], "~dnorm(0,0.00001)", sep = "", collapse = "")
                  }
                  modcoef <- c("beta0~dnorm(0,0.00001)", modcoef)
                }

                if (coefPrior == "Cauchy") {
                  modcoef <- rep(NA, length(betanames))
                  for (i in 1:length(betanames)) {
                    modcoef[i] <- paste(betanames[i], "~dt(0,0.16,1)", sep = "", collapse = "")
                  }
                  modcoef <- c("beta0~dt(0,0.01,1)", modcoef)
                }

                subjint <- c("for(i in 1:Subj){", "bi[i]~dnorm(0,taubi)", "}#i", "taubi<-pow(sdbi,-2)", "sdbi~dunif(0,10)")


                # additional nodes to compare categorical coefficients
                if ("categorical" %in% covariate.typesb) {
                  addbetanodes <- list()
                  nodetemp <- list()
                  for (i in 1:nvars) {
                    if (covariate.typesb[i] == "categorical") {
                      # to provide names, need to get mapping matrix, column of betas and column of corresponding name
                      betanames2 <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      covarnames2 <- paste0(covariatesb[i], "_", colnames(covsb[[i]]$mapping))
                      namemap <- data.frame(betanames2, covarnames2)

                      nodetemp[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      nodemat <- expand.grid(paste(nodetemp[[i]]), paste(nodetemp[[i]]))
                      nodemat <- nodemat[nodemat[, 1] != nodemat[, 2], ]

                      nodemat[, 3] <- 0
                      for (j in 1:nrow(nodemat)) {
                        if (as.numeric(strsplit(as.character(nodemat[j, 1]), "_")[[1]][2]) > as.numeric(strsplit(as.character(nodemat[j, 2]), "_")[[1]][2])) {
                          nodemat[j, 3] <- 1
                        }
                      }
                      nodemat <- nodemat[nodemat[, 3] == 1, ]
                      nodemat <- nodemat[order(nodemat[, 1]), ]
                      nodemat <- nodemat[, c(1, 2)]

                      nodemat <- as.data.frame(nodemat)
                      colnames(nodemat) <- c("beta1", "beta2")
                      nodemat$cov1 <- NA
                      nodemat$cov2 <- NA
                      for (j in 1:nrow(nodemat)) {
                        nodemat$cov1[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta1[j]])
                        nodemat$cov2[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta2[j]])
                      }

                      # need to subtract columns
                      addbetanodes[[i]] <- rep(NA, nrow(nodemat))
                      for (j in 1:nrow(nodemat)) {
                        addbetanodes[[i]][j] <- paste0(nodemat$cov1[j], "_minus_", nodemat$cov2[j], "<-", nodemat$beta1[j], "-", nodemat$beta2[j])
                      }
                    }
                  }
                  addbetanodes <- paste(unlist(addbetanodes))
                }



                if (multiSampsPerSubj == TRUE) {
                  if (any(covariate.typesb == "categorical") == TRUE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, addbetanodes, subjint, "})"))))
                  }
                  if (any(covariate.typesb == "categorical") == FALSE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, subjint, "})"))))
                  }
                }
                if (multiSampsPerSubj == FALSE) {
                  if (any(covariate.typesb == "categorical") == TRUE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, addbetanodes, modcoef, "})"))))
                  }
                  if (any(covariate.typesb == "categorical") == FALSE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, modcoef, "})"))))
                  }
                }

            }  # typeOfZero=='censored'





        if (typeOfZero == "true")
            {

                # covariates and coefficients for marginalized model part
                betanames <- list()
                for (i in 1:nvars) {
                  if (covariate.typesb[i] %in% c("binary", "continuous")) {
                    betanames[[i]] <- paste0("beta", i)
                  }
                  if (covariate.typesb[i] == "categorical") {
                    betanames[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                  }
                }
                betanames <- paste(unlist(betanames))

                covarnames <- list()
                for (i in 1:nvars) {
                  if (covariate.typesb[i] %in% c("binary", "continuous")) {
                    covarnames[[i]] <- paste0("*", covariatesb[i], "[", coef.indices$index[coef.indices$level == covariate.levelsb[i]], "]")
                  }
                  if (covariate.typesb[i] == "categorical") {
                    covarnames[[i]] <- paste0("*", covariatesb[i], "_", colnames(covsb[[i]]$mapping), "[", coef.indices$index[coef.indices$level ==
                      covariate.levelsb[i]], "]")
                  }
                }
                covarnames <- paste(unlist(covarnames))

                modterms <- rep(NA, length(covarnames))
                for (i in 1:length(covarnames)) {
                  modterms[i] <- paste0(betanames[i], covarnames[i])
                }
                modterms <- c("beta0", modterms)
                linpred <- paste(modterms, sep = "", collapse = " + ")



                # covariates and coefficients for binary model part
                if (is.null(covariatesForBinary) == TRUE) {
                  covsa <- covsb
                  covariatesForBinary <- covariatesb
                }
                if (is.null(covariatesForBinary) == FALSE) {
                  covstokeep <- rep(NA, length(covs))
                  for (i in 1:length(covs)) {
                    covstokeep[i] <- ifelse(covs[[i]]$covariate %in% covariatesForBinary, 1, 0)
                  }
                  covstokeep <- which(covstokeep == 1)
                  covsa <- covs[covstokeep]
                }
                # need to get correct order of covariatesa
                covariatesForBinary <- factor(covariatesForBinary, levels = ordcovnames)
                covariatesForBinary <- sort(covariatesForBinary)
                covariatesForBinary <- as.character(covariatesForBinary)
                covariatesa <- covariatesForBinary[covstokeep]
                covariate.typesa <- covariate.types[covstokeep]
                covariate.levelsa <- covariate.levels[covstokeep]
                nvarsa <- length(covsa)

                alphanames <- list()
                for (i in 1:nvarsa) {
                  if (covariate.typesa[i] %in% c("binary", "continuous")) {
                    alphanames[[i]] <- paste0("alpha", i)
                  }
                  if (covariate.typesa[i] == "categorical") {
                    alphanames[[i]] <- paste0("alpha", i, "_", seq(1:ncol(covsa[[i]]$mapping)))
                  }
                }
                alphanames <- paste(unlist(alphanames))

                covarnamesa <- list()
                for (i in 1:nvarsa) {
                  if (covariate.typesa[i] %in% c("binary", "continuous")) {
                    covarnamesa[[i]] <- paste0("*", covariatesa[i], "[", coef.indices$index[coef.indices$level == covariate.levelsa[i]], "]")
                  }
                  if (covariate.typesa[i] == "categorical") {
                    covarnamesa[[i]] <- paste0("*", covariatesa[i], "_", colnames(covsa[[i]]$mapping), "[", coef.indices$index[coef.indices$level ==
                      covariate.levelsa[i]], "]")
                  }
                }
                covarnamesa <- paste(unlist(covarnamesa))

                modtermsa <- rep(NA, length(covarnamesa))
                for (i in 1:length(covarnamesa)) {
                  modtermsa[i] <- paste0(alphanames[i], covarnamesa[i])
                }
                modtermsa <- c("alpha0", modtermsa)
                linpreda <- paste(modtermsa, sep = "", collapse = " + ")



                # body of the model
                if (multiSampsPerSubj == TRUE) {
                  modbod <- rep(NA, 12)

                  modbod[1] <- "for(i in 1:Subj){"
                  modbod[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                  modbod[3] <- "for(k in (Raster[j]+1):Raster[j+1]){"
                  modbod[4] <- "zeros[k]~dpois(phi[k])"
                  modbod[5] <- "phi[k]<- -loglik[k] + C "
                  modbod[6] <- paste("linbin[k]<-", linpreda, " + bi1[i] + v[j]", sep = "", collapse = "")
                  modbod[7] <- "binprob[k]<-ilogit(linbin[k])"
                  modbod[8] <- paste("mu[k]<-", linpred, " + bi2[i] + PC[k] - log(binprob[k]) - pow(sigma", errorvar.index, ",2)/2", sep = "", collapse = "")
                  modbod[9] <- paste("loglik[k]<-llcal.true(y=y[k], binprob=binprob[k], sigma=sigma", errorvar.index, ", mu=mu[k])", sep = "", collapse = "")
                  modbod[10] <- "} #i (raster loop)"
                  modbod[11] <- "} #j (sample loop)"
                  modbod[12] <- "} #i (subject loop)"
                }

                if (multiSampsPerSubj == FALSE) {
                  modbod <- rep(NA, 10)

                  modbod[1] <- "for(j in 1:Sample){"
                  modbod[2] <- "for(k in (Raster[j]+1):Raster[j+1]){"
                  modbod[3] <- "zeros[k]~dpois(phi[k])"
                  modbod[4] <- "phi[k]<- -loglik[k] + C "
                  modbod[5] <- paste("linbin[k]<-", linpred, " + v[j]", sep = "", collapse = "")
                  modbod[6] <- "binprob[k]<-ilogit(linbin[k])"
                  modbod[7] <- paste("mu[k]<-", linpred, " + PC[k] - log(binprob[k]) - pow(sigma", errorvar.index, ",2)/2", sep = "", collapse = "")
                  modbod[8] <- paste("loglik[k]<-llcal.true(y=y[k], binprob=binprob[k], sigma=sigma", errorvar.index, ", mu=mu[k])", sep = "", collapse = "")
                  modbod[9] <- "} #k (raster loop)"
                  modbod[10] <- "} #j (sample loop)"
                }


                # PC component for each sample and spatial covariate level
                if (is.null(rastpervar) == FALSE) {
                  modsampPC <- rep(NA, nrow(rastpervar) * 5)
                  for (i in 1:nrow(rastpervar)) {
                    modsampPC[(i - 1) * 5 + 1] <- ""
                    modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", rastpervar$sample[i], ", Spatial level=", rastpervar$cov.level[i], sep = "", collapse = "")
                    modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", rastpervar$start[i], ":", rastpervar$stop[i], "){", sep = "", collapse = "")
                    if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==0){
                      modsampPC[(i - 1) * 5 + 4]<-paste('PC[n]<-rk',which(spatvarlevels == rastpervar$cov.level[i]),'[n]', sep = "", collapse = "")
                    }
                    if(recStructures[which(spatvarlevels == rastpervar$cov.level[i])]==1){
                      modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i],
                                                          ",1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]],Kmat[n,1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", rastpervar$sample[i], "]])", sep = "", collapse = "")
                    }
                    modsampPC[(i - 1) * 5 + 5] <- "}"
                  }
                }

                if (is.null(rastpervar) == TRUE) {
                  if(recStructures==0){
                    modsampPC<-''
                  } else {
                    modsampPC <- rep(NA, nsamps * 5)
                    for (i in 1:nsamps) {
                      modsampPC[(i - 1) * 5 + 1] <- ""
                      modsampPC[(i - 1) * 5 + 2] <- paste("# Sample ", samplelabs[i], sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 3] <- paste("for(n in ", (cNRastPerSamp[i] + 1), ":", cNRastPerSamp[i + 1], "){", sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 4] <- paste("PC[n]<-inprod(xc", "[", samplelabs[i], ",1:nsup", "[", samplelabs[i], "]],Kmat[n,1:nsup", "[",
                                                          samplelabs[i], "]])", sep = "", collapse = "")
                      modsampPC[(i - 1) * 5 + 5] <- "}"
                    }
                  }
                }


                # PC priors
                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == TRUE & latentVarianceLevel == "sample") {
                  if(recStructures==0){
                    PCpriors<-''
                  }
                  if(recStructures==1){
                    if(latentVarianceLevel=='sample'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)",
                                    "sdxc[j]~dunif(0,10)", "}#j","}#i")
                    }
                    if(latentVarianceLevel=='subject'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[i])", "}#l",
                                    "}#j","tauxc[i]<-pow(sdxc[i],-2)","sdxc[i]~dunif(0,10)", "}#i")
                    }
                    if(latentVarianceLevel=='overall'){
                      PCpriors <- c("for(i in 1:Subj){","for(j in (Sample[i]+1):Sample[i+1]){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l", "}#j","}#i","tauxc<-pow(sdxc,-2)",
                                    "sdxc~dunif(0,10)")
                    }

                  }
                }
                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == TRUE) {
                  if(recStructures==0){
                    PCpriors<-''
                  }
                  if(recStructures==1){
                    if(latentVarianceLevel=='sample'){
                      PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc[j])", "}#l", "tauxc[j]<-pow(sdxc[j],-2)", "sdxc[j]~dunif(0,10)",
                                    "}#j")
                    }
                    if(latentVarianceLevel=='overall'){
                      PCpriors <- c("for(j in 1:Sample){", "for(l in 1:nsup[j]){", "xc[j,l]~dnorm(0,tauxc)", "}#l","}#j", "tauxc<-pow(sdxc,-2)", "sdxc~dunif(0,10)","}#j")
                    }

                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "sample") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 5))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 5 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 5] <- paste("}#l")
                      PCpriors[(i - 1) * 5 + 6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                            rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, "}#j", "}#i")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                        PCpriorstemp[5] <- paste("}#l")
                        PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                        PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[8]<-"}#j"
                        PCpriorstemp[9]<-"}#i"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "subject") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 3))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
                    }
                    PCpriors <- c(PCpriors, "}#j")

                    PCpriors2 <- rep(NA, 2 * nvarlevels)
                    for (i in 1:nvarlevels) {
                      PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                             rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
                      PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, PCpriors2, "}#i")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[i])", sep = "", collapse = "")
                        PCpriorstemp[5] <- "}#l"
                        PCpriorstemp[6] <- "}#j"
                        PCpriorstemp[7] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[i],-2)", sep = "", collapse = "")
                        PCpriorstemp[8] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[i]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[9] <- "}#i"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == TRUE & is.null(spatialVar) == FALSE & latentVarianceLevel == "overall") {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (2 + nvarlevels * 3))
                    PCpriors[1] <- "for(i in 1:Subj){"
                    PCpriors[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"

                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 5] <- paste("}#l")
                    }
                    PCpriors <- c(PCpriors, "}#j",'}#i')

                    PCpriors2 <- rep(NA, 2 * nvarlevels)
                    for (i in 1:nvarlevels) {
                      PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                          rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                      PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, PCpriors2)
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 9)
                        PCpriorstemp[1] <- "for(i in 1:Subj){"
                        PCpriorstemp[2] <- "for(j in (Sample[i]+1):Sample[i+1]){"
                        PCpriorstemp[3] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                        PCpriorstemp[5] <- "}#l"
                        PCpriorstemp[6] <- "}#j"
                        PCpriorstemp[7] <- "}#i"
                        PCpriorstemp[8] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                 rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                        PCpriorstemp[9] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='sample') {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (1 + nvarlevels * 5))
                    PCpriors[1] <- "for(j in 1:Sample){"
                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 5 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 4] <- paste("}#l")
                      PCpriors[(i - 1) * 5 + 5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                            rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                      PCpriors[(i - 1) * 5 + 6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                    }
                    PCpriors <- c(PCpriors, "}#j")
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 7)
                        PCpriorstemp[1] <- "for(j in 1:Sample){"
                        PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), "[j])", sep = "", collapse = "")
                        PCpriorstemp[4] <- paste("}#l")
                        PCpriorstemp[5] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                    rastpervar$cov.level[i]), "[j],-2)", sep = "", collapse = "")
                        PCpriorstemp[6] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "[j]~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[7]<-"}#j"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }

                if (multiSampsPerSubj == FALSE & is.null(spatialVar) == FALSE & latentVarianceLevel=='overall') {
                  if(sum(recStructures)==length(recStructures)){
                    PCpriors <- rep(NA, (1 + nvarlevels * 5))
                    PCpriors[1] <- "for(j in 1:Sample){"
                    for (i in 1:nvarlevels) {
                      PCpriors[(i - 1) * 3 + 2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                                       latentVarianceLevel], "]){", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                               rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                      PCpriors[(i - 1) * 3 + 4] <- paste("}#l")
                    }
                    PCpriors<-c(PCpriors,"}#j")
                    PCpriors2<-rep(NA, nvarlevels*2)
                    PCpriors2[(i - 1) * 2 + 1] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                        rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                    PCpriors2[(i - 1) * 2 + 2] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                    PCpriors<-c(PCpriors,PCpriors2)
                  }
                  if(any(recStructures==0)){
                    PCpriors<-''
                    for(i in 1:nvarlevels){
                      if(recStructures[i]==0){
                        PCpriorstemp <- rep(NA, 6)
                        PCpriorstemp[1]<-paste0('for(j in 1:',nrow(eval(parse(text = paste0('rastpervar',i)))),'){')
                        PCpriorstemp[2]<-paste0('for(k in samprast',i,'[(j-1)*2+1]:samprast',i,'[j*2]){')
                        PCpriorstemp[3]<-paste0('rk',i,'[k]~dnorm(0,taurk',i,')')
                        PCpriorstemp[4]<-'}#k'
                        PCpriorstemp[5]<-'}#j'
                        PCpriorstemp[6]<-paste0('taurk',i,'<-pow(sdrk',i,',-2)')
                        PCpriorstemp[7]<-paste0('sdrk',i,'~dunif(0,10)')
                      }
                      if(recStructures[i]==1){
                        PCpriorstemp <- rep(NA, 7)
                        PCpriorstemp[1] <- "for(j in 1:Sample){"
                        PCpriorstemp[2] <- paste("for(l in 1:nsup", which(spatvarlevels == rastpervar$cov.level[i]), "[", coef.indices$index[coef.indices$level ==
                                                                                                                                               latentVarianceLevel], "]){", sep = "", collapse = "")
                        PCpriorstemp[3] <- paste("xc", which(spatvarlevels == rastpervar$cov.level[i]), "[j,l]~dnorm(0,tauxc", which(spatvarlevels ==
                                                                                                                                       rastpervar$cov.level[i]), ")", sep = "", collapse = "")
                        PCpriorstemp[4] <- "}#l"
                        PCpriorstemp[5]<- "}#j"
                        PCpriorstemp[6] <- paste("tauxc", which(spatvarlevels == rastpervar$cov.level[i]), "<-pow(sdxc", which(spatvarlevels ==
                                                                                                                                 rastpervar$cov.level[i]), ",-2)", sep = "", collapse = "")
                        PCpriorstemp[7] <- paste("sdxc", which(spatvarlevels == rastpervar$cov.level[i]), "~dunif(0,10)", sep = "", collapse = "")
                        PCpriorstemp[7]<-"}#j"
                      }
                      PCpriors<-c(PCpriors,PCpriorstemp)
                    }
                    PCpriors<-PCpriors[-1]
                  }
                }


                # error variance prior
                if (errorVarianceLevel == "sample" & multiSampsPerSubj == TRUE) {
                  errlines <- c("for(i in 1:Subj){", "for(j in (Sample[i]+1):Sample[i+1]){", "sigma[j]~dunif(0,10)", "}#j", "}#i")
                }
                if (errorVarianceLevel == "sample" & multiSampsPerSubj == FALSE) {
                  errlines <- c("for(j in 1:Sample){", "sigma[j]~dunif(0,10)", "}#j")
                }
                if (errorVarianceLevel == "subject" & multiSampsPerSubj == TRUE) {
                  errlines <- c("for(i in 1:Subj){", "sigma[i]~dunif(0,10)", "}#i")
                }
                if (errorVarianceLevel == "overall") {
                  errlines <- c("sigma~dunif(0,10)")
                }
                if (errorVarianceLevel == "subject" & multiSampsPerSubj == FALSE)
                  stop("If multiSampsPerSubj==FALSE then you cannot estimate the error variance on the subject level.
             If multiSampsPerSubj==FALSE is correct, then change the error variance to sample or overall.")

                # intercepts
                if (multiSampsPerSubj == TRUE) {
                  sampint <- c("for(i in 1:Subj){", "for(j in (Sample[i]+1):Sample[i+1]){", "v[j]~dnorm(0,tauv)", "}#j", "}#i")
                  subjint <- c("for(i in 1:Subj){", "bi1[i]~dnorm(0,taubi1)", "bi2[i]~dnorm(meanbi2[i],precbi2)", "meanbi2[i]<-(sdbi2/sdbi1)*rho*bi1[i]",
                    "}#i", "precbi2<-taubi2/(1-rho*rho)", "taubi1<-pow(sdbi1,-2)", "sdbi1~dunif(0,10)", "taubi2<-pow(sdbi2,-2)", "sdbi2~dunif(0,10)",
                    "rho~dunif(-1,1)", "tauv<-pow(sdv,-2)", "sdv~dunif(0,10)")
                }

                if (multiSampsPerSubj == FALSE) {
                  sampint <- c("for(j in 1:Sample){", "v[j]~dnorm(0,tauv)", "}", "tauv<-pow(sdv,-2)", "sdv~dunif(0,10)")
                }


                # model coefficients
                if (coefPrior == "sdunif") {
                  modcoef <- list()
                  for (i in 1:nvars) {
                    if (covariate.typesb[i] %in% c("binary", "continuous")) {
                      modcoef[[i]] <- rep(NA, 3)
                      modcoef[[i]][1] <- paste("beta", i, "~dnorm(0,taub", i, ")", sep = "", collapse = "")
                      modcoef[[i]][2] <- paste("taub", i, "<-pow(sdb", i, ",-2)", sep = "", collapse = "")
                      modcoef[[i]][3] <- paste("sdb", i, "~dunif(0,10)", sep = "", collapse = "")

                    }
                    if (covariate.typesb[i] == "categorical") {
                      modcoef[[i]] <- rep(NA, 3 * ncol(covsb[[i]]$mapping))
                      for (j in 1:ncol(covsb[[i]]$mapping)) {
                        modcoef[[i]][(j - 1) * 3 + 1] <- paste("beta", i, "_", j, "~dnorm(0,taub", i, "_", j, ")", sep = "", collapse = "")
                        modcoef[[i]][(j - 1) * 3 + 2] <- paste("taub", i, "_", j, "<-pow(sdb", i, "_", j, ",-2)", sep = "", collapse = "")
                        modcoef[[i]][(j - 1) * 3 + 3] <- paste("sdb", i, "_", j, "~dunif(0,10)", sep = "", collapse = "")
                      }

                      # unused code for same variance for each beta modcoef[[i]]<-rep(NA,(ncol(covsb[[i]]$mapping)+2)) for(j in 1:ncol(covsb[[i]]$mapping)){
                      # modcoef[[i]][j]<-paste('beta',i,'_',j,'~dnorm(0,taub',i,')', sep='', collapse='') }
                      # modcoef[[i]][ncol(covsb[[i]]$mapping)+1]<-paste('taub',i,'<-pow(sdb',i,',-2)', sep='', collapse = '')
                      # modcoef[[i]][ncol(covsb[[i]]$mapping)+2]<-paste('sdb',i,'~dunif(0,10)', sep='', collapse = '')
                    }
                  }
                  modcoef <- c("beta0~dnorm(0,taub0)", "taub0<-pow(sdb0,-2)", "sdb0~dunif(0,10)", paste0(unlist(modcoef)))

                  modcoefa <- list()
                  for (i in 1:nvarsa) {
                    if (covariate.typesa[i] %in% c("binary", "continuous")) {
                      modcoefa[[i]] <- rep(NA, 3)
                      modcoefa[[i]][1] <- paste("alpha", i, "~dnorm(0,taua", i, ")", sep = "", collapse = "")
                      modcoefa[[i]][2] <- paste("taua", i, "<-pow(sda", i, ",-2)", sep = "", collapse = "")
                      modcoefa[[i]][3] <- paste("sda", i, "~dunif(0,10)", sep = "", collapse = "")

                    }
                    if (covariate.typesa[i] == "categorical") {
                      modcoefa[[i]] <- rep(NA, 3 * ncol(covsa[[i]]$mapping))
                      for (j in 1:ncol(covsa[[i]]$mapping)) {
                        modcoefa[[i]][(j - 1) * 3 + 1] <- paste("alpha", i, "_", j, "~dnorm(0,taua", i, "_", j, ")", sep = "", collapse = "")
                        modcoefa[[i]][(j - 1) * 3 + 2] <- paste("taua", i, "_", j, "<-pow(sda", i, "_", j, ",-2)", sep = "", collapse = "")
                        modcoefa[[i]][(j - 1) * 3 + 3] <- paste("sda", i, "_", j, "~dunif(0,10)", sep = "", collapse = "")
                      }

                      # unused code for same variance for each beta modcoefa[[i]]<-rep(NA,(ncol(covsb[[i]]$mapping)+2)) for(j in 1:ncol(covsb[[i]]$mapping)){
                      # modcoefa[[i]][j]<-paste('alpha',i,'_',j,'~dnorm(0,taua',i,')', sep='', collapse='') }
                      # modcoefa[[i]][ncol(covsb[[i]]$mapping)+1]<-paste('taua',i,'<-pow(sda',i,',-2)', sep='', collapse = '')
                      # modcoefa[[i]][ncol(covsb[[i]]$mapping)+2]<-paste('sda',i,'~dunif(0,10)', sep='', collapse = '')
                    }
                  }
                  modcoefa <- c("alpha0~dnorm(0,taua0)", "taua0<-pow(sdb0,-2)", "sda0~dunif(0,10)", paste0(unlist(modcoefa)))

                  modcoef <- c(modcoefa, modcoef)
                }

                if (coefPrior == "dnorm") {
                  modcoef <- rep(NA, length(betanames))
                  for (i in 1:length(betanames)) {
                    modcoef[i] <- paste(betanames[i], "~dnorm(0,0.00001)", sep = "", collapse = "")
                  }
                  modcoef <- c("beta0~dnorm(0,0.00001)", modcoef)

                  modcoefa <- rep(NA, length(alphanames))
                  for (i in 1:length(alphanames)) {
                    modcoefa[i] <- paste(alphanames[i], "~dnorm(0,0.00001)", sep = "", collapse = "")
                  }
                  modcoefa <- c("alpha0~dnorm(0,0.00001)", modcoefa)

                  modcoef <- c(modcoefa, modcoef)
                }

                if (coefPrior == "Cauchy") {
                  modcoef <- rep(NA, length(betanames))
                  for (i in 1:length(betanames)) {
                    modcoef[i] <- paste(betanames[i], "~dt(0,0.16,1)", sep = "", collapse = "")
                  }
                  modcoef <- c("beta0~dt(0,0.01,1)", modcoef)

                  modcoefa <- rep(NA, length(alphanames))
                  for (i in 1:length(alphanames)) {
                    modcoefa[i] <- paste(alphanames[i], "~dt(0,0.16,1)", sep = "", collapse = "")
                  }
                  modcoefa <- c("alpha0~dt(0,0.01,1)", modcoefa)

                  modcoef <- c(modcoefa, modcoef)
                }

                # additional nodes to compare categorical coefficients
                if ("categorical" %in% covariate.typesb) {
                  addbetanodes <- list()
                  nodetemp <- list()
                  for (i in 1:nvars) {
                    if (covariate.typesb[i] == "categorical") {
                      # to provide names, need to get mapping matrix, column of betas and column of corresponding name
                      betanames2 <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      covarnames2 <- paste0(covariatesb[i], "_", colnames(covsb[[i]]$mapping))
                      namemap <- data.frame(betanames2, covarnames2)

                      nodetemp[[i]] <- paste0("beta", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      nodemat <- expand.grid(paste(nodetemp[[i]]), paste(nodetemp[[i]]))
                      nodemat <- nodemat[nodemat[, 1] != nodemat[, 2], ]

                      nodemat[, 3] <- 0
                      for (j in 1:nrow(nodemat)) {
                        if (as.numeric(strsplit(as.character(nodemat[j, 1]), "_")[[1]][2]) > as.numeric(strsplit(as.character(nodemat[j, 2]), "_")[[1]][2])) {
                          nodemat[j, 3] <- 1
                        }
                      }
                      nodemat <- nodemat[nodemat[, 3] == 1, ]
                      nodemat <- nodemat[order(nodemat[, 1]), ]
                      nodemat <- nodemat[, c(1, 2)]

                      nodemat <- as.data.frame(nodemat)
                      colnames(nodemat) <- c("beta1", "beta2")
                      nodemat$cov1 <- NA
                      nodemat$cov2 <- NA
                      for (j in 1:nrow(nodemat)) {
                        nodemat$cov1[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta1[j]])
                        nodemat$cov2[j] <- as.character(namemap$covarnames2[namemap$betanames2 == nodemat$beta2[j]])
                      }

                      # need to subtract columns
                      addbetanodes[[i]] <- rep(NA, nrow(nodemat))
                      for (j in 1:nrow(nodemat)) {
                        addbetanodes[[i]][j] <- paste0(nodemat$cov1[j], "_minus_", nodemat$cov2[j], "<-", nodemat$beta1[j], "-", nodemat$beta2[j])
                      }
                    }
                  }
                  addbetanodes <- paste(unlist(addbetanodes))


                  addalphanodes <- list()
                  nodetemp <- list()
                  for (i in 1:nvars) {
                    if (covariate.typesa[i] == "categorical") {
                      # to provide names, need to get mapping matrix, column of alphas and column of corresponding name
                      alphanames2 <- paste0("alpha", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      covarnames2 <- paste0(covariatesa[i], "_", colnames(covsb[[i]]$mapping))
                      namemap <- data.frame(alphanames2, covarnames2)

                      nodetemp[[i]] <- paste0("alpha", i, "_", seq(1:ncol(covsb[[i]]$mapping)))
                      nodemat <- expand.grid(paste(nodetemp[[i]]), paste(nodetemp[[i]]))
                      nodemat <- nodemat[nodemat[, 1] != nodemat[, 2], ]

                      nodemat[, 3] <- 0
                      for (j in 1:nrow(nodemat)) {
                        if (as.numeric(strsplit(as.character(nodemat[j, 1]), "_")[[1]][2]) > as.numeric(strsplit(as.character(nodemat[j, 2]), "_")[[1]][2])) {
                          nodemat[j, 3] <- 1
                        }
                      }
                      nodemat <- nodemat[nodemat[, 3] == 1, ]
                      nodemat <- nodemat[order(nodemat[, 1]), ]
                      nodemat <- nodemat[, c(1, 2)]

                      nodemat <- as.data.frame(nodemat)
                      colnames(nodemat) <- c("alpha1", "alpha2")
                      nodemat$cov1 <- NA
                      nodemat$cov2 <- NA
                      for (j in 1:nrow(nodemat)) {
                        nodemat$cov1[j] <- as.character(namemap$covarnames2[namemap$alphanames2 == nodemat$alpha1[j]])
                        nodemat$cov2[j] <- as.character(namemap$covarnames2[namemap$alphanames2 == nodemat$alpha2[j]])
                      }

                      # need to subtract columns
                      addalphanodes[[i]] <- rep(NA, nrow(nodemat))
                      for (j in 1:nrow(nodemat)) {
                        addalphanodes[[i]][j] <- paste0(nodemat$cov1[j], "_minus_", nodemat$cov2[j], "<-", nodemat$alpha1[j], "-", nodemat$alpha2[j])
                      }
                    }
                  }
                  addalphanodes <- paste(unlist(addalphanodes))
                }


                if (multiSampsPerSubj == TRUE) {
                  if (any(covariate.typesb == "categorical") == TRUE | any(covariate.typesa == "categorical") == TRUE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, sampint, subjint, modcoef, addbetanodes,
                      "})"))))
                  }
                  if (any(covariate.typesb == "categorical") == FALSE & any(covariate.typesa == "categorical") == FALSE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, sampint, subjint, modcoef, "})"))))
                  }
                }
                if (multiSampsPerSubj == FALSE) {
                  if (any(covariate.typesb == "categorical") == TRUE | any(covariate.typesa == "categorical") == TRUE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, sampint, modcoef, addbetanodes, "})"))))
                  }
                  if (any(covariate.typesb == "categorical") == FALSE & any(covariate.typesa == "categorical") == FALSE) {
                    assign("mymod", eval(parse(text = c("nimbleCode({", modbod, modsampPC, PCpriors, errlines, sampint, modcoef, "})"))))
                  }
                }
            }  # typeOfZero=='true'
    }


    outlist <- list(model = mymod,
                    covariates = covariatesb, covariateTypes=covariate.typesb, covariateLevels=covariate.levelsb)

    if(typeOfZero=='true'){
      outlist[['covariatesForBinary']]<-covariatesa
      outlist[['covariatesBinaryTypes']]<-covariate.typesa
      outlist[['covariatesBinaryLevels']]<-covariate.levelsa
    }
    outlist[['coefPrior']]<-coefPrior
    outlist[['multiSampsPerSubj']]<-multiSampsPerSubj
    outlist[['errorVarianceLevel']]<-errorVarianceLevel
    outlist[['latentVarianceLevel']]<-latentVarianceLevel
    outlist[['typeOfZero']]<-typeOfZero

    class(outlist) <- "PCModelList"
    return(outlist)

}







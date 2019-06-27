# imagingPC: Analysis of imaging mass spectrometry (IMS) data using a process convolution (PC) approach

In this R package, we employ a discrete PC approach in which a zero-centered latent process is convolved with a smoothing kernel function to account for spatial information.  For computational efficiency with large datasets, this package implements a semivariogram-based approach to estimate and fix the smoothing kernel function.  In the context of a Bayesian mixed model framework, this package writes and fits models to estimate the latent process at a limited set of locations called cupport istes while also incorporating covariates of interest.  Furthermore, this package incorporates the PC approach into left-censored and marginalized two-part models to account for different zero-generating processes.  This package was designed for IMS data, but the methods can be extended to imaging data collected over a regular grid.

The imagingPC package requires the R packages geoR, plyr, coda, ggplot2, grid, gridExtra, cowplot, and nimble.  The last package, nimble, also requires the installation of Rtools to write C++ code.

The research for which this R package was developed will be published in multiple papers in peer-reveiwed journals.  As soon as those papers are published, I will reference those publications to provide insight into and justification for the methods we use.

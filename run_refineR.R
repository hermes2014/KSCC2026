# Install refineR package from CRAN, if it is not installed.
pkg <- "refineR"
if (!requireNamespace(pkg, quietly = TRUE)) {
  install.packages(pkg)
}


# Load refineR package
library(refineR) 


# To open the help page and documentation with all important 
# information about the functions and their arguments 
?refineR


# Prefiltering, cleaning, and possible partitioning of the extracted 
# data from the laboratory information system is not included here and
# has to be done in advance. 
# To showcase the refineR application, an exemplary dataset is used 
# containing only a vector of (prefiltered) numeric test results. 
# Replace this by importing or loading your own dataset. 
data(testcase1)
data <- testcase1


# Run refineR estimation 
fit <- findRI(Data=data)


# Run refineR estimation with bootstrapping, by setting the argument 
# “NBootstrap” to a value > 0. In our experience, using a minimum 
# number of 200 bootstrap iterations leads to reasonable results.
# Please note that the computation time is increased considerably when # bootstrapping is activated. 
fit.bs <- findRI(Data=data, NBootstrap=200)


# Print summary of estimated model
print(fit)
 


# Print summary of estimated model with bootstrapping
print(fit.bs)

 

# Compute reference interval estimates for specified percentiles 
# (RIperc, default: c(0.025,0.975)) with confidence intervals for 
# specified confidence region (CIprop, default: 0.95) and specify 
# if you want to use the full data estimate as point estimate or the 
# median or mean from the bootstrap samples by setting pointEst to 
# "fullDataEst" (Default), "medianBS", or "meanBS". 
# We recommend to use the median of all bootstrap samples (medianBS) 
# here. 
getRI(fit.bs, RIperc=c(0.025, 0.975), CIprop=0.95, pointEst="medianBS")  

 
# Default plot function 
x11()
plot(fit)
 
# Optional: Plot model in Box-Cox transformed domain. 
# This step should only be carried out to assess the goodness of fit. # Please note that confidence intervals cannot be shown in this domain # and the reference interval is shown on the transformed scale, which # does not correspond to the original concentration scale.

# plot(fit, Scale = ”transformed”)

# Additional parameters that could be useful for the plot function
# (in addition to RIperc, CIprop and pointEst): 
# xlim, xlab, title,(ylab, ylim) 
# (xlim and xlab are just exemplary values and should be adjusted for 
# your data)

plot(fit.bs, RIperc=c(0.025, 0.975), CIprop=0.95, pointEst="medianBS",    
     xlim=c(0,40), xlab="Concentration [xx/yy]", title="Testset")
 
# Check visual representation of the result. If distribution appears
# to be skewed and shifted away from zero ( e.g. like the AST case),  # perform the estimation using the alternative model ‘modBoxCox’.
# This way, the two-parameter Box-Cox transformation will be utilized # and an optimal shift will be identified (the option ‘modBoxCoxFast’ # is a faster approximation of the shift with less accurate results)
fit.mbc <- findRI(Data=data, model="modBoxCox")

# Print summary of estimated model
print(fit.mbc)

 


# Default plot function 
plot(fit.mbc)

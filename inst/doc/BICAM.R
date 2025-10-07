## ----initial, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----packageload, message=FALSE, warning=FALSE--------------------------------
library(BTIME)


## ----datasetup----------------------------------------------------------------

# dfmif <- read.csv('~/mIF_Test_Data.csv', check.names = FALSE) # Load in data
# dat <- dfmif |> dplyr::select(suid,total,stage,M1,M2,M3) # Construct data to match setup of image below


## ----modelrun, echo=TRUE, results='hide'--------------------------------------

##################### Unstructured Model #####################
# posterior_unstr <- BICAM(dat = dat,
#                          M = 3,
#                          adapt = 2000,
#                          burn = 100,
#                          it = 100)

##################### Exponential Decay Model #####################

# WILL GET ERROR IF 'dis' MATRIX IS NOT SUPPLIED. 'dis' MUST ALSO BE INVERTIBLE.

# dis <- matrix(c(0,1,1,
#                  1,0,2,
#                  1,2,0),nrow = 3)

# posterior_expdecay <- BICAM(dat = dat,
#                          M = 3,
#                          adapt = 2000,
#                          burn = 100,
#                          it = 100,
#                          model = "ExpDecay",
#                          dis = dis)

##################### Tree Covariance Structure Model #####################

# WILL GET ERROR IF 'tree' MATRIX IS NOT SUPPLIED. 'tree' MUST ALSO BE INVERTIBLE.

# tree <- matrix(c(1,1,1,
#                 1,2,1,
#                 1,1,2),nrow = 3)

# posterior_tree <- BICAM(dat = dat,
#                          M = 3,
#                          adapt = 2000,
#                          burn = 100,
#                          it = 100,
#                          model = "Tree",
#                          tree = tree)

########### Scaled Tree

# posterior_treescaled <- BICAM(dat = dat,
#                          M = 3,
#                          adapt = 2000,
#                          burn = 100,
#                          it = 100,
#                          model = "TreeScaled",
#                          tree = tree)


##################### Multi-Level Tree Covariance Structure Model #####################

# WILL GET ERROR IF 'treelevels' MATRIX IS NOT SUPPLIED.

# treelevels <- list(matrix(c(1,1,1,
#                             1,1,1,
#                             1,1,1),nrow = 3),
#                    matrix(c(0,0,0,
#                             0,1,0,
#                             0,0,1),nrow = 3))

# posterior_tree <- BICAM(dat = dat,
#                          M = 3,
#                          adapt = 2000,
#                          burn = 100,
#                          it = 100,
#                          model = "TreeLevels",
#                          treelevels = treelevels)


## ----post1, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[1]]


## ----post2, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[2]]


## ----post3, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[3]]


## ----post4, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[4]]


## ----post5, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[5]]


## ----post6, echo=TRUE, results='hide'-----------------------------------------
# posterior_unstr[[6]]



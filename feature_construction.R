library(R.matlab)
library(pracma)
library(neuralnet)
library(emdist)

# craft the features that we use for classification

pindx <- NULL
for(i in c(2,3,5,7)) {  # get indices of the primes (as reference indices)
  pindx <- c(pindx, which(mnist$test$labels == i)[1])  # get the first index of a digit that is the prime (= reference image)
}

features <- matrix(data = rep(0, 10*4), nrow = 10)

# first feature: is the number strictly positive (indicator 0/1)
#features[,1] <- c(0,rep(1,9))  # only 0 is not positive; all other numbers are

## only for debugging purposes, we may put the number directly into the first feature, to later verify
## (when looking at point_pairs) that the primes and composites are put close to each other accordingly 
#features[,1] <- 0:9 
## -----

# second feature: is the number an integer
features[,2] <- rep(1,10) # all numbers are integers (but practically, a dot with commas could be recognized, making this indicator a "0")
# third feature: is it a prime?
features[1+c(2,3,5,7),3] <- 1  # only these numbers are primes; the "+1" is only to get the right index (starting from 1, rather than 0)
# last feature: earth mover distance (2D) between a digit and the set of primes
dindx <- NULL
for(d in 0:9) {  # iterate over all digits
  dindx <- which(mnist$test$labels == d)[2]  # retrieve another instance of the digit (including the primes themselves)
  # compute the minimum earth-mover distance to the primes
  digit <- matrix(mnist$test$images[dindx,], nrow=28)
  f <- Inf
  for(j in pindx) {   # iterate over the reference prime digit pictures 
    p <- matrix(mnist$test$images[j,], nrow=28)   # retrieve the prime reference digit picture
    f <- min(f, emd2d(digit,p, max.iter = 1000)) # get the earth mover distance between the digit picture and the reference prime
  }
  features[d + 1, 4] <- f
}
rownames(features) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
colnames(features) <- c("positive", "integer", "nontrivial_divisor", "similarity")

## determine the "small" and "large" distances
## NOTE: the variable "features" here is the same as "y" in the Octave code; since parts of the Octave implementation
## were taken over to here, we changed only the syntax, but not the variable naming (to remain closer to the Octave part)
y <- features
m <- 10
## BEGIN of what has been copied from Octave ################################################
## determine the closest and farest pair of points
d_min <- Inf
d_max <- -Inf
for(i in 1:m) {
  for(j in (i+1):m) {
    if (j > m) break;
    d = norm(y[i,] - y[j,], type = "2")
    d_min = min(d_min, d)
    d_max = max(d_max, d)
  }
}

## set values for points to be "close" or "far" based on their current
## geometric scattering. The choices here are arbitrary and have no mutual relation other than being "small" and "large"
smallDistance = 0.005 * d_min;  ## half of 1% of the minimum distance
largeDistance = 200 * d_max;  ## large multiple of the maximum distance
## END of what has been copied from Octave ################################################

h <- nchoosek(m,2)

## we use inverse distance weighting as the interpolation method instead
## since the usual packages all work only on 2D grids, we use a direct implementation here
idw_helper <- function(x, points, values) {
  N <- nrow(points)
  
  # look up if the point is one of the known points
  i <- which(apply(points, 1, function(p) return(all(abs(p - x) < 1e-7))))  # avoid numerical issues by == comparison
  
  if (length(i) == 0) {
    # point not found; now we spatially interpolate
    rho <- 5    # arbitrary choice (larger than the dimension)
    A <- 0
    B <- 0
    for(i in 1:N) {
      p <- points[i,]
      w <- 1 / norm(x - p, type = "2")^rho
      A <- A + w * values[i]
      B <- B + w
    }
    return(A / B)
  }
  else {
    return(values[i])  # return the point as found
  }
}

source("dbscan_customDist.R")

point_pairs <- matrix(rep(0,h*8), nrow = h, ncol=8)
colnames(point_pairs) <- c("f1", "f2", "f3", "f4", "g1", "g2", "g3", "g4") #, "dist")

####################################################
# desired classification, case #1: recognize primes
class_label <- features[,3] # let the classification be exactly whether the digit is a prime
class_label_case1 <- class_label   # save it for writing it to the Octave file

source("compute_pairwise_distances.R")

# here is our crafted distance, based on interpolation (as an alternative "baseline") to the eps-semimetric
# OBSERVE that the distances as desired are put in *exactly here*
crafted_distance <- function(x,y) { return(idw_helper(c(x,y), point_pairs, distances)) }

dbscan_result_unmodified <- dbscan_customDist(y, 2, mean(range(dist(y))))  # run with Euclidean distances and mean average distance between features
# use-case 1 (not intended, but demonstrated): accuracy can improve
dbscan_result_case1 <- dbscan_customDist(y, 2, 0.5*(smallDistance + largeDistance), crafted_distance)
# result: primes are exactly recognized

####################################################
# desired classification, case #2: let us include 1 as a prime number
#                     0 1 2 3 4 5 6 7 8 9
class_label <- c(0,1,1,1,0,1,0,1,0,0)
class_label_case2 <- class_label   # save it for writing it to the Octave file
source("compute_pairwise_distances.R")
crafted_distance <- function(x,y) { return(idw_helper(c(x,y), point_pairs, distances)) }
dbscan_result_case2 <- dbscan_customDist(y, 2, 0.5*(smallDistance + largeDistance), crafted_distance)

# here is our crafted distance, based on interpolation (as an alternative "baseline") to the eps-semimetric
# OBSERVE that the distances as desired are put in *exactly here*

####################################################
# desired classification, case #3: let us declare exactly all odd numbers as prime (which is incorrect for 1, 2 and 9)
#                     0 1 2 3 4 5 6 7 8 9
class_label <- c(0,1,0,1,0,1,0,1,0,1)
class_label_case3 <- class_label   # save it for writing it to the Octave file
source("compute_pairwise_distances.R")
crafted_distance <- function(x,y) { return(idw_helper(c(x,y), point_pairs, distances)) }
dbscan_result_case3 <- dbscan_customDist(y, 2, 0.5*(smallDistance + largeDistance), crafted_distance)


#fire test-cases against crafted_distance to make it recognizable as *not* being a metric!
test_metric <- function(nTests) {
  message("testing positive definiteness")
  for(i in 1:nTests) {
    x <- 5*rnorm(n = 4)
    y <- 5*rnorm(n = 4)
    d <- crafted_distance(x,y)  # should be > 0
    if (d < 0) {
      message(" positive definiteness failed on", x , ", ", y)
    }
    d <- crafted_distance(x,x)  # should give zero
    if (d > 0) {
      message(" positive definiteness failed on d(x,x)")
    }
  }
    
  message("testing symmetry")
  for(i in 1:nTests) {
    x <- 5*rnorm(n = 4)
    y <- 5*rnorm(n = 4)
    d1 <- crafted_distance(x,y)
    d2 <- crafted_distance(y,x)
    if (d1 != d2) {
      message(" symmetry failed on", x , ", ", y)
    }
  }
  
  message("testing the triangle inequality")
  for(i in 1:nTests) {
    x <- 5*rnorm(n = 4)
    y <- 5*rnorm(n = 4)
    z <- 5*rnorm(n = 4)
    dxy <- crafted_distance(x,y)
    dxz <- crafted_distance(x,z)
    dzy <- crafted_distance(z,y)
    if (dxy > dxz + dzy) {
      message(" triangle inequality failed on", x , ", ", y, ", ", z)
    }
  }
}

##### now, for k-means...
## since this was subject to strong roundoff errors anyway, we spared testing it


# now, to get this feature vector into MATLAB/Octave
writeMat("for_octave_script.mat", y = features, classlabels1 = class_label_case1, classlabels2 = class_label_case2, classlabels3 = class_label_case3)  # ...save it as a MAT file
# and in MATLAB/Octave, invoke: load("features.mat") to get the variable "y" accordingly, and the other class labels for all 3 cases
# ATTENTION: works well with manipulation_dbscan_demo3.m
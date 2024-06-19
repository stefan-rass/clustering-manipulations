library(dslabs)

mnist <- read_mnist()  # load from http source; may take *quite some time*
i <- 5
image(1:28, 1:28, matrix(mnist$test$images[i,], nrow=28)[ , 28:1], col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
mnist$test$labels[i]  # label of the picture

# convert image into a matrix:
matrix(mnist$test$images[i,], nrow=28)
# extract features

which(mnist$test$labels == 0)  # get all indices for which the digit (label) is "0"
which(mnist$test$labels %in% c(2,3,5,7))  # indices of digits that are primes {2,3,5,7}

for(d in 0:9) {
  i <- which(mnist$test$labels == d)[2]
  image(1:28, 1:28, matrix(mnist$test$images[i,], nrow=28)[ , 28:1], col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
  # display the pictures; we use the 2nd "instance" of each digit as reference in "feature_construction.r"
}


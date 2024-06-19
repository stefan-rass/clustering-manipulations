# helper code to compute pairs of vectors with desired distances based on (predefined) class labels
k <- 0
distances <- rep(0, h)  ## distances = zeros([h 1]);
## set the distance data vector to be exactly the pairwise distances reflecting the classification
for(i in 1:m) {
  for(j in (i+1):m) {
    if (j > m) break;
    k <- k + 1
    point_pairs[k,1:4] <- y[i,]
    point_pairs[k,5:8] <- y[j,]
    if (class_label[i] == class_label[j]) {
      distances[k] <- smallDistance
    } 
    else {
      distances[k] <- largeDistance
    }
  }
}

# variation of DBSCAN with a custom metric "custom_metric"
# added comments to relate this code to the Pseudocode as given on https://en.wikipedia.org/wiki/DBSCAN

dbscan_customDist <- function(X,minpts,EPS, custom_metric = function(x,y) { sqrt(sum((x-y)^2)) }) {
  C <- 0;
  nr <- dim(X)[1]
  assignments <- zeros(dim(X)[1],1);
  clustered <- zeros(dim(X)[1],1);
  for(i in 1:dim(X)[1]) {
    if(clustered[i] == 1) {
      next;
    }
    clustered[i] <- 1;    ## mark P as visited
    isneighbour <- NULL;
    neighbourcount <- 0;
    for(j in 1:dim(X)[1]) { ## N <- D.regionQuery(P, eps)  "isneighbor" <- N here
      #dist <- sqrt(sum((X[i,]-X[j,]).^2));
      dist <- custom_metric(X[i,],X[j,]);
      if(dist < EPS) {
        neighbourcount <- neighbourcount + 1;
        isneighbour <- c(isneighbour,j);
      }
    }
    if(neighbourcount<minpts) {  ## if sizeof(N) < MinPts
      next;  ## mark P as NOISE (not explicitly done)
    } 
    else {
      C <- C + 1;   ## C <- next cluster; now "expandCluster(P, N, C, eps, MinPts)" begins
      assignments[i] <- C;  ## add P to cluster C
      for(k in isneighbour) { ## for each point P' in N
        if (clustered[k]==0) { ## if P' is not visited
          clustered[k] <- 1  ## mark P' as visited
          b_isneighbour <- NULL
          b_neighbourcount <- 0
          for (j in 1:dim(X)[1]) { ## N' <- D.regionQuery(P', eps); "b_isneighbour" <- N' here
            #dist <- sqrt(sum((X[k,]-X[j,]).^2));
            dist <- custom_metric(X[k,],X[j,])
            if(dist<EPS) {
              b_neighbourcount <- b_neighbourcount + 1
              b_isneighbour <- c(b_isneighbour,j)
            }
          }
          if (b_neighbourcount >= minpts) { ## if sizeof(N') >= MinPts
            isneighbour <- c(isneighbour,b_isneighbour);  ## N <- N joined with N'
          }
        }
        ## the next line from the pseudocode is missing in this implementation
        ##if P' is not yet member of any cluster
        assignments[k] <- C;  ## add P' to cluster C
        ## unmark P' as NOISE if necessary => not necessary, since we never mark as noise
      }
    }
  }
  return(assignments)
}
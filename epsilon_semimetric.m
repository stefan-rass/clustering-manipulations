## Change (as of October 2nd, 2025): adapted the noise injection during the mapping of y's into the high-dimensional space, to make the method also work when there is a rank deficiency in the feature vector matrix
##       Status of the change: works as intended

## Parameters
## x, y: two points in the (low-dimensional) space
## s: scaling factor for the (deterministic) noise
## A: matrix to define the quadratic form (must be square)
function d = epsilon_semimetric(x, y, s, A)
  h = rows(A);  # just for efficiency, since we need this value many times hereafter
  ex = h - max(size(x));
  ## initialize the pseudorandom noise
  rand("seed", prod(x)*sum(y));  ## let the seed depend on *both* x and y
  #zi = [x, s * rand([1 ex])];  # this version fails if the feature vector matrix has a rank deficiency; the line below fixes this problem
  zi = [x, zeros([1 ex])] + s * rand([1 h]);
  rand("seed", prod(y)*sum(x));  ## let the seed depend on *both* x and y
  #zj = [y, s * rand([1 ex])];  # this version fails if the feature vector matrix has a rank deficiency; the line below fixes this problem
  zj = [y, zeros([1 ex])] + s * rand([1 h]);
  rand("seed", "reset");  ## avoid unwanted side-effects (for other routines)
  v = zi - zj;
  d = sqrt(v*A*v');
endfunction

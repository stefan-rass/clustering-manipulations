## Parameters
## x, y: two points in the (low-dimensional) space
## s: scaling factor for the (deterministic) noise
## A: matrix to define the quadratic form (must be square)
function d = epsilon_semimetric(x, y, s, A)
  ex = rows(A) - max(size(x));
  ## initialize the pseudorandom noise
  rand("seed", prod(x)*sum(y));  ## let the seed depend on *both* x and y
  zi = [x, s * rand([1 ex])];
  rand("seed", prod(y)*sum(x));  ## let the seed depend on *both* x and y
  zj = [y, s * rand([1 ex])];
  rand("seed", "reset");  ## avoid unwanted side-effects (for other routines)
  v = zi - zj;
  d = sqrt(v*A*v');
endfunction

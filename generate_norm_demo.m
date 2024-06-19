##Experiment: Choice of randomly scattered points, and independently random
##            distances assigned between them. Code to verify if the constructed
##            quadratic form delivers the desired distances
##Conducted on: August 24th, 2022
##Description: see the paper, version of October 6th, 2022
##Result: Successful

clear all ## avoid side effects
rand("seed", "reset");
randn("seed", "reset");
##rand("seed", 456789); ## make results reproducibe as they appear in the paper
##randn("seed", 456789); ## attention: the Gaussian random number generator needs to be seeded as well

m = 10 %50  ## number of points
ell = 4   ## number of dimensions
epsilon = 0.1
h = nchoosek(m,2) % number of difference vectors to consider
y = randn([m ell]);
s = epsilon / sqrt(h);

% pick random distances to craft our metric for;
% these need to be set accordingly to how far we wish to have the data points
% we leave it with random distances here, for the verification of the method
% to work with arbitrary choices (other experiments will make the choice of distance more systematically)
distances = 1 + 2*rand([h 1])

################################################
# here is the main construction of the metric

ex = nchoosek(m,2) - ell;   ## number of dimensions to pad with noise
y_prime = [y zeros([rows(y) ex])];  % embedding into the larger space
M = zeros([h,h]);
k = 0; ## just a counter to enumerate the pairs
for i = 1:length(y)
  for j = (i+1):length(y)
    k = k + 1;
    % pick a fresh point in the proximity of the given data point to compute the difference
    zi = [y(i,:), s * rand([1 ex])];
    assert(norm(y_prime(i,:) - zi, 2) <= epsilon)   % y_prime is just y with zero-padding for equal dimensions
    zj = [y(j,:), s * rand([1 ex])];
    assert(norm(y_prime(j,:) - zj, 2) <= epsilon)
    M(:,k) = (zi - zj)';
  endfor
endfor

## M should have all linearly independent rows; ...print just for verification
disp("rank of M is ")
rank(M)

% construct the quadratic form taking the desired values
idx = 2:h;
B = null(M(:,idx)')';
A_1 = B'*B;
lambda = M(:,1)' * A_1 * M(:,1);
A = A_1 / lambda * distances(1)^2;
disp([1 norm(A_1, 2) lambda])
for i = 1:(h-1)  % iterate and exclude another vector each time
  idx(i) = i;
  B = null(M(:,idx)')';
  A_i = B'*B;  % actually, it is A_(i+1), but this would not be a valid variable identifier :-)
  lambda = M(:,i+1)' * A_i * M(:,i+1);
  A_i = A_i / lambda * distances(i+1)^2;
  A = A + A_i;
endfor

crafted_norm = @(x) (sqrt(x'*A*x));

% scaled version to check if we cannot at least get a scalar multiple
lambda_max = max(eig(A));
A_scaled = A / lambda_max;
scaled_crafted_norm = @(x) (sqrt(x'*A_scaled*x));

# end of constructions of norms (metric is constructed in the other experimental scripts)
################################################

% test if the quadratic form-norm delivers the sought distances
disp(['point pair  ' 'desired  ' 'computed distance' ])
for i = 1:h
  disp([i distances(i) crafted_norm(M(:,i))])
endfor

% test if the quadratic form-norm delivers the sought distances
disp(['point pair  ' 'desired_scaled ' 'computed distance with scaling' ])
for i = 1:h
  disp([i distances(i)/sqrt(lambda_max) scaled_crafted_norm(M(:,i))])
endfor

% verify that the scaled norm is <= 2-norm (...trivial...)
for trials = 1:1000
  y = randn([h 1]);  % ATTENTION: need to draw points from larger-dimensional space!
  assert(scaled_crafted_norm(y) <= norm(y,2))
endfor
disp('finished (all checks ok)');

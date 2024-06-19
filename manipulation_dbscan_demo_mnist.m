##Experiment: Manipulation of DBSCAN with the MNIST Dataset prepared in the R scripts
##Conducted on: April 3rd, 2024
##Description: see the paper, version of October 6th, 2022
##  April 2nd, 2024: bugfix in the loops (line 56 and others)
##Result: Successful

clear all
load('for_octave_script.mat')   # load class labels and data from R scripts

rand("seed", "reset"); ## make sure we always have fresh random points and classes
randn("seed", "reset");
##rand("seed", 456789); ## make results reproducibe as they appear in the paper
##randn("seed", 456789); ## attention: the Gaussian random number generator needs to be seeded as well

m = 10
ell = 4
y = 10*randn([m ell]);  # OBSERVE that we may even *overwrite* the features, i.e., entirely ignore them, but the classification still works
number_of_classes = 2;

disp(y);

## assign random classifications
#class_label = randi(number_of_classes, m, 1);

# take the classes as desired in the experiments (from the paper)
# Experiment 2: classlabels1: only primes
# Experiment 3: classlabels2: primes incl. the number 1
# Experiment 4: classlabels3: all odd numbers
#class_label = classlabels1 + 1 # "+1" to avoid 0 as a class (would cause errors with indexing below)
#class_label = classlabels2 + 1 # "+1" to avoid 0 as a class (would cause errors with indexing below)
#class_label = classlabels3 + 1 # "+1" to avoid 0 as a class (would cause errors with indexing below)

## we let the point y(i) have class "class_label(i)" hereafter

## determine the closest and farest pair of points
d_min = Inf;
d_max = -Inf;
for i = 1:m
  for j = (i+1):m
    d = norm(y(i, :) - y(j, :), 2);
    d_min = min(d_min, d);
    d_max = max(d_max, d);
  endfor
endfor

## set values for points to be "close" or "far" based on their current
## geometric scattering. The choices here are arbitrary and have no mutual relation other than being "small" and "large"
smallDistance = 0.005 * d_min;  ## half of 1% of the minimum distance
largeDistance = 200 * d_max;  ## large multiple of the maximum distance

## assign the virtual neighbors in a distance that is "small"
epsilon =  smallDistance
ex = nchoosek(m,2) - ell;  ## dimensions to extend for the embedding in the larger space
y_prime = [y zeros([rows(y) ex])];  % embedding into the larger space

h = nchoosek(m,2) ## number of difference vectors to consider
M = zeros([h,h]);
distances = zeros([h 1]);
s = 0.9 * epsilon / sqrt(ex + ell);  ## scaling factor for the noise, so that the noise-points are less than epsilon far from the original points

################################################
# here is the main construction of the metric

k = 0;
for i = 1:size(y)(1)  % iterate over all rows of y
  for j = (i+1):size(y)(1)
    k = k + 1;
    % pick a fresh point in the proximity of the given data point to compute the difference
    rand("seed", prod(y(i,:))*sum(y(j,:)))  ## let the seed depend on *both* y(i) and y(j)
    zi = [y(i,:), s * rand([1 ex])];  ## now, the sequence is pseudorandom
    assert(norm(y_prime(i,:) - zi, 2) <= epsilon)   % y_prime is just y with zero-padding for equal dimensions

    rand("seed", prod(y(j,:))*sum(y(i,:)))  ## let the seed depend on *both* y(j) and y(i)
    zj = [y(j,:), s * rand([1 ex])];
    assert(norm(y_prime(j,:) - zj, 2) <= epsilon)
    M(:,k) = (zi - zj)';

    ## set the desired distance to be "small" if the two points belong to the
    ## same class, and "large" if they shall belong to different classes
    if (class_label(i) == class_label(j))
      distances(k) = smallDistance;
    else ## different classes, so the points should be separated widely
      distances(k) = largeDistance;
    endif
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

crafted_norm = @(y) (sqrt(y'*A*y));

% test if the quadratic form-norm delivers the sought distances
disp(['point pair  ' 'desired  ' 'computed distance' ])
for i = 1:h
  disp([i distances(i) crafted_norm(M(:,i))])
endfor

% scaled version to check if we cannot at least get a scalar multiple
lambda_max = max(eig(A));
A_scaled = A / lambda_max;
scaled_crafted_norm = @(y) (sqrt(y'*A_scaled*y));
% test if the quadratic form-norm delivers the sought distances
disp(['point pair  ' 'desired_scaled ' 'computed distance with scaling' ])
for i = 1:h
  disp([i distances(i)/sqrt(lambda_max) scaled_crafted_norm(M(:,i))])
endfor

crafted_metric = @(y,y_prime) (epsilon_semimetric(y, y_prime, s, A));

# end of constructions of metric and norms
################################################


intra_cluster_distances = zeros([1 number_of_classes]) - Inf;  ## for max-aggegation
inter_cluster_distances = zeros([number_of_classes number_of_classes]) + Inf;   ## for min-aggegation
k = 0;
for i = 1:(size(y)(1))
  for j = (i+1):(size(y)(1))
    k = k + 1;

    #d = epsilon_semimetric(y(i,:), y(j, :), s, A);   ## A ... unscaled matrix
    d = crafted_metric(y(i,:), y(j,:));

    ci = class_label(i);
    cj = class_label(j);
    if ci == cj
      intra_cluster_distances(ci) = max(intra_cluster_distances(ci), d);
     else
      icd = inter_cluster_distances(ci, cj);
      inter_cluster_distances(ci, cj) = min([icd d]);
      inter_cluster_distances(cj, ci) = min([icd d]);
    endif
  endfor
endfor
disp('-----------------------------');
disp('analysis of data with noise');
disp('unscaled distances inside a cluster');
disp(intra_cluster_distances);
disp('distances between clusters');
disp(inter_cluster_distances);

disp('middle distance')
disp(0.5*(smallDistance + largeDistance))

#[assignments,C] = dbscan(X,minpts,EPS,custom_metric)
[assignments,C] = dbscan_customDist(y, 2, 0.5*(smallDistance + largeDistance), crafted_metric);
disp([assignments, class_label]);
## for the verification of a bijection between the class labels assigned by DBSCAN and the ones we desired,
## we need the set of labels to be identical. DBSCAN assigns class labels counting from zero, while we assign from number 1
## -> we correct this if necessary by shifting all labels by +1
if (min(assignments) == 0)  ## just correct in case that DBSCAN has started with class label 0, whereas we started with class label 1
  assignments += 1;
endif

X = [assignments, class_label]; # should be "the same" ("isomorphic")
## we check if there is a bijection between "assignments" and "class_label"
Y = unique(X, "rows"); ## remove all duplicate rows
## check if we can permute the left and right columns (thinking of these as pre-images and images)
## to give the same set. If so, then obviously we have pi_1(y) == pi_2(f(y)), and f(y) = pi_2^(-1)(pi_1(y)) is a permutation and hence bijective
assert(!any(sort(Y(:,1)) != sort(Y(:,2))));

disp('cluster assignment as desired, up to isomorphism');

##Experiment: Manipulation of k-means
##Conducted on: September 29, 2022
##Description: see the paper, version of October 6th, 2022
##Result: Successful in many, but not all cases. Specifically,
##        the matrix M contains often NaN values. This effect is reproducible
##        by seeding the PRNG with the value 12345
##Conjectured explanation: numerical issues with the required linear independence
##Status: overall SUCCESSFUL
## - Manipulation works in all cases where there were no numerical errors
##   but does so only if we provide the cluster centers pre-computed as the
##   average of the desired cluster points
## - experiment FAILED, if we provide *other* cluster centers, e.g., any (arbitrary)
##   neighbor z_{i,j} of a point y_j. Explanation: k-Means in that case will frequently
##   update the cluster centers, thus "fiddling" with the design of distances that we want
## - Conjecture from above unverified (as of October 6th, 2022)

clear all
rand("seed", "reset"); ## make sure we always have fresh random points and classes
randn("seed", "reset");
##rand("seed", 456789); ## make results reproducibe as they appear in the paper
##randn("seed", 456789); ## attention: the Gaussian random number generator needs to be seeded as well
## note that seeding with the value "12345" will produce the numeric issues also reported in the paper

m = 10 %50
ell = 2 %4
y = 10*randn([m ell]);
number_of_classes = 3;

disp(y);

## assign random classifications
class_label = randi(number_of_classes, m, 1);
## we let the point y(i) have class C(i) hereafter

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

## in foresight of k-means, we determine a "virtual" cluster center to define the points as being "close" to the cluster-center
## while the cluster centers themselves should be far apart from one another
k_means_centers = zeros([number_of_classes ell]);
k_means_cluster_counts = zeros([number_of_classes 1]);
for i = 1:length(y)
  ci = class_label(i);
  k_means_centers(ci,:) += y(i,:);
  k_means_cluster_counts(ci)++;
endfor
for i = 1:number_of_classes
  k_means_centers(i,:) /= k_means_cluster_counts(i);
endfor
## now, we add the cluster centers to the data points, and set the distances accordingly
x_original = y;  ## make a "backup-copy" to restore the original data later on
y = [k_means_centers; y]; ## put the new points first, to test with the ell-means default of taking the first ell points as cluster centers (which they are, if we put it like this)
class_label = [(1:number_of_classes)'; class_label]; ## label the centers accordingly
m += number_of_classes; ## now we have a few more points to consider

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
for i = 1:length(y)
  for j = (i+1):length(y)
    k = k + 1;
    % pick a fresh point in the proximity of the given data point to compute the difference

    % pick the point in *dependence* of the respective neighbor that it refers to,
    % by using a bivariate quadratic multivariate form.

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

crafted_metric = @(y,y_prime) (epsilon_semimetric(y, y_prime, s, A_scaled));  ## scaled version

# end of constructions of metric and norms
################################################

intra_cluster_distances = zeros([1 number_of_classes]) - Inf;  ## for max-aggegation
inter_cluster_distances = zeros([number_of_classes number_of_classes]) + Inf;   ## for min-aggegation
k = 0;
for i = 1:length(y)
  for j = (i+1):length(y)
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

#[assignment,centers] = kmeans(X,ell, centers = 0, maxIter = 20000, custom_metric)
## pre-computed cluster centers; this will make k-means "not update/change" the cluster centers over iterations
[assignment,centers] = kmeans_customDist(y, number_of_classes, centers = k_means_centers, maxIter = 100, crafted_metric)
[assignment, class_label]
assert(!any(assignment != class_label));
disp('all class labels assigned as desired');


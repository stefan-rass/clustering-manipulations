# Scripts for the paper
## Metricizing the Euclidean Space towards Desired Distance Relations in Point Clouds
to appear in: IEEE Transactions on Information Forensics and Security

Recommended citation: Stefan Rass, Sandra König, Shahzad Ahmad, und Maksim Goman, "Metricizing the Euclidean Space Toward Desired Distance Relations in Point Clouds", IEEE Trans. Inf. Forensics Secur., Vol. 19, pp. 7304–7319, 2024, doi: 10.1109/TIFS.2024.3420246.

BibTeX: 

```
@article{rass_metricizing_2024,
	title = {Metricizing the {Euclidean} {Space} {Toward} {Desired} {Distance} {Relations} in {Point} {Clouds}},
	volume = {19},
	copyright = {https://creativecommons.org/licenses/by/4.0/legalcode},
	issn = {1556-6013, 1556-6021},
	url = {https://ieeexplore.ieee.org/document/10574843/},
	doi = {10.1109/TIFS.2024.3420246},
	urldate = {2024-08-05},
	journal = {IEEE Transactions on Information Forensics and Security},
	author = {Rass, Stefan and König, Sandra and Ahmad, Shahzad and Goman, Maksim},
	year = {2024},
	pages = {7304--7319},
}
```


---

* `generate_norm_demo.m`
  Verification of the construction of a norm to give desired distances 
* `manipulation_kmeans_demo.m`
  Experiment: Manipulation of the k-means algorithm
* `manipulation_dbscan_demo.m`
  Experiment: Manipulation of the DBSCAN algorithm
  Note that the version `manipulation_dbscan_demo_mnist.m` is different herein only in containing three different pre-desired classifications that we used to test the manipulation on the MNIST dataset (in Figure 4); details are explained below
* `epsilon_semimetric.m`
  Implementation of the epsilon semimetric, based on the quadratic form and noise magnitude for automatic embedding
* `dbscan_customDist.m`
  Adapted version of DBSCAN with custom distance function to be provided externally. 
* `kmeans_customDist.m`
  Adapted version of k-keans with custom distance function to be provided externally
* `kmeans_example_run.log` and `dbscan_example_run.log`
  example executions from running the code, showing all outputs
  
The Implementations of DBSCAN and k-Means have been adapted from the Implementation found here: https://github.com/devil1993/DBSCAN
The difference is only in 2 positions: (i) replacing the Euclidean distance by a function that can be provided externally, and (ii) adding a parameter to the procedure to externally supply the metric to be used.

The call-sequence on the Octave/Matlab prompt to run the examples is

`> manipulation_dbscan_demo`

`> manipulation_kmeans_demo`

The latter was observed to (frequently) fail due to numeric roundoff errors (which cause troubles in the internal eigenvalue computation routines). 

The code has been tested to run with Octave, version 9.1.0

## for the MNIST Experiments
These were more convenient to run in the R environment (www.r-project.org), which used the following files:
* `compute_pairwise_distances.R`:
  helper function to compute the actual distances between feature points
* `dbscan_customDist.R`
  an R implementation of the above dbscan_customDist.m
* `MNIST_image_extraction.R`
  a script to extract the features and do the "base-line" comparison with the inverse distance weighting (indicated in Experiment 5); source this file and look into the output variables "dbscan_result_case1", "dbscan_result_case2" and "dbscan_result_case3", corresponding to the results of Exp2, Exp3 and Exp4 in Figure 4 (and confirming them).
  On the R-prompt, you may invoke
  
`> test_metric(n)`

  to run, e.g., n random tests of the inverse distance weighting to be recognizable as *not* being a metric (this verification is analytically assured for the construction in the paper, except for the allowed additive error in the triangle inequality)
* `MNIST_image_extraction.R`
  script to load the images from the MNIST database
* `manipulation_dbscan_demo_mnist.m`
  This file has, in lines 29, 30 and 31, the desired classifications of the numbers 0, 1, ..., 9 encoded, to make DBSCAN produce exactly these desired outputs. Simply "include" the desired classification by uncommenting it (and commenting out the others) directly in the file. The actual values are produced with the supplementary R scripts and loaded from their output (the file "for_octave_script.mat", not included in the repo as being auxiliary intermediate only)
  
On the R prompt, one can source the following files in this order:

`> source("MNIST_image_extraction")`

`> source("MNIST_image_extraction")`

and, then finally, on the Octave prompt, invoke

`> manipulation_dbscan_demo_mnist`

to let the manipulation run with the desired classifications. We emphasize that the whole feature extraction in the R scripts is a typical and in fact crucial part of any AI data engineering, but the *actual feature values are entirely irrelevant* in this attack: Indeed, the Octave script actually overwrites the features with random values, and still accomplishes the desired classification. Only the "classlabels" variable, produced inside the R-scripts, is relevant for the Octave code.

## Remarks and Updates
Theorem 3 in the paper needs the unstated addtional hypothesis that the input points from $Y$ are all linearly independent. This is easy to accomplish a priori, by adding pseudo-features that are drawn from an absolutely continuous distribution (alluding to Lemma 1 whose proof argument also shows that the resulting vectors will be linearly independent). Alternatively, leaving the theorem in its hypothesis and claim unchanged, one can adapt the proof as follows: we create the noisy vectors as $z_{i,j}=y'_i+f(y_i,y_j)$, i.e., adding an $h$-dimensional noise vector to the canonic embedding $y'_i=(y_i,0)\in\mathbb{R}^h$ instead of just appending the noise to the original data. The proof, literally, remains intact upon changing the argument accordingly, and *likewise remains Theorem 3 correct as stated, under this change to its proof.*
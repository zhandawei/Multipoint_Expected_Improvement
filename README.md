## The Multipoint Expected Improvement (a.k.a *q*EI) Criterion 

* This is a MATLAB implementation of the multipoint expected improvement algorithm[1,2,3]. I am interesting in the multipoint expected improvement criterion, but I only find a implementation in R DiceOptim package[4] and could not find a MATLAB implementation. So, I made this implementation.

* I referred to the *qEI.R* function in the DiceOptim package[4] when doing the coding, which computes the analytical expression of the multipoint expected improvement criterion using the formula in [3].

* Specifically, I used the *qsimvnv* code of Alan Genz [5] for the numerical computation of multivariate normal distribution values, and the modified Cholesky algorithm[6] for the Cholesky decomposition in this code.

* Also, I have made some small changes to the *predictor* funtion of the DACE toolbox[7] to output the covaricance matrix of the inputs.

* I have verified this code against the Monte Carlo method.

* Feel free to contact me if you find any error in this code.



### reference
1. Schonlau, M. (1997). Computer experiments and global optimization, University of Waterloo.
2. Ginsbourger, D., R. Le Riche and L. Carraro (2010). Kriging Is Well-Suited to Parallelize Optimization. Computational Intelligence in Expensive Optimization Problems. Y. Tenne and C.-K. Goh, Springer Berlin Heidelberg. 2: 131-162.
3. Chevalier, C. and D. Ginsbourger (2013). Fast Computation of the Multi-Points Expected Improvement with Applications in Batch Selection. Learning and Intelligent Optimization. G. Nicosia and P. Pardalos, Springer Berlin Heidelberg. 7997: 59-69.
4. Roustant, O., D. Ginsbourger and Y. Deville (2012). "DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization." Journal of Statistical Software 51(1): 1-55. [Links](https://cran.r-project.org/web/packages/DiceOptim/index.html) 
5. Alan Genz. Numerical Computation of Multivariate Normal Probabilities, in J. of Computational and Graphical Stat., 1(1992), pp. 141-149, WSU Math, PO Box 643113, Pullman, WA 99164-3113. [Links](http://www.math.wsu.edu/faculty/genz/software/software.html)
6. S. H. Cheng and N. J. Higham. A modified Cholesky algorithm based on a symmetric indefinite factorization. SIAM J. Matrix Anal. Appl., 19(4):1097-1110, 1998. [Github links](https://github.com/higham/modified-cholesky)
7. S. N. Lophaven, H. B. Nielsen, and J. Sodergaard, DACE - A MATLAB Kriging Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical Modelling, Technical University of Denmark, 2002.

# UQ_benchmark_models
This is a collection of benchmark models that I have used for uncertainty 
quantification (UQ) purposes. They are parametric models, the parameters of 
which can be easily recasted as random variables (RVs). They have been 
predominantly used to test own UQ software, as well as to perform UQ studies 
and develop surrogate models with UQLab, SPINTERP, Chaospy, OpenTURNS, and the 
MATLAB-Sparse-Grids-Kit. UQ studies using the models provided in this repository 
are available in the following works:

@article{loukrezis2019assessing,
author  = {Dimitrios Loukrezis and Ulrich  Römer and Herbert  De Gersem},
title   = {Assessing the performance of Leja and Clenshaw-Curtis collocation for 
           computational electromagnetics with random input data},
journal = {International Journal for Uncertainty Quantification},
issn    = {2152-5080},
year    = {2019},
volume  = {9},
number  = {1},
pages   = {33--57}
}

@article{loukrezis2019approximation,
author = {{Loukrezis}, Dimitrios and {De Gersem}, Herbert},
title  = "{Approximation and Uncertainty Quantification of Stochastic Systems 
          with Arbitrary Input Distributions using Weighted Leja 
		  Interpolation}",
journal = {arXiv e-prints},
year    = "2019",
eid     = {arXiv:1904.07709},
}

@phdthesis{loukrezis2019adaptive,
author = {Dimitrios Loukrezis},
title  = {Adaptive Approximations for High-DImensional Uncertainty 
          Quantification in Stochastic Parametric Electromagnetic Field 
		  Simulations}, 
school = {TU Darmstadt},
year   = {2019}		  
}
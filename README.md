# multiBeyn
Contour integral solver for multiparameter eigenvalue problems

We implement a contour solver for two and three parameter eigenvalue problems. 
The code in this repository is based on the paper:

[1] Emil Graf and Alex Townsend, "A Contour Method for Multiparameter Eigenvalue Problems," 2026.

The folder biBeyn gives three versions of the algorithm, one for analytic problems, one for quadratic problems, and one for linear problems, 
along with four test files.
The folder triBeyn gives a single proof-of-concept implementation for linear three parameter eigenvalue problems, and a single small random test.

# Files needed

We have included the functions legpts and besselroots from chebfun, due to the University of Oxford and the chebfun developers, http://www.chebfun.org/.

The examples in test_SL and test_ARMA11_2 require files from MultiParEig, https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig, Bor Plestenjak, 2025. To run these examples, first download and install MultiParEig. Both ARMA examples are based on code from MultiParEig, B. Plestenjak and A. Muhic, University of Ljubljana, P. Holoborodko, Advanpix LLC.

# Tests

The test files reproduce the numerical experiments from [1]. Originally the experiments were performed on an AMD Ryzen 9 5950x using MATLAB r2023a. Due to variations in the order of operations, MATLAB's parfor is not strictly deterministic, so the experiments are not precisely reproducible, but are as close as possible.

With 12-16 workers, all experiments take at least several hours to run. The longest is test_ARMA11_2, which takes ~36 hours.

# Delay-Differential Equations

We sample a segment of the critical-delay curve for a DDE with two independent delays. This is an analytic, nonpolynomial two parameter eigenvalue problem. See section 7.1 of [1].

# Multiparameter Sturm-Liouville Problems

We find eigenvalues of a two parameter Sturm-Liouville problem given in section 7.2 of [1] for a range of large oscillation numbers. This is a linear two parameter eigenvalue problem.

# ARMA(1,1) Model

We find critical points for an ARMA(1,1) model, for a long sequence from Fatigue Crack Propagation Benchmark, GDR 3651 FATACRACK, https://doi.org/10.5281/zenodo.1478472,
and a short sequence from M.E.Hochstenbach, T.Kosir, B.Plestenjak, "On the Solution of Rectangular Multiparameter Eigenvalue Problems," arXiv 2212.01867. See section 7.3 of [1]. This is a quadratic two parameter eigenvalue problem.

#  triBeyn

We give a proof-of-concept test for triBeyn on a generic three parameter linear eigenvalue problem. See section 7.4 of [1].






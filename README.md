# RamsesPatchApproxMGSolver
Author: Hans A. Winther, University of Oxford (2015)

An implementation of the fast approximate method of [Winther &amp; Ferreira 2015](https://arxiv.org/abs/1403.6492) to perform N-body simulations of modified gravity models with screening on the chameleon type (works for all m(a), beta(a) models [Brax et al. 2012](https://arxiv.org/abs/1203.4812)). 

# Installation and use

- Place the patch in the patch/ folder of [RAMSES](https://bitbucket.org/rteyssie/ramses/overview) and use the Makefile included. Alternatively add DEFINES += -DMGSOLVER ; PATCH = ../patch/mgsolver and user\_defined\_functions.o to POISSONOBJS in your makefile.

- Specify the functions defining your model in user\_defined\_functions.f90 or compile with -DFOFR or -DSYMMETRON to use already implemented models.

- Two example models have been implemented: 1) [Hu-Sawicky f(R) gravity](https://arxiv.org/abs/0705.1158) 2) [Symmetron model](https://arxiv.org/abs/1001.4525).

- NB: This version of the patch has not been properly tested yet so use it carefully.

# Namelist parameters

- The parameters are specified in POISSON\_PARAMS in the namelist file

- mg\_solver\_active (logical; false) Use the modified version if true. Standard RAMSES otherwise.

- include\_screening (logical; false) If true use the screening method when solving. Otherwise we just solve the linear Klein-Gordon equation for the scalar-field.

- epsilon\_mg (real; 1d-5) Convergence criterion for the MG solver (1d-4 is used for the standard Poisson solver in RAMSES)

- n\_fofr, fofr0 (real; 1.0, 1d-5) Only if compiled with -DFOFR. [Hu-Sawicky f(R) gravity](https://arxiv.org/abs/0705.1158) parameters

- assb, lambda\_0, beta\_0 (real; 0.5, 1.0, 1.0) Only if compiled with -DSYMMETRON. [Symmetron](https://arxiv.org/abs/1001.4525) parameters

- myparam (real array; 0.0, 0.0, 0.0, 0.0, 0.0) User-defined parameters that can be used in user\_defined\_functions.f90

- eps\_to\_start\_solving (real; 0.001) The solver is only activated if the maximum |Geff(k,a)/G - 1| > eps\_to\_start\_solving

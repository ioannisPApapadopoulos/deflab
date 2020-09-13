# README #

### Synopsis ###

Deflation in MATLAB

This is a MATLAB library that implements the deflation method of Farrell, Birkisson and Funke (2015) doi:10.1137/140984798 

Deflation is used to discover multiple solutions of nonlinear problems. In most iterative methods, such as Newton's method, the user is required to pick an initial guess. In order to discover all the solutions of the nonlinear problem, the user often has to resort to randomly selecting different initial guesses in hopes of converging to different solutions.

With deflation the first solution is discovered as normal. However, once it has been discovered, the solution is then deflated away. The iterative algorithm can be run again - from the same initial guess - to discover subsequent solutions. Deflation will not allow the iterative algorithm to converge to an already discovered solution nor will it remove any solutions that have not yet been discovered. However, convergence is only guaranteed with a globally convergent method. 

This library also includes implementations of Benson and Munson's rsls and ssls active-set solvers (doi:10.1080/10556780500065382) and Hintermuller, Ito and Kunisch's primal-dual active set strategy (doi:10.1137/S1052623401383558)

### Dependencies ###

This code has been written in MATLAB.

### Code examples ###

The easiest way to learn how to use deflation is to examine the examples in examples/. Start with examples/parabola.

### Contributors ###

Ioannis P. A. Papadopoulos (ioannis.papadopoulos@maths.ox.ac.uk)

### License ###

GNU LGPL, version 3.

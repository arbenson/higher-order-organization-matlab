Higher-order organization of complex networks
--------

This is experimental Matlab code for the methods and some of the examples in

Higher-order Organization of Complex Networks.
Austin R. Benson, David F. Gleich, and Jure Leskovec.
Science, vol. 353, no. 6295, pp. 163-166, 2016.

If you use this code in a publication, please cite this paper.

The implementations used here do not use the most efficient triangle enumeration
algorithms.  Instead, they use the simplest matrix computations.  Thus, the code
here may not scale to super large networks.  The
[C++ code](http://snap.stanford.edu/higher-order/) uses all of the efficient
counting techniques.

See the [project page](http://snap.stanford.edu/higher-order/) for more
information about the methods.

Dependencies
--------
* [MatlabBGL](https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/) for graph processing in Matlab.
  This is package is required to run the examples, but is not necessary to use the rest of the code, such
  as the MotifAdjacency() function.

Examples
--------
* `celegans_example.m` reproduces some of the results in Figure 2 of the paper.
* `foodweb_example.m` reproduces some of the results in Section S7.1 of the paper.


[![Build Status](https://travis-ci.org/riklund/MScComputationSoftware.png)](https://travis-ci.org/riklund/MScComputationSoftware)


#MScComputationSoftware

A program used for computations for my Master's Thesis project. It basically solves the Schrödinger equation in a many-body particle basis expansion. This is done by constructing large matrices and finding their eigenvalues. 

![Resonance States Illustration](images/resonance.png?raw=true)
The figure shows an illustration of solutions to the Schrödinger equation, with a marked resonance state.


A general documentation of the purpose and functionality of this software, as well as the underlying mathematical principles, can be found in my thesis (http://publications.lib.chalmers.se/records/fulltext/197684/197684.pdf). It is recommended to read this thesis before looking at the code, in order to get a sense of what the code actually does.

Note also that for large matrices, it may be better to instruct the software to save the matrix to disk and use another tool (such as https://github.com/riklund/MatrixUtils) to find the eigenvalues, in order to limit RAM usage. Further, note that the matrix needs to be stored in RAM memory at some point, so the number of basis states is ultimately limited by the capacity of your computer.

##Prerequisities
This software depends on a number of external libraries, some of which are bundeled with it and some of which are not. Most importantly, the following libraries must be installed in your PATH: 

* gfortran
* Intel MKL
* libconfig++

The MKL library can be substituted for LAPACK + LAPACKE libraries with minor modifications to the makefile and some headers.


This software also bundles with the following libraries:

* fparser4.5.1 (released under LGPL license)
* ARPACK (released under a modified BSD license)



##Installation
First make sure that the prerequisites are installed. Then clone the repository (recursively, the subrepositories are needed as well) and issue

```bash
./setup.sh
make
```



##Usage
The software is used by simply issuing 

```bash
./Compute
```

As an optional parameter, the flag 

```bash
--configFile (configuration filename)
```

can be used to load a different config file than the default "config.conf"


##Auxiliary programs
Bundled with this software is a large number of additional programs, some more and some less useful, for data analysis and for parallelizing computations. Documentation is included when not self-descriptive.


##Known issues

* Automatic creation of configuration files is not fully implemented. To overcome this issue, use an old config file as template when creating new ones.
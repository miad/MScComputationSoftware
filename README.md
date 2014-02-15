[![Build Status](https://travis-ci.org/riklund/MScComputationSoftware.png)](https://travis-ci.org/riklund/MScComputationSoftware)


#MScComputationSoftware

A program used for computations for my Master's Thesis project. It basically solves the Schrödinger equation in a many-body particle basis expansion. This is done by constructing large matrices and finding their eigenvalues. 

The functionality of the software is best described in my upcoming thesis.

Note that this software is currently under construction, and the results from its computations are currently under verification. Also, a general documentation of the functionality is written in conjunction with my thesis, and currently not available.

##Prerequisities
This software depends on a number of external libraries, some of which are bundeled with it and some of which are not. Most importantly, the following packages needs to be installed in the latest stable Debian versions for the software to even compile:

* liblapack-dev
* liblapacke
* libarpack2-dev
* libblas-dev
* libconfig++-dev

This software also bundles with the following libraries, released under a LGPL license:

* fparser4.5.1



##Installation
First make sure that the prerequisites are installed. Then clone the repository and issue

	  ./setup.sh
	  make




##Usage
The software is used by simply issuing 

	./Compute

As an optional parameter, the flag 

   --configFile (configuration filename)

can be used to load a different config file than the default config.conf


##Auxiliary programs
Bundled with this software is a large number of additional programs useful in data analysis and for parallelizing computations. Documentation is sometimes included. 


##Known issues

* Documentation is not yet complete.
* Automatic creation of configuration files is not fully implemented. Workaround: use a template file.
* Program output not yet verified.
* Repository cleanup pending.
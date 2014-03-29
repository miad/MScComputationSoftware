#!/bin/bash
function Main()
{
	if [ "$1" != "--skip-requirement-check" ]
	then
		CheckRequirements
	fi
	InitRLlib
	InitFParser
	InitArpack
	if [ -f "makefile" ]
	then
		chmod a+w makefile
	fi
	cp makefile.in makefile
	#To prevent accidental modification of makefile instead of makefile.in
	chmod a-w makefile 
	echo "Setup completed without errors."
}

function InitArpack()
{
	cd ARPACK_LIB
	cp ARmake.inc.in ARmake.inc
	chmod u+w ARmake.inc
	sed -i -e "s@home = .*ARPACK_LIB@home = $(pwd)@g" ARmake.inc
	chmod u-w ARmake.inc
	make -j 200 lib
	if [ "$?" -ne "0" ]
	then
		echo "Error in ARPACK compilation."
		exit 15
	fi
	cd ..
}


#Install fparser files in correct directories for build.
function InitFParser()
{
	FPARSERDIR="fparser4.5.1"
	EXTRALIBS="extra_libs"
	echo "Unpacking and preparing fparser..."
	mkdir -p $FPARSERDIR
	if [ ! -f "$EXTRALIBS/$FPARSERDIR.zip" ]
	then
		echo "Could not find $EXTRALIBS/$FPARSERDIR.zip, critical error, exiting."
		exit 9
	fi

	if [ -d "$FPARSERDIR" ]
	then
		read -p "Directory $FPARSERDIR exists, overwrite? (y/n)" -n 1 -r
		echo -en "\n"
		if [[ $REPLY =~ ^[yY]$ ]]
		then
			rm -rf "$FPARSERDIR"
		fi
	fi

	if [ ! -d "$FPARSERDIR" ]
	then
		unzip $EXTRALIBS/$FPARSERDIR.zip -d $FPARSERDIR
	fi

	ln -fs ../$FPARSERDIR/extrasrc include/extrasrc

	FPOBJECTS=(fparser fpoptimizer fpconfig)
	for OBJ in "${FPOBJECTS[@]}"
	do
		touch $FPARSERDIR/$OBJ.cc
		touch $FPARSERDIR/$OBJ.hh
		ln -fs ../$FPARSERDIR/$OBJ.cc src/$OBJ.cc
		ln -fs ../$FPARSERDIR/$OBJ.hh include/$OBJ.hh
	done
}

#Run setup script on RLlib
function InitRLlib()
{
	if [ ! -f "RLlib/setup.sh" ]
	then
		echo "Could not find RLlib/setup.sh, fatal error."
		exit 4
	fi
	cd RLlib
	echo -n "RLlib : "
	./setup.sh
	if [ "$?" -ne "0" ]
	then
		echo "RLlib setup returned nonzero exit code, fatal error."
		exit 5
	fi
	cd ..
}


function CheckRequirements()
{
	#Check if we have a makefile in the current directory. If not, we are obviously not in the project directory.
	if [ ! -f "makefile.in" ]
	then
		echo "Fatal error: makefile.in not found."
		RETURNCODE=1
		exit $RETURNCODE
	fi
	AssurePackageInstalled liblapack-dev
	AssurePackageInstalled libarpack2-dev
	AssurePackageInstalled libblas-dev
	AssurePackageInstalled liblapacke-dev
	AssurePackageInstalled libconfig++-dev
}


#Check if a package is installed and terminate the script if it isn't
function AssurePackageInstalled()
{
	if [[ "$(dpkg -s $1 2>&1 | grep 'Status: install ok installed' | wc -l)" -ne "1" ]]
	then 
		echo "Fatal error: This program requires the package $1. Install it and try again."
		RETURNCODE=2
		exit $RETURNCODE
	fi
}







#Run the main function.
Main "$*"
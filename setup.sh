#!/bin/bash
function Main()
{
	CheckRequirements
	InitRLlib
	if [ -f "makefile" ]
	then
		chmod a+w makefile
	fi
	cp makefile.in makefile
	#To prevent accidental modification of makefile instead of makefile.in
	chmod a-w makefile 
	echo "Setup completed without errors."
}

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
}

function AssurePackageInstalled()
{
	if [[ "$(dpkg -s $1 2>&1 | grep 'Status: install ok installed' | wc -l)" -ne "1" ]]
	then 
		echo "Fatal error: This program requires the package $1. Install it and try again."
		RETURNCODE=2
		exit $RETURNCODE
	fi
}

Main
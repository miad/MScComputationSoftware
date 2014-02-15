#!/bin/bash
function Main()
{
	CheckRequirements
	chmod a+w makefile
	cp makefile.in makefile
	#To prevent accidental modification of makefile instead of makefile.in
	chmod a-w makefile 
	echo "Setup completed without errors."
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
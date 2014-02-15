#!/bin/bash
#Runs a program when more memory than 41 GB becomes available.
while [[ "`ps -u riklund | grep Compute | wc`" -gt "2" ]]
do
	sleep 60
done
$*
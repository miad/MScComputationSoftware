#!/bin/bash
#Runs a program when more memory than 41 GB becomes available.
while [[ "`grep MemFree /proc/meminfo | awk '{print $2}'`" -lt "42991616" ]]
do
	sleep 60
done
$*
#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: ./checkForJobFailures.sh <version>"
    exit 0
fi

version=$1

if [ ! -d $version ]
then
    echo "Subdirectory $version does not exist."
    exit 0
fi

echo "Checking for STL exceptions and CMS exceptions..."
echo "-------------------------------------------------"
grep -r xception $version/LSFJOB_*
grep -r xception $version/*.txt
echo

echo "Checking for segfaults..."
echo "-------------------------------------------------"
grep -r egmentation $version/LSFJOB_*
grep -r iolation $version/LSFJOB_*
grep -r egmentation $version/*.txt
grep -r iolation $version/*.txt
echo

echo "Checking for file copy failures..."
echo "-------------------------------------------------"
grep -r "Error accessing" $version/LSFJOB_*
grep -r quota $version/LSFJOB_*
grep -r "Error accessing" $version/*.txt
grep -r quota $version/*.txt
echo

echo "Checking for job timeouts..."
echo "-------------------------------------------------"
grep -r Cputime $version/LSFJOB_*
grep -r Killed $version/LSFJOB_*
grep -r Cputime $version/*.txt
grep -r Killed $version/*.txt
echo

exit 0

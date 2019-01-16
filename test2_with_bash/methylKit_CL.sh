#!/bin/sh

#  methylKit_CL.sh
#  
#
#  Created by Faraz Ahmed on 1/16/19.
#  


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ "$1" = "help" ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash methylKit_CL.sh <nameOfExperiment> hg38.bed hg38_CPGs.bed"
    echo ""
    echo " Dependencies: hg38.bed, hg38_CPGs.bed, methylKit.R, inputFile.csv and *.cov files"
    echo "--------------------------------------------------------------------------------------"
    echo ""

fi

if [ -e "hg38.bed" ] && [ -e "hg38_CPGs.bed" ] && [ -e "methylKit.R" ] && [ "inputFile.csv" ]; then

        ### Execution of methylKit.R via command line (MACOS terminal)

        if [ ! -z "$1" ] && [ ! -z "$2" ] && [ ! -z "$3" ] && [ "$1" != "help" ]; then

        RScript methylKit.R $1 $2 $3

        fi
else
    echo ""
    echo " Please make sure all of the following files are present in the current directory/folder:"
    echo "  hg38.bed, hg38_CPGs.bed, methylKit.R, inputFile.csv and all bismark coverage files"
    echo ""

fi

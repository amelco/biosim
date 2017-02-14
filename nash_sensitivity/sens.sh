#!/bin/bash
#
###########################################################
##
## Sensitivity analysis
##
## 1. Change of all input parameters at a given percentage 
## 2. Compare to biomass output (the whole period) by Nash and Sutcliff model eficiency 
## assuming the results from the default parameter set as "measured data"

## Begin of script

if [ $# -eq 0 ]; then
  echo -n "Please specify the percentual change of all parameters: "
	read delta
else
	delta=$1
fi
	pdelta=$(echo "$delta/100" | bc -l)

# number of parameters to be changed
np=23

# Beginning to write sensitivity output to a file
#echo "    -=== Sensitivity Summary ===-    "      > result.sns
#echo ""                                          >> result.sns
#echo "Number of parameters to be changed: $np"   >> result.sns
#echo "Percentage of variation: $delta%"          >> result.sns


### changes ###

for (( i=1; i<=$np; i++ )); do
  # changing input file params.in
  case "$i" in
       1) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="beta";;
       2) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="ftr" ;;
       3) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="radext" ;;
       4) awk -v delta=$pdelta '{if(NR==8) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Am" ;;
       5) awk -v delta=$pdelta '{if(NR==9) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Aini" ;;
       6) awk -v delta=$pdelta '{if(NR==10) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Nini" ;;
       7) awk -v delta=$pdelta '{if(NR==11) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Nmax" ;;
       8) awk -v delta=$pdelta '{if(NR==12) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Bini" ;;
       9) awk -v delta=$pdelta '{if(NR==13) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="tau" ;;
      10) awk -v delta=$pdelta '{if(NR==14) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="rdg" ;;
      11) awk -v delta=$pdelta '{if(NR==15) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="rho" ;;
#   [4-9]) awk -v delta=$pdelta -v line=$i '{if(NR==line+4) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Am" ;;
#  1[0-1]) awk -v delta=$pdelta -v line=$i '{if(NR==line+4) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Aini" ;;
      12) awk -v delta=$pdelta '{if(NR==16) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="part_root" ;;
      13) awk -v delta=$pdelta '{if(NR==16) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="part_shoot" ;;
      14) awk -v delta=$pdelta '{if(NR==17) {printf "%12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="bioloss" ;;
      15) awk -v delta=$pdelta '{if(NR==20) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="kcmin" ;;
      16) awk -v delta=$pdelta '{if(NR==20) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="kcfull" ;;
			# Soil input file
      17) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="Or1" ;;
      18) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="Os1" ;;
      19) awk -v delta=$pdelta '{if(NR==3) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="Or2" ;;
      20) awk -v delta=$pdelta '{if(NR==3) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="Os2" ;;
      21) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f,%12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="z1" ;;
      22) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f,%12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="z2" ;;
      23) awk -v delta=$pdelta '{if(NR==7) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' soil.in.default > soil.in ; var="Wini" ;;

     *) echo $i;;
  esac

	# running the model
	./bmg

	# copying results and renaming according to the variable number
	echo "bio_${var}-p${delta}.out"
	filebio=`echo ""`
	fileswb=`echo ""`
	cp bio.out "results/bio_${var}-p${delta}.out"
	cp swb.out "results/swb_${var}-p${delta}.out"


done


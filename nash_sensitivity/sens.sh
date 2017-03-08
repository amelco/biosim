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
  echo ""
	echo "If the first parameter is 0, then the script will compute Nash coefficients."
	echo "OBS.: You have run all percentual variations before compute Nash coefficients."
	echo ""
  #echo -n "Please specify the percentual change of all parameters: "
	#read delta
	echo -n "Type all percentual changes: "
	read -a delta
else
	#delta=$1
	echo "Execute this script without any arguments"
	exit
fi
#pdelta=$(echo "$delta/100" | bc -l)
# number of parameters to be changed
np=23

#if [ $delta -ne 0 ]; then
if [ ${delta[0]} -ne 0 ]; then
  
  ## storing the default result ##
  cp params.in.default params.in
  cp soil.in.default soil.in
 	./bmg &> /dev/null
  #./bmg
  cp bio.out results/bio_default.out
  cp swb.out results/swb_default.out
  #cd results
  #awk '{print $1,$2,$4}' bio_default.out > bdef.tmp
  #cp bdef.tmp biomass.sns
  #cd ..
  
  ### changing parameter values ###
	j=0
for (( c=0; c<${#delta[@]}; c++)); do  
  echo ""
  pdelta=$(echo "${delta[$j]}/100" | bc -l)
	echo "pdelta = ${pdelta}"
	for (( i=1; i<=$np; i++ )); do
    # changing input file params.in
    case "$i" in
         1) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="beta";;
         2) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="ftr" ;;
         3) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="radext" ;;
         4) awk -v delta=$pdelta '{if(NR==8) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Am" ;;
         5) awk -v delta=$pdelta '{if(NR==9) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Aini" ;;
         6) awk -v delta=$pdelta '{if(NR==10) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Nini" ;;
         7) awk -v delta=$pdelta '{if(NR==11) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Nmax" ;;
         8) awk -v delta=$pdelta '{if(NR==12) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="Bini" ;;
         9) awk -v delta=$pdelta '{if(NR==13) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="tau" ;;
        10) awk -v delta=$pdelta '{if(NR==14) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="rdg" ;;
        11) awk -v delta=$pdelta '{if(NR==15) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="rho" ;;
        12) awk -v delta=$pdelta '{if(NR==16) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="part_root" ;;
        13) awk -v delta=$pdelta '{if(NR==16) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="part_shoot" ;;
        14) awk -v delta=$pdelta '{if(NR==17) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' params.in.default > params.in ; var="bioloss" ;;
        15) awk -v delta=$pdelta '{if(NR==20) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' params.in.default > params.in ; var="kcmin" ;;
        16) awk -v delta=$pdelta '{if(NR==20) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' params.in.default > params.in ; var="kcfull" ;;
  			# Soil input file
        17) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="Or1" ;;
        18) awk -v delta=$pdelta '{if(NR==2) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="Os1" ;;
        19) awk -v delta=$pdelta '{if(NR==3) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="Or2" ;;
        20) awk -v delta=$pdelta '{if(NR==3) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="Os2" ;;
        21) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f %12.4f\n", $1+$1*delta,$2} else {print $0}}' soil.in.default > soil.in ; var="z1" ;;
        22) awk -v delta=$pdelta '{if(NR==5) {printf "%12.4f %12.4f\n", $1,$2+$2*delta} else {print $0}}' soil.in.default > soil.in ; var="z2" ;;
        23) awk -v delta=$pdelta '{if(NR==7) {printf "%12.4f\n", $1+$1*delta} else {print $0}}' soil.in.default > soil.in ; var="Wini" ;;
  
       *) echo $i;;
    esac
  
  	# running the model
  	./bmg &> /dev/null
  
  	# copying results and renaming according to the variable number
  	echo "${i} bio_${var}-p${delta[$j]}.out"
  	cp bio.out "results/bio_${var}-p${delta[$j]}.out"
  	cp swb.out "results/swb_${var}-p${delta[$j]}.out"
  #	cd results
  #	awk '{print $4}' bio_${var}-p${delta}.out > out.tmp
  #	cp biomass.sns t.tmp
  #	paste -d' ' t.sns out.tmp > biomass.sns
  
  done
	j=`expr $j+1`
done

  #echo ""
  #echo "Now run statsANDplots.gpt with GNUPLOT" 
  #echo ""

else
  #clear
  echo ""
	echo "-== Computing Nash coeficcients ==-"
	echo ""
	echo -n "Type all percentual changes: "
	read -a p
	cd results

  # preparing nash output file
	echo -e "param" > biomass_nash.sns
	for ((g=0; g<${#p[@]}; g++)); do
	  echo -e "${p[$g]}" >> biomass_nash.sns
	done
	cp biomass_nash.sns n.tmp

  echo "Creating sns files..."
	for f in bio_*-p${p[1]}.out; do
	  # getting name of parameters
	  i=`expr index "$f" -`
		tam=$((${#f} - 10 - ${#p[1]}))
		par_name=${f:4:$tam}
		#echo $par_nameaa

	  # creating tmp file that stores acumulated columns
	  awk 'BEGIN{print "ano","SQD","B_def"} NR>1{print $1, $2, $4}' bio_default.out > tmp

		echo -e -n "\n  Processing parameter ${par_name}..."

		# loop through all variations
		for ((g=0; g<${#p[@]}; g++)); do
      #echo "bio_${par_name}-p${p[$g]}.out"
		  awk -v perc=${p[$g]} 'BEGIN{print "B_"perc} NR>1{print $4}' bio_${par_name}-p${p[$g]}.out > 1.tmp
			paste -d' ' tmp 1.tmp > biomass_${par_name}.sns
			cp biomass_${par_name}.sns tmp
		done
		# removing tmp file
		rm tmp
	  rm 1.tmp

		echo -e -n " \t\tCalculating efficiency..."

	  # Calculating Nash model efficiency e
		echo -e "${par_name}" > 2.tmp
		for ((g=0; g<${#p[@]}; g++)); do
		  avg=`awk 'BEGIN{sum=0;i=0}NR>1{sum=sum+$3;i++}END{printf "%8.2f", sum/i}' biomass_${par_name}.sns`
	    e=`awk -v avg=$avg -v j=$g 'BEGIN{num=0;den=0}{num=num+($(j+3)-$3)**2;den=den+($3-avg)**2}END{print 1-num/den}' biomass_${par_name}.sns`
		  echo -e "$e" >> 2.tmp
		done
		paste -d'\t' n.tmp 2.tmp > biomass_nash.sns
		cp biomass_nash.sns n.tmp
	done
	rm *.tmp
fi

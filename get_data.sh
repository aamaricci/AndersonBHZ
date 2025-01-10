#!/bin/bash
shopt -s expand_aliases
ARGUMENTS=${@}
PATH_TO_ME=$(dirname `realpath ${0}`)
SCRIPT=${0##*/}
LOGFILE=log.${SCRIPT%.sh}
>$LOGFILE
exec >  >(tee -a $LOGFILE)
exec 2> >(tee -a $LOGFILE >&2)
echo "Executing: $SCRIPT from $PATH_TO_ME"

source ~/.aliasrc
source ~/.aliasosxrc
source ~/.functionsrc
source ~/.osx_functionsrc

HERE=`pwd`			# Working directory
WLIST=list_wdis			# List of disorder strenght values

rm z2VSwdis.dat EgapVSwdis.dat PSzPgapVSwdis.dat


while read WDIS
do
    WDIR=WDIS$WDIS
    if [ -d $WDIR ];then
	cd $WDIR

	read ID IDUM < list_idum
	cp IDUM_$IDUM/used.inputABHZ.conf .
	$HOME/.bin/get_data_Abhz_2d 
	echo $WDIR
	#
	echo "get Z2"
	cat IDUM_*/z2.dat|mean  > tmp1
	cat IDUM_*/z2.dat|gmean > tmp2
	cat IDUM_*/z2.dat|sdev  > tmp3
	paste tmp1 tmp3 tmp2 > z2.dat
	#
	echo "get Egap"	
	cat IDUM_*/Egap.dat|mean  > tmp1
	cat IDUM_*/Egap.dat|gmean > tmp2
	cat IDUM_*/Egap.dat|sdev  > tmp3
	paste tmp1 tmp3 tmp2 > Egap.dat
	rm tmp*
	#
	echo "get PSzPgap"
	cat IDUM_*/PSzPgap.dat|mean   > tmp1
	cat IDUM_*/PSzPgap.dat|gmean  > tmp2
	cat IDUM_*/PSzPgap.dat|sdev   > tmp3
	paste tmp1 tmp3 tmp2 > PSzPgap.dat
	#
	rm tmp*
	cd $HERE
    fi
    cat $WDIR/z2.dat|trans $WDIS \#1 \#2 >> z2VSwdis.dat
    cat $WDIR/Egap.dat|trans $WDIS \#1 \#2 >> EgapVSwdis.dat
    cat $WDIR/PSzPgap.dat|trans $WDIS \#1 \#2 >> PSzPgapVSwdis.dat
done < $WLIST

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
WLIST=list_wdis.used		# List of disorder strenght values
ldir |cut -d"S" -f2 |tee $WLIST


alias mgs="python -c 'import numpy as np; import sys;from numpy import *; file=sys.stdin;x=np.asarray(loadtxt(file,unpack=True));print(np.mean(x),np.std(x),np.exp(np.mean(np.log(np.abs(x)))))'"

while read WDIS
do
    WDIR=WDIS$WDIS
    if [ -d $WDIR ];then
	cd $WDIR
	pwd 
	ldir|cut -dU -f2 |tee list_u.used

	rm *mgs*.dat
	
	while read U;do
	    UDIR=U$U
	    if [ -d $UDIR ];then
		cd $UDIR
		ldir |cut -d_ -f2 |tee list_idum.used
		pwd 
		read IDUM < list_idum.used
		cp IDUM_$IDUM/used.inputABHZ.conf .
		$HOME/.bin/get_data_Abhz_2d idumFILE=list_idum.used
		cd ..
	    fi
	done < list_u.used


	echo "get Z2: $WDIS"
	while read U;do
	    UDIR=U$U
	    cat $UDIR/IDUM_*/z2.dat|mgs |awk -v w=$WDIS -v u=$U '{ print u,$1,$2,$3,w}' 
	done < list_u.used |tee z2VSuVSw.mgs

	echo "get Egap: $WDIS"	
	while read U;do
	    UDIR=U$U
	    cat $UDIR/IDUM_*/Egap.dat|mgs |awk -v w=$WDIS -v u=$U '{ print u,$1,$2,$3,w}' 
	done < list_u.used |tee EgapVSuVSw.mgs

	echo "get PSzPgap: $WDIS"
	while read U;do
	    UDIR=U$U
	    cat $UDIR/IDUM_*/PSzPgap.dat|mgs |awk -v w=$WDIS -v u=$U '{ print u,$1,$2,$3,w}' 
	done < list_u.used |tee PSzPgapVSuVSw.mgs
		
	cd $HERE
    fi
done < $WLIST

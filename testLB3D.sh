#!/bin/bash
#=========================================================================
#
# Copyright 1999-2012, Owners retain copyrights to their respective works.
#
# This file is part of lb3d.
#
# lb3d is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# lb3d is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with lb3d. If not, see <http://www.gnu.org/licenses/>.
#
#=========================================================================


# LB3DTEST.sh 1.0 2011/12/29
# Sebastian Schmieschek
#
# This scripts automates simple input output relation tests of the program LB3D
# By Comparing provided refrence HDF5 output to output of a compiled snapshot
# It will loop through the subdirectories <TESTNAME> of ${TESTDIR} and look for a <TESTNAME>.cfg file
#
# The <TESTNAME>.cfg may contain the following data
#
# FLAGS=-DFLAG1 -DFLAG2          	LB3D Compilerflags (see manual) for the test
# INPUTFILE=template.in        		Name of the LB3D input-file (options see manual) used for the test
# OUTFILES=od wd sur dir vel arr	Prefixed of LB3D output (see manual) to be compared
# REFDIR=./reference            	Location of reference output files (HDF5) relative to ${TESTDIR}
# OUTDIR=./output            	Location of output files (HDF5) relative to ${TESTDIR}
# TIMESTEP=100			 	Timestep to compare
#
# A binary lbe.<TESTNAME> is build in ${LB3DDIR}/code/ according to defines.MACHINE using FLAGS
# Stdout and stderr are written to ${TESTDIR}/<TESTNAME>/${STARTTIME}_Compilation.<TESTNAME>.[out|err]
#
# lbe.<TESTNAME> is copied to ${TESTDIR}/<TESTNAME>/ and executed using ${INPUT} as input-file
# Stdout and stderr are written to ${TESTDIR}/<TESTNAME>/${STARTTIME}_Execution.<TESTNAME>.[out|err]
#


# Settings
#####################################################################

MACHINE="LINUX64_TEST" # Which architecture to compile for
LB3DDIR=$(pwd)
LB3DCONFIG='configLB3D.sh'
TESTDIR=${LB3DDIR}/test
HDF5TOOLSPATH="/usr/bin"
MPIPREFIX="mpirun -n 1" # Or leave empty for serial execution. 
FILETYPES="hdf xdr"


# Initialisation 
#####################################################################
# What time is it?
STARTTIME=$(date +%Y%m%d-%H%M%S)
# Where am I
RUNDIR=$(pwd)
# Who am I 
MYNAME=$0

TESTS=''

echo -e "Starting ${MYNAME} in ${RUNDIR} on ${STARTTIME}..\n"
echo "---------------------------------------------------------------------"

# Get a list of Tests
#####################################################################
echo "Checking for tests in ${TESTDIR}.."
echo "---------------------------------------------------------------------"
cd ${TESTDIR}
for i in *
do
	if [ -d $i ]
	then
		if [ -e "./${i}/${i}.cfg" ]
		then
			echo -e "\tFound test \"${i}\" .."
			TESTS="${TESTS} ${i}"
		fi
	fi
done
if [[ ${TESTS} == "" ]] 
then
	echo "Error: Found no test. Exiting."
	exit 1
else
	echo -e "Will run tests: ${TESTS}.\n"
fi


# Check LB3DDIR
#####################################################################
echo "Checking LB3D in ${LB3DDIR}.."
echo "---------------------------------------------------------------------"
if [ -d ${LB3DDIR} ]
then
	if [ -e ${LB3DDIR}/${LB3DCONFIG} ]
	then
		echo -e "Will compile LB3D in ${LB3DDIR}.\n"
	else
		echo "Error: ${LB3DCONFIG} not found in ${LB3DDIR}. Exiting."
		exit 1
	fi
else
	echo "Error: ${LB3DDIR} is not a directory. Exiting."
	exit
fi


# Check H5DIFF
#####################################################################
echo "Checking h5diff in ${HDF5TOOLSPATH}.."
echo "---------------------------------------------------------------------"
if [ -d ${HDF5TOOLSPATH} ] 
then
	if [ -e ${HDF5TOOLSPATH}/h5diff ]
	then
		H5DIFFBIN=${HDF5TOOLSPATH}/h5diff
	else
		echo "Warning: ${HDF5TOOLSPATH}/h5diff does not exist."
	fi
else
	echo "Warning: ${HDF5TOOLSPATH} is not a diretory."
fi
if [ "${H5DIFFBIN}" == "" ]
then
	if [ -e $(which h5diff) ]
	then
		H5DIFFBIN=$(which h5diff)
	else
		echo "Error: h5diff not found. Exiting."
		exit 1
	fi
fi
echo -e "Using ${H5DIFFBIN}.\n"


# Run tests
#####################################################################
for TEST in ${TESTS}
do
if [ "${TEST}" == "template" ]
then 
	echo "Ignoring template."
else 
	cd ${TESTDIR}/${TEST}/
# Read test config
#####################################################################
	echo "Reading ${TESTDIR}/${TEST}/${TEST}.cfg.."
echo "---------------------------------------------------------------------"
	while read LINE
	do
		KEY=$(echo ${LINE} | cut -d '=' -f '1') 	
		VAL=$(echo ${LINE} | cut -d '=' -f '2')
			
		case ${KEY} in
		FLAGS)
			if [ "${VAL}" == "" ]
			then
				echo "Warning: No Compilerflags given."	
			else
				echo "Using Compilerflags ${VAL}."
				FLAGS=${VAL}
			fi
			;;
		REFDIR)
			if [ "${VAL}" == "" ]
			then
				echo "Warning: No location of reference given, using ./ ."
				REFDIR=./
			else
				if [ -d ${TESTDIR}/${TEST}/${VAL} ] 
				then
					echo "Using reference in ${VAL}."
					REFDIR=${VAL}
				else
					echo "Error: ${TESTDIR}/${TEST}/${VAL} is not a directory. Exiting."
					exit 1
				fi
			fi
		 	;;
		OUTDIR)
			if [ "${VAL}" == "" ]
			then
				echo "Warning: No location of output given, using ./ ."
				OUTDIR=./
			else
				if [ -d ${TESTDIR}/${TEST}/${VAL} ] 
				then
					echo "Using output in ${VAL}."
					OUTDIR=${VAL}
				else
					echo "Error: ${TESTDIR}/${TEST}/${VAL} is not a directory. Exiting."
					exit 1
				fi
			fi
			;;
		INPUTFILE)
			if [ "${VAL}" == "" ]
			then    
				echo "Error: No input-file given. Exiting."
				exit
			else    
			        if [ -e ${TESTDIR}/${TEST}/${VAL} ]
				then
					echo "Using input-file ${VAL}."
					INPUTFILE=${VAL}
				else
					echo "Error: Input-file ${TESTDIR}/${TEST}/${VAL} not found. Exiting."
					exit
				fi
			fi      
			;;
		OUTFILES)
			if [ "${VAL}" == "" ]
			then    
				echo "Error: No output to compare given. Exiting."
				exit
			else    
			        echo "Will compare output of ${VAL}."
				OUTFILES=${VAL}
			fi      
			;;
		TIMESTEP)
			if [ "${VAL}" == "" ]
			then    
				echo "Error: No timstep for comparison given. Exiting."
				exit
			else    
			        echo "Will compare at timestep ${VAL}."
				TIMESTEP=${VAL}
			fi      
			;;
		esac
	done < "${TEST}.cfg"
echo "---------------------------------------------------------------------"
# Evaluate input-file 
#####################################################################
	

# Compile LB3D
#####################################################################
	cd ${LB3DDIR} 

	echo "Configuring.."
	./${LB3DCONFIG} CONFIG ${MACHINE} ${FLAGS} 1> ${TESTDIR}/${TEST}/logs/${STARTTIME}_config.${TEST}.out 2> ${TESTDIR}/${TEST}/logs/${STARTTIME}_config.${TEST}.err

echo "---------------------------------------------------------------------"
	echo "Compiling.."
	./${LB3DCONFIG} MKLB3D 1> ${TESTDIR}/${TEST}/logs/${STARTTIME}_compile.${TEST}.out 2> ${TESTDIR}/${TEST}/logs/${STARTTIME}_compile.${TEST}.err
	
echo "---------------------------------------------------------------------"
	echo "Cleaning up ${OUTDIR}.."
	rm ${TESTDIR}/${TEST}/${OUTDIR}/*h5
	rm ${TESTDIR}/${TEST}/${OUTDIR}/*xdr

echo "---------------------------------------------------------------------"

	cp ${LB3DDIR}/code/lb3d ${TESTDIR}/${TEST}/lb3d.${TEST}
	cd ${TESTDIR}/${TEST}/
	for FTYPE in ${FILETYPES}
	do
		echo "Executing ${TEST} writing ${FTYPE}-output.."
		${MPIPREFIX} ./lb3d.${TEST} -f ${INPUTFILE} -d "dump_${FTYPE}.din" 1> ${TESTDIR}/${TEST}/logs/${STARTTIME}_execute.${TEST}.${FTYPE}.out 2> ${TESTDIR}/${TEST}/logs/${STARTTIME}_execute.${TEST}.${FTYPE}.err
		
echo "---------------------------------------------------------------------"
	TIMESTRING=$(printf t%08d ${TIMESTEP})
	FAILED=0
	for OTYPE in ${OUTFILES}
	do
		case ${FTYPE} in
		hdf)
			echo "Comparing ${OTYPE}.."
# Diff'ing HDF
			${H5DIFFBIN} $(ls ./${OUTDIR}/${OTYPE}_*${TIMESTRING}*.h5) $(ls ${REFDIR}/${OTYPE}_*${TIMESTRING}*.h5) > ${TESTDIR}/${TEST}/logs/${STARTTIME}_h5diff.out
		
			if [ "$(cat ${TESTDIR}/${TEST}/logs/${STARTTIME}_h5diff.out)" == "" ]

			then
				echo hdf passed.
			else 
				FAILED=1
				echo hdf failed.
			fi
 			;;
		xdr)
# Diff'ing XDR
			diff $(ls ./${OUTDIR}/${OTYPE}_*${TIMESTRING}*.xdr) $(ls ${REFDIR}/${OTYPE}_*${TIMESTRING}*.xdr) > ${TESTDIR}/${TEST}/logs/${STARTTIME}_xdrdiff.out

                	if [ "$(cat ${TESTDIR}/${TEST}/logs/${STARTTIME}_xdrdiff.out)" == "" ]

                	then
                        	echo xdr passed.
                	else
                        	FAILED=1
                        	echo xdr failed.
                	fi
			;;
		esac
	done
		echo "---------------------------------------------------------------------"
	done
	if [ ${FAILED} -eq 0 ] 
	then
		echo "Test ${TEST} passed."
	else
		echo "Test ${TEST} failed."
	fi
 echo "---------------------------------------------------------------------"
fi
done


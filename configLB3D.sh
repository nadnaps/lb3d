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


# This script sets up a Makefile suitable for the platform the code
# is to be compiled on.
#
# If you try to compile on some unsupported platform foo, then create
# a file called "defines.foo" with the appropriate flags for the compilers,
# then type "./LB3DCONFIG.sh foo".
#

######################################################################
# Variable definition
######################################################################
MYDIR=$(pwd)
MKMFDIR="${MYDIR}/utils/scripts"
ME=$0

DIRS="code/xdrf utils/obstacles" 
CODEDIRS="code "
DOCDIRS=""

# Uncomment to rebuild documentation
#DOCDIRS="doc/doxygen doc/manual"

######################################################################
# Function definition
######################################################################

function getFlags()
{
for CODEDIR in ${CODEDIRS}; do
        cd ${CODEDIR}
        FLAGLIST=$(ls *.F90 *.h| xargs grep '^#e\?l\?ifn\?\(def\)\?' | sed 's/.*#e\?l\?ifn\?\(def\)\? \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq)
        for FLAG in ${FLAGLIST}; do
                echo -D${FLAG}
        done
        cd ..
done
}

function usage()
{
echo Syntax is $0 OPTION PLATFORM COMPILERFLAGS
echo
echo OPTION is one of:
echo CONFIG
echo MAKE
echo CLEAN
echo MKLB3D
echo CLEANLB3D
echo
echo PLATFORM is one of:
for i in defines.* ; do
        new=`echo $i|sed s/defines.//`
	if [ $new != 'TEMPLATE' ] 
	then
        	echo -n $new " "	
	fi
done
echo 
echo
echo  COMPILERFLAGS is one or more of:
getFlags
exit -1

}

######################################################################
# Check commandline arguments 
######################################################################

if  expr $# '<' 1 > /dev/null  ; then
        usage
fi

# Check if option is valid

OPTION=$1
shift

if [ ! "${OPTION}" = "CONFIG" ] ; then
        if [ ! "${OPTION}" = "MAKE" ] ; then
        if [ ! "${OPTION}" = "MKLB3D" ] ; then
	        if [ ! "${OPTION}" = "CLEAN" ]; then
	        if [ ! "${OPTION}" = "CLEANLB3D" ]; then
	        echo OPTION contains an unknown value.
                exit -1
	        fi
		fi
	fi
        fi
fi


######################################################################
# Execute CONFIG
######################################################################

if [ "${OPTION}" = "CONFIG" ] ; then

# Check if platform is valid

PLATFORM=$1
shift
DEFFILE=defines.$PLATFORM

	if [ $PLATFORM == 'TEMPLATE' ]; then
		echo TEMPLATE is not a platform.
	        echo
                usage
	fi
        if [ ! -f $DEFFILE ]; then
	        echo I do not know about platform \"$PLATFORM\".
	        echo
                usage
        fi

# Check if compilerflags are valid

VALIDFLAGS=$(getFlags)
while [[ $1 ]]; do
ISVALID=0
for VALIDFLAG in ${VALIDFLAGS};do
        if [ $1 == ${VALIDFLAG} ]; then
                ISVALID=1
                MAKEFFLAGS="$MAKEFFLAGS $1"
        fi
done
        if [ ${ISVALID} == 0 ]; then
                echo I do not know about compilerflag \"$1\"
                echo 
                usage
        fi
        shift
done

echo 
echo "Creating Makefiles..."
echo

for i in $DIRS $DOCDIRS; do
    echo $i
    cat $MYDIR/$DEFFILE > $MYDIR/$i/Makefile
    echo -e "MAKEFFLAGS=$MAKEFFLAGS" >> $MYDIR/$i/Makefile
    echo "XDRLIB=-L${MYDIR}/code/xdrf" >> $MYDIR/$i/Makefile
    cat $MYDIR/$i/Makefile.STUB >> $MYDIR/$i/Makefile
done

for i in $CODEDIRS ; do
    echo $i
    EXECUTABLE=lb3d
    TEMPLATE=$MYDIR/$i/Makefile.template
    OUTFILE=$MYDIR/$i/Makefile

    # change command line syntax for preprocessor switches if LB3D is
    # to be compiled by an IBM compiler
    if echo $MAKEFFLAGS |grep '[^ ]' > /dev/null
    then
	case $PLATFORM in
	    JUGENE2 | HUYGENS | JUGENE | SP2MPI)
		CODEMAKEFFLAGS=-WF,`echo $MAKEFFLAGS |sed 's/ /,/g'`
		;;
	    *)
		CODEMAKEFFLAGS=$MAKEFFLAGS
		;;
	esac
    fi

    echo "MAKEFFLAGS=$CODEMAKEFFLAGS" > $TEMPLATE
    echo "XDRLIB=-L${MYDIR}/code/xdrf" >> $TEMPLATE
    cat $MYDIR/$DEFFILE >> $TEMPLATE

    cd $MYDIR/$i
    $MKMFDIR/mkmf.pl -m $OUTFILE -t $TEMPLATE -p $EXECUTABLE
    cd $MYDIR
done

fi

######################################################################
# Execute MAKE 
######################################################################

if [ "${OPTION}" = "MAKE" ] ; then
echo
echo "Makefiles copied. Will compile now...."
echo

for i in $DIRS $CODEDIRS $DOCDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make
    cd $MYDIR
done
fi

######################################################################
# Execute MKLB3D
######################################################################

if [ "${OPTION}" = "MKLB3D" ] ; then
echo
echo "Makefiles copied. Will compile lbe now...."
echo

for i in code/xdrf $CODEDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make
    cd $MYDIR
done
fi

######################################################################
# Execute CLEAN
######################################################################

if [ "${OPTION}" = "CLEAN" ] ; then
echo
echo "Cleaning up...."
echo

for i in $DIRS $CODEDIRS $DOCDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make clean
    cd $MYDIR
done
fi

######################################################################
# Execute CLEANLB3D
######################################################################

if [ "${OPTION}" = "CLEANLB3D" ] ; then
echo
echo "Cleaning up...."
echo

for i in code/xdrf $CODEDIRS; do
    echo --------- $i ---------
    cd $MYDIR/$i
    make clean
    cd $MYDIR
done
fi



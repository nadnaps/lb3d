#!/bin/bash

# extract the current revision
REV=`grep LBEREVISION ../../code/lbe_revision.h |sed 's/#.* "//' |sed 's/"//'`

# extract all preprocessor switches---I'm sure there's a more elegant way...
DEFS=`cat ../../code/Makefile.template \
    |grep -v '^#' \
    |sed ':a;N;$!ba;s/\\\n/ /g' \
    |grep -- '-D' \
    |sed 's/[= \n\t,][ \n\t]*/\n/g' \
    |grep '^-D' \
    |sed 's/^-D//' \
    |while read i; do echo -n "$i "; done; echo -ne "\b"`
DEFS=${DEFS% *}

cat Doxyfile.template \
    |sed "s/__PROJECT_NUMBER__/\"r$REV - $DEFS\"/" \
    |sed "s/__PREDEFINED__/$DEFS/" \
    >Doxyfile

doxygen 

if [ $? -ne 0 ]
then
    echo "==========================================================="
    echo "Doxygen failed."
    echo "Look at doxygen-warnings.log or/and rerun with QUIET = NO ."
    echo "==========================================================="
    rm -f doxygen.completed
else
    echo "==============================="
    echo "Doxygen completed successfully."
    echo "==============================="
    touch doxygen.completed
fi

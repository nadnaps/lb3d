#/bin/bash
# This script searches all .F90 and .h files for #ifdef <FLAG> statements and reports all <FLAG>s found.

CODEDIR="../../code"
REPORTFILE="flags.txt"

cd $CODEDIR 

# Write two temp files - one with just <FLAG>, the other with <FLAG> <TOTAL COUNT>
ls *.F90 *.h | xargs grep '^#e\?l\?ifn\?\(def\)\?' | sed 's/.*#e\?l\?ifn\?\(def\)\? \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq > flags.tmp.total
ls *.F90 *.h | xargs grep '^#e\?l\?ifn\?\(def\)\?' | sed 's/.*#e\?l\?ifn\?\(def\)\? \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > flags.tmp.total.cnt

# Mostly the same as above, except we only look inside the lbe_detect_flags() subroutine
cat lb3d_io.F90 | awk '/subroutine lbe_detect_flags()/,/end subroutine lbe_detect_flags/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq > flags.tmp.lbe_detect_flags
cat lb3d_io.F90 | awk '/subroutine lbe_detect_flags()/,/end subroutine lbe_detect_flags/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > flags.tmp.lbe_detect_flags.cnt

# Mostly the same as above, except we only look inside the lbe_more_metadata_phdf5() subroutine
cat lb3d_io_hdf5.F90 | awk '/subroutine lbe_more_metadata_phdf5()/,/end subroutine lbe_more_metadata_phdf5/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq > flags.tmp.lbe_more_metadata_phdf5
cat lb3d_io_hdf5.F90 | awk '/subroutine lbe_more_metadata_phdf5()/,/end subroutine lbe_more_metadata_phdf5/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > flags.tmp.lbe_more_metadata_phdf5.cnt

# Add <FLAG> 0 entries for the missing flags
cat flags.tmp.total flags.tmp.lbe_detect_flags | sort | uniq -u | sed 's/\(.*\)/\1 0/' >> flags.tmp.lbe_detect_flags.cnt
cat flags.tmp.total flags.tmp.lbe_more_metadata_phdf5 | sort | uniq -u | sed 's/\(.*\)/\1 0/' >> flags.tmp.lbe_more_metadata_phdf5.cnt

# Sort the files so the entries exactly match the order of the main file
sort flags.tmp.lbe_detect_flags.cnt > flags.tmp.lbe_detect_flags.cnttotal
sort flags.tmp.lbe_more_metadata_phdf5.cnt > flags.tmp.lbe_more_metadata_phdf5.cnttotal

# Just some text written to the report
echo "# Current compiler flags:" > $REPORTFILE
echo "# NAME total stdout hdf5" >> $REPORTFILE
echo "# " >> $REPORTFILE

# Paste the lines of all files, line by line, and retain only <FLAG> <TOTAL COUNT> <IO> <IO_HDF5>
paste flags.tmp.total.cnt flags.tmp.lbe_detect_flags.cnttotal flags.tmp.lbe_more_metadata_phdf5.cnttotal | awk '{print$1, $2, $4, $6}' >> flagsreport.tmp

# Cleanup
rm -f flags.tmp.*

# Append data to final report
cat flagsreport.tmp >> $REPORTFILE
echo "# " >> $REPORTFILE
echo "# Checking for errors: " >> $REPORTFILE

echo "Creating $REPORTFILE ..."

# Line by line, check if the data is what we want
cat flagsreport.tmp | while read line; do
  varn=`echo $line | awk '{print$1}'`
  vart=`echo $line | awk '{print$2}'`
  ldf=`echo $line | awk '{print$3}'`
  lmmp=`echo $line | awk '{print$4}'`
  # Total number of entries outside the output routines
  vart=$(( $vart - $ldf - $lmmp ))
  # We need exactly one entry for <IO>
  if [ $ldf -ne 1 ]; then
    echo "# WARNING: Flag $varn not reported exactly once in lbe_detect_flags()." | tee -a $REPORTFILE
  fi
  # We need exactly one entry for <IO_HDF5>
  if [ $lmmp -ne 1 ]; then
    echo "# WARNING: Flag $varn not reported exactly once in lbe_more_metadata_phdf5()." | tee -a $REPORTFILE
  fi
  # We need at least one entry outside those two output routines
  if [ $vart -lt 1 ]; then
   echo "# WARNING: Flag $varn is reported but unused." | tee -a $REPORTFILE
  fi
done

# More text output to the report
echo "# " >> $REPORTFILE
echo "# Done!" >> $REPORTFILE

# Cleanup
rm -f flagsreport.tmp

cd ..

mv $CODEDIR/$REPORTFILE .


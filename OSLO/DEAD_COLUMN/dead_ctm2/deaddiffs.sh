#! /local/gnu/bin/bash

#Purpose: Check if the Oslo CTM2 version of DEAD
#is up to date with respect to Charlie's developer's version.
#If there are important changes, you should consider changing the files
#in ctm2_src/DUST_DUST/

#Zender updated CVS version of dust model
dir1=/mn/hox/u1/alfgr/match_dst/dst/
#Oslo CTM2 version of dust model
dir2=/mn/hox/u1/alfgr/ctm2_src/DUST_DUST/

#Go to the ctm2_src/DUST_DUST/ directory
cd $dir1

#Get the list of files (which should equal the DEAD model)
filelist=$(ls)

#Go back to the CTM2 dust directory
cd $dir2

#Print the list of files
echo LIST OF FILES IN ZENDER VERSION AT $dir1
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo $filelist

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo NOW YOU ARE AT $(pwd)

#Go through the list of files
for file in $filelist
do
echo DIFFS IN $dir1$file $file   #Print what files you are looking at
diff $dir1$file ./$file            #Print the differences with DEAD
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
done


#!/bin/bash

##### to deactivate codemods functionality, uncomment nmodlines=0 below #####

cmF90="${1}/codemods.F90"
dcmF90="${1}/dummy_codemods.F90"
cmF90c="${1}/codemods.F90-compare"

master_exits=`git show-ref master`
if [ -n "$master_exists" ]; then
  gitref="master"
else
  gitref=`git rev-parse --abbrev-ref HEAD`
fi

gitmaster=`git rev-parse $gitref 2>/dev/null`
if [ -n "$gitmaster" ]; then
  # compare local code with (local) master
  # (the branch hash might still change, e.g. through squashing)
  echo "Determining code changes relative to git master..." >&2
  echo "If this step takes very long, merging may be indicated:" >&2
  echo "> git checkout master" >&2
  echo "> git pull" >&2
  echo "> git checkout [user]/master-pending" >&2
  echo "> git merge master" >&2
  echo "(backup/commit local changes beforehand)" >&2

  diffcmd="git diff $gitref"
else
  filediff=`diff $dcmF90 $cmF90 2>&1`
  if [ -n "$filediff" ]; then
    # if the code is *newly* unmodified, restore/touch dummy_codemods.F90
    cp $dcmF90 $cmF90
    touch $cmF90
  fi
  exit
fi

nmodlines=`$diffcmd $1 | wc -l`
# only construct codemods.F90 from diffcmd if there are modifications;
# otherwise, use dummy file
nmodlines=0 # uncomment to use dummy_codemods.F90 (no codemods.dat output)
if [ $nmodlines -lt 1 ]; then
  filediff=`diff $dcmF90 $cmF90 2>&1`
  if [ -n "$filediff" ]; then
    # if the code is *newly* unmodified, restore/touch dummy_codemods.F90
    cp $dcmF90 $cmF90
    touch $cmF90
  fi
  exit
fi

# select length of longest modification line for array declaration; need to
# loop to account for string replacements in every line; otherwise, could use
#longestlength=`$diffcmd $1 | awk '{if (length > x) {x=length}}END{print x}'`
longestlength=0
while read -r line; do
  length=`echo "$line" | sed 's/%/+perc;+/g' | sed 's/"/+quot;+/g' | wc -m`
  if [ $length -gt $longestlength ]; then longestlength=$length; fi
done <<< "`$diffcmd $1`"

# store last codemods.F90 for comparison
if [ -f $cmF90 ]; then mv $cmF90 $cmF90c; fi

# construct F90 code: module, variable definitions
printf '!!!!!!!!!!!!!!!!! DO NOT ADD THIS FILE TO VERSION CONTROL !!!!!!!!!!!!!!!!!\n' > $cmF90
printf '! This module is created dynamically by tools/codemods.sh upon compiling. !\n' >> $cmF90
printf '! When running GENE, it outputs codemods.dat containing the git diff info.!\n' >> $cmF90
printf '!!!!!!!!!!!!!!!!! DO NOT ADD THIS FILE TO VERSION CONTROL !!!!!!!!!!!!!!!!!\n' >> $cmF90
printf 'MODULE codemods_mod\n' >> $cmF90
printf '  USE discretization\n' >> $cmF90
printf '  USE par_in\n' >> $cmF90
printf '  USE par_other\n' >> $cmF90
printf '  USE file_io\n' >> $cmF90
printf '  IMPLICIT NONE\n' >> $cmF90
printf '  PUBLIC::write_codemods\nCONTAINS\n' >> $cmF90
printf '  SUBROUTINE write_codemods\n' >> $cmF90
printf '    INTEGER::codemodfile,nlines,lineind\n' >> $cmF90
printf '    CHARACTER(LEN=8)::filestat="replace",filepos="rewind"\n' >> $cmF90
printf '    CHARACTER(LEN=' >> $cmF90
echo $longestlength >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf '),DIMENSION(' >> $cmF90
echo $nmodlines >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf ')::modstring\n' >> $cmF90
printf '    CHARACTER(LEN=150)::codemodfilename\n' >> $cmF90
printf '    nlines=' >> $cmF90
echo $nmodlines >> $cmF90

# extract diffcmd for src/, substitute special characters; note: with awk,
# lines>80 are broken into 60 char segments to not split the newline substitute
#$diffcmd $1 | sed 's/%/+perc;+/g' | sed 's/"/+quot;+/g' | sed ':a;N;$!ba;s#\n#+newl;+" // \&\n      \& "#g' | awk 'BEGIN{p=60}{while(length()>p){if(length()>80){printf substr($0,1,p) "\" // &\n      & \"";$0=substr($0,p+1)}else{printf substr($0,1,80);$0=substr($0,81)}};print}' >> $cmF90
j=1
while read -r line; do
  printf '    modstring(' >> $cmF90
  echo $j >> $cmF90
  # remove newline at end of string
  printf %s "$(< $cmF90)" > $cmF90
  printf ')="' >> $cmF90
  echo "$line" | sed 's/%/+perc;+/g' | sed 's/"/+quot;+/g' | awk \
    'BEGIN{p=60}{while(length()>p){if(length()>80){printf substr($0,1,p) "\" // &\n      & \"";$0=substr($0,p+1)}else{printf substr($0,1,80);$0=substr($0,81)}};print}' >> $cmF90
  # remove newline at end of string
  printf %s "$(< $cmF90)" > $cmF90
  printf '"\n' >> $cmF90
  j=`expr $j + 1`
done <<< "`$diffcmd $1`"

# write to GENE output file codemods.dat
printf '    IF (mype.EQ.0) THEN\n' >> $cmF90
printf '      CALL get_unit_nr(codemodfile)\n' >> $cmF90
printf '      codemodfilename=TRIM(diagdir)//"/codemods"//TRIM(file_extension)\n' >> $cmF90
printf '      OPEN(codemodfile,FILE=codemodfilename,FORM="formatted",STATUS=filestat,POSITION=filepos)\n' >> $cmF90
# substitute quotation marks, newlines back into text
printf '      CALL reinsert_chars(modstring)\n' >> $cmF90

# add header with GENE version info
if [ -n "$gitmaster" ]; then
  printf '      WRITE(codemodfile,"(A)") "Code modifications relative to:"\n' >> $cmF90
  printf '      WRITE(codemodfile,"(A)") "GIT master "//TRIM(git_master)\n' >> $cmF90
  printf '      WRITE(codemodfile,"(A)") "GENE "//release\n' >> $cmF90
else
  printf '      WRITE(codemodfile,"(A)") "Code modifications: SVN revision "//TRIM(svn_rev)//" of GENE "//release\n' >> $cmF90
fi
printf '      DO lineind=1,nlines\n' >> $cmF90
printf '        WRITE(codemodfile,"(A)") TRIM(modstring(lineind))\n' >> $cmF90
printf '      END DO\n' >> $cmF90
printf '      CALL FLUSH(codemodfile)\n' >> $cmF90
printf '      CLOSE(codemodfile)\n' >> $cmF90
printf '    END IF\n' >> $cmF90
printf '  END SUBROUTINE write_codemods\n' >> $cmF90

# subroutine for reinsertion of special characters into string
printf '  SUBROUTINE reinsert_chars(inoutstring)\n' >> $cmF90
printf '    CHARACTER(LEN=' >> $cmF90
echo $longestlength >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf '),DIMENSION(' >> $cmF90
echo $nmodlines >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf ')::inoutstring\n' >> $cmF90
printf '    CHARACTER(LEN=' >> $cmF90
echo $longestlength >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf ')::line\n' >> $cmF90
printf '    INTEGER::ind,nlines,lineind\n' >> $cmF90
printf '    nlines=' >> $cmF90
echo $nmodlines >> $cmF90
# remove newline at end of string
printf %s "$(< $cmF90)" > $cmF90
printf '\n    DO lineind=1,nlines\n' >> $cmF90
printf '      line=inoutstring(lineind)\n' >> $cmF90
printf '      DO\n' >> $cmF90
printf '        ind = INDEX(line,"+quot;+")\n' >> $cmF90
printf '        IF (ind.EQ.0) EXIT\n' >> $cmF90
printf "        line = line(:ind-1)//'\"'//line(ind+7:)\n" >> $cmF90
printf '      END DO\n' >> $cmF90
#printf '      DO\n' >> $cmF90
#printf "        ind = INDEX(inoutstr,'+newl;+')\n" >> $cmF90
#printf '        IF (ind.EQ.0) EXIT\n' >> $cmF90
#printf '        inoutstr = inoutstr(:ind-1)//NEW_LINE(" ")//inoutstr(ind+7:)\n' >> $cmF90
#printf '      END DO\n' >> $cmF90
printf '      DO\n' >> $cmF90
printf '        ind = INDEX(line,"+perc;+")\n' >> $cmF90
printf '        IF (ind.EQ.0) EXIT\n' >> $cmF90
printf "        line = line(:ind-1)//'%%'//line(ind+7:)\n" >> $cmF90
printf '      END DO\n' >> $cmF90
printf '      inoutstring(lineind)=line\n' >> $cmF90
printf '    END DO\n' >> $cmF90
printf '  END SUBROUTINE reinsert_chars\n' >> $cmF90
printf 'END MODULE codemods_mod' >> $cmF90

if [ -f $cmF90c ]; then
  filediff=`diff $cmF90 $cmF90c`
  if [ -z "$filediff" ]; then
    # if there are no changes to codemods.F90, use stored version for older timestamp
    mv $cmF90c $cmF90
  else
    rm $cmF90c
  fi
fi

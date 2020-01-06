#!/bin/bash

# This scripts gets all the GITTA_USER/... branches from the gene-dev.git
# repo on gitta and pushes them to the new GITLAB_USER fork of GITLAB_REPO_NAME.
# If GITTA_USER and GITLAB_USER differ, both have to be given on the command line.

# The GITTA_USER/master and the GITTA_USER/master-pending branch are not copied,
# because the name would then collide with the overall master branch of the fork.

# All "-pending" suffixes are deleted and only the basenames are used as new branch
# names.

startdir=$(pwd)
LOCAL_CLONE=gene_gitlab_user_clone
GITLAB_REPO_NAME=gene-dev

if [ $# -eq 1 ]; then
    GITTA_USER=$1
    GITLAB_USER=$GITTA_USER
else
    if [ $# -eq 2 ]; then
	GITTA_USER=$1
	GITLAB_USER=$2
    else
	echo "Usage: get_user_branches.sh GITTA_USERID [GITLAB_USERID]"
	exit 1
    fi
fi

skip=no

if [ "$skip" == "no" ]; then
    URL=git@gitlab.mpcdf.mpg.de:${GITLAB_USER}/${GITLAB_REPO_NAME}.git

    PREFIX=${LOCAL_CLONE}
    while [ -d ${LOCAL_CLONE} ]; do
	counter=$(( $counter + 1 ))
	LOCAL_CLONE=$(printf "${PREFIX}_%2.2d" $counter)
    done

    # first check if the remote repository already exists
    git ls-remote -h $URL &>/dev/null
    if [ $? -ne 0 ]; then
	echo "--------------------------------------------------------"
	echo "-- Please fork the GENE/gene-dev repository first     --"
	echo "-- via the webinterface https://gitlab.mpcdf.mpg.de   --"
	echo "-- Then restart this script to get the user branches. --"
	echo "--------------------------------------------------------"
	exit 1
    fi
    git clone $URL -o gitlab ${LOCAL_CLONE}
    cd ${LOCAL_CLONE}

    # just to make it difficult to push something accidentally
    git config push.default nothing

    # add gitta as remote and change default fetch behaviour
    git remote add gitta GENE@gitta.rzg.mpg.de:gene-dev.git
    git config remote.gitta.fetch +refs/heads/${GITTA_USER}/*:refs/remotes/gitta/*

    # fetch all user branches from gitta
    git fetch gitta +refs/heads/${GITTA_USER}/*:refs/remotes/gitta/*
fi

git checkout master

# get a list of all branches of the gitlab remote
gitlab_branchlist=$(git branch --list -a *gitlab*)

# now push all user branches
# first get a list of all branches on gitta
branchlist=$(git branch --list -a *gitta*)
for b in $branchlist; do
    branchname=$(basename $b)
    echo ""
    echo "branch : $b -> basename = $branchname"
    if [ "$branchname" != "master" ] && [ "$branchname" != "master-pending" ]; then
	# we skip the user/master branch to not collide with the overall master branch

	# from the pending branches we delete the pending suffix
	bname_no_pending=$(sed 's/^\([a-zA-Z_0-9.-]*\)-pending$/\1/' <<< $branchname)

	# check if there is a name collision with the branches in the gitlab-main remote
	skip_this_branch=no
	for gb in $gitlab_branchlist; do
	    #echo "gb = $gb"
	    # exclude the lines which contain HEAD or -> from the test
	    if [ `grep "\->" <<< $gb` ]; then
		continue
	    fi
	    gitlab_branchname=$(basename $gb)
	    if [ "$gitlab_branchname" == "HEAD" ]; then
		continue
	    fi
	    #echo "Testing for name conflict with $gitlab_branchname"
	    if [ "$bname_no_pending" == "$gitlab_branchname" ]; then
		# Name conflict: Now check if there are differences
		git diff --quiet --exit-code $b $gb
		if [ $? -ne 0 ]; then
		    # check if the two branches are in a fast-forward relation to each other
		    rl_gitta_gitlab=$(git rev-list ${b}..$gb)
		    rl_gitlab_gitta=$(git rev-list ${gb}..${b})
		    if [ -z "${rl_gitta_gitlab}" ] && [ -n "${rl_gitlab_gitta}" ]; then
			# gitlab is an ancestor of gitta, hence push the newer version
			echo "Branch $gb on gitlab is an ancestor of branch $g on gitta"
			echo "hence we merge the gitta branch into the gitlab branch"
			echo "and push it to gitlab."
			git checkout ${gitlab_branchname}
			git merge $b
			git push gitlab ${gitlab_branchname}:refs/heads/${gitlab_branchname}
		    elif [ -n "${rl_gitta_gitlab}" ] && [ -z "${rl_gitlab_gitta}" ]; then
			# gitta is an ancestor of gitlab, just skip as the version
			# on gitlab is already newer and contains all of the gitta
			# revisions
			skip_this_branch=yes
		    else
			echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
			echo "Name conflict for branch $bname_no_pending"
			echo "Skipping this branch, you have to get it manually."
			echo ""
			echo "The easiest way is to rename the branch with"
			echo "  git push gitlab ${b}:refs/heads/<new name of the branch>"
			echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
			skip_this_branch=yes
		    fi
		fi
		break
	    fi
	done
	if [ "$skip_this_branch" == "no" ]; then
	    echo "git push gitlab +$b:refs/heads/${bname_no_pending}"
	    git push gitlab +${b}:refs/heads/${bname_no_pending}
	fi
    fi
done

cd $startdir
rm -Rf ${LOCAL_CLONE}

#!/bin/bash

echo "Number of arguments is $#"
if [ "$#" -eq 1 ]; then
    GITTA_USERNAME=$1
    GITLAB_USERNAME=${GITTA_USERNAME}
elif [ "$#" -eq 2 ]; then
    GITTA_USERNAME=$1
    GITLAB_USERNAME=$2
else
    echo "You have to pass your gitta username [and your gitlab username]."
    exit 1
fi

echo "Using GITTA_USERNAME=${GITTA_USERNAME} and GITLAB_USERNAME=${GITLAB_USERNAME}."

delete_old_remote_branches=yes

# check for existing remotes
remotes=$(git remote)
gitlab_exist=0
gitlab_main_exist=0
for r in ${remotes}; do
    if [ "$r" == "gitlab" ]; then
	gitlab_exist=1
    fi
    if [ "$r" == "gitlab-main" ]; then
	gitlab_main_exist=1
    fi
done

# Add the two new remotes
if [ ${gitlab_exist} -eq 0 ]; then
    echo "Adding gitlab as remote"
    git remote add gitlab git@gitlab.mpcdf.mpg.de:${GITLAB_USERNAME}/gene-dev.git
fi

if [ ${gitlab_main_exist} -eq 0 ]; then
    echo "Adding gitlab-main as remote"
    git remote add gitlab-main git@gitlab.mpcdf.mpg.de:GENE/gene-dev.git
fi

#git config remote.gitlab.push refs/heads/master:refs/heads/master

git fetch gitlab
git fetch gitlab-main

#branchlist=$(git branch --list ${GITTA_USERNAME}/*)
branchlist=$(git branch | cut -c 3-)

#echo "branchlist=$branchlist"
echo "++++++ Loop over all local branches +++++++++++++++++++++++++++++++++++"
for b in $branchlist; do
    echo "Working on branch $b"
    branchname=$(basename $b)
    username=$(dirname $b)
    # from the pending branches we delete the pending suffix
    bname_no_pending=$(sed 's/^\([a-zA-Z_0-9.-]*\)-pending$/\1/' <<< $branchname)
    #echo "$b -> branchname = $branchname, user=$username"

    # We skip all local branches from other users
    if [ "$username" == "${GITTA_USERNAME}" ]; then
	#echo "branch : $b -> basename = $branchname"
	# we skip the user/master branch to not collide with the overall master branch
	if [ "$branchname" != "master" ] && [ "$branchname" != "master-pending" ]; then
	    # First, rename the branch without username and without -pending
	    if [ "${bname_no_pending}" == "titan" ]; then
		# titan-pending has to be moved to gpu branch
		echo "Rename branch $b to gpu"
		target_bname="gpu"
	    else
		target_bname=${bname_no_pending}
	    fi
	    git branch -m $b ${target_bname}
	    git config branch.${target_bname}.remote gitlab
	    git branch --set-upstream-to=gitlab/${target_bname} ${target_bname}
	    git config branch.${target_bname}.merge refs/heads/${target_bname}
	fi
    else
	# username of branch is not the actual username
	if [ -n "$username" ] && [ "$username" != "." ]; then
	    # username of branch is another user's gitta username
	    echo "---------------------------------------------------------------------------"
	    echo "  Branch $b is from a different user. To get this branch again from gitlab,"
	    echo "  you have to first make sure the username $username is also valid on gitlab"
	    echo "  and then setup a new remote:"
	    echo "   git remote add gitlab-$username git@gitlab.mpcdf.mpg.de:$username/gene-dev"
	    echo "  The repository of the other user must also be accessible for you."
	    echo ""
	    echo "  Then rename the branch to get rid of the -pending suffix:"
	    echo "   git mv $b $username/${bname_no_pending}"
	    echo ""
	    echo "  Finally, you have to set the upstream of your branch to the new repository with:"
	    echo "   git branch --set-upstream-to=gitlab-$username/${bname_no_pending} $username/${bname_no_pending}"
	    echo "---------------------------------------------------------------------------"
	else
	    # username of branch is empty or .
	    if [ "$branchname" != "master" ] && \
		[ "$branchname" != "master-pending" ] && \
		[ "$branchname" != "master_temp" ]; then
		branch_remote=$(git config --get branch.$b.remote)
		if [ "x$?" == "x0" ]; then
		    # we have a remote of the branch
		    if [ "${branch_remote}" == "gitta" ]; then
			echo "change remote from gitta to gitlab"
			git config branch.$b.remote gitlab
			branch_merge=$(git config --get branch.$b.merge)
			echo "make sure the upstream branch $branch_merge does exist on gitlab"
		    else
			echo "remote is not gitta (it is $branch_remote), hence no changes applied."
		    fi
		else
		    echo "Skipping the branch $b, as it does not have an upstream."
		fi
	    fi
	fi
    fi
    echo ""
done

# Finally, redirect master
# read from gitlab-main
# this seems to be a bit confusing
#git config branch.master.remote gitlab-main
git config branch.master.remote gitlab
# write to gitlab
#git config branch.master.pushRemote gitlab
git config branch.master.merge refs/heads/master
git config push.default current

# what to do with the USER/master-pending branch
# we merge it into the master branch and delete it
git checkout master
if [ `git rev-parse --verify --quiet ${GITTA_USERNAME}/master-pending` ]; then
    git merge ${GITTA_USERNAME}/master-pending
    git branch -d ${GITTA_USERNAME}/master-pending
fi


if [ "$delete_old_remote_branches" == "yes" ]; then
    branchlist=$(git branch -r --list | grep "gitta/" | cut -c 3-)
    for b in $branchlist; do
	if [ ! `grep HEAD <<<$b` ] && [ "$b" != "->" ]; then
	    git branch -d -r $b
	fi
    done
fi

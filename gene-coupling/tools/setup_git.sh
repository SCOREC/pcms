#!/bin/bash
#configures local git working copy for usage with gene-dev.git repository
for arg in $@; do
    if [ "$arg" == "--replace-all" ]; then
	replace_all=True
	echo ""
	echo "  All hooks and config entries will be replaced by this script."
	echo "  This is safe, if you did not install manually any hooks or"
	echo "  changed your config file."
	echo ""
    else
	replace_all=False
    fi
done

pushdefault=`git config --global --get push.default`
if [ -z "$pushdefault" ]; then
    git config --global push.default upstream
else
    if [ "${replace_all}" == "True" ];then
	echo "push.default is kept: $pushdefault"
    fi
fi

# set the pager to a less which understands color codings
value=`git config --global --get core.pager`
if [ -z "$value" ]; then
    git config --global core.pager 'less -FRMX'
else
    if [ "${replace_all}" == "True" ];then
	echo "core.pager is kept: $value"
    fi
fi

# ascii art for the log command
value=`git config --get alias.lg`
if [ -z "$value" ]; then
    git config --global alias.lg "log --graph --pretty=tformat:'%Cred%h%Creset -%C(yellow)%d%Creset %aN - %s %Cgreen(%cr)%Creset' --abbrev-commit --date=relative"
fi

# an alias for getting the root path of the repository
value=`git config --global --get alias.root`
if [ -z "$value" ]; then
    git config --global --add alias.root '!pwd'
fi

# ask to set the editor
value=`git config --global --get core.editor`
if [ -z "$value" ]; then
    echo "Please set your favorite editor for editing commit messages with the command"
    echo " git config --global core.editor <your editor>"
fi

# -------------------------------------------
# here comes the repo specific configuration

# activate the default pre-commit hook
# which checks for whitespace problems
basedir=`git root`
if [ ! -e "$basedir/.git/hooks/pre-commit" ] || [ "${replace_all}" == "True" ]; then
    if [ -e "$basedir/.git/hooks/pre-commit.sample" ]; then
	cp $basedir/.git/hooks/pre-commit.sample $basedir/.git/hooks/pre-commit
    else
	cp $basedir/tools/git_hooks/pre-commit $basedir/.git/hooks/pre-commit
    fi
fi

# hook for preparing the commit message
if [ ! -e "$basedir/.git/hooks/prepare-commit-msg" ] || [ "${replace_all}" == "True" ]; then
    cp $basedir/tools/git_hooks/prepare-commit-msg $basedir/.git/hooks/prepare-commit-msg
else
    echo "There is already a prepare-commit-msg hook installed. No changes are applied."
fi

# copy the keyword testing script into the commit-msg hook
if [ ! -e "$basedir/.git/hooks/commit-msg" ] || [ "${replace_all}" == "True" ]; then
    cp $basedir/tools/git_hooks/commit-msg $basedir/.git/hooks/commit-msg
    cp $basedir/tools/git_hooks/check_message $basedir/.git/hooks/check_message
else
    echo "There is already a commit-msg hook installed. No changes are applied."
fi


cd $basedir
remote=gitlab
remote_url=$(git config --get remote.${remote}.url)
if [ $? -ne 0 ]; then
    #perhaps the repo was cloned with remote named origin
    remote_url=$(git config --get remote.origin.url)
    if [ $? -eq 0 ]; then
	remote=origin
	echo "The usual advice is to use the gitlab repository and name the linked remote"
	echo " with the name gitlab. This can be achieved by cloning with the -o gitlab"
	echo " command line option."
    else
	echo "Cannot determine the name of the remote. Neither origin nor $remote."
	exit 1
    fi
fi

# add the main gitlab repository as the gitlab-main remote
remotes=$(git remote)
gitlab_main_exists=0
for r in $remotes; do
    if [ "$r" == "gitlab-main" ]; then
	gitlab_main_exists=1
    fi
    #echo $r
done
if [ ${gitlab_main_exists} -ne 1 ]; then
    git remote add gitlab-main git@gitlab.mpcdf.mpg.de:GENE/gene-dev.git
fi
git fetch gitlab-main

# get everything from gitlab-main
git fetch gitlab-main

# set the read and write behaviour of the master branch
# read from gitlab-main
git config branch.master.remote gitlab-main
# write to gitlab
git config branch.master.pushRemote gitlab
git config branch.master.merge refs/heads/master
git config push.default current

git checkout master

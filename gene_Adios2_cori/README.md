<!--
PLEASE NOTE: The instructions below are only relevant for the development
repository and may be safely ignored by gitta users.
Information on installing and running GENE can be found in
./INSTALL
./doc/gene.pdf (*)
./doc/tutorial.pdf (*)
(*) requires 'gmake doc' or 'make doc' first
-->
How to use this repository?
===========================

Create a fork (one-time action!)
-------------------------

If new to this repository, you first need to create a [fork](/../forks/new)
(one time only). Please visit the website of the repository and click
**[here](/../forks/new)** or find the corresponding button on top.

*Only one fork per user* is permitted.
By default, a fork only contains the `master` branch,
the `tags` (releases) and some long-term development projects involving more
than one developer. Instructions on how to
[add new branches](README.md#createswitch-branches)
or
[import branches from gitta](README.md#import-your-gitta-user-branches-one-time-action)
can be found below.

Upload your public SSH key(s)
-----------------------------

For `git+ssh` features, you need to upload your public SSH key(s)
**[here](/../../../profile/keys)**. Please visit the following
**[help page](/../../../help/ssh/README)** for further instructions on generating
and managing SSH keys.

Create a local clone
-------------------

A local clone of your own fork can be obtained via

`git clone git@gitlab.mpcdf.mpg.de:$GITLAB_USERID/gene-dev.git -o gitlab <your local clone>`

**Please run the `./tools/setup_git.sh` script** hereafter to set up some hooks,
remotes, aliases and default settings for working with git and the `gene-dev`
repository. This script in particular adds a remote `gitlab-main` which points
to the main development repository `GENE/gene-dev.git` on `gitlab.mpcdf.mpg.de`.

Note that changes to `gitlab-main` can only be triggered via
**[Merge Requests](/../merge_requests)**
([see instructions below](README.md#upload-and-incorporate-changes-into-the-main-development-merge-request)).

Create/switch branches
----------------------

For working on one of the development branches, a local branch has to be created
first via

`git checkout <local branchname>`

If `<local branchname>` can unambiguously mapped to a remote tracking branch,
the latter is automatically set as upstream. If none or more than one remote
branch are found with this name, the tracking branch has to be provided
explicitly via

`git checkout -b <local branchname> --track=<remote>/<remote branchname>`

This is, for instance, always necessary for the `gitlab-main` branches with
more than one developer (`minigene`, `gpu`, etc) as they exist both on
`gitlab-main` and in your `gitlab` fork.

Synchronize your local clone
------------------------------
First, download the actual state of `gitlab-main` into the local remote
tracking branches `refs/remotes/gitlab-main/*` via

`git fetch gitlab-main`

which should create no harm.
Hereafter, you have three different possibilities on how to concile your
changes with the changes from master. But before continuing, switch to your
local branch with

`git checkout <your local branch>`


### Fast-forward update
This is the easiest way but is only possible if only you did changes in
your working copy and there are no changes on the public gitlab-main/master.
Then you can just make

`git merge gitlab-main/<remote branch>`

### Rebasing
If there are changes on the master branch of gitlab-main and you have some
commits locally, you can first try if a rebasing works. Rebasing means that
your local commits are replayed on top of the gitlab-main/master and leads
to a nice linear history. Your local commits are rewritten (Hashes will change
and also the commit date), so never do a rebase on commits that you already
published in a way. You cannot disturb the gitlab-main/master as it is
protected, so you cannot do any harm to this official branch, but this is not
true for the public development branches like gene3D or refactoring.
How works rebasing:

`git rebase gitlab-main/<remote branch>`

If there are conflicts, you have to resolve them and

`git rebase --continue`

### Merging
Sometimes you have many commits and each of them changes lot of things, then a
rebase will become very inconvenient as you have to resolve the conflicts for
each commit. Then it might be a better idea to do a merge, where you have to
resolve the conflicts only once.
It is done with:

`git merge gitlab-main/<remote branch>`


Upload and incorporate changes into the main development (Merge request)
------------------------------------------------------------------

First, commit your local changes as usual via

`git commit -m "KEYWORD: meaningful message" <files>`

and make sure to synchronize your local and remote branches as explained
above. Then run

`git push`

In order to incorporate your changes into the main development branch,
you need to issue a **[Merge request](/../merge_requests)**
which is currently only possible through the web interface.

**[Click here](/../merge_requests)**
or the corresponding button on top and select the branch of your fork that
shall be merged into `gitlab-main`.

One of the senior developers will then take care and initiate the `merge`
(perhaps after some discussion, which is now facilitated via the gitlab
web interface).

Please make sure to first incorporate all actual changes from
`gitlab-main` into the local branch. The **Merge request** will be rejected
otherwise.

Further information
-------------------

Please check the **[Wiki](/../wikis/home)** for further information.



Moving from gitta to gitlab
===========================

The following section is only relevant to developers with previous
write access to `gitta`.

Import your gitta USER branches (one-time action!)
--------------------------------------------------

On `gitta`, every developer had a `USERID/master`, a `USERID/master-pending`
and typically a number of `USERID/...-pending` branches. In order to import
these branches to `gitlab`, you may run

`./tools/get_user_branches.sh GITTA_USERID GITLAB_USERID`

or

`./tools/get_user_branches.sh GITTA_USERID`

if your USERIDs are identical in an up-to-date clone.

The script will then create a temporary clone of your fork,
add `gitta` as a remote, and fetch all your `gitta` branches associated with
your USERID. Hereafter, all of your USERID branches will be pushed to
your `gitlab` fork. The following rules apply:

1. `USERID/master` and `USERID/master-pending` are ignored and __NOT__ copied,
to avoid conflicts with the master branch from the gitlab-main repository.

2. All `USERID/XXX-pending` branches are pushed without the -pending suffix.
This suffix is not anymore needed and only the basename is used in the new
forked repository.

3. In case of conflicting branch names between your `gitlab` fork and
a `USERID/XXX-pending` branch on gitta (most likely for `minigene-pending`,
`gpu-pending` and `gene3D-pending`) a fast-forward push from the `gitta` branch
to the `gitlab` branch will be attempted. If not successful, a warning will be
displayed containing instructions on how to proceed.


Modifying an existing local gitta clone to work with gitlab
-----------------------------------------------------------
First make sure all of your local changes are committed to your local
repository and your local clone is up-to-date with `gitta`.

Then run

`./tools/switch_from_gitta_to_gitlab.sh GITTA_USERID [GITLAB_USERID]`

in your local `gitta` clone. This script adds `gitlab-main` as remote and
fetches all the branches into the local remote-tracking branches.

It furthermore adapts your local branches to the new naming scheme and
sets the upstream to the respective `gitlab` branches, if possible. It
is important to check the output of this script and the local
branches carefully to ensure everything worked as expected.

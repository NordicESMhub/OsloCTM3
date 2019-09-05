# Oslo CTM3
A global chemical transport model developed at MetOs

version: See MODEL_VERSION in cmn_size.F90
======================================================================
When referring to Oslo CTM3, use "Oslo CTM3 vX.Y".
======================================================================

User manual
The user manual is located in the file manual_osloctm3.pdf.
Or can be found on git-hub.

+++++++++++
Quick start
+++++++++++

1. If you are NOT planning to change code, aka development, "fork" from master on git-hub and clone your fork from git-hub:

> git clone git@github.com:<github_username>/OsloCTM3.git

2. Otherwise create your own development branch from master (on git-hub) and give it a descriptive (not too long) name. Clone that branch.

This can also be done in a terminal. Log-in to HPC machine
Clone master from git-hub:

> git clone github.com:NordicESMhub/OsloCTM3.git

Make your local branch and give it a descriptive name, e.g. username and project abbreviation:

> git branch <branch_name>

Push your branch to remote (git-hub):

> git push -u origin <branch_name>

+++
Before compiling
+++

If you are running on abel, you should do

> module purge

to unload all automatically load modules and do

> module load netcdf.intel/4.3.3.1

This will load all necessary dependencies.

Export the path of the working copy of YOUR Oslo CTM3:

> export CTM3_DIR=$HOME/OsloCTM3

and do the same with your work directory:

> export WORK=/work/users/<YOUR_USER_NAME>

Export your notur project number:

> export PROJECT=nnXXXXk

Set an alias for the job queue on abel:

> alias squeue='squeue -lA ${PROJECT}'

Tip: 
Since these steps have to be repeated every time you log in,
it is wise to put the commands into a bash script (e.g. setpaths) 
that you can "source" after each log-in.

+++
Compile
+++

> cd $CTM3_DIR

> make -j

+++
Before running
+++

Make sure you are member of the group "cic-hpc" and have access to

> ls /work/projects/cicero/ctm_input/

Then export the path to forcing etc.:

> export CTM3_INPUT=/work/projects/cicero/ctm_input/Indata_CTM3

Make sure you are member of the group "gf-ozone" and have access to

> ls /projects/researchers/researchers01/sfalk/input/ctm_input/

Then export the path:

> export CTM_USR_INPUT=/projects/researchers/researchers01/sfalk/input/ctm_input/

+++
To run the example ("c3run_example.job") in c3run
+++

Open c3run_example.job in an editor of your choice (e.g., emacs, vi, gedit)
and get familiar with it. 
You do not have to change much to run the example.
In line 8, replace $PROJECT with your project number:

> #SBATCH --account=$PROJECT

Save this change but remember do not commit it to git-hub!
Now change to your work directory:

> cd $WORK

Create a new directory - it should have the same name as is given in line 5 of c3run_example.job - and change to it:

> mkdir C3RUN_example

> cd C3RUN_example

Finally send the model to the batch system of abel:

> sbatch $CTM3_DIR/c3run/c3run_example.job

Wait for about 15-30 min.

+++
During the run
+++

You can check the progress with:

> squeue

+++
After the run
+++

Look at the results.
For a quick few on results, you can use ncview on abel

> module load ncview

or panoply (https://www.giss.nasa.gov/tools/panoply/download/).


Have fun!


+++
HOW TO git
+++

1. Fork
Sometimes it will be necessary to "sync" your fork to the original master.
Details about that process can be found in github help (https://help.github.com/en/articles/syncing-a-fork).
Here we briefly summarize the steps:

Change into the working directory of your OsloCTM3

> cd $CTM3_DIR

Fetch changes from the "original" master:

> git fetch upstream

Switch to your "local" master:

> git checkout master

Merge the "original" master into your "local" master:

> git merge upstream/master

2. Branch
Sometimes it will be necessary to "merge" changes in the master branch beck into your own branch. Here we will give a brief summary of the steps:

Change into the working directory of your OsloCTM3

> cd $CTM3_DIR

Change to "master" branch:

> git checkout master

Get the changes from remote:

> git pull

Change back into your own branch:

> git checkout <branch_name>

Merge the changes into your branch:

> git merge master

*Unfortunately, the automatic merging may cause "conflicts" in case your own changes clash with the changes in master. Which means you have to go through the files that cause this manually to "resolve" the "conflicts".*

After merging is done you can push the updated branch to remote:

> git push

# Oslo CTM3
A global chemical transport model developed at MetOs

~~~
version: See MODEL_VERSION in cmn_size.F90
When referring to Oslo CTM3, use "Oslo CTM3 vX.Y".
~~~

## User manual
The user manual is located in the file manual_osloctm3.pdf.
Or can be found on git-hub.

---
## Quick start

Log in to ABEL/SAGA and create a directory:
~~~
> mkdir models
~~~
You are going to clone the Oslo CTM3 into that directory.
Change into the directory:
~~~
> cd models
~~~

> WARNING: The code changes and configuration for SAGA will only be merged into the master after December 1!!

### 1. Fork

If you are NOT planning to change code, aka development, 
"fork" from master on git-hub (online in your browser). 
Clone your fork:
~~~
> git clone git@github.com:<github_username>/OsloCTM3.git
~~~
Configure the remote repository for your fork.
List the current configured remote repository:
~~~
> git remote -v
~~~
It will be most likely empty.
Specify a new remote upstream repository that will be synced with the fork:
~~~
> git remote add upstream https://github.com/NordicESMhub/OsloCTM3.git
~~~
Check configuration again:
~~~
> git remote -v
~~~
### 2. Branch

Otherwise create your own development branch from master (on git-hub) and give it a descriptive (not too long) name, e.g. use your git-hub username and an abbreviation for your project. Clone that branch:
~~~
> git clone github.com:NordicESMhub/<username_project>.git
~~~
Branching can also be done in terminal. Log-in to HPC machine (abel).
Clone master from git-hub:
~~~
> git clone github.com:NordicESMhub/OsloCTM3.git
~~~
Make your local branch and give it a descriptive name, e.g. username and project abbreviation:
~~~
> git branch <username_project>
~~~
Push your branch to remote (git-hub):
~~~
> git push -u origin <username_project>
~~~
## Before compiling

You have to set up your environment. 
You should always unload all automatically loaded modules first:
~~~
> module purge
~~~
Export the path of your Oslo CTM3 working directory <username_project>:
~~~
> export CTM3_ROOT=${HOME}/OsloCTM3/<username_project>
~~~
Export your notur project number:
~~~
> export PROJECT=nnXXXXk
~~~
Set an alias for the job queue on abel:
~~~
> alias squeue='squeue -lA ${PROJECT}'
~~~

### 1. ABEL
Load all necessary dependencies:
~~~
> module load netcdf.intel/4.3.3.1
~~~
Export your work directory:
~~~
> export WORK=/work/users/${USER}
~~~

### 2. SAGA
Load all necessary dependencies:
~~~
> module load netCDF-Fortran/4.4.4-intel-2018b
~~~
Export your work directory:
~~~
> export WORK=/cluster/work/users/${USER}
~~~
Export the input directory:
~~~
> export CICERO=/cluster/projects/nn9188k/OsloCTM3
~~~
You have to set:
~~~
> export NETCDF_ROOT=$EBROOTNETCDFMINFORTRAN
~~~

> Tip: 
> Since these steps have to be repeated every time you log in,
> it is wise to put the commands into a bash script (e.g. setpaths) 
> that you can "source" after each log-in (example further below).

## Compile

Change to the working director of your Oslo CTM3:
~~~
> cd $CTM3_ROOT
~~~
Compile the model (in parallel, using all available resources):
~~~
> make -j
~~~
In case the compilation was successful, you will find the Oslo CTM3 executable
"osloctm3" in your $CTM3_ROOT.

## Before running

### 1. ABEL
Make sure you are member of the group "cic-hpc" and have access to
~~~
> ls /work/projects/cicero/ctm_input/
~~~
Then export the path to forcing etc.:
~~~
> export CTM3_INPUT=/work/projects/cicero/ctm_input/Indata_CTM3
~~~
Make sure you are member of the group "gf-ozone" and have access to
~~~
> ls /projects/researchers/researchers01/sfalk/input/ctm_input/
~~~
Then export the path:
~~~
> export CTM_USR_INPUT=/projects/researchers/researchers01/sfalk/input/ctm_input/
~~~

### 2. SAGA
Make sure you are member of the group "nn9188k" and have access to
~~~
ls $CICERO
~~~
Then export the path to forcing etc.:
~~~
> export CTM3_INPUT=${CICERO}/Indata_CTM3
~~~
You can use your own forcings by setting:
~~~
> export CTM_USR_INPUT=<your_forcing_files>
~~~

## Run the example ("c3run_example.job") in c3run

Change to the run script directory:
~~~
> cd $CTM3_ROOT/c3run
~~~
Open c3run_example.job in an editor of your choice (e.g., emacs, vi, gedit)
and get familiar with it. 
You do not have to change much to run the example.
In line 8 ("#SBATCH --account=$PROJECT"), replace $PROJECT with your project number. 
You can also execute:
~~~
> sed -i 's^$PROJECT^'$PROJECT'^g' c3run_example.job
~~~
Now change to your work directory:
~~~
> cd $WORK
~~~
Create a new directory - it should always have the same name as is given in line 5 of c3run_example.job - and change to it:
~~~
> mkdir C3RUN_example
> cd C3RUN_example
~~~
Finally send the model run to the batch system on ABEL/SAGA:
~~~
> sbatch $CTM3_ROOT/c3run/c3run_example.job
~~~
Wait for about 15-30 min.

## During the run

You can check the progress with:
~~~
> squeue
~~~
## After the run

Look at the results.
For a quick few on results, you can use ncview on abel
~~~
> module load ncview
~~~
or panoply (https://www.giss.nasa.gov/tools/panoply/download/).

~~~~~~~~~~~~~
Have fun!
~~~~~~~~~~~~~

---
# HOW TO 

## bash

As mentioned above, you do not want to execute all export commands by hand every time you log in to abel. We therefore create a bash script which will do the job.
Create a directory in your home directory on abel:
~~~
> mkdir ${HOME}/bin
~~~
Change into the directory:
~~~
> cd ${HOME}/bin
~~~
Create a new file:
~~~
> emacs setpaths &
~~~
Add the lines and don't forget to change the names in <>!
Save the file (EMACS: CTRL-x-s).
### 1. ABEL
~~~~~~~~~~~~~
#! /bin/bash
echo "To set work environment do set_up <opt>"
# Export the current project number
export PROJECT=<nnXXXXk>
# alias for squeue
alias squeue='squeue -lA ${PROJECT}'
# Export the work directory
export WORK=/work/users/$USER
# Export data storage on abel
export ASTRA=/projects/researchers/researchers01/sfalk
# Export CICERO directory
export CICERO=/work/projects/cicero/ctm_input

# Load modules and setup
module purge
echo "Settings for OsloCTM3git"
export CTM3_INPUT=${CICERO}/Indata_CTM3
export CTM3_USR_INPUT=$ASTRA/input_data/ctm_input/
export CTM3_ROOT=${HOME}/OsloCTM3/<username_project>
module load netcdf.intel/4.3.3.1
~~~~~~~~~~~~~

### 2. SAGA
~~~~~~~~~~~~~
#! /bin/bash
echo "To set work environment do set_up <opt>"
# Export the current project number
export PROJECT=<nnXXXXk>
# alias for squeue
alias squeue='squeue -lA ${PROJECT}'
# Export the work directory
export WORK=/cluster/work/users/${USER}
# Export CICERO directory
export CICERO=/cluster/projects/nn9188k/OsloCTM3

# Load modules and setup
module purge
echo "Settings for OsloCTM3git"
module load netCDF-Fortran/4.4.4-intel-2018b
export CTM3_INPUT=${CICERO}/Indata_CTM3
export CTM3_USR_INPUT=<your_forcing_files>
export CTM3_ROOT=${WORKSPACE}/OsloCTM3
export NETCDF_ROOT=$EBROOTNETCDFMINFORTRAN
~~~~~~~~~~~~~

## git

### 1. Fork

Sometimes it will be necessary to "sync" your fork to the original master.
Details about that process can be found in git-hub help 
(https://help.github.com/en/articles/syncing-a-fork).
Here we briefly summarize the steps:

Change into the working directory of your OsloCTM3:
~~~
> cd $CTM3_ROOT
~~~
Fetch changes from the "original" master:
~~~
> git fetch upstream
~~~
Switch to your "local" master:
~~~
> git checkout master
~~~
Merge the "original" master into your "local" master:
~~~
> git merge upstream/master
~~~
### 2. Branch

Sometimes it will be necessary to "merge" changes in the master branch beck into your own branch. Here we will give a brief summary of the steps:

Change into the working directory of your OsloCTM3
~~~
> cd $CTM3_ROOT
~~~
Change to "master" branch:
~~~
> git checkout master
~~~
Get the changes from remote:
~~~
> git pull
~~~
Change back into your own branch:
~~~
> git checkout <branch_name>
~~~
Merge the changes into your branch:
~~~
> git merge master
~~~

> Remark: 
> The automatic merging may cause "conflicts" in case your own changes clash 
> with the changes in master. Which means you have to go through the files 
> that cause this manually to "resolve" the "conflicts".

After merging is done you can push the updated branch to remote:
~~~
> git push
~~~

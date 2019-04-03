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

If you are not planning to change code, aka development, clone the master branch on git-hub:

> git clone git@github.com:NordicESMhub/OsloCTM3.git

Otherwise create your own development branch (fork) from master (on git-hub) and give it a descriptive (not too long) name. Clone that branch (fork).

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

> export INPUT_OCTM3=/work/projects/cicero/ctm_input/Indata_CTM3

Make sure you are member of the group "gf-ozone" and have access to

> ls /projects/researchers/researchers01/sfalk/input/ctm_input/

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


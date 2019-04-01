# Oslo CTM3
# - A global chemical transport model developed at MetOs

version: See MODEL_VERSION in cmn_size.F90
======================================================================
When refering to Oslo CTM3, use "Oslo CTM3 vX.Y".
======================================================================

User manual
The user manual is located in the file manual_osloctm3.pdf.
Or can be found on git-hub.

+++++++++++
Quick start
+++++++++++

+++
Before compiling
+++
If you are running on abel, you should do

> module purge

to unload all automatically load modules and do

> module load netcdf.intel/4.3.3.1

This will load all necessery dependencies.

Export the path of the working copy of your OsloCTM3:

> export CTM3_DIR=$HOME/<OsloCTM3>

and do the same with your work directory:

> export WORK=/work/users/<YOUR_USER_NAME>

Export your notur project number:

> export PROJECT=nnXXXX

Set an alias for the job queue on abel:

> alias squeue='squeue -lA ${PROJECT}'

Tip: 
Since these steps have to be repeated every time you log in,
it is wise to put the commants into a bash script (e.g. setpaths) 
that you can "source" after each login.

+++
Compile
+++
> cd $CTM3_DIR
> make -j

+++
Before running
+++
Make sure you are member of the group "cic-hpc"
and have the right to access

> ls /work/projects/cicero/ctm_input/

Then export the path to forcing etc.:

> export INPUT_OCTM3=/work/projects/cicero/ctm_input/Indata_CTM3

Make sure you are member of the group "gf-ozone"
and have access to

> /projects/researchers/researchers01/sfalk/input/ctm_input/

+++
To run the example ("c3run_example.job") in c3run
+++

> cd $WORK
> mkdir C3RUN_example
> cd C3RUN_example

Open c3run_example.job in an editor of your choice (e.g., emacs, vi, gedit)
and get familiar with it. 
You should not have to change anything to run the example.

Finally send the model to the batch system

> sbatch $CTM3_DIR/c3run/c3run_example.job

Wait for about 15-30 min.

+++
During the run
+++
You can check the progess with

> squeue

+++
After the run
+++

Look at the results.
You can use ncview on abel

> module load ncview

or panoply for a quick few on results.


Have fun!


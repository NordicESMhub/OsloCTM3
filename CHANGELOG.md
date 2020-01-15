*******************************
# Oslo CTM3 v1.1
~~~
Developer: 	Stefanie Falk
		Amund Søvde Haslerud
Organization: 	MetOs
Date:		2018-2019
~~~
*******************************
## Changes
	v1.0 -> v1.1
### Update of dry deposition scheme -> mOSaic
	- Mosaic approach accounting for different vegetation types
	  based on Simpson et al. (2012).
	- Published in the journal of Geoscience Model Development
	  Stefanie Falk and Amund Søvde Haslerud
	  Geosci. Model Dev., 12, 4705–4728, 2019
          https://doi.org/10.5194/gmd-12-4705-2019
	- For details of the implemented mOSaic scheme refer to the 
	  publication above.
	- The drydeposition velocities have changed.
### Code and structure has been tidied up
	- The routine for computing a local Monin-Obukhov length
	  has been moved from pbl_mixing.f90 to utilities_oslo.f90. 
	  The same routine is now used throughout the code.
	- The runscript (c3run) has been extended for easier use.
	- Modifications related to run on sigma2 SAGA.
	- Only the files in tables which are in use are copied to
	  WORK directory.
	- Hardcoded paths have been moved to .job.
	- The README has been extended regarding "quick start".
	- Easy switch between old and new dry deposition scheme.
	- Slightly more flexible output handling.


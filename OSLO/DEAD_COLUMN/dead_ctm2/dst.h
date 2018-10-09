! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dst.h,v 1.1 2003/04/15 14:41:35 alfgr Exp $ -*-f90-*-

! Purpose: dst.h sets and describes all tokens used in dust parameterization

! Usage: 
! #include <dst.h> /* Dust preprocessor tokens */ 

#ifndef DST_H
#define DST_H

#ifdef DST

! dst.h MUST have access to params.h to work
! CCM files all include #include params.h, but MATCH files do not
! The following three lines ensure access to params.h for MATCH files
#ifndef PARAMS_H
#include <params.h>
#endif /* not PARAMS_H */ 

! Host model:
! Tokens must be set to enable the correct output routines
! CCM activates the 5-argument CCM outfld() call
! When CCM is not defined, the MATCH outfld() call is used
! BXM activates the extensive 1-D netCDF output calls
! Normally CCM, and BXM are set in the Makefile rather than here
! BXM should imply CCM as well so CCM hooks get tested in BXM
!#define CCM
!#define BXM

! Set radiative mode:
! DST_RAD activates passing of dust to radiation routines
! If DST_RAD is not defined then the code should be bit-for-bit with CCM
! In this case the dust transport will occur, but optical depth is the only radiative property archived
#define DST_RAD
#ifdef DST_RAD
! When DST_RAD is defined, either DST_FDB or DST_FRC must also be defined
! DST_FDB implements radiative feedback mode
! DST_FDB calls the radiation routines only once
! The DST_FDB dynamics simulation will differ from CCM
! The other alternative when DST_RAD is defined is to use forcing mode, DST_FRC
! DST_FRC calls the radiation routines with and without dust
! The difference is the diagnostic radiative forcing due to dust
! Archiving all the DST_FRC fields involves a lot of bookkeeping
! The DST_FRC dynamics simulation should be bit-for-bit with CCM
#define DST_FRC
!#define DST_FDB
#endif /* not DST_RAD */ 

! Set convective transport mode:
! DST_TRN_CNV calls conv_tran after wet deposition rather than in deep convection
! DST_TRN_CNV applies to CCM 3.X runs only, since MATCH does the former by default
! DST_TRN_CNV improves realism of transport processes so should generally be defined
#define DST_TRN_CNV

! Set land surface mode:
! DST_LSM activates passing of LSM data into CCM
! Currently, DST_LSM is required by the mobilization scheme
#define DST_LSM

! Set chemistry mode:
! DST_CHM activates input and processing of chemistry arrays
#define DST_CHM

! Set debugging modes:
! DST_DBG activates non-vectorized sanity checks and verbose output:
!c++alfgr
! DST_DBG is set in Oslo CTM2 Makefile if wanted
!c--alfgr
! 1. Tendency procedures (e.g., dst_dps_wet()) will check for unreasonable mixing ratios after each adjustment
! 2. Surface soil properties at the debug gridpoint are sometimes printed
! 3. qneg3() is called after each source/sink routine
! Production runs are either on CRAY or SGI nowadays so default is no DST_DBG on these machines
!#if (!defined CRAY) && (!defined SGI)
!#define DST_DBG
!#endif /* not CRAY or SGI */
! DST_MSS_BDG activates global mass budget output into netCDF file
! The diagnostic mass budget requires running global averages on many quantities each timestep
! This can be very expensive, but makes checking for mass balance much easier
!#if (!defined CRAY) && (!defined SGI)
!#define DST_MSS_BDG
!#endif /* not CRAY or SGI */

! Automated check for erroneous token combinations.
! No editing should need to be done beneath this line.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dynamic source scheme requires interactive LSM data
#if (defined DST_DYN) && (!defined DST_LSM)
      "ERROR dst.h defines DST_DYN but not DST_LSM"
#endif /* not (DST_DYN && DST_TMS) */ 

#ifdef DST_RAD
! Exactly one of DST_FDB and DST_FRC must be defined when DST_RAD is defined 
#if (defined DST_FDB) && (defined DST_FRC)
      "ERROR dst.h defines both DST_FDB and DST_FRC"
#endif /* not (DST_FDB && DST_FRC) */ 
#if (!defined DST_FDB) && (!defined DST_FRC)
      "ERROR dst.h defines neither DST_FDB nor DST_FRC"
#endif /* (DST_FDB && DST_FRC) */ 
#endif /* not DST_RAD */ 

#endif /* not DST */

#endif /* not DST_H */


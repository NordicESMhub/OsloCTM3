! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/params.h,v 1.3 2007/04/11 08:35:28 alfgr Exp $ -*-f90-*- 

! Purpose: params.h sets all tokens used in dust parameterization

! Usage: 
! #include <params.h> /* Preprocessor tokens */ 

! DST_NBR is number of dust tracers
! DST_IDX_SRT is array index of first dust component in constituent array
! Provide defaults for DST_IDX_SRT and DST_NBR even when DST is not defined
! This permits files in dst to compile for non-dust control runs
! CCM and MATCH differ significantly in their enumeration of advected species
! CCM guarantees that water vapor is species 1, thus DST_IDX_SRT may be 2
! MATCH (with prognostic cloudwater) sets condensed water to species PCNST 
! MATCH carries water vapor in a separate array (shadv), thus DST_IDX_SRT is 1

#ifndef PARAMS_H
#define PARAMS_H

! Define Box Model parameters if appropriate, else default to global models
#ifdef BXM

! Begin Box Model section
#define PLAT 1 
#define PLON 1 
#define PLEV 1
#define DST_NBR 4
#define DST_IDX_SRT 2
#define PCNST 1+DST_NBR

! End Box Model section
#else /* not BXM */ 

! Define CCM parameters if appropriate, else default to MATCH/MOZART
#ifdef CCM 

! Begin CCM section
#ifdef DST
#define DST_NBR 4 
#define DST_IDX_SRT 2
#define PCNST 1+DST_NBR
#else /* not DST */
#define DST_IDX_SRT 1
#define DST_NBR 1 
#define PCNST 1
#endif /* not DST */

#define PNATS 0
#define PLEV 18
#define PLEVR 18
#define POZLEV 23

! End CCM section
#else /* not CCM */ 

! Begin MATCH/MOZART section
#ifdef DST
#ifdef OC_M7
   !++alfgr M7 not standard: 2
#define DST_NBR 2 
#else
   ! standard: 8
#define DST_NBR 8
#endif /* OC_M7 */
#define DST_IDX_SRT 1
#define PCNST DST_NBR
#else /* not DST */
#define DST_IDX_SRT 1
#define DST_NBR 1 
#define PCNST 2
#endif /* not DST */

#define PNATS  0

!++alfgr
#ifdef L40
#ifdef COLLAPSE
#define PLEV 37
#else
#define PLEV 40
#endif /* COLLAPSE */
#endif

#ifdef L60
#ifdef COLLAPSE
#define PLEV 57
#else
#define PLEV 60
#endif /* COLLAPSE */
#endif
!--alfgr

#define CALC_ETADOT
#define DI_VDIFF
#define DI_CONV_CCM
#define DI_CLOUD_PHYS

! End MATCH/MOZART section
#endif /* not CCM */ 

! Define Resolution parameters for global models
! NB: MATCH/MOZART do not need PTRM, PTRN, or PTRK
#ifdef T5
#define PLON 16
#define PLAT 8
#define PTRM 5
#define PTRN 5
#define PTRK 5
#endif /* not T5 */

#ifdef HT21
#define PLON 64
#define PLAT 32
#define PTRM 21
#define PTRN 21
#define PTRK 21
#endif /* not T21 */

#ifdef T31
#define PLON 96
#define PLAT 48
#define PTRM 31
#define PTRN 31
#define PTRK 31
#endif /* not T31 */

#ifdef HT42
#define PLON 128
#define PLAT 64
#define PTRM 42
#define PTRN 42
#define PTRK 42
#endif /* not T42 */

#ifdef T62
#define PLON 192
#define PLAT 94
#define PTRM 62
#define PTRN 62
#define PTRK 62
#endif /* not T62 */

#ifdef T63
#define PLON 192
#define PLAT 96
#define PTRM 63
#define PTRN 63
#define PTRK 63
#endif /* not T63 */

#ifdef H1x1
#define PLON 360
#define PLAT 180
#define PTRM 319
#define PTRN 319
#define PTRK 319
#endif /* not T319 */

! End Global Model section
#endif /* not BXM */ 

#endif /* not PARAMS_H */ 


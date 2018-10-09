! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstcmnini.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Routines to initialize modules variables throughout the dust model, 
! currently including these modules
! dstchm.F90, dstdbg.F90, dstodx.F90, dstscv.F90
! This module simply initializes all those modules whose initialization 
! is host model independent, and thus one call (dst_msc_cmn_ini()) suffices
! to initialize these modules.
! The order of initialization of other modules is host-model dependent
! and the calls to initialize them must therefore be made explicitly in
! the host model

! Usage:
! use dstcmnini ! [mdl] Module initialization

! dst.h needed for DST_CHM token
!#include <dst.h> /* Dust preprocessor tokens */

module dstcmnini ! [mdl] Module initialization
  implicit none

contains
  
subroutine dst_msc_cmn_ini()
  ! Purpose: Initialize parameter-dependent fields in miscellaneous common blocks
  ! dst_msc_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  use dstblm,only:dst_blm_cmn_ini ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
  use dstdbg,only:dst_dbg_cmn_ini ! [mdl] Debugging information for dust model
  use dstodx,only:dst_odx_cmn_ini ! [mdl] Optical depth information
  use dstchm,only:dst_chm_cmn_ini ! [mdl] Chemical properties of dust
  use dstscv,only:dst_scv_cmn_ini ! [mdl] Aerosol scavenging properties
  use dstsltsbl,only:dst_slt_sbl_cmn_ini ! [mdl] Saltation sandblasting physics
  ! Following three modules used to be block data which had to be initialized in main()
  ! fxm: I assume this is no longer necessary
  !  use dstlsm ! Initialize LSM block data
  !  use dstpsd ! Initialize aerosol block data
  !  use dstrad ! Initialize radiation block data
  implicit none
  ! Input
  ! Output
  ! Input/Output
  ! Local
  ! Main code
  
  ! Initialize boundary layer meteorology (checked ok for CTM3)
  call dst_blm_cmn_ini()
  
  ! Initialize debugging information (checked ok for CTM3)
  call dst_dbg_cmn_ini()
  
  ! Initialize optical depth information (checked ok for CTM3)
  call dst_odx_cmn_ini()
  
  ! Initialize wet scavenging information (checked ok for CTM3)
  call dst_scv_cmn_ini()

  ! Initialize saltation sandblasting physics (checked ok for CTM3)
  call dst_slt_sbl_cmn_ini()
  
#ifdef DST_CHM
  ! Initialize dust chemistry (checked ok for CTM3)
  call dst_chm_cmn_ini()
#endif /* not DST_CHM */
  
  return
end subroutine dst_msc_cmn_ini                       ! end dst_msc_cmn_ini()

end module dstcmnini ! [mdl] Module initialization

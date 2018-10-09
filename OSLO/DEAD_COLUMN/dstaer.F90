! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstaer.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: dstaer.F90 contains common blocks which store important 
! microphysical properties of the mineral dust aerosol. 
! Variables in dstaer are used each timestep by source and sink routines
! Related, but strictly diagnostic, microphysical quantities are stored in dstpsd.F90

! Most of these variables are initialized in dst_psd_ini()
! dst_psd_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
! dstaer.F90 MUST have access to dstgrd.F90 and to pmgrid.F90 to work

! Usage:
! use dstaer ! [mdl] Aerosol microphysical properties

module dstaer             ! [mdl] Aerosol microphysical properties
  use dead_precision             ! [mdl] Precision r8, i8, ...
  use dstgrd,only:dst_nbr,dst_src_nbr ! [mdl] Dust grid sizes
  use pmgrid,only:plond,plon ! [mdl] Spatial resolution parameters 
  implicit none
  save                      ! [stt] Changes to common variables are sticky
  
  ! dst_aer_cmn is computed at run time, initialized in dst_psd_ini()
  real(r8) dmt_vwr(dst_nbr) ! [m] Mass weighted diameter resolved
  real(r8) dns_aer(dst_nbr) ! [kg m-3] Particle density
  real(r8) mss_frc_src(dst_src_nbr) ! [frc] Mass fraction of source distribution
  real(r8) mss_frc_trn_dst_src(dst_nbr) ! [frc] Fraction of transported dust mass at source
  real(r8) ovr_src_snk_mss(dst_src_nbr,dst_nbr) ! [frc] Overlap of src with snk
  real(r8) ovr_src_snk_mss_ttl ! [frc] Total transported mass fraction of dust flux
  real(r8) stk_crc(dst_nbr) ! [frc] Correction to Stokes settling velocity
  real(r8) sfc_spc_rsl(dst_nbr) ! [m2 kg-1] Specific surface area resolved 
  
end module dstaer

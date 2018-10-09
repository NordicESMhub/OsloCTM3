! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstdbg.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Debugging module for dust model
! Initialized in dst_dbg_cmn_ini() which is called by CCM:inti() 

! Usage: 
! use dstdbg ! [mdl] Debugging information for dust model

module dstdbg             ! [mdl Debugging information for dust model
  use precision             ! [mdl] Precision r8, i8, ...
  implicit none
  save                      ! [stt] Changes to common variables are sticky
  
  ! Debug info initialized in dst_dbg_cmn_ini()
  integer lat_dbg           ! [idx] Latitude  index for diagnostics
  integer lev_dbg           ! [idx] Level     index for diagnostics
  integer lon_dbg           ! [idx] Longitude index for diagnostics
  integer time_dbg          ! [idx] Time      index for diagnostics
  
  ! Physics info initialized in dst_dbg_cmn_ini()
  real(r8) flx_mss_sgn      ! [kg m-2 s-1] Minimum significant flux of dust
  real(r8) flx_mss_mxm      ! [kg m-2 s-1] Maximum plausible flux of dust
  real(r8) q_dst_mxm        ! [kg kg-1] Maximum plausible mass mixing ratio of dust
  real(r8) q_dst_sgn        ! [kg kg-1] Minimum significant mass mixing ratio of dust
contains
  
  subroutine dst_dbg_cmn_ini()
    ! Purpose: Initialize parameter-dependent fields debugging common block
    ! dst_dbg_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
    use precision ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    ! Output
    ! Input/Output
    ! Local
    integer idx               ! [idx] Counting index 
    ! Main code
    
    ! All source/sink routines compare new dust mixing ratio to q_dst_mxm
    ! q_dst > q_dst_mxm indicates unbelievably high dust emissions
    ! Model should abort until veracity can be established
    ! q_dst_mxm should exceed greatest recorded atmospheric dust concentration
    q_dst_mxm=1.0e-3          ! [kg kg-1] Maximum plausible mass mixing ratio of dust
    if(q_dst_mxm < 1.0e-7) stop 'dst_dbg_cmn_ini() reports q_dst_mxm too small, will trip for reasonable dust concentrations'
    
    ! Source/sink routines compare new dust mixing ratio to q_dst_sgn
    ! q_dst < q_dst_sgn indicates unbelievably low dust emissions
    ! probably due to numerical artifacts like rounding, underflow, order of operations
    ! q_dst_sgn = 0.0 is fine, but numeric artifacts will not be eliminated
    ! This can result in underflows from operations involving very small values
    ! q_dst_sgn should not exceed about 1.0e-20 kg kg-1
    ! WARNING: Model is free to (and does) eliminate concentrations lower than q_dst_sgn
    q_dst_sgn=1.0e-30_r8         ! [kg kg-1] Minimum plausible mass mixing ratio of dust
    if(q_dst_sgn > 1.0e-20_r8) stop 'dst_dbg_cmn_ini() reports q_dst_sgn too large, fixers may cause significant mass imbalances'
    
    ! Threshold flux indicating presence of "significant" dust
    ! Evaporative source routine prints warning when evaporative dust flux 
    ! exceeds both the available precipitating flux and flx_mss_sgn
    flx_mss_sgn=1.0e-30_r8       ! [kg m-2 s-1] Minimum significant flux of dust
    
    ! Source routine compares new dust flux to flx_mss_mxm
    ! flx_mss_vrt_dst_ttl > flx_mss_mxm indicates suspiciously high dust fluxes
    flx_mss_mxm=100.0e-9      ! [kg m-2 s-1] Maximum plausible flux of dust
    if(flx_mss_mxm < 1.0e-9) stop 'dst_dbg_cmn_ini() reports flx_mss_mxm too small, will trip for reasonable mass fluxes'
    
    ! Initialize common block used for debugging
#ifdef BXM
    lat_dbg=1                 ! [idx] Latitude  index for diagnostics
    lev_dbg=1                 ! [idx] Level     index for diagnostics
    lon_dbg=1                 ! [idx] Longitude index for diagnostics
    time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not BXM */
#ifdef T5
    ! T5 lon = i = 13 = 112.5W and lat = j = 3 = 31.7S (land): Chile
    lat_dbg=3                 ! [idx] Latitude  index for diagnostics
    lev_dbg=11                ! [idx] Level     index for diagnostics
    lon_dbg=13                ! [idx] Longitude index for diagnostics
    time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T5 */
#ifdef T42
    ! T42 lon = i = 99 = 84.3W and lat = j = 37 = 12.5N (land): Panama isthmus
    ! T42 lon = i = 98 = 87.1W and lat = j = 37 = 12.5N (ocean): Pacific
    ! T42 lon = i = 115 = 39.8W and lat = j = 61 = 79.5N (ocean): Greenland
    ! T42 lon = i = 125 = 11.25W and lat = j = 40 = 20.9N (land): Desert at center of Mauritania
    ! T42 lon = i = 106 = 64.6W and lat = j = 22 = -29.3S (land): near Atacama
    ! fxm: pass in coordinate arrays to simplify index identification
    ! lat_dbg=vec_val2idx(lat,lat_nbr,20.9) ! [idx] Latitude  index for diagnostics
    ! lon_dbg=vec_val2idx(lon,lon_nbr,348.75) ! [idx] Longitude index for diagnostics
  lat_dbg=40                ! [idx] Latitude  index for diagnostics
  lev_dbg=9                 ! [idx] Level     index for diagnostics
  lon_dbg=125               ! [idx] Longitude index for diagnostics
  time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T42 */ 
#ifdef T62
  ! T62 lon = i = 159 = 63W and lat = j = 29 = 35.0S (land): Pasture in Argentina
  ! T62 lon = i = 187 = 11.25W and lat = j = 58 = 20.0N (land): Desert at center of Mauritania
  ! T62 lon = i = 57 = 105E and lat = j = 66 = 35.23N (land): Desert near Takla Makan
  lat_dbg=29                ! [idx] Latitude  index for diagnostics
  lev_dbg=9                 ! [idx] Level     index for diagnostics
  lon_dbg=159               ! [idx] Longitude index for diagnostics
  time_dbg=0                ! [idx] Time      index for diagnostics
#endif /* not T62 */ 
  
  return
end subroutine dst_dbg_cmn_ini                       ! end dst_dbg_cmn_ini()

end module dstdbg

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstsfc.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Surface fields needed by dust parameterization
! dstsfc.F90 needs to work independently of host model (MATCH, CCM, CSM, MOZART)
! dstsfc.F90 contains interface routines to access 2-D surface variables stored in common blocks

! Usage: 
! use dstsfc ! [mdl] 2-D surface fields on PLON x PLAT grid

! params.h needed for PLON and PLAT
!//CTM3: These are inclued in pmgrid.F90
!#include <params_dust.h>

!// ------------------------------------------------------------------
module dstsfc ! [mdl] 2-D surface fields on PLON x PLAT grid
  !// ------------------------------------------------------------------
  use dead_precision ! [mdl] Precision r8, i8, ...
  use pmgrid ! [mdl] Spatial resolution parameters
  use dstgrd,only:bln_nbr ! [mdl] Dust grid sizes
  implicit none
  save ! [stt] Changes to common variables are sticky
  private ! [stt] Symbols are private unless individually qualified as public
!  public::dst_sfc_set ! [sbr] Box model method to set surface tvbds common block
  public::flx_rad_sfc_get ! [sbr] Return latitude slice of surface radiative fluxes
  public::flx_rad_sfc_set ! [sbr] Set latitude slice of surface radiative fluxes
  public::sfc_typ_get ! [sbr] Return latitude slice of LSM surface type
  public::soi_txt_get ! [sbr] Return latitude slice of soil texture
  public::tpt_gnd_soi_get ! [sbr] Return latitude slice of soil temperature and ground temperature
  public::tpt_gnd_soi_set ! [sbr] Set latitude slice of surface and soil temperatures
  public::vwc_sfc_get ! [sbr] Return latitude slice of volumetric water content
  public::vwc_sfc_set ! [sbr] Set latitude slice of surface volumetric water content
  
  ! Time-independent surface information supplied by boundary dataset
  ! Initialized in CCM:dst/dst_tibds_ini()
  integer,public::sfc_typ(PLON,PLAT) ! [idx] LSM surface type (0..28)
  real(r8),public::lnd_frc_dry(PLON,PLAT) ! [frc] Land surface that is not lake or wetland (by area)
  real(r8),public::mbl_bsn_fct(PLON,PLAT) ! [frc] Erodibility factor
  real(r8),public::mss_frc_CaCO3(PLON,PLAT) ! [frc] Mass fraction of CaCO3
  real(r8),public::mss_frc_cly(PLON,PLAT) ! [frc] Mass fraction of clay
  real(r8),public::mss_frc_snd(PLON,PLAT) ! [frc] Mass fraction of sand
  real(r8),public::sfc_frc_bln(PLON,PLAT,bln_nbr) ! [frc] Fraction of 4 available soil types
  ! Time-varying surface information supplied by CCM
  ! Variables initialized in CCM:physics/tphysbc()
  ! fxm: Probably should store radiative fluxes elsewhere
  real(r8) flx_LW_dwn_sfc(PLON,PLAT) ! [W m-2] Longwave downwelling flux at surface
  real(r8) flx_SW_abs_sfc(PLON,PLAT) ! [W m-2] Solar flux absorbed by ground
  
  ! NB: Use PLON,PLAT tokens instead of plon,plat
  ! params.h was dropped from LSM after CCM3.5, but we re-include it for PLON,PLAT tokens
  ! In LSM, plon and plat are passed as variables into lsmdrv()
  ! In CCM, plon and plat are parameters in every routine with pmgrid.h or prgrid.h
  ! dst_sfc needs to work in both LSM and CCM so use PLON,PLAT
  
  ! Time-varying surface information supplied by LSM
  ! Variables initialized in CCM:lsm/lsmdrv()
  logical sfc_set_sgs_mbl(PLON,PLAT) ! [flg] Are surface properties already set?
  real(r8) tpt_gnd(PLON,PLAT)   ! [K] Ground temperature
  real(r8) tpt_soi(PLON,PLAT)   ! [K] Soil temperature
  real(r8) vwc_sfc(PLON,PLAT)   ! [m3 m-3] Volumetric water content

  !// ------------------------------------------------------------------

  
contains
  

  !// ------------------------------------------------------------------
  subroutine flx_rad_sfc_set( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       flx_LW_dwn_sfc_in,   & ! I [W m-2] Longwave downwelling flux at surface
       flx_SW_abs_sfc_in)   ! I [W m-2] Solar flux absorbed by ground
    !// ------------------------------------------------------------------
    ! Set latitude slice of surface radiative fluxes
    ! flx_rad_sfc_set() is called by CCM:physics/tphysbc(), MATCH:physlic()
    !//
    !// Rewritten for Oslo CTM3, uses only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in) :: lat_idx, & ! I [idx] Latitude index
                          lon_idx    ! I [idx] Longitude index
    real(r8),intent(in)::flx_LW_dwn_sfc_in ! I [W m-2] Longwave downwelling flux at surface
    real(r8),intent(in)::flx_SW_abs_sfc_in ! I [W m-2] Solar flux absorbed by ground
    !// ------------------------------------------------------------------

    ! Main Code
    flx_LW_dwn_sfc(lon_idx,lat_idx) = flx_LW_dwn_sfc_in
    flx_SW_abs_sfc(lon_idx,lat_idx) = flx_SW_abs_sfc_in

    !// ------------------------------------------------------------------
  end subroutine flx_rad_sfc_set                       ! end flx_rad_sfc_set()
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine flx_rad_sfc_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       flx_LW_dwn_sfc_out,  & ! O [W m-2] Longwave downwelling flux at surface
       flx_SW_abs_sfc_out)  ! O [W m-2] Solar flux absorbed by ground
    !// ------------------------------------------------------------------
    ! Return latitude slice of surface radiative fluxes
    ! flx_rad_sfc_get() is called by dst_mbl()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in) :: lat_idx, & ! I [idx] Latitude index
                          lon_idx    ! I [idx] Longitude index
    ! Output
    real(r8),intent(out)::flx_LW_dwn_sfc_out ! O [W m-2] Longwave downwelling flux at surface
    real(r8),intent(out)::flx_SW_abs_sfc_out ! O [W m-2] Solar flux absorbed by ground
    !// ------------------------------------------------------------------

    ! Main code
    flx_LW_dwn_sfc_out = flx_LW_dwn_sfc(lon_idx,lat_idx)
    flx_SW_abs_sfc_out = flx_SW_abs_sfc(lon_idx,lat_idx)

    !// ------------------------------------------------------------------
  end subroutine flx_rad_sfc_get                       ! end flx_rad_sfc_get()
  !// ------------------------------------------------------------------




  !// ------------------------------------------------------------------
  subroutine soi_txt_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       lnd_frc_dry_out,     & ! O [frc] Dry land fraction
       mbl_bsn_fct_out,     & ! O [frc] Erodibility factor
       mss_frc_CaCO3_out,   & ! O [frc] Mass fraction of CaCO3
       mss_frc_cly_out,     & ! O [frc] Mass fraction of clay
       mss_frc_snd_out,     & ! O [frc] Mass fraction of sand
       sfc_frc_bln_out)       ! O [frc] Fraction of 4 available soil types
    !// ------------------------------------------------------------------
    ! Return latitude slice of soil texture
    ! soi_txt_get() is called by dst_mbl()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    real(r8),parameter :: mss_frc_cly_glb=0.20 ! [kg kg-1] Ad hoc globally uniform clay mass fraction

    ! Input
    integer,intent(in) :: lat_idx, & ! I [idx] Latitude index
                          lon_idx    ! I [idx] Longitude index

    ! Output
    real(r8),intent(out)::lnd_frc_dry_out ! O [frc] Dry land fraction
    real(r8),intent(out)::mbl_bsn_fct_out ! O [frc] Erodibility factor
    real(r8),intent(out)::mss_frc_CaCO3_out ! O [frc] Mass fraction of CaCO3
    real(r8),intent(out)::mss_frc_cly_out ! O [frc] Mass fraction of clay
    real(r8),intent(out)::mss_frc_snd_out ! O [frc] Mass fraction of sand
    real(r8),intent(out)::sfc_frc_bln_out(bln_nbr) ! [frc] Fraction of 4 available soil types
    !// ------------------------------------------------------------------

    ! Main code
    lnd_frc_dry_out = lnd_frc_dry(lon_idx,lat_idx)
    mbl_bsn_fct_out = mbl_bsn_fct(lon_idx,lat_idx)
    ! fxm: CaCO3 currently has missing value of 1.0e36 which causes problems
    if (mss_frc_CaCO3(lon_idx,lat_idx) <= 1.0_r8) then
       mss_frc_CaCO3_out = mss_frc_CaCO3(lon_idx,lat_idx)
    else
       mss_frc_CaCO3_out = 0.0_r8
    endif
    ! fxm: Temporarily set mss_frc_cly used in mobilization to globally uniform SGS value of 0.20
    ! Put excess mass fraction into sand
    ! mss_frc_cly_out(lon_idx)=mss_frc_cly(lon_idx,lat_idx)
    ! mss_frc_snd_out(lon_idx)=mss_frc_snd(lon_idx,lat_idx)+mss_frc_cly_dlt
    mss_frc_cly_out = mss_frc_cly_glb
    mss_frc_snd_out = mss_frc_snd(lon_idx,lat_idx)+mss_frc_cly(lon_idx,lat_idx)-mss_frc_cly_glb


    sfc_frc_bln_out(:)=sfc_frc_bln(lon_idx,lat_idx,:)


    !// ------------------------------------------------------------------
  end subroutine soi_txt_get                       ! end soi_txt_get()
  !// ------------------------------------------------------------------
  



  !// ------------------------------------------------------------------
  subroutine sfc_typ_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       sfc_typ_out)           ! O [idx] LSM surface type (0..28)
    !// ------------------------------------------------------------------
    ! Return latitude slice of LSM surface type
    ! sfc_typ_get() is called by dst_mbl(), dst_dps_dry()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in) :: lat_idx, & ! I [idx] Latitude index
                          lon_idx    ! I [idx] Longitude index
    ! Output
    integer,intent(out) :: sfc_typ_out ! O [idx] LSM surface type (0..28)
    !// ------------------------------------------------------------------
    ! Main code
    sfc_typ_out = sfc_typ(lon_idx,lat_idx)

    !// ------------------------------------------------------------------
  end subroutine sfc_typ_get                       ! end sfc_typ_get()
  !// ------------------------------------------------------------------
  


  !// ------------------------------------------------------------------
  subroutine tpt_gnd_soi_set( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       tpt_gnd_in,          & ! I [K] Ground temperature
       tpt_soi_in)            ! I [K] Soil temperature
    !// ------------------------------------------------------------------
    ! Set latitude slice of surface and soil temperatures
    ! tpt_gnd_soi_set() is called by MATCH:physlic()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in)::lat_idx, & ! I [idx] Latitude index
                        lon_idx    ! I [idx] Longitude index
    real(r8),intent(in)::tpt_gnd_in ! I [K] Ground temperature
    real(r8),intent(in)::tpt_soi_in ! I [K] Soil temperature
    !// ------------------------------------------------------------------

    ! Main Code
    tpt_gnd(lon_idx,lat_idx) = tpt_gnd_in
    tpt_soi(lon_idx,lat_idx) = tpt_soi_in

    !// ------------------------------------------------------------------
  end subroutine tpt_gnd_soi_set                       ! end tpt_gnd_soi_set()
  !// ------------------------------------------------------------------


  !// ------------------------------------------------------------------
  subroutine tpt_gnd_soi_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       tpt_gnd_out,         & ! O [K] Ground temperature
       tpt_soi_out)         ! O [K] Soil temperature
    !// ------------------------------------------------------------------
    ! Return latitude slice of soil temperature and ground temperature
    ! tpt_gnd_soi_get() is called by dst_mbl()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in)::lat_idx, & ! I [idx] Latitude index
                        lon_idx    ! I [idx] Longitude index
    ! Output
    real(r8),intent(out)::tpt_gnd_out ! O [K] Ground temperature
    real(r8),intent(out)::tpt_soi_out ! O [K] Soil temperature
    !// ------------------------------------------------------------------

    tpt_gnd_out = tpt_gnd(lon_idx,lat_idx)
    tpt_soi_out = tpt_soi(lon_idx,lat_idx)

    !// ------------------------------------------------------------------
  end subroutine tpt_gnd_soi_get                       ! end tpt_gnd_soi_get()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine vwc_sfc_set( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       vwc_sfc_in)          ! I [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------
    ! Set latitude slice of surface volumetric water content
    ! vwc_sfc_set() is called by MATCH:physlic()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in)::lat_idx, & ! I [idx] Latitude index
                        lon_idx    ! I [idx] Longitude index
    real(r8),intent(in)::vwc_sfc_in ! I [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------

    ! Main Code
    vwc_sfc(lon_idx,lat_idx) = vwc_sfc_in ! [m3 m-3]

    !// ------------------------------------------------------------------
  end subroutine vwc_sfc_set                       ! end vwc_sfc_set()
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine vwc_sfc_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
       vwc_sfc_out)         ! O [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------
    ! Return latitude slice of volumetric water content
    ! vwc_sfc_get() is called by dst_mbl()
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    integer,intent(in)::lat_idx, & ! I [idx] Latitude index
                        lon_idx    ! I [idx] Longitude index
    ! Output
    real(r8),intent(out)::vwc_sfc_out ! O [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------

    ! Main code
    vwc_sfc_out = vwc_sfc(lon_idx,lat_idx)

    !// ------------------------------------------------------------------
  end subroutine vwc_sfc_get                       ! end vwc_sfc_get()
  !// ------------------------------------------------------------------


  !// ------------------------------------------------------------------
end module dstsfc
!// ------------------------------------------------------------------

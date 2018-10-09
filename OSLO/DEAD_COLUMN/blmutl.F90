! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/blmutl.F90,v 1.3 2003/08/29 11:46:53 alfgr Exp $

! Purpose: blmutl.F90 contains driver routines for grid-scale boundary layer meteorology
! Usage:
! use blmutl ! [mdl] Boundary layer meteorology driver

module blmutl ! [mdl] Boundary layer meteorology driver
  use dead_precision ! [mdl] Precision r8, i8, ...
  !  use dstmblutl ! [mdl] Mobilization utilities
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  private::blm_ice ! [sbr] Compute boundary layer properties over sea ice
  private::blm_lnd ! [sbr] Compute boundary layer properties over land
  private::blm_ocn ! [sbr] Compute boundary layer properties over ocean
  private::mno_stb_crc_heat_uns_get ! [fnc] Stability correction for heat, unstable case
  private::mno_stb_crc_mmn_uns_get ! [fnc] Stability correction for momentum, unstable case
  private::rgh_mmn_get ! [sbr] Set roughness length
!  private::rgh_zpd_get ! [sbr] Set roughness length and zero plane displacement
  private::zpd_get ! [sbr] Set zero plane displacement
  public::blm_glb ! [sbr] Solve boundary layer meteorology on global scale
  public::oro_is_ice ! [fnc] Point is > 50% sea ice
  public::oro_is_lnd ! [fnc] Point is > 50% land
  public::oro_is_ocn ! [fnc] Point is > 50% ocean
  public::snw_frc_get ! [fnc] Convert snow depth to fractional snow cover
  public::wnd_rfr_get ! [fnc] Interpolate wind speed to reference height
  public::xch_cff_mmn_ocn_ntr_get ! [frc] Neutral 10 m drag coefficient over ocean

contains
  
  !// ------------------------------------------------------------------
  logical function oro_is_ocn(oro_val)
    ! Purpose: True if > 50% ocean
    implicit none
    real(r8),intent(in) :: oro_val ! [frc] Orography
    oro_is_ocn = nint(oro_val)==0
  end function oro_is_ocn
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  logical function oro_is_lnd(oro_val)
    ! Purpose: True if > 50% land
    implicit none
    real(r8),intent(in) :: oro_val ! [frc] Orography
    oro_is_lnd = nint(oro_val)==1
  end function oro_is_lnd
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  logical function oro_is_ice(oro_val)
    ! Purpose: True if > 50% sea ice
    implicit none
    real(r8),intent(in) :: oro_val ! [frc] Orography
    oro_is_ice = nint(oro_val)==2
  end function oro_is_ice
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function mno_stb_crc_heat_uns_get( & ! [fnc] Stability correction for heat, unstable case
       sml_fnc_mmn_uns_rcp) ! I [frc] Reciprocal of Monin-Obukhov similarity function
    ! Purpose: Stability correction for heat, unstable case
    ! Given the reciprocal of the Monin-Obukhov similarity function 
    ! (usually called phi) for momentum in an unstable atmosphere, return the 
    ! stability correction factor for heat, usually called psi
    ! References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
    ! Currently this function is BFB with CCM:dom/flxoce()
    implicit none
    real(r8),intent(in) :: sml_fnc_mmn_uns_rcp ! I [frc] Reciprocal of Monin-Obukhov similarity function
    mno_stb_crc_heat_uns_get = & ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
         2.0_r8*log((1.0_r8+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0_r8)
  end function mno_stb_crc_heat_uns_get
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function mno_stb_crc_mmn_uns_get( & ! [fnc] Stability correction for momentum, unstable case
       sml_fnc_mmn_uns_rcp) ! I [frc] Reciprocal of Monin-Obukhov similarity function
    ! Purpose: Stability correction for heat, unstable case
    ! Given the reciprocal of the Monin-Obukhov similarity function 
    ! (usually called phi) for momentum in an unstable atmosphere, return the 
    ! stability correction factor for momentum, usually called psi 
    ! References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
    ! Currently this function is BFB with CCM:dom/flxoce()
    implicit none
    real(r8),intent(in) :: sml_fnc_mmn_uns_rcp ! I [frc] Reciprocal of Monin-Obukhov similarity function
    mno_stb_crc_mmn_uns_get = & ! [frc] Ary88 p. 167, Bru82 p. 71, SeP97 p. 869, Bon96 p. 52, BKL97 p. F1, LaP81 p. 325
         log((1.0_r8+sml_fnc_mmn_uns_rcp*(2.0_r8+sml_fnc_mmn_uns_rcp))*(1.0_r8+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/8.0_r8) &
         -2.0_r8*atan(sml_fnc_mmn_uns_rcp)+1.571_r8
  end function mno_stb_crc_mmn_uns_get
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function xch_cff_mmn_ocn_ntr_get(wnd_10m_ntr)
    ! Purpose: [frc] Neutral 10 m drag coefficient over ocean
    implicit none
    real(r8),intent(in) :: wnd_10m_ntr ! [m s-1] Wind speed at 10 m
    xch_cff_mmn_ocn_ntr_get = 0.0027_r8/wnd_10m_ntr + 0.000142_r8 + 0.0000764_r8*wnd_10m_ntr ! [frc] LaP82 CCM:dom/flxoce(), NOS97 p. I2
  end function xch_cff_mmn_ocn_ntr_get
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine rgh_mmn_get( & ! [sbr] Set roughness length
       oro,                 & ! I [frc] Orography
       rgh_mmn,             & ! O [m] Roughness length momentum
       sfc_typ,             & ! I [idx] LSM surface type (0..28)
       snw_frc,             & ! I [frc] Fraction of surface covered by snow
       wnd_10m)             ! I [m s-1] 10 m wind speed
    !// ------------------------------------------------------------------
    ! Purpose: Set roughness length
    ! NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
    ! NB: Currently all lakes are treated as unfrozen
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    use dstlsm ! [mdl] LSM data
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    real(r8),parameter :: rgh_mmn_ice_lak=0.04 ! [m] Roughness length over frozen lakes Bon96 p. 59
    real(r8),parameter :: rgh_mmn_ice_lnd=0.05 ! [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
    real(r8),parameter :: rgh_mmn_ice_ocn=0.0005 ! [m] Roughness length over sea ice BKL97 p. F-3
    real(r8),parameter :: rgh_mmn_lak_wrm=0.001 ! [m] Roughness length over unfrozen lakes Bon96 p. 59
    real(r8),parameter :: rgh_mmn_snw=0.04 ! [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
    real(r8),parameter :: wnd_min_dps=1.0 ! [m s-1] Minimum windspeed for momentum exchange
    ! Input
    integer, intent(in)  :: sfc_typ        ! [idx] LSM surface type (0..28)
    real(r8), intent(in) :: oro, &         ! [frc] Orography
                            snw_frc, &     ! [frc] Fraction of surface covered by snow
                            wnd_10m        ! [m s-1] 10 m wind speed
    ! Output
    real(r8), intent(out) :: rgh_mmn       ! [m] Roughness length momentum

    ! Local
    integer ice_nbr           ! [nbr] Number of sea ice points
    integer lnd_nbr           ! [nbr] Number of land points
    integer ocn_nbr           ! [nbr] Number of ocean points
    integer pln_typ_idx       ! [idx] Plant type index
    integer sgs_idx           ! [idx] Surface sub-gridscale index
    real(r8) rlm_crr              ! [m] Roughness length of current sub-gridscale
    real(r8) wnd_10m_bnd          ! [m s-1] Bounded wind speed at 10 m
    real(r8) xch_cff_mmn_ocn_ntr  ! [frc] Neutral 10 m drag coefficient over ocean
    !// ------------------------------------------------------------------
    ! Main Code
    
    ! Initialize
    rgh_mmn = 0.0_r8               ! [m]
    
    ! Initialize ocean, sea-ice, and land vectors
    ocn_nbr=0
    if (oro_is_ocn(oro)) ocn_nbr = 1

    ice_nbr=0
    if (oro_is_ice(oro)) ice_nbr = 1

    lnd_nbr=0
    if (oro_is_lnd(oro)) lnd_nbr = 1
    
    ! Ocean point
    if (ocn_nbr == 1) then
       ! Convert wind speed to roughness length over ocean
       wnd_10m_bnd = max(wnd_min_dps,wnd_10m) ! [m s-1]
       ! Approximation: neutral 10 m wind speed unavailable, use 10 m wind speed
       xch_cff_mmn_ocn_ntr = xch_cff_mmn_ocn_ntr_get(wnd_10m_bnd) ! [frc]
       rgh_mmn = 10.0_r8*exp(-cst_von_krm/sqrt(xch_cff_mmn_ocn_ntr)) ! [m] BKL97 p. F-4, LaP81 p. 327 (14) 
    endif
    
    ! Sea ice point
    if (ice_nbr == 1) then
       rgh_mmn = snw_frc * rgh_mmn_snw + (1.0_r8 - snw_frc)*rgh_mmn_ice_ocn ! [m] Bon96 p. 59
    endif
    
    ! Land point
    if (lnd_nbr == 1) then
       ! Store surface blend for current gridpoint
       if (sfc_typ == 0) then ! Inland lake
          ! fxm: Add temperature input and so ability to discriminate warm from frozen lakes here
          rgh_mmn = rgh_mmn_lak_wrm ! [m] Bon96 p. 59
       else if(sfc_typ == 1) then ! Land ice
          rgh_mmn = snw_frc * rgh_mmn_snw + (1.0_r8 - snw_frc)*rgh_mmn_ice_lnd ! [m] Bon96 p. 59
       else                   ! Normal land
          do sgs_idx=1,3
             ! Bare ground is pln_typ=14, ocean is pln_typ=0
             pln_typ_idx=pln_typ(sfc_typ,sgs_idx)
             if (pln_typ_idx == 14) then ! Bare ground
                rlm_crr = snw_frc * rgh_mmn_snw + (1.0_r8 - snw_frc)*rgh_mmn_ice_lnd ! [m] Bon96 p. 59 (glacial ice is same as bare ground)
             else if (pln_typ_idx > 0) then ! Regular plant type
                rlm_crr = snw_frc * rgh_mmn_snw + (1.0_r8 - snw_frc)*z0mvt(pln_typ_idx) ! [m] Bon96 p. 59
             else             ! Presumably ocean snuck through
                stop 'dst: ERROR: pln_typ_idx == 0 in rgh_mmn_get()'
             endif            ! endif
             rgh_mmn = & ! [m]
                  rgh_mmn + &
                  pln_frc(sfc_typ,sgs_idx)* & ! [frc]
                  rlm_crr     ! [m]
          end do               ! end loop over number of possible plant types
       endif                  ! endif normal land

       !++alfgr (Stop iteration for wnd_frc from going out of hand)
       !I choose to limit roughness lengths to 70 cm
       rgh_mmn = min(0.7_r8, rgh_mmn)
       !--alfgr

    endif
    
    ! Sanity check
    ! do lon_idx=1,plon
    ! if (rgh_mmn <= 0.0_r8) then
    ! write(6,'(a,3(a,es8.1),a)') &
    ! 'rgh_mmn_get: ','oro = ',oro,', rgh_mmn = ',rgh_mmn,' m, wnd_10m = ',wnd_10m,' m s-1'
    ! endif                  ! endif err
    ! end do                     ! end loop over ocn lon

    !// ------------------------------------------------------------------
  end subroutine rgh_mmn_get                       ! end rgh_mmn_get()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine snw_frc_get( &   ! [fnc] Convert snow depth to fractional snow cover
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       snw_frc              & ! O [frc] Fraction of surface covered by snow
       )
    !// ------------------------------------------------------------------
    ! Purpose: Convert equivalent liquid water snow depth to fractional snow cover
    ! Use snow thickness -> fraction algorithm of Bon96
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstblm, only: dns_H2O_lqd_std ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! Note disparity in bulk snow density between CCM and LSM
    ! WiW80 p. 2724, 2725 has some discussion of bulk snow density
    ! real(r8),parameter::dns_H2O_snw_gnd_LSM=250.0 ! [kg m-3] Bulk density of snow, CCM:lsm/snoconi.F
    real(r8),parameter :: dns_H2O_snw_gnd_std=100._r8; ! [kg m-3] Standard bulk density of snow on ground WiW80 p. 2724, 2725, CCM:physics/tsinti()
    real(r8),parameter :: snw_hgt_thr=0.05_r8 ! [m] Geometric snow thickness for 100% coverage, CCM:lsm/snoconi.F 

    ! Input
    real(r8),intent(in) :: snw_hgt_lqd   ! [m] Equivalent liquid water snow depth
    ! Output
    real(r8),intent(out) :: snw_frc      ! [frc] Fraction of surface covered by snow
    ! Local
    real(r8) :: snw_hgt         ! [m] Geometric bulk thickness of snow
    real(r8) :: hgt_lqd_snw_cnv ! [frc] Conversion factor from liquid water depth to geometric snow thickness
    !// ------------------------------------------------------------------
    ! Main code
    ! Initialize some variables
    hgt_lqd_snw_cnv = dns_H2O_lqd_std/dns_H2O_snw_gnd_std ! [frc] Conversion factor from liquid water depth to geometric snow thickness
    ! Fractional snow cover
    snw_hgt = snw_hgt_lqd * hgt_lqd_snw_cnv ! [m] NB: CCM and LSM seem to disagree on this
    snw_frc = min(snw_hgt/snw_hgt_thr,1.0_r8) ! [frc]
    !// ------------------------------------------------------------------
  end subroutine snw_frc_get                       ! end snw_frc_get()
  !// ------------------------------------------------------------------


  
  !// ------------------------------------------------------------------
  subroutine wnd_rfr_get( & ! [fnc] Interpolate wind speed to reference height
       flg_oro,             & ! I [flg] Orography flag
       hgt_mdp,             & ! I [m] Midpoint height above surface
       hgt_rfr,             & ! I [m] Reference height
       hgt_zpd,             & ! I [m] Zero plane displacement
       mno_lng,             & ! I [m] Monin-Obukhov length
       wnd_frc,             & ! I [m s-1] Surface friction velocity
       wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
       wnd_min,             & ! I [m s-1] Minimum windspeed
       wnd_rfr)             ! O [m s-1] Wind speed at reference height
    !// ------------------------------------------------------------------
    ! Purpose: Interpolate wind speed at given height to wind speed at reference height
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    integer,parameter::rfr_hgt_idx=1 ! Named index for lower (target) height
    integer,parameter::gcm_hgt_idx=2 ! Named index for upper (known) height
    ! Input
    logical,intent(in) :: flg_oro    ! I [flg] Orography flag
    real(r8),intent(in) :: hgt_mdp, &       ! I [m] Midlayer height above surface
                           hgt_rfr, &              ! I [m] Reference height
                           hgt_zpd, &       ! I [m] Zero plane displacement
                           mno_lng, &       ! I [m] Monin-Obukhov length
                           wnd_frc, &       ! I [m s-1] Surface friction velocity
                           wnd_mdp, &       ! I [m s-1] Surface layer mean wind speed
                           wnd_min              ! I [m s-1] Minimum windspeed
    ! Output
    real(r8),intent(out) :: wnd_rfr       ! O [m s-1] Wind speed at reference height

    ! Local
    integer lvl_idx           ! Stability computation loop index
    real(r8) mno_stb_crc_mmn(2) ! [frc] Monin-Obukhov stability correction momentum
    real(r8) mno_stb_prm(2) ! [frc] Monin-Obukhov stability parameter 
    real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
    real(r8) tmp2                 ! Term in stability correction computation
    real(r8) tmp3                 ! Term in stability correction computation
    real(r8) tmp4                 ! Term in stability correction computation
    real(r8) wnd_crc_fct   ! [frc] Wind correction factor
    real(r8) hgt_rfr_rcp      ! [m-1] Reciprocal of reference height
    !// ------------------------------------------------------------------
    ! Initialize
    hgt_rfr_rcp = 1.0_r8/hgt_rfr ! [m-1]
    ! Initialize output
    wnd_rfr = wnd_min           ! [m s-1]
    
    ! Compute horizontal wind speed at reference height
    if (flg_oro .and. hgt_zpd < hgt_rfr) then
       ! Code uses notation of Bon96 p. 50, where lvl_idx=1 is 10 m ref. hgt, lvl_idx=2 is atm. hgt
       mno_stb_prm(rfr_hgt_idx) = min((hgt_rfr - hgt_zpd)/mno_lng, 1.0_r8) ! [frc]
       mno_stb_prm(gcm_hgt_idx) = min((hgt_mdp - hgt_zpd)/mno_lng, 1.0_r8) ! [frc]
       do lvl_idx=1,2
          if (mno_stb_prm(lvl_idx) < 0.0_r8) then
             sml_fnc_mmn_uns_rcp = (1.0_r8 - 16.0_r8*mno_stb_prm(lvl_idx))**0.25_r8
             tmp2 = log((1.0_r8 + sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0_r8)
             tmp3 = log((1.0_r8+sml_fnc_mmn_uns_rcp)/2.0_r8)
             mno_stb_crc_mmn(lvl_idx) = 2.0_r8*tmp3+tmp2-2.0_r8*atan(sml_fnc_mmn_uns_rcp)+1.5707963
          else             ! not stable
             mno_stb_crc_mmn(lvl_idx)=-5.0_r8*mno_stb_prm(lvl_idx)
          endif            ! stable
       end do               ! end loop over lvl_idx
       tmp4 = log((hgt_mdp - hgt_zpd)/(hgt_rfr-hgt_zpd))
       ! Correct neutral stability assumption
       wnd_crc_fct = tmp4-mno_stb_crc_mmn(gcm_hgt_idx)+mno_stb_crc_mmn(rfr_hgt_idx) ! [frc]
       wnd_rfr = wnd_mdp - wnd_frc*cst_von_krm_rcp*wnd_crc_fct ! [m s-1]
       wnd_rfr = max(wnd_rfr,wnd_min) ! [m s-1]
       ! Use neutral stability assumption
       ! wnd_rfr=wnd_mdp-wnd_frc*cst_von_krm_rcp*tmp4
    endif                  ! endif (hgt_zpd < hgt_rfr.and.flg_oro)
    
    !// ------------------------------------------------------------------
  end subroutine wnd_rfr_get                       ! end wnd_rfr_get()
  !// ------------------------------------------------------------------




  !// ------------------------------------------------------------------
  subroutine zpd_get( &       ! [sbr] Set zero plane displacement
       hgt_zpd,             & ! O [m] Zero plane displacement
       oro,                 & ! I [frc] Orography
       sfc_typ,             & ! I [idx] LSM surface type (0..28)
       snw_frc)               ! I [frc] Fraction of surface covered by snow
    !// ------------------------------------------------------------------
    ! Purpose: Set zero plane displacement 
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstlsm, only: pln_typ, pln_frc, zpdvt! [mdl] LSM data
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! Input
    integer  :: sfc_typ       ! [idx] LSM surface type (0..28)
    real(r8) :: oro, &        ! [frc] Orography
                snw_frc       ! [frc] Fraction of surface covered by snow
    ! Output
    real(r8) :: hgt_zpd       ! [m] Zero plane displacement

    ! Local
    integer lnd_nbr           ! [nbr] Number of land points
    integer pln_typ_idx       ! [idx] Plant type index
    integer sgs_idx           ! [nbr] Surface sub-gridscale index
    !// ------------------------------------------------------------------
    ! Main Code
    
    ! Initialize array
    ! Zero plane displacement is identically 0.0_r8 everywhere except land
    hgt_zpd = 0.0_r8         ! [m]
    
    ! Land ahoy!
    lnd_nbr=0
    if (oro_is_lnd(oro)) lnd_nbr = 1


    ! Land point
    if (lnd_nbr == 1) then
       ! Store surface blend for current gridpoint
       if (sfc_typ == 0) then ! Inland lake
          hgt_zpd = 0.0_r8 ! [m]
       else if(sfc_typ == 1) then ! Land ice
          hgt_zpd = 0.0_r8 ! [m]
       else                   ! Normal land
          do sgs_idx=1,3
             ! Bare ground is pln_typ=14, ocean is pln_typ=0
             pln_typ_idx = pln_typ(sfc_typ,sgs_idx)
             hgt_zpd = & ! [m]
                  hgt_zpd + &
                  pln_frc(sfc_typ,sgs_idx)* & ! [frc]
                  zpdvt(pln_typ_idx) ! [m]
          end do               ! end loop over number of possible plant types
       endif                  ! endif normal land
    endif !// if (lnd_nbr == 1) then

    !// ------------------------------------------------------------------
  end subroutine zpd_get                       ! end zpd_get()
  !// ------------------------------------------------------------------




  !// ------------------------------------------------------------------
  subroutine blm_ocn( & ! [sbr] Compute boundary layer properties over ocean
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_ocn,             & ! I [flg] Ocean flag
       hgt_mdp,             & ! I [m] Midlayer height above surface
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
       tpt_mdp,             & ! I [K] Midlayer temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
       wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
       mno_lng,             & ! O [m] Monin-Obukhov length
       rgh_mmn,             & ! O [m] Roughness length momentum
       wnd_frc)             ! O [m s-1] Surface friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
    ! compute the boundary layer exchange properties over open ocean
    ! blm_ocn() is called by blm_glb()
    ! Routine uses specified surface temperature rather than solving energy balance equation for new Ts
    ! Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
    ! Coding: Zender
    
    ! Notes on algorithm:
    ! In this routine surface temperature tpt_sfc = Sea surface temperature SST
    ! Routine uses sfc suffix rather than SST simply for consistency with land and sea ice routines
    ! Suffix mdp quantity evaluated at height hgt_mdp
    ! Suffix sfc quantity evaluated at (sea) surface temperature
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstblm ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    integer,parameter::itr_max=5 ! [nbr] Maximum number of iterations for tpt_gnd loop
    real(r8),parameter::dlt_nbr_ntr_10m_cff=0.0346_r8 ! [frc] Coefficient for neutral 10 m Dalton number CCM:dom/flxoce() LaP82 p. 477
    real(r8),parameter::eps_max=1.0d-5 ! [frc] Relative accuracy for convergence
    real(r8),parameter::hgt_rfr_LaP81=10.0_r8 ! [m] Reference height for turbulent flux parameterization of LaP81
    real(r8),parameter::ssh_H2O_sln_crc=0.98_r8 ! [frc] Salinity correction to saturation specific humidity of H2O LaP81 p. 328 CCM:dom/flxoce()
    real(r8),parameter::stn_nbr_ntr_10m_cff_stb=0.0180_r8 ! [frc] Coefficient for neutral, stable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
    real(r8),parameter::stn_nbr_ntr_10m_cff_uns=0.0327_r8 ! [frc] Coefficient for neutral, unstable, 10 m Stanton number CCM:dom/flxoce() LaP82 p. 476
    real(r8),parameter::wnd_min_dps=1.0_r8 ! [m s-1] Minimum windspeed used for mobilization
    ! Input
    logical, intent(in) :: flg_ocn    ! I [flg] Ocean flag
    real(r8), intent(in) :: dns_mdp, &       ! I [kg m-3] Midlayer density
                            hgt_mdp, &       ! I [m] Midlayer height above surface
                            prs_mdp, &       ! I [Pa] Pressure
                            q_H2O_vpr, &     ! I [kg kg-1] Specific humidity
                            tpt_mdp, &       ! I [K] Midlayer temperature
                            tpt_ptn_mdp, &   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
                            tpt_sfc, &       ! I [K] Surface temperature
                            wnd_mrd_mdp, &   ! I [m s-1] Surface layer meridional wind speed
                            wnd_znl_mdp      ! I [m s-1] Surface layer zonal wind speed
    ! Input/Output
    ! Output
    real(r8), intent(out) :: mno_lng, &       ! O [m] Monin-Obukhov length
                             rgh_mmn, &       ! O [m] Roughness length momentum
                             wnd_frc          ! O [m s-1] Surface friction velocity
    ! Local
    integer itr_idx           ! [idx] Counting index
    real(r8) eps_crr              ! [frc] Current relative accuracy
    real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
    real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
    real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
    real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
    real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
    real(r8) mno_stb_prm          ! [frc] Monin-Obukhov stability parameter 
    real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
    real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
    real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
    real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
    real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
    real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
    real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
    real(r8) tpt_dlt              ! [K] Temperature change
    real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
    real(r8) tpt_scl              ! [K] Temperature scale
    real(r8) vpr_scl              ! [kg kg-1] Moisture scale
    real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
    real(r8) wnd_mdp              ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
    real(r8) wnd_rfr_ntr          ! [m s-1] Neutral 10 m wind speed
    real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
    real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
    real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
    real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
    real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
    real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
    real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
    real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange

#ifdef UNDEFINED
    ! Variables used in diagnostics
    real(r8) flx_LW_up_sfc        ! [W m-2] Longwave upwelling flux at surface
    real(r8) flx_ltn              ! [W m-2] Latent heat flux to atmosphere
    real(r8) flx_sns              ! [W m-2] Sensible heat flux to atmosphere
    real(r8) wnd_str              ! [kg m-1 s-2] Wind stress
    real(r8) wnd_str_mrd          ! [kg m-1 s-2] Meridional wind stress
    real(r8) wnd_str_znl          ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */
    
    !  Next three variables needed only for 10 m wind speed and mid-layer drag coefficient
    !  real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
    !  real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height LaP81 p. 327 (14)
    !  real(r8) xch_cff_mmn(plond)   ! [frc] Drag coefficient at mid-layer
    !// ------------------------------------------------------------------
    
    ! Main Code
    
    ! Initialize variables which are independent of stability iteration
    if (flg_ocn) then
          
       ! Midlayer wind speeds
       wnd_mdp =   & ! [m s-1] Surface layer mean wind speed
            sqrt(wnd_znl_mdp*wnd_znl_mdp +  &
            wnd_mrd_mdp*wnd_mrd_mdp) 
       wnd_mdp_bnd = max(wnd_mdp,wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded
          
       tpt_ptn_vrt_mdp = tpt_ptn_mdp*(1.0_r8+eps_H2O_rcp_m1*q_H2O_vpr) ! [K] Midlayer Virtual potential temperature
       tpt_bnd_cls = tpt_bnd_cls_get(tpt_sfc) ! [C]
       ! Liquid sea water supercools befored freezing to sea ice (CCM:dom/parsst.h/tsice = -1.7999 C)
       ! Assume gridpoint liquid sea water thus use heat of vaporization not sublimation
       svp_H2O_sfc = svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
       ltn_heat_trn = ltn_heat_vpr_H2O_std ! [J kg-1]
       ssh_H2O_sfc = eps_H2O*svp_H2O_sfc/(prs_mdp-one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
       ! Correct freshwater saturation specific humidity for salinity effects
       ssh_H2O_sfc = ssh_H2O_sfc*ssh_H2O_sln_crc ! [kg kg-1] Saturation specific humidity of H2O at surface
          
       tpt_dlt = tpt_ptn_mdp-tpt_sfc ! [K] Temperature change
       q_H2O_vpr_dlt = q_H2O_vpr-ssh_H2O_sfc ! [kg kg-1] Humidity change
       hgt_rat_log = log(hgt_mdp/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
       spc_heat_mst_air = spc_heat_dry_air*(1.0_r8+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
          
       ! Stable if tpt_ptn_mdp > tpt_sfc
       stb_val = 0.5_r8+sign(0.5_r8,tpt_dlt) ! [flg] 1.0 if stable, 0.0 if unstable
       ! Initial guess for roots of neutral exchange coefficients: z/L=0 and u10n=vmag
       xch_cff_mmn_ntr_sqrt = sqrt(xch_cff_mmn_ocn_ntr_get(wnd_mdp_bnd)) ! [frc] Squareroot of neutral 10 m drag coefficient
       xch_cff_heat_ntr_sqrt = (1.0_r8-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
       xch_cff_vpr_ntr_sqrt = dlt_nbr_ntr_10m_cff ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
          
#ifdef UNDEFINED
       if (.true.) then
          write (6,'(a)') 'Ocean fluxes:'
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               'itr','mno_lng','mno_stb','wnd_frc','  U10N ','  CDN  ','   CD  ',' H(atm)','   L   ',' LW(up)','  eps  '
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               '   ','   m   ','  frc  ',' m s-1 ',' m s-1 ',' x 1000',' x 1000',' W m-2 ',' W m-2 ',' W m-2 ','  frc  '
       endif               ! endif dbg
#endif /* not UNDEFINED */
          
       ! Iteration loop
       do itr_idx=1,itr_max
             
          ! Save old friction speed for convergence diagnostic
          if (itr_idx == 1) then
             eps_crr = eps_max+1.0_r8 ! [frc] Current relative accuracy
             wnd_frc_old = 0.0_r8 ! [m s-1]
          else 
             wnd_frc_old = wnd_frc ! [m s-1]
          endif            ! endif
             
          if (itr_idx == 1) then
             ! First iteration: Use estimated roots of neutral exchange coefficients
             wnd_frc = xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          else
             ! Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
             wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          endif            ! endif first iteration
             
          ! Compute stability parameter at midlayer and evaluate stability corrections  
          ! Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
          mno_stb_prm = & ! [frc] Monin-Obukhov stability parameter 
                  cst_von_krm*grv_sfc*hgt_mdp &
                  *(tpt_scl/tpt_ptn_vrt_mdp+vpr_scl/(1.0_r8/eps_H2O_rcp_m1+q_H2O_vpr)) &
                  /(wnd_frc*wnd_frc)
          ! Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
          mno_stb_prm = sign(min(abs(mno_stb_prm),10.0_r8),mno_stb_prm) ! [frc] Monin-Obukhov stability parameter
          stb_val = 0.5_r8+sign(0.5_r8,mno_stb_prm) ! [flg] 1.0 if stable, 0.0 if unstable
          sml_fnc_mmn_uns_rcp = max(sqrt(abs(1.0_r8-16.0_r8*mno_stb_prm)),1.0_r8) ! [frc] BKL97 p. F1, LaP81 p. 325
          sml_fnc_mmn_uns_rcp = sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_mmn = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_heat = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
             
          ! Shift old neutral 10 m exchange coefficient to measurement height and stability 
          xch_cff_mmn_sqrt = xch_cff_mmn_ntr_sqrt/(1.0_r8+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
             
          ! Define neutral 10 m wind speed 
          wnd_rfr_ntr = wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt ! [m s-1] Neutral 10 m wind speed
             
          ! Update neutral 10 m exchange coefficients
          xch_cff_mmn_ntr_sqrt = sqrt(xch_cff_mmn_ocn_ntr_get(wnd_rfr_ntr)) ! [frc] Squareroot of neutral 10 m drag coefficient
          xch_cff_vpr_ntr_sqrt = dlt_nbr_ntr_10m_cff ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
          xch_cff_heat_ntr_sqrt = (1.0_r8-stb_val)*stn_nbr_ntr_10m_cff_uns+stb_val*stn_nbr_ntr_10m_cff_stb ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
             
          ! Shift old neutral 10 m exchange coefficients to measurement height and stability
          xch_cff_mmn_sqrt = xch_cff_mmn_ntr_sqrt/(1.0_r8+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
          xch_cff_heat_sqrt = xch_cff_heat_ntr_sqrt/(1.0_r8+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
          xch_cff_vpr_sqrt = xch_cff_vpr_ntr_sqrt/(1.0_r8+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
             
#ifdef UNDEFINED
          if (.true.) then
             eps_crr = abs((wnd_frc-wnd_frc_old)/wnd_frc) ! Relative convergence
             mno_lng = hgt_mdp/mno_stb_prm ! [m] Monin-Obukhov length
             wnd_str = dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
             flx_sns = -spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
             flx_ltn = -ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
             flx_LW_up_sfc = cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
             write (6,'(i3,1x,f9.3,1x,f8.3,1x,f7.4,1x,f7.3,1x,f7.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.6)') &
                  itr_idx,mno_lng,mno_stb_prm, &
                  wnd_frc,wnd_rfr_ntr, &
                  1000.0_r8*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, &
                  1000.0_r8*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt, &
                  flx_sns,flx_ltn,flx_LW_up_sfc,eps_crr
          endif            ! endif dbg
#endif /* not UNDEFINED */
             
       end do               ! end loop over itr
          
       ! Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
       wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
       tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
       vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 
          
       ! Compute surface stress components
       !        wnd_str=dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
       !        wnd_str_znl=-wnd_str*wnd_znl_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
       !        wnd_str_mrd=-wnd_str*wnd_mrd_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
          
       ! Compute heat flux components at current surface temperature
       ! Define positive latent and sensible heat as upwards into atmosphere 
       !        flx_sns=-spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
       !        flx_ltn=-ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
       !        flx_LW_up_sfc=cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
          
       ! Additional diagnostics
       mno_lng = hgt_mdp/mno_stb_prm ! [m] Monin-Obukhov length
       rgh_mmn = hgt_rfr_LaP81*exp(-cst_von_krm/xch_cff_mmn_ntr_sqrt) ! [m] BKL97 p. F-4, LaP81 p. 327 (14) 
          
       !       xch_cff_mmn=xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt ! [frc] Drag coefficient at mid-layer
          
    endif                  ! endif flg_ocn

    
    ! Initialize hgt_zpd to 0.0 m for ocean
    !  hgt_zpd=0.0_r8               ! [m]
    ! Interpolate midlayer wind speed to 10 m
    !  call wnd_rfr_get( &
    !      flg_ocn,             & ! I [flg] Ocean flag
    !      hgt_mdp,             & ! I [m] Midpoint height above surface
    !      hgt_rfr_LaP81,       & ! I [m] Reference height for turbulent flux parameterization of LaP81
    !      hgt_zpd,             & ! I [m] Zero plane displacement
    !      mno_lng,             & ! I [m] Monin-Obukhov length
    !      wnd_frc,             & ! I [m s-1] Surface friction velocity
    !      wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
    !      wnd_min_dps,         & ! I [m s-1] Minimum windspeed
    !      wnd_rfr)             ! O [m s-1] Wind speed at reference height
    

    !// ------------------------------------------------------------------
  end subroutine blm_ocn                       ! end blm_ocn()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine blm_ice( & ! [sbr] Compute boundary layer properties over ice
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_ice,             & ! I [flg] Sea ice flag
       hgt_mdp,             & ! I [m] Midlayer height above surface
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
       rgh_mmn,             & ! I [m] Roughness length momentum
       tpt_mdp,             & ! I [K] Midlayer temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
       wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
       mno_lng,             & ! O [m] Monin-Obukhov length
       wnd_frc)             ! O [m s-1] Surface friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
    ! compute the boundary layer exchange properties over sea ice
    ! blm_ice() is called by blm_glb()
    ! Routine uses specified surface temperature rather than solving energy balance equation for new Ts
    ! Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
    ! Coding: Zender
    
    ! Notes on algorithm:
    ! Suffix mdp quantity evaluated at height hgt_mdp
    ! Suffix sfc quantity evaluated at (sea) surface temperature
    ! The roughness length of sea ice depends on many factors
    ! A globally uniform roughness length for sea ice is foolish but necessary
    ! CCM1/2/3 used z0m=0.04 m (CCM:dom/parpbl), based on a glacial ice value
    ! LSM also adopted this value for land ice
    ! This value is more appropriate for ridged, multi-year ice
    ! The NCAR Oceanography section used z0m=0.05 m (BKL97 p. F-3)
    ! These values lead to large drag coefficients and excessive ice extent off of Antarctica
    ! CSM later adopted z0m=0.0005 m, appropriate for new-formed, seasonal ice
    ! This value greatly improved sea-ice dynamics in CSM, so we adopt it
    ! This blm_ice() routine reads roughness length as an input parameter
    ! for algorithmic reasons, and checks that its value agrees with CCM sea ice
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstblm ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    integer,parameter::itr_max=5 ! [nbr] Maximum number of iterations for tpt_gnd loop
    real(r8),parameter::eps_max=1.0d-5 ! [frc] Relative accuracy for convergence
    real(r8),parameter::hgt_rfr_LaP81=10.0_r8 ! [m] Reference height for turbulent flux parameterization of LaP81
    real(r8),parameter::rgh_mmn_ice_ocn=0.0005_r8 ! [m] Roughness length over sea ice BKL97 p. F-3
    real(r8),parameter::wnd_min_dps=1.0_r8 ! [m s-1] Minimum windspeed used for mobilization
    ! Input
    logical, intent(in) :: flg_ice           ! I [flg] Sea ice flag
    real(r8), intent(in) :: dns_mdp, &       ! I [kg m-3] Midlayer density
                            hgt_mdp, &       ! I [m] Midlayer height above surface
                            prs_mdp, &       ! I [Pa] Pressure
                            q_H2O_vpr, &     ! I [kg kg-1] Specific humidity
                            rgh_mmn, &       ! I [m] Roughness length momentum
                            tpt_mdp, &       ! I [K] Midlayer temperature
                            tpt_ptn_mdp, &   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
                            tpt_sfc, &       ! I [K] Surface temperature
                            wnd_mrd_mdp, &   ! I [m s-1] Surface layer meridional wind speed
                            wnd_znl_mdp      ! I [m s-1] Surface layer zonal wind speed
    ! Input/Output
    ! Output
    real(r8), intent(out) :: mno_lng, &       ! O [m] Monin-Obukhov length
                             wnd_frc          ! O [m s-1] Surface friction velocity
    ! Local
    integer itr_idx           ! [idx] Counting index
    real(r8) eps_crr              ! [frc] Current relative accuracy
    real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
    real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
    real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
    real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
    real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
    real(r8) mno_stb_prm   ! [frc] Monin-Obukhov stability parameter 
    real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
    real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
    real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
    real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
    real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
    real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
    real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
    real(r8) tpt_dlt              ! [K] Temperature change
    real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
    real(r8) tpt_scl       ! [K] Temperature scale
    real(r8) vpr_scl       ! [kg kg-1] Moisture scale
    real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
    real(r8) wnd_mdp       ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
    real(r8) wnd_rfr_ntr   ! [m s-1] Neutral 10 m wind speed
    real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
    real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
    real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
    real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
    real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
    real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
    real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
    real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange
    
#ifdef UNDEFINED
    ! Variables used in diagnostics
    real(r8) flx_LW_up_sfc ! [W m-2] Longwave upwelling flux at surface
    real(r8) flx_ltn       ! [W m-2] Latent heat flux to atmosphere
    real(r8) flx_sns       ! [W m-2] Sensible heat flux to atmosphere
    real(r8) wnd_str       ! [kg m-1 s-2] Wind stress
    real(r8) wnd_str_mrd   ! [kg m-1 s-2] Meridional wind stress
    real(r8) wnd_str_znl   ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */
    !// ------------------------------------------------------------------
   
    ! Main Code
    
    ! Initialize variables which are independent of stability iteration
    if (flg_ice) then
          
       ! Midlayer wind speeds
       wnd_mdp =   & ! [m s-1] Surface layer mean wind speed
            sqrt(wnd_znl_mdp*wnd_znl_mdp+  &
            wnd_mrd_mdp*wnd_mrd_mdp) 
       wnd_mdp_bnd = max(wnd_mdp,wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded
          
       tpt_ptn_vrt_mdp = tpt_ptn_mdp*(1.0_r8+eps_H2O_rcp_m1*q_H2O_vpr) ! [K] Midlayer Virtual potential temperature
       tpt_bnd_cls = tpt_bnd_cls_get(tpt_sfc) ! [C]
       ! Assume gridpoint is sea ice thus use heat of sublimation not vaporization (contrary to CSM)
       ! NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
       svp_H2O_sfc = svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
       ltn_heat_trn = ltn_heat_sbl_H2O_std ! [J kg-1]
       ssh_H2O_sfc = eps_H2O*svp_H2O_sfc/(prs_mdp-one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
       ! NB: No salinity correction for sea ice (contrary to CSM)
          
       tpt_dlt = tpt_ptn_mdp-tpt_sfc ! [K] Temperature change
       q_H2O_vpr_dlt = q_H2O_vpr-ssh_H2O_sfc ! [kg kg-1] Humidity change
       hgt_rat_log = log(hgt_mdp/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
       spc_heat_mst_air = spc_heat_dry_air*(1.0_r8+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
          
#ifdef UNDEFINED
       if (.true.) then
          write (6,'(a)') 'Sea ice fluxes:'
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               'itr','mno_lng','mno_stb','wnd_frc','  U10N ','  CDN  ','   CD  ',' H(atm)','   L   ',' LW(up)','  eps  '
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               '   ','   m   ','  frc  ',' m s-1 ',' m s-1 ',' x 1000',' x 1000',' W m-2 ',' W m-2 ',' W m-2 ','  frc  '
       endif               ! endif dbg
#endif /* not UNDEFINED */
          
       ! Iteration loop
       do itr_idx=1,itr_max
             
          ! Save old friction speed for convergence diagnostic
          if (itr_idx == 1) then
             eps_crr = eps_max+1.0_r8 ! [frc] Current relative accuracy
             wnd_frc_old = 0.0_r8 ! [m s-1]
          else 
             wnd_frc_old = wnd_frc ! [m s-1]
          endif            ! endif
             
          ! Roots of neutral, 10 m exchange coefficients
          xch_cff_mmn_ntr_sqrt = cst_von_krm/log(hgt_rfr_LaP81/rgh_mmn) ! [frc] Squareroot of neutral 10 m drag coefficient
          xch_cff_heat_ntr_sqrt = xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
          xch_cff_vpr_ntr_sqrt = xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
             
          if (itr_idx == 1) then
             ! First iteration: Use estimated roots of neutral exchange coefficients
             wnd_frc = xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          else
             ! Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
             wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          endif            ! endif first iteration
             
          ! Compute stability parameter at midlayer and evaluate stability corrections  
          ! Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
          mno_stb_prm = & ! [frc] Monin-Obukhov stability parameter 
               cst_von_krm*grv_sfc*hgt_mdp &
               *(tpt_scl/tpt_ptn_vrt_mdp+vpr_scl/(1.0_r8/eps_H2O_rcp_m1+q_H2O_vpr)) &
               /(wnd_frc*wnd_frc)
          ! Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
          mno_stb_prm = sign(min(abs(mno_stb_prm),10.0_r8),mno_stb_prm) ! [frc] Monin-Obukhov stability parameter
          stb_val = 0.5_r8+sign(0.5_r8,mno_stb_prm) ! [flg] 1.0 if stable, 0.0 if unstable
          sml_fnc_mmn_uns_rcp = max(sqrt(abs(1.0_r8-16.0_r8*mno_stb_prm)),1.0_r8) ! [frc] BKL97 p. F1, LaP81 p. 325
          sml_fnc_mmn_uns_rcp = sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_mmn = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_heat = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
             
          ! Shift old neutral 10 m exchange coefficients to measurement height and stability
          xch_cff_mmn_sqrt=xch_cff_mmn_ntr_sqrt/(1.0_r8+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
          xch_cff_heat_sqrt=xch_cff_heat_ntr_sqrt/(1.0_r8+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
          xch_cff_vpr_sqrt=xch_cff_vpr_ntr_sqrt/(1.0_r8+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
             
#ifdef UNDEFINED
          if (.true.) then
             eps_crr=abs((wnd_frc-wnd_frc_old)/wnd_frc) ! Relative convergence
             wnd_rfr_ntr=wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt ! [m s-1] Neutral 10 m wind speed (by definition)
             mno_lng=hgt_mdp/mno_stb_prm ! [m] Monin-Obukhov length
             wnd_str=dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
             flx_sns=-spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
             flx_ltn=-ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
             flx_LW_up_sfc=cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
             write (6,'(i3,1x,f9.3,1x,f8.3,1x,f7.4,1x,f7.3,1x,f7.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,f8.6)') &
                  itr_idx,mno_lng,mno_stb_prm, &
                  wnd_frc,wnd_rfr_ntr, &
                  1000.0_r8*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, &
                  1000.0_r8*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt, &
                  flx_sns,flx_ltn,flx_LW_up_sfc,eps_crr
          endif            ! endif dbg
#endif /* not UNDEFINED */
             
       end do               ! end loop over itr
          
       ! Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
       wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
       tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
       vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 
          
       ! Compute surface stress components
       !        wnd_str=dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
       !        wnd_str_znl=-wnd_str*wnd_znl_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
       !        wnd_str_mrd=-wnd_str*wnd_mrd_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
          
       ! Compute heat flux components at current surface temperature
       ! Define positive latent and sensible heat as upwards into atmosphere 
       !        flx_sns=-spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
       !        flx_ltn=-ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
       !        flx_LW_up_sfc=cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
       
       ! Additional diagnostics
       mno_lng = hgt_mdp/mno_stb_prm ! [m] Monin-Obukhov length
    endif                  ! endif flg_ice
    

    !// ------------------------------------------------------------------
  end subroutine blm_ice                       ! end blm_ice()
  !// ------------------------------------------------------------------


  
  !// ------------------------------------------------------------------
  subroutine blm_lnd( & ! [sbr] Compute boundary layer properties over land
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_lnd,             & ! I [flg] Land flag
       hgt_mdp,             & ! I [m] Midlayer height above surface
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
       rgh_mmn,             & ! I [m] Roughness length momentum
       tpt_mdp,             & ! I [K] Midlayer temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
       wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
       mno_lng,             & ! O [m] Monin-Obukhov length
       wnd_frc)             ! O [m s-1] Surface friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
    ! compute the boundary layer exchange properties over land
    ! blm_lnd() is called by blm_glb()
    ! Routine uses specified surface temperature rather than solving energy balance equation for new Ts
    ! Theory and algorithms: Large and Pond (1981,1982) CCM:dom/flxoce(), Bonan (1996) CCM:lsm/surtem()
    ! Coding: Zender
    
    ! Notes on algorithm:
    ! Suffix mdp quantity evaluated at height hgt_mdp
    ! Suffix sfc quantity evaluated at (sea) surface temperature
    ! Currently, routine is virtually identical to blm_ice()
    ! fxm: Using blm_ice()/blm_ocn()-type physics is unjustifiable over
    ! heterogeneous land surface types because of varying transpiration
    ! resistances over heterogeneous land surface types
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstblm ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    integer,parameter  :: itr_max=5      ! [nbr] Maximum number of iterations for tpt_gnd loop
    real(r8),parameter :: eps_max=1.0e-5 ! [frc] Relative accuracy for convergence
    real(r8),parameter :: hgt_rfr_LaP81=10.0_r8 ! [m] Reference height for turbulent flux parameterization of LaP81
    real(r8),parameter :: wnd_min_dps=1.0_r8 ! [m s-1] Minimum windspeed used for mobilization

    ! Input
    logical, intent(in) :: flg_lnd          ! I [flg] Land flag
    real(r8), intent(in) :: dns_mdp, &      ! I [kg m-3] Midlayer density
                            hgt_mdp, &      ! I [m] Midlayer height above surface
                            prs_mdp, &      ! I [Pa] Pressure
                            q_H2O_vpr, &    ! I [kg kg-1] Specific humidity
                            rgh_mmn, &      ! I [m] Roughness length momentum
                            tpt_mdp, &      ! I [K] Midlayer temperature
                            tpt_ptn_mdp, &  ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
                            tpt_sfc, &      ! I [K] Surface temperature
                            wnd_mrd_mdp, &  ! I [m s-1] Surface layer meridional wind speed
                            wnd_znl_mdp     ! I [m s-1] Surface layer zonal wind speed
    ! Input/Output
    ! Output
    real(r8), intent(out) :: mno_lng, &       ! O [m] Monin-Obukhov length
                             wnd_frc          ! O [m s-1] Surface friction velocity
    ! Local
    integer itr_idx           ! [idx] Counting index
    real(r8) eps_crr              ! [frc] Current relative accuracy
    real(r8) hgt_rat_log          ! [frc] Log of ratio of local to reference height
    real(r8) ltn_heat_trn         ! [J kg-1] Latent heat of sublimation or evaporation
    real(r8) mno_stb_crc_heat     ! [frc] Monin-Obukhov stability correction heat
    real(r8) mno_stb_crc_mmn      ! [frc] Monin-Obukhov stability correction momentum
    real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
    real(r8) mno_stb_prm          ! [frc] Monin-Obukhov stability parameter 
    real(r8) ntp_fct              ! [frc] Interpolation factor in reference height temperature calculation
    real(r8) q_H2O_vpr_dlt        ! [kg kg-1] Humidity change
    real(r8) spc_heat_mst_air     ! [J kg-1 K-1] Specific heat of moist air
    real(r8) ssh_H2O_sfc          ! [kg kg-1] Saturation specific humidity of H2O at surface
    real(r8) stb_val              ! [flg] 1.0 if stable, 0.0 if unstable
    real(r8) svp_H2O_sfc          ! [Pa] Saturation vapor pressure over planar condensed water at surface
    real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
    real(r8) tpt_dlt              ! [K] Temperature change
    real(r8) tpt_ptn_vrt_mdp      ! [K] Midlayer virtual potential temperature
    real(r8) tpt_scl              ! [K] Temperature scale
    real(r8) vpr_scl              ! [kg kg-1] Moisture scale
    real(r8) wnd_frc_old          ! [m s-1] Surface friction velocity old
    real(r8) wnd_mdp              ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_mdp_bnd          ! [m s-1] Surface layer mean wind speed bounded
    real(r8) wnd_rfr_ntr          ! [m s-1] Neutral 10 m wind speed
    real(r8) xch_cff_heat_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
    real(r8) xch_cff_heat_sqrt    ! [frc] Squareroot of mid-layer Stanton number for heat exchange
    real(r8) xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m drag coefficient
    real(r8) xch_cff_mmn_sqrt     ! [frc] Squareroot of mid-layer drag coefficient
    real(r8) xch_cff_vpr_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
    real(r8) xch_cff_vpr_sqrt     ! [frc] Squareroot of mid-layer Dalton number for vapor exchange
    real(r8) xpn_heat_fct         ! [frc] Factor of heat exchange
    real(r8) xpn_mmn_fct          ! [frc] Factor of momentum exchange
    
#ifdef UNDEFINED
    ! Variables used in diagnostics
    real(r8) flx_LW_up_sfc ! [W m-2] Longwave upwelling flux at surface
    real(r8) flx_ltn       ! [W m-2] Latent heat flux to atmosphere
    real(r8) flx_sns       ! [W m-2] Sensible heat flux to atmosphere
    real(r8) wnd_str       ! [kg m-1 s-2] Wind stress
    real(r8) wnd_str_mrd   ! [kg m-1 s-2] Meridional wind stress
    real(r8) wnd_str_znl   ! [kg m-1 s-2] Zonal wind stress
#endif /* not UNDEFINED */
    !// ------------------------------------------------------------------
    
    ! Main Code
    
    ! Initialize variables which are independent of stability iteration
    if (flg_lnd) then

       ! Midlayer wind speeds
       wnd_mdp =   & ! [m s-1] Surface layer mean wind speed
            sqrt(wnd_znl_mdp*wnd_znl_mdp +  &
            wnd_mrd_mdp*wnd_mrd_mdp) 
       wnd_mdp_bnd = max(wnd_mdp,wnd_min_dps) ! [m s-1] Surface layer mean wind speed bounded

       tpt_ptn_vrt_mdp = tpt_ptn_mdp * (1.0_r8 + eps_H2O_rcp_m1*q_H2O_vpr) ! [K] Midlayer Virtual potential temperature
       ! Saturation vapor pressure of water at ground temperature
       ! NB: Prognostic models must sometimes use heat of vaporization everywhere to conserve energy
       tpt_bnd_cls = tpt_bnd_cls_get(tpt_sfc) ! [C]
       if (tpt_bnd_cls > 0.0_r8) then 
          svp_H2O_sfc = svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
          ltn_heat_trn = ltn_heat_vpr_H2O_std ! [J kg-1]
       else
          svp_H2O_sfc = svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
          ltn_heat_trn = ltn_heat_sbl_H2O_std ! [J kg-1]
       endif               ! endif frozen
       ssh_H2O_sfc = eps_H2O*svp_H2O_sfc/(prs_mdp - one_mns_eps_H2O*svp_H2O_sfc) ! [kg kg-1] Saturation specific humidity of H2O at surface
       ! NB: No salinity correction for land
          
       tpt_dlt = tpt_ptn_mdp - tpt_sfc ! [K] Temperature change
       q_H2O_vpr_dlt = q_H2O_vpr - ssh_H2O_sfc ! [kg kg-1] Humidity change
       hgt_rat_log = log(hgt_mdp/hgt_rfr_LaP81) ! [frc] Log of ratio of local to reference height
       spc_heat_mst_air = spc_heat_dry_air*(1.0_r8+cp_vpr_rcp_cp_dry_m1*ssh_H2O_sfc) ! [J kg-1 K-1] Specific heat of moist air
          
#ifdef UNDEFINED
       if (.true.) then
          write (6,'(a)') 'Land fluxes:'
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               'itr','mno_lng','mno_stb','wnd_frc','  U10N ','  CDN  ','   CD  ',' H(atm)','   L   ',' LW(up)','  eps  '
          write (6,'(a3,1x,a9,1x,a8,1x,a7,1x,a7,1x,a7,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
               '   ','   m   ','  frc  ',' m s-1 ',' m s-1 ',' x 1000',' x 1000',' W m-2 ',' W m-2 ',' W m-2 ','  frc  '
       endif               ! endif dbg
#endif /* not UNDEFINED */
          
       ! Iteration loop
       do itr_idx=1,itr_max
             
          ! Save old friction speed for convergence diagnostic
          if (itr_idx == 1) then
             eps_crr = eps_max+1.0_r8 ! [frc] Current relative accuracy
             wnd_frc_old = 0.0_r8 ! [m s-1]
          else 
             wnd_frc_old = wnd_frc ! [m s-1]
          endif            ! endif
             
          ! Roots of neutral, 10 m exchange coefficients
          xch_cff_mmn_ntr_sqrt = cst_von_krm/log(hgt_rfr_LaP81/rgh_mmn) ! [frc] Squareroot of neutral 10 m drag coefficient
          xch_cff_heat_ntr_sqrt = xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Stanton number for heat exchange
          xch_cff_vpr_ntr_sqrt = xch_cff_mmn_ntr_sqrt ! [frc] Squareroot of neutral 10 m Dalton number for vapor exchange
             
          if (itr_idx == 1) then
             ! First iteration: Use estimated roots of neutral exchange coefficients
             wnd_frc = xch_cff_mmn_ntr_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_ntr_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_ntr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          else
             ! Subsequently: Use current roots of (non-neutral) exchange coefficients at measurement height
             wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity
             tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale
             vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale
          endif            ! endif first iteration
             
          ! Compute stability parameter at midlayer and evaluate stability corrections  
          ! Stable if tpt_ptn_mdp > tpt_sfc or if mno_stb_prm > 0
          mno_stb_prm = & ! [frc] Monin-Obukhov stability parameter 
                  cst_von_krm*grv_sfc*hgt_mdp &
                  *(tpt_scl/tpt_ptn_vrt_mdp+vpr_scl/(1.0_r8/eps_H2O_rcp_m1+q_H2O_vpr)) &
                  /(wnd_frc*wnd_frc)
          ! Ensure |z/L| < 10.0 because similarity function asymptotes to (-z/L)^(-1/3) for z/L << -1, Ary88 p. 166
          mno_stb_prm = sign(min(abs(mno_stb_prm),10.0_r8),mno_stb_prm) ! [frc] Monin-Obukhov stability parameter
          stb_val = 0.5_r8+sign(0.5_r8,mno_stb_prm) ! [flg] 1.0 if stable, 0.0 if unstable
          sml_fnc_mmn_uns_rcp = max(sqrt(abs(1.0_r8-16.0_r8*mno_stb_prm)),1.0_r8) ! [frc] BKL97 p. F1, LaP81 p. 325
          sml_fnc_mmn_uns_rcp = sqrt(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_mmn = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_mmn_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
          mno_stb_crc_heat = -5.0_r8*mno_stb_prm*stb_val+(1.0_r8-stb_val)*mno_stb_crc_heat_uns_get(sml_fnc_mmn_uns_rcp) ! [frc] BKL97 p. F1, LaP81 p. 325
             
          ! Shift old neutral 10 m exchange coefficients to measurement height and stability
          xch_cff_mmn_sqrt = xch_cff_mmn_ntr_sqrt/(1.0_r8+xch_cff_mmn_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_mmn)) ! [frc] Squareroot of mid-layer drag coefficient LaP81 p. 327 (15), LaP82 p. 466 (10)
          xch_cff_heat_sqrt = xch_cff_heat_ntr_sqrt/(1.0_r8+xch_cff_heat_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Stanton number for heat exchange LaP82 p. 466 (10)
          xch_cff_vpr_sqrt = xch_cff_vpr_ntr_sqrt/(1.0_r8+xch_cff_vpr_ntr_sqrt/cst_von_krm*(hgt_rat_log-mno_stb_crc_heat)) ! [frc] Squareroot of mid-layer Dalton number for vapor exchange LaP82 p. 466 (10)
             
#ifdef UNDEFINED
          if (.true.) then
             eps_crr = abs((wnd_frc-wnd_frc_old)/wnd_frc) ! Relative convergence
             wnd_rfr_ntr = wnd_mdp_bnd*xch_cff_mmn_sqrt/xch_cff_mmn_ntr_sqrt ! [m s-1] Neutral 10 m wind speed (by definition)
             mno_lng = hgt_mdp/mno_stb_prm ! [m] Monin-Obukhov length
             wnd_str = dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
             flx_sns = -spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
             flx_ltn = -ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
             flx_LW_up_sfc = cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
             write (6,'(i3,1x,f9.3,1x,f8.3,1x,f7.4,1x,f7.3,1x,f7.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.6)') &
                  itr_idx,mno_lng,mno_stb_prm, &
                  wnd_frc,wnd_rfr_ntr, &
                  1000.0_r8*xch_cff_mmn_ntr_sqrt*xch_cff_mmn_ntr_sqrt, &
                  1000.0_r8*xch_cff_mmn_sqrt*xch_cff_mmn_sqrt, &
                  flx_sns,flx_ltn,flx_LW_up_sfc,eps_crr
          endif            ! endif dbg
#endif /* not UNDEFINED */
             
       end do               ! end loop over itr
          
       ! Update wnd_frc, tpt_scl, and vpr_scl using updated exchange coefficients at measurement height
       wnd_frc = xch_cff_mmn_sqrt*wnd_mdp_bnd ! [m s-1] Surface friction velocity 
       tpt_scl = xch_cff_heat_sqrt*tpt_dlt ! [K] Temperature scale 
       vpr_scl = xch_cff_vpr_sqrt*q_H2O_vpr_dlt ! [kg kg-1] Moisture scale 
          
       ! Compute surface stress components
       !        wnd_str=dns_mdp*wnd_frc*wnd_frc ! [kg m-1 s-2] Wind stress
       !        wnd_str_znl=-wnd_str*wnd_znl_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Zonal wind stress
       !        wnd_str_mrd=-wnd_str*wnd_mrd_mdp/wnd_mdp_bnd ! [kg m-1 s-2] Meridional wind stress
          
       ! Compute heat flux components at current surface temperature
       ! Define positive latent and sensible heat as upwards into atmosphere 
       !        flx_sns=-spc_heat_mst_air*wnd_str*tpt_scl/wnd_frc ! [W m-2] Sensible heat flux to atmosphere
       !        flx_ltn=-ltn_heat_trn*wnd_str*vpr_scl/wnd_frc ! [W m-2] Latent heat flux to atmosphere
       !        flx_LW_up_sfc=cst_Stefan_Boltzmann*(tpt_sfc**4.0_r8) ! [W m-2] Longwave upwelling flux to atmosphere
          
       ! Additional diagnostics
       mno_lng = hgt_mdp /mno_stb_prm ! [m] Monin-Obukhov length
    endif                  ! endif flg_lnd

    !// ------------------------------------------------------------------
  end subroutine blm_lnd                       ! end blm_lnd()
  !// ------------------------------------------------------------------
  



  !// ------------------------------------------------------------------
  subroutine blm_glb(       & ! [sbr] Solve boundary layer meteorology on global scale
       dns_mdp,             & ! I [kg m-3] Midlayer density
       hgt_mdp,             & ! I [m] Midlayer height above surface
       oro,                 & ! I [frc] Orography
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
       sfc_typ,             & ! I [idx] LSM surface type (0..28)
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       tpt_mdp,             & ! I [K] Midlayer temperature
       tpt_ptn_mdp,         & ! I [K] Potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
       wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
       hgt_zpd,             & ! O [m] Zero plane displacement
       mno_lng,             & ! O [m] Monin-Obukhov length
       rgh_mmn,             & ! O [m] Roughness length momentum
       snw_frc,             & ! O [frc] Fraction of surface covered by snow
       wnd_frc,             & ! O [m s-1] Surface friction velocity
       wnd_rfr)             ! O [m s-1] Wind speed at reference height
    !// ------------------------------------------------------------------
    ! Purpose: Given vectors of state variables and surface type, 
    ! divide vectors in land, ocean, and sea-ice components and call the routines 
    ! to compute the boundary layer exchange properties
    ! blm_glb() is called by dstdpsdry()
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    real(r8),parameter::hgt_rfr=10.0_r8   ! [m] Reference height for deposition processes
    real(r8),parameter::wnd_min_dps=1.0_r8 ! [m s-1] Minimum windspeed for momentum exchange
    ! Input
    integer, intent(in)  :: sfc_typ     ! I [idx] LSM surface type (0..28)
    real(r8), intent(in) :: oro         ! I [frc] Orography
    real(r8), intent(in) :: dns_mdp     ! I [kg m-3] Midlayer density
    real(r8), intent(in) :: hgt_mdp     ! I [m] Midlayer height above surface
    real(r8), intent(in) :: prs_mdp     ! I [Pa] Pressure
    real(r8), intent(in) :: q_H2O_vpr   ! I [kg kg-1] Specific humidity
    real(r8), intent(in) :: snw_hgt_lqd ! I [m] Equivalent liquid water snow depth
    real(r8), intent(in) :: tpt_mdp     ! I [K] Midlayer temperature
    real(r8), intent(in) :: tpt_ptn_mdp ! I [K] Midlayer local potential temperature
    real(r8), intent(in) :: tpt_sfc     ! I [K] Surface temperature
    real(r8), intent(in) :: wnd_mrd_mdp ! I [m s-1] Surface layer meridional wind speed
    real(r8), intent(in) :: wnd_znl_mdp ! I [m s-1] Surface layer zonal wind speed
    ! Input/Output
    ! Output
    real(r8), intent(out) :: hgt_zpd   ! O [m] Zero plane displacement
    real(r8), intent(out) :: mno_lng   ! O [m] Monin-Obukhov length
    real(r8), intent(out) :: rgh_mmn   ! O [m] Roughness length momentum
    real(r8), intent(out) :: snw_frc   ! O [frc] Fraction of surface covered by snow
    real(r8), intent(out) :: wnd_frc   ! O [m s-1] Surface friction velocity
    real(r8), intent(out) :: wnd_rfr   ! O [m s-1] Wind speed at reference height
    ! Local
    !integer lon_idx           ! [idx] Counting index for lon
    logical flg_ice    ! [flg] Sea ice flag
    logical flg_lnd    ! [flg] Land flag
    logical flg_ocn    ! [flg] Ocean flag
    logical flg_true   ! [flg] True flag
    real(r8) wnd_mdp       ! [m s-1] Surface layer mean wind speed
    !// ------------------------------------------------------------------
    ! Main Code
    
    ! Construct land, ocean, and sea-ice flags
    if (oro_is_lnd(oro)) then
       flg_lnd = .true.
    else 
       flg_lnd = .false.
    endif                  ! endif land
    if (oro_is_ocn(oro)) then
       flg_ocn = .true.
    else 
       flg_ocn = .false.
    endif                  ! endif ocean
    if (oro_is_ice(oro)) then
       flg_ice = .true.
    else 
       flg_ice = .false.
    endif                  ! endif ice
#ifdef DST_DBG
    ! Sanity check
    if (.not.flg_lnd.and..not.flg_ocn.and..not.flg_ice) stop 'Invalid surface type in blm_glb()'
#endif /* not DST_DBG */

    
    !// Midlayer wind speed
    flg_true = .true. ! [flg] True flag
    wnd_mdp =  &      ! [m s-1] Surface layer mean wind speed
         sqrt(wnd_znl_mdp * wnd_znl_mdp +  &
            wnd_mrd_mdp * wnd_mrd_mdp) 

    
    !// Fraction of surface covered by snow (checked ok for CTM3)
    call snw_frc_get( &
         snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
         snw_frc)             ! O [frc] Fraction of surface covered by snow

    
    !// Zero plane displacement (checked ok for CTM3)
    call zpd_get( &
         hgt_zpd,             & ! O [m] Zero plane displacement
         oro,                 & ! I [frc] Orography
         sfc_typ,             & ! I [idx] LSM surface type (0..28)
         snw_frc)               ! I [frc] Fraction of surface covered by snow
    

    !alfgr: For high vertical resolutions, zpd can be higher than hgt_mdp (not physical)
    !Dirty fix: (limit displacement height to be 70% of midpoint height)
    hgt_zpd = min(0.7*hgt_mdp, hgt_zpd)


    ! Roughness length (checked ok for CTM3)
    call rgh_mmn_get( &
         oro,                 & ! I [frc] Orography
         rgh_mmn,             & ! O [m] Roughness length momentum
         sfc_typ,             & ! I [idx] LSM surface type (0..28)
         snw_frc,             & ! I [frc] Fraction of surface covered by snow
         wnd_mdp)             ! I [m s-1] Surface layer mean wind speed

    !// Find surface properties over land (checked ok for CTM3)
    call blm_lnd( &
         dns_mdp,             & ! I [kg m-3] Midlayer density
         flg_lnd,             & ! I [flg] Land flag
         hgt_mdp,             & ! I [m] Midlayer height above surface
         prs_mdp,             & ! I [Pa] Pressure
         q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
         rgh_mmn,             & ! I [m] Roughness length momentum
         tpt_mdp,             & ! I [K] Midlayer temperature
         tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
         tpt_sfc,             & ! I [K] Surface temperature
         wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
         wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
         mno_lng,             & ! O [m] Monin-Obukhov length
         wnd_frc)             ! O [m s-1] Surface friction velocity
    
    !// Find surface properties over ocean (checked ok for CTM3)
    call blm_ocn( &
         dns_mdp,             & ! I [kg m-3] Midlayer density
         flg_ocn,             & ! I [flg] Ocean flag
         hgt_mdp,             & ! I [m] Midlayer height above surface
         prs_mdp,             & ! I [Pa] Pressure
         q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
         tpt_mdp,             & ! I [K] Midlayer temperature
         tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
         tpt_sfc,             & ! I [K] Surface temperature
         wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
         wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
         mno_lng,             & ! O [m] Monin-Obukhov length
         rgh_mmn,             & ! O [m] Roughness length momentum
         wnd_frc)             ! O [m s-1] Surface friction velocity
    
    !// Find surface properties over ice (checked ok for CTM3)
    call blm_ice( &
         dns_mdp,             & ! I [kg m-3] Midlayer density
         flg_ice,             & ! I [flg] Sea ice flag
         hgt_mdp,             & ! I [m] Midlayer height above surface
         prs_mdp,             & ! I [Pa] Pressure
         q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
         rgh_mmn,             & ! I [m] Roughness length momentum
         tpt_mdp,             & ! I [K] Midlayer temperature
         tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
         tpt_sfc,             & ! I [K] Surface temperature
         wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
         wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
         mno_lng,             & ! O [m] Monin-Obukhov length
         wnd_frc)             ! O [m s-1] Surface friction velocity
    
    ! Interpolate midlayer wind speed to 10 m (checked ok for CTM3)
    call wnd_rfr_get( &
         flg_true,            & ! I [flg] True flag
         hgt_mdp,             & ! I [m] Midpoint height above surface
         hgt_rfr,             & ! I [m] Reference height for deposition processes
         hgt_zpd,             & ! I [m] Zero plane displacement
         mno_lng,             & ! I [m] Monin-Obukhov length
         wnd_frc,             & ! I [m s-1] Surface friction velocity
         wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
         wnd_min_dps,         & ! I [m s-1] Minimum windspeed used for deposition 
         wnd_rfr)             ! O [m s-1] Wind speed at reference height
    

    !// ------------------------------------------------------------------
  end subroutine blm_glb                       ! end blm_glb()
  !// ------------------------------------------------------------------
  
end module blmutl

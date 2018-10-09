! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstblm.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: dstblm.F90 contains boundary layer meteorology driver and utilities 
! required for sub-gridscale, non-vegetated land surfaces
! Usage:
! use dstblm ! [mdl] Boundary layer meteorology for non-vegetated land surfaces

module dstblm ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  save ! [stt] Changes to common variables are sticky
  ! private ! [stt] Symbols are private unless individually qualified as public
  private::dsvpdt_H2O_ice_PrK78_fst_scl ! [fnc] Derivative of saturation vapor pressure over planar ice water
  private::dsvpdt_H2O_lqd_PrK78_fst_scl ! [fnc] Derivative of saturation vapor pressure over planar liquid water
  public::svp_H2O_ice_PrK78_fst_scl ! [fnc] Saturation vapor pressure over planar ice water
  public::svp_H2O_lqd_PrK78_fst_scl ! [fnc] Saturation vapor pressure over planar liquid water
  public::tpt_bnd_cls_get ! [fnc] Given T [K], Return -50 < T [C] < 50 C
  
  ! Physical constants for boundary layer meteorology
  ! These variables are initialized in dstcmnini:dst_blm_cmn_ini()
  real(r8) cp_vpr_rcp_cp_dry_m1 ! 0.83745 [frc] Constant for moist specific heat IrG81 p. 77
  real(r8) cst_Stefan_Boltzmann ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
  real(r8) dns_H2O_lqd_std      ! [kg m-3] Density of liquid water
  real(r8) ltn_heat_fsn_H2O_std ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
  real(r8) ltn_heat_sbl_H2O_std ! [J kg-1] Latent heat of sublimation of H2O, standard
  real(r8) ltn_heat_vpr_H2O_std ! [J kg-1] Latent heat of vaporization of H2O, standard
  real(r8) one_mns_eps_H2O      ! (0.378) [frc] Constant for saturation specific humidity
  real(r8) spc_heat_H2O_ice_vlm ! [J m-3 K-1] Volumetric specific heat capacity of ice water
  real(r8) spc_heat_H2O_lqd_vlm ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
  real(r8) spc_heat_dry_air     ! (1005.0) [J kg-1 K-1] IrG81 p. 25
  real(r8) tpt_frz_pnt          ! [K] Kelvin--Celsius scale offset
  
contains
  
  real(r8) function dsvpdt_H2O_lqd_PrK78_fst_scl( & ! [fnc] Derivative of saturation vapor pressure over planar liquid water
       tpt_cls) ! I [C] Temperature in celsius
    ! Purpose: Derivative of saturation vapor pressure over planar liquid water
    ! Coefficients for d(SVP)/dT over planar liquid water
    real(r8),parameter::c0=4.438099984e-01
    real(r8),parameter::c1=2.857002636e-02
    real(r8),parameter::c2=7.938054040e-04
    real(r8),parameter::c3=1.215215065e-05
    real(r8),parameter::c4=1.036561403e-07
    real(r8),parameter::c5=3.532421810e-10
    real(r8),parameter::c6=-7.090244804e-13
    real(r8),intent(in)::tpt_cls ! [C] Temperature in celsius
    dsvpdt_H2O_lqd_PrK78_fst_scl= & ! [Pa] 
         100.0_r8*(c0+tpt_cls*(c1+tpt_cls*(c2+tpt_cls*(c3+tpt_cls*(c4+tpt_cls*(c5+tpt_cls*c6))))))
  end function dsvpdt_H2O_lqd_PrK78_fst_scl

  real(r8) function dsvpdt_H2O_ice_PrK78_fst_scl( & ! [fnc] Derivative of saturation vapor pressure over planar ice water
       tpt_cls) ! I [C] Temperature in celsius
    ! Purpose: Derivative of saturation vapor pressure over planar ice water
    ! Coefficients for d(SVP)/dT over planar ice water
    real(r8),parameter::d0=5.030305237e-01
    real(r8),parameter::d1=3.773255020e-02
    real(r8),parameter::d2=1.267995369e-03
    real(r8),parameter::d3=2.477563108e-05
    real(r8),parameter::d4=3.005693132e-07
    real(r8),parameter::d5=2.158542548e-09
    real(r8),parameter::d6=7.131097725e-12
    real(r8),intent(in)::tpt_cls ! [C] Temperature in celsius
    dsvpdt_H2O_ice_PrK78_fst_scl= & ! [Pa] 
         100.0_r8*(d0+tpt_cls*(d1+tpt_cls*(d2+tpt_cls*(d3+tpt_cls*(d4+tpt_cls*(d5+tpt_cls*d6))))))
  end function dsvpdt_H2O_ice_PrK78_fst_scl

  real(r8) function svp_H2O_lqd_PrK78_fst_scl( & ! [fnc] Saturation vapor pressure over planar liquid water
       tpt_cls) ! I [C] Temperature in celsius
    ! Purpose: Saturation vapor pressure over planar liquid water
    ! Coefficients for SVP over planar liquid water
    ! Inline statement functions for thermodynamics from Lowe and Ficke (1974) as reported in PrK78 p. 625. Range of validity is -50 C < T < 50 C.
    real(r8),parameter::a0=6.107799961
    real(r8),parameter::a1=4.436518521e-01
    real(r8),parameter::a2=1.428945805e-02
    real(r8),parameter::a3=2.650648471e-04
    real(r8),parameter::a4=3.031240396e-06
    real(r8),parameter::a5=2.034080948e-08
    real(r8),parameter::a6=6.136820929e-11
    real(r8),intent(in)::tpt_cls ! [C] Temperature in celsius
    svp_H2O_lqd_PrK78_fst_scl= & ! [Pa] 
         100.0_r8*(a0+tpt_cls*(a1+tpt_cls*(a2+tpt_cls*(a3+tpt_cls*(a4+tpt_cls*(a5+tpt_cls*a6))))))
  end function svp_H2O_lqd_PrK78_fst_scl

  real(r8) function svp_H2O_ice_PrK78_fst_scl( & ! [fnc] Saturation vapor pressure over planar ice water
       tpt_cls) ! I [C] Temperature in celsius
    ! Purpose: Saturation vapor pressure over planar ice water
    ! Coefficients for SVP over planar ice water
    real(r8),parameter::b0=6.109177956
    real(r8),parameter::b1=5.034698970e-01
    real(r8),parameter::b2=1.886013408e-02
    real(r8),parameter::b3=4.176223716e-04
    real(r8),parameter::b4=5.824720280e-06
    real(r8),parameter::b5=4.838803174e-08
    real(r8),parameter::b6=1.838826904e-10
    real(r8),intent(in)::tpt_cls ! [C] Temperature in celsius
    svp_H2O_ice_PrK78_fst_scl= & ! [Pa] 
         100.0_r8*(b0+tpt_cls*(b1+tpt_cls*(b2+tpt_cls*(b3+tpt_cls*(b4+tpt_cls*(b5+tpt_cls*b6))))))
  end function svp_H2O_ice_PrK78_fst_scl

  real(r8) function tpt_bnd_cls_get( & ! [fnc] Given T [K], Return -50 < T [C] < 50 C
       tpt) ! I [K] Temperature
    ! Purpose: Given temperature T in [K], return bounded temperature in [C], 
    ! i.e., -50 < T [C] < 50 C
    real(r8),intent(in)::tpt ! I [K] Temperature
    tpt_bnd_cls_get=min(50.0_r8,max(-50.0_r8,(tpt-tpt_frz_pnt))) ! [C] Temperature bounded celsius
  end function tpt_bnd_cls_get
    
  subroutine dst_blm_cmn_ini()
    ! Purpose: Initialize boundary layer meteorology module
    ! dst_blm_cmn_ini() is called by dst_msc_cmn_ini()
    ! dst_msc_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun() 
    use dstcst,only:eps_H2O,gas_cst_dry_air,gas_cst_unv,mmw_H2O ! [mdl] Physical constants for dust routines
    implicit none
    ! Parameters for boundary layer meteorology
    real(r8),parameter::cst_Stefan_Boltzmannx=5.67032e-8 ! [W m-2 K-4] Stefan-Boltzmann constant
    real(r8),parameter::dns_H2O_lqd_stdx=1000.0 ! [kg m-3] Density of liquid water
    real(r8),parameter::ltn_heat_fsn_H2O_stdx=0.3336e06 ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard CCM:lsm/phyconi.F 
    real(r8),parameter::ltn_heat_sbl_H2O_stdx=2.8440e06 ! [J kg-1] Latent heat of sublimation of H2O, standard CCM:lsm/phyconi.F
    real(r8),parameter::ltn_heat_vpr_H2O_stdx=2.5104e06 ! [J kg-1] Latent heat of vaporization of H2O, standard CCM:lsm/phyconi.F
    real(r8),parameter::spc_heat_H2O_ice_vlmx=2.094e06 ! [J m-3 K-1] Volumetric specific heat capacity of ice water CCM:lsm/phyconi.F
    real(r8),parameter::spc_heat_H2O_lqd_vlmx=4.188e06 ! [J m-3 K-1] Volumetric specific heat capacity of liquid water CCM:lsm/phyconi.F
    real(r8),parameter::tpt_frz_pntx=273.15 ! [K] Kelvin--Celsius scale offset Bol80
    ! Local
    real(r8) gas_cst_H2O          ! (461.65) [J kg-1 K-1] Gas constant of H2O
    real(r8) spc_heat_H2O_vpr     ! (1850.0) [J kg-1 K-1] IrG81 pp. 77, 245
    
    ! Parameters for boundary layer meteorology
    cst_Stefan_Boltzmann=cst_Stefan_Boltzmannx ! (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant
    dns_H2O_lqd_std=dns_H2O_lqd_stdx ! [kg m-3] Density of liquid water
    ltn_heat_fsn_H2O_std=ltn_heat_fsn_H2O_stdx ! [J kg-1] Latent heat of fusion of H2O at 0 C, standard
    ltn_heat_sbl_H2O_std=ltn_heat_sbl_H2O_stdx ! [J kg-1] Latent heat of sublimation of H2O, standard
    ltn_heat_vpr_H2O_std=ltn_heat_vpr_H2O_stdx ! [J kg-1] Latent heat of vaporization of H2O, standard
    one_mns_eps_H2O=1.0_r8-eps_H2O   ! (0.378) [frc] Constant for saturation specific humidity
    spc_heat_H2O_ice_vlm=spc_heat_H2O_ice_vlmx ! [J m-3 K-1] Volumetric specific heat capacity of ice water
    spc_heat_H2O_lqd_vlm=spc_heat_H2O_lqd_vlmx ! [J m-3 K-1] Volumetric specific heat capacity of liquid water
    tpt_frz_pnt=tpt_frz_pntx  ! [K] Kelvin--Celsius scale offset
    
    ! Derived 
    spc_heat_dry_air=7.0_r8*gas_cst_dry_air/2.0_r8 ! (1005.0) [J kg-1 K-1]
    gas_cst_H2O=gas_cst_unv/mmw_H2O ! (461.65) [J kg-1 K-1] Gas constant of H2O
    spc_heat_H2O_vpr=4.0_r8*gas_cst_H2O ! (1850.0) [J kg-1 K-1] Specific heat of H2O vapor IrG81 pp. 77, 245
    cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0_r8 ! (0.83745) [frc] Constant for moist specific heat IrG81 p. 77
    
    ! Sanity checks
    if (abs(273.15_r8-tpt_frz_pnt)/273.15_r8 > 1.0e-4_r8) stop 'dst_blm_cmn_ini(): tpt_frz_pnt error'
    if (abs(1004.697_r8-spc_heat_dry_air)/1004.697_r8 > 2.0e-4_r8) then
       print*, 'dst_blm_cmn_ini(): spc_heat_dry_air error',spc_heat_dry_air,gas_cst_dry_air
       print*,(abs(1004.697_r8-spc_heat_dry_air)/1004.697_r8),1004.697_r8
       print*,r8
       stop
    endif
    if (abs(461.52_r8-gas_cst_H2O)/461.52_r8 > 1.0e-4_r8) stop 'dst_blm_cmn_ini(): gas_cst_H2O error'
    if (abs(1846.08_r8-spc_heat_H2O_vpr)/1846.08_r8 > 1.0e-4_r8) stop 'dst_blm_cmn_ini(): spc_heat_H2O_vpr error'
    if (abs(0.83745_r8-cp_vpr_rcp_cp_dry_m1)/0.83745_r8 > 1.0e-3_r8) then
       write (6,'(a,f11.6)') 'cp_vpr_rcp_cp_dry_m1 = ',cp_vpr_rcp_cp_dry_m1
       stop 'dst_blm_cmn_ini(): cp_vpr_rcp_cp_dry_m1 error'
    endif
    
    return
  end subroutine dst_blm_cmn_ini
  
  

  !// ------------------------------------------------------------------
  subroutine hyd_prp_get( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       mss_frc_cly,         & ! I [frc] Mass fraction clay 
       mss_frc_snd,         & ! I [frc] Mass fraction sand
       vwc_dry,             & ! O [m3 m-3] Dry volumetric water content (no E-T)
       vwc_opt,             & ! O [m3 m-3] E-T optimal volumetric water content
       vwc_sat)             ! O [m3 m-3] Saturated volumetric water content
    !// ------------------------------------------------------------------
    ! Purpose: Determine hydrologic properties from soil texture
    ! NB: All I/O for this routine is time-invariant
    ! Thus, the hydrologic properties coulde be computed once at initialization
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
    ! Commons
    ! Input
    logical,intent(in)::flg_mbl    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::mss_frc_cly   ! I [frc] Mass fraction clay
    real(r8),intent(in)::mss_frc_snd   ! I [frc] Mass fraction sand
    ! Output
    real(r8),intent(out)::vwc_dry       ! O [m3 m-3] Dry volumetric water content (no E-T)
    real(r8),intent(out)::vwc_opt       ! O [m3 m-3] E-T optimal volumetric water content 
    real(r8),intent(out)::vwc_sat       ! O [m3 m-3] Saturated volumetric water content (sand-dependent)
    ! Local
    integer idx_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    real(r8) smp_xpn_b     ! [frc] Exponent "b" for smp (clay-dependent)
    real(r8) smp_sat       ! [mm H2O] Saturated soil matric potential (sand-dependent)
    !// ------------------------------------------------------------------
    ! Main Code
    ! Initialize outputs
    vwc_dry = 0.0_r8   ! [m3 m-3] Dry volumetric water content (no E-T)
    vwc_opt = 0.0_r8   ! [m3 m-3] E-T optimal volumetric water content 
    vwc_sat = 0.0_r8   ! [m3 m-3] Saturated volumetric water content (sand-dependent)

    
    ! Time-invariant soil hydraulic properties
    ! See Bon96 p. 98, implemented in CCM:lsm/lsmtci()
    if (flg_mbl) then
       smp_xpn_b = & ! [frc] Exponent "b" for smp (clay-dependent)
            2.91_r8+0.159_r8*mss_frc_cly*100.0_r8
       ! NB: Adopt convention that matric potential is positive definite
       smp_sat =   & ! [mm H2O] Saturated soil matric potential (sand-dependent)
            10.0_r8*(10.0_r8**(1.88_r8-0.0131_r8*mss_frc_snd*100.0_r8))
       vwc_sat =   & ! [m3 m-3] Saturated volumetric water content (sand-dependent)
            0.489_r8-0.00126_r8*mss_frc_snd*100.0_r8
       vwc_dry =   & ! [m3 m-3] Dry volumetric water content (no E-T)
            vwc_sat* &
            (316230.0_r8/smp_sat)**(-1.0_r8/smp_xpn_b)
       vwc_opt=   & ! [m3 m-3] E-T optimal volumetric water content
            vwc_sat* &
            (158490.0_r8/smp_sat)**(-1.0_r8/smp_xpn_b)
    endif                  ! endif flg_mbl

    !// ------------------------------------------------------------------
  end subroutine hyd_prp_get                       ! end hyd_prp_get()
  !// ------------------------------------------------------------------


  !// ------------------------------------------------------------------
  subroutine cnd_trm_soi_get( &
       cnd_trm_soi,         & ! O [W m-1 K-1] Soil thermal conductivity
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       lvl_dlt,             & ! O [m] Soil layer thickness
       mss_frc_cly,         & ! I [frc] Mass fraction clay 
       mss_frc_snd,         & ! I [frc] Mass fraction sand
       tpt_soi,             & ! I [K] Soil temperature
       vwc_dry,             & ! I [m3 m-3] Dry volumetric water content (no E-T)
       vwc_opt,             & ! I [m3 m-3] E-T optimal volumetric water content
       vwc_sat,             & ! I [m3 m-3] Saturated volumetric water content
       vwc_sfc)             ! I [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------
    ! Purpose: Thermal properties of soil
    ! NB: Currently this routine is optimized for ground without snow-cover
    ! Although snow thickness is read in, it is not currently used
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
    real(r8),parameter::cnd_trm_H2O_ice=2.2 ! [W m-1 K-1] Thermal conductivity of ice water CCM:lsm/phyconi
    real(r8),parameter::cnd_trm_H2O_lqd=0.6 ! [W m-1 K-1] Thermal conductivity of liquid water CCM:lsm/phyconi
    real(r8),parameter::cnd_trm_snw=0.34 ! [W m-1 K-1] Thermal conductivity of snow CCM:lsm/snoconi, Bon96 p. 77
    real(r8),parameter::lvl_dlt_sfc=0.1 ! [m] Soil layer thickness, top layer
    real(r8),parameter::tpt_dlt=0.5    ! [K] Temperature range of mixed phase soil
    ! Input
    logical,intent(in)::flg_mbl    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::mss_frc_cly   ! I [frc] Mass fraction clay 
    real(r8),intent(in)::mss_frc_snd   ! I [frc] Mass fraction sand
    ! real(r8),intent(in)::snw_hgt       ! I [m] Geometric bulk thickness of snow
    real(r8),intent(in)::tpt_soi       ! I [K] Soil temperature
    real(r8),intent(in)::vwc_dry       ! I [m3 m-3] Dry volumetric water content (no E-T)
    real(r8),intent(in)::vwc_opt       ! I [m3 m-3] E-T optimal volumetric water content
    real(r8),intent(in)::vwc_sat       ! I [m3 m-3] Saturated volumetric water content
    real(r8),intent(in)::vwc_sfc       ! I [m3 m-3] Volumetric water content
    ! Output
    real(r8),intent(out)::cnd_trm_soi   ! O [W m-1 K-1] Soil thermal conductivity
    real(r8),intent(out)::lvl_dlt       ! O [m] Soil layer thickness
    ! Local
    real(r8) cnd_trm_soi_dry ! [W m-1 K-1] Thermal conductivity of dry soil
    real(r8) cnd_trm_soi_frz ! [W m-1 K-1] Soil thermal conductivity, frozen
    real(r8) cnd_trm_soi_sld ! [W m-1 K-1] Thermal conductivity of soil solids
    real(r8) cnd_trm_soi_wrm ! [W m-1 K-1] Soil thermal conductivity, unfrozen
    real(r8) ltn_heat_fsn_vlm ! [J m-3] Volumetric latent heat of fusion
    real(r8) snw_hgt_bnd          ! [m] Bounded geometric bulk thickness of snow
    !// ------------------------------------------------------------------
    ! Main Code
    ! Initialize arrays
    lvl_dlt = lvl_dlt_sfc    ! [m] Soil layer thickness
    cnd_trm_soi = 0.0_r8        ! [W m-1 K-1] Soil thermal conductivity
    

    if (flg_mbl) then
       ltn_heat_fsn_vlm = vwc_sfc*ltn_heat_fsn_H2O_std*dns_H2O_lqd_std ! [J m-3] Volumetric latent heat of fusion
       cnd_trm_soi_sld = & ! [W m-1 K-1] Thermal conductivity of soil solids CCM:lsm/lsmtci() Bon96 p. 77
            (8.80_r8*mss_frc_snd+2.92_r8*mss_frc_cly) &
            /(mss_frc_snd+mss_frc_cly)
       cnd_trm_soi_dry = 0.15_r8 ! [W m-1 K-1] Thermal conductivity of dry soil CCM:lsm/lsmtci() Bon96 p. 77
       cnd_trm_soi_wrm = & ! [W m-1 K-1] Soil thermal conductivity, unfrozen
            +cnd_trm_soi_dry &
            +(cnd_trm_soi_sld**(1.0_r8-vwc_sat) &
            *(cnd_trm_H2O_lqd**vwc_sfc)-cnd_trm_soi_dry) &
            *vwc_sfc/vwc_sat
       cnd_trm_soi_frz = & ! [W m-1 K-1] Soil thermal conductivity, frozen
            +cnd_trm_soi_dry &
            +(cnd_trm_soi_sld**(1.0_r8-vwc_sat) &
            *(cnd_trm_H2O_ice**vwc_sfc)-cnd_trm_soi_dry) &
            *vwc_sfc/vwc_sat
       if (tpt_soi < tpt_frz_pnt-tpt_dlt) then
          cnd_trm_soi = cnd_trm_soi_frz ! [W m-1 K-1] Soil thermal conductivity
       endif               ! endif
       if (tpt_soi >= tpt_frz_pnt-tpt_dlt.and.tpt_soi <= tpt_frz_pnt+tpt_dlt) then
          cnd_trm_soi = & ! [W m-1 K-1] Soil thermal conductivity
               +cnd_trm_soi_frz &
               +(cnd_trm_soi_frz-cnd_trm_soi_wrm) &
               *(tpt_soi-tpt_frz_pnt+tpt_dlt) &
               /(2.0_r8*tpt_dlt)
       endif               ! endif
       if (tpt_soi > tpt_frz_pnt+tpt_dlt) then
          cnd_trm_soi = cnd_trm_soi_wrm ! [W m-1 K-1] Soil thermal conductivity
       endif               ! endif
          
       ! ! Blend snow into first soil layer
       ! Snow is not allowed to cover dust mobilization regions
       ! snw_hgt_bnd=min(snw_hgt(lon_idx),1.0_r8) ! [m] Bounded geometric bulk thickness of snow
       ! lvl_dlt_snw(lon_idx)=lvl_dlt(lon_idx)+snw_hgt_bnd ! O [m] Soil layer thickness including snow Bon96 p. 77
       ! cnd_trm_soi(lon_idx)= & ! [W m-1 K-1] Soil thermal conductivity Bon96 p. 77
       !         cnd_trm_snw*cnd_trm_soi(lon_idx)*lvl_dlt_snw(lon_idx) &
       !         /(cnd_trm_snw*lvl_dlt(lon_idx)+cnd_trm_soi(lon_idx)*snw_hgt_bnd)
    endif                  ! endif flg_mbl

    !// ------------------------------------------------------------------
  end subroutine cnd_trm_soi_get                       ! end cnd_trm_soi_get()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine trn_fsh_vpr_soi_atm_get( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       tpt_soi,             & ! I [K] Soil temperature
       tpt_soi_frz,         & ! I [K] Temperature of frozen soil
       trn_fsh_vpr_soi_atm, & ! O [frc] Transfer efficiency of vapor from soil to atmosphere
       vwc_dry,             & ! I [m3 m-3] Dry volumetric water content (no E-T)
       vwc_opt,             & ! I [m3 m-3] E-T optimal volumetric water content
       vwc_sfc)             ! I [m3 m-3] Volumetric water content
    !// ------------------------------------------------------------------
    ! Purpose: Compute factor describing effects of soil texture and moisture on vapor transfer between soil and atmosphere
    ! Taken from Bon96 p. 59, CCM:lsm/surphys
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
    real(r8),parameter::trn_fsh_vpr_soi_atm_frz=0.01 ! [frc] Transfer efficiency of vapor from frozen soil to atmosphere CCM:lsm/surphy()
    ! Input
    logical,intent(in)::flg_mbl    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::tpt_soi       ! I [K] Soil temperature
    real(r8),intent(in)::tpt_soi_frz          ! I [K] Temperature of frozen soil
    real(r8),intent(in)::vwc_dry       ! I [m3 m-3] Dry volumetric water content (no E-T)
    real(r8),intent(in)::vwc_opt       ! I [m3 m-3] E-T optimal volumetric water content
    real(r8),intent(in)::vwc_sfc       ! I [m3 m-3] Volumetric water content
    ! Output
    real(r8),intent(out)::trn_fsh_vpr_soi_atm ! O [frc] Transfer efficiency of vapor from soil to atmosphere
    !// ------------------------------------------------------------------
    
    ! The trn_fsh_vpr_soi_atm efficiency factor attempts to tie soil texture and 
    ! moisture properties to the vapor conductance of the soil-atmosphere system.
    ! When the soil temperature is sub-freezing, the conductance describes the 
    ! resistance to vapor sublimation (or deposition) and transport through the 
    ! open soil pores to the atmosphere.
    ! For warm soils, vapor transfer is most efficient at the optimal VWC for E-T
    ! Thus when vwc_sfc = vwc_opt, soil vapor transfer is perfectly efficient 
    ! (trn_fsh_vpr_soi_atm = 1.0) so the soil does not contribute any resistance  
    ! to the surface vapor transfer.
    ! When vwc_sfc > vwc_opt, the soil has an excess of moisture and, again,
    ! vapor transfer is not limited by soil characteristics.
    ! In fact, according to Bon96 p. 98, vwc_dry is only slightly smaller than
    ! vwc_opt, so trn_fsh_vpr_soi_atm is usually either 0 or 1 and intermediate
    ! efficiencies occur over only a relatively small range of VWC.
    ! When vwc_sfc < vwc_dry, the soil matrix is subsaturated and acts as a 
    ! one-way sink for vapor through osmotic and capillary potentials.
    ! In this case trn_fsh_vpr_soi_atm = 0, which would cause the surface resistance
    ! rss_vpr_sfc to blow up, but this is guarded against and rss_sfc_vpr 
    ! is set to ~1.0e6*rss_aer_vpr instead. 
    ! Note that this formulation does not seem to allow vapor transfer from
    ! the atmosphere to the soil when vwc_sfc < vwc_dry, even when 
    ! e_atm > esat(Tg).
    ! Air at the apparent sink for moisture is has vapor pressure e_sfc
    ! e_atm = Vapor pressure of ambient air at z = hgt_mdp
    ! e_sfc = Vapor pressure at apparent sink for moisture at z = zpd + rgh_vpr 
    ! e_gnd = Vapor pressure at air/ground interface temperature 
    ! Air at the soil interface is assumed saturated, i.e., e_gnd = esat(Tg)
    
    ! Main code
    if (flg_mbl) then
       if (tpt_soi > tpt_soi_frz) then
          trn_fsh_vpr_soi_atm = & ! [frc] Transfer efficiency of vapor from soil to atmosphere
                  min(max(vwc_sfc-vwc_dry,0.0_r8)/(vwc_opt-vwc_dry),1.0_r8) ! CCM:lsm/surphys Bon96 p. 59
       else
          trn_fsh_vpr_soi_atm = trn_fsh_vpr_soi_atm_frz ! [frc] Bon96 p. 59
       endif
    endif                  ! endif flg_mbl

    !// ------------------------------------------------------------------
  end subroutine trn_fsh_vpr_soi_atm_get                       ! end trn_fsh_vpr_soi_atm_get()
  !// ------------------------------------------------------------------


  !// ------------------------------------------------------------------
  subroutine blm_mbl( &
       cnd_trm_soi,         & ! I [W m-1 K-1] Soil thermal conductivity
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_LW_dwn_sfc,      & ! I [W m-2] Longwave downwelling flux at surface
       flx_SW_abs_sfc,      & ! I [W m-2] Solar flux absorbed by ground
       hgt_mdp,             & ! I [m] Midlayer height above surface
       hgt_zpd,             & ! I [m] Zero plane displacement
       lvl_dlt,             & ! I [m] Soil layer thickness
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
       rgh_mmn,             & ! I [m] Roughness length momentum
       tpt_mdp,             & ! I [K] Midlayer temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
       tpt_soi,             & ! I [K] Soil temperature
       trn_fsh_vpr_soi_atm, & ! I [frc] Transfer efficiency of vapor from soil to atmosphere
       !wnd_mrd_mdp,         & ! I [m s-1] Surface layer meridional wind speed
       !wnd_znl_mdp,         & ! I [m s-1] Surface layer zonal wind speed
       wnd_mdp,             & !I [m s-1] Surface layer wind speed
       mno_lng,             & ! O [m] Monin-Obukhov length
       tpt_gnd,             & ! I/O [K] Ground temperature
       wnd_frc)             ! O [m s-1] Surface friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given meteorology at reference level (i.e., GCM layer midpoint)
    ! compute the boundary layer exchange properties
    ! Routine is optimized for dust source regions: dry, bare, uncovered land
    ! Theory and algorithms: Bonan (1996) CCM:lsm/surtem()
    ! Notes on algorithm:
    ! Suffix mdp quantity evaluated at height hgt_mdp
    ! Suffix gnd quantity evaluated at ground
    ! Relationships between various temperatures: Bon96 p. 55 and Fig. 16 p. 57:
    ! T1 = Skin temperature = Soil temperature of 1st layer (top 10 cm)
    ! Tg = "Ground" temperature = Air temperature at z = rgh_heat
    ! Ts = "Surface" temperature = Air temperature at z=zpd+rgh_heat
    ! Ta = "Aerodynamic" temperature = Air temperature at z=zpd+rgh_mmn
    ! Te = Emitting temperature = (FLWup/sigma)**0.25
    ! Tgcm = Ambient temperature = Air temperature at z=hgt_mdp
    ! For bare ground, Ts = Tg = Air temperature at z=rgh_heat (Bon96 p. 55)
    ! For bare ground, Tg = potential temperture defined relative to local surface pressure
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
    integer,parameter::itr_max=10     ! [nbr] Maximum number of iterations for tpt_gnd loop
    integer,parameter::sgn_chg_ctr_max=4 ! [nbr] Maximum number of sign changes in stability parameter
    real(r8),parameter::eps_dbz=1.0e-6 ! [frc] Prevents division by zero
    real(r8),parameter::sfc_ems_gnd=0.96 ! [frc] Surface emissivity of bare ground CCM:lsm/snoconi()
    real(r8),parameter::wnd_min_mbl=1.0 ! [m s-1] Minimum windspeed used for mobilization
    ! Input
    logical,intent(in)::flg_mbl    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::cnd_trm_soi   ! I [W m-1 K-1] Soil thermal conductivity
    real(r8),intent(in)::dns_mdp       ! I [kg m-3] Midlayer density
    real(r8),intent(in)::flx_LW_dwn_sfc ! I [W m-2] Longwave downwelling flux at surface
    real(r8),intent(in)::flx_SW_abs_sfc ! I [W m-2] Solar flux absorbed by ground
    real(r8),intent(in)::hgt_mdp       ! I [m] Midlayer height above surface
    real(r8),intent(in)::hgt_zpd              ! I [m] Zero plane displacement
    real(r8),intent(in)::lvl_dlt       ! I [m] Soil layer thickness
    real(r8),intent(in)::prs_mdp       ! I [Pa] Pressure
    real(r8),intent(in)::q_H2O_vpr     ! I [kg kg-1] Specific humidity
    real(r8),intent(in)::rgh_mmn              ! I [m] Roughness length momentum
    real(r8),intent(in)::trn_fsh_vpr_soi_atm ! I [frc] Transfer efficiency of vapor from soil to atmosphere
    real(r8),intent(in)::tpt_mdp       ! I [K] Midlayer temperature
    real(r8),intent(in)::tpt_ptn_mdp   ! I [K] Midlayer local potential temperature (relative to surface pressure, not to 1000 mb)
    real(r8),intent(in)::tpt_soi       ! I [K] Soil temperature
    !real(r8),intent(in)::wnd_mrd_mdp   ! I [m s-1] Surface layer meridional wind speed
    !real(r8),intent(in)::wnd_znl_mdp   ! I [m s-1] Surface layer zonal wind speed
    real(r8),intent(in):: wnd_mdp       !I [m s-1] Surface layer wind speed 
    ! Input/Output
    real(r8),intent(inout)::tpt_gnd       ! I/O [K] Ground temperature
    ! Output
    real(r8),intent(out)::mno_lng       ! O [m] Monin-Obukhov length
    real(r8),intent(out)::wnd_frc       ! O [m s-1] Surface friction velocity
    ! Local
    integer itr_idx           ! [idx] Counting index
    integer sgn_chg_ctr ! [nbr] Number of sign changes in stability parameter
    real(r8) cff_xch_mmn   ! [frc] Exchange coefficient for momentum transfer
    real(r8) cnd_heat_sfc_mdp     ! [m s-1] Sensible heat conductance surface air to midlayer air
    real(r8) cnd_vpr_gnd_sfc      ! [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
    real(r8) cnd_vpr_sfc_mdp      ! [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57
    real(r8) cnd_vpr_ttl          ! [m s-1] Sum of conductances
    real(r8) cst_psych     ! [Pa K-1] Psychrometric constant
    real(r8) dsvpdt_H2O_gnd ! [Pa K-1] Derivative of saturation vapor pressure over planar condensed water, ground
    real(r8) flx_LW_net    ! [W m-2] Net longwave flux to atmosphere
    real(r8) flx_LW_net_cff_a ! [W m-2] a in FLWnet = a + b*T^4 Bon96 p. 45
    real(r8) flx_LW_net_cff_b ! [W m-2 K-4] b in FLWnet = a + b*T^4 Bon96 p. 45
    real(r8) flx_LW_up_sfc ! [W m-2] Longwave upwelling flux at surface
    real(r8) flx_ltn_evp   ! [W m-2] Evaporation flux to atmosphere
    real(r8) flx_ltn_evp_cff_a ! [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
    real(r8) flx_ltn_evp_cff_b ! [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
    real(r8) flx_ltn_fct          ! Factor in vapor flux calculations
    real(r8) flx_sns_atm   ! [W m-2] Sensible heat flux to atmosphere
    real(r8) flx_sns_atm_cff_a ! [W m-2] a in SH = a + b*Ts
    real(r8) flx_sns_atm_cff_b ! [W m-2 K-1] b in SH = a + b*Ts
    real(r8) flx_sns_atm_fct      ! Factor in sensible heat computation
    real(r8) flx_sns_atm_tmp      ! [W m-2] Temporary sensible heat flux
    real(r8) flx_sns_atm_vrt_tmp  ! [W m-2] Temporary virtual sensible heat flux Bon96 p. 49
    real(r8) flx_sns_gnd   ! [W m-2] Sensible heat flux to soil
    real(r8) flx_sns_gnd_cff_a ! [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
    real(r8) flx_sns_gnd_cff_b ! [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
    real(r8) flx_vpr_tmp          ! [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
    real(r8) ltn_heat_trn  ! [J kg-1] Latent heat of sublimation or evaporation
    real(r8) mno_dnm              ! Denominator of Monin-Obukhov length Bon96 p. 49
    real(r8) mno_stb_crc_heat ! [frc] Monin-Obukhov stability correction heat
    real(r8) mno_stb_crc_heat_crr ! Undamped correction factor heat [frc]
    real(r8) mno_stb_crc_mmn ! [frc] Monin-Obukhov stability correction momentum
    real(r8) mno_stb_crc_mmn_crr  ! Undamped correction factor momentum [frc]
    real(r8) sml_fnc_mmn_uns_rcp  ! [frc] Reciprocal of similarity function for momentum, unstable atmosphere
    real(r8) mno_stb_crc_tmp2     ! Term in stability correction computation
    real(r8) mno_stb_crc_tmp3     ! Term in stability correction computation
    real(r8) mno_stb_crc_tmp4     ! Term in stability correction computation
    real(r8) mno_stb_crc_tmp5     ! Term in stability correction computation
    real(r8) mno_stb_prm   ! [frc] Monin-Obukhov stability parameter 
    real(r8) mno_stb_prm_old ! [frc] Monin Obukhov stability parameter old
    real(r8) msv_gnd       ! [frc] Ground emissivity
    real(r8) nrg_bdg       ! [W m-2] Surface energy budget
    real(r8) nrg_bdg_dlt          ! [W m-2 K-1] Temperature derivative of surface energy budget
    real(r8) ppr_H2O_cnp   ! [Pa] Canopy vapor pressure of H2O
    real(r8) ppr_H2O_cnp_cff_a ! [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
    real(r8) ppr_H2O_cnp_cff_b ! [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55
    real(r8) ppr_H2O_mdp   ! [Pa] Ambient vapor pressure of H2O
    real(r8) rgh_heat      ! [m] Roughness length heat
    real(r8) rss_aer_fct   ! [s m-1] Term in resistance calculation
    real(r8) rss_aer_heat  ! [s m-1] Aerodynamic resistance to heat transfer
    real(r8) rss_aer_heat_fct ! [frc] Term in resistance calculation
    real(r8) rss_aer_mmn   ! [s m-1] Aerodynamic resistance to momentum transfer
    real(r8) rss_aer_mmn_fct ! [frc] Term in resistance calculation
    real(r8) rss_aer_vpr   ! [s m-1] Aerodynamic resistance to vapor transfer
    real(r8) rss_sfc_vpr   ! [s m-1] Surface resistance to vapor transfer
    real(r8) svp_H2O_gnd   ! [Pa] Saturation vapor pressure over planar condensed water at ground
    real(r8) tpt_bnd_cls          ! [C] Temperature bounded celsius
    real(r8) tpt_gnd_dlt   ! [K] Ground temperature adjustment
    real(r8) tpt_vrt       ! [K] Virtual temperature
    !real(r8) wnd_mdp       ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_mdp_bnd   ! [m s-1] Surface layer mean wind speed bounded
    !// ------------------------------------------------------------------
    
    ! Main Code
    
    ! Initialize variables which normally would be available in LSM
    if (flg_mbl) then
       ppr_H2O_mdp = q_H2O_vpr*prs_mdp/(eps_H2O+one_mns_eps_H2O*q_H2O_vpr) ! [Pa] Ambient vapor pressure of H2O
       msv_gnd = sfc_ems_gnd ! [frc] Ground emissivity
       ! NB: Initial guess for stability parameter depends on last Monin-Obukhov length
       ! mno_stb_prm(lon_idx)=min((hgt_mdp(lon_idx)-hgt_zpd)/mno_lng(lon_idx),1.0_r8) ! [frc]
       mno_stb_prm = min(-(hgt_mdp-hgt_zpd)/15.0_r8,1.0_r8) ! [frc]
       rgh_heat = rgh_mmn ! [m] Roughness length heat
       tpt_vrt = tpt_mdp*(1.0_r8+eps_H2O_rcp_m1*q_H2O_vpr) ! [K] Virtual temperature
    endif                  ! endif flg_mbl

    
    ! Initialize variables which are independent of stability iteration
    if (flg_mbl) then
       ! Zero temperature adjustments from last timestep
       tpt_gnd_dlt = 0.0_r8 ! [K] Change in ground temperature
       ! Latent heat of water transformation
       if (tpt_mdp > tpt_frz_pnt) then 
          ltn_heat_trn = ltn_heat_vpr_H2O_std ! [J kg-1]
       else 
          ltn_heat_trn = ltn_heat_sbl_H2O_std ! [J kg-1]
       endif               ! endif
       ! Psychrometric constant for transformations of water vapor at surface
       cst_psych = spc_heat_dry_air*prs_mdp/(eps_H2O*ltn_heat_trn) ! [Pa K-1] 
    endif                  ! endif flg_mbl

    
    if (flg_mbl) then
       ! Midlayer wind speeds
       !wnd_mdp(lon_idx)=   & ! [m s-1] Surface layer mean wind speed
       !     sqrt(wnd_znl_mdp(lon_idx)*wnd_znl_mdp(lon_idx)+  &
       !     wnd_mrd_mdp(lon_idx)*wnd_mrd_mdp(lon_idx)) 
       wnd_mdp_bnd = max(wnd_mdp,wnd_min_mbl) ! [m s-1] Surface layer mean wind speed bounded
       ! Miscellaneous
       sgn_chg_ctr = 0 ! [nbr] Number of sign changes in stability parameter
       mno_stb_prm_old = 0.0_r8 ! [frc] Monin Obukhov stability parameter old
       ! Stability-independent components of aerodynamic resistance calculations
       rss_aer_fct = 1.0_r8/(cst_von_krm*cst_von_krm*wnd_mdp_bnd) ! [s m-1]
       rss_aer_mmn_fct = log((hgt_mdp-hgt_zpd)/rgh_mmn) ! [frc]
       rss_aer_heat_fct = log((hgt_mdp-hgt_zpd)/rgh_heat) ! [frc]
    endif                  ! endif flg_mbl


    if (flg_mbl) then
       ! F(LW net) = -msv*Fdwn + msv*sigma*Tg^4 = a + b*Tg^4 Bon96 p. 45
       flx_LW_net_cff_a = -msv_gnd*flx_LW_dwn_sfc ! [W m-2] a in FLWnet = a + b*Tg^4 Bon96 p. 45
       flx_LW_net_cff_b = msv_gnd*cst_Stefan_Boltzmann ! [W m-2 K-4] b in FLWnet = a + b*Tg^4 Bon96 p. 45
          
       ! F(sns heat dwn into soil) = 2*k*(Tg-T1)/dz = a + b*Tg Bon96 p. 64
       flx_sns_gnd_cff_b = 2.0_r8*cnd_trm_soi/lvl_dlt ! [W m-2 K-1] b in Fgnd = a + b*Tg Bon96 p. 64
       flx_sns_gnd_cff_a = -flx_sns_gnd_cff_b*tpt_soi ! [W m-2] a in Fgnd = a + b*Tg Bon96 p. 64
    endif                  ! endif flg_mbl

    
    ! Iteration loop
    if (flg_mbl) then
      do itr_idx=1,itr_max
       
       if (flg_mbl) then
          ! Stability functions computed as in Bon96 p. 52
          if (mno_stb_prm < 0.0_r8) then
             sml_fnc_mmn_uns_rcp = (1.0_r8-16.0_r8*mno_stb_prm)**0.25_r8
             mno_stb_crc_tmp2 = log((1.0_r8+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0_r8)
             mno_stb_crc_tmp3 = log((1.0_r8+sml_fnc_mmn_uns_rcp)/2.0_r8)
             mno_stb_crc_mmn_crr = 2.0_r8*mno_stb_crc_tmp3+mno_stb_crc_tmp2-2.0_r8*atan(sml_fnc_mmn_uns_rcp)+1.5707963_r8 ! [frc]
             mno_stb_crc_heat_crr = 2.0_r8*mno_stb_crc_tmp2 ! [frc]
          else             ! not stable
             mno_stb_crc_mmn_crr = -5.0_r8*mno_stb_prm ! [frc]
             mno_stb_crc_heat_crr = mno_stb_crc_mmn ! [frc]
          endif            ! not stable
             
          ! Filter stability corrections to reduce numerical ping-pong
          if (itr_idx == 1) then
             mno_stb_crc_mmn = mno_stb_crc_mmn_crr ! [frc]
             mno_stb_crc_heat = mno_stb_crc_heat_crr ! [frc]
          else
             mno_stb_crc_mmn = 0.5*(mno_stb_crc_mmn_crr+mno_stb_crc_mmn) ! [frc]
             mno_stb_crc_heat = 0.5*(mno_stb_crc_heat_crr+mno_stb_crc_heat) ! [frc]
          endif            ! endif first iteration
          
          ! Aerodynamic resistance between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
          mno_stb_crc_tmp4 = rss_aer_mmn_fct-mno_stb_crc_mmn ! [frc]
          mno_stb_crc_tmp5 = rss_aer_heat_fct-mno_stb_crc_heat ! [frc]
          rss_aer_mmn = max(rss_aer_fct*mno_stb_crc_tmp4*mno_stb_crc_tmp4,1.0_r8) ! [s m-1]
          rss_aer_heat = max(rss_aer_fct*mno_stb_crc_tmp4*mno_stb_crc_tmp5,1.0_r8) ! [s m-1]
          ! Resistances are equal because rgh_heat = rgh_vpr
          rss_aer_vpr = rss_aer_heat ! [s m-1]
             
          ! Exchange coefficients between z(gcm) and zpd+z0m, zpd+z0h, and zpd+z0w
          ! Exchange coefficients are dimensionless, inversely proportional to wind speed 
          cff_xch_mmn = 1.0_r8/(rss_aer_mmn*wnd_mdp_bnd) ! [frc]
          ! Friction velocity
          wnd_frc = wnd_mdp_bnd*sqrt(cff_xch_mmn) ! [m s-1]
       endif               ! endif flg_mbl
       

       if (flg_mbl) then
          ! Saturation vapor pressure of water at ground temperature
          tpt_bnd_cls = tpt_bnd_cls_get(tpt_gnd) ! [C]
          if (tpt_bnd_cls > 0.0_r8) then 
             svp_H2O_gnd = svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
             dsvpdt_H2O_gnd = dsvpdt_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa K-1]
          else
             svp_H2O_gnd = svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
             dsvpdt_H2O_gnd = dsvpdt_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa K-1]
          endif            ! endif frozen
       endif               ! endif flg_mbl


       ! Surface resistance to vapor transfer
       if (itr_idx == 1) then ! NB: Bon96 only computes this on the first iteration
          if (flg_mbl) then
             rss_sfc_vpr = & ! [s m-1] Surface resistance to vapor transfer Bon96 p. 56, Fgr. 17 p. 57
                  rss_aer_vpr *(1.0_r8-trn_fsh_vpr_soi_atm)/ &
                  max(trn_fsh_vpr_soi_atm,eps_dbz)
             ! Set minimum rss_sfc_vpr so that in bare ground case es = svp(tg)
             rss_sfc_vpr = max(rss_sfc_vpr,eps_dbz)
          endif            ! endif flg_mbl
       endif                  ! endif first iteration
       
       ! Heat conductances from ground to ambient air
       if (flg_mbl) then
          ! Conductances are dimensional, exact inverses of resistances
          ! Coefficients for sensible heat flux SH = a + b*Ts
          cnd_heat_sfc_mdp = 1.0_r8/rss_aer_heat ! [m s-1] Sensible heat conductance surface air to midlayer air Bon96 p. 60, Fig. 16 p. 57
          flx_sns_atm_fct = dns_mdp*spc_heat_dry_air*cnd_heat_sfc_mdp ! Bon96 p. 55
          flx_sns_atm_cff_a = -tpt_ptn_mdp*flx_sns_atm_fct ! [W m-2] a in SH = a + b*Ts Bon96 p. 55, 69, 70
          flx_sns_atm_cff_b = flx_sns_atm_fct ! [W m-2 K-1] b in SH = a + b*Ts Bon96 p. 55, 69, 70
       endif               ! endif flg_mbl


       if (flg_mbl) then
          ! Vapor conductances from ground to ambient air
          flx_ltn_fct = dns_mdp*spc_heat_dry_air/cst_psych ! 
          cnd_vpr_sfc_mdp = 1.0_r8/rss_aer_vpr ! [m s-1] Water vapor conductance surface air to midlayer air Bon96 p. 60, Fgr. 17 p. 57
          cnd_vpr_gnd_sfc = 1.0_r8/rss_sfc_vpr ! [m s-1] Water vapor conductance ground to surface air Bon96 p. 60, Fgr. 17 p. 57
          cnd_vpr_ttl = cnd_vpr_sfc_mdp+cnd_vpr_gnd_sfc ! [m s-1] Sum of conductances

          ! Coefficients for canopy vapor pressure e(cnp) = a + b*e(Ts)
          ppr_H2O_cnp_cff_a = ppr_H2O_mdp*cnd_vpr_sfc_mdp/cnd_vpr_ttl ! [Pa] a in e(cnp) = a + b*e(Ts) Bon96 p. 55
          ppr_H2O_cnp_cff_b = cnd_vpr_gnd_sfc/cnd_vpr_ttl ! [frc] b in e(cnp) = a + b*e(Ts) Bon96 p. 55
             
          ! Coefficients for evaporation heat flux LHE = a + b*e(Ts) 
          flx_ltn_evp_cff_a = -flx_ltn_fct*(ppr_H2O_mdp-ppr_H2O_cnp_cff_a)*cnd_vpr_sfc_mdp ! [W m-2] a in LHE = a + b*e(Ts) Bon96 p. 55
          flx_ltn_evp_cff_b = flx_ltn_fct*ppr_H2O_cnp_cff_b*cnd_vpr_sfc_mdp ! [W m-2 Pa-1] b in LHE = a + b*e(Ts) Bon96 p. 55
       endif               ! endif flg_mbl

       
       ! Evaluate fluxes for current tpt_gnd
       if (flg_mbl) then
          flx_sns_gnd = flx_sns_gnd_cff_a+flx_sns_gnd_cff_b*tpt_gnd ! [W m-2] Sensible heat flux to soil
          flx_sns_atm = flx_sns_atm_cff_a+flx_sns_atm_cff_b*tpt_gnd ! [W m-2] Sensible heat flux to atmosphere
          flx_ltn_evp = flx_ltn_evp_cff_a+flx_ltn_evp_cff_b*svp_H2O_gnd ! [W m-2] Evaporation flux to atmosphere
          flx_LW_net = flx_LW_net_cff_a+flx_LW_net_cff_b*(tpt_gnd**4.0_r8) ! [W m-2] Net longwave flux to atmosphere
          nrg_bdg = & ! [W m-2] Total energy budget at surface
               +flx_SW_abs_sfc &
               -flx_LW_net &
               -flx_ltn_evp &
               -flx_sns_atm &
               -flx_sns_gnd
          nrg_bdg_dlt=     & ! [W m-2 K-1] Temperature derivative of surface energy budget
               -4.0_r8*flx_LW_net_cff_b*(tpt_gnd**3.0_r8) &
               -flx_ltn_evp_cff_b*dsvpdt_H2O_gnd &
               -flx_sns_atm_cff_b &
               -flx_sns_gnd_cff_b
          tpt_gnd_dlt = -nrg_bdg/nrg_bdg_dlt ! [K] Newton-Raphson temperature adjustment
       endif               ! endif flg_mbl


       if (flg_mbl) then
          ! Adjust temperatures 
          tpt_gnd = tpt_gnd+tpt_gnd_dlt ! [K]
             
          ! Adjust canopy vapor pressure
          tpt_bnd_cls = tpt_bnd_cls_get(tpt_gnd) ! [C]
          if (tpt_bnd_cls > 0.0_r8) then
             ppr_H2O_cnp = ppr_H2O_cnp_cff_a+ppr_H2O_cnp_cff_b*svp_H2O_lqd_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
          else
             ppr_H2O_cnp = ppr_H2O_cnp_cff_a+ppr_H2O_cnp_cff_b*svp_H2O_ice_PrK78_fst_scl(tpt_bnd_cls) ! [Pa]
          endif            ! endif frozen
          
          ! Sensible heat flux
          flx_sns_atm_tmp = & ! [W m-2] Temporary sensible heat flux
               -(tpt_ptn_mdp-tpt_gnd) &
               *dns_mdp*spc_heat_dry_air &
               /rss_aer_heat
          ! Following step approximates psi_h(z0m/L) = 0 
          ! Monin-Obukhov stability parameter mno_stb_prm for next iteration
          flx_vpr_tmp = & ! [kg m-2 s-1] Temporary water vapor flux Bon96 p. 56
               ! 19990501 fxm: this appears to be latent heat flux, not water vapor flux
               -(ppr_H2O_mdp-ppr_H2O_cnp) &
               *dns_mdp*spc_heat_dry_air &
               /(cst_psych*rss_aer_vpr)
          flx_sns_atm_vrt_tmp = & ! [W m-2] Virtual sensible heat flux Bon96 p. 49
               +flx_sns_atm_tmp &
               +eps_H2O_rcp_m1*spc_heat_dry_air*tpt_mdp*flx_vpr_tmp &
               /ltn_heat_trn
          mno_dnm =         & ! Denominator of Monin-Obukhov length Bon96 p. 49
               +cst_von_krm &
               *(grv_sfc/tpt_vrt) &
               *flx_sns_atm_vrt_tmp &
               /(dns_mdp*spc_heat_dry_air)
          ! Set denominator of Monin-Obukhov length to minimum value if vapor and heat fluxes equal 0.0
          if (abs(mno_dnm) <= eps_dbz) mno_dnm=eps_dbz
          mno_lng = -1.0_r8*(wnd_frc**3.0_r8)/mno_dnm ! [m] Monin-Obukhov length Bon96 p. 49
          ! Stability functions only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
          mno_stb_prm = min((hgt_mdp-hgt_zpd)/mno_lng,1.0_r8) ! [frc] Monin Obukhov stability parameter
             
          ! Accumulate number of times stability parameter changes sign 
          if (mno_stb_prm_old*mno_stb_prm < 0.0_r8) sgn_chg_ctr=sgn_chg_ctr+1
          ! Zero stability parameter if it has changed sign too many times
          if (sgn_chg_ctr >= sgn_chg_ctr_max) then
             mno_stb_prm = 0.0_r8 ! [frc]
             ! Zero stability corrections for consistency with stability parameter
             mno_stb_crc_mmn = 0.0_r8 ! [frc]
             mno_stb_crc_heat = 0.0_r8 ! [frc]
          endif           ! endif
          mno_stb_prm_old = mno_stb_prm ! [frc]
       endif               ! endif flg_mbl
       
      end do                     ! end loop over itr
    endif !// if (flg_mbl) then

    !// ------------------------------------------------------------------
  end subroutine blm_mbl                       ! end blm_mbl()
  !// ------------------------------------------------------------------


end module dstblm

! Purpose: dstmbl.F90 controls mineral dust mobilization processes
! Usage:
! use dstmbl ! [mdl] Mobilization driver
!// ------------------------------------------------------------------
!// Rewritten as a box model for Oslo CTM3.
!//
!// Amund Sovde, February 2015
!// Also rewritten to take dust mass per grid box as input, instead
!// of mass mixing ratio. The old method of using mixing ratio and
!// calculated air mass and density created small (?) inconsistencies
!// between flux diagnostics and the actual change in tracer mass.
!// ------------------------------------------------------------------

! Requires dst.h for DST_MSS_BDG
!#include <dst.h> /* Dust preprocessor tokens */

module dstmbl ! [mdl] Mobilization driver
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_mbl ! [sbr] Driver for aerosol mobilization

contains
  
  !// ------------------------------------------------------------------
  subroutine dst_mbl(lchnk,ncol,obuf, &
       doy,                 & ! I [day] Day of year [1.0..366.0)
       hgt_mdp,             & ! I [m] Midpoint height above surface
       lat_idx,             & ! I [idx] Model latitude index
       lon_idx,             & ! I [idx] Model longitude index
       lat_rdn,             & ! I [rdn] Latitude
       oro,                 & ! I [frc] Orography
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       prs_sfc,             & ! I [Pa] Surface pressure
       q_H2O_vpr,           & ! I [kg kg-1] Water vapor mixing ratio
       m_dst,               & ! I/O [kg] Dust mass
       airm,                & ! I [kg] Air mass
       volu,                & ! I [m3] Air volume
       area,                & ! I [m2] Grid box area
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp,             & ! I [K] Temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
       wnd_mrd_mdp,         & ! I [m s-1] Meridional wind component
       wnd_znl_mdp,         & ! I [m s-1] Zonal wind component
       production)            ! O [kg] Dust produced
    !// ------------------------------------------------------------------
    ! Purpose: Driver for aerosol mobilization
    ! dst_mbl() is called by CCM:physics/tphysac(), MATCH:src/physlic()
    ! NB: dst_mbl() is designed to require only single layer surface fields 
    ! Eliminating multi-layer fields allows easier implementation in LSM
    !//
    !// Rewritten for Oslo CTM3, calculates only value for (lon_idx,lat_idx).
    !//
    !// Also rewritten to take dust mass per grid box as input, instead
    !// of mass mixing ratio. The old method of using mixing ratio and
    !// calculated air mass and density created small (?) inconsistencies
    !// between flux diagnostics and the actual change in tracer mass.
    !//
    !// Amund Sovde, February 2015, October 2009
    !// ------------------------------------------------------------------
    ! [mdl] Boundary layer meteorology driver
    use blmutl,only: wnd_rfr_get, snw_frc_get
    ! [mdl] Aerosol microphysical properties
    use dstaer,only: dmt_vwr,dns_aer, mss_frc_src, mss_frc_trn_dst_src, &
                     ovr_src_snk_mss, ovr_src_snk_mss_ttl
    ! [mdl] Mass budget diagnostics
    !use dstbdg,only:bdg_aal,bdg_gam_wet_2d
    ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use dstblm,only: blm_mbl, cnd_trm_soi_get, hyd_prp_get, tpt_frz_pnt, &
                     trn_fsh_vpr_soi_atm_get
    ! [mdl] Physical constants for dust routines
    use dstcst
    ! [mdl] Debugging information for dust model
    use dstdbg
    ! [mdl] Dust grid sizes
    use dstgrd,only: dst_nbr, bln_nbr, dst_src_nbr
    ! [mdl] Mobilization utilities
    use dstmblutl
    ! [mdl] Mass budget utilities
    use dstmssutl,only: dst_add_nbr
    ! [mdl] Nomenclature for outfld()
    use dstnm
    ! [mdl] 2-D surface fields on PLON x PLAT grid
    use dstsfc,only: sfc_typ_get, soi_txt_get, vwc_sfc_get, tpt_gnd_soi_get, flx_rad_sfc_get
    ! [mdl] Saltation sandblasting physics
    use dstsltsbl
    ! [mdl] Time-varying boundary data sets
    use dsttvbds,only:dst_tvbds_get
!#ifndef BXM
    ! [mdl] History/archiving
    use dead_history,only: outfld_1
!#endif /* BXM */
    ! [mdl] String manipulation
    use sng_mdl,only: ftn_strnul
    ! [mdl] Weibull wind speed distribution
    use wbl_mdl,only: wbl_wnd, wnd_mdp_bin_get
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! [m] Reference height for mobilization processes
    real(r8),parameter::hgt_rfr=10.0_r8
    ! [m] Zero plane displacement for erodible surfaces
    real(r8),parameter::hgt_zpd_mbl=0.0_r8
#ifndef AlG01
    ! [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
    real(r8),parameter::rgh_mmn_mbl=100.0e-6_r8
#else
    ! [m] Roughness length momentum for erodible surface (AlG01 formulation)
    real(r8),parameter::rgh_mmn_mbl=3000.0e-6_r8
#endif
    ! [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
    ! fxm: rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
    real(r8),parameter::rgh_mmn_smt=33.3e-6_r8
    real(r8),parameter::wnd_min_mbl=1.0_r8 ! [m s-1] Minimum windspeed used for mobilization 

    ! Input
    integer,intent(in) :: lat_idx, & ! I [idx] Model latitude index
                          lon_idx, & ! I [idx] Model longitude index
                          lchnk, &   ! I [id] Chunk identifier
                          ncol       ! I [nbr] Number of atmospheric columns
    real(r8),intent(in) :: doy, &         ! I [day] Day of year [1.0..366.0)
                           hgt_mdp, &     ! I [m] Midlayer height above surface
                           lat_rdn, &     ! I [rdn] Latitude
                           obuf(*), &     ! I [ptr] Output buffer
                           oro, &         ! I [frc] Orography
                           prs_dlt, &     ! I [Pa] Pressure thickness
                           prs_mdp , &    ! I [Pa] Pressure
                           prs_sfc, &     ! I [Pa] Surface pressure
                           q_H2O_vpr, &   ! I [kg kg-1] Water vapor mixing ratio
                           snw_hgt_lqd, & ! I [m] Equivalent liquid water snow depth
                           tm_adj, &      ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
                           tpt_mdp, &     ! I [K] Temperature
                           tpt_ptn_mdp, & ! I [K] Midlayer local potential temperature
                           wnd_mrd_mdp, & ! I [m s-1] Meridional wind component
                           wnd_znl_mdp, & ! I [m s-1] Zonal wind component
                           airm, &        ! I [kg] Air mass @ sfc
                           volu, &        ! I [m3] Air volume @ sfc
                           area           ! I [m2] Grid box area

    ! Input/Output
    real(r8),intent(inout) :: m_dst(dst_nbr) ! I/O [kg] Dust mass
    ! Output
    real(r8), intent(out) :: production(dst_nbr) ! O [kg] Dust mass produced

    ! Local Output
    real(r8) :: dst_slt_flx_rat_ttl, &      ! [m-1] Ratio of vertical dust flux to streamwise mass flux
                flx_mss_hrz_slt_ttl, &      ! [kg m-1 s-1] Vertically integrated streamwise mass flux
                flx_mss_vrt_dst(dst_nbr), & ! [kg m-2 s-1] Vertical mass flux of dust
                flx_mss_vrt_dst_ttl, &      ! [kg m-2 s-1] Total vertical mass flux of dust
                frc_thr_ncr_drg, &          ! [frc] Threshold friction velocity increase from roughness
                frc_thr_ncr_wtr, &          ! [frc] Threshold friction velocity increase from moisture
                hgt_zpd, &                  ! [m] Zero plane displacement
                lnd_frc_mbl, &              ! [frc] Bare ground fraction
                mno_lng, &                  ! [m] Monin-Obukhov length
                snw_frc, &                  ! [frc] Fraction of surface covered by snow
                trn_fsh_vpr_soi_atm, &      ! [frc] Transfer efficiency of vapor from soil to atmosphere
                wnd_frc, &                  ! [m s-1] Friction velocity
                wnd_frc_slt, &              ! [m s-1] Saltating friction velocity
                wnd_frc_thr_slt, &          ! [m s-1] Threshold friction velocity for saltation
                wnd_mdp, &                  ! [m s-1] Surface layer mean wind speed
                wnd_rfr, &                  ! [m s-1] Wind speed at reference height
                wnd_rfr_thr_slt             ! [m s-1] Threshold 10 m wind speed for saltation

    ! Local
    logical :: flg_CaCO3, &        ! [flg] Activate CaCO3 tracer
               flg_mbl             ! [flg] Mobilization candidates
    character(80) :: fl_out        ! [sng] Name of netCDF output file
    integer :: m, &                ! [idx] Counting index for species
               mbl_nbr, &          ! [nbr] Number of mobilization candidates
               sfc_typ             ! [idx] LSM surface type (0..28)
    real(r8) :: cnd_trm_soi, &     ! [W m-1 K-1] Soil thermal conductivity
                dns_mdp, &         ! [kg m-3] Midlayer density
                flx_LW_dwn_sfc, &  ! [W m-2] Longwave downwelling flux at surface
                flx_SW_abs_sfc, &  ! [W m-2] Solar flux absorbed by ground
                lnd_frc_dry, &     ! [frc] Dry land fraction
                mbl_bsn_fct , &    ! [frc] Erodibility factor
                mss_frc_CaCO3, &   ! [frc] Mass fraction of CaCO3
                lvl_dlt, &         ! [m] Soil layer thickness
                mpl_air, &         ! [kg m-2] Air mass path in layer
                mss_frc_cly, &     ! [frc] Mass fraction of clay
                mss_frc_snd, &     ! [frc] Mass fraction of sand
                sfc_frc_bln(bln_nbr), & ! [frc] fraction of each soil type
                tm_dlt, &               ! [s] Mobilization timestep
                tpt_gnd, &              ! [K] Ground temperature
                tpt_soi, &              ! [K] Soil temperature
                tpt_soi_frz, &          ! [K] Temperature of frozen soil
                !tpt_vrt_mdp , &         ! [K] Midlayer virtual temperature
                vai_dst, &         ! [m2 m-2] Vegetation area index, one-sided
                vwc_dry, &         ! [m3 m-3] Dry volumetric water content (no E-T)
                vwc_opt, &         ! [m3 m-3] E-T optimal volumetric water content 
                vwc_sat, &         ! [m3 m-3] Saturated volumetric water content (sand-dependent)
                vwc_sfc, &         ! [m3 m-3] Volumetric water content
                gwc_sfc            ! [kg kg-1] Gravimetric water content

    !// In case we need mixing ratio
    real(r8) :: q_dst(dst_nbr) ! [kg kg-1] Dust mixing ratio

    ! [frc] Global mass flux tuning factor (a posteriori)
    real(r8) :: flx_mss_fdg_fct

    ! ++alfgr variables needed for Weibull code
    integer,parameter :: wnd_mdp_nbr=5 ! [nbr] Number of discrete wind speeds
    integer :: wnd_mdp_idx ! [idx] Counter for wind speed discretization
    real(r8),parameter :: wnd_frc_rsl=0.95_r8    ! Fraction of wind PDF to resolve
    real(r8) :: wnd_mdp_wgt(wnd_mdp_nbr), &    ! Weighting for winds
                wnd_mdp_bin, &                 ! Wind in the actual wind bin
                wnd_mdp_min, &                 ! Mininum wind speed we are interested in
                wnd_mdp_max, &                 ! Maximum wind speed we are interested in
                ovr_src_snk_mss_wbn(dst_src_nbr,dst_nbr), &   ! Overlap mass factors for wbin
                ovr_src_snk_mss_add(dst_src_nbr,dst_nbr), &   ! Overlap mass factors (added up)
                ovr_src_snk_mss_ttl_wbn, &    ! Total overlap between source and sink from wbin
                ovr_src_snk_mss_ttl_add, &    ! Total overlap between source and sink (added up)
                flx_mss_hrz_slt_ttl_wbn, &    ! Horizontal flux due to a wind bin    
                flx_mss_vrt_dst_ttl_wbn, &    ! Vertical flux due to a wind bin
                mss_frc_src_wbn(dst_src_nbr), &     ! Mass fraction source due to this wind bin
                mss_frc_src_add(dst_src_nbr), &     ! Mass fraction source (added up)
                mss_frc_trn_dst_src_wbn(dst_nbr), & ! Mass fraction of source transported from wbin
                mss_frc_trn_dst_src_add(dst_nbr)    ! Mass fraction of source transported (added up)
    ! --alfgr variables needed for Weibull winds

    ! GCM diagnostics
    real(r8) :: q_dst_tnd_mbl(dst_nbr), &  ! [kg kg-1 s-1] Dust tendency due to gravitational settling
                q_dst_tnd_mbl_ttl          ! [kg kg-1 s-1] Total dust tendency due to gravitational settling
    !// ------------------------------------------------------------------
    
    ! Main Code
    ! Timesplit if desired
    tm_dlt = tm_adj ! [s] (default CCM: 2*dt, MATCH: dt)

    ! Assume water in soil freezes at 0 C
    tpt_soi_frz = tpt_frz_pnt ! [K] Temperature of frozen soil

    ! Global mass flux tuning factor (a posteriori, from dstcst.F90)
    flx_mss_fdg_fct = flx_mss_fdg_fct0

    ! Initialize output fluxes and tendencies
    q_dst_tnd_mbl(:)   = 0.0_r8 ! [kg kg-1 s-1]
    q_dst_tnd_mbl_ttl  = 0.0_r8 ! [kg kg-1 s-1]
    flx_mss_vrt_dst(:) = 0.0_r8 ! [kg m-2 s-1]
    flx_mss_vrt_dst_ttl= 0.0_r8 ! [kg m-2 s-1]
    frc_thr_ncr_wtr = 0.0_r8 ! [frc]
    wnd_rfr         = 0.0_r8 ! [m s-1]
    wnd_frc         = 0.0_r8 ! [m s-1]
    wnd_frc_slt     = 0.0_r8 ! [m s-1]
    wnd_frc_thr_slt = 0.0_r8 ! [m s-1]
    wnd_rfr_thr_slt = 0.0_r8 ! [m s-1]
    hgt_zpd         = hgt_zpd_mbl ! [m]
    flx_mss_hrz_slt_ttl = 0.0_r8

    ! Wind PDF initialization
    ! Initialize mass fraction of source
    mss_frc_src_wbn(:) = mss_frc_src(:)
    mss_frc_src_add(:) = 0.0_r8 ! Must be zero because it is updated and weighted
    ! Initialize overlap of source and sink
    ovr_src_snk_mss_wbn(:,:) = ovr_src_snk_mss(:,:)
    ovr_src_snk_mss_add(:,:) = 0.0_r8 ! Must start as zero because it is updated and weighted
    ! Initialize total overlap src/sink
    ovr_src_snk_mss_ttl_wbn = ovr_src_snk_mss_ttl
    ovr_src_snk_mss_ttl_add = 0.0_r8 ! Must start as zero because it is updated and weighted
    mss_frc_trn_dst_src_wbn(:) = mss_frc_trn_dst_src(:)
    mss_frc_trn_dst_src_add(:) = 0.0_r8 ! Must start as zero because it is updated and weighted      
    
    ! Compute required derived fields
    if(tpt_mdp > 350.0_r8) stop 'ERROR: dst_mbl() reports tpt_mdp(i) > 350.0'
    ! Midlayer virtual temperature
    !//tpt_vrt_mdp = tpt_mdp*(1.0_r8 + eps_H2O_rcp_m1 * q_H2O_vpr) ! [K]

    !// Density at center of gridbox, now from CTM air mass and volume.
    !//dns_mdp = prs_mdp/(tpt_vrt_mdp*gas_cst_dry_air) ! [kg m-3]
    dns_mdp = airm / volu

    ! Approximate surface virtual temperature (uses midlayer moisture)
    ! tpt_vrt_sfc=tpt_sfc(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
    ! Surface density
    ! dns_sfc(i)=prs_sfc(i)/(tpt_vrt_sfc*gas_cst_dry_air) ! [kg m-3]

    !// Mass of air currently in gridbox
    !//mpl_air = prs_dlt*grv_sfc_rcp ! [kg m-2]
    mpl_air = airm / area

    ! Mean surface layer horizontal wind speed 
    wnd_mdp = sqrt(wnd_znl_mdp*wnd_znl_mdp + wnd_mrd_mdp*wnd_mrd_mdp)
       
    !write(6,*)'WIND SPEED AT MIDPOINT',wnd_mdp(i)

    
    ! Gather input variables from common blocks
    ! Surface type (checked ok for CTM3)
    call sfc_typ_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         sfc_typ)               ! O [idx] LSM surface type (0..28)


    ! Soil texture and dry land fraction (checked ok for CTM3)
    call soi_txt_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         lnd_frc_dry,         & ! O [frc] Dry land fraction
         mbl_bsn_fct,         & ! O [frc] Erodibility factor
         mss_frc_CaCO3,       & ! O [frc] Mass fraction of CaCO3
         mss_frc_cly,         & ! O [frc] Mass fraction of clay
         mss_frc_snd,         & ! O [frc] Mass fraction of sand
         sfc_frc_bln)           ! O [frc] Fraction of four soil types


    ! Soil moisture (checked ok for CTM3)
    call vwc_sfc_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         vwc_sfc)               ! O [m3 m-3] Volumetric water content


    ! Surface and soil temperature (checked ok for CTM3)
    call tpt_gnd_soi_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         tpt_gnd,             & ! O [K] Ground temperature
         tpt_soi)               ! O [K] Soil temperature


    ! Radiative fluxes at surface (checked ok for CTM3)
    call flx_rad_sfc_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         flx_LW_dwn_sfc,      & ! O [W m-2] Longwave downwelling flux at surface
         flx_SW_abs_sfc)        ! O [W m-2] Solar flux absorbed by ground


    ! Time varying boundary datasets: Vegetation area index (checked ok for CTM3)
    call dst_tvbds_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         vai_dst)               ! O [m2 m-2] Vegetation area index, one-sided


    ! All variables from common blocks have been assembled
    ! Use these to compute derived fields


    ! Fraction of surface covered by snow (checked ok for CTM3)
    call snw_frc_get( &
         snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
         snw_frc)               ! O [frc] Fraction of surface covered by snow


    ! Fraction of each gridcell suitable for dust mobilization (checked ok for CTM3)
    call lnd_frc_mbl_get( &
         doy,                 & ! I [day] Day of year [1.0..366.0)
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         flg_mbl,             & ! O [flg] Mobilization candidate flag
         lat_rdn,             & ! I [rdn] Latitude
         lnd_frc_dry,         & ! I [frc] Dry land fraction
         lnd_frc_mbl,         & ! O [frc] Bare ground fraction
         mbl_nbr,             & ! O [flg] Number of mobilization candidates
         oro,                 & ! I [frc] Orography
         sfc_typ,             & ! I [idx] LSM surface type (0..28)
         snw_frc,             & ! I [frc] Fraction of surface covered by snow
         tpt_soi,             & ! I [K] Soil temperature
         tpt_soi_frz,         & ! I [K] Temperature of frozen soil
         vai_dst)               ! I [m2 m-2] Vegetation area index, one-sided
    
    ! Much ado about nothing
    if (mbl_nbr == 0) goto 737

    
    ! Hydrologic properties (checked ok for CTM3)
    call hyd_prp_get(         & ! NB: These properties are time-invariant
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         mss_frc_cly,         & ! I [frc] Mass fraction clay 
         mss_frc_snd,         & ! I [frc] Mass fraction sand
         vwc_dry,             & ! O [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt,             & ! O [m3 m-3] E-T optimal volumetric water content
         vwc_sat)               ! O [m3 m-3] Saturated volumetric water content

    
    ! Soil thermal conductivity and layer thickness (checked ok for CTM3)
    call cnd_trm_soi_get( &
         cnd_trm_soi,         & ! O [W m-1 K-1] Soil thermal conductivity
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         lvl_dlt,             & ! O [m] Soil layer thickness
         mss_frc_cly,         & ! I [frc] Mass fraction clay 
         mss_frc_snd,         & ! I [frc] Mass fraction sand
         tpt_soi,             & ! I [K] Soil temperature
         vwc_dry,             & ! I [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt,             & ! I [m3 m-3] E-T optimal volumetric water content
         vwc_sat,             & ! I [m3 m-3] Saturated volumetric water content
         vwc_sfc)               ! I [m3 m-3] Volumetric water content


    ! Transfer efficiency of vapor from soil to atmosphere (checked ok for CTM3)
    call trn_fsh_vpr_soi_atm_get( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         tpt_soi,             & ! I [K] Soil temperature
         tpt_soi_frz,         & ! I [K] Temperature of frozen soil
         trn_fsh_vpr_soi_atm, & ! O [frc] Transfer efficiency of vapor from soil to atmosphere
         vwc_dry,             & ! I [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt,             & ! I [m3 m-3] E-T optimal volumetric water content
         vwc_sfc)               ! I [m3 m-3] Volumetric water content


    ! Surface exchange properties over erodible surfaces (checked ok for CTM3)
    call blm_mbl( &
         cnd_trm_soi,         & ! I [W m-1 K-1] Soil thermal conductivity
         dns_mdp,             & ! I [kg m-3] Midlayer density
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         flx_LW_dwn_sfc,      & ! I [W m-2] Longwave downwelling flux at surface
         flx_SW_abs_sfc,      & ! I [W m-2] Solar flux absorbed by ground
         hgt_mdp,             & ! I [m] Midlayer height above surface
         hgt_zpd_mbl,         & ! I [m] Zero plane displacement
         lvl_dlt,             & ! I [m] Soil layer thickness
         prs_mdp,             & ! I [Pa] Pressure
         q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
         rgh_mmn_mbl,         & ! I [m] Roughness length momentum
         tpt_mdp,             & ! I [K] Midlayer temperature
         tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
         tpt_soi,             & ! I [K] Soil temperature
         trn_fsh_vpr_soi_atm, & ! I [frc] Transfer efficiency of vapor from soil to atmosphere
         wnd_mdp,             & ! I [m s-1] wind in layer midpoint
         mno_lng,             & ! O [m] Monin-Obukhov length
         tpt_gnd,             & ! I/O [K] Ground temperature
         wnd_frc)               ! O [m s-1] Surface friction velocity


    ! Interpolate midlayer wind speed to 10 m (checked ok for CTM3)
    call wnd_rfr_get( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         hgt_mdp,             & ! I [m] Midpoint height above surface
         hgt_rfr,             & ! I [m] Reference height for mobilization processes
         hgt_zpd,             & ! I [m] Zero plane displacement
         mno_lng,             & ! I [m] Monin-Obukhov length
         wnd_frc,             & ! I [m s-1] Surface friction velocity
         wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
         wnd_min_mbl,         & ! I [m s-1] Minimum windspeed
         wnd_rfr)               ! O [m s-1] Wind speed at reference height


    ! Get the weighting factors and winds we want to calculate for (checked ok for CTM3)
    call wbl_wnd( &
         hgt_mdp, & ! I [m] Height of layer midpoint
         hgt_rfr, & ! I [m] Reference height
         wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
         wnd_mdp, & ! I [m s-1] Wind speed at midpoint
         wnd_mdp_min, & ! O [m s-1] Min wind speed at midpoint
         wnd_mdp_max, & ! O [m s-1] Max wind speed at midpoint 
         wnd_mdp_nbr, & ! I [nbr] Number of wind speeds to calculate 
         wnd_mdp_wgt, & ! O [frc] Wind speed bin weight
         wnd_rfr)       ! I [m s-1] Wind speed at reference height



    ! Now, do the loop on wnd_mdp
    do wnd_mdp_idx=1,wnd_mdp_nbr

       if(wnd_mdp_nbr > 1) then  ! More than one wind bin

          ! Get mean wind speed in this bin (checked ok for CTM3)
          call wnd_mdp_bin_get( &
               wnd_mdp_bin, & ! O [m s-1] wind speed in bin
               wnd_mdp_idx, & ! I [nbr] The bin we are interested in
               wnd_mdp_max, & ! I [m s-1] largest wind speed which we calculate
               wnd_mdp_min, & ! I [m s-1] smallest wind speed which we calculate
               wnd_mdp_nbr)   ! I [nbr] Weibull wind number
    
          !write(6,*)'  '
          !write(6,*)'WIND SPEED',wnd_mdp_bin,'@ height',hgt_mdp,'weight',wnd_mdp_wgt(:,wnd_mdp_idx)
          
          ! Get wind friction speed for this wind speed (checked ok for CTM3)
          call blm_mbl( &
               cnd_trm_soi,         & ! I [W m-1 K-1] Soil thermal conductivity
               dns_mdp,             & ! I [kg m-3] Midlayer density
               flg_mbl,             & ! I [flg] Mobilization candidate flag
               flx_LW_dwn_sfc,      & ! I [W m-2] Longwave downwelling flux at surface
               flx_SW_abs_sfc,      & ! I [W m-2] Solar flux absorbed by ground
               hgt_mdp,             & ! I [m] Midlayer height above surface
               hgt_zpd_mbl,         & ! I [m] Zero plane displacement
               lvl_dlt,             & ! I [m] Soil layer thickness
               prs_mdp,             & ! I [Pa] Pressure
               q_H2O_vpr,           & ! I [kg kg-1] Specific humidity
               rgh_mmn_mbl,         & ! I [m] Roughness length momentum
               tpt_mdp,             & ! I [K] Midlayer temperature
               tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
               tpt_soi,             & ! I [K] Soil temperature
               trn_fsh_vpr_soi_atm, & ! I [frc] Transfer efficiency of vapor from soil to atmosphere
               wnd_mdp_bin,         & ! I [m s-1] wind speed in layer midpoint
               mno_lng,             & ! O [m] Monin-Obukhov length
               tpt_gnd,             & ! I/O [K] Ground temperature
               wnd_frc)               ! O [m s-1] Surface friction velocity

       else ! Only one wind speed
          ! Set weights to unity then continue as before using mean quantities
          wnd_mdp_wgt(:)=1.0_r8 ! [frc] Wind speed bin weight
       endif ! endif 


       ! Factor by which surface roughness increases threshold friction velocity (checked ok for CTM3)
       call frc_thr_ncr_drg_get( &
            frc_thr_ncr_drg, & ! O [frc] Factor by which surface roughness increases threshold friction velocity
            rgh_mmn_mbl, & ! I [m] Roughness length momentum for erodible surfaces
            rgh_mmn_smt) ! I [m] Smooth roughness length
    

       ! Convert volumetric water content to gravimetric water content (checked ok for CTM3)
       call vwc2gwc( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            gwc_sfc, & ! O [kg kg-1] Gravimetric water content
            vwc_sat, & ! I [m3 m-3] Saturated volumetric water content (sand-dependent)
            vwc_sfc) ! I [m3 m-3] Volumetric water content


       ! Factor by which soil moisture increases threshold friction velocity (checked ok for CTM3)
       call frc_thr_ncr_wtr_get( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            frc_thr_ncr_wtr, & ! O [frc] Factor by which moisture increases threshold friction velocity
            mss_frc_cly, & ! I [frc] Mass fraction of clay
            gwc_sfc)             ! I [kg kg-1] Gravimetric water content


       ! Threshold friction velocity for saltation over dry, bare, smooth ground
       ! fxm: Use surface density not midlayer density (checked ok for CTM3)
       call wnd_frc_thr_slt_get( &
            dmt_vwr, & ! I [m] Mass weighted diameter resolved 
            dns_aer, & ! I [kg m-3] Particle density
            dns_mdp, & ! I [kg m-3] Midlayer density
            wnd_frc_thr_slt)     ! O [m s-1] Threshold friction velocity for saltation
    
       ! Adjust threshold friction velocity to account for moisture and roughness
       wnd_frc_thr_slt =      & ! [m s-1] Threshold friction velocity for saltation
            wnd_frc_thr_slt   & ! [m s-1] Threshold for dry, flat ground
            * frc_thr_ncr_wtr & ! [frc] Adjustment for moisture
            * frc_thr_ncr_drg   ! [frc] Adjustment for roughness

    
       ! Threshold saltation wind speed
       if (flg_mbl) then
          wnd_rfr_thr_slt = & ! [m s-1] Threshold 10 m wind speed for saltation
               wnd_rfr * wnd_frc_thr_slt / wnd_frc
       endif !// endif flg_mbl


       ! Saltation increases friction speed by roughening surface (checked ok for CTM3)
       call wnd_frc_slt_get( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            wnd_frc, & ! I [m s-1] Surface friction velocity
            wnd_frc_slt, & ! O [m s-1] Saltating friction velocity
            wnd_rfr, & ! I [m s-1] Wind speed at reference height
            wnd_rfr_thr_slt) ! I [m s-1] Threshold 10 m wind speed for saltation


#ifdef AlG01
       ! Horizontal streamwise mass flux consistent with AlG01 formulation (checked ok for CTM3)
       call flx_mss_hrz_slt_ttl_AlG01_get( &
            flx_mss_hrz_slt_ttl_wbn, & ! O [kg m-1 s-1] Saltation flux
            flx_mss_hrz_slt_ttl_lut, & ! I [kg m-1 s-1] Look up table for saltation flux
            sfc_frc_bln, & ! I [frc] Fraction of soil blends
            wnd_frc_slt) ! I [m s-1] Wind friction speed 
#else  /* !AlG01 */
       ! Horizontal streamwise mass flux for old "bulk" formulation (checked ok for CTM3)
       call flx_mss_hrz_slt_ttl_Whi79_get( &
            dns_mdp, & ! I [kg m-3] Midlayer density
            flg_mbl, & ! I [flg] Mobilization candidate flag
            flx_mss_hrz_slt_ttl_wbn, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
            wnd_frc_slt, & ! I [m s-1] Saltating friction velocity
            wnd_frc_thr_slt) ! I [m s-1] Threshold friction speed for saltation
#endif /* !AlG01 */


       ! Apply land surface and vegetation limitations and global tuning factor
       flx_mss_hrz_slt_ttl_wbn = flx_mss_hrz_slt_ttl_wbn & ! [kg m-2 s-1]
            * lnd_frc_mbl    & ! [frc] Bare ground fraction
            * mbl_bsn_fct    & ! [frc] Erodibility factor
            * flx_mss_fdg_fct  ! [frc] Global mass flux tuning factor (empirical)


#ifdef AlG01
       !// (checked ok for CTM3)
       call mss_frc_src_AlG01_get( &
            wnd_frc, & ! I [m s-1] wind friction velocity (Should we use wnd_frc_slt here ??)
            mss_frc_src_lut, & ! I [frc] mass fraction source look up table
            mss_frc_src_wbn, & ! O [frc] mass fraction of source (looked up)
            sfc_frc_bln) ! I [frc] weighting of soil types 

       !write(6,*)'ALG01 wnd_frc ',wnd_frc
       !write(6,*)'AlG01 mss_frc_src',mss_frc_src_wbn


       !// (checked ok for CTM3)
       call ovr_src_snk_mss_AlG01_get( &
            mss_frc_src_wbn, & ! I [frc] mass fraction of source modes (looked up)
            ovr_src_snk_mss_wbn, & ! O [frc] mass overlap source and sink
            mss_frc_trn_dst_src_wbn, & ! O [frc] fraction of transported dust in each transport bin
            ovr_src_snk_mss_ttl_wbn) ! O [frc] total fraction of produced dust transported

       
       !// (checked ok for CTM3)
       call flx_mss_vrt_dst_ttl_AlG01_get( &
            dst_slt_flx_rat_ttl_lut, & ! I [frc] look up table for alpha (function of wind and soil type)
            dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
            flx_mss_hrz_slt_ttl_wbn, & ! I [kg m s-1] horizontal dust flux
            flx_mss_vrt_dst_ttl_wbn, & ! O [kg m-2 s-1] total vertical dust flux
            wnd_frc, & ! I [m s-1] wind friction speed (fxm: should be wnd_frc_slt?)
            sfc_frc_bln) ! I [frc] fraction of soil type
#else
       
       ! Vertical dust mass flux (checked ok for CTM3)
       call flx_mss_vrt_dst_ttl_MaB95_get( &
            dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
            flg_mbl, & ! I [flg] Mobilization candidate flag
            flx_mss_hrz_slt_ttl_wbn, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
            flx_mss_vrt_dst_ttl_wbn, & ! O [kg m-2 s-1] Total vertical mass flux of dust
            mss_frc_cly)         ! I [frc] Mass fraction clay 
       
#endif /* !AlG01 */


       !write(6,*) 'Bin horizontal flux',flx_mss_hrz_slt_ttl_wbn(i)
       !write(6,*) 'Bin vertical flux  ',flx_mss_vrt_dst_ttl_wbn(i)
       ! Add up total horizontal flux (update no matter what)
       ! NB: Some horizontal fluxes do not give vertical fluxes!
       ! Weight horizontal flux by wind speed PDF
       flx_mss_hrz_slt_ttl_wbn = &
            flx_mss_hrz_slt_ttl_wbn*wnd_mdp_wgt(wnd_mdp_idx)
       ! Add up the horizontal saltation soil flux
       flx_mss_hrz_slt_ttl = &
            flx_mss_hrz_slt_ttl+flx_mss_hrz_slt_ttl_wbn
       ! Modify microphysical variables iff flx_mss_vrt_dst_ttl > 0.0
       if (flx_mss_vrt_dst_ttl_wbn > 0.0_r8) then
          ! Weight vertical flux by wind speed PDF
          flx_mss_vrt_dst_ttl_wbn = &
                  flx_mss_vrt_dst_ttl_wbn*wnd_mdp_wgt(wnd_mdp_idx)
          ! Add up mass fraction of source modes
          mss_frc_src_add(:)=mss_frc_src_add(:)+ &    
               mss_frc_src_wbn(:)*flx_mss_vrt_dst_ttl_wbn
          ! Add up transported mass fraction of source
          mss_frc_trn_dst_src_add(:)=mss_frc_trn_dst_src_add(:)+ &
               mss_frc_trn_dst_src_wbn(:)*flx_mss_vrt_dst_ttl_wbn
          ! Add up mass overlap fraction
          ovr_src_snk_mss_add(:,:)= ovr_src_snk_mss_add(:,:)+ &
               ovr_src_snk_mss_wbn(:,:)*flx_mss_vrt_dst_ttl_wbn
          ! Add up total mass overlap fraction
          ovr_src_snk_mss_ttl_add = ovr_src_snk_mss_ttl_add+ &
               ovr_src_snk_mss_ttl_wbn*flx_mss_vrt_dst_ttl_wbn
          flx_mss_vrt_dst_ttl = flx_mss_vrt_dst_ttl+flx_mss_vrt_dst_ttl_wbn
       endif ! endif flux is zero
    enddo ! end loop over pdf


    ! Normalize total mass fractions
    ! write(6,*)'diagnostics for sum of wind speeds'
    ! Normalize sums by total vertical dust flux
    if (flx_mss_vrt_dst_ttl > 0.0_r8) then
       mss_frc_src_add(:) = mss_frc_src_add(:)/flx_mss_vrt_dst_ttl
       mss_frc_trn_dst_src_add(:) = mss_frc_trn_dst_src_add(:)/flx_mss_vrt_dst_ttl
       ovr_src_snk_mss_add(:,:) = ovr_src_snk_mss_add(:,:)/flx_mss_vrt_dst_ttl 
       ovr_src_snk_mss_ttl_add = ovr_src_snk_mss_ttl_add/flx_mss_vrt_dst_ttl
       ! Diagnose final alpha as (vertical flux)/(horizontal flux)
       dst_slt_flx_rat_ttl = flx_mss_vrt_dst_ttl/flx_mss_hrz_slt_ttl
    else ! else flux==0.0
       ! Must reset mass overlap factors when flux is zero
       ! Only matters for MaB95 formulation, not for AlG01
       mss_frc_src_add(:) = mss_frc_src_wbn(:)
       ovr_src_snk_mss_add(:,:) = ovr_src_snk_mss_wbn(:,:)
       ovr_src_snk_mss_ttl_add = ovr_src_snk_mss_ttl_wbn
       mss_frc_trn_dst_src_add(:) = mss_frc_trn_dst_src_wbn(:)
    endif ! endif flux==0.0
       

    ! Partition vertical mass flux into transport bins (checked ok for CTM3)
    call flx_mss_vrt_dst_prt( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         flx_mss_vrt_dst,     & ! O [kg m-2 s-1] Vertical mass flux of dust
         flx_mss_vrt_dst_ttl, & ! I [kg m-2 s-1] Total vertical mass flux of dust
         ovr_src_snk_mss_add)   ! I [frc] The overlap fractions added up 


    ! Mask dust mass flux by tracer mass fraction at source (checked ok for CTM3)
    flg_CaCO3=.false.         ! [flg] Activate CaCO3 tracer
    if (flg_CaCO3) then
       call flx_mss_CaCO3_msk( &
            dmt_vwr,          & ! I [m] Mass weighted diameter resolved 
            flg_mbl,          & ! I [flg] Mobilization candidate flag
            flx_mss_vrt_dst,  & ! I/O [kg m-2 s-1] Vertical mass flux of dust
            mss_frc_CaCO3,    & ! I [frc] Mass fraction of CaCO3
            mss_frc_cly,      & ! I [frc] Mass fraction of clay
            mss_frc_snd)      ! I [frc] Mass fraction of sand
    endif                     ! endif flg_CaCO3



    !// Fluxes are known, so adjust mixing ratios
    if (flg_mbl) then
       do m=1,dst_nbr
          !// Update dust mass [kg]
          m_dst(m) = m_dst(m) + flx_mss_vrt_dst(m) * area * tm_dlt

          !// Not necessary to calculate q-tendencies anymore,
          !// we only need mass path tendencies.
          !q_dst_tnd_mbl(m) = & ! [kg kg-1 s-1]
          !     flx_mss_vrt_dst(m)/mpl_air
          !q_dst(m) = q_dst(m)+ & ! [kg kg-1]
          !     tm_dlt * q_dst_tnd_mbl(m)
!#ifdef DST_DBG
          !if (q_dst(m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,es8.1,a)')  &
          !     'dst_mbl: lat = ',lat_idx,' q_dst(',lon_idx,',',m,') = ',q_dst(m),' kg kg-1'
!#endif /*!DST_DBG */

       end do !// end loop over cst

#ifdef DST_DBG
       if (flx_mss_vrt_dst_ttl > flx_mss_mxm) write(6,'(a,i2,a,i3,a,es8.1,a)')  &
            'dst_mbl: lat = ',lat_idx,' flx_mss_vrt_dst_ttl(',lon_idx,') = ',flx_mss_vrt_dst_ttl,' kg m-2 s-1'
#endif /*!DST_DBG */
    endif                  ! endif flg_mbl


    ! Integrate over all size categories for diagnostic output
    !call dst_add_nbr(q_dst_tnd_mbl,q_dst_tnd_mbl_ttl)

    ! Recompute flx_mss_vrt_dst_ttl as sum of transported dust
    call dst_add_nbr(flx_mss_vrt_dst, & ! I [kg m-2 s-1] Vertical mass flux of dust
         flx_mss_vrt_dst_ttl) ! O [kg m-2 s-1] Total vertical mass flux of dust

    ! Jump to here when no points are mobilization candidates
737 continue

    !// Save production term of all dst_nbr species  [kg]
    do m = 1, dst_nbr
       production(m) = flx_mss_vrt_dst(m) * area * tm_dlt
    end do

    !// Put to history fields
    call outfld_1('DSTSFMBL',flx_mss_vrt_dst_ttl,lon_idx,lat_idx,obuf)
!    call outfld_1('BSN_FCT',mbl_bsn_fct,lon_idx,lat_idx,obuf)
!    call outfld_1('DNS_MDP',dns_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('SNW_FRC',snw_frc,lon_idx,lat_idx,obuf)
!    call outfld_1('FLX_SWA',flx_SW_abs_sfc,lon_idx,lat_idx,obuf)
!    call outfld_1('FLX_LWD',flx_LW_dwn_sfc,lon_idx,lat_idx,obuf)
!    call outfld_1('FRC_WET',frc_thr_ncr_wtr,lon_idx,lat_idx,obuf)
!    call outfld_1('HGT_MDP',hgt_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('HGT_ZPD',hgt_zpd,lon_idx,lat_idx,obuf)
!    call outfld_1('LND_DRY',lnd_frc_dry,lon_idx,lat_idx,obuf)
!    call outfld_1('MNO_LNG',mno_lng,lon_idx,lat_idx,obuf)
!    call outfld_1('ORO',oro,lon_idx,lat_idx,obuf)
!    call outfld_1('PRS_DLT',prs_dlt,lon_idx,lat_idx,obuf)
!    call outfld_1('PRS_MDP',prs_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('PRS_SFC',prs_sfc,lon_idx,lat_idx,obuf)
!    call outfld_1('Q_H2O',q_h2o_vpr,lon_idx,lat_idx,obuf)
!    call outfld_1('RGH_MMN',rgh_mmn_mbl,lon_idx,lat_idx,obuf)
!    call outfld_1('SNW_FRC',snw_frc,lon_idx,lat_idx,obuf)
!    call outfld_1('TPT_GND',tpt_gnd,lon_idx,lat_idx,obuf)
!    call outfld_1('TPT_MDP',tpt_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('TPT_PTN',tpt_ptn_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('WND_FRC_MBL',wnd_frc,lon_idx,lat_idx,obuf)
!    call outfld_1('WND_RFR_MBL',wnd_rfr,lon_idx,lat_idx,obuf)
!    call outfld_1('WNDZNLMDP',wnd_znl_mdp,lon_idx,lat_idx,obuf)
!    call outfld_1('WNDMRDMDP',wnd_mrd_mdp,lon_idx,lat_idx,obuf)

     !// Could be interesting; need to check dead_history.f90 to
     !// see which are included.
!    call outfld_1('WND_FRCS',wnd_frc_slt,lon_idx,lat_idx,obuf)
!    call outfld_1('WND_FRCT',wnd_frc_thr_slt,lon_idx,lat_idx,obuf)
!    call outfld_1('WND_RFRT',wnd_rfr_thr_slt,lon_idx,lat_idx,obuf)
!    call outfld_1('GWC_SFC',gwc_sfc,lon_idx,lat_idx,obuf)
!    call outfld_1('VAI_DST',vai_dst,lon_idx,lat_idx,obuf)
!    call outfld_1('VWC_SFC',vwc_sfc,lon_idx,lat_idx,obuf)


    !// ------------------------------------------------------------------
  end subroutine dst_mbl                       ! end dst_mbl()
  !// ------------------------------------------------------------------
  
 
end module dstmbl

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstmbl.F90,v 1.3 2007/04/11 08:35:28 alfgr Exp $

! Purpose: dstmbl.F90 controls mineral dust mobilization processes
! Usage:
! use dstmbl ! [mdl] Mobilization driver

! Requires dst.h for DST_MSS_BDG
#include <dst.h> /* Dust preprocessor tokens */

module dstmbl ! [mdl] Mobilization driver
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_mbl ! [sbr] Driver for aerosol mobilization

contains
  
  subroutine dst_mbl(lchnk,ncol,obuf, &
       doy,                 & ! I [day] Day of year [1.0..366.0)
       hgt_mdp,             & ! I [m] Midpoint height above surface
       lat_idx,             & ! I [idx] Model latitude index
       lat_rdn,             & ! I [rdn] Latitude
       oro,                 & ! I [frc] Orography
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       prs_sfc,             & ! I [Pa] Surface pressure
       q_H2O_vpr,           & ! I [kg kg-1] Water vapor mixing ratio
       q_dst,               & ! I/O [kg kg-1] Dust mixing ratio
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp,             & ! I [K] Temperature
       tpt_ptn_mdp,         & ! I [K] Midlayer local potential temperature
       wnd_mrd_mdp,         & ! I [m s-1] Meridional wind component
       wnd_znl_mdp)         ! I [m s-1] Zonal wind component
    ! Purpose: Driver for aerosol mobilization
    ! dst_mbl() is called by CCM:physics/tphysac(), MATCH:src/physlic()
    ! NB: dst_mbl() is designed to require only single layer surface fields 
    ! Eliminating multi-layer fields allows easier implementation in LSM
    use blmutl,only:wnd_rfr_get,snw_frc_get ! [mdl] Boundary layer meteorology driver
    use dstaer,only:dmt_vwr,dns_aer,mss_frc_src,mss_frc_trn_dst_src,ovr_src_snk_mss,ovr_src_snk_mss_ttl ! [mdl] Aerosol microphysical properties
    use dstbdg,only:bdg_aal,bdg_gam_wet_2d ! [mdl] Mass budget diagnostics
    use dstblm,only:blm_mbl,cnd_trm_soi_get,hyd_prp_get,tpt_frz_pnt,trn_fsh_vpr_soi_atm_get ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use dstcst ! [mdl] Physical constants for dust routines
    use dstdbg ! [mdl] Debugging information for dust model
    use dstgrd,only:dst_nbr,bln_nbr,dst_src_nbr ! [mdl] Dust grid sizes
    use dstmblutl ! [mdl] Mobilization utilities
    use dstmssutl,only:dst_add_lon ! [mdl] Mass budget utilities
    use dstnm ! [mdl] Nomenclature for outfld()
    use dstsfc,only:sfc_typ_get,soi_txt_get,vwc_sfc_get,tpt_gnd_soi_get,flx_rad_sfc_get ! [mdl] 2-D surface fields on PLON x PLAT grid
    use dstsltsbl  ! [mdl] Saltation sandblasting physics
    use dsttvbds,only:dst_tvbds_get ! [mdl] Time-varying boundary data sets
#ifndef BXM
    use histout,only:outfld ! [mdl] History/archiving
#endif /* BXM */
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    use precision ! [mdl] Precision r8, i8, ...
    use sng_mdl,only:ftn_strnul ! [mdl] String manipulation
    use wbl_mdl,only:wbl_wnd,wnd_mdp_bin_get ! [mdl] Weibull wind speed distribution
    implicit none
    ! Parameters
    real(r8),parameter::flx_mss_fdg_fct=7.0e-4*1500./1454. ! [frc] Global mass flux tuning factor (a posteriori)
    real(r8),parameter::hgt_rfr=10.0 ! [m] Reference height for mobilization processes
    real(r8),parameter::hgt_zpd_mbl=0.0 ! [m] Zero plane displacement for erodible surfaces
#ifndef AlG01
    real(r8),parameter::rgh_mmn_mbl=100.0e-6 ! [m] Roughness length momentum for erodible surfaces MaB95 p. 16420, GMB98 p. 6205
#else
    real(r8),parameter::rgh_mmn_mbl=3000.0e-6 ! [m] Roughness length momentum for erodible surface (AlG01 formulation)
#endif
    ! fxm: rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
    real(r8),parameter::rgh_mmn_smt=33.3e-6 ! [m] Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
    real(r8),parameter::wnd_min_mbl=1.0 ! [m s-1] Minimum windspeed used for mobilization 
    ! Input
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::doy ! I [day] Day of year [1.0..366.0)
    real(r8),intent(in)::hgt_mdp(plond) ! I [m] Midlayer height above surface
    real(r8),intent(in)::lat_rdn ! I [rdn] Latitude
    real(r8),intent(in)::obuf(*) ! I [ptr] Output buffer
    real(r8),intent(in)::oro(plond) ! I [frc] Orography
    real(r8),intent(in)::prs_dlt(plond) ! I [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond) ! I [Pa] Pressure
    real(r8),intent(in)::prs_sfc(plond) ! I [Pa] Surface pressure
    real(r8),intent(in)::q_H2O_vpr(plond) ! I [kg kg-1] Water vapor mixing ratio
    real(r8),intent(in)::snw_hgt_lqd(plond) ! I [m] Equivalent liquid water snow depth
    real(r8),intent(in)::tm_adj ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
    real(r8),intent(in)::tpt_mdp(plond) ! I [K] Temperature
    real(r8),intent(in)::tpt_ptn_mdp(plond) ! I [K] Midlayer local potential temperature
    real(r8),intent(in)::wnd_mrd_mdp(plond) ! I [m s-1] Meridional wind component
    real(r8),intent(in)::wnd_znl_mdp(plond) ! I [m s-1] Zonal wind component
    ! Input/Output
    real(r8),intent(inout)::q_dst(plond,dst_nbr) ! I/O [kg kg-1] Dust mixing ratio
    ! Local Output
    real(r8) dst_slt_flx_rat_ttl(plond) ! [m-1] Ratio of vertical dust flux to streamwise mass flux
    real(r8) flx_mss_hrz_slt_ttl(plond) ! [kg m-1 s-1] Vertically integrated streamwise mass flux
    real(r8) flx_mss_vrt_dst(plond,dst_nbr) ! [kg m-2 s-1] Vertical mass flux of dust
    real(r8) flx_mss_vrt_dst_ttl(plond) ! [kg m-2 s-1] Total vertical mass flux of dust
    real(r8) frc_thr_ncr_drg(plond) ! [frc] Threshold friction velocity increase from roughness
    real(r8) frc_thr_ncr_wtr(plond) ! [frc] Threshold friction velocity increase from moisture
    real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement
    real(r8) lnd_frc_mbl(plond)   ! [frc] Bare ground fraction
    real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
    real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
    real(r8) trn_fsh_vpr_soi_atm(plond) ! [frc] Transfer efficiency of vapor from soil to atmosphere
    real(r8) wnd_frc(plond)       ! [m s-1] Friction velocity
    real(r8) wnd_frc_slt(plond)   ! [m s-1] Saltating friction velocity
    real(r8) wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
    real(r8) wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
    ! Local
    logical flg_CaCO3         ! [flg] Activate CaCO3 tracer
    logical flg_mbl(plond)    ! [flg] Mobilization candidates
    character(80) fl_out      ! [sng] Name of netCDF output file
    integer i                 ! [idx] Counting index for longitude
    integer m                 ! [idx] Counting index for species
    integer mbl_nbr           ! [nbr] Number of mobilization candidates
    integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
    real(r8) cnd_trm_soi(plond)   ! [W m-1 K-1] Soil thermal conductivity
    real(r8) dns_mdp(plond)       ! [kg m-3] Midlayer density
    real(r8) flx_LW_dwn_sfc(plond) ! [W m-2] Longwave downwelling flux at surface
    real(r8) flx_SW_abs_sfc(plond) ! [W m-2] Solar flux absorbed by ground
    real(r8) lnd_frc_dry(plond)   ! [frc] Dry land fraction
    real(r8) mbl_bsn_fct(plond)   ! [frc] Erodibility factor
    real(r8) mss_frc_CaCO3(plond) ! [frc] Mass fraction of CaCO3
    real(r8) lvl_dlt(plond)       ! [m] Soil layer thickness
    real(r8) mpl_air(plond)       ! [kg m-2] Air mass path in layer
    real(r8) mss_frc_cly(plond)   ! [frc] Mass fraction of clay
    real(r8) mss_frc_snd(plond)   ! [frc] Mass fraction of sand
    real(r8) sfc_frc_bln(plond,bln_nbr) ! [frc] fraction of each soil type
    real(r8) tm_dlt               ! [s] Mobilization timestep
    real(r8) tpt_gnd(plond)       ! [K] Ground temperature
    real(r8) tpt_soi(plond)       ! [K] Soil temperature
    real(r8) tpt_soi_frz          ! [K] Temperature of frozen soil
    real(r8) tpt_vrt_mdp          ! [K] Midlayer virtual temperature
    real(r8) vai_dst(plond)       ! [m2 m-2] Vegetation area index, one-sided
    real(r8) vwc_dry(plond)       ! [m3 m-3] Dry volumetric water content (no E-T)
    real(r8) vwc_opt(plond)       ! [m3 m-3] E-T optimal volumetric water content 
    real(r8) vwc_sat(plond)       ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    real(r8) vwc_sfc(plond)       ! [m3 m-3] Volumetric water content
    real(r8) gwc_sfc(plond)       ! [kg kg-1] Gravimetric water content
    ! ++alfgr variables needed for Weibull code
    integer,parameter::wnd_mdp_nbr=5 ! [nbr] Number of discrete wind speeds
    integer::wnd_mdp_idx ! [idx] Counter for wind speed discretization
    real(r8),parameter::wnd_frc_rsl=0.95_r8 ! Fraction of wind PDF to resolve
    real(r8)::wnd_mdp_wgt(plond,wnd_mdp_nbr) ! Weighting for winds
    real(r8)::wnd_mdp_bin(plond) ! Wind in the actual wind bin
    real(r8)::wnd_mdp_min(plond) ! Mininum wind speed we are interested in
    real(r8)::wnd_mdp_max(plond) ! Maximum wind speed we are interested in
    real(r8)::ovr_src_snk_mss_wbn(plond,dst_src_nbr,dst_nbr) ! Overlap mass factors for wbin
    real(r8)::ovr_src_snk_mss_add(plond,dst_src_nbr,dst_nbr) ! Overlap mass factors (added up)
    real(r8)::ovr_src_snk_mss_ttl_wbn(plond) ! Total overlap between source and sink from wbin
    real(r8)::ovr_src_snk_mss_ttl_add(plond) ! Total overlap between source and sink (added up)
    real(r8)::flx_mss_hrz_slt_ttl_wbn(plond) ! Horizontal flux due to a wind bin    
    real(r8)::flx_mss_vrt_dst_ttl_wbn(plond) ! Vertical flux due to a wind bin
    real(r8)::mss_frc_src_wbn(plond,dst_src_nbr) ! Mass fraction source due to this wind bin
    real(r8)::mss_frc_src_add(plond,dst_src_nbr) ! Mass fraction source (added up)
    real(r8)::mss_frc_trn_dst_src_wbn(plond,dst_nbr) ! Mass fraction of source transported from wbin
    real(r8)::mss_frc_trn_dst_src_add(plond,dst_nbr) ! Mass fraction of source transported (added up)
    ! --alfgr variables needed for Weibull winds

    ! GCM diagnostics
    real(r8) q_dst_tnd_mbl(plond,dst_nbr) ! [kg kg-1 s-1] Dust tendency due to gravitational settling
    real(r8) q_dst_tnd_mbl_ttl(plond) ! [kg kg-1 s-1] Total dust tendency due to gravitational settling
    
    ! Main Code
    ! Timesplit if desired
    tm_dlt=tm_adj ! [s] (default CCM: 2*dt, MATCH: dt)
    ! Assume water in soil freezes at 0 C
    tpt_soi_frz=tpt_frz_pnt ! [K] Temperature of frozen soil
    
    ! Initialize output fluxes and tendencies
    q_dst_tnd_mbl(:,:)=0.0_r8 ! [kg kg-1 s-1]
    q_dst_tnd_mbl_ttl(:)=0.0_r8 ! [kg kg-1 s-1]
    flx_mss_vrt_dst(:,:)=0.0_r8 ! [kg m-2 s-1]
    flx_mss_vrt_dst_ttl(:)=0.0_r8 ! [kg m-2 s-1]
    frc_thr_ncr_wtr(:)=0.0_r8 ! [frc]
    wnd_rfr(:)=0.0_r8 ! [m s-1]
    wnd_frc(:)=0.0_r8 ! [m s-1]
    wnd_frc_slt(:)=0.0_r8 ! [m s-1]
    wnd_frc_thr_slt(:)=0.0_r8 ! [m s-1]
    wnd_rfr_thr_slt(:)=0.0_r8 ! [m s-1]
    hgt_zpd(:)=hgt_zpd_mbl ! [m]

    ! Wind PDF initialization
    ! Initialize mass fraction of source
    mss_frc_src_wbn(:,:)=mss_frc_src(:,:)
    mss_frc_src_add(:,:)=0.0_r8 ! Must be zero because it is updated and weighted
    ! Initialize overlap of source and sink
    ovr_src_snk_mss_wbn(:,:,:)=ovr_src_snk_mss(:,:,:)
    ovr_src_snk_mss_add(:,:,:)=0.0_r8 ! Must start as zero because it is updated and weighted
    ! Initialize total overlap src/sink
    ovr_src_snk_mss_ttl_wbn(:)=ovr_src_snk_mss_ttl(:)
    ovr_src_snk_mss_ttl_add(:)=0.0_r8 ! Must start as zero because it is updated and weighted
    mss_frc_trn_dst_src_wbn(:,:)=mss_frc_trn_dst_src(:,:)
    mss_frc_trn_dst_src_add(:,:)=0.0_r8 ! Must start as zero because it is updated and weighted      
    
    ! Compute required derived fields
    do i=1,plon
       ! Stop occasional haywire model runs 
       if(tpt_mdp(i) > 350.0_r8) stop 'ERROR: dst_mbl() reports tpt_mdp(i) > 350.0'
       ! Midlayer virtual temperature
       tpt_vrt_mdp=tpt_mdp(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
       ! Density at center of gridbox
       dns_mdp(i)=prs_mdp(i)/(tpt_vrt_mdp*gas_cst_dry_air) ! [kg m-3]
       ! Approximate surface virtual temperature (uses midlayer moisture)
       ! tpt_vrt_sfc=tpt_sfc(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
       ! Surface density
       ! dns_sfc(i)=prs_sfc(i)/(tpt_vrt_sfc*gas_cst_dry_air) ! [kg m-3]
       ! Mass of air currently in gridbox
       mpl_air(i)=prs_dlt(i)*grv_sfc_rcp ! [kg m-2]
       ! Mean surface layer horizontal wind speed 
       wnd_mdp(i)=sqrt(wnd_znl_mdp(i)*wnd_znl_mdp(i)+wnd_mrd_mdp(i)*wnd_mrd_mdp(i))
       
       !write(6,*)'WIND SPEED AT MIDPOINT',wnd_mdp(i)
    end do                    ! end loop over lon
    
    ! Gather input variables from common blocks
    ! Surface type
    call sfc_typ_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         sfc_typ)             ! O [idx] LSM surface type (0..28)
    ! Soil texture and dry land fraction
    call soi_txt_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         lnd_frc_dry,         & ! O [frc] Dry land fraction
         mbl_bsn_fct,         & ! O [frc] Erodibility factor
         mss_frc_CaCO3,       & ! O [frc] Mass fraction of CaCO3
         mss_frc_cly,         & ! O [frc] Mass fraction of clay
         mss_frc_snd,         & ! O [frc] Mass fraction of sand
         sfc_frc_bln)         ! O [frc] Fraction of four soil types
    ! Soil moisture
    call vwc_sfc_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         vwc_sfc)             ! O [m3 m-3] Volumetric water content
    ! Surface and soil temperature
    call tpt_gnd_soi_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         tpt_gnd,             & ! O [K] Ground temperature
         tpt_soi)             ! O [K] Soil temperature
    ! Radiative fluxes at surface
    call flx_rad_sfc_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         flx_LW_dwn_sfc,      & ! O [W m-2] Longwave downwelling flux at surface
         flx_SW_abs_sfc)      ! O [W m-2] Solar flux absorbed by ground
    ! Time varying boundary datasets: Vegetation area index
    call dst_tvbds_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         vai_dst)             ! O [m2 m-2] Vegetation area index, one-sided
    
    ! All variables from common blocks have been assembled
    ! Use these to compute derived fields
    ! Fraction of surface covered by snow
    call snw_frc_get( &
         snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
         snw_frc)             ! O [frc] Fraction of surface covered by snow
    
    ! Fraction of each gridcell suitable for dust mobilization
    call lnd_frc_mbl_get( &
         doy,                 & ! I [day] Day of year [1.0..366.0)
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
         vai_dst)             ! I [m2 m-2] Vegetation area index, one-sided
    
    ! Much ado about nothing
    if (mbl_nbr == 0) goto 737
    
    ! Hydrologic properties
    call hyd_prp_get(         & ! NB: These properties are time-invariant
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         mss_frc_cly,         & ! I [frc] Mass fraction clay 
         mss_frc_snd,         & ! I [frc] Mass fraction sand
         vwc_dry,             & ! O [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt,             & ! O [m3 m-3] E-T optimal volumetric water content
         vwc_sat)             ! O [m3 m-3] Saturated volumetric water content
    
    ! Soil thermal conductivity and layer thickness
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
         vwc_sfc)             ! I [m3 m-3] Volumetric water content
    
    ! Transfer efficiency of vapor from soil to atmosphere
    call trn_fsh_vpr_soi_atm_get( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         tpt_soi,             & ! I [K] Soil temperature
         tpt_soi_frz,         & ! I [K] Temperature of frozen soil
         trn_fsh_vpr_soi_atm, & ! O [frc] Transfer efficiency of vapor from soil to atmosphere
         vwc_dry,             & ! I [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt,             & ! I [m3 m-3] E-T optimal volumetric water content
         vwc_sfc)             ! I [m3 m-3] Volumetric water content
    
    ! Surface exchange properties over erodible surfaces
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
         wnd_frc)             ! O [m s-1] Surface friction velocity
    
    ! Interpolate midlayer wind speed to 10 m
    call wnd_rfr_get( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         hgt_mdp,             & ! I [m] Midpoint height above surface
         hgt_rfr,             & ! I [m] Reference height for mobilization processes
         hgt_zpd,             & ! I [m] Zero plane displacement
         mno_lng,             & ! I [m] Monin-Obukhov length
         wnd_frc,             & ! I [m s-1] Surface friction velocity
         wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
         wnd_min_mbl,         & ! I [m s-1] Minimum windspeed
         wnd_rfr)             ! O [m s-1] Wind speed at reference height

    ! Get the weighting factors and winds we want to calculate for
    call wbl_wnd( &
         hgt_mdp, & ! I [m] Height of layer midpoint
         hgt_rfr, & ! I [m] Reference height
         wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
         wnd_mdp, & ! I [m s-1] Wind speed at midpoint
         wnd_mdp_min, & ! O [m s-1] Min wind speed at midpoint
         wnd_mdp_max, & ! O [m s-1] Max wind speed at midpoint 
         wnd_mdp_nbr, & ! I [nbr] Number of wind speeds to calculate 
         wnd_mdp_wgt, & ! O [frc] Wind speed bin weight
         wnd_rfr) ! I [m s-1] Wind speed at reference height

    ! Now, do the loop on wnd_mdp
    do wnd_mdp_idx=1,wnd_mdp_nbr

       if(wnd_mdp_nbr > 1) then  ! More than one wind bin

          ! Get mean wind speed in this bin
          call wnd_mdp_bin_get( &
               wnd_mdp_bin, & ! O [m s-1] wind speed in bin
               wnd_mdp_idx, & ! I [nbr] The bin we are interested in
               wnd_mdp_max, & ! I [m s-1] largest wind speed which we calculate
               wnd_mdp_min, & ! I [m s-1] smallest wind speed which we calculate
               wnd_mdp_nbr) ! I [nbr] Weibull wind number
    
          !write(6,*)'  '
          !write(6,*)'WIND SPEED',wnd_mdp_bin,'@ height',hgt_mdp,'weight',wnd_mdp_wgt(:,wnd_mdp_idx)
          
          ! Get wind friction speed for this wind speed
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
               wnd_frc)             ! O [m s-1] Surface friction velocity

       else ! Only one wind speed
          ! Set weights to unity then continue as before using mean quantities
          wnd_mdp_wgt(:,:)=1.0_r8 ! [frc] Wind speed bin weight
       endif ! endif 

       ! Factor by which surface roughness increases threshold friction velocity 
       call frc_thr_ncr_drg_get( &
            frc_thr_ncr_drg, & ! O [frc] Factor by which surface roughness increases threshold friction velocity
            rgh_mmn_mbl, & ! I [m] Roughness length momentum for erodible surfaces
            rgh_mmn_smt) ! I [m] Smooth roughness length
    
       ! Convert volumetric water content to gravimetric water content
       call vwc2gwc( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            gwc_sfc, & ! O [kg kg-1] Gravimetric water content
            vwc_sat, & ! I [m3 m-3] Saturated volumetric water content (sand-dependent)
            vwc_sfc) ! I [m3 m-3] Volumetric water content

       ! Factor by which soil moisture increases threshold friction velocity 
       call frc_thr_ncr_wtr_get( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            frc_thr_ncr_wtr, & ! O [frc] Factor by which moisture increases threshold friction velocity
            mss_frc_cly, & ! I [frc] Mass fraction of clay
            gwc_sfc)             ! I [kg kg-1] Gravimetric water content
    
       ! Threshold friction velocity for saltation over dry, bare, smooth ground
       ! fxm: Use surface density not midlayer density
       call wnd_frc_thr_slt_get( &
            dmt_vwr, & ! I [m] Mass weighted diameter resolved 
            dns_aer, & ! I [kg m-3] Particle density
            dns_mdp, & ! I [kg m-3] Midlayer density
            wnd_frc_thr_slt)     ! O [m s-1] Threshold friction velocity for saltation
    
       ! Adjust threshold friction velocity to account for moisture and roughness
       do i=1,plon
          wnd_frc_thr_slt(i)= & ! [m s-1] Threshold friction velocity for saltation
               wnd_frc_thr_slt(i)* & ! [m s-1] Threshold for dry, flat ground
               frc_thr_ncr_wtr(i)* & ! [frc] Adjustment for moisture
               frc_thr_ncr_drg(i) ! [frc] Adjustment for roughness
       end do                    ! end loop over lon
    
       ! Threshold saltation wind speed
       do i=1,plon
          if (flg_mbl(i)) then
             wnd_rfr_thr_slt(i)= & ! [m s-1] Threshold 10 m wind speed for saltation
                  wnd_rfr(i)*wnd_frc_thr_slt(i)/wnd_frc(i)
          endif                  ! endif flg_mbl
       end do                    ! end loop over lon
    
       ! Saltation increases friction speed by roughening surface
       call wnd_frc_slt_get( &
            flg_mbl, & ! I [flg] Mobilization candidate flag
            wnd_frc, & ! I [m s-1] Surface friction velocity
            wnd_frc_slt, & ! O [m s-1] Saltating friction velocity
            wnd_rfr, & ! I [m s-1] Wind speed at reference height
            wnd_rfr_thr_slt) ! I [m s-1] Threshold 10 m wind speed for saltation

#ifdef AlG01
       ! Horizontal streamwise mass flux consistent with AlG01 formulation
       call flx_mss_hrz_slt_ttl_AlG01_get( &
            flx_mss_hrz_slt_ttl_wbn, & ! O [kg m-1 s-1] Saltation flux
            flx_mss_hrz_slt_ttl_lut, & ! I [kg m-1 s-1] Look up table for saltation flux
            sfc_frc_bln, & ! I [frc] Fraction of soil blends
            wnd_frc_slt) ! I [m s-1] Wind friction speed 
#else  /* !AlG01 */
       ! Horizontal streamwise mass flux for old "bulk" formulation
       call flx_mss_hrz_slt_ttl_Whi79_get( &
            dns_mdp, & ! I [kg m-3] Midlayer density
            flg_mbl, & ! I [flg] Mobilization candidate flag
            flx_mss_hrz_slt_ttl_wbn, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
            wnd_frc_slt, & ! I [m s-1] Saltating friction velocity
            wnd_frc_thr_slt) ! I [m s-1] Threshold friction speed for saltation
#endif /* !AlG01 */
       
       ! Apply land surface and vegetation limitations and global tuning factor
       do i=1,plon
          flx_mss_hrz_slt_ttl_wbn(i)=flx_mss_hrz_slt_ttl_wbn(i) & ! [kg m-2 s-1]
               *lnd_frc_mbl(i) & ! [frc] Bare ground fraction
               *mbl_bsn_fct(i) & ! [frc] Erodibility factor
               *flx_mss_fdg_fct  ! [frc] Global mass flux tuning factor (empirical)
       end do ! end loop over lon

#ifdef AlG01
       call mss_frc_src_AlG01_get( &
            wnd_frc, & ! I [m s-1] wind friction velocity (Should we use wnd_frc_slt here ??)
            mss_frc_src_lut, & ! I [frc] mass fraction source look up table
            mss_frc_src_wbn, & ! O [frc] mass fraction of source (looked up)
            sfc_frc_bln) ! I [frc] weighting of soil types 

       write(6,*)'ALG01 wnd_frc ',wnd_frc
       write(6,*)'AlG01 mss_frc_src',mss_frc_src_wbn

       call ovr_src_snk_mss_AlG01_get( &
            mss_frc_src_wbn, & ! I [frc] mass fraction of source modes (looked up)
            ovr_src_snk_mss_wbn, & ! O [frc] mass overlap source and sink
            mss_frc_trn_dst_src_wbn, & ! O [frc] fraction of transported dust in each transport bin
            ovr_src_snk_mss_ttl_wbn) ! O [frc] total fraction of produced dust transported
       
       call flx_mss_vrt_dst_ttl_AlG01_get( &
            dst_slt_flx_rat_ttl_lut, & ! I [frc] look up table for alpha (function of wind and soil type)
            dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
            flx_mss_hrz_slt_ttl_wbn, & ! I [kg m s-1] horizontal dust flux
            flx_mss_vrt_dst_ttl_wbn, & ! O [kg m-2 s-1] total vertical dust flux
            wnd_frc, & ! I [m s-1] wind friction speed (fxm: should be wnd_frc_slt?)
            sfc_frc_bln) ! I [frc] fraction of soil type
#else
       
       ! Vertical dust mass flux
       call flx_mss_vrt_dst_ttl_MaB95_get( &
            dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
            flg_mbl, & ! I [flg] Mobilization candidate flag
            flx_mss_hrz_slt_ttl_wbn, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
            flx_mss_vrt_dst_ttl_wbn, & ! O [kg m-2 s-1] Total vertical mass flux of dust
            mss_frc_cly)         ! I [frc] Mass fraction clay 
       
#endif /* !AlG01 */

       do i=1,plon
          !write(6,*) 'Bin horizontal flux',flx_mss_hrz_slt_ttl_wbn(i)
          !write(6,*) 'Bin vertical flux  ',flx_mss_vrt_dst_ttl_wbn(i)
          ! Add up total horizontal flux (update no matter what)
          ! NB: Some horizontal fluxes do not give vertical fluxes!
          ! Weight horizontal flux by wind speed PDF
          flx_mss_hrz_slt_ttl_wbn(i)= &
               flx_mss_hrz_slt_ttl_wbn(i)*wnd_mdp_wgt(i,wnd_mdp_idx)
          ! Add up the horizontal saltation soil flux
          flx_mss_hrz_slt_ttl(i)= &
               flx_mss_hrz_slt_ttl(i)+flx_mss_hrz_slt_ttl_wbn(i)
          ! Modify microphysical variables iff flx_mss_vrt_dst_ttl > 0.0
          if (flx_mss_vrt_dst_ttl_wbn(i) > 0.0_r8) then
             ! Weight vertical flux by wind speed PDF
             flx_mss_vrt_dst_ttl_wbn(i)= &
                  flx_mss_vrt_dst_ttl_wbn(i)*wnd_mdp_wgt(i,wnd_mdp_idx)
             ! Add up mass fraction of source modes
             mss_frc_src_add(i,:)=mss_frc_src_add(i,:)+ &    
                  mss_frc_src_wbn(i,:)*flx_mss_vrt_dst_ttl_wbn(i)
             ! Add up transported mass fraction of source
             mss_frc_trn_dst_src_add(i,:)=mss_frc_trn_dst_src_add(i,:)+ &
                  mss_frc_trn_dst_src_wbn(i,:)*flx_mss_vrt_dst_ttl_wbn(i)
             ! Add up mass overlap fraction
             ovr_src_snk_mss_add(i,:,:)= ovr_src_snk_mss_add(i,:,:)+ &
                  ovr_src_snk_mss_wbn(i,:,:)*flx_mss_vrt_dst_ttl_wbn(i)
             ! Add up total mass overlap fraction
             ovr_src_snk_mss_ttl_add(i)=ovr_src_snk_mss_ttl_add(i)+ &
                  ovr_src_snk_mss_ttl_wbn(i)*flx_mss_vrt_dst_ttl_wbn(i)
             flx_mss_vrt_dst_ttl(i)=flx_mss_vrt_dst_ttl(i)+flx_mss_vrt_dst_ttl_wbn(i)
          endif ! endif flux is zero
       enddo ! end loop over longitude
    enddo ! end loop over pdf

    ! Normalize total mass fractions
    do i=1,plon
       ! write(6,*)'diagnostics for sum of wind speeds'
       ! Normalize sums by total vertical dust flux
       if (flx_mss_vrt_dst_ttl(i) > 0.0_r8) then
          mss_frc_src_add(i,:)=mss_frc_src_add(i,:)/flx_mss_vrt_dst_ttl(i) 
          mss_frc_trn_dst_src_add(i,:)=mss_frc_trn_dst_src_add(i,:)/flx_mss_vrt_dst_ttl(i)
          ovr_src_snk_mss_add(i,:,:)=ovr_src_snk_mss_add(i,:,:)/flx_mss_vrt_dst_ttl(i) 
          ovr_src_snk_mss_ttl_add(i)=ovr_src_snk_mss_ttl_add(i)/flx_mss_vrt_dst_ttl(i)
          ! Diagnose final alpha as (vertical flux)/(horizontal flux)
          dst_slt_flx_rat_ttl(i)=flx_mss_vrt_dst_ttl(i)/flx_mss_hrz_slt_ttl(i)
       else ! else flux==0.0
          ! Must reset mass overlap factors when flux is zero
          ! Only matters for MaB95 formulation, not for AlG01
          mss_frc_src_add(i,:)=mss_frc_src_wbn(i,:)
          ovr_src_snk_mss_add(i,:,:)=ovr_src_snk_mss_wbn(i,:,:)
          ovr_src_snk_mss_ttl_add(i)=ovr_src_snk_mss_ttl_wbn(i)
          mss_frc_trn_dst_src_add(i,:)=mss_frc_trn_dst_src_wbn(i,:)
       endif ! endif flux==0.0
       
    enddo ! end loop over longitude

    ! Partition vertical mass flux into transport bins
    call flx_mss_vrt_dst_prt( &
         flg_mbl,             & ! I [flg] Mobilization candidate flag
         flx_mss_vrt_dst,     & ! O [kg m-2 s-1] Vertical mass flux of dust
         flx_mss_vrt_dst_ttl, & ! I [kg m-2 s-1] Total vertical mass flux of dust
         ovr_src_snk_mss_add)   ! I [frc] The overlap fractions added up 

    ! Mask dust mass flux by tracer mass fraction at source
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
    
    ! Fluxes are known, so adjust mixing ratios
    do i=1,plon               ! NB: Inefficient loop order
       if (flg_mbl(i)) then
          do m=1,dst_nbr
             q_dst_tnd_mbl(i,m)= & ! [kg kg-1 s-1]
                  flx_mss_vrt_dst(i,m)/mpl_air(i)
             q_dst(i,m)=q_dst(i,m)+ & ! [kg kg-1]
                  tm_adj*q_dst_tnd_mbl(i,m)
#ifdef DST_DBG
             if (q_dst(i,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,es8.1,a)')  &
                  'dst_mbl: lat = ',lat_idx,' q_dst(',i,',',m,') = ',q_dst(i,m),' kg kg-1'
#endif /*!DST_DBG */
          end do              ! end loop over cst
#ifdef DST_DBG
          if (flx_mss_vrt_dst_ttl(i) > flx_mss_mxm) write(6,'(a,i2,a,i3,a,es8.1,a)')  &
               'dst_mbl: lat = ',lat_idx,' flx_mss_vrt_dst_ttl(',i,') = ',flx_mss_vrt_dst_ttl(i),' kg m-2 s-1'
#endif /*!DST_DBG */
       endif                  ! endif flg_mbl
    end do                    ! end loop over lon
    
    ! Integrate over all size categories for diagnostic output
    call dst_add_lon(q_dst_tnd_mbl,q_dst_tnd_mbl_ttl)
    ! Recompute flx_mss_vrt_dst_ttl as sum of transported dust
    call dst_add_lon(flx_mss_vrt_dst, & ! I [kg m-2 s-1] Vertical mass flux of dust
         flx_mss_vrt_dst_ttl) ! O [kg m-2 s-1] Total vertical mass flux of dust

    ! Jump to here when no points are mobilization candidates
737 continue
    
#if 0
    if (lat_idx == lat_dbg) then
       write (6,'(a,1(a,f9.5,a),19(a,es9.2,a))') &
            'dst: Diagnostics at (lat_dbg,lon_dbg) set in dst_dbg_cmn_ini() :', &
            'doy = ',doy,' day, ', &
            'flx_LW_dwn_sfc = ',flx_LW_dwn_sfc(lon_dbg),' W m-2, ', &
            'flx_SW_abs_sfc = ',flx_SW_abs_sfc(lon_dbg),' W m-2, ', &
            'dst_slt_flx_rat_ttl = ',dst_slt_flx_rat_ttl(lon_dbg),' kg m-1 s-1, ', &
            'flx_mss_hrz_slt_ttl = ',flx_mss_hrz_slt_ttl(lon_dbg),' kg m-1 s-1, ', &
            'frc_thr_ncr_drg = ',frc_thr_ncr_drg(lon_dbg),' frc, ', &
            'frc_thr_ncr_wtr = ',frc_thr_ncr_wtr(lon_dbg),' frc, ', &
            'gwc_sfc = ',gwc_sfc(lon_dbg),' kg kg-1, ', &
            'lnd_frc_mbl = ',lnd_frc_mbl(lon_dbg),' frc, ', &
            'mbl_bsn_fct = ',mbl_bsn_fct(lon_dbg),' frc, ', &
            'mss_frc_CaCO3 = ',mss_frc_CaCO3(lon_dbg),' frc, ', &
            'mno_lng = ',mno_lng(lon_dbg),' m, ', &
            'rgh_mmn = ',rgh_mmn_mbl,' m, ', &
            'tpt_gnd = ',tpt_gnd(lon_dbg),' K, ', &
            'tpt_ptn_mdp = ',tpt_ptn_mdp(lon_dbg),' K, ', &
            'tpt_soi = ',tpt_soi(lon_dbg),' K, ', &
            'vai_dst = ',vai_dst(lon_dbg),' m2 m-2, ', &
            'vwc_sfc = ',vwc_sfc(lon_dbg),' m3 m-3, ', &
            'wnd_frc = ',wnd_frc(lon_dbg),' m s-1, ', &
            'wnd_frc_slt = ',wnd_frc_slt(lon_dbg),' m s-1, ', &
            'wnd_frc_thr_slt = ',wnd_frc_thr_slt(lon_dbg),' m s-1, '
    endif                     ! endif dbg
#endif /* !0 */
    
#ifdef BXM
    ! Update netCDF file
    fl_out='aer.nc' ! [sng] Name of netCDF output file
    call ftn_strnul(fl_out)
    call mbl2nc(             &
         dst_slt_flx_rat_ttl, & ! I [m-1] Ratio of vertical dust flux to streamwise mass flux
         fl_out,              & ! I [sng] Name of netCDF output file
         flx_mss_hrz_slt_ttl, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
         flx_mss_vrt_dst,     & ! I [kg m-2 s-1] Vertical mass flux of dust
         flx_mss_vrt_dst_ttl, & ! I [kg m-2 s-1] Total vertical mass flux of dust
         frc_thr_ncr_drg,     & ! I [frc] Threshold friction velocity increase from roughness
         frc_thr_ncr_wtr,     & ! I [frc] Threshold friction velocity increase from moisture
         gwc_sfc,             & ! I [kg kg-1] Gravimetric water content
         hgt_zpd_mbl,         & ! I [m] Zero plane displacement height
         lat_idx,             & ! I [idx] Model latitude index
         lnd_frc_mbl,         & ! I [frc] Bare ground fraction
         mno_lng,             & ! I [m] Monin-Obukhov length
         rgh_mmn_mbl,         & ! I [m] Roughness length momentum
         snw_frc,             & ! I [frc] Fraction of surface covered by snow
         vai_dst,             & ! I [m2 m-2] Vegetation area index, one-sided
         wnd_frc,             & ! I [m s-1] Friction velocity
         wnd_frc_slt,         & ! I [m s-1] Saltating friction velocity
         wnd_frc_thr_slt,     & ! I [m s-1] Threshold friction velocity for saltation
         wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
         wnd_rfr,             & ! I [m s-1] Wind speed at reference height
         wnd_rfr_thr_slt,     & ! I [m s-1] Threshold 10 m wind speed for saltation
         mss_frc_src_add,         & ! I [frc] Mass fraction of source distribution 
         mss_frc_trn_dst_src_add, & ! I [frc] Fraction of transported dust mass at source
         ovr_src_snk_mss_add,     & ! I [frc] Overlap of src with snk 
         ovr_src_snk_mss_ttl_add  ) ! I [frc] Total transported mass fraction of dust flux 

#endif /*!BXM */
    
#ifdef CCM
    call outfld('BSN_FCT',mbl_bsn_fct,ncol,lchnk)
    call outfld('FRC_WET',frc_thr_ncr_wtr,ncol,lchnk)
    call outfld('GWC_SFC',gwc_sfc,ncol,lchnk)
    call outfld('LND_MBL',lnd_frc_mbl,ncol,lchnk)
    call outfld('TPT_GND',tpt_gnd,ncol,lchnk)
    call outfld('VAI_DST',vai_dst,ncol,lchnk)
    call outfld('VWC_SFC',vwc_sfc,ncol,lchnk)
    call outfld('WND_FRC',wnd_frc,ncol,lchnk)
    call outfld('WND_FRCS',wnd_frc_slt,ncol,lchnk)
    call outfld('WND_FRCT',wnd_frc_thr_slt,ncol,lchnk)
    call outfld('WND_RFR',wnd_rfr,ncol,lchnk)
    call outfld('WND_RFRT',wnd_rfr_thr_slt,ncol,lchnk)
    call outfld('DSTSFMBL',flx_mss_vrt_dst_ttl,ncol,lchnk)
    do m=1,dst_nbr
       call outfld(flx_mss_mbl_sfc_nm(m),flx_mss_vrt_dst(1,m),ncol,lchnk)
    end do                    ! end loop over cst
#else /*!CCM */
    call outfld('BSN_FCT',mbl_bsn_fct,plond,lat_idx,obuf)
    call outfld('FRC_WET',frc_thr_ncr_wtr,plond,lat_idx,obuf)
    call outfld('GWC_SFC',gwc_sfc,plond,lat_idx,obuf)
    call outfld('LND_MBL',lnd_frc_mbl,plond,lat_idx,obuf)
    call outfld('TPT_GND',tpt_gnd,plond,lat_idx,obuf)
    call outfld('VAI_DST',vai_dst,plond,lat_idx,obuf)
    call outfld('VWC_SFC',vwc_sfc,plond,lat_idx,obuf)
    call outfld('WND_FRC',wnd_frc,plond,lat_idx,obuf)
    call outfld('WND_FRCS',wnd_frc_slt,plond,lat_idx,obuf)
    call outfld('WND_FRCT',wnd_frc_thr_slt,plond,lat_idx,obuf)
    call outfld('WND_RFR',wnd_rfr,plond,lat_idx,obuf)
    call outfld('WND_RFRT',wnd_rfr_thr_slt,plond,lat_idx,obuf)
    call outfld('DSTSFMBL',flx_mss_vrt_dst_ttl,plond,lat_idx,obuf)
    do m=1,dst_nbr
       call outfld(flx_mss_mbl_sfc_nm(m),flx_mss_vrt_dst(1,m),plond,lat_idx,obuf)
    end do                    ! end loop over cst
!++alfgr
    call outfld('SNW_FRC',snw_frc,plond,lat_idx,obuf)
    call outfld('LND_MBL',lnd_frc_mbl,plond,lat_idx,obuf)
    call outfld('WNDZNLMDP',wnd_znl_mdp,plond,lat_idx,obuf)
    call outfld('WNDMRDMDP',wnd_mrd_mdp,plond,lat_idx,obuf)
    call outfld('DSTSFMBL',flx_mss_vrt_dst_ttl,plond,lat_idx,obuf)
    call outfld('FLX_SWA',flx_SW_abs_sfc,plond,lat_idx,obuf)
    call outfld('FLX_LWD',flx_LW_dwn_sfc,plond,lat_idx,obuf)
    call outfld('WND_FRC_MBL',wnd_frc,plond,lat_idx,obuf)
    call outfld('TPT_MDP',tpt_mdp,plond,lat_idx,obuf)
    call outfld('Q_H2O',q_h2o_vpr,plond,lat_idx,obuf)
!--alfgr
#endif /*!CCM */
#ifdef DST_MSS_BDG
    call bdg_aal('dst_sf_mbl',lat_idx,flx_mss_vrt_dst_ttl)
    ! NB: Only call to bdg_gam_wet_2D is made here
    call bdg_gam_wet_2D('dst_ss_mbl',lat_idx,prs_dlt,q_dst_tnd_mbl_ttl)
#endif /*!DST_MSS_BDG */
    return
  end subroutine dst_mbl                       ! end dst_mbl()
  
  subroutine mbl2nc(             &
       dst_slt_flx_rat_ttl, & ! I [m-1] Ratio of vertical dust flux to streamwise mass flux
       fl_out,              & ! I [sng] Name of netCDF output file
       flx_mss_hrz_slt_ttl, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
       flx_mss_vrt_dst,     & ! I [kg m-2 s-1] Vertical mass flux of dust
       flx_mss_vrt_dst_ttl, & ! I [kg m-2 s-1] Total vertical mass flux of dust
       frc_thr_ncr_drg,     & ! I [frc] Threshold friction velocity increase from roughness
       frc_thr_ncr_wtr,     & ! I [frc] Threshold friction velocity increase from moisture
       gwc_sfc,             & ! I [kg kg-1] Gravimetric water content
       hgt_zpd,             & ! I [m] Zero plane displacement height
       lat_idx,             & ! I [idx] Model latitude index
       lnd_frc_mbl,         & ! I [frc] Bare ground fraction
       mno_lng,             & ! I [m] Monin-Obukhov length
       rgh_mmn,             & ! I [m] Roughness length momentum
       snw_frc,             & ! I [frc] Fraction of surface covered by snow
       vai_dst,             & ! I [m2 m-2] Vegetation area index, one-sided
       wnd_frc,             & ! I [m s-1] Friction velocity
       wnd_frc_slt,         & ! I [m s-1] Saltating friction velocity
       wnd_frc_thr_slt,     & ! I [m s-1] Threshold friction velocity for saltation
       wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
       wnd_rfr,             & ! I [m s-1] Wind speed at reference height
       wnd_rfr_thr_slt,     & ! I [m s-1] Threshold 10 m wind speed for saltation
       mss_frc_src,         & ! I [frc] Mass fraction of source distribution 
       mss_frc_trn_dst_src, & ! I [frc] Fraction of transported dust mass at source
       ovr_src_snk_mss,     & ! I [frc] Overlap of src with snk 
       ovr_src_snk_mss_ttl) ! I [frc] Total transported mass fraction of dust flux 
    ! Purpose: Output aerosol mobilization properties to netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstctl ! [mdl] Control variables, routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='mbl2nc' ! [sng] Subroutine name
    ! Input
    character(len=*) fl_out      ! [sng] Name of netCDF output file
    integer  lat_idx        ! I [idx] Model latitude index
    real(r8) dst_slt_flx_rat_ttl(plond) ! [m-1] Ratio of vertical dust flux to streamwise mass flux
    real(r8) flx_mss_hrz_slt_ttl(plond) ! [kg m-1 s-1] Vertically integrated streamwise mass flux
    real(r8) flx_mss_vrt_dst(plond,dst_nbr) ! [kg m-2 s-1] Vertical mass flux of dust
    real(r8) flx_mss_vrt_dst_ttl(plond) ! [kg m-2 s-1] Total vertical mass flux of dust
    real(r8) frc_thr_ncr_drg(plond) ! [frc] Threshold friction velocity increase from roughness
    real(r8) frc_thr_ncr_wtr(plond) ! [frc] Threshold friction velocity increase from moisture
    real(r8) gwc_sfc(plond)       ! [kg kg-1] Gravimetric water content
    real(r8) hgt_zpd              ! [m] Zero plane displacement
    real(r8) lnd_frc_mbl(plond)   ! [frc] Bare ground fraction
    real(r8) mno_lng(plond)       ! [m] Monin-Obukhov length
    real(r8) rgh_mmn              ! [m] Roughness length momentum
    real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
    real(r8) vai_dst(plond)       ! [m2 m-2] Vegetation area index, one-sided
    real(r8) wnd_frc(plond)       ! [m s-1] Friction velocity
    real(r8) wnd_frc_slt(plond)   ! [m s-1] Saltating friction velocity
    real(r8) wnd_frc_thr_slt(plond) ! [s m-1] Threshold friction velocity for saltation
    real(r8) wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
    real(r8) wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
    real(r8) wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
    real(r8) mss_frc_src(plond,dst_src_nbr) ! [frc] Mass fraction of source distribution
    real(r8) mss_frc_trn_dst_src(plond,dst_nbr) ! [frc] Fraction of transported dust mass at source
    real(r8) ovr_src_snk_mss(plond,dst_src_nbr,dst_nbr) ! [frc] Overlap of src with snk
    real(r8) ovr_src_snk_mss_ttl(plond) ! [frc] Total transported mass fraction of dust flux
    ! Output
    ! Local
    ! File metadata and dimension IDs
    integer dim_lon_sz_src_sz_time(4) ! [enm] Dimension IDs
    integer dim_lon_sz_src_time(3) ! [enm] Dimension IDs
    integer dim_lon_sz_time(3) ! [enm] Dimension IDs
    integer dim_lon_time(2)   ! [enm] Dimension IDs
    integer srt_lon_time(2)   ! Starting index array
    integer srt_lon_sz_src_sz_time(4) ! Starting index array
    integer srt_lon_sz_src_time(3) ! Starting index array
    integer cnt_lon_time(2)   ! Count array
    integer srt_lon_sz_time(3) ! Starting index array
    integer cnt_lon_sz_time(3) ! Count array
    integer cnt_lon_sz_src_sz_time(4) ! Count array
    integer cnt_lon_sz_src_time(3) ! Count array
    integer fll_mode_old      ! Old fill mode
    integer lon_dim_id        ! [enm] Dimension ID for lon
    integer time_dim_id       ! [enm] Dimension ID for time
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer sz_dim_id         ! [enm] Dimension ID for sz
    integer sz_src_dim_id     ! [enm] Dimension ID for sz_src
    integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
    ! Variable IDs
    integer dst_slt_flx_rat_ttl_id ! [enm] Variable ID
    integer flx_mss_hrz_slt_ttl_id ! [enm] Variable ID
    integer flx_mss_vrt_dst_id ! [enm] Variable ID
    integer flx_mss_vrt_dst_ttl_id ! [enm] Variable ID
    integer frc_thr_ncr_drg_id ! [enm] Variable ID
    integer frc_thr_ncr_wtr_id ! [enm] Variable ID
    integer gwc_sfc_id        ! [enm] Variable ID
    integer hgt_zpd_id        ! [enm] Variable ID
    integer lnd_frc_mbl_id    ! [enm] Variable ID
    integer mno_lng_id        ! [enm] Variable ID
    integer rgh_mmn_id        ! [enm] Variable ID
    integer snw_frc_id        ! [enm] Variable ID
    integer vai_dst_id        ! [enm] Variable ID
    integer wnd_frc_id        ! [enm] Variable ID
    integer wnd_frc_slt_id    ! [enm] Variable ID
    integer wnd_frc_thr_slt_id ! [enm] Variable ID
    integer wnd_mdp_id        ! [enm] Variable ID
    integer wnd_rfr_id        ! [enm] Variable ID
    integer wnd_rfr_thr_slt_id ! [enm] Variable ID
    integer mss_frc_src_id         ! [enm] Variable ID
    integer mss_frc_trn_dst_src_id ! [enm] Variable ID
    integer ovr_src_snk_mss_id ! [enm] Variable ID
    integer ovr_src_snk_mss_ttl_id ! [enm] Variable ID
    ! Variable data
    
    ! Begin netCDF output routines
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=rcd+nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
    rcd=rcd+nf90_redef(nc_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'sz',sz_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'sz_src',sz_src_dim_id)
    ! Add global attributes
    ! Define dimension IDs
    ! Assemble ID and count vectors for each multidimensional combination of dimensions
    dim_lon_time=(/lon_dim_id,time_dim_id/)
    srt_lon_time=(/1,nstep/)  ! Starting index array
    cnt_lon_time=(/plon,1/)   ! Count array
    
    dim_lon_sz_time=(/lon_dim_id,sz_dim_id,time_dim_id/)
    srt_lon_sz_time=(/1,1,nstep/)
    cnt_lon_sz_time=(/plon,dst_nbr,1/)

    dim_lon_sz_src_sz_time=(/lon_dim_id,sz_src_dim_id,sz_dim_id,time_dim_id/)
    srt_lon_sz_src_sz_time=(/1,1,1,nstep/)
    cnt_lon_sz_src_sz_time=(/plon,dst_src_nbr,dst_nbr,1/)

    dim_lon_sz_src_time=(/lon_dim_id,sz_src_dim_id,time_dim_id/)
    srt_lon_sz_src_time=(/1,1,nstep/)
    cnt_lon_sz_src_time=(/plon,dst_src_nbr,1/)
    
    if (nstep == 1) then
       ! Variable definitions
       rcd=rcd+nf90_def_var(nc_id,'frc_thr_ncr_drg',nf90_float,dim_lon_time,frc_thr_ncr_drg_id)
       rcd=rcd+nf90_def_var(nc_id,'dst_slt_flx_rat_ttl',nf90_float,dim_lon_time,dst_slt_flx_rat_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_hrz_slt_ttl',nf90_float,dim_lon_time,flx_mss_hrz_slt_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_vrt_dst',nf90_float,dim_lon_sz_time,flx_mss_vrt_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_vrt_dst_ttl',nf90_float,dim_lon_time,flx_mss_vrt_dst_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'frc_thr_ncr_wtr',nf90_float,dim_lon_time,frc_thr_ncr_wtr_id)
       rcd=rcd+nf90_def_var(nc_id,'gwc_sfc',nf90_float,dim_lon_time,gwc_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'hgt_zpd_mbl',nf90_float,hgt_zpd_id)
       rcd=rcd+nf90_def_var(nc_id,'lnd_frc_mbl',nf90_float,dim_lon_time,lnd_frc_mbl_id)
       rcd=rcd+nf90_def_var(nc_id,'mno_lng_mbl',nf90_float,dim_lon_time,mno_lng_id)
       rcd=rcd+nf90_def_var(nc_id,'rgh_mmn_mbl',nf90_float,rgh_mmn_id)
       rcd=rcd+nf90_def_var(nc_id,'snw_frc',nf90_float,dim_lon_time,snw_frc_id)
       rcd=rcd+nf90_def_var(nc_id,'vai_dst',nf90_float,dim_lon_time,vai_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_frc_mbl',nf90_float,dim_lon_time,wnd_frc_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_frc_slt',nf90_float,dim_lon_time,wnd_frc_slt_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_frc_thr_slt',nf90_float,dim_lon_time,wnd_frc_thr_slt_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_mdp',nf90_float,dim_lon_time,wnd_mdp_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_rfr_mbl',nf90_float,dim_lon_time,wnd_rfr_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_rfr_thr_slt',nf90_float,dim_lon_time,wnd_rfr_thr_slt_id)
       rcd=rcd+nf90_def_var(nc_id,'mss_frc_src',nf90_float,dim_lon_sz_src_time,mss_frc_src_id)
       rcd=rcd+nf90_def_var(nc_id,'mss_frc_trn_dst_src',nf90_float,dim_lon_sz_time,mss_frc_trn_dst_src_id)
       rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_mss',nf90_float,dim_lon_sz_src_sz_time,ovr_src_snk_mss_id)
       rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_mss_ttl',nf90_float,dim_lon_time,ovr_src_snk_mss_ttl_id)
       ! Add english text descriptions
       rcd=rcd+nf90_put_att(nc_id,dst_slt_flx_rat_ttl_id,'long_name','Ratio of vertical dust flux to streamwise mass flux')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_hrz_slt_ttl_id,'long_name','Vertically integrated streamwise mass flux')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_id,'long_name','Vertical mass flux of dust')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_ttl_id,'long_name','Total vertical mass flux of dust')
       rcd=rcd+nf90_put_att(nc_id,frc_thr_ncr_drg_id,'long_name','Threshold friction velocity increase from roughness')
       rcd=rcd+nf90_put_att(nc_id,frc_thr_ncr_wtr_id,'long_name','Threshold friction velocity increase from moisture')
       rcd=rcd+nf90_put_att(nc_id,gwc_sfc_id,'long_name','Gravimetric water content')
       rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'long_name','Zero plane displacement height')
       rcd=rcd+nf90_put_att(nc_id,lnd_frc_mbl_id,'long_name','Bare ground fraction')
       rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'long_name','Monin-Obukhov length')
       rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'long_name','Roughness length momentum')
       rcd=rcd+nf90_put_att(nc_id,snw_frc_id,'long_name','Fraction of surface covered by snow')
       rcd=rcd+nf90_put_att(nc_id,vai_dst_id,'long_name','Vegetation area index, one-sided')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'long_name','Friction velocity')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_slt_id,'long_name','Saltating friction velocity')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_thr_slt_id,'long_name','Threshold friction velocity for saltation')
       rcd=rcd+nf90_put_att(nc_id,wnd_mdp_id,'long_name','Surface layer mean wind speed')
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'long_name','Wind speed at reference height')
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_thr_slt_id,'long_name','Threshold 10 m wind speed for saltation')
       rcd=rcd+nf90_put_att(nc_id,mss_frc_src_id,'long_name','Mass fraction of source distribution')
       rcd=rcd+nf90_put_att(nc_id,mss_frc_trn_dst_src_id,'long_name','Fraction of transported dust mass at source')
       rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_id,'long_name','Mass overlap of src dist. i with sink bin j')
       rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_ttl_id,'long_name','Total transported mass fraction of dust flux')
       ! Add units
       rcd=rcd+nf90_put_att(nc_id,dst_slt_flx_rat_ttl_id,'units','meter-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_hrz_slt_ttl_id,'units','kilogram meter-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_vrt_dst_ttl_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,frc_thr_ncr_drg_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,frc_thr_ncr_wtr_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,gwc_sfc_id,'units','kilogram kilogram-1')
       rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,lnd_frc_mbl_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,snw_frc_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,vai_dst_id,'units','meter2 meter-2')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_slt_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_thr_slt_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,wnd_mdp_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_thr_slt_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_mss_ttl_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,mss_frc_src_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,mss_frc_trn_dst_src_id,'units','fraction')
    else                      ! endif nstep == 1
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,'dst_slt_flx_rat_ttl',dst_slt_flx_rat_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_hrz_slt_ttl',flx_mss_hrz_slt_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_vrt_dst',flx_mss_vrt_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_vrt_dst_ttl',flx_mss_vrt_dst_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'frc_thr_ncr_drg',frc_thr_ncr_drg_id)
       rcd=nf90_wrp_inq_varid(nc_id,'frc_thr_ncr_wtr',frc_thr_ncr_wtr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'gwc_sfc',gwc_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'hgt_zpd_mbl',hgt_zpd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'lnd_frc_mbl',lnd_frc_mbl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'mno_lng_mbl',mno_lng_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rgh_mmn_mbl',rgh_mmn_id)
       rcd=nf90_wrp_inq_varid(nc_id,'snw_frc',snw_frc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vai_dst',vai_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_frc_mbl',wnd_frc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_frc_slt',wnd_frc_slt_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_frc_thr_slt',wnd_frc_thr_slt_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_mdp',wnd_mdp_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_rfr_mbl',wnd_rfr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_rfr_thr_slt',wnd_rfr_thr_slt_id)
       rcd=nf90_wrp_inq_varid(nc_id,'mss_frc_src',mss_frc_src_id)
       rcd=nf90_wrp_inq_varid(nc_id,'mss_frc_trn_dst_src',mss_frc_trn_dst_src_id)
       rcd=nf90_wrp_inq_varid(nc_id,'ovr_src_snk_mss',ovr_src_snk_mss_id)
       rcd=nf90_wrp_inq_varid(nc_id,'ovr_src_snk_mss_ttl',ovr_src_snk_mss_ttl_id)
    endif                     ! endif nstep /= 1
    ! End define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    if (nstep == 1) then
       rcd=rcd+nf90_put_var(nc_id,hgt_zpd_id,hgt_zpd)
       rcd=rcd+nf90_put_var(nc_id,rgh_mmn_id,rgh_mmn)
    endif                     ! endif nstep == 1
    rcd=rcd+nf90_put_var(nc_id,dst_slt_flx_rat_ttl_id,dst_slt_flx_rat_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,gwc_sfc_id,gwc_sfc,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_hrz_slt_ttl_id,flx_mss_hrz_slt_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_vrt_dst_id,flx_mss_vrt_dst,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_vrt_dst_ttl_id,flx_mss_vrt_dst_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,frc_thr_ncr_drg_id,frc_thr_ncr_drg,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,frc_thr_ncr_wtr_id,frc_thr_ncr_wtr,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,lnd_frc_mbl_id,lnd_frc_mbl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,mno_lng_id,mno_lng,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,snw_frc_id,snw_frc,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,vai_dst_id,vai_dst,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_frc_id,wnd_frc,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_frc_slt_id,wnd_frc_slt,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_frc_thr_slt_id,wnd_frc_thr_slt,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_mdp_id,wnd_mdp,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_rfr_id,wnd_rfr,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_rfr_thr_slt_id,wnd_rfr_thr_slt,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_src_id,mss_frc_src,start=srt_lon_sz_src_time,count=cnt_lon_sz_src_time)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_trn_dst_src_id,mss_frc_trn_dst_src,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_mss_id,ovr_src_snk_mss,start=srt_lon_sz_src_sz_time,count=cnt_lon_sz_src_sz_time)
    rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_mss_ttl_id,ovr_src_snk_mss_ttl,start=srt_lon_time,count=cnt_lon_time)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
    if (nstep == 1) then
       write (6,'(a,a43,a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': Initialized mobilization data archive in ',fl_out(1:ftn_strlen(fl_out))
    endif                     ! endif
  end subroutine mbl2nc                       ! end mbl2nc()
  
end module dstmbl

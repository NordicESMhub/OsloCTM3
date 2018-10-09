! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstdpswet.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: dstdpswet.F contains routines which represent mineral dust 
! wet deposition tendencies

! Usage
! use dstdpswet ! [mdl] Wet deposition driver

! Implementation is designed with two goals in mind:
! 1. Make scheme work with CCM3 clouds for Paleo-CSM
! 2. Make scheme work with PJR's clouds for CCM4 and MATCH
! PJR's prognostic scheme contains a fuller set of wet deposition predictors than CCM3 cloud scheme 
! In addition to q_H2O_cnd2pcp_tnd, which is available in CCM3, PJR predicts q_H2O_wtr2pcp_ttl_tnd and q_H2O_pcp2vpr_tnd
! These are necessary for "exact" representation of wet deposition
! Assuming a simple scavenging ratio is sufficient to parameterize the 
! wet deposition process, the q_H2O_wtr2pcp_ttl_tnd, q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd tendencies suffice to predict dust rain-out rates.
! However, vanilla CCM3 only provides dqcond, the net vapor tendency due to condensavtion and evaporation
! CCM3 does not provide condensation and evaporation tendencies separately
! Thus we map available CCM3 tendencies onto tendencies of same name that PJR provides
! Thus q_H2O_wtr2pcp_ttl_tnd and q_H2O_pcp2vpr_tnd in CCM3 do not mean quite the same thing as they do in CCM4
! However, this strategy has the advantage that once inputs are defined in either
! version of the model (CCM3 or CCM4), dust prediction can use the same code.

! Cloud nomenclature:
! cld_frc = [frc] Local total cloud fraction
! = Horizontal fraction of gridcell covered by cloud
! = MATCH:cldt
! cld_frc_cnv = [frc] Convective cloud fraction
! = Horizontal fraction of gridcell covered by convective cloud
! = CCM/MATCH:cldc
! cld_vlm = [frc] Cloud volume
! = Fraction of gridcell volume occupied by rain or cloud water
! = Is larger than cld_frc when rain from above passes through gridcell
! = MATCH:cldv

! Tracer nomenclature:
! q_H2O_vpr = [kg kg-1] Water vapor mixing ratio
! = CCM/MATCH:q

! q_H2O_cnd = [kg kg-1] Condensed H2O mixing ratio
! = Cloud water amount
! = Sum of liquid and ice water
! = MATCH:as1 if diagnosed, cwat iff archived

! q_H2O_cnd_cnv = [kg kg-1] Condensed H2O mixing ratio in convective clouds
! = fxm
! = MATCH:conicw,ql in conv_ccm_pjr()

! q_H2O_pcp_lqd = ! [kg kg-1] Rainwater mixing ratio
! = Mixing ratio of liquid precipitation (for wet deposition purposes)
! = MATCH:rain

! q_H2O_vpr2cnd_tnd = [kg kg-1 s-1] H2O vapor to condensate tendency
! q_H2O_cnd_tnd = [kg kg-1 s-1] Net H2O condensate formation tendency
! = Condensation - evaporation 
! = This is net condensation minus evaporation, so may be positive or negative 
! = Positive when condensation exceeds evaporation
! = Negative when evaporation exceeds condensation
! = Evaporation here refers to in-cloud evaporation, and does not include
! = evaporation from falling rain which is counted separately in q_H2O_pcp2vpr
! = MATCH:cme but MATCH sometimes uses this as an in-cloud rate
! = MATCH:conds

! q_H2O_pcp2vpr_tnd = [kg kg-1 s-1] H2O precipitation to vapor tendency
! = H2O stratiform precipitation to vapor tendency
! = Evaporation of non-local, stratiform rain into local gridcell
! = Rate of evaporation of falling stratiform precipitation
! = Does not include any evaporation from in-cloud water
! = Evaporation from convective rain (convective downdrafts) is not included (it may be accounted for directly in convection routines, I do not know, but it is not included in this term)
! = Positive definite
! = MATCH:evapr,evaps

! q_H2O_cnd2pcp_tnd = [kg kg-1 s-1] Condensed H2O to precipitation tendency
! = Local water to stratiform precipitation tendency
! = Stratiform precipitation tendency
! = CCM 3.6: vpr2pcp as predicted by physics/cond()
! = CCM 3.7+: stratiform precipitation is taken directly from condensate
! = Production tendency only, does not include local evaporation of falling rain
! = Positive definite
! = MATCH:prain,precs

! q_H2O_vpr2pcp_cnv_tnd = [kg kg-1 s-1] H2O vapor to convective precipitation tendency
! = Precipitation formation tendency from convection
! = Production tendency only, does not include local evaporation of falling rain
! = Positive definite
! = CCM/MATCH:cmfdqr

! q_H2O_wtr2pcp_ttl_tnd = [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
! = Local water (vapor + condensate) to total (convective + stratiform) precipitation tendency
! = Sum of convective and stratiform precipitation production tendencies
! = Positive definite
! = q_H2O_cnd2pcp_tnd + q_H2O_vpr2pcp_cnv_tnd

! Precipitation flux nomenclature:
! pcp_flx_sfc = [kg m-2 s-1] Total precipitation reaching surface
! = CCM:prect (but in mass flux units not m s-1)

! pcp_flx_str = [kg m-2 s-1] Stratiform precipitation flux
! = CCM:precl (but in mass flux units not m s-1)

! pcp_flx_cnv = [kg m-2 s-1] Convective precipitation flux
! = CCM:precc (but in mass flux units not m s-1)

! wtr2pcp_flx = [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
! = Total gridcell precipitation formation mass flux
! = Sum of vapor->precipitation and condensate->precipitation
! = Until 20000413 wtr2pcp_flx was mislabeled as vpr2pcp_flx but was defined correctly (=vpr2pcp_flx+cnd2pcp_flx)
! = CCM/MATCH: q_H2O_wtr2pcp_ttl_tnd*mpl_air

! Total tendency equations in terms of the above:
! tpt_mdp = tpt_mdp + &
!    latvap/cpair*(q_H2O_vpr2cnd_tnd-q_H2O_pcp2vpr)*tm_dlt +
!    rmelt(i,k)*tm_dlt
! q_H2O_vpr = q_H2O_vpr - &
!    (q_H2O_vpr2cnd_tnd - q_H2O_pcp2vpr)*tm_dlt
! q_H2O_cnd = q_H2O_cnd + &
!    (q_H2O_vpr2cnd_tnd - q_H2O_cnd2pcp_tnd(i,k))*tm_dlt
! pcp_flx_str(i) = pcp_flx_str(i) + &
!    (q_H2O_cnd2pcp_tnd(i,k)-q_H2O_pcp2vpr)*prs_dlt/grv_mean_sfc

! Requires dst.h for BXM, CCM, DST_MSS_BDG tokens
#include <dst.h> /* Dust preprocessor tokens */

module dstdpswet ! [mdl] Wet deposition driver
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::cld_dgn ! [sbr] Cloud diagnostics
  public::dst_dps_wet ! [sbr] Wet deposition driver
  public::dst_dps_wet_old ! [sbr] Wet deposition driver, old
  
contains
  
  subroutine dst_dps_wet(lchnk,ncol,obuf, & ! [sbr] Wet deposition driver
       cld_frc,             & ! I [frc] Local total cloud fraction
       cld_frc_cnv,         & ! I/O [frc] Convective cloud fraction
       cld_vlm,             & ! I [frc] Cloud volume
       frc_trc_trn_cnv_ptn, & ! O [frc] Interstitial tracer fraction
       lat_idx,             & ! I [idx] Model latitude index
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_cnd,           & ! I [kg kg-1] Condensed H2O mixing ratio
       q_H2O_cnd_cnv,       & ! I [kg kg-1] Condensed H2O mixing ratio in convective clouds
       q_H2O_cnd2pcp_tnd,   & ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
       q_H2O_pcp2vpr_tnd,   & ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
       q_H2O_vpr2pcp_cnv_tnd, & ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
       q_dst,               & ! I/O [kg kg-1] Dust mixing ratio
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp)             ! I [K] Temperature
    ! q_H2O_vpr2cnd_tnd,  & ! I [kg kg-1 s-1] H2O vapor to condensate tendency
    ! Purpose: Simulate wet deposition processes 
    ! dst_dps_wet() is called by CCM:physics/aphys(), MATCH:src/physlic()
    ! Algorithm based on Barth, Rasch, and Kiehl (2000) (BRK00)
    use dstbdg,only:bdg_aal,bdg_gam_wet ! [mdl] Mass budget diagnostics
    use dstblm,only:tpt_frz_pnt ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use dstcst ! [mdl] Physical constants for dust routines
    use dstdbg ! [mdl] Debugging information for dust model
    use dstgrd,only:dst_nbr ! [mdl] Dust grid sizes
    use dstmssutl,only:dst_add_lon_lev,dst_add_lon ! [mdl] Mass budget utilities
    use dstnm ! [mdl] Nomenclature for outfld()
    use dstscv ! [mdl] Aerosol scavenging properties
#ifndef BXM
    use histout,only:outfld ! [mdl] History/archiving
#endif /* BXM */
    use pmgrid ! [mdl] Spatial resolution parameters
    use precision ! [mdl] Precision r8, i8, ...
    use sng_mdl,only:ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Parameters
    real(r8),parameter::pcp_flx_str_min=1.0e-12 ! [kg m-2 s-1] Minimum stratiform precipitation flux for determining evaporated fraction
    real(r8),parameter::q_H2O_cnd_min=1.0e-12 ! [kg kg-1] Minimum condensed water for determining condensation to precipitation conversion fraction
    real(r8),parameter::cld_frc_cll_min=1.0e-5 ! [frc] Minimum cloud fraction for collision scavenging for determining potential scavenging
    real(r8),parameter::q_dst_cld_min=1.0e-30 ! [kg kg-1] Epsilon to prevent divide by zero for interstitial fraction
    real(r8),parameter::q_trc_tnd_scv_eps=1.0e-30 ! [kg kg-1 s-1] Epsilon to prevent to prevent divide by zero in denominators of scavenged fractions
    
    ! Input
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
    real(r8),intent(in)::cld_frc(plond,plevp) ! I [frc] Local total cloud fraction
    real(r8),intent(in)::cld_vlm(plond,plev) ! I [frc] Cloud volume
    real(r8),intent(in)::q_H2O_cnd(plond,plev) ! I [kg kg-1] Condensed H2O mixing ratio
    real(r8),intent(in)::q_H2O_cnd_cnv(plond,plev) ! I [kg kg-1] Condensed H2O mixing ratio in convective clouds
    real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
    real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure
    real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
    real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
    ! Output
    real(r8),intent(out)::frc_trc_trn_cnv_ptn(plond,plev,dst_nbr) ! O [frc] Interstitial tracer fraction
    ! Input/Output
    ! Fields besides q_dst are only altered when they are physically inconsistent
    real(r8),intent(inout)::cld_frc_cnv(plond,plev) ! I/O [frc] Convective cloud fraction
    real(r8),intent(inout)::q_H2O_pcp2vpr_tnd(plond,plev) ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
    real(r8),intent(inout)::q_H2O_cnd2pcp_tnd(plond,plev) ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
    real(r8),intent(inout)::q_dst(plond,plev,dst_nbr) ! I/O [kg dust/kg moist air]
    ! Local Output
    real(r8) flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
    real(r8) flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
    real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Total precipitation reaching surface
    real(r8) prect2(plond)        ! [m s-1] Total precipitation reaching surface (diagnostic, should be same as PRECT)
    real(r8) q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
    real(r8) q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
    real(r8) q_dst_tnd_ncl(plond,plev,dst_nbr) ! O [kg kg-1 s-1] Nucleation scavenging tendency (positive definite)
    real(r8) q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
    real(r8) q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
    real(r8) q_dst_tnd_wet(plond,plev,dst_nbr) ! [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
    real(r8) spc_xsx_ncl_scv(dst_nbr) ! [m2 kg-1] Specific cross section for nucleation scavenging, GBS98, BJG93
    ! Local
    character(80) fl_out      ! [sng] Name of netCDF output file
    integer i                 ! [idx] lon index
    integer k                 ! [idx] lev index
    integer m                 ! [idx] cst index 
    real(r8) cld_frc_cll_cnv      ! [frc] Cloud fraction for collision scavenging by convective rain
    real(r8) cld_frc_cll_str      ! [frc] Cloud fraction for collision scavenging by stratiform rain
    real(r8) cld_frc_eff          ! [frc] Fraction of dust removable by precip
    real(r8) cld_frc_ncl_cnv      ! [frc] Cloud fraction for nucleation scavenging in convective clouds
    real(r8) cld_frc_ncl_str      ! [frc] Cloud fraction for nucleation scavenging in stratiform clouds
    real(r8) cnd_frc_ice          ! [frc] Ice fraction of condensate
    real(r8) dst2pcp_flx          ! [kg m-2 s-1] Local removal of dust by precip (q_H2O_wtr2pcp_ttl_tnd*q_dst*z_scv)
    real(r8) flx_mss_evp          ! [kg m-2 s-1] Local source of dust due to evaporation of precipitation from above (pcpfrc_evp*flx_mss_pcp_sfc)
    real(r8) flx_mss_scv_abv_cnv(plond,dst_nbr) ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
    real(r8) flx_mss_scv_abv_str(plond,dst_nbr) ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
    real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
    real(r8) pcp2vpr_flx          ! [kg m-2 s-1] Local conversion of precip to vapor (q_H2O_pcp2vpr_tnd*mpl_air)
    real(r8) pcp_flx_cnv(plond)   ! [kg m-2 s-1] Convective precipitation flux
    real(r8) pcp_flx_str(plond)   ! [kg m-2 s-1] Stratiform precipitation flux
    real(r8) pcp_str_frc_evp      ! [frc] Stratiform precipitation from above that evaporates locally (pcp2vpr_flx/pcp_flx_sfc)
    real(r8) ptn_frc_scv_cll_cnv  ! [frc] Potential tracer fraction removed by collision scavenging in convective precipitation
    real(r8) ptn_frc_scv_cll_str  ! [frc] Potential tracer fraction removed by collision scavenging in stratiform precipitation
    real(r8) ptn_frc_scv_ncl_cnv  ! [frc] Potential tracer fraction removed by nucleation scavenging in convective precipitation
    real(r8) ptn_frc_scv_ncl_str  ! [frc] Potential tracer fraction removed by nucleation scavenging in stratiform precipitation
    real(r8) q_H2O_cnd2pcp_frc    ! [frc] Fraction of condensate converted to (stratiform) precipitation
    real(r8) q_H2O_wtr2pcp_ttl_tnd(plond,plev) ! [kg kg-1 s-1] Local water to total precipitation tendency
    real(r8) q_trc_tnd_scv_cll_cnv ! [kg kg-1 s-1] Convective rain-collision scavenging tracer flux
    real(r8) q_trc_tnd_scv_cll_str ! [kg kg-1 s-1] Stratiform rain-collision scavenging tracer tendency
    real(r8) q_trc_tnd_scv_cnv    ! [kg kg-1 s-1] Total convective scavenging tracer tendency
    real(r8) q_trc_tnd_scv_ncl_cnv ! [kg kg-1 s-1] Convective nucleation scavenging tracer tendency
    real(r8) q_trc_tnd_scv_ncl_str ! [kg kg-1 s-1] Stratiform nucleation scavenging tracer tendency
    real(r8) q_trc_tnd_scv_str    ! [kg kg-1 s-1] Total stratitform scavenging tracer tendency
    real(r8) scv_ncl_frc_cnv      ! [frc] Nucleation scavenging fraction in convective clouds
    real(r8) scv_ncl_frc_str      ! [frc] Nucleation scavenging fraction in stratiform clouds
    real(r8) scv_tnd_fct          ! Factor limiting scavenging tendencies to amount present [frc]
    real(r8) tm_dlt               ! [s] Wet deposition timestep
    real(r8) tpt_cls              ! [C] Temperature
    real(r8) wtr2pcp_flx          ! [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
    ! real(r8) cld_frc_max_abv_cnv  ! [frc] Maximum convective cloud fraction at or above this level
    ! real(r8) cld_frc_max_abv_str  ! [frc] Maximum stratiform cloud fraction at or above this level
    ! real(r8) mpc_trc_abv          ! [kg m-2] Column integrated tracer path
    
    ! Main Code
    ! Timesplit if desired
    tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
    
    ! Initialize fluxes and tendencies
    q_dst_tnd_pcp(:,:,:)=0.0_r8         ! [kg kg-1 s-1] Scavenging tendency
    q_dst_tnd_evp(:,:,:)=0.0_r8         ! [kg kg-1 s-1] Evaporation tendency
    flx_mss_pcp_sfc(:,:)=0.0_r8       ! [kg m-2 s-1]
    
    ! Specific cross-sections adopted were defined empirically by Balkanski
    ! fxm: specific cross-sections should use rain and aerosol size information
    ! Set size and species dependent scavenging cross-sections
    do m=1,dst_nbr
       ! NB: RBK00 use 0.1 for all first order scavenging coefficients, based on Balkanski
       spc_xsx_ncl_scv(m)=0.10_r8 ! [m2 kg-1] Specific cross section for nucleation scavenging, csz 50% of value used for hygrophilic aerosols
       ! scv_cff_mss_avg_pcp_nrm_cnv(m)=0.05_r8 ! [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, convective csz 50% of value used for hygrophilic aerosols
       ! scv_cff_mss_avg_pcp_nrm_str(m)=0.05_r8 ! [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, stratiform csz 50% of value used for hygrophilic aerosols
    end do                    ! end loop over cst
    
    ! Sanity checks
    ! Alter fields that are physically inconsistent
    ! fxm: Put this loop in an #ifdef DST_DBG block once the routines responsible are fixed
    do k=1,plev
       do i=1,plon
          ! fxm: q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd may be < 0.0, though I do not understand why, so just make them positive definite
          q_H2O_pcp2vpr_tnd(i,k)=max(q_H2O_pcp2vpr_tnd(i,k),0.0_r8) ! [kg kg-1 s-1] 
          q_H2O_cnd2pcp_tnd(i,k)=max(q_H2O_cnd2pcp_tnd(i,k),0.0_r8) ! [kg kg-1 s-1] 
          ! if (q_H2O_pcp2vpr_tnd(i,k) < 0.0_r8) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_pcp2vpr_tnd(',i,',',k,') = ',q_H2O_pcp2vpr_tnd(i,k)
          ! if (q_H2O_cnd2pcp_tnd(i,k) < 0.0_r8) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_cnd2pcp_tnd(',i,',',k,') = ',q_H2O_cnd2pcp_tnd(i,k)
          ! 20000612: MATCH sets max(cld_frc=0.999) but allows cld_frc_cnv=1.0
          cld_frc_cnv(i,k)=min(cld_frc_cnv(i,k),cld_frc(i,k))
          ! if (cld_frc(i,k) < cld_frc_cnv(i,k)) then
          ! write(6,'(a)') 'dst: ERROR dst_dps_wet() reports cloud fraction error:'
          ! write(6,'(a,i2,a,i3,a,i2,a,es15.9,a,es15.9)') &
          !     'dst_pcp: lat = ',lat_idx,' cld_frc(',i,',',k,') = ',cld_frc(i,k) &
          !     ,' < cld_frc_cnv = ',cld_frc_cnv(i,k)
          ! endif               ! endif
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    ! 100 format(a,i2,a,i3,a,i2,a,es8.1) 
    
    ! Compute necessary derived fields
    do k=1,plev
       do i=1,plon
          ! Mass of air currently in gridbox
          mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
          ! Total precipitation formation
          q_H2O_wtr2pcp_ttl_tnd(i,k)= & ! [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
               q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    ! Keep track of amount of stratiform precip and tracer falling from above
    pcp_flx_str(:)=0.0_r8           ! [kg m-2 s-1] Stratiform precipitation flux
    pcp_flx_cnv(:)=0.0_r8           ! [kg m-2 s-1] Convective precipitation flux
    flx_mss_scv_abv_str(:,:)=0.0_r8   ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
    flx_mss_scv_abv_cnv(:,:)=0.0_r8   ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
    
    ! Notes on algorithm:
    ! This algorithm originated by PJR and MKB for highly soluble aerosols (sulfate)
    ! For some species, e.g., HNO3, any tracer within cloud volume is in cloud water
    ! For other species, e.g., mineral dust, amount in cloud water varies
    ! Convective clouds: We do not know cloud water amount 
    ! Assume all convective cloud water (and thus all tracer in convective cloud water) falls out each time step
    ! Stratiform clouds: Fraction of cloud water converted to precipitation defines amount of tracer which is scavenged
    ! All clouds: Scavenging refers to "removal from gridcell"
    ! Removal from gridcell means gridcell aerosol mass/number decreases
    ! Some would include CCN activation (without precipitation) as scavenging
    ! While activation means aerosol no longer competes as CCN, or competes
    ! for H2O condensation as haze particle, activated droplets are still present
    ! in gridcell so aerosol cannot be considered "removed" from mass budget.
    ! Thus vapor condensation does not contribute to nucleation scavenging
    ! because condensate must precipitate before aerosol is removed from model.
    
    ! Invert loops (lon is inner dimension) process from top down on each longitude
    ! Surface fluxes accumulate as local scavenging and evaporation modify falling precipitation
    ! Top-down structure of algorithm based on CCM2 cond() routine of J. J. Hack
    do k=1,plev               ! NB: Inefficient loop order
       do i=1,plon
          
          tpt_cls=tpt_mdp(i,k)-tpt_frz_pnt ! [C] Temperature
          
          ! Ice fraction is ignored in dust wet deposition, but is used by other aerosols (e.g., SO4)
          cnd_frc_ice=max(0.0_r8,min(-tpt_cls*0.05_r8,1.0_r8)) ! [frc] Ice fraction of condensate
          
          ! Convert local condensation and evaporation tendencies to fluxes
          wtr2pcp_flx=q_H2O_wtr2pcp_ttl_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Local conversion of water to precip
          pcp2vpr_flx=q_H2O_pcp2vpr_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Local conversion of precip to vapor
          
          ! Compute fraction of stratiform rain from above that evaporates in this gridcell
          ! This fraction will be applied to tracer in stratitform rain from above to compute
          ! tracer released into this gridcell during evaporation of stratiform rain.
          ! Use safety limits because pcp_str_frc_evp often equals, e.g., 1.0000000000000007
          pcp_str_frc_evp=min(pcp2vpr_flx/max(pcp_flx_str_min,pcp_flx_str(i)),1.0_r8) ! [frc]
          
          ! Compute fraction of condensate converted to (stratiform) precipitation
          ! This fraction determines what fraction of tracer in stratiform clouds is potentially nucleation scavenged
          ! [frc] Fraction of condensate converted to (stratiform) precipitation
          q_H2O_cnd2pcp_frc=max(0.0_r8,min(1.0_r8, &
               q_H2O_cnd2pcp_tnd(i,k)*tm_dlt/max(q_H2O_cnd(i,k),q_H2O_cnd_min)))
          ! fxm: PJR and MCB decided not to add condensate formed in current timestep to denominator, why?
          ! q_H2O_cnd2pcp_tnd(i,k)*tm_dlt/max(q_H2O_cnd(i,k)+q_H2O_vpr2cnd_tnd(i,k)*tm_dlt,q_H2O_cnd_min)))
          
          ! Cloud fraction properties for each process are independent of aerosol:
          ! MCB and PJR considered multiple candidates for each process, but chose
          ! to use cld_vlm for all SO4 processes
          
          ! Cloud used in convective nucleation scavenging should be local gridcell amount
          ! MCB and PJR considered cld_frc, cld_frc_cnv, and cld_vlm for SO4
          cld_frc_ncl_cnv=cld_frc_cnv(i,k) ! [frc] Cloud fraction for nucleation scavenging in convective clouds
          ! Cloud in convective rain-collision scavenging depends on local and higher cloud
          ! MCB and PJR considered max(cld_frc,cld_frc_abv_cnv), max(cld_frc_cnv,cld_frc_abv_cnv), max(cld_vlm,cld_frc_abv_cnv), and simply cld_vlm for SO4
          cld_frc_cll_cnv=cld_vlm(i,k) ! [frc] Cloud fraction for collision scavenging by convective rain
          ! Cloud used in stratiform nucleation scavenging should be local gridcell amount
          ! MCB and PJR considered cld_frc-cld_frc_cnv, cld_frc, and cld_vlm for SO4
          cld_frc_ncl_str=max(0.0_r8,cld_frc(i,k)-cld_frc_cnv(i,k)) ! [frc] Cloud fraction for nucleation scavenging in stratiform clouds
          ! Cloud in stratiform rain-collision scavenging depends on local and higher cloud
          ! MCB and PJR considered max(cld_frc,cld_frc_abv_str), max(cld_vlm,cld_frc_abv_str), and simply cld_vlm for SO4
          cld_frc_cll_str=cld_vlm(i,k) ! [frc] Cloud fraction for collision scavenging by stratiform rain
          
          ! Scavenging processes may be species and size dependent
          do m=1,dst_nbr
             
             ! Scavenging in convective clouds:
             ! Nucleation scavenging (SeP97 p. 1027):
             ! In-cloud scavenging occurs through nucleation scavenging and interstitial aerosol collection by cloud droplets (as opposed to raindrops)
             ! Nucleation scavenging is growth (activation) of CCN to cloud drops
             ! Textbook definition is not currently useful for model implementation
             ! As described above, model defines in-cloud scavenging as _removal_ 
             ! of aerosol from local gridcell by locally formed precipitation
             ! Thus activation and collection in cloud droplets _do not_ contribute
             ! to model in-cloud scavenging unless/until cloud droplets rain out.

             ! Interstitial aerosol collection by cloud droplets is inefficient unless aerosol or droplet is larger than about 20 microns in diameter
             ! However, removal of interstitial aerosol is likely to be much more efficient in convective updrafts than in the stratiform region
             ! For dust, interstitial collection by droplets is rare because large dust particles settle out quickly, before precipitation has a chance to act

             ! Following representation of nucleation scavenging involves assumptions:
             ! 1. Aerosol activation contributes to model in-cloud scavenging only in gridpoints/timesteps when local cloud droplets rain out
             ! 2. In-cloud re-evaporation (de-activation) does not occur (because we do not separately track aerosol mass in cloud droplets)
             ! 3. Collection of interstitial aerosol by cloud droplets does not occur (again, because we do not separately track aerosol mass in cloud droplets)
             ! 4. Only precipitation collects (thus removes) interstitial aerosol
             ! 5. Remaining interstitial aerosol is susceptible to convective redistribution
             
             ! Convective precipitation formation during timestep times specific cross-section
             ! for nucleation scavenging yields fraction of homogeneously distributed tracer
             ! in liquid cloud that is potentially nucleation scavenged this timestep
             ! Species-dependence is incorporated through spc_xsx_ncl_scv factor
             ! [frc] Potential tracer fraction removed by nucleation scavenging in convective precipitation
             ptn_frc_scv_ncl_cnv=max(0.0_r8,min(1.0_r8, &
                  q_H2O_vpr2pcp_cnv_tnd(i,k)*mpl_air(i,k)*spc_xsx_ncl_scv(m)*tm_dlt))
             
             ! [kg kg-1 s-1] Convective nucleation scavenging tracer tendency
             ! Allow dust to nucleate ice, thus allow snow to remove dust
             q_trc_tnd_scv_ncl_cnv=cld_frc_ncl_cnv*ptn_frc_scv_ncl_cnv*q_dst(i,k,m)*(1.0-cnd_frc_ice*frc_ice_scv(m))/tm_dlt
             
             ! Convective rain-collision scavenging (SeP97 p. 1021):
             ! Proportional to cross-sectional area of gridcell swept by precipitation 
             ! Collision efficiency E is probability that collision occurs within
             ! cross-sectional area swept out by path of falling raindrop.
             ! For rain and hygroscopic aerosol, collision efficiency is same
             ! as collection efficiency because sticking efficiency is ~ 1.0. 
             ! E is usually << 1 for atmospheric aerosol due to streamlines diverging around around falling raindrop
             ! Very small (D < 0.1 um) particles are only collected by Brownian diffusion
             ! Somewhere near D ~ 1.0 um particles exceed critical Stokes number and 
             ! inertial impaction collection becomes very efficient 
             ! Particles D > 10.0 um are also collected by interception
             ! In between (0.1 < D < 1.0 um) is the "Greenfield gap" where E << 1
             ! A "reasonable" raindrop radius is ~ 1 mm RoY94 p. 179
             
             ! Area flux = pcp_flx [kg m-2 s-1] * pi rds^2       [m2 H2O] * dt [s]
             ! ------- ------------   -------------- --------  
             ! dns_H2O [kg m-3 H2O]   (4/3) pi rds^3 [m3 H2O]
             
             ! = 3*rds*pcp_flx*dt [m2 H2O m-2]
             ! ----------------
             ! 4
             
             ! [frc] Potential tracer fraction removed by collision scavenging in convective precipitation
             ptn_frc_scv_cll_cnv=max(0.0_r8,min(1.0_r8, &
                  pcp_flx_cnv(i)*scv_cff_mss_avg_pcp_nrm_cnv(m)*tm_dlt/max(cld_frc_cll_cnv,cld_frc_cll_min)))
             
             ! [kg kg-1 s-1] Convective rain-collision scavenging tracer tendency
             q_trc_tnd_scv_cll_cnv=cld_frc_cll_cnv*ptn_frc_scv_cll_cnv*q_dst(i,k,m)/tm_dlt
             
             ! [kg kg-1 s-1] Total convective scavenging tracer tendency
             q_trc_tnd_scv_cnv=q_trc_tnd_scv_ncl_cnv+q_trc_tnd_scv_cll_cnv
             
             ! [frc] Nucleation scavenging fraction in convective clouds
             scv_ncl_frc_cnv=q_trc_tnd_scv_ncl_cnv/(q_trc_tnd_scv_cnv+q_trc_tnd_scv_eps) 
             
             ! Scavenging in stratiform clouds:
             
             ! Nucleation scavenging (SeP97 p. 1027):
             ! This is handled slightly differently than in convective case
             ! First we determine fraction of local condensate that precipitates
             ! We assume same fraction of tracer is removed from liquid portion of cloud
             ! Species-dependence is currently not incorporated, i.e., there is no
             ! factor to account for scavenging efficiency for in-cloud nucleation
             
             ! [kg kg-1 s-1] Stratiform nucleation scavenging tracer tendency
             ! Allow dust to nucleate ice, thus allow snow to remove dust
             q_trc_tnd_scv_ncl_str=cld_frc_ncl_str*q_H2O_cnd2pcp_frc*q_dst(i,k,m)*(1.0-cnd_frc_ice*frc_ice_scv(m))/tm_dlt
             
             ! Stratiform rain-collision scavenging (SeP97 p. 1021):
             
             ! [frc] Potential tracer fraction removed by collisions scavenging in stratiform precipitation
             ptn_frc_scv_cll_str=max(0.0_r8,min(1.0_r8, &
                  pcp_flx_str(i)*scv_cff_mss_avg_pcp_nrm_str(m)*tm_dlt/max(cld_frc_cll_str,cld_frc_cll_min)))
             
             ! [kg kg-1 s-1] Stratiform rain-collision scavenging tracer tendency
             q_trc_tnd_scv_cll_str=cld_frc_cll_str*ptn_frc_scv_cll_str*q_dst(i,k,m)/tm_dlt
             
             ! [kg kg-1 s-1] Total stratiform scavenging tracer tendency
             q_trc_tnd_scv_str=q_trc_tnd_scv_ncl_str+q_trc_tnd_scv_cll_str
             
             ! [frc] Nucleation scavenging fraction in stratiform clouds
             scv_ncl_frc_str=q_trc_tnd_scv_ncl_str/(q_trc_tnd_scv_str+q_trc_tnd_scv_eps) 
             
             ! Ensure sum of convective and stratiform scavenging tendencies does not exceed amount present
             scv_tnd_fct=q_dst(i,k,m)/max(tm_dlt*(q_trc_tnd_scv_cnv+q_trc_tnd_scv_str),q_trc_tnd_scv_eps) ! [frc]
             if (scv_tnd_fct < 1.0_r8) then
                ! Total stratiform scavenging tracer tendency
                q_trc_tnd_scv_str=q_trc_tnd_scv_str*scv_tnd_fct ! [kg kg-1 s-1]
                ! Total convective scavenging tracer tendency
                q_trc_tnd_scv_cnv=q_trc_tnd_scv_cnv*scv_tnd_fct ! [kg kg-1 s-1]
                ! It is now "safe" to use these tendencies to adjust tracer mixing ratio
                ! However, it is only "safe" to machine precision since order of operations can still cause problems
             endif            ! endif
             
             ! Nomenclature: 
             ! Scavenging tendency refers to removal by nucleation and collision scavenging [positive definite]
             ! Evaporation tendency refers to source by evaporation of stratiform precipitation from above [positive definite]
             ! Wet deposition tendency is evaporation tendency (source) minus scavenging tendency (sink) [positive when source > sink]
             ! [kg kg-1 s-1] Scavenging tendency
             q_dst_tnd_pcp(i,k,m)=q_trc_tnd_scv_str+q_trc_tnd_scv_cnv 
             
             ! [kg kg-1 s-1] Evaporation tendency
             q_dst_tnd_evp(i,k,m)=pcp_str_frc_evp*flx_mss_scv_abv_str(i,m)/mpl_air(i,k)
             
             ! Tracer not removed within cloud is interstitial and subject to convective transport
             ! fxm: restrict to truly convective fraction and use entrainment factor beta to allow for "vacuuming"
             frc_trc_trn_cnv_ptn(i,k,m)= & ! [frc] Interstitial tracer fraction
                  1.0-max(0.0_r8,min(1.0_r8, &
                  q_dst_tnd_pcp(i,k,m)*tm_dlt/max(cld_vlm(i,k)*q_dst(i,k,m),q_dst_cld_min)))
             
             ! Diagnose nucleation scavenging tendency for possible archival (superfluous, could be deleted)
             q_dst_tnd_ncl(i,k,m)= & ! [kg kg-1 s-1] Nucleation scavenging tendency
                  q_trc_tnd_scv_cnv*scv_ncl_frc_cnv+q_trc_tnd_scv_str*scv_ncl_frc_str
             
             ! Diagnose wet deposition (evaporation minus scavenging) tendency for possible archival (superfluous, could be deleted)
             q_dst_tnd_wet(i,k,m)=q_dst_tnd_evp(i,k,m)-q_dst_tnd_pcp(i,k,m) ! [kg kg-1 s-1]
             
             ! Tendencies are known, so adjust mixing ratios
             ! max() clause is necessary because result cannot be guarranteed positive otherwise
             q_dst(i,k,m)=max(0.0_r8,q_dst(i,k,m)+tm_dlt*q_dst_tnd_wet(i,k,m)) ! [kg kg-1] Dust mixing ratio
             
             ! Accumulate aerosol-dependent vertical integrals
             ! All tracer in stratiform rain may potentially be evaporated at every level
             flx_mss_scv_abv_str(i,m)=flx_mss_scv_abv_str(i,m)*(1.0_r8-pcp_str_frc_evp)+q_trc_tnd_scv_str*mpl_air(i,k) ! [kg m-2 s-1] Tracer mass flux scavenged from above, stratiform
             ! Tracer in convective rain increases monotonically on way down--convective evaporation is not accounted for
             flx_mss_scv_abv_cnv(i,m)=flx_mss_scv_abv_cnv(i,m)+q_trc_tnd_scv_cnv*mpl_air(i,k) ! [kg m-2 s-1] Tracer mass flux scavenged from above, convective
             
          end do              ! end loop over cst
          
          ! Accumulate aerosol-independent vertical integrals
          ! Precipitation falling into cell below is precipitation that entered this cell from above, 
          ! plus local source of precipitation, minus what was evaporated. 
          pcp_flx_str(i)=pcp_flx_str(i)+(q_H2O_cnd2pcp_tnd(i,k)-q_H2O_pcp2vpr_tnd(i,k))*mpl_air(i,k) ! [kg m-2 s-1] Stratiform precipitation flux
          ! Evaporation of convective precipitation is not allowed here
          pcp_flx_cnv(i)=pcp_flx_cnv(i)+q_H2O_vpr2pcp_cnv_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1] Convective precipitation flux
          
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
#ifdef DST_DBG
    ! Sanity check
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             if (q_dst(i,k,m) < 0.0_r8 .or. q_dst(i,k,m) > q_dst_mxm) then
                write(6,'(a)') 'dst: ERROR dst_dps_wet() reports mixing ratio too high or too low:'
                write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1,a,2(a,es8.1,a))') &
                     'dst_pcp: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m),' kg kg-1' &
                     ,', q_dst_tnd_evp = ',q_dst_tnd_evp(i,k,m),' kg kg-1 s-1' &
                     ,', q_dst_tnd_pcp = ',q_dst_tnd_pcp(i,k,m),' kg kg-1 s-1'
                stop
             endif            ! endif err
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
#endif /* not DST_DBG */
    
    ! Diagnostic totals
    do i=1,plon
       pcp_flx_sfc(i)=pcp_flx_str(i)+pcp_flx_cnv(i) ! [kg m-2 s-1] Total precipitation reaching surface
       ! Convert [kg m-2 s-1] to [m s-1]
       prect2(i)=pcp_flx_sfc(i)*0.001 ! [m s-1] Total precipitation reaching surface (same as PRECT)
    end do                    ! end loop over lon
    do m=1,dst_nbr
       do i=1,plon
          ! Tracer reaching surface in precipitation
          flx_mss_pcp_sfc(i,m)=flx_mss_scv_abv_str(i,m)+flx_mss_scv_abv_cnv(i,m) ! [kg m-2 s-1]
       end do                 ! end loop over cst
    end do                    ! end loop over lon
    
    ! Integrate over all size categories and levels for diagnostic output
    call dst_add_lon_lev(q_dst_tnd_evp,q_dst_tnd_evp_ttl)
    call dst_add_lon_lev(q_dst_tnd_pcp,q_dst_tnd_pcp_ttl)
    call dst_add_lon(flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl)
    
#ifdef BXM
    ! Update netCDF file
    fl_out='aer.nc'           ! [sng] Name of netCDF output file
    call ftn_strnul(fl_out)
    call dpswet2nc(             &
         cld_frc,             & ! I [frc] Local total cloud fraction
         cld_frc_cnv,         & ! I [frc] Convective cloud fraction
         cld_vlm,             & ! I [frc] Cloud volume
         fl_out,              & ! I [sng] Name of netCDF output file
         flx_mss_pcp_sfc,     & ! I [kg m-2 s-1] Dust reaching surface in precipitation
         flx_mss_pcp_sfc_ttl, & ! I [kg m-2 s-1] Total dust reaching surface in precipitation
         frc_trc_trn_cnv_ptn,         & ! I [frc] Interstitial tracer fraction
         lat_idx,             & ! I [idx] Model latitude index
         pcp_flx_sfc,         & ! I [kg m-2 s-1] Total precipitation reaching surface
         q_H2O_cnd,           & ! I [kg kg-1] Condensed H2O mixing ratio
         q_H2O_cnd2pcp_tnd,   & ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
         q_H2O_pcp2vpr_tnd,   & ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
         q_H2O_vpr2pcp_cnv_tnd, & ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
         q_dst_tnd_evp,       & ! I [kg kg-1 s-1] Evaporation tendency
         q_dst_tnd_evp_ttl,   & ! I [kg kg-1 s-1] Total evaporation tendency
         q_dst_tnd_ncl,       & ! I [kg kg-1 s-1] Nucleation scavenging tendency
         q_dst_tnd_pcp,       & ! I [kg kg-1 s-1] Scavenging tendency
         q_dst_tnd_pcp_ttl,   & ! I [kg kg-1 s-1] Total scavenging tendency
         q_dst_tnd_wet,       & ! I [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
         scv_cff_mss_avg_pcp_nrm_cnv, & ! I [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, convective
         scv_cff_mss_avg_pcp_nrm_str, & ! I [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, stratiform
         spc_xsx_ncl_scv)     ! I [m2 kg-1] Specific cross section for nucleation scavenging
#endif /* not BXM */
#ifdef CCM
    call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,ncol,lchnk)
    call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,ncol,lchnk)
    call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,ncol,lchnk)
    do m=1,dst_nbr
       call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),ncol,lchnk)
    end do                    ! end loop over cst
#else /* not CCM */
    call outfld('PRECT2',prect2,plond,lat_idx,obuf)
    call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,obuf)
    call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,obuf)
    call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,obuf)
    do m=1,dst_nbr
       call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,obuf)
    end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
    call bdg_aal('dst_sf_pcp',lat_idx,flx_mss_pcp_sfc_ttl)
    call bdg_gam_wet('dst_ss_pcp',lat_idx,prs_dlt,q_dst_tnd_pcp_ttl)
    call bdg_gam_wet('dst_ss_evp',lat_idx,prs_dlt,q_dst_tnd_evp_ttl)
#endif /* not DST_MSS_BDG */
    return
  end subroutine dst_dps_wet                       ! end dst_dps_wet()
  
  subroutine cld_dgn( & ! [sbr] Cloud diagnostics
       cld_frc,             & ! I [frc] Local total cloud fraction
       cld_vlm,             & ! O [frc] Cloud volume
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_cnd2pcp_tnd,   & ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
       q_H2O_cnd_tnd,       & ! I [kg kg-1 s-1] Net H2O condensate formation tendency
       q_H2O_pcp_lqd,       & ! O [kg kg-1] Rainwater mixing ratio
       q_H2O_pcp2vpr_tnd,   & ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
       q_H2O_vpr2pcp_cnv_tnd, & ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
       tpt_mdp)             ! I [K] Temperature
    ! Purpose: Diagnose statistics of cloud properties
    ! Same as MATCH:src/cloud/clddiag()
    ! cld_dgn() is called by MATCH:physlic() and CCM:aphys()
    use dstblm,only:dns_H2O_lqd_std,tpt_frz_pnt ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use dstcst ! [mdl] Physical constants for dust routines
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    real(r8),parameter::rds_pcp=1000.0 ! [m] Radius of precipitation RoY94 p. 179
    real(r8),parameter::q_H2O_pcp_lqd_min=1.0e-14 ! [kg kg-1] Minimum rainwater mixing ratio
    real(r8),parameter::pcp_flx_sfc_min=1.0e-30 ! [kg m-2 s-1] Minimum precipitation rate to prevent divide by zero
    ! Input
    real(r8),intent(in)::cld_frc(plond,plevp) ! I [frc] Local total cloud fraction
    real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure
    real(r8),intent(in)::q_H2O_cnd2pcp_tnd(plond,plev) ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
    real(r8),intent(in)::q_H2O_cnd_tnd(plond,plev) ! I [kg kg-1 s-1] Net H2O condensate formation tendency
    real(r8),intent(in)::q_H2O_pcp2vpr_tnd(plond,plev) ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
    real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
    real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
    ! Output
    real(r8),intent(out)::cld_vlm(plond,plev) ! O [frc] Cloud volume
    real(r8),intent(out)::q_H2O_pcp_lqd(plond,plev) ! O [kg kg-1] Rainwater mixing ratio
    ! Local
    integer i                 ! [idx] Longitude index
    integer k                 ! [idx] Level index
    real(r8) cld_frc_pcp_wgt(plond) ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
    real(r8) dns_mdp              ! [kg m-3] Midlayer density
    real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Precipitation rate 
    real(r8) pcp_flx_pst_sfc(plond) ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes
    real(r8) pcp_flx_lcl          ! [kg m-2 s-1] Local net precipitation production rate (+'ve when formation exceeds evaporation)
    real(r8) pcp_flx_lcl_pst      ! [kg m-2 s-1] Positive definite local net precipitation production rate
    real(r8) vlc_grv_fct          ! [kg1/2 m-1/2 s-1] Factor in terminal velocity formula SeP97 p. 467 (8.44)
    real(r8) vlc_grv_pcp          ! [m s-1] Fallspeed of precipitation
    
    ! Collision scavenging depends on horizontal fraction of local gridcell through which rain from above falls
    ! This ought to be true for collision scavenging by both stratiform and convective rain
    ! Nucleation scavenging should depend solely on local cloud fraction
    ! This ought to be true for nucleation scavenging by both stratiform and convective clouds
    ! Note that pre-CCM4, no information about is known about vertical correlation functions between stratiform cloud at different levels
    
    ! In developing the SO4 wet deposition scheme, PJR and MKB settled on using a single cloud fraction,
    ! called cld_vlm, for all wet scavenging processes (nucleation and collision)
    ! They define this factor as the maximum of local cloud amount and
    ! 
    ! sum above of (cloud*positive precipitation production)   sum real precipitation from above
    ! ------------------------------------------------------ x --------------------------------------------
    ! sum above of (      positive precipitation production)   sum positive precipitation from above
    ! 
    ! According to this definition, cld_vlm 
    ! 1. increases monotonically from the upper troposphere to the lower troposphere
    ! 2. is unity in and everywhere beneath the first socked-in layer
    
    vlc_grv_fct=1.94_r8*2.13_r8*sqrt(dns_H2O_lqd_std*rds_pcp*grv_sfc*2.7e-4_r8) ! [kg1/2 m-1/2 s-1] Factor in terminal velocity formula SeP97 p. 467 (8.44)
    do i=1,plon
       pcp_flx_sfc(i)=0.0_r8     ! [kg m-2 s-1] Precipitation rate
       cld_frc_pcp_wgt(i)=0.0_r8 ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
       pcp_flx_pst_sfc(i)=pcp_flx_sfc_min ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes
    end do                    ! end loop over lon
    
    do k=1,plev
       do i=1,plon
          cld_vlm(i,k)=       & ! [frc] Cloud volume
               max(min(1.0_r8, &
               cld_frc_pcp_wgt(i)/pcp_flx_pst_sfc(i) &
               )*pcp_flx_sfc(i)/pcp_flx_pst_sfc(i), &
               cld_frc(i,k))
          pcp_flx_lcl=        & ! [kg m-2 s-1] Local net precipitation production rate (+'ve when formation exceeds evaporation)
               prs_dlt(i,k)*grv_sfc_rcp* &
               (q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)-q_H2O_pcp2vpr_tnd(i,k))
          pcp_flx_sfc(i)=pcp_flx_sfc(i)+pcp_flx_lcl ! [kg m-2 s-1] Precipitation rate
          
          ! Positive definite version of pcp_flx_lcl used to weight contribution of layer cloud to precipitation-weighted cloud fraction
          pcp_flx_lcl_pst=max(pcp_flx_lcl,0.0_r8) ! [kg m-2 s-1] Positive definite local net precipitation production rate
          cld_frc_pcp_wgt(i)=cld_frc_pcp_wgt(i)+cld_frc(i,k)*pcp_flx_lcl_pst ! [kg m-2 s-1] Integral of precipitation-weighted cloud fraction
          pcp_flx_pst_sfc(i)=pcp_flx_pst_sfc(i)+pcp_flx_lcl_pst ! [kg m-2 s-1] Sum of positive definite precipitation production fluxes
          
          q_H2O_pcp_lqd(i,k)=0.0_r8 ! [kg kg-1] Rainwater mixing ratio
          if (tpt_mdp(i,k) > tpt_frz_pnt) then
             ! fxm: Should use virtual temperature in density calculation
             dns_mdp=prs_mdp(i,k)/(gas_cst_dry_air*tpt_mdp(i,k)) ! [kg m-3] Midlayer density
             vlc_grv_pcp=vlc_grv_fct/sqrt(dns_mdp) ! [m s-1] Fallspeed of precipitation
             q_H2O_pcp_lqd(i,k)=pcp_flx_sfc(i)/(dns_mdp*vlc_grv_pcp) ! [kg kg-1] Rainwater mixing ratio
             if (q_H2O_pcp_lqd(i,k) < q_H2O_pcp_lqd_min) q_H2O_pcp_lqd(i,k)=0.0_r8 ! [kg kg-1] Rainwater mixing ratio
          endif               ! endif not freezing
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    return
  end subroutine cld_dgn
  
  subroutine dst_dps_wet_old(lchnk,ncol,obuf, &
       cld_frc,             & ! I [frc] Local total cloud fraction
       cld_frc_cnv,         & ! I [frc] Convective cloud fraction
       frc_trc_trn_cnv_ptn, & ! O [frc] Interstitial tracer fraction
       lat_idx,             & ! I [idx] Model latitude index
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       q_H2O_cnd2pcp_tnd,   & ! I/O [kg kg-1 s-1] Condensed H2O to precipitation tendency
       q_H2O_pcp2vpr_tnd,   & ! I/O [kg kg-1 s-1] H2O precipitation to vapor tendency
       q_H2O_vpr2pcp_cnv_tnd, & ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
       q_dst,               & ! I/O [kg kg-1] Dust mixing ratio
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp)             ! I [K] Temperature
    ! Simulate wet deposition processes 
    ! dst_dps_wet_old() is called by CCM:physics/aphys(), MATCH:src/physlic()
    ! Current processes considered are:
    ! 1. Nucleation scavenging
    ! 2. Evaporation of dust in stratiform rain evaporation
    ! Algorithm based on Tegen and Fung 1995
    ! Assumes dust is homogeneously distributed throughout a gridbox
    ! Scavenging, however, should only be allowed to occur in that fraction of a gridbox which is covered by clouds
    ! fxm: 20000525 Slight differences in diagnosed precipitation between this routine and host model
    ! might be due to applying stratiform evaporation tendency to both stratiform and convective rain
    use dstbdg,only:bdg_aal,bdg_gam_wet ! [mdl] Mass budget diagnostics
    use dstcst ! [mdl] Physical constants for dust routines
    use dstdbg ! [mdl] Debugging information for dust model
    use dstgrd,only:dst_nbr ! [mdl] Dust grid sizes
    use dstmssutl,only:dst_add_lon_lev,dst_add_lon ! [mdl] Mass budget utilities
    use dstnm ! [mdl] Nomenclature for outfld()
    use dstscv,only:z_scv ! [mdl] Aerosol scavenging properties
#ifndef BXM
    use histout,only:outfld ! [mdl] History/archiving
#endif /* BXM */
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    real(r8),parameter::q_dst_cld_min=1.0e-30 ! [kg kg-1] Epsilon to prevent divide by zero for interstitial fraction
    ! Input
    integer,intent(in)::lat_idx ! Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::obuf(*)  ! Output buffer
    real(r8),intent(in)::cld_frc(plond,plevp) ! [frc] Local total cloud fraction
    real(r8),intent(in)::cld_frc_cnv(plond,plev) ! [frc] Convective cloud fraction
    real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
    real(r8),intent(in)::prs_dlt(plond,plev) ! [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond,plev) ! [Pa] Midlayer pressure
    real(r8),intent(in)::tpt_mdp(plond,plev) ! [K] Temperature
    real(r8),intent(in)::tm_adj   ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
    ! Output
    real(r8),intent(out)::frc_trc_trn_cnv_ptn(plond,plev,dst_nbr) ! O [frc] Interstitial tracer fraction
    ! Input/Output
    real(r8),intent(inout)::q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
    real(r8),intent(inout)::q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
    real(r8),intent(inout)::q_dst(plond,plev,dst_nbr) ! [kg dust/kg moist air]
    ! Local Output
    ! Local
    integer i                 ! [idx] lon index
    integer k                 ! [idx] lev index
    integer m                 ! [idx] cst index 
    real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
    real(r8) tm_dlt               ! [s] Wet deposition timestep
    real(r8) q_H2O_wtr2pcp_ttl_tnd(plond,plev) ! [kg kg-1 s-1] Local water to total precipitation tendency
    ! GCM diagnostics
    real(r8) cld_frc_eff          ! [frc] Fraction of dust removable by precip
    real(r8) dst2pcp_flx          ! [kg m-2 s-1] Local removal of dust by precip (q_H2O_wtr2pcp_ttl_tnd*q_dst*z_scv)
    real(r8) flx_mss_evp          ! [kg m-2 s-1] Local source of dust due to evaporation of precipitation from above (pcpfrc_evp*flx_mss_pcp_sfc)
    real(r8) flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
    real(r8) flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
    real(r8) pcp2vpr_flx          ! [kg m-2 s-1] Local conversion of precip to vapor (q_H2O_pcp2vpr_tnd*mpl_air)
    real(r8) pcp_flx_sfc(plond)   ! [kg m-2 s-1] Total precipitation reaching surface
    real(r8) prect2(plond)        ! [m s-1] Total precipitation reaching surface (diagnostic, should be same as PRECT)
    real(r8) pcp_frc_evp          ! [frc] Precipitation from above that evaporates locally (pcp2vpr_flx/pcp_flx_sfc)
    real(r8) q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
    real(r8) q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
    real(r8) q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
    real(r8) q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
    real(r8) wtr2pcp_flx          ! [kg m-2 s-1] Local conversion of water to precip (q_H2O_wtr2pcp_ttl_tnd*mpl_air)
    
    ! Main Code
    ! Timesplit if desired
    tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
    
    ! Initialize fluxes and tendencies
    q_dst_tnd_pcp(:,:,:)=0.0_r8         ! [kg kg-1 s-1]
    q_dst_tnd_evp(:,:,:)=0.0_r8         ! [kg kg-1 s-1]
    flx_mss_pcp_sfc(:,:)=0.0_r8       ! [kg m-2 s-1]
    
    ! Sanity check
    do k=1,plev
       do i=1,plon
          ! q_H2O_pcp2vpr_tnd, and q_H2O_cnd2pcp_tnd may be < 0.0, though I do not understand why
          q_H2O_pcp2vpr_tnd(i,k)=max(q_H2O_pcp2vpr_tnd(i,k),0.0_r8) ! [kg kg-1 s-1] 
          q_H2O_cnd2pcp_tnd(i,k)=max(q_H2O_cnd2pcp_tnd(i,k),0.0_r8) ! [kg kg-1 s-1] 
          ! if (q_H2O_pcp2vpr_tnd(i,k) < 0.0_r8) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_pcp2vpr_tnd(',i,',',k,') = ',q_H2O_pcp2vpr_tnd(i,k)
          ! if (q_H2O_cnd2pcp_tnd(i,k) < 0.0_r8) write(6,100) 'dst: dst_dps_wet: lat = ',lat_idx,' q_H2O_cnd2pcp_tnd(',i,',',k,') = ',q_H2O_cnd2pcp_tnd(i,k)
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    ! 100 format(a,i2,a,i3,a,i2,a,es8.1) 
    
    ! Compute necessary derived fields
    do k=1,plev
       do i=1,plon
          ! Mass of air currently in gridbox
          mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
          ! Total precipitation formation
          q_H2O_wtr2pcp_ttl_tnd(i,k)= & ! [kg kg-1 s-1] H2O vapor + condensate conversion to precipitation tendency
               q_H2O_cnd2pcp_tnd(i,k)+q_H2O_vpr2pcp_cnv_tnd(i,k)
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    ! Scavenging of dust by precipitation (wet deposition)
    ! Scavenging ratio relates precip tendency to dust removal tendency
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             ! Compute gridpoint dust tendency for wet deposition
             ! Currently, application of cloud fraction to dust in rain is rather haphazard
             ! Needs to be consistent with physical picture of a cloud and rainshaft
             ! fxm: Dust should be last thing to evaporate in a raindrop since it is a CCN
             ! Thus, dust should only be evaporated from virga
             ! Using cld_frc_eff = cld_frc(i,k) (dynamic cloud fraction) results in lifetime of clay of 30 days
             ! In CCM, using cld_frc_eff = 0.85 results in precipitation lifetime of clay of 12-13 days for dynamic source scheme, 19-20 days for TOMS scheme 
             ! In MATCH, using cld_frc_eff = 0.85 results in precipitation lifetime of clay of 4--5 days for dynamic source scheme, ??? days for TOMS scheme 
             ! For a given model, The DSS/TOMS timescales differ because TOMS scheme injects more dust into subtropical easterlies which undergo fewer precipitation events (??)
             ! Actual lifetime of clay should be 14--21 days
             ! fxm: should use cld_frc_eff=f(cld_frc), should be fixed by PJR wet dps routine
             cld_frc_eff=0.85_r8 ! [frc] Until dstmch18
             ! cld_frc_eff=cld_frc(i,k) ! [frc] As of dstmch18
             q_dst_tnd_pcp(i,k,m)= & ! [kg kg-1 s-1]
                  cld_frc_eff*z_scv(m)*q_dst(i,k,m)*q_H2O_wtr2pcp_ttl_tnd(i,k)
             ! Do not scavenge more than available dust
             if (tm_dlt*q_dst_tnd_pcp(i,k,m) < q_dst(i,k,m)) then
                ! Tendencies are known, so adjust mixing ratios
                q_dst(i,k,m)=q_dst(i,k,m)- &
                     tm_dlt*q_dst_tnd_pcp(i,k,m) ! [kg kg-1]
             else             ! endif physical
                ! Only bother printing a warning when it concerns significant amounts dust
                if (q_dst(i,k,m) > q_dst_sgn) write (6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
                     'dst: pcp_fxr: lat = ',lat_idx,' q_dst_tnd_pcp(',i,',',k,',',m,') = ',q_dst_tnd_pcp(i,k,m)
                ! Reset to physical value
                q_dst_tnd_pcp(i,k,m)=q_dst(i,k,m)/tm_dlt ! [kg kg-1 s-1]
                q_dst(i,k,m)=0.0_r8 ! [kg kg-1]
             endif            ! endelse unphysical
             
             ! Tracer not removed within cloud is interstitial and subject to convective transport
             ! fxm: multiply cld_frc_cnv by entrainment factor beta to allow for "vacuuming"
             frc_trc_trn_cnv_ptn(i,k,m)= & ! [frc] Interstitial tracer fraction
                  1.0-max(0.0_r8,min(1.0_r8, &
                  q_dst_tnd_pcp(i,k,m)*tm_dlt/max(cld_frc_cnv(i,k)*q_dst(i,k,m),q_dst_cld_min)))
             
#ifdef DST_DBG
             if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
                  'dst_pcp: lat = ',lat_idx,' qdst(',i,',',k,',',m,') = ',q_dst(i,k,m)
#endif /* not DST_DBG */
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    ! Initialize surface precipitation with top level precipitation
    ! Assume no evaporation in top level
    do i=1,plon
       pcp_flx_sfc(i)=q_H2O_wtr2pcp_ttl_tnd(i,1)*mpl_air(i,1) ! [kg m-2 s-1]
    end do                    ! end loop over lon
    do m=1,dst_nbr
       do i=1,plon
          flx_mss_pcp_sfc(i,m)=q_dst_tnd_pcp(i,1,m)*mpl_air(i,1) ! [kg m-2 s-1]
          q_dst_tnd_evp(i,1,m)=0.0_r8 ! [kg kg-1 s-1]
       end do                 ! end loop over lon
    end do                    ! end loop over cst
    ! Invert loops (lon is inner dimension) process from top down on each longitude
    ! Surface fluxes accumulate as local washout and evaporation modify precipitation on its way down
    ! Algorithm is based on JJH's cond() routine
    do k=2,plev               ! NB: level starts with 2 (may generate harmless warning when plev=1, i.e., in BXM)
       do i=1,plon
          ! Convert local condensation and evaporation tendencies to fluxes
          wtr2pcp_flx=q_H2O_wtr2pcp_ttl_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1]
          pcp2vpr_flx=q_H2O_pcp2vpr_tnd(i,k)*mpl_air(i,k) ! [kg m-2 s-1]
          
          ! Compute fraction of rain from above that evaporates in this gridcell
          ! This fraction will be applied to dust in rain from above to compute
          ! dust released into this cell due to evaporation.
          if (pcp_flx_sfc(i) > 0.0_r8) then
             ! Use safety limiter because pcp_frc_evp often equals, e.g., 1.0000000000000007
             pcp_frc_evp=min(pcp2vpr_flx/pcp_flx_sfc(i),1.0_r8) ! [frc]
#ifdef DST_DBG
             if (pcp_frc_evp > 1.0.or.pcp_frc_evp < 0.0_r8) write(6,'(a,i2,a,i3,a,i2,a,es8.1)')  &
                  'dst: pcp_evp: lat = ',lat_idx,' pcp_frc_evp(',i,',',k,') = ',pcp_frc_evp
#endif /* not DST_DBG */
          else
             pcp_frc_evp=0.0_r8  ! [frc]
          endif               ! endif pcp > 0
          
          ! Precipitation falling into cell below is precipitation entering cell from above, 
          ! plus local source of precipitation, minus what was evaporated. 
          pcp_flx_sfc(i)=max( & ! [kg m-2 s-1]
               pcp_flx_sfc(i)+wtr2pcp_flx- &
               pcp2vpr_flx &
               ,0.0_r8)
          do m=1,dst_nbr
             dst2pcp_flx=q_dst_tnd_pcp(i,k,m)*mpl_air(i,k) ! [kg m-2 s-1]
             ! Evaporate dust into cell based on its concentration in rain from above
             flx_mss_evp=pcp_frc_evp*flx_mss_pcp_sfc(i,m) ! [kg m-2 s-1]
             q_dst_tnd_evp(i,k,m)=flx_mss_evp/mpl_air(i,k) ! [kg kg-1 s-1]
             ! Do not let more dust be released by evaporation than is falling through gridcell
             if (q_dst_tnd_evp(i,k,m) > flx_mss_pcp_sfc(i,m)/mpl_air(i,k)) then
                if (flx_mss_evp > flx_mss_sgn) write (6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1,a,es8.1)')  &
                     'dst: evp_fxr: lat = ',lat_idx,' flx_mss_evp(',i,',',k,',',m,') = ',flx_mss_evp, &
                     ' flx_mss_pcp_sfc = ',flx_mss_pcp_sfc(i,m)
                ! Reset to physical value
                q_dst_tnd_evp(i,k,m)=flx_mss_pcp_sfc(i,m)/(mpl_air(i,k)*tm_dlt) ! [kg kg-1 s-1]
             endif            ! endif unphysical
             ! Tendencies are known, so adjust mixing ratios
             q_dst(i,k,m)=q_dst(i,k,m)+ &
                  tm_dlt*q_dst_tnd_evp(i,k,m) ! [kg kg-1]
#ifdef DST_DBG
             if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
                  'dst_evp: lat = ',lat_idx,' qdst(',i,',',k,',',m,') = ',q_dst(i,k,m)
#endif /* not DST_DBG */
             
             ! Wet dust falling into next cell down is wet dust entering cell from above, 
             ! plus local source of wet dust, minus what was evaporated. 
             flx_mss_pcp_sfc(i,m)=max( &
                  flx_mss_pcp_sfc(i,m)+dst2pcp_flx- &
                  flx_mss_evp &
                  ,0.0_r8)
          end do              ! end loop over cst
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    ! Get rid of unphysical values
    ! 20000415: Better to allow qneg() routines to perform these checks and fixes
    if (.false.) then
       do m=1,dst_nbr
          do k=1,plev
             do i=1,plon
                q_dst(i,k,m)=max(q_dst(i,k,m),0.0_r8) ! [kg kg-1]
                if (q_dst(i,k,m) < q_dst_sgn) q_dst(i,k,m)=0.0_r8 ! [kg kg-1]
             end do           ! end loop over lon
          end do              ! end loop over lev
       end do                 ! end loop over cst
    endif                     ! endif false
    
    ! Integrate over all size categories and levels for diagnostic output
    call dst_add_lon_lev(q_dst_tnd_evp,q_dst_tnd_evp_ttl)
    call dst_add_lon_lev(q_dst_tnd_pcp,q_dst_tnd_pcp_ttl)
    call dst_add_lon(flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl)
    
#ifdef CCM
    call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,ncol,lchnk)
    call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,ncol,lchnk)
    call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,ncol,lchnk)
    do m=1,dst_nbr
       call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),ncol,lchnk)
    end do                    ! end loop over cst
#else /* not CCM */
    do i=1,plon
    ! Convert [kg m-2 s-1] to [m s-1]
       prect2(i)=pcp_flx_sfc(i)*0.001 ! [m s-1] Total precipitation reaching surface (same as PRECT)
    end do                    ! end loop over lon
    call outfld('PRECT2',prect2,plond,lat_idx,obuf)
    call outfld('DSTSFPCP',flx_mss_pcp_sfc_ttl,plond,lat_idx,obuf)
    call outfld('DSTSSPCP',q_dst_tnd_pcp_ttl,plond,lat_idx,obuf)
    call outfld('DSTSSEVP',q_dst_tnd_evp_ttl,plond,lat_idx,obuf)
    do m=1,dst_nbr
       call outfld(flx_mss_pcp_sfc_nm(m),flx_mss_pcp_sfc(1,m),plond,lat_idx,obuf)
    end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
    call bdg_aal('dst_sf_pcp',lat_idx,flx_mss_pcp_sfc_ttl)
    call bdg_gam_wet('dst_ss_pcp',lat_idx,prs_dlt,q_dst_tnd_pcp_ttl)
    call bdg_gam_wet('dst_ss_evp',lat_idx,prs_dlt,q_dst_tnd_evp_ttl)
#endif /* not DST_MSS_BDG */
    return
  end subroutine dst_dps_wet_old                       ! end dst_dps_wet_old()
  
  subroutine dpswet2nc(             &
       cld_frc,             & ! I [frc] Local total cloud fraction
       cld_frc_cnv,         & ! I [frc] Convective cloud fraction
       cld_vlm,             & ! I [frc] Cloud volume
       fl_out,              & ! I [sng] Name of netCDF output file
       flx_mss_pcp_sfc,     & ! I [kg m-2 s-1] Dust reaching surface in precipitation
       flx_mss_pcp_sfc_ttl, & ! I [kg m-2 s-1] Total dust reaching surface in precipitation
       frc_trc_trn_cnv_ptn, & ! I [frc] Interstitial tracer fraction
       lat_idx,             & ! I [idx] Model latitude index
       pcp_flx_sfc,         & ! I [kg m-2 s-1] Total precipitation reaching surface
       q_H2O_cnd,           & ! I [kg kg-1] Condensed H2O mixing ratio
       q_H2O_cnd2pcp_tnd,   & ! I [kg kg-1 s-1] Condensed H2O to precipitation tendency
       q_H2O_pcp2vpr_tnd,   & ! I [kg kg-1 s-1] H2O precipitation to vapor tendency
       q_H2O_vpr2pcp_cnv_tnd, & ! I [kg kg-1 s-1] H2O vapor to convective precipitation tendency
       q_dst_tnd_evp,       & ! I [kg kg-1 s-1] Evaporation tendency
       q_dst_tnd_evp_ttl,   & ! I [kg kg-1 s-1] Total evaporation tendency
       q_dst_tnd_ncl,       & ! I [kg kg-1 s-1] Nucleation scavenging tendency
       q_dst_tnd_pcp,       & ! I [kg kg-1 s-1] Scavenging tendency
       q_dst_tnd_pcp_ttl,   & ! I [kg kg-1 s-1] Total scavenging tendency
       q_dst_tnd_wet,       & ! I [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
       scv_cff_mss_avg_pcp_nrm_cnv, & ! I [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, convective
       scv_cff_mss_avg_pcp_nrm_str, & ! I [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, stratiform
       spc_xsx_ncl_scv)     ! I [m2 kg-1] Specific cross section for nucleation scavenging
    ! Purpose: Output time-varying wet deposition properties to netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstctl ! [mdl] Control variables, routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='dpswet2nc' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    real(r8),intent(in)::cld_frc(plond,plevp) ! [frc] Cloud fraction
    real(r8),intent(in)::cld_frc_cnv(plond,plev) ! [frc] Convective cloud fraction
    real(r8),intent(in)::cld_vlm(plond,plev) ! [frc] Cloud volume
    real(r8),intent(in)::q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
    real(r8),intent(in)::q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
    real(r8),intent(in)::q_H2O_cnd(plond,plev) ! [kg kg-1] Condensed H2O mixing ratio
    real(r8),intent(in)::q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
    real(r8),intent(in)::spc_xsx_ncl_scv(dst_nbr) ! [m2 kg-1] Specific cross section for nucleation scavenging
    real(r8),intent(in)::scv_cff_mss_avg_pcp_nrm_cnv(dst_nbr) ! [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, convective
    real(r8),intent(in)::scv_cff_mss_avg_pcp_nrm_str(dst_nbr) ! [m2 kg-1] Mass mean scavenging coefficient, precipitation normalized, stratiform
    real(r8),intent(in)::flx_mss_pcp_sfc(plond,dst_nbr) ! [kg m-2 s-1] Dust reaching surface in precipitation
    real(r8),intent(in)::flx_mss_pcp_sfc_ttl(plond) ! [kg m-2 s-1] Total dust reaching surface in precipitation
    real(r8),intent(in)::frc_trc_trn_cnv_ptn(plond,plev,dst_nbr) ! O [frc] Interstitial tracer fraction
    real(r8),intent(in)::pcp_flx_sfc(plond) ! [kg m-2 s-1] Total precipitation reaching surface
    real(r8),intent(in)::q_dst_tnd_evp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Evaporation tendency
    real(r8),intent(in)::q_dst_tnd_evp_ttl(plond,plev) ! [kg kg-1 s-1] Total evaporation tendency
    real(r8),intent(in)::q_dst_tnd_ncl(plond,plev,dst_nbr) ! O [kg kg-1 s-1] Nucleation scavenging tendency
    real(r8),intent(in)::q_dst_tnd_pcp(plond,plev,dst_nbr) ! [kg kg-1 s-1] Scavenging tendency
    real(r8),intent(in)::q_dst_tnd_pcp_ttl(plond,plev) ! [kg kg-1 s-1] Total scavenging tendency
    real(r8),intent(in)::q_dst_tnd_wet(plond,plev,dst_nbr) ! [kg kg-1 s-1] Wet deposition (evaporation minus scavenging) tendency
    ! Output
    ! Local
    ! File metadata and dimension IDs
    integer cnt_lon_lev_sz_time(4) ! Count array
    integer cnt_lon_levp_sz_time(4) ! Count array
    integer cnt_lon_lev_time(3) ! Count array
    integer cnt_lon_levp_time(3) ! Count array
    integer cnt_lon_sz_time(3) ! Count array
    integer cnt_lon_time(2)   ! Count array
    integer dim_lon_lev_sz_time(4) ! [enm] Dimension IDs
    integer dim_lon_levp_sz_time(4) ! [enm] Dimension IDs
    integer dim_lon_lev_time(3) ! [enm] Dimension IDs
    integer dim_lon_levp_time(3) ! [enm] Dimension IDs
    integer dim_lon_sz_time(3) ! [enm] Dimension IDs
    integer dim_lon_time(2)   ! [enm] Dimension IDs
    integer fll_mode_old      ! Old fill mode
    integer lev_dim_id        ! [enm] Dimension ID for lev
    integer levp_dim_id       ! [enm] Dimension ID for levp
    integer lon_dim_id        ! [enm] Dimension ID for lon
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer srt_lon_lev_sz_time(4) ! Starting index array
    integer srt_lon_levp_sz_time(4) ! Starting index array
    integer srt_lon_lev_time(3) ! Starting index array
    integer srt_lon_levp_time(3) ! Starting index array
    integer srt_lon_sz_time(3) ! Starting index array
    integer srt_lon_time(2)   ! Starting index array
    integer sz_dim_id         ! [enm] Dimension ID for sz
    integer time_dim_id       ! [enm] Dimension ID for time
    ! Variable IDs
    integer cld_frc_cnv_id    ! [enm] Variable ID
    integer cld_frc_id        ! [enm] Variable ID
    integer cld_vlm_id        ! [enm] Variable ID
    integer q_H2O_cnd_id      ! [enm] Variable ID
    integer q_H2O_cnd2pcp_tnd_id ! [enm] Variable ID
    integer q_H2O_pcp2vpr_tnd_id ! [enm] Variable ID
    integer q_H2O_vpr2pcp_cnv_tnd_id ! [enm] Variable ID
    integer spc_xsx_ncl_scv_id ! [enm] Variable ID
    integer scv_cff_mss_avg_pcp_nrm_cnv_id ! [enm] Variable ID
    integer scv_cff_mss_avg_pcp_nrm_str_id ! [enm] Variable ID
    integer flx_mss_pcp_sfc_id ! [enm] Variable ID
    integer flx_mss_pcp_sfc_ttl_id ! [enm] Variable ID
    integer frc_trc_trn_cnv_ptn_id ! [enm] Variable ID
    integer pcp_flx_sfc_id    ! [enm] Variable ID
    integer q_dst_tnd_evp_id  ! [enm] Variable ID
    integer q_dst_tnd_evp_ttl_id ! [enm] Variable ID
    integer q_dst_tnd_ncl_id  ! [enm] Variable ID
    integer q_dst_tnd_pcp_id  ! [enm] Variable ID
    integer q_dst_tnd_pcp_ttl_id ! [enm] Variable ID
    integer q_dst_tnd_wet_id  ! [enm] Variable ID
    ! Variable data
    
    ! Begin netCDF output routines
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
    rcd=rcd+nf90_redef(nc_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'sz',sz_dim_id)
    
    rcd=nf90_wrp_inq_dimid(nc_id,'levp',levp_dim_id)
    
    ! Add global attributes
    ! Define dimension IDs
    ! Assemble ID and count vectors for each multidimensional combination of dimensions
    dim_lon_time=(/lon_dim_id,time_dim_id/)
    cnt_lon_time=(/plon,1/)
    srt_lon_time=(/1,nstep/)
    
    dim_lon_sz_time=(/lon_dim_id,sz_dim_id,time_dim_id/)
    cnt_lon_sz_time=(/plon,dst_nbr,1/)
    srt_lon_sz_time=(/1,1,nstep/)
    
    dim_lon_lev_time=(/lon_dim_id,lev_dim_id,time_dim_id/)
    cnt_lon_lev_time=(/plon,plev,1/)
    srt_lon_lev_time=(/1,1,nstep/)
    
    dim_lon_lev_sz_time=(/lon_dim_id,lev_dim_id,sz_dim_id,time_dim_id/)
    cnt_lon_lev_sz_time=(/plon,plev,dst_nbr,1/)
    srt_lon_lev_sz_time=(/1,1,1,nstep/)
    
    dim_lon_levp_time=(/lon_dim_id,levp_dim_id,time_dim_id/)
    cnt_lon_levp_time=(/plon,plevp,1/)
    srt_lon_levp_time=(/1,1,nstep/)
    
    dim_lon_levp_sz_time=(/lon_dim_id,levp_dim_id,sz_dim_id,time_dim_id/)
    cnt_lon_levp_sz_time=(/plon,plevp,dst_nbr,1/)
    srt_lon_levp_sz_time=(/1,1,1,nstep/)
    
    if (nstep == 1) then
       ! Variable definitions
       rcd=rcd+nf90_def_var(nc_id,'cld_frc',nf90_float,dim_lon_levp_time,cld_frc_id)
       rcd=rcd+nf90_def_var(nc_id,'cld_frc_cnv',nf90_float,dim_lon_lev_time,cld_frc_cnv_id)
       rcd=rcd+nf90_def_var(nc_id,'cld_vlm',nf90_float,dim_lon_lev_time,cld_vlm_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp_sfc',nf90_float,dim_lon_sz_time,flx_mss_pcp_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp_sfc_ttl',nf90_float,dim_lon_time,flx_mss_pcp_sfc_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'frc_trc_trn_cnv_ptn',nf90_float,dim_lon_lev_sz_time,frc_trc_trn_cnv_ptn_id)
       rcd=rcd+nf90_def_var(nc_id,'pcp_flx_sfc',nf90_float,dim_lon_time,pcp_flx_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'q_H2O_cnd',nf90_float,dim_lon_lev_time,q_H2O_cnd_id)
       rcd=rcd+nf90_def_var(nc_id,'q_H2O_cnd2pcp_tnd',nf90_float,dim_lon_lev_time,q_H2O_cnd2pcp_tnd_id)
       rcd=rcd+nf90_def_var(nc_id,'q_H2O_pcp2vpr_tnd',nf90_float,dim_lon_lev_time,q_H2O_pcp2vpr_tnd_id)
       rcd=rcd+nf90_def_var(nc_id,'q_H2O_vpr2pcp_cnv_tnd',nf90_float,dim_lon_lev_time,q_H2O_vpr2pcp_cnv_tnd_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_evp',nf90_float,dim_lon_lev_sz_time,q_dst_tnd_evp_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_evp_ttl',nf90_float,dim_lon_lev_time,q_dst_tnd_evp_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_ncl',nf90_float,dim_lon_lev_sz_time,q_dst_tnd_ncl_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_pcp',nf90_float,dim_lon_lev_sz_time,q_dst_tnd_pcp_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_pcp_ttl',nf90_float,dim_lon_lev_time,q_dst_tnd_pcp_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_wet',nf90_float,dim_lon_lev_sz_time,q_dst_tnd_wet_id)
       rcd=rcd+nf90_def_var(nc_id,'scv_cff_mss_avg_pcp_nrm_cnv',nf90_float,sz_dim_id,scv_cff_mss_avg_pcp_nrm_cnv_id)
       rcd=rcd+nf90_def_var(nc_id,'scv_cff_mss_avg_pcp_nrm_str',nf90_float,sz_dim_id,scv_cff_mss_avg_pcp_nrm_str_id)
       rcd=rcd+nf90_def_var(nc_id,'spc_xsx_ncl_scv',nf90_float,sz_dim_id,spc_xsx_ncl_scv_id)
       ! Add english text descriptions
       rcd=rcd+nf90_put_att(nc_id,cld_frc_cnv_id,'long_name','Convective cloud fraction')
       rcd=rcd+nf90_put_att(nc_id,cld_frc_id,'long_name','Cloud fraction')
       rcd=rcd+nf90_put_att(nc_id,cld_vlm_id,'long_name','Cloud volume')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_id,'long_name','Dust reaching surface in precipitation')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_ttl_id,'long_name','Total dust reaching surface in precipitation')
       rcd=rcd+nf90_put_att(nc_id,frc_trc_trn_cnv_ptn_id,'long_name','Interstitial tracer fraction')
       rcd=rcd+nf90_put_att(nc_id,pcp_flx_sfc_id,'long_name','Total precipitation reaching surface')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd2pcp_tnd_id,'long_name','Condensed H2O to precipitation tendency')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd_id,'long_name','Condensed H2O mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_pcp2vpr_tnd_id,'long_name','H2O precipitation to vapor tendency')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,'long_name','H2O vapor to convective precipitation tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_id,'long_name','Evaporation tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_ttl_id,'long_name','Total evaporation tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_ncl_id,'long_name','Nucleation scavenging tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_id,'long_name','Scavenging tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_ttl_id,'long_name','Total scavenging tendency')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_wet_id,'long_name','Wet deposition (evaporation minus scavenging) tendency')
       rcd=rcd+nf90_put_att(nc_id,spc_xsx_ncl_scv_id,'long_name','Specific cross section for nucleation scavenging')
       rcd=rcd+nf90_put_att(nc_id,scv_cff_mss_avg_pcp_nrm_cnv_id,'long_name', &
            'Mass mean scavenging coefficient, precipitation normalized, convective')
       rcd=rcd+nf90_put_att(nc_id,scv_cff_mss_avg_pcp_nrm_str_id,'long_name', &
            'Mass mean scavenging coefficient, precipitation normalized, stratiform')
       ! Add units
       rcd=rcd+nf90_put_att(nc_id,cld_frc_cnv_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,cld_frc_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,cld_vlm_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_pcp_sfc_ttl_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,frc_trc_trn_cnv_ptn_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,pcp_flx_sfc_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd2pcp_tnd_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_cnd_id,'units','kilogram kilogram-1')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_pcp2vpr_tnd_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_evp_ttl_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_ncl_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_pcp_ttl_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_wet_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,scv_cff_mss_avg_pcp_nrm_cnv_id,'units','meter2 kilogram-1')
       rcd=rcd+nf90_put_att(nc_id,scv_cff_mss_avg_pcp_nrm_str_id,'units','meter2 kilogram-1')
       rcd=rcd+nf90_put_att(nc_id,spc_xsx_ncl_scv_id,'units','meter2 kilogram-1')
    else                      ! endif nstep == 1
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,'cld_frc',cld_frc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'cld_frc_cnv',cld_frc_cnv_id)
       rcd=nf90_wrp_inq_varid(nc_id,'cld_vlm',cld_vlm_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_pcp_sfc',flx_mss_pcp_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_pcp_sfc_ttl',flx_mss_pcp_sfc_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'frc_trc_trn_cnv_ptn',frc_trc_trn_cnv_ptn_id)
       rcd=nf90_wrp_inq_varid(nc_id,'pcp_flx_sfc',pcp_flx_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_H2O_cnd',q_H2O_cnd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_H2O_cnd2pcp_tnd',q_H2O_cnd2pcp_tnd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_H2O_pcp2vpr_tnd',q_H2O_pcp2vpr_tnd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_H2O_vpr2pcp_cnv_tnd',q_H2O_vpr2pcp_cnv_tnd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_evp',q_dst_tnd_evp_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_evp_ttl',q_dst_tnd_evp_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_ncl',q_dst_tnd_ncl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_pcp',q_dst_tnd_pcp_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_pcp_ttl',q_dst_tnd_pcp_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_wet',q_dst_tnd_wet_id)
    endif                     ! endif nstep /= 1
    ! End define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    if (nstep == 1) then
       rcd=rcd+nf90_put_var(nc_id,spc_xsx_ncl_scv_id,spc_xsx_ncl_scv)
       rcd=rcd+nf90_put_var(nc_id,scv_cff_mss_avg_pcp_nrm_cnv_id,scv_cff_mss_avg_pcp_nrm_cnv)
       rcd=rcd+nf90_put_var(nc_id,scv_cff_mss_avg_pcp_nrm_str_id,scv_cff_mss_avg_pcp_nrm_str)
    endif                     ! endif nstep == 1
    rcd=rcd+nf90_put_var(nc_id,cld_frc_cnv_id,cld_frc_cnv,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,cld_frc_id,cld_frc,start=srt_lon_levp_time,count=cnt_lon_levp_time)
    rcd=rcd+nf90_put_var(nc_id,cld_vlm_id,cld_vlm,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_sfc_id,flx_mss_pcp_sfc,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_sfc_ttl_id,flx_mss_pcp_sfc_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,frc_trc_trn_cnv_ptn_id,frc_trc_trn_cnv_ptn,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,pcp_flx_sfc_id,pcp_flx_sfc,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,q_H2O_cnd2pcp_tnd_id,q_H2O_cnd2pcp_tnd,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_H2O_cnd_id,q_H2O_cnd,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_H2O_pcp2vpr_tnd_id,q_H2O_pcp2vpr_tnd,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_H2O_vpr2pcp_cnv_tnd_id,q_H2O_vpr2pcp_cnv_tnd,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_evp_id,q_dst_tnd_evp,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_evp_ttl_id,q_dst_tnd_evp_ttl,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_ncl_id,q_dst_tnd_ncl,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_pcp_id,q_dst_tnd_pcp,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_pcp_ttl_id,q_dst_tnd_pcp_ttl,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_wet_id,q_dst_tnd_wet,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
    if (nstep == 1) then
       write (6,'(a,a45,a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': Initialized wet deposition data archive in ',fl_out(1:ftn_strlen(fl_out))
    endif                     ! endif
  end subroutine dpswet2nc                       ! end dpswet2nc()
  
end module dstdpswet

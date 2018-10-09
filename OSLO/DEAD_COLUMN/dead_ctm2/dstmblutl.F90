! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstmblutl.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: dstmblutl.F contains routines to compute mineral dust source tendencies
! Usage:
! use dstmblutl ! [mdl] Mobilization utilities

module dstmblutl ! [mdl] Mobilization utilities
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  private::tm_2_idx_wgt ! [fnc] Day of year --> bounding indices, weights
  private::dst_lsm_ini ! [fnc] Erodible fraction of gridcell
  private::vai_lsm_get ! [fnc] Update time-varying surface fields

contains
  
  subroutine wnd_frc_thr_slt_get( &
       dmt_aer,             & ! I [m] Particle diameter (Currently not used)
       dns_aer,             & ! I [kg m-3] Particle density (Currently not used)
       dns_mdp,             & ! I [kg m-3] Midlayer density
       wnd_frc_thr_slt)     ! O [m s-1] Threshold friction velocity for saltation
    ! Purpose: Compute dry threshold friction velocity for saltation
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    ! Parameters
    real(r8),parameter::dmt_slt_opt=75.0e-6 ! [m] Optimal diameter for saltation, IvW82 p. 117 Fgr. 8, Pye87 p. 31, MBA97 p. 4388, SRL96 (2)
    real(r8),parameter::dns_slt=2650.0 ! [kg m-3] Density of optimal saltation particles, MBA97 p. 4388 
    ! Input
    real(r8),intent(in)::dns_aer(dst_nbr)     ! [kg m-3] Particle density
    real(r8),intent(in)::dmt_aer(dst_nbr)     ! [m] Particle diameter
    real(r8),intent(in)::dns_mdp(plond)       ! [kg m-3] Midlayer density
    ! Output
    real(r8),intent(out)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    ! Local
    integer lon_idx           ! [idx] Counting index
    real(r8) ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) dns_fct              ! Density ratio factor for saltation calculation
    real(r8) icf_fct              ! Interparticle cohesive forces factor for saltation calculation
    real(r8) tmp1                 ! Factor in saltation computation
    ! Main Code
    ! Initialize some variables
    ! MaB95 pzn. for Re*t(D_opt) circumvents iterative solution
    ryn_nbr_frc_thr_prx_opt=0.38_r8+1331.0_r8*(100.0_r8*dmt_slt_opt)**1.56_r8 ! [frc] "B" MaB95 p. 16417 (5)
    ! Given Re*t(D_opt), compute time independent factors contributing to u*t
    icf_fct=1.0_r8+6.0e-07_r8/(dns_slt*grv_sfc*(dmt_slt_opt**2.5_r8)) ! [frc] IvW82 p. 115 (6) MaB95 p. 16417 (4) Interparticle cohesive forces
    dns_fct=dns_slt*grv_sfc*dmt_slt_opt ! IvW82 p. 115 (6) MaB95 p. 16417 (4)
    if (ryn_nbr_frc_thr_prx_opt < 0.03_r8) then
       stop 'dst: wnd_frc_thr_slt_get() reports ryn_nbr_frc_thr_prx_opt < 0.03'
    else if (ryn_nbr_frc_thr_prx_opt < 10.0_r8) then
       ryn_nbr_frc_thr_opt_fnc=-1.0_r8+1.928_r8*(ryn_nbr_frc_thr_prx_opt**0.0922_r8) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (6)
       ryn_nbr_frc_thr_opt_fnc=0.1291_r8*0.1291_r8/ryn_nbr_frc_thr_opt_fnc ! [frc] 
    else 
       ryn_nbr_frc_thr_opt_fnc=1.0_r8-0.0858_r8*exp(-0.0617_r8*(ryn_nbr_frc_thr_prx_opt-10.0_r8)) ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (7)
       ryn_nbr_frc_thr_opt_fnc=0.120_r8*0.120_r8*ryn_nbr_frc_thr_opt_fnc*ryn_nbr_frc_thr_opt_fnc! [frc]
    endif                     ! endif
    ! This method minimizes the number of square root computations performed
    tmp1=sqrt(icf_fct*dns_fct*ryn_nbr_frc_thr_opt_fnc)
    do lon_idx=1,plon
       wnd_frc_thr_slt(lon_idx)=tmp1/sqrt(dns_mdp(lon_idx)) ! [m s-1] Threshold friction velocity for saltation dry ground 
    end do                    ! end loop over lon
    return
  end subroutine wnd_frc_thr_slt_get                       ! end wnd_frc_thr_slt_get()
  
  subroutine wnd_rfr_thr_slt_get( &
       wnd_frc,             & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt,     & ! I [m s-1] Threshold friction velocity for saltation
       wnd_mdp,             & ! I [m s-1] Surface layer mean wind speed
       wnd_rfr,             & ! I [m s-1] Wind speed at reference height
       wnd_rfr_thr_slt)     ! O [m s-1] Threshold 10 m wind speed for saltation
    ! Purpose: Compute threshold 10 m wind speed
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    real(r8),intent(in)::wnd_frc(plond)       ! [m s-1] Friction velocity
    real(r8),intent(in)::wnd_frc_thr_slt(plond) ! [m s-1] Threshold friction velocity for saltation
    real(r8),intent(in)::wnd_mdp(plond)       ! [m s-1] Surface layer mean wind speed
    real(r8),intent(in)::wnd_rfr(plond)       ! [m s-1] Wind speed at reference height
    ! Output
    real(r8),intent(out)::wnd_rfr_thr_slt(plond) ! [m s-1] Threshold 10 m wind speed for saltation
    ! Local
    integer lon_idx           ! [idx] Counting index
    ! Main Code
    ! Compute threshold horizontal wind speed at reference height
    do lon_idx=1,plon
       ! A more complicated procedure would recompute mno_lng for wnd_frc_thr,
       ! and then integrate vertically from rgh_mmn+hgt_zpd to hgt_rfr
       ! wnd_crc_fct is (1/k)*[ln(z-D)/z0 - psi(zeta2) + psi(zeta1)]
       wnd_rfr_thr_slt(lon_idx)=wnd_frc_thr_slt(lon_idx)*wnd_rfr(lon_idx)/wnd_frc(lon_idx) ! [m s-1]
    end do                    ! end loop over lon
    return
  end subroutine wnd_rfr_thr_slt_get                       ! end wnd_rfr_thr_slt_get()
  
  subroutine vwc2gwc( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       gwc_sfc,             & ! O [kg kg-1] Gravimetric water content
       vwc_sat,             & ! I [m3 m-3] Saturated volumetric water content (sand-dependent)
       vwc_sfc)             ! I [m3 m-3] Volumetric water content
    ! Purpose: Convert volumetric water content to gravimetric water content
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstblm,only:dns_H2O_lqd_std ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    implicit none
    ! Parameters
    real(r8),parameter::dns_prt_sfc=2650.0 ! [kg m-3] Dry density of soil particles (excluding pores)
    ! Input
    logical,intent(in)::flg_mbl(plond) ! [flg] Mobilization candidate flag
    real(r8),intent(in)::vwc_sat(plond) ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    real(r8),intent(in)::vwc_sfc(plond) ! [m3 m-3] Volumetric water content
    ! Output
    real(r8),intent(out)::gwc_sfc(plond) ! [kg kg-1] Gravimetric water content
    ! Local
    integer lon_idx           ! [idx] Counting index
    real(r8) dns_blk_dry(plond) ! [kg m-3] Bulk density of dry surface soil
    ! Main Code
    ! Initialize output
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
       ! Assume volume of air pores when dry equals saturated VWC
       ! This implies air pores are completely filled by water in saturated soil
          dns_blk_dry(lon_idx)=dns_prt_sfc*(1.0-vwc_sat(lon_idx)) ! [kg m-3] Bulk density of dry surface soil
          gwc_sfc(lon_idx)=vwc_sfc(lon_idx)*dns_H2O_lqd_std/dns_blk_dry(lon_idx) ! O [kg kg-1] Gravimetric water content
       endif ! endif flg_mbl
    enddo ! end loop over lon
    return
  end subroutine vwc2gwc
  
  subroutine frc_thr_ncr_wtr_get( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       frc_thr_ncr_wtr,     & ! O [frc] Factor by which moisture increases threshold friction velocity
       mss_frc_cly,         & ! I [frc] Mass fraction of clay
       gwc_sfc)             ! I [kg kg-1] Gravimetric water content
    ! Purpose: Compute factor by which soil moisture increases threshold friction velocity 
    ! This parameterization is based on FMB99
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::mss_frc_cly(plond) ! [frc] Mass fraction of clay
    real(r8),intent(in)::gwc_sfc(plond) ! [kg kg-1] Gravimetric water content
    ! Output
    real(r8),intent(out)::frc_thr_ncr_wtr(plond) ! [frc] Factor by which moisture increases threshold friction velocity
    ! Local
    integer lon_idx           ! [idx] Counting index
    real(r8) gwc_thr(plond)       ! [kg kg-1] Threshold gravimetric water content
    ! Main Code
    ! Initialize output
    frc_thr_ncr_wtr(:)=1.0_r8 ! [frc] Factor by which moisture increases threshold friction velocity
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
          ! Adjust threshold velocity for inhibition by moisture
          ! frc_thr_ncr_wtr(lon_idx)=exp(22.7_r8*vwc_sfc(lon_idx)) ! [frc] SRL96
          ! Compute threshold soil moisture based on clay content
          ! gwc_thr=mss_frc_cly*(0.17_r8+0.14_r8*mss_frc_cly) ! [m3 m-3] FMB99 p. 155 (14)
          ! fxm: 19991105 remove factor of mss_frc_cly from gwc_thr to improve large scale behavior
          gwc_thr(lon_idx)=0.17_r8+0.14_r8*mss_frc_cly(lon_idx) ! [m3 m-3] 
          if (gwc_sfc(lon_idx) > gwc_thr(lon_idx)) &
               frc_thr_ncr_wtr(lon_idx)=sqrt(1.0_r8+1.21_r8*(100.0_r8*(gwc_sfc(lon_idx)-gwc_thr(lon_idx)))**0.68_r8) ! [frc] FMB99 p. 155 (15)
       endif ! endif flg_mbl
    enddo ! end loop over lon
    ! Uncomment following line to remove all dependence on gwc_sfc
    ! frc_thr_ncr_wtr(lon_idx)=1.0_r8 ! [frc]
    return
  end subroutine frc_thr_ncr_wtr_get
  
  subroutine frc_thr_ncr_drg_get( &
       frc_thr_ncr_drg,     & ! O [frc] Factor by which surface roughness increases threshold friction velocity
       rgh_mmn_mbl,         & ! I [m] Roughness length momentum for erodible surfaces
       rgh_mmn_smt)         ! I [m] Smooth roughness length
    ! Purpose: Compute factor by which surface roughness increases threshold friction velocity
    ! This parameterization is based on MaB95 and GMB98
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    real(r8),intent(in)::rgh_mmn_mbl ! [m] Roughness length momentum for erodible surfaces
    real(r8),intent(in)::rgh_mmn_smt ! [m] Smooth roughness length
    ! Output
    real(r8),intent(out)::frc_thr_ncr_drg(plond) ! [frc] Factor by which roughness increases threshold friction velocity
    ! Local
    integer lon_idx           ! [idx] Counting index
    real(r8) wnd_frc_fsh_frc      ! [frc] Efficient fraction of wind friction
    real(r8) wnd_frc_fsh_frc_rcp  ! [frc] Reciprocal of wnd_frc_fsh_frc
    ! Main Code
    ! Adjust threshold velocity for inhibition by roughness elements
    wnd_frc_fsh_frc=          & ! [frc] MaB95 p. 16420, GMB98 p. 6207
         +1.0_r8-log(rgh_mmn_mbl/rgh_mmn_smt)/log(0.35_r8*((0.1/rgh_mmn_smt)**0.8_r8))
    if (wnd_frc_fsh_frc <= 0.0_r8.or.wnd_frc_fsh_frc > 1.0_r8) stop  &
         'dst: ERROR frc_thr_ncr_drg_get() reports wnd_frc_fsh_frc out of range'
    wnd_frc_fsh_frc_rcp=1.0_r8/wnd_frc_fsh_frc ! [frc] 
    frc_thr_ncr_drg=wnd_frc_fsh_frc_rcp ! [frc]
    ! fxm: 19991012 Set frc_thr_ncr_drg=1.0, equivalent to assuming mobilization takes place at smooth roughness length
    frc_thr_ncr_drg=1.0_r8
    return
  end subroutine frc_thr_ncr_drg_get                       ! end frc_thr_ncr_drg_get()
  
  subroutine wnd_frc_slt_get( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       wnd_frc,             & ! I [m s-1] Surface friction velocity
       wnd_frc_slt,         & ! O [m s-1] Saltating friction velocity
       wnd_rfr,             & ! I [m s-1] Wind speed at reference height
       wnd_rfr_thr_slt)     ! I [m s-1] Threshold 10 m wind speed for saltation
    ! Purpose: Compute the saltating friction velocity
    ! Saltation increases friction speed by roughening surface, AKA "Owen's effect"
    ! This acts as a positive feedback to the friction speed
    ! GMB98 parameterized this feedback in terms of 10 m windspeeds
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::wnd_frc(plond)       ! I [m s-1] Surface friction velocity
    real(r8),intent(in)::wnd_rfr(plond)       ! I [m s-1] Wind speed at reference height
    real(r8),intent(in)::wnd_rfr_thr_slt(plond) ! I [m s-1] Threshold 10 m wind speed for saltation
    ! Output
    real(r8),intent(out)::wnd_frc_slt(plond)   ! O [m s-1] Saltating friction velocity
    ! Local
    integer lon_idx           ! [idx] Counting index
    real(r8) wnd_rfr_dlt          ! [m s-1] Reference windspeed excess over threshold
    real(r8) wnd_frc_slt_dlt      ! [m s-1] Friction velocity increase from saltation
    ! Main Code
    ! Compute saltating friction velocity, accounting for "Owen's effect"
    wnd_frc_slt=wnd_frc       ! [m s-1] Saltating friction velocity
    do lon_idx=1,plon
       if (flg_mbl(lon_idx).and.wnd_rfr(lon_idx) >= wnd_rfr_thr_slt(lon_idx)) then
          ! Saltation roughens the boundary layer, AKA "Owen's effect"
          ! GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
          ! GMB98 p. 6209 (12) has u* in cm s-1 and U, Ut in m s-1, personal communication, D. Gillette, 19990529
          ! With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
          ! Increase in friction velocity due to saltation varies as square of 
          ! difference between reference wind speed and reference threshold speed 
          wnd_rfr_dlt=wnd_rfr(lon_idx)-wnd_rfr_thr_slt(lon_idx)
          wnd_frc_slt_dlt=0.003_r8*wnd_rfr_dlt*wnd_rfr_dlt ! [m s-1] Friction velocity increase from saltation GMB98 p. 6209
          wnd_frc_slt(lon_idx)=wnd_frc(lon_idx)+wnd_frc_slt_dlt ! [m s-1] Saltating friction velocity
       endif                  ! endif wnd_frc_mbl > wnd_frc_thr_slt
    end do                    ! end loop over lon
    return
  end subroutine wnd_frc_slt_get                       ! end wnd_frc_slt_get()
  
  subroutine flx_mss_CaCO3_msk( &
       dmt_vwr,             & ! I [m] Mass weighted diameter resolved 
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_mss_vrt_dst_CaCO3, & ! I/O [kg m-2 s-1] Vertical mass flux of dust on input, CaCO3 on output
       mss_frc_CaCO3,       & ! I [frc] Mass fraction of CaCO3
       mss_frc_cly,         & ! I [frc] Mass fraction of clay
       mss_frc_snd)         ! I [frc] Mass fraction of sand
    ! Purpose: Mask dust mass flux by CaCO3 mass fraction at source
    ! Theory: Uses soil CaCO3 mass fraction from Global Soil Data Task, 1999 (Sch99)
    ! Uses size dependent apportionment of CaCO3 from Claquin et al, 1999 (CSB99)
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    real(r8),parameter::dmt_cly_max=2.0e-6 ! [m] Maximum diameter of Clay soil texture CSB99 p. 22250
    real(r8),parameter::dmt_slt_max=50.0e-6 ! [m] Maximum diameter of Silt soil texture CSB99 p. 22250
    real(r8),parameter::dns_CaCO3=2950.0 ! [kg m-3] Density of CaCO3 http://www.ssc.on.ca/mandm/calcit.htm
    ! Input
    logical,intent(in)::flg_mbl(plond) ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::dmt_vwr(dst_nbr) ! I [m] Mass weighted diameter resolved
    real(r8),intent(in)::mss_frc_CaCO3(plond) ! I [frc] Mass fraction of CaCO3
    real(r8),intent(in)::mss_frc_cly(plond) ! I [frc] Mass fraction of clay
    real(r8),intent(in)::mss_frc_snd(plond) ! I [frc] Mass fraction of sand
    ! Output
    ! Input/Output
    ! Fluxes refer to dust on input and to CaCO3 on output
    real(r8),intent(inout)::flx_mss_vrt_dst_CaCO3(plond,dst_nbr) ! I/O [kg m-2 s-1] Vertical mass flux of dust on input, CaCO3 on output
    integer rcd               ! [rcd] Return success code
    ! Local
    integer m                 ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index for lon
    real(r8) mss_frc_slt(plond)   ! [frc] Mass fraction of silt
    real(r8) mss_frc_CaCO3_sz_crr ! [frc] Fraction of soil CaCO3 in size bin
    real(r8) mss_frc_CaCO3_cly    ! [frc] Fraction of CaCO3 in clay
    real(r8) mss_frc_CaCO3_slt    ! [frc] Fraction of CaCO3 in silt
    real(r8) mss_frc_CaCO3_snd    ! [frc] Fraction of CaCO3 in sand
    ! Main Code
    ! Initialize output
    do m=1,dst_nbr
       do lon_idx=1,plon
          ! Simple technique is to mask dust mass by tracer mass fraction
          ! Model transports (hence conserves) CaCO3 rather than total dust itself
          ! Method assumes source, transport, and removal processes are linear with tracer mass
          if (flg_mbl(lon_idx)) then
             ! 20000320: Currently this is only process  in dust model requiring mss_frc_slt
             mss_frc_slt(lon_idx)= & ! [frc] Mass fraction of silt
                  max(0.0_r8,1.0_r8-mss_frc_cly(lon_idx)-mss_frc_snd(lon_idx))
             ! CSB99 showed that CaCO3 is not uniformly distributed across sizes
             ! There is more CaCO3 per unit mass of silt than per unit mass of clay
             mss_frc_CaCO3_cly=max(0.0_r8,-0.045_r8+0.5_r8*min(0.5_r8,mss_frc_cly(lon_idx))) ! [frc] Fraction of CaCO3 in clay CSB99 p. 22249 Figure 1b
             mss_frc_CaCO3_slt=max(0.0_r8,-0.175_r8+1.4_r8*min(0.5_r8,mss_frc_slt(lon_idx))) ! [frc] Fraction of CaCO3 in silt CSB99 p. 22249 Figure 1a
             mss_frc_CaCO3_snd=1.0_r8-mss_frc_CaCO3_cly-mss_frc_CaCO3_snd ! [frc] Fraction of CaCO3 in sand CSB99 p. 22249 Figure 1a
             ! Set CaCO3 fraction of total CaCO3 for each transport bin 
             if (dmt_vwr(m) < dmt_cly_max) then
                ! Transport bin carries Clay
                mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_cly ! [frc] Fraction of soil CaCO3 in size bin
             else if (dmt_vwr(m) < dmt_slt_max) then ! endif size bin carries clay
                ! Transport bin carries Silt
                mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_slt ! [frc] Fraction of soil CaCO3 in size bin
             else             ! endif size bin carries silt
                ! Transport bin carries Sand
                mss_frc_CaCO3_sz_crr=mss_frc_CaCO3_snd ! [frc] Fraction of soil CaCO3 in size bin
             endif            ! endif size bin carries sand
#ifdef DST_DBG
             ! Sanity check
             if (mss_frc_CaCO3_sz_crr < 0.0_r8 .or. mss_frc_CaCO3_sz_crr > 1.0_r8) stop &
                  'dst: flx_mss_CaCO3_msk() mss_frc_CaCO3_sz_crr < 0.0.or.mss_frc_CaCO3_sz_crr > 1.0'
             if (mss_frc_CaCO3(lon_idx) < 0.0_r8 .or. mss_frc_CaCO3(lon_idx) > 1.0_r8) stop &
                  'dst: flx_mss_CaCO3_msk() mss_frc_CaCO3 < 0.0.or.mss_frc_CaCO3 > 1.0'
#endif /* not DST_DBG */
             
             ! Convert dust flux to CaCO3 flux
             flx_mss_vrt_dst_CaCO3(lon_idx,m)=flx_mss_vrt_dst_CaCO3(lon_idx,m) & ! [kg m-2 s-1]
                  *mss_frc_CaCO3(lon_idx) & ! [frc] Mass fraction of CaCO3 (at this location)
                  ! 20020925 fxm: Remove size dependence of CaCO3
                  *1.0_r8
                  !                  *mss_frc_CaCO3_sz_crr ! [frc] Fraction of soil CaCO3 in size bin
          endif               ! endif 
       end do                  ! end loop over lon
    end do                     ! end loop over sz
    
    return
  end subroutine flx_mss_CaCO3_msk                       ! end flx_mss_CaCO3_msk()
  
  subroutine flx_mss_hrz_slt_ttl_Whi79_get( &
       dns_mdp,             & ! I [kg m-3] Midlayer density
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
       wnd_frc,             & ! I [m s-1] Surface friction velocity
       wnd_frc_thr_slt)     ! I [m s-1] Threshold friction speed for saltation
    ! Purpose: Compute vertically integrated streamwise mass flux of particles
    ! Theory: Uses method proposed by White (1979)
    ! fxm: use surface air density not midlayer density
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    ! Parameters
    real(r8),parameter::cst_slt=2.61 ! [frc] Saltation constant Whi79 p. 4648, MaB97 p. 16422 
    ! Input
    logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::dns_mdp(plond)       ! I [kg m-3] Midlayer density
    real(r8),intent(in)::wnd_frc(plond)       ! I [m s-1] Surface friction velocity
    real(r8),intent(in)::wnd_frc_thr_slt(plond) ! I [m s-1] Threshold friction speed for saltation
    ! Output
    real(r8),intent(out)::flx_mss_hrz_slt_ttl(plond) ! O [kg m-1 s-1] Vertically integrated streamwise mass flux
    integer rcd               ! [rcd] Return success code
    ! Local
    real(r8) wnd_frc_rat          ! [frc] Ratio of wind friction threshold to wind friction
    integer lon_idx           ! [idx] Counting index for lon
    ! Main Code
    ! Initialize output
    flx_mss_hrz_slt_ttl=0.0_r8   ! [kg m-1 s-1]
    
    do lon_idx=1,plon
       if (flg_mbl(lon_idx).and.wnd_frc(lon_idx) > wnd_frc_thr_slt(lon_idx)) then
          wnd_frc_rat=wnd_frc_thr_slt(lon_idx)/wnd_frc(lon_idx) ! [frc]
          flx_mss_hrz_slt_ttl(lon_idx)= & ! [kg m-1 s-1] 
               cst_slt*dns_mdp(lon_idx)*(wnd_frc(lon_idx)**3.0_r8)* &
               (1.0_r8-wnd_frc_rat)*(1.0_r8+wnd_frc_rat)*(1.0_r8+wnd_frc_rat)/grv_sfc ! Whi79 p. 4648 (19), MaB97 p. 16422 (28)
       endif                  ! endif 
    end do                     ! end loop over lon
    
    return
  end subroutine flx_mss_hrz_slt_ttl_Whi79_get                       ! end flx_mss_hrz_slt_ttl_Whi79_get()
  
  subroutine flx_mss_vrt_dst_ttl_MaB95_get( &
       dst_slt_flx_rat_ttl, & ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_mss_hrz_slt_ttl, & ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
       flx_mss_vrt_dst_ttl, & ! O [kg m-2 s-1] Total vertical mass flux of dust
       mss_frc_cly)         ! I [frc] Mass fraction clay 
    ! Purpose: Diagnose total vertical mass flux of dust from vertically integrated streamwise mass flux
    ! Theory: Uses clay-based method proposed by Marticorena & Bergametti (1995)
    ! Their parameterization is based only on data for mss_frc_cly < 0.20
    ! For clayier soils, dst_slt_flx_rat_ttl may behave dramatically differently
    ! Whether this behavior changes when mss_frc_cly > 0.20 is unknown
    ! Anecdotal evidence suggests vertical flux decreases for mss_frc_cly > 0.20
    ! Thus we use min[mss_frc_cly,0.20] in MaB95 parameterization
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    ! Input
    logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::flx_mss_hrz_slt_ttl(plond) ! I [kg m-1 s-1] Vertically integrated streamwise mass flux
    real(r8),intent(in)::mss_frc_cly(plond)   ! I [frc] Mass fraction clay 
    ! Output
    real(r8),intent(out)::dst_slt_flx_rat_ttl(plond) ! O [m-1] Ratio of vertical dust flux to streamwise mass flux
    real(r8),intent(out)::flx_mss_vrt_dst_ttl(plond) ! O [kg m-2 s-1] Total vertical mass flux of dust
    integer rcd               ! [rcd] Return success code
    ! Local
    integer lon_idx           ! [idx] Counting index for lon
    real(r8) mss_frc_cly_vld      ! [frc] Mass fraction clay limited to 0.20
    real(r8) ln10                 ! [frc] Natural log of 10
    
    ! Initialize some variables
    ln10=log(10.0_r8)            ! [frc] Natural log of 10
    do lon_idx=1,plon
       if (flg_mbl(lon_idx)) then
          ! 19990603: fxm: Dust production is EXTREMELY sensitive to this parameter, which changes flux by 3 orders of magnitude in 0.0 < mss_frc_cly < 0.20
          mss_frc_cly_vld=min(mss_frc_cly(lon_idx),0.2_r8) ! [frc]
          dst_slt_flx_rat_ttl(lon_idx)= & ! [m-1]
               100.0_r8*exp(ln10*(13.4_r8*mss_frc_cly_vld-6.0_r8)) ! MaB95 p. 16423 (47)
          flx_mss_vrt_dst_ttl(lon_idx)=flx_mss_hrz_slt_ttl(lon_idx)*dst_slt_flx_rat_ttl(lon_idx) ! [kg m-1 s-1] 
       endif                  ! endif flg_mbl
    end do                     ! end loop over lon
    
    return
  end subroutine flx_mss_vrt_dst_ttl_MaB95_get ! end flx_mss_vrt_dst_ttl_MaB95_get()

  !****************************************************************
  !April 2002: 3 new subroutines for AlG01 parameterization
  !***************************************************************

  subroutine mss_frc_src_AlG01_get(&
       wnd_frc,  &        ! I [m s-1] wind friction velocity
       mss_frc_src_lut, & ! I [frc] mass fraction source look up table
       mss_frc_src, &     ! O [frc] mass fraction of source (looked up)
       sfc_frc_bln )    ! I [frc] weighting of soil types 
    ! Purpose find mass fraction of each source mode from the three
    ! source modes of Alfaro/Gomes JGR 2001.
    ! Code by Alf Grini, UCI/UiO, 2002
    use dstgrd
    use pmgrid,only:plon,plond
    implicit none
    
    ! OUTPUT
    real(r8), intent(out):: mss_frc_src(plond,dst_src_nbr)  ! O [frc] mass fraction of source modes
    ! INPUT
    real(r8), intent(in) :: mss_frc_src_lut(wnd_frc_nbr,bln_nbr,dst_src_nbr) ! I Look up table for mss_frc_src
    real(r8), intent(in) :: wnd_frc(plond)            ! I [m/s] wind friction speed
    real(r8), intent(in) :: sfc_frc_bln(plond,bln_nbr) ! [frc] weighting of soil types
    !LOCAL
    integer         :: lon_idx     ! [idx] counting variable for longitude
    integer         :: src_idx     ! [idx] counting variable for source mode
    integer         :: soil_idx    ! [idx] counting variable for soil type
    integer         :: wnd_frc_idx ! [idx] u* to look up values for [1-100 (cm/s)]
    
    ! INITIALIZE
    mss_frc_src(:,:)=0.d0
    
    do src_idx=1,dst_src_nbr
       do soil_idx=1,bln_nbr
          do lon_idx=1,plon
             
             !Getting the right wind friction speed to look up 1-100               
             wnd_frc_idx=min(100,int(100.d0*wnd_frc(lon_idx)+0.5d0))
             
             !Summing up mss_frc_src getting sfc_frc_bln contribution from each soiltype
             mss_frc_src(lon_idx,src_idx) = & 
                  sfc_frc_bln(lon_idx,soil_idx) &
                  *mss_frc_src_lut(wnd_frc_idx,soil_idx,src_idx) &   !Contains the mass fraction of each soil 
                  + mss_frc_src(lon_idx,src_idx)                     !mss_frc_src is three first elements in 
             !look up table
          enddo  !Loop on lon
       enddo     !Loop on soil type
    enddo        !Loop on source mode
    
    return
  end subroutine mss_frc_src_AlG01_get
  
  subroutine ovr_src_snk_mss_AlG01_get( &
       mss_frc_src, &         ! I [frc] mass fraction of source modes (looked up)
       ovr_src_snk_mss,&      ! O [frc] mass overlap source and sink
       mss_frc_trn_dst_src, & ! O [frc] fraction of transported dust in each transport bin
       ovr_src_snk_mss_ttl)   ! O [frc] total fraction of produced dust transported
    
    !Purpose Calculate the overlap between source and sink based on mass fraction of each source mode.
    !Code by Alf Grini, UCI/UiO, 2002
    use pmgrid,only:plon,plond        !longitude number
    use dstgrd                        !dst dimensions and sizes
    use dstpsd,only:ovr_src_snk_frc  !Fractional overlap source and sink
    
    implicit none
    
    ! OUTPUT
    real(r8), intent(out):: ovr_src_snk_mss(plond,dst_src_nbr,dst_nbr) ! O [frc] mass overlap source sink
    real(r8), intent(out):: ovr_src_snk_mss_ttl(plond)                 ! O [frc] total overlap source and sink
    real(r8), intent(out):: mss_frc_trn_dst_src(plond,dst_nbr)         ! O [frc] fracion of transported dust in sink bins
    ! INPUT
    real(r8), intent(in) :: mss_frc_src(plond,dst_src_nbr)             ! I [frc] mass fraction of each source mode
    !LOCAL
    integer              :: lon_idx              ![idx] index for longitude
    integer              :: src_idx              ![idx] index for source modes
    integer              :: snk_idx              ![idx] index for dust sinks
    
    ovr_src_snk_mss_ttl(:)=0.0_r8   ! [frc]
    mss_frc_trn_dst_src(:,:)=0.0_r8 ! [frc] Fraction of transported dust mass at source
    do snk_idx=1,dst_nbr
       do src_idx=1,dst_src_nbr
          do lon_idx=1,plon
             ovr_src_snk_mss(lon_idx,src_idx,snk_idx)= & ! [frc]
                  ovr_src_snk_frc(src_idx,snk_idx)*mss_frc_src(lon_idx,src_idx) ! [frc]
             mss_frc_trn_dst_src(lon_idx,snk_idx)= & ! [frc] Fraction of transported dust mass at source
                  mss_frc_trn_dst_src(lon_idx,snk_idx)+ovr_src_snk_mss(lon_idx,src_idx,snk_idx)
             ovr_src_snk_mss_ttl(lon_idx)= & ! [frc]
                  ovr_src_snk_mss_ttl(lon_idx)+ovr_src_snk_mss(lon_idx,src_idx,snk_idx)
          enddo                ! end loo  over lon
       end do                  ! end loop over src
    end do                     ! end loop over snk
    ! Convert fraction of mobilized mass to fraction of transported mass
    do lon_idx=1,plon
       mss_frc_trn_dst_src(lon_idx,:)= & ! [frc] Fraction of transported dust mass at source
            mss_frc_trn_dst_src(lon_idx,:) &
            /max(1.e-20,ovr_src_snk_mss_ttl(lon_idx))
    enddo
    
    return
  end subroutine ovr_src_snk_mss_AlG01_get
  
  subroutine flx_mss_vrt_dst_ttl_AlG01_get( &
       dst_slt_flx_rat_ttl_lut, &  ! I [frc] look up table for alpha (function of wind and soil type)
       dst_slt_flx_rat_ttl, &      ! O [m-1] vertical to horizontal dust flux (alpha)
       flx_mss_hrz_slt_ttl, &      ! I [kg m s-1] horizontal dust flux
       flx_mss_vrt_dst_ttl, &      ! O [kg m-2 s-1] total vertical dust flux
       wnd_frc, &                  ! I [m/s] wind friction speed
       sfc_frc_bln)             ! I [frc] fraction of soil types
    
    !Purpose: Calculate the vertical flux of dust based on horizontal flux of soil.
    !We look up the ratio based on parameterization of Alfaro/Gomes 2001
    !Code by Alf Grini, UCI/UiO, 2002
    
    use pmgrid,only:plond,plon
    use dstgrd
    
    implicit none
    
    ! OUTPUT
    real(r8), intent(out) :: dst_slt_flx_rat_ttl(plond) ! O [m-1] Total horizontal soil flux
    real(r8), intent(out) :: flx_mss_vrt_dst_ttl(plond) ! O [kg m-2 s-1] total vertical dust mass flux
    ! INPUT
    real(r8), intent(in)  :: flx_mss_hrz_slt_ttl(plond)                         ! I [kg m-1 s-1] total horizontal mass flux of soil
    real(r8), intent(in)  :: dst_slt_flx_rat_ttl_lut(wnd_frc_nbr,bln_nbr) ! I [m-1] ratio of vertical to horizontal flux
    real(r8), intent(in)  :: wnd_frc(plond)                   ! I [m/s] wind friction speed
    real(r8), intent(in)  :: sfc_frc_bln(plond,bln_nbr) ! [frc]
    !LOCAL
    integer               :: soil_idx              ![idx] index for soil type
    integer               :: wnd_frc_idx           ![idx] index for wind friction speed
    integer               :: lon_idx               ![idx] index for longitude
    
    ! INITIALIZE
    dst_slt_flx_rat_ttl(:)=0.d0
    
    do soil_idx=1,bln_nbr
       do lon_idx=1,plon
          wnd_frc_idx=min(100,int(100.d0*wnd_frc(lon_idx)+0.5d0)) !wind friction speed to look up
          
          dst_slt_flx_rat_ttl(lon_idx)= &                         !ratio equals
               dst_slt_flx_rat_ttl(lon_idx) &                     !old ratio 
               +sfc_frc_bln(lon_idx,soil_idx) &                   !plus contribution from each soiltype
               *dst_slt_flx_rat_ttl_lut(wnd_frc_idx,soil_idx)     !And their respective ratio at that wind speed
          
          
       enddo  !loop on longitude
    enddo     !loop on latitude
    
    do lon_idx=1,plon
       flx_mss_vrt_dst_ttl(lon_idx)= &           !Vertical dust flux equals
            flx_mss_hrz_slt_ttl(lon_idx) &       !horizontal soil flux times
            *dst_slt_flx_rat_ttl(lon_idx)        !Ratio just looked up above
    enddo
    
    return
  end subroutine flx_mss_vrt_dst_ttl_AlG01_get

  !***************************************************************
  !January 2002 : Another new AlG01 subroutine
  subroutine flx_mss_hrz_slt_ttl_AlG01_get(   &
       flx_mss_hrz_slt_ttl                    & !O [kg m-1 s-1] horizontal saltating flux
       ,flx_mss_hrz_slt_ttl_lut               & !I [kg m-1 s-1] lut for horizontal flux
       ,sfc_frc_bln                           & !I [frc] lut for soil blends
       ,wnd_frc                               & !I [m s-1] wind friction speed
       )

    !**************************************************************************************'
    !Purpose: Look up "complex" horizontal saltation flux. The complex flux is calculated offline
    !in the sltsbl module which is now accessible by cvs at dust
    !The complex flux takes into account the real microphysics of saltation where different
    !soil sizes interact with wind friction speed and soil size distributions.
    !Lookuptables are initialized in dst_slt_sbl_cmn_ini
    !*************************************************************************************

    !Author: Alf Grini: alf.grini@geofysikk.uio.no

    use precision
    use pmgrid, only:plon,plond
    use dstgrd       !We need dust grid sizes; Here in particular: bln_nbr
    implicit none
    
    !Input
    real(r8), intent(in)            :: flx_mss_hrz_slt_ttl_lut(wnd_frc_nbr,bln_nbr) ![kg m-1 s-1]Look up table for horizontal flux
    real(r8), intent(in)            :: sfc_frc_bln(plond,bln_nbr)                   ![frc] Surface fraction of soil blends
    real(r8), intent(in)            :: wnd_frc(plond)                               ![m s-1] Wind friction speed

    !Output
    real(r8), intent(out)           :: flx_mss_hrz_slt_ttl(plond)                   ![kg m-1 s-1] Horizontal saltation flux

    !Local
    integer                         :: bln_idx    !Counter for soil blends
    integer                         :: i          !Counter for longitude
    integer                         :: wnd_frc_idx ![idx] Integer for wind friction speed (Actually wnd_frc in cm/s)

    !Initialize:
    flx_mss_hrz_slt_ttl(:)=0.0_r8
    
    do bln_idx=1,bln_nbr
       do i=1,plon

          !Get number from 1-100 describing wind friction speed
          wnd_frc_idx=min(100,nint(wnd_frc(i)*100.0_r8))  

          !Add up horizontal saltation flux
          flx_mss_hrz_slt_ttl(i)=flx_mss_hrz_slt_ttl(i)      & !Horizontal flux equals
               +flx_mss_hrz_slt_ttl_lut(wnd_frc_idx,bln_idx) & !The looked up horizontal flux for a soil type
               *sfc_frc_bln(i,bln_idx)                        !weighted by the fraction of that soil type at this longitude
          !write(6,*)'looking up ',i,bln_idx
          !write(6,*)flx_mss_hrz_slt_ttl(i)
          !write(6,*)flx_mss_hrz_slt_ttl_lut(wnd_frc_idx,bln_idx)
          !write(6,*)sfc_frc_bln(i,bln_idx)
          
       enddo
    enddo
    
    !write(6,*)'lut  :',flx_mss_hrz_slt_ttl_lut
    !stop
    
    end subroutine flx_mss_hrz_slt_ttl_AlG01_get

  !**************************************************************
  !END OF SUBROUTINES USING ALG01 PARAMETERIZATION
  !**************************************************************
  
  subroutine flx_mss_vrt_dst_prt( &
       flg_mbl,             & ! I [flg] Mobilization candidate flag
       flx_mss_vrt_dst,     & ! O [kg m-2 s-1] Vertical mass flux of dust
       flx_mss_vrt_dst_ttl, & ! I [kg m-2 s-1] Total vertical mass flux of dust
       ovr_src_snk_mss)
    ! Purpose: Partition total vertical mass flux of dust into transport bins
    ! Theory: 
    ! ++alfgr Edited 28.01.03 to use an overlap fraction calculated in dstmbl
    ! and not the overlap fraction which is saved in the dstaer module
    use dstgrd ! [mdl] Dust grid sizes
    !++alfgr use dstaer,only:ovr_src_snk_mss ! [mdl] Aerosol microphysical properties
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    logical,intent(in)::flg_mbl(plond)    ! I [flg] Mobilization candidate flag
    real(r8),intent(in)::flx_mss_vrt_dst_ttl(plond) ! I [kg m-2 s-1] Total vertical mass flux of dust
    real(r8),intent(in)::ovr_src_snk_mss(plond,dst_src_nbr,dst_nbr) ![frc] Overlap fractions 
    ! Output
    real(r8),intent(out)::flx_mss_vrt_dst(plond,dst_nbr) ! O [kg m-2 s-1] Vertical mass flux of dust
    integer rcd               ! [rcd] Return success code
    ! Local
    integer lon_idx           ! [idx] Counting index for lon
    integer src_idx           ! [idx] Counting index for src
    integer snk_idx           ! [idx] Counting index for snk
    integer snk_nbr           ! [nbr] Dimension size
    
    ! Initialize some variables
    flx_mss_vrt_dst=0.0_r8       ! [frc]
    ! Recall that snk_nbr=dst_nbr
    do lon_idx=1,plon         ! NB: Inefficient loop order
       if (flg_mbl(lon_idx)) then
          do snk_idx=1,dst_nbr
             do src_idx=1,dst_src_nbr
                flx_mss_vrt_dst(lon_idx,snk_idx)= & ! [kg m-2 s-1] 
                     flx_mss_vrt_dst(lon_idx,snk_idx)+ &
                     ovr_src_snk_mss(lon_idx,src_idx,snk_idx)* &
                     flx_mss_vrt_dst_ttl(lon_idx) 
             end do            ! end loop over src
          end do               ! end loop over snk
       endif                  ! endif flg_mbl
    end do                     ! end loop over lon
    
    return
  end subroutine flx_mss_vrt_dst_prt                       ! end flx_mss_vrt_dst_prt()
  
  subroutine vai_lsm_get( & ! [fnc] Update time-varying surface fields
       doy,                 & ! I [day] Day of year [1.0..366.0)
       lat_dgr,             & ! I [dgr] Latitude
       sfc_typ,             & ! I [idx] LSM surface type (0..28)
       vai_lsm)             ! O [m2 m-2] Vegetation area index, one-sided
    ! Purpose: Update time-varying surface fields to current timestep
    ! Routine currently computes latitude slice of VAI based on surface type, date, and latitude
    ! Algorithm is from CCM:lsm/phenol()
    ! vai_lsm_get() is not currently called, but would be called by CCM:dynamics/advnce(), MATCH:src/main()
    ! Note that routine currently returns leaf area index, not total vegetation area index = leaf + stem
    use pmgrid                ! [mdl] Spatial resolution parameters
    use dstlsm                ! [mdl] LSM data
    implicit none
    ! Parameters
    ! Input
    integer,intent(in)::sfc_typ(plond)    ! [idx] LSM surface type (0..28)
    real(r8),intent(in)::lat_dgr              ! [dgr] Latitude
    real(r8),intent(in)::doy                  ! [day] Day of year [1.0..366.0)
    ! Output
    real(r8),intent(out)::vai_lsm(plond)        ! [m2 m-2] Vegetation area index, one-sided
    ! Input/Output
    ! Local
    integer idx_mth_glb       ! [idx] Interpolation month, future
    integer idx_mth_lub       ! [idx] Interpolation month, past
    integer lon_idx           ! [idx] Counting index
    integer pln_typ_idx       ! [idx] Plant type index
    integer sfc_typ_idx       ! [idx] Surface type index
    integer sgs_idx           ! Surface sub-gridscale index
    logical::flg_SH_adj=.true. ! [flg] Add 182 days in southern hemisphere
    real(r8) wgt_glb          ! [frc] Interpolation weight
    real(r8) wgt_lub          ! [frc] Interpolation weight
    
    ! Main Code
    ! Sanity check
    if (doy < 1.0_r8 .or. doy >= 366.0_r8) stop 'dst: ERROR doy < 1.0.or.doy >= 366.0 in vai_lsm_get()'      
    ! Initialize array
    vai_lsm(:)=0.0_r8               ! [m2 m-2]
    
    call tm_2_idx_wgt(  & ! [fnc] Day of year --> bounding indices, weights
         doy,                 & ! I [day] Day of year [1.0..366.0)
         flg_SH_adj,          & ! I [flg] Add 182 days in southern hemisphere
         lat_dgr,             & ! I [dgr] Latitude
         idx_mth_glb,         & ! O [idx] Interpolation month, future
         idx_mth_lub,         & ! O [idx] Interpolation month, past
         wgt_glb,             & ! O [frc] Interpolation weight
         wgt_lub)             ! O [frc] Interpolation weight
    
    do lon_idx=1,plon
       ! Store surface blend of current gridpoint
       sfc_typ_idx=sfc_typ(lon_idx)
       do sgs_idx=1,3
          ! Non-vegetated surfaces have pln_typ=14
          pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
          ! Current VAI is weighted past VAI...
          vai_lsm(lon_idx)=   & ! [m2 m-2]
               vai_lsm(lon_idx)+ & ! [m2 m-2]
               pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc]
               tai(pln_typ_idx,idx_mth_glb)* & ! [m2 m-2]
               wgt_glb        ! [frc] 
          ! ...plus weighted future VAI
          vai_lsm(lon_idx)=   & ! [m2 m-2]
               vai_lsm(lon_idx)+ & ! [m2 m-2]
               pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc] 
               tai(pln_typ_idx,idx_mth_lub)* & ! [m2 m-2] 
               wgt_lub        ! [frc]
       end do                 ! end loop over number of possible plant types
    end do                    ! end loop over lon
    return 
  end subroutine vai_lsm_get
  
  subroutine dst_lsm_ini( & ! [fnc] Erodible fraction of gridcell
       pln_frc,             & ! I
       pln_typ,             & ! I
       lnd_frc_dry,         & ! I
       sfc_typ,             & ! I
       sfc_dst_mbl_frc)     ! O
    ! Purpose: LSM vegetation blends and surface type fields determine time-invariant fraction of each gridcell available for deflation
    ! dst_lsm_ini() is not currently called, but used to be called by CCM:control/initext(), MATCH:src/inirun()
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    integer,intent(in)::pln_typ(0:28,3)   ! [idx] LSM plant type (1..14 = nbr_LSM_pln_typ)
    integer,intent(in)::sfc_typ(plon,plat) ! [idx] LSM surface type 
    
    real(r8),intent(in)::pln_frc(0:28,3)      ! [frc] LSM fractional coverage of corresponding plant type (sums to 1)
    real(r8),intent(in)::lnd_frc_dry(plon,plat)   ! [frc] Land surface that is not lake or wetland (by area)
    ! Output
    real(r8),intent(out)::sfc_dst_mbl_frc(plon,plat) ! [frc] Non-lake, non-wetland area
    ! Local
    integer sgs_idx           ! Surface sub-gridscale index
    integer lat_idx           ! Latitude index
    integer lon_idx           ! Longitude index
    integer sfc_typ_idx       ! Surface type index
    
    ! Set the fraction of the surface available for mobilization by dust
    do lat_idx=1,plat
       do lon_idx=1,plon
          ! LSM surface type described in Bon96 Table 5, p. 12
          ! Exclude mobilization from beneath tall canopies, lakes, wetlands, etc
          
          ! Get rid of ocean and sea ice points
          ! if (nint(oro(i)) /= 1) flx_mss_mbl_sfc(i,m)=0.0_r8
          
          ! Store surface blend of current gridpoint
          sfc_typ_idx=sfc_typ(lon_idx,lat_idx)
          ! NB: plant type (pln_typ) and plant fraction (pln_frc) arrays MUST be dimensioned
          ! (0:28,3) in order to be correctly indexed by sfc_typ.
          ! Dimensioning them (29,3) will fail unless you index them by sfc_typ+1
          ! This is because sfc_typ=0 (Ocean) is allowed by LSM
          
          ! Initialize the array
          sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 
          
          ! Include shrubs, grasses, crops, and bare ground (TeF94 p. 22900)
          do sgs_idx=1,3
             if (pln_typ(sfc_typ_idx,sgs_idx) >= 6.and. &
                  pln_typ(sfc_typ_idx,sgs_idx) <= 14) then
                sfc_dst_mbl_frc(lon_idx,lat_idx)= & ! [frc]
                     sfc_dst_mbl_frc(lon_idx,lat_idx)+ &
                     pln_frc(sfc_typ_idx,sgs_idx)
             endif            ! endif bare ground
          end do               ! end loop over number of possible plant types
          
          ! Exclude ocean
          if (sfc_typ_idx == 0) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 ! [frc]
          
          ! Exclude ice and sea ice
          if (sfc_typ_idx == 1) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 ! [frc]
          
          ! Exclude warm crops (mostly Europe)
          if (sfc_typ_idx == 26) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 ! [frc]
          
          ! Exclude wetlands
          if (sfc_typ_idx == 27) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 ! [frc]
          if (sfc_typ_idx == 28) sfc_dst_mbl_frc(lon_idx,lat_idx)=0.0_r8 ! [frc]
          
          ! Land fraction suitable for dust mobilization is known for each surface gridbox
          ! However, many gridboxes are part land and part water
          ! Renormalize sfc_dst_mbl_frc to account for fractional lakes, wetlands, and ocean
          ! This removes as candidate source regions any areas which are wet year-round
          ! GBB uses same procedure to normalize land fluxes to actual grid-cell area in CCM:lsm/lsmini.F
          sfc_dst_mbl_frc(lon_idx,lat_idx)= & ! [frc]
               sfc_dst_mbl_frc(lon_idx,lat_idx)* & ! [frc]
               lnd_frc_dry(lon_idx,lat_idx) ! [frc]
          ! sfc_dst_mbl_frc now contains the (time-invariant) surface fraction of each gridcell considered suitable for dust mobilization. 
          
          ! Sanity check
          if ((sfc_dst_mbl_frc(lon_idx,lat_idx) > 1.0_r8).or. &
               (sfc_dst_mbl_frc(lon_idx,lat_idx) < 0.0_r8)) then
             write(6,'(a,i3,a,i2,a,f6.4)')  &
                  'dst: dst_lsm_ini(): sfc_dst_mbl_frc(', &
                  lon_idx,',',lat_idx,') = ', &
                  sfc_dst_mbl_frc(lon_idx,lat_idx)
          endif               ! endif out of bounds error
          
       end do                 ! end loop over lon
    end do                    ! end loop over lat
    
    return
  end subroutine dst_lsm_ini
  
  subroutine tm_2_idx_wgt(  & ! [fnc] Day of year --> bounding indices, weights
       doy,                 & ! I [day] Day of year [1.0..366.0)
       flg_SH_adj,          & ! I [flg] Add 182 days in southern hemisphere
       lat_dgr,             & ! I [dgr] Latitude
       idx_mth_glb,         & ! O [idx] Interpolation month, future
       idx_mth_lub,         & ! O [idx] Interpolation month, past
       wgt_glb,             & ! O [frc] Interpolation weight
       wgt_lub)             ! O [frc] Interpolation weight
    ! Purpose: Given day of year (doy), return bounding indices and weights 
    ! These values are used to retrieve bracketing boundary value data stored
    ! at mid-monthly points and to interpolate to time-interpolate them to doy
    ! Algorithm is from CCM:lsm/phenol()
    ! tm_2_idx_wgt() is called by lnd_frc_mbl_get(), CCM:dynamics/advnce(), MATCH:src/main()
    ! Usage:
    ! Set flg_SH_adj = true to stagger output indices and weights by 6 months
    ! Use this when indices/weights will be applied to 12-time step arrays
    ! which are valid only for NH points and SH is assumed 6 months out of phase
    ! This occurs, for example, with LSM LAI/VAI data
    ! Note that this only affects SH output, i.e., where lat_dgr < 0.0
    ! Currently routine is only used for LSM LAI/VAI so flg_SH_adj should be true
    implicit none
    ! Parameters
    integer,parameter::mpy=12 ! Resolution of external data
    integer,parameter::dpy=365 ! Model days per year
    ! Input
    logical,intent(in)::flg_SH_adj ! [flg] Add 182 days in southern hemisphere
    real(r8),intent(in)::lat_dgr  ! [dgr] Latitude
    real(r8),intent(in)::doy      ! [day] Day of year [1.0..366.0)
    ! Output
    integer,intent(out)::idx_mth_glb ! [idx] Interpolation month, future
    integer,intent(out)::idx_mth_lub ! [idx] Interpolation month, past
    real(r8),intent(out)::wgt_glb ! [frc] Interpolation weight
    real(r8),intent(out)::wgt_lub ! [frc] Interpolation weight
    ! Input/Output
    ! Local
    integer iday_NH           ! [day] Day number Northern Hemisphere [1..365]
    integer iday_SH           ! [day] iday_NH shifted 6 mth for SH [1..365]
    integer idoy              ! [day] Current day of year [1..365]
    integer imoy              ! [mth] Current month of year [1..12]
    integer lon_idx           ! [idx] Counting index
    real(r8) day_grow             ! [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
    real(r8) dom                  ! [day] Current day of month [1.0..32.0)
    real(r8) mth_grow             ! [mth] Months since Jan (NH) or Jul (SH) (1.0..12.0)
    integer dpm(mpy)          ! [day] Days per month
    integer fdom_doy(mpy)     ! [day] Day of year of first day of month 
    integer ldom_doy(mpy)     ! [day] Day of year of last day of month
    ! NB: Initializing named data statements seems to be non-standard in F90
    data dpm     / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
    data fdom_doy/  1, 32, 60, 91,121,152,182,213,244,274,305,335/
    data ldom_doy/ 31, 59, 90,120,151,181,212,243,273,304,334,365/
    save dpm
    save fdom_doy
    save ldom_doy
    
    ! Main Code
    ! Sanity check
    if (doy < 1.0.or.doy >= 366.0_r8) stop 'dst: ERROR doy < 1.0.or.doy >= 366.0 in tm_2_idx_wgt()'                  
    
    ! Initialize output to bogus values
    idx_mth_glb=-2147483647 ! [idx] Interpolation month, future
    idx_mth_lub=-2147483647 ! [idx] Interpolation month, past
    wgt_glb=1.0e36_r8 ! [frc] Interpolation weight
    wgt_lub=1.0e36_r8 ! [frc] Interpolation weight
    
    ! Convert doy to moy and dom
    imoy=1 ! [mth] Current month of year [1..12]
    idoy=int(doy) ! [day] Current day of year [1..365]
    do while(ldom_doy(imoy) < idoy.and.imoy < mpy)
       imoy=imoy+1 ! [mth] Current month of year [1..12]
    end do                    ! end while
    dom=doy-fdom_doy(imoy)+1.0_r8 ! [day] Current day of month [1.0..32.0)
    ! Sanity check
    if (imoy < 1.or.imoy > mpy) stop 'dst: ERROR imoy < 1.or.imoy > mpy in tm_2_idx_wgt()'
    if (dom < 1.0_r8.or.dom > 32.0_r8) stop 'dst: ERROR dom < 1.0.or.dom > 32.0 in tm_2_idx_wgt()'
    
    ! Compute iday_NH and iday_SH:
    ! iday_NH = day of current year since Jan 0 [1..365]
    iday_NH=int(doy)          ! [day] Day of year [1..365]
    
    day_grow=iday_NH          ! [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
    
    ! Adjust for growing season offset in SH when arrays are relative to NH
    if (flg_SH_adj) then
       if (lat_dgr < 0.0_r8) then
          ! iday_SH = iday_NH shifted 6 mth for SH: 1->183, 183->365, 184->1, 365->182
          iday_SH=mod(iday_NH-1+dpy/2,dpy)+1 ! [day] iday_NH shifted 6 mth for SH [1..365]
          day_grow=iday_SH    ! [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
       endif                  ! endif flg_SH_adj
    end if                    ! endif
    mth_grow=12.0_r8*(day_grow-0.5_r8)/dpy ! [mth] Months since Jan (NH) or Jul (SH) (1.0..12.0)
    
    ! Compute cyclic indices into monthly data idx_mth_glb,idx_mth_lub
    ! Indices and weights for interpolating time-varying fields, e.g., LSM VAI
    idx_mth_glb=mth_grow+0.5_r8 ! [idx] Interpolation month, future
    idx_mth_lub=idx_mth_glb+1 ! [idx] Interpolation month, past
    wgt_glb=(idx_mth_glb+0.5_r8)-mth_grow ! [frc] Interpolation weight
    wgt_lub=1.0_r8-wgt_glb ! [frc] Interpolation weight
    if (idx_mth_glb < 1) idx_mth_glb=mpy ! [frc] Interpolation weight
    if (idx_mth_lub > mpy) idx_mth_lub=1 ! [frc] Interpolation weight
    
    return 
  end subroutine tm_2_idx_wgt
  
  subroutine lnd_frc_mbl_get( &
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
    ! Purpose: Return fraction of each gridcell suitable for dust mobilization
    ! lnd_frc_mbl_get() is called by dst_mbl()
    ! The date is used to obtain the time-varying vegetation cover
    ! Routine currently computes latitude slice of VAI based on surface type, date, and latitude
    ! LAI/VAI algorithm is from CCM:lsm/phenol() Bon96
    ! The LSM data are mid-month values, i.e., valid on the 15th of the month
    use blmutl,only:oro_is_lnd ! [mdl] Boundary layer meteorology driver
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstlsm ! [mdl] LSM data
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Parameters
    real(r8),parameter::vai_mbl_thr=0.30 ! [m2 m-2] VAI threshold quenching dust mobilization
    ! Input
    integer,intent(in)::sfc_typ(plond) ! I [idx] LSM surface type (0..28)
    real(r8),intent(in)::doy      ! I [day] Day of year [1.0..366.0)
    real(r8),intent(in)::lat_rdn  ! I [rdn] Latitude
    real(r8),intent(in)::lnd_frc_dry(plond) ! I [frc] Dry land fraction
    real(r8),intent(in)::oro(plond) ! I [frc] Orography
    real(r8),intent(in)::snw_frc(plond) ! I [frc] Fraction of surface covered by snow
    real(r8),intent(in)::tpt_soi(plond) ! [K] Soil temperature
    real(r8),intent(in)::tpt_soi_frz ! [K] Temperature of frozen soil
    real(r8),intent(in)::vai_dst(plond) ! I [m2 m-2] Vegetation area index, one-sided
    ! Output
    integer,intent(out)::mbl_nbr ! O [nbr] Number of mobilization candidates
    logical,intent(out)::flg_mbl(plond) ! O [flg] Mobilization candidates
    real(r8),intent(out)::lnd_frc_mbl(plond) ! O [frc] Bare ground fraction
    ! Input/Output
    ! Local
    real(DBLKIND) pi       ! [frc] 3
    integer idx_idx           ! [idx] Counting index
    integer idx_mth_glb       ! [idx] Interpolation month, future
    integer idx_mth_lub       ! [idx] Interpolation month, past
    integer lnd_idx(plond)    ! [idx] Longitude index array (land)
    integer lnd_nbr           ! [nbr] Number of land points
    integer lon_idx           ! [idx] Counting index for longitude
    integer pln_typ_idx       ! [idx] Plant type index
    integer sfc_typ_idx       ! [idx] Surface type index
    integer sgs_idx           ! [idx] Surface sub-gridscale index
    logical::flg_vai_tvbds=.true. ! [flg] Use VAI data in time-varying boundary dataset
    logical::flg_SH_adj=.true. ! [flg] Add 182 days in southern hemisphere
    real(r8) lat_dgr              ! [dgr] Latitude
    real(r8) pln_frc_mbl          ! [frc] "Bare ground" fraction of sub-gridscale cell
    real(r8) pln_frc_sgs          ! [frc] Plant fraction of current sub-gridscale cell
    real(r8) vai_sgs              ! [m2 m-2] Leaf + stem area index, one-sided
    real(r8) wgt_glb              ! [frc] Interpolation weight
    real(r8) wgt_lub              ! [frc] Interpolation weight
    ! Main Code
    ! Sanity check
    if (doy < 1.0_r8.or.doy >= 366.0_r8) stop 'dst: ERROR doy < 1.0.or.doy >= 366.0 in lnd_frc_mbl_get()'
    if (vai_mbl_thr <= 0.0_r8) stop 'dst: ERROR vai_mbl_thr <= 0.0 in lnd_frc_mbl_get()'
    ! Initialize outputs
    pi=4.0_DBLKIND*atan(1.0_DBLKIND) ! [frc] 3
    lat_dgr=180.0_r8*lat_rdn/pi  ! [dgr] Latitude
    mbl_nbr=0                 ! [nbr] Number of mobilization candidates
    do lon_idx=1,plon         ! NB: plond
       flg_mbl(lon_idx)=.false. ! [flg] Mobilization candidates
    end do                     ! end loop over lon
    lnd_frc_mbl(:)=0.0_r8           ! [frc] Bare ground fraction
    
    ! Land ahoy!
    lnd_nbr=0
    do lon_idx=1,plon
       if (oro_is_lnd(oro(lon_idx))) then
          lnd_nbr=lnd_nbr+1
          lnd_idx(lnd_nbr)=lon_idx
       end if                 ! endif
    end do                     ! end loop over lon
    
    ! Much ado about nothing
    if (lnd_nbr == 0) return
    
    if (.not.flg_vai_tvbds) then
       ! LSM monthly prescribed VAI should always be offset for SH
       flg_SH_adj=.true.      ! [flg] Add 182 days in southern hemisphere
       ! Obtain bounding indices, weights for current timestep and latitude
       call tm_2_idx_wgt(     & ! [fnc] Day of year --> bounding indices, weights
            doy,              & ! I [day] Day of year [1.0..366.0)
            flg_SH_adj,       & ! I [flg] Add 182 days in southern hemisphere
            lat_dgr,          & ! I [dgr] Latitude
            idx_mth_glb,      & ! O [idx] Interpolation month, future
            idx_mth_lub,      & ! O [idx] Interpolation month, past
            wgt_glb,          & ! O [frc] Interpolation weight
            wgt_lub)          ! O [frc] Interpolation weight
    endif                     ! not flg_vai_tvbds
    
    ! Land points
    do idx_idx=1,lnd_nbr
       lon_idx=lnd_idx(idx_idx)
       ! Store surface blend of current gridpoint
       sfc_typ_idx=sfc_typ(lon_idx)
       if (sfc_typ_idx <= 1.or. & ! Inland lakes and land ice
            sfc_typ_idx >= 27.or. & ! Wetlands
            tpt_soi(lon_idx) < tpt_soi_frz) then ! Frozen soil
          lnd_frc_mbl(lon_idx)=0.0_r8 ! [frc] Bare ground fraction
       else                   ! Normal, unfrozen land
          if (flg_vai_tvbds) then
             
             ! "bare ground" fraction of current gridcell decreases
             ! linearly from 1.0 to 0.0 as VAI increases from 0.0 to vai_mbl_thr
             lnd_frc_mbl(lon_idx)= & ! O [frc] Bare ground fraction
                  1.0_r8-min(1.0_r8,min(vai_dst(lon_idx),vai_mbl_thr)/vai_mbl_thr)
             
          else                ! not flg_vai_tvbds
             
             do sgs_idx=1,3
                ! Non-vegetated surfaces have pln_typ=14
                pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx) ! [idx]
                pln_frc_sgs=pln_frc(sfc_typ_idx,sgs_idx) ! [frc]
                ! Current VAI is weighted past VAI plus weighted future VAI
                vai_sgs=      & ! [m2 m-2] Leaf + stem area index, one-sided
                     +tai(pln_typ_idx,idx_mth_glb)*wgt_glb &
                     +tai(pln_typ_idx,idx_mth_lub)*wgt_lub
                ! "bare ground" fraction of current sub-gridscale cell decreases
                ! linearly from 1.0 to 0.0 as VAI increases from 0.0 to vai_mbl_thr
                pln_frc_mbl=  & ! [frc] "Bare ground" fraction of sub-gridscale cell
                     1.0_r8-min(1.0_r8,min(vai_sgs,vai_mbl_thr)/vai_mbl_thr)
                ! if (dbg_lvl == dbg_crr) then
                ! write(6,'(a,3(a,f6.4,a))') 'lnd_frc_mbl_get(): ',
                ! $                 'vai_sgs = ',vai_sgs,' m-2 m-2, ',
                ! $                 'pln_frc = ',pln_frc_sgs,', ',
                ! $                 'pln_frc_mbl = ',pln_frc_mbl,''
                ! endif            ! endif dbg
                lnd_frc_mbl(lon_idx)=lnd_frc_mbl(lon_idx)+pln_frc_sgs*pln_frc_mbl ! [frc]
             end do            ! end loop over number of possible plant types
          endif               ! not flg_vai_tvbds
       endif                  ! endif normal land
       
       ! Adjust for factors which constrain entire gridcell
       lnd_frc_mbl(lon_idx)=   &
            lnd_frc_mbl(lon_idx) & ! [frc] Bare ground fraction
            *lnd_frc_dry(lon_idx) & ! Account for wetlands, lake, ocean, ice
            *(1.0_r8-snw_frc(lon_idx)) ! Account for snow coverage
       
       if ((lnd_frc_mbl(lon_idx) > 1.0_r8).or.(lnd_frc_mbl(lon_idx) < 0.0_r8)) then
          write(6,'(a,i3,a,f15.8)')  &
               'dst: lnd_frc_mbl_get(): lnd_frc_mbl(', &
               lon_idx,') = ',lnd_frc_mbl(lon_idx)
          stop
       endif                  ! endif out of bounds error
       
       if (lnd_frc_mbl(lon_idx) > 0.0_r8) then
          flg_mbl(lon_idx)=.true. ! [flg] Mobilization candidates
          mbl_nbr=mbl_nbr+1   ! [nbr] Number of mobilization candidates
       endif                  ! endif mobilization candidate
       
    end do                    ! end loop over land
    
    return 
  end subroutine lnd_frc_mbl_get
  
end module dstmblutl

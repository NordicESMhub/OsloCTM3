! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/aernvr.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: aernvr.F90 contains common blocks which store environmental properties
! needed for aerosol simulations

! These variables are initialized in aer_nvr_cmn_ini()
! aer_nvr_cmn_ini() is called by aer()
! aernvr.F90 MUST have access to pmgrid.F90 to work

! Usage:
! use aernvr ! [mdl] Aerosol environmental properties

module aernvr             ! [mdl] Aerosol environmental properties
  use precision ! [mdl] Precision r8, i8, ...
  use pmgrid ! [mdl] Spatial resolution parameters
  use dstgrd ! [mdl] Dust grid sizes
  implicit none
  save ! [stt] Changes to common variables are sticky
  
  ! Environmental properties needed for simulations, initialized in aer_nvr_cmn_ini()
  ! Coordinate grid
  real(r8) lat(plat)        ! [dgr] Latitude
  real(r8) lev(plev)        ! [dgr] Midlayer pressure
  real(r8) levp(plevp)      ! [dgr] Interface pressure
  real(r8) lon(plon)        ! [dgr] Longitude
  
  ! aer_nvr_cmn is initialized in aer_nvr_cmn_ini()
  integer lchnk ! [id] Chunk identifier
  integer ncol ! [nbr] Number of atmospheric columns
  integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
  
  real(r8) asp_rat_lps ! [frc] Ellipsoidal aspect ratio
  real(r8) cld_frc(plond,plevp) ! [frc] Cloud fraction
  real(r8) cld_frc_cnv(plond,plev) ! [frc] Convective cloud fraction
  real(r8) dns_mdp(plond,plev) ! [kg m-3] Midlayer density
  real(r8) doy              ! [day] Day of year [1.0..366.0)
  real(r8) flx_LW_dwn_sfc(plond) ! [W m-2] Longwave downwelling flux at surface
  real(r8) flx_SW_abs_sfc(plond) ! [W m-2] Solar flux absorbed by ground
  real(r8) hgt_mdp(plond,plev) ! [m] Midlayer height above surface
  real(r8) hgt_ntf(plond,plevp) ! [m] Interface height above surface
  real(r8) lat_dgr(plat)    ! [dgr] Latitude
  real(r8) lat_rdn(plat)    ! [rdn] Latitude
  real(r8) lnd_frc_dry(plond) ! [frc] Dry land fraction
  real(r8) mbl_bsn_fct(plond) ! [frc] Erodibility factor
  real(r8) mno_lng(plond)   ! [m] Monin-Obukhov length
  real(r8) mpl_air(plond,plev) ! [kg m-2] Air mass path in layer
  real(r8) mss_cnc_dst(plond,plev) ! [kg m-3] Mass concentration of dust
  real(r8) mss_frc_CaCO3(plond) ! [frc] Mass fraction CaCO3 
  real(r8) mss_frc_cly(plond) ! [frc] Mass fraction clay 
  real(r8) mss_frc_slt(plond) ! [frc] Mass fraction silt 
  real(r8) mss_frc_snd(plond) ! [frc] Mass fraction sand 
  real(r8) obuf(1)          ! Output buffer
  real(r8) oro(plond)       ! [frc] Orography
  real(r8) prs_dlt(plond,plev) ! [Pa] Pressure thickness
  real(r8) prs_mdp(plond,plev) ! [Pa] Midlayer pressure 
  real(r8) prs_ntf(plond,plevp) ! [Pa] Interface pressure
  real(r8) q_H2O_cnd(plond,plev) ! [kg kg-1] Condensed H2O mixing ratio
  real(r8) q_H2O_cnd2pcp_tnd(plond,plev) ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
  real(r8) q_H2O_cnd_cnv(plond,plev) ! [kg kg-1] Condensed H2O mixing ratio in convective clouds
  real(r8) q_H2O_cnd_pcp(plond,plev) ! [kg kg-1] H2O precipitation mixing ratio
  real(r8) q_H2O_cnd_tnd(plond,plev) ! [kg kg-1 s-1] Net H2O condensate formation tendency
  real(r8) q_H2O_pcp_lqd(plond,plev) ! [kg kg-1] Rain water mixing ratio
  real(r8) q_H2O_pcp2vpr_tnd(plond,plev) ! [kg kg-1 s-1] H2O precipitation to vapor tendency
  real(r8) q_H2O_vpr(plond,plev) ! [kg kg-1] Water vapor mixing ratio
  real(r8) q_H2O_vpr2pcp_cnv_tnd(plond,plev) ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
  real(r8) q_cst(plond,plev,pcnst) ! [kg kg-1] Full constituent array
  real(r8) q_dst(plond,plev,dst_nbr) ! [kg kg-1] Dust mixing ratio
  real(r8) rxrc_chm(plond,plev,chm_nbr) ! [s-1] Pseudo first order rate coefficients
  real(r8) sfc_frc_bln(plond,bln_nbr) ! [frc] Fraction of 4 available soil types
  real(r8) snw_hgt_lqd(plond) ! [m] Equivalent liquid water snow depth
  real(r8) tm_adj           ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
  real(r8) tpt_gnd(plond)   ! [K] Ground temperature
  real(r8) tpt_ice(plond)   ! [K] Ice temperature
  real(r8) tpt_mdp(plond,plev) ! [K] Temperature
  real(r8) tpt_ptn(plond,plev) ! [K] Potential temperature
  real(r8) tpt_sfc(plond)   ! [K] Surface temperature
  real(r8) tpt_soi(plond)   ! [K] Soil temperature
  real(r8) tpt_sst(plond)   ! [K] Sea surface temperature
  real(r8) tpt_vrt(plond,plev) ! [K] Virtual temperature
  real(r8) vai_dst(plond)   ! [m2 m-2] Vegetation area index, one-sided
  real(r8) vmr_HNO3_gas(plond,plev) ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
  real(r8) vmr_NO3_aer(plond,plev) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
  real(r8) vmr_SO4_aer(plond,plev) ! [mlc mlc-1] Particulate SO4 volume mixing ratio
  real(r8) vmr_chm(plond,plev,chm_nbr) ! [mlc mlc-1] Chemical species volume mixing ratios
  real(r8) vwc_sfc(plond)   ! [m3 m-3] Volumetric water content
  real(r8) wnd_frc(plond)   ! [m s-1] Friction velocity
  real(r8) wnd_mdp(plond)   ! [m s-1] Surface layer mean wind speed
  real(r8) wnd_mrd_mdp(plond) ! [m s-1] Meridional wind component
  real(r8) wnd_znl_mdp(plond) ! [m s-1] Zonal wind component
  
contains
  
  subroutine aer_nvr_cmn_cmd_ln_dfl() ! [fnc] Initialize aer_nvr_cmn command-line fields
    ! Initialize fields in aerosol environmental properties common block aer_nvr_cmn
    ! that may be superceded by command line input
    ! aer_nvr_cmn_cmd_ln_dfl() is called by aer
    ! aer_nvr_cmn_cmd_ln_dfl() must be called before aer_nvr_cmn_ini()
    ! Some initializations in aer_nvr_cmn_ini() depend on aer_nvr_cmn_cmd_ln_dfl()
    ! Command line overrides of aer_nvr_cmn_cmd_ln_dfl() initializations are handled in main()
    use precision             ! [mdl] Precision r8, i8, ...
    use dstgrd                ! [mdl] Dust grid sizes
    use dstblm                ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use pmgrid                ! [mdl] Spatial resolution parameters
    use dstcst                ! [mdl] Physical constants for dust routines
    use dstchm                ! [mdl] Chemical properties of dust
    implicit none
    ! Parameters
    ! Local
    integer i                 ! [idx] Counting index for longitude
    integer j                 ! [idx] Counting index for latitude
    integer k                 ! [idx] Counting index for level
    integer m                 ! [idx] Counting index for constituent
    real(DBLKIND) pi          ! [frc] 3
    ! Main code
    ! Initialize defaults
    pi=4.0_DBLKIND*atan(1.0_DBLKIND)        ! [frc] 3
    
    ! Set environmental pressure and temperature for settling speed
    ! fxm: 19990901 coordinate doy with nbdate,ndcur, etc.
    doy=135.0_r8                 ! [day] Day of year [1.0..366.0)
    do j=1,plat
       lat_dgr(j)=40.0_r8        ! [dgr]
       lat(j)=lat_dgr(j)      ! [dgr]
       lat_rdn(j)=pi*lat_dgr(j)/180.0_r8 ! [rdn]
    end do                    ! end loop over lat
    do i=1,plon
       lon(i)=15.0_r8*(i-1)      ! [dgr]
    end do                    ! end loop over lon
    do k=1,plevp              ! Loop ends at plevp
       levp(k)=101325.0_r8-1000.0_r8*(plevp-k) ! [Pa]
    end do                    ! end loop over lev
    do k=1,plev
       lev(k)=0.5_r8*(levp(k)+levp(k+1)) ! [Pa]
    end do                    ! end loop over lev
    
    do k=1,plevp              ! Loop ends at plevp
       do i=1,plon
          hgt_ntf(i,k)=100.0_r8-10.0_r8*(k-1) ! [m] Interface height above surface
          prs_ntf(i,k)=levp(k) ! [Pa] Interface pressure
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    do k=1,plev
       do i=1,plon
          q_H2O_vpr(i,k)=1.76721e-2_r8 ! [kg kg-1] Water vapor mixing ratio
          prs_mdp(i,k)=lev(k) ! [Pa] Midlayer pressure 
          tpt_mdp(i,k)=300.0_r8+3.0_r8*(k-1) ! [K] Temperature
          hgt_mdp(i,k)=0.5_r8*(hgt_ntf(i,k)+hgt_ntf(i,k+1)) ! [m] Midlayer height above surface
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    do i=1,plon
       wnd_znl_mdp(i)=10.0_r8    ! [m s-1] Zonal wind component
       wnd_mrd_mdp(i)=0.0_r8     ! [m s-1] Meridional wind component
       tpt_gnd(i)=300.0_r8       ! [K] Ground temperature
       tpt_soi(i)=297.0_r8      ! [K] Soil temperature
       ! fxm: tpt_frz_pnt has not been initialized yet
       ! tpt_ice(i)=tpt_frz_pnt ! [K] Ice temperature
       tpt_ice(i)=273.15_r8      ! [K] Ice temperature
       tpt_sst(i)=297.0_r8       ! [K] Sea surface temperature
       oro(i)=1.0_r8             ! [frc] Orography
       sfc_typ(i)=2           ! [idx] LSM surface type (0..28)
       vai_dst(i)=0.0_r8         ! [m2 m-2] Vegetation area index, one-sided
    end do                    ! end loop over lon
    
    ! Set from command line
    asp_rat_lps=1.0_r8 ! [frc] Ellipsoidal aspect ratio

    return
  end subroutine aer_nvr_cmn_cmd_ln_dfl ! end aer_nvr_cmn_cmd_ln_dfl()
  
  subroutine aer_nvr_cmn_ini()
    ! Initialize aerosol environmental properties common block aer_nvr_cmn
    ! aer_nvr_cmn_ini() is called by aer
    use dstblm                ! [mdl] Boundary layer meteorology for non-vegetated land surfaces
    use dstchm                ! [mdl] Chemical properties of dust
    use dstcst                ! [mdl] Physical constants for dust routines
    use dstgrd                ! [mdl] Dust grid sizes
    use dstpsd,only:dst_psd_ini,asp_rat_lps_set ! [mdl] Dust particle size distributions
    use pmgrid                ! [mdl] Spatial resolution parameters
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Parameters
    ! Local
    integer i                 ! [idx] Counting index for longitude
    integer j                 ! [idx] Counting index for latitude
    integer k                 ! [idx] Counting index for level
    integer m                 ! [idx] Counting index for constituent
    ! Main code
    ! Initialize defaults
    
    ! Variables needed for wet deposition parameterization
    do k=1,plev
       do i=1,plon
          cld_frc(i,k)=0.0_r8    ! [frc] Cloud fraction
          cld_frc_cnv(i,k)=0.0_r8 ! [frc] Convective cloud fraction
          q_H2O_cnd(i,k)=0.0_r8  ! [kg kg-1] Condensed H2O mixing ratio
          q_H2O_cnd_cnv(i,k)=0.0_r8 ! [kg kg-1] Condensed H2O mixing ratio in convective clouds
          q_H2O_cnd2pcp_tnd(i,k)=1.0e-9_r8 ! [kg kg-1 s-1] Condensed H2O to precipitation tendency
          q_H2O_cnd_pcp(i,k)=0.0_r8 ! [kg kg-1] H2O precipitation mixing ratio
          q_H2O_cnd_tnd(i,k)=0.0_r8 ! [kg kg-1 s-1] Net H2O condensate formation tendency
          q_H2O_pcp2vpr_tnd(i,k)=1.0e-9_r8 ! [kg kg-1 s-1] H2O precipitation to vapor tendency
          q_H2O_vpr2pcp_cnv_tnd(i,k)=1.0e-9_r8 ! [kg kg-1 s-1] H2O vapor to convective precipitation tendency
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    do i=1,plon
       flx_LW_dwn_sfc(i)=350.0_r8 ! [W m-2] Longwave downwelling flux at surface
       flx_SW_abs_sfc(i)=450.0_r8 ! [W m-2] Solar flux absorbed by ground
       snw_hgt_lqd(i)=0.0_r8     ! [m] Equivalent liquid water snow depth
    end do                    ! end loop over lon
    do i=1,plon
       if (nint(oro(i)) == 0) then
          lnd_frc_dry(i)=0.0_r8  ! [frc] Dry land fraction
          mbl_bsn_fct(i)=0.0_r8  ! [frc] Erodibility factor
          mss_frc_CaCO3(i)=0.0_r8 ! [frc] Mass fraction of CaCO3
          mss_frc_cly(i)=0.0_r8  ! [frc] Mass fraction of clay 
          mss_frc_snd(i)=0.0_r8  ! [frc] Mass fraction of sand
          tpt_sfc(i)=tpt_sst(i) ! [K] Surface temperature
          vwc_sfc(i)=1.0e30_r8   ! [m3 m-3] Volumetric water content
       else if (nint(oro(i)) == 1) then
          lnd_frc_dry(i)=1.0_r8  ! [frc] Dry land fraction
          mbl_bsn_fct(i)=1.0_r8  ! [frc] Erodibility factor
          mss_frc_CaCO3(i)=0.05_r8 ! [frc] Mass fraction of CaCO3
          mss_frc_cly(i)=0.19_r8 ! [frc] Mass fraction of clay 
          mss_frc_snd(i)=0.777_r8 ! [frc] Mass fraction of sand
          tpt_sfc(i)=tpt_gnd(i) ! [K] Surface temperature
          vwc_sfc(i)=0.03_r8     ! [m3 m-3] Volumetric water content
       else if (nint(oro(i)) == 2) then
          lnd_frc_dry(i)=0.0_r8  ! [frc] Dry land fraction
          mbl_bsn_fct(i)=0.0_r8  ! [frc] Erodibility factor
          mss_frc_CaCO3(i)=0.0_r8 ! [frc] Mass fraction of CaCO3
          mss_frc_cly(i)=0.0_r8  ! [frc] Mass fraction of clay 
          mss_frc_snd(i)=0.0_r8  ! [frc] Mass fraction of sand
          tpt_sfc(i)=tpt_ice(i) ! [K] Surface temperature
          vwc_sfc(i)=1.0_r8      ! [m3 m-3] Volumetric water content
       else 
          stop
       endif
    end do                    ! end loop over lon

    ! Soil blend from Chatenet
    do i=1,plon
       sfc_frc_bln(i,1)=0.375_r8 ! ASS (dgm=125um) Alumino-silicated sand
       sfc_frc_bln(i,2)=0.625_r8 ! FS  (dgm=210um) Fine sand
       sfc_frc_bln(i,3)=0.0_r8 ! SS  (dgm=520um) Salts
       sfc_frc_bln(i,4)=0.0_r8 ! CS  (dgm=610um) Coarse Sand
       !sfc_frc_bln(i,1)=0.0_r8 ! ASS (dgm=125um) Alumino-silicated sand
       !sfc_frc_bln(i,2)=0.1_r8 ! FS  (dgm=210um) Fine sand
       !sfc_frc_bln(i,3)=0.0_r8 ! SS  (dgm=520um) Salts
       !sfc_frc_bln(i,4)=0.9_r8 ! CS  (dgm=610um) Coarse Sand
    enddo
    
    do i=1,plon
       mss_frc_slt(i)=1.0_r8-mss_frc_cly(i)-mss_frc_snd(i) ! [frc] Mass fraction silt
       wnd_mdp(i)=sqrt(wnd_znl_mdp(i)*wnd_znl_mdp(i)+wnd_mrd_mdp(i)*wnd_mrd_mdp(i)) ! [m s-1] Surface layer mean wind speed
    end do                    ! end loop over lon
    
    do k=1,plev
       do i=1,plon
          prs_dlt(i,k)=prs_ntf(i,k+1)-prs_ntf(i,k) ! [Pa] Pressure thickness
          tpt_ptn(i,k)=tpt_mdp(i,k)*(prs_ntf(i,plevp)/prs_mdp(i,k))**kappa_dry_air ! [K] Virtual temperature
          tpt_vrt(i,k)=tpt_mdp(i,k)*(1.0_r8+eps_H2O_rcp_m1*q_H2O_vpr(i,k)) ! [K] Virtual temperature
          dns_mdp(i,k)=prs_mdp(i,k)/(tpt_vrt(i,k)*gas_cst_dry_air) ! [kg m-3]
          mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2] Air mass path in layer
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             ! Initializing aerosol mixing ratio to zero ensures diagnostics in dst_mss_bdg() file are synchronized
             q_dst(i,k,m)=0.0_r8 ! [kg kg-1] Dust mixing ratio
             ! q_dst(i,k,m)=1.0e-9_r8 ! [kg kg-1] Dust mixing ratio
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    do k=1,plev
       do i=1,plon
          ! Initialize gaseous and particulate VMRs with PEM-WEST-A data from ZhC99 p. 358 Tbl. 3
          ! Initial value = 0.0 ppbv for short-lived species
          ! Gas VMRs
          vmr_chm(i,k,idx_H2O2_gas)=0.0e-09_r8 ! [mlc mlc-1] Gaseous H2O2 volume mixing ratio
          vmr_chm(i,k,idx_HNO3_gas)=0.05e-09_r8 ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
          vmr_chm(i,k,idx_HO2_gas)=0.0e-09_r8 ! [mlc mlc-1] Gaseous HO2 volume mixing ratio
          vmr_chm(i,k,idx_N2O5_gas)=0.0e-09_r8 ! [mlc mlc-1] Gaseous N2O5 volume mixing ratio
          vmr_chm(i,k,idx_NO3_gas)=0.0e-09_r8 ! [mlc mlc-1] Gaseous NO3 volume mixing ratio
          vmr_chm(i,k,idx_O3_gas)=30.0e-09_r8 ! [mlc mlc-1] Gaseous O3 volume mixing ratio
          vmr_chm(i,k,idx_OH_gas)=0.0e-09_r8 ! [mlc mlc-1] Gaseous OH volume mixing ratio
          vmr_chm(i,k,idx_SO2_gas)=1.0e-09_r8 ! [mlc mlc-1] Gaseous SO2 volume mixing ratio
          ! Particulate VMRs
          vmr_SO4_aer(i,k)=0.6e-9_r8 ! [mlc mlc-1] Particulate SO4 volume mixing ratio
          vmr_NO3_aer(i,k)=0.0e-9_r8 ! [mlc mlc-1] Particulate NO3 volume mixing ratio
          vmr_chm(i,k,idx_SO4_aer)=vmr_SO4_aer(i,k) ! [mlc mlc-1] Particulate SO4 volume mixing ratio
          vmr_chm(i,k,idx_NO3_aer)=vmr_NO3_aer(i,k) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    
    mss_cnc_dst=sum(q_dst,dim=3)*dns_mdp ! [kg m-3] Mass concentration of dust
    
    do m=1,pcnst
       do k=1,plev
          do i=1,plon
#ifdef CCM
             if (m == 1) then 
#else /* !CCM */
             if (m == pcnst) then 
#endif /* !CCM */
                ! CCM, MATCH adhere to convention that water vapor is stored as specific humidity [kg H2O vapor kg-1 moist air]
                q_cst(i,k,m)=q_H2O_vpr(i,k) ! [kg kg-1] Water vapor mixing ratio
             else
                ! CCM, MATCH adhere to convention that non-water vapor tracers are stored as dry mass mixing ratios [kg tracer kg-1 dry air]
                ! Generally speaking, non-water tracers are carried as dry mass mixing ratios everywhere but the physics routines
                q_cst(i,k,m)=q_dst(i,k,m-dst_idx_srt+1) ! [kg kg-1] Dust mixing ratio
             endif            ! endif
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    call asp_rat_lps_set(asp_rat_lps) ! [frc] Ellipsoidal aspect ratio

    return
  end subroutine aer_nvr_cmn_ini            ! end aer_nvr_cmn_ini()
  
end module aernvr

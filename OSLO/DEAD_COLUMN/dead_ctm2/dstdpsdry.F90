! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstdpsdry.F90,v 1.2 2003/07/11 14:01:26 alfgr Exp $

! Purpose: dstdpsdry.F controls mineral dust dry deposition processes

! Usage:
! use dstdpsdry ! [mdl] Dry deposition driver

! Requires dst.h for DST_MSS_BDG
#include <dst.h> /* Dust preprocessor tokens */

module dstdpsdry ! [mdl] Dry deposition driver
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_dps_dry

contains
  
  subroutine dst_dps_dry(lchnk,ncol,obuf, &
       hgt_mdp,             & ! I [m] Midpoint height above surface
       lat_idx,             & ! I [idx] Model latitude index
       mno_lng,             & ! O [m] Monin-Obukhov length
       oro,                 & ! I [frc] Orography
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Midlayer pressure
       q_H2O_vpr,           & ! I [kg kg-1] Water vapor mixing ratio
       q_dst,               & ! I/O [kg kg-1] Dust mixing ratio
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp,             & ! I [K] Temperature
       tpt_ptn_mdp,         & ! I [K] Potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_frc,             & ! O [m s-1] Surface friction velocity
       wnd_mrd_mdp,         & ! I [m s-1] Meridional wind component
       wnd_rfr,             & ! O [m s-1] Wind speed at reference height
       wnd_znl_mdp)         ! I [m s-1] Zonal wind component
    ! Driver for aerosol dry deposition
    ! dst_dps_dry() is called by CCM:physics/tphysac(), MATCH:src/physlic()
    use blmutl,only:blm_glb ! [mdl] Boundary layer meteorology driver
    use dpsdryutl,only:rss_aer_get,stk_crc_get,vlc_grv_get,rss_lmn_get ! [mdl] Dry deposition utilities
    use dstaer ! [mdl] Aerosol microphysical properties
    use dstbdg,only:bdg_aal,bdg_gam_wet ! [mdl] Mass budget diagnostics
    use dstcst ! [mdl] Physical constants for dust routines
    use dstdbg ! [mdl] Debugging information for dust model
    use dstgrd ! [mdl] Dust grid sizes
    use dstmssutl,only:dst_add_lon_lev,dst_add_lon ! [mdl] Mass budget utilities
    use dstnm ! [mdl] Nomenclature for outfld()
    use dstsfc,only:sfc_typ_get ! [mdl] 2-D surface fields on PLON x PLAT grid
#ifndef BXM
    use histout,only:outfld ! [mdl] History/archiving
#endif /* BXM */
    use pmgrid ! [mdl] Spatial resolution parameters
    use precision ! [mdl] Precision r8, i8, ...
    use sng_mdl,only:ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Parameters
    real(r8),parameter::hgt_rfr=10.0 ! [m] Reference height for deposition processes
    real(r8),parameter::wnd_min_dps=1.0 ! [m s-1] Minimum windspeed used for deposition 
    ! Input
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::hgt_mdp(plond) ! I [m] Midlayer height above surface
    real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
    real(r8),intent(in)::oro(plond) ! I [frc] Orography
    real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure 
    real(r8),intent(in)::q_H2O_vpr(plond,plev) ! I [kg kg-1] Water vapor mixing ratio
    real(r8),intent(in)::snw_hgt_lqd(plond) ! I [m] Equivalent liquid water snow depth
    real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
    real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
    real(r8),intent(in)::tpt_ptn_mdp(plond) ! I [K] Potential temperature
    real(r8),intent(in)::tpt_sfc(plond) ! I [K] Surface temperature
    real(r8),intent(in)::wnd_mrd_mdp(plond) ! I [m s-1] Meridional wind component
    real(r8),intent(in)::wnd_znl_mdp(plond) ! I [m s-1] Zonal wind component
    ! Input/Output
    real(r8),intent(inout)::q_dst(plond,plev,dst_nbr) ! I/O [kg kg-1] Dust mixing ratio
    ! Output
    real(r8),intent(out)::mno_lng(plond) ! O [m] Monin-Obukhov length
    real(r8),intent(out)::wnd_frc(plond) ! O [m s-1] Friction velocity
    real(r8),intent(out)::wnd_rfr(plond) ! O [m s-1] Wind speed at reference height
    ! Local Output
    real(r8) rss_aer(plond)       ! [s m-1] Aerodynamic resistance
    real(r8) rss_lmn(plond,dst_nbr) ! [s m-1] Quasi-laminar layer resistance
    real(r8) rss_trb(plond,dst_nbr) ! [s m-1] Resistance to turbulent deposition
    real(r8) shm_nbr(plond,dst_nbr) ! [frc] Schmidt number
    real(r8) stk_nbr(plond,dst_nbr) ! [frc] Stokes number
    real(r8) vlc_dry(plond,plev,dst_nbr) ! [m s-1] Total dry deposition velocity
#ifdef BXM
    real(r8) vlc_grv(plond,plev,dst_nbr) ! [m s-1] Settling velocity
#endif /* not BXM */
    real(r8) vlc_trb(plond,dst_nbr) ! [m s-1] Turbulent deposition velocity
    ! Local
    character(80) fl_out      ! [sng] Name of netCDF output file
    integer i                 ! [idx] Counting index
    integer k                 ! [idx] lev index
    integer m                 ! [idx] Counting index
    integer sfc_typ(plond)    ! [idx] LSM surface type (0..28)
    real(r8) dns_mdp(plond,plev)  ! [kg m-3] Midlayer density
    real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
    real(r8) rgh_mmn(plond)       ! [m] Roughness length momentum
    real(r8) snw_frc(plond)       ! [frc] Fraction of surface covered by snow
    real(r8) tm_dlt               ! [s] Dry deposition timestep
    real(r8) tpt_vrt              ! [K] Virtual temperature
    real(r8) hgt_zpd(plond)       ! [m] Zero plane displacement height
    ! GCM diagnostics
    real(r8) flx_mss_dry(plond,0:plev,dst_nbr) ! [kg m-2 s-1] Flux due to settling and turbulence (+'ve downwards)
    real(r8) flx_mss_dry_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to dry deposition (+'ve downwards)
    real(r8) flx_mss_dry_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to dry deposition (+'ve downwards)
    real(r8) flx_mss_trb_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to turbulent mix-out (+'ve downwards)
    real(r8) flx_mss_trb_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to turbulent mix-out (+'ve downwards)
    real(r8) flx_mss_grv_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to gravitational settling (+'ve downwards)
    real(r8) flx_mss_grv_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to gravitational settling (+'ve downwards)
    real(r8) mpl_dst(plond,plev,dst_nbr) ! [kg m-2] Dust mass path in layer
    real(r8) q_dst_tnd_dry(plond,plev,dst_nbr) ! [kg kg-1 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
    real(r8) q_dst_tnd_dry_ttl(plond,plev) ! [kg kg-1 s-1] Total dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
    ! Main Code
    ! Initialize
    
    ! Timesplit if desired
    tm_dlt=tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)
    
    ! Initialize fluxes and tendencies
    flx_mss_dry(:,:,:)=0.0_r8           ! [kg m-2 s-1] NB: Vertical level starts at 0
    q_dst_tnd_dry(:,:,:)=0.0_r8         ! [kg kg-1 s-1]
    flx_mss_dry_sfc(:,:)=0.0_r8       ! [kg m-2 s-1]
    flx_mss_trb_sfc(:,:)=0.0_r8       ! [kg m-2 s-1]
    flx_mss_grv_sfc(:,:)=0.0_r8       ! [kg m-2 s-1]
    
    ! Compute necessary derived fields
    do k=1,plev
       do i=1,plon
          tpt_vrt=tpt_mdp(i,k)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i,k)) ! [K] Virtual temperature
          dns_mdp(i,k)=prs_mdp(i,k)/(tpt_vrt*gas_cst_dry_air) ! [kg m-3]
          mpl_air(i,k)=prs_dlt(i,k)*grv_sfc_rcp ! [kg m-2]
       end do                 ! end loop over lon
    end do                    ! end loop over lev
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             ! Mass of dust currently in gridbox
             mpl_dst(i,k,m)=q_dst(i,k,m)*mpl_air(i,k) ! [kg m-2]
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    ! Gather necessary variables
    ! Surface type
    call sfc_typ_get( &
         lat_idx,             & ! I [idx] Latitude index
         plon,                & ! I [frc] Dimension size (normally plon)
         sfc_typ)             ! O [idx] LSM surface type (0..28)
    
    call blm_glb(             & ! Solve boundary layer meteorology on global scale
         dns_mdp(1,plev),     & ! I [kg m-3] Midlayer density
         hgt_mdp,             & ! I [m] Midlayer height above surface
         oro,                 & ! I [frc] Orography
         prs_mdp(1,plev),     & ! I [Pa] Pressure
         q_H2O_vpr(1,plev),   & ! I [kg kg-1] Specific humidity
         sfc_typ,             & ! I [idx] LSM surface type (0..28)
         snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
         tpt_mdp(1,plev),     & ! I [K] Midlayer temperature
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
    
    ! Aerodynamic resistance
    call rss_aer_get( &
         hgt_mdp,             & ! I [m] Midlayer height above surface
         hgt_zpd,             & ! I [m] Zero plane displacement height
         mno_lng,             & ! I [m] Monin-Obukhov length
         rgh_mmn,             & ! I [m] Roughness length momentum
         rss_aer,             & ! O [s m-1] Aerodynamic resistance
         wnd_frc)             ! I [m s-1] Surface friction velocity
    
    ! Gravitational settling velocity
    call vlc_grv_get( &
         dmt_vwr,             & ! I [m] Mass weighted diameter resolved
         dns_aer,             & ! I [kg m-3] Particle density
         prs_mdp,             & ! I [Pa] Pressure
         stk_crc,             & ! I [frc] Correction to Stokes settling velocity
         tpt_mdp,             & ! I [K] Temperature
         vlc_dry)             ! O [m s-1] Total dry deposition velocity
    
    ! Quasi-laminar layer resistance
    call rss_lmn_get( &
         dmt_vwr,             & ! I [m] Mass weighted diameter resolved
         dns_aer,             & ! I [kg m-3] Particle density
         dns_mdp(1,plev),     & ! I [kg m-3] Midlayer density
         prs_mdp(1,plev),     & ! I [Pa] Pressure
         rss_lmn,             & ! O [s m-1] Quasi-laminar layer resistance
         shm_nbr,             & ! O [frc] Schmidt number
         stk_crc,             & ! I [frc] Correction to Stokes settling velocity
         stk_nbr,             & ! O [frc] Stokes number
         tpt_mdp(1,plev),     & ! I [K] Temperature
         wnd_frc)             ! I [m s-1] Friction velocity    
    
    ! Lowest layer: Turbulent + Gravitational deposition
    do m=1,dst_nbr
       do i=1,plon
          rss_trb(i,m)=rss_aer(i)+rss_lmn(i,m)+rss_aer(i)*rss_lmn(i,m)*vlc_dry(i,plev,m) ! [s m-1] Resistance to turbulent deposition
          vlc_trb(i,m)=1.0_r8/rss_trb(i,m) ! [m s-1] Turbulent deposition velocity
          vlc_dry(i,plev,m)=vlc_dry(i,plev,m)+vlc_trb(i,m) ! [m s-1] Total dry deposition velocity
       end do                 ! end loop over lon
    end do                     ! end loop over sz
    
#ifdef DST_DBG
    ! Sanity checks
    do i=1,plon
       if (wnd_frc(i) < 0.0_r8) then 
          write(6,'(a,i2,a,i3,a,es8.1)')   &
               'dst_dps_dry: ERROR lat = ',lat_idx,' wnd_frc(',i,') = ',wnd_frc(i)
          write(6,*)'rgh_mmn is ',rgh_mmn(i)
          call abort
       endif                  ! endif
    end do                    ! end loop over lon
    do m=1,dst_nbr
       do i=1,plon
          if (vlc_trb(i,m) < 0.0_r8) then 
             write(6,'(a,i2,a,i3,a,i2,a,es8.1)')   &
                  'dst_dps_dry: ERROR lat = ',lat_idx,' vlc_trb(',i,',',m,') = ',vlc_trb(i,m)
             call abort
          endif               ! endif
          if (vlc_dry(i,plev,m) < 0.0_r8) then 
             write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
                  'dst_dps_dry: ERROR lat = ',lat_idx,' vlc_dry(',i,',',plev,',',m,') = ',vlc_dry(i,plev,m)
             call abort
          endif               ! endif
       end do                 ! end loop over lon
    end do                     ! end loop over sz
#endif /* not DST_DBG */
    
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             ! Downward flux of bin m from level k due to settling and turbulent mixout
             flx_mss_dry(i,k,m)=vlc_dry(i,k,m)*dns_mdp(i,k)*q_dst(i,k,m) ! [kg m-2 s-1]
             ! Do not let more dust fall than is available
             if (tm_dlt*flx_mss_dry(i,k,m) > mpl_dst(i,k,m)) then
                ! #ifdef DST_DBG
                ! Courant condition has probably been violated if this code is seen
                ! This is natural and expected if large particles and long timesteps are used
                ! Only bother printing a warning when it concerns significant amounts dust
                ! if (q_dst(i,k,m) > q_dst_sgn) write (6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)') &
                ! 'dst: dry_fxr: lat = ',lat_idx,' flx_mss_dry(',i,',',k,',',m,') = ',flx_mss_dry(i,k,m)
                ! #endif /* not DST_DBG */ 
                ! Reset to physical value
                flx_mss_dry(i,k,m)=mpl_dst(i,k,m)/tm_dlt ! [kg m-2 s-1]
             endif
             q_dst_tnd_dry(i,k,m)= & ! [kg kg-1 s-1]
                  (flx_mss_dry(i,k,m)-flx_mss_dry(i,k-1,m)) & ! NB: level 0 must equal 0.0 
                  /mpl_air(i,k)
             ! Fluxes are known, so adjust mixing ratios
             q_dst(i,k,m)=q_dst(i,k,m)- & ! [kg kg-1]
                  tm_dlt*q_dst_tnd_dry(i,k,m)
#ifdef DST_DBG
             if (q_dst(i,k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
                  'dst_dry: lat = ',lat_idx,' q_dst(',i,',',k,',',m,') = ',q_dst(i,k,m)
#endif /* not DST_DBG */
             ! Get rid of unphysical values
             ! Final bit of result is sensitive to order of arithmetic operations
             ! Prevent exceedingly small negative concentrations which may appear in long runs
             if (q_dst(i,k,m) < 0.0_r8) q_dst(i,k,m)=0.0_r8 ! [kg kg-1]
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    do m=1,dst_nbr
       do i=1,plon
          flx_mss_dry_sfc(i,m)=flx_mss_dry(i,plev,m) ! [kg m-2 s-1]
          flx_mss_trb_sfc(i,m)=vlc_trb(i,m)*flx_mss_dry_sfc(i,m)/vlc_dry(i,plev,m) ! [kg m-2 s-1]
          flx_mss_grv_sfc(i,m)=flx_mss_dry_sfc(i,m)-flx_mss_trb_sfc(i,m) ! [kg m-2 s-1]
       end do                 ! end loop over lon
    end do                    ! end loop over cst
    
    ! Integrate over all size categories for diagnostic output
    call dst_add_lon_lev(q_dst_tnd_dry,q_dst_tnd_dry_ttl)
    call dst_add_lon(flx_mss_dry_sfc,flx_mss_dry_sfc_ttl)
    call dst_add_lon(flx_mss_trb_sfc,flx_mss_trb_sfc_ttl)
    call dst_add_lon(flx_mss_grv_sfc,flx_mss_grv_sfc_ttl)
    
#ifdef BXM
    ! Diagnostic vlc_grv
    do m=1,dst_nbr
       ! fxm: Following line may generate warnings on some compilers
       ! This loop will not be executed when plev=1, and is here for completeness only
       do k=1,plev-1          ! NB: plev-1 
          do i=1,plon
             vlc_grv(i,k,m)=vlc_dry(i,k,m) ! [m s-1] Settling velocity
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    do m=1,dst_nbr
       do i=1,plon
          vlc_grv(i,plev,m)=vlc_dry(i,plev,m)-vlc_trb(i,m) ! [m s-1] Settling velocity
       end do                 ! end loop over lon
    end do                     ! end loop over sz
    ! Update netCDF file
    fl_out='aer.nc'           ! [sng] Name of netCDF output file
    call ftn_strnul(fl_out)
    call dpsdry2nc(             &
         fl_out,              & ! i [sng] Name of netCDF output file
         flx_mss_dry_sfc,     & ! I [kg m-2 s-1] Surface flux due to dry deposition (+'ve downwards)
         flx_mss_dry_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to dry deposition (+'ve downwards)
         flx_mss_grv_sfc,     & ! I [kg m-2 s-1] Surface flux due to gravitational settling (+'ve downwards)
         flx_mss_grv_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to gravitational settling (+'ve downwards)
         flx_mss_trb_sfc,     & ! I [kg m-2 s-1] Surface flux due to turbulent mix-out (+'ve downwards)
         flx_mss_trb_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to turbulent mix-out (+'ve downwards)
         hgt_zpd,             & ! I [m] Zero plane displacement height
         lat_idx,             & ! I [idx] Model latitude index
         mno_lng,             & ! I [m] Monin-Obukhov length
         q_dst_tnd_dry,       & ! I [kg kg-1 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
         q_dst_tnd_dry_ttl,   & ! I [kg kg-1 s-1] Total dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
         rgh_mmn,             & ! I [m] Roughness length momentum
         rss_aer,             & ! I [s m-1] Aerodynamic resistance
         rss_lmn,             & ! I [s m-1] Quasi-laminar layer resistance
         rss_trb,             & ! I [s m-1] Resistance to turbulent deposition
         shm_nbr,             & ! I [frc] Schmidt number
         stk_nbr,             & ! I [frc] Stokes number
         vlc_dry,             & ! I [m s-1] Total dry deposition velocity
         vlc_grv,             & ! I [m s-1] Settling velocity
         vlc_trb,             & ! I [m s-1] Turbulent deposition velocity
         wnd_frc,             & ! I [m s-1] Friction velocity
         wnd_rfr)             ! I [m s-1] Wind speed at reference height
#endif /* not BXM */
#ifdef CCM
    call outfld('DSTSFDRY',flx_mss_dry_sfc_ttl,ncol,lchnk)
    call outfld('DSTSFTRB',flx_mss_trb_sfc_ttl,ncol,lchnk)
    call outfld('DSTSFGRV',flx_mss_grv_sfc_ttl,ncol,lchnk)
    call outfld('DSTSSDRY',q_dst_tnd_dry_ttl,ncol,lchnk)
    do m=1,dst_nbr
       call outfld(flx_mss_dry_sfc_nm(m),flx_mss_dry_sfc(1,m),ncol,lchnk)
       call outfld(flx_mss_trb_sfc_nm(m),flx_mss_trb_sfc(1,m),ncol,lchnk)
       call outfld(flx_mss_grv_sfc_nm(m),flx_mss_grv_sfc(1,m),ncol,lchnk)
    end do                    ! end loop over cst
#else /* not CCM */
    call outfld('DSTSFDRY',flx_mss_dry_sfc_ttl,plond,lat_idx,obuf)
    call outfld('DSTSFTRB',flx_mss_trb_sfc_ttl,plond,lat_idx,obuf)
    call outfld('DSTSFGRV',flx_mss_grv_sfc_ttl,plond,lat_idx,obuf)
    call outfld('DSTSSDRY',q_dst_tnd_dry_ttl,plond,lat_idx,obuf)
    do m=1,dst_nbr
       call outfld(flx_mss_dry_sfc_nm(m),flx_mss_dry_sfc(1,m),plond,lat_idx,obuf)
       call outfld(flx_mss_trb_sfc_nm(m),flx_mss_trb_sfc(1,m),plond,lat_idx,obuf)
       call outfld(flx_mss_grv_sfc_nm(m),flx_mss_grv_sfc(1,m),plond,lat_idx,obuf)
    end do                    ! end loop over cst
#endif /* not CCM */
#ifdef DST_MSS_BDG
    call bdg_aal('dst_sf_dry',lat_idx,flx_mss_dry_sfc_ttl)
    call bdg_aal('dst_sf_trb',lat_idx,flx_mss_trb_sfc_ttl)
    call bdg_aal('dst_sf_grv',lat_idx,flx_mss_grv_sfc_ttl)
    call bdg_gam_wet('dst_ss_dry',lat_idx,prs_dlt,q_dst_tnd_dry_ttl)
#endif /* not DST_MSS_BDG */
    return
  end subroutine dst_dps_dry                       ! end dst_dps_dry()
  
  subroutine dpsdry2nc(             &
       fl_out,              & ! I Name of netCDF output file
       flx_mss_dry_sfc,     & ! I [kg m-2 s-1] Surface flux due to dry deposition (+'ve downwards)
       flx_mss_dry_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to dry deposition (+'ve downwards)
       flx_mss_grv_sfc,     & ! I [kg m-2 s-1] Surface flux due to gravitational settling (+'ve downwards)
       flx_mss_grv_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to gravitational settling (+'ve downwards)
       flx_mss_trb_sfc,     & ! I [kg m-2 s-1] Surface flux due to turbulent mix-out (+'ve downwards)
       flx_mss_trb_sfc_ttl, & ! I [kg m-2 s-1] Total surface flux due to turbulent mix-out (+'ve downwards)
       hgt_zpd,             & ! I [m] Zero plane displacement height
       lat_idx,             & ! I [idx] Model latitude index
       mno_lng,             & ! I [m] Monin-Obukhov length
       q_dst_tnd_dry,       & ! I [kg kg-1 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
       q_dst_tnd_dry_ttl,   & ! I [kg kg-1 s-1] Total dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
       rgh_mmn,             & ! I [m] Roughness length momentum
       rss_aer,             & ! I [s m-1] Aerodynamic resistance
       rss_lmn,             & ! I [s m-1] Quasi-laminar layer resistance
       rss_trb,             & ! I [s m-1] Resistance to turbulent deposition
       shm_nbr,             & ! I [frc] Schmidt number
       stk_nbr,             & ! I [frc] Stokes number
       vlc_dry,             & ! I [m s-1] Total dry deposition velocity
       vlc_grv,             & ! I [m s-1] Settling velocity
       vlc_trb,             & ! I [m s-1] Turbulent deposition velocity
       wnd_frc,             & ! I [m s-1] Friction velocity
       wnd_rfr)             ! I [m s-1] Wind speed at reference height
    ! Purpose: Output time-varying dry deposition properties to netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstctl ! [mdl] Control variables, routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='dpsdry2nc' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    real(r8),intent(in)::flx_mss_dry_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to dry deposition (+'ve downwards)
    real(r8),intent(in)::flx_mss_dry_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to dry deposition (+'ve downwards)
    real(r8),intent(in)::flx_mss_grv_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to gravitational settling (+'ve downwards)
    real(r8),intent(in)::flx_mss_grv_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to gravitational settling (+'ve downwards)
    real(r8),intent(in)::flx_mss_trb_sfc(plond,dst_nbr) ! [kg m-2 s-1] Surface flux due to turbulent mix-out (+'ve downwards)
    real(r8),intent(in)::flx_mss_trb_sfc_ttl(plond) ! [kg m-2 s-1] Total surface flux due to turbulent mix-out (+'ve downwards)
    real(r8),intent(in)::hgt_zpd(plond) ! [m] Zero plane displacement
    real(r8),intent(in)::mno_lng(plond) ! [m] Monin-Obukhov length
    real(r8),intent(in)::q_dst_tnd_dry(plond,plev,dst_nbr) ! [kg kg-1 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
    real(r8),intent(in)::q_dst_tnd_dry_ttl(plond,plev) ! [kg kg-1 s-1] Total dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
    real(r8),intent(in)::rgh_mmn(plond) ! [m] Roughness length momentum
    real(r8),intent(in)::rss_aer(plond) ! [s m-1] Aerodynamic resistance
    real(r8),intent(in)::rss_lmn(plond,dst_nbr) ! [s m-1] Quasi-laminar resistance
    real(r8),intent(in)::rss_trb(plond,dst_nbr) ! [s m-1] Resistance to turbulent deposition
    real(r8),intent(in)::shm_nbr(plond,dst_nbr) ! [frc] Schmidt number
    real(r8),intent(in)::stk_nbr(plond,dst_nbr) ! [frc] Stokes number
    real(r8),intent(in)::vlc_dry(plond,plev,dst_nbr) ! [m s-1] Total dry deposition velocity
    real(r8),intent(in)::vlc_grv(plond,plev,dst_nbr) ! [m s-1] Settling velocity
    real(r8),intent(in)::vlc_trb(plond,dst_nbr) ! [m s-1] Turbulent deposition velocity
    real(r8),intent(in)::wnd_frc(plond) ! [m s-1] Friction velocity
    real(r8),intent(in)::wnd_rfr(plond) ! [m s-1] Wind speed at reference height
    ! Output
    ! Local
    ! File metadata and dimension IDs
    integer cnt_lon_lev_sz_time(4) ! Count array
    integer cnt_lon_lev_time(3) ! Count array
    integer cnt_lon_sz_time(3) ! Count array
    integer cnt_lon_time(2)   ! Count array
    integer dim_lon_lev_sz_time(4) ! [enm] Dimension IDs
    integer dim_lon_lev_time(3) ! [enm] Dimension IDs
    integer dim_lon_sz_time(3) ! [enm] Dimension IDs
    integer dim_lon_time(2)   ! [enm] Dimension IDs
    integer fll_mode_old      ! Old fill mode
    integer lev_dim_id        ! [enm] Dimension ID for lev
    integer lon_dim_id        ! [enm] Dimension ID for lon
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer srt_lon_lev_sz_time(4) ! Starting index array
    integer srt_lon_lev_time(3) ! Starting index array
    integer srt_lon_sz_time(3) ! Starting index array
    integer srt_lon_time(2)   ! Starting index array
    integer sz_dim_id         ! [enm] Dimension ID for sz
    integer time_dim_id       ! [enm] Dimension ID for time
    ! Variable IDs
    integer flx_mss_dry_sfc_id ! [enm] Variable ID
    integer flx_mss_dry_sfc_ttl_id ! [enm] Variable ID
    integer flx_mss_grv_sfc_id ! [enm] Variable ID
    integer flx_mss_grv_sfc_ttl_id ! [enm] Variable ID
    integer flx_mss_trb_sfc_id ! [enm] Variable ID
    integer flx_mss_trb_sfc_ttl_id ! [enm] Variable ID
    integer hgt_zpd_id        ! [enm] Variable ID
    integer mno_lng_id        ! [enm] Variable ID
    integer q_dst_tnd_dry_id  ! [enm] Variable ID
    integer q_dst_tnd_dry_ttl_id ! [enm] Variable ID
    integer rgh_mmn_id        ! [enm] Variable ID
    integer rss_aer_id        ! [enm] Variable ID
    integer rss_lmn_id        ! [enm] Variable ID
    integer rss_trb_id        ! [enm] Variable ID
    integer shm_nbr_id        ! [enm] Variable ID
    integer stk_nbr_id        ! [enm] Variable ID
    integer vlc_dry_id        ! [enm] Variable ID
    integer vlc_grv_id        ! [enm] Variable ID
    integer vlc_trb_id        ! [enm] Variable ID
    integer wnd_frc_id        ! [enm] Variable ID
    integer wnd_rfr_id        ! [enm] Variable ID
    ! Variable data

    ! Begin netCDF output routines
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
    rcd=rcd+nf90_redef(nc_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'sz',sz_dim_id)
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
    
    if (nstep == 1) then
       ! Variable definitions
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_dry_sfc',nf90_float,dim_lon_sz_time,flx_mss_dry_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_dry_sfc_ttl',nf90_float,dim_lon_time,flx_mss_dry_sfc_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_grv_sfc',nf90_float,dim_lon_sz_time,flx_mss_grv_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_grv_sfc_ttl',nf90_float,dim_lon_time,flx_mss_grv_sfc_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_trb_sfc',nf90_float,dim_lon_sz_time,flx_mss_trb_sfc_id)
       rcd=rcd+nf90_def_var(nc_id,'flx_mss_trb_sfc_ttl',nf90_float,dim_lon_time,flx_mss_trb_sfc_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'hgt_zpd_dps',nf90_float,dim_lon_time,hgt_zpd_id)
       rcd=rcd+nf90_def_var(nc_id,'mno_lng_dps',nf90_float,dim_lon_time,mno_lng_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_dry',nf90_float,dim_lon_lev_sz_time,q_dst_tnd_dry_id)
       rcd=rcd+nf90_def_var(nc_id,'q_dst_tnd_dry_ttl',nf90_float,dim_lon_lev_time,q_dst_tnd_dry_ttl_id)
       rcd=rcd+nf90_def_var(nc_id,'rgh_mmn_dps',nf90_float,dim_lon_time,rgh_mmn_id)
       rcd=rcd+nf90_def_var(nc_id,'rss_aer',nf90_float,dim_lon_time,rss_aer_id)
       rcd=rcd+nf90_def_var(nc_id,'rss_lmn',nf90_float,dim_lon_sz_time,rss_lmn_id)
       rcd=rcd+nf90_def_var(nc_id,'rss_trb',nf90_float,dim_lon_sz_time,rss_trb_id)
       rcd=rcd+nf90_def_var(nc_id,'shm_nbr',nf90_float,dim_lon_sz_time,shm_nbr_id)
       rcd=rcd+nf90_def_var(nc_id,'stk_nbr',nf90_float,dim_lon_sz_time,stk_nbr_id)
       rcd=rcd+nf90_def_var(nc_id,'vlc_dry',nf90_float,dim_lon_lev_sz_time,vlc_dry_id)
       rcd=rcd+nf90_def_var(nc_id,'vlc_grv',nf90_float,dim_lon_lev_sz_time,vlc_grv_id)
       rcd=rcd+nf90_def_var(nc_id,'vlc_trb',nf90_float,dim_lon_sz_time,vlc_trb_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_frc_dps',nf90_float,dim_lon_time,wnd_frc_id)
       rcd=rcd+nf90_def_var(nc_id,'wnd_rfr_dps',nf90_float,dim_lon_time,wnd_rfr_id)
       ! Add english text descriptions
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'long_name','Wind speed at reference height')
       rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'long_name','Roughness length momentum')
       rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'long_name','Zero plane displacement height')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_dry_sfc_id,'long_name','Surface flux due to dry deposition')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_dry_sfc_ttl_id,'long_name','Total surface flux due to dry deposition')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_trb_sfc_id,'long_name','Surface flux due to turbulent mix-out')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_trb_sfc_ttl_id,'long_name','Total surface flux due to turbulent mix-out')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_grv_sfc_id,'long_name','Surface flux due to gravitational settling')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_grv_sfc_ttl_id,'long_name','Total surface flux due to gravitational settling')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_dry_id,'long_name','Dust tendency due to settling and turbulence')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_dry_ttl_id,'long_name',' Total dust tendency due to settling and turbulence')
       rcd=rcd+nf90_put_att(nc_id,rss_aer_id,'long_name','Aerodynamic resistance')
       rcd=rcd+nf90_put_att(nc_id,rss_lmn_id,'long_name','Laminar resistance')
       rcd=rcd+nf90_put_att(nc_id,rss_trb_id,'long_name','Resistance for turbulent deposition')
       rcd=rcd+nf90_put_att(nc_id,shm_nbr_id,'long_name','Schmidt number')
       rcd=rcd+nf90_put_att(nc_id,stk_nbr_id,'long_name','Stokes number')
       rcd=rcd+nf90_put_att(nc_id,vlc_dry_id,'long_name','Total dry deposition velocity')
       rcd=rcd+nf90_put_att(nc_id,vlc_grv_id,'long_name','Gravitational settling velocity')
       rcd=rcd+nf90_put_att(nc_id,vlc_trb_id,'long_name','Turbulent deposition velocity')
       rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'long_name','Monin-Obukhov length')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'long_name','Friction velocity')
       ! Add units
       rcd=rcd+nf90_put_att(nc_id,wnd_rfr_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,mno_lng_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,wnd_frc_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,hgt_zpd_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,'units','meter')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_dry_sfc_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_dry_sfc_ttl_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_trb_sfc_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_trb_sfc_ttl_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_grv_sfc_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,flx_mss_grv_sfc_ttl_id,'units','kilogram meter-2 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_dry_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,q_dst_tnd_dry_ttl_id,'units','kilogram kilogram-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,rss_aer_id,'units','second meter-1')
       rcd=rcd+nf90_put_att(nc_id,rss_lmn_id,'units','second meter-1')
       rcd=rcd+nf90_put_att(nc_id,rss_trb_id,'units','second meter-1')
       rcd=rcd+nf90_put_att(nc_id,shm_nbr_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,stk_nbr_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,vlc_dry_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,vlc_grv_id,'units','meter second-1')
       rcd=rcd+nf90_put_att(nc_id,vlc_trb_id,'units','meter second-1')
    else                      ! endif nstep == 1
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_dry_sfc',flx_mss_dry_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_dry_sfc_ttl',flx_mss_dry_sfc_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_trb_sfc',flx_mss_trb_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_trb_sfc_ttl',flx_mss_trb_sfc_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_grv_sfc',flx_mss_grv_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'flx_mss_grv_sfc_ttl',flx_mss_grv_sfc_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_dry',q_dst_tnd_dry_id)
       rcd=nf90_wrp_inq_varid(nc_id,'q_dst_tnd_dry_ttl',q_dst_tnd_dry_ttl_id)
       rcd=nf90_wrp_inq_varid(nc_id,'mno_lng_dps',mno_lng_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_frc_dps',wnd_frc_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rgh_mmn_dps',rgh_mmn_id)
       rcd=nf90_wrp_inq_varid(nc_id,'hgt_zpd_dps',hgt_zpd_id)
       rcd=nf90_wrp_inq_varid(nc_id,'wnd_rfr_dps',wnd_rfr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rss_aer',rss_aer_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rss_lmn',rss_lmn_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rss_trb',rss_trb_id)
       rcd=nf90_wrp_inq_varid(nc_id,'shm_nbr',shm_nbr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'stk_nbr',stk_nbr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vlc_dry',vlc_dry_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vlc_grv',vlc_grv_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vlc_trb',vlc_trb_id)
    endif                     ! endif nstep /= 1
    ! End define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    rcd=rcd+nf90_put_var(nc_id,flx_mss_dry_sfc_id,flx_mss_dry_sfc,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_dry_sfc_ttl_id,flx_mss_dry_sfc_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_grv_sfc_id,flx_mss_grv_sfc,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_grv_sfc_ttl_id,flx_mss_grv_sfc_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_trb_sfc_id,flx_mss_trb_sfc,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_trb_sfc_ttl_id,flx_mss_trb_sfc_ttl,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,hgt_zpd_id,hgt_zpd,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,mno_lng_id,mno_lng,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_dry_id,q_dst_tnd_dry,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,q_dst_tnd_dry_ttl_id,q_dst_tnd_dry_ttl,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rgh_mmn_id,rgh_mmn,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,rss_aer_id,rss_aer,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,rss_lmn_id,rss_lmn,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,rss_trb_id,rss_trb,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,shm_nbr_id,shm_nbr,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,stk_nbr_id,stk_nbr,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,vlc_dry_id,vlc_dry,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,vlc_grv_id,vlc_grv,start=srt_lon_lev_sz_time,count=cnt_lon_lev_sz_time)
    rcd=rcd+nf90_put_var(nc_id,vlc_trb_id,vlc_trb,start=srt_lon_sz_time,count=cnt_lon_sz_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_frc_id,wnd_frc,start=srt_lon_time,count=cnt_lon_time)
    rcd=rcd+nf90_put_var(nc_id,wnd_rfr_id,wnd_rfr,start=srt_lon_time,count=cnt_lon_time)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
    if (nstep == 1) then
       write (6,'(a,a45,a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': Initialized dry deposition data archive in ',fl_out(1:ftn_strlen(fl_out))
    endif                     ! endif
  end subroutine dpsdry2nc                       ! end dpsdry2nc()
  
end module dstdpsdry

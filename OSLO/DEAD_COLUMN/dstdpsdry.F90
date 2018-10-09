! Purpose: dstdpsdry.F controls mineral dust dry deposition processes
! Usage:
! use dstdpsdry ! [mdl] Dry deposition driver
!
!// ------------------------------------------------------------------
!// Rewritten as a column model for Oslo CTM3
!//
!// Amund Sovde, February 2015
!// Also rewritten to take dust mass per grid box as input, instead
!// of mass mixing ratio. The old method of using mixing ratio and
!// calculated air mass and density created small (?) inconsistencies
!// between flux diagnostics and the actual change in tracer mass.
!// ------------------------------------------------------------------
!
! Requires dst.h for DST_MSS_BDG
!#include <dst.h> /* Dust preprocessor tokens */

module dstdpsdry ! [mdl] Dry deposition driver
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_dps_dry

contains
  
  !// ------------------------------------------------------------------
  subroutine dst_dps_dry(lchnk,ncol,obuf, &
       hgt_mdp,             & ! I [m] Midpoint height above surface
       lat_idx,             & ! I [idx] Model latitude index
       lon_idx,             & ! I [idx] Model longitude index
       mno_lng,             & ! O [m] Monin-Obukhov length
       oro,                 & ! I [frc] Orography
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Midlayer pressure
       q_H2O_vpr,           & ! I [kg kg-1] Water vapor mixing ratio
       m_dst,               & ! I/O [kg] Dust mass
       airm,                & ! I [kg] Air mass
       volu,                & ! I [m3] Air volume
       area,                & ! I [m2] Grid box area
       snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       tpt_mdp,             & ! I [K] Temperature
       tpt_ptn_mdp,         & ! I [K] Potential temperature
       tpt_sfc,             & ! I [K] Surface temperature
       wnd_frc,             & ! O [m s-1] Surface friction velocity
       wnd_mrd_mdp,         & ! I [m s-1] Meridional wind component
       wnd_rfr,             & ! O [m s-1] Wind speed at reference height
       wnd_znl_mdp,         & ! I [m s-1] Zonal wind component
       drydep)                ! O [kg] Dry deposited
    !// ------------------------------------------------------------------
    ! Driver for aerosol dry deposition
    ! dst_dps_dry() is called by CCM:physics/tphysac(), MATCH:src/physlic()
    use blmutl, only: blm_glb ! [mdl] Boundary layer meteorology driver
    use dpsdryutl, only: rss_aer_get, vlc_grv_get, rss_lmn_get ! [mdl] Dry deposition utilities
    use dstaer ! [mdl] Aerosol microphysical properties
! not used in CTM3
!    use dstbdg,only:bdg_aal,bdg_gam_wet ! [mdl] Mass budget diagnostics
    use dstcst ! [mdl] Physical constants for dust routines
    use dstdbg, only: q_dst_mxm ! [mdl] Debugging information for dust model
    use dstgrd ! [mdl] Dust grid sizes
    use dstmssutl, only: dst_add_lev, dst_add_nbr ! [mdl] Mass budget utilities
    !use dstnm ! [mdl] Nomenclature for outfld()
    use dstsfc, only: sfc_typ_get ! [mdl] 2-D surface fields on PLON x PLAT grid
    use dead_history, only: outfld_1 ! [mdl] History/archiving
    use pmgrid ! [mdl] Spatial resolution parameters
    !use sng_mdl,only:ftn_strnul ! [mdl] String manipulation
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Also rewritten to take dust mass per grid box as input, instead
    !// of mass mixing ratio. The old method of using mixing ratio and
    !// calculated air mass and density created small (?) inconsistencies
    !// between flux diagnostics and the actual change in tracer mass.
    !//
    !// Amund Sovde, February 2015, October 2009
    !// ------------------------------------------------------------------
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    real(r8),parameter :: hgt_rfr=10.0    ! [m] Reference height for deposition processes
    real(r8),parameter :: wnd_min_dps=1.0 ! [m s-1] Minimum windspeed used for deposition 

    ! Input
    integer,intent(in) :: lat_idx, &  ! I [idx] Model latitude index
                          lon_idx, &  ! I [idx] Model longitude index
                          lchnk, &    ! I [id] Chunk identifier
                          ncol        ! I [nbr] Number of atmospheric columns
    real(r8),intent(in) :: obuf(*)            ! I [ptr] Output buffer
    real(r8),intent(in) :: hgt_mdp, &         ! I [m] Midlayer height above surface
                           oro, &             ! I [frc] Orography
                           prs_dlt(plev), &   ! I [Pa] Pressure thickness
                           prs_mdp(plev), &   ! I [Pa] Midlayer pressure 
                           q_H2O_vpr(plev), & ! I [kg kg-1] Water vapor mixing ratio
                           snw_hgt_lqd, &     ! I [m] Equivalent liquid water snow depth
                           tm_adj, &          ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
                           tpt_mdp(plev), &   ! I [K] Temperature
                           tpt_ptn_mdp, &     ! I [K] Potential temperature
                           tpt_sfc, &         ! I [K] Surface temperature
                           wnd_mrd_mdp, &     ! I [m s-1] Meridional wind component
                           wnd_znl_mdp, &     ! I [m s-1] Zonal wind component
                           airm(plev), &      ! I [kg] Air mass
                           volu(plev), &      ! I [m3] Air volume
                           area               ! I [m2] Grid box area

    ! Input/Output
    real(r8),intent(inout) :: m_dst(plev,dst_nbr) ! I/O [kg] Dust mass

    ! Output
    real(r8),intent(out) :: mno_lng, &  ! O [m] Monin-Obukhov length
                            wnd_frc, &  ! O [m s-1] Friction velocity
                            wnd_rfr  ! O [m s-1] Wind speed at reference height

    real(r8),intent(out) :: drydep(dst_nbr) ! I/O [kg] Dry deposited

    ! Local Output
    real(r8) :: rss_aer, &               ! [s m-1] Aerodynamic resistance
                rss_lmn(dst_nbr), &      ! [s m-1] Quasi-laminar layer resistance
                rss_trb(dst_nbr), &      ! [s m-1] Resistance to turbulent deposition
                shm_nbr(dst_nbr), &      ! [frc] Schmidt number
                stk_nbr(dst_nbr), &      ! [frc] Stokes number
                vlc_dry(plev,dst_nbr), & ! [m s-1] Total dry deposition velocity
!#ifdef BXM
!                vlc_grv(plev,dst_nbr), & ! [m s-1] Settling velocity
!#endif /* not BXM */
                vlc_trb(dst_nbr)   ! [m s-1] Turbulent deposition velocity

    ! Local
    character(len=80) :: fl_out   ! [sng] Name of netCDF output file
    integer :: k, &               ! [idx] lev index
               m, &               ! [idx] Counting index
               sfc_typ            ! [idx] LSM surface type (0..28)
    real(r8) :: dns_mdp(plev), &  ! [kg m-3] Midlayer density
                mpl_air(plev), &  ! [kg m-2] Air mass path in layer
                rgh_mmn, &        ! [m] Roughness length momentum
                snw_frc, &        ! [frc] Fraction of surface covered by snow
                tm_dlt, &         ! [s] Dry deposition timestep
                tpt_vrt, &        ! [K] Virtual temperature
                hgt_zpd, &        ! [m] Zero plane displacement height
                zarea, &          ! [m-2] inverse area
                zvolu(plev)       ! [m-3] inverse volume

    real(r8) :: q_dst(plev,dst_nbr) ! [kg kg-1] Dust mass mixing ratio

    ! GCM diagnostics
    real(r8) :: flx_mss_dry(0:plev,dst_nbr), & ! [kg m-2 s-1] Flux due to settling and turbulence (+'ve downwards)
                flx_mss_dry_sfc(dst_nbr), &    ! [kg m-2 s-1] Surface flux due to dry deposition (+'ve downwards)
                flx_mss_dry_sfc_ttl, &         ! [kg m-2 s-1] Total surface flux due to dry deposition (+'ve downwards)
                flx_mss_trb_sfc(dst_nbr), &    ! [kg m-2 s-1] Surface flux due to turbulent mix-out (+'ve downwards)
                flx_mss_trb_sfc_ttl, &         ! [kg m-2 s-1] Total surface flux due to turbulent mix-out (+'ve downwards)
                flx_mss_grv_sfc(dst_nbr), &    ! [kg m-2 s-1] Surface flux due to gravitational settling (+'ve downwards)
                flx_mss_grv_sfc_ttl, &         ! [kg m-2 s-1] Total surface flux due to gravitational settling (+'ve downwards)
                mpl_dst(plev,dst_nbr), &       ! [kg m-2] Dust mass path in layer
                m_dst_tnd_dry(plev,dst_nbr), & ! [kg m-2 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
                q_dst_tnd_dry(plev,dst_nbr), & ! [kg kg-1 s-1] Dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
                q_dst_tnd_dry_ttl(plev)        ! [kg kg-1 s-1] Total dust tendency due to settling and turbulence (+'ve when a sink to the gridcell)
    !// ------------------------------------------------------------------
    ! Main Code
    ! Initialize
    
    !// Timesplit if desired
    tm_dlt = tm_adj             ! [s] (default CCM: 2*dt, MATCH: dt)

    !// Initialize fluxes and tendencies
    flx_mss_dry(:,:)   = 0.0_r8       ! [kg m-2 s-1] NB: Vertical level starts at 0
    q_dst_tnd_dry(:,:) = 0.0_r8       ! [kg kg-1 s-1]
    flx_mss_dry_sfc(:) = 0.0_r8       ! [kg m-2 s-1]
    flx_mss_trb_sfc(:) = 0.0_r8       ! [kg m-2 s-1]
    flx_mss_grv_sfc(:) = 0.0_r8       ! [kg m-2 s-1]

    !// Inverse area [m-2] and volume [m-3]
    zarea = 1._r8 / area
    do k=1,plev
       zvolu(k) = 1._r8 / volu(k)
    end do

    !// Compute necessary derived fields
    do k=1,plev
       !// Corrected air path and density by using CTM values
       !//tpt_vrt = tpt_mdp(k)*(1.0_r8 + eps_H2O_rcp_m1 * q_H2O_vpr(k)) ! [K] Virtual temperature
       !//dns_mdp(k) = prs_mdp(k)/(tpt_vrt*gas_cst_dry_air) ! [kg m-3] Air density
       !//mpl_air(k) = prs_dlt(k)*grv_sfc_rcp               ! [kg m-2] Air path
       dns_mdp(k) = airm(k) * zvolu(k)    ! [kg m-3] Air density
       mpl_air(k) = airm(k) * zarea       ! [kg m-2] Air path
    end do !// end loop over lev

    do m=1,dst_nbr
       do k=1,plev
          ! Mass of dust currently in gridbox
          !//mpl_dst(k,m) = q_dst(k,m)*mpl_air(k) ! [kg m-2]
          mpl_dst(k,m) = m_dst(k,m) * zarea      ! [kg m-2]
          !// Dust mass mixing ratio [kg kg-1] (not really necessary)
          q_dst(k,m) = m_dst(k,m) / airm(k)
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    !// Gather necessary variables
    !// Surface type (checked ok for CTM3)
    call sfc_typ_get( &
         lat_idx,             & ! I [idx] Latitude index
         lon_idx,             & ! I [idx] Longitude index
         sfc_typ)               ! O [idx] LSM surface type (0..28)
    
    !// Solve boundary layer meteorology on global scale (checked ok for CTM3)
    call blm_glb(             & 
         dns_mdp(plev),     & ! I [kg m-3] Midlayer density
         hgt_mdp,             & ! I [m] Midlayer height above surface
         oro,                 & ! I [frc] Orography
         prs_mdp(plev),     & ! I [Pa] Pressure
         q_H2O_vpr(plev),   & ! I [kg kg-1] Specific humidity
         sfc_typ,             & ! I [idx] LSM surface type (0..28)
         snw_hgt_lqd,         & ! I [m] Equivalent liquid water snow depth
         tpt_mdp(plev),     & ! I [K] Midlayer temperature
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


    !// Aerodynamic resistance (checked ok for CTM3)
    call rss_aer_get( &
         hgt_mdp,             & ! I [m] Midlayer height above surface
         hgt_zpd,             & ! I [m] Zero plane displacement height
         mno_lng,             & ! I [m] Monin-Obukhov length
         rgh_mmn,             & ! I [m] Roughness length momentum
         rss_aer,             & ! O [s m-1] Aerodynamic resistance
         wnd_frc)             ! I [m s-1] Surface friction velocity

    !// Gravitational settling velocity (checked ok for CTM3)
    call vlc_grv_get( &
         dmt_vwr,             & ! I [m] Mass weighted diameter resolved
         dns_aer,             & ! I [kg m-3] Particle density
         prs_mdp,             & ! I [Pa] Pressure
         stk_crc,             & ! I [frc] Correction to Stokes settling velocity
         tpt_mdp,             & ! I [K] Temperature
         vlc_dry)             ! O [m s-1] Total dry deposition velocity
    !// Here the gravitational+turbulent velocity (vlc_dry(plev,m)+vlc_trb(m) above)
    !// is used to transport dust downwards.
    !vlc_dry(:,:) = 0._r8


    !// Quasi-laminar layer resistance (checked ok for CTM3)
    call rss_lmn_get( &
         dmt_vwr,             & ! I [m] Mass weighted diameter resolved
         dns_aer,             & ! I [kg m-3] Particle density
         dns_mdp(plev),     & ! I [kg m-3] Midlayer density
         prs_mdp(plev),     & ! I [Pa] Pressure
         rss_lmn,             & ! O [s m-1] Quasi-laminar layer resistance
         shm_nbr,             & ! O [frc] Schmidt number
         stk_crc,             & ! I [frc] Correction to Stokes settling velocity
         stk_nbr,             & ! O [frc] Stokes number
         tpt_mdp(plev),     & ! I [K] Temperature
         wnd_frc)             ! I [m s-1] Friction velocity    

    ! Lowest layer: Turbulent + Gravitational deposition
    do m=1,dst_nbr
       !// [s m-1] Resistance to turbulent deposition
       rss_trb(m) = rss_aer &
                    + rss_lmn(m) &
                    + rss_aer * rss_lmn(m) * vlc_dry(plev,m)
       !// [m s-1] Turbulent deposition velocity
       vlc_trb(m) = 1.0_r8 / rss_trb(m)
       !// [m s-1] Total dry deposition velocity
       vlc_dry(plev,m)=vlc_dry(plev,m)+vlc_trb(m)
    end do !// end loop over sz


#ifdef DST_DBG
    ! Sanity checks
    if (wnd_frc < 0.0_r8) then 
       write(6,'(a,i2,a,i3,a,es8.1)')   &
            'dst_dps_dry: ERROR lat = ',lat_idx,' wnd_frc(',lon_idx,') = ',wnd_frc
       write(6,*)'rgh_mmn is ',rgh_mmn
       call abort
    endif                  ! endif
    do m=1,dst_nbr
       if (vlc_trb(m) < 0.0_r8) then 
          write(6,'(a,i2,a,i3,a,i2,a,es8.1)')   &
               'dst_dps_dry: ERROR lat = ',lat_idx,' vlc_trb(',lon_idx,',',m,') = ',vlc_trb(m)
          call abort
       endif               ! endif
       if (vlc_dry(plev,m) < 0.0_r8) then 
          write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
               'dst_dps_dry: ERROR lat = ',lat_idx,' vlc_dry(',lon_idx,',',plev,',',m,') = ',vlc_dry(plev,m)
          call abort
       endif               ! endif
    end do                     ! end loop over sz
#endif /* not DST_DBG */



    !// Gravitational + turbulent settling
    do m=1,dst_nbr
       do k=1,plev
          ! Downward flux of bin m from level k due to settling and turbulent mixout

          !//flx_mss_dry(k,m) = vlc_dry(k,m)*dns_mdp(k)*q_dst(k,m) ! [kg m-2 s-1]
          !// Mass flux: v[m s-1] * mass[kg] / volume[m3]
          flx_mss_dry(k,m) = vlc_dry(k,m) * m_dst(k,m) * zvolu(k) ! [kg m-2 s-1]

          ! Do not let more dust fall than is available
          if (tm_dlt*flx_mss_dry(k,m) > mpl_dst(k,m)) then
             ! Courant condition has probably been violated if this code is seen
             ! This is natural and expected if large particles and long timesteps are used
             ! Only bother printing a warning when it concerns significant amounts dust
             ! if (q_dst(i,k,m) > q_dst_sgn) write (6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)') &
             ! 'dst: dry_fxr: lat = ',lat_idx,' flx_mss_dry(',i,',',k,',',m,') = ',flx_mss_dry(i,k,m)
             ! Reset to physical value (i.e. remove all that is in the box)
             flx_mss_dry(k,m) = mpl_dst(k,m) / tm_dlt ! [kg m-2 s-1]
          endif

          !// Find tendency
          !// Assumption: Mass only fall into level below, not past it.
          !// NB: flx_mss_dry is defined from 0:plev, and for k=0 it is 0.0
          !//q_dst_tnd_dry(k,m) = & ! [kg kg-1 s-1]
          !//     (flx_mss_dry(k,m) - flx_mss_dry(k-1,m)) / mpl_air(k)
          !// Tendency is mass out minus mass coming in from above
          m_dst_tnd_dry(k,m) = & ! [kg m-2 s-1]
               (flx_mss_dry(k,m) - flx_mss_dry(k-1,m))

          !// Fluxes are known, so adjust mixing ratios [kg kg-1]
          !//q_dst(k,m) = q_dst(k,m) - q_dst_tnd_dry(k,m) * tm_dlt

          !// Have flux [kg m-2 s-1]; adjust dust mass [kg]
          m_dst(k,m) = m_dst(k,m) - m_dst_tnd_dry(k,m) * area * tm_dlt

          !// Skip test on q_dst (which is not calculated)
!#ifdef DST_DBG
          !if (q_dst(k,m) > q_dst_mxm) write(6,'(a,i2,a,i3,a,i2,a,i2,a,es8.1)')  &
          !     'dst_dry: lat = ',lat_idx,' q_dst(',lon_idx,',',k,',',m,') = ',q_dst(k,m)
!#endif /* not DST_DBG */

          ! Get rid of unphysical values
          ! Final bit of result is sensitive to order of arithmetic operations
          ! Prevent exceedingly small negative concentrations which may appear in long runs
          if (m_dst(k,m) < 0.0_r8) m_dst(k,m) = 0.0_r8 ! [kg kg-1]
       end do !// end loop over lev
    end do !// end loop over cst



    do m=1,dst_nbr
       flx_mss_dry_sfc(m) = flx_mss_dry(plev,m) ! [kg m-2 s-1]
       !// Note that dry deposition to ground is usually larger than
       !// the change in m_dst, because some dust has fallen into
       !// bottom layer from above.

       flx_mss_trb_sfc(m) = vlc_trb(m)*flx_mss_dry_sfc(m)/vlc_dry(plev,m) ! [kg m-2 s-1]
       flx_mss_grv_sfc(m) = flx_mss_dry_sfc(m)-flx_mss_trb_sfc(m) ! [kg m-2 s-1]

       !// Save for global CTM scavenging diagnostic [kg]
       drydep(m) = flx_mss_dry_sfc(m) * area * tm_dlt

    end do                    ! end loop over cst



    !// Integrate over all size categories for diagnostic output
    call dst_add_lev(q_dst_tnd_dry,q_dst_tnd_dry_ttl)
    call dst_add_nbr(flx_mss_dry_sfc,flx_mss_dry_sfc_ttl)
    call dst_add_nbr(flx_mss_trb_sfc,flx_mss_trb_sfc_ttl)
    call dst_add_nbr(flx_mss_grv_sfc,flx_mss_grv_sfc_ttl)


    !// Skip BXM and CCM, do only not-CCM

    !// Dry dep flux [kg m-2 s-1]
    call outfld_1('DSTSFDRY',flx_mss_dry_sfc_ttl,lon_idx,lat_idx,obuf)

    !// The following are not in outfld_1:
    !call outfld_1('DSTSFTRB',flx_mss_trb_sfc_ttl,lon_idx,lat_idx,obuf)
    !call outfld_1('DSTSFGRV',flx_mss_grv_sfc_ttl,lon_idx,lat_idx,obuf)
    !call outfld_1('DSTSSDRY',q_dst_tnd_dry_ttl,lon_idx,lat_idx,obuf)
    !do m=1,dst_nbr
    !   call outfld_1(flx_mss_dry_sfc_nm(m),flx_mss_dry_sfc(1,m),lon_idx,lat_idx,obuf)
    !   call outfld_1(flx_mss_trb_sfc_nm(m),flx_mss_trb_sfc(1,m),lon_idx,lat_idx,obuf)
    !   call outfld_1(flx_mss_grv_sfc_nm(m),flx_mss_grv_sfc(1,m),lon_idx,lat_idx,obuf)
    !end do                    ! end loop over cst

    !// skip for now
!#ifdef DST_MSS_BDG
    !call bdg_aal('dst_sf_dry',lat_idx,flx_mss_dry_sfc_ttl)
    !call bdg_aal('dst_sf_trb',lat_idx,flx_mss_trb_sfc_ttl)
    !call bdg_aal('dst_sf_grv',lat_idx,flx_mss_grv_sfc_ttl)
    !call bdg_gam_wet('dst_ss_dry',lat_idx,prs_dlt,q_dst_tnd_dry_ttl)
!#endif /* not DST_MSS_BDG */

    !// ------------------------------------------------------------------
  end subroutine dst_dps_dry                       ! end dst_dps_dry()
  !// ------------------------------------------------------------------
  
  
end module dstdpsdry

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dsttibds.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Read time-invariant boundary data sets into dust model 
! These routines depend heavily on netCDF I/O

! use dsttibds ! [mdl] Time-invariant boundary data sets

! params.h needed for resolution tokens
!//CTM3 rather uses pmgrid module
!#include <params_dust.h> /* Preprocessor tokens */ 

module dsttibds ! [mdl] Time-invariant boundary data sets
  use pmgrid, only: PCNST, PLON, PLAT, PLEV, PLOND, IPARW, JPARW, &
       XDGRD, YDGRD, XDEDG, YDEDG, AREAXY, IMAP, JMAP, IDGRD, JDGRD
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  save ! [stt] Changes to common variables are sticky

!//CTM3
  !// Read-in is done in native CTM resolution, so I have updated this
  !// routine accordingly.
  !// Amund Sovde, February 2015
  
!//CTM3 rather uses pmgrid module
  !integer,parameter::pcnst= PCNST ! [nbr] number of advected constituents
  !integer,parameter::plon= PLON ! [nbr] number of longitudes
  !integer,parameter::plat= PLAT ! [nbr] number of latitudes
  !integer,parameter::plev= PLEV ! [nbr] number of vertical levels
  !integer,parameter::plond= plon ! [nbr] slt extended domain longitude
  
contains 
  
  subroutine dst_tibds_ini( &
       fl_in, &                ! I file name
       mbl_name)               ! I mobilisation name
    ! Initialize time invariant land surface characteristics
    ! dst_tibds_ini() is called by CCM:control/initext(), MATCH:src/inirun()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dstsfc ! [mdl] 2-D surface fields on PLON x PLAT grid
    use dstgrd,only:bln_nbr ! [mdl] Dust grid sizes
    use regridding, only: E_GRID
    implicit none
#ifdef SCCM
    ! The SCCM header file latlon.h contains two floats, datalat and datalon, 
    ! which specify the lat and lon of the SCCM profile.
#include <latlon.h>
#endif /* not SCCM */
    ! Parameters
    !integer,parameter::lat_nbr_max=plat ! NB: This must change to enable SCCM to work
    !integer,parameter::lon_nbr_max=plon
    ! Input
    character(len=*),intent(in)::fl_in    ! netCDF input file
    character(len=*),intent(in)::mbl_name ! mobilisation name
    ! Output
    ! Input/Output
    ! Local
    integer area_id           ! [enm] Dimension ID for area
    integer bln_dmn_id        ! [enm] Dimension ID for blend
    integer bln_nbr_lcl       ! [nbr] Number of soil blends in file
    integer cnt_lon_lat(2) ! [nbr] Hyperslab size
    integer cnt_lon_lat_bln(3) ! [nbr] Hyperslab size
    integer dst_idx           ! [idx] Counting index
    integer idx               ! [idx] Counting index
    integer lat_dmn_id        ! [enm] Dimension ID for lat
    integer lat_id            ! [enm] Variable ID
    integer lat_idx           ! [idx] Counting index
    integer lat_nbr           ! [nbr] Dimension size
    integer latedg_id            ! [enm] Variable ID
    integer lnd_frc_dry_id    ! [enm] Variable ID
    integer lon_dmn_id        ! [enm] Dimension ID for lon
    integer lon_id            ! [enm] Variable ID
    integer lon_idx           ! [idx] Counting index
    integer lon_nbr           ! [nbr] Dimension size
    integer lonedg_id            ! [enm] Variable ID
    integer mbl_bsn_fct_id    ! [enm] Variable ID
    integer mss_frc_CaCO3_id  ! [enm] Variable ID
    integer mss_frc_cly_id    ! [enm] Variable ID
    integer mss_frc_snd_id    ! [enm] Variable ID
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer sfc_typ_id        ! [enm] Variable ID
    integer sfc_frc_bln_id    ! [enm] Variable ID
    integer srt_lon_lat(2) ! [idx] Starting offsets
    integer srt_lon_lat_bln(3) ! [idx] Starting offsets
    logical flg_err           ! [flg] Bad input data
    
    ! Entire lat and lon arrays are retrieved from disk
    real(r8) lat(JPARW)     ! Entire lat and lon arrays are retrieved from disk
    real(r8) lon(IPARW)
    real(r8) latedg(JPARW+1)! Entire lat and lon edge arrays are retrieved from disk
    real(r8) lonedg(IPARW+1)
    real(r8) area(IPARW,JPARW)
    real(r8) r2Dnative(IPARW,JPARW)
    real(r8) r3Dnative(IPARW,JPARW,bln_nbr)
    integer  i2Dnative(IPARW,JPARW)
    integer  stypecount(29), stmp(29), stypefreq(29),six(idgrd,jdgrd), n1,n2,ix,jx,smax,stot

#ifdef SCCM
    integer lat_idx_sccm      ! SCCM nearest latitude index into boundary data
    integer lon_idx_sccm      ! SCCM nearest longitude index into boundary data
    real(r8) dst_crr              ! Current least distance from gridpoint to target
    real(r8) dst_new              ! Distance of current gridpoint to target
#endif /* not SCCM */ 
    
    ! Initialize default values
    srt_lon_lat=(/1,1/)
    srt_lon_lat_bln=(/1,1,1/)
    !// CTM3: Data will be read in native resolution
    !cnt_lon_lat=(/plon,plat/)
    !cnt_lon_lat_bln=(/plon,plat,bln_nbr/)
    cnt_lon_lat=(/IPARW,JPARW/)
    cnt_lon_lat_bln=(/IPARW,JPARW,bln_nbr/)
    
    ! Read in land surface types from high resolution grid
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=rcd+nf90_open(fl_in,nf90_nowrite,nc_id)
    if (rcd /= nf90_noerr) then
       print*,'Cannot read file '//trim(fl_in)
       call nf90_err_exit(rcd,"dst_tibds_ini()")
    endif
    write(6,"(a29,1x,a)") "dsttibds: Opened netCDF input file",fl_in
    
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"bln",bln_dmn_id)
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr),"dst_tibds_ini(): inquire_dim lat")
    if (lat_nbr /= JPARW .and. lat_nbr /= 1) stop "lat_nbr /= JPARW .and. lat_nbr /= 1"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr),"dst_tibds_ini(): inquire_dim lon")
    if (lon_nbr /= IPARW .and. lon_nbr /= 1) stop "lon_nbr /= IPARW .and. lon_nbr /= 1"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bln_dmn_id,len=bln_nbr_lcl),"dst_tibds_ini(): inquire_dim bln")
    if (bln_nbr_lcl /= bln_nbr) stop "bln_nbr_lcl /= bln_nbr"
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"lon",lon_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_id)
    rcd=nf90_wrp_inq_varid(nc_id,"area",area_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lon_grd",lonedg_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lat_grd",latedg_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lnd_frc_dry",lnd_frc_dry_id)
    rcd=nf90_wrp_inq_varid(nc_id,mbl_name,mbl_bsn_fct_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_CaCO3",mss_frc_CaCO3_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_snd",mss_frc_snd_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_cly",mss_frc_cly_id)
    rcd=nf90_wrp_inq_varid(nc_id,"sfc_frc_bln",sfc_frc_bln_id)
    rcd=nf90_wrp_inq_varid(nc_id,"sfc_typ",sfc_typ_id)
    !// Get area
    rcd=rcd+nf90_get_var(nc_id,area_id,area)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var area")
    rcd=rcd+nf90_get_var(nc_id,latedg_id,latedg)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var latedg")
    rcd=rcd+nf90_get_var(nc_id,lonedg_id,lonedg)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var lonedg")
    ! Get data
    rcd=rcd+nf90_get_var(nc_id,lat_id,lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var lat")
    rcd=rcd+nf90_get_var(nc_id,lon_id,lon)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var lon")
    ! In SCCM, find nearest lat and lon indices into boundary data
#ifdef SCCM
    !//CTM3
    !//Not set up to to this when reading from native resolution.
    stop 'SCCM in dsttidbds.F90'
    
    lat_idx_sccm=1
    dst_crr=abs(datalat-lat(lat_idx_sccm))
    do lat_idx=1,lat_nbr
       dst_new=abs(datalat-lat(lat_idx))
       if(dst_new < dst_crr) then 
          lat_idx_sccm=lat_idx
          dst_crr=dst_new
       endif
    end do                    ! end loop over lat
    
    lon_idx_sccm=1
    dst_crr=abs(datalon-lon(lon_idx_sccm))
    do lon_idx=1,lon_nbr
       dst_new=abs(datalon-lon(lon_idx))
       if(dst_new < dst_crr) then 
          lon_idx_sccm=lon_idx
          dst_crr=dst_new
       endif
    end do                    ! end loop over lon
    
    ! Overwrite hyperslab indices with coordinates of SCCM column
    srt_lon_lat=(/lon_idx_sccm,1at_idx_sccm/)
#endif /* not SCCM */ 
    
    !// Read 2D lnd_frc_dry
    rcd=rcd+nf90_get_var(nc_id,lnd_frc_dry_id,r2Dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara lnd_frc_dry")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             lnd_frc_dry(lon_idx,lat_idx) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,lnd_frc_dry,XDEDG,YDEDG,PLON,PLAT,1)
       lnd_frc_dry(:,:) = lnd_frc_dry(:,:) / areaxy(:,:)
    end if

    !// Read 2D mobilisation map: Code names this mbl_bsn_fct, but
    !// it can be a different field defined by mbl_name above.
    rcd=rcd+nf90_get_var(nc_id,mbl_bsn_fct_id,r2Dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mbl_bsn_fct")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             mbl_bsn_fct(lon_idx,lat_idx) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,mbl_bsn_fct,XDEDG,YDEDG,PLON,PLAT,1)
       mbl_bsn_fct(:,:) = mbl_bsn_fct(:,:) / areaxy(:,:)
    end if


    !// Read 2D mss_frc_CaCO3
    rcd=rcd+nf90_get_var(nc_id,mss_frc_CaCO3_id,r2Dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_CaCO3")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             mss_frc_CaCO3(lon_idx,lat_idx) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,mss_frc_CaCO3,XDEDG,YDEDG,PLON,PLAT,1)
       mss_frc_CaCO3(:,:) = mss_frc_CaCO3(:,:) / areaxy(:,:)
    end if


    !// Read 2D mss_frc_snd
    rcd=rcd+nf90_get_var(nc_id,mss_frc_snd_id,r2Dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_snd")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             mss_frc_snd(lon_idx,lat_idx) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,mss_frc_snd,XDEDG,YDEDG,PLON,PLAT,1)
       mss_frc_snd(:,:) = mss_frc_snd(:,:) / areaxy(:,:)
    end if

    !// Read 2D mss_frc_cly
    rcd=rcd+nf90_get_var(nc_id,mss_frc_cly_id,r2Dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_cly")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             mss_frc_cly(lon_idx,lat_idx) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,mss_frc_cly,XDEDG,YDEDG,PLON,PLAT,1)
       mss_frc_cly(:,:) = mss_frc_cly(:,:) / areaxy(:,:)
    end if

    !// Read 3D sfc_frc_bln
    rcd=rcd+nf90_get_var(nc_id,sfc_frc_bln_id,r3Dnative,start=srt_lon_lat_bln,count=cnt_lon_lat_bln)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara sfc_frc_bln")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do idx=1,bln_nbr
          do lat_idx=1,lat_nbr
             do lon_idx=1,lon_nbr
                sfc_frc_bln(lon_idx,lat_idx,idx) = r3Dnative(lon_idx,lat_idx,idx)
             end do
          end do
       end do
    else
       !// Interpolate; must multiply with area
       do idx=1,bln_nbr
          r2Dnative(:,:) = r3Dnative(:,:,idx) * area(:,:)
          call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,sfc_frc_bln(:,:,idx),XDEDG,YDEDG,PLON,PLAT,1)
          sfc_frc_bln(:,:,idx) = sfc_frc_bln(:,:,idx) / areaxy(:,:)
       end do
    end if

    !// Read 2D sfc_typ (this is an integer)
    rcd=rcd+nf90_get_var(nc_id,sfc_typ_id,i2dnative,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara sfc_typ")
    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             sfc_typ(lon_idx,lat_idx) = i2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Not possible to interpolate surface types; need to pick one.
       do lat_idx = 1, PLAT
          do lon_idx = 1, PLON
             stypecount(:) = 0
             six(:,:) = 0
             do jx = 1, jdgrd
                do ix = 1, idgrd
                   !// Index is 0:28, so fortran needs 1:29
                   idx = i2Dnative(imap(ix,lon_idx),jmap(jx,lat_idx))
                   six(ix,jx) = idx ! sfc types in window
                   stypecount(idx+1) = stypecount(idx+1) + 1 ! count #
                end do
             end do
             !// sort type frequencies
             stypefreq(:) = 0
             stmp(:) = stypecount(:)
             stot = 0
             do n1 = 1, 29 !// Loop through surface types 0-28 (start on 1)
                smax = 0
                do n2 = 1, 29 !// Find index for type occuring most frequently
                   !// There may be several of same freq; if so, pick the first
                   if (stmp(n2) .gt. smax) smax = n2
                end do
                !// Keep track of most frequent surface type.
                stypefreq(n1) = smax
                !// Need to zero this type so it is not counted next time
                if (smax .gt. 0) then
                   stmp(smax) = 0 !// This entry is counted
                   stot = stot + 1
                else
                   exit !// No more data
                end if
             end do

             !// How many different entries
             if (stot .eq. 1) then
                sfc_typ(lon_idx,lat_idx) = six(1,1) !// All are the same
             else if (stot .gt. 1) then
                !// Pick the largest. Could consider area, but they are not that different.
                sfc_typ(lon_idx,lat_idx) = stypefreq(1) - 1
             else
                print*,'STOT zero'
                stop
             end if
             !print*,lon_idx,lat_idx,stot,sfc_typ(lon_idx,lat_idx)

          end do !// do lon_idx=1,PLON
       end do !// do lat_idx=1,PLAT
    end if

    rcd=rcd+nf90_close(nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"Unable to close file "//fl_in)
    write (6,"(a25,1x,a)") "dsttibds: Ingested netCDF file",fl_in

    
    ! Sanity check (can save lots of wasted time if boundary dataset is corrupt)
    !// These are all on CTM grid (possibly degraded)
    flg_err=.false.           ! [flg] Bad input data
    do lat_idx=1,plat
       do lon_idx=1,plon
          if (sfc_typ(lon_idx,lat_idx) < 0.or.sfc_typ(lon_idx,lat_idx) > 28) flg_err=.true. ! [idx] LSM surface type (0..28)
          if (lnd_frc_dry(lon_idx,lat_idx) < 0.0.or.lnd_frc_dry(lon_idx,lat_idx) > 1.0_r8) flg_err=.true. ! [frc] Dry land fraction
          if (mbl_bsn_fct(lon_idx,lat_idx) < 0.0) flg_err=.true. ! [frc] Erodibility factor
          if (mss_frc_CaCO3(lon_idx,lat_idx) < 0.0.or.mss_frc_CaCO3(lon_idx,lat_idx) > 1.0_r8) flg_err=.true. ! [frc] Mass fraction of CaCO3
          if (mss_frc_cly(lon_idx,lat_idx) < 0.0.or.mss_frc_cly(lon_idx,lat_idx) > 1.0_r8) flg_err=.true. ! [frc] Mass fraction of clay
          if (mss_frc_snd(lon_idx,lat_idx) < 0.0.or.mss_frc_snd(lon_idx,lat_idx) > 1.0_r8) flg_err=.true. ! [frc] Mass fraction of sand
          ! fxm: put this in loop up to bln_nbr
          if (sfc_frc_bln(lon_idx,lat_idx,1) < 0.0.or.sfc_frc_bln(lon_idx,lat_idx,1) > 1.0_r8) flg_err=.true. ![frc] frac. ASS
          if (sfc_frc_bln(lon_idx,lat_idx,2) < 0.0.or.sfc_frc_bln(lon_idx,lat_idx,2) > 1.0_r8) flg_err=.true. ![frc] frac. FS
          if (sfc_frc_bln(lon_idx,lat_idx,3) < 0.0.or.sfc_frc_bln(lon_idx,lat_idx,3) > 1.0_r8) flg_err=.true. ![frc] frac. SS
          if (sfc_frc_bln(lon_idx,lat_idx,4) < 0.0.or.sfc_frc_bln(lon_idx,lat_idx,4) > 1.0_r8) flg_err=.true. ![frc] frac. CS
          if (flg_err) then
             write(6,"(a,2(a4,i3,a4,f9.4,a2))")  &
                  "dsttibds: ERROR dst_tibds_ini() sanity check failed at ", &
                  "lat(",lat_idx,") = ",ydgrd(lat_idx),", ", &
                  "lon(",lon_idx,") = ",xdgrd(lon_idx),", "
             write(6,'(a,i3)') 'sfc_typ: ',sfc_typ(lon_idx,lat_idx)
             write(6,'(a,f16.12)') 'lnd_frc_dry: ',lnd_frc_dry(lon_idx,lat_idx)
             write(6,'(a,es20.12)') 'mbl_bsn_fct: ',mbl_bsn_fct(lon_idx,lat_idx)
             write(6,'(a,es20.12)') 'mss_frc_CaCO3: ',mss_frc_CaCO3(lon_idx,lat_idx)
             write(6,'(a,es20.12)') 'mss_frc_cly: ',mss_frc_cly(lon_idx,lat_idx)
             write(6,'(a,es20.12)') 'mss_frc_sbd: ',mss_frc_snd(lon_idx,lat_idx)
             write(6,'(a,es20.12)') 'sfc_frc_bln 1: ',sfc_frc_bln(lon_idx,lat_idx,1)
             write(6,'(a,es20.12)') 'sfc_frc_bln 2: ',sfc_frc_bln(lon_idx,lat_idx,2)
             write(6,'(a,es20.12)') 'sfc_frc_bln 3: ',sfc_frc_bln(lon_idx,lat_idx,3)
             write(6,'(a,es20.12)') 'sfc_frc_bln 4: ',sfc_frc_bln(lon_idx,lat_idx,4)
             stop 
          endif               ! endif err
       end do                 ! end loop over lon
    end do                    ! end loop over lat
    return
  end subroutine dst_tibds_ini                       ! end dst_tibds_ini()
  
end module dsttibds

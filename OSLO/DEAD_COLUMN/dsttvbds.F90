! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dsttvbds.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: dstvbds.F90 contains subroutines for initializing, reading,
! interpolating, time-varying boundary data sets into the dust model. 
! These datasets are currently assumed be annual cycle datasets 
! containing twelve monthly averages.
! These routines depend heavily on netCDF I/O

! Usage:
! use dsttvbds ! [mdl] Time-varying boundary data sets

! params.h required for horizontal resolution tokens
!#include <params_dust.h> /* Preprocessor tokens */ 

module dsttvbds ! [mdl] Time-varying boundary data sets
  use pmgrid, only: PLAT, PLON, IPARW, JPARW, YDEDG, XDEDG, AREAXY
  use dead_precision ! [mdl] Precision r8, i8, ...
  use dstdbg ! [mdl] Debugging information for dust model
  implicit none
  save ! [stt] Changes to common variables are sticky
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_tvbds_ini ! [sbr] Initialize annual cycle boundary data file
  public::dst_tvbds_ntp ! [sbr] Interpolate annual cycle data to current time value
  public::dst_tvbds_close ! [sbr] Close annual cycle boundary data file
  public::dst_tvbds_get ! [sbr] Return specifed latitude slice of boundary data
  
  ! Time-varying boundary data sets
  ! fxm: these variables should be private
  integer,public,parameter::lat_nbr=PLAT
  integer,public,parameter::lon_nbr=PLON
  integer,parameter::time_nbr=12
  
  ! fxm: initialize fl_nm somewhere to improve diagnostics, currently fl_nm is never set
  character(80) fl_nm       ! [sng] netCDF file containing TVBDS data
  
  ! Variables initialized in dst_tvbds_ini()
  integer cnt_lon_lat_time(3) ! [nbr] Hyperslab size
  integer idx_glb_dsk       ! [idx] Index of most recent past time slice in input file
  integer idx_glb_ram       ! [idx] Index of most recent past time slice in memory
  integer idx_lub_dsk       ! [idx] Index of nearest future time slice in input file
  integer idx_lub_ram       ! [idx] Index of nearest future time slice in memory
  integer nc_id             ! [id] File ID
  integer srt_lon_lat_time(3) ! [idx] Starting offsets
  
  ! Variables initialized in dst_tvbds_ntp() and dst_tvbds_ini()
  real(r8),public::vai_dst(lon_nbr,lat_nbr) ! [m2 m-2] Vegetation area index, one-sided, interpolated fxm:should be private but BXM needs it
  real(r8) vai_dst_bnd(lon_nbr,lat_nbr,2) ! [m2 m-2] Vegetation area index, one-sided, boundary data
#ifdef TOMS
  real(r8),public::src_str(lon_nbr,lat_nbr) ! [frc] Source strength
  real(r8) src_str_bnd(lon_nbr,lat_nbr,2) ! [frc] Source strength, boundary data
#endif /* !TOMS */
  real(r8) time(time_nbr) ! [day] Time coordinate (day of year)
  
contains
  
  subroutine dst_tvbds_ini(fl_in,time_doy) ! I
    ! Purpose: Initialize annual cycle boundary data file
    ! NB: This routine opens fl_in and leaves it open for the duration of the run
    ! Remember to use dst_tvbds_close() to close it
    ! dst_tvbds_ini() is called by CCM:control/initext(), MATCH:src/main.F:inirun()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
    use regridding, only: E_GRID
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="dst_tvbds_ini" ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_in       ! I [sng] netCDF input file
    real(r8),intent(in)::time_doy         ! I [day] Day of year
    ! Output (common block dstbnd)
    ! Input/Output
    ! Local
    integer area_id           ! [enm] Dimension ID for area
    integer lat_dmn_id        ! [enm] Dimension ID for lat
    integer lat_nbr_in        ! [nbr] Dimension size
    integer latedg_id            ! [enm] Variable ID
    integer lon_dmn_id        ! [enm] Dimension ID for lon
    integer lon_nbr_in        ! [nbr] Dimension size
    integer lonedg_id            ! [enm] Variable ID
    integer rcd               ! [rcd] Return success code
    integer time_dmn_id       ! [enm] Dimension ID for time
    integer time_id           ! [enm] Variable ID
#ifdef TOMS
    integer src_str_id        ! [enm] Variable ID
#endif /* !TOMS */
    integer vai_dst_id        ! [enm] Variable ID
    integer time_nbr_in       ! [nbr] Dimension size

    ! Entire lat and lon arrays are retrieved from disk
    real(r8) latedg(JPARW+1)! Entire lat and lon edge arrays are retrieved from disk
    real(r8) lonedg(IPARW+1)
    real(r8) area(IPARW,JPARW)
    real(r8) r2Dnative(IPARW,JPARW)
    real(r8) r2Drun(PLON,PLAT)
    integer lat_idx,lon_idx

    ! Read in land surface types from high resolution grid
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
    write(6,"(a29,1x,a)") "dsttvbds: Opened netCDF input file",fl_in
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    rcd=nf90_wrp_inq_varid(nc_id,"area",area_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lon_grd",lonedg_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lat_grd",latedg_id)

    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr_in),"dst_tvbds_ini(): inquire_dim lat")
    if (JPARW /= lat_nbr_in) stop "IPARW /= lat_nbr_in"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr_in),"dst_tvbds_ini(): inquire_dim lon")
    if (IPARW /= lon_nbr_in) stop "lon_nbr /= lon_nbr_in"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr_in),"dst_tvbds_ini(): inquire_dim time")
    if (time_nbr /= time_nbr_in) stop "time_nbr /= time_nbr_in"
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"time",time_id)
#ifdef TOMS
    rcd=nf90_wrp_inq_varid(nc_id,"src_str_mdl",src_str_id)
#endif /* !TOMS */
    rcd=nf90_wrp_inq_varid(nc_id,"vai_ttl_clm",vai_dst_id)
    ! Get time coordinate
    rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time),"get_var time")

    ! Verify time coordinate is monotonic increasing
    if(.not.mnt_ncr_chk(time,time_nbr)) stop "dst_tvbds_ini(): time coordinate not monotonic increasing"

    ! Get indices of bounding time slices
    idx_glb_dsk=idx_glb_get(time,time_nbr,time_doy)
    idx_lub_dsk=mod(idx_glb_dsk,time_nbr)+1
    
    ! Initialize hyperslab indices for disk file
    srt_lon_lat_time=(/1,1,1/)
    !// CTM3: Data will be read in native resolution
    !cnt_lon_lat_time=(/lon_nbr,lat_nbr,1/)
    cnt_lon_lat_time=(/IPARW,JPARW,1/)
    
    ! Initialize hyperslab indices for variables in RAM
    idx_glb_ram=1
    idx_lub_ram=2

    rcd=nf90_get_var(nc_id,area_id,area)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var area")
    rcd=nf90_get_var(nc_id,latedg_id,latedg)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var latedg")
    rcd=nf90_get_var(nc_id,lonedg_id,lonedg)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var lonedg")

    ! Get glb records
    srt_lon_lat_time(3)=idx_glb_dsk
    rcd=nf90_wrp(nf90_get_var(nc_id,vai_dst_id,r2dnative(:,:), &
         start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var vai_dst")
#ifdef TOMS
    rcd=nf90_wrp(nf90_get_var(nc_id,src_str_id,src_str_bnd(:,:,idx_glb_ram), &
         start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var src_str")
#endif /* !TOMS */

    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             vai_dst_bnd(lon_idx,lat_idx,idx_glb_ram) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,r2drun,XDEDG,YDEDG,PLON,PLAT,1)
       vai_dst_bnd(:,:,idx_glb_ram) = r2drun(:,:) / areaxy(:,:)
    end if


    ! Get lub records
    srt_lon_lat_time(3)=idx_lub_dsk
    rcd=nf90_wrp(nf90_get_var(nc_id,vai_dst_id,r2dnative(:,:), &
         start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var vai_dst")
#ifdef TOMS
    rcd=nf90_wrp(nf90_get_var(nc_id,src_str_id,src_str_bnd(:,:,idx_lub_ram), &
         start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var src_str")
#endif /* !TOMS */

    !// Use native or degrade?
    if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             vai_dst_bnd(lon_idx,lat_idx,idx_lub_ram) = r2Dnative(lon_idx,lat_idx)
          end do
       end do
    else
       !// Interpolate; must multiply with area
       r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
       call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,r2drun,XDEDG,YDEDG,PLON,PLAT,1)
       vai_dst_bnd(:,:,idx_lub_ram) = r2drun(:,:) / areaxy(:,:)
    end if

    
    write (6,"(3(a,f9.4))") "dsttvbds: dst_tvbds_ini(): Day of year = ",time_doy, &
         ": read data for days",time(idx_glb_dsk)," and ",time(idx_lub_dsk)

    return
  end subroutine dst_tvbds_ini                       ! end dst_tvbds_ini()
  
  subroutine dst_tvbds_ntp(time_doy) ! I
    ! Purpose: Interpolate annual cycle data to current time value
    ! Input new data as necessary
    ! dst_tvbds_ntp() is called by CCM:dynamics/advnce(), MATCH:src/main()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use regridding, only: E_GRID
    implicit none
    ! Parameters
    ! Input
    real(r8),intent(in)::time_doy             ! I [day] Day of year
    ! Output (common block dstbnd)
    ! Input/Output
    ! Local
#ifdef TOMS
    integer src_str_id        ! [enm] Variable ID
#endif /* !TOMS */
    integer vai_dst_id        ! [enm] Variable ID
    integer area_id           ! [enm] Dimension ID for area
    integer latedg_id            ! [enm] Variable ID
    integer lonedg_id            ! [enm] Variable ID
    integer dat_nbr           ! Array size passed to interpolation routine
    integer idx_glb_dsk_old   ! [idx] Index of most recent past time slice in input file, old
    integer idx_lub_dsk_old   ! [idx] Index of nearest future time slice in input file, old
    integer rcd               ! [rcd] Return success code
    real(r8) crd_max          ! Maximum coordinate value for interpolation routine
    real(r8) crd_min          ! Minimum coordinate value for interpolation routine
    real(r8) crd_ntp          ! Arrival coordinate for interpolation routine
    real(r8) time_dlt         ! [day] Time interval between glb and lub data
    real(r8) time_dlt_1st     ! [day] Time interval between glb and current time

    ! Entire lat and lon arrays are retrieved from disk
    real(r8) latedg(JPARW+1)! Entire lat and lon edge arrays are retrieved from disk
    real(r8) lonedg(IPARW+1)
    real(r8) area(IPARW,JPARW)
    real(r8) r2Dnative(IPARW,JPARW)
    real(r8) r2Drun(PLON,PLAT)
    integer lat_idx,lon_idx

    ! Main code
    ! Vet time
    if (time_doy < 1.0_r8.or.time_doy > 366.25_r8) stop "dsttvbds: dst_tvbds_ntp() reports time_doy out of bounds"
    
    ! Save current indices of boundary time slices
    idx_lub_dsk_old=idx_lub_dsk ! [idx] Index of most recent past time slice in input file, old
    idx_glb_dsk_old=idx_glb_dsk ! [idx] Index of nearest future time slice in input file, old

    ! Get new indices of boundary time slices
    idx_glb_dsk=idx_glb_get(time,time_nbr,time_doy) ! [idx] Index of nearest future time slice in input file
    idx_lub_dsk=mod(idx_glb_dsk,time_nbr)+1 ! [idx] Index of nearest future time slice in input file
    
    ! Is model time still bounded by dataset times?
    if (idx_glb_dsk /= idx_glb_dsk_old.or.idx_lub_dsk /= idx_lub_dsk_old) then
       
       ! Get variable IDs (defensive programming, dataset may have changed...)
#ifdef TOMS
       rcd=nf90_wrp_inq_varid(nc_id,"src_str_mdl",src_str_id)
#endif /* !TOMS */
       rcd=nf90_wrp_inq_varid(nc_id,"vai_ttl_clm",vai_dst_id)

       ! Get dimension IDs
       rcd=nf90_wrp_inq_varid(nc_id,"area",area_id)
       rcd=nf90_wrp_inq_varid(nc_id,"lon_grd",lonedg_id)
       rcd=nf90_wrp_inq_varid(nc_id,"lat_grd",latedg_id)

       
       if (idx_lub_dsk /= idx_lub_dsk_old) then
          ! Read in new "future" data, replace old "past" data
          
          idx_glb_ram=idx_lub_ram ! Swap pointers to "past" and "future" slices
          idx_lub_ram=mod(idx_glb_ram,2)+1 ! "future" index is whatever past is not
          
          srt_lon_lat_time(3)=idx_lub_dsk
          rcd=nf90_wrp(nf90_get_var(nc_id,vai_dst_id,r2Dnative(:,:), &
               start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var vai_dst")
          !// Use native or degrade?
          if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
             do lat_idx=1,lat_nbr
                do lon_idx=1,lon_nbr
                   vai_dst_bnd(lon_idx,lat_idx,idx_lub_ram) = r2Dnative(lon_idx,lat_idx)
                end do
             end do
          else
             !// Get stuff for interpolation
             rcd=nf90_get_var(nc_id,area_id,area)
             if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var area")
             rcd=nf90_get_var(nc_id,latedg_id,latedg)
             if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var latedg")
             rcd=nf90_get_var(nc_id,lonedg_id,lonedg)
             if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ini():get_var lonedg")

             !// Interpolate; must multiply with area
             r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
             call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,r2drun,XDEDG,YDEDG,PLON,PLAT,1)
             vai_dst_bnd(:,:,idx_lub_ram) = r2drun(:,:) / areaxy(:,:)
          end if

#ifdef TOMS
          rcd=nf90_wrp(nf90_get_var(nc_id,src_str_id,src_str_bnd(:,:,idx_lub_ram), &
               start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var src_str")
#endif /* !TOMS */

          write (6,"(2(a,f9.4))") "dsttvbds: dst_tvbds_ntp(): Day of year = ",time_doy,": read data for day ",time(idx_lub_dsk)
          
          ! Update glb
          ! NB: This is defensive programming---this code should never need to be executed
          if (idx_glb_dsk /= idx_lub_dsk_old) then
             srt_lon_lat_time(3)=idx_glb_dsk
             rcd=nf90_wrp(nf90_get_var(nc_id,vai_dst_id,r2Dnative(:,:), &
                  start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var vai_dst")
             !// Use native or degrade?
             if (PLON.eq.IPARW .and. PLAT.eq.JPARW) then
                do lat_idx=1,lat_nbr
                   do lon_idx=1,lon_nbr
                      vai_dst_bnd(lon_idx,lat_idx,idx_glb_ram) = r2Dnative(lon_idx,lat_idx)
                   end do
                end do
             else
                !// Interpolate; must multiply with area
                r2Dnative(:,:) = r2Dnative(:,:) * area(:,:)
                call E_GRID(r2Dnative,lonedg,latedg,IPARW,JPARW,r2drun,XDEDG,YDEDG,PLON,PLAT,1)
                vai_dst_bnd(:,:,idx_glb_ram) = r2drun(:,:) / areaxy(:,:)
             end if

#ifdef TOMS
             rcd=nf90_wrp(nf90_get_var(nc_id,src_str_id,src_str_bnd(:,:,idx_glb_ram), &
                  start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var src_str")
#endif /* !TOMS */
             if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tvbds_ntp():get_vara vai_dst")
             write (6,"(2(a,f9.4))") "dsttvbds: dst_tvbds_ntp(): WARNING Day of year = ",time_doy, &
                  ": read data for day ",time(idx_glb_dsk)
             stop "dsttvbds: dst_tvbds_ntp() executing code that should not be reached"
          endif               ! endif reading new glb
       endif                  ! endif reading new lub
    endif                     ! endif updating either glb or lub
    
    ! Renormalize coordinate distance to unity to facilitate interpolation 
    if (time(idx_lub_dsk) < time(idx_glb_dsk)) then ! i.e., Dec 15 -- Jan 15
       time_dlt=365.0_r8 - time(idx_glb_dsk) + time(idx_lub_dsk)
       if (time_doy <= time(idx_lub_dsk)) then ! i.e., Jan 1 -- Jan 15
          time_dlt_1st=365.0_r8 - time(idx_glb_dsk) + time_doy
       else
          time_dlt_1st=time_doy - time(idx_glb_dsk) ! i.e., Dec 15 -- Dec 31
       end if                 ! endif
    else                      ! i.e., Jan 15 -- Dec 15
       time_dlt=time(idx_lub_dsk) - time(idx_glb_dsk)
       time_dlt_1st=time_doy - time(idx_glb_dsk)
    end if                    ! endif
    
    ! Interpolate to current time
    dat_nbr=lat_nbr*lon_nbr
    crd_min=0.0_r8
    crd_max=1.0_r8
    crd_ntp=time_dlt_1st/time_dlt
    if (crd_ntp < 0.0_r8.or.crd_ntp > 1.0_r8) stop "dsttvbds: dst_tvbds_ntp() reports crd_ntp out of bounds"
    
    call ntp_arr(dat_nbr,vai_dst_bnd(1,1,idx_glb_ram),vai_dst_bnd(1,1,idx_lub_ram), & ! I
         crd_min,crd_max,crd_ntp, & ! I
         vai_dst)             ! O
#ifdef TOMS
    call ntp_arr(dat_nbr,src_str_bnd(1,1,idx_glb_ram),src_str_bnd(1,1,idx_lub_ram), & ! I
         crd_min,crd_max,crd_ntp, & ! I
         src_str)             ! O
#endif /* !TOMS */
    
#if 0
    ! Sanity check on boundary data interpolation
    write (6,"(a,9(a,i3,a),4(a,es9.2,a))") &
         "dsttvbds: Diagnostics at from dst_tvbds_ntp(): ", &
         "lon_dbg = ",lon_dbg," [idx], ", &
         "lat_dbg = ",lat_dbg," [idx], ", &
         "lon_nbr = ",lon_nbr," [nbr], ", &
         "lat_nbr = ",lat_nbr," [nbr], ", &
         "time_nbr = ",time_nbr," [nbr], ", &
         "idx_glb_dsk = ",idx_glb_dsk," [idx], ", &
         "idx_lub_dsk = ",idx_lub_dsk," [idx], ", &
         "idx_glb_ram = ",idx_glb_ram," [idx], ", &
         "idx_lub_ram = ",idx_lub_ram," [idx], ", &
         "time_doy = ",time_doy," [day], ", &
         "vai_dst_bnd(1) = ",vai_dst_bnd(lon_dbg,lat_dbg,idx_glb_ram)," [m2 m-2], ", &
         "vai_dst_bnd(2) = ",vai_dst_bnd(lon_dbg,lat_dbg,idx_lub_ram)," [m2 m-2], ", &
         "vai_dst = ",vai_dst(lon_dbg,lat_dbg)," [m2 m-2], "
#endif /* !0 */

    return
  end subroutine dst_tvbds_ntp                       ! end dst_tvbds_ntp()
  
  subroutine dst_tvbds_close()
    ! Purpose: Close annual cycle boundary data file
    ! dst_tvbds_close() is called by CCM:dynamics/eul/stepon(), MATCH:main()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    implicit none
    ! Local
    integer rcd               ! [rcd] Return success code
    ! Main Code
    rcd=nf90_close(nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"Unable to close tvbds netCDF file")
    write (6,"(a25,1x,a)") "dsttvbds: Closed tvbds netCDF file"
    return
  end subroutine dst_tvbds_close                       ! end dst_tvbds_close()
  
  !// ------------------------------------------------------------------
  subroutine dst_tvbds_get( &
       lat_idx,             & ! I [idx] Latitude index
       lon_idx,             & ! I [idx] Longitude index
#ifdef TOMS
       src_str_out,         & ! O [frc] Source strength, current
#endif /* !TOMS */
       vai_dst_out)         ! O [m2 m-2] Vegetation area index, one-sided, current
    !// ------------------------------------------------------------------
    ! Purpose: Return specifed latitude slice of surface boundary data
    ! dst_tvbds_get() is called by dst_mbl() and used by lnd_frc_mbl_get()
    ! NB: dst_tvbds_get() does not convert arrays from netCDF
    ! dataset resolution (lat_nbr = plat X lon_nbr = plon) into latitude 
    ! slices dimensioned according to model resolution (i.e., set 
    ! lon_nbr_out = plon)
    !//
    !// Rewritten for Oslo CTM3, fetches only value for (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! Input
    integer,intent(in)::lat_idx, &        ! I [idx] Latitude index
                        lon_idx           ! I [nbr] Dimension size (normally plon)
    ! Output
#ifdef TOMS
    real(r8),intent(out)::src_str_out ! O [frc] Source strength, current
#endif /* !TOMS */
    real(r8),intent(out)::vai_dst_out ! O [m2 m-2] Vegetation area index, one-sided, current
    !// ------------------------------------------------------------------

    ! Copying data from buffer to I/O variable ensures integrity of buffer
    ! However, vai information is read only so copying it wastes time
    ! fxm: Use fortran pointers instead
#ifdef TOMS
    src_str_out = src_str(lon_idx,lat_idx) ! [frc] Source strength
#endif /* !TOMS */
    vai_dst_out = vai_dst(lon_idx,lat_idx) ! [m2 m-2] Vegetation area index, one-sided

    !// ------------------------------------------------------------------
  end subroutine dst_tvbds_get                       ! end dst_tvbds_get()
  !// ------------------------------------------------------------------


  
  subroutine ntp_arr(dat_nbr,dat_min,dat_max, & ! I
       crd_min,crd_max,crd_ntp, & ! I
       dat_out)             ! O
    ! Purpose: Interpolate the input arrays to the given coordinate
    ! ntp_arr() is called by dst_tvbds_ntp()
    ! Input
    integer,intent(in)::dat_nbr           ! I [nbr] Dimension size
    real(r8),intent(in)::crd_max              ! I [frc] Maximum coordinate value for interpolation routine
    real(r8),intent(in)::crd_min              ! I [frc] Minimum coordinate value for interpolation routine
    real(r8),intent(in)::crd_ntp              ! I [frc] Arrival coordinate for interpolation routine
    real(r8),intent(in)::dat_max(dat_nbr)     ! I [frc] Array values at crd_max
    real(r8),intent(in)::dat_min(dat_nbr)     ! I [frc] Array values at crd_min
    ! Output
    real(r8),intent(out)::dat_out(dat_nbr)     ! O [frc] Array values at crd_ntp
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) crd_dlt              ! [frc] Interval between min and max coordinates
    real(r8) crd_dlt_1st          ! [frc] Interval between min and arrival coordinates
    real(r8) crd_fct              ! [frc] Factor crd_dlt_1st/crd_dlt
    
    if (crd_max <= crd_min) then
       write (6,"(a)") "dsttvbds: ERROR bnt_ntp() reports crd_max <= crd_min"
       stop
    endif                     ! endif err
    
    ! Compute interpolation weight outside main loop
    crd_dlt=crd_max-crd_min ! [frc] Interval between min and max coordinates
    crd_dlt_1st=crd_ntp-crd_min ! [frc] Interval between min and arrival coordinates
    crd_fct=crd_dlt_1st/crd_dlt ! [frc] Factor crd_dlt_1st/crd_dlt
    
    do idx=1,dat_nbr
       dat_out(idx)=dat_min(idx)+crd_fct*(dat_max(idx)-dat_min(idx))
    end do                    ! end loop over idx
    
    return
  end subroutine ntp_arr                       ! end ntp_arr()
  
  integer function idx_glb_get(dat,dat_nbr,val)
    ! Purpose: Return index of largest data point (in dat) which is less than specified value (val), 
    ! i.e., find index of greatest lower bound of specified value.
    ! Data array is assumed to be cyclic so routine should return dat_nbr when val < dat(1)
    ! idx_glb_get() is called by dst_tvbds_ini() and dst_tvbds_ntp()
    ! Input
    integer,intent(in)::dat_nbr           ! I [nbr] Dimension size
    real(r8),intent(in)::dat(dat_nbr)     ! I [frc] Monotonically increasing input array
    real(r8),intent(in)::val              ! I [frc] Value to bound
    ! Local
    integer dat_idx           ! [idx] Counting index
    
    if (val < dat(1).or.val >= dat(dat_nbr)) then
       idx_glb_get=dat_nbr
       return
    endif                     ! endif
    do dat_idx=2,dat_nbr
       if (val < dat(dat_idx)) then
          idx_glb_get=dat_idx-1
          return
       endif                  ! endif
    end do                    ! end loop over dat
    write (6,"(a)") "dsttvbds: ERROR idx_glb_get() unable to find greatest lower bound"
    stop
  end function idx_glb_get                       ! end idx_glb_get()
  
end module dsttvbds

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dsttibds.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Read time-invariant boundary data sets into dust model 
! These routines depend heavily on netCDF I/O

! use dsttibds ! [mdl] Time-invariant boundary data sets

! params.h needed for resolution tokens
#include <params.h> /* Preprocessor tokens */ 

module dsttibds ! [mdl] Time-invariant boundary data sets
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  save ! [stt] Changes to common variables are sticky
  
  integer,parameter::pcnst= PCNST ! [nbr] number of advected constituents
  integer,parameter::plon= PLON ! [nbr] number of longitudes
  integer,parameter::plat= PLAT ! [nbr] number of latitudes
  integer,parameter::plev= PLEV ! [nbr] number of vertical levels
  integer,parameter::plond= plon ! [nbr] slt extended domain longitude
  
contains 
  
  subroutine dst_tibds_ini( &
       fl_in)               ! I
    ! Initialize time invariant land surface characteristics
    ! dst_tibds_ini() is called by CCM:control/initext(), MATCH:src/inirun()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use precision ! [mdl] Precision r8, i8, ...
    use dstsfc ! [mdl] 2-D surface fields on PLON x PLAT grid
    use dstgrd,only:bln_nbr ! [mdl] Dust grid sizes
    implicit none
#ifdef SCCM
    ! The SCCM header file latlon.h contains two floats, datalat and datalon, 
    ! which specify the lat and lon of the SCCM profile.
#include <latlon.h>
#endif /* not SCCM */
    ! Parameters
    integer,parameter::lat_nbr_max=plat ! NB: This must change to enable SCCM to work
    integer,parameter::lon_nbr_max=plon
    ! Input
    character(len=*),intent(in)::fl_in ! netCDF input file
    ! Output
    ! Input/Output
    ! Local
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
    integer lnd_frc_dry_id    ! [enm] Variable ID
    integer lon_dmn_id        ! [enm] Dimension ID for lon
    integer lon_id            ! [enm] Variable ID
    integer lon_idx           ! [idx] Counting index
    integer lon_nbr           ! [nbr] Dimension size
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
    real(r8) lat(lat_nbr_max)     ! Entire lat and lon arrays are retrieved from disk
    real(r8) lon(lon_nbr_max)
    
#ifdef SCCM
    integer lat_idx_sccm      ! SCCM nearest latitude index into boundary data
    integer lon_idx_sccm      ! SCCM nearest longitude index into boundary data
    real(r8) dst_crr              ! Current least distance from gridpoint to target
    real(r8) dst_new              ! Distance of current gridpoint to target
#endif /* not SCCM */ 
    
    ! Initialize default values
    srt_lon_lat=(/1,1/)
    srt_lon_lat_bln=(/1,1,1/)
    cnt_lon_lat=(/plon,plat/)
    cnt_lon_lat_bln=(/plon,plat,bln_nbr/)
    
    ! Read in land surface types from high resolution grid
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=rcd+nf90_open(fl_in,nf90_nowrite,nc_id)
    if (rcd /= nf90_noerr) then
       print*,'Cannot read file '//trim(fl_in)
       call nf90_err_exit(rcd,"dst_tibds_ini()")
    endif
    write(6,"(a29,1x,a)") "dst: Opened netCDF input file",fl_in
    
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"bln",bln_dmn_id)
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr),"dst_tibds_ini(): inquire_dim lat")
    if (lat_nbr /= lat_nbr_max.and.lat_nbr /= 1) stop "lat_nbr /= lat_nbr_max.and.lat_nbr /= 1"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr),"dst_tibds_ini(): inquire_dim lon")
    if (lon_nbr /= lon_nbr_max.and.lon_nbr /= 1) stop "lon_nbr /= lon_nbr_max.and.lon_nbr /= 1"
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bln_dmn_id,len=bln_nbr_lcl),"dst_tibds_ini(): inquire_dim bln")
    if (bln_nbr_lcl /= bln_nbr) stop "bln_nbr_lcl /= bln_nbr"
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"lon",lon_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lnd_frc_dry",lnd_frc_dry_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mbl_bsn_fct",mbl_bsn_fct_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_CaCO3",mss_frc_CaCO3_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_snd",mss_frc_snd_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mss_frc_cly",mss_frc_cly_id)
    rcd=nf90_wrp_inq_varid(nc_id,"sfc_frc_bln",sfc_frc_bln_id)
    rcd=nf90_wrp_inq_varid(nc_id,"sfc_typ",sfc_typ_id)
    ! Get data
    rcd=rcd+nf90_get_var(nc_id,lat_id,lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var lat")
    rcd=rcd+nf90_get_var(nc_id,lon_id,lon)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_var lon")
    ! In SCCM, find nearest lat and lon indices into boundary data
#ifdef SCCM
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
    
    rcd=rcd+nf90_get_var(nc_id,lnd_frc_dry_id,lnd_frc_dry,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara lnd_frc_dry")
    rcd=rcd+nf90_get_var(nc_id,mbl_bsn_fct_id,mbl_bsn_fct,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mbl_bsn_fct")
    rcd=rcd+nf90_get_var(nc_id,mss_frc_CaCO3_id,mss_frc_CaCO3,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_CaCO3")
    rcd=rcd+nf90_get_var(nc_id,mss_frc_snd_id,mss_frc_snd,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_snd")
    rcd=rcd+nf90_get_var(nc_id,mss_frc_cly_id,mss_frc_cly,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara mss_frc_cly")
    rcd=rcd+nf90_get_var(nc_id,sfc_frc_bln_id,sfc_frc_bln,start=srt_lon_lat_bln,count=cnt_lon_lat_bln)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara sfc_frc_bln")
    rcd=rcd+nf90_get_var(nc_id,sfc_typ_id,sfc_typ,start=srt_lon_lat,count=cnt_lon_lat)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"dst_tibds_ini():get_vara sfc_typ")
    
    rcd=rcd+nf90_close(nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,"Unable to close file "//fl_in)
    write (6,"(a25,1x,a)") "dst: Ingested netCDF file",fl_in
    
    ! Sanity check (can save lots of wasted time if boundary dataset is corrupt)
    flg_err=.false.           ! [flg] Bad input data
    do lat_idx=1,lat_nbr
       do lon_idx=1,lon_nbr
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
             write (6,"(a,2(a4,i3,a4,f9.4,a2))")  &
                  "dst: ERROR dst_tibds_ini() sanity check failed at ", &
                  "lat(",lat_idx,") = ",lat(lat_idx),", ", &
                  "lon(",lon_idx,") = ",lon(lon_idx),", "
             stop 
          endif               ! endif err
       end do                 ! end loop over lon
    end do                    ! end loop over lat
    return
  end subroutine dst_tibds_ini                       ! end dst_tibds_ini()
  
end module dsttibds

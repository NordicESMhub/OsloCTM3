!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// CTM averages
!//=========================================================================
module averages
  !//-----------------------------------------------------------------------
  !// MODULE: averages
  !// DESCRIPTION: Routines to make tracer averages (mostly not used).
  !//
  !// Important
  !//   The main work on making averages is done in oc_diagnostics.f90.
  !//   Here are some of the basic routines for making simple averages,
  !//   and also some print-to-screen routines.
  !//
  !// Contains
  !//   subroutine AVG_WRT_NC4
  !//   subroutine AVG_ADD2
  !//   subroutine AVG_CLR2
  !//   subroutine AVG_ADD2_H2O
  !//   subroutine AVG_P1
  !//   subroutine AVG_WRT2_OLD
  !//-----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'averages.f90'
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine AVG_WRT_NC4
    !//---------------------------------------------------------------------
    !// Writes average for STT and XSTT, as netCDF4 files.
    !// Includes tracer info and grid info.
    !//
    !// This routine is a bit long, since it does all the netcdf
    !// calls necessary. Sorry about that.
    !//
    !// Amund Sovde Haslerud, May 2017
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, r4, i8
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR, IPARW, JPARW, &
         LPARW, IDBLK, JDBLK, MPBLK, TNAMELEN
    use cmn_ctm, only: JYEAR, JMON, JDATE, IDAY, &
         XDGRD, YDGRD, XDEDG, YDEDG, ETAA, ETAB, AREAXY, NTM, &
         ZGRD, ZEDG, ETAA, ETAB
    use cmn_met, only: METTYPE, metCYCLE, metREVNR, MET_ROOT, MPATH1,MPATH2
    use cmn_diag, only: NRAVG, STTAVG, AIRAVG, DVAVG, ZHAVG, PSFCAVG, &
         RUNTITLE, NDAY1, JYEAR1, JMON1, JDATE1, &
         metTypeInfo, resolutionInfo, nc4deflate_global, nc4shuffle_global
    use cmn_chem, only: TMASS, TNAME
    use cmn_parameters, only: M_AIR
    use cmn_oslo, only: chem_idx, Xchem_idx, XTMASS, trsp_idx, XTNAME, &
         TEMPAVG, H2OAVG, QAVG, AMAVG, LMTROPAVG, H2OAVG_LMT1, &
         RESULTDIR, XSTTAVG
    use netcdf
    use ncutils, only: handle_error
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer :: I, J, L, N, ioerr
    integer, dimension(6) :: start_time, end_time
    real(r8)  :: ZNRAVG
    character(len=80) ::  FN_AVG
    character(LEN=8) :: datestamp, datestamp0
    real(r4) :: r4xyz(IPAR,JPAR,LPAR) !// To put out reverse-order XSTTAVG
    !//---------------------------------------------------------------------
    !// Version number
    integer, parameter :: version = 8
    !// Update version 8: File format is netCDF4
    !//---------------------------------------------------------------------

    !// netCDF variables
    !//---------------------------------------------------------------------
    integer :: &
         status, &          !Status for netcdf file 0=OK
         ncid, &            !file id for output netcdf file
         lon_dim_id, &      !Dimension id for longitude
         lon_id, &          !Variable id for longitude
         lat_dim_id, &      !Dimension id for latitude
         lat_id, &          !Variable id for latitude
         lev_dim_id, &      !Dimension id for level
         lev_id, &          !Variable id for level
         ilon_dim_id, &      !Dimension id for longitude interstices
         ilon_id, &          !Variable id for longitude interstices
         ilat_dim_id, &      !Dimension id for latitude interstices
         ilat_id, &          !Variable id for latitude interstices
         ilev_dim_id, &     !Dimension id for level interstices
         ilev_id, &         !Variable id for level interstices
         ihya_dim_id, &     !Dimension id for hybrid A on interstices
         ihya_id, &         !Variable id for hybrid A on interstices
         ihyb_dim_id, &     !Dimension id for hybrid B on interstices
         ihyb_id, &         !Variable id for hybrid B on interstices
         ncomps_dim_id, &     !Dimension id for NPAR+NOTRPAR
         tracer_name_len_dim_id, & !Dimension id for tname charater length
         date_size_dim_id, & !Dimension id for date sizes
         dim_lon_lat_id(2), &
         dim_lon_lat_lev_id(3), &
         dim_ilev_lon_lat_id(3), &
         nravg_id, &                !ID for nNumber of accumulations in average
         version_id, &              !ID for file version number
         start_time_id, end_time_id, & !IDs for start/end dates for average
         start_day_id, end_day_id, &   !IDs for start/end day (NDAY)
         npar_id, notrpar_id, &        !IDs for scalar output NPAR/NOTRPAR
         tracer_idx_id, &     !ID for all tracer number
         tracer_molw_id, & !ID for all tracer molecular weights
         tracer_name_id, & !ID for all tracer names
         tracer_transported_id, & !ID for transport flag
         native_lon_id, &          !Variable id for native longitude size
         native_lat_id, &          !Variable id for native latitude size
         native_lev_id, &          !Variable id for native level size
         areaxy_id, &
         psfc_id, &
         air_id, volume_id, airdens_id, temperature_id, &
         height_id, lmtrop_id, h2o_id, q_id, &
         comps_id(NPAR+NOTRPAR)
    character(len=TNAMELEN), dimension(NPAR+NOTRPAR) :: tracer_name
    real(r8), dimension(NPAR+NOTRPAR) :: tracer_molw
    integer, dimension(NPAR+NOTRPAR) :: tracer_idx
    integer, dimension(NPAR+NOTRPAR) :: tracer_transported
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'AVG_WRT_NC4'
    !//---------------------------------------------------------------------

    ZNRAVG = 1._r8 / real(NRAVG,r8)
    write(6,'(A,I5,F9.5)') 'avg_wrt_nc ',NRAVG,ZNRAVG
    !// The transported species
    do N = 1, NTM
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                STTAVG(I,J,L,N) = ZNRAVG * STTAVG(I,J,L,N)
             end do
          end do
       end do
    end do
    !// H2O
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             H2OAVG(I,J,L) = ZNRAVG * H2OAVG(I,J,L)
          end do
       end do
    end do
    !// Q
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             QAVG(I,J,L) = ZNRAVG * QAVG(I,J,L)
          end do
       end do
    end do
    !// Air
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             AIRAVG(I,J,L) = ZNRAVG * AIRAVG(I,J,L)
          end do
       end do
    end do
    !// DV
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             DVAVG(I,J,L) = ZNRAVG * DVAVG(I,J,L)
          end do
       end do
    end do
    !// Air density
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             AMAVG(I,J,L) = ZNRAVG * AMAVG(I,J,L)
          end do
       end do
    end do
    !// Temperature
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             TEMPAVG(I,J,L) = ZNRAVG * TEMPAVG(I,J,L)
          end do
       end do
    end do
    !// Surface pressure
    do J = 1,JPAR
       do I = 1,IPAR
          PSFCAVG(I,J) = ZNRAVG * PSFCAVG(I,J)
       end do
    end do
    !// Height
    do J = 1,JPAR
       do I = 1,IPAR
          do L = 1,LPAR+1
             ZHAVG(L,I,J) = ZNRAVG * ZHAVG(L,I,J)
          end do
       end do
    end do
    !// LMTROP
    do J = 1,JPAR
       do I = 1,IPAR
          LMTROPAVG(I,J) = ZNRAVG * LMTROPAVG(I,J)
       end do
    end do
    !// H2O mixing ratio at LMTROP+1
    do J = 1,JPAR
       do I = 1,IPAR
          H2OAVG_LMT1(I,J) = ZNRAVG * H2OAVG_LMT1(I,J)
       end do
    end do

    !// Non-transported species if present
    if (NOTRPAR .gt. 0) then
       do J = 1,JPAR
          do I = 1,IPAR
             do N = 1,NOTRPAR
                do L = 1,LPAR
                   XSTTAVG(L,N,I,J) = ZNRAVG * XSTTAVG(L,N,I,J)
                end do
             end do
          end do
       end do
    end if



    !// Make list of all tracer names
    do N = 1, NPAR
       tracer_name(N) = TNAME(N)
    end do
    do N = 1, NOTRPAR
       tracer_name(NPAR + N) = XTNAME(N)
    end do

    !// List of component numbers
    do N = 1, NPAR
       tracer_idx(N) = chem_idx(N)
    end do
    do N = 1, NOTRPAR
       tracer_idx(NPAR + N) = Xchem_idx(N)
    end do

    !// List of flags for transported (1) or non-transported (0)
    do N = 1, NPAR
       tracer_transported(N) = 1
    end do
    do N = 1, NOTRPAR
       tracer_transported(NPAR + N) = 0
    end do


    !// Make list of their molecular mass
    do N = 1, NPAR
       tracer_molw(N) = TMASS(N)
    end do
    do N = 1, NOTRPAR
       tracer_molw(NPAR + N) = XTMASS(N)
    end do


    !//---------------------------------------------------------------------
    !// Write netcdf4 file.
    !//---------------------------------------------------------------------
    !// File name is avgsavSTARTDATE-ENDDATE.nc, where STARTDATE is taken
    !// from JYEAR1,JMON1,JDATE1,00UTC and ENDDATE is taken from current
    !// time, which is JYEAR,JMON,JDATE,00UTC (already updated by calendar).
    !// Format of the dates are YYYYMMDDHHMM, where HHMM always is 0000.
    !//---------------------------------------------------------------------
    !// Start of average period at 00:00UTC
    start_time = (/JYEAR1,JMON1,JDATE1,0,0,0/)
    write(datestamp0(1:8),'(i4.4,2i2.2)') JYEAR1,JMON1,JDATE1
    !// End of average period at 00:00UTC
    end_time = (/JYEAR,JMON,JDATE,0,0,0/)
    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE

    !//---------------------------------------------------------------------
    !// File name
    FN_AVG = trim(RESULTDIR)//'avgsav_'//datestamp0//'_'//datestamp//'.nc'
    !// Open netCDF4 file for writing
    status=nf90_create(path=FN_AVG,cmode=nf90_netcdf4,ncid=ncid)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title','Oslo CTM3 averages')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    status=nf90_put_att(ncid,nf90_global,'model_info', &
         'Oslo CTM3 is a 3D Chemical Transport Model')
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': model_info')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology',metTypeInfo)
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': driving_meteorology')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology_path',MPATH1)
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': driving_meteorology_path')
    status=nf90_put_att(ncid,nf90_global,'resolution_info',resolutionInfo)
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': resolution_info')

    status=nf90_put_att(ncid,nf90_global,'runtitle',trim(runtitle))
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': runtitle')
    status=nf90_put_att(ncid,nf90_global,'contact_info','For errors, contact CICERO')
    if (status/=nf90_noerr) call handle_error(status, f90file//':'//subr//': contact_info')

    !// Define lat/lon/lev
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')
    status = nf90_def_dim(ncid,"lev",LPAR,lev_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev dim')

    !// Define ilat/ilon/ilev
    status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':define ilon dim')
    status = nf90_def_dim(ncid,"ilev",LPAR+1,ilev_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev dim')

    !// Define NCOMPS = NPAR+NOTRPAR
    status = nf90_def_dim(ncid,"NCOMPS",NPAR+NOTRPAR,ncomps_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NCOMPS dim')

    !// Define length of tracer name string (use TNAME for this)
    status = nf90_def_dim(ncid,"tracer_name_len",TNAMELEN,tracer_name_len_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,"date_size",size(start_time),date_size_dim_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')

     
    !// Defining the combined id for a 2d field (lon,lat)
    dim_lon_lat_id(1) = lon_dim_id
    dim_lon_lat_id(2) = lat_dim_id
    !// Defining the combined id for a 3d field (lon,lat,lev)
    dim_lon_lat_lev_id(1) = lon_dim_id
    dim_lon_lat_lev_id(2) = lat_dim_id
    dim_lon_lat_lev_id(3) = lev_dim_id
    !// Defining the combined id for a reversed 3d field (ilev,lon,lat)
    dim_ilev_lon_lat_id(1) = ilev_dim_id
    dim_ilev_lon_lat_id(2) = lon_dim_id
    dim_ilev_lon_lat_id(3) = lat_dim_id


    !// Defining the lon variable
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    !// Attributes
    status = nf90_put_att(ncid,lon_id,'units','degrees east')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lon_id,'description','Value at grid box center.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lon')

    !// Defining the lat variable
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    !// Attributes
    status = nf90_put_att(ncid,lat_id,'units','degrees north')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lat_id,'description','Value at grid box center.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lat')

    !// Defining the lev variable
    status = nf90_def_var(ncid,"lev",nf90_double,lev_dim_id,lev_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev variable')
    !// Attributes
    status = nf90_put_att(ncid,lev_id,'units','hPa')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lev')
    status = nf90_put_att(ncid,lev_id,'description', &
         'Pressure at mass-weighted grid box center (mean of pressures at upper and lower edges), assuming surface at 1000hPa.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lev')


    !// Defining ilon variable (lon on interstices)
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    !// Attributes
    status = nf90_put_att(ncid,ilon_id,'units','degrees east')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilon_id,'description','Value at eastern edge of grid box.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilon')

    !// Defining ilat variable (lat on interstices)
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    !// Attributes
    status = nf90_put_att(ncid,ilat_id,'units','degrees north')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilat_id,'description','Value at southern edge of grid box.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilat')

    !// Defining ilev variable (lev on interstices)
    status = nf90_def_var(ncid,"ilev",nf90_double,ilev_dim_id,ilev_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev variable')
    !// Attributes
    status = nf90_put_att(ncid,ilev_id,'units','hPa')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilev')
    status = nf90_put_att(ncid,ilev_id,'description', &
         'Pressure at lower edge of grid box, assuming surface at 1000hPa.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilev')

    !// Defining hybrid sigma coordinate A
    status = nf90_def_var(ncid,"ihya",nf90_double,ilev_dim_id,ihya_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihya variable')
    !// Attributes
    status = nf90_put_att(ncid,ihya_id,'units','hPa')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihya')
    status = nf90_put_att(ncid,ihya_id,'description', 'Sigma hybrid coordinate A.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihya')
    status = nf90_put_att(ncid,ihya_id,'usage','p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihya')

    !// Defining hybrid sigma coordinate B
    status = nf90_def_var(ncid,"ihyb",nf90_double,ilev_dim_id,ihyb_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihyb variable')
    !// Attributes
    status = nf90_put_att(ncid,ihyb_id,'units','1')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihyb')
    status = nf90_put_att(ncid,ihyb_id,'description', 'Sigma hybrid coordinate B.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihyb')
    status = nf90_put_att(ncid,ihyb_id,'usage','p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihyb')

    !// Info about native size IPARW
    status = nf90_def_var(ncid,"IPARW",nf90_int,native_lon_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPARW variable')
    status = nf90_put_att(ncid,native_lon_id,'description', &
         'Meteorological data native longitudinal resolution')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPARW')
    !// Info about native size JPARW
    status = nf90_def_var(ncid,"JPARW",nf90_int,native_lat_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPARW variable')
    status = nf90_put_att(ncid,native_lat_id,'description', &
         'Meteorological data native latitudinal resolution')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPARW')
    !// Info about native size LPARW
    status = nf90_def_var(ncid,"LPARW",nf90_int,native_lev_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPARW variable')
    status = nf90_put_att(ncid,native_lev_id,'description', &
         'Meteorological data vertical resolution')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPARW')

    !// Number of time steps accumulated - NRAVG
    status = nf90_def_var(ncid,"NRAVG",nf90_int,nravg_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NRAVG variable')
    status = nf90_put_att(ncid,nravg_id,'description', &
         'Number of time steps accumulated')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description NRAVG')

    !// Start time for accumulating data - START_TIME
    status = nf90_def_var(ncid,"START_TIME",nf90_int,date_size_dim_id,start_time_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_TIME variable')
    status = nf90_put_att(ncid,start_time_id,'description', &
         'Start date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_TIME')

    !// End time for accumulating data - END_TIME
    status = nf90_def_var(ncid,"END_TIME",nf90_int,date_size_dim_id,end_time_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_TIME variable')
    status = nf90_put_att(ncid,end_time_id,'description', &
         'End date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_TIME')

    !// Start day for accumulating data - START_DAY (NDAY)
    status = nf90_def_var(ncid,"START_DAY",nf90_int,start_day_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_DAY variable')
    status = nf90_put_att(ncid,start_day_id,'description', &
         'Start day for accumulating data (model counter NDAY).')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_DAY')

    !// End day for accumulating data - END_DAY
    status = nf90_def_var(ncid,"END_DAY",nf90_int,end_day_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_DAY variable')
    status = nf90_put_att(ncid,end_day_id,'description', &
         'End day for accumulating data (model counter NDAY).')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_DAY')

    !// File version number
    status = nf90_def_var(ncid,"VERSION",nf90_int,version_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define VERSION variable')
    status = nf90_put_att(ncid,version_id,'description', &
         'Output file version number')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description VERSION')

    !// NPAR and NOTRPAR - they could be handy
    status = nf90_def_var(ncid,"NPAR",nf90_int,npar_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NPAR variable')
    status = nf90_put_att(ncid,npar_id,'description', &
         'Number of transported species in simulation.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description NPAR')
    status = nf90_def_var(ncid,"NOTRPAR",nf90_int,notrpar_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NOTRPAR variable')
    status = nf90_put_att(ncid,notrpar_id,'description', &
         'Number of non-transported species in simulation.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description NOTRPAR')


    !// All the components
    status = nf90_def_var(ncid,"tracer_idx",nf90_int,ncomps_dim_id,tracer_idx_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_idx variable')
    status = nf90_put_att(ncid,tracer_idx_id,'description', &
         'ID numbers for species used in Oslo CTM3.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_idx')

    !// Molecular masses of transported components (r8)
    status = nf90_def_var(ncid,"tracer_molweight",nf90_double,ncomps_dim_id,tracer_molw_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_molweight variable')
    status = nf90_put_att(ncid,tracer_molw_id,'units','g/mol')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_molweight')
    status = nf90_put_att(ncid,tracer_molw_id,'description', &
         'Molecular weights.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_molweight')

    !// Tracer names
    status = nf90_def_var(ncid,"tracer_name",nf90_char, &
         (/tracer_name_len_dim_id,ncomps_dim_id/),tracer_name_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Tracer/species names.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Transported or not?
    status = nf90_def_var(ncid,"tracer_transported",nf90_int,ncomps_dim_id,tracer_transported_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_transported variable')
    status = nf90_put_att(ncid,tracer_transported_id,'units','1/0')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_transported')
    status = nf90_put_att(ncid,tracer_transported_id,'description', &
         'Flag if tracer is transported or not.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_transported')



    !// Grid area (r8), deflate netcdf4
    status = nf90_def_var(ncid,"gridarea",nf90_double,dim_lon_lat_id,areaxy_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'units','m2')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')

    !// Average surface pressure (r4), deflate netcdf4
    status = nf90_def_var(ncid,"Psfc",nf90_float,dim_lon_lat_id,psfc_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define Psfc variable')
    status = nf90_def_var_deflate(ncid,psfc_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate Psfc variable')
    status = nf90_put_att(ncid,psfc_id,'units','hPa')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit Psfc')


    !// Air (r4), deflate netcdf4
    status = nf90_def_var(ncid,"AIR",nf90_float,dim_lon_lat_lev_id,air_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define AIR variable')
    status = nf90_def_var_deflate(ncid,air_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate AIR variable')
    status = nf90_put_att(ncid,air_id,'units','kg')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit AIR')

    !// Volume (r4), deflate netcdf4
    status = nf90_def_var(ncid,"volume",nf90_float,dim_lon_lat_lev_id,volume_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define volume variable')
    status = nf90_def_var_deflate(ncid,volume_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate volume variable')
    status = nf90_put_att(ncid,volume_id,'units','m3')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit volume')

    !// Air density (r4), deflate netcdf4
    status = nf90_def_var(ncid,"air_density",nf90_float,dim_lon_lat_lev_id,airdens_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define air_density variable')
    status = nf90_def_var_deflate(ncid,airdens_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate air_density variable')
    status = nf90_put_att(ncid,airdens_id,'units','molec/cm3')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit air_density')

    !// Tempterature (r4), deflate netcdf4
    status = nf90_def_var(ncid,"temperature",nf90_float,dim_lon_lat_lev_id,temperature_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define temperature variable')
    status = nf90_def_var_deflate(ncid,temperature_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate temperature variable')
    status = nf90_put_att(ncid,temperature_id,'units','K')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit temperature')

    !// Height (r4), deflate netcdf4
    status = nf90_def_var(ncid,"height",nf90_float,dim_ilev_lon_lat_id,height_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define height variable')
    status = nf90_def_var_deflate(ncid,height_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate height variable')
    status = nf90_put_att(ncid,height_id,'units','m')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units height')
    status = nf90_put_att(ncid,height_id,'info','Reverse order (ilev,lon,lat)')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute info height')

    !// LMTROP, deflate netcdf4
    status = nf90_def_var(ncid,"LMTROP",nf90_float,dim_lon_lat_id,lmtrop_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LMTROP variable')
    status = nf90_def_var_deflate(ncid,lmtrop_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate LMTROP variable')
    status = nf90_put_att(ncid,lmtrop_id,'units','Uppermost model level in troposphere.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit LMTROP')

    !// H2O (r4), deflate netcdf4
    status = nf90_def_var(ncid,"H2O_all",nf90_float,dim_lon_lat_lev_id,h2o_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define H2O_all variable')
    status = nf90_def_var_deflate(ncid,h2o_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate H2O_all variable')
    status = nf90_put_att(ncid,h2o_id,'units','kg')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units H2O_all')
    status = nf90_put_att(ncid,h2o_id,'info','H2O from metdata Q, but overwritten in the stratosphere by calculated H2O when available.')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units H2O_all')

    !// Q (r4), deflate netcdf4
    status = nf90_def_var(ncid,"Q",nf90_float,dim_lon_lat_lev_id,q_id)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define Q variable')
    status = nf90_def_var_deflate(ncid,q_id,nc4shuffle,1,nc4deflate)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate Q variable')
    status = nf90_put_att(ncid,q_id,'units','kg(H2O)/kg(air)')
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units Q')
    status = nf90_put_att(ncid,q_id,'description', &
            'Specific humidity retrieved from meteorological data.')
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description QAVG')

    !// Component by component (r4)
    do N = 1, NPAR+NOTRPAR
       status = nf90_def_var(ncid,trim(tracer_name(N)),nf90_float,dim_lon_lat_lev_id,comps_id(N))
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define '//trim(tracer_name(N))//' variable')
       status = nf90_def_var_deflate(ncid,comps_id(N),nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate '//trim(tracer_name(N))//' variable')
       status = nf90_put_att(ncid,comps_id(N),'units','kg of species')
       if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units '//trim(tracer_name(N)))
    end do



    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------
    
 
    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')
    status = nf90_put_var(ncid,lev_id,ZGRD)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lev')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')
    status = nf90_put_var(ncid,ilev_id,ZEDG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilev')

    status = nf90_put_var(ncid,ihya_id,ETAA)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihya')
    status = nf90_put_var(ncid,ihyb_id,ETAB)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihyb')

    !// Info about native size
    status = nf90_put_var(ncid,native_lon_id,IPARW)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPARW')
    status = nf90_put_var(ncid,native_lat_id,JPARW)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPARW')
    status = nf90_put_var(ncid,native_lev_id,LPARW)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPARW')


    !// Number of time steps added
    status = nf90_put_var(ncid,NRAVG_id,NRAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NRAVG')
    !// Start time
    status = nf90_put_var(ncid,start_time_id,start_time)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_TIME')
    !// End time
    status = nf90_put_var(ncid,end_time_id,end_time)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_TIME')
    !// Start day (NDAY counter)
    status = nf90_put_var(ncid,start_day_id,NDAY1)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_DAY')
    !// End day (NDAY counter)
    status = nf90_put_var(ncid,end_day_id,IDAY-1)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_DAY')

    !// File version number
    status = nf90_put_var(ncid,version_id,version)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting VERSION')

    !// NPAR and NOTRPAR
    status = nf90_put_var(ncid,npar_id,NPAR)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NPAR')
    status = nf90_put_var(ncid,notrpar_id,NOTRPAR)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NOTRPAR')

    !// All the transported components
    status = nf90_put_var(ncid,tracer_idx_id,tracer_idx)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':putting tracer_idx')

    !// Molecular masses of transported components (r8)
    status = nf90_put_var(ncid,tracer_molw_id,tracer_molw)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_molw')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,tracer_name)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Tracer transport flag
    status = nf90_put_var(ncid,tracer_transported_id,tracer_transported)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_transported')


    !// Grid area (r8)
    status = nf90_put_var(ncid,areaxy_id,AREAXY)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AREAXY')

    !// Average surface pressure (r4)
    status = nf90_put_var(ncid,psfc_id,PSFCAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting PSFCAVG')

    !// Average gridbox air mass (r4)
    status = nf90_put_var(ncid,air_id,AIRAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AIRAVG')

    !// Average gridbox volume (r4)
    status = nf90_put_var(ncid,volume_id,DVAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting DVAVG')

    !// Average air density (r4)
    status = nf90_put_var(ncid,airdens_id,AMAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AMAVG')

    !// Average temperature (r4)
    status = nf90_put_var(ncid,temperature_id,TEMPAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting TEMPAVG')

    !// Average heights (r4)
    status = nf90_put_var(ncid,height_id,ZHAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ZHAVG')

    !// Average LMTROP (r4)
    status = nf90_put_var(ncid,lmtrop_id,LMTROPAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LMTROPAVG')

    !// Average H2O (trop+strat) (r4)
    status = nf90_put_var(ncid,h2o_id,H2OAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting H2OAVG')

    !// Average Q (r4)
    status = nf90_put_var(ncid,q_id,QAVG)
    if (status/=nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting QAVG')

    !// Component by component (r4) - transported first
    do N = 1, NPAR
       status = nf90_put_var(ncid,comps_id(N),STTAVG(:,:,:,N))
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting '//trim(TNAME(N)))
    end do

    !// Component by component - not-transported (r4)
    if (NOTRPAR .gt. 0) then
       do N = 1, NOTRPAR
          !// Back from reversed order
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   r4xyz(I,J,L) = XSTTAVG(L,N,I,J)
                end do
             end do
          end do
          status = nf90_put_var(ncid,comps_id(NPAR+N),r4xyz)
          if (status/=nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting '//trim(XTNAME(N)))
       end do
    end if !// if (NOTRPAR .gt. 0) then


    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status/=nf90_noerr) call handle_error(status,'close file: '//trim(fn_avg))

    write(6,'(a,2i8)') &
         f90file//':'//subr//': averages: '//trim(FN_AVG),NRAVG,IDAY-1
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !//---------------------------------------------------------------------
  end subroutine AVG_WRT_NC4
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine AVG_ADD2_H2O(BTT,AIRB,MP)
    !//---------------------------------------------------------------------
    !// Routine for H2O 3D averages, when H2O is not calculated in
    !// chemistry. Works in IJ-blocks.
    !//
    !// Ole Amund Sovde, March 2010
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, &
         MPBLK, LOSLOCSTRAT
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TMASSMIX2MOLMIX
    use cmn_met, only: Q
    use cmn_parameters, only: M_AIR
    use cmn_oslo, only: H2OAVG, LMTROP, trsp_idx, QAVG, H2OAVG_LMT1
    use strat_h2o, only: STR_H2O
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: AIRB
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT
    integer, intent(in) :: MP

    !// Locals
    integer :: I, J, L, II, JJ, LL
    !//---------------------------------------------------------------------

    !// Q from metdata
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             QAVG(I,J,L) = QAVG(I,J,L) + Q(I,J,L)
          end do
       end do
    end do

    !// Tropospheric H2O from metdata (stratosphere fetched below)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = 1, LMTROP(I,J)
             H2OAVG(I,J,L) = H2OAVG(I,J,L) + Q(I,J,L)*AIRB(L,II,JJ)
          end do
       end do
    end do

    !// Stratospheric H2O from BTT or STR_H2O
    if (trsp_idx(114).gt.0) then
       !// Aircraft H2O is included in tracer 114
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = LMTROP(I,J)+1, LPAR
                H2OAVG(I,J,L) = H2OAVG(I,J,L) + BTT(L,trsp_idx(114),II,JJ)
             end do
             !// Lower stratosphere H2O mixing ratio
             LL = LMTROP(I,J)+1
             H2OAVG_LMT1(I,J) = H2OAVG_LMT1(I,J) &
                  + BTT(LL,trsp_idx(114),II,JJ) &
                    / AIRB(LL,II,JJ) * TMASSMIX2MOLMIX(trsp_idx(114))
          end do
       end do
    else
       !// Aircraft H2O must be added separately
       if (LOSLOCSTRAT) then
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = LMTROP(I,J)+1, LPAR
                   H2OAVG(I,J,L) = H2OAVG(I,J,L) + STR_H2O(L,II,JJ,MP)
                end do
                !// Lower stratosphere H2O mixing ratio
                LL = LMTROP(I,J)+1
                H2OAVG_LMT1(I,J) = H2OAVG_LMT1(I,J) &
                     + STR_H2O(LL,II,JJ,MP) / AIRB(LL,II,JJ) * M_AIR/18._r8
             end do
          end do
       else
          !// Tropospheric chemistry only
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = LMTROP(I,J)+1, LPAR
                   H2OAVG(I,J,L) = H2OAVG(I,J,L) + Q(I,J,L)*AIRB(L,II,JJ)
                end do
                !// Lower stratosphere H2O mixing ratio; fetch from Q
                LL = LMTROP(I,J)+1
                H2OAVG_LMT1(I,J) = H2OAVG_LMT1(I,J) &
                     + Q(I,J,LL)*M_AIR/18._r8
             end do
          end do
       end if

    end if !// if (trsp_idx(114).gt.0) then

    !//---------------------------------------------------------------------
  end subroutine AVG_ADD2_H2O
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine AVG_ADD2
    !//---------------------------------------------------------------------
    !// Improved adding routine for 3D averages. Treats STT, XSTT and PSFC.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rAvg
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: NTM, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, STT, AIR
    use cmn_diag, only: STTAVG, AIRAVG, DVAVG, ZHAVG, PSFCAVG, NRAVG
    use cmn_met, only: T, P, ZOFLE
    use cmn_oslo, only: TEMPAVG, AMAVG, DV_IJ, LMTROP, &
          LMTROPAVG, XSTT, XSTTAVG
    use cmn_parameters, only: AVOGNR, M_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: I, J, L, N, II, JJ, MP
    !// --------------------------------------------------------------------

    NRAVG = NRAVG + 1
    !// The transported species
    do N = 1,NTM
       do L = 1,LPAR
          do J = 1,JPAR
             do I = 1,IPAR
                STTAVG(I,J,L,N) = STTAVG(I,J,L,N) + real(STT(I,J,L,N),rAvg)
             end do
          end do
       end do
    end do
    !// Air
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             AIRAVG(I,J,L) = AIRAVG(I,J,L) + real(AIR(I,J,L),rAvg)
          end do
       end do
    end do
    !// Volume: DV_IJ(LPAR,IDBLK,JDBLK,MPBL) calculated in oc_main.f
    do L = 1,LPAR
       do MP = 1,MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                DVAVG(I,J,L) = DVAVG(I,J,L) + DV_IJ(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    !// Air density: DV_IJ(LPAR,IDBLK,JDBLK,MPBL) calculated in oc_main.f
    do L = 1, LPAR
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                AMAVG(I,J,L) = AMAVG(I,J,L) &
                     + real(AIR(I,J,L)*1.e-3_r8 &
                           /(M_AIR*DV_IJ(L,II,JJ,MP)) * AVOGNR, rAvg)
             end do
          end do
       end do
    end do
    !// Temperature
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             TEMPAVG(I,J,L) = TEMPAVG(I,J,L) + real(T(I,J,L),rAvg)
          end do
       end do
    end do
    !// Surface pressure
    do J = 1, JPAR
       do I = 1, IPAR
          PSFCAVG(I,J) = PSFCAVG(I,J) + real(P(I,J),rAvg)
       end do
    end do
    !// ZOFLE (change to I,J,L)
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LPAR+1
             ZHAVG(L,I,J) = ZHAVG(L,I,J) + real(ZOFLE(L,I,J), rAvg)
          end do
       end do
    end do
    !// LMTROP
    do J = 1, JPAR
       do I = 1, IPAR
          LMTROPAVG(I,J) = LMTROPAVG(I,J) + real(LMTROP(I,J), rAvg)
       end do
    end do

    !// The non-transported species
    if (NOTRPAR.gt.0) then
       do J = 1, JPAR
          do I = 1, IPAR
             do N = 1, NOTRPAR
                do L = 1, LPAR
                  XSTTAVG(L,N,I,J) = XSTTAVG(L,N,I,J) + real(XSTT(L,N,I,J),rAvg)
                end do
             end do
          end do
       end do
    end if

    !// --------------------------------------------------------------------
  end subroutine AVG_ADD2
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine AVG_CLR2
    !// --------------------------------------------------------------------
    !// Improved routine for clearing averages.
    !//
    !// Ole Amund Sovde, October 2013
    !//    Added QAVG and H2OAVG_LMT1.
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rAvg
    use cmn_size, only: NOTRPAR
    use cmn_ctm, only: NTM, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, GMTAU, &
         IDAY, JDAY, JYEAR, JMON, JDATE, TMON
    use cmn_diag, only: STTAVG, AIRAVG, DVAVG, ZHAVG, PSFCAVG, NRAVG, &
         TAU1, NDAY1, JDAY1, JYEAR1, JMON1, JDATE1, TMON1
    use cmn_oslo, only: TEMPAVG, H2OAVG, AMAVG, LMTROPAVG, &
         QAVG, H2OAVG_LMT1, XSTTAVG
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    NRAVG = 0
    !// The transported species
    STTAVG(:,:,:,:) = 0._rAvg
    !// Air
    AIRAVG(:,:,:)   = 0._rAvg
    !// Volume
    DVAVG(:,:,:)    = 0._rAvg
    !// Surface pressure
    PSFCAVG(:,:)    = 0._rAvg
    !// Temperature
    TEMPAVG(:,:,:)  = 0._rAvg
    !// H2O
    H2OAVG(:,:,:)   = 0._rAvg
    !// Q
    QAVG(:,:,:)     = 0._rAvg
    !// Air density
    AMAVG(:,:,:)    = 0._rAvg
    !// Height
    ZHAVG(:,:,:)    = 0._rAvg
    !// LMTROP
    LMTROPAVG(:,:)  = 0._rAvg
    !// H2O mixing ratio at LMTROP+1
    H2OAVG_LMT1(:,:)= 0._rAvg

    !// The non-transported species
    if (NOTRPAR.gt.0) then
       XSTTAVG(:,:,:,:) = 0._r8
    end if


    !// Set time of AVERAGE initialization
    TAU1  = GMTAU
    NDAY1 = IDAY
    JDAY1 = JDAY
    JYEAR1= JYEAR
    JMON1 = JMON
    JDATE1 = JDATE
    TMON1 = TMON

    !// --------------------------------------------------------------------
  end subroutine AVG_CLR2
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine AVG_P1
    !//---------------------------------------------------------------------
    !// Processes/prints 0-D & 1-D average tracer mixing ratios for
    !// stdout!
    !// This uses the TENDENCY boxes but the calendar for AVGS
    !// assumes that STTAVG and AIRAVG are non-zero! and correctly
    !// summed/averaged by a call to AVG_WRT immediately before.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: IDAY, GMTAU,JDATE,TMON,JYEAR, NTM, &
         TALT, TLAT, TLNG
    use cmn_chem, only: TNAME
    use cmn_diag, only: STTAVG, AIRAVG, NTBPAR, NDAY1, &
         TAU1,JDATE1,TMON1,JYEAR1, LBGTA, &
         NBOXD, IBOXD, JBOXD, KBOXD, TBOXD
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    real(r8) :: DT_HR
    integer :: I, II, J, K, K2, L, M, N, IB0, IB1, JB0, JB1, LB0, LB1
    character(len=80) :: BTITLE
    real(r8) :: BGTBX(NTBPAR)
    real(r8) :: BGTBXL(LPAR,NTBPAR),BGTBXJ(JPAR,NTBPAR),BGTBXI(IPAR,NTBPAR)
    real(r8) :: BGTAX(NTBPAR)
    real(r8) :: BGTAXL(LPAR,NTBPAR),BGTAXJ(JPAR,NTBPAR),BGTAXI(IPAR,NTBPAR)
    !//---------------------------------------------------------------------
    
    DT_HR = 24._r8 * real(IDAY-NDAY1,r8)
    write(BTITLE(1:25),'(a25)') '1-D Average Profiles(m/m)'
    write(BTITLE(26:45),'(f8.3,i3,1x,a3,i5)') TAU1,JDATE1,TMON1,JYEAR1
    write(BTITLE(46:50),'(a5)') ' -to-'
    write(BTITLE(51:70),'(f8.3,i3,1x,a3,i5)') GMTAU,JDATE,TMON,JYEAR
    write(BTITLE(71:80),'(f10.2)') DT_HR

    !//---calculate 0-D and 1-D mixing ratio profiles for tendency boxes
    do N=1,NTM
      if (LBGTA(N)) then
        write(6,'(a10,a80)') TNAME(N), BTITLE
        BGTBX(:)    = 0.0_r8
        BGTBXL(:,:) = 0.0_r8
        BGTBXJ(:,:) = 0.0_r8
        BGTBXI(:,:) = 0.0_r8
        BGTAX(:)    = 0.0_r8
        BGTAXL(:,:) = 0.0_r8
        BGTAXJ(:,:) = 0.0_r8
        BGTAXI(:,:) = 0.0_r8

        do K = 1,NBOXD
          IB0 = IBOXD(1,K)
          IB1 = IBOXD(2,K)
          JB0 = JBOXD(1,K)
          JB1 = JBOXD(2,K)
          LB0 = KBOXD(1,K)
          LB1 = KBOXD(2,K)

          !// ---tracer average mixing ratio (0-D & 1-D) in box
          do L = LB0,LB1
            do J = JB0,JB1
              do II = IB0,IB1
                I = mod(II-1,IPAR) + 1
                BGTBX(   K) = BGTBX(   K) + STTAVG(I,J,L,N)
                BGTBXL(L,K) = BGTBXL(L,K) + STTAVG(I,J,L,N)
                BGTBXJ(J,K) = BGTBXJ(J,K) + STTAVG(I,J,L,N)
                BGTBXI(I,K) = BGTBXI(I,K) + STTAVG(I,J,L,N)
                BGTAX(   K) = BGTAX(   K) + AIRAVG(I,J,L)
                BGTAXL(L,K) = BGTAXL(L,K) + AIRAVG(I,J,L)
                BGTAXJ(J,K) = BGTAXJ(J,K) + AIRAVG(I,J,L)
                BGTAXI(I,K) = BGTAXI(I,K) + AIRAVG(I,J,L)
              end do
            end do
          end do

          BGTBX(K) = BGTBX(K)/BGTAX(K)
          do L = LB0,LB1
            BGTBXL(L,K) = BGTBXL(L,K)/BGTAXL(L,K)
          end do
          do J = JB0,JB1
            BGTBXJ(J,K) = BGTBXJ(J,K)/BGTAXJ(J,K)
          end do
          do II = IB0,IB1
            I = mod(II-1,IPAR) + 1
            BGTBXI(I,K) = BGTBXI(I,K)/BGTAXI(I,K)
          end do
        end do

        !//---print 0-D & 1-D mixing ratio profiles for each box K
        do K = 1,NBOXD,10
          K2 = min(K+9, NBOXD)
          write(6,'(2i4,10(1x,a9))') K,K2, (TBOXD(M),M=K,K2)
          write(6,'(a8, 1p,10e10.3)')'    0-D ',(BGTBX(M),M=K,K2)

          write(6,'(a10,a80)') TNAME(N), BTITLE
          do L = LPAR, 1, -1
            write(6,'(2a4,1p,10e10.3)')' L  ',TALT(L),(BGTBXL(L,M),M=K,K2)
          end do
          write(6,'(a10,a80)') TNAME(N), BTITLE
          do J = JPAR, 1, -1
            write(6,'(2a4,1p,10e10.3)')' J  ',TLAT(J),(BGTBXJ(J,M),M=K,K2)
          end do
          write(6,'(a10,a80)') TNAME(N), BTITLE
          do I = 1, IPAR
            write(6,'(2a4,1p,10e10.3)')' I  ',TLNG(I),(BGTBXI(I,M),M=K,K2)
          end do
        end do

      end if
    end do
        
    !//---------------------------------------------------------------------
  end subroutine AVG_P1
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine AVG_WRT2_OLD
    !//---------------------------------------------------------------------
    !// Writes average for STT and XSTT.
    !// Includes tracer info and grid info.
    !//
    !// Ole Amund Sovde, October 2013, ..., March 2010
    !//    See version history at the end of variable declarations.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: JYEAR, JMON, JDATE, IDAY, &
         XDGRD, YDGRD, XDEDG, YDEDG, ETAA, ETAB, AREAXY, NTM
    use cmn_diag, only: NRAVG, STTAVG, AIRAVG, DVAVG, ZHAVG, PSFCAVG, &
         RUNTITLE, NDAY1
    use cmn_chem, only: TMASS, TNAME
    use cmn_parameters, only: M_AIR
    use cmn_oslo, only: chem_idx, Xchem_idx, XTMASS, trsp_idx, XTNAME, &
         TEMPAVG, H2OAVG, QAVG, AMAVG, LMTROPAVG, H2OAVG_LMT1, &
         RESULTDIR, XSTTAVG
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer :: I, J, L, N,ioerr
    real(r8)  :: ZNRAVG
    character(len=80) ::  FN_AVG
    character(LEN=8) :: datestamp
    logical :: fnr_ok
    integer :: ifnr

    integer :: put_h2o !// flag to put out extra H2O and Q averages
    !//---------------------------------------------------------------------
    !// Version number
    integer, parameter :: version=7
    !// Update version 7: Xchem_idx, XTMASS and XTNAME are now put out before
    !//                   STT. ADDH2O is no longer in use; H2OAVG and QAVG
    !//                   are now put out as standard fields, along with
    !//                   a flag to possibly turn them off later.
    !// Extras in version 6: LMTROP
    !// Extras in version 5: ZOFLE
    !// Extras in version 4: DV
    !// Extras in version 3: TNAME, XTNAME
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'AVG_WRT2'
    !//---------------------------------------------------------------------

    ZNRAVG = 1._r8 / real(NRAVG,r8)
    write(6,'(A,I5,F9.5)') 'avg_wrt ',NRAVG,ZNRAVG
    !// The transported species
    do N = 1, NTM
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                STTAVG(I,J,L,N) = ZNRAVG * STTAVG(I,J,L,N)
             end do
          end do
       end do
    end do
    !// H2O
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             H2OAVG(I,J,L) = ZNRAVG * H2OAVG(I,J,L)
          end do
       end do
    end do
    !// Q
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             QAVG(I,J,L) = ZNRAVG * QAVG(I,J,L)
          end do
       end do
    end do
    !// Air
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             AIRAVG(I,J,L) = ZNRAVG * AIRAVG(I,J,L)
          end do
       end do
    end do
    !// DV
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             DVAVG(I,J,L) = ZNRAVG * DVAVG(I,J,L)
          end do
       end do
    end do
    !// Air density
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             AMAVG(I,J,L) = ZNRAVG * AMAVG(I,J,L)
          end do
       end do
    end do
    !// Temperature
    do L = 1,LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             TEMPAVG(I,J,L) = ZNRAVG * TEMPAVG(I,J,L)
          end do
       end do
    end do
    !// Surface pressure
    do J = 1,JPAR
       do I = 1,IPAR
          PSFCAVG(I,J) = ZNRAVG * PSFCAVG(I,J)
       end do
    end do
    !// Height
    do J = 1,JPAR
       do I = 1,IPAR
          do L = 1,LPAR+1
             ZHAVG(L,I,J) = ZNRAVG * ZHAVG(L,I,J)
          end do
       end do
    end do
    !// LMTROP
    do J = 1,JPAR
       do I = 1,IPAR
          LMTROPAVG(I,J) = ZNRAVG * LMTROPAVG(I,J)
       end do
    end do
    !// H2O mixing ratio at LMTROP+1
    do J = 1,JPAR
       do I = 1,IPAR
          H2OAVG_LMT1(I,J) = ZNRAVG * H2OAVG_LMT1(I,J)
       end do
    end do

    !// Non-transported species if present
    if (NOTRPAR .gt. 0) then
       do J = 1,JPAR
          do I = 1,IPAR
             do N = 1,NOTRPAR
                do L = 1,LPAR
                   XSTTAVG(L,N,I,J) = ZNRAVG * XSTTAVG(L,N,I,J)
                end do
             end do
          end do
       end do
    end if

    !// YYYYMMDD for filenames, taken from JYEAR,JMON,JDATE
    !// They are already updated by calendar.
    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
    FN_AVG = trim(RESULTDIR)//'avgsav'//datestamp//'.dta'

    !// Find file number
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do
    !// Open file
    open(ifnr,file=FN_AVG,form='UNFORMATTED',status='new',iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Error opening save-average file: '//trim(FN_AVG)
       stop 'STOP in '//subr
    end if
    !// Header
    write(ifnr) RUNTITLE
    !// Resolution (add one for additional H2O?)
    write(ifnr) IPAR,JPAR,LPAR,NPAR
    !// Average info
    write(ifnr) NRAVG,NDAY1,IDAY - 1,version

    !// All the transported components
    write(ifnr) chem_idx
    !// Molecular masses of transported components (r8)
    write(ifnr) TMASS
    !// Tracer names
    write(ifnr) TNAME

    !// Non-transported components
    write(ifnr) NOTRPAR
    if (NOTRPAR .gt. 0) then
       !// All the non-transported components
       write(ifnr) Xchem_idx
       !// Molecular masses of non-transported components
       write(ifnr) XTMASS
       !// Tracer names
       write(ifnr) XTNAME
    end if

    !// Grid area (r8)
    write(ifnr) AREAXY
    !// Longitude centers and edges (r8)
    write(ifnr) XDGRD,XDEDG
    !// Latitude centers and edges (r8)
    write(ifnr) YDGRD,YDEDG
    !// Hybrid sigma/pressure coordinates (r8)
    write(ifnr) ETAA,ETAB
    !// Average surface pressure (rAvg)
    write(ifnr) PSFCAVG

    !// Air (rAvg)
    write(ifnr) AIRAVG
    !// Volume (rAvg)
    write(ifnr) DVAVG
    !// Air density (rAvg)
    write(ifnr) AMAVG
    !// Tempterature (rAvg)
    write(ifnr) TEMPAVG
    !// Height (rAvg)
    write(ifnr) ZHAVG
    !// LMTROP
    write(ifnr) LMTROPAVG
    !// H2O (rAvg) and Q (rAvg) if wanted
    put_h2o = 1
    write(ifnr) put_h2o      !// Flag default output; could use LOSLOCTROP
    if (put_h2o .eq. 1) then !// Test on the same flag; could use LOSLOCTROP
       write(ifnr) H2OAVG
       write(ifnr) QAVG
    end if

    !// Component by component (rAvg)
    do N = 1, NPAR
       write(ifnr) chem_idx(N)
       write(ifnr) STTAVG(:,:,:,N)
    end do

    !// Write non-transported species if present
    if (NOTRPAR .gt. 0) then
       !// The tracer masses
       write(ifnr) XSTTAVG
    end if

    close(ifnr)
    write(6,'(a,2i8)') &
         'wrote STT+AIR+XSTT avgs: '//trim(FN_AVG),NRAVG,IDAY-1

    !//---------------------------------------------------------------------
  end subroutine AVG_WRT2_OLD
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
end module averages
!//=========================================================================

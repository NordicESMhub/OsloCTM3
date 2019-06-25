!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// CTM diagnostics
!//=========================================================================
module diagnostics_general
  !//-----------------------------------------------------------------------
  !// MODULE: diagnostics_general
  !// DESCRIPTION: Routines for diagnostics, developed along with the Oslo
  !//              chemistry routines.
  !// ----------------------------------------------------------------------
  !// Diagnostics for the Oslo CTM3.
  !// Contains:
  !//   subroutine diag_ohch4n2o_init
  !//   subroutine init_lifetime
  !//   subroutine du_columns
  !//   subroutine write_ducolumns
  !//   subroutine init_daily_diag
  !//   subroutine daily_diag_output
  !//   subroutine nops_diag
  !//   subroutine mp_diag
  !//   subroutine TBGT_2FILE
  !//   subroutine REPORTS_CHEMISTRY
  !//   subroutine sumup_burden_and_lifetime
  !//   subroutine report_burden_and_lifetime
  !//   subroutine report_ch4n2o
  !//   subroutine diag_burden_snapshot
  !//   subroutine write_snapshot
  !//   subroutine write_snapshot_the
  !//   subroutine ch4n2o_burden
  !//   subroutine ch4_loss3
  !//   subroutine n2o_loss3
  !//   subroutine diags_tpset
  !//   subroutine report_negO3
  !//   subroutine tnd_emis_daily
  !//   subroutine tnd_emis2file
  !//   subroutine save_chemPL
  !//   subroutine chembud_output
  !//
  !// Ole Amund Sovde, February 2015, August 2013, February 2012,
  !//              October 2008
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, MPBLK, IDBLK, JDBLK
  use cmn_parameters, only: M_AIR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Levels for tropopause based on O3=150ppbv
  integer, dimension(IPAR,JPAR) :: tp_o3, tp_hpa
  !// O3-tp limit as mass mixing ratio
  real(r8),parameter :: o3lim = 150.e-9_r8 * 48._r8 / M_AIR
  !// Pressure limit tp [hPa]
  real(r8),parameter :: preslim = 200._r8

  !// For instant average OH and lifetime of CH4 (and CH3CCl3?)
  real(r8), dimension(MPBLK) :: &
       RKERN_CH4_TP,   RKERN_CH4_TP_OH, &
       RKERN_CH4_LTOP, RKERN_CH4_LTOP_OH, &
       RKERN_CH4_HPA,  RKERN_CH4_HPA_OH, &
       RKERN_CO_TP,    RKERN_CO_TP_OH, &
       RKERN_CO_LTOP,  RKERN_CO_LTOP_OH, &
       RKERN_CO_HPA,   RKERN_CO_HPA_OH, &
       RKERN_SPIV,     RKERN_SPIV_OH, &
       RKERN_LAWR,     RKERN_LAWR_OH, &
       RKERN_LAWR_AIR, RKERN_LAWR_AIR_OH, &
       RKERN_LAWR_M,   RKERN_LAWR_kM, &
       RLOSS_CH4_TP, RLOSS_CH4_HPA, RLOSS_CH4_LTOP, &
       RTOT_CH4
  !// Accumulated values
  real(r8)  :: accum_oh_ch4_tp, accum_oh_ch4_hpa, accum_oh_ch4_ltop, &
       accum_oh_co_tp, accum_oh_co_hpa, accum_oh_co_ltop, &
       accum_ch4_tp, accum_ch4_hpa, accum_ch4_ltop, num_accum, &
       accum_oh_spiv, accum_oh_lawr, accum_oh_lawr_air, accum_ch4_lawr

  !// For total CH4/N2O lifetime
  !// 1: sfc-L60, 2: sfc-xxxhPa, 3: sfc-TP, 4: sfc-150ppbvO3
  real(r8), dimension(MPBLK,4) :: &
       GLINS_N2O, &          !// Instant N2O (kg)
       GLOSS_N2O_MON, &   !// Total loss (kg) per month
       GMASS_N2O_MON, &   !// Average mass per month
       GLOSS_N2O_RTOT, &  !// Total loss (kg) running total
       GMASS_N2O_RTOT     !// Average mass running total
  real(r8), dimension(MPBLK,4) :: &
       GLINS_CH4, &        !// Instant CH4 (kg)
       GLOSS_CH4_MON, &    !// Total loss (kg) monthly
       GMASS_CH4_MON, &    !// Average mass monthly
       GLOSS_CH4_RTOT, &   !// Total loss (kg) running total
       GMASS_CH4_RTOT      !// Average mass running total
  real(r8) :: num_ch4n2o, num_ch4n2o_rtot

  !// Diagnose chemistry loss
  integer, parameter :: nchemdiag = 10
  real(r8), dimension(nchemdiag,NPAR,LPAR,IPAR,JPAR) :: CHEMLOSSMASS
  real(r8), dimension(nchemdiag,NPAR,LPAR,IPAR,JPAR) :: CHEMPRODMASS
  real(r8), dimension(LPAR,IPAR,JPAR) :: OxCHEMLOSSMASS
  real(r8), dimension(LPAR,IPAR,JPAR) :: OxCHEMPRODMASS
  !// Only save prod and loss for selected components
  integer,parameter :: ncPL = 4
  integer,dimension(ncPL), parameter :: compsPL = (/13,46,113,114/)

  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='diagnostics_general.f90'
  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public diag_ohch4n2o_init, TBGT_2FILE, &
       init_daily_diag, daily_diag_output, &
       nops_diag, mp_diag, REPORTS_CHEMISTRY, sumup_burden_and_lifetimes, &
       write_snapshot, &
       ch4n2o_burden, &
       tnd_emis2file, &
       nchemdiag, save_chemPL, save_chemOxPL, chembud_output, init_chembud
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine diag_ohch4n2o_init()
    !// --------------------------------------------------------------------
    !// Called from input_oslo.f90
    !//
    !// Ole Amund Sovde, February 2011, May 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Accumulated number in average
    num_accum    = 0._r8

    !// For instant OH aveages (CH4-kernel)
    RKERN_CH4_TP(:)     = 0._r8
    RKERN_CH4_TP_OH(:)  = 0._r8
    RKERN_CH4_LTOP(:)   = 0._r8
    RKERN_CH4_LTOP_OH(:)= 0._r8
    RKERN_CH4_HPA(:)    = 0._r8
    RKERN_CH4_HPA_OH(:) = 0._r8
    accum_oh_ch4_tp   = 0._r8
    accum_oh_ch4_hpa  = 0._r8
    accum_oh_ch4_ltop = 0._r8

    !// For instant OH aveages (CO-kernel)
    RKERN_CO_TP(:)     = 0._r8
    RKERN_CO_TP_OH(:)  = 0._r8
    RKERN_CO_LTOP(:)   = 0._r8
    RKERN_CO_LTOP_OH(:)= 0._r8
    RKERN_CO_HPA(:)    = 0._r8
    RKERN_CO_HPA_OH(:) = 0._r8
    accum_oh_co_tp   = 0._r8
    accum_oh_co_hpa  = 0._r8
    accum_oh_co_ltop = 0._r8

    !// For instant OH aveages (Spivakovsky-kernel)
    RKERN_SPIV(:)     = 0._r8
    RKERN_SPIV_OH(:)  = 0._r8
    accum_oh_spiv     = 0._r8

    !// For instant OH aveages (Lawrence 2001)
    RKERN_LAWR(:)     = 0._r8
    RKERN_LAWR_OH(:)  = 0._r8
    RKERN_LAWR_AIR(:)    = 0._r8
    RKERN_LAWR_AIR_OH(:) = 0._r8
    RKERN_LAWR_M(:)   = 0._r8
    RKERN_LAWR_kM(:)  = 0._r8
    accum_oh_lawr     = 0._r8
    accum_oh_lawr_air = 0._r8
    accum_ch4_lawr    = 0._r8

    !// For instant CH4 OH-lifetime
    RLOSS_CH4_TP(:)   = 0._r8
    RLOSS_CH4_HPA(:)  = 0._r8
    RLOSS_CH4_LTOP(:) = 0._r8
    RTOT_CH4(:)   = 0._r8
    accum_ch4_tp  = 0._r8
    accum_ch4_hpa = 0._r8
    accum_ch4_ltop= 0._r8


    !// New CH4 and N2O lifetime
    GLOSS_N2O_MON(:,:) = 0._r8
    GMASS_N2O_MON(:,:) = 0._r8
    GLOSS_CH4_MON(:,:) = 0._r8
    GMASS_CH4_MON(:,:) = 0._r8
    GLOSS_CH4_RTOT(:,:) = 0._r8
    GMASS_CH4_RTOT(:,:) = 0._r8
    num_ch4n2o = 0._r8
    num_ch4n2o_rtot = 0._r8

    !// --------------------------------------------------------------------
  end subroutine diag_ohch4n2o_init
  !// ----------------------------------------------------------------------


      
  !// ----------------------------------------------------------------------
  subroutine init_lifetime()
    !// --------------------------------------------------------------------
    !// Initialize CH4 and N2O accumulated losses and burdens.
    !//
    !// Ole Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'init_lifetime'
    !// --------------------------------------------------------------------
    GMASS_N2O_MON(:,:) = 0._r8
    GMASS_CH4_MON(:,:) = 0._r8
    GLOSS_N2O_MON(:,:) = 0._r8
    GLOSS_CH4_MON(:,:) = 0._r8
    num_ch4n2o = 0._r8
    write(6,'(a)') f90file//':'//subr// &
         ': Initialized CH4/N2O lifetime diagnostics'
    !// --------------------------------------------------------------------
  end subroutine init_lifetime
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine du_columns(BTT, NMET, MP)
    !// --------------------------------------------------------------------
    !// Calculates instantaneous dobson columns (O3) for troposphere and
    !// stratosphere each NMET.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, AREAXY
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: LMTROP, dobson_snapshot, dobson_snapshot_ts, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: NMET, MP
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT

    !// Locals
    real(r8)  :: o3_col, o3_col_ts
    real(r8)  :: zo3
    integer :: I,J,L, II, JJ
    !// --------------------------------------------------------------------

    !// No O3; return
    if (trsp_idx(1) .lt. 0) return

    !// Conversion factors from mass [kg] to DU
    !// TMASS is g/mol
    !// AVOGNR is molecules/mol
    !// 1kg = 1000g
    !// 1DU = 2.69d20 molec/m2
    !// Combine kg- and DU-conversion to factor 1.d0/2.69d17
    ZO3 = AVOGNR / TMASS(trsp_idx(1)) / 2.69e17_r8


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Sum up mass
        o3_col = BTT(1,trsp_idx(1),II,JJ)
        do L = 1, LMTROP(I,J)
           o3_col = o3_col + BTT(L,trsp_idx(1),II,JJ)
        end do
        o3_col_ts = o3_col !// Tropospheric column
        do L = LMTROP(I,J)+1, LPAR
           o3_col = o3_col + BTT(L,trsp_idx(1),II,JJ)
        end do

        !// Convert from kg to DU
        dobson_snapshot(I,J,NMET)    = o3_col * ZO3 / AREAXY(I,J)
        dobson_snapshot_TS(I,J,NMET) = o3_col_ts * ZO3 / AREAXY(I,J)

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine du_columns
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine write_ducolumns_old(nday, ndayi)
    !// --------------------------------------------------------------------
    !// Write dobson columns to file. Done at the end of each day; puts out
    !// all NMETs for that day.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: NRMETD
    use cmn_chem, only: TNAME
    use utilities, only: get_free_fileid
    use cmn_oslo, only: dobson_snapshot, dobson_snapshot_ts, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: nday, ndayi

    !// File number for dobson values
    integer :: ifnr, M
    character(len=80) :: filename
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'write_ducolumns_old'
    !// --------------------------------------------------------------------

    !// If O3 not included, return
    if (trsp_idx(1) .lt. 0) return
    if (trim(TNAME(trsp_idx(1))) .ne. 'O3') return

    !// File name
    filename = 'dobson_nmet.dta'

    !// Find file number to use
    ifnr = get_free_fileid()

    if (NDAY .eq. NDAYI) then
       !// Check the size of dobson array
       if (size(dobson_snapshot,3) .lt. NRMETD) then
          write(6,'(a,i7)') f90file//':'//subr// &
               ': Wrong size of dobson_snapshot:', size(dobson_snapshot,3)
          write(6,'(a,i7)') f90file//':'//subr// &
               ': NRMETD: ',NRMETD
          stop
       end if

       !// Open file
       open(ifnr, file=filename, form='unformatted')
       !// Horizontal resolution
       write(ifnr) IPAR, JPAR
    else
       !// Open for append
       open(ifnr, file=filename, form='unformatted', position='append')
    end if

    !// Write data for this NDAY
    do M = 1, NRMETD
       write(ifnr) nday, M
       write(ifnr) real(dobson_snapshot(:,:,M), r4)
       write(ifnr) real(dobson_snapshot_ts(:,:,M), r4)
    end do
    !// Close until next day
    close(ifnr)

    write(6,'(a,i7)') '* Wrote file '//trim(filename), NDAY

    !// --------------------------------------------------------------------
  end subroutine write_ducolumns_old
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine write_ducolumns(nday, ndayi)
    !// --------------------------------------------------------------------
    !// Write dobson columns to file. Done at the end of each day; puts out
    !// all NMETs for that day.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//    Output as netCDF4.
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NRMETD, JYEAR, JMON, JDATE, XDGRD, YDGRD, XDEDG, YDEDG
    use cmn_chem, only: TNAME
    use cmn_diag, only: nc4deflate_global, nc4shuffle_global
    use cmn_oslo, only: dobson_snapshot, dobson_snapshot_ts, trsp_idx
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: nday, ndayi

    !// Local variables
    character(len=80) :: filename
    character(len=8) :: datestamp
    integer, dimension(NRMETD) :: time_hrs_since_start
    integer :: &
         N, &
         lat_dim_id, lon_dim_id, time_dim_id, &
         lat_id, lon_id, time_id, &
         ilat_dim_id, ilon_dim_id, &
         ilat_id, ilon_id, &
         dutot_id, dutrp_id, &
         start_year_id, start_month_id, start_date_id, &
         nlats, nlons, nsteps, status, ncid, &
         srt_lon_lat_time(3), &
         cnt_lon_lat_time(3)

    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'write_ducolumns'
    !// --------------------------------------------------------------------

    !// If O3 not included, return
    if (trsp_idx(1) .lt. 0) return
    if (trim(TNAME(trsp_idx(1))) .ne. 'O3') return

    !// File name
    filename = 'dobson_nmet.nc'

    if (NDAY .eq. NDAYI) then
       !// Check the size of dobson array
       if (size(dobson_snapshot,3) .lt. NRMETD) then
          write(6,'(a,i7)') f90file//':'//subr// &
               ': Wrong size of dobson_snapshot:', size(dobson_snapshot,3)
          write(6,'(a,i7)') f90file//':'//subr// &
               ': NRMETD: ',NRMETD
          stop
       end if

       !// initialize number of time steps stored
       nsteps = 0

       !// Create file
       write(6,'(a)') f90file//':'//subr//': creating file: '//trim(filename)
       status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': in creating file')

       !//File headers
       status=nf90_put_att(ncid,nf90_global,'title', &
            '3h O3 column (DU) from Oslo CTM3')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': file header')

       !// Define spatial dimensions (lat, lon)
       status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lat dim')
       status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lon dim')
       !// Define spatial dimensions (ilat, ilon)
       status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define ilat dim')
       status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define ilon dim')
       !// Define time dimension unlimited
       status = nf90_def_dim(ncid,"time",nf90_unlimited,time_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define time dim')

       !// Define the lon/lat
       status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lon variable')
       status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define lat variable')
       !// Define ilon/ilat
       status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define ilon variable')
       status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define ilat variable')
       !// Define time
       status = nf90_def_var(ncid,"time",nf90_double,time_dim_id,time_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define time variable')

       !// Putting attributes to lon/lat variables
       status = nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units lon')
       status = nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units lat')
       !// Putting attributes to ilon/ilat variables
       status = nf90_put_att(ncid,ilon_id,'units','degree_east')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units ilon')
       status = nf90_put_att(ncid,ilat_id,'units','degree_north')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units ilat')

       !// Putting attributes to time variable, defined to be
       !// hours since model start.
       write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
       status = nf90_put_att(ncid,time_id,'units', &
            'hours since model start '//datestamp)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units time')

       status = nf90_def_var(ncid,"start_year",nf90_int,start_year_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define start_year variable')
       status = nf90_def_var(ncid,"start_month",nf90_int,start_month_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define start_month variable')
       status = nf90_def_var(ncid,"start_date",nf90_int,start_date_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define start_date variable')


       !// Define total O3 column variables
       status = nf90_def_var(ncid,'o3col_tot',nf90_float, &
            (/lon_dim_id,lat_dim_id,time_dim_id/), dutot_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define o3col_tot variable')
       !// Deflate netcdf4
       status = nf90_def_var_deflate(ncid,dutot_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define o3col_tot variable deflate')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,dutot_id,'units','DU')        
       if (status .ne. nf90_noerr) call handle_error(status,&
            f90file//':'//subr//': attribute units o3col_tot')
       status = nf90_put_att(ncid,dutot_id, 'long_name', &
            'Total O3 column [DU]')
       if (status .ne. nf90_noerr) call handle_error(status,&
            f90file//':'//subr//': attribute long_name o3col_tot')

       !// Define tropospheric O3 column variables
       status = nf90_def_var(ncid,'o3col_trp',nf90_float, &
            (/lon_dim_id,lat_dim_id,time_dim_id/), dutrp_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define o3col_trp variable')
       !// Deflate netcdf4
       status = nf90_def_var_deflate(ncid,dutrp_id,nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define o3col_trp variable deflate')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,dutrp_id,'units','DU')        
       if (status .ne. nf90_noerr) call handle_error(status,&
            f90file//':'//subr//': attribute units o3col_trp')
       status = nf90_put_att(ncid,dutrp_id, 'long_name', &
            'Tropospheric O3 column [DU]')
       if (status .ne. nf90_noerr) call handle_error(status,&
            f90file//':'//subr//': attribute long_name o3col_trp')

       !// -----------------------------------------------------------------
       !// End definition mode
       status = nf90_enddef(ncid)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': end defmode')
       !// -----------------------------------------------------------------

       !// Put lon/lat
       status = nf90_put_var(ncid,lon_id,XDGRD)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lon')
       status = nf90_put_var(ncid,lat_id,YDGRD)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lat')
       !// Put ilon/ilat
       status = nf90_put_var(ncid,ilon_id,XDEDG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lon')
       status = nf90_put_var(ncid,ilat_id,YDEDG)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting lat')

       !// Start date info
       status = nf90_put_var(ncid,start_year_id,JYEAR)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting start_year')
       status = nf90_put_var(ncid,start_month_id,JMON)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting start_month')
       status = nf90_put_var(ncid,start_date_id,JDATE)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting start_date')

       !// -----------------------------------------------------------------
       !// file will be closed after end of this if-else statement
       !// -----------------------------------------------------------------

    else

       !// Not first time step - open existing file for append
       status = nf90_open(filename, nf90_write, ncid)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': open existing file')

       !// Inquire dimension IDs
       status = nf90_inq_dimid(ncid,"lat",lat_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting lat')
       status = nf90_inq_dimid(ncid,"lon",lon_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting lon')
       status = nf90_inq_dimid(ncid,"time",time_dim_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': getting time')

       !// no need to check size dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq lat dim')
       if (nlats .ne. JPAR) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': wrong resolution: on file vs JPAR: ',nlats,JPAR
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq lon dim')
       if (nlons .ne. IPAR) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': wrong resolution: on file vs IPAR: ',nlons,IPAR
          stop
       end if
       !// Check time step
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq time dim')

       if (nsteps.le.0) then
          write(6,'(a)') f90file//':'//subr// &
               ': file reports already added nsteps = ',nsteps
          stop
       end if

       !// Get variable IDs for O3-columns
       status = nf90_inq_varid(ncid,'o3col_tot',dutot_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq varid dutot_id')
       status = nf90_inq_varid(ncid,'o3col_trp',dutrp_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq varid dutrp_id')

       !// Get variable ID for time
       status = nf90_inq_varid(ncid,"time",time_id) 
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': inq varid time_id')
    end if



    !// Put data to file:
    !// --------------------------------------------------------------------
    !// Defining how far to count for each time a data set is added
    cnt_lon_lat_time = (/IPAR, JPAR,  NRMETD/)
    !// Defining where to start adding the new time step
    srt_lon_lat_time = (/1, 1, nsteps + 1/)

    !// time since model start
    !//   nsteps is the previous time step added
    !//   nsteps + 1 is start index of the new time steps added 
    do N = 1, NRMETD
       time_hrs_since_start(N) = &
            real(nsteps + N - 1, r8)*24._r8 / real(NRMETD, r8)
    end do
    status = nf90_put_var(ncid,time_id,time_hrs_since_start, &
         start = (/ nsteps + 1 /), count = (/ NRMETD /))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting time')

    status = nf90_put_var(ncid,dutot_id,dobson_snapshot(:,:,:), &
            start=srt_lon_lat_time, count=cnt_lon_lat_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting o3col_tot')

    status = nf90_put_var(ncid,dutrp_id,dobson_snapshot_ts(:,:,:), &
            start=srt_lon_lat_time, count=cnt_lon_lat_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting o3col_trp')

    !// --------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file')
    !// --------------------------------------------------------------------

    write(6,'(a,i7)') f90file//':'//subr// &
         ': wrote file '//trim(filename), NDAY

    !// --------------------------------------------------------------------
  end subroutine write_ducolumns
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine init_daily_diag(JYEAR, JMON, JDAY, NDAY, NDAYI)
    !// --------------------------------------------------------------------
    !// Do daily initializations of diagnostics.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_diag, only: TROPMASS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR, JMON, JDAY, NDAY, NDAYI
    !// --------------------------------------------------------------------

    !// set TROPMASS every day
    TROPMASS(:,:,:) = 0._r8

    !// --------------------------------------------------------------------
  end subroutine init_daily_diag
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine daily_diag_output(DDYEAR, DDMON, DDDATE, DDAY, NDAY, NDAYI)
    !// --------------------------------------------------------------------
    !// Do daily diagnostics at the end of the day, before calendar update.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use diagnostics_scavenging, only: scav_diag_collect_daily, &
         scav_diag_2fileB
    use cmn_ctm, only: LDLYSCAV
    !use emissions_megan, only: megan_report
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: DDYEAR, DDMON, DDDATE, DDAY, NDAY, NDAYI
    !// --------------------------------------------------------------------

    !// Write dobson columns to file
    call write_ducolumns(nday, ndayi)

    !// Collect wet scavenging diagnose
    call scav_diag_collect_daily(DDAY)

    !// Daily output for wet deposition maps
    if (LDLYSCAV(2)) then
       call scav_diag_2fileB(DDYEAR,DDMON,DDDATE,NDAY)
    end if

    !// Print MEGAN report (can be added temporarily if needed)
    !call megan_report(NDAY, NDAYI)

    !// Collect total emissions this day
    call tnd_emis_daily(DDAY)


    !// --------------------------------------------------------------------
  end subroutine daily_diag_output
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine nops_diag(JYEAR, JMON, JDATE, NDAY, NMET, NOPS, NDAYI, LNEWM)
    !// --------------------------------------------------------------------
    !// Do diagnostics at the beginning of each NOPS.
    !// Not parallelized!
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LBCOC, LOSLOCTROP, LOSLOCSTRAT
    use cmn_sfc, only: LDDEPmOSaic, VGSTO3, VRAO3, VRBO3, VRCO3
    use bcoc_oslo, only: bcsnow_nmet_output, bcsnow_nmet_output_nc
    use ch4routines, only: reportsfcch4
    use satelliteprofiles_mls, only: satprofs_mls_master
    use verticalprofiles_stations2, only: vprofs_master
    !use atom, only: atom_master
    use caribic2, only: caribic2_master
    use hippo, only: hippo_master
    use troccinox_fal, only: troccifal_master
    use troccinox_geo, only: troccigeo_master
    use troccinox_ban, only: trocciban_master
    use diagnostics_scavenging, only: scav_diag_put_gsto, scav_diag_put_drydepvelo
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI
    logical, intent(in) :: LNEWM
    !// --------------------------------------------------------------------

    !// Update "tropopauses" for diagnostics (not LMTROP)
    call diags_tpset()

    !// BCsnow output
    if (LBCOC) call bcsnow_nmet_output_nc(NDAY,NMET,NOPS)

    !// Process satellite profiles
    !if (LOSLOCTROP .or. LOSLOCSTRAT) &
    !     call satprofs_mls_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI)

    !// Process station profiles
    call vprofs_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI)

    !// Process flight data caribic (10 second data)
    !if (LOSLOCTROP .or. LOSLOCSTRAT) &
    !     call caribic2_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI,LNEWM)

    !// Process flight data hippo (10 second data)
    !// Do this for BC only, but can be set up to do other species
    !// as well.
    if (LBCOC) &
         call hippo_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI,LNEWM)

    !if (LBCOC) &
         !call atom_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI,LNEWM)

    !// Process flight data for troccinox
    !if (LOSLOCTROP .or. LOSLOCSTRAT) then
    !   call troccifal_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI)
    !   call troccigeo_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI)
    !   call trocciban_master(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI)
    !end if

    !// Surface CH4
    if (LOSLOCTROP .and. NOPS .eq. 1) call reportsfcch4()

    !// Snapshots on theta levels (mode=1 is NH)
    !call write_snapshot_the(nday,nmet,nops,1)

    !// Snapshot from stomatal conductance
    if (LDDEPmOSaic) then
       call scav_diag_put_gsto(NMET,VGSTO3)
       call scav_diag_put_drydepvelo(NMET,VRAO3,VRBO3,VRCO3)
    end if
    !// --------------------------------------------------------------------
  end subroutine nops_diag
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine mp_diag(BTT,NDAY,NMET,NOPS,NSUB,CCYC,MP)
    !// --------------------------------------------------------------------
    !// Do diagnoses in parallel, i.e. for a IJ-block.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: IDBLK, JDBLK, NPAR, LPAR
    use cmn_sfc, only: LDDEPmOSaic, VGSTO3
    use diagnostics_scavenging, only: scav_diag_put_fsto
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NMET, NOPS, NSUB, CCYC, MP
    real(r8), intent(in)  :: BTT(LPAR,NPAR,IDBLK,JDBLK)
    !// --------------------------------------------------------------------

    !// Dobson columns
    if (NOPS.eq.1. .and. NSUB.eq.1. .and. CCYC.eq.1) then
       call du_columns(BTT,NMET,MP)
      
    end if
    if (NSUB.eq.1. .and. CCYC.eq.1) then
       !// Snapshot from stomatal flux
       if (LDDEPmOSaic) call scav_diag_put_fsto(BTT,NMET,VGSTO3,MP)
    end if
    !// --------------------------------------------------------------------
  end subroutine mp_diag
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TBGT_2FILE(INIT)
    !// --------------------------------------------------------------------
    !// Processes horizontal (column sums) 2-D tendency budgets for
    !// specified boxes.
    !// Dumps to file
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_size, only: NTDPAR
    use cmn_ctm, only: NTM, AIR, STT, GMTAU, JDATE, JMON, JYEAR, &
         MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, YDGRD, IDAY
    use cmn_diag, only: LBGA2, LBGT2, NTND, STTTND, STTTN0, &
         TAU0, JDATE0, JMON0, JYEAR0, NDAY0, TLDIAG
    use utilities, only: get_free_fileid
    use cmn_oslo, only: CONVWASHOUT, chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  INIT
    !// Locals
    real(r8) :: DT_HR
    integer :: I,J,II,JJ,L,M,MP,M1,N, ifnr
    !// For file stuff
    character(len=80) :: BTITLE, filenameij,filenamejl
    !// Map (IJ)
    real(r8) :: BGTIJ(IPAR,JPAR,NTDPAR+1,NPAR)
    real(r8) :: BGTIJA(IPAR,JPAR),BGTIJT(IPAR,JPAR)
    !// Zonal mean
    real(r8) :: BGTJL(JPAR,LPAR,NTDPAR+1,NPAR)
    real(r8) :: BGTJLA(JPAR,LPAR),BGTJLT(JPAR,LPAR)
    !// For total global output
    real(r8) :: GLOB_BUDGET(NTDPAR+1,NPAR)
    !// Names of diagnosed processes
    character(len=6), dimension(NTDPAR+1) :: PTITLE
    !// Version number
    integer, parameter :: VERSION = 2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tbgt_2file'
    !// --------------------------------------------------------------------

    DT_HR = 24._r8 * real(IDAY-NDAY0, r8)
    M1 = NTND+1               !// Convective washout

    BGTIJ(:,:,:,:) = 0._r8
    BGTIJA(:,:) = 0._r8
    BGTIJT(:,:) = 0._r8

    BGTJL(:,:,:,:) = 0._r8
    BGTJLA(:,:) = 0._r8
    BGTJLT(:,:) = 0._r8

    PTITLE(1:NTND)  = TLDIAG(1:NTND)
    PTITLE(M1) = 'CNV WO'

    !// Files to write to
    filenameij = 'ctm3_budgets_ij.dta'
    filenamejl = 'ctm3_budgets_jl.dta'
    !// Find non-used file number for input file
    ifnr = get_free_fileid()

    !// Initialize?
    if (INIT .eq. 0) then
       !// IJ-budgets
       write(6,'(a)') f90file//':'//subr// &
            ': Initializing budget file '//trim(filenameij)
       open(ifnr,file=filenameij,form='unformatted')
       write(ifnr) NTM, M1,IPAR,JPAR,LPAR, VERSION
       write(ifnr) LBGA2, LBGT2(1:NTM)
       close(ifnr)
       !// JL-budgets
       write(6,'(a)') f90file//':'//subr// &
            ': Initializing budget file '//trim(filenamejl)
       open(ifnr,file=filenamejl,form='unformatted')
       write(ifnr) NTM, M1,IPAR,JPAR,LPAR, VERSION
       write(ifnr) LBGA2, LBGT2(1:NTM)
       close(ifnr)
       !// Skip diagnostics when initializing
       return
    end if


    !// AIR (N=0) calculate 2-D budget tendencies = ONLY total inst mass
    if (LBGA2) then
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                BGTIJA(I,J) = BGTIJA(I,J) + AIR(I,J,L)
                BGTJLA(J,L) = BGTJLA(J,L) + AIR(I,J,L)
             end do
          end do
       end do
    end if

    !// Calculate 2-D tendency budgets (if .true.) for all boxes for tracer N
    do N = 1, NTM
       if (LBGT2(N)) then

          !// total instantaneous tracer mass in box
          do L = 1, LPAR
             do J = 1, JPAR
                do I = 1, IPAR
                   BGTIJT(I,J) = BGTIJT(I,J) + STT(I,J,L,N)
                   BGTJLT(J,L) = BGTJLT(J,L) + STT(I,J,L,N)
                end do
             end do
          end do
       end if
    end do

    !// Acccumulate all 2-D tendencies 1:NTND for the box
    do M = 1, NTND
       do N = 1, NTM
          if (LBGT2(N)) then
             do L = 1, LPAR
                do J = 1, JPAR
                   do I = 1, IPAR
                      BGTIJ(I,J,M,N) = BGTIJ(I,J,M,N) + STTTND(I,J,L,N,M)
                      BGTJL(J,L,M,N) = BGTJL(J,L,M,N) + STTTND(I,J,L,N,M)
                   end do
                end do
             end do
          end if
       end do
    end do

    !// Sum up global budget numbers (STTTN0)
    GLOB_BUDGET(:,:) = 0._r8
    do M = 1, NTND
       do N = 1, NTM
          !// Switch indices
          GLOB_BUDGET(M,N) = STTTN0(N,M)
       end do
    end do

    !// Convective washout
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do N = 1, NTM
            if (LBGT2(N)) then
              do L = 1, LPAR
                BGTIJ(I,J,M1,N) = BGTIJ(I,J,M1,N) + CONVWASHOUT(L,N,II,JJ,MP)
                BGTJL(J,L,M1,N) = BGTJL(J,L,M1,N) + CONVWASHOUT(L,N,II,JJ,MP)
                !// Global budget
                GLOB_BUDGET(M1,N) = GLOB_BUDGET(M1,N) &
                                    + CONVWASHOUT(L,N,II,JJ,MP)
              end do
            end if
          end do
        end do
      end do
    end do

    !// open file
    write(6,'(a)') f90file//':'//subr//': Writing to '//trim(filenameij)
    !// File is already generated
    open(ifnr,file=filenameij,form='unformatted',position='append')


    !// Write
    write(ifnr) real(TAU0, r4),JDATE0,JMON0,JYEAR0
    write(ifnr) real(GMTAU, r4),JDATE,JMON,JYEAR
    write(ifnr) DT_HR
    write(ifnr) PTITLE
    if (LBGA2) write(ifnr) real(BGTIJA, r4)
    write(ifnr) real(BGTIJT, r4)
    do N = 1, NTM
       if (LBGT2(N)) then
          write(ifnr) chem_idx(N)
          write(ifnr) real(BGTIJ(:,:,1:M1,N), r4)
       end if
    end do
    write(ifnr) GLOB_BUDGET !// All tracers, real(r8)
    close(ifnr)

    !// open file
    write(6,'(a)') f90file//':'//subr//': Writing to '//trim(filenamejl)
    !// File is already generated
    open(ifnr,file=filenamejl,form='unformatted',position='append')

    !// Write
    write(ifnr) real(TAU0, r4),JDATE0,JMON0,JYEAR0
    write(ifnr) real(GMTAU, r4),JDATE,JMON,JYEAR
    write(ifnr) DT_HR
    write(ifnr) PTITLE
    if (LBGA2) write(ifnr) real(BGTJLA, r4)
    write(ifnr) real(BGTJLT, r4)
    do N = 1, NTM
       if (LBGT2(N)) then
          write(ifnr) chem_idx(N)
          write(ifnr) real(BGTJL(:,:,1:M1,N), r4)
       end if
    end do
    write(ifnr) GLOB_BUDGET !// All tracers, real(r8)
    close(ifnr)

    !// --------------------------------------------------------------------
  end subroutine TBGT_2FILE
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine REPORTS_CHEMISTRY(NDAY,NMET,NOPS,TAU, NDAYI,LNEWM)
    !// --------------------------------------------------------------------
    !// Prints out report for the chemistry & physics.
    !// Mainly for global reports on IJ-blocks.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LOSLOCTROP
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS, NDAYI
    real(r8), intent(in)  :: TAU
    logical, intent(in) :: LNEWM
    !// --------------------------------------------------------------------

    !// Print TP info
    if (NOPS .eq. 1) then
       write(6,'(a,i5,2i3,1x,2i3)') '* Tropopause ',NDAY, &
            NMET,NOPS,minval(LMTROP),maxval(LMTROP)
    end if

    !// Print number of negative ozone corrected in tropchem.
    call report_negO3()

    if (LOSLOCTROP .and. NOPS.eq.1) then
       !// Print out CH4 lifetime and OH burden (3hr avg and running mean)
       call report_burden_and_lifetime(NDAY,NMET,NOPS,TAU, &
            NDAY.ne.NDAYI.and.LNEWM)
       !// Print out CH4 and N2O lifetime (instant/monthly/running mean)
       call report_ch4n2o(NDAY,NMET,NOPS,TAU, &
            NDAY.ne.NDAYI.and.LNEWM)
    end if

    if (NOPS.eq.1) then
       !// Burden diagnostics (mass)
       call diag_burden_snapshot(nday,nmet,nops)
    end if

    !// --------------------------------------------------------------------
  end subroutine REPORTS_CHEMISTRY
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine sumup_burden_and_lifetimes(BTT,AIRB,BVOL,BTEM,MP)
    !// --------------------------------------------------------------------
    !// Find average OH concentration based on CH4-kernel and CO-kernel.
    !// Also find OH-lifetime based on OH and reaction rate with CH4.
    !// Domains for calculation:
    !//   - surface up to LMTROP
    !//   - air mass below 200 hPa
    !//   - from level 1:LPAR-1
    !//
    !// Sums up the numbers in each IJ-block.
    !// Will calculate zero values if components are not included.
    !//
    !// Different from CTM2 method, which made no sense for OH.
    !// CTM2 also calculated wrong CH4 lifetime (did not use total burden).
    !//
    !// Ole Amund Sovde, May 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK, LOSLOCTROP
    use cmn_ctm, only: ETAA, ETAB, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         YDGRD, YGRD
    use cmn_chem, only: TMASS
    use cmn_met, only: P
    use cmn_parameters, only: AVOGNR, M_AIR, MINTEMP
    use cmn_oslo, only: LMTROP, trsp_idx, Xtrsp_idx, XTMASS, XSTT
    use chem_oslo_rates, only: r_oh_ch4
    use utilities_oslo, only: rate3B
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: AIRB, BTEM
    real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK), intent(in) :: BVOL
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT

    !// Local Parameters
    real(r8) :: &
         MOLEC_OH, CONC_OH, &
         MOLEC_CH4, CONC_CH4, MASS_CH4, &
         CONC_AIR, &
         !// OH average CH4-kernel
         KERN_CH4_TP,  KERN_CH4_TP_OH, &
         KERN_CH4_LTOP,KERN_CH4_LTOP_OH, &
         KERN_CH4_HPA, KERN_CH4_HPA_OH, &
         !// OH average CO-kernel
         KERN_CO_TP,  KERN_CO_TP_OH, &
         KERN_CO_LTOP,KERN_CO_LTOP_OH, &
         KERN_CO_HPA, KERN_CO_HPA_OH, &
         !// Spivakovsky kernel (Spivakovsky et al, 2000, JGR 105(D7), 8931)
         KERN_SPIV, KERN_SPIV_OH, pspiv, &
         !// Lawrence et al, 2001, ACP, doi:10.5194/acp-1-37-2001
         KERN_LAWR, KERN_LAWR_OH, plawr, &
         KERN_LAWR_AIR, KERN_LAWR_AIR_OH, &
         KERN_LAWR_M, KERN_LAWR_kM, &
         !// CH4 lifetime
         LOSS_CH4_TP, LOSS_CH4_HPA, LOSS_CH4_LTOP, & !// CH4 lifetime
         TOT_CH4, &
         !// Other variables
         RDUM, RVCM3, RCOOH, R2, P1, P2, PFRAC
    integer :: I,J,L,TI,II,JJ
    !// --------------------------------------------------------------------

    !// If not tropospheric chemistry, return
    if (.not. LOSLOCTROP) return

    !// Initialize
    KERN_CH4_TP      = 0._r8
    KERN_CH4_TP_OH   = 0._r8
    KERN_CH4_LTOP    = 0._r8
    KERN_CH4_LTOP_OH = 0._r8
    KERN_CH4_HPA     = 0._r8
    KERN_CH4_HPA_OH  = 0._r8

    KERN_CO_TP      = 0._r8
    KERN_CO_TP_OH   = 0._r8
    KERN_CO_LTOP    = 0._r8
    KERN_CO_LTOP_OH = 0._r8
    KERN_CO_HPA     = 0._r8
    KERN_CO_HPA_OH  = 0._r8

    KERN_SPIV    = 0._r8
    KERN_SPIV_OH = 0._r8

    KERN_LAWR    = 0._r8
    KERN_LAWR_OH = 0._r8
    KERN_LAWR_AIR    = 0._r8
    KERN_LAWR_AIR_OH = 0._r8
    KERN_LAWR_M  = 0._r8
    KERN_LAWR_kM = 0._r8

    LOSS_CH4_TP   = 0._r8
    LOSS_CH4_HPA  = 0._r8
    LOSS_CH4_LTOP = 0._r8
    TOT_CH4 = 0._r8

    !// Factor for converting from mass (kg) to molecules
    !// molec = mass(kg) / (1d-3kg/g g/mol) * Na
    !//       = mass(kg) * 1d3 * Na
    RDUM = 1.e3_r8 * AVOGNR
      
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Pressure for Spivakovsky kernel
      if (abs(YDGRD(J)) .gt. 32._r8) then
        pspiv = 200._r8
      else
        pspiv = 100._r8
      end if
      !// Pressure for Lawrence calculation
      plawr = 300._r8 - 215._r8 * (cos(YGRD(J))**2)

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        do L = 1, LPAR-1
          !// Temperature index for reaction rate
          TI   = Nint(BTEM(L,II,JJ)) - MINTEMP

          !// Volume [cm3]
          RVCM3 = 1.e6_r8 * BVOL(L,II,JJ,MP)


          !// Get OH molecules
          if (Xtrsp_idx(40) .gt. 0) then
             MOLEC_OH = XSTT(L,Xtrsp_idx(40),I,J) &
                        / XTMASS(Xtrsp_idx(40)) * RDUM
          else if (trsp_idx(40) .gt. 0) then
             MOLEC_OH = BTT(L,trsp_idx(40),II,JJ) &
                        / TMASS(trsp_idx(40)) * RDUM
          else
             !// OH is not included
             MOLEC_OH = 0._r8
          end if
          !// Get OH concentration [molec/cm3]
          CONC_OH  = MOLEC_OH / RVCM3

          !// CH4 molecules
          if (trsp_idx(46) .gt. 0) then
             MOLEC_CH4 = BTT(L,trsp_idx(46),II,JJ) &
                         / TMASS(trsp_idx(46)) * RDUM
             MASS_CH4 = BTT(L,trsp_idx(46),II,JJ)
          else if (Xtrsp_idx(46) .gt. 0) then
             MOLEC_CH4 = XSTT(L,Xtrsp_idx(46),I,J) &
                         / XTMASS(Xtrsp_idx(46)) * RDUM
             MASS_CH4 = XSTT(L,Xtrsp_idx(46),I,J)
          else
             !// CH4 is not included
             MOLEC_CH4 = 0._r8
             MASS_CH4 = 0._r8
          end if
          !// CH4 concentration [molec/cm3]
          CONC_CH4  = MOLEC_CH4 / RVCM3


          !// OH burden will be calculated with CH4 kernel and CO kernel,
          !// as described in subroutine header.

          !// Reaction rate for CO+OH (IUPAC)
          !RCOOH = AIRB(L,II,JJ) * RDUM / (M_AIR * RVCM3)
          !RCOOH = 1.44d-13*(1._r8 + RCOOH/4.0d19)

          !// JPL
          CONC_AIR = AIRB(L,II,JJ) * RDUM / (M_AIR * RVCM3)
          RCOOH = rate3B(110, BTEM(L,II,JJ)/300._r8, CONC_AIR, &
               5.9e-33_r8, 1.4_r8, 1.1e-12_r8, -1.3_r8, 0.6_r8, 0)
          R2 = rate3B(111, BTEM(L,II,JJ)/300._r8, CONC_AIR, &
               1.5e-13_r8, -0.6_r8, 2.9e9_r8, -6.1_r8, 0.6_r8, 1)
          RCOOH = RCOOH + R2



          !// Below LMTROP
          if (L .lt. LMTROP(I,J)) then
             !// OH average CH4-kernel
             KERN_CH4_TP = KERN_CH4_TP + AIRB(L,II,JJ) * r_oh_ch4(TI)
             KERN_CH4_TP_OH = KERN_CH4_TP_OH &
                  + AIRB(L,II,JJ) * r_oh_ch4(TI) * CONC_OH
             !// OH average CO-kernel
             KERN_CO_TP = KERN_CO_TP + AIRB(L,II,JJ) * RCOOH
             KERN_CO_TP_OH = KERN_CO_TP_OH &
                  + AIRB(L,II,JJ) * RCOOH * CONC_OH

             !// OH-loss of CH4 (sum up molecules, not concentration)
             LOSS_CH4_TP = LOSS_CH4_TP + r_oh_ch4(TI) * CONC_OH * MOLEC_CH4
          end if


          !// Up to xxxhPa pressure level
          P1 = ETAA(L) + ETAB(L) * P(I,J)
          P2 = ETAA(L+1) + ETAB(L+1) * P(I,J)
          !// Find fraction to add
          if (P2 .ge. preslim) then
             !// Add whole gridbox
             PFRAC = 1._r8
          else if (p1.gt.preslim .and. p2.lt.preslim) then
             !// Fraction between P1 and preslim
             PFRAC = (P1 - preslim)/(P1 - P2)
          else
             !// Add nothing (will not occur in the next if-statement)
             PFRAC = 0._r8
          end if

          !// With PFRAC, we can add as long as box bottom pressure is
          !// greater than preslim:
          if (P1 .ge. preslim) then
             !// OH average CH4-kernel
             KERN_CH4_HPA = KERN_CH4_HPA &
                  + AIRB(L,II,JJ) * r_oh_ch4(TI) * PFRAC
             KERN_CH4_HPA_OH = KERN_CH4_HPA_OH &
                  + AIRB(L,II,JJ) * r_oh_ch4(TI) * CONC_OH * PFRAC
             !// OH average CO-kernel
             KERN_CO_HPA = KERN_CO_HPA &
                  + AIRB(L,II,JJ) * RCOOH * PFRAC
             KERN_CO_HPA_OH = KERN_CO_HPA_OH &
                  + AIRB(L,II,JJ) * RCOOH * CONC_OH * PFRAC

             !// OH-loss of CH4 (sum up molecules, not concentration)
             LOSS_CH4_HPA = LOSS_CH4_HPA &
                  + r_oh_ch4(TI) * CONC_OH * MOLEC_CH4 * PFRAC
          end if

          !// Spivakovsky loss uses pressure pspiv
          if (P2 .ge. pspiv) then
             PFRAC = 1._r8
          else if (p1.gt.pspiv .and. p2.lt.pspiv) then
             PFRAC = (P1 - pspiv) / (P1 - P2)
          else
             PFRAC = 0._r8
          end if
          if (P1 .ge. pspiv) then
             KERN_SPIV = KERN_SPIV + AIRB(L,II,JJ) * PFRAC
             KERN_SPIV_OH = KERN_SPIV_OH + AIRB(L,II,JJ) * CONC_OH * PFRAC
          end if

          !// Lawrence
          if (P2 .ge. plawr) then
             PFRAC = 1._r8
          else if (p1.gt.plawr .and. p2.lt.plawr) then
             PFRAC = (P1 - plawr) / (P1 - P2)
          else
             PFRAC = 0._r8
          end if
          if (P1 .ge. plawr) then
             KERN_LAWR = KERN_LAWR + r_oh_ch4(TI) * PFRAC
             KERN_LAWR_OH = KERN_LAWR_OH + r_oh_ch4(TI) * CONC_OH * PFRAC
             KERN_LAWR_AIR = KERN_LAWR_AIR + AIRB(L,II,JJ) * PFRAC
             KERN_LAWR_AIR_OH = KERN_LAWR_AIR_OH &
                                + AIRB(L,II,JJ) * CONC_OH * PFRAC
             !// For calculating CH4 lifetime
             KERN_LAWR_M = KERN_LAWR_M + MASS_CH4 * PFRAC
             KERN_LAWR_kM = KERN_LAWR_kM + MASS_CH4 * r_oh_ch4(TI) * PFRAC
          end if


          !// Up to LTOP
          !// OH average CH4-kernel
          KERN_CH4_LTOP = KERN_CH4_LTOP + AIRB(L,II,JJ) * r_oh_ch4(TI)
          KERN_CH4_LTOP_OH = KERN_CH4_LTOP_OH &
               + AIRB(L,II,JJ) * r_oh_ch4(TI) * CONC_OH
          !// OH average CO-kernel
          KERN_CO_LTOP = KERN_CO_LTOP + AIRB(L,II,JJ) * RCOOH
          KERN_CO_LTOP_OH = KERN_CO_LTOP_OH + AIRB(L,II,JJ) * RCOOH * CONC_OH

          !// Up to LTOP: OH-loss of CH4 (sum up molecules, not concentration)
          LOSS_CH4_LTOP = LOSS_CH4_LTOP &
               + r_oh_ch4(TI) * CONC_OH * MOLEC_CH4
          !// Total number of CH4 molecules (i.e. burden)
          TOT_CH4 = TOT_CH4 + MOLEC_CH4

        end do
      end do
    end do

    !// Save values in MP-arrays

    !// Tropospheric
    RKERN_CH4_TP(MP) = RKERN_CH4_TP(MP) + KERN_CH4_TP
    RKERN_CH4_TP_OH(MP) = RKERN_CH4_TP_OH(MP) + KERN_CH4_TP_OH
    RKERN_CO_TP(MP) = RKERN_CO_TP(MP) + KERN_CO_TP
    RKERN_CO_TP_OH(MP) = RKERN_CO_TP_OH(MP) + KERN_CO_TP_OH
    RLOSS_CH4_TP(MP) = RLOSS_CH4_TP(MP) + LOSS_CH4_TP

    !// Up to some pressure level
    RKERN_CH4_HPA(MP) = RKERN_CH4_HPA(MP) + KERN_CH4_HPA
    RKERN_CH4_HPA_OH(MP) = RKERN_CH4_HPA_OH(MP) + KERN_CH4_HPA_OH
    RKERN_CO_HPA(MP) = RKERN_CO_HPA(MP) + KERN_CO_HPA
    RKERN_CO_HPA_OH(MP) = RKERN_CO_HPA_OH(MP) + KERN_CO_HPA_OH
    RLOSS_CH4_HPA(MP) = RLOSS_CH4_HPA(MP) + LOSS_CH4_HPA

    RKERN_SPIV(MP) = RKERN_SPIV(MP) + KERN_SPIV
    RKERN_SPIV_OH(MP) = RKERN_SPIV_OH(MP) + KERN_SPIV_OH

    RKERN_LAWR(MP) = RKERN_LAWR(MP) + KERN_LAWR
    RKERN_LAWR_OH(MP) = RKERN_LAWR_OH(MP) + KERN_LAWR_OH
    RKERN_LAWR_AIR(MP) = RKERN_LAWR_AIR(MP) + KERN_LAWR_AIR
    RKERN_LAWR_AIR_OH(MP) = RKERN_LAWR_AIR_OH(MP) + KERN_LAWR_AIR_OH
    RKERN_LAWR_M(MP) = RKERN_LAWR_M(MP) + KERN_LAWR_M
    RKERN_LAWR_kM(MP) = RKERN_LAWR_kM(MP) + KERN_LAWR_kM

    !// Whole atmosphere
    RKERN_CH4_LTOP(MP) = RKERN_CH4_LTOP(MP) + KERN_CH4_LTOP
    RKERN_CH4_LTOP_OH(MP) = RKERN_CH4_LTOP_OH(MP) + KERN_CH4_LTOP_OH
    RKERN_CO_LTOP(MP) = RKERN_CO_LTOP(MP) + KERN_CO_LTOP
    RKERN_CO_LTOP_OH(MP) = RKERN_CO_LTOP_OH(MP) + KERN_CO_LTOP_OH
    RLOSS_CH4_LTOP(MP) = RLOSS_CH4_LTOP(MP) + LOSS_CH4_LTOP
    RTOT_CH4(MP) = RTOT_CH4(MP) + TOT_CH4



    !// --------------------------------------------------------------------
  end subroutine sumup_burden_and_lifetimes
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine report_burden_and_lifetime(NDAY,NMET,NOPS,TAU,zero)
    !// --------------------------------------------------------------------
    !// Print out average OH and lifetimes.
    !// Sums up contributions from all MP-blocks.
    !// Values are initialized in diag_init called from input_oslo.f90.
    !//
    !// Does both OH and lifetimes from sumup_burden_and_lifetimes,
    !// and LOSS-based lifetimes from ch4_loss3 and n2o_loss3.
    !//
    !// Ole Amund Sovde, October 2008, May 2012
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NRMETD, NROPSM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    real(r8), intent(in)  :: TAU
    logical, intent(in) :: zero
    !// Locals
    real(r8) :: RTMP
    character(len=3) :: cpres
    real(r8), parameter :: sec_yr = 365._r8*86400._r8
    !// Flag to print out more numbers for different domains.
    !// These are in general not so useful, perhaps even misleading.
    logical, parameter :: verbosePrint=.false.
    !// --------------------------------------------------------------------


    !// Accumulate first
    !// --------------------------------------------------------------------
    !// Tropospheric
    accum_oh_ch4_tp = accum_oh_ch4_tp &
         + sum(RKERN_CH4_TP_OH)/sum(RKERN_CH4_TP)
    accum_oh_co_tp = accum_oh_co_tp &
         + sum(RKERN_CO_TP_OH)/sum(RKERN_CO_TP)
    accum_ch4_tp = accum_ch4_tp &
         + sum(RTOT_CH4)/(sum(RLOSS_CH4_TP)*sec_yr)

    !// Below preslim
    accum_oh_ch4_hpa = accum_oh_ch4_hpa &
         + sum(RKERN_CH4_HPA_OH)/sum(RKERN_CH4_HPA)
    accum_oh_co_hpa = accum_oh_co_hpa &
         + sum(RKERN_CO_HPA_OH)/sum(RKERN_CO_HPA)
    accum_ch4_hpa = accum_ch4_hpa &
         + sum(RTOT_CH4)/(sum(RLOSS_CH4_HPA)*sec_yr)

    !// Spivakovsky
    accum_oh_spiv = accum_oh_spiv &
         + sum(RKERN_SPIV_OH)/sum(RKERN_SPIV)

    !// Lawrence
    RTMP = sum(RKERN_LAWR_OH) / sum(RKERN_LAWR) ! Avg OH
    accum_oh_lawr = accum_oh_lawr + RTMP
    accum_oh_lawr_air = accum_oh_lawr_air &
         + sum(RKERN_LAWR_AIR_OH)/sum(RKERN_LAWR_AIR)
    accum_ch4_lawr = accum_ch4_lawr &
         + sum(RKERN_LAWR_M)/(RTMP * sum(RKERN_LAWR_kM) * sec_yr)


    !// Global
    accum_oh_ch4_ltop = accum_oh_ch4_ltop &
         + sum(RKERN_CH4_LTOP_OH)/sum(RKERN_CH4_LTOP)
    accum_oh_co_ltop = accum_oh_co_ltop &
         + sum(RKERN_CO_LTOP_OH)/sum(RKERN_CO_LTOP)
    accum_ch4_ltop = accum_ch4_ltop &
         + sum(RTOT_CH4)/(sum(RLOSS_CH4_LTOP)*sec_yr)

    !// Increase number of accumulations
    num_accum = num_accum + 1._r8

    write(cpres(1:3),'(i3.3)') nint(preslim)

    !// OH from CH4 kernel
    !// ------------------------------------------------------------------
    write(6,'(A,f5.2,a4)') &
         '* Snapshot and average OH [molec/cm^3] at ',tau,' UTC'
    if (verbosePrint) then
       !// Up to LMTROP
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH (CH4) SFC-LMTROP        = ', &
            sum(RKERN_CH4_TP_OH)/sum(RKERN_CH4_TP),NDAY,NMET,NOPS, &
            accum_oh_ch4_tp/num_accum
       !// Up to xxxhPa
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH (CH4) SFC-'//cpres//'hPa        = ', &
            sum(RKERN_CH4_HPA_OH)/sum(RKERN_CH4_HPA),NDAY,NMET,NOPS, &
            accum_oh_ch4_hpa/num_accum
    end if
    !// Up to LTOP
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH (CH4) SFC-LTOP          = ', &
         sum(RKERN_CH4_LTOP_OH)/sum(RKERN_CH4_LTOP),NDAY,NMET,NOPS, &
         accum_oh_ch4_ltop/num_accum


    !// OH from CO kernel
    !// ------------------------------------------------------------------
    if (verbosePrint) then
       !// Up to LMTROP
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH (CO) SFC-LMTROP         = ', &
            sum(RKERN_CO_TP_OH)/sum(RKERN_CO_TP),NDAY,NMET,NOPS, &
            accum_oh_co_tp/num_accum
       !// Up to xxxhPa
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH (CO) SFC-'//cpres//'hPa         = ', &
            sum(RKERN_CO_HPA_OH)/sum(RKERN_CO_HPA),NDAY,NMET,NOPS, &
            accum_oh_co_hpa/num_accum
    end if
    !// Up to LTOP
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH (CO) SFC-LTOP           = ', &
         sum(RKERN_CO_LTOP_OH)/sum(RKERN_CO_LTOP),NDAY,NMET,NOPS, &
         accum_oh_co_ltop/num_accum

    !// OH from Spivakovsky kernel
    !// ------------------------------------------------------------------
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH (Spivakovsky)           = ', &
         sum(RKERN_SPIV_OH)/sum(RKERN_SPIV),NDAY,NMET,NOPS, &
         accum_oh_spiv/num_accum

    !// OH from Lawrence kernel
    !// ------------------------------------------------------------------
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH (Lawrence, kCH4)        = ', &
         sum(RKERN_LAWR_OH)/sum(RKERN_LAWR),NDAY,NMET,NOPS, &
         accum_oh_lawr/num_accum
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH (Lawrence, airmass)     = ', &
         sum(RKERN_LAWR_AIR_OH)/sum(RKERN_LAWR_AIR),NDAY,NMET,NOPS, &
         accum_oh_lawr_air/num_accum


    !// CH4 OH-lifetime
    !// ------------------------------------------------------------------
    write(6,'(A,f5.2,a4)') &
         '* Snapshot and average OH-lifetime of CH4 [years] at ',tau,' UTC'
    if (verbosePrint) then
       !// Up to LMTROP
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH-lifetime CH4 SFC-LMTROP = ', &
            sum(RTOT_CH4)/(sum(RLOSS_CH4_TP)*sec_yr), &
            NDAY,NMET,NOPS,accum_ch4_tp/num_accum
       !// Up to xxxhPa
       write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
            '  OH-lifetime CH4 SFC-'//cpres//'hPa = ', &
            sum(RTOT_CH4)/(sum(RLOSS_CH4_HPA)*sec_yr), &
            NDAY,NMET,NOPS,accum_ch4_hpa/num_accum
    end if
    !// Up to LTOP
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH-lifetime CH4 SFC-LTOP   = ', &
         sum(RTOT_CH4)/(sum(RLOSS_CH4_LTOP)*sec_yr), &
         NDAY,NMET,NOPS,accum_ch4_ltop/num_accum
    !// Lawrence
    RTMP = sum(RKERN_LAWR_OH) / sum(RKERN_LAWR) ! Avg OH [molec/cm3]
    write(6,'(A,es12.6,1x,i4,2(1x,i2),1x,es9.3)') &
         '  OH-lifetime CH4 Lawrence   = ', &
         sum(RKERN_LAWR_M)/(RTMP * sum(RKERN_LAWR_kM) * sec_yr), &
         NDAY,NMET,NOPS,accum_ch4_lawr/num_accum


    !// reset
    !// ------------------------------------------------------------------
    RKERN_CH4_TP(:)     = 0._r8
    RKERN_CH4_TP_OH(:)  = 0._r8
    RKERN_CH4_LTOP(:)   = 0._r8
    RKERN_CH4_LTOP_OH(:)= 0._r8
    RKERN_CH4_HPA(:)    = 0._r8
    RKERN_CH4_HPA_OH(:) = 0._r8

    RKERN_CO_TP(:)     = 0._r8
    RKERN_CO_TP_OH(:)  = 0._r8
    RKERN_CO_LTOP(:)   = 0._r8
    RKERN_CO_LTOP_OH(:)= 0._r8
    RKERN_CO_HPA(:)    = 0._r8
    RKERN_CO_HPA_OH(:) = 0._r8

    RKERN_SPIV(:)      = 0._r8
    RKERN_SPIV_OH(:)   = 0._r8

    RKERN_LAWR(:)     = 0._r8
    RKERN_LAWR_OH(:)  = 0._r8
    RKERN_LAWR_AIR(:)    = 0._r8
    RKERN_LAWR_AIR_OH(:) = 0._r8
    RKERN_LAWR_M(:)   = 0._r8
    RKERN_LAWR_kM(:)  = 0._r8

    RLOSS_CH4_TP(:)   = 0._r8
    RLOSS_CH4_HPA(:)  = 0._r8
    RLOSS_CH4_LTOP(:) = 0._r8
    RTOT_CH4(:)   = 0._r8

    !// --------------------------------------------------------------------
  end subroutine report_burden_and_lifetime
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine report_ch4n2o(NDAY,NMET,NOPS,TAU,zero)
    !// --------------------------------------------------------------------
    !// Print out N2O and CH4 burdens and lifetimes, based on
    !// chemical loss from ch4_loss3 and n2o_loss3.
    !//
    !// Ole Amund Sovde, February 2016, May 2012, October 2008
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NRMETD, NROPSM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    real(r8), intent(in) :: TAU
    logical, intent(in) :: zero
    !// Locals
    real(r8) :: RTMP
    !// --------------------------------------------------------------------

    !// For new N2O/CH4 lifetimes
    !// ---------------------------------------------------------------------
    !// For accumulated mass, divide by num_ch4n2o.
    !// Accumulated loss convert to kg/yr, first by
    !//    kg/acc.time * 1/(num_ch4n2o * dtops) => kg/s
    !// and then multiply by 31536000 to get kg/yr
    !//   accloss / (num_ch4n2o * dtops) * 31536000
    !// Lifetime is then
    !//   (accmass / num_ch4n2o) / (accloss / (num_ch4n2o * dtops) * 31536000)
    !// num_ch4n2o cancels, and dtops = 86400/(nropsm*nrmetd), so we
    !// multiply accmass/accloss with:
    RTMP = (86400._r8/real(NROPSM*NRMETD, r8)) / 31536000._r8
    !// Skip printout of lifetimes; these can be calculated from
    !// accumulated mass and loss.
    !//
    !// RAVG values are running averages for each month; i.e. it
    !// is reset every month.
    !//
    write(6,'(a,4es13.5)') 'RAVG.CH4.mass: ', &
         sum(GMASS_CH4_MON(:,1)) / num_ch4n2o, &
         sum(GMASS_CH4_MON(:,2)) / num_ch4n2o, &
         sum(GMASS_CH4_MON(:,3)) / num_ch4n2o, &
         sum(GMASS_CH4_MON(:,4)) / num_ch4n2o
    write(6,'(a,4es13.5)') 'RAVG.CH4.loss: ', &
         sum(GLOSS_CH4_MON(:,1)) / (RTMP*num_ch4n2o), &
         sum(GLOSS_CH4_MON(:,2)) / (RTMP*num_ch4n2o), &
         sum(GLOSS_CH4_MON(:,3)) / (RTMP*num_ch4n2o), &
         sum(GLOSS_CH4_MON(:,4)) / (RTMP*num_ch4n2o)
    write(6,'(a,4es13.5)') 'RAVG.CH4.life: ', &
         sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,1))) * RTMP, &
         sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,2))) * RTMP, &
         sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,3))) * RTMP, &
         sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,4))) * RTMP

    !// Monthly average burden and accumulated loss
    if (zero) then
       write(6,'(a,f10.0)') 'MONTHLY AVG timestep', num_ch4n2o
       write(6,'(a,4es13.5)') 'MONTHLY N2O life AVG:', &
            sum(GMASS_N2O_MON(:,1)) / (sum(GLOSS_N2O_MON(:,1))) * RTMP, &
            sum(GMASS_N2O_MON(:,1)) / (sum(GLOSS_N2O_MON(:,2))) * RTMP, &
             -99._r8, &          !// it is infinity
             sum(GMASS_N2O_MON(:,1))/ (sum(GLOSS_N2O_MON(:,4))) * RTMP
       write(6,'(a,4es13.5)') 'MONTHLY CH4 life AVG:', &
            sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,1))) * RTMP, &
            sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,2))) * RTMP, &
            sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,3))) * RTMP, &
            sum(GMASS_CH4_MON(:,1)) / (sum(GLOSS_CH4_MON(:,4))) * RTMP

       write(6,'(a,4es13.5)') 'MONTHLY N2O mass AVG:', &
            sum(GMASS_N2O_MON(:,1)) / num_ch4n2o, &
            sum(GMASS_N2O_MON(:,2)) / num_ch4n2o, &
            sum(GMASS_N2O_MON(:,3)) / num_ch4n2o, &
            sum(GMASS_N2O_MON(:,4)) / num_ch4n2o
       write(6,'(a,4es13.5)') 'MONTHLY CH4 mass AVG:', &
            sum(GMASS_CH4_MON(:,1)) / num_ch4n2o, &
            sum(GMASS_CH4_MON(:,2)) / num_ch4n2o, &
            sum(GMASS_CH4_MON(:,3)) / num_ch4n2o, &
            sum(GMASS_CH4_MON(:,4)) / num_ch4n2o
       !// Monthly loss [kg/month], not scaled to kg/yr
       write(6,'(a,4es13.5)') 'MONTHLY CH4 loss    :', &
            sum(GLOSS_CH4_MON(:,1)), &
            sum(GLOSS_CH4_MON(:,2)), &
            sum(GLOSS_CH4_MON(:,3)), &
            sum(GLOSS_CH4_MON(:,4))
       !// Stratospheric lifetime, calculated as
       !//    1 / (1/total_lifetime - 1/trop_lifetime)
       write(6,'(a,3es13.5)') 'MONTHLY CH4 STRAT LIFE AVG:       ', &
            !// Above preslim
            1._r8 / (1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,1))) * RTMP ) &
                   - 1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,2))) * RTMP ) ), &
            !// Above TP
            1._r8 / (1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,1))) * RTMP ) &
                  - 1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,3))) * RTMP ) ), &
            !// Above O3lim
            1._r8 / (1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,1))) * RTMP ) &
                  - 1._r8/( sum(GMASS_CH4_MON(:,1)) &
                           / (sum(GLOSS_CH4_MON(:,4))) * RTMP ) )

       !// Accumulate running means
       GLOSS_CH4_RTOT(:,:) = GLOSS_CH4_RTOT(:,:) + GLOSS_CH4_MON(:,:)
       GMASS_CH4_RTOT(:,:) = GMASS_CH4_RTOT(:,:) + GMASS_CH4_MON(:,:)
       GLOSS_N2O_RTOT(:,:) = GLOSS_N2O_RTOT(:,:) + GLOSS_N2O_MON(:,:)
       GMASS_N2O_RTOT(:,:) = GMASS_N2O_RTOT(:,:) + GMASS_N2O_MON(:,:)
       num_ch4n2o_rtot = num_ch4n2o_rtot + num_ch4n2o

       !// Print out running means
       write(6,'(a,f10.0)') 'RUNNING AVG timestep', num_ch4n2o_rtot
       write(6,'(a,4es13.5)') 'RUNNING N2O life AVG:', &
            sum(GMASS_N2O_RTOT(:,1)) / (sum(GLOSS_N2O_RTOT(:,1))) * RTMP, &
            sum(GMASS_N2O_RTOT(:,1)) / (sum(GLOSS_N2O_RTOT(:,2))) * RTMP, &
            -99._r8, &          !// it is infinity
            sum(GMASS_N2O_RTOT(:,1)) / (sum(GLOSS_N2O_RTOT(:,4))) * RTMP
       write(6,'(a,4es13.5)') 'RUNNING CH4 life AVG:', &
            sum(GMASS_CH4_RTOT(:,1)) / (sum(GLOSS_CH4_RTOT(:,1))) * RTMP, &
            sum(GMASS_CH4_RTOT(:,1)) / (sum(GLOSS_CH4_RTOT(:,2))) * RTMP, &
            sum(GMASS_CH4_RTOT(:,1)) / (sum(GLOSS_CH4_RTOT(:,3))) * RTMP, &
            sum(GMASS_CH4_RTOT(:,1)) / (sum(GLOSS_CH4_RTOT(:,4))) * RTMP
       write(6,'(a,4es13.5)') 'RUNNING N2O mass AVG:', &
            sum(GMASS_N2O_RTOT(:,1)) / num_ch4n2o_rtot, &
            sum(GMASS_N2O_RTOT(:,2)) / num_ch4n2o_rtot, &
            sum(GMASS_N2O_RTOT(:,3)) / num_ch4n2o_rtot, &
            sum(GMASS_N2O_RTOT(:,4)) / num_ch4n2o_rtot

       write(6,'(a,4es13.5)') 'RUNNING CH4 mass AVG:', &
            sum(GMASS_CH4_RTOT(:,1)) / num_ch4n2o_rtot, &
            sum(GMASS_CH4_RTOT(:,2)) / num_ch4n2o_rtot, &
            sum(GMASS_CH4_RTOT(:,3)) / num_ch4n2o_rtot, &
            sum(GMASS_CH4_RTOT(:,4)) / num_ch4n2o_rtot
       write(6,'(a,4es13.5)') 'RUNNING CH4 loss    :', &
            sum(GLOSS_CH4_RTOT(:,1)), &
            sum(GLOSS_CH4_RTOT(:,2)), &
            sum(GLOSS_CH4_RTOT(:,3)), &
            sum(GLOSS_CH4_RTOT(:,4))

       write(6,'(a,3es13.5)') 'RUNNING AVG CH4 STRAT LIFE:       ', &
            !// Above preslim
            1._r8 / (1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,1))) * RTMP ) &
                  - 1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,2))) * RTMP ) ), &
            !// Above TP
            1._r8 / (1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,1))) * RTMP ) &
                  - 1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,3))) * RTMP ) ), &
            !// Above O3lim
            1._r8 / (1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,1))) * RTMP ) &
                  - 1._r8/( sum(GMASS_CH4_RTOT(:,1)) &
                           / (sum(GLOSS_CH4_RTOT(:,4))) * RTMP ) )

       !// Re-initialize
       call init_lifetime()
    end if

    !// --------------------------------------------------------------------
  end subroutine report_ch4n2o
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine diag_burden_snapshot(nday,nmet,nops)
    !// --------------------------------------------------------------------
    !// Diagnose trop/strat burdens.
    !// Only selected species are included, such as O3, but you
    !// can add your own following the same structure.
    !//
    !// Currently only carried out for NOPS=1.
    !//
    !// Ole Amund Sovde, November 2014, January 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: LOSLOCTROP
    use cmn_ctm, only: STT, AIR, ETAA, ETAB
    use cmn_chem, only: TNAME
    use cmn_met, only: P
    use cmn_oslo, only: LMTROP, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    integer, intent(in) :: nday,nmet,nops

    !// Locals
    integer :: I,J,L,K,cid,trnr
    real(r8) :: tot_mass_trp,tot_mass_str,tot_mass_trp2,tot_mass_str2, &
         tot_mass_trp3,tot_mass_str3,pfrac,P1,P2
    character(len=3) :: cpres

    integer, parameter :: nc=1
    integer, parameter, dimension(nc) :: comps = (/ 1 /)
    !// --------------------------------------------------------------------

    !// Only calculate each NMET
    if (nops .ne. 1) return

    !// But in case output should be done every NOPS, update
    !// tropopause data every NMET
    if (nops .eq. 1) then
       !// Find tp_o3 and tp_hpa each NMET
       if (trsp_idx(1) .gt. 0) then
         tp_o3(:,:) = LMTROP(:,:) !// Initialize
         do J = 1, JPAR
           do I = 1, IPAR
             do L = 10, LPAR
                if (STT(I,J,L,trsp_idx(1))/AIR(I,J,L) .gt. o3lim) then
                   tp_o3(I,J) = L-1
                   exit
                end if
             end do
           end do
         end do
       else
         tp_o3(:,:) = LMTROP(:,:)
       end if
       tp_hpa(:,:) = 1 !// Initialize
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                if ( (etaa(L+1) + p(i,j)*etab(L+1)) .lt. preslim) exit
                tp_hpa(I,J) = L
             end do
          end do
       end do
    end if

    !// Print header if instantaneous masses are to be put out
    do K = 1, nc
       CID = comps(K)
       TRNR = trsp_idx(CID)
       if (TRNR .gt. 0) then
          write(cpres(1:3),'(i3.3)') nint(preslim)
          write(6,'(a,3f15.5)')'* Budgets below/above '// &
               'LMTROP/150ppbv(O3)/'//cpres//'hPa'
          exit !// Exit loop to start putting out values
       end if
    end do

    !// Put out tracers
    do K = 1, nc
       CID = comps(K)
       TRNR = trsp_idx(CID)
       if (TRNR .gt. 0) then
          tot_mass_trp  = 0._r8
          tot_mass_str  = 0._r8
          tot_mass_trp2 = 0._r8
          tot_mass_str2 = 0._r8
          tot_mass_trp3 = 0._r8
          tot_mass_str3 = 0._r8
          do L = 1, LPAR
           do J = 1, JPAR
            do I = 1, IPAR
             if (L .le. tp_o3(I,J)) then
               tot_mass_trp2 = tot_mass_trp2 + STT(I,J,L,TRNR)
             else
               tot_mass_str2 = tot_mass_str2 + STT(I,J,L,TRNR)
             end if
             if (L .le. LMTROP(I,J)) then
               tot_mass_trp = tot_mass_trp + STT(I,J,L,TRNR)
             else
               tot_mass_str = tot_mass_str + STT(I,J,L,TRNR)
             end if
             if (L .le. tp_hpa(I,J)) then
               tot_mass_trp3 = tot_mass_trp3 + STT(I,J,L,TRNR)
             else if (L .eq. tp_hpa(I,J)+1) then
               !// Split: P1 below preslim, P2 above
               P1 = etaa(L) + p(i,j) * etab(L)
               P2 = etaa(L+1) + p(i,j) * etab(L+1)
               PFRAC = (P1 - preslim) / (P1 - P2)
               tot_mass_trp3 = tot_mass_trp3 + STT(I,J,L,TRNR) * pfrac
               tot_mass_str3 = tot_mass_str3 + STT(I,J,L,TRNR) *(1._r8-pfrac)
             else
               tot_mass_str3 = tot_mass_str3 + STT(I,J,L,TRNR)
             end if
            end do
           end do
          end do
          write(6,'(a,3f13.5)')'  Trop '//trim(TNAME(TRNR))//'[Tg]: ', &
          tot_mass_trp*1.e-9_r8, tot_mass_trp2*1.e-9_r8, tot_mass_trp3*1.e-9_r8
          write(6,'(a,3f13.5)')'  Strt '//trim(TNAME(TRNR))//'[Tg]: ', &
          tot_mass_str*1.e-9_r8, tot_mass_str2*1.e-9_r8, tot_mass_str3*1.e-9_r8
       end if
    end do !// do K = 1, nc

    !// --------------------------------------------------------------------
  end subroutine diag_burden_snapshot
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine write_snapshot(nday,nmet,nops,mode)
    !// --------------------------------------------------------------------
    !// Writes out snapshots of selected species/metdata for
    !// northern hemisphere only.
    !// Called in pmain at the beginning of each NOPS.
    !//
    !// Ole Amund Sovde, June 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: NROPSM, NRMETD, ETAA, ETAB, XDGRD, YDGRD, &
         XDEDG, YDEDG, STT, AIR, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         JYEAR, JMON, JDATE
    use cmn_met, only: P, T, PVU
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx
    use strat_h2o, only: LOLD_H2OTREATMENT, str_h2o
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NMET, NOPS, MODE
    !// Locals
    character(len=80) ::  daily_file
    character(len=8)  ::  datestamp
    integer :: ifnr, JSTART, JEND, MP, II, JJ, I,J,L, N, NHOUR
    real(r4) :: dumr4(ipar,jpar,lpar)
    integer, parameter :: NsnapComps = 2
    integer, dimension(NsnapComps), parameter :: &
         snapComps = (/1, 107/)
    !// --------------------------------------------------------------------

    !// Only put out every hour, check NOPS vs NROPSM
    NHOUR = nint(real(NROPSM)/(real(24)/real(NRMETD)))
    if ( mod(NOPS-1, NHOUR) .ne. 0) return

    !// Which region to snapshot is controlled by MODE
    if (MODE .eq. 1) then
       !// Northern hemisphere
       JSTART = JPAR/2 + 1
       JEND   = JPAR
    else if (MODE .eq. 2) then
       !// Southern hemisphere
       JSTART = 1
       JEND   = JPAR/2
    else
       !// GLOBAL
       JSTART = 1
       JEND   = JPAR
    end if

    !// Get file number
    ifnr = get_free_fileid()

    !// File name
    write(datestamp(1:8),'(i4.4,2i2.2)') jyear,jmon,jdate
    daily_file='snapshots_'//datestamp//'.dta'

    if (nmet.eq.1 .and. nops.eq.1) then
       !// Start new file at the beginning of the day
       open(ifnr, file=daily_file, form='unformatted', status='unknown')
       write(ifnr) ipar,jpar,lpar,jstart,jend, NRMETD, NROPSM, NsnapComps
       write(ifnr) etaa, etab
       write(ifnr) xdgrd, xdedg
       write(ifnr) ydgrd, ydedg
    else
       !// Open for append on the rest of the day
       open(ifnr,file=daily_file,form='unformatted',&
            status='unknown',access='append')
    end if
    !// Write
    write(ifnr) NDAY,NMET,NOPS
    do N = 1, NsnapComps
       write(ifnr) snapComps(N)
       if (snapComps(N) .eq. 114 .and. LOLD_H2OTREATMENT) then
          !// Get str_h2o
          do MP = 1, MPBLK
             do J = MPBLKJB(MP), MPBLKJE(MP)
                JJ    = J - MPBLKJB(MP) + 1
                do I = MPBLKIB(MP), MPBLKIE(MP)
                   II    = I - MPBLKIB(MP) + 1
                   do L = 1, LPAR
                      dumr4(i,j,l) = str_h2o(L,II,JJ,MP)
                   end do
                end do
             end do
          end do
          write(ifnr) dumr4(:,JSTART:JEND,:)
       else
          write(ifnr) real(STT(:,JSTART:JEND,:,trsp_idx(snapComps(N))), r4)
       end if
    end do !// do N = 1, NsnapComps
    write(ifnr) real(AIR(:,JSTART:JEND,:), r4)
    write(ifnr) real(pvu(:,:,JSTART:JEND), r4)
    write(ifnr) real(T(:,JSTART:JEND,:), r4)
    write(ifnr) real(p(:,JSTART:JEND), r4)
    close(ifnr)
    !// --------------------------------------------------------------------
  end subroutine write_snapshot
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine write_snapshot_the(nday,nmet,nops,mode)
    !// --------------------------------------------------------------------
    !// Writes out number mixing ratio snapshots of selected species/metdata
    !// on theta levels defined in physics_oslo.f90.
    !//
    !// Called from subroutine nops_diag in this module.
    !//
    !// Ole Amund Sovde, February 2016
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: NROPSM, NRMETD, ETAA, ETAB, XDGRD, YDGRD, &
         XDEDG, YDEDG, STT, AIR, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         JYEAR, JMON, JDATE
    use cmn_chem, only: TMASSMIX2MOLMIX
    use cmn_met, only: P, T, PVU
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx
    use physics_oslo, only: NTHE, pvthe, pvtheta, theqlat, IJLfield2ThetaLvs
    use strat_h2o, only: LOLD_H2OTREATMENT, str_h2o
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NMET, NOPS, MODE
    !// Locals
    character(len=80) ::  daily_file
    character(len=8)  ::  datestamp
    integer :: ifnr, JSTART, JEND, MP, II, JJ, I,J,L, N, NHOUR, TRNR
    real(r8) :: r8tmp(ipar,jpar,lpar)
    !// Define components to put out
    integer, parameter :: NsnapComps = 2
    integer, dimension(NsnapComps), parameter :: &
         snapComps = (/1, 107/)
    real(r8) :: the_stt(IPAR,JPAR,NTHE)
    !// --------------------------------------------------------------------

    !// Only put out every hour, check NOPS vs NROPSM
    NHOUR = nint(real(NROPSM)/(real(24)/real(NRMETD)))
    if ( mod(NOPS-1, NHOUR) .ne. 0) return

    !// Which region to snapshot is controlled by MODE
    if (MODE .eq. 1) then
       !// Northern hemisphere
       JSTART = JPAR/2 + 1
       JEND   = JPAR
    else if (MODE .eq. 2) then
       !// Southern hemisphere
       JSTART = 1
       JEND   = JPAR/2
    else
       !// GLOBAL
       JSTART = 1
       JEND   = JPAR
    end if

    !// Get file number
    ifnr = get_free_fileid()


    !// File name
    write(datestamp(1:8),'(i4.4,2i2.2)') jyear,jmon,jdate
    daily_file='theta_snapshots_'//datestamp//'.dta'

    if (nmet.eq.1 .and. nops.eq.1) then
       !// Start new file at the beginning of the day
       open(ifnr, file=daily_file, form='unformatted', status='unknown')
       write(ifnr) ipar,jpar,lpar,jstart,jend, NTHE, NRMETD, NROPSM, NsnapComps
       write(ifnr) etaa, etab
       write(ifnr) xdgrd, xdedg
       write(ifnr) ydgrd, ydedg
       write(ifnr) pvthe
    else
       !// Open for append on the rest of the day
       open(ifnr,file=daily_file,form='unformatted',&
            status='unknown',access='append')
    end if

    !// Write time info
    write(ifnr) NDAY,NMET,NOPS

    !// For each species find values on theta levels
    do N = 1, NsnapComps
       TRNR = trsp_idx(snapComps(N)) !// transport number

       write(ifnr) snapComps(N)
       if (snapComps(N) .eq. 114 .and. LOLD_H2OTREATMENT) then
          !// Get str_h2o
          do MP = 1, MPBLK
             do J = MPBLKJB(MP), MPBLKJE(MP)
                JJ    = J - MPBLKJB(MP) + 1
                do I = MPBLKIB(MP), MPBLKIE(MP)
                   II    = I - MPBLKIB(MP) + 1
                   do L = 1, LPAR
                      r8tmp(i,j,l) = str_h2o(L,II,JJ,MP)
                   end do
                end do
             end do
          end do
          !// Special treatment for quasi-static H2O
          r8tmp(:,:,:) = r8tmp(:,:,:)/AIR(:,:,:) * M_AIR / 18._r8
       else
          !// Convert tracer to number mixing ratio
          r8tmp(:,:,:) = STT(:,:,:,TRNR) / AIR(:,:,:) * TMASSMIX2MOLMIX(TRNR)
       end if
       !// Interpolate to theta levels
       call IJLfield2ThetaLvs(r8tmp, the_stt, JSTART, JEND)
       !// Write tracer
       write(ifnr) real(the_stt(:,JSTART:JEND,:), r4)
    end do !// do N = 1, NsnapComps

    !// Write PVU and equivalent latitude on theta levels
    write(ifnr) real(pvtheta(:,JSTART:JEND,:), r4)
    write(ifnr) real(theqlat(:,JSTART:JEND,:), r4)
    close(ifnr)
    !// --------------------------------------------------------------------
  end subroutine write_snapshot_the
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4n2o_burden(BTT,ADDMASS,DT,MP)
    !// --------------------------------------------------------------------
    !// Find CH4 and N2O burden.
    !// Note that routine is called NCYCLES times during time step. We
    !// sum up mass (kg) and loss (kg/s) equal times, so a lifetime can
    !// be calculated as the ratio of the two.
    !//
    !// Ole Amund Sovde, October 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: LOSLOCTROP, LOSLOCSTRAT
    use cmn_ctm, only: ETAA, ETAB, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: P
    use cmn_oslo, only: trsp_idx, LMTROP, METHANEMIS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT
    real(r8), intent(in) :: DT
    logical, intent(in) :: ADDMASS

    !// Local Parameters
    real(r8) :: RLOSS(4), PRES(LPAR), frac
    integer :: I,J,L,II,JJ
    !// --------------------------------------------------------------------

    !// If not tropospheric chemistry, return
    if (.not. LOSLOCTROP) return

    !// If not time for adding mass, return
    if (.not. ADDMASS) return

    !// Find CH4 first, then N2O

    !// Find accumulated CH4 before chemistry
    GLINS_CH4(MP,:) = 0._r8
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Whole domain
        do L = 1, LPAR-1
           GLINS_CH4(MP,1) = GLINS_CH4(MP,1) + BTT(L,trsp_idx(46),II,JJ)
        end do

        !// Upper edge of boxes
        PRES(:) = ETAA(2:LPAR+1) + P(I,J)*ETAB(2:LPAR+1)
        !// Up to preslim
        do L = 1, LPAR-1
           !// Exit when box top is above preslim, interpolate below
           if (PRES(L) .lt. preslim) exit
           GLINS_CH4(MP,2) = GLINS_CH4(MP,2) + BTT(L,trsp_idx(46),II,JJ)
        end do
        !// Add fraction between PRES(L-1) and preslim
        if (L .le. LPAR) then
           frac = (PRES(L-1)-preslim)/(PRES(L-1)-PRES(L))
           GLINS_CH4(MP,2) = GLINS_CH4(MP,2) + BTT(L,trsp_idx(46),II,JJ) * frac
        end if

        !// Up to TP
        do L = 1, LMTROP(I,J)
           GLINS_CH4(MP,3) = GLINS_CH4(MP,3) + BTT(L,trsp_idx(46),II,JJ)
        end do

        !// Up to 150ppbv O3
        do L = 1, TP_O3(I,J)
           GLINS_CH4(MP,4) = GLINS_CH4(MP,4) + BTT(L,trsp_idx(46),II,JJ)
        end do

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// Save the total amount of CH4
    GMASS_CH4_MON(MP,:) = GMASS_CH4_MON(MP,:) + GLINS_CH4(MP,:)

    !// Add up number of accumulations
    if (MP .eq. 1) num_ch4n2o = num_ch4n2o + 1._r8

    !// Find N2O if stratchem is included
    if (LOSLOCSTRAT) then
       !// Find accumulated CH4 before chemistry
       GLINS_N2O(MP,:) = 0._r8
       !// Loop over latitude (J is global, JJ is block)
       do J = MPBLKJB(MP),MPBLKJE(MP)
         JJ    = J - MPBLKJB(MP) + 1

         !// Loop over longitude (I is global, II is block)
         do I = MPBLKIB(MP),MPBLKIE(MP)
           II    = I - MPBLKIB(MP) + 1

           !// Whole domain
           do L = 1, LPAR-1
              GLINS_N2O(MP,1) = GLINS_N2O(MP,1) + BTT(L,trsp_idx(107),II,JJ)
           end do

           !// Upper edge of boxes
           PRES(:) = ETAA(2:LPAR+1) + P(I,J)*ETAB(2:LPAR+1)
           !// Up to preslim
           do L = 1, LPAR-1
              !// Exit when box top is above preslim, interpolate below
              if (PRES(L) .lt. preslim) exit
              GLINS_N2O(MP,2) = GLINS_N2O(MP,2) + BTT(L,trsp_idx(107),II,JJ)
           end do
           !// Add fraction between PRES(L-1) and preslim
           if (L .le. LPAR) then
              frac = (PRES(L-1)-preslim)/(PRES(L-1)-PRES(L))
              GLINS_N2O(MP,2) = GLINS_N2O(MP,2) + BTT(L,trsp_idx(107),II,JJ) * frac
           end if

           !// Up to TP
           do L = 1, LMTROP(I,J)
              GLINS_N2O(MP,3) = GLINS_N2O(MP,3) + BTT(L,trsp_idx(107),II,JJ)
           end do

           !// Up to 150ppbv O3
           do L = 1, TP_O3(I,J)
              GLINS_N2O(MP,4) = GLINS_N2O(MP,4) + BTT(L,trsp_idx(107),II,JJ)
           end do
         end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

       !// Save the total amount of N2O
       GMASS_N2O_MON(MP,:) = GMASS_N2O_MON(MP,:) + GLINS_N2O(MP,:)

    end if !// if (LOSLOCSTRAT) then

    !// --------------------------------------------------------------------
  end subroutine ch4n2o_burden
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4_loss3(TRNR,CHEMLOSS,DV,PUEDG, I,J,II,JJ,MP,L1,L2,LMT)
    !// --------------------------------------------------------------------
    !// Accumulate CH4 lost in chemistry, convert to units [kg].
    !//
    !// Ole Amund Sovde, February 2015
    !// --------------------------------------------------------------------
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: TRNR,I,J,II,JJ,MP,L1,L2,LMT
    real(r8), intent(in) :: CHEMLOSS(LPAR)
    real(r8), intent(in) :: DV(LPAR)
    real(r8), intent(in) :: PUEDG(LPAR)

    !// Locals
    real(r8)  :: rfac, lostmass, RLOSS(4)
    integer :: L
    !// --------------------------------------------------------------------

    rfac = 1.e6_r8 / AVOGNR * TMASS(TRNR) * 1.e-3_r8

    RLOSS(:) = 0._r8 !// Loss diag for 4 vertical domains
    do L = L1, L2

       !// Total loss in entry 1 [molec/cm3, accumulated for this step]
       !// molecules/cm3 > kg/gridbox
       lostmass = CHEMLOSS(L) * DV(L) * rfac

       RLOSS(1) = RLOSS(1) + lostmass

       !// Up to preslim
       if (PUEDG(L) .ge. preslim) RLOSS(2) = RLOSS(2) + lostmass
       !// Add fraction between PRES(L-1) and preslim
       if (PUEDG(L) .lt. preslim) then
          if (PUEDG(L-1).gt.preslim) RLOSS(2) = RLOSS(2) &
               + lostmass * (PUEDG(L-1)-preslim)/(PUEDG(L-1)-PUEDG(L))
       end if

       !// Up to TP
       if (L .le. LMT) RLOSS(3) = RLOSS(3) + lostmass

       !// Up to 150ppbv O3
       if (L .le. TP_O3(I,J)) RLOSS(4) = RLOSS(4) + lostmass

    end do

    !// Save to be put out each avgsav file
    GLOSS_CH4_MON(MP,:) = GLOSS_CH4_MON(MP,:) + RLOSS(:)

    !// --------------------------------------------------------------------
  end subroutine ch4_loss3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine n2o_loss3(TRNR,CHEMLOSS,DV,PUEDG, I,J,II,JJ,MP,L1,L2,LMT)
    !// --------------------------------------------------------------------
    !// Accumulate N2O lost in chemistry, convert to units [kg].
    !//
    !// Ole Amund Sovde, February 2015
    !// --------------------------------------------------------------------
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: TRNR,I,J,II,JJ,MP,L1,L2,LMT
    real(r8), intent(in) :: CHEMLOSS(LPAR), DV(LPAR), PUEDG(LPAR)

    !// Locals
    real(r8)  :: rfac, lostmass, RLOSS(4)
    integer :: L
    !// --------------------------------------------------------------------

    rfac = 1.e6_r8 / AVOGNR * TMASS(TRNR) * 1.e-3_r8

    RLOSS(:) = 0._r8 !// Loss diag for 4 vertical domains
    do L = L1, L2

       !// Total loss in entry 1 [molec/cm3, accumulated for this step]
       !// molecules/cm3 > kg/gridbox
       lostmass = CHEMLOSS(L) * DV(L) * rfac
       RLOSS(1) = RLOSS(1) + lostmass

       !// Up to preslim
       if (PUEDG(L) .ge. preslim) RLOSS(2) = RLOSS(2) + lostmass
       !// Add fraction between PRES(L-1) and preslim
       if (PUEDG(L) .lt. preslim) then
          if (PUEDG(L-1).gt.preslim) RLOSS(2) = RLOSS(2) &
               + lostmass * (PUEDG(L-1)-preslim)/(PUEDG(L-1)-PUEDG(L))
       end if

       !// Up to TP
       if (L .le. LMT) RLOSS(3) = RLOSS(3) + lostmass

       !// Up to 150ppbv O3
       if (L .le. TP_O3(I,J)) RLOSS(4) = RLOSS(4) + lostmass

    end do !// do L = L1, L2

    !// Save to be put out each avgsav file
    GLOSS_N2O_MON(MP,:) = GLOSS_N2O_MON(MP,:) + RLOSS(:)

    !// --------------------------------------------------------------------
  end subroutine n2o_loss3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine diags_tpset()
    !// --------------------------------------------------------------------
    !// Set levels of 150ppbv O3 and xxxhPa.
    !//
    !// Ole Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    use cmn_ctm, only: ETAA, ETAB, STT, AIR
    use cmn_met, only: P
    use cmn_oslo, only: LMTROP, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: i,j,l, LSTART
    !// --------------------------------------------------------------------

    !// Find tp_o3 and tp_hpa
    LSTART = 1
    do L = 1, LPAR
       if ((etaa(L+1) + minval(p)*etab(L+1)) .lt. (600._r8)) exit
       !// Start at the level where top pressure is > 600hPa.
       LSTART = L
    end do
    if (trsp_idx(1) .gt. 0) then
       tp_o3(:,:) = LSTART !LMTROP(:,:)  !// Initialize
       do J = 1, JPAR
          do I = 1, IPAR
             do L = LSTART, LPAR
                if (STT(I,J,L,trsp_idx(1))/AIR(I,J,L) .gt. o3lim) then
                   tp_o3(I,J) = L-1
                   exit
                end if
             end do
          end do
       end do
    else
       tp_o3(:,:) = 1
    end if

    !// TP defined using pressure limit preslim
    LSTART = 1
    do L = 1, LPAR
       if ((etaa(L+1) + minval(p)*etab(L+1)) .lt. (preslim+100.)) exit
       LSTART = L
    end do
    tp_hpa(:,:) = LSTART
    do J = 1, JPAR
       do I = 1, IPAR
          do L = LSTART, LPAR
             if ( (etaa(L+1) + p(i,j)*etab(L+1)) .lt. preslim) exit
             tp_hpa(I,J) = L
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine diags_tpset
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine report_negO3()
    !// --------------------------------------------------------------------
    !// Prints out number of negative O3 reported by OSLOCHEM.
    !// This is accumulated over the sub-stepping loop NSUB.
    !//
    !// Ole Amund Sovde, December 2013
    !// --------------------------------------------------------------------
    use cmn_oslo, only: TROPCHEMnegO3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8) :: totNEGS
    integer :: L
    !// --------------------------------------------------------------------

    !// Number of negatives
    totNEGS = sum(TROPCHEMnegO3)
    if (totNEGS .eq. 0) return

    !// Print out for each level (probably only L1
    do L = 1, LPAR
       totNEGS = sum(TROPCHEMnegO3(L,:))
       if (totNEGS .eq. 0) cycle
       write(6,'(a,i2,a3,i5)') &
            '* Correcting tropchem NEGATIVE O3 [LEV/NUMBER]: ',L, &
            ' / ',nint(totNEGS)
    end do
    TROPCHEMnegO3(:,:) = 0._r8

    !// --------------------------------------------------------------------
  end subroutine report_negO3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tnd_emis_daily(DDAY)
    !// --------------------------------------------------------------------
    !// Saves accumulated emissions. Applies only when emissions are
    !// treated in chemistry.
    !// Called from daily_diag_output.
    !//
    !// Ole Amund Sovde, September 2014
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TNAME
    use cmn_oslo, only: chem_idx, DIAGEMIS_IJ, emisTotalsDaily, emisTotalsOld
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: DDAY
    !// Locals
    integer :: I,J,II,JJ,N,L,MP
    real(r8) :: DTDIAG, ZTG, totemis(NPAR)
    !// --------------------------------------------------------------------
    !// Sum up total emissions accumulated so far
    totemis(:) = 0._r8
    do MP = 1, MPBLK
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II   = I - MPBLKIB(MP) + 1
             do N = 1, NPAR
                do L = 1, LPAR
                   totemis(N) = totemis(N) + DIAGEMIS_IJ(L,N,II,JJ,MP)
                end do
             end do
          end do
       end do
    end do
    !// Diagnose emissions of this day
    emisTotalsDaily(:,DDAY) = totemis(:) - emisTotalsOld(:)
    !// Update the values to be subtracted next day. Note that
    !// when DIAGEMIS_IJ is set to zero in TBGT_G (p-bdgts.f),
    !// emisTotalsOld also has to be set to zero.
    emisTotalsOld(:) = totemis(:)
    !// --------------------------------------------------------------------
  end subroutine tnd_emis_daily
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tnd_emis2file(NDAY)
    !// --------------------------------------------------------------------
    !// Write daily emission totals. This is written to file, not
    !// to std.out, in diagnostic/budget print-out section of p-main.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Output as netCDF4.
    !// Ole Amund Sovde, September 2014
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_size, only: TNAMELEN
    use cmn_ctm, only: ETAA, ETAB, AREAXY, MPBLKIB, MPBLKIE, MPBLKJB, &
         MPBLKJE, JYEAR, JMON, JDATE, XDGRD, YDGRD, XDEDG, YDEDG, &
         ZGRD, ZEDG
    use cmn_chem, only: TNAME, TMASS
    use cmn_diag, only: NDAY0,JYEAR0,JMON0,JDATE0, &
         nc4deflate_global, nc4shuffle_global
    use utilities, only: get_free_fileid
    use cmn_oslo, only: chem_idx, emisTotalsDaily, DIAGEMIS_IJ
    use netcdf
    use ncutils, only: handle_error
    use emissions_megan, only: megan_emis2file
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    !// Locals
    integer :: fnr
    character(len=80) :: filename
    character(len=4) :: CYEAR
    character(len=8) :: datestamp, datestamp0
    integer :: I,J,II,JJ,N,L,MP, LMAX, Ncount
    real(r8) :: DTDIAG, ZTG, totemis(NPAR), emis3d(IPAR,JPAR,LPAR)
    integer, dimension(6) :: start_time, end_time
    character(len=TNAMELEN), dimension(NPAR) :: tracer_name

    integer :: &
         lat_dim_id, lon_dim_id, lev_dim_id, time_dim_id, &
         lat_id, lon_id, lev_id, time_id, &
         ilat_dim_id, ilon_dim_id, ilev_dim_id, &
         ilat_id, ilon_id, ilev_id, &
         ihya_dim_id, ihyb_dim_id, &
         ihya_id, ihyb_id, &
         areaxy_id, &
         date_size_dim_id, &
         start_time_id, end_time_id, & !IDs for start/end dates for average
         start_day_id, end_day_id, &   !IDs for start/end day (NDAY)
         tracer_name_len_dim_id, tracer_name_id, &
         emcomp_dim_id, chemid_id, molw_id, &
         status, ncid
    integer, dimension(NPAR) :: comps_id

    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tnd_emis2file'
    !// --------------------------------------------------------------------

    !// Get file number
    fnr = get_free_fileid()

    !// YEAR may have been updated by CALENDR before calling
    !// so we check if the new day is 1 Jan
    if (JMON.eq.1 .and. JDATE.eq.1) then
       write(CYEAR(1:4),'(i4.4)') JYEAR - 1
    else
       write(CYEAR(1:4),'(i4.4)') JYEAR
    end if

    !// Write daily total emissions
    !// ------------------------------------------------------------
    !// Open file for this year
!    filename = 'emis_daily_totals_'//CYEAR//'.dta'
!    open(fnr,file=filename,form='unformatted')
!    write(fnr) NPAR, JYEAR
!    !// Write tracer info
!    write(fnr) TNAME
!    write(fnr) chem_idx
!    write(fnr) TMASS
!    !// Write total emissions
!    write(fnr) emisTotalsDaily
!    !// Close file
!    close(fnr)

    !// Also print total accumulated emissions to screen
    !// ------------------------------------------------------------
    totemis(:) = 0._r8
    do MP = 1, MPBLK
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II   = I - MPBLKIB(MP) + 1
             do N = 1, NPAR
                do L = 1, LPAR
                   totemis(N) = totemis(N) + DIAGEMIS_IJ(L,N,II,JJ,MP)
                end do
             end do
          end do
       end do
    end do

    !// Print out totals for the diagnostic period
    DTDIAG = real(NDAY+1 - NDAY0, r8)
    ZTG = 365._r8/(DTDIAG) * 1.e-9_r8
    write(6,'(a)') '* Total emissions this period '// &
                     '[Tg/year] + DT(diagnose) [days]:'

    !// Find number of components emitted, retrieve names and
    !// print totals to screen.
    Ncount = 0
    do N = 1, NPAR
       if (totemis(N) .gt. 0._r8) then
          Ncount = Ncount + 1
          write(6,'(a,1x,a10,i4,i4,es20.10,f6.0)') 'Emission:', &
               TNAME(N),chem_idx(N),N,totemis(N)*ZTG,DTDIAG
          !// Save tracer name for emitted species
          tracer_name(Ncount) = TNAME(N)
       end if
    end do

    !//---------------------------------------------------------------------
    !// Write daily total emissions - netCDF4
    !//---------------------------------------------------------------------
    !// Open file for this year
    filename = 'emis_daily_totals_'//CYEAR//'.nc'
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')

    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title', &
         'Daily total emissions in Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    !//---------------------------------------------------------------------

    !// Define spatial dimensions (npar)
    status = nf90_def_dim(ncid,"NPAR",NPAR,emcomp_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NPAR dim')
    !// Define spatial dimensions (time)
    status = nf90_def_dim(ncid,"time",366,time_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define time dim')

    !// Define length of tracer name string
    status = nf90_def_dim(ncid,"tracer_name_len", &
         TNAMELEN, tracer_name_len_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')

    !// Year of data
    status = nf90_def_var(ncid,"year",nf90_int,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define year variable')

    !// Tracer names
    status = nf90_def_var(ncid,"tracer_name",nf90_char, &
         (/tracer_name_len_dim_id,emcomp_dim_id/),tracer_name_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Tracer names')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Component IDs
    status = nf90_def_var(ncid,"chem_idx",nf90_int, &
         (/emcomp_dim_id/),chemid_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define chem_idx variable')
    status = nf90_put_att(ncid,chemid_id,'description', &
         'Component ID in Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description chem_idx')

    !// Tracer mass (molecular weight)
    status = nf90_def_var(ncid,"tracer_molweight",nf90_double, &
         (/emcomp_dim_id/),molw_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_molweight variable')
    status = nf90_put_att(ncid,molw_id,'description', &
         'Tracer mlecular weights')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_molweight')
    status = nf90_put_att(ncid,molw_id,'units','g/mol')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_molweight')

    !// Define emissions variable
    status = nf90_def_var(ncid,"emisTotalsDaily", nf90_double, &
         (/emcomp_dim_id, time_dim_id/), comps_id(1))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define emisTotalsDaily variable')
    status = nf90_def_var_deflate(ncid,comps_id(1),nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate emisTotalsDaily variable')
    status = nf90_put_att(ncid,comps_id(1),'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit emisTotalsDaily')

    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Tracer daily emission totals
    status = nf90_put_var(ncid,comps_id(1),emisTotalsDaily)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting emisTotalsDaily')

    !// Year of data
    status = nf90_put_var(ncid,time_id,JYEAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting year')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,TNAME)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Component ID
    status = nf90_put_var(ncid,chemid_id,chem_idx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting chem_idx')

    !// Tracer molecular weights
    status = nf90_put_var(ncid,molw_id,TMASS)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_molweight')

    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------


    !//---------------------------------------------------------------------


    !//---------------------------------------------------------------------
    !// Write 3D emission data to file - netCDF4
    !//---------------------------------------------------------------------
    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
    write(datestamp0(1:8),'(i4.4,2i2.2)') JYEAR0,JMON0,JDATE0
    start_time = (/JYEAR0,JMON0,JDATE0,0,0,0/)
    end_time   = (/JYEAR, JMON, JDATE, 0,0,0/)

    filename = 'emis_accumulated_3d_'//datestamp0//'_'//datestamp//'.nc'
    write(6,'(a)') f90file//':'//subr//': creating file: '//trim(filename)
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')

    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title', &
         'Accumulated emissions in Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    !//---------------------------------------------------------------------

    !// Define spatial dimensions (lat, lon, lev)
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')
    status = nf90_def_dim(ncid,"lev",LPAR,lev_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev dim')
    !// Define spatial dimensions (ilat, ilon, ilev)
    status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon dim')
    status = nf90_def_dim(ncid,"ilev",LPAR+1,ilev_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev dim')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,"date_size",size(start_time),date_size_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')

    !// Define tracer names
    status = nf90_def_dim(ncid,"emcomps",Ncount,emcomp_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define emcomps dim')

    !// Define length of tracer name string
    status = nf90_def_dim(ncid,"tracer_name_len", &
         TNAMELEN, tracer_name_len_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')
    !//---------------------------------------------------------------------

    !// Define the lon/lat/lev
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    status = nf90_def_var(ncid,"lev",nf90_double,lev_dim_id,lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev variable')

    !// Define ilon/ilat/ilev
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    status = nf90_def_var(ncid,"ilev",nf90_double,ilev_dim_id,ilev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev variable')

    !// Putting attributes to lon/lat/lev variables
    status = nf90_put_att(ncid,lon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lev')
    !// Putting attributes to ilon/ilat variables
    status = nf90_put_att(ncid,ilon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilev')

    !// Defining hybrid sigma coordinate A
    status = nf90_def_var(ncid,"ihya",nf90_double,ilev_dim_id,ihya_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihya variable')
    !// Attributes
    status = nf90_put_att(ncid,ihya_id,'units','hPa')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihya')
    status = nf90_put_att(ncid,ihya_id,'description', &
         'Sigma hybrid coordinate A.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihya')
    status = nf90_put_att(ncid,ihya_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihya')

    !// Defining hybrid sigma coordinate B
    status = nf90_def_var(ncid,"ihyb",nf90_double,ilev_dim_id,ihyb_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihyb variable')
    !// Attributes
    status = nf90_put_att(ncid,ihyb_id,'units','1')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'description', &
         'Sigma hybrid coordinate B.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihyb')

    !// Need info on accumulated time period
    status = nf90_def_var(ncid,"delta_time",nf90_double,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define time variable')
    status = nf90_put_att(ncid,time_id,'units', &
            'accumulation period given in seconds')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units time_delta')

    !// Start date for accumulating data - START_TIME
    status = nf90_def_var(ncid,"START_TIME", nf90_int, &
         date_size_dim_id, start_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_TIME variable')
    status = nf90_put_att(ncid,start_time_id,'description', &
         'Start date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_TIME')

    !// End date for accumulating data - END_TIME
    status = nf90_def_var(ncid,"END_TIME", nf90_int, &
         date_size_dim_id, end_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_TIME variable')
    status = nf90_put_att(ncid,end_time_id,'description', &
         'End date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_TIME')

    !// Tracer names of emitted species
    status = nf90_def_var(ncid,"tracer_name",nf90_char, &
         (/tracer_name_len_dim_id,emcomp_dim_id/),tracer_name_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Name of emitted tracer/species.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Grid area (r8), deflate netcdf4
    status = nf90_def_var(ncid,"gridarea", nf90_double, &
         (/lon_dim_id, lat_dim_id/), areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')

    !// Define data component by component
    do N = 1, NPAR
       if (totemis(N) .gt. 0._r8) then
          status = nf90_def_var(ncid,trim(TNAME(N)), nf90_float, &
               (/lon_dim_id, lat_dim_id, lev_dim_id/), comps_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define '//trim(TNAME(N))//' variable')
          status = nf90_def_var_deflate(ncid, comps_id(N), &
               nc4shuffle, 1, nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate '//trim(TNAME(N))//' variable')
          status = nf90_put_att(ncid,comps_id(N),'units','accumulated kg of species')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute units '//trim(TNAME(N)))
       end if !// if (totemis(N) .gt. 0._r8) then
    end do

    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')
    status = nf90_put_var(ncid,lev_id,ZGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lev')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')
    status = nf90_put_var(ncid,ilev_id,ZEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilev')

    status = nf90_put_var(ncid,ihya_id,ETAA)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihya')
    status = nf90_put_var(ncid,ihyb_id,ETAB)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihyb')

    !// delta-time in seconds
    DTDIAG = real(NDAY + 1 - NDAY0, r8) * 84600._r8
    status = nf90_put_var(ncid,time_id,DTDIAG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NRAVG')
    !// Start time
    status = nf90_put_var(ncid,start_time_id,start_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_TIME')
    !// End time
    status = nf90_put_var(ncid,end_time_id,end_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_TIME')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,tracer_name(1:Ncount))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Grid area
    status = nf90_put_var(ncid, areaxy_id, AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting gridarea')

    !// Component by component
    do N = 1, NPAR
       emis3d(:,:,:) = 0._r8
       if (totemis(N) .gt. 0._r8) then
          do MP = 1, MPBLK
            do J = MPBLKJB(MP), MPBLKJE(MP)
              JJ   = J - MPBLKJB(MP) + 1
              do I = MPBLKIB(MP), MPBLKIE(MP)
                II   = I - MPBLKIB(MP) + 1
                do L = 1, LPAR
                  emis3d(I,J,L) = DIAGEMIS_IJ(L,N,II,JJ,MP)
                end do
              end do
            end do
          end do
          status = nf90_put_var(ncid, comps_id(N),emis3d)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting '//trim(TNAME(N)))
       end if !// if (totemis(N) .gt. 0._r8) then
    end do !// do N = 1, NPAR

    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------


    !// Keep track of MEGAN emissions
    call megan_emis2file(NDAY)

    !// --------------------------------------------------------------------
  end subroutine tnd_emis2file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine init_chembud()
    !// --------------------------------------------------------------------
    !// Initialize CHEMLOSS/CHEMPROD/OxCHEMLOS/OxCHEMPROD
    !//
    !// Amund Sovde Haslerud, February 2017
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Re-initialize
    CHEMLOSSMASS(:,:,:,:,:) = 0._r8
    CHEMPRODMASS(:,:,:,:,:) = 0._r8
    OxCHEMLOSSMASS(:,:,:) = 0._r8
    OxCHEMPRODMASS(:,:,:) = 0._r8
    !// --------------------------------------------------------------------
  end subroutine init_chembud
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine save_chemPL(CHEMLOSS,CHEMPROD,DV,PUEDG, I,J,II,JJ,MP,L1,L2,LMT)
    !// --------------------------------------------------------------------
    !// Accumulate chemistry losses, convert to units [kg].
    !//
    !// Ole Amund Sovde, 2014
    !// --------------------------------------------------------------------
    use cmn_size, only: TRACER_ID_MAX
    use cmn_chem, only: TMASS
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J,II,JJ,MP,L1,L2,LMT
    real(r8), intent(in) :: CHEMLOSS(nchemdiag,TRACER_ID_MAX,LPAR)
    real(r8), intent(in) :: CHEMPROD(nchemdiag,TRACER_ID_MAX,LPAR)
    real(r8), intent(in) :: DV(LPAR)    !// Volume [m3]
    real(r8), intent(in) :: PUEDG(LPAR) !// Box top pressure [hPa]
    !// Locals
    real(r8)  :: rfac, lostmass, gainedmass, CLOSS(LPAR)
    integer :: L, M, N, TRNR, TRID
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'save_chemPL'
    !// --------------------------------------------------------------------

    do L = L1, L2
       do N = 1, ncPL
          TRID = compsPL(N)
          TRNR = trsp_idx(TRID)
          !// Skip if not included
          if (TRNR <= 0) cycle

          !// Factor to convert from molec/cm3 to kg/gridbox
          rfac = 1.e6_r8 / AVOGNR * TMASS(TRNR) * 1.e-3_r8 * DV(L)

          do M = 1, nchemdiag
             !// Save all losses, but convert from molec/cm3 to
             !// kg/gridbox. Numbers are accumulated for current time step.
             lostmass = CHEMLOSS(M,TRID,L) * rfac
             CHEMLOSSMASS(M,TRNR,L,I,J) = CHEMLOSSMASS(M,TRNR,L,I,J) &
                  + lostmass
          end do !// do M = 1, nchemdiag
          do M = 1, nchemdiag
             !// Save all prods, but convert from molec/cm3 to
             !// kg/gridbox. Numbers are accumulated for current time step.
             gainedmass = CHEMPROD(M,TRID,L) * rfac
             CHEMPRODMASS(M,TRNR,L,I,J) = CHEMPRODMASS(M,TRNR,L,I,J) &
                  + gainedmass
          end do !// do M = 1, nchemdiag
       end do !// do N = 1, nco
    end do !// do L = L1, L2

    !// N2O loss specifically (for lifetime calculation)
    TRID = 107
    TRNR = trsp_idx(TRID)
    if (TRNR > 0 .and. L2 >= LMT) then
       !// Only do this for stratosphere (L2 >= LMT)
       CLOSS(:) = CHEMLOSS(1,TRID,:)
       call n2o_loss3(TRNR,CLOSS,DV,PUEDG,I,J,II,JJ,MP,L1,L2,LMT)
    end if

    !// CH4 loss specifically (for lifetime calculation)
    TRID = 46
    TRNR = trsp_idx(TRID)
    if (TRNR > 0) then
       CLOSS(:) = CHEMLOSS(1,TRID,:)
       call ch4_loss3(TRNR,CLOSS,DV,PUEDG,I,J,II,JJ,MP,L1,L2,LMT)
    end if

    !// --------------------------------------------------------------------
  end subroutine save_chemPL
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine save_chemOxPL(OxCHEMLOSS,OxCHEMPROD,DV,I,J,L1,L2)
    !// --------------------------------------------------------------------
    !// Accumulate Ox chemistry losses, convert to units [kg].
    !//
    !// Amund Sovde Haslerud, February 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: AVOGNR
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I,J,L1,L2
    real(r8), intent(in) :: OxCHEMLOSS(LPAR)
    real(r8), intent(in) :: OxCHEMPROD(LPAR)
    real(r8), intent(in) :: DV(LPAR)    !// Volume [m3]
    !// Locals
    real(r8)  :: rfac, lostmass, gainedmass
    integer :: L
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'save_chemOxPL'
    !// --------------------------------------------------------------------

    do L = L1, L2

       !// Factor to convert from molec/cm3 to kg/gridbox
       rfac = 1.e6_r8 / AVOGNR * 48._r8 * 1.e-3_r8 * DV(L)

       !// Save loss, but convert from molec/cm3 to
       !// kg/gridbox. Numbers are accumulated for current time step.
       lostmass = OxCHEMLOSS(L) * rfac
       OxCHEMLOSSMASS(L,I,J) = OxCHEMLOSSMASS(L,I,J) &
            + lostmass
       !// Save prod, but convert from molec/cm3 to
       !// kg/gridbox. Numbers are accumulated for current time step.
       gainedmass = OxCHEMPROD(L) * rfac
       OxCHEMPRODMASS(L,I,J) = OxCHEMPRODMASS(L,I,J) &
            + gainedmass

    end do !// do L = L1, L2

    !// --------------------------------------------------------------------
  end subroutine save_chemOxPL
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine chembud_output(YEAR,MONTH,DATE,NDAY)
    !// --------------------------------------------------------------------
    !// Write chemistry budgets to file. Called from pmain.
    !// 1. Burden   (kg/gridbox)
    !// 2. CHEMLOSS (kg/gridbox)
    !// 3. CHEMPROD (kg/gridbox)
    !//
    !// To calculate molec/cm3, these fields must be combined with
    !// volume, so budgets and averages should use same calendar set up.
    !//
    !// Amund Sovde Haslerud, February 2017, February 2016, February 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: STT, AREAXY
    use cmn_diag, only: NDAY0
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR, MONTH, DATE, NDAY
    !// Locals
    integer :: I,J,L, M, N, TRNR, TRID, NCOUNT, NINCID(nchemdiag)
    integer :: ifnr
    character(len=80) :: filename
    character(len=8) :: datestamp
    real(r8) :: RTMP(IPAR,JPAR,LPAR,nchemdiag)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'chembud_output'
    !// --------------------------------------------------------------------

    !// Find file number
    ifnr = get_free_fileid()

    !// YYYYMMDD for filenames, taken from JYEAR,JMON,JDATE, are
    !// the arguments of this subroutine.
    !// They are already updated by calendar.
    write(datestamp(1:8),'(i4.4,2i2.2)') YEAR,MONTH,DATE
    filename = 'chemistryBPL_'//datestamp//'.dta'

    !// Total number of tracers diagnosed and their IDs
    NCOUNT = 0
    NINCID(:) = 0
    do N = 1, ncPL
       TRID = compsPL(N)
       if (trsp_idx(TRID) .gt. 0) then
          NINCID(N) = TRID
          NCOUNT = NCOUNT + 1
       end if
    end do

    !// If no tracers were diagnosed, skip making the file
    if (NCOUNT .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': No chemical prod/loss diagnosed, skipping '//trim(filename)
    else
       !// Write file
       open(ifnr,file=filename,form='unformatted')
       write(ifnr) IPAR,JPAR,LPAR,NDAY0,NDAY
       write(ifnr) AREAXY

       write(ifnr) NCOUNT
       write(ifnr) NINCID(1:NCOUNT)

       do N = 1, ncPL
          TRID = compsPL(N)
          TRNR = trsp_idx(TRID)
          !// Skip if tracer not included
          if (TRNR .le. 0) cycle

          !// Write component number
          write(ifnr) TRID
          !// Burden
          write(ifnr) real(STT(:,:,:,TRNR), r4) !// Write as float
          !// Lost mass in chemistry
          do L = 1, LPAR
             do J = 1, JPAR
                do I = 1, IPAR
                   do M = 1, nchemdiag
                      RTMP(I,J,L,M) = CHEMLOSSMASS(M,TRNR,L,I,J)
                   end do
                end do
             end do
          end do
          do M = 1, nchemdiag
             if (sum(RTMP(:,:,:,M)) .eq. 0._r8) then
                write(ifnr) M, 0 !// Flag 0 denotes no data
                write(ifnr) -1   !// To save space, write -1
             else
                write(ifnr) M, 1 !// Flag 1 denotes data
                write(ifnr) real(RTMP(:,:,:,M), r4) !// Write as float
             end if
          end do
          !// Produced mass in chemistry
          do L = 1, LPAR
             do J = 1, JPAR
                do I = 1, IPAR
                   do M = 1, nchemdiag
                      RTMP(I,J,L,M) = CHEMPRODMASS(M,TRNR,L,I,J)
                   end do
                end do
             end do
          end do
          do M = 1, nchemdiag
             if (sum(RTMP(:,:,:,M)) .eq. 0._r8) then
                write(ifnr) M, 0 !// Flag 0 denotes no data
                write(ifnr) -1   !// To save space, write -1
             else
                write(ifnr) M, 1 !// Flag 1 denotes data
                write(ifnr) real(RTMP(:,:,:,M), r4) !// Write as float
             end if
          end do

       end do
       close(ifnr)
       write(6,'(a)') f90file//':'//subr// &
            ': Wrote chemical prod/loss/burden: '//trim(filename)
    end if

    !// Ox Prod and Loss
    filename = 'chemistryOxPL_'//datestamp//'.dta'

    !// Write file
    open(ifnr,file=filename,form='unformatted')
    write(ifnr) IPAR,JPAR,LPAR,NDAY0,NDAY
    write(ifnr) AREAXY

    !// Produced mass in chemistry
    write(ifnr) real(OxCHEMPRODMASS, r4) !// Write as float
    write(ifnr) real(OxCHEMLOSSMASS, r4) !// Write as float
    close(ifnr)
    write(6,'(a)') f90file//':'//subr// &
         ': Wrote chemical prod/loss: '//trim(filename)


    !// Re-initialize
    call init_chembud()

    !// --------------------------------------------------------------------
  end subroutine chembud_output
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine chembud_output_nc(NDAY)
    !// --------------------------------------------------------------------
    !// Write chemistry production and loss to file. Called from pmain.
    !// Should also put out volume and molecular weights.
    !//
    !// Units are kg/gridbox
    !//
    !// To calculate molec/cm3, these fields must be combined with
    !// volume, so budgets and averages should use same calendar set up.
    !//
    !// Amund Sovde Haslerud, October 2017, February 2016, February 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_ctm, only: AREAXY, XDGRD, YDGRD, ZGRD, XDEDG, YDEDG, ZEDG, &
         ETAA, ETAB, JYEAR, JMON, JDATE
    use cmn_diag, only: NDAY0, JYEAR0, JMON0, JDATE0, &
         nc4deflate_global, nc4shuffle_global
    use cmn_chem, only: TNAME
    use cmn_oslo, only: trsp_idx
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    !// Locals
    integer :: I,J,L, M, N, TRNR, TRID, NCOUNT
    integer, dimension(6) :: start_time, end_time
    character(len=80) :: filename, varname
    character(len=8) :: datestamp, datestamp0
    real(r8) :: RTMP(IPAR,JPAR,LPAR,nchemdiag)
    real(r8) :: DTDIAG
    logical :: Lox

    integer :: &
         lat_dim_id, lon_dim_id, lev_dim_id, nchemdiag_dim_id, time_dim_id, &
         lat_id, lon_id, lev_id, nchemdiag_id, time_id, &
         ilat_dim_id, ilon_dim_id, ilev_dim_id, &
         ilat_id, ilon_id, ilev_id, &
         ihya_dim_id, ihyb_dim_id, &
         ihya_id, ihyb_id, &
         areaxy_id, &
         date_size_dim_id, &
         start_time_id, end_time_id, & !IDs for start/end dates for average
         status, ncid
    integer, dimension(ncPL) :: comploss_id, compprod_id
    integer :: oxloss_id, oxprod_id

    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'chembud_output_nc'
    !// --------------------------------------------------------------------


    !// filename: start-date to end date (both 00UTC)
    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
    write(datestamp0(1:8),'(i4.4,2i2.2)') JYEAR0,JMON0,JDATE0

    !// Start/end of period at 00:00UTC
    start_time = (/JYEAR0,JMON0,JDATE0,0,0,0/)
    end_time = (/JYEAR,JMON,JDATE,0,0,0/)

    filename = 'chemistryPL_'//datestamp0//'_'//datestamp//'.nc'

    !// Total number of tracers diagnosed and their IDs
    NCOUNT = 0
    do N = 1, ncPL
       TRID = compsPL(N)
       if (trsp_idx(TRID) .gt. 0) then
          NCOUNT = NCOUNT + 1
       end if
    end do

    if (sum(OxCHEMPRODMASS) .eq. 0._r8 .and. sum(OxCHEMLOSSMASS) .eq. 0._r8) &
         Lox = .false.

    !// If no tracers were diagnosed, skip making the file
    if (NCOUNT .eq. 0 .and. .not. Lox) then
       write(6,'(a)') f90file//':'//subr// &
            ': No chemical prod/loss diagnosed, skipping '//trim(filename)
       return
    end if

    !// Write file

    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')

    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title', &
         'Chemistry budget terms from Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': title')
    !//---------------------------------------------------------------------

    !// Define spatial dimensions (lat, lon, lev)
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')
    status = nf90_def_dim(ncid,"lev",LPAR,lev_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev dim')
    status = nf90_def_dim(ncid,"nchemdiag",nchemdiag_id,nchemdiag_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define nchemdiag dim')


    !// Define ilon/ilat/ilev
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    status = nf90_def_var(ncid,"ilev",nf90_double,ilev_dim_id,ilev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev variable')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,"date_size",size(start_time),date_size_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')

    !//---------------------------------------------------------------------

    !// Define the lon/lat/lev
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    status = nf90_def_var(ncid,"lev",nf90_double,lev_dim_id,lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev variable')

    !// Define ilon/ilat
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    !// Define time
    status = nf90_def_var(ncid,"time",nf90_double,time_dim_id,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define time variable')

    !// Putting attributes to lon/lat/lev variables
    status = nf90_put_att(ncid,lon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lev')

    !// Putting attributes to ilon/ilat variables
    status = nf90_put_att(ncid,ilon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilev')

    !// Defining hybrid sigma coordinate A
    status = nf90_def_var(ncid,"ihya",nf90_double,ilev_dim_id,ihya_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihya variable')
    !// Attributes
    status = nf90_put_att(ncid,ihya_id,'units','hPa')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihya')
    status = nf90_put_att(ncid,ihya_id,'description', &
         'Sigma hybrid coordinate A.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihya')
    status = nf90_put_att(ncid,ihya_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihya')

    !// Defining hybrid sigma coordinate B
    status = nf90_def_var(ncid,"ihyb",nf90_double,ilev_dim_id,ihyb_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihyb variable')
    !// Attributes
    status = nf90_put_att(ncid,ihyb_id,'units','1')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'description', &
         'Sigma hybrid coordinate B.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihyb')

    !// Info on accumulated time period
    status = nf90_def_var(ncid,"delta_time",nf90_double,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define time variable')
    status = nf90_put_att(ncid,time_id,'units', &
         'accumulation period given in seconds')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units time_delta')

    !// Start date for accumulating data - START_TIME
    status = nf90_def_var(ncid,"START_TIME", nf90_int, &
         date_size_dim_id, start_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_TIME variable')
    status = nf90_put_att(ncid,start_time_id,'description', &
         'Start date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_TIME')

    !// End date for accumulating data - END_TIME
    status = nf90_def_var(ncid,"END_TIME", nf90_int, &
         date_size_dim_id, end_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_TIME variable')
    status = nf90_put_att(ncid,end_time_id,'description', &
         'End date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_TIME')

    !// Grid area (r8), deflate netcdf4
    status = nf90_def_var(ncid,"gridarea", nf90_double, &
         (/lon_dim_id, lat_dim_id/), areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')

    !// For each tracer put out loss and prod
    do N = 1, ncPL
       !// Only put out if included
       if (trsp_idx(compsPL(N)) .le. 0) cycle

       varname = trim(TNAME(trsp_idx(compsPL(N))))//'_LOSS'
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/lon_dim_id, lat_dim_id, lev_dim_id, nchemdiag_dim_id/), &
            comploss_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define '//trim(varname)//' variable')
       status = nf90_def_var_deflate(ncid, comploss_id(N), &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate '//trim(varname)//' variable')
       status = nf90_put_att(ncid,comploss_id(N),'units','accumulated kg of species')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))

       varname = trim(TNAME(trsp_idx(compsPL(N))))//'_PROD'
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/lon_dim_id, lat_dim_id, lev_dim_id, nchemdiag_dim_id/), &
            compprod_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define '//trim(varname)//' variable')
       status = nf90_def_var_deflate(ncid, compprod_id(N), &
            nc4shuffle, 1, nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate '//trim(varname)//' variable')
       status = nf90_put_att(ncid,compprod_id(N),'units','accumulated kg of species')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))

    end do !// do N = 1, ncPL

    !// Ox
    status = nf90_def_var(ncid, 'Ox_LOSS', nf90_float, &
         (/lon_dim_id, lat_dim_id, lev_dim_id/), oxloss_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define Ox_LOSS variable')
    status = nf90_def_var_deflate(ncid, oxloss_id, &
         nc4shuffle, 1, nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate Ox_LOSS variable')
    status = nf90_put_att(ncid,oxloss_id,'units','accumulated kg of species')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units Ox_LOSS')

    status = nf90_def_var(ncid, 'Ox_PROD', nf90_float, &
         (/lon_dim_id, lat_dim_id, lev_dim_id/), oxprod_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define Ox_PROD variable')
    status = nf90_def_var_deflate(ncid, oxprod_id, &
         nc4shuffle, 1, nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate Ox_PROD variable')
    status = nf90_put_att(ncid,oxprod_id,'units','accumulated kg of species')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units Ox_PROD')


    !// --------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !// --------------------------------------------------------------------

    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')
    status = nf90_put_var(ncid,lev_id,ZGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lev')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')
    status = nf90_put_var(ncid,ilev_id,ZEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilev')

    status = nf90_put_var(ncid,ihya_id,ETAA)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihya')
    status = nf90_put_var(ncid,ihyb_id,ETAB)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihyb')

    !// delta-time in seconds
    DTDIAG = real(NDAY + 1 - NDAY0, r8) * 84600._r8
    status = nf90_put_var(ncid,time_id,DTDIAG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NRAVG')
    !// Start time
    status = nf90_put_var(ncid,start_time_id,start_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_TIME')
    !// End time
    status = nf90_put_var(ncid,end_time_id,end_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_TIME')

    !// Grid area
    status = nf90_put_var(ncid, areaxy_id, AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting gridarea')

    !// For each tracer put out loss and prod
    do N = 1, ncPL

       !// Only put out if included
       if (trsp_idx(compsPL(N)) .le. 0) cycle

       varname = trim(TNAME(N))//'_LOSS'

       !// Lost mass in chemistry
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                do M = 1, nchemdiag
                   RTMP(I,J,L,M) = CHEMLOSSMASS(M,TRNR,L,I,J)
                end do
             end do
          end do
       end do

       status = nf90_put_var(ncid, comploss_id(N), RTMP)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting '//trim(varname))

       !// Produced mass in chemistry
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                do M = 1, nchemdiag
                   RTMP(I,J,L,M) = CHEMPRODMASS(M,TRNR,L,I,J)
                end do
             end do
          end do
       end do

       status = nf90_put_var(ncid, compprod_id(N), RTMP)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting '//trim(varname))

    end do


    !// Ox_LOSS
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             RTMP(I,J,L,1) = OxCHEMLOSSMASS(L,I,J)
          end do
       end do
    end do
    status = nf90_put_var(ncid, oxloss_id, RTMP(:,:,:,1))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting Ox_LOSS')

    !// Ox_PROD
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             RTMP(I,J,L,1) = OxCHEMPRODMASS(L,I,J)
          end do
       end do
    end do
    status = nf90_put_var(ncid, oxprod_id, RTMP(:,:,:,1))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting Ox_PROD')


    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr// &
         ': Wrote chemical prod/loss: '//trim(filename)



    !// Re-initialize
    call init_chembud()

    !// --------------------------------------------------------------------
  end subroutine chembud_output_nc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module diagnostics_general
!//=========================================================================

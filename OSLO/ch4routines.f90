!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, February 2017
!//=========================================================================
!// CH4 routines.
!//=========================================================================
module ch4routines
  !// ----------------------------------------------------------------------
  !// MODULE: ch4routines
  !// DESCRIPTION: Routines for setting CH4 at surface and for soil uptake
  !//              when using emissions.
  !//
  !// Due to its lifetime, CH4 is usually not emitted, but set at the
  !// surface. This routine reads CH4 surface fields for this, but also
  !// for setting CH4 soil uptake in case CH4 emissions are used.
  !//
  !// Fixed surface CH4
  !// Default is CH4 from from the HYMN project (currently fixed to 2003),
  !// where CH4 was emitted at the surface and as good as in steady state.
  !// You can also use the old versions for POET and RETRO emissions.
  !//
  !// Main routine is oc_update_ch4surface.
  !//
  !// Also a routine for setting 3D field, to "boost" up a
  !// previously slightly low CH4.
  !//
  !// If pre-industrial simulations are to be carried out, CH4 must
  !// be set in other ways.
  !//
  !//
  !// Emitted CH4
  !// Emissions are treated as other emissions, but you need to set
  !// parameter METHANEMIS=.true. in oc_globalvariables.f.
  !// Here you find routines for setting the soil uptake (dry deposition),
  !// which accompagnies emissions.
  !//
  !// Contains:
  !//   subroutine update_ch4surface
  !//   subroutine ch4surface_hymn
  !//   subroutine ch4surface_scale_hymn
  !//   subroutine ch4surface_poet
  !//   subroutine ch4surface_retro
  !//   subroutine READCH4
  !//   subroutine setch4sfc
  !//   subroutine set_ch4_stt
  !//   subroutine updateSOILUPTAKEbousquet
  !//   subroutine read_ch4sfc4soiluptake
  !//   subroutine read_ch4bousquet
  !//   subroutine ch4drydep_bousquet
  !//   subroutine reportsfcch4
  !//
  !// Amund Sovde, October 2013
  !//   Included routines for CH4 surface uptake.
  !// Amund Sovde, March 2011
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IDBLK, JDBLK, MPBLK
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Manually set the type of CH4 field to use:
  !// 1: HYMN (T42, will be interpolated to model resolution)
  !// 2: POET
  !// 3: RETRO
  integer, parameter :: CH4TYPE = 1

  !// Fields for CH4 soil uptake (i.e. dry deposition process in the model)
  real(r8) :: CH4SOILUPTAKE(IDBLK,JDBLK,MPBLK)
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'ch4routines.f90'
  !// ----------------------------------------------------------------------
  !// All variables are to be saved. This should be default for global
  !// variables, although Fortran claims it may not be.
  save
  private
  public update_ch4surface, setch4sfc, set_ch4_stt, scale_ch4_stt, &
       updateSOILUPTAKEbousquet, ch4drydep_bousquet, reportsfcch4
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine update_ch4surface()
    !// --------------------------------------------------------------------
    !// Controls which surface CH4 field to use.
    !//
    !// Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_oslo, only: trsp_idx, Xtrsp_idx, CH4FIELD, METHANEMIS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='update_ch4surface'
    !// --------------------------------------------------------------------

    !// Skip if METHANEMIS == .true.
    if (METHANEMIS) then
       CH4FIELD(:,:) = -1._r8
       return
    end if

    !// Check if CH4 is included
    if (.not. (trsp_idx(46) .gt. 0 .or. Xtrsp_idx(46) .gt. 0)) then
       CH4FIELD(:,:) = 0._r8
       write(6,'(a)') f90file//':'//subr// &
            ': CH4 not included, no need to set surface CH4'
       return
    end if

    !// Update CH4FIELD
    if (CH4TYPE .eq. 1) then
       !// HYMN (default)
       call ch4surface_hymn()
       !// Scale hymn data - uses MYEAR for scaling.
       !call ch4surface_scale_hymn()
    else if (CH4TYPE .eq. 2) then
       !// POET
       call ch4surface_poet(2000)
    else if (CH4TYPE .eq. 3) then
       !// RETRO
       call ch4surface_retro(2000)
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': No such surface field type for CH4: ',CH4TYPE
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine update_ch4surface
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4surface_hymn()
    !// --------------------------------------------------------------------
    !// Read surface CH4 from HYMN files and put field into CH4FIELD.
    !// Input is T42.
    !//
    !// Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR
    use cmn_ctm, only:  areaxy, xdedg, ydedg, JMON, TMON
    use cmn_met, only: MYEAR
    use cmn_parameters, only: a0, CPI180, ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use cmn_oslo, only: CH4FIELD
    use ncutils, only: readnc_3d_from4d
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// Local parameters
    integer,parameter :: IT42=128, JT42=64
    integer :: CH4YEAR, J, I
    real(r8) :: INCH8(IT42,JT42),T42R8XY(IT42,JT42), dx
    real(r8) :: T42AREA, T42XDEDG(IT42+1), T42YDEDG(JT42+1)
    real(r8) :: WGAULAT(JT42),WGAUWT(JT42), WYEDG1P1(JT42+1)
    real(r8) :: YEDG1P1(JT42+1)

    character(len=100) :: FILENAME

    character(len=4) :: CYEAR
    character(len=3),parameter :: CH4LAB='CH4'
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='ch4surface_hymn'
    !// --------------------------------------------------------------------

    !// Fields only available for 2003-2005
    if (MYEAR .lt. 2003) then
       CH4YEAR = 2003
    else if (MYEAR .gt. 2005) then
       CH4YEAR = 2005
    else
       CH4YEAR = MYEAR
    end if

    !// Keep fixed at 2003
    CH4YEAR = 2003

    Write(CYEAR,'(I4)') CH4YEAR

    write(6,'(a)') f90file//':'//subr//': Reading CH4 surface field'

    !// Data are on mixing ratio, T42 resolution
    INCH8(:,:) = 0._r8
    FILENAME = 'Indata_CTM3/ch4_fields/CH4_2D_'//CYEAR//'_'//TMON//'.nc'

    !// Read the file (INCH8 must be R8)
    !// Month to get is '1', since there is only one month on each file.
    call readnc_3d_from4d(filename,'lon',IT42,'lat',JT42,'lev',1, &
         'time', 1, CH4LAB, INCH8)

    if (IPAR .ne. IT42) then
       !// Interpolate - must multiply by area before interpolating
       !// Find input edges
       !// Zonal
       dx = 360._r8/real(IT42,r8)
       T42XDEDG(1) = -0.5*dx
       do J=2,IT42+1
          T42XDEDG(J) = T42XDEDG(J-1) + dx
       end do


       !// Meridional (gaussian, as in p-grid.f)
       call GAUSST2(JT42,WGAULAT,WGAUWT)
       WYEDG1P1(1) = -1._r8
       do J = 2,JT42/2
          WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
       end do
       WYEDG1P1(JT42/2+1) = 0._r8
       do J = JT42/2+2,JT42+1
          WYEDG1P1(J) = -WYEDG1P1(JT42+2-J)
       end do
       !// T42 Y-edged (degrees)
       do J = 1,JT42+1
          T42YDEDG(J) = ZPI180 * asin(WYEDG1P1(J))
       end do

       !// Grid box areas, multiply mmr-field with area
       do J=1,JT42
          T42AREA =  A0*A0 * CPI180*(T42XDEDG(2)-T42XDEDG(1)) &
               * (sin(CPI180*T42YDEDG(J+1)) - sin(CPI180*T42YDEDG(J)))

          T42R8XY(:,J) = INCH8(:,J) * T42AREA
       end do


       Call E_GRID(T42R8XY,T42XDEDG,T42YDEDG,IT42,JT42, &
            CH4FIELD,XDEDG,YDEDG,IPAR,JPAR,1)

       !// Finally divide by current area
       CH4FIELD(:,:) = CH4FIELD(:,:)/AREAXY(:,:)

    else

       !// Keep T42 field
       do J = 1, JPAR
          do I = 1, IPAR
             CH4FIELD(I,J) = INCH8(I,J)
          end do
       end do

    end if

    write(6,'(a)') '* Surface CH4 updated to HYMN fields.'

    !// --------------------------------------------------------------------
  end subroutine ch4surface_hymn
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4surface_scale_hymn()
    !// --------------------------------------------------------------------
    !// Scale HYMN surface CH4 according to year, based on ESRL Global
    !// Monitoring Division marine global annual CH4.
    !// http://www.esrl.noaa.gov/gmd/ccgg/trends_ch4/
    !//
    !// Routine must be called right after CH4FIELD has been set.
    !//
    !// Amund Sovde Haslerud, January 2017
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    use cmn_oslo, only: CH4FIELD
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input

    !// Local parameters
    integer, parameter :: nObs = 33
    real(r8), dimension(nObs), parameter :: ANNUAL_CH4 = &
      (/ 1644.56_r8, 1657.39_r8, 1669.79_r8, 1682.08_r8, 1693.18_r8, &
         1704.10_r8, 1714.00_r8, 1724.64_r8, 1735.35_r8, 1736.38_r8, &
         1741.84_r8, 1748.73_r8, 1751.02_r8, 1754.36_r8, 1765.34_r8, &
         1772.00_r8, 1773.00_r8, 1771.04_r8, 1772.59_r8, 1776.92_r8, &
         1776.94_r8, 1773.96_r8, 1774.69_r8, 1781.14_r8, 1786.77_r8, &
         1793.21_r8, 1798.65_r8, 1803.01_r8, 1808.23_r8, 1813.32_r8, &
         1822.49_r8, 1833.99_r8, 1842.99_r8 /)
    integer, dimension(nObs), parameter :: ANNUAL_YEAR = &
       (/ 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, &
          1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, &
          2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, &
          2014, 2015, 2016 /)

    integer :: YSC, YCH4
    real(r8) :: SCALE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='ch4surface_scale_hymn'
    !// --------------------------------------------------------------------

    !// Scale to observed year, but not outside the range
    if (MYEAR .lt. ANNUAL_YEAR(1)) then
       YSC = 1 !// Use the first year (1984)
    else if (MYEAR .gt. maxval(ANNUAL_YEAR)) then
       YSC = nObs !// Use the last year (2015)
    else
       !// Scale to year
       do YSC = 1, nObs
          if (MYEAR .eq. ANNUAL_YEAR(YSC)) then
             exit
          end if
       end do
    end if

    !// Scale from year (2003)
    do YCH4 = 1, nObs
       if (2003 .eq. ANNUAL_YEAR(YCH4)) then
          exit
       end if
    end do

    SCALE = ANNUAL_CH4(YSC)/ANNUAL_CH4(YCH4)
    CH4FIELD(:,:) = CH4FIELD(:,:) * SCALE

    write(6,'(a24,i4,a4,i4,a3,f8.4)') 'CH4@surface scaled from ',&
         ANNUAL_YEAR(YCH4),' to ',ANNUAL_YEAR(YSC),' : ',scale

    !// --------------------------------------------------------------------
  end subroutine ch4surface_scale_hymn
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine ch4surface_poet(YEAR)
    !// --------------------------------------------------------------------
    !// Read CH4 as in POET emissions.
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR
    use cmn_ctm, only:  areaxy, xdedg, ydedg, JMON
    use cmn_parameters, only: a0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: CH4FIELD
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR

    !// Local parameters
    integer, parameter :: IEMIS=360, JEMIS=180
    integer :: NRCOMP, I,J,L,M,N,IWIDTH
    real(r8) :: EMISDATA(IEMIS,JEMIS), &
         CTMDATA(IPAR,JPAR)

    real(r8) :: XBEDGE(IEMIS+1), YBEDGE(JEMIS+1), XYBOX(JEMIS), &
         DELXEMIS,DELYEMIS

    character(len=80) :: FILNAVN, PATH
    character(len=4) :: CYEAR 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='ch4surface_poet'
    !// --------------------------------------------------------------------

    write(CYEAR,'(I4)') YEAR
  
    write(6,'(a)') f90file//':'//subr//': Reading CH4 surface field (POET)'

    EMISDATA(:,:) = 0._R8
    PATH='Indata_CTM3/'
    IWIDTH = Len_Trim(PATH)

      
    FILNAVN='ch4_'//CYEAR


    Call READCH4(EMISDATA,FILNAVN,PATH,JMON,IEMIS,JEMIS,IWIDTH)


    !// UCI style interpolation

    !// Edges given in degrees
    DELXEMIS = 360._r8/real(IEMIS,r8)
    Do I=1,IEMIS+1
       XBEDGE(I) = DELXEMIS*real(I-1,r8)
    End Do

    DELYEMIS = 180._r8 / real(JEMIS,r8)
    YBEDGE(1) = -90._r8
    Do J=1,JEMIS
       YBEDGE(J+1) = YBEDGE(J) + DELYEMIS
    End Do

    !// dataset can not be per area - need to convert from /cm2

    !// multiply with area [m2]
    do J = 1, JEMIS
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    do J = 1, JEMIS
       EMISDATA(:,J) = EMISDATA(:,J)*XYBOX(J)
    end do

    call E_GRID(EMISDATA,XBEDGE,YBEDGE,IEMIS,JEMIS, &
         CTMDATA,XDEDG,YDEDG,IPAR,JPAR,1)


    !// Convert back to per area
    CTMDATA(:,:) = CTMDATA(:,:)/AREAXY(:,:)


    CH4FIELD(:,:) = CTMDATA(:,:) * 1.e-9_r8


    !// --------------------------------------------------------------------
  end subroutine ch4surface_poet
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4surface_retro(YEAR)
    !// --------------------------------------------------------------------
    !// Read CH4 as in RETRO emissions.
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR
    use cmn_ctm, only:  areaxy, xdedg, ydedg, JMON
    use cmn_parameters, only: CPI180
    use cmn_oslo, only: CH4FIELD
    use regridding, only: E_GRID_Y
    use ncutils, only: readnc_1d, readnc_1d_from2d
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: YEAR

    integer, parameter :: LATDIM=41
    integer :: J, PICKMONTH
    real(r8) :: RCH4_IN(LATDIM), RLAT_IN(LATDIM),RCH4_OUT(JPAR)
    real(r8) :: RDEDG(LATDIM+1), delta

    character(len=80) :: FILENAME
    character(len=40) :: FIELDNAME, dimname_lat, dimname_time
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='ch4surface_retro'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Reading CH4 surface field (RETRO)'

    RCH4_IN(:) = 0._r8


    FILENAME = 'Indata_CTM3/ch4_fields/ch4_surface_retro_1800-2100.nc'
    fieldname='ch4'
    dimname_lat ='lat'
    dimname_time='time'
   
    !// Pick month (conatins monthly info for 1800-2100)
    PICKMONTH = 12 * (YEAR - 1800) + JMON

    !// Read latitude dimension
    call readnc_1d(RLAT_IN, dimname_lat, dimname_lat, LATDIM,FILENAME)

    !// Extract 1D-field from netcdf 2D-field
    call readnc_1d_from2d(RCH4_IN, FIELDNAME, dimname_lat, LATDIM, &
         dimname_time, PICKMONTH, FILENAME)


    !// Get RDEDG (edge of grids)
    RDEDG(1) = RLAT_IN(1)
    RDEDG(2) = RLAT_IN(1) + 0.5_r8*(RLAT_IN(2)-RLAT_IN(1))
    do J = 3, LATDIM
       RDEDG(J) = RDEDG(J-1) + (RLAT_IN(J)-RLAT_IN(J-1))
    end do
    RDEDG(LATDIM+1) = RLAT_IN(LATDIM)

    !// Interpolate 1-D
    !// Must multiply with area (original field is vmr)
    !// Since the field is a zonal mean, the horizontal extent is 2PI,
    !// and the area of each latitudinal slice is
    !//       A0*A0*CPI180*2*PI *
    !//           (sin(CPI180*YBDEDG(J+1)) - sin(CPI180*YBDEDG(J)))
    !// A0*A0*CPI180*2*PI is skipped since we divide by it after
    !// interpolation.
    do J = 1, LATDIM
       delta = (sin(CPI180*RDEDG(J+1)) - sin(CPI180*RDEDG(J)))
       RCH4_IN(J) = RCH4_IN(J)*delta
    end do

    call E_GRID_Y(RCH4_IN,RDEDG,LATDIM, RCH4_OUT,YDEDG,JPAR)

    !// Divide by new area
    do J = 1, JPAR
       delta =  (sin(CPI180*YDEDG(J+1)) - sin(CPI180*YDEDG(J)))
       RCH4_OUT(J) = RCH4_OUT(J)/delta
    end do

    !// Put field into CH4FIELD
    do J = 1, JPAR
       CH4FIELD(:,J) = RCH4_OUT(J) * 1.e-9_r8
    end do


    !// --------------------------------------------------------------------
  end subroutine ch4surface_retro
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  Subroutine READCH4 (EMISDATA, FILNAVN, PATH, NUMTIMES,IEMIS,JEMIS,IWIDTH)
    !// --------------------------------------------------------------------
    !// Old routine to read surface CH4 for POET.
    !// --------------------------------------------------------------------
    Implicit None
    !// --------------------------------------------------------------------
    !// Input parameters
    Integer, intent(in) :: NUMTIMES,IEMIS,JEMIS,IWIDTH
    Character(len=*),intent(in)  :: FILNAVN
    Character(len=*), intent(in) :: PATH
    !// Output parameters
    Real(r8), intent(out) :: EMISDATA(IEMIS,JEMIS)

    !// Values declarations
    real(r8) :: out(JEMIS,NUMTIMES), sgrid,tgrid,sgrid2,tgrid2

    integer ::     numLevels,ilevel,igrid,jgrid,itime

    !// Div.
    Integer :: I,J,K,M,YR,MNTH
    character(len=120) :: datafilename

    logical :: fex
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='READCH4'
    !// --------------------------------------------------------------------

    !// Initialize
    numLevels=1 ! Only one level, eventual veritcal allocation done later 


    datafilename = trim(PATH)//trim(FILNAVN)
    inquire(file=datafilename,exist=fex)
    if (.not. fex) then
       write(6,'(a)') f90file//':'//subr//': Cannot find file '//trim(datafilename)
       stop
    end if
    write(6,'(a)') f90file//':'//subr//': Reading '//trim(datafilename)
    open (150,file=datafilename,status='old',form='formatted')


    do itime=1,NUMTIMES    !  (360x180=64800 gridboxes)
       !// ... read in all levels and times for one grid
       !// read (150, 103)
       !// 103 format (/)         
       read (150,*) yr,mnth, (out(jgrid,itime),jgrid=1,JEMIS)
    end do

    Do I=1,IEMIS
       EMISDATA(I,:)=out(:,NUMTIMES)
    end do

    close(150)

    !// --------------------------------------------------------------------
  end Subroutine READCH4
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine setch4sfc(BTT,BTTBCK,AIRB,MP)
    !// --------------------------------------------------------------------
    !// Set CH4 at surface to keep it constant.
    !// CH4FIELD is read in in routine ...
    !//
    !// Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TMOLMIX2MASSMIX
    use cmn_oslo, only: trsp_idx, Xtrsp_idx, XTMOLMIX2MASSMIX, &
         CH4FIELD, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT, BTTBCK
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: AIRB
    integer, intent(in) :: MP
    !// Locals
    integer :: I,II,J,JJ, TNR, XTNR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='setch4sfc'
    !// --------------------------------------------------------------------

    !// Return if CH4 is not included
    if (trsp_idx(46).lt.0 .and. Xtrsp_idx(46).lt.0) return

    TNR = trsp_idx(46)
    XTNR = Xtrsp_idx(46)


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Set ch4 (mass) from vmr
          if (TNR .gt. 0) then
             BTT(1,TNR,II,JJ) = CH4FIELD(I,J)*AIRB(1,II,JJ) &
                                * TMOLMIX2MASSMIX(TNR)
             !// Resetting the diagnostics also
             BTTBCK(1,TNR,II,JJ) =  BTT(1,TNR,II,JJ)
          else if (XTNR .gt. 0) then
             XSTT(1,XTNR,I,J) = CH4FIELD(I,J)*AIRB(1,II,JJ) &
                                * XTMOLMIX2MASSMIX(XTNR)
          end if

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine setch4sfc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_ch4_stt()
    !// --------------------------------------------------------------------
    !// If CH4 in general is low, use this routine to boost up the CH4
    !// to monthly average from HYMN results.
    !// Horizontal resolution on input is T42, will interpolate to current
    !// resolution. So far, there is no vertical interpolation, so this
    !// works only for L60.
    !//
    !// Amund Sovde Haslerud, February 2017
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: STT, AIR, areaxy, xdedg, ydedg, JMON
    use cmn_chem, only: TMOLMIX2MASSMIX
    use cmn_met, only: MYEAR
    use cmn_parameters, only: a0, CPI180, ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use cmn_oslo, only: CH4FIELD, trsp_idx, METHANEMIS
    use ncutils, only: GET_NETCDF_VAR_1D, GET_NETCDF_VAR_3D, GET_NETCDF_VAR_4D
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Local parameters
    integer,parameter :: IT42=128, JT42=64, LT42=60
    integer :: CH4YEAR, I, J, L
    real(r8) :: dx, CTMR8(IPAR,JPAR)
    real(r8) :: TOTAREA, hymnmean

    character(len=100) :: FILENAME

    character(len=4) :: CYEAR

    integer :: nLev, nLon, nLat, nTime

    real(r8), parameter :: meansfc = 0._r8 ! Skip scaling
    !real(r8), parameter :: meansfc = 1682.08e-9_r8 ! Observed 1990
    real(r8), allocatable, dimension(:) :: inLon, inLat, inLev, inTime
    real(r8), allocatable, dimension(:) :: xbedge, ybedge, xybox
    real(r8), allocatable, dimension(:) :: WGAULAT, WGAUWT, WYEDG1P1, YEDG1P1
    real(r8), allocatable, dimension(:,:) :: INR8XY
    real(r8), allocatable, dimension(:,:,:) :: INR8XYZ
    real(r8), allocatable, dimension(:,:,:,:) :: INR8XYZT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='set_ch4_stt'
    !// --------------------------------------------------------------------

    !// Skip if METHANEMIS == .true.
    if (METHANEMIS) return

    !// Check if CH4 is included
    if (trsp_idx(46) .le. 0) return

    !// Check vertical resolution; only works for L60
    if (LPAR .ne. LT42) then
       write(6,'(a)') f90file//':'//subr//': not set up for LPAR=/=LT42'
       stop
    end if

    !// Fields only available for 2003-2005
    !if (MYEAR .lt. 2003) then
    !   CH4YEAR = 2003
    !else if (MYEAR .gt. 2005) then
    !   CH4YEAR = 2005
    !else
    !   CH4YEAR = MYEAR
    !end if

    !// Keep fixed at 2003
    CH4YEAR = 2003

    write(CYEAR,'(I4)') CH4YEAR

    write(6,'(a)') f90file//':'//subr//': January CH4 3D field'


    !// 3D file for January
    FILENAME = 'Indata_CTM3/ch4_fields/CH4_3D_'//CYEAR//'_JAN.nc'

    !// Get data from dataset
    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( filename, 'lon', inLon  )
    call get_netcdf_var_1d( filename, 'lat', inLat  )
    call get_netcdf_var_1d( filename, 'lev', inLev )
    call get_netcdf_var_1d( filename, 'time', inTime )

    nLon = SIZE( inLon )
    nLat = SIZE( inLat )
    nLev = SIZE( inLev )
    nTime = SIZE( inTime )

    !// Deallocate all local variables
    IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
    IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
    IF ( ALLOCATED(inLev) ) DEALLOCATE(inLev)
    IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)

    ALLOCATE(INR8XYZT(nlon,nlat,nlev,nTime))
    ALLOCATE(INR8XYZ(nlon,nlat,nlev), INR8XY(nlon,nlat))

    !// Get dataset
    !call get_netcdf_var_3d( filename, 'CH4', inR8XYZ, nlon,nlat,nlev )
    call get_netcdf_var_4d( filename, 'CH4', inR8XYZT, nlon,nlat,nlev,nTime )

    INR8XYZ(:,:,:) = INR8XYZT(:,:,:,1)
    DEALLOCATE(INR8XYZT)

    !// Calculate grid area of dataset
    ALLOCATE(xbedge(nlon+1), ybedge(nlat+1), WGAULAT(nlat), WGAUWT(nlat), &
         WYEDG1P1(nLat+1), YEDG1P1(nLat+1), xybox(nlat+1))



    if (IPAR .ne. nLon) then
       !// Interpolate - must multiply by area before interpolating
       !// Find input edges
       !// Zonal
       dx = 360._r8/real(nLon,r8)
       XBEDGE(1) = -0.5*dx
       do J = 2, nLon+1
          XBEDGE(J) = XBEDGE(J-1) + dx
       end do
       !// Meridional (gaussian, as in p-grid.f)
       call GAUSST2(nLat, WGAULAT, WGAUWT)
       WYEDG1P1(1) = -1._r8
       do J = 2, nLat/2
          WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
       end do
       WYEDG1P1(nLat/2+1) = 0._r8
       do J = nLat/2+2, nLat+1
          WYEDG1P1(J) = -WYEDG1P1(nLat+2-J)
       end do
       !// Input Y-edged (degrees)
       do J = 1, nLat+1
          YBEDGE(J) = ZPI180 * asin(WYEDG1P1(J))
       end do

       !// Grid box areas
       do J = 1, nLat
          XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
               * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
       end do
       TOTAREA = sum(XYBOX) * real(nLon, r8)


       !// Grid box areas, multiply mmr-field with area
       hymnmean = 0._r8
       do J = 1, nLat
          do I = 1, nLon
             hymnmean = hymnmean + inR8XYZ(I,J,1) * XYBOX(J) / TOTAREA
          end do
       end do

       do L = 1, LPAR
          do J = 1, nLat
             do I = 1, nLon
                inR8XY(I,J) = inR8XYZ(I,J,L) * XYBOX(J)
             end do
          end do
          call E_GRID(INR8XY, XBEDGE, YBEDGE, nLon, nLat, &
               CTMR8, XDEDG, YDEDG, IPAR, JPAR,1)

          !// Finally divide by current area
          STT(:,:,L,trsp_idx(46)) = CTMR8(:,:)/AREAXY(:,:) &
                               * AIR(:,:,L) * TMOLMIX2MASSMIX(trsp_idx(46))
       end do

    else

       !// Resolution of data is the same as model run.
       !// Set STT (mass)
       hymnmean = sum(inR8XYZ(:,:,1) * AREAXY(:,:)) / sum(AREAXY)

       do L = 1, LPAR
         do J = 1, JPAR
           do I = 1, IPAR
             STT(I,J,L,trsp_idx(46)) = inR8XYZ(I,J,L)*AIR(I,J,L) &
                                       * TMOLMIX2MASSMIX(trsp_idx(46))
           end do
         end do
       end do

    end if

    write(6,'(a,2f12.2)') '* Overriding CH4 min/max (ppbv):',&
         minval(inR8XYZ(:,:,:))*1.e9_r8,maxval(inR8XYZ(:,:,:))*1.e9_r8

    if (meansfc .gt. 0._r8) then
       STT(:,:,:,trsp_idx(46)) = STT(:,:,:,trsp_idx(46)) * meansfc / hymnmean
       write(6,'(a,2es16.4)') '* Scaling 3D CH4 (ppbv) from/to:', &
            hymnmean, meansfc
    end if

    !// Deallocate
    deallocate (inR8XYZ, inR8XY, xbedge, ybedge, WGAULAT, WGAUWT, &
         WYEDG1P1, YEDG1P1, xybox)


    !// --------------------------------------------------------------------
  end subroutine set_ch4_stt
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine scale_ch4_stt()
    !// --------------------------------------------------------------------
    !// Scale up/down 3D CH4. Currently set up for using marine global
    !// average CH4, together with annual observations of ESRL Global
    !// Monitoring Division.
    !//
    !// Amund Sovde Haslerud, February 2017
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: STT, AIR, areaxy
    use cmn_chem, only: TMASSMIX2MOLMIX
    use cmn_sfc, only: LSMASK
    use cmn_oslo, only: CH4FIELD, trsp_idx, METHANEMIS
    use ncutils, only: GET_NETCDF_VAR_1D, GET_NETCDF_VAR_3D, GET_NETCDF_VAR_4D
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: I, J
    real(r8) :: sttmean, totarea
    !// Observed 1990: ESRL Global Monitoring Division marine global annual CH4.
    real(r8), parameter :: meansfc = 1682.08e-9_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='scale_ch4_stt'
    !// --------------------------------------------------------------------

    sttmean = 0._r8
    totarea = 0._r8

    do J = 1, JPAR
       do I = 1, IPAR
          !// Make a marine average
          if (LSMASK(I,J,1) .gt. 0.5_r8) then
             sttmean = sttmean + STT(I,J,1,trsp_idx(46)) / AIR(I,J,1) &
                  * AREAXY(I,J) * TMASSMIX2MOLMIX(trsp_idx(46))
             totarea = totarea + AREAXY(I,J)
          end if
       end do
    end do
    sttmean = sttmean / totarea

    STT(:,:,:,trsp_idx(46)) = STT(:,:,:,trsp_idx(46)) * meansfc / sttmean
    
    write(6,'(a,2es16.4)') f90file//':'//subr//': Scaling 3D CH4 (ppbv) from/to:', &
         sttmean, meansfc
    !// --------------------------------------------------------------------
  end subroutine scale_ch4_stt
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine updateSOILUPTAKEbousquet(LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Update soil uptake of CH4, from Bousquet data (2013, GAME project).
    !// The uptake is calculated from kg/h for the meteorological year,
    !// to 1/s by the use of surface values of CH4, taken from Oslo CTM2
    !// model simulation during the HYMN project (year 2003 values).
    !// This is an approximation, but we need uptake to be 1/s, not kg.
    !//
    !// Called from update_drydepvariables in drydeposition_oslo.f90.
    !//
    !// Amund Sovde, March 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, MPBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, JMON
    use cmn_met, only: MYEAR
    use cmn_parameters, only: AVOGNR
    use utilities, only: get_free_fileid
    use cmn_oslo, only: METHANEMIS, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    logical, intent(in) :: LNEW_MONTH

    !// Locals
    real(r8)  :: CTMXY(IPAR,JPAR), CH4SFC(IPAR,JPAR), uscale
    integer :: MP, I, J, II, JJ
    character(len=80) :: filename
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='updateSOILUPTAKEbousquet'
    !// --------------------------------------------------------------------

    !// No soil uptake if no CH4 emissions?
    !// May probably include it, but wait until later (needs testing).
    if (.not. METHANEMIS) then
       CH4SOILUPTAKE(:,:,:) = 0._r8
       return
    end if

    !// Skip if CH4 is not included
    if (trsp_idx(46) .le. 0) return

    !// Only update field every month
    if (.not.LNEW_MONTH) return

    !// This file was originally filed under /projects/geo/EMIS/GAME/
    filename = 'fch4.ref.mask11.nc'
    write(6,'(a)') f90file//':'//subr// &
         ': Updating CH4 soil uptake for VDEP: '//trim(filename)

    !// Read a certain month from netcdf file
    !// ---------------------------------------------------------
    !// Use current month (JMON) and meteorological year 2003 (matching
    !// HYMN year). Hence, this value will not change from year to year.
    !// Output values from this routine are POSITIVE for soil uptake.
    call read_ch4bousquet(filename, 'soils', JMON, 2003, CTMXY)


    !// Read CTM surface values of CH4 at surface. Use CTM3 runs where
    !// CH4 was set to HYMN 2003; data are given in kg.
    !// This is only to get velocities in 1/s, not kg/s, so there is no need
    !// to scale down or up for pre-industrial or future runs.
    call read_ch4sfc4soiluptake(CH4SFC)

    !// Scaling variables to scale up the loss rate
    uscale = 1._r8
    !// Unit scaling: 1/h to 1/s: 1/(3600)
    uscale = uscale / 3600._r8


    !// Set CH4SOILUPTAKE (units are kg/h, convert to 1/s)
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ   = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II   = I - MPBLKIB(MP) + 1
             CH4SOILUPTAKE(II,JJ,MP) = CTMXY(I,J) / CH4SFC(I,J) * uscale
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine updateSOILUPTAKEbousquet
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_ch4sfc4soiluptake(CH4SFC)
    !// --------------------------------------------------------------------
    !// Read surface CH4 from HYMN, to be used together with Bousquet
    !// surface uptake, to calculate uptake velocity for CH4.
    !//
    !// Amund Sovde, May 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR,JPAR
    use cmn_ctm, only: xdedg, ydedg, JMON, TMON, AREAXY
    use cmn_parameters, only: a0, CPI180, ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use ncutils, only: get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    real(r8), intent(out) :: CH4SFC(IPAR,JPAR)

    integer,parameter :: IT42=128, JT42=64
    integer :: J, I
    real(r8) :: INCH8(IT42,JT42,12),T42R8XY(IT42,JT42), dx
    real(r8) :: T42AREA, T42XDEDG(IT42+1), T42YDEDG(JT42+1)
    real(r8) :: WGAULAT(JT42),WGAUWT(JT42), WYEDG1P1(JT42+1)
    real(r8) :: YEDG1P1(JT42+1)

    character(len=100) :: FILENAME
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='read_ch4sfc4soiluptake'
    !// --------------------------------------------------------------------

    !// Read CH4 for 2003
    !// Data are on mixing ratio, T42 resolution
    FILENAME = 'Indata_CTM3/ch4_fields/ctm3_ch4sfc_2003.nc'

    !// Read the file (INCH8 must be R8)
    !// Month to get is '1', since there is only one month on each file.
    !call readnc_3d_from4d(filename,'lon',IT42,'lat',JT42,'lev',1, &
    !     'time', 1, 'CH4', INCH8)
    call get_netcdf_var_3d( filename, 'CH4SFC',INCH8, IT42,JT42,12)

    if (IPAR .ne. IT42) then
       !// Interpolate
       !// Find input edges
       !// Zonal
       dx = 360._r8/real(IT42,r8)
       T42XDEDG(1) = -0.5*dx
       do J=2,IT42+1
          T42XDEDG(J) = T42XDEDG(J-1) + dx
       end do

       !// Meridional (gaussian, as in p-grid.f)
       call GAUSST2(JT42,WGAULAT,WGAUWT)
       WYEDG1P1(1) = -1._r8
       do J = 2,JT42/2
          WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
       end do
       WYEDG1P1(JT42/2+1) = 0._r8
       do J = JT42/2+2,JT42+1
          WYEDG1P1(J) = -WYEDG1P1(JT42+2-J)
       end do
       !// T42 Y-edged (degrees)
       do J = 1,JT42+1
          T42YDEDG(J) = ZPI180 * asin(WYEDG1P1(J))
       end do

       !// Input is kg per grid box, so there is no need to multiply by
       !// area before interpolating
       T42R8XY(:,:) = INCH8(:,:,JMON) !// kg CH4

       Call E_GRID(T42R8XY,T42XDEDG,T42YDEDG,IT42,JT42, &
            CH4SFC,XDEDG,YDEDG,IPAR,JPAR,1)

    else

       !// Keep T42 field [kg]
       do J = 1, JPAR
          do I = 1, IPAR
             CH4SFC(I,J) = INCH8(I,J,JMON)
          end do
       end do

    end if

    !// --------------------------------------------------------------------
  end subroutine read_ch4sfc4soiluptake
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_ch4bousquet(filename,datasetname, mon,year, CTMXY)
    !// --------------------------------------------------------------------
    !// Reads a dataset 'soils' from Bousquet file of CH4 emissions and
    !// soil uptake.
    !// Routine changes negative signs to positive.
    !//
    !// Amund Sovde, March 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, JMON
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use emisutils_oslo, only: get_xyedges
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: mon, year
    character(len=*), intent(in) :: filename, datasetname
    !// Output
    real(r8), intent(out) :: CTMXY(IPAR,JPAR)

    !// Locals
    integer :: getY, getM, sTime, nTime, nLon, nLat, I, J
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    real(r8), allocatable, dimension(:) :: xbedge, ybedge, xybox
    real(r8), allocatable, dimension(:,:,:) :: allindata
    real(r8), allocatable, dimension(:,:) :: R8XY
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='read_ch4bousquet'
    !// --------------------------------------------------------------------

    !// Which year to fetch
    if (year .eq. 9999) then
       !// Use meteorological year
       getY = MYEAR
    else
       getY = year
    end if

    !// Which month?
    if (mon .eq. 99) then
       !// Use current month
       getM = JMON
    else
       getM = mon
    end if

    !// Get data from dataset
    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( filename, 'longitude',  inLon  )
    call get_netcdf_var_1d( filename, 'latitude',  inLat  )
    call get_netcdf_var_1d( filename, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    !// Select months to get: Find first entry for year getY
    sTime = -1
    do I=1,nTime
       if (inTime(I) .ge. real(getY)) then
          sTime = I
          exit
       end if
    end do
    if (sTime .le. 0) then
       write(6,'(a,i4,a)') f90file//':'//subr// &
            ': Cannot find year ',getY,' in dataset:'//trim(filename)
       write(6,'(a,i4,a)') '           sTime ',sTime
       stop
    end if

    sTime = sTime + getM - 1

    write(6,'(a,i5,i3,i5)') '* Reading bousquet soil deposition:',getY,getM,stime

    !// Test time
    if (sTime .gt. ntime) then
       write(6,'(a)') f90file//':'//subr//': Problem reading dataset'
       print*,'  stime/ntime:',sTime,nTime
       stop
    end if

    !// Emission variable, expected units: kg/m2/time for each month
    ALLOCATE(allindata(nlon,nlat,ntime), R8XY(nlon,nlat))
    call get_netcdf_var_3d( filename, datasetname, allindata, nlon,nlat,ntime )

    !// Fetch R8XY, switch signs from negative to positive
    R8XY(:,:) = -min(allindata(:,:,sTime),0._r8)

    !// Deallocate all local variables
    IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
    IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
    IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)
    IF ( ALLOCATED(allindata) ) DEALLOCATE(allindata)

    !// Calculate grid area of dataset
    ALLOCATE(xbedge(nlon), ybedge(nlat+1), xybox(nlat+1))
    !// Type 1 (starts at 180W,90S) for this dataset
    call get_xyedges(nLon,nLat,XBEDGE,YBEDGE,1)

    !// Grid box areas
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    !// Data is kg/m2/h, so we multiply with area before interpolation
    do J = 1, nLat
       do I = 1, nLon
          R8XY(I,J) = R8XY(I,J) * XYBOX(J)
       end do
    end do

    !// Interpolate
    call E_GRID(R8XY(:,:),XBEDGE,YBEDGE,nLon,nLat,CTMXY,XDEDG,YDEDG,IPAR,JPAR,1)

    !// Deallocate rest of variables
    DEALLOCATE(R8XY, xbedge, ybedge, xybox)

    print*,'* Done read_ch4bousquet'

    !// --------------------------------------------------------------------
  end subroutine read_ch4bousquet
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ch4drydep_bousquet(VDEP, BTT, DZ, MP)
    !// --------------------------------------------------------------------
    !// Sets up VDEP for CH4, based on Bousquet data.
    !// Applies only to Bousquet data, which has after read-in
    !// has units kg/s. The unit is here converted to m/s.
    !//
    !// Called from setdrydep in drydeposition_oslo.f90
    !//
    !// Amund Sovde, March 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IDBLK, JDBLK, IPAR, JPAR, LPAR, NPAR
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_oslo, only: trsp_idx, METHANEMIS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8),  intent(in) :: BTT(LPAR,NPAR,IDBLK,JDBLK)
    real(r8),  intent(in) :: DZ(IDBLK,JDBLK)
    integer, intent(in) :: MP

    !// Output
    real(r8), intent(inout) :: VDEP(NPAR,IPAR,JPAR)

    !// Locals
    integer :: NTR, I,II,J,JJ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='ch4drydep_bousquet'
    !// --------------------------------------------------------------------

    !// If no CH4 emissions, skip setting VDEP for CH4. It has already been
    !// initialised to zero, so we return:

    !// You may want to comment out this if your purpose is to create new
    !// scaling factors; then you may want the soil uptake as part of the loss.
    if (.not. METHANEMIS) return

    !// Transport number of CH4
    NTR = trsp_idx(46)
    if (NTR .le. 0) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Need to change 1/s to m/s
          !// DZ was calculated in calling routine
          VDEP(NTR,I,J) = CH4SOILUPTAKE(II,JJ,MP) * DZ(II,JJ)

       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine ch4drydep_bousquet
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine reportsfcch4()
    !// --------------------------------------------------------------------
    !// Print out global average surface mixing ratio of CH4 (ppb) and
    !// total CH4 burden (Tg).
    !//
    !// Amund Sovde, May 2014
    !// --------------------------------------------------------------------
    use cmn_ctm, only: STT, AIR, AREAXY
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8) :: total
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='reportsfcch4'
    !// --------------------------------------------------------------------
    if (trsp_idx(46) .le. 0) return
    total = sum(STT(:,:,1,trsp_idx(46))/AIR(:,:,1)*AREAXY(:,:)) &
           /sum(AREAXY)*28.97e9_r8/16._r8
    write(6,'(a,2f9.2)') '* Average CH4@sfc[ppb] & burden[Tg]: ', &
         total, sum(STT(:,:,:,trsp_idx(46)))*1.e-9_r8
    !// --------------------------------------------------------------------
  end subroutine reportsfcch4
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
end module ch4routines
!//=========================================================================

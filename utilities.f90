!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Amund Sovde Haslerud, August 2017
!//=========================================================================
!// Utilities for running the CTM.
!//=========================================================================
module utilities
  !//-----------------------------------------------------------------------
  !// MODULE: utilities
  !// DESCRIPTION: Utilities for running the CTM.
  !//
  !// Contains:
  !//   subroutine write_log
  !//   subroutine model_info
  !//   subroutine calendar
  !//   subroutine is_leap
  !//   subroutine get_soldecdis
  !//   subroutine CALENDR_OLD
  !//   subroutine CALENDL
  !//   subroutine LOCSZA
  !//   subroutine LCM
  !//   subroutine ctmExitC
  !//   subroutine ctmExitL
  !//   subroutine ctmExitIJL
  !//   integer function get_free_fileid
  !//   subroutine get_dinm
  !//   subroutine check_btt
  !//   subroutine adjust_moments
  !//   real(r8) function moninobukhov_length 
  !//
  !// Stefanie Falk, May 2019
  !// Amund Sovde Haslerud, August 2017
  !//   Split CALENDR into calendar and get_soldecdis.
  !// Ole Amund Sovde, March 2015
  !//   To f90 based on p-utils.f.
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !//-----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'utilities.f90'
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------


contains


  !//-----------------------------------------------------------------------
  subroutine write_log(ISW, START_TIME, END_TIME, TOT_TIME)
    !//---------------------------------------------------------------------
    !// Write info about start/end of run
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPARW,JPARW,LPARW, IPAR,JPAR,LPAR, &
         MODEL, MODEL_VERSION
    use cmn_ctm, only: JDATE, TMON, JYEAR
    use cmn_oslo, only: DINM
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ISW, START_TIME(8), END_TIME(8)
    real(r8), intent(in) :: TOT_TIME
    !// Locals
    integer :: MM_USED, I
    character(len=80) :: BAR
    !//---------------------------------------------------------------------

    !// Make a bar of '>>>>'
    do I = 0, 19
       BAR((i*4+1):(i+1)*4) = '>>>>'
    end do

    !// ISW=0: start
    if (ISW .eq. 0) then
       !//---version title & resolution
       write(6,'(a)') BAR
       write(6,'(a)') MODEL
       write(6,'(a,1x,a)') ' version: '//trim(MODEL_VERSION), &
            '(based on UCI CTM qcode 7.1 1/2015)'
       write(6,'(a)') BAR
       write(6,'(A,5I4)') 'run starts at YR/Mon/DD/HH/MM/ ',&
            START_TIME(1), START_TIME(2), START_TIME(3), &
            START_TIME(5), START_TIME(6)
       write(6,'(a)') BAR
       !// Print out more specific output.
       call model_info()
       write(6,'(a)') BAR
    end if

    !//ISW=1: end of run
    if (ISW .eq. 1) then
       write(6,'(a)') BAR
       write(6,'(a)') MODEL
       write(6,'(a)') ' version: '//trim(MODEL_VERSION)
       write(6,'(a)') BAR
       write(6,'(A,5I4)') 'run ends at YR/Mon/DD/HH/MM/ ',&
            END_TIME(1), END_TIME(2), END_TIME(3), &
            END_TIME(5), END_TIME(6)


       !// Minutes at end minus minutes at start
       MM_USED = nint(TOT_TIME / 60._r8)

       write(6,'(A,I9)') 'total run time (min)', MM_USED
       write(6,'(A,I6,2x,A3,I5)') 'end of simulation: ',JDATE,TMON,JYEAR
       write(6,'(a)') BAR

    end if

    !//---------------------------------------------------------------------
  end subroutine write_log
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine model_info()
    !//---------------------------------------------------------------------
    !// Write info about model run.
    !//---------------------------------------------------------------------
    use cmn_size, only: IPARW,JPARW,LPARW, IPAR,JPAR,LPAR, MPBLK, &
         MODEL, MODEL_VERSION, &
         LOSLOCHEM, LOSLOCTROP, LOSLOCSTRAT, LSULPHUR, LBCOC, LSALT, &
         LDUST, LNITRATE, LSOA, LEMISDEP_INCHEM, LE90, LLINOZ, &
         NPAR, NPAR_TROP, NPAR_STRAT, NPAR_SUL, NPAR_NITRATE, &
         NPAR_SALT, NPAR_DUST, NPAR_BC, NPAR_OM, NPAR_E90, NPAR_LINOZ, &
         NPAR_SOA
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer :: IOMP_THREADS, omp_get_num_threads
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'model_info'
    !//---------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': model specifics:'
    write(6,'(a,3I5,5X,3I5)') ' Native resolution:      ', IPARW,JPARW,LPARW
    write(6,'(a,3I5,5X,3I5)') ' Model resolution:       ', IPAR,JPAR,LPAR

    write(6,'(a)') ' Included modules:'
    if (LOSLOCHEM) then
       write(6,'(a,l1)') '   Tropospheric chemistry:           ', LOSLOCTROP
       write(6,'(a,l1)') '   Stratospheric chemistry:          ', LOSLOCSTRAT
       write(6,'(a,l1)') '   Sulphur module:                   ', LSULPHUR
       write(6,'(a,l1)') '   Black/organic carbon (BCOC):      ', LBCOC
       write(6,'(a,l1)') '   Sea salt (SALT):                  ', LSALT
       write(6,'(a,l1)') '   Nitrate module:                   ', LNITRATE
       write(6,'(a,l1)') '   Mineral dust module:              ', LDUST
       write(6,'(a,l1)') '   Secondary organic aerosols (SOA): ', LSOA
       write(6,'(a,l1)') '   Emis & dep inside chemistry?:     ', LEMISDEP_INCHEM
    else
       write(6,'(a,l1)') '   Oslo chemistry:                   ', LOSLOCHEM
    end if

    write(6,'(a,i3)')    ' Total # tracers (NPAR): ', NPAR
    if (LOSLOCTROP)  write(6,'(10x,a,i3)')    'NPAR_TROP:     ', NPAR_TROP
    if (LOSLOCSTRAT) write(6,'(10x,a,i3)')    'NPAR_STRAT:    ', NPAR_STRAT
    if (LSULPHUR)    write(6,'(10x,a,i3)')    'NPAR_SUL:      ', NPAR_SUL
    if (LNITRATE)    write(6,'(10x,a,i3)')    'NPAR_NITRATE:  ', NPAR_NITRATE
    if (LSALT)       write(6,'(10x,a,i3)')    'NPAR_SALT:     ', NPAR_SALT
    if (LDUST)       write(6,'(10x,a,i3)')    'NPAR_DUST:     ', NPAR_DUST
    if (LBCOC)       write(6,'(10x,a,i3)')    'NPAR_BC:       ', NPAR_BC
    if (LBCOC)       write(6,'(10x,a,i3)')    'NPAR_OM:       ', NPAR_OM
    if (LSOA)        write(6,'(10x,a,i3)')    'NPAR_SOA:      ', NPAR_SOA
    if (LE90)        write(6,'(10x,a,i3)')    'NPAR_E90:      ', NPAR_E90
    if (LLINOZ)      write(6,'(10x,a,i3)')    'NPAR_LINOZ:    ', NPAR_LINOZ

    !// OpenMP information about threads
    IOMP_THREADS = 1 !// Initialize
!$omp parallel
!$omp master
!$      IOMP_THREADS = omp_get_num_threads()
!$omp end master
!$omp end parallel
    write(6,'(a,i3)')    ' Number of CPUs (OMP_NUM_THREADS): ', IOMP_THREADS
    write(6,'(a,i7)')    ' Number of IJ parallel blocks: ', MPBLK

    !//---------------------------------------------------------------------
  end subroutine model_info
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine calendar(IYEAR,IDAY, LLPYR,LFIXMET,MYEAR, &
                      JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET, &
                      JYEAR_NEXT,JDAY_NEXT,JMON_NEXT,JDATE_NEXT)
    !//---------------------------------------------------------------------
    !//--- compute the year/month/day counting from day# IDAY of year IYEAR
    !//--- note that IYEAR is the reference year and IDAY >>365 for long runs
    !//---     JYEAR = current year
    !//---     JDAY = day of year (1:365 or 366)
    !//---     LLPYR = allow for leap year
    !//---     LFIXMET = recycle met fields(no leap yr)
    !//---     JMON = no. of month (1:12)
    !//---     TMON = 3-char label of month
    !//---     JDATE = day of the month (1:31)
    !//---     LYEAR = .true. = current year is being treated as a leap year
    !//---     TMET = 3-char label for day of year (Feb 29 = '901')
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input/output
    integer, intent(in)      :: IYEAR,IDAY
    logical, intent(in)      :: LLPYR,LFIXMET
    integer, intent(out)     :: JYEAR,JDAY,JMON,JDATE,MYEAR
    integer, intent(out)     :: JYEAR_NEXT,JDAY_NEXT,JMON_NEXT,JDATE_NEXT
    logical, intent(out)     :: LYEAR
    character(len=3), intent(out) :: TMON,TMET

    !// Locals
    integer :: KDAY,KYEAR, JMET,I, DIY
    real(r8) :: G,SOLLNG,SINDEC
    logical :: TMPLOG

    !// Parameters
    !// Month of day
    integer, parameter, dimension(365) ::  JMOFD = [(1,I=1,31), &
         (2,I=32,59),     (3,I=60,90),    (4,I=91,120), (5,I=121,151), &
         (6,I=152,181),  (7,I=182,212),  (8,I=213,243), (9,I=244,273), &
         (10,I=274,304), (11,I=305,334), (12,I=335,365) ]

    integer, parameter, dimension(366) ::  JMOFDL= [(1,I=1,31), &
         (2,I=32,60),     (3,I=61,91),    (4,I=92,121), (5,I=122,152), &
         (6,I=153,182),  (7,I=183,213),  (8,I=214,244), (9,I=245,274), &
         (10,I=275,305), (11,I=306,335), (12,I=336,366) ]

    !// Number of days in previous months
    integer, parameter, dimension(13)  :: JDOFM = &
         [0,31,59,90,120,151,181,212,243,273,304,334,365]
    integer, parameter, dimension(13)  :: JDOFML = &
         [0,31,60,91,121,152,182,213,244,274,305,335,366]

    character(len=3), parameter, dimension(12) :: AMON = ['JAN','FEB', &
         'MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    !//---------------------------------------------------------------------

    KDAY = IDAY   !// Current IDAY
    KYEAR = IYEAR !// This is reference year

    !// Leap year?
    I = 0
    do while (I .eq. 0)

       !// New leap year calculation includes the special century years:
       !// Note that now you can run LFIXMET with LLPYR.
       if (LLPYR) then
          call is_leap(KYEAR, LYEAR)
       else
          !// Do not take leap years into account
          LYEAR = .false.
       end if

       if (LYEAR) then
          DIY = 366
       else
          DIY = 365
       end if

       if (KDAY .gt. DIY) then
          KDAY = KDAY - DIY
          KYEAR = KYEAR + 1
       else
          I = 1 !// Done with iteration; got to current year
       end if
    end do


    !//---have reached current year
    JYEAR = KYEAR
    JDAY = KDAY
    if (LFIXMET)  then
       MYEAR = IYEAR
    else
       MYEAR = JYEAR
    end if
    if (LYEAR) then
       JMON = JMOFDL(JDAY)
       JDATE = JDAY - JDOFML(JMON)
       if (JDAY .gt. 60) then
          JMET = JDAY - 1
       else if (JDAY .eq. 60) then
          JMET = 901
       else
          JMET = JDAY
       end if
    else
       JMON = JMOFD(JDAY)
       JDATE = JDAY - JDOFM(JMON)
       JMET = JDAY
    end if
    TMON = AMON(JMON)

    write (TMET(1:3),'(I3.3)')  JMET


    !// Calculate for next day + 1
    !//---------------------------------------------------------------------
    !// See above for comments
    KDAY = IDAY + 1
    KYEAR = IYEAR

    I = 0
    do while (I .eq. 0)
       if (LLPYR) then
          call is_leap(KYEAR, TMPLOG)
       else
          TMPLOG = .false.
       end if

       if (TMPLOG) then
          DIY = 366
       else
          DIY = 365
       end if

       if (KDAY .gt. DIY) then
          KDAY = KDAY - DIY
          KYEAR = KYEAR + 1
       else
          I = 1 !// Done
       end if
    end do


    !//---have reached current year
    JYEAR_NEXT = KYEAR
    JDAY_NEXT = KDAY
    if (TMPLOG) then
       JMON_NEXT = JMOFDL(JDAY_NEXT)
       JDATE_NEXT = JDAY_NEXT - JDOFML(JMON_NEXT)
    else
       JMON_NEXT = JMOFD(JDAY_NEXT)
       JDATE_NEXT = JDAY_NEXT - JDOFM(JMON_NEXT)
    end if

    write(6,'(A,I6,2(i4,i3,i3,i5))') 'calendar: IDAY/DAY/NEXT: ', IDAY, &
         JDAY,JDATE,JMON,JYEAR, JDAY_NEXT,JDATE_NEXT,JMON_NEXT,JYEAR_NEXT
    !//---------------------------------------------------------------------
  end subroutine calendar
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine is_leap(KYEAR,LYEAR)
    !//---------------------------------------------------------------------
    !// Calculates whether KYEAR is a leap year.
    !//
    !// Amund Sovde Haslerud, January 2017
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: KYEAR
    logical, intent(out) :: LYEAR
    !//---------------------------------------------------------------------

    !// New leap year calculation includes the special century years:
    if (mod(KYEAR,4) .ne. 0) then
       LYEAR = .false.
    else if (mod(KYEAR,100) .ne. 0) then
       LYEAR = .true.
    else if (mod(KYEAR,400) .ne. 0) then
       LYEAR = .false.
    else
       LYEAR = .true.
    end if

    !//---------------------------------------------------------------------
  end subroutine is_leap
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  subroutine get_soldecdis(JDAY,SOLDEC,SOLDIS)
    !//---------------------------------------------------------------------
    !// Calculate from JDAY the solar declination and distance to the Sun.
    !// Taken from old CALENDR routine.
    !//
    !// Amund Sovde Haslerud, January 2017
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: CPI180, ZPI180
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: JDAY
    real(r8),  intent(out)     :: SOLDEC,SOLDIS

    real(r8) :: G, SOLLNG, SINDEC

    !//--------------sun-earth data from US Naval Obs-------------------------
    !//---small shifts from year-to-year, just use year 2000
    !//---  e(deg) = 23.439 - 0.00013 * (Y - 2000)  (obliquity)
    !//---Earth-sun distance astronomical units (AU)
    !//---  g(deg) = 356.543 + 0.98560028*JDAY (deg, JDAY=day of year=1=Jan 1.5)
    !//---  R = 1.00014 - 0.01671 cos g - 0.00014 cos 2g
    !//---  q(deg) = 279.473 + 0.98564736*JDAY
    !//---  L = q + 1.915 sin g + 0.020 sin 2g  (L= apparent eclipt long of Sun)
    !//--- sin d = sin e * sin L (solar declination)
    !real(r8), parameter :: CPI180 = 0.01745329252_r8
    real(r8), parameter :: COBLIQ = 0.3977725_r8
    !real(r8), parameter :: ZPI180 = 1._r8/CPI180
    !//---------------------------------------------------------------------
    !//---------------------------------------------------------------------

    !//---solar declination and distance to sun: from US Naval Obs
    G = CPI180*(356.543_r8 + 0.98560028_r8 * real(JDAY, r8))
    SOLDIS = 1.00014_r8 - 0.01671_r8*cos(G) - 0.00014_r8*cos(G+G)
    SOLLNG = 279.473_r8 + 0.98564736_r8 * real(JDAY, r8) &
                       + 1.915_r8*sin(G) + 0.02_r8*sin(G+G)
    SINDEC = COBLIQ * sin(SOLLNG*CPI180)
    SOLDEC = asin(SINDEC)*ZPI180

    !//---------------------------------------------------------------------
  end subroutine get_soldecdis
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine CALENDR_OLD (IYEAR,IDAY, LLPYR,LFIXMET,MYEAR, &
                      JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
    !//---------------------------------------------------------------------
    !//--- compute the year/month/day counting from day# IDAY of year IYEAR
    !//--- note that IYEAR is the reference year and IDAY >>365 for long runs
    !//---     JYEAR = current year
    !//---     JDAY = day of year (1:365 or 366)
    !//---     LLPYR = allow for leap year
    !//---     LFIXMET = recycle met fields(no leap yr)
    !//---     JMON = no. of month (1:12)
    !//---     TMON = 3-char label of month
    !//---     JDATE = day of the month (1:31)
    !//---     LYEAR = .true. = current year is being treated as a leap year
    !//---     TMET = 3-char label for day of year (Feb 29 = '901')
    !//---     SOLDEC = solar declination (degrees)
    !//---     SINDEC = sin(solar declination)
    !//---     SOLDIS = distance to sun (A.U.)
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input/output
    integer, intent(in)      :: IYEAR,IDAY
    logical, intent(in)      :: LLPYR,LFIXMET
    integer, intent(out)     :: JYEAR,JDAY,JMON,JDATE,MYEAR
    logical, intent(out)     :: LYEAR
    real(r8),  intent(out)     :: SOLDEC,SOLDIS
    character(len=3), intent(out) :: TMON,TMET

    !// Locals
    integer :: KDAY,KYEAR, JMET,I, DIY
    real(r8) :: DFYR, G,SOLLNG,SINDEC

    !// Parameters
    integer, parameter, dimension(365) ::  JMOFD = [(1,I=1,31), &
         (2,I=32,59),     (3,I=60,90),    (4,I=91,120), (5,I=121,151), &
         (6,I=152,181),  (7,I=182,212),  (8,I=213,243), (9,I=244,273), &
         (10,I=274,304), (11,I=305,334), (12,I=335,365) ]

    integer, parameter, dimension(366) ::  JMOFDL= [(1,I=1,31), &
         (2,I=32,60),     (3,I=61,91),    (4,I=92,121), (5,I=122,152), &
         (6,I=153,182),  (7,I=183,213),  (8,I=214,244), (9,I=245,274), &
         (10,I=275,305), (11,I=306,335), (12,I=336,366) ]

    integer, parameter, dimension(13)  :: JDOFM = &
         [0,31,59,90,120,151,181,212,243,273,304,334,365]
    integer, parameter, dimension(13)  :: JDOFML = &
         [0,31,60,91,121,152,182,213,244,274,305,335,366]

    character(len=3), parameter, dimension(12) :: AMON = ['JAN','FEB', &
         'MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

    !//--------------sun-earth data from US Naval Obs-------------------------
    !//---small shifts from year-to-year, just use year 2000
    !//---  e(deg) = 23.439 - 0.00013 * (Y - 2000)  (obliquity)
    !//---Earth-sun distance astronomical units (AU)
    !//---  g(deg) = 356.543 + 0.98560028*JDAY (deg, JDAY=day of year=1=Jan 1.5)
    !//---  R = 1.00014 - 0.01671 cos g - 0.00014 cos 2g
    !//---  q(deg) = 279.473 + 0.98564736*JDAY
    !//---  L = q + 1.915 sin g + 0.020 sin 2g  (L= apparent eclipt long of Sun)
    !//--- sin d = sin e * sin L (solar declination)
    real(r8), parameter :: CPI180 = 0.01745329252_r8
    real(r8), parameter :: COBLIQ = 0.3977725_r8
    real(r8), parameter :: ZPI180 = 1._r8/CPI180
    !//---------------------------------------------------------------------

    KDAY = IDAY   !// Current IDAY
    KYEAR = IYEAR !// This is reference year

    !// Leap year?
    I = 0
    do while (I .eq. 0)
       !// Old leap year only valid 1904 - 2096 because assumed to be
       !// divisible by 4.
       !LYEAR = (mod(KYEAR,4).eq.0) .and. LLPYR .and. .not.LFIXMET

       !// New leap year calculation includes the special century years:
       !// Note that now you can run LFIXMET with LLPYR.
       if (LLPYR) then
          if (mod(KYEAR,4) .ne. 0) then
             LYEAR = .false.
          else if (mod(KYEAR,100) .ne. 0) then
             LYEAR = .true.
          else if (mod(KYEAR,400) .ne. 0) then
             LYEAR = .false.
          else
             LYEAR = .true.
          end if
       else
          !// Do not take leap years into account
          LYEAR = .false.
       end if

       if (LYEAR) then
          DIY = 366
       else
          DIY = 365
       end if
       if (KDAY .gt. DIY) then
          KDAY = KDAY - DIY
          KYEAR = KYEAR + 1
       else
          I = 1 !// Done with iteration; got to current year
       end if
    end do


    !//---have reached current year
    JYEAR = KYEAR
    JDAY = KDAY
    if (LFIXMET)  then
       MYEAR = IYEAR
    else
       MYEAR = JYEAR
    end if
    if (LYEAR) then
       JMON = JMOFDL(JDAY)
       JDATE = JDAY - JDOFML(JMON)
       DFYR = real(JDAY, r8) / 366._r8
       if (JDAY .gt. 60) then
          JMET = JDAY - 1
       else if (JDAY .eq. 60) then
          JMET = 901
       else
          JMET = JDAY
       end if
    else
       JMON = JMOFD(JDAY)
       JDATE = JDAY - JDOFM(JMON)
       DFYR = real(JDAY, r8) / 365._r8
       JMET = JDAY
    end if
    TMON = AMON(JMON)

    write (TMET(1:3),'(I3.3)')  JMET

    !//---solar declination and distance to sun: from US Naval Obs
    G = CPI180*(356.543_r8 + 0.98560028_r8 * real(JDAY, r8))
    SOLDIS = 1.00014_r8 - 0.01671_r8*cos(G) - 0.00014_r8*cos(G+G)
    SOLLNG = 279.473_r8 + 0.98564736_r8 * real(JDAY, r8) &
                       + 1.915_r8*sin(G) + 0.02_r8*sin(G+G)
    SINDEC = COBLIQ * sin(SOLLNG*CPI180)
    SOLDEC = asin(SINDEC)*ZPI180

    write(6,'(A,I6)') 'CALENDR: IDAY: ', IDAY
    !//---------------------------------------------------------------------
  end subroutine CALENDR_OLD
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine CALENDL (JYEAR,JDAY,LYEAR, NDAY366, LD,LD2)
    !//---------------------------------------------------------------------
    !//--- look up calendar array NDAYD for 0 or 1 and set LD = false or true
    !//---     JYEAR = current year
    !//---     JDAY = day of year (1:365 or 1:366)
    !//---     LYEAR = .true. = current year is being treated as a leap year
    !//---     NDAY366(1:366) = input cal. for triggering diagnostics or saves
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input/output
    integer, intent(in)     :: JYEAR,JDAY, NDAY366(366)
    logical, intent(in)     :: LYEAR
    logical, intent(out)    :: LD,LD2
    !// Local
    integer :: I
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CALENDL'
    !//---------------------------------------------------------------------

    if (JYEAR .lt. 1900) then
       !// Probably due to leap year
       write(6,'(a)') f90file//':'//subr// &
            ': Calendar problems with JYEAR < 1900'
       stop 'STOP in '//subr
    end if
    if (LYEAR) then
       I = JDAY
    else
       if (JDAY .gt. 59) then
          I = JDAY + 1
       else
          I = JDAY
       end if
    end if
    !//---use NDAY366 for triggering two logicals (e.g., unf & stdout writes)
    LD  = NDAY366(I) .gt. 0
    LD2 = NDAY366(I) .gt. 1

    !//---------------------------------------------------------------------
  end subroutine CALENDL
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine LOCSZA (UTT,XDEG,YDEG,SOLDEC,SOLDIS, COSSZA,SOLFX)
    !//---------------------------------------------------------------------
    !//---calculates COSSZA = cosine of local Solar Zenith Angle
    !//              SOLFX  = solar flux factor (from distance to sun)
    !//---for:   UTT  = Universal Time (in hours, UTT=0 = midnight GMT)
    !//          XDEG = longitude (degrees)
    !//          YDEG = latitude (degrees)
    !//---assumes CALENDR called, giving:
    !//           SOLDEC = solar declination (degrees)
    !//           SOLDIS = distance to sun (A.U.)
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input/output
    real(r8), intent(in)  :: UTT,XDEG,YDEG,SOLDEC,SOLDIS
    real(r8), intent(out) :: COSSZA,SOLFX
    !// Locals
    real(r8)  SINDEC, COSDEC, SINLAT, COSLAT, LST 
    real(r8),parameter :: CPI180 = 3.141592653589793_r8 /180._r8
    real(r8),parameter :: CPIHRS = 3.141592653589793_r8 / 12._r8
    !//---------------------------------------------------------------------

    !//---solar declination
    SINDEC = sin(SOLDEC*CPI180)
    COSDEC = cos(SOLDEC*CPI180)
    !//---latitude
    SINLAT = sin(YDEG*CPI180)
    COSLAT = cos(YDEG*CPI180)
        
    !//---local solar time (hr), defined as 0.00 at local noon
    LST = (UTT - 12._r8) + XDEG*0.066666666_r8
    !//---cos of local solar zenith angle
    COSSZA = COSDEC*COSLAT*cos(LST*CPIHRS) + SINDEC*SINLAT
    !//---solar flux factor (distance to sun, 0.034 is avg for +-.0167)
    SOLFX  = SOLDIS**(-2)

    !//---------------------------------------------------------------------
  end subroutine LOCSZA
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine LCM (NN1,NN2,LCM12)
    !//---------------------------------------------------------------------
    !// NLCM = least common multiple of N1 & N2 by Euclidean algorithm for GCD
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input/output
    integer, intent(in)  :: NN1, NN2
    integer, intent(out) :: LCM12
    !// Locals
    integer :: N1, N2, N3, GCD
    !//---------------------------------------------------------------------

    N1 = NN1
    N2 = NN2
    N3 = mod(N1, N2)
    do while (N3.ne.0)
       N1 = N2
       N2 = N3
       N3 = mod(N1, N2)
    end do
    GCD = N2
    LCM12 = NN1*NN2/GCD

    !//---------------------------------------------------------------------
  end subroutine LCM
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine ctmExitC(MESSAG)
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=*),intent(in)::   MESSAG
    !//---------------------------------------------------------------------
    write(6,'(a)')  MESSAG
    stop
    !//---------------------------------------------------------------------
  end subroutine ctmExitC
  !//-----------------------------------------------------------------------


   
  !//-----------------------------------------------------------------------
  subroutine ctmExitL(MESSAG,LABEL)
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=*),intent(in)::   MESSAG, LABEL
    !//---------------------------------------------------------------------
    write(6,'(2a)')  MESSAG,LABEL
    stop
    !//---------------------------------------------------------------------
  end subroutine ctmExitL
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine ctmExitIJL(MESSAG,I,J,L)
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=*),intent(in)::   MESSAG
    integer, intent(in):: I,J,L
    !//---------------------------------------------------------------------
    write(6,'(a)')  MESSAG
    write(6,'(3i10)') I,J,L
    stop
    !//---------------------------------------------------------------------
  end subroutine ctmExitIJL
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  integer function get_free_fileid()
    !// --------------------------------------------------------------------
    !// Get a free file number to use for opening files.
    !//
    !// Ole Amund Sovde, July 2011
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: file_nr
    logical :: file_io
    !// --------------------------------------------------------------------
    file_nr = 20      !// Starting number
    file_io = .true.  !// Assume file is opened
    do while (file_io)
       file_nr = file_nr + 1
       inquire(file_nr,opened=file_io)
    end do
    !// file_nr is not opened and can be used
    get_free_fileid = file_nr
    !// --------------------------------------------------------------------
  end function get_free_fileid
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_dinm(LYEAR)
    !// --------------------------------------------------------------------
    !// Set up number of days in month for regular year or leap year.
    !//
    !// Ole Amund Sovde, November 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    logical, intent(in) :: LYEAR

    !// Parameters
    integer, dimension(12), parameter :: daysinmonth = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    integer, dimension(12), parameter :: daysinmonth_leap = &
         (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    !// --------------------------------------------------------------------

    if (LYEAR) then
       DINM(:) = daysinmonth_leap(:)
    else
       DINM(:) = daysinmonth(:)
    end if

    !// --------------------------------------------------------------------
  end subroutine get_dinm
  !// ----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine CFRMIN(T, CFR, CW, CI, EPS)
    !//---------------------------------------------------------------------
    !// Put the minimum limit to cloud fraction.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in)    :: T
    real(r8), intent(in)    :: EPS
    real(r8), intent(in)    :: CW, CI ! water/ice ratio averaged over box
    !// Input/Output
    real(r8), intent(inout) :: CFR    ! cloud fraction
    !// Parameters
    real(r8), parameter  :: TIW   = 245._r8  ! Liquid above, ice below
    real(r8), parameter  :: CWMIN = 1.e-5_r8 ! in-cloud water ratio limit
    real(r8), parameter  :: CIMIN = 1.e-6_r8 ! in-cloud ice ratio limit
    !//---------------------------------------------------------------------

    !// Remember that this routine is called when CFR<EPS,
    !// it does not limit any cloud fractions.
    if (T .gt. TIW .and. CW .gt. (CWMIN * EPS)) then
       !// There is liquid cloud
       CFR = EPS
    else if (T .le. TIW .and. CI .gt. (CIMIN * EPS)) then
       !// There is ice cloud
       CFR = EPS
    else
       !// There is no cloud
       CFR = 0._r8
    end if

    !//---------------------------------------------------------------------
  end subroutine CFRMIN
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine CIWMIN(T, CW, CI)
    !//---------------------------------------------------------------------
    !// Put the minimum limits to in-cloud water and ice ratios (kg/kg).
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in)    :: T
    !// Input/Output
    real(r8), intent(inout) :: CW, CI ! in-cloud water/ice ratio
    !// Parameters
    real(r8), parameter  :: TIW   = 245._r8
    real(r8), parameter  :: CWMIN = 1.e-5_r8 ! in-cloud water ratio limit
    real(r8), parameter  :: CIMIN = 1.e-6_r8 ! in-cloud ice ratio limit
    !//---------------------------------------------------------------------

    if (T .gt. TIW) then
       CW = max(CW, CWMIN)
    else
       CI = max(CI, CIMIN)
    end if

    !//---------------------------------------------------------------------
  end subroutine CIWMIN
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine check_btt(BTT,MP,MESSAGE)
    !// --------------------------------------------------------------------
    !// Putting BTT (transported tracers) and XSTT (non-transported tracers)
    !// into ZC_LOCAL, to avoid striding in chemistry integration routine.
    !// ONLY moves the column between L_START and L_END, so that this
    !// routine is applicable to both the troposphere and the stratosphere.
    !//
    !// Can consider only putting tropospheric components into ZC_LOCAL for
    !// tropospheric chemistry (stratospheric chemistry will also need some
    !// tropospheric components), but I doubut it will be any faster.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TNAME
    use cmn_oslo, only: LMTROP, chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    character(len=*), intent(in) :: MESSAGE

    !// Locals
    integer :: TRACER_ID, N, L
    integer :: I,J, II,JJ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'check_btt'
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        do N = 1, NPAR
          !// Generating vertical arrays
          do L = 1, LPAR
            if (BTT(L,N,II,JJ) .lt. 0._r8) then
              !// Only print out for other species than SO, which is
              !// allowed to be negative.
              if (chem_idx(N).ne.110) then
                write(6,'(a,i4,es13.5,i3,6i4,i3)') f90file//':'//subr// &
                     ': NEGATIVE '//trim(TNAME(N)), chem_idx(N), &
                     BTT(L,N,II,JJ),L,N,I,J,II,JJ,MP, LMTROP(I,J)
                write(6,'(a)') MESSAGE
                stop 'STOP in '//subr
              end if
            end if
            if (BTT(L,N,II,JJ) .ne. BTT(L,N,II,JJ)) then
              !// NANs!
               write(6,'(a,i4,es13.5,i3,6i4,i3)') f90file//':'//subr// &
                    ': NAN '//trim(TNAME(N)), chem_idx(N), &
                   BTT(L,N,II,JJ),L,N,I,J,II,JJ,MP, LMTROP(I,J)
              write(6,'(a)') MESSAGE
              stop 'STOP in '//subr
            end if
          end do
        end do
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine check_btt
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine check_stt(MESSAGE)
    !// --------------------------------------------------------------------
    !// Check STT
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, NPAR, IPAR, JPAR
    use cmn_ctm, only: STT
    use cmn_chem, only: TNAME
    use cmn_oslo, only: LMTROP, chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: MESSAGE

    !// Locals
    integer :: TRACER_ID, N, L
    integer :: I,J 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'check_stt'
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do N = 1, NPAR
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (STT(I,J,L,N) .lt. 0._r8) then
                   if (chem_idx(N).ne.110) then
                      write(6,'(a,i4,es13.5,5i4)') f90file//':'//subr// &
                           ': NEGATIVE '//trim(TNAME(N)), chem_idx(N), &
                           STT(I,J,L,N),L,N,I,J, LMTROP(I,J)
                      write(6,'(a)') MESSAGE
                      stop 'STOP in '//subr
                   end if
                end if
                if (STT(I,J,L,N) .ne. STT(I,J,L,N)) then
                   !// NANs!
                   write(6,'(a,i4,es13.5,5i4)') f90file//':'//subr// &
                        ': NAN '//trim(TNAME(N)), chem_idx(N), &
                        STT(I,J,L,N),L,N,I,J, LMTROP(I,J)
                   write(6,'(a)') MESSAGE
                   stop 'STOP in '//subr
                end if
             end do
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine check_stt
  !// ----------------------------------------------------------------------





  !// ----------------------------------------------------------------------
  subroutine adjust_moments(BTT,BTTBCK,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,MP)
    !// --------------------------------------------------------------------
    !// Adjust moments if tracer has been reduced.
    !// This is standard treatment for UCI CTM, to reduce the moments
    !// when tracer has been reduced. It should hinder the possibility
    !// that moments contain (and hence will try to move) more mass than
    !// is actually in the gridbox.
    !//
    !// Amund Sovde, September 2014
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom, r8
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK):: BTT, BTTBCK
    !// Input/output
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ

    !// Locals
    integer :: II,JJ,N,L
    real(r8) :: frac
    !// --------------------------------------------------------------------

    do JJ = 1, JDBLK
       do II = 1, IDBLK
          do N = 1, NPAR
             do L = 1, LPAR
                if ( BTT(L,N,II,JJ) .lt. BTTBCK(L,N,II,JJ)) then
                   frac = BTT(L,N,II,JJ) / BTTBCK(L,N,II,JJ)
                   BZT(L,N,II,JJ) = BZT(L,N,II,JJ) * real(frac, rMom)
                   BZZ(L,N,II,JJ) = BZZ(L,N,II,JJ) * real(frac, rMom)
                   BXZ(L,N,II,JJ) = BXZ(L,N,II,JJ) * real(frac, rMom)
                   BYZ(L,N,II,JJ) = BYZ(L,N,II,JJ) * real(frac, rMom)
                   BXT(L,N,II,JJ) = BXT(L,N,II,JJ) * real(frac, rMom)
                   BXX(L,N,II,JJ) = BXX(L,N,II,JJ) * real(frac, rMom)
                   BYT(L,N,II,JJ) = BYT(L,N,II,JJ) * real(frac, rMom)
                   BYY(L,N,II,JJ) = BYY(L,N,II,JJ) * real(frac, rMom)
                   BXY(L,N,II,JJ) = BXY(L,N,II,JJ) * real(frac, rMom)
                end if
             end do
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine adjust_moments
  !// ----------------------------------------------------------------------
!// ----------------------------------------------------------------------
  real(r8) function moninobukhov_length(SFCD,SFCT,USTAR,SFCS)
    !// --------------------------------------------------------------------
    !// Description: 
    !//  Calculate the Monin-Obukhov length which is used pbl-mixing.f90, 
    !//  seasalt.f90, fallingaerosols.f90, drydeposition.f90.
    !//  
    !//  It is directly taken from pbl-mixing.f90 (e.g. Garratt, 1992)
    !//
    !// History: 
    !//  Stefanie Falk, Mai 2019
    !// --------------------------------------------------------------------
    use cmn_parameters, only: G0, cp_air, VONKARMAN
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: &
         SFCD,SFCT,USTAR,SFCS
    !// Local variables
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'moninobukhov_length'
    real(r8), parameter :: cp_k_g = cp_air/(VONKARMAN*G0)
    !// Formerly cp_k_g = 256._r8 had been used
    !// --------------------------------------------------------------------
    moninobukhov_length = (-1._r8)*cp_k_g*SFCD*SFCT*(USTAR**3)/SFCS 
    !// Make sure we get physical values
    if (moninobukhov_length .eq. 0._r8) then
       write(6,'(a,2i5,4es16.6)') f90file//':'//subr//': MO_LEN is 0! -> 0.001: '// &
            'SFCD,SFCT,USTAR,SFCS',SFCD,SFCT,USTAR,SFCS
       moninobukhov_length = 1.e-3_r8
    end if
    if (moninobukhov_length .ne. moninobukhov_length) then
       write(6,'(a,2i5,4es16.6)') f90file//':'//subr//': MO_LEN is 0! -> 1000: '// &
            'SFCD,SFCT,USTAR,SFCS',SFCD,SFCT,USTAR,SFCS
       moninobukhov_length = 1.e3_r8
    end if
    return
  end function moninobukhov_length
  !// ----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module utilities
!//=========================================================================

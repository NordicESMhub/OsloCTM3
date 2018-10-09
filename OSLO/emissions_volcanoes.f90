!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Volcanoe emissions.
!//=========================================================================
module emissions_volcanoes
  !// ----------------------------------------------------------------------
  !// MODULE: emissions_volcanoes
  !// DESCRIPTION: Module setting up volcanoe emissions of SO2.
  !//
  !// HTAP data are given for year 1979-2010.
  !// You specify the path and year in the STV section of the emission
  !// file (Ltracer_emis). If year is 9999, the meteorological year
  !// will be used.
  !//
  !// Contains:
  !//   subroutine init_volcPATH
  !//   subroutine read_volcEMIS
  !//   subroutine add_volcEMIS
  !//
  !// Amund Sovde, March 2014
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  integer, parameter :: events_max = 426800
  integer :: events_tot

  !// Indices for volcanoe emissions
  integer,dimension(events_max) :: volc_ii, volc_jj, volc_mp
  !// Emissions of SO2
  real(r8),dimension(events_max)  :: volc_emis_so2
  !// Volcanoe elevation [m.a.s.l.]
  real(r8),dimension(events_max)  :: volc_elev
  !// Cloud column height [m.a.s.l.]
  real(r8),dimension(events_max)  :: volc_cch

  !// Path to the files
  character(len=100) :: volc_path
  !// Year for emissions
  integer :: volc_year

  !// Flag for using these emissions
  logical :: LHTAP_VOLC=.false.
  logical :: LACOM_VOLC=.false.

  !// Keep track of which events belonging to which days
  real(r8), dimension(events_max) :: volc_mm, volc_dd !// May not be necessary
  integer, dimension(31,12) :: volc_event_start, volc_event_end

  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public init_volcPATH, read_volcEMIS, add_volcEMIS
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine init_volcPATH(pathname,year, INFMT)
    !// --------------------------------------------------------------------
    !// Set volcanoe emission path and year of emissions.
    !//
    !// Amund Sovde, February 2014
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: pathname
    integer, intent(in) :: YEAR, INFMT
    !// --------------------------------------------------------------------
    volc_path = trim(pathname)
    if (INFMT .eq. 553) then
       LHTAP_VOLC = .true.
    else if (INFMT .eq. 554) then
       LACOM_VOLC = .true.
    else
       write(6,'(a,i5)') &
           '*** emissions_volcanoes.f90: No such INFMT for volcanoes:',INFMT
    end if
    if (YEAR.eq.9999) then
       !// Use meteorological year
       volc_year = MYEAR
    else
       !// Force year of emissions
       volc_year = YEAR
    end if
    write(6,'(a)') '* Volcanic path/file: '//trim(volc_path)
    !// --------------------------------------------------------------------
  end subroutine init_volcPATH
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_volcEMIS()
    !// --------------------------------------------------------------------
    !// Call read-in of emissions.
    !//
    !// Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LSULPHUR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Is SO2 included at all?
    if (.not. LSULPHUR) then
       LHTAP_VOLC = .false.
       LACOM_VOLC = .false.
       return
    end if

    if (LHTAP_VOLC) then
       call read_volcEMIS_HTAP()
    else if (LACOM_VOLC) then
       call read_volcEMIS_ACOM()
    end if
    !// --------------------------------------------------------------------
  end subroutine read_volcEMIS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_volcEMIS_HTAP()
    !// --------------------------------------------------------------------
    !// Read volcanoe emissions for current year.
    !// Routine needs to be called at model start and then once a year.
    !//
    !// Amund Sovde, February 2014
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LSULPHUR
    use cmn_ctm, only: XDEDG, YDEDG, all_mp_indices
    use ncutils, only: get_netcdf_var_1d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    character(len=120) :: filename
    character(len=4)   :: cyear
    integer            :: nEvents, iii,jjj, lastD,lastM, inY,inM,inD
    integer :: i,j,n, tmpDate
    real(r8)  :: temp_lat, temp_lon
    !real(r8)  :: ctmlonb(IPAR+1)

    real(r8), allocatable, dimension(:) :: &
         inLon, inLat, inEmis, inElev, inCch, inDate
    !// --------------------------------------------------------------------

    !// Initialise
    volc_ii(:) = 0
    volc_jj(:) = 0
    volc_emis_so2(:) = 0._r8
    volc_event_start(:,:) = 0
    volc_event_end(:,:)   = 0

    !// Is SO2 included at all?
    if (.not. LSULPHUR) then
       LHTAP_VOLC = .false.
       return
    end if

    !// Skip if not included
    if (.not. LHTAP_VOLC) return

    !// Fetch data from file
    write(cyear(1:4),'(i4.4)') volc_year
    filename = trim(volc_path)//'volc_so2_'//cyear//'.nc'
    write(6,'(a)') '* Reading volcanic SO2 emissions: '//trim(filename)

    call get_netcdf_var_1d( filename, 'lon',  inLon )
    nEvents  = SIZE( inLon  )

    if (nEvents .gt. events_max) then
       print*,'emissions_volcanoes.f90: nEvents > events_max'
       print*,'nEvents:   ',nEvents
       print*,'events_max:',events_max
       stop
    end if
    call get_netcdf_var_1d( filename, 'lat', inLat )
    call get_netcdf_var_1d( filename, 'date', inDate )
    call get_netcdf_var_1d( filename, 'so2', inEmis )
    call get_netcdf_var_1d( filename, 'elevation', inElev )
    call get_netcdf_var_1d( filename, 'cloud_column_height', inCch )


    !// To keep track of days/month
    lastM = -1
    lastD = -1

    do n = 1, nEvents

       !// get year/month/day info
       tmpDate = nint(inDate(n))
       inY = tmpDate/10000
       inM = (tmpDate - inY*10000)/100
       inD = (tmpDate - inY*10000 - inM*100)

       !// Keep track of start/end of events for each day
       if (inD .ne. lastD) then
          !// we have a new day
          volc_event_start(inD,inM) = n
          volc_event_end(inD,inM)   = n !// This works as an initialisation
          lastD = inD
       else
          !// Working on the same day as before; keep increasing the end event
          volc_event_end(inD,inM)   = n
       end if


       !// Get I/J index from lat/lon boundaries
       !// Initialize
       iii = -1
       jjj = -1
       temp_lat = inLat(n)
       temp_lon = inLon(n)

       !// Find grid boxes J-index
       do j = 1, JPAR
          if ( temp_lat.ge.YDEDG(j) .and. (temp_lat.lt.YDEDG(j+1)) ) then
             !// We have y-index
             jjj = j
             exit
          end if
       end do
       if (temp_lat.ge.YDEDG(jpar+1)) jjj = jpar !// special for NP
       if (temp_lat.lt.YDEDG(1)) jjj = 1         !// special for SP

       !// XDEDGE is defined as ~0 to ~360
       !ctmlonb(:) = XDEDG(:)
       !ctmlonb(IPAR/2+2:IPAR+1) = XDEDG(IPAR/2+2:IPAR+1) + 360._r8
 
       !// Find grid boxes I-index
       !// Longitude on degrees in the range <ctmlonb(1),ctmlonb(ipar+1)]
       if (temp_lon .lt. XDEDG(1)) temp_lon = temp_lon + 360._r8
       do I = 1, IPAR
          if ( (temp_lon .gt. XDEDG(i)) .and. &
               (temp_lon .le. XDEDG(i+1)) ) then
             !// We have x-index
             iii = i
             exit
          end if
       end do

       if (iii.lt.0 .or. jjj.lt.0) then
          write(6,'(a,2i5)') &
               '* emissions_volcanoes.f90: ERROR: III/JJJ is wrong:',iii,jjj
          stop
       end if

       !// Assign i/j through ii/jj/mp
       volc_ii(n) = all_mp_indices(1,iii,jjj)
       volc_jj(n) = all_mp_indices(2,iii,jjj)
       volc_mp(n) = all_mp_indices(3,iii,jjj)

       !// Emissions (convert from kt(SO2)/d to kg(SO2)/s
       volc_emis_so2(n) = inEmis(n) * 1.e6_r8 / 86400._r8
       !// Volcanoe elevation [m]
       volc_elev(n)     = inElev(n)
       !// Volcanoe plume heigh [m]
       volc_cch(n)     = inCch(n)

       if (volc_elev(n) .gt. volc_cch(n)) then
          !// Correcting wrong values (may do this the other way around?)
          volc_elev(n) = volc_cch(n)
       end if

    end do !// do n = 1, nEvents

    !// Deallocate all local variables
    if ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
    if ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
    if ( ALLOCATED(inDate) ) DEALLOCATE(inDate)
    if ( ALLOCATED(inEmis) ) DEALLOCATE(inEmis)
    if ( ALLOCATED(inElev) ) DEALLOCATE(inElev)
    if ( ALLOCATED(inCch) ) DEALLOCATE(inCch)

    !// --------------------------------------------------------------------
  end subroutine read_volcEMIS_HTAP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_volcEMIS_ACOM()
    !// --------------------------------------------------------------------
    !// Read volcanoe emissions.
    !// Routine needs to be called at model start and then once a year.
    !//
    !// Amund Sovde, March 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LSULPHUR
    use cmn_ctm, only: XDEDG, YDEDG, all_mp_indices
    use ncutils, only: get_netcdf_var_1d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: nEvents, iii,jjj
    integer :: I,J,N
    real(r8)  :: temp_lat, temp_lon
    !real(r8)  :: ctmlonb(IPAR+1)
    real(r8)  :: inEmis, inElev, inCch

    logical :: fnr_ok
    integer :: ifnr, io_err
    real(r8), parameter :: dtyear = 365._r8*24._r8*3600._r8
    !// --------------------------------------------------------------------

    !// Initialise
    volc_ii(:) = 0
    volc_jj(:) = 0
    volc_emis_so2(:) = 0._r8
    volc_event_start(:,:) = 0
    volc_event_end(:,:)   = 0

    !// Is SO2 included at all?
    if (.not. LSULPHUR) then
       LACOM_VOLC = .false.
       return
    end if

    !// Skip if not included
    if (.not. LACOM_VOLC) return

    !// Fetch data from file
    !// Will use file name as is; it does not change over years since
    !// emissions are continuous.
    write(6,'(a)') '* Reading volcanic SO2 emissions: '//trim(volc_path)

    !// File number for reading emission files
    fnr_ok = .true.
    ifnr = 7
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do

    open(ifnr,file=trim(volc_path),form='formatted',status='old')
    !// read header
    read(ifnr,*)
    read(ifnr,*)
    read(ifnr,*)

    !// Number of (continuous) events
    nEvents = 0
    volc_event_start(:,:) = 1
    volc_event_end(:,:)   = 1


    !// Read data (as events)
    do N=1,99999
       read(ifnr,'(2(1x,i4),3(1x,e10.4))',iostat=io_err) &
            I, J, inEmis, inElev, inCch
       if (io_err.ne.0) exit

       !// File starts at -180 E, 90 S, -179 E, 90 S etc.
       !// Get lat/lon
       temp_lon = real(i, r8) + 0.5_r8 ! e.g. -179.5W
       if (temp_lon .lt. 0._r8) temp_lon = temp_lon + 360._r8
       temp_lat = real(j, r8) + 0.5_r8 ! e.g. -89.5S

       !// Get I/J index from lat/lon boundaries
       !// Initialize
       iii = -1
       jjj = -1

       !// Find grid boxes J-index
       do j = 1, JPAR
          if ( temp_lat.ge.YDEDG(j) .and. (temp_lat.lt.YDEDG(j+1)) ) then
             !// We have y-index
             jjj = j
             exit
          end if
       end do
       if (temp_lat.ge.YDEDG(jpar+1)) jjj = jpar !// special for NP
       if (temp_lat.lt.YDEDG(1)) jjj = 1         !// special for SP

       !// convert XDEDGE to 0-360
       !ctmlonb(:) = XDEDG(:)
       !ctmlonb(IPAR/2+2:IPAR+1) = XDEDG(IPAR/2+2:IPAR+1) + 360._r8
 
       !// Find grid boxes I-index
       !// Longitude on degrees in the range <ctmlonb(1),ctmlonb(ipar+1)]
       if (temp_lon .lt. XDEDG(1)) temp_lon = temp_lon + 360._r8
       do I = 1, IPAR
          if ( (temp_lon .gt. XDEDG(i)) .and. &
               (temp_lon .le. XDEDG(i+1)) ) then
             !// We have x-index
             iii = i
             exit
          end if
       end do

       if (iii.lt.0 .or. jjj.lt.0) then
          write(6,'(a,2i5)') &
               '* emissions_volcanoes.f90: ERROR: III/JJJ is wrong:',iii,jjj
          stop
       end if

       !// Assign i/j through ii/jj/mp
       volc_ii(n) = all_mp_indices(1,iii,jjj)
       volc_jj(n) = all_mp_indices(2,iii,jjj)
       volc_mp(n) = all_mp_indices(3,iii,jjj)

       !// Emissions (convert from kg/yr to kg/s
       volc_emis_so2(n) = inEmis / dtyear
       !// According to the emission file, AEROCOM emissions are to be
       !// distributed linearly between elevation and plume altitude.
       !// The routine add_volcEMIS, however, distributes in the
       !// upper 1/3 of the plume. We stick to that anyway.
       !// Volcanoe elevation [m]
       volc_elev(n)     = inElev
       !// Volcanoe plume heigh [m]
       volc_cch(n)     = inCch

       if (volc_elev(n) .gt. volc_cch(n)) then
          !// Correcting wrong values (may do this the other way around?)
          volc_elev(n) = volc_cch(n)
       end if

       !// Save number of events
       nEvents = N

    end do !// do n = 1, 99999

    !// Since emissions are continuous, all days have emissions from
    !// event 1 to nEvents.
    volc_event_end(:,:)   = nEvents

    print*,'number of events ACOM volcanic emissions',nevents

    !// --------------------------------------------------------------------
  end subroutine read_volcEMIS_ACOM
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine add_volcEMIS(BX,DT,DAY,MONTH,MP)
    !// --------------------------------------------------------------------
    !// Routine to add volcanic SO2 emissions.
    !// Routine is called from either emis4chem_oslo or SOURCE.
    !// IMPORTANT:
    !//   Note that the arguments are different things in the
    !//   two routines!:
    !//     emis4chem_oslo: BX is BEMIS [kg/s] and needs DT=1.d0
    !//     SOURCE:         BX is BTT [kg] and needs DT to be the time step.
    !//
    !// Amund Sovde, March 2015, February 2014
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, MPBLK, LSULPHUR
    use cmn_ctm, only:  MPBLKIB, MPBLKJB, all_mp_indices
    use cmn_met, only: ZOFLE
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: DAY, MONTH, MP
    real(r8), intent(in)  :: DT
    !// In/Out: Emissions (emis4chem_oslo) or tracer array (SOURCE)
    real(r8), intent(inout)  :: BX(LPAR,NPAR,IDBLK,JDBLK)

    !// Local
    integer :: N_SO2, L,I,J, II,JJ, LV1, LV2, N
    real(r8) :: frac, DZ, ZH, topoH, plumeH, topZ3bot
    !// --------------------------------------------------------------------

    !// Skip volcanic emissions?
    if (.not. LSULPHUR) return
    if (.not. (LHTAP_VOLC .or. LACOM_VOLC)) return

    N_SO2 = trsp_idx(72) !// Transport number for SO2

    !// Loop through all events of the day
    do n = volc_event_start(DAY,MONTH), volc_event_end(DAY,MONTH)

       !// Check for this MP-block
       if (volc_mp(n) .eq. MP) then

          !// MP-block indices
          II = volc_ii(n)
          JJ = volc_jj(n)
          !// Global indices
          I = II - 1 + MPBLKIB(MP)
          J = JJ - 1 + MPBLKJB(MP)


          if (volc_elev(n) .ge. volc_cch(n)) then
             !// Volcanoe altitude is the same as the cloud height it
             !// generates. In other words it is non-eruptive, but degassing.
             !// In this case, place emissions into the model level of the
             !// crater elevation.

             !// Check volcanoe altitude vs model levels
             LV1 = 1
             do L = 1, LPAR
                if (volc_elev(n) .ge. ZOFLE(L,I,J) .and. &
                     volc_elev(n) .lt. ZOFLE(L+1,I,J)) then
                   LV1 = L !// Model level where volcanoe is located
                   exit
                end if
             end do
             BX(LV1,N_SO2,II,JJ) = BX(LV1,N_SO2,II,JJ) &
                  + volc_emis_so2(n) * DT
          else
             !// Eruptive volcanoe. Assume emissions are placed in the
             !// levels spanning the top 1/3 of the volcano plume.
             !// IMPORTANT: We only look for the model levels, and do not care
             !//            about interpolation to get exact thickness of the
             !//            plume. Emissions are to be equally distributed in
             !//            the model levels containing the top 1/3 of the
             !//            plume.
             plumeH = volc_cch(n)
             !// Bottom of top 1/3
             topZ3bot = plumeH - (plumeH - volc_elev(n)) / 3._r8

             !// Locate model levels where we find the top 1/3 of plume.
             if (topZ3bot .lt. ZOFLE(1,I,J)) then
                !// In case bottom of top 1/3 is below model surface
                LV1 = 1
             else
                lv1 = -1
                do L = 1, LPAR
                   if (topZ3bot .ge. ZOFLE(L,I,J) .and. &
                        topZ3bot .lt. ZOFLE(L+1,I,J)) then
                      !// This is the level where top 1/3 of volcanoe
                      !// plume starts
                      lv1 = L
                      exit
                   end if
                end do
                if (lv1 .lt. 0) then
                   write(6,'(a,f9.2)') &
                        'emissions_volcanoes.f90: Volcanoe LV1 is not set',topZ3bot
                   stop
                end if
                !// Limit LV1 just in case
                if (LV1 .gt. LPAR) LV1 = LPAR
             end if


             lv2 = lv1 !// In case plumeH < ZOFLE(1,I,J)
             do L = lv1, LPAR
                if (plumeH .ge. ZOFLE(L,I,J) .and. &
                     plumeH .lt. ZOFLE(L+1,I,J)) then
                   !// This is the level of the plume top
                   lv2 = L
                   exit
                end if
             end do
             if (lv2 .lt. 0) then
                print*,'Volcanoe LV2 is not set',plumeH
                stop
             end if
             !// Limit LV2 just in case
             if (LV2 .gt. LPAR) LV2 = LPAR


             !// Thickness of model layers spanning top 1/3 of volcanoe plume
             dZ = ZOFLE(LV2+1,I,J) - ZOFLE(LV1,I,J)
             do L = LV1, LV2
                frac = (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) / DZ
                BX(L,N_SO2,II,JJ) = BX(L,N_SO2,II,JJ) &
                     + volc_emis_so2(n) * DT * frac
                !write(6,'(a4,1x,2i3,2f8.4,4f8.1)') 'VOLC', &
                !     n,L,volc_emis_so2(n) * DT, frac,&
                !     ZOFLE(L,I,J),ZOFLE(L+1,I,J),volc_elev(n),volc_cch(n)
             end do

          end if !// if (volc_elev(n).le. volc_cch(n)) then

       end if !// if (volc_mp(n) .eq. MP) then

    end do !// do n = volc_event_start(DAY,MONTH), volc_event_end(DAY,MONTH)

    !// --------------------------------------------------------------------
  end subroutine add_volcEMIS
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module emissions_volcanoes
!//=========================================================================

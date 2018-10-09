!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// CARIBIC II output from Oslo CTM3
!//=========================================================================
module caribic2
  !// ----------------------------------------------------------------------
  !// MODULE: caribic2
  !// DESCRIPTION: Routines for putting out 10-sec CARIBIC data.
  !//
  !// Ole Amund Sovde, March 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, &
       LOSLOCSTRAT
  use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
       AIR, STT, AREAXY, NROPSM, NRMETD
  use cmn_met, only: ZOFLE, P, T, Q, BLH_CUR, BLH_NEXT, PVU
  use cmn_oslo, only: LMTROP, XSTT, trsp_idx, Xtrsp_idx
  use physics_oslo, only: NTHE, pvthe, theqlat
  use strat_h2o, only: sumH2
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// STT components
  !// O3, CO, C2H4, C2H6, HCHO, CH3CHO, NO, NO2, CH4, ACETONE,SO2, N2O, NOy
  !//  1,  6,    7,    8,   13,     14, 43,  44,  46,      50, 72, 107, 147
  integer,parameter :: ntracer=14
  integer,dimension(ntracer),parameter :: out_comps = &
       (/ 1,4,6,7,8,13,14,43,44,46,50,107,124,147 /)
  real(r8), dimension(ntracer)           :: mole_mass
  !character(len=10), dimension(ntracer):: out_names

  !// Version
  integer, parameter :: version=2

  integer                      :: n_event, cur_event
  integer,dimension(6)         :: idate_event
  integer                      :: itau_event
  integer :: file_id

  type observation
     integer,dimension(6)       :: idate     !// date of observation
     integer                    :: itau      !// seconds since reference time
     character(len=24)          :: ident     !// observation identifier
     real(r4)                   :: lat       !// latitude of observation
     real(r4)                   :: lon       !// longitude of observation
     real(r4)                   :: height    !// height of observation
     integer                    :: eventnumber !// event number
     integer                    :: ii        !// longitudinal grid box index
     integer                    :: jj        !// latitudinal grid box index
     integer                    :: nb_ii     !// neighbor lon. grid box index
     integer                    :: nb_jj     !// neighbor lat. grid box index
     real(r4)                   :: nb_xfrac
     real(r4)                   :: nb_yfrac
     integer                    :: ctm_time  !// CTM time [sec] since midnight
     real(r8),dimension(4)      :: areaxy    !// Surface area
     integer,dimension(2)       :: ctm_itau  !// CTM itau
     !// For each field, store data from the two closest NOPS
     !// 1: prior to measurement, 2: after
     real(r8),dimension(4,2)     :: psfc          !// Surface pressure
     real(r8),dimension(4,2)     :: blh           !// Boundary layer height
     real(r8),dimension(LPAR,ntracer,4,2) :: mass !// Tracer mass
     real(r8),dimension(LPAR,4,2) :: h2o          !// H2O (not part of STT)
     real(r8),dimension(LPAR,4,2) :: temperature  !// Temperature
     real(r8),dimension(LPAR,4,2) :: airmass      !// Air mass
     !real(r8),dimension(LPAR,4,2) :: cldfr        !// Cloud fraction (v3)
     real(r8),dimension(LPAR+1,4,2):: zoflev !// Height of box bottoms (1 is topography)
     real(r8), dimension(NTHE,4,2) :: eqlat        !// Equivalent latitude
     real(r8),dimension(LPAR,4,2)  :: pvu          !// Potential vorticity units
     integer,dimension(4,2)      :: tph_lev !// Uppermost troposphere grid box index
  end type observation


  !// Max observations per day
  integer,parameter :: maxOPD = 9000
  type(observation), dimension(maxOPD) :: observations

  !// Date information for the day in progress
  integer :: curyear, curdate, curmonth

  logical :: nodata, endoffile
  logical, parameter :: verbose=.false.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'caribic2.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public caribic2_master
  !// ----------------------------------------------------------------------

contains

  !// --------------------------------------------------------------------


  !// --------------------------------------------------------------------
  subroutine get_new_events(JYEAR,JMON,JDATE,NDAY,NDAYI, NMET,NOPS,LNEWM)
    !// --------------------------------------------------------------------
    !// Find observations for this day.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: XDEDG, YDEDG, YDGRD, XDGRD, LFIXMET
    use cmn_chem, only: TMASS
    use cmn_met, only: MYEAR
    use cmn_oslo, only: XTMASS
    use dateconv, only: idate2itau
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: JYEAR,JMON,JDATE,NDAY,NDAYI, NMET,NOPS
    logical, intent(in) :: LNEWM

    !// Locals
    real(r8)  :: ctmlonb(IPAR+1),ctmlon(IPAR),deltalon,deltalat
    logical :: file_io, linread, okdata
    character(len=80) :: filename
    integer           :: iii, jjj, nii, njj, nbb, nmp, N, i,j
    integer           :: hh,mm,ss
    integer :: file_err

    integer :: temp_idate(6),temp_itau, scount
    integer :: now_date(6),now_itau,now_itau_p
    real(r8)   :: nbxfrac,nbyfrac
    integer :: temp_pres,temp_ident, NO, NM, hour, minute, secs
    real(r8)  :: temp_lon,temp_lat, time, dtmet, dtops
    character(len=4) :: cyear
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_new_events'
    !// --------------------------------------------------------------------

    temp_idate(6)=0
    !// Initialize number of events
    n_event = 0
    !// Initialize current event
    cur_event = 1

    !// Open file, do this every month for simplicity
    if (LNEWM) then
       !// Check you have the tracers available
       do N = 1, ntracer
          if ( .not. (trsp_idx(out_comps(n)).gt.0 .or. &
                     Xtrsp_idx(out_comps(n)).gt.0) ) then
             write(6,'(a,i5,a)') '* '//f90file//':'//subr// &
                  ': Tracer number ',out_comps(n), &
                  ' is not included in the simulation!'
             stop
          endif
       enddo

       !// Get tracer molecular masses
       do N = 1, ntracer
          if ( trsp_idx(out_comps(n)).gt.0) then
             mole_mass(N) = TMASS(trsp_idx(out_comps(n)))
          else
             mole_mass(N) = XTMASS(Xtrsp_idx(out_comps(n)))
          endif
       enddo

       !// Meteorological year
       write(cyear,'(i4.4)') MYEAR

       !// The flight list
       filename='Indata_CTM3/observations/caribic_'//cyear//'_10sec.dat'

       write(6,'(a)') '* '//f90file//':'//subr//': Checking file'


       !// Find a file number
       file_id=40
       file_io=.true.
       do while (file_io)
          file_id=file_id+1
          inquire(file_id,opened=file_io)
       enddo

       nodata=.false.
       endoffile=.false.

       open(file_id,file=filename,form='formatted',status='old',iostat=file_err)
       if (file_err .ne. 0) then
          write(6,'(a,a)') '*** Could not open file: ',trim(filename)
          write(6,'(a)')   '    Will not simulate caribic flights!'
          return
       else
          write(6,'(a)') '* Reading file: '//trim(filename)
       endif

       !// Input file is open, read header
       read(file_id,*)

    endif

    !// Check for end of file
    if (endoffile) return

    write(6,'(a)') '* '//f90file//':'//subr//': Initializing caribic measurements'

    !// Time stamp now
    now_date = (/JYEAR,JMON,JDATE,0,0,0/)
    call idate2itau(now_date,now_itau)


    !// Read off observations prior to model time (spin-off)
    temp_itau = 0
    scount = 0
    do while (temp_itau .lt. now_itau)
       read(file_id,'(i4,i2,i2,1x,i2,i2,i2,f10.4,f9.4,i7,i4)',iostat=file_err) &
            temp_idate(1),temp_idate(2),temp_idate(3), &
            temp_idate(4),temp_idate(5),temp_idate(6), &
            temp_lon,temp_lat,temp_pres,temp_ident
       !// Count spinoffs
       scount = scount + 1

       !// If end of file...
       if (file_err .ne. 0) then
          print*,'Read error1',file_err,temp_itau, now_itau
          if (.not.endoffile) print*,'End of caribic input file'
          endoffile=.true.
          nodata = .true.
          close(file_id)
          exit
       endif


       !// Find tau
       call idate2itau(temp_idate,temp_itau)

       if (verbose) write(6,'(a,2(i4,2i2.2,1x,3i2.2,1x),2i16)') 'SKIP ', &
            now_date,temp_idate, now_itau,temp_itau

    enddo

    !// If spinup read successfully, step back to previous line
    if (scount.gt.0 .and. .not.endoffile) backspace(file_id)

    !// Skip calculations when no data is available
    if (nodata) return

    !// Time steps
    dtmet = 24._r8/real(NRMETD, r8)
    dtops = dtmet/real(NROPSM, r8)

    !// longitude boundary
    ctmlonb(1) = XDEDG(1)
    do i = 2, ipar+1
       if (XDEDG(i) .lt. 0._r8) then
          ctmlonb(i) = XDEDG(i) + 360._r8
       else
          ctmlonb(i) = XDEDG(i)
       end if
    end do
    !// longitude center
    do i = 1, ipar
       if (XDGRD(i) .lt. 0._r8) then
          ctmlon(i) = XDGRD(i) + 360._r8
       else
          ctmlon(i) = XDGRD(i)
       end if
    end do

    do NM = 1,NRMETD
       do NO = 1,NROPSM

          time   = real(NM-1, r8)*dtmet + real(NO-1, r8)*dtops
          hour   = floor(time)
          minute = floor((time - hour) * 60._r8)
          secs   = nint((time - hour)*60._r8 - minute) * 60._r8

          !// Time stamp for this NOPS and the next
          now_date = (/ JYEAR,JMON,JDATE,hour,minute,secs /)
          call idate2itau(now_date,now_itau)
          !// Next NOPS in seconds
          now_itau_p = now_itau + nint(dtops*3600._r8)

          okdata=.true.
          do while (okdata)
             read(file_id,'(i4,2i2,1x,3i2,f10.4,f9.4,i7,i4)',iostat=file_err) &
                  temp_idate(1),temp_idate(2),temp_idate(3),&
                  temp_idate(4),temp_idate(5),temp_idate(6),&
                  temp_lon,temp_lat,temp_pres,temp_ident

             if (file_err .ne. 0) then
                print*,'Read error2',file_err
                print*,temp_itau, now_itau
                if (.not.endoffile) print*,'End of caribic input file'
                endoffile = .true.
                close(file_id)
                exit
             endif


             !// Find tau for observation
             call idate2itau(temp_idate,temp_itau)

             !// Save observation data for this NOPS
             if (temp_itau .ge. now_itau .and. temp_itau .lt. now_itau_p) then


                !// Increase number of observations
                n_event = n_event + 1
                if (n_event .gt. maxOPD) then
                   write(6,'(a)') '* '//f90file//':'//subr//': n_event > maxOPD'
                   stop
                endif

                write(6,'(a,2(i4,2i2.2,1x,3i2.2,1x),2i16,i6)') ' OK ', &
                     now_date,temp_idate, now_itau,temp_itau,n_event




                !// Find grid info


                !// Initialize
                iii = -1
                jjj = -1

                !// Find grid boxes J-index
                do j = 1, JPAR
                   if ( temp_lat.ge.YDEDG(j) .and. &
                        (temp_lat.lt.YDEDG(j+1)) ) then
                      !// We have y-index
                      jjj=j
                      exit
                   endif
                enddo
                if (temp_lat.ge.YDEDG(jpar+1)) jjj=jpar !// special for NP
                if (temp_lat.lt.YDEDG(1)) jjj=1         !// special for SP


                !// Find grid boxes I-index
                !// Longitude on degrees in the range
                !// <ctmlonb(1),ctmlonb(ipar+1)]
                if (temp_lon .lt. ctmlonb(1)) temp_lon = temp_lon + 360._r8
                do I = 1, IPAR
                   if ( (temp_lon .gt. ctmlonb(i)) .and. &
                        (temp_lon .le. ctmlonb(i+1)) ) then
                      !// We have x-index
                      iii=i
                      exit
                   endif
                enddo

                if (iii.lt.0 .or. jjj.lt.0) then
                   write(6,'(a,2i5)') '* '//f90file//':'//subr// &
                        ': ERROR: III/JJJ is wrong:',iii,jjj
                   stop
                endif


                !// Find neighbour boxes: closest i-box
                if (temp_lon .gt. ctmlon(iii)) then
                   !// east of grid box center
                   nii = iii + 1
                else
                   !// west of grid box center
                   nii = iii - 1
                endif
                !// check the boundaries
                if (nii .eq. (ipar+1)) nii = 1
                if (nii .eq. 0) nii = ipar


                !// Find neighbour boxes: closest j-box
                if (temp_lat .gt. YDGRD(jjj)) then
                   njj = jjj + 1
                else
                   njj = jjj - 1
                endif
                !// Check the boundaries
                if (njj .eq. 0) njj = 1         !// special for SP
                if (njj .eq. jpar+1) njj = jpar !// special for NP


                !// Find fraction of neighbour boxes
                deltalon = ctmlon(2) - ctmlon(1)
                if (jjj .gt. 1) then
                   deltalat = YDGRD(jjj) - YDGRD(jjj-1)
                else
                   deltalat = YDGRD(2) - YDGRD(1) !// special for SP
                endif
                nbyfrac = abs(temp_lat - YDGRD(jjj))/deltalat
                nbxfrac = abs(temp_lon - ctmlon(iii))/deltalon

                !// Check fractions
                if (nbxfrac .gt. 1._r8) then
                   write(6,'(a)') '* '//f90file//':'//subr// &
                        ': Something wrong with lon:'
                   print*,nbxfrac,temp_lon,xdgrd(iii),ctmlon(iii),deltalon
                   print*,iii,nii
                   stop
                endif
                if (nbyfrac .gt. 1._r8) then
                   write(6,'(a)') '* '//f90file//':'//subr// &
                        ': Something wrong with lat:'
                   print*,nbyfrac,temp_lat,YDGRD(jjj)
                   print*,jjj,njj
                   stop
                endif

                !// Update measurement array

                !// Fill the observation info
                write(observations(n_event)%ident(1:24),'(i24)') temp_ident
                observations(n_event)%lat = temp_lat
                observations(n_event)%lon = temp_lon
                observations(n_event)%height = temp_pres
                observations(n_event)%itau = temp_itau    !// itau for observation
                observations(n_event)%idate = temp_idate
                observations(n_event)%eventnumber = n_event
                observations(n_event)%ctm_time    = 0 !// CTM seconds of day
                observations(n_event)%ctm_itau(:) = 0 !// CTM itau for both NOPS

                observations(n_event)%jj = jjj      !// CTM j-index
                observations(n_event)%ii = iii      !// CTM i-index
                observations(n_event)%nb_ii = nii  !// CTM closest i-neighbor
                observations(n_event)%nb_jj = njj  !// CTM closest j-neighbor
                observations(n_event)%nb_xfrac = nbxfrac !// Fractional distance x-direction
                observations(n_event)%nb_yfrac = nbyfrac !// Fractional distance y-direction
                observations(n_event)%areaxy(1)   = AREAXY(iii,jjj)
                observations(n_event)%areaxy(2)   = AREAXY(nii,jjj)
                observations(n_event)%areaxy(3)   = AREAXY(iii,njj)
                observations(n_event)%areaxy(4)   = AREAXY(nii,njj)

                observations(n_event)%mass(:,:,:,:) = -999._r8 ! tracer mass
                observations(n_event)%psfc(:,:) = 0._r8
                observations(n_event)%blh(:,:) = 0._r8
                observations(n_event)%temperature(:,:,:) = 0._r8
                observations(n_event)%airmass(:,:,:) = 0._r8
                !observations(n_event)%cldfr(:,:,:) = 0._r8
                observations(n_event)%zoflev(:,:,:) = 0._r8
                observations(n_event)%eqlat(:,:,:) = 0._r8
                observations(n_event)%pvu(:,:,:) = 0._r8
                observations(n_event)%tph_lev(:,:) = 0

             else
                if (temp_itau .lt. now_itau) then
                   write(6,'(a)') '* '//f90file//':'//subr// &
                        ': temp_itau < now_itau'
                   write(6,'(a)') '  Should have been read off at start'
                   stop
                endif
                if (verbose) write(6,'(a,2(i4,2i2.2,1x,3i2.2,1x),2i16)') ' NOT ', &
                     now_date,temp_idate, now_itau,temp_itau
                !// Data for the next model hour
                okdata = .false.
                !// Step back to previous line for next hour
                if (temp_itau .ge. now_itau_p) backspace(file_id)

             endif

          enddo !// do while (okdata)

          !// End of file (read error), exit loop
          if (endoffile) exit

       enddo !// do NO = 1,NROPSM

       !// End of file (read error), exit loop
       if (endoffile) exit

    enddo !// do NM=1,NRMETD

    if (file_err .ne. 0) then
       write(6,'(a,i5)') '* '//f90file//':'//subr// &
            ': Read error3',file_err
       write(6,'(a,2i20)') '  temp_itau, now_itau: ',temp_itau, now_itau
       write(6,'(a)') '  End of '//f90file//' input file?'
       return
    endif

    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'
    write(6,'(a,i5,a)')'CARIBIC: ',n_event,' events from file.'
    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'


    !// --------------------------------------------------------------------
  end subroutine get_new_events
  !// --------------------------------------------------------------------





  !// --------------------------------------------------------------------
  subroutine caribic_output(JYEAR,JMON,JDATE, NDAY, NMET, NOPS)
    !// --------------------------------------------------------------------
    !// Produces output for the given NOPS and the previous NOPS.
    !//
    !// Ole Amund Sovde, January 2011
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    use dateconv, only: idate2itau
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input variables
    integer,intent(in) :: JYEAR,JMON,JDATE, NDAY, NMET, NOPS

    !// Local variables
    integer :: nev, nbb, i, j, L, n, NN, tops
    integer :: max_itau, min_itau, temp_itau, &
         hour, minute, secs
    integer :: now_date(6),now_itau,now_itau_m,now_itau_p
    real(r8) :: dtmet, dtops, time
    real(r8) :: mass_CH4, mass_H2, H2OTMP(LPAR), blh_nops, frac_next

    logical :: LFIRST, LLAST

    !// For converting H2O from concentration to mass
    real(r8),parameter :: C2MH2O  = 18._r8/M_AIR, &
         VMR2MMRH2 = 2._r8/M_AIR, &
         MW_H2 = 2._r8, MW_CH4 = 16._r8, MW_H2O = 18._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'caribic_output'
    !// --------------------------------------------------------------------

    !// Time steps for NMET and NOPS
    dtmet = 24._r8/real(NRMETD, r8)
    dtops = dtmet/real(NROPSM, r8)
    frac_next = real(NOPS-1, r8) / real(NROPSM) !// For BLH temp.interp.


    time   = real(NMET-1, r8)*dtmet + real(NOPS-1, r8)*dtops
    hour   = floor(time)
    minute = floor((time - hour) * 60._r8)
    secs   = nint((time - hour)*60._r8 - minute) * 60._r8

    !// Time stamp for this NOPS and the next
    now_date = (/ JYEAR,JMON,JDATE,hour,minute,secs /)
    call idate2itau(now_date,now_itau)
    !// Next NOPS in seconds
    now_itau_p = now_itau + nint(dtops*3600._r8)
    !// Previous NOPS in seconds
    now_itau_m = now_itau - nint(dtops*3600._r8)

    !// Will save data for the two closest model hours (plus minus one NOPS).


    !// Check for first time step of day
    if (NMET .eq. 1 .and. NOPS.eq.1) then
       LFIRST = .true.
    else
       LFIRST = .false.
    endif

    !// Last time step is 00 the next day. Only when
    !// routine is called from output routine.
    if (NMET .eq. (NRMETD+1) .and. NOPS.eq.1) then
       LLAST = .true.
    else
       LLAST = .false.
    endif



    !// 00: Only store this time step
    !// 01-23: Store data for measurements done the last hour


    do nev = 1, n_event

      !// If observations this hour, save output
      !// If not, loop on
      temp_itau = observations(nev)%itau

      !// Store for data this NOPS (1) or for data measured last NOPS (2)
      do tops = 1, 2
         
        !// Save data for this NOPS; skip for NMET==9, which
        !// is only for catching the last time step of the day
        !// TOPS=1: Save for this NOPS. This is doe for all NOPS.
        !// TOPS=2: Save for NOPS+1. This is carried out at the next NOPS, and
        !//         is therefore also done for NOPS=1/NMET=9 to finish last day
        !//         before writing to file.
        if ( (tops.eq.1 .and. (.not. LLAST)) .or. &
             (tops.eq.2 .and. (.not. LFIRST)) ) then

          !// Possible routes:
          !// tops:1 LLAST:false: Put out data for this NOPS
          !// tops:2 LFIRST:false: Put out data data measured last NOPS
          !//        which should not be done for NOPS=1
          !// other: Nothing to be done

          if (tops .eq. 1) then
            !// Start/end time for observations to be saved this NOPS
            !// i.e. observations between this NOPS and the next.
            min_itau = now_itau
            max_itau = now_itau_p
          else if (tops .eq. 2) then
            !// Start/end time for observations to be saved at NOPS+1
            !// (relative to tops=1), i.e. the same observations that
            !// were stored by tops=2 at (previous NOPS).
            min_itau = now_itau_m
            max_itau = now_itau
          else
            write(6,'(a,i5)') '* '//f90file//':'//subr// &
                  ': tops wrong',tops
            stop
          endif

          !// Check the time
          if (temp_itau .ge. min_itau .and. &
              temp_itau .lt. max_itau) then
            !// Save CTM time stamps
            if (tops.eq.1) then
               observations(nev)%ctm_time = hour*3600 + minute*60 + secs
               observations(nev)%ctm_itau(1) = now_itau
               observations(nev)%ctm_itau(2) = now_itau_p
            endif

            if (verbose) print*,'CARIBIC',nev,tops,lfirst,llast,min_itau,&
                 temp_itau,max_itau

            !// Loop through neighboring boxes
            !// 1: iii/jjj, 2: nii/jjj, 3: iii/njj, 4: nii/njj
            do nbb = 1, 4
              if (nbb.eq.1) then
                 i = observations(nev)%ii
                 j = observations(nev)%jj
              elseif (nbb.eq.2) then
                 i = observations(nev)%nb_ii
                 j = observations(nev)%jj
              elseif (nbb.eq.3) then
                 i = observations(nev)%ii
                 j = observations(nev)%nb_jj
              elseif (nbb.eq.4) then
                 i = observations(nev)%nb_ii
                 j = observations(nev)%nb_jj
              else
                 write(6,'(a,i5)') '* '//f90file//':'//subr// &
                  ': NBB is wrong!',nbb
                 stop
              endif

              blh_nops =   (1._r8 - frac_next) * BLH_CUR(i,j) &
                                  + frac_next * BLH_NEXT(i,j)

              

              !// Calculate H2O mass
              H2OTMP(:) = Q(I,J,:) * AIR(I,J,:)
              if (LOSLOCSTRAT) then
                if (trsp_idx(114) .gt. 0) then
                  Do L = LMTROP(I,J)+1, LPAR
                    !// Stratospheric H2O
                    H2OTMP(L) = STT(I,J,L,trsp_idx(114))
                  End Do
                else
                  !// Need to calculate H2O from sum of H2
                  Do L = LMTROP(I,J)+1, LPAR
                    !// Get CH4
                    if (trsp_idx(46).gt.0) then
                      mass_CH4 = STT(I,J,L,trsp_idx(46))
                    else if (Xtrsp_idx(46).gt.0) then
                      mass_CH4 = XSTT(L,Xtrsp_idx(46),I,J)
                    else
                      write(6,'(a)') '* '//f90file//':'//subr// &
                            ': CH4 not available'
                      stop
                    endif
                    !// Get H2
                    if (trsp_idx(113).gt.0) then
                      mass_H2 = STT(I,J,L,trsp_idx(113))
                    else if (Xtrsp_idx(113).gt.0) then
                      mass_H2 = XSTT(L,Xtrsp_idx(113),I,J)
                    else
                      write(6,'(a)') '* '//f90file//':'//subr// &
                            ': H2 not available'
                      stop
                    endif
                    !// Use moles to calculate stratospheric H2O mass
                    ! vmr: molec/molecair, molec = mass /M * Na
                    H2OTMP(L) = ( sumH2 &
                       * VMR2MMRH2 * AIR(I,J,L) / MW_H2 & !// Sum H2 [kmol/gbox]
                       - 2._r8 * mass_CH4 / MW_CH4 &      !// -2*CH4 [kmol/gbox]
                       - mass_H2  / MW_H2  &              !// - H2 [kmol/gbox]
                       ) * MW_H2O !// From kmol to kg
                  End Do
                endif
              endif

              !// Save the variables (position tops)
              observations(nev)%psfc(nbb,tops)           = p(i,j)
              !// Boundary layer height
              observations(nev)%blh(nbb,tops)            = blh_nops
              !// Tracers
              do n = 1, ntracer
                 !// Transported
                 NN = trsp_idx(out_comps(n))
                 if (NN.gt.0) then
                    observations(nev)%mass(:,n,nbb,tops) = STT(i,j,:,NN)
                 else
                    !// Non-transported if not transprorted (checked in read-in)
                    NN = Xtrsp_idx(out_comps(n))
                    observations(nev)%mass(:,n,nbb,tops) = XSTT(:,NN,i,j)
                 endif
              enddo

              !// H2O
              observations(nev)%h2o(:,nbb,tops)          = H2OTMP(:)
              !// Temperature
              observations(nev)%temperature(:,nbb,tops)  = T(i,j,:)
              !// Air mass
              observations(nev)%airmass(:,nbb,tops)      = AIR(i,j,:)
              !// Cloud fraction
              !observations(nev)%cldfr(:,nbb,tops)        = CLDFR(i,j,:)
              !// Height of box bottom
              observations(nev)%zoflev(:,nbb,tops)       = ZOFLE(:,i,j)
              !// Equivalent latitudes on NTHE theta levels
              observations(nev)%eqlat(:,nbb,tops)        = theqlat(i,j,:)
              !// Potential voritcity units (PVU)
              observations(nev)%pvu(:,nbb,tops)          = pvu(:,i,j)
              !// Uppermost troposphere grid box index
              observations(nev)%tph_lev(nbb,tops)        = LMTROP(i,j)

            enddo !// do nbb = 1, 4

          endif !// if (temp_itau ...

        endif !// if ( (tops.eq.1 .and. ...
      enddo !// do tops = 1, 2
    enddo !// do nev = 1, n_event




    !// --------------------------------------------------------------------
  end subroutine caribic_output
  !// --------------------------------------------------------------------



  !// --------------------------------------------------------------------
  subroutine caribic_data_to_file(NDAY)
    !// --------------------------------------------------------------------
    !// Write collected flight path data to file.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: ETAA, ETAB
    use cmn_oslo, only: RESULTDIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY

    !// Locals
    integer :: error
    character(LEN=160) :: filename
    character(len=3) :: cday
    logical :: fnr_ok
    integer :: ifnr, n
    !// --------------------------------------------------------------------

    !// If no data, then no file
    if (n_event .eq. 0) return

    !// File name
    write(cday(1:3),'(i3.3)') NDAY
    filename = trim(RESULTDIR)//'caribic2_day'//cday//'.dta'


    !// Find file number to use
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    enddo
    !// Open file
    open(ifnr,file=filename,form='unformatted')
    !// Write info
    write(ifnr) curyear,curmonth,curdate,version !// Time & version
    write(ifnr) lpar,ntracer,n_event, nthe,ipar,jpar !// Resolution
    write(ifnr) etaa, etab             !// Sigma coordinates
    write(ifnr) out_comps              !// Components (STT)
    write(ifnr) mole_mass              !// Molecular masses of components
    !write(ifnr) out_names              !// Tracer names of components
    write(ifnr) pvthe                  !// Theta levels for eqlat
    do n = 1, n_event
       if (verbose) print*,'CARIBIC to file; writing event',N
       write(ifnr) observations(n)%idate
       write(ifnr) observations(n)%itau
       write(ifnr) observations(n)%ident
       write(ifnr) observations(n)%lat
       write(ifnr) observations(n)%lon
       write(ifnr) observations(n)%height
       write(ifnr) observations(n)%ii
       write(ifnr) observations(n)%jj
       write(ifnr) observations(n)%nb_ii
       write(ifnr) observations(n)%nb_jj
       write(ifnr) observations(n)%nb_xfrac
       write(ifnr) observations(n)%nb_yfrac

       write(ifnr) observations(n)%ctm_time
       write(ifnr) observations(n)%ctm_itau
       write(ifnr) real(observations(n)%psfc, r4)
       write(ifnr) real(observations(n)%blh, r4)
       write(ifnr) real(observations(n)%areaxy, r4)
       write(ifnr) real(observations(n)%mass, r4)
       write(ifnr) real(observations(n)%h2o, r4)
       write(ifnr) real(observations(n)%temperature, r4)
       write(ifnr) real(observations(n)%airmass, r4)
       write(ifnr) real(observations(n)%zoflev, r4)
       write(ifnr) real(observations(n)%eqlat, r4)
       write(ifnr) real(observations(n)%pvu, r4)
       write(ifnr) observations(n)%tph_lev
    enddo


    !// --------------------------------------------------------------------
  end subroutine caribic_data_to_file
  !// --------------------------------------------------------------------




  !// --------------------------------------------------------------------
  subroutine caribic2_master(JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI,LNEWM)
    !// --------------------------------------------------------------------
    !// Process CARIBIC profiles. Called outside parallell region.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI
    logical, intent(in) :: LNEWM
    !// --------------------------------------------------------------------

    if (NDAY.ne.NDAYI .and. NMET.eq.1 .and. NOPS.eq.1) then
       !// Finish end of last day
       call caribic_output(curyear,curmonth,curdate,NDAY-1, NRMETD+1, 1)

       !// Write to file (NDAY-1 since we save previous day)
       call caribic_data_to_file(NDAY-1)

    endif

    !// Initialize output for vertical profiles from CARIBIC
    if (nmet.eq.1 .and. nops.eq.1) then
       !// Initialize
       call get_new_events(JYEAR,JMON,JDATE,NDAY,NDAYI, NMET,NOPS,LNEWM)

       !// Set date for output file/info
       curyear = JYEAR
       curmonth= JMON
       curdate = JDATE
    endif

    !// Process profiles
    call caribic_output(JYEAR,JMON,JDATE,NDAY, NMET, NOPS)


    !// --------------------------------------------------------------------
  end subroutine caribic2_master
  !// --------------------------------------------------------------------




  !// --------------------------------------------------------------------
end module caribic2
!// --------------------------------------------------------------------

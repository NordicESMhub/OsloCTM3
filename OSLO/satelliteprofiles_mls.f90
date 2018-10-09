!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Satellite profiles MLS
!//=========================================================================
module satelliteprofiles_mls
  !// ----------------------------------------------------------------------
  !// MODULE: satelliteprofiles_mls
  !// DESCRIPTION: Satellite (MLS) measurement locations are fed into the
  !//              CTM on a daily basis. Spatial interpolation to satellite
  !//              location is done on-line, to save output file size.
  !// ----------------------------------------------------------------------
  !// VERSION 3
  !// 
  !// Satellite (e.g. MLS) measurement locations are fed into the CTM on
  !// a daily basis. Spatial interpolation to satellite location is done
  !// on-line, to save output file size.
  !//
  !// Also put out is the fractional distance to grid box centers as
  !// explained below.
  !//
  !// Input is read each day, and output is saved at the end of each day.
  !//
  !// Data for the ***two closest hours*** are retrieved, and because of this
  !// the routine is outside MP-blocks (for simplicity).
  !//
  !// VERSION 2 does horizontal interpolation on-line to save output file size.
  !//
  !// Ole Amund Sovde, May 2010
  !// ----------------------------------------------------------------------
  !//
  !// How to interpolate the data:
  !// The easiest is to interpolate linearly from model grid to station
  !// location. However, since the pressure may differ in neighboring boxes
  !// at least below the layer where ETAB is non-zero, vertical interpolation
  !// should be carried out first.
  !//
  !// Having the grids (iii,jjj) of the station and the closest neighboring
  !// grids ((iii,njj),(nii,jjj),(nii,njj)), as well as the fractional
  !// distance from the grid center to the location, we can do the
  !// interpolation.
  !//
  !// Method:
  !//  |---------|---------|
  !//  |         |         |
  !//  |    3    |    4    |   - Interpolate in x-direc. for ind_j and ind_nbj
  !//  |         |         |   - Use the results and interpolate in y-direction
  !//  |---------|---------|   - 1: iii,jjj
  !//  |      s  |         |     2: nii,jjj
  !//  |    1    |    2    |     3: iii,njj
  !//  |         |         |     4: nii,njj
  !//  |---------|---------|
  !//       --- distance (degrees) from to station to grid center: ds_x
  !//                              (ds_y for latitude)
  !//       ----------- distance between grid centers: dx (dy for latitude)
  !//
  !// A station (s) in grid 1, with neighbours 2, 3 and 4 is a distance ds_x
  !// from the grid center of grid 1. Calling the grid center distance dx, the
  !// fractional distance from the grid center is ds_x/dx.
  !//
  !// In a linear interpolation, the fractional distance ds_x/dx gives the
  !// fraction of the contribution of the neighbouring box to the interpolated
  !// value.
  !//   xfrac = ds_x/dx
  !//   yfrac = ds_y/dy
  !// The contribution from the grid box itself is therefore (1 - ds_x/dx).
  !// Using this method, we get from a field f
  !//   tmpf_jjj = f(iii,jjj)*(1 - xfrac) + f(nii,jjj)*xfrac
  !//   tmpf_njj = f(iii,njj)*(1 - xfrac) + f(nii,njj)*xfrac
  !// and then the interpolated is found by using these in the y-direction:
  !//   value =  tmpf_jjj*(1-yfrac) + tmpf_njj*yfrac
  !//
  !// The shortest way is to do everything at once
  !//   value = (f(iii,jjj)*(1 - xfrac) + f(nii,jjj)*xfrac) * (1-yfrac)
  !//         + (f(iii,njj)*(1 - xfrac) + f(nii,njj)*xfrac) * yfrac
  !//
  !//
  !// Ole Amund Sovde, February 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, &
       LOSLOCSTRAT
  use cmn_ctm, only: AREAXY, NRMETD, NROPSM, STT, AIR, &
       MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, ETAA, ETAB
  use cmn_met, only: ZOFLE, P, T, Q, BLH_CUR, BLH_NEXT, PVU
  use cmn_oslo, only: XSTT, trsp_idx, Xtrsp_idx, LMTROP
  use physics_oslo, only: NTHE, pvthe, theqlat
  use strat_h2o, only: sumH2
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// STT components
  !// MLS: O3, HNO3, CO, HO2, OH, SO2, CH3Cl, N2O, HCl, H2O, ClO, HOCl, BrO, NOy, T
  !// MLS:  1,    4,  6,  21, 40,  72,   106, 107, 111, 114, 133, 134,  139, 147
  integer,parameter                    :: ntracer = 13
  integer,dimension(ntracer),parameter :: stt_components=(/1,4,6,21,40,106,107,111,124,133,134,139,147/)
  real(r8), dimension(ntracer)           :: mole_mass

  !// Version
  integer, parameter :: version=3
  !// 3: Some variables are now double precision (etaa, etab, mole_mass)
  real(r8) :: version_mls

  !// Station data
  !// --------------------------------------------------------------------
  !// Store hourly data for 24 hours
  type measurement
     integer                    :: date        !// YYYYMMDD
     integer                    :: time        !// HHMMSS
     integer                    :: time_sec    !// seconds of day
     integer                    :: time_sec_org !// seconds of day
     real(r4)                   :: lat         !// latitude of profile
     real(r4)                   :: lon         !// longitude of profile
     integer                    :: ii        !// longitudinal grid box index
     integer                    :: jj        !// latitudinal grid box index
     integer                    :: nb_ii     !// neighbor lon. grid box index
     integer                    :: nb_jj     !// neighbor lat. grid box index
     real(r4)                   :: nb_xfrac
     real(r4)                   :: nb_yfrac
     real(r8)                   :: areaxy        !// Surface area
     integer,dimension(2)       :: ctmtime      !// CTM time [HHMMSS]
     integer,dimension(2)       :: ctmtime_sec  !// CTM time [sec]
     !// For each field, store data from the two closest NOPS
     !// 1: prior to measurement, 2: after
     real(r8),dimension(2)      :: psfc          !// Surface pressure
     real(r8),dimension(2)      :: blh           !// Boundary layer height
     real(r8),dimension(LPAR,ntracer,2)  :: mass !// Tracer mass
     real(r8),dimension(LPAR,2)  :: h2o          !// H2O (not part of STT)
     real(r8),dimension(LPAR,2)  :: temperature  !// Temperature
     real(r8),dimension(LPAR,2)  :: airmass      !// Air mass
     real(r8),dimension(LPAR+1,2):: zoflev       !// Height of box bottoms (1 is topography)
     real(r8),dimension(NTHE,2)  :: eqlat        !// Equivalent latitude
     real(r8),dimension(LPAR,2)  :: pvu          !// Potential vorticity units
     real(r8),dimension(2)       :: tph_pres     !// Tropopause pressure
  end type measurement


  !// Max profile measurements per day
  integer,parameter    :: nsp_max=3500

  !// Hourly data for each station
  type(measurement), dimension(nsp_max) :: measurements

  !// Areas of all neighbor grid boxes are needed
  real(r8), dimension(4,nsp_max) :: area_gridboxes

  !// The actual number of profiles
  integer              :: nsatprofs

  !// Date information for the day in progress
  integer :: fileyear, filedate, filemonth

  !// Number of stations per MP-block
  !integer,dimension(MPBLK)  :: mpstat
  !// Corresponding profile and neighbor position
  !// Assume roughly evenly distributed profiles around the globe, so
  !// nsp_max should be ok for these (if all profiles in one MP-block
  !// we would need nsp_max*4).
  !integer,dimension(nsp_max,MPBLK)  :: mpnsp, mpnbb

  !// Has the arrays been initialized?
  logical :: LINITIALIZED = .false.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'satelliteprofiles_mls.f90'
  !// ----------------------------------------------------------------------
  !// All variables are to be saved.
  save
  !// All is private
  private
  !// except
  public satprofs_mls_master
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine get_satprofiles_mls(NDAY, NMET, NOPS)
    !// --------------------------------------------------------------------
    !// For given locations save model profiles for the two closest operator
    !// split steps (NOPS). Save the four closest columns.
    !// Write to file at the end of the day.
    !//
    !// Input is mass.
    !// Carried out outside IJ-blocks.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input variables
    integer,intent(in) :: NDAY, NMET, NOPS

    !// Local variables
    integer :: nsp, nbb, i, j, L, n, NN, tops, &
         ii, jj, nii, njj
    integer :: sectime, sectime_pnops, sectime_mnops,mintime, maxtime, &
         ctmtime, ctmtime_pnops, hhmmsstime, hhmmsstime_p
    real(r8) :: time_now, time_next, time_nmet, time_nops
    real(r8) :: mass_CH4, mass_H2, H2O_TMP(4,LPAR), &
         vs, vp(LPAR), vpp1(LPAR+1), tp(4), &
         xfrac, yfrac, t1, t2, t3, frac_next

    logical :: LFIRST, LLAST

    !// For converting H2O from concentration to mass
    real(r8),parameter :: &
         C2MH2O  = 18._r8/M_AIR, &
         VMR2MMRH2 = 2._r8/M_AIR, &
         MW_H2 = 2._r8, MW_CH4 = 16._r8, MW_H2O = 18._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_satprofiles_mls'
    !// --------------------------------------------------------------------

    !// Arrays must have been initialized
    if (.not.LINITIALIZED) then
       if (NMET.eq.1.and.NOPS.eq.1) &
            write(6,'(a)') '* '//f90file//':'//subr// &
            ': No satellite measurements to simulate!'
       return
    endif


    !// Time steps for NMET and NOPS
    time_nmet = 24._r8/real(NRMETD, r8)
    time_nops = time_nmet/real(NROPSM, r8)
    frac_next = real(NOPS-1, r8) / real(NROPSM) !// For BLH temp.interp.

    !// Time goes from 000000 to 235959, depending on NROPSM.
    !// When NROPSM=3, time goes to 230000.
    time_now = real(NMET - 1, r8)*time_nmet + real(NOPS-1, r8)*time_nops

    !// Will save data for the two closest model hours (plus minus one hour).

    !// CTM time in seconds
    sectime = nint(time_now*3600._r8)
    sectime_pnops = nint( (time_now + time_nops) * 3600._r8 )
    sectime_mnops = nint( (time_now - time_nops) * 3600._r8 )

    !// CTM time in HHMMSS
    t1 = nint(time_now)
    t2 = nint((time_now - real(t1, r8))*60._r8)
    t3 = nint(((time_now - real(t1, r8))*60._r8 - real(t2, r8))*60._r8)
    hhmmsstime = t1*10000 + t2*100 + t3
    t1 = nint(time_now + time_nops)
    t2 = nint((time_now + time_nops - real(t1, r8))*60._r8)
    t3 = nint(((time_now + time_nops - real(t1, r8))*60._r8 &
                                     - real(t2, r8))*60._r8)
    hhmmsstime_p = t1*10000 + t2*100 + t3

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


    do nsp = 1, nsatprofs

      !// Measurements between NOPS (tnops=1) and NOPS+1 (tnops=2)
      !// are to be stored.
      do tops = 1, 2
         
        !// NOPS (tops=1): done each time step, but not for NMET==NRMETD+1.
        !//                NMET==NRMETD+1 is only for catching the last
        !//                time step of the day.
        !// NOPS+1 (tops=2): Cannot be stored for NOPS==1, so it must be stored
        !//                  at the next NOPS. Hence, at the next NOPS we must
        !//                  store data for the _previous_ NOPS.
        !//                  This is also done for NMET==NRMETD+1.
        !// Flags: LLAST: NMET==NRMETD+1 and NOPS==1
        !//        LFIRST: NOPS==1 and NMET==1
        if ( (tops.eq.1 .and. (.not. LLAST)) .or. &
             (tops.eq.2 .and. (.not. LFIRST)) ) then

          !// Possible routes:
          !// tops:1 LLAST:false: Put out data for this NOPS
          !// tops:2 LFIRST:false: Put out data data measured last NOPS
          !//        which should not be done for NOPS=1
          !// other: Nothing to be done

          if (tops .eq. 1) then
            !// Time interval NOPS:NOPS+1
            mintime = sectime
            maxtime = sectime_pnops
          else if (tops .eq. 2) then
            !// Time interval NOPS-1:NOPS, since we are storing data for the
            !// previous NOPS (NOPS has been increased, and we need the NOPS+1
            !// values for the last NOPS).
            !// I.e. for a given NOPS: sectime_pnops(NOPS) == sectime(NOPS+1)
            !//                        sectime_mnops(NOPS+1) == sectime(NOPS)
            mintime = sectime_mnops
            maxtime = sectime
          else
            write(6,'(a,i5)') '* '//f90file//':'//subr// &
                  ': tops wrong',tops
            stop
          endif

          !// Check the time
          if (measurements(nsp)%time_sec .ge. mintime .and. &
              measurements(nsp)%time_sec .lt. maxtime) then
            !// Save CTM time stamps
            if (tops.eq.1) then
               measurements(nsp)%ctmtime_sec(1) = sectime
               measurements(nsp)%ctmtime_sec(2) = sectime_pnops
               measurements(nsp)%ctmtime(1) = hhmmsstime
               measurements(nsp)%ctmtime(2) = hhmmsstime_p
            endif

            !// 1. Find values for neighbor columns
            !// 2. Interpolate horizontally
            ii  = measurements(nsp)%ii
            jj  = measurements(nsp)%jj
            nii = measurements(nsp)%nb_ii
            njj = measurements(nsp)%nb_jj
            xfrac = measurements(nsp)%nb_xfrac
            yfrac = measurements(nsp)%nb_yfrac



            !// Surface pressure
            tp(1) = P( ii, jj)
            tp(2) = P(nii, jj)
            tp(3) = P( ii,njj)
            tp(4) = P(nii,njj)
            !// Interpolate x- and y-direction.
            vs = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            measurements(nsp)%psfc(tops)           = vs

            !// Boundary layer height
            !// Interpolate in time between BLH_CUR and BLH_NEXT
            !// When NOPS runs from 1 to 3 (NROPSM=3):
            !//   NOPS=1: weight 1 for current, 0 for next.
            !//   NOPS=2: weight 2/3 for current, 1/3 for next.
            !//   NOPS=3: weight 1/3 for current, 2/3 for next.
            !//   Fraction for current weight is:
            !//     frac_cur = 1 - (NOPS - 1)/NROPSM
            !//   Fraction for next weight is:
            !//     frac_next = (NOPS - 1)/NROPSM
            tp(1) = (1._r8 - frac_next) * BLH_CUR( ii, jj) &
                           + frac_next * BLH_NEXT( ii, jj)
            tp(2) = (1._r8 - frac_next) * BLH_CUR(nii, jj) &
                            + frac_next * BLH_NEXT(nii, jj)
            tp(3) = (1._r8 - frac_next) * BLH_CUR( ii,njj) &
                            + frac_next * BLH_NEXT( ii,njj)
            tp(4) = (1._r8 - frac_next) * BLH_CUR(nii,njj) &
                            + frac_next * BLH_NEXT(nii,njj)
            !// Interpolate x- and y-direction.
            vs = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            measurements(nsp)%blh(tops)            = vs


            !// Calculate H2O mass
            do nbb = 1, 4
               if (nbb.eq.1) then
                  i = measurements(nsp)%ii
                  j = measurements(nsp)%jj
               elseif (nbb.eq.2) then
                  i = measurements(nsp)%nb_ii
                  j = measurements(nsp)%jj
               elseif (nbb.eq.3) then
                  i = measurements(nsp)%ii
                  j = measurements(nsp)%nb_jj
               elseif (nbb.eq.4) then
                  i = measurements(nsp)%nb_ii
                  j = measurements(nsp)%nb_jj
               else
                  write(6,'(a,i5)') '* '//f90file//':'//subr// &
                       ': NBB is wrong!',nbb
                  stop
               endif
               !// Set from meteorology
               do L = 1, LPAR
                  H2O_TMP(nbb,L) = Q(i,j,L) * AIR(i,j,L)
               enddo
               !// Calculate in stratosphere
               if (LOSLOCSTRAT) then
                 if (trsp_idx(114) .gt. 0) then
                   !// H2O is transported (use only stratosphere)
                   do L = LMTROP(I,J)+1, LPAR
                     H2O_TMP(nbb,L) = STT(I,J,L,trsp_idx(114))
                   end do
                 else
                   !// Calculate stratospheric H2O from sum of H2
                   do L = LMTROP(I,J)+1, LPAR
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
                     H2O_TMP(nbb,L) = ( sumH2 * &
                          VMR2MMRH2 * AIR(I,J,L) / MW_H2 & !sumH2 [kmol/gbox]
                          - 2._r8 * mass_CH4 / MW_CH4 &  ! -2*CH4 [kmol/gbox]
                          - mass_H2  / MW_H2  &          ! - H2 [kmol/gbox]
                          ) * MW_H2O !// From kmol to kg
                   end do
                 endif
               endif !// if (LOSLOCSTRAT) then
            enddo !// do nbb = 1, 4

            do L = 1, LPAR
               tp(1) = H2O_TMP(1,L)
               tp(2) = H2O_TMP(2,L)
               tp(3) = H2O_TMP(3,L)
               tp(4) = H2O_TMP(4,L)
               !// Interpolate.
               vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                       + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo !// do L = 1, LPAR
            measurements(nsp)%h2o(:,tops)          = vp(:)


            !// Tracers
            do n = 1, ntracer
               do L = 1, LPAR
                  !// Transported
                  NN = trsp_idx(stt_components(n))
                  if (NN.gt.0) then
                     tp(1) = STT( ii, jj,L,NN)
                     tp(2) = STT(nii, jj,L,NN)
                     tp(3) = STT( ii,njj,L,NN)
                     tp(4) = STT(nii,njj,L,NN)
                  else
                     !// Non-transported (checked in read-in)
                     NN = Xtrsp_idx(stt_components(n))
                     tp(1) = XSTT(L,NN, ii, jj)
                     tp(2) = XSTT(L,NN,nii, jj)
                     tp(3) = XSTT(L,NN, ii,njj)
                     tp(4) = XSTT(L,NN,nii,njj)
                  endif
                  !// Interpolate.
                  vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                          + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
               enddo !// do L = 1, LPAR
               measurements(nsp)%mass(:,n,tops) = vp(:)
            enddo !// do n = 1, ntracer


            !// Temperature
            do L = 1, LPAR
               tp(1) = T( ii, jj,L)
               tp(2) = T(nii, jj,L)
               tp(3) = T( ii,njj,L)
               tp(4) = T(nii,njj,L)
               !// Interpolate.
               vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                       + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo !// do L = 1, LPAR
            measurements(nsp)%temperature(:,tops)  = vp(:)


            !// Air mass
            do L = 1, LPAR
               tp(1) = AIR( ii, jj,L)
               tp(2) = AIR(nii, jj,L)
               tp(3) = AIR( ii,njj,L)
               tp(4) = AIR(nii,njj,L)
               !// Interpolate.
               vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                       + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo
            measurements(nsp)%airmass(:,tops)      = vp(:)


            !// Height of box bottom
            do L = 1, LPAR+1
               tp(1) = ZOFLE(L, ii, jj)
               tp(2) = ZOFLE(L,nii, jj)
               tp(3) = ZOFLE(L, ii,njj)
               tp(4) = ZOFLE(L,nii,njj)
               !// Interpolate.
               vpp1(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                        + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo
            measurements(nsp)%zoflev(:,tops)       = vpp1(:)


            !// Equivalent latitudes on NTHE theta levels
            do L = 1, NTHE
               tp(1) = theqlat( ii, jj,L)
               tp(2) = theqlat(nii, jj,L)
               tp(3) = theqlat( ii,njj,L)
               tp(4) = theqlat(nii,njj,L)
               !// Interpolate.
               vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                       + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo
            measurements(nsp)%eqlat(1:NTHE,tops)    = vp(1:NTHE)


            !// Potential voritcity units (PVU)
            do L = 1, LPAR
               tp(1) = PVU(L, ii, jj)
               tp(2) = PVU(L,nii, jj)
               tp(3) = PVU(L, ii,njj)
               tp(4) = PVU(L,nii,njj)
               !// Interpolate.
               vp(L) = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                       + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            enddo
            measurements(nsp)%pvu(:,tops)          = vp(:)


            !// Tropopause height pressure (bottom edge of LMTROP+1)
            tp(1) = etaa( LMTROP( ii, jj)+1 ) + etab( LMTROP( ii, jj)+1 )*p( ii, jj)
            tp(2) = etaa( LMTROP(nii, jj)+1 ) + etab( LMTROP(nii, jj)+1 )*p(nii, jj)
            tp(3) = etaa( LMTROP( ii,njj)+1 ) + etab( LMTROP( ii,njj)+1 )*p( ii,njj)
            tp(4) = etaa( LMTROP(nii,njj)+1 ) + etab( LMTROP(nii,njj)+1 )*p(nii,njj)
            !// Interpolate.
            vs = ( (tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8 - yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
            measurements(nsp)%tph_pres(tops)        = vs


          endif !// if (measurements(nsp)%sectime ...

        endif !// if ( (tops.eq.1 .and. ...
      enddo !// do tops = 1, 2
    enddo !// do nsp = 1, nsatprofs



    !// --------------------------------------------------------------------
  end subroutine get_satprofiles_mls
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine initialize_satprofiles_mls(JYEAR,JMON,JDATE)
    !// --------------------------------------------------------------------
    !// Initializes the satellite list based on a input list.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: XDEDG, YDEDG, YDGRD, XDGRD
    use cmn_chem, only: TMASS
    use cmn_oslo, only: XTMASS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer,intent(in)  :: jyear, jmon, jdate

    !// Local variables
    integer           :: i,j,l, n, m
    integer           :: iii, jjj, nii, njj, nsp
    integer           :: hh,mm,ss
    integer :: file_nr,file_err
    real(r8)  :: ctmlonb(IPAR+1),ctmlon(IPAR),deltalon,deltalat
    real(r8)  :: tp(4),rtmp
    logical :: file_io, linread
    character(len=80) :: filename
    character(len=8) :: cdate

    !// Profile data
    integer,dimension(nsp_max)  :: sat_time, sat_date,sat_sectime
    integer,dimension(nsp_max)  :: ind_i, ind_j     ! grid box
    integer,dimension(nsp_max)  :: ind_nbi,ind_nbj  ! closest neighbor boxes
    real(r8),dimension(nsp_max)   :: vlats,vlons      ! lat/lon
    real(r8),dimension(nsp_max)   :: nbxfrac,nbyfrac  ! fraction of neighbor box
    real(r8),dimension(nsp_max)   :: axi              ! interpolated grid area
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'initialize_satprofiles_mls'
    !// --------------------------------------------------------------------

    write(*,'(a)')'* satelliteprofiles_mls.f90: Initializing daily MLS '// &
         'satellite measurements'

    !// Check you have the tracers available
    do N = 1, ntracer
       if ( .not. (trsp_idx(stt_components(n)).gt.0 .or. &
            Xtrsp_idx(stt_components(n)).gt.0) ) then
          write(6,'(a,i5,a)') '* '//f90file//':'//subr// &
               ': Tracer number ',stt_components(n), &
               ' is not included in the simulation!'
          stop
       endif
    enddo

    !// Get tracer molecular masses
    do N = 1, ntracer
       if ( trsp_idx(stt_components(n)).gt.0) then
          mole_mass(N) = TMASS(trsp_idx(stt_components(n)))
       else
          mole_mass(N) = XTMASS(Xtrsp_idx(stt_components(n)))
       endif
    enddo

    !// The profiles list
    write(cdate(1:8),'(i4.4,2i2.2)') jyear,jmon,jdate 
    filename='MLSv3_LOCATIONS/MLSv3_locations_'//cdate//'.dat'

    !// Find a file number
    file_nr=20
    file_io=.true.
    do while (file_io)
       file_nr=file_nr+1
       inquire(file_nr,opened=file_io)
    enddo
    open(file_nr,file=filename,form='formatted',status='old',iostat=file_err)
    if (file_err .ne. 0) then
       write(*,'(a,a)') '*** Could not open file: ',trim(filename)
       write(*,'(a)')   '    Will not simulate satellite profiles!'
       nsatprofs = 0
       LINITIALIZED = .false.
       return
    endif

    !// Read header
    read(file_nr,*)
    read(file_nr,'(9x,f5.2)',iostat=file_err) version_mls
    if (file_err .ne. 0) then
       !// Failsafe for version 2 without header
       version_mls=2.2
       backspace(file_nr)
       backspace(file_nr)
    endif

    !// Read the station list
    nsp = 0
    linread = .true.
    do while (linread)
       nsp = nsp + 1 !// Count the profiles
       read(file_nr,'(i8,3i2,f9.4,1x,f9.4)',iostat=file_err) &
            sat_date(nsp),hh,mm,ss, vlons(nsp),vlats(nsp)
       if (file_err .eq. 0) then
          !// Time [hh:mm:ss]
          sat_time(nsp) = hh*10000 + mm*100 + ss
          !// Time [seconds of day]
          sat_sectime(nsp) = hh*3600 + mm*60 + ss
          !// Allow 10 seconds past midnight; treated as 235959 below.
          if ((sat_sectime(nsp) - 864000).gt. 10) then
             write(6,'(a,i5)') '* '//f90file//':'//subr// &
                  ': Problem with time:',nsp
             write(6,'(a,i10)')'  sat_sectime',sat_sectime(nsp)
             write(6,'(a,i2.2,a1,i2.2,a1,i2.2)')'  hh:mm:ss',hh,':',mm,':',ss
             stop
          endif
          !// If data was read, update the number of stations
          nsatprofs = nsp
       else
          linread = .false.
       endif
    enddo
    !// Check number of stations
    if (nsatprofs .gt. nsp_max) then
       write(6,'(a)') '* '//f90file//':'//subr// &
            ': Too many measurements on file!'
       write(6,'(a,i4)') '   max_stations: ',nsp_max
       write(6,'(a,i4)') '   nrofstations: ',nsatprofs
       stop
    endif

    !// Close file
    close(file_nr)

    write(*,'(a,i4)') '  Number of satprofs read from file: ',nsatprofs
    !// --------------------------------------------------------------------

    !// Get longitude boundary
    ctmlonb(1) = XDEDG(1)
    do i = 2, ipar+1
       if (XDEDG(i) .lt. 0._r8) then
          ctmlonb(i) = XDEDG(i) + 360._r8
       else
          ctmlonb(i) = XDEDG(i)
       end if
    end do
    !// Get longitude center
    do i = 1, ipar
       if (XDGRD(i) .lt. 0._r8) then
          ctmlon(i) = XDGRD(i) + 360._r8
       else
          ctmlon(i) = XDGRD(i)
       end if
    end do

    !// Initialize
    iii = -1
    jjj = -1

    !// Find grid boxes J-index
    do nsp = 1, nsatprofs
       do j = 1, JPAR
          if ( (vlats(nsp).ge.YDEDG(j)) .and. (vlats(nsp).lt.YDEDG(j+1)) ) then
             !// We have y-index
             jjj=j
             exit
          endif
       enddo
       if (vlats(nsp).ge.YDEDG(jpar+1)) jjj=jpar !// special for NP
       if (vlats(nsp).lt.YDEDG(1)) jjj=1         !// special for SP
       !// Update array
       ind_j(nsp) = jjj
    enddo

    !// Find grid boxes I-index
    do nsp = 1, nsatprofs
       !// Longitude on degrees in the range <ctmlonb(1),ctmlonb(ipar+1)]
       if (vlons(nsp) .gt. ctmlonb(ipar+1)) vlons(nsp) = vlons(nsp) - 360._r8
       if (vlons(nsp) .le. ctmlonb(1)) vlons(nsp) = vlons(nsp) + 360._r8
       do I = 1, IPAR
          if ( (vlons(nsp) .ge. ctmlonb(i)) .and. &
               (vlons(nsp) .lt. ctmlonb(i+1)) ) then
             !// We have x-index
             iii=i
             exit
          endif
       enddo
       !// Update array
       ind_i(nsp)=iii
    enddo

    if (iii.lt.0 .or. jjj.lt.0) then
       write(6,'(a,2i5)') '* '//f90file//':'//subr// &
            ': ERROR: III/JJJ is wrong:',iii,jjj
       stop
    endif

    !// Find neighbour boxes: closest i-box
    do nsp = 1, nsatprofs
       !// Double check overlap 0E
       if (ind_i(nsp) .eq. 1 .and. vlons(nsp) .gt. ctmlonb(ipar) &
            .and. vlons(nsp) .le. ctmlonb(ipar+1)) then
          !// This should not be necessary if the above works
          !// Neighbor grid box is just to the west of box 1
          rtmp = vlons(nsp) - 360._r8
          write(6,'(a,i7,f12.3)') f90file//':'//subr//': correcting vlons', &
               nsp,vlons(nsp)
       else
          rtmp = vlons(nsp)
       endif
       if (rtmp .ge. ctmlon(ind_i(nsp))) then
          !// east of grid box center
          nii = ind_i(nsp) + 1
       else
          !// west of grid box center
          nii = ind_i(nsp) - 1
       endif
       !// check the "boundaries"
       if (nii .eq. (ipar+1)) nii = 1
       if (nii .eq. 0) nii = ipar
       !// update array
       ind_nbi(nsp)=nii
    enddo

    !// Find neighbour boxes: closest j-box
    do nsp = 1, nsatprofs
       if (vlats(nsp) .gt. YDGRD(ind_j(nsp))) then
          njj = ind_j(nsp) + 1
       else
          njj = ind_j(nsp) - 1
       endif
       !// Check the boundaries
       if (njj .eq. 0) njj=1         !// special for SP
       if (njj .eq. jpar+1) njj=jpar !// special for NP
       !// Update array
       ind_nbj(nsp) = njj
    enddo


    !// Find fraction of neighbour boxes
    deltalon = ctmlon(2) - ctmlon(1)
    do nsp = 1, nsatprofs
       iii = ind_i(nsp)
       jjj = ind_j(nsp)
       if (jjj .gt. 1) then
          deltalat = YDGRD(jjj) - YDGRD(jjj-1)
       else
          deltalat = YDGRD(2) - YDGRD(1) ! special for SP
       endif
       nbyfrac(nsp) = abs(vlats(nsp) - YDGRD(jjj))/deltalat
       nbxfrac(nsp) = abs(vlons(nsp) - ctmlon(iii))/deltalon
       !// Check fractions
       if (nbxfrac(nsp) .gt. 1._r8) then
          write(6,'(a)') '* '//f90file//':'//subr// &
               ': Something wrong with lon:'
          print*,nbxfrac(nsp),nbyfrac(nsp)
          print*,vlons(nsp),ctmlon(iii),XDGRD(iii)
          print*,iii,ind_nbi(nsp)
          stop
       endif
       if (nbyfrac(nsp) .gt. 1._r8) then
          write(6,'(a)') '* '//f90file//':'//subr// &
               ': Something wrong with lat:'
          print*,nbxfrac(nsp),nbyfrac(nsp)
          print*,vlats(nsp),YDGRD(jjj)
          print*,jjj,ind_nbj(nsp)
          stop
       endif
    enddo

    !// Find mean grid box area
    do nsp = 1, nsatprofs
       tp(1) = AREAXY(ind_i(nsp),  ind_j(nsp))
       tp(2) = AREAXY(ind_nbi(nsp),ind_j(nsp))
       tp(3) = AREAXY(ind_i(nsp),  ind_nbj(nsp))
       tp(3) = AREAXY(ind_nbi(nsp),ind_nbj(nsp))
       axi(nsp) = &
            ( (tp(1)*(1._r8 - nbxfrac(nsp)) + &
            + tp(2)*nbxfrac(nsp)) * (1._r8-nbyfrac(nsp)) &
            + (tp(3)*(1._r8 - nbxfrac(nsp)) + tp(4)*nbxfrac(nsp)) * nbyfrac(nsp))
    enddo


    !// --------------------------------------------------------------------

    !// Put into measurements array
    do nsp = 1, nsatprofs
       measurements(nsp)%date     = sat_date(nsp)    !// profile date
       measurements(nsp)%time     = sat_time(nsp)    !// profile time [hhmmss]
       !// Profile time in seconds of day. When comparing this to the MLS time
       !// stamp, you must subtract the seconds between 1/1-93 and your year.
       !// For a given year and counting the number of leap days from 1/1-93,
       !// subtract:   previous_secs=(365L*(year-1993)+leaps)*3600L*24L
       !// For other satellites this may be different.
       measurements(nsp)%time_sec_org = sat_sectime(nsp) !// profile time [sec]
       !// Unfortunately some data are for a few seconds over 240000.
       !// These are treated as 235959 in sat_sectime (which is not stored).
       measurements(nsp)%time_sec = min(86399,sat_sectime(nsp)) !// profile time [sec]
       measurements(nsp)%lat      = vlats(nsp)   !// Station latitude
       measurements(nsp)%lon      = vlons(nsp)   !// Station longitude
       measurements(nsp)%ii       = ind_i(nsp)   !// CTM i-index
       measurements(nsp)%jj       = ind_j(nsp)   !// CTM j-index
       measurements(nsp)%nb_ii    = ind_nbi(nsp) !// CTM closest i-neighbor
       measurements(nsp)%nb_jj    = ind_nbj(nsp) !// CTM closest j-neighbor
       measurements(nsp)%nb_xfrac = nbxfrac(nsp) !// Fractional distance x-direction
       measurements(nsp)%nb_yfrac = nbyfrac(nsp) !// Fractional distance y-direction
       measurements(nsp)%areaxy   = axi(nsp)
       measurements(nsp)%ctmtime(:)     = 0
       measurements(nsp)%ctmtime_sec(:) = 0
       !// Initialize the profiles
       measurements(nsp)%psfc(:)          = 0._r8
       measurements(nsp)%blh(:)           = 0._r8
       measurements(nsp)%mass(:,:,:)      = 0._r8
       measurements(nsp)%h2o(:,:)         = 0._r8
       measurements(nsp)%temperature(:,:) = 0._r8
       measurements(nsp)%airmass(:,:)     = 0._r8
       measurements(nsp)%zoflev(:,:)      = 0._r8
       measurements(nsp)%eqlat(:,:)       = 0._r8
       measurements(nsp)%pvu(:,:)         = 0._r8
       measurements(nsp)%tph_pres(:)      = 0._r8
    enddo

    !// We are initialized
    LINITIALIZED = .true.

    !// --------------------------------------------------------------------
  end subroutine initialize_satprofiles_mls
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine satprofiles_mls_to_file()
    !// --------------------------------------------------------------------
    !// Write satellite profiles to file.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: RESULTDIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    character(LEN=80) :: filename
    character(LEN=8) :: datestamp
    logical :: fnr_ok
    integer :: ifnr, nsp
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'satprofiles_mls_to_file'
    !// --------------------------------------------------------------------

    !// If no data, then no file
    if (nsatprofs .eq. 0) return

    !// Write file for previous day
    write(datestamp(1:8),'(i4.4,2i2.2)') fileyear,filemonth,filedate
    filename = trim(RESULTDIR)//'mls_profiles_'//datestamp//'.dta'

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
    write(ifnr) fileyear,filemonth,filedate,version,version_mls  !// Date & version
    write(ifnr) lpar,ntracer,nsatprofs, nthe,ipar,jpar !// Resolution
    write(ifnr) etaa, etab             !// Sigma coordinates
    write(ifnr) stt_components         !// Components (STT)
    write(ifnr) mole_mass              !// Molecular masses of components
    write(ifnr) pvthe                  !// Theta levels for eqlat
    do nsp = 1, nsatprofs
       write(ifnr) measurements(nsp)%time     !// HHMMSS
       write(ifnr) measurements(nsp)%time_sec_org  !// sec, the original (see init)
       write(ifnr) measurements(nsp)%ctmtime  !// HHMMSS
       write(ifnr) measurements(nsp)%ctmtime_sec  !// sec
       write(ifnr) measurements(nsp)%lat
       write(ifnr) measurements(nsp)%lon
       write(ifnr) measurements(nsp)%ii
       write(ifnr) measurements(nsp)%jj
       write(ifnr) measurements(nsp)%nb_ii
       write(ifnr) measurements(nsp)%nb_jj
       write(ifnr) measurements(nsp)%nb_xfrac
       write(ifnr) measurements(nsp)%nb_yfrac
       write(ifnr) real(measurements(nsp)%areaxy, r4)
       write(ifnr) real(measurements(nsp)%psfc, r4)
       write(ifnr) real(measurements(nsp)%blh, r4)
       write(ifnr) real(measurements(nsp)%mass, r4)
       write(ifnr) real(measurements(nsp)%h2o, r4)
       write(ifnr) real(measurements(nsp)%temperature, r4)
       write(ifnr) real(measurements(nsp)%airmass, r4)
       write(ifnr) real(measurements(nsp)%zoflev, r4)
       write(ifnr) real(measurements(nsp)%eqlat, r4)
       write(ifnr) real(measurements(nsp)%pvu, r4)
       write(ifnr) real(measurements(nsp)%tph_pres, r4)
    enddo
    close(ifnr)

    !// Done for now; Re-initialize counter and flag
    write(6,'(a,i6)') '* '//f90file//':'//subr// &
         ': Wrote satellite profiles: ',nsatprofs
    nsatprofs = 0
    LINITIALIZED = .false.

    !// --------------------------------------------------------------------
  end subroutine satprofiles_mls_to_file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine satprofs_mls_master(JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI)
    !// --------------------------------------------------------------------
    !// Process satellite profiles. Called outside parallell region.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI
    !// --------------------------------------------------------------------

    if (NDAY.ne.NDAYI .and. NMET.eq.1 .and. NOPS.eq.1) then
       !// Finish end of last day
       call get_satprofiles_mls(NDAY-1, NRMETD+1, 1)

       !// Write to file
       call satprofiles_mls_to_file()

    endif

    !// Initialize output for vertical profiles from satellites
    if (nmet.eq.1 .and. nops.eq.1) then
       !// Initialize
       call initialize_satprofiles_mls(JYEAR,JMON,JDATE)

       !// Set date for output file/info
       fileyear = JYEAR
       filemonth= JMON
       filedate = JDATE
    endif

    !// Process profiles
    call get_satprofiles_mls(NDAY, NMET, NOPS)

    !// --------------------------------------------------------------------
  end subroutine satprofs_mls_master
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module satelliteprofiles_mls
!//=========================================================================

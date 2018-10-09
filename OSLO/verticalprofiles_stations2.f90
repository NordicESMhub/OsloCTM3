!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Vertical profiles at selected stations.
!//=========================================================================
module verticalprofiles_stations2
  !// ----------------------------------------------------------------------
  !// MODULE: verticalprofiles_stations2
  !// DECRIPTION: Produce vertical profiles at selected stations.
  !// ----------------------------------------------------------------------
  !// VERSION 2
  !//
  !// Reads a list of station names and locations. Spatial interpolation to
  !// satellite location is done on-line, to save output file size.
  !//
  !// Output is stored for all NOPS steps and saved at the end of each day.
  !//
  !// VERSION 2 does horizontal interpolation on-line to save output file size.
  !//
  !// Ole Amund Sovde, May 2010
  !// ----------------------------------------------------------------------
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
  !//  |         |         |   - Use the results and interpolate in y-direc.
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
  !// from the grid center of grid 1. Calling the grid center distance dx,
  !// the fractional distance from the grid center is ds_x/dx.
  !//
  !// In a linear interpolation, the fractional distance ds_x/dx gives the
  !// fraction of the contribution of the neighbouring box to the
  !// interpolated value.
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
  !// Ole Amund Sovde, February 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, &
       LOSLOCSTRAT
  use cmn_ctm, only: AREAXY, STT, AIR, ETAA, ETAB, NROPSM, NRMETD, &
       MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
  use cmn_met, only: ZOFLE, P, T, Q, BLH_CUR, BLH_NEXT, PVU
  use cmn_oslo, only: XSTT, trsp_idx, Xtrsp_idx, LMTROP
  use physics_oslo, only: NTHE, pvthe, theqlat
  use strat_h2o, only: sumH2
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// STT components
  integer,parameter                    :: ntracer = 16
  integer,dimension(ntracer),parameter :: stt_components=(/1,4,6,7,8,9,13,15,17,20,21,41,42,43,44,48/)
  !integer,parameter                    :: ntracer = 6
  !integer,dimension(ntracer),parameter :: stt_components=(/240,241,242,243,244,245/)
  real(r8), dimension(ntracer)         :: mole_mass

  !// Number of diagnoses per day. Must be >= NRMETD*NROPSM
  integer, parameter :: ndiag=96

  !// Version number
  integer, parameter :: version=2

  !// Station data
  !// ----------------------------------------------------------------------
  !// Store hourly data for 24 hours
  type measurement
     character(len=30)      :: locname     !// station name
     character(len=3)       :: loccode     !// station code
     real(r4)               :: lat         !// latitude of station
     real(r4)               :: lon         !// longitude of station
     real(r4)               :: alt         !// altitude of station
     integer                :: ii          !// longitudinal grid box index
     integer                :: jj          !// latitudinal grid box index
     integer                :: nb_ii       !// neighbor lon. grid box index
     integer                :: nb_jj       !// neighbor lat. grid box index
     real(r4)               :: nb_xfrac    !// fractional distance x-direction
     real(r4)               :: nb_yfrac    !// fractional distance y-direction
     real(r8)                 :: areaxy      !// Surface area
     real(r8),dimension(ndiag)     :: blh           !// Boundary layer height
     real(r8),dimension(ndiag)     :: psfc          !// Surface pressure
     real(r8),dimension(LPAR,ntracer,ndiag) :: mass !// Tracer mass
     real(r8),dimension(LPAR,ndiag)  :: h2o         !// H2O (not part of STT)
     real(r8),dimension(LPAR,ndiag)  :: temperature !// Temperature
     real(r8),dimension(LPAR,ndiag)  :: airmass     !// Air mass
     real(r8),dimension(LPAR+1,ndiag):: zoflev !// Height of box bottoms (1 is topography)
     real(r8),dimension(NTHE,ndiag)  :: eqlat  !// Equivalent latitude
     real(r8),dimension(LPAR,ndiag)  :: pvu    !// Potential vorticity units
     real(r8),dimension(ndiag)       :: tph_pres     !// Tropopause pressure
  end type measurement

  !// Max measurements per day
  integer, parameter    :: max_stations=250

  !// Hourly data for each station
  type(measurement), dimension(max_stations) :: measurements
  !// The actual number of stations
  integer :: nrofstations

  !// Date information for the day in progress
  integer :: fileyear, filedate, filemonth

  !// Has the arrays been initialized?
  logical :: LINITIALIZED = .false.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: &
       f90file = 'verticalprofiles_stations2.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public vprofs_master
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine vprof_stations(NDAY, NMET, NOPS)
    !// --------------------------------------------------------------------
    !// Save vertical profiles of tracer mass and meteorology data for
    !// pre-defined station locations. Profiles are saved each operator
    !// split loop (NOPS).
    !//
    !// Carried out outside IJ-blocks.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input variables
    integer,intent(in) :: NDAY, NMET, NOPS

    !// Local variables
    integer :: NDSTEP !// Diagnose step (1 to NROPSM*NRMETD)
    integer :: nst, nbb, i, j, L, n, NN, ii,jj,nii,njj

    !// For calculating H2O in stratosphere
    real(r8) :: mass_CH4, mass_H2, H2O_TMP(4,LPAR), &
         vs, vp(LPAR), vpp1(LPAR+1), tp(4), &
         xfrac, yfrac, frac_next

    !// For converting H2O from concentration to mass
    real(r8),parameter :: &
         C2MH2O  = 18._r8 / M_AIR, &
         VMR2MMRH2 = 2._r8 / M_AIR, &
         MW_H2 = 2._r8, MW_CH4 = 16._r8, MW_H2O = 18._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'vprof_stations'
    !// --------------------------------------------------------------------

    !// Arrays must have been initialized
    if (.not.LINITIALIZED) then
       if (NMET.eq.1 .and. NOPS.eq.1) &
            write(6,'(a)') f90file//':'//subr// &
            ': Not initialized; skipping.'
       return
    end if

    !// Diagnose this NOPS
    NDSTEP = (NMET - 1)*NROPSM + (NOPS - 1) + 1
    frac_next = real(NOPS-1, r8) / real(NROPSM) !// For BLH temp.interp.

    !// Save output for each station
    do nst = 1, nrofstations

       !// 1. Find values for neighbor columns
       !// 2. Interpolate horizontally
       ii  = measurements(nst)%ii
       jj  = measurements(nst)%jj
       nii = measurements(nst)%nb_ii
       njj = measurements(nst)%nb_jj
       xfrac = measurements(nst)%nb_xfrac
       yfrac = measurements(nst)%nb_yfrac



       !// Surface pressure
       tp(1) = P( ii, jj)
       tp(2) = P(nii, jj)
       tp(3) = P( ii,njj)
       tp(4) = P(nii,njj)
       !// Interpolate x- and y-direction.
       vs = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
           + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       measurements(nst)%psfc(NDSTEP)           = vs

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
       tp(1) = (1._r8 - frac_next) * BLH_CUR( ii, jj) + frac_next * BLH_NEXT( ii, jj)
       tp(2) = (1._r8 - frac_next) * BLH_CUR(nii, jj) + frac_next * BLH_NEXT(nii, jj)
       tp(3) = (1._r8 - frac_next) * BLH_CUR( ii,njj) + frac_next * BLH_NEXT( ii,njj)
       tp(4) = (1._r8 - frac_next) * BLH_CUR(nii,njj) + frac_next * BLH_NEXT(nii,njj)
       !// Interpolate x- and y-direction.
       vs = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
           + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       measurements(nst)%blh(NDSTEP)            = vs




       !// Calculate H2O mass
       do nbb = 1, 4
          if (nbb.eq.1) then
             i = measurements(nst)%ii
             j = measurements(nst)%jj
          else if (nbb.eq.2) then
             i = measurements(nst)%nb_ii
             j = measurements(nst)%jj
          else if (nbb.eq.3) then
             i = measurements(nst)%ii
             j = measurements(nst)%nb_jj
          else if (nbb.eq.4) then
             i = measurements(nst)%nb_ii
             j = measurements(nst)%nb_jj
          else
             write(6,'(a,i5)') f90file//':'//subr// &
                  ': NBB is wrong - stopping!',nbb
             stop
          end if
          !// Set from meteorology
          do L = 1, LPAR
             H2O_TMP(nbb,L) = Q(i,j,L) * AIR(i,j,L)
          end do
          !// Calculate in stratosphere
          if (LOSLOCSTRAT) then
             if (trsp_idx(114) .gt. 0) then
                !// H2O is transported (use only stratosphere)
                Do L = LMTROP(I,J)+1, LPAR
                   H2O_TMP(nbb,L) = STT(I,J,L,trsp_idx(114))
                End Do
             else
                !// Calculate stratospheric H2O from sum of H2
                Do L = LMTROP(I,J)+1, LPAR
                   !// Get CH4
                   if (trsp_idx(46).gt.0) then
                      mass_CH4 = STT(I,J,L,trsp_idx(46))
                   else if (Xtrsp_idx(46).gt.0) then
                      mass_CH4 = XSTT(L,Xtrsp_idx(46),I,J)
                   else
                      write(6,'(a)') f90file//':'//subr// &
                            ': CH4 not available - stopping'
                      stop
                   end if
                   !// Get H2
                   if (trsp_idx(113).gt.0) then
                      mass_H2 = STT(I,J,L,trsp_idx(113))
                   else if (Xtrsp_idx(113).gt.0) then
                      mass_H2 = XSTT(L,Xtrsp_idx(113),I,J)
                   else
                      write(6,'(a)') f90file//':'//subr// &
                            ': H2 not available - stopping'
                      stop
                   end if
                   !// Use moles to calculate stratospheric H2O mass
                   ! vmr: molec/molecair, molec = mass /M * Na
                   H2O_TMP(nbb,L) = ( sumH2 * &
                        VMR2MMRH2 * AIR(I,J,L) / MW_H2 & !// Sum H2 [kmol/gbox]
                        - 2._r8 * mass_CH4 / MW_CH4 &    !// -2*CH4 [kmol/gbox]
                        - mass_H2  / MW_H2  &            !// - H2 [kmol/gbox]
                        ) * MW_H2O !// From kmol to kg
                End Do
             end if
          end if !// if (LOSLOCSTRAT) then
       end do !// do nbb = 1, 4

       do L = 1, LPAR
          tp(1) = H2O_TMP(1,L)
          tp(2) = H2O_TMP(2,L)
          tp(3) = H2O_TMP(3,L)
          tp(4) = H2O_TMP(4,L)
          !// Interpolate x direction. Area weight for J-direction.
          vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
               + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do !// do L = 1, LPAR
       measurements(nst)%h2o(:,NDSTEP)          = vp(:)




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
                !// Non-transported if not transprorted (checked in read-in)
                NN = Xtrsp_idx(stt_components(n))
                tp(1) = XSTT(L,NN, ii, jj)
                tp(2) = XSTT(L,NN,nii, jj)
                tp(3) = XSTT(L,NN, ii,njj)
                tp(4) = XSTT(L,NN,nii,njj)
             end if
             !// Interpolate x direction. Area weight for J-direction.
             vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                    + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
          end do !// do L = 1, LPAR
          measurements(nst)%mass(:,n,NDSTEP) = vp(:)
       end do !// do n = 1, ntracer




       !// Temperature
       do L = 1, LPAR
          tp(1) = T( ii, jj,L)
          tp(2) = T(nii, jj,L)
          tp(3) = T( ii,njj,L)
          tp(4) = T(nii,njj,L)
          !// Interpolate x direction. Area weight for J-direction.
          vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do !// do L = 1, LPAR
       measurements(nst)%temperature(:,NDSTEP)  = vp(:)


       !// Air mass
       do L = 1, LPAR
          tp(1) = AIR( ii, jj,L)
          tp(2) = AIR(nii, jj,L)
          tp(3) = AIR( ii,njj,L)
          tp(4) = AIR(nii,njj,L)
          !// Interpolate x direction. Area weight for J-direction.
          vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do
       measurements(nst)%airmass(:,NDSTEP)      = vp(:)



       !// Height of box bottom
       do L = 1, LPAR+1
          tp(1) = ZOFLE(L, ii, jj)
          tp(2) = ZOFLE(L,nii, jj)
          tp(3) = ZOFLE(L, ii,njj)
          tp(4) = ZOFLE(L,nii,njj)
          !// Interpolate x direction. Area weight for J-direction.
          vpp1(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do
       measurements(nst)%zoflev(:,NDSTEP)       = vpp1(:)


       !// Equivalent latitudes on NTHE theta levels
       do L = 1, NTHE
          tp(1) = theqlat( ii, jj,L)
          tp(2) = theqlat(nii, jj,L)
          tp(3) = theqlat( ii,njj,L)
          tp(4) = theqlat(nii,njj,L)
          !// Interpolate x direction. Area weight for J-direction.
          vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do
       measurements(nst)%eqlat(1:NTHE,NDSTEP)    = vp(1:NTHE)


       !// Potential voritcity units (PVU)
       do L = 1, LPAR
          tp(1) = PVU(L, ii, jj)
          tp(2) = PVU(L,nii, jj)
          tp(3) = PVU(L, ii,njj)
          tp(4) = PVU(L,nii,njj)
          !// Interpolate x direction. Area weight for J-direction.
          vp(L) = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
                 + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       end do
       measurements(nst)%pvu(:,NDSTEP)          = vp(:)


       !// Tropopause height pressure (bottom edge of LMTROP+1)
       tp(1) = etaa( LMTROP( ii, jj)+1 ) + etab( LMTROP( ii, jj)+1 )*p( ii, jj)
       tp(2) = etaa( LMTROP(nii, jj)+1 ) + etab( LMTROP(nii, jj)+1 )*p(nii, jj)
       tp(3) = etaa( LMTROP( ii,njj)+1 ) + etab( LMTROP( ii,njj)+1 )*p( ii,njj)
       tp(4) = etaa( LMTROP(nii,njj)+1 ) + etab( LMTROP(nii,njj)+1 )*p(nii,njj)
       !// Interpolate.
       vs = ((tp(1)*(1._r8 - xfrac) + tp(2)*xfrac) * (1._r8-yfrac) &
           + (tp(3)*(1._r8 - xfrac) + tp(4)*xfrac) * yfrac )
       measurements(nst)%tph_pres(NDSTEP)        = vs


    end do !// do nst = 1, nrofstations


    !// --------------------------------------------------------------------
  end subroutine vprof_stations
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine initialize_stations()
    !// --------------------------------------------------------------------
    !// Initializes the station output based on a input list.
    !// The input list is hard coded in variable filename.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: XDEDG, YDEDG, YDGRD, XDGRD
    use cmn_chem, only: TMASS
    use cmn_oslo, only: XTMASS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Local variables
    integer           :: i,j,l, n, m
    integer           :: iii, jjj, nii, njj, nst, nbb

    integer :: file_nr,file_err
    real(r8)  :: ctmlonb(IPAR+1),ctmlon(IPAR),deltalon,deltalat
    real(r8)  :: tp(4),rtmp
    logical :: file_io, linread
    character(len=80) :: filename

    !// Station data
    integer,dimension(max_stations)  :: ind_i, ind_j     ! grid box
    integer,dimension(max_stations)  :: ind_nbi,ind_nbj  ! closest neighbor boxes
    real(r8),dimension(max_stations)   :: vlats,vlons,valt   ! lat/lon/alt
    real(r8),dimension(max_stations)   :: nbxfrac,nbyfrac ! fraction of neighbor box
    character(len=30),dimension(max_stations) :: locname ! name of stations
    character(len=3),dimension(max_stations)  :: loccode ! station code
    real(r8),dimension(max_stations)   :: axi              ! interpolated grid area

    !// Temporaries for read-in
    character(len=30) :: tname
    character(len=3)  :: tcode
    real(r8)            :: tlat,tlon,talt
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'initialize_stations'
    !// --------------------------------------------------------------------

    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'
    write(6,'(a)') f90file//':'//subr//': Initializing hourly station output'

    !// Check you have the tracers available
    do N = 1, ntracer
       if ( .not. (trsp_idx(stt_components(n)).gt.0 &
                   .or. Xtrsp_idx(stt_components(n)).gt.0) ) then
          write(6,'(a,i5,a)') f90file//':'//subr// &
               ': Tracer number ',stt_components(n), &
               ' is not included in the simulation!'
          stop
       end if
    end do

    !// Get tracer molecular masses
    do N = 1, ntracer
       if ( trsp_idx(stt_components(n)).gt.0) then
          mole_mass(N) = TMASS(trsp_idx(stt_components(n)))
       else
          mole_mass(N) = XTMASS(Xtrsp_idx(stt_components(n)))
       end if
    end do

    !// Check array size
    if (NDIAG .lt. NROPSM*NRMETD) then
       write(6,'(a)') f90file//':'//subr//': NDIAG is too small!'
       write(6,'(a,i3,a,i3)') '    NDIAG:',NDIAG,'NROPSM*NRMETD: ',NROPSM*NRMETD
       stop
    end if


    !// The station list
    filename='stationlist_verticalprofiles.dat'

    !// Find a file number
    file_nr=20
    file_io=.true.
    do while (file_io)
       file_nr=file_nr+1
       inquire(file_nr,opened=file_io)
    end do
    open(file_nr,file=filename,form='formatted',status='old',iostat=file_err)
    if (file_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Could not open file: '//trim(filename)
       write(6,'(a)')   '    Will not simulate station vertical profiles!'
       nrofstations=0
       LINITIALIZED = .false.
       return
    end if

    !// Read the station list
    nst = 0
    linread = .true.
    write(6,'(a)') f90file//':'//subr//': Reading '//trim(filename)
    do while (linread)
       !// Read line
       read(file_nr,'(a3,1x,a30,1x,f8.3,1x,f8.3,1x,f8.1)',iostat=file_err) &
            tcode,tname,tlon,tlat,talt
       if (file_err .eq. 0) then
          !// A line at the end
          if (tcode .eq. '---') then
             exit
          else
             !// If data was read, update the number of stations
             nst = nst + 1    !// Count the stations
             if (nst .le. max_stations) then
                !// Only save for nst.le.max_stations; will stop below if
                !// nst>max_stations.
                loccode(nst) = tcode
                locname(nst) = tname
                vlons(nst) = tlon
                vlats(nst) = tlat
                valt(nst)  = talt
             end if
             write(6,'(i4,1x,a3,1x,a30,1x,f8.3,1x,f8.3,1x,f8.1)') &
                  nst,tcode,tname,tlon,tlat,talt
          end if
       else
          linread = .false.
       end if
    end do

    !// Set number of stations
    nrofstations = nst

    !// Check number of stations
    if (nrofstations .gt. max_stations) then
       write(6,'(a)') f90file//':'//subr//': Too many stations on file!'
       write(6,'(a,i4)') '   max_stations: ',max_stations
       write(6,'(a,i4)') '   nrofstations: ',nrofstations
       stop
    end if


    !// Close file
    close(file_nr)

    write(6,'(a,i4)') '  Number of stations read from file: ',nrofstations
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
    do nst=1,nrofstations
       do j=1,JPAR
          if ( (vlats(nst).ge.YDEDG(j)) .and. (vlats(nst).lt.YDEDG(j+1)) ) then
             !// We have y-index
             jjj=j
             exit
          end if
       end do
       if (vlats(nst).ge.YDEDG(jpar+1)) jjj=jpar !// special for NP
       if (vlats(nst).lt.YDEDG(1)) jjj=1         !// special for SP
       !// Update array
       ind_j(nst) = jjj
    end do

    !// Find grid boxes I-index
    do nst=1,nrofstations
       !// Longitude on degrees in the range <ctmlonb(1),ctmlonb(ipar+1)]
       if (vlons(nst) .gt. ctmlonb(ipar+1)) vlons(nst) = vlons(nst) - 360._r8
       if (vlons(nst) .le. ctmlonb(1)) vlons(nst) = vlons(nst) + 360._r8
       do I=1,IPAR
          if ( (vlons(nst) .gt. ctmlonb(i)) .and. &
               (vlons(nst) .le. ctmlonb(i+1)) ) then
             !// We have x-index
             iii=i
             exit
          end if
       end do
       !// Update array
       ind_i(nst)=iii
    end do

    if (iii.lt.0 .or. jjj.lt.0) then
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': ERROR: III/JJJ is wrong:',iii,jjj
       stop
    end if

    !// Find neighbour boxes: closest i-box
    do nst=1,nrofstations
       !// Double check overlap 0E
       if (ind_i(nst) .eq. 1 .and. vlons(nst) .gt. ctmlonb(ipar) &
            .and. vlons(nst) .le. ctmlonb(ipar+1)) then
          !// This should not be necessary if the above works
          !// Neighbor grid box is just to the west of box 1
          rtmp = vlons(nst) - 360._r8
          write(6,'(a,i7,f12.3)') f90file//':'//subr//': correcting vlons', &
               nst,vlons(nst)
       else
          rtmp = vlons(nst)
       end if
       if (rtmp .ge. ctmlon(ind_i(nst))) then
          !// east of grid box center
          nii = ind_i(nst) + 1
       else
          !// west of grid box center
          nii = ind_i(nst) - 1
       end if
       !// check the "boundaries"
       if (nii .eq. (ipar+1)) nii = 1
       if (nii .eq. 0) nii = ipar
       !// update array
       ind_nbi(nst)=nii
    end do

    !// Find neighbour boxes: closest j-box
    do nst=1,nrofstations
       if (vlats(nst) .gt. YDGRD(ind_j(nst))) then
          njj = ind_j(nst) + 1
       else
          njj = ind_j(nst) - 1
       end if
       !// Check the boundaries
       if (njj .eq. 0) njj=1         !// special for SP
       if (njj .eq. jpar+1) njj=jpar !// special for NP
       !// Update array
       ind_nbj(nst) = njj
    end do


    !// Find fraction of neighbour boxes (can be north or south)
    deltalon = ctmlon(2) - ctmlon(1)
    write(6,'(a)') &
         'Stn   Lon     Grd.cnt iii nbi xfrac  Lat    Grd.cnt jjj nbj yfrac'
    do nst=1,nrofstations
       iii = ind_i(nst)
       jjj = ind_j(nst)
       if (jjj .gt. 1) then
          deltalat = YDGRD(jjj) - YDGRD(jjj-1)
       else
          !// special for SP
          deltalat = YDGRD(2) - YDGRD(1)
       end if
       !// For SP ind_j == ind_nbj, so the fraction really does not matter
       nbyfrac(nst) = abs(vlats(nst) - YDGRD(jjj))/deltalat
       !// Not absolute value for xfrac yet
       nbxfrac(nst) = (vlons(nst) - ctmlon(iii))/deltalon

       !// Double check the fractions
       if (nbyfrac(nst).gt.1. .or. nbxfrac(nst).gt.1.) then
          write(6,'(a)') f90file//':'//subr//': Check nbxfrac/nbyfrac'
          print*,nbxfrac(nst),nbyfrac(nst)
          print*,iii
          print*,ctmlon
          print*,vlons(nst),ctmlon(iii)
          stop
       end if
       if (nbxfrac(nst) .lt. 0._r8) then
          !// Neighbor is to the left
          if ((ind_i(nst) - ind_nbi(nst)) .ne. 1) then
             !// Check if we have 1/ipar
             if (.not. (ind_i(nst).eq.1 .and. ind_nbi(nst).eq.ipar)) then
                write(6,'(a,f12.3)') f90file//':'//subr// &
                     ': nbxfrac 1:',nbxfrac(nst)
                write(6,'(a,i5)') 'ind_i:   ',ind_i(nst)
                write(6,'(a,i5)') 'ind_nbi: ',ind_nbi(nst)
                stop
             end if
          end if
       else
          !// Neighbor is to the right
          if ((ind_nbi(nst) - ind_i(nst)) .ne. 1) then
             !// Check ipar/1
             if (.not. (ind_i(nst).eq.ipar .and. ind_nbi(nst).eq.1)) then
                write(6,'(a,f12.3)') f90file//':'//subr// &
                     ': nbxfrac 2:',nbxfrac(nst)
                write(6,'(a,i5)') 'ind_i:   ',ind_i(nst)
                write(6,'(a,i5)') 'ind_nbi: ',ind_nbi(nst)
                stop
             end if
          end if
       end if


       !// Override with absolute values for xfrac
       nbxfrac(nst) = abs(nbxfrac(nst))

       write(6,'(i4,1x,2f8.3,2i4,f6.3, 2f8.3,2i4,f6.3)') &
            nst,vlons(nst),ctmlon(iii),ind_i(nst),ind_nbi(nst),nbxfrac(nst), &
                vlats(nst),YDGRD(jjj),ind_j(nst),ind_nbj(nst),nbyfrac(nst)
    end do



    !// Find mean grid box area
    do nst=1,nrofstations
       tp(1) = AREAXY(ind_i(nst),  ind_j(nst))
       tp(2) = AREAXY(ind_nbi(nst),ind_j(nst))
       tp(3) = AREAXY(ind_i(nst),  ind_nbj(nst))
       tp(3) = AREAXY(ind_nbi(nst),ind_nbj(nst))
       axi(nst) = &
            ((tp(1)*(1._r8 - nbxfrac(nst)) + tp(2)*nbxfrac(nst)) * (1._r8-nbyfrac(nst)) &
           + (tp(3)*(1._r8 - nbxfrac(nst)) + tp(4)*nbxfrac(nst)) * nbyfrac(nst))
    end do

    !// --------------------------------------------------------------------

    !// Put into measurements array
    do nst=1,nrofstations
       measurements(nst)%locname  = locname(nst) !// Station name
       measurements(nst)%loccode  = loccode(nst) !// Station code
       measurements(nst)%lat      = vlats(nst)   !// Station latitude
       measurements(nst)%lon      = vlons(nst)   !// Station longitude
       measurements(nst)%alt      = valt(nst)    !// Station altitude
       measurements(nst)%ii       = ind_i(nst)   !// CTM i-index
       measurements(nst)%jj       = ind_j(nst)   !// CTM j-index
       measurements(nst)%nb_ii    = ind_nbi(nst) !// CTM closest i-neighbor
       measurements(nst)%nb_jj    = ind_nbj(nst) !// CTM closest j-neighbor
       measurements(nst)%nb_xfrac = nbxfrac(nst) !// Fractional distance x-direction
       measurements(nst)%nb_yfrac = nbyfrac(nst) !// Fractional distance y-direction
       measurements(nst)%areaxy   = axi(nst)
       !// Initialize the hourly data to zero
       measurements(nst)%psfc(:)          = 0._r8
       measurements(nst)%blh(:)           = 0._r8
       measurements(nst)%mass(:,:,:)      = 0._r8
       measurements(nst)%h2o(:,:)         = 0._r8
       measurements(nst)%temperature(:,:) = 0._r8
       measurements(nst)%airmass(:,:)     = 0._r8
       measurements(nst)%zoflev(:,:)      = 0._r8
       measurements(nst)%eqlat(:,:)       = 0._r8
       measurements(nst)%pvu(:,:)         = 0._r8
       measurements(nst)%tph_pres(:)      = 0._r8
    end do

    !// We are initialized
    LINITIALIZED = .true.

    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'

    !// --------------------------------------------------------------------
  end subroutine initialize_stations
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine vprofs_to_file()
    !// --------------------------------------------------------------------
    !// Write vertical profiles to file.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: RESULTDIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    character(LEN=160) :: filename
    character(LEN=8) :: datestamp
    logical :: fnr_ok
    integer :: ifnr, nst, actual_diag
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'vprofs_to_file'
    !// --------------------------------------------------------------------

    !// If no profiles are generated, no files are generated ...
    if (.not.LINITIALIZED) return

    write(datestamp(1:8),'(i4.4,2i2.2)') fileyear,filemonth,filedate
    filename = trim(RESULTDIR)//'hourly_station_vprof_'//datestamp//'.dta'

    !// Find file number to use
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do

    !// We have diagnosed NROPSM*NRMETD steps
    actual_diag = NROPSM*NRMETD

    !// Open file
    open(ifnr,file=filename,form='unformatted')

    !// Write info
    write(ifnr) fileyear,filemonth,filedate,version  !// Time & version
    write(ifnr) lpar,ntracer,nrofstations,nthe, actual_diag,ipar,jpar !// Resolution
    write(ifnr) etaa, etab             !// Sigma coordinates
    write(ifnr) stt_components         !// Components (STT)
    write(ifnr) mole_mass              !// Molecular masses of components
    write(ifnr) pvthe                  !// Theta levels for eqlat
    do nst = 1, nrofstations
       write(ifnr) measurements(nst)%locname
       write(ifnr) measurements(nst)%loccode
       write(ifnr) measurements(nst)%lat
       write(ifnr) measurements(nst)%lon
       write(ifnr) measurements(nst)%alt
       write(ifnr) measurements(nst)%ii
       write(ifnr) measurements(nst)%jj
       write(ifnr) measurements(nst)%nb_ii
       write(ifnr) measurements(nst)%nb_jj
       write(ifnr) measurements(nst)%nb_xfrac
       write(ifnr) measurements(nst)%nb_yfrac
       write(ifnr) real(measurements(nst)%areaxy, r4)
       write(ifnr) real(measurements(nst)%psfc(1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%blh(1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%mass(:,:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%h2o(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%temperature(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%airmass(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%zoflev(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%eqlat(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%pvu(:,1:actual_diag), r4)
       write(ifnr) real(measurements(nst)%tph_pres(1:actual_diag), r4)
    end do
    close(ifnr)

    !// No need to re-initialize; all hours will be overwritten until next time.
    write(6,'(a,i6)') f90file//':'//subr// &
         ': Wrote station profiles: ',nrofstations

    !// --------------------------------------------------------------------
  end subroutine vprofs_to_file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine vprofs_master(JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI)
    !// --------------------------------------------------------------------
    !// Master routine for station profiles output.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JYEAR, JMON, JDATE,NDAY,NMET,NOPS, NDAYI
    !// --------------------------------------------------------------------

    if (NDAY.ne.NDAYI .and. NMET.eq.1 .and. NOPS.eq.1) then

       !// Write data of previous day to file
       call vprofs_to_file()

    else if (NDAY.eq.NDAYI .and. NMET.eq.1 .and. NOPS.eq.1) then

       !// Initialize output for vertical profiles from satellites
       call initialize_stations()

    end if


    !// Process profiles every NOPS
    call vprof_stations(NDAY, NMET, NOPS)


    !// Set date for file name/info
    if (NMET.eq.1 .and. NOPS.eq.1) then
       fileyear = JYEAR
       filemonth= JMON
       filedate = JDATE
    end if

    !// --------------------------------------------------------------------
  end subroutine vprofs_master
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module verticalprofiles_stations2
!//=========================================================================

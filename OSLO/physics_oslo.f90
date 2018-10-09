!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, September 2017
!//=========================================================================
!// Routines for Oslo physics.
!//=========================================================================
module physics_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: physics_oslo
  !// DESCRIPTION: Routines for controlling Oslo physics.
  !// 
  !// Contains:
  !//   subroutine update_physics
  !//   subroutine defineTP
  !//   subroutine defineTP_IJ
  !//   subroutine tp_pvu_ij
  !//   subroutine tp_dtdz_ij
  !//   subroutine tp_e90_ij
  !//   subroutine tp_o3_150
  !//   subroutine get_pvu
  !//   subroutine ijlw2lij
  !//   subroutine theta_pv
  !//   subroutine theta_eqlat
  !//   subroutine IJLfield2ThetaLvs
  !//   subroutine metdata_ij
  !//   subroutine set_blh_ij
  !//   subroutine check_lmtrop
  !//
  !// Amund Sovde Haslerud, September 2017
  !//   Updated for other tropopause definitions.
  !// Ole Amund Sovde, February 2016, September 2008
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  logical,parameter,private :: verbose=.false.

  !// Calculate equivalent latitudes on specific theta levels
  integer,parameter :: NTHE=18
  real(r8), dimension(NTHE),parameter :: &
       pvthe= (/400._r8, 425._r8, 450._r8, 475._r8, 500._r8, 525._r8, &
                550._r8, 575._r8, 600._r8, 650._r8, 700._r8, 750._r8, &
                800._r8, 850._r8, 900._r8, 950._r8,1000._r8,1100._r8/)
  real(r8),dimension(IPAR,JPAR,NTHE) :: pvtheta, theqlat
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'physics_oslo.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public update_physics, get_pvu, ijlw2lij, metdata_ij, &
       NTHE, pvthe, pvtheta, theqlat, IJLfield2ThetaLvs, set_blh_ij
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine update_physics(NDAYI,NDAY,NMET,NOPS,LNEWM)
    !// --------------------------------------------------------------------
    !// Update physics needed in Oslo chemistry.
    !// Also initialises E90-tracer and its tropopause definition.
    !// Done outside parallell loop, but inside NOPS-loop.
    !//
    !// Amund Sovde Haslerud, September 2017
    !//   Updated to handle different LMTROP definitions.
    !// Ole Amund Sovde, November 2014, October 2008
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NRMETD, JMON
    use cmn_sfc, only: LANDUSE_YEAR, LAI_YEAR, ZOI_YEAR
    use grid, only: read_landSurfTypeFrac, read_LAI, read_ZOI
    use utilities_oslo, only: init_e90, tpause_e90
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAYI, NDAY, NMET, NOPS
    logical, intent(in) :: LNEWM
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'update_physics'
    !// --------------------------------------------------------------------
    
    !// Initialize E90 tracer. The routine init_e90 needs LMTROP,
    !// so it must be set from PVU when E90 needs to be
    !// initialised. This has to be done before defineTP, which 
    !// should be able to use the tropopause from E90 tracer.
    if (NDAY.eq.NDAYI .and.NMET.eq.1.and.NOPS.eq.1) then
       call check_lmtrop()
       !// Initialize E90
       call init_e90()
       !// Update E90-tropopause (also done in MP-block loop)
       call tpause_e90()
    end if

    !// Update physical properties each meteorological time step.
    if (NOPS.eq.1) then
       !// Calculate equivalent latitudes on given potential temperature
       !// surfaces. This is done in two steps:
       !// Find PVU on the pre-defined theta levels pvthe
       call theta_pv()
       !// Find equivalent latitude on theta levels pvthe
       call theta_eqlat()

       !// Update tropopause height; needed for some diags and for chemistry,
       !// so we do this outside IJ-blocks.
       call defineTP(NDAY,NMET)

    end if

    if (JMON .eq. 1 .and. LNEWM) then
       !// Set landSurfTypeFrac every 1 Jan, if needed.
       if (LANDUSE_YEAR .eq. 9999) call read_landSurfTypeFrac()
       !// Possible update of LAI
       if (LAI_YEAR .eq. 9999) call read_LAI()
       !// Possible update of ZOI
       if (ZOI_YEAR .eq. 9999) call read_ZOI()

    end if

    !// --------------------------------------------------------------------
  end subroutine update_physics
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine defineTP(NDAY,NMET)
    !// --------------------------------------------------------------------
    !// Define tropopause. Needs to be done for every met-field.
    !// Done outside IJ-block.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Added other options.
    !// Ole Amund Sovde, March 2010
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: ETAA, ETAB, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: P, T, ZOFLE, PVU
    use cmn_chem, only: Ne90, LPAUZTOP
    use cmn_oslo, only: LMTROP, TP_TYPE
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    integer, intent(in):: NDAY, NMET
    !// ------------------------------------------------------------------
    !// Locals
    integer :: MP
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'defineTP'
    !// --------------------------------------------------------------------

    !// Only new tropopause at the beginning of a meteorological time step.
    do MP = 1, MPBLK

       if (TP_TYPE .eq. 1) then
          !// TP_TYPE=1: PVU-based tropopause
          call tp_pvu_ij(LMTROP,P,ETAA,ETAB,T,PVU,ZOFLE,IPAR,JPAR,LPAR, &
               IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else if (TP_TYPE .eq. 2) then
          !// TP_TYPE=2: Tropopause based on lapse rate (-dT/dz)
          call tp_dtdz_ij(LMTROP,P,ETAA,ETAB,T,ZOFLE,IPAR,JPAR,LPAR, &
               IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else if (TP_TYPE .eq. 3) then
          !// TP_TYPE=3: Tropopause based on E90 tracer
          if (Ne90 .le. 0) then
             write(6,'(a,i5)') f90file//':'//subr// &
               ': E90 tracer is not included - cannot use it '// &
               'to calculate tropopause!'
             stop
          end if
          call tp_e90_ij(LMTROP,P,ETAA,ETAB,LPAUZTOP,IPAR,JPAR,LPAR, &
               MPBLK, MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else
          write(6,'(a,i5)') f90file//':'//subr//': '// &
               'NO SUCH TP_TYPE AVAILABLE:',TP_TYPE
          stop
       end if

    end do !// do MP = 1, MPBLK

    write(6,'(a,i5,i3,2(a4,1x,2i3))') '* Tropopause updated',NDAY,NMET, &
         ' SH:',minval(LMTROP(:,1:JPAR/2)), &
                maxval(LMTROP(:,1:JPAR/2)), &
         ' NH:',minval(LMTROP(:,(JPAR/2+1):JPAR)), &
                maxval(LMTROP(:,(JPAR/2+1):JPAR))

    !// --------------------------------------------------------------------
  end subroutine defineTP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine defineTP_IJ(NDAY,NMET,NOPS,NSUB,MP)
    !// --------------------------------------------------------------------
    !// Define tropopause. Needs to be done for every met-field.
    !// Will do this inside IJ-block.
    !// Climatological based tropopause then needs to have a predefined
    !// field.
    !//
    !// Not in use.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: ETAA, ETAB, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: LPAUZTOP
    use cmn_met, only: P, T, ZOFLE, PVU
    use cmn_oslo, only: LMTROP, TP_TYPE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in):: NDAY,NMET,NOPS,NSUB,MP
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'defineTP_IJ'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr// &
         ': Routine is not in use - stopping'
    stop

    !// Only new tropopause at the beginning of a meteorological time step.
    if (NOPS.eq.1 .and. NSUB.eq.1) then

       if (TP_TYPE .eq. 1) then
          !// TP_TYPE=1: PVU-based tropopause
          call tp_pvu_ij(LMTROP,P,ETAA,ETAB,T,PVU,ZOFLE,IPAR,JPAR,LPAR, &
               IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else if (TP_TYPE .eq. 2) then
          !// TP_TYPE=2: Tropopause based on lapse rate
          call tp_dtdz_ij(LMTROP,P,ETAA,ETAB,T,ZOFLE,IPAR,JPAR,LPAR, &
               IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else if (TP_TYPE .eq. 3) then
          !// TP_TYPE=3: Tropopause based on E90 tracer
          call tp_e90_ij(LMTROP,P,ETAA,ETAB,LPAUZTOP,IPAR,JPAR,LPAR, &
               MPBLK, MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)

       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': NO SUCH TP_TYPE AVAILABLE:',TP_TYPE
          stop
       end if

    end if !// if (NOPS.eq.1 .and. NSUB.eq.1) then

    !// --------------------------------------------------------------------
  end subroutine defineTP_IJ
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tp_pvu_ij(LMTROP,PSFC,ETAA,ETAB,T,PVU,ZOFL,IPAR,JPAR, &
       LPAR,IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)
    !// --------------------------------------------------------------------
    !// Get the tropopause height from PV-data. PV is given by the almost
    !// all of the IFS meteorological data, and when it is not given, it is
    !// calculated from wind data.
    !//
    !// The upper limit of the tropopause is potential temperature of 380K
    !// at the grid box upper edge. A lower edge is set to 5km.
    !//
    !// Ole Amund Sovde, Noveber 2015, September 2009
    !// --------------------------------------------------------------------
    use cmn_parameters, only: R_AIR, cp_air
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLK,MP
    integer,dimension(MPBLK), intent(in) :: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    real(r8), dimension(IPAR,JPAR,LPAR), intent(in) :: T
    !// PVU is reversed
    real(r8), dimension(LPAR,IPAR,JPAR), intent(in) :: PVU
    real(r8), dimension(LPAR+1,IPAR,JPAR), intent(in) :: ZOFL
    real(r8), dimension(IPAR,JPAR), intent(in) :: PSFC
    real(r8), dimension(LPAR+1), intent(in)    :: ETAA, ETAB
    !// Output
    integer, intent(out) :: LMTROP(IPAR,JPAR)

    !// Parameters
    real(r8), parameter :: pvlimit   = 2.5_r8 !// pvu-limit for tropopause [PVU]
    real(r8), parameter :: tplim_min = 5.e3_r8  !// min tropopause height [m]
    real(r8), parameter :: rocp = -R_AIR / cp_air !// -R/cp
    logical, parameter :: verbose=.false.

    !// Local variables
    integer :: I, J, L, II, JJ, Lminpvu
    real(r8) :: &
         pvminval, &        !// min PVU in the column
         pvmaxval, &        !// min PVU in the column
         theta_hi, &        !// potential temperature [K] upper box edge
         pres              !// upper edge pressure
    integer :: L380     !// Uppermost level for theta_hi <= 380K
    integer :: LMT      !// LMTROP for (I,J)
    integer :: LMIN,LMAX !// Lowest/highest found LMT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tp_pvu_ij'
    !// --------------------------------------------------------------------

    !// Calculations are done inside IJ-block!
    !//   I,J  : the actual (global) indices of longitude/latitude.
    !//   II,JJ: the indices of the IJ-block.
 
    !// Initialize tropopause height
    do J = MPBLKJB(MP), MPBLKJE(MP)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          LMTROP(I,J) = 0 !// Initialize tropopause height
       end do
    end do

    LMIN = 1000 !// Keep track of the lowest LMT in the IJ-block
    LMAX = 0    !// Keep track of the highest LMT in the IJ-block

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Will calculate tropopause level LMT for this (I,J)
          LMT = 0

          !// Find the level for 380K, based on potential temperature at top
          !// of boxes: THETA = T * (p0/p)^(r/cp) = T * (p/p0)^(-r/cp)
          L380 = 0
          !// Loop from above due to THETA structure.
          !// Start at LPAR-1 since T is in grid center
          do L = LPAR-1, 2, -1
             !// Pressure on upper edge of grid box [hPa]
             pres = ETAA(L+1) + ETAB(L+1) * PSFC(I,J)
             !// T is box center value
             theta_hi = (T(I,J,L) + T(I,J,L+1))*0.5_r8 *(pres*1.e-3_r8)**rocp
             if (theta_hi .le. 380._r8) then
                !// Upper boundary is less than 380K; set this as uppermost
                !// troposphere.
                L380 = L
                exit
             end if
          end do


          !// Check absolute PVU value from top and downwards, and where
          !// the PVU is less than pvlimit, we have the tropopause level.

          !// Find minimum PVU in column (only down to level 2)
          pvminval = 1000._r8
          do L = L380, 2, -1
             if (pvminval .ge. abs(PVU(L,I,J))) then
                pvminval = abs(PVU(L,I,J))
                Lminpvu = L
             end if
          end do
          if (pvminval .le. pvlimit) then
             !// Somewhere PVU is smaller than the limit. This leaves two
             !// options;
             !// 1. PVU > pvlimit above, so we can find PVU = pvlimit
             !// 2. PVU is always < pvlimit (low latitudes). Need to use 380K
             !// Start at LPAR-1 since T is in grid center
             do L = L380, 2, -1
                !// Checking regime ( PVU < pvlimit & theta_hi(L) <= 380K)
                if (abs(PVU(L,I,J)) .le. pvlimit .and. &
                     L .le. L380) then 
                   !// We got the height of the uppermost troposphere level
                   LMT = L
                   exit
                else
                   !// Keep looping downwards; do nothing unless we reach
                   !// minimum allowed TP height.
                   if (abs(PVU(L,I,J)) .gt. pvlimit .and. &
                      ZOFL(L,I,J) .le. tplim_min) then
                      !// Bottom of box is below tplim_min.
                      !// Let tplim_min be lowest value, set TP to this layer
                      LMT = L
                      exit
                   end if
                end if
             end do
             if (LMT .le. 1) then
                print*,'LMTROP<1',I,J,II,JJ,MP,pvminval,minval(abs(PVU(:,I,J)))
                do L = LPAR, 1, -1
                   print*,L,I,J,abs(PVU(L,I,J)),L380
                end do
                stop
             end if
          else
             !// PVU larger than pvminval, check lowest PVU value: if it is
             !// below tplim_min, then we choose TP to be at tplim_min.
             !// If the lowest PVU value is above tplim_min, we use the height
             !// where PVU is lowest.
             !// Do not need the 380K test here, since the TP should be low in
             !// this case, but we still loop from L380 since that is the max
             !// allowed TP height.
             write(6,'(a,f7.3,a,f7.3)')'*** PVU larger than ',pvlimit, &
                  '; trying to use',pvminval
             do L = L380, 2, -1
                if (ZOFL(L,I,J) .le. tplim_min) then
                   !// Lower than minimum height, setting the TP level at
                   !// tplim_min.
                   write(6,'(a,2i3,f7.3,f8.2,i3,i4,i3)') &
                        '    TP1X: TP < tplim_min',L,Lminpvu, &
                        abs(PVU(L,I,J)), ZOFL(L,i,j),MP,I,J
                   LMT = L
                   exit
                else
                   !// Above minimum height. Check for minimum value, or go to
                   !// level below.
                   if (abs(PVU(L,I,J)) .eq. pvminval) then
                      !// Have minimum value. PVU-minimum is between min
                      !// (tplim_min) and max height (L380). I doubt that PVU
                      !// will have minimum above 380K, but stick to that
                      !// maximum value.
                      LMT = L
                      write(6,'(a,2i3,f7.3,f8.2,i3,i4,i3)') &
                           '    TP2X:',L,Lminpvu, &
                           abs(PVU(L,I,J)), ZOFL(L,i,j),MP,I,J
                      exit
                   else
                      !// Keep looping downwards. The first test on ZOFL will
                      !// catch us if we fall below tplim_min.
                   end if
                end if

             end do !// do L = L380, 2, -1

          end if !// if (pvminval .le. pvlimit) then

          !// Set the calculated LMTROP
          LMTROP(I,J) = LMT
          !// Keep track of the lowest/highest LMT in this IJ-block
          if (LMT .lt. LMIN) LMIN = LMT
          if (LMT .gt. LMAX) LMAX = LMT

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    if (LMIN .eq. 0) then
       write(6,'(a,2i5)') f90file//':'//subr//': LMTROP problem? '// &
            'Min. level zero. MP/LMAX: ',MP,LMAX
       stop
    end if

    if (verbose) write(6,'(A,I2,1x,I2,A3,I2)') &
         f90file//':'//subr//': Tropopause updated from PVU: MIN/MAX: ', &
         MP, LMIN, LMAX

    !// --------------------------------------------------------------------
  end subroutine tp_pvu_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tp_dtdz_ij(LMTROP,PSFC,ETAA,ETAB,T,ZOFL,IPAR,JPAR, &
       LPAR,IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)
    !// --------------------------------------------------------------------
    !// Get the tropopause height from temperature data.
    !// 
    !// Tropopause is where -dT/dz falls below 2K/km.
    !//
    !// Amund Sovde Haslerud, August 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: R_AIR, cp_air
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLK,MP
    integer,dimension(MPBLK), intent(in) :: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    real(r8), dimension(IPAR,JPAR,LPAR), intent(in) :: T
    real(r8), dimension(LPAR+1,IPAR,JPAR), intent(in) :: ZOFL
    real(r8), dimension(IPAR,JPAR), intent(in) :: PSFC
    real(r8), dimension(LPAR+1), intent(in)    :: ETAA, ETAB
    !// Output
    integer, intent(out) :: LMTROP(IPAR,JPAR)

    !// Parameters
    real(r8), parameter :: gama_tp = 2._r8 !// 2K/km
    real(r8), parameter :: tpp_limit = 50._r8  !// min TP pressure [hPa]
    real(r8), parameter :: alt_limit = 5.e3_r8  !// min TP height [m]
    logical, parameter :: verbose=.false.

    !// Local variables
    integer :: I, J, L, II, JJ
    real(r8),dimension(LPAR) :: vtem, vpres
    real(r8) :: &
         dT, dz, lapse, minpres
    integer :: LMT      !// LMTROP for (I,J)
    integer :: LMT0,LMT1
    integer :: LMIN,LMAX !// Lowest/highest found LMT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tp_dtdz_ij'
    !// --------------------------------------------------------------------

    !// Initialize tropopause height
    do J = MPBLKJB(MP), MPBLKJE(MP)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          LMTROP(I,J) = 0 !// Initialize tropopause height
       end do
    end do

    LMIN = 1000 !// Keep track of the lowest LMT in the IJ-block
    LMAX = 0    !// Keep track of the highest LMT in the IJ-block
    minpres = 10000._r8 !// minimum TP pressure in IJ-block

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Calculate tropopause level LMT for this (I,J)
          LMT = 0

          !// Vertical temperature
          vtem(:) = T(I,J,:)

          !// Pressure grid box tops
          do L = 1, LPAR
             vpres(L) = etaa(L+1) + psfc(i,j) * etab(L+1)
          end do

          !// Find lowest allowed tropopause level LMT0
          LMT0 = 1
          do L = 1, LPAR-1
             !// Check altitude of grid box top
             !// (i.e. where top is higher than alt_limit)
             if (ZOFL(L+1,I,J) - ZOFL(1,I,J) .ge. alt_limit) then
                LMT0 = L
                exit
             end if
          end do !// do L = 1, LPAR-1

          !// Initialise LMT
          LMT = LMT0
          do L = LMT0, LPAR-1

             !// Calculate lapse rate
             !// From L to L+1, but should check from L-1 to L also.

             !// delta T [K]
             dT = vtem(L+1) - vtem(L)
             !// delta Z [km]
             dz = 1.e-3_r8 * (ZOFL(L+1,I,J) - ZOFL(L,I,J))
             !// Lapse rate is -dT/dz
             lapse = -dT / dz

             if (lapse .le. 2._r8) then
                LMT = L
                exit !// exit loop
             end if

          end do !// do L = LMT0, LPAR-1


          !// Do not allow TP above some limit
          !// Make sure to also check if 
          if (vpres(LMT) .le. tpp_limit) then
             LMT1 = LMT
             do L = LMT1, LMT0, -1
                !// Use the first level where pressure is > tpp_limit
                if (vpres(L) .ge. tpp_limit) then
                   LMT = L
                   exit
                end if
             end do
             write(6,'(a,5i5,2es10.2)') f90file//':'//subr// &
                  ': TPP too low - adjusted!',i,j,LMT1,LMT,vpres(LMT)
          end if


          !// Set the calculated LMTROP
          LMTROP(I,J) = LMT
          !// Keep track of the lowest/highest LMT in this IJ-block
          if (LMT .lt. LMIN) LMIN = LMT
          if (LMT .gt. LMAX) LMAX = LMT
          !// Keep track of the lowest TP pressure in this IJ-block
          if (vpres(LMT) .le. minpres) minpres = vpres(LMT)


       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)



    if (verbose) write(6,'(A,I2,1x,I2,A3,I2,f6.1)') &
         f90file//':'//subr//' Tropopause updated: MIN/MAX: ', &
         MP, LMIN, LMAX, minpres

    !// --------------------------------------------------------------------
  end subroutine tp_dtdz_ij
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine tp_e90_ij(LMTROP,PSFC,ETAA,ETAB,LPAUZTOP,IPAR,JPAR,LPAR, &
       MPBLK, MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)
    !// --------------------------------------------------------------------
    !// Get the tropopause height from LPAUZTOP.
    !//
    !// Amund Sovde Haslerud, August 2017
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR,JPAR,LPAR,MPBLK,MP
    integer,dimension(MPBLK), intent(in) :: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    real(r8), dimension(IPAR,JPAR), intent(in) :: PSFC
    real(r8), dimension(LPAR+1), intent(in)    :: ETAA, ETAB
    integer, dimension(IPAR,JPAR), intent(in)  :: LPAUZTOP
    !// Output
    integer, intent(out) :: LMTROP(IPAR,JPAR)

    !// Local variables
    integer :: I, J, LMT
    real(r8) :: ptop, minpres
    integer :: LMIN, LMAX !// Lowest/highest found LMT
    logical, parameter :: verbose=.false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tp_e90_ij'
    !// --------------------------------------------------------------------

    LMIN = 1000 !// Keep track of the lowest LMT in the IJ-block
    LMAX = 0    !// Keep track of the highest LMT in the IJ-block
    minpres = 10000._r8 !// minimum TP pressure in IJ-block

    do J = MPBLKJB(MP), MPBLKJE(MP)
       do I = MPBLKIB(MP), MPBLKIE(MP)

          LMT = LPAUZTOP(I,J)

          LMTROP(I,J) = LMT

          !// Pressure grid box tops
          ptop = etaa(LMT+1) + psfc(I,J) * etab(LMT+1)

          !// Keep track of the lowest/highest LMT in this IJ-block
          if (LMT .lt. LMIN) LMIN = LMT
          if (LMT .gt. LMAX) LMAX = LMT
          !// Keep track of the lowest TP pressure in this IJ-block
          if (ptop .le. minpres) minpres = ptop

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    if (verbose) write(6,'(A,I2,1x,I2,A3,I2,f6.1)') &
         f90file//':'//subr//' Tropopause updated: MIN/MAX: ', &
         MP, LMIN, LMAX, minpres

    !// --------------------------------------------------------------------
  end subroutine tp_e90_ij
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine tp_o3_150(LMTROP150,IPAR,JPAR,LPAR,NPAR,MTCO3,STT,AIR)
    !// --------------------------------------------------------------------
    !// Find tropopause based on where O3<=150ppbv.
    !// Currently not for MP-blocks.
    !//
    !// Ole Amund Sovde, February 2011
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR,JPAR,LPAR,NPAR,MTCO3
    real(r8), intent(in) :: STT(IPAR,JPAR,LPAR,NPAR),AIR(IPAR,JPAR,LPAR)
    !// Output
    integer, intent(out) :: LMTROP150(IPAR,JPAR)
    !// Local variables
    integer :: I, J, L
    !// O3-tp limit as mass mixing ratio
    real(r8),parameter :: o3lim = 150.e-9_r8 * 48._r8 / M_AIR
    !// --------------------------------------------------------------------

    if (MTCO3 .le. 0) then
       LMTROP150(:,:) = 1
       return
    end if

    !// Initialize
    LMTROP150(:,:) = 10
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 10, LPAR
             if (STT(I,J,MTCO3,L) .gt. AIR(I,J,L)*o3lim) then
                LMTROP150(I,J) = L-1
                exit
             end if
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine tp_o3_150
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_pvu()
    !// --------------------------------------------------------------------
    !// Calculate PV from met-fields. This is only done for IPAR,JPAR,LPAR
    !// resolution.
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: ETAA, ETAB, YGRD
    use cmn_met, only: P, T, UMS, VMS, PVU
    use cmn_parameters, only: G0, CPI, A0, R_AIR, CP_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: I,J,L
    real(r8) :: &
         theta_lo,theta_hi, & ! Pot. temp. boundaries [K]
         deltap, &                ! difference in pressure [hPa]
         deltav,deltau, &         ! difference in velocity [m/s]
         deltax,deltay, &         ! difference in extention [m]
         cor                      ! coriolis parameter
    !// --------------------------------------------------------------------

    !// Initialize
    PVU(:,:,:) = 0._r8

    do L = 2, LPAR-1 !// due to theta_lo and theta_hi
      do J = 1, JPAR
        do I = 1, IPAR

          ! potential temperature at lower boundary of layer L [K]
          theta_lo = (T(I,J,L-1) + T(I,J,L)) * 0.5_r8 &
               * ( P(I,J) / ( ETAA(L) + ETAB(L) * P(I,J) ) )**(r_air/cp_air)
          ! potential temperature at upper boundary of layer L [K]
          theta_hi = (T(I,J,L) + T(I,J,L+1)) * 0.5_r8 &
               * ( P(I,J) / ( ETAA(L+1) + ETAB(L+1) * P(I,J) ) )**(r_air/cp_air)
          ! pressure difference across layer L [hPa]
          !dp = ETAA(L)-ETAA(L+1) + (ETAB(L)-ETAB(L+1)) * P(I,J)
          deltap = (ETAA(L+1) - ETAA(L) + (ETAB(L+1) - ETAB(L)) * P(I,J))
          ! Coriolis parameter at J [1/s]
          cor = 2._r8 * 7.292e-5_r8 * sin(YGRD(J))
          ! v diff between I+1 and I-1 [m/s]
          if (I .eq. 1) then
             deltav = VMS(L,I+1,J) - VMS(L,IPAR,J)
          else if (I .eq. IPAR) then
             deltav = VMS(L,1,J) - VMS(L,I-1,J)
          else
             deltav = VMS(L,I+1,J) - VMS(L,I-1,J)
          end if
          ! dx [m]
          deltax = 2._r8 * 2._r8 * CPI * A0 * cos(YGRD(J)) / real(IPAR, r8)
          ! u diff between I+1 and I-1 [m/s]
          if (J .eq. 1) then
             deltau = UMS(L,I,J+1) - UMS(L,I,J)
          else if (J .eq. JPAR) then
             deltau = UMS(L,I,J) - UMS(L,I,J-1)
          else
             deltau = UMS(L,I,J+1) - UMS(L,I,J-1)
          end if
          ! dy [m]
          if (J .eq. 1) then
             deltay = A0 * (YGRD(J+1) - (YGRD(J)))
          else if (J .eq. JPAR) then
             deltay = A0 * (YGRD(J) - (YGRD(J-1)))
          else   
             deltay = A0 * (YGRD(J+1) - (YGRD(J-1)))
          end if


          !// Calculate PV (1.d-2 to get dp in Pa).
          !// Multiply by 1.d6 to get PVU.
          PVU(L,I,J) = -(deltav/deltax - deltau/deltay + cor) * g0 &
               * (theta_hi - theta_lo)/deltap*1.e-2_r8 * 1.e6_r8

        end do ! do I = 1, IPAR
      end do ! do J = 1, JPAR
    end do ! do L = 2, LPAR-1

    !// Set PV at surface = first layer above
    PVU(1,:,:) = PVU(2,:,:)
    !// And top = the one below
    PVU(LPAR,:,:) = PVU(LPAR-1,:,:)

    !// --------------------------------------------------------------------
  end subroutine get_pvu
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ijlw2lij(FLDW,LMAP,XLMMAP,IM,JM,LM,IPARW, &
       JPARW,LPARW,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,JDGRD,FLD)
    !// --------------------------------------------------------------------
    !// Put FLDW(IPARW,JPARW,LPARW) on a degraded grid, and convert to
    !// FLD(LPAR,IPAR,JPAR).
    !//
    !// Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use regridding, only: TRUNG8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IM,JM,LM,IPARW,JPARW,LPARW,IDGRD,JDGRD
    real(r8), dimension(IPARW,JPARW,LPARW), intent(in)  :: FLDW
    real(r8), intent(in) :: XLMMAP(LPARW+1)
    integer,intent(in)  :: LMAP(LM+1)
    integer,intent(in)  :: IMAP(IDGRD,IM), JMAP(JDGRD,JM)
    real(r8), intent(in)  :: ZDEGI(IDGRD,IM), ZDEGJ(JDGRD,JM)
    !// Output
    real(r8),dimension(LM,IM,JM), intent(out)  :: FLD
    !// Locals
    integer  :: I,J,L,LL
    real(r8) :: VDEG(IPARW,JPARW,LM) !// Vertically degraded
    real(r8) :: stdArray(IM,JM,LM)
    !// --------------------------------------------------------------------


    !// Initialize
    VDEG(:,:,:)     = 0._r8 !// Degraded vertically
    stdArray(:,:,:) = 0._r8 !// Degraded vertically AND horizontally
    FLD(:,:,:)      = 0._r8 !// Output field

    !// Collapsing layers
    do L = 1,LM !// Remember that LM = LPAR
       do LL = LMAP(L),LMAP(L+1)-1
          do J = 1,JPARW
             do I = 1,IPARW
                VDEG(I,J,L) = VDEG(I,J,L) + FLDW(I,J,LL)*XLMMAP(LL)
             end do
          end do
       end do
    end do

    !// Put field VDEG onto stdArray, degrading or not
    call TRUNG8(VDEG,stdArray,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD, &
         JDGRD,IPARW,JPARW,IM,JM,LM,1)

    !// Reverse to LIJ
    do L = 1, LM
       do J = 1, JM
          do I = 1, IM
             FLD(L,I,J) = stdArray(I,J,L)
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine ijlw2lij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine theta_pv()
    !// --------------------------------------------------------------------
    !// Interpolate PV to pre-defined theta levels.
    !//
    !// Must be called outside parallell loop, since it uses openmp.
    !//
    !// Ole Amund Sovde, December 2009 - February 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR
    use cmn_ctm, only: ETAA, ETAB
    use cmn_met, only: T, P, PVU
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    real(r8)  :: theta(LPAR), pres
    integer :: the_low(NTHE), the_hig(NTHE)
    integer :: I,J,L, TH
    !// --------------------------------------------------------------------

    !// Loop through columns; in parallel
!$omp parallel &
!$omp  private(L,I,J,TH,theta,pres,the_low,the_hig) &
!$omp  shared(ETAA,ETAB,P,PVU,T,pvtheta) &
!$omp  default(NONE)
!$omp do
    do j = 1, JPAR
       do i = 1, IPAR
          !// Find potential temperature
          do l=1,LPAR
             pres = 0.5_r8 * ( ETAA(L) + ETAA(L+1) &
                             + P(I,J)*(ETAB(L) + ETAB(L+1)) )
             theta(L) = T(I,J,L) * (pres * 1.e-3_r8)**(-0.2859_r8)
          end do

          !// Find model layers closest to the specified theta levels
          the_low(:) = 1
          the_hig(:) = 2
          do TH = NTHE, 1, -1
             do L = LPAR-1, 1, -1
                if (theta(L) .lt. pvthe(TH)) then
                   the_low(TH) = L   !// CTM level below pvthe(TH)
                   the_hig(TH) = L+1 !// CTM level above pvthe(TH)
                   exit !// exit L-loop
                end if
             end do
          end do !// do TH = NTHE, 1, -1

          !// Interpolate pvu between the_low:the_hig
          !// p_t = pv_l * (t_h - t) / (t_h-t_l)  + pv_l * (t-t_l) / (t_h-t_l)
          !//     + pv_h * (t - t_l) / (t_h-t_l)  - pv_l * (t-t_l) / (t_h-t_l)
          !//     = pv_l + (pv_h  - pv_l) * (t - t_l)/(t_h - t_l)
          do TH = 1, NTHE
             pvtheta(I,J,TH) = pvu(the_low(TH),I,J) &
                  + (pvu(the_hig(TH),I,J) - pvu(the_low(TH),I,J)) &
                    / (theta(the_hig(TH)) - theta(the_low(TH))) &
                    * (pvthe(TH) - theta(the_low(TH)))
          end do !// do TH = 1, NTHE

       end do !// do i = 1, IPAR
    end do !// do j = 1, JPAR
!$omp end do
!$omp end parallel
    !// --------------------------------------------------------------------
  end subroutine theta_pv
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine theta_eqlat()
    !// --------------------------------------------------------------------
    !// Find equivalent latitude at pre-defined theta levels, from
    !// interpolated PV on the pre-defined theta levels.
    !//
    !// Must be called outside parallell loop, since it uses openmp.
    !//
    !// Ole Amund Sovde, December 2009 - February 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: AREAXY
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Number of PV-steps to calculate per hemisphere
    integer,parameter :: antsteps=128
    real(r8)          :: zstep,pvbins
    integer :: TH, I, J, AI, arrayindex, JP2
    integer :: eqlat_indices(IPAR,JPAR)  ! the indices for finding eqlat
    real(r8) :: &
         maxpv,minpv, &
         pvstep,totarea,curpv,curarea, &
         pvarea(antsteps), &
         equivlat(antsteps)
    !// --------------------------------------------------------------------
    !// Initialize eqlat
    theqlat(:,:,:) = 0._r8
    JP2 = JPAR/2
    pvbins = real(antsteps-1, r8)
    zstep = 1._r8/pvbins

!$omp parallel &
!$omp  private(TH,I,J,AI,maxPV,minPV,curPV,PVstep,totarea, &
!$omp          PVarea,curArea,eqlat_indices,equivLat,arrayIndex) &
!$omp  shared(JP2,pvbins,zstep,theqlat,pvtheta,AREAXY) &
!$omp  default(NONE)
!$omp do schedule(static,1)
    !// Loop over theta levels
    do TH = 1, NTHE

       !// Northern Hemisphere first
       !// ------------------------------------------------------------------
       maxPV = maxval(pvtheta(:,(JP2+1):JPAR,TH))
       minPV = minval(pvtheta(:,(JP2+1):JPAR,TH))

       PVstep  = (maxPV-minPV) * zstep !// the PV step
       totArea = 0._r8                  !// total area


       PVarea(:) = 0._r8
       eqlat_indices(:,:) = -1
       ! need the ctm2 areas areay(:) to calculate the PV area

       do J = JP2+1, JPAR
          do I = 1, IPAR

             !// Current PVU
             curPV = pvtheta(I,J,TH)

             if (curPV .ge. minPV) then
                !// NH index of PVU
                arrayIndex = int(pvbins / (maxPV-minPV)*(curPV-minPV)) + 1
                if (arrayindex .lt. 1) then
                   print*,'AI WRONG NH',arrayindex,maxpv,minpv,curpv, &
                        (maxPV-minPV)*(curPV-minPV)
                   stop
                end if
                curArea = AREAXY(i,j)          !// grid box area
                totArea = totArea + curArea    !// summing up aera
                !// find area enclosing each PVU line defined by
                !// minPV,maxPV and PVstep
                do AI = arrayIndex, 1, -1
                   PVarea(AI) = PVarea(AI) + curArea
                end do
                !// The index corresponding to the eqlat
                eqlat_indices(I,J) = arrayIndex
             end if
          end do
       end do

       !// Normalize to totArea
       PVarea(:) = PVarea(:)/totArea

       !// Finally we retrieve the equivalent latitude
       do AI = 1, antSteps
          equivLat(AI) = 57.29577951_r8 * asin(1._r8 - PVarea(AI))
       end do
       do J = JP2+1, JPAR
          do I = 1, IPAR
             if (eqlat_indices(I,J) .gt. 0) then
                theqlat(I,J,TH) = equivLat(eqlat_indices(I,J))
             else
                theqlat(I,J,TH) = 0._r8
             end if
          end do
       end do



       !// Southern Hemisphere
       !// ------------------------------------------------------------------
       !// Note that for SH the minimum PV is at the South Pole
       maxPV = maxval(pvtheta(:,1:JP2,TH))
       minPV = minval(pvtheta(:,1:JP2,TH))

       PVstep  = (maxPV-minPV) * zstep !// the PV step
       totArea = 0._r8                  !// total area

       PVarea(:) = 0._r8
       eqlat_indices(:,:) = -1

       do J = 1, JP2
          do I = 1, IPAR
             !// Current PVU
             curPV = pvtheta(I,J,TH)
             if (curPV .ge. minPV) then
                !// SH index of PVU
                arrayIndex = int(pvbins / (maxPV-minPV)*(maxPV-curPV)) + 1
                if (arrayindex .lt. 1) then
                   print*,'AI WRONG SH',arrayindex,maxpv,minpv,curpv, &
                        (maxPV-minPV),(curPV-minPV),(curPV-minPV)/(maxPV-minPV)
                   stop
                end if
                curArea = AREAXY(I,J)          !// grid box area
                totArea = totArea + curArea    !// summing up aera
                !// find area enclosing each PVU line defined by
                !// minPV,maxPV and PVstep
                do AI=arrayIndex,1,-1
                   PVarea(AI) = PVarea(AI) + curArea
                end do
                !// the index corresponding to the eqlat
                eqlat_indices(I,J) = arrayIndex
             end if
          end do
       end do


       !// Normalize to totArea
       PVarea(:) = PVarea(:)/totArea

       !// finally we retrieve the equivalent latitude
       do AI = 1, antSteps
          equivLat(AI)=-57.29577951_r8 * asin(1._r8 - PVarea(AI))
       end do
       do J = 1, JP2
          do I = 1, IPAR
             if (eqlat_indices(I,J) .gt. 0) then
                theqlat(I,J,TH) = equivLat(eqlat_indices(I,J))
             else
                theqlat(I,J,TH) = 0._r8
             end if
          end do
       end do

    end do !// do TH = 1, NTHE
!$omp end do
!$omp end parallel

    !// --------------------------------------------------------------------
  end subroutine theta_eqlat
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine IJLfield2ThetaLvs(indata, outdata, JSTART, JEND)
    !// --------------------------------------------------------------------
    !// Interpolate IJL-field to pre-defined theta levels.
    !//
    !// Must be called outside parallell loop, since it uses openmp.
    !//
    !// Ole Amund Sovde, February 2016
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR
    use cmn_ctm, only: ETAA, ETAB
    use cmn_met, only: T, P, PVU
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JSTART, JEND
    real(r8), intent(in) :: indata(IPAR,JPAR,LPAR)
    !// Output
    real(r8), intent(out) :: outdata(IPAR,JPAR,NTHE)
    !// Locals
    real(r8)  :: theta(LPAR), pres
    integer :: the_low(NTHE), the_hig(NTHE)
    integer :: I,J,L, TH
    !// --------------------------------------------------------------------

    !// Initialise
    outdata(:,:,:) = 0._r8

    !// Loop through columns; in parallel
!$omp parallel &
!$omp  private(L,I,J,TH,theta,pres,the_low,the_hig) &
!$omp  shared(JSTART, JEND, ETAA,ETAB,P,T,indata,outdata) &
!$omp  default(NONE)
!$omp do
    do J = JSTART, JEND
       do I = 1, IPAR
          !// Find potential temperature
          do L = 1, LPAR
             pres = 0.5_r8 * ( ETAA(L) + ETAA(L+1) &
                             + P(I,J)*(ETAB(L) + ETAB(L+1)) )
             theta(L) = T(I,J,L) * (pres * 1.e-3_r8)**(-0.2859_r8)
          end do

          !// Find model layers closest to the specified theta levels
          the_low(:) = 1
          the_hig(:) = 2
          do TH = NTHE, 1, -1
             do L = LPAR-1, 1, -1
                if (theta(L) .lt. pvthe(TH)) then
                   the_low(TH) = L   !// CTM level below pvthe(TH)
                   the_hig(TH) = L+1 !// CTM level above pvthe(TH)
                   exit !// exit L-loop
                end if
             end do
          end do !// do TH = NTHE, 1, -1

          !// Interpolate tracer F between the_low:the_hig
          !// I_h = indata for upper theta, I_l = indata for lower theta
          !// I_t = I_l * (t_h - t) / (t_h-t_l)  + I_l * (t-t_l) / (t_h-t_l)
          !//     + I_h * (t - t_l) / (t_h-t_l)  - I_l * (t-t_l) / (t_h-t_l)
          !//     = I_l + (I_h  - I_l) * (t - t_l)/(t_h - t_l)
          do TH = 1, NTHE
             outdata(I,J,TH) = indata(I,J,the_low(TH)) &
                  + (indata(I,J,the_hig(TH)) - indata(I,J,the_low(TH))) &
                    / (theta(the_hig(TH)) - theta(the_low(TH))) &
                    * (pvthe(TH) - theta(the_low(TH)))
          end do !// do TH = 1, NTHE

       end do !// do i = 1, IPAR
    end do !// do j = JSTART, JEND
!$omp end do
!$omp end parallel
    !// --------------------------------------------------------------------
  end subroutine IJLfield2ThetaLvs
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine metdata_ij(AIRB, BTEM, MP)
    !// --------------------------------------------------------------------
    !// Sets private array of air molecular density AIRMOLEC_IJ and
    !// box volume DV_IJ. They are used in converting to concentration.
    !//
    !// Ole Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, AREAXY
    use cmn_met, only: ZOFLE
    use cmn_parameters, only: M_AIR, AVOGNR, R_AIR, G0
    use cmn_oslo, only: DV_IJ, AIRMOLEC_IJ
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), dimension(LPAR,IDBLK,JDBLK),intent(in) :: AIRB, BTEM
    !// Locals
    !// Pressure variables, thickness
    real(r8) :: P1, P2, DZ

    !// For looping
    integer :: I, J, L, II, JJ

    !// Some helping constants
    real(r8),parameter ::  &
         AIR2AM = AVOGNR * 1.e-3_r8 / M_AIR
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Generating vertical arrays
          do L = 1, LPAR

             !// It is now possible to use the global array ZOFLE
             DZ = ZOFLE(L+1,I,J) - ZOFLE(L,I,J)

             !// Box volume [m3]
             DV_IJ(L,II,JJ,MP) = AREAXY(I,J) * DZ

             !// Calculate air density [molec/cm^3]
             AIRMOLEC_IJ(L,II,JJ,MP) = &
                  AIRB(L,II,JJ) / DV_IJ(L,II,JJ,MP) * AIR2AM

          end do !// do L = 1, LPAR

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine metdata_ij
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine set_blh_ij(NMET, NOPS, NSUB, CCYC, dtchem2, &
       dtmet, nmetTimeIntegrated, MP)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR
    use cmn_ctm, only: NROPSM, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: BLH, BLH_CUR, BLH_NEXT, LBLH, ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NMET, NOPS, NSUB, CCYC, MP
    real(r8), intent(in) :: dtmet, nmetTimeIntegrated, dtchem2
    !// Locals
    real(r8) :: frac_next, ZBOT, ZMID, ZTOP
    integer :: I,J,L, LL
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'set_blh_ij'
    !// --------------------------------------------------------------------

    !// BLH from metfields is instantaneous.
    !// Here we interpolate between this meteorolgical time step and the
    !// next. This can be done in two ways:, either to the start of each
    !// NOPS or halfway between them.
    !// Whole nops, using NROPSM=3:
    !//   NOPS=1: weight 1 for current, 0 for next.
    !//   NOPS=2: weight 2/3 for current, 1/3 for next.
    !//   NOPS=3: weight 1/3 for current, 2/3 for next.
    !//   Fraction for current weight is:
    !//     frac_cur = 1 - (NOPS - 1)/NROPSM
    !// Halfway between:
    !//   NOPS=1: weight 5/6 for current, 1/6 for next.
    !//   NOPS=2: weight 1/2 for current, 1/2 for next.
    !//   NOPS=3: weight 1/6 for current, 5/6 for next.
    !//   Fraction for current weight is:
    !//     frac_cur = (NROPSM - (NOPS-1) + 1/2) / NROPSM
    !//            = 1 - (NOPS - 0.5) / NROPSM
    !// Since CNVBLD (pbl_mixing.f90) is done several times (NSUB) during
    !// a NOPS, it should interpolate every time step (using
    !// BLH_CUR and BLH_NEXT). Here we set BLH at each NOPS - this BLH is
    !// only used for output (profiles and such) where instantaneous values
    !// are given for each hour (NOPS).
    !// To save one operation, we calculate frac_next:
    !// Halfway:
    !// frac_next = (real(NOPS, r8) - 0.5_r8) / real(NROPSM, r8)
    !// Whole step:
    !// frac_next = (real(NOPS, r8) - 1._r8) / real(NROPSM, r8)


    !// Integrate linearly in time; CCYC time step is dtchem2, and we
    !// need the fraction at half of that. For this we use the
    !// nmetTimeIntegrated and dtmet:
    frac_next = (nmetTimeIntegrated + 0.5 * dtchem2) / dtmet


    do J = MPBLKJB(MP), MPBLKJE(MP)
       do I = MPBLKIB(MP), MPBLKIE(MP)

          BLH(I,J) = (1._r8 - frac_next) * BLH_CUR(I,J) + frac_next * BLH_NEXT(I,J)

          !// Find model level for this
          LL = 1
          ZTOP = ZOFLE(2,I,J) - ZOFLE(1,I,J)
          ZBOT = 0._r8
          do L = 2, LPAR
             ZTOP = ZOFLE(L+1,I,J) - ZOFLE(1,I,J)
             ZBOT = ZOFLE(L,I,J) - ZOFLE(1,I,J)
             ZMID = 0.5_r8 * (ZTOP + ZBOT)
             if (ZMID .gt. BLH(I,J)) then
                LL = L - 1
                exit
             end if
          end do
          LBLH(I,J) = LL
          if (LL .eq. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  'LBLH is 1!!!!',BLH(I,J),ZBOT,ZMID,ZTOP
             stop
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine set_blh_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_blh(NOPS)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR
    use cmn_ctm, only: NROPSM
    use cmn_met, only: BLH, BLH_CUR, BLH_NEXT, LBLH, ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NOPS
    real(r8) :: frac_next, ZBOT, ZMID, ZTOP
    integer :: I,J,L, LL
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'set_blh'
    !// --------------------------------------------------------------------

    !// BLH from metfields is instantaneous.
    !// Here we interpolate between this meteorolgical time step and the
    !// next. This can be done in two ways:, either to the start of each
    !// NOPS or halfway between them.
    !// Whole nops, using NROPSM=3:
    !//   NOPS=1: weight 1 for current, 0 for next.
    !//   NOPS=2: weight 2/3 for current, 1/3 for next.
    !//   NOPS=3: weight 1/3 for current, 2/3 for next.
    !//   Fraction for current weight is:
    !//     frac_c = 1 - (NOPS - 1)/NROPSM
    !// Halfway between:
    !//   NOPS=1: weight 5/6 for current, 1/6 for next.
    !//   NOPS=2: weight 1/2 for current, 1/2 for next.
    !//   NOPS=3: weight 1/6 for current, 5/6 for next.
    !//   Fraction for current weight is:
    !//     frac_cur = (NROPSM - (NOPS-1) + 1/2) / NROPSM
    !//            = 1 - (NOPS - 0.5) / NROPSM
    !// Since CNVBLD (pbl_mixing.f90) is done several times (NSUB) during
    !// a NOPS, it should interpolate every time step (using
    !// BLH_CUR and BLH_NEXT). Here we set BLH at each NOPS - this BLH is
    !// only used for output (profiles and such) where instantaneous values
    !// are given for each hour (NOPS).
    !// To save one operation, we calculate frac_next:
    !// Halfway:
    !// frac_next = (real(NOPS, r8) - 0.5_r8) / real(NROPSM, r8)
    !// Whole step:
    frac_next = (real(NOPS, r8) - 1._r8) / real(NROPSM, r8)

    BLH(:,:) = (1._r8 - frac_next) * BLH_CUR(:,:) + frac_next * BLH_NEXT(:,:)

    !// Find model level for this
    do J = 1, JPAR
       do I = 1, IPAR
          LL = 1
          ZTOP = ZOFLE(2,I,J) - ZOFLE(1,I,J)
          ZBOT = 0._r8
          do L = 2, LPAR
             ZTOP = ZOFLE(L+1,I,J) - ZOFLE(1,I,J)
             ZBOT = ZOFLE(L,I,J) - ZOFLE(1,I,J)
             ZMID = 0.5_r8 * (ZTOP + ZBOT)
             if (ZMID .gt. BLH(I,J)) then
                LL = L - 1
                exit
             end if
          end do
          LBLH(I,J) = LL
          if (LL .eq. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  'LBLH is 1!!!!',BLH(I,J),ZBOT,ZMID,ZTOP
             stop
          end if
       end do
    end do
    !// --------------------------------------------------------------------
  end subroutine set_blh
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine check_lmtrop()
    !// --------------------------------------------------------------------
    !// Must check if LMTROP is set before E90 initialisation.
    !//
    !// Amund Sovde Haslerud, September 2017
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: STT, AIR, ETAA, ETAB, &
         MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: Ne90, TMASSMIX2MOLMIX
    use cmn_met, only: P,T,PVU,ZOFLE
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: MP
    real(r8) :: rtmp
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'check_lmtrop'
    !// --------------------------------------------------------------------

    !// Skip this if E90-tracer is not included
    if (Ne90 .le. 0) return

    !// Is perhaps E90 already set?
    rtmp = maxval(STT(:,:,:,Ne90) / AIR(:,:,:)) * TMASSMIX2MOLMIX(Ne90)
    if (rtmp .gt. 1.e-12_r8) return

    if (maxval(LMTROP) .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
           ': LMTROP will be initialised from PVU before E90 initialisation'
       !// Set LMTROP
       do MP = 1, MPBLK
          !// PVU-based tropopause
          call tp_pvu_ij(LMTROP,P,ETAA,ETAB,T,PVU,ZOFLE,IPAR,JPAR,LPAR, &
               IDBLK,JDBLK,MPBLK,MPBLKIB,MPBLKIE,MPBLKJB,MPBLKJE,MP)
       end do
    end if



    !// --------------------------------------------------------------------
  end subroutine check_lmtrop
  !// ----------------------------------------------------------------------

  !// ------------------------------------------------------------------
end module physics_oslo
!//=========================================================================

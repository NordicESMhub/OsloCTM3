!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routines for setting up grid.
!//=========================================================================
module grid
  !//-----------------------------------------------------------------------
  !// MODULE: grid
  !// DESCRIPTION: Routines for setting up grid.
  !//              This is the UCI p-grid.f, but slightly modified
  !//              to CTM3 f90.
  !//
  !// Contains
  !//   subroutine SET_GRID
  !//   subroutine LABELG
  !//   subroutine DIAGBLK
  !//   subroutine GAUSST2
  !//   subroutine DBLDBL
  !//   subroutine AIRSET
  !//   subroutine SURF_IN
  !//   subroutine gridICAO
  !//
  !//-----------------------------------------------------------------------
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'grid.f90'
  !// ----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains



  !//-----------------------------------------------------------------------
  subroutine SET_GRID(GM0000,JMPOLR)
    !//---------------------------------------------------------------------
    !// General grid set up, global only, handles EC and GISS grids
    !// Latitude (J=1:JPAR) grid:
    !// LGAUGRD = .true. - use Gaussian grid for latitude, else uniform grid
    !// JMPOLR (regualr grid only)
    !//         =0 =>polar box same delta-lat as others,
    !//         =1 =>polar box is halfsize (eg GISS JPAR=46)
    !// Longitude (I=1:IM) grid:   always regular, uniform spacing
    !// GM0000 = I-coord of the Greenwich Merid: 
    !/           1.5 => mid-box[1] on GM
    !//          IPARW/2+1 => left-box[1] at Dateline
    !// Currently for EC, GM0000=1.5;   for GISS4x5, GM0000=37.0 & JMPOLR=1
    !//
    !// Pressure/Altitude (L=1:LM) grid:  based on eta coords
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPARW, JPARW, IPAR, JPAR, LGAUGRD, &
         NTRUNW, NMMAX, IDGRD, JDGRD
    use cmn_ctm, only: LDEG, IMEPZ, IDTLN, WGLYE, WGLYG, ALP, ALPV, &
         DD, SS, TRIG, IFAX, TRIGU, AREAXYW, AREAXY, &
         YGRD, YDGRD, XGRD, XDGRD, YEDG, YDEDG, XEDG, XDEDG, &
         YDGRDW, XDGRDW, YDEDGW, XDEDGW, &
         ZDEGI, ZDEGJ, IMAP, JMAP, &
         DISTX, DISTY, &
         AREAG, TLAT, TLATE, TLNG, TLNGE
    use cmn_met, only: metTYPE, PMEAN, PMEANW
    use cmn_parameters, only: CPI, C2PI, CPI180, ZPI180, A0
    use utilities, only: ctmExitC
    use ncutils, only: get_netcdf_var_1d
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8),  intent(in):: GM0000
    integer, intent(in):: JMPOLR
    !// Locals
    integer :: I,J,IMB2,JMB2
    real(r8) :: DEL_JD,DEL_I,DEL_ID

    real(r8) :: WGAULAT(JPARW),WGAUWT(JPARW), WYEDG1P1(JPARW+1)
    real(r8) :: YEDG1P1(JPAR+1)
    real(r8), dimension(JPARW) :: YGRDW
    real(r8), dimension(JPARW+1) :: YEDGW
    real(r8), dimension(IPARW) :: XGRDW
    real(r8), dimension(IPARW+1) :: XEDGW

    !// lat/lon from NorESM, real native resolution
    real(r8), allocatable, dimension(:) :: inLat, inLatE
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='set_grid'
    !// --------------------------------------------------------------------

    write(6,'(a)') '--------------------------------------------'// &
         '---------------------------'
    write(6,'(a)') f90file//':'//subr//': Initialising horizontal grid'

    JMB2 = JPARW/2
    IMB2 = IPARW/2
    if (mod(IPARW,2) .ne. 0) &
         call ctmExitC(f90file//':'//subr//': IPARW must be even')
    !// check EPZ fit:
    do J = 1, JPAR
       if (mod(IPARW, IMEPZ(J)) .ne. 0) &
            call ctmExitC(f90file//':'//subr//': IMEPZ/IM misfit')
    end do


    !// LONGITUDE grid is always regular and set by IPARW and GM000
    !//---------------------------------------------------------------------
    DEL_ID = 360._r8 / real(IPARW, r8)
    do I = 1, IPARW+1
       XDEDGW(I) = (real(I, r8) - GM0000) * DEL_ID
    end do
    do I = 1, IPARW
       XDGRDW(I) = 0.5_r8 * (XDEDGW(I) + XDEDGW(I+1)) 
    end do
    do I = 1, IPARW+1
       XEDGW(I) = XDEDGW(I) * CPI180
    end do
    do I = 1, IPARW
       XGRDW(I) = XDGRDW(I) * CPI180
    end do
    !// Will skip negative values for labeling (produces 170W instead of 190E)
    !// Not important.
    !// for computation XEDGW must be montonic, for labeling shift to
    !// neg. longit.
    !do I = 1, IPARW+1
    !   if (XDEDGW(I) .gt. 180._r8) XDEDGW(I) = XDEDGW(I) - 360._r8
    !end do
    !do I = 1, IPARW
    !   if (XDGRDW(I) .gt. 180._r8) XDGRDW(I) = XDGRDW(I) - 360._r8
    !end do
    !// set dateline box
    IDTLN = mod(int(GM0000) + IPAR/2 - 1, IPAR) + 1


    !// LATITUDE grid can be Gauss-pt, or regular, incl. half-size boxes
    !// at pole,
    !// WYEDG1P1 = sin(latitude) of box edges (-1 to +1) [temporary]
    !//---------------------------------------------------------------------

    if (mod(JPARW, 2) .ne. 0)  &
         call ctmExitC(f90file//':'//subr//': JPARW must be even')


    !// Meteorology type (EC or other)
    !//---------------------------------------------------------------------
    if (METTYPE(1:6) .eq. 'norESM')  then
       !// How to use the grids
       !// Most data are given on lat/lon parameters, which should
       !// be grid centers. However, lat starts at 90S and ends at 90N.
       !// This is because the polar cap is treated as a single box
       !// in transport, so each polar pie piece has the extention of
       !// a half grid box.
       !//
       !// All variables on lat/lon are given as grid center values,
       !// e.g. rain and convective fluxes.
       !//
       !// U and V are also given on lat/lon, while US and VS are
       !// given on slat/slon (the staggered values). However,
       !// US is (lon/slat) and VS is (slon,lat), so this is not
       !// what we want.
       !//
       !// The polar pies will be very small, so I will combine it with
       !// the polar-most grid box. Hence, if NorESM has 96 latitudes
       !// the CTM will use 94.
       !// First read NorESM native data
       call get_netcdf_var_1d('norESMgrid.nc', 'lat', inLat)
       call get_netcdf_var_1d('norESMgrid.nc', 'ilat', inLatE)
       !// New edges from pole to pole
       YDEDGW(1) = -90._r8
       !// Skip the next-to-polar-pie box
       YDEDGW(2:JPARW) = inLatE(2:JPARW)
       YDEDGW(JPARW+1) = 90._r8
       do J = 1,JPARW+1
          YEDGW(J) = YDEDGW(J) * CPI180
          WYEDG1P1(J) = sin(YEDGW(J))
          write(6,'(i4,f12.5)') J,YDEDGW(J)
       end do

       !// Skip the next-to-polar-pie box
       !// and calculate new center value
       YDGRDW(1) = 0.5_r8 * (-90._r8 + inLatE(2))
       YDGRDW(2:JPARW-1) = inLat(3:JPARW) !// Skip the next-to-polar-pie box
       YDGRDW(JPARW) = 0.5_r8 * (inLatE(JPARW) + 90._r8)

       do J = 1,JPARW
          YGRDW(J) = YDGRDW(J) * CPI180
          WGAULAT(J) = sin(YGRDW(J))
          write(6,'(i4,f12.5)') J,YDGRDW(J)
       end do

       do J = 1,JPARW+1
          WGLYE(J) = Asin(WYEDG1P1(J))
       enddo
       do J = 1,JPARW
          WGLYG(J) = Asin(WGAULAT(J))
       enddo

       deallocate(inLat, inLatE)
       !// Now two questions remain:
       !// 1. What to do with winds
       !//    U needs to be shifted (U(I) is to be flux into box I)
       !//      At poles, the grid point shifts halfway towards meridional
       !//      staggered location.
       !//    V: In CTM3, V(J) is wind (converted to flux) into grid box J.
       !//       It is by default set to zero at the poles???
       !//       NorESM-V thus has to be shifted polewards, so that
       !//       V(1) = inV(1) !// wind into pie piece
       !//       do J=2,JPAR
       !//          V(J) = 0.5d0 * (inV(J-1) + inV(J))
       !//       end do
       !// 2. What to do with grid center polar cap values for other
       !//    species?
       !//    I will use the polar pie pieces as though they are given
       !//    in their new centers.
    else if (METTYPE(1:9) .eq. 'ECMWF_IFS' .or. &
         METTYPE(1:10) .eq. 'ECMWF_oIFS') then
       !// EC grids
       !//------------------------------------------------------------------
       if (.not. LGAUGRD) then

          !// Regular latitude grid
          DEL_JD = 180._r8 / real(JPARW - JMPOLR, r8)
          YDEDGW(1) = -90._r8
          YDEDGW(2) = -90._r8 + DEL_JD / real(1 + JMPOLR, r8)
          do J = 3, JPARW
             YDEDGW(J) = YDEDGW(2) + real(J-2, r8) * DEL_JD
          end do
          YDEDGW(JPARW+1) = 90._r8
          do J = 1, JPARW+1
             YEDGW(J) = YDEDGW(J)*CPI180
             WYEDG1P1(J) = sin(YDEDGW(J)*CPI180)
          end do
          do J = 1, JPARW
             YGRDW(J) = asin(0.5_r8 * (WYEDG1P1(J) + WYEDG1P1(J+1)))
             YDGRDW(J) = ZPI180 * YGRDW(J)
          end do
          do J = 1, JPARW+1
             WGLYE(J) = Asin(WYEDG1P1(J))
          end do
          do J = 1, JPARW
             WGLYG(J) = YGRDW(J)
          end do

       else

          !// Gaussian latitude grid
          call GAUSST2(JPARW, WGAULAT, WGAUWT)

          WYEDG1P1(1) = -1._r8
          do J = 2, JMB2
             WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
          end do
          WYEDG1P1(JMB2+1) = 0._r8
          do J = JMB2+2, JPARW+1
             WYEDG1P1(J) = -WYEDG1P1(JPARW+2-J)
          end do

          do J = 1, JPARW+1
             WGLYE(J) = Asin(WYEDG1P1(J))
          end do
          do J = 1, JPARW
             WGLYG(J) = Asin(WGAULAT(J))
          end do

          do J = 1, JPARW+1
             YEDGW(J) = asin(WYEDG1P1(J))
             YDEDGW(J) = ZPI180 * YEDGW(J)
          end do
          do J = 1, JPARW
             YGRDW(J) = asin(WGAULAT(J))
             YDGRDW(J) = ZPI180 * YGRDW(J)
          end do

       end if !// if (.not. LGAUGRD) then



       !// EC: set up Spectral to Grid transform functions: Assoc-Leg Polys
       do J = 1, JMB2
          call LEGGEN(ALP(1,J),YGRDW(J),NTRUNW,NMMAX)
       end do
       do J = 1, JMB2
          call LEGGEN(ALPV(1,J),YEDGW(J+1),NTRUNW,NMMAX)
       end do
       !// EC: transforms from Vorticiy and Divergence to U and V
       call FFT_DDSS (NTRUNW,A0,NMMAX,DD,SS)
       !// EC: FFT transforms, use grid
       call FFT_SET (TRIG,IFAX,IPARW)
       !// EC: trig functions for U (shift to boundaries)
       DEL_I = C2PI / real(IPARW, r8)
       do I = 1, IMB2
          TRIGU(2*I-1) = cos(DEL_I*(real(I, r8)-0.5_r8))
          TRIGU(2*I  ) = sin(DEL_I*(real(I, r8)-0.5_r8))
       end do

    else
       write(6,'(a)') f90file//':'//subr// &
            ': Unknown metTYPE: '//trim(metTYPE)
       stop
    end if !// if (METTYPE(1:2) .eq. 'EC')  then

    if (.false.) then
       !// REGULAR grid can be found as this
       !//------------------------------------------------------------------
       DEL_JD = 180._r8 / real(JPARW - JMPOLR, r8)
       YDEDGW(1) = -90._r8
       YDEDGW(2) = -90._r8 + DEL_JD/real(1 + JMPOLR, r8)
       do J = 3, JPARW
          YDEDGW(J) = YDEDGW(2) + real(J-2, r8) * DEL_JD
       end do
       YDEDGW(JPARW+1) = 90._r8
       do J = 1, JPARW+1
          YEDGW(J) = YDEDGW(J) * CPI180
          WYEDG1P1(J) = sin(YDEDGW(J) * CPI180)
       end do
       do J = 1,JPARW
          YGRDW(J) = asin(0.5_r8*(WYEDG1P1(J)+WYEDG1P1(J+1)))
          YDGRDW(J) = ZPI180 * YGRDW(J)
       end do
    end if !// if (METTYPE(1:2) .eq. 'EC')  then


    !// Area of grid boxes
    !//---------------------------------------------------------------------
    do J = 1,JPARW
       do I = 1,IPARW
          AREAXYW(I,J)= A0*A0 * (XEDGW(I+1) - XEDGW(I)) &
                              * (WYEDG1P1(J+1) - WYEDG1P1(J))
       end do
    end do

    !// Define lat-lon of ctm grid and weightings if horizontal
    !// resolution degraded
    !//---------------------------------------------------------------------
    call DBLDBL(YGRD,YDGRD,XGRD,XDGRD,YEDG,YDEDG,XEDG,XDEDG, &
                YGRDW,YDGRDW,XGRDW,XDGRDW,YEDGW,YDEDGW,XEDGW,XDEDGW, &
                ZDEGI,ZDEGJ,IMAP,JMAP, &
                IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD,LDEG)

    YEDG1P1(1) = -1._r8
    do J = 2, JPAR
       YEDG1P1(J) = sin(YEDG(J))
    end do
    YEDG1P1(JPAR+1) = 1._r8

    !// DIAGNOSTICS AND WEIGHTING FACTORS
    !//---------------------------------------------------------------------
    !//>>>>>this new gridding does not yet deal with nested, doubled grids!

    !// AREAS:  use the IPAR,JPAR rather than IPARW,JPARW in case of re-grid
    do J = 1, JPAR
       DISTY(J) = A0*(YEDG(J+1)-YEDG(J))
    end do
    do J = 1, JPAR
       DISTX(J) = A0 * (XEDG(2) - XEDG(1)) * 0.5_r8 &
                     * (cos(YEDG(J)) + cos(YEDG(J+1)))
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          AREAXY(I,J) = A0*A0 * (XEDG(I+1) - XEDG(I)) &
                              * (YEDG1P1(J+1) - YEDG1P1(J))
       end do
    end do

    !// Global total area
    AREAG = sum(AREAXY)

    !// Labels
    !//---------------------------------------------------------------------
    do J = 1, JPAR
       call LABELG(YDGRD(J),TLAT(J),1)
    end do
    do J = 1, JPAR+1
       call LABELG(YDEDG(J),TLATE(J),1)
    end do
    do I = 1, IPAR
       call LABELG(XDGRD(I),TLNG(I),2)
    end do
    do I = 1, IPAR+1
       call LABELG(XDEDG(I),TLNGE(I),2)
    end do

    !// Write-out grid data
    !//---------------------------------------------------------------------
    !// will need to shift absolute global grid to actual window calc later!
    write(6,*) '  ------------------grid----------------------'
    write(6,'(a,1p,e15.8)') ' AREA globe(m^2)    : ', AREAG
    write(6,'(A)') 'Area of grid boxes (m^2) 1:JPAR '
    write(6,'(1P10E10.3)') (AREAXY(1,J),J=1,JPAR)

    write(6,'(a)') ' MID-POINT of GRID BOXES'
    write(6,'(a,i5)') '     J=1:',JPAR
    write(6,'(5(i5,f9.5,f6.1,a5))') (J, YGRD(J), YDGRD(J), TLAT(J), J=1,JPAR)
    write(6,'(a,i5)') '     I=1:',IPAR
    write(6,'(5(i5,f9.5,f6.1,a5))') (I, XGRD(I), XDGRD(I), TLNG(I), I=1,IPAR)
    write(6,'(a,i8)') '  IDATELINE=', IDTLN
    write(6,*) ' EDGES of GRID'
    write(6,'(a,i5)') '     J=1:',JPAR+1
    write(6,'(5(i5,f9.5,f6.1,a5))') (J, YEDG(J), YDEDG(J), TLATE(J), J=1,JPAR+1)
    write(6,'(a,i5)') '     I=1:',IPAR+1
    write(6,'(5(i5,f9.5,f6.1,a5))') (I, XEDG(I), XDEDG(I), TLNGE(I), I=1,IPAR+1)
    write(6,'(a)') ' IMEPZ = polar box averaging over met fields'
    write(6,'(8(i6,a4,i4))') (J, TLAT(J), IMEPZ(J), J=1,JPAR)

    !//---------------------------------------------------------------------
  end subroutine SET_GRID
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SET_GRID_VERT()
    !//---------------------------------------------------------------------
    !// General grid set up, global only and vertical only.
    !// Uses hybrid coordinates (ETAAW and ETABW), along with
    !// PMEANW and PMEAN.
    !//
    !// Separated out vertical part from SET_GRID, since mean pressure
    !// needed horizontal part om SET_GRID and the vertical part needed
    !// mean pressure.
    !//
    !// Ole Amund Sovde, March 2016
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPARW, JPARW, LPARW, IPAR, JPAR, LPAR, &
         IDGRD, JDGRD
    use cmn_ctm, only: AREAXY, &
         ETAAW, ETABW, ETAA, ETAB, &
         LMMAP, XLMMAP, XYZA, XYZB, XYA, XYB, PIJL, &
         ZEDG, ZGRD, &
         TALT, TALTE
    use cmn_met, only: PMEAN, PMEANW
    use cmn_parameters, only: G0
    use ncutils, only: get_netcdf_var_1d
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    !// Locals
    integer :: I,J,L,LW
    real(r8) :: SUMXYZ, AIRGLB
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='set_grid_vert'
    !// --------------------------------------------------------------------

    write(6,'(a)') '--------------------------------------------'// &
         '---------------------------'
    write(6,'(a)') f90file//':'//subr//': Initialising vertical grid'

    !// PRESSURE/ALTITUDE grid, remap the global wind grid (W)
    !//---------------------------------------------------------------------
    !// if LPAR < LPARW (met field layers) use LMMAP to collapse
    !// met-layers on CTM L
    !// Calculate new ETAA & ETAB for CTM(1:LM) with remapped vertical
    !// grid(1:LPARW)
    do LW = 1, LPARW+1
       ETAAW(LW) = 0.01_r8 * ETAAW(LW)
    end do
    ETAA(1) = ETAAW(1)
    ETAB(1) = ETABW(1)
    do LW = 1, LPARW
       L = LMMAP(LW)
       ETAA(L+1) = ETAAW(LW+1)
       ETAB(L+1) = ETABW(LW+1)
    end do
    !// LPAR is top anyway
    !LM = LMMAP(LPARW)

    !// Calculate weighting for LW=1:LPARW to merged CTM layer L=1:LM
    do LW = 1, LPARW
       L = LMMAP(LW)
       XLMMAP(LW) = ((ETAAW(LW) - ETAAW(LW+1))  &
            + 1.e3_r8 * (ETABW(LW) - ETABW(LW+1)) ) / ((ETAA(L) - ETAA(L+1)) &
            + 1.e3_r8 * (ETAB(L) - ETAB(L+1)) )
    end do


    !// Air mass factors for each grid box:
    !// air mass(I,J,L) in kg = XYZA() + XYZB()*surface pressure (mbar!!)
    SUMXYZ  = 0._r8
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             XYZA(I,J,L) = (ETAA(L) - ETAA(L+1)) * AREAXY(I,J) * 1.e2_r8 / G0
             XYZB(I,J,L) = (ETAB(L) - ETAB(L+1)) * AREAXY(I,J) * 1.e2_r8 / G0
             SUMXYZ  = SUMXYZ + XYZB(I,J,L)
          end do
       end do
    end do
    do J = 1, JPAR
       do I = 1,IPAR
          XYZA(I,J,LPAR+1) = ETAA(LPAR+1) * AREAXY(I,J) * 1.e2_r8 / G0
          XYZB(I,J,LPAR+1) = ETAB(LPAR+1) * AREAXY(I,J) * 1.e2_r8 / G0
       end do
    end do
    !// Column air mass factors XYA, XYB do not include > model top
    do J = 1, JPAR
       do I = 1, IPAR
          XYA(I,J)   = 0._r8
          XYB(I,J)   = 0._r8
          do L = 1, LPAR
             XYA(I,J) = XYA(I,J) + XYZA(I,J,L)
             XYB(I,J) = XYB(I,J) + XYZB(I,J,L) 
          end do
       end do
    end do


    AIRGLB = 0._r8
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             AIRGLB  = AIRGLB + XYZA(I,J,L) + XYZB(I,J,L) * PMEAN(I,J)
          end do
       end do
    end do

    !// Create new 3-D P-grid for interpolating sources, etc. using PMEAN
    !// PIJL(I,J,L=1:LM+1) is pressure at edges
    do L = 1, LPAR+1
       do J = 1, JPAR
          do I = 1,IPAR
             PIJL(I,J,L) =  ETAA(L) + ETAB(L) * PMEAN(I,J)
          end do
       end do
    end do

    do L = 1, LPAR+1
       ZEDG(L) = ETAA(L) + ETAB(L) * 1000._r8 
       call LABELG(ZEDG(L),TALTE(L),3)
    end do
    do L = 1, LPAR
       ZGRD(L) = 0.5_r8 * (ZEDG(L) + ZEDG(L+1))
       call LABELG(ZGRD(L),TALT(L),3)
    end do

    write(6,*) '  --------------vertical grid----------------------'
    write(6,'(a,1p,e15.8)') ' AIR mass wet (kg)  : ', AIRGLB
    write(6,'(a)') ' MID-POINT of GRID BOXES'
    write(6,'(a)')  '   LW   wting  ==> L (lbl)'
    write(6,'(i5,f10.5,i5,a5)') &
         (LW, XLMMAP(LW),LMMAP(LW),TALT(LMMAP(LW)), LW=1,LPARW)
    write(6,*) ' EDGES of GRID'
    write(6,'(a,i5)') '     L=1:',LPAR+1
    write(6,'(1x,a1,i2,a2,a5,f9.3,a4)') & 
         ('L',L-1,'.5', TALTE(L), ZEDG(L),' hPa', L=1,LPAR+1)

    !//---------------------------------------------------------------------
  end subroutine SET_GRID_VERT
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine set_mean_psfc(infile)
    !//---------------------------------------------------------------------
    !// Reads 2D annual mean surface pressure from file (PMEANW) to be
    !// used as topography (through PMEAN).
    !//
    !// Reads netCDF file and interpolates to model resolution.
    !//
    !// Ole Amund Sovde, March 2016
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPARW, JPARW, IPAR, JPAR, &
         NTRUNW, NMMAX, IDGRD, JDGRD
    use cmn_ctm, only: ZDEGI, ZDEGJ, IMAP, JMAP, &
         XDGRDW, YDGRDW, XDEDGW, YDEDGW, AREAXYW, AREAXY, &
         PMEANG, AREAG
    use cmn_met, only: PMEAN, PMEANW
    use regridding, only: E_GRID, TRUNG8
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=*), intent(in) :: infile
    integer :: I, J, I1, J1
    !// Mean surface pressure
    real(r8), allocatable :: inPMEAN(:,:),inAREA(:,:), inXDEDG(:),inYDEDG(:)
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='set_mean_psfc'
    !//---------------------------------------------------------------------

    write(6,'(a)') '--------------------------------------------'// &
         '---------------------------'
    write(6,'(a)') f90file//':'//subr// &
         ': Read mean surf press from: '//trim(INFILE)
    !// Fetch data from netcdf file
    !// lon/lat for interfaces
    call get_netcdf_var_1d( trim(infile), 'ilon', inXDEDG )
    call get_netcdf_var_1d( trim(infile), 'ilat', inYDEDG )
    I1 = size(inXDEDG)-1
    J1 = size(inYDEDG)-1
    !// surface pressure and area
    allocate( inPMEAN(I1, J1), inAREA(I1, J1) )
    call get_netcdf_var_2d( trim(infile), 'Pmean', inPMEAN, I1, J1 )
    call get_netcdf_var_2d( trim(infile), 'gridarea', inAREA, I1, J1 )
    !// Interpolate to native resolution (IPARW,JPARW)
    !// Must multiply by area and divide by area after interpolation
    inPMEAN(:,:) = inPMEAN(:,:) * inAREA(:,:)
    call E_GRID(inPmean,inXDEDG,inYDEDG,I1,J1, &
                PMEANW,XDEDGW,YDEDGW,IPARW,JPARW,1)
    deallocate( inXDEDG, inYDEDG, inPMEAN, inAREA)
    PMEANW(:,:) = PMEANW(:,:) / AREAXYW(:,:)

    !// Find PMEAN
    call TRUNG8(PMEANW, PMEAN, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)

    write(6,'(a,2f12.5)') 'Min/max native mean pressure [hPa]: ', &
         minval(pmeanw),maxval(pmeanw)
    write(6,'(a,2f12.5)') 'Min/max mean pressure [hPa]: ', &
         minval(pmean),maxval(pmean)

    !// Mean global surface pressure
    PMEANG = 0._r8
    do J = 1, JPAR
       do I = 1, IPAR
          PMEANG = PMEANG + PMEAN(I,J) * AREAXY(I,J)
       end do
    end do
    PMEANG = PMEANG / AREAG
    write(6,'(a,f15.8)') 'Global mean surface pressure [hPa]: ',PMEANG

    !//---------------------------------------------------------------------
  end subroutine set_mean_psfc
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine LABELG(X,TITLX,NDEX)
    !//---------------------------------------------------------------------
    !// Generate char*4 titles for lat(NDEX=1), long(NDEX=2) or vert(NDEX=3)
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) ::  X
    integer, intent(in) :: NDEX
    !// Output
    character(len=4), intent(out) :: TITLX
    !// Locals
    character(len=1), parameter, dimension(2)  :: TITLNS =['N','S']
    character(len=1), parameter, dimension(2)  :: TITLEW =['E','W']
    character(len=1), parameter, dimension(12) :: &
         TITLNN =[' ','1','2','3','4','5','6','7','8','9','0','-']
    real(r8) :: XABS
    integer :: I1000,I0100,I0010,I0001,I100,I010,I001,I10,I01,IXABS
    !//---------------------------------------------------------------------

    !// latitude (NDEX=1)
    XABS = abs(X)
    if (NDEX .eq. 1) then
       IXABS = int(XABS + 0.5_r8)
       I10 = IXABS / 10
       I01 = IXABS - 10 * I10
       if (I01 .eq. 0) I01 = 10
       I10 = min(10, I10)
       TITLX(1:1) = TITLNN(1)
       TITLX(2:2) = TITLNN(I10+1)
       TITLX(3:3) = TITLNN(I01+1)
       if (X .ge. 0._r8) then
          TITLX(4:4) = TITLNS(1)
       else
          TITLX(4:4) = TITLNS(2)
       end if
    end if

    !// longitude (NDEX=2)
    if (NDEX .eq. 2) then
       IXABS = int(XABS + 0.5_r8)
       I100 = IXABS / 100
       I010 = (IXABS - 100 * I100) / 10
       I001 = IXABS - 100 * I100 - 10 * I010
       if (I100 .gt. 0 .AND. I010 .eq. 0) I010 = 10
       if (I001 .eq. 0) I001 = 10
       I100 = min(10, I100)
       TITLX(1:1) = TITLNN(I100+1)
       TITLX(2:2) = TITLNN(I010+1)
       TITLX(3:3) = TITLNN(I001+1)
       if (X .ge. 0._r8) then
          TITLX(4:4) = TITLEW(1)
       else
          TITLX(4:4) = TITLEW(2)
       end if
    end if

    !// pressure (NDEX=3) - nearest mbar only for now
    if (NDEX .eq. 3) then
       if (XABS .gt. 1._r8) then
          IXABS = int(XABS + 0.5_r8)
          I1000 = IXABS / 1000
          I0100 = (IXABS - 1000 * I1000) / 100
          I0010 = (IXABS - 1000 * I1000 - 100 * I0100) / 10
          I0001 =  IXABS - 1000 * I1000 - 100 * I0100 - 10 * I0010
          I1000 = min(I1000, 10)
          if (I1000 .gt. 0 .AND. I0100 .eq. 0) I0100 = 10
          if (I0100 .gt. 0 .AND. I0010 .eq. 0) I0010 = 10
          if (I0001 .eq. 0) I0001 = 10
          TITLX(1:1) = TITLNN(I1000+1)
          TITLX(2:2) = TITLNN(I0100+1)
          TITLX(3:3) = TITLNN(I0010+1)
          TITLX(4:4) = TITLNN(I0001+1)
       else
          TITLX(1:1) = TITLNN(1)
          TITLX(2:2) = TITLNN(1)
          TITLX(3:3) = TITLNN(1)
          XABS = 10._r8 * XABS
          TITLX(3:3) = TITLNN(12)
          if(XABS .lt. 1._r8) then
             XABS = 10._r8 * XABS
             TITLX(2:2) = TITLNN(12)
             if(XABS .lt. 1._r8) then
                XABS = 10._r8 * XABS
                TITLX(1:1) = TITLNN(12)
             end if
          end if
          I0001 = int(XABS + 0.50_r8)
          I0001 = min(I0001, 9)
          TITLX(4:4) = TITLNN(I0001+1)
       end if
    end if

    !//---------------------------------------------------------------------
  end subroutine LABELG
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine LABELG2(X,TITLX,NDEX)
    !//---------------------------------------------------------------------
    !// Generate char*4 titles for lat(NDEX=1), long(NDEX=2) or vert(NDEX=3)
    !// A bit more elegant?
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) ::  X
    integer, intent(in) :: NDEX
    !// Output
    character(len=4), intent(out) :: TITLX
    !// Locals
    character(len=1), parameter, dimension(2)  :: TITLNS =['N','S']
    character(len=1), parameter, dimension(2)  :: TITLEW =['E','W']
    real(r8) :: XABS
    integer :: IABS
    !//---------------------------------------------------------------------

    !// latitude (NDEX=1)
    XABS = abs(X)
    if (NDEX .eq. 1) then
       write(TITLX(1:3),'(I3)') nint(XABS)
       if (X .ge. 0._r8) then
          TITLX(4:4) = TITLNS(1)
       else
          TITLX(4:4) = TITLNS(2)
       end if
    end if

    !// longitude (NDEX=2)
    if (NDEX .eq. 2) then
       if (X .ge. 0._r8 .and. X .le. 180._r8) then
          !// Latitude is E
          write(TITLX(1:3),'(I3)') nint(XABS)
          TITLX(4:4) = TITLEW(1)
       else if (X .gt. 180._r8) then
          !// Change to negative for W
          write(TITLX(1:3),'(I3)') nint(abs(X - 360._r8))
          TITLX(4:4) = TITLEW(2)
       else
          !// Latitude is W
          write(TITLX(1:3),'(I3)') nint(XABS)
          TITLX(4:4) = TITLEW(2)
       end if
    end if

    !// pressure (NDEX=3) - nearest mbar only for now
    if (NDEX .eq. 3) then
       if (XABS .gt. 1._r8) then
          write(TITLX(1:4),'(I4)') nint(XABS)
       else
          !// Altarnative
          TITLX(1:3) = '  -'
          if(XABS .lt. 1._r8) then
             XABS = 10._r8 * XABS
             TITLX(2:2) = '-'
             if(XABS .lt. 1._r8) then
                XABS = 10._r8 * XABS
                TITLX(1:1) = '-'
             end if
          end if
          IABS = min(nint(XABS), 9)
          write(TITLX(4:4),'(i1)') IABS

       end if
    end if

    !//---------------------------------------------------------------------
  end subroutine LABELG2
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine DIAGBLK(XLNG,YLAT,XDEDG,YDEDG, IBOX,JBOX,IPAR,JPAR,IDTLN)
    !//---------------------------------------------------------------------
    !// Finds grid box start/end indices for box diagnose.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  IPAR, JPAR, IDTLN
    real(r8), intent(in) ::   XDEDG(IPAR+1), YDEDG(JPAR+1)
    !// Input/Output
    real(r8), intent(inout) :: XLNG(2), YLAT(2)
    !// Output
    integer, intent(out) :: IBOX(2), JBOX(2)
    !// Locals
    integer :: I, J1,J2,I1,I2
    !//---------------------------------------------------------------------

    YLAT(1) = min(YDEDG(JPAR+1), max(YDEDG(1), YLAT(1)))
    YLAT(2) = min(YDEDG(JPAR+1), max(YDEDG(1), YLAT(2)))

    J1 = 2
    do while (YLAT(1) .gt. YDEDG(J1))
       J1 = J1 + 1
    end do
    J1 = J1 - 1
    J2 = J1 + 1
    do while (YLAT(2) .gt. YDEDG(J2))
       J2 = J2 + 1
    end do
    J2 = J2 - 1

    JBOX(1) = J1
    JBOX(2) = J2

    !// Note that IDTLN should be the I-box either across the dateline:  
    !//   XDEDG(IDTLN) = +179. & XDEDG(IDTLN+1) = -179.
    !// or the first box on the east side of dateline:
    !//   XDEDG(IDTLN) = +180. & XDEDG(IDTLN+1) = -178.

    !// Then, XDEDG(IDTLN+1) is the smallest (most negative) longitude.
    !// and XDEDG(IDTLN) is the largest longitude value.
    if (XLNG(1) .gt. XDEDG(IDTLN)) XLNG(1) = -180._r8
    if (XLNG(2) .gt. XDEDG(IDTLN)) XLNG(2) = -180._r8

    I1 = IDTLN + 1
    I  = mod(I1-1, IPAR) + 1
    do while (XLNG(1) .gt. XDEDG(I))
       I1 = I1 + 1
       I  = mod(I1-1, IPAR) + 1
    end do
    I1 = I1 - 1
    I2 = IDTLN + 1
    I  = mod(I2-1, IPAR) + 1
    do while (XLNG(2) .gt. XDEDG(I))
       I2 = I2 + 1
       I  = mod(I2-1, IPAR) + 1
    end do
    I2 = I2 - 1

    !// The sequence I1:I2 is monotonic, but may cross into I2>IPAR
    IBOX(1) = mod(I1-1, IPAR) + 1
    IBOX(2) = mod(I2-1, IPAR) + 1
    if (IBOX(2) .lt. IBOX(1))  IBOX(2) = IPAR + IBOX(2)

    !//---------------------------------------------------------------------
  end subroutine DIAGBLK
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine DIAG_LTSTN(XLNG,YLAT,XDEDG,YDEDG,GM0_LT,LTM1,LTM2, &
                        N,ISTN,JSTN,IPAR,JPAR,IDTLN,NOPSTL,LTSTN)
    !//---------------------------------------------------------------------
    !// Finds grid box indices for station diagnose.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  IPAR, JPAR, IDTLN, NOPSTL, N
    real(r8), intent(in) ::   XDEDG(IPAR+1),YDEDG(JPAR+1),GM0_LT(IPAR+1)
    !// Input/Output
    real(r8), intent(inout) :: XLNG, YLAT, LTM1, LTM2
    !// Output
    integer, intent(out) ::  ISTN, JSTN
    logical, intent(out) ::  LTSTN(NOPSTL)
    !// Locals
    real(r8) :: DTSTEP, LTGRD1, LTGRD2, DIFT1, DIFT2
    integer :: I, J, K, M, II, K1, K2
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='DIAG_LTSTN'
    !// --------------------------------------------------------------------

    !// check for sign change at dateline (eg, box[1] from +175. to -175.
    if (XLNG .gt. XDEDG(IDTLN)) XLNG = -180._r8

    II = IDTLN + 1
    I  = mod(II-1, IPAR) + 1
    do while (XLNG .gt. XDEDG(I))
       II = II + 1
       I  = mod(II-1, IPAR) + 1
    end do
    II = II - 1
    ISTN = mod(II-1, IPAR) + 1

    YLAT = min(YDEDG(JPAR+1), max(YDEDG(1),YLAT) )
    J = 2
    do while  (YLAT .gt. YDEDG(J))
       J = J + 1
    end do
    JSTN = J - 1

    DTSTEP = 24._r8 / real(NOPSTL, r8)

    if (LTM2 .lt. LTM1)  LTM2 = LTM2 + 24._r8
    I  = ISTN
    K  = 0
    do M = 1, NOPSTL
       LTGRD1 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I)), 24._r8)
       LTGRD2 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
       if (LTGRD2 .lt. LTGRD1)   LTGRD2 = LTGRD2 + 24._r8

       if (LTM1 .le. LTGRD1 .and. LTGRD2 .le. LTM2) then
          LTSTN(M) = .true.
          K  = K + 1
          !//write(6,'(A,5I5,2(3X,2f6.2),F9.2)') 'NOPSTL/N/K/I/M', &
          !//     NOPSTL,N,K,I,M,LTM1,LTM2,LTGRD1,LTGRD2,XLNG
       else if (LTM2 .ge. 24._r8)  then
          if (LTM1 .le. LTGRD1 + 24._r8 .and. LTGRD2 + 24._r8 .le. LTM2) then
             LTSTN(M) = .true.
             K  = K + 1
             !//write(6,'(A,5I5,2(3X,2f6.2),F9.2)') 'NOPSTL/N/K/I/M', &
             !//     NOPSTL,N,K,I,M,LTM1,LTM2,LTGRD1,LTGRD2,XLNG
          else
             LTSTN(M) = .false.
          end if
       else
          LTSTN(M) = .false.
       end if
    end do

    K2  = 0
    if (K .eq. 0) then
       do M = 1, NOPSTL
          K1  = M + 1
          LTGRD1 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
          LTGRD2 = dmod((real(K1-1, r8) * DTSTEP + GM0_LT(I)), 24._r8)
          if (LTGRD2 .lt. LTGRD1)   LTGRD2 = LTGRD2 + 24._r8

          if (LTGRD1 .le. LTM1 .and. LTM2 .le. LTGRD2) then
             K2  = M
             DIFT1 = LTM1 - LTGRD1
             DIFT2 = LTGRD2 - LTM2
          end if
       end do

       if (DIFT2 .gt. DIFT1) then
          write(6,'(A,I4,A,2F6.2,A,2F6.2)') 'STN', N, &
               ' Local time diag. interval is changed from ', &
               LTM1,LTM2,'  to ', &
               dmod((real(K2-1, r8) * DTSTEP + GM0_LT(I)), 24._r8), LTM2
          LTSTN(K2) = .true.
          LTM1  = dmod((real(K2-1, r8) * DTSTEP + GM0_LT(I)), 24._r8)

       else

          write(6,'(A,I4,A,2F6.2,A,2F6.2)') 'STN', N, &
               ' Local time diag. interval is changed from ', &
               LTM1,LTM2,'  to ', &
               LTM1,dmod((real(K2, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
          LTSTN(K2+1) = .true.
          LTM2  = dmod((real(K2, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
       end if
    end if
    !// write(6,'(I2,24(1X,4L1))') N,(LTSTN(K),K=1,NOPSTL)

    if (K .eq. 0 .and. K2 .eq. 0)  then
       write(6,'(A,I5)') f90file//':'//subr// &
            ': check local time diag. interval at STN',N
       write(6,'(A)') '*** job killed ***'
       stop
    end if

    !//---------------------------------------------------------------------
  end subroutine DIAG_LTSTN
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine DIAG_LTGL(GM0_LT,LTM1,LTM2,IPAR,NOPSTL,NRGLTD,LTGBL)
    !//---------------------------------------------------------------------
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  IPAR, NOPSTL
    real(r8), intent(in) ::   GM0_LT(IPAR+1)
    !// Input/Output
    real(r8), intent(inout) :: LTM1, LTM2
    !// Ouput
    integer, intent(out) ::  NRGLTD
    logical, intent(out) ::  LTGBL(IPAR,100)
    !// Locals
    real(r8) :: DTSTEP, LTGRD1, LTGRD2, GINTVL, LTM2P
    integer :: I, K, M, MINTVL, K2
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='DIAG_LTGL'
    !// --------------------------------------------------------------------

    !// check for sign change at dateline (eg, box[1] from +175. to -175.

    if (LTM2 .lt. LTM1)  LTM2 = LTM2 + 24._r8
    LTGBL(:,:) = .false.
    DTSTEP = 24._r8 / real(NOPSTL, r8)
    GINTVL = GM0_LT(3) - GM0_LT(2)    ! time interval of one grid box
    LTM2   = max(LTM2, LTM1 + DTSTEP - GINTVL)
    MINTVL = int( (LTM2 - LTM1) / GINTVL ) + 1
    LTM2   = LTM1 + (MINTVL + 0.2_r8) * GINTVL
    NRGLTD = max(int((LTM2 - LTM1 + 0.99_r8 * GINTVL) / DTSTEP), 1)
    write(6,'(A,E17.9)') ' Time interval of one grid',GINTVL
    write(6,'(A,I3)') ' Number of recordings per day', NRGLTD
    write(6,'(A,2X,2F9.4)') &
         ' <<<<< adjusted local time interval for 3-D diag.',LTM1,LTM2

    do I = 1, IPAR
       K  = 0
       M  = 0
       do while (K .lt. NRGLTD .and. M .lt. NOPSTL)
          M = M + 1
          LTGRD1 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I)), 24._r8)
          LTGRD2 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
          if (LTGRD2 .lt. LTGRD1)   LTGRD2 = LTGRD2 + 24._r8

          if (LTM1 .lt. LTGRD2 .and. LTGRD2 .le. LTM2) then
             LTGBL(I,M) = .true.
             K  = K + 1
             !//write(6,'(A,4I5,2(3X,2f6.2),F9.2)') 'NOPSTL/K/I/M', &
             !//     NOPSTL,K,I,M,LTM1,LTM2,LTGRD1,LTGRD2
          else if (LTM2 .ge.  24._r8) then
             if (LTM1 .lt. LTGRD2 + 24._r8 .and. LTGRD2 + 24._r8 .le. LTM2) then
                LTGBL(I,M) = .true.
                K  = K + 1
                !//write(6,'(A,4I5,2(3X,2f6.2),F9.2)') 'NOPSTL/K/I/M', &
                !//     NOPSTL,K,I,M,LTM1,LTM2,LTGRD1,LTGRD2
             else
                LTGBL(I,M) = .false.
             end if
          else
             LTGBL(I,M) = .false.
          end if
       end do

       if (K .lt. NRGLTD) then
          LTM2P = LTM2 + 0.99_r8 * GINTVL
          K2  = 0
          M   = 0
          do while (K2 .lt. NRGLTD .and. M .lt. NOPSTL)
             M = M + 1
             LTGRD1 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I)), 24._r8)
             LTGRD2 = dmod((real(M-1, r8) * DTSTEP + GM0_LT(I+1)), 24._r8)
             if (LTGRD2 .lt. LTGRD1)   LTGRD2 = LTGRD2 + 24._r8

             if (LTM1 .lt. LTGRD2 .and. LTGRD2 .le. LTM2P) then
                LTGBL(I,M) = .true.
                K2  = K2 + 1
                !//write(6,'(A,4I5,2(3X,2f6.2),F9.2)') 'NOPSTL/K2/I/M', &
                !//     NOPSTL,K2,I,M,LTM1,LTM2P,LTGRD1,LTGRD2
             else if (LTM2P .ge. 24._r8) then
                if (LTM1 .lt. LTGRD2 + 24._r8 .and. LTGRD2 + 24._r8 .le. LTM2P) then
                   LTGBL(I,M) = .true.
                   K2  = K2 + 1
                   !//write(6,'(A,4I5,2(3X,2f6.2),F9.2)') 'NOPSTL/K2/I/M', &
                   !//     NOPSTL,K2,I,M,LTM1,LTM2P,LTGRD1,LTGRD2
                else
                   LTGBL(I,M) = .false.
                end if
             else
                LTGBL(I,M) = .false.
             end if
          end do
          if (K2 .eq. NRGLTD) then
             write(6,'(A,I3,A,2X,2F9.4)') ' At I=',I, &
                  ' Local time interval is expanded',LTM1,LTGRD2
          else
             write (6,'(A,I3,A,I3,A)') f90file//':'//subr//': At I=',I, &
                  ' Number of counts', NRGLTD,' cannot be reached'
             write(6,'(A)') '*** job killed ***'
             stop
          end if
       end if
    end do

    !//---------------------------------------------------------------------
  end subroutine DIAG_LTGL
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine GAUSST2(NG,XP,WT)
    !//---------------------------------------------------------------------
    !// Calculate NG Gaussian quadrature points (XP) & weights (WT)
    !// for interval (X1,X2)  from GISS-Lacis code.
    !// tested against EC version GAUAW, simpler, same to 1.d-13 or better)
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  NG
    !// Output
    real(r8), dimension(NG), intent(out) :: XP, WT
    !// Parameters
    real(r8), parameter :: PI = 3.141592653589793_r8
    real(r8), parameter :: PS = 1.013211836423378e-1_r8
    real(r8), parameter :: DXL = 1.e-16_r8
    real(r8), parameter :: X1 = -1._r8,  X2 = 1._r8
    !// Locals
    real(r8) :: XMID,XHAF,DNG,DN,DM,DI,C,Z,ZZ,PN,PNI,PNJ,PNK,DX,X
    integer :: NN,N2,N,I,J
    logical :: do130
    !//---------------------------------------------------------------------

    XMID = (X2+X1) / 2._r8
    XHAF = (X2-X1) / 2._r8
    DNG = NG
    NN = NG / 2
    N2 = NN * 2

    !//if (N2 .eq. NG) goto 110
    if (N2 .ne. NG) then

       XP(NN+1) = XMID
       WT(NN+1) = 1._r8
       if (NG .lt. 2) return !// Done

       PN = 1._r8

       !N = 0
       !100 N = N+2
       !    DN = N
       !    DM = DN - 1.d0
       !    PN = PN * (DM / DN)
       !    if (N .lt. N2) goto 100
       !// The 100-loop tests on N<N2. Since N=0 at start, it can be
       !// replaced with a do while loop.
       !// Better to use do N = 2, N2, 2 since N2 is always >=2.
       !do while (N .lt. N2)
       !   N = N + 2
       !   DN = N
       !   DM = DN - 1.d0
       !   PN = PN * (DM / DN)
       !end do
       do N = 2, N2, 2
          DN = real(N, r8)
          DM = DN - 1._r8
          PN = PN * (DM / DN)
       end do

       WT(NN+1) = 2._r8 * XHAF / (DNG * PN)**2
    end if !// if (N2 .ne. NG) then

    !// 110: we get straight here if N2.eq.NG
    I = 0
    C = PI / sqrt(DNG * (DNG + 1._r8) + 0.5_r8 - PS) / 105._r8

    !// 120 starts at I=1, with increments of 1 and goes on as long
    !// as I<NN, so that in the final round I=NN.
    !// NN is fixed throughout the loop.
    do I = 1, NN
       DI = I
       Z = PS / (4._r8 * DI - 1._r8)**2
       ZZ = (105._r8 + Z * (210._r8 - Z &
                          * (2170._r8 - Z * (105812._r8 - 12554474._r8 * Z))))
       X = DCOS(ZZ * C * (DI - 0.25_r8))

       !// 130 sets N = 1, but loops as long as abs(DX).gt.DXL
       do130 = .true.
       do while (do130)
          !  130 N = 1
          DM = 1._r8
          PNI = 1._r8
          PNJ = X

          !// 140 loops from N=2 to N.eq.NG
          ! 140 N = N+1
          do N = 2, NG
             DN = N
             PNK = ((DM + DN) * X * PNJ - DM * PNI) / DN
             PNI = PNJ
             PNJ = PNK
             DM = DN
          end do !// do N = 2, NG
          !if (N.lt.NG) goto 140

          DX = PNJ * (1._r8 - X * X) / DNG / (PNI -X * PNJ)
          X = X - DX
          !if (abs(DX).gt.DXL) goto 130
          if (abs(DX) .le. DXL) do130 = .false.
       end do

       J = NG + 1 - I
       XP(I) = XMID - XHAF * X
       XP(J) = XMID + XHAF * X
       WT(I) = 2._r8 * XHAF * (1._r8 - X * X) / (DNG * PNI)**2
       WT(J) = WT(I)
    end do !// do I = 1, NN (120-loop)

    !//---------------------------------------------------------------------
  end subroutine GAUSST2
  !//-----------------------------------------------------------------------







  !//-----------------------------------------------------------------------
  subroutine DBLDBL(YGRD,YDGRD,XGRD,XDGRD,YEDG,YDEDG,XEDG,XDEDG, &
                 YGRDW,YDGRDW,XGRDW,XDGRDW,YEDGW,YDEDGW,XEDGW,XDEDGW, &
                 ZDEGI,ZDEGJ,IMAP,JMAP, &
                 IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD,LDEG)
    !//---------------------------------------------------------------------
    !// Define lat-lon of ctm grid and weightings if horizontal resolution
    !// degraded.
    !// For example, take 1x1 met data and run at 2x2 or 2x3
    !//  IDGRD    Number of high-res longitude grid boxes in each low-res box
    !//  JDGRD    Number of high-res latitude ...
    !//  IMAP     Index of high-res longitude grid box for each IMAPN
    !//  ZDEGI    Fraction of low-res box accounted for by high-res box
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: CPI, ZPI180
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  IPARW,JPARW,IPAR,JPAR,IDGRD,JDGRD
    real(r8), intent(in)  ::  YGRDW(JPARW), YDGRDW(JPARW), &
                            YEDGW(JPARW+1), YDEDGW(JPARW+1), &
                            XGRDW(IPARW), XDGRDW(IPARW), &
                            XEDGW(IPARW+1), XDEDGW(IPARW+1)
    real(r8), intent(out) ::  YGRD(JPAR), YDGRD(JPAR), &
                            YEDG(JPAR+1), YDEDG(JPAR+1), &
                            XGRD(IPAR), XDGRD(IPAR), &
                            XEDG(IPAR+1), XDEDG(IPAR+1)
    !// Ouput
    real(r8), intent(out)  :: ZDEGI(IDGRD,IPAR), ZDEGJ(JDGRD,JPAR)
    integer, intent(out) :: IMAP(IDGRD,IPAR), JMAP(JDGRD,JPAR)
    logical, intent(out) :: LDEG

    !// Locals
    integer :: I, J, K, JW
    real(r8) :: Y1P1W(361), Y1P1(181), ZDY
    !//---------------------------------------------------------------------

    if (IDGRD .eq. 1 .and. JDGRD .eq. 1) then
       !// For native resolution no degradation required
       LDEG = .false.
       ZDEGI(:,:) = 1._r8
       ZDEGJ(:,:) = 1._r8

       do I = 1, IPAR
          IMAP(1,I) = I
          XGRD(I)   = XGRDW(I)
          XDGRD(I)  = XDGRDW(I)
          XEDG(I)   = XEDGW(I)
          XDEDG(I)  = XDEDGW(I)
       end do
       XEDG(IPAR+1)  = XEDGW(IPAR+1)
       XDEDG(IPAR+1) = XDEDGW(IPAR+1)

       do J = 1, JPAR
          JMAP(1,J) = J
          YGRD(J)   = YGRDW(J)
          YDGRD(J)  = YDGRDW(J)
          YEDG(J)   = YEDGW(J)
          YDEDG(J)  = YDEDGW(J)
       end do
       YEDG(JPAR+1)  = YEDGW(JPAR+1)
       YDEDG(JPAR+1) = YDEDGW(JPAR+1)

    else
       !// Otherwise, degrade
       LDEG = .true.
       write(6,'(a,2i5)') 'Degrade Longitude resolution:', IPARW,IPAR
       write(6,'(a,2i5)') 'Degrade Latitude  resolution:', JPARW,JPAR

       do I = 1, IPAR
          do K = 1, IDGRD
             IMAP(K,I)  = K + (I-1) * IDGRD
             ZDEGI(K,I) = 1._r8 / real(IDGRD, r8)
          end do
       end do
       do J = 1, JPAR
          do K = 1, JDGRD
             JMAP(K,J)  = K + (J-1) * JDGRD
          end do
       end do

       do I = 1, IPAR
          XEDG(I)  = XEDGW(IMAP(1,I))
          XDEDG(I) = XDEDGW(IMAP(1,I))
       end do
       XEDG(IPAR+1)  = XEDGW(IPARW+1)
       XDEDG(IPAR+1) = XDEDGW(IPARW+1)
       do I = 1, IPAR
          XGRD(I)  = 0.5_r8 * (XEDG(I) + XEDG(I+1))
          XDGRD(I) = ZPI180 * XGRD(I)
       end do

       do J = 1, JPAR
          YEDG(J)  = YEDGW(JMAP(1,J))
          YDEDG(J) = YDEDGW(JMAP(1,J))
       end do
       YEDG(JPAR+1)  = YEDGW(JPARW+1)
       YDEDG(JPAR+1) = YDEDGW(JPARW+1)
       do J = 1, JPAR
          YGRD(J)  = 0.5_r8 * (YEDG(J) + YEDG(J+1))
          YDGRD(J) = ZPI180 * YGRD(J)
       end do

       Y1P1(1) = -1._r8
       do J = 2, JPAR
          Y1P1(J) = sin(YEDG(J))
       end do
       Y1P1(JPAR+1) = 1._r8
       Y1P1W(1) = -1._r8
       do J = 2, JPARW
          Y1P1W(J) = sin(YEDGW(J))
       end do
       Y1P1W(JPARW+1) = 1._r8

       do J = 1, JPAR
          ZDY = 1._r8 / (Y1P1(J+1) - Y1P1(J))
          do K = 1, JDGRD
             JW = JMAP(K,J)
             ZDEGJ(K,J)  = (Y1P1W(JW+1) - Y1P1W(JW)) * ZDY
          end do
       end do

       !// Output mapping
       write(6,'(a)') 'Longitude re-mapping: '
       do I = 1,3
          write(6,1000) &
               I,IDGRD,(IMAP(K,I),ZDEGI(K,I),K=1,IDGRD)
       end do
       do I = IPAR-1, IPAR
          write(6,1000) &
               I,IDGRD,(IMAP(K,I),ZDEGI(K,I),K=1,IDGRD)
       end do
       write(6,'(a)') 'Latitude  re-mapping: '
       do J = 1,3
          write(6,1000) &
               J,JDGRD,(JMAP(K,J),ZDEGJ(K,J),K=1,JDGRD)
       end do
       do J = JPAR-1, JPAR
          write(6,1000) &
               J,JDGRD,(JMAP(K,J),ZDEGJ(K,J),K=1,JDGRD)
       end do

1000 format(2i3,6(i4,' (',f6.4,')'))

       !// Done setting up degraded arrays.
    end if

    !//---------------------------------------------------------------------
  end subroutine DBLDBL
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine AIRSET(DTMET)
    !//---------------------------------------------------------------------
    !// AIRSET sets up the air mass
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: AIR, XYZA, XYZB, XYB
    use cmn_met, only: U, V, P, Q
    use cmn_parameters, only: G0
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: DTMET
    !// Locals
    integer :: I,J,L
    real(r8) :: PSTART(IPAR,JPAR), DAIRDT(IPAR,JPAR), G100, &
         SUMAD0,SUMAW0,AIRWET,AIRDRY,AIRH2O
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='airset'
    !// --------------------------------------------------------------------
    !// UNITS OF AIR MASS AND TRACER = (kg)
    !// Air mass (kg) is given by area (m^2) * pressure thickness (Pa) / g0
    !// P(I,J)    = surf pressure (Pa) averaged in extended zone.
    !//---------------------------------------------------------------------
    !// INITIALIZE
    !//   AIRD() = air mass at beginning of time step
    !//   SUMAD0 = initialized dry-air mass (kg), this quantity is assumed
    !//            to be conserved throughout the run and is used to
    !//            reset AIRD()
    SUMAD0 = 0._r8
    SUMAW0 = 0._r8
    G100   = 100._r8 / G0

    !// Assume that all values relocated to degraded grid (IPAR,JPAR,LPAR)
    !// average  Q(I,J,L) and P(I,J) over extended polar zones.

    DAIRDT(:,:) = 0._r8
    do L = 1, LPAR
      do J = 1, JPAR
        do I = 1, IPAR-1
          DAIRDT(I,J)  = DAIRDT(I,J) + (U(I,J,L) - U(I+1,J,L)) * G100
        end do
        DAIRDT(IPAR,J) = DAIRDT(IPAR,J) + (U(IPAR,J,L) - U(1,J,L)) * G100
      end do
    end do
    do L = 1, LPAR
      do J = 1, JPAR-1
        do I = 1, IPAR
          DAIRDT(I,J) = DAIRDT(I,J) + (V(I,J,L) - V(I,J+1,L)) * G100
        end do
      end do
      do I = 1, IPAR
        DAIRDT(I,JPAR) = DAIRDT(I,JPAR) + V(I,JPAR,L) * G100
      end do
    end do

    do J = 1, JPAR
      do I = 1, IPAR
        PSTART(I,J) = P(I,J) - (DTMET * DAIRDT(I,J)) / XYB(I,J)
        !// pressure calculated by (DTMET*DAIRDT(I,J)- XYA(I,J)) / XYB(I,J)
        !// doesn't count air mass above LPAR+1.
      end do
    end do

    do L = 1, LPAR
      do J = 1, JPAR
        do I = 1, IPAR
          !//AIRWET      = XYZA(I,J,L) + P(I,J) * XYZB(I,J,L)
          AIRWET      = XYZA(I,J,L) + PSTART(I,J) * XYZB(I,J,L)
          AIRH2O      = Q(I,J,L) * AIRWET
          AIR(I,J,L)  = AIRWET - AIRH2O
          SUMAD0      = SUMAD0 + AIR(I,J,L)
          SUMAW0      = SUMAW0 + AIRH2O
        end do
      end do
    end do
    AIRDRY  = SUMAD0
    AIRWET  = SUMAD0 + SUMAW0
    write(6,'(A,1P,2E12.5)') f90file//':'//subr// &
         ' Initialize Air: dry+water(kg):',SUMAD0,SUMAW0

    !//---------------------------------------------------------------------
  end subroutine AIRSET
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SURF_IN(INFIL1, INFIL2, INFIL3, LISLSCP2)
    !//---------------------------------------------------------------------
    !// Read GEIA 1 x 1 degree land fraction datasets, remap to CTM grid
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDGRD, YDGRD, XDEDG, YDEDG, PLAND, AREAXY
    use cmn_sfc, only: landSurfTypeFrac, ZOI, LAI
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) ::  INFIL1,INFIL2,INFIL3
    logical, intent(in)       ::  LISLSCP2

    !// Parameters
    integer, parameter :: I1x1 = 360,  J1x1 = 180
    integer, parameter :: NR8XY=1

    !// Locals
    real(r8) :: XBEDGE(I1x1+1), YBEDGE(J1x1+1), XYBOX(J1x1)
    real(r8) :: ZOIX(I1x1,J1x1,12), LAIX(I1x1,J1x1,12)
    integer :: LS_TYPX(I1x1,J1x1)

    integer :: I, J, L, M, LENF, fnr, ioerr
    integer :: I4DATA(I1x1,J1x1)
    real(r8) :: R8DATA(I1x1,J1x1), R8XY(IPAR,JPAR,NR8XY), RTMP,RTMP2

    character(len=70) :: TITLE
    character(len=100) :: FNAME
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='surf_in'
    !// --------------------------------------------------------------------

    !// Initialize to be able to distinguish between LISLSCP2/ISLSCP88
    landSurfTypeFrac(:,:,:) = 0._r8

    !// Get a free file id
    fnr = get_free_fileid()

    if (.not.LISLSCP2) then

       !// Typical 1x1 grid with (1,1)-box lower-left edges at (180W, 90S)
       !// Area calculation must be consistent with AREAXY
       do I = 1, I1x1+1
          XBEDGE(I) = (real(I, r8) - 181._r8)
       end do
       do J = 1, J1x1+1
          YBEDGE(J) = (real(J, r8) - 91._r8)
       end do
       do J = 1, J1x1
          XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
                      * (sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
       end do

       !// Open and read land-fraction data (in %)
       open(fnr, file=INFIL1, status='old', form='formatted', iostat=ioerr)
       if (ioerr .ne. 0) then
          write(6,'(a)') f90file//':'//subr//': error opening: '//trim(INFIL1)
          stop
       end if
       read (fnr, '(A)') TITLE
         write(6,'(4a)') ' >>read: ',TITLE,' from file:',INFIL1
       read(fnr,'(20i4)', iostat=ioerr) ((I4DATA(I,J),I=1,I1x1),J=J1x1,1,-1)
       if (ioerr .ne. 0) then
          write(6,'(a,i7)') f90file//':'//subr// &
               ': error reading: '//trim(INFIL1),ioerr
          stop
       end if
       close(fnr)

       do J = 1, J1x1
          do I = 1, I1x1
             R8DATA(I,J) = 1.e-2_r8 * real(I4DATA(I,J), r8) * XYBOX(J)
          end do
       end do

       !// Note that this method of calculating Land Fraction is NEW & PRECISE
       !// it averages the areas, not the % fractions!
       call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
            R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)

       do J = 1, JPAR
          do I = 1, IPAR
             PLAND(I,J) = R8XY(I,J,1) / AREAXY(I,J)
          end do
       end do


       !//------------------------------------------------------------------
       !// read ISLSCP-88 Land Surface Data:  
       !// 16 ISLSCP land surface types mapped onto first 9 of the Olson 
       !// abbreviated Land Surface types (10=wetland & 11=urban not used)
       !// LS_TYP = index of dominant land surface type of 1x1 grid (1:9=MXLS)
       !// ZOIX = monthly roughness length (m) for each 1x1 grid
       !// LAIX = monthly leaf area index  (#) for each 1x1 grid
       !// ZOL  = annual  roughness length (m) for each 1x1 grid (type=>ZO)
       !//------------------------------------------------------------------
       open(fnr, FILE=INFIL2, form='formatted', status='old', iostat=ioerr)
       if (ioerr .ne. 0) then
          write(6,'(a)') f90file//':'//subr//': error opening: '//trim(INFIL2)
          stop
       end if
         write (6,'(2a)') 'Reading vegetation types from: ',INFIL2

       read(fnr,'(a)', iostat=ioerr) TITLE
       if (ioerr .ne. 0) then
          write(6,'(a,i7)') f90file//':'//subr// &
               ': error2 reading TITLE: '//trim(INFIL2),ioerr
          stop
       end if
         write(6,'(a)') TITLE

       do M = 1, 14
          read(fnr,*, iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error reading header: '//trim(INFIL2),ioerr,m
             stop
          end if
       end do

       read(fnr,'(36i2)', iostat=ioerr) LS_TYPX
       if (ioerr .ne. 0) then
          write(6,'(a,i7)') f90file//':'//subr// &
               ': error2 reading LS_TYPX: '//trim(INFIL2),ioerr
          stop
       end if

       do M = 1, 12
          read(fnr,*, iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error reading header2: '//trim(INFIL2),ioerr,m
             stop
          end if
          read(fnr,'(18f7.4)', iostat=ioerr) ((ZOIX(I,J,M), I=1,I1x1),J=1,J1x1)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error reading ZOIX: '//trim(INFIL2),ioerr,m
             stop
          end if
       end do

       do M = 1, 12
          read(fnr,*, iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error reading header3: '//trim(INFIL2),ioerr,m
             stop
          end if
          read(fnr,'(18f5.2)', iostat=ioerr) ((LAIX(I,J,M), I=1,I1x1),J=1,J1x1)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error reading LAIX: '//trim(INFIL2),ioerr,m
             stop
          end if
       end do
       close(fnr)


       !//------------------------------------------------------------------
       !// Set up Land-surface type fractions:
       !// ISLSCP88 uses 100% def on 1x1, only first 9 type of Olson veg class
       !// 1=water     2=snow/ice    3=decid fo    4=conif fo     5=agricult
       !// 6=shrub/gr  7=tropc fo    8=tundra      9=desert
       !//------------------------------------------------------------------
       do L = 1, 9
         do J = 1, J1x1
           do I = 1, I1x1
             if (LS_TYPX(I,J) .eq. L) then
               R8DATA(I,J) = XYBOX(J)
             else
               R8DATA(I,J) = 0._r8
             end if
           end do
         end do

         call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
              R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)

         do J = 1, JPAR
           do I = 1, IPAR
             landSurfTypeFrac(L,I,J) = R8XY(I,J,1)/AREAXY(I,J)
           end do
         end do
       end do

       !// Surface roughness (m)
       do M = 1, 12
         do J = 1, J1x1
           do I = 1, I1x1
             R8DATA(I,J) = ZOIX(I,J,M) * XYBOX(J)
           end do
         end do
         call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
              R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
         do J = 1, JPAR
           do I = 1, IPAR
             ZOI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
           end do
         end do
       end do

       !// Leaf Area Index (dimensionless)
       do M = 1, 12
         do J = 1, J1x1
           do I = 1, I1x1
              R8DATA(I,J) = LAIX(I,J,M) * XYBOX(J)
           end do
         end do
         call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
              R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
         do J = 1, JPAR
           do I = 1, IPAR
             LAI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
           end do
         end do
       end do

    else

       !//------------------------------------------------------------------
       !// ISLSCP2 1x1 grid with (1,1)-box lower-left edges at (180W,90N)
       !//------------------------------------------------------------------
       do I = 1, I1x1+1
          XBEDGE(I) = (real(I, r8) - 181._r8)
       end do
       do J = 1, J1x1+1
          !// Make it easier for interpolation: Data are read in reverse below
          !//YBEDGE(J) = (91._r8 - real(J, r8))
          YBEDGE(J) = (real(J, r8) - 91._r8)
       end do
       do J = 1, J1x1
          XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
               * abs(sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
       end do

       !//------------------------------------------------------------------
       !// Open and read ISLSCP2 MODIS land fraction and type data (unit: %)
       !// into
       !// PLAND   = land fraction (m2/m2)
       !// LS_FRAC = land surface type (up to 18, which defined as following)
       !// Note: type 0 represents water, which is used as land fraction
       !//       and stored in PLAND
       !//------------------------------------------------------------------
       !// The following 17 IGBP classes are used.
       !//  0=Water Bodies                        1=Evergreen Needleleaf Forests
       !//  2=Evergreen Broadleaf Forests         3=Deciduous Needleleaf Forests
       !//  4=Deciduous Broadleaf Forests         5=Mixed Forests
       !//  6=Closed Shrublands                   7=Open Shrublands
       !//  8=Woody Savannas                      9=Savannas
       !// 10=Grasslands                         11=Permanent Wetlands
       !// 12=Croplands                          13=Urban and Built-Up
       !// 14=Cropland/Natural Vegetation Mosaic 15=Permanent Snow and Ice
       !// 16=Barren or Sparsely Vegetated       17=Unclassified
       !//------------------------------------------------------------------
       FNAME = INFIL1
       LENF = len(trim(FNAME))-5
       do M = 1, 18
          R8DATA(:,:) = 0._r8
          write(FNAME(LENF:LENF+1),'(i2.2)') M-1
          write(6,'(a)') 'Reading MODIS landcover file:'
          write(6,'(a)') trim(FNAME)
          open(fnr, file=trim(FNAME), status='old', form='formatted', iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error opening: '//trim(FNAME),m,ioerr
             stop
          end if
          do J = J1x1, 1, -1
             read(fnr,*, iostat=ioerr) (R8DATA(I,J),I=1,I1x1)
             if (ioerr .ne. 0) then
                write(6,'(a,3i7)') f90file//':'//subr// &
                     ': error reading R8DATA: '//trim(FNAME),m,J,ioerr
                stop
             end if
          end do
          close(fnr)

          !// Change from % to fraction
          R8DATA(:,:) = R8DATA(:,:) * 1.e-2_r8
          if (M .eq. 1) then
             do J = 1, J1x1
               do I = 1, I1x1
                 R8DATA(I,J) = (1._r8 - R8DATA(I,J)) * XYBOX(J)
               end do
             end do
             call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
                  R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
             do J = 1, JPAR
               do I = 1, IPAR
                 PLAND(I,J) = R8XY(I,J,1) / AREAXY(I,J)
               end do
             end do
          else
             do J = 1, J1x1
               do I = 1, I1x1
                 R8DATA(I,J) = R8DATA(I,J)*XYBOX(J)
               end do
             end do
             call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
                   R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
             do J = 1, JPAR
               do I = 1, IPAR
                 landSurfTypeFrac(M-1,I,J) = R8XY(I,J,1) / AREAXY(I,J)
               end do
             end do
          end if
       end do

       !// Correct PLAND if small differences
       do J = 1, JPAR
         do I = 1, IPAR
           RTMP = real(sum(landSurfTypeFrac(:,I,J)), r8)
           if (RTMP .gt. 1._r8) then
              !// Scale down
              landSurfTypeFrac(:,I,J) = landSurfTypeFrac(:,I,J) / RTMP
              RTMP = 1._r8
           end if
         
           if (PLAND(I,J) .gt. 0._r8) then
              RTMP2 = PLAND(I,J) - RTMP
              if (abs(RTMP2) .lt. 1.e-12_r8) PLAND(I,J) = RTMP
           end if
         end do
       end do


       !// Open and read ISLSCP2 FASIR surface roughness length into
       !// ZOI = monthly roughness length (m) for each CTM grid.
       !// Data range (0.001,8), permanent ice = -77,
       !// no data over land = -88, water = -99
       FNAME = INFIL2
       LENF = len(trim(FNAME)) - 9
       do M = 1, 12
          write(FNAME(LENF:LENF+5),'(i4,i2.2)') 1998,M ! using data in 1998
          write(6,*) 'Reading FASIR surface roughness length file:'
          write(6,*) FNAME

          open(fnr, file=FNAME, status='old', form='formatted',iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error opening: '//trim(FNAME),m,ioerr
             stop
          end if

          do J = J1x1, 1, -1
             read(fnr,*, iostat=ioerr) (R8DATA(I,J),I=1,I1x1)
             if (ioerr .ne. 0) then
                write(6,'(a,3i7)') f90file//':'//subr// &
                     ': error reading R8DATA: '//trim(FNAME),m,J,ioerr
                stop
             end if
          end do
          close(fnr)
          do J = 1, J1x1
            do I = 1, I1x1
              if (R8DATA(I,J) .lt. 0) R8DATA(I,J) = 0._r8 ! flagged point -> zero
              R8DATA(I,J) = R8DATA(I,J) * XYBOX(J)
            end do
          end do
          call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
               R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)

          do J = 1, JPAR
            do I = 1, IPAR
              ZOI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
            end do
          end do
       end do

       !// open and read ISLSCP2 FASIR leaf area index into
       !// LAI = monthly leaf area index for each CTM grid
       !// data range (0,2.497), permanent ice = -77,
       !// no data over land = -88, water = -99
       FNAME = INFIL3
       LENF = len(trim(FNAME)) - 9
       do M = 1, 12
          write(FNAME(LENF:LENF+5),'(i4,i2.2)') 1998,M  ! using data in 1998
          write(6,*) 'Reading FASIR leaf area index file:'
          write(6,*) FNAME

          open(fnr, file=FNAME, status='old', form='formatted',iostat=ioerr)
          if (ioerr .ne. 0) then
             write(6,'(a,2i7)') f90file//':'//subr// &
                  ': error opening: '//trim(FNAME),m,ioerr
             stop
          end if

          do J = J1x1, 1, -1
             read(fnr,*, iostat=ioerr) (R8DATA(I,J),I=1,I1x1)
             if (ioerr .ne. 0) then
                write(6,'(a,3i7)') f90file//':'//subr// &
                     ': error reading R8DATA: '//trim(FNAME),m,J,ioerr
                stop
             end if
          end do
          close(fnr)
          do J = 1, J1x1
            do I = 1, I1x1
              if (R8DATA(I,J) .lt. 0) R8DATA(I,J) = 0._r8 ! flagged point -> zero
              R8DATA(I,J) = R8DATA(I,J) * XYBOX(J)
            end do
          end do
          call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
               R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
          do J = 1, JPAR
            do I = 1, IPAR
              LAI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
            end do
          end do
       end do

    end if
      
    !//---------------------------------------------------------------------
  end subroutine SURF_IN
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_landSurfTypeFrac()
    !// --------------------------------------------------------------------
    !// Reads surface input for CTM3:
    !//   1. Plant function types (landSurfTypeFrac)
    !//
    !// Additional:
    !//   Land fraction
    !//
    !// Based on SURF_IN.
    !//
    !// Amund Sovde Haslerud, November 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, PLAND, AREAXY
    use cmn_sfc, only: landSurfTypeFrac, LANDUSE_IDX, LANDUSE_YEAR, &
         ZOI, LAI, NVGPAR, fileLandSurfTypeFrac
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use ncutils, only: get_netcdf_var_1d, get_netcdf_dim, readnc_3d_from4d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: getYear, nLon, nLat, nPft, nTime, t_idx, I, J, M
    real(r8) :: RTMP, RTMP2
    real(r8), dimension(IPAR, JPAR, 1) :: R8XY
    real(r8), allocatable, dimension(:) :: XBEDGE, YBEDGE, inTime, XYBOX
    real(r8), allocatable, dimension(:,:) :: R8DATA
    real(r8), allocatable, dimension(:,:,:) :: inPFT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_landSurfTypeFrac'
    !// --------------------------------------------------------------------

    !// landSurfTypeFrac
    landSurfTypeFrac(:,:,:) = 0._r8

    !// Input needed:
    !// 1. File name, or path
    !// 2. Year to read (specific or 9999 for metdata)
    if (LANDUSE_YEAR .eq. 9999) then
       getYear = MYEAR
    else
       getYear = LANDUSE_YEAR
    end if

    !// All datasets should have the same format. The original UCI CTM
    !// MODIS files have been converted to netCDF4 format to be easy to read.


    !// Note that the categories differ from LANDUSE_IDX 2 and 3.


    !// Need to find time index for dataset
    if (LANDUSE_IDX .eq. 2) then
       !//------------------------------------------------------------------
       !// ISLSCP2
       !// Only available for year 2001 (generated from Oct 2000 - Oct 2001)
       !//------------------------------------------------------------------
       if (.not. (LANDUSE_YEAR .eq. 9999 .or. LANDUSE_YEAR .eq. 2001)) then
          write(6,'(a)') f90file//':'//subr// &
               ': LANDUSE_YEAR must be 9999 or 2001 for MODIS (LANDUSE_IDX = 2)'
          stop
       end if
       !// Override the year to get - there is only one.
       getYear = 2001
    end if

    !// Read grid edges and time info from file
    call get_netcdf_var_1d(fileLandSurfTypeFrac, 'ilongitude',  XBEDGE)
    call get_netcdf_var_1d(fileLandSurfTypeFrac, 'ilatitude',  YBEDGE)
    call get_netcdf_var_1d(fileLandSurfTypeFrac, 'time', inTime  )

    nLon  = size(XBEDGE) - 1
    nLat  = size(YBEDGE) - 1
    nTime = size(inTime)

    !// Need the last dimension (pft)
    call get_netcdf_dim(fileLandSurfTypeFrac, 'pft', nPft  )

    !// Allocate PFT, area size and array for regridding
    allocate( inPFT(nLon,nLat,nPft), XYBOX(nLat), R8DATA(nLon,nLat) )


    !// Which year to fetch
    do t_idx = 1, nTime
       if (nint(inTime(t_idx)) .eq. getYear) exit
    end do

    if (t_idx .gt. nTime) then
       write(6,'(a)') f90file//':'//subr// &
            ': Time index not found on file '//trim(fileLandSurfTypeFrac)
       stop
    end if
    write(6,'(a,2i5)') f90file//':'//subr// &
         ': getting year and index on file:',nint(inTime(t_idx)),t_idx


    !// Read dataset (unit is %)
    call readnc_3d_from4d(fileLandSurfTypeFrac, 'lon',nLon,'lat',nLat, &
       'pft', nPft, 'time', t_idx, 'PFT_PCT', inPFT)

    !// Change from % to fraction
    inPFT(:,:,:) = inPFT(:,:,:) * 1.e-2_r8


    !// Generate area of grid boxes
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
            * abs(sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
    end do

    if (LANDUSE_IDX .eq. 2) then
       !//------------------------------------------------------------------
       !// ISLSCP2 MODIS land fraction and type data
       !//------------------------------------------------------------------
       !// Input CTM
       !// index index Land cover type
       !//  0    N/A   Water Bodies
       !/   1     1    Evergreen Needleleaf Forests
       !//  2     2    Evergreen Broadleaf Forests
       !//  3     3    Deciduous Needleleaf Forests
       !//  4     4    Deciduous Broadleaf Forests
       !//  5     5    Mixed Forests
       !//  6     6    Closed Shrublands
       !//  7     7    Open Shrublands
       !//  8     8    Woody Savannas
       !//  9     9    Savannas
       !// 10    10    Grasslands
       !// 11    11    Permanent Wetlands
       !// 12    12    Croplands
       !// 13    13    Urban and Built-Up
       !// 14    14    Cropland/Natural Vegetation Mosaic
       !// 15    15    Permanent Snow and Ice
       !// 16    16    Barren or Sparsely Vegetated
       !// 17    17    Unclassified
       !//------------------------------------------------------------------
       if (nPft .gt. NVGPAR) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': nPft / NVGPAR mismatch: ',nPft,NVGPAR
          stop
       end if

       do M = 1, nPft

          do J = 1, nLat
             do I = 1, nLon
                R8DATA(I,J) = inPFT(I,J,M) * XYBOX(J)
             end do
          end do

          call E_GRID(R8DATA, XBEDGE, YBEDGE, nLon, nLat, &
               R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)

          if (M .eq. 1) then
             !// PLAND from ocean coverage
             do J = 1, JPAR
                do I = 1, IPAR
                   PLAND(I,J) = max(0._r8, 1._r8 - R8XY(I,J,1) / AREAXY(I,J) )
                end do
             end do
          else
             !// landSurfTypeFrac
             do J = 1, JPAR
                do I = 1, IPAR
                   landSurfTypeFrac(M-1,I,J) = R8XY(I,J,1) / AREAXY(I,J)
                end do
             end do
          end if
       end do !// do M = 1, nPft

       !// Possible correction to landSurfTypeFrac and PLAND
       do J = 1, JPAR
          do I = 1, IPAR
             RTMP = sum(landSurfTypeFrac(:,I,J))
             if (RTMP .gt. 1._r8) then
                !// Scale down if larger than 1. Unchanged if within [0-1].
                landSurfTypeFrac(:,I,J) = landSurfTypeFrac(:,I,J) / RTMP
                RTMP = 1._r8 !// New total of landSurfTypeFrac
             end if

             !// For MODIS data, there could be some small inconsistencies
             !// in PLAND vs sum of landSurfTypeFrac.
             if (PLAND(I,J) .gt. 0._r8) then
                RTMP2 = PLAND(I,J) - RTMP
                if (abs(RTMP2) .lt. 1.e-12_r8) PLAND(I,J) = RTMP
             end if
          end do
       end do


    else if (LANDUSE_IDX .eq. 3) then
       !//------------------------------------------------------------------
       !// CLM4 PFTs
       !//------------------------------------------------------------------
       !// Note the order of CLM land use types, starting on barren land
       !// and then MEGAN types - where type 2 and 3 follows Gunther et al,
       !// 2012 (G2012, doi:10.5194/gmd-5-1471-2012), not
       !// original MEGAN code. CTM3 follows G2012.
       !//
       !// While the CLM dataset contains two types of Crop, only one
       !// contains data. This makes no difference to MEGAN, which have
       !// same settings for these (Crop2 is Corn in MEGAN).
       !//
       !// Structured into CTM with barren land at the end:
       !// Input CTM
       !// index index Land cover type
       !//  1    17    Barren land
       !//  2     1    Needleaf evergreen temperate tree
       !//  3     2    Needleaf evergreen boreal tree
       !//  4     3    Needleaf deciduous boreal tree
       !//  5     4    Broadleaf evergreen tropical tree
       !//  6     5    Broadleaf evergreen temperate tree
       !//  7     6    Broadleaf deciduous tropical tree
       !//  8     7    Broadleaf deciduous temperate tree
       !//  9     8    Broadleaf deciduous boreal tree
       !// 10     9    Broadleaf evergreen temperate shrub
       !// 11    10    Broadleaf deciduous temperate shrub
       !// 12    11    Broadleaf deciduous boreal shrub
       !// 13    12    Arctic C3 grass (cold)
       !// 14    13    C3 grass (cool)
       !// 15    14    C4 grass (warm)
       !// 16    15    Crop1
       !// 17    16    Crop2
       !//------------------------------------------------------------------
       do M = 1, nPft

          do J = 1, nLat
             do I = 1, nLon
                R8DATA(I,J) = inPFT(I,J,M) * XYBOX(J)
             end do
          end do

          call E_GRID(R8DATA, XBEDGE, YBEDGE, nLon, nLat, &
               R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)

          !// landSurfTypeFrac
          if (M .eq. 1) then
             !// Put barren land at the end
             do J = 1, JPAR
                do I = 1, IPAR
                   landSurfTypeFrac(nPft,I,J) = R8XY(I,J,1) / AREAXY(I,J)
                end do
             end do
          else
             do J = 1, JPAR
                do I = 1, IPAR
                   landSurfTypeFrac(M-1,I,J) = R8XY(I,J,1) / AREAXY(I,J)
                end do
             end do
          end if
       end do !// do M = 1, nPft


       !// Possible correction to landSurfTypeFrac
       do J = 1, JPAR
          do I = 1, IPAR
             RTMP = sum(landSurfTypeFrac(:,I,J))
             if (RTMP .gt. 1._r8) then
                !// Scale down if larger than 1. Unchanged if within [0-1].
                landSurfTypeFrac(:,I,J) = landSurfTypeFrac(:,I,J) / RTMP
             end if
          end do
       end do

       !// PLAND is set from total landSurfTypeFrac
       do J = 1, JPAR
          do I = 1, IPAR
             PLAND(I,J) = sum(landSurfTypeFrac(:,I,J))
          end do
       end do

    else
       write(6,'(a)') f90file//':'//subr// &
            ': No such LANDUS_IDX available',LANDUSE_IDX
       stop
    end if



       
    !// At the end
    deallocate( XBEDGE, YBEDGE, inTime, XYBOX, inPFT, R8DATA )

    write(6,'(a)') f90file//':'//subr//': landSurfTypeFrac is set'


    !// --------------------------------------------------------------------
  end subroutine read_landSurfTypeFrac
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine read_LAI()
    !// --------------------------------------------------------------------
    !// Reads surface input for CTM3:
    !//   1. Monthly leaf area index from ISLSCP2 FASIR adjusted.
    !//
    !// Based on SURF_IN.
    !//
    !// Amund Sovde Haslerud, December 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, PLAND, AREAXY
    use cmn_sfc, only: LAI, LAI_YEAR, fileLAI
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d, &
         get_netcdf_dim, readnc_3d_from4d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: getYear, nLon, nLat, nTime, t_idx, I, J, M, status
    real(r8), dimension(IPAR, JPAR, 1) :: R8XY
    real(r8), allocatable, dimension(:) :: XBEDGE, YBEDGE, inTimeYear, XYBOX
    real(r8), allocatable, dimension(:,:) :: R8DATA
    real(r8), allocatable, dimension(:,:,:) :: inLAI
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_LAI'
    !// --------------------------------------------------------------------

    !// LAI
    LAI(:,:,:) = 0._r8

    !// Input needed:
    !// 1. File name, or path
    !// 2. Year to read (specific or 9999 for metdata)
    if (LAI_YEAR .eq. 9999) then
       getYear = MYEAR
    else
       getYear = LAI_YEAR
    end if

    !// Need to find time index for dataset
    !//---------------------------------------------------------------------
    !// ISLSCP2
    !// Only available for year 1982-1998
    !//---------------------------------------------------------------------
    if (LAI_YEAR .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': LAI read from climatology'
    else if (LAI_YEAR .eq. 9999) then
       write(6,'(a)') f90file//':'//subr// &
            ': LAI will be read from MYEAR'
       if (MYEAR .lt. 1982 .or. MYEAR .gt. 1998) then
          write(6,'(a)') f90file//':'//subr// &
               ': LAI_YEAR must be 0000, 9999 or in the range of 1982-1998.'
          stop
       end if
    end if

    !// Read grid edges and time info from file
    call get_netcdf_var_1d(fileLAI, 'ilongitude',  XBEDGE)
    call get_netcdf_var_1d(fileLAI, 'ilatitude',  YBEDGE)

    nLon  = size(XBEDGE) - 1
    nLat  = size(YBEDGE) - 1

    !// Allocate LAI, area size and array for regridding
    allocate( inLAI(nLon,nLat,12), XYBOX(nLat), R8DATA(nLon,nLat) )


    if (LAI_YEAR .eq. 0) then
       !// Use climatology
       call get_netcdf_var_3d(fileLAI,'CLIM_LAI', inLAI, nlon,nlat,12 )

    else
       !// Use meteorological year or specific year (getYear)
       call get_netcdf_var_1d(fileLAI, 'time_year', inTimeYear  )
       nTime = size(inTimeYear)

       !// Which year to fetch
       do t_idx = 1, nTime
          if (nint(inTimeYear(t_idx)) .eq. getYear) exit
       end do

       if (t_idx .gt. nTime) then
          write(6,'(a)') f90file//':'//subr// &
               ': Time index not found on file '//trim(fileLAI)
          stop
       end if
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': getting year and index on file:',nint(inTimeYEAR(t_idx)),t_idx
       !// Read dataset
       call readnc_3d_from4d(fileLAI, 'lon',nLon,'lat',nLat, &
            'time_month', 12, 'time_year', t_idx, 'LAI', inLAI)

    end if


    !// Generate area of grid boxes
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
            * abs(sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
    end do


    !// Interpolate LAI for each month
    do M = 1, 12
       do J = 1, nLat
          do I = 1, nLon
             ! flagged point -> zero
             R8DATA(I,J) = max(inLAI(I,J,M) * XYBOX(J), 0._r8)
          end do
       end do

       call E_GRID(R8DATA, XBEDGE, YBEDGE, nLon, nLat, &
            R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)

       do J = 1, JPAR
          do I = 1, IPAR
             LAI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
          end do
       end do
    end do !// do M = 1, 12


    !// At the end
    deallocate( XBEDGE, YBEDGE, XYBOX, inLAI, R8DATA)
    if (allocated(inTimeYear)) deallocate( inTimeYear )

    write(6,'(a,2es16.6)') f90file//':'//subr//': LAI is set, min/max: ', &
         minval(LAI),maxval(LAI)

    !// --------------------------------------------------------------------
  end subroutine read_LAI
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_ZOI()
    !// --------------------------------------------------------------------
    !// Reads surface input for CTM3:
    !//   1. Monthly leaf area index from ISLSCP2 FASIR adjusted.
    !//
    !// Based on SURF_IN.
    !//
    !// Amund Sovde Haslerud, December 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, PLAND, AREAXY
    use cmn_sfc, only: ZOI, ZOI_YEAR, fileZOI
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d, &
         get_netcdf_dim, readnc_3d_from4d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: getYear, nLon, nLat, nTime, t_idx, I, J, M
    real(r8), dimension(IPAR, JPAR, 1) :: R8XY
    real(r8), allocatable, dimension(:) :: XBEDGE, YBEDGE, inTimeYear, XYBOX
    real(r8), allocatable, dimension(:,:) :: R8DATA
    real(r8), allocatable, dimension(:,:,:) :: inZOI
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_ZOI'
    !// --------------------------------------------------------------------

    !// ZOI (Z0)
    ZOI(:,:,:) = 0._r8

    !// Input needed:
    !// 1. File name, or path
    !// 2. Year to read (specific or 9999 for metdata)
    if (ZOI_YEAR .eq. 9999) then
       getYear = MYEAR
    else
       getYear = ZOI_YEAR
    end if

    !// Need to find time index for dataset
    !//---------------------------------------------------------------------
    !// ISLSCP2
    !// Only available for year 1982-1998
    !//---------------------------------------------------------------------
    if (ZOI_YEAR .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Z0 read from climatology'
    else if (ZOI_YEAR .eq. 9999) then
       write(6,'(a)') f90file//':'//subr// &
            ': Z0 read from MYEAR'
       if (MYEAR .lt. 1982 .or. MYEAR .gt. 1998) then
          write(6,'(a)') f90file//':'//subr// &
               ': ZOI_YEAR must be 0000, 9999 or in the range of 1982-1998.'
          stop
       end if
    end if

    !// Read grid edges and time info from file
    call get_netcdf_var_1d(fileZOI, 'ilongitude',  XBEDGE)
    call get_netcdf_var_1d(fileZOI, 'ilatitude',  YBEDGE)

    nLon  = size(XBEDGE) - 1
    nLat  = size(YBEDGE) - 1

    !// Allocate ZOI, area size and array for regridding
    allocate( inZOI(nLon,nLat,12), XYBOX(nLat), R8DATA(nLon,nLat) )


    if (ZOI_YEAR .eq. 0) then
       !// Use climatology
       call get_netcdf_var_3d(fileZOI,'CLIM_Z0', inZOI, nlon,nlat,12 )

    else
       !// Use meteorological year or specific year (getYear)
       call get_netcdf_var_1d(fileZOI, 'time_year', inTimeYear  )
       nTime = size(inTimeYear)

       !// Which year to fetch
       do t_idx = 1, nTime
          if (nint(inTimeYear(t_idx)) .eq. getYear) exit
       end do

       if (t_idx .gt. nTime) then
          write(6,'(a)') f90file//':'//subr// &
               ': Time index not found on file '//trim(fileZOI)
          stop
       end if
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': getting year and index on file:',nint(inTimeYEAR(t_idx)),t_idx
       !// Read dataset
       call readnc_3d_from4d(fileZOI, 'lon',nLon,'lat',nLat, &
            'time_month', 12, 'time_year', t_idx, 'Z0', inZOI)

    end if


    !// Generate area of grid boxes
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
            * abs(sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
    end do

    !// Interpolate Z0 for each month
    do M = 1, 12
       do J = 1, nLat
          do I = 1, nLon
             ! flagged point -> zero
             R8DATA(I,J) = max(inZOI(I,J,M) * XYBOX(J), 0._r8)
          end do
       end do

       call E_GRID(R8DATA, XBEDGE, YBEDGE, nLon, nLat, &
            R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)

       do J = 1, JPAR
          do I = 1, IPAR
             ZOI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
          end do
       end do
    end do !// do M = 1, 12


    !// At the end
    deallocate( XBEDGE, YBEDGE, XYBOX, inZOI, R8DATA )
    if (allocated(inTimeYear)) deallocate( inTimeYear )

    write(6,'(a,2es16.6)') f90file//':'//subr//': ZOI is set, min/max: ', &
         minval(ZOI),maxval(ZOI)

    !// --------------------------------------------------------------------
  end subroutine read_ZOI
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine read_ZOI_LAI(fileZOI, fileLAI)
    !// --------------------------------------------------------------------
    !// Reads surface input for CTM3:
    !//   2. Roughness lenght (ZOI)
    !//   3. Leaf area index (LAI)

    !// ZOI
    !// Could it be calculated from PFT characteristics?
    !// 1/10 of canopy height, and 0.0002m for sea, 0.005 for barren.
    !// Not sure how to specify big cities.


    !// LAI

       !// Open and read ISLSCP2 FASIR surface roughness length into
       !// ZOI = monthly roughness length (m) for each CTM grid.
       !// Data range (0.001,8), permanent ice = -77,
       !// no data over land = -88, water = -99
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, AREAXY
    use cmn_sfc, only: ZOI, LAI
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) ::  fileZOI, fileLAI

    !// Parameters
    integer, parameter :: I1x1 = 360,  J1x1 = 180
    integer, parameter :: NR8XY=1

    !// Locals
    real(r8) :: XBEDGE(I1x1+1), YBEDGE(J1x1+1), XYBOX(J1x1)
    real(r8) :: ZOIX(I1x1,J1x1,12), LAIX(I1x1,J1x1,12)
    integer :: LS_TYPX(I1x1,J1x1)

    integer :: I, J, L, M, LENF, fnr, ioerr
    integer :: I4DATA(I1x1,J1x1)
    real(r8) :: R8DATA(I1x1,J1x1), R8XY(IPAR,JPAR,NR8XY), RTMP,RTMP2

    character(len=70) :: TITLE
    character(len=160) :: FNAME

    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='surf_in'
    !// --------------------------------------------------------------------

    fnr = get_free_fileid()
    do I = 1, I1x1+1
       XBEDGE(I) = (real(I, r8) - 181._r8)
    end do
    do J = 1, J1x1+1
       !// Make it easier for interpolation: Data are read in reverse below
       !//YBEDGE(J) = (91._r8 - real(J, r8))
       YBEDGE(J) = (real(J, r8) - 91._r8)
    end do
    do J = 1, J1x1
       XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
            * abs(sin(CPI180 * YBEDGE(J+1)) - sin(CPI180 * YBEDGE(J)))
    end do

    FNAME = fileZOI
    print*,trim(fname)
    LENF = len(trim(FNAME)) - 9
    do M = 1, 12
       write(FNAME(LENF:LENF+5),'(i4,i2.2)') 1998,M ! using data in 1998
       write(6,'(a)') f90file//':'//subr// &
            ': Reading FASIR surface roughness length file: '//trim(FNAME)

       open(fnr, file=FNAME, status='old', form='formatted',iostat=ioerr)
       if (ioerr .ne. 0) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': error opening: '//trim(FNAME),m,ioerr
          stop
       end if

       do J = J1x1, 1, -1
          read(fnr,*, iostat=ioerr) (R8DATA(I,J),I=1,I1x1)
          if (ioerr .ne. 0) then
             write(6,'(a,3i7)') f90file//':'//subr// &
                  ': error reading R8DATA: '//trim(FNAME),m,J,ioerr
             stop
          end if
       end do
       close(fnr)
       do J = 1, J1x1
          do I = 1, I1x1
             if (R8DATA(I,J) .lt. 0) R8DATA(I,J) = 0._r8 ! flagged point -> zero
             R8DATA(I,J) = R8DATA(I,J) * XYBOX(J)
          end do
       end do
       call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
            R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)

       do J = 1, JPAR
          do I = 1, IPAR
             ZOI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
          end do
       end do
    end do

    !// open and read ISLSCP2 FASIR leaf area index into
    !// LAI = monthly leaf area index for each CTM grid
    !// data range (0,2.497), permanent ice = -77,
    !// no data over land = -88, water = -99
    FNAME = fileLAI
    LENF = len(trim(FNAME)) - 9
    do M = 1, 12
       write(FNAME(LENF:LENF+5),'(i4,i2.2)') 1998,M  ! using data in 1998
       write(6,'(a)') f90file//':'//subr// &
            ': Reading FASIR leaf area index file: '//trim(FNAME)

       open(fnr, file=FNAME, status='old', form='formatted',iostat=ioerr)
       if (ioerr .ne. 0) then
          write(6,'(a,2i7)') f90file//':'//subr// &
               ': error opening: '//trim(FNAME),m,ioerr
          stop
       end if

       do J = J1x1, 1, -1
          read(fnr,*, iostat=ioerr) (R8DATA(I,J),I=1,I1x1)
          if (ioerr .ne. 0) then
             write(6,'(a,3i7)') f90file//':'//subr// &
                  ': error reading R8DATA: '//trim(FNAME),m,J,ioerr
             stop
          end if
       end do
       close(fnr)
       do J = 1, J1x1
          do I = 1, I1x1
             if (R8DATA(I,J) .lt. 0) R8DATA(I,J) = 0._r8 ! flagged point -> zero
             R8DATA(I,J) = R8DATA(I,J) * XYBOX(J)
          end do
       end do
       call E_GRID(R8DATA, XBEDGE, YBEDGE, I1x1, J1x1, &
            R8XY, XDEDG, YDEDG, IPAR, JPAR, NR8XY)
       do J = 1, JPAR
          do I = 1, IPAR
             LAI(I,J,M) = R8XY(I,J,1) / AREAXY(I,J)
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine read_ZOI_LAI
  !// ----------------------------------------------------------------------





  !//-----------------------------------------------------------------------
  subroutine gridICAO (YDEDG,XDEDG,PIJL,LCAO,IPAR,JPAR,LPAR,NCAO, &
                       ICAO,JCAO,FNAME)
    !//---------------------------------------------------------------------
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR, JPAR, LPAR, NCAO  ! NCAO=100  stop at 50 kft
    real(r8), intent(in)  :: YDEDG(JPAR+1), XDEDG(IPAR+1), &
                           PIJL(IPAR,JPAR,LPAR+1)
    character(len=80), intent(in) ::  FNAME
    !// Output
    integer, intent(out) :: LCAO(NCAO,IPAR,JPAR), &
                            JCAO(180), ICAO(360)

    !// Local 
    real(r8) :: PCAO(NCAO+1), YCAO(180), XCAO(360)
    real(r8) :: YJ(JPAR+1), XI(IPAR+2), PX(LPAR+1)
    integer :: K,N,L,LL,J,JJ,I,II, fnr, ioerr
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr='gridICAO'
    !// --------------------------------------------------------------------

    !// Get a free file id
    fnr = get_free_fileid()

    !// ICAO emisisons input with gridded indices: IC, JC, LC
    !// are put into CTM grid (I,J,L):
    !//   I = ICAO(IC),  J = JCAO(JC),  L = LCAO(LC,I,J)
    !// these transfer arrays are calculated once with the initial gridding

    !// set up the ICAO grid:
    open(fnr, file=trim(FNAME), status='OLD', iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a,i7)') f90file//':'//subr// &
            ': error opening: '//trim(FNAME)
       stop
    end if
    read(fnr,*)
    read(fnr,'(4x,f10.3)') PCAO
    close (fnr)

    write(6,*) 'finished reading PCAO',PCAO(NCAO+1)
    do J = 1,180
       YCAO(J) = float(J) - 90.5_r8
    end do
    do I = 1,360
       XCAO(I) = float(I) - 0.5_r8
    end do

    !// calculate JCAO(1:180) = J-index for each latitude ICAO box
    !// given YJ(1:JPAR+1) = edge of latitude boxes (deg): -90.0 to +90.0
    do J = 1, JPAR+1
       YJ(J) = YDEDG(J)
    end do
    K = 1
    N = 1
    do J = 1, JPAR
       do while (YCAO(K) .lt. YJ(J+1) .and. N .le. 180)
          JCAO(K) = J
          N = N + 1
          K = min(N, 180)
       end do
    end do
    JCAO(:) = max(1, min(JPAR, JCAO(:)))

    !//do K = 1, 180
    !//   JJ = JCAO(K)
    !//   write(6,'(i4,f9.3,i6,2f9.3)') K,YCAO(K),JJ,YJ(JJ),YJ(JJ+1)
    !//end do

    !// calculate ICAO(1:360) = I-index for each longitude ICAO box
    !// given XI(1:IPAR+1) = edge of longitude boxes (deg) monoton.increasing
    !// really only meant for EC-like fields with first box centered on GMT
    do I = 1, IPAR+1
       XI(I) = XDEDG(I)
    end do
    do I = 2, IPAR+1
       if (XI(I) .lt. 0._r8) then
          XI(I) = XI(I) + 360._r8
       end if
    end do
    XI(IPAR+2) = XI(IPAR+1) + XI(2) - XI(1)
    K = 1
    N = 1
    do I = 1, IPAR+1
       do while (XCAO(K) .lt. XI(I+1) .and. N .le. 360)
          ICAO(K) = I
          if (I .gt. IPAR) then
             ICAO(K) = 1
          end if
          N = N + 1
          K = min(N, 360)
       end do
    end do
    ICAO(:) = max(1, min(IPAR, ICAO(:)))

    !//do K = 1, 360
    !//   II = ICAO(K)
    !//   write(6,'(i4,f9.3,i6,2f9.3)') K,XCAO(K),II,XI(II),XI(II+1)
    !//end do

    !// LCAO(1:110,I,J) is the layer L=1:LPAR to inject ICAO emissions
    !// use annual mean surf pressure to determine CTM layers: PIJL(I,J,L)
    do J = 1, JPAR
       do I = 1, IPAR

          do L = 1, LPAR+1
             PX(L) = PIJL(I,J,L)
          end do
          K = 1
          N = 1
          do L = 2, LPAR+1
             do while (PCAO(K+1) .ge. PX(L) .and. N .le.100)
                LCAO(K,I,J) = L - 1
                N = N + 1
                K = min(N, 100)
             end do
          end do
          LCAO(:,I,J) = max(1, min(LPAR, LCAO(:,I,J)))

          !//do K = 1, NCAO
          !//   LL = LCAO(K,I,J)
          !//   write(6,'(i4,2f9.3,i4,2f9.3)') K,PCAO(K),PCAO(K+1), &
          !//        LL,PX(LL),PX(LL+1)
          !//end do

       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine gridICAO
  !//-----------------------------------------------------------------------

 !// ----------------------------------------------------------------------
  subroutine read_growing_season()
    !// --------------------------------------------------------------------
    !// Reads surface input for CTM3:
    !//   1. Preprocessed daily GDAY -> GDAY_MAP.
    !//   2. Preprocessed GLEN -> GLEN_MAP
    !//
    !// Only called during initialization. For runs exceeding one year
    !// a restart and change of file in *.inp is currently needed!
    !// Based on SURF_IN.
    !//
    !// Stefanie Falk, August 2018
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_ctm, only: JYEAR, JMON, JDAY, LYEAR, JDATE, NRMETD, &
         ZDEGI, ZDEGJ, IMAP, JMAP
    use cmn_size, only: IPARW, JPARW, &
         IPAR, JPAR, IDGRD, JDGRD
    use cmn_sfc, only: GDAY_MAP, GLEN_MAP, fileGSMAP
        use regridding, only: TRUNG8
    use ncutils, only: get_netcdf_var_2d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    logical, parameter :: VERBOSE = .false.
    !// Local parameters
    integer :: totdays = 365
    !// To be read from file
    logical :: fex
    !//---------------------------------------------------------------------
    !// Allocatable arrays - double precision
    real(r8), dimension(:,:), allocatable :: W2D, R8XY
    real(r8), dimension(:,:,:), allocatable :: W3D, R8XYZ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_GROWING_SEASON'
    !// --------------------------------------------------------------------
    if (LYEAR) then
       totdays = 366
    end if
    !// Allocate
    allocate( GDAY_MAP(IPAR,JPAR,totdays) )
    !// Initialize
    GDAY_MAP(:,:,:) = 0._r8
    GLEN_MAP(:,:) = 0._r8
    !// Allocate 2D arrays - native resolution
    allocate( W2D(IPARW,JPARW) )
    !// Allocate 2D arrays - window resolution (IPAR/JPAR)
    allocate( R8XY(IPAR,JPAR) )
    !// Allocate 3D arrays - native resolution
    allocate( W3D(IPARW,JPARW,totdays) )
    !// Allocate 3D arrays - window resolution (IPAR/JPAR)
    allocate( R8XYZ(IPAR,JPAR,totdays) )
    
    !// Check if files exist
    inquire(FILE=trim(fileGSMAP), exist=fex)
    if (.not. fex) then
       write(6,'(a)') f90file//':'//subr// &
            ': No such file: '//trim(fileGSMAP)
       stop
    end if
    !// Growing season parameters GDAY and GLEN
    !// --------------------------------------------------------------------
    call get_netcdf_var_2d(fileGSMAP, 'GLEN',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    GLEN_MAP(:,:) = max(0._r8, R8XY(:,:))  !// Limit to positive values
    if (verbose) call gotData('2da','GLEN')
    
    call get_netcdf_var_3d(fileGSMAP, 'GDAY',W3D, IPARW, JPARW, totdays)
    call TRUNG8(W3D, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, totdays, 1)
    GDAY_MAP(:,:,:) = max(0._r8, R8XYZ(:,:,:))  !// Limit to positive values
    if (verbose) call gotData('3da','GDAY')
       
    !// --------------------------------------------------------------------
  end subroutine read_growing_season
  !// ----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module grid
!//=========================================================================

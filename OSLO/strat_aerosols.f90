!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Stratospheric aerosols.
!//=========================================================================
module strat_aerosols
  !// ----------------------------------------------------------------------
  !// MODULE: strat_aerosols
  !// DESCRIPTION: Routines for setting stratospheric aerosols in the
  !//              Oslo CTM. Aerosols are made by Considine at NASA.
  !//
  !// Contains
  !//   subroutine update_strat_backaer
  !//   subroutine update_ba
  !//   subroutine ctm2_set_partarea
  !//   subroutine ctm2_set_partarea5
  !//   subroutine meri_interpol
  !//
  !// Ole Amund Sovde, February 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: JPAR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
      
  !// Vertical resolution of backaer file
  integer, parameter :: LM_S = 27
  !// Resolution of SAD data
  integer, parameter :: &
       startyear = 1979, &
       endyear = 1999, &
       JMSAD = 16, &  ! # of lats (J=1:80S-70S, J=2:70S-60S,..., J=16:70N-80N)
       LMSAD = 36     ! # of layers  (L=1:5-6km, L=2:6-7km, ..., L=36:40-41km)
  !// SAGE grid: /-75 : 10 : 75/
  real(r8), dimension(JMSAD), parameter :: SAGElats = &
          (/ -75._r8,-65._r8,-55._r8,-45._r8,-35._r8,-25._r8,-15._r8,-5._r8, &
               5._r8, 15._r8, 25._r8, 35._r8, 45._r8, 55._r8, 65._r8, 75._r8 /)

  real(r8) :: &
       CTMsad(JPAR,LM_S), & ! Raw data for SA area densities @ LM_S layers
       SAalt(LM_S), &       ! Heights where SA field is defined
       SApres(LM_S)         ! Std.Pres. where SA field is defined

  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private update_ba, ctm2_set_partarea, ctm2_set_partarea5, meri_interpol
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine update_strat_backaer(NDAY,LHET,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Update straospheric aerosols when heterogeneous chemistry is included.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    logical, intent(in) :: LHET !// =LAER.or.LPSC
    logical, intent(in) :: LNEW_MONTH
    !// --------------------------------------------------------------------

    !// Is there a new dataset for background aerosols?
    if (LNEW_MONTH .and. LHET) call update_ba(NDAY)


    !// Adjust zonal means to current meteorology
    call ctm2_set_partarea(NDAY,LHET)

    !// --------------------------------------------------------------------
  end subroutine update_strat_backaer
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_ba(NDAY)
    !// --------------------------------------------------------------------
    !// Reads SAD zonal distribution and interpolates to CTM meridional
    !// grid (lat, SAD heights). Will be interpolated to height in
    !// subroutine ctm2_set_partarea.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    use cmn_ctm, only: JMON, YDEDG
    use cmn_parameters, only: CPI, CPI180
    use regridding, only: E_GRID_Y
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY

    !// Locals
    !// SAD input data
    real(r8) :: SAGEsad(JMSAD,LMSAD), YDEDGSAD(JMSAD+1),delta
    character(len=3) :: MonthChar

    !// For interpolation
    real(r8) :: r8in(jmsad), r8out(jpar)

    !// Indices etc
    integer :: I,J,L,N,Y,COUNTER,SADYEAR,ifnr,tyear

    !// Text strings specifying year and resolution:
    character(len=4)  :: YEAR
    character(len=80) :: FileName
    !// Logical for existence of input file
    logical :: fnr_ok

    !// Parameters
    logical, parameter ::  LOLDINTERP = .true.
    !// --------------------------------------------------------------------

    !// Which year?
    !// The dataset contains 1979-1999, use the closest year
    SADYEAR = min(max(startyear,MYEAR),endyear)
    write(YEAR,'(I4)') SADYEAR

    !// Read original data
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do
    filename = 'Indata_CTM3/backaer_monthly/sad_1979-1999_new.dat'
    open(ifnr, File=filename, Status='Old', Form='Formatted')

    !// Read intro:
    do N = 1, 162
       read(ifnr,*)
    end do

    !// Read data: square microns per cubic centimeter
    !// Read off pre SADYEAR
    do Y = 1, (SADYEAR - startyear)
       do N = 1, 12
          read(ifnr,'(a3,I5)') MonthChar,tYear
          do L = 1, LMSAD
             read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=1, 6)
             read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=7,12)
             read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=13,JMSAD)
          end do
       end do
    end do
    !// Read until this month
    do N = 1, JMON
       read(ifnr,'(a3,I5)') MonthChar,tYear
       do L = 1, LMSAD
          read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=1, 6)
          read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=7,12)
          read(ifnr,'(6e12.4)') (SAGEsad(J,L),J=13,JMSAD)
       end do
    end do
    write(6,'(A,A3,X,I4)') '* Updating background SAD for ',MonthChar,MYear

    close(ifnr)

    !// Interpolate to latitudes
    !// ---------------------------------------------------------------------
    CTMsad(:,:) = 0._r8
    if (LOLDINTERP) then
       !// Use old linearly interpolation
       do L = 1, LM_S
          call meri_interpol(SAGEsad(:,L),CTMsad(:,L))
       end do
    else
       !// Use CTM3 style interpolation. Produces less continuous gradients
       !// when going from a coarse grid to finer grid.
       delta = 160._r8 / real(JMSAD, r8) !// From 80S to 80N
       do J = 1, JMSAD+1
          YDEDGSAD(J) = real(J-1, r8) * delta - 80._r8
       end do
       !// Extend to poles
       YDEDGSAD(1)       = -90._r8
       YDEDGSAD(JMSAD+1) =  90._r8

       do L = 1, LM_S
          !// Must multiply with area (original field is um2/cm3)
          !// Since the field is a zonal mean, the horizontal extent is 2PI,
          !// and the area of each latitudinal slice is
          !//       A0*A0*CPI180*2*PI *
          !//           (sin(CPI180*YDEDGSAD(J+1)) - sin(CPI180*YDEDGSAD(J)))
          !// A0*A0*CPI180*2*PI is skipped since we divide by it after
          !// interpolation.
          do J = 1, JMSAD
             delta = (sin(CPI180 * YDEDGSAD(J+1)) - sin(CPI180 * YDEDGSAD(J)))
             r8in(J) = SAGEsad(J,L) * delta
          end do
          call e_grid_y(r8in,YDEDGSAD,JMSAD,r8out,YDEDG,JPAR)
          !// Divide by new area
          do J = 1, JPAR
             delta = (sin(CPI180 * YDEDG(J+1)) - sin(CPI180 * YDEDG(J)))
             CTMsad(J,L) = r8out(J) / delta
          end do
       end do
    end if !// if (LOLDINTERP) then

    write(6,'(a,i2)') '* Updated stratospheric aerosols SAD for month',JMON

    !// --------------------------------------------------------------------
  end subroutine update_ba
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ctm2_set_partarea(NDAY, LHET)
    !// --------------------------------------------------------------------
    !// Calculates surface area densities for sulfate aerosols are taken
    !// from SAGE II monthly means. Reads original satellite data.
    !// Routine is a combination of old oc_set_partarea and the interpolation
    !// routines made by Michael Gauss.
    !//
    !// Routine assumes that SAD data are means on geographical latitude
    !// (i.e. not equivalent latitude). Then we must also use geographical
    !// latitude.
    !//
    !// To be called outside parallell loop.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: ETAA, ETAB, JMON
    use cmn_met, only: P
    use cmn_oslo, only: PARTAREA
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    logical, intent(in) :: LHET

    !// For interpolation
    real(r8) :: pres, tmpr8
    !// Indices etc
    integer :: I,J,L,COUNTER
    !// --------------------------------------------------------------------

    if (.not. LHET) then
       !// No heterogeneous chemistry; set SAD to zero.
       PARTAREA(:,:,:) = 0._r8
       return
    end if

    !// Interpolate to heights
    !// ------------------------------------------------------------------
    !// Now, SA surface area field is stored in CTMsad(JPAR,LM_S) at JM CTM
    !// latitudes and LM_S SA levels.
    !// Currently these are z*=5.5,6.5,...,31.5 km (LM_S=27)
    do L = 1, LM_S
       tmpr8  = 4.5_r8 + real(L, r8)
       SAalt(L)  = tmpr8                        !// km
       SApres(L) = exp(-tmpr8/7._r8)*1000._r8   !// hPa
    end do

    !// Assume SAD data are means on geographical latitude (i.e.
    !// not equivalent latitude). Then we must also use geographical
    !// latitude.

    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LPAR

             !// Find pressure
             pres = 0.5_r8 * (ETAA(L  ) + ETAB(L  )*P(I,J) &
                            + ETAA(L+1) + ETAB(L+1)*P(I,J))

             !// Do interpolations in pressure (i.e. mass weighted)
             if (pres.gt.SApres(1) .or. pres.lt.SApres(LM_S)) then
                PartArea(L,I,J) = 0._r8  !// out of range
             else
                if (sum(CTMsad(J,:)) .gt. 0._r8) then
                   Counter = 1
                   do while (SApres(Counter).ge.pres)
                      Counter = Counter + 1
                   end do
                   PartArea(L,I,J) = &
                        ( CTMsad(J,Counter-1) * (SApres(Counter) - pres) &
                          + CTMsad(J,Counter) * (pres - SApres(Counter-1)) &
                        ) / (SApres(Counter) - SApres(Counter-1))
                else
                   !// Nothing to interpolate
                   PartArea(L,I,J) = 0._r8
                end if
             end if

          end do !// do L = 1, LPAR
       end do !// do I = 1, IPAR
    end do !// do J = 1, JPAR


    !// Convert from microns2/cm3 into cm2/cm3
    PartArea(:,:,:) = PartArea(:,:,:)*1.e-08_r8
    write(6,'(A,I4,a3,i2)')'  Background SAD (PartArea) updated at '// &
         'NDAY/MONTH: ',NDAY,' / ',JMON
    write(6,'(A,es12.6,1x,es12.6)')'  Min/max value [um2/cm3] ', &
          minval(PartArea(:,:,:))*1.e8_r8,maxval(PartArea(:,:,:))*1.e8_r8

    !// --------------------------------------------------------------------
  end subroutine ctm2_set_partarea
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ctm2_set_partarea5(NDAY, LHET)
    !// --------------------------------------------------------------------
    !// Calculates surface area densities for sulfate aerosols are taken
    !// from SAGE II monthly means. Reads original satellite data.
    !//
    !// This routine is not correct to use with SAD data from Considine,
    !// since those data are on geographical latitude. When a dataset
    !// on equivalent latitude is available, this would be the routine
    !// to use, after doing appropriate re-gridding of data.
    !//
    !// To be called outside parallell loop.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, LPAR
    use cmn_ctm, only: ETAA, ETAB, JMON, YDEDG
    use cmn_met, only: P, T
    use cmn_oslo, only: PARTAREA
    use physics_oslo, only: NTHE, pvthe, theqlat
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    logical, intent(in) :: LHET

    !// For interpolation
    real(r8) :: pres, theta, eqlat, tmpr8
    integer :: TH, the_low, the_hig

    !// Indices etc
    integer :: I,J,L,COUNTER, Jeq,JSLIM,JNLIM

    !// Limit the use of equivalent latitude polarward of this latitude
    !// (use southern hemisphere for quick calculation)
    real(r8), parameter :: eqlim = -30._r8
    !// ------------------------------------------------------------------

    if (.not. LHET) then
       !// No heterogeneous chemistry; set SAD to zero.
       PARTAREA(:,:,:) = 0._r8
       return
    end if

    !// Interpolate to heights
    !// ------------------------------------------------------------------
    !// Now, SA surface area field is stored in CTMsad(JPAR,LM_S) at JM CTM
    !// latitudes and LM_S SA levels.
    !// Currently these are z*=5.5,6.5,...,31.5 km (LM_S=27)
    do L = 1, LM_S
       tmpr8  = 4.5_r8 + real(L, r8)
       SAalt(L)  = tmpr8                      !// km
       SApres(L) = exp(-tmpr8/7._r8)*1000._r8 !// hPa
    end do

    !// Generate mapping from CTM equivalent latitude to SA latitude:
    !// 1. Find eqlat on theta surfaces
    !// 2. For each grid box interpolate from theta to find eqlat
    !//    For theta 
    !// 3. use YDEDG to find corresponding SA latitude

    if (maxval(theqlat) .eq. 0._r8) then
       print*,'ctm2_set_partarea5: Eq.latitude has not been calculated!'
       stop
    end if

    !// Will use equivalent latitude poleward of |eqlim| degrees
    !// Find J-boxes to define these boxes; initialize first:
    JSLIM = 1
    JNLIM = JPAR
    do J = 1, JPAR
       if (ydedg(J) .lt. eqlim) then
          JSLIM = J            !// Box north of southern limit
          JNLIM = JPAR - J + 1 !// Box south of northern limit
       else
          exit
       end if
    end do

    !// Simplified vertical interpolation:
    !// OsloCTM3: Loop has changed order!
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LPAR

             !// Find pressure
             pres = 0.5_r8 * (ETAA(L  ) + ETAB(L  )*P(I,J) &
                            + ETAA(L+1) + ETAB(L+1)*P(I,J))

             !// Approximately height [km]
             !ZStar = -7._r8*log(pres*1.d-3)

             !// Select which latitude to use; CTM or the
             !// equivalent latitude box corresponding to the CTM latitude
             if (J .ge. JSLIM .and. J .le. JNLIM) then
                !// At low latitudes, use CTM grid
                Jeq = J
             else
                !// Find equivalent latitude from NTHE-levels above/below CTM
                theta = T(i,j,l) * (pres*1.e-3_r8)**(-0.2859_r8)
                if (theta .lt. pvthe(1)) then
                   !// Use eqlat for the lowest theta level
                   theta = pvthe(1)
                   the_low = 1
                   the_hig = 2
                else if (theta .gt. pvthe(NTHE)) then
                   !// Use eqlat for highest theta level
                   theta = pvthe(NTHE)
                   the_low = NTHE-1
                   the_hig = NTHE
                else
                   !// Initialize
                   the_low = 1
                   the_hig = 2
                   do TH = NTHE-1,1,-1
                      if (pvthe(TH) .lt. theta) then
                         the_low = TH     !// Theta surface below CTM
                         the_hig = TH + 1 !// Theta surface above CTM
                         exit             !// exit L-loop
                      end if
                   end do
                end if

                !// The interpolated equivalent latitude, on CTM theta value
                eqlat = (theqlat(I,J,the_hig)*(theta - pvthe(the_low)) &
                     + theqlat(I,J,the_low)*(pvthe(the_hig) - theta) ) &
                     / ( pvthe(the_hig) - pvthe(the_low) )
                !// Find equivalent latitude box corresponding to the CTM
                !// latitude
                do Jeq = 1, JPAR
                   if ( eqlat .ge. ydedg(Jeq) .and. &
                        eqlat .lt. ydedg(Jeq+1) ) exit
                end do
             end if !// if (J .ge. JSLIM .and. J .le. JNLIM) then


             !// Do interpolations in pressure (i.e. mass weighted)
             if (pres.gt.SApres(1) .or. pres.lt.SApres(LM_S)) then
                PartArea(L,I,J) = 0._r8  !// out of range
             else
                if (sum(CTMsad(Jeq,:)) .gt. 0._r8) then
                   Counter=1
                   do while (SApres(Counter).ge.pres)
                      Counter = Counter + 1
                   end do
                   PartArea(L,I,J)= &
                        ( CTMsad(Jeq,Counter-1) * (SApres(Counter) - pres) &
                        + CTMsad(Jeq,Counter) * (pres - SApres(Counter-1)) &
                        ) / (SApres(Counter) - SApres(Counter-1))
                else
                   !// Nothing to interpolate
                   PartArea(L,I,J) = 0._r8
                end if
             end if

          end do !// do L = 1, LPAR
       end do !// do I = 1, IPAR
    end do !// do J = 1, JPAR

    !// Convert from microns2/cm3 into cm2/cm3
    PartArea(:,:,:) = PartArea(:,:,:)*1.e-08_r8
    write(6,'(A,I4,a3,i2)')'  Background SAD (PartArea) updated at '// &
         'NDAY/MONTH: ',NDAY,' / ',JMON
    write(6,'(A,es12.6,1x,es12.6)')'  Min/max value [um2/cm3] ', &
         minval(PartArea(:,:,:))*1.e8_r8,maxval(PartArea(:,:,:))*1.e8_r8

    !// --------------------------------------------------------------------
  end subroutine ctm2_set_partarea5
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine meri_interpol(Field_1,Field_2)
    !// --------------------------------------------------------------------
    !// Interpolate linearly in latitude.
    !// Written by Michael Gauss; updated for CTM2/CTM3.
    !//
    !// Ole Amund Sovde, February 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: YDGRD
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/output
    real(r8), intent(in)  :: Field_1(JMSAD)
    real(r8), intent(out) :: Field_2(JPAR)

    integer :: J_1, J_2, N1, N2
    real(r8) :: Deviation
    !// --------------------------------------------------------------------

    !// Grid definitions?

    !// Initialize Field_2
    Field_2(:) = 0._r8

    !// Interpolation in meridional direction
    !// Field_inter(JMSAD) --> Field_2(JPAR):
    do J_2 = 1, JPAR

       !// In grid 1 find southward neighbour, SAGElats(N1), and
       !// northward neighbour, SAGElats(N2), of Lats_2(J_2)

       !// The following method is not very elegant, but compatible
       !// with many different grid definitions
       Deviation = 1000._r8
       N2 = -1
       do J_1 = 1, JMSAD
          if ( (SAGElats(J_1) - YDGRD(J_2)) .lt. Deviation &
               .and. (SAGElats(J_1) - YDGRD(J_2)) .gt. 0._r8 ) then
             Deviation = SAGElats(J_1) - YDGRD(J_2)
             N2 = J_1
          end if
       end do

       !// If grid 2 goes further north than grid 1:
       if (Deviation .eq. 1000) then
          Field_2(J_2) = Field_1(JMSAD)

       !// If grid 2 goes further south than grid 1:
       else if (N2.eq.1) then
          Field_2(J_2) = Field_1(1)

       else
          if (N2 .lt. 0) then
             print*,'* meri_interpol: N2 negative',N2
             stop
          end if
          N1 = N2 - 1
          !// Calculate Field_2(J_2) as a linear combination
          Field_2(J_2) = &
               Field_1(N1) * (SAGElats(N2) - YDGRD(J_2)) &
                  / (SAGElats(N2) - SAGElats(N1)) &
               + Field_1(N2) * (YDGRD(J_2) - SAGElats(N1)) &
                  / (SAGElats(N2) - SAGElats(N1))

       end if

    end do

    !// --------------------------------------------------------------------
  end subroutine meri_interpol
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module strat_aerosols
!//=========================================================================

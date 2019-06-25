!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Amund Sovde Haslerud, September 2017
!//=========================================================================
!// Utilities for running the CTM, but more specific to Oslo CTM3.
!//=========================================================================
module utilities_oslo
  !// ----------------------------------------------------------------------
  !// Utilities for the Oslo chemistry.
  !// Try to avoid using modules here, this module is supposed to be
  !// compiled before other modules, as they may need the routines here.
  !//
  !// Contains:
  !//   subroutine gotoZC_IJ
  !//   subroutine backfromZC_IJ
  !//   subroutine ZC_MASS2CONC
  !//   subroutine ZC_CONC2MASS
  !//   real(r8) function troe
  !//   real(r8) function rate3B
  !//   subroutine get_chmcycles
  !//   subroutine us76_atmosphere
  !//   real(r8) function h2o_sat
  !//   subroutine source_e90
  !//   subroutine decay_e90
  !//   subroutine init_e90
  !//   subroutine tpause_e90
  !//   subroutine tpauseb_e90
  !//   subroutine tpauseb_o3
  !//   subroutine SZA_PN
  !//   subroutine stringUpCase
  !//   subroutine landfrac2mosaic
  !//   subroutine GROWSEASON
  !//   subroutine MAPPED_GROWSEASON
  !//
  !// Stefanie Falk, March - August 2018
  !//   Moved GROWSEASON to utilities.
  !//   Created mapping of daily preprocessed growseason.
  !// Amund Sovde Haslerud, September 2017
  !//   Some updates to use different tropopause definitions.
  !// Ole Amund Sovde, November 2015, October 2008
  !// 
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'utilities_oslo.f90'
  !// ----------------------------------------------------------------------
  public
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
contains
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  subroutine gotoZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,L_START,L_END, &
       IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
       NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)
    !// --------------------------------------------------------------------
    !// Putting BTT (transported tracers) and XSTT (non-transported tracers)
    !// into ZC_LOCAL, to avoid striding in chemistry integration routine.
    !// ONLY moves the column between L_START and L_END, so that this
    !// routine is applicable to both the troposphere and the stratosphere.
    !//
    !// Can consider only putting tropospheric components into ZC_LOCAL for
    !// tropospheric chemistry (stratospheric chemistry will also need some
    !// tropospheric components), but I doubt it will be any faster.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR, JPAR, LPAR, IDBLK, JDBLK, MP, &
                           NPAR, NOTRPAR, TRACER_ID_MAX
    real(r8), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(r8), intent(in), dimension(LPAR,NOTRPAR,IPAR,JPAR):: XSTT
    integer, intent(in) :: I,J, II,JJ, L_START,L_END
    integer,dimension(IDBLK*JDBLK),intent(in) :: MPBLKIB, MPBLKIE
    integer, intent(in) :: chem_idx(NPAR), Xchem_idx(NOTRPAR)
    !// Output
    real(r8), intent(out), dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL

    !// Locals
    integer :: TRACER_ID, N, L, K
    logical, parameter :: LOCALDEBUG = .false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gotoZC_IJ'
    !// --------------------------------------------------------------------

    !// Initialize (may want to remove it for efficiency)
    ZC_LOCAL(:,:) = 0._r8


    !// Transported species
    do N = 1,NPAR
       TRACER_ID = chem_idx(N)
       !// When tracer list was read, it was checked that TRACER_ID
       !// should be ok.
       do L = L_START, L_END
          ZC_LOCAL(TRACER_ID,L) = BTT(L,N,II,JJ)
          if (LOCALDEBUG) then
             !// Should not be necessary - rather test on BTT earlier.
             if (ZC_LOCAL(TRACER_ID,L).ne.ZC_LOCAL(TRACER_ID,L)) then
                write(6,'(a,5i5)') f90file//':'//subr// &
                     ': Tracer is NAN:',TRACER_ID,L,II,JJ,N,L_START, L_END
                do K = 1, LPAR
                   write(6,'(a,i5,es16.6)') '  L, BTT(L)',K,BTT(K,N,II,JJ)
                end do
                stop
             end if
          end if
       end do
    end do

    !// Non-transported species
    do N = 1,NOTRPAR
       TRACER_ID = Xchem_idx(N)
       !// When tracer list was read, it was checked that TRACER_ID
       !// should be ok.
       do L = L_START, L_END
          ZC_LOCAL(TRACER_ID,L) = XSTT(L,N,I,J)
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine gotoZC_IJ
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine backfromZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,L_START, &
       L_END,IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
       NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)
    !// --------------------------------------------------------------------
    !// Putting ZC_LOCAL back into BTT (transported tracers) and XSTT
    !// (non-transported tracers). See gotoZC_IJ.
    !// ONLY moves the column between L_START and L_END, so that this
    !// routine is applicable to both the troposphere and the stratosphere.
    !//
    !// Ole Amund Sovde, November 2014, October 2008
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPAR, JPAR, LPAR, IDBLK, JDBLK, MP, &
                           NPAR, NOTRPAR, TRACER_ID_MAX
    real(r8), intent(in), dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL
    integer, intent(in) :: I,J, II,JJ, L_START,L_END
    integer,dimension(IDBLK*JDBLK),intent(in) :: MPBLKIB, MPBLKIE
    integer, intent(in) :: chem_idx(NPAR), Xchem_idx(NOTRPAR)
    !// Output
    real(r8), intent(out), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(r8), intent(out), dimension(LPAR,NOTRPAR,IPAR,JPAR):: XSTT

    !// Locals
    integer :: TRACER_ID, N, L
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'backfromZC_IJ'
    !// --------------------------------------------------------------------


    !// Transported species back to BTT (no striding in BTT)
    do N = 1,NPAR
       TRACER_ID = chem_idx(N)
       !// When tracer list was read, it was checked that TRACER_ID
       !// should be ok.
       do L = L_START, L_END
          BTT(L,N,II,JJ) = ZC_LOCAL(TRACER_ID,L)
       end do
    end do

    !// Non-transported species
    do N = 1,NOTRPAR
       TRACER_ID = Xchem_idx(N)
       !// When tracer list was read, it was checked that TRACER_ID
       !// should be ok.
       do L = L_START, L_END
          XSTT(L,N,I,J) = ZC_LOCAL(TRACER_ID,L)
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine backfromZC_IJ
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ZC_MASS2CONC(ZC_LOCAL,L_START,L_END,LPAR, &
          NPAR,NOTRPAR,TRACER_ID_MAX, trsp_idx, Xtrsp_idx, &
          TMASS,XTMASS,DV,AVOGNR)
    !// --------------------------------------------------------------------
    !// Converting ZC_LOCAL from mass to concentration.
    !// ONLY converts the column between L_START and L_END, so that this
    !// routine is applicable to both the troposphere and the stratosphere.
    !// 
    !// The conversion uses molecular weight (MW [g/mol]), volume (DV [m3]),
    !// and Avogadros number (Na):
    !// ZC(molec/cm3) = ZC(kg) * 1/(MW * kg/(1e3g)) * Na / (DV * 1.e6cm3/m3)
    !//               = ZC(kg) * 1e3/MW * Na / DV / 1e6
    !//               = ZC(kg) * 1e-3/MW * Na / DV
    !//
    !// Ole Amund Sovde, October 2008, March 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LPAR, NPAR, NOTRPAR, TRACER_ID_MAX
    real(r8), intent(in) :: DV(LPAR) !// Volume of box
    integer, intent(in) :: L_START, L_END
    integer, intent(in) :: trsp_idx(TRACER_ID_MAX), Xtrsp_idx(TRACER_ID_MAX)
    real(r8), intent(in) :: TMASS(NPAR), XTMASS(NOTRPAR)
    real(r8), intent(in) :: AVOGNR

    !// Input/Output
    real(r8), intent(inout), dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL

    !// Locals
    integer :: TRACER_ID, N, XN, L
    real(r8) :: RDUM
    logical, parameter :: LOCALDEBUG = .false.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'ZC_MASS2CONC'
    !// --------------------------------------------------------------------

    !// Loop through chemical species
    !// While this produce some striding in ZC_LOCAL, it reduces the number
    !// of IF statements. Old routine checked on species for all
    !// tracers and levels.
    do TRACER_ID = 1, TRACER_ID_MAX

       N  = trsp_idx(TRACER_ID)  !// Transported index, if transported
       XN = Xtrsp_idx(TRACER_ID) !// Non-transported index, if not transported

       if (N .gt. 0) then
          !// Transported species
          RDUM = 1.e-3_r8 * AVOGNR / TMASS(N)   !// See above for comments
       else if (XN .gt. 0) then
          !// Non-transported species
          RDUM = 1.e-3_r8 * AVOGNR / XTMASS(XN) !// See above for comments
       else
          !// Species not in use, go to next species
          cycle
       end if


       !// Loop through partial column
       do L = L_START,L_END
          !// kg * 1d-3molec/kg / cm3 = molec/cm3
          ZC_LOCAL(TRACER_ID,L) = ZC_LOCAL(TRACER_ID,L) * RDUM / DV(L)
          if (LOCALDEBUG) then
             if (ZC_LOCAL(TRACER_ID,L) .ne. ZC_LOCAL(TRACER_ID,L)) then
                write(6,'(a,4i5,f6.1,es16.4)') f90file//':'//subr// &
                     ': ', TRACER_ID, L, N, XN, RDUM, DV(L)
                stop
             end if
          end if !// if (LOCALDEBUG) then
 
       end do !// do L = L_START, L_END
    end do !// do TRACER_ID = 1, TRACER_ID_MAX

    !// --------------------------------------------------------------------
  end subroutine ZC_MASS2CONC
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ZC_CONC2MASS(ZC_LOCAL,L_START,L_END,LPAR, &
           NPAR,NOTRPAR,TRACER_ID_MAX, trsp_idx, Xtrsp_idx, &
           TMASS,XTMASS,DV,AVOGNR)
    !// --------------------------------------------------------------------
    !// Converting ZC_LOCAL from concentration to mass.
    !// ONLY converts the column between L_START and L_END, so that this
    !// routine is applicable to both the troposphere and the stratosphere.
    !// 
    !// The conversion uses molecular weight (MW [g/mol]), volume (DV [m3]),
    !// and Avogadros number (Na):
    !// ZC(molec/cm3) = ZC(kg) * 1/(MW * kg/(1e3g) * Na / (DV * 1.e6cm3/m3)
    !//               = ZC(kg) * 1e3/MW * Na / DV / 1e6
    !//               = ZC(kg) * 1e-3/MW * Na / DV
    !//
    !// ZC(kg) = ZC(molec/cm3) * 1e3 * DV * MW / Na
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LPAR, NPAR, NOTRPAR, TRACER_ID_MAX
    real(r8), intent(in) :: DV(LPAR) !// Volume of box
    integer, intent(in) :: L_START,L_END
    integer, intent(in) :: trsp_idx(TRACER_ID_MAX), Xtrsp_idx(TRACER_ID_MAX)
    real(r8), intent(in) :: TMASS(NPAR), XTMASS(NOTRPAR)
    real(r8), intent(in) :: AVOGNR
    !// Input/Output
    real(r8), intent(inout), dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL

    !// Locals
    integer :: TRACER_ID, N, XN, L
    real(r8) :: RDUM
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'ZC_CONC2MASS'
    !// --------------------------------------------------------------------

    !// Loop through chemical species
    !// While this produce some striding in ZC_LOCAL, it reduces the number
    !// of IF statements. Old routine checked on species for all
    !// tracers and levels.
    do TRACER_ID = 1, TRACER_ID_MAX

       N  = trsp_idx(TRACER_ID)  !// Transported index, if transported
       XN = Xtrsp_idx(TRACER_ID) !// Non-transported index, if not transported

       if (N .gt. 0) then
          !// Transported species
          RDUM = 1000._r8 * TMASS(N) / AVOGNR   !// See above for comments
       else if (XN .gt. 0) then
          !// Non-transported species
          RDUM = 1000._r8 * XTMASS(XN) / AVOGNR !// See above for comments
       else
          !// Species not in use, go to next species
          cycle
       end if

       !// Loop through partial column
       do L = L_START, L_END

          !// kg * 1d-3molec/kg / cm3 = molec/cm3
          ZC_LOCAL(TRACER_ID,L) = ZC_LOCAL(TRACER_ID,L) * RDUM * DV(L)

       end do !// do L = L_START, L_END

    end do !// do TRACER_ID = 1, TRACER_ID_MAX

    !// --------------------------------------------------------------------
  end subroutine ZC_CONC2MASS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function TROE (KZERO,KINF,FC)
    !// --------------------------------------------------------------------
    !// Calculates ratecoeff. for three body reactions according to
    !// expressions developed by troe. see f.ex. atkinson et al. 1989
    !//
    !// This is the old tropospheric calculation of 3-body reactions,
    !// and is included for historical reasons. It is replaced by
    !// rate3B, which is more flexible.
    !// Amund Sovde, December 2014
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O parameters
    real(r8), intent(in) :: KZERO,KINF,FC
    !// Local parameters
    real(r8) :: RHLP1,RHLP2,RHLP3
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'TROE'
    !// --------------------------------------------------------------------
    write(6,'(a)') f90file//':'//subr//': Not to be used! Use rate3B instead'
    stop
    RHLP1 = KZERO / KINF
    RHLP2 = KZERO * KINF / (KZERO + KINF)
    RHLP3 = 1._r8 / (1._r8 + log10(RHLP1) * log10(RHLP1))
    TROE  = RHLP2 * (FC**RHLP3)
    !// --------------------------------------------------------------------
  end function TROE
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function rate3B(INUM, TZ300,AIRDENS,kzero300,Nzero, &
                         kinf300,Minf, FC,ACT)
    !// --------------------------------------------------------------------
    !// Three-body reaction rates. This is very much the same routine
    !// as the stratospheric code BODY3_BOX, but written a bit more
    !// cleanly.
    !// INUM:    Reaction number for debug use
    !// TZ300:   T/300
    !// AIRDENS: molec/cm3
    !// kzero300: k_o^300 in JPL
    !// kzero300: k_o^300 in JPL
    !// Nzero:    n in JPL
    !// Minf:     m in JPL
    !// FC:     constant, usually 0.6d0
    !// ACT:    0: Regular calculation (k_f)
    !//         1: Calculated as an activation channel (k_f^ca)
    !//
    !// Ole Amund Sovde, December 2014
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: &
         TZ300, AIRDENS, kzero300, Nzero, kinf300, Minf, FC
    integer, intent(in) :: INUM, ACT
    !// Local variables
    real(r8) :: kzero, kinf, B, D, RX
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'rate3B'
    !// --------------------------------------------------------------------

    if (kinf300 .eq. 0._r8) then
       !// No high-pressure limit is given; reaction should
       !// only apply for low pressure.
       RX = kzero300 * TZ300**(-Nzero)
    else
       !// Assume both low and high-pressure limit is given
       kzero = kzero300 * TZ300**(-Nzero)
       kinf  = kinf300  * TZ300**(-Minf)
       B = kzero * AIRDENS / kinf
       D = 1._r8 / ( 1._r8 + log10(B) * log10(B) )
       RX = kzero / ( 1._r8 + B ) * (FC**D)
    end if

    if (ACT .eq. 0) then
       !// Standard 3-body reaction; multiply by air density.
       RATE3B = RX * AIRDENS
    else if (ACT .eq. 1) then
       !// Activation channel; skip multiplication by AIRDENS
       !// Check the JPL reactions if you are in doubt.
       RATE3B = RX
    else
       write(6,'(a, i7)') f90file//':'//subr// &
            ': ACT must be 0 or 1! INUM:', INUM
       stop
    end if

    !// --------------------------------------------------------------------
  end function rate3B
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_chmcycles(CHMCYCLES, NRCHEM, NROPSM)
    !// --------------------------------------------------------------------
    !// Internal UiO loop for chemistry/boundary layer mixing
    !// Number of cycles per NOPS when NROPSM==3.
    !//
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NRCHEM, NROPSM
    !// Output
    integer, intent(out) :: CHMCYCLES
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_chmcycles'
    !// --------------------------------------------------------------------

    if (NRCHEM .eq. 1) then
       CHMCYCLES = 4
       !// Simple test to take NROPSM into account,
       if (NROPSM .ge. 6 .and. NROPSM .lt. 12) CHMCYCLES = 2
       if (NROPSM .ge. 12) CHMCYCLES = 1
    else if (NRCHEM .eq. 2) then
       CHMCYCLES = 2
       !// Simple test to take NROPSM into account
       if (NROPSM .ge. 6) CHMCYCLES = 1
    else
       CHMCYCLES = 1
    end if

    !// --------------------------------------------------------------------
  end subroutine get_chmcycles
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine US76_Atmosphere(alt, sigma, delta, theta)
    !// --------------------------------------------------------------------
    ! PURPOSE: Compute the properties of the 1976 standard atmosphere to 86 km.
    ! AUTHOR: Ralph Carmichael, Public Domain Aeronautical Software
    ! NOTE: If alt > 86, the values returned will not be correct, but they will
    !   not be too far removed from the correct values for density.
    !   The reference document does not use the terms pressure and temperature
    !   above 86 km.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// A R G U M E N T S
    !// --------------------------------------------------------------------
    REAL(r8),INTENT(IN)::  alt    ! geometric altitude, km.
    REAL(r8),INTENT(OUT):: sigma  ! density/sea-level standard density
    REAL(r8),INTENT(OUT):: delta  ! pressure/sea-level standard pressure
    REAL(r8),INTENT(OUT):: theta  ! temperature/sea-level standard temperature
    !// --------------------------------------------------------------------
    !// L O C A L   C O N S T A N T S
    !// --------------------------------------------------------------------
    REAL(r8),PARAMETER:: REARTH = 6369.0_r8 ! radius of the Earth (km)
    REAL(r8),PARAMETER:: GMR = 34.163195_r8 ! hydrostatic constant
    INTEGER,PARAMETER:: NTAB = 8 ! number of entries in the defining tables
    !// --------------------------------------------------------------------
    !// L O C A L   V A R I A B L E S
    !// --------------------------------------------------------------------
    INTEGER :: i,j,k         ! counters
    REAL(r8) :: h            ! geopotential altitude (km)
    REAL(r8) :: tgrad, tbase ! temperature gradient and base temp of this layer
    REAL(r8) :: tlocal       ! local temperature
    REAL(r8) :: deltah       ! height above base of this layer
    !// --------------------------------------------------------------------
    !//  L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )
    !// --------------------------------------------------------------------
    REAL(r8),DIMENSION(NTAB),PARAMETER:: htab = &
         (/   0.0_r8,  11.0_r8,   20.0_r8,   32.0_r8,   47.0_r8, &
             51.0_r8,  71.0_r8,   84.852_r8/)
    REAL(r8),DIMENSION(NTAB),PARAMETER:: ttab = &
         (/288.15_r8, 216.65_r8, 216.65_r8, 228.65_r8, 270.65_r8, &
           270.65_r8, 214.65_r8, 186.946_r8/)
    REAL(r8),DIMENSION(NTAB),PARAMETER:: ptab = &
         (/1.0_r8, 2.233611e-1_r8, 5.403295e-2_r8, 8.5666784e-3_r8, &
         1.0945601e-3_r8, 6.6063531e-4_r8, 3.9046834e-5_r8, 3.68501e-6_r8/)
    REAL(r8),DIMENSION(NTAB),PARAMETER:: gtab = &
         (/-6.5_r8, 0.0_r8, 1.0_r8, 2.8_r8, 0.0_r8, -2.8_r8, -2.0_r8, 0.0_r8/)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'US76_Atmosphere'
    !// --------------------------------------------------------------------

    ! convert geometric to geopotential altitude
    h = alt * REARTH / (alt + REARTH)

    i = 1
    j = NTAB                    ! setting up for binary search
    DO
       k = (i + j) / 2              ! integer division
       IF (h .lt. htab(k)) THEN
          j = k
       ELSE
          i = k
       END IF
       IF (j .le. i+1) EXIT
    END DO

    tgrad = gtab(i)             ! i will be in 1...NTAB-1
    tbase = ttab(i)
    deltah = h - htab(i)
    tlocal = tbase + tgrad * deltah
    theta = tlocal / ttab(1)      ! temperature ratio

    IF (tgrad .eq. 0._r8) THEN    ! pressure ratio
       delta = ptab(i) * EXP(-GMR * deltah / tbase)
    ELSE
       delta = ptab(i) * (tbase / tlocal)**(GMR / tgrad)
    END IF

    sigma = delta / theta         ! density ratio

    !// --------------------------------------------------------------------
  end subroutine US76_Atmosphere
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function h2o_sat(T)
    !// --------------------------------------------------------------------
    !// Returns saturated value of H2O [molec/cm3] at temperature T.
    !// For T<273.15K: From Marti and Mauersberger, GRL Vol 20, No 5,
    !// 363-366, 1993.
    !// --------------------------------------------------------------------
    use cmn_parameters, only: KBOLTZ, AVOGNR, TK_0C
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: T
    !// Locals
    real(r8) :: AMP0, E, TC
    !// --------------------------------------------------------------------
    AMP0 = KBOLTZ * T ! mb * cm3/molec
    if (T .lt. TK_0C) then
       H2O_SAT = 10._r8**(12.537_r8 - 2663.5_r8/T) / (AMP0 * 100._r8)
    else
       TC = T - TK_0C
       E  = 6.1_r8 * 10._r8**(7.45_r8 * TC / (235._r8 + TC))
       H2O_SAT = 217.e-6_r8 * E / T / 18._r8 * AVOGNR
    end if
    !// --------------------------------------------------------------------
  end function h2o_sat
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine source_e90(BTT,DTSRCE,MP)
    !// --------------------------------------------------------------------
    !// Adds source to E90 tracer.
    !//
    !// Ole Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLK, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, &
         AREAXY
    use cmn_chem, only: Ne90
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/Output
    integer, intent(in) :: MP
    real(r8),  intent(inout) ::  BTT(LPAR,NPAR,IDBLK,JDBLK)
    real(r8),  intent(in) :: DTSRCE

    !// Locals
    real(r8) :: ZTOTAREA    !// Total area
    real(r8) :: ADDTC       !// Emitted mass
    integer :: I,II,J,JJ

    !// Global emission rate kg/s
    real(r8), parameter :: g_rate = 1.0605e12_r8/31536000_r8
    !// --------------------------------------------------------------------

    !// Check for E90-tracer
    if (Ne90 .le. 0) return

    !// Increase source at surface
    ZTOTAREA = 1._r8 / sum(AREAXY)

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ   = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II   = I - MPBLKIB(MP) + 1
          ADDTC = g_rate * AREAXY(I,J) * ZTOTAREA * DTSRCE
          BTT(1,Ne90,II,JJ) = BTT(1,Ne90,II,JJ) + ADDTC
       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine source_e90
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine decay_e90 (BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,DTCHEM,MP)
    !// --------------------------------------------------------------------
    !// Provides a simple e-fold decay of species throughout the model
    !// domain. Applies for the E90 tracer, i.e. uses 90 days e-fold decay.
    !// Take from UCI routine DECAY.
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: AREAXY
    use cmn_chem, only: Ne90
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: BTT
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    real(r8), intent(in)::  DTCHEM
    integer, intent(in):: MP

    !// Locals
    integer :: II,JJ,L
    real(r8)  :: F1L
    !// e-fold decay over DECTIM days
    real(r8), parameter :: DECTIM = 90._r8
      !// ------------------------------------------------------------------

    !// Check for E90-tracer
    if (Ne90 .le. 0) return

    if (DECTIM .gt. 0._r8) then
       F1L = exp(-DTCHEM / (86400._r8 * DECTIM))
       do JJ = 1, JDBLK
          do II = 1, IDBLK
             do L = 1, LPAR
                BTT(L,Ne90,II,JJ) = BTT(L,Ne90,II,JJ) * F1L
                BZT(L,Ne90,II,JJ) = BZT(L,Ne90,II,JJ) * F1L
                BZZ(L,Ne90,II,JJ) = BZZ(L,Ne90,II,JJ) * F1L
                BXZ(L,Ne90,II,JJ) = BXZ(L,Ne90,II,JJ) * F1L
                BYZ(L,Ne90,II,JJ) = BYZ(L,Ne90,II,JJ) * F1L
                BXT(L,Ne90,II,JJ) = BXT(L,Ne90,II,JJ) * F1L
                BXX(L,Ne90,II,JJ) = BXX(L,Ne90,II,JJ) * F1L
                BYT(L,Ne90,II,JJ) = BYT(L,Ne90,II,JJ) * F1L
                BYY(L,Ne90,II,JJ) = BYY(L,Ne90,II,JJ) * F1L
                BXY(L,Ne90,II,JJ) = BXY(L,Ne90,II,JJ) * F1L
             end do
          end do
       end do
    end if

    !// --------------------------------------------------------------------
  end subroutine decay_e90
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine init_e90()
    !// --------------------------------------------------------------------
    !// Initialize e90. Set constant in troposphere, assume e-fold above.
    !// If E90 needs to be initialised, LMTROP is first set from PVU.
    !// This is because LMTROP should be calculated after E90 tropopause
    !// is calculated (in physics_oslo.f90), which would fail if E90 was not
    !// set.
    !//
    !// Amund Sovde Haslerud, September 2017
    !//   Added possible initialisation of LMTROP.
    !// Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: STT, AIR
    use cmn_chem, only: Ne90, TMASSMIX2MOLMIX, TMOLMIX2MASSMIX
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: I, J, L
    real(r8) :: rtmp
    real(r8), parameter :: init_ppbv = 110.e-9_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'init_e90'
    !// --------------------------------------------------------------------

    !// Check for E90-tracer
    if (Ne90 .le. 0) return

    !// Is perhaps E90 already set?
    rtmp = maxval(STT(:,:,:,Ne90) / AIR(:,:,:)) * TMASSMIX2MOLMIX(Ne90)
    if (rtmp .gt. 1.e-12_r8) return

    if (maxval(LMTROP) .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
           ': Why is notLMTROP initialised from PVU before E90 '// &
           'initialisation?'
       stop
    end if


    do J = 1, JPAR
       do I = 1, IPAR
          !// Set 110 ppb up to LMTROP-2
          do L = 1, LMTROP(I,J)-2
             STT(I,J,L,Ne90) = init_ppbv * AIR(I,J,L) * TMOLMIX2MASSMIX(NE90)
          end do
          !// Exponential above
          do L = LMTROP(I,J)-1, LPAR
             STT(I,J,L,Ne90) = init_ppbv * AIR(I,J,L) * TMOLMIX2MASSMIX(NE90) &
                  * exp(-real(L-LMTROP(I,J), r8) / real(LPAR-LMTROP(I,J), r8))
          end do
       end do
    end do

    write(6,'(a)') f90file//':'//subr//': Initialized E90 tracer'

    !// --------------------------------------------------------------------
  end subroutine init_e90
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tpause_e90()
    !// --------------------------------------------------------------------
    !// Initialize LSTRATAIS_E90, a 3D logical array where stratopsheric
    !// air is T and tropospheric air is F, based on the E90 tracer.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Changed names.
    !// Ole Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: STT, AIR
    use cmn_chem, only: Ne90, TMOLMIX2MASSMIX, &
         E90VMR_TP, LSTRATAIR_E90, LPAUZTOP
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: I,J,L
    real(r8) :: tpval
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tpause_e90'
    !// --------------------------------------------------------------------

    !// Check for E90-tracer
    if (Ne90 .le. 0) return

    !// E90 vmr value for tropopause, converted to mmr
    tpval = E90VMR_TP * TMOLMIX2MASSMIX(Ne90)

    !// Initialize LSTRATAIR_E90
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             LSTRATAIR_E90(L,I,J) = STT(I,J,L,Ne90) .lt. tpval*AIR(I,J,L)
          end do
       end do
    end do

    !// Fint max height level of tropospheric air
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LPAR
             if (.not.LSTRATAIR_E90(L,I,J)) LPAUZTOP(I,J) = L
          end do
       end do
    end do

    write(6,'(a)') f90file//':'//subr//': Initialized LSTRATAIR_E90'

    !// --------------------------------------------------------------------
  end subroutine tpause_e90
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tpauseb_e90(BTT,AIRB,LSTRATAIR_E90B,MP)
    !// --------------------------------------------------------------------
    !// Define the stratosphere based on E90 tracer.
    !// LSTRATAIR_E90(L,I,J) is a 3D logical, .true. = stratospheric
    !// called within the IJ-blocks.
    !// MPSPLIT splits it into private array LSTRATAIR_E90B.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Changed names.
    !// Ole Amund Sovde, February 2012
    !// --..----------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: Ne90, TMOLMIX2MASSMIX, E90VMR_TP, LPAUZTOP
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/output
    real(r8), intent(in) ::  BTT(LPAR,NPAR,IDBLK,JDBLK)
    real(r8), intent(in) ::  AIRB(LPAR,IDBLK,JDBLK)
    integer, intent(in) ::  MP
    logical, intent(out) :: LSTRATAIR_E90B(LPAR,IDBLK,JDBLK)

    !// Locals
    real(r8) :: tpval
    integer :: I, J, L, II, JJ
    !// ..------------------------------------------------------------------
    character(len=*), parameter :: subr = 'toauseb_e90'
    !// --------------------------------------------------------------------

    !// Check for E90-tracer
    if (Ne90 .le. 0) return

    !// E90 vmr value for tropopause, converted to mmr
    tpval = E90VMR_TP * TMOLMIX2MASSMIX(Ne90)

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = 1, LPAR
             LSTRATAIR_E90B(L,II,JJ) = &
                  BTT(L,Ne90,II,JJ) .lt. tpval * AIRB(L,II,JJ)
             !// Fint max height level of tropospheric air
             if (.not.LSTRATAIR_E90B(L,II,JJ)) LPAUZTOP(I,J) = L
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine tpauseb_e90
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tpauseb_o3(BTT,AIRB,LSTRATAIR_O3,MP)
    !// --------------------------------------------------------------------
    !// Define the stratosphere based on O3 isopleths.
    !// Note that LSTRATAIR_O3(L,I,J) is a 3D logical, .true. = stratospheric
    !// called within the IJ-blocks.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Changed names.
    !// Ole Amund Sovde, February 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TMOLMIX2MASSMIX, O3iso1, O3iso2
    use cmn_met, only: ZOFLE
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/output
    real(r8), intent(in) ::  BTT(LPAR,NPAR,IDBLK,JDBLK)
    real(r8), intent(in) ::  AIRB(LPAR,IDBLK,JDBLK)
    integer, intent(in) ::  MP
    logical, intent(out) :: LSTRATAIR_O3(LPAR,IPAR,JPAR,2)

    !// Locals
    real(r8) :: tpval
    integer :: I, J, L, II, JJ
    !// ------------------------------------------------------------------
    character(len=*), parameter :: subr = 'toauseb_o3'
    !// --------------------------------------------------------------------

    !// Check for O3
    if (trsp_idx(1) .le. 0) return

    !// Select O3 isopleth 1 as chemical tropopause (mmr)
    tpval = O3iso1 * TMOLMIX2MASSMIX(trsp_idx(1))

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = 1, LPAR
             LSTRATAIR_O3(L,I,J,1) = &
                  (ZOFLE(L,I,J) .gt. 4000._r8) .and. &
                   BTT(L,trsp_idx(1),II,JJ) .ge. tpval*AIRB(L,II,JJ)
          end do
          !// Loop down from top and set possible strat boxes
          do L = LPAR, 1, -1
             if (LSTRATAIR_O3(L,I,J,1)) then
                exit
             else
                LSTRATAIR_O3(L,I,J,1) = .true.
             end if
          end do
          !// Repeat to set whole tropopause to false
          !// Will mess up tropopause folds; limit to tropics!
          !if (abs(YDGRD(J)).lt. 30._r8) then
          !   do L = LPAR, 1, -1
          !      if (.not.LSTRATAIR_O3(L,I,J,1)) exit
          !   end do
          !   LSTRATAIR_O3(1:L,I,J,1) = .false.
          !end if

       end do
    end do

    !// Select O3 isopleth 2 as chemical tropopause 
    tpval = O3iso2 * TMOLMIX2MASSMIX(trsp_idx(1))

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = 1, LPAR
             LSTRATAIR_O3(L,I,J,2) = &
                  (ZOFLE(L,I,J) .gt. 4000._r8) .and. &
                  BTT(L,trsp_idx(1),II,JJ) .ge. tpval*AIRB(L,II,JJ)
          end do
          !// Loop down from top and set possible strat boxes
          do L = LPAR, 1, -1
             if (LSTRATAIR_O3(L,I,J,2)) then
                exit
             else
                LSTRATAIR_O3(L,I,J,2) = .true.
             end if
          end do
          !// Repeat to set whole tropopause to false
          !// Will mess up tropopause folds; limit to tropics!
          !if (abs(YDGRD(J)).lt. 30._r8) then
          !   do L = LPAR, 1, -1
          !      if (.not.LSTRATAIR_O3(L,I,J,2)) exit
          !   end do
          !   LSTRATAIR_O3(1:L,I,J,2) = .false.
          !end if

       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine tpauseb_o3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine SZA_PN(GMTIME,NDAY,YGRDJ,XGRDI, SZA,PN)
    !// --------------------------------------------------------------------
    !// Calculate solar zenith angle (SZA) and whether it is polar
    !// night (PN) or not.
    !//
    !// Based on Version FJX6.6.
    !// GMTIME = UT for when J-values are wanted 
    !//   (for implicit solver this is at the end of the time step)
    !// NDAY   = integer day of the year (used for solar lat and declin)
    !// YGRDJ  = laitude (radians) for grid (I,J)
    !// XGDRI  = longitude (radians) for grid (I,J)
    !//
    !// SZA = solar zenith angle in degrees
    !// COSSZA = U0 = cos(SZA)
    !// PN = polar night (1) not polar night (0)
    !// --------------------------------------------------------------------
    use cmn_parameters, only: CPI
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) ::   GMTIME,YGRDJ,XGRDI
    integer, intent(in) ::  NDAY
    !// Output
    real(r8), intent(out) ::  SZA
    integer, intent(out)::  PN
    !// Locals
    real(r8) :: PI, PI180, LOCT, minSZA
    real(r8) :: SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT,COSZ,COSSZA
    integer :: I
    !// --------------------------------------------------------------------

    PI     = CPI
    PI180  = PI/180._r8
    SINDEC = 0.3978_r8 * sin(0.9863_r8 * (real(NDAY, r8) - 80._r8) * PI180)
    SOLDEK = asin(SINDEC)
    COSDEC = cos(SOLDEK)
    SINLAT = sin(YGRDJ)
    SOLLAT = asin(SINLAT)
    COSLAT = cos(SOLLAT)

    !// Find if we have polar night:
    !// Loop through 24 hours and check for minimum SZA
    minSZA = 100._r8
    do I = 1, 24
       LOCT   = ((real(I, r8) * 15._r8) - 180._r8) * PI180 + XGRDI
       COSSZA = COSDEC * COSLAT * cos(LOCT) + SINDEC * SINLAT
       minSZA    = min(minSZA, acos(COSSZA) / PI180)
    end do
    if (minSZA .ge. 90._r8) then
       PN = 1
    else
       PN = 0
    end if

    !// Find SZA
    LOCT   = (((GMTIME) * 15._r8) - 180._r8) * PI180 + XGRDI
    COSSZA = COSDEC * COSLAT * cos(LOCT) + SINDEC * SINLAT
    SZA    = acos(COSSZA) / PI180

    !// --------------------------------------------------------------------
  end subroutine SZA_PN
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine stringUpCase(UCSTR,INSTR)
    !// --------------------------------------------------------------------
    !// Put a string into upper case letters (English alphabet).
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in) :: INSTR
    !// Output
    character(len=*),intent(out) :: UCSTR

    !// Locals
    character(len=*),parameter :: LC='abcdefghijklmnopqrstuvwxyz'
    character(len=*),parameter :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i, n,strlen
    !// --------------------------------------------------------------------

    UCSTR = INSTR
    strlen = min(len(INSTR),len(UCSTR))

    do i = 1, strlen
       !// Find location of letter in lower case constant string
       n = index( LC, INSTR(i:i))
       !// If current substring is a lower case letter, make it upper case,
       !// otherwise keep the input letter.
       if ( n .ne. 0 ) then
          UCSTR(i:i) = UC(n:n)
       else
          UCSTR(i:i) = INSTR(i:i)
       end if
    end do

    !// --------------------------------------------------------------------
  end subroutine stringUpCase
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine landfrac2mosaic(LFmosaic,nLFmosaic,LFin,nLFin, lat, LANDUSE_IDX)
    !// --------------------------------------------------------------------
    !// Convert land type fractions to EMEP categories, to be used
    !// in dry deposition scheme.
    !// Used two times in drydeposition_oslo.f90 and one time in
    !// bcoc_oslo.f90.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: nLFmosaic, nLFin, LANDUSE_IDX
    real(r8), intent(in) :: lat
    real(r8), dimension(nLFin), intent(in) :: LFin
    real(r8), dimension(nLFmosaic), intent(out) :: LFmosaic
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'landfrac2mosaic'
    !// --------------------------------------------------------------------

    LFmosaic(:) = 0._r8

    !// Will set land categories, and then estimate ocean, since there
    !// may be slight inconsistencies between 1-PLAND and 1 - land
    !// categories. There should not be, but I have seen that
    !// 1 - land categories can be a tiny negative number.

    !// Vegetation type fractions (Fraction Land-use)
    if (LANDUSE_IDX .eq. 2) then
       !// ----------------------------------------------------------------
       !// ISLSCP2 MODIS land fraction and type data
       !// Note: Type 0 represents water and is not included in LS_FRAC.
       !//       It has been used to calculate land fraction PLAND.
       !//  0=Water Bodies                      1=Evergreen Needleleaf Forests
       !//  2=Evergreen Broadleaf Forests       3=Deciduous Needleleaf Forests
       !//  4=Deciduous Broadleaf Forests       5=Mixed Forests
       !//  6=Closed Shrublands                 7=Open Shrublands
       !//  8=Woody Savannas                    9=Savannas
       !// 10=Grasslands                       11=Permanent Wetlands
       !// 12=Croplands                        13=Urban and Built-Up
       !// 14=Cropland/Natural Vegetation Mosaic 15=Permanent Snow and Ice
       !// 16=Barren or Sparsely Vegetated       17=Unclassified
       !// ----------------------------------------------------------------
       !// mOSaic categories (reduced from Simpson et al., 2012)
       !//  1. Needleleaftree (temperated/boreal)
       !//  2. Deciduoustree (temperated/boral)
       !//  3. Needleleaftree (mediterranean)
       !//  4. Broadleaftree (mdeiterranean)
       !//  5. Crops <a,b,c>
       !//  6. Moorland (savanna++)
       !//  7. Grassland
       !//  8. Scrubs (med.)
       !//  9. Wetlands
       !//  10. Tundra
       !//  11. Desert
       !//  12. Water
       !//  13. Urban
       !//  14. Ice/Snow
       !// MODIS -> mOSaic
       !// 1: 1(90-45S, 45-90N),3, 0.5*5
       !// 2: 4, 0.5*5
       !// 3: 1 (45S-45N)
       !// 4: 2
       !// 5: 12,14
       !// 6: - [6,7,8,9]
       !// 7: 10
       !// 8: 6,7,8,9 [-]
       !// 9: 11
       !// 10: 16(90S-60S,60N-90N)
       !// 11: 16(60S-60N), 17
       !// 12: 1.d0 - PLAND(I,J)
       !// 13: 13
       !// 14: 15 (but will be treated seperately under snow/ice)
       !// ----------------------------------------------------------------
       LFmosaic(1) = LFin(3)+0.5_r8*LFin(5)
       LFmosaic(2) = LFin(4)+0.5_r8*LFin(5)
       LFmosaic(3) = 0._r8
       if (lat.lt.-45 .or. lat.gt.45) then
          LFmosaic(1) = LFmosaic(1)+LFin(1)
       else
          LFmosaic(3) = LFmosaic(3)+LFin(1)
       end if
       LFmosaic(4) = LFin(2)
       LFmosaic(5) = LFin(12) + LFin(14)
       !LFmosaic(6) = 0._r8
       LFmosaic(6) = sum(LFin(6:9))
       LFmosaic(7) = LFin(10)
       !LFmosaic(8) = sum(LFin(6:9))
       LFmosaic(8) = 0._r8
       LFmosaic(9) = LFin(11)

       if (lat.lt.-60._r8 .or. lat.gt.60._r8) then
          LFmosaic(10) = LFin(16)
          LFmosaic(11) = 0._r8
       else
          LFmosaic(10) = 0._r8
          LFmosaic(11) = LFin(16)
       end if
       LFmosaic(11) = LFmosaic(11) + LFin(17)
       LFmosaic(13) = LFin(13)
       LFmosaic(14) = LFin(15)
       !// Set Category "ocean" at the end
       LFmosaic(12) = max(0._r8, 1._r8 - (sum(LFmosaic(1:11)) + sum(LFmosaic(13:))))

    else if (LANDUSE_IDX .eq. 3) then

       !// --------------------------------------------------------------------
       !// mOSaic categories (reduced from Simpson et al., 2012)
       !//  1. Needleleaftree (temperated/boreal)
       !//  2. Deciduoustree (temperated/boral)
       !//  3. Needleleaftree (mediterranean)
       !//  4. Broadleaftree (mediterranean)
       !//  5. Crops <a,b,c>
       !//  6. Moorland (savanna++)
       !//  7. Grassland
       !//  8. Scrubs (med.)
       !//  9. Wetlands
       !//  10. Tundra
       !//  11. Desert
       !//  12. Water
       !//  13. Urban
       !//  14. Ice/Snow
       !// CLM categories
       !// CLM  CTM indices
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
       !// CLM -> mOSaic
       !// 1: 1 (90-45S, 45-90N) 2,3
       !// 2: 5,7,8
       !// 3: 1 (45S-45N)
       !// 4: 4,6
       !// 5: 15,16
       !// 6: - [14]
       !// 7: 13,14 [12,13]
       !// 8: 9,10
       !// 9: -
       !// 10: 11,12 [11]
       !// 11: 17
       !// 12: Ocean (1-PLAND)
       !// 13: -
       !// 14: -
       !// --------------------------------------------------------------------
       LFmosaic(1) = LFin(2)+ LFin(3)
       LFmosaic(2) = LFin(5) + sum(LFin(7:8))
       LFmosaic(3) = 0._r8
       if ((lat.lt.-45) .or. (lat.gt.45)) then
          LFmosaic(1) = LFmosaic(1)+LFin(1)
       else
          LFmosaic(3) = LFmosaic(3)+LFin(1)
       end if
       LFmosaic(4) = LFin(4) + LFin(6)
       LFmosaic(5) = sum(LFin(15:16))
       !LFmosaic(6) = 0._r8
       LFmosaic(6) = LFin(14)
       !LFmosaic(7) = LFin(13) + LFin(14)
       LFmosaic(7) = sum(LFin(12:13))
       LFmosaic(8) = LFin(9) + LFin(10)
       LFmosaic(9) = 0._r8
       !LFmosaic(10) = LFin(11) + LFin(12)
       LFmosaic(10) = LFin(11)
       if ((lat.gt.-60) .and. (lat.lt.60)) then
          LFmosaic(11) = LFin(17)
          LFmosaic(14) = 0._r8
       else
          LFmosaic(14) = LFin(17)
          LFmosaic(11) = 0._r8
       end if
       LFmosaic(13) = 0._r8
      
       !// Ocean may not be fully compatible with 1-PLAND:
       LFmosaic(12) = max(0._r8, 1._r8 - (sum(LFmosaic(1:11)) + sum(LFmosaic(13:))))

    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': unknown LANDUSE_IDX ',LANDUSE_IDX
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine landfrac2mosaic
  !// ----------------------------------------------------------------------
!// ----------------------------------------------------------------------
  subroutine set_vegetation_height(start_height,lat, out_height)
    !// --------------------------------------------------------------------
    !// MODIFY forest vegetation height dependent on latitude.
    !//
    !// Based on the modification polwards from 60deg North 
    !// in Simpson et al. 2012.
    !// Extended southward to the tropics and usage on southern hemesphere.
    !//
    !// Additional ESTIMATION: 
    !// Average Tropical forest is about 2-times higher (e.g. 40m) than 
    !// midlatitude forest. Decrease height linearly to 20m at 40deg.
    !// Constant height thereafter till 60deg. 
    !// Valide for both, northern and southern, hemisphere.
    !//
    !// Not sure if this is realistic, though.
    !//
    !// Used two times in drydeposition_oslo.f90 and one time in
    !// bcoc_oslo.f90.
    !//
    !// Stefanie Falk, March 2018
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: start_height, lat
    real(r8), intent(out) :: out_height
    !// Locals
    real(r8) :: rate
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_vegetation_height'
    !// --------------------------------------------------------------------
    rate = 0.05_r8

    if (abs(lat) .le. 10._r8) then
       out_height = start_height*2._r8
    else if (abs(lat) .le. 40._r8) then
       out_height = (1._r8-(abs(lat)-10._r8)*rate/3._r8)*start_height*2._r8
    else if (abs(lat) .le. 60._r8) then
       !// Forest vegetation height
       out_height = start_height
    else if (abs(lat) .lt. 74._r8) then
       !// From 60-74 latitude, linearly decreasing with rate
       out_height = (1._r8-(abs(lat)-60._r8)*rate)*start_height
    else
       !// High-latitude forest
       out_height = start_height*3._r8/10._r8
    end if

    if (out_height .lt. 0._r8) then
       write(6,'(a,i5)') f90file//':'//subr// &
            ': There is something wrong with your vegetation height!'
       write(6,'(a)') 'Start (m): '
       write(6,*) start_height
       write(6,'(a)') 'Stop (m): '
       write(6,*), out_height
       stop
    end if
    !// --------------------------------------------------------------------
  end subroutine set_vegetation_height
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  subroutine MAPPED_GROWSEASON ( JDAY, I, J, GDAY, GLEN )
    !// --------------------------------------------------------------------
    !  DESCRIPTION
    !     This function uses a preprocessed map of grow season GDAY and GLEN.
    !     The preprocessed data was created using the 5 days/5 degC
    !     criterum to estimate the begin and end of growing season
    !     from the OpenIFS meteorological data of the same year
    !     for latitudes between 45 and 85deg north.
    !     All others are basically the same as computed with GROWSEASON.
    !
    !// Stefanie Falk, August 2018
    !//   Based on GROWSEASON (below).
    !//   Enable usage of preprocessed GDAY/GLEN
    !// --------------------------------------------------------------------
    use cmn_sfc, only: GDAY_MAP, GLEN_MAP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)  :: JDAY
    integer, intent(in)  :: I, J
    integer, intent(out) :: GDAY, GLEN
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'MAPPED_GROWSEASON'
    !// --------------------------------------------------------------------
    GLEN = GLEN_MAP(I,J)
    GDAY = GDAY_MAP(I,J,JDAY)
    
  end subroutine MAPPED_GROWSEASON
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  subroutine GROWSEASON ( JDAY, LAT, GDAY, GLEN )
    !// --------------------------------------------------------------------
    !  DESCRIPTION
    !    This function computes the day of the growing season
    !    corresponding to the given day of year.
    !
    !  HISTORY:
    !    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
    !               Variation of growing season depends on latitude
    !               (Guenther)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    use cmn_ctm, only: LYEAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)  :: JDAY
    real(r8), intent(in) :: LAT
    integer, intent(out) :: GDAY, GLEN
    !// Locals
    integer  :: &
         GSEASON_START, GSEASON_END, &
         MAY31,NOV1,DEC31
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'GROWSEASON'
    !// --------------------------------------------------------------------
           
    !// Days of year for important dates, taking leap year into account
    if (LYEAR) then
       MAY31 = 152
       NOV1  = 306
       DEC31 = 366
    else
       MAY31 = 151
       NOV1  = 305
       DEC31 = 365
    end if
    
    if ( LAT .le. 23.0_r8 .and. LAT .ge. -23.0_r8 ) then
       ! tropical regions, year round
       GSEASON_START = 1
       GSEASON_END = DEC31
          
       GDAY = JDAY - GSEASON_START + 1
       GLEN = GSEASON_END - GSEASON_START + 1
          
    else if ( LAT .lt. -23.0_r8 ) then
       ! southern hemisphere
       if ( LAT .lt. -60.0_r8 ) then
          ! antarctic start = 0 end = 0, no growing
          GDAY = 0
          GLEN = 0
       else
          ! southern hemisphere temperate, NOV, DEC, JAN-MAY
          if (JDAY .ge. NOV1 .AND. JDAY .le. DEC31 ) THEN
             GSEASON_START = NOV1
             GSEASON_END   = DEC31
             
             GDAY = JDAY - GSEASON_START + 1
             !// Total length of growing season (Nov-May)
             GLEN = DEC31-NOV1 + 1 + MAY31
          else if (JDAY .ge. 1 .and. JDAY .le. MAY31) then
             GSEASON_START = 1
             GSEASON_END   = MAY31
             !// Add Nov+Dec
             GDAY = JDAY - GSEASON_START + 1 + 61
             GLEN = 61 + GSEASON_END - GSEASON_START + 1
             !// Total length of growing season (Nov-May)
             GLEN = DEC31-NOV1 + 1 + MAY31
          else
             GDAY = 0
             GLEN = 0
          end if
          
       end if

    else if ( LAT .gt. 23._r8 ) then

       ! northern hemisphere
       if ( LAT .gt. 65.0_r8 ) then
          ! arctic start = 0 end = 0, no growing season
          GDAY = 0
          GLEN = 0
       else
          ! northern hemisphere temperate
          ! start= (lat-23)*4.5            189
          ! end = 365 -((lat-23)*3.3)      226
          
          GSEASON_START = INT( (LAT - 23.0_r8) * 4.5_r8 )
          GSEASON_END   = DEC31 - INT( (LAT - 23.0_r8) * 3.3_r8 )
          
          if (JDAY .ge. GSEASON_START .and. JDAY .le. GSEASON_END) then
             GDAY = JDAY - GSEASON_START + 1
          else
             GDAY = 0
          end if
          GLEN = GSEASON_END - GSEASON_START + 1
       end if
    else
       write(6,'(a,es16.6)') f90file//':'//subr//': Invalid LAT ',LAT
       stop
    end if
    
    !// --------------------------------------------------------------------
  end subroutine GROWSEASON
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
end module utilities_oslo
!//=========================================================================

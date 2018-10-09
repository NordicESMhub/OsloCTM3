!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, November 2017
!//=========================================================================
!// PSC microphysics
!//=========================================================================
module psc_microphysics
  !// ----------------------------------------------------------------------
  !// MODULE: psc_microphysics
  !// DESCRIPTION: Routines for PSC chemistry (stratosphere).
  !//
  !// Parameters and distributions are so far constant, since the mean radius
  !// is constant. If the MEAN RADIUS is to be CHANGED, this will have to be
  !// revised!
  !//
  !// Contains:
  !//   subroutine oslochem_psc
  !//   subroutine get_surface_area
  !//   subroutine get_psc12_sad
  !//   subroutine psc_box
  !//   subroutine get_mw:
  !//   subroutine CARS
  !//   function H2O_SAT
  !//   function WSED
  !//   SUBROUTINE SED
  !//   FUNCTION AM_BIN
  !//   FUNCTION X_CARS
  !//   FUNCTION HENRIC
  !//   FUNCTION RHO_S_CAR
  !//   FUNCTION RHO_N_CAR
  !//   SUBROUTINE LOGN
  !//   SUBROUTINE MOM
  !//   SUBROUTINE DIS_LN
  !//   SUBROUTINE sps_SURF
  !//   SUBROUTINE sps_WSASAS
  !//   FUNCTION sps_ROSAS
  !//   SUBROUTINE PSC_diagnose
  !//   subroutine set_psc_constants
  !//   function getTnat
  !//
  !// Amund Sovde Haslerud, November 2017
  !//   Cleaned up and made more efficient. Included max height ZHMAX for
  !//   calculations.
  !// Ole Amund Sovde, October 2008
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR
  use cmn_parameters, only: KBOLTZ, AVOGNR, R_ATM, LDEBUG, CPI
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Number of bins in log-normal distribution
  !// Originally, there were 40 bins, but that made the code rather slow.
  !// To reduce computing time for microphysics, the number of bins
  !// are reduced to 15.
  integer, parameter  :: N_B = 15   ! Number of radii

  !// Flags for turning on/off SAD/PSC1/PSC2
  logical, parameter :: &
       LAEROSOL = .true., &  ! Calculate SAD for backround aerosols
       LPSC = .true., &      ! Calculate PSC SAD
       LPSC2 = .true.        ! Specifically also for PSC2

  real(r8), parameter  :: &
       R0_SA  = 0.07e-4_r8, &       !// mean radii SAD [cm]
       R0_NAT = 0.4e-4_r8, &        !// :          NAT [cm]
       R0_ICE = 5.0e-4_r8, &        !// :          ICE [cm]
       S_SA   = 1.9_r8, &           !// standard deviations
       S_NAT  = 1.6_r8, &           !// :
       S_ICE  = 1.8_r8, &           !// :
       TICE_MIN = 170._r8, &     !// Minimum value of TICE (see TICE for reference)
       ZHMAX = 30000._r8         !// Max allowed height for PSC existence [m]

  real(r8) :: &          !// To be calculated once in set_psc_constants()
       F_SAD(N_B), &     !// Mass distribution SAD
       F_NAT(N_B), &     !// Mass distribution NAT
       F_ICE(N_B), &     !// Mass distribution ICE
       R_A(N_B), &       !// Radii of size distribution [cm]
       DR(N_B), &        !// Difference between the radii [cm]
       VOL_PAR(N_B),&    !// Volume of particles in bin (4/3*pi*r^3dr)
       RA, &             !// 1. momentum of log-normal distribution, SAD
       RNAT, &           !// 1. momentum of log-normal distribution, NAT
       RICE, &           !// 1. momentum of log-normal distribution, ICE
       S2_SA, &          !// 2. momentum of log-normal distribution, SAD
       S2_NAT, &         !// 2. momentum of log-normal distribution, NAT
       S2_ICE, &         !// 2. momentum of log-normal distribution, ICE
       S3_SA, &          !// 3. momentum of log-normal distribution, SAD
       S3_NAT, &         !// 3. momentum of log-normal distribution, NAT
       S3_ICE            !// 3. momentum of log-normal distribution, ICE


  !// Coefficients for the calculation of the concentration of a
  !// binary solution of HNO3 and H2SO4 (Carslaw et al.,GRL,22,1877-1880,1995)
  real(r8)  :: AK_CARS(7,2)

  !// Coefficients for the calculation of Henry's law constant
  !// of HNO3 and H2SO4 (Carslaw et al.,GRL,22,1877-1880,1995)
  real(r8)  :: Q_CARS(10,2)

  !// PSC 3D variables
  real(r8), dimension(LPAR,IPAR,JPAR) :: &
       PSC1, PSC2, & !// PSC1 and PSC2 surface area densities
       VAER, &       !// Volume of frozen aerosols
       SPS_PARTAREA  !// PARTAREA after PSC modification
 
  !// Logical to tell if SAT or not
  logical, dimension(LPAR,IPAR,JPAR) :: SAT


  !// Coefficients for calculating sedimentation velocities in function
  !// wsed (Kasten, 1968: Falling speed of aerosol particles, J. Appl.
  !// Meteorol. 7, pp. 944-947,
  !// doi:10.1175/1520-0450(1968)007<0944:FSOAP>2.0.CO;2 )
  integer, parameter  :: wsed_LL = 31   !// Values for 31 levels
  integer, parameter  :: wsed_NN = 8    !// Size bins
  real(r8) :: wsed_AMN(wsed_NN)         !// Exponents for wsed_WF
  real(r8) :: wsed_RF(wsed_NN)          !// Particle radius [microns]
  real(r8) :: wsed_WF(wsed_LL,wsed_NN)  !// Sedimentation velocities 
  real(r8) :: wsed_ZF(wsed_LL)          !// Height of levels
  !// Lower limit for treating VAER
  real(r8), parameter :: vaer_lim = 1.e-15_r8
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'psc_microphysics.f90'
  !// ----------------------------------------------------------------------
  save !// All variables are to be saved.
  private
  public oslochem_psc, get_psc12_sad, set_psc_constants, LPSC, LAEROSOL, &
       psc_diagnose
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_psc(BTT,BTEM,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Calculates formation and evolution of PSCs/STS.
    !//
    !// Will calculate PSC1 and PSC2 areas.
    !//
    !// Amund Sovde Haslerud, 2017
    !//   Optimised and cleaned up.
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LOSLOCSTRAT
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, ETAA, ETAB
    use cmn_chem, only: TMASS
    use cmn_met, only: ZOFLE, P
    use cmn_parameters, only: M_AIR
    use cmn_oslo, only: LMTROP, PARTAREA, trsp_idx, Xtrsp_idx, DV_IJ, &
         AIRMOLEC_IJ, XSTT
    use strat_h2o, only: set_d_h2o, LOLD_H2OTREATMENT, STR_H2O
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In
    integer, intent(in) :: MP
    real(r8), intent(in)  :: DTCHM
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM
    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout):: BTT

    !// Locals

    !// Column arrays ------------------------------------------------------
    !// box volume [m], temperature [K]
    real(r8), dimension(LPAR) :: DV, TEMP, ZDV

    !// Air density [molecules / cm3]
    real(r8) :: AIR_MOLEC

    !// Tracer concentrations [molecules / cm3]
    real(r8) :: C_H2SO4_gas !// Column not needed
    real(r8), dimension(LPAR) :: &
         C_HNO3_gas, C_H2O_gas, C_NOy, C_HNO3_iceliq, C_H2O_ice, C_NOx, &
         C_HNO3_ice, C_HNO3_liq

    !// Old tracer concentrations for diagnostics [molecules / cm3]
    real(r8),dimension(LPAR) :: &
         old_HNO3_gas, old_H2O_gas, &
         old_HNO3_ice, old_HNO3_iceliq

    !// SAD, PSC1, PSC2 columns
    real(r8), dimension(LPAR) :: sSAD, sPSC1, sPSC2, sVAER, &
         sVICE, sVAER_OLD, sRHO, & ! STS density
         n_particles, &            ! number of particles
         sTNAT, &                  ! TNAT
         sTICE                     ! TICE
    logical, dimension(LPAR) :: sSAT

    !// Altitude center [km] and upper edge of box[m]
    real(r8), dimension(LPAR)   :: ALT, ZTOP
    !// Altitude lower edge of box [cm]
    real(r8), dimension(LPAR+1) :: ALTST

    !// --------------------------------------------------------------------

    !// Inverse molecular masses multiplied with Avogadros number
    real(r8) :: Z_HNO3, Z_NOy, Z_H2O, IZ_HNO3, IZ_NOy, IZ_H2O

    !// For looping
    integer :: I,J,L,II,JJ

    !// End layer of tropospherer and start layer of stratosphere
    integer :: LMT, LOWERL

    !// Highest level below ZHMAX
    integer :: LMAX

    !// Transport indices
    integer :: N_HNO3, N_NOy, N_HNO3s, N_H2O, N_H2Os

    !// For H2SO4 calculation
    real(r8) :: PRES, H2SO4_MASMR, SAREA

    !// Average NAT size
    real(r8) :: accNATsize, countNAT

    real(r8), dimension(LPAR)   :: dH2Osed
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'oslochem_psc'
    !//---------------------------------------------------------------------

    !// If no stratospheric chemistry, skip PSC routine
    if (.not. LOSLOCSTRAT) return

    !// Skip if .not.LPSC
    if (.not.LPSC) return

    !// Get molecular masses for transforming to concentration
    !// Can be fetched from get_mw, hard code for now;
    N_HNO3  = trsp_idx(  4)
    N_H2O   = trsp_idx(114)
    N_HNO3s = trsp_idx(124)
    N_H2Os  = trsp_idx(125) ! can be non-transported
    N_NOy   = trsp_idx(147)

    Z_HNO3  = AVOGNR * 1.e-3_r8 / TMASS(N_HNO3)
    Z_NOy   = AVOGNR * 1.e-3_r8 / TMASS(N_NOy)
    Z_H2O   = AVOGNR * 1.e-3_r8 / 18._r8
    !// Inverse coefficients
    IZ_HNO3 = 1._r8 / Z_HNO3
    IZ_NOy  = 1._r8 / Z_NOy
    IZ_H2O  = 1._r8 / Z_H2O


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1


        !// Height of grid box tops and the uppermost box below ZHMAX
        LMAX = 1
        do L = 1, LPAR
           ZTOP(L) = ZOFLE(L+1,I,J) - ZOFLE(1,I,J)
           if (ZTOP(L) .le. ZHMAX) LMAX = L
        end do

        !// Startpoint of PSC calculations (lowermost stratosphere)
        LMT = LMTROP(I,J)
        LOWERL = LMTROP(I,J) + 1

        if (LMAX .le. LOWERL) then
           write(6,'(a,2i5)') f90file//':'//subr// &
                ': LMAX < LMTROP+1',LMAX,LOWERL
           stop
        end if


        !// Set variables in L=1,LPAR (LMAX)
        !// Formation is calculated in stratosphere only,
        !// but sedimentation needs values below also.
        !//-----------------------------------------------------------------
        !// Temperature
        TEMP(1:LMAX) = BTEM(1:LMAX,II,JJ)
        !// Volume [m3]
        DV(1:LMAX) = DV_IJ(1:LMAX,II,JJ,MP)
        !// Specific volume [m3]
        ZDV(1:LMAX) = 1._r8 / DV(1:LMAX)

        !// Tracer concentrations [molec/cm^3]
        C_HNO3_gas(1:LMAX)    = BTT(1:LMAX,N_HNO3, II,JJ) * ZDV(1:LMAX) * Z_HNO3
        !// iceliq is combined tracer, and works as input to psc_box.
        !// There it may get altered slightly, but the main thing psc_box
        !// does is to calculate HNO3_liq and HNO3_ice, which are
        !// treated by sedimentation before HNO3_iceliq is set to their sum.
        C_HNO3_iceliq(1:LMAX) = BTT(1:LMAX,N_HNO3s,II,JJ) * ZDV(1:LMAX) * Z_HNO3
        C_NOy(1:LMAX)         = BTT(1:LMAX,N_NOy,  II,JJ) * ZDV(1:LMAX) * Z_NOy
        if (N_H2Os .gt. 0) then
           C_H2O_ice(1:LMAX)  = BTT(1:LMAX,N_H2Os, II,JJ) * ZDV(1:LMAX) * Z_H2O
        else
           C_H2O_ice(1:LMAX)  = 0._r8 !// Old treatment needs zero
        end if
        if (N_H2O .gt. 0) then
           C_H2O_gas(1:LMAX)  = BTT(1:LMAX,N_H2O, II,JJ) * ZDV(1:LMAX) * Z_H2O
        else
           C_H2O_gas(1:LMAX)  = STR_H2O(1:LMAX,II,JJ,MP) * ZDV(1:LMAX) * Z_H2O
        end if
        !// Keep track of NOx
        !// Solid HNO3 is not included in NOy, so NOx is just
        !// NOx = NOy - HNO3
        C_NOx(1:LMAX) = max(0._r8, C_NOy(1:LMAX) - C_HNO3_gas(1:LMAX))

        !// Save previous values for diagnostics
        old_HNO3_gas(1:LMAX)    = C_HNO3_gas(1:LMAX)
        old_HNO3_iceliq(1:LMAX) = C_HNO3_iceliq(1:LMAX)
        old_H2O_gas(1:LMAX)     = C_H2O_gas(1:LMAX)

        !// Initialise HNO3_ice and HNO3_liq (important for sedimentation)
        C_HNO3_ice(1:LMAX) = 0._r8
        C_HNO3_liq(1:LMAX) = 0._r8


        !// Set whole column for some variables
        sSAD(:) = PARTAREA(:,I,J) !// Need to set whole column
        n_particles(:) = 0._r8    !// Maybe not necessary, but keep it
        !// It could be enough to set these in 1:LMT only
        sRHO(:) = 1.62_r8         !// Need troposphere for sedimentation
        sVAER_OLD(:) = 0._r8      !// Need troposhere zero for sedimentation
        dH2Osed(:) = 0._r8      !// Sedimented H2O

        !// Initialise LOWERL:LMAX
        !//-----------------------------------------------------------------
        !// Volume of aerosols
        sVAER(LOWERL:LMAX) = VAER(LOWERL:LMAX,I,J)


        !// Initialise troposphere
        !//-----------------------------------------------------------------
        !// No PSCs (or particle volumes) in troposphere.
        sPSC1(1:LMT) = 0._r8
        sPSC2(1:LMT) = 0._r8
        sSAT(1:LMT) = .false.
        sVAER(1:LMT) = 0._r8
        sVICE(1:LMT) = 0._r8
        !// TNAT and TICE should be initialised. SHould not necessary,
        !// they are only used in LOWERL:LMAX.
        sTNAT(1:LMT) = 0._r8
        sTICE(1:LMT) = 0._r8


        !// Do evaporation in troposphere (probably better? to do after sed)
        !//-----------------------------------------------------------------
        !// Gaseous species; evaporate solid species
        C_HNO3_gas(1:LMT) = C_HNO3_gas(1:LMT) + C_HNO3_iceliq(1:LMT)
        C_HNO3_iceliq(1:LMT) = 0._r8

        !// No need to change NOy, it is set afterwards from NOx+new_HNO3_gas

        !// No solid H2O in troposphere
        C_H2O_gas(1:LMT) = C_H2O_gas(1:LMT) + C_H2O_ice(1:LMT)
        C_H2O_ice(1:LMT) = 0._r8


        !// Need to calculate average NAT size for sedimentation
        accNATsize = 0._r8
        countNAT = 0._r8


        do L = LOWERL, LMAX

           !// Cannot check on T > 240K here:
           !// There may have been particles before which needs to evaporate

           !// Air molecules
           AIR_MOLEC = AIRMOLEC_IJ(L,II,JJ,MP)

           !// H2SO4 concentration
           !// H2SO4 mass mixing ratio is calculated from aerosol surface
           !// density area, assuming the aerosols to be described by a
           !// log-normal distribution, and the particles in equilibrium
           !// with water. Alternatively, it can be set maually.

           !// Calculated from aerosol surface density area:
           !// Box center pressure [hPa]
           PRES = 0.5_r8 * ( (ETAA(L) + ETAA(L+1)) &
                          + (ETAB(L) + ETAB(L+1)) * P(I,J) )
           !// Calculate H2SO4 mass mixing ratio
           CALL sps_SURF(TEMP(L),PRES,C_H2O_GAS(L) / AIR_MOLEC,sSAD(L),H2SO4_MASMR)
           !// Convert to concentration.
           !// The CARS-routine ignores negative values of H2SO4, but we set
           !// negative values to zero just in case.
           C_H2SO4_gas = max(H2SO4_MASMR * M_AIR / 98._r8 * AIR_MOLEC, 0._r8)

           !// Alternatively H2SO4(L) can be set manually, e.g. with a
           !// typical vol. mix. rat. of 0.5ppbv:
           !//   IF (L.GE.LOWERL) THEN
           !//     C_H2SO4_gas = 0.5d-9 * AIR_MOLEC
           !//   ELSE
           !//     C_H2SO4_gas = 0._r8
           !//   END IF

           !// Note that calculations are only done for LOWERL:LMAX
           !// HNO3_iceliq is not used in sedimentation, but after
           !// HNO3_ice and HNO3_liq have been sedimented, HNO3_iceliq
           !// is set to their sum before going back to tracer 124.
           call psc_box(TEMP(L), sSAD(L), &
                C_HNO3_gas(L), C_H2O_gas(L), C_H2SO4_gas, &        ! Gases
                C_HNO3_iceliq(L), C_HNO3_ice(L), C_HNO3_liq(L), &  ! Solid
                C_H2O_ice(L), &
                sVAER(L), sSAT(L), accNATsize, countNAT, &
                sVICE(L), sVAER_OLD(L), sRHO(L), &
                n_particles(L), sTNAT(L), sTICE(L), &
                I,J,L)

        end do !// do L = LOWERL, LMAX


        !// Sedimentation for column
        !// ----------------------------------------------------------------
        !// If no particles, go to next grid box
        if (sum(sVAER(1:LMAX)) .eq. 0._r8 .and. &
            sum(sVICE(1:LMAX)) .eq. 0._r8) then
           if (sum(sVAER_OLD(1:LMAX)) .eq. 0._r8) then
              !// There were no PSCs before, no need to update tracers
              write(6,'(a)') f90file//':'//subr// &
                   ': None before, none now'
              PSC1(:,I,J) = 0._r8
              PSC2(:,I,J) = 0._r8
              cycle
           else
              !// There has been evaporation of all particles.
              !// set gasesous species and skip sedimentation.
              BTT(1:LMAX,N_HNO3s,II,JJ) = 1.e-20_r8
              BTT(1:LMAX,N_HNO3, II,JJ) = C_HNO3_gas(1:LMAX) * DV(1:LMAX) * IZ_HNO3
              BTT(1:LMAX,N_NOy,  II,JJ) = &
                      (C_NOx(1:LMAX) + C_HNO3_gas(1:LMAX)) * DV(1:LMAX) * IZ_NOy
              if (N_H2O .gt. 0) then
                 BTT(1:LMAX,N_H2O,II,JJ) = &
                         max(C_H2O_gas(1:LMAX) * DV(1:LMAX) * IZ_H2O, 1.e-20_r8)
              else
                 STR_H2O(1:LMAX,II,JJ,MP) = &
                         max(C_H2O_gas(1:LMAX) * DV(1:LMAX) * IZ_H2O, 1.e-20_r8)
              end if
              if (N_H2Os .gt. 0) then
                 BTT(1:LMAX,N_H2Os,II,JJ) = 1.e-20_r8
              else
                 XSTT(1:LMAX,Xtrsp_idx(125),I,J) = 1.e-20_r8
              end if

              PSC1(:,I,J) = 0._r8
              PSC2(:,I,J) = 0._r8

              write(6,'(a)') f90file//':'//subr// &
                   ': All particles were evaporated'
              cycle
           end if
        end if

        !// Set altitudes
        ALT(1)   = ZTOP(1) * 0.5e-3_r8  !// Center altitude [km]
        ALTST(1) = 0._r8                !// Edge altitude [cm]
        do L = 2, LPAR
          ALT(L)   = (ZTOP(L) + ZTOP(L-1)) * 0.5e-3_r8
          ALTST(L) = ZTOP(L-1) * 100._r8
        end do
        ALTST(LPAR+1) = ZTOP(LPAR) * 100._r8


        call sedimentation_driver(ALT, ALTST, accNATsize, countNAT, DTCHM, &
             sVAER, sVICE, sVAER_OLD, sRHO, C_H2O_ice, &
             C_HNO3_ice, C_HNO3_liq, &
             dH2Osed, LMAX)

        !// Total not gaseous HNO3 after sedimentation
        !// ----------------------------------------------------------------
        C_HNO3_iceliq(1:LMAX) = C_HNO3_liq(1:LMAX) + C_HNO3_ice(1:LMAX)

        !// In case of negative / zero values
        do L = 1, LMAX
           if (C_HNO3_iceliq(L) .lt. 0._r8) then
              write(6,'(a,4i5,2es12.3)') f90file//':'//subr// &
                   ': NEG C_HNO3_iceliq (MP,I,J,L)',MP,I,J,L, &
                   C_HNO3_iceliq(L), C_HNO3_ice(L)
              stop
              C_HNO3_iceliq(L) = max(C_HNO3_iceliq(L), 0._r8)
           end if

           if (C_H2O_ice(L) .lt. 0._r8) then
              write(6,'(a,4i5,2es12.3)') f90file//':'//subr// &
                   ': NEG C_H2O_ice (MP,I,J,L)',MP,I,J,L, &
                      C_H2O_ice(L)
              stop
              C_H2O_ice(L)  = max(C_H2O_ice(L), 0._r8)
           end if
        end do

        if (LDEBUG) then
           do L = 1, LMAX
              if (C_HNO3_gas(L) .ne. C_HNO3_gas(L)) then
                 write(6,'(a,4i5,6es12.3)') f90file//':'//subr// &
                      ': PSC problem (MP,I,J,L)',MP,I,J,L, &
                      old_HNO3_gas(L),C_HNO3_gas(L),DV(L),TEMP(L), &
                      old_H2O_gas(L),C_H2O_gas(L)
                 stop
              end if
           end do !// do L = 1, LMAX
        end if !// if (LDEBUG) then







        do L = LOWERL, LMAX
           !// Calculate surface area
           !// Only in stratosphere
           call get_surface_area(TEMP(L), sTNAT(L), sTICE(L), &
                sSAT(L), sVAER(L), sVICE(L), &
                sPSC1(L), sPSC2(L), sSAD(L), &
                n_particles(L), I,J,L)

           !// "EQULIBRIUM APPROACH" TO VOLUME
           !// -------------------------------------------------------------
           !// As in the latest SCTM1 - but what is really happening here?
           !// 17/9-2004: The answer from Sergei:
           !// This an idea to combine non-equilibrium and equilibrium
           !// approach. VOL_OLD is frozen aerosol (SAT) and its amount is
           !// permanent until it reaches melting temperature. The rest HNO3
           !// is a subject for equilibrium treatment like in original
           !// Carslaw code. Their sum is combined to evaluate surface area
           !// but then they devided for futher treatment.
           !// However, try to comment this line and let me know the result.
           !// 
           !// If STS is frozen, it is assumed composition of NAT, and the
           !// VOLA is then amount of frozen aerosol, which is assumed
           !// constant. 
           !// 2011: Threw this back in, since PSC area density grows too fast
           !// and creates a lot more O3 loss than observed. Ideally, it should
           !// not be used since particles should grow.
           if (sSAT(L)) sVAER(L) = sVAER_OLD(L)
           !// -------------------------------------------------------------

        end do



        !// UPDATE AEROSOL AREA - whole column
        do L = 1, LPAR
           !// AR is read from satellite data, and is used for H2SO4
           !// calculations. A new SAD is calculated from the H2SO4, so AR
           !// should not be overwritten and then used the next time step.
           !// To solve this, sps_PARTAREA is updated from satellite data
           !// each time step before entering this subroutine. We can then
           !// update AR to use in the heterogeneous chemistry.
           if (sSAD(L) .ne. sSAD(L)) then
              write(6,'(a,4i5)') f90file//':'//subr// &
                   ': sSAD(L) NAN (MP,I,J,L)',MP,I,J,L
              stop
           end if

           SPS_PARTAREA(L,I,J) = sSAD(L)
        end do

        !// UPDATE NEW VALUES and calculate parameters / psc rates
        !// Troposphere: solid HNO3/H2O should evaporate
        C_HNO3_gas(1:LMT)    = C_HNO3_gas(1:LMT) + C_HNO3_iceliq(1:LMT)
        C_HNO3_iceliq(1:LMT) = 0._r8
        C_H2O_gas(1:LMT)     = C_H2O_gas(1:LMT) + C_H2O_ice(1:LMT)
        C_H2O_ice(1:LMT)     = 0._r8


        !// Put back masses
        BTT(1:LMAX,N_HNO3s,II,JJ) = &
             max(C_HNO3_iceliq(1:LMAX) * DV(1:LMAX) * IZ_HNO3, 1.e-20_r8)

        BTT(1:LMAX,N_HNO3,II,JJ) = &
             max(C_HNO3_gas(1:LMAX) * DV(1:LMAX) * IZ_HNO3,1.e-20_r8)

        BTT(1:LMAX,N_NOy,II,JJ) = &
             max((C_NOx(1:LMAX) + C_HNO3_gas(1:LMAX)) &
                  * DV(1:LMAX) * IZ_NOy, 1.e-20_r8)

        if (N_H2O .gt. 0) then
           BTT(1:LMAX,N_H2O,II,JJ) = &
                max(C_H2O_gas(1:LMAX) * DV(1:LMAX) * IZ_H2O,1.e-20_r8)
        else
           if (minval(C_H2O_gas(1:LMAX)) .le. 0._r8) then
              do L = 1, LMAX
                 write(6,'(a,4i5,5es12.3)') f90file//':'//subr// &
                      ': C_H2O_gas <= 0 (MP,I,J,L)',MP,I,J,L, &
                      old_H2O_gas(L),C_H2O_gas(L), &
                      C_H2O_ice(L),DV(L),TEMP(L)
              end do
              stop
           end if
           STR_H2O(1:LMAX,II,JJ,MP) = &
                max(C_H2O_gas(1:LMAX) * DV(1:LMAX) * IZ_H2O, 1.e-20_r8)
        end if

        if (N_H2Os .gt. 0) then
           BTT(1:LMAX,N_H2Os,II,JJ) = &
                max(C_H2O_ice(1:LMAX) * DV(1:LMAX) * IZ_H2O,1.e-20_r8)
        else
           XSTT(1:LMAX,Xtrsp_idx(125),I,J) = &
                max(C_H2O_ice(1:LMAX) * DV(1:LMAX) * IZ_H2O,1.e-20_r8)
        end if

        !// PSC surface area density (not only LOWERL:LMAX, to get zero in
        !// troposphere).
        PSC1(1:LMAX,I,J) = sPSC1(1:LMAX)
        PSC2(1:LMAX,I,J) = sPSC2(1:LMAX)

        !// Set zero PSCs above LMAX
        PSC1(LMAX+1:LPAR,I,J) = 0._r8
        PSC2(LMAX+1:LPAR,I,J) = 0._r8

        !// Send sedimented diagnose to strath2o.f90
        call set_d_h2o(II,JJ,MP,LMAX,dH2Osed)

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine oslochem_psc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_surface_area(TEMP, TNAT, TICE, &
       SAT, VAER, VICE, PSC1_L, PSC2_L, SAD_L, &
       n_particles, I,J,L)
    !// --------------------------------------------------------------------
    !// Calculates surface area density of PSC1/PSC2/SAD based on
    !// the number of particles calculated and their volumes.
    !//
    !// Taken from the original PSC_1d code, but calculates only for a
    !// given level.
    !//
    !// Amund Sovde Haslerud, November 2017
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: TEMP, TNAT, TICE, VAER, VICE
    real(r8), intent(inout) :: n_particles
    logical, intent(in) :: SAT
    integer, intent(in) :: I,J,L
    real(r8), intent(inout) :: SAD_L
    real(r8), intent(out) :: PSC1_L, PSC2_L

    !// Locals
    real(r8) :: SURFACE, S2_NAT_loc, S3_NAT_loc, r_particle, &
         RNAT_loc, F_NAT_loc(N_B)
    integer :: METHOD
    real(r8), parameter :: PI43 = CPI * 4._r8/3._r8  ! Pi*4/3
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_surface_area'
    !//---------------------------------------------------------------------


    !// Calculate surface area of PSC1, PSC2 and SAT
    !// --------------------------------------------------------------------
    !// PSCI
    !// --------------------------------------------------------------------
    if (TEMP .gt. 220._r8) then
       PSC1_L = 0._r8
       PSC2_L = 0._r8
       !// Do not touch SAD
       return
    else

       S2_NAT_loc = 0._r8
       S3_NAT_loc = 0._r8

       !// Use TMELT to distinguish between NAT and SA.
       if (SAT .and. VAER .gt. vaer_lim) then
          !// Treat particles as NAT
          !// n_particles has already been set, and is not sedimented.
          !// (Maybe it should.)
          !// We need to find S2_NAT_loc and S3_NAT_loc from
          !// each of the model levels.
          !//
          !// If n_particles would be sedimented, the minimum value must be
          !// set once again.
          !//
          !// Seems I did not get this right:
          !if (n_particles .le. 0._r8) then
          !   write(6,'(a,3i5,6es12.3)') f90file//':'//subr// &
          !       ': n_particles == 0: ', I, J, L, n_particles, VAER,TEMP, SAD_L
          !   stop
          !end if
          !// n_particles may somehow be
          !// zero while VAER is close to vaer_lim, and SAT true while
          !// SAD_L is zero. In other words, tiny NAT with no sulphate.
          !// Could either skip and set SURFACE=0 or set n_particles=2,
          !// Try the latter:
          if (n_particles .le. 0._r8) n_particles = 2._r8
          
          r_particle = (VAER / (PI43 * n_particles))**(1._r8 / 3._r8) * 0.595238_r8
          call DIS_LN(N_B,R_A,F_NAT_loc,r_particle,S_NAT, &
                  RNAT_loc,S2_NAT_loc,S3_NAT_loc)

          if (S3_NAT_loc .gt. 0._r8) then
             SURFACE = 3._r8 * VAER * S2_NAT_loc / S3_NAT_loc ! cm(2)/cm(3)
          else
             SURFACE = 0._r8
          end if

          METHOD = 1
       else if (.not.SAT .and. VAER .gt. vaer_lim) then

          !// Only STS or SA
          SURFACE = 3._r8 * VAER * S2_SA / S3_SA ! cm(2)/cm(3)
          METHOD = 2
       else

          !// No NAT, no STS, no SA
          !// May change the background aerosol area. Should not have
          !// a large impact.
          SURFACE = 0._r8
          METHOD = 0
       end if


       if (SURFACE.ne.SURFACE) then
          write(6,'(a,3i5,6es12.3)') f90file//':'//subr// &
               ': SURFACE NAN ', I, J, L, SURFACE,VAER,S2_NAT_loc, &
               S3_NAT_loc,TEMP, SAD_L
          stop
       end if
          
       if (Temp .le. TNAT) then
          !// PSC (can be NAT or STS), no SAD
          PSC1_L = SURFACE
          SAD_L  = 0._r8
       else
          !// particle is ternary or binary solution, but not supercooled
          SAD_L  = SURFACE
          !// i.e. no PSCs
          PSC1_L = 0._r8
       end if
    end if

    !// PSCII
    !// --------------------------------------------------------------------
    if (Temp .le. TICE .and. Temp .lt. 210._r8) then
       !// Used to be 190K instead of TICE
       !// But we limit to 210 in case of overshooting when H2O is
       !// transported, which causes problems in stratchem for low
       !// tropopauses.
       PSC2_L = 3._r8 * VICE * S2_ICE / S3_ICE ! cm(2)/cm(3)
    else
       PSC2_L = 0._r8
    end if

    !// --------------------------------------------------------------------
  end subroutine get_surface_area
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine sedimentation_driver(ALT, ALTST, accNATsize, countNAT, TSTEP, &
       VOLA, VICE, VOL_OLD, RHO, H2Os, HNO3I, HNO3L, &
       dH2Osed, LMAX)
    !// --------------------------------------------------------------------
    !// SEDIMENTATION
    !// HNO3 in STS/NAT and on PSC2, and volumes of PSC1 and PSC2 are
    !// redistributed vertically by sedimentation.
    !// The sedimentation is done through the strato- and troposphere,
    !// by the modified routine sedimentation.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR), intent(in) :: ALT, RHO
    real(r8), dimension(LPAR+1), intent(in) :: ALTST
    real(r8), dimension(LPAR), intent(inout) :: VOLA, VICE, VOL_OLD, &
         H2Os, HNO3I, HNO3L, dH2Osed
    real(r8), intent(in) :: accNatsize, countNAT, TSTEP
    integer, intent(in) :: LMAX

    real(r8) :: &
         AVRNAT, &
         F_NAT_loc(N_B), &
         S2_NAT_loc, &
         S3_NAT_loc, &
         RNAT_loc, &
         VU_NAT, &      ! Some volumes*dr
         VU_ICE, &      ! Some volumes*dr
         VUF_ICE, &     ! Volume fraction (?), F_ICE(JI)*VOL_PAR/VU_ICE
         VUF_NAT        ! Volume fraction (?), VUF_NAT=F_NAT(JI)*VOL_PAR/VU_NAT


    !// SEDIMENTATION PARAMETERS
    real(r8), dimension(LPAR) :: &
         W, &                 ! Sedimentation velocity
         SR, SRS, SRT         ! Sedimented values

    integer :: &
         JI, &         ! Index for going through radii
         L             ! Index for going through layers
    !// --------------------------------------------------------------------
    real(r8), dimension(2), parameter :: &
         ro   = (/1.62_r8, 0.928_r8/), & ! Densities
         AMUR = (/117._r8, 18._r8/)      ! Molecular weights, HNO3*3H2O(?), H2O
    real(r8), parameter :: &
         ssatH2O = 1._r8, &        ! To get overcooled value of H2O
         PI = CPI, &               ! Pi
         PI43 = CPI * 4._r8/3._r8  ! Pi*4/3
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'sedimentation_driver'
    !//---------------------------------------------------------------------



    !// Do sedimentation only if there are PSC1a or PSC2
    !// STS is sedimented, but perhaps should not. It is done because
    !// HNO3L is sedimented, and it is not separated into NAT and STS.
    !// It is therefore sedimented as NAT.
    if (.not. (maxval(VOLA) .gt. vaer_lim .or. maxval(VICE) .gt. 0._r8)) then
       dH2Osed(:) = 0._r8
       return
    end if

    if (countNAT .gt. 0._r8) then
       !// Average radius
       AVRNAT = accNATsize / countNAT

       !// Local log-normal distribution
       !// For small AVRNAT, this may end up in zero F_NAT,RNAT,S2_NAT,
       !// and S3_NAT. This must be checked before sedimentation, where we
       !// divide by VU_NAT.
       call DIS_LN(N_B, R_A, F_NAT_loc, AVRNAT, S_NAT, &
            RNAT_loc, S2_NAT_loc, S3_NAT_loc)

    else
       !// No PSC1 (no VOLA)
       S3_NAT_loc = 0._r8
       S2_NAT_loc = 0._r8
       F_NAT_loc(:) = 0._r8
    end if


    VU_NAT = PI43 * S3_NAT_loc
    VU_ICE = PI43 * S3_ICE


    !// Loop through size bins for ICE first, then NAT/STS
    if (VU_ICE .gt. 0._r8) then

       do JI = 1, N_B
          !VOL_PAR = PI43 * R_A(JI)**3._r8 * DR(JI)

          !// ICE SEDIMENTATION
          !// --------------------------------------------------------------
          VUF_ICE = F_ICE(JI) * VOL_PAR(JI) / VU_ICE

          !// Cycle if zero (check for .le. just in case)
          if (VUF_ICE .le. 0._r8) cycle

          !// Initialize (may not be necessary)
          !SRS(:)= 0._r8
          !SRT(:)= 0._r8
          !SR(:) = 0._r8
          !W(:)  = 0._r8
          !// Only sediment from LMAX and downwards
          do L = 1, LMAX
             W(L) = WSED(R_A(JI), ALT(L), ro(2))
             if (LDEBUG) then
                if (W(L) .lt. 0._r8) then
                   write(6,'(a,2i4,4es14.6)') f90file//':'//subr// &
                        ': CRITICAL: ICE: W(L) LT 0',L,JI,W(L),ro(2), &
                        ALT(L),R_A(JI)
                   stop
                end if
             end if !// if (LDEBUG) then
             SRS(L) = H2Os(L) * VUF_ICE
             SRT(L) = HNO3I(L) * VUF_ICE
             SR(L)  = VICE(L) * VUF_ICE
          end do
          !// Do the sedimentation
          call sedimentation2(VICE,SR,ALTST,W,TSTEP,1,LMAX,1,JI)

          !// Also find sedimented H2Os
          dH2Osed(1:LMAX) = - H2Os(1:LMAX)
          call sedimentation2(H2Os,SRS,ALTST,W,TSTEP,1,LMAX,2,JI)
          dH2Osed(1:LMAX) = dH2Osed(1:LMAX) + H2Os(1:LMAX)

          call sedimentation2(HNO3I,SRT,ALTST,W,TSTEP,1,LMAX,3,JI)

       end do ! do JI = 1, N_B
    end if !// if (VU_ICE .gt. 0._r8) then


    !// STS SEDIMENTATION
    !// ----------------------------------------------------------------
    if (vu_nat .gt. 0._r8) then

       do JI = 1, N_B
          !VOL_PAR = PI43 * R_A(JI)**3._r8 * DR(JI)

          VUF_NAT = F_NAT_loc(JI) * VOL_PAR(JI) / VU_NAT

          !// Cycle if zero (check for .le. just in case)
          if (VUF_NAT .le. 0._r8) cycle

          do L = 1, LMAX
             W(L) = WSED(R_A(JI), ALT(L), RHO(L))
             if (LDEBUG) then
                if (W(L) .lt. 0._r8) then
                   write(6,'(a,2i4,4es14.6)') f90file//':'//subr// &
                        ': CRITICAL: STS: W(L) LT 0',L,JI,W(L),RHO(L), &
                        ALT(L),R_A(JI)
                   stop
                end if
             end if !// if (LDEBUG) then
             SRS(L) = HNO3L(L) * VUF_NAT
             !// Added SRT as in last SCTM1 version
             SRT(L) = VOL_OLD(L) * VUF_NAT
             SR(L)  = VOLA(L) * VUF_NAT
          end do
          !// Do the sedimentation
          call sedimentation2(VOLA,SR,ALTST,W,TSTEP,1,LMAX,4,JI)
          call sedimentation2(VOL_OLD,SRT,ALTST,W,TSTEP,1,LMAX,5,JI)
          call sedimentation2(HNO3L,SRS,ALTST,W,TSTEP,1,LMAX,6,JI)
          if (LDEBUG) then
             do L = 1, LMAX
                if (VOLA(L) .lt. 0._r8 .OR. VOLA(L) .NE. VOLA(L)) then
                   write(6,'(a,2i4,es14.6)') f90file//':'//subr// &
                        ':CRITICAL: VOLA after sed',L,JI,VOLA(L)
                   stop
                end if
             end do
          end if !// if (LDEBUG) then

       end do ! do JI = 1, N_B
    end if !// if (vu_nat .gt. 0._r8) then

    !// --------------------------------------------------------------------
  end subroutine sedimentation_driver
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine psc_box(TEMP, SAD, &                  ! in
       HNO3_gas, H2O_gas, H2SO4, &
       HNO3_iceliq, HNO3_ice, HNO3_liq, H2O_ice, &
       VAER_L, SAT, accNATsize, countNAT, &        ! inout
       VICE, VAER_L_OLD, RHO, &                    ! out
       n_particles, TNAT, TICE, &
       I,J,L)
    !// --------------------------------------------------------------------
    !// Calculates amount of PSCs, and new values for solid/liquid/gaseous
    !// HNO3 and solid H2O. Volumes of PSCs are
    !// calculated, as is the aerosol surface aera.
    !//
    !// Solid H2O is calculated from saturation value of gaseous H2O, and the
    !// solid H2O is calculated every time. If there is need of using H2O_ice in
    !// e.g. the chemistry, one have to take that into account when
    !// calculating new H2O_ice values. At this moment, H2O is calculated from
    !// AIR_MOLEC, CH4 and H2 in the stratosphere, and then modified here
    !// (each time).
    !//
    !// The surface area density (SAD) is calculated here for temperatures
    !// 205 < T < 220. SAD may be used instead of PARTAREA in calculating
    !// -PARBK constants, but one should be aware that the psps_psc routine
    !// uses PARTAREA to calculate H2SO4, so that indirectly SAD is calculated
    !// based on PARTAREA.
    !//
    !// Input:
    !//   SAD     - Aerosol surface area [cm2/cm3]
    !//   TEMP    - Temperature          [K]
    !//   H2SO4   - total H2SO4       [molec/cm3]
    !//
    !// In/out:
    !//   HNO3_gas : gaseous HNO3             [molec/cm3]
    !//   H2O_gas  : gaseous H2O              [molec/cm3]
    !//   HNO3_iceliq : solid (liq+ice) HNO3  [molec/cm3]

    !//
    !//
    !// OUTPUT:
    !//
    !//   HNO3_iceliq   - solid HNO3        [molec/cm3]
    !//             old value is overwritten
    !//   H2O_ice    - solid H2O         [molec/cm3]
    !//   PSC1    - PSC1 surface area [cm2/cm3]
    !//   PSC2    - PSC2 surface area [cm2/cm3]
    !// --------------------------------------------------------------------
    use utilities_oslo, only: H2O_SAT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: &
         I,J,L             ! The grid indices

    real(r8), intent(in) :: &
         H2SO4, &        ! H2SO4 input
         TEMP, &         ! Temperature
         SAD             ! Area of aerosols

    !// Input/output
    real(r8), intent(inout) :: &
         HNO3_gas, &     ! HNO3 gaseous
         HNO3_iceliq, &  ! HNO3 solid, will be altered
         H2O_gas, &      ! H2O gaseous
         H2O_ice, &      ! H2O solid, not input yet, only output
         VAER_L, &       ! Volume of aerosols, will be altered
         accNATsize, &   ! Accumulated NAT size, may be increased
         countNAT        ! Accumulated count of NATs, may be increased

    logical, intent(inout) :: &
         SAT             ! True if STS has frozen to NAT

    !// Output
    real(r8), intent(out) :: &
         HNO3_ice, &     ! Formed solid HNO3, ice
         HNO3_liq, &     ! Formed solid HNO3, liquid
         VICE, &         ! 
         VAER_L_OLD, &   !
         RHO, &          ! Density of STS
         n_particles, &
         TNAT, &         ! TNAT; where NAT can start forming.
         TICE            ! TICE; where PSC2 forms


    !// Locals
    real(r8) :: &
         VICE2           ! Volume of HNO3_ice, the coating of PSC2 (scalar)


    !// Other variables
    real(r8) :: &
         TMELT, TSAT, &  ! Temperatures giving melting or freezing.
         HNO3_ambient, & ! Total HNO3 (solid + gaseous)
         H2O_ambient, &  ! Total H2O (solid + gaseous)
         HNO3torr, &     ! Pressure of HNO3
         H2Osat, &       ! Saturated H2O
         H2Otorr, &      ! Pressure of H2O
         H2Oss, &        ! Overcooled value of H2O
         AMP0, &         ! k*T (Boltzman's constant * Temperature)
         AMP, &          ! k*T*H2O
         r_particle, &
         WS, &           ! Weight fraction of H2SO4 in the aerosol, calculated
                         ! in subroutine CARS, but not used here.
         WN              ! same, but for HNO3


    !// --------------------------------------------------------------------
    real(r8), dimension(2), parameter :: &
         ro   = (/1.62_r8, 0.928_r8/), & ! Densities
         AMUR = (/117._r8, 18._r8/)      ! Molecular weights, HNO3*3H2O(?), H2O
    real(r8), parameter :: &
         ssatH2O = 1._r8, &        ! To get overcooled value of H2O
         PI = CPI, &               ! Pi
         PI43 = CPI * 4._r8/3._r8  ! Pi*4/3
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'psc_box'
    !//---------------------------------------------------------------------

    !// Routine is called only for stratospheric air. Calling routine should
    !// do initialisations.


    !// --------------------------------------------------------------------
    !// PARTICLE SIZE BINS
    !// --------------------------------------------------------------------
    !// They are already set up. If subject to change, they
    !// are set further below.

    !// --------------------------------------------------------------------
    !// MOMENTUMS OF LOG-NORMAL DISTRIBUTION
    !// --------------------------------------------------------------------
    !// Global/fixed values are already set up. If subject to change, they
    !// are set further below.

    !// --------------------------------------------------------------------
    !// CALCULATE PSC/SAD FORMATION / MELTING in stratosphere
    !// AND EVENTUALLY MELTING in troposphere
    !// --------------------------------------------------------------------
    !// Comment by Ole Amund Sovde, 15. November 2004
    !// Use L=1,LM instead of L=LOWERL,LM for melting and for sedimentation,
    !// but calculate formation from LOWERL and upwards only.
    !// --------------------------------------------------------------------
    HNO3_ambient = HNO3_gas + HNO3_iceliq   ! Total HNO3 [molec/cm3]
    H2O_ambient = H2O_gas + H2O_ice          ! Total H2O [molec/cm3]
    H2Osat  = H2O_SAT(TEMP)                ! SATURATED H2O
    H2Oss   = ssatH2O * H2Osat             ! OVERCOOLED VALUE OF H2O
    !// H2O: find solid H2O
    H2O_ice = MAX(0._r8, H2O_ambient - H2Oss)  ! SOLID H2O (ice)
    H2O_gas = H2O_ambient - H2O_ice               ! WATER VAPOR

    AMP0 = KBOLTZ * Temp                   ! k*T [mbcm3/molec]
    AMP  = AMP0 * H2O_gas                     ! k*T*H2O (mb)

    !// 760Torr=1013.25mb, so we multiply with 760./1013.25
    H2Otorr  = AMP * 760._r8/1013.25_r8                ! pressure of H2O
    HNO3torr = HNO3_ambient * AMP0 * 760._r8/1013.25_r8 ! pressure of HNO3

    !// Temperature for STS freezing. Found no reference, although I suspect
    !// it is freezing of one of the H2O/H2SO4/HNO3 mixtures described e.g.
    !// by Fox et al, Science, vol 267, pp 351, 1995, where
    !//   log10(K) = A/T + B
    !// I think this is close to the MixH:SA4 line in their Figure 4,
    !// where SAT (SA4) is in equlibrium with a solution of H2SO4/HNO3/H2O
    !// stochiometrically: H2SO4x4H2O in equlibrium with
    !// H2SO4xHNO3x5H2O (i.e.
    !// SA4 with one extra HNO3 and H2O). This means n=2 in their data.
    !// Then we have k=4 (SA4), and l=5 from the MixH, and combined they
    !// give the equlibirium K(l-k) = pHNO3*pH2O.
    !// According to Fox et al, the values should be
    !//    A = -9229.19162 and B = 38.06445
    !// Close, and maybe a cigar; the calculated temperature is slightly
    !// higher; 0 - 0.5K, lowest departure at low temperatures.
    !// Reading off the figure may get you closer to the values used here:
    TSAT    = -9384._r8 / (log10(HNO3torr * H2Otorr) - 39._r8)

    !// Temperature for SAT melting:
    !// Zhang, R. et al.,J. Phys. Chem. Vol 97, pp 7351, 1993
    !//   Physical Chemistry of the H2SO4 binary system at Low
    !//   Temperatures: Stratospheric Implications
    TMELT   = 3236._r8 / (11.502_r8 - log10(H2Otorr))

    !// Find TNAT (from Hanson & Mauersberger 1988)
    TNAT = getTnat(HNO3torr, H2Otorr)
    !// TNAT can also be calculated from (no reference available)
    !//   TNAT(L) = -12414.3_r8/(DLOG10(HNO3torr)
    !//             +3._r8*LOG10(H2Otorr) - 45.68_r8)

    !// TICE, temperature for freezing,
    !// from Marti and Mauersberger, GRL 20, No 5, 363-366, 1993
    TICE = 2663.5_r8 / (12.537_r8 - log10(AMP * 100._r8))
    !// As in CARS, assume a supersaturation corresponding to 3K
    !// (See Considine et al., JGR 105, D3, 3955-3973, 2000)
    !// and a minimum of TICE_MIN (187K).
    TICE = max(TICE - 3._r8, TICE_MIN)

    !// We always have TSAT < TMELT, and so far I have seen TSAT < TNAT,
    !// so that once TSAT is reached, STS freezes to NAT, and we keep
    !// SAT=.TRUE. until T > TMELT
    if (TEMP .le. TSAT)  SAT = .TRUE.
    if (TEMP .gt. TMELT) SAT = .FALSE.

    !// Number of particles;
    n_particles  = SAD / (4._r8 * PI * S2_SA) ! Number of SAD particles
    !   PARL(L) = PSC1(L)/(4._r8*PI*S2_NAT)    ! Number of PSC1 particles
    !   PARI(L) = PSC2(L)/(4._r8*PI*S2_ICE)    ! Number of PSC2 particles

    if (.not.SAT) then
       !// Temperature is too high for solid PSCs. We may still get some STS.
       !//   PAR(L) = PAR(L)+PARL(L)+PARI(L)
       HNO3_gas    = HNO3_ambient ! All is gaseous
       HNO3_iceliq = 0._r8        ! No solid, but may be formed as STS
       VAER_L_OLD  = 0._r8        ! Volume is set to zero
    else
       !// Temperature is low enough for STS to freeze. Already formed
       !// STS freeze, otherwise what is formed now will freeze.
       !// Keep the incoming gaseous HNO3, incoming solid HNO3 stays solid
       VAER_L_OLD = VAER_L       ! Volume of already frozen particles
       !TNAT(L)=TNAT(L) + 3._r8   ! Not used
    end if

    !// Calculate PSC formation in the stratosphere
    !// --------------------------------------------------------------

    ! Initialize weight fraction of H2SO4 and HNO3 in aerosol
    WS = 0._r8
    WN = 0._r8

    !// WATER VAPOR CAN NOT EXCEED ITS SATURATED VALUE
    !// --------------------------------------------------------------
    if (LPSC2) then
       VICE = AMUR(2) * H2O_ice / (AVOGNR * ro(2)) ! Volume of ice
       !+4.*PI*PARI(L)*S3_ICE/3.   
    else
       !// Do not change H2O_gas, keep input value until CARS. After
       !// CARS, H2O_gas is set to H2O_ambient.
       VICE = 0._r8
    end if
    VICE2   = 0._r8

    !// TERNARY SOLUTION H2SO4/HNO3/H2O
    !// --------------------------------------------------------------
    ! VAER_L  - VOLUME OF TERNARY SOLUTION
    ! HNO3_liq - LIQIUD HNO3
    ! HNO3_ice - SOLID HNO3
    call CARS(Temp, H2O_gas, HNO3_gas, H2SO4, VAER_L, &
               RHO, HNO3_liq, HNO3_ice, VICE2, I,J,L)

    if (LDEBUG) then
       if (VAER_L .lt. 0._r8 .or. VAER_L.ne.VAER_L) &
            write(6,'(a,i5,10es12.3)') f90file//':'//subr// &
            ': VAER_L wrong after cars',L,VAER_L,RHO, &
            VAER_L_OLD, HNO3_liq, HNO3_ice, HNO3_gas, H2SO4, &
            Temp,VICE2
    end if

    !// --------------------------------------------------------------
    !// SO FAR
    !// We have now new values for gaseous species, and the amount
    !// of liquid/ice HNO3 formed from the (old) gaseous phase.
    !// VOLA(L) is the volume of these new particles, so we need to
    !// add this to the already existing volume, VOL_OLD(L).
    !// If all is assumed evaporated, VOLA(L)=0._r8, VOL_OLD(L)=0._r8
    !// --------------------------------------------------------------



    !// VOLA AS IN SCTM1
    !// --------------------------------------------------------------
    !// Although, in total, I do not understand this way of seeing
    !// the volume growth. See comment from Sergei at the bottom
    !// of this routine.
    !//   VOLA(L)=MAX1(VOLA(L),VOL_OLD(L))
    VAER_L = VAER_L + VAER_L_OLD
    if (LPSC2) then
       !// Add the ice HNO3 (HNO3 coating)
       VICE = VICE + VICE2
    else
       !// Leave out the PSC2; both H2O and ice HNO3
       VICE = 0._r8  ! This one is actually set at the top
       HNO3_gas = HNO3_gas + HNO3_ice
       HNO3_ice = 0._r8
       H2O_gas = H2O_ambient
       H2O_ice = 0._r8
    end if

    !// Allow NAT to grow (Calculate NAT size, sps07)
    if (VAER_L .le. vaer_lim) then
       VAER_L = 0._r8
    else
       !// If there were no particles before formation, PAR will be 0.
       !// To avoid division by 0, PAR is set to 2 when this happens.
       if (n_particles .eq. 0._r8) n_particles = 2._r8
       r_particle = (VAER_L / (PI43 * n_particles))**(1._r8 / 3._r8) * 0.595238_r8
       accNATsize = accNATsize + r_particle
       countNAT = countNAT + 1._r8
    end if


    !// PRESUMEABLY UNNECESSARY TESTS
    !// -----------------------------------------------------------------
    if (LDEBUG) then
       !// Do some testing in case of negative values. I don't
       !// think this ever happens, but it is a part of the testing.
       if (HNO3_liq .lt. 0._r8) then
          write(6,'(a,3i4,5es11.3,f7.2)') &
               f90file//':'//subr//': HNO3_liq<0! (-> 1e-20)',&
               I,J,L, HNO3_liq, HNO3_gas, &
               HNO3_ice, HNO3_iceliq, HNO3_ambient, TEMP
          stop
          HNO3_liq = 1.e-20_r8
       end if
       if (HNO3_iceliq .lt. 0._r8) then
          write(6,'(a,3i4,5es11.3,f7.2)') &
               f90file//':'//subr//': HNO3_iceliq<0! Should never happen',&
               I,J,L, HNO3_iceliq, HNO3_gas, &
               HNO3_ice, HNO3_liq, HNO3_ambient, TEMP
          stop
       end if
       if (HNO3_ice .lt. 0._r8) then
          write(6,'(a,3i4,5es11.3,f7.2)') &
               f90file//':'//subr//': HNO3_ice<0! (-> 1e-20)',&
               I,J,L,HNO3_ice, HNO3_gas, &
               HNO3_liq, HNO3_iceliq, HNO3_ambient, TEMP
          stop
          HNO3_ice = 1.e-20_r8
       end if
    end if !// if (LDEBUG) then



    !// LIQUID HNO3
    !// -----------------------------------------------------------------
    !// Add the old solid particles to the new liquid part
    HNO3_liq = HNO3_liq + HNO3_iceliq

    !// GASEOUS HNO3
    !// -----------------------------------------------------------------
    HNO3_gas = HNO3_ambient - HNO3_liq - HNO3_ice

    !// Testing the gaseous value.
    !// For small HNO3AMBIENT, GHNO3 may become slightly negative...
    if (HNO3_gas .lt. 0._r8) then
       write(6,'(a,3i4,5es11.3,f7.2)') &
            f90file//':'//subr//': HNO3_gas<0! (-> 1e-20)',&
            I,J,L,HNO3_gas, &
            HNO3_liq, HNO3_ice, HNO3_iceliq, HNO3_ambient, TEMP
       HNO3_gas = 1.e-20_r8
    end if


    !// --------------------------------------------------------------------
  end subroutine psc_box
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine get_psc12_sad(I, J, P1, P2, SAD)
    !// --------------------------------------------------------------------
    !// Retrieve column values of PSC1 and PSC2.
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I, J
    !// Output
    real(r8), intent(out) :: P1(LPAR), P2(LPAR), SAD(LPAR)
    !// --------------------------------------------------------------------
    !// Retrieve column values of PSC1
    P1(:) = PSC1(:,I,J)
    !// Retrieve column values of PSC2
    P2(:) = PSC2(:,I,J)
    !// Retrieve column values of sps_PARTAREA (i.e. background SAD)
    SAD(:) = sps_PARTAREA(:,I,J)
    !// --------------------------------------------------------------------
  end subroutine get_psc12_sad
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_mw(trsp_idx, Xtrsp_idx, TMASS, XTMASS, TRACER_ID_MAX, &
       NPAR, NOTRPAR, LPAR, Z_HNO3, Z_NOy, Z_HNO3s, Z_H2Os)
    !// --------------------------------------------------------------------
    !// Fetches coefficients for mass to concentration conversion. Checks both
    !// trsp_idx and Xtrsp_idx.
    !// 
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Array sizes
    integer, intent(in) :: TRACER_ID_MAX, NPAR, NOTRPAR, LPAR
    !// Transport number and non-transport number
    integer, dimension(TRACER_ID_MAX), intent(in) :: trsp_idx, Xtrsp_idx
    !// Tracer data:   mol.wt.
    real(r8), dimension(NPAR), intent(in) :: TMASS
    !// Tracer data:   mol.wt.
    real(r8), dimension(NOTRPAR), intent(in) ::  XTMASS
    !// Output
    real(r8), intent(out) :: Z_HNO3, Z_NOy, Z_HNO3s, Z_H2Os
    !// --------------------------------------------------------------------

    !// Get molecular masses for transforming to concentration
    if (trsp_idx(4).gt.0) then
       Z_HNO3 = AVOGNR * 1.e-3_r8 / TMASS(trsp_idx(4))
    else
       if (Xtrsp_idx(4).gt.0) then
          Z_HNO3 = AVOGNR * 1.e-3_r8 / XTMASS(Xtrsp_idx(4))
       else
          write(*,'(a)') 'Microphysics without HNO3!'
          stop
       end if
    end if
    if (trsp_idx(147).gt.0) then
       Z_NOy = AVOGNR * 1.e-3_r8 / TMASS(trsp_idx(147))
    else
       if (Xtrsp_idx(147).gt.0) then
          Z_NOy = AVOGNR * 1.e-3_r8 / XTMASS(Xtrsp_idx(147))
       else
          write(*,'(a)') 'Microphysics without NOy_str!'
          stop
       end if
    end if
    if (trsp_idx(124).gt.0) then
       Z_HNO3s = AVOGNR * 1.e-3_r8 / TMASS(trsp_idx(124))
    else
       if (Xtrsp_idx(124).gt.0) then
          Z_HNO3s = AVOGNR * 1.e-3_r8 / XTMASS(Xtrsp_idx(124))
       else
          write(*,'(a)') 'Microphysics without HNO3s!'
          stop
       end if
    end if
    if (trsp_idx(125).gt.0) then
       Z_H2Os = AVOGNR * 1.e-3_r8 / TMASS(trsp_idx(125))
    else
       if (Xtrsp_idx(125).gt.0) then
          Z_H2Os = AVOGNR * 1.e-3_r8 / XTMASS(Xtrsp_idx(125))
       else
          write(*,'(a)') 'Microphysics without H2Os!'
          stop
       end if
    end if

    !// --------------------------------------------------------------------
  end subroutine get_mw
  !// ----------------------------------------------------------------------






  !// ----------------------------------------------------------------------
  subroutine CARS(T,H2O,HNO3,H2SO4,VA,RHO,HNO3A,HNO3I,VI, I,J,L)
    !// --------------------------------------------------------------------
    !// CARSLAW ET AL., GRL,22,1877-1880,1995.
    !//
    !// INPUT:  
    !//   T     - TEMPERATURE (K)
    !//   H2O   - TOTAL H2O,   CM(-3)
    !//   HNO3  - TOTAL HNO3,  CM(-3)
    !//   H2SO4 - TOTAL H2SO4, CM(-3)
    !//
    !// OUTPUT:
    !//   VA    - VOLUME OF AEROSOL PER UNIT VOLUME OF AIR, 
    !//   RHO   - TERNARY SOLUTION MASS DENSITY, KG/M(3)
    !//   PN    - EQUILIBRIUM VAPOUR PRESSURE OF HNO3 OVER THE DROPLET, ATM
    !//   HNO3A - LIQUID HNO3, CM(-3)
    !//   HNO3I - SOLID HNO3, CM(-3)
    !//
    !//   VI    - VOLUME OF ICE HNO3 (T < TMIN)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: T, H2O, HNO3, H2SO4
    integer, intent(in) :: I,J,L
    !// Output
    real(r8), intent(out):: VA, RHO, HNO3A, HNO3I, VI

    !// Locals
    real(r8) :: &
         AMNB,AMSB,ANB,ASB,TT,ANBT,ASBT,AMSB2,AMNB2,AMNB_AMSB, &
         A2B3,A2C, &
         A,B,C, &! constants
         PN0,PN,PHI,PI,AMS,AMN,RHOS,RHON
    real(r8) :: &
         lphno3s, &  ! log10(saturation pressure of hno3)
         pHNO3s, &   ! saturation pressure of hno3
         HNO3s, &    ! saturated conc. of hno3
         AMP0,AMP, & ! kB*T, kB*T*H2O
         pH2O, &     ! partial pressure of H2O
         PH, &       ! Molarity of H2SO4
         HNSB,ZNAM,HNNB, & ! variables
         AM, &
         WS, &       ! Weight fraction of H2SO4 in the aerosol
         WN, &       ! Weight fraction of HNO3 in the aerosol
         PW,PW_MB, & ! partial pressure H2O in atm and mb
         TMIN,TICE   ! Temperatures for minimum freezing temperature
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CARS'
    !//---------------------------------------------------------------------

    ! Initializing values (they might not be initialized otherwise)
    PN0 = HNO3 * R_ATM / AVOGNR * T * 1.e6_r8 ! HNO3 partial pressure [atm]
    PH  = H2SO4 * 1.e6_r8 / AVOGNR            ! H2SO4 [mol/m3]
    RHO = 1.62_r8                  ! Initialize RHO
    VI    = 0._r8
    VA    = 0._r8
    TMIN = 0._r8

    !// WATER PARTIAL PRESSURE
    !// ----------------------
    PW    = H2O * R_ATM / AVOGNR * T * 1.e6_r8  ! pressure [atm]
    PW_MB = PW * 1013._r8            ! pressure [mb]

    !// VALID ONLY FOR 2.E-5MB < PW < 2.E-3MB
    !//            MAX(TICE-3K,185K) < T < 240K
    !//              STS-caclulation only up to 215K.
    !// --------------------------------------------------------------------
    if (PW_MB .gt. 2.e-5_r8 .and. PW_MB .lt. 2.e-3_r8 &
         .and. T .le. 240._r8) GOTO 3
    !// Temperature is too high for PSC, set values to zero and return.
    !// --------------------------------------------------------------------
8   VA    = 0._r8    ! This is if temperature is too high for PSC
    VI    = 0._r8    !   : ice volume
    HNO3A = 0._r8    !   : liquid
    HNO3I = 0._r8    !   : ice
    GOTO 2          ! And then go to return

    !// Temperature is .LE. 240K
    !// ------------------------
    !// Calculate TICE from Marti and Mauersberger,
    !// GRL Vol 20, No 5, 363-366, 1993
3   TICE = 2663.5_r8 / (12.537_r8 - log10(KBOLTZ * T * H2O * 100._r8))
    !// Assume a supersaturation corresponding to 3K
    !// (See Considine et al., JGR 105, D3, 3955-3973, 2000)
    TMIN = max(TICE - 3._r8, TICE_MIN)    ! MIN of ice temp and 187K

    if (T .LT. 215._r8) GOTO 9             ! It's colder than 215K
    !// Carslaw et al, 1995: No ternary solution above 215K:
    !// 215K < T <= 240K, which means no frozen or liquid HNO3
    !// Only binary solution of H2SO4/H2O is possible
    !// -------------------------------------------------------
    HNO3I = 0._r8                          ! NO FROZEN HNO3
10  AMS   = AM_BIN(T,PW,AK_CARS(1,2))     ! (KG H2O)(-1)
    AMN   = 0._r8                          ! No HNO3 in solution
    PN    = PN0                           ! No HNO3 in solution
    !// We may still have H2SO4 in binary solution with H2O;
    !// The volume of sulphate aerosols are calculated at nr 7.
    GOTO 7

    !// Here temperature is lower than 220K, but how low is it? If it's
    !// not low enough to form HNO3 ice as PSC2 coating, go to 1
    !// ---------------------------------------------------------------
9   if (T .ge. TMIN) GOTO 1           ! Not cold enough to get ice.

    !// T < TMIN: We have HNO3 as PSC2 coating, as HNO3 ice
    !// ---------------------------------------------------
    !// Saturation pressure of HNO3 is calculated according to
    !// Hanson and Mauersberger, GRL Vol 15, No 8, 855-858,1988
6   AM = -2.7836_r8 - 0.00088_r8 * T
    B  = 38.9855_r8 - 11397._r8 / T + 0.009179_r8 * T
    !// This was the original, one coefficient is wrong:
    !// B  = 36.1047d0 - 11397.d0/T + 0.009179d0*T
    !// We also might get here from further down in the code, where
    !// we have liquid HNO3/H2O binary solution
    AMP0  = KBOLTZ * T                       ! [mbcm3/molec]
    AMP   = AMP0 * H2O                       ! [mb]
    !!!pH2O  = AMP/1013.25_r8                  ! [atm], should be torr
    !!!pHNO3s= AM*(LOG10(pH2O)+2.8808_r8)+B              ! sat.pres. of HNO3
    !!!HNO3s = MIN(HNO3,10._r8**(pHNO3S)/AMP0*1013.25_r8) ! SATURATED HNO3
    !!!According to Hanson and Mauersberger, pressure should be in Torr.
    !!!which gives: log10(pTorr) = log10(pAtm) + 2.8808
    !!!Still, we have here that pHNO3s=log10(pressure given in Torr).
    !!!To get molec/cm3, we have to convert to mb to use AMP0:
    !!!HNO3s = (10._r8**pHNO3s)/(760Torr/atm)*(1013.25mb/atm)/AMP0
    !!!A more understandable code would be:

    pH2O   = AMP / 1013.25_r8 * 760._r8  ! [mb]/[mb/atm]*[Torr/atm] = [Torr]
    lpHNO3s= AM * (log10(pH2O)) + B      ! log10(HNO3 sat.pres.)
                                         !     = log10(pHNO3s/Torr)
    pHNO3s = 10._r8**(lpHNO3s)           ! HNO3 sat. pres. [Torr]
    HNO3s  = min(HNO3, &
         pHNO3s * 1013.25_r8 / (760._r8 * AMP0)) ! SATURATED HNO3
                      ! [Torr]*[mb/atm]/[Torr/atm]/[mbcm3/molec]=[molec/cm3]
    if (T .lt. TMIN) then
       ! T < TMIN :: PSC2 coating, HNO3 ice
       HNO3A = 0._r8                   ! ALL HNO3 IS FROZEN
       HNO3I = HNO3 - hno3s            ! SOLID HNO3
                                         ! (total HNO3 - saturated HNO3)
       !VA=117._r8*HNO3I/(AVOGNR*0.928_r8)  ! VOLUME OF SOLID HNO3
       VI = 117._r8 * HNO3I / (AVOGNR * 0.928_r8)  ! VOLUME OF SOLID HNO3
       ! PSC1 will not be formed now.
       VA = 0._r8                         ! NO VOLUME OF AEROSOLS
    else
       ! T > TMIN :: we have binary solution of HNO3/H2O
       HNO3I = 0._r8                   ! NO FROZEN HNO3
       HNO3A = HNO3 - HNO3s            ! SOLID (LIQUID) HNO3
       RHO   = 1.62_r8                 ! DENSITY OF SOLID HNO3
       VA    = 117._r8 * HNO3A / (AVOGNR * RHO) ! VOLUME OF SOLID (LIQUID) HNO3
       VI    = 0._r8
    end if

    !// Finnished, all is ice or liquid, go to return.
    GOTO 2

1   HNO3I = 0._r8
    VI    = 0._r8
    !// TMIN <= T <= 220K: no HNO3 ice.
    !// Check for binary/ternary solutions of H2SO4/HNO3 with H2O?
    !// -------------------------------------------------------------

    !// CONCENTRATION OF H2SO4 IN THE BINARY SOLUTION WITH H2O
    !// ------------------------------------------------------
    AMSB = AM_BIN(T, PW, AK_CARS(1,2)) ! (KG H2O)(-1)

    !// H2SO4 MOL/M(3)
    !// --------------
    PH = H2SO4 * 1.e6_r8 / AVOGNR
    IF (PH) 6,6,5
    !// If we have H2SO4, we continue to calculate binary solution,
    !// otherwise we go to nr 6 to calculate liquid HNO3 in binary
    !// solution.

5   continue

    !// H2SO4 PARTIAL PRESSURE
    !// ----------------------
    TT = R_ATM * T * PH

    !// HNO3 PARTIAL PRESSURE
    !// ---------------------
    PN0 = HNO3 * R_ATM / AVOGNR * T * 1.e6_r8      ! ATM

    !// HENRY LAW CONSTANTS
    !// -------------------
    HNNB = HENRIC(T, PW, Q_CARS(1,1))  ! MOL*KG(-1)*ATM(-1)
    HNSB = HENRIC(T, PW, Q_CARS(1,2))  ! MOL*KG(-1)*ATM(-1)

    !// CONCENTRATION OF HNO3 IN THE BINARY SOLUTION WITH H2O
    !// -----------------------------------------------------
    AMNB = AM_BIN(T, PW, AK_CARS(1,1)) ! (KG H2O)(-1)

    !// If no HNO3/H2O binary solution, check if there is H2SO4/H2O solution.
    if (AMNB .le. 0._r8) GOTO 10
    !// We have HNO3/H2O, what about H2SO4?
    ANB = HNNB * AMNB
    ASB = HNSB * AMSB
    ANBT = TT * ANB
    ASBT = TT * ASB
    AMSB2 = AMSB * AMSB
    AMNB2 = AMNB * AMNB
    AMNB_AMSB = AMNB - AMSB
    A = ((ANBT - ASBT) * AMNB - 2._r8 * AMNB2 * AMSB + AMNB * AMSB2 &
         + (ANB - ASB) * AMSB * PN0) / (AMNB2 - AMNB * AMSB)
    B = AMSB * (-2._r8 * ANBT + ASBT + AMSB * AMNB - HNNB * AMSB * PN0) &
         / AMNB_AMSB
    C = ANBT * AMSB2 / AMNB_AMSB
    A2B3 = A * A - 3._r8 * B
    A2C  = -2._r8*A*A*A + 9._r8*A*B - 27._r8*C
    PHI  = atan(sqrt(4._r8 * A2B3 * A2B3 * A2B3 - A2C * A2C) / A2C)
    PI = 2._r8 * ASIN(1._r8)
    if (PHI .lt. 0._r8) PHI = PHI + PI
    AMS = -(A + 2._r8 * sqrt(A2B3) * COS((PI + PHI) / 3._r8)) / 3._r8
    !// We have no H2SO4: calculate liquid HNO3 at nr 6
    if (AMS .le. 0._r8) GOTO 6

    !// We have ternary solution HNO3/H2SO4/H2O
    AMN  = MAX(0._r8, AMNB * (1._r8 - AMS / AMSB))
    PN   = AMN * (AMN + AMS) / (HNNB * AMN + HNSB * AMS)
7   ZNAM = 1._r8 / (1._r8 + AMS * 0.09812_r8 + AMN * 0.06303_r8)
    WS   = AMS * 0.09812_r8 * ZNAM
    WN   = AMN * 0.06303_r8 * ZNAM
    RHOS = RHO_S_CAR(T, AMS)  ! density of binary H2SO4/H2O
    RHON = RHO_N_CAR(T, AMN)  ! density of binary HNO3/H2O
    RHO  = 1.e-3_r8 / ((1._r8 / RHOS) * AMS / (AMN + AMS) &
         + (1._r8 / RHON) * AMN / (AMN + AMS))
    ! Added test to find possible bug
    if (RHO .ne. RHO) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': RHO is NAN ',I,J,L
       write(6,'(a,3es12.3)') ': RHO/RHOS/RHON',RHO,RHOS,RHON
       write(6,'(a,3es12.3)') ': AMS/AMN',AMS,AMN
       write(6,'(a,3es12.3)') ': T/TMIN/TICE',T,TMIN,TICE
       write(6,'(a,3es12.3)') ': PN/PN0/PH',PN0,PN, PH
       write(6,'(a,3es12.3)') ': HNSB/HNNB',HNSB,HNNB
       write(6,'(a,3es12.3)') ': AMNB/AMSB',AMNB,AMSB
       write(6,'(a,3es12.3)') ': A/B/C',A,B,C
       write(6,'(a,3es12.3)') ': WS/WN',WS,WN
       write(6,'(a,3es12.3)') ': H2O',H2O
       stop
       print*,'RESCUE: RHO=1.62'
       RHO = 1.62_r8
    end if
    VA   = PH * 0.09812e-3_r8 / (WS * RHO)
    if (.false.) then
       !// Added test to find possible bug
       if (VA.ne.VA) then
          print*,'WARNING: VA in CARS',VA,PH,WS,RHO,RHON,RHOS,AMN,AMS, &
               PN,PN0,T,HNSB,HNNB,AMNB,AMSB,WN
          print*,'RESCUE: VA=0.0'
          stop
          !VA=0.0
          !PN=PN0
       end if
    end if
    HNO3A = max(0._r8,(PN0 - PN) * AVOGNR / (R_ATM * T * 1.e6_r8))
2   return
    !// --------------------------------------------------------------------
  end subroutine CARS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function WSED(Rcm,Z,RO)
    !// --------------------------------------------------------------------
    !// FALLING SPEED OF AEROSOL PARTICLES
    !// F.KASTEN, J.APPL.MET.,7,944-947,1968.
    !//
    !// March 2007: Modified for psps_psc module in the Oslo CTM2.
    !//
    !// INPUT:
    !//   Rcm  - RADIUS OF PARTICLE, cm
    !//   Z    - ALTITUDE, km
    !//   RO   - MASS DENSITY OF AEROSOL PARTICLES, gm/cm(3)
    !// OUTPUT:
    !//   WSED - FALLING SPEED, cm/sec
    !//
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: Rcm,Z,RO
    !// Locals
    integer :: II,KK, IgotNN
    !real(r8) :: ZF(wsed_LL)
    real(r8) :: ZXP,ZXM,YXP,YXM,ZN,R
    !// --------------------------------------------------------------------

    R = Rcm * 1.e4_r8 ! convert to microns
    II = 1
    do while(wsed_ZF(II) + 0.1_r8 .lt. Z)
       II = II + 1
       if (II .gt. wsed_LL) then !GOTO 2
          !// Too high for calculation
          WSED = 0._r8
          return
       end if
    end do
    if (II .eq. 1) II = II + 1
    ZXP = wsed_ZF(II) - Z
    ZXM = Z - wsed_ZF(II-1)
    KK = 1
    IgotNN = 0
    do while (wsed_RF(KK) + 0.0001_r8 .lt. R)
       KK = KK + 1
       if (KK .ge. wsed_NN) then !GOTO 4
          !// Larger than given in look-up table
          IgotNN = 1
          exit
       end if
    end do

    !// Regular calculation
    if (IgotNN .eq. 0) then
       if (KK .eq. 1) KK = KK + 1
       YXP = (wsed_RF(KK) - R) * wsed_AMN(KK - 1)
       YXM = (R - wsed_RF(KK - 1)) * wsed_AMN(KK)
       ZN  = 1._r8/((wsed_ZF(II) - wsed_ZF(II-1))*(wsed_RF(KK) - wsed_RF(KK-1)))
       WSED= ((wsed_WF(II-1, KK-1) * ZXP + wsed_WF(II,KK-1) * ZXM) * YXP &
            + (wsed_WF(II-1, KK) * ZXP + wsed_WF(II,KK) * ZXM) * YXM) * ZN
    else
       !// KK >= wset_NN
       YXP = log10(wsed_RF(KK) / R)
       YXM = log10(R / wsed_RF(KK-1))
       ZN  = 1._r8 / ((wsed_ZF(II) - wsed_ZF(II-1)) &
            * (log10(wsed_RF(KK) / wsed_RF(KK-1))))
       WSED= 10._r8**(((log10(wsed_WF(II-1, KK-1) * wsed_AMN(KK-1)) * ZXP &
            + log10(wsed_WF(II,KK-1) * wsed_AMN(KK-1)) * ZXM) * YXP &
            + (log10(wsed_WF(II-1,KK) * wsed_AMN(KK))*ZXP &
            + log10(wsed_WF(II,KK) * wsed_AMN(KK)) * ZXM) * YXM) * ZN)
    end if
    WSED = RO * WSED
    !// --------------------------------------------------------------------
  end function WSED
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine sedimentation2(CONC,SED,Z,W,DT,LOWERL,UPPERL,NR,NB)
    !// --------------------------------------------------------------------
    !// Based on sedimentation routine from Sergei, modified to better fit
    !// the Oslo CTM2.
    !//
    !// Sedimentation of an atmospheric particle with concentration CONC
    !// in solid state with falling speed W. The sedimentation starts at
    !// LOWERL and then each layer above is treated. Sedimentation from
    !// LOWERL (to LOWERL-1) is assumed to rain out. In the Oslo CTM2
    !// LOWERL is 1, so that stratospheric values are sedimented into the
    !// troposphere (where they may eventually evaporate).
    !//
    !// Input:
    !//        LOWERL    - The tropopause level (start of stratosphere)
    !//        Z(LPAR+1)   - Straggled altitude values, cm
    !//        W         - Falling speed, cm/sec
    !//        DT        - Time step, sec
    !//        SED(LPAR)   - Sedimented concentration, cm(-3)
    !//        (Just for testing:)
    !//        NR        - Number of the called sedimentation (see PSC_1d)
    !//        NB        - The size bin
    !//
    !// Output:
    !//        CONC(L)   - Implicit concentration, cm(-3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: &
         LOWERL, &           ! The lowermost box for sedimentation
         UPPERL, &           ! The uppermost box for sedimentation
         NR, NB              ! Test variables
    real(r8), intent(in) :: &
         Z(LPAR+1), &        ! Straggled altitude values [cm]
         W(LPAR), &          ! Falling speed [cm/sec]
         DT, &               ! Time step [s]
         SED(LPAR)           ! Sedimented concentration [molec/cm3]
    real(r8), intent(inout) :: &
         CONC(LPAR)          ! Implicit concentration [molec/cm3]

    !// Local variables
    integer :: &
         L, L2               ! Looping indices
    real(r8) :: &
         E, &                ! Fraction of box sedimented
         HSED                ! The amount sedimented from the box
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'sedimentation2'
    !//---------------------------------------------------------------------

    do L = LOWERL, UPPERL
       L2 = L
       ! It might be a good idea to check if CONC .GT. 0.d0 before
       ! doing the sedimentation.
       if (CONC(L) .gt. 0._r8) then

          ! First find the fraction of the box that is sedimented
          if (Z(L+1) .gt. Z(L)) then
             E = W(L) * DT / (Z(L+1) - Z(L))
          else
             write(6,'(a,3i4,2es16.6)') f90file//':'//subr// &
                  ': Seriously wrong: Z(L+1) <= Z(L): L/NR/NB/Z(L+1),Z(L)',&
                  L,NR,NB,Z(L+1),Z(L)
             stop
             !E = 0._r8
          end if

          ! Test if we sediment more than we have
          if (SED(L) .gt. CONC(L)) then
             write(6,'(a,3i4,2es16.6)') f90file//':'//subr// &
                  ': TOO MUCH SEDIMENTATION', &
                  L, NR, NB, SED(L), CONC(L)
             !SED(L) = CONC(L)
          end if

          HSED = min(SED(L), CONC(L))

          if (E .lt. 0._r8) then
             write(6,'(a,2i4,5es16.6)') f90file//':'//subr// &
                  ': E IS NEGATIVE!',L,L2,E,W(L),(Z(L+1) - Z(L)),Z(L+1),Z(L)
             stop
          end if

          !// If E > 1, the sedimentation goes on through more than one layer
          do while (E .gt. 1._r8 .and. L2 .gt. LOWERL)
             ! The new concentration is the old minus the sedimented value
             CONC(L2) = CONC(L2) - HSED

             !// For L2 .GT. LOWERL, the sedimented value is added in the
             !// box below.
             CONC(L2-1) = CONC(L2-1) + HSED

             ! Then update E value (fraction of boxes to sediment through)
             E  = E - 1._r8

             ! And in case E still .GT. 1.d0, sediment to the layer below
             L2 = L2 - 1

             ! Keep going until E .LE. 1.d0 .OR. L2.EQ.LOWERL
          end do

          ! After the above loop we now have THREE possibilities:
          !   1.) E .LE. 1.d0 .and. L2 .GT. LOWERL
          !   2.) L2 .EQ. LOWERL .and. E .GT. 1.d0
          !   3.) L2 .EQ. LOWERL .and. E .LE. 1.d0
          ! - CASE 1 -
          ! E parts of the sedimentation is sedimented to the box below.
          ! - CASE 2 -
          ! Since L2.EQ.LOWERL, the sedimented value is assumed to rain
          ! out. This is done below, but because this is computed as for
          ! case 3, we set E=1.d0 (all of SED is sedimented out).
          ! - CASE 3 -
          ! The particles fall only E parts of the distance to the box
          ! below, so the sedimented value is E*SED. This value is also
          ! assumed to rain out, because L2.EQ.LOWERL.

          ! Need to check case 2:
          if (L2 .eq. LOWERL .and. E .gt. 1._r8) then
             ! This is the CASE 2.
             E = 1._r8
          end if

          ! This now applies for all CASEs
          HSED = HSED * E            ! Amount sedimented from the box
          CONC(L2) = CONC(L2) - HSED ! Out of layer L2
          if (L2 .gt. LOWERL) then
             ! This only applies to CASE 1!
             ! Sedimentation is added to the below level, but for the
             ! lowermost level (L2.EQ.LOWERL), HSED is assumed to vanish
             ! (rain out). In the Oslo CTM2 LOWERL=1, so this ensures the
             ! calculation to stay inside array bounds.
             CONC(L2-1) = CONC(L2-1) + HSED ! Into layer L2-1
          end if

       else if (CONC(L) .lt. 0._r8) then
          ! So, the concentration is negative. This is wrong.
          write(6,'(a,3i4,es16.6)') f90file//':'//subr// &
               ': CONC < 0, setting to 0', L,NR,NB,CONC(L)
          CONC(L) = 0._r8
       end if ! IF (CONC(L).GT.0._r8) THEN

    end do

    !// --------------------------------------------------------------------
  end subroutine sedimentation2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine SED(Y,S,L,Z,LL,W,DT,NR,NB)
    !// --------------------------------------------------------------------
    !// SEDIMENTATION OF ATMOSPHERIC PARTICLE Y IN SOLID STATE WITH SPEED W
    !//
    !// INPUT:
    !//        L       - NUMBER OF ALTITUDE LEVELS
    !//        Z(LL)   - STRAGGLED ALTITUDE VALUES, CM
    !//        LL=L+1  - STRAGGLED GRID
    !//        W       - FALLING SPEED, CM/SEC
    !//        DT      - TIME STEP, SEC
    !//        S       - SEDIMENTED CONCENTRATION OF Y, cm(-3)
    !//
    !// OUTPUT:
    !//        Y(L)    - IMPLICIT CONCENTRATION OF Y, cm(-3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Inout/Ouput
    integer :: L,LL,NR,NB
    real(r8) ::  Y(L),S(L),Z(LL),W(L),DT
    !// Locals
    integer :: I,II
    real(r8) :: E,HSED
    !// --------------------------------------------------------------------
    do I = 1, L
       II=I
       E=W(I)*DT/(Z(I+1)-Z(I)) ! E parts of box is sedimented
       !if (E.gt.1._r8) then
       !print*,'E is greater than 1',E,I,Y(I),S(I)
       !end if
       !CAStest > in case we sediment more than we have
       if (S(I) .gt. Y(I)) then
          print*, '****** TOO MUCH SEDIMENTATION ******',S(I),Y(I),I, &
               NR,NB,Y(I+1),Y(I+2)
          S(I)=Y(I)
       end if
       !CAStest <
3      if (E .le. 1._r8) GOTO 2 ! In kept in box > to 2
       Y(II)=Y(II)-S(I)        ! New conc: old - sedimented
       ! for II=1, Y(II) is the new sfc conc.
       if (II .eq. 1) cycle !// Go to next I
       ! Following is for the L>1
       Y(II-1)=Y(II-1)+S(I)    ! Into level below: old + sedimented
       E=E-1._r8                ! sedimented down to the
       II=II-1                 ! box below
       GOTO 3                  ! until E<1
2      HSED = S(I) * E         ! How much that is sedimented.
       Y(II) = Y(II) - HSED    ! From box II to either
       if (II .eq. 1) cycle    ! sfc, go to next I
       Y(II-1) = Y(II-1) + HSED ! or the layer below
    end do
    !// --------------------------------------------------------------------
  end subroutine SED
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function AM_BIN(T,P,AK)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    real(r8), intent(in) :: T,P,AK(7)
    !// Local variables
    real(r8)             :: X
    !// --------------------------------------------------------------------
    X = X_CARS(T,P,AK)
    AM_BIN = MAX(55.51_r8 * X/(1._r8 - X), 0._r8)
    !// --------------------------------------------------------------------
  end function AM_BIN
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function X_CARS(T,P,AK)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    real(r8), intent(in) :: T,P,AK(7)
    !// Local variables
    real(r8)             :: ASUM1,ASUM2
    !// --------------------------------------------------------------------
    ASUM1 = AK(1) + AK(2)/T
    ASUM2 = 2._r8 * (AK(3) + AK(4)/T)
    X_CARS = (-ASUM1 - sqrt(ASUM1*ASUM1 - 2._r8*ASUM2 &
         * (AK(5) + AK(6)/T + AK(7)*log(T) - log(P)) ) ) / ASUM2
    !// --------------------------------------------------------------------
  end function X_CARS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function HENRIC(T,P,Q)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    real(r8), intent(in) :: T,P,Q(10)
    !// Local variables
    real(r8)             :: TR,TR2,PR,PR2,AHEN
    !// --------------------------------------------------------------------
    PR   = LOG(P) + 18.4_r8
    TR   = 1.e4_r8 * (1._r8 / T - 1._r8 / 230._r8)
    TR2  = TR*TR
    PR2  = PR*PR
    AHEN = Q(1)+Q(2)*TR2 &
         + (Q(3)+Q(4)*TR+Q(5)*TR2+Q(6)*TR2*TR)*PR &
         + (Q(7)+Q(8)*TR+Q(9)*TR2)*PR2 &
         + Q(10)*TR*PR2*PR
    HENRIC = exp(AHEN)
    !// --------------------------------------------------------------------
  end function HENRIC
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function RHO_S_CAR(T,AM)
    !// --------------------------------------------------------------------
    !// DENSITY OF H2SO4 IN A BINARY SOLUTION, KG/M(3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    real(r8), intent(in) :: AM,T
    !// Local variables
    real(r8)             :: T2
    !// --------------------------------------------------------------------
    T2 = T*T
    RHO_S_CAR = 1000._r8 + (123.64_r8 - 5.6e-4_r8 * T2) * AM &
         - (29.54_r8 - 1.81e-4_r8*T2) * AM**1.5_r8 &
         + (2.343_r8 - 1.487e-3_r8*T - 1.324e-5_r8*T2) * AM*AM
    !// --------------------------------------------------------------------
  end function RHO_S_CAR
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function RHO_N_CAR(T,AM)
    !// --------------------------------------------------------------------
    !// DENSITY OF HNO3 IN A BINARY SOLUTION, KG/M(3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    real(r8), intent(in) :: AM,T
    ! Local variables
    real(r8)             :: T2
    !// --------------------------------------------------------------------
    T2 = T*T
    RHO_N_CAR = 1000._r8 + (85.11_r8 - 5.04e-4_r8*T2) * AM &
         - (18.96_r8 - 1.427e-4_r8*T2) * AM**1.5_r8 &
         + (1.458_r8 - 1.198e-3_r8*T - 9.703e-6_r8*T2) * AM*AM
    !// --------------------------------------------------------------------
  end function RHO_N_CAR
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine LOGN(R,N,S,R0,F)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    integer, intent(in) :: N
    real(r8), intent(in)  :: R(N),S,R0
    real(r8), intent(out) :: F(N)
    !// Local variables
    real(r8)            :: PI,SL2,R0R
    integer             :: I
    !// --------------------------------------------------------------------
    PI  = CPI
    SL2 = LOG(S)**2
    R0R = 1._r8 / R0
    do I = 1, N
       F(I) = 1._r8 / (SQRT(2._r8 * PI * SL2) * R(I)) &
            * exp(-(log(R(I) * R0R))**2 / (2._r8*SL2))
    end do
    !// --------------------------------------------------------------------
  end subroutine LOGN
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine MOM(R,N,F,S1,S2,S3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// I/O variables
    integer, intent(in)     :: N
    real(r8), intent(in)    :: R(N),F(N)
    real(r8), intent(inout) :: S1,S2,S3
    ! Local variables
    real(r8)               :: DR, R2
    integer                :: I
    !// --------------------------------------------------------------------
    S1 = 0._r8
    S2 = 0._r8
    S3 = 0._r8
    do I = 1, N
       IF (I.EQ.1) THEN
          DR = R(2) - R(1)
       ELSE IF (I.EQ.N) THEN
          DR = R(N) - R(N-1)
       ELSE
          DR = 0.5_r8*(R(I+1) - R(I-1))
       END IF
       R2 = R(I)*R(I)
       S1 = S1 + F(I)*R(I)*DR
       S2 = S2 + F(I)*R2*DR
       S3 = S3 + F(I)*R2*R(I)*DR
    end do
    !// --------------------------------------------------------------------
  end subroutine MOM
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DIS_LN(N,R,F,R0,S,S1,S2,S3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)     :: N
    real(r8), intent(in)    :: R(N),R0,S
    real(r8), intent(inout) :: S1,S2,S3
    real(r8), intent(out) :: F(N)
    !// --------------------------------------------------------------------
    call LOGN(R,N,S,R0,F)
    call MOM(R,N,F,S1,S2,S3)
    !// --------------------------------------------------------------------
  end subroutine DIS_LN
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  !// In addition, if one wants to use aerosol surface area to
  !// calculate H2SO4 mixing ratio.
  !// ----------------------------------------------------------------------
  subroutine sps_SURF(TAIR,PMEAN,H2O,AREA,MASMIX_H2SO4)
    !// --------------------------------------------------------------------
    !// Subroutine to calculate H2SO4 mass mixing ratio from aerosol
    !// surface area, assuming the aerosols to be described by a
    !// log-normal distribution, and the particles in equilibrium with
    !// the water.
    !//
    !// Input:
    !//     TAIR:    (K)            Ambient air temperature
    !//     PMEAN:   (mb)           Ambient air pressure
    !//     H2O:     (mol/mol)      Volume mixing ratio of water vapor
    !//     AREA:    (cm**2/cm**3)  Particle surface area volume of air
    !//
    !// Output:
    !//     MASMIX_H2SO4            Mass mix.rat. of H2SO4 (kg H2SO4/ kg air)
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR, R_UNIV
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: &
         TAIR, &       ! Ambient air temperature [K]
         PMEAN, &      ! Ambient air pressure [mb]
         H2O, &        ! Volume mixing ratio of water vapor
         AREA          ! Particle surface area per volume of air [cm2/cm3]
    !// Output
    real(r8), intent(out) :: &
         MASMIX_H2SO4  ! Mass mixing ratio of H2SO4 [kg H2SO4/ kg air]

    !// Parameters
    !// Median radius and geometric standard deviation
    !// of assumed log-normal size distribution:
    real(r8), parameter :: &
         R_MEDIAN = 0.0725e-6_r8, & ! Median radius [m]
         GSTDEV   = 1.86_r8         ! standard deviation
    !//Physical constants:
    real(r8), parameter :: &
         MAIR     = M_AIR * 1.e-3_r8, & ! Molar weight of dry air (kg/mole)
         RGAS     = R_UNIV              ! Universal gas constant (J/(mole K))

    !// Local variables
    real(r8) :: &
         PAIR, &                 ! Air pressure
         PPWV, &                 ! H2O partial pressure
         WSAS, DENSITY, &
         RHOAIR, VOLUME, FACTOR
    !// --------------------------------------------------------------------

    !// For particles with area less than 1.d-20 cm2/cm3, set H2SO4=0
    if (AREA .lt. 1.e-20_r8) then
       MASMIX_H2SO4 = 0._r8
       return
    end if

    !// Air and H2O pressure given in [Pa]
    PAIR = PMEAN * 100._r8 ! Change [mb] into [Pa]
    PPWV = H2O * PAIR      ! H2O-pressure [Pa]

    !// Factor for log-normal distribution
    FACTOR = R_MEDIAN * exp(2.5_r8*(log(GSTDEV))**2)/3._r8

    !// Particles in equilibrium with water vapor;
    !// Calculate sulfuric acid weightfraction: 
    call sps_WSASAS(TAIR,PPWV,WSAS)

    !// H2SO4 density (kg/m3):
    DENSITY = sps_ROSAS(TAIR,WSAS)

    !// Air density (kg/m3):
    RHOAIR = (MAIR * PAIR) / (RGAS * TAIR)

    !// Convert surface area density (m**2/m**3 air)
    !// into volume density (m**3/ kg air):
    !// In the Oslo CTM, AREA is given as cm2/cm3. Here we convert to
    !// SI units (m2/m3): AREA = (AREA*1.0d+6)*1.0d-4 = AREA*1.0d+2
    VOLUME = (AREA * 100._r8) * FACTOR / RHOAIR

    !// Then the mass mixing ratio of H2SO4
    MASMIX_H2SO4 = WSAS * DENSITY * VOLUME

    !// --------------------------------------------------------------------
  end subroutine sps_SURF
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine sps_WSASAS(TAIR,PPWV,WSAS)
    !// --------------------------------------------------------------------
    !// This subroutine calculates the sulfuric acid mass fraction in 
    !// sulfuric acid solution with a plane surface in equilibrium with 
    !// the ambient water vapor partial pressure. 
    !//
    !// Source: D.R. Hanson, A.R. Ravishankara & S. Solomon, JGR 99,3615,1994
    !//         Numerical fit to data Steele and Hamill (1981)
    !//
    !// Input/output variables:
    !// REAL TAIR,PPWV,WSAS
    !//
    !// Input:       
    !//     TAIR:    K         Temperature of ambient air [TAIR.GT. 190 K]
    !//     PPWV:    Pa        Partial pressure of ambient water vapor 
    !//
    !// Output:
    !//     WSAS:              mass fraction of sulfuric acid. [0.1;1]
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: &
         TAIR, &           ! Temperature of ambient air [K]
         PPWV              ! Partial pressure of ambient H2O [Pa]
    !// Output
    real(r8), intent(out) :: &
         WSAS              ! mass fraction of sulfuric acid. [0.1;1]

    !// Local variables
    real(r8) :: LNPPWV,T
    real(r8), parameter :: &
         A0 = -14.458_r8, A1 = 0.19988_r8, A2 = 0.62456_r8, &
         B0 = 3565.0_r8, B1 = 44.777_r8, B2=1.3204_r8
    !// --------------------------------------------------------------------
    LNPPWV = log(PPWV / 100._r8)
    T      = max(TAIR, 190._r8)
    WSAS = ((A0 + A2*LNPPWV)*T + B0)/(B1 + B2*LNPPWV - A1*T)
    WSAS = WSAS / 100._r8
    !// --------------------------------------------------------------------
  end subroutine sps_WSASAS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function sps_ROSAS(TAIR,WSA)
    !// --------------------------------------------------------------------
    !// Density of liquid sulfuric acid solution.
    !//
    !// Source: John H.Perry (ed.):Chemical Engineers' Handbook,
    !//         McGraw-Hill, New York 1963, p. 3-79 & 3-
    !//
    !// 0 to 100 % has been fitted with a polynomium of two variables
    !// of order 5 in W and lineary in T. Fit quality better than 0.5 
    !//
    !// Input:  TAIR: Temperature  (K)
    !//         WSA:  Weight fraction of H2SO4  [0;1] 
    !// Output: Density of sulfuric acid solution  (kg/m**3)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: TAIR,WSA
    !// Locals
    integer :: I,J
    real(r8) :: C(6), A(6), B(2,6)
    !// --------------------------------------------------------------------
    
    A(1:6) = (/ 1.00190e+03_r8, 5.50496e+02_r8, 1.54093e+03_r8, &
               -4.89219e+03_r8, 7.56555e+03_r8, -3.92739e+03_r8/)

    B(1,1:6) = (/ 1.98378e+01_r8, 1.02256e+03_r8, -1.48665e+03_r8, &
                  -7.24651e+02_r8, 3.68348e+03_r8, -2.22159e+03_r8 /)
    B(2,1:6) = (/ -6.97011e-02_r8, -3.59886_r8, 5.24992_r8, &
                  2.54047_r8, -1.29355e+01_r8, 7.80553_r8 /)

    do I = 1, 6
       C(I) = A(I) + B(1,I) + B(2,I) * TAIR
    end do
    sps_ROSAS = C(1) + &
         WSA*(C(2) + WSA*(C(3) + WSA*(C(4) + WSA*(C(5) + WSA*C(6)))))
    !// --------------------------------------------------------------------
  end function sps_ROSAS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine PSC_diagnose(NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// By Ole Amund Sovde, 25. February 2005
    !// Checking number of boxes with PSCs and where cold aerosols may
    !// possibly form.
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_met, only: T
    use cmn_oslo, only: LMTROP
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: &
         NDAY,NMET,NOPS

    !// Locals
    integer :: &
         NROPSC1, NROPSC2, NROAER ! Counters for PSC and SA boxes


    real(r8), parameter :: &
         TCRIT1 = 205._r8, &       ! Crit. temp. PSC1
         TCRIT2 = 190._r8, &       ! Crit. temp. PSC2
         TAER = 220._r8            ! Crit. temp. Cold aerosol

    integer :: I,J,L,COUNTER ! Loop counters
    integer :: imax,jmax,lmax
    real(r8) :: maxPSC1
    !// --------------------------------------------------------------------

    !// Initialization:
    NROPSC1 = 0
    NROPSC2 = 0
    NROAER = 0

    !// This does not mean PSCs / aerosols are formed in all these boxes,
    !// only that the temperature is cold enough.

    !// Check for PSCs on/off...
    if (LPSC) then
       ! Checking the stratosphere
       maxPSC1 = 0._r8
       do J = 1, JPAR
          do I = 1, IPAR
             do L = LMTROP(I,J)+1, LPAR-1  ! Layers where PSCs are allowed
                ! Check temp vs critical temperature
                if (PSC2(L,I,J).GT.0._r8) NROPSC2=NROPSC2+1
                if (PSC1(L,I,J).GT.0._r8) NROPSC1=NROPSC1+1
                if (PSC1(L,I,J) .gt. maxPSC1) then
                   maxPSC1 = PSC1(L,I,J)
                   imax = I
                   jmax = J
                   lmax = L
                end if
             end do
          end do
       end do
    end if

    !// Check for AEROSOLSs on/off...
    if (LAEROSOL) then
       do J = 1, JPAR
          do I = 1, IPAR
             do L = LMTROP(I,J)+1, LPAR-1  ! Layers where PSCs are allowed
                ! Check temp vs critical temperature
                if (T(I,J,L) .lt. TAER .and. T(I,J,L) .GT. TCRIT1) then
                   NROAER=NROAER+1
                end if
             end do
          end do
       end do
    end if

    write(6,'(a2,i5,2i3,A14,I6,A2,I6,A7,es12.3,3i4)') &
         'At',NDAY,NMET,NOPS,' PSC1/PSC2 in', &
         NROPSC1,'/',NROPSC2,' boxes.',maxPSC1,imax,jmax,lmax
    write(6,'(A27,I6,A7)') 'Cold aerosols may form in ', &
         NROAER,' boxes.'
    !// --------------------------------------------------------------------
  end subroutine PSC_diagnose
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_psc_constants()
    !// --------------------------------------------------------------------
    !// Defining coefficients for PSC calculations. These used to be set for
    !// each i,j,l, but now we set them once.
    !// --------------------------------------------------------------------
    use cmn_oslo, only: trsp_idx
    use strat_h2o, only: LOLD_H2OTREATMENT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: ji
    real(r8) :: tmin,tmax,drmean, delta
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'set_psc_constants'
    !//---------------------------------------------------------------------

    !// Coefficients for subroutine CARS
    !// Coefficients for the calculation of a binary solution of
    !// HNO3 and H2SO4
    !// --------------------------------------------------------------------
    AK_CARS(:,1) = (/ &
         -3.9136e1_r8, 6.3584e3_r8, 8.329e1_r8, -1.7650e4_r8, 1.9853e2_r8, &
         -1.1948e4_r8,-2.8469e1_r8 /)
    AK_CARS(:,2) = (/ &
         -2.1661e1_r8, 2.7242e3_r8, 5.181e1_r8, -1.5732e4_r8,  4.7004e1_r8, &
         -6.9690e3_r8,-4.6183_r8 /)

    !// Coefficients for the calculation of Henry's law constant of
    !// HNO3 and H2SO4
    !// --------------------------------------------------------------------
    Q_CARS(:,1) = (/ &
         14.5734_r8,   0.0615994_r8,-1.14895_r8,  0.691693_r8, &
         -0.098863_r8, 0.0051579_r8, 0.123472_r8,-0.115574_r8, &
         0.0110113_r8, 0.0097914_r8 /)
    Q_CARS(:,2) = (/ &
         14.47_r8,     0.0638795_r8, &
         -3.29597_r8,  1.778224_r8, -0.223244_r8, 0.0086486_r8, &
          0.536695_r8,-0.335164_r8,  0.0265153_r8,0.0157550_r8/)

    !// Coefficients for sedimentation (Kasten, 1968)
    !// Exponents for wsed_WF
    ! ----------------------------------------------------------------
    wsed_AMN(:) = (/ 1.e-6_r8, 1.e-5_r8, 1.e-5_r8, 1.e-4_r8, &
                     1.e-3_r8, 1.e-2_r8, 1.e-1_r8, 1._r8 /)

    !// Particle radius (size bins) [microns]
    !// --------------------------------------------------------------
    wsed_RF(:)  = (/ 0.003_r8, 0.01_r8, 0.03_r8, 0.1_r8, &
                     0.3_r8, 1._r8, 3._r8, 10._r8 /)

    !// Falling speed for every 2km (wsed_WF*wsed_AMN = cm/sec)
    !// for each size bin
    ! ----------------------------------------------------------------
    wsed_WF(:,1) =(/ &
         4.1150_r8,5.1727_r8,6.5809_r8,8.4837_r8,11.096_r8,14.747_r8, &
         20.005_r8,27.337_r8,37.362_r8,51.065_r8,69.792_r8,95.383_r8, &
         130.00_r8,176.66_r8,239.37_r8,323.44_r8,435.86_r8,586.05_r8, &
         782.70_r8,1038.3_r8,1368.8_r8,1793.3_r8,2335.9_r8,3025.5_r8, &
         3893.3_r8,4988.9_r8,6391.9_r8,8177.8_r8,10492._r8,13510._r8, &
         17463._r8/)
    wsed_WF(:,2) =(/ &
         1.4283_r8,1.7825_r8,2.2537_r8,2.8901_r8,3.7630_r8,4.9831_r8, &
         6.7369_r8,9.1809_r8,12.522_r8,17.090_r8,23.332_r8,31.862_r8, &
         43.399_r8,58.953_r8,79.856_r8,107.88_r8,145.35_r8,195.41_r8, &
         260.96_r8,346.17_r8,456.33_r8,597.82_r8,778.68_r8,1008.5_r8, &
         1297.8_r8,1663.0_r8,2130.7_r8,2726.0_r8,3497.3_r8,4503.5_r8, &
         5821.0_r8/)
    wsed_WF(:,3) =(/ &
         4.8096_r8,5.8811_r8,7.3058_r8,9.2287_r8,11.864_r8,15.544_r8, &
         20.814_r8,28.140_r8,38.159_r8,51.858_r8,70.582_r8,96.165_r8, &
         130.77_r8,177.43_r8,240.13_r8,324.20_r8,436.61_r8,586.79_r8, &
         783.42_r8,1039.0_r8,1369.5_r8,1794.0_r8,2336.5_r8,3026.1_r8, &
         3893.9_r8,4989.6_r8,6392.5_r8,8178.4_r8,10492._r8,13511._r8, &
         17463._r8/)
    wsed_WF(:,4) =(/ &
         2.3182_r8,2.6796_r8,3.1600_r8,3.8083_r8,4.6970_r8,5.9377_r8, &
         7.6918_r8,10.115_r8,13.440_r8,17.995_r8,24.227_r8,32.743_r8, &
         44.268_r8,59.812_r8,80.706_r8,108.72_r8,146.19_r8,196.23_r8, &
         261.76_r8,346.95_r8,457.10_r8,598.57_r8,779.42_r8,1009.3_r8, &
         1298.5_r8,1663.8_r8,2131.4_r8,2726.7_r8,3498.1_r8,4504.3_r8, &
         5821.7_r8/)
    wsed_WF(:,5) =(/ &
         1.4008_r8,1.5224_r8,1.6778_r8,1.8818_r8,2.1565_r8,2.5363_r8, &
         3.0490_r8,3.7427_r8,4.7117_r8,6.0547_r8,7.9057_r8,10.441_r8, &
         13.882_r8,18.533_r8,24.790_r8,33.186_r8,44.418_r8,59.418_r8, &
         79.065_r8,104.61_r8,137.64_r8,180.08_r8,234.32_r8,303.27_r8, &
         390.04_r8,499.61_r8,639.90_r8,818.50_r8,1049.9_r8,1351.8_r8, &
         1747.0_r8/)
    wsed_WF(:,6) =(/ &
         1.3188_r8,1.3890_r8,1.4725_r8,1.5736_r8,1.6985_r8,1.8575_r8, &
         2.0299_r8,2.2259_r8,2.5077_r8,2.9122_r8,3.4878_r8,4.2890_r8, &
         5.3987_r8,6.9179_r8,8.9785_r8,11.757_r8,15.484_r8,20.463_r8, &
         26.993_r8,35.491_r8,46.488_r8,60.619_r8,78.688_r8,101.66_r8, &
         130.58_r8,167.10_r8,213.86_r8,273.40_r8,350.54_r8,451.17_r8, &
         582.93_r8/)
    wsed_WF(:,7) =(/ &
         1.1264_r8,1.1738_r8,1.2279_r8,1.2903_r8,1.3635_r8,1.4510_r8, &
         1.5235_r8,1.5776_r8,1.6521_r8,1.7559_r8,1.9027_r8,2.1034_r8, &
         2.3936_r8,2.8074_r8,3.3872_r8,4.1873_r8,5.2775_r8,6.7402_r8, &
         8.6733_r8,11.201_r8,14.482_r8,18.706_r8,24.114_r8,30.993_r8, &
         39.660_r8,50.615_r8,64.642_r8,82.506_r8,105.65_r8,135.85_r8, &
         175.38_r8/)
    wsed_WF(:,8) =(/ &
         1.2280_r8,1.2745_r8,1.3264_r8,1.3848_r8,1.4509_r8,1.5269_r8, &
         1.5769_r8,1.5942_r8,1.6183_r8,1.6515_r8,1.6972_r8,1.7489_r8, &
         1.8237_r8,1.9311_r8,2.0855_r8,2.3077_r8,2.6247_r8,3.0564_r8, &
         3.6503_r8,4.4497_r8,5.5071_r8,6.8847_r8,8.6620_r8,10.934_r8, &
         13.808_r8,17.453_r8,22.124_r8,28.080_r8,35.799_r8,45.868_r8, &
         59.050_r8 /)

    do JI = 1, wsed_LL
       wsed_ZF(JI) = 2._r8 * real(JI - 1, r8) ! km height
    end do
    

    !// The log-normal distribution
    ! ----------------------------------------------------------------
    ! Particle size bins
    ! ----------------------------------------------------------------
    !// Particle radius (cm)
    !// Originally calculated as: R_A(JI)=1.d-4*10**(0.08d0*(JI-24))
    !// by using N_B=40, this gave the range of R_A:
    !//   I=1:  (0.08*(1-24)  = -1.84
    !//   I=40: (0.08*(40-24) =  1.28
    tmin = -1.84_r8
    tmax = 1.28_r8
    !// We keep min/max of that range, and adjust to a smaller number
    !// of bins, by using a fixed log-interval.
    !// However, we must make this compatible with a different N_B:
    !// The original drmean was 0.08 for N_B=40:
    !//   drmean = (tmax - tmin) / (real(N_B, r8) - 1._r8)
    drmean = (tmax - tmin) / 39._r8
    !// Scale to use for N_B != 40, so that R_A covers the range
    !// [1.d-4*10**(tmin), 1.d-4*10**(tmax)]
    delta =  39._r8 / real(N_B-1, r8)
    !// Adjust to N_B: The old distribution asumed the exponential term
    !// to be zero at JI-0.6N_B, and is not needed now; we only use the
    !// fixed interval (equivalent methods):
    write(6,'(a)') f90file//':'//subr//': Setting particle R_A'
    do JI = 1, N_B
       !// Radius [cm] from tmin to tmax with fixed log-intervals
       R_A(JI) = 1.e-4_r8 * 10._r8**(tmin + drmean*(real(JI-1,r8)*delta))
    end do
    DR(1) = R_A(2) - R_A(1)
    do JI = 2, N_B - 1
       DR(JI) = 0.5_r8 * (R_A(JI+1) - R_A(JI-1))   ! d(PARTICLE RADIUS)
    end do
    DR(N_B) = R_A(N_B) - R_A(N_B-1)

    !// Volume of particles in bin
    do JI = 1, N_B
       VOL_PAR(JI) = CPI * 4._r8/3._r8 * R_A(JI)**3._r8 * DR(JI)
    end do

    write(6,'(a)') '      Bin    R_A             DR              VOL_PAR'
    do JI = 1, N_B
       write(6,'(a,i3,3es16.6)') '  Bin:',JI,R_A(JI),DR(JI),VOL_PAR(JI)
    end do

    !// Momentums of log-normal distribution
    !// --------------------------------------------------------------
    !// Ole Amund Sovde, March 2007: F_SAD is not used in the Oslo CTM2, but
    !// S2_SA and S3_SA are.
    call DIS_LN(N_B,R_A,F_SAD,R0_SA, S_SA, RA,  S2_SA, S3_SA )
    call DIS_LN(N_B,R_A,F_NAT,R0_NAT,S_NAT,RNAT,S2_NAT,S3_NAT)
    call DIS_LN(N_B,R_A,F_ICE,R0_ICE,S_ICE,RICE,S2_ICE,S3_ICE)

    !// Initialize
    PSC1(:,:,:)  = 0._r8
    PSC2(:,:,:)  = 0._r8
    VAER(:,:,:)  = 0._r8
    SAT(:,:,:)   = .false.

    if (.not.( (trsp_idx(  4).gt.0) .and. (trsp_idx(124).gt.0) .and. &
         (trsp_idx(147).gt.0) ) ) then
       write(6,'(a)') '* oc_microphysics.f:'
       write(6,'(a)') '  Heterogeneous chemistry requires transported'// &
            'HNO3, HNO3s and NOy!, and not H2Os'
       stop
    end if
    if (.not.LOLD_H2OTREATMENT &
         .and. (trsp_idx(125).lt.0 .or. trsp_idx(114).lt.0) ) then
       write(6,'(a)') '* oc_microphysics.f:'
       write(6,'(a)') '  Heterogeneous chemistry requires '// &
            'transported H2O and H2Os'
       stop
    end if
    if (LOLD_H2OTREATMENT &
         .and. (trsp_idx(125).gt.0 .or. trsp_idx(114).gt.0) ) then
       write(6,'(a)') '* oc_microphysics.f:'
       write(6,'(a)') '  Heterogeneous chemistry requires '// &
            'non-transported H2O and H2Os'
       stop
    end if

    write(6,'(a)') '* PSC calculations initialized'
    !// --------------------------------------------------------------------
  end subroutine set_psc_constants
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function getTnat(phno3,ph2o)
    !// --------------------------------------------------------------------
    !// Calculate TNAT from Hanson & Mauersberger 1988:
    !// Pressures are in Torr.
    !//
    !// log(phno3) = m(T)log(ph2o) + b(T)
    !//   m(T) = -2.7836d0 - 0.00088d0*T               = a1 + a2*T
    !//   b(T) = 38.9855d0 - 11397.d0/T + 0.009179d0*T = a3 + a4/T + a5*T
    !// where
    !//   a1 = -2.7836d0 ; a2 = -0.00088d0 ; a3 = 38.9855d0
    !//   a4 = -11397.d0 ; a5 = 0.009179d0
    !//
    !// Re-arranging
    !//   log(phno3) = (a1 + a2*T)log(ph2o) + a3 + a4/T + a5*T
    !//
    !//   (a5 + a2*log(ph2o))*T**2 + (a1*log(ph2o)
    !//                             + a3 - log(phno3))*T + a4 = 0
    !//
    !//   aa = a5 + a2*log(ph2o)
    !//   bb = a1*log(ph2o) + a3 - log(phno3)
    !//   cc = a4
    !//
    !// Solution
    !//   Tnat = -bb/2aa +- sqrt(bb**2 - 4aacc)/2aa
    !//
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8),intent(in)   :: phno3,ph2o ! torr

    !// locals
    real(r8),parameter :: &
         a1 = -2.7836_r8, &
         a2 = -0.00088_r8, &
         a3 = 38.9855_r8, &
         cc = -11397._r8, &  ! cc = a4
         a5 = 0.009179_r8
    real(r8)  :: tnat,aa,bb,d2,t1,t2
    !// --------------------------------------------------------------------

    aa = a5 + a2 * log10(ph2o)
    bb = a1 * log10(ph2o) + a3 - log10(phno3)
               
    d2 = bb * bb - 4._r8 * aa * cc
    if (d2 .ge. 0._r8) then
       !// Need to calculate both solutions for the 2. order polynom
       !// Remember that aa may be negative
       t1 = ( - bb + sqrt(d2) ) * 0.5_r8 / aa
       t2 = ( - bb - sqrt(d2) ) * 0.5_r8 / aa
       tnat = max(t1,t2)
       if (tnat .lt. 0._r8) tnat = 0._r8 ! no PSCs.
    else
       !// complex number, treat as no PSCs
       tnat = 0._r8              ! no PSCs.
    end if

    !// return TNAT
    getTnat = tnat
    !// --------------------------------------------------------------------
  end function getTnat
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module psc_microphysics
!//=========================================================================

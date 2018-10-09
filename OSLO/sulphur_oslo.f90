!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Sulphur stuff for Oslo chemistry.
!//=========================================================================
module sulphur_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: sulphur_oslo
  !// DESCRIPTION: Sulphur stuff for Oslo chemistry.
  !// ----------------------------------------------------------------------
  !// Contains:
  !//   subroutine TCRATE_CONST_S
  !//   subroutine TCRATE_TP_S_IJ
  !//
  !// Amund Sovde, September 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR,JPAR
  use cmn_parameters, only: MINTEMP, TEMPRANGE, AVOGNR, R_ATM
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Parameters for sulphur chemistry
  real(r8), parameter :: HPLUS = 3.16e-5_r8 !// Assume pH=4,5

  real(r8), dimension(TEMPRANGE) :: &
!       R4071a, &
!       R4074, &
       R4077, &
       R3877, &
       R4076

  !// Fields for DMS short term variations
  real(r8), dimension(IPAR,JPAR,12) :: DMSseaconc
  !// ----------------------------------------------------------------------
  save
  public
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine TCRATE_CONST_S(LSULPHUR)
    !// --------------------------------------------------------------------
    !// Find constant and temperature dependent reaction rates for
    !// the sulphur chemistry.
    !//
    !//
    !// Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input parameters
    logical, intent(in) :: LSULPHUR

    !// Locals
    real(r8) :: ZTEM,TEM, FC
    integer :: I
    !// --------------------------------------------------------------------

    if (.not. LSULPHUR) then
!       R4071a(:) = 0._r8
!       R4074(:)  = 0._r8
       R4077(:)  = 0._r8
       R3877(:)  = 0._r8
       R4076(:)  = 0._r8
    else
       do I = 1, TEMPRANGE

          TEM  = real(I + MINTEMP, r8)
          ZTEM = 1._r8 / TEM

!// OH + DMS --> H2O + CH3SCH2 (abstraction)                JPL (1997), I15
!          R4071a(I) = 1.2e-11_r8 * exp(-260._r8 * ZTEM)

!// H2S + OH --> --> SO2                                    JPL (1997)
!          R4074(I) = 6.0e-12_r8 * exp(-75._r8 * ZTEM)

!// COS + OH --> --> SO2                                    JPL (1997)
          R4077(I) = 1.1e-13_r8 * exp(-1200._r8 * ZTEM)

!// O3P + COS --> CO + SO --> SO2                           JPL (1997)
          R3877(I) = 2.1e-11_r8 * exp(-2200._r8 * ZTEM)

!// OH + CS2 --> --> SO2 + COS                see note I13, JPL (1997)
          R4076(I) = 0._r8

       end do
    end if


    !// --------------------------------------------------------------------
  end subroutine TCRATE_CONST_S
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TCRATE_TP_S_IJ(LM,LMTROP,TEMP, AIR_MOLEC, &
       H2O_MOLEC, PMIDL, LSULPHUR, LAT, AREA, &
       LWEPAR, CLDFR, CLDLWC, CLDIWC, RAIN,DV, AIR, &
       R4071b, RTOT4072,RAQ0172, RAQ1572, RAQ1772, &
       RCATSO2 )
    !// --------------------------------------------------------------------
    !// Find sulphur reaction rates
    !// 1. dependent on T and p.
    !// 2. aquous rates in cloud droplets
    !// 3. catalytic oxidation of SO2 by metals inside clouds.
    !// 4. sub cloud scavenging rates of SO2 and SO4 (to be treated as
    !//    loss in chemistry)
    !//    HAS BEEN REMOVED June 2012. Subcloud scavenging is covered
    !//    by the new scavenging scheme.
    !//
    !// Amund Sovde, June 2012, August 2009
    !// --------------------------------------------------------------------
    use utilities_oslo, only: rate3B
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LM, LMTROP, LWEPAR
    real(r8), intent(in) :: &
         AIR_MOLEC(LM), &
         H2O_MOLEC(LM), &
         PMIDL(LM), &
         CLDFR(LWEPAR), &
         CLDLWC(LWEPAR), &
         CLDIWC(LWEPAR), &
         RAIN(LWEPAR), &
         DV(LM), &
         AIR(LM), &
         TEMP(LM)
    real(r8), intent(in) :: LAT, AREA
    logical, intent(in) :: LSULPHUR
    !// Output
    real(r8), intent(out) ::  &
         R4071b(LM),  RTOT4072(LM), &
         RAQ0172(LM), RAQ1572(LM), RAQ1772(LM),RCATSO2(LM)

    !// Local variables
    real(r8) :: &
         KZERO,KINF, &
         ZTEM,TEM,TZ300,TEMPFAC, &
         O2,DZ, &
         HTMP1, HTMP2, HH, RDUM, RAQDUM,CWVOL,CWVOL_FRAC, &
         EFF_SO2, FAQSO2, FAQO3, FAQH2O2, FAQHO2NO2, K_H2O2_HSO3, &
         AQFIX_HO2NO2_SO2, AQFIX_O3_SO3
    integer :: L,TI,L4km

    !// Parameters
    real(r8), parameter :: INV298 = 1._r8 / 298._r8
    !// Modified universal gas constant, scaled from [J/(mole*K)]
    !// to [atm*L/(mol*K)]
    real(r8), parameter :: RGAS_MOD = R_ATM * 1000._r8
    !// Subcloud scavenging of SO2 Rate= 2.61d-5*RAIN for RAIN kg/m2/hr
    real(r8), parameter :: RMARTIN = 2.61e-5_r8
    !// Subcloud scavenging of SO4
    real(r8), parameter :: RBERGE = 0.52_r8    !// C*Em, where C=5.2 and Em=0.1
    !// --------------------------------------------------------------------

    !// Initialize (always, they have to be zero when not using sulphur
    !// chemistry).
    R4071b(:)   = 0._r8
    RTOT4072(:) = 0._r8
    RAQ0172(:)  = 0._r8
    RAQ1572(:)  = 0._r8
    RAQ1772(:)  = 0._r8
    RCATSO2(:)  = 0._r8

    !// If no sulfur, return
    if (.not. LSULPHUR) return


    !// Gas phase first, then aquous
    do L = 1, LMTROP

       TEM   = TEMP(L)
       ZTEM  = 1._r8 / TEM
       TZ300 = TEM / 300._r8
       O2    = 0.21_r8 * AIR_MOLEC(L)

!// OH + DMS --> CH3S(OH)CH3 --> 0.75 SO2 + 0.25 MSA (addition) Chin et al., 
       R4071b(L) = (O2 * 1.7e-42_r8 * exp(7810._r8 * ZTEM)) &
            / (1._r8 + (O2 * 5.5e-31_r8 * exp(7460._r8 * ZTEM)))

!// SO2 + OH --> --> HOSO4                           Chin et al. (1996)
       !KZERO = 3.0D-31*((300._R8*ZTEM)**(3.3_R8))*AIR_MOLEC(L)
       !KINF = 1.5D-12
       !RTOT4072(L) = TROE(KZERO,KINF,0.6_R8)
       RTOT4072(L) = rate3B(301, TZ300, AIR_MOLEC(L), 3.e-31_r8, 3.3_r8, &
            1.5e-12_r8, 0._r8, 0.6_r8, 0)
       !// JPL 2010
       !RTOT4072(L) = rate3B(301, TZ300, AIR_MOLEC(L), 3.3e-31_r8, 4.3_r8, &
       !     1.6e-12_r8, 0._r8, 0.6_r8, 0)

       !// Will be modified below

    end do



    !// Modify/set rates due to cloud water/ice
    !// --------------------------------------------------------------------
    !// HO2NO2(aq) + HSO3(aq) <--> 2H(aq) + SO4(aq) + NO3(aq)
    AQFIX_HO2NO2_SO2 = 3.1e5_r8  ! roed bok, tabell 3.10

    !// O3(aq) + SO3(aq) <--> SO4(aq) + O2(aq)
    AQFIX_O3_SO3 = 1.8e4_r8 * (HPLUS**(-0.4_r8)) ! L/mol, MOLLER, 1980

    do L = 1, min(LMTROP, LWEPAR)

       if (CLDFR(L) .gt. 0.05_r8) then ! required minimum for aq. chem
          !// Find volume fraction of cloud water
          !// Volume: Total water mass divided by density (1000kg/m3)
          CWVOL = (CLDLWC(L) + CLDIWC(L)) * AIR(L) * 1.e-3_r8
          !// Fraction of cloud volume
          CWVOL_FRAC = CWVOL / (CLDFR(L) * DV(L))

          if (CWVOL_FRAC .gt. 1.e-8_r8) then
             RAQDUM = 1000._r8 * CLDFR(L) / (AVOGNR * CWVOL_FRAC)
             TEMPFAC = 1._r8 / TEMP(L) - INV298

             !// Find fraction of dissolved tracers

             !// SO2(g) <--> SO2(aq)
             !// A regular Henry expression [mol/(atm*L)], which is to be
             !// modified
             HTMP1 = 1.23_r8 * exp(3020._r8 * TEMPFAC) * RGAS_MOD * TEMP(L)
             !// SO2(aq) <--> HSO3(aq) + H(aq)
             HTMP2 = 1.23e-2_r8 * exp(2010._r8 * TEMPFAC)

             !// Calculate efficient Henry's coefficient
             EFF_SO2 = HTMP1 * (1._r8 + (HTMP2 / HPLUS))

             !// Dissolved SO2: Fraction of dissolved SO2
             RDUM = EFF_SO2 * CWVOL_FRAC
             FAQSO2 = RDUM / (RDUM + 1._r8)


             !// Dissolved O3: O3(g) <--> O3(aq)
             RDUM = CWVOL_FRAC &
                  * 1.13e-2_r8 * exp(2300._r8 * TEMPFAC) * RGAS_MOD * TEMP(L)
             FAQO3 = RDUM / (RDUM + 1._r8)

             !// Dissolved H2O2: H2O2(g) <--> H2O2(aq)
             RDUM = CWVOL_FRAC &
                  * 7.1e4_r8 * exp(6800._r8 * TEMPFAC) * RGAS_MOD * TEMP(L)
             FAQH2O2 = RDUM / (RDUM + 1._r8)

             !// Dissolved HO2NO2: HO2NO2(g) <--> HO2NO2(aq)
             RDUM = CWVOL_FRAC &
                  * 1.2e4_r8 * exp(6900._r8 * TEMPFAC) * RGAS_MOD * TEMP(L)
             FAQHO2NO2 = RDUM / (RDUM + 1._r8)

             !// Modify rates
!// SO2 + OH --> --> H2SO4                           Chin et al. (1996)
             RTOT4072(L) = RTOT4072(L) &
                  * (((1._r8 - FAQSO2) * CLDFR(L)) + 1._r8 - CLDFR(L))


             !// Set rates for aquous reactions
!// O3(aq) + SO3(aq) <--> SO4(aq) + O2(aq)
             RAQ0172(L) = RAQDUM * AQFIX_O3_SO3 * FAQSO2 * FAQO3

!// H2O2(aq) + HSO3(aq) <--> H(aq) + SO4(aq) + H2O
             !// RAQH2O2_SO2(I,J,L)=8.3E5 ! L/mol, MARTIN & DAMSCHEN, 1981
             K_H2O2_HSO3 = (8.0e4_r8 * exp(-3650._r8 * TEMPFAC)) &
                                       / (0.1_r8 + HPLUS)
             RAQ1572(L) = RAQDUM * K_H2O2_HSO3 * HTMP1 &
                  * FAQSO2 * FAQH2O2 / EFF_SO2

!// HO2NO2(aq) + HSO3(aq) <--> 2H(aq) + SO4(aq) + NO3(aq)
             RAQ1772(L) = RAQDUM * AQFIX_HO2NO2_SO2 * HTMP1 &
                  * FAQSO2 * FAQHO2NO2 / EFF_SO2

          end if
       end if !// if (CLDFR(L) .gt. 0.05) then 

       !// Catalytic oxidation of SO2 by metals inside clouds
       !// -----------------------------------------------------------------
       if (LAT .le. 45._r8) then
          RCATSO2(L) = 2.77e-6_r8 * CLDFR(L) !// lifetime 100hrs south of 45N
       else
          RCATSO2(L) = 5.55e-6_r8 * CLDFR(L) !// lifetime 50hrs
       end if
    end do

    !// --------------------------------------------------------------------
  end subroutine TCRATE_TP_S_IJ
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module sulphur_oslo
!//=========================================================================

!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// Constant parameters for CTM run.
!//=========================================================================
module CMN_PARAMETERS
  !//-----------------------------------------------------------------------
  !// MODULE: CMN_PARAMETERS
  !// DESCRIPTION: Constant parameters for CTM run.
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

  !// Physical constants
  !//-----------------------------------------------------------------------
  !// Earth radius [m]
  real(r8), parameter :: A0        = 6371000._r8
  !// Gravitational constant [m/s^2]
  real(r8), parameter :: G0        = 9.80665_r8
  !// pi & 2pi
  real(r8), parameter :: CPI       = 3.141592653589793_r8
  real(r8), parameter :: C2PI      = 2._r8*CPI
  !// Conversions to / from radians
  real(r8), parameter :: CPI180    = CPI/180._r8
  real(r8), parameter :: ZPI180    = 1._r8/CPI180
  !// Von Karman constant
  real(r8), parameter :: VONKARMAN = 0.4_r8
  !// Similarity profiles (from Kansas experiment)
  real(r8), parameter :: CSIM_BETA      = 4.7_r8
  real(r8), parameter :: CSIM_GAMMA     = 15._r8


  !// Chemistry constants
  !//-----------------------------------------------------------------------
  !// Conversion atm to Pa [Pa/atm]
  real(r8), parameter :: atm2Pa = 101325._r8
  real(r8), parameter :: Pa2atm = 1._r8 / atm2Pa
  !// Numbers of kcal in a J
  real(r8), parameter :: J2kcal = 4186.8_r8

  !// Molecular mass for air [g/mol]
  real(r8), parameter :: M_AIR  = 28.97_r8

  !// Universal gas constant [J/(K*mol)] = [m3*Pa/(K*mol)]
  real(r8), parameter :: R_UNIV = 8.31446_r8
  !// Gas constant in "chemical units" [m3*atm/(K*mol)] (~8.205d-5)
  real(r8), parameter :: R_ATM  = R_UNIV * Pa2atm
  !// Specific gas constant for air [J/(K*kg)] (~287.d0)
  real(r8), parameter :: R_AIR  = R_UNIV / M_AIR * 1.e3_r8
  !// Specific gas constant for water vapor [J/(K*kg)] (~461.d0)
  real(r8), parameter :: R_H2O  = R_UNIV / 18.01528_r8 * 1.e3_r8

  !// Avogadro's number [molecules/mol]
  real(r8), parameter :: AVOGNR = 6.022149e23_r8
  !// Boltzmann's const [J/(K*molecules)]
  real(r8), parameter :: BOLTZMANN = 1.38063e-23_r8
  !// Boltzmann's const [mb*cm3/(K*molecules)]
  real(r8), parameter :: KBOLTZ = 1.e4_r8 * BOLTZMANN

  !// Specific heat of dry air at constant pressure [J/(K*kg)
  real(r8), parameter :: cp_air = 1004._r8
  !// Latent heat of vaporization at 0C [J/kg] (Stull, 1988)
  real(r8), parameter :: Lv_0C = 2.501e6_r8
  !// Gradient of Lv between 0C and 100C [(J/kg)/K] (Stull, 1988)
  real(r8), parameter :: dLv_dT = -0.00237e6_r8
  !// Temperature [K] at 0C
  real(r8), parameter :: TK_0C = 273.15_r8
  !// Saturation vapor pressure of H2O [Pa]
  real(r8), parameter :: es_0C = 611.2_r8

  !// Temporal constants
  real(r8), parameter :: secDay  = 86400._r8
  real(r8), parameter :: secYear = 31536000._r8

  !// Limits on temperature for temperature dependent rates
  integer, parameter :: MINTEMP    = 150
  integer, parameter :: MAXTEMP    = 350
  integer, parameter :: TEMPRANGE  = MAXTEMP - MINTEMP


  !// Logical operator for including more extensive debugging
  logical, parameter :: LDEBUG=.true.

  !//-----------------------------------------------------------------------
end module CMN_PARAMETERS
!//=========================================================================

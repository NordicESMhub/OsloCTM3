!//=========================================================================
!// Oslo CTM3 v2015.01
!//=========================================================================
!// Ole Amund Sovde, April 2015
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
  !//   - subroutine oslochem_psc
  !//   subroutine get_psc12_sad
  !//   - subroutine get_mw:
  !//   - SUBROUTINE PSC_1d: Main column model.
  !//   - SUBROUTINE CARS
  !//   - FUNCTION H2O_SAT
  !//   - FUNCTION WSED
  !//   - subroutine sedimentation
  !//   - SUBROUTINE SED
  !//   - FUNCTION AM_BIN
  !//   - FUNCTION X_CARS
  !//   - FUNCTION HENRIC
  !//   - FUNCTION RHO_S_CAR
  !//   - FUNCTION RHO_N_CAR
  !//   - SUBROUTINE LOGN
  !//   - SUBROUTINE MOM
  !//   - SUBROUTINE DIS_LN
  !//   - SUBROUTINE sps_SURF
  !//   - SUBROUTINE sps_WSASAS
  !//   - FUNCTION sps_ROSAS
  !//   - SUBROUTINE PSC_diagnose
  !//   - subroutine set_psc_constants
  !//   - function getTnat
  !//
  !// Amund Sovde, October 2008
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Flags for turning on/off SAD/PSC1/PSC2
  logical, parameter :: &
       LAEROSOL = .false., &
       LPSC = .false.
  !// ----------------------------------------------------------------------
  public
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_psc(BTT,BTEM,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Calculates formation and evolution of PSCs/STS.
    !//
    !// Will calculate PSC1 and PSC2 areas.
    !// 
    !// Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In
    integer, intent(in) :: MP
    real(r8), intent(in)  :: DTCHM
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM
    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in):: BTT
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine oslochem_psc
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
    P1(:) = 0._r8
    !// Retrieve column values of PSC2
    P2(:) = 0._r8
    !// Retrieve column values of sps_PARTAREA (i.e. background SAD)
    SAD(:) = 0._r8
    !// --------------------------------------------------------------------
  end subroutine get_psc12_sad
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_psc_constants()
    !// --------------------------------------------------------------------
    !// DUMMY.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine set_psc_constants
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module psc_microphysics
!//=========================================================================

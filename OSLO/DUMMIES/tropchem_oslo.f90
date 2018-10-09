!//=========================================================================
!// Oslo CTM3 v2015.01
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// DUMMY.
!//=========================================================================
module tropchem_oslo
  !// ----------------------------------------------------------------------
  !// DUMMY routines for tropospheric chemistry.
  !//
  !// Contains:
  !//   - subroutine oslochem_trop
  !//
  !// Amund Sovde, September 2009
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_trop(BTT,BJV,BAIR_MOLEC,BVOL,BTEM,AIRB,BEMIS,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// DUMMY
    !//
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, IDBLK, JDBLK, NPAR
    use cmn_fjx, only: JVN_
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DTCHM
    integer, intent(in) :: MP
    real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK),intent(in) :: BJV
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BAIR_MOLEC,BVOL
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM,AIRB
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BEMIS
    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine oslochem_trop
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module tropchem_oslo
!//=========================================================================

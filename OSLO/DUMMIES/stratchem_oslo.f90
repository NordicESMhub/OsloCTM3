!//=========================================================================
!// Oslo CTM3 v2015.01
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// DUMMY.
!//=========================================================================
module stratchem_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: stratchem_oslo
  !// DECRIPTION: DUMMY.
  !//
  !// Contains:
  !//   subroutine oslochem_strat
  !//   subroutine read_oslo2d
  !//   subroutine update_strat_boundaries
  !//   subroutine set_fam_in_trop
  !//
  !// Amund Sovde, September 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_strat(BTT,BJV,BAIR_MOLEC,BVOL,BTEM,BEMIS,BTTBCK,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// DUMMY
    !//
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, IDBLK, JDBLK, MPBLK, NPAR
    use cmn_fjx, only: JVN_
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DTCHM
    integer, intent(in) :: MP
    real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK,MPBLK),intent(in) :: BJV
    real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK), intent(in) :: BAIR_MOLEC, BVOL
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BEMIS

    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTTBCK
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine oslochem_strat
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_oslo2d2(LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// DUMMY
    !//
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    logical, intent(in) :: LNEW_MONTH
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine read_oslo2d2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_strat_boundaries(BTT,BTTBCK,AIRB,BTEM,DTADV,MP)
    !// --------------------------------------------------------------------
    !// Update stratosphere before transport, when stratospheric
    !// chemistry is not included.
    !// 
    !//
    !// Amund Sovde, June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, IDBLK, JDBLK, NPAR, LPARW, LOSLOCTROP
    use strat_o3noy_clim, only: update_stratNOY, update_stratNOX
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTTBCK
    real(r8), dimension(LPAR,IDBLK,JDBLK),intent(in)  :: AIRB, BTEM
    real(r8), intent(in)  :: DTADV
    integer, intent(in) :: MP
    !// --------------------------------------------------------------------

    if (LOSLOCTROP) then
      !// Update NOy in stratosphere when running without stratospheric
      !// chemistry.
      if (LPARW.le.40) then
        call update_stratNOX(BTT,BTTBCK,AIRB,DTADV,MP)
      else
        call update_stratNOY(BTT,BTTBCK,AIRB,MP)
      end if
    end if

    return

    !// --------------------------------------------------------------------
  end subroutine update_strat_boundaries
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_fam_in_trop(BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// DUMMY.
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR,IDBLK,JDBLK, NPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTTBCK
    integer, intent(in) :: MP
    !// --------------------------------------------------------------------
    return
    !// --------------------------------------------------------------------
  end subroutine set_fam_in_trop
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module stratchem_oslo
!//=========================================================================

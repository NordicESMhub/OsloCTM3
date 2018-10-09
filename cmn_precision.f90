!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Precision parameters for Oslo CTM3.
!//=========================================================================
module cmn_precision
  !//-------------------------------------------------------------------------
  !// MODULE: cmn_precision
  !//
  !// DESCRIPTION: Precision parameters for Oslo CTM3.
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------
  !// Standard precisions
  integer, parameter :: r8 = selected_real_kind(15) !// 8 byte real
  integer, parameter :: r4 = selected_real_kind(6)  !// 4 byte real

  integer, parameter :: i8 = selected_int_kind(15) !// 8 byte integer

  !// Precision used for moments (e.g. SUT, SVT)
  integer, parameter :: rMom = selected_real_kind(6) !// 4 byte real

  !// Precision used for averages (e.g. STTAVG)
  integer, parameter :: rAvg = selected_real_kind(6) !// 4 byte real

  !// Precision used for 5D tendencies array (STTTND)
  integer, parameter :: rTnd = selected_real_kind(6) !// 4 byte real
  !//-----------------------------------------------------------------------
end module cmn_precision
!//=========================================================================

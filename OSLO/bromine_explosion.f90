!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Marit Kollstuen, November 2019
!//=========================================================================
!// Routines for the bromine explosion implementation.
!//=========================================================================
module bromine_explosion
  !// ----------------------------------------------------------------------
  !// MODULE: bcoc_oslo
  !// DESCRIPTION: Routines for the Arctic/Antarctic bromine explosion events.

  use cmn_precision, only: r8
  
  implicit none

contains

  !//-----------------------------------------------------------------------
  !subroutine be_snow_seaice_ij


  !// ----------------------------------------------------------------------
  subroutine be_getspringsummer(YDGRD, JDAY, &
       springday,springend,summermid,summerend,KDAY)
    !// --------------------------------------------------------------------
    !// Taken from bcoc_oslo.f90 by Amund Søvde, 14.11.19
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in)   :: ydgrd
    integer, intent(in)  :: JDAY
    integer, intent(out) :: KDAY,springday,springend,summermid,summerend
    !// --------------------------------------------------------------------

    !// Melting start days NH/SH
    !// NH: 136 16. May, SH: 320 16. November
    if (YDGRD .ge. 0._r8) then
       !// NH
       springday = 136
       springend = 172
       summermid = 219
       summerend = 244
       KDAY = JDAY
    else
       !// SH
       springday = 320
       springend = 354
       !// Adjust summermid and KDAY so days are monotonic
       summermid =  36 + 365
       summerend =  60 + 365
       if (JDAY .le. 60) then
          KDAY = JDAY + 365
       else
          KDAY = JDAY
       end if
    end if

    !// --------------------------------------------------------------------
  end subroutine be_getspringsummer
  !// ----------------------------------------------------------------------


end module bromine_explosion
!//=========================================================================

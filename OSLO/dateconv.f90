!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Date conversion.
!//=========================================================================
module dateconv
  !// ----------------------------------------------------------------------
  !// MODULE: dateconv
  !// DESCRIPTION: Contains routines for date/time calculation.
  !//
  !// Contains
  !//   subroutine itau2idate
  !//   subroutine idate2itau
  !//   function julianday
  !//   subroutine calendarday
  !//
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  ! model times are assumed to be supplied in seconds relative to
  ! a reference year (in this example 1995 (1 Jan, 00:00 UTC))
  integer, parameter :: iyear0 = 1995  ! reference year
  !// ----------------------------------------------------------------------
 
contains

  !// ----------------------------------------------------------------------
  subroutine itau2idate(itaux,idatex)
!-----------------------------------------------------------------------
!**** itau2idate
!
!     purpose
!     -------
!     calculate date from given time in seconds relative to reference year
!
!     parameters
!     ----------
!     on input : itaux contains date/time in seconds
!     on output: idatex contains date in year,month,day,hour,min,sec
!
!     dependencies
!     ------------
!     iyear0 is the reference year for the calculation
!
!     externals
!     ---------
!     subroutines: calendardate
!     funtions:    julianday
!-----------------------------------------------------------------------
 
    implicit none
    integer itaux,idatex(6)
    integer julian,idayy
!
! compute time (hour,min,sec) and number of days
!
    idatex(6)=mod(itaux,60)
    idatex(5)=mod(itaux/60,60)
    idatex(4)=mod(itaux/3600,24)
    idayy=itaux/86400
!
! real calendar
!
    julian=julianday(1,1,iyear0)+idayy
    call calendardate(julian,idatex(2),idatex(3),idatex(1))

  end subroutine itau2idate


  subroutine idate2itau(idatex,itaux)
!-----------------------------------------------------------------------
!**** idate2itau
!
!     purpose
!     -------
!     calculate time in seconds since 1 Jan 1995
!
!     parameters
!     ----------
!     on input : idatex contains date in year,month,day,hour,min,sec
!     on output: itaux contains date/time in seconds rel. to reference year
!
!     dependencies
!     ------------
!     iyear0 is the reference year for the calculation
!
!     externals
!     ---------
!     funtions: julianday
!-----------------------------------------------------------------------

    implicit none
    integer idatex(6),itaux
    integer idaysec
!
! compute the seconds the day is old
! 
    idaysec=idatex(6)+idatex(5)*60+idatex(4)*3600

    itaux=86400*(julianday(idatex(2),idatex(3),idatex(1))-julianday(1,1,iyear0)) &
         +idaysec

  end subroutine idate2itau


  function julianday(mm,id,iyyy)
!-----------------------------------------------------------------------
!**** julianday
!
!     purpose
!     -------
!     calculate julian day from given date
!
!     parameters
!     ----------
!     on input : mm, id, iyyy contain month, day and year
!     on output: julidanday contains the julian day
!
!     reference
!     ---------
!     J. Meeuws, "Astronomical formulea for calculators" 19xx
!-----------------------------------------------------------------------
    implicit none

    integer mm, id, iyyy, julianday
    integer igreg, jy, jm, ja
    parameter (igreg=15+31*(10+12*1582))
!
! handle dates before 0 AD
!
    if (iyyy.eq.0) stop ' ERROR invalid year 0 AD'
    if (iyyy.lt.0) iyyy=iyyy+1
!
! calculate julian day from date in gregorian calendar
!
    if (mm.gt.2) then
       jy=iyyy
       jm=mm+1
    else
       jy=iyyy-1
       jm=mm+13
    endif
    julianday=int(365.25*jy)+int(30.6001*jm)+id+1720995
!
! handle julian calendar
!
    if (id+31*(mm+12*iyyy).ge.igreg) then
       ja=int(0.01*jy)
       julianday=julianday+2-ja+int(0.25*ja)
    endif
    return
  end function julianday


  subroutine calendardate(julian,mm,id,iyyy)
!-----------------------------------------------------------------------
!*** calendardate
!
!     purpose
!     -------
!     calculate date from given julian day
!
!     parameters
!     ----------
!     on input : julian contains the julian day
!     on output: mm, id, iyyy contain month, day and year
!
!     reference
!     ---------
!     J. Meeuws, "Astronomical formulea for calculators" 19xx
!-----------------------------------------------------------------------
    implicit none
    
    integer julian, mm, id, iyyy
    integer igreg, jalpha, ja, jb, jc, jd, je
    parameter (igreg=2299161)
!
! handle gregorian and julian date
!
    if (julian.ge.igreg)then
       jalpha=int(((julian-1867216)-0.25)/36524.25)
       ja=julian+1+jalpha-int(0.25*jalpha)
    else
       ja=julian
    endif
    jb=ja+1524
    jc=int(6680.+((jb-2439870)-122.1)/365.25)
    jd=365*jc+int(0.25*jc)
    je=int((jb-jd)/30.6001)
    id=jb-jd-int(30.6001*je)
    mm=je-1
    if (mm.gt.12) mm=mm-12
    iyyy=jc-4715
    if (mm.gt.2) iyyy=iyyy-1
!
! handle dates before 0 AD
!
    if (iyyy.le.0) iyyy=iyyy-1

    return
  end subroutine calendardate

  !// ----------------------------------------------------------------------
end module dateconv
!//=========================================================================

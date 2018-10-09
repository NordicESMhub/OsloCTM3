!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Cloud routines.
!//=========================================================================
module cloudjx
  !//-----------------------------------------------------------------------
  !// MODULE: cloudjx
  !// DESCRIPTION: Cloud parameters, variables and routines.
  !//
  !// 
  !//
  !// Contains
  !//   !subroutine cloudw
  !//   !subroutine ranset
  !//
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, r4
  use cmn_size, only: IPAR, JPAR, LPAR
  use cmn_fjx, only: NQD_, CBIN_, ICA_
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  !// Parameters
  integer, parameter ::  NRAN_ = 10007  ! dimension for random number
  !// Fetch from fast-JX
  !integer, parameter ::  CBIN_ = 20     ! # of quantized cld fration bins
  !integer, parameter ::  ICA_  = 20000  ! # of Indep Colm Atmospheres

  !// Variables
  !//-----------------------------------------------------------------------
  real(r4)  :: RAN4(NRAN_)       ! Random number set
  integer :: IRAN0

  integer :: CLDFLAG, NRANDO  ! # of ICA's used for LCDFLAG(5)
  !      1 : Clear sky J's
  !      2 : Averaged cloud cover
  !      3 : cloud-fract**3/2, then average cloud cover
  !      4 : Average direct solar beam over all ICAs, invert to get clouds
  !      5 : Random select NRANDO ICA's from all(Independent Column Atmos.)
  !      6 : Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
  !      7 : Use all (up to 4) QCAs (average clouds within each Q-bin)
  !      8 : Calcluate J's for ALL ICAs (up to 20,000 per cell!)

  !// Cloud2 variables:
  !// Fractional Cloud data/storage
  logical :: LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
  real(r8)  :: CLDSTORE(LPAR+1,NQD_,IPAR,JPAR)
  real(r8)  :: SWSTORE(NQD_,IPAR,JPAR)
  integer :: TYPSTORE(LPAR+1,IPAR,JPAR)
  integer :: RANSEED
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'cloudjx.f90'
  !//-----------------------------------------------------------------------
  private
  public NQD_, NRAN_, CBIN_, ICA_, RAN4, IRAN0, CLDFLAG, NRANDO, &
       LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ, &
       CLDSTORE, SWSTORE, TYPSTORE, RANSEED, &
       cloud_init
  save
  !//-----------------------------------------------------------------------

! The new cloud-jx will be placed here. Until f90 update is done, this
! routine will only contain parameters and variables for cloud2.

contains

  !//-----------------------------------------------------------------------
  subroutine cloud_init()
    !//---------------------------------------------------------------------
    !// Sets up necessary variables for cloud routines.
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    call RANSET(NRAN_, RAN4, RANSEED)
    !//---------------------------------------------------------------------
  end subroutine cloud_init
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine RANSET(ND, RAN4L, ISTART)
    !//---------------------------------------------------------------------
    !  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND) 
    !     program RAN3 from Press, based on Knuth
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer,intent(in) :: ND
    !// Input/Output
    integer,intent(inout) :: ISTART
    !// Output
    real(r4), intent(out)   :: RAN4L(ND)

    !// Parameters
    integer, parameter ::  MBIG=1000000000
    integer, parameter ::  MSEED=161803398
    integer, parameter ::  MZ=0
    !real(r4) , parameter ::  FAC=1.e-9_r4
    real(r8), parameter ::  FAC=1.e-9_r8

    !// Locals
    integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
    real(r8) :: RAN8(ND)
    !//---------------------------------------------------------------------
    !---initialization and/or fix of ISEED < 0 -----------------------------
    write(6,'(a,i10)') 'Initialize RAN4 in RANSET. RANSEED: ',ISTART
    MJ = MSEED - abs(ISTART)
    MJ = mod(MJ,MBIG)
    MA(55) = MJ
    MK = 1
    do I=1,54
       II = mod(21*I,55)
       MA(II) = MK
       MK = MJ-MK
       if (MK.lt.MZ) then
          MK=MK+MBIG
       end if
       MJ = MA(II)
    end do
    do K=1,4
       do I=1,55
          MA(I)=MA(I)-MA(1+MOD(I+30,55))
          if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
          end if
       end do
    end do
    INEXT = 0
    INEXTP = 31
    ISTART = 1
    !---generate next ND pseudo-random numbers -----------------------------
    do J=1,ND
       INEXT = mod(INEXT,55) +1
       INEXTP = mod(INEXTP,55) +1
       MJ = MA(INEXT) - MA(INEXTP)
       if (MJ .lt. MZ) then
          MJ=MJ+MBIG
       end if
       MA(INEXT) = MJ
       RAN8(J) = real(MJ, r8) * FAC
    end do
    !// Convert to float
    RAN4L = real(RAN8, r4)
    !//---------------------------------------------------------------------
  end subroutine RANSET
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
end module cloudjx
!//=========================================================================

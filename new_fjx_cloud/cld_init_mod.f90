!---------------------------------------------------------------------------
!        UCI CTM core p-7.1 (1/2013)
!---------------------------------------------------------------------------
! Adjusted to CTM3
! Amund Sovde, February 2015
!---------------------------------------------------------------------------
!MODULE: CLD_INIT
!
!DESCRIPTION: calculate R effective in liquid and ice cloud
!                         liquid water and ice water path
!
!=======================================================================
! IPARW      = met fields global longitude dimension
! JPARW      = met fields global latitude dimension
! IPAR       = window longitude dimension
! JPAR       = window latitude dimension
! LPAR       = window vertical dimension
! LWEPAR     = wet convection vertical dimension
! IDGRD      = degrade parameter (IPARW/IPAR)
! JDGRD      = degrade parameter (JPARW/JPAR)
! LDEG       = .TRUE. if horizontal degradation used
! G0         = 9.80665d0
!=======================================================================
!//-------------------------------------------------------------------------
MODULE CLD_INIT_MOD
  !//-----------------------------------------------------------------------
  USE CMN_SIZE_MOD, ONLY : IPARW, JPARW, IPAR, JPAR, LPAR, LWEPAR,  &
                             IDGRD, JDGRD, NDGRD
  USE CMN_CTM_MOD,  ONLY : LDEG, IMAP, JMAP, ZDEGI, ZDEGJ
  USE CMN_PARAMETERS,  ONLY : G0, CPI180
  !//-----------------------------------------------------------------------
  IMPLICIT NONE
  !//-----------------------------------------------------------------------

  !-------data to set up the random number sequence for use in cloud-JX
  INTEGER, PARAMETER :: NRAN_ = 10007  ! dimension for random number
  REAL*4  :: RAN4(NRAN_)       ! Random number set
  INTEGER :: IRAN0

  !-----------------------------------------------------------------------

  INTEGER :: ICAIJ(IPAR,JPAR)  ! Total number of ICA at I,J
  INTEGER :: NCLDIJ(IPAR,JPAR) ! Total number of PHOTO-JX called at I,J

  REAL*8, DIMENSION(LWEPAR,IPAR,JPAR) :: &
       CLDFRSTORE      & ! Cloud fraction
       ,WLCSTORE        & ! Liquid water content
       ,REFFLSTORE      & ! R effective in liguid cloud
       ,REFFISTORE      & ! R effective in ice cloud
       ,LWPSTORE        & ! Liguid water path
       ,IWPSTORE          ! Ice water path

  !// PUBLIC SUBROUTINES:
  !
  PUBLIC  :: CLOUDW, RANSET
  SAVE

  !//-----------------------------------------------------------------------
CONTAINS
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  SUBROUTINE CLOUDW(CLDFRW,CLDIWCW,CLDLWCW,PW,ETAA,ETAB,ZOFLEW,  &
                        MYEAR,JDAY,GMTAU)
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------

    ! Input parameter
    real*8, intent(in), dimension(IPARW,JPARW,LWEPAR) :: &
         CLDFRW    & ! Cloud Fraction [0,1] - midpoint value
         ,CLDIWCW   & ! Cloud Ice Water Content [Kg/Kg] - midpoint value
         ,CLDLWCW     ! Cloud Liquid Water Content [Kg/Kg] - midpoint value

    real*8, intent(in), dimension(IPARW,JPARW) :: &
         PW          ! Surface pressure (hPa)

    real*8, intent(in), dimension(LPAR+1) :: &
         ETAA      & ! Edge of ETA pressure coordinates
         ,ETAB        !     (etaa + etab*P_sfc)

    real*8, intent(in) :: &
         ZOFLEW(LPAR+1,IPARW,JPARW)   ! Edge of vertical grid (m)

    real*8, intent(in)  :: GMTAU         ! time in a day
    integer, intent(in) :: JDAY          ! day in a year
    integer, intent(in) :: MYEAR         ! year
    real*8, parameter   :: G100 = 100.d0/G0
    real*8, dimension(LWEPAR+1) :: PLWR, ZLWR
    real*8, dimension(LWEPAR) :: PDEL, CLF, WLC, WIC
    real*8, dimension(NDGRD)  :: CLFN, WLCN, WICN, LWPN, IWPN,  &
                                 FACTR, DZN, PMIDN
    real*8 ::  IWCC, PMID, DP, P1, P2, DZ, F1, CF, ZCLF

    integer :: I, J, IX, JX, II, JJ, L, LTOP, N
    !//---------------------------------------------------------------------

    LTOP  = LWEPAR

    if (.not. LDEG)  then

      do J = 1,JPAR
        do I = 1,IPAR

          if (maxval(CLDFRW(I,J,:)) .le. 0.005d0) then
             CLDFRSTORE(:,I,J) = 0.d0
             WLCSTORE(:,I,J) = 0.d0
             IWPSTORE(:,I,J) = 0.d0
             REFFISTORE(:,I,J) = 0.d0
             LWPSTORE(:,I,J) = 0.d0
             REFFLSTORE(:,I,J) = 0.d0
             cycle
          end if

          do L = 1,LTOP+1
             PLWR(L) = ETAA(L) + ETAB(L)*PW(I,J)
             ZLWR(L) = ZOFLEW(L,I,J)
          end do
          do L = 1,LTOP
             CF  = CLDFRW(I,J,L)
             PDEL(L) = PLWR(L) - PLWR(L+1)
             if (CF .gt. 0.005d0) then
                CLF(L) = CF
                WLC(L) = CLDLWCW(I,J,L) / CF
                WIC(L) = CLDIWCW(I,J,L) / CF
             else
                CLF(L) = 0.d0
                WLC(L) = 0.d0
                WIC(L) = 0.d0
             end if
          end do

          !---based on liquid and ice water mass ratio AVG_D IN CLOUD
          !   (divide by CF above)
          !---WL = liquid water ratio (kg-liq/kg-air), WI = ice water ratio
          do L = 1,LTOP

             CLDFRSTORE(L,I,J)  = CLF(L)
             WLCSTORE(L,I,J) = WLC(L)

!---ice clouds
             if (WIC(L) .gt. 1.d-12) then
                IWPSTORE(L,I,J) = 1000.d0*WIC(L)*PDEL(L)*G100     ! g/m2
                IWCC = IWPSTORE(L,I,J) / (ZLWR(L+1)-ZLWR(L))      ! g/m3
                REFFISTORE(L,I,J) = 164.d0 * IWCC**0.23d0
             else
                IWPSTORE(L,I,J) = 0.d0
                REFFISTORE(L,I,J) = 0.d0
             end if

!---water clouds
             if (WLC(L) .gt. 1.d-12) then
                PMID = 0.5d0 * (PLWR(L)+PLWR(L+1))
                F1   = 0.005d0 * (PMID - 610.d0)
                F1   = min(1.d0, max(0.d0, F1))
                LWPSTORE(L,I,J) = 1000.d0*WLC(L)*PDEL(L)*G100     ! g/m2
                REFFLSTORE(L,I,J) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
             else
                LWPSTORE(L,I,J) = 0.d0
                REFFLSTORE(L,I,J) = 0.d0
             end if

          end do

        end do
      end do

    else  ! degrid

      do J = 1,JPAR
        do I = 1,IPAR

          do L = 1,LTOP

            N  = 0
            do JX = 1,JDGRD
              do IX = 1,IDGRD
                N  = N + 1
                II = IMAP(IX,I)
                JJ = JMAP(JX,J)
                P1 = ETAA(L) + ETAB(L)*PW(II,JJ)
                P2 = ETAA(L+1) + ETAB(L+1)*PW(II,JJ)
                DP = P1 - P2
                FACTR(N) = ZDEGI(IX,I) * ZDEGJ(JX,J)
                PMIDN(N) = FACTR(N) * 0.5d0*(P1+P2)
                DZN(N) = FACTR(N) * (ZOFLEW(L+1,II,JJ)-ZOFLEW(L,II,JJ))

                CF  = CLDFRW(II,JJ,L)
                if (CF .gt. 0.005d0) then
                   CLFN(N) = FACTR(N) * CF
                   WLCN(N) = CLDLWCW(II,JJ,L)/CF
                   WICN(N) = CLDIWCW(II,JJ,L)/CF
                else
                   CLFN(N) = 0.d0
                   WLCN(N) = 0.d0
                   WICN(N) = 0.d0
                end if
!---ice clouds
                if (WICN(N) .gt. 1.d-12) then
                   IWPN(N) = 1000.d0*WICN(N)*DP*G100     ! g/m2
                else
                   IWPN(N) = 0.d0
                end if
!---water clouds
                if (WLCN(N) .gt. 1.d-12) then
                   LWPN(N) = 1000.d0*WLCN(N)*DP*G100     ! g/m2
                else
                   LWPN(N) = 0.d0
                end if
              end do
            end do

            CLDFRSTORE(L,I,J) = 0.d0
            WLCSTORE(L,I,J) = 0.d0
            IWPSTORE(L,I,J) = 0.d0
            REFFISTORE(L,I,J) = 0.d0
            LWPSTORE(L,I,J) = 0.d0
            REFFLSTORE(L,I,J) = 0.d0

            DZ   = sum(DZN) / sum(FACTR)
            PMID = sum(PMIDN) / sum(FACTR)
            CF   = sum(CLFN) / sum(FACTR)
            if (CF .gt. 0.005d0) then
              CLDFRSTORE(L,I,J) = CF
              ZCLF = 1.d0 / sum(CLFN)
              CF   = 0.d0
              do N = 1,NDGRD
                 IWPSTORE(L,I,J) = IWPSTORE(L,I,J) + IWPN(N)*CLFN(N)
                 LWPSTORE(L,I,J) = LWPSTORE(L,I,J) + LWPN(N)*CLFN(N)
                 WLCSTORE(L,I,J) = WLCSTORE(L,I,J) + WLCN(N)*CLFN(N)
                 CF   = CF + WICN(N)*CLFN(N)
              end do
              IWPSTORE(L,I,J) = IWPSTORE(L,I,J) * ZCLF
              LWPSTORE(L,I,J) = LWPSTORE(L,I,J) * ZCLF
              WLCSTORE(L,I,J) = WLCSTORE(L,I,J) * ZCLF
              CF   = CF * ZCLF

              if (CF .gt. 1.d-12) then
                 IWCC = IWPSTORE(L,I,J) / DZ                ! g/m3
                 REFFISTORE(L,I,J) = 164.d0 * IWCC**0.23d0
              end if
              if (WLCSTORE(L,I,J) .gt. 1.d-12) then
                 F1   = 0.005d0 * (PMID - 610.d0)
                 F1   = min(1.d0, max(0.d0, F1))
                 REFFLSTORE(L,I,J) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
              end if
            end if

          end do  ! L loop

        end do
      end do
    end if
    !//---------------------------------------------------------------------
  END SUBROUTINE CLOUDW
  !//---------------------------------------------------------------------

!-----------------------------------------------------------------------

!  fastJX version 7.2 (f90) - Prather notes (Jan 2013)

!---calculation of cloud optical depth in FJX-72 !!!
!---    assumes that clouds are 100% if in layer

!   IWP = ice water path (in layer, in cloud) in g/m**2
!   LWP = liquid water path (in layer, in cloud) in g/m**2
!   REFFI = effective radius of ice cloud (microns)
!   REFFL = effective radius of liquid cloud (microns)

!>>>>method for calculating cloud OD (in layer) used by FJX core or
!outside
!>>>>FJX core needs only the _WP and the REFF_
!>>>> note that FJX can use corrrect Q's and density updates if need be.

!   ODL = LWP * 0.75 * 2.10 / REFFL
!   ODI = IWP * 0.75 * 2.00 / (REFFI * 0.917)

!>>>R-effective determined by main code, not FJX

!   REFF determined by user - some recommendations below (from Neu &
!   Prather)
!     REFFI is a simplle function of ice water content IWC (g/m3, 0.0001
!     to 0.1)
!          IWC = IWP / delta-Z (of layer in m, approx OK)
!          REFFI = 50. * (1. + 8.333 * IWC)
!   prefer Heymsfield++ 2003 JAM, log-log fit ext(/m) vs. IWC, Fig B1a,
!   p.1389
!              EXT (/m) = 1.7e-3 * (IWC/0.1)**0.77
!          REFFI = 164. * IWC**0.23     (33 microns at 0.001 --- 164 at
!          1.0)

!          REFFL is a simple function of pressure (PCLD):
!            FACTOR = (PCLD - 610.) / 200.
!            FACTOR = min(1.0, max(0.0, FACTOR)
!          REFFL = 9.60*FACTOR + 12.68*(1.-FACTOR)
!-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  SUBROUTINE RANSET (ND,RAN4L,ISTART)
    !//---------------------------------------------------------------------
    !  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND) 
    !     program RAN3 from Press, based on Knuth
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, parameter ::  MBIG=1000000000
    integer, parameter ::  MSEED=161803398
    integer, parameter ::  MZ=0
    real*4 , parameter ::  FAC=1.e-9
    integer,intent(in)    :: ND
    real*4, intent(out)   :: RAN4L(ND)
    integer,intent(inout) :: ISTART
    integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
!---initialization and/or fix of ISEED < 0
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
!---generate next ND pseudo-random numbers
    do J=1,ND
       INEXT = mod(INEXT,55) +1
       INEXTP = mod(INEXTP,55) +1
       MJ = MA(INEXT) - MA(INEXTP)
       if (MJ .lt. MZ) then
          MJ=MJ+MBIG
       end if
       MA(INEXT) = MJ
       RAN4L(J) = MJ*FAC
    end do

    !//---------------------------------------------------------------------
  END SUBROUTINE RANSET
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
END MODULE CLD_INIT_MOD
!//-------------------------------------------------------------------------

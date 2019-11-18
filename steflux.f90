!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// STE budget
!//=========================================================================
module steflux
  !//-----------------------------------------------------------------------
  !// MODULE: steflux
  !// DESCRIPTION: Routines to calculate horizontal fluxes and STE.
  !//
  !// Contains
  !//   subroutine CHEMFLUX
  !//   subroutine CHEMFLUX_E90
  !//   subroutine dumpuvflux
  !//   subroutine DUMPTRMASS
  !//   subroutine DUMPTRMASS_E90
  !//   subroutine SAVETRMASS
  !//   subroutine STEBGT_CLR
  !//   subroutine STEBGT_WRITE
  !//   subroutine ctm3_pml
  !//   subroutine ctm3_o3scav
  !//   subroutine chemflux_setup
  !//-----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'steflux.f90'
  !//-----------------------------------------------------------------------
  private
  public CHEMFLUX, CHEMFLUX_E90, dumpuvflux, DUMPTRMASS, DUMPTRMASS_E90, &
       SAVETRMASS, STEBGT_CLR, STEBGT_WRITE, ctm3_pml, ctm3_o3scav, &
       chemflux_setup
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine CHEMFLUX(DTADV,L,QFU,QFV,N,RX,CHMFLX)
    !//---------------------------------------------------------------------
    !// Saves the horizontal flux of specified tracer (O3) in troposphere,
    !// but only for tropospheric to tropospheric boxes.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, NPAR
    use cmn_ctm, only: ALFA, BETA
    use cmn_met, only: ZOFLE
    use cmn_chem, only: N_STE, TMASSMIX2MOLMIX
    use cmn_diag, only: LFLXDG
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    INTEGER, intent(in) :: L
    REAL(r8),  intent(in) :: DTADV
    REAL(r8),  intent(in) :: QFU(IPAR+1,JPAR,NPAR,2)
    REAL(r8),  intent(in) :: QFV(IPAR,JPAR+1,NPAR,2)
    !// Output
    REAL(r8), intent(out)::  CHMFLX(IPAR+1,JPAR+1,2)
    INTEGER :: I, J, II, N

    !// Locals
    REAL(r8) ::  QU,QV,Q0F,Q1F,F0,XF,QF0,RX
    logical :: LTROPAIR
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CHEMFLUX'
    !//---------------------------------------------------------------------

    !// N_STE may be either N_LZ or CTM3 O3.
    if (LFLXDG .and. N .ne. N_STE) then
       write(6,'(a)') f90file//':'//subr// &
            ': This program is only '// &
            'written to collect O3 STE fluxes so far'
       stop
    end if

    F0 = RX / TMASSMIX2MOLMIX(N) !// RX is mole/mole vmr for o3strat
    !// U-flux
    do J = 1, JPAR
       do I = 1, IPAR + 1
          QU  = ALFA(I,J,L) * DTADV  !// Air mass flux out
          Q0F = QFU(I,J,N,1)  !// tracer N flux out
          Q1F = QFU(I,J,N,2)  !// tracer N 1 moment flux out

          !// Define air below 4km as tropospheric
          if (I .eq.  IPAR+1) then
             LTROPAIR = 0.5_r8 * (ZOFLE(L,IPAR,J) + ZOFLE(L,1,J)) .lt. 4000._r8
          else if (I .eq. 1) then
             LTROPAIR = 0.5_r8 * (ZOFLE(L,IPAR,J) + ZOFLE(L,1,J)) .lt. 4000._r8
          else
             LTROPAIR = 0.5_r8 * (ZOFLE(L,I,J) + ZOFLE(L,I-1,J)) .lt. 4000._r8
          end if

          !// Count out-of-box flux from troposphere to troposphere:
          !// Note that F0 = mass mixing ratio of isopleth mixing ratio
          if ((F0 * ABS(QU) .ge. (ABS(Q0F) + ABS(Q1F))) .or. LTROPAIR) then
             !// Pure troposphere
             !// Maximum mass mixing ratio in flux is
             !//   ( ABS(Q0F)+ABS(Q1F) ) / ABS(QU)
             !// so if this mass mixing ratio is smaller than F0, air
             !// is tropospheric. Hence, we count the whole flux.
             XF  = 1._r8
             QF0 = Q0F
          else if (F0 * ABS(QU) .lt. (ABS(Q0F) - ABS(Q1F)) ) then
             !// Pure stratospheric: count no flux
             !// Minimum mass mixing ratio in flux is
             !//   ( ABS(Q0F)-ABS(Q1F) ) / ABS(QU)
             !// so if this is larger than F0, air is stratospheric,
             !// and we count no flux.
             XF  = 0._r8
             QF0 = 0._r8
          ELSE
             !// Partly trop/strat
             !// Fraction of tropospheric air can be calculated by
             !// the mixing ratio in the flux. For mean values, if
             !// ABS(Q0F)/ABS(QU) = F0 we have "tropopause".
             !// To find F0 around Q0F, we recognize that F0 must be some
             !// distance X between abs(Q0F)-abs(Q1F) and abs(Q0F)+abs(Q1F):
             !//   abs(Q0F)-abs(Q1F) + X*(2*abs(Q1F)) = F0*QU
             !//   X = (F0*QU - abs(Q0F) + abs(Q1F))/(2*abs(Q1F))
             !// Then we limit abs(Q1F) with 1.D-10*ABS(Q0F) in case Q1F
             !// is approximately zero.
             XF = (F0*ABS(QU) - ABS(Q0F) + ABS(Q1F)) / &
                      (2._r8 * MAX(ABS(Q1F), 1.e-10_r8*ABS(Q0F)))
             !// Now, how to use this fraction:
             !// Tropospheric ozone flux is then the integral over the area
             !// bound by the lower base (q0f-q1f), and upper base (f0*qu),
             !// with height xf:
             QF0 = (F0 * abs(QU) + ABS(Q0F) - ABS(Q1F)) / 2._r8 * XF
             QF0 = SIGN(QF0, Q0F) 
          end if
          CHMFLX(I,J,1) = QF0
       end do
    end do

    !// V-flux
    do J = 1, JPAR + 1
       do I = 1, IPAR

          QV  = BETA(I,J,L) *  DTADV

          Q0F = QFV(I,J,N,1)
          Q1F = QFV(I,J,N,2)

          !// Define air below 4km as tropospheric
          if (J .eq. JPAR + 1) then
             LTROPAIR = ZOFLE(L,I,JPAR) .lt. 4000._r8
          else if (J .eq. 1) then
             LTROPAIR = ZOFLE(L,I,1) .lt. 4000._r8
          else
             LTROPAIR = 0.5_r8*(ZOFLE(L,I,J) + ZOFLE(L,I,J-1)) .lt. 4000._r8
          end if

          !// For comments, see U-flux comments.
          if ( (F0 * ABS(QV) .ge. (ABS(Q0F) + ABS(Q1F))) .or. LTROPAIR) then
             XF  = 1._r8
             QF0 = Q0F
          else if (F0 * ABS(QV) .lt. (ABS(Q0F) - ABS(Q1F)) ) then
             XF  = 0._r8
             QF0 = 0._r8
          else
             XF  = (F0 * ABS(QV) - ABS(Q0F) + ABS(Q1F)) / &
                  (2._r8 * MAX(ABS(Q1F), 1.e-10_r8*ABS(Q0F)))
             QF0 = (F0 * ABS(QV) + ABS(Q0F) - ABS(Q1F)) / 2._r8 * XF

             QF0 = SIGN(QF0,Q0F)             
          end if
          CHMFLX(I,J,2) = QF0
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine CHEMFLUX
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine CHEMFLUX_E90(DTADV,L,QFU,QFV,CHMFLX, NN)
    !//---------------------------------------------------------------------
    !// Finds horizontal fluxes within the troposphere based
    !// on e90-tracer.
    !// NN is tracer number, and should be either N_STE or N_LZ.
    !//
    !// Ole Amund Sovde, February 2012
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, NPAR
    use cmn_ctm, only: ALFA, BETA
    use cmn_met, only: ZOFLE
    use cmn_chem, only: Ne90, N_STE, N_LZ, E90VMR_TP, TMASSMIX2MOLMIX
    use cmn_diag, only: LFLXDG
    use utilities, only: ctmExitC
   !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: L, NN
    real(r8), intent(in)  :: DTADV
    real(r8), intent(in)  :: QFU(IPAR+1,JPAR,NPAR,2)
    real(r8), intent(in)  :: QFV(IPAR,JPAR+1,NPAR,2)
    !// Output
    real(r8), intent(out) :: CHMFLX(IPAR+1,JPAR+1,2)

    !// Locals
    integer :: I, J, II
    real(r8) :: QU,QV,Q0F,Q1F,F0,XF,QF0,RX,Y
    real(r8) :: Q0F_E90,Q1F_E90
    logical :: LTROPAIR
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CHEMFLUX_E90'
    !//---------------------------------------------------------------------


    !// Find flux using e90 tracer
    if (Ne90 .le. 0) & 
         CALL ctmExitC(f90file//':'//subr//': This routine needs e90 tacer')
    !// Find flux using e90 tracer or linoz tracer
    if (.not.(NN .eq. N_STE .or. NN .eq. N_LZ)) &
         CALL ctmExitC(f90file//':'//subr//': This routine needs N_STE/linoz')

    !// mass mixing ratio of e90 tracer at tropopause
    F0 = E90VMR_TP / TMASSMIX2MOLMIX(Ne90)   !RX is mole/mole vmr for o3strat


    !// trapezoid y axis lower (left) base |q0f|-|q1f|
    !//                 upper base (right) |q0f|+|q1f|
    !// trapezoid x axis  from 0 to 1
    !// U-flux
    do J = 1, JPAR
       do I = 1, IPAR + 1
          QU      = ALFA(I,J,L) * DTADV !// Air mass flux out
          Q0F_E90 = QFU(I,J,Ne90,1)  !// tracer N_E90 flux out
          Q1F_E90 = QFU(I,J,Ne90,2)  !// tracer N_E90 1 moment flux out
          Q0F     = QFU(I,J,NN,1)    !// O3 tracer flux out
          Q1F     = QFU(I,J,NN,2)    !// O3 tracer 1 moment flux out

          !// Define air below 4km as tropospheric
          if (I .eq. IPAR + 1) then
             LTROPAIR = 0.5_r8 * (ZOFLE(L,IPAR,J) + ZOFLE(L,1,J)) .lt. 4000._r8
          else if (I .eq. 1) then
             LTROPAIR = 0.5_r8 * (ZOFLE(L,IPAR,J) + ZOFLE(L,1,J)) .lt. 4000._r8
          else
             LTROPAIR = 0.5_r8 * (ZOFLE(L,I,J) + ZOFLE(L,I-1,J)) .lt. 4000._r8
          end if

          !// To calculate the flux, we must decide whether the air is
          !// tropospheric, stratospheric or both. To do this, we apply the
          !// mean flux of E90 tracer out of the grid box (Q0F_E90) and the
          !// first moment (Q1F_E90).
          !// Figuratively, Q0F_E90 and Q1F_E90 set up a trapezoid, with left
          !// top point at ABS(Q0F_E90)+ABS(Q1F_E90) ant right top point at
          !// ABS(Q0F_E90)-ABS(Q1F_E90).
          !// The signs are opposite of the original CHEMFLUX, where the
          !// gradient in O3 is used: E90 and O3 has opposite gradients.
          !//
          !// Stratospheric air is found when the mass mixing ratio of the
          !// E90-tracer is smaller than F0 given above. We have to compare
          !// this with the mass mixing ratio in the flux, and for simplicity
          !// we compare masses:
          !// Tropospheric air:
          !//   Mass in left corner of trapezoid is higher than the
          !//   stratospheric limit F0*ABS(QU):
          !//     F0*ABS(QU) < (ABS(Q0F_E90)-ABS(Q1F_E90))
          !// Stratospheric air:
          !//   Mass in right corner of trapezoid is lower than the
          !//   stratospheric limit F0*ABS(QU):
          !//     F0*ABS(QU) >= (ABS(Q0F_E90)+ABS(Q1F_E90))
          !// Partly tropospheric and partly stratospheric air:
          !//   We need to locate the tropopause first. Having the range of
          !//   the trapezoid x-axis [0,1], the tropopause lies within, at
          !//   a point X. X is therefore the fraction of tropospheric air,
          !//   and lies between abs(Q0F_E90)+abs(Q1F_E90) and
          !//   abs(Q0F_E90)-abs(Q1F_E90):
          !//     abs(Q0F_E90)+abs(Q1F_E90) - X*(2*abs(Q1F_E90)) = F0*QU
          !//     X = -(F0*QU - abs(Q0F_E90) - abs(Q1F_E90))/(2*abs(Q1F_E90))
          !//   Stratospheric fraction is therefore (perhaps more easily
          !//   derived by looking from right to left on the trapezoid):
          !//     XF = (1 - X)
          !//        = (F0*QU - abs(Q0F_E90) + abs(Q1F_E90))/(2*abs(Q1F_E90))
          !//   To avoid dividing by zero, we limit the divisor abs(Q1F_E90)
          !//   with 1.D-10*ABS(Q0F_E90), in case Q1F_E90 is approximately
          !//   zero.
          !//   Again: Note that when using Q0F and Q1F for O3, XF is the
          !//          tropospheric fraction because its gradient is opposite
          !//          of E90.
          !//
          !// We now have XF, and then calculate the O3 flux at (1-XF), which
          !// we call Y. This equation is based on the fact that the O3
          !// gradient is of opposite sign as for E90, and we interpolate
          !// linearly from the left corner of the ozone trapezoid (not E90,
          !// thus opposite signs).
          !//   Y = ABS(Q0F)-ABS(Q1F) + (1-XF) * (2*ABS(Q1F))
          !//
          !// Having the flux (ABS(Q0F)-ABS(Q1F)) at the left corner, and
          !// Y at (1-XF), we need to integrate over the ozone trapezoid
          !// (i.e. finding the area bound by x=[0,1-XF] and
          !// y=[(ABS(Q0F)-ABS(Q1F)),Y]:
          !//   QF0 = (ABS(Q0F) - ABS(Q1F) + Y) * (1.0d0 - XF) / 2.0d0
          !//
          !// If you are not familiar with the fluxes, you may find it strange
          !// that we need to integrate. Q0F is the mean flux out of the box,
          !// and if Q1F=0, the tropospheric part is just (1-XF)*Q0F, which
          !// is actually the area of a trapezoid with Q1F=0.
          !// Given a gradient Q1F, there is less O3 flux to the left
          !// (ABS(Q0F)-ABS(Q1F)) and more to the right (Y), integration
          !// is necessary.

          if ( (F0 * ABS(QU) .lt. (ABS(Q0F_E90) - ABS(Q1F_E90)) ) &
               .or. LTROPAIR ) then
             !// Pure troposphere:
             XF  = 1._r8
             QF0 = Q0F  !// Count O3 flux
          else if (F0 * ABS(QU) .ge. (ABS(Q0F_E90) + ABS(Q1F_E90))) then
             !// Pure stratosphere
             XF  = 0._r8
             QF0 = 0._r8 !// Zero O3 flux
          else
             !// Partly trop/strat
             !// Find fraction of strat air (e90 gradient is opposite of O3)
             XF = (F0 * ABS(QU) - ABS(Q0F_E90) + ABS(Q1F_E90)) / &
                      (2._r8 * MAX(ABS(Q1F_E90), 1.e-10_r8*ABS(Q0F_E90)))
             !// XF is where f0*abs(qu)= tracer flux in between the lower bound
             !// and the upper bound for e90. Solve XF to get the dividing
             !// point of the box for stratospheric/tropospheric for e90.
             !// For e90: 1 to xf is stratospheric part on the right of
             !//  trapezoid.
             !//
             !// The x-axis in the trapezoid (with upper and lower bound of
             !// tracer fluxes) is the fraction of air mass.
             !// We should solve the ozone flux at (1-XF), which
             !// is named Y here. (0 to (1-XF) is tropospheric air.)
             !// And then we integrate the area from (Q0F-Q1F)
             !// to Y (another trapezoid inside the big one).
             Y = (1 - XF) * (2._r8*ABS(Q1F)) + (ABS(Q0F) - ABS(Q1F))
             !// Tropospheric ozone flux is then the integral over the area
             !// bound by the lower base (q0f-q1f), and upper base y,
             !// with height (1-xf).
             QF0 = (ABS(Q0F) - ABS(Q1F) + Y) * (1.0_r8 - XF) / 2._r8
             QF0 = SIGN(QF0, Q0F) 
          end if
          CHMFLX(I,J,1) = QF0
       end do
    end do

    !// V-flux
    do J = 1, JPAR+1
       do I = 1,IPAR

          QV  = BETA(I,J,L) *  DTADV

          Q0F_E90 = QFV(I,J,Ne90,1)
          Q1F_E90 = QFV(I,J,Ne90,2)
          Q0F     = QFV(I,J,NN,1)
          Q1F     = QFV(I,J,NN,2)

          !// Define air below 4km as tropospheric
          if (J .eq. JPAR + 1) then
             LTROPAIR = ZOFLE(L,I,JPAR) .lt. 4000._r8
          else if (J .eq. 1) then
             LTROPAIR = ZOFLE(L,I,1) .lt. 4000._r8
          else
             LTROPAIR = 0.5_r8 * (ZOFLE(L,I,J) + ZOFLE(L,I,J-1)) .lt. 4000._r8
          end if
            
          !// For comments, see U-flux comments.
          if ( (F0 * ABS(QV) .lt. (ABS(Q0F_E90) - ABS(Q1F_E90)) ) &
               .or. LTROPAIR) then
             !// Pure troposphere
             XF  = 1._r8
             QF0 = Q0F
          else if (F0 * ABS(QV) .ge. (ABS(Q0F_E90) + ABS(Q1F_E90))) then
             !// Pure stratosphere
             XF  = 0._r8
             QF0 = 0._r8
          else
             XF  = (F0 * ABS(QV) - ABS(Q0F_E90) + ABS(Q1F_E90)) / &
                    (2._r8 * MAX(ABS(Q1F_E90), 1.e-10_r8*ABS(Q0F_E90)))
             !// See comments on U-flux for info on the following:
             !QF0= Q0F * (1._r8-XF)
             Y = (1 - XF) * (2._r8*ABS(Q1F)) + (ABS(Q0F)-ABS(Q1F))
             QF0 = (ABS(Q0F) - ABS(Q1F) + Y) * (1._r8 - XF) / 2._r8
             QF0 = SIGN(QF0, Q0F) 
          end if

          CHMFLX(I,J,2)= QF0

       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine CHEMFLUX_E90
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine dumpuvflux(uflux,vflux,k) ! for flux less than x ppb
    !//---------------------------------------------------------------------
    !// Accumulate horizontal fluxes into 3D arrays. Accumulates over
    !// time period defined by STE flux calendar.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_diag, only: UFLX, VFLX
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: uflux(ipar+1,jpar,lpar)
    real(r8), intent(in) :: vflux(ipar,jpar+1,lpar)
    integer, intent(in):: k
    !// Locals
    integer :: i,j,l
    !//---------------------------------------------------------------------
    !// note uflx(x,y,1)  sum of uflx below isopleth value 1

    !// dump flux every time after the flux is computed
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR+1
             uflx(I,J,K) = uflx(I,J,K) + uflux(I,J,L)
          end do
       end do
    end do
      
    do L = 1, LPAR
       do J = 1, JPAR+1
          do I = 1, IPAR
             vflx(I,J,K) = vflx(I,J,K) + vflux(I,J,L)
          end do
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine dumpuvflux
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine DUMPTRMASS(LO3PAUZ,K)    
    !//---------------------------------------------------------------------
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: NRMETD, STT
    use cmn_chem, only: N_STE
    use cmn_diag, only: TROPM
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: K
    logical, intent(in) :: LO3PAUZ(LPAR,IPAR,JPAR,2)
    !// Locals
    real(r8) :: PPB, ZNRMETD
    integer :: I,J,L
    !//---------------------------------------------------------------------

    !// CTM3: Have changed N_LZ to N_STE
    if (N_STE .le. 0)  return

    !// Get daily average to avoid diurnal changes and small variability.
    !// This routine has to be placed in met step (8 times a day).
    !// TROPM is reset every day, in SAVETRMASS
    ZNRMETD = 1._r8 / real(NRMETD, r8)
    do L = 1, LPAR
       do J = 1, JPAR
          do I = 1, IPAR
             if (.not.LO3PAUZ(L,I,J,K)) &
                  TROPM(I,J,K) = TROPM(I,J,K) + STT(I,J,L,N_STE) * ZNRMETD
          end do
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine DUMPTRMASS
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine DUMPTRMASS_E90(K,N_TRACER)
    !//---------------------------------------------------------------------
    !// Save daily average of tropospheric O3 mass, using e90-tracer,
    !// through LSTRATAIR_E90.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: NRMETD, STT
    use cmn_chem, only: Ne90, LSTRATAIR_E90
    use cmn_diag, only: TROPM
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: K, N_TRACER
    !// Locals
    real(r8) :: PPB, ZNRMETD
    integer :: I,J,L
    !//---------------------------------------------------------------------

    if (N_TRACER .le. 0)  return  ! *** exit
    if (Ne90 .le. 0)  return  ! *** exit

    !// Get daily average to avoid diurnal changes and small variability.
    !// This routine has to be placed in met step (8 times a day).
    !// TROPM is reset every day, in SAVETRMASS
    ZNRMETD = 1._r8 / real(NRMETD,r8)
    do L = 1, LPAR
       do J = 1,JPAR
          do I = 1,IPAR
             if (.not.LSTRATAIR_E90(L,I,J)) then
                TROPM(I,J,K) = TROPM(I,J,K) + STT(I,J,L,N_TRACER) * ZNRMETD
             end if
          end do
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine DUMPTRMASS_E90
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine SAVETRMASS(N)    
    !//---------------------------------------------------------------------
    !// Saves daily average TROPM in TROPMASS and clears TROPM.
    !// First timestep sets TROPMASS_0 = TROPM
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   Removed loops.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_chem, only: N_STE, N_LZ
    use cmn_diag, only: TROPM, TROPMASS, TROPMASS_0
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in):: N
    !//---------------------------------------------------------------------

    !// Skip if relevant tracers are not included
    if (N_STE .le. 0 .and. N_LZ .le. 0) return

    !// dump daily average to tropmass
    !// clear tropm
    if (N .eq. 0) then
       !// Initialise
       TROPMASS_0(:,:,:) = TROPM(:,:,:)
    end if

    if (N .eq. 1) then
       !// Save current TROPM
       TROPMASS(:,:,:) = TROPM(:,:,:)
       !// Reset TROPM
       TROPM(:,:,:) = 0._r8
    end if

    !//---------------------------------------------------------------------
  end subroutine SAVETRMASS
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine STEBGT_CLR(NDAY,N)
    !//---------------------------------------------------------------------
    !// Clear STE budget arrays; save TROPMASS_0.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: GMTAU, JDATE, JMON, JYEAR
    use cmn_chem, only: N_STE
    use cmn_diag, only: O3PBLSINK, N2OPBLSINK, UFLX, VFLX, O3PML, &
         TROPMASS, TROPMASS_0, O3WSCAV, &
         IDAY0_STE, ITAU0_STE, JDATE0_STE, JMON0_STE, JYEAR0_STE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,N
    !// Locals
    integer :: I, J
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'STEBGT_CLR'
    !//---------------------------------------------------------------------

    if (N_STE.eq.0) return  ! *** exit

    !// save before clear

    if (N .eq. 0) then
       O3PBLSINK(:,:,:)  = 0._r8
       N2OPBLSINK(:,:) = 0._r8
       UFLX(:,:,:) = 0._r8
       VFLX(:,:,:) = 0._r8
       O3PML(:,:,:) = 0.0_r8
       TROPMASS_0(:,:,:) = 0.0_r8
       TROPMASS(:,:,:) = 0.0_r8
       O3WSCAV(:,:,:) = 0.0_r8
       !//RTN2O(:,:) = 0.0_r8
       IDAY0_STE = NDAY
       ITAU0_STE = int(GMTAU)
       JDATE0_STE = JDATE
       JMON0_STE =  JMON
       JYEAR0_STE = JYEAR
    end if

    if (N .eq. 1) then
       O3PBLSINK(:,:,:)  = 0._r8
       N2OPBLSINK(:,:) = 0._r8
       UFLX(:,:,:) = 0._r8
       VFLX(:,:,:) = 0._r8
       O3PML(:,:,:) = 0.0_r8
       O3WSCAV(:,:,:) = 0.0_r8
       !//RTN2O(:,:) = 0._r8

       !// Set TROPMASS_0 to current TROPMASS
       TROPMASS_0(:,:,:) = TROPMASS(:,:,:)
       !// Zero current TROPMASS
       TROPMASS(:,:,:) = 0.0_r8
       !// Date stuff
       IDAY0_STE = NDAY
       ITAU0_STE = int(GMTAU)
       JDATE0_STE = JDATE
       JMON0_STE =  JMON
       JYEAR0_STE = JYEAR
       write(6,'(2(A,I6))') f90file//':'//subr// &
            ': jdate0_STE=',jdate0_STE,', jmon0_STE',jmon
    end if

    !//---------------------------------------------------------------------
  end subroutine STEBGT_CLR
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine STEBGT_WRITE_OLD
    !//---------------------------------------------------------------------
    !// Write STE budget to file.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LLINOZ
    use cmn_ctm, only: IDAY, AREAXY, JYEAR, JMON, JDATE, GMTAU
    use cmn_chem, only: N_STE, o3iso1, o3iso2, E90VMR_TP
    use cmn_diag, only: IDAY0_STE, itau0_STE, jmon0_STE, jdate0_STE, &
         jyear0_STE, NFLX, UFLX, VFLX, O3PML, O3WSCAV, &
         O3PBLSINK, TROPMASS, TROPMASS_0
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Locals
    character(len=28) :: filnm
    character(len=4) :: fryrlab, toyrlab
    character(len=6) :: frdaylab, todaylab
    integer :: itau,ia,ib,i,j,nc,k,NNFX, fnr
    real(r8) :: xflx, yflx, dyear,zdyear, stesum, pmlsum, scavsum
    real(r8) :: ste(ipar,jpar,NFLX), hflx(ipar,jpar), dm(ipar,jpar,NFLX)
    !//---------------------------------------------------------------------

    if (N_STE .le. 0) return  ! *** exit

    dyear = real(IDAY - IDAY0_STE, r8) / 365.0_r8
    write(6,'(2(A,I6,5X),A,1PE12.5)') &
         'iday=',iday,'iday0_STE',iday0_STE,'dyear=',dyear     
    zdyear = 1._r8 / dyear

    do nc = 1, NFLX ! for 2 diff o3 isopleths + e90 + linoz
       
       do j = 1,jpar
          do i = 1,ipar
             xflx = uflx(i+1,j,nc) - uflx(i,j,nc)
             yflx = vflx(i, j+1,nc) - vflx(i,j,nc)
             hflx(i,j) = xflx + yflx
          end do
       end do

       !// This is just to smooth down the noise of flux divergence.
       !// do 1-2-1 filtering for horizontal flux  for k=1,2 twice.
       do k = 1, 2

          do j = 1, jpar
             do i = 1, ipar
                ib = mod(i, ipar) + 1
                ia = mod(i-2+ipar, ipar) + 1
                hflx(i,j) = 0.25_r8 * hflx(ib,j) &
                          + 0.5_r8 * hflx(i,j) &
                          + 0.25_r8 *hflx(ia,j)
             end do
          end do

          do i = 1,ipar
             do j = 6, jpar - 5         
                hflx(i,j) = 0.25_r8 * hflx(i,j-1) &
                          + 0.5_r8 * hflx(i,j) &
                          + 0.25_r8 * hflx(i,j+1)
             end do
          end do

       end do  ! end of k-loop 1-2-1 filter

       !// the last 5 grid points in the polar region 
       !do i = 1,ipar
       !   do j = 1,5
       !      hflx(i,jpar-j+1) = hflx(i,jpar-5)
       !      hflx(i,j) = hflx(i, 6)
       !   end do
       !end do

       !// Change in tropospheric mass
       do i = 1, ipar
          do j = 1, jpar
             dm(i,j,nc) = tropmass(i,j,nc) - tropmass_0(i,j,nc)
          end do
       end do

       do j = 1, jpar
          do i = 1, ipar
             ib = mod(i, ipar) + 1
             ia = mod(i-2+ipar, ipar) + 1
             dm(i,j,nc) = 0.25_r8 * dm(ia,j,nc) &
                        + 0.5_r8 * dm(i,j,nc) &
                        + 0.25_r8 * dm(ib,j,nc)
          end do
       end do

       do j = 2, jpar-1
          do i = 1, ipar
             ib = mod(i,ipar) + 1
             ia = mod(i-2+ipar, ipar) + 1
             dm(i,j,nc) = 0.25_r8 * dm(i,j-1,nc) &
                        + 0.5_r8 * dm(i,j,nc) &
                        + 0.25_r8 * dm(i,j+1,nc)
          end do
       end do

       !// Polar boxes
       do i = 1, ipar
          dm(i,jpar,nc) = dm(i,jpar-1,nc)
          dm(i,1,nc) = dm(i,2,nc)
       end do

       stesum = 0._r8
       pmlsum = 0._r8
       scavsum = 0._r8

       !// Do the trop mass o3 budget of ste by add and subtract all terms
       do j = 1, jpar
          do i = 1, ipar  
             ste(i,j,nc) = dm(i,j,nc) + hflx(i,j) &
                           - o3pml(i,j,nc) - o3pblsink(i,j,nc) &
                           - o3wscav(i,j,nc)

             stesum  = stesum + ste(i,j,nc)
             pmlsum  = pmlsum + o3pml(i,j,nc)
             scavsum = scavsum + o3wscav(i,j,nc)

             !// Convert to grams/m2/year from kg
             ste(i,j,nc) = ste(i,j,nc) * 1.e3_r8 / areaxy(i,j) * zdyear
          end do
       end do

       !// print to std.out
       if (nc .le. 2) then
          write(6,'(A,I3,1P,2(A,ES14.6),A)') &
               'isopleth ', nc, '   total ste=', stesum*1.e-9_r8, &
               ' Tg. PML from CTM3 O3',pmlsum*1.e-9_r8, ' Tg'
       else if (nc .eq. 3) then
          write(6,'(A,ES14.6,A,ES14.6,A)') &
               'below e90 tpz, total ste=', stesum*1.e-9_r8, &
               ' Tg. PML from CTM3 O3',pmlsum*1.e-9_r8, ' Tg'
       else if (nc .eq. 4) then
          if (LLINOZ) write(6,'(A,ES14.6,A,ES14.6,A)') &
               'LZ below e90 tpz, total ste=', stesum*1.e-9_r8, &
               ' Tg. PML from LINOZ',pmlsum*1.e-9_r8, ' Tg'
       end if

       if (nc .eq. 4 .and. LLINOZ) then
          write(6,'(A,i1,A,5(1X,es10.3))') &
               'STE',nc,'/PML/HFLX/DM/PSINK:', &
               stesum*1.e-9_r8*zdyear, &
               pmlsum*1.e-9_r8*zdyear, &
               sum(hflx)*1.e-9_r8*zdyear, &
               sum(dm(:,:,nc))*1.e-9_r8*zdyear, &
               sum(o3pblsink(:,:,nc))*1.e-9_r8*zdyear
       else
          write(6,'(A,i1,A,5(1X,es10.3))') &
               'STE',nc,'/PML/HFLX/DM/SCAV: ', &
               stesum*1.e-9_r8*zdyear, &
               pmlsum*1.e-9_r8*zdyear, &
               sum(hflx)*1.e-9_r8*zdyear, &
               sum(dm(:,:,nc))*1.e-9_r8*zdyear, &
               scavsum*1.e-9_r8*zdyear
       end if

    end do

    !// write out to files and rezero variables
    itau = int(gmtau)   
    write(fryrlab,'(i4.4)') jyear0_STE
    write(frdaylab,'(3i2.2)') jmon0_STE,jdate0_STE,itau0_STE
    print*,'frdaylab=',frdaylab
    write(toyrlab,'(i4.4)') jyear
    write(todaylab,'(3i2.2)')jmon,jdate,itau
    print*,'todaylab=',todaylab
    write(filnm,'(a28)') &
         'stea_'//fryrlab//frdaylab//'z_'//toyrlab//todaylab//'z'

    write(6,'(a,a28)') 'write in files named: ',filnm

    !// Get free file id
    fnr = get_free_fileid()

    open(fnr,file=filnm,form='unformatted',status='unknown')
    write(fnr) NFLX
    write(fnr) o3iso1, o3iso2, E90VMR_TP, E90VMR_TP
    do nc = 1, NFLX
       write(fnr) uflx(:,:,nc), vflx(:,:,nc), &
            tropmass_0(:,:,nc), tropmass(:,:,nc), o3pml(:,:,nc), &
            o3pblsink, ste(:,:,nc), dm(:,:,nc)
    end do
    close(fnr)

    !//---------------------------------------------------------------------
  end subroutine STEBGT_WRITE_OLD
  !//-----------------------------------------------------------------------




  !//-----------------------------------------------------------------------
  subroutine STEBGT_WRITE
    !//---------------------------------------------------------------------
    !// Write STE budget to file.
    !//
    !// Amund Sovde Haslerud, August 2017
    !//   To netCDF4. Skip filtering (better to post process).
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LLINOZ
    use cmn_ctm, only: IDAY, AREAXY, JYEAR, JMON, JDATE, &
         XDGRD, YDGRD, XDEDG, YDEDG
    use cmn_chem, only: N_STE, o3iso1, o3iso2, E90VMR_TP
    use cmn_diag, only: IDAY0_STE, jmon0_STE, jdate0_STE, &
         jyear0_STE, NFLX, UFLX, VFLX, O3PML, O3WSCAV, &
         O3PBLSINK, TROPMASS, TROPMASS_0, &
         RUNTITLE, metTypeInfo, resolutionInfo, &
         nc4deflate_global, nc4shuffle_global
    use cmn_met, only: METTYPE, MET_ROOT, MPATH1,MPATH2
    use netcdf
    use ncutils, only: handle_error
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Locals
    character(len=80) :: filename
    character(len=8) :: datestamp0, datestamp1
    integer :: nc, k
    real(r8) :: dyear, zdyear
    real(r8), dimension(NFLX) :: &
         stesum, pmlsum, scavsum, hflxsum, dmsum, pblsinksum
    real(r8) :: ste(ipar,jpar,NFLX), hflx(ipar,jpar,NFLX), dm(ipar,jpar,NFLX)

    !// netCDF variables
    !//---------------------------------------------------------------------
    integer :: &
         status, &           ! Status for netcdf file 0=OK
         ncid, &             ! file id for output netcdf file
         lon_dim_id, &       ! Dimension id for longitude
         lat_dim_id, &       ! Dimension id for latitude
         ilon_dim_id, &      ! Dimension id for longitude interstices
         ilat_dim_id, &      ! Dimension id for latitude interstices
         nflx_dim_id, &      ! Dimension id for NFLX
         date_size_dim_id, & ! Dimension id for date sizes
         lon_id, &           ! Variable id for longitude
         lat_id, &           ! Variable id for latitude
         ilon_id, &          ! Variable id for longitude interstices
         ilat_id, &          ! Variable id for latitude interstices
         start_time_id, end_time_id, & ! IDs for start/end time
         ndays_id, &
         areaxy_id, &        ! grid area
         isopleths_id, &     ! isopleths
         uflx_id, vflx_id, & ! flux
         tropmass0_id, tropmass_id, &
         pml_id, sink_id, wetscav_id, ste_id, dm_id
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'STEBGT_WRITE'
    !//---------------------------------------------------------------------

    if (N_STE .le. 0) return  ! *** exit

    dyear = real(IDAY - IDAY0_STE, r8) / 365.0_r8
    write(6,'(2(a,i6),a,es12.5)') f90file//':'//subr// &
         ': iday=',iday,', iday0_STE',iday0_STE,', dyear=',dyear     
    zdyear = 1._r8 / dyear

    !// For NFLX (2 diff o3 isopleths + e90 + linoz) calculate hflx and dm.

    !// Change caused by horizontal flux
    hflx(:,:,:) = uflx(2:ipar+1,:,:) - uflx(1:ipar,:,:) &
                  + vflx(:,2:jpar+1,:) - vflx(:,1:jpar,:)
    !// Change in tropospheric mass
    dm(:,:,:) = tropmass(:,:,:) - tropmass_0(:,:,:)

    !// Original UCI code uses a 1-2-1 filter method (0.25:0.5:0.25 weight)
    !// to reduce noise on hflx (zonal direction then meridional direction,
    !// repeated / done twice) and dm (zonal direction then meridional
    !// direction, done once).
    !// This should be done in post processing. Here we calculate global STE
    !// without filters. The filters are described in the manual.

    !// Find STE
    ste(:,:,:) = dm(:,:,:) + hflx(:,:,:) &
                 - o3pml(:,:,:) - o3pblsink(:,:,:) - o3wscav(:,:,:)

    !// print to std.out part 1
    do nc = 1, NFLX

       !// Sums in Tg
       stesum(nc)     = sum(ste(:,:,nc)) * 1.e-9_r8
       pmlsum(nc)     = sum(o3pml(:,:,nc)) * 1.e-9_r8
       scavsum(nc)    = sum(o3wscav(:,:,nc)) * 1.e-9_r8
       hflxsum(nc)    = sum(hflx(:,:,nc)) * 1.e-9_r8
       dmsum(nc)      = sum(dm(:,:,nc)) * 1.e-9_r8
       pblsinksum(nc) = sum(o3pblsink(:,:,nc)) * 1.e-9_r8

       if (nc .le. 2) then
          write(6,'(A,I3,1P,2(A,ES14.6),A)') &
               'isopleth ', nc, '   total ste=', stesum(nc), &
               ' Tg. PML from CTM3 O3',pmlsum(nc), ' Tg'
       else if (nc .eq. 3) then
          write(6,'(A,ES14.6,A,ES14.6,A)') &
               'below e90 tpz, total ste=', stesum(nc), &
               ' Tg. PML from CTM3 O3',pmlsum(nc), ' Tg'
       else if (nc .eq. 4) then
          if (LLINOZ) write(6,'(A,ES14.6,A,ES14.6,A)') &
               'LZ below e90 tpz, total ste=', stesum(nc), &
               ' Tg. PML from LINOZ',pmlsum(nc), ' Tg'
       end if
    end do !// do nc = 1, NFLX

    !// print to std.out part 2
    do nc = 1, NFLX
       write(6,'(a3,i1,a,6(1x,es10.3))') &
            'STE',nc,'/PML/HFLX/DM/WSCV/PBL: ', &
            stesum(nc) * zdyear, &
            pmlsum(nc) * zdyear, &
            hflxsum(nc) * zdyear, &
            dmsum(nc) * zdyear, &
            scavsum(nc) * zdyear, &
            pblsinksum(nc) * zdyear
    end do !// do nc = 1, NFLX


    !// Save to file
    !//---------------------------------------------------------------------
    write(datestamp0(1:8),'(i4.4,2i2.2)') jyear0_STE,jmon0_STE,jdate0_STE
    write(datestamp1(1:8),'(i4.4,2i2.2)') jyear,jmon,jdate
    filename = 'ste_'//datestamp0//'_'//datestamp1//'.nc'
    !//---------------------------------------------------------------------
    !// Open netCDF4 file for writing
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title', &
         'Oslo CTM3 Stratosphere-Troposphere Exchange (STE)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    status=nf90_put_att(ncid,nf90_global,'model_info', &
         'Oslo CTM3 is a 3D Chemical Transport Model')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': model_info')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology',metTypeInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology_path',MPATH1)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology_path')
    status=nf90_put_att(ncid,nf90_global,'resolution_info',resolutionInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': resolution_info')

    status=nf90_put_att(ncid,nf90_global,'runtitle',trim(runtitle))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': runtitle')
    status=nf90_put_att(ncid,nf90_global,'contact_info', &
         'For errors, contact CICERO')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': contact_info')

    !// Define dimensions
    !//---------------------------------------------------------------------
    !// Define lat/lon/lev
    status = nf90_def_dim(ncid,'lat',JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,'lon',IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')

    !// Define ilat/ilon/ilev
    status = nf90_def_dim(ncid,'ilat',JPAR+1,ilat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,'ilon',IPAR+1,ilon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//':define ilon dim')


    !// Define STE fluxes
    status = nf90_def_dim(ncid,'NFLX',NFLX,nflx_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NFLX dim')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,'date_size',6,date_size_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')


    !// Define variables
    !//---------------------------------------------------------------------
    !// Defining the lon variable
    status = nf90_def_var(ncid,'lon',nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    status = nf90_put_att(ncid,lon_id,'units','degrees east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lon_id,'description','Value at grid box center.')

    !// Defining the lat variable
    status = nf90_def_var(ncid,'lat',nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    status = nf90_put_att(ncid,lat_id,'units','degrees north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lat_id,'description','Value at grid box center.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description lat')


    !// Defining ilon variable (lon on interstices)
    status = nf90_def_var(ncid,'ilon',nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_put_att(ncid,ilon_id,'units','degrees east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilon_id,'description', &
         'Value at eastern edge of grid box.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilon')

    !// Defining ilat variable (lat on interstices)
    status = nf90_def_var(ncid,'ilat',nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    status = nf90_put_att(ncid,ilat_id,'units','degrees north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilat_id,'description', &
         'Value at southern edge of grid box.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ilat')


    !// Start time for accumulating data - START_TIME
    status = nf90_def_var(ncid, 'START_TIME', nf90_int, &
         date_size_dim_id, start_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_TIME variable')
    status = nf90_put_att(ncid,start_time_id,'description', &
         'Start date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_TIME')

    !// End time for accumulating data - END_TIME
    status = nf90_def_var(ncid,'END_TIME', nf90_int, &
         date_size_dim_id, end_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_TIME variable')
    status = nf90_put_att(ncid,end_time_id,'description', &
         'End date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_TIME')

    !// Put out NDAYS and IYEAR to indicate what kind of spinup
    status = nf90_def_var(ncid,"NDAYS",nf90_int,ndays_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NDAYS variable')
    status = nf90_put_att(ncid,ndays_id,'description', &
         'Numper of days since model start (NDAY-NDAYI+1)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description NDAYS')

    !// O3 isopleths
    status = nf90_def_var(ncid, 'isopleths', nf90_double, &
         nflx_dim_id, isopleths_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define isopleths variable')
    status = nf90_put_att(ncid,isopleths_id,'description', &
         'Isopleth values (vmr) for flux calculations')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description isopleths')
    status = nf90_put_att(ncid,isopleths_id,'description1',&
         'NFLX=1: STE through O3 isopleth 1')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description1 isopleths')
    status = nf90_put_att(ncid,isopleths_id,'description2',&
         'NFLX=2: STE through O3 isopleth 2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description2 isopleths')
    status = nf90_put_att(ncid,isopleths_id,'description3',&
         'NFLX=3: STE using E90 tropopause (E90 isopleth 3)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description3 isopleths')
    status = nf90_put_att(ncid,isopleths_id,'description4',&
         'NFLX=4: STE using E90 tropopause (E90 isopleth 4)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description4 isopleths')


    !// gridarea
    status = nf90_def_var(ncid, 'gridarea', nf90_double, &
         (/lon_dim_id, lat_dim_id/), areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')


    !// UFLX
    status = nf90_def_var(ncid, 'UFLX', nf90_float, &
         (/ilon_dim_id,lat_dim_id,nflx_dim_id/), uflx_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define UFLX variable')
    status = nf90_def_var_deflate(ncid,uflx_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate UFLX variable')
    status = nf90_put_att(ncid,uflx_id,'longname','Zonal flux of O3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname VFLX')
    status = nf90_put_att(ncid,uflx_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit UFLX')
    status = nf90_put_att(ncid,uflx_id,'description', &
         'In 2D, for latitude j and longitude i: '//&
         'hflx(j,i) = uflx(j,i+1) - uflx(j,i) '// &
                   '+ vflx(j+1,i) - vflx(j,i)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description UFLX')

    !// VFLX
    status = nf90_def_var(ncid, 'VFLX', nf90_float, &
         (/lon_dim_id,ilat_dim_id,nflx_dim_id/), vflx_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define VFLX variable')
    status = nf90_def_var_deflate(ncid,vflx_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate VFLX variable')
    status = nf90_put_att(ncid,vflx_id,'longname','Meridional flux of O3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname VFLX')
    status = nf90_put_att(ncid,vflx_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit VFLX')
    status = nf90_put_att(ncid,vflx_id,'description', &
         'In 2D, for latitude j and longitude i: '//&
         'hflx(j,i) = uflx(j,i+1) - uflx(j,i) '// &
                   '+ vflx(j+1,i) - vflx(j,i)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description VFLX')

    !// TROPMASS_0
    status = nf90_def_var(ncid, 'TROPMASS_0', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), tropmass0_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define TROPMASS_0 variable')
    status = nf90_def_var_deflate(ncid,tropmass0_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate TROPMASS_0 variable')
    status = nf90_put_att(ncid,tropmass0_id,'longname', &
         'Tropospheric mass of O3 at START_TIME')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname TROPMASS_0')
    status = nf90_put_att(ncid,tropmass0_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit TROPMASS_0')

    !// TROPMASS
    status = nf90_def_var(ncid, 'TROPMASS', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), tropmass_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define TROPMASS variable')
    status = nf90_def_var_deflate(ncid,tropmass_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate TROPMASS variable')
    status = nf90_put_att(ncid,tropmass_id,'longname', &
         'Tropospheric mass of O3 at END_TIME')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname TROPMASS')
    status = nf90_put_att(ncid,tropmass_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit TROPMASS')

    !// Prod minus loss (PML)
    status = nf90_def_var(ncid, 'PML', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), pml_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define PML variable')
    status = nf90_def_var_deflate(ncid,pml_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate PML variable')
    status = nf90_put_att(ncid,pml_id,'longname', &
         'Chemical production minus loss in troposphere.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname PML')
    status = nf90_put_att(ncid,pml_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit PML')

    !// PBL sink
    !// (drydep, should be zero for Oslo chemistry, but not for LINOZ diag)
    status = nf90_def_var(ncid, 'PBL_SINK', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), sink_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define PBL_SINK variable')
    status = nf90_def_var_deflate(ncid,sink_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate PBL_SINK variable')
    status = nf90_put_att(ncid,sink_id,'longname', &
         'Tropospheric (PBL) sink of O3 (dry deposition).')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname PBL_SINK')
    status = nf90_put_att(ncid,sink_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit PBL_SINK')

    !// WETSCAV (wet scavenging)
    status = nf90_def_var(ncid, 'WETSCAV', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), wetscav_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define WETSCAV variable')
    status = nf90_def_var_deflate(ncid,wetscav_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate WETSCAV variable')
    status = nf90_put_att(ncid,wetscav_id,'longname', &
         'Tropospheric wet scavenged O3.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname WETSCAV')
    status = nf90_put_att(ncid,wetscav_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit WETSCAV')

    !// STE
    status = nf90_def_var(ncid, 'STE', nf90_float, &
         (/lon_dim_id,lat_dim_id,nflx_dim_id/), ste_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define STE variable')
    status = nf90_def_var_deflate(ncid,ste_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate STE variable')
    status = nf90_put_att(ncid,ste_id,'longname', &
         'STE of O3.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute longname STE')
    status = nf90_put_att(ncid,ste_id,'unit','kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit STE')

    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------
    
 
    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')

    !// Start time
    status = nf90_put_var(ncid, start_time_id, &
         (/jyear0_STE, jmon0_STE, jdate0_STE, 0, 0, 0/) )
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_TIME')
    !// End time
    status = nf90_put_var(ncid, end_time_id, &
         (/jyear, jmon, jdate, 0, 0, 0/) )
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_TIME')
    !// NDAYS
    status = nf90_put_var(ncid,ndays_id,IDAY - IDAY0_STE)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NDAYS')

    !// O3 isopleths
    status = nf90_put_var(ncid, isopleths_id, &
         (/o3iso1, o3iso2, E90VMR_TP, E90VMR_TP/) )
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting isopleths')

    !// gridarea
    status = nf90_put_var(ncid, areaxy_id, AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting gridarea')

    !// UFLX
    status = nf90_put_var(ncid, uflx_id, uflx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting UFLX')

    !// VFLX
    status = nf90_put_var(ncid, vflx_id, vflx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting VFLX')

    !// TROPMASS_0
    status = nf90_put_var(ncid, tropmass0_id, tropmass_0)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting TROPMASS_0')

    !// TROPMASS
    status = nf90_put_var(ncid, tropmass_id, tropmass)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting TROPMASS')

    !// PML
    status = nf90_put_var(ncid, pml_id, O3PML)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting PML')

    !// PBL_SINK
    status = nf90_put_var(ncid, sink_id, o3pblsink)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting PBL_SINK')

    !// WETSCAV
    status = nf90_put_var(ncid, wetscav_id, o3wscav)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting WETSCAV')

    !// STE
    status = nf90_put_var(ncid, ste_id, ste)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting STE')

    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))

    write(6,'(a)') &
         f90file//':'//subr//': wrote '//trim(filename)
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !//---------------------------------------------------------------------
  end subroutine STEBGT_WRITE
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine ctm3_pml(BTT,BTTBCK,AIRB,LSTRATAIR_E90B,LO3PAUZ,MP)
    !//---------------------------------------------------------------------
    !// From BTT and BTTBCK, calculate chemical
    !// tendency (PML) below pre-defined O3-surfaces.
    !//
    !// Skip uppermost layers where O3 may be less than the isopleths.
    !//
    !// Ole Amund Sovde, February 2012
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: N_STE
    use cmn_diag, only: O3PML, LFLXDG
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8),dimension(LPAR,NPAR,IDBLK,JDBLK),intent(in) :: BTT, BTTBCK
    real(r8),dimension(LPAR,IDBLK,JDBLK),intent(in) :: AIRB
    integer, intent(in) :: MP
    logical, intent(in) :: LSTRATAIR_E90B(LPAR,IDBLK,JDBLK)
    logical, intent(in) :: LO3PAUZ(LPAR,IPAR,JPAR,2)
    !// Locals
    real(r8) :: PMLNET
    integer :: I,J,II,JJ,L
    !//---------------------------------------------------------------------

    if (.not.LFLXDG) return

    !// Skip if N_STE is not set
    if (N_STE .le. 0) return

    !// Calculate chemical tendency
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ = J - MPBLKJB(MP) + 1

       do I = MPBLKIB(MP), MPBLKIE(MP)
          II = I - MPBLKIB(MP) + 1

          !// Trop/strat air is defined by LO3PAUZ and LSTRATAIR_E90

          do L = 1, LPAR

             !// Find chemical tendency
             PMLNET = BTT(L,N_STE,II,JJ) - BTTBCK(L,N_STE,II,JJ)

             !// O3 isopleth 1
             if (.not.LO3PAUZ(L,I,J,1)) O3PML(I,J,1) = PMLNET + O3PML(I,J,1)

             !// O3 isopleth 2
             if (.not.LO3PAUZ(L,I,J,2)) O3PML(I,J,2) = PMLNET + O3PML(I,J,2)

             !// Count PML in e90 troposphere
             if (.not.LSTRATAIR_E90B(L,II,JJ)) O3PML(I,J,3) = PMLNET + O3PML(I,J,3)

             !// Linoz tracer is calculated in separate routine

          end do

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    !//---------------------------------------------------------------------
  end subroutine ctm3_pml
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine ctm3_o3scav(BTT,BTTBCK,AIRB,LSTRATAIR_E90B,LO3PAUZ,MP)
    !//---------------------------------------------------------------------
    !// From BTT and BTTBCK, calculate net scavenged O3
    !// below pre-defined O3-surfaces and below e90 tropopause.
    !//
    !// O3WSCAV is also used in p-cnvw_oc.f
    !//
    !// Ole Amund Sovde, March 2012
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: N_STE
    use cmn_diag, only: LFLXDG, O3PML, O3WSCAV
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8),dimension(LPAR,NPAR,IDBLK,JDBLK),intent(in) :: BTT, BTTBCK
    real(r8),dimension(LPAR,IDBLK,JDBLK),intent(in) :: AIRB
    integer, intent(in) :: MP
    logical, intent(in) :: LSTRATAIR_E90B(LPAR,IDBLK,JDBLK)
    logical, intent(in) :: LO3PAUZ(LPAR,IPAR,JPAR,2)
    !// Locals
    real(r8) :: WNET
    integer :: I,J,II,JJ,L
    !//---------------------------------------------------------------------

    !// Skip if STE is not to be calculated
    if (.not.LFLXDG) return

    !// Skip if N_STE is not set
    if (N_STE .le. 0) return

    !// Calculate stratospheric chemical tendency
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ = J - MPBLKJB(MP) + 1

       do I = MPBLKIB(MP),MPBLKIE(MP)
          II = I - MPBLKIB(MP) + 1

          do L = 1, LPAR

             !// Find chemical tendency
             WNET = BTT(L,N_STE,II,JJ) - BTTBCK(L,N_STE,II,JJ)

             if (.not.LO3PAUZ(L,I,J,1)) O3WSCAV(I,J,1) = WNET + O3WSCAV(I,J,1)
             if (.not.LO3PAUZ(L,I,J,2)) O3WSCAV(I,J,2) = WNET + O3WSCAV(I,J,2)

             if (.not.LSTRATAIR_E90B(L,II,JJ)) O3WSCAV(I,J,3) = WNET + O3WSCAV(I,J,3)

          end do

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    !//---------------------------------------------------------------------
  end subroutine ctm3_o3scav
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine chemflux_setup()
    !//---------------------------------------------------------------------
    !// Reads tables/polar_o3loss.dat
    !// to set variables E90VMR_TP, O3iso1 and O3iso2
    !// when LINOZ is *not* used.
    !//
    !// Ole Amund Sovde, February 2012
    !//---------------------------------------------------------------------
    use cmn_size, only: LLINOZ
    use cmn_chem, only: E90VMR_TP, O3iso1, O3iso2, INFILE_POLAR_O3LOSS
    use utilities, only: get_free_fileid
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=80) :: TITL1, INITFI
    integer :: fnr
    !//---------------------------------------------------------------------
    !// Just in case
    if (LLINOZ) return

    !// Read parameters from Linoz-file
    INITFI = INFILE_POLAR_O3LOSS

    !// Get free file id
    fnr = get_free_fileid()

    open(fnr,file=INITFI,status='old',form ='formatted')
    read(fnr,*)  TITL1
    read(fnr,*)  TITL1
    read(fnr,*)  TITL1
    read(fnr,*)  TITL1
    read(fnr,'(8e10.3)')  E90VMR_TP, O3iso1, O3iso2
    write(6,'(a,es10.3)') &
         ' Recommended tropopause threshold mix. ratio for E90 tracer =',E90VMR_TP
    write(6,'(a,2es10.3)') &
         ' STEFLUX isopleth 1/2 vmr =',O3iso1,O3iso2
    close(fnr)
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    !//---------------------------------------------------------------------
  end subroutine chemflux_setup
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module steflux
!//=========================================================================

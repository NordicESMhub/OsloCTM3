C---(p-vect3.f)----v 6.1  (mjp 7/2011)
c---new SOM Lim=3 allows small over/undershoot see C3MAX & C3MIN
c
c---subroutines:  QVECT3
c Version: qcode_61a; 20120221
c
c-----------------------------------------------------------------------
c---Second-Order Moments (SOM) advection code:  M.J. Prather (2007)
c---   This subroutine is written for the X advection step, but commutation
c---   of the order of the moments allows it to be use for Y and Z steps.
C---
c---   This code was initally written in intrinsic vector code (no DO loops)
c---   but this runs much more slowly on most computers (AIX,Intel,...)

c-----------------------------------------------------------------------
      subroutine QVECT3 (QM,QU0,NQ,LIMTR, MSTEP,
     &         QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ, Q0F,Q1F)
c-----------------------------------------------------------------------

c--Computes SOM advection for a single cyclic pipe of tracer in the X-direction
c-----  =QU(1)=> BOX[1] =QU(2)=> ... =QU(NQ)=> BOX[NQ] =QU(1)=> BOX[1] ...
c----LIMTR=0  no limits, neg. tracer allowed - OK
c----LIMTR=1  pos. definite, simple SOM fix (per Prather 1986):adjust QXT & QXX
c----LIMTR=2  pos. definite & monotonic (ca. 1994-95): adjust QXT,QXX,QXY,QXZ
c----LIMTR=3  start with LIMTR=2 and use min-max in parents to limit daughter
c----             new v6.1 allows for +-3% over/undershoot of max/min
c--NOTE the order in which the Limiters are applied:
c---     LIMTR = 1 & 2 are applied ONLY to the X-direction and just before
c---      the advection step.  The moments coming out of QVECT3 are NOT limited
c---     LIMTR = 3 applies LIMTR=2 before step, and the min-max criteria after.

c---Q0F, Q1F are the tracer & 1st-moments fluxes across boundary [J-1] => [J]
c---        carried by the air mass QU
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use utilities, only: ctmExitC
      implicit none
c-----------------------------------------------------------------------
c input parameters
      integer, intent(in) ::
     &     NQ,             ! # of boxes, assumed to be cyclic
     &     LIMTR           ! LIMITR #
      real(r8),  intent(in), dimension(NQ+1) ::
     &     QU0             ! air mass flux to be moved [I-1]->[I] in adv step
c input/output parameters
      real(r8),  intent(inout), dimension(NQ+1) ::
     &     QM,             ! air mass in box at start,
     &     QTT,            ! tracer mass in box [I]
     &     QXT, QXX,       ! 1st & 2nd moments of tracer in X-direction
     &     QXY, QXZ,       ! coupled 1st moments including X component
     &     QYZ, QYT, QYY, QZT, QZZ     !** other SOMs without X
c output parameters
      integer, intent(out) ::
     &      MSTEP          ! no. of multi-steps need to complete step with flux
c                            QU0 and not exceed exceed CFL limit (no neg. mass!)
      real(r8),  intent(out), dimension(NQ+1) ::
     &      Q0F,Q1F        ! computed tracer flux & x-moments from [I-1] to [I]

c local parameters
      integer J, NJ, JJ, JJM1, JP1, JJP1, M

      real(r8), dimension(2050) ::
     &   QMFIN, QM1, QU

      real(r8), dimension(1025) ::
     &   FTT, FXX, FYY, FZZ, FXT, FYT, FZT, FXY, FXZ, FYZ, EPS
     &  ,A0M, A0TT, A0XT, A0XX, A0YT, A0YY, A0ZT, A0ZZ, A0XY, A0XZ, A0YZ
     &  ,B0M, B0TT, B0XT, B0XX, B0YT, B0YY, B0ZT, B0ZZ, B0XY, B0XZ, B0YZ
     &  ,A1M, A1TT, A1XT, A1XX, A1YT, A1YY, A1ZT, A1ZZ, A1XY, A1XZ, A1YZ
     &  ,B1M, B1TT, B1XT, B1XX, B1YT, B1YY, B1ZT, B1ZZ, B1XY, B1XZ, B1YZ
     &  ,ABXMX, ABXMN, ABYMX, ABYMN, ABZMX, ABZMN,   UAB

      real(r8) :: EPS1T, EPSSQ, EPS1S, F0, F1, F2, ABSCAL,
     &     QMMIN, QM1A, QM1B

      logical LGC_AB(1025), LGC_LX, LGC_L2

c---some constants
      real(r8), parameter :: C000 = 1.e-30_r8
      real(r8), parameter :: C033 = 1.0_r8/3.0_r8
      real(r8), parameter :: C034 = 0.33334_r8
      real(r8), parameter :: C149 = 1.49999_r8

c---CXYZ used to limit cross moments
      real(r8), parameter :: CXYZ = 1.5_r8

c---C3MAX & C3MIN used to allow limited under/overshoot in Limiter=3
      real(r8), parameter :: C3MAX = 1.01_r8
      real(r8), parameter :: C3MIN = 0.99_r8

c-----------------------------------------------------------------------
C---REQUIRE that NQ be even (easy if adding non-cyclic columns, but not if true
C---  odd-numbered, cyclic solution) - cannot do separation otherwise.
c---N.B. QVECT3 CANNOT DO CYCLIC PIPEFLOW WITH AN ODD NUMBER OF BOXES!

      if (mod(NQ,2) .ne. 0) call ctmEXITC('NQ must be even in QVECT3')
      NJ = NQ/2

      if (NJ .gt. 1024)
     &     call ctmEXITC('NQ out of bounds(1024) in QVECT3')

C--Assume that always cyclical and add box #1 as #NQ+1.
C--    this is easy even if LCYCLE = .F., just a zero connecting velocity.

C--Now first prepare a column of 'advects' A0 => B0
C--  if U>0:  A1 & B1 have all moments stored
C--  if U<0:  switch A & B, reverse sign of X moments
C---  NB first-order X-moments switch sign for reverse advection


c---QnF() are the fluxes (+ 1st & 2nd moments) of tracer carried by QU()
c---     currently it is not practical to calculated Q2F.
      do J = 1,NQ+1
        Q0F(J) = 0._r8
        Q1F(J) = 0._r8
      enddo

C---UCI CTM uses a new approach to determining the time-step limits.
c---   (1) a GLOBAL limit is derived based only on divergence
c---     i.e., maximum time step 'dt-div' such that the sequence of
c---     operator-split advection (W-V-U in our case) does not remove
c---     more than 95% (a parameter) of the mass in any grid box.
c---     This is our global-divergence limiter
c---   (2) for each pipe-flow advection in QVECT3, check if traditional
c---     CFL limit is exceeded (e.g., c * delta-t > delta-x).
c---     If it is, then each pipe determines its own number of multi-steps
c---     (an integer number of smaller dt) and does a quick internal looping

c---This new combination CFL limiter greatly reduces the number of advection
c---     calls as well as the cost of each step (if done internally in QVECT3).
c---     The multi-stepping works well in jet and polar regions, limiting the
c---     use of a short time to only those levels and latitudes.

C--FIRST check that for global time-step dt-div, the divergence limit is met
          QMMIN = 1._r8
        do J=1,NQ
          QMFIN(J) = QM(J) + (QU0(J) - QU0(J+1))
          QMMIN = min(QMMIN,QM(J),QMFIN(J))
        enddo
         if (QMMIN .le. 0._r8) then
            do J=1,NQ
               write(*,'(i3,5es13.5)') J,QMFIN(J),QM(J),QU0(J),QU0(J+1),
     &              QM(J)+(QU0(J)-QU0(J+1))
            enddo
           call ctmEXITC(' NEG air mass begin-end in QVECT3')
         endif
C--SECOND check for need to multi-step
        do J=2,NQ,2
         QM1(J-1) = QM(J-1) - QU0(J)
         QM1(J)   = QM(J)   + QU0(J)
        enddo
         QMMIN = 1._r8
        do J=1,NQ
         QM1A = QM1(J)/QMFIN(J)
         QM1B = QM1(J)/QM(J)
         QMMIN = min (QM1A, QM1B, QMMIN)
        enddo
      if (QMMIN .gt. 0.01_r8) then
c---no multi-stepping needed to assure non-negative intermediate steps
        MSTEP = 1
       do J=1,NQ+1
        QU(J) = QU0(J)
       enddo
      else
c---must multi-step to avoid negative intermediate steps
        MSTEP = int((1._r8 - QMMIN)/0.99_r8) + 1
       do J=1,NQ+1
        QU(J) = QU0(J) / real(MSTEP, r8)
       enddo
      endif


c>>>>>>begin primary multi-stepping loop
      do M = 1,MSTEP

C---#1 NOW load the EVEN pairs:
      do J = 1,NJ
        LGC_AB(J) = QU(J+J) .ge. 0._r8
        UAB(J)    = abs(QU(J+J))
      enddo

C---  load the EVEN pair of Q--_s into A--_s and B--_s
      do J = 1,NJ
        JJ    = J + J
        JJM1  = J + J - 1
        if (LGC_AB(J))  then
          A0M(J)  = QM(JJM1)
          A0TT(J) = QTT(JJM1)
          A0XX(J) = QXX(JJM1)
          A0YT(J) = QYT(JJM1)
          A0YY(J) = QYY(JJM1)
          A0ZT(J) = QZT(JJM1)
          A0ZZ(J) = QZZ(JJM1)
          A0YZ(J) = QYZ(JJM1)
          A0XT(J) = QXT(JJM1)
          A0XY(J) = QXY(JJM1)
          A0XZ(J) = QXZ(JJM1)
C---  same for the B_s
          B0M(J)  = QM(JJ)
          B0TT(J) = QTT(JJ)
          B0XX(J) = QXX(JJ)
          B0YT(J) = QYT(JJ)
          B0YY(J) = QYY(JJ)
          B0ZT(J) = QZT(JJ)
          B0ZZ(J) = QZZ(JJ)
          B0YZ(J) = QYZ(JJ)
          B0XT(J) = QXT(JJ)
          B0XY(J) = QXY(JJ)
          B0XZ(J) = QXZ(JJ)
         else
          A0M(J)  = QM(JJ)
          A0TT(J) = QTT(JJ)
          A0XX(J) = QXX(JJ)
          A0YT(J) = QYT(JJ)
          A0YY(J) = QYY(JJ)
          A0ZT(J) = QZT(JJ)
          A0ZZ(J) = QZZ(JJ)
          A0YZ(J) = QYZ(JJ)
          A0XT(J) =-QXT(JJ)
          A0XY(J) =-QXY(JJ)
          A0XZ(J) =-QXZ(JJ)
C---  same for the B_s
          B0M(J)  = QM(JJM1)
          B0TT(J) = QTT(JJM1)
          B0XX(J) = QXX(JJM1)
          B0YT(J) = QYT(JJM1)
          B0YY(J) = QYY(JJM1)
          B0ZT(J) = QZT(JJM1)
          B0ZZ(J) = QZZ(JJM1)
          B0YZ(J) = QYZ(JJM1)
          B0XT(J) =-QXT(JJM1)
          B0XY(J) =-QXY(JJM1)
          B0XZ(J) =-QXZ(JJM1)
        endif
      enddo

C---LIMTR = 1:  positive definite, reset Sxx to allow maximum Sx (1.5*So)
      if (LIMTR.eq.1) then
        do J = 1,NJ
            F0 = A0TT(J)
            F1 = min(+C149*F0, max(-C149*F0, A0XT(J)))
            F2 = min((F0+F0-C034*abs(F1)), max(abs(F1)-F0, A0XX(J)))
          A0XT(J) = F1
          A0XX(J) = F2
          A0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XY(J)))
          A0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XZ(J)))
        enddo

C--------do NOT impose limits on downwind B boxes since they are not transported
C       B0XT = min(+C149*B0TT, max(-C149*B0TT,B0XT))
C       B0XX = min((B0TT+B0TT-C034*abs(B0XT)), max(abs(B0XT)-B0TT, B0XX))
C       B0XY = min(+C149*B0TT,max(-C149*B0TT,B0XY))
C       B0XZ = min(+C149*B0TT,max(-C149*B0TT,B0XZ))

      endif

C---LIMTR = 2:  positive and monotonic in box (<X>-dim only)
c---LIMTR = 3:  min-max, needs to have positive-monotonic also
      if (LIMTR.ge.2) then
        do J = 1,NJ
          LGC_L2  = A0XX(J) .gt. 0._r8
          LGC_LX  = A0XT(J) .gt. 0._r8
            F0      = A0TT(J)
         if (LGC_L2) then
c---A0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(A0XX(J), abs(A0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(A0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---A0XX(J) < 0:  curved down at ends, allow if A0XT(J) < A0TT(J)
            F1 = min(+F0, max(-F0, A0XT(J)))
            F2 = max(A0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          A0XT(J) = F1
          A0XX(J) = F2
          A0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XY(J)))
          A0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XZ(J)))
C---no need to impose limits on downwind B since it is not transported
        enddo
      endif



c---LIMTR = 3:  min-max, force positive-monotonic on B boxes also
c---             find and store min/max of both A & B boxes
      if (LIMTR.ge.3) then
        do J = 1,NJ
          LGC_L2  = B0XX(J) .gt. 0._r8
          LGC_LX  = B0XT(J) .gt. 0._r8
            F0      = B0TT(J)
         if (LGC_L2) then
c---B0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(B0XX(J), abs(B0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(B0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---B0XX(J) < 0:  curved down at ends, allow if B0XT(J) < B0TT(J)
            F1 = min(+F0, max(-F0, B0XT(J)))
            F2 = max(B0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          B0XT(J) = F1
          B0XX(J) = F2
          B0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XY(J)))
          B0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XZ(J)))
        enddo

C---LIMTR = 3 min/max in A-- & B-- boxes occurs at end points if monotonic
       do J = 1,NJ
       ABXMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J))/B0M(J)  )
       ABXMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J))/B0M(J)  )
       ABYMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J)+abs(A0XY(J)))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J)+abs(B0XY(J)))/B0M(J)  )
       ABYMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J)-abs(A0XY(J)))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J)-abs(B0XY(J)))/B0M(J)  )
       ABZMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J)+abs(A0XZ(J)))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J)+abs(B0XZ(J)))/B0M(J)  )
       ABZMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J)-abs(A0XZ(J)))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J)-abs(B0XZ(J)))/B0M(J)  )
       enddo

      endif

C---Advect from A0 ==[UAB]==> to B0  for even pairs [1:2] [3:4] [5:6]
      do J = 1,NJ
        EPS(J) = UAB(J) / A0M(J)
        EPS1T  = 1._r8 - EPS(J)
        EPS1S  = EPS1T * EPS1T
        EPSSQ  = EPS(J) * EPS(J)
c------Seperate left A0 box to 2 parts: F_s = section tranferred to B0 box
        FYY(J)  = EPS(J) * A0YY(J)
        FZZ(J)  = EPS(J) * A0ZZ(J)
        FYZ(J)  = EPS(J) * A0YZ(J)
        FYT(J)  = EPS(J) * (A0YT(J) + EPS1T * A0XY(J))
        FZT(J)  = EPS(J) * (A0ZT(J) + EPS1T * A0XZ(J))
        FXY(J)  = EPSSQ * A0XY(J)
        FXZ(J)  = EPSSQ * A0XZ(J)
        FXX(J)  = EPSSQ * EPS(J) * A0XX(J)
        FXT(J)  = EPSSQ * (A0XT(J) + 3._r8*EPS1T*A0XX(J))
        FTT(J)  = EPS(J) * (A0TT(J) + EPS1T*A0XT(J)
     &                   + EPS1T*(EPS1T-EPS(J))*A0XX(J) )
      enddo
c---FLUX vector: FTT(J) is mass transferred [J+J-1] => [J+J], FXT, FXX moments
      do J = 1,NJ
        if (LGC_AB(J))  then
          Q0F(J+J) = Q0F(J+J) + FTT(J)
          Q1F(J+J) = Q1F(J+J) + FXT(J)
        else
          Q0F(J+J) = Q0F(J+J) - FTT(J)
          Q1F(J+J) = Q1F(J+J) - FXT(J)
        endif
      enddo


C---order of calculation below is important for A0XT which is overwritten!
      do J = 1,NJ
        EPS1T  = 1._r8 - EPS(J)
        EPS1S  = EPS1T * EPS1T
        A0YY(J) = A0YY(J) - FYY(J)
        A0ZZ(J) = A0ZZ(J) - FZZ(J)
        A0YZ(J) = A0YZ(J) - FYZ(J)
        A0YT(J) = A0YT(J) - FYT(J)
        A0ZT(J) = A0ZT(J) - FZT(J)
        A0XY(J) = A0XY(J) * EPS1S
        A0XZ(J) = A0XZ(J) * EPS1S
        A0XT(J) = (A0XT(J) - 3._r8*EPS(J)*A0XX(J)) * EPS1S
        A0XX(J) = A0XX(J) * EPS1S * EPS1T
        A0TT(J) = A0TT(J) - FTT(J)
        A0M (J) = A0M(J) - UAB(J)

        B0M (J) = B0M(J) + UAB(J)
        EPS(J)  = UAB(J) / B0M(J)
      enddo

C combine transferred moments (F_s) into box B0 - ORDER important !
      do J = 1,NJ
        EPS1T   = 1._r8 - EPS(J)
        EPS1S   = EPS1T * EPS1T
        EPSSQ   = EPS(J) * EPS(J)
        B0ZZ(J) = B0ZZ(J) + FZZ(J)
        B0YY(J) = B0YY(J) + FYY(J)
        B0YZ(J) = B0YZ(J) + FYZ(J)
        B0XY(J) = EPS(J)*FXY(J) + EPS1T*B0XY(J)
     &            + 3._r8*(EPS(J)*B0YT(J) - EPS1T*FYT(J) )
        B0XZ(J) = EPS(J)*FXZ(J) + EPS1T*B0XZ(J)
     &            + 3._r8*(EPS(J)*B0ZT(J) - EPS1T*FZT(J))
        B0YT(J) = B0YT(J) + FYT(J)
        B0ZT(J) = B0ZT(J) + FZT(J)
        B0XX(J) = FXX(J)*EPSSQ + B0XX(J)*EPS1S
     &            + 5._r8*( EPS(J)*EPS1T*(B0XT(J) - FXT(J))
     &            + (EPS1T-EPS(J))*(EPS1T*FTT(J)-EPS(J)*B0TT(J)) )
        B0XT(J) = EPS(J)*FXT(J) + EPS1T*B0XT(J)
     &            + 3._r8*(EPS(J)*B0TT(J) - EPS1T*FTT(J))
        B0TT(J) = B0TT(J) + FTT(J)
      enddo

      if (LIMTR.ge.3) then
C---LIMTR = 3 & 4:  force daughter box B to be monotonic and bounded by min/max
        do J = 1,NJ
          LGC_L2  = B0XX(J) .gt. 0._r8
          LGC_LX  = B0XT(J) .gt. 0._r8
            F0      = B0TT(J)
         if (LGC_L2) then
c---B0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(B0XX(J), abs(B0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(B0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---B0XX(J) < 0:  curved down at ends, allow if B0XT(J) < B0TT(J)
            F1 = min(+F0, max(-F0, B0XT(J)))
            F2 = max(B0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          B0XT(J) = F1
          B0XX(J) = F2
          B0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XY(J)))
          B0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XZ(J)))
c---limit daughter box B by min/max of parent A & B boxes
C---use near-zero C0000 to avoid zero divide
c---rescale for X-max
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABXMX(J)*B0M(J)-B0TT(J))
     &                       /max(abs(B0XT(J))+B0XX(J), C000)))
          B0XT(J) = ABSCAL*B0XT(J)
          B0XX(J) = ABSCAL*B0XX(J)
c---rescale for X-min
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABXMN(J)*B0M(J))
     &                       /max(abs(B0XT(J))-B0XX(J), C000)))
          B0XT(J) = ABSCAL*B0XT(J)
          B0XX(J) = ABSCAL*B0XX(J)
c---rescale for XY moment
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABYMX(J)*B0M(J)-B0TT(J))
     &             /max(abs(B0XT(J))+B0XX(J)+abs(B0XY(J)), C000)))
          B0XY(J) = ABSCAL*B0XY(J)
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABYMN(J)*B0M(J))
     &             /max(abs(B0XT(J))-B0XX(J)+abs(B0XY(J)), C000)))
          B0XY(J) = ABSCAL*B0XY(J)
c---rescale for XZ moment
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABZMX(J)*B0M(J)-B0TT(J))
     &             /max(abs(B0XT(J))+B0XX(J)+abs(B0XZ(J)), C000)))
          B0XZ(J) = ABSCAL*B0XZ(J)
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABZMN(J)*B0M(J))
     &             /max(abs(B0XT(J))-B0XX(J)+abs(B0XZ(J)), C000)))
          B0XZ(J) = ABSCAL*B0XZ(J)
        enddo

      endif

C---Now re-mask back onto A1-- AND B1-- arrays
      do J = 1,NJ
        if (LGC_AB(J))  then
          A1M(J)  = A0M(J)
          A1TT(J) = A0TT(J)
          A1XX(J) = A0XX(J)
          A1YT(J) = A0YT(J)
          A1YY(J) = A0YY(J)
          A1ZT(J) = A0ZT(J)
          A1ZZ(J) = A0ZZ(J)
          A1YZ(J) = A0YZ(J)
          A1XT(J) = A0XT(J)
          A1XY(J) = A0XY(J)
          A1XZ(J) = A0XZ(J)
C---
          B1M(J)  = B0M(J)
          B1TT(J) = B0TT(J)
          B1XX(J) = B0XX(J)
          B1YT(J) = B0YT(J)
          B1YY(J) = B0YY(J)
          B1ZT(J) = B0ZT(J)
          B1ZZ(J) = B0ZZ(J)
          B1YZ(J) = B0YZ(J)
          B1XT(J) = B0XT(J)
          B1XY(J) = B0XY(J)
          B1XZ(J) = B0XZ(J)
        else
          A1M(J)  = B0M(J)
          A1TT(J) = B0TT(J)
          A1XX(J) = B0XX(J)
          A1YT(J) = B0YT(J)
          A1YY(J) = B0YY(J)
          A1ZT(J) = B0ZT(J)
          A1ZZ(J) = B0ZZ(J)
          A1YZ(J) = B0YZ(J)
          A1XT(J) =-B0XT(J)
          A1XY(J) =-B0XY(J)
          A1XZ(J) =-B0XZ(J)
C---
          B1M(J)  = A0M(J)
          B1TT(J) = A0TT(J)
          B1XX(J) = A0XX(J)
          B1YT(J) = A0YT(J)
          B1YY(J) = A0YY(J)
          B1ZT(J) = A0ZT(J)
          B1ZZ(J) = A0ZZ(J)
          B1YZ(J) = A0YZ(J)
          B1XT(J) =-A0XT(J)
          B1XY(J) =-A0XY(J)
          B1XZ(J) =-A0XZ(J)
        endif
c       if (A1M(J).lt.0._r8 .or. B1M(J).lt.0._r8)  then
c         write(6,'(A,4I5,1P4E15.7)')  'qvect3',LIMTR,MSTEP,M,J
c    +           ,A1M(J),B1M(J),UAB(J)
c         do I = 1,NQ
c           write(6,'(I4,1P3E15.7)') I,QM(I),QU(I),QU0(I)
c         enddo
c         write(6,*) QMMIN
c         call ctmEXITC('*** negative air mass in QVECT3 ***')
c       endif
      enddo




C---#2 NOW we need to repeat the whole process with ODD pairs:
C---   this wraps box(NQ) to box(1), need to have QU(NQ+1) defined

      QU(NQ+1)  = QU(1)

      do J = 1,NJ
        LGC_AB(J) = QU(J+J+1) .ge. 0._r8
        UAB(J)    = abs(QU(J+J+1))
      enddo

C---Now shift A,B down by one, B(1)=>A(1), A(2)=>B(1),...B(N)=>A(N), A(1)=>B(N)
C---Now create vector set of A0 & B0 quantities reflecting sign of QU
C---  NB first-order X-moments switch sign for reverse advection
      do J = 1,NJ-1
        JP1  = J + 1
        if (LGC_AB(J))  then
          A0M(J)  = B1M(J)
          A0TT(J) = B1TT(J)
          A0XX(J) = B1XX(J)
          A0YT(J) = B1YT(J)
          A0YY(J) = B1YY(J)
          A0ZT(J) = B1ZT(J)
          A0ZZ(J) = B1ZZ(J)
          A0YZ(J) = B1YZ(J)
          A0XT(J) = B1XT(J)
          A0XY(J) = B1XY(J)
          A0XZ(J) = B1XZ(J)
C---
          B0M(J)  = A1M(JP1)
          B0TT(J) = A1TT(JP1)
          B0XX(J) = A1XX(JP1)
          B0YT(J) = A1YT(JP1)
          B0YY(J) = A1YY(JP1)
          B0ZT(J) = A1ZT(JP1)
          B0ZZ(J) = A1ZZ(JP1)
          B0YZ(J) = A1YZ(JP1)
          B0XT(J) = A1XT(JP1)
          B0XY(J) = A1XY(JP1)
          B0XZ(J) = A1XZ(JP1)
         else
          A0M(J)  = A1M(JP1)
          A0TT(J) = A1TT(JP1)
          A0XX(J) = A1XX(JP1)
          A0YT(J) = A1YT(JP1)
          A0YY(J) = A1YY(JP1)
          A0ZT(J) = A1ZT(JP1)
          A0ZZ(J) = A1ZZ(JP1)
          A0YZ(J) = A1YZ(JP1)
          A0XT(J) =-A1XT(JP1)
          A0XY(J) =-A1XY(JP1)
          A0XZ(J) =-A1XZ(JP1)
C---
          B0M(J)  = B1M(J)
          B0TT(J) = B1TT(J)
          B0XX(J) = B1XX(J)
          B0YT(J) = B1YT(J)
          B0YY(J) = B1YY(J)
          B0ZT(J) = B1ZT(J)
          B0ZZ(J) = B1ZZ(J)
          B0YZ(J) = B1YZ(J)
          B0XT(J) =-B1XT(J)
          B0XY(J) =-B1XY(J)
          B0XZ(J) =-B1XZ(J)
        endif
      enddo
      if (LGC_AB(NJ))  then
        A0M(NJ)  = B1M(NJ)
        A0TT(NJ) = B1TT(NJ)
        A0XX(NJ) = B1XX(NJ)
        A0YT(NJ) = B1YT(NJ)
        A0YY(NJ) = B1YY(NJ)
        A0ZT(NJ) = B1ZT(NJ)
        A0ZZ(NJ) = B1ZZ(NJ)
        A0YZ(NJ) = B1YZ(NJ)
        A0XT(NJ) = B1XT(NJ)
        A0XY(NJ) = B1XY(NJ)
        A0XZ(NJ) = B1XZ(NJ)
C---
        B0M(NJ)  = A1M(1)
        B0TT(NJ) = A1TT(1)
        B0XX(NJ) = A1XX(1)
        B0YT(NJ) = A1YT(1)
        B0YY(NJ) = A1YY(1)
        B0ZT(NJ) = A1ZT(1)
        B0ZZ(NJ) = A1ZZ(1)
        B0YZ(NJ) = A1YZ(1)
        B0XT(NJ) = A1XT(1)
        B0XY(NJ) = A1XY(1)
        B0XZ(NJ) = A1XZ(1)
       else
        A0M(NJ)  = A1M(1)
        A0TT(NJ) = A1TT(1)
        A0XX(NJ) = A1XX(1)
        A0YT(NJ) = A1YT(1)
        A0YY(NJ) = A1YY(1)
        A0ZT(NJ) = A1ZT(1)
        A0ZZ(NJ) = A1ZZ(1)
        A0YZ(NJ) = A1YZ(1)
        A0XT(NJ) =-A1XT(1)
        A0XY(NJ) =-A1XY(1)
        A0XZ(NJ) =-A1XZ(1)
C---
        B0M(NJ)  = B1M(NJ)
        B0TT(NJ) = B1TT(NJ)
        B0XX(NJ) = B1XX(NJ)
        B0YT(NJ) = B1YT(NJ)
        B0YY(NJ) = B1YY(NJ)
        B0ZT(NJ) = B1ZT(NJ)
        B0ZZ(NJ) = B1ZZ(NJ)
        B0YZ(NJ) = B1YZ(NJ)
        B0XT(NJ) =-B1XT(NJ)
        B0XY(NJ) =-B1XY(NJ)
        B0XZ(NJ) =-B1XZ(NJ)
      endif


C---LIMTR = 1:  positive definite, reset Sxx to allow maximum Sx (1.5*So)
      if (LIMTR.eq.1) then
        do J = 1,NJ
            F0 = A0TT(J)
            F1 = min(+C149*F0, max(-C149*F0, A0XT(J)))
            F2 = min((F0+F0-C034*abs(F1)), max(abs(F1)-F0, A0XX(J)))
          A0XT(J) = F1
          A0XX(J) = F2
          A0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XY(J)))
          A0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XZ(J)))
        enddo
      endif

C---LIMTR = 2:  positive and monotonic in box (<X>-dim only)
c---LIMTR = 3:  min-max, needs to have positive-monotonic also
      if (LIMTR.ge.2) then
        do J = 1,NJ
          LGC_L2  = A0XX(J) .gt. 0._r8
          LGC_LX  = A0XT(J) .gt. 0._r8
            F0      = A0TT(J)
         if (LGC_L2) then
c---A0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(A0XX(J), abs(A0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(A0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---A0XX(J) < 0:  curved down at ends, allow if A0XT(J) < A0TT(J)
            F1 = min(+F0, max(-F0, A0XT(J)))
            F2 = max(A0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          A0XT(J) = F1
          A0XX(J) = F2
          A0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XY(J)))
          A0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, A0XZ(J)))
C---no need to impose limits on downwind B since it is not transported
        enddo
      endif


c---LIMTR = 3:  min-max, force positive-monotonic on B boxes also
c---             find and store min/max of both A & B boxes
      if (LIMTR.ge.3) then
        do J = 1,NJ
          LGC_L2  = B0XX(J) .gt. 0._r8
          LGC_LX  = B0XT(J) .gt. 0._r8
            F0      = B0TT(J)
         if (LGC_L2) then
c---B0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(B0XX(J), abs(B0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(B0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---B0XX(J) < 0:  curved down at ends, allow if B0XT(J) < B0TT(J)
            F1 = min(+F0, max(-F0, B0XT(J)))
            F2 = max(B0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          B0XT(J) = F1
          B0XX(J) = F2
          B0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XY(J)))
          B0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XZ(J)))
        enddo

C---LIMTR = 3 min/max in A-- & B-- boxes occurs at end points if monotonic
       do J = 1,NJ
       ABXMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J))/B0M(J)  )
       ABXMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J))/B0M(J)  )
       ABYMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J)+abs(A0XY(J)))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J)+abs(B0XY(J)))/B0M(J)  )
       ABYMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J)-abs(A0XY(J)))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J)-abs(B0XY(J)))/B0M(J)  )
       ABZMX(J)= max((A0TT(J)+abs(A0XT(J))+A0XX(J)+abs(A0XZ(J)))/A0M(J),
     &             (B0TT(J)+abs(B0XT(J))+B0XX(J)+abs(B0XZ(J)))/B0M(J)  )
       ABZMN(J)= min((A0TT(J)-abs(A0XT(J))+A0XX(J)-abs(A0XZ(J)))/A0M(J),
     &             (B0TT(J)-abs(B0XT(J))+B0XX(J)-abs(B0XZ(J)))/B0M(J)  )
       enddo
      endif

C---Advect from A0 ==[UAB]==> to B0  for even pairs [1:2] [3:4] [5:6]
      do J = 1,NJ
        EPS(J) = UAB(J) / A0M(J)
        EPS1T  = 1._r8 - EPS(J)
        EPS1S  = EPS1T * EPS1T
        EPSSQ  = EPS(J) * EPS(J)
c------Seperate left A0 box to 2 parts: F_s = section tranferred to B0 box
        FYY(J)  = EPS(J) * A0YY(J)
        FZZ(J)  = EPS(J) * A0ZZ(J)
        FYZ(J)  = EPS(J) * A0YZ(J)
        FYT(J)  = EPS(J) * (A0YT(J) + EPS1T * A0XY(J))
        FZT(J)  = EPS(J) * (A0ZT(J) + EPS1T * A0XZ(J))
        FXY(J)  = EPSSQ * A0XY(J)
        FXZ(J)  = EPSSQ * A0XZ(J)
        FXX(J)  = EPSSQ * EPS(J) * A0XX(J)
        FXT(J)  = EPSSQ * (A0XT(J) + 3._r8*EPS1T*A0XX(J))
        FTT(J)  = EPS(J) * (A0TT(J) + EPS1T*A0XT(J)
     &                   + EPS1T*(EPS1T-EPS(J))*A0XX(J) )
      enddo
c---FLUX vector: mass transferred A<=>B i.e., [J+J] <=> [J+J+1]
      do J = 1,NJ
        if (LGC_AB(J))  then
          Q0F(J+J+1) = Q0F(J+J+1) + FTT(J)
          Q1F(J+J+1) = Q1F(J+J+1) + FXT(J)
        else
          Q0F(J+J+1) = Q0F(J+J+1) - FTT(J)
          Q1F(J+J+1) = Q1F(J+J+1) - FXT(J)
        endif
      enddo

C---order is important for A0XT!
      do J = 1,NJ
        EPS1T  = 1._r8 - EPS(J)
        EPS1S  = EPS1T * EPS1T
        A0YY(J) = A0YY(J) - FYY(J)
        A0ZZ(J) = A0ZZ(J) - FZZ(J)
        A0YZ(J) = A0YZ(J) - FYZ(J)
        A0YT(J) = A0YT(J) - FYT(J)
        A0ZT(J) = A0ZT(J) - FZT(J)
        A0XY(J) = A0XY(J) * EPS1S
        A0XZ(J) = A0XZ(J) * EPS1S
        A0XT(J) = (A0XT(J) - 3._r8*EPS(J)*A0XX(J)) * EPS1S
        A0XX(J) = A0XX(J) * EPS1S * EPS1T
        A0TT(J) = A0TT(J) - FTT(J)
        A0M (J) = A0M(J) - UAB(J)

        B0M (J) = B0M(J) + UAB(J)
        EPS(J)  = UAB(J) / B0M(J)
      enddo

C combine transferred moments (F_s) into box B0 - ORDER important !
      do J = 1,NJ
        EPS1T   = 1._r8 - EPS(J)
        EPS1S   = EPS1T * EPS1T
        EPSSQ   = EPS(J) * EPS(J)
        B0ZZ(J) = B0ZZ(J) + FZZ(J)
        B0YY(J) = B0YY(J) + FYY(J)
        B0YZ(J) = B0YZ(J) + FYZ(J)
        B0XY(J) = EPS(J)*FXY(J) + EPS1T*B0XY(J)
     &            + 3._r8*(EPS(J)*B0YT(J) - EPS1T*FYT(J) )
        B0XZ(J) = EPS(J)*FXZ(J) + EPS1T*B0XZ(J)
     &            + 3._r8*(EPS(J)*B0ZT(J) - EPS1T*FZT(J))
        B0YT(J) = B0YT(J) + FYT(J)
        B0ZT(J) = B0ZT(J) + FZT(J)
        B0XX(J) = FXX(J)*EPSSQ + B0XX(J)*EPS1S
     &            + 5._r8*( EPS(J)*EPS1T*(B0XT(J) - FXT(J))
     &            + (EPS1T-EPS(J))*(EPS1T*FTT(J)-EPS(J)*B0TT(J)) )
        B0XT(J) = EPS(J)*FXT(J) + EPS1T*B0XT(J)
     &            + 3._r8*(EPS(J)*B0TT(J) - EPS1T*FTT(J))
        B0TT(J) = B0TT(J) + FTT(J)
      enddo


      if (LIMTR.ge.3) then
C---LIMTR = 3:  force daughter box B to be monotonic and bounded by min/max
        do J = 1,NJ
          LGC_L2  = B0XX(J) .gt. 0._r8
          LGC_LX  = B0XT(J) .gt. 0._r8
            F0      = B0TT(J)
         if (LGC_L2) then
c---B0XX(J) > 0:  do not allow reversal of curvature
            F2 = min(B0XX(J), abs(B0XT(J))*C033, 0.5_r8*F0)
            F1 = min(F0+F2, abs(B0XT(J)))
           if (.not. LGC_LX) then
            F1 = -F1
           endif
         else
c---B0XX(J) < 0:  curved down at ends, allow if B0XT(J) < B0TT(J)
            F1 = min(+F0, max(-F0, B0XT(J)))
            F2 = max(B0XX(J), -F0+abs(F1), -abs(F1)*C033)
         endif
          B0XT(J) = F1
          B0XX(J) = F2
          B0XY(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XY(J)))
          B0XZ(J) = min(+CXYZ*F0, max(-CXYZ*F0, B0XZ(J)))
c---limit daughter box B by min/max of parent A & B boxes
C---use near-zero C0000 to avoid zero divide
c---rescale for X-max
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABXMX(J)*B0M(J)-B0TT(J))
     &                       /max(abs(B0XT(J))+B0XX(J), C000)))
          B0XT(J) = ABSCAL*B0XT(J)
          B0XX(J) = ABSCAL*B0XX(J)
c---rescale for X-min
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABXMN(J)*B0M(J))
     &                       /max(abs(B0XT(J))-B0XX(J), C000)))
          B0XT(J) = ABSCAL*B0XT(J)
          B0XX(J) = ABSCAL*B0XX(J)
c---rescale for XY moment
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABYMX(J)*B0M(J)-B0TT(J))
     &             /max(abs(B0XT(J))+B0XX(J)+abs(B0XY(J)), C000)))
          B0XY(J) = ABSCAL*B0XY(J)
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABYMN(J)*B0M(J))
     &             /max(abs(B0XT(J))-B0XX(J)+abs(B0XY(J)), C000)))
          B0XY(J) = ABSCAL*B0XY(J)
c---rescale for XZ moment
        ABSCAL = max(0._r8, min(1._r8, (C3MAX*ABZMX(J)*B0M(J)-B0TT(J))
     &             /max(abs(B0XT(J))+B0XX(J)+abs(B0XZ(J)), C000)))
          B0XZ(J) = ABSCAL*B0XZ(J)
        ABSCAL = max(0._r8, min(1._r8, (B0TT(J)-C3MIN*ABZMN(J)*B0M(J))
     &             /max(abs(B0XT(J))-B0XX(J)+abs(B0XZ(J)), C000)))
          B0XZ(J) = ABSCAL*B0XZ(J)
         enddo

      endif

C---  Now put the A--_s and B--_s back to the pair of Q--_s
      do J = 1,NJ
        JJ    = J + J
        JJP1  = J + J + 1
        if (LGC_AB(J))  then
          QM(JJ)  = A0M(J)
          QTT(JJ) = A0TT(J)
          QXX(JJ) = A0XX(J)
          QYT(JJ) = A0YT(J)
          QYY(JJ) = A0YY(J)
          QZT(JJ) = A0ZT(J)
          QZZ(JJ) = A0ZZ(J)
          QYZ(JJ) = A0YZ(J)
          QXT(JJ) = A0XT(J)
          QXY(JJ) = A0XY(J)
          QXZ(JJ) = A0XZ(J)
C---  same for the B_s
          QM(JJP1)  = B0M(J)
          QTT(JJP1) = B0TT(J)
          QXX(JJP1) = B0XX(J)
          QYT(JJP1) = B0YT(J)
          QYY(JJP1) = B0YY(J)
          QZT(JJP1) = B0ZT(J)
          QZZ(JJP1) = B0ZZ(J)
          QYZ(JJP1) = B0YZ(J)
          QXT(JJP1) = B0XT(J)
          QXY(JJP1) = B0XY(J)
          QXZ(JJP1) = B0XZ(J)
         else
          QM(JJP1)  = A0M(J)
          QTT(JJP1) = A0TT(J)
          QXX(JJP1) = A0XX(J)
          QYT(JJP1) = A0YT(J)
          QYY(JJP1) = A0YY(J)
          QZT(JJP1) = A0ZT(J)
          QZZ(JJP1) = A0ZZ(J)
          QYZ(JJP1) = A0YZ(J)
          QXT(JJP1) =-A0XT(J)
          QXY(JJP1) =-A0XY(J)
          QXZ(JJP1) =-A0XZ(J)
C---  same for the B_s
          QM(JJ)  = B0M(J)
          QTT(JJ) = B0TT(J)
          QXX(JJ) = B0XX(J)
          QYT(JJ) = B0YT(J)
          QYY(JJ) = B0YY(J)
          QZT(JJ) = B0ZT(J)
          QZZ(JJ) = B0ZZ(J)
          QYZ(JJ) = B0YZ(J)
          QXT(JJ) =-B0XT(J)
          QXY(JJ) =-B0XY(J)
          QXZ(JJ) =-B0XZ(J)
        endif
      enddo
      QM(1)  = QM(1+NQ)
      QTT(1) = QTT(1+NQ)
      QXT(1) = QXT(1+NQ)
      QXX(1) = QXX(1+NQ)
      QYT(1) = QYT(1+NQ)
      QYY(1) = QYY(1+NQ)
      QZT(1) = QZT(1+NQ)
      QZZ(1) = QZZ(1+NQ)
      QXY(1) = QXY(1+NQ)
      QXZ(1) = QXZ(1+NQ)
      QYZ(1) = QYZ(1+NQ)

      Q0F(1) = Q0F(1+NQ)
      Q1F(1) = Q1F(1+NQ)

      enddo
c<<<<<<end primary multi-stepping loop

c---The flux and moments of tracer -- Q0F(I), Q1F(I) -- are carried with
c---    air mass QU(I) from box[I] to box [I-1].
c---Some applications at UCI require the flux for air with tracer < F0 (kg/kg)
c---    With multi-stepping, we can compute the 1st-moment assuming it is the
c---    same for all of the steps  (this is done below).
c---    1) calculate the fraction of QU with tracer abundance < F0 (kg/kg)
c---         XF0 = [ F0*QU - Q0F + abs(Q1F) ] / [ 2 * abs(Q1F) ]
c---            make sure that Q1F is at least +-1.d-10*Q0F for division
c---            place limits on XF0 of [0, 1]
c---    2) calculate total trace mass with abundance < F0
c---         QF0 = [ F0*QU + Q0F - abs(Q1F) ]/2 * XF0
c---             place limits on QF0 of [0, Q0F]
c---          (or put upper lim on F0*QU of Q0F+abs(Q1F))

c---adjust the 'flux' boxes Q0/Q1/Q2 if multi-step, just assume a single slope:
c---     need to increase first moment, second moment not reliable.
        do J=1,NQ+1
c         Q1F(J) = min(Q0F(J), max(-Q0F(J), real(MSTEP)*Q1F(J))) !wrong
          Q1F(J)  =SIGN(MAX(ABS(Q0F(J)),
     &          real(MSTEP,r8)*ABS(Q1F(J))),Q1F(J))
        enddo

      return
      end

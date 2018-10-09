c-----------------------------------------------------------------------
c---(p-dyn2.f)---(p-code 6.1b 2/2012): new over-the-pole flow version 4
c
c         U-V tracer advection parallel over layers L=1:LM
c         W   tracer adv parallel over IJ-blocks
c         W   tracer adv also combines large-scale W & convective W
c
c---subroutines:  DYN2UL, DYN2VL, DYN2W_OC, QLIMIT2, POLES1, POLES2
c      ***have dropped non-L, parallel-over-tracers versions**


c-----------------------------------------------------------------------
      subroutine DYN2UL(DTALFA,AIRUV,ALFA2D,L,QFU,USTEP)
c-----------------------------------------------------------------------
c       WEST-to-EAST ADVECTION OF TRACE COMPOUNDS using S.O.M.
c-----------------------------------------------------------------------
c  new CTM q-code 6.0b (5/2011) - OMP parallelized over LM(layer)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, NPAR
      use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW,
     &     NTM, LMTSOM
      implicit none
C-----------------------------------------------------------------------

      real(r8),  intent(in) :: DTALFA
      real(r8),  intent(in) :: ALFA2D(IPAR+1,JPAR)
      integer, intent(in) :: L                     ! Layer index
      real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
      real(r8),  intent(out) :: QFU(IPAR+1,JPAR,NPAR,2)
      integer, intent(out) :: USTEP                ! sum of NSTEP

      real(r8), dimension(IPAR+1) ::
     &     QM,            ! air mass in box at start,
     &     QU,            ! air mass flux moved [I-1]->[I] in adv. step
     &     QTT,           ! tracer mass in box [I]
     &     QXT,QYT,QZT,   ! 1st moments of tracer in U, V, W direction
     &     QXX,QYY,QZZ,   ! 2nd moments of tracer in U, V, W direction
     &     QXY,QYZ,QXZ,   ! cross-moments of tracer
     &     Q0F,Q1F        ! computed tracer flux from [I-1] to [I]
      real(r8)   AIRUVIJ(IPAR,JPAR)
      integer :: NQ       ! length of vector pipe for advection (assumed cyclic)
      integer :: NSTEP    ! #multi-steps needed for local CFL, ret by QVECT3
      integer :: I,J,N
      integer :: IM, JM   !Locals for CTM3
c-----------------------------------------------------------------------
c>>>check on negative airmass for debug ONLY
c     do J = 1,JM
c       do I = 1,IM
c         if (AIRUV(I,J) .lt. 0._r8)  then
c           call EXITIJL (' AIRTR < 0 in DYN2UL', I,J,L)
c         endif
c       enddo
c     enddo

      USTEP = 0

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR

      do N = 1,NTM
C---V 6.0b has over-the-pole flow, combining J=1:2 and JM-1:JM
        do J = 2,JM-1

c---transfer STT(), SXT(), ... into local vector for piped flow
          do I = 1,IM
            QM(I)  = AIRUV(I,J)
            QU(I)  = ALFA2D(I,J) * DTALFA
            QTT(I) = STT(I,J,L,N)
            QXT(I) = SUT(I,J,L,N)
            QYT(I) = SVT(I,J,L,N)
            QZT(I) = SWT(I,J,L,N)
            QXX(I) = SUU(I,J,L,N)
            QYY(I) = SVV(I,J,L,N)
            QZZ(I) = SWW(I,J,L,N)
            QXY(I) = SUV(I,J,L,N)
            QXZ(I) = SUW(I,J,L,N)
            QYZ(I) = SVW(I,J,L,N)
          enddo
            QU(IM+1) = QU(1)
            NQ = IM
c         write(6,'(A,2I4)') 'call QVECT3 in DYN2UL N/J',N,J
c-----------------------------------------------------------------------
          call QVECT3(QM,QU,NQ,LMTSOM, NSTEP,
     &           QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ, Q0F,Q1F)
          USTEP = USTEP + NSTEP
c-----------------------------------------------------------------------
c---transfer Q__() back to STT() & moments
          do I = 1,IM
            STT(I,J,L,N) = QTT(I)
            SUT(I,J,L,N) = QXT(I)
            SVT(I,J,L,N) = QYT(I)
            SWT(I,J,L,N) = QZT(I)
            SUU(I,J,L,N) = QXX(I)
            SVV(I,J,L,N) = QYY(I)
            SWW(I,J,L,N) = QZZ(I)
            SUV(I,J,L,N) = QXY(I)
            SUW(I,J,L,N) = QXZ(I)
            SVW(I,J,L,N) = QYZ(I)
            AIRUVIJ(I,J) = QM(I)
            QFU(I,J,N,1) = Q0F(I)
            QFU(I,J,N,2) = Q1F(I)
          enddo
            QFU(IM+1,J,N,1) = Q0F(IM+1)
            QFU(IM+1,J,N,2) = Q1F(IM+1)

        enddo               ! J loop
      enddo                 ! N tracer loop

c--reset air mass in layer L
        do J = 2,JM-1
          do I = 1,IM
            AIRUV(I,J) = AIRUVIJ(I,J)
          enddo
        enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine DYN2VL(DTBETA,AIRUV,BETA2D,L,QFV,VSTEP)
c-----------------------------------------------------------------------
c       SOUTH-to-NORTH ADVECTION OF TRACE COMPOUNDS using S.O.M.
c-----------------------------------------------------------------------
c  new CTM p-code 6.0b (5/2011) - OMP parallelized over L(layer)
c  over-the-pole flow changed:
c        boxes 1&2 JM&JM-1 combined, all ALFAs used, BETA(J=1 & JM+1)=0
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, NPAR
      use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW,
     &     NTM, LMTSOM
      implicit none
C-----------------------------------------------------------------------

      real(r8),  intent(in) :: DTBETA
      real(r8),  intent(in) :: BETA2D(IPAR,JPAR+1)
      integer, intent(in) :: L                      ! Layer index
      real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
      real(r8),  intent(out) :: QFV(IPAR,JPAR+1,NPAR,2)
      integer, intent(out) :: VSTEP                ! sum of NSTEP

      real(r8), dimension(2*JPAR+1) ::
     &     QM,            ! air mass in box at start,
     &     QU,            ! air mass flux moved [I-1]->[I] in adv. step
     &     QTT,           ! tracer mass in box [I]
     &     QXT,QYT,QZT,   ! 1st moments of tracer in U, V, W direction
     &     QXX,QYY,QZZ,   ! 2nd moments of tracer in U, V, W direction
     &     QXY,QYZ,QXZ,   ! cross-moments of tracer
     &     Q0F,Q1F        ! computed tracer flux from [I-1] to [I]
      real(r8) :: AIRUVIJ(IPAR,JPAR)
      real(r8) :: STTJ1(2),STTJM(2)
      integer :: NQ       ! length of vector pipe for advection (assumed cyclic)
      integer :: NSTEP    ! #multi-steps needed for local CFL, ret by QVECT3
      integer :: I,J,N,II,JJ,IQD2
      integer :: IM, JM   !Locals for CTM3
c-----------------------------------------------------------------------
c---check on negative airmass for debug ONLY
c     do J = 1,JM
c       do I = 1,IM
c         if (AIRUV(I,J) .lt. 0._r8)  then
c           call EXITIJL (' AIRTR < 0 in DYN2VL', I,J,L)
c         endif
c       enddo
c     enddo

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR

c---N-S flow is periodic, over-the-pole,  combines opposite meridions:
      IQD2  = IM/2
      NQ    = JM+JM-4
      VSTEP = 0

      do N = 1,NTM
        do I = 1,IQD2
           II = I + IQD2

c-----------------------------------------------------------------------
c---start at SPole and follow up merdion stripe at I, using boxes 2:JM-1
          do J = 1,JM-2
              JJ = J+1
            QM(J)  = AIRUV(I,JJ)
            QU(J)  = BETA2D(I,JJ)*DTBETA
            QTT(J) = STT(I,JJ,L,N)
            QXT(J) = SUT(I,JJ,L,N)
            QYT(J) = SVT(I,JJ,L,N)
            QZT(J) = SWT(I,JJ,L,N)
            QXX(J) = SUU(I,JJ,L,N)
            QYY(J) = SVV(I,JJ,L,N)
            QZZ(J) = SWW(I,JJ,L,N)
            QXY(J) = SUV(I,JJ,L,N)
            QXZ(J) = SUW(I,JJ,L,N)
            QYZ(J) = SVW(I,JJ,L,N)
          enddo
c---continue down (sign reversal on U & V moments, QU) merdion stripe II
c---                        using boxes JM-1:2
          do J = JM-1,JM+JM-4
              JJ = JM+JM-2 - J
            QM(J) = AIRUV(II,JJ)
            QU(J)  = -BETA2D(II,JJ+1)*DTBETA
            QTT(J) =  STT(II,JJ,L,N)
            QXT(J) = -SUT(II,JJ,L,N)
            QYT(J) = -SVT(II,JJ,L,N)
            QZT(J) =  SWT(II,JJ,L,N)
            QXX(J) =  SUU(II,JJ,L,N)
            QYY(J) =  SVV(II,JJ,L,N)
            QZZ(J) =  SWW(II,JJ,L,N)
            QXY(J) = +SUV(II,JJ,L,N)
            QXZ(J) = -SUW(II,JJ,L,N)
            QYZ(J) = -SVW(II,JJ,L,N)
          enddo

            QU(NQ+1) = QU(1)

c-----------------------------------------------------------------------
          call QVECT3(QM,QU,NQ,LMTSOM, NSTEP,
     &           QTT,QYT,QYY,QYZ,QXY,QZT,QZZ,QXT,QXX,QXZ, Q0F,Q1F)
          VSTEP = VSTEP + NSTEP
c-----------------------------------------------------------------------

c---replace meridion stripe up I
          do J = 1,JM-2
             JJ = J+1
            STT(I,JJ,L,N) = QTT(J)
            SUT(I,JJ,L,N) = QXT(J)
            SVT(I,JJ,L,N) = QYT(J)
            SWT(I,JJ,L,N) = QZT(J)
            SUU(I,JJ,L,N) = QXX(J)
            SVV(I,JJ,L,N) = QYY(J)
            SWW(I,JJ,L,N) = QZZ(J)
            SUV(I,JJ,L,N) = QXY(J)
            SUW(I,JJ,L,N) = QXZ(J)
            SVW(I,JJ,L,N) = QYZ(J)
            AIRUVIJ(I,JJ) = QM(J)
c---Q0F, Q1F, Q2F[J] are fluxes from [J-1]=>[J] of STT, SXT, SXX
c     with ver 6.0b they are ONLY defined for J=2,JM, not 1,JM+1
            QFV(I,JJ,N,1)  = Q0F(J)
            QFV(I,JJ,N,2)  = Q1F(J)
          enddo
            QFV(I,JM,N,1)  = Q0F(JM-1)
            QFV(I,JM,N,2)  = Q1F(JM-1)

c---replace meridion stripe down II (note reverse of Y moments)
          do J = JM-1,JM+JM-4
              JJ = JM+JM-2 - J
            STT(II,JJ,L,N) =  QTT(J)
            SUT(II,JJ,L,N) = -QXT(J)
            SVT(II,JJ,L,N) = -QYT(J)
            SWT(II,JJ,L,N) =  QZT(J)
            SUU(II,JJ,L,N) =  QXX(J)
            SVV(II,JJ,L,N) =  QYY(J)
            SWW(II,JJ,L,N) =  QZZ(J)
            SUV(II,JJ,L,N) = +QXY(J)
            SUW(II,JJ,L,N) = -QXZ(J)
            SVW(II,JJ,L,N) = -QYZ(J)
            AIRUVIJ(II,JJ) =  QM(J)
            QFV(II,JJ+1,N,1) = -Q0F(J)
            QFV(II,JJ+1,N,2) = -Q1F(J)
          enddo
            QFV(II,2,N,1)  = -Q0F(JM+JM-3)
            QFV(II,2,N,2)  = -Q1F(JM+JM-3)

        enddo               ! I loop
      enddo                 ! N tracer loop

c--reset air mass in layer L
        do J = 1,JM
          do I = 1,IM
            AIRUV(I,J) = AIRUVIJ(I,J)
          enddo
        enddo

      return
      end



c-----------------------------------------------------------------------
      subroutine DYN2W_OC(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ
     &                ,AIRB,GAMACB,GAMAB,DTGAMA,MP,WSTEP)
c-----------------------------------------------------------------------
c       VERTICAL(UPWARD) ADVECTION OF TRACE COMPOUNDS using S.O.M.
c-----------------------------------------------------------------------
c  new CTM p-code 5.3 (7/2007) - OMP parallelized over IJ-blocks
c
c  Amund Sovde, October 2008
c  Instead of several vertical pipes (I,J),(I+1,J),(I+2,J),... combined
c  into a long pipe, IMDIV components are stacked. In this way,
c  all components in the same long pipe need the same internal time
c  stepping, so that no column (pipe-part) forces other columns (parts of
c  the pipe) to be slower than necessary. Also, the striding is reduced.
c
c  I think this still allowes QVECT3 to be vectorized.
c  The vertical advection works on private IJ-block arrays BTT, BXT, ...
c  the large-scale vertical w's (GAMAB) and the convective mesoscale w
c     (GAMACB) are combined into a single advective step
c-----------------------------------------------------------------------
      use cmn_precision, only: r8, rMom
      use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, IMDIV
      use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE,
     &     NTM, LMTSOM
      implicit none
C-----------------------------------------------------------------------

      real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
      real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: 
     &     BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
      real(r8), intent(inout), dimension(LPAR,IDBLK,JDBLK) :: AIRB
      real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: GAMACB,GAMAB
      real(r8), intent(in) :: DTGAMA
      integer,intent(in) :: MP
      integer, intent(out) :: WSTEP                ! sum of NSTEP

      !// Use another name; NPARDIV
      integer, parameter :: NPARDIV=IMDIV
      real(r8), dimension(NPARDIV*LPAR+1) ::
     &     QM,            ! air mass in box at start,
     &     QU,            ! air mass flux moved [I-1]->[I] in adv. step
     &     QTT,QTT2,           ! tracer mass in box [I]
     &     QXT,QYT,QZT,   ! 1st moments of tracer in U, V, W direction
     &     QXX,QYY,QZZ,   ! 2nd moments of tracer in U, V, W direction
     &     QXY,QYZ,QXZ,   ! cross-moments of tracer
     &     Q0F,Q1F        ! computed tracer flux from [I-1] to [I]
      integer :: NQ       ! length of vector pipe for advection (assumed cyclic)
      integer :: NSTEP    ! #multi-steps needed for local CFL, ret by QVECT3
      integer :: I,J,L,II,JJ,NB,M,ML,N,NREST,NXTRA
      integer :: LM !Locals for CTM3
c-----------------------------------------------------------------------
      character(len=*), parameter :: subr = 'DYN2W_OC'
c-----------------------------------------------------------------------


      WSTEP = 0

c Set for local use (IM, JM, LM are removed from CTM3)
      LM = LPAR

      II   = MPBLKIE(MP)-MPBLKIB(MP)+1
      JJ   = MPBLKJE(MP)-MPBLKJB(MP)+1
c---combine IMDIV columns of W advection into one:  NQ = LM * IMDIV
c---if (mod(II,IMDIV).ne.0) then TROUBLE, but check once at startup

      !// Combines NPARDIV tracers in the column to a long pipe:
      !// NQ = LM*NPARDIV.

      NQ   = LM * NPARDIV

      !// If NPAR (transported components) is not a multiple of NPARDIV,
      !// we get a rest which must be treated together with dummy tracers.
      NREST = mod(NPAR,NPARDIV)
      !// To transport all components, we need NXTRA dummy tracers to fill
      !// up the long pipe:
      NXTRA = NPARDIV-NREST

c>>>this check on negative airmass SHOULD be turned off after debug mode
c     do J = 1,JJ
c       do I = 1,II
c         do L = 1,LM
c           if (AIRB(L,I,J) .lt. 0._r8)  then
c               call EXITIJL (' AIRB < 0 in DYN2W', I,J,L)
c           endif
c         enddo
c       enddo
c     enddo

c---begin major loop over J's and I's
      do J = 1,JJ
        do I = 1,II

c---combine the vertical pipes,
c---need to ensure that GAMA(1,I,J)=0 so that pipes do NOT connect
          do M = 1,NPARDIV
            !// All work on same I, i.e. the same air columns are stacked
            do L = 1,LM
              ML      = L + LM*(M-1)
              QU(ML)  = (GAMAB(L,I,J)+GAMACB(L,I,J))*DTGAMA
            enddo
          enddo
              QU(ML+1)  = 0._r8
c---loop over tracers N, need NXTRA dummy tracers; stack NPARDIV tracers
c---in the long pipe
          do NB = 1,NTM+NXTRA,NPARDIV
            do M = 1,NPARDIV
              !// The actual tracer number
              N  = NB + M - 1
              !// Fill the last pipe with the same tracer.
              if (N .gt. NTM) N = NTM
              do L = 1,LM
                ML      = L + LM*(M-1)
C---must re-init the airmass for each tracer
                QM(ML)  = AIRB(L,I,J)
                QTT(ML) = BTT(L,N,I,J)
                QXT(ML) = BXT(L,N,I,J)
                QYT(ML) = BYT(L,N,I,J)
                QZT(ML) = BZT(L,N,I,J)
                QXX(ML) = BXX(L,N,I,J)
                QYY(ML) = BYY(L,N,I,J)
                QZZ(ML) = BZZ(L,N,I,J)
                QXY(ML) = BXY(L,N,I,J)
                QXZ(ML) = BXZ(L,N,I,J)
                QYZ(ML) = BYZ(L,N,I,J)
              enddo
            enddo
            QTT2=QTT
C-----------------------------------------------------------------------
            call QVECT3(QM,QU,NQ,LMTSOM, NSTEP,
     &           QTT,QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY, Q0F,Q1F)

            WSTEP = WSTEP + NSTEP
C-----------------------------------------------------------------------

            do M = 1,NPARDIV
              N  = NB + M - 1
              !// Fill up except the NXTRA for the last last pipe
              if (N .le. NTM) then
                do L = 1,LM
                  ML = L + LM*(M-1)
                  BTT(L,N,I,J) = QTT(ML)
                  BXT(L,N,I,J) = QXT(ML)
                  BYT(L,N,I,J) = QYT(ML)
                  BZT(L,N,I,J) = QZT(ML)
                  BXX(L,N,I,J) = QXX(ML)
                  BYY(L,N,I,J) = QYY(ML)
                  BZZ(L,N,I,J) = QZZ(ML)
                  BXY(L,N,I,J) = QXY(ML)
                  BXZ(L,N,I,J) = QXZ(ML)
                  BYZ(L,N,I,J) = QYZ(ML)
                enddo
              endif
            enddo !// M-loop
          enddo  !// NB loop

c---reset/update air mass in boxes after last tracer N
          !// Same column of air is transported for all parts of the
          !// long pipe; pick the first part of it
          do L = 1,LM
            ML = L
            AIRB(L,I,J) = QM(ML)
          enddo

        enddo       !// I loop
      enddo       !// J loop

      return
      end


c-----------------------------------------------------------------------
      subroutine QLIMIT2 (ETT,EXT,EXX,EXY,EXZ)
c-----------------------------------------------------------------------
c   quick SOM limiter only for LIM=2 (pos+mono) in the X direction
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      real(r8), intent(in) ::    ETT
      real(r8), intent(inout) :: EXT,EXX,EXY,EXZ

      logical :: LGC_LX,LGC_L2
      real(r8) :: F0,F1,F2

      real(r8), parameter :: C033 = 0.333333333333_r8
      real(r8), parameter :: CXYZ = 1.0_r8     ! limit XY-moments, 1.0 or 1.5

      LGC_L2  = EXX .gt. 0._r8
      LGC_LX  = EXT .gt. 0._r8
         F0 = ETT
      if (LGC_L2) then
c---EXX > 0:  do not allow reversal of curvature
         F2 = min(EXX, abs(EXT)*C033, 0.5_r8*F0)
         F1 = min(F0+F2, abs(EXT))
        if (.not. LGC_LX) then
         F1 = -F1
        endif
      else
c---EXX < 0:  curved down at ends, allow if EXT < ETT
         F1 = min(+F0, max(-F0, EXT))
         F2 = max(EXX, -F0+abs(F1), -abs(F1)*C033)
      endif

      EXT = F1
      EXX = F2
      EXY = min(+CXYZ*F0, max(-CXYZ*F0, EXY))
      EXZ = min(+CXYZ*F0, max(-CXYZ*F0, EXZ))

      return
      end


c-----------------------------------------------------------------------
       subroutine POLES1(AIRUV,ALFA2D,BETA2D,L)
c-----------------------------------------------------------------------
c     Combine polar-pie box with next lower latitude (J=1,2  =JM-1,JM)
c           using S.O.M. (Prather, Sovde)
c-----------------------------------------------------------------------
c  new CTM p-code 6.0b (5/2011) - OMP parallelized over LM(layer)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR
      use cmn_ctm, only: STT, NTM,
     &     SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW, ALFA,BETA
      implicit none
C-----------------------------------------------------------------------

       real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
       integer, intent(in) :: L                     ! Layer index
       real(r8),  intent(out) :: ALFA2D(IPAR+1,JPAR),BETA2D(IPAR,JPAR+1)

       real(r8) :: FTT,FUT,FVT,FWT,FUU,FVV,FWW,FUV,FUW,FVW
       real(r8) ::
     &      FTT1,FUT1,FVT1,FWT1,FUU1,FVV1,FWW1,FUV1,FUW1,FVW1,Q1,EPS1
       real(r8) ::
     &      FTT2,FUT2,FVT2,FWT2,FUU2,FVV2,FWW2,FUV2,FUW2,FVW2,Q2,EPS2
       real(r8) :: F0,F1,F2
       integer :: I,J,JJ,N
       integer :: IM, JM        !Locals for CTM3
c-----------------------------------------------------------------------
c---for alternative version 3 of polar flow:
c---   combine two next-to-pole boxes into one, do U-V transport,
c---   split back in 2 boxes when complete

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR

      do N = 1,NTM

       do J = 1,JM-1,JM-2
       do I = 1,IM
         Q1 = AIRUV(I,J)
         FTT1 = STT(I,J,L,N)
         FUT1 = SUT(I,J,L,N)
         FVT1 = SVT(I,J,L,N)
         FWT1 = SWT(I,J,L,N)
         FUU1 = SUU(I,J,L,N)
         FVV1 = SVV(I,J,L,N)
         FWW1 = SWW(I,J,L,N)
         FUV1 = SUV(I,J,L,N)
         FUW1 = SUW(I,J,L,N)
         FVW1 = SVW(I,J,L,N)
         Q2 = AIRUV(I,J+1)
         FTT2 = STT(I,J+1,L,N)
         FUT2 = SUT(I,J+1,L,N)
         FVT2 = SVT(I,J+1,L,N)
         FWT2 = SWT(I,J+1,L,N)
         FUU2 = SUU(I,J+1,L,N)
         FVV2 = SVV(I,J+1,L,N)
         FWW2 = SWW(I,J+1,L,N)
         FUV2 = SUV(I,J+1,L,N)
         FUW2 = SUW(I,J+1,L,N)
         FVW2 = SVW(I,J+1,L,N)
c---SOMADD subroutine
c---orientation:  box1 on Left, box2 on Right, positive V is to Right
         EPS1 = Q1/(Q1+Q2)
         EPS2 = 1._r8 - EPS1
         FTT     = FTT1 + FTT2
         FVT     = EPS1*FVT1 + EPS2*FVT2 + 3._r8*(EPS1*FTT2 - EPS2*FTT1)
         FUV     = EPS1*FUV1 + EPS2*FUV2 + 3._r8*(EPS1*FUT2 - EPS2*FUT1)
         FVW     = EPS1*FVW1 + EPS2*FVW2 + 3._r8*(EPS1*FWT2 - EPS2*FWT1)
         FVV     = EPS1*EPS1*FVV1 + EPS2*EPS2*FVV2 +
     &        5._r8*(EPS1*EPS2*(FVT2-FVT1)
     &                 + (EPS1-EPS2)*(EPS1*FTT2-EPS2*FTT1))
         FUT     = FUT1 + FUT2
         FWT     = FWT1 + FWT2
         FUU     = FUU1 + FUU2
         FWW     = FWW1 + FWW2
         FUW     = FUW1 + FUW2
c---
        if (J.eq.1) then
          JJ = 2
        else
          JJ = JM-1
        endif
         STT(I,JJ,L,N) = FTT
         SUT(I,JJ,L,N) = FUT
         SVT(I,JJ,L,N) = FVT
         SWT(I,JJ,L,N) = FWT
         SUU(I,JJ,L,N) = FUU
         SVV(I,JJ,L,N) = FVV
         SWW(I,JJ,L,N) = FWW
         SUV(I,JJ,L,N) = FUV
         SUW(I,JJ,L,N) = FUW
         SVW(I,JJ,L,N) = FVW

       enddo
       enddo

      enddo  ! n loop

c---new v 6.0b version of over-the-pole flow combines J=1:2 & JM-1:JM
c---        DYN2UL & DYN2VL only do transport for boxes 2:JM-1
c---  ALFAs must be combined from J=1:2 inot J=2 because DYN2U/V does J=2:JM-1

      do J = 3,JM-2
        do I = 1,IM
          ALFA2D(I,J) = ALFA(I,J,L)
          BETA2D(I,J) = BETA(I,J,L)
        enddo
      enddo
      do I = 1,IM
        AIRUV(I,2) = AIRUV(I,2) + AIRUV(I,1)
        AIRUV(I,JM-1) = AIRUV(I,JM-1) + AIRUV(I,JM)
        ALFA2D(I,2) = ALFA(I,2,L) + ALFA(I,1,L)
        ALFA2D(I,JM-1) = ALFA(I,JM-1,L) + ALFA(I,JM,L)
        BETA2D(I,JM-1) = BETA(I,JM-1,L)
        BETA2D(I,2)  = 0._r8
        BETA2D(I,JM) = 0._r8
      enddo
        ALFA2D(IM+1,2) = ALFA2D(1,2)
        ALFA2D(IM+1,JM-1) = ALFA2D(1,JM-1)

       return
       end


c-----------------------------------------------------------------------
       subroutine POLES2(AIRUV,QFU,QFV,DTADV,L)
c-----------------------------------------------------------------------
c     Split extended polar-pie box back into two (J=1 & 2  J=JM-1 & JM)
c     Use S.O.M. to split, need to know final air mass in J=1 & J=JM
c-----------------------------------------------------------------------
c  new CTM p-code 6.0b (5/2011) - OMP parallelized over LM(layer)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, NPAR
      use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW,
     &     AIR, ALFA, BETA, NTM, AREAXY
      implicit none
C-----------------------------------------------------------------------

       real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
       real(r8),  intent(inout) :: QFU(IPAR+1,JPAR,NPAR,2)
       real(r8),  intent(inout) :: QFV(IPAR,JPAR+1,NPAR,2)
       real(r8),  intent(in) :: DTADV
       integer, intent(in) :: L                     ! Layer index

       real(r8) :: FTT,FUT,FVT,FWT,FUU,FVV,FWW,FUV,FUW,FVW
       real(r8) ::
     &      FTT1,FUT1,FVT1,FWT1,FUU1,FVV1,FWW1,FUV1,FUW1,FVW1,Q1,EPS1
       real(r8) ::
     &      FTT2,FUT2,FVT2,FWT2,FUU2,FVV2,FWW2,FUV2,FUW2,FVW2,Q2,EPS2
       real(r8) :: F0,F1,F2, AIR1,AIR2
       integer :: I,J,JJ,N
       integer :: IM, JM        !Locals for CTM3
c-----------------------------------------------------------------------

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR

c---calculate new air mass in J=1 & J=JM (currently merged into J=2 & J=JM-1)
       do I = 1,IM
         AIRUV(I,1) = AIR(I,1,L) + DTADV*
     &      (-BETA(I,2,L) + ALFA(I,1,L) - ALFA(I+1,1,L))
         AIRUV(I,JM) = AIR(I,JM,L) + DTADV*
     &      (+BETA(I,JM,L) + ALFA(I,JM,L) - ALFA(I+1,JM,L))
       enddo

       do N = 1,NTM

       do J = 1,JM-1,JM-2
       do I = 1,IM

        if (J.eq.1) then
          JJ = 2
          EPS1 = AIRUV(I,1)/AIRUV(I,2)
          EPS2 = 1._r8 - EPS1
        else
          JJ = JM-1
          EPS2 = AIRUV(I,JM)/AIRUV(I,JM-1)
          EPS1 = 1._r8 - EPS2
        endif

         FTT = STT(I,JJ,L,N)
         FUT = SUT(I,JJ,L,N)
         FVT = SVT(I,JJ,L,N)
         FWT = SWT(I,JJ,L,N)
         FUU = SUU(I,JJ,L,N)
         FVV = SVV(I,JJ,L,N)
         FWW = SWW(I,JJ,L,N)
         FUV = SUV(I,JJ,L,N)
         FUW = SUW(I,JJ,L,N)
         FVW = SVW(I,JJ,L,N)
c---must apply LIMTR=2 (positive monotonic) before splitting box
         F0  = FTT
         F1  = FVT
         F2  = FVV
           if (F2 .gt. 0._r8) then
               F2  = min(F2, abs(F1)/3._r8, 0.5_r8*F0)
            if (F1 .lt. 0._r8) then
               F1  = max(-(F0+F2), F1)
            else
               F1  = min(+(F0+F2), F1)
            endif
           else
               F1  = min(+F0 ,max(-F0, F1))
               F2  = max(F2, -F0+abs(F1), -abs(F1)/3._r8)
           endif
         FVT = F1
         FVV = F2
         FUV = min(+F0, max(-F0, FUV))
         FVW = min(+F0, max(-F0, FvW))
c-----
         FTT1 = EPS1*(FTT - EPS2*(FVT + (EPS1-EPS2)*FVV))
         FVT1 = EPS1*EPS1*(FVT - 3._r8*EPS2*FVV)
         FVV1 = EPS1*EPS1*EPS1*FVV
         FUT1 = EPS1*(FUT - EPS2*FUV)
         FWT1 = EPS1*(FWT - EPS2*FVW)
         FUV1 = EPS1*EPS1*FUV
         FVW1 = EPS1*EPS1*FVW
         FUU1 = EPS1*FUU
         FWW1 = EPS1*FWW
         FUW1 = EPS1*FUW
c-----
         FTT2 = EPS2*(FTT + EPS1*(FVT + (EPS1-EPS2)*FVV))
         FVT2 = EPS2*EPS2*(FVT + 3._r8*EPS1*FVV)
         FVV2 = EPS2*EPS2*EPS2*FVV
         FUT2 = EPS2*(FUT + EPS1*FUV)
         FWT2 = EPS2*(FWT + EPS1*FVW)
         FUV2 = EPS2*EPS2*FUV
         FVW2 = EPS2*EPS2*FVW
         FUU2 = EPS2*FUU
         FWW2 = EPS2*FWW
         FUW2 = EPS2*FUW
c---
         STT(I,J,L,N) = FTT1
         SUT(I,J,L,N) = FUT1
         SVT(I,J,L,N) = FVT1
         SWT(I,J,L,N) = FWT1
         SUU(I,J,L,N) = FUU1
         SVV(I,J,L,N) = FVV1
         SWW(I,J,L,N) = FWW1
         SUV(I,J,L,N) = FUV1
         SUW(I,J,L,N) = FUW1
         SVW(I,J,L,N) = FVW1
c---
         STT(I,J+1,L,N) = FTT2
         SUT(I,J+1,L,N) = FUT2
         SVT(I,J+1,L,N) = FVT2
         SWT(I,J+1,L,N) = FWT2
         SUU(I,J+1,L,N) = FUU2
         SVV(I,J+1,L,N) = FVV2
         SWW(I,J+1,L,N) = FWW2
         SUV(I,J+1,L,N) = FUV2
         SUW(I,J+1,L,N) = FUW2
         SVW(I,J+1,L,N) = FVW2

       enddo
       enddo

       enddo

c---partition the QFU & QFV fluxes into the polar-pie boxes (J=1 & JM)
       do N = 1,NTM
        do I = 1,IM
           EPS1 = AIRUV(I,1)/AIRUV(I,2)
         QFV(I,1,N,1) = QFV(I,2,N,1)
         QFV(I,1,N,2) = QFV(I,2,N,2)
         QFV(I,2,N,1) = EPS1*QFV(I,3,N,1)
         QFV(I,2,N,2) = EPS1*QFV(I,3,N,2)
           EPS2 = AIRUV(I,JM)/AIRUV(I,JM-1)
         QFV(I,JM+1,N,1) = QFV(I,JM,N,1)
         QFV(I,JM+1,N,2) = QFV(I,JM,N,2)
         QFV(I,JM,N,1) = EPS2*QFV(I,JM-1,N,1)
         QFV(I,JM,N,2) = EPS2*QFV(I,JM-1,N,2)
        enddo
c---it is hard to partition QFU, just used fixed EPS1 from AREAXY
           EPS1 = AREAXY(1,1)/(AREAXY(1,1)+AREAXY(1,2))
           EPS2 = 1._r8 - EPS1
        do I = 1,IM+1
         QFU(I,1,N,1) = EPS1*QFU(I,2,N,1)
         QFU(I,2,N,1) = EPS2*QFU(I,2,N,1)
         QFU(I,1,N,2) = EPS1*QFU(I,2,N,2)
         QFU(I,2,N,2) = EPS2*QFU(I,2,N,2)
         QFU(I,JM,N,1) = EPS1*QFU(I,JM-1,N,1)
         QFU(I,JM-1,N,1) = EPS2*QFU(I,JM-1,N,1)
         QFU(I,JM,N,2) = EPS1*QFU(I,JM-1,N,2)
         QFU(I,JM-1,N,2) = EPS2*QFU(I,JM-1,N,2)
        enddo
       enddo

c---reset air mass for J=2 & J=JM-1 after N-loops
       do I = 1,IM
         AIRUV(I,2) = AIRUV(I,2) - AIRUV(I,1)
         AIRUV(I,JM-1) = AIRUV(I,JM-1) - AIRUV(I,JM)
       enddo

       return
       end

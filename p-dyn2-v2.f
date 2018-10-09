c-----------------------------------------------------------------------
c---(p-dyn2.f)---(p-code 6.1b 2/2012):
c---   6/2011: fixed moment sign error in DYN2V !!
c---    Alternate 0b - retain X-moments in the polar-pie boxes (5/2011)
c
c              U-V tracer advection parallel over layers L=1:LM
c              W   tracer adv parallel over IJ-blocks
c              W   tracer adv also combines large-scale W & convective W
c
c---subroutines:  DYN2UL, DYN2VL, DYN2W_OC, QLIMIT2
c             ***have dropped non-L, parallel-over-tracers versions**


c-----------------------------------------------------------------------
      subroutine DYN2UL(DTALFA,AIRUV,L,QFU,USTEP)
c-----------------------------------------------------------------------
c       WEST-to-EAST ADVECTION OF TRACE COMPOUNDS using S.O.M.
c-----------------------------------------------------------------------
c  new CTM p-code 5.3 (7/2007) - OMP parallelized over LM(layer)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, NPAR
      use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW,
     &     NTM, LMTSOM
      implicit none
C-----------------------------------------------------------------------

      real(r8),  intent(in) :: DTALFA
      real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
      integer, intent(in) :: L                     ! Layer index
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
      real(r8) :: AIRUVIJ(IPAR,JPAR)
      integer :: NQ       ! length of vector pipe for advection (assumed cyclic)
      integer :: NSTEP    ! #multi-steps needed for local CFL, ret by QVECT3
      integer :: I,J,N
      integer :: IM, JM   !Locals for CTM3
c-----------------------------------------------------------------------

c>>>this check on negative airmass SHOULD be turned off after debug mode
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
        do J = 1,JM

c---transfer STT(), SXT(), ... into local vector for piped flow
          do I = 1,IM
            QM(I)  = AIRUV(I,J)
            QU(I)  = ALFA(I,J,L) * DTALFA
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
        do J = 1,JM
          do I = 1,IM
            AIRUV(I,J) = AIRUVIJ(I,J)
          enddo
        enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine DYN2VL(DTBETA,AIRUV,L,QFV,VSTEP)
c-----------------------------------------------------------------------
c       SOUTH-to-NORTH ADVECTION OF TRACE COMPOUNDS using S.O.M.
c-----------------------------------------------------------------------
c  new CTM p-code 5.3 (7/2007) - OMP parallelized over N(tracer)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, NPAR
      use cmn_ctm, only: STT, NTM,
     &     SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW
      implicit none
C-----------------------------------------------------------------------

      real(r8),  intent(in) :: DTBETA
      real(r8),  intent(inout) :: AIRUV(IPAR,JPAR)
      integer, intent(in) :: L                      ! Layer index
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
      integer :: I,J,N,II,JJ,IQD2,JMT2
      integer :: IM, JM   !Locals for CTM3
c-----------------------------------------------------------------------

c Set for local use (IM, JM, LM are removed from CTM3)
      IM = IPAR
      JM = JPAR

c---N-S flow is periodic, over-the-pole,  combines opposite meridions:
      IQD2  = IM/2            ! I loop is half of IM
      JMT2  = JM + JM         ! J loop is twice JM, but
      NQ    = JMT2            !   Alt 2 - leaves polar pie-wedges
      VSTEP = 0
c-----------------------------------------------------------------------

c>>>this check on negative airmass SHOULD be turned off after debug mode
c     do J = 1,JM
c       do I = 1,IM
c         if (AIRUV(I,J) .lt. 0._r8)  then
c           call EXITIJL (' AIRTR < 0 in DYN2VL', I,J,L)
c         endif
c       enddo
c     enddo

      do N = 1,NTM
        do I = 1,IQD2
           II = I + IQD2

c-----------------------------------------------------------------------
c---start at SPole and follow up merdion stripe at I
          do J = 1,JM
            QM(J)  = AIRUV(I,J)
            QU(J)  = BETA(I,J,L)*DTBETA
            QTT(J) = STT(I,J,L,N)
            QXT(J) = SUT(I,J,L,N)
            QYT(J) = SVT(I,J,L,N)
            QZT(J) = SWT(I,J,L,N)
            QXX(J) = SUU(I,J,L,N)
            QYY(J) = SVV(I,J,L,N)
            QZZ(J) = SWW(I,J,L,N)
            QXY(J) = SUV(I,J,L,N)
            QXZ(J) = SUW(I,J,L,N)
            QYZ(J) = SVW(I,J,L,N)
          enddo
c---continue down (sign reversal on U & V moments, QU) merdion stripe II
          do J = JM+1,JMT2
             JJ = JMT2+1 - J
            QM(J) = AIRUV(II,JJ)
            QU(J)  = -BETA(II,JJ+1,L)*DTBETA
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
            QU(JMT2+1) = QU(1)

c         write(6,'(A,3I4)') 'call QVECT3 in DYN2VL N/I/II',N,I,II
c-----------------------------------------------------------------------
          call QVECT3(QM,QU,NQ,LMTSOM, NSTEP,
     &           QTT,QYT,QYY,QYZ,QXY,QZT,QZZ,QXT,QXX,QXZ, Q0F,Q1F    )
          VSTEP = VSTEP + NSTEP
c-----------------------------------------------------------------------

c---replace meridion stripe up I
          do J = 1,JM
            STT(I,J,L,N) = QTT(J)
            SUT(I,J,L,N) = QXT(J)
            SVT(I,J,L,N) = QYT(J)
            SWT(I,J,L,N) = QZT(J)
            SUU(I,J,L,N) = QXX(J)
            SVV(I,J,L,N) = QYY(J)
            SWW(I,J,L,N) = QZZ(J)
            SUV(I,J,L,N) = QXY(J)
            SUW(I,J,L,N) = QXZ(J)
            SVW(I,J,L,N) = QYZ(J)
            AIRUVIJ(I,J) = QM(J)
          enddo
          do J = 1,JM+1
            QFV(I,J,N,1)  = Q0F(J)
            QFV(I,J,N,2)  = Q1F(J)
          enddo
c---replace meridion stripe down II (note reverse of Y moments)
          do J = JM+1,JMT2
            JJ = JMT2+1 - J
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
          enddo
          do J = JM+1,JMT2+1
            JJ = JMT2+1 - J
            QFV(II,JJ+1,N,1) = -Q0F(J)
            QFV(II,JJ+1,N,2) = -Q1F(J)
          enddo

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
     &     QTT,           ! tracer mass in box [I]
     &     QXT,QYT,QZT,   ! 1st moments of tracer in U, V, W direction
     &     QXX,QYY,QZZ,   ! 2nd moments of tracer in U, V, W direction
     &     QXY,QYZ,QXZ,   ! cross-moments of tracer
     &     Q0F,Q1F        ! computed tracer flux from [I-1] to [I]
      integer :: NQ       ! length of vector pipe for advection (assumed cyclic)
      integer :: NSTEP    ! #multi-steps needed for local CFL, ret by QVECT3
      integer :: I,J,L,II,JJ,NB,M,ML,N,NREST,NXTRA
      integer :: LM !Locals for CTM3
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

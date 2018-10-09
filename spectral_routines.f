!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Collection of routines used for converting spectral fields to
!// grid point data.
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
C-----------------------------------------------------------------------
C Routines for converting spectral data to gridpoint data.
C These are not put in f90 module, because their arguments typically
C have different array structure than the calling routine.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      Subroutine SPE2GP ( SPE,NRNM,LUV,GPT,ID,JD,LD, ALPGAU,NMMAX,
     &                    FFTAR,WORK,TRIG,IFAX,NTP1,FFTPTS,NRFFT  )
C-----------------------------------------------------------------------
C Do the spectral to gridpoint transform. First do an inverse Legendre
C transform (ILT) then do a Fast Fourier Transform (FFT).
C
C Original version    J.K.Sundet January 1995
C-----------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      use utilities, only: ctmExitC
      Implicit NONE
C-----------------------------------------------------------------------
C I/O parameters
      Integer
     &     IFAX(10),
     &     NMMAX,NRNM,NTP1,
     &     FFTPTS,NRFFT,
     &     ID,JD,LD
      real(r8)
     &     GPT(ID,JD,LD),       ! gaussian gridpt values or fluxes
     &     FFTAR(FFTPTS,NRFFT), ! fourier to gridpoint transf. array
     &     WORK(FFTPTS*NRFFT),  ! work area for FFT
     &     ALPGAU(NMMAX,JD/2),
     &     TRIG(ID)
      Real(r4)
     &     SPE(2,NRNM/2,LD)     ! spectral coeffs
      Logical
     &     LUV                  ! .TRUE. if U or V
C Local parameters
      real(r8)
     &     SPEWRK(2,51680),     ! 51680 = T319 size (NRNM/2) Obs:NRNM is 2*NMMAX
     &     SUMRES,SUMIMS,       ! symetric real and img sums
     &     SUMREA,SUMIMA        ! anti-symetric real and img sums
      Integer
     &     NGROUP,              ! number of *(32) latitudes
     &     NRES,                ! residual number of latitudes
     &     IGRP,                ! NGROUP counter
     &     ILIM,                ! NTP1
     &     IMLIM,               ! NTP1
     &     IN,                  ! n Max(dimension) for given m
     &     JBACK,               ! counter keeping track of where in NGROUP
     &     JNH,JSH,JLAT,        ! latitude counter
     &     I,J,ILONG,JF,L,INL,  ! loop counters
     &     IMN,IMP,             ! m,n counters
     &     IXTRA                ! an extra counter if U or V
C-----------------------------------------------------------------------
C Initialize
      ILIM  = NTP1
      IMLIM = NTP1
      If (.NOT.LUV) Then
        IXTRA = 0
      Else
        IXTRA = 1
      Endif

      NRES = Mod(JD,NRFFT)
      NGROUP = JD/NRFFT
      If (NRES.NE.0) Then
        Write(*,*) ' * Only programed for multiples of 32 ! *'
        WRITE(*,*) NGROUP, NRES
        Call ctmExitC('=================SP2GP================')
      Endif

C Do the transform in two operations, multiply the spec coeffs
C with appropriate functions. Then sum the ones (odd and even) to
C find the symetric and anti-symetric functions (rel. to equator).

      Do 109 L=1,LD

C For all pairs of NRFFT latitudes
        Do 105 IGRP=1,NGROUP
          Do J=1,NRFFT
            Do I=1,FFTPTS
              FFTAR(I,J) = 0._r8
              WORK(I + (J-1)*FFTPTS) = 0._r8
            End Do
          End Do
          If (IGRP.EQ.1) Then
            JBACK = 0
          Else
            JBACK = ((IGRP-1)*NRFFT)/2
          Endif
C Real(r4) --> real(r8)
          Do I=1,NRNM/2
            SPEWRK(1,I) = SPE(1,I,L)
            SPEWRK(2,I) = SPE(2,I,L)
          End Do
          Do J=1,NRFFT/2
            IMP = 0
            IMN = 0
C For each zonal wavenumber sum the coeffs*A-Leg coeffs.
            Do ILONG=1,IMLIM
              SUMRES = 0._r8
              SUMIMS = 0._r8
              SUMREA = 0._r8
              SUMIMA = 0._r8
              IN  = ILIM - ILONG + (1+IXTRA)
C Find even (symetric) part
              Do JF=1,IN,2
                SUMRES = SUMRES +ALPGAU(IMP+JF,J+JBACK)*SPEWRK(1,IMN+JF)
                SUMIMS = SUMIMS +ALPGAU(IMP+JF,J+JBACK)*SPEWRK(2,IMN+JF)
              Enddo
C Find odd (anti-symetric) part
              Do JF=2,IN,2
                SUMREA = SUMREA +ALPGAU(IMP+JF,J+JBACK)*SPEWRK(1,IMN+JF)
                SUMIMA = SUMIMA +ALPGAU(IMP+JF,J+JBACK)*SPEWRK(2,IMN+JF)
              Enddo
C For the southern hemisphere row, the legendre functions are
C the complex conjugates of the corresponding northern row -
C hence the juggling with the signs.
              FFTAR(2*ILONG-1,J*2-1) = SUMRES + SUMREA
              FFTAR(2*ILONG,  J*2-1) = SUMIMS + SUMIMA
              FFTAR(2*ILONG-1,J*2  ) = SUMRES - SUMREA
              FFTAR(2*ILONG,  J*2)   = SUMIMS - SUMIMA
C
              IMP = IMP + IN + (1-IXTRA)
              IMN = IMN + IN
            Enddo
          Enddo
          
C Do multple FFT tranforms - currently 32, 
          CALL FFT_99(FFTAR,WORK,TRIG,IFAX,1,FFTPTS,ID,NRFFT)
          
C The latitudes are arranged: odd NH, even SH; rearrange
          Do J=1,NRFFT/2
            JNH  = J*2 - 1
            JSH  = J*2
            JLAT = JD - (J+JBACK) + 1
            Do INL=1,ID
              GPT(INL,J+JBACK,L) = FFTAR(INL,JNH)
              GPT(INL,JLAT,L)    = FFTAR(INL,JSH)
            Enddo
          Enddo
 105    Continue
 109  Continue
C
      Return
      End
C
C----------------------------------------------------------------------
      Subroutine ZD2UV(ZSPE,DSPE,U,V,DD,SS,ID,JD,LD,ALP,ALPV,NMMAX,NRNM,
     &                 FFTAR,WORK,TRIG,  IFAX,NTRUNW,NRFFT,FFTPTS)
C----------------------------------------------------------------------
C Take care of the transform to horizontal wind from vorticity and
C divergence fields.
C----------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      Implicit NONE
C----------------------------------------------------------------------
C I/O parameters
      Integer
     &     IFAX(10),
     &     ID,JD,LD,
     &     NRNM,NMMAX,NTRUNW,
     &     NRFFT,FFTPTS
      real(r8)
     &     U(ID,JD,LD),
     &     V(ID,JD,LD),
     &     FFTAR(FFTPTS,NRFFT),
     &     WORK(FFTPTS*NRFFT),
     &     ALP(NMMAX,JD/2),
     &     ALPV(NMMAX,JD/2),
     &     DD(NMMAX),
     &     SS(NMMAX),
     &     TRIG(ID)
cccccccccccccccccc     ,   &     TRIGU(ID)
      Real(r4)
     &     ZSPE(2,NRNM/2,LD),
     &     DSPE(2,NRNM/2,LD)
C Local parameters
      Real(r4), Allocatable ::
     &        USPE(:,:,:), VSPE(:,:,:)
      Logical
     &     LUV
C-------------------------------------------------------------------
      Allocate (
     &     VSPE(2,51680,LD),   ! 51680 = T319 size (NRNM/2) Obs:NRNM is 2*NMMAX
     &     USPE(2,51680,LD) )  ! 5885 = T106 NMMAX=(106+1)*(106+4)/2

C Initialize
      LUV = .TRUE.

C Find the spectral coeffs for U and V from the vorticity
C and divergence coeffs
      Call UVCOEF(ZSPE,DSPE,DD,SS,NRNM,NMMAX,LD,NTRUNW,USPE,VSPE)

C Transform the fields to GridPoints
      Call SPE2GP(USPE,NMMAX*2,LUV, U,ID,JD,LD, ALP,NMMAX,
     &            FFTAR,WORK,TRIG,IFAX,(NTRUNW+1),FFTPTS,NRFFT)
      Call SPE2GP(VSPE,NMMAX*2,LUV, V,ID,JD,LD, ALPV,NMMAX,
     &            FFTAR,WORK,TRIG,IFAX,(NTRUNW+1),FFTPTS,NRFFT)
C
      Deallocate ( VSPE, USPE )
c
      Return
      End
C
C
C----------------------------------------------------------------------
      Subroutine UVCOEF(ZSPE,DSPE,DD,SS,NRNM,NMMAX,LD,NTRUNW,USPE,VSPE)
C----------------------------------------------------------------------
C Instead of using the vorticity and divergence we use the spectral 
C coeffisients for u and v when the fluxes are calculated. It's also
C possible to derive the transport directly from the divergence and 
C vorticity, but the approach chosen is more accessible and cleaner.
C
C Method: The spectral coeffisients for u and v are related to vorticity
C         and divergence in the following way;
C
C   u(n,m)=+dd(n,m)*vor(n-1,m)+i*ss(n,m)*div(n,m)-dd(n+1,m)*vor(n+1,m)
C          and
C   v(n,m)=-dd(n,m)*div(n-1,m)+i*ss(n,m)*vor(n,m)-dd(n+1,m)*div(n+1,m)
C 
C Original version:    J. K. Sundet, November 1994
C----------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      Implicit NONE
C----------------------------------------------------------------------
C I/O parameters
      Integer
     &     NRNM,NMMAX,
     &     LD,NTRUNW
      real(r8)
     &     DD(NMMAX),
     &     SS(NMMAX)
      Real(r4)
     &     USPE(2,NMMAX,LD),
     &     VSPE(2,NMMAX,LD),
     &     ZSPE(2,NRNM/2,LD),
     &     DSPE(2,NRNM/2,LD)
C Local parameters
      real(r8), Allocatable ::
     &        ZWRK(:,:), DWRK(:,:)
      Integer
     &     JL,JN,I,L,         ! loop counters
     &     INU,INM,           ! coeffs counters
     &     NTP1
C----------------------------------------------------------------------
      Allocate (
     &     DWRK(2,51680),        ! 51680 = T319 NMMAX
     &     ZWRK(2,51680) )       ! 5885 = T106 NMMAX=(106+1)*(106+4)/2
C Initialize      
      NTP1 = NTRUNW + 1
      ZWRK(:,:) = 0._r8
      DWRK(:,:) = 0._r8

C For all levels
      Do L=1,LD

C Real(r4) --> real(r8)
         Do I=1,NRNM/2
           ZWRK(1,I) = ZSPE(1,I,L)
           ZWRK(2,I) = ZSPE(2,I,L)
           DWRK(1,I) = DSPE(1,I,L)
           DWRK(2,I) = DSPE(2,I,L)
         End Do

        INU = 1
        INM = 1
        Do JL=1,NTP1
C For  n = m
          If (JL.EQ.1) Then
            USPE(1,INU,L) = -DD(INU+1)*ZWRK(1,INM+1)
            USPE(2,INU,L) = -DD(INU+1)*ZWRK(2,INM+1)
            VSPE(1,INU,L) = +DD(INU+1)*DWRK(1,INM+1)
            VSPE(2,INU,L) = +DD(INU+1)*DWRK(2,INM+1)
          Elseif (JL.EQ.NTP1) Then
            USPE(1,INU,L) = -SS(INM)*DWRK(2,INM)
            USPE(2,INU,L) = +SS(INM)*DWRK(1,INM)
            VSPE(1,INU,L) = -SS(INM)*ZWRK(2,INM)
            VSPE(2,INU,L) = +SS(INM)*ZWRK(1,INM)
          Else
            USPE(1,INU,L) = -SS(INM)*DWRK(2,INM)-DD(INU+1)*ZWRK(1,INM+1)
            USPE(2,INU,L) = +SS(INM)*DWRK(1,INM)-DD(INU+1)*ZWRK(2,INM+1)
            VSPE(1,INU,L) = -SS(INM)*ZWRK(2,INM)+DD(INU+1)*DWRK(1,INM+1)
            VSPE(2,INU,L) = +SS(INM)*ZWRK(1,INM)+DD(INU+1)*DWRK(2,INM+1)
          Endif
          INM = INM + 1
          INU = INU + 1
C For  m < n < NT
          Do JN=(JL+1),NTRUNW
            USPE(1,INU,L) = +DD(INU)*ZWRK(1,INM-1)
     &                     -SS(INM)*DWRK(2,INM)
     &                     -DD(INU+1)*ZWRK(1,INM+1)
            USPE(2,INU,L) = +DD(INU)*ZWRK(2,INM-1)
     &                     +SS(INM)*DWRK(1,INM)
     &                     -DD(INU+1)*ZWRK(2,INM+1)
            VSPE(1,INU,L) = -DD(INU)*DWRK(1,INM-1)
     &                     -SS(INM)*ZWRK(2,INM)
     &                     +DD(INU+1)*DWRK(1,INM+1)
            VSPE(2,INU,L) = -DD(INU)*DWRK(2,INM-1)
     &                     +SS(INM)*ZWRK(1,INM)
     &                     +DD(INU+1)*DWRK(2,INM+1)
            INM = INM + 1
            INU = INU + 1
          Enddo
C For  n = NT
          If (JL.LT.NTP1) Then
            USPE(1,INU,L) = +DD(INU)*ZWRK(1,INM-1)
     &                     -SS(INM)*DWRK(2,INM)
            USPE(2,INU,L) = +DD(INU)*ZWRK(2,INM-1)
     &                     +SS(INM)*DWRK(1,INM)
            VSPE(1,INU,L) = -DD(INU)*DWRK(1,INM-1)
     &                     -SS(INM)*ZWRK(2,INM)
            VSPE(2,INU,L) = -DD(INU)*DWRK(2,INM-1)
     &                     +SS(INM)*ZWRK(1,INM)
            INU = INU + 1
            INM = INM + 1
          Endif
C For  n = NTP1
          USPE(1,INU,L) = +DD(INU)*ZWRK(1,INM-1)
          USPE(2,INU,L) = +DD(INU)*ZWRK(2,INM-1)
          VSPE(1,INU,L) = -DD(INU)*DWRK(1,INM-1)
          VSPE(2,INU,L) = -DD(INU)*DWRK(2,INM-1)
          INU = INU + 1
        Enddo
      Enddo
c
      Deallocate ( ZWRK, DWRK )
C
      Return
c-----------------------------------------------------------------------
      End
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
      subroutine FFT_99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT)
c-----------------------------------------------------------------------
c  MULTIPLE FAST REAL PERIODIC TRANSFORM OF LENGTH N  PERFORMED BY 
c     REMOVING REDUNDANT OPERATIONS FROM COMPLEX TRANSFORM LENGTH N
c
c     A     = ARRAY CONTAINING INPUT & OUTPUT DATA
c     WORK  = WORK AREA OF SIZE (N+1)*MIN(LOT,64)
c     TRIGS = PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
c     IFAX  = PREVIOUSLY PREPARED LIST OF FACTORS OF N
c     INC   = INCREMENT WITHIN EACH DATA VECTOR (INC=1 = CONSEC STORED DATA)
c     JUMP  = INCREMENT BETWEEN THE START OF EACH DATA VECTOR
c     N     = LENGTH OF THE DATA VECTORS
c     LOT   = NUMBER OF DATA VECTORS
c assumed for this version: ISIGN = +1 TRANSFORM SPECTRAL=>GRIDPOINT
c ordering of coeff's
c         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
c         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS ReqUIRED
c ordering of data
c         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) 
c         N must be factor into 2 x 3 x 5 's (odd OK)
c
c transforms:
c     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
c         where C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
c     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
c               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
c
c Adapted to CTM J.K. Sundet Jan. 1996
c-----------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      use utilities, only: ctmExitIJL
      implicit none
      integer, intent(inout)      ::  INC,JUMP,N,LOT
      integer, dimension(10), intent(inout)      :: IFAX
      real(r8), dimension(JUMP*LOT), intent(inout) :: A,WORK
      real(r8), dimension(N), intent(inout)        :: TRIGS

      integer
     &     IOC,                 ! = JUMP
     &     NFAX,                ! number of FFT factors
     &     NX,NBLOX,NVEX,       ! loop limits
     &     ISTART,IA,LA,        ! counters
     &     I,J,K,NB,II,JJ,      ! counters
     &     IGO,                 ! determines what half-plane to sample points
     &     IFAC,IERR,           ! = IFAX(K) and error check
     &     IBASE,JBASE,IX       ! counters
      logical  LNODD, LFAXODD
c-----------------------------------------------------------------------

c--Check that SET_FFT is done:
      if (IFAX(10).ne.N) call FFT_SET(TRIGS,IFAX,N)

      IOC     = JUMP
      NFAX    = IFAX(1)
      LFAXODD = mod(NFAX,2).ne.0
      if (mod(N,2).eq.1) then
        NX    = N
        LNODD = .TRUE.
      else
        NX    = N+1
        LNODD = .FALSE.
      endif
      NBLOX = (LOT-1)/64 + 1
      NVEX  = LOT - (NBLOX-1)*64
      ISTART = 1
      
c---Spectral-to-gridpoint transform (ISIGN=+1 in old notation)
      do NB = 1,NBLOX
          IA = ISTART
          I  = ISTART
        do J=1,NVEX
          A(I+INC) = 0.5_r8*A(I)
          I        = I + JUMP
        enddo
        if (.not.LNODD) then
            I = ISTART + N*INC
          do J=1,NVEX
            A(I) = 0.5_r8*A(I)
            I    = I + JUMP
          enddo
        endif
         IA = ISTART + INC
         LA = 1
         IGO = +1

c---Sample gridpoints on both half-planes
        do K=1,NFAX
          IFAC = IFAX(K+1)
          IERR = -1
         if (IGO.eq.+1) then
            CALL FFT_RPASS(A(IA),A(IA+LA*INC),WORK(1),WORK(IFAC*LA+1),
     &           TRIGS,INC,1,JUMP,NX,NVEX,N,IFAC,IOC,LA,IERR)
            IGO = -1
         else
            CALL FFT_RPASS(WORK(1),WORK(LA+1),A(IA),A(IA+IFAC*LA*INC),
     &           TRIGS,1,INC,NX,JUMP,NVEX,N,IFAC,IOC,LA,IERR)
            IGO = +1
         endif
          LA = IFAC*LA
          IA = ISTART
          if (IERR.ne.0) then
      call ctmEXITIJL('>>>>error FFT_RPASS: IERR,NVEX,IFAC',
     &            IERR,NVEX,IFAC)
          endif
        enddo

c---If necessary, copy results back to A()
        if (LFAXODD) then
          IBASE = 1
          JBASE = IA
          do JJ=1,NVEX
            I = IBASE
            J = JBASE
           do II=1,N
              A(J) = WORK(I)
              I = I + 1
              J = J + INC
           enddo
            IBASE = IBASE + NX
            JBASE = JBASE + JUMP
          enddo
        endif

c---Fill in zero's at end
          IX = ISTART + N*INC
        do J=1,NVEX
          A(IX)     = 0._r8
          A(IX+INC) = 0._r8
          IX = IX + JUMP
        enddo
          ISTART = ISTART + NVEX*JUMP
          NVEX = 64

      enddo

      return
      end



c-----------------------------------------------------------------------
      subroutine FFT_RPASS(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,
     &     IOC,LA,IERR)
c-----------------------------------------------------------------------
c---performs one pass thru data as part of multiple real FFT 
c     A    = FIRST REAL INPUT VECTOR   eqUIVALENCE B(1) WITH A (LA*INC1+1)
c     C    = FIRST REAL OUTPUT VECTOR  eqUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)
c     TRIGS=  ECALCULATED LIST OF SIneS & COSIneS
c     INC1 =  ADDRESSING INCREMENT FOR A
c     INC2 =  ADDRESSING INCREMENT FOR C
c     INC3 =  INCREMENT BETWEEN INPUT VECTORS A
c     INC4 =  INCREMENT BETWEEN OUTPUT VECTORS C
c     LOT  =  NUMBER OF VECTORS
c     N    =  LENGTH OF THE VECTORS
c     IFAC =  CURRENT FACTOR OF N
c     LA   =  PRODUCT OF PREVIOUS FACTORS
c     IERR =  ERROR INDICATOR:
c                0 - PASS COMPLETED WITHOUT ERROR
c                1 - LOT GREATER THAN 64
c                2 - IFAC NOT CATERED FOR
c                3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
c  Adapted to CTM:    J. K. Sundet   July 1994
c------------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      implicit none
      integer, intent(in)  :: INC1,INC2,INC3,INC4,LOT,N,IFAC,IOC,LA
      integer, intent(out) :: IERR
      real(r8), dimension(IOC*LOT), intent(in)  :: A,B
      real(r8), dimension(IOC*LOT), intent(out) :: C,D
      real(r8), dimension(N), intent(in) :: TRIGS

      real(r8), dimension(64) :: A10,A11,A20,A21,B10,B11,B20,B21
      real(r8)  C1,C2,C3,C4,C5,S1,S2,S3,S4,S5,
     &     FN,FIFAC,FLA,FINC1,FINC2,FM
      integer  I,J,K,L,M,IJK, KB,KC,KD,KE,KF
      integer  IA,IB,IC,ID,IE,IG,  JA,JB,JC,JD,JE,JF,JG,JH
      integer  IBASE,JBASE,IINK,JINK,JUMP,KSTOP,IHLP
      logical  LAEQM

      real(r8),parameter:: SIN36=0.587785252292473137_r8
      real(r8),parameter:: SSIN36=2._r8*SIN36
      real(r8),parameter:: SIN45=0.707106781186547524_r8
      real(r8),parameter:: SSIN45=2._r8*SIN45
      real(r8),parameter:: SIN60=0.866025403784438597_r8
      real(r8),parameter:: SSIN60=2._r8*SIN60
      real(r8),parameter:: SIN72=0.951056516295153531_r8
      real(r8),parameter:: SSIN72=2._r8*SIN72
      real(r8),parameter:: QRT5 =0.559016994374947451_r8
      real(r8),parameter:: QQRT5 =2._r8*QRT5

C---------------------------------------------------------------------
      if (LOT.GT.64) then
        IERR = 1
        return
      endif

      IBASE= 0
      JBASE= 0
      IERR = 0
      FN   = real(N,r8)
      FIFAC= real(IFAC,r8)
      FLA  = real(LA,r8)
      FINC1= real(INC1,r8)
      FINC2= real(INC2,r8)
      FM   = FN/FIFAC
      M    = int(FM)
      IINK = int(FLA*FINC1)
      JINK = int(FLA*FINC2)
      JUMP = int((FIFAC-1._r8)*real(JINK,r8))
      KSTOP= int((FN-FIFAC)/(2._r8*FIFAC))
      LAEQM= LA.eq.M

c---defactorize the FFT coeff's, determined by IFAC
c---IFAC=2
      if (IFAC.eq.2) then

        IA = 1
        IB = IA + int((2._r8*FM - FLA)*FINC1)
        JA = 1
        JB = JA + JINK
       if (.not.LAEQM) then
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=A(IA+I)+A(IB+I)
              C(JB+J)=A(IA+I)-A(IB+I)
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo

          IA   = IA+IINK
          IINK = 2*IINK
          IB   = IB-IINK
          IBASE= 0
          JBASE= JBASE+JUMP
          JUMP = 2*JUMP+JINK
         if (IA.ne.IB) then
          do K=LA,KSTOP,LA
            KB=K+K
            C1=TRIGS(KB+1)
            S1=TRIGS(KB+2)
            IBASE=0
           do L=1,LA
             I=IBASE
             J=JBASE
            do IJK=1,LOT
              C(JA+J)=A(IA+I)+A(IB+I)
              D(JA+J)=B(IA+I)-B(IB+I)
              C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
              D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
              I=I+INC3
              J=J+INC4
            enddo
             IBASE=IBASE+INC1
             JBASE=JBASE+INC2
           enddo
            IA=IA+IINK
            IB=IB-IINK
            JBASE=JBASE+JUMP
          enddo

              if(IA.GT.IB) then
               return
              endif

         endif
          IBASE = 0 
         do L=1,LA
           I=IBASE
           J=JBASE
          do IJK=1,LOT
             C(JA+J)=A(IA+I)
             C(JB+J)=-B(IA+I)
             I=I+INC3
             J=J+INC4
          enddo
           IBASE=IBASE+INC1
           JBASE=JBASE+INC2
         enddo
       else
         do L=1,LA
           I=IBASE
           J=JBASE
          do IJK=1,LOT
             C(JA+J)=2._r8*(A(IA+I)+A(IB+I))
             C(JB+J)=2._r8*(A(IA+I)-A(IB+I))
             I=I+INC3
             J=J+INC4
          enddo
           IBASE=IBASE+INC1
           JBASE=JBASE+INC2
         enddo
       endif

c---IFAC=3
      elseif(IFAC.eq.3) then

        IA = 1
        IB = IA+int((2._r8*FM-FLA)*FINC1)
        IC = IB
        JA = 1
        JB = JA+JINK
        JC = JB+JINK
       if (.not.LAEQM) then
         do L=1,LA
            I=IBASE
            J=JBASE
          do IJK=1,LOT
             C(JA+J)=A(IA+I)+A(IB+I)
             C(JB+J)=(A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))
             C(JC+J)=(A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))
             I=I+INC3
             J=J+INC4
           enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         enddo
          IA   = IA+IINK
          IINK = 2*IINK
          IB   = IB+IINK
          IC   = IC-IINK
          JBASE= JBASE+JUMP
          JUMP = 2*JUMP+JINK
         if (IA.ne.IC) then
           do K=LA,KSTOP,LA
             KB=K+K
             KC=KB+KB
             C1=TRIGS(KB+1)
             S1=TRIGS(KB+2)
             C2=TRIGS(KC+1)
             S2=TRIGS(KC+2)
             IBASE=0
            do L=1,LA
               I=IBASE
               J=JBASE
             do IJK=1,LOT
               C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
               D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
               C(JB+J)= C1*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I))) - 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S1*((B(IA+I)-0.5_r8*(B(IB+I)-B(IC+I))) +
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JB+J)= S1*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I))) - 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C1*((B(IA+I)-0.5_r8*(B(IB+I)-B(IC+I))) + 
     &              (SIN60*(A(IB+I)-A(IC+I))))
               C(JC+J)= C2*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I))) + 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S2*((B(IA+I)-0.5_r8*(B(IB+I)-B(IC+I))) -
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JC+J)= S2*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I))) + 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C2*((B(IA+I)-0.5_r8*(B(IB+I)-B(IC+I))) - 
     &              (SIN60*(A(IB+I)-A(IC+I))))
               I = I + INC3
               J = J + INC4
             enddo
              IBASE = IBASE + INC1
              JBASE = JBASE + INC2
            enddo
             IA  = IA + IINK
             IB  = IB + IINK
             IC  = IC - IINK
             JBASE = JBASE + JUMP
           enddo

              if (IA.GT.IC) then
               return
              endif

         endif
          IBASE = 0
         do L=1,LA
            I=IBASE
            J=JBASE
          do IJK=1,LOT
              C(JA+J)=A(IA+I)+A(IB+I)
              C(JB+J)=(0.5_r8*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
              C(JC+J)=-(0.5_r8*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
              I=I+INC3
              J=J+INC4
          enddo
           IBASE=IBASE+INC1
           JBASE=JBASE+INC2
         enddo
       else
         do L=1,LA
           I=IBASE
           J=JBASE
          do IJK=1,LOT
            C(JA+J)=2._r8*(A(IA+I)+A(IB+I))
            C(JB+J)=(2._r8*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
            C(JC+J)=(2._r8*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
            I=I+INC3
            J=J+INC4
          enddo
           IBASE=IBASE+INC1
           JBASE=JBASE+INC2
         enddo
       endif

c---IFAC=4
      elseif(IFAC.eq.4) then

        IA = 1
        IB = IA+int((2._r8*FM-FLA)*FINC1)
        IC = IB+int(2._r8*FM*FINC1)
        ID = IB
        JA = 1
        JB = JA+JINK
        JC = JB+JINK
        JD = JC+JINK
        if (.not.LAEQM) then
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
              C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
              C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
              C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo

          IA   = IA+IINK
          IINK = 2*IINK
          IB   = IB+IINK
          IC   = IC-IINK
          ID   = ID-IINK
          JBASE= JBASE+JUMP
          JUMP = 2*JUMP+JINK
          if(IB.ne.IC) then
            do K=LA,KSTOP,LA
              KB=K+K
              KC=KB+KB
              KD=KC+KB
              C1=TRIGS(KB+1)
              S1=TRIGS(KB+2)
              C2=TRIGS(KC+1)
              S2=TRIGS(KC+2)
              C3=TRIGS(KD+1)
              S3=TRIGS(KD+2)
              IBASE=0
             do L=1,LA
               I=IBASE
               J=JBASE
               do IJK=1,LOT
                 C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
                 D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
                 C(JC+J)= C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     &                -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
                 D(JC+J)= S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     &                +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
                 C(JB+J)= C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     &                -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
                 D(JB+J)= S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     &                +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
                 C(JD+J)= C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     &                -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
                 D(JD+J)= S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     &                +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
                 I=I+INC3
                 J=J+INC4
               enddo
               IBASE=IBASE+INC1
               JBASE=JBASE+INC2
             enddo
              IA=IA+IINK
              IB=IB+IINK
              IC=IC-IINK
              ID=ID-IINK
              JBASE=JBASE+JUMP
            enddo

              if(IB.GT.IC) then
               return
              endif

          endif
          IBASE = 0
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=A(IA+I)+A(IB+I)
              C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
              C(JC+J)=B(IB+I)-B(IA+I)
              C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        else
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=2._r8*((A(IA+I)+A(IC+I))+A(IB+I))
              C(JB+J)=2._r8*((A(IA+I)-A(IC+I))-B(IB+I))
              C(JC+J)=2._r8*((A(IA+I)+A(IC+I))-A(IB+I))
              C(JD+J)=2._r8*((A(IA+I)-A(IC+I))+B(IB+I))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        endif

c---IFAC=5
      elseif(IFAC.eq.5) then
        IA = 1
        IB = IA+int((2._r8*FM-FLA)*FINC1)
        IC = IB+int(2._r8*FM*FINC1)
        ID = IC
        IE = IB
        JA = 1
        JB = JA+JINK
        JC = JB+JINK
        JD = JC+JINK
        JE = JD+JINK
        if(.not.LAEQM) then
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
              C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I))) + 
     &          QRT5*(A(IB+I)-A(IC+I)))-(SIN72*B(IB+I)+SIN36*B(IC+I))
              C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I))) - 
     &          QRT5*(A(IB+I)-A(IC+I)))-(SIN36*B(IB+I)-SIN72*B(IC+I))
              C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I))) - 
     &          QRT5*(A(IB+I)-A(IC+I)))+(SIN36*B(IB+I)-SIN72*B(IC+I))
              C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I))) + 
     &          QRT5*(A(IB+I)-A(IC+I)))+(SIN72*B(IB+I)+SIN36*B(IC+I))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo

           IA   = IA+IINK
           IINK = 2*IINK
           IB   = IB+IINK
           IC   = IC+IINK
           ID   = ID-IINK
           IE   = IE-IINK
           JBASE= JBASE+JUMP
           JUMP = 2*JUMP+JINK
          if(IB.ne.ID) then
            do K=LA,KSTOP,LA
              KB=K+K
              KC=KB+KB
              KD=KC+KB
              KE=KD+KB
              C1=TRIGS(KB+1)
              S1=TRIGS(KB+2)
              C2=TRIGS(KC+1)
              S2=TRIGS(KC+2)
              C3=TRIGS(KD+1)
              S3=TRIGS(KD+2)
              C4=TRIGS(KE+1)
              S4=TRIGS(KE+2)
              IBASE=0
             do L=1,LA
                I=IBASE
                J=JBASE
              do IJK=1,LOT
                A10(IJK)= +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))+ 
     &              (A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))
                A20(IJK)= -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))+
     &             (A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))
                B10(IJK)= +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))+
     &             (B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))
                B20(IJK)= -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))+
     &             (B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))
                A11(IJK)=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
                A21(IJK)=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
                B11(IJK)=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
                B21(IJK)=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))
                C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
                D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
                C(JB+J)=C1*(A10(IJK)-A11(IJK))-S1*(B10(IJK)+B11(IJK))
                D(JB+J)=S1*(A10(IJK)-A11(IJK))+C1*(B10(IJK)+B11(IJK))
                C(JE+J)=C4*(A10(IJK)+A11(IJK))-S4*(B10(IJK)-B11(IJK))
                D(JE+J)=S4*(A10(IJK)+A11(IJK))+C4*(B10(IJK)-B11(IJK))
                C(JC+J)=C2*(A20(IJK)-A21(IJK))-S2*(B20(IJK)+B21(IJK))
                D(JC+J)=S2*(A20(IJK)-A21(IJK))+C2*(B20(IJK)+B21(IJK))
                C(JD+J)=C3*(A20(IJK)+A21(IJK))-S3*(B20(IJK)-B21(IJK))
                D(JD+J)=S3*(A20(IJK)+A21(IJK))+C3*(B20(IJK)-B21(IJK))
                I=I+INC3
                J=J+INC4
              enddo
                IBASE=IBASE+INC1
                JBASE=JBASE+INC2
             enddo
              IA=IA+IINK
              IB=IB+IINK
              IC=IC+IINK
              ID=ID-IINK
              IE=IE-IINK
              JBASE=JBASE+JUMP
            enddo

              if(IB.GT.ID) then
               return
              endif

          endif

          IBASE = 0
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
              C(JB+J)=(QRT5*(A(IA+I)-A(IB+I)) + 
     &             (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &             -(SIN36*B(IA+I)+SIN72*B(IB+I))
              C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I)) + 
     &             (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &             -(SIN36*B(IA+I)+SIN72*B(IB+I))
              C(JC+J)=(QRT5*(A(IA+I)-A(IB+I)) - 
     &             (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &             -(SIN72*B(IA+I)-SIN36*B(IB+I))
              C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I)) - 
     &             (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &             -(SIN72*B(IA+I)-SIN36*B(IB+I))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        else
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=2._r8*(A(IA+I)+(A(IB+I)+A(IC+I)))
              C(JB+J)=(2._r8*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             +QQRT5*(A(IB+I)-A(IC+I))) - 
     &             (SSIN72*B(IB+I)+SSIN36*B(IC+I))
              C(JC+J)=(2._r8*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             -QQRT5*(A(IB+I)-A(IC+I))) -
     &             (SSIN36*B(IB+I)-SSIN72*B(IC+I))
              C(JD+J)=(2._r8*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             -QQRT5*(A(IB+I)-A(IC+I))) + 
     &             (SSIN36*B(IB+I)-SSIN72*B(IC+I))
              C(JE+J)=(2._r8*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             +QQRT5*(A(IB+I)-A(IC+I))) + 
     &             (SSIN72*B(IB+I)+SSIN36*B(IC+I))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        endif

c---IFAC=6
      elseif (IFAC.eq.6) then

        IHLP= int(2._r8*FM*FINC1)
        IA  = 1
        IB  = IA+int((2._r8*FM-FLA)*FINC1)
        IC  = IB+IHLP
        ID  = IC+IHLP
        IE  = IC
        IG  = IB
        JA  = 1
        JB  = JA+JINK
        JC  = JB+JINK
        JD  = JC+JINK
        JE  = JD+JINK
        JF  = JE+JINK
        if (.not.LAEQM) then 
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
              C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
              C(JB+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))
     &             -(SIN60*(B(IB+I)+B(IC+I)))
              C(JF+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))
     &             +(SIN60*(B(IB+I)+B(IC+I)))
              C(JC+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))
     &             -(SIN60*(B(IB+I)-B(IC+I)))
              C(JE+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))
     &             +(SIN60*(B(IB+I)-B(IC+I)))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo

          IA   = IA+IINK
          IINK = 2*IINK
          IB   = IB+IINK
          IC   = IC+IINK
          ID   = ID-IINK
          IE   = IE-IINK
          IG   = IG-IINK
          JBASE= JBASE+JUMP
          JUMP = 2*JUMP+JINK
          if (IC.ne.ID) then
            do K=LA,KSTOP,LA
              KB=K+K
              KC=KB+KB
              KD=KC+KB
              KE=KD+KB
              KF=KE+KB
              C1=TRIGS(KB+1)
              S1=TRIGS(KB+2)
              C2=TRIGS(KC+1)
              S2=TRIGS(KC+2)
              C3=TRIGS(KD+1)
              S3=TRIGS(KD+2)
              C4=TRIGS(KE+1)
              S4=TRIGS(KE+2)
              C5=TRIGS(KF+1)
              S5=TRIGS(KF+2)
              IBASE=0
              do L=1,LA
                I=IBASE
                J=JBASE
                do IJK=1,LOT
                  A11(IJK)= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IG+I))
                  A20(IJK)=(A(IA+I)+A(ID+I))-0.5*A11(IJK)
                  A21(IJK)=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IG+I)))
                  B11(IJK)= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IG+I))
                  B20(IJK)=(B(IA+I)-B(ID+I))-0.5*B11(IJK)
                  B21(IJK)=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IG+I)))
                  C(JA+J)=(A(IA+I)+A(ID+I))+A11(IJK)
                  D(JA+J)=(B(IA+I)-B(ID+I))+B11(IJK)
                  C(JC+J)=C2*(A20(IJK)-B21(IJK))-S2*(B20(IJK)+A21(IJK))
                  D(JC+J)=S2*(A20(IJK)-B21(IJK))+C2*(B20(IJK)+A21(IJK))
                  C(JE+J)=C4*(A20(IJK)+B21(IJK))-S4*(B20(IJK)-A21(IJK))
                  D(JE+J)=S4*(A20(IJK)+B21(IJK))+C4*(B20(IJK)-A21(IJK))
                  A11(IJK)=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IG+I))
                  B11(IJK)=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IG+I))
                  A20(IJK)=(A(IA+I)-A(ID+I))-0.5*A11(IJK)
                  A21(IJK)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IG+I)))
                  B20(IJK)=(B(IA+I)+B(ID+I))+0.5*B11(IJK)
                  B21(IJK)=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IG+I)))
                  C(JD+J)=C3*((A(IA+I)-A(ID+I))+A11(IJK)) -
     &                 S3*((B(IA+I)+B(ID+I))-B11(IJK))
                  D(JD+J)=S3*((A(IA+I)-A(ID+I))+A11(IJK)) + 
     &                 C3*((B(IA+I)+B(ID+I))-B11(IJK))
                  C(JB+J)=C1*(A20(IJK)-B21(IJK))-S1*(B20(IJK)-A21(IJK))
                  D(JB+J)=S1*(A20(IJK)-B21(IJK))+C1*(B20(IJK)-A21(IJK))
                  C(JF+J)=C5*(A20(IJK)+B21(IJK))-S5*(B20(IJK)+A21(IJK))
                  D(JF+J)=S5*(A20(IJK)+B21(IJK))+C5*(B20(IJK)+A21(IJK))

                  I=I+INC3
                  J=J+INC4
                enddo
                IBASE=IBASE+INC1
                JBASE=JBASE+INC2
              enddo
              IA=IA+IINK
              IB=IB+IINK
              IC=IC+IINK
              ID=ID-IINK
              IE=IE-IINK
              IG=IG-IINK
              JBASE=JBASE+JUMP
            enddo

              if (IC.GT.ID) then
               return
              endif

          endif

          IBASE = 0
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
              C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
              C(JB+J)=(SIN60*(A(IA+I)-A(IC+I))) - 
     &             (0.5*(B(IA+I)+B(IC+I))+B(IB+I))
              C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I))) -
     &             (0.5*(B(IA+I)+B(IC+I))+B(IB+I))
              C(JC+J)=SIN60*(B(IC+I)-B(IA+I)) +
     &             (0.5*(A(IA+I)+A(IC+I))-A(IB+I))
              C(JE+J)=SIN60*(B(IC+I)-B(IA+I)) - 
     &             (0.5*(A(IA+I)+A(IC+I))-A(IB+I))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        else
          do L=1,LA
            I=IBASE
            J=JBASE
            do IJK=1,LOT
              C(JA+J)=(2._r8*(A(IA+I)+A(ID+I)))
     &              +(2._r8*(A(IB+I)+A(IC+I)))
              C(JD+J)=(2._r8*(A(IA+I)-A(ID+I)))
     &             -(2._r8*(A(IB+I)-A(IC+I)))
              C(JB+J)=(2._r8*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &             -(SSIN60*(B(IB+I)+B(IC+I)))
              C(JF+J)=(2._r8*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &             +(SSIN60*(B(IB+I)+B(IC+I)))
              C(JC+J)=(2._r8*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     &             -(SSIN60*(B(IB+I)-B(IC+I)))
              C(JE+J)=(2._r8*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     &             +(SSIN60*(B(IB+I)-B(IC+I)))
              I=I+INC3
              J=J+INC4
            enddo
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
          enddo
        endif

c---IFAC=8
      elseif (IFAC.eq.8) then

          if (.not.LAEQM) then
           IERR = 3
            return
          endif

         IHLP= int(2._r8*FLA*FINC1)
         IA  = 1
         IB  = IA+LA*INC1
         IC  = IB+IHLP
         ID  = IC+IHLP
         IE  = ID+IHLP
         JA  = 1
         JB  = JA+JINK
         JC  = JB+JINK
         JD  = JC+JINK
         JE  = JD+JINK
         JF  = JE+JINK
         JG  = JF+JINK
         JH  = JG+JINK
       do L=1,LA
         I=IBASE
         J=JBASE
        do IJK=1,LOT
          C(JA+J)=2._r8*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JE+J)=2._r8*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JC+J)=2._r8*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
          C(JG+J)=2._r8*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
          C(JB+J)=2._r8*((A(IA+I)-A(IE+I))-B(IC+I))
     &         +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JF+J)=2._r8*((A(IA+I)-A(IE+I))-B(IC+I))
     &         -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JD+J)=2._r8*((A(IA+I)-A(IE+I))+B(IC+I))
     &         -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          C(JH+J)=2._r8*((A(IA+I)-A(IE+I))+B(IC+I))
     &         +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          I=I+INC3
          J=J+INC4
        enddo
         IBASE = IBASE + INC1
         JBASE = JBASE + INC2
       enddo

c---IFAC failure
      else
        IERR = 2
        return

      endif

      return
      end




c-----------------------------------------------------------------------
      subroutine FFT_DDSS (NT,A0,NMMAX,DD,SS)
c-----------------------------------------------------------------------
c Calculate the constant arrays needed for evaluating spectral coeffisients 
c     for u and v when the voritcity and divergence is given.
c
c   DD(N,M) = -A0 * sqrt( (N*N-M*M)/(4*N*N-1) ) / N  "divergence" factors  
c   SS(N,M) = -A0 * M / ( N*(N+1) )                  "vorticity" factors
c
c  Revised:  J. K. Sundet, January 1994 (Method from ECMWF manual)
c----------------------------------------------------------------------
      use cmn_precision, only: r8, r4
      implicit none
      integer, intent(in)                  :: NT, NMMAX
      real(r8), intent(in)                   :: A0
      real(r8), dimension(NMMAX),intent(out) :: DD, SS

      real(r8) :: ZM,ZN
      integer :: NTP1, IMD,IMN, JM,JN
c-----------------------------------------------------------------------
       NTP1 = NT+1
       IMD  = 1
       IMN  = 1
      do JM=1,NTP1
          ZM = JM - 1._r8
        do JN=JM,NTP1
           ZN = JN - 1._r8
          if (JN.GT.1) then
            DD(IMD) = -sqrt((ZN*ZN - ZM*ZM)/(4._r8*ZN*ZN - 1._r8))/ZN*A0
            SS(IMN) = -ZM*A0/(ZN*(ZN + 1._r8))
          else
            DD(IMD) = 0._r8
            SS(IMN) = 0._r8
          endif
            IMN = IMN + 1
            IMD = IMD + 1
        enddo
          ZN  = ZN + 1._r8
          DD(IMD) = -sqrt((ZN*ZN - ZM*ZM)/(4._r8*ZN*ZN - 1._r8))/ZN*A0
          IMD = IMD + 1
      enddo
      return
      end


c------------------------------------------------------------------------
      subroutine FFT_SET (TRIGS,IFAX,N)
c------------------------------------------------------------------------
c---computes factors of n & trigonometric functions required by FFT991
c
c    Adapted to CTM, March 1996: J.K. Sundet
C------------------------------------------------------------------------
      use cmn_precision, only: r8
      use utilities, only: ctmExitIJL
      implicit none
      integer,intent(in)                  :: N
      integer,dimension(10),intent(out)   :: IFAX
      real(r8), dimension(N), intent(out)   :: TRIGS

      real(r8) :: ANGLE,DEL,FN
      integer :: JFAX(10),  I,J,K,NHL,NU,  IFAC,NFAX,II1,II2
c---LFAX are the factorization for FFTs, notes say the order is 8,6,5,4,3,2?
      integer, parameter, dimension(7) :: LFAX = [6,8,5,4,3,2,1]

      real(r8), parameter ::  CPI = 3.141592653589793_r8

C-----------------------------------------------------------------------
        FN  = real(N, r8)
        NHL = int(FN*0.5_r8 - 1._r8)
        DEL = 2._r8*CPI/FN
      do K=0,NHL
        ANGLE        = real(K, r8)*DEL
        TRIGS(2*K+1) = cos(ANGLE)
        TRIGS(2*K+2) = sin(ANGLE)
      enddo

c---this section of factorization is a mess, only updated exit calls
C Find factors of N (8,6,5,4,3,2; only one 8 allowed)
         NU = N
         K  = 0
         I  = 1
      do 100 II1=1,999
        if (NU.GT.1) then
          IFAC = LFAX(I)
          if (mod(NU,IFAC).eq.0) then
            K = K + 1
            JFAX(K) = IFAC
            if (IFAC.ne.8) then
C Apart from 8, there might be other factors that appears more that once
              do II2=1,999
                NU = NU/IFAC
                if (mod(NU,IFAC).eq.0) then
                  K = K + 1
                  JFAX(K) = IFAC
                else
                  goto 100
                endif
              enddo
            else
C Eight is encountered
              if(K.GT.1) then
                JFAX(1) = 8
                JFAX(K) = 6
              endif
              NU = NU/IFAC
            endif
          endif
C There are only six factors available
          I = I + 1
          if (I.gt.7) then
          call ctmEXITIJL('>>>>FFT_SET: N has illegal factors',N,N,N)
          endif
          goto 100
        elseif (NU.eq.1) then
C Reverse order of factors
          NFAX    = K
          IFAX(1) = NFAX
          do J=1,NFAX
            IFAX(NFAX+2-J) = JFAX(J)
          enddo
          IFAX(10) = N
          goto 99
C NFAC.LT.1
        else
          call ctmEXITIJL('>>>>FFT_SET: N has illegal factors',N,N,N)
        endif

  100 continue
   99 continue

      return
      end


c----------------------------------------------------------------------
      subroutine LEGGEN (PLEG,PLAT,NT,NMMAX)
c----------------------------------------------------------------------
c   generate Legendre functions, rewritten from SPLEG1 (ECMWF model) 
c       NT = CTM truncation
c       NMMAX = max number of spectral coeffs = (T+1)*(T+4)/2
c       PLAT  = Latitude in radians
c
c   Adapted to CTM:    J. K. Sundet   July 1994
c----------------------------------------------------------------------
      use cmn_precision, only: r8
      use utilities, only: ctmExitC
      implicit none
      integer,intent(in)      ::  NT,NMMAX
      real(r8), intent(in)      ::  PLAT
      real(r8), dimension(NMMAX), intent(out) :: PLEG

c  for T159, NMMAX=(159+1)*(159+4)/2 = 13040, for T319, = 51680
      real(r8), dimension(51680) :: ZHLP1, ZHLP2, ZHLP3
      real(r8) ZSIN,ZCOS,ZF1M,ZRE1,ZF2M,ZM,ZN,ZE1,ZE2
      integer NT1,J1M,JLM,JM,JCN,NPC
     
      if (NMMAX.gt.51680) then
        call ctmExitC('>>>>LEGGEN not enough scratch space T>319')
      endif

        NT1 = NT + 1
        ZSIN   = sin(PLAT)
        ZCOS   = sqrt(1._r8 - ZSIN*ZSIN)
        JLM    = 2
        PLEG(1)= 1._r8
        ZF1M   = sqrt(3._r8)
        PLEG(2)= ZF1M*ZSIN
        NPC    = 2
      do JM = 2,NT1
        ZM       = JM - 1._r8
        ZHLP1(JM)= sqrt(2._r8*ZM + 3._r8)
        ZHLP2(JM)= 1._r8/sqrt(2._r8*ZM)
        NPC      = NPC + 1
      enddo
        ZHLP1(1) = sqrt(3._r8)

      do JM = 1,NT1
        J1M = JM - 1
        ZM  = J1M
        ZRE1= ZHLP1(JM)
        ZE1 = 1._r8/ZRE1
       if (J1M.ne.0) then
          ZF2M = ZF1M*ZCOS*ZHLP2(JM)
          ZF1M = ZF2M*ZRE1
          JLM  = JLM + 1
          PLEG(JLM)= ZF2M
          JLM  = JLM + 1
          PLEG(JLM)= ZF1M*ZSIN
          NPC  = NPC + 1
       endif
       do JCN = J1M+2,NT1
         ZN   = JCN
         ZHLP3(JCN)= sqrt((4._r8*ZN*ZN - 1._r8)/(ZN*ZN - ZM*ZM))
         NPC  = NPC + 1
       enddo
       do JCN = J1M+2,NT1
         ZE2  = ZHLP3(JCN)
         JLM  = JLM + 1
         PLEG(JLM)= ZE2*( ZSIN*PLEG(JLM-1) - ZE1*PLEG(JLM-2) )
         ZE2  = 1._r8/ZE2
         ZE1  = ZE2
         NPC  = NPC + 1
       enddo
      enddo

      return
      end




c---------------------------------------------------------------------
      subroutine GAUSST(NG,XP,WT)
c---------------------------------------------------------------------
c---calculate NG Gaussian quadrature points (XP) & weights (WT)
c---     for interval (X1,X2)  from GISS-Lacis code.
c---tested against EC version GAUAW, simpler, same to 1.d-13 or better)
      use cmn_precision, only: r8

      implicit none
      integer, intent(in) ::  NG
      real(r8), dimension(NG), intent(out) :: XP, WT

      real(r8), parameter :: PI = 3.141592653589793_r8
      real(r8), parameter :: PS = 1.013211836423378e-01_r8
      real(r8), parameter :: DXL = 1.e-16_r8
      real(r8), parameter :: X1 = -1._r8,  X2 = 1._r8

      real(r8) XMID,XHAF,DNG,DN,DM,DI,C,Z,ZZ,PN,PNI,PNJ,PNK,DX,X
      integer NN,N2,N,I,J

      stop 'Rather use grid.f90: GAUSST2'

        XMID = (X2+X1)/2._r8
        XHAF = (X2-X1)/2._r8
        DNG = NG
        NN = NG/2
        N2 = NN*2
      if (N2.eq.NG) goto 110

        XP(NN+1) = XMID
        WT(NN+1) = 1._r8
      if (NG.lt.2) return

        PN = 1._r8
        N = 0
  100 N = N+2
        DN = N
        DM = DN-1._r8
        PN = PN*(DM/DN)
      if (N.lt.N2) goto 100

        WT(NN+1) = 2._r8*XHAF/(DNG*PN)**2
  110 I = 0
        C = PI/sqrt(DNG*(DNG+1._r8)+0.5_r8-PS)/105._r8
  120 I = I+1
        DI = I
        Z = PS/(4._r8*DI-1._r8)**2
        ZZ = (105._r8+Z*(210._r8-Z*(2170._r8-Z
     &                *(105812._r8-12554474._r8*Z))))
        X = DCOS(ZZ*C*(DI-0.25_r8))
  130 N = 1
        DM = 1._r8
        PNI = 1._r8
        PNJ = X
  140 N = N+1
        DN = N
        PNK = ((DM+DN)*X*PNJ-DM*PNI)/DN
        PNI = PNJ
        PNJ = PNK
        DM = DN
      if (N.lt.NG) goto 140

        DX = PNJ*(1._r8-X*X)/DNG/(PNI-X*PNJ)
        X = X-DX
      if (abs(DX).gt.DXL) goto 130

        J = NG+1-I
        XP(I) = XMID-XHAF*X
        XP(J) = XMID+XHAF*X
        WT(I) = 2._r8*XHAF*(1._r8-X*X)/(DNG*PNI)**2
        WT(J) = WT(I)
      if (I.lt.NN) goto 120

      return
      end


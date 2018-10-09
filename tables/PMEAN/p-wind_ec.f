c-----------------------------------------------------------------------
c---(pwind_ec.f)----generic CTM shell from UCIrvine (p-code 5.5, 1/2008)


c------ PWIND_E19 reads met data from ECMWF Integrated Fcst System
c---subroutines:  WIND,SPE2GP,ZD2UV,UVCOEF,TRUNGR,ZEROW


c      subroutine FFT_99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT)
c      subroutine FFT_RPASS(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,
c                            IOC,LA,IERR)
c      subroutine CFRMIN(T,CFR,CW,CI,EPS)

c-----------------------------------------------------------------------
      Subroutine WIND (ILOOP,NDAY,NMET)
C-----------------------------------------------------------------------
C Read input data. ECMWF data has three files, one for spectral 
C      data (.b01), one for 2-d gridpoint data (.b03) and one for 
C      3-d gridpoint data (.b02).
C Spectral transformation is done in this routine, and sumation of
C      physical fields from T63 ---> T42, T21 or other resolution.
C Allow several different formats for 2-d gridpoint files by labelling
C      the record numbers we need (we need N15R records)
C
C-----------------------------------------------------------------------
C      Implicit NONE is in cmn_h.f
C-----------------------------------------------------------------------
      Include 'cmn_h.f'
      Include 'cmn_w.f'
      Include 'cmn_s.f'
C-----------------------------------------------------------------------
      Integer, Intent(in) ::  ILOOP, NDAY, NMET
C-----------------------------------------------------------------------
C I/O parameters:    ILOOP      ! <=0, init, > 0 ordinary run
c                    NMET       ! the dynamical timestep
c                    NDAY       ! day number, should match TMET
C-----------------------------------------------------------------------
c Resolution degradation
      Integer   N15R
C
      Integer   NINTW
C     NINTW should not be a parameter, since we use this routine for
C     ERA-40 data as well as IFS, where it should be 6. Test for the
C     ERA-40 is done below according to logical ERA, and NINTW is set to
C     3 or 6 for IFS and ERA-40, respectively.
      Logical, Parameter ::  ERA=.false., VERBOSE=.false.

      Real*8, Parameter  ::  EPS=0.01d0 ! minimum cloud fraction

C Local parameters
      Real*8, Allocatable ::
     &        FFTAR(:,:), WORK(:), VTMP(:,:), UTMP(:), PW(:,:)
     &       ,TW(:,:,:), CLDFRW(:,:,:), CLDIWCW(:,:,:), CLDLWCW(:,:,:)
     &       ,ZOFLEW(:,:,:),CDETU(:,:,:),CDETD(:,:,:)
      Real*8
     &     VEDGE,BAND,DETA,DETB,DELP,
     &     ZCOS,DT,ZDT,SFTD,LV,ESAT,SDEN,
     &     MZERO,MZICE,RZERO,RZICE,PT,PB,
     &     DELZ,LWC,IWC,ODCW,ODCI,TFACT,QMIN,
     &     DMASS,ZAREA,PSRF,
     &     POFLE(LPAR+1)
      Real*4, Dimension(:,:), Allocatable :: 
     &        GRDAT, EWSS, NSSS
      Real*4, Allocatable ::
     &        GRDATA(:,:,:)
      Real*4, Dimension(:,:), Allocatable :: 
     &        GRDATHI, SPWK5, SPWK6
      Real*4, Allocatable ::
     &        GRDATAHI(:,:,:), PSPEHI(:)
      Real*4  RSUM
      Real*8, Dimension(:,:,:), Allocatable :: 
     &        WK3DT, WK3DU, WK3DV
      Integer
     &     LMAP(LPAR+1)
      Integer
     &     ISEC(16),
     &     NREAD,NSTAGE,
     &     I,J,II,JJ,L,LL,
     &     FLD_ID, IRAN0
      Logical  LUV
Cmga increased characters to be more flexible      Character*15
      Character*80
     &     FNAME1,FNAME2,FNAME3

      Integer NF, NrOfFields
C To check if RAIN is split into CONVRAIN/LSRAIN or not
      Logical found216 ! if field 216 is found we have only PRECLS
C If 217 and 218 are found instead of 216 we have CONVRAIN/LSRAIN:
      Real*8, Dimension(:,:,:), Allocatable ::
     &     CONVRAIN, LSRAIN
Cmga--^

      Save   LMAP,N15R

      Allocate (
     &          FFTAR(FFTPTS,NRFFT), WORK(FFTPTS*NRFFT),
     &          VTMP(IPARW,JPARW), UTMP(IPARW),
     &          WK3DU(IPARW,JPARW,LPARW), WK3DV(IPARW,JPARW,LPARW),
     &          WK3DT(IPARW,JPARW,LPARW),
     &          ZOFLEW(LPAR+1,IPARW,JPARW),
     &          PW(IPARW,JPARW), TW(IPARW,JPARW,LPAR),
     &          CLDFRW(IPARW,JPARW,LWEPAR), CLDIWCW(IPARW,JPARW,LWEPAR),
     &          CLDLWCW(IPARW,JPARW,LWEPAR),
     &          CDETU(IPAR,JPAR,LWEPAR),CDETD(IPAR,JPAR,LWDPAR),
     &          CONVRAIN(IPAR,JPAR,LPARW), LSRAIN(IPAR,JPAR,LPARW) )
      Allocate (
     &          GRDAT(IPAR,JPAR), GRDATA(IPAR,JPAR,LPARW),
     &          EWSS(IPAR,JPAR), NSSS(IPAR,JPAR) )
      Allocate (
     &          GRDATHI(IPARW,JPARW), GRDATAHI(IPARW,JPARW,LPARW),
     &          SPWK5(NRNM,LPARW), SPWK6(NRNM,LPARW),
     &          PSPEHI(NRNM) )
c
c  Record numbers of desired 2-d grid point data, corresponding to
c               SSHF,SLHF,U10m,V10m,T2m,TDw2m,EWSS,NSSS,BLH, SA, MSLP
C-----------------------------------------------------------------------
C NEED
C Should check that the date contained in ISEC is OK compared
C      to the model date ...
C NEED

Cmga --v
      found216=.false.
Cmga --^

      If (.not. ERA) Then
        NINTW = 3 ! IFS: 3 hours per meteorological timestep
        NREAD = 24/NRMETD
      Else 
        NINTW = 6 ! ERA-40: 6 hours timestep
        NREAD = 24/NRMETD*(NINTW/3)
      Endif

C ERROR check
c     NREAD  = 24/NRMETD
      NSTAGE = NREAD/NINTW
      DT     = 3600.D0 * dble(NINTW)
      ZDT    = 1.d0 / DT
      If (NSTAGE*NINTW.NE.NREAD)
     +    Call EXITC('INCONSISTENT TIMING OF DATA')

C Initialize
      LUV  = .FALSE.
C Define minimum humidity QMIN in kg/kg
      QMIN = 3.d-6*18.d0/29.d0

c---locate the position of random number sequence based on year/day/hour
      IRAN0 = 1 + 3*(JYEAR-2000) + 7*JDAY + 11*nint(GMTAU)

C Initial step - setup dry airmass (P() and Q())
      If (ILOOP .le. 0)  then

        If (Mod(NREAD,NINTW).NE.0)
     +      Call EXITC('INCONSISTENT TIME INTERVAL')
        If (IM.NE.IPAR .OR. JM.NE.JPAR .OR. LM.NE.LPAR)
     +      Call EXITC('INCONSISTENT DIMENSIONS')

C Set up level weightings if vertical resolution degraded
        do LL = LPARW+1,1,-1
          L   = LMMAP(LL)
          LMAP(L) = LL
        enddo
c
        FNAME1 = ''
        FNAME2 = ''
        FNAME3 = ''
        write (MPATH2(1:4),'(i4)') MYEAR
        Write (MFILE3(3:5),'(A3)')  TMET
        FNAME1 = trim(MPATH1)//trim(MPATH2)//trim(MFILE3)//'.b01'

C print out info, open files and read in data, transform and sum up
        Write(6,1000) FNAME1,JMON,NDAY
        Open (13,FILE=trim(FNAME1),FORM='UNFORMATTED'
     +           ,STATUS='OLD',ERR=94)
C Pressure field
        Read(13) FLD_ID
        If (FLD_ID.NE.152) Call EXITC('ERROR, input Data PSPE')
          Read(13) PSPEHI

C Pressure; transform from spectral to grid point data
          Call SPE2GP(PSPEHI,NRNM,LUV,PW,IPARW,JPARW,1, ALP,NMMAX,
     &                FFTAR,WORK,TRIG,IFAX,(NTRUNW+1),FFTPTS,NRFFT)

C save high-res pressure Log(Ps) --> hPa
          Do J = 1,JPARW
            Do I = 1,IPARW
              PW(I,J) = Exp(PW(I,J))*1.d-2
            Enddo
          Enddo
        if (ldeg) then
          Call TRUNG8(PW,P,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,
     &                JDGRD,IPARW,JPARW,IM,JM,1,1)
        else
          Do J = 1,JM
          Do I = 1,IM
            P(I,J) = PW(I,J)
          Enddo
          Enddo
        endif
         
C Initialization is done - these steps are for all time steps
      Else

C Clear all key arrays
        Call ZEROW(U,V,T,Q,CDRY,CWETN,CWETE,CWETD,CENTU,CENTD,
     &             CDETU,CDETD,
     &             PRECMC,PRECLS,ODW,ODI,OD,CLDFR,CNVCF,STRCF,CLDLWC,
     &             CLDIWC, SLH,SHF,SMF,SFT,SFQ,SFU,SFV,BLH,
     &             MSLP,SA,IPAR,JPAR,LPAR,LDPAR,LWEPAR,LWDPAR)

C New day? If so, open new files
        If (NMET.EQ.1) Then

          write (MPATH2(1:4),'(i4)') MYEAR
          Write (MFILE3(3:5),'(A3)')  TMET
          FNAME1 = trim(MPATH1)//trim(MPATH2)//trim(MFILE3)//'.b01'
          Write(6,1000) FNAME1,JMON,NDAY

          Close(13)
          Open (13,FILE=trim(FNAME1),FORM='UNFORMATTED'
     +           ,STATUS='OLD',ERR=94)

        Endif
C
C -- SPECTRAL DATA --
C Read in spectral data - error check only on Pressure since it 
C  is assumed that the ordering of fields is fine
C
C Pressure, Temperature, Vorticity, Divergence, (Vertical Velocity ignored)
        Read(13) FLD_ID
        If (FLD_ID.NE.152) Call EXITC('ERROR, input Data PSPE')

        Read(13) PSPEHI
C Pressure; transform from spectral to grid point data
        Call SPE2GP(PSPEHI,NRNM,LUV,PW,IPARW,JPARW,1,ALP,
     &           NMMAX,FFTAR,WORK,TRIG,IFAX,(NTRUNW+1),FFTPTS,NRFFT)
C save high-res pressure Log(Ps) --> hPa
        Do J = 1,JPARW
          Do I = 1,IPARW
            PW(I,J) = Exp(PW(I,J))*1.d-2
          Enddo
        Enddo

        If (ldeg) then
          Call TRUNG8(PW,P,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,
     &                JDGRD,IPARW,JPARW,IM,JM,1,1)
        Else
          Do J = 1,JM
          Do I = 1,IM
            P(I,J) = PW(I,J)
          Enddo
          Enddo
        Endif

C Temperature
        Read(13) FLD_ID
        If (FLD_ID.NE.130) Call EXITC('ERROR, input Data T')
        Read(13) SPWK5
C Pressure; transform from spectral to grid point data
        Call SPE2GP(SPWK5,NRNM,LUV,WK3DT,IPARW,JPARW,LPARW,ALP,
     &           NMMAX,FFTAR,WORK,TRIG,IFAX,(NTRUNW+1),FFTPTS,NRFFT)
        TW(:,:,:)  = 0.d0
        Do L = 1,LM
          Do LL = LMAP(L),LMAP(L+1)-1
            Do J = 1,JPARW
              Do I = 1,IPARW
                TW(I,J,L) = TW(I,J,L) + WK3DT(I,J,LL)*XLMMAP(LL)
              End Do
            End Do
          End Do
        End Do
        if (ldeg) then
          Call TRUNG8(TW,T,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,
     &                JDGRD,IPARW,JPARW,IM,JM,LM,1)
        else
          Do L = 1,LM
            Do J = 1,JM
              Do I = 1,IM
                T(I,J,L) = TW(I,J,L)
              End Do
            End Do
          End Do
        endif

C Vorticity and Divergence - find U and V
        Read(13) FLD_ID
        If (FLD_ID.NE.138) Call EXITC('ERROR, input Data U')
        Read(13) SPWK5

        Read(13) FLD_ID
        If (FLD_ID.NE.155) Call EXITC('ERROR, input Data V')
        Read(13) SPWK6

C       VERTICAL VELOCITY [Pa/s] - PHONY READ IN, NOT FOR YEAR 2001-->
C-----------------------------------------------------------------------
Cmga:
C W is included in Jostein's old 1997/2000/2001 data, but not in the
C new 2001-2006 data sets. To make this routine compatible with both
C data sets read first the FLD_ID. If it is 135 (as in the old data)
C then print 'Phoney read 3D-field' and read RSUM. In the other case,
C i.e. FLD_ID.ne.135, we have read the FLD_ID for the next field
C and need to jump back one line: Backspace(13)!
C-----------------------------------------------------------------------
        Read(13,END=100) FLD_ID
        If (FLD_ID.eq.135) then
c          print*,'PWIND: Phoney read in W:',FLD_ID
          Read(13) RSUM
        Else
c          print*,'PWIND: No W in this data set! FLD_ID now:',FLD_ID
c          print*,'PWIND: Backing up one line.'
          Backspace(13)
        End If
 100    Continue

      Endif

      Deallocate (
     &          FFTAR, WORK, VTMP, UTMP, GRDAT, GRDATA, EWSS, NSSS,
     &          GRDATHI, GRDATAHI, SPWK5, SPWK6,
     &          PSPEHI, WK3DU, WK3DV, WK3DT, ZOFLEW, PW, TW,
     &          CLDFRW, CLDIWCW,CLDLWCW,
     &          CONVRAIN, LSRAIN, CDETU,CDETD )

      Return
 1000 Format(1X,A50,4I5)
C      
 91   Call EXITC('============INCONSISTENT TIME INTERVAL==============')
 93   Call EXITC('============INCONSISTENT TIMING OF DATA=============')
 94   Call EXITC('================Error on .b01 file==================')
 95   Call EXITC('================Error on .b02 file==================')
 96   Call EXITC('================Error on .b03 file==================')
C      
      End
C
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
      Implicit NONE
C-----------------------------------------------------------------------
C I/O parameters
      Integer
     &     IFAX(10),
     &     NMMAX,NRNM,NTP1,
     &     FFTPTS,NRFFT,
     &     ID,JD,LD
      Real*8
     &     GPT(ID,JD,LD),       ! gaussian gridpt values or fluxes
     &     FFTAR(FFTPTS,NRFFT), ! fourier to gridpoint transf. array
     &     WORK(FFTPTS*NRFFT),  ! work area for FFT
     &     ALPGAU(NMMAX,JD/2),
     &     TRIG(ID)
      Real*4
     &     SPE(2,NRNM/2,LD)     ! spectral coeffs
      Logical
     &     LUV                  ! .TRUE. if U or V
C Local parameters
      Real*8
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
        Call EXITC('=================SP2GP================')
      Endif

C Do the transform in two operations, multiply the spec coeffs
C with appropriate functions. Then sum the ones (odd and even) to
C find the symetric and anti-symetric functions (rel. to equator).

      Do 109 L=1,LD

C For all pairs of NRFFT latitudes
        Do 105 IGRP=1,NGROUP
          Do J=1,NRFFT
            Do I=1,FFTPTS
              FFTAR(I,J) = 0.d0
              WORK(I + (J-1)*FFTPTS) = 0.d0
            End Do
          End Do
          If (IGRP.EQ.1) Then
            JBACK = 0
          Else
            JBACK = ((IGRP-1)*NRFFT)/2
          Endif
C Real*4 --> Real*8
          Do I=1,NRNM/2
            SPEWRK(1,I) = SPE(1,I,L)
            SPEWRK(2,I) = SPE(2,I,L)
          End Do
          Do J=1,NRFFT/2
            IMP = 0
            IMN = 0
C For each zonal wavenumber sum the coeffs*A-Leg coeffs.
            Do ILONG=1,IMLIM
              SUMRES = 0.d0
              SUMIMS = 0.d0
              SUMREA = 0.d0
              SUMIMA = 0.d0
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
c
c
C---------------------------------------------------------------------
      Subroutine TRUNG8(SPHI,SPLO,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,JDGRD,
     &                  IPARW,JPARW,IM,JM,LM,LV )
C---------------------------------------------------------------------
C Average gridpoint data to go with truncation of spectral fields
C---------------------------------------------------------------------
      Implicit NONE
C---------------------------------------------------------------------
C I/O Variables
      Integer IPARW, JPARW, IM, JM, LM, LV
      Integer IDGRD,JDGRD,IMAP(IDGRD,IM),JMAP(JDGRD,JM)
      Real*8  SPHI(IPARW,JPARW,LM,LV),SPLO(IM,JM,LM,LV)
      Real*8  ZDEGI(IDGRD,IM),ZDEGJ(JDGRD,IM)

c Local Variables
      Integer I,J,K,L,ix,jx
c
c  Zero output array
      SPLO(:,:,:,:) = 0.0d0
c
      do K=1,LV
        do L=1,LM
          do J=1,JM
            do I=1,IM
              do jx=1,jdgrd
                do ix=1,idgrd
                  SPLO(I,J,L,K) = SPLO(I,J,L,K) +
     &                                 SPHI(imap(ix,i),jmap(jx,j),L,K)
     &                                          *zdegi(ix,i)*zdegj(jx,j)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      RETURN
      END
c
c
C---------------------------------------------------------------------
      Subroutine TRUNGR(SPHI,SPLO,ZDEGI,ZDEGJ,IMAP,JMAP,IDGRD,JDGRD,
     &                  IPARW,JPARW,IM,JM,LM,LV )
C---------------------------------------------------------------------
C Average gridpoint data to go with truncation of spectral fields
C---------------------------------------------------------------------
      Implicit NONE
C---------------------------------------------------------------------
C I/O Variables
      Integer IPARW, JPARW, IM, JM, LM, LV
      Integer IDGRD,JDGRD,IMAP(IDGRD,IM),JMAP(JDGRD,JM)
      Real*4  SPHI(IPARW,JPARW,LM,LV),SPLO(IM,JM,LM,LV)
      Real*8  ZDEGI(IDGRD,IM),ZDEGJ(JDGRD,IM)
c Local Variables
      Integer I,J,K,L,ix,jx
c
c  Zero output array
      SPLO(:,:,:,:) = 0.0d0
c
      do K=1,LV
        do L=1,LM
          do J=1,JM
            do I=1,IM
              do jx=1,jdgrd
                do ix=1,idgrd
                  SPLO(I,J,L,K) = SPLO(I,J,L,K) +
     &                                 SPHI(imap(ix,i),jmap(jx,j),L,K)
     &                                          *zdegi(ix,i)*zdegj(jx,j)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      RETURN
      END


c-----------------------------------------------------------------------
      Subroutine ZEROW(U,V,T,Q,CDRY,CWETN,CWETE,CWETD,CENTU,CENTD,
     &                CDETU,CDETD,
     &                PRECMC,PRECLS,ODW,ODI,OD,CLDFR,CNVCF,STRCF,CLDLWC,
     &                CLDIWC, SLH,SHF,SMF,SFT,SFQ,SFU,SFV,BLH,
     &                MSLP,SA,IM,JM,LM,LD,LWE,LWD)
c-----------------------------------------------------------------------
c  Old routine from pwind_g9 to zero all met field variables
c-----------------------------------------------------------------------
      implicit none
      integer  IM,JM,LM,LD,LWE,LWD,K
      real*8   U(*),V(*),T(*),Q(*)
      real*8   CDRY(*),CWETN(*),CWETE(*),CWETD(*),CENTU(*),CENTD(*)
      real*8   CDETU(*),CDETD(*)
      real*8   PRECMC(*),PRECLS(*),ODW(*),ODI(*),OD(*),CLDFR(*)
      real*8   CNVCF(*),STRCF(*),CLDLWC(*),CLDIWC(*)
      real*8   SLH(*),SHF(*),SMF(*),SFT(*),SFQ(*),SFU(*),SFV(*),BLH(*)
      real*8   MSLP(*),SA(*)
c
      do K = 1,IM*JM*LM
        U(K)     = 0.d0
        V(K)     = 0.d0
        Q(K)     = 0.d0 
        T(K)     = 0.d0
      enddo
      do K = 1,IM*JM*LD
        CDRY(K)  = 0.d0
      enddo
      do K = 1,IM*JM*LWE
        CWETN(K) = 0.d0
        CWETE(K) = 0.d0
        CENTU(K) = 0.d0
        CDETU(K) = 0.d0
        PRECMC(K)= 0.d0
        PRECLS(K)= 0.d0
        ODW(K)   = 0.d0
        ODI(K)   = 0.d0
        OD(K)    = 0.d0
        CLDFR(K) = 0.d0
        CNVCF(K) = 0.d0
        STRCF(K )= 0.d0
        CLDLWC(K)= 0.d0
        CLDIWC(K)= 0.d0
      enddo
      do K = 1,IM*JM*LWD
        CWETD(K) = 0.d0
        CENTD(K) = 0.d0
        CDETD(K) = 0.d0
      enddo
      do K = 1,IM*JM
        SLH(K)   = 0.d0
        SHF(K)   = 0.d0
        SMF(K)   = 0.d0
        SFT(K)   = 0.d0
        SFQ(K)   = 0.d0
        SFU(K)   = 0.d0
        SFV(K)   = 0.d0
        BLH(K)   = 0.d0
        MSLP(K)  = 0.d0
        SA(K)    = 0.d0
      enddo
c
      return
      end


C-----------------------------------------------------------------------
C
C    Summary of EC fields
C
C    *.b01  Spectral Data read on unit 13 - one record per met record
C
C            PSPE     Surface Pressure        (nrnm)
C            Temp     Temperature             (nrnm,lpar)
C            U        Vorticity               (nrnm,lpar)
C            V        Divergence              (nrnm,lpar)
C            W        Vertical Velocity       (nrnm,lpar)   [not used]
C
C    *.b02  Gridpoint Data read on unit 14 - one record per day (recl varies!)
C
C            Q        Water Vapour                       (ipar,jpar,lpar)
C            UpN      Mass Flux updrafts [Kg/(m^2*s)]    (ipar,jpar,*)
C            DnN      Mass Flux downdrafts               (ipar,jpar,*)
C            UpE      Mass Flux updrafts entrainment     (ipar,jpar,*)
C            DnE      Mass Flux downdrafts entrainment   (ipar,jpar,*)
C            Rain     Rainfall                           (ipar,jpar,*)
C            CLW      Cloud Liquid Water Content [Kg/Kg] (ipar,jpar,*)
C            CIWC     Cloud Ice Water Content [Kg/Kg]    (ipar,jpar,*)
C            CLDFR    Cloud Fraction [0,1]               (ipar,jpar,*)
C
C    *.b03  Gridpoint Data read on unit 15 - one record per variable
C
C            TCW      Total Column Water             (ipar,jpar)   [not used]
C            LSP      Large Scale Precip.            (ipar,jpar)   [not used]
C            CP       Convective Precip.             (ipar,jpar)   [not used]
C          [ BLD      Boundary Layer Dissipation     (ipar,jpar)  T63 Only ]
C            SSHF     Surface Sensible Heat Flux     (ipar,jpar)
C            SLHF     Surface Latent Heat FLux SLHF  (ipar,jpar)
C            U10M     10m U wind component           (ipar,jpar)
C            V10M     10m V wind component           (ipar,jpar)
C            T2M      2m Temperature                 (ipar,jpar)
C            TDEW2M   2m Dewpoint Temperature        (ipar,jpar)
C            EWSS     E/W Surface Stress             (ipar,jpar)
C            NSSS     N/S Surface Stress             (ipar,jpar)
C          [ RUNOFF   Runoff                         (ipar,jpar)  T63 Only ]
C            CLDTYP   Accumulated cloud type         (ipar,jpar)   [not used]
C            BLH      Boundary Layer Height          (ipar,jpar)
C
C-----------------------------------------------------------------------


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
      implicit none
      integer, intent(inout)      ::  INC,JUMP,N,LOT
      integer, dimension(10), intent(inout)      :: IFAX
      real*8, dimension(JUMP*LOT), intent(inout) :: A,WORK
      real*8, dimension(N), intent(inout)        :: TRIGS

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
          A(I+INC) = 0.5d0*A(I)
          I        = I + JUMP
        enddo
        if (.not.LNODD) then
            I = ISTART + N*INC
          do J=1,NVEX
            A(I) = 0.5d0*A(I)
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
      call EXITIJL('>>>>error FFT_RPASS: IERR,NVEX,IFAC',IERR,NVEX,IFAC)
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
          A(IX)     = 0.d0
          A(IX+INC) = 0.d0
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
      implicit none
      integer, intent(in)  :: INC1,INC2,INC3,INC4,LOT,N,IFAC,IOC,LA
      integer, intent(out) :: IERR
      real*8, dimension(IOC*LOT), intent(in)  :: A,B
      real*8, dimension(IOC*LOT), intent(out) :: C,D
      real*8, dimension(N), intent(in) :: TRIGS

      real*8, dimension(64) :: A10,A11,A20,A21,B10,B11,B20,B21
      real*8  C1,C2,C3,C4,C5,S1,S2,S3,S4,S5, FN,FIFAC,FLA,FINC1,FINC2,FM
      integer  I,J,K,L,M,IJK, KB,KC,KD,KE,KF
      integer  IA,IB,IC,ID,IE,IG,  JA,JB,JC,JD,JE,JF,JG,JH
      integer  IBASE,JBASE,IINK,JINK,JUMP,KSTOP,IHLP
      logical  LAEQM

      real*8,parameter:: SIN36=0.587785252292473137d0, SSIN36=2.d0*SIN36
      real*8,parameter:: SIN45=0.707106781186547524d0, SSIN45=2.d0*SIN45
      real*8,parameter:: SIN60=0.866025403784438597d0, SSIN60=2.d0*SIN60
      real*8,parameter:: SIN72=0.951056516295153531d0, SSIN72=2.d0*SIN72
      real*8,parameter:: QRT5 =0.559016994374947451d0, QQRT5 =2.d0*QRT5

C---------------------------------------------------------------------
      if (LOT.GT.64) then
        IERR = 1
        return
      endif

      IBASE= 0
      JBASE= 0
      IERR = 0
      FN   = real(N)
      FIFAC= real(IFAC)
      FLA  = real(LA)
      FINC1= real(INC1)
      FINC2= real(INC2)
      FM   = FN/FIFAC
      M    = int(FM)
      IINK = int(FLA*FINC1)
      JINK = int(FLA*FINC2)
      JUMP = int((FIFAC-1.D0)*real(JINK))
      KSTOP= int((FN-FIFAC)/(2.D0*FIFAC))
      LAEQM= LA.eq.M

c---defactorize the FFT coeff's, determined by IFAC
c---IFAC=2
      if (IFAC.eq.2) then

        IA = 1
        IB = IA + int((2.D0*FM - FLA)*FINC1)
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
             C(JA+J)=2.D0*(A(IA+I)+A(IB+I))
             C(JB+J)=2.D0*(A(IA+I)-A(IB+I))
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
        IB = IA+int((2.D0*FM-FLA)*FINC1)
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
               C(JB+J)= C1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I))) - 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S1*((B(IA+I)-0.5d0*(B(IB+I)-B(IC+I))) +
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JB+J)= S1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I))) - 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C1*((B(IA+I)-0.5d0*(B(IB+I)-B(IC+I))) + 
     &              (SIN60*(A(IB+I)-A(IC+I))))
               C(JC+J)= C2*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I))) + 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S2*((B(IA+I)-0.5d0*(B(IB+I)-B(IC+I))) -
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JC+J)= S2*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I))) + 
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C2*((B(IA+I)-0.5d0*(B(IB+I)-B(IC+I))) - 
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
              C(JB+J)=(0.5d0*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
              C(JC+J)=-(0.5d0*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
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
            C(JA+J)=2.D0*(A(IA+I)+A(IB+I))
            C(JB+J)=(2.D0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
            C(JC+J)=(2.D0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
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
        IB = IA+int((2.D0*FM-FLA)*FINC1)
        IC = IB+int(2.D0*FM*FINC1)
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
              C(JA+J)=2.D0*((A(IA+I)+A(IC+I))+A(IB+I))
              C(JB+J)=2.D0*((A(IA+I)-A(IC+I))-B(IB+I))
              C(JC+J)=2.D0*((A(IA+I)+A(IC+I))-A(IB+I))
              C(JD+J)=2.D0*((A(IA+I)-A(IC+I))+B(IB+I))
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
        IB = IA+int((2.D0*FM-FLA)*FINC1)
        IC = IB+int(2.D0*FM*FINC1)
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
              C(JA+J)=2.D0*(A(IA+I)+(A(IB+I)+A(IC+I)))
              C(JB+J)=(2.D0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             +QQRT5*(A(IB+I)-A(IC+I))) - 
     &             (SSIN72*B(IB+I)+SSIN36*B(IC+I))
              C(JC+J)=(2.D0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             -QQRT5*(A(IB+I)-A(IC+I))) -
     &             (SSIN36*B(IB+I)-SSIN72*B(IC+I))
              C(JD+J)=(2.D0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &             -QQRT5*(A(IB+I)-A(IC+I))) + 
     &             (SSIN36*B(IB+I)-SSIN72*B(IC+I))
              C(JE+J)=(2.D0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
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

        IHLP= int(2.D0*FM*FINC1)
        IA  = 1
        IB  = IA+int((2.D0*FM-FLA)*FINC1)
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
              C(JA+J)=(2.D0*(A(IA+I)+A(ID+I)))+(2.D0*(A(IB+I)+A(IC+I)))
              C(JD+J)=(2.D0*(A(IA+I)-A(ID+I)))-(2.D0*(A(IB+I)-A(IC+I)))
              C(JB+J)=(2.D0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &             -(SSIN60*(B(IB+I)+B(IC+I)))
              C(JF+J)=(2.D0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &             +(SSIN60*(B(IB+I)+B(IC+I)))
              C(JC+J)=(2.D0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     &             -(SSIN60*(B(IB+I)-B(IC+I)))
              C(JE+J)=(2.D0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
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

         IHLP= int(2.D0*FLA*FINC1)
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
          C(JA+J)=2.D0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JE+J)=2.D0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JC+J)=2.D0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
          C(JG+J)=2.D0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
          C(JB+J)=2.D0*((A(IA+I)-A(IE+I))-B(IC+I))
     &         +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JF+J)=2.D0*((A(IA+I)-A(IE+I))-B(IC+I))
     &         -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JD+J)=2.D0*((A(IA+I)-A(IE+I))+B(IC+I))
     &         -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          C(JH+J)=2.D0*((A(IA+I)-A(IE+I))+B(IC+I))
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


!//--------------------------------------------------------------------------
      subroutine fluxfilter(cdetu,cdetd)
!//--------------------------------------------------------------------------
      !// Filter convective mass fluxes.
      !// Based on filtering in p-wind_ec.f.
      !//
      !// Amund Sovde, March 2010
!//--------------------------------------------------------------------------
      include 'cmn_h.f'
      include 'cmn_w.f'
!//--------------------------------------------------------------------------
      !// Input
      real*8, intent(inout) :: cdetu(ipar,jpar,lwepar)
      real*8, intent(inout) :: cdetd(ipar,jpar,lwdpar)

      !// Local variables
      real*8 :: sden
      integer :: KCTOP(ipar,jpar), KDBOT(ipar,jpar), i,j,l,LL
!//--------------------------------------------------------------------------

c Remove noise and wrong sign convections
c Noise filter for CWETE
c   the min cwete is -2.3e-5 kg/s/m^2, the neg fluxes are often seen as
c  'dipoles' with an equally small positive updraft flux in a neighboring layer
      do l = 1,lwepar
         do j = 1,jm
            do i = 1,im
               sden = cwete(i,j,l) / areaxy(i,j)
               if (sden .lt. 2.3d-5) cwete(i,j,l) = 0.d0
            enddo
         enddo
      enddo

      if (MODEL(1:6) .eq. 'ECT319')  then
      !// Add more filter for boxes where the box above and below have no flux 
        do l = 2,lwepar-1
          do j = 1,jm
            do i = 1,im
              if (cwete(i,j,l-1).eq.0.d0 .and. cwete(i,j,l).gt.0.d0
     &             .and.cwete(i,j,l+1).eq.0.d0) then
                 !// Some boxes have flux, but no flux below. Filter these
                 !// if smaller than 2.d-4
                sden = CWETE(i,j,l) / areaxy(i,j)
                if (sden .lt. 2.d-4) cwete(i,j,l) = 0.d0
              endif
            enddo
          enddo
        enddo
      endif

c  Noise filter for CWETD
c   the max cwetd is 7.1e-6 kg/s/m^2, the + fluxes are often seen as 'dipoles'
c   with an equally small negative downdraft flux in a neighboring layer
      do l = 1,lwdpar
         do j = 1,jm
            do i = 1,im
               sden = cwetd(i,j,l) / areaxy(i,j)
               if (sden .gt. -7.1d-6) cwetd(i,j,l) = 0.d0
            enddo
         enddo
      enddo

      if (MODEL(1:6) .eq. 'ECT319')  then
      !// Add more filter for boxes where the box above and below have no flux 
        do l = 2,lwdpar-1
          do j = 1,jm
            do i = 1,im
              if (cwetd(i,j,l-1).eq.0.d0 .and. cwetd(i,j,l).lt.0.d0
     &             .and.cwetd(i,j,l+1).eq.0.d0) then
                sden = cwetd(i,j,l) / areaxy(i,j)
                if (sden .gt. -1.d-5) cwetd(i,j,l) = 0.d0
              endif
            enddo
          enddo
        enddo
      endif

        KCTOP(:,:) = LWEPAR - 1

      if (MODEL(1:6) .eq. 'ECT319')  then   ! N160 fixes

c Noise filter for CDETU
        do l = 1,lwepar
          do j = 1,jm
            do i = 1,im
              sden = cdetu(i,j,l) / areaxy(i,j)
              if (sden .lt. 2.3d-5) cdetu(i,j,l) = 0.d0
            enddo
          enddo
        enddo

c Noise filter for CDETD
        do l = 1,lwdpar
          do j = 1,jm
            do i = 1,im
               sden = cdetd(i,j,l) / areaxy(i,j)
               if (sden .lt. 7.1d-6) cdetd(i,j,l) = 0.d0
            enddo
          enddo
        enddo


C more filter for updrafts/downdrafts

      !// Find top of convection
        KCTOP(:,:) = 1
        do j = 1,jm
          do i = 1,im
            do l = lwepar-1,1,-1
               !// Loop downwards
               if (cwete(i,j,l).gt. 0.d0) then
                  KCTOP(i,j) = L
                  exit
               else
                  !// No flux; no detrainment
                  cdetu(i,j,l) = 0.d0
               endif
            enddo
          enddo
        enddo

      !// Find bottom of downdrafts
        KDBOT(:,:) = LWDPAR-1
        do j = 1,jm
          do i = 1,im
            do l = 1,lwdpar-1
               !// Loop downwards
               if (cwetd(i,j,l+1).lt. 0.d0) then
                  !// Mass comes in at top
                  KDBOT(i,j) = L
                  exit
               else
                  !// No flux in at top; no detrainment
                  cdetd(i,j,l) = 0.d0
               endif
            enddo
          enddo
        enddo

c Detrainment cannot be larger than flux in at bottom
        do j = 1,jm
          do i = 1,im
            do l = 1,KCTOP(I,J)
              !// Loop downwards; cdetu and cwete are non-negative
              !// due to filtering above
              if (cdetu(i,j,l) .gt. cwete(i,j,l)) then
                 !// CDETU cannot be larger than the flux in!
                 if (cwete(i,j,l).eq.0.d0) then
                    CDETU(I,J,L) = 0.d0
                 else
                    CDETU(i,j,l) = CWETE(i,j,l)
                 endif
              endif
            enddo
          enddo
        enddo

c Detrainment cannot be larger than flux in at top
        do j = 1,jm
          do i = 1,im
            !// Make sure top is zero
            CDETD(i,j,LWDPAR) = 0.d0
            do L = LWDPAR-1,KDBOT(I,J),-1
              !// Loop downwards; cdetd is non-negative and cwetd is
              !// non-positive due to filtering above
              if (CDETD(i,j,l) .gt. -CWETD(i,j,l+1)) then
                 if (CWETD(i,j,l+1) .lt. 0.d0) then
                    CDETD(i,j,l) = -CWETD(i,j,l+1)
                 else
                    CDETD(i,j,l) = 0.d0
                 endif
              endif
            enddo
          enddo
        enddo

      endif   ! N160 fixes


C FINALLY
c Build entrainment - updrafts (detrainment is positive)
      do j = 1,jm
        do i = 1,im
          !do l = 1,LWEPAR
          do l = 1,KCTOP(I,J)
            !// F(L) + E = F(L+1) + D
            sden = cwete(i,j,l+1) - cwete(i,j,l) + cdetu(i,j,l)
            !// Only allow positive entrainment
            centu(i,j,l) = max(0.d0,sden)
          enddo
        enddo
      enddo
c Build entrainment - downdrafts (flux is negative)
      do l = 1,lwdpar-1
        do j = 1,jm
          do i = 1,im
            !// -F(L+1) + E = -F(L) + D
            sden = cwetd(i,j,l+1) - cwetd(i,j,l) + cdetd(i,j,l)
            !// Only allow positive entrainment
            centd(i,j,l) = max(0.d0,sden)
          enddo
        enddo
      enddo
!//--------------------------------------------------------------------------
      return
      end
!//--------------------------------------------------------------------------

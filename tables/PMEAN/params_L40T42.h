c-----------------------------------------------------------------------
c     parameter file for the UCI CTM core p-5.6 (9/2009)
c-----------------------------------------------------------------------
c  ECMWF Integrated Forecast System, T21/T42/T63/1x1, 40(37)/60 layers
c   T21 Resolution (IPARW =  64, JPARW = 32) 
c           can be nested in T42 or native Gauss grid
c   T42 Resolution (IPARW = 128, JPARW = 64)
c   T63 Resolution (IPARW = 192, JPARW = 96)
c   1x1 Resolution (IPARW = 360, JPARW =180)
c  T319 Resolution (IPARW = 640, JPARW =320)
c           can run native 1x1 or upward multiples:  2x2, 2x3, 3x3, ...

c---T42_L40 run at full resolution and L37 or L40, Gaussian grid N32
      character*10, parameter ::  MODEL = 'EC_T42_L40'
c---1x1_L40 run at full resolution and L37 or L40, non-Gaussian grid n90
c     character*10, parameter ::  MODEL = 'EC_1x1_L40'
c---T42_L60 run at full resolution and L57 or L60, Gaussian grid N32
c     character*10, parameter ::  MODEL = 'EC_T42_L60'
c--T319_L60 run at full resolution and L57 or L60, Gaussian grid N160
c     character*10, parameter ::  MODEL = 'ECT319_L60'

c---resolution of input met fields
      integer, parameter :: IPARW=128, JPARW=64,  LPARW=40   !T42L40 EC grid
c     integer, parameter :: IPARW=360, JPARW=180, LPARW=40   !1X1L40 EC grid
c     integer, parameter :: IPARW=128, JPARW=64,  LPARW=60   !T42L60 EC grid
c     integer, parameter :: IPARW=640, JPARW=320, LPARW=60   !319L60 EC grid

C-----------------------------------------------------------------------
C  Spectral to Grid transformations (Sundet 11/98)
c          NSPGT  = spectral grid type:  3=Quad:T*, 2=Lin:TL*
c          NTRUNW = spectral truncation  (supplied fields)
c          NRNM   = number of spectral coeffs.
c          NMMAX  = number of spectral coeff pairs for wind
      integer, parameter :: NSPGT = 3-(IPARW/300)     ! T42
      integer, parameter :: NTRUNW = (IPARW-1)/NSPGT  ! T42
c     integer, parameter :: NTRUNW = 319              ! 1x1  and T319
      integer, parameter :: NRNM = (NTRUNW+1)*(NTRUNW+2)
      integer, parameter :: NMMAX = (NTRUNW+1)*(NTRUNW+4)/2

C  FFT transforms
      integer, parameter :: NRFFT = 32            ! T42
      integer, parameter :: FFTPTS = IPARW+10     ! T42
c     integer, parameter :: NRFFT = 4             ! 1X1
c     integer, parameter :: FFTPTS = 3*(NTRUNW+2) ! 1X1
c     integer, parameter :: NRFFT = 32            ! T319
c     integer, parameter :: FFTPTS = 3*(NTRUNW+2) ! T319

c---useful constants
      real*8, parameter ::  A0 = 6371000.d0
      real*8, parameter ::  G0 = 9.80665d0
      real*8, parameter ::  CPI    = 3.141592653589793d0
      real*8, parameter ::  C2PI    = 2.d0*CPI
      real*8, parameter ::  CPI180 = CPI/180.D0
      real*8, parameter ::  ZPI180 = 1.d0/CPI180

c---resolution of CTM run (consistent with met fields and degrade)
      integer, parameter :: IPAR=128, JPAR=64, LPAR=37   !**  ECT42L40 grid
c     integer, parameter :: IPAR=360, JPAR=180, LPAR=37  !**  EC1x1L40 grid
c     integer, parameter :: IPAR=128, JPAR=64, LPAR=57   !**  ECT42L60 grid
c     integer, parameter :: IPAR=180, JPAR=90,  LPAR=37  !**  EC2x2 grid
c     integer, parameter :: IPAR=640, JPAR=320, LPAR=57  !**  ECT319L60 grid
      integer, parameter :: LWEPAR=34, LWDPAR=30, LDPAR=1
      integer, parameter :: LCONVM=34
      integer, parameter :: IMDIV=2

c---degrade params to reduce horizontal resolution of read-in met fields
      integer, parameter :: IDGRD=IPARW/IPAR, JDGRD=JPARW/JPAR
     &                     ,NDGRD=IDGRD*JDGRD

c---OpenMP blocks (careful to ensure exact mulltiples)
      integer, parameter ::  MPIPAR=8, MPJPAR=4   !** ECT42 grid
c     integer, parameter ::  MPIPAR=18, MPJPAR=9  !** EC1x1 grid
c     integer, parameter ::  MPIPAR=10, MPJPAR=5  !** EC2x2 grid
c     integer, parameter ::  MPIPAR=32, MPJPAR=16 !** ECT319 grid
      integer, parameter ::  MPBLK=MPIPAR*MPJPAR
      integer, parameter ::  IDBLK=IPAR/MPIPAR, JDBLK=JPAR/MPJPAR

c---max number of tracers (CHEM sets NTM = no. of used tracers):
      integer, parameter ::  NPAR=28  ! if ASAD is on, NPAR=JITR+JITE+JIIN

c---Dimension for 2-D emission tables
c     integer, parameter ::  NE2DS=1     ! 1: mass only, 6: mass plus moments
      integer, parameter ::  NE2DS=6     ! 1: mass only, 6: mass plus moments

c---diagnostics:
      integer, parameter ::  NTDPAR =   9  ! processes for 3D Tendency Budgets
      integer, parameter ::  NTBPAR =  50  ! boxes for tendencies
      integer, parameter ::  NSBPAR =  50  ! stations for time series(all L)
      integer, parameter ::  NSTPAR =  96  ! steps/day of series(=NRMETD*NROPSM)
c---KSTMAX limits disk space from 3-D unf writes (1 T42 tracer = 1.2 Mbytes)
      integer, parameter ::  KSTMAX = 999

c---Parameters for ASAD
c
c  Master Chemistry Switch
      logical, parameter ::  LPASAD=.true.    ! ASAD switch
c
      integer, parameter ::
     &  JITR = 24, ! Number of TR tracers: implicit chem. and transport
     &  JITE = 2,  ! Number of TR tracers: explicit chem. and transport
     &  JISS = 0,  ! Number of SS species: steady state, no transport
     &  JITS = 5,  ! Number of TS species: implicit chem., no transport
     &  JICT = 4,  ! Number of CT/CF species:  constant, no transport
     &  JIIN = 2,  ! Number of IN tracers: no chemistry, transport
     &  JPSPEC = JITR+JITE+JISS+JITS+JICT, ! Number of species: species in ASAD
     &  JPCTR  = JITR+JITS,                ! Number of species: families in ASAD
     &  JPBG   = JITR+JITE+JITS+JICT,      ! Number of species: diag in ASAD
     &  JITE1  = max(JITE,1),              ! in case jite=0
     &  JITS1  = max(JITS,1)               ! in case jits=0

      integer, parameter ::
     &  JPBK = 54,      ! Number of bimolecular reactions. May be 0.
     &  JPTK = 10,      ! Number of trimolecular reactions.
     &  JPPJ = 20,      ! Number of photolysis reactions.
     &  JPHK = 0,       ! Number of heterogeneous reactions.
     &  JPNR = JPBK+JPTK+JPPJ+JPHK    ! Total number of reactions.

      integer, parameter ::
     &   JPNL = LPAR-3  ! Maximum number of levels for chemistry    (up to lpar)


c---Cloud Cover parameters
c
      integer, parameter ::  NQD_  = 4
      integer, parameter ::  NRAN_ = 10007  ! dimension for random number
      integer, parameter ::  CBIN_ = 20     ! # of quantized cld fration bins
      integer, parameter ::  ICA_  = 20000  ! # of Indep Colm Atmospheres

c-----------------------------------------------------------------------
c---Fast-JX parameters v5.3 (prather 6/05)
c
      integer, parameter :: WX_=18 ! Maximum number of wavelength bins
      integer, parameter :: X_=64  ! Maximum # of species requiring J-values
      integer, parameter :: A_=40  ! no. of Aerosol/cloud Mie sets (input data)
      integer, parameter :: W_=8   ! no. of Wavelength bins: 
                                   ! =18 std, =12 trop only, =8 trop-quick
c
c   L2_ = 2*L1_ = 2*LPAR+2 = no. levels in the basic Fast-JX grid (mid-level)
c   N_  = no. of levels in Mie scattering arrays
c       = 4*LPAR + 5 + 2*sum(JADDLV)
c   M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
c   M2_ = 2*M_ = 8, replaces MFIT
c
      integer, parameter :: L1_=LPAR+1, L2_=2*LPAR+2
      integer, parameter :: N_=501, M_=4, M2_=2*M_
c
c--- 4 Gauss pts = 8-stream
      real*8, dimension(M_), parameter  :: EMU = [0.06943184420297d0,
     &   0.33000947820757d0, 0.66999052179243d0, 0.93056815579703d0]
      real*8, dimension(M_), parameter  :: WT  = [0.17392742256873d0,
     &   0.32607257743127d0, 0.32607257743127d0, 0.17392742256873d0]
c
c   SZAMAX  Solar zenith angle cut-off, above which to skip calculation
c   ZZHT    scale height (cm) used above top of CTM ZHL(LPAR+1)
c   MASFAC  Conversion factor for pressure to column density
c   MASFZC  Mass factor  delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      real*8, parameter :: SZAMAX=98.0d0
      real*8, parameter :: ZZHT=5.d5
      real*8, parameter :: RAD=6375.d5     ! Radius of Earth (cm)
      real*8, parameter :: MASFAC=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
c

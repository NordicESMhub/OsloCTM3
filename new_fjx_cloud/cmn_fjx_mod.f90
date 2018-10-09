!---------------------------------------------------------------------------
!    'cmn_FJX_mod.f90'  for fast-JX code v 7.0+ (prather 9/12)
!---------------------------------------------------------------------------
! Adjusted to CTM3.
! Amund Sovde, February 2015
!---------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!---------------------------------------------------------------------------
!MODULE: CMN_FJX_MOD   & CMN_CLD_MOD
!
!DESCRIPTION: CMN_FJX contains fast-JX variables
!
!---------------------------------------------------------------------------
!//-------------------------------------------------------------------------
MODULE CMN_FJX_MOD
  !//-----------------------------------------------------------------------
  USE CMN_SIZE_MOD, ONLY : IPAR, JPAR, LPAR
  !//-----------------------------------------------------------------------
  IMPLICIT NONE
  PUBLIC
  !//-----------------------------------------------------------------------
  ! JXL_: vertical(levels) dim for J-values computed within fast-JX
  INTEGER, PARAMETER ::  JXL_=100, JXL1_=JXL_+1
  ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid
  ! (mid-level)
  INTEGER, PARAMETER ::  JXL2_=2*JXL_+2
  ! WX_  = dim = no. of wavelengths in input file
  INTEGER, PARAMETER ::  WX_=18
  ! X_   = dim = max no. of X-section data sets (input data)
  INTEGER, PARAMETER ::  X_=72
  ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
  INTEGER, PARAMETER ::  A_=41
  ! C_   = dim = no. of cld-data sets (input data)
  INTEGER, PARAMETER ::  C_=16
  ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
  INTEGER, PARAMETER ::  W_=18     ! W_=8, 12 or 18
  ! N_  = no. of levels in Mie scattering arrays
  !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
  INTEGER, PARAMETER ::  N_=601
  ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
  INTEGER, PARAMETER ::  M_=4
  ! M2_ = 2*M_ = 8, replaces MFIT
  INTEGER, PARAMETER ::  M2_=2*M_
  !//-----------------------------------------------------------------------

  ! 4 Gauss pts = 8-stream
  REAL*8, DIMENSION(M_), PARAMETER  ::  &
       EMU = [.06943184420297d0, .33000947820757d0, &
              .66999052179243d0, .93056815579703d0]
  REAL*8, DIMENSION(M_), PARAMETER  ::  &
       WT  = [.17392742256873d0, .32607257743127d0, &
              .32607257743127d0, .17392742256873d0]
  !//-----------------------------------------------------------------------


  ! MASFAC: Conversion factor for pressure to column density
  REAL*8, PARAMETER   ::  &
       MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
  ! ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
  REAL*8, PARAMETER   :: ZZHT = 5.d5
  ! RAD: Radius of Earth (cm)
  REAL*8, PARAMETER   :: RAD = 6375.d5
  ! ATAU: heating rate (factor increase from one layer to the next)
  REAL*8, PARAMETER   :: ATAU = 1.120d0
  ! ATAU0: minimum heating rate
  REAL*8, PARAMETER   :: ATAU0 = 0.010d0
  ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
  INTEGER, PARAMETER  :: JTAUMX = (N_ - 4*JXL_)/2

  !---- Variables in file 'FJX_spec.dat' (RD_XXX)

  ! WBIN: Boundaries of wavelength bins
  REAL*8 :: WBIN(WX_+1)
  ! WL: Centres of wavelength bins - 'effective wavelength'
  REAL*8 :: WL(WX_)
  ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
  REAL*8 :: FL(WX_)

  REAL*8 :: QO2(WX_,3)   ! QO2: O2 cross-sections
  REAL*8 :: QO3(WX_,3)   ! QO3: O3 cross-sections
  REAL*8 :: Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

  ! QQQ: Supplied cross sections in each wavelength bin (cm2)
  REAL*8 :: QQQ(WX_,3,X_)
  ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
  REAL*8 :: QRAYL(WX_+1)
  ! TQQ: Temperature for supplied cross sections
  REAL*8 :: TQQ(3,X_)
  ! LQQ = 1, 2, or 3 to determine interpolation with T or P
  INTEGER :: LQQ(X_)

  ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*6 :: TITLEJX(X_)
  ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*1 :: SQQ(X_)

  !---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)

  ! QAA: Aerosol scattering phase functions
  REAL*8 :: QAA(5,A_)
  ! WAA: 5 Wavelengths for the supplied phase functions
  REAL*8 :: WAA(5,A_)
  ! PAA: Phase function: first 8 terms of expansion
  REAL*8 :: PAA(8,5,A_)
  ! RAA: Effective radius associated with aerosol type
  REAL*8 :: RAA(A_)
  ! SAA: Single scattering albedo
  REAL*8 :: SAA(5,A_)
  ! DAA: density (g/cm^3)
  REAL*8 :: DAA(A_)
  ! NAA: Number of categories for scattering phase functions
  INTEGER :: NAA

  !---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)

  ! QCC: Cloud scattering phase functions
  REAL*8 :: QCC(5,C_)
  ! WCC: 5 Wavelengths for supplied phase functions
  REAL*8 :: WCC(5,C_)
  ! PCC: Phase function: first 8 terms of expansion
  REAL*8 :: PCC(8,5,C_)
  ! RCC: Effective radius associated with cloud type
  REAL*8 :: RCC(C_)
  ! SCC: Single scattering albedo
  REAL*8 :: SCC(5,C_)
  ! DCC: density (g/cm^3)
  REAL*8 :: DCC(C_)
  ! NCC: Number of categories for cloud scattering phase functions
  INTEGER :: NCC

  !---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)

  ! WMM: U Michigan aerosol wavelengths
  REAL*8 :: WMM(6)
  ! UMAER: U Michigan aerosol data sets
  integer, parameter :: n_umset = 44
  REAL*8 :: UMAER(3,6,21,n_umset)

  !---- Variables in file 'atmos_std.dat' (RD_PROF)

  integer, parameter ::  LREF=51   ! layer dim. in reference profiles
  integer, parameter ::  JREF=18   ! latitude dim. in reference profiles

  ! T and O3 reference profiles
  REAL*8, DIMENSION(LREF,JREF,12) :: TREF, OREF

  INTEGER :: NJX, NW1, NW2

  !-----------------------------------------------------------------------
  INTEGER, PARAMETER ::  &
       AN_=25  & ! no of separate aerosols per layer (needs NDX for each)
       ,L_=LPAR     &     ! L_ = number of CTM layers
       ,L1_=L_+1    &
       ,L2_=2*L_+2        ! no. levels in the Fast-JX grid that
                          ! includes both layer edges and layer mid-points

  !-----------------------------------------------------------------------
  ! variables used to map fast-JX J's onto CTM J's
  !-----------------------------------------------------------------------

  INTEGER, PARAMETER :: JVN_=19   ! total number of J-values in ratj.d

  REAL*8  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
  INTEGER JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
  INTEGER NRATJ         ! number of Photolysis reactions in CTM chemistry,
                            ! derived here NRATJ must be .le. JVN_
  INTEGER JO1D          ! index of O3(1D) reaction (not directly used)
  CHARACTER*6 JVMAP(JVN_) ! label of J-value used to match with fast-JX J's
  CHARACTER*50 JLABEL(JVN_) ! label of J-value used in the main chem model

  !-----------------------------------------------------------------------
  ! mass, O3, & T in atmosphere above top of CTM
  !   assumes that PTOP is same for all (I,J) = ETAA(LPAR+1)
  !   pressure levels for O3/T climatology (at 2 km in z*)
  !-----------------------------------------------------------------------
  real*8, dimension(JPAR,12) :: TOPM, TOP3, TOPT

  !-----------------------------------------------------------------------
  ! aerosol
  !-----------------------------------------------------------------------
  real*8, dimension(IPAR,JPAR,L1_) :: AER1P, AER2P  ! aerosol path (g/m2)
  integer,dimension(IPAR,JPAR,L1_) :: AER1N, AER2N  ! aerosol index type

  !-----------------------------------------------------------------------
  ! Cloud Cover parameters
  !-----------------------------------------------------------------------
  integer, parameter :: CBIN_ = 20     ! # of quantized cloud fration bins
  integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
  integer, parameter :: L3RG = 3       ! max-ran overlap groups
  INTEGER, PARAMETER :: NQD_ = 4       ! # of cloud fraction bins (4)
  
  !//-----------------------------------------------------------------------
END MODULE CMN_FJX_MOD
!//-------------------------------------------------------------------------

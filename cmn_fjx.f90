!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// fast-JX parameters and variables.
!//=========================================================================
module cmn_fjx
  !//-----------------------------------------------------------------------
  !// MODULE: cmn_fjx
  !// DESCRIPTION: fast-JX parameters and variables.
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: JPPJ, W_chem, IPAR, JPAR, LPAR
  use cmn_parameters, only: AVOGNR, M_AIR, G0
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------

  !// FASTJ specific input files
  character(len=80) :: INFILE_FJX_SPEC, INFILE_FJX_SCAT, INFILE_FJX_AERO
  character(len=80) :: INFILE_FJX_JS, INFILE_FJX_O1D, INFILE_FJX_CLIM

  !// Fast-JX parameters v5.3 (prather 6/05)

  !// Not yet compatible with 7.1


  !// WX_  = dim = no. of wavelengths in input file
  integer, parameter ::  WX_=18

  !// X_   = dim = max no. of X-section data sets (input data)
  !integer, parameter ::  X_=72 ! jx7.1
  integer, parameter ::  X_=64

  !// A_   = dim = no. of Aerosol/cloud Mie sets (input data)
  integer, parameter ::  A_=41


  !// W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
  integer, parameter ::  W_= W_chem     ! W_=8, 12 or 18

  !// N_  = no. of levels in Mie scattering arrays
  !//     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
  !integer, parameter ::  N_=601 ! jx7.1
  integer, parameter ::  N_=501

  !// M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
  integer, parameter ::  M_=4

  !// M2_ = 2*M_ = 8, replaces MFIT
  integer, parameter ::  M2_=2*M_

  !// C_   = dim = no. of cld-data sets (input data)
  !integer, parameter ::  C_=16



  !// JXL_: vertical(levels) dim for J-values computed within fast-JX
  !integer, parameter ::  JXL_=100, JXL1_=JXL_+1 !jx7.1

  !// JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid
  !// (mid-level)
  !integer, parameter ::  JXL2_=2*JXL_+2 !jx7.1





  !// O3, T climatology
  !//-----------------------------------------------------------------------
  integer, parameter ::  LREF=51   ! layer dim. in reference profiles
  integer, parameter ::  JREF=18   ! latitude dim. in reference profiles

  ! T and O3 reference profiles
  real(r8), dimension(LREF,JREF,12) :: TREF, OREF


  !//-----------------------------------------------------------------------


  ! 4 Gauss pts = 8-stream
  real(r8), dimension(M_), parameter  ::  &
       EMU = [.06943184420297_r8, .33000947820757_r8, &
              .66999052179243_r8, .93056815579703_r8]
  real(r8), dimension(M_), parameter  ::  &
       WT  = [.17392742256873_r8, .32607257743127_r8, &
              .32607257743127_r8, .17392742256873_r8]
  !//-----------------------------------------------------------------------

  !// MASFAC: Conversion factor for pressure to column density
  real(r8), parameter   ::  &
       MASFAC = 100._r8 * AVOGNR/(M_AIR * G0 * 10._r8)

  !// ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
  real(r8), parameter   :: ZZHT = 5.e5_r8

  !// RAD: Radius of Earth (cm)
  real(r8), parameter   :: RAD = 6375.e5_r8

  !// SZAMAX  Solar zenith angle cut-off, above which to skip calculation
  real(r8), parameter   :: SZAMAX=98.0_r8 ! jx6.7


  !// ATAU: heating rate (factor increase from one layer to the next)
  !real(r8), parameter   :: ATAU = 1.120_r8 !jx7.1
  real(r8) :: ATAU

  !// ATAU0: minimum heating rate
  !real(r8), parameter   :: ATAU0 = 0.010_r8 !jx7.1
  real(r8) :: ATAU0

  !// JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
  !integer, parameter  :: JTAUMX = (N_ - 4*JXL_)/2 !jx7.1
  integer :: JTAUMX

  !//-----------------------------------------------------------------------


  !// Spectral data ---- Variables in file 'FJX_spec.dat' (RD_XXX)
  !//-----------------------------------------------------------------------

  !// WBIN: Boundaries of wavelength bins
  real(r8) :: WBIN(WX_+1)
  !// WL: Centres of wavelength bins - 'effective wavelength'
  real(r8) :: WL(WX_)

  !// FL: Solar flux incident on top of atmosphere (cm-2.s-1)
  real(r8) :: FL(WX_)

  real(r8) :: QO2(WX_,3)   ! QO2: O2 cross-sections
  real(r8) :: QO3(WX_,3)   ! QO3: O3 cross-sections
  real(r8) :: Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

  !// QQQ: Supplied cross sections in each wavelength bin (cm2)
  !real(r8) :: QQQ(WX_,3,X_) ! jx7.1
  real(r8) :: QQQ(WX_,2,X_)

  !// QRAYL: Rayleigh parameters (effective cross-section) (cm2)
  !real(r8) :: QRAYL(WX_+1) ! jx7.1
  real(r8) :: QRAYL(WX_)

  !// TQQ: Temperature for supplied cross sections
  real(r8) :: TQQ(3,X_)

  !// LQQ = 1, 2, or 3 to determine interpolation with T or P
  !integer :: LQQ(X_) ! jx7.1

  !// TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
  !character(len=6) :: TITLEJX(X_) !jx7.1
  character(len=7) :: TITLEJ(X_)

  !// SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
  !character(len=1) :: SQQ(X_) !jx7.1


  !// Scattering aerosols ---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)
  !//-----------------------------------------------------------------------

  !// QAA: Aerosol scattering phase functions
  real(r8) :: QAA(5,A_)
  !// WAA: 5 Wavelengths for the supplied phase functions
  real(r8) :: WAA(5,A_)
  !// PAA: Phase function: first 8 terms of expansion
  real(r8) :: PAA(8,5,A_)
  !// RAA: Effective radius associated with aerosol type
  real(r8) :: RAA(A_)
  !// SAA: Single scattering albedo
  real(r8) :: SAA(5,A_)
  !// DAA: density (g/cm^3)
  real(r8) :: DAA(A_)
  !// NAA: Number of categories for scattering phase functions
  integer :: NAA


  !// Scattering clouds ---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)
  !//-----------------------------------------------------------------------

  !// QCC: Cloud scattering phase functions
  !real(r8) :: QCC(5,C_)
  !// WCC: 5 Wavelengths for supplied phase functions
  !real(r8) :: WCC(5,C_)
  !// PCC: Phase function: first 8 terms of expansion
  !real(r8) :: PCC(8,5,C_)
  !// RCC: Effective radius associated with cloud type
  !real(r8) :: RCC(C_)
  !// SCC: Single scattering albedo
  !real(r8) :: SCC(5,C_)
  !// DCC: density (g/cm^3)
  !real(r8) :: DCC(C_)
  !// NCC: Number of categories for cloud scattering phase functions
  !integer :: NCC


  !// Scattering U Michigan ---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)
  !//-----------------------------------------------------------------------
  !// n_umset: U Michigan aerosol
  integer, parameter :: n_umset = 44
  !// WMM: U Michigan aerosol wavelengths
  real(r8) :: WMM(6)
  !// UMAER: U Michigan aerosol data sets
  !real(r8) :: UMAER(3,6,21,n_umset) !jx7.1
  real(r8) :: UMAER(3,5,21,n_umset)


  !// CTM layers
  !//-----------------------------------------------------------------------
  integer, parameter :: &
        L_ = LPAR, &      ! L_ = number of CTM layers
        L1_= L_+1, &
        L2_= 2*L_+2       ! no. levels in the Fast-JX grid that
                          ! includes both layer edges and layer mid-points


  !// variables used to map fast-JX J's onto CTM J's
  !-----------------------------------------------------------------------
  integer, parameter :: JVN_=JPPJ   ! total number of J-values in ratj.d
  real(r8)  :: JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
  real(r8)  :: JFACTAO !jx6.7
  integer :: JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
  integer :: NRATJ         ! number of Photolysis reactions in CTM chemistry,
                           ! derived here NRATJ must be .le. JVN_
  !integer :: JO1D          ! index of O3(1D) reaction (not directly used)
  integer :: JINDO         ! index of O3(1D) reaction (not directly used)
  !// label of J-val used to match with fast-JX J's
  !character(len=6) :: JVMAP(JVN_) !jx7.1
  !// label of J-val used in the main chem model
  !character(len=50) :: JLABEL(JVN_) !jx7.1
  character(len=7) :: JLABEL(JVN_), JLABELO

  logical :: LTROP !jx6.7

  integer :: NJVAL, NW1, NW2 !jx6.7


  !// mass, O3, & T in atmosphere above top of CTM
  !//   assumes that PTOP is same for all (I,J) = ETAA(LPAR+1)
  !//   pressure levels for O3/T climatology (at 2 km in z*)
  !-----------------------------------------------------------------------
  !real(r8), dimension(JPAR,12) :: TOPM, TOP3, TOPT ! jx7.1 ????
  real(r8), dimension(IPAR,JPAR) :: TOPT, TOPM, TOP3
  real(r8), dimension(IPAR,JPAR) :: LPART, LPARM, LPAR3

  !// Mass and O3 taken from STT and possibly climatology
  real(r8), dimension(LPAR,IPAR,JPAR) :: DMS, DO3 !jx6.7


  !// aerosol
  !-----------------------------------------------------------------------
  real(r8), dimension(IPAR,JPAR,L1_) :: AER1P, AER2P  ! aerosol path (g/m2)
  integer,dimension(IPAR,JPAR,L1_) :: AER1N, AER2N  ! aerosol index type

  !// Cloud Cover parameters
  !-----------------------------------------------------------------------
  integer, parameter :: CBIN_ = 20     ! # of quantized cloud fration bins
  integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
  integer, parameter :: L3RG = 3       ! max-ran overlap groups
  integer, parameter :: NQD_ = 4       ! # of cloud fraction bins (4)




  !//-----------------------------------------------------------------------
end module cmn_fjx
!//=========================================================================

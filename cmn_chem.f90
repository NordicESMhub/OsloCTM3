!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
module CMN_CHEM
  !//-----------------------------------------------------------------------
  !// MODULE: CMN_CHEM_MOD
  !// DESCRIPTION: CMN_CHEM_MOD contains emission/chemistry variables for CTM
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only : IPAR, JPAR, LPAR, NPAR, IPARW, JPARW, LPARW, &
       LWEPAR, LWDPAR, LDPAR, NE2DS, TNAMELEN
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  save
  !//-----------------------------------------------------------------------
  !// CTM3: Commented out linoz flux and sink/mass
  !//     LITSRC is reverse order
  !//     TCNION is replaced by SCViceFR and SCV_T258

  !// Tracer specific input files
  character(len=80) :: INFILE_T, INFILE_WET, INFILE_DRY
  character(len=80) :: INFILE_EMIS, INFILE_POLAR_O3LOSS
  character(len=80) :: INFILE_RES, INFILE_MEGAN, INFILE_LIGHTNING

  !// Basic chemistry & init info:
  logical, dimension(NPAR) ::  LZONE        ! for init, scaling of tracer

  !// Tracer specific info
  real(r8), dimension(NPAR) :: &
       TMASS,           & ! molecular weight
       TMASSMIX2MOLMIX, & ! conversion from kg/kg to mole/mole
       TMOLMIX2MASSMIX    ! conversion from mole/mole to kg/kg
  character(len=TNAMELEN) :: TNAME(NPAR)   ! tracer names



  !// Emissions data:
  !//-----------------------------------------------------------------------
  !// 2D stored as pattern with scale factor for each species
  !// Note that different patterns, including monthly stored as new table
  integer, parameter :: ETPAR=1800         !// dim. no. of 2-D tables
  integer :: NE2TBL                        !// actual no. of 2-D tables
  integer :: NM2TBL(ETPAR)   !// Table of months the emissions apply for
  integer :: NY2TBL(ETPAR)   !// Table of years emissions apply for
  real(r8), dimension(NPAR,ETPAR) :: E2STBL  !// scale factor for each N
  real(r8), dimension(IPAR,JPAR,NE2DS,ETPAR):: E2DS !// CTM gridded emission w/XY
  logical, dimension(NPAR,ETPAR) :: E2LTBL !// .T. = species N uses tbl E

  !// Emissions data:  3D 
  integer, parameter :: E3PAR=48
  integer :: NE3TBL
  integer :: NM3TBL(E3PAR) !// Month table
  integer :: NY3TBL(E3PAR) !// Year table
  real(r8), dimension(NPAR,E3PAR) :: E3STBL
  !// Need to change order of E3DS to avoid striding (IPAR,JPAR,LPAR,E3PAR)
  real(r8), dimension(LPAR,IPAR,JPAR,E3PAR):: E3DSNEW  !// no XY moments for 3D
  logical, dimension(NPAR,E3PAR) :: E3LTBL !// .T. = species N uses tbl E

  !// Lightning:
  real(r8) :: &
       !// CTM3: reversed order instead of UCI LITSRC(IPAR,JPAR,LWEPAR)
       LITSRC(LWEPAR,IPAR,JPAR), &  !// Lightning source (% of annual total)
       LITFAC(NPAR)  !// Total mass of lightn. source for each species (kg/yr)
  integer :: &
       NLIT(NPAR), & !// Index of species affected by lightning
       NEMLIT        !// Number of species affected by lightning emissions

  !// LinoZ data
  integer, parameter::  NLINZL = 25
  integer, parameter::  NLINZJ = 18
  integer, parameter::  NLINZT = 7
  real(r8), dimension(JPAR,LPAR,NLINZT)::   LZPML0, LZPML1, LZPML2
  real(r8), dimension(JPAR,LPAR,12,NLINZT):: LZTBL0, LZTBL1, LZTBL2
  real(r8), dimension(JPAR,12)::  LZO3COL         !O3 colm above CTM
  real(r8), dimension(JPAR)::  O3TOPLZ, TACTLZ    ! O3-colm & PSC activ T
  real(r8), dimension(LPAR)::  TPSCLZ, ZPSCLZ     ! PSC: loss freq & SZA
  real(r8) :: O3SFCLZ    !  Surface O3 abundance (mole/mole) for Linoz O3
  real(r8) :: ZN2OSFC, ZCH4SFC
  real(r8) :: O3TAULZ    !  First-order lifetime ==> O3SFCLZ
  real(r8) :: E90VMR_TP    !  Tropopause threshold mix. ratio for E90 tracer
  real(r8) :: O3iso1, O3iso2  !  Tropopause threshold mix. ratio
  real(r8) :: CHLORLZ    !  Chlorine loading (ppb)
  real(r8) :: T195LZ     !  threshold temperature (K) to invoke PSC loss param.
  integer :: LZMIN      !  lowest CTM layer with Linoz data

  integer :: N_LZ       !  tracer number (index) that is using Linoz chem
  integer :: N_STE      !  tracer number (index) for STE
  integer :: Ne90       !  tracer number (index) of e090
  integer :: N_O3       !  tracer number (index) of O3
  integer :: N_NO       !  tracer number (index) of NO
  integer :: N_NO2      !  tracer number (index) of NO2
  integer :: LZLBO3     !  level of boundary layer to reset O3 to O3SFCLZ
  logical, dimension(LPAR,IPAR,JPAR) ::  LSTRATAIR_E90 !inst 3D .T. = strat air
  integer, dimension(IPAR,JPAR) ::  LPAUZTOP !max L of tropospheric air


  !// Tracer transport indices
  !//-----------------------------------------------------------------------
!  INTEGER :: &
!       N_LZ,      & ! tracer number (index) for Linoz chem O3strat
!       Ne90,      & ! tracer number (index) of e090
!       N_O3,      & ! tracer number (index) of O3
!       N_N2O,     & ! tracer number (index) of N2O
!       N_NOY,     & ! tracer number (index) of NOY
!       N_CH4,     & ! tracer number (index) of CH4
!       N_CH4x,    & ! tracer number (index) of CH4x
!       N_NO,      & ! tracer number (index) of NO
!       N_NO2        ! tracer number (index) of NO2

     

  !// Scavenging:  wet dep, dry dep, conv, washout/rainout
  logical :: TCAER(NPAR)
  real(r8), dimension(NPAR)   :: TCWETL     !// convective loss fract
  real(r8), dimension(NPAR)   :: TCHENA     !// Henry constant A
  real(r8), dimension(NPAR)   :: TCHENB     !// Henry constant B
  real(r8), dimension(NPAR)   :: TCKAQA     !// Flag for hard-coded Henry
  real(r8), dimension(NPAR)   :: TCKAQB     !// Henry or kinetically limited
  real(r8), dimension(NPAR)   :: SCV_RETEFF !// ice soluble retention coefficient
  real(r8), dimension(NPAR)   :: SCViceFR   !// gridbox fraction soluble in ice
  integer, dimension(NPAR)  :: SCV_T258   !// ice scavenging below 258K

  real(r8), dimension(LPAR,NPAR) :: TRWETL  !// washout vs. L

  !//-----------------------------------------------------------------------
end module CMN_CHEM
!//=========================================================================

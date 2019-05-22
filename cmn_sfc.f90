!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Surface variables and parameters (not meteorological variables).
!//=========================================================================
module CMN_SFC
  !//-----------------------------------------------------------------------
  !// MODULE: CMN_SFC
  !//
  !// DESCRIPTION: Contains surface variables for CTM.
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, NPAR
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  save
  !//-----------------------------------------------------------------------

  !// Surface and vegetation data:  fraction in each LandSurface type
  integer, parameter ::  NVGPAR = 17    ! number of vegetation types
  !real(r8), dimension(IPAR,JPAR,NVGPAR) :: LS_FRAC
  real(r8), dimension(NVGPAR,IPAR,JPAR) :: landSurfTypeFrac
  !// Stomatal resistance (conductance is inverse of resistance)
  real(r8), dimension(NVGPAR,IPAR,JPAR) :: StomRes
  real(r8), dimension(IPAR,JPAR,12)     :: ZOI, LAI  !rougness, leaf-area
  !// Land use dataset: 2=ISLSCP2, 3=CLM4 PFTs
  integer :: LANDUSE_IDX
  !// The year of the wanted land use dataset
  integer :: LANDUSE_YEAR
  !// File name to read from
  character(len=160) :: fileLandSurfTypeFrac

  !// LAI file name and year
  character(len=160) :: fileLAI
  integer :: LAI_YEAR
  !// Z0 file name and year
  character(len=160) :: fileZOI
  integer :: ZOI_YEAR
  
  !// DRYDEP mOSaic parameter switch
  logical :: LDDEPmOSaic
  !// Parameters and land use type from EMEP
  real(r8), dimension(28,16) :: DDEP_PAR
  
  !// Displacement height [m] for given land type fraction (LS_FRAC)
  real(r8) :: ZPDVT_C3(NVGPAR)

  !// Land-sea masks (ocean=0, land=1, lake=2, small island=3, ice shelf=4
  real(r8), dimension(IPAR,JPAR,5) :: LSMASK

  !// Dry deposition
  real(r8), dimension(NPAR,IPAR,JPAR) :: VDEP      !// dep veloc, mean
  real(r8), dimension(NPAR,3)         :: VDEPIN    !// UCI dep vel for land/ocean/cryo
  !// Stomatal conductance
  real(r8), dimension(IPAR,JPAR)      :: VGSTO3    !// dep veloc, stomata

  !// Growing season
  logical :: LGSMAP
  character(len=160) :: fileGSMAP
  integer, dimension(:,:,:), allocatable  :: GDAY_MAP  !// Daily map of number of growing days
  integer, dimension(IPAR,JPAR)           :: GLEN_MAP  !// Map of length of the growing season
  
  !//-----------------------------------------------------------------------
end module CMN_SFC
!//=========================================================================

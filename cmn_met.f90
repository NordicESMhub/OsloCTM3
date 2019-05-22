!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Meteorological variables and parameters.
!//=========================================================================
module CMN_MET
  !//-----------------------------------------------------------------------
  !// MODULE: CMN_MET
  !//
  !// DESCRIPTION: Contains meteorological variables for CTM.
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IPARW, JPARW, LPARW, &
       LWEPAR, LWDPAR, LDPAR, IDBLK, JDBLK, MPBLK
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  save
  !//-----------------------------------------------------------------------

  !// Met-file names
  character(len=80) :: MET_ROOT   !// Root where metdata resides
  character(len=80) :: MPATH1, MPATH2, MFILE3, PPFDFILE
  character(len=120):: PPFDPATH

  !// Metdata type (e.g. ECMWF_IFS/ECMWF_oIFS/ECMWF_oIFSnc4),
  !// cycle and revision number
  character(len=16) :: metTYPE
  integer :: metCYCLE
  integer :: metREVNR


  !// Meteorological year
  integer :: MYEAR

  !// Resolution as strings
  character(len=3) :: VRES, VRESW !// Vertical resolution
  character(len=4) :: HnativeRES  !// Horizontal resolution metdata


  !// Winds(U & V), water(Q), temperature(T), surf press(P)
  real(r8), dimension(IPAR,JPAR,LPAR) ::  U, V, Q, T
  !//real(r8), dimension(LPAR,IPAR,JPAR) ::  Q_LIJ, T_LIJ
  real(r8), dimension(IPAR,JPAR) ::       P, MSLP

  !// Annual mean surface pressure (hPa)
  real(r8), dimension(IPAR,JPAR) :: PMEAN
  real(r8), dimension(IPARW,JPARW) :: PMEANW

  !// Convective fluxes & cloud properties
  real(r8), dimension(IPAR,JPAR,LWEPAR) :: &
       CWETN, &   !// Convective non-entraining updraft flux [kg/s] (not used)
       CWETE, &   !// Convective updraft flux [kg/s] (also covers non-entr.)
       CENTU, &   !// Convective entrainment into updrafts [kg/s]
       PRECLS, &  !// Precipitation, large scale [kg/s]
       PRECCNV, & !// Precipitation, convective [kg/s]
       CLDFR, &   !// Cloud fraction [0,1]
       CLDLWC, &  
       CLDIWC
  !// Convective fluxes & cloud properties
  !real(r8), dimension(LWEPAR,IPAR,JPAR) :: &
  !     PRECLS_LIJ, &  !// Precipitation, large scale [kg/s]
  !     PRECCNV_LIJ, & !// Precipitation, convective [kg/s]


  !// Downdrafts
  real(r8), dimension(IPAR,JPAR,LWDPAR) :: &
       CWETD, & !// Convective downdraft flux [kg/s]
       CENTD    !// Convective entrainment into downdrafts [kg/s]



  !// 3D fields set up for column treatment (LPAR,IPAR,JPAR)
  real(r8), dimension(LPAR,IPAR,JPAR) :: &
       PVU, &  !// Potential vorticity
       UMS, &  !// Winds [m/s] (used in physics routines, not in transport)
       VMS     !// Winds [m/s] (used in physics routines, not in transport)


  !// Surface meteorology 
  real(r8), dimension(IPAR,JPAR)  :: &
       SA, &   !// Surface albedo [0,1] (forecast albedo)
       SLH, &
       SHF, &
       SMF, &
       SFT, &
       SFQ, &
       SFU, &
       SFV, &
       !SPRECLS, &    !// Surface value of PRECLS_LIJ [kg/s]
       BLH, &        !// Boundary layer height [m], see LBLH for model level.
       BLH_CUR, &    !// Boundary layer height [m] for current time step
       BLH_NEXT, &   !// Boundary layer height [m] for next time step
       USTR, &
       CI, &         !// Sea ice cover [0,1]
       SD, &         !// Snow depth [m water equivalent]
       PhotActRad, & !// Photosynthetically active radiation [W/m2]
       PPFD, &       !// Photosynthetically active radiation - deaccumulated [W/m2]
       SWVL1, &      !// Soil volumetric water level 1 0-7 cm [m3/m3]
       SWVL3, &      !// Soil volumetric water level 3 28-100 cm [m3/m3]
       STL1          !// Soil temperature level 1 [K]

  integer, dimension(IPAR,JPAR) :: &
       LBLH          !// Uppermost model level for PBL

  real(r8), dimension(IDBLK,JDBLK,MPBLK) :: &
       SMLT, &       !// Snow melt [m water equivalent per seconds]
       ES, &         !// Snow evaporation [m water equivalent per seconds]
       SNFL, &       !// Snow fall [m water equivalent per seconds]
       PBL_KEDDY, &  !// KEDDY (KH) for PBL
       MO_LENGTH, &  !// Monin-Obukhov (MO_LENGTH) from PBL
       PRANDTLL1     !// Prandtl for surface layer




  !// Altitudes
  real(r8), dimension(LPAR+1,IPAR,JPAR) :: ZOFLE



  !//-----------------------------------------------------------------------
end module CMN_MET
!//=========================================================================

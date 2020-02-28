!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// Main variables used for running Oslo CTM3.
!//=========================================================================
module CMN_CTM
  !//-------------------------------------------------------------------------
  !// MODULE: CMN_CTM
  !//
  !// DESCRIPTION: CMN_CTM contains major variables for running the CTM.
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, rMom
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IPARW, JPARW, &
       LPARW, IDGRD, JDGRD, MPBLK, NMMAX, NRMETD
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  save
  !//-----------------------------------------------------------------------


  !// Grids: coords for lat & lng, radians & deg, mid & edges, areas
  !//-----------------------------------------------------------------------
  real(r8), dimension(JPAR)    :: &
       YGRD,   & ! latitude mid-point of grid boxes in radians
       YDGRD,  & ! latitude mid-point of grid boxes in degree
       DISTY     ! latitudinal distance of grid boxes in km

  real(r8), dimension(JPAR+1)  :: &
       YEDG,   & ! latitude edges of grid in radians
       YDEDG,  & ! latitude edges of grid in degree
       DISTX     ! longitudinal distance of grid in km

  real(r8), dimension(IPAR)    :: &
       XGRD,   & ! longitude mid-point of grid boxes in radians
       XDGRD     ! longitude mid-point of grid boxes in degree

  real(r8), dimension(IPAR+1)  :: &
       XEDG,   & ! longitude edges of grid in radians
       XDEDG     ! longitude edges of grid in degree

  real(r8), dimension(LPAR)    :: &
       ZGRD    ! center of vertical grid (hPa)

  real(r8), dimension(LPAR+1)  :: &
       ZEDG,       & ! edge of vertical grid (hPa)
       ETAA, ETAB    ! edge of ETA coords: (etaa + etab*P_sfc)

  real(r8), dimension(IPAR,JPAR) :: &
       AREAXY        ! area of grid boxes (m^2)
  real(r8), dimension(IPARW,JPARW) :: &
       AREAXYW       ! Native area of grid boxes (m^2)

  !// Grid: ETA coords of the met fields and LMMAP remapping to CTM grid
  real(r8), dimension(LPARW+1) :: ETAAW, ETABW, XLMMAP
  integer, dimension(LPARW+1) :: LMMAP
  !// Grid: native lat/lon
  real(r8), dimension(JPARW)    :: &
       YDGRDW     ! latitude mid-point of grid boxes in degree
  real(r8), dimension(JPARW+1)    :: &
       YDEDGW     ! latitude edges of grid in degree
  real(r8), dimension(IPARW)    :: &
       XDGRDW     ! longitude mid-point of grid boxes in degree
  real(r8), dimension(IPARW+1)  :: &
       XDEDGW     ! longitude edges of grid in degree

  !// Grid mappings:  includes regular & degraded grid mappings (lower res)
  logical :: LDEG              ! .TRUE. if horizontal degradation used
  integer :: IMAP(IDGRD,IPAR)  ! Index of high-res longitude grid box
  integer :: JMAP(JDGRD,JPAR)  ! Index of high-res latitude grid box
  !// Fraction of low-res box accounted for by high-res box
  real(r8)  :: ZDEGI(IDGRD,IPAR), ZDEGJ(JDGRD,JPAR)

  !// Latitude variables
  real(r8)  ::  WGLYE(JPARW+1), WGLYG(JPARW)

  !// Strings for lat/lon/height printouts
  character(len=4) :: &
       TLAT(JPAR), TLATE(JPAR+1), &
       TLNG(IPAR), TLNGE(IPAR+1), &
       TALT(LPAR), TALTE(LPAR+1)

  !// Grid data: areas, pressure at edges of box (PIJL) for PMEAN
  real(r8), dimension(IPAR,JPAR,LPAR+1) :: XYZA, XYZB, PIJL
  real(r8), dimension(IPAR,JPAR) :: XYA, XYB
  !// Grid data: annual mean ht(m), land fract
  real(r8), dimension(IPAR,JPAR) :: PALTD, PLAND
  !// Grid data: map for averaging met field (p,T,Q) over extended polar zones
  integer, dimension(JPAR)  :: IMEPZ
  !// Total area and mean surface pressure
  real(r8) :: AREAG, PMEANG



  !// Transport mass fluxes (U=Alfa, V=Beta, W=Gama)
  !//-----------------------------------------------------------------------
  real(r8), dimension(IPAR+1,JPAR,LPAR) :: ALFA
  real(r8), dimension(IPAR,JPAR+1,LPAR) :: BETA
  real(r8), dimension(IPAR,JPAR,LPAR+1) :: GAMA

  !// Air Mass and Tracer Mass and Moments
  !//-----------------------------------------------------------------------
  real(r8), dimension(IPAR,JPAR,LPAR) :: &
       AIR,   &
       AIRX  ! dry-air mass expected, based on PCTM and eta-levels

  real(r8), dimension(IPAR,JPAR,LPAR,NPAR) :: STT ! Tracer Mass

  real(rMom), dimension(IPAR,JPAR,LPAR,NPAR) :: &
       SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW  ! Tracer Mass Moments


  character(len=10) :: ANAME               ! air
!  character(len=12) :: METTYPE             ! label metdata, e.g. ECMWF_IFS


  !// Calendar stuff
  !//-----------------------------------------------------------------------
  integer :: IYEAR, IDAY, JDAY, JYEAR, JMON, JDATE
  integer :: JDAY_NEXT, JYEAR_NEXT, JMON_NEXT, JDATE_NEXT !// next day
  character(len=3) :: TMON ! 3-char month (JAN, FEB, ...)
  character(len=3) :: TMET ! 3-char day of year (001, 002, ...)
  real(r8) :: GMTAU ! Greenwich time
  real(r8) :: modelTimeIntegrated ! Seconds since model start [s]

  !// Indices
  !//-----------------------------------------------------------------------

  !// IM/JM/LM should be EXACTLY IPAR/JPAR/LPAR, so they are removed.
  !// NTM is just NPAR, but let us assume you can compile the CTM and
  !// run with fewer tracers (NTM<=NPAR).
  !INTEGER :: IM      ! window longitude dimension
  !INTEGER :: JM      ! window latitude dimension
  !INTEGER :: LM      ! window vertical dimension

  !// Number of transported tracers as variable
  integer :: NTM

  !// Number of met fields per day (NRMETD). This is now a parameter
  !// defined in cmn_size.F90, but is accessible through cmn_ctm.f90.
  !// It should be a parameter because it is fixed to the chosen set
  !// of meteorological data, and because it makes NMET diagnostics
  !// easier.

  !// Number of operator splits per met fields
  integer :: NROPSM
  !// Number of chem calculations per operator split
  integer :: NRCHEM
  !// Seconday Order Moments (SOM) limiter index
  integer :: LMTSOM


  !// Current run data
  !//-----------------------------------------------------------------------
  real(r8) :: SOLDEC  ! solar declination (degrees)
  real(r8) :: SOLDIS  ! distance to sun (A.U.)
  real(r8) :: CFLLIM  ! CFL limiter
  logical :: LLPYR   ! Allow leap year (true: looks for day=901 for Feb 29)
  logical :: LYEAR   ! Actual flag for when a year is leap year
  logical :: LFIXMET ! Recycle meteorological year
  logical :: LCONT   ! Continue from restart file
  logical :: LNCR    ! Switch between old (.sav) and new (.nc) restart files
  integer :: START_AVG  !// For init when LCONT=F (0:zero, 1: avgsav file)

  integer :: IDTLN   ! dateline index
  integer :: NSCX    ! Scavenging scheme
  integer :: NDPX    ! Deposition scheme
  integer :: NBLX    ! Boundary layer mixing scheme
  
  !// Diagnostic output
  !//-----------------------------------------------------------------------
  logical, dimension(13) :: LDUMP3HRS ! Switch on 3 hourly output of defined tracers
  logical :: LSTOM1HRS ! Switch on 1 hourly output of GstO3 and FstO3
  ! Switch on daily total/2d scavening output channels (burden, large-scale,...)
  logical, dimension(7) :: LDLYSCAV  
  
  !// OpenMP Blocks  (MPBLK: total number of OpenMP blocks)
  !//-----------------------------------------------------------------------
  integer, dimension(MPBLK) :: &
       MPBLKIB,   & ! beginning longitude index in main domain
       MPBLKIE,   & ! end longitude index in main domain
       MPBLKJB,   & ! beginning latitude index in main domain
       MPBLKJE      ! end latitude index in main domain

  !// Mapping I,J to II,JJ,MP
  integer, dimension(3,IPAR,JPAR) :: all_mp_indices


  !// (cmn_s.f)---spectral to grid transforms --  Sundet 11/98
  !// spectral transform data
  !//-----------------------------------------------------------------------
  real(r8) :: SS(NMMAX)            ! weights for vorticity to gen U and V
  real(r8) :: DD(NMMAX)            ! weights for divergence to gen U and V
  real(r8) :: ALP(NMMAX,JPARW/2)   ! asocciate legendre polys for grid
  real(r8) :: ALPV(NMMAX,JPARW/2)  ! asocciate legendre polys for boundary
  real(r8) :: TRIG(IPARW)          ! trig. functions for FFT at grid
  real(r8) :: TRIGU(IPARW)         ! trig. functions for FFT at boundary
  integer :: IFAX(10)            ! FFT factors 


  !//-----------------------------------------------------------------------
end module CMN_CTM
!//=========================================================================

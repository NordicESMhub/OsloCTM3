!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// Size parameters and flags for CTM, determined by ifdefs in Makefile.
!//=========================================================================
module CMN_SIZE
  !//-----------------------------------------------------------------------
  !// MODULE: CMN_SIZE
  !// DESCRIPTION: Contains size parameters for CTM, and also
  !//              flags defining modules to include.
  !//
  !// This module contains C-code (ifdefs) to set up resolution based
  !// on compiling tokens.
  !//-----------------------------------------------------------------------
  !// GRID PARAMETERS
  !//   IPARW      = met fields global longitude dimension
  !//   JPARW      = met fields global latitude dimension
  !//   LPARW      = met fields max number of sigma levels 
  !//   IPAR       = window longitude dimension
  !//   JPAR       = window latitude dimension
  !//   LPAR       = window vertical dimension
  !//   LWEPAR     = updraft / wet convection vertical dimension
  !//   LWDPAR     = downdraft / dry convection vertical dimension
  !//
  !//   For nested grids: IPAR<IPARW, JPAR<JPARW.
  !//
  !//   LGAUGRD    = logical parameter for Gaussian grid
  !//
  !// OpenMP block dimension (careful to ensure exact mulltiples)
  !//   MPIPAR, MPJPAR
  !//
  !// Tracer parameters
  !//   NPAR       = Number of transported tracers
  !//   NOTRPAR    = Number of non-transported tracers
  !//
  !// Spectral to Grid transformations (Sundet 11/98)
  !//   NSPGT     = spectral grid type:  3=Quad:T*, 2=Lin:TL*
  !//   NTRUNW    = spectral truncation  (supplied fields)
  !//   NRNM      = number of spectral coeffs.
  !//   NMMAX     = number of spectral coeff pairs for wind
  !//
  !// FFT transforms
  !//   NRFFT, FFTPTS
  !//-----------------------------------------------------------------------
  implicit none
  public
  !//-----------------------------------------------------------------------
  character(len=9), parameter ::  MODEL = 'Oslo CTM3'
  character(len=4), parameter ::  MODEL_VERSION = 'v1.1'

  !// Metdata: Number of meteorological time steps per day
  !//-----------------------------------------------------------------------
  integer, parameter :: NRMETD     = 8

  !// Metdata: Native horizontal resolution
  !//-----------------------------------------------------------------------
#if defined(HNATIVET42)
  integer, parameter :: IPARW      = 128
  integer, parameter :: JPARW      = 64
  logical, parameter :: LGAUGRD    = .TRUE.

  integer, parameter :: NSPGT   = 3-(IPARW/300)
  integer, parameter :: NTRUNW  = (IPARW-1)/NSPGT
  integer, parameter :: NRNM    = (NTRUNW+1)*(NTRUNW+2)
  integer, parameter :: NMMAX   = (NTRUNW+1)*(NTRUNW+4)/2
  integer, parameter :: NRFFT   = 32
  integer, parameter :: FFTPTS  = IPARW+10

#elif defined(HNATIVET159)
  integer, parameter :: IPARW      = 320
  integer, parameter :: JPARW      = 160
  logical, parameter :: LGAUGRD    = .TRUE.

  integer, parameter :: NSPGT   = 3-(IPARW/300)
  integer, parameter :: NTRUNW  = 159
  integer, parameter :: NRNM    = (NTRUNW+1)*(NTRUNW+2)
  integer, parameter :: NMMAX   = (NTRUNW+1)*(NTRUNW+4)/2
  integer, parameter :: NRFFT   = 32
  integer, parameter :: FFTPTS  = 3*(NTRUNW+2)

#elif defined(HNATIVET319)
  integer, parameter :: IPARW      = 640
  integer, parameter :: JPARW      = 320
  logical, parameter :: LGAUGRD    = .TRUE.

  integer, parameter :: NSPGT   = 3-(IPARW/300)
  integer, parameter :: NTRUNW  = 319
  integer, parameter :: NRNM    = (NTRUNW+1)*(NTRUNW+2)
  integer, parameter :: NMMAX   = (NTRUNW+1)*(NTRUNW+4)/2
  integer, parameter :: NRFFT   = 32
  integer, parameter :: FFTPTS  = 3*(NTRUNW+2)

#endif

  !// Metdata: Native vertical resolution
  !//-----------------------------------------------------------------------
#if defined(VNATIVEL40)
  integer, parameter :: LPARW    = 40
  integer, parameter :: LWEPARW  = 38
#elif defined(VNATIVEL60)
  integer, parameter :: LPARW    = 60
  integer, parameter :: LWEPARW  = 40
#endif



  !// Model run resolution (multiples of native horizontal resolution)
  !//-----------------------------------------------------------------------
#if defined(HORIGINAL)
  !// Original, i.e. same as native, resolution
  integer, parameter :: IPAR       = IPARW
  integer, parameter :: JPAR       = JPARW

#elif defined(HTWO)
  !// Put 2 native boxes into each lat/lon direction, combining
  !// 2x2 boxes in total
  integer, parameter :: IPAR       = IPARW/2
  integer, parameter :: JPAR       = JPARW/2

#elif defined(HFOUR)
  !// Put 4 native boxes into each lat/lon direction, combining
  !// 4x4 boxes in total
  integer, parameter :: IPAR       = IPARW/4
  integer, parameter :: JPAR       = JPARW/4

#endif /* HORIGINAL */

  !// Collapsed vertical or not?
  !// IMDIV: number of combined pipes in vertical advection. It can
  !// speed up the code in vectorized machines, but will otherwise
  !// slow the computer due to striding for IMDIV>1. In general
  !// this cost more than the gain in a longer pipe.
  !// IMDIV must be even number for L37/L57 (due to qvect3).
#if defined(COLLAPSE)
  integer, parameter :: LPAR       = LPARW - 3
  integer, parameter :: IMDIV      = 2
#else
  integer, parameter :: LPAR       = LPARW
  integer, parameter :: IMDIV      = 1
#endif /* COLLAPSE */

  !// Convective parameters (MUST be even number)
  !//-----------------------------------------------------------------------
  integer, parameter :: LWDPAR  = 30
  integer, parameter :: LDPAR   = 1
#if defined(COLLAPSE)
  integer, parameter :: LWEPAR  = 34
  integer, parameter :: LCONVM  = 34
#else
  integer, parameter :: LWEPAR  = LWEPARW
  integer, parameter :: LCONVM  = LWEPARW
#endif


  !// MP-block sizes (check HNATIVEx against Hy)
  !//-----------------------------------------------------------------------
#if defined(HNATIVET42)
#if defined(HORIGINAL)
  integer, parameter :: MPIPAR  = 2
  integer, parameter :: MPJPAR  = JPAR
#else
  integer, parameter :: MPIPAR  = 1
  integer, parameter :: MPJPAR  = JPAR
#endif

#elif defined(HNATIVET159)
#if defined(HORIGINAL)
  !// With dynamic scheduling, MPIPAR of 4 and 8 are rather similar
  integer, parameter :: MPIPAR  = 8
  integer, parameter :: MPJPAR  = JPAR
#elif defined(HTWO)
  !// With dynamic scheduling, MPIPAR of 2 and 4 are rather similar
  integer, parameter :: MPIPAR  = 2
  integer, parameter :: MPJPAR  = JPAR
#else
  integer, parameter :: MPIPAR  = 1
  integer, parameter :: MPJPAR  = JPAR
#endif

#elif defined(HNATIVET319)
#if defined(HORIGINAL)
  !// It may be that MPIPAR 16 works well. Not tested.
  integer, parameter :: MPIPAR  = 8
  integer, parameter :: MPJPAR  = JPAR
#elif defined(HTWO)
  integer, parameter :: MPIPAR  = 8
  integer, parameter :: MPJPAR  = JPAR
#elif defined(HFOUR)
  integer, parameter :: MPIPAR  = 2
  integer, parameter :: MPJPAR  = JPAR
#else
  integer, parameter :: MPIPAR  = 1
  integer, parameter :: MPJPAR  = JPAR
#endif

#endif

  !// OpenMP blocks
  integer, parameter :: MPBLK   = MPIPAR*MPJPAR   ! number of OpenMP blocks
  integer, parameter :: IDBLK   = IPAR/MPIPAR  ! block longitude dimension
  integer, parameter :: JDBLK   = JPAR/MPJPAR  ! block latitude dimension

  !// degrade parameters to reduce horizontal resolution of met fields
  integer, parameter :: IDGRD   = IPARW/IPAR
  integer, parameter :: JDGRD   = JPARW/JPAR
  integer, parameter :: NDGRD   = IDGRD*JDGRD



  !// Oslo CTM3 chemistry parameters: Tracers and module parameters
  !//-----------------------------------------------------------------------
  integer, parameter :: TNAMELEN = 10 !// Length of TNAME/XTNAME strings
#if defined(OSLOCHEM)
  integer, parameter :: MAX_NOTRPAR=15 !// Max tracers not transported

  !// Find number of components fromdifferent CTM3 packages
  !// E90-tracer?
#if defined(E90)
  integer, parameter :: NPAR_E90      = 1
  logical, parameter :: Le90          = .true.
#else
  integer, parameter :: NPAR_E90      = 0
  logical, parameter :: Le90          = .false.
#endif /* E90 */

  !// LINOZ tracer?
#if defined(LINOZ)
  integer, parameter :: NPAR_LINOZ    = 1
  logical, parameter :: LLINOZ        = .true.
#else
  integer, parameter :: NPAR_LINOZ    = 0
  logical, parameter :: LLINOZ        = .false.
#endif /* LINOZ */

  !// Lightning generation tracer
#ifdef LITGEN
      integer, parameter :: NPAR_LIT=1
#else
      integer, parameter :: NPAR_LIT=0
#endif /* LITGEN */

  !// Tropospheric chemistry parameters
#if defined(TROPCHEM)
  integer, parameter :: NPAR_TROP     = 39
  integer, parameter :: NOTRPAR_TROP  = 7
  logical, parameter :: LOSLOCTROP    = .true.
#else
  integer, parameter :: NPAR_TROP     = 0
  integer, parameter :: NOTRPAR_TROP  = 0
  logical, parameter :: LOSLOCTROP    = .false.
#endif /* TROPCHEM */

  !// Stratospheric chemistry parameters
#if defined(STRATCHEM)
  integer, parameter :: NPAR_STRAT    = 38
  integer, parameter :: NOTRPAR_STRAT = 7
  logical, parameter :: LOSLOCSTRAT   = .true.
#else
  integer, parameter :: NPAR_STRAT    = 0
  integer, parameter :: NOTRPAR_STRAT = 0
  logical, parameter :: LOSLOCSTRAT   = .false.
#endif /* STRATCHEM */

  !// Sulphur parameters
#if defined(SULPHUR)
  integer, parameter :: NPAR_SUL      = 5
  logical, parameter :: LSULPHUR      = .true.
#else
  integer, parameter :: NPAR_SUL      = 0
  logical, parameter :: LSULPHUR      = .false.
#endif /* SULPHUR */

  !// Black & Organic carbon parameters
#if defined(BCOC)
  integer, parameter :: NPAR_BC       = 6
  integer, parameter :: NPAR_OM       = 8
  logical, parameter :: LBCOC         = .true.
#else
  integer, parameter :: NPAR_BC       = 0
  integer, parameter :: NPAR_OM       = 0
  logical, parameter :: LBCOC         = .false.
#endif /* BCOC */

  !// Secondary organic carbon parameters
#if defined(SOA)
  integer, parameter :: NPAR_SOA      = 56
  logical, parameter :: LSOA          = .true.
#else
  integer, parameter :: NPAR_SOA      = 0
  logical, parameter :: LSOA          = .false.
#endif /* SOA */

  !// Sea salt parameters
#if defined(SEASALT)
  integer, parameter :: NPAR_SALT     = 8
  logical, parameter :: LSALT         = .true.
#else
  integer, parameter :: NPAR_SALT     = 0
  logical, parameter :: LSALT         = .false.
#endif /* SEASALT */

  !// Mineral DUST parameters
#if defined(DUST)
  integer, parameter :: NPAR_DUST     = 8
  logical, parameter :: LDUST         = .true.
#else
  integer, parameter :: NPAR_DUST     = 0
  logical, parameter :: LDUST         = .false.
#endif /* DUST */

  !// Nitrate parameters
#if defined(NITRATE)
  integer, parameter :: NPAR_NITRATE  = 5
  logical, parameter :: LNITRATE      = .true.
#else
  integer, parameter :: NPAR_NITRATE  = 0
  logical, parameter :: LNITRATE      = .false.
#endif /* NITRATE */

  !// General Oslo CTM3 parameters
  !// max number of transported tracers
  integer, parameter ::  NPAR    = NPAR_TROP + NPAR_STRAT + NPAR_SUL &
                                   + NPAR_BC + NPAR_OM &
                                   + NPAR_SALT + NPAR_DUST &
                                   + NPAR_NITRATE + NPAR_E90 + NPAR_LINOZ &
                                   + NPAR_SOA + NPAR_LIT
  !// max number of non-transported tracers
  integer, parameter ::  NOTRPAR = NOTRPAR_TROP + NOTRPAR_STRAT
  !// max number of all tracers
  integer, parameter ::  TOT_TRACER    = NPAR + NOTRPAR
  integer, parameter ::  TRACER_ID_MAX = 300    ! max value of tracer number
  logical, parameter ::  LOSLOCHEM     = .true. ! Oslo chemistry switch



#if defined(STRATCHEM)
  !// Number of species in Oslo chemistry (strat+trop) requiring J-values
  integer, parameter :: JPPJ          = 51 !// Number of photolysis reactionsW
  integer, parameter :: W_chem        = 18
#else
  !// Number of species in Oslo chemistry (trop only) requiring J-values
  integer, parameter :: JPPJ          = 51 !// Number of photolysis reactions
  integer, parameter :: W_chem        = 8
#endif /* STRATCHEM */


  !// Treatment of emissions and deposition; Oslo chemistry
#if defined(EMISDEPCHEM)
  !// Emissions and deposition inside chemistry
  logical, parameter :: LEMISDEP_INCHEM =.true.
  !// Dimension for 2-D emission tables
  !// 1: mass only, 6: mass plus moments
  integer, parameter :: NE2DS           = 1
#else
  !// Emissions and deposition as separate processes
  logical, parameter :: LEMISDEP_INCHEM = .false.
  !// Dimension for 2-D emission tables
  !// 1: mass only, 6: mass plus moments
  integer, parameter :: NE2DS           = 6
#endif /* EMISDEPCHEM */

#else   /* OSLOCHEM */
  logical, parameter :: LOSLOCHEM       = .false. ! Oslo chemistry switch
  integer, parameter :: JPPJ            = 20 !// Number of photolysis reactions
  !// Treatment of emissions and deposition; UCI style
  !// Emissions and deposition as separate processes
  logical, parameter :: LEMISDEP_INCHEM = .false.
  !// Dimension for 2-D emission tables
  !// 1: mass only, 6: mass plus moments
  integer, parameter :: NE2DS           = 6
  integer, parameter :: NPAR            = 1
#endif /* OSLOCHEM */


  !// Diagnostics
  integer, parameter ::  NTDPAR =   8  ! processes for 3D Tendency Budgets
  integer, parameter ::  NTBPAR =  50  ! boxes for tendencies
  integer, parameter ::  NSBPAR =  50  ! stations for time series(all L)
  integer, parameter ::  NSTPAR =  96  ! steps/day of series(=NRMETD*NROPSM)


  !//-----------------------------------------------------------------------
end module CMN_SIZE
!//=========================================================================

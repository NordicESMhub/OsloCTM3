!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// CTM global variables (should probably be put elsewhere)
!//=========================================================================
module cmn_oslo
  !//-----------------------------------------------------------------------
  !// MODULE: cmn_oslo
  !//
  !// DESCRIPTION: Module defining global variables for Oslo chemistry.
  !//
  !// Some variables are of size (LPAR,IPAR,JPAR), others are
  !// (LPAR,IDBLK,JDBLK,MPBLK). The latter is to reduce striding for an
  !// IJ-block. For each IJ-block (there are MPBLK of them), there will
  !// be as good as no striding.
  !//
  !// Covers
  !//   - physical variables
  !//   - emission variables
  !//   - chemistry variables
  !//   - diagnose variables
  !//   - help variables
  !//   - tracer related variables
  !//
  !// Amund Sovde, October 2008 - July 2010
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, rAvg
  use cmn_size, only: IPAR, JPAR, LPAR, LWEPAR, IDBLK, JDBLK, MPBLK, &
       NPAR, NOTRPAR, NPAR_DUST, TRACER_ID_MAX, TNAMELEN
  use cmn_sfc, only: NLCAT
  use cmn_chem, only: ETPAR
  use cmn_fjx, only: JVN_
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------
  !// All variables are to be saved. This should be default for global
  !// variables, although Fortran claims it may not be.
  save
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  !// To reduce striding for an IJ-block, most of the larger arrays defined
  !// here are not on the form (IPAR,JPAR) but (IDBLK,JDBLK,MPBLK).
  !// This is to reduce striding in each IJ-block.
  !//-----------------------------------------------------------------------


  !// --- Physical variables -----------------------------------------------
  !// Uppermost level of troposphere
  integer, dimension(IPAR,JPAR) :: LMTROP
  !// Select method to define LMTROP
  !//   TP_TYPE=1: PVU-based tropopause.
  !//   TP_TYPE=2: Based on -dT/dz.
  !//   TP_TYPE=3: Based on E90 tracer (LPAUZTOP).
  integer, parameter :: TP_TYPE = 1
  !// When no stratospheric chemistry is used, calculate tropospheric
  !// chemitstry LVS2ADD2TC levels above LMTROP.
  integer, parameter :: LVS2ADD2TC = 3

  !// Background particle surface area density
  real(r8), dimension(LPAR,IPAR,JPAR) :: PARTAREA


  !// --- For convective washout -------------------------------------------
  !// Flag for type of dissolved calculation (set in wetset)
  integer, dimension(NPAR) :: TCCNVHENRY
  !// Flag for convective plume temperatures (258K and 273K)
  integer, dimension(2,IDBLK,JDBLK,MPBLK) :: LELEVTEMP

  !// Fraction: convective precipitation / liquid water in convective plume
  real(r8), dimension(LWEPAR,IDBLK,JDBLK,MPBLK) :: QFRAC
  !// Excess convective rain to be treated as large scale rain
  !real(r8), dimension(LWEPAR,IDBLK,JDBLK,MPBLK) :: LSEXTRA

  !// Fraction: volume concentration of elevator liquid water mass
  !//           i.e. elevator liquid water volume / elevator volume
  real(r8), dimension(LWEPAR,IDBLK,JDBLK,MPBLK) :: LW_VOLCONC



  !// --- For air ----------------------------------------------------------
  !// Air molecules defined every meteorological time step
  real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK) :: AIRMOLEC_IJ

  !// Grid box volume defined every meteorological time step
  real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK) :: DV_IJ


  !// --- Tracer related variables -----------------------------------------
  !// Array for mapping chemical ids to transport numbers in STT, and reverse
  integer, dimension(TRACER_ID_MAX) :: trsp_idx
  integer, dimension(NPAR)          :: chem_idx
  !// Array for mapping chemical ids to no-transport numbers in XSTT,
  !// and reverse
  integer, dimension(TRACER_ID_MAX) :: Xtrsp_idx
  integer, dimension(NOTRPAR)       :: Xchem_idx

  !// Name array for non-transported tracers
  character(len=TNAMELEN), dimension(NOTRPAR) :: XTNAME
  !// Tracer data:   mol.wt. & conversion from kg/kg to mole/mole
  real(r8), dimension(NOTRPAR) :: XTMASS
  real(r8), dimension(NOTRPAR) :: XTMASSMIX2MOLMIX
  real(r8), dimension(NOTRPAR) :: XTMOLMIX2MASSMIX

  !// --- Non-transported species (should be out-phased) -------------------
  !// Non-transported species defined in extra STT
  real(r8), dimension(LPAR,NOTRPAR,IPAR,JPAR) :: XSTT

  !// Average species defined in extra STT:
  !// To be removed when the chemistry core diagnostics work.
  real(rAvg), dimension(LPAR,NOTRPAR,IPAR,JPAR)  :: XSTTAVG

  !// --- Chemical stuff ---------------------------------------------------
  !// J-values for chemistry (NPHPAR,LPAR,IDBLK,JDBLK,MPBLK).
  real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK,MPBLK) :: JVAL_IJ
  !// Control J-values in internal chemistry/boundary layer mixing loop
  logical :: LJCCYC

  !// Oslo 2D arrays upper / lower boundaries
  real(r8) :: &
       STT_2D_LB(JPAR,10,47), &!// Lower 10 layers
       STT_2D_LT(JPAR,64), &   !// Uppermost layer
       O3_UPBND(JPAR)          !// Flux of O3 from stratosphere

  !// Count the number of corrected negative ozone after tropchem.
  !// Generally only happens in surface level.
  real(r8) :: TROPCHEMnegO3(LPAR,MPBLK)

  !// Aerosol surface conversion of N2O5 to HNO3
  real(r8) :: PR42HET(LPAR,JPAR)

  !// Fields for CH4 surface treatment (for fixed field)
  real(r8) :: CH4FIELD(IPAR,JPAR)


  !// --- Emissions --------------------------------------------------------
  !// Emissions for chemistry treatment [kg/s]
  real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK,MPBLK) :: EMIS_IJ
  !// Accumulated emissions [kg] to be used in STTTND and STTTN0
  real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK,MPBLK) :: DIAGEMIS_IJ
  real(r8), dimension(NPAR,366) :: emisTotalsDaily
  real(r8), dimension(NPAR)     :: emisTotalsOld

  !// Decide whether we want to use methane emissions or methane surface
  !// fields based on monthly means
  logical, parameter ::  METHANEMIS = .false.

  !// --- Additions for 2D emissions ---------------------------------------
  !// Number of categories, given in E2D_CATNAMES
  !// 1: AGR: Agricultural
  !// 2: ENE: Energy
  !// 3: IND: Industrial
  !// 4: TRA: Traffic
  !// 5: RCO: Residential combustion, other
  !// 6: SLV: Solvent use
  !// 7: WAS: Waste burning
  !// 8: SHI: Ship emissions
  !// 9: BIO: Biogenic emissions
  !//10: OCN: Oceanic emissions
  !//11: NAT: Natural sources
  !//12: EXF: Extraction distribution of fossil fuels
  !//13: POW: Power generation (if not ENE is used)
  !//14: AWB: Agricultural waste burning emissions
  !//15: RES: Residential (if not RCO is used)
  integer, parameter :: NECAT = 16
  character(len=3),dimension(NECAT), parameter :: ECATNAMES= &
       (/'AGR','ENE','IND','TRA','RCO','WAS','SLV','WAS','SHI','BIO', &
         'OCE','NAT','EXF','POW','AWB','RES'/)

  !// Table mapping 2D emission table to category number
  integer, dimension(ETPAR) :: E2CTBL


  !// Scalings for 2D emissions - horizontal and vertical
  !// ----------------------------------------------------------------------
  !// Local hour scaling is done in two ways:
  !// A. Using predefined scalings for each local hour.
  !// B. By calculating factors for each grid box and scaling with
  !//    average values e.g. for the month.
  !// For simplicity, these are treated separately, but it is not
  !// possible to apply both to a dataset.
  !//
  !// Vertical scaling is used to distribute surface emissions on the
  !// lowermost 4 levels.

  ! This is only implemented for Oslo chemistry in emisdep4chem.f90,
  ! not in source_uci.f90.

  !// A. Table mapping 2D emission table to diurnal variation index
  integer, dimension(ETPAR) :: E2LocHourTBL

  !// Pre-defined entries of diurnal scaling based on local hour
  !// 0: None; use dataset as is, together with global scalings.
  !// 1: RETRO variations (TNO)
  integer, parameter :: E2LocHourRETRO = 1
  !// 2: +50% from 8am to 7pm, -50% from 8pm to 7am (e.g. for BCOC)
  integer, parameter :: E2LocHour5050  = 2

  !// Max number of LocHour types (i.e. not dependent on 2D field input)
  integer, parameter :: NE2LocHourVARS = 2
  !// Table of LocHour (diurnal) variations
  real(r8), dimension(24,NECAT,NE2LocHourVARS) :: E2LocHourSCALE


  !// B. Table mapping 2D emission table to horizontal 2D variation index
  integer, dimension(ETPAR) ::   E22dTBL

  !// Pre-defined entries of diurnal scaling from 2D fields
  !// 1: Scale with both sunlight and temperature
  integer, parameter :: E22dSunTemp = 1
  !// 2: Scale with sunlight
  integer, parameter :: E22dSun     = 2
  !// 3: Scale with heating degree day (HDD)
  integer, parameter :: E22dHDD     = 3

  !// Max number of 2D variation types
  integer, parameter :: NE22dVARS = 3
  !// Table of 2D variations
  real(r8), dimension(IDBLK,JDBLK,MPBLK,NE22dVARS) :: E22dSCALE
  !// Keep track of how many 2D scalings are used
  logical, dimension(NE22dVARS) :: E2L2dTBL


  !// Table mapping 2D emission to vertical distribution on
  !// several layers (defined by NE2vertLVS and hardcoded in code)
  integer, dimension(ETPAR) ::   E2vertTBL
  !// Max number of 2D variation types
  integer, parameter :: NE2vertVARS = 3
  !// Number of layers to distribute on
  integer, parameter :: NE2vertLVS  = 4
  !// Table of vertical variations
  real(r8), dimension(NE2vertLVS,NE2vertVARS) :: E2vertSCALE




  !// --- Specific 3D emissions --------------------------------------------
  !// Forest fires - only some layers
  integer :: FF_TYPE, FF_YEAR, NEFIR ! Type, year, # of components emitted
  character(len=120) :: FF_PATH
  integer, parameter :: EPAR_FIR = 23    !// Number of components
  integer, parameter :: EPAR_FIR_LM = 33 !// Number of layers
  integer,dimension(EPAR_FIR) :: ECOMP_FIR
  !// Emission array used in model
  real(r8), dimension(EPAR_FIR_LM,EPAR_FIR,IDBLK,JDBLK,MPBLK) :: EMIS_FIR
  !// CTM components to be emitted
  character(len=10), dimension(EPAR_FIR) :: FF_CNAMES
  !// Biomass burning prefix for file name
  character(len=20), dimension(EPAR_FIR) :: FF_BNAMES
  !// Biomass burning dataset name
  character(len=20), dimension(EPAR_FIR) :: FF_VARNAME
  !// Global scaling of dataset
  real(r8), dimension(EPAR_FIR) :: FF_SCALE
  !// Partitions: (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)
  integer, parameter :: EPAR_NPARTS=6
  real(r8), dimension(EPAR_NPARTS, EPAR_FIR) :: FF_PARTITIONS

  !// GFED4 emissions factors
  integer, parameter :: GFED4_NC = 41 !// Number of GFED4 species
  integer, parameter :: GFED4_NP = 6  !// For now: must equal EPAR_NPARTS
  real(r8), dimension(GFED4_NP,GFED4_NC) :: GFED4_EF !// emission factors
  character(len=16), dimension(GFED4_NC) :: GFED4_NM !// name of species


  !// Lightning - see LITSRC in cmn_chem.f90 and lightning.f90

  !// Aircraft - see emissions_aircraft.f90

  !// Ship - treated as monthly 2D surface emissions


  !// --- Global diagnose variables ----------------------------------------
  !// Diagnose convective washout
  real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK,MPBLK) :: CONVWASHOUT

  !// Diagnose column of species [DU] each NMET
  real(r8), dimension(IPAR,JPAR,24) :: dobson_snapshot, dobson_snapshot_ts

  !// Average temperature
  real(rAvg), dimension(IPAR,JPAR,LPAR) :: TEMPAVG
  !// Average H2O
  real(rAvg), dimension(IPAR,JPAR,LPAR) :: H2OAVG
  !// Average H2O mixing ratio just above LMTROP
  real(rAvg), dimension(IPAR,JPAR) :: H2OAVG_LMT1
  !// Average Q
  real(rAvg), dimension(IPAR,JPAR,LPAR) :: QAVG
  !// Average air density (AIR_MOLEC)
  real(rAvg), dimension(IPAR,JPAR,LPAR) :: AMAVG
  !// Average LMTROP
  real(rAvg), dimension(IPAR,JPAR) :: LMTROPAVG


  !// --- Diagnose global scavenging (LS/CNV/DRY) --------------------------
  !// Total scavenged within each IJ-block [kg]
  real(r8), dimension(NPAR,MPBLK) :: &
       SCAV_LS, SCAV_CN, SCAV_DD
  !// Accumulated burden -> average burden [kg]
  real(r8), dimension(NPAR,MPBLK) :: SCAV_BRD
  !// Total scavenged for each process
  real(r8), dimension(NPAR,4,366) :: SCAV_DIAG !// LS, CNV, DRY and BRD
  !// Maps (surface) of each scavenging process
  real(r8), dimension(NPAR,IDBLK,JDBLK,MPBLK) :: &
       SCAV_MAP_WLS, SCAV_MAP_WCN, SCAV_MAP_DRY

  !//--- Diagnose stomatal ozone uptake ------------------------------------
  !// Stomatal conductance: daily average
  real(r8), dimension(IPAR,JPAR,NLCAT) :: GSTO3_AVG
  !// Stomatal flux: daily average
  real(r8), dimension(IPAR,JPAR,NLCAT) :: FSTO3_AVG
  !// Dry deposition velocities
  real(r8), dimension(IPAR,JPAR,NLCAT) :: VRAO3_AVG
  real(r8), dimension(IPAR,JPAR,NLCAT) :: VRBO3_AVG
  real(r8), dimension(IPAR,JPAR,NLCAT) :: VRCO3_AVG
  
  !// --- Help variables ---------------------------------------------------
  !// Days of month
  integer, dimension(12) :: DINM

  !// Flag initialized components
  integer :: ZEROINIT(NPAR), XZEROINIT(NOTRPAR)

  !// Path for result files
  character(len=160) :: RESULTDIR

  !// Radii of dust bins
  real(r8) :: dustbinsradii(NPAR_DUST)


  !//-----------------------------------------------------------------------
end module cmn_oslo
!//=========================================================================

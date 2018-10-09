!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// Diagnostic variables for Oslo CTM3.
!//=========================================================================
module CMN_DIAG
  !// ----------------------------------------------------------------------
  !// MODULE: CMN_DIAG_MOD
  !// DESCRIPTION: CMN_DIAG_MOD contains diagnostic variables for CTM
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, r4, rTnd, rAvg
  use cmn_size, only : IPAR, JPAR, LPAR, NPAR, MPBLK, NTDPAR, NTBPAR, &
       NSBPAR, NSTPAR
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  save
  !//-----------------------------------------------------------------------

  !// Title read from input file and printed to some diag files.
  character(len=80) :: RUNTITLE
  !// Metdata info for diagnostic files
  character(len=80) :: metTypeInfo
  character(len=120) :: resolutionInfo

  !// Calendar for CONTinuation = write a restart file
  integer, dimension(366)              :: JDO_C

  !// Tendencies
  !//-----------------------------------------------------------------------
  integer, dimension(366)      :: JDO_T !// Calendar for tendencies
  integer                      :: NTND  !// Actual # of tendency processes
  logical :: LBGA1, LBGT1(NPAR) !// Controlling 1D output (Air+Tracers)
  logical :: LBGA2, LBGT2(NPAR) !// Controlling 2D output (Air+Tracers)
  character(len=6), dimension(NTDPAR) :: TLDIAG !// Names of processes
  !// Process indices in tendency array
  integer :: NTND_SOURCE, NTND_BNDLYR, NTND_DRYDEP, NTND_CHEM, &
             NTND_LSSCAV, NTND_CNSCAV, NTND_WADV, NTND_UVADV

  !// Tendency arrays
  real(r8) :: STTBCK(IPAR,JPAR,LPAR,NPAR)          !// STT before process
  real(r8) :: STTTN0(NPAR,NTDPAR)                  !// Global tendency
  real(rTnd) :: STTTND(IPAR,JPAR,LPAR,NPAR,NTDPAR) !// 3D tendency

  !// Boxes for tendencies (& 1-D avg profiles): 
  integer :: NBOXD
  character(len=15), dimension(NTBPAR) :: TBOXD
  integer, dimension(2,NTBPAR) :: IBOXD, JBOXD, KBOXD
  integer, dimension(2,NTBPAR,MPBLK) :: IBOXDMP, JBOXDMP
  logical :: LBOXDMP(NTBPAR,MPBLK)


  !// netCDF4 - global packing parameters
  !//-----------------------------------------------------------------------
  integer, parameter :: nc4deflate_global = 9
  integer, parameter :: nc4shuffle_global = 1


  !// mixing ratio output
  !//-----------------------------------------------------------------------
  integer, parameter :: NDTRMX = 1
  logical            :: LTRMX(NPAR), LMXDG
  integer, dimension(366) :: JDO_X
  integer            :: NTRMX
  real(r4)             :: STTMX(IPAR,JPAR,LPAR,NDTRMX)
  character(len=10), dimension(NDTRMX) :: TTRMX

  !// Flux diag?
  !//-----------------------------------------------------------------------
  logical :: LFLXDG

  !// Averages (STT and more)
  !//-----------------------------------------------------------------------
  integer, dimension(366)      :: JDO_A
  logical :: LBGTA(NPAR)
  integer :: NRAVG
  real(rAvg) :: &
       STTAVG(IPAR,JPAR,LPAR,NPAR),&
       AIRAVG(IPAR,JPAR,LPAR), &
       DVAVG(IPAR,JPAR,LPAR), &
       ZHAVG(LPAR+1,IPAR,JPAR), &
       PSFCAVG(IPAR,JPAR)

  !// Local time tracer diag.
  real(r8) :: GM0_LT(IPAR+1), LTGBL1(5), LTGBL2(5)
  real(r8), dimension(NSBPAR)       :: LTLAT,LTLNG ,LTSTN1,LTSTN2
  integer, dimension(NSBPAR)      :: ILTX,JLTX, NCNTLTS, LTSTID
  integer :: NBOXLT, NRGBLT, NRGLTD(5), LTGLID(5), NCNTDAY
  character(len=10), dimension(NSBPAR) :: TLTAX, TLTRST
  character(len=10) :: TLTRGL(5)
  logical :: LTSTNSV(100,NSBPAR), LTGBLSV(IPAR,100,5)
  real(r4) :: STTLTS(LPAR,NSBPAR), STTLTG(IPAR,JPAR,LPAR,5)

  !// Series of tracer profile at stations 
  integer :: NBOXS
  integer, dimension(366)              :: JDO_S
  logical :: LBGTS(NPAR)
  character(len=10), dimension(NSBPAR) :: TSTAX
  real(r8), dimension(NSBPAR)            :: STLAT,STLNG
  integer, dimension(NSBPAR)           :: ISTA, JSTA
  real(r4), dimension(LPAR,NSTPAR,NSBPAR,NPAR) :: STTS4

  !// date records for tendencies (0) and averages (1):
  integer :: NDAY0,JDAY0,JYEAR0,JMON0,JDATE0
  integer :: NDAY1,JDAY1,JYEAR1,JMON1,JDATE1
  real(r8) :: TAU0, TAU1
  character(len=3) :: TMON0, TMON1

  !// multiple U and V steps in QVECT, and accumulated NADV
  integer :: USTEP(LPAR),VSTEP(LPAR),WSTEP(MPBLK),MNADV

  !// Linoz tracer flux diagnostics
  integer, parameter :: NFLX = 4
  real(r8)  :: UFLX(IPAR+1,JPAR,NFLX), VFLX(IPAR,JPAR+1,NFLX)
  integer :: IDAY0_STE, ITAU0_STE, JDATE0_STE, JMON0_STE, JYEAR0_STE
  !// Linoz boundary sink and trop ozone mass
  real(r8) :: &
       TROPMASS(IPAR,JPAR,NFLX), TROPMASS_0(IPAR,JPAR,NFLX), &
       TROPM(IPAR,JPAR,NFLX), O3PML(IPAR,JPAR,NFLX), &
       O3WSCAV(IPAR,JPAR,NFLX), &
       O3PBLSINK(IPAR,JPAR,NFLX),N2OPBLSINK(IPAR,JPAR)

  !//-----------------------------------------------------------------------
end module CMN_DIAG
!//=========================================================================

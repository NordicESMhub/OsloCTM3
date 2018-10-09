c---(cmn_h.f)---main (head) common blocks p-5.5 (9/2007)-------

      implicit none
      include 'params.h'

c---Title for cmn blks
      character*40, parameter::   TITCMN = 
     &    ' UCI CTM: common blocks     p-5.5 9/2007 '

c---Met Field Model title set in params.h: char*10, param::MODEL='EC_T42_L40'
c---Titles of run and labels for tracer names, grid, dates, met fields ids...
      character*80 RTITLE
      character*10 TNAME(NPAR),ANAME
      character*4  TLAT(JPAR),TLATE(JPAR+1),TLNG(IPAR),TLNGE(IPAR+1)
      character*4  TALT(LPAR),TALTE(LPAR+1)
      character*3  TMON,TMET

c---Air Mass and Tracer Mass and Moments
      real*8, dimension(IPAR,JPAR,LPAR)       ::  AIR, AIRX
      real*8, dimension(IPAR,JPAR,LPAR,NPAR)  ::  STT
      real*4, dimension(IPAR,JPAR,LPAR,NPAR)  ::  SUT,SVT,SWT, SUU,SVV,
     &               SWW, SUV,SUW,SVW
      
c---Indices:   NTM = number of transported tracers
      integer        IM,JM,LM,NTM, NRMETD,NROPSM,NRCHEM, LMTSOM

c---OpenMP Blocks
      integer, dimension(MPBLK) :: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE

c---Calendar data (integer)
      integer        IYEAR,IDAY,JDAY,JYEAR,JMON,JDATE,MYEAR
c---Current run data: UT (0. to 24. hr), solar decl, solar dist.,CFL lim
      real*8         GMTAU, SOLDEC,SOLDIS, CFLLIM

c---Tracer data:   mol.wt. & conversion from kg/kg to mole/mole
      real*8, dimension(NPAR) ::  TMASS, TMMVV

c---Logical switches for leap year, Gaussian grid, recycling met fields
      logical LLPYR, LGAUGRD, LFIXMET, LCONT

c---Grid data: areas, pressure at edges of box (PIJL) for PMEAN
      real*8, dimension(IPAR,JPAR,LPAR+1) :: XYZA, XYZB, PIJL
      real*8, dimension(IPAR,JPAR) :: XYA, XYB, AREAXY
c---Grid data: annual mean surf press, ht(m), land fract, ...
      real*8, dimension(IPAR,JPAR) :: PMEAN, PALTD, PLAND
c---Grid data: coords for lat & lng, radians & deg, mid & edges 
      real*8, dimension(JPAR)      :: YGRD, YDGRD, DISTY
      real*8, dimension(JPAR+1)    :: YEDG, YDEDG, DISTX
      real*8, dimension(IPAR)      :: XGRD, XDGRD
      real*8, dimension(IPAR+1)    :: XEDG, XDEDG
c---Grid:  Z__ = center/edge of vertical grid (km), ETA coord's
      real*8, dimension(LPAR)      :: ZGRD
      real*8, dimension(LPAR+1)    :: ZEDG, ETAA, ETAB
c---Grid: ETA coords of the met fields and LMMAP remapping to CTM grid
      real*8, dimension(LPARW+1)   :: ETAAW, ETABW, XLMMAP
      integer,dimension(LPARW+1)   :: LMMAP
c---Grid data:  map for averaging met field (p,T,Q) over extended polar zones
      integer, dimension(JPAR)  :: IMEPZ
c---Indices:  dateline; params to scaveng, deposition, boundary layer      
      integer     IDTLN,   NSCX,NDPX,NBLX
c---Diagnostics: total area and mean surface pressure
      real*8                     AREAG, PMEANG

c-----------------------------------------------------------------------
      common/CCH_R8/ AIR,AIRX,STT, 
     &               TMASS,TMMVV,GMTAU,CFLLIM,SOLDEC,SOLDIS
      common/CCH_R4/ SUT,SVT,SWT,SUU,SVV,SWW,SUV,SUW,SVW
      common/CCH_I/ IM,JM,LM,NTM, NRMETD,NROPSM,NRCHEM, LMTSOM,
     &              IYEAR,IDAY,JDAY,JYEAR,JMON,JDATE,MYEAR,LMMAP,
     &              IMEPZ,IDTLN, NSCX,NDPX,NBLX,
     &              MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
      common/CCH_L/ LLPYR,LGAUGRD,LFIXMET,LCONT
      common/CCH_T/ RTITLE,TNAME,ANAME,
     &              TLAT,TLATE,TLNG,TLNGE,TALT,TALTE,TMET,TMON
      
      common/CCG_R8/ XYZA,XYZB,PIJL,XYA,XYB,AREAXY,PMEAN,PALTD,PLAND,
     &     YGRD,YDGRD,YEDG,YDEDG,DISTY, XGRD,XDGRD,XEDG,XDEDG,DISTX,
     &     ZGRD,ZEDG,ETAA,ETAB,ETAAW,ETABW,XLMMAP, AREAG,PMEANG
c-----------------------------------------------------------------------


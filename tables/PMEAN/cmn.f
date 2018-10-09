c---(cmn_h.f)---main (head) common blocks p-5.5 (9/2007)-------

      implicit none
      include 'params.h'

c---Met Field Model title set in params.h: char*10, param::MODEL='EC_T42_L40'
c---Titles of run and labels for tracer names, grid, dates, met fields ids...
      character*80 RTITLE
      character*10 TNAME(NPAR),ANAME
      character*4  TLAT(JPAR),TLATE(JPAR+1),TLNG(IPAR),TLNGE(IPAR+1)
      character*4  TALT(LPAR),TALTE(LPAR+1), TMET
      character*3  TMON

c---Air Mass and Tracer Mass and Moments
      real*8, dimension(IPAR,JPAR,LPAR)       ::  AIR, AIRX
      real*8, dimension(IPAR,JPAR,LPAR,NPAR)  ::  STT
      real*4, dimension(IPAR,JPAR,LPAR,NPAR)  ::  SUT,SVT,SWT, SUU,SVV,
     &               SWW, SUV,SUW,SVW
      
c---Indices:   NTM = number of transported tracers
      integer        IM,JM,LM,NTM, NRMETD

c---Calendar data (integer)
      integer        IYEAR,IDAY,JDAY,JYEAR,JMON,JDATE,MYEAR
c---Current run data: UT (0. to 24. hr), solar decl, solar dist.,CFL lim
      real*8         GMTAU

c---Tracer data:   mol.wt. & conversion from kg/kg to mole/mole
      real*8, dimension(NPAR) ::  TMASS, TMMVV

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
c---Diagnostics: total area and mean surface pressure
      real*8   AREAG, PMEANG
      logical  LLPYR, LGAUGRD, LFIXMET

c-----------------------------------------------------------------------
      common/CCH_R8/ AIR,AIRX,STT, 
     &               TMASS,TMMVV,GMTAU
      common/CCH_R4/ SUT,SVT,SWT,SUU,SVV,SWW,SUV,SUW,SVW
      common/CCH_I/ IM,JM,LM,NTM, NRMETD,
     &              IYEAR,IDAY,JDAY,JYEAR,JMON,JDATE,MYEAR,LMMAP,
     &              IMEPZ
      common/CCH_T/ RTITLE,TNAME,ANAME,
     &              TLAT,TLATE,TLNG,TLNGE,TALT,TALTE,TMET,TMON
      common/CCH_L/ LLPYR,LGAUGRD,LFIXMET
      
      common/CCG_R8/ XYZA,XYZB,PIJL,XYA,XYB,AREAXY,PMEAN,PALTD,PLAND,
     &     YGRD,YDGRD,YEDG,YDEDG,DISTY, XGRD,XDGRD,XEDG,XDEDG,DISTX,
     &     ZGRD,ZEDG,ETAA,ETAB,ETAAW,ETABW,XLMMAP, AREAG,PMEANG
c-----------------------------------------------------------------------

c---Met-file names ---------------------------------------
      character*40 MPATH1, MPATH2, MFILE3

c---Met Fields: Winds(U & V), water(Q), temperature(T), surf press(P)
      real*8, dimension(IPAR,JPAR,LPAR) ::  Q, T
      real*8, dimension(IPAR,JPAR) :: P, AREAXYW
      common/CCW_R8/   Q,T,P, AREAXYW
      common/CCW_T/    MPATH1,MPATH2,MFILE3

c---(cmn_s.f)---spectral to grid transforms --  Sundet 11/98

c---spectral transform data
      real*8   SS(NMMAX)            ! weights for vorticity to gen U and V
      real*8   DD(NMMAX)            ! weights for divergence to gen U and V
      real*8   ALP(NMMAX,JPARW/2)   ! asocciate legendre polys for grid
      real*8   ALPV(NMMAX,JPARW/2)  ! asocciate legendre polys for boundary
      real*8   TRIG(IPARW)          ! trig. functions for FFT at grid
      real*8   TRIGU(IPARW)         ! trig. functions for FFT at boundary
      integer  IFAX(10)           ! FFT factors 

c---grid data:  includes regular & degraded grid mappings (lower res)
      real*8   WGLYE(JPARW+1), WGLYG(JPARW)
      logical  LDEG                   ! .TRUE. if horizontal degradation used
      integer  IMAP(IDGRD,IPAR), JMAP(JDGRD,JPAR)
      real*8   ZDEGI(IDGRD,IPAR), ZDEGJ(JDGRD,JPAR)

c-----------------------------------------------------------------------
      common/CCS_R8/   SS, DD, ALP, ALPV, TRIG, TRIGU, 
     &                 WGLYE,WGLYG,ZDEGI,ZDEGJ
      common/CCS_I/    IFAX, IMAP,JMAP
      common/CCS_L/    LDEG

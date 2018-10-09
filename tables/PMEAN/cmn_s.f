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


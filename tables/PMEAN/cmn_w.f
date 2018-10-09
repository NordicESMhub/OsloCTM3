c---(cmn_w.f)---met-field & related common blocks p-5.5 (9/2007)-------

c---Met-file names ---------------------------------------
      character*40 MPATH1, MPATH2, MFILE3

c---Met Fields: Winds(U & V), water(Q), temperature(T), surf press(P)
      real*8, dimension(IPAR,JPAR,LPAR) ::  U, V, Q, T
      real*8, dimension(IPAR,JPAR) ::       P, MSLP
      real*8, dimension(IPARW,JPARW) :: PMEANW, AREAXYW

c---Met Fields: Convective fluxes & cloud properties
      real*8, dimension(IPAR,JPAR,LWEPAR) :: CWETN,CWETE,CENTU,PRECMC,
     &         PRECLS,CLDFR, CNVCF, STRCF, CLDLWC,CLDIWC, ODW,ODI,OD
      real*8, dimension(IPAR,JPAR,LWDPAR) ::  CWETD,CENTD
      real*8, dimension(IPAR,JPAR,LDPAR)  ::  CDRY

c---Surface meteorology 
      real*8, dimension(IPAR,JPAR)  ::  SA,SLH,SHF,SMF,SFT,SFQ,SFU,SFV,
     &                                  BLH,USTR

c---Derived: Mass fluxes (U=Alfa, V=Beta, W=Gama)
      real*8, dimension(IPAR+1,JPAR,LPAR) :: ALFA
      real*8, dimension(IPAR,JPAR+1,LPAR) :: BETA
      real*8, dimension(IPAR,JPAR,LPAR+1) :: GAMA

c---Fractional Cloud data/storage
c       number of cld fraction bins:  NQD_=4  defined in params.h
      logical LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
      real*8  CLDSTORE(LPAR+1,NQD_,IPAR,JPAR),SWSTORE(NQD_,IPAR,JPAR)
      integer TYPSTORE(LPAR+1,IPAR,JPAR), RANSEED
     
c---Surface and vegetation data:  fraction in each LandSurface type
      integer, parameter ::  NVGPAR = 17    ! number of vegetation types
      real*8, dimension(IPAR,JPAR,NVGPAR) :: LS_FRAC
      real*8, dimension(IPAR,JPAR,12)     :: ZOI, LAI  !rougness, leaf-area
      real*8, dimension(IPAR,JPAR,NPAR)   :: VDEP      !dep veloc, mean

c---Random number set
      real*4 RAN4(NRAN_)
      
c---Altitudes and TOPs
      real*8  ZOFLE(LPAR+1,IPAR,JPAR)
      real*8, dimension(IPAR,JPAR) :: TOPT, TOPM, TOP3

c-----------------------------------------------------------------------
      common/CCW_R8/   U,V,Q,T,P,MSLP, CWETN,CWETE,CENTU,PRECMC,PRECLS,
     &      CLDFR,CNVCF,STRCF,CLDLWC,CLDIWC,ODW,ODI,OD,CWETD,CENTD,CDRY,
     &      SA,SLH,SHF,SMF,SFT,SFQ,SFU,SFV,BLH,USTR,  CLDSTORE,SWSTORE,
     &      ALFA,BETA,GAMA,  LS_FRAC,ZOI,LAI, VDEP,
     &      ZOFLE, PMEANW, AREAXYW, TOPT,TOPM,TOP3
      common/CCW_R4/   RAN4
      common/CCW_I/    TYPSTORE, RANSEED
      common/CCW_L/    LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
      common/CCW_T/    MPATH1,MPATH2,MFILE3




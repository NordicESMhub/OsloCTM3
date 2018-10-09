c---(p-phot.f)-------generic CTM shell from UCIrvine (p-code 5.6, 6/08)
c     CTM3: Modified for CTM3; PHOTOL removed, -NAER fixed
c           Current code revised to JX ver 6.7
c
c           Amund Sovde, December 2012
c ---------------------------------------------------------------------
c
c version 6.7  (should output same J-values as ver 6.6 EXCEPT for VOCs)
c       New JPL-2010 update complete through all VOCs (see notes in FJX_spec.dat)
c           Most important change is incorporation of Stern-Volmer pressure-dependent
c           and wavelength-dependent quantum yields into the Temperature interpolation.
c           Acetone now has only one pair of X sections for each pathway.
c       Redo of mapping of (external) chemical model J's onto fastJX J's
c           see examples of splitting a J with fixed branching q-yields
c           see how chem model labels are not tied by the cross section labels in fastJX
c       Changes in FX_spec.dat make it incompatible with earlier versions, although
c           all of the Xsection data have the same format.
c       Now the number of X sections read in equals the number of J's that fast-JX calculates.
c       As before, all J's except for O2, O3, O3(1D) have only pairs of data at different T.
c
c version 6.6  (should output same J-values as ver 6.5)
c       Redo of the 4x4 matrix inversions from subroutine to inline code (faster)
c       Small changes in read-in of cross sections to allow for a reduced set.
c           Previous X-section data required "Acetone" to be last and included
c            extra tables for special pressure dependence formulation.
c           Now, IF acetone is included as the last X-scetion, it reads the tables.
c       N.B. *** It appears that with the new, updated solar fluxes, the J-NO is too large.
c           The J-O2 and O3 steady-state are imporved, but the NOy:N2O is too small.
c           No obvious fix at present, possibly due lack of thermopsheric NO absorption.
c
c version 6.5
c >>>> J-values are now averaged over CTM layer (and not just mid-layer point)
c        This is important when cloud optical depth is large (~1).
c
c version 6.4
c >>>> allows for shortened, speeded up troposphere versions<<<<<<
c     STD:           W_=18
c        identical results to v-6.2 if cloud OD is consistent
c     TROP-ONLY:     W_=12
c        collapses the wavelength bins from 18 to 12 (5-6-7-8 & 11-18) 
c        drops many 'stratospheric' cross-sections (denoted by 'x' in 2nd title)
c        allows use of single standard spectral data set:  FJX_spec.dat
c        results close to W_=18, largest difference is J-O2 (<1% in 13-18 km!!)
c        This is recommended as accurate for troposphere only calculations.
c     TROP-QUICK:    W_=8
c        reverts to original fast-J 7-bins (12-18) plus 1 scaled UV (5) for J-O2
c        errors in 12-18 km range for J-O2, high sun are 10%, worse above.
c     ***Photolysis of O2 in the upper tropical troposphere is an important 
c        source of O3.  It needs to be included in tropospheric runs.  
c        TROP-ONLY is recommended, W_=8 is a quick fix if speed essential.
c     
c     Major rewrite of code to minimize calls and allow better vector-type ops.
c     loop over wavelengths internal to Mie soln.
c     Driven by profiling of CTM code, may still need optimization.
c     Wavelengths can be contracted to W_=12 (trop only) and strat-only 
c        X-sections are dropped.  With parm W_=18, the std fast-JX is retrieved.
c     Many call eliminated and summations in BLKSLV and GEN_ID are explicit
c     GEN_ID replaces GEN and calculates all matrix coeff's (1:L_) at once
c     RD_XXX changed to collapse wavelengths & x-sections to Trop-only:
c           WX_ = 18 (params.h) should match the JX_spec.dat wavelengths
c           W_ = 12 (Trop-only) or 18 (std) is set in (params.h).
c       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
c
c version 6.3
c   revise cloud/aerosol OD & wavelength properties for CTM link:
c         OPTICL is new sub for cloud optical properties, but it
c              now starts with cloud OD and cloud NDX
c         OPTICA & OPTICM are new subs to convert aerosol path (g/m2) to OD
c              A is std UCI scat data 
c              M is U Michigan data tables for aerosols, includes Rel Hum effect
c
c version 6.2 corrects a long-standing problem at SZA > 89 degrees.
c   In prior versions the ray-tracing of the path (and air-mass functions)
c   back to the sun was done at the edges of the CTM layers (it was developed
c   for the grid-point J-value code at Harvard/GISS/UCI).  This left the 
c   interpolation to the mid-layer (needed for J's) open.  The prior method
c   gave irregular fluctuations in the direct solar beam at mid-layer for
c   large SZA > 88.  This is now corrected with exact ray-tracing from
c   the mid-pt of each CTM layer.  For small SZA, there is no effective 
c   difference, for large SZA, results could be erratic.
c
c   6.2 fix should be easy if you have migrated to v6.1, else some minor
c   caution may be needed:
c      replace sub SPHERE with SPHERE2, AMF2 report factors for mid and egdes.
c      replace sub OPMIE with new OPMIE, this uses the new AMF2 correctly.
c      replace sub PHOTOJ with new PHOTOJ, this just hands off AMF2 from
c            SPHERE2 to OPMIE.
c
c version 6.1 adds
c      6.1b simplifies calling sequences feeds solar factor, albedo, to PHOTOJ
c         and read LAT, LNG directly.  No substantive changes.
c      new read-in of scat data for clouds/aerosols to allow for UMich data
c      This has required substantial rewrite of some of the core subroutines:
c         OPMIE is now called for each wavelength and without aersol/cloud data
c              all subs below OPMIE are unchanged
c         OPTICD & OPTICM are new subs to convert path (g/m2) to OD and phase fn
c              D is std UCI scat data (re-ordered for clouds 1st)
c              M is U Michigan data tables for aerosols, includes Rel Hum effect
c         PHOTOJ now assembles the aerosol data (better for CTM implementation)
c      This version can reproduce earlier versions exactly, but the test input 
c         is changed from OD and NDX to PATH (g/m2) and NDX.
c version 6.0 adds
c      new 200-nm scattering data so that stratospheric aerosols can be done!
c version 5.7
c     adds the new flux diagnostics (including heating rates)
c        accurate fluxes for spherical atmos and SZA > 90 !
c     recommend geometric delta-tau factor from 1.18 to 1.12 for more accurate
c        heating rates (but more layers!)
c     tuned and corrected to be almost flux conserving (1.e-5), except
c        deep clouds, where diffusive flux is created (1.e-4)
c     still needs to return to the original 1970-code for the block-tri solution
c        after extensive profiling with F95 and 'modern' versions
c        it was found that they are much more expensive!!!
c     corrects typo in JAC(2000) fast-J paper on I+ (reflected from l.b.):
c        I+(lb) = refl/(1+refl) * (4*Integ[j(lb)*mu*dmu] + mu0*Fdirect(lb))
c version 5.6 adds
c      clean up problems with thick clouds does correct solar attenuation
c        into cloud sub-layers and into the mid-point of the CTM level
c      New calculated upward and downward FLUXES at each wavelength at TOP/BOT
c      Correct deposition of solar flux in each CTM layer (spherical)
c        awaits new diagnostics of the h's for heating rates.
c      back to old matrix solver (UCI blocksolver and matinv-4)
c version 5.5 adds
c      new code for generating and solving the block tri-diagonal scattering
c           problem.  Uses single call to GEM and general 4x4 block-tri solver.
c version 5.3c adds
c      calculates reflected UV-vis solar energy (relative to 1.0)
c      new solar spectrum (J-O2 increases in strat by 10%, J-NO by 15+%)
c
c version 5.3b changes include:
c      new data files for specral Xsection and mie-scattering.
c      add sub-layers (JXTRA) to thick cloud/aerosol layers,
c           sets up log-spaced sub-layers of increasing thickness ATAU
c      correction 'b' does massive clean up of the linking code,
c           now the only subroutine that has access to CTM arrays is PHOTOJ
c           Also, the access to the cmn_JVdat.f is 'read-only' after init.
c           This should enable safe openMP/MPI coding.
c
c common files and what they mean:
c   cmn_h.f    Main CTM header file: CTM grid, time-of-day, parameters.
c   cmn_w.f    CTM 3-D arrays (clds, aersl, O3, p, T)
c   cmn_jv.f   Xsects, Mie, etc., all for 1-D J-values, std atmospheres.
c   pam_mie.f  dimension for mie code variables.
c
c<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
c
c      photol:   Update the photolysis rate constants
c
c      PHOTOJ(UTIME,ILNG,JLAT,SZA,U0,FREFL,ODCLD,NCLDX,ZPJ)
c              Gateway to fast-JX, Update the photolysis rates
c              COMMON BLOCKS: cmn_h.f, cmn_w.f, cmn_jv.f
c
c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
c  N.B. all these need access to cmn_jv.f, but do NOT write into it.
c
c     OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
c        UCI CLOUD data sets, calculate scattering properties, 
c           scales OD at 600 nm to JX wavelength 
c              COMMON BLOCKS: cmn_jv.f
c
c     OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L)
c     OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)
c       UC Irvine & U Michigan aerosol data sets, generates fast-JX data formats
c              COMMON BLOCKS: cmn_jv.f
c
c     JRATET(PPJ,TTJ,FFF, VALJL):  Calculate J-value, called by PTOTOJ.
c              COMMON BLOCKS: cmn_h.f, cmn_jv.f
c
c     JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,DTAUX,POMEGAX,JXTRA)
c              print out atmosphere used in J-value calc.
c              >>>will be superseded by CTM routines
c              COMMON BLOCKS: params.h
c
c      FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
c
c     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
c              calc SZA and Solar Flux factor for given lat/lon/UT
c
c     SPHERE2(U0,RAD,ZHL,ZZHT,AMF2,L1_):  
c              calculate spherical geometry, air-mass factors (v 6.2)
c
c      EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
c              add sub-layers (JXTRA) to thick cloud/aerosol layers
c
c      OPMIE (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,
c    &    JXTRA,FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
c              calculate mean intensity (actinic) at each CTM levels
c              calculate fluxes and deposition (heating rates)
c              COMMON BLOCKS: cmn_h.f, cmn_jv.f, pam_mie.f
c
c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
c    SUBROUTINES:
c      MIESCT (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,ND)
c            include 'params.h' = dimension parameters
c
c      BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,M,N,ND)
c            include 'params.h' = dimension parameters
c
c      GEN (POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,B,CC,AA,A,H,C1
c            ,N,ND,ID)
c            include 'params.h' = dimension parameters
c
c      LEGND0 (X,PL,N)
c
c      GAUSSP (N,XPT,XWT)
c
c      EFOLD  (F0, F1, N, F)
c
c
c
C-----------------------------------------------------------------------
      subroutine PHOTOJ(UTIME,ILNG,JLAT,SZA,U0,FREFL
     &                 ,ODCLD,TEML,NCLDX,ZPJ,ZPJO)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
c
c  PHOTOJ is the gateway to fast-JX calculations:
c        only access to CTM 3-D GLOBAL arrays
c        sets up the 1-D column arrays for calculating J_s
c
C--------------------------------
c     AVGF   Attenuation of beam at each level for each wavelength
c     FFF    Actinic flux at each desired level
c     XQO2   Absorption cross-section of O2
c     XQO3   Absorption cross-section of O3
C--------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: LPAR
      use cmn_ctm, only: IDAY, YGRD, XGRD, ETAA, ETAB
      use cmn_fjx, only: JVN_, L1_, L2_, W_, L_, N_, X_, 
     &     JTAUMX, SZAMAX, ZZHT, RAD,
     &     DMS, DO3, TOPT, TOPM, TOP3, LPART, LPARM, LPAR3, WL, FL,
     &     QRAYL, RAA, DAA,
     &     TQQ, QO2, QO3, NRATJ, JIND, JFACTA, JINDO, JFACTAO,
     &     ATAU, ATAU0
      use cmn_met, only: P, ZOFLE, Q, SA
      use aerosols2fastjx, only: NUM_AER_TYPES, AERPATH, AER_TYPES_USED,
     &     LJV_AEROSOL
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(in)  ::  UTIME, ODCLD(LPAR+1), TEML(LPAR)
      integer,intent(in)  ::  ILNG, JLAT, NCLDX(LPAR+1)
c
      real(r8), intent(out) ::  ZPJ(LPAR,JVN_)   ! J_s indexed to CTM chemistry
      real(r8), intent(out) ::  ZPJO(LPAR)  ! J_O3(1D) indexed to CTM chemistry
      real(r8), intent(out) ::  FREFL, U0, SZA   ! fraction of energy reflected
c-----------------------------------------------------------------------

c--------key amtospheric data needed to solve plane-parallel J
      real(r8),  dimension(L1_+1) :: PPJ, TTJ,DDJ,ZZJ,ZHL
      integer, dimension(L2_+1) :: JXTRA
      
      real(r8), dimension(W_)       :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
      real(r8), dimension(L_, W_)   :: AVGF, FJFLX
      real(r8), dimension(L1_,W_)   :: DTAUX, FLXD
      real(r8), dimension(8,L1_,W_) :: POMEGAX
      real(r8), dimension(W_,L1_)   ::  FFX
      real(r8), dimension(W_,8)     ::  FFXNET
c---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real(r8)  FLXJ(L1_),FFX0,FXBOT,FABOT

      real(r8)  QL(L1_),ODABS,ODRAY
      real(r8)  RFLECT,SOLF,FREFS,FREFI
      real(r8)  AMF2(2*L1_+1,2*L1_+1)
c------------key SCATTERING arrays for clouds+aerosols------------------
      real(r8)  OPTX(5),SSAX(5),SLEGX(8,5)
      real(r8)  ODL(5,L1_),SSA(5,L1_),SLEG(8,5,L1_)
      real(r8)  OD600(L1_),  PATH,RH
c------------key arrays AFTER solving for J_s---------------------------
      real(r8)  FFF(W_,LPAR),VALJ(X_)
      real(r8)  FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)
      real(r8)  VALJL(LPAR,X_)             ! 2-D array of J_s returned by JRATET

      integer  I,J,K,L,M,KMIE,KW,NAER,NDCLD, RATIO(W_)
      real(r8)   XQO3,XQO2,DTAUC, WAVE, FLINT, TTT


c---wavelength index flags for trop-only:
c      integer, parameter, dimension(18) ::  LTWVL =
c     & [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
c---    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18

C-----------------------------------------------------------------------

      ZPJ(:,:) = 0._r8
      ZPJO(:)  = 0._r8
      FFF(:,:) = 0._r8
      FREFL = 0._r8
      FREFI = 0._r8
      FREFS = 0._r8
C-----------------------------------------------------------------------
      call SOLARZ (UTIME,IDAY,YGRD(JLAT),XGRD(ILNG), SZA,U0,SOLF)
C-----------------------------------------------------------------------
c---  SOLF = 1._r8   ! this needs to be dropped to include 6.7% annual cycle

c     write(6,'(A,I5,1PE14.5)')
c    &   ' IDAY/ Solar Flux factor:', IDAY,SOLF

c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
c                        or         99.                  80 km
      if (SZA .gt. SZAMAX) goto 99

c---load the amtospheric column data
      !// Set 1:LPAR-1, treat uppermost layer and above separately below
      do L = 1,L1_ -1
        PPJ(L) = ETAA(L) + ETAB(L)*P(ILNG,JLAT)
        TTJ(L) = TEML(L)
        DDJ(L) = DMS(L,ILNG,JLAT)
        ZZJ(L) = DO3(L,ILNG,JLAT)
        ZHL(L) = 100._r8 * ZOFLE(L,ILNG,JLAT)
        QL(L)  = Q(ILNG,JLAT,L)
      enddo

        QL(L1_)  = 0._r8
        PPJ(L1_) = ETAA(L1_)
        PPJ(L1_+1) = 0._r8
        ZHL(L1_)   = 100._r8 * ZOFLE(L1_,ILNG,JLAT)
        ZHL(L1_+1) = ZHL(L1_) + ZZHT

        !// CTM3: We do no chemistry in the uppermost layer, and should set
        !// T,Mass,O3 from climatology as well (L1_ = LPAR+1).
        !TTJ(L1_-1) = TEML(L_-1) !LPART(ILNG,JLAT)
        TTJ(L1_-1) = LPART(ILNG,JLAT)
        DDJ(L1_-1) = DMS(L1_-1,ILNG,JLAT) !should be equal to LPARM(ILNG,JLAT)
        !ZZJ(L1_-1) = DO3(L,ILNG,JLAT) !LPAR3(ILNG,JLAT)
        ZZJ(L1_-1) = LPAR3(ILNG,JLAT)

        !// Set values for T/Mass/O3 above model (L1_ = LPAR+1)
        TTJ(L1_) = TOPT(ILNG,JLAT)
        DDJ(L1_) = TOPM(ILNG,JLAT)
        ZZJ(L1_) = TOP3(ILNG,JLAT)


c---calculate spherical weighting functions (AMF: Air Mass Factor)
C-----------------------------------------------------------------------
      call SPHERE2 (U0,ZHL,AMF2)
C-----------------------------------------------------------------------

c---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8)
c---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      OD600(:) = 0._r8
      ODL(:,:) = 0._r8 
      SSA(:,:) = 0._r8
      SLEG(:,:,:) = 0._r8

      do L = 1,L1_

c---cloud in layer:  if valid cloud index (now 4:13)
         NDCLD = NCLDX(L)
c---cloud PATH here is optical depth: not g/m2
         PATH  = ODCLD(L)
         RH = QL(L)
       if (PATH .gt. 0._r8) then
         !// PATH is already given in correct units (i.e. optical depth)
         !// Test on NDCLD is done in OPTICL
         call OPTICL (OPTX,SSAX,SLEGX,  PATH,NDCLD)
         do K=1,5
            ODL(K,L)  = ODL(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
           do I=1,8
            SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
           enddo
         enddo
       endif

c---use ODL of clouds (not aerosols) at 600 nm to determine added layers
       OD600(L) = ODL(4,L)

c---aerosols in layer: check aerosol index
c---this uses data from climatology OR from current CTM (STT of aerosols)

c---FIND useful way to sum over different aerosol types!
c       do M = 1,2
c        if (M.eq.1) then
c           NAER = AER1N(ILNG,JLAT,L)
c           PATH = AER1P(ILNG,JLAT,L)
c        else
c           NAER = AER2N(ILNG,JLAT,L)
c           PATH = AER2P(ILNG,JLAT,L)
c        endif
      if (LJV_AEROSOL) then
       do M = 1, NUM_AER_TYPES
          NAER = AER_TYPES_USED(M)
          !// If aerosol type is not included, then go to next aerosol type
          !// This should not happen when looping through AER_TYPES_USED
          if (NAER.eq.0) cycle
          !// Get aerosol path [g/m2]
          PATH = AERPATH(L,M,ILNG,JLAT)

C---subroutines OPTICA & OPTICM return the same information:
C---  optical depth (OPTX), single-scat albedo (SSAX) and phase fn (SLEGX(8))
C---subs have slightly different inputs:
C---  PATH is the g/m2 in the layer, NAER in the cloud/aerosol index
C---  UMich aerosols use relative humidity (RH)

        if (PATH .gt. 0._r8) then
         if (NAER .gt.0) then
           call OPTICA (OPTX,SSAX,SLEGX,  PATH,RH, NAER)
         else
           NAER = -NAER
           call OPTICM (OPTX,SSAX,SLEGX,  PATH,RH,NAER)
         endif
         do K=1,5
            ODL(K,L) = ODL(K,L)  + OPTX(K)
            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
          do I=1,8
            SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
          enddo
         enddo
        endif
       enddo
      endif !// if (.not.LJV_AEROSOL) then

       do K=1,5
         if (ODL(K,L) .gt. 0._r8) then
            SSA(K,L) = SSA(K,L)/ODL(K,L)
          do I=1,8
           SLEG(I,K,L) = SLEG(I,K,L)/ODL(K,L)
          enddo
         endif
       enddo

      enddo !// do L = 1,L1_

c---can add aerosol OD at 600 nm to determine added layers, but not done yet
c---    OD600(L) = OD(4,L) + aerosol OD_s

c---when combining with Rayleigh and O2-O3 abs, remember the SSA and 
c---  phase fn SLEG are weighted by ODL and ODL*SSA, respectively.

c---Given the aerosol+cloud ODL/layer in visible (600 nm) calculate how to add 
C       additonal levels at top of clouds (now uses log spacing)
C-----------------------------------------------------------------------
      call EXTRAL(OD600,L1_,L2_,N_,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------

c---set surface reflectance
        RFLECT = max(0._r8,min(1._r8,SA(ILNG,JLAT)))
        do K = 1,W_
          RFL(K) = RFLECT
        enddo

C---Loop over all wavelength bins to calc mean actinic flux AVGF(L)

      do K = 1,W_

        WAVE = WL(K)
C---Pick nearest Mie wavelength to get scattering properites------------
                               KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255._r8 ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355._r8 ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500._r8 ) KMIE=4
        if( WAVE .gt. 800._r8 ) KMIE=5


c---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
c---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
        do L = 1,L1_
         TTT     = TTJ(L)
         XQO3 = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2)
     &                      ,QO3(K,1),QO3(K,2),QO3(K,3))

         XQO2 = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1)
     &                      ,QO2(K,1),QO2(K,2),QO2(K,3))

         ODABS = XQO3*ZZJ(L) + XQO2*DDJ(L)*0.20948_r8
         ODRAY = DDJ(L)*QRAYL(K)

         DTAUX(L,K) = ODL(KMIE,L) + ODABS + ODRAY

         do I=1,8
           POMEGAX(I,L,K) = SLEG(I,KMIE,L)*ODL(KMIE,L)
         enddo
           POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0_r8*ODRAY
           POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5_r8*ODRAY
         do I=1,8
           POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
         enddo
        enddo

      enddo

C-----------------------------------------------------------------------
      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
     &        AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
C-----------------------------------------------------------------------

      do K = 1,W_

c direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar=negative by convention)
c----     also at bottom (DN), does not include diffuse reflected flux.
        FLXUP(K) =  FJTOP(K)
        DIRUP(K) = -FLXD0(K)
        FLXDN(K) = -FJBOT(K)
        DIRDN(K) = -FSBOT(K)

        do L = 1,LPAR
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
        enddo
          FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
          FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
          FREFS = FREFS + SOLF*FL(K)/WL(K)

c---for each wavelength calculate the flux budget/heating rates:
c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L)]
c            but for spherical atmosphere!
c  FJFLX(L) = diffuse flux across top of layer L

c---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
c---     need special fix at top and bottom:
c---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
         FABOT = (1._r8-RFL(K))*(FJBOT(K)+FSBOT(K))
         FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
         FLXJ(1) = FJFLX(1,K) - FXBOT
       do L=2,L_
         FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
       enddo
         FLXJ(L_+1) = FJTOP(K) - FJFLX(L_,K)
c---calculate net flux deposited in each CTM layer (direct & diffuse):
         FFX0 = 0._r8
       do L=1,L1_
         FFX(K,L) = FLXD(L,K) - FLXJ(L)
         FFX0 = FFX0 + FFX(K,L)
       enddo

c  NB: the radiation level ABOVE the top CTM level is included in these budgets
c      these are the flux budget/heating terms for the column:
c  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
c  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
c  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
c  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
c  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
c  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
c       these are surface fluxes to compare direct vs. diffuse:
c  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
c  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

      FFXNET(K,1) = FLXD0(K)
      FFXNET(K,2) = FSBOT(K)
      FFXNET(K,3) = FLXD0(K) + FSBOT(K)
      FFXNET(K,4) = FJTOP(K)
      FFXNET(K,5) = FFX0
      FFXNET(K,6) = FABOT
      FFXNET(K,7) = FSBOT(K)
      FFXNET(K,8) = FJBOT(K)

c-----------------------------------------------------------------------
      enddo       ! end loop over wavelength K

          FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
          FREFI = FREFI/FREFS

c---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

C-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJL)
C-----------------------------------------------------------------------

c---map the J-values from fast-JX onto CTM ones (use JIND & JFACTA)
      do L = 1,LPAR
        do J = 1,NRATJ
          if (JIND(J).gt.0) then 
            ZPJ(L,J) = VALJL(L,JIND(J))*JFACTA(J)
          else
            ZPJ(L,J) = 0._r8
          endif
        enddo
        if (JINDO .gt. 0) then
          ZPJO(L) = VALJL(L,JINDO)*JFACTAO
        else
          ZPJO(L) = 0._r8
        endif
      enddo

c---diagnostics that are NOT returned to the CTM code

C-----------------------------------------------------------------------
c     write(6,*)'fast-JX-(6.1)----PHOTOJ internal print: Atmosphere----'
c---used last called values of DTAUX and POMEGAX, should be 600 nm

c      call JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,
c     &            DTAUX(1,W_),POMEGAX(1,1,W_),JXTRA)

C---PRINT SUMMARY of mean intensity, flux, heating rates:
c     write(6,*)
c     write(6,*)'fast-JX(6.1)----PHOTOJ internal print: Mean Intens----'
c     write(6,'(a,5f10.4)')
c    & ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',
c    &  RFLECT,SZA,U0,FREFI,FREFL

c     write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
c     write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
c     write(6,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin'
c     do L = LPAR,1,-1
c      do K=NW1,NW2
c       RATIO(K) = (1.d5*FFF(K,L)/FL(K))
c      enddo
c       write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
c     enddo

c     write(6,*)
c     write(6,*)'fast-JX(6.1)----PHOTOJ internal print: Net Fluxes----'
c     write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
c     write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
c     write(6,*) ' ---NET FLUXES--- '
c     write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
c     write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
c     write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
c     write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
c     write(6,*) ' ---SRF FLUXES--- '
c     write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
c     write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)
c     write(6,'(2a)') '  ---NET ABS per layer:       10000=Fsolar',
c    & '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
c     do L = LPAR,1,-1
c      do K=NW1,NW2
c       RATIO(K) = 1.d5*FFX(K,L)
c      enddo
c       write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
c     enddo

C-----------------------------------------------------------------------
      
   99 continue
      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<end CTM-fastJX linking subroutines<<<<<<<<<<<<<<



c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<

c------------------------------------------------------------------------------
      subroutine OPTICL (OPTD,SSALB,SLEG, RDCLD,NDCLD)
c------------------------------------------------------------------------------
cc     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
c---set CLOUD fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
c---  this now requires the cloud layer OD rather than water path
c---  assumed CLOUD Optical Depth given with type for optical properties (v-6.3)
c
c 04 W_H01 (H1/Deir)GAMMA:r-m=0.1/alf=2 n=1.335   reff=0.250___G=.0940_rho=1.000
c 05 W_H04 (H1/Deir)GAMMA:r-m=0.4/alf=2 n=1.335   reff=1.000___G=1.508_rho=1.000
c 06 W_H40 (H1/Deir)GAMMA:r-m=4.0/alf=2 n=1.335   reff=10.00___G=146.4_rho=1.000
c 07 W_C02 (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
c 08 W_C04 (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
c 09 W_C08 (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
c 10 W_C13 (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
c 11 W_L06 (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
c 12 Ice-Hexagonal (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
c 13 Ice-Irregular (Mishchencko)                  reff=50.00___G=999.9_rho=0.917
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: QAA, SAA, PAA
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(out)::    OPTD(5)    ! optical depth of layer
      real(r8), intent(out)::    SSALB(5)   ! single-scattering albedo
      real(r8), intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real(r8), intent(in)::     RDCLD      ! optical depth of cloud layer @600 nm
      integer,intent(inout)::  NDCLD      ! index of cloud layer:  4:13

      integer I,J
      real(r8)  XTINCT, REFF,RHO

c---default cloud type C1, Reff = 12 microns
      if (NDCLD .gt. 13 .or. NDCLD .lt.4) then
         NDCLD = 9
      endif

c--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = RDCLD * QAA(J,NDCLD)/QAA(4,NDCLD)
         SSALB(J) = SAA(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PAA(I,J,NDCLD)
        enddo
      enddo

      return
      end


c------------------------------------------------------------------------------
      subroutine OPTICA (OPTD,SSALB,SLEG, PATH,RELH,L)
c------------------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
c---for the UCI aerosol data sets, calculates optical properties at fast-JX's
c              std 5 wavelengths:200-300-400-600-999nm
c---UCI aersols optical data  v-6.1:
c
c   14 S-Bkg LOGN:r=.090 s=.600 reff=.221  n=1.514/1.473/1.459/1.448/1.435  
c   15 S-Vol LOGN:r=.080 s=.800 reff=.386  n=1.514/1.473/1.459/1.448/1.435  
c   16 UT-sulfate LOGN:r=0.05 s=.693 n=1.44           reff=0.166___rho=1.769
c   17 UT-sulfate LOGN:r=0.05 s=.693 n=1.46           reff=0.166___rho=1.769
c   18 UM-SULFate LOGN:r=.050 s=.642 n=1.53           reff=0.140___rho=1.769
c   19 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i     reff=0.140___rho=1.500
c   20 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i     reff=0.150___rho=1.500
c   21 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i  reff=0.149___rho=1.230
c   22 UM-FF04(%BC) LOGN:r=.050 s=.642 n=1.541+0.02i  reff=0.140___rho=1.212
c   23 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i  reff=0.140___rho=1.230
c   24 Mdust .15 (R.V. Martin generated phase fns)
c   25 Mdust .25 (R.V. Martin generated phase fns)
c   26 Mdust 0.4 (R.V. Martin generated phase fns)
c   27 Mdust 0.8 (R.V. Martin generated phase fns)
c   28 Mdust 1.5 (R.V. Martin generated phase fns)
c   29 Mdust 2.5 (R.V. Martin generated phase fns)
c   30 Mdust 4.0 (R.V. Martin generated phase fns)

c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: NAA, QAA, SAA, PAA, RAA, DAA
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(out)::    OPTD(5)    ! optical depth of layer
      real(r8), intent(out)::    SSALB(5)   ! single-scattering albedo
      real(r8), intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real(r8), intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real(r8), intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     L       ! index of cloud/aerosols

      integer I,J
      real(r8)  XTINCT, REFF,RHO

      !if (L .gt. NAA .or. L .lt. 4) then
      if (L .gt. NAA .or. L .lt. 3) then
         write(6,*) ' aerosol index out-of-range: L/NAA',L,NAA
         L = 18
      endif

      !if (L .lt. 14) then
      if (L.ne.3 .and. L .lt. 14) then
         write(6,*) ' aerosol as cloud: L/NAA',L,NAA
      endif

         REFF = RAA(L)
         RHO = DAA(L)
      do J=1,5
c---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75_r8*QAA(J,L)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,L)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,L)
       enddo
      enddo

      return
      end


c------------------------------------------------------------------------------
      subroutine OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)
c------------------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c     Modified for CTM3 with increased number of datasets (N_UMSET)
C-----------------------------------------------------------------------
c   UC Irvine & U Michigan aerosol data sets, generates fast-JX data formats
c---NB Approximates the Legendre expansion(L) of the scattering phase fn 
c---           as (2*L+1)*g**L
c---UMAER(I,J,K,L):
c   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
c   J=1:6 = [200, 300, 400, 550, 600 , 1000 nm]
c   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
c   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC), 
c                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]

c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: n_umset, UMAER
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(out)::    OPTD(5)    ! optical depth of layer
      real(r8), intent(out)::    SSALB(5)   ! single-scattering albedo
      real(r8), intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
      real(r8), intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real(r8), intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     L       ! index of cloud/aerosols

      integer KR,J
      real(r8)  R,FRH, GCOS, XTINCT

c---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
c---interpolate in Relative Humidity
c---extrapolate pahse fn from first term (g)

      if (L .gt. n_umset) then
         write(6,*) ' UM aer index too large: L',L
         L = 1
      endif

      R = 100._r8*min(1._r8, max(.01_r8, RELH))
         KR = (R/5._r8)
         KR = max(0, min(19, KR)) + 1
        if (KR.lt.20) then
         FRH = 0.20_r8*(R - 5._r8*real(KR-1,r8))
        else
         FRH = 0.25_r8*(R - 5._r8*real(KR-1,r8))
        endif

      do J=1,5
       SSALB(J) = UMAER(1,J,KR,L)*(1._r8-FRH) + UMAER(1,J,KR+1,L)*FRH

       XTINCT = UMAER(3,J,KR,L)*(1._r8-FRH) + UMAER(3,J,KR+1,L)*FRH
       OPTD(J) = PATH*XTINCT
       
       GCOS   = UMAER(2,J,KR,L)*(1._r8-FRH) + UMAER(2,J,KR+1,L)*FRH
       SLEG(1,J) =  1._r8
       SLEG(2,J) =  3._r8*GCOS
       SLEG(3,J) =  5._r8*GCOS**2
       SLEG(4,J) =  7._r8*GCOS**3
       SLEG(5,J) =  9._r8*GCOS**4
       SLEG(6,J) = 11._r8*GCOS**5
       SLEG(7,J) = 13._r8*GCOS**6
       SLEG(8,J) = 15._r8*GCOS**7
  
      enddo
      return
      end


c-----------------------------------------------------------------------
      subroutine SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c
c     GMTIME = UT for when J-values are wanted 
c           (for implicit solver this is at the end of the time step)
c     NDAY   = integer day of the year (used for solar lat and declin)
c     YGRDJ  = laitude (radians) for grid (I,J)
c     XGDRI  = longitude (radians) for grid (I,J)
c
c     SZA = solar zenith angle in degrees
c     COSSZA = U0 = cos(SZA)
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      real(r8), intent(in) ::   GMTIME,YGRDJ,XGRDI
      integer, intent(in) ::  NDAY
      real(r8), intent(out) ::  SZA,COSSZA,SOLFX

      real(r8)  PI, PI180, LOCT
      real(r8)  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
c
      PI     = 3.141592653589793_r8
      PI180  = PI/180._r8
      SINDEC = 0.3978_r8*sin(0.9863_r8*(real(NDAY,r8)-80._r8)*PI180)
      SOLDEK = asin(SINDEC)
      COSDEC = cos(SOLDEK)
      SINLAT = sin(YGRDJ)
      SOLLAT = asin(SINLAT)
      COSLAT = cos(SOLLAT)
c
      LOCT   = (((GMTIME)*15._r8)-180._r8)*PI180 + XGRDI
      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/PI180

c     write(6,*) ' XGRDI,YGRDJ',XGRDI,YGRDJ
c     write(6,*) ' LOCT (rad)',LOCT
c     write(6,*) ' SINDEC,COSDEC', SINDEC,COSDEC
c     write(6,*) ' SINLAT,COSLAT', SINLAT,COSLAT
c     write(6,*) ' COS, SZA',COSSZA,SZA


      SOLFX  = 1._r8-(0.034_r8*cos(real(NDAY-186,r8)*2._r8*PI/365._r8))

      return
      end


c-----------------------------------------------------------------------
      subroutine SPHERE2(GMU,ZHL,AMF2)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c
c----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
c  This new AMF2 does each of the half-layers of the CTM separately,
c     whereas the original, based on the pratmo code did the whole layers
c     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
c  Since fast-JX is meant to calculate the intensity at the mid-layer, the
c     solar beam at low sun (interpolated between layer edges) was incorrect.
c  This new model does make some approximations of the geometry of the layers:
c     the CTM layer is split evenly in mass (good) and in height (approx).
c
c  Calculation of spherical geometry; derive tangent heights, slant path
c  lengths and air mass factor for each layer. Not called when
c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
c  beam (where tangent height is below altitude J-value desired at).
C-----------------------------------------------------------------------
c in:
c     GMU     = MU0 = cos(solar zenith angle)
c     RAD     radius of Earth mean sea level (cm)
c     ZHL(L)  height (cm) of the bottome edge of CTM level L
c     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1)
c     L1_     dimension of CTM = levels +1 (L+1 = above-CTM level)
c out:
c     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: L1_, RAD
      implicit none
c-----------------------------------------------------------------------
      real(r8), intent(in)  ::   GMU,ZHL(L1_+1)
      real(r8), intent(out) ::   AMF2(2*L1_+1,2*L1_+1)

c     RZ      Distance from centre of Earth to each point (cm)
c     RQ      Square of radius ratios
c     SHADHT  Shadow height for the current SZA
c     XL      Slant path between points

      integer  I, J, K, II, L2
      real(r8)   XMU1,XMU2,XL,DIFF,SHADHT,RZ(L1_+1)
      real(r8)   RZ2(2*L1_+1),RQ2(2*L1_+1)
c
c--- must have top-of-atmos (NOT top-of-CTM) defined
c      ZHL(L1_+1) = ZHL(L1_) + ZZHT

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1_+1
        RZ(II)   = RAD + ZHL(II)
      enddo

c---calculate heights for edges of split CTM-layers
      L2 = 2*L1_
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5_r8*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1_+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      enddo


c---shadow height for SZA > 90
      if (GMU .lt. 0.0_r8)  then
        SHADHT = RZ2(1)/dsqrt(1.0_r8-GMU**2)
      else
        SHADHT = 0._r8
      endif

c---up from the surface calculating the slant paths between each level
c---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0._r8

      do 16 J = 1,2*L1_+1

c  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
c  Ascend from layer J calculating AMF2s
        XMU1 = abs(GMU)
        do I = J,2*L1_
          XMU2     = dsqrt(1.0_r8 - RQ2(I)*(1.0_r8-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
c--fix above top-of-atmos (L=L1_+1), must set DTAU(L1_+1)=0
          AMF2(2*L1_+1,J) = 1._r8
c
c  Twilight case - Emergent Beam, calc air mass factors below layer
        if (GMU .ge. 0.0_r8) goto 16

c  Descend from layer J 
          XMU1       = abs(GMU)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0_r8-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0._r8)   ! filter
c  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0_r8)  then
            XMU2      = sqrt(1.0_r8 - (1.0_r8-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2._r8*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
c  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0_r8
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue
      return
      end
c
c
C-----------------------------------------------------------------------
      subroutine EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c
c
c    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
c    this version sets up log-spaced sub-layers of increasing thickness ATAU
c
c     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
c        This can be just cloud or cloud+aerosol, it is used only to set
c        the number in levels to insert in each layer L
c        Set for log-spacing of tau levels, increasing top-down.
c
c     N.B. the TTAU, etc calculated here are NOT used elsewhere

c---The log-spacing parameters have been tested for convergence and chosen
c---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
c---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100 
c---  ATAU = 1.12 now recommended for more -accurate heating rates (not J_s)
C-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: L_
      implicit none
c-----------------------------------------------------------------------
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real(r8),  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real(r8),  intent(in) ::  ATAU,ATAU0
      integer, intent(out) ::  JXTRA(L2X+1)   !number of sub-layers to be added
c
      integer JTOTL,I,L,L2 
      real(r8)  TTAU(2*L_+3),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
c
C---Reinitialize arrays
      TTAU(:)  = 0._r8
      JXTRA(:) = 0
c
c---combine these edge- and mid-layer points into grid of size:
c---              L2X+1 = 2*L1X+1 = 2*LPAR+3
c---calculate column optical depths above each level, TTAU(1:L2X+1)
c---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
c
c---Divide thick layers to achieve better accuracy in the scattering code
c---In the original fast-J, equal sub-layers were chosen, this is wasteful
c---and this new code (ver 5.3) uses log-scale:  
c---        Each succesive layer (down) increase thickness by ATAU > 1
c---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
c---        4 sub-layers with ODs = 1 - 2 - 4 - 8
c---The key parameters are:
c---        ATAU = factor increase from one layer to the next
c---        ATAUMN = the smallest OD layer desired
c---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
c---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1._r8
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0_r8
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5_r8 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
c---Now compute the number of log-spaced sub-layers to be added in
c---   the interval TTAU(L2) > TTAU(L2+1)
c---The objective is to have successive TAU-layers increasing by factor ATAU >1
c---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5_r8)))
        endif
      enddo

c---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      return
      end


c-----------------------------------------------------------------------
      real(r8) FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
c-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c
c  Three-point linear interpolation function
c-----------------------------------------------------------------------

      use cmn_precision, only: r8
      implicit none
      real(r8)   TINT,T1,T2,T3,F1,F2,F3
      if (TINT .le. T2)  then
        if (TINT .le. T1)  then
          FLINT = F1
        else
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        endif
      else
        if (TINT .ge. T3)  then
          FLINT = F3
        else
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        endif
      endif
      return
      end


c-----------------------------------------------------------------------
      subroutine JRATET(PPJ,TTJ,FFF, VALJL)
c-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c
c in:
c        PPJ(L1_+1) = pressure profile at edges
c        TTJ(L1_) = = temperatures at mid-level
c        FFF(K=1:NW, L=1:LPAR) = mean actinic flux 
c out:
c        VALJ(LPAR,X_)  LPAR = no of levels
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: LPAR
      use cmn_fjx, only: L1_, W_, X_, TQQ, QO2, QO3, Q1D, NJVAL, QQQ
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(in)  ::  PPJ(L1_+1),TTJ(L1_)
      real(r8), intent(inout)  ::  FFF(W_,LPAR)
      real(r8), intent(out) ::  VALJL(LPAR,X_)

      real(r8)  FLINT             ! external function for X-sections
      real(r8)  VALJ(X_)          ! temp for call J_s at one L
      real(r8)  QO2TOT(W_), QO3TOT(W_), QO31DY(W_), QO31D, QQQT, TFACT
      real(r8)  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,
     &     QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
c
      do L = 1,LPAR    ! master loop over layer = L

c---need temperature and density (for some quantum yields):
c---in this case the Pressures PPJ are defined at the boundaries,
c---                Temperatures in the middle of each layer
          TT   = TTJ(L)
         if (L .eq. 1) then
          PP = PPJ(1)
         else
          PP  = (PPJ(L)+PPJ(L+1))*0.5_r8
         endif
          DD = 7.24e18_r8*PP/TT

c---if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is too weak
c        and does not represent the decay of 215.5-221.5 nm sunlight.
        if (PP .gt. 100._r8) then
          if (W_ .eq. 18) then
            J=11
            FFF(J,L) = 0._r8
          elseif (W_ .eq. 12) then
            J=5
            FFF(J,L) = 0._r8
          endif
        endif

        do J = 1,NJVAL
          VALJ(J) = 0._r8
        enddo

c First treat O2+hv, total O3+hv and fraction of O3+hv to O1D
c     i.e. QO2TOT, QO3TOT and QO31D
        if (TT .le. TQQ(2,1))  then
          if (TT .le. TQQ(1,1))  then
            TFACT = 0._r8
          else
            TFACT = (TT -TQQ(1,1))/(TQQ(2,1) -TQQ(1,1))
          endif
          do K = 1,W_
           QO2TOT(K) = QO2(K,1) + (QO2(K,2) - QO2(K,1))*TFACT
          enddo
        else
          if (TT .ge. TQQ(3,1))  then
            TFACT = 1._r8
          else
            TFACT = (TT -TQQ(2,1))/(TQQ(3,1) -TQQ(2,1))
          endif
          do K = 1,W_
           QO2TOT(K) = QO2(K,2) + (QO2(K,3) - QO2(K,2))*TFACT
          enddo
        endif
  
        if (TT .le. TQQ(2,2))  then
          if (TT .le. TQQ(1,2))  then
            TFACT = 0._r8
          else
            TFACT = (TT -TQQ(1,2))/(TQQ(2,2) -TQQ(1,2))
          endif
          do K = 1,W_
           QO3TOT(K) = QO3(K,1) + (QO3(K,2) - QO3(K,1))*TFACT
          enddo
        else
          if (TT .ge. TQQ(3,2))  then
            TFACT = 1._r8
          else
            TFACT = (TT -TQQ(2,2))/(TQQ(3,2) -TQQ(2,2))
          endif
          do K = 1,W_
           QO3TOT(K) = QO3(K,2) + (QO3(K,3) - QO3(K,2))*TFACT
          enddo
        endif
  
        if (TT .le. TQQ(2,3))  then
          if (TT .le. TQQ(1,3))  then
            TFACT = 0._r8
          else
            TFACT = (TT -TQQ(1,3))/(TQQ(2,3) -TQQ(1,3))
          endif
          do K = 1,W_
           QO31DY(K) = Q1D(K,1) + (Q1D(K,2) - Q1D(K,1))*TFACT
          enddo
        else
          if (TT .ge. TQQ(3,3))  then
            TFACT = 1._r8
          else
            TFACT = (TT -TQQ(2,3))/(TQQ(3,3) -TQQ(2,3))
          endif
          do K = 1,W_
           QO31DY(K) = Q1D(K,2) + (Q1D(K,3) - Q1D(K,2))*TFACT
          enddo
        endif

C     Set the first 3 J-values
        do K = 1,W_
           QO31D  = QO31DY(K)*QO3TOT(K)
          VALJ(1) = VALJ(1) + QO2TOT(K)*FFF(K,L)
          !// CTM3 does not use total O3; it can be separated here
          !// (subtracting QO31D*FFF(K,L) in VAL(2)), or in the calling routine.
          VALJ(2) = VALJ(2) + QO3TOT(K)*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo


C     Now set the rest of the J-values
        do J = 4,NJVAL

          if (TQQ(2,J) .gt. TQQ(1,J)) then
           TFACT = max(0._r8,
     &            min(1._r8,(TT-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J))))
          else
           TFACT = 0._r8
          endif

          do K = 1,W_
            QQQT    = QQQ(K,1,J) + (QQQ(K,2,J) - QQQ(K,1,J))*TFACT
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
          enddo

        enddo

c>>>>   all these Stern-Volmer q-ylds incoprporated into the T interpolation
c>>>>Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])
c>>>>Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])
c>>>>Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])
c>>>>Acetone is a special case:   (as per Blitz et al GRL, 2004)

c----Load array of J-values in native order, need to be indexed/scaled
c    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
        do J=1,NJVAL
          VALJL(L,J) = VALJ(J)
        enddo

      enddo    ! master loop over L=1,LPAR
      return
      end


c-----------------------------------------------------------------------
      subroutine JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,DTAUX,POMEGAX,JXTRA)
c-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: L1_, L2_, ZZHT
      implicit none
c-----------------------------------------------------------------------
c--------key amtospheric data needed to solve plane-parallel J----------
      real(r8), intent(in), dimension(L1_+1) :: TTJ,DDJ,ZZJ,ZHL
      real(r8), intent(in), dimension(L1_+1) :: PPJ
      integer,intent(in), dimension(L2_+1) :: JXTRA
      real(r8), intent(in)                 :: DTAUX(L1_),POMEGAX(8,L1_)
c-----------------------------------------------------------------------
      integer  I,J,K,L
      real(r8)   COLO2,COLO3,ZKM,DELZ,ZTOP

      write(6,'(4a)') '   L z(km)     p      T   ',
     & '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',
     & '  g(cos) CTM lyr=>'

          COLO2 = 0._r8
          COLO3 = 0._r8
          ZTOP = ZHL(L1_) + ZZHT

        do L = L1_,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948_r8  
          COLO3 = COLO3 + ZZJ(L)
          DELZ = ZTOP-ZHL(L)
          ZTOP = ZHL(L)
          ZKM = ZHL(L)*1.e-5_r8

      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')
     &      L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ,
     &      COLO2,COLO3,DTAUX(L),POMEGAX(1,L),POMEGAX(2,L)/3._r8,
     &      JXTRA(L+L),JXTRA(L+L-1)

        enddo

      return
      end


C-----------------------------------------------------------------------
      subroutine OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
     &                  FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: L_, L1_, W_, L2_, M2_, N_, ATAU0
      implicit none
c-----------------------------------------------------------------------

      real(r8), intent(in) ::   DTAUX(L1_,W_),POMEGAX(8,L1_,W_)
      real(r8), intent(in) ::   AMF2(2*L1_+1,2*L1_+1)
      real(r8), intent(in) ::   U0,RFL(W_)
      integer, intent(in) ::  JXTRA(L2_+1)
      real(r8), intent(out) :: 
     &     FJACT(L_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_)
      real(r8), intent(out) ::  FJFLX(L_,W_),FLXD(L1_,W_),FLXD0(W_)
c
      integer JNDLEV(L_),JNELEV(L1_)
      integer JADDLV(L2_+1),JADDTO(L2_+1),L2LEV(L2_+1)
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
      integer LZ0,LZ1,LZMID
      real(r8)   SUMT,SUMJ

      real(r8)  DTAU(L1_+1,W_),POMEGAJ(M2_,L2_+1,W_),TTAU(L2_+1,W_)
      real(r8)  FTAU2(L2_+1,W_),POMEGAB(M2_,W_)
      real(r8)  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
      real(r8), dimension(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX
c--- variables used in mie code-----------------------------------------
      real(r8), dimension(W_)         :: FJT,FJB 
      real(r8), dimension(N_,W_)      :: FJ,FZ,ZTAU
      real(r8), dimension(M2_,N_,W_) :: POMEGA
      real(r8), dimension(2*L1_,W_)   :: FLXD2

C  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
c in:    
c     DTAUX(1:L1_,1:W_) = optical depth of each layer
c     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s abledo)
c     U0  = cos (SZA)
c     RFL(1:W_) = Lambertian albedo of surface
c     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
c        AMF2 now does both edges and middle of CTM layers
c     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
c out:
c     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
c  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
c     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
c     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
c     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
c     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
C        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
c     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
c        this should take into account sphericity, and is not just = mu0
c     FLXD0(1:W_) = sum of solar flux deposited in atmos
c        does NOT include flux on lower surface, does NOT mean absorbed!
C-----------------------------------------------------------------------
c
c     DTAU     Local optical depth of each CTM level
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU2     Attenuation of solar beam
c     POMEGAJ  Scattering phase function
c
c---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the 
c   factor increase from sub-layer to sub-layer 
c
C---------------------SET UP FOR MIE CODE-------------------------------

c-----------------wavelength independent--------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c
c  The following mapping holds for JADDLV=0
c        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
c        Top:       TTAU(L2_)  ==> ZTAU(3)
c        Infinity:     0.0     ==> ZTAU(1)
c        index: 2*(L2_+1-L2)+1 ==> LZ
c
c  Mie scattering code only used from surface to level L2_
C------------------------------------------------------------------------
c
C------------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
c  to be incremented linearly, and the flux fz to be attenuated top-down 
c    (avoiding problems where lower level fluxes are zero).
C------------------------------------------------------------------------
c
c  Ascend through atmosphere transposing grid and adding extra points
c  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
c  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
c    because we need to insert the intermediate layers (even LZ) for the 
c    asymmetric scattering code.


c  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
c    order, expanded, doubled-level scatter grid. 
c    Note that we need to deal with the expansion by JADD levels (L2L).
c      These JADDLV levels are skipped and need to be interpolated later.
c    Note that only odd LZ levels are filled, 

C----------------------re-grid data---------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
c           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
c     JADDLV(L2)  Number of new levels actually added at each wavelength
c            where JADDLV = 0 when there is effectively no FTAU2 
c     JADDTO(L2)   Total number of new levels to add to and above level (L2)
c     JNDLEV(L) = L2 index that maps on CTM mid-layer L
c
c---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
c---    JADDLV is taken from JXTRA, which is based on visible OD.
c---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
c---these should be fixed for all wavelengths to lock-in the array sizes
      do L2 = 1,L2_,1
        JADDLV(L2) = JXTRA(L2)
      enddo
        JADDTO(L2_+1) = 0
      do L2 = L2_,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

c---expanded grid now included CTM edge and mid layers plus expanded 
c---    grid to allow for finer delta-tau at tops of clouds.
c---    DIM of new grid = L2_ + JADDTO(1) + 1

c---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
c     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2_+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

c---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
c---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,L_
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(L_+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2_ + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        write(6,'(a,2i9)') ' overflow of scatter arrays:',ND,N_
        stop
      endif

c----------------begin wavelength dependent set up------------------------------

C---Reinitialize arrays
      ZTAU(:,:)     = 0._r8
      FZ(:,:)       = 0._r8
      POMEGA(:,:,:) = 0._r8

      do K=1,W_

C---Set up optical depth DTAU(L)
       do L = 1,L1_
        DTAU(L,K) = DTAUX(L,K)
       enddo
        DTAU(L1_+1,K) = 0._r8

c---Define the total scattering phase fn for each CTM layer L=1:L_+1
c---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
C---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1_
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)
        enddo
       enddo

C---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
c---      at the middle & edges of the CTM layers L=1:2*L1_+1
c---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
c---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0._r8
        FTAU2(L2_+1,:) = 1.0_r8
       do LL = 1,2*L1_+1
         L = (LL+1)/2
        if (AMF2(LL,LL) .gt. 0.0_r8) then
           XLTAU = 0.0_r8
         do II = 1,2*L1_+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5_r8*DTAU(I,K)*AMF2(II,LL)
         enddo
         if (XLTAU .lt. 76._r8) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         endif
        endif
       enddo

c---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
c---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0._r8
       do LL = 1,2*L1_
        if (AMF2(LL,LL) .gt. 0._r8) then 
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        endif
       enddo
        if (AMF2(1,1) .gt. 0._r8) then 
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0._r8
        endif

       do LL = 2,2*L1_,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

c---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
c---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0._r8
       if (AMF2(2*L1_,2*L1_) .gt. 0._r8) then
        do L=1,L1_
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        enddo
       endif

C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the JNDLEV array.
C------------------------------------------------------------------------
c---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
c---    points (1:L_) plus 1 for the mid point of added top layer.
c---combine these edge- and mid-layer points into grid of size:
c---              L2_+1 = 2*L1_+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2_+1)
c---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD

        TTAU(L2_+1,K) = 0.0_r8
       do L2 = L2_,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5_r8 * DTAU(L,K)
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

c----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0._r8) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1._r8+RFL(K))
       else
        ZFLUX(K) = 0._r8
       endif

c  Calculate scattering properties, level centres then level boundaries
c ***be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index***
       do L2 = L2_,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        enddo
       enddo
c---lower boundary value is set (POMEGAJ(I,1), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2_+1,K) = POMEGAJ(I,L2_,K)
       enddo
c---now have POMEGAJ filled at even points from L2=3:L2_-1
c---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2_-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + 
     &           POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        enddo
       enddo

C---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
c---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        enddo
       enddo

c   Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2_             ! L2 = index of CTM edge- and mid-layers
         L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
         LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
         L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

        if (L22 .gt. 0) then
          TAUBTM(K) = TTAU(L2,K)
          TAUTOP(K) = TTAU(L2+1,K)
          FBTM(K)   = FTAU2(L2,K)
          FTOP(K)   = FTAU2(L2+1,K)
         do I = 1,M2_
          POMEGAB(I,K) = POMEGAJ(I,L2,K)
         enddo

c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
         ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

c---fraction from TAUBTM=>TAUTOP
          ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
c---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg), 
c---else scale by TAU
          if (U0 .gt. 0.02_r8) then
            FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
          else
            if (FBTM(K) .lt. 1.e-32_r8) then
              FZ(LZZ,K) = 0._r8
            else    
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            endif
          endif
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) + 
     &               ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          enddo
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          enddo
         enddo
        endif
       enddo

c   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5_r8*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5_r8*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        enddo
       enddo
      
      enddo  ! wavelength loop!

C-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
C-----------------------------------------------------------------------

c---Move mean intensity from scatter array FJ(LZ=1:ND) 
c---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_

C FJX 6.6: Use mean intensity averaged throughout layer instead of mid-layer

       if (.false.) then
c---mean intensity:  4*<I> + solar at mid-layer
        do L = 1,L_
          L2L = JNDLEV(L)
          LZ  = ND+2 - 2*L2L
          FJACT(L,K) = 4._r8*FJ(LZ,K) + FZ(LZ,K)
        enddo

c---mean diffuse flux: after if-block; calculations are the same for both
c---                   methods.

      else
c---mean intensity averaged throughout layer:
c---   average (tau-weighted) the odd-points from: NELEV(L-1) to NELEV(L)
       do L = 1,L_
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        endif
         SUMJ = (4._r8*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K))
     &         +(4._r8*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ
     &          +(4._r8*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        enddo
         FJACT(L,K) = SUMJ/SUMT
       enddo

      endif !// mean intensity

c---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
c---      average (tau-wtd) the h's just above and below the L-edge
        do L = 1,L_
          L2L = JNELEV(L)
          LZ  = ND+2 - 2*L2L
          FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
          FJFLX(L,K)=
     &         4._r8*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1._r8-FJFLX0))
       enddo


c---diffuse fluxes reflected at top, incident at bottom 
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)

      enddo  ! wavelength loop!

      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<end core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<


c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
C-----------------------------------------------------------------------
      subroutine MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: M2_, N_, W_, M_, EMU
      implicit none
c-----------------------------------------------------------------------

c--- expect parameters M_, N_ in params.h ------------------------------
c
      integer, intent(in) ::  ND
      real(r8), intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_)
     &                       ,RFL(W_),U0,ZFLUX(W_)
      real(r8), intent(out) ::  FJ(N_,W_),FJT(W_),FJB(W_)

      real(r8)  PM(M_,M2_),PM0(M2_)
      integer I, IM
C-----------------------------------------------------------------------
C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol_n of inhomogeneous Rayleigh scattering atmosphere. 
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
C-----------------------------------------------------------------------
      do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
         PM(I,IM) = PM0(IM)
       enddo
      enddo

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25_r8*PM0(IM)
       enddo

c---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

      return
      end


C-----------------------------------------------------------------------
      subroutine LEGND0 (X,PL,N)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
C---Calculates ORDINARY Legendre fns of X (real) 
C---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      use cmn_precision, only: r8
      implicit none
      integer, intent(in) :: N
      real(r8), intent(in)  :: X
      real(r8), intent(out) :: PL(N)
c
      integer I
      real(r8)  DEN
C---Always does PL(2) = P[1]
        PL(1) = 1._r8
        PL(2) = X
        do I = 3,N
         DEN = real(I-1,r8)
         PL(I) = PL(I-1)*X*(2._r8-1.0/DEN) - PL(I-2)*(1._r8-1._r8/DEN)
        enddo
      return
      end


C-----------------------------------------------------------------------
      subroutine BLKSLV
     &     (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
C  Sets up and solves the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
c  This goes back to the old, dumb, fast version 5.3
c
c Updated for Oslo CTM3 18 Aug 2010, when MP was in Oslo.
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: M2_, N_, W_, M_, WT, EMU
      implicit none
c-----------------------------------------------------------------------
c--- expect parameters M_, N_ in parm_MIE.f------------------------------

      integer, intent(in) ::  ND
      real(r8), intent(in)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_)
     &                       ,PM(M_,M2_),PM0(M2_)
     &                       ,RFL(W_),ZFLUX(W_)
      real(r8), intent(out) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_)

      real(r8), dimension(M_,N_,W_)    ::  A,C,H,   RR

      real(r8), dimension(M_,M_,N_,W_) ::  B,AA,CC,  DD
      real(r8), dimension(M_,M_) ::  E
      real(r8)  SUMB,SUMBX,SUMT
      integer I, J, K, L


      do K = 1,W_
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K),
     &     PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),
     &             A(1,1,K),H(1,1,K),C(1,1,K), ND)
      enddo

      do K = 1,W_
C-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1._r8/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1._r8/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1._r8/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1._r8/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K)
     &                  -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K)
     &              +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1._r8/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1._r8/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1._r8/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1._r8/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K)
     &              + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

C---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K)
     &     + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K)
     &     + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K)
     &     - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K)
     &     - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1._r8/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1._r8/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1._r8/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1._r8/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K)
     &              +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

C-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K)
     &    + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K)
     &    + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

C----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2)
     &          + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2)
     &          + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

c---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
c---FJBOT = scaled diffuse flux onto surface:
c---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
c---SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2)
     &      + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2)
     &      + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4._r8*SUMB*RFL(K)/(1.0_r8 + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4._r8*SUMT
       FJBOT(K) = 4._r8*SUMB - SUMBX

      enddo

      return
      end


C-----------------------------------------------------------------------
      subroutine GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0
     &              ,B,CC,AA,A,H,C,  ND)
C-----------------------------------------------------------------------
c     Version FJX6.7 (checked December 2012)
C-----------------------------------------------------------------------
C  Generates coefficient matrices for the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_fjx, only: M2_, N_, M_, EMU, WT
      implicit none
c-----------------------------------------------------------------------
c--- expect parameters M_, N_ in params.h ------------------------------
      integer, intent(in) ::  ND
      real(r8), intent(in)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      real(r8), intent(in)  ::  ZFLUX,RFL
      real(r8), intent(in),dimension(N_) :: FZ,ZTAU

      real(r8), intent(out),dimension(M_,M_,N_) ::  B,AA,CC
      real(r8), intent(out),dimension(M_,N_) ::  A,C,H

      integer I, J, K, L1,L2,LL
      real(r8)  SUM0, SUM1, SUM2, SUM3
      real(r8)  DELTAU, D1, D2, SURFAC
c
      real(r8), dimension(M_,M_) :: S,T,U,V,W
c---------------------------------------------

C---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
       do I = 1,M_
        SUM0 = 
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = 
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = 
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = 
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5_r8*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5_r8*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = 
     &   POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = 
     &   POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = 
     &   POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = 
     &   POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5_r8*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5_r8*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0_r8
         T(I,I)   = T(I,I)   + 1.0_r8
         V(I,I)   = V(I,I)   + 1.0_r8
         B(I,I,L1)= B(I,I,L1) + 1.0_r8

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) 
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) 
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25_r8*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0_r8*D2*C(I,L1)
         A(I,L1) = 0.0_r8
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

c------------intermediate points:  can be even or odd, A & C diagonal
c---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*(
     &     POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4)
     &   + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = 
     &     POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4)
     &    +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0_r8
        enddo
       enddo

c---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*(
     &     POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3)
     &   + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 = 
     &     POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3)
     &    +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0_r8
        enddo
       enddo

C---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 = 
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 = 
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 = 
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 = 
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5_r8*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5_r8*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 = 
     &   POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 = 
     &   POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 = 
     &   POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 = 
     &   POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5_r8*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5_r8*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0_r8
         T(I,I)   = T(I,I)   + 1.0_r8
         V(I,I)   = V(I,I)   + 1.0_r8
         B(I,I,L1)= B(I,I,L1) + 1.0_r8

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) 
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) 
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

C------------lower boundary, 2nd-order, A-matrix is full (AA)
         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25_r8*DELTAU
         SURFAC = 4.0_r8*RFL/(1.0_r8 + RFL)
       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
          SUM1 = SURFAC*SUM0
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) - 2.0_r8*D2*C(I,L1) + SUM0*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0_r8
       enddo

      return
      end



C-----------------------------------------------------------------------
      subroutine EFOLD (F0, F1, N, F)
C-----------------------------------------------------------------------
C     ***not used in fast-J, part of original scattering code
C---Speciality subroutine for calculating consistent exp(-tau/mu0)
C---  values on the tau grid so that photons are conserved.
C---  ***only works for plane-parallel, NOT psuedo-spherical atmos
C
C---  calculate the e-fold between two boundaries, given the value
C---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
C---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
C---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
C---           alternatively, could define A = ln[F0/F1]
C---  let X = A*x, d2F/dX2 = F
C---  assume equal spacing (not necessary, but makes this easier)
C---      with N-1 intermediate points (and N layers of thickness dX = A/N)
C---
C---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
C---      let D = 1 / dX*dX:
C
C  1  |   1        0        0        0        0        0   |    | F0 |
C     |                                                    |    | 0  |
C  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
C     |                                                    |    | 0  |
C  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
C     |                                                    |    | 0  |
C     |   0        0       -D      2D+1      -D        0   |    | 0  |
C     |                                                    |    | 0  |
C  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
C     |                                                    |    | 0  |
C N+1 |   0        0        0        0        0        1   |    | F1 |
C      
C-----------------------------------------------------------------------
C  Advantage of scheme over simple attenuation factor: conserves total
C  number of photons - very useful when using scheme for heating rates.
C  Disadvantage: although reproduces e-folds very well for small flux
C  differences, starts to drift off when many orders of magnitude are
C  involved.
C-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer, intent(in) :: N
      real(r8), intent(in)  :: F0,F1
      real(r8), intent(out) :: F(101)     !F(N+1)  
      integer I
      real(r8) A,DX,D,DSQ,DDP1, B(101),R(101)
C
      if (F0 .eq. 0._r8) then
        do I = 1,N
          F(I)=0._r8
        enddo
        return
      elseif (F1.eq.0._r8) then
        A = log(F0/1.e-250_r8)
      else
        A = log(F0/F1)
      endif
      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1._r8
      B(2) = DDP1
      R(2) = +D*F0
      do I = 3,N
        B(I) = DDP1 - DSQ/B(I-1)
        R(I) = +D*R(I-1)/B(I-1)
      enddo
      F(N+1) = F1
      do I = N,2,-1
        F(I) = (R(I) + D*F(I+1))/B(I)
      enddo
      F(1) = F0
      return
      end

c<<<<<<<<<<<<<<<<<<<<<<<<<end core scattering subroutines<<<<<<<<<<<<<<<

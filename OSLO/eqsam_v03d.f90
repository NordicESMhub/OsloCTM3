module eqsam_v03d

  implicit none
  
  public
  
contains
  
subroutine eqsam_v03d_sub(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
use cmn_precision, only: r8
use cmn_parameters, only: AVOGNR, R_ATM, R_UNIV, J2kcal
implicit none
!___________________________________________________________________________________________________________________________________
!      Written by Swen Metzger 3/11/99. Modified 2002, 2003.
!
!      Department of Atmospheric Chemistry, Max-Planck-Institute for Chemistry.
!      email: metzger@mpch-mainz.mpg.de
!      http://www.mpch-mainz.mpg.de/~metzger
!
!      COPYRIGHT 1999-2003
!
!      purpose
!      -------
!      EQSAM (EQuilibrium Simplified Aerosol Model) is a new and computationally efficient thermodynamic
!      aerosol composition model that allows to calculate the gas/aerosol equilibrium partitioning,
!      including aerosol water, sufficiently fast and accurate for global (or even regional) modeling.
!      EQSAM is based on a number of parameterizations, including single solute molalities and activity 
!      coefficients (AC). The thermodynamic framework (domains and subdomains, internally mixed aerosols) 
!      is the same as of more sophisticated thermodynamic equilibrium models (EQMs), e.g. of ISORROPIA 
!      (Nenes et al., 1998). Details are given in the references below (and the references therein).
!
!      The main assumption on which EQSAM/EQMs are based is thermodynamical and chemical equilibrium. 
!      From this assumption it directly follows that the aerosol water activity (aw) equals the ambient 
!      relative humidity (RH), if the water vapor pressure is sufficiently larger than the partial vapor
!      pressure of the aerosol compounds. This is approximately true for tropospheric aerosols. Given the 
!      large amount of water vapor present, water vapor and aerosol water equilibrate relatively faster 
!      compared to all other aerosol compounds. This is subsequently also true for single aerosol compounds.
!      The water activity of single solutes must also equal RH under this assumption. Therefore, the so 
!      called ZSR-relation is (and can be) used to calculate the aerosol associated water mass (simply
!      from the sum of all water mass fractions that are derived from measured single solute molalities). 
!
!      In contrast to other EQMs, EQSAM utilizes the fact that the RH fixes the water activity 
!      (under the above assumptions) and the consequence that any changes in RH also causes changes in 
!      the aerosol water mass and, hence, aerosol activity (including activity coefficients). Thus, an decrease
!      (increase) in RH decrease (increases) the aerosol water mass (and water activity). This can change the
!      aerosol composition, e.g. due to condensation (evaporation/crystallization), because the vapor pressure 
!      above the aerosol reduces (increases). In turn, a vapor pressure reduction (increase) due to changes
!      in the aerosol composition is compensated by an associated condensation (evaporation) of water vapor 
!      to maintain the aerosol molality to remain constant (because aw=RH). Furthermore, the aerosol water 
!      mainly depends on the aerosol mass and the type of solute, so that parameterizations of single solute 
!      molalities and activity coefficients can be defined, only depending on the type of solute and RH.
!      The advantage of using such parameterizations is that the entire aerosol equilibrium composition 
!      can be solved analytically, i.e. non-iteratively, which considerably reduces the amount of CPU time 
!      that is usually need for aerosol thermodynamic calculations (especially if an EQM is incorporated in
!      an aerosol dynamical model that is in turn embedded in a high resolution regional or global model).
!
!      However, EQSAM should still be regarded as a starting point for further developments. There is still 
!      room for improvements. For instance, this code is not yet numerically optimized (vectorized) and a 
!      number of improvements with respect to an explicit treatment of additional equilibrium reactions,
!      missing (or only implicit) dissociation, and a basic parameterization of the water uptake. 
!      
!      Note that EQSAM was originally developed to calculate the gas/aerosol equilibrium partitioning of the 
!      ammonium-sulfate-nitrate-water system for climate models, excluding solid compounds. 
!      This version (eqsam_v03d.f90) is extended with respect to sea salt. Solids/hysteresis are treated in a 
!      simplified manner. Results of a box model comparison with ISORROPIA will be available from the web page.
!      Please also note that the water uptake is based on additional (unpublished) parameterizations for single 
!      solute molalities, which are derived from tabulated measurements used in ISORROPIA. Note further that 
!      this extended version (eqsam_v03d.f90) is not yet published. A publication is in progress.
!
! ToDo:
!     Split ion-pairs into ions for water parameterizations (since info is actually available)
!     Include uptake/dissociation of NH3, HNO3, HCl (mainly to get pH right at near neutral conditions)
!     Extension to K+,Ca++,Mg++, CO2/(CO3)2--/HCO3-,SOA,etc.. (maybe not)
!     Vectorization. Translation of hardcoded formulas in array syntax.
!     I/O Interface and program structure clean up.
!     EQSAM info webpage.
!
! Version History:
!
!  eqsam_v03d.f90 (MPI-CH, June 2003): 
!   - gama parameterizations now according to Metzger 2002 (JGR Appendix)
!   - improved pH calculations (still restricted to strong acids)
!   - removed bug that lead to too high nitrate formation at dry and cold regions (UT/LS) 
!   - removed bug in solid/hysteresis calculations 
!     (both bugs introduced in eqsam_v03b.f90 by cleaning up eqsam_v02a.f90)
!   
!  eqsam_v03c.f90 (MPI-CH, April 2003):
!   - more accurate paramterizations of single solute molalities (Na, Cl species)
!   - cleanded up RHD subdomain structure
!   - improved water uptake (Na, Cl species)
!
!  eqsam_v03b.f90 (MPI-CH, March 2003): 
!                 System extended to HCl,Cl-/Na+. 
!                 Parameterization (fit) of additional HNO3 uptake removed.
!                 Instead, complete analytical solution of equilibrium reactions, based on the AC-RH relationship.
!  eqsam_v03.f90  (IMAU, October 1999): 
!                 Test version (included in TM3).
!  eqsam_v02a.f90 (IMAU, April 2000):
!                 Box model version.
!  eqsam_v02.f90  (IMAU, October 1999):
!                 TM3 version.
!                 Version including solids and additional HNO3 uptake on acidic aerosols (parameterized).
!  eqsam_v01b.f90 (MPI-CH, January 2003):
!                 Same as eqsam_v01a.f90 (additional lines though uncommented for test purposes only).
!  eqsam_v01a.f90 (IMAU, April 2000):
!                 Box model version.
!  eqsam_v01.f90  (IMAU, October 1999):
!                 TM3 version.
!                 First and most basic version (without solids) for better vectorization (for global modeling).
!                 System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, H2O 
!                 based on equilibrium / internal mixture assumption / aw=rh / ZSR-relation
!                 parameterization of activcity coefficients (AC), i.e. an AC-RH relationship
!
!      
!      interface
!      ---------
!      call  eqsam_v03d(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
!      yi = input  array (imax, nca)
!      yo = output array (imax, nco)
!      imax = max loop (e.g. time steps)
!      nca >= 11
!      nc0 >= 35
!      iopt = 1 metastable 
!      iopt = 2 solids 
!      iopt = 3 hysteresis (metastable/solids) for online calculations
!      iopt = 31 hysteresis lower branch 
!      iopt = 32 hysteresis upper branch 
!      ipunit = I/O unit (can be skipped)
!      in = array        (can be skipped)
!         
!      method
!      ------
!      equilibrium / internal mixture assumption / aw=rh
!      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O 
!              (K+,Ca++,Mg++)
!      external
!      --------
!      program    eqmd.f90    (driver only needed for the box model version)
!      subroutine gribio.f90  (provides diagnostics output in grib/binary/ascii format)
!      
!      references
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000.
!         http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis, 
!         GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001102, 2002
!         http://www.agu.org/journals/jd/jd0216/2001JD001102/index.html
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld, 
!         GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001103, 2002.
!         http://www.agu.org/journals/jd/jd0216/2001JD001103/index.html
!___________________________________________________________________________________________________________________________________
real(r8),parameter                        :: RH_HIST_DW=1.50_r8          ! mean value for mixture of wet (2) and dry (1) gridboxes (needed for HYSTERESIS)
real(r8),parameter                        :: T0=298.15_r8, T1=298.0_r8, &
                                             AVO=AVOGNR, &
                                             R=R_ATM,  &                ! in cu.m*atm/deg/mole (82.0567e-6_r8)
                                             r_kcal  = R_UNIV/J2kcal ! Ideal gas constant [1.986e-3_r8 kcal K-1.mole-1]
real(r8),parameter                        :: RHMAX=0.99_r8, RHMIN=0.0001_r8 ! restrict to max / min RH
real(r8),parameter                        :: MWNH4=18._r8, MWSO4=96._r8, MWNO3=62._r8, MWCl=35.5_r8 ! mole mass of species considered
real(r8),parameter                        :: MWNa=23._r8,MWCa=40.1_r8,MWN=14._r8, MWS=32.1_r8
real(r8),parameter                        :: MWH20=55.51_r8*18.01_r8, ZERO=0._r8
real(r8),parameter                        :: GF1=0.25_r8, GF2=0.50_r8, GF3=0.40_r8, GF4=1.00_r8, K=2._r8           ! exponents of AC-RH functions
!______________________________________________
integer,parameter                         :: NPAIR=10 
!
integer                                   :: ii,il,IHYST
integer,intent(in)                        :: nca,nco,imax,loop,ipunit
integer,intent(inout)                     :: iopt
!______________________________________________
integer,dimension(6),intent(in)           :: in
!______________________________________________
real(r8)                                  :: T0T,TT,RH,PX,RHD,KAN,KAC,ZIONIC,RH_HIST,GAMA,GG,GF,GFN
real(r8)                                  :: X00,X01,X02,X03,X04,X05,X08,X09,X10,X11
real(r8)                                  :: X0,X1,X2,X3,X4,X5,X6,XK10,XK6
real(r8)                                  :: ZFLAG,ZKAN,ZKAC,PH,COEF,HPLUS,AKW,XKW,MOLAL
real(r8)                                  :: TNH4,TSO4,TNO3,TNa,TCl,TPo,TCa,TMg
real(r8)                                  :: PNH4,PSO4,PNO3,PCl,PNa,GNO3,GNH3,GSO4,GHCl
real(r8)                                  :: ASO4,ANO3,ANH4,ACl,ANa,SNH4,SSO4,SNO3,SCl,SNa
real(r8)                                  :: WH2O,PM,PMs,PMt,RINC,DON,RATIONS,GR,NO3P,NH4P
!_______________________________________________
real(r8),dimension(imax,nca),intent(in)   :: yi
real(r8),dimension(imax,nco),intent(out)  :: yo
real(r8),dimension(8)                     :: w1,w2
real(r8),dimension(8)                     :: RHDA,RHDE,RHDX,RHDZ    ! RHD / MRHD arrays for different aerosol types
real(r8),dimension(NPAIR)                 :: M0,MW,NW,ZW            ! arrays of ion pairs
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
! mole mass of the salt solute
DATA MW(1:NPAIR)/ 58.5_r8, 142._r8,  88._r8, 132._r8, 80._r8, 53.5_r8, 98._r8, 115._r8, 120._r8, 247._r8/
! square of max. dissocation number (not consistent)
DATA NW(1:NPAIR)/  2._r8,    2.5_r8,  2.5_r8,  2.5_r8, 3.5_r8, 1._r8,   4.5_r8,  2._r8,   2._r8,   2.5_r8/
! exponents of water activity functions/
DATA ZW(1:NPAIR)/  0.67_r8,  1._r8,   1._r8,   1._r8,  1._r8,  1._r8,   0.5_r8,  1._r8,   1._r8,   1._r8/
!
! RHD / MRHD values as of ISORROPIA / SCAPE (T=298.15K)
DATA RHDA(1:8)/0.32840_r8, 0.4906_r8, 0.6183_r8, 0.7997_r8, 0.67500_r8, 0.5000_r8, 0.4000_r8, 0.0000_r8/
! Temp. coeff.
DATA RHDE(1:8)/-1860.0_r8, -431.0_r8, 852.00_r8, 80.000_r8, 262.000_r8, 3951.0_r8, 384.00_r8, 0.0000_r8/
!___________________________________________________________________________________________________________________________________
IHYST=2
IF(IOPT.EQ.31) THEN      ! SOLID HYSTORY
   IHYST=1
   IOPT=3
ELSE IF(IOPT.EQ.32) THEN  ! WET   HISTORY
   IHYST=2
   IOPT=3
END IF

!write(ipunit,*)'eqsam_v03d ...'
!print*,'                                                          '
!print*,'              EQuilibrium Simplified Aerosol Model (EQSAM)'
!print*,'                        for global modeling               ' 
!print*,'                                 by                       '
!print*,'                         Swen Metzger, MPI-CH             '
!print*,'                         Copyright 1999-2003              '
!print*,'                    >> metzger@mpch-mainz.mpg.de <<       '
!print*,'                     last change:  04. June, 2003         '
!print*,'                           (version 3.0d)                 '
!print*,'                   gas/aerosol calculations assuming      '
!print*,'                  System: NH3,NH4+/H2SO4+,HSO4-,SO4--     '
!print*,'                      HNO3,NO3-, HCl,Cl-/Na+, H2O         '
!if(iopt.eq.1) then
!print*,'                         metastable aeorsols              '
!elseif(iopt.eq.2) then
!print*,'                            solid aeorsols                '
!elseif(iopt.eq.3) then
!print*,'                             hysteresis                   '
!print*,'                         (metastable/solids)              '
!if(IHYST.eq.1) then
!print*,'                            solid hystory                 '
!elseif(IHYST.eq.2) then
!print*,'                             wet hystory                  '
!endif
!endif
!print*,'                                                          '
!print*,'loop over ',loop,' data sets'
!print*,'   '
!___________________________________________________________________________________________________________________________________
! init/reset
yo=0._r8
w1=0._r8
w2=0._r8
!___________________________________________________________________________________________________________________________________
do il=1,loop

! get old relative humidity to calculate aerosol hysteresis (online only)

   RH_HIST = 2._r8                                                                                   ! WET HISTORY (DEFAULT)
   IF(IHYST.EQ.1 .OR. IOPT.EQ.2)  RH_HIST = 1._r8                                                      ! SET TO SOLIDS

!  meteorology
   TT = yi(il,1)                    ! T                      [K]
   RH = yi(il,2)                    ! RH                     [0-1]
   PX = yi(il,11)                   ! p                      [hPa]
!
! gas+aerosol:
   w1(1) = yi(il,6)                 ! Na+ (ss  + xsod) (a)   [umol/m^3 air]
   w1(2) = yi(il,4)                 ! H2SO4    + SO4-- (p)   [umol/m^3 air]
   w1(3) = yi(il,3)                 ! NH3  (g) + NH4+  (p)   [umol/m^3 air]
   w1(4) = yi(il,5)                 ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
   w1(5) = yi(il,7)                 ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
   w1(6) = yi(il, 8)                ! K+   (p) from Dust     [umol/m^3 air]
   w1(7) = yi(il, 9)                ! Ca++ (p) from Dust     [umol/m^3 air]
   w1(8) = yi(il,10)                ! Mg++ (p) from Dust     [umol/m^3 air]
!______________________________________________

   zflag=1._r8

   w1=w1*1.0e-6_r8                     ! [mol/m^3 air]

   TNa   = w1(1)                    ! total input sodium   (g+p) 
   TSO4  = w1(2)                    ! total input sulfate  (g+p) 
   TNH4  = w1(3)                    ! total input ammonium (g+p)
   TNO3  = w1(4)                    ! total input nitrate  (g+p) 
   TCl   = w1(5)                    ! total input chloride (g+p) 
   TPo   = w1(6)                    ! total input potasium (g+p) 
   TCa   = w1(7)                    ! total input calcium  (g+p)
   TMg   = w1(8)                    ! total input magnesium(g+p)

! SULFATE RICH

      if((w1(1)+w1(3)+w1(6)+2._r8*(w1(7)+w1(8))).le.(2._r8*w1(2))) then
          zflag=3._r8
      end if

! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1

      if((w1(1)+w1(3)+w1(6)+2._r8*(w1(7)+w1(8))).le.w1(2)) then
          zflag=4._r8
      end if

! SULFATE NEUTRAL CASE

      if((w1(1)+w1(3)+w1(6)+2._r8*(w1(7)+w1(8))).gt.(2._r8*w1(2))) then
          zflag=2._r8
      end if

! SULFATE POOR AND CATION POOR CASE

      if((w1(1)+w1(6)+2._r8*(w1(7)+w1(8))).gt.(2._r8*w1(2))) then       
          zflag=1._r8
      end if

      IF ( RH .LT. RHMIN ) RH=RHMIN
      IF ( RH .GT. RHMAX ) RH=RHMAX

! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs

      RHDX(:)=RHDA(:)*exp(RHDE(:)*(1._r8/TT-1._r8/T0))
      RHDZ(:)=RHDX(:)
      
! ACCOUNT FOR VARIOUS AMMOMIUM/SODIUM SULFATE SALTS ACCORDING TO MEAN VALUE AS OF ISORROPIA

      GG=2._r8                         ! (Na)2SO4 / (NH4)2SO4 IS THE PREFFERED SPECIES FOR SULFATE DEFICIENT CASES
      IF(ZFLAG.EQ.3._r8) THEN
         IF(RH.LE.RHDZ(7)) THEN       ! ACCOUNT FOR MIXTURE OF (NH4)2SO4(s) & NH4HSO4(s) & (NH4)3H(SO4)2(s) 
            GG=1.677_r8                !                        (Na)2SO4 &  NaHSO4
!           GG=1.5
         ELSE IF(RH.GT.RHDZ(7).AND.RH.LE.RHDZ(5)) THEN ! MAINLY (Na)2SO4 / (NH4)2SO4(s) & (NH4)3H(SO4)2(s)
            GG=1.75_r8
!           GG=1.5
         ELSE IF(RH.GE.RHDZ(5)) THEN   ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- & HSO4-
            GG=1.5_r8                  !  (Na)2SO4 &  NaHSO4
         END IF
      END IF
      IF(ZFLAG.EQ.4._r8) GG=1._r8       ! IF SO4 NEUTRALIZED, THEN ONLY AS NaHSO4 / NH4HSO4(S) OR  HSO4- / H2SO4

      RHD=RH
      IF(IOPT.EQ.2.OR.RH_HIST.LT.RH_HIST_DW) THEN   ! GET RHD FOR SOLIDS / HYSTERESIS
!
! GET LOWEST DELIQUESCENCE RELATIVE HUMIDITIES ACCORDING TO THE CONCENTRATION DOMAIN (APROXIMATION) 
! BASED ON RHD / MRHD ISORROPIA/SCAPE
!
      w2(:)=1._r8
      do ii=1,8
         if(w1(ii).le.1.e-12_r8) w2(ii)=0._r8  ! skip compound in RHD calculation if value is concentration is zero or rather small
      end do

! GET LOWEST RHD ACCORDING TO THE CONCENTRATION DOMAIN

! zflag=1. (cation rich)  ...
! 1. sea salt      aerosol          : RHDX(1)=MgCl2
! 2. mineral dust  aerosol          : RHDX(2)=Ca(NO3)2
!
! zflag=2. (sulfate neutral) ...
! 3. ammonium + nitrate             : RHDX(3)= NH4NO3
! 4. ammonium + sulfate             : RHDX(4)=(NH4)2SO4        
! 5. ammonium + sulfate mixed salt  : RHDX(5)=(NH4)3H(SO4)2, (NH4)2SO4        
! 6. ammonium + nitrate  + sulfate  : RHDX(6)=(NH4)2SO4, NH4NO3, NA2SO4, NH4CL
!
! zflag=3. (sulfate poor) ...
! 7. ammonium + sulfate  (1:1,1.5)  : RHDX(7)= NH4HSO4
!
! zflag=4. (sulfate very poor) ...
! 8. sulfuric acid                  : RHDX(8)= H2SO4       

   IF(ZFLAG.EQ.1._r8)THEN

      RHD=W2(1)+W2(5)                     ! Na+  dependency
      IF(RHD.EQ.0._r8)  RHDX(1)=1._r8 
      RHD=W2(6)+W2(7)+W2(8)               ! K+/Ca++/Mg++ dependency (incl. ss)
      IF(RHD.EQ.0._r8)  RHDX(2)=1._r8

      RHD=MINVAL(RHDX(1:2))

   ELSE IF(ZFLAG.EQ.2._r8)THEN

      RHD=W2(3)*W2(4)                     ! NH4+ & NO3-  dependency
      IF(RHD.EQ.0._r8)  RHDX(3)=1._r8
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(GG.NE.2._r8)   RHD=0.               ! account only for pure (NH4)2SO4
      IF(RHD.EQ.0._r8)  RHDX(4)=1._r8
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0._r8)  RHDX(5)=1._r8
      RHD=W2(2)+W2(3)+W2(4)+W2(5)         ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL dependency
      IF(RHD.EQ.0._r8)  RHDX(6)=1._r8

!     RHD=MINVAL(RHDX(3:4))
      RHD=MINVAL(RHDX(3:6))

   ELSE IF(ZFLAG.EQ.3._r8)THEN

      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0._r8)  RHDX(7)=1._r8
      RHD=RHDX(7)

   ELSE IF(ZFLAG.EQ.4._r8)THEN

      RHD=W2(2)                           ! H2SO4 dependency (assume no dry aerosol)
      IF(RHD.EQ.0._r8)  RHDX(8)=1._r8

      RHD=RHDX(8)

   END IF ! ZFLAG
   END IF ! SOLIDS

! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

      M0(:) = ((NW(:)*MWH20/MW(:)*(1._r8/RH-1._r8)))**ZW(:)

! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

      T0T=T0/TT
      COEF=1._r8+LOG(T0T)-T0T

! EQUILIBRIUM CONSTANT NH4NO3(s) <==> NH3(g) + HNO3(g) [atm^2] (ISORROPIA)

!Changed 20111116, RBS:      XK10 = 5.746e-17
      XK10 = 3.000e-17_r8
      XK10= XK10 * EXP(-74.38_r8*(T0T - 1._r8) + 6.120_r8*COEF)
      KAN = XK10 / (R*TT) / (R*TT)

! EQUILIBRIUM CONSTANT  NH4CL(s) <==> NH3(g) + HCL(g) [atm^2] (ISORROPIA)

      XK6  = 1.086e-16_r8
      XK6 = XK6 * EXP(-71.00_r8*(T0T - 1._r8) + 2.400_r8*COEF)
      KAC = XK6 / (R*TT) / (R*TT)

!
! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER H2O <==> H(aq) + OH(aq) [mol^2/kg^2] (ISORROPIA)

      XKW  = 1.010e-14_r8
      XKW = XKW *EXP(-22.52_r8*(T0T - 1._r8) + 26.920_r8*COEF)

! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

      GAMA=0._r8
      IF(RH.GE.RHD) GAMA=(RH**ZFLAG/(1000._r8/ZFLAG*(1._r8-RH)+ZFLAG))
      GAMA = GAMA**GF1            ! ONLY GAMA TYPE OF NH4NO3, NaCl, etc. NEEDED SO FAR

      GAMA=0._r8
      GFN=K*K                      ! K=2, i.e. condensation of 2 water molecules per 1 mole ion pair
      GF=GFN*GF1                   ! = GFN[=Nw=4] * GF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
                                   ! ONLY GAMA TYPE OF NH4NO3, NH4Cl, etc. needed so far

      IF(RH.GE.RHD) GAMA=RH**GF/((GFN*MWH20*(1._r8/RH-1._r8)))**GF1

      GAMA = min(GAMA,1._r8)        ! FOCUS ON 0-1 SCALE
      GAMA = max(GAMA,0._r8)
      GAMA = (1._r8-GAMA)**K        ! transplate into aqueous phase equillibrium and account for 
                                   ! enhanced uptake of aerosol precursor gases with increasing RH
                                   ! (to match the results of ISORROPIA)


! CALCULATE RHD DEPENDENT EQ: IF RH <  RHD => NH4NO3(s) <==> NH3 (g) + HNO3(g) (ISORROPIA)
!                             IF RH >> RHD => HNO3  (g)   -> NO3 (aq)

      X00  = MAX(ZERO,MIN(TNa,GG*TSO4))       ! MAX SODIUM   SULFATE
      X0   = MAX(ZERO,MIN(TNH4,GG*TSO4-X00))  ! MAX AMMOMIUM SULFATE
      X01  = MAX(ZERO,MIN(TNa-X00, TNO3))     ! MAX SODIUM   NITRATE
      X1   = MAX(ZERO,MIN(TNH4-X0,TNO3-X01))  ! MAX AMMOMIUM NITRATE
!
      X02  = MAX(ZERO,MIN(TNa-X01-X00,TCl))   ! MAX SODIUM   CHLORIDE
      X03  = MAX(ZERO,MIN(TNH4-X0-X1,TCl-X02))! MAX AMMOMIUM CHLORIDE

      X2   = MAX(TNH4-X1-X0-X03,ZERO)         ! INTERIM RESIDUAL NH3
      X3   = MAX(TNO3-X1-X01,ZERO)            ! INTERIM RESIDUAL HNO3
      X04  = MAX(TSO4-(X0+X00)/GG,ZERO)       ! INTERIM RESIDUAL H2SO4
      X05  = MAX(TCl-X03-X02,ZERO)            ! INTERIM RESIDUAL HCl
!     X06  = MAX(TNa-X02-X01-X00,ZERO)        ! INTERIM RESIDUAL Na (should be zero for electro-neutrality in input data)
!
      ZKAN=2._r8
      IF(RH.GE.RHD) ZKAN=ZKAN*GAMA

      X4   = X2 + X3
      X5   = SQRT(X4*X4+KAN*ZKAN*ZKAN)
      X6   = 0.5_r8*(-X4+X5)
      X6   = MIN(X1,X6)
      
      GHCl = X05                              ! INTERIM RESIDUAl HCl
      GNH3 = X2 + X6                          ! INTERIM RESIDUAl NH3
      GNO3 = X3 + X6                          ! RESIDUAl HNO3
      GSO4 = X04                              ! RESIDUAl H2SO4
      PNa  = X02 + X01 + X00                  ! RESIDUAl Na (neutralized)

      ZKAC=2._r8
      IF(RH.GE.RHD) ZKAC=ZKAC*GAMA

      X08   = GNH3 + GHCl
      X09   = SQRT(X08*X08+KAC*ZKAC*ZKAC)
      X10   = 0.5_r8*(-X08+X09)
      X11   = MIN(X03,X10)

      GHCl = GHCl + X11                       ! RESIDUAL HCl
      GNH3 = GNH3 + X11                       ! RESIDUAL NH3

! GO SAVE ...

      IF(GHCl.LT.0._r8)   GHCl=0._r8
      IF(GSO4.LT.0._r8)   GSO4=0._r8
      IF(GNH3.LT.0._r8)   GNH3=0._r8
      IF(GNO3.LT.0._r8)   GNO3=0._r8
      IF(PNa.LT.0._r8)    PNa=0._r8
      IF(GSO4.GT.TSO4) GSO4=TSO4
      IF(GNH3.GT.TNH4) GNH3=TNH4
      IF(GNO3.GT.TNO3) GNO3=TNO3
      IF(GHCl.GT.TCl)  GHCl=TCl
      IF(PNa.GT.TNa)   PNa=TNa
!     IF(PNa.LT.TNa)   print*,il,' PNa.LT.TNa => no electro-neutrality in input data! ',PNa,TNa


! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4>1, TEN BRINK, ET AL., 1996, ATMOS ENV, 24, 4251-4261)

!     IF(TSO4.EQ.ZERO.AND.TNO3.GT.ZERO.OR.TNO3/TSO4.GE.1.) RHD=RH

!     IF(IOPT.EQ.2.AND.RH.LT.RHD.OR.IOPT.EQ.2.AND.RH_HIST.LT.RH_HIST_DW) THEN        ! SOLIDS / HYSTERESIS
      IF(RH_HIST.EQ.1._r8.AND.RH.LT.RHD) THEN        ! SOLIDS / HYSTERESIS

       ! EVERYTHING DRY, ONLY H2SO4 (GSO4) REMAINS IN THE AQUEOUSE PHASE

         ANH4 = 0._r8
         ASO4 = 0._r8
         ANO3 = 0._r8
         ACl  = 0._r8
         ANa  = 0._r8

      ELSE  !  SUPERSATURATED SOLUTIONS NO SOLID FORMATION

         ASO4 = TSO4 - GSO4
         ANH4 = TNH4 - GNH3
         ANO3 = TNO3 - GNO3
         ACl  = TCl  - GHCl
         ANa  = PNa

      END IF ! SOLIDS / HYSTERESIS

! CALCULATE AEROSOL WATER [kg/m^3(air)]
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
      IF(ZFLAG.EQ.1._r8) WH2O = ASO4/M0( 2) + ANO3/M0(3) + ACl/M0(6)  !sulfate poor
      IF(ZFLAG.EQ.2._r8) WH2O = ASO4/M0( 9) + ANO3/M0(5) + ACl/M0(6)  !sulfate neutral
      IF(ZFLAG.EQ.3._r8) WH2O = ASO4/M0( 8) + ANO3/M0(5) + ACl/M0(6)  !sulfate rich
      IF(ZFLAG.EQ.4._r8) WH2O = ASO4/M0( 8) + GSO4/M0(7)              !sulfate very rich

! CALCULATE AQUEOUS PHASE PROPERTIES

!     PH    = 9999.
      PH    = 7._r8
      MOLAL = 0._r8
      HPLUS = 0._r8
      ZIONIC= 0._r8

      IF(WH2O.GT.0._r8) THEN            

         ! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER

         AKW=XKW*RH*WH2O*WH2O                         ! H2O <==> H+ + OH- with kw [mol^2/kg^2]
         AKW=AKW**0.5_r8                               ! [OH-] = [H+] [mol]

         ! Calculate hydrogen molality [mol/kg], i.e. H+ of the ions:
         !                                   Na+, NH4+, NO3-, Cl-, SO4--, HH-SO4- [mol/kg(water)]
         !                                   with [OH-] = kw/[H+]

         HPLUS = (-ANa/WH2O-ANH4/WH2O+ANO3/WH2O+ACl/WH2O+GG*ASO4/WH2O+GG*GSO4/WH2O+ & 
            SQRT(( ANa/WH2O+ANH4/WH2O-ANO3/WH2O-ACl/WH2O-GG*ASO4/WH2O-GG*GSO4/WH2O)**2+XKW/AKW*WH2O))/2._r8

         ! Calculate pH

         !PH=-ALOG10(HPLUS)                             ! aerosol pH 
         !++alfgr changed to
         PH=-1._r8*log10(HPLUS)

         ! Calculate ionic strength [mol/kg]

         ZIONIC=0.5_r8*(ANa+ANH4+ANO3+ACl+ASO4*GG*GG+GSO4*GG*GG+XKW/AKW*WH2O*WH2O)
         ZIONIC=ZIONIC/WH2O                                           ! ionic strength [mol/kg]
!        ZIONIC=min(ZIONIC,200.0)                                     ! limit for output
!        ZIONIC=max(ZIONIC,0.0)

      END IF ! AQUEOUS PHASE 
!
!-------------------------------------------------------
! calculate diagnostic output consistent with other EQMs ...
!
      ASO4 = ASO4 + GSO4                                     ! assuming H2SO4 remains aqueous

      TNa   = TNa  * 1.e6_r8                                    ! total input sodium   (g+p)  [umol/m^3]
      TSO4  = TSO4 * 1.e6_r8                                    ! total input sulfate  (g+p)  [umol/m^3]
      TNH4  = TNH4 * 1.e6_r8                                    ! total input ammonium (g+p)  [umol/m^3]
      TNO3  = TNO3 * 1.e6_r8                                    ! total input nitrate  (g+p)  [umol/m^3]
      TCl   = TCl  * 1.e6_r8                                    ! total input chloride (g+p)  [umol/m^3]
      TPo   = TPo  * 1.e6_r8                                    ! total input potasium (g+p)  [umol/m^3]
      TCa   = TCa  * 1.e6_r8                                    ! total input calcium  (g+p)  [umol/m^3]
      TMg   = TMg  * 1.e6_r8                                    ! total input magnesium(g+p)  [umol/m^3]
!
! residual gas:
      GNH3 = GNH3 * 1.e6_r8                                     ! residual NH3
      GSO4 = GSO4 * 1.e6_r8                                     ! residual H2SO4
      GNO3 = GNO3 * 1.e6_r8                                     ! residual HNO3
      GHCl = GHCl * 1.e6_r8                                     ! residual HCl

! total particulate matter (neutralized)
      PNH4=TNH4-GNH3                                         ! particulate ammonium   [umol/m^3]
      PNO3=TNO3-GNO3                                         ! particulate nitrate    [umol/m^3]
      PCl =TCl -GHCl                                         ! particulate chloride   [umol/m^3]
      PNa =TNa                                               ! particulate sodium     [umol/m^3]
      PSO4=TSO4                                              ! particulate sulfate    [umol/m^3]

! liquid matter
      ASO4 = ASO4 * 1.e6_r8                                     ! aqueous phase sulfate  [umol/m^3] 
      ANH4 = ANH4 * 1.e6_r8                                     ! aqueous phase ammonium [umol/m^3]
      ANO3 = ANO3 * 1.e6_r8                                     ! aqueous phase nitrate  [umol/m^3]
      ACl  = ACl  * 1.e6_r8                                     ! aqueous phase chloride [umol/m^3]
      ANa  = ANa  * 1.e6_r8                                     ! aqueous phase sodium   [umol/m^3]

! solid matter
      SNH4=PNH4-ANH4                                         ! solid phase ammonium   [umol/m^3]
      SSO4=PSO4-ASO4                                         ! solid phase sulfate    [umol/m^3]
      SNO3=PNO3-ANO3                                         ! solid phase nitrate    [umol/m^3]
      SCl =PCl -ACl                                          ! solid phase chloride   [umol/m^3]
      SNa =PNa -ANa                                          ! solid phase sodium     [umol/m^3]

! GO SAVE ...

      IF(SNH4.LT.0._r8)   SNH4=0._r8
      IF(SSO4.LT.0._r8)   SSO4=0._r8
      IF(SNO3.LT.0._r8)   SNO3=0._r8
      IF(SCl.LT.0._r8)    SCl=0._r8
      IF(SNa.LT.0._r8)    SNa=0._r8

      PM=SNH4+SSO4+SNO3+SNH4+SCl+SNa+ANH4+ASO4+ANO3+ACl+ANa  ! total PM [umol/m^3]
      PMs=SNH4*MWNH4+SSO4*MWSO4+SNO3*MWNO3+SCl*MWCl+SNa*MWNa ! dry particulate matter (PM)     [ug/m^3]
      PMt=PMs+ANH4*MWNH4+ASO4*MWSO4+ANO3*MWNO3+ACl*MWCl+  &
          ANa*MWNa                                           ! total (dry + wet) PM, excl. H20 [ug/m^3]

      WH2O = WH2O * 1.e9_r8                                     ! convert aerosol water from [kg/m^3] to [ug/m^3]
      IF(WH2O.LT.1.e-3_r8) WH2O=0._r8

! UPDATE HISTORY RH FOR HYSTERESIS (ONLINE CALCULATIONS ONLY)

      RH_HIST=2._r8                                           ! wet
      IF(WH2O.EQ.0.) RH_HIST=1._r8                            ! dry

      RINC = 1._r8
      IF(PMt.GT.0._r8)   RINC = (WH2O/PMt + 1._r8)**(1._r8/3._r8)  ! approx. radius increase due to water uptake
      IF(RINC.EQ.0._r8)  RINC = 1._r8

      RATIONS = 0._r8
      IF(PSO4.GT.0._r8) RATIONS = PNO3/PSO4                   ! nitrate / sulfate mol ratio

      GR = 0._r8
      IF(GNO3.GT.0._r8) GR = GNH3/GNO3                        ! gas ratio = residual NH3 / residual HNO3   [-]

      DON = 0._r8
      IF((PNO3+2._r8*PSO4).GT.0._r8) DON = 100._r8*PNH4/(PNO3+2._r8*PSO4)! degree of neutralization by ammonia : ammonium / total nitrate + sulfate  [%]

      NO3P = 0._r8
      IF(TNO3.GT.0._r8) NO3P = 100._r8*PNO3/TNO3               ! nitrate  partitioning = nitrate  / total nitrate    [%]

      NH4P = 0._r8
      IF(TNH4.GT.0._r8) NH4P = 100._r8*PNH4/TNH4               ! ammonium partitioning = ammonium / total ammonium   [%]
!
! store aerosol species for diagnostic output:
!______________________________________________________________
! Input values:
      yo(il, 1) = TT   - 273.15_r8                                     ! T                                    [degC]
      yo(il, 2) = RH   * 100.00_r8                                     ! RH                                      [%]
      yo(il, 3) = TNH4                                                ! total input ammonium (g+p)       [umol/m^3]
      yo(il, 4) = TSO4                                                ! total input sulfate  (g+p)       [umol/m^3]
      yo(il, 5) = TNO3                                                ! total input nitrate  (g+p)       [umol/m^3]
      yo(il, 6) = TNa                                                 ! total input sodium   (p)         [umol/m^3]
      yo(il,33) = TCl                                                 ! total input chloride (g+p)       [umol/m^3]
      yo(il, 7) = TPo                                                 ! total input potasium (p)         [umol/m^3]
      yo(il,34) = TCa                                                 ! total input calcium  (p)         [umol/m^3]
      yo(il,35) = TMg                                                 ! total input magnesium(p)         [umol/m^3]
      yo(il,25) = PX                                                  ! atmospheric pressure                  [hPa]
! Output values:
      yo(il, 8) = GHCL                                                ! residual HCl   (g)               [umol/m^3]
      yo(il, 9) = GNO3                                                ! residual HNO3  (g)               [umol/m^3]
      yo(il,10) = GNH3                                                ! residual NH3   (g)               [umol/m^3]
      yo(il,11) = GSO4                                                ! residual H2SO4 (aq)              [umol/m^3]
      yo(il,12) = WH2O                                                ! aerosol Water  (aq)                [ug/m^3]
      yo(il,13) = PH                                                  ! aerosol pH                            [log]
      yo(il,14) = ZFLAG                                               ! concnetration domain [1=SP,2=SN,3=SR,4=SVR]
      yo(il,15) = PM                                                  ! total particulate matter         [umol/m^3]
      yo(il,16) = SNH4                                                ! solid ammonium (s)               [umol/m^3]
      yo(il,17) = SNO3                                                ! solid nitrate  (s)               [umol/m^3]
      yo(il,18) = SSO4                                                ! solid sulfate  (s)               [umol/m^3]
      yo(il,19) = PNH4                                                ! particulate ammonium (p=a+s)     [umol/m^3]
      yo(il,20) = PNO3                                                ! particulate nitrate  (p=a+s)     [umol/m^3]
      yo(il,21) = PSO4                                                ! particulate sulfate  (p=a+s)     [umol/m^3]
      yo(il,22) = RATIONS                                             ! mol ratio Nitrate/Sulfate (p)           [-]
      yo(il,23) = GAMA                                                ! activity coefficient (e.g. NH4NO3)           [-]
      yo(il,24) = ZIONIC                                              ! ionic strength (aq)                [mol/kg]
      yo(il,26) = PMt                                                 ! total PM (liquids & solids)        [ug/m^3]
      yo(il,27) = PMs                                                 ! total PM (solid)                   [ug/m^3]
      yo(il,28) = RINC                                                ! radius increase (H2O/PMt+1)**(1/3)      [-]
      yo(il,29) = SCl                                                 ! solid chloride (s)               [umol/m^3]
      yo(il,30) = SNa                                                 ! solid sodium (s)                 [umol/m^3]
      yo(il,31) = PCl                                                 ! particulate chloride (p=a+s)     [umol/m^3]
      yo(il,32) = PNa                                                 ! particulate sodium (p=a+s)       [umol/m^3]
end do
!

end subroutine eqsam_v03d_sub


end module eqsam_v03d

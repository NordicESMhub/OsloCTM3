c-----------------------------------------------------------------------
c---(p-linoz.f)    p-code 5.5a (Prather 5/2008 - fix colm & interp)
c---Linoz chemistry subs:  
c        LNZ_INIT, LNZ_SET, LNZ_PML, LNZ_SETO3, INT_LAT, INT_SOM, INT_MID
c        DECAY, TPAUSEB, TPAUSEG
c
c-----------------------------------------------------------------------
      subroutine LNZ_INIT
c-----------------------------------------------------------------------
c---read/init Linoz data & other strat chem tables from pratmo box model
c-----------------------------------------------------------------------

c---converts standard stratospheric tables (1:NLINZT) from PRATMO chemistry model
c---  STRTBL is converted in SCTM with vertical moments ___0, ___1, ___2
c---  O3 colm above top CTM layer calc as LZO3COL(J,M) (DU)

c---LINOZ tables are pre-computed as STRTBL(25,18,12,7):
c     NLINZT = 1:7 = no of std Linoz tables (can be expanded)
c     NLINZL = 1:25 = pressure levels from 58 km down to 10 km by 2 km
c     NLINZJ = 1:18 = std latitudes 85S to 85N, 10 deg intervals
c     M      = 1:12 = month
c---alternative table no's are chemistries, eg, loss freq for N2O, CFCs, ...
c     these can be inclded if NT > 7
c---note that the 1st and 2nd order moments are not used now in Linoz
c     but have been important in tracer studies (N2O, CFCs, ...)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: JPAR, LPAR
      use cmn_ctm, only: ZEDG, YDGRD
      use cmn_chem, only: NLINZJ, NLINZL, NLINZT,
     &     O3SFCLZ, ZN2OSFC, ZCH4SFC, O3TAULZ, LZLBO3, E90VMR_TP,
     &     O3iso1, O3iso2, CHLORLZ, T195LZ,
     &     LZTBL0, LZTBL1, LZTBL2, TPSCLZ, ZPSCLZ, TACTLZ,
     &     LZMIN, LZO3COL
      use cmn_parameters, only: CPI180
      implicit none
c-----------------------------------------------------------------------

      character(len=80) TITL1, INITFI

      integer J,K,L,M,N,NL

      real(r8)  P0L(LPAR+1),YSTRT(NLINZJ)

      real(r8) XPSD,XPSLM1,XPSL,P0,P1,P2,F0,F1,F2,PS(NLINZL+6),
     &     F(NLINZL+5)

      real(r8) STRTBL(NLINZL,NLINZJ,12,NLINZT),TPARM(NLINZL,JPAR,12)
      real(r8) TLCAR(NLINZL)
      
c-----------------------------------------------------------------------

c---Assumes that primary input routine has already read in:
c    LSTRAT = .true. = do Linoz chemistry
c    N_LZ = tracer number for Linoz O3 (nominally = 1, must be >0)
c    O3SFCLZ = O3 value to reset lower boundary (mole fraction = v/v)
c    O3TAULZ = e-fold time (s) to force O3=>O3SFCLZ  (if < 300 sec, inst. reset)
c    LZLBO3 = L value of the top of Boudnary Layer to reset Linoz O3 to O3SFCLZ
c    E90VMR_TP = O3 value that determeines chemical tropopause (mole fraction)

c---Linoz coefficients: 7 tables, 12 months, 18 latitudes, 25 layers
c---z* layers are arranged in 2 km intervals from 10 km to 58 km.
      INITFI='tables/linoz2004_2006jpl'
      open(8,file=INITFI,status='old',form ='formatted')
        read(8,1001)  TITL1
          write(6,1001) TITL1
        do N = 1,NLINZT
         read(8,1001) TITL1
         do M = 1,12
          do J = 1,18
            read(8,1003) (STRTBL(K,J,M,N),K=NLINZL,1,-1)
          enddo
         enddo
        enddo
      close(8)

c---data for PSC-parameterized O3 loss: 
c---Cl loading for Polar O3 loss, scaled to 1987 figures (WMO'98 Ch.11)
c---      PSC loss = (CHLORLZ)**2 / (2.75)**2
c---      e.g., 2004 = 3.427, 1987 = 2.750
      INITFI='tables/polar_o3loss.dat'
      open(8,file=INITFI,status='old',form ='formatted')
        read(8,1001)  TITL1
         write(6,1001) TITL1
        read(8,1002)  O3SFCLZ, ZN2OSFC, ZCH4SFC
         write(6,1004) ' Surface O3Strat concentraion =',O3SFCLZ
         write(6,1004) ' Surface ZN2OSFC concentraion =',ZN2OSFC
         write(6,1004) ' Surface ZCH4SFC concentraion =',ZCH4SFC
        read(8,1002)  O3TAULZ
         write(6,1004)' First-order lifetime to attain o3ssfc =',O3TAULZ
        read(8,'(I10)') LZLBO3
        read(8,1002)  E90VMR_TP, O3iso1, O3iso2
         write(6,1004)
     &    ' Recommended tropopause threshold mix. ratio =',E90VMR_TP
         write(*,'(a,2es10.3)')
     &        ' STEFLUX isopleth 1/2 vmr =',O3iso1,O3iso2
C---default CHLORLZ value, can/should be reset for appropriate year
        read(8,1002)  CHLORLZ
         write(6,1004) ' Cl loading =',CHLORLZ
        read(8,1002)  T195LZ
         write(6,1004) ' T PSC threshold =',T195LZ
        read(8,1001)  TITL1
         write(6,1001) TITL1
        read(8,1003) (TLCAR(K),K=NLINZL,1,-1)
      close(8)

 1001 format (a80)
 1002 format (8e10.3)
 1003 format (20x,6e10.3/(8e10.3))
 1004 format (a,1p,e11.3)

c---set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c---scan downward from 58 km to 10 km in 2 km intervals, constant >58 km
c---assume STRTBL defined at z* level and linearly interpolated to z* +-2km 
c---no. of 2-km layers from highest alt to surface (lowest 5 layers not defined)
c---
c---PS(1) = 58 km, ..., PS(25) = 10 km, PS(30=NL) = 0 km
      NL = NLINZL + 5
        XPSD       = 10._r8 ** (0.125_r8)
        PS(1)      = 1000._r8 * 10._r8**(-58._r8 / 16._r8)
      do L = 2,NL
        PS(L)    = PS(L-1) * XPSD
      enddo

c---pratmo table latitudes
        YSTRT(1) = -85._r8
      do J = 2,NLINZJ
        YSTRT(J) = YSTRT(J-1) + 10._r8
      enddo

c---lowest layer to allow Linoz chemistry: table values below 10 km are zero.
        LZMIN = 1
      do L = 1,LPAR
        if (ZEDG(L) .gt. 237._r8)  LZMIN = L
      enddo

c---CTM layer edges down from top
      do K = 1,LPAR+1
        P0L(K)   = ZEDG(LPAR+2-K)
      enddo

c---loop over parameter tables
      do N = 1,NLINZT

c---interpolate vs. latitude, from pratmo's STRTBL to CTM's TPARM(temporary)
        call INT_LAT(NLINZL,NLINZJ, JPAR,JPAR, YSTRT,YDGRD,
     &                         STRTBL(1,1,1,N),TPARM)

        do M = 1,12
        do J = 1,JPAR

c---load the profile value in F(:), 0-2-4-...-58 km (zero lowest L)
          do K = 1,NLINZL
            F(K) = TPARM(K,J,M)
          enddo
            K=NLINZL
            F(K+1) = F(K)
            F(K+2) = F(K)
            F(K+3) = F(K)
            F(K+4) = F(K)
            F(K+5) = F(K)

          if (N .ne. 3) then   
c---now integrate with moments the values over CTM layers L
           do L = 1,LPAR
             P1 = P0L(L+1)
             P2 = P0L(L)
             call INT_SOM(P1,P2,F0,F1,F2,PS,F,NL)
             LZTBL0(J,L,M,N) = F0
             LZTBL1(J,L,M,N) = F1
             LZTBL2(J,L,M,N) = F2
           enddo
          else
c---for O3 column, interpolate (not integrate) to the mid-point of each CTM L
           do L = 1,LPAR
             P0 = 0.5_r8*(P0L(L+1)+P0L(L))
             call INT_MID (P0,F0,PS,F,NL)
             LZTBL0(J,L,M,N) = F0
           enddo
             P0 = P0L(1)
             call INT_MID(P0,F0,PS,F,NL)
             LZO3COL(J,M) = F0
          endif

        enddo
        enddo

      enddo

c---set up the PSC loss parameters: (Cariolle et al. 1990): 
      do K = 1,NLINZL
        F(K) = TLCAR(K)
      enddo
        K=NLINZL
        F(K+1) = F(K)
        F(K+2) = F(K)
        F(K+3) = F(K)
        F(K+4) = F(K)
        F(K+5) = F(K)
      do L = 1,LPAR
         P1 = P0L(L+1)
         P2 = P0L(L)
        call INT_SOM(P1,P2,F0,F1,F2,PS,F,NL)
         TPSCLZ(L) = -F0
c---needs to be scaled by CHLORLZ squared: [Cl]**2 / [Cl @ 1987]**2
      enddo

c---maximum SZA for PSC loss (= tangent ht as sunset)
c---calc approx ht of mid layer =F0 (km) and estimate sunset SZA
      do L = 1,LPAR
        F0 = max(16._r8*log10(1000._r8/(0.5_r8*(ZEDG(L)+ZEDG(L+1)))),
     &           0._r8)
        F1 = (90._r8 + sqrt(F0)) * CPI180
        ZPSCLZ(L) = cos(F1)
      enddo

c---activation temperature: <195 K, and poleward of 40 deg only
      do J=1,JPAR
        if (abs(YDGRD(J)) .gt. 40._r8) then
          TACTLZ(J) = T195LZ
        else
          TACTLZ(J) = 0._r8
        endif
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine LNZ_SET (MONTH,NYR)
c-----------------------------------------------------------------------
c---set up Linoz (& other strat) table data for each day/month/year 
c-----------------------------------------------------------------------
c---LINOZ tables
c  1- ozone (Logan climatology), v/v
c  2- Temperature climatology, K                     
c  3- Column ozone climatology (Logan) above box mid-pt, DU
c  4- ozone (P-L) for climatological ozone, v/v/s
c  5- d(P-L) / dO3, 1/s
c  6- d(P-L) / dT, v/v/s/K
c  7- d(P-L) / d(column O3), v/v/s/DU
c  8- loss freq, 1/s for strat tracers (N2O, ...)
c-----------------------------------------------------------------------
      use cmn_size, only: JPAR, LPAR
      use cmn_chem, only: NLINZT, LZTBL0, LZTBL1, LZTBL2,
     &     LZPML0, LZPML1, LZPML2, O3TOPLZ, LZO3COL
      use utilities, only: ctmExitC
      implicit none
c-----------------------------------------------------------------------
      integer, intent(in):: MONTH,NYR

      integer J,L,N
c-----------------------------------------------------------------------

      if (NYR .lt. 1900) then
         call ctmExitC('>>>Linoz stop--cannot do YR<1900')
      endif


c---LINOZ parameter tables (1:NLINZT) are top-down, filled to CTM's L=1:LPAR
c---NB these could be interp to DAY/MONTH/YEAR
      do N = 1,NLINZT
       do L = 1,LPAR
        do J = 1,JPAR
          LZPML0(J,L,N) = LZTBL0(J,L,MONTH,N)
          LZPML1(J,L,N) = LZTBL1(J,L,MONTH,N)  ! moments not used by Linoz
          LZPML2(J,L,N) = LZTBL2(J,L,MONTH,N)  ! moments not used by Linoz
        enddo
       enddo
      enddo

c---climatology:  O3 column above model top for DAY/MONTH/YEAR
        do J = 1,JPAR
          O3TOPLZ(J) = LZO3COL(J,MONTH)
        enddo

c---polar O3 loss rates (Cariolle et al. 1990): just passed thru
c---   loss freq TPSCLZ(1:LM), SZA limit ZPSCLZ(1:LM), & activ T TLACT(1:JM)

c---need to reset CHLORLZ here if is to be changed over the simulation

      return
      end



c-----------------------------------------------------------------------
      subroutine LNZ_SETO3
c-----------------------------------------------------------------------
c--initialize Linoz O3 (N=N_LZ) based on supplied climatology
c-----------------------------------------------------------------------
      use cmn_precision, only: r8, rMom
      use cmn_size, only: IPAR, JPAR, LPAR, LLINOZ
      use cmn_ctm, only:
     &     AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW
      use cmn_chem, only: N_LZ, LZPML0, LZPML1, LZPML2, O3SFCLZ, LZMIN,
     &     LZONE, TMASSMIX2MOLMIX
      implicit none
c-----------------------------------------------------------------------
      integer I,J,L,LR

      if (LLINOZ .and. N_LZ.ne.0) then
c---zero moments - could use supplied vertical moments
        SUT(:,:,:,N_LZ)=0._r8
        SVT(:,:,:,N_LZ)=0._rMom
        SWT(:,:,:,N_LZ)=0._rMom
        SUU(:,:,:,N_LZ)=0._rMom
        SVV(:,:,:,N_LZ)=0._rMom
        SWW(:,:,:,N_LZ)=0._rMom
        SUV(:,:,:,N_LZ)=0._rMom
        SUW(:,:,:,N_LZ)=0._rMom
        SVW(:,:,:,N_LZ)=0._rMom
c---set Linoz O3 (N_LZ) to climatology above LZMIN, to O3SFCLZ below
c---assumes LNZ_SET has been called to set LZPML0 for the MONTH/DAY
        do L = LZMIN,LPAR
          LR = LPAR+1-L
          do J=1,JPAR
            do I=1,IPAR
             STT(I,J,L,N_LZ) = LZPML0(J,LR,1)*AIR(I,J,L)
     &                         / TMASSMIX2MOLMIX(N_LZ)
             SWT(I,J,L,N_LZ) = LZPML1(J,LR,1)*AIR(I,J,L)
     &                         / TMASSMIX2MOLMIX(N_LZ)
             SWW(I,J,L,N_LZ) = LZPML2(J,LR,1)*AIR(I,J,L)
     &                         / TMASSMIX2MOLMIX(N_LZ)
            enddo
          enddo
        enddo
        do L = 1,LZMIN-1
          do J=1,JPAR
            do I=1,IPAR
             STT(I,J,L,N_LZ) = O3SFCLZ*AIR(I,J,L)/TMASSMIX2MOLMIX(N_LZ)
            enddo
          enddo
        enddo
c---do not rescale the LINOZ tracer with LZONE 
        LZONE(N_LZ) = .false.
      endif

      print*, '* LINOZ tracer is initialized'
      return
      end



c-----------------------------------------------------------------------
      subroutine LNZ_PML(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ
     &                ,AIRB,BTEM,LSTRATAIR_E90B,UTTAU,DTCHEM,MP,
     &                LZLBO3X, N)
c-----------------------------------------------------------------------
c---Linoz = linearize P-L for stratospheric ozone based on tables from
c   the PRATMO model using climatological T, O3, Month
c-----------------------------------------------------------------------
      use cmn_precision, only: r8, rMom
      use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
      use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE,
     &     AREAXY, XDGRD, YDGRD, SOLDEC, SOLDIS
      use cmn_chem, only: N_LZ, O3SFCLZ, O3TAULZ, CHLORLZ, LZPML0,
     &     LZMIN, TPSCLZ, O3TOPLZ, TACTLZ, ZPSCLZ, TMASSMIX2MOLMIX
      use cmn_diag, only: O3PML, LFLXDG, O3PBLSINK
      use utilities, only: LOCSZA
      implicit none
c-----------------------------------------------------------------------
      real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: BTT
      real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK)::
     &         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
      real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: AIRB, BTEM
      real(r8), intent(in) :: UTTAU, DTCHEM
      integer, intent(in):: MP,N,LZLBO3X
      logical, intent(in):: LSTRATAIR_E90B(LPAR,IDBLK,JDBLK)

c---DOBF = 0.5*6.022d23*1000/(48.*1.d4*2.687d16) = weighting for half-layer
      real(r8),parameter ::  DOBF = 23345.45_r8
      real(r8),parameter ::  HRSEC2 = 0.5_r8/3600._r8

      real(r8)  AIRFAC,CLIMO3,CLIMPML,DERO3,DERTMP,DERCO3,DCO3,DTMP
      real(r8)  UTSZA,COSSZA,SOLFX, FXYZP(LPAR),COLO3(LPAR)
      real(r8)  STTNEW,STTOLD,PMLNET, FXYZL,FXYZS,CLPSC

      real(r8)  ZO3SSFC, SSO3
      integer I,J,II,JJ,L,LR

c---------------------set ups needed------------------------------------
c---Linoz applied to tracer = N:  nominally N=1 for strat or trop+Linoz runs
c---note that Lower Boundary is reset (L=1:LZLBO3) so LZLBO3=0 for trop-O3 runs
c
c---for ASAD, call LINOZ (... 1,4) for O3strat, and (... n,0) for true O3
c
c---LZMIN =lowest layer to allow Linoz or any strat-table chemistry: 
c---       i.e., no table values below 10 km in pressure alt (237.14 hPa)
c
c---O3TOPLZ = ozone column above top of model (climatology, see STRATL)
c
c---LZPML0(J,LR,1:7) tables loaded for this day/month/year (currently just month)
c       7 tables, each a function of 
c         month (12), 
c         latitude (18: -85 to 85 in 10 deg. increments)
c         altitude (25: z*=10-58 km in 2 km increments).
c  1- ozone (Logan climatology), v/v                 [CLIMO3]
c  2- Temperature climatology, K                     
c  3- Column ozone climatology (Logan) integrated above box, DU
c  4- ozone (P-L) for climatological ozone, v/v/s    [CLIMPML]
c  5- d(P-L) / dO3, 1/s                              [DERO3]
c  6- d(P-L) / dT, v/v/s/K                           [DERTMP]
c  7- d(P-L) / d(column O3), v/v/s/DU                [DERCO3]
c-----------------------------------------------------------------------

c---UT time (hrs) for PSC chemistry uses mid-point in time step 
c---(for J-values, UTSZA must be the end point in the time step!)
          UTSZA = UTTAU + HRSEC2*DTCHEM

c e-fold time for O3 at lower boundary (L=1:LZLBO3): O3TAULZ < 300s => just rest
          if (O3TAULZ .gt. 300._r8) then
            FXYZS = 1._r8 - exp(-DTCHEM/O3TAULZ)
          else
            FXYZS = 1._r8
          endif
          ZO3SSFC = O3SFCLZ/TMASSMIX2MOLMIX(N_LZ)


c---enhanced polar O3 loss (Cariolle 1990) with loss freq TPSCLZ(LR) < 0 
              FXYZP(:) = 0._r8
          if (CHLORLZ .gt. 0._r8) then
              CLPSC = (CHLORLZ/2.75_r8)**2
            do L = LZMIN,LPAR
              LR = LPAR+1-L
              FXYZP(LR) = 1._r8 - exp(TPSCLZ(LR)*CLPSC*DTCHEM)
            enddo
          endif

c-----------------------------------------------------------------------
c---major loop pair over IJ-block
c-----------------------------------------------------------------------
      do J = MPBLKJB(MP),MPBLKJE(MP)
         JJ     = J - MPBLKJB(MP) + 1

       do I = MPBLKIB(MP),MPBLKIE(MP)
          II     = I - MPBLKIB(MP) + 1

c---need SZA if doing enhanced, PSC Cl-driven O3 loss
        if (CHLORLZ .gt. 0._r8) then
         call LOCSZA(UTSZA,XDGRD(I),YDGRD(J),SOLDEC,SOLDIS,COSSZA,SOLFX)
        endif

c---ozone column above mid-layer for each CTM Level (in Dobson Units)
           COLO3(LPAR) = O3TOPLZ(J) + BTT(LPAR,N,II,JJ)*DOBF/AREAXY(I,J)
         do L = LPAR-1,1,-1
           COLO3(L) = COLO3(L+1) +
     &       (BTT(L,N,II,JJ)+BTT(L+1,N,II,JJ))*DOBF/AREAXY(I,J)
         enddo
c-----------------------------------------------------------------------
c---LINOZ chemistry over altitudes LZMIN:LM, boundary over 1:LZLBO3
        do L = LZMIN,LPAR
        if (LSTRATAIR_E90B(L,II,JJ)) then

          LR = LPAR+1-L
          STTOLD = BTT(L,N,II,JJ)
          STTNEW = STTOLD
          AIRFAC=AIRB(L,II,JJ)/TMASSMIX2MOLMIX(N)
c---climatological ozone (v/v = mole fraction mixing ratio)
          CLIMO3  = LZPML0(J,LR,1)*AIRFAC
c---climatological P-L            
          CLIMPML = LZPML0(J,LR,4)
c---partial derivative: d(P-L)/dO3 < 0 
          DERO3   = LZPML0(J,LR,5)
c---partial derivative: d(P-L)/dT
          DERTMP  = LZPML0(J,LR,6)
c---partial derivative: d(P-L)/dcol-O3
          DERCO3  = LZPML0(J,LR,7)
c---differences from climatology
          DTMP = BTEM(L,II,JJ) - LZPML0(J,LR,2)
          DCO3 = COLO3(L) - LZPML0(J,LR,3)

c---steady-state ozone:  can be negative in lower strat, but timescale is long
          SSO3 = CLIMO3 - AIRFAC*(CLIMPML+DCO3*DERCO3+DTMP*DERTMP)/DERO3
c---use DERO3 (<0) timescale to decay to Steady-State.
          PMLNET = (SSO3 - STTOLD)*(1._r8 - exp(DERO3*DTCHEM))
          STTNEW = STTOLD + PMLNET
          
c---PSC-activ loss depends on  T (incl. <195K & lat>40), sun above horizon(L)
          if (CHLORLZ .gt. 0._r8) then
           if (BTEM(L,II,JJ) .lt. TACTLZ(J)) then
            if (COSSZA .gt. ZPSCLZ(L)) then
              PMLNET = -STTOLD*FXYZP(LR)
              STTNEW = STTOLD + PMLNET
            endif
           endif
          endif

          if (LFLXDG .and. N.eq.N_LZ) then
            AIRFAC = AIRB(L,II,JJ)/TMASSMIX2MOLMIX(N)
c            if (BTT(L,N,II,JJ)/AIRFAC .le. O3iso1 )
c     &                    O3PML(I,J,1) = PMLNET + O3PML(I,J,1)
c            if (BTT(L,N,II,JJ)/AIRFAC .le. O3iso2 )
c     &                    O3PML(I,J,2) = PMLNET + O3PML(I,J,2)
c Rather check Linoz STE through e90 surface
            if (.not.LSTRATAIR_E90B(L,II,JJ))
     &           O3PML(I,J,4) = PMLNET + O3PML(I,J,4)
          endif

          BTT(L,N,II,JJ) = STTNEW
c---reduce moments if O3 is reduced
          if (STTNEW .lt. STTOLD) then
            FXYZL = STTNEW/STTOLD
            BXT(L,N,II,JJ) = BXT(L,N,II,JJ) * FXYZL
            BYT(L,N,II,JJ) = BYT(L,N,II,JJ) * FXYZL
            BZT(L,N,II,JJ) = BZT(L,N,II,JJ) * FXYZL
            BXX(L,N,II,JJ) = BXX(L,N,II,JJ) * FXYZL
            BYY(L,N,II,JJ) = BYY(L,N,II,JJ) * FXYZL
            BZZ(L,N,II,JJ) = BZZ(L,N,II,JJ) * FXYZL
            BXY(L,N,II,JJ) = BXY(L,N,II,JJ) * FXYZL
            BXZ(L,N,II,JJ) = BXZ(L,N,II,JJ) * FXYZL
            BYZ(L,N,II,JJ) = BYZ(L,N,II,JJ) * FXYZL
          endif 

        endif
        enddo

c---reset O3 at lower boundary (L=1:LZLBO3) using e-fold time O3TAULZ
        if (N .eq. N_LZ)  then
          do L = 1,LZLBO3X
            STTOLD = BTT(L,N,II,JJ)
            PMLNET = (ZO3SSFC*AIRB(L,II,JJ) - STTOLD) * FXYZS
            STTNEW = STTOLD + PMLNET
c            if (LFLXDG)  O3PBLSINK(I,J) = PMLNET + O3PBLSINK(I,J)
c Store in FLX diagnose 4
            if (LFLXDG)  O3PBLSINK(I,J,4) = PMLNET + O3PBLSINK(I,J,4)
            BTT(L,N,II,JJ) = STTNEW
c---reduce moments if O3 is reduced
            if (STTNEW .lt. STTOLD) then
              FXYZL = STTNEW/STTOLD
              BXT(L,N,II,JJ) = BXT(L,N,II,JJ) * FXYZL
              BYT(L,N,II,JJ) = BYT(L,N,II,JJ) * FXYZL
              BZT(L,N,II,JJ) = BZT(L,N,II,JJ) * FXYZL
              BXX(L,N,II,JJ) = BXX(L,N,II,JJ) * FXYZL
              BYY(L,N,II,JJ) = BYY(L,N,II,JJ) * FXYZL
              BZZ(L,N,II,JJ) = BZZ(L,N,II,JJ) * FXYZL
              BXY(L,N,II,JJ) = BXY(L,N,II,JJ) * FXYZL
              BXZ(L,N,II,JJ) = BXZ(L,N,II,JJ) * FXYZL
              BYZ(L,N,II,JJ) = BYZ(L,N,II,JJ) * FXYZL
            endif 
          enddo
        endif

       enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine INT_LAT(KPRA,JPRA,JCTM,JPAR, XPRA,XCTM, TPRA,TCTM)
c-----------------------------------------------------------------------
c---interpolate from pratmo standard tables (85S to 85N) -> CTM J-grid
c---   table TPRA with lat grid XPRA ==> table TCTM with lat grid XCTM
c---assume tables are for 12 months
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer, intent(in)::  KPRA,JPRA,JCTM,JPAR
      real(r8), intent(in)::   XPRA(JPRA),XCTM(JPAR),TPRA(KPRA,JPRA,12)
      real(r8), intent(out)::  TCTM(KPRA,JPAR,12)

      integer  I,II,J,K,M
      real(r8)   CNST1,CNST2
c----------------------------------------------------------------------
        J = 2 
c---this logic works only if CTM grid finer than 10 deg
       do I = 1,JCTM
        if (XCTM(I) .gt. XPRA(1))  then
          if (XCTM(I) .lt. XPRA(JPRA))  then
              CNST1 = (XPRA(J) - XCTM(I)) / (XPRA(J) - XPRA(J-1))
              CNST2 = 1._r8 - CNST1
            do M = 1,12
            do K = 1,KPRA
              TCTM(K,I,M) = CNST1*TPRA(K,J-1,M) + CNST2*TPRA(K,J,M)
            enddo
            enddo
             II = min(I+1,JCTM)
             if (XCTM(II) .gt. XPRA(J))  J = min(JPRA,J+1)
          else
c---CTM latitude > +85
            do M = 1,12
            do K = 1,KPRA
              TCTM(K,I,M) = TPRA(K,JPRA,M)
            enddo
            enddo
          endif
        else
c---CTM latitude < -85
            do M = 1,12
            do K = 1,KPRA
              TCTM(K,I,M) = TPRA(K,1,M)
            enddo
            enddo
        endif
       enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine INT_MID (P0,F0,PS,F,NL)
c-----------------------------------------------------------------------
c---For a CTM model level at mid-pt pressure = P0, 
c---    interpolates the value F0 from F linearly in p on the std grid PS
c---    assumes  PS(1) > PS(2) > PS(3) ... > PS(NL+1)
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer, intent(in)::   NL
      real(r8), intent(in) ::   P0,PS(NL+1),F(NL)
      real(r8), intent(out)::   F0

      integer  I
      real(r8)   XB
c-----------------------------------------------------------------------
        F0 = 0._r8

      do I = NL,1,-1
        if (PS(I) .lt. P0) then
c-------- have found:  PS(I) >= P0 > PS(I+1)--------------
          XB = (P0-PS(I))/(PS(I+1)-PS(I))
          XB = min( 1._r8, max( 0._r8, XB))
          F0 = F(I) + XB*(F(I+1) - F(I))
          goto 2
        endif
      enddo
          F0 = F(1)
    2 continue

      return
      end


c-----------------------------------------------------------------------
      subroutine INT_SOM (P1,P2,F0,F1,F2,PS,F,NL)
c-----------------------------------------------------------------------
c---For a CTM model level bounded by pressure P1 > P2, 
c---    integrates (p-avg) the value F0 from F on the std (2-km) grid PS
c---    assumes  top=PS(1) < PS(2) < PS(3) ... > PS(30) = 1000 mb
c---NOTE reverse order in P's
c---Assume that the quantity is constant over range halfway to layer above/below
c---    and calculate box edges from P=0 to P=1000

c---For a CTM model level between pressure range P1 > P2 (decreasing up)
c---calculate the SOM Z-moments of the loss freq at std z* (log-p) intervals
c--------  the pressure levels BETWEEN z* values are:
c                         PS(i) > PS(i+1) bounds z*(i)
c-------- The MOMENTS for a square-wave or 'bar': F(x)=F0  b<=x<=c, =0.0 else
c-----     S0 =   f0 (x)                      [from x=b to x=c]
c-----     S1 = 3 f0 (x^2 - x)                [from x=b to x=c]
c-----     S2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      implicit none
      integer, intent(in)::   NL
      real(r8), intent(in) ::   P1,P2,PS(NL),F(NL)
      real(r8), intent(out)::   F0,F1,F2

      integer  I
      real(r8)   XB,XC,PC,PB,SGNF0, PF1,PF2
c-----------------------------------------------------------------------
        F0 = 0._r8
        F1 = 0._r8
        F2 = 0._r8
      do I = 1,NL
       if (I.eq.1) then
         PF1 = 0._r8
       else
         PF1 = 0.5_r8*(PS(I-1)+PS(I))
       endif
       if (I.eq.NL) then
         PF2 = 1000._r8
       else
         PF2 = 0.5_r8*(PS(I)+PS(I+1))
       endif
        PC   = min(P1,PF2)
        PB   = max(P2,PF1)
        if (PC .gt. PB)  then
C--- have condition:  P1 .ge. PC .gt. PB .ge. P2 
C---      and          0 .le. XB .lt. XC .le. 1
          XC = (PC-P2)/(P1-P2)
          XB = (PB-P2)/(P1-P2)

c-------- assume that the quantity, F, is constant over interval [XLO,XUP],
c--------   F0: (c-b),   F1: 6((c2-c)-(b2-b)),  F2: 5((2c3-3c2+c)-(2b3-3b2+b))
c-------- calculate its contribution to the moments in the interval [0,1]
          F0 = F0 +F(I) *(XC -XB)
          F1 = F1 +F(I) *3._r8 *((XC *XC -XC) - (XB *XB -XB))
          F2 = F2 +F(I) *5._r8 *
     &         ((XC+XC-1._r8)*(XC*XC -XC) - (XB+XB-1._r8)*(XB*XB -XB))
        endif
      enddo

c---limiter on Z-moments: force monotonicity (tables can be + or -)
        SGNF0 = sign(1._r8, F0) 
        F0 = abs(F0)
      if (F2 .gt. 0._r8)  then
c-------- do.not.allow reversal of curvature: F2 > 0 -------------------
        F2   = min(F2, abs(F1)*0.333333_r8, 0.5_r8*F0)
        if (F1 .lt. 0._r8)  then
          F1 = max(-(F0+F2), F1)
         else
          F1 = min(+(F0+F2), F1)
        endif
      else
c-------- F2 < 0 = curved down at ends, allow if F1 < F0 ---------------
        F1  = min(F0, max(-F0, F1))
        F2  = max(F2, (abs(F1)-F0), (-abs(F1)*0.333333_r8))
      endif

        F0 = SGNF0 * F0
        F1 = SGNF0 * F1
        F2 = SGNF0 * F2

      return
      end


c-----------------------------------------------------------------------
      subroutine DECAY (BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ
     &                 ,DTCHEM,MP)
c-----------------------------------------------------------------------
c---provides a simple e-fold decay of species throughout the model domain
c---supplements ASAD chemistry for Rn-222, other labeled tracers, ...
c---only acts on tracers with an 'e' as first char in their name. 
c------------------------------------------------------------------------
      use cmn_precision, only: r8, rMom
      use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
      use cmn_ctm, only: NTM, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
      use cmn_chem, only: TNAME
      implicit none
c-----------------------------------------------------------------------
      real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: BTT
      real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK)::
     &         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
      real(r8), intent(in)::  DTCHEM
      integer, intent(in):: MP

      integer I,J,L,N
      real(r8)  F1L,DECTIM
c------------------------------------------------------------------------

      !// CTM3: Not in use (see oc_utilities.f: decay_e90)
      stop 'p-linoz: DECAY not in use for CTM3!'

      do N = 1,NTM
        if (TNAME(N)(1:1).eq.'e') then
          read(TNAME(N),'(1x,f3.0)') DECTIM
          if (DECTIM .gt. 0._r8) then
            F1L = exp(-DTCHEM/(86400._r8*DECTIM))
            do J = 1,MPBLKJE(MP)-MPBLKJB(MP)+1
              do I = 1,MPBLKIE(MP)-MPBLKIB(MP)+1
                do L = 1,LPAR
                 BTT(L,N,I,J) = BTT(L,N,I,J)*F1L
                 BZT(L,N,I,J) = BZT(L,N,I,J)*F1L
                 BZZ(L,N,I,J) = BZZ(L,N,I,J)*F1L
                 BXZ(L,N,I,J) = BXZ(L,N,I,J)*F1L
                 BYZ(L,N,I,J) = BYZ(L,N,I,J)*F1L
                 BXT(L,N,I,J) = BXT(L,N,I,J)*F1L
                 BXX(L,N,I,J) = BXX(L,N,I,J)*F1L
                 BYT(L,N,I,J) = BYT(L,N,I,J)*F1L
                 BYY(L,N,I,J) = BYY(L,N,I,J)*F1L
                 BXY(L,N,I,J) = BXY(L,N,I,J)*F1L
                enddo
              enddo
            enddo
          endif
        endif
        if (TNAME(N).eq.'O3f') then
          DECTIM = 4._r8  ! Synoz: 4 day efold to zero at L=1
          F1L = exp(-DTCHEM/(86400._r8 * DECTIM))
          do J = 1,MPBLKJE(MP)-MPBLKJB(MP)+1
            do I = 1,MPBLKIE(MP)-MPBLKIB(MP)+1
              BTT(1,N,I,J) = BTT(1,N,I,J)*F1L
              BZT(1,N,I,J) = BZT(1,N,I,J)*F1L
              BZZ(1,N,I,J) = BZZ(1,N,I,J)*F1L
              BXZ(1,N,I,J) = BXZ(1,N,I,J)*F1L
              BYZ(1,N,I,J) = BYZ(1,N,I,J)*F1L
              BXT(1,N,I,J) = BXT(1,N,I,J)*F1L
              BXX(1,N,I,J) = BXX(1,N,I,J)*F1L
              BYT(1,N,I,J) = BYT(1,N,I,J)*F1L
              BYY(1,N,I,J) = BYY(1,N,I,J)*F1L
              BXY(1,N,I,J) = BXY(1,N,I,J)*F1L
            enddo
          enddo
        endif

      enddo

      return
      end


      
c-----------------------------------------------------------------------
      subroutine TPAUSEB (BTT,AIRB,LSTRATAIR_E90B,MP)
c-----------------------------------------------------------------------
c define the tropopause level based on Linoz O3 (N=N_LZ)
c note that LSTRATAIR_E90(L,I,J) is a 3D logical, .true. = stratospheric
c called within the IJ-blocks, ONLY re-set of LSTRATAIR_E90 during the run
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
      use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
      use cmn_chem, only: Ne90, TNAME, E90VMR_TP, TMASSMIX2MOLMIX
      implicit none
c-----------------------------------------------------------------------
      real(r8),  Intent(in) ::  BTT(LPAR,NPAR,IDBLK,JDBLK)
      real(r8),  Intent(in) ::  AIRB(LPAR,IDBLK,JDBLK)
      Integer, Intent(in) ::  MP
      Logical, Intent(out) :: LSTRATAIR_E90B(LPAR,IDBLK,JDBLK)

      real(r8)   ZPAUSE
      integer  I, J, L

      !// CTM3: Not in use (see oc_utilities.f: tpauseb_e90)
      stop 'in TPAUSEB'

c Select O3 isopleth as chemical tropopause 
      ZPAUSE = E90VMR_TP/TMASSMIX2MOLMIX(Ne90)

      do J = 1,MPBLKJE(MP)-MPBLKJB(MP)+1
      do I = 1,MPBLKIE(MP)-MPBLKIB(MP)+1
        do L = 1,LPAR
          LSTRATAIR_E90B(L,I,J) =
     &          BTT(L,Ne90,I,J) .lt. ZPAUSE*AIRB(L,I,J)
        enddo
      enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine TPAUSEG
c-----------------------------------------------------------------------
c define the tropopause level based on Linoz O3 (N=N_LZ)
c note that LSTRATAIR_E90(I,J,L) is a 3D logical, .true. = stratospheric
c ***called/used during setup, based on init Linoz O3 or restart values
c-----------------------------------------------------------------------
      use cmn_precision, only: r8
      use cmn_size, only: IPAR, JPAR, LPAR
      use cmn_ctm, only: AIR, STT
      use cmn_chem, only: Ne90, E90VMR_TP, LSTRATAIR_E90,TMASSMIX2MOLMIX
      implicit none
c-----------------------------------------------------------------------
      real(r8)   ZPAUSE
      integer  I, J, L

      !// CTM3: Not in use (see oc_utilities.f: tpause_e90)
      stop 'in TPAUSEG'

      LSTRATAIR_E90(:,:,:) = .false.

      if (Ne90 .gt. 0) then
c Select O3 isopleth as chemical tropopause 
      ZPAUSE = E90VMR_TP/TMASSMIX2MOLMIX(Ne90)

       do J=1,JPAR
       do I=1,IPAR
         do L=1,LPAR
           LSTRATAIR_E90(L,I,J) = STT(I,J,L,Ne90) .lt. ZPAUSE*AIR(I,J,L)
         enddo
       enddo
       enddo

      endif
      return
      end

c---(p-cloud2.f)----new (10/2008) revised off-line cloud fraction generator
c---(p-code 6.0a)
!// Amund Sovde, November 2010
!// I have tried to include more comments, they start with '!//'.
!
!-----------------------------------------------------------------------
      subroutine CLOUD(CLDFRW,CLDIWCW,CLDLWCW,PW,TW,ETAA,ETAB,AREAXYW
     &                ,ZOFLEW,ZDEGI,ZDEGJ,IMAP,JMAP
     &                ,SWSTORE,CLDSTORE,TYPSTORE,RAN4,IRAN0
     &                ,LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)
!-----------------------------------------------------------------------

      use cmn_precision, only: r8, r4
      use cmn_size, only: IPARW, JPARW, IPAR, JPAR, LPAR, LWEPAR,
     &     IDGRD, JDGRD, NDGRD
      use cloudjx, only: NQD_, CBIN_, NRAN_
      use utilities, only: ctmExitC

      implicit none

! Input parameter
      real(r8), intent(in)   :: PW(IPARW,JPARW), TW(IPARW,JPARW,LPAR)
     &                       ,ETAA(LPAR+1), ETAB(LPAR+1)
     &                       ,ZOFLEW(LPAR+1,IPARW,JPARW)
      real(r8), intent(in), dimension(IPARW,JPARW,LWEPAR) ::
     &                        CLDFRW, CLDIWCW, CLDLWCW
      real(r8), intent(in)   :: AREAXYW(IPARW,JPARW)
     &                       ,ZDEGI(IDGRD,IPAR),ZDEGJ(JDGRD,JPAR)
      real(r4), intent(in)   :: RAN4(NRAN_)
      integer, intent(in)  :: IRAN0,  IMAP(IDGRD,IPAR),JMAP(JDGRD,JPAR)
      logical, intent(in)  :: LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ

! Output parameter
      real(r8), intent(out)  :: CLDSTORE(LPAR+1,NQD_,IPAR,JPAR)
     &                       ,SWSTORE(NQD_,IPAR,JPAR)
      integer, intent(out) :: TYPSTORE(LPAR+1,IPAR,JPAR)

! local variable
      real(r8), dimension(LPAR,NQD_) :: TAUQCA
      real(r8), dimension(NQD_)      :: WTQCA
      integer,dimension(LPAR)      :: JXQCA
c     logical LCLDRAN  ! .T. means that up to 4 ICAs will be picked and loaded
c     in the QCA arrays, can be used, eg, 3 different ones every hour.

      real(r8)  RANNUM(NQD_), GAREA(NDGRD), AWT(NDGRD)
      real(r8)  PSRF(NDGRD), TEM(LPAR,NDGRD), ZL(LPAR+1,NDGRD)
     &       ,CLDFN(LPAR,NDGRD), BOXIWR(LPAR,NDGRD), BOXLWR(LPAR,NDGRD)
      integer KBOX, I,J, IX,JX, ILNG, JLAT, L, N, IRAN
     &       ,LTOP
!
!-----------------------------------------------------------------------
!
      LTOP  = LWEPAR
      CLDSTORE(:,:,:,:)  = 0._r8
      SWSTORE(:,:,:)   = 0._r8
      TYPSTORE(:,:,:)  = 1

      CLDFN(:,:)  = 0._r8
      BOXLWR(:,:) = 0._r8
      BOXIWR(:,:) = 0._r8
      RANNUM(:) = -1._r8
!
      do J = 1,JPAR
      do I = 1,IPAR

        IRAN = IRAN0 + 2*I + 5*J  ! multiply of yr, day & time are 3, 7 & 11
        do N = 1,3
          IRAN = mod (IRAN, NRAN_) + 1
          RANNUM(N) = RAN4(IRAN)
          if (RANNUM(N).lt.0._r8 .or. RANNUM(N).gt.1._r8)
     &      call ctmExitC('*** check RANSET routine ***')
        enddo
        !// For each column, which in a horizontally degraded grid consists
        !// of KBOX native meteorology columns, we find the properties.
        !// For non-degraded run, IDGRD=1, JDGRD=1 and KBOX=1.
        KBOX   = 0
        do JX = 1,JDGRD
        do IX = 1,IDGRD
          KBOX = KBOX + 1
          ILNG = IMAP(IX,I)
          JLAT = JMAP(JX,J)
          GAREA(KBOX) = AREAXYW(ILNG,JLAT)
          PSRF(KBOX)  = PW(ILNG,JLAT)
          AWT(KBOX)   = ZDEGI(IX,I)*ZDEGJ(JX,J)
          do L = 1,LPAR
            TEM(L,KBOX)  = TW(ILNG,JLAT,L)
            ZL(L,KBOX)   = ZOFLEW(L,ILNG,JLAT)
          enddo
            ZL(LPAR+1,KBOX) = ZOFLEW(LPAR+1,ILNG,JLAT)
          do L = 1,LWEPAR
            CLDFN(L,KBOX)  = CLDFRW(ILNG,JLAT,L)
            BOXIWR(L,KBOX) = CLDIWCW(ILNG,JLAT,L)
            BOXLWR(L,KBOX) = CLDLWCW(ILNG,JLAT,L)
          enddo
        enddo
        enddo
        !// Make quadrature independent atmospheres based on these data
        call QUADCA(LTOP,RANNUM,CLDFN,BOXIWR,BOXLWR,TEM,ZL
     &             ,PSRF,ETAA,ETAB,GAREA,AWT,TAUQCA,JXQCA,WTQCA
     &             ,LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)
        do N = 1,NQD_
          SWSTORE(N,I,J) = WTQCA(N)
          do L = 1,LWEPAR
            CLDSTORE(L,N,I,J) = TAUQCA(L,N)
          enddo
        enddo
        do L = 1,LWEPAR
          TYPSTORE(L,I,J) = JXQCA(L)
        enddo

      enddo
      enddo

      return
      end
!
c
c-----------------------------------------------------------------------
      subroutine QUADCA(LTOP,RANNUM,CLDFN,BOXIWR,BOXLWR,T,ZOFL
     &                 ,PSRF,ETAA,ETAB,GAREA,AWT,TAUQCA,JXQCA,WTQCA
     &                 ,LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)
c-----------------------------------------------------------------------
c      
c---generates the 4 quadrature ICAs (indpendent column atmospheres):
c     1 = clear sky,          total OD  <= 0.5
c     2 = cirrus-hazy,  0.5 < total OD  <= 4.0
c     3 = stratus,      4.0 < total OD  <= 30.
c     4 = cumulus,      30. < total OD  <999.  (drops >999.)
c
c     TAUQCA(L,1;4) = visible optical depth (600nm) in each layer 
c     JXQCA(L) = index for optical proerties, fast-JX indices:
c            4   W_H01 (H1/Deirm.)  R eff= 0.25 G=.094
c            5   W_H04 (H1/Deirm.)  R eff= 1.00 G=1.51
c            6   W_H40 (H1/Deirm.)  R eff=10.0  G=146.4
c            7   W_C02 (C1/Deirm.)  R eff= 3.0  G=19.5
c            8   W_C04 (C1/Deirm.)  R eff= 6.0  G=78.2
c            9   W_C08 (C1/Deirm.)  R eff= 12.  G=301.
c           10   W_C13 (C1/Deirm.)  R eff= 20.  G=472.
c           11   W_L06 (W/Lacis)    R eff= 10.  G=194.
c           12   Ice-H = hexagonal ice cloud (Mishchenko)
c           13   Ice-I = irregular ice cloud (Mishchenko)
c     WTQCA(1:4) = weighting for ICA (WT=0.0 ==> no ICA for that type)
c
c---can use random number sequence to select probabilistically likey ICAs
c     recommend that an intial sequence of pseudo-random numbers be 
c     calculated in p-MAIN in the same way every time (eg, 50,000 numbers)
c     and then recycled in the same way based on day, year, and IJ.
c---LCLDAVG = .true. means average clouds over box
c---LCLDRAN = .true. means that up to 4 ICAs will be picked and loaded
c     in the QCA arrays, can be used, eg, 3 different ones every hour.
      
      use cmn_precision, only: r8
      use cmn_size, only: LPAR, NDGRD
      use cmn_parameters, only: G0
      use cloudjx, only: NQD_

      implicit none

      real(r8), parameter  :: G100 = 100._r8/G0

      integer, intent(in)  :: LTOP
      logical,intent(in)   :: LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ
      real(r8), intent(in)   :: RANNUM(NQD_)
      real(r8), intent(in)   :: PSRF(NDGRD), GAREA(NDGRD), AWT(NDGRD)
     &                       ,ZOFL(LPAR+1,NDGRD)
      real(r8), intent(in), dimension(LPAR,NDGRD) ::
     &                        CLDFN,BOXIWR,BOXLWR,T
      real(r8), intent(in),dimension(LPAR+1) :: ETAA,ETAB

      real(r8), intent(out),dimension(LPAR,NQD_) :: TAUQCA
      integer,intent(out),dimension(LPAR)      :: JXQCA
      real(r8), intent(out),dimension(NQD_)      :: WTQCA

      real(r8), dimension(LPAR)  :: P1000,POFL,DELP,DELZ,VOL,AIR,IWR,LWR 
     &                           ,TAU, TAUL, TAUI, TT, CLDF
      real(r8), dimension(LPAR+1):: POFLE
      real(r8), dimension(LPAR,NDGRD) :: TAUIN, TAULN, TAUN
      real(r8)  REFF,QEXT, IWC
      integer L,JXTYPE, N
      integer JXCLD(LPAR,NDGRD)

c---set up profile of pressures, altitudes and volumes:
c--- R = 287.*(1-Q) + 461.5*Q -- assume 0.5% bl w.v.==> R = 288.
c--- delta-z (m) = dln(P) * R * T / g   where R/g = 288/9.81 = 29.36
c--- VOLume (m^3)
      do N = 1,NDGRD

        !// Loop trhough all native columns covered by a degraded horizontal
        !// grid. For non-degraded grid, NDGRD is 1.

        POFLE(1)  = PSRF(N)
        do L = 2,LPAR+1
          POFLE(L)  = ETAA(L) + ETAB(L)*PSRF(N)
          POFL(L-1) = 0.5_r8*(POFLE(L-1)+POFLE(L))
          DELP(L-1) = POFLE(L-1) - POFLE(L)
          DELZ(L-1) = ZOFL(L,N) - ZOFL(L-1,N)
          VOL(L-1) = GAREA(N)*DELZ(L-1)            ! m3 in box
          AIR(L-1) = DELP(L-1)*G100*GAREA(N)       ! kg in box
        enddo

c---P1000 = mid-layer pressure=POFL @1000 hPa>>>>fix to match old interp by L
        do L=1,LTOP
          P1000(L) = 0.5_r8*(ETAA(L)+ETAA(L+1))
     &              +0.5e3_r8*(ETAB(L)+ETAB(L+1))
        enddo

        JXCLD(:,N)  = 1
        do L=1,LTOP
c---TAUIN/LN is the in-cloud optical depth, not average over box.
          TAUIN(L,N) = 0._r8
          TAULN(L,N) = 0._r8

c---BOXLWR = grid-box average liquid water ratio (kg-liq/kg-air-in-whole-box)
c---BOXIWR = grid-box avg ice water ratio (kg/kg)
c---LWR = liquid water ratio (kg-liq/kg-air in-cloud), IWR = ice water ratio
          LWR(L) = 0._r8
          IWR(L) = 0._r8
          if (CLDFN(L,N) .gt. 0._r8)  then

c---ICE cloud OD from non-linear fit to obs:
            if (BOXIWR(L,N) .gt. 0._r8) then
              IWR(L) = BOXIWR(L,N)/CLDFN(L,N)
              IWC = 1.e3_r8*IWR(L)*AIR(L)/VOL(L)         ! IWC = g/m3 ???CF factor
c---  old Sun & Shine 1995 formula, but the IWR should be IWC (g/m3)
c           TAUIN(L,N) = ( 1.d3*IWR(L)/
c    &                 (3.06d-2 + 2.5481d-1*1.d3*IWR(L)) )*
c    &                (7._r8*dlog(POFLE(L)/POFLE(L+1)) )
c           TAUIN(L,N) = ( 1.d3*IWC/
c    &                 (3.06d-2 + 2.5481d-1*1.d3*IWC) )*
c    &                (7._r8*dlog(POFLE(L)/POFLE(L+1)) )
c---  new Heymsfield 2003, log-log fit ext(/m) to Fig B1(a), p.1389
c           IWC = 1.d3*IWR(L)*AIR(L)/VOL(L)         ! IWC = g/m3
              TAUIN(L,N) = DELZ(L) * 9.e-3_r8 * (IWC)**0.752_r8
              if (T(L,N) .ge. 233.15_r8) then
                JXCLD(L,N) = 13  ! ice irreg
              else
                JXCLD(L,N) = 12  ! ice hexag (cold)
              endif
            endif

c---LIQ cloud OD based on R-eff and Q-ext in optical (600 nm)
            if (BOXLWR(L,N) .gt. 0._r8) then
              LWR(L) = BOXLWR(L,N)/CLDFN(L,N)
c---   tau = Qext / (4/3)*Reff(microns)*0.001    *   LWR*delta-P*100/g
              call OD_LIQ (P1000(L),REFF,QEXT,JXTYPE)
              TAULN(L,N) = 0.75d3*LWR(L)*DELP(L)*G100*QEXT/REFF
              !// Overwrite cloud type if liquid cloud dominates over ice
              if (TAULN(L,N).gt.TAUIN(L,N))  JXCLD(L,N) = JXTYPE
            endif

          endif

c---sum the layer ODs, let the liquid cloud JX type prevail (last set).
          TAUN(L,N) = TAULN(L,N) + TAUIN(L,N)
        enddo

      enddo

c-----------------------------------------------------------------------
      if (.not. LCLDQMD .or. .not.LCLDRANA) then
         !// LCLDQMD  = Mid point of quadrature cloud cover ICAs
         !// LCLDQMN  = Mean quadrature cloud cover ICAs
         !// LCLDRANA = Random selected from all cloud cover ICAs
         !// LCLDRANQ = Random selected from 4 Mean quadrature cloud cover ICAs
         !// I.e. for LCLDAVG, LCLDQMN or LCLDRANQ
        JXQCA(:) = 1 !// index for optical proerties, fast-JX indices
        if (NDGRD .eq. 1)  then
          !// Native horizontal resolution
          do L=1,LTOP
            JXQCA(L) =  JXCLD(L,1)
          enddo
        else
          !// Degraded horizontal resolution
          do L = 1,LTOP
            CLDF(L) = 0._r8
            TAUL(L) = 0._r8
            TAUI(L) = 0._r8
            TT(L) = 0._r8
            !// Make average out of the native columns, weighting them
            !// with AWT, the fraction of the simulated grid area (if it
            !// covers 4 native boxes, it will be 0.25).
            do N = 1,NDGRD
              CLDF(L) = CLDF(L) + CLDFN(L,N)*AWT(N)
              TAUL(L) = TAUL(L) + TAULN(L,N)*AWT(N)
              TAUI(L) = TAUI(L) + TAUIN(L,N)*AWT(N)
              TT(L) = TT(L) + T(L,N)*AWT(N)
            enddo
            !// Check again if ice dominates in the average.
            if (CLDF(L) .gt. 0._r8)  then
              if (TAUI(L) .gt. TAUL(L))  then
                if (TT(L) .ge. 233.15_r8)  then
                  JXQCA(L) = 13
                else
                  JXQCA(L) = 12
                endif
              endif
            endif
          enddo
        endif
      endif
c-----------------------------------------------------------------------
      if (LCLDAVG) then
        WTQCA(:) = 0._r8
        TAUQCA(:,:) = 0._r8
        WTQCA(1) = 1._r8     !// Will only use one column
        !// And sum up all optical depths (OD) weighted by cloud
        !// fraction (CF), i.e. make an average so the cloud covers
        !// the whole grid area.
        if (NDGRD .eq. 1)  then
          !// Native horizontal resolution
          do L=1,LTOP
            TAUQCA(L,1) = TAUN(L,1)*CLDFN(L,1)
          enddo
        else
          !// Degraded horizontal resolution
          do L = 1,LTOP
            !// Must use TAUI and TAUL generated after averaging above.
            TAU(L) = TAUI(L) + TAUL(L)
            !// Make average by summing up OD weighted by CF
            TAUQCA(L,1) = TAU(L)*CLDF(L)
          enddo
        endif

      else

c---now generate quadrature column atmospheres (QCAs)
c    if only WTQCA(1) .ne. 0 then have single atmos (AVG or RAN) or all clear

        call CLDQUAD(CLDFN,TAUN,JXCLD,LTOP,RANNUM,TAUQCA,JXQCA,WTQCA
     &              ,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)

      endif

c---now to store the TAUQCA,JXQCA and WTQCA

      return
      end


c-----------------------------------------------------------------------
      subroutine OD_LIQ (PLVL,REFF,QEXT,JXTYPE)
c-----------------------------------------------------------------------
c  using PLEVL = P(L)/P-surf, set R-eff(microns) and Q-ext for liquid clds
c       > 810 hPa use marine stratus:  Reff =  9.60,  Qext = 2.10 (@ 600 nm)
c       < 610 hPa use marine cumulus:  Reff = 12.68,  Qext = 2.08
c  JXTYPE = fast-JX index for scattering and scaling to other wavelengths:
c            4   W_H01 (H1/Deirm.)  R eff= 0.25 G=.094
c            5   W_H04 (H1/Deirm.)  R eff= 1.00 G=1.51
c            6   W_H40 (H1/Deirm.)  R eff=10.0  G=146.4
c            7   W_C02 (C1/Deirm.)  R eff= 3.0  G=19.5
c            8   W_C04 (C1/Deirm.)  R eff= 6.0  G=78.2
c            9   W_C08 (C1/Deirm.)  R eff= 12.  G=301.
c           10   W_C13 (C1/Deirm.)  R eff= 20.  G=472.
c           11   W_L06 (W/Lacis)    R eff= 10.  G=194.
c           12   Ice-H = hexagonal ice cloud (Mishchenko)
c           13   Ice-I = irregular ice cloud (Mishchenko)
c  
      use cmn_precision, only: r8
      implicit none
      real(r8), intent(in)  :: PLVL
      real(r8), intent(out) :: REFF, QEXT
      integer, intent(out):: JXTYPE
      real(r8)   ::  F1

      if (PLVL .gt. 810._r8) then
        REFF = 9.60_r8
        QEXT = 2.10_r8
      elseif (PLVL .lt. 610._r8) then
        REFF = 12.68_r8
        QEXT = 2.08_r8
      else
        F1 = 0.005_r8 * (PLVL - 610._r8)
        REFF = 12.68_r8*(1._r8-F1) + F1*9.60_r8
        QEXT =  2.08_r8*(1._r8-F1) + F1*2.10_r8
      endif
        JXTYPE = 9 
      return
      end


c-----------------------------------------------------------------------
      subroutine CLDQUAD(CLDF,TAU,JXT,LTOP,RANNUM,TAUQCA,JXQCA,WTQCA
     &                  ,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)
c-----------------------------------------------------------------------
c---input is cloud fraction CLDF(1:LTOP), in-cloud optical depth TAU(1:LTOP),
c        and Mie/fast-JX type of cloud JXT(1:LTOP).  
c   generates Indep Colm Atmospheres (ICA) from max-ran overlap criteria
c        LCLDRAN = uses RANNUM(1:3) to select an ICA from cum prob of total OD
c   else generate 4 quadrature colm atmos QCA per Neu et al JGR 2007.
c---output is optical depth in each layer TAUQCA(1:LTOP,1:4), cloud type 
c        JXQCA(1:LTOP), and weight for each atmosphere WTQCA(1:4)
c---note the QCA are in order and hence can have WTQCA(3)=0 while WTQCA(4)=.5

      use cmn_precision, only: r8
      use cmn_size, only: LPAR, NDGRD
      use cloudjx, only: NQD_, CBIN_, ICA_

      implicit none

c---need to reconcile these, put in params.h and remove excess from params.h
      real(r8), parameter  :: OD_QUAD(NQD_) =
     &                         [0.5_r8, 4.0_r8, 30._r8, 1.e9_r8]

      real(r8), intent(in),dimension(LPAR,NDGRD) :: CLDF,TAU
      integer,intent(in),dimension(LPAR,NDGRD) :: JXT
      integer,intent(in) :: LTOP
      logical,intent(in) :: LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
      real(r8), intent(in) :: RANNUM(NQD_)

      real(r8), intent(out),dimension(LPAR,NQD_) :: TAUQCA
      integer,intent(out),dimension(LPAR)      :: JXQCA
      real(r8), intent(out),dimension(NQD_)      :: WTQCA

      real(r8)  TAUSUM(LPAR), TAUBAK(LPAR,NQD_)
      real(r8)  ZNDGRD,GWT,ODCOL,WTCOL,WTSUM
      integer I, II, IG, J, K, L, M, N, KBOX
     &       ,NCOLTL, NRICA
      logical LQMEAN

      integer, dimension(NQD_)  :: NDXQCA
      integer, dimension(NDGRD) ::  NRG, NCOL
      integer, dimension(LPAR,NDGRD) :: GBOT,GTOP,GNR,NCLDF
      integer, dimension(LPAR,CBIN_+1,NDGRD) :: GFNR
      real(r8),  dimension(LPAR,CBIN_+1,NDGRD) :: GFCOL

c---Independent Column Atmosphere data, only need total OD and wts
      real(r8),  dimension(NDGRD*ICA_):: OCOL,WCOL,OCDFS
      integer, dimension(NDGRD*ICA_):: OSORT
c-----------------------------------------------------------------------
      WTQCA(:) = 0._r8
      TAUQCA(:,:) = 0._r8
      LQMEAN = LCLDQMN .or. LCLDRANQ

      call ICANR(CLDF,GFCOL,NRG,GNR,GBOT,GTOP,NCLDF,GFNR,NCOL,LTOP)
      !// NRG: number of max-overlapping groups
      !// GNR: Number of CF bins for each overlap group
      !// NCLDF: Index of CF bin (1:CBIN_) for a given layer
      !// GFNR: Cloud fraction bin for each overlap group (decreasing order)
      !// GFCOL: Its cloud fraction
      !// NCOL: Number of ICAs

c     WTSUM = 0._r8
      K = 0
      do KBOX = 1,NDGRD
        !// Loop through all ICAs
        do I = 1,NCOL(KBOX)
          K  = K + 1
          WTCOL = 1._r8
          ODCOL = 0._r8
          TAUSUM(:) = 0._r8
          II = I
          !// Loop through number of max-overlapping groups
          do N=1,NRG(KBOX)
            !// Index of CF
            IG = mod(II-1, GNR(N,KBOX)) + 1
            II = (II-1)/GNR(N,KBOX) + 1
            !// CF interval, i.e. weight
            GWT = GFCOL(N,IG,KBOX) - GFCOL(N,IG+1,KBOX)
            WTCOL = WTCOL*GWT
            do L=GBOT(N,KBOX),GTOP(N,KBOX)
              if (NCLDF(L,KBOX) .ge. GFNR(N,IG,KBOX)) then
                !// CF bin .ge. CF bin for overlap group
                !// Sum up OD and TAU
                ODCOL = ODCOL + TAU(L,KBOX)
                TAUSUM(L) = TAU(L,KBOX)
              endif
            enddo
          enddo
          WCOL(K) = WTCOL !// Weight of each ICA
          OCOL(K) = ODCOL !// OD column of each ICA
c         WTSUM = WTSUM + WTCOL
          if (LQMEAN) then
            !// Add (weighted) to the corresponding quadrature column (QCA),
            !// i.e. making means for the 4 QCA
            M = 1
            do while (ODCOL .gt. OD_QUAD(M))
              M = M + 1
            enddo
            do L = 1,LTOP
              TAUQCA(L,M) = TAUQCA(L,M) + WTCOL*TAUSUM(L)
            enddo
            WTQCA(M) = WTQCA(M) + WTCOL
          endif
        enddo
      enddo
      NCOLTL = K

      if (LQMEAN) then
        if (NDGRD .gt. 1) then
          !// Degraded horizontal resolution
          ZNDGRD = 1._r8 / real(NDGRD, r8)
          do N = 1,NQD_
            do L = 1,LTOP 
              TAUQCA(L,N) = ZNDGRD * TAUQCA(L,N)
            enddo
            WTQCA(N) = ZNDGRD * WTQCA(N)
          enddo
        endif

        do M = 1,NQD_
          if (WTQCA(M) .gt. 0._r8) then
            do L = 1,LTOP 
              TAUQCA(L,M) = TAUQCA(L,M) / WTQCA(M)
            enddo
          endif
        enddo

        if (LCLDQMN)  return  ! *** LCLDQMN=T  exit

        do N = 1,NQD_
          do L = 1,LTOP
            TAUBAK(L,N) = TAUQCA(L,N)
          enddo
        enddo
        TAUQCA(:,:) = 0._r8

        OCDFS(1) = WTQCA(1)
        do N = 2,NQD_
          OCDFS(N) = OCDFS(N-1) + WTQCA(N)
        enddo
        WTQCA(:) = 0._r8
        N  = 1
        !// Random pick of one of the 4 mean quadrature cloud cover ICAs 
        do while (N.lt.4)
          I  = 1
          do while (OCDFS(I) .lt. RANNUM(N))
            I  = I + 1
          enddo
          NDXQCA(N) = I
          WTQCA(N) = 1._r8
          NRICA = N
          N = N+1
c         write(6,'(A,2I4,2F9.5)') 'NDXQCA',NDXQCA(NRICA),NRICA
c    &           ,RANNUM(NRICA),OCDFS(I)
        enddo
        do N = 1,NRICA
          if (WTQCA(N) .gt. 0._r8) then
            do L = 1,LTOP  
              TAUQCA(L,N) = TAUBAK(L,NDXQCA(N))
            enddo 
          endif
        enddo
        return  ! *** LCLDRANQ=T  exit
      endif


c---choice of NQD_ quadrature QCAs based on break points or random from CDF
      NDXQCA(:) = 0
      JXQCA(:)  = 1
      if (LCLDRANA) then
c---select random ICA based on cumulative prob fn of wts, up to 4 (NQD_)
        OSORT(1) = 1
        OCDFS(1) = WCOL(1)
        do I = 2,NCOLTL
          OCDFS(I) = OCDFS(I-1) + WCOL(I)
          OSORT(I) = I
        enddo
        N  = 1
        do while(N .lt. 4)
          I = 1
          do while (OCDFS(I) .lt. RANNUM(N))
            I = I + 1
          enddo
          NDXQCA(N) = I
          WTQCA(N) = 1._r8
          NRICA = N
          N = N+1
c         write(6,'(A,2I4,2F9.5)') 'NDXQCA',NDXQCA(NRICA),NRICA
c    &           ,RANNUM(NRICA),OCDFS(I)
        enddo
      else
        call QUADMD(OCOL,WCOL,OD_QUAD,WTQCA,NDXQCA,OSORT,NCOLTL)
        NRICA  = NQD_
      endif

c---locate representive ICAs and put into QCAs:  TAUQCA & JXQCA
      do J = 1,NRICA
        I = NDXQCA(J)
        if (I.gt.0 .and. I.le.NCOLTL) then
          II = OSORT(I)
          KBOX = 1
          K    = NCOL(KBOX)
          do while (II .gt. K)
            KBOX = KBOX + 1
            K    = K + NCOL(KBOX)
          enddo
          N   = 0
          do K = 1,KBOX-1
            N = N + NCOL(K)
          enddo
          II  = II - N
          do N = 1,NRG(KBOX)
            IG = mod(II-1, GNR(N,KBOX)) + 1
            II = (II-1)/GNR(N,KBOX) + 1
            do L=GBOT(N,KBOX),GTOP(N,KBOX)
              if (NCLDF(L,KBOX) .ge. GFNR(N,IG,KBOX)) then
                TAUQCA(L,J) = TAU(L,KBOX)
                JXQCA(L) = JXT(L,KBOX)
              endif
            enddo
          enddo
        endif
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine QUADMD(OCOL,WCOL,OD_QUAD,WTQCA,NDXQCA,OSORT,NCOLTL)
c-----------------------------------------------------------------------

      use cmn_precision, only: r8
      use cmn_size, only: NDGRD
      use cloudjx, only: NQD_, ICA_

      implicit none

c---need to reconcile these, put in params.h and remove excess from params.h
      integer, intent(in) :: NCOLTL
      real(r8), intent(in)  :: OD_QUAD(NQD_)
      real(r8), intent(in), dimension(NDGRD*ICA_):: WCOL, OCOL

      real(r8), intent(out)  :: WTQCA(NQD_)
      integer, intent(out) :: NDXQCA(NQD_)
      integer, intent(out) :: OSORT(NDGRD*ICA_)

      integer   NQ1(NQD_),NQ2(NQD_), I, N, N1,N2

c---Independent Column Atmosphere data, only need total OD and wts
      real(r8),  dimension(NDGRD*ICA_):: OCOLS, OCDFS
      real(r8)   WWW4(NQD_)
c-----------------------------------------------------------------------
      WTQCA(:) = 0._r8

c---sort all the Indep Column Atmos (ICAs) in order of increasing OD
      if (NCOLTL .eq. 1)  then
        OSORT(1) = 1
        OCOLS(1) = OCOL(1)
      else
        call HEAPSRT (NCOLTL,OCOL,OCOLS,OSORT,NDGRD*ICA_)
      endif
        OCDFS(1) = WCOL(OSORT(1))
      do I = 2,NCOLTL
        OCDFS(I) = OCDFS(I-1) + WCOL(OSORT(I))
      enddo

c---choice of NQD_ quadrature QCAs based on break points or random from CDF
        NDXQCA(:) = 0
c---find beginning/end of quad range, note NQ2 < NQ1 means nothing in that range
        I = 1
        do N = 1,NQD_
          do while (OCOLS(I) .lt. OD_QUAD(N) .and. I.le.NCOLTL)
            I = I+1
          enddo
          NQ2(N) = I-1
        enddo
        NQ1(1) = 1
        do N = 2,NQD_
          NQ1(N) = NQ2(N-1)+1
        enddo

c---define ICA wts from cum prob, pick middle ICA as representative
        NDXQCA(:) = 0
        WTQCA(:) = 0
        do N = 1,NQD_
          N1 = NQ1(N)
          N2 = NQ2(N)
          if (N2 .ge. N1) then
            NDXQCA(N) = (N1+N2)/2
            if (N1 .gt. 1) then
              WTQCA(N) = OCDFS(N2)-OCDFS(N1-1)
            else
              WTQCA(N) = OCDFS(N2)
            endif
          endif
        enddo
c---cuts off total OD > OD_QUAD(NQD_), must renormalize ICA weights
        if (OCDFS(NQ2(NQD_)) .gt. 0._r8)  then
          do N = 1,NQD_
            WTQCA(N) = WTQCA(N)/OCDFS(NQ2(NQD_))
          enddo
        endif

      return
      end


c-----------------------------------------------------------------------
      subroutine ICANR(CLDF,GFCOL,NRG,GNR,GBOT,GTOP,NCLDF,GFNR,NCOL
     &                ,LTOP)
c-----------------------------------------------------------------------
c---evaluate number of ICAs: NCOL

c---input is cloud fraction CLDF(1:LTOP)

      use cmn_precision, only: r8
      use cmn_size, only: LPAR, NDGRD
      use cloudjx, only: NQD_, CBIN_, ICA_

      implicit none

c---need to reconcile these, put in params.h and remove excess from params.h
      real(r8), parameter  :: GRP_BRK = 0.049_r8

      real(r8), intent(in),dimension(LPAR,NDGRD) :: CLDF
      integer,intent(in) :: LTOP

      integer, intent(out) ::  NRG(NDGRD), NCOL(NDGRD)
      integer, intent(out), dimension(LPAR,NDGRD) :: GBOT,GTOP,GNR,NCLDF
      integer, intent(out), dimension(LPAR,CBIN_+1,NDGRD) :: GFNR
      real(r8),  intent(out), dimension(LPAR,CBIN_+1,NDGRD) :: GFCOL

      real(r8)   FBIN
      integer, dimension(NDGRD) ::  NRGX, NCOLX
      integer, dimension(LPAR)  :: GMIN,GMAX
      integer, dimension(CBIN_) :: NSAME
      integer  I,L,LL,N,NC, NGRBRK, KBOX
c-----------------------------------------------------------------------
      FBIN = CBIN_
      if (GRP_BRK .gt. 0._r8) then
        NGRBRK = min(int(GRP_BRK*FBIN) + 1, CBIN_)
      else
        NGRBRK = 0
      endif
c       write(6,*) 'NGRP_break',NGRBRK

c---quantize cloud fractions into bins, CBIN_=40 => N=0=[0%], N=1=[0.001-2.5%],
c        N=2=[2.5-5.0%], ... N=39=[95.0-97.5%], N=40=[97.5-100%]
      do KBOX = 1,NDGRD

        do L = 1,LTOP
          !// NCLDF is index for the bins, 0 is cloud free
          if (CLDF(L,KBOX) .gt. 0._r8) then
            NCLDF(L,KBOX) = min(int(CLDF(L,KBOX)*FBIN) + 1, CBIN_)
          else
            NCLDF(L,KBOX) = 0
          endif
        enddo

c---search from bottom to top, finding 1st level in group with cloud fraction 
c    .ge. threshold, and then first level above that at which the cld fraction
c    is .lt. threshold. NRG = number of such groups.
        L = 1
        NRG(KBOX) = 0
        do while (L.lt.LTOP)
          if (NCLDF(L,KBOX) .gt. NGRBRK) then
            NRG(KBOX) = NRG(KBOX)+1
            GMIN(NRG(KBOX)) = L
            do LL = L+1,LTOP
c  look for first layer to drop below CLDF threshold = NGRBRK
              if (NCLDF(LL,KBOX) .le. NGRBRK) then
                GMAX(NRG(KBOX)) = LL
                goto 11
              endif
            enddo
            GMAX(NRG(KBOX)) = LTOP
   11       continue
            L = GMAX(NRG(KBOX))+1
          else
            L = L+1
          endif
        enddo

c---assign levels between maximum overlap groups to group above.
        if (NRG(KBOX) .eq. 0) then
          !// Cloud free
          NRG(KBOX) = 1
          GBOT(1,KBOX) = 1
          GTOP(1,KBOX) = LTOP
        else
          !// There are clouds
          GBOT(1,KBOX) = 1        !// From surface
          GTOP(1,KBOX) = GMAX(1)  !//   to top of first group
          do N=2,NRG(KBOX)        !//  etc...
            GBOT(N,KBOX) = min(GTOP(N-1,KBOX)+1, GMIN(N))
            GTOP(N,KBOX) = GMAX(N)
          enddo
          GTOP(NRG(KBOX),KBOX) = LTOP
        endif

c---for each max-overlap group calculate number of unique cloud fractions
        do N=1,NRG(KBOX)
          NSAME(:) = 0
          do L=GBOT(N,KBOX),GTOP(N,KBOX)
            if (NCLDF(L,KBOX) .gt. 0) then
              !// Flag whether we have cloud in this cloud fraction bin
              NSAME(NCLDF(L,KBOX)) = 1
            endif
          enddo

c---sort cloud fractions in deceasing order for each max-overlap group
c---  note that largest bin N=CBIN_ (eg, 97.5%) will be treated as 100%
          !// GFNR: Cloud fraction bin for first overlap group
          !// GFCOL: cloud fraction
          GFNR(N,1,KBOX) = CBIN_
          GFCOL(N,1,KBOX) = 1.0_r8
          NC=1
          do I=CBIN_-1,1,-1
            if(NSAME(I) .gt. 0) then
              NC = NC+1
              !// Cloud fraction bin for this overlap group
              GFNR(N,NC,KBOX) = I
              !// Its cloud fraction (for I=20,CBIN_=20: 0.95)
              GFCOL(N,NC,KBOX) = (float(GFNR(N,NC,KBOX))-0.5_r8) / FBIN
            endif
          enddo
          GNR(N,KBOX) = NC  !// Number of CF bins for this overlap group (N)
          GFNR(N,NC+1,KBOX) = 0      !// Its CF bin
          GFCOL(N,NC+1,KBOX) = 0._r8 !// and CF
        enddo

c---number of unique columns in group:  if too many ICAs, drop upper groups!
        NCOL(KBOX) = 1
        do N=1,NRG(KBOX)
          NCOL(KBOX) = NCOL(KBOX)*GNR(N,KBOX)
          if (NCOL(KBOX) .le. ICA_) then
            NCOLX(KBOX) = NCOL(KBOX)
            NRGX(KBOX) = N
          endif
        enddo
        if (NCOL(KBOX) .gt. ICA_) then
          write(6,'(a,5i8)') ' NCOL greater than ICA_'
     &           ,NCOL(KBOX),ICA_,NCOLX(KBOX),NRG(KBOX),NRGX(KBOX)
          NCOL(KBOX) = NCOLX(KBOX)
          NRG(KBOX) = NRGX(KBOX)
        endif

      enddo    ! KBOX = 1,NDGRD

      return
      end


c-----------------------------------------------------------------------   
      subroutine HEAPSRT (N,A,AX,IX,ND)
c-----------------------------------------------------------------------
c  classic heapsort, sorts real*8 array A(N) into ascending order, 
c     places sorted array AX(N):   AX(1) .le. AX(N)
c     returns indexing IX(N) that records the location of A in sequence:
c           A(IX(J)) ==> AX(J), s.t. IX(1) = orig location of smallest A
c                           and IX(N) = original loc. of largest A
c
      use cmn_precision, only: r8
      implicit none
      integer, intent(in)  :: N, ND
      real(r8), dimension(ND),intent(in)  :: A

      real(r8), dimension(ND),intent(out) :: AX
      integer,dimension(ND),intent(out) :: IX

      integer :: I,J,L,IR,IA
      real(r8) :: RA

      do I=1,N
        IX(I) = I
        AX(I) = A(I)
      enddo
        L  = N/2+1
        IR = N
   10 continue
      if (L .gt. 1) then
         L = L-1
         RA = AX(L)
         IA = IX(L)
      else
         RA = AX(IR)
         IA = IX(IR)
         AX(IR) = AX(1)
         IX(IR) = IX(1)
         IR = IR-1
        if (IR .eq. 1) then
          AX(1) = RA
          IX(1) = IA
          return
        endif
      endif
        I = L
        J = L+L
   20 continue
      if (J .le. IR) then
        if (J .lt. IR) then
          if(AX(J) .lt. AX(J+1)) then
            J = J+1
          endif
        endif
        if (RA .lt. AX(J)) then
          AX(I) = AX(J)
          IX(I) = IX(J)
          I = J
          J = J+J
        else
          J = IR+1
        endif
        goto 20
      endif
        AX(I) = RA
        IX(I) = IA
      goto 10
      end



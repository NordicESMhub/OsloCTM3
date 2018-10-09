!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, June 2015
!//=========================================================================
!// Large scale scavenging routines.
!//=========================================================================
module scavenging_largescale_uci
  !// ----------------------------------------------------------------------
  !// MODULE: scavenging_largescale_uci
  !// DESCRIPTION: Large scale wet scavenging (Neu ...)
  !//              Converted from UCI p-scav.f and modified for CTM3.
  !//
  !// Contains
  !//   subroutine WASHO
  !//   subroutine WASH1
  !//   subroutine WASH2
  !//   subroutine getHstar
  !//   subroutine DISGAS2
  !//   subroutine WASHGAS2
  !//   subroutine DISGAS
  !//   subroutine HENRYS
  !//   subroutine HENRYSallT
  !//   subroutine RAINGAS
  !//   subroutine WASHGAS
  !//   subroutine DIAMEMP
  !//   subroutine GAMMAX
  !//   subroutine WETSET_CTM3
  !//
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Parameters
  real(r8), parameter  :: TICE   = 273._r8
  real(r8), parameter  :: TSOLID = 258._r8 ! Also referred to as TMIX
  !// ----------------------------------------------------------------------
  character(len=*), parameter :: f90file = 'scavenging_largescale_uci.f90'
  !// ----------------------------------------------------------------------
  !//---(p-scav.f)    p-code 5.6 (Neu, Tang, Prather Aug 2008)
  !//Version: qcode_56d; 20090318
  !//     CTM3: modified VDEP
  !//         SCAV_FMT is longer
  !//         TCCNVHENRY is read
  !// ----------------------------------------------------------------------
  private
  public WASHO, WETSET_CTM3, getHstar
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine WASHO(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
       AIRB,BTEM,DELT,MP)
    !// --------------------------------------------------------------------
    !// ---p-scav 5.6 (2008  J. Neu, MJP)
    !// ---called from pmain to calculate rainout and washout of tracers
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: NSCX
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: AIRB, BTEM
    real(r8), intent(in) :: DELT
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'WASHO'
    !//---------------------------------------------------------------------

    !// NSCX = 0 = simpel e-fold vs. ht per tracer
    if (NSCX .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': DO NOT USE THIS ROUTINE in Oslo CTM3'
       write(6,'(a)') '    : NSCX should be 1'
       stop
       call WASH1 (BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,DELT,MP)
    end if

    !//NSCX = 1 = new water/ice scavenging (J.Neu 8/2008)
    if (NSCX .eq. 1) then
       call WASH2 (BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, & 
            AIRB,BTEM,DELT,MP)
    end if

    !// --------------------------------------------------------------------
  end subroutine WASHO
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine WASH1(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,DTSCAV,MP)
    !// --------------------------------------------------------------------
    !//---p-scav 5.6 (2008)
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LCONVM, LWEPAR, LPAR
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, NTM
    use cmn_chem, only: TRWETL
    use cmn_met, only: ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: DTSCAV
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    !// Locals
    integer :: I,II,J,JJ,L,N,LTOP,LZ
    real(r8) :: FXYW1
    !// --------------------------------------------------------------------

    !// NSCX = 0 = simpel e-fold vs. ht per tracer
    !// --------------------------------------------------------------------

    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ   = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II   = I - MPBLKIB(MP) + 1

        !// ---globally uniform "rainout" loss rate (/s) = e-fold vs. height
        !//    per tracer
        !// ---TRWETL(LZ=1:9) stores rate from 1 to 8 km, and then above
        !//    8 km (LZ=9)
        !// ---applied to both aerosols and gases when NSCX = 0
        !// ---top level for rain/wash out scavenging 
        !//    (possibly include LPAUZ?, hence inside IJ)
        LTOP = min(LCONVM, LWEPAR, LPAR-1)
        do L = 1, LTOP
          LZ = max(1,min(9,&
               (int(0.5e-3_r8 * (ZOFLE(L,I,J) + ZOFLE(L+1,I,J))))))

          do N = 1, NTM
            if (TRWETL(LZ,N) .gt. 0._r8) then
              FXYW1 = exp(-DTSCAV * TRWETL(LZ,N))
              BTT(L,N,II,JJ) = BTT(L,N,II,JJ) * FXYW1
              BZT(L,N,II,JJ) = BZT(L,N,II,JJ) * FXYW1
              BZZ(L,N,II,JJ) = BZZ(L,N,II,JJ) * FXYW1
              BXZ(L,N,II,JJ) = BXZ(L,N,II,JJ) * FXYW1
              BYZ(L,N,II,JJ) = BYZ(L,N,II,JJ) * FXYW1
              BXT(L,N,II,JJ) = BXT(L,N,II,JJ) * FXYW1
              BXX(L,N,II,JJ) = BXX(L,N,II,JJ) * FXYW1
              BYT(L,N,II,JJ) = BYT(L,N,II,JJ) * FXYW1
              BYY(L,N,II,JJ) = BYY(L,N,II,JJ) * FXYW1
              BXY(L,N,II,JJ) = BXY(L,N,II,JJ) * FXYW1
            end if
          end do

        end do

      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine WASH1
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine WASH2(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
    AIRB,BTEM,DTSCAV,MP)
    !// --------------------------------------------------------------------
    !// ---p-scav 5.6 (2008  J. Neu, MJP, Q. Tang)
    !//---called from WASHO (from MAIN) to calc rainout and washout of tracers
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LCONVM, LWEPAR, LPAR
    use cmn_ctm, only: AREAXY, ETAA, ETAB, &
         MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, NTM
    use cmn_chem, only: TRWETL, TCWETL, TCHENA, TCHENB, TCKAQA, &
         TCKAQB, SCV_RETEFF, SCViceFR,  SCV_T258, TMASS
    use cmn_met, only: ZOFLE, P, PRECLS, CLDLWC, CLDIWC, CLDFR
    use utilities, only: CIWMIN
    use cmn_oslo, only: TCCNVHENRY, chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: DTSCAV
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: AIRB, BTEM
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ

    !// Locals
    integer :: I,J,L,N, LTOP, II, JJ
    real(r8), dimension(LPAR+1) :: POFLE
    real(r8), dimension(LPAR) :: DELZ,POFL
    real(r8), dimension(LPAR) :: RLS, RMC, CLWC,CIWC,CFR,TEM
    real(r8), dimension(LPAR) :: QM
    real(r8), dimension(LPAR,NPAR) :: QTT, QTTNEW
    real(r8), dimension(LPAR,NPAR) :: QTTnoSCAV
    real(r8) :: GAREA, FXYW1 
    integer :: LWSHTYP(NPAR),LICETYP(NPAR)

    real(r8) :: CLWX,CFXX
    real(r8) :: RNEW, RPRECIP, DELTARIMEMASS, DELTARIME
    real(r8) :: MASSLOSS
    real(r8) :: DOR, DNEW, DEMP, COLEFFSNOW, RHOSNOW
    real(r8) :: WEMP, REMP, RRAIN, RWASH
    real(r8) :: QTPRECIP, QTCXA, QTAX
    real(r8) :: QTDISCF, QTDISRIME, QTDISCXA
    real(r8) :: QTDISSTAR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'WASH2'
    !//---------------------------------------------------------------------

    !// Cloud overlap nomenclature
    !// --------------------------------------------------------------------
! TOP=above | 1 = rain under cloud  | 2 = rain from clear | (3) = no rain      |
!-------------------------------------------------------------------------------
! lyr L     | 4=cld core(under rain) | 5=cld(not/rain) | 6=rain thru clr | (7) |
!-------------------------------------------------------------------------------
! BOT=below |  8 = rain from under cloud in lyr L    | 9 = rain (ambient)| (10)|
!-------------------------------------------------------------------------------
!
! notes:
!       4 = cloud under rainfall (including cloud above and ambient above)
!       5 = cloud with no rain coming in, can initiate rain.
!       2 & 6 & 9 precip falling thru clear sky (aka ambient)
!           in 6, the area and precip shrink (per km)
!
! rain = all forms of precipitation (kg/s) from each sub-box
!
!  FR1:FR9 = FRaction of grid box area  (0->1)
!  PR1:PR9 = PRecipitation thru bottom of region (kg/s)
!  QT____ = tracer quantity in sub-box, or precip over time step (kg)


    !// Local tracer arrays: QT__ = kg of tracer in region
    !//     NET = net loss from region (1,2,3) in layer L
    !//     TOP = amount coming in top region (1,3)
    !//     PRC = amount in precip regions (0,1,2) 0=1+2
    !//     RIM = amount picked up riming, all cloud regions
    !//     EVP = evaporative losses/gains in cloud
    !//           (1 = all cloud) & ambient (3)
    real(r8), dimension(NPAR) :: QTTOP1, QTTOP2, QTTOP3
    real(r8), dimension(NPAR) :: QTTOP4, QTTOP6
    real(r8), dimension(NPAR) :: QTBOT7, QTBOT8, QTBOT9
    real(r8), dimension(NPAR) :: QTT4, QTT5, QTT6
    real(r8), dimension(NPAR) :: QTNET4, QTNET5, QTNET6
    real(r8), dimension(NPAR) :: QTPRC0, QTPRC4, QTPRC5
    real(r8), dimension(NPAR) :: QTWSH4, QTWSH6, QTRIM4, QTEVP4, QTEVP6
    real(r8), dimension(NPAR) :: QTEVP4X, QTEVP4Y, QTEVP6X, QTEVP6Y

    real(r8), dimension(LPAR) :: CEVAP6, DRLS

    !// Cloud fractions
    real(r8) :: FR1,FR2,FR3, FR4,FR5,FR6, FR7,FR8,FR9
    !// Precipitation rates
    real(r8) :: PR1,PR2,PR3, PR4,PR5,PR6, PR7,PR8,PR9

    real(r8) :: W1TO4,W1TO6, W2TO4,W2TO6, W3TO4,W3TO6
    real(r8) :: D1,D2,D3, D4,D5,D6, D7,D8,D9
    real(r8) :: PR4T, D4T, PR6T, D6T, FR6T
    real(r8) :: PR6TMP

    !// For Henry expressions
    real(r8) :: HSTAR
    integer, parameter :: HCALLER = 1 !// caller flag for getHstar

    !// Parameters
    real(r8), parameter  :: &
         TMIX       = TSOLID, &
         TFROZ      = 240._r8, &
         EVAPRATE   = 0.25_r8, &   !fraction/km
         CFMIN      = 0.025_r8, &
         CFCOREMIN  = 0.1_r8, &
         DMIN       = 0.1_r8, &       !mm
         VOLPOW     = 1._r8/3._r8, &
         RHORAIN    = 1000._r8, &     !kg/m3
         RHOSNOWFIX = 100._r8, &      !kg/m3
         COLEFFRAIN = 0.7_r8, &
         COLEFFAER  = 0.05_r8, &
         QTTMIN     = 1.e-20_r8    ! Min QTT to scavenge (tiny QTT may -> NAN)
    !// --------------------------------------------------------------------


    !// --------------------------------------------------------------------
    !// NSCX = 1 = new water/ice scavenging (J.Neu 8/2008)
    !// --------------------------------------------------------------------

    !// top level for rain/wash out scavenging (possibly include LPAUZ(I,J)?)
    LTOP = min(LCONVM, LWEPAR, LPAR-1)    ! possibly min (..., LPAUZB(I,J)) ?
    QTT(:,:) = 0._r8
    QTTNEW(:,:) = 0._r8
    QTTnoSCAV(:,:) = 0._r8

    !//---set up flags for each species:
    do N = 1, NTM
       LWSHTYP(N) = 0
       LICETYP(N) = 0
       !//-----check if soluble
       if(TCHENA(N) .gt. 0._r8) then
          !// check whether mass-limited or henry's law
          if(TCKAQB(N) .gt. 0._r8) then
             LWSHTYP(N) = 1
          else
             LWSHTYP(N) = 2
          end if
          !// check whether soluble in ice
          if(SCViceFR(N) .gt. 0._r8) then
             !// Could eventually check on SCV_RETEFF(N) also, but not
             !// until it is listed in input file.
             LICETYP(N) = 1
          else
             LICETYP(N) = 2
          end if
       end if
       !// CTM3: possibility of overriding large scale scavenging, so that
       !//       only convective scavenging is calculated
       if (TCCNVHENRY(N) .eq. 5) then
          LWSHTYP(N) = 0
          LICETYP(N) = 0
       else if (TCCNVHENRY(N) .eq. 6) then
          LWSHTYP(N) = 3 !// Skip removal on rain
          LICETYP(N) = 1 !// Include removal on ice
       else if (TCCNVHENRY(N) .eq. 7) then
          LWSHTYP(N) = 3 !// Skip removal on rain
          LICETYP(N) = 1 !// Include removal on ice
       end if
    end do


    !// Loop through each column
    !// --------------------------------------------------------------------
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ   = J - MPBLKJB(MP) + 1
    do I = MPBLKIB(MP),MPBLKIE(MP)
      II   = I - MPBLKIB(MP) + 1
      !// ------------------------------------------------------------------

      !// set ups from (I,J) grid:
      RLS(:) = 0._r8
      RMC(:) = 0._r8
      CFXX   = 0._r8
      GAREA  = AREAXY(I,J)

      do L = 1, LTOP
        RLS(L)    = PRECLS(I,J,L)        ! kg/s thru bottom of layer L
        !RMC(L)    = PRECCNV(I,J,L) !// Convective rain not used
        CLWC(L)   = CLDLWC(I,J,L)      ! kg/kg in box; LWC/CF=in-cloud
        CIWC(L)   = CLDIWC(I,J,L)
        CFR(L)    = CLDFR(I,J,L)
        TEM(L)    = BTEM(L,II,JJ)
        QM(L)     = AIRB(L,II,JJ)
      end do
      RLS(LTOP+1) = 0._r8


      !// approx height (m) of layer edges (ignore Q)
      POFLE(1)  = P(I,J)
      do L = 2, LPAR+1
        POFLE(L) = ETAA(L) + ETAB(L) * POFLE(1)
        !// R = 287.*(1-Q) + 461.5*Q -- assume 0.5% bl w.v.==> R = 288.
        !// delta-z (m) = dln(P) * R * T / g   where R/g = 288/9.81 = 29.36 
        DELZ(L-1)= ZOFLE(L,I,J) - ZOFLE(L-1,I,J)
      end do
      do L = 1, LPAR
        POFL(L) = 0.5_r8 * (POFLE(L) + POFLE(L+1))
        !fract loss of amb.precip.
        CEVAP6(L) = EVAPRATE * DELZ(L) * 1.e-3_r8
      end do


      do N = 1, NTM
        !// check if soluble
        if(LWSHTYP(N) .gt. 0) then
          do L = 1,LPAR
            !// Allow scavening of a fraction of grid box
            if (TCWETL(N) .eq. 1._r8) then
              QTTnoSCAV(L,N) = 0._r8
              QTT(L,N)    = BTT(L,N,II,JJ)
              QTTNEW(L,N) = BTT(L,N,II,JJ)
            else
              !// Keep out fraction not subject for scavenging
              QTTnoSCAV(L,N) = BTT(L,N,II,JJ) * (1._r8 - TCWETL(N))
              QTT(L,N)    = BTT(L,N,II,JJ) - QTTnoSCAV(L,N)
              QTTNEW(L,N) = BTT(L,N,II,JJ) - QTTnoSCAV(L,N)
            end if
          end do
        end if
      end do !// do N = 1, NTM

      !// initialize input to zero above top layer (LTOP)
      QTBOT7(:) = 0._r8
      QTBOT8(:) = 0._r8
      QTBOT9(:) = 0._r8

      D7 = 0._r8
      D8 = 0._r8
      D9 = 0._r8

      PR7 = 0._r8
      PR8 = 0._r8
      PR9 = 0._r8

      FR7 = 0._r8
      FR8 = 0._r8
      FR9 = 0._r8

      W1TO4 = 0._r8
      W1TO6 = 0._r8
      W2TO4 = 0._r8
      W2TO6 = 0._r8
      W3TO4 = 0._r8
      W3TO6 = 0._r8

      do L = LTOP, 1, -1

        QTTOP1(:) = QTBOT7(:)
        QTTOP2(:) = QTBOT8(:)
        QTTOP3(:) = QTBOT9(:)
        D1  = D7
        D2  = D8
        D3  = D9
        PR1 = PR7
        PR2 = PR8
        PR3 = PR9
        FR1 = FR7
        FR2 = FR8
        FR3 = FR9

        PR4 = 0._r8
        FR4 = 0._r8
        PR5 = 0._r8
        FR5 = 0._r8
        PR6 = 0._r8
        FR6 = 0._r8

        D4 = 0._r8
        D5 = 0._r8
        D6 = 0._r8

        QTTOP4(:) = 0._r8
        QTTOP6(:) = 0._r8

        QTPRC0(:) = 0._r8
        QTPRC4(:) = 0._r8
        QTPRC5(:) = 0._r8

        QTRIM4(:) = 0._r8
        QTWSH4(:) = 0._r8
        QTEVP4(:) = 0._r8
        QTEVP4X(:) = 0._r8
        QTEVP4Y(:) = 0._r8

        QTWSH6(:) = 0._r8
        QTEVP6(:) = 0._r8
        QTEVP6X(:) = 0._r8
        QTEVP6Y(:) = 0._r8

        QTDISCF = 0._r8
        QTDISRIME = 0._r8
        QTDISCXA = 0._r8

        RPRECIP = 0._r8
        DELTARIMEMASS = 0._r8
        DELTARIME = 0._r8
        DOR = 0._r8
        DNEW = 0._r8

        if (RLS(L).le.0._r8 .and. RLS(L+1).le.0._r8) then
          !// No rain out at bottom and no rain in at top
          PR7 = 0._r8
          PR8 = 0._r8
          PR9 = 0._r8
          FR7 = 0._r8
          FR8 = 0._r8
          FR9 = 0._r8
          D7  = 0._r8
          D8  = 0._r8
          D9  = 0._r8
          QTBOT7(:) = 0._r8
          QTBOT8(:) = 0._r8
          QTBOT9(:) = 0._r8
        else if (RLS(L).le.0._r8 .and. RLS(L+1).gt.0._r8) then
          !// No rain out at bottom, but rain in at top
          PR7 = 0._r8
          PR8 = 0._r8
          PR9 = 0._r8
          FR7 = 0._r8
          FR8 = 0._r8
          FR9 = 0._r8
          D7  = 0._r8
          D8  = 0._r8
          D9  = 0._r8
          QTBOT7(:) = 0._r8
          QTBOT8(:) = 0._r8
          QTBOT9(:) = 0._r8
          do N = 1, NTM
            if (LWSHTYP(N) .gt. 0) then 
              if (QTTOP1(N) + QTTOP2(N) + QTTOP3(N) .lt. 0._r8) then
                QTTNEW(L,N) = QTT(L,N)
              else
                QTTNEW(L,N) = QTT(L,N) + QTTOP1(N) + QTTOP2(N) + QTTOP3(N)
              end if
            end if
          end do
        else
          !// Rain out at bottom and rain in at top
          ! Determine if CFMIN is needed
          if (CFR(L).gt.0._r8 .and. CFR(L).lt.CFMIN) then
            CFR(L) = CFMIN
            CLWC(L) = CLWC(L) / CFR(L)
            CIWC(L) = CIWC(L) / CFR(L)
     
            call CIWMIN(TEM(L),CLWC(L),CIWC(L))
     
            CLWC(L) = CLWC(L)*CFR(L)
            CIWC(L) = CIWC(L)*CFR(L)
          end if

          ! Determine if CFCOREMIN is needed
          PR6TMP = (PR1 + PR2 + PR3) * (1._r8 - CEVAP6(L))
          if (RLS(L) .ge. PR6TMP) then
            CFXX = max(CFR(L), CFMIN)
            CLWC(L) = CLWC(L) / CFXX
            CIWC(L) = CIWC(L) / CFXX

            call CIWMIN(TEM(L),CLWC(L),CIWC(L))

            CLWC(L) = CLWC(L) * CFXX
            CIWC(L) = CIWC(L) * CFXX
          else
            CFXX = CFR(L)
          end if
!          if (DRLS(L) .gt. 0._r8) then
!            CFXX = max(CFR(L),CFCOREMIN)
!            CLWC(L) = CLWC(L)/CFXX
!            CIWC(L) = CIWC(L)/CFXX

!            call CIWMIN(TEM(L),CLWC(L),CIWC(L))

!            CLWC(L) = CLWC(L)*CFXX
!            CIWC(L) = CIWC(L)*CFXX
!          else
!            CFXX = CFR(L)
!          endif

          !// Determine the weighting coeffients and fractions of
          !// sub-boxes and D4,D6
          if (FR1+FR2+FR3 .gt. 0._r8) then
            if (FR1+FR2+FR3 .gt. CFXX) then
              FR4 = CFXX
              FR6 = FR1 + FR2 + FR3 - CFXX
              FR5 = 0._r8
              if(FR1.eq.0._r8 .and. FR2.ne.0._r8 .and. FR3.ne.0._r8) then
                W1TO4 = 0._r8
                W1TO6 = 0._r8
                if (FR4 .le. FR2) then
                  W2TO4 = FR4 / FR2 
                  W2TO6 = 1._r8 - W2TO4
                  W3TO4 = 0._r8
                  W3TO6 = 1._r8
                else
                  W2TO4 = 1._r8
                  W2TO6 = 0._r8
                  W3TO6 = FR6 / FR3
                  W3TO4 = 1._r8 - W3TO6
                end if
              else if(FR2.eq.0._r8.and.FR1.ne.0._r8.and.FR3.ne.0._r8) then
                W2TO4 = 0._r8
                W2TO6 = 0._r8
                if (FR4 .le. FR1) then
                  W1TO4 = FR4 / FR1
                  W1TO6 = 1._r8 - W1TO4
                  W3TO4 = 0._r8
                  W3TO6 = 1._r8
                else
                  W1TO4 = 1._r8
                  W1TO6 = 0._r8
                  W3TO6 = FR6 / FR3
                  W3TO4 = 1._r8 - W3TO6
                end if
              else if(FR3.eq.0._r8.and.FR1.ne.0._r8.and.FR2.ne.0._r8) then
                W3TO4 = 0._r8
                W3TO6 = 0._r8
                if (FR4 .le. FR1) then
                  W1TO4 = FR4 / FR1
                  W1TO6 = 1._r8 - W1TO4
                  W2TO4 = 0._r8
                  W2TO6 = 1._r8
                else
                  W1TO4 = 1._r8
                  W1TO6 = 0._r8
                  W2TO6 = FR6 / FR2
                  W2TO4 = 1._r8 - W2TO6
                end if
              else if (FR1.eq.0._r8 .and. FR2.eq.0._r8) then
                W1TO4 = 0._r8
                W1TO6 = 0._r8
                W2TO4 = 0._r8
                W2TO6 = 0._r8
                W3TO4 = FR4 / FR3
                W3TO6 = 1._r8 - W3TO4
              else if (FR1.eq.0._r8 .and. FR3.eq.0._r8) then
                W1TO4 = 0._r8
                W1TO6 = 0._r8
                W2TO4 = FR4 / FR2 
                W2TO6 = 1._r8 - W2TO4
                W3TO4 = 0._r8
                W3TO6 = 0._r8
              else if (FR2.eq.0._r8 .and. FR3.eq.0._r8) then
                W1TO4 = FR4 / FR1 
                W1TO6 = 1._r8 - W1TO4
                W2TO4 = 0._r8 
                W2TO6 = 0._r8
                W3TO4 = 0._r8
                W3TO6 = 0._r8
              else if (FR4 .le. FR1) then
                W1TO4 = FR4 / FR1
                W1TO6 = 1._r8 - W1TO4
                W2TO4 = 0._r8
                W2TO6 = 1._r8
                W3TO4 = 0._r8
                W3TO6 = 1._r8
              else if (FR4 .le. FR1+FR2) then
                W1TO4 = 1._r8
                W1TO6 = 0._r8
                W2TO4 = (FR4 - FR1) / FR2
                W2TO6 = 1._r8 - W2TO4
                W3TO4 = 0._r8
                W3TO6 = 1._r8
              else
                W1TO4 = 1._r8
                W1TO6 = 0._r8
                W2TO4 = 1._r8
                W2TO6 = 0._r8
                W3TO6 = FR6 / FR3
                W3TO4 = 1._r8 - W3TO6
              end if
            else if ( (FR1 + FR2 + FR3) .eq. CFXX) then
              FR4 = FR1 + FR2 + FR3
              FR5 = 0._r8 !aviod CFXX-(FR1+FR2+FR3) when they are equal
              FR6 = 0._r8
              W1TO4 = 1._r8
              W1TO6 = 0._r8
              W2TO4 = 1._r8
              W2TO6 = 0._r8
              W3TO4 = 1._r8
              W3TO6 = 0._r8
            else
              FR4 = FR1 + FR2 + FR3
              FR5 = CFXX - (FR1 + FR2 + FR3)
              FR6 = 0._r8
              W1TO4 = 1._r8
              W1TO6 = 0._r8
              W2TO4 = 1._r8
              W2TO6 = 0._r8
              W3TO4 = 1._r8
              W3TO6 = 0._r8
            end if
            PR4 = W1TO4*PR1 + W2TO4*PR2 + W3TO4*PR3
            PR6 = W1TO6*PR1 + W2TO6*PR2 + W3TO6*PR3
            if (PR4 .gt. 0._r8) then
              D4 = (W1TO4*PR1*D1 + W2TO4*PR2*D2 + W3TO4*PR3*D3) / PR4
            else
              D4 = 0._r8
            end if
            if (PR6 .gt. 0._r8) then
              D6 = (W1TO6*PR1*D1 + W2TO6*PR2*D2 + W3TO6*PR3*D3) / PR6
            else
              D6 = 0._r8
            end if
          else
            FR4 = 0._r8
            FR5 = CFXX
            FR6 = 0._r8
            D4  = 0._r8
            D6  = 0._r8
            PR4 = 0._r8
            PR6 = 0._r8
            W1TO4 = 0._r8
            W1TO6 = 0._r8
            W2TO4 = 0._r8
            W2TO6 = 0._r8
            W3TO4 = 0._r8
            W3TO6 = 0._r8
          end if

          do N = 1, NTM
            QTTOP4(N) = W1TO4*QTTOP1(N) + W2TO4*QTTOP2(N) + W3TO4*QTTOP3(N)
            QTTOP6(N) = W1TO6*QTTOP1(N) + W2TO6*QTTOP2(N) + W3TO6*QTTOP3(N)
          end do

          !// Get TOP values
          PR4T = PR4
          D4T  = D4
          PR6T = PR6
          D6T  = D6
          FR6T = FR6
          if (CFXX .gt. 0._r8) then
            CLWX = CLWC(L) + CIWC(L)
          else
            CLWX = 0._r8
          end if

!-----Evaporate ambient precip and decrease area-------------------------
!-----If ice, diam=diam falling from above  If rain, diam=4mm (not used)
!-----Evaporate tracer contained in evaporated precip
!-----Can't evaporate more than we start with-----------------------------
!-----Don't do washout until we adjust ambient precip to match Rbot if needed
!------(after RNEW if statements)

          FR6 = max(0._r8, FR6T - CEVAP6(L)*FR6T)
          !// Do not use: PR6 = max(0._r8,PR6T*(1._r8-CEVAP6(L))**2)   !kg/s
          PR6 = max(0._r8, PR6T - CEVAP6(L)*PR6T)   !kg/s
          if (FR6 .gt. 0._r8) then
            if (TEM(L) .lt. TICE) then
              D6 = D6T      !mm
            else
              D6 = 4.0_r8    !mm - not necessary
            end if
          else
            D6 = 0._r8
          end if
          if (PR6T .gt. 0._r8) then
            do N = 1, NTM
              QTEVP6X(N) = (PR6T-PR6)/PR6T*QTTOP6(N)
            end do
          end if

          RNEW = RLS(L) - PR4T - PR6     !kg/s

          !// if RNEW>0, there is growth and/or new precip formation
          !if(RNEW.gt.0._r8 .and. CFXX.gt.0._r8) then
          if (RNEW .gt. 0._r8) then

            !// Min cloudwater requirement for cloud with new precip
            !// CWMIN is only needed for new precip formation - do not
            !// need for RNEW<0
            !// CLWX = max((CLWC(L)+CIWC(L)),CWMIN)

!-------------------------ICE-----------------------------------------
!-----For ice and mixed phase, grow precip in old cloud by riming
!-----Use only portion of cloudwater in old cloud fraction
!-----and rain above old cloud fraction
!-----COLEFF from Lohmann and Roeckner (1996), Loss rate from Rotstayn (1997)
            if (TEM(L) .lt. TICE) then
              COLEFFSNOW = exp(2.5e-2_r8 * (TEM(L) - 273._r8))

              if(TEM(L) .le. TFROZ) then
                RHOSNOW = RHOSNOWFIX
              else
                RHOSNOW = 0.303_r8 * (TEM(L) - TFROZ) * RHOSNOWFIX
              end if
              if(FR4 .gt. 0._r8) then
                if(D4T .gt. 0._r8) then
                  DELTARIMEMASS = &
                       (CLWX * QM(L)) * (FR4 / CFXX) &
                       * (1._r8 - exp((-COLEFFSNOW / (D4T * 1.e-3_r8)) &
                          * (PR4T / (2._r8 * RHOSNOW * GAREA * FR4)) * DTSCAV))
                else
                  DELTARIMEMASS = 0._r8
                end if
              else
                DELTARIMEMASS = 0._r8
              end if
              !// Increase in precip rate due to riming (kg/s):
              !// Limit to total increase in R in cloud
              if (FR4 .gt. 0._r8) then
                DELTARIME = min(RNEW, DELTARIMEMASS / DTSCAV)
              else
                DELTARIME = 0._r8
              end if
              !// Find diameter of rimed precip, must be at least .1mm
              if (PR4T .gt. 0._r8) then
                DOR = max(DMIN,(((PR4T + DELTARIME) / PR4T)**VOLPOW) * D4T)
              else
                DOR = 0._r8
              end if
!-----If there is some in-cloud precip left, we have new precip formation
!-----Will be spread over whole cloud fraction 
!------Calculate precip rate in old and new cloud fractions
              RPRECIP = RNEW - DELTARIME                    !kg/s
              !// Calculate precip rate in old and new cloud fractions
              PR4 = PR4T + DELTARIME + RPRECIP * FR4 / CFXX !kg/s
              PR5 = RPRECIP * FR5 / CFXX                    !kg/s

!-----Find diameter of new precip from empirical relation using Rprecip
!----in given area of box- use density of water, not snow, to convert kg/s
!----to mm/s -> as given in Field and Heymsfield, J. Atm. Sci., 60, 544-560,
!----           2003, doi:10.1175/1520-0469(2003)060%3C0544:AASOIC%3E2.0.CO;2

!-----Also calculate diameter of mixed precip,D4, from empirical relation
!-----using total R in FR4 - this will give larger particles than averaging DOR and
!-----DNEW in the next level
!-----DNEW and D4 must be at least .1mm
              if (RPRECIP .gt. 0._r8) then
                !// New rain is formed
                WEMP = CLWX * QM(L) / (GAREA * CFXX * DELZ(L))
                !mm/s local
                REMP = RPRECIP / (GAREA * CFXX * RHORAIN * 1.e-3_r8)

                call DIAMEMP(WEMP,REMP,DNEW)

                DNEW = max(DMIN, DNEW)

                !// in sub-box 5
                if (FR5 .gt. 0._r8) then
                  D5 = DNEW
                else
                  D5 = 0._r8
                end if
              else
                D5 = 0._r8
              end if

              !// in sub-box 4
              if (FR4 .gt. 0._r8) then
                WEMP = CLWX * QM(L) / (GAREA * CFXX * DELZ(L)) !kg/m3
                !mm/s local
                REMP = PR4 / (GAREA * FR4 * RHORAIN * 1.e-3_r8)

                call DIAMEMP(WEMP,REMP,DEMP)

                D4 = ((PR4T + DELTARIME) / PR4) * DOR &
                       + (RPRECIP * FR4 / CFXX / PR4) * DNEW
                D4 = max(DEMP, D4)
                D4 = max(DMIN, D4)
              else
                D4 = 0._r8
              end if


              !// ICE SCAVENGING
              !// --------------------------------------------------------
!------For ice, rainout only hno3/aerosols using new precip
!------Tracer dissolved given by Kaercher and Voigt (2006) for T<258K
!-----For T>258K, use Henry's Law with Retention coefficient
!------Rain out in whole CF
              if (RPRECIP .gt. 0._r8) then
                do N = 1, NTM
                  if (QTT(L,N) .gt. QTTMIN) then
                    if (LICETYP(N) .eq. 1) then
                      RRAIN = RPRECIP !kg/s

                      !// Calculate effective Henry expression
                      call getHstar(chem_idx(N), TCHENA(N), &
                           TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                      !// Calculate dissolved tracer
                      call DISGAS2(chem_idx(N),CLWX, CFXX, TMASS(N), HSTAR, &
                           TEM(L), POFL(L), QM(L), &
                           QTT(L,N)*CFXX, SCV_RETEFF(N), SCV_T258(N), &
                           SCViceFR(N), QTDISCF)

                      !// Rain out
                      call RAINGAS(RRAIN, DTSCAV, CLWX, CFXX, &
                           QM(L), QTT(L,N), QTDISCF, QTPRC0(N))

                      !// Prec.rate in old cloud
                      QTPRC4(N) = FR4 * QTPRC0(N) / CFXX
                      !// Prec.rate of new prec.
                      QTPRC5(N) = FR5 * QTPRC0(N) / CFXX
                    else if(LICETYP(N) .eq. 2) then
                      QTPRC4(N) = 0._r8 !// No tracer in old ice prec.
                      QTPRC5(N) = 0._r8 !// No tracer in new ice prec.
                    end if
                  end if
                end do
              end if !// if (RPRECIP .gt. 0._r8) then

              !// For ice, accretion removal for hno3 and aerosols is 
              !// propotional to riming, no accretion removal for gases
              !// remove only in mixed portion of cloud
              !// Limit DELTARIMEMASS to RNEW*DTSCAV for ice - evaporation
              !// of rimed ice to match
              !// RNEW precip rate would result in HNO3 escaping from
              !// ice (no trapping) 
              if (DELTARIME .gt. 0._r8) then
                do N = 1, NTM
                  if (QTT(L,N) .gt. QTTMIN) then
                    if (LICETYP(N) .eq. 1) then

!!!! this looks odd below, discontinuous at TFROZ
                      if(TEM(L) .le. TFROZ) then
                        RHOSNOW = RHOSNOWFIX
                      else
                        RHOSNOW = 0.303_r8 * (TEM(L) - TFROZ) * RHOSNOWFIX
                      end if

                      !// Calculate effective Henry expression
                      call getHstar(chem_idx(N), TCHENA(N), &
                           TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                      !// This process is riming, but since this is ice,
                      !// we apply SCViceFR, so that not all tracer in
                      !// cloud may rime.
                      call DISGAS2(chem_idx(N), CLWX * FR4 / CFXX, FR4, TMASS(N), &
                          HSTAR, TEM(L), POFL(L), &
                          QM(L), QTT(L,N)*FR4 - QTPRC4(N), &
                          SCV_RETEFF(N), SCV_T258(N), SCViceFR(N), &
                          QTDISRIME)

                      QTDISSTAR = QTDISRIME

                      QTRIM4(N) = (QTT(L,N)*FR4 - QTPRC4(N)) &
                             * (1._r8 - exp(-COLEFFSNOW / (D4T * 1.e-3_r8) &
                             * PR4T / (2._r8 * RHOSNOW * GAREA * FR4) &
                             * QTDISSTAR / (QTT(L,N)*FR4 - QTPRC4(N))*DTSCAV))
                      QTRIM4(N) = min(QTRIM4(N), &
                                      DELTARIME * DTSCAV / &
                                      (CLWX * QM(L) * FR4 / CFXX) * QTDISSTAR)

                    else if (LICETYP(N) .eq. 2) then
                      QTRIM4(N) = 0._r8 !// No tracer to be removed by ice
                    end if
                  end if
                end do
              end if !// if (DELTARIME .gt. 0._r8) then
              !// For ice, no washout in interstitial cloud air
              QTWSH4(:) = 0._r8     !! not necessary?
              QTEVP4(:) = 0._r8     !! not necessary?


           else !// if(TEM(L) .lt. TICE) then

             !// RAIN
             !// --------------------------------------------------------
             !// For rain, accretion increases rain rate but diameter
             !// remains constant.
             !// Diameter is 4mm (not used)
             if (FR4 .gt. 0._r8) then
               DELTARIMEMASS = &
                    CLWX * QM(L) * FR4 / CFXX &
                    * (1._r8 - exp(-0.24_r8 * COLEFFRAIN &
                            * ((PR4T / (GAREA * FR4))**0.75_r8) * DTSCAV))
             else
               DELTARIMEMASS = 0._r8
             end if
             !// Increase in precip rate due to riming (kg/s):
             !// Limit to total increase in R in cloud
             if (FR4 .gt. 0._r8) then
               DELTARIME = min(RNEW, DELTARIMEMASS / DTSCAV)
             else
               DELTARIME = 0._r8
             end if
             !// If there is some in-cloud precip left, we have new
             !// precip formation 
             RPRECIP = RNEW - DELTARIME

             PR4 = PR4T + DELTARIME + RPRECIP * FR4 / CFXX  !kg/s
             PR5 = RPRECIP * FR5 / CFXX                     !kg/s

             if (FR4 .gt. 0._r8) then
               D4 = 4._r8
             else
               D4 = 0._r8
             end if
             if (FR5 .gt. 0._r8) then
               D5 = 4._r8
             else
               D5 = 0._r8
             end if

             !// RAIN SCAVENGING
             !// --------------------------------------------------------
             !// For rain, rainout both hno3/aerosols and gases using
             !// new precip
             if (RPRECIP .gt. 0._r8) then
               !// RPRECIP is new precipitation (RNEW-DELTARIME)
               do N = 1, NTM
                 if (LWSHTYP(N) .eq. 3) then
                   !// No new precip for type 3
                   QTPRC4(N) = 0._r8
                   QTPRC5(N) = 0._r8
                 else if (TCHENA(N) .gt. 0._r8 .and. &
                      QTT(L,N).gt.QTTMIN) then
                   !// There is tracer which can be dissolved

                   RRAIN = RPRECIP           !kg/s

                   !// Calculate effective Henry expression
                   call getHstar(chem_idx(N), TCHENA(N), &
                        TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                   !// T>273: Instead of SCViceFR, send in 1._r8
                   call DISGAS2(chem_idx(N), CLWX, CFXX, TMASS(N), HSTAR, &
                        TEM(L), POFL(L), QM(L), &
                        QTT(L,N)*CFXX, SCV_RETEFF(N), SCV_T258(N), &
                        1._r8, QTDISCF)

                   call RAINGAS(RRAIN, DTSCAV, CLWX, CFXX, &
                        QM(L), QTT(L,N), QTDISCF, QTPRC0(N))

                   QTPRC4(N) = FR4 * QTPRC0(N) / CFXX
                   QTPRC5(N) = FR5 * QTPRC0(N) / CFXX
                 end if
               end do
             end if !// if(RPRECIP .gt. 0._r8) then

!--------For rain, accretion removal is propotional to riming
!--------caclulate for hno3/aerosols and gases
!-------Remove only in mixed portion of cloud
!-------Limit DELTARIMEMASS to RNEW*DTSCAV
             if (DELTARIME .gt. 0._r8) then

               do N = 1, NTM
                 if (LWSHTYP(N).eq.3) then
                   !// No riming for type 3
                   QTRIM4(N) = 0._r8
                 else if (TCHENA(N) .gt. 0._r8 .and. &
                      QTT(L,N) .gt. QTTMIN) then

                   !// Calculate effective Henry expression
                   call getHstar(chem_idx(N), TCHENA(N), &
                        TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                   !// There is tracer which can be dissolved
                   !// T>273: Instead of SCViceFR, send in 1._r8
                   call DISGAS2(chem_idx(N), CLWX * FR4 / CFXX, FR4, TMASS(N), &
                        HSTAR, TEM(L), POFL(L), &
                        QM(L), QTT(L,N)*FR4 - QTPRC4(N), &
                        SCV_RETEFF(N), SCV_T258(N), 1._r8, QTDISRIME)

                   QTDISSTAR = QTDISRIME

                   QTRIM4(N) = (QTT(L,N)*FR4 - QTPRC4(N)) &
                        * (1._r8 - exp(-0.24_r8 * COLEFFRAIN &
                                 * (PR4T / (GAREA*FR4))**0.75_r8 &
                                 * QTDISSTAR &
                                   / (QTT(L,N)*FR4 - QTPRC4(N)) * DTSCAV))

                   QTRIM4(N) = min(QTRIM4(N), &
                                   DELTARIME * DTSCAV / &
                                     (CLWX * QM(L) * FR4 / CFXX) * QTDISSTAR)
                 else
                   QTRIM4(N)=0._r8
                 end if
               end do
             end if !// if(DELTARIME .gt. 0._r8) then


             !// For rain, washout gases and HNO3/aerosols using rain from
             !// above old cloud
             !// Washout for HNO3/aerosols is only on non-dissolved
             !// portion, impaction-style
             !// Washout for gases is on non-dissolved portion, limited by
             !// QTTOP+QTRIME
             if (PR4T .gt. 0._r8) then

               do N = 1, NTM
                 !// Available for possible washout
                 QTPRECIP = FR4*QTT(L,N) - QTPRC4(N) - QTRIM4(N)

                 if(LWSHTYP(N).eq.1 .and. QTT(L,N).gt.QTTMIN) then
                   !// Impact-style
                   if(QTPRECIP .gt. 0._r8) then
                     QTWSH4(N) = QTPRECIP &
                          * (1._r8 - exp(-0.24_r8 * COLEFFAER &
                          * ((PR4T / (GAREA*FR4))**0.75_r8) * DTSCAV)) !local
                   else
                     QTWSH4(N) = 0._r8
                   end if
                   QTEVP4(N) = 0._r8
                 else if (LWSHTYP(N).eq.2 .and. QTT(L,N).gt.QTTMIN) then
                   RWASH = PR4T                    !kg/s local
                   if (QTPRECIP .gt. 0._r8) then

                     !// Calculate effective Henry expression
                     call getHstar(chem_idx(N), TCHENA(N), &
                           TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                     !// CTM3: Send in SCV_RETEFF(N)
                     call WASHGAS2(chem_idx(N), RWASH, FR4, DTSCAV, &
                          QTTOP4(N)+QTRIM4(N), &
                          HSTAR, TEM(L), &
                          POFL(L), QM(L), QTPRECIP, SCV_RETEFF(N), &
                          QTWSH4(N), QTEVP4(N))
                   else
                     QTWSH4(N) = 0._r8
                     QTEVP4(N) = 0._r8
                   end if
                 else
                   !// Other LWSHTYP where rain washout is not to be done
                   QTWSH4(N) = 0._r8
                   QTEVP4(N) = 0._r8
                 end if
               end do
             end if !// if (PR4T .gt. 0._r8) then

           end if !// if(TEM(L) .lt. TICE) then

           !// -------------------------------------------------------------
           !// End of RNEW > 0
           !          elseif (RNEW.gt.0._r8 .and. CFXX.eq.0._r8) then
           !            PR6 = RLS(L)
           !            do N = 1,NTM
           !              QTEVP6X(N) = QTTOP6(N)*(PR6T-PR6)/PR6T
           !            enddo
           !// -------------------------------------------------------------
         else               !// if (RNEW .gt. 0._r8) then
           !// -------------------------------------------------------------
           !// If RNEW<O, confine precip to area of cloud above-----------
           !// FR4 does not require a minimum (could be zero if R(L).le.what
           !// evaporated in ambient)
           !//
           !// Put rain into cloud up to PR8 so that we evaporate
           !// from ambient first
           !// Adjust ambient to try to match RLS(L)
           !// If no cloud, PR9=R(L)
           if (PR6T .gt. 0._r8) then
             if (RLS(L) .le. PR4T) then
               !// Less in metdata than comes in on top; adjust PR4, set PR6=0
               PR6 = 0._r8
               FR6 = 0._r8
               PR4 = RLS(L)
             else
               !// Some rain falls into ambient
               PR6 = RLS(L)-PR4T
               PR4 = PR4T
               !// Do not use: FR6 = sqrt(PR6/PR6T)*FR6T
               FR6 = PR6 / PR6T * FR6T  ! shrink FR6 at the same rate as PR6?
               do N = 1, NTM
                 QTEVP6X(N) = (PR6T-PR6)/PR6T*QTTOP6(N)
               end do
             end if
           else
             !// No rain through ambient (clear) air
             PR6 = 0._r8
             FR6 = 0._r8
             PR4 = RLS(L) !// LS rain from metdata
           end if

           !// IN-CLOUD EVAPORATION/WASHOUT
           !// If precip out cloud bottom is 0, evaporate everything
           !// If there is no cloud, QTTOP1(N)=0, so nothing happens
           if (PR4.le.0.0 .or. PR4T.le.0._r8) then
             QTEVP4(:) = QTTOP4(:)
             PR4 = 0._r8
             D4 = 0._r8
           else
             !// If rain out the bottom of the cloud is >0 (but .le. PR8):
             !// For ice, decrease particle size,
             !// no washout
             !// no evap for non-ice gases (b/c there is nothing in ice)
             !// T<Tmix,release hno3 & aerosols:
             !//    release is amount dissolved in ice mass released
             !// T>Tmix, hno3&aerosols are incorporated into ice structure:
             !//    do not release
             !// For rain, assume full evaporation of some raindrops
             !//    proportional evaporation for all species 
             !//    washout for gases using Rbot 
             !//    impact washout for hno3/aerosol portion in gas phase

             if (TEM(L) .lt. TICE) then
               D4 = ((PR4 / PR4T)**VOLPOW) * D4T
               do N = 1, NTM
                 QTWSH4(N) = 0._r8
                 if (LICETYP(N).eq.1 .and. QTT(L,N).gt.QTTMIN) then
                   if (TEM(L) .lt. TMIX) then
                     !// T < 258K
                     MASSLOSS = (PR4T - PR4) * DTSCAV

                     !// note-QTT doesn't matter b/c T<258K
                     !// CTM3: Only KV06 needs to be done; Henrys law
                     !// is not necessary.
                     !// I wonder if there is something wrong/weird with
                     !// this: should we really send in QTT?
                     if (SCV_T258(N) .eq. 1) then

                       !// Calculate effective Henry expression
                       call getHstar(chem_idx(N), TCHENA(N), &
                             TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                       call DISGAS2(chem_idx(N), MASSLOSS/QM(L), FR4, TMASS(N), &
                            HSTAR, TEM(L), POFL(L), &
                            QM(L), QTT(L,N), SCV_RETEFF(N), SCV_T258(N), &
                            SCViceFR(N), QTEVP4(N))

!not sure how to fix this call above, QTEVP4(N) is now limited by QTT(L)
!as in std ratio.
                       QTEVP4(N) = 1._r8 / (1._r8/QTEVP4(N) - 1._r8/QTT(L,N))
                       QTEVP4(N) = min(QTTOP4(N), QTEVP4(N))
                     else if (SCV_T258(N) .eq. 2) then

                       !// Assume solved in ice for T<258K
                       !// Calculate effective Henry expression
                       call getHstar(chem_idx(N), TCHENA(N), &
                             TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                       !// Routine needs QTT in sub-box as input!
                       call DISGAS2(chem_idx(N), MASSLOSS/QM(L), FR4, TMASS(N), &
                         HSTAR, TEM(L), POFL(L), &
                         QM(L), QTT(L,N)*FR4, SCV_RETEFF(N), SCV_T258(N), &
                         SCViceFR(N), QTEVP4(N))
                       QTEVP4(N) = 1._r8 / (1._r8/QTEVP4(N) - 1._r8/QTT(L,N))
                       QTEVP4(N) = min(QTTOP4(N),QTEVP4(N))

                     else
                       !// Aerosols etc.: May include other treatments here
                       QTEVP4(N) = 0._r8
                     end if
                   else
                     !// TMIX < T
                     QTEVP4(N) = 0._r8
                   end if
                 else if(LICETYP(N).eq.2 .and. QTT(L,N).gt.QTTMIN) then
                   QTEVP4(N) = 0._r8
                 end if
               end do
             else      !// if (TEM(L) .lt. TICE) then
               !// T > 273K
               D4 = 4._r8
               do N = 1, NTM
                 !// QTEVP4X is first part of evaporation, the part that
                 !// needs to evaporate due to the change (PR4T-PR4)
                 QTEVP4X(N) = (PR4T - PR4) / PR4T * QTTOP4(N)
                 QTCXA = FR4 * QTT(L,N)

                 if (LWSHTYP(N).eq.1 .and. QTT(L,N).gt.QTTMIN) then
                   !// Mass limited evaporation (e.g. HNO3/aerosols)

                   !// Calculate effective Henry expression
                   call getHstar(chem_idx(N), TCHENA(N), &
                         TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                   call DISGAS2(chem_idx(N), CLWX*(FR4/CFXX), FR4, TMASS(N), &
                        HSTAR, TEM(L), POFL(L), &
                        QM(L), QTCXA, SCV_RETEFF(N), SCV_T258(N), &
                        1._r8, QTDISCXA)

!>>>>>>>not sure how to fix this call above, QTEVP4(N) is now limited by QTT(L)
!>>>>>>>as in std ratio.
                   QTDISCXA = 1._r8 / (1._r8/QTDISCXA - 1._r8/QTCXA)

                   if ((QTCXA - QTDISCXA) .gt. 0._r8) then
                     !// Need to wash out more
                     QTWSH4(N) = (QTCXA - QTDISCXA) &
                          * (1._r8 - exp(-0.24_r8 * COLEFFAER &
                          * ((PR4/(GAREA*FR4))**0.75_r8)*DTSCAV)) !local

                   else
                     QTWSH4(N) = 0._r8
                   end if
                   !// Do not allow evaporation of tracer
                   QTEVP4Y(N) = 0._r8
                 else if (LWSHTYP(N).eq.2 .and. QTT(L,N).gt.QTTMIN) then
                   RWASH = PR4

                   !// Calculate effective Henry expression
                   call getHstar(chem_idx(N), TCHENA(N), &
                         TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                   !// CTM3: Send in SCV_RETEFF(N)
                   call WASHGAS2(chem_idx(N), RWASH, FR4, DTSCAV, QTTOP4(N), HSTAR, &
                        TEM(L), POFL(L), QM(L), &
                        !//Old: QTCXA-QTDISCXA,QTWSH4(N),QTEVP4Y(N))  !???
                        QTCXA, SCV_RETEFF(N), QTWSH4(N), QTEVP4Y(N))
                 else if (LWSHTYP(N).eq.3 .and. QTT(L,N).gt.QTTMIN) then
                   !// Allow evaporation 4X, but not 4Y. There is no QTWSH4
                   !// in rain for this case, so there is only old rain to
                   !// evaporate.
                   QTEVP4Y(N)= 0._r8
                   QTWSH4(N) = 0._r8
                 end if
                 QTEVP4(N) = QTEVP4X(N) + QTEVP4Y(N)
               end do
             end if
           end if !// if (PR4.le.0.0 .or. PR4T.le.0._r8) then

         end if      !// End RNEW if statements


         !// AMBIENT WASHOUT
         !// Ambient precip is finalized - if it is rain, washout
         !// no ambient washout for ice, since gases are in vapor phase
         if (PR6 .gt. 0._r8) then
           if (TEM(L) .ge. TICE) then
             do N = 1, NTM
               QTAX = FR6 * QTT(L,N)
               if (LWSHTYP(N).eq.1 .and. QTT(L,N).gt.QTTMIN) then
                 QTWSH6(N) = QTAX * (1._r8 - exp(-0.24_r8 * COLEFFAER &
                      * ((PR6/(GAREA*FR6))**0.75_r8)*DTSCAV))  !local

                 QTEVP6Y(N) = 0._r8
               else if (LWSHTYP(N).eq.2 .and. QTT(L,N).gt.QTTMIN) then
                 RWASH = PR6  !kg/s local

                 !// Calculate effective Henry expression
                 call getHstar(chem_idx(N), TCHENA(N), &
                      TCHENB(N), TCKAQA(N), TEM(L), HCALLER, HSTAR)

                 call WASHGAS2(chem_idx(N), RWASH, FR6, DTSCAV, QTTOP6(N), HSTAR, &
                     TEM(L), POFL(L), QM(L), QTAX, &
                     SCV_RETEFF(N), QTWSH6(N), QTEVP6Y(N))
               else if (LWSHTYP(N).eq.3) then
                 !// No ambient wash out
                 QTEVP6Y(:)=0._r8
                 QTWSH6(:)=0._r8
               end if
             end do
           else
             QTEVP6Y(:)=0._r8
             QTWSH6(:)=0._r8
           end if
         else
           QTEVP6Y(:)=0._r8
           QTWSH6(:)=0._r8
         end if
         do N = 1, NTM
           QTEVP6(N) = QTEVP6X(N) + QTEVP6Y(N)
         end do

         !// Store PRs and QTs for next layer
         PR7 = PR4                 !kg/s
         PR8 = PR5
         PR9 = PR6

         FR7 = FR4
         FR8 = FR5
         FR9 = FR6

         D7  = D4
         D8  = D5
         D9  = D6

         !// Net loss can not exceed QTT in each region

         do N = 1, NTM
           if (LWSHTYP(N) .gt. 0) then     !  Solubility check
             QTNET4(N) = QTPRC4(N) + QTRIM4(N) + QTWSH4(N) - QTEVP4(N)
             QTNET4(N) = min(QTT(L,N)*FR4,QTNET4(N))

             QTNET5(N) = QTPRC5(N)
             QTNET5(N) = min(QTT(L,N)*FR5,QTNET5(N))

             QTNET6(N) = QTWSH6(N) - QTEVP6(N)
             QTNET6(N) = min(QTT(L,N)*FR6,QTNET6(N))
           
             QTTNEW(L,N) = QTT(L,N) - QTNET4(N) - QTNET5(N) - QTNET6(N)
             !// Make sure small noise is reduced.
             if (QTTNEW(L,N) .lt. 0._r8 .and. QTTNEW(L,N) .gt. -1.e-10_r8) &
                  QTTNEW(L,N) = 1.e-20_r8

             !// Tracer masses to send down to next level: rain above
             QTBOT7(N) = QTTOP4(N) + QTNET4(N)
             !// Tracer masses to send down to next level: no rain above, only new precip
             QTBOT8(N) = QTNET5(N)
             !// Tracer masses to send down to next level: ambient
             QTBOT9(N) = QTTOP6(N) + QTNET6(N)
           end if
         end do

       end if   !End RLS(L)<=0._r8 or QTT(L,N)<=0._r8 statements

       !// END SCAVENGING
       !// -----------------------------------------------------------------

       do N = 1, NTM
          if ( (QTTNEW(L,N).lt.0._r8 .and. QTT(L,N).gt.0._r8) &
               .or. (QTTNEW(L,N).ne.QTTNEW(L,N)) ) then
            write(6,'(a)') f90file//':'//subr// &
                 ': QTTNEW is negative or NAN!!!!!!'
            print*,'N=',N,'L=',L,'I=',I,'J=',J,'M=',MP
            print*,'QTT=',QTT(L,N)
            print*,'QTTNEW=',QTTNEW(L,N)
            print*,'RLS=',RLS(L),'DRLS=',DRLS(L)
            print*,'RNEW=',RNEW/GAREA
            print*,'PR4=',PR4,'PR8(4X)=',PR4T
            print*,'PR5=',PR5,'PR6(6X)=',PR6T
            print*,'PR9(6)=',PR6
            print*,'RPRECIP=',RPRECIP
            print*,'DELTARIME=',DELTARIME
            print*,'RRAIN=',RRAIN
            print*,'CFXX=',CFXX
            print*,'FR4=',FR4,'FR8(4X)=',FR4
            print*,'FR5=',FR5,'FR6(6X)=',FR6T,'FR9(6)=',FR6
            print*,'QTPRC4=',QTPRC4(N),'QTRIM4=',QTRIM4(N)
            print*,'QTWSH4=',QTWSH4(N),'QTEVP4=',QTEVP4(N)
            print*,'QTPRC5=',QTPRC5(N),'QTWSH6=',QTWSH6(N)
            print*,'QTEVP6=',QTEVP6(N)
            print*,'QT123=',QTTOP1(N),QTTOP2(N),QTTOP3(N)
            print*,'QTTOP4=',QTTOP4(N),'QTTOP6=',QTTOP6(N)
            print*,'SUM_PR=',PR8+PR9+PR7
            print*,'PR7=',PR7,'PR8=',PR8,'PR9=',PR9
            print*,'SUM_QT=',QTBOT8(N)+QTBOT9(N)+QTBOT7(N)
            print*,'FR456=',FR4,FR5,FR6
            print*,'FR123=',FR1,FR2,FR3
            print*,'SUM456=',FR4+FR5+FR6
            print*,'SUM123=',FR1+FR2+FR3
            print*,'W4',W1TO4,W2TO4,W3TO4
            print*,'W6',W1TO6,W2TO6,W3TO6
            print*,'!!!!!!!'
            if(QTTNEW(L,N).ne.QTTNEW(L,N)) print*,' WASH2 NAN'
            stop 'stop in scavenging_largescale_uci.f90: WASH2'
          end if

        end do

      end do !// L-loop

      !// reload new tracer mass and down-scale moments if tracer lost: 
      do L = 1, LTOP
        do N = 1, NTM
          if (LWSHTYP(N) .gt. 0) then     !  Solubility check
            BTT(L,N,II,JJ) = QTTNEW(L,N) + QTTnoSCAV(L,N)
 
            if (QTTNEW(L,N) .lt. QTT(L,N)) then
              !FXYW1 = QTTNEW(L,N)/QTT(L,N)
              FXYW1 = (QTTNEW(L,N) + QTTnoSCAV(L,N))/ &
                                (QTT(L,N) + QTTnoSCAV(L,N))
              BZT(L,N,II,JJ) = BZT(L,N,II,JJ) * FXYW1
              BZZ(L,N,II,JJ) = BZZ(L,N,II,JJ) * FXYW1
              BXZ(L,N,II,JJ) = BXZ(L,N,II,JJ) * FXYW1
              BYZ(L,N,II,JJ) = BYZ(L,N,II,JJ) * FXYW1
              BXT(L,N,II,JJ) = BXT(L,N,II,JJ) * FXYW1
              BXX(L,N,II,JJ) = BXX(L,N,II,JJ) * FXYW1
              BYT(L,N,II,JJ) = BYT(L,N,II,JJ) * FXYW1
              BYY(L,N,II,JJ) = BYY(L,N,II,JJ) * FXYW1
              BXY(L,N,II,JJ) = BXY(L,N,II,JJ) * FXYW1
            end if
          end if
        end do
      end do

      !// ------------------------------------------------------------------
    end do !// I block
    end do !// J block
    !// --------------------------------------------------------------------

    !// --------------------------------------------------------------------
  end subroutine WASH2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine KV06 (CLWX,CFX,MOLMASS,TM,QM,QT,QTDIS)
    !// --------------------------------------------------------------------
    !// For T below solid phase, only calculate buried HNO3
    !// (MU = mole fraction) from empirical fit to
    !// Kaercher and Voigt (2006)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: CLWX,CFX      !cloud water,cloud fraction 
    real(r8), intent(in) :: MOLMASS     !molecular mass of tracer
    real(r8), intent(in) :: TM          !temperature of box (K)
    real(r8), intent(in) :: QM          !air mass in box (kg)
    real(r8), intent(in) :: QT          !tracer in sub-box (kg)
    real(r8), intent(out) :: QTDIS      !tracer dissolved in aqueous phase 
 
    real(r8)  MUEMP,TDIS1
    real(r8)  TMASS,AMASS,WMASS
    !real(r8), parameter  :: TSOLID = 258._r8 ! Set in module
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'KV06'
    !//---------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Do not use KV06 - STOP'
    stop

    AMASS = QM*CFX       ! kg of air in sub-box
    WMASS = CLWX*QM      ! kg of liquid water in sub-box
    TMASS = QT           ! kg of tracer in sub-box

    !// Treatment for T below solid phase
    if (TM .lt. TSOLID) then
       !// For T below solid phase, only calculate buried HNO3
       !// (MU = mole fraction) from empirical fit to
       !// Kaercher and Voigt (2006)
       MUEMP = exp(-14.2252_r8 + 0.155704_r8*TM - 7.1929e-4_r8*(TM**2))
       TDIS1 = MUEMP*(MOLMASS/18._r8)*WMASS
       QTDIS = TDIS1*TMASS/(TDIS1 + TMASS)
    else
       QTDIS = 0._r8
    end if

    !// --------------------------------------------------------------------
  end subroutine KV06
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getHstar(chem_idx,KHA,KHB,HFLAG,TM,HCALLER,HSTAR)
    !// --------------------------------------------------------------------
    !// From Henry constants and the flag HFLAG, the effective HSTAR is
    !// calculated.
    !//
    !// KHA & KHB are Henrys Law coefficients (not extrapolated below 273K)
    !// Skipping the use of UCI KAQ, which was a pre-calculated (vs. pH)
    !// enhancement due to dissociation/polymers.
    !// Rather use HFLAG to set these here, or enhance KHA directly in
    !// input table.
    !//
    !// Effective Henry's Law constant:
    !/           H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)   :: chem_idx   !Tracer number
    real(r8), intent(in)  :: KHA,KHB    !Henry's Law coeffs A*exp(-B/T)
    real(r8), intent(in)  :: TM         !Temperature (K)
    real(r8), intent(in)  :: HFLAG      !Flag to denote hard-coded Henry
    integer, intent(in)   :: HCALLER    !Flag for calling routine 1:LS,2:CNV
    real(r8), intent(out) :: HSTAR      !HSTAR
 
    real(r8), parameter  :: INV298 = 1._r8/298._r8
    real(r8), parameter  :: HPLUS  = 3.16e-5_r8
    real(r8) :: HTEMP, TEMPFAC, HTMP1, HTMP2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'getHstar'
    !//---------------------------------------------------------------------

    if (HCALLER .eq. 1) then
       !// Called from LS scavenging
       HTEMP = max(TM, TICE)
    else if (HCALLER .eq. 2) then
       !// Called from CNV scavenging
       HTEMP = TM
    else
       write(6,'(a)') f90file//':'//subr// &
            ': HCALLER unknown - need to set HTEMP'
       stop
    end if

    TEMPFAC = 1._r8/HTEMP - INV298

    if (HFLAG .gt. 0._r8) then
       !// Hard coded effective HSTAR - Large scale scavenging

       !// Check for hard coded chemical id. If not hardcoded, then stop.
       if (chem_idx .eq. 72) then
          !// SO2

          !// SO2(g) <--> SO2(aq)
          !// A regular Henry expression, which is to be modified
          !// In the convective routine, this is multiplied
          !// by RGAS_MOD*HTEMP
          HTMP1 = KHA * exp(KHB * TEMPFAC)

          !// SO2(aq) <--> HSO3(aq) + H(aq)
          HTMP2 = 1.23e-2_r8 * exp(2.01e3_r8 * TEMPFAC)

          !// Calculate efficient Henry's coefficient
          HSTAR = HTMP1 * (1._r8 + (HTMP2/HPLUS))
       else if (chem_idx .eq. 4) then
          !// HNO3
          !// Regular expression, but enhanced
          HSTAR = 3.80e+2_r8 * KHA * exp(KHB * TEMPFAC)
       else if (chem_idx .eq. 124) then
          !// HNO3s
          !// Regular expression, but enhanced
          HSTAR = 3.80e+2_r8 * KHA * exp(KHB * TEMPFAC)
       else if (chem_idx .eq. 13) then
          !// CH2O
          !// Regular expression, but enhanced
          HSTAR = 1.97_r8 * KHA * exp(KHB * TEMPFAC)
       else
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Hard coded Henry not defined for component ID ',chem_idx
       end if

    else

       !// Regular expression
       HSTAR = KHA * exp(KHB * TEMPFAC)

    end if


    !// --------------------------------------------------------------------
  end subroutine getHstar
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getQTDIS(chem_idx,HSTAR,TM,PR, AMASS,WMASS,TMASS,RETEFF, TDIS)
    !// --------------------------------------------------------------------
    !// This routine works on a sub-box, i.e. the masses are for
    !// the precipitating part of the gridbox.
    !//
    !// HSTAR is effective Henry expression (pre-calculated)
    !// Effective Henry's Law constant:
    !/           H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: chem_idx     !chemistry id (for debugging)
    real(r8), intent(in) :: HSTAR       !Henry's Law expression
    real(r8), intent(in) :: TM, PR      !temperature (K) and pressure (hPa) 
    real(r8), intent(in) :: WMASS       !liquid water mass in box (kg)
    real(r8), intent(in) :: AMASS       !air mass in box (kg)
    real(r8), intent(in) :: TMASS       !tracer in box (kg)
    real(r8), intent(in) :: RETEFF      !retention coefficient
    real(r8), intent(out) :: TDIS       !tracer in aqueous phase (kg)
 
    real(r8) :: TDIS1
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'getQTDIS'
    !//---------------------------------------------------------------------
    !// Fraction dissolved is
    !// TDIS1 = HSTAR * 29.e-6_r8 * PR * WMASS / AMASS
    !// And multiplied by TMASS we get the tracer dissolved:
    TDIS1 = HSTAR * 29.e-6_r8 * PR * WMASS * TMASS / AMASS

    if (TM .lt. TSOLID) then
       !// For true ice, do not allow disllved gases
       TDIS = 0._r8
    else if (TM .lt. TICE) then
       !// For super cooled liquid, allow for retention efficiency < 1
       !// The multiplication of TMASS transforms the fraction dissolved
       !// to the mass dissollved.
       TDIS = RETEFF * TDIS1 * TMASS / (RETEFF * TDIS1 + TMASS)
    else
       !// The multiplication of TMASS transforms the fraction dissolved
       !// to the mass dissollved.
       TDIS = TDIS1 * TMASS / (TDIS1 + TMASS)
    end if

    !// --------------------------------------------------------------------
  end subroutine getQTDIS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getQTDISallT(HSTAR,TM,PR, AMASS,WMASS,TMASS,RETEFF, TDIS)
    !// --------------------------------------------------------------------
    !// This routine works on a sub-box, i.e. the masses are for
    !// the precipitating part of the gridbox.
    !//
    !// HSTAR is effective Henry expression (pre-calculated)
    !// Effective Henry's Law constant:
    !/           H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: HSTAR       !Henry's Law expression
    real(r8), intent(in) :: TM, PR      !temperature (K) and pressure (hPa) 
    real(r8), intent(in) :: WMASS       !liquid water mass in box (kg)
    real(r8), intent(in) :: AMASS       !air mass in box (kg)
    real(r8), intent(in) :: TMASS       !tracer in box (kg)
    real(r8), intent(in) :: RETEFF      !retention coefficient
    real(r8), intent(out) :: TDIS       !tracer in aqueous phase (kg)
 
    real(r8) :: TDIS1
    !// --------------------------------------------------------------------

    TDIS1 = HSTAR * 29.e-6_r8 * PR * WMASS * TMASS / AMASS

    !// allow dissolved mass for all ice and water, i.e. assume
    !// ice to be similar as for super cooled liquid, with
    !// a retention efficiency < 1
    if (TM .lt. TICE) then
       TDIS = RETEFF * TDIS1 * TMASS / (RETEFF * TDIS1 + TMASS)
    else
       TDIS = TDIS1 * TMASS / (TDIS1 + TMASS)
    end if

    !// --------------------------------------------------------------------
  end subroutine getQTDISallT
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DISGAS2(chem_idx,CLWX,CFX,MOLMASS,HSTAR,TM,PR,QM,QT, &
       RETEFF,IT258K,SCViceFR,QTDIS)
    !// --------------------------------------------------------------------
    !// Input is grid box values, except tracer mass which is sub-box
    !// tracer mass. Air mass and liquid water is grid box values and
    !// here we calculate sub-gridbox values (what is available for
    !// scavenging).
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: CLWX,CFX      !cloud water,cloud fraction 
    real(r8), intent(in) :: MOLMASS     !molecular mass of tracer
    real(r8), intent(in) :: HSTAR       !Henry's Law HSTAR
    real(r8), intent(in) :: TM          !temperature of box (K)
    real(r8), intent(in) :: PR          !pressure of box (hPa)
    real(r8), intent(in) :: QM          !air mass in box (kg)
    real(r8), intent(in) :: QT      !tracer in sub-box (kg) *** OBS: SUB BOX ***
    real(r8), intent(in) :: RETEFF      !retention coefficient for ice
    real(r8), intent(in) :: SCViceFR    !fraction of sub-box available for scav.
    integer, intent(in) :: IT258K       !flag for ice treatment below 258K
    integer, intent(in) :: chem_idx     !chemical tracer id (for debugging)
    real(r8), intent(out) :: QTDIS      !tracer dissolved in aqueous phase 
 
    !// Locals
    real(r8)  MUEMP,TDIS1
    real(r8)  TMASS,AMASS,WMASS
    !// Parameters
    !real(r8), parameter  :: TICE   = 273._r8 ! Set in module
    !real(r8), parameter  :: TSOLID = 258._r8 ! Set in module
    !// --------------------------------------------------------------------

    !// Convert values to sub-box
    !// WMASS is total liquid water, so it must be calculated from
    !// gridbox mass (not multiplied by CFX).
    AMASS = QM*CFX   ! kg of air in sub-box
    WMASS = CLWX*QM  ! kg of liquid water in sub-box
    TMASS = QT       ! kg of tracer in sub-box

    !// IMPORTANT for T<273K
    !// The tracer disolved in aqueous phase is QTDIS (kg), so if
    !// only a fraction SCViceFR of the sub-box is available for ice
    !// scavenging it should be multiplied with SCViceFR. This will be
    !// carried out at the end for TM<TICE.

    !// Calculate QTDIS, usually from HSTAR

    if (TM .ge. TSOLID) then

       !// Treatment above TSOLID
       call getQTDIS(chem_idx,HSTAR,TM,PR, AMASS,WMASS,TMASS,RETEFF, QTDIS)

    else

       if (IT258K .eq. 1) then
          !// For T below solid phase, only calculate buried HNO3
          !// (MU = mole fraction) from empirical fit to
          !// Kaercher and Voigt (2006)
          MUEMP = exp(-14.2252_r8 + 0.155704_r8*TM - 7.1929e-4_r8*(TM**2))
          TDIS1 = MUEMP*(MOLMASS/18._r8)*WMASS
          QTDIS = TDIS1*TMASS/(TDIS1 + TMASS)
       else if (IT258K .eq. 2) then
          !// Assume Henry's law to apply also below 258K, using
          !// retention coefficient.
          call getQTDISallT(HSTAR,TM,PR,AMASS,WMASS,TMASS,RETEFF, QTDIS)
       else if (IT258K .eq. 3) then
          !// No removal below 258K, but sets retention coefficient
          !// to 1 for 258K<T<273K.
          QTDIS = 0._r8
       else if (IT258K .eq. 4) then
          !// Assume Henry's law to apply also below 258K, using
          !// retention coefficient = 1.
          !// This is the same as IT258K=2; when RETEFF is eventually read
          !// from scavenging file, this can be removed and IT258K=2 can
          !// be used.
          call getQTDISallT(HSTAR,TM,PR,AMASS,WMASS,TMASS,RETEFF, QTDIS)
       else
          !// May include other treatments here
          QTDIS = 0._r8
       end if

    end if

    !// Adjust to reduce the mass dissolved:
    if (TM .lt. TICE) QTDIS = QTDIS * SCViceFR

    !// --------------------------------------------------------------------
  end subroutine DISGAS2
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine WASHGAS2(chem_idx,RWASH,BOXF,DTSCAV,QTRTOP,HSTAR_IN,TM,PR,QM, &
       QT,RETEFF, QTWASH,QTEVAP)
    !// --------------------------------------------------------------------
    !// In this routine the QM is gridbox air mass, while QTWASH is the
    !// sub-box tracer mass available for washout.
    !//
    !// For most gases below-cloud washout assume Henry-Law equilib with
    !// precip, assumed that precip is liquid.
    !// If frozen, do not call this sub.
    !// Since solubility is moderate, fraction of box with rain does not
    !// matter.
    !//NB this code does not consider the aqueous dissociation (eg, C-q) 
    !// that makes uptake of HNO3 and H2SO4 so complete.  To do so would
    !// require that we keep track of the pH of the falling rain.
    !// THUS the Henry's Law coefficient KHA needs to be enhanced to incldue
    !// this!
    !// ALSO the possible formation of other soluble species from eg, CH2O,
    !// H2O2, can be considered with enhanced values of KHA.
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: chem_idx ! chemistry id (for debugging)
    real(r8), intent(in)  :: RWASH   ! precip leaving bottom of box (kg/s)
    real(r8), intent(in)  :: BOXF    ! fraction of box with washout
    real(r8), intent(in)  :: DTSCAV  ! time step (s)
    real(r8), intent(in)  :: QTRTOP  ! tracer-T in rain entering top of box 
                                     !  over time step (kg)
    real(r8), intent(in)  :: HSTAR_IN ! Henry's Law expression
    real(r8), intent(in)  :: TM      ! temperature of box (K)
    real(r8), intent(in)  :: PR      ! pressure of box (hPa)
    real(r8), intent(in)  :: QT      ! tracer in sub-box (kg)
    real(r8), intent(in)  :: QM      ! air mass in box (kg)
    real(r8), intent(in)  :: RETEFF  ! retention coeff, for scavenging on ice
    real(r8), intent(out) :: QTWASH  ! tracer picked up by precip (kg)
    real(r8), intent(out) :: QTEVAP  ! tracer evaporated from precip (kg)
      
    real(r8), parameter ::  INV298 = 1._r8/298._r8
    real(r8) :: HSTAR, AMASS, FWASH, QTMAX, QTDIF
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'WASHGAS2'
    !//---------------------------------------------------------------------

    AMASS = QM * BOXF

    !// effective Henry's Law constant:
    !//    H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000

    if (TM .le. TSOLID) then
       !// Set HSTAR=0 for T < TSOLID
       HSTAR = 0._r8
    else if  (TM .lt. TICE) then
       !// Use retention coefficient for TSOLID<T<TICE
       HSTAR = RETEFF * HSTAR_IN
    else
       HSTAR = HSTAR_IN
    end if

    !// effective tracer washout frequency (1/s):
    FWASH = HSTAR * 29.e-6_r8 * PR * RWASH / AMASS

    !// equilib amount of T (kg) in rain thru bottom of box over time step
    QTMAX = QT * FWASH * DTSCAV

    if (QTMAX .gt. QTRTOP) then
       !// more of tracer T can go into rain
       QTDIF = min (QT, QTMAX - QTRTOP)
       QTWASH = QTDIF * (1._r8 - exp(-DTSCAV * FWASH))
       QTEVAP = 0._r8
    else
       !// too much of T in rain, must degas/evap T
       QTWASH = 0._r8
       QTEVAP = QTRTOP - QTMAX
    end if
     
    !// --------------------------------------------------------------------
  end subroutine WASHGAS2
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine DISGAS(CLWX,CFX,MOLMASS,KHA,KHB,CQA,TM,PR,QM,QT, &
       RETEFF,IT258K,SCViceFR,QTDIS)
    !// --------------------------------------------------------------------
    !// Input is grid box values, except tracer mass which is sub-box
    !// tracer mass. Air mass and liquid water is grid box values and
    !// here we calculate sub-gridbox values (what is available for
    !// scavenging).
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: CLWX,CFX      !cloud water,cloud fraction 
    real(r8), intent(in) :: MOLMASS     !molecular mass of tracer
    real(r8), intent(in) :: KHA,KHB,CQA !Henry's Law coeffs A*exp(-B/T)
    real(r8), intent(in) :: TM          !temperature of box (K)
    real(r8), intent(in) :: PR          !pressure of box (hPa)
    real(r8), intent(in) :: QM          !air mass in box (kg)
    real(r8), intent(in) :: QT      !tracer in sub-box (kg) *** OBS: SUB BOX ***
    real(r8), intent(in) :: RETEFF      !retention coefficient for ice
    real(r8), intent(in) :: SCViceFR    !fraction of sub-box available for scav.
    integer, intent(in):: IT258K      !flag for ice treatment below 258K
    real(r8), intent(out) :: QTDIS      !tracer dissolved in aqueous phase 
 
    !// Locals
    real(r8)  MUEMP,TDIS1
    real(r8)  TMASS,AMASS,WMASS
    !// Parameters
    !real(r8), parameter  :: TICE   = 273._r8 ! Set in module
    !real(r8), parameter  :: TSOLID = 258._r8 ! Set in module
    !// --------------------------------------------------------------------

    !// Convert values to sub-box
    !// WMASS is total liquid water, so it must be calculated from
    !// gridbox mass (not multiplied by CFX).
    AMASS = QM*CFX   ! kg of air in sub-box
    WMASS = CLWX*QM  ! kg of liquid water in sub-box
    TMASS = QT       ! kg of tracer in sub-box

    !// IMPORTANT for T<273K
    !// The tracer disolved in aqueous phase is QTDIS (kg), so if
    !// only a fraction SCViceFR of the sub-box is available for ice
    !// scavenging it should be multiplied with SCViceFR. This will be
    !// carried out at the end for TM<TICE.

    !// Treatment above TSOLID
    if (TM .ge. TSOLID) then
       call HENRYS (KHA,KHB,CQA,TM,PR, AMASS,WMASS,TMASS,RETEFF, QTDIS)
    else
       !// Treatment for T below solid phase
       !if (TM .lt. TSOLID) then

       if (IT258K .eq. 1) then
          !// For T below solid phase, only calculate buried HNO3
          !// (MU = mole fraction) from empirical fit to
          !// Kaercher and Voigt (2006)
          MUEMP = exp(-14.2252_r8 + 0.155704_r8*TM - 7.1929e-4_r8*(TM**2))
          TDIS1 = MUEMP*(MOLMASS/18._r8)*WMASS
          QTDIS = TDIS1*TMASS/(TDIS1 + TMASS)
       else if (IT258K .eq. 2) then
          !// Assume Henry's law to apply also below 258K, using
          !// retention coefficient.
          call HENRYSallT(KHA,KHB,CQA,TM,PR, &
                         AMASS,WMASS,TMASS,RETEFF, QTDIS)
       else if (IT258K .eq. 3) then
          !// No removal below 258K, but sets retention coefficient
          !// to 1 for 258K<T<273K (done in WETSET_CTM3).
          QTDIS = 0._r8
       else if (IT258K .eq. 4) then
          !// Assume Henry's law to apply also below 258K, using
          !// retention coefficient = 1.
          !// This is the same as IT258K=2; when RETEFF is eventually read
          !// from scavenging file, this can be removed and IT258K=2 can
          !// be used.
          call HENRYSallT(KHA,KHB,CQA,TM,PR, &
                         AMASS,WMASS,TMASS,RETEFF, QTDIS)
       else
          !// May include other treatments here
          QTDIS = 0._r8
       end if

    end if

    !// Adjust to reduce the mass dissolved:
    if (TM .lt. TICE) QTDIS = QTDIS * SCViceFR

    !// --------------------------------------------------------------------
  end subroutine DISGAS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine HENRYS(KHA,KHB,KAQ,TM,PR, AMASS,WMASS,TMASS,RETEFF, TDIS)
    !// --------------------------------------------------------------------
    !// This routine works on a sub-box, i.e. the masses are for
    !// the precipitating part of the gridbox.
    !//
    !// KHA & KHB are the Henrys Law coefficients (not extrapolated below 273K)
    !// KAQ is a pre-calculated (vs. pH) enhancement due to
    !// dissociation/polymers.
    !// Effective Henry's Law constant:
    !/           H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: KHA,KHB,KAQ !Henry's Law coeffs A*exp(-B/T)
    real(r8), intent(in) :: TM, PR      !temperature (K) and pressure (hPa) 
    real(r8), intent(in) :: WMASS       !liquid water mass in box (kg)
    real(r8), intent(in) :: AMASS       !air mass in box (kg)
    real(r8), intent(in) :: TMASS       !tracer in box (kg)
    real(r8), intent(in) :: RETEFF      !retention coefficient
    real(r8), intent(out) :: TDIS       !tracer in aqueous phase (kg)
 
    real(r8), parameter  :: INV298 = 1._r8/298._r8
    !real(r8), parameter  :: TICE   = 273._r8 ! Set in module
    !real(r8), parameter  :: TSOLID = 258._r8 ! Set in module
    real(r8) :: HTEMP, HSTAR, TDIS1
    !// --------------------------------------------------------------------

    HTEMP = max(TM, TICE)
    HSTAR = KAQ * KHA * exp(KHB/HTEMP - KHB*INV298)
    TDIS1 = HSTAR * 29.e-6_r8 * PR * WMASS * TMASS / AMASS

    !// for true ice, do not allow disllved gases
    if (TM .lt. TSOLID) then
       TDIS = 0._r8
       !// for super cooled liquid, allow for retention efficiency < 1
    else if (TM .lt. TICE) then
       TDIS = RETEFF * TDIS1 * TMASS / (RETEFF*TDIS1 + TMASS)
    else
       TDIS = TDIS1 * TMASS / (TDIS1 + TMASS)
    end if

    !// --------------------------------------------------------------------
  end subroutine HENRYS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine HENRYSallT(KHA,KHB,KAQ,TM,PR, AMASS,WMASS,TMASS, RETEFF,TDIS)
    !// --------------------------------------------------------------------
    !// This routine works on a sub-box, i.e. the masses are for
    !// the precipitating part of the gridbox.
    !//
    !// Same as HENRYS, but uses retention coefficient also below 258K.
    !// KHA & KHB are the Henrys Law coefficients (not extrapolated below 273K)
    !// KAQ is a pre-calculated (vs. pH) enhancement due to
    !// dissociation/polymers.
    !// Effective Henry's Law constant:
    !//          H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: KHA,KHB,KAQ !Henry's Law coeffs A*exp(-B/T)
    real(r8), intent(in) :: TM, PR      !temperature (K) and pressure (hPa) 
    real(r8), intent(in) :: WMASS       !liquid water mass in box (kg)
    real(r8), intent(in) :: AMASS       !air mass in box (kg)
    real(r8), intent(in) :: TMASS       !tracer in box (kg)
    real(r8), intent(in) :: RETEFF      !retention coefficient
    real(r8), intent(out) :: TDIS       !tracer in aqueous phase (kg)
 
    real(r8), parameter  :: INV298 = 1._r8/298._r8
    !real(r8), parameter  :: TICE   = 273._r8 ! Set in module
    real(r8) :: HTEMP, HSTAR, TDIS1
    !// --------------------------------------------------------------------

    HTEMP = max(TM, TICE)
    HSTAR = KAQ * KHA * exp(KHB/HTEMP - KHB*INV298)
    TDIS1 = HSTAR * 29.e-6_r8 * PR * WMASS * TMASS / AMASS

    !// allow dissolved mass for all ice and water, i.e. assume
    !// ice to be similar as for super cooled liquid, with
    !// a retention efficiency < 1
    if (TM .lt. TICE) then
       TDIS = RETEFF * TDIS1 * TMASS / (RETEFF*TDIS1 + TMASS)
    else
       TDIS = TDIS1 * TMASS / (TDIS1 + TMASS)
    end if

    !// --------------------------------------------------------------------
  end subroutine HENRYSallT
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine RAINGAS(RRAIN,DTSCAV,CLWX,CFX,QM,QT,QTDIS,QTRAIN)
    !// --------------------------------------------------------------------
    !// In contrast to DISGAS and WASHGAS, this routine use actual
    !// gridbox mass of QT.
    !//
    !// New trace-gas rainout from large-scale precip with two time scales,
    !// one based on precip formation from cloud water and one based on 
    !// Henry's Law solubility: correct limit for delta-t
    !//
    !// NB this code does not consider the aqueous dissociation (eg, C-q) 
    !// that makes uptake of HNO3 and H2SO4 so complete.  To do so would
    !// require that we keep track of the pH of the falling rain.
    !// THUS the Henry's Law coefficient KHA needs to be enhanced to
    !// incldue this!
    !// ALSO the possible formation of other soluble species from eg, CH2O,
    !// H2O2, can be considered with enhanced values of KHA.
    !//
    !// Does NOT now use RMC (moist conv rain) but could, assuming 30% coverage
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: RRAIN       !new rain formation in box (kg/s)
    real(r8), intent(in) :: DTSCAV      !time step (s)
    real(r8), intent(in) :: CLWX,CFX    !cloud water and cloud fraction
    real(r8), intent(in) :: QM          !air mass in box (kg)
    real(r8), intent(in) :: QT          !tracer in box (kg) 
    real(r8), intent(in) :: QTDIS       !tracer in aqueous phase (kg) 
    real(r8), intent(out):: QTRAIN      !tracer picked up by new rain  

    real(r8)   QTLFRQ                   ! loss freq  (1/s)
    real(r8)   WMASS, TMASS
    !// set a maximum loss freq here of 166 sec = 6 m/s fall vel over 1 km
    real(r8), parameter  :: LFREQ0 = 6.e-3_r8    !max loss freq (/s) 
    real(r8), parameter  :: QMAX  = 0.99_r8    !max amount removed in step
    !// --------------------------------------------------------------------

    WMASS = CLWX*QM      ! kg of liquid water in sub-box
    TMASS = QT*CFX        ! kg of tracer in sub-box

    !// tracer loss frequency (1/s) within a cloud (can be sub-box):
    QTLFRQ = (RRAIN * QTDIS) / (WMASS * TMASS)
    QTLFRQ = min (QTLFRQ, LFREQ0)
    !// amount of tracer leaving in precip:
    QTRAIN = TMASS * (1._r8 - exp(-DTSCAV*QTLFRQ))
    QTRAIN = min(QTRAIN, QMAX * TMASS)

    !// --------------------------------------------------------------------
  end subroutine RAINGAS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine WASHGAS(RWASH,BOXF,DTSCAV,QTRTOP,KHA,KHB,KAQ,TM,PR,QM, &
       QT,RETEFF, QTWASH,QTEVAP)
    !// --------------------------------------------------------------------
    !// In this routine the QM is gridbox air mass, while QTWASH is the
    !// sub-box tracer mass available for washout.
    !//
    !// For most gases below-cloud washout assume Henry-Law equilib with
    !// precip, assumed that precip is liquid.
    !// If frozen, do not call this sub.
    !// Since solubility is moderate, fraction of box with rain does not
    !// matter.
    !//NB this code does not consider the aqueous dissociation (eg, C-q) 
    !// that makes uptake of HNO3 and H2SO4 so complete.  To do so would
    !// require that we keep track of the pH of the falling rain.
    !// THUS the Henry's Law coefficient KHA needs to be enhanced to incldue
    !// this!
    !// ALSO the possible formation of other soluble species from eg, CH2O,
    !// H2O2, can be considered with enhanced values of KHA.
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in)  :: RWASH   ! precip leaving bottom of box (kg/s)
    real(r8), intent(in)  :: BOXF    ! fraction of box with washout
    real(r8), intent(in)  :: DTSCAV  ! time step (s)
    real(r8), intent(in)  :: QTRTOP  ! tracer-T in rain entering top of box 
                                     !  over time step (kg)
    real(r8), intent(in)  :: KHA,KHB,KAQ ! Henry's Law coeffs A*exp(-B/T)
    real(r8), intent(in)  :: TM      ! temperature of box (K)
    real(r8), intent(in)  :: PR      ! pressure of box (hPa)
    real(r8), intent(in)  :: QT      ! tracer in sub-box (kg)
    real(r8), intent(in)  :: QM      ! air mass in box (kg)
    real(r8), intent(in)  :: RETEFF  ! retention coeff, for scavenging on ice
    real(r8), intent(out) :: QTWASH  ! tracer picked up by precip (kg)
    real(r8), intent(out) :: QTEVAP  ! tracer evaporated from precip (kg)
      
    real(r8), parameter ::  INV298 = 1._r8/298._r8
    !real(r8), parameter ::  TICE = 273._r8 ! Set in module
    !real(r8), parameter ::  TSOLID = 258._r8 ! Set in module
    real(r8) :: HSTAR, HTEMP, AMASS, FWASH, QTMAX, QTDIF
    !// --------------------------------------------------------------------

    AMASS = QM * BOXF

    !// effective Henry's Law constant:
    !//    H* = moles-T / liter-precip / press(atm-T)
    !// p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
    if (TM .gt. TSOLID) then
       HTEMP = max(TM, TICE)
       HSTAR = KAQ * KHA * exp(KHB/HTEMP - KHB*INV298)
    else
       HSTAR = 0._r8
    end if
    if  (TM .lt. TICE) then
       HSTAR = RETEFF * HSTAR
    end if
      
    !// effective tracer washout frequency (1/s):
    FWASH = HSTAR * 29.e-6_r8 * PR * RWASH / AMASS

    !// equilib amount of T (kg) in rain thru bottom of box over time step
    QTMAX = QT * FWASH * DTSCAV

    if (QTMAX .gt. QTRTOP) then
       !// more of tracer T can go into rain
       QTDIF = min (QT, QTMAX - QTRTOP)
       QTWASH = QTDIF * (1._r8 - exp(-DTSCAV * FWASH))
       QTEVAP = 0._r8
    else
       !// too much of T in rain, must degas/evap T
       QTWASH = 0._r8
       QTEVAP = QTRTOP - QTMAX
    end if
     
    !// --------------------------------------------------------------------
  end subroutine WASHGAS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DIAMEMP(CWATER,RRATE,DIAM)
    !// --------------------------------------------------------------------
    !// empirical fit of precip DIAMeter (mm) to cloud water dens & rain rate
    !// as given in Field and Heymsfield, J. Atm. Sci., 60, 544-560, 2003,
    !// doi:10.1175/1520-0469(2003)060%3C0544:AASOIC%3E2.0.CO;2
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in)  :: CWATER    !  kg/m3 
    real(r8), intent(in)  :: RRATE     !  kg/m2/s = mm/s
    real(r8), intent(out) :: DIAM  

    real(r8) :: RX,WX,THETA,PHI,ETA,BETA1,ALPHA,BEXP
    real(r8) :: GAMTHETA,GAMBETA1
    !// --------------------------------------------------------------------

    RX = RRATE * 3600._r8               !mm/hr
    WX = CWATER * 1.e3_r8                !g/m3

    if (RX .gt. 0.04_r8) then
       THETA = exp(-1.43_r8 * log10(7._r8 * RX)) + 2.8_r8
    else
       THETA = 5._r8
    end if
    BETA1  = 1._r8 + 0.6105_r8 * THETA       ! 0.6105 = 1/(1+.638)
    call GAMMAX(THETA,GAMTHETA)
    call GAMMAX(BETA1,GAMBETA1)
    PHI   = 0.1_r8 * RRATE
    ETA   = exp(3.01_r8 * THETA - 10.5_r8)
    ALPHA = exp(4._r8 * BETA1 - 18.0_r8)

    BEXP = -1._r8 / (0.3895_r8 * THETA - 1._r8)
    DIAM = ( (WX * ETA * GAMTHETA) &
              / (1.e6_r8 * ALPHA * PHI * GAMBETA1) )**BEXP
    DIAM = 10._r8 * DIAM

    !// --------------------------------------------------------------------
  end subroutine DIAMEMP
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine GAMMAX(XX,GA)
    !// --------------------------------------------------------------------
    !// Returns the value GAMMA(XX) using ln gamma algorithm
    !// bounded by 1.d+30 to +8.d+30  (forces XX > 0)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: XX
    real(r8), intent(out):: GA
    real(r8) :: X,Y,TMP,SER,GAMMLN
    integer :: J
    real(r8), dimension(6), parameter:: COF = &
         [76.18009172947146_r8, -86.50532032941677_r8, &
          24.01409824083091_r8,  -1.231739572450155_r8, &
          0.1208650973866179e-2_r8,  -0.5395239384953e-5_r8]
    real(r8), parameter:: STP = 2.5066282746310005_r8
    real(r8), parameter:: SER0 = 1.000000000190015_r8
    !// --------------------------------------------------------------------

    X = min(30._r8, max(1.e-30_r8, XX))
    Y = X
    TMP = X + 5.5_r8
    TMP =(X + 0.5_r8) * log(TMP) - TMP
    SER = SER0
    do J = 1, 6
       Y = Y + 1._r8
       SER = SER + COF(J) / Y
    end do
    GAMMLN = TMP + log(STP * SER / X)
    GA = exp(GAMMLN)

    !// --------------------------------------------------------------------
  end subroutine GAMMAX
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine WETSET_CTM3(INFILE)
    !// --------------------------------------------------------------------
    !// Sets wet deposition/scavenging for WASH2.
    !// Needs to be called in setup once.
    !//
    !// Modified for CTM3. Reads retention coefficients for all transported
    !// tracers.
    !//
    !// Ole Amund Sovde, March-May 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_ctm, only: NSCX, NDPX, NTM
    use cmn_chem, only: TCAER, TRWETL, TCWETL, TCHENA, TCHENB, TCKAQA, &
         TCKAQB, SCV_RETEFF, SCViceFR,  SCV_T258, TNAME
    use cmn_oslo, only: TCCNVHENRY
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), intent(in) ::  INFILE

    character(len=80) :: TITLE, SCAV_FMT, OUT_FMT
    character(len=10) :: XTRACR
    logical :: LXX_AER
    real(r8) :: XXSOLU,XXHA,XXHB, XXQA,XXQB,XXSCVFR
    real(r8) :: XXION
    integer :: XXCHN, XX258
    integer :: N_TRACERS, I,J,N,NN,LSURF
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'WETSET_CTM3'
    !//---------------------------------------------------------------------

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// Cannot use this treatment for NSCX==0
    if (NSCX .ne. 1) then
       write(6,'(a)') f90file//':'//subr// &
            ': Cannot use WETSET_CTM3 to set WASH1 lifetimes - STOP'
       stop
    end if

    !// Initialize
    TCAER(:) =  .false.
    TCWETL(:) = 0._r8
    TCCNVHENRY(:) = 0         !// for CTM3
    TCHENA(:) = 0._r8
    TCHENB(:) = 0._r8
    TCKAQA(:) = 0._r8
    TCKAQB(:) = 0._r8
    !// CTM3 variable retention coefficients and treatment below 258K
    SCV_RETEFF(:) = 0.5_r8 !// Assume 0.5 as default
    SCViceFR(:) = 0._r8    !// Fraction available for ice scavenging.
    SCV_T258(:) = 0
    LXX_AER = .false.

    !// Zero out wet scavenging table for WASH1
    TRWETL(:,:) = 0._r8

    !// Open scavenging data file
      write(6,*) ' open file: '//trim(INFILE)
    open(8,file=INFILE,status='old',form='formatted')
    read(8,'(a60)') TITLE !// Title 1
      write(6,*) trim(TITLE)
    read(8,'(a60)') TITLE !// Title 2
      write(6,*) trim(TITLE)
    read(8,'(a)') SCAV_FMT
    read(8,SCAV_FMT) N_TRACERS
    !// Print out table values, so output can be checked
    out_fmt='(a,a10,1x,2a4,1x,a4,a4,a8,1x,a5,1x,a4,1x,a4,1x,a6,a4)'
      write(6,out_fmt) '  match: ', 'tracername', &
           '#scv','#ntm','SOLU',' CHN','KHA    ', &
           'KHB  ','KAQA    ','KAQB','ISCVFR','ICE'
    !// Format for numbers
    out_fmt='(a,a10,1x,2i4,1x,f4.2,i3,1x,es8.1,1x,f5.0,1x,f4.0,1x,f4.0,1x,es7.1,i3)'

    do N = 1, N_TRACERS
      read(8,SCAV_FMT) I,XTRACR,XXSOLU,XXCHN,XXHA,XXHB, XXQA,XXQB,XXSCVFR,XX258
      !// Break do-loop at -99
      if (I .eq. -99) exit
      !// Search for current tracers:
      do NN = 1, NTM
        if (TNAME(NN) .eq. XTRACR) then
          write(6,out_fmt) '  match: ', XTRACR,I,NN,XXSOLU,XXCHN, &
               XXHA,XXHB, XXQA,XXQB,XXSCVFR,XX258
          TCAER(NN) = LXX_AER
          TCWETL(NN) = XXSOLU
          TCCNVHENRY(NN) = XXCHN !// for CTM3
          !// Check WETL against CHN
          if ((XXSOLU.eq.0._r8) .and. (XXCHN.ne.0)) then
             print*,'Stopping: Check SOLU against CHN'
             print*,'XXSOLU:',XXSOLU
             print*,'XXCHN: ',XXCHN
             stop
          end if
          if ((XXSOLU.eq.0._r8) .and. (XXHA.gt.0)) then
             print*,'WARNING: Wet savenging due to Henrys law zeroed'
             print*,'XXSOLU:',XXSOLU
             print*,'XXHA/XXHB: ',XXHA,XXHB
          end if
          TCHENA(NN) = XXHA
          TCHENB(NN) = XXHB
          TCKAQA(NN) = XXQA
          TCKAQB(NN) = XXQB
          !// Fraction of grid box available for ice scavenging
          SCViceFR(NN)   = max(min(XXSCVFR, 1._r8), 0._r8)
          !// Ice uptake retention coefficient
          SCV_RETEFF(NN) = 0.5_r8
          if (XX258 .eq. 3 .or. XX258 .eq. 4) then
             !// 3: No removal below 258K, 4: remove as Henry's law below 258K
             SCV_RETEFF(NN) = 1._r8
          end if
          !// Treatment below 258K
          SCV_T258(NN) = XX258
        end if
      end do
    end do
    close(8)

    !// --------------------------------------------------------------------
  end subroutine WETSET_CTM3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module scavenging_largescale_uci
!//=========================================================================

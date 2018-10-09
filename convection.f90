!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Convection routine.
!//=========================================================================
module convection
  !// ----------------------------------------------------------------------
  !// MODULE: convection
  !// DESCRIPTION: Routine for convection and convective removal.
  !//              This is the UCI p-cnvw.f, modified for CTM3 and converted
  !//              to f90.
  !// ---(p-cnvw.f)-----UCIrvine CTM  p-code 5.5 (9/2007)
  !// Version: qcode_56d; 20090318
  !// CTM3: Modified to CTM2 convective rainout.
  !// wet(moist) convection code
  !// Contains:
  !//   subroutine CONVW_OSLO
  !//   subroutine ADJFLX2
  !//   subroutine ADJFLX
  !//   subroutine QCNVW2_OSLO
  !//   subroutine QCNVW2 (UCI)
  !//   subroutine SCAV_UPD (UCI)
  !//
  !// To f90 for CTM3: Ole Amund Sovde, May 2015
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'convection.f90'
  !// ----------------------------------------------------------------------
  private
  public CONVW_OSLO
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine CONVW_OSLO(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                         AIRB,BTEM,GAMACB,LPAUZB,DTCONV,NOPS,NSUB,MP)
    !// --------------------------------------------------------------------
    !// NEW convection code, CFL limit has been dropped, extra mass flux
    !// is taken from next above layer *IF* needed, similar effect as
    !// subsidence.
    !// New sub ADJFLX needed to adjust updraft fluxes & check for conv.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LWEPAR, LWDPAR, LCONVM
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, NTM
    use cmn_chem, only: TCAER, TCHENA, N_STE, TMASSMIX2MOLMIX, &
         O3iso1, O3iso2, LPAUZTOP
    use cmn_diag, only: O3WSCAV
    use cmn_met, only: CWETN, CWETE, CWETD, CENTU, CENTD, PRECCNV, ZOFLE
    use cmn_oslo, only: QFRAC, LW_VOLCONC, TCCNVHENRY, CONVWASHOUT, &
         LELEVTEMP
    use cnv_oslo, only: elevator_fractions, liquid_fractions
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NOPS,NSUB,MP
    real(r8), intent(in) :: DTCONV
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM
    logical, dimension(LPAR,IDBLK,JDBLK), intent(in) :: LPAUZB
    !// Inout/Output
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(inout) :: AIRB
    !// Output
    real(r8), intent(out), dimension(LPAR,IDBLK,JDBLK)   :: GAMACB
    !// --------------------------------------------------------------------
    !// Locals
    integer I,J,L,N,LE,LDN, II,JJ
    real(r8), dimension(LPAR+1) :: &
         FLUX_N, FLUX_E, FLUX_D, ENT_D, ENT_U
    real(r8), dimension(LPAR) ::    TEM, QMOLD, QM, &
         QTT, QXT, QXX, QXY, QXZ, QYT, QYY, QYZ, QZT, QZZ
    real(r8) :: DELZ, SUM_N, SUM_E, SUM_D, SUM_S

    logical :: LAER, LIJ

    !// CTM3 wetloss
    real(r8), dimension(LWEPAR) :: CRAIN, TWETLN, CNVLOST
    real(r8), dimension(LWEPAR,NPAR) :: CNV_WETL
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CONVW_OSLO'
    !// --------------------------------------------------------------------

    LE      = min(LCONVM,LWEPAR)
    LDN     = min(LE,LWDPAR)
    if (mod(LE,2) .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': LE is odd; not allowed!'
       stop 'STOP in '//subr
    end if

    GAMACB(:,:,:) = 0._r8

    !// --------------------------------------------------------------------
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ   = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
         II   = I - MPBLKIB(MP) + 1

        !// Initialize
        do L = 1,LPAR
          QM(L)    = AIRB(L,II,JJ) !// Air
          QMOLD(L) = AIRB(L,II,JJ)
          TEM(L)   = BTEM(L,II,JJ) !// Temperature
        end do
        do L = 1,LPAR+1
          FLUX_N(L) = 0._r8        !// Non-entraining flux updrafts (kg/s)
          FLUX_D(L) = 0._r8        !// Flux downdrafts (kg/s)
          FLUX_E(L) = 0._r8        !// Entraining flux updrafts (kg/s)
          ENT_U(L)  = 0._r8        !// Entrainment into updrafts (kg/s)
          ENT_D(L)  = 0._r8        !// Entrainment into downdrafts (kg/s)
        end do

!// Convective mass fluxes => mass moved per time step, conv-CFL dropped
        !// Set updraft fluxes/entrainment
        do L = 1,LE
          FLUX_N(L) = CWETN(I,J,L) * DTCONV !// ECMWF-field = 0.
          FLUX_E(L) = CWETE(I,J,L) * DTCONV
          ENT_U(L)  = CENTU(I,J,L) * DTCONV
          CRAIN(L) = PRECCNV(I,J,L) !// kg/s
        end do
        !// This is normally not done, since LE is equal to LWEPAR
        do L = LE+1, LWEPAR
          CRAIN(L) = 0._r8
        end do
        !// Downdraft flux/entrainment
        do L = 1,LDN
          FLUX_D(L) = CWETD(I,J,L) * DTCONV
          ENT_D(L)  = CENTD(I,J,L) * DTCONV
        end do


!// Adjust fluxes if need be to take full DTCONV time step
        call ADJFLX2(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D,QM,LE, &
             SUM_N,SUM_E,SUM_D,SUM_S, LIJ)

        !// No convection, go to next I,J
        !// ----------------------------------------------------------------
        if (.not.LIJ) cycle


        !// CTM3: Calculate fractions for wet loss:
        !//   QFRAC:      Fraction of convective rain [kg/s] to liquid
        !//               cloud water [kg/s]. Only done for NOPS=1 and
        !//               NSUB=1, and are thus constant through
        !//               meteorological time step.
        !//   LW_VOLCONC: Fraction of liquid cloud water volume to
        !//               total elevator volume. Only done for NOPS=1
        !//               and NSUB=1.
        !//   CNV_WETL:   Fraction of tracers lost in convective rainout.
        if (NOPS.eq.1 .and. NSUB.eq.1) then
          !// QFRAC and LW_VOLCONC only depend on meteorology.
          !// IMPORTANT: ENT_U and FLUX_E are now [kg], while
          !//            CRAIN is [kg/s]. The routine will handle
          !//            the different units.
          call elevator_fractions(I,J,II,JJ,MP, LE, &
               BTEM, CRAIN, ENT_U, FLUX_E, DTCONV)
        end if

        !// Calculate liquid fractions (all time steps) and sets
        !// CNV_WETL(L,N), which are used instead of TCWETL further
        !// down in the code.
        call liquid_fractions(I,J,II,JJ,MP, BTEM, CNV_WETL)

        !// Possible override to turn off convective washout
        !// is to set CNV_WETL(:,:) = 0._r8

        !// QCNVW2 is refined with respect to CNV_WETL and TCWETL, and
        !// the CTM2 method of removal replaces SCAV_UPD

        do N = 1,NTM
          !// CTM3 wetloss parameter
          !TWETLN = TCWETL(N)
          TWETLN(:) = CNV_WETL(:,N)
          LAER   = TCAER(N)
          do L = 1,LPAR
            QM(L)  = QMOLD(L)
            QTT(L) = BTT(L,N,II,JJ)
            QXT(L) = BXT(L,N,II,JJ)
            QYT(L) = BYT(L,N,II,JJ)
            QZT(L) = BZT(L,N,II,JJ)
            QXX(L) = BXX(L,N,II,JJ)
            QYY(L) = BYY(L,N,II,JJ)
            QZZ(L) = BZZ(L,N,II,JJ)
            QXY(L) = BXY(L,N,II,JJ)
            QXZ(L) = BXZ(L,N,II,JJ)
            QYZ(L) = BYZ(L,N,II,JJ)
          end do


          call QCNVW2_OSLO(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D, ZOFLE(1,I,J), &
               QM,QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ, &
               SUM_N,SUM_E,SUM_D,SUM_S, TWETLN,CNVLOST,CRAIN,TCHENA(N), &
               LAER,LE)

          do L = 1,LE
            BTT(L,N,II,JJ) = QTT(L)
            BXT(L,N,II,JJ) = QXT(L)
            BYT(L,N,II,JJ) = QYT(L)
            BZT(L,N,II,JJ) = QZT(L)
            BXX(L,N,II,JJ) = QXX(L)
            BYY(L,N,II,JJ) = QYY(L)
            BZZ(L,N,II,JJ) = QZZ(L)
            BXY(L,N,II,JJ) = QXY(L)
            BXZ(L,N,II,JJ) = QXZ(L)
            BYZ(L,N,II,JJ) = QYZ(L)
          end do

          !// Sum up convective loss (save as negative values)
          !// Note that at lower altitudes, CNVLOST may be negative due
          !// to evaporation.
          do L = 1, LE
            CONVWASHOUT(L,N,II,JJ,MP) = CONVWASHOUT(L,N,II,JJ,MP) &
                  - CNVLOST(L)
          end do
          !// Save scavenging of O3
          if (N .eq. N_STE) then
            do L = 1,LE
              if (BTT(L,N_STE,II,JJ) * TMASSMIX2MOLMIX(N_STE) &
                        .le. AIRB(L,II,JJ) * O3iso1 ) &
                   O3WSCAV(I,J,1) = -cnvlost(L) + O3WSCAV(I,J,1)
              if (BTT(L,N_STE,II,JJ) * TMASSMIX2MOLMIX(N_STE) &
                        .le. AIRB(L,II,JJ) * O3iso2 ) &
                   O3WSCAV(I,J,2) = -cnvlost(L) + O3WSCAV(I,J,2)
              if (L.le.LPAUZTOP(I,J) .and. (.not.LPAUZB(L,II,JJ))) &
                   O3WSCAV(I,J,3) = -cnvlost(L) + O3WSCAV(I,J,3)
            end do
          end if
        end do

!// Caclulate the conv-subsidence w, combine with large-scale w later
        !// The subsidence mass equals the mass lost by convection, i.e.
        !// no net change in mass due to convection after subsidence.
        do L = LE, 2, -1
          GAMACB(L,II,JJ) = GAMACB(L+1,II,JJ) + (QMOLD(L) - QM(L))
        end do
        do L = 1, LE
          AIRB(L,II,JJ) = QM(L)                      !// New mass
          GAMACB(L,II,JJ) = GAMACB(L,II,JJ) / DTCONV !// To w (kg/s)
        end do

      end do     ! II block
    end do       ! JJ block

    !// --------------------------------------------------------------------
  end subroutine CONVW_OSLO
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ADJFLX(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D,QM,NQ, &
       SUM_N,SUM_E,SUM_D,SUM_S, LIJ)
    !// --------------------------------------------------------------------
    !// Adjust fluxes so that no more than a fraction CNST of the grid
    !// box mass is moved each time step.
    !// Downdraft flux is also limited to a fraction UDFRAC of air
    !// detrained by updrafts.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Inout
    integer, intent(in) :: NQ
    real(r8), intent(in), dimension(LPAR) :: QM
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR+1) :: &
         FLUX_N, FLUX_E, FLUX_D, ENT_U, ENT_D
    !// Output
    real(r8), intent(out) :: SUM_N, SUM_E, SUM_D, SUM_S
    logical, intent(out) :: LIJ
    !// --------------------------------------------------------------------
    !// Locals
    real(r8), dimension(101) :: &
         FLUX_N0, FLUX_E0, FLUX_D0, ENT_U0, ENT_D0
    real(r8), dimension(100) :: QMCL, DQM, TOTQM
    real(r8) :: DFLUX, FLXD_max
    integer :: L

    !// Max fraction allowed to be removed from grid box
    real(r8), parameter :: CNST = 0.85_r8
    !// Limit downdraft flux to max fraction UDFRAC of the air detrained
    !// by updrafts.
    real(r8), parameter :: UDFRAC = 0.25_r8
    !// --------------------------------------------------------------------

    !// Initialize
    LIJ = .false.
    SUM_N      = 0._r8
    SUM_E      = 0._r8
    SUM_D      = 0._r8
    SUM_S      = 0._r8

    !// Sum up the fluxes
    do L = 1,NQ
       SUM_N    = SUM_N + FLUX_N(L)
       SUM_E    = SUM_E + FLUX_E(L) + ENT_U(L)
       SUM_D    = SUM_D - FLUX_D(L) + ENT_D(L)
    end do
    SUM_S = SUM_N + SUM_E + SUM_D

    if (SUM_S .lt. 1.e-30_r8) return   ! if no convection, exit


    LIJ = .true.

    !// Save mass and fluxes before adjustment
    do L = 1,LPAR
       QMCL(L)    = QM(L)
       FLUX_E0(L) = FLUX_E(L)
       FLUX_D0(L) = FLUX_D(L)
       FLUX_N0(L) = FLUX_N(L)
       ENT_U0(L)  = ENT_U(L)
       ENT_D0(L)  = ENT_D(L)
       DQM(L)    = 0._r8      !// Increase in elevator due to net in-flux
    end do

    if (SUM_N .ge. 1.e-10_r8) then
      !// Do non-entraining updrafts first 
      FLUX_N(1)     = 0._r8
      do L = 1,NQ
        !// Calculate net flux upwards
        DFLUX         = FLUX_N0(L+1) - FLUX_N(L)
        if (DFLUX .gt. 0._r8)  then
          !// More out on top than in at bottom; make sure no more than CNST
          !// is moved. First do flux:
          FLUX_N(L+1) = min(DFLUX, CNST*QMCL(L)) + FLUX_N(L)
          !// then adjust air mass
          QMCL(L)     = QMCL(L) - (FLUX_N(L+1) - FLUX_N(L))
        else
          !// Less out on top than in at bottom; net increase of mass inside
          !// box.
          !// Not necessary: FLUX_N(L+1) = FLUX_N0(L+1)
          !// Increase in air mass in convective plume
          DQM(L)      = DQM(L) - DFLUX
        end if
      end do
    end if !// non-entraining updrafts

    if (SUM_E .ge. 1.e-10_r8) then
      !// Next do entraining updrafts
      do L = 1,NQ    
        !// Limit entrainment to CNST of grid box mass
        ENT_U(L)    = min(ENT_U0(L), CNST*QMCL(L))
      end do
      !// Make sure no mass comes from the surface
      FLUX_E(1)     = 0._r8

      do L = 1,NQ
        !// Calculate net flux upwards
        DFLUX         = FLUX_E0(L+1) - FLUX_E(L)
        if (DFLUX .gt. 0._r8)  then
          !// More out on top than in at bottom; make sure no more than CNST
          !// is moved. First do flux:
          FLUX_E(L+1) = min(DFLUX, CNST*QMCL(L)) + FLUX_E(L)
          !// then adjust air mass
          QMCL(L)     = QMCL(L) - (FLUX_E(L+1)-FLUX_E(L))
        else
          !// Less out on top than in at bottom; net increase of mass inside
          !// box.
          !// Not necessary: FLUX_E(L+1) = FLUX_E0(L+1)
          !// Increase in air mass in convective plume
          DQM(L)      = DQM(L) - DFLUX
        end if
      end do
    end if !// Entraining updrafts

    if (abs(SUM_D) .ge. 1.e-10_r8) then
      !// Last do downdrafts
      do L = 1,NQ
        !// Total air; ambient + convective plume
        TOTQM(L) = QMCL(L) + DQM(L)
        !// Entrain max what is available
        ENT_D(L) = min(ENT_D0(L), TOTQM(L))
      end do
      !// No flux in at the top
      FLUX_D(NQ+1) = 0._r8

      do L = NQ,1,-1        ! from top down    
        !// Limit influx to CNST of total mass
        !// Also so that no more than 25% of the detrained upward flux
        !// (not from the detrained diag, but from the divergence of
        !// the upward flux left in the layer):
        !//
        !// Check for net downward flux
        DFLUX = max(-FLUX_D0(L) + FLUX_D(L+1), 0._r8)
        if (DFLUX .gt. 1.e-10_r8) then !0.d0) then
           !// There is net downward flux; allow max 25% of detrained
           !//  upward flux. Only allow F(L+1)-F(L) < 0 (i.e. detrainment)
           FLXD_max = - min( 0._r8, UDFRAC * (FLUX_N(L+1) - FLUX_N(L) &
                                             + FLUX_E(L+1) - FLUX_E(L)) )
        else
           !// No flux, no need to limit to detrained value.
           FLXD_max = 0._r8
        end if
        !// Limit flux
        FLXD_max = min(FLXD_max, - FLUX_D0(L) + FLUX_D(L+1), CNST*TOTQM(L))
        FLUX_D(L) = FLUX_D(L+1) - FLXD_max
        !// Old treatment
        !// FLUX_D(L) = FLUX_D(L+1) &
        !//       - min(-FLUX_D0(L)+FLUX_D(L+1), CNST*TOTQM(L))
      end do

    end if !// Downdrafts

    !// --------------------------------------------------------------------
  end subroutine ADJFLX
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ADJFLX2(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D,QM,NQ, &
       SUM_N,SUM_E,SUM_D,SUM_S, LIJ)
    !// --------------------------------------------------------------------
    !// Adjust fluxes so that no more than a fraction CNST of the grid
    !// box mass is moved each time step.
    !// Downdraft flux is also limited to a fraction UDFRAC of air
    !// detrained by updrafts.
    !// According to ECMWF (Peter Bechtold) downdraft and detrainment
    !// should be scaled with the scaling resulting from updraft flux
    !// adjustment. This is available in ZMFS.
    !//
    !// Updated from UCI adjflx.
    !//
    !// Ole Amund Sovde, March 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NQ
    real(r8), intent(in), dimension(LPAR) :: QM
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR+1) :: &
         FLUX_N, FLUX_E, FLUX_D, ENT_U, ENT_D
    !// Output
    real(r8), intent(out) :: SUM_N, SUM_E, SUM_D, SUM_S
    logical, intent(out) :: LIJ
    !// --------------------------------------------------------------------
    real(r8), dimension(101) :: &
         FLUX_N0, FLUX_E0, FLUX_D0, ENT_U0, ENT_D0
    real(r8), dimension(100) :: QMCL, DQM, TOTQM
    real(r8) :: DFLUX, FLXD_max, ZMFS(LPAR)
    integer :: L

    !// Max fraction allowed to be removed from grid box
    real(r8), parameter :: CNST = 0.8_r8
    !// Limit downdraft flux to max fraction UDFRAC of the air detrained
    !// by updrafts.
    real(r8), parameter :: UDFRAC = 0.25_r8
    !// --------------------------------------------------------------------

    LIJ = .false.
    SUM_N      = 0._r8
    SUM_E      = 0._r8
    SUM_D      = 0._r8
    SUM_S      = 0._r8
    !// scale for entraining updrafts
    ZMFS(:) = 1._r8
    !// Sum up the fluxes
    do L = 1,NQ
       SUM_N    = SUM_N + FLUX_N(L)
       SUM_E    = SUM_E + FLUX_E(L) + ENT_U(L)
       SUM_D    = SUM_D - FLUX_D(L) + ENT_D(L)
    end do
    SUM_S    = SUM_N + SUM_E + SUM_D

    if (SUM_S .lt. 1.e-30_r8) return  ! if no convection, exit

    LIJ = .true.

    !// Save mass and fluxes before adjustment
    do L = 1, LPAR
       QMCL(L)    = QM(L)
       FLUX_E0(L) = FLUX_E(L)
       FLUX_D0(L) = FLUX_D(L)
       FLUX_N0(L) = FLUX_N(L)
       ENT_U0(L)  = ENT_U(L)
       ENT_D0(L)  = ENT_D(L)
       DQM(L)    = 0._r8      !// Increase in elevator due to net in-flux
    end do

    if (SUM_N .ge. 1.e-10_r8) then
      !// Do non-entraining updrafts first 
      FLUX_N(1)     = 0._r8
      do L = 1,NQ
        !// Calculate net flux upwards
        DFLUX         = FLUX_N0(L+1) - FLUX_N(L)
        if (DFLUX .gt. 0._r8)  then
          !// More out on top than in at bottom; make sure no more than CNST
          !// is moved. First do flux:
          FLUX_N(L+1) = min(DFLUX, CNST*QMCL(L)) + FLUX_N(L)
          !// then adjust air mass
          QMCL(L)     = QMCL(L) - (FLUX_N(L+1) - FLUX_N(L))
        else
          !// Less out on top than in at bottom; net increase of mass inside
          !// box.
          !// Not necessary: FLUX_N(L+1) = FLUX_N0(L+1)
          !// Increase in air mass in convective plume
          DQM(L)      = DQM(L) - DFLUX
        end if
      end do
    end if !// Non-entraining updrafts

    if (SUM_E .ge. 1.e-10_r8) then
      !// Next do entraining updrafts
      do L = 1,NQ    
        !// Limit entrainment to CNST of grid box mass
        ENT_U(L)    = min(ENT_U0(L), CNST*QMCL(L))
      end do
        !// Make sure no mass comes from the surface
        FLUX_E(1)     = 0._r8
      do L = 1,NQ
          !// Calculate net flux upwards
          DFLUX         = FLUX_E0(L+1) - FLUX_E(L)
        if (DFLUX .gt. 0._r8)  then
          !// More out on top than in at bottom; make sure no more than CNST
          !// is moved. First do flux:
          FLUX_E(L+1) = min(DFLUX, CNST*QMCL(L)) + FLUX_E(L)
          !// Scale for updraft flux
          ZMFS(L) = min(ZMFS(L),min(DFLUX, CNST*QMCL(L))/DFLUX)
          !// then adjust air mass
          QMCL(L)     = QMCL(L) - (FLUX_E(L+1)-FLUX_E(L))
        else
          !// Less out on top than in at bottom; net increase of mass inside
          !// box.
          !// Not necessary: FLUX_E(L+1) = FLUX_E0(L+1)
          !// Increase in air mass in convective plume
          DQM(L)      = DQM(L) - DFLUX
        end if
      end do
    end if !// Entraining updrafts

    if (abs(SUM_D) .ge. 1.e-10_r8) then
      !// Last do downdrafts
      do L = 1,NQ
        !// Total air; ambient + convective plume
        TOTQM(L) = QMCL(L) + DQM(L)
        !// Entrain max what is available
        !ENT_D(L) = min(ENT_D0(L), TOTQM(L))
        !// Scale down detrainment to scale factor from updrafts
        ENT_D(L) = min(ENT_D0(L)*ZMFS(L), TOTQM(L))
      end do
        !// No flux in at the top
        FLUX_D(NQ+1) = 0._r8
      do L = NQ,1,-1        ! from top down    
        !// Limit influx to CNST of total mass
        !// Also so that no more than 25% of the detrained upward flux
        !// (not from the detrained diag, but from the divergence of
        !// the upward flux left in the layer):
        !//
        !// Scale down downdrafts to scale factor from updrafts
        FLUX_D0(L) = FLUX_D0(L)*ZMFS(L)
        !// Check for net downward flux coming in
        DFLUX = max(-FLUX_D(L+1), 0._r8)
        !DFLUX = max(-FLUX_D0(L)+FLUX_D(L+1),0._r8)
        if (DFLUX .gt. 0._r8) then
           !// There is downward flux coming in; allow max 25% of detrained
           !// upward flux. Only allow F(L+1)-F(L) < 0 (i.e. detrainment)
           FLXD_max = max( 0._r8, - UDFRAC * (FLUX_N(L+1) - FLUX_N(L) &
                                            + FLUX_E(L+1) - FLUX_E(L)) )
        else
           !// No flux, no need to limit to detrained value.
           FLXD_max = 0._r8
        end if
        !// Limit flux; remember that values are positive
        FLXD_max = min(FLXD_max, -FLUX_D0(L) + FLUX_D(L+1), CNST*TOTQM(L))
        FLUX_D(L) = FLUX_D(L+1) - FLXD_max
      end do

    end if !// Downdrafts

    !// --------------------------------------------------------------------
  end subroutine ADJFLX2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine QCNVW2_OSLO(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D,ZOFL, &
       QM,QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ, &
       SUM_N,SUM_E,SUM_D,SUM_S, TWETLN,CNVLOST,CRAIN,HENRY_A, &
       LAER,NQ)
    !// --------------------------------------------------------------------
    !// Version: qcode_56d; 20090318
    !// CTM3: Modified to CTM2 convective rainout.
    !//       TWETLN changed to dimension(LPAR) and tracers are scavenged
    !//       by removing this fraction from the plume traveling upwards.
    !//
    !// Ole Amund Sovde, May 2009
    !// --------------------------------------------------------------------
    !// New version drops the CFL criteria and just removes up to 99% of
    !// air mass in box, see ADJFLX.
    !// Also w* subsidence has ben moved out and combine with W-V-U sequence
    !// Processes:
    !//    1 = N = non-entraining updrafts
    !//    2 = E = entraining updrafts
    !//    3 = D = downdrafts
    !// Explicitly calculates transport in a column due to wet convection.
    !// Includes entraining updrafts/downdrafts, nonentraining updrafts.
    !// Calls scavenging routine SCAV_UPD for both updraft types.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, LWEPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in), dimension(LPAR+1) :: &
         FLUX_N, FLUX_E, FLUX_D, ENT_U, ENT_D, ZOFL
    real(r8), intent(in)  :: SUM_N, SUM_E, SUM_D, SUM_S, HENRY_A
    real(r8), intent(in), dimension(LWEPAR)  :: TWETLN, CRAIN
    logical, intent(in) :: LAER
    integer, intent(in) :: NQ
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR) :: &
         QM, QTT, QXT, QXX, QYT, QYY, QZT, QZZ, QXY, QXZ, QYZ
    !// Output
    real(r8), intent(out), dimension(LWEPAR) :: CNVLOST
    !// --------------------------------------------------------------------
    !// Locals
    real(r8) :: LQM(LPAR), LQTT(LPAR), DQM(LPAR,3)
    real(r8) :: DQTT(LPAR,3), DQMNT(9,LPAR,3), LQMNT(9,LPAR), PQMNT(9)
    real(r8) :: DFLUX,DTT,DMNT,PQM,PQTT,DFRACT
    integer :: K,L,M 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'QCNVW2_OSLO'
    !// --------------------------------------------------------------------

    CNVLOST(:) = 0._r8 !// Initialize what is lost

    if (SUM_S .lt. 1.e-30_r8) return  ! if no convection, exit

    !// LQ__ = air/tracer mass remaining as "ambient/environment"
    do L = 1, NQ
       LQM(L)     = QM(L)  !// Air mass
       LQTT(L)    = QTT(L) !// Tracer mass
       LQMNT(1,L) = QXT(L) !// Tracer moments
       LQMNT(2,L) = QYT(L) !//       :
       LQMNT(3,L) = QXX(L) !//       :
       LQMNT(4,L) = QYY(L) !//       :
       LQMNT(5,L) = QXY(L) !//       :
       LQMNT(6,L) = QZT(L) !//       :
       LQMNT(7,L) = QXZ(L) !//       :
       LQMNT(8,L) = QYZ(L) !//       :
       LQMNT(9,L) = QZZ(L) !//       :
    end do

    !// DQ__ = air/tracer parcels left behind by convection (1=N, 2=E, 3=D)
    do M = 1, 3
       do L = 1, NQ 
          DQM(L,M)   = 0._r8     !// Air mass
          DQTT(L,M)  = 0._r8     !// Tracer mass
          do K = 1,9
             DQMNT(K,L,M) = 0._r8 !// Tracer moments
          end do
       end do  ! L
    end do  ! M

!// ------------------------------------------------------------------------
!// #1 = N = non-entraining updrafts
!// ------------------------------------------------------------------------
    if (SUM_N .ge. 1.e-10_r8) then

       !// PQ__ = temporary storage of "plume" values in updraft/downdraft
       !// These are moved up in updrafts, and down in downdrafts
       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1, 9
          PQMNT(K) = 0._r8
       end do

       do  L = 1, NQ
         DFLUX = FLUX_N(L+1) - FLUX_N(L)          !// Net flux upwards

         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, bring ambient air into the updraft core
           if (DFLUX .gt. LQM(L)) DFLUX = LQM(L) !// Limit to grid box air mass
           DFRACT  = DFLUX / LQM(L)            !// Fraction of air moved
           PQM     = PQM + DFLUX               !// Increase plume air mass
           LQM(L)  = LQM(L) - DFLUX            !// Decrease ambient air mass
           DTT     = LQTT(L) * DFRACT          !// Tracer mass to be moved
           PQTT    = PQTT + DTT                !// Increase plume tracer mass
           LQTT(L) = LQTT(L) - DTT             !// Decrease ambient tracer mass
           do K = 1, 9
             DMNT = LQMNT(K,L) * DFRACT        !// Change in moments
             PQMNT(K) = PQMNT(K) + DMNT        !// Increase plume moments
             LQMNT(K,L) = LQMNT(K,L) - DMNT    !// Decrease ambient moments
           end do
         else if (DFLUX .lt. 0._r8) then
           !// DFLUX < 0, dump updraft core air into ambient
           if (PQM .gt. 0._r8)  then
             DFLUX     = -DFLUX                !// Change sign for simlicity
             if (DFLUX .gt. PQM)  DFLUX = PQM  !// Limit to plume air mass
             DFRACT    = DFLUX / PQM           !// Fraction of air moved
             DQM(L,1)  = DQM(L,1) + DFLUX      !// Air mass dumped to ambient
             PQM       = PQM - DFLUX           !// Decrease plume air mass
             DTT       = PQTT * DFRACT         !// Tracer mass to be dumped
             DQTT(L,1) = DQTT(L,1) + DTT       !// Tracer mass dumped to ambient
             PQTT      = PQTT - DTT            !// Decrease plume tracer mass
             do K = 1,9
               DMNT = PQMNT(K) * DFRACT           !// Change in moments
               PQMNT(K) = PQMNT(K) - DMNT         !// Decrease plume moments
               DQMNT(K,L,1) = DQMNT(K,L,1) + DMNT !// Increase ambient moments
             end do
           end if
         end if
         !// DFLUX = 0, fluxes are same, updraft simply passes through
       end do !// do  L = 1, NQ

    end if !// Non-entraining updrafts

!// ------------------------------------------------------------------------
!// #2 = E = entraining updrafts
!// ------------------------------------------------------------------------
    if (SUM_E .ge. 1.e-10_r8) then
       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1, 9
          PQMNT(K) = 0._r8
       end do

       do L = 1, NQ
         !// if entraining mass flux exists, first entrain air from ambient
         DFLUX = 0._r8
         if (ENT_U(L) .ge. 1.e-10_r8) then
           !// There is detrainment
           if (ENT_U(L) .gt. LQM(L))  then
             DFLUX  = LQM(L)                 !// Limit to grid box air mass
           else
             DFLUX  = ENT_U(L)
           end if
           DFRACT  = DFLUX / LQM(L)          !// Fraction of air to be entr.
           PQM     = PQM + DFLUX             !// Increase plume air mass
           LQM(L)  = LQM(L) - DFLUX          !// Decrease ambient air mass
           DTT     = LQTT(L) * DFRACT        !// Tracer mass to be entrained
           PQTT    = PQTT + DTT              !// Increase plume tracer mass
           LQTT(L) = LQTT(L) - DTT           !// Decrease ambient tracer mass
           do K = 1, 9
             DMNT = LQMNT(K,L) * DFRACT      !// Change in moments
             PQMNT(K) = PQMNT(K) + DMNT      !// Increase plume moments
             LQMNT(K,L) = LQMNT(K,L) - DMNT  !// Decrease ambient moments
           end do
         end if !// if (ENT_U(L) .ge. 1.e-10_r8) then

         !// Extra entrainment to balance fluxes
         DFLUX = FLUX_E(L+1) - (FLUX_E(L) + DFLUX)
         !// proceed as for non-entraining, assume full mixing in updraft core
         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, bring more ambient air into the updraft core
           if (DFLUX .gt. LQM(L)) DFLUX = LQM(L)  !// See the lines above
           DFRACT  = DFLUX / LQM(L)               !// for comments
           PQM     = PQM + DFLUX
           LQM(L)  = LQM(L) - DFLUX
           DTT     = LQTT(L) * DFRACT
           PQTT    = PQTT + DTT
           LQTT(L) = LQTT(L) - DTT
           do K = 1, 9
             DMNT = LQMNT(K,L) * DFRACT
             PQMNT(K) = PQMNT(K) + DMNT
             LQMNT(K,L) = LQMNT(K,L) - DMNT
           end do
         else if (DFLUX .lt. 0._r8) then
           !// DFLUX < 0, dump updraft core air into ambient
           !// (Detrainment is needed to balance fluxes)
           if (PQM .gt. 0._r8)  then
             DFLUX     = -DFLUX                !// Change sign for simlicity
             if (DFLUX .gt. PQM)  DFLUX = PQM  !// Limit to plume air mass
             DFRACT    = DFLUX / PQM           !// Fraction of air moved
             DQM(L,2)  = DQM(L,2) + DFLUX      !// Air mass dumped to ambient
             PQM       = PQM - DFLUX           !// Decrease plume air mass
             DTT       = PQTT * DFRACT         !// Tracer mass to be dumped
             DQTT(L,2) = DQTT(L,2) + DTT       !// Tracer mass dumped to ambient
             PQTT      = PQTT - DTT            !// Decrease plume tracer mass
             do K = 1, 9
               DMNT = PQMNT(K) * DFRACT           !// Change in moments
               PQMNT(K) = PQMNT(K) - DMNT         !// Decrease plume moments
               DQMNT(K,L,2) = DQMNT(K,L,2) + DMNT !// Increase ambient moments
             end do
           end if
         end if

         !// Oslo CTM2 removal method
         !// ------------------------
         !// PQTT is the mass traveling upwards; it may have grown (DFLUX>0)
         !// or it may have been dumped to the ambient (DFLUX<0), but after
         !// those processes, we remove the tracers in the plume according to
         !// the CTM3 wetloss fraction.
         !// The removed mass is diagnosed as PQTT*TWETLN(L) before the
         !// actual removal. After the downdrafts process (#3), CNVLOST is
         !// transported downwards and evaporation is calculated.
         if (L.le.NQ) then
           CNVLOST(L) = PQTT * TWETLN(L)
           PQTT = PQTT * (1._r8 - TWETLN(L))
         end if

         !// -DFLUX = 0, fluxes are same, updraft simply passes through
      end do !// do L = 1, NQ

    end if !// Entraining updrafts

!// ------------------------------------------------------------------------
!// Evaporation of falling rain done after #3
!// ------------------------------------------------------------------------

!// ------------------------------------------------------------------------
!// #3 = D = downdrafts
!// ------------------------------------------------------------------------
    if (abs(SUM_D) .ge. 1.e-10_r8) then

       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1, 9
          PQMNT(K) = 0._r8
       end do

       do L = NQ, 1, -1
      
         DFLUX = 0._r8
         if (ENT_D(L) .ge. 1.e-10_r8) then
           !// if downdraft mass flux exists, first entrain air from ambient
           if (ENT_D(L) .gt. LQM(L))  then
             DFLUX  = LQM(L)                 !// Limit to grid box air mass
           else
             DFLUX  = ENT_D(L)
           end if
           DFRACT  = DFLUX / LQM(L)          !// Fraction of mass to entrain
           PQM     = PQM + DFLUX             !// Increase plume air mass
           LQM(L)  = LQM(L) - DFLUX          !// Decrease ambient air mass
           DTT     = LQTT(L) * DFRACT        !// Tracer mass to be entrained
           PQTT    = PQTT + DTT              !// Increase plume tracer mass
           LQTT(L) = LQTT(L) - DTT           !// Decrease ambient tracer mass
           do K = 1, 9
             DMNT = LQMNT(K,L) * DFRACT      !// Change in moments
             PQMNT(K) = PQMNT(K) + DMNT      !// Increase plume moments
             LQMNT(K,L) = LQMNT(K,L) - DMNT  !// Decrease ambient moments
           end do
         end if !// if (ENT_D(L) .ge. 1.e-10_r8) then 

         !// Extra entrainment to balance flux
         DFLUX = FLUX_D(L+1) - (FLUX_D(L) + DFLUX)
         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, entrain air into downdraft
           if (DQM(L,1) .gt. 0._r8)  then
             !// pull first from air dumped at L by non-entraining updraft
             if (DFLUX .lt. DQM(L,1)) then
               DFRACT    = DFLUX / DQM(L,1)
               PQM       = PQM + DFLUX
               DQM(L,1)  = DQM(L,1) - DFLUX
               DTT       = DQTT(L,1) * DFRACT
               PQTT      = PQTT + DTT
               DQTT(L,1) = DQTT(L,1) - DTT
               do K = 1, 9
                 DMNT = DQMNT(K,L,1) * DFRACT
                 PQMNT(K) = PQMNT(K) + DMNT
                 DQMNT(K,L,1) = DQMNT(K,L,1) - DMNT
               end do
               cycle !// Done with this level (go to 34, which was continue)
             else
               !// Pull all dumped at L by non-entraining updraft
               PQM       = PQM + DQM(L,1)
               DFLUX     = DFLUX - DQM(L,1)
               DQM(L,1)  = 0._r8
               PQTT      = PQTT + DQTT(L,1) 
               DQTT(L,1) = 0._r8
               do K = 1, 9
                 PQMNT(K) = PQMNT(K) + DQMNT(K,L,1)  
                 DQMNT(K,L,1) = 0._r8
               end do
             end if
           end if

           if (DQM(L,2) .gt. 0._r8)  then
             !// if not enough, pull from air dumped at L by entraining updraft
             if (DFLUX .lt. DQM(L,2)) then
               DFRACT    = DFLUX / DQM(L,2)
               PQM       = PQM + DFLUX
               DQM(L,2)  = DQM(L,2) - DFLUX
               DTT       = DQTT(L,2) * DFRACT
               PQTT      = PQTT + DTT
               DQTT(L,2) = DQTT(L,2) - DTT
               do K = 1, 9
                 DMNT = DQMNT(K,L,2) * DFRACT
                 PQMNT(K) = PQMNT(K) + DMNT
                 DQMNT(K,L,2) = DQMNT(K,L,2) - DMNT
               end do
               cycle !// Done with this level (go to 34, which was continue)
             else
               PQM       = PQM + DQM(L,2)  
               DFLUX     = DFLUX - DQM(L,2)
               DQM(L,2)  = 0._r8
               PQTT      = PQTT + DQTT(L,2) 
               DQTT(L,2) = 0._r8
               do K = 1, 9
                 PQMNT(K) = PQMNT(K) + DQMNT(K,L,2)  
                 DQMNT(K,L,2) = 0._r8
               end do
             end if
           end if

           if (DFLUX .gt. 0._r8) then
             !// if still not enough, pull remainder from ambient
             DFRACT  = DFLUX / LQM(L)
             PQM     = PQM + DFLUX
             LQM(L)  = LQM(L) - DFLUX
             DTT     = LQTT(L) * DFRACT
             PQTT    = PQTT + DTT
             LQTT(L) = LQTT(L) - DTT 
             do K = 1, 9
               DMNT = LQMNT(K,L) * DFRACT
               PQMNT(K) = PQMNT(K) + DMNT
               LQMNT(K,L) = LQMNT(K,L) - DMNT
             end do
           end if

         else if (DFLUX .lt. 0._r8) then
           !// DFLUX < 0, dump downdraft core air into ambient
           if (PQM .lt. 1.e-10_r8) then
              write(6,'(a)') f90file//':'//subr// &
                   ': Downdrafts PQM < 1.e-10: should not be possible'
              stop 'STOP in '//subr
           end if
           DFLUX     = -DFLUX                 !// Change sign for simlicity
           if (DFLUX .gt. PQM)  DFLUX = PQM   !// Limit to plume air mass
           DFRACT    = DFLUX / PQM            !// Fraction of air moved
           PQM       = PQM - DFLUX            !// Decrease plume air mass
           DQM(L,3)  = DQM(L,3) + DFLUX       !// Air mass dumped to ambient
           DTT       = PQTT * DFRACT          !// Tracer mass to be dumped
           DQTT(L,3) = DQTT(L,3) + DTT        !// Tracer mass dumped to ambient
           PQTT      = PQTT - DTT             !// Decrease plume tracer mass
           do K = 1, 9
             DMNT = PQMNT(K) * DFRACT           !// Change in moments
             PQMNT(K) = PQMNT(K) - DMNT         !// Decrease plume moments
             DQMNT(K,L,3) = DQMNT(K,L,3) + DMNT !// Increase ambient moments
           end do
         end if
       end do !// do L = NQ, 1, -1

    end if !// Downdrafts
!// ------------------------------------------------------------------------
!// #4 = convective scavenging from updraft masses (DQTT(L,1:2)
!// ------------------------------------------------------------------------

!// UCI method of scavenging
!//    if (TWETLN .gt. 0._r8) then
!//      call SCAV_UPD (TWETLN,LAER,DQTT,DQMNT,ZOFL, NQ)
!//    end if

!// ------------------------------------------------------------------------
!// CTM3: Re-evaporation of falling rain
!// ------------------------------------------------------------------------
    !// Two options for evaporation of rain drops:
    !// 1. Easily dissolved species stay inside drops until all
    !//    is evaporated. Species are dumped in layer where this
    !//    happens.
    !// 2. Less soluble species can evaporate as long as there
    !//    is evaporation going on. Assume evaporation linearly
    !//    dependent on the fraction of rain evaporated, i.e.
    !//    calculated from rain out at bottom of grid box vs rain
    !//    in on top.

    !// A possible 3 would be a mix of the two, assuming that during
    !// evaporation some drops evaporate fully and others do not.
    !// Option 1 means that all drops get smaller when evaporation
    !// occurs.

    !// Skip adjusting moments. Assume evaporation is
    !// evenly distributed in grid box.

    !// CNVLOST is only accounting for what was removed during
    !// updraft. When this falls to the ground, there may be some evaporation
    !// and therefore CNVLOST must be adjusted for evaporation
    !// when going downwards

    if (HENRY_A .gt. 1._r8) then
       !// Option 1: Assume easily dissolved if Henry constant at
       !//           298K is greater than 1.
       PQTT = 0._r8           !// Elevator mass falling from top of plume
       do L = NQ, 1, -1
         !// Rain falling at level L
         PQTT = PQTT + CNVLOST(L)
         !// Go to next level if no rainout
         if (PQTT .eq. 0._r8) cycle

         if (CRAIN(L) .eq. 0._r8) then
           !// Dump all in ambient air
           LQTT(L) = LQTT(L) + PQTT
           !// CNVLOST accounts for what was lost above, summed up
           !// by PQTT down to this level L. To let CNVLOST account
           !// for what is evaporated, CNVLOST(L) should now become
           !// negative.
           !// All that was lost from top to L is evaporated here; PQTT
           !// already contains CNVLOST(L), and all is dumped.
           CNVLOST(L) = -PQTT
           !// CNVLOST(L:NQ) = 0._r8 !// must be wrong?: CNVLOST(L) - PQTT
           !// No tracer in plume, reset and keep looping downwards
           !// in case there is more rain below
           PQTT = 0._r8
         end if
       end do
    else
       !// Option 2
       PQTT = 0._r8           !// Elevator mass falling from top of plume
       do L = NQ, 1, -1
         !// Rain falling at level L
         PQTT = PQTT + CNVLOST(L)
         !// Go to next level if no rainout
         if (PQTT .eq. 0._r8) cycle

         if (L .lt. LWEPAR) then
           if (CRAIN(L) .lt. CRAIN(L+1)) then
             !// Rain out is smaller than rain in; find fractional
             !// mass to evaporate
             DTT = PQTT * (CRAIN(L+1) - CRAIN(L)) / CRAIN(L+1)
             !// New amount of tracer dissolved in falling rain
             PQTT = max(0._r8, PQTT - DTT)
             !// Update ambient tracer mass
             LQTT(L) = LQTT(L) + DTT
             !// Remove this part from CNVLOST, i.e. evaporate.
             CNVLOST(L) = CNVLOST(L) - DTT
             !// Skip adjusting moments. Assume evaporation is
             !// evenly distributed in grid box.
           end if
         end if
       end do
    end if !// if (TCHENA .gt. 1._r8) then

!// ------------------------------------------------------------------------
!// #5 = recombine up/downdraft parcels (DQ's) with remaining ambient (LQ's)
!// ------------------------------------------------------------------------
    do M = 1, 3
      do L = 1, NQ
        LQM(L) = LQM(L) + DQM(L,M)
        LQTT(L) = LQTT(L) + DQTT(L,M)
        !// note that only XY moments are added, all vertical moments
        !// in DQ__ are lost
        do K = 1, 5
          LQMNT(K,L) = LQMNT(K,L) + DQMNT(K,L,M)
        end do ! K
      end do ! L
    end do ! M
    do L = 1, NQ
      QM(L)    = LQM(L)
      QTT(L)   = LQTT(L)
      QXT(L)   = LQMNT(1,L)
      QYT(L)   = LQMNT(2,L)
      QXX(L)   = LQMNT(3,L)
      QYY(L)   = LQMNT(4,L)
      QXY(L)   = LQMNT(5,L)
      QZT(L)   = LQMNT(6,L)
      QXZ(L)   = LQMNT(7,L)
      QYZ(L)   = LQMNT(8,L)
      QZZ(L)   = LQMNT(9,L)
    end do

    !// --------------------------------------------------------------------
  end subroutine QCNVW2_OSLO
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine QCNVW2(FLUX_N,FLUX_E,FLUX_D,ENT_U,ENT_D,ZOFL, &
       QM,QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ, &
       SUM_N,SUM_E,SUM_D,SUM_S, TWETLN,LAER,NQ)
    !// --------------------------------------------------------------------
    !// New version drops the CFL criteria and just removes up to 99% of
    !// air mass in box, see ADJFLX.
    !// Also w* subsidence has ben moved out and combine with W-V-U sequence
    !// Processes:
    !//      1 = N = non-entraining updrafts
    !//      2 = E = entraining updrafts
    !//      3 = D = downdrafts
    !// Explicitly calculates transport in a column due to wet convection.
    !// Includes entraining updrafts/downdrafts, nonentraining updrafts.
    !// Calls scavenging routine SCAV_UPD for both updraft types.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, LWEPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in), dimension(LPAR+1) :: &
         FLUX_N, FLUX_E, FLUX_D, ENT_U, ENT_D, ZOFL
    real(r8), intent(in)  :: SUM_N, SUM_E, SUM_D, SUM_S
    real(r8), intent(in)  :: TWETLN
    logical, intent(in) :: LAER
    integer, intent(in) :: NQ
    !// Input/Output
    real(r8), intent(inout), dimension(LPAR) :: &
         QM, QTT, QXT, QXX, QYT, QYY, QZT, QZZ, QXY, QXZ, QYZ
    !// Locals
    real(r8) :: LQM(LPAR), LQTT(LPAR), DQM(LPAR,3)
    real(r8) :: DQTT(LPAR,3), DQMNT(9,LPAR,3), LQMNT(9,LPAR), PQMNT(9)
    real(r8) :: DFLUX,DTT,DMNT,PQM,PQTT,DFRACT
    integer :: K,L,M 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'QCNVW2'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr// &
         ': Should use QCNVW2_OSLO instead!'
    stop 'STOP in '//subr

    if (SUM_S .LT. 1.e-30_r8) return   ! if no convection, exit

    !// LQ__ = air/tracer mass remaining as "ambient/environment"
    do L = 1, NQ
       LQM(L)     = QM(L)
       LQTT(L)    = QTT(L)
       LQMNT(1,L) = QXT(L)
       LQMNT(2,L) = QYT(L)
       LQMNT(3,L) = QXX(L)
       LQMNT(4,L) = QYY(L)
       LQMNT(5,L) = QXY(L)
       LQMNT(6,L) = QZT(L)
       LQMNT(7,L) = QXZ(L)
       LQMNT(8,L) = QYZ(L)
       LQMNT(9,L) = QZZ(L)
    end do

    !// DQ__ = air/tracer parcels left behind by convection (1=N, 2=E, 3=D)
    do M = 1, 3
       do L = 1, NQ 
          DQM(L,M)   = 0._r8
          DQTT(L,M)  = 0._r8
          do K = 1,9
             DQMNT(K,L,M) = 0._r8
          end do
       end do  ! L
    end do  ! M

!// ------------------------------------------------------------------------
!// #1 = N = non-entraining updrafts
!// ------------------------------------------------------------------------
    if (SUM_N .ge. 1.e-10_r8) then

       !// PQ__ = temporary storage of "plume" values in updraft/downdraft
       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1,9
          PQMNT(K) = 0._r8
       end do

       do  L = 1,NQ
         DFLUX = FLUX_N(L+1) -FLUX_N(L)
         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, bring ambient air into the updraft core
           if (DFLUX .gt. LQM(L))  DFLUX = LQM(L)
           DFRACT  = DFLUX / LQM(L)
           PQM     = PQM + DFLUX
           LQM(L)  = LQM(L) - DFLUX
           DTT     = LQTT(L) * DFRACT
           PQTT    = PQTT + DTT
           LQTT(L) = LQTT(L) - DTT
           do K = 1,9
             DMNT = LQMNT(K,L) * DFRACT
             PQMNT(K) = PQMNT(K) + DMNT
             LQMNT(K,L) = LQMNT(K,L) - DMNT
           end do
         else if (DFLUX .lt. 0._r8) then
           !// DFLUX < 0, dump updraft core air into ambient
           if (PQM .gt. 0._r8)  then
             DFLUX     = -DFLUX
             if (DFLUX .gt. PQM)  DFLUX = PQM
             DFRACT    = DFLUX / PQM
             DQM(L,1)  = DQM(L,1) + DFLUX
             PQM       = PQM - DFLUX
             DTT       = PQTT * DFRACT
             DQTT(L,1) = DQTT(L,1) + DTT
             PQTT      = PQTT - DTT
             do K = 1,9
               DMNT = PQMNT(K) * DFRACT
               PQMNT(K) = PQMNT(K) - DMNT
               DQMNT(K,L,1) = DQMNT(K,L,1) + DMNT
             end do
           end if
         end if
         !// DFLUX = 0, fluxes are same, updraft simply passes through
       end do
    end if !// Non-entraining updrafts

!// ------------------------------------------------------------------------
!// #2 = E = entraining updrafts
!// ------------------------------------------------------------------------
    if (SUM_E .ge. 1.e-10_r8) then

       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1,9
          PQMNT(K) = 0._r8
       end do

       do L = 1,NQ    
         !// if entraining mass flux exists, first entrain air from ambient
         DFLUX = 0._r8
         if (ENT_U(L) .ge. 1.e-10_r8) then
           if (ENT_U(L) .gt. LQM(L))  then
             DFLUX  = LQM(L)
           else
             DFLUX  = ENT_U(L)
           end if
           DFRACT  = DFLUX / LQM(L)
           PQM     = PQM + DFLUX
           LQM(L)  = LQM(L) - DFLUX
           DTT     = LQTT(L) * DFRACT
           PQTT    = PQTT + DTT
           LQTT(L) = LQTT(L) - DTT
           do K = 1,9
             DMNT = LQMNT(K,L) * DFRACT
             PQMNT(K) = PQMNT(K) + DMNT
             LQMNT(K,L) = LQMNT(K,L) - DMNT
           end do
         end if

         DFLUX = FLUX_E(L+1) - (FLUX_E(L) + DFLUX)
         !// proceed as for non-entraining, assume full mixing in updraft core
         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, bring more ambient air into the updraft core
           if (DFLUX .gt. LQM(L))  DFLUX=LQM(L)
           DFRACT  = DFLUX / LQM(L)
           PQM     = PQM + DFLUX
           LQM(L)  = LQM(L) - DFLUX
           DTT     = LQTT(L) * DFRACT
           PQTT    = PQTT + DTT
           LQTT(L) = LQTT(L) - DTT
           do K = 1,9
             DMNT = LQMNT(K,L) * DFRACT
             PQMNT(K) = PQMNT(K) + DMNT
             LQMNT(K,L) = LQMNT(K,L) - DMNT
           end do
         else if (DFLUX .lt. 0._r8) then
           !// DFLUX < 0, dump updraft core air into ambient
           if (PQM .gt. 0._r8)  then
             DFLUX     = -DFLUX
             if (DFLUX .gt. PQM)  DFLUX=PQM
             DFRACT    = DFLUX / PQM
             DQM(L,2)  = DQM(L,2) + DFLUX
             PQM       = PQM - DFLUX
             DTT       = PQTT * DFRACT
             DQTT(L,2) = DQTT(L,2) + DTT
             PQTT      = PQTT - DTT
             do K = 1,9
               DMNT = PQMNT(K) * DFRACT
               PQMNT(K) = PQMNT(K) - DMNT
               DQMNT(K,L,2) = DQMNT(K,L,2) + DMNT
             end do
           end if
         end if
         !// DFLUX = 0, fluxes are same, updraft simply passes through
       end do
    end if !// Entraining updrafts

!// ------------------------------------------------------------------------
!// #3 = D = downdrafts
!// ------------------------------------------------------------------------
    if (abs(SUM_D) .ge. 1.e-10_r8) then

       PQM  = 0._r8
       PQTT = 0._r8
       do K = 1,9 
          PQMNT(K) = 0._r8
       end do

       do L = NQ, 1, -1
      
         DFLUX = 0._r8
         if (ENT_D(L) .ge. 1.e-10_r8)  then
           !// if downdraft mass flux exists, first entrain air from ambient
           if (ENT_D(L) .gt. LQM(L))  then
             DFLUX  = LQM(L)
           else
             DFLUX  = ENT_D(L)
           end if
           DFRACT  = DFLUX / LQM(L)
           PQM     = PQM + DFLUX
           LQM(L)  = LQM(L) - DFLUX
           DTT     = LQTT(L) * DFRACT
           PQTT    = PQTT + DTT
           LQTT(L) = LQTT(L) - DTT
           do K = 1, 9
             DMNT = LQMNT(K,L) * DFRACT
             PQMNT(K) = PQMNT(K) + DMNT
             LQMNT(K,L) = LQMNT(K,L) - DMNT
           end do
         end if

         DFLUX = FLUX_D(L+1) - (FLUX_D(L) + DFLUX)
         if (DFLUX .gt. 0._r8) then
           !// DFLUX > 0, entrain air into downdraft
           if (DQM(L,1) .gt. 0._r8)  then
             !// pull first from air dumped at L by non-entraining updraft
             if (DFLUX .lt. DQM(L,1)) then
               DFRACT    = DFLUX / DQM(L,1)
               PQM       = PQM + DFLUX
               DQM(L,1)  = DQM(L,1) - DFLUX
               DTT       = DQTT(L,1) * DFRACT
               PQTT      = PQTT + DTT
               DQTT(L,1) = DQTT(L,1) - DTT
               do K = 1, 9
                 DMNT = DQMNT(K,L,1) * DFRACT
                 PQMNT(K) = PQMNT(K) + DMNT
                 DQMNT(K,L,1) = DQMNT(K,L,1) - DMNT
               end do
               cycle !// Done with this level (go to 34, which was continue)
             else
               PQM       = PQM + DQM(L,1)
               DFLUX     = DFLUX - DQM(L,1)
               DQM(L,1)  = 0._r8
               PQTT      = PQTT + DQTT(L,1) 
               DQTT(L,1) = 0._r8
               do K = 1,9
                 PQMNT(K) = PQMNT(K) + DQMNT(K,L,1)  
                 DQMNT(K,L,1) = 0._r8
               end do
             end if
           end if
           if (DQM(L,2) .gt. 0._r8)  then
             !// if not enough, pull from air dumped at L by entraining updraft
             if (DFLUX .lt. DQM(L,2)) then
               DFRACT    = DFLUX / DQM(L,2)
               PQM       = PQM + DFLUX
               DQM(L,2)  = DQM(L,2) - DFLUX
               DTT       = DQTT(L,2) * DFRACT
               PQTT      = PQTT + DTT
               DQTT(L,2) = DQTT(L,2) - DTT
               do K = 1, 9
                 DMNT = DQMNT(K,L,2) * DFRACT
                 PQMNT(K) = PQMNT(K) + DMNT
                 DQMNT(K,L,2) = DQMNT(K,L,2) - DMNT
               end do
               cycle !// Done with this level (go to 34, which was continue)
             else
               PQM       = PQM + DQM(L,2)  
               DFLUX     = DFLUX - DQM(L,2)
               DQM(L,2)  = 0._r8
               PQTT      = PQTT + DQTT(L,2) 
               DQTT(L,2) = 0._r8
               do K = 1, 9
                 PQMNT(K) = PQMNT(K) + DQMNT(K,L,2)  
                 DQMNT(K,L,2) = 0._r8
               end do
             end if
           end if
           if (DFLUX .gt. 0._r8) then
             !// if still not enough, pull remainder from ambient
             DFRACT  = DFLUX / LQM(L)
             PQM     = PQM + DFLUX
             LQM(L)  = LQM(L) - DFLUX
             DTT     = LQTT(L) * DFRACT
             PQTT    = PQTT + DTT
             LQTT(L) = LQTT(L) - DTT 
             do K = 1, 9
               DMNT = LQMNT(K,L) * DFRACT
               PQMNT(K) = PQMNT(K) + DMNT
               LQMNT(K,L) = LQMNT(K,L) - DMNT
             end do
           end if
         else if (DFLUX.lt.0._r8) then
           !// DFLUX < 0, dump downdraft core air into ambient
           if (PQM .lt. 1.e-10_r8) then
              write(6,'(a)') f90file//':'//subr// &
                   ': Downdrafts PQM < 1.e-10: should not be possible'
              stop 'STOP in '//subr
           end if
           DFLUX     = -DFLUX
           if (DFLUX .gt. PQM)  DFLUX = PQM
           DFRACT    = DFLUX / PQM
           PQM       = PQM - DFLUX
           DQM(L,3)  = DQM(L,3) + DFLUX
           DTT       = PQTT * DFRACT
           DQTT(L,3) = DQTT(L,3) + DTT
           PQTT      = PQTT - DTT
           do K = 1, 9
             DMNT = PQMNT(K) * DFRACT
             PQMNT(K) = PQMNT(K) - DMNT
             DQMNT(K,L,3) = DQMNT(K,L,3) + DMNT
           end do
         end if
       end do !// do L = NQ, 1, -1

    end if !// Downdrafts

!// ------------------------------------------------------------------------
!// #4 = convective scavenging from updraft masses (DQTT(L,1:2)
!// ------------------------------------------------------------------------
    if (TWETLN .gt. 0._r8) then
       call SCAV_UPD (TWETLN,LAER,DQTT,DQMNT,ZOFL, NQ)
    end if

!// ------------------------------------------------------------------------
!// #5 = recombine up/downdraft parcels (DQ's) with remaining ambient (LQ's)
!// ------------------------------------------------------------------------
    do M = 1, 3
       do L = 1, NQ
          LQM(L) = LQM(L) + DQM(L,M)
          LQTT(L) = LQTT(L) + DQTT(L,M)
          !// note that only XY moments are added, all vertical moments
          !// in DQ__ are lost
          do K = 1, 5
             LQMNT(K,L) = LQMNT(K,L) + DQMNT(K,L,M)
          end do ! K
       end do ! L
    end do ! M
    do L = 1,NQ
       QM(L)    = LQM(L)
       QTT(L)   = LQTT(L)
       QXT(L)   = LQMNT(1,L)
       QYT(L)   = LQMNT(2,L)
       QXX(L)   = LQMNT(3,L)
       QYY(L)   = LQMNT(4,L)
       QXY(L)   = LQMNT(5,L)
       QZT(L)   = LQMNT(6,L)
       QXZ(L)   = LQMNT(7,L)
       QYZ(L)   = LQMNT(8,L)
       QZZ(L)   = LQMNT(9,L)
    end do

    !// --------------------------------------------------------------------
  end subroutine QCNVW2
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine SCAV_UPD (FSCAV,LAER,DQTT,DQMNT,ZOFL,NQ)
    !// --------------------------------------------------------------------
    !// scavenge from convective updrafts:  DQTT(L,1:2)
    !// remove a fixed fraction of tracer: FSCAV 
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)    :: FSCAV
    logical, intent(in)   :: LAER
    real(r8), intent(in)    :: ZOFL(LPAR+1)
    integer, intent(in)   :: NQ
    !// Input/Output
    real(r8), intent(inout) :: DQTT(LPAR,3), DQMNT(9,LPAR,3)
    !// Locals
    integer :: K,L,M
    real(r8) :: FSC, FRAC
    !// --------------------------------------------------------------------

    !// scavenge from convective updrafts:  DQTT(L,1:2)
    !// note that aerosols also use FSCAV, but reduced by 2 below 2.6 km
    if (FSCAV .gt. 0._r8) then
       do L = 1, NQ
          FSC = FSCAV
          if (LAER) then
             if ((ZOFL(L+1)-ZOFL(1)).lt.2600._r8) then
                FSC = 0.5_r8*FSCAV
             end if
          end if
          FRAC = 1._r8 - FSC
          do M = 1, 2
             !// reduce tracer mass and moments in updrafts
             DQTT(L,M) = DQTT(L,M) * FRAC
             do K= 1, 5
                DQMNT(K,L,M) = DQMNT(K,L,M) * FRAC
             end do
          end do
       end do
    end if

    !// --------------------------------------------------------------------
  end subroutine SCAV_UPD
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module convection
!//=========================================================================

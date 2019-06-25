!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Gravitational settling of aerosols.
!//=========================================================================
module fallingaerosols
  !// ----------------------------------------------------------------------
  !// MODULE: fallingaerosols
  !// DESCRIPTION: Routine for transporting aerosols downwards due to
  !//              gravitation, calculated using the second order
  !//              moments (SOM) scheme of Prather.
  !//
  !// Amund Sovde, February - March 2013
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8, rMom
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Falling velocities
  !// Fritz Kasten, Journal of Applied Meteorology, vol. 7, p. 944, 1968
  integer, parameter :: NRAD = 8, NLEV = 54
  real(r8) :: k68_rad(NRAD)
  real(r8) :: k68_z(NLEV),k68_v(NRAD,NLEV)
  !// ----------------------------------------------------------------------
  save
  private
  public aerosolsettling, readkasten1968
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine aerosolsettling(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
       AIRB,BTEM,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)
    !// --------------------------------------------------------------------
    !//
    !// Amund Sovde, February 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, NPAR, LSULPHUR, &
         LDUST, LSALT, LNITRATE, LBCOC, LSOA, NPAR_DUST, &
         NPAR_BC, NPAR_OM
    use cmn_ctm, only: ETAA, ETAB
    use cmn_oslo, only: SCAV_MAP_DRY, trsp_idx, dustbinsradii
    use bcoc_oslo, only: bc_trsp_idx, BC_IDS, OM_IDS, bcsnow_diag_dd
    use dust_oslo, only: dust_trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in)  ::  AIRB, BTEM
    real(r8), intent(in)  :: DTGAMA
    integer, intent(in) :: NDAY, NMET, NOPS, NSUB, MP

    real(r8), dimension(LPAR,IDBLK,JDBLK) :: GAMAAER
    real(r8), dimension(IDBLK,JDBLK,NPAR) :: BTTonGRND
    real(r8) :: radius, rho, pres
    real(r8) :: radiusL(LPAR), rhoL(LPAR)
    integer :: TRID, N, II, JJ, L
    integer :: NDEPPED(NPAR)
    !// --------------------------------------------------------------------

    !// Initialize
    BTTonGRND(:,:,:) = 0._r8
    NDEPPED(:) = 0

    !// Sulphur particles (must be removed from stratloss)
    if (LSULPHUR) then
       !// Sediment SO4
       TRID = trsp_idx(73)
       if (TRID .gt. 0) then
          NDEPPED(TRID) = 1 !// Flag the deposited component
          radius = 0.08e-6_r8 ! [m]
          call getGAMAAER(radius, DTGAMA, AIRB, GAMAAER,MP)
          call GRAV_SETN(TRID,BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
               BTTonGRND(:,:,TRID),AIRB,GAMAAER,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)
       end if
       !// Sediment MSA
       TRID = trsp_idx(75)
       if (TRID .gt. 0) then
          NDEPPED(TRID) = 1 !// Flag the deposited component
          radius = 0.08e-6_r8 ! [m]
          call getGAMAAER(radius, DTGAMA, AIRB, GAMAAER,MP)
          call GRAV_SETN(TRID,BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
               BTTonGRND(:,:,TRID),AIRB,GAMAAER,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)
       end if
    end if

    !// DUST (must be removed from dust treatment, in dstdpsdry.F90)
    if (LDUST) then
       do N = 1, NPAR_DUST
          !// Transport id
          TRID = dust_trsp_idx(N)
          NDEPPED(TRID) = 1 !// Flag the deposited component
          !// Get radius for this bin (dstpsd.F90)
          radiusL(:) = dustbinsradii(N)
          !// Set density [kg/m3]
          rhoL(:) = 2.6e3_r8
          !// Get falling speed
          !call oc_getGAMAAER(radius, DTGAMA, AIRB, GAMAAER,MP)
          call getGAMAAER_NS(radiusL, rhoL, DTGAMA, BTEM, AIRB, GAMAAER,MP,1)
          !// Do transport
          call GRAV_SETN(TRID,BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
               BTTonGRND(:,:,TRID),AIRB,GAMAAER,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)

       end do
    end if

    !// SALT (must be removed from salt treatment)

    !// BCOC (must be removed from bcoc treatment)
    if (LBCOC) then
!       do N = 1, NPAR_BCOC
!          !// BCOC_IDS gives ID number
!          !// Transport id
!          TRID = bc_trsp_idx(N)
!          NDEPPED(TRID) = 1 !// Flag the deposited component
!          !// Get radius for this bin (use effective radius from FJX)
!          if (BCOC_IDS(N).le.243) then
!             radiusL(:) = 0.08e-6_r8
!          else
!             radiusL(:) = 0.039e-6_r8
!          end if
!          !// Increase radius as a parameterization to allow for coagulation.
!          !do L=1,LPAR
!          !   pres = 0.5_r8 * (ETAA(L)+ETAA(L+1) + 1000._r8*(ETAB(L)+ETAB(L+1)))
!          !   !// Double radius at 200hPa and upwards.
!          !   if (pres .lt. 200._r8) then
!          !      radiusL(L) = radiusL(L) * 2._r8
!          !   else
!          !      radiusL(L) = radiusL(L) * (1._r8 + (1000._r8-pres)/800._r8)
!          !   end if
!          !end do
!          !// Set density [kg/m3]
!          if (BCOC_IDS(N).eq.242 .or. BCOC_IDS(N).eq.243 .or. &
!               BCOC_IDS(N).eq.246 .or. BCOC_IDS(N).eq.247) then
!             rhoL(:) = 1.35e3_r8
!          else
!             rhoL(:) = 1.e3_r8
!          end if
!          !// Get falling speed
!          !call oc_getGAMAAER(radius, DTGAMA, AIRB, GAMAAER,MP)
!          call getGAMAAER_NS(radiusL, rhoL, DTGAMA, BTEM, AIRB, GAMAAER,MP,0)
!          !// Do transport
!          call GRAV_SETN(TRID,BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
!               BTTonGRND(:,:,TRID),AIRB,GAMAAER,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)
!
!          !// Diagnose dry dep on snow
!          call bcsnow_diag_dd(BTTonGRND(:,:,TRID),BCOC_IDS(N),MP)
!       end do
    end if

    !// SOA (must be removed from stratloss)

    !// NITRATE  (must be removed from nitrate treatment)


    !// Add BTTonGND to drydep diagnoses
    do N = 1, NPAR
       if (NDEPPED(N).eq.0) cycle !// Not deposited
       do JJ = 1, JDBLK
          do II = 1, IDBLK
             SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) &
                                        + BTTonGRND(II,JJ,N)
          end do
       end do
    end do


    !// --------------------------------------------------------------------
  end subroutine aerosolsettling
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine readkasten1968()
    !// --------------------------------------------------------------------
    !// Read in the Kasten (1968) data.
    !// --------------------------------------------------------------------
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=80) :: cline,filename
    integer :: N,L, fnr, io_err
    !// --------------------------------------------------------------------

    !// Open file
    filename = 'Indata_CTM3/kasten1968_fallingvelocities.dat'
    fnr = get_free_fileid()
    open(fnr,file=filename,status='old',form='formatted',iostat=io_err)

    !// Check for table file
    if (io_err .ne. 0) then
       print*,'No such file: '//trim(filename)
       stop
    end if

    !// Read header
    do N=1,5
       read(fnr,'(a)') cline
    end do
    read(fnr,'(4x,8(4x,f6.3))') k68_rad(:)

    !// Read height and velocities
    do L=1, NLEV
       read(fnr,'(f5.0,8(1x,e9.4))') k68_z(L), k68_v(:,L)
    end do

    !// Close file
    close(fnr)

    !// --------------------------------------------------------------------
  end subroutine readkasten1968
  !// ----------------------------------------------------------------------






  !// ----------------------------------------------------------------------
  subroutine getGAMAAER(radius, DTGRAV, AIRB, GAMAAER,MP)
    !// --------------------------------------------------------------------
    !// The vertical velocity of particles are found based on their radius.
    !// Assumes fixed radius for the whole IJ-block, but can in principle
    !// be set up using an array containing 3D (IJ-block) radii of aerosols.
    !//
    !// Amund Sovde, February 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIE, MPBLKIB, MPBLKJE, MPBLKJB, CFLLIM
    use cmn_met, only: ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in)  :: radius, DTGRAV
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK)
    !// Output
    real(r8), intent(out) :: GAMAAER(LPAR,IDBLK,JDBLK)
    !// Locals
    real(r8) :: QM(LPAR), QMO(LPAR), QU(LPAR), ZH(LPAR), DZ(LPAR)
    real(r8) :: wvel,whmin,whmax,rfrac,hfrac,mflux
    integer :: N, I, J, II, JJ, L, LRmin, LRmax, Hmin, Hmax
    !// --------------------------------------------------------------------

    !// Find closest radii (k68_rad(LRmin) <= radius < k68_rad(LRmax))
    !// in Kasten data
    LRmin = 1
    do N = 1, NRAD
       if (radius .ge. k68_rad(N)) LRmin = N
    end do
    LRmax = N+1

    !// If radius is smaller than k68_rad(1), set equal to LRmin
    if (radius .lt. k68_rad(1)) LRmax = LRmin
    !// If radius is bigger than k68_rad(NRAD), set equal to LRmin
    if (radius .ge. k68_rad(NRAD)) LRmax = LRmin
    !// Sanity check
    if (LRmax .gt. NRAD) then
       print*,'shouldnt happen, LRmin,LRmax,radius',LRmin,LRmax,radius
       stop
    end if

    !// Fractional distance from LRmax (linear interpolation)
    !// Should consider logarithmic interpolation
    if (LRmax .eq. LRmin) then
       rfrac = 0._r8
    else
       rfrac = (radius - k68_rad(LRmin))/(k68_rad(LRmax)-k68_rad(LRmin))
    end if


    !// Loop on columns
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Air mass
        QM(:) = AIRB(:,II,JJ)
        !// Assume ZOFLE to be similar to standard atmosphere for now.
        !// Future: Calculate height in std atm
        !ZH(:) = 1.e-3_r8 * ZOFLE(:,I,J) !// km
        ZH(:) = 0.5e-3_r8 * (ZOFLE(1:LPAR,I,J) + ZOFLE(2:LPAR+1,I,J))
        !// Level thickness [m]
        DZ(:) = ZOFLE(2:LPAR+1,I,J) - ZOFLE(1:LPAR,I,J)

        !// Will include velocity at ground, to calculate dry deposition.
        !// If you want separate routine for that, this must be changed here
        !// and in transport routine below. GAMAAER(1,II,JJ) = 0._r8


        do L = LPAR, 1, -1
          !// Pick level from k68
          Hmin = 1
          do N = 1, NLEV
            if (ZH(L) .ge. k68_z(N)) then
               Hmin = N
               exit
            end if
          end do
          Hmax = Hmin + 1
          if (Hmax .gt. NLEV) Hmax = Hmin

          !// get velocity
          if (Hmax .eq. Hmin) then
            !// below k68_z(1) or above k68_z(NLEV); no interpolation of
            !// height, only radius
            wvel = k68_v(Hmin,LRmin)*(1._r8-rfrac) + k68_v(Hmin,LRmax)*rfrac
          else
            !// Interpolate in height and radius
            whmin = k68_v(Hmin,LRmin)*(1._r8-rfrac) + k68_v(Hmin,LRmax)*rfrac
            whmax = k68_v(Hmax,LRmax)*(1._r8-rfrac) + k68_v(Hmax,LRmax)*rfrac
            !// Find velocity of grid box center
            if (ZH(L) .gt. k68_z(Hmin)) then
               hfrac = (ZH(L) - k68_z(Hmin)) / (k68_z(Hmax)-k68_z(Hmin))
            else
               hfrac = 0._r8
            end if
            wvel = whmin * (1._r8-hfrac) + whmax*hfrac
          end if

          !// Find air mass flux for this velocity. Airmass is not
          !// to be transported, but VECT3 uses the airmass flux as
          !// velocity.
          !//   Unit: 1.e-2_r8[m/cm] * w[cm/s] / dz[m] * mass[kg] = kg/s
          mflux = 1.e-2_r8 * wvel / DZ(L) * QM(L)

          !// Limit flux to CFLLIM of gridbox mass:
          !// Have to change sign to negative for falling particles.
          if (mflux*DTGRAV .gt. CFLLIM*QM(L)) then
             QU(L) = -CFLLIM*QM(L)/DTGRAV
          else
             QU(L) = -mflux
          end if

        end do

        !// Store the flux [kg/s]
        GAMAAER(:,II,JJ) = QU(:)

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine getGAMAAER
  !// ----------------------------------------------------------------------







  !// ----------------------------------------------------------------------
  subroutine GRAV_SETN(NTR,BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
       BTTonGRND,AIRB,GAMAAER,DTGAMA,NDAY,NMET,NOPS,NSUB,MP)
    !// --------------------------------------------------------------------
    !// Routine for transporting a tracer downwards by a settling velocity.
    !// Based on DYN2W_OC.
    !//
    !// Will calculate mass deposited on surface/ground!
    !//
    !// Amund Sovde, February 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIE, MPBLKIB, MPBLKJE, MPBLKJB, LMTSOM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    real(r8), dimension(IDBLK,JDBLK), intent(out) :: BTTonGRND
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: GAMAAER
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in)  ::  AIRB
    real(r8), intent(in)  :: DTGAMA
    integer, intent(in) :: NTR, NDAY, NMET, NOPS, NSUB, MP

    real(r8), dimension(LPAR+3) :: & ! add ground + possible 1 for collapse
         QM,          &  ! air mass in box at start,
         QU,          &  ! air mass flux moved [I-1]->[I] in adv. step
         QTT,         &  ! tracer mass in box [I]
         QXT,QYT,QZT, &  ! 1st moments of tracer in U, V, W direction
         QXX,QYY,QZZ, &  ! 2nd moments of tracer in U, V, W direction
         QXY,QYZ,QXZ, &  ! cross-moments of tracer
         Q0F,Q1F         ! computed tracer flux from [I-1] to [I]
    integer :: NQ        ! length of vector pipe for advection (assumed cyclic)
    integer :: NSTEP     ! #multi-steps needed for local CFL, ret by QVECT3
    integer :: L,II,JJ
    !// --------------------------------------------------------------------

    !// This routine does transport on only one species.
    !// In case of odd numbers of layers, when NPARDIV =/= 1, the pipe
    !// needs to be increased by 1. When transporting several tracers,
    !// this should rather be done as in DYN2W_OC, where NPARDIV tracers
    !// in the column are stacked to form a long pipe of even size.
    if (mod(LPAR,2) .eq. 0) then
      NQ = LPAR+2 !// Must add 2 to account for surface
    else
      NQ = LPAR+1 !// Must add 1 to account for surface+collapse
    end if


    !// begin major loop over J's and I's
    do JJ = 1,JDBLK
      do II = 1,IDBLK

        !// Need to ensure that GAMAAER(1,I,J)=0 so that pipes do NOT connect
        QU(1) = 0._r8 !// Ground
        QU(LPAR+1)  = 0._r8 !// Top when collapsed
        QU(LPAR+2)  = 0._r8 !// Top when not collapsed
        QU(LPAR+3)  = 0._r8 !// Top when not collapsed
        do L = 1,LPAR
          QU(L+1)  = GAMAAER(L,II,JJ)*DTGAMA
        end do

        if (NTR .gt. NPAR) then
          print*,'*** OC_GRAV_SETN: NTR > NTM; very wrong: STOPPING'
          stop
        end if

        !// QVECT3 does not allow air mass to be zero after transport.
        !// To avoid this, we set it to 1 at ground
        QM(1)  = 1._r8
        QTT(1) = 0._r8
        QXT(1) = 0._r8
        QYT(1) = 0._r8
        QZT(1) = 0._r8
        QXX(1) = 0._r8
        QYY(1) = 0._r8
        QZZ(1) = 0._r8
        QXY(1) = 0._r8
        QXZ(1) = 0._r8
        QYZ(1) = 0._r8
        !// In case of odd LM, add 1 in LM+1 (will not affect results; GAMA=0)
        !// In case of even LM, add 1 in LM+2 (will not affect results; GAMA=0)
        !// And Q-s are of size NQ+1:
        do L=1,3
           QTT(LPAR+L) = 0._r8
           QXT(LPAR+L) = 0._r8
           QYT(LPAR+L) = 0._r8
           QZT(LPAR+L) = 0._r8
           QXX(LPAR+L) = 0._r8
           QYY(LPAR+L) = 0._r8
           QZZ(LPAR+L) = 0._r8
           QXY(LPAR+L) = 0._r8
           QXZ(LPAR+L) = 0._r8
           QYZ(LPAR+L) = 0._r8
        end do
        QM(LPAR+1)  = 1._r8
        QM(LPAR+2)  = 1._r8
        QM(LPAR+3)  = 0._r8

        !// Fetch air and tracer
        do L = 1,LPAR
          !// must re-init the airmass for each tracer
          QM(L+1)  = AIRB(L,II,JJ)
          QTT(L+1) = BTT(L,NTR,II,JJ)
          QXT(L+1) = BXT(L,NTR,II,JJ)
          QYT(L+1) = BYT(L,NTR,II,JJ)
          QZT(L+1) = BZT(L,NTR,II,JJ)
          QXX(L+1) = BXX(L,NTR,II,JJ)
          QYY(L+1) = BYY(L,NTR,II,JJ)
          QZZ(L+1) = BZZ(L,NTR,II,JJ)
          QXY(L+1) = BXY(L,NTR,II,JJ)
          QXZ(L+1) = BXZ(L,NTR,II,JJ)
          QYZ(L+1) = BYZ(L,NTR,II,JJ)
        end do

        !// ----------------------------------------------------------
        call QVECT3(QM,QU,NQ,LMTSOM, NSTEP, &
              QTT,QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY, Q0F,Q1F)
        !WSTEP = WSTEP + NSTEP
        !// ----------------------------------------------------------

        do L = 1,LPAR
          !// Fetch new masses. QTT(1) is mass deposited on ground.
          BTT(L,NTR,II,JJ) = QTT(L+1)
          BXT(L,NTR,II,JJ) = QXT(L+1)
          BYT(L,NTR,II,JJ) = QYT(L+1)
          BZT(L,NTR,II,JJ) = QZT(L+1)
          BXX(L,NTR,II,JJ) = QXX(L+1)
          BYY(L,NTR,II,JJ) = QYY(L+1)
          BZZ(L,NTR,II,JJ) = QZZ(L+1)
          BXY(L,NTR,II,JJ) = QXY(L+1)
          BXZ(L,NTR,II,JJ) = QXZ(L+1)
          BYZ(L,NTR,II,JJ) = QYZ(L+1)
        end do

        !// Diagnose mass deposited on surface
        BTTonGRND(II,JJ) = QTT(1)

        !// Do *NOT* update air mass, it is not to be transported
        !// when aerosols fall due to gravity.

      end do       !// I loop
    end do       !// J loop

    !// --------------------------------------------------------------------
  end subroutine GRAV_SETN
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine fallspeed(LPAR,diam,rhoL,pres,tem,wvel)
    !// --------------------------------------------------------------------
    !// Purpose: Given size and density of particle, and column
    !// thermodynamic profile, compute terminal fall speed.
    !// Based on DEAD/DUST code.
    !//
    !// diam: Particle diameter [m]
    !//       Can in principle vary with height.
    !// rho:  Particle density [kg m-3]
    !//       Can in principle vary with height.
    !// pres: Pressure [Pa]
    !// tem:  Temperature [K]
    !// wvel: Settling velocity [m s-1]
    !//
    !// Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR, R_UNIV, G0, CPI
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LPAR
    real(r8), intent(in) :: &
         rhoL(LPAR), & ! I [kg m-3] Particle density
         diam(LPAR), & ! I [m] Particle diameter
         pres(LPAR), & ! I [Pa] Pressure 
         tem(LPAR)     ! I [K] Temperature
    !// Output
    real(r8), intent(out) :: wvel(LPAR)   ! O [m s-1] Settling velocity

    !// Locals
    integer :: L     ! [idx] Counting index
    real(r8) :: &
         mfp_atm, &  ! [m] Mean free path of air
         slp_crc, &  ! [frc] Slip correction factor
         vsc_dyn_atm ! [kg m-1 s-1] Dynamic viscosity of air

    !// Parameters
    real(r8), parameter :: &
         mmw_dry_air = M_AIR*1e-3_r8, & ! [kg mol-1] Mean mol. weight of dry air
         gas_cst_unv = R_UNIV           ! [J mol-1 K-1] Universal gas constant
    !// --------------------------------------------------------------------

    !// Size-independent thermokinetic properties
    do L = 1, LPAR
       !// [kg m-1 s-1] Dynamic viscosity of air
       !// [kg m-1 s-1] RoY94  p. 102
       vsc_dyn_atm = 1.72e-5_r8 * ((tem(L)/273._r8)**1.5_r8) &
                             * 393._r8/(tem(L) + 120._r8)
       !// [m] Mean free path of air SeP97 p. 455
       mfp_atm = 2._r8*vsc_dyn_atm / &
            (pres(L) * sqrt(8._r8*mmw_dry_air / (cpi * gas_cst_unv * tem(L))))

       !// [frc] Slip correction factor SeP97 p. 464
       slp_crc = 1._r8 + 2._r8*mfp_atm*(1.257_r8 &
                   + 0.4_r8*exp(-1.1_r8 * diam(L) / (2._r8*mfp_atm))) / diam(L)
       !// [m s-1] Stokes settling velocity SeP97 p. 466
       wvel(L) = (1._r8/18._r8) * diam(L)**2 * rhoL(L) &
                              * g0 * slp_crc / vsc_dyn_atm
       !// [m s-1] Corrected Stokes settling velocity (depends on
       !// particle size etc.)
       !// Skipping this for now; in CTM3 DUST code, this is about 0.995,
       !// which we approximate to 1.
       !wvel(L) = wvel(L) * stk_crc
    end do

    !// --------------------------------------------------------------------
  end subroutine fallspeed
  !// ----------------------------------------------------------------------
 


  !// ----------------------------------------------------------------------
  subroutine turbfallspeed(diam, rho, pres, tem, wvel, &
       wnd_frc, air_rho, hgt_mdp, hgt_zpd, mno_lng, rgh_mmn, tvel)
    !// --------------------------------------------------------------------
    !// Purpose: Given size and density of particle, and column
    !// thermodynamic profile, compute 
    !// Based on DEAD/DUST code.
    !//
    !// diam: Mass weighted particle diameter [m]
    !// rho:  Particle density [kg m-3]
    !// pres: Pressure [Pa]
    !// tem:  Temperature [K]
    !// wvel: Calculated settling velocity at surface (will be used to
    !//       calculate turbulent settling)
    !// wnd_frc: friction velocity [m s-1] (i.e. USTAR)
    !// air_rho: air density [kg m-3]
    !// hgt_mdp: Midlayer height above surface [m]
    !// hgt_zpd: Zero plane displacement height [m]
    !// mno_lng: Monin-Obukhov length [m]
    !// rgh_mmn: Roughness length momentum [m]
    !// tvel: Turbulent settling velocity [m s-1]
    !//
    !// Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR, R_UNIV, G0, CPI, BOLTZMANN
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: &
         rho, &     ! [kg m-3] Particle density
         diam, &    ! [m] Particle diameter
         pres, &    ! [Pa] Pressure at surface
         tem, &     ! [K] Temperature at surface
         wvel, &    ! [m s-1] Gravitational settling velocity
         wnd_frc, & ! [m s-1] Friction velocity
         air_rho,&  ! [kg m-3] Air density
         hgt_mdp, & ! [m] Midlayer height above surface
         hgt_zpd, & ! [m] Zero plane displacement height
         mno_lng, & ! [m] Monin-Obukhov length
         rgh_mmn    ! [m] Roughness length momentum

    !// Output
    real(r8), intent(out) :: tvel   ! O [m s-1] Turbulent settling velocity


    !// Locals
    integer :: L              ! [idx] Counting index
    real(r8) :: &
         mfp_atm, &     ! [m] Mean free path of air
         slp_crc, &     ! [frc] Slip correction factor
         vsc_dyn_atm, & ! [kg m-1 s-1] Dynamic viscosity of air
         vsc_knm_atm, & ! [m2 s-1] Kinematic viscosity of air
         vlc_grv, &     ! [m s-1] Stokes settling velocity 
         stk_nbr, &     ! [frc] Stokes number SeP97 p. 965
         dff_aer, &     ! [m2 s-1] Brownian diffusivity of particle
         shm_nbr, &     ! [frc] Schmidt number
         shm_nbr_xpn, & ! [frc] Surface-dependent exponent
         rss_lmn, &     ! [s m-1] Laminar resistance
         rss_trb, &     ! [s m-1] Resistance to turbulent deposition
         tmp            ! Temporary variable

    real(r8) :: &
         eta_sqr_rlm, &   ! Eta squared at roughness height
         eta_sqr_ctm, &   ! Eta squared at model layer height
         eta_rlm, &       ! Eta at roughness height
         eta_ctm, &       ! Eta at model layer height
         nmr_rlm, &       ! Numerator
         dnm_ctm, &       ! Denominator
         tmp4, &          ! Correction to neutral atmosphere factor
         tmp5, &          ! Neutral atmosphere factor
         mno_prm_rlm, &   ! [frc] Monin-Obukhov parameter at roughness height
         mno_prm_ctm, &   ! [frc] Monin-Obukhov parameter at model layer height
         rss_aer          ! [s m-1] Aerodynamical resistance

    !// Parameters
    real(r8), parameter :: &
         mmw_dry_air = M_AIR*1e-3_r8, & ![kg mol-1] Mean mol. weight dry air
         gas_cst_unv = R_UNIV, &        ![J mol-1 K-1] Universal gas constant
         cst_Boltzmann = BOLTZMANN, &   ![J K-1] Boltzmann's constant
         shm_nbr_xpn_lnd = -2._r8/3._r8, & ![frc] Exponent for aerosol-diffusion
                                           !dependence on Schmidt number
                                           !over land
         shm_nbr_xpn_ocn = -0.5_r8, &  ![frc] Exponent for aerosol-diffusion
                                       !dependence on Schmidt number over ocean
         cst_von_krm = 0.4_r8          ![fraction] Von Karman constant
    !// --------------------------------------------------------------------

    !// Size-independent thermokinetic properties
    !// [kg m-1 s-1] Dynamic viscosity of air
    !// [kg m-1 s-1] RoY94  p. 102
    vsc_dyn_atm = 1.72e-5_r8 * ((tem/273._r8)**1.5_r8) &
                           * 393._r8/(tem + 120._r8)
    !// [m] Mean free path of air SeP97 p. 455
    mfp_atm = 2._r8*vsc_dyn_atm / &
             (pres * sqrt(8._r8*mmw_dry_air / (cpi * gas_cst_unv * tem)))
    !// [m2 s-1] Kinematic viscosity of air
    vsc_knm_atm = vsc_dyn_atm/air_rho

    !// [frc] Slip correction factor SeP97 p. 464
    slp_crc = 1._r8 + 2._r8*mfp_atm*(1.257_r8 &
                   + 0.4_r8*exp(-1.1_r8 * diam / (2._r8*mfp_atm))) / diam

    !// [m s-1] Stokes settling velocity SeP97 p. 466
    vlc_grv = (1._r8/18._r8) * diam**2 * rho * g0 * slp_crc / vsc_dyn_atm
    !// [m s-1] Corrected Stokes settling velocity (depends on particle
    !// size etc.)
    !// Skipping this for now; in CTM3 DUST code, this is about 0.995,
    !// which we approximate to 1.
    !vlc_grv = vlc_grv * stk_crc

    !// [frc] Stokes number SeP97 p. 965
    stk_nbr = vlc_grv * wnd_frc**2 / (g0 * vsc_knm_atm)
    !// [m2 s-1] Brownian diffusivity of particle SeP97 p. 474
    dff_aer = cst_Boltzmann * tem * slp_crc &
               / (3._r8 * cpi * vsc_dyn_atm * diam)
    !// [frc] Schmidt number SeP97 p. 972
    shm_nbr = vsc_knm_atm / dff_aer
    !// [frc] Surface-dependent exponent for aerosol-diffusion dependence
    !// on Schmidt number          
    shm_nbr_xpn = shm_nbr_xpn_lnd
    !// fxm: Turning this on dramatically reduces deposition velocity in low
    !// wind regimes.
    !// Schmidt number exponent is -2/3 over solid surfaces and
    !// -1/2 over liquid surfaces SlS80 p. 1014
    !// if (oro(i) == 0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else
    !// shm_nbr_xpn=shm_nbr_xpn_lnd ! [frc] Surface-dependent exponent for
    !// aerosol-diffusion dependence on Schmidt number

    !// Calculate laminar resistance
    tmp = shm_nbr**shm_nbr_xpn + 10._r8**(-3._r8/stk_nbr)
    rss_lmn = 1._r8/(tmp * wnd_frc) !// [s m-1] SeP97 p. 972, 965
 
    !// Calculate aerodynamical resistance
    !// Compute stability parameter: Maximum value is 1 because stability
    !// correction function is valid only for zeta < 1, e.g.,
    !// Bon96 p. 52, Bru82 p. 71, SeP97 p. 963
    mno_prm_rlm = min(rgh_mmn/mno_lng, 1._r8) ! [frc]
    mno_prm_ctm = min((hgt_mdp - hgt_zpd) / mno_lng,1._r8) ! [frc]
    if (mno_lng < 0._r8) then
       !// Difference between unstable corrections
       eta_sqr_rlm = sqrt(1._r8 - 16._r8 * mno_prm_rlm)
       eta_sqr_ctm = sqrt(1._r8 - 16._r8 * mno_prm_ctm)
       eta_rlm = sqrt(eta_sqr_rlm)
       eta_ctm = sqrt(eta_sqr_ctm)
       nmr_rlm = (eta_sqr_rlm + 1._r8) * (eta_rlm + 1._r8)**2
       dnm_ctm = (eta_sqr_ctm + 1._r8) * (eta_ctm + 1._r8)**2
       tmp4 = log(nmr_rlm/dnm_ctm) + 2._r8*(atan(eta_ctm)-atan(eta_rlm)) ! [frc]
    else
       ! Difference between stable corrections
       tmp4 = 5._r8 * (mno_prm_ctm - mno_prm_rlm) ! [frc]
    end if
    tmp5 = log((hgt_mdp - hgt_zpd)/rgh_mmn)
    rss_aer = (tmp4 + tmp5) / (cst_von_krm * wnd_frc) ! [s m-1] Bon96 p. 54
 

    !// [s m-1] Resistance to turbulent deposition
    rss_trb = rss_aer + rss_lmn + rss_aer * rss_lmn * wvel
    !// [m s-1] Turbulent deposition velocity
    tvel = 1._r8 / rss_trb

    !// --------------------------------------------------------------------
  end subroutine turbfallspeed
  !// ----------------------------------------------------------------------
 



  !// ----------------------------------------------------------------------
  subroutine getGAMAAER_NS(radiusL, rhoL, DTGRAV, BTEM, AIRB, GAMAAER,MP,TURB)
    !// --------------------------------------------------------------------
    !// The vertical velocity of particles are found based on their radius.
    !// Assumes fixed radius for the whole IJ-block, but can in principle
    !// be set up using an array containing 3D (IJ-block) radii of aerosols.
    !//
    !// Amund Sovde, June 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIE, MPBLKIB, MPBLKJE, MPBLKJB, CFLLIM, &
         ETAA, ETAB, PLAND, JMON
    use cmn_met, only: ZOFLE, P, USTR, SHF, SFT
    use cmn_sfc, only: ZOI, NVGPAR, landSurfTypeFrac, ZPDVT_C3
    use utilities, only: moninobukhov_length
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP, TURB
    real(r8), intent(in)  :: radiusL(LPAR),rhoL(LPAR), DTGRAV
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK), BTEM(LPAR,IDBLK,JDBLK)
    !// Output
    real(r8), intent(out) :: GAMAAER(LPAR,IDBLK,JDBLK)
    !// Locals
    real(r8) :: QM(LPAR), QMO(LPAR), QU(LPAR), PRES(LPAR), TEM(LPAR), DZ(LPAR)
    real(r8) :: wvel(LPAR),tvel,mflux,adens, mo_len,zo,SFMU,SFNU, hgt_zpd
    integer :: N, I, J, II, JJ, L, LRmin, LRmax, Hmin, Hmax
    !// --------------------------------------------------------------------

    !// Loop on columns
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Air mass
        QM(:) = AIRB(:,II,JJ)
        !// Temperature
        TEM(:) = BTEM(:,II,JJ)
        !// Pressure [Pa]
        PRES(:) = 0.5e2_r8 * (ETAA(1:LPAR) + P(I,J)*ETAB(1:LPAR) &
                         + ETAA(2:LPAR+1) + P(I,J)*ETAB(2:LPAR+1))
        !// Level thickness [m]
        DZ(:) = ZOFLE(2:LPAR+1,I,J) - ZOFLE(1:LPAR,I,J)

        !// Calculate fall speed [m/s]
        call fallspeed(LPAR,2._r8*radiusL,rhoL,pres,tem,wvel)

        !// Calculate additional settling due to turbulent mixing in
        !// lowermost layer. Use surface rather than midpoint:
        !adens = QM(1) / (AREAXY(I,J)*DZ(1)) !// [kg/m3]
        adens = 1._r8 * P(I,J) / (SFT(I,J) * 287._r8) !// [kg/m3]

        !// Calculate M-O length (differs from DUST code)
        MO_LEN = moninobukhov_length(adens,SFT(I,J),USTR(I,J),SHF(I,J))

        !// Surface roughness [m] and displacement height [m]
        !// Similar to DEAD calculation, but with CTM3 metdata as input.
        if (PLAND(I,J) .ge. 0.25_r8) then
           ZO = ZOI(I,J,JMON)
           hgt_zpd = 0._r8
        else
           !// Correct water surf roughness for wind/waves:
           !// Absol.visc. 6.2d-8*T (lin fit:-30C to +40C):
           SFMU  = 6.2e-8_r8*SFT(I,J)
           !// Kinematic visc(nu) = mu/density  (m*m/s)
           SFNU  = SFMU / adens
           ZO = min(0.135_r8 * SFNU/USTR(I,J) &
                       + 1.83e-3_r8 * USTR(I,J)**2, 2.e-3_r8)
           hgt_zpd = 0._r8
           do L = 1, NVGPAR
              hgt_zpd = hgt_zpd + landSurfTypeFrac(L,I,J) * ZPDVT_C3(L)
           end do
        end if
          

        !// Include turbulent fall speed at surface?
        if (TURB .eq. 1) then
           call turbfallspeed(2._r8*radiusL(1), rhoL(1),pres(1), tem(1), &
                wvel(1), USTR(I,J), adens, dZ(1)*0.5_r8, hgt_zpd, mo_len, &
                zo, tvel)
           !// Add turbulent velocity to total velocity
           wvel(1) = wvel(1) + tvel
        end if

        do L = LPAR, 1, -1
           !// Find air mass flux for this velocity. Airmass is not
           !// to be transported, but VECT3 uses the airmass flux as
           !// velocity.
           !//   Unit:  w[m/s] / dz[m] * mass[kg] = kg/s
           mflux = wvel(L) / DZ(L) * QM(L)

           !// Limit flux to CFLLIM of gridbox mass:
           !// Have to change sign to negative for falling particles.
           if (mflux*DTGRAV .gt. CFLLIM*QM(L)) then
              QU(L) = -CFLLIM*QM(L)/DTGRAV
           else
              QU(L) = -mflux
           end if

        end do !// do L = LPAR, 1, -1

        !// Store the flux [kg/s]
        GAMAAER(:,II,JJ) = QU(:)

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine getGAMAAER_NS
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module fallingaerosols
!//=========================================================================

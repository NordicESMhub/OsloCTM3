!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Sea salt module.
!//=========================================================================
module seasalt
  !// ----------------------------------------------------------------------
  !// MODULE: seasalt
  !// DESCRIPTION: Module for SEA SALT package.
  !//
  !// Most of the SEA SALT module is kept private.
  !// Master routines are global.
  !//
  !// Contains SALT variables and the following routines:
  !//   subroutine seasalt_init
  !//   subroutine get_seasaltscheme
  !//   subroutine get_rsalt80um
  !//   subroutine seasalt_getflux
  !//   subroutine seasalt_master
  !//   subroutine dry2wet
  !//   subroutine growth_factor
  !//   subroutine falling
  !//   subroutine drydeppart
  !//   subroutine addtogether
  !//   subroutine saltbdg2file
  !//   subroutine seasaltbudget_wdep
  !//   subroutine vdep_n_coarse
  !//   subroutine seasalt_emis
  !//
  !// References (by Alf Grini)
  !//   J. W. Fitzgerald, Approximation formula for the equilibrium
  !//   size of an aerosol particle as a function of its dry 
  !//   size and composition and the ambient relative humidity, 1975
  !//   J. Appl. Met,14,pp 1044-1049
  !//
  !// Ole Amund Sovde, June 2015, September 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK, NPAR_SALT
  use cmn_parameters, only: G0, CPI, R_UNIV, TK_0C, BOLTZMANN, cp_air
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Parameters
  !// Max number of salt tracers
  integer, parameter :: max_SALT_comps = NPAR_SALT
  !// Number of salt bins
  integer, parameter :: Nsaltbins      = NPAR_SALT
  !// Lowest diam dry bin limit (m)
  real(r8), parameter  :: Dsaltdrymin = 0.03e-6_r8
  !// Highest diam dry bin limit (m)
  real(r8), parameter  :: Dsaltdrymax = 25.0e-6_r8
  !// Density [kg/m3]
  real(r8), parameter  :: rhosaltdry  = 2200._r8

  logical, parameter :: &
       LRHFIXED = .false. !// Fixed RH or not

  !// Scheme for sea salt production
  !// 1: Monahan et al (Oceanic Whitecaps, 1986, Reidel)
  !//    and Smith (QJRMS, 1993, doi:10.1002/qj.49711951211)
  !// 2: Mårtensson et al (JGR, 2003, doi:10.1029/2002JD002263)
  !// 3: Gantt et al (GMD, 2015, doi:10.5194/gmd-8-619-2015)
  !// 4: Witek et al (JGR, 2016, doi:10.1002/2015JD023726)
  integer, parameter :: SeaSaltScheme = 4
  character(len=80), dimension(4), parameter :: SeaSaltSchemeName = &
       (/'Monahan et al 1986 & Smith 1993', &
         'Mårtensson et al 2003', &
         'Gantt et al 2015 (Gahn et al 2003 & SST adjustment of Jaegle 2011)', &
         'Witek et al 2016 (Sofiev et al 2011 & SST adjustment of Jaegle 2011)' /)

  !// Variables
  !// List of all transport numbers used for SALT components
  integer   :: salt_trsp_idx(Nsaltbins)
  !// Total numbers of components in use
  integer   :: salt_idx_tot

  real(r8), dimension(Nsaltbins+1) :: &
       Dsaltdry, & !// Dry diameters of bins (m)
       Dsalt80, &  !// Diameters at 80% RH (m)
       rsalt80um   !// Radius at 80%RH in um
  real(r8), dimension(Nsaltbins) :: &
       saltbinsradii !// Mean radii at 80%RH (m)
  real(r8), dimension(LPAR) :: &
       RHfixed

  !// Emission array: Calculated from wind, i.e. every time meteorology is
  !// updated. Unit is kg/s.
  real(r8), dimension(NPAR_SALT, IDBLK, JDBLK, MPBLK) :: seasalt_flux


  !// Budget variables
  !// ----------------------------------------------------------------------
  !// Flag to turn on 3D diagnostics. If false, only total budgets
  !// will be put out to budget file.
  logical, parameter :: LSALTDIAG3D = .false.

  !// 2D fields for burden and budget processes; all salt tracers
  real(r8), dimension(Nsaltbins,ipar,jpar) :: &
       salt_production, &
       salt_drydep, &
       salt_lscv, &
       salt_cscv, &
       salt_burden

  real(r8) :: saltbudgetacctime   !// Accumulated time span (s)
  !// File for averages
  character(len=80), parameter :: saltbdgfile='ctm3sltbudget.nc'
  integer :: saltaverages_written !// Number of averages written to file

  !// Local debug
  logical, parameter :: Ldebug_salt = .TRUE.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'seasalt.f90'
  !// ----------------------------------------------------------------------
  save    !// All variables are to be saved.
  private !// Keep most of the module private
  !// And specify what will be public
  public max_SALT_comps, salt_trsp_idx, Nsaltbins, &
       seasalt_init, seasalt_master, saltbdg2file, seasaltbudget_lscav, &
       saltbinsradii, seasalt_emis, get_rsalt80um, get_seasaltscheme, &
       seasalt_getflux, emissions_seasalt_total
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine seasalt_init()
    !// --------------------------------------------------------------------
    !// Initialize SALT simulations. Only figures out the transport numbers
    !// and indices for SALT components.
    !//
    !// Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, LNITRATE
    use cmn_chem, only: TNAME
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Local variables
    integer :: N
    character(len=10) :: tracername !// Name of tracer

    real(r8) :: &
         logbin1, &    !// Logarithm of lowest bin length (log(m))
         logbinmax, &  !// Logarithm of highest bin length (log(m))
         totlogstep, & !// Logarithm of whole bin length
         logstep, &    !// Logarithm of one bin length
         logbin        !// Logarithm of arbitrary bin length
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='seasalt_init'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Sea SALT initialisation'

    !// Initialize number of SALT components
    salt_idx_tot = 0

    !// Get transport numbers for salt tracers
    do N = 1, NPAR
       tracername = TNAME(N)
       if (tracername(1:4) .eq. 'SALT') then
          salt_idx_tot = salt_idx_tot + 1
          salt_trsp_idx(salt_idx_tot) = N
       end if
    end do

    write(6,'(a,i5)') 'Sea salt scheme index: ',SeaSaltScheme
    write(6,'(a)')    'Sea salt scheme name:  '// &
         trim(SeaSaltSchemeName(SeaSaltScheme))

    !// Check if you read any SALT tracers at all
    if (salt_idx_tot .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': YOU HAVE NO COMPONENTS STARTING WITH SALT!'
       stop
    else if (salt_idx_tot .gt. Nsaltbins) then
       write(6,'(a,2i7)') f90file//':'//subr// &
            ': Too many SALT tracers found ', salt_idx_tot, Nsaltbins
       stop
    else
       write(6,'(a,i7)') f90file//':'//subr// &
            ': Number of SALT tracers/bins found ',salt_idx_tot
    end if


    !// Set salt lowest and highest diameter
    Dsaltdry(1) = Dsaltdrymin
    Dsaltdry(Nsaltbins+1) = Dsaltdrymax

    !// logarithmically even-spaced Diameters
    logbin1   = log(Dsaltdry(1))
    logbinmax = log(Dsaltdry(Nsaltbins+1))
    totlogstep= logbinmax - logbin1
    logstep   = totlogstep / real(Nsaltbins, r8)

    logbin = logbin1
    do N = 2,Nsaltbins+1
       logbin = logbin+logstep
       Dsaltdry(N) = exp(logbin)
    end do


    !// Fixed RH at 80% if LRHfixed is true 
    RHfixed(:)  = 0.8_r8
    !// Diameters at 80 % RH
    Dsalt80(:)  = 2._r8 * Dsaltdry(:)  !Fitzgerald 1975
    !// Radius in um at 80 % RH
    rsalt80um(:)= 0.5e6_r8 * Dsalt80(:)
    do N = 1, Nsaltbins
       saltbinsradii(N) = 0.5_r8 * (Dsalt80(N) + Dsalt80(N+1))
    end do


    write(6,'(a)') '  Bin    Dsalt80    Dsaltdry'
    do N = 1, Nsaltbins
       write(6,'(2x,i2,5x,es9.3,2x,es9.3)') N,Dsaltdry(N),Dsaltdry(N)
    end do


    !// Budget stuff
    saltbudgetacctime = 0._r8
    saltaverages_written = 0

    salt_production(:,:,:) = 0._r8
    salt_drydep(:,:,:)     = 0._r8
    salt_lscv(:,:,:)       = 0._r8
    salt_cscv(:,:,:)       = 0._r8
    salt_burden(:,:,:)     = 0._r8

    !// Check nitrate tracers
    if (LNITRATE) then
       if (trsp_idx(63) .le. 0 .or. trsp_idx(65) .le. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Nitrate is included, but tracers are missing!'
          write(6,'(a,i5)') ' trsp_idx(63) NH4coarse',trsp_idx(63)
          write(6,'(a,i5)') ' trsp_idx(65) NO3coarse',trsp_idx(65)
          stop
       end if
    end if !// if (LNITRATE) then

    write(6,'(a)') f90file//':'//subr//': Sea SALT initialised'

    !// --------------------------------------------------------------------
  end subroutine seasalt_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_rsalt80um(rsalt80,NB)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NB
    real(r8), intent(out) :: rsalt80(NB+1)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_rsalt80um'
    !// --------------------------------------------------------------------
    if (NB .ne. NPAR_SALT) then
       write(6,'(a,2i7)') f90file//':'//subr// &
            ': NB.ne.NPAR_SALT: ',NB, NPAR_SALT
       stop
    end if
    rsalt80(1:NB+1) = rsalt80um(1:NB+1)
    !// --------------------------------------------------------------------
  end subroutine get_rsalt80um
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_seasaltscheme(NN)
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(out) :: NN
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_seasaltscheme'
    !// --------------------------------------------------------------------
    NN = SeaSaltScheme
    !// --------------------------------------------------------------------
  end subroutine get_seasaltscheme
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_getflux(ssProd,NB,II,JJ,MP)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NB,II,JJ,MP
    real(r8), intent(out) :: ssProd(NB)
    !// --------------------------------------------------------------------
    ssProd(1:NB) = seasalt_flux(1:NB,II,JJ,MP)
    !// --------------------------------------------------------------------
  end subroutine seasalt_getflux
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_master(AIRB, BTEM, BDV, BTT, QHUM, ZOI, JMON, DTCHM, MP)
    !// --------------------------------------------------------------------
    !// Integrate Oslo Chemistry package for SEA SALT aerosols using the QSSA
    !// integrator.
    !//
    !// Will loop through IJ-block and integrate vertically.
    !//
    !// Note that sea salt is calculated for all (I,J), but with zero
    !// production over land and over sea ice.
    !//
    !// Changes from Oslo CTM2 to Oslo CTM3:
    !//   - Treating a column instead of longitude and column
    !//   - Roughness length not read from file; it is available in the CTM
    !//     core (ZOI)
    !//
    !// To CTM3: Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, LNITRATE
    use cmn_ctm, only: PLAND, AREAXY, &
       MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_met, only: SFT, SHF, USTR, UMS, VMS, CI, ZOFLE
    use cmn_sfc, only: LSMASK
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// In
    integer, intent(in)   :: MP, JMON
    real(r8), intent(in)  :: DTCHM
    real(r8), intent(in)  :: BTEM(LPAR,IDBLK,JDBLK), &
                             AIRB(LPAR,IDBLK,JDBLK), &
                             BDV(LPAR,IDBLK,JDBLK), &
                             QHUM(IPAR,JPAR,LPAR), &
                             ZOI(IPAR,JPAR,12)

    !// In/Out
    real(r8), intent(inout) :: BTT(LPAR,NPAR,IDBLK,JDBLK)

    !// Locals
    real(r8), dimension(LPAR, Nsaltbins) :: SALTT   !// Local salt array
    real(r8)  :: Dsaltwet(LPAR,Nsaltbins+1), &
               rhosaltwet(LPAR,Nsaltbins), &
               sloss(LPAR,Nsaltbins), &
               prodseasalt(Nsaltbins)
    !// Local meteorology
    real(r8), dimension(LPAR) :: T_COL, AIR_COL, Q_COL, ALPHA_COL
    real(r8) :: AXY, Z0_LOC, mass, massflux, DZ, losscoarse
    real(r8) :: bcol(Nsaltbins)

    real(r8), dimension(Nsaltbins) :: DRYDEP

    !// Indices
    integer :: I,J,II,JJ, M
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='seasalt_master'
    !// --------------------------------------------------------------------



    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Save production for M7, otherwise do the following:

          !// Set up column variables
          !// --------------------------------------------------------------
          !// Salt masses
          do M = 1, salt_idx_tot
             SALTT(:,M) = BTT(:,salt_trsp_idx(M),II,JJ)
          end do

          !// Temperature
          T_COL(:) = BTEM(:,II,JJ)

          !// Air
          AIR_COL(:) = AIRB(:,II,JJ)

          !// Get inverse volume
          Q_COL(:) = QHUM(I,J,:)

          !// Get inverse volume
          ALPHA_COL(:) = 1._r8 / BDV(:,II,JJ)

          !// Grid box height
          DZ = (ZOFLE(2,I,J) - ZOFLE(1,I,J))

          !// Grid box area
          AXY = AREAXY(I,J)

          !// Z0 can change in this code. Initialize with CTM3 ZOI
          Z0_LOC = ZOI(I,J,JMON)


          !// Treat salt processes
          !// --------------------------------------------------------------
          !// Dry to wet conversion. Output: Dsaltwet, rhosaltwet
          call dry2wet(LRHfixed, RHfixed, Q_COL, AIR_COL, ALPHA_COL, T_COL, &
                       Dsaltwet, rhosaltwet)


          !// Get production of particles at 80 % RH (S=0.8)
          !// Production is now calculated in a separate routine, called
          !// each from pmain (update_emis_ij) each time meteorology is
          !// updated. The emission array (seasalt_flux) can then be used
          !// by other routines as well.
          prodseasalt(:) = seasalt_flux(:,II,JJ,MP)


          !// Get the falling velocities for layer 2 and above
          call FALLING( T_COL, ALPHA_COL, AXY, Dsaltwet, rhosaltwet, sloss)

          !// Get the roughness lengths for sea surface
          call z0_get(USTR(I,J), PLAND(I,J), T_COL(1), Z0_LOC)

          !// Get the dry deposition velocities, correcting the falling
          !// velocity for turbulent mix-out and laminar resistance in the
          !// lowest layer
          call DRYDEPPART(T_COL(1), AIR_COL(1), ALPHA_COL(1), AXY, &
               Dsaltwet(1,:), rhosaltwet(1,:), SHF(I,J), Z0_LOC, USTR(I,J), &
               SFT(I,J), sloss(1,:))


          !// Final steps for salt
          !// --------------------------------------------------------------
          !// Adding together.
          !// BUDGET: Also sums up production and loss for diagnostics. At the
          !// bottom of oc_salt_master the time steps are summed up, so an
          !// average of the diagnostics can be generated.
          call addtogether(SALTT, prodseasalt, sloss, AXY, dtchm, DRYDEP, I,J)

          !// Save for global CTM scavenging diagnostic [kg]
          do M = 1, salt_idx_tot
             call add2drydepdiag(II,JJ,MP,salt_trsp_idx(m),DRYDEP(m))
          end do

          !// Add up for salt drydep budget, all tracers [kg]
          salt_drydep(:,I,J) = salt_drydep(:,I,J) + DRYDEP(:)
          !// Add up for salt production budget, all tracers [kg]
          salt_production(:,I,J) = salt_production(:,I,J) &
                                   + prodseasalt(:) * DTCHM


          !// Put the masses back in BTT
          do M = 1, salt_idx_tot
             BTT(:,salt_trsp_idx(M),II,JJ) = SALTT(:,M)
          end do

          !// BUDGET: Accumulate column burden for average [kg]
          do M = 1, salt_idx_tot
             bcol(M) = sum(SALTT(:,M))
          end do
          !// Save burden [kg*s]
          !// Burden is weighted with DTCHM to allow for a possible
          !// varying DTCHM (between NOPS, not varying between MP-blocks).
          salt_burden(:,I,J) = salt_burden(:,I,J) + bcol(:) * DTCHM
          

          !// Get a mass-weighted loss term for the coarse mode aerosols
          massflux = 0._r8  !kg/sec (loss out of box)
          mass = 0._r8      !kg     (concentration in box)     
          do M = 1, salt_idx_tot
             !//Note "L" is replaced by "1" since we are using first layer
             massflux = massflux + SALTT(1,M) * sloss(1,M)
             mass     = mass + SALTT(1,M)
          end do !// do M = 1,salt_idx_tot
          if (mass .gt. 0._r8) then
             losscoarse = massflux / mass
          else
             losscoarse = 1.e-5_r8  !Set low loss if no mass
          end if !// Check on mass.gt.0.0

          !// Save losscoarse for NITRATE species
          if (LNITRATE) call vdep_n_coarse(losscoarse,DZ,I,J)

       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)

    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// Save total time (seconds) for average
    !// Only do this for one IJ-block, to count the time only once.
    if (MP .eq. 1) saltbudgetacctime = saltbudgetacctime + DTCHM


    !// --------------------------------------------------------------------
  end subroutine seasalt_master
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DRY2WET( Lrhknown, RH, Q, AIR, ZDV, T, Dpartwet, rhopartwet )
    !// --------------------------------------------------------------------
    !// PURPOSE: 
    !// Calculate the density and diameter of wet sea salt particles 
    !//
    !// Theory: J. W. Fitzgerald, Approximation formula for the equilibrium
    !// size of an aerosol particle as a function of its dry 
    !// size and composition and the ambient relative humidity, 1975
    !// J. Appl. Met,14,pp 1044-1049
    !//
    !// Author for Oslo CTM2: Alf Grini
    !//
    !// Routine strides a little for L,N loops, but it is probably ok due to
    !// small array sizes.
    !//
    !// To CTM3: Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// INPUT
    real(r8), intent(in) :: Q(LPAR)   !// Water vapor mixing ratio (kg/kg)
    real(r8), intent(in) :: T(LPAR)   !// Temperature (K)
    real(r8), intent(in) :: AIR(LPAR) !// mass of air in grid cell (kg)
    real(r8), intent(in) :: ZDV(LPAR) !// Inverse volume (m-3)
    real(r8), intent(in) :: RH(LPAR)  !// Relative humidity (0-1) (if known)
    logical, intent(in)  :: LRHknown  !// is rel.humidity already known (?)

    !// OUTPUT
    real(r8), intent(out) :: Dpartwet(LPAR,Nsaltbins+1) !// Diameter of wet particles (m)
    real(r8), intent(out) :: rhopartwet(LPAR,Nsaltbins) !// Density of wet particles (kg/m3)

    !// LOCAL
    integer ::  I,L,N        !Counting indexes for longitudes, layers and bins
    real(r8) :: Dbindry      !radius of dry aerosol (m)
    real(r8) :: Dbinwet      !radius of wet aerosol (m)
    real(r8) :: Vpartwet     !Volume of wet aerosol (m3)
    real(r8) :: Vpartdry     !Volume of dry aerosol
    real(r8) :: Vwater       !Volume of water
    real(r8) :: Newmass      ! ??
    real(r8) :: Pvapor       !Partial pressure of water vapor (Pa)
    real(r8) :: P0           !Saturation partial pressure of water vapor (Pa)
    real(r8) :: TC           !Temperature in degrees Celcius
    real(r8) :: weightvapor  !Mass of water vapor
    real(r8) :: S            !Relative humidity (0-1)
    logical :: LSUL          !Are we calculating on sulphate (?) 
    real(r8) :: rwet         !wet radius (m)
    real(r8) :: rdry         !Dry radium (m)

    !// LOCAL PARAMETERS
    logical, parameter :: Ldebug_slt_d2w=.FALSE.  ! Write debug to screen
    real(r8), parameter :: &
         a0 = 6.1078_r8, &      !
         a1 = 4.4365e-1_r8, &   !Coefficients used to calculate
         a2 = 1.4289e-2_r8, &   !Saturation partial pressure of water
         a3 = 2.6507e-4_r8, &   !at a certain temperature
         a4 = 3.0312e-6_r8, &   !Seinfeld & Pandis, atmos. chem. and physics 
         a5 = 2.0341e-8_r8, &   !pp 781, table 15.2
         a6 = 6.1368e-11_r8, &  !
         rhowater = 1000._r8    !Density of water (kg/m3)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='DRY2WET'
    !// --------------------------------------------------------------------
    
    !// Initializing, not calculating sulphate here
    LSUL = .FALSE.           

    !// Loop through column
    DO L = 1, LPAR

       !// Check on if we already know RH in all grid cells
       IF (LRHKNOWN) THEN
          S = RH(L)
       ELSE
          !// WE DON'T KNOW WHAT RH WE ARE CALCULATING FOR
          TC = T(L) - TK_0C
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: TC',TC
          !// Saturation water pressure               
          P0 = a0 + a1*TC + a2*TC*TC &
               + a3*TC**3 + a4*TC**4 &
               + a5*TC**5 + a6*TC**6
          P0 = P0 * 100._r8  !mbar-->*Pa*(mbar/bar)-->Pa

          !// weight of H2O (in kg)
          weightvapor = Q(L)*AIR(L)


          !// Partial pressure of H2O using ideal gas (pV=nRT)
          pvapor = (weightvapor * 1.e3_r8/18.02_r8) * R_UNIV * T(L) * ZDV(L)
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Pvapor',pvapor

          !// Relative humidity
          S = Pvapor/P0      !RELATIVE HUMIDITY
          IF (S.gt.1._r8) then
             S = 1._r8
          ELSE IF (S.lt.0._r8) then
             S = 0.01_r8
          END IF
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Relative humidity ',I,L,S

          !// Don't do calculation id T (degrees C) < -50 
          IF (TC.lt.-50._r8) S = 0._r8

       END IF  !Check on known relative humidity



       !// Loop on all bin limits
       DO N = 1, Nsaltbins + 1

          !// Set dry RADIUS
          rdry = 0.5_r8*Dsaltdry(N)

          !// Calculating growth factor (new diameters)
          call growth_factor(S, rdry, rwet, LSUL)

          !// Set wet DIAMETER
          Dpartwet(L,N) = 2._r8*rwet

       END DO  !Loop on bins



       !// Loop on all bins
       DO N = 1, Nsaltbins

          !// Representative diameter of "wet" bin
          Dbinwet = 0.5_r8*(Dpartwet(L,N) + Dpartwet(L,N+1))
          IF(Ldebug_slt_d2w)write(6,*)subr//':debug: Dbinwet',N,Dbinwet

          !// Representative Diameter of "dry" bin
          Dbindry = 0.5_r8*(Dsaltdry(N) + Dsaltdry(N+1))
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Dbindry',N,Dbindry

          !// Volume of wet particle (4/3*pi*r^3)
          Vpartwet = 4._r8 / 3._r8 * CPI * (0.5_r8 * Dbinwet)**3
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Vpartwet',Vpartwet

          !// Volume of dry particle
          Vpartdry = 4._r8 / 3._r8 * CPI * (0.5_r8 * Dbindry)**3
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Vpartdry',Vpartdry

          !// Volume of pure water, ignoring changement of volume due to mixing
          Vwater = Vpartwet - Vpartdry
          IF (Ldebug_slt_d2w) write(6,*)subr//':debug: Vwater',Vwater

          !// density of wet particles, ignoring changement of volume due
          !// to mixing
          rhopartwet(L,N) = (Vpartdry*rhosaltdry + Vwater*rhowater) / Vpartwet

          IF (Ldebug_slt_d2w) THEN
             write(6,*)subr//':debug: RH, rhopartwet',S,rhopartwet(L,N)
             write(6,*)subr//':debug: rhosaltdry',rhosaltdry
          END IF
               
       END DO !// DO N = 1, Nsaltbins

    END DO !// DO L=1,LPAR


    !// --------------------------------------------------------------------
  end subroutine dry2wet
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine GROWTH_FACTOR(S,Rd,R,LSUL)
    !// --------------------------------------------------------------------
    !// FITZGERALD (1975)
    !//
    !// From CTM2 to CTM3: Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: Rd   !// Dry radius
    logical, intent(in)   :: LSUL !// With or without sulphur
    !// Input/Output
    real(r8), intent(inout)  :: S    !// Supersaturation
    !// Output
    real(r8), intent(out)    :: R    !// Wet radius
    !// Locals
    real(r8) :: BETA,ALFA,TAU
    real(r8), parameter :: SMIN=0.3_r8, SMID=0.81_r8, SMAX = 0.995_r8
    !// --------------------------------------------------------------------

    IF (S .GT. SMAX)  S = SMAX

    IF(S.GT.SMIN .AND. S.LT.SMID) THEN
       !BETA = EXP(0.00077_r8*0.81_r8/(1.009_r8-0.81_r8))
       !ALFA = 1.2_r8*EXP(0.066_r8*0.81_r8/(1.058_r8-0.81_r8))
       BETA = EXP(0.00077_r8 * SMID / (1.009_r8 - SMID))
       ALFA = 1.2_r8*EXP(0.066_r8 * SMID / (1.058_r8 - SMID))
       IF (.NOT.LSUL) ALFA = ALFA * 1.35_r8
       !R = ((0.81_r8 - S)*Rd + (S - 0.3_r8)*ALFA*Rd**BETA)/0.51_r8
       R = ((SMID - S) * Rd + (S - SMIN) * ALFA * Rd**BETA) / (SMAX-SMIN)
    ELSE IF (S.GE.SMID .AND. S.LE.SMAX) THEN
       BETA = EXP(0.00077_r8 * S / (1.009_r8 - S))
       IF (S .LE. 0.97_r8) THEN
          TAU = 1.058_r8
       ELSE
          TAU = 1.058_r8 - (0.0155_r8*(S - 0.97_r8) / (1.02_r8 - S**1.4_r8))
       END IF
       ALFA = 1.2_r8 * EXP(0.066_r8 * S / (TAU - S))
       IF (.NOT.LSUL) ALFA = ALFA * 1.35_r8
       R = ALFA * Rd**BETA
    ELSE
       !// S .le. SMIN
       R = Rd
    END IF

    !// --------------------------------------------------------------------
  end subroutine growth_factor
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine FALLING( TEM, ZDV, AREAXY, Dpart, rhopart, sloss)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Calculate the falling velocity of an aerosol in all other layers
    !// than the first layer. Turbulent mix out is not taken into account
    !// here, only sedimentation. 
    !//
    !// Author: Alf Grini
    !//
    !// To CTM3: Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: TEM(LPAR)               !// Temp. of grid cell (K) 
    real(r8), intent(in) :: ZDV(LPAR)               !// Inv. vol. of box (m-3)
    real(r8), intent(in) :: AREAXY                  !// Area of grid cell (m2)
    real(r8), intent(in) :: Dpart(LPAR,Nsaltbins+1) !// Aerosol diam limits (m)
    real(r8), intent(in) :: rhopart(LPAR,Nsaltbins) !// Aerosol density (kg/m3)

    !// Output
    real(r8), intent(out) :: &
         sloss(LPAR,Nsaltbins)  !// Loss of aerosols due to falling (1/s)

    !// Local
    integer :: Ninner       !// Counter for sub bin calculation
    integer :: deltar       !// Size of bin (m) in sub bin calculation
    real(r8) :: height     !// Height of grid box (m)
    real(r8) :: viscair    !// Viscosity of air (kg m ???) 
    real(r8) :: rbin       !// Radius of bin (m)
    real(r8) :: vst        !// Falling velocity (m/s)
    real(r8) :: vstavg     !//    averaged over sub bins (if any) (m/s) 
    integer :: i,l,n                       ! Counting indexes for lon,lev,bin

    !// Local parameters
    integer, parameter :: Nsubsteps = 1        !// Number of sub bin loops
    real(r8), parameter  :: &
         visc0air = 1.75e-5_r8, & !// Viscosity of air at temperature T0visc
         T0visc = TK_0C, &        !// Temperature at which visocity=visc0 (K)  
         nviscair = 0.65_r8       !// Power law dependence viscosity on T
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='FALLING'
    !// --------------------------------------------------------------------

    !// Loop through bins (tracers)
    DO N = 1, Nsaltbins

       !// Loop through levels
       DO L = 1, LPAR

          !// Size of sub grid step in falling velocities (m)
          deltar = 0.5_r8 * (Dpart(L,N+1) - Dpart(L,N)) / real(Nsubsteps, r8)

          !// Viscocity of air at ambient temperature
          viscair = visc0air * (TEM(L) / T0visc)**nviscair

          !// Height of grid cell 
          HEIGHT = 1._r8 / (AREAXY * ZDV(L))

          !// Initializing falling velocity before averaging over sub grid 
          vstavg = 0._r8

          !// Check on if we will do sub grid falling calculations
          IF (Nsubsteps .eq. 1) THEN  !Will not do sub grid calculations

             rbin = 0.5_r8 * sqrt(Dpart(L,N) * Dpart(L,N+1))

             vst = 1._r8 / 18._r8 * (2._r8*rbin) * (2._r8*rbin) &
                   * rhopart(L,N) * G0 / viscair

             vstavg = vstavg + vst

          ELSE             !Will Average over sub bin falling velocities

             !// Lowest limit of bins
             rbin = 0.5_r8 * Dpart(L,N)

             DO Ninner = 1, Nsubsteps

                !// Radius of size class
                rbin = rbin + deltar
                !// Falling velocity ingoring the cunningham correction factor
                vst = 1._r8 / 18._r8 * (2._r8*rbin) * (2._r8*rbin) &
                     * rhopart(L,N) * G0 / viscair

                !// Adding the new falling velocity to the old one
                vstavg = vstavg + vst

             END DO !// Loop on sub grid falling velocities

          END IF !// DO sub grid calculation or not

          sloss(L,N) = (vstavg / real(Nsubsteps, r8)) / HEIGHT

          if (LDEBUG_SALT) then
             if (sloss(L,N) .ne. sloss(L,N)) then
                write(6,'(a,2i7)') f90file//':'//subr// &
                     ': SALT FALLING (sloss) is NAN at L,N: ',L,N
                write(6,'(a,i7)') 'Nsubsteps: ',Nsubsteps
                write(6,'(a,es20.12)') 'vstavg: ',vstavg
                write(6,'(a,es20.12)') 'vst:    ',vst
                write(6,'(a,es20.12)') 'HEIGHT: ',HEIGHT
                write(6,'(a,es20.12)') 'deltar: ',deltar
                write(6,'(a,es20.12)') 'rbin:   ',rbin
                write(6,'(a,es20.12)') 'DPART(L,N):   ',DPART(L,N)
                write(6,'(a,es20.12)') 'DPART(L,N+1): ',DPART(L,N+1)
                write(6,'(a,es20.12)') 'rhopart(L,N): ',rhopart(L,N)
                write(6,'(a,es20.12)') 'viscair: ',viscair
                stop
             end if
          end if

       END DO !// Loop on levels

    END DO !// Loop on bins

    !// --------------------------------------------------------------------
  end subroutine FALLING
  !// ----------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine DRYDEPPART(TEM, AIR, ZDV, AREAXY, Dpart, ROP, SSHF, Z0, USTR, &
       SFT, LOSS)
    !// ------------------------------------------------------------------
    !// PURPOSE: 
    !// Calculate drydeposition as done by Seinfeld and Pandis (1998).
    !// The variable SSHF(Surface sensible heat flux) is a meteorological
    !// field and is used as turbulent heat flux.
    !//
    !// ONLY to be used at model surface!
    !//
    !// Author: Alf Grini
    !//
    !// Theory:
    !//   J.H. Seinfeld and S.N. Pandis: 
    !//   Atmospheric chemistry and physics, John Wiley and Sons,
    !//   New York, 1998, Chapter 9.4, pp 962-972
    !//
    !// To CTM3: Ole Amund Sovde, September 2009
    !// ------------------------------------------------------------------
    use utilities, only: moninobukhov_length
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// INPUT
    real(r8), intent(in) :: &
         TEM, &       !Temperature (K)
         AIR, &       !Weight of air in cell (kg)
         ZDV, &       !Inverse of volume (m-3)
         AREAXY, &             !area of grid cell (m2)
         Dpart(Nsaltbins+1), & !Diameter of aerosols (m)
         rop(Nsaltbins), &     !Density of particles (kg/m3)
         SSHF, &               !Surface sensible heat flux (W/m2)
         z0, &                 !Roughness length (m)
         ustr, &               !Friction wind speed (m/s)
         SFT                   !Surface temperature (K)

    !// OUTPUT
    real(r8), intent(out) :: &
         loss(Nsaltbins)  ! Loss from lowest layer due to dry deposition (1/s)

    !// LOCAL
    integer :: I           !Counting variable for longitudes
    integer :: N           !Counting variable for bins
    integer :: Ninner      !Couting variable for loop inside a bin
    real(r8) :: vdep        !dry deposition velocity (m/s)
    real(r8) :: ra          !Aerodynamic resistance (s/m) 
    real(r8) :: rb          !laminar layer resistance (s/m)
    real(r8) :: Sc          !Schmidt number (UNITS ??)
    real(r8) :: St          !Stokes number (UNITS ??)
    real(r8) :: z           !Height of middle of layer (m)
    real(r8) :: vg          !Setteling velocity due to gravity only (m/s)
    real(r8) :: etacoeff    !Parameter describing stability (check Seinfeld)
    real(r8) :: tanetacoeff !Parameter describing stability (check Seinfeld)
    real(r8) :: Dbin        !Diameter of bin in question
    real(r8) :: DeltaD      !length of a size bin (m)
    real(r8) :: kji,kji0    !Parameters describing stability (check Seinfeld)
    real(r8) :: eta, eta0   !Parameters describing stability (check Seinfeld)
    real(r8) :: LMO         !Monin Obukhov length (m)
    real(r8) :: Trel        !Relaxation time (s)
    real(r8) :: Dpdiff      !Particle diffusivity in air (m2/s) (?)
    real(r8) :: roair       !Density of air (kg/m3)
    real(r8) :: viscair     !Viscocity of air (kg/(m*s))
    real(r8) :: ustr_loc    !local ustr
    !// LOCAL PARAMETERS
    integer, parameter :: Nsubsteps = 1   !Max inner loops inside a bin
    real(r8), parameter :: &
         vonkarman = 0.4_r8, & !Von karman constant
         T0visc = TK_0C, &     !Standard temperature for air viscosity
         kboltz = BOLTZMANN, & !Boltzman's constant (Pa*m3/(K*molecule))
         nviscair = 0.65_r8, & !Parameter to get air viscosity @ different temp
         visc0air = 1.75e-5_r8, & !Standard air viscosity (UNITS ???)
         Cp = cp_air           !Heat capacity of air (m2 s-2 K-1)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='DRYDEPPART'
    !// --------------------------------------------------------------------

    !// LOSS TERM IS ZERO INITIALLY
    LOSS(:) = 0._r8

    DO N = 1, Nsaltbins


       deltaD=(Dpart(N+1) - Dpart(N)) / real(Nsubsteps, r8)

       DO Ninner = 1, Nsubsteps

          !// initialize vdep
          vdep = 0._r8

          !// CALCULATING PARAMETERS NECESSARY FOR THE FURTHER CALCUATIONS
          !// --------------------------------------------------------------
          !// RADIUS OF "ACTIVE" BIN
          IF (Ninner .eq. 1) THEN
             Dbin = sqrt(Dpart(N) * Dpart(N+1))
          ELSE
             Dbin = Dbin + deltaD
          END IF
          !// VISCOCITY OF AIR
          viscair = visc0air * (TEM / T0visc)**nviscair

          !// HEIGHT IN MIDDLE OF BOX (ONLY FOR SURFACE LAYER!!!)
          z = 0.5_r8 / (ZDV * AREAXY)

          !// AIR DENISITY
          roair = AIR * ZDV

          !// PARTICLE DIFFUSIVITY
          Dpdiff = kboltz * TEM / (3._r8 * CPI * viscair * Dbin)

          !// STOKES FALLING VELOCITY
          vg = 1._r8 / 18._r8 * Dbin*Dbin * rop(N) * G0 / viscair

          !// PARTICLE RELAXATION TIME
          Trel = vg / G0

          !// STOKES NUMBER
          St = vg * ustr * ustr / (G0 * viscair / roair)

          !// SCHMIDT NUMBER
          Sc = (viscair / roair) / Dpdiff


          !// CALCULATING MICROMETEOROLOGICAL PARAMETERS
          !// ------------------------------------------------------------------
          !// Ensure non-zero ustr
          if (ustr .gt. 0._r8) then
             ustr_loc = ustr
          else
             ustr_loc = 5.e-3_r8
          end if
          !// MONIN OBUKHOV LENGTH
          LMO = moninobukhov_length(roair,SFT,ustr_loc,SSHF)

          !// Limiting abs(z/L) to be between approx (-1 and 1) 
          IF (LMO.gt.0 .and. LMO.lt.z) LMO = z
          IF (LMO.lt.0 .and. LMO.gt.(-1._r8*z)) LMO = -1._r8 * z
          !// Cannot allow zero LMO (z is always positive)
          if (LMO.eq.0._r8) LMO = z

          !// KJI (repeat limit)
          KJI = max(min(z/LMO, 1._r8), -1._r8)

          !// KJI0
          KJI0 = z0 / LMO



          !// CHECKING FOR STABLE/UNSTABLE/NEUTRAL REGIME
          !// ------------------------------------------------------------------
          IF (KJI.gt.0._r8 .and. KJI.le.1._r8) THEN

             !// STABLE REGIME; ra for stable
             ra = 1._r8 / (vonkarman * ustr_loc) &
                  * (log(z/z0) + 4.7_r8 * (kji - kji0))
                  
          ELSE IF (KJI.eq.0._r8) THEN 

             !// NEUTRAL REGIME; ra for neutral
             ra = 1._r8 / (vonkarman * ustr_loc) * log(z / z0)

          ELSE !IF (KJI.lt.0._r8 .AND. KJI.ge.-1._r8) THEN

             !// UNSTABLE REGIME; ra for unstable
             !// ETA
             ETA = (1._r8 - 15._r8 * Kji)**0.25_r8
             !// ETA0
             ETA0 = (1._r8 - 15._r8 * Kji0)**0.25_r8
             etacoeff = ((eta0**2 + 1._r8) * (eta0 + 1._r8)**2) &
                                  / ((eta**2 + 1._r8) * (eta + 1._r8)**2)

             tanetacoeff = atan(eta) - atan(eta0)

             ra = 1._r8 / (vonkarman * ustr_loc) * &
                  (log(z / z0) + log(etacoeff) + 2._r8 * tanetacoeff)

             IF (.not.(KJI.lt.0._r8 .AND. KJI.ge.-1._r8)) THEN
                !// This should not be possible
                write(6,'(a,i7)') f90file//':'//subr// &
                     ': WARNING unstable regime, check KJI: ', KJI
                write(6,'(a,4es20.12)') '  ra, ustr_loc, z, z0: ', &
                     ra, ustr_loc, z, z0
                write(6,'(a,3es20.12)') '  eta,eta0,tanetacoeff: ', &
                     eta,eta0,tanetacoeff
             end if

          END IF !// Check on turbulence (stable/unstable)


          !// CALCULATING LAMINAR LAYER RESISTANCE
          !// ------------------------------------------------------------------
          rb = 1._r8 / (ustr_loc * (Sc**(-0.667_r8) + 10._r8**(-3._r8/St)))


          !// CALCULATING TOTAL DRYDEP-VELOCITY FOR A BIN
          !// ------------------------------------------------------------------
          vdep = 1._r8 / (ra + rb + ra * rb * vg) + vg + vdep

          IF(vdep.lt.0._r8)THEN
             write(6,'(a,es20.12)') f90file//':'//subr// &
                  ': NEGATIVE vdep: ',vdep
             write(6,*)'N,Ninner',N,Ninner
             write(6,*)'ra = ',ra, log(z/z0), 4.7_r8*(kji-kji0)
             write(6,*)'rb=  ',rb
             write(6,*)'sshf ',sshf
             write(6,*)'ustr',ustr_loc
             write(6,*)'z0',z0
             write(6,*)'z',z
             write(6,*)'KJI',KJI
             write(6,*)'1/ZDV',1._r8/ZDV
             stop
          END IF

          !// DEBUG
          if (LDEBUG_SALT) then
             if (vdep .ne. vdep) then
                write(6,'(a,es20.12)') f90file//':'//subr// &
                     ': vdep NAN: ',vdep
                write(6,*)'N,Ninner',N,Ninner
                write(6,*)'lmo= ',LMO
                write(6,*)'rb=  ',rb
                write(6,*)'sshf ',sshf
                write(6,*)'ustr',ustr_loc
                write(6,*)'z0',z0
                write(6,*)'z',z
                write(6,*)'KJI',KJI
                write(6,*)'DV/st/sc',1._r8/ZDV,St,Sc
                IF (KJI.gt.0._r8 .and. KJI.le.1._r8) THEN
                   !// STABLE REGIME; ra for stable
                   write(6,*)'ra stable = ',ra, log(z/z0), 4.7_r8*(kji-kji0)
                ELSE IF (KJI.eq.0._r8) THEN 
                   !// NEUTRAL REGIME; ra for neutral
                   write(6,*)'ra neutral = ',ra, log(z/z0)
                ELSE IF (KJI.lt.0._r8 .AND. KJI.ge.-1._r8) THEN
                   !// UNSTABLE REGIME; ra for unstable
                   write(6,*)'ra unstable = ',ra, log(z/z0), log(etacoeff), tanetacoeff
                ELSE
                   write(6,*)'ra unknown = ',ra, z, z0,eta,eta0
                END IF !// Check on turbulence (stable/unstable)
                stop
             end if
          end if


          !// FINISHED WITH INNER LOOP
       END DO               !Ninner

       !// AVERAGING THE VELOCITIES FOUND
       vdep = vdep / real(Nsubsteps, r8)

       !write(6,*)'D,vdep',sqrt(Dpart(N)*Dpart(N+1))*1.e6,'um',vdep*100.,'cm/s'


       !// CALCULATING LOSS TERM
       LOSS(N) = vdep / (2._r8 * Z) !(Loss speed devided by height of layer)

    END DO !// DO N = 1, Nsaltbins

    !// --------------------------------------------------------------------
  end subroutine DRYDEPPART
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine z0_get(USTR, PLAND, T, Z0)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Modify z0 for momentum (and thus here, drydep). If we have Ice or
    !// sea-ice
    !//
    !// Theory:
    !//   Stull, R.B. : An introduction to boundary layer meteorology
    !//   Kluwer academic publishers, Dortrecht, 1994
    !//   pp: 380-381, eq: 9.7.2.c) and table 9.6
    !//
    !//  Author: Alf Grini
    !//
    !// To CTM3: Ole Amund Sovde, September 2009
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    real(r8), intent(in)      :: ustr  !Wind friction speed (m/s)
    real(r8), intent(in)      :: T     !Temperature(K)
    real(r8), intent(in)      :: PLAND !Percentage of land (frc)

    !// In/out
    real(r8), intent(inout) :: z0  !Roughness length (m) (modified at sea)

    !// Local
    real(r8), parameter :: &
         Tfreeze = TK_0C, &  !Freezing point (K)
         z0seaice = 1.e-4_r8 !m (approximately, Stull pp 380) 
    integer :: I             !Counting variable for longitude
    !// ------------------------------------------------------------------

    !// Modify z0 read in initsalt only if we have sea or sea_ice

    IF (PLAND.lt.0.5_r8 .AND. T.gt.Tfreeze) THEN
       !// Unfrozen ocean, use Charnock relationship
       !// Stull, pp381, eq 9.7.2.c

       z0 = 0.016_r8 * ustr * ustr / G0
            
    ELSE IF (PLAND.lt.0.5_r8 .AND. T.le.Tfreeze) THEN
       !// Frozen ocean, set fixed z0
       z0 = z0seaice
    END IF

    !// Don't let z0 be lower than that of sea/ice
    z0 = max(z0seaice, z0)

    !// --------------------------------------------------------------------
  end subroutine z0_get
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine addtogether(SALTT,prodseasalt, sloss,areaxy,dtime,drydep,I,J)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Add together sea salt to get sea salt concentrations knowing the
    !// respective production and loss terms.
    !//
    !// Author: Alf Grini
    !//
    !// Changes from Oslo CTM2 to Oslo CTM3
    !//   - Only loops over NSALTBINS and vertical!
    !//   - Adds to budget diagnostics: production and dry deposition
    !//
    !// Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: prodseasalt(Nsaltbins) !// Prod. of seasalt (kg/s)
    real(r8), intent(in) :: sloss(LPAR,Nsaltbins)  !// Loss due to falling/mix out (1/s)
    real(r8), intent(in) :: areaxy            !// Area of grid square (m^2)
    real(r8), intent(in) :: dtime             !// Timestep (s)
    integer, intent(in) :: I,J                !// Global indices

    !// Input/Output
    !// Salt species (kg)
    real(r8), intent(inout) :: SALTT(LPAR,Nsaltbins)
    real(r8), intent(out)   :: DRYDEP(Nsaltbins) !// Dry deposited (kg)

    !// Local
    real(r8) :: &
         Cold, &                 !// Old concentration before timestep (kg)
         Cnew, &                 !// New concentration after timestep (kg)
         lost(LPAR,Nsaltbins), & !// Lost from layer (kg)
         stau, &               !// Lifetime with respect to loss term (s)
         sprod, &              !// Production in kg/s
         prod_dgn, &           !// Production to output(kg/m2/s)
         ddep_dgn              !// Drydep to output (kg/m2/s)
    integer :: L, &        !// Couting variable for levels
               N           !// Counting variable for saltbins
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='addtogether'
    !// --------------------------------------------------------------------

    drydep(:) = 0._r8

    !// Start loop on saltbins
    do N = 1, Nsaltbins

       !// Loop on Layers: Going from the top and down !!
       do L = LPAR, 1, -1

          !// Nothing is lost at this level yet
          lost(L,N) = 0._r8

          IF (L.eq.LPAR) THEN
             !// No production in the uppermost layer
             sprod = 0._r8
          ELSE IF (L.eq.1) THEN
             !// Production is what is lost from layer above 
             !// And what is produced from sea surface
             sprod = lost(L+1,N)/dtime + prodseasalt(N)
          ELSE !//IF (L.lt.LPAR .AND. L.gt.1) THEN
             !// Production is what is lost from layer above
             sprod = lost(L+1,N)/dtime
          END IF !// Check on what layer we are in

          !// Using the QSSA solver to calculate the new concentration
          !// and how much is lost. Save old value first:
          Cold = SALTT(L,N)


          if (sloss(L,N) .eq. 0._r8) then
             !// QSSA without loss
             Cnew = Cold + sprod*dtime
          else
             !// Have loss term, calculate lifetime
             stau = 1._r8 / sloss(L,N)

             !// DEBUG
             if (sloss(L,N) .lt. 0._r8) then
                write(6,'(a,2i5,es20.12)') f90file//':'//subr// &
                     ': sloss NEGATIVE: ',L,N,sloss(L,N)
                stop
             end if

             !// Solve as Euler forward or as QSSA
             IF (stau .gt. 100._r8*dtime) THEN
                !// Solve Euler forward equation
                Cnew = Cold + (sprod - sloss(L,N)*Cold)*dtime
             ELSE  
                !// Solve QSSA equation
                Cnew = sprod / sloss(L,N) + &
                     (Cold - sprod/sloss(L,N)) * exp(-sloss(L,N)*dtime)
             END IF !// Choice on QSSA

          end if !// if (sloss(L,N) .eq. 0._r8) then

          !// DEBUG
          if (LDEBUG_SALT) then
             if (Cnew .ne. Cnew) then
                write(6,'(a,2i5,4es20.12)') f90file//':'//subr// &
                     ': Cnew NAN',N,L,Cold, sprod, sloss(L,N),dtime
                stop
             end if
          end if

          !// Update STT, new production and loss terms
          SALTT(L,N) = Cnew

          !// Some of the produced salt may have been lost (i.e. sedimented
          !// downwards), so the total lost salt is approximated assuming
          !// linear production:
          !//    (Cnew - Cold) = PROD - LOSS
          !//    LOSS = PROD + Cold - Cnew
          !// where LOSS is positive; LOSS can NOT be allowed to become negative
          lost(L,N) = max(0._r8, sprod*dtime + Cold - Cnew)

       end do !// DO L = LPAR, 1, -1
  
    end do !// DO N = 1, Nsaltbins


    !// Initialize production (sum of all bins)
    prod_dgn = 0._r8

    !// Initialize drydep (sum of all bins)
    ddep_dgn = 0._r8

    !// Add up total production in kg/m2
    do N = 1, Nsaltbins
       prod_dgn = prod_dgn + prodseasalt(N) / AREAXY * dtime !kg/m2/s * s
    end do

    do N = 1, Nsaltbins
       ddep_dgn = ddep_dgn + lost(1,N) / AREAXY !kg/m2
       !// Save drydep [kg] for global diagnose
       drydep(N) = lost(1,N)
    end do

    !// --------------------------------------------------------------------
  end subroutine addtogether
  !// ----------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine add2drydepdiag(ii,jj,mp,n,ddmass)
    !// ------------------------------------------------------------------
    !// Saves SALT tracer mass lost to drydep.
    !//
    !// Ole Amund Sovde, September 2015
    !// ------------------------------------------------------------------
    use cmn_ctm, only: all_mp_indices
    use cmn_oslo, only: SCAV_MAP_DRY, SCAV_DD
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    integer, intent(in) :: ii,jj,mp,n
    real(r8), intent(in) :: ddmass
    !// ------------------------------------------------------------------
    SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) + ddmass
    !// Save totals
    SCAV_DD(N,MP) = SCAV_DD(N,MP) + ddmass
    !// ------------------------------------------------------------------
  end subroutine add2drydepdiag
  !// ------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine saltbdg2file(NDAY, NDAYI, NDAY0)
    !// --------------------------------------------------------------------
    !// Writes out SALT budget averages, more comprehensive than
    !// the CTM2 saltbdg2d: Can also write budgets for all dust tracers.
    !// Depends on core diagnostics, as it fetches the diagnostics of wet
    !// scavenging from the STTTND array.
    !// Also writes totals and lifetimes to screen.
    !//
    !// Ole Amund Sovde, March 2016
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, NTDPAR, MODEL
    use cmn_ctm, only: XDGRD, YDGRD, IYEAR, AREAXY, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_diag, only: STTTND, NTND_LSSCAV
    use cmn_oslo, only: CONVWASHOUT
    use netcdf
    use ncutils, only:  handle_error, handle_err
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI, NDAY0
    
    !// Locals
    integer :: I, J, II, JJ, MP, L, N, K     !// Indices
    real(r8) :: davgdays, davgtnd  !// Time durations
    real(r8) :: burden,lscv,cscv,prod,ddep
    real(r8) :: ZAREA(IPAR,JPAR)
    integer :: strLen
    !// Total budgets (summed up over SALT tracers)
    real(r8), dimension(IPAR,JPAR) :: &
         sltTotBdgProd, &
         sltTotBdgDdep, &
         sltTotBdgLscv, &
         sltTotBdgCscv, & 
         sltTotBdgBrdn

    !// --------------------------------------------------------------------
    !// netCDF: Dimensions ID and counters
    integer :: &
         lat_dim_id, &             !// Dimension ID for latitude
         lon_dim_id, &             !// Dimension ID for longitude
         nbr_dim_id, &             !// Dimension ID for number of tracers
         time_dim_id, &            !// Dimension ID for time
         timespan_dim_id, &        !// Dimension ID for timespan
         str_len_id, &             !// Dimension ID for string length
         time_id, &                !// Variable ID for time
         timespan_id, &            !// Variable ID for average time span (day)
         lon_id, &                 !// Variable ID for longitude
         lat_id, &                 !// Variable ID for latitude
         prod_id, ddep_id, lscv_id, cscv_id, brdn_id, &
         tot_prod_id, tot_ddep_id, tot_lscv_id, tot_cscv_id, tot_brdn_id, &
         gridarea_id, &            !// Variable ID for gridarea
         saltnames_id              !// Variable ID for saltnames
    integer :: &
         srt_time(1), &            !// starting point for time array
         dim_lon_lat_time(3), &    !// Dimension ID for processes
         srt_lon_lat_time(3), &    !// Start array for lon/lat/time
         cnt_lon_lat_time(3), &       !// Counting array for lon/lat/time
         dim_nbr_lon_lat_time(4), &    !// Dimension ID for processes
         srt_nbr_lon_lat_time(4), &    !// Start array for nbr/lon/lat/time
         cnt_nbr_lon_lat_time(4)       !// Counting array for nbr/lon/lat/time

    !// Other locals
    integer :: &
         ncid, &                   !// File ID for nc file
         status, &                 !// Error status for nc file
         nlons, &                  !// Number of longitudes found in file
         nlats, &                  !// Number of latitudes found in file
         nbrs, &                   !// Number of tracers found in file
         nsteps                    !// Number of steps found in file 
    character(len=4)      :: cyear                  !// Year in character*4
    character(len=3)      :: cday                   !// Day in character*3
    character(len=80)     :: time_label             !// Label for variable "time"
    real(r8)              :: time                   !// Time in this timestep
    character(len=:), allocatable :: saltnames(:)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'saltbdg2file'
    !// --------------------------------------------------------------------

    !// Get budgets from STTTND (6=LSSCAV, 8=C_SCAV)
    !// These tendencies are summed up (kg/gridbox) from NDAY0 to NDAY
    davgtnd = 86400._r8 * real(NDAY+1 - NDAY0, r8)

    !// Check the time vs time for other diagnoses
    davgdays = saltbudgetacctime
    if (davgdays .ne. davgtnd) then
       write(6,'(a,2f12.3,2i8)') f90file//':'//subr//': '// &
            'wrong dt',davgdays,davgtnd, NDAY0, NDAY
       stop
    end if

    !// Large scale wash out [kg] (accumulated from NDAY0 to NDAY)
    salt_lscv(:,:,:) = 0._r8
    do K = 1, Nsaltbins
       N = salt_trsp_idx(K)
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                !// Save as positive values: STTTND is BTT-BTTBCK.
                salt_lscv(K,I,J) = salt_lscv(K,I,J) &
                     - STTTND(I,J,L,N,NTND_LSSCAV)
             end do
          end do
       end do
    end do

    !// Convective wash out [kg] (accumulated from NDAY0 to NDAY)
    salt_cscv(:,:,:) = 0._r8
    do K = 1, Nsaltbins
      N = salt_trsp_idx(K)
      do MP = 1, MPBLK
        do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = 1, LPAR
              !// Save as positive values: CONVWASHOUT is neg for removal.
               salt_cscv(K,I,J) = salt_cscv(K,I,J) &
                    - CONVWASHOUT(L,N,II,JJ,MP)
            end do
          end do
        end do
      end do
    end do

    !// Generate average; i.e. divide by time duration of accumulation,
    !// which is davgtnd: [kg] -> [kg/s]
    salt_production(:,:,:) = salt_production(:,:,:) / davgtnd
    salt_drydep(:,:,:) = salt_drydep(:,:,:) / davgtnd
    salt_lscv(:,:,:) = salt_lscv(:,:,:) / davgtnd
    salt_cscv(:,:,:) = salt_cscv(:,:,:) / davgtnd
    !// For burden, this means [kg*s] -> [kg]
    salt_burden(:,:,:) = salt_burden(:,:,:) / davgtnd

    !// Average prod/ddep/lscv/cscv/burden during period (already divided by dt)
    prod   = sum(salt_production)
    ddep   = -sum(salt_drydep)
    lscv   = -sum(salt_lscv)
    cscv   = -sum(salt_cscv)
    burden = sum(salt_burden)

    write(6,'(a,6es10.2)') 'SALT AVG P/D/LS/CN[kg/s]/B[kg]', &
         prod, ddep, lscv, cscv, burden
    write(6,'(a,2f11.3)')  'SALT LIFE B/P, B/L [days]', &
         burden/prod / 86400._r8,&
         -burden/(ddep+lscv+cscv) / 86400._r8

    !// Prepare for output
    !// --------------------------------------------------------------------

    !// For totals tendencies [kg/s] and brdn [kg]
    do J = 1, JPAR
       do I = 1, IPAR
          sltTotBdgProd(I,J) = sum(salt_production(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          sltTotBdgDdep(I,J) = sum(salt_drydep(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          sltTotBdgLscv(I,J) = sum(salt_lscv(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          sltTotBdgCscv(I,J) = sum(salt_cscv(:,I,J))
       end do
    end do
    do J = 1, JPAR
       do I = 1, IPAR
          sltTotBdgBrdn(I,J) = sum(salt_burden(:,I,J))
       end do
    end do

    !// Divide by area [kg] -> [kg/m2]
    ZAREA(:,:) = 1._r8 / AREAXY(:,:)
    do N = 1, nsaltbins
       salt_production(N,:,:) = salt_production(N,:,:) * ZAREA(:,:)
       salt_drydep(N,:,:) = salt_drydep(N,:,:) * ZAREA(:,:)
       salt_lscv(N,:,:) = salt_lscv(N,:,:) * ZAREA(:,:)
       salt_cscv(N,:,:) = salt_cscv(N,:,:) * ZAREA(:,:)
       salt_burden(N,:,:) = salt_burden(N,:,:) * ZAREA(:,:)
    end do
    !// For totals and brdn
    sltTotBdgProd(:,:) = sltTotBdgProd(:,:) * ZAREA(:,:)
    sltTotBdgDdep(:,:) = sltTotBdgDdep(:,:) * ZAREA(:,:)
    sltTotBdgLscv(:,:) = sltTotBdgLscv(:,:) * ZAREA(:,:)
    sltTotBdgCscv(:,:) = sltTotBdgCscv(:,:) * ZAREA(:,:)
    sltTotBdgBrdn(:,:) = sltTotBdgBrdn(:,:) * ZAREA(:,:)

    !// The total time span for the average [days]
    davgdays = davgdays / 86400._r8

    !// Write next item to file (it is initialized as 0)
    saltaverages_written = saltaverages_written + 1


    !// --------------------------------------------------------------------
    !// Generate netCDF file
    status = nf90_noerr

    !// First time routine is called, we need to define variables
    if (saltaverages_written .eq. 1) then
       !// Create file
       write(6,'(a)') f90file//':'//subr//': creating '//trim(saltbdgfile)
       !// Clobber means you can overwrite existing data
       status=nf90_create(path=saltbdgfile,cmode=nf90_clobber,ncid=ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error opening '//trim(saltbdgfile))

       !// File headers. The group should agree on these so that output from
       !// our model always has the same headers.
       status = nf90_put_att(ncid,nf90_global,'title', &
            trim(MODEL)//' salt averaged fields')
       if (status .ne. nf90_noerr) call handle_error(status,'attributing title')

       !// Define dimensions
       status=nf90_def_dim(ncid,'lon',IPAR,lon_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim lon')
       status=nf90_def_dim(ncid,'lat',JPAR,lat_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim lat')
       status=nf90_def_dim(ncid,'Nsaltbins',nsaltbins,nbr_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim Nsaltbins')
       strLen = len(TNAME(salt_trsp_idx(1))) !// Length of TNAMEs
       status=nf90_def_dim(ncid,'str_len',strLen,str_len_id)
       if (status .ne. nf90_noerr) call &
            handle_error(status,'defining dim strLen')
       status=nf90_def_dim(ncid,'time',nf90_unlimited,time_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'defining dim time')

       !// Defining the combined id for a field (lon /lat /time)
       dim_lon_lat_time(1)=lon_dim_id
       dim_lon_lat_time(2)=lat_dim_id
       dim_lon_lat_time(3)=time_dim_id

       dim_nbr_lon_lat_time(1)=nbr_dim_id
       dim_nbr_lon_lat_time(2)=lon_dim_id
       dim_nbr_lon_lat_time(3)=lat_dim_id
       dim_nbr_lon_lat_time(4)=time_dim_id

       !// Defining the lon/lat/time-variable
       status = nf90_def_var(ncid,'lon',nf90_double,lon_dim_id,lon_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining lon')
       status = nf90_def_var(ncid,'lat',nf90_double,lat_dim_id,lat_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining lat')
       status = nf90_def_var(ncid,'time',nf90_double,time_dim_id,time_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining time')
       !// Defining average time span
       status=nf90_def_var(ncid,'timespan',nf90_double,time_dim_id,timespan_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining timespan')
       !// Defining salt names
       status = nf90_def_var(ncid,'saltname',nf90_char, &
                                  (/str_len_id,nbr_dim_id/),saltnames_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltname')
       !// Defining grid box area
       status=nf90_def_var(ncid,'gridarea',nf90_double, &
                                (/lon_dim_id,lat_dim_id/),gridarea_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining gridarea')

       !// Putting attributes to /lon/lat/ variables
       !// Don't change the units. It is important to keep these for GRADS
       status=nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit lon')
       status=nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit lat')

       write(cyear,'(I4)')IYEAR
       write(cday,'(I3)')NDAYI
       time_label='days since '//cyear//': day number '//trim(cday)

       !// Putting units attribute to time variable
       status=nf90_put_att(ncid,time_id,'units',time_label)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit time')
       !// Average time span
       status=nf90_put_att(ncid,timespan_id,'units','Average time span (days)')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit timespan')
       !// Gridbox area
       status=nf90_put_att(ncid,gridarea_id,'units','Grid box area [m2]')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit gridarea')
       !// Salt names
       status=nf90_put_att(ncid,saltnames_id,'units','Salt tracer names')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltnames')

       !// Define the tracer field variables for totals and their units
       status = nf90_def_var(ncid,'saltprod',nf90_float,dim_lon_lat_time,tot_prod_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltprod')
       status = nf90_put_att(ncid,tot_prod_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltprod')
       status = nf90_put_att(ncid,tot_prod_id,'longname','Total salt production')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname saltprod')

       status = nf90_def_var(ncid,'saltddep',nf90_float,dim_lon_lat_time,tot_ddep_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltddep')
       status = nf90_put_att(ncid,tot_ddep_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltddep')
       status = nf90_put_att(ncid,tot_ddep_id,'longname','Total salt drydep')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname saltddep')

       status = nf90_def_var(ncid,'saltwdep',nf90_float,dim_lon_lat_time,tot_lscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltwdep')
       status = nf90_put_att(ncid,tot_lscv_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltwdep')
       status = nf90_put_att(ncid,tot_lscv_id,'longname', &
            'Total salt largescale scavenged')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname saltwdep')

       status = nf90_def_var(ncid,'saltcnvw',nf90_float,dim_lon_lat_time,tot_cscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltcnvw')
       status = nf90_put_att(ncid,tot_cscv_id,'units','kg/m2/s')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltcnvw')
       status = nf90_put_att(ncid,tot_cscv_id,'longname', &
            'Total salt convective scavenged')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname saltcnvw')

       status = nf90_def_var(ncid,'saltbrdn',nf90_float,dim_lon_lat_time,tot_brdn_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error defining saltbrdn')
       status = nf90_put_att(ncid,tot_brdn_id,'units','kg/m2')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error unit saltbrdn')
       status = nf90_put_att(ncid,tot_brdn_id,'longname','Total salt burden')
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error longname saltbrdn')


       !// Define the tracer field variables for all SALT tracers
       if (LSALTDIAG3D) then
          status = nf90_def_var(ncid,'salt_prod',nf90_float, &
               dim_nbr_lon_lat_time, prod_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining salt_prod')
          status = nf90_put_att(ncid,prod_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit salt_prod')
          status = nf90_put_att(ncid,prod_id,'longname', &
               'Salt production, all salt tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname salt_prod')

          status = nf90_def_var(ncid,'salt_ddep',nf90_float, &
               dim_nbr_lon_lat_time, ddep_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining salt_ddep')
          status = nf90_put_att(ncid,ddep_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit salt_ddep')
          status = nf90_put_att(ncid,ddep_id,'longname', &
               'Salt drydep, all salt tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname salt_ddep')

          status = nf90_def_var(ncid,'salt_lscv',nf90_float, &
               dim_nbr_lon_lat_time, lscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining salt_lscv')
          status = nf90_put_att(ncid,lscv_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit salt_lscv')
          status = nf90_put_att(ncid,lscv_id,'longname', &
               'Salt largescale scavenged, all salt tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname salt_lscv')

          status = nf90_def_var(ncid,'salt_cscv',nf90_float, &
               dim_nbr_lon_lat_time, cscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining salt_cscv')
          status = nf90_put_att(ncid,cscv_id,'units','kg/m2/s')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error unit salt_cscv')
          status = nf90_put_att(ncid,cscv_id,'longname', &
               'Salt convective scavenged, all salt tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname salt_cscv')

          status = nf90_def_var(ncid,'salt_brdn',nf90_float, &
               dim_nbr_lon_lat_time, brdn_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error defining salt_brdn')
          status = nf90_put_att(ncid,brdn_id,'units','kg/m2')
          if (status .ne. nf90_noerr)&
               call handle_error(status,'error unit salt_brdn')
          status = nf90_put_att(ncid,brdn_id,'longname', &
               'Salt burden, all salt tracers')
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error longname salt_brdn')
       end if !// if (LSALTDIAG3D) then

       !// End definition mode
       status = nf90_enddef(ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error nf90_enddef')
 
       !// Putting the lon/lat/areaxy variables
       status = nf90_put_var(ncid,lon_id,XDGRD)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting lon')
       status = nf90_put_var(ncid,lat_id,YDGRD)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting lat')
       status = nf90_put_var(ncid,gridarea_id,AREAXY)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting gridarea')
       !// salt names
       allocate( character(strLen) :: saltnames(nsaltbins) )
       do N = 1, nsaltbins
          saltnames(N) = TNAME(salt_trsp_idx(N))
       end do
       status = nf90_put_var(ncid,saltnames_id,saltnames)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltnames')
       deallocate( saltnames )
       write(6,'(a)') f90file//':'//subr//': created '//trim(saltbdgfile)

    else !// THE FILE HAS BEEN USED BEFORE
       !// Open the existing file
       write(6,'(a)') f90file//':'//subr//': opening '//trim(saltbdgfile)
       status = nf90_open(saltbdgfile, nf90_write, ncid)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error in open')
       !// Inquire dimension ids
       status = nf90_inq_dimid(ncid,'lat',lat_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting lat_id')
       status = nf90_inq_dimid(ncid,'lon',lon_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting lon_id')
       status = nf90_inq_dimid(ncid,'Nsaltbins',nbr_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting nbr_id')
       status = nf90_inq_dimid(ncid,'time',time_dim_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting time_id')
     
       !// Inquire dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim lat')
       if (nlats .ne. JPAR)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': JPAR mismatch on existing file: ',nlats,JPAR
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim lon')
       if (nlons .ne. IPAR)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': IPAR mismatch on existing file: ',nlons,IPAR
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,nbr_dim_id,len=nbrs)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim nbr')
       if (nbrs .ne. nsaltbins)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': nbrs mismatch on existing file: ',nbrs,nsaltbins
          stop
       endif
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting dim time')
       if (nsteps+1 .ne. saltaverages_written)then
          write(6,'(a,2i7)')f90file//':'//subr// &
               ': nsteps mismatch on existing file: ', &
               nsteps, saltaverages_written
          stop
       endif

       !// Get ids for totals
       status = nf90_inq_varid(ncid,'saltprod',tot_prod_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_prod_id')
       status = nf90_inq_varid(ncid,'saltddep',tot_ddep_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_ddep_id')
       status = nf90_inq_varid(ncid,'saltwdep',tot_lscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_lscv_id')
       status = nf90_inq_varid(ncid,'saltcnvw',tot_cscv_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_cscv_id')
       status = nf90_inq_varid(ncid,'saltbrdn',tot_brdn_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting tot_brdn_id')

       !// Get ids for all salt bins
       if (LSALTDIAG3D) then
          status = nf90_inq_varid(ncid,'salt_prod',prod_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting prod_id')
          status = nf90_inq_varid(ncid,'salt_ddep',ddep_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting ddep_id')
          status = nf90_inq_varid(ncid,'salt_lscv',lscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting lscv_id')
          status = nf90_inq_varid(ncid,'salt_cscv',cscv_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting cscv_id')
          status = nf90_inq_varid(ncid,'salt_brdn',brdn_id)
          if (status .ne. nf90_noerr) &
               call handle_error(status,'error getting brdn_id')
       end if !// if (LSALTDIAG3D) then

       !// Get variable id for time
       status=nf90_inq_varid(ncid,'time',time_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting time_id')
       status = nf90_inq_varid(ncid,'timespan',timespan_id)
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error getting timespan_id')

    end if !// Has the file already been in use or not ?

    !// --------------------------------------------------------------------
    !// Defining how far to count for each time a data set is added
    !// Adding 2D fields
    cnt_lon_lat_time = (/IPAR , JPAR,  1/)
    !// Defining where to start adding the new time step
    srt_lon_lat_time = (/1, 1, saltaverages_written/)

    !// Adding 3D fields
    cnt_nbr_lon_lat_time = (/nsaltbins, IPAR , JPAR,  1/)
    !// Defining where to start adding the new time step
    srt_nbr_lon_lat_time = (/1, 1, 1, saltaverages_written/)


    srt_time(1) = saltaverages_written    !// Start value for new time step
    time = real(NDAY - NDAYI + 1, r8)     !// Time in r8 format (days)
  
    status = nf90_put_var(ncid,time_id,time,start=srt_time)
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting time')

    !// Also put out average time span
    status = nf90_put_var(ncid,timespan_id,davgdays,start=srt_time)
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting avg time span')


    !// Put tracer fields for totals
    status = nf90_put_var(ncid, tot_prod_id, sltTotBdgProd(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltprod')

    status = nf90_put_var(ncid, tot_ddep_id, sltTotBdgDdep(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltddep')

    status = nf90_put_var(ncid, tot_lscv_id, sltTotBdgLscv(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltwdep')

    status = nf90_put_var(ncid, tot_cscv_id, sltTotBdgCscv(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltcnvw')

    status = nf90_put_var(ncid, tot_brdn_id, sltTotBdgBrdn(:,:), &
         start=srt_lon_lat_time, &   !// starting point for writing
         count=cnt_lon_lat_time )    !// Counts how many bytes written
    if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting saltbrdn')


    !// Put fields for all tracers
    if (LSALTDIAG3D) then
       !// Put out salt_prod
       status = nf90_put_var(ncid, prod_id, &  !// File id, variable id
            salt_production(:,:,:), &       !// Salt produced
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting salt_prod')

       !// Put out salt_ddep
       status = nf90_put_var(ncid, ddep_id, &  !// File id, variable id
            salt_drydep(:,:,:), &           !// Salt drydepped
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting salt_ddep')

       !// Put out salt_lscv
       status = nf90_put_var(ncid, lscv_id, &  !// File id, variable id
            salt_lscv(:,:,:), &             !// Salt ls scav
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting salt_lscv')

       !// Put out salt_cscv
       status = nf90_put_var(ncid, cscv_id, &  !// File id, variable id
            salt_cscv(:,:,:), &             !// Salt cnv scav
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting salt_cscv')

       !// Put out salt_brdn
       status = nf90_put_var(ncid, brdn_id, &  !// File id, variable id
            salt_burden(:,:,:), &           !// Salt burden
            start=srt_nbr_lon_lat_time, &   !// starting point for writing
            count=cnt_nbr_lon_lat_time )    !// Counts how many bytes written
       if (status .ne. nf90_noerr) &
            call handle_error(status,'error putting salt_brdn')
    end if !// if (LSALTDIAG3D) then
  
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         subr//':close:'//trim(saltbdgfile))

    write(6,'(a)') f90file//':'//subr//': updated '//trim(saltbdgfile)


    !// Re-initialize
    saltbudgetacctime = 0._r8

    salt_production(:,:,:) = 0._r8
    salt_drydep(:,:,:)     = 0._r8
    salt_lscv(:,:,:)       = 0._r8
    salt_cscv(:,:,:)       = 0._r8
    salt_burden(:,:,:)     = 0._r8

    !// --------------------------------------------------------------------
  end subroutine saltbdg2file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasaltbudget_lscav(BTT,BTTBCK,DT,MP)
    !// --------------------------------------------------------------------
    !// Possible routine to diagnose large scale scavenging of SALT
    !// from BTT and BTTBCK.
    !// No need for this, the numbers are diagnosed by STTTND.
    !//
    !// If you include it, you have to make a similar routine for
    !// convective scavenging.
    !//
    !// Ole Amund Sovde, July 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR, IDBLK, JDBLK
    use cmn_ctm, only: AREAXY, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: DT
    real(r8), intent(in) :: BTT(LPAR,NPAR,IDBLK,JDBLK), &
                            BTTBCK(LPAR,NPAR,IDBLK,JDBLK)
    !// Locals
    integer :: I,II,J,JJ, N, L, TRNR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='seasaltbudget_lscav'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Not use yet! - STOP'
    stop

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Diagnose sea salt removed by large scale wet scavenging.
          do N = 1, Nsaltbins
             TRNR = salt_trsp_idx(Nsaltbins)
             do L = 1, LPAR
                !// Accumulate scavenged [kg]
                salt_lscv(N,I,J) = salt_lscv(N,I,J) &
                     + (BTTBCK(L,TRNR,II,JJ) - BTT(L,TRNR,II,JJ))
             end do
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine seasaltbudget_lscav
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine vdep_n_coarse(losscoarse,DZ,I,J)
    !// --------------------------------------------------------------------
    !// Set drydep velocities of NH4coarse and NO3coarse based on sea salt
    !// fall velocities.
    !//
    !// Transport numbers for NH4coarse and NO3coarse are checked by
    !// salt_init.
    !//
    !// Ole Amund Sovde, March 2016, September 2009
    !// --------------------------------------------------------------------
    use cmn_sfc, only: VDEP
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: losscoarse, DZ
    integer, intent(in) :: I, J
    !// Locals
    integer :: TRSPNR
    !// --------------------------------------------------------------------

    !// NH4coarse (63)
    TRSPNR = trsp_idx(63)
    !// From 1/s to m/s
    VDEP(TRSPNR,I,J) = losscoarse * DZ

    !// NO3coarse (65)
    TRSPNR = trsp_idx(65)
    !// From 1/s to m/s
    VDEP(TRSPNR,I,J) = losscoarse * DZ

    !// --------------------------------------------------------------------
  end subroutine vdep_n_coarse
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_emis(NMET, NOPS, NSUB, CCYC, MP)
    !// --------------------------------------------------------------------
    !// Calculate sea salt emissions (kg/s) for this MP block, after
    !// meteorology has been updated.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_ctm, only: PLAND, AREAXY, &
       MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_met, only: SFT, SHF, USTR, SFU, SFV, CI
    use cmn_sfc, only: LSMASK
    use seasaltprod, only: seasalt_production, &
         seasalt_production_martensson03, seasalt_production_gantt15, &
         seasalt_production_witek16
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NMET, NOPS, NSUB, CCYC, MP

    !// Locals
    real(r8) :: &
         prodseasalt, &        !// production / flux [kg/s]
         rbinStart, rbinEnd, & !// start and end of radius bin [um]
         wind10, &  !// Wind speed at 10m
         AXY, &     !// Area [m2]
         fsaltwater !// Grid fraction covered by salt water [0:1]

    !// Indices
    integer :: I,J,II,JJ, N, UPDATEEMIS
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='seasalt_emis'
    !// --------------------------------------------------------------------

    !// Only update sea salt if meteorology is updated.
    UPDATEEMIS = NOPS * NSUB * CCYC
    if (UPDATEEMIS .ne. 1) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Grid fraction of salt water ocean
          fsaltwater = 1._r8 - PLAND(I,J)
          if (fsaltwater .lt. 1.e-10_r8) fsaltwater = 0._r8
          !// Check land-sea mask if this could be freshwater.
          if (LSMASK(I,J,3) .gt. 0._r8) then
             !// Scale down fsaltwater to match fraction of ocean to
             !// total water areas
             fsaltwater = fsaltwater &
                  * LSMASK(I,J,1) / (LSMASK(I,J,1) + LSMASK(I,J,3))
          end if

          !// Sea ice adjustment: Assume fraction (i.e. production)
          !// proportional to (1-CI)
          fsaltwater = fsaltwater * (1._r8 - CI(I,J))

          !// Zero production if no sea water
          if (fsaltwater .le. 0._r8) then
             seasalt_flux(:,II,JJ,MP) = 0._r8
             !// Go to next (I,J)
             cycle
          end if


          !// Wind speed at 10m
          wind10 = sqrt(SFU(I,J)**2 + SFV(I,J)**2)

          !// Grid box area
          AXY = AREAXY(I,J)


          !// Find sea salt flux (kg/s) for each bin
          do N = 1, NPAR_SALT

             rbinStart = rsalt80um(N)
             rbinEnd   = rsalt80um(N+1)

             if (SeaSaltScheme .eq. 1) then
                !// Default is calling Monahan86 production
                call seasalt_production(wind10, fsaltwater, AXY, &
                     rbinStart, rbinEnd, rhosaltdry, I, J, prodseasalt)
             else if (SeaSaltScheme .eq. 2) then
                !// Mårtensson et al, 2003
                call seasalt_production_martensson03(wind10, fsaltwater, AXY, &
                     rbinStart, rbinEnd, rhosaltdry, SFT(I,J), I, J, &
                     prodseasalt)
             else if (SeaSaltScheme .eq. 3) then
                !// Gantt et al, 2015
                call seasalt_production_gantt15(wind10, fsaltwater, AXY, &
                     rbinStart, rbinEnd, rhosaltdry, SFT(I,J), I, J, &
                     prodseasalt)
             else if (SeaSaltScheme .eq. 4) then
                !// Witek et al, 2016
                call seasalt_production_witek16(wind10, fsaltwater, AXY, &
                     rbinStart, rbinEnd, rhosaltdry, SFT(I,J), I, J, &
                     prodseasalt)
             else
                !// Not defined
                write(6,'(a,i3)') f90file//':'//subr// &
                     ': production scheme not defined', SeaSaltScheme
                stop
             end if

             !// Save production [kg/s] in array
             seasalt_flux(N,II,JJ,MP) = prodseasalt

          end do !// do N = 1, NPAR_SALT

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine seasalt_emis
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emissions_seasalt_total()
    !// --------------------------------------------------------------------
    !// Print out current total sea salt flux as [Tg/yr].
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_parameters, only: secYear
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    write(6,'(a,es12.3)') 'Instant sea salt from ocean [Tg/yr]: ',&
         sum(seasalt_flux) * secYear * 1.e-9_r8
    !// --------------------------------------------------------------------
  end subroutine emissions_seasalt_total
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module seasalt
!//=========================================================================

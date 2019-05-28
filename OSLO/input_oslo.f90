!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Further initialization of CTM.
!//=========================================================================
module input_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: input_oslo
  !// DESCRIPTION: Reads input to initialize CTM beyond routine INPUT.
  !//
  !// Contains:
  !//   subroutine init_oslo
  !//   subroutine info_qssa
  !//   subroutine read_tracer_list
  !//   subroutine tracer_specific_input
  !//   subroutine set_resultdir
  !//   subroutine read_lsmask
  !//
  !// Ole Amund Sovde, September 2008
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'input_oslo.f90'
  !// ----------------------------------------------------------------------
  public
  private info_qssa, tracer_specific_input, &
       set_resultdir, read_lsmask
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine clear_oslo_variables()
    !// --------------------------------------------------------------------
    !// Clear variables used for chemistry.
    !//
    !// Ole Amund Sovde, November 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: NOTRPAR
    use cmn_chem, only: TNAME, TMASS, TMOLMIX2MASSMIX, TMASSMIX2MOLMIX
    use cmn_oslo, only: trsp_idx, chem_idx, Xtrsp_idx, Xchem_idx, &
         XTMASS, XTNAME, XTMASSMIX2MOLMIX, XTMOLMIX2MASSMIX
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Transport index array
    trsp_idx(:)   = -99
    !// Reverse transport index array
    chem_idx(:) = -99
    !// Transported tracer names
    TNAME(:) = ' '
    !// Transported tracer mass
    TMASS(:) = 0._r8
    !// Conversion factors
    TMOLMIX2MASSMIX(:) = 0._r8
    TMASSMIX2MOLMIX(:) = 0._r8

    if (NOTRPAR .gt. 0) then
       !// NO-Transport index array
       Xtrsp_idx(:) = -99
       !// Reverse NO-transport index array
       Xchem_idx(:) = -99
       !// Non-transported tracer names
       XTNAME(:) = ' '
       !// Non-transported tracer mass
       XTMASS(:) = 0._r8
       !// Conversion factors
       XTMOLMIX2MASSMIX(:) = 0._r8
       XTMASSMIX2MOLMIX(:) = 0._r8
    end if
    !// --------------------------------------------------------------------
  end subroutine clear_oslo_variables
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine init_oslo(NDAY, DTCHM, CHMCYCLES)
    !// --------------------------------------------------------------------
    !// Initialize oslo chemistry. Tracer related data is read, and constant
    !// chemical reaction rates are set.
    !// Dry deposition parameters are read and some diagnostics are also
    !// initialized.
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LOSLOCTROP, LOSLOCSTRAT, LSULPHUR, LBCOC, &
         LSALT, LNITRATE, LDUST, LSOA, LEMISDEP_INCHEM, &
         TRACER_ID_MAX, IPARW, LPARW, LPAR, IPAR, JPAR
    use cmn_ctm, only: AREAXY, AREAG
    use cmn_chem, only: TNAME
    use cmn_sfc, only: ZPDVT_C3, LANDUSE_IDX
    use cmn_chem, only: INFILE_WET, INFILE_DRY, E2DS, E3DSNEW
    use fastjx, only: PHOT_IN
    use omp, only: get_all_mpind
    use scavenging_largescale_uci, only: WETSET_CTM3
    use scavenging_drydep_uci, only: DRYSET_CTM3
    !// --------------------------------------------------------------------
    use cmn_oslo, only: trsp_idx, chem_idx, Xtrsp_idx, Xchem_idx, XTNAME, &
         CH4FIELD, O3_UPBND
    use bcoc_oslo, only: bcoc_init
    use chem_oslo_rates, only: TCRATE_CONST2
    use ch4routines, only: update_ch4surface
    use sulphur_oslo, only: TCRATE_CONST_S, DMSseaconc
    use psc_microphysics, only: set_psc_constants
    use drydeposition_oslo, only: drydepinit
    use diagnostics_general, only: diag_ohch4n2o_init, init_chembud
    use diagnostics_scavenging, only: scav_diag_init
    use seasalt, only: seasalt_init
    use dust_oslo, only: dust_init
    use fallingaerosols, only: readkasten1968
    use nitrate, only: nitrate_init
    use soa_oslo, only: soa_init
    use aerosols2fastjx, only: initialize_tropaerosols4fjx
    use emissions_oslo, only: emis_input
    use emissions_ocean, only: emissions_ocean_organiccarbon_init
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, CHMCYCLES
    real(r8), intent(in) :: DTCHM     ! chemistry time step in seconds
    !// Locals
    integer :: J
    real(r8) :: RDUM
    !// ------------------------------------------------------------------
    character(len=*), parameter :: subr = 'init_oslo'
    !// ------------------------------------------------------------------

    !// Initialize grid mapping from I,J to II,JJ,MP
    call get_all_mpind()

    !// Emission arrays
    E2DS(:,:,:,:) = 0._r8
    E3DSNEW(:,:,:,:) = 0._r8
    DMSseaconc(:,:,:) = 0._r8

    !// Initialize tracer info
    call tracer_specific_input()

    !// Load CH4 field; can specify year
    if (trsp_idx(46).gt. 0) Call update_ch4surface()


    !// Vegetation is read in subroutine input (initialize.f90)

    !// Read additional land-sea mask to separate ocean/lakes
    call read_lsmask()

    !// Set chemical constant rates (chem_oslo_rates.f90)
    call TCRATE_CONST2()

    !// Initialize/set chemical constant rates for sulphur reactions
    !// (sulphur.f90)
    call TCRATE_CONST_S(LSULPHUR)

    !// Henry coefficients are calculated on the fly.


    !//  Stratospheric initializations
    if (LOSLOCSTRAT) then
       !// Regular rates dependent on temperature or constants
       !// are set in TCRATE_CONST2, covering both trop and strat.

       !// PSC parameters (psc_microphysics.f90)
       call set_psc_constants()

    else if (LOSLOCTROP) then
       !// Initialize stuff for only tropospheric chemistry
       !// Set the O3 flux at the model top for L40 data
       O3_UPBND(:) = 0._r8
       if (LPARW .eq. 40) then
          !// set ozone flux to 450Tg/year = 4.5E11 Kg/year, constant flux !
          do J = 1, JPAR
             O3_UPBND(J) = 4.5e11_r8/(365._r8*86400._r8)*AREAXY(1,J)/AREAG
          end do
       end if
    end if

    !// Initialize photochemistry
    call PHOT_IN

    !// Set wet scavenging parameters
    call WETSET_CTM3 (INFILE_WET)

    !// UCI dry deposition (will be overwritten by Oslo treatment)
    call DRYSET_CTM3 (INFILE_DRY)
    !// CTM3 dry deposition rates (used by setdrydep in drydeposition_oslo.f90)
    call drydepinit()


    !// JV arrays to account for tropospheric aerosols
    call initialize_tropaerosols4fjx()

    !// Read Kasten (1968) table for gravitational settling
    !// Not really necessary if fallingaerosols.f90 is not used.
    call readkasten1968()

    !// Global variables for OH and lifetimes
    call diag_ohch4n2o_init()

    !// Initialize chemistry budget arrays
    call init_chembud()

    !// Read emissions with new routine for the CTM3
    call emis_input()

    !// BCOC initializations
    if (LBCOC) call bcoc_init()

    !// SALT initializations
    if (LSALT) call seasalt_init()

    !// Organic matter from ocean initializations (needs salt)
    if (LBCOC) call emissions_ocean_organiccarbon_init()

    !// DUST initializations (needs emis_input to be carried out)
    if (LDUST) call dust_init(NDAY,DTCHM/real(CHMCYCLES, r8))

    !// NITRATE initializations
    if (LNITRATE) call nitrate_init(trsp_idx, TRACER_ID_MAX, LSULPHUR, LSALT)

    !// SOA initializations
    if (LSOA) call soa_init()

    !// Diagnose scavenging (wet/dry) - init
    call scav_diag_init(1)

    !// Look for resultdir_file
    call set_resultdir()

    !// displacement height (m) matching MODIS (taken from DEAD and adjusted
    !// to CTM3 vegetation)
    ZPDVT_C3(:) = 0._r8
    if (LANDUSE_IDX .eq. 2) then
       ZPDVT_C3(1:17) = (/11.39_r8, 23.45_r8, 9.38_r8, 13.40_r8, &
            12.06_r8, 0.34_r8, 0.34_r8, 0.34_r8, 0.34_r8, 0.34_r8, 0.34_r8, &
            0.34_r8, 0.34_r8, 0.34_r8, 0._r8, 0._r8, 0._r8/)
    else if (LANDUSE_IDX .eq. 3) then
       !// CLM4-PFT
       !// Typically 2/3 of the height of obstacles.
       !// Using CLM4 ztop:
       ZPDVT_C3(1:17) = 2._r8/3._r8 &
            * (/17.0_r8, 17.0_r8, 14.0_r8, 35.0_r8, 35.0_r8, 18.0_r8, &
                20.0_r8, 20.0_r8,  0.5_r8,  0.5_r8,  0.5_r8,  0.5_r8, &
                 0.5_r8,  0.5_r8,  0.5_r8,  0.5_r8,  0.0_r8/)
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': No such LANDUSE_IDX available',LANDUSE_IDX
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine init_oslo
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine info_qssa()
    !// --------------------------------------------------------------------
    !// Print out info about oslo chemistry.
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    write(6,'(a)') 'Basic QSSA chemistry timesteps are:'
    write(6,'(a)') '   time step    = DTCHM'
    write(6,'(a)') '   steady state = DTCHM/10'
    write(6,'(a)') '   euler        = DTCHM*10'
    write(6,'(a)') '   test         = DTCHM'

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    !// --------------------------------------------------------------------
  end subroutine info_qssa
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine read_tracer_list(fname)
    !// --------------------------------------------------------------------
    !// Read tracer list.
    !//
    !// Ole Amund Sovde, November 2015, September 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LOSLOCTROP, LOSLOCSTRAT, LSULPHUR, LBCOC, &
         LSALT, LNITRATE, LDUST, LEMISDEP_INCHEM, LE90, LLINOZ, &
         NPAR, NOTRPAR, IDGRD, JDGRD, NDGRD, TNAMELEN
    use cmn_ctm, only: LCONT, START_AVG, NTM
    use cmn_chem, only: TNAME, TMASS, N_LZ, Ne90, N_STE, TMASSMIX2MOLMIX, &
         TMOLMIX2MASSMIX
    use cmn_met, only: HnativeRES,VRES,VRESW
    use cmn_parameters, only: M_AIR
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx, chem_idx, Xtrsp_idx,Xchem_idx, &
         XTMASS,XTMASSMIX2MOLMIX,XTMOLMIX2MASSMIX, XTNAME
    use psc_microphysics, only: LPSC, LAEROSOL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: FNAME

    !// Locals
    integer :: TRACER_ID, IERR, IDUM, fnr, &
         INP, INP2, I
    real(r8) :: MASS
    character(len=TNAMELEN) :: SNAM
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_tracer_list'
    !// ------------------------------------------------------------------

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    write(6,'(a)') f90file//':'//subr//': Reading tracer list: '//trim(FNAME)
    write(6,'(a)') '- - - Running with Oslo chemistry - - -'


    !// Initialize
    INP      = 0
    INP2     = 0
    TNAME(:) = ' '
    XTNAME(:)= ' '
    !// Linoz/E90/STEFLUX
    N_LZ = 0
    Ne90 = 0
    N_STE= 0

    !// Get a free file number
    fnr = get_free_fileid()
    open(fnr,file=FNAME,status='old',form ='formatted')
    read(fnr,*) !// read header
    read(fnr,*) !// read dummy line

    if (LOSLOCSTRAT) then
       write(6,'(a,l1)') 'PSC heterogeneous chemistry: ', LPSC
       write(6,'(a,l1)') 'Aerosol heterogeneous chemistry: ', LAEROSOL
    end if

    !// No need to read number of tracers here, when defined in params.F

    write(6,'(a,i3)') 'Transported tracers:     ', NPAR
    write(6,'(a,i3)') 'Non-transported tracers: ', NOTRPAR

    write(6,'(a)') 'Horizontal resolution: '//trim(HnativeRES)
    write(6,'(a,3i3)') '  window multiplier:   ',IDGRD,JDGRD,NDGRD
    write(6,'(a)') 'Vertical resolution:   '//VRES//'/'//VRESW

    do i=1,4
       read(fnr,*) !// read dummy line
    end do



    !// First read transported species
    !// ------------------------------------------------------------------
    write(6,'(a)') '- - - - - - - - - - - - - - - - - - - - - - - -'
    write(6,'(A)') ' NR  ID NAME        Mw                         '
    write(6,'(a)') '- - - Transported species - - - - - - - - - - -'
    do I = 1,NPAR
       read(fnr,*,iostat=ierr) TRACER_ID,SNAM,MASS
       if (ierr .ne. 0) then
          write(6,'(a)') f90file//':'//subr//': Problems reading '//trim(fname)
          write(6,'(a,i5)') '* Check NPAR vs tracer list, NPAR: ',NPAR
          write(6,'(a,2i5)') '* Trying to read tracer #/ID ',I,TRACER_ID
          write(6,'(a,i5)') '* NPAR is larger than listed:',NPAR
          stop
       end if
       INP = INP + 1
       trsp_idx(TRACER_ID) = I  ! Mapping from tracer id to transport number
       chem_idx(I) = TRACER_ID  ! Mapping from trsp.nr. (MTC) to tracer id
       TMASS(I) = MASS          ! M_tracer
       TMASSMIX2MOLMIX(I) = M_AIR / MASS  ! M_air/M_tracer
       TMOLMIX2MASSMIX(I) = MASS / M_AIR  ! M_tracer/M_air
       TNAME(I) = SNAM          ! Tracer name
       !// Check if LINOZ is to be used
       if (trim(SNAM) .eq. 'O3LINOZ') then
          N_LZ = I             !// Linoz tracer
          !N_STE = I            !// STE tracer (set to N_LZ if LINOZ)
          if (.not.LLINOZ) then
             write(6,'(a)') f90file//':'//subr// &
                  ': LINOZ tracer in '//trim(fname)//' but not in Makefile'
             stop
          end if
       end if
       !// Check if E90-tracer is to be used
       if (trim(SNAM) .eq. 'E090') then
          Ne90 = I             !// E90-tracer
          if (.not.Le90) then
             write(6,'(a)') f90file//':'//subr// &
                  'E90 tracer in '//trim(fname)//' but not in Makefile'
             stop
          end if
       end if

       write(6,'(I3,1x,I3,1x,A,1x,f7.3)') I,TRACER_ID,SNAM,MASS
    end do

    !// STEFLUX: Set STE tracer number if not LINOZ
    !if (N_LZ .le. 0 .and. trsp_idx(1).gt.0) N_STE = trsp_idx(1)
    !// Will use O3 as STE tracer, 
    if (trsp_idx(1).gt.0) N_STE = trsp_idx(1)

    !// Check input file for correct number of transported species
    read(fnr,'(i3)',iostat=ierr) IDUM
    if (idum .ne. -2) then
       write(6,'(a)') f90file//':'//subr//': Problems reading '//trim(fname)
       write(6,'(a)') '* Too many transported species in tracer list'
       write(6,'(a,i5)') '* Check NPAR vs tracer list, NPAR: ',NPAR
       stop
    end if

    !// Read non-transported species
    !// ------------------------------------------------------------------
    write(6,'(a)') '- - - Non-transported species - - - - - - - - -'
    do I = 1, NOTRPAR
       read(fnr,*,iostat=ierr) TRACER_ID,SNAM,MASS
       if (ierr .ne. 0) then
          write(6,'(a)') f90file//':'//subr//': Problems reading '//trim(fname)
          write(6,'(a,i5)') '* Check NOTRPAR vs tracer list ',NOTRPAR
          write(6,'(a,2i5)') '* Trying to read tracer #/ID ',I,TRACER_ID
          write(6,'(a)') '* NOTRPAR is larger than listed!'
          stop
       end if
       INP2 = INP2 + 1
       Xtrsp_idx(TRACER_ID) = I ! Mapping from tracer id to NO-transport number
       Xchem_idx(I) = TRACER_ID ! Mapping from non-trsp. nr. (XMTC) to tracer id
       XTMASS(I) = MASS         ! M_tracer
       XTMASSMIX2MOLMIX(I) = M_AIR / MASS  ! M_air/M_tracer
       XTMOLMIX2MASSMIX(I) = MASS / M_AIR  ! M_tracer/M_air
       XTNAME(I) = SNAM         ! Tracer name
       write(6,'(I3,1x,I3,1x,A10,1x,f7.3)') I,TRACER_ID,SNAM,MASS
    end do
    write(6,'(a)') '- - - - - - - - - - - - - - - - - - - - - - - -'

    !// Check if all tracers are read
    read(fnr,'(i3)',iostat=ierr) IDUM
    if (idum .ne. -3) then
       write(6,'(a)') f90file//':'//subr//': Problems reading '//trim(fname)
       write(6,'(a)') '* Too many non-transported species in tracer list'
       write(6,'(a,i5)') '* Check NOTRPAR vs tracer list ',NOTRPAR
       stop
    end if

    !// Check number of species
    NTM = INP
    if (NTM .ne. NPAR) then
       write(6,'(a)') '* Mis-match: # transported tracers vs NPAR'
       write(6,'(a)') '* Tracer list: '//trim(fname)
       write(6,'(a,2i5)') '  NTM,NPAR: ',NTM,NPAR
       stop
    end if
    if (INP2 .ne. NOTRPAR) then
       write(6,'(a)') '* Mis-match: # non-transported tracers vs NOTRPAR'
       write(6,'(a)') '* Tracer list: '//trim(fname)
       write(6,'(a,2i5)') '  INP2,NOTRPAR: ',INP2,NOTRPAR
       stop
    end if

    !// Double check that all species are set
    do I = 1,NPAR
       TRACER_ID = chem_idx(I)
       if (TRACER_ID .le. 0) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': VERY VERY WRONG (N,ID)',I,chem_idx(I)
          stop
       end if
    end do
    do I = 1,NOTRPAR
       TRACER_ID = Xchem_idx(I)
       if (TRACER_ID .le. 0) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': VERY VERY WRONG (N,XID)',I,Xchem_idx(I)
          stop
       end if
    end do

    !// Chemistry diagnostics
    !// --------------------------------------------------------------------
    !// Here diagnostics for specific reactions was placed in UCI-CTM.
    !// This is not an option in CTM3.
    close(fnr)

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine read_tracer_list
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine tracer_specific_input()
    !// --------------------------------------------------------------------
    !// Read in the tracer-specific run instructions and some data sets.
    !// Based on UCI CHEM_IN.
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LOSLOCHEM, LOSLOCSTRAT, LOSLOCTROP, LBCOC, &
         LSALT, LDUST, LLINOZ, JPAR, LPAR, NPAR
    use cmn_ctm, only: ANAME, LCONT, JMON, IYEAR
    use cmn_chem, only: INFILE_T
    use steflux, only: chemflux_setup
    use utilities, only: ctmExitC
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'tracer_specific_input'
    !// --------------------------------------------------------------------

    !// Wrong place if not Oslo chemistry is chosen
    if (.not. LOSLOCHEM) call ctmExitC('LOSLOCHEM is not defined!')

    !// Name of AIR
    ANAME = 'AIR'

    !// Test Oslo chemistry switches
    if (LOSLOCSTRAT .and. .not.LOSLOCTROP) then
       write(6,'(a)') f90file//':'//subr// &
            ': Oslo chemistry does not allow only stratospheric chemistry!'
       write(6,'(a)') 'Check Makefile and parameters.'
       stop
    end if
    if (.not.(LOSLOCSTRAT .or. LOSLOCTROP)) then
       !// May use stand-alone modules...
       if ( .not. (LBCOC .or. LSALT .or. LDUST) ) then
          write(6,'(a)') '    LOSLOCSTRAT/LOSLOCTROP not set!'
          write(6,'(a)') '    And no stand-alone modules!'
       end if
    end if

    !// Print info about chemistry/qssa
    call info_qssa()

    !// In case of Linoz
    if (LLINOZ)  then
       call LNZ_INIT
       call LNZ_SET(JMON,IYEAR)
       call LNZ_SETO3
       !call TPAUSEG !// Will be set in pmain, using CTM3 routine
    else
       !// Needs 3 parameters from Linoz files
       !// O3TPALZ: Mixing ratio for finding e90 tropopause
       !// O3iso1/O3iso2: O3 isopleths for flux calculations
       call chemflux_setup()
    end if


    write(6,'(a,l2)') ' LCONTinuation =',LCONT

    !// --------------------------------------------------------------------
  end subroutine tracer_specific_input
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_resultdir()
    !// --------------------------------------------------------------------
    !// Reads file resultdir_file to get path to where to save result
    !// files. Default is './'.
    !//
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_oslo, only: RESULTDIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(100) ::  file_name
    integer :: N
    logical :: fex
    !// --------------------------------------------------------------------

    !// Set result dir to be used for saving diag. Default is './'
    file_name = 'resultdir_file'
    inquire(file=file_name,exist=fex)
    if (fex) then
       N = 20
       fex=.true.
       do while (fex)
          N = N + 1
          inquire(N,opened=fex)
       end do
       open(N,file=file_name,form='formatted')
       read(N,'(A)') RESULTDIR
       close(N)
       if (RESULTDIR .eq. '') RESULTDIR = '.'
       RESULTDIR = trim(RESULTDIR)//'/'
    else
       print*,'* set_resultdir: No resultdir_file'
       RESULTDIR = './'
    end if
    !// --------------------------------------------------------------------
  end subroutine set_resultdir
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_lsmask()
    !// --------------------------------------------------------------------
    !// Read land-sea mask. The field values are:
    !// 0. Ocean
    !// 1. Land
    !// 2. Lake
    !// 3. Small Islands
    !// 4. Ice Shelf
    !// 
    !// Ole Amund Sovde, October 2014
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, AREAXY
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use netcdf
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    use cmn_sfc, only: LSMASK
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=80) :: FILENAME
    integer           :: nLon, nLat, I,J, N
    real(r8), dimension(IPAR,JPAR) :: R8CTM
    real(r8), allocatable, dimension(:) :: inLon, inLat, xbedge,ybedge,xybox
    real(r8), allocatable, dimension(:,:) :: IXY
    real(r8), allocatable, dimension(:,:,:) :: R8XYZ
    real(r8), allocatable, dimension(:,:) :: R8XY
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_lsmask'
    !// ------------------------------------------------------------------

    filename = 'Indata_CTM3/landsea.nc'

    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( filename, 'lon',  inLon  )
    call get_netcdf_var_1d( filename, 'lat',  inLat  )
    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )

    allocate(xbedge(nLon+1), ybedge(nLat+1), xybox(nLat), &
         IXY(nLon,nLat), R8XYZ(nLon,nLat,5), R8XY(nLon,nLat))

    !// Read LSMASK
    call get_netcdf_var_2d( filename, 'LSMASK', IXY, nlon,nlat )

    !// Set map fractions
    R8XYZ(:,:,:) = 0._r8
    do J = 1, nlat
       do I = 1, nlon
          if (IXY(I,J) .eq. 0) then
             R8XYZ(I,J,1) = 1._r8
          else if (IXY(I,J) .eq. 1) then
             R8XYZ(I,J,2) = 1._r8
          else if (IXY(I,J) .eq. 2) then
             R8XYZ(I,J,3) = 1._r8
          else if (IXY(I,J) .eq. 3) then
             R8XYZ(I,J,4) = 1._r8
          else if (IXY(I,J) .eq. 4) then
             R8XYZ(I,J,5) = 1._r8
          end if
       end do
    end do

    !// Interpolate to model grid
    !// Set up grid for indata
    do I = 1,nLon
       XBEDGE(I) = inLon(I) - 0.5_r8*abs(inLon(2)-inLon(1))
    end do
    XBEDGE(nLon+1) = XBEDGE(1) + 360._r8
    do J = 1,nLat
       YBEDGE(J) = inLat(J) - 0.5_r8*abs(inLat(2)-inLat(1))
    end do
    YBEDGE(nLat+1) = 90._r8
    do J=1,nLat
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * abs(sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    do N = 1, 5
       do J = 1, nLat
          do I = 1, nLon
             R8XY(I,J) = R8XYZ(I,J,N)*XYBOX(J)
          end do
       end do
       call E_GRID(R8XY,XBEDGE,YBEDGE,nLon,nLat, &
            R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)
       do J = 1,JPAR
          do I = 1,IPAR
             LSMASK(I,J,N) = R8CTM(I,J)/AREAXY(I,J)
          end do
       end do
    end do

    !// Deallocate all local variables
    deallocate(inLon, inLat, xbedge,ybedge, IXY, R8XYZ, R8XY)

    write(6,'(a)') f90file//':'//subr//': land-sea mask is read'
    !// --------------------------------------------------------------------
  end subroutine read_lsmask
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module input_oslo
!//=========================================================================

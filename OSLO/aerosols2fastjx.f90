!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Further initialization of CTM.
!//=========================================================================
module aerosols2fastjx
  !// ----------------------------------------------------------------------
  !// MODULE: aerosols2fastjx
  !// DESCRIPTION: Allows CTM3 tropospheric aerosols to be used in Fast-JX.
  !//
  !// In this module we define all the aerosols to include in Fast-JX
  !// calculations. The aerosols can be retrieved from the tracer array
  !// or from climatologies.
  !//
  !// IMPORTANT:
  !// Fast-JX does not use the tracer fields (STT) directly:
  !// Aerosols with the same optical properties will be summed up to
  !// give total aerosol paths for each type. Fast-JX will then use
  !// the total path.
  !//
  !// A climatology of aerosols produced by the Oslo CTM2, is used as a
  !// basis, and aerosols not included in the simulation can be taken from
  !// there.
  !//
  !// Contains
  !//   - subroutine initialize_tropaerosols
  !//   - subroutine set_aer4fjx
  !//   - subroutine get_tropaerosols
  !//   - subroutine update_tropaerosols
  !//   - subroutine set_aer4fjx_ctm2
  !//
  !// Amund Sovde, May-July 2010, June 2012
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR
  use cmn_ctm, only: JDATE, XDEDG, YDEDG
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Switch to turn on/off aerosol effect on Fast-JX
  logical,parameter :: LJV_AEROSOL = .true.

  !// Number of aerosols to be included in Fast-JX
  integer, parameter :: NUM_AER = 1
!  integer, parameter :: NUM_AER = 51

  !// Tracer IDs for the included aerosols
  !// These will be included from STT if available, or from climatology
  !// APM defined by CLIM_ID below.
  !//   Possible tracers
  !//     73 Sulphate
  !//    211-218 Mineral dust
  !//    240-247 BCOC
  !//     64 NO3fine
  !// List all possible in a comment here?
  integer, dimension(NUM_AER), parameter :: AER_TRACER_ID = (/-1/)
!  integer, dimension(NUM_AER), parameter :: AER_TRACER_ID = &
!       (/ 73, &                                     !// Sulphate
!          64, &                                     !// NO3fine
!          201, 202, 203, 204, 205, 206, 207, 208, & !// SALT
!          211, 212, 213, 214, 215, 216, 217, 218, & !// DUST
!          230, 231, 232, 233, 234, 235, &           !// OM
!          240, 241, 242, 243, 244, 245, &           !// BC
!          165, 166, 167, 168, 169, 170, 171, 172, 173, 174, &     !// SOA
!          175, 176, 177, 178, 179, 182, 183, 188, 189, 190, 191 /)!// SOA

  !// Type of aerosol: Specific for each tracer in AER_TRACER_ID, and must
  !//                  match entry in FJX_scat.dta
  integer, dimension(NUM_AER) :: AERTYPE


  !// Number of aerosol TYPES (some aerosols have the same optical properties).
  !// Aerosols will be treated according to optical types, for which PATHs
  !// will be summed up and treated with the same optical properties.
  integer, parameter :: MAX_AER_TYPES = 25

  !// Keep track of the actual number of types
  integer :: NUM_AER_TYPES

  !// Aerosol types in use
  integer, dimension(MAX_AER_TYPES) :: AER_TYPES_USED

  !// Paths [g/m2] of the aerosol types (set from STT/BTT or climatology)
  real(r8), dimension(LPAR+1,MAX_AER_TYPES,IPAR,JPAR) :: AERPATH

  !// Aerosols to be read from climatology
  integer, dimension(NUM_AER) :: AERFROMCLIM


  !// Climatology variables
  !// ----------------------------------------------------------------------
  !// Number of aerosol types to be fetched from climatology.
  !// The climatology may contain more aerosols than this.
  integer, parameter :: NUM_CLIM = 2

  !// List aerosol IDs to be read from file
  integer,dimension(NUM_CLIM),parameter :: CLIM_ID = &
       (/ 254, 255 /)

  !// Save data for two months
  real(r8), dimension(2,LPAR,NUM_CLIM,IPAR,JPAR) :: APM

  !// For time interpolation between the two months, updated at 00UTC unless
  !// you change the code.
  real(r8) :: aerdtfrac

  !// ----------------------------------------------------------------------
  !// All variables are to be saved. This should be default for global
  !// variables, although Fortran claims it may not be.
  save
  !// ----------------------------------------------------------------------
  !// All is private
  private
  !// except
  public update_tropaerosols4fjx, set_aer4fjx, set_aer4fjx_ctm2, &
       initialize_tropaerosols4fjx, NUM_AER_TYPES, AERPATH, &
       AER_TYPES_USED, LJV_AEROSOL
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine initialize_tropaerosols4fjx()
    !// --------------------------------------------------------------------
    !// Initialize arrays for tropopsheric aerosols to be used in Fast-JX
    !// scattering/absorption.
    !//
    !// Amund Sovde, July 2010
    !// --------------------------------------------------------------------
    use cmn_chem, only: TNAME
    use cmn_oslo, only: chem_idx, trsp_idx
    use bcoc_oslo, only: BC_JVAER_TYPE, OM_JVAER_TYPE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: NT, NA, NC, TRID, AERINC
    logical :: foundtype
    !// --------------------------------------------------------------------

    !// Initialize
    AERPATH(:,:,:,:) = 0._r8
    AERTYPE(:) = 0
    AERFROMCLIM(:) = 0

    NUM_AER_TYPES = 0
    AER_TYPES_USED(:) = 0

    write(*,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    write(*,'(a)') '* Aerosols to be used in Fast-JX (in aerosols2fastjx.f90):'
    write(*,'(a)') '  Initializing...'

    !// Count number of included aerosols
    AERINC = 0

    if (.not.LJV_AEROSOL) then
       !// Check tracers for aerosols
       write(*,'(a)') '  No aerosols will be treated in Fast-JX!'
       write(*,'(a)') '  aerosols2fastjx: LJV_AEROSOL = .false.'
       write(*,'(a71)') '--------------------------------------------'// &
            '---------------------------'
       return
    end if

    !// Check all transported tracers
    do NT = 1, NPAR

       !// Tracer ID
       TRID = chem_idx(NT)

       !// Check whether tracer should be used in Fast-JX
       do NA = 1, NUM_AER

          if (TRID .eq. AER_TRACER_ID(NA)) then
             !// Aerosol is included in the simulation
             !// Set aerosol type, hardcoded.
             select case(TRID)
             case (64)
                !// 64: NOfine (depends on relative humidity: UMaer)
                AERTYPE(NA) = -44
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (73)
                !// 73: Sulphate (depends on relative humidity: UMaer)
                AERTYPE(NA) = -35
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (165:179, 182:183, 188:191)
                !// SOA (use OCFF)
                AERTYPE(NA) = -34
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (201:208)
                !// 201:208: SALT
                if (TNAME(NT).eq.'SALT01') AERTYPE(NA) = -36
                if (TNAME(NT).eq.'SALT02') AERTYPE(NA) = -37
                if (TNAME(NT).eq.'SALT03') AERTYPE(NA) = -38
                if (TNAME(NT).eq.'SALT04') AERTYPE(NA) = -39
                if (TNAME(NT).eq.'SALT05') AERTYPE(NA) = -40
                if (TNAME(NT).eq.'SALT06') AERTYPE(NA) = -41
                if (TNAME(NT).eq.'SALT07') AERTYPE(NA) = -42
                if (TNAME(NT).eq.'SALT08') AERTYPE(NA) = -43
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (211:218)
                !// 211:218: DUST
                if (TNAME(NT).eq.'DUST01') AERTYPE(NA) = 34
                if (TNAME(NT).eq.'DUST02') AERTYPE(NA) = 35
                if (TNAME(NT).eq.'DUST03') AERTYPE(NA) = 36
                if (TNAME(NT).eq.'DUST04') AERTYPE(NA) = 37
                if (TNAME(NT).eq.'DUST05') AERTYPE(NA) = 38
                if (TNAME(NT).eq.'DUST06') AERTYPE(NA) = 39
                if (TNAME(NT).eq.'DUST07') AERTYPE(NA) = 40
                if (TNAME(NT).eq.'DUST08') AERTYPE(NA) = 41
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (230:235)
                !// 230:235: OM
                AERTYPE(NA) = OM_JVAER_TYPE(TRID-229)
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (240:245)
                !// 240:245: BC
                AERTYPE(NA) = BC_JVAER_TYPE(TRID-239)
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case (254)
                !// 254: BC1 from M7 simulation
                AERTYPE(NA) = 3
                write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID, &
                     ': STT will be used, TYPE: ',AERTYPE(NA)
                AERINC = AERINC + 1
             case default
                !// do nothing
                write(*,'(4x,a,i3,a)') 'Aerosol number: ',TRID, &
                ': not defined in initialize_tropaerosols in '// &
                'aerosols2fastjx.f90'
             end select
             !// Go to next tracer
             exit
          end if !// if (TRID .eq. AER_TRACER_ID(NA)) then

       end do !// do NA = 1, NUM_AER

    end do !// do N = 1, NPAR


    !// Check which tracers to read from climatology
    !// --------------------------------------------------------------------
    do NA = 1, NUM_AER

       !// If aerosol is included in STT, skip to next tracer
       if (.not.(AERTYPE(NA) .eq. 0)) cycle

       !// otherwise, we check the climatology

       !// Tracer ID to check for
       TRID = AER_TRACER_ID(NA)

       if (NUM_AER.eq.1 .and. TRID.eq.-1) then
          !// CTM2 climatology of BC
          AERTYPE(1) = 3
       else
          !// Is the tracer ID included in the climatology? If so, use it
          !// for the tracer in Fast-JX.
          do NC = 1, NUM_CLIM

             if (TRID .eq. CLIM_ID(NC)) then
                !// Aerosol is included in the climatology
                !// Use index in climatology as flag (positive number)
                AERFROMCLIM(NA) = NC
                !// Set aerosol type, hardcoded.
                select case(TRID)
                case (254)
                   !// 254: BC1 from M7 simulation
                   AERTYPE(NA) = 3
                   write(*,'(4x,a,i3,a,i3)') 'Aerosol number: ',TRID,&
                        ': Climatology will be used, TYPE: ',AERTYPE(NA)
                case default
                   !// do nothing
                   write(*,'(4x,a,i3,a)') 'Skipping aerosol ',TRID, &
                        ': in J-calculations.'
                end select
                !// Go to next climatology aerosol
                exit
             end if
          end do !// do NC = 1, NUM_CLIM
       end if !// if (NUM_AER.eq.1 .and. TRID.eq.-1) then

    end do !// do NA = 1, NUM_AER


    !// Print out info on tracers
    !// --------------------------------------------------------------------
    write(*,'(a)') '  Aerosols to be used in Fast-JX are initialized as:'
    do NA = 1, NUM_AER
       TRID = AER_TRACER_ID(NA)
       if (TRID .le. 0) cycle

       write(*,'(4x,A,3x,I3,A,L1,A,L1)') &
            'Tracer ID / transported / from climatology: ',&
            AER_TRACER_ID(NA),' / ',trsp_idx(TRID).gt.0,' / ',AERFROMCLIM(NA).gt.0
    end do
    if (NUM_AER .eq. 1 .and. AER_TRACER_ID(1).eq.-1) then
       write(*,'(4x,A)') &
            'Climatology of BC will be used in FastJX'
    end if


    !// Find aerosol types based on optical properties
    !// --------------------------------------------------------------------
    AER_TYPES_USED(:) = 0
    !// Must check both STT tracers and the ones from climatology
    !// Set first type, if present
    if (AERTYPE(1).ne.0) then
       AER_TYPES_USED(1) = AERTYPE(1)
       NUM_AER_TYPES = 1
    else
       NUM_AER_TYPES = 0
    end if
    !// Loop through the rest
    do NA = 2, NUM_AER
       !// Only check nonzero types (i.e. either in STT or in climatology)
       if (AERTYPE(NA) .eq. 0) cycle
       !// Loop through AER_TYPES_USED to find new entry
       foundtype = .false.
       do NC = 1, NUM_AER_TYPES
          if (AER_TYPES_USED(NUM_AER_TYPES) .eq. AERTYPE(NA)) then
             foundtype = .true.
             exit
          end if
       end do
       if (.not.foundtype) then
          NUM_AER_TYPES = NUM_AER_TYPES + 1
          if (NUM_AER_TYPES .gt. MAX_AER_TYPES) then
             write(*,'(a)') '* NUM_AER_TYPES > MAX_AER_TYPES: Check aerosols2fastjx.f90'
             write(*,'(a,2i3)') '  NUM_AER_TYPES/MAX_AER_TYPES: ',NUM_AER_TYPES,MAX_AER_TYPES
             stop
          end if
          AER_TYPES_USED(NUM_AER_TYPES) = AERTYPE(NA)
       end if
    end do
    write(*,'(a)') '  The aerosols tracers will be grouped in types in Fast-JX'
    write(*,'(a,i3)') '  Aerosols types used in Fast-JX: ',NUM_AER_TYPES
    do NC = 1, NUM_AER_TYPES
       write(*,'(4x,A,2(1x,I3))') 'Aerosol type: ', NC, AER_TYPES_USED(NC)
    end do
    write(*,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine initialize_tropaerosols4fjx
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_aer4fjx(BTT,MP)
    !// --------------------------------------------------------------------
    !// Calculate tropopsheric aerosol paths for Fast-JX absorption and
    !// scattering. If an aerosol is transported (is in trsp_idx), the path is
    !// calculated from aerosol mass. If not, it is fetched from climatology
    !// if available.
    !// If climatology is not available, path and type are already zero, and
    !// will not be treated in Fast-JX.
    !//
    !// Sums up tracer masses when optical properties are equal.
    !//
    !// Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY
    use cmn_oslo, only: chem_idx, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in):: BTT

    !// Locals
    integer :: I,J,L,II,JJ
    integer :: NT, NA, NC, NUSED
    real(r8) :: ZAREA, PATH(LPAR)
    real(r8), parameter :: ZEROPATH = 1.e-20_r8
    !// --------------------------------------------------------------------

    !// Check for switch
    if (.not.LJV_AEROSOL) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Inverse grid box area [m^{-2}]
          ZAREA = 1._r8/AREAXY(I,J)

          !// Loop through AER_TYPES_USED and sum up for AERPATH
          do NUSED = 1, NUM_AER_TYPES

             !// Initialize path for this type
             PATH(:) = 0._r8

             !// Loop through aerosol tracers for Fast-JX
             do NA = 1, NUM_AER

                !// Skip non-included tracers
                if (AER_TRACER_ID(NA) .le. 0) cycle

                !// Sum up for this type?
                if (AER_TYPES_USED(NUSED) .ne. AERTYPE(NA)) cycle

                !// Get transport number
                NT = trsp_idx(AER_TRACER_ID(NA))
                !// Possible number in the saved climatology
                NC = AERFROMCLIM(NA)

                if (NT .gt. 0) then
                   !// Tracer is included in simulation. Convert
                   !// tracer mass to path [g/m2]. 1.d3 to change from kg->g
                   !AERPATH(1:LPAR,NA,I,J) = 1.e3_r8 * BTT(:,NT,II,JJ) * ZAREA
                   do L = 1, LPAR
                      !// Sum up
                      PATH(L) = PATH(L) + 1.e3_r8 * BTT(L,NT,II,JJ) * ZAREA
                   end do
                else if (NC .gt. 0) then
                   !// Get aerosol PATH from climatology.
                   !// Aerosol type was set in initialize_tropaerosols.
                   !// APM is already [g/m2].
                   do L = 1, LPAR
                      PATH(L) = PATH(L) + APM(1,L,NC,I,J)*aerdtfrac + &
                           APM(2,L,NC,I,J)*(1._r8 - aerdtfrac)
                   end do
                else
                   !// Not necessary to do anything; AERPATH and AERTYPE
                   !// are already zero.
                end if


             end do !// do NA = 1, NUM_AER

             !// Use calculated total PATH of aerosol type
             do L = 1, LPAR
                if (PATH(L) .le. ZEROPATH) then
                   !// Skip very small values
                   AERPATH(L,NUSED,I,J) = 0._r8
                else
                   AERPATH(L,NUSED,I,J) = PATH(L)
                end if
             end do !// do L = 1, LPAR

          end do !// do NUSED = 1, NUM_AER_TYPES

       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine set_aer4fjx
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_clim_tropaerosols4fjx(MONTH)
    !// --------------------------------------------------------------------
    !// Read tropospheric aerosols from file. Levels 1-32, which are the
    !// levels where L40 and L60 are roughly similar.
    !// Handles collapsed layers, but reads data from non-collapsed CTM levels.
    !//
    !// Amund Sovde, May 2010
    !// --------------------------------------------------------------------
    use netcdf
    use cmn_precision, only: r4
    use cmn_ctm, only: XLMMAP, LMMAP, AREAXY
    use regridding, only: E_GRID
    use ncutils, only: handle_error, handle_err
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    integer, intent(in) :: MONTH

    !// Locals
    integer :: I, J, L, LL, N, APAR, LMAP(LPAR+1)
    integer :: MON1, MON2

    !// File treatment
    character(len=80) :: filename
    real(r4), dimension(:,:,:,:),allocatable :: aer_in !// In-field
    real(r8), dimension(:,:),allocatable :: aer2d      !// For interpolation
    real(r8), dimension(:),allocatable :: &
         aer_xdedg, aer_ydedg !// For interpolation
    integer, dimension(:),allocatable :: aer_ids

    !// netcdf integers
    integer                  :: lon_dim_id   !// Id for longitude dimension
    integer                  :: lon_id       !// Id for variable longitude
    integer                  :: lat_dim_id   !// Id for latitude dimension
    integer                  :: lat_id       !// Id for latitude
    integer                  :: lev_dim_id   !// Id for level dimension
    integer                  :: lev_id       !// Id for level
    integer                  :: time_dim_id  !// Id for time dimension
    integer                  :: time_id      !// Id for time
    integer                  :: naer_dim_id  !// Id for num_aerosols dimension
    integer                  :: naer_id      !// Id for number of aerosols
    integer                  :: field_id     !// Variable id for field
    integer                  :: nlons        !// Longitudes in file
    integer                  :: nlats        !// Latitudes in file
    integer                  :: nlevs        !// Levels in file
    integer                  :: nsteps       !// Timesteps avaiable in file
    integer                  :: naer         !// Number of aerosols in file
    integer                  :: status       !// status of process (0=OK)
    integer                  :: ncid         !// file id 

    !// Name of the climatology field
    character(len=10) :: CLIM_NAME
    character(len=3) :: CNUM
    !// Are the tracers we want in the climatology?
    integer,dimension(NUM_CLIM) :: AERinCLIM

    !// For interpolation
    logical :: LINTERP
    real(r8), dimension(IPAR, JPAR) :: R8XY
    real(r8) :: maxcol1,maxcol2,meancol1,meancol2

    real(r8), parameter :: ZEROPATH = 1.e-20_r8
    !// --------------------------------------------------------------------

    !// File to be read is netcdf
    filename = 'Indata_CTM3/tropospheric_aerosols.nc'
    status = nf90_noerr  !Status is 0 and should be kept that way !!
    write(*,'(a)') '* Reading tropospheric aerosols for Fast-JX'
    write(*,'(a)') '  File: '//trim(filename)

    !// Open the existing file
    status = nf90_open(filename, nf90_nowrite, ncid)
    if (status /= nf90_noerr) call handle_error(status,'error in open')

    !// Inquire dimension ids
    status = nf90_inq_dimid(ncid,'Latitudes',lat_dim_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_dimid: Latitudes')

    status = nf90_inq_dimid(ncid,'Longitudes',lon_dim_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_dimid: Longitudes')

    status = nf90_inq_dimid(ncid,'Levels',lev_dim_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_dimid: Levels')

    status = nf90_inq_dimid(ncid,'Months',time_dim_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_dimid: Months')

    status = nf90_inq_dimid(ncid,'Num_aerosols',naer_dim_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_dimid: Num_aerosols')

    !// Inquire dimensions

    !// Latitude
    status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
    if (status /= nf90_noerr) call handle_error(status,'inquire_dimension: lat_dim_id')
    if (nlats == JPAR) then
       LINTERP = .false.
    else
       LINTERP = .true.
    end if

    !// Longitude
    status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
    if (status /= nf90_noerr) call handle_error(status,'inquire_dimension: lon_dim_id')
    if (nlons /= IPAR) then
       LINTERP = .true.
    end if

    !// Levels
    status = nf90_Inquire_Dimension(ncid,lev_dim_id,len=nlevs)
    if (status /= nf90_noerr) call handle_error(status,'inquire_dimension: lev_dim_id')

    !// Times
    status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
    if (status /= nf90_noerr) call handle_error(status,'inquire_dimension: time_dim_id')
    if (nsteps /= 12) then
       write(*,'(a)') '* Not 12 months of data in '//trim(filename)
       stop
    end if

    !// Number of aerosols
    status = nf90_Inquire_Dimension(ncid,naer_dim_id,len=naer)
    if (status /= nf90_noerr) call handle_error(status,'inquire_dimension: naer_dim_id')


    !// Allocate
    allocate( aer_in(nlons,nlats,nlevs,2), aer2d(nlons,nlats), &
         aer_xdedg(nlons+1), aer_ydedg(nlats+1),aer_ids(naer) )

    !// Get grid size
    status = nf90_inq_varid(ncid,'XDEDG',field_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_varid: XDEDG')
    status = nf90_get_var( ncid,field_id,aer_xdedg )

    status = nf90_inq_varid(ncid,'YDEDG',field_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_varid: YDEDG')
    status = nf90_get_var( ncid,field_id,aer_ydedg )

    status = nf90_inq_varid(ncid,'Aerosol_IDs',field_id)
    if (status /= nf90_noerr) call handle_error(status,'inq_varid: Aerosol_IDs')
    status = nf90_get_var( ncid,field_id,aer_ids )


    !// Check which of the wanted aerosols that are in climatology
    AERinCLIM(:) = 0
    do N = 1, NUM_CLIM
       do L = 1, NAER
          if (CLIM_ID(N) .eq. AER_IDS(L)) then
             AERinCLIM(N) = 1
             exit
          end if
       end do
    end do
    if (NUM_CLIM .lt. NAER) then
       write(*,'(a)') '  * There are more aerosols in climatology than NUM_CLIM'//&
                      ' (not a problem)'
       write(*,'(a,i3)') '    NUM_CLIM (saved in run): ',NUM_CLIM
       write(*,'(a,i3)') '    NAER (on file):          ',NAER
       write(*,'(a,20i4)') '    CLIM_ID (saved in run): ',CLIM_ID
       write(*,'(a,20i4)') '    AER_ID (on file):       ',AER_IDS
    end if


    !// Get data for this month and the next
    mon1 = MONTH
    if (MONTH .lt. 12) then
       mon2 = MONTH + 1
    else
       mon2 = 1
    end if


    !// May handle collapsed layers
    !// Number of levels in CTM; less than nlevs when collapsing layers:
    APAR = maxval(LMMAP(1:nlevs))
    !// Each full (non-collapsed) level starting points are given in LMAP
    !// LMAP(1)=1, and if layers 1:3 are collapsed, LMAP(2) = 4.
    do LL = LPAR+1,1,-1
       L   = LMMAP(LL)
       LMAP(L) = LL    !// The start levels (non-collapsed) of each CTM level
    end do

    !// Initialize
    APM(:,:,:,:,:) = 0._r8


    !// Get all defined aerosol types
    do N = 1, NUM_CLIM

       !// Aerosol not in the climatology; cycle
       if (AERinCLIM(N) .eq. 0) then
          write(*,'(a,i3,a)') '  * Aerosol ID ',CLIM_ID(N),&
               ' is not in climatology; will skip!'
          cycle
       end if

       !// Get variable ID
       write(CNUM(1:3),'(i3)') CLIM_ID(N)
       CLIM_NAME = 'Tracer'//CNUM
       status = nf90_inq_varid(ncid,trim(CLIM_NAME),field_id)
       if (status /= nf90_noerr) call handle_error(status,'inq_varid: CLIM_ID: '//CLIM_NAME)
       
       !// Get the variable month 1
       status = nf90_get_var( ncid,field_id,aer_in(:,:,:,1), &
            start=(/1, 1, 1, mon1/), count=(/nlons, nlats, nlevs, 1/) )
       if (status /= nf90_noerr) call handle_error(status,'get_var: month 1')

       !// Get the variable month 2
       status = nf90_get_var( ncid,field_id,aer_in(:,:,:,2), &
            start=(/1, 1, 1, mon2/), count=(/nlons, nlats, nlevs, 1/) )
       if (status /= nf90_noerr) call handle_error(status,'get_var: month 2')

       write(*,'(a)') '  Got variable: '//trim(CLIM_NAME)


       !// Interpolate?
       if (LINTERP) then
         !// CTM layers corresponding to NLEVS
         do L = 1, APAR

           !// Loop over all full levels for each (collapsed) CTM levels
           do LL = LMAP(L),LMAP(L+1)-1
             
             !// Interpolate month 1
             aer2d(:,:) = aer_in(:,:,L,1)
             call E_GRID(aer2d, aer_xdedg, aer_ydedg, nlons, nlats, &
                         R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
             !// Put into CTM array. Make average if collapsed layers.
             do J = 1, JPAR
               do I = 1, IPAR
                 APM(1,L,N,I,J) = APM(1,L,N,I,J) + R8XY(I,J)*XLMMAP(LL)
               end do
             end do

             !// Interpolate month 2
             aer2d(:,:) = aer_in(:,:,L,2)
             call E_GRID(aer2d, aer_xdedg, aer_ydedg, nlons, nlats, &
                         R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
             !// Put into CTM array
             do J = 1, JPAR
               do I = 1, IPAR
                 APM(2,L,N,I,J) = APM(2,L,N,I,J) + R8XY(I,J)*XLMMAP(LL)
               end do
             end do

           end do !// do LL = LMAP(L),LMAP(L+1)-1
         end do !// do L = 1, LPAR
       else
          !// This month
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  APM(1,L,N,I,J) = APM(1,L,N,I,J) + aer_in(I,J,L,1)*XLMMAP(LL)
                end do
              end do
              do L = APAR+1,lpar
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  APM(1,L,N,I,J) = 0._r8
                end do
              end do
            end do
          end do
          !// Next month
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  APM(2,L,N,I,J) = APM(2,L,N,I,J) + aer_in(I,J,L,2)*XLMMAP(LL)
                end do
              end do
              do L = APAR+1,lpar
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  APM(2,L,N,I,J) = 0._r8
                end do
              end do
            end do
          end do
       end if

       !// Set zero for value lower than ZEROPATH
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                if (APM(1,L,N,I,J) .le. ZEROPATH) APM(1,L,N,I,J) = 0._r8
                if (APM(2,L,N,I,J) .le. ZEROPATH) APM(2,L,N,I,J) = 0._r8
            end do
          end do
       end do

       !// Check min/max path values
       maxcol1=0.
       maxcol2=0.
       meancol1=0.
       meancol2=0.
       do J = 1, JPAR
          do I = 1, IPAR
           if (maxcol1 .lt. sum(APM(1,:,N,I,J))) maxcol1=sum(APM(1,:,N,I,J))
           if (maxcol2 .lt. sum(APM(2,:,N,I,J))) maxcol2=sum(APM(2,:,N,I,J))
           meancol1 = meancol1 + sum(APM(1,:,N,I,J))*AREAXY(I,J)
           meancol2 = meancol2 + sum(APM(2,:,N,I,J))*AREAXY(I,J)
          end do
       end do

       meancol1=meancol1/sum(AREAXY)
       meancol2=meancol2/sum(AREAXY)
       write(*,'(4x,a,2es12.4)') 'Max/mean column [g/m2] month 1:',&
            maxcol1,meancol1
       write(*,'(4x,a,2es12.4)') 'Max/mean column [g/m2] month 2:',&
            maxcol2,meancol2


    end do !// do N = 1, NUM_CLIM

    !// Closing file
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
  
    !// Deallocate
    deallocate(aer_in, aer2d, aer_xdedg, aer_ydedg,aer_ids)


    !// --------------------------------------------------------------------
  end subroutine get_clim_tropaerosols4fjx
  !// ----------------------------------------------------------------------



  !// --------------------------------------------------------------------
  subroutine update_tropaerosols4fjx(MONTH, LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Set tropopsheric aerosols for Fast-JX scattering.
    !// Updates monthly means every month, and interpolates linearly in time
    !// each day at 00UTC.
    !// Called from update_chemistry in main_oslo.f90.
    !//
    !// Amund Sovde, May 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MONTH
    logical, intent(in) :: LNEW_MONTH
    !// Locals
    integer :: I,J,L,N
    real(r8) :: dtfrac
    !// --------------------------------------------------------------------

    !// Skip aerosol effect on Fast-JX?
    if (.not.LJV_AEROSOL) return

    !// Read new climatology?
    !// Possible future read-in for tropospheric aerosols.
    !if (LNEW_MONTH) call get_clim_tropaerosols4fjx(MONTH)

    !// Interpolate linearly in time between two monthly means
    aerdtfrac = 1._r8 - real(jdate-1, r8) / DINM(MONTH)
    dtfrac = 1._r8 - real(jdate-1, r8) / DINM(MONTH)

   !// --------------------------------------------------------------------
  end subroutine update_tropaerosols4fjx
  !// ----------------------------------------------------------------------





  !// ----------------------------------------------------------------------
  subroutine set_aer4fjx_ctm2(MP)
    !// --------------------------------------------------------------------
    !// Set tropopsheric aerosols for Fast-JX absorbtion as in CTM2.
    !// Only to be used when NUM_AER==1.
    !//
    !// Amund Sovde, July 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, ETAA, ETAB
    use cmn_met, only: P, ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP

    !// Locals
    integer :: I,J,L,II,JJ,K          !// Indices
    integer, parameter :: LVL=51, &   !// Number of levels in standard p
                          NTROP_LEV=8 !// Number of levels used
    real(r8),dimension(LVL+1) :: BREF, PSTD
    real(r8) :: dtfrac,PJ(LPAR+2), B0, XC, PC, PB, DLOGP
    !// --------------------------------------------------------------------


    if (NUM_AER .ne. 1) then
       print*,'* aerosols2fastjx.f90: set_aer4fjx_ctm2: NUM_AER != 1'
       stop
    end if

    !// Set aerosol type
    AERTYPE(1) = 3

    !// Approx. Black Carbon up to 10 km. Surface 200 ng/m3 (Liousse et al)
    !// Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-9 g/m2/m)
    !// Below we interpolate BREF into CTM grid.
    !// Assume 10 ng/m3, with extinction coefficient 10g/m2.
    BREF(1:NTROP_LEV)       = 10._r8 * 1.0e-9_r8
    !Do I = NTROP_LEV+1,LVL+1
    BREF(NTROP_LEV+1:LVL+1) = 0._r8
    !End Do
    
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          do L=1,LPAR+1
             PJ(L)  = ETAA(L) + ETAB(L)*P(I,J)
          end do
          PJ(LPAR+2) = 0._r8

          PSTD(1) = Max(PJ(1),1000._r8)
          PSTD(2) = 1000._r8*10._r8**(-1._r8/16._r8)
          DLOGP   = 10._r8**(-1._r8/8._r8)
          do L = 3, LVL
             PSTD(L) = PSTD(L-1) * DLOGP
          end do
          PSTD(LVL+1) = 0._r8

          !// with mass (pressure) weighting, assuming constant mixing ratio and
          !// temperature half a layer on either side of the point supplied.
          do L = 1,LPAR
             B0 = 0._r8
             do K = 1,LVL
                PC = Min(PJ(L),PSTD(K))
                PB = Max(PJ(L+1),PSTD(K+1))
                if (PC.GT.PB) Then
                   XC = (PC-PB)/(PJ(L)-PJ(L+1))
                   B0 = B0 + BREF(K)*XC
                end if
             end do
             !// Currently use soot as absorbing only (type 3).
             !// Type 3 gives extinction of 5m2/g, so to get CTM2 values, we
             !// need to change B0 (which was multiplied by 10):
             !AER1P(I,J,L) = B0*2._r8*(ZOFLE(L+1,I,J) - ZOFLE(L,I,J))
             !AER1N(I,J,L) = 3
             AERPATH(L,1,I,J) = B0*2._r8*(ZOFLE(L+1,I,J) - ZOFLE(L,I,J))
          end do !// do L = 1,LPAR
          AERPATH(LPAR+1,1,I,J) = 0._r8

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine set_aer4fjx_ctm2
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
end module aerosols2fastjx
!//=========================================================================

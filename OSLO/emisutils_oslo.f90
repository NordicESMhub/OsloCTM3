!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, April 2017
!//=========================================================================
!// Emission utilities.
!//=========================================================================
module emisutils_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: emisutils_oslo
  !// DESCRIPTION: Contains utilities to read emissions into CTM.
  !//
  !// Contains:
  !//   subroutine reademis_2d
  !//   subroutine reademis_3d
  !//   subroutine reademis_stv
  !//   subroutine emisinterp2d
  !//   subroutine emis_setscaling
  !//   subroutine scale_area
  !//   subroutine get_xyedges
  !//   subroutine set_diurnal_scalings
  !//   subroutine set_vertical_scalings
  !//   subroutine emis_setscaling_2dfields
  !//   subroutine emis_diag
  !//
  !//   subroutine ctm2_rdmolec2
  !//   subroutine iiasa_read2d
  !//   subroutine gfed_read3d
  !//   subroutine volc_read1
  !//   subroutine read_retrobbh
  !//   subroutine read_bcoc_bond_2d
  !//   subroutine read_aerocom_2d
  !//   subroutine gfed4_rd
  !//   subroutine gfed4_rd_novert
  !//   subroutine gfed4_rd_novert_daily
  !//   subroutine read_ceds_cicero
  !//   subroutine ceds_biomass_burning
  !//   subroutine ceds_biomass_burning_novert
  !//   subroutine gfed4_rd_daily
  !//   subroutine getmeganmonthly
  !//
  !// Amund Sovde Haslerud, March 2017,
  !//                       November 2015, June 2013, November 2011, 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// String length of CATLIST
  integer, parameter::catlistlen=3
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='emisutils_oslo.f90'
  !// ----------------------------------------------------------------------
  public
  private scale_area, ctm2_rdmolec2, iiasa_read2d, &
       volc_read1, read_retrobbh, emisinterp2d
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine reademis_2d(INFILE,INFMT,INRES,INTYP,INCAT,INMON,INYEAR, &
       ESCENYEAR,INSCAL_LH, INSCAL_2D, INSCAL_V, &
       INVAR, IRES,JRES,MRES,MTBL,  ES_NT,ES_TRNR,ES_USCALE)
    !// --------------------------------------------------------------------
    !// Based on UCI RDEM1x1.
    !// What is done:
    !// 1. Reads 2D data into R8XY. Based on INFMT code.
    !// 2. Possibly scales with area (depending on INTYP)
    !// 3. Interpolates into emission arrays (E2DS).
    !// 4. Assigns dataset to the tracers specified to use it.
    !//    Also assign scaling factors for dataset and tracers.
    !// Arguments:
    !//    XBEDGE(1:IB+1) & YBEDGE(1:JB+1) are em grid edges in degrees
    !//    XYBOX(1:IB,1:JB) = box area (m2)
    !//    INFMT = Read-in format codes. Numbers <500 starts at (180W,90S),
    !//            while >=500 starts at (0E,90S). We try to list all
    !//            options here, but you should check code comments also.
    !//            401 = netCDF 12 month dataset, e.g. MEGAN.
    !//            402 = netCDF annual dataset, e.g. ECLIPSE.
    !//            403 = netCDF for many months, used for GAME CH4.
    !//            501 = POET emissions 2000, annual data, ascii.
    !//            502 = POET emissions 2000, 12 month data, ascii.
    !//            551 = GEIA H2S emissions, 12 month data, ascii.
    !//            552 = Older IIASA emissions such as SO2, annual data, ascii.
    !//            601 = netCDF 12 month data, e.g. RETRO.
    !//            602 = netCDF annual data, e.g. EDGARv4.
    !//            701 = BCOC 2D emissions from Tami Bond, annual, ascii.
    !//            702 = SOA AEROCOM emissions, 12 month data, ascii.
    !//    INRES = char*5 description (eg, '1x1')
    !//    INTYP = 0: leave emissions as is
    !//            1: scale by area (cm2)
    !//            2: scale by area (m2)
    !//    INCAT = Emission category name.
    !//    INMON = Which month to apply dataset (99 is for all months if
    !//            dataset is annual).
    !//    INYEAR= Which year to apply dataset (9999 is for all years).
    !//            Dataset is used when INYEAR matches MYEAR (i.e. not JYEAR).
    !//    INSCAL_LH= Diurnal (local hour) variation to apply to dataset?
    !//    INSCAL_2D= Diurnal (local hour) 2D based variation
    !//    INSCAL_V = Vertical distribution on lowermost layers?
    !//    INVAR = variable name for netCDF files.
    !//    MTBL  = Emission table number so far.
    !//    ES_NT    = Number of tracers to use the current dataset.
    !//    ES_TRNR  = The transport numbers of these tracers.
    !//    ES_USCALE= Scaling factors to apply for each tracer.
    !// --------------------------------------------------------------------
    use cmn_size, only: LSOA, NPAR
    use cmn_parameters, only: A0, CPI180
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// In
    character(len=*), intent(in) :: INFILE, INRES, INCAT, INVAR
    integer, intent(in)  :: INFMT, INTYP, INMON,INYEAR, ESCENYEAR, &
                            INSCAL_LH, INSCAL_2D, INSCAL_V, &
                            IRES,JRES,MRES, ES_NT
    integer, intent(in)  :: ES_TRNR(NPAR)
    real(r8), intent(in) :: ES_USCALE(NPAR)
    
    !// Out
    integer, intent(inout) :: MTBL

    !// Locals
    real(r8) :: XBEDGE(IRES+1), YBEDGE(JRES+1), XYBOX(JRES)
    real(r8) :: SUM1x1, SUMCTM, rpostot, rtot
    real(r8) :: R8XY(IRES,JRES,MRES)

    integer :: I, J, M, N, io_err, NSETS

    character(len=70) :: TITLE
    integer :: efnr, ii, jj, k, grid_type

    character(len=60) :: LONG,LAT,TIME
    !// Categories on file
    integer :: MIDX(MRES)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'reademis_2d'
    !// --------------------------------------------------------------------


    !// SET UP TYPICAL EMISSION GRIDS & AREAS
    !// --------------------------------------------------------------------
    if (INFMT .ge. 500) then
       !// Grid with lower-left ((1,1)-box) edges at (0E, 90S)
       grid_type = 1
    else
       !// Grid with lower-left ((1,1)-box) edges at (180W, 90S)
       grid_type = 2
    end if
    !// Other emisison grids may be needed eventually

    !// Set up grid
    call get_xyedges(IRES,JRES,XBEDGE,YBEDGE,grid_type)

    !// Grid box areas (all zonal boxes are of same size)
    do J = 1, JRES
       XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    !// --------------------------------------------------------------------

    !// File number for reading emission files
    efnr = get_free_fileid()

    !// Initialize 
    MIDX(:) = 0         !// Dataset numbers to be defined
    R8XY(:,:,:) = 0._r8 !// Emission dataset

    !// DETERMINE TYPE AND READ IN EMISSION DATA
    !// --------------------------------------------------------------------

    if (INFMT.eq.401) then
       !// netCDF file, 2D emissions
       !// Applies for e.g. MEGAN. Starts at 180W,90S.
       !// Contains 12 months of data field.
       !// May contain more datasets, but retrieves only specified dataset.

       !// Read data, return 3D fields (IRES,JRES,MRES)
       call nc_getdata_sfc(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.402) then
       !// netCDF file, 2D emissions
       !// Applies for ECLIPSE, which starts at 180W,90S

       !// Files contain annual emissions.
       !// Find these by standard netCDF commands.
       !// Retrieve one dataset.

       !// Read data, return 2D fields (IRES,JRES,1)
       call nc_getdata_sfc_2d(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.403) then
       !// netCDF file, 2D emissions, for many months, dimension 720x360
       !// Applies for GAME CH4

       !// Read data, return 3D fields (IRES,JRES,12)
       call nc_getch4emis_game(INFILE,INVAR,INYEAR,R8XY,IRES,JRES,MRES,NSETS)

       !// Scaling must be applied: Use only positive values, but scale to
       !// total of all values. Area-weighted!
       do m = 1, nsets
         rpostot = 0._r8
         rtot = 0._r8
         do J = 1, JRES
           do I = 1, IRES
             rtot = rtot + R8XY(i,j,m) * xybox(J)
             if (R8XY(i,j,m) .ge. 0._r8) &
                  rpostot = rpostot + R8XY(i,j,m) * xybox(J)
           end do
         end do
         R8XY(:,:,m) = max(R8XY(:,:,m), 0._r8) * rtot / rpostot
         write(6,'(a,f16.10)') '* Scaling dataset to positive values', &
              rtot/rpostot
      end do

    else if (INFMT.eq.404) then
       !// netCDF file, 2D emissions
       !// Applies for CTM3 MEGAN-accumulated files.
       !// Retrieve one dataset.
       !// Will overwrite XBEDGE and YBEDGE

       !// Read data, return 2D fields (IRES,JRES,1)
       call getmeganmonthly(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS, &
            XBEDGE, YBEDGE)

       !// Grid box areas (all zonal boxes are of same size)
       do J = 1, JRES
          XYBOX(J) =  A0*A0 * CPI180 * (XBEDGE(2) - XBEDGE(1)) &
               * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
       end do

    else if (INFMT.eq.460) then
       !// CEDS on CICERO format, i.e. sector datasets have separate
       !// names and data are on the format (LON,LAT,TIME).
       !// Need additional YEAR info from emission list.
       call read_ceds_cicero(INFILE, INVAR, INYEAR, R8XY, &
            IRES, JRES, MRES, NSETS, ESCENYEAR)


    else if (INFMT.eq.501) then
       !// POET 2D emissions, unit: molecule/cm2/s
       !// Contains annual data, 1 data set
       NSETS = 1

       !// Read data
       call ctm2_rdmolec2(R8XY,INFILE,efnr,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.502) then
       !// POET 2D emissions
       !// Contains monthly data, 12 data sets
       NSETS = 12

       !// Read data
       call ctm2_rdmolec2(R8XY,INFILE,efnr,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.551) then
       !// 2D emissions, unit: molecule/cm2/s
       !// Mainly used for H2S emissions.
       !// Contains 12 months of data, only one category
       NSETS = 12

       !// Open file
       open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
       if (io_err .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Problems with file: '//trim(infile)
          write(6,'(a,i12)') 'io_err: ',io_err
          stop
       end if
       do M = 1, NSETS
          !// Read all data
          read(efnr,'(I4)') N   !// month number on file
          if (N .ne. M) then
             write(6,'(a)') f90file//':'//subr// &
                  ': Problems with file: '//trim(infile)
             stop
          end if
          read(efnr,'(10E12.6)') ((R8XY(I,J,M),I=1,IRES),J=1,JRES)
       end do
       !// Close file
       close(efnr)

    else if (INFMT.eq.552) then
       !// 2D emissions, unit: kg fmm/y
       !// Used for older IIASA emissions such as SO2.
       !// Contains annual data, only one category
       NSETS = 1
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          write(6,'(a,i12)') 'INFMT: ',INFMT
          stop
       end if

       call iiasa_read2d(R8XY,INFILE,efnr,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.601) then
       !// netCDF file, 2D emissions
       !// Applies for e.g. RETRO. Starts at 0E,90S
       !// Contains 12 months of data field.
       !// May contain more datasets, but retrieves only specified dataset.

       !// Read data, return 3D fields (IRES,JRES,MRES)
       call nc_getdata_sfc(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.602) then
       !// netCDF file, 2D emissions
       !// Applies for e.g. EDGAR4.2, dimension 3600x1800

       !// Files contain annual emissions.
       !// Find these by standard netCDF commands.
       !// Retrieve one dataset.

       !// Read data, return 2D fields (IRES,JRES,1)
       call nc_getdata_sfc_2d(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.701) then
       !// BCOC 2D emissions from Tami Bond, unit on file: Gg/y
       !// Contains annual data, only one category
       !// First grid box in TBond data is located at 180W->179W & 90S->89S,
       !// but -180:0 is swapped to 181:360 during read-in. So grid_type 1
       !// is correct.
       NSETS = 1
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          write(6,'(a,i12)') 'INFMT: ',INFMT
          stop
       end if

       call read_bcoc_bond_2d(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS)

    else if (INFMT.eq.702) then
       !// Read data from AEROCOM, e.g. Secondary orcanic aerosols.
       !// 12 months, unit on file is kg/month.

       !// Test if SOA is read; Not to be used with SOA scheme.
       N = len(trim(INFILE))
       if (INFILE(N-6:N).eq.'SOA.1x1' .and. LSOA) then
          write(6,'(a)') f90file//':'//subr// &
               ': Should not use AEROCOM SOA emissions when '// &
               'running SOA scheme!'
          stop
       end if

       NSETS = 12
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          write(6,'(a,i12)') 'INFMT: ',INFMT
          stop
       end if

       call read_aerocom_2d(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS)

    else
       !// Unknown format to read!
       write(6,'(a,i5)') f90file//':'//subr// &
            ': Format not defined!: ',INFMT
       write(6,'(a)') '  INFILE: '//trim(infile)
       stop
    end if

    !// --------------------------------------------------------------------

    !// Scale with area
    call scale_area(R8XY,NSETS,XYBOX,IRES,JRES,MRES,INTYP)


    !// Interpolate and put into emission table (E2DS), assigning a number
    !// to each dataset. The numbers are stored in MIDX(:).
    !// - Routine will convert from monthly to annual dataset if
    !//   all the monthly sets are equal.
    !// - Returns NSETS=0 if all datasets are zero.
    call emisinterp2d(R8XY,IRES,JRES,MRES,NSETS,MTBL,MIDX, &
         XBEDGE, YBEDGE, INMON, INYEAR, INCAT)


    !// If NSETS > 0, the emissions are read.
    !// Assign dataset to tracers and also the scaling factors.
    if (NSETS .gt. 0) then
       call emis_setscaling(MRES,NSETS,INSCAL_LH,INSCAL_2D,INSCAL_V, &
            MTBL,MIDX,INCAT, ES_NT,ES_TRNR,ES_USCALE,2)
    end if

    !// --------------------------------------------------------------------
  end subroutine reademis_2d
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine reademis_3d(INFILE,INFMT,INRES,INTYP,INCAT,INMON,INYEAR,&
       INVAR, IRES,JRES,MRES,MTBL, ES_NT,ES_TRNR,ES_USCALE)
    !// --------------------------------------------------------------------
    !// Based on UCI RDEM_3d.
    !//
    !// Reads 3D data and interpolates into emission arrays.
    !//    XBEDGE(1:IB+1) & YBEDGE(1:JB+1) are em grid edges in degrees
    !//    XYBOX(1:IB,1:JB) = box area (m2)
    !//    INFMT = Read-in format codes. Numbers <500 starts at (180W,90S),
    !//            while >=500 starts at (0E,90S). We try to list all
    !//            options here, but you should check code comments also.
    !//            401 = netCDF 12 month dataset, e.g. MEGAN, distributes
    !//                  on 4 lowermost levels depending on category.
    !//            402 = netCDF annual dataset, e.g. ECLIPSE, distributes
    !//                  on 4 lowermost levels depending on category.
    !//            651 = Old GFED 3D emissions, 12 month data.
    !//                  *NOT IN USE*
    !//            652 = Reads IIASA 2D emissions, annual data, and
    !//                  distributed on levels 4-7.
    !//            653 = Volcanic 3D emissions (SO2), annual data, ascii.
    !//            703 = GFED BCOC with separate vertical distribution,
    !//                  12 month data, ascii.
    !//                  *OUTDATED, SHOULD BE REMOVED*
    !//    INRES = char*5 description (eg, '1x1')
    !//    INTYP = 0: leave emissions as is
    !//            1: scale by area (cm2)
    !//            2: scale by area (m2)
    !//    INCAT = Emission category name.
    !//    INMON = Which month to apply dataset (99 is for all months if
    !//            dataset is annual).
    !//    INYEAR= Which year to apply dataset (9999 is for all years).
    !//            Dataset is used when INYEAR matches MYEAR (i.e. not JYEAR).
    !//    INVAR = variable name for netCDF files.
    !//    MTBL  = Emission table number so far.
    !//    ES_NT    = Number of tracers to use the current dataset.
    !//    ES_TRNR  = The transport numbers of these tracers.
    !//    ES_USCALE= Scaling factors to apply for each tracer.
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR
    use cmn_ctm, only: XDEDG, YDEDG, LMMAP
    use cmn_chem, only: E3DSNEW, NY3TBL, NM3TBL, E3PAR, ETPAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// In
    character(len=*), intent(in) :: INFILE, INRES, INCAT, INVAR
    integer, intent(in)       :: INFMT, INTYP, INMON,INYEAR, &
                                 IRES,JRES,MRES, ES_NT
    integer, intent(in)       :: ES_TRNR(NPAR)
    real(r8), intent(in)        :: ES_USCALE(NPAR)
    
    !// Out
    integer, intent(inout) :: MTBL

    !// Locals
    real(r8) ::  XBEDGE(IRES+1), YBEDGE(JRES+1), XYBOX(JRES)
    real(r8) ::  SUM1x1, SUMCTM,ioffset,joffset,factor,rdum
    !// 3D field IRES,JRES,L1x1,MRES (S is sets/months)
    integer, parameter :: L1x1 = 6
    real(r8) ::  R8XYZ_IN(IRES,JRES,L1x1,MRES),R8XY(IRES,JRES,MRES)
    !// 3D field IRES,JRES,LPAR,MRES
    real(r8) ::  R8XYZ(IRES,JRES,LPAR,MRES)
    !// 2D temporary field for interpolation (NO moments for 3D)
    real(r8) ::  R8CTM(IPAR,JPAR)

    integer ::  I, J, L, M, N, NN, ML, io_err, NLEVS,NSETS

    character(len=70) :: TITLE
    integer :: efnr,ii,jj,k, grid_type
    real(r8) :: rsum, vdist(4)
    integer :: MIDX(MRES)

    character(len=60) :: LONG,LAT,TIME
    real(r8),parameter :: Z3=1._r8/3._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'reademis_3d'
    !// --------------------------------------------------------------------


    !// SET UP TYPICAL EMISSION GRIDS & AREAS
    !// --------------------------------------------------------------------
    if (INFMT .ge. 500) then
       !// Grid with lower-left ((1,1)-box) edges at (0E, 90S)
       grid_type = 1
    else
       !// Grid with lower-left ((1,1)-box) edges at (180W, 90S)
       grid_type = 2
    end if
    !// Other emisison grids may be needed here

    !// Set up grid
    call get_xyedges(IRES,JRES,XBEDGE,YBEDGE,grid_type)

    !// Grid box areas
    do J = 1, JRES
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do
    !// --------------------------------------------------------------------

    !// File number for reading emission files
    efnr = get_free_fileid()

    !// Initialize 
    MIDX(:) = 0


    !// Possible category vertical scalings for the 4 lowest model layers.
    !// These are only used for some read-in formats (INFMT), to
    !// distribute chimney-gases vertically according to category.
    !// E.g. VOCs are perhaps not emitted through combustion
    !// processes, and should apply a dfferent scaling.
    !// Levels, globally averaged: 1: 0-16m, 2: 16-41m, 3: 41-77m, 4: 77-128m
    if (INCAT .eq. 'POW') then
       vdist(:) = (/0.25_r8,0.125_r8,0.125_r8,0.5_r8/)
    else if (INCAT .eq. 'INC') then
       vdist(:) = (/0.25_r8,0.125_r8,0.125_r8,0.5_r8/)
    else if (INCAT .eq. 'RES') then
       vdist(:) = (/0.6_r8,0.35_r8,0.05_r8,0._r8/)
    else if (INCAT .eq. 'SHI') then
       vdist(:) = (/0.5_r8,0.25_r8,0.25_r8,0._r8/)
    else if (INCAT .eq. 'AWB') then
       vdist(:) = (/0.6_r8,0.3_r8,0._r8,0.1_r8/)
    else if (INCAT .eq. 'WAS') then
       vdist(:) = (/0.9_r8,0.1_r8,0._r8,0._r8/)
    else
       vdist(1)   = 1._r8
       vdist(2:4) = 0._r8
    end if


    !// DETERMINE TYPE AND READ IN EMISSION DATA
    !// --------------------------------------------------------------------
    if (INFMT.eq.401) then
       !// Read 2D-format code 401 and distribute vertically in the lowest
       !// 4 layers, ala CTM2. File is netCDF file, 2D emissions, conatining
       !// 12 months of data.
       !// Applies for e.g. ECLIPSE monthly data. Starts at 90S,180W.

       !// Read data, return 3D fields (IRES,JRES,MRES)
       call nc_getdata_sfc(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

       R8XYZ(:,:,:,1:NSETS) = 0._r8
       do M = 1, NSETS
          R8XYZ(:,:,LMMAP( 1),M) = R8XYZ(:,:,LMMAP( 1),M) + vdist(1)*R8XY(:,:,M)
          R8XYZ(:,:,LMMAP( 2),M) = R8XYZ(:,:,LMMAP( 2),M) + vdist(2)*R8XY(:,:,M)
          R8XYZ(:,:,LMMAP( 3),M) = R8XYZ(:,:,LMMAP( 3),M) + vdist(3)*R8XY(:,:,M)
          R8XYZ(:,:,LMMAP( 4),M) = R8XYZ(:,:,LMMAP( 4),M) + vdist(4)*R8XY(:,:,M)
       end do

       !// Interpolate and put into arrays
       do M=1,NSETS

          !// Update 3D table index
          MTBL = MTBL + 1
          if (MTBL .gt. E3PAR) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': MTBL > E3PAR: ', MTBL, E3PAR
             write(6,'(a,i5)') '    INFMT:',INFMT
             stop
          end if
          !// Save table index for this dataset
          MIDX(M) = MTBL

          !// Loop over each layer
          do L=1,LMMAP(4)
             !// Interpolate into R8CTM (no moments, only mean field
             call E_GRID(R8XYZ(:,:,L,M),XBEDGE,YBEDGE,IRES,JRES, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)
             !// Update emission table
             E3DSNEW(L,:,:,MTBL) = E3DSNEW(L,:,:,MTBL) + R8CTM(:,:)
          end do

          !// Table applies for all months in year
          if (NSETS.ne.12) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': WRONG INMON and NSETS',INMON,NSETS
             stop
          end if
          NM3TBL(MTBL) = M
          NY3TBL(MTBL) = INYEAR

       end do !// do M=1,NSETS

    else if (INFMT.eq.402) then
       !// Read 2D-format code 402 and distribute vertically in the lowest
       !// 4 layers, ala CTM2. File is netCDF file, 2D emissions, conatining
       !// annual data.
       !// Applies for e.g. ECLIPSE monthly data. Starts at 90S,180W.

       !// Read data, return 3D fields (IRES,JRES)
       call nc_getdata_sfc_2d(INFILE,INVAR,R8XY(:,:,1),IRES,JRES,MRES,NSETS)

       R8XYZ(:,:,:,:) = 0._r8
       R8XYZ(:,:,LMMAP( 1),1) = R8XYZ(:,:,LMMAP( 1),1) + vdist(1)*R8XY(:,:,1)
       R8XYZ(:,:,LMMAP( 2),1) = R8XYZ(:,:,LMMAP( 2),1) + vdist(2)*R8XY(:,:,1)
       R8XYZ(:,:,LMMAP( 3),1) = R8XYZ(:,:,LMMAP( 3),1) + vdist(3)*R8XY(:,:,1)
       R8XYZ(:,:,LMMAP( 4),1) = R8XYZ(:,:,LMMAP( 4),1) + vdist(4)*R8XY(:,:,1)

       !// Update 3D table index
       MTBL = MTBL + 1
       if (MTBL .gt. E3PAR) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': MTBL > E3PAR: ', MTBL, E3PAR
          write(6,'(a,i5)') '    INFMT: ',INFMT
          stop
       end if

       !// Loop over each layer
       do L = 1, LMMAP(4)
          !// Interpolate into R8CTM (no moments, only mean field
          call E_GRID(R8XYZ(:,:,L,1),XBEDGE,YBEDGE,IRES,JRES, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)
          !// Update emission table
          E3DSNEW(L,:,:,MTBL) = E3DSNEW(L,:,:,MTBL) + R8CTM(:,:)
       end do

       !// Table applies for all months in year
       if (NSETS .ne. 1) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': WRONG INMON and NSETS',INMON,NSETS
          stop
       end if
       NM3TBL(MTBL) = INMON
       NY3TBL(MTBL) = INYEAR

    else if (INFMT.eq.651) then

       !// GFED 3D emissions, unit: kg/month
       !// Contains monthly data, only one category.
       !// Data are given in 3D and put into the corresponding CTM levels
       NSETS = 12
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          stop
       end if

       !// Read 2D data
       call gfed_read3d(R8XYZ_IN,INFILE,efnr,IRES,JRES,L1x1,MRES,NSETS,NLEVS)

       !// Move into CTM vertical grid
       !// Biomass burning is distributed in 6 height intervals
       !// Layer A    0 -  100 m
       !// Layer B  100 -  500 m
       !// Layer C  500 - 1000 m
       !// Layer D 1000 - 2000 m
       !// Layer E 2000 - 3000 m
       !// Layer F 3000 - 6000 m
       R8XYZ(:,:,:,1:NSETS) = 0._r8
       do M = 1, NSETS
          R8XYZ(:,:,LMMAP( 1),M) = R8XYZ(:,:,LMMAP( 1),M) + 0.25_r8*R8XYZ_IN(:,:,1,M)
          R8XYZ(:,:,LMMAP( 2),M) = R8XYZ(:,:,LMMAP( 2),M) + 0.25_r8*R8XYZ_IN(:,:,1,M)
          R8XYZ(:,:,LMMAP( 3),M) = R8XYZ(:,:,LMMAP( 3),M) + 0.25_r8*R8XYZ_IN(:,:,1,M)
          R8XYZ(:,:,LMMAP( 4),M) = R8XYZ(:,:,LMMAP( 4),M) + 0.25_r8*R8XYZ_IN(:,:,1,M)
          R8XYZ(:,:,LMMAP( 5),M) = R8XYZ(:,:,LMMAP( 5),M) + 0.25_r8*R8XYZ_IN(:,:,2,M)
          R8XYZ(:,:,LMMAP( 6),M) = R8XYZ(:,:,LMMAP( 6),M) + 0.25_r8*R8XYZ_IN(:,:,2,M)
          R8XYZ(:,:,LMMAP( 7),M) = R8XYZ(:,:,LMMAP( 7),M) + 0.25_r8*R8XYZ_IN(:,:,2,M)
          R8XYZ(:,:,LMMAP( 8),M) = R8XYZ(:,:,LMMAP( 8),M) + 0.25_r8*R8XYZ_IN(:,:,2,M)
          R8XYZ(:,:,LMMAP( 9),M) = R8XYZ(:,:,LMMAP( 9),M) +   Z3 * R8XYZ_IN(:,:,3,M)
          R8XYZ(:,:,LMMAP(10),M) = R8XYZ(:,:,LMMAP(10),M) +   Z3 * R8XYZ_IN(:,:,3,M)
          R8XYZ(:,:,LMMAP(11),M) = R8XYZ(:,:,LMMAP(11),M) +   Z3 * R8XYZ_IN(:,:,3,M)
          R8XYZ(:,:,LMMAP(12),M) = R8XYZ(:,:,LMMAP(12),M) + 0.25_r8*R8XYZ_IN(:,:,4,M)
          R8XYZ(:,:,LMMAP(13),M) = R8XYZ(:,:,LMMAP(13),M) + 0.25_r8*R8XYZ_IN(:,:,4,M)
          R8XYZ(:,:,LMMAP(14),M) = R8XYZ(:,:,LMMAP(14),M) + 0.25_r8*R8XYZ_IN(:,:,4,M)
          R8XYZ(:,:,LMMAP(15),M) = R8XYZ(:,:,LMMAP(15),M) + 0.25_r8*R8XYZ_IN(:,:,4,M)
          R8XYZ(:,:,LMMAP(16),M) = R8XYZ(:,:,LMMAP(16),M) +   Z3 * R8XYZ_IN(:,:,5,M)
          R8XYZ(:,:,LMMAP(17),M) = R8XYZ(:,:,LMMAP(17),M) +   Z3 * R8XYZ_IN(:,:,5,M)
          R8XYZ(:,:,LMMAP(18),M) = R8XYZ(:,:,LMMAP(18),M) +   Z3 * R8XYZ_IN(:,:,5,M)
          R8XYZ(:,:,LMMAP(19),M) = R8XYZ(:,:,LMMAP(19),M) + 0.20_r8*R8XYZ_IN(:,:,6,M)
          R8XYZ(:,:,LMMAP(20),M) = R8XYZ(:,:,LMMAP(20),M) + 0.20_r8*R8XYZ_IN(:,:,6,M)
          R8XYZ(:,:,LMMAP(21),M) = R8XYZ(:,:,LMMAP(21),M) + 0.20_r8*R8XYZ_IN(:,:,6,M)
          R8XYZ(:,:,LMMAP(22),M) = R8XYZ(:,:,LMMAP(22),M) + 0.20_r8*R8XYZ_IN(:,:,6,M)
          R8XYZ(:,:,LMMAP(23),M) = R8XYZ(:,:,LMMAP(23),M) + 0.20_r8*R8XYZ_IN(:,:,6,M)
       end do


       !// Interpolate and put into arrays
       do M = 1, NSETS

          !// Update 3D table index
          MTBL = MTBL + 1
          if (MTBL .gt. E3PAR) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': MTBL > E3PAR: ', MTBL, E3PAR
             write(6,'(a,i5)') '    INFMT: ',INFMT
             stop
          end if
          !// Save table index for this dataset
          MIDX(M) = MTBL

          !// Loop over each layer
          do L = 1, LMMAP(23)
             !// Interpolate into R8CTM (no moments, only mean field
             call E_GRID(R8XYZ(:,:,L,M),XBEDGE,YBEDGE,IRES,JRES, &
                  R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)
             !// Update emission table
             E3DSNEW(L,:,:,MTBL) = E3DSNEW(L,:,:,MTBL) + R8CTM(:,:)
          end do


          !// Table applies for all months in year
          if (NSETS.ne.12) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': WRONG INMON and NSETS',INMON,NSETS
             stop
          end if
          NM3TBL(MTBL) = M
          NY3TBL(MTBL) = INYEAR

       end do !// do M=1,NSETS

    else if (INFMT.eq.652) then
       !// 3D emissions, unit: kg/y
       !// Contains annual data, only one category.
       !// Data are given in 2D and put into the lowermost levels of the model.
       NSETS = 1
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          stop
       end if

       !// Read 2D data (surface, distribute verticallt below)
       call iiasa_read2d(R8XYZ,INFILE,efnr,IRES,JRES,MRES,NSETS)

       !// Interpolate into R8CTM (no moments, only mean field
       call E_GRID(R8XYZ(:,:,1,1),XBEDGE,YBEDGE,IRES,JRES, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,1)

       !// Update 3D table index
       MTBL = MTBL + 1
       if (MTBL .gt. E3PAR) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': MTBL > E3PAR: ', MTBL, E3PAR
          write(6,'(a,i5)') '    INFMT: ',INFMT
          stop
       end if
       !// Save table index for this dataset
       MIDX(1) = MTBL

       !// Put into 3D
       if (LPARW.eq.40.or.LPARW.eq.60) then
          !// 15% in layer 4  80-130m;  35% in layer 5 130-200m
          !// 40% in layer 6 200-280m;  10% in layer 7 280-390m
          E3DSNEW(LMMAP(4),:,:,MTBL) = E3DSNEW(LMMAP(4),:,:,MTBL) + 0.15_r8*R8CTM(:,:)
          E3DSNEW(LMMAP(5),:,:,MTBL) = E3DSNEW(LMMAP(5),:,:,MTBL) + 0.35_r8*R8CTM(:,:)
          E3DSNEW(LMMAP(6),:,:,MTBL) = E3DSNEW(LMMAP(6),:,:,MTBL) + 0.40_r8*R8CTM(:,:)
          E3DSNEW(LMMAP(7),:,:,MTBL) = E3DSNEW(LMMAP(7),:,:,MTBL) + 0.10_r8*R8CTM(:,:)
       else
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': Not coded for this verticalresolution yet',IPAR,JPAR
          stop
       end if

       !// Table applies for all months in year
       if (INMON .ne. 99) then
          write(6,'(a,i5,a,i5)') f90file//':'//subr// &
               ': WRONG INMON ',inmon,' for 3D code ',infmt
          stop
       end if
       NM3TBL(MTBL) = INMON
       NY3TBL(MTBL) = INYEAR

    else if (INFMT.eq.653) then

       !// Volcanic 3D emissions, unit: kg/y
       !// Contains annual data, only one category.
       !// Data are given in 3D and put into the corresponding CTM levels
       NSETS = 1
       if (IRES.ne.360) then
          write(6,'(a)') f90file//':'//subr// &
               ': Wrong resolutuion: '//trim(infile)
          stop
       end if

       !// Read volcano data, interpolate ant put into emission tables
       call volc_read1(INFILE,efnr,IRES,JRES,MRES,NSETS,MTBL,MIDX, &
            XBEDGE, YBEDGE, XDEDG, YDEDG)

       !// Table applies for all months in year
       if (INMON.ne.99) then
          write(6,'(a,i5,a,i5)') f90file//':'//subr// &
               ': WRONG INMON ',inmon,' for 3D code ',infmt
          stop
       end if
       NM3TBL(MTBL) = INMON
       NY3TBL(MTBL) = INYEAR


    else
       !// Unknown format to read!
       write(6,'(a,i5)') f90file//':'//subr// &
            ': Format not defined!: ',INFMT
       stop
    end if


    !// Emissions are ok, now to the scaling factors
    !// In the future, NSETS may be zero if dataset is zero
    !// LSCAL_LH not set up yet: Use 0.
    !// LSCAL_2D not allowed: Use 0.
    !// LSCAL_V is of course not allowed. Use 0.
    call emis_setscaling(MRES,NSETS,0,0,0,MTBL,MIDX,INCAT, &
         ES_NT,ES_TRNR,ES_USCALE,3)

    !// --------------------------------------------------------------------
  end subroutine reademis_3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine reademis_stv(INFILE,INFMT,INFACT,INRES,INYEAR,INTYP,INCAT, &
       IRES,JRES,MRES,INVAR)
    !// --------------------------------------------------------------------
    !// Reads files for short term variations. Typically 2D fields which are
    !// modified later based on meteorology. Routine based on OC_RDEMxxx.
    !//
    !// Reads 2D data into R8XY and interpolates into specific arrays.
    !//    XBEDGE(1:IRES+1) & YBEDGE(1:JRES+1) are em grid edges in degrees
    !//    XYBOX(1:IRES,1:JRES) = box area (m2)
    !//Currently:
    !//    INFMT > 5??  : Fields with grid type 1
    !//    INRES = char*5 description (eg, '1x1')
    !//    INTYP = 0, leave emissions as is, = 1 scale by area (cm2),
    !//            = 2 scale by area (m2)
    !//    INVAR = Input field
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG, PLAND
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    use dust_oslo, only: dust_set_mbl_name_ff
    use sulphur_oslo, only: DMSseaconc
    use emissions_megan, only: megan_input
    use emissions_volcanoes, only: init_volcPATH
    use emissions_ocean, only: emissions_ocean_organiccarbon_setpath
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: INFILE, INRES, INCAT, INVAR
    integer, intent(in)       :: INFMT, INYEAR, INTYP, IRES,JRES,MRES
    real(r8), intent(in)      :: INFACT

    !// Locals
    real(r8) ::  XBEDGE(IRES+1), YBEDGE(JRES+1), XYBOX(JRES)
    real(r8) ::  SUM1x1, SUMCTM
    real(r8) ::  R8XY(IRES,JRES,MRES)
    integer ::  I, J, L, M, N, NN, ML, io_err, NLEVS
    character(len=70) :: TITLE
    integer :: efnr,ii,jj,k, grid_type, NSETS
    character(len=60) :: LONG,LAT,TIME
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'reademis_stv'
    !// --------------------------------------------------------------------


    !// SET UP TYPICAL SURFACE GRIDS & AREAS
    !// --------------------------------------------------------------------
    if (INFMT .ge. 500) then
       !// Grid with lower-left ((1,1)-box) edges at (0E, 90S)
       grid_type = 1
    else
       !// Grid with lower-left ((1,1)-box) edges at (180W, 90S)
       grid_type = 2
    end if
    !// Other emisison grids may be needed here

    !// Set up grid
    call get_xyedges(IRES,JRES,XBEDGE,YBEDGE,grid_type)

    !// Grid box areas
    do J = 1, JRES
       XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
            * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
    end do

    !// --------------------------------------------------------------------

    !// File number for reading emission files
    efnr = get_free_fileid()


    if (INFMT.eq.551) then
       !// 2D field DMSseaconc [nM m2/L] Kettle et al.
       if (INVAR .ne. 'DMSseaconc') then
          write(6,'(a)') f90file//':'//subr// &
               ': Field not DMSseaconc: '//trim(INVAR)
          stop
       end if
       if (INTYP .ne. 2) then
          write(6,'(a)') f90file//':'//subr// &
               ': DMSseaconc INTYP should be 2!'
          stop
       end if
       !// Initialize
       NSETS = 12
       R8XY(:,:,:) = 0._r8


       !// Open file
       open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
       if (io_err .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Problems with file '//trim(infile)
          write(6,'(a,i4)') '    io_err:',io_err
          stop
       end if
       do M = 1, NSETS
          !// Read all data
          read(efnr,'(I4)') N   !// month number on file
          if (N .ne. M) then
             write(6,'(a)') f90file//':'//subr// &
                  ': Problems with file '//trim(infile)
             write(6,'(a,2(1x,i4))') '    N .ne. M:',N,M
             stop
          end if
          read(efnr,'(10E12.6)') ((R8XY(I,J,M),I=1,IRES),J=1,JRES)
       end do
       !// Close file
       close(efnr)

       !// Scale with area: yes, input is per volume
       call scale_area(R8XY,NSETS,XYBOX,IRES,JRES,MRES,INTYP)

       !// Change from nmol m2/L to kg/m: *1.d-9*Mw(DMS)*1.d3/(1.d3L/m3)
       R8XY(:,:,:) = R8XY(:,:,:) * 62.e-9_r8
       !// This value will be multiplied with a velocity to create a flux.
       !// Taken care of in the source treatment each operator time step.
 


       !// Interpolate dataset into R8CTM
!$omp parallel private (M) &
!$omp          shared (NSETS,R8XY,IRES,JRES,XBEDGE,YBEDGE) &
!$omp          shared (DMSseaconc,XDEDG,YDEDG) &
!$omp          default(NONE)
!$omp do
       do M = 1, NSETS
          call E_GRID(R8XY(:,:,M),XBEDGE,YBEDGE,IRES,JRES, &
               DMSseaconc(:,:,M),XDEDG,YDEDG,IPAR,JPAR,1)
          write(6,'(a,1x,i2,2(1x,es20.13))') '    * Interpolated month', &
               M,sum(R8XY(:,:,M)),sum(DMSseaconc(:,:,M))
       end do
!$omp end do
!$omp end parallel

       !// No need to scale with ocean fraction; assume 1x1 is good enough to
       !// resolve coasts and that the mass put into current resolution is
       !// correct.
       !// Final unit is kg/m, to be multiplied with a velocity m/s to get the
       !// source flux kg/s.

    else if (INFMT.eq.552) then
       !// 2D field DMSseaconc for 12 months (Lana etal 2011, GBC)
       if (INVAR .ne. 'DMSseaconc') then
          write(6,'(a)') f90file//':'//subr// &
               ': Field not DMSseaconc: '//trim(INVAR)
          stop
       end if
       if (INTYP .ne. 2) then
          write(6,'(a)') f90file//':'//subr// &
               ': DMSseaconc INTYP should be 2!'
          stop
       end if
       !// Initialize
       NSETS = 12
       R8XY(:,:,:) = 0._r8

       !// Read data [nM], return 3D fields (IRES,JRES,MRES)
       call nc_getdata_sfc(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)

       !// Scale with area: yes, input is per volume
       call scale_area(R8XY,NSETS,XYBOX,IRES,JRES,MRES,INTYP)

       !// Change from nmol/L*m2 to kg/m: *1.d-9*Mw(DMS)*1.d3/(1.d3L/m3)
       R8XY(:,:,:) = R8XY(:,:,:) * 62.e-9_r8
       !// This value will be multiplied with a velocity to create a flux.
       !// Taken care of in the source treatment each operator time step.
 
       !// To calculate a flux [kg/s], the final unit must be kg/m. Before
       !// chemistry, we multiply with a velocity [m/s] we get the source
       !// flux kg/s.
       !// Hence, we do not divide by area after interpolation.

       !// Interpolate dataset into R8CTM
!$omp parallel private (M) &
!$omp          shared (NSETS,R8XY,IRES,JRES,XBEDGE,YBEDGE) &
!$omp          shared (DMSseaconc,XDEDG,YDEDG) &
!$omp          default(NONE)
!$omp do
       do M = 1, NSETS
          call E_GRID(R8XY(:,:,M),XBEDGE,YBEDGE,IRES,JRES, &
               DMSseaconc(:,:,M),XDEDG,YDEDG,IPAR,JPAR,1)
          write(6,'(a,1x,i2,2(1x,es20.13))') '    * Interpolated month', &
               M, sum(R8XY(:,:,M)), sum(DMSseaconc(:,:,M))
       end do
!$omp end do
!$omp end parallel

       !// No need to scale with ocean fraction; assume 1x1 is good enough to
       !// resolve coasts and that the mass put into current resolution is
       !// correct.

    else if (INFMT.eq.553) then
       !// Read path for volcanic SO2 emissions. These emissions are
       !// from HTAP, and are given as daily emissions. Therefore we
       !// only set the path and year here, while reading the data later.
       call init_volcPATH(INFILE,INYEAR,INFMT)

    else if (INFMT.eq.554) then
       !// Read path for volcanic SO2 emissions. These emissions are
       !// from AEROCOM (kg/yr), which will be converted to kg/s.
       call init_volcPATH(INFILE,INYEAR,INFMT)

    else if (INFMT.eq.555) then
       !// Include organic carbon emissions from ocean.
       !// Sets path to chlorophyll A data, and emission flag if species
       !// is included.
       call emissions_ocean_organiccarbon_setpath(INFILE,INYEAR)
       !// Will do further initialisation in input_oslo.f90

    else if (INFMT.eq.556) then
       !// Set name of mobilisation map and fudge factor to read
       !// from dust input.
       call dust_set_mbl_name_ff(INVAR, INFACT)

    else if (INFMT.eq.557) then
       !// MEGAN emissions: set true and send in path
       call megan_input(.true., INFILE)

    else
       !// Unknown format to read!
       write(6,'(a,i5)') f90file//':'//subr// &
            ': Format not defined!: ',INFMT
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine reademis_stv
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emisinterp2d(R8XY,IRES,JRES,MRES,NSETS,MTBL,MIDX, &
   XBEDGE, YBEDGE,INMON,INYEAR,CATNAME)
    !// --------------------------------------------------------------------
    !// Interpolates field R8XY(IRES,JRES,MRES) for 
    !// sets 1:NSETS (i.e. R8XY(IRES,JRES,1:NSETS)), using its
    !// boundary information XBEDGE and YBEDGE.
    !// Uses UCI E_GRID, and works for any IRES, at least >= 360.
    !//
    !// This routine is parallelized over NSETS.
    !//
    !// Ole Amund Sovde, March 2012, August 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR
    use cmn_ctm, only: XDEDG, YDEDG
    use cmn_chem, only: NE2DS, E2DS, NY2TBL, NM2TBL, ETPAR
    use cmn_oslo, only: ECATNAMES
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IRES,JRES,MRES,INMON,INYEAR
    real(r8), intent(in) :: R8XY(IRES,JRES,MRES)
    real(r8), intent(in) :: XBEDGE(IRES+1), YBEDGE(JRES+1)
    character(len=*), intent(in) :: CATNAME
    !// Input/Output
    integer, intent(inout) :: MTBL, NSETS
    !// Output
    integer, intent(out) :: MIDX(MRES)

    !// Locals
    real(r8)  :: R8CTM(IPAR,JPAR,NE2DS)
    real(r8)  :: SUM1x1, SUMCTM
    integer :: C1, S1, NN, I,J, NUMSETS,MPRIV
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emisinterp2d'
    !// --------------------------------------------------------------------

    write(6,'(a,i3)') subr//': Interpolating # of datasets: ',NSETS

    !// Check the sets; maybe all sets are equal: use one field.
    !// If set is zero, skip it.
    numsets = NSETS
    sumctm = 0._r8
    if (numsets .ne. 1) then
       do S1 = 2, NSETS
          !// There is very little chance that the NSETS fields
          !// are different but end up with zero difference
          sumctm = sumctm + sum(R8XY(:,:,S1-1) - R8XY(:,:,S1))
          if (sumctm .ne. 0._r8) exit
       end do
       if (sumctm .eq. 0._r8 ) then
          if (sum(R8XY(:,:,1)) .le. 0._r8) then
             !// All fields are zero
             write(6,'(a)') 'ZERO; All datasets are zero: Skipping'
             numsets = 0
          else
             !// All sets for this category are equal
             write(6,'(a)') 'CHANGE to ANNUAL! '
             numsets = 1
          end if
       end if
    end if

    !// Do interpolation. It is parallellized over datasets.
    !// --------------------------------------------------------------------
!$omp parallel private (S1,MPRIV,NN,R8CTM,SUM1x1,SUMCTM) &
!$omp          shared (R8XY,MTBL,numsets,IRES,JRES,XBEDGE,YBEDGE) &
!$omp          shared (INMON,INYEAR,NM2TBL,NY2TBL,XDEDG,YDEDG,MIDX) &
!$omp          shared (E2DS) &
!$omp          default(NONE)
!$omp do
    do S1 = 1, numsets

       !// Interpolate dataset into R8CTM
       if (IPAR .ne. IRES .or. JPAR .ne. JRES) then
          call E_GRID(R8XY(:,:,S1),XBEDGE,YBEDGE,IRES,JRES, &
               R8CTM,XDEDG,YDEDG,IPAR,JPAR,NE2DS)
       else
          !// No need to interpolate
          do J = 1, JPAR
             do I = 1, IPAR
                R8CTM(I,J,1) = R8XY(I,J,S1)
             end do
          end do
       end if

       !// Check sums before/after
       SUM1x1 = sum(R8XY(:,:,S1))
       SUMCTM = sum(R8CTM(:,:,1))

       !// Check totals before/after
       if (abs(sum1x1-sumctm)/sum1x1 .gt. 1.e-3_r8) then
          write(6,'(a)') f90file//':'//subr// &
               ': Mismatch TOTALS before/after interpolation'
          write(6,'(a,2e20.12)') '  * TOTALS before/after: ',SUM1x1,SUMCTM
          write(6,'(a,2e20.12)') '    DIFF (%)',(sum1x1-sumctm)/sum1x1*100._r8
          stop
       end if


       !// Emission dataset is now interpolated.
       !// Assign new table index. MTBL is not updated until after the
       !// parallel S1-loop, so we call the dataset number is MPRIV:
       MPRIV = MTBL + S1
       if (MPRIV .gt. ETPAR) then
          write(6,'(a,2i5)') f90file//':'//subr//': MTBL > ETPAR',MPRIV,ETPAR
          stop
       end if
       !// Save table index for this dataset
       MIDX(S1) = MPRIV


       !// Put dataset into 2D-emis table (entries NE2DS>1 are moments)
       do NN = 1, NE2DS
          E2DS(:,:,NN,MPRIV) = R8CTM(:,:,NN)
       end do
  
       !// Assign which month table applies for
       if (numsets .eq. 1) then
          NM2TBL(MPRIV) = INMON  !// Month defined by emission list
       else
          NM2TBL(MPRIV) = S1     !// Month defined by dataset number
       end if
       !// Assign which year table applies for
       NY2TBL(MPRIV) = INYEAR

    end do !// do S1 = 1, numsets
!$omp end do
!$omp end parallel

    !// Increase emission table index
    MTBL = MTBL + numsets
    !// Return number of assigned datasets
    NSETS = numsets

    !// --------------------------------------------------------------------
  end subroutine emisinterp2d
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine emis_setscaling(MRES,NSETS, INSCAL_LH, INSCAL_2D, INSCAL_V, &
       MTBL,MIDX, CATNAME, ES_NT,ES_TRNR,ES_USCALE,DIM)
    !// --------------------------------------------------------------------
    !// Assigns datasets to tracers, and also sets scalings for tracers
    !// using the emission datasets.
    !// Some other table properties are also set.
    !//
    !// This routine is called from both reademis_2d and reademis_3d, and
    !// does slightly different things for 2d and 3d. Here is what is
    !// assigned:
    !//
    !// 2D:
    !//   E2LTBL: Table in use for tracer/dataset
    !//   E2STBL: Scaling factor for tracer/dataset
    !//   E2LocHourTBL: Local hour variation number for tracer/dataset
    !//   E22dTBL: Diurnal 2d variation number for tracer/dataset
    !//   E2vertTBL: Vertical distribution into several layers
    !//   E2CTBL: Category number for dataset
    !// 3D:
    !//   E3LTBL: Table in use for tracer/dataset
    !//   E3STBL: Scaling factor for tracer/dataset
    !//
    !// Ole Amund Sovde, November 2015, March 2012, August 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_chem, only: TNAME, E2DS, E2LTBL, E2STBL, &
         E3DSNEW, E3LTBL, E3STBL
    use cmn_oslo, only: NECAT, ECATNAMES, E2CTBL, E2LocHourTBL, E22dTBL, &
         NE2vertVARS, E2vertTBL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MRES, ES_NT, NSETS, &
         INSCAL_LH, INSCAL_2D, INSCAL_V, MTBL, MIDX(MRES)
    integer, intent(in) :: ES_TRNR(NPAR)
    character(len=*), intent(in) :: CATNAME
    real(r8), intent(in) :: ES_USCALE(NPAR)
    integer, intent(in)  :: DIM

    !// Locals
    integer :: I,J,N,M, NC, TN, CATidx, S1
    real(r8) :: ESCALE, ADD2EDIAG
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_setscaling'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr// &
         ': scaling properties for tracers vs datasets'

    !// 1. First assign properties only to dataset
    !// 2. Then assign properties of tracers vs dataset

    !// Diagnose categories; sum up per category (for both 2d and 3d).
    !// Find category index
    do CATidx = 1, NECAT
       if (CATNAME .eq. ECATNAMES(CATidx)) exit
    end do
    if (CATidx .gt. NECAT) then
       !// Wrong species in emission list
       write(6,'(a)') f90file//':'//subr// &
            ': UNKNOWN emission category: '//CATNAME
       write(6,'(a)') '  Must update global parameter ECATNAMES'
       stop
    end if


    !// 1. Assign properties only to dataset
    !// --------------------------------------------------------------------
    do S1 = 1, NSETS

       !// Index of the emission table for this set/month (S1)
       M = MIDX(S1)

       if (M .le. 0) cycle !// Failsafe: Go to next S1 if M == 0

       !// Check M (in case it has not been done before)
       if (M .gt. MTBL) then
          write(6,'(a,2i5)') f90file//':'//subr// &
               ': MTBL too small? ', M, MTBL
          stop
       end if

       !// 2D or 3D field
       if (DIM .eq. 2) then

          !// Set category number for dataset.
          if (E2CTBL(M) .gt. 0) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': Seriously wrong: category already set!',S1,N
             write(6,'(a,2i5)') 'E2CTBL(M), CATidx:', E2CTBL(M), CATidx
             stop
          end if
          E2CTBL(M) = CATidx

          !// Set flag for diurnal (local hour) variation on dataset
          !// (i.e. all tracers will use it)
          if (E2LocHourTBL(M) .gt. 0) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': Seriously wrong: E2LocHourTBL already set!',S1,N
             write(6,'(a,2i5)') 'E2LocHourTBL(M), INSCAL_LH:',&
                  E2LocHourTBL(M), INSCAL_LH
             write(6,'(a)') 'This should NEVER happen!'
             stop
          end if
          E2LocHourTBL(M) = INSCAL_LH

          !// Set 2D-based scaling flags
          if (E22dTBL(M) .gt. 0) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': Seriously wrong: E22dTBL already set!',S1,N
             write(6,'(a,2i5)') 'E22dTBL(M), INSCAL_2D:',&
                  E22DTBL(M), INSCAL_2D
             write(6,'(a)') 'This should NEVER happen!'
             stop
          end if
          E22dTBL(M) = INSCAL_2D

          !// Set E2vertTBL
          if (INSCAL_V .gt. NE2vertVARS) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': INSCAL_V > NE2vertVARS: ',INSCAL_V,NE2vertVARS
             stop
          end if
          E2vertTBL(M) = INSCAL_V

       else if (DIM .eq. 3) then

          !// Set flag for diurnal (local hour) variation on dataset
          !// (i.e. all tracers will use it)
          if (INSCAL_LH .ne. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  ': INSCAL_LH /= 0 not allowed for 3D yet!'
             stop
          end if
          if (INSCAL_2D .ne. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  ': INSCAL_2D /= 0 not allowed for 3D!'
             stop
          end if
          if (INSCAL_V .ne. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  ': INSCAL_V /= 0 not allowed for 3D!'
             stop
          end if

       else
          write(6,'(a,i5)') f90file//':'//subr//':A: No such DIM: ',DIM
       end if

    end do !// do S1 = 1,NSETS


    !// 2. Assign tracer-dataset properties
    !// --------------------------------------------------------------------
    do N = 1, ES_NT
       !// Tracer transport number
       TN = ES_TRNR(N)
       !// Initialize scaling factor
       ESCALE = ES_USCALE(N)
       !// Write info
       write(6,'(a,es10.4)') '  * Scaling '//trim(TNAME(TN))// &
            ' for '//CATNAME//' with ',ESCALE

       !// Find all datasets and mark the scalings
       do S1 = 1, NSETS

          !// Index of the emission table for this set/month (S1)
          M = MIDX(S1)

          if (M .le. 0) cycle !// Failsafe: Go to next S1 if M == 0

          !// 2D or 3D field
          if (DIM .eq. 2) then

             if (E2LTBL(TN,M)) then
                write(6,'(a,3i5)') f90file//':'//subr// &
                     ': E2LTBL ALREADY SET!!!',N,TN,M
                stop
             end if
             !// Link tracer TN to table M
             E2LTBL(TN,M) = .true.
             !// Set global scale for tracer TN using table M
             E2STBL(TN,M) = ESCALE
             !// Diagnostic: Total emitted by dataset*scaling
             ADD2EDIAG = sum(E2DS(:,:,1,M)) * ESCALE

          else if (DIM .eq. 3) then

             !// Link tracer TN to table M
             if (E3LTBL(TN,M)) then
                write(6,'(a,3i5)') f90file//':'//subr// &
                     ': E3LTBL ALREADY SET!!!',N,TN,M
                stop
             end if
             E3LTBL(TN,M) = .true.
             !// Set scale for tracer TN using table M
             E3STBL(TN,M) = ESCALE
             !// Diagnostic: Total emitted by dataset*scaling
             ADD2EDIAG = sum(E3DSNEW(:,:,:,M)) * ESCALE

          else
             write(6,'(a,i5)') f90file//':'//subr// &
                  ':B: No such DIM: ',DIM
             stop
          end if

          !// Print out how much was emitted of each tracer/process
          write(6,'(a,a,1x,a,i3,i5,es10.3,es12.5)') &
               '  NAM:CAT:S:M:ESCALE:Tg/mon: ', &
               TNAME(TN),CATNAME,S1,M,ESCALE, ADD2EDIAG*2628000.e-9_r8

       end do !// do S1 = 1,NSETS
    end do !// do N = 1, ES_NT


    !// --------------------------------------------------------------------
  end subroutine emis_setscaling
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine scale_area(R8XY, NSETS, XYBOX, IX,JX,SX,INTYP)
    !// --------------------------------------------------------------------
    !// Scale R8XY with area (xybox is m2). Done before interpolation,
    !// since interpolation routine needs field to be per grid box.
    !//
    !// Ole Amund Sovde, March 2012, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In
    integer, intent(in) :: IX,JX,SX,INTYP,NSETS
    real(r8), intent(in) :: XYBOX(JX)
    !// In/Out
    real(r8), intent(inout) :: R8XY(IX,JX,SX)

    !// Locals
    integer :: I,J,S
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'scale_area'
    !// --------------------------------------------------------------------

    !// INTYP == 0: no area scaling
    if (INTYP .eq. 0) return

    !// Do area scaling
    if (INTYP .eq. 1) then
       !// indata given per cm2
       do S = 1, NSETS
         do j = 1, JX
           do i = 1, IX
              R8XY(I,J,S) = R8XY(I,J,S) * XYBOX(J) * 1.e4_r8
           end do
         end do
       end do
    else if (INTYP .eq. 2) then
       !// indata given per m2
       do S = 1, NSETS
         do j = 1, JX
           do i = 1, IX
              R8XY(I,J,S) = R8XY(I,J,S) * XYBOX(J)
           end do
         end do
        end do
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': INTYP not defined',INTYP
       stop 'STOP in '//subr
    end if

    !// --------------------------------------------------------------------
  end subroutine scale_area
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_xyedges(IRES,JRES,XBEDGE,YBEDGE,START)
    !// --------------------------------------------------------------------
    !// Set up edges for a certain resolution.
    !// START: 1 = Lower-left edge at (0,90S), 2= lower-left edges at
    !//            (180W, 90S)
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IRES,JRES
    integer, intent(in) :: START
    !// Output
    real(r8), intent(out) :: XBEDGE(IRES+1), YBEDGE(JRES+1)

    !// Locals
    real(r8) :: ioffset, joffset, factor
    integer :: I, J
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'get_xyedges'
    !// --------------------------------------------------------------------

    !// Check resolution
    if (IRES.eq.360) then
       
       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 1._r8
          joffset = 91._r8
          factor = 1._r8
       else if (START .eq. 2) then
          !// Grid with (1,1)-box lower-left edges at (180W, 90S)
          ioffset = 181._r8
          joffset = 91._r8
          factor = 1._r8
       else
          write(6,'(a,i5)') f90file//':'//subr//': Unknown map: ',START
          stop
       end if

    else if (IRES.eq.720) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 0.5_r8
          joffset = 90.5_r8
          factor = 0.5_r8
       else if (START .eq. 2) then
          !// Standard 1x1 grid with (1,1)-box lower-left edges at (180W, 90S)
          ioffset = 180.5_r8
          joffset = 90.5_r8
          factor = 0.5_r8
       else
          write(6,'(a,i5)') f90file//':'//subr//': Unknown map: ',START
          stop
       end if

    else if (IRES.eq.1440) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 0.25_r8
          joffset = 90.25_r8
          factor = 0.25_r8
       else if (START .eq. 2) then
          !// Standard 1x1 grid with (1,1)-box lower-left edges at (180W, 90S)
          ioffset = 180.25_r8
          joffset = 90.25_r8
          factor = 0.25_r8
       else
          write(6,'(a,i5)') f90file//':'//subr//': Unknown map: ',START
          stop
       end if

    else if (IRES.eq.3600) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 0.1_r8
          joffset = 90.1_r8
          factor = 0.1_r8
       else if (START .eq. 2) then
          !// Standard 1x1 grid with (1,1)-box lower-left edges at (180W, 90S)
          ioffset = 180.1_r8
          joffset = 90.1_r8
          factor = 0.1_r8
       else
          write(6,'(a,i5)') f90file//':'//subr//': Unknown map: ',START
          stop
       end if

    else if (IRES.eq.160) then

       if (START .eq. 1) then
          !// Lower-left edge at (0,90S)
          ioffset = 1.125_r8
          joffset = 90.1_r8
          factor = 1.125_r8
       else if (START .eq. 2) then
          !// Standard 1x1 grid with (1,1)-box lower-left edges at (180W, 90S)
          ioffset = 181.125_r8
          joffset = 91.125_r8
          factor = 1.125_r8
       else
          write(6,'(a,i5)') f90file//':'//subr//': Unknown map: ',START
          stop
       end if

    else
       write(6,'(a,2i7)') f90file//':'//subr// &
       ': wrong resolution: ',IRES,JRES
       stop
    end if

    do I = 1, IRES+1
       XBEDGE(I) = (real(I, r8) * factor - ioffset)
    end do
    do J = 1, JRES+1
       YBEDGE(J) = (real(J, r8) * factor - joffset)
    end do

    !// --------------------------------------------------------------------
  end subroutine get_xyedges
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_diurnal_scalings()
    !// --------------------------------------------------------------------
    !// Simple hourly diurnal variations may be imposed on category datasets
    !// (monthly). The variations apply for all datasets for a given category.
    !// 
    !// Defining a type of diurnal variations is done in the emission input
    !// file, where also all choices are listed.
    !//
    !// When a new diurnal variation is applied, it must be hardcoded into
    !// this routine.
    !//
    !// Ole Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_oslo, only: E2LocHourSCALE, NECAT, ECATNAMES, &
         E2LocHourRETRO, E2LocHour5050
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Local parameters
    integer :: EC
    !// Diurnal profiles TNO, uncertain but better than nothing.
    !// Local time is index - 1, i.e. starting at 00 local time.
    real(r8), dimension(24) :: &
         POW = (/0.79_r8, 0.72_r8, 0.72_r8, 0.71_r8, 0.74_r8, 0.80_r8, &
                 0.92_r8, 1.08_r8, 1.19_r8, 1.22_r8, 1.21_r8, 1.21_r8, &
                 1.17_r8, 1.15_r8, 1.14_r8, 1.13_r8, 1.10_r8, 1.07_r8, &
                 1.04_r8, 1.02_r8, 1.02_r8, 1.01_r8, 0.96_r8, 0.88_r8/), &
         INC = (/0.75_r8, 0.75_r8, 0.78_r8, 0.82_r8, 0.88_r8, 0.95_r8, &
                 1.02_r8, 1.09_r8, 1.16_r8, 1.22_r8, 1.28_r8, 1.30_r8, &
                 1.22_r8, 1.24_r8, 1.25_r8, 1.16_r8, 1.08_r8, 1.01_r8, &
                 0.95_r8, 0.90_r8, 0.85_r8, 0.81_r8, 0.78_r8, 0.75_r8/), &
         RES = (/0.38_r8, 0.36_r8, 0.36_r8, 0.36_r8, 0.37_r8, 0.50_r8, &
                 1.19_r8, 1.53_r8, 1.57_r8, 1.56_r8, 1.35_r8, 1.16_r8, &
                 1.07_r8, 1.06_r8, 1.00_r8, 0.98_r8, 0.99_r8, 1.12_r8, &
                 1.41_r8, 1.52_r8, 1.39_r8, 1.35_r8, 1.00_r8, 0.42_r8/), &
         SOL = (/0.50_r8, 0.35_r8, 0.20_r8, 0.10_r8, 0.10_r8, 0.20_r8, &
                 0.75_r8, 1.25_r8, 1.40_r8, 1.50_r8, 1.50_r8, 1.50_r8, &
                 1.50_r8, 1.50_r8, 1.50_r8, 1.50_r8, 1.50_r8, 1.40_r8, &
                 1.25_r8, 1.10_r8, 1.00_r8, 0.90_r8, 0.80_r8, 0.70_r8/), &
         TRA = (/0.19_r8, 0.09_r8, 0.06_r8, 0.05_r8, 0.09_r8, 0.22_r8, &
                 0.86_r8, 1.84_r8, 1.86_r8, 1.41_r8, 1.24_r8, 1.20_r8, &
                 1.32_r8, 1.44_r8, 1.45_r8, 1.59_r8, 2.03_r8, 2.08_r8, &
                 1.51_r8, 1.06_r8, 0.74_r8, 0.62_r8, 0.61_r8, 0.44_r8/), &
         AGR = (/0.60_r8, 0.60_r8, 0.60_r8, 0.60_r8, 0.60_r8, 0.65_r8, &
                 0.75_r8, 0.90_r8, 1.10_r8, 1.25_r8, 1.45_r8, 1.60_r8, &
                 1.75_r8, 1.75_r8, 1.70_r8, 1.55_r8, 1.35_r8, 1.1_r8, &
                 0.90_r8, 0.75_r8, 0.65_r8, 0.60_r8, 0.60_r8, 0.60_r8/)
    !// --------------------------------------------------------------------


    !// Initialize scalings to 1
    E2LocHourSCALE(:,:,:) = 1._r8

    !// Scalings for RETRO (ESCAL_LH = 1)
    !// --------------------------------------------------------------------
    !// If no RETRO categories are found, the variations will not be
    !// initialized. This is of no problem, as the CTM3 cannot use
    !// non-defined categories anyway.
    do EC = 1, NECAT
       if (ECATNAMES(EC).eq.'POW') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = POW(:)
       else if (ECATNAMES(EC).eq.'ENE') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = POW(:)
       else if (ECATNAMES(EC).eq.'INC') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = INC(:)
       else if (ECATNAMES(EC).eq.'IND') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = INC(:)
       else if (ECATNAMES(EC).eq.'RES') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = RES(:)
       else if (ECATNAMES(EC).eq.'RCO') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = RES(:)
       else if (ECATNAMES(EC).eq.'SOL') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = SOL(:)
       else if (ECATNAMES(EC).eq.'SLV') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = SOL(:)
       else if (ECATNAMES(EC).eq.'TRA') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = TRA(:)
       else if (ECATNAMES(EC).eq.'AGR') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = AGR(:)
       else if (ECATNAMES(EC).eq.'AWB') then
          E2LocHourSCALE(:,EC,E2LocHourRETRO) = AGR(:)
       end if
    end do

    !// Scalings for 5050 (ESCAL_LH = 2)
    !//   50% reduction from 8pm to 7am
    !//   50% increase from 8am to 7pm
    !//   50% reduction from 8pm to 7am
    !// --------------------------------------------------------------------
    E2LocHourSCALE( 1: 7,:,E2LocHour5050) = 0.5_r8
    E2LocHourSCALE( 8:19,:,E2LocHour5050) = 1.5_r8
    E2LocHourSCALE(20:24,:,E2LocHour5050) = 0.5_r8

    !// --------------------------------------------------------------------
  end subroutine set_diurnal_scalings
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_vertical_scalings()
    !// --------------------------------------------------------------------
    !// Simple vertical scalings for 2D emissions.
    !// Factors for distributing surface emissions into the lowermost
    !// model levels.
    !//
    !// Ole Amund Sovde, November 2015
    !// --------------------------------------------------------------------
    use cmn_oslo, only: NE2vertLVS, E2vertSCALE, E2vertTBL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Vertical scalings
    E2vertSCALE(:,:) = 0._r8
    E2vertTBL(:) = 0
    !// Average layer thickness [m] for ECMWF L1-L4:
    !// 0-16 / 16-41 / 41-77 / 77-128

    !// Scaling typically used for power and industrial combustion
    E2vertSCALE(1:NE2vertLVS,1) = (/ 0.25_r8, 0.125_r8, 0.125_r8, 0.5_r8 /)

    !// Scaling typically used for residential heating
    E2vertSCALE(1:NE2vertLVS,2) = (/ 0.60_r8, 0.200_r8, 0.200_r8, 0.0_r8 /)

    !// Scaling typically used for ship emissions
    E2vertSCALE(1:NE2vertLVS,3) = (/ 0.30_r8, 0.400_r8, 0.300_r8, 0.0_r8 /)

    !// --------------------------------------------------------------------
  end subroutine set_vertical_scalings
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emis_setscaling_2dfields()
    !// --------------------------------------------------------------------
    !// Update 2D scaling datasets.
    !//
    !// Ole Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, MPBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, JMON, &
         XDEDG, YDEDG, AREAXY
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: E22dSCALE, NE22dVARS, E2L2dTBL
    use ncutils, only: get_netcdf_var_3d, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    !// Output

    !// Locals
    integer :: N, I,J,II,JJ, MP
    integer,parameter :: I1x1=360, J1x1=180
    real(r8) :: RIN(I1x1,J1x1,12), RIN2(I1x1,J1x1,12),EBOX(I1x1,J1x1)
    real(r8) :: XBEDGE(I1x1+1), YBEDGE(J1x1+1), XYBOX(J1x1)
    real(r8) :: RRUN(IPAR,JPAR)

    character(len=120) :: INFILE
    character(len=4) :: CYEAR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_setscaling_2dfields'
    !// --------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr// &
         ': Updating 2D fields for emission scalings'
    do N = 1, NE22dVARS
       if (E2L2dTBL(N)) then
          select case(N)
          case(1)
             !// Read scaling for daylight temperatures
             !// read dataset 1
             infile = 'Indata_CTM3/emissions_scalings_t_light.nc'
             call get_netcdf_var_3d( infile,'tempscale', RIN, I1x1,J1x1,12 )
          case(2)
             !// read scaling for daylight
             infile = 'Indata_CTM3/emissions_scalings_t_light.nc'
             call get_netcdf_var_3d( infile,'lightscale', RIN, I1x1,J1x1,12 )
          case(3)
             !// read scaling for domestic biofuel
             write(CYEAR(1:4),'(i4.4)') MYEAR
             infile = 'Indata_CTM3/emissions_scalings_hddfrac_'//CYEAR//'.nc'
             write(6,'(a)') 'Reading: '//trim(infile)
             call get_netcdf_var_3d( infile,'hddfrac', RIN, I1x1,J1x1,12 )
             !// The new RIN is monthly fraction of annual total. When the input
             !// emissions are kg/s, we need to multiply by 12 to get correct
             !// scaling.
             RIN(:,:,:) = RIN(:,:,:) * 12._r8
          case default
             !// No such dataset
             write(6,'(a,i3)') f90file//':'//subr// &
                  ': No such 2DSCAL dataset: ',N
             write(6,'(a)') '  STOP in '//subr
             stop
          end select

          call get_xyedges(I1x1,J1x1,XBEDGE,YBEDGE,1)

          !// Grid box areas
          do J=1,J1x1
             XYBOX(J) =  A0*A0 * CPI180*(XBEDGE(2)-XBEDGE(1)) &
                  * (sin(CPI180*YBEDGE(J+1)) - sin(CPI180*YBEDGE(J)))
          end do
          !// Multiply to get per area
          do J = 1, J1x1
             EBOX(:,J) = RIN(:,J,JMON)*XYBOX(J)
          end do

          !// Interpolate
          call E_GRID(EBOX(:,:),XBEDGE,YBEDGE,I1x1,J1x1, &
                      RRUN,XDEDG,YDEDG,IPAR,JPAR,1)

          !// Divide by area
          RRUN(:,:) = RRUN(:,:) / AREAXY(:,:)

          !// Put into array
          do MP = 1, MPBLK
            do J = MPBLKJB(MP), MPBLKJE(MP)
              JJ   = J - MPBLKJB(MP) + 1
              do I = MPBLKIB(MP), MPBLKIE(MP)
                II   = I - MPBLKIB(MP) + 1
                E22dSCALE(II,JJ,MP,N) = RRUN(I,J)
              end do
            end do
          end do

          write(6,'(a,i3)') '* Updated 2D scalings: ',N

       end if
    end do
    write(6,'(a)') f90file//':'//subr// &
         ': Done updating 2D fields for emission scalings'
    write(6,'(a)') '-------------------------------------------------' &
         //'----------------------'

    !// --------------------------------------------------------------------
  end subroutine emis_setscaling_2dfields
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emis_diag(BEMIS,ORG_DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Diagnose emissions using STTTND and STTTN0 (cmn_d.f).
    !// Because routine is called in IJ-block, we need to fetch it
    !// into p-main. Then it added to STTTND and STTTN0.
    !//
    !// Routine is called ONLY every NOPS.
    !//
    !// Ole Amund Sovde, September 2012
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, NRCHEM
    use cmn_oslo, only: DIAGEMIS_IJ
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)  :: MP                         !// IJ-block index
    real(r8), intent(in) :: ORG_DTCHM                  !// Org. chem time step
    real(r8), intent(in) :: BEMIS(LPAR,NPAR,IDBLK,JDBLK) !// Grid box emissions

    !// Locals
    real(r8) :: DT
    integer :: I,J,II,JJ,N,L
    !// --------------------------------------------------------------------

    DT = ORG_DTCHM * real(NRCHEM,r8)

    !// Add BEMIS to DIAGEMIS_IJ
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ   = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II   = I - MPBLKIB(MP) + 1
          do N = 1, NPAR
             do L = 1, LPAR
                DIAGEMIS_IJ(L,N,II,JJ,MP) = DIAGEMIS_IJ(L,N,II,JJ,MP) &
                     + BEMIS(L,N,II,JJ) * DT
             end do
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine emis_diag
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine ctm2_rdmolec2(R8XY,FILENAME,EFNR,IEMIS,JEMIS, &
       MRES,NTIMES)
    !// --------------------------------------------------------------------
    !// Reads CTM2 POET format. Based on CTM2 readmolec.
    !//   NTIMES: Number of months in file.
    !// 
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IEMIS,JEMIS,MRES,NTIMES
    character(len=*), intent(in) :: FILENAME
    !// Output
    Real(r8), intent(out)       :: R8XY(IEMIS,JEMIS,MRES)

    !// Locals
    real(r8) :: out(NTIMES),  sgrid,tgrid,sgrid2,tgrid2
    integer :: ilevel,igrid,jgrid,itime
    Integer :: I,J,K,M
    !// --------------------------------------------------------------------

    open (efnr,file=FILENAME,status='unknown',form='formatted')
    read (efnr, 103)
103 format (//)

    do I = 1, 64800    !  (360x180=64800 gridboxes)
       !// ... read in all levels and times for one grid

       read (efnr,*) sgrid,sgrid2,tgrid,tgrid2, &
            (out(itime),itime=1,NTIMES)

       !// Adjust grid so that we start at Greenwich and South Pole
       igrid = INT(sgrid)
       jgrid = INT(tgrid)

       if (igrid .lt. 0) Then
          igrid = 360 - abs(igrid) + 1
       else
          igrid = igrid + 1
       end if
       if (jgrid .lt. 0) Then
          jgrid = 90 - abs(jgrid) + 1
       else
          jgrid = 90 + jgrid + 1   
       end if

       !// move data into array
       do K = 1, NTIMES
          R8XY(IGRID,JGRID,K) = out(K)
       end do
      
    end do

    close(efnr)

    !// --------------------------------------------------------------------
  end subroutine ctm2_rdmolec2
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine iiasa_read2d(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Read annual IIASA 2D data and put into R8XY.
    !// NSETS is assumed to be 1 for now.
    !//
    !// Ole Amund Sovde, August 2009, March 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IRES,JRES,MRES,NSETS
    character(len=*), intent(in) :: INFILE
    !// Output
    real(r8), intent(out)       :: R8XY(IRES,JRES,MRES)

    !// Locals
    real(r8)  :: RDUM
    integer :: I,J,K, II,JJ, io_err
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'iiasa_read2d'
    !// --------------------------------------------------------------------

    !// Open file
    open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems with file: '//trim(infile)
       write(6,'(a,i5)') 'io_err',io_err
       stop
    end if

    !// Initialize
    R8XY(:,:,:) = 0._r8

    read(efnr,*)
    read(efnr,*)
    read(efnr,*)

    do K = 1, 64800   ! 360x180 max (actually 967)
       read(efnr,'(2(1x,i4),1(1x,e10.4))',iostat=io_err) I,J,RDUM
       if (io_err.ne.0) exit
       !// File starts at -180 E, 90 S, -179 E, 90 S etc.
       !// Turn the array from -180->180 into 1->360 & -90->90 into 1->180
       if (I .lt. 0) then       ! W of Greenwich
          II = I + 361
       else
          II = I + 1
       end if                 ! if (I .lt. 0)
       JJ = J + 91          ! 1->91
       R8XY(II,JJ,1) = RDUM

    end do

    !// Close file
    close(efnr)
  
    !// --------------------------------------------------------------------
  end subroutine iiasa_read2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_oldNH3(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS,ANTR)
    !// --------------------------------------------------------------------
    !// Read monthly NH3 emissions used in CTM2, with one exception:
    !// Biomass burning is not included in CTM3-files, because it is
    !// included in GFEDv3.
    !//
    !// Ole Amund Sovde, March 2012, November 2011
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IRES,JRES,MRES,NSETS,ANTR
    character(len=*), intent(in) :: INFILE
    !// Output
    real(r8), intent(out)       :: R8XY(IRES,JRES,MRES)

    !// Locals
    real(r8)  :: RDUM
    integer :: I,J,M,io_err
    character(len=80) :: filename
    character(len=3), parameter, dimension(12) :: AMON = ['JAN','FEB', &
         'MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_oldNH3'
    !// --------------------------------------------------------------------

    !// Initialize
    R8XY(:,:,1:NSETS) = 0._r8

    do M = 1, NSETS

       if (ANTR .gt. 0) then
          filename = trim(INFILE)//'NH3emis_noBB_'//AMON(M)//'.1x1'
       else
          !filename = trim(INFILE)//'NH3NOMAN'//AMON(M)//'.1x1'
          filename = trim(INFILE)//'NH3emis_NAT_noBB_'//AMON(M)//'.1x1'
       end if

       !// Open file
       open(efnr,file=filename,status='old',form='formatted',iostat=io_err)
       if (io_err .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Problems with file: '//trim(infile)
          write(6,'(a,i5)') 'io_err',io_err
          stop
       end if

       !// Reading in file in kg(N)/month pr gridbox
       do while (.true.)
          read(efnr,*,iostat=io_err) I,J,R8XY(I,J,M)
          if (io_err .ne. 0) exit
       end do

       !// Close file
       close (efnr)

    end do !// do M = 1, NSETS

    !// --------------------------------------------------------------------
  end subroutine read_oldNH3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine gfed_read3d(R8XYZ_IN,INFILE,efnr,IRES,JRES,L1x1,MRES,NSETS,NLEVS)
    !// --------------------------------------------------------------------
    !// Read annual GFED 3D data and put into R8XYZ_IN.
    !// NSETS is assumed to be 1 for now.
    !//
    !// Ole Amund Sovde, August 2009, March 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IRES,JRES,MRES,L1x1,NSETS
    character(len=*), intent(in) :: INFILE
    !// Output
    real(r8), intent(out)       :: R8XYZ_IN(IRES,JRES,L1x1,MRES)
    integer, intent(out)      :: NLEVS

    !// Locals
    real(r8)  :: RDUM(MRES)
    integer :: I,J,K,L,M, io_err
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed_read3d'
    !// --------------------------------------------------------------------

    !// Open file
    open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems with file: '//trim(infile)
       write(6,'(a,i5)') 'io_err',io_err
       stop
    end if

    !// Initialize
    R8XYZ_IN(:,:,:,1:NSETS)=0._r8
    !// Only one level
    NLEVS = 1

    read(efnr,*)
    read(efnr,*)
    read(efnr,*)

    do K = 1, 99999   ! 360x180 max (actually 967)
       read(efnr,'(3(1x,i4),12(1x,e10.4))',iostat=io_err) &
            I, J, L, (RDUM(M), M=1,NSETS)
       if (io_err .ne. 0) exit
       !// File starts at -180 E, 90 S, -179 E, 90 S etc.
       !// Turn the array from -180->180 into 1->360 & -90->90 into 1->180
       if (I .lt. 0) then     ! W of Greenwich
          I = I + 361
       else
          I = I + 1
       end if                 ! if (I .lt. 0)
       J = J + 91             ! 1->91
       R8XYZ_IN(I,J,L,1:NSETS) = RDUM(1:NSETS)

    end do

    !// Close file
    close(efnr)
  
    !// --------------------------------------------------------------------
  end subroutine gfed_read3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine volc_read1(INFILE,efnr,IRES,JRES,MRES,NSETS,MTBL,MIDX, &
       XBEDGE, YBEDGE, XDEDG, YDEDG)
    !// --------------------------------------------------------------------
    !// Read annual volcano 3D data and put into emission array.
    !// Only one category; natural.
    !// Update MTBL and MIDX.
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: ETAA, ETAB
    use cmn_chem, only: E3DSNEW, E3PAR, ETPAR
    use cmn_met, only: PMEAN, T, Q
    use regridding, only: E_GRID
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IRES,JRES,MRES,NSETS
    character(len=*), intent(in) :: INFILE
    real(r8), intent(in)        :: XBEDGE(IRES+1), YBEDGE(JRES+1)
    real(r8), intent(in)        :: XDEDG(IPAR+1), YDEDG(JPAR+1)
    !// Output
    integer, intent(inout)    :: MTBL
    integer, intent(out)      :: MIDX(MRES)

    !// Locals
    real(r8) :: VMASS1x1(IRES,JRES), VMINH1x1(IRES,JRES), VMAXH1x1(IRES,JRES)
    real(r8) :: VMASS(IPAR,JPAR), VMINH(IPAR,JPAR), VMAXH(IPAR,JPAR)
    real(r8)  :: RDUM,v1,v2,v3,dvdh,dz,zh(lpar+1)
    integer :: I,J,K,L,Lmin,Lmax, io_err
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'volc_read1'
    !// --------------------------------------------------------------------

    !// Open file
    open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems with file: '//trim(infile)
       write(6,'(a,i5)') 'io_err',io_err
       stop
    end if

    !// Read off header
    read(efnr,*)
    read(efnr,*)
    read(efnr,*)

    !// Read data
    do K = 1, 99999
       read(efnr,'(2(1x,i4),3(1x,e10.4))',iostat=io_err) I,J,v1,v2,v3
       if (io_err.ne.0) exit

       !// File starts at -180 E, 90 S, -179 E, 90 S etc.
       !// Turn the array from -180->180 into 1->360 & -90->90 into 1->180
       if (I .lt. 0) then     !// W of Greenwich
          I = I + 361
       else
          I = I + 1
       end if                 !// if (I .lt. 0)
       J = J + 91             !// 1->91
 
       !// Mass into 2D array
       vmass1x1(i,j) = v1
       !// Min height into 2D array
       vminh1x1(i,j) = v2
       !// Max height into 2D array
       vmaxh1x1(i,j) = v3
    end do

    close(efnr)

    !// Update 3D table index
    MTBL = MTBL + 1
    if (MTBL .gt. E3PAR) then
       write(6,'(a,2i5)') f90file//':'//subr// &
            ': MTBL > E3PAR', MTBL, E3PAR
       stop
    end if
    !// Save table index for this dataset (only one)
    MIDX(1) = MTBL

    !// Interpolate all arrays
    call E_GRID(vmass1x1,XBEDGE,YBEDGE,IRES,JRES, &
                  vmass,XDEDG,YDEDG,IPAR,JPAR,1)
    call E_GRID(vminh1x1,XBEDGE,YBEDGE,IRES,JRES, &
                  vminh,XDEDG,YDEDG,IPAR,JPAR,1)
    call E_GRID(vmaxh1x1,XBEDGE,YBEDGE,IRES,JRES, &
                  vmaxh,XDEDG,YDEDG,IPAR,JPAR,1)

    !// Find Z and DZ for each T42 column and calculate emissions
    do j = 1, jpar
       do i = 1, ipar
          if (vmass(i,j).gt.0._r8) then

             !// Height of lower edges
             zh(1) = 0._r8
             do L = 1, LPAR
                dz = - 29.26586_r8 * T(I,J,L) * (1._r8 + 0.6_r8*Q(I,J,L)) &
                    *log((etaa(L+1)+etab(L+1)*pmean(i,j)) &
                          / (etaa(L)+etab(L)*pmean(i,j)))
                zh(L+1) = zh(L) + dz
             end do

             !// Gradient for emissions with height
             dvdh = vmass(i,j) / (vmaxh(i,j) - vminh(i,j))
             do L = 1, LPAR
                if (vminh(i,j).lt.zh(l+1)) exit
             end do
             Lmin = min(L,LPAR)
             !// Lmin+1 is the first level confined in the volcano plume,
             !// and the plume starts below this level (i.e. in layer Lmin).

             !// Volcano top is therefore in layer Lmin+1 or above. This test
             !// uses Lmax=Lmin+1 when if-test fails.
             !do L=LPAR,Lmin+2,-1
             !   if (vmaxh(i,j).ge.zh(l)) exit
             !end do
             do L = Lmin, LPAR-1
                if (zh(l+1) .gt. vmaxh(i,j)) exit
             end do
             Lmax = min(L,LPAR)

             !// Put results directly into E3DS
             !// Bottom volcano layer, from vminh to zh(Lmin+1) into layer Lmin
             !E3DS(i,j,Lmin,MTBL) = E3DS(i,j,Lmin,MTBL) &
             E3DSNEW(Lmin,i,j,MTBL) = E3DSNEW(Lmin,i,j,MTBL) &
                    + dvdh * (zh(Lmin+1) - vminh(i,j))
             if ((zh(Lmin+1) - vminh(i,j)).lt.0._r8) then
                write(6,'(a,4i5,3f9.2)') f90file//':'//subr// &
                     ':A:',lmin,i,j,mtbl,dvdh,zh(Lmin+1),vminh(i,j)
                stop
             end if
             !// Iside the volcano plume, layer thickness is zh(L+1) - zh(L)
             do L = Lmin+1, (Lmax-1)
                !E3DS(i,j,L,MTBL) = E3DS(i,j,L,MTBL) &
                E3DSNEW(L,i,j,MTBL) = E3DSNEW(L,i,j,MTBL) &
                    + dvdh * (zh(L+1) - zh(L))
                if ((zh(L+1) - zh(L)).lt.0._r8) then
                   write(6,'(a,4i5,3f9.2)') f90file//':'//subr// &
                        ':B:',l,i,j,mtbl,dvdh,zh(L+1),zh(l)
                   stop
                end if
             end do
             !// Top, from zh(Lmax) to top of plume vmaxh
             !E3DS(i,j,Lmax,MTBL) = E3DS(i,j,Lmax,MTBL) &
             E3DSNEW(Lmax,i,j,MTBL) = E3DSNEW(Lmax,i,j,MTBL) &
                    + dvdh * (vmaxh(i,j) - zh(Lmax))
             if ((vmaxh(i,j) - zh(Lmax)).lt.0._r8) then
                write(6,'(a,4i5,3f9.2)') f90file//':'//subr// &
                ':C:',lmax,i,j,mtbl,dvdh,vmaxh(i,j),zh(Lmax)
                stop
             end if
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine volc_read1
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_retrobbh(INFILE,IHLF,JHLF,LMAX,BBH,LPARW,LMMAP)
    !// --------------------------------------------------------------------
    !// Read biomass burning height distribution.
    !// Flips array so it starts on (0E,90S) instead of (180W,90N).
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)       :: IHLF,JHLF,LMAX,LPARW,LMMAP(LPARW)
    character(len=*), intent(in) :: INFILE
    real(r8), intent(out)       :: BBH(IHLF,JHLF,LMAX)

    integer, parameter :: LHLF=13 !// Levels of input data
    integer :: I,J, L, nLon, nLat, nLev
    real(r8) :: BFRAC(IHLF,JHLF,LHLF),BBFRAC(IHLF,JHLF,LHLF), vsum
    character(len=60) :: FIELDNAME
    real(r8), allocatable, dimension(:) :: inLon, inLat, inLev
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_retrobbh'
    !// --------------------------------------------------------------------
    FIELDNAME = 'altitude_fraction'

    !call readnc_field(BBFRAC,FIELDNAME,LONG,LAT,LEVEL,INFILE,IHLF,JHLF,LHLF)

    !// Check resolution (latitude/longitude)
    !// This routine allocates inLon/inLat
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'level',  inLev  )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nLev  = SIZE( inLev  )

    deallocate( inLon, inLat, inLev )

    if (nLon .ne. IHLF .or. nLat.ne.JHLF .or. nLev.ne.LHLF) then
       write(6,'(a,3i5)') f90file//':'//subr// &
            ': wrong resolution ', nLon, nLat, nLev
       stop
    end if

    !// Get data, should be resolution 720x360x13
    call get_netcdf_var_3d(infile, fieldname, BBFRAC, nLon, nLat, nLev)


    !// Flip array so it starts on (0E,90S)
    do I = 1, IHLF/2
       BFRAC(IHLF/2+I,:,:)=BBFRAC(I,:,:)
    end do
    do I = (IHLF/2)+1, IHLF
       BFRAC(I-IHLF/2,:,:)=BBFRAC(I,:,:)         
    end do
    do J = 1, JHLF
       BBFRAC(:,JHLF-J+1,:)=BFRAC(:,J,:)
    end do

    !// Initialize output field
    BBH(:,:,:)=0._r8

    !// Do vertical adjustments
    BBH(:,:,LMMAP( 1)) = BBH(:,:,LMMAP( 1)) + 0.3_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 2)) = BBH(:,:,LMMAP( 2)) + 0.2_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 3)) = BBH(:,:,LMMAP( 3)) + 0.1_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 4)) = BBH(:,:,LMMAP( 4)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 5)) = BBH(:,:,LMMAP( 5)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 6)) = BBH(:,:,LMMAP( 6)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 7)) = BBH(:,:,LMMAP( 7)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 8)) = BBH(:,:,LMMAP( 8)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP( 9)) = BBH(:,:,LMMAP( 9)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP(10)) = BBH(:,:,LMMAP(10)) + 0.05_r8*BBFRAC(:,:,1)     
    BBH(:,:,LMMAP(11)) = BBH(:,:,LMMAP(11)) + 0.05_r8*BBFRAC(:,:,1)
    BBH(:,:,LMMAP(12)) = BBH(:,:,LMMAP(12)) + 0.25_r8*BBFRAC(:,:,2)
    BBH(:,:,LMMAP(13)) = BBH(:,:,LMMAP(13)) + 0.25_r8*BBFRAC(:,:,2)
    BBH(:,:,LMMAP(14)) = BBH(:,:,LMMAP(14)) + 0.25_r8*BBFRAC(:,:,2)
    BBH(:,:,LMMAP(15)) = BBH(:,:,LMMAP(15)) + 0.25_r8*BBFRAC(:,:,2)
    BBH(:,:,LMMAP(16)) = BBH(:,:,LMMAP(16)) + 0.5_r8*BBFRAC(:,:,3)
    BBH(:,:,LMMAP(17)) = BBH(:,:,LMMAP(17)) + 0.5_r8*BBFRAC(:,:,3)
    BBH(:,:,LMMAP(18)) = BBH(:,:,LMMAP(18)) + (1._r8/3._r8)*BBFRAC(:,:,4)
    BBH(:,:,LMMAP(19)) = BBH(:,:,LMMAP(19)) + (1._r8/3._r8)*BBFRAC(:,:,4)
    BBH(:,:,LMMAP(20)) = BBH(:,:,LMMAP(20)) + (1._r8/3._r8)*BBFRAC(:,:,4)
    BBH(:,:,LMMAP(21)) = BBH(:,:,LMMAP(21)) + 0.5_r8*BBFRAC(:,:,5)
    BBH(:,:,LMMAP(22)) = BBH(:,:,LMMAP(22)) + 0.5_r8*BBFRAC(:,:,5)
    BBH(:,:,LMMAP(23)) = BBH(:,:,LMMAP(23)) + 0.5_r8*BBFRAC(:,:,6)
    BBH(:,:,LMMAP(24)) = BBH(:,:,LMMAP(24)) + 0.5_r8*BBFRAC(:,:,6)
    BBH(:,:,LMMAP(25)) = BBH(:,:,LMMAP(25)) + BBFRAC(:,:,7)
    BBH(:,:,LMMAP(26)) = BBH(:,:,LMMAP(26)) + 0.5_r8*BBFRAC(:,:,8)
    BBH(:,:,LMMAP(27)) = BBH(:,:,LMMAP(27)) + 0.5_r8*BBFRAC(:,:,8)
    BBH(:,:,LMMAP(28)) = BBH(:,:,LMMAP(28)) + BBFRAC(:,:,9)
    BBH(:,:,LMMAP(29)) = BBH(:,:,LMMAP(29)) + BBFRAC(:,:,10)
    BBH(:,:,LMMAP(30)) = BBH(:,:,LMMAP(30)) + BBFRAC(:,:,11)
    BBH(:,:,LMMAP(31)) = BBH(:,:,LMMAP(31)) + BBFRAC(:,:,12)
    BBH(:,:,LMMAP(32)) = BBH(:,:,LMMAP(32)) + 0.5_r8*BBFRAC(:,:,13)
    BBH(:,:,LMMAP(33)) = BBH(:,:,LMMAP(33)) + 0.5_r8*BBFRAC(:,:,13)


    !// Check BB
    do J = 1, JHLF
       do I = 1, IHLF
          if (BBH(I,J,1) .eq. 0._r8 .and. sum(BBH(I,J,:)).gt.0._r8) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': BBH is wrong'
             do L = 1, LMAX
                write(6,'(i5,es12.3)') L, BBH(I,J,L)
             end do
             stop
          end if
       end do
    end do

    !// Will adjust vertical sum so it is either 0 or 1.
    !// (Due to single precision on file, the sum in double precision
    !//  is often not exactly 1.)
    do J = 1, JHLF
       do I = 1, IHLF
          vsum = sum(BBH(I,J,:))
          if (vsum .gt. 0._r8 .and. vsum .ne. 1._r8) &
               BBH(I,J,:) = BBH(I,J,:) / vsum
       end do
    end do


    !// --------------------------------------------------------------------
  end subroutine read_retrobbh
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_bcoc_bond_2d(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Read annual 2D data for BCOC (produced by Tami Bond), and put into
    !// layer 1 of R8XY. NSETS is therefore 1.
    !//
    !// Reads both VERSION 4 (bc1-files) and 5 (BC1EN-files).
    !//
    !// Flips dataset to grid_type 1.
    !//
    !// Ole Amund Sovde, September 2009, March 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)       :: EFNR,IRES,JRES,MRES,NSETS
    character(len=*), intent(in) :: INFILE
    !// Output
    real(r8), intent(out)       :: R8XY(IRES,JRES,MRES)

    !// Locals
    real(r8)  :: RDUM, Input10(10), RTEMP(IRES,JRES)
    integer :: IDEG, JDEG, I1, I2, J, N, io_err, VERSION
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'volc_read1'
    !// --------------------------------------------------------------------

    !// Open file
    open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems with file: '//trim(infile)
       write(6,'(a,i5)') 'io_err',io_err
       stop
    end if

    !// Initialize
    R8XY(:,:,1:NSETS)=0._r8

    !// Read header info (9 lines)
    read(efnr,'(/)')
    read(efnr,'(14x,i1)') version
    write(6,'(a,i5)') ' * Data from Tami Bond, VERSION',version
    read(efnr,'(/////)')

    !// Weird looping, could simplify but decided to keep the main loops.
    do JDEG = -89, 90
       !// For CTM 1x1 grid, starts at 1
       J = JDEG + 90
       do IDEG = -17, 18
         !// Read 10 longitudes at a time
         if (VERSION .eq. 4) then
            read(efnr,'(10(F8.3,1X))') (Input10(N),N=1,10)
         else
            read(efnr,'(10(F10.3,1X))') (Input10(N),N=1,10)
         end if
         I1 = (IDEG * 10) - 9
         I2 = (IDEG * 10)
         !// For CTM 1x1 grid, shift -179:0 to 181:360
         if (I1.lt.1) I1 = I1 + 360
         if (I2.lt.1) I2 = I2 + 360
         R8XY(I1:I2,J,1) = Input10(:)
       end do
    end do

    !// Close file
    close(efnr)

    !// --------------------------------------------------------------------
  end subroutine read_bcoc_bond_2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_aerocom_2d(R8XY,INFILE,EFNR,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Read Secondary Organic aerosol data for 12 months, used in the BCOC
    !// module. Data are from AEROCOM, and should not be used when a more
    !// sophisticated SOA scheme is used. NSETS=12.
    !//
    !// Flips dataset to grid_type 1.
    !//
    !// Ole Amund Sovde, September 2009, March 2012
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    Integer, intent(in)       :: EFNR,IRES,JRES,MRES,NSETS
    Character(len=*), intent(in) :: INFILE
    !// Output
    Real(r8), intent(out)       :: R8XY(IRES,JRES,MRES)

    !// Locals
    real(r8)  :: RDUM, Input12(12), RTEMP(IRES,JRES)
    integer :: IBOX, JBOX, I, J, K, M, io_err
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'volc_aerocom_2d'
    !// --------------------------------------------------------------------

    !// Open file
    open(efnr,file=INFILE,status='old',form='formatted',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems with file: '//trim(infile)
       write(6,'(a,i5)') 'io_err',io_err
       stop
    end if

    !// Initialize
    R8XY(:,:,1:NSETS) = 0._r8

    !// Read header info (3 lines)
    read(efnr,'(//)')

    !// Read data
    do K = 1, 64800   ! 360x180 max (actually 16440)
       read(efnr,'(2(1x,i4),12(1x,e10.4))',iostat=io_err) &
            I,J,(Input12(M),M=1,12)
       if (io_err .ne. 0) exit
       !// Data starts at 180E,90S, assume edge values
       !// J=-90 is south edge of gridbox 1 (add 91)
       !// I=-180 is west edge of gridbox 181 (add 361)
       !//  = 0 is west of gridbox 1 (add 1)
       if (I .lt. 0) then
          IBOX = I + 361
       else
          IBOX = I + 1
       end if
       JBOX = J + 91
       do M = 1, 12
          R8XY(IBOX,JBOX,M) = Input12(M)
       end do

    end do

    !// Close file
    close(efnr)

    !// --------------------------------------------------------------------
  end subroutine read_aerocom_2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine nc_getdata_sfc(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Reads surface emissions for 1 or 12 months, from netCDF file, 
    !// fetching variable name INVAR.
    !//
    !// Based on routines from Chris D. Holmes, 2012.
    !//
    !// Ole Amund Sovde, March 2012
    !// --------------------------------------------------------------------
    use netcdf
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
 
    !// Input
    character(len=*), intent(in)  :: infile   !// Name of netCDFfile
    character(len=*), intent(in)  :: invar    !// Name of input variable
    integer,intent(in)   :: IRES     !// Number of longitudes
    integer,intent(in)   :: JRES     !// Number of latitudes
    integer,intent(in)   :: MRES     !// Number of sets (months) per category

    !// Output
    real(r8), intent(out) :: R8XY(IRES,JRES,MRES)  !// 3D field
    integer,intent(out)   :: NSETS

    !// Locals
    integer            :: status, ncid, nDims, nVars, nAtts
    integer            :: nLon, nLat, nLev, nTime
    character(len=25)  :: units !CDH caution using fixed length
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nc_getdata_sfc'
    !// --------------------------------------------------------------------


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    !// Deallocate all local variables
    IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
    IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
    IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)

    if (nLon.ne.IRES .or. nLat.ne.JRES) then
       write(6,'(a)') f90file//':'//subr// &
            ': Horizontal resolution on file does not match'// &
            ' specified resolution'
       write(6,'(a,2i5)') '  Specified:',IRES,JRES
       write(6,'(a,2i5)') '  On file:  ',nLon,nLat
       stop
    end if
    if (nTime.eq.1) then
       write(6,'(a)') f90file//':'//subr//': One dataset'
    else if (nTime.eq.12) then
       write(6,'(a)') f90file//':'//subr//': 12 datasets'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': not 1 or 12 datasets - do not know what to do'
       stop
    end if
       
    !// Emission variable, expected units: kg/m2/s for each month
    call get_netcdf_var_3d( infile, invar, R8XY, nlon,nlat,ntime )

    !// Number of datasets
    NSETS = nTime

    !// --------------------------------------------------------------------
  end subroutine nc_getdata_sfc
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  subroutine nc_getdata_sfc_2d(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Reads surface emissions for 1 time step, from netCDF file, 
    !// fetching variable name INVAR.
    !//
    !// Based on routines from Chris D. Holmes, 2012.
    !//
    !// Ole Amund Sovde, September 2012
    !// --------------------------------------------------------------------
    use netcdf
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
 
    !// Input
    character(len=*), intent(in)  :: infile   !// Name of netCDFfile
    character(len=*), intent(in)  :: invar    !// Name of input variable
    integer,intent(in) :: IRES     !// Number of longitudes
    integer,intent(in) :: JRES     !// Number of latitudes
    integer,intent(in) :: MRES     !// Number of sets (months) per category

    !// Output
    real(r8), intent(out) :: R8XY(IRES,JRES,MRES)  !// 3D field
    integer,intent(out)   :: NSETS

    !// Locals
    integer            :: status, ncid, nDims, nVars, nAtts
    integer            :: nLon, nLat, nLev, nTime, I,J
    character(len=25)  :: units !CDH caution using fixed length
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nc_getdata_sfc_2d'
    !// ------------------------------------------------------------------


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = 1

    !// Deallocate all local variables
    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)
    if ( allocated(inTime) ) deallocate(inTime)

    if (nLon.ne.IRES .or. nLat.ne.JRES) then
       write(6,'(a)') f90file//':'//subr// &
            ': Horizontal resolution on file does not match'// &
            ' specified resolution'
       write(6,'(a,2i5)') '  Specified:',IRES,JRES
       write(6,'(a,2i5)') '  On file:  ',nLon,nLat
       stop
    end if
    if (nTime.eq.1) then
       write(6,'(a)') f90file//':'//subr//': One dataset'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': not 1 or 12 datasets - do not know what to do'
       stop
    end if
       

    !// Emission variable, expected units: kg/m2/s for each month
    call get_netcdf_var_2d( infile, invar, R8XY(:,:,1), nlon,nlat )

    !// Check for nans
    do J = 1, nlat
       do I = 1, nlon
          if (R8XY(I,J,1) .ne. R8XY(I,J,1)) R8XY(I,J,1) = 0._r8
       end do
    end do

    !// Number of datasets
    NSETS = nTime

    !// --------------------------------------------------------------------
  end subroutine nc_getdata_sfc_2d
  !// ----------------------------------------------------------------------





  !// ----------------------------------------------------------------------
  subroutine nc_getch4emis_game(INFILE,INVAR,INYEAR,R8XY,IRES,JRES,MRES,NSETS)
    !// --------------------------------------------------------------------
    !// Reads CH4 surface emissions for GAME for 12 months, from netCDF file.
    !// Contains 1984 to 2009, use INMON to select.
    !//
    !//
    !// Ole Amund Sovde, February 2013
    !// --------------------------------------------------------------------
    use cmn_met, only: MYEAR
    use netcdf
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
 
    !// Input
    character(len=*), intent(in)  :: infile   !// Name of netCDFfile
    character(len=*), intent(in)  :: invar    !// Name of input variable
    integer,intent(in) :: INYEAR   !// Year to fetch
    integer,intent(in) :: IRES     !// Number of longitudes
    integer,intent(in) :: JRES     !// Number of latitudes
    integer,intent(in) :: MRES     !// Number of sets (months) per category

    !// Output
    real(r8), intent(out) :: R8XY(IRES,JRES,MRES)  !// 3D field
    integer,intent(out)   :: NSETS

    !// Locals
    integer            :: nLon, nLat, nLev, nTime, sTime, M, getY
    character(len=25)  :: units !CDH caution using fixed length
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    real(r8), allocatable, dimension(:,:,:) :: allindata
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nc_getch4emis_game'
    !// --------------------------------------------------------------------

    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'longitude',  inLon  )
    call get_netcdf_var_1d( infile, 'latitude',  inLat  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    !// Which year to fetch
    if (INYEAR .eq. 9999) then
       !// Use meteorological year
       getY = MYEAR
    else
       getY = INYEAR
    end if


    !// Select months to get: Find first entry
    sTime = -1
    do M = 1, nTime
       if (inTime(M) .ge. getY) then
          sTime = M
          exit
       end if
    end do
    if (sTime .le. 0) then
       write(6,'(a,i4,a)') f90file//':'//subr// &
            ': Cannot find year ',getY,' in dataset:'//trim(infile)
       stop
    end if

    if (nLon.ne.IRES .or. nLat.ne.JRES) then
       write(6,'(a)') f90file//':'//subr// &
            ': Horizontal resolution on file does not match'// &
            ' specified resolution'
       write(6,'(a,2i5)') '  Specified:',IRES,JRES
       write(6,'(a,2i5)') '  On file:  ',nLon,nLat
       stop
    end if

    write(6,'(a,i5,f12.3)') '  Month index (sTime) and time (inTime(sTime)):',sTime,inTime(sTime)

    !// Test time sTime:sTime+11
    if (sTime+11 .gt. ntime) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problem reading emission dataset'
       write(6,'(a,2i7)') '  stime/ntime:',sTime,nTime
       stop
    end if


    !// Emission variable, expected units: kg/m2/s for each month
    allocate(allindata(nlon,nlat,ntime))
    call get_netcdf_var_3d( infile, invar, allindata, nlon,nlat,ntime )


    !// Fetch R8XY
    do M = 1, 12
       R8XY(:,:,M) = allindata(:,:,sTime+M-1)
    end do



    !// Number of datasets
    NSETS = 12

    !// Deallocate all local variables
    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)
    if ( allocated(inTime) ) deallocate(inTime)
    if ( allocated(allindata) ) deallocate(allindata)

    !// --------------------------------------------------------------------
  end subroutine nc_getch4emis_game
  !// ----------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine gfed4_rd(JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This reads monthly emissions of dry matter (DM) from GFEDv4,
    !// for the six partitions:
    !// SAVA: Savanna
    !// BORF: Boreal forrest
    !// TEMF: Temperature forest
    !// DEFO: Deforestation
    !// PEAT: Peat
    !// AGRI: Agricultural wastet burning
    !//
    !// Input resolution is 0.25x0.25 degrees, which is first interpolated
    !// to 0.5x0.5 degrees, then distributed according
    !// to RETRO vertical distribution (also 0.5x0.5 degrees), before
    !// interpolated to CTM3 resolution.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// The routine is hardcoded to read separate emission array (EMIS_FIR).
    !// This is necessary due to all scalings for NMVHCs, and because the
    !// emission fields need to be independent of emission arrays
    !// E2DS and E3DS to allow for more frequent updating.
    !//
    !// Variables/parameters defined in oc_globalvariables.f:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Ole Amund Sovde, February 2012, updated June 2013
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_SCALE, FF_PARTITIONS, EPAR_NPARTS
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,M,N,N_EMIS,II,JJ,MP, io_err

    !// Dimensions
    integer, parameter :: IRES=360 ,JRES=180, IHLF=720, JHLF=360, LHLF=13
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: HXBEDGE(IHLF+1), HYBEDGE(JHLF+1),HXYBOX(JHLF)
    real(r8)   :: QXYBOX(JQRT)

    real(r8) :: BBH_FIR(IHLF,JHLF,EPAR_FIR_LM), HLFDUM(IHLF,JHLF)
    real(r8) :: R8Q(IQRT,JQRT)

    real(r8) :: RDUM(IHLF,JHLF)
    character(len=20) :: COMP
    character(len=2) :: CMONTH
    character(len=4) :: CYEAR
    !// Partitions: Agriculture, deforestation, forest, peat, savanna, woodland
    integer, parameter :: NPARTS = 6
    integer :: partitions(NPARTS)
    real(r8)  :: DM_PARTS(IHLF,JHLF,NPARTS)
    character(len=4) :: CPARTS(NPARTS)

    !// Interpolated array
    real(r8)  :: RXY8(IPAR,JPAR,EPAR_FIR_LM)

    !// File variables
    logical :: file_status
    integer :: efnr, nLon, nLat
    character(len=120) :: INFILE
    integer :: NTRNR
    real(r8) :: EDATA(IHLF,JHLF,EPAR_FIR),tscale, scalefac

    !// GFED4 emission factors
    character(len=16) :: VAR
    real(r8) :: f1,f2,f3,f4,f5,f6
    integer, parameter :: NG4 = 41
    integer :: GFED4_N
    real(r8) :: GFED4_EF(NPARTS,NG4)
    character(len=16) :: GFED4_NM(NG4)

    real(r8), allocatable, dimension(:) :: QXBEDGE, QYBEDGE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed4_rd'
    !// ------------------------------------------------------------------

    write(6,'(a)') subr//': Biomass burning (GFEDv4)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// DM emissions from each partition, to be read from file
    DM_PARTS(:,:,:) = 0._r8

    !// ------------------------------------------------------------------

    !// Select partitions to use (default is all)
    CPARTS = (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)

    !// Get an unused file number
    file_status = .true.
    efnr = 7
    do while (file_status)
       efnr = efnr + 1
       inquire(efnr,opened=file_status)
    end do
    
    !// Get emission factors for all GFED4 species
    infile = trim(FF_PATH)//'GFED4_Emission_Factors.txt'
    open(efnr,file=infile,form='formatted',action='read')
    do I = 1, 17
       read(efnr,*) ! Reading header
    end do
    do I = 1, NG4
       !// Read emission factors
       read(efnr,*) VAR,f1,f2,f3,f4,f5,f6
       GFED4_EF(1:6,I) = (/ f1,f2,f3,f4,f5,f6 /)
       GFED4_NM(I) = VAR
    end do
    close(efnr)

    write(6,'(a)') 'GFED4 emission factors'
    write(6,'(a16,6a9)') 'Component',CPARTS
    do I = 1, NG4
       write(6,'(a16,6f9.3)')GFED4_NM(I),GFED4_EF(1:6,I)
    end do


    !// Read DM from partitions
    !// ------------------------------------------------------------------
    !// Month to read
    write(CMONTH(1:2),'(i2.2)') JMONTH
    !// Year to read
    if (FF_YEAR .ne. 9999) then
       write(CYEAR,'(i4.4)') FF_YEAR
    else
       write(CYEAR,'(i4.4)') MYEAR !// Use meteorological year
    end if

    !// Filename
    infile = trim(FF_PATH)//'GFED4.1s_'//CYEAR//'_'//CMONTH//'.nc'

    !// Inquire file
    inquire(file=infile,exist=file_status)

    if (.not. file_status) then
       write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
            trim(infile)
       !// Fall-back to 1997/2016
       if (MYEAR .gt. 2016) then
          infile = trim(FF_PATH)//'GFED4.1s_2016_'//CMONTH//'.nc'
       else if (MYEAR .lt. 1997) then
          infile = trim(FF_PATH)//'GFED4.1s_1997_'//CMONTH//'.nc'
       else
          stop
       end if
       !// In case of fall-back
       write(6,'(a)') 'Will try to use file: '//trim(infile)
       inquire(file=infile,exist=file_status)
       if (.not. file_status) then
          write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
               trim(infile)
          stop
       end if
    end if


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lone', QXBEDGE )
    call get_netcdf_var_1d( infile, 'late', QYBEDGE )

    nLon  = SIZE( QXBEDGE ) - 1
    nLat  = SIZE( QYBEDGE ) - 1

    if (nLon .ne. IQRT .or. nLat.ne.JQRT) then
       write(6,'(a,2i5)') subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do


    !// Set up halfdegrees (Grid type is 1, start at 0E, 90S)
    !// (Actually, GFED starts at 180W, 90N, but it is swicthed after read-in.)
    call get_xyedges(IHLF,JHLF,HXBEDGE,HYBEDGE,1)

    !// HLF Grid box areas
    do J = 1, JHLF
       HXYBOX(J) =  A0*A0 * CPI180*(HXBEDGE(2) - HXBEDGE(1)) &
            * (sin(CPI180*HYBEDGE(J+1)) - sin(CPI180*HYBEDGE(J)))
    end do



    !// Read DM for six categories; save data on 0.5x0.5 degree
    !// For each species, add all partitions together in 0.5x0.5.
    !// Finally convert to model resolution.
    do N = 1, NPARTS

       !// Emission variable, expected units: kg/m2/s for each month
       !// GFED4 on netcdf has been converted to start at 0E,90S
       call get_netcdf_var_2d( infile, 'DM_'//CPARTS(N),R8Q(:,:), nlon,nlat )

       !// multiply by area for interpolation
       do J = 1, JQRT
          R8Q(:,J) = R8Q(:,J) * QXYBOX(J)
       end do

       !// From 0.25 to 0.5 degree
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, HLFDUM, &
                 HXBEDGE,HYBEDGE,IHLF,JHLF,1)

       !// divide by area after interpolation
       do J = 1, JHLF
          HLFDUM(:,J) = HLFDUM(:,J) / HXYBOX(J)
       end do

       !// Save DM for each partition
       DM_PARTS(:,:,N) = HLFDUM(:,:)

    end do !// do N = 1, NPARTS

    !// Deallocate allocated variables
    if ( allocated(QXBEDGE) ) deallocate(QXBEDGE)
    if ( allocated(QYBEDGE) ) deallocate(QYBEDGE)


    !// Get vertical distribution
    !// ------------------------------------------------------------------
    INFILE = trim(FF_PATH)//'bb_altitudes_aggregated.0.5x0.5.nc'
    write(6,'(a)') ' * Reading BBH '//trim(INFILE)

    !// Read the BBH data
    call read_retrobbh(INFILE,IHLF,JHLF,EPAR_FIR_LM,BBH_FIR,LPARW,LMMAP)

    !// Check BB
    do J = 1, JHLF
       do I = 1, IHLF
          if (BBH_FIR(I,J,1) .eq.0._r8 .and. sum(BBH_FIR(I,J,:)).gt.0._r8) then
             write(6,'(a,2i5)') f90file//':'//subr// &
                  ': BBH STRANGE'
             do L = 1, EPAR_FIR_LM
                write(6,'(i5,es12.3)') L, BBH_FIR(I,J,L)
             end do
             stop
          end if
       end do
    end do



    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    !// Initialize 3D HLF array (which will be interpolated later)
    EDATA(:,:,:) = 0._r8

    do M = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(M) .le. 0) cycle

       !// Find corresponding GFED4 name. Some use same as CTM3, others
       !// are different and must be overrided below.
       COMP = trim(FF_BNAMES(M))
       !// Additional scaling factor
       scalefac = FF_SCALE(M)

       !// Select partitions to use (default is all)
       !// 'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'
       partitions(:) = FF_PARTITIONS(:,M)


       !// Find the GFED4 component number, if available.
       GFED4_N = -1
       do I = 1, NG4
          if (trim(GFED4_NM(I)) .eq. trim(COMP)) then
             GFED4_N = I
             exit
          end if
       end do

       !// If not available, skip component
       if (GFED4_N .le. 0) cycle

       !// Sum up partitions for this species
       RDUM(:,:) = 0._r8
       do N = 1, NPARTS

          if (partitions(N) .eq. 0) cycle
          !// Emission factors are g/kgDM, so we convert to kg, and also
          !// apply the scalefac.
          tscale =  GFED4_EF(N,GFED4_N) * 1.e-3_r8 * scalefac

          RDUM(:,:) = RDUM(:,:) + DM_PARTS(:,:,N) * tscale

       end do

       !// Convert to kg/s
       do J = 1, JHLF
          EDATA(:,J,M) = RDUM(:,J) * HXYBOX(J)
       end do

    end do !// do N_EMIS = 1, NEFIR



    !// Loop over components
    do N_EMIS = 1, NEFIR
       !// If component is not included, go to next N_EMIS
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       NTRNR = ECOMP_FIR(N_EMIS)


       !// Distribute with height
!$omp parallel private (L,I,J,HLFDUM) &
!$omp          shared (EDATA,BBH_FIR,RXY8,LMMAP,HXBEDGE,HYBEDGE, &
!$omp                  XDEDG,YDEDG,N_EMIS) &
!$omp          default(NONE)
!$omp do
       do L = 1, EPAR_FIR_LM
          HLFDUM(:,:) = EDATA(:,:,N_EMIS) * BBH_FIR(:,:,LMMAP(L))
          if (maxval(HLFDUM) .gt. 0._r8) then
             if (L.eq.1) then
                !// Check problems with BBH on half degree vs emissions on
                !// one degree
                do J = 1, JHLF
                  do I = 1, IHLF
                     !// Put everyting in layer 1
                     if (BBH_FIR(I,J,1) .eq. 0._r8 .and. &
                          HLFDUM(I,J) .gt. 0._r8) then
                        HLFDUM(I,J) = EDATA(I,J,N_EMIS)
                     end if
                  end do
               end do
            end if

            !// Interpolate
            call E_GRID(HLFDUM,HXBEDGE,HYBEDGE,IHLF,JHLF, RXY8(:,:,L), &
                 XDEDG,YDEDG,IPAR,JPAR,1)
         else
            RXY8(:,:,L) = 0._r8
         end if
      end do !// do L = 1, EPAR_FIR_LM
!$omp end do
!$omp end parallel


      !// Put into EMIS_FIR
!$omp parallel private (MP,I,J,L,II,JJ) &
!$omp         shared (RXY8,EMIS_FIR,MPBLKJB,MPBLKJE,MPBLKIB,MPBLKIE,N_EMIS) &
!$omp         default(NONE)
!$omp do
      do MP = 1, MPBLK
        !// Loop over latitude (J is global, JJ is block)
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1

          !// Loop over longitude (I is global, II is block)
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1

            do L = 1, EPAR_FIR_LM
               EMIS_FIR(L,N_EMIS,II,JJ,MP) = RXY8(I,J,L)
            end do
          end do
        end do
      end do

!$omp end do
!$omp end parallel
      tscale = 0._r8
      do MP = 1, MPBLK
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = 1, EPAR_FIR_LM
               tscale = tscale + EMIS_FIR(L,N_EMIS,II,JJ,MP)
            end do
          end do
        end do
      end do
      write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ',FF_CNAMES(N_emis),NTRNR,&
           tscale*2628000.e-9_r8,' Tg/mon'

    end do !// do do N_EMIS = 1, NEFIR

    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'

    !// ------------------------------------------------------------------
  end subroutine gfed4_rd
  !// ------------------------------------------------------------------





  !// ------------------------------------------------------------------
  subroutine gfed4_init()
    !// ------------------------------------------------------------------
    !//
    !// Amund Sovde Haslerud, December 2017
    !// ------------------------------------------------------------------
    use cmn_oslo, only: FF_PATH, GFED4_NC, GFED4_NP, GFED4_EF, GFED4_NM
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input

    !// Locals
    integer :: I

    !// Partitions
    character(len=4), dimension(GFED4_NP) :: CPARTS

    !// File variables
    logical :: file_status
    integer :: efnr, io_err
    character(len=120) :: INFILE

    !// GFED4 emission factors
    character(len=16) :: VAR
    real(r8) :: f1,f2,f3,f4,f5,f6
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed4_init'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Initialising GFED4 biomass burning'

    !// Get an unused file number
    file_status = .true.
    efnr = 7
    do while (file_status)
       efnr = efnr + 1
       inquire(efnr,opened=file_status)
    end do

    CPARTS = (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)
    
    !// Get emission factors for all GFED4 species
    infile = trim(FF_PATH)//'GFED4_Emission_Factors.txt'
    open(efnr,file=infile,form='formatted',action='read',iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': Error opening '//trim(infile)
       stop
    end if
    do I = 1, 17
       read(efnr,*) ! Reading header
    end do
    do I = 1, GFED4_NC
       !// Read emission factors
       read(efnr,*) VAR,f1,f2,f3,f4,f5,f6
       GFED4_EF(1:6,I) = (/ f1,f2,f3,f4,f5,f6 /)
       GFED4_NM(I) = VAR
    end do
    close(efnr)

    write(6,'(a)') 'GFED4 emission factors'
    write(6,'(a16,6a9)') 'Component',CPARTS
    do I = 1, GFED4_NC
       write(6,'(a16,6f9.3)') GFED4_NM(I),GFED4_EF(1:6,I)
    end do

    write(6,'(a)') f90file//':'//subr//': initialised.'

    !// ------------------------------------------------------------------
  end subroutine gfed4_init
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine gfed4_rd_novert(JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This reads monthly emissions of dry matter (DM) from GFEDv4,
    !// for the six partitions:
    !// SAVA: Savanna
    !// BORF: Boreal forrest
    !// TEMF: Temperature forest
    !// DEFO: Deforestation
    !// PEAT: Peat
    !// AGRI: Agricultural wastet burning
    !//
    !// Input resolution is 0.25x0.25 degrees, interpolated directly
    !// into model resolution, with no vertical distribution.
    !// This is supposed to be combined with distributing evenly within
    !// PBL.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// The routine is hardcoded to read separate emission array (EMIS_FIR).
    !// This is necessary due to all scalings for NMVHCs, and because the
    !// emission fields need to be independent of emission arrays
    !// E2DS and E3DS to allow for more frequent updating.
    !//
    !// Variables/parameters defined in oc_globalvariables.f:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Ole Amund Sovde, February 2012, updated June 2013
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_SCALE, FF_PARTITIONS, EPAR_NPARTS
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,M,N,N_EMIS,II,JJ,MP, io_err

    !// Dimensions
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: QXYBOX(JQRT)

    real(r8) :: R8Q(IQRT,JQRT)
    real(r8) :: RDUM(IQRT,JQRT)

    character(len=20) :: COMP
    character(len=2) :: CMONTH
    character(len=4) :: CYEAR

    !// Partitions: Agriculture, deforestation, forest, peat, savanna, woodland
    integer, parameter :: NPARTS = 6
    integer :: partitions(NPARTS)
    real(r8)  :: DM_PARTS(IQRT,JQRT,NPARTS)
    character(len=4) :: CPARTS(NPARTS)

    !// Interpolated array
    real(r8)  :: RXY8(IPAR,JPAR)

    !// File variables
    logical :: file_status
    integer :: efnr, nLon, nLat
    character(len=120) :: INFILE
    real(r8) :: tscale, scalefac

    !// GFED4 emission factors
    character(len=16) :: VAR
    real(r8) :: f1,f2,f3,f4,f5,f6
    integer, parameter :: NG4 = 41
    integer :: GFED4_N
    real(r8) :: GFED4_EF(NPARTS,NG4)
    character(len=16) :: GFED4_NM(NG4)

    real(r8), allocatable, dimension(:) :: QXBEDGE, QYBEDGE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed4_rd_novert'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Biomass burning (GFEDv4)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// DM emissions from each partition, to be read from file
    DM_PARTS(:,:,:) = 0._r8

    !// ------------------------------------------------------------------

    !// Select partitions to use (default is all)
    CPARTS = (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)

    !// Get an unused file number
    file_status = .true.
    efnr = 7
    do while (file_status)
       efnr = efnr + 1
       inquire(efnr,opened=file_status)
    end do
    
    !// Get emission factors for all GFED4 species
    infile = trim(FF_PATH)//'GFED4_Emission_Factors.txt'
    open(efnr,file=infile,form='formatted',action='read')
    do I = 1, 17
       read(efnr,*) ! Reading header
    end do
    do I = 1, NG4
       !// Read emission factors
       read(efnr,*) VAR,f1,f2,f3,f4,f5,f6
       GFED4_EF(1:6,I) = (/ f1,f2,f3,f4,f5,f6 /)
       GFED4_NM(I) = VAR
    end do
    close(efnr)

    write(6,'(a)') 'GFED4 emission factors'
    write(6,'(a16,6a9)') 'Component',CPARTS
    do I = 1, NG4
       write(6,'(a16,6f9.3)')GFED4_NM(I),GFED4_EF(1:6,I)
    end do


    !// Read DM from partitions
    !// ------------------------------------------------------------------
    !// Month to read
    write(CMONTH(1:2),'(i2.2)') JMONTH
    !// Year to read
    if (FF_YEAR .ne. 9999) then
       write(CYEAR,'(i4.4)') FF_YEAR
    else
       write(CYEAR,'(i4.4)') MYEAR !// Use meteorological year
    end if

    !// Filename
    infile = trim(FF_PATH)//'GFED4.1s_'//CYEAR//'_'//CMONTH//'.nc'

    !// Inquire file
    inquire(file=infile,exist=file_status)

    if (.not. file_status) then
       write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
            trim(infile)
       !// Fall-back to 1997/2016
       if (MYEAR .gt. 2016) then
          infile = trim(FF_PATH)//'GFED4.1s_2016_'//CMONTH//'.nc'
       else if (MYEAR .lt. 1997) then
          infile = trim(FF_PATH)//'GFED4.1s_1997_'//CMONTH//'.nc'
       else
          stop
       end if
       !// In case of fall-back
       write(6,'(a)') 'Will try to use file: '//trim(infile)
       inquire(file=infile,exist=file_status)
       if (.not. file_status) then
          write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
               trim(infile)
          stop
       end if
    end if


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lone', QXBEDGE )
    call get_netcdf_var_1d( infile, 'late', QYBEDGE )

    nLon  = SIZE( QXBEDGE ) - 1
    nLat  = SIZE( QYBEDGE ) - 1

    if (nLon .ne. IQRT .or. nLat.ne.JQRT) then
       write(6,'(a,2i5)') f90file//':'//subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do



    !// Read DM for six categories; save data on 0.5x0.5 degree
    !// For each species, add all partitions together in 0.5x0.5.
    !// Finally convert to model resolution.
    do N = 1, NPARTS

       !// Emission variable, expected units: kg/m2/s for each month
       !// GFED4 on netcdf has been converted to start at 0E,90S
       call get_netcdf_var_2d( infile, 'DM_'//CPARTS(N),R8Q(:,:), nlon,nlat )

       !// Save DM for each partition
       DM_PARTS(:,:,N) = R8Q(:,:)

    end do !// do N = 1, NPARTS



    !// Get vertical distribution
    !// ------------------------------------------------------------------
    !// None is used


    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    do N_EMIS = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       !// Find corresponding GFED4 name. Some use same as CTM3, others
       !// are different and must be overrided below.
       COMP = trim(FF_BNAMES(N_EMIS))
       !// Additional scaling factor
       scalefac = FF_SCALE(N_EMIS)

       !// Select partitions to use (default is all)
       !// 'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'
       partitions(:) = FF_PARTITIONS(:,N_EMIS)


       !// Find the GFED4 component number, if available.
       GFED4_N = -1
       do I = 1, NG4
          if (trim(GFED4_NM(I)) .eq. trim(COMP)) then
             GFED4_N = I
             exit
          end if
       end do

       !// If not available, skip component
       if (GFED4_N .le. 0) cycle

       !// Sum up partitions for this species
       R8Q(:,:) = 0._r8
       do N = 1, NPARTS

          if (partitions(N) .eq. 0) cycle
          !// Emission factors are g/kgDM, so we convert to kg, and also
          !// apply the scalefac.
          tscale =  GFED4_EF(N,GFED4_N) * 1.e-3_r8 * scalefac

          R8Q(:,:) = R8Q(:,:) + DM_PARTS(:,:,N) * tscale

       end do

       !// Convert to kg/s
       do J = 1, JQRT
          R8Q(:,J) = R8Q(:,J) * QXYBOX(J)
       end do


       !// Interpolate
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, RXY8(:,:), &
                 XDEDG,YDEDG,IPAR,JPAR,1)

       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1

                EMIS_FIR(1,N_EMIS,II,JJ,MP) = RXY8(I,J)
             end do
          end do
       end do

       tscale = sum(RXY8)
       write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ', &
            FF_CNAMES(N_emis),ECOMP_FIR(N_EMIS),&
            tscale*2628000.e-9_r8,' Tg/mon'

    end do !// do do N_EMIS = 1, NEFIR

    !// Deallocate allocated variables
    if ( allocated(QXBEDGE) ) deallocate(QXBEDGE)
    if ( allocated(QYBEDGE) ) deallocate(QYBEDGE)


    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'

    !// ------------------------------------------------------------------
  end subroutine gfed4_rd_novert
  !// ------------------------------------------------------------------




  !// ------------------------------------------------------------------
  subroutine gfed4_rd_novert_daily(JDATE,JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This routine reads monthly emissions of dry matter (DM) from GFEDv4,
    !// along with daily fractions of the monthly totals.
    !// To be accurate, the emissions are read every day - although they
    !// could have been read once per month and kept in original resolution.
    !//
    !// Six partitions:
    !// SAVA: Savanna
    !// BORF: Boreal forrest
    !// TEMF: Temperature forest
    !// DEFO: Deforestation
    !// PEAT: Peat
    !// AGRI: Agricultural wastet burning
    !//
    !// Input resolution is 0.25x0.25 degrees, interpolated directly
    !// into model resolution, with no vertical distribution.
    !// This is supposed to be combined with distributing evenly within
    !// PBL.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// The routine is hardcoded to read separate emission array (EMIS_FIR).
    !// This is necessary due to all scalings for NMVHCs, and because the
    !// emission fields need to be independent of emission arrays
    !// E2DS and E3DS to allow for more frequent updating.
    !//
    !// Variables/parameters defined in oc_globalvariables.f:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Ole Amund Sovde, February 2012, updated June 2013
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_SCALE, FF_PARTITIONS, EPAR_NPARTS, &
         GFED4_EF, GFED4_NM, GFED4_NC, GFED4_NP, DINM
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JDATE, JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,M,N,N_EMIS,II,JJ,MP

    !// Dimensions
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: QXYBOX(JQRT)

    real(r8), dimension(IQRT,JQRT) :: R8Q, DFRAC
    real(r8) :: RDUM(IQRT,JQRT)

    character(len=20) :: COMP
    character(len=2) :: CMONTH, CDATE
    character(len=4) :: CYEAR

    !// Partitions: Agriculture, deforestation, forest, peat, savanna, woodland
    integer, dimension(GFED4_NP) :: partitions
    real(r8), dimension(IQRT,JQRT,GFED4_NP)  :: DM_PARTS
    character(len=4), dimension(GFED4_NP) :: CPARTS

    !// Interpolated array
    real(r8), dimension(IPAR,JPAR) :: RXY8

    !// File variables
    logical :: file_status
    integer :: efnr, nLon, nLat
    character(len=120) :: INFILE
    real(r8), dimension(IQRT,JQRT,EPAR_FIR) :: EDATA
    real(r8), dimension(EPAR_FIR) :: ESUM
    real(r8) :: tscale, scalefac

    !// GFED4 emission factors
    character(len=16) :: VAR
    real(r8) :: f1,f2,f3,f4,f5,f6
    integer :: GFED4_N

    real(r8), allocatable, dimension(:) :: QXBEDGE, QYBEDGE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed4_rd_novert_daily'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Biomass burning (GFEDv4)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// DM emissions from each partition, to be read from file
    DM_PARTS(:,:,:) = 0._r8

    !// Select partitions to use (default is all)
    CPARTS = (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)

    !// Read DM from partitions
    !// ------------------------------------------------------------------
    !// Date and month to read
    write(CDATE(1:2),'(i2.2)') JDATE
    write(CMONTH(1:2),'(i2.2)') JMONTH
    !// Year to read
    if (FF_YEAR .ne. 9999) then
       write(CYEAR,'(i4.4)') FF_YEAR
    else
       write(CYEAR,'(i4.4)') MYEAR !// Use meteorological year
    end if

    !// Filename
    infile = trim(FF_PATH)//'GFED4.1s_'//CYEAR//'_'//CMONTH//'.nc'

    !// Inquire file
    inquire(file=infile,exist=file_status)

    if (.not. file_status) then
       write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
            trim(infile)
       !// Fall-back to 2003/2016
       if (MYEAR .gt. 2016) then
          infile = trim(FF_PATH)//'GFED4.1s_2016_'//CMONTH//'.nc'
       else if (MYEAR .lt. 2003) then
          infile = trim(FF_PATH)//'GFED4.1s_2003_'//CMONTH//'.nc'
       else
          stop
       end if
       !// In case of fall-back
       write(6,'(a)') 'Will try to use file: '//trim(infile)
       inquire(file=infile,exist=file_status)
       if (.not. file_status) then
          write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
               trim(infile)
          stop
       end if
    end if


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lone', QXBEDGE )
    call get_netcdf_var_1d( infile, 'late', QYBEDGE )

    nLon  = SIZE( QXBEDGE ) - 1
    nLat  = SIZE( QYBEDGE ) - 1

    if (nLon .ne. IQRT .or. nLat.ne.JQRT) then
       write(6,'(a,2i5)') f90file//':'//subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do



    !// Read DM for six categories; save data on 0.5x0.5 degree
    !// For each species, add all partitions together in 0.5x0.5.
    !// Finally convert to model resolution.
    do N = 1, GFED4_NP

       !// Emission variable, expected units: kg/m2/s for each month
       !// GFED4 on netcdf has been converted to start at 0E,90S
       call get_netcdf_var_2d( infile, 'DM_'//CPARTS(N),R8Q(:,:), nlon,nlat )

       !// Save DM for each partition
       DM_PARTS(:,:,N) = R8Q(:,:)

    end do !// do N = 1, GFED4_NP



    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    !// Initialize 3D array (which will be interpolated later)
    EDATA(:,:,:) = 0._r8

    do M = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(M) .le. 0) cycle

       !// Find corresponding GFED4 name. Some use same as CTM3, others
       !// are different and must be overrided below.
       COMP = trim(FF_BNAMES(M))
       !// Additional scaling factor
       scalefac = FF_SCALE(M)

       !// Select partitions to use (listed in emissions list)
       !// 'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'
       partitions(:) = FF_PARTITIONS(:,M)


       !// Find the GFED4 component number, if available.
       GFED4_N = -1
       do I = 1, GFED4_NC
          if (trim(GFED4_NM(I)) .eq. trim(COMP)) then
             GFED4_N = I
             exit
          end if
       end do

       !// If not available, skip component
       if (GFED4_N .le. 0) cycle

       !// Sum up partitions for this species
       RDUM(:,:) = 0._r8
       do N = 1, GFED4_NP

          if (partitions(N) .eq. 0._r8) cycle
          !// Emission factors are g/kgDM, so we convert to kg, and also
          !// apply the scalefac.
          tscale =  GFED4_EF(N,GFED4_N) * 1.e-3_r8 * scalefac

          RDUM(:,:) = RDUM(:,:) + DM_PARTS(:,:,N) * tscale

       end do

       !// Convert to kg/s
       do J = 1, JQRT
          EDATA(:,J,M) = RDUM(:,J) * QXYBOX(J)
       end do

    end do !// do N_EMIS = 1, NEFIR


    !// Include daily fraction of monthly totals
    !// ------------------------------------------------------------------
    infile = trim(FF_PATH)//'daily/'//CYEAR// &
         '/fraction_emissions_'//CYEAR//CMONTH//CDATE//'_nc4.nc'
    write(6,'(a)') f90file//':'//subr//': Reading '//trim(INFILE)
    !// Emission variable, expected units: kg/m2/s for each month
    call get_netcdf_var_2d(infile, 'Fraction_of_Emissions', RDUM, IQRT,JQRT)

    !// IMPORTANT: Must scale by days in month, since emissions
    !//            are given in units kg/s.
    !// Skip possible negatives.
    RDUM(:,:) = max(0._r8, RDUM(:,:)) * real(DINM(JMONTH),r8)
    !// Switch from -180:180 to 0:360, as GFED4 monthly fields
    do J = 1, JQRT
       DFRAC(1:IQRT/2,J) = RDUM(IQRT/2+1:IQRT,J)
       DFRAC(IQRT/2+1:IQRT,J) = RDUM(1:IQRT/2,J)
    end do


    !// Regrid - Loop over components
    !// ------------------------------------------------------------------
    ESUM(:) = 0._r8
!$omp parallel private (MP,I,J,II,JJ,RXY8,R8Q) &
!$omp          shared (N_EMIS,NEFIR,EDATA, QXBEDGE,QYBEDGE, XDEDG,YDEDG, &
!$omp                  ECOMP_FIR,ESUM,EMIS_FIR, DFRAC, &
!$omp                  MPBLKJB,MPBLKJE,MPBLKIB,MPBLKIE) &
!$omp          default(NONE)
!$omp do
    do N_EMIS = 1, NEFIR
       !// If component is not included, go to next N_EMIS
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       R8Q(:,:) = EDATA(:,:,N_EMIS) * DFRAC(:,:)

       !// Interpolate
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, RXY8(:,:), &
                 XDEDG,YDEDG,IPAR,JPAR,1)

       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1

                EMIS_FIR(1,N_EMIS,II,JJ,MP) = RXY8(I,J)
             end do
          end do
       end do

       ESUM(N_EMIS) = sum(RXY8)

    end do !// do do N_EMIS = 1, NEFIR
!$omp end do
!$omp end parallel

    !// Write out estimated totals
    do N_EMIS = 1, NEFIR
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle
       write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ',&
            FF_CNAMES(N_emis),ECOMP_FIR(N_EMIS),&
            ESUM(N_EMIS) * 86400.e-9_r8,' Tg/day'
    end do !// do do N_EMIS = 1, NEFIR


    !// Deallocate allocated variables
    if ( allocated(QXBEDGE) ) deallocate(QXBEDGE)
    if ( allocated(QYBEDGE) ) deallocate(QYBEDGE)


    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'

    !// ------------------------------------------------------------------
  end subroutine gfed4_rd_novert_daily
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine read_ceds_cicero(INFILE, INVAR, INYEAR, R8XY, &
       IRES, JRES, MRES, NSETS, ESCENYEAR)
    !// ----------------------------------------------------------------
    !// Read CEDS emission files generated at CICERO. I.e. not the
    !// original CEDS files - those files were 4D where all sectors were
    !// included, which was tricky to include in CTM3.
    !//
    !// Amund Sovde Haslerud, March 2017
    !// ----------------------------------------------------------------
    use netcdf
    use ncutils, only: get_netcdf_var_1d, readnc_sub3d_from3d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
 
    !// Input
    character(len=*), intent(in)  :: infile   !// Name of netCDFfile
    character(len=*), intent(in)  :: invar    !// Name of input variable
    integer, intent(in)         :: INYEAR   !// Year to fetch
    integer, intent(in)         :: ESCENYEAR   !// Year to fetch !MTL, not use INYEAR anymore
    integer, intent(in)         :: IRES     !// Number of longitudes
    integer, intent(in)         :: JRES     !// Number of latitudes
    integer, intent(in)         :: MRES     !// Number of sets (months) per category
    !// Output
    real(r8), intent(out)       :: R8XY(IRES,JRES,MRES)  !// 3D field
    integer, intent(out)        :: NSETS    !// Number of sets (months) per category

    !// Locals
    integer            :: nLon, nLat, nLev, nTime, sTime, M, getY
    character(len=25)  :: units !CDH caution using fixed length
    character(len=10)  :: svar
    integer, dimension(3) :: srt_entry, cnt_entry
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    real(r8), allocatable, dimension(:,:,:) :: sub_indata
    !// ------------------------------------------------------------------
    character(len=*), parameter :: subr = 'read_ceds_cicero'
    !// ------------------------------------------------------------------

    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    !// Will fetch emission data for year ESCENYEAR
    getY = ESCENYEAR

    !// Select months to get: Find first entry of year getY;
    !// inTime is given as days since 1750, and starts at mid-month
    !// (15 Jan for each year).
    sTime = -1
    do M = 1, nTime
       if (inTime(M) .ge. (getY - 1750)*365 + 15) then
          sTime = M
          exit
       end if
    end do
    if (sTime .le. 0) then
       write(6,'(a,i4,a)') f90file//':'//subr// &
            ': Cannot find year ',getY,' in file: '//trim(infile)
       stop
    end if

    if (nLon.ne.IRES .or. nLat.ne.JRES) then
       write(6,'(a)') f90file//':'//subr// &
            ': Horizontal resolution mismatch'
       write(6,'(a)')       '  File: '//trim(infile)
       write(6,'(a,2(i5))') '  Specified:',IRES,JRES
       write(6,'(a,2(i5))') '  On file:  ',nLon,nLat
       stop
    end if

    write(6,'(a,i5,f12.3)') '  Month index (sTime) and time (inTime(sTime)):',sTime,inTime(sTime)

    !// Test time sTime:sTime+11
    if (sTime+11 .gt. ntime) then
       write(6,'(a,2i7)') f90file//':'//subr// &
            ': stime/ntime:',sTime,nTime
       stop
    end if

    !// Will fetch 12 months
    NSETS = 12

    !// Emission variable, expected units: kg/m2/s for each month
    !// Will retrieve NSETS datasets (months)
    allocate( sub_indata(nLon,nLat,NSETS) )

    !// Start entry of dataset
    srt_entry = (/ 1, 1, stime /)
    !// Get NSETS datasets
    cnt_entry = (/ nLon, nLat, NSETS /)

    call readnc_sub3d_from3d(INFILE, 'lon', nlon, 'lat', nlat,'time',nTime,&
         invar, srt_entry, cnt_entry, sub_indata)

    R8XY(:,:,1:NSETS) = sub_indata(:,:,:)

    !// Deallocate all local variables
    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)
    if ( allocated(inTime) ) deallocate(inTime)
    if ( allocated(sub_indata) ) deallocate(sub_indata)




    !// ----------------------------------------------------------------
  end subroutine read_ceds_cicero
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine ceds_biomass_burning(JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This reads monthly emissions of CEDS biomass burning.
    !//
    !// Input resolution is 0.25x0.25 degrees, which is first interpolated
    !// to 0.5x0.5 degrees, then distributed according
    !// to RETRO vertical distribution (also 0.5x0.5 degrees), before
    !// interpolated to CTM3 resolution.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// Variables info:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Amund Sovde Haslerud, January 2017
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_VARNAME, FF_SCALE
    use ncutils, only: get_netcdf_var_1d, readnc_2d_from3d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,N_EMIS,II,JJ,MP

    !// Dimensions
    integer, parameter :: IRES=360 ,JRES=180, IHLF=720, JHLF=360, LHLF=13
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: HXBEDGE(IHLF+1), HYBEDGE(JHLF+1),HXYBOX(JHLF)
    real(r8)   :: QXBEDGE(IQRT+1), QYBEDGE(JQRT+1),QXYBOX(JQRT)

    real(r8) :: BBH_FIR(IHLF,JHLF,EPAR_FIR_LM), HLFDUM(IHLF,JHLF)
    real(r8) :: R8Q(IQRT,JQRT), RDUM2(IQRT,JQRT)

    !// Interpolated array
    real(r8)  :: RXY8(IPAR,JPAR,EPAR_FIR_LM)

    !// File variables
    logical :: ex
    integer :: nLon, nLat, nTime, sTime
    character(len=200) :: INFILE
    character(len=13) :: infile_daterange
    integer :: NTRNR, getY, startY
    real(r8) :: EDATA(IHLF,JHLF,EPAR_FIR), tscale, mmday, sumA, sumB

    real(r8),dimension(12) :: midmonth = (/15.5_r8, 45.0_r8, 74.5_r8, 105.0_r8, 135.5_r8, 166.0_r8, 196.5_r8, 227.5_r8, 258.0_r8, 288.5_r8, 319.0_r8, 349.5_r8 /)

    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'ceds_biomass_burning'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Biomass burning (CEDS)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// ------------------------------------------------------------------

    !// Largely follow GFED4 read-in.

    !// Year to read
    if (FF_YEAR .ne. 9999) then
       getY = FF_YEAR
    else
       getY = MYEAR !// Use meteorological year
    end if


    !// Open file and read the data...
    if (getY .ge. 1850) then 
       startY = 1850
       infile_daterange = '185001-201512'
    else
       startY = 1750
       infile_daterange = '175001-184912'
    endif



    !// Read file resolution
    infile = trim(FF_PATH)//trim(FF_BNAMES(1))//'-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_'//trim(infile_daterange)//'.nc'


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'longitude',  inLon  )
    call get_netcdf_var_1d( infile, 'latitude',  inLat  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    if (nLon .ne. IQRT .or. nLat .ne. JQRT) then
       write(6,'(a,2i5)') f90file//':'//subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if
    
    !// Find timestep to read
    write(6,'(a,i5)') '  month to fetch: ', JMONTH 
    write(6,'(a,2i5)') '  year to fetch, start year on file:', getY, startY

    mmday = midmonth(JMONTH) + real(getY - startY, r8) * 365._r8
    sTime = -1
    do I = 1, nTime
       if (inTime(I) .eq. mmday) then
          sTime = I
          exit
       end if
    end do

    if (sTime .eq. -1) then 
       write(6,'(a)') f90file//':'//subr//': Wrong inTime index'
       stop
    else 
       write(6,'(a,i5,f12.3)') &
            '  Month index (sTime) and time (inTime(sTime)):',sTime,inTime(sTime)
    end if

    deallocate( inLon, inLat, inTime )
  
    !// Set up grid (start 0E,90S, all input data will be flipped to that.)
    call get_xyedges(nLon,nLat,QXBEDGE,QYBEDGE,1)

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do


    !// Set up halfdegrees grid (start 0E,90S, all input data will be flipped to that.)
    call get_xyedges(IHLF,JHLF,HXBEDGE,HYBEDGE,1)

    !// HLF Grid box areas
    do J = 1, JHLF
       HXYBOX(J) =  A0*A0 * CPI180*(HXBEDGE(2) - HXBEDGE(1)) &
            * (sin(CPI180*HYBEDGE(J+1)) - sin(CPI180*HYBEDGE(J)))
    end do


    !// Get vertical distribution
    !// ------------------------------------------------------------------
    INFILE = trim(FF_PATH)//'bb_altitudes_aggregated.0.5x0.5.nc'
    write(6,'(a)') '  Reading BBH '//trim(INFILE)

    !// Read the BBH data
    call read_retrobbh(INFILE,IHLF,JHLF,EPAR_FIR_LM,BBH_FIR,LPARW,LMMAP)


    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    !// Initialize 3D HLF array (which will be interpolated later)
    EDATA(:,:,:) = 0._r8

    do N_EMIS = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       !// Skip if scaling factor is zero
       if (FF_SCALE(N_EMIS) .eq. 0._r8) cycle

       !// File name, version 1-2
       infile = trim(FF_PATH)//trim(FF_BNAMES(N_EMIS))//'-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_'//trim(infile_daterange)//'.nc'

       inquire(file=infile,exist=ex)
       if (.not. ex) then
          write(6,'(a)') f90file//':'//subr//': no such file: '//trim(infile)
          stop
       end if

       write(6,'(a)') '  Reading '//trim(FF_VARNAME(N_EMIS))//' from '//trim(INFILE)

       call readnc_2d_from3d(INFILE, 'longitude', nlon, 'latitude', nlat,'time', &
            sTime, FF_VARNAME(N_EMIS), R8Q)

       !// Check for and remove missing values
       do J = 1, JQRT
          do I = 1, IQRT
             if (R8Q(I,J) .ge. 1.e+20_r8) then
                R8Q(I,J) = 0._r8
             end if
          end do
       end do

       RDUM2(:,:) = 0._r8
       !// Flip so field starts at (0E,90S)
       do J = 1, JQRT
          RDUM2(:,J) = R8Q(:,JQRT+1-J)
       end do
       do I = 1, IQRT/2
          R8Q(IQRT/2+I,:) = RDUM2(I,:)
       end do
       do I = (IQRT/2)+1, IQRT
          R8Q(I-IQRT/2,:) = RDUM2(I,:)         
       end do



       !// multiply by area for interpolation, kg/m2/s --> kg/s
       do J = 1, JQRT
          R8Q(:,J) = R8Q(:,J) * QXYBOX(J)
       end do


       HLFDUM(:,:) = 0._r8
       !// From 0.25 to 0.5 degree
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, HLFDUM, &
                 HXBEDGE,HYBEDGE,IHLF,JHLF,1)

       sumB = sum(R8Q)
       sumA = sum(HLFDUM)
       if (abs(sumB-sumA)/sumB .gt. 1.e-5_r8) then
          write(6,'(a,es20.12,es20.12)') f90file//':'//subr// &
               ':  Wrong sum before/after interp. ', sumB, sumA
          stop
       end if

       EDATA(:,:,N_EMIS) = HLFDUM(:,:)

    end do !// do N_EMIS = 1, NEFIR

    !// Distribute vertically and interpolate to model resolution
    !// ------------------------------------------------------------------
    do N_EMIS = 1, NEFIR

      !// If component is not included, go to next N_EMIS
      if (ECOMP_FIR(N_EMIS) .le. 0) cycle

      NTRNR = ECOMP_FIR(N_EMIS)

       !// Distribute with height
!$omp parallel private (L,I,J,HLFDUM) &
!$omp          shared (EDATA,BBH_FIR,RXY8,LMMAP,HXBEDGE,HYBEDGE, &
!$omp                  XDEDG,YDEDG,N_EMIS) &
!$omp          default(NONE)
!$omp do
      do L = 1, EPAR_FIR_LM
         HLFDUM(:,:) = EDATA(:,:,N_EMIS)
         !// Seems BBH_FIR may be zero for some EDATA. If so, we put
         !// emissions at the surface.
         if (maxval(HLFDUM) .gt. 0._r8) then
            HLFDUM(:,:) = HLFDUM(:,:) * BBH_FIR(:,:,LMMAP(L))
            if (L .eq. 1) then
               do J = 1, JHLF
                  do I = 1, IHLF
                     if (sum(BBH_FIR(I,J,:)) .eq. 0._r8) &
                           HLFDUM(I,J) = EDATA(I,J,N_EMIS)
                  end do
               end do
            end if
            !// No need to check the other levels. Sum of BBH_FIR in vertical
            !// is always 0 or 1.

            !// Interpolate
            call E_GRID(HLFDUM,HXBEDGE,HYBEDGE,IHLF,JHLF, RXY8(:,:,L), &
                 XDEDG,YDEDG,IPAR,JPAR,1)

         else
            RXY8(:,:,L) = 0._r8
         end if
      end do !// do L = 1, EPAR_FIR_LM
!$omp end do
!$omp end parallel



      !// Put into EMIS_FIR and scale it
!$omp parallel private (MP,I,J,L,II,JJ) &
!$omp         shared (RXY8,EMIS_FIR,MPBLKJB,MPBLKJE,MPBLKIB,MPBLKIE,N_EMIS,FF_SCALE) &
!$omp         default(NONE)
!$omp do
      do MP = 1, MPBLK
        !// Loop over latitude (J is global, JJ is block)
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1

          !// Loop over longitude (I is global, II is block)
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1

            do L = 1, EPAR_FIR_LM
               EMIS_FIR(L,N_EMIS,II,JJ,MP) = RXY8(I,J,L) * FF_SCALE(N_EMIS)
            end do
          end do
        end do
      end do

!$omp end do
!$omp end parallel
      tscale = 0._r8
      do MP = 1, MPBLK
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = 1, EPAR_FIR_LM
               tscale = tscale + EMIS_FIR(L,N_EMIS,II,JJ,MP)
            end do
          end do
        end do
      end do
      write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ',FF_CNAMES(N_emis),NTRNR,&
           tscale*2628000.e-9_r8,' Tg/mon'

    end do !// do do N_EMIS = 1, NEFIR

    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'




    !// ------------------------------------------------------------------
  end subroutine ceds_biomass_burning
  !// ----------------------------------------------------------------------


  !// ------------------------------------------------------------------
  subroutine ceds_biomass_burning_novert(JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This reads monthly emissions of CEDS biomass burning.
    !//
    !// Input resolution is 0.25x0.25 degrees, which is first interpolated
    !// to 0.5x0.5 degrees, then distributed according
    !// to RETRO vertical distribution (also 0.5x0.5 degrees), before
    !// interpolated to CTM3 resolution.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// Variables info:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Amund Sovde Haslerud, January 2017
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_VARNAME, FF_SCALE
    use ncutils, only: get_netcdf_var_1d, readnc_2d_from3d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,N_EMIS,II,JJ,MP

    !// Dimensions
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: QXBEDGE(IQRT+1), QYBEDGE(JQRT+1),QXYBOX(JQRT)

    real(r8) :: R8Q(IQRT,JQRT), RDUM2(IQRT,JQRT)

    !// Interpolated array
    real(r8)  :: RXY8(IPAR,JPAR,1)

    !// File variables
    logical :: ex
    integer :: nLon, nLat, nTime, sTime
    character(len=200) :: INFILE
    character(len=13) :: infile_daterange
    integer :: NTRNR, getY, startY
    real(r8) :: mmday, sumA, sumB

    real(r8),dimension(12) :: midmonth = (/15.5_r8, 45.0_r8, 74.5_r8, 105.0_r8, 135.5_r8, 166.0_r8, 196.5_r8, 227.5_r8, 258.0_r8, 288.5_r8, 319.0_r8, 349.5_r8 /)

    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'ceds_biomass_burning_novert'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Biomass burning (CEDS)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// ------------------------------------------------------------------

    !// Largely follow GFED4 read-in.

    !// Year to read
    if (FF_YEAR .ne. 9999) then
       getY = FF_YEAR
    else
       getY = MYEAR !// Use meteorological year
    end if


    !// Open file and read the data...
    if (getY .ge. 1850) then 
       startY = 1850
       infile_daterange = '185001-201512'
    else
       startY = 1750
       infile_daterange = '175001-184912'
    endif



    !// Read file resolution
    infile = trim(FF_PATH)//trim(FF_BNAMES(1))//'-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_'//trim(infile_daterange)//'.nc'


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'longitude',  inLon  )
    call get_netcdf_var_1d( infile, 'latitude',  inLat  )
    call get_netcdf_var_1d( infile, 'time', inTime )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )
    nTime = SIZE( inTime )

    if (nLon .ne. IQRT .or. nLat .ne. JQRT) then
       write(6,'(a,2i5)') f90file//':'//subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if
    
    !// Find timestep to read
    write(6,'(a,i5)') '  month to fetch: ', JMONTH 
    write(6,'(a,2i5)') '  year to fetch, start year on file:', getY, startY

    mmday = midmonth(JMONTH) + real(getY - startY, r8) * 365._r8
    sTime = -1
    do I = 1, nTime
       if (inTime(I) .eq. mmday) then
          sTime = I
          exit
       end if
    end do

    if (sTime .eq. -1) then 
       write(6,'(a)') f90file//':'//subr//': Wrong inTime index'
       stop
    else 
       write(6,'(a,i5,f12.3)') &
            '  Month index (sTime) and time (inTime(sTime)):', &
            sTime,inTime(sTime)
    end if

    deallocate( inLon, inLat, inTime )
  
    !// Set up grid (start 0E,90S, all input data will be flipped to that.)
    call get_xyedges(nLon,nLat,QXBEDGE,QYBEDGE,1)

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do


    !// No vertical distribution
    !// ------------------------------------------------------------------


    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    do N_EMIS = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       !// Skip if scaling factor is zero
       if (FF_SCALE(N_EMIS) .eq. 0._r8) cycle

       !// File name, version 1-2
       infile = trim(FF_PATH)//trim(FF_BNAMES(N_EMIS))//'-em-biomassburning_input4MIPs_emissions_CMIP_VUA-CMIP-BB4CMIP6-1-2_gn_'//trim(infile_daterange)//'.nc'

       inquire(file=infile,exist=ex)
       if (.not. ex) then
          write(6,'(a)') f90file//':'//subr//': no such file: '//trim(infile)
          stop
       end if

       write(6,'(a)') '  Reading '//trim(FF_VARNAME(N_EMIS))//' from '//trim(INFILE)

       call readnc_2d_from3d(INFILE, 'longitude', nlon, 'latitude', nlat,'time', &
            sTime, FF_VARNAME(N_EMIS), R8Q)

       !// Check for and remove missing values
       do J = 1, JQRT
          do I = 1, IQRT
             if (R8Q(I,J) .ge. 1.e+20_r8) then
                R8Q(I,J) = 0._r8
             end if
          end do
       end do

       RDUM2(:,:) = 0._r8
       !// Flip so field starts at (0E,90S)
       do J = 1, JQRT
          RDUM2(:,J) = R8Q(:,JQRT+1-J)
       end do
       do I = 1, IQRT/2
          R8Q(IQRT/2+I,:) = RDUM2(I,:)
       end do
       do I = (IQRT/2)+1, IQRT
          R8Q(I-IQRT/2,:) = RDUM2(I,:)         
       end do



       !// multiply by area for interpolation, kg/m2/s --> kg/s
       do J = 1, JQRT
          R8Q(:,J) = R8Q(:,J) * QXYBOX(J)
       end do

       !// From 0.25 to model resolution
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, RXY8(:,:,1), &
                 XDEDG,YDEDG,IPAR,JPAR,1)


       sumB = sum(R8Q)
       sumA = sum(RXY8)
       if (abs(sumB-sumA)/sumB .gt. 1.e-5_r8) then
          write(6,'(a,es20.12,es20.12)') f90file//':'//subr// &
               ':  Wrong sum before/after interp. ', sumB, sumA
          stop
       end if

       !// Put into EMIS_FIR and scale it [kg/s]
       EMIS_FIR(1,N_EMIS,II,JJ,MP) = RXY8(I,J,1) * FF_SCALE(N_EMIS)

       write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ', &
            FF_CNAMES(N_emis),NTRNR,&
            sumB * 2628000.e-9_r8,' Tg/mon'

    end do !// do N_EMIS = 1, NEFIR

    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'


    !// ------------------------------------------------------------------
  end subroutine ceds_biomass_burning_novert
  !// ----------------------------------------------------------------------





  !// ------------------------------------------------------------------
  subroutine gfed4_rd_daily(JDATE,JMONTH,JYEAR)
    !// ------------------------------------------------------------------
    !// This reads monthly emissions of dry matter (DM) from GFEDv4,
    !// for the six partitions:
    !// SAVA: Savanna
    !// BORF: Boreal forrest
    !// TEMF: Temperature forest
    !// DEFO: Deforestation
    !// PEAT: Peat
    !// AGRI: Agricultural wastet burning
    !//
    !// Input resolution is 0.25x0.25 degrees, which is first interpolated
    !// to 0.5x0.5 degrees, then distributed according
    !// to RETRO vertical distribution (also 0.5x0.5 degrees), before
    !// interpolated to CTM3 resolution.
    !//
    !// Interpolation is PARALLELLIZED!
    !//
    !// Not yet configured for daily fraction of monthly emissions.
    !// 
    !// The routine is hardcoded to read separate emission array (EMIS_FIR).
    !// This is necessary due to all scalings for NMVHCs, and because the
    !// emission fields need to be independent of emission arrays
    !// E2DS and E3DS to allow for more frequent updating.
    !//
    !// Variables/parameters defined in oc_globalvariables.f:
    !//   EPAR_FIR:    Max number of components in EMIS_FIR
    !//   NEFIR:       The number of components in EMIS_FIR
    !//   EPAR_FIR_LM: Number of layers in emission set, starting from surface.
    !//   ECOMP_FIR:   Transport numbers of the components in EMIS_FIR
    !//   EMIS_FIR:    Size (LPAR,EPAR_FIR,IDBLK,JDBLK,MPBLK)
    !//
    !// Amund Sovde Haslerud, December 2017
    !//   Added daily fractions.
    !// Ole Amund Sovde, February 2012, updated June 2013
    !// ------------------------------------------------------------------
    use cmn_size, only: IPAR, JPAR, LPAR, LPARW, NPAR, MPBLK
    use cmn_ctm, only: LMMAP, XDEDG, YDEDG, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use cmn_oslo, only: EPAR_FIR, NEFIR, EPAR_FIR_LM, ECOMP_FIR, EMIS_FIR, &
         FF_PATH, FF_YEAR, METHANEMIS, &
         FF_CNAMES, FF_BNAMES, FF_SCALE, FF_PARTITIONS, &
         GFED4_NC, GFED4_NP, GFED4_EF, GFED4_NM, DINM
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JDATE, JMONTH, JYEAR

    !// Locals
    integer :: I,J,L,M,N,N_EMIS,II,JJ,MP, io_err

    !// Dimensions
    integer, parameter :: IRES=360 ,JRES=180, IHLF=720, JHLF=360, LHLF=13
    integer, parameter :: IQRT=1440, JQRT=720
    !// XY grid
    real(r8)   :: HXBEDGE(IHLF+1), HYBEDGE(JHLF+1),HXYBOX(JHLF)
    real(r8)   :: QXYBOX(JQRT)

    real(r8) :: BBH_FIR(IHLF,JHLF,EPAR_FIR_LM)
    real(r8), dimension(IQRT,JQRT) :: R8Q, QDUM, DFRAC
    real(r8), dimension(IHLF,JHLF) :: HDUM

    character(len=20) :: COMP
    character(len=2) :: CDATE, CMONTH
    character(len=4) :: CYEAR
    !// Partitions: Agriculture, deforestation, forest, peat, savanna, woodland
    integer :: partitions(GFED4_NP)
    real(r8)  :: DM_PARTS(IHLF,JHLF,GFED4_NP)
    character(len=4) :: CPARTS(GFED4_NP)

    !// Interpolated array
    real(r8)  :: RXY8(IPAR,JPAR,EPAR_FIR_LM)

    !// File variables
    logical :: file_status
    integer :: efnr, nLon, nLat
    character(len=120) :: INFILE
    real(r8) :: EDATA(IHLF,JHLF),tscale, scalefac

    integer :: GFED4_N

    real(r8), allocatable, dimension(:) :: QXBEDGE, QYBEDGE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'gfed4_rd_daily'
    !// ------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': Biomass burning (GFEDv4)'

    !// Initialize emissions
    EMIS_FIR(:,:,:,:,:) = 0._r8

    !// DM emissions from each partition, to be read from file
    DM_PARTS(:,:,:) = 0._r8

    !// ------------------------------------------------------------------

    !// Select partitions to use (default is all)
    CPARTS = (/'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'/)


    !// Read DM from partitions
    !// ------------------------------------------------------------------
    !// Date and month to read
    write(CDATE(1:2),'(i2.2)') JDATE
    write(CMONTH(1:2),'(i2.2)') JMONTH
    !// Year to read
    if (FF_YEAR .ne. 9999) then
       write(CYEAR,'(i4.4)') FF_YEAR
    else
       write(CYEAR,'(i4.4)') MYEAR !// Use meteorological year
    end if


    !// Include daily fraction of monthly totals
    !// ------------------------------------------------------------------
    infile = trim(FF_PATH)//'daily/'//CYEAR// &
         '/fraction_emissions_'//CYEAR//CMONTH//CDATE//'_nc4.nc'
    write(6,'(a)') f90file//':'//subr//': Reading '//trim(INFILE)
    !// Emission variable, expected units: kg/m2/s for each month
    call get_netcdf_var_2d(infile, 'Fraction_of_Emissions', QDUM, IQRT,JQRT)

    !// IMPORTANT: Must scale by days in month, since emissions
    !//            are given in units kg/s.
    !// Skip possible negatives.
    QDUM(:,:) = max(0._r8, QDUM(:,:)) * real(DINM(JMONTH),r8)
    !// Switch from -180:180 to 0:360, as GFED4 monthly fields
    do J = 1, JQRT
       DFRAC(1:IQRT/2,J) = QDUM(IQRT/2+1:IQRT,J)
       DFRAC(IQRT/2+1:IQRT,J) = QDUM(1:IQRT/2,J)
    end do

    write(6,'(a,2es16.6)') f90file//':'//subr//': DFRAC min/max',minval(DFRAC),maxval(DFRAC)

    !// Filename
    infile = trim(FF_PATH)//'GFED4.1s_'//CYEAR//'_'//CMONTH//'.nc'

    !// Inquire file
    inquire(file=infile,exist=file_status)

    if (.not. file_status) then
       write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
            trim(infile)
       !// Fall-back to 2003/2016
       if (MYEAR .gt. 2016) then
          infile = trim(FF_PATH)//'GFED4.1s_2016_'//CMONTH//'.nc'
       else if (MYEAR .lt. 2003) then
          infile = trim(FF_PATH)//'GFED4.1s_2003_'//CMONTH//'.nc'
       else
          stop
       end if
       !// In case of fall-back
       write(6,'(a)') 'Will try to use file: '//trim(infile)
       inquire(file=infile,exist=file_status)
       if (.not. file_status) then
          write(6,'(a)') f90file//':'//subr//': file does not exist: '// &
               trim(infile)
          stop
       end if
    end if


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lone', QXBEDGE )
    call get_netcdf_var_1d( infile, 'late', QYBEDGE )

    nLon  = SIZE( QXBEDGE ) - 1
    nLat  = SIZE( QYBEDGE ) - 1

    if (nLon .ne. IQRT .or. nLat.ne.JQRT) then
       write(6,'(a,2i5)') f90file//':'//subr//': Wrong lon/lat size',nLon,nLat
       write(6,'(a,2i5)') '  Should be:',IQRT,JQRT
       stop
    end if

    !// QRT Grid box areas
    do J = 1, JQRT
       QXYBOX(J) =  A0*A0 * CPI180 * (QXBEDGE(2) - QXBEDGE(1)) &
            * (sin(CPI180*QYBEDGE(J+1)) - sin(CPI180*QYBEDGE(J)))
    end do


    !// Set up halfdegrees (Grid type is 1, start at 0E, 90S)
    !// (Actually, GFED starts at 180W, 90N, but it is swicthed after read-in.)
    call get_xyedges(IHLF,JHLF,HXBEDGE,HYBEDGE,1)

    !// HLF Grid box areas
    do J = 1, JHLF
       HXYBOX(J) =  A0*A0 * CPI180*(HXBEDGE(2) - HXBEDGE(1)) &
            * (sin(CPI180*HYBEDGE(J+1)) - sin(CPI180*HYBEDGE(J)))
    end do




    !// Read DM for six categories; save data on 0.5x0.5 degree
    !// For each species, add all partitions together in 0.5x0.5.
    !// Finally convert to model resolution.
    do N = 1, GFED4_NP

       !// Emission variable, expected units: kg/m2/s for each month
       !// GFED4 on netcdf has been converted to start at 0E,90S
       call get_netcdf_var_2d( infile, 'DM_'//CPARTS(N),R8Q(:,:), nlon,nlat )

       !// Multiply by daily fraction of month
       R8Q(:,:) = R8Q(:,:) * DFRAC(:,:)
       !// multiply by area for interpolation
       do J = 1, JQRT
          R8Q(:,J) = R8Q(:,J) * QXYBOX(J)
       end do

       !// From 0.25 to 0.5 degree
       call E_GRID(R8Q,QXBEDGE,QYBEDGE,IQRT,JQRT, HDUM, &
                 HXBEDGE,HYBEDGE,IHLF,JHLF,1)

       !// divide by area after interpolation
       do J = 1, JHLF
          HDUM(:,J) = HDUM(:,J) / HXYBOX(J)
       end do

       !// Save DM for each partition
       DM_PARTS(:,:,N) = HDUM(:,:)

    end do !// do N = 1, GFED4_NP

    !// Deallocate allocated variables
    if ( allocated(QXBEDGE) ) deallocate(QXBEDGE)
    if ( allocated(QYBEDGE) ) deallocate(QYBEDGE)


    !// Get vertical distribution
    !// ------------------------------------------------------------------
    INFILE = trim(FF_PATH)//'bb_altitudes_aggregated.0.5x0.5.nc'
    write(6,'(a)') ' * Reading BBH '//trim(INFILE)

    !// Read the BBH data
    call read_retrobbh(INFILE,IHLF,JHLF,EPAR_FIR_LM,BBH_FIR,LPARW,LMMAP)

    !// Originally, there was a check whether BBH(I,J,1) was 0 but non-zero
    !// above. This never happened, so I removed the test.


    !// Collect emissions of specified partitionings, for each of
    !// the emitted components.
    !// ------------------------------------------------------------------

    !// Loop through componens and interpolate
    EDATA(:,:) = 0._r8

    do N_EMIS = 1, NEFIR

       !// Skip component if not included
       if (ECOMP_FIR(N_EMIS) .le. 0) cycle

       !// Find corresponding GFED4 name. Some use same as CTM3, others
       !// are different and must be overrided below.
       COMP = trim(FF_BNAMES(N_EMIS))
       !// Additional scaling factor
       scalefac = FF_SCALE(N_EMIS)

       !// Select partitions to use (default is all)
       !// 'SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'
       partitions(:) = FF_PARTITIONS(:,N_EMIS)


       !// Find the GFED4 component number, if available.
       GFED4_N = -1
       do I = 1, GFED4_NC
          if (trim(GFED4_NM(I)) .eq. trim(COMP)) then
             GFED4_N = I
             exit
          end if
       end do

       !// If not available, skip component
       if (GFED4_N .le. 0) cycle

       !// Sum up partitions for this species
       HDUM(:,:) = 0._r8
       do N = 1, GFED4_NP

          if (partitions(N) .eq. 0._r8) cycle
          !// Emission factors are g/kgDM, so we convert to kg, and also
          !// apply the scalefac.
          tscale =  GFED4_EF(N,GFED4_N) * 1.e-3_r8 * scalefac

          HDUM(:,:) = HDUM(:,:) + DM_PARTS(:,:,N) * tscale

       end do

       !// Convert to kg/s
       do J = 1, JHLF
          EDATA(:,J) = HDUM(:,J) * HXYBOX(J)
       end do



       !// Distribute with height
       do L = 1, EPAR_FIR_LM
          HDUM(:,:) = EDATA(:,:) * BBH_FIR(:,:,LMMAP(L))
          if (maxval(HDUM) .gt. 0._r8) then
             if (L.eq.1) then
                !// Check problems with BBH on half degree vs emissions on
                !// one degree
                do J = 1, JHLF
                  do I = 1, IHLF
                     !// Put everyting in layer 1
                     if (BBH_FIR(I,J,1) .eq. 0._r8 .and. &
                          HDUM(I,J) .gt. 0._r8) then
                        HDUM(I,J) = EDATA(I,J)
                     end if
                  end do
               end do
            end if

            !// Interpolate
            call E_GRID(HDUM,HXBEDGE,HYBEDGE,IHLF,JHLF, RXY8(:,:,L), &
                 XDEDG,YDEDG,IPAR,JPAR,1)
         else
            RXY8(:,:,L) = 0._r8
         end if
      end do !// do L = 1, EPAR_FIR_LM

      do MP = 1, MPBLK
        !// Loop over latitude (J is global, JJ is block)
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1

          !// Loop over longitude (I is global, II is block)
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1

            do L = 1, EPAR_FIR_LM
               EMIS_FIR(L,N_EMIS,II,JJ,MP) = RXY8(I,J,L)
            end do
          end do
        end do
      end do

      tscale = 0._r8
      do MP = 1, MPBLK
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = 1, EPAR_FIR_LM
               tscale = tscale + EMIS_FIR(L,N_EMIS,II,JJ,MP)
            end do
          end do
        end do
      end do
      write(6,'(a,a10,1x,i3,es12.5,a)') ' * Total FF ',FF_CNAMES(N_emis),ECOMP_FIR(N_EMIS),&
           tscale*86400.e-9_r8,' Tg/day'

    end do !// do do N_EMIS = 1, NEFIR

    write(6,'(a)') f90file//':'//subr//': Forest fires emissions are updated'

    !// ------------------------------------------------------------------
  end subroutine gfed4_rd_daily
  !// ------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getmeganmonthly(INFILE,INVAR,R8XY,IRES,JRES,MRES,NSETS,&
       XBEDGE,YBEDGE)
    !// --------------------------------------------------------------------
    !// Reads MEGAN accumulated surface emissions,
    !// fetching variable name INVAR.
    !//
    !// Amund Sovde Haslerud, January 2018
    !// --------------------------------------------------------------------
    use netcdf
    use ncutils, only: get_netcdf_var_0d, get_netcdf_var_1d, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
 
    !// Input
    character(len=*), intent(in)  :: infile   !// Name of netCDFfile
    character(len=*), intent(in)  :: invar    !// Name of input variable
    integer,intent(in) :: IRES     !// Number of longitudes
    integer,intent(in) :: JRES     !// Number of latitudes
    integer,intent(in) :: MRES     !// Number of sets (months) per category

    !// Output
    real(r8), intent(out) :: R8XY(IRES,JRES,MRES)  !// 3D field
    integer,intent(out)   :: NSETS
    real(r8), intent(out) :: XBEDGE(IRES+1), YBEDGE(JRES+1)

    !// Locals
    real(r8) :: delta_time
    integer            :: status, ncid
    integer            :: nLon, nLat, I, J
    real(r8), allocatable, dimension(:) :: inLon, inLat, inLonE, inLatE
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'nc_getdata_sfc_2d'
    !// ------------------------------------------------------------------


    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'ilon',  inLonE  )
    call get_netcdf_var_1d( infile, 'ilat',  inLatE  )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )

    !// Set XBEDGE and YBEDGE
    XBEDGE(:) = inLonE(:)
    YBEDGE(:) = inLatE(:)

    !// Deallocate all local variables
    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)
    if ( allocated(inLonE) ) deallocate(inLonE)
    if ( allocated(inLatE) ) deallocate(inLatE)

    !// Get delta_time
    call get_netcdf_var_0d( infile, 'delta_time',  delta_time  )

    write(6,'(a,es16.6)') f90file//':'//subr// &
         ': delta_time: ',delta_time

    if (nLon.ne.IRES .or. nLat.ne.JRES) then
       write(6,'(a)') f90file//':'//subr// &
            ': Horizontal resolution on file does not match'// &
            ' specified resolution'
       write(6,'(a,2i5)') '  Specified:',IRES,JRES
       write(6,'(a,2i5)') '  On file:  ',nLon,nLat
       stop
    end if

       

    !// Emission variable, expected units: kg/m2/s for each month
    call get_netcdf_var_2d( infile, invar, R8XY(:,:,1), nlon,nlat )

    !// Check for nans
    do J = 1, nlat
       do I = 1, nlon
          if (R8XY(I,J,1) .ne. R8XY(I,J,1)) R8XY(I,J,1) = 0._r8
       end do
    end do

    !// Make field kg/s
    R8XY(:,:,:) = R8XY(:,:,:) / delta_time

    !// Number of datasets
    write(6,'(a)') f90file//':'//subr//': One dataset'
    NSETS = 1

    !// --------------------------------------------------------------------
  end subroutine getmeganmonthly
  !// ----------------------------------------------------------------------






  !// ----------------------------------------------------------------------
end module emisutils_oslo
!//=========================================================================

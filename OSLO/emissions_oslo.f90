!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Søvde Haslerud, January 2017
!//=========================================================================
!// Read in emissions.
!//=========================================================================
module emissions_oslo
  !// ----------------------------------------------------------------------
  !// Routines for reading in emissions.
  !// Contains:
  !//   subroutine emis_input: Reading in emis data.
  !//   subroutine update_emis: update short term emissions.
  !//   subroutine emis_prop: properties for the emission dataset
  !//   subroutine emis_unit: calculate scale factor for dataset
  !//   subroutine emis_define_2dscaling: define scaling flags
  !//
  !// Amund Søvde Haslerud, January 2017
  !// Ole Amund Sovde, November 2015, August - September 2009
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'emissions_oslo.f90'
  !// ----------------------------------------------------------------------
  public
  private emis_prop, emis_unit, emis_define_2dscaling
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine emis_input()
    !// --------------------------------------------------------------------
    !// Read in emission information for CTM3. This is read from the file
    !// Ltracer_emis.inp, and consists of several sections.
    !// 1. How to treat biomass burning.
    !// 2. Which aircraft inventory to use.
    !// 3. Lightning NOx emissions.
    !// 4. 2D (surface) emissions, monthly. Can put diurnal cycles on top.
    !// 5. 3D emissions, monthly.
    !// 6. Emissions on shorter term basis.
    !//
    !// Basic info about each dataset are read, e.g. read-in format code
    !// (e.g. 501), scaling factors, category name, info on diurnal variations
    !//  etc. See manual for example from a 2D dataset.
    !//
    !// The structure of the input file is explained at the end of the
    !// input file.
    !//
    !// Ole Amund Sovde, April 2015, August - September 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_ctm, only: NTM, NPAR
    use cmn_chem, only: INFILE_EMIS, NEMLIT, NLIT, LITFAC, TNAME, &
         TMASS, NE2TBL, E2LTBL, E2STBL, NM2TBL, NY2TBL, E2DS, &
         NE3TBL, E3LTBL, E3STBL, NM3TBL, NY3TBL, E3DSNEW
    use cmn_parameters, only: AVOGNR
    use utilities, only: get_free_fileid
    use cmn_oslo, only: E2LocHourTBL, E2CTBL, &
         FF_TYPE, FF_YEAR, FF_PATH, EMIS_FIR, &
         ECOMP_FIR, E22dTBL, E22dSCALE, E2L2dTBL, METHANEMIS, &
         EPAR_FIR, FF_CNAMES, FF_BNAMES, EPAR_NPARTS, FF_PARTITIONS, &
         FF_SCALE, NEFIR, FF_VARNAME
    use emisutils_oslo, only: reademis_2d, reademis_3d, reademis_stv, &
         set_diurnal_scalings, set_vertical_scalings, gfed4_init
    use sulphur_oslo, only: DMSseaconc
    use emissions_aircraft, only: AirScen, AirScenYear, AirEmisPath
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// File variables
    integer :: ifnr

    !// Indices for loops etc.
    integer :: C, M,  EC,NC, io_err, MTBL_OLD, MTBL,NTR

    !// Data to be read from emission input file
    real(r8)  :: EFACT, EMFACM, ESCALE
    real(r8) :: partitions(EPAR_NPARTS)
    integer :: ETAG, EYEAR, EMON, ETYP, EUNIT, ESCENYEAR
    integer :: EDIURN, ESCAL_LH, ESCAL_2D
    integer :: EVERT
    character(len=3)   :: ECAT !// Category name for data
    character(len=3)   :: ERES !// Resolution for data
    character(len=10)  :: ENAME
    character(len=20)  :: BNAME, EVAR
    character(len=50)  :: EDATASET
    character(len=120) :: ELINE, ELINE2
    character(len=200) :: EFILE


    !// input dataset variables
    integer, parameter :: MRES = 12   !// Max number of data sets/months
    integer :: MIDX(MRES)   !// Emission table index for set/cat
    integer :: &
         NSETS, &           !// Number of sets on file
         N1, N2, &          !// Indices
         IEMIS,JEMIS, &     !// Input resolution, set from ERES
         J, N               !// For lightning read-in

    !// Tracer scalings
    integer  :: ES_TRNR(NPAR)     !// Transport numbers for these tracers
    integer  :: ES_NT             !// Number of tracers to use the dataset
    real(r8) :: USCALE            !// Unit scaling for tracers
    real(r8) :: ES_USCALE(NPAR)   !// Total scalings for tracers

    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_input'
    !// --------------------------------------------------------------------


    !// Initialize diurnal variations
    call set_diurnal_scalings()

    !// Initialize vertical scalings
    call set_vertical_scalings()

    !// Find non-used file number for input file
    ifnr = get_free_fileid()

    !// Open emission file
    write(6,'(a)') f90file//':'//subr// &
         ': Reading emissions from '//trim(infile_emis)
    open(ifnr, file=INFILE_EMIS, status='old', form='formatted', iostat=io_err)
    if (io_err .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': emis_input: Error reading ' &
            //trim(INFILE_EMIS)
       write(6,'(a,i5)') '    netCDF error code:',io_err
       stop
    end if

    !// Read the whole file and print to screen - this is to have a more
    !// easy way to inspect what the run used.
    write(6,'(a)') '  First print the whole list to std.out.'
    write(6,'(a)') '  Then the datasets are read one by one.'
    do while (.true.)
       read(ifnr,'(a200)') EFILE
       write(6,'(a)') trim(EFILE)
       if (EFILE(1:5) .eq. '-END-') exit
    end do
    !// Close and re-open file to start at beginning.
    close(ifnr)
    open(ifnr, file=INFILE_EMIS, status='old', form='formatted', iostat=io_err)


    !// main header
    read(ifnr,*) ELINE
    write(6,'(a)') ELINE
    !// --------------------------------------------------------------------

    !// General information
    do J = 1,6
       read(ifnr,*) ELINE
    end do


    !// Read in aircraft emissions
    !// --------------------------------------------------------------------
    read(ifnr,*) !// skip header
    read(ifnr,*) ELINE, AirScen, AirScenYear
    read(ifnr,*) ELINE, AirEmisPath

    write(6,'(a)') 'Aircraft emissions scenario: '//trim(AirScen)
    write(6,'(a)') 'Aircraft emissions year: '//trim(AirScen)
    write(6,'(a)') 'Aircraft emissions path: '//trim(AirEmisPath)

    !// Read in lightning emissions
    !// --------------------------------------------------------------------
    read(ifnr,*) !// skip line
    read(ifnr,'(i3)') NEMLIT
    NLIT(:)  = 0
    LITFAC(:) = 0._r8
    J = 0
    do M = 1,NEMLIT
       read(ifnr,*) ENAME, EMFACM, EFACT
       do N = 1,NTM
          if (ENAME.eq.TNAME(N) .and. EFACT.ne.0._r8) then
             J  = J + 1
             NLIT(J)   = N
             !// Scale LITFAC from Tg to kg
             LITFAC(J) = EFACT * TMASS(N) * 1.e9_r8 / EMFACM
          end if
       end do
    end do
    NEMLIT = J


    !// 2D emission tables
    !// --------------------------------------------------------------------
    !// The tables tell us which tracers use the different datasets, and
    !// which month or year they are to be used.

    !// Initialize
    NE2TBL = 0            !// no. of 2-D tables read
    E2LTBL(:,:) = .false. !// .T. = species N uses tbl E
    E2STBL(:,:) = 0._r8   !// scale factor for each (N,E)
    NM2TBL(:) = 0         !// month table for dataset (E)
    NY2TBL(:) = 0         !// year table for dataset (E)
    MTBL = 0              !// emis. table idx (the # of emis. datasets)
    E2DS(:,:,:,:) = 0._r8 !// The emissions

    !// Initialize CTM3 specific tables
    E2CTBL(:) = 0 !// Category number for 2D datasets

    !// Index for diurnal variation for 2D datasets
    E2LocHourTBL(:) = 0

    !// CTM3 2D scaling for 2D emissions.
    E22dTBL(:) = 0             !// Index for hourly 2D variation datasets
    E22dSCALE(:,:,:,:) = 0._r8 !// 2D scale (IDBLK,JDBLK,MPBLK,NE22dVARS)
    E2L2dTBL(:) = .false.      !// Flag which 2d sets are used.


    !// The following procedure is:
    !// 1. Read file name of dataset
    !// 2. Read info about dataset
    !// 3. Read info about tracers using this dataset
    !// 4. Send this info into read-in routine
    !//    a. Dataset is read
    !//    b. Dataset is scaled
    !//    c. Set scaling factors for the whole dataset
    !//    d. Set tracer specific scaling factors for this dataset.

    !// Read header for 2D
    do J = 1,4
       read(ifnr,*) ELINE
    end do
    !// 1. File name of the first dataset (or 'end2D')
    read(ifnr,*) EFILE

    !// Read 2D data
    !// --------------------------------------------------------------------
    do while (EFILE .ne. 'end2D')

       !// 2. Read info about the dataset.
       read(ifnr,*) ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EVERT, EDATASET, ESCENYEAR
       !// What was read:
       !// ETAG:   Read-in format code
       !// EFACT:  Dataset scaling factor
       !// ERES:   Resolution
       !// EMON:   Month for which dataset applies
       !// EYEAR:  Year for which dataset applies
       !// ECAT:   Category on file
       !// ETYP:   Area weighted type for dataset
       !// EUNIT:  Unit flag/number for dataset
       !// EDIURN: Flag for imposing a local hour/diurnal variation
       !//         Will set either ESCAL_LH for local hour or
       !//         ESCAL_2D for 2D field.
       !// EVERT:  Flag for distributing emissions in a few layers.
       !// EDATASET: Dataset name (especilly for netCDF files)
       !// ESCENYEAR: Dataset/scenario year, i.e. the dataset fetched
       !//            from file. Included to make it possible to use
       !//            a specific year of emissions for another meteorological
       !//            year (EYEAR).
       write(6,'(a)') ' 2D emis file:'
       write(6,'(a)') trim(EFILE)
       write(6,'(i5,es14.6,a5,2i5,1x,a3,1x,i3,1x,3i3,1x,a40,1x,i4)') &
            ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EVERT, EDATASET, ESCENYEAR

       !// Set up ESCAL_LH and ESCAL_2D based on EDIURN
       call emis_define_2dscaling(EDIURN, ESCAL_LH, ESCAL_2D)

       !// Keep track of table index befor reading new datasets
       MTBL_OLD = MTBL

       !// Find emission dataset properties; resolution IEMIS,JEMIS
       call emis_prop(ERES,IEMIS,JEMIS)


       !// 3. Read in which tracers that will use this dataset, along
       !//    with the scalings to apply.
       !// -----------------------------------------------------------------
       ES_NT = 0             !// Number of tracers to use the dataset
       ES_TRNR(:) = 0        !// Transport number of these tracers
       ES_USCALE(:)  = 1._r8 !// Scaling factor e.g. due to units

       !// Read component name and component scaling factor (or 'xxx 0')
       read(ifnr,*) ENAME, ESCALE

       !// Check for a component using the dataset ('xxx' = no more species)
       !// -----------------------------------------------------------------
       do while (ENAME .ne. 'xxx')

          !// Find transport number
          do NTR = 1, NPAR
             if (trim(ENAME) .eq. trim(TNAME(NTR))) exit
          end do
          if (NTR .gt. NPAR) then
             !// Wrong species in emission list
             write(6,'(a)') '* Trying to emit a non-included species: ' &
                  //trim(ENAME)
             write(6,'(a)') '  * Skipping emission for tracer '//trim(ENAME)
             !// Skip it; read next tracer to use the data
             read(ifnr,*) ENAME, ESCALE
             cycle
          end if

          !// Check if CH4 is emitted or fixed at surface
          if (trim(ENAME).eq.'CH4' .and. .not.METHANEMIS) then
             write(6,'(a)') f90file//':'//subr// &
                  ': Cannot emit CH4 when fixed at surface!'
             stop 'STOP in '//subr
          end if


          !// Keep track how many species will use this 2D dataset, and also
          !// their transport numbers:
          ES_NT = ES_NT + 1
          ES_TRNR(ES_NT) = NTR


          !// Component scaling
          !// --------------------------------------------------------------
          !// If dataset is not in units of molecules, but e.g. kg, and
          !// tracer is different from dataset tracer, ESCALE have to
          !// reflect the conversion due to molecular masses.
          !// For molecules as input, emis_unit scales to correct
          !// molecular mass.
          call emis_unit(NTR,EUNIT,NPAR,TMASS,ENAME,AVOGNR,USCALE)

          write(6,'(a,es15.8)') '    scale input units: ', USCALE
          write(6,'(a,es15.8)') '    scale with ESCALE: ', ESCALE
          write(6,'(a,es15.8)') '    scale with EFACT:  ', EFACT

          !// Set total scaling factor for this component, i.e. combine
          !// unit scaling (USCALE), component scaling factor (ESCALE) and
          !// dataset scaling factor (EFACT):
          ES_USCALE(ES_NT) = USCALE * EFACT * ESCALE

          !// Read next tracer to use the data
          read(ifnr,*) ENAME, ESCALE

       end do !// do while (ENAME .ne. 'xxx')


       !// 4. Read emission dataset and assign to species.
       !// -----------------------------------------------------------------
       if (ES_NT .gt. 0) then
          !// Now we know:
          !//   - which dataset to read
          !//   - which species to apply it to
          !//   - scaling factors
          !//
          !// The following routine
          !//   - reads the dataset
          !//   - assigns values to the emission tables
          !//   - assigns values to diurnal/2d scaling factors
          call reademis_2d(EFILE,ETAG,ERES,ETYP,ECAT,EMON,EYEAR, ESCENYEAR, &
               ESCAL_LH, ESCAL_2D, EVERT, &
               EDATASET, IEMIS,JEMIS,MRES,MTBL, ES_NT,ES_TRNR,ES_USCALE)


          !// Print out the start and end indices for the datasets.
          !// This is only to be able to check that everything works.
          if (MTBL .gt. MTBL_OLD) &
               write(6,'(a,2i5)')'  * 2D emis idx start/end', MTBL_OLD+1,MTBL


          if (minval(E2DS(:,:,1,(MTBL_OLD+1):MTBL)) .lt. 0._r8) then
             write(6,'(a,es12.6)') f90file//':'//subr//': NEG. 2d emissions!', &
                  minval(E2DS(:,:,1,(MTBL_OLD+1):MTBL))
             write(6,'(a)') '*** CHECK EMISSIONS!!! ***'
             stop
          end if
       end if !// if (ES_NT .gt. 0) then

       !// Next dataset file name (or end2D)
       read(ifnr,*) EFILE

    end do !// do while (EFILE .ne. 'end2D')
    !// --------------------------------------------------------------------
    !// Save number of 2D tables
    NE2TBL = MTBL
    !// --------------------------------------------------------------------



    !// 3D emission tables
    !// --------------------------------------------------------------------
    !// Initialize
    NE3TBL = 0               !// no. of 2-D tables
    E3LTBL(:,:) = .false.    !// .T. = species N uses tbl E
    E3STBL(:,:) = 0._r8      !// scale factor for each N
    NM3TBL(:) = 0            !// month table
    NY3TBL(:) = 0            !// year table
    MTBL=0                   !// emission table index
    E3DSNEW(:,:,:,:) = 0._r8 !// The emissions

    !// Read header for 3D
    do J = 1,3
       read(ifnr,*) ELINE
    end do
    read(ifnr,*) EFILE !// first 3D file or 'end3D'
    !// 3D data
    !// --------------------------------------------------------------------
    do while (EFILE .ne. 'end3D')

       !// 1. First read info about the dataset.
       !// 2. Then read which species to let the dataset apply for.

       !// Info on how to read emissions
       read(ifnr,*) ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EDATASET
       write(6,'(a)') ' 3D emis file:'
       write(6,'(a)') trim(EFILE)
       write(6,'(i5,es14.6,a5,2i5,1x,a3,1x,i3,1x,2i3,1x,a20)') &
            ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EDATASET

       !// Keep track of table index befor reading new datasets
       MTBL_OLD = MTBL

       !// Find emission properties; resolution IEMIS,JEMIS
       call emis_prop(ERES,IEMIS,JEMIS)

       !// Read in tracers to use this dataset, and their scalings
       ES_NT = 0             !// Number of tracers to use the dataset
       ES_TRNR(:) = 0        !// Transport number of these tracers
       ES_USCALE(:)  = 1._r8 !// Default scaleing due to units (for each tracer)

       !// Read component name and number of user defined category scalings
       !// (May be 'xxx 0')
       read(ifnr,*) ENAME, ESCALE

       !// Check for a component using the dataset
       !// -----------------------------------------------------------------
       do while (ENAME .ne. 'xxx')
          !// Find transport number
          do NTR = 1, NPAR
             if (ENAME .eq. TNAME(NTR)) exit
          end do
          if (NTR .gt. NPAR) then
             !// Wrong species in emission list
             write(6,'(a)') '* Trying to emit a non-included species: ' &
                  //trim(ENAME)
             write(6,'(a)') '  * Skipping emission dataset!'
             !// Skip it; read next tracer to use the data
             read(ifnr,*) ENAME, ESCALE
             cycle
          end if

          !// Keep track how many species will use this 3D dataset, and also
          !// their transport numbers:
          ES_NT = ES_NT + 1
          ES_TRNR(ES_NT) = NTR

          !// Component scaling
          !// --------------------------------------------------------------
          !// If dataset is not in units of molecules, but e.g. kg, and
          !// tracer is different from dataset tracer, ESCALE have to
          !// reflect the conversion due to molecular masses.
          !// For molecules as input, emis_unit scales to correct
          !// molecular mass.
          call emis_unit(NTR,EUNIT,NPAR,TMASS,ENAME,AVOGNR,ES_USCALE(ES_NT))

          !// Also scale with ESCALE (component scaling factor) and
          !// EFACT (dataset scaling factor)
          ES_USCALE(ES_NT) = ES_USCALE(ES_NT) * EFACT * ESCALE

          write(6,'(a,es15.8)') '    scaling with ESCALE: ',ESCALE
          write(6,'(a,es15.8)') '    scaling with EFACT: ',EFACT

          !// Read tracer to use the data
          read(ifnr,*) ENAME, ESCALE

       end do !// do while (ENAME .ne. 'xxx')


       !// Check if there are included species to use dataset
       if (ES_NT .gt. 0) then
          !// Now we know:
          !//   - which dataset to read
          !//   - which species to apply it to
          !//   - scaling factors
          !//
          !// The following routine
          !//   - reads the dataset
          !//   - assigns values to the emission tables
          call reademis_3d(EFILE,ETAG,ERES,ETYP,ECAT,EMON,EYEAR,&
               EDATASET, IEMIS,JEMIS,MRES,MTBL, ES_NT,ES_TRNR,ES_USCALE)

          !// The MTBL is now updated
          write(6,'(a,2i5)')'  * 3D emis idx start/end',MTBL_OLD+1,MTBL

          if (minval(E3DSNEW(:,:,:,(MTBL_OLD+1):MTBL)).lt.0._r8) then
             write(6,'(a,es12.6)') f90file//':'//subr//': NEG. 3d emissions!', &
                  minval(E3DSNEW(:,:,:,(MTBL_OLD+1):MTBL))
             write(6,'(a)') '*** CHECK EMISSIONS!!! ***'
             stop
          end if
       end if !// if (ES_NT .gt. 0) then

       !// Next file name (or 3D end) and description/comment
       read(ifnr,*) EFILE

    end do !// do while (EFILE.ne.'end3D')
    !// --------------------------------------------------------------------
    !// Save number of 3D tables
    NE3TBL = MTBL
    !// --------------------------------------------------------------------


    !// Short term variations
    !// --------------------------------------------------------------------
    !// Read in datasets to be used for short term variation. Typically a
    !// 2D field modified by meteorology. This requires hard coding.
    !// Therefore, only one tracer can use each dataset.

    !// Initialize the known datasets
    DMSseaconc(:,:,:) = 0._r8

    !// Read header for STV
    do J = 1,2
       read(ifnr,*) ELINE
    end do
    read(ifnr,*) EFILE !// first STV file or '-END-'
    !// STV data
    !// --------------------------------------------------------------------
    do while (EFILE.ne.'-END-')

       !// Important:
       !// Short term variations can apply for only one component at
       !// a time.

       !// Important 2:
       !// However, this section can also turn on some groups of emissions,
       !// such as oceanic organic carbon aerosols. These may apply to several
       !// species, and to allow this, we also accept the dummy tracer name
       !// '---'.

       !// Info on how to read emissions
       read(ifnr,*) ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
             EUNIT, EDIURN, EDATASET
       write(6,'(a)') ' 3D STV file:'
       write(6,'(a)') trim(EFILE)
       write(6,'(i5,es14.6,a5,2i5,1x,a3,1x,i3,1x,2i3,1x,a20)') &
            ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EDATASET

       if (ERES .eq. '1x1') then
          IEMIS = 360
          JEMIS = 180
       else if (ERES.eq.'HLF') then
          !// 0.5x0.5 degrees resolution (HLF = half)
          IEMIS = 720
          JEMIS = 360
       else if (ERES.eq.'QRT') then
          !// 0.25x0.25 degrees resolution (QRT = quarter)
          IEMIS = 1440
          JEMIS = 720
       else if (ERES .eq. 'PPP') then
          !// Point sources, reading is done elsewhere, but IEMIS/JEMIS needs
          !// values; use 360/180.
          IEMIS = 360
          JEMIS = 180
       else
          write(6,'(a)') f90file//':'//subr// &
               ': Unknown resolution: '//ERES
          stop
       end if

       read(ifnr,*) ENAME

       !// Check dataset vs name
       !if (trim(EDATASET) .eq. 'DMSseaconc' .and. .not.(ENAME.eq.'DMS')) then
       !   write(6,'(a)') '* STV: Should not use DMSseaconc with other' &
       !        //' species than DMS!'
       !   stop
       !end if

       if (trim(ENAME) .eq. '---') then
          NTR = -1
       else
          !// Find transport number
          do NTR = 1, NPAR
             if (ENAME .eq. TNAME(NTR)) exit
          end do
       end if
       if (NTR .gt. NPAR) then
          !// Wrong species in emission list
          write(6,'(a)') '* STV: Trying to emit a non-included species: ' &
               //trim(ENAME)
          write(6,'(a)') '  * Skipping emission dataset!'
       else
          !// Read the STV data, interpolate and put them in their arrays.
          !// Will also be carried out if NTR == -1.
          call reademis_stv(EFILE,ETAG,EFACT,ERES,EYEAR,ETYP,ECAT,&
               IEMIS,JEMIS, MRES,EDATASET)
       end if

       !// Next file name (or end) and description/comment
       read(ifnr,*) EFILE
    end do !// do while (EFILE.ne.'-END-')


    !// Print number of emission tables
    write(6,'(a,i6)') '** Number of 2D emission tables:',NE2TBL
    write(6,'(a,i6)') '** Number of 3D emission tables:',NE3TBL



    !// Forest fires - NOT YET READY
    !// --------------------------------------------------------------------
    !// Initialize
    FF_TYPE = 0
    EMIS_FIR(:,:,:,:,:) = 0._r8
    ECOMP_FIR(:) = 0


    !// Forest fires
    !// --------------------------------------------------------------------
    !// Read FF header
    read(ifnr,*) ELINE
    if (trim(ELINE) .ne. '-BBB-') then
       print*,'No BBB in emission list'
       stop
    end if
    do J = 1, 2
       read(ifnr,*) ELINE
    end do
    read(ifnr,*) EFILE !// Read FF path

    if (trim(EFILE) .ne. 'endFF') then

       FF_PATH = trim(EFILE)

       !// Info on how to read 
       read(ifnr,*) ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EDATASET !, ESCENYEAR
       write(6,'(a)') ' Forest fires path: '
       write(6,'(a)') trim(FF_PATH)
       !write(6,'(i5,es14.6,a5,2i5,1x,a3,1x,i3,1x,3i3,1x,a40,1x,i4)') &
       write(6,'(i5,es14.6,a5,2i5,1x,a3,1x,i3,1x,2i3,1x,a20)') &
            ETAG, EFACT, ERES, EMON, EYEAR, ECAT, ETYP, &
            EUNIT, EDIURN, EDATASET !, ESCENYEAR

       if (ETAG .eq. 567 .or. ETAG .eq. 569 .or. ETAG .eq. 570 .or. &
            ETAG .eq. 571 .or. ETAG .eq. 572 .or. ETAG .eq. 573) then
          !// GFEDv4
          if (ETAG .eq. 567) then
             FF_TYPE = 4 !// Old version of vertical distribution
          else if (ETAG .eq. 573) then
             FF_TYPE = 3 !// Old version of vertical distribution,
                         !// and use daily fractions.
          else if (ETAG .eq. 569) then
             FF_TYPE = 6 !// Distribute according to air mass in BLH
          else if (ETAG .eq. 570) then
             FF_TYPE = 7 !// Distribute according to air mass in BLH
                         !// and use diurnal fractions.
          else if (ETAG .eq. 571) then
             FF_TYPE = 8 !// Distribute according to altitude in BLH
          else if (ETAG .eq. 572) then
             FF_TYPE = 9 !// Distribute according to altitude in BLH
                         !// and use diurnal fractions.
          end if
          !// Year to read (9999 will use meteorological year)
          FF_YEAR = EYEAR

          !// Read info line
          read(ifnr,*) ELINE

          !// Read all that components that should be read
          NEFIR = 0 !// Counter
          do N = 1, EPAR_FIR
             read(ifnr,*) ENAME, BNAME, EFACT, partitions

             if (trim(ENAME) .ne. '---') then
                !// Find transport number
                do NTR = 1, NPAR
                   if (trim(ENAME) .eq. trim(TNAME(NTR))) exit
                end do
             else
                exit !// Done reading
             end if

             if (NTR .gt. NPAR) then
                !// Wrong species in emission list
                write(6,'(a)') '* '//trim(ENAME)//' not included - skipping'
             else if (NTR .gt. 0) then
                NEFIR = NEFIR + 1
                if (NEFIR .gt. EPAR_FIR) then
                   write(6,'(a,2i5)') f90file//':'//subr// &
                        ': GFEDv4: Number of species > EPAR_FIR',NEFIR,EPAR_FIR
                   stop
                end if

                !// Save transport numbers
                ECOMP_FIR(NEFIR) = NTR

                FF_CNAMES(NEFIR) = ENAME
                FF_BNAMES(NEFIR) = BNAME
                FF_SCALE(NEFIR) = EFACT
                FF_PARTITIONS(:,NEFIR) = partitions(:)
                write(6,'(a10,a20,f7.4,6f6.2)') &
                     FF_CNAMES(NEFIR), FF_BNAMES(NEFIR), FF_SCALE(NEFIR), &
                     FF_PARTITIONS(:,NEFIR)

             end if

          end do !// do N = 1, EPAR_FIR
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Number of emitted species (NEFIR): ',NEFIR

          !// Also initialise GFED4 emission factors
          call gfed4_init()

       else if (ETAG .eq. 568 .or. ETAG .eq. 574) then
          !// CEDS
          if (ETAG .eq. 568) then
             FF_TYPE = 5  !// CEDS distributed with RETRO vertical
          else
             FF_TYPE = 51 !// CEDS distributed in BLH (altitude-weighted)
          end if

          !// Year to read (9999 will use meteorological year)
          FF_YEAR = EYEAR

          !// Read info line
          read(ifnr,*) ELINE

          !// Read all that components that should be read
          NEFIR = 0 !// Counter
          do N = 1, EPAR_FIR
             read(ifnr,*) ENAME, BNAME, EVAR, EFACT

             if (trim(ENAME) .ne. '---') then
                !// Find transport number
                do NTR = 1, NPAR
                   if (trim(ENAME) .eq. trim(TNAME(NTR))) exit
                end do
             else
                exit !// Done reading
             end if

             if (NTR .gt. NPAR) then
                !// Wrong species in emission list
                write(6,'(a)') '* '//trim(ENAME)//' not included - skipping'
             else if (NTR .gt. 0) then
                NEFIR = NEFIR + 1
                if (NEFIR .gt. EPAR_FIR) then
                   write(6,'(a,2i5)') f90file//':'//subr// &
                        ': CEDS: Number of species > EPAR_FIR',NEFIR,EPAR_FIR
                   stop
                end if

                !// Save transport numbers
                ECOMP_FIR(NEFIR) = NTR

                FF_CNAMES(NEFIR) = ENAME
                FF_BNAMES(NEFIR) = BNAME ! change to FF_PREFIX
                FF_VARNAME(NEFIR) = EVAR
                FF_SCALE(NEFIR) = EFACT
                write(6,'(a10,a20,a20,f7.4)') &
                     FF_CNAMES(NEFIR), FF_BNAMES(NEFIR), &
                     FF_VARNAME(NEFIR), FF_SCALE(NEFIR)

             end if

          end do !// do N = 1, EPAR_FIR
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Number of emitted species (NEFIR): ',NEFIR

       end if
    end if !// if (trim(EFILE) .ne. 'endFF') then

    !// Close file
    close(ifnr)

    !// --------------------------------------------------------------------
  end subroutine emis_input
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_emis(NDAY,NMET,NOPS,NDAYI,LNEWM,DT)
    !// --------------------------------------------------------------------
    !// Update short term emissions, using separate emission arrays (i.e.
    !// not E2DS or E3DS).
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LBCOC
    use cmn_ctm, only: JDAY, JMON, JDATE, JYEAR
    use emisutils_oslo, only: gfed4_rd, gfed4_rd_daily, ceds_biomass_burning, &
         gfed4_rd_novert, gfed4_rd_novert_daily, &
         emis_setscaling_2dfields, ceds_biomass_burning_novert
    use cmn_oslo, only: FF_TYPE
    use emissions_aircraft, only: aircraft_emis_master
    use emissions_megan, only: megan_update_metdata
    use emissions_ocean, only: emissions_ocean_getChlA
    use emissions_volcanoes, only: read_volcEMIS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS, NDAYI
    logical, intent(in) :: LNEWM
    real(r8), intent(in)  :: DT
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'update_emis'
    !// --------------------------------------------------------------------

    !// MEGAN emissions - only each NMET
    if (NOPS .eq. 1) call megan_update_metdata(NDAY, NMET, NDAYI)

    !// Forest fires - monthly data, possibly diurnal data
    if (FF_TYPE .eq. 3) then
       !// GFEDv4 daily, old vertical
       if (NOPS .eq. 1 .and. NMET .eq. 1) &
            call gfed4_rd_daily(JDATE,JMON,JYEAR)
    else if (FF_TYPE .eq. 4) then
       !// GFEDv4 monthly, old vertical - see also FF_TYPE 6
       if (LNEWM) call gfed4_rd(JMON,JYEAR)
    else if (FF_TYPE .eq. 5) then
       !// CEDS
       if (LNEWM) call ceds_biomass_burning(JMON,JYEAR)
    else if (FF_TYPE .eq. 51) then
       !// CEDS distributed in BLH
       if (LNEWM) call ceds_biomass_burning_novert(JMON,JYEAR)
    else if (FF_TYPE .eq. 6) then
       !// GFEDv4 no vertical scaling in read-in, use BLH p for that.
       if (LNEWM) call gfed4_rd_novert(JMON,JYEAR)
    else if (FF_TYPE .eq. 7) then
       !// GFEDv4, daily, no vertical scaling in read-in, use BLH p for that.
       if (NOPS .eq. 1 .and. NMET .eq. 1) &
            call gfed4_rd_novert_daily(JDATE,JMON,JYEAR)
    else if (FF_TYPE .eq. 8) then
       !// GFEDv4 no vertical scaling in read-in, use BLH Z for that.
       if (LNEWM) call gfed4_rd_novert(JMON,JYEAR)
    else if (FF_TYPE .eq. 9) then
       !// GFEDv4, daily, no vertical scaling in read-in, use BLH Z for that.
       if (NOPS .eq. 1 .and. NMET .eq. 1) &
            call gfed4_rd_novert_daily(JDATE,JMON,JYEAR)
    else if (FF_TYPE .ne. 0) then
       write(6,'(a,i2,a)') f90file//':'//subr//': FF_TYPE = ',FF_TYPE,&
            ' is not valid:'
       stop
    end if


    !// Aircraft emissions
    call aircraft_emis_master(NMET,NOPS,LNEWM,DT)

    !// Set 2D scalings for e.g. isoprene and monoterpenes
    if (LNEWM) call emis_setscaling_2dfields()

    !// Volcanoe SO2 emissions to be read at start and every year
    if ((JDAY.eq.1 .and. JMON.eq.1) .or. (NDAY.eq.NDAYI)) then
       if (NMET.eq.1 .and. NOPS.eq.1) then
          call read_volcEMIS()
       end if
    end if

    !// Update chlorophyll A
    if (LBCOC .and. LNEWM) call emissions_ocean_getChlA(LNEWM)

    !// --------------------------------------------------------------------
  end subroutine update_emis
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_emis_ij(NMET,NOPS,NSUB,CCYC,MP)
    !// --------------------------------------------------------------------
    !// Update short term emissions, using separate emission arrays (i.e.
    !// not E2DS or E3DS).
    !//
    !// Ole Amund Sovde, August 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: LSALT, LBCOC
    use cmn_precision, only: r8
    use emissions_ocean, only: emissions_ocean_organiccarbon
    use seasalt, only: seasalt_emis
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NMET, NOPS, NSUB, CCYC, MP
    !// --------------------------------------------------------------------

    if (LSALT) call seasalt_emis(NMET, NOPS, NSUB, CCYC, MP)

    if (LBCOC) call emissions_ocean_organiccarbon(NMET, NOPS, NSUB, CCYC, MP)
    
    !// --------------------------------------------------------------------
  end subroutine update_emis_ij
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine emis_prop(ERES,IEMIS,JEMIS)
    !// --------------------------------------------------------------------
    !// Sets the following properties of emission data:
    !//   - Resolution (IEMIS,JEMIS) from ERES
    !//
    !// Ole Amund Sovde, March 2012, August 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: ERES
    !// Output
    integer, intent(out)     :: IEMIS,JEMIS
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_prop'
    !// --------------------------------------------------------------------

    !// What is the dimension
    if (ERES.eq.'1x1') then
       !// 1x1 degrees resolution
       IEMIS = 360
       JEMIS = 180
    else if (ERES.eq.'HLF') then
       !// 0.5x0.5 degrees resolution (HLF = half)
       IEMIS = 720
       JEMIS = 360
    else if (ERES.eq.'QRT') then
       !// 0.25x0.25 degrees resolution (QRT = quarter)
       IEMIS = 1440
       JEMIS = 720
    else if (ERES.eq.'ZP1') then
       !// 0.1x0.1 degrees resolution (ZP1 = zero point 1)
       IEMIS = 3600
       JEMIS = 1800
    else if (ERES.eq.'x22') then
       !// T159N80 in 2x2 degraded resolution
       IEMIS = 160
       JEMIS = 80
    else
       write(6,'(a)') f90file//':'//subr// &
            ': Unknown resoltuion: '//ERES
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine emis_prop
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emis_unit(NTR,EUNIT,NPAR,TMASS,ENAME,AVOGNR,USCALE)
    !// --------------------------------------------------------------------
    !// Sets the scaling factor based on unit flag in emission list.
    !// - Scales from input units to kg/s
    !// - Not accounted for automatically:
    !//   Difference in molecular weights of dataset tracer
    !//   vs tracer for which the dataset is applied.
    !// - Area scalings (e.g. 1/cm2 to 1/m2 are treated elsewhere by ETYP.
    !//
    !// More unit flags can be added when needed. So far we have:
    !// 0: No scaling; input field is kg/s or kg/area/s
    !// 1: Input data is molec/s or molec/area/s
    !// 2: Input data is kg/y or kg/area/y
    !// 3: Input data is kg/month or kg/area/month
    !// 
    !//
    !// Ole Amund Sovde, March 2012, August - September 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)         :: NTR, EUNIT, NPAR
    real(r8), intent(in)        :: TMASS(NPAR),AVOGNR
    character(len=*),intent(in) :: ENAME
    !// Output
    real(r8), intent(out)       :: USCALE

    !// Locals
    real(r8), parameter :: SECYR  = 31536000._r8 !// 365-day
    real(r8), parameter :: SECMON = 2628000._r8  !// SECYR/12
    character(len=50) :: cscaleinfo
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_unit'
    !// --------------------------------------------------------------------

    !// IMPORTANT!!!
    !// This scaling assumes dataset contains same tracer as it applies for.
    !// In case they differ, a scaling based on molecular weights must
    !// be applied. The calling routine (emis_input) explains this.

    !// Initialize
    USCALE = 1._r8
    if (EUNIT .eq. 0) then
       !// kg/s
       USCALE = 1._r8
       cscaleinfo = ', indata is kg/s, no scaling necessary: '
    else if (EUNIT .eq. 1) then
       !// molec/s (ETYP takes care of 1/area)
       !// mass = molec /AVOGNR*Mw*1.d-3kg/g
       !// Uses the tracer mass of species NTR, not for dataset species.
       !// Scaling due to differences in molecular masses must be done
       !// in the SCALING factor in the emission file.
       USCALE = 1.e-3_r8 * TMASS(NTR) / AVOGNR
       cscaleinfo = ', from molec/s to kg/s, for '//trim(ENAME)//':'
    else if (EUNIT .eq. 2) then
       !// kg/y (ETYP takes care of 1/area)
       !// 1/y to 1/s: 1/(3600*24*365=31536000)
       USCALE = 1._r8 / SECYR
       cscaleinfo = ', from kg/year to kg/s: '
    else if (EUNIT .eq. 3) then
       !// kg/month (ETYP takes care of 1/area)
       !// 1/month to 1/s: 1/(3600*24*365/12=2628000)
       USCALE = 1._r8 / SECMON
       cscaleinfo = ', from kg/month to kg/s: '
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': No unit specified:',EUNIT
       stop
    end if


    !// Does not scale kg of input with tracer molecular mass. Use
    !// ESCALE for that.
    write(6,'(a,i2,a)') '  * EUNIT = ',EUNIT,trim(cscaleinfo)

    !// --------------------------------------------------------------------
  end subroutine emis_unit
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emis_define_2dscaling(EDIURN,E2LocHour,E22dField)
    !// --------------------------------------------------------------------
    !// Check diurnal scaling flag and define whether LocHour or 2d field
    !// is used for scaling purposes.
    !//
    !// Ole Amund Sovde, November 2015
    !// --------------------------------------------------------------------
    use cmn_oslo, only: E2LocHourRETRO, E2LocHour5050, &
         E22dSunTemp, E22dSun, E22dHDD, NE22dVARS, E2L2dTBL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: EDIURN
    integer, intent(out) :: E2LocHour, E22dField
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'emis_define_2dscaling'
    !// --------------------------------------------------------------------

    !// Return if no scaling is set
    if (EDIURN .eq. 0) then
       E2LocHour = 0
       E22dField = 0
       return
    end if

    if (EDIURN .eq. 1) then
       E2LocHour = E2LocHourRETRO
       E22dField = 0
    else if (EDIURN .eq. 2) then
       E2LocHour = E2LocHour5050
       E22dField = 0
    else if (EDIURN .eq. 3) then
       E2LocHour = 0
       E22dField = E22dSunTemp
       E2L2dTBL(E22dSunTemp) = .true. !// Keep track of 2d dataset used
    else if (EDIURN .eq. 4) then
       E2LocHour = 0
       E22dField = E22dSun
       E2L2dTBL(E22dSun) = .true.     !// Keep track of 2d dataset used
    else if (EDIURN .eq. 5) then
       E2LocHour = 0
       E22dField = E22dHDD
       E2L2dTBL(E22dHDD) = .true.     !// Keep track of 2d dataset used
    else
       write(6,'(a, i5)') f90file//':'//subr// &
            ': Not defined EDIURN: ',EDIURN
       stop 'in '//subr
    end if

    !// --------------------------------------------------------------------
  end subroutine emis_define_2dscaling
  !// ----------------------------------------------------------------------


  !// ------------------------------------------------------------------
end module emissions_oslo
!//=========================================================================

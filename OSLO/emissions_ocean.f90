!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, June 2015
!//=========================================================================
!// Oceanic emissions.
!//=========================================================================
module emissions_ocean
  !// ----------------------------------------------------------------------
  !// MODULE: emissions_ocean
  !// DESCRIPTION: Contains utilities to calculate oceanic emissions,
  !//              mainly primary organic aerosols (POA), also called
  !//              organic carbon.
  !//
  !// Contains:
  !//
  !// Ole Amund Sovde, June 2015
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IDBLK, JDBLK, MPBLK
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// In case where sea salt module is turned off, these emissions needs
  !// to define their own sea salt flux.
  !// To make it possible to run organic aerosol with sea salt module,
  !// N_SS_BINS_STANDALONE must be as large as the defined NPAR_SALT.
  integer, parameter :: N_SS_BINS_STANDALONE = 8
  !// Lowest/highest diam dry bin limit (m)
  real(r8), parameter  :: Dsaltdrymin = 0.03e-6_r8
  real(r8), parameter  :: Dsaltdrymax = 25.0e-6_r8
  !// Density [kg/m3]
  real(r8), parameter  :: rhosaltdry  = 2200._r8

  real(r8), dimension(N_SS_BINS_STANDALONE+1) :: &
       rsalt80, & !// Radius at 80%RH in um
       Dsaltdry   !// Dry diameters of bins (m)

  !// Scheme for sea salt production for STANDALONE parameterisation
  !// 1: Monahan et al (Oceanic Whitecaps, 1986, Reidel)
  !//    and Smith (QJRMS, 1993, doi:10.1002/qj.49711951211)
  !// 2: Mårtensson et al (JGR, 2003, doi:10.1029/2002JD002263)
  !// 3: Gantt et al (GMD, 2015, doi:10.5194/gmd-8-619-2015)
  !// 4: Witek et al (JGR, 2016, doi:10.1002/2015JD023726)
  integer, parameter :: SeaSaltScheme = 4

  real(r8) :: scaleEPOA !// Scale to match EPOA to observed values.

  integer :: N_SS_BINS

  !// Flag to include these emissions (set T if correct tracer name exists)
  logical :: LOCEAN_OCARB
  !// Year for which chlorophyll is used
  integer :: YEAR_CHLA
  !// Path to chlorophyll a data
  character(len=80) :: chlaPATH

  !// Emission array: Emissions of OC kg/s.
  real(r8), dimension(IDBLK, JDBLK, MPBLK) :: EPOAseaspray

  !// Chlorophyll A (mg/m3)
  real(r8), dimension(IDBLK, JDBLK, MPBLK) :: chlorophyllA

  !// Emitted primary organic carbon (POA) from ocean
  integer, parameter :: NPAR_ocnPOA = 2
  character(len=10), dimension(NPAR_ocnPOA), parameter :: &
       ocnPOAnames = (/ 'omOCNfob', 'omOCNfil' /)
  !// Fraction of emissions (hydrophil / hydrophob)
  !// According to Gantt et al (2012, ACP) emissions should be hydrophobic,
  !// following Decesari et al (2007, Env.Sci.Tech.) and
  !// Facchini et al (2008, Geophys.Res.Lett.).
  !// We use 80/20 for now.
  real(r8), dimension(NPAR_ocnPOA), parameter :: &
       ocnPOA_emisFraction = (/ 0.8_r8, 0.2_r8 /)
  !// Transport numbers
  integer, dimension(NPAR_ocnPOA) :: ocnPOA_trsp_idx
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'emissions_ocean.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public EPOAseaspray, emissions_ocean_organiccarbon_setpath, &
       emissions_ocean_organiccarbon_init, emissions_ocean_getChlA, &
       emissions_ocean_organiccarbon, emissions_ocean_total, &
       add_oceanOCemis
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine emissions_ocean_organiccarbon_setpath(inpath,inyear)
    !// --------------------------------------------------------------------
    !// Sets path and year for chlorophyll dataset.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LSALT, NPAR_SALT, NPAR
    use cmn_chem, only: TNAME
    use seasalt, only: get_rsalt80um, get_seasaltscheme
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: inpath
    integer, intent(in) :: inyear

    !// Locals
    real(r8) :: &
         logbin1, &    !// Logarithm of lowest bin length (log(m))
         logbinmax, &  !// Logarithm of highest bin length (log(m))
         totlogstep, & !// Logarithm of whole bin length
         logstep, &    !// Logarithm of one bin length
         logbin        !// Logarithm of arbitrary bin length
    integer :: N,NN
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='emissions_ocean_organiccarbon_setpath'
    !// ------------------------------------------------------------------

    !// Initialise emission array
    EPOAseaspray(:,:,:) = 0._r8
    !// Transport numbers of oceanic carbon species
    ocnPOA_trsp_idx(:) = -1

    !// Should include oceanic carbon emissions if species are
    !// included.
    do NN = 1, NPAR_ocnPOA
       do N = 1, NPAR
          if (trim(TNAME(N)) .eq. ocnPOAnames(NN)) ocnPOA_trsp_idx(NN) = N
       end do
    end do

    if (minval(ocnPOA_trsp_idx) .lt. 1) then
       write(6,'(a)') f90file//':'//subr// &
            ': Not all ocnPOA species'// &
            ' are included: skipping oceanic emissions.'
       LOCEAN_OCARB = .false.
       return
    end if

    LOCEAN_OCARB = .true.
    !// Year to read chlorophyll A data for
    YEAR_CHLA = INYEAR
    !/ Path to chlorophyll A data
    chlaPATH = inpath
    !// Make sure path ends with /
    n = len(trim(chlaPATH))
    if (chlaPATH(n:n) .ne. '/') chlaPATH = trim(chlaPATH)//'/'
    
    !// --------------------------------------------------------------------
  end subroutine emissions_ocean_organiccarbon_setpath
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emissions_ocean_organiccarbon_init()
    !// --------------------------------------------------------------------
    !// Initialise variables for using oceanic emissions.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LSALT, NPAR_SALT, NPAR
    use cmn_chem, only: TNAME
    use seasalt, only: get_rsalt80um, get_seasaltscheme
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    real(r8) :: &
         logbin1, &    !// Logarithm of lowest bin length (log(m))
         logbinmax, &  !// Logarithm of highest bin length (log(m))
         totlogstep, & !// Logarithm of whole bin length
         logstep, &    !// Logarithm of one bin length
         logbin        !// Logarithm of arbitrary bin length
    integer :: N
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='emissions_ocean_organiccarbon_init'
    !// ------------------------------------------------------------------

    if (.not. LOCEAN_OCARB) return
    
    !// Sea salt data needed for calculating oceanic organic carbon
    if (LSALT) then
       !// Check if SeaSaltScheme matches seasalt.f90
       call get_seasaltscheme(N)
       if (N .ne. SeaSaltScheme) then
          write(6,'(a)') f90file//':'//subr// &
               ': SeaSaltScheme does not match SeaSaltScheme in SALT routine'
          stop
       end if

       !// Use data from sea salt module
       N_SS_BINS = NPAR_SALT
       !// Check for array size
       if (NPAR_SALT .gt. N_SS_BINS_STANDALONE) then
          write(6,'(a)') f90file//':'//subr// &
               ': NPAR_SALT is larger '// &
               'than N_SS_BINS_STANDALONE: will not work'
          stop
       end if
       !// Get rsalt80_ssm(1:NPAR_SALT) from sea salt module
       rsalt80(:) = 0._r8
       call get_rsalt80um(rsalt80,N_SS_BINS)
 
    else
       !// Stand-alone version: Must to calculate sea salt itself
       N_SS_BINS = N_SS_BINS_STANDALONE
       !// Calculate these based on N_SS_BINS, as for sea salt code.
       !// Use r80 as in sea salt code to generate rsalt80(:N_SS_BINS)
       !// Set salt lowest and highest diameter
       Dsaltdry(1) = Dsaltdrymin
       Dsaltdry(N_SS_BINS+1) = Dsaltdrymax
       !// logarithmically even-spaced Diameters
       logbin1   = log(Dsaltdry(1))
       logbinmax = log(Dsaltdry(N_SS_BINS+1))
       totlogstep= logbinmax - logbin1
       logstep   = totlogstep / real(N_SS_BINS, r8)
       logbin = logbin1
       do N = 2, N_SS_BINS + 1
          logbin = logbin + logstep
          Dsaltdry(N) = exp(logbin)
       end do
       !// Diameters at 80% RH is 2x dry diameter (Fitzgerald 1975)
       !// and radius at 80% RH is half of diameter, so in [um]:
       rsalt80(:)= 1.e6_r8 * Dsaltdry(:)
    end if

    !// Check rsalt
    if (maxval(rsalt80) .le. 0._r8) then
       write(6,'(a)') f90file//':'//subr// &
            ': rsalt80 <= 0: Check if it is initialised'
       if (LSALT) write(6,'(a)') 'You are retrieving it from SALT'
       stop
    end if

    !// Set scaleEPOA. Should be calculated from climatological averages.
    !// Note that scaling will depend slightly on resolution.
    if (SeaSaltScheme .eq. 1) then
       scaleEPOA = 0.35_r8  !// First guess
    else if (SeaSaltScheme .eq. 2) then
       scaleEPOA = 0.35_r8  !// First guess
    else if (SeaSaltScheme .eq. 3) then
       scaleEPOA = 0.5_r8 !// About right
    else if (SeaSaltScheme .eq. 4) then
       scaleEPOA = 0.57_r8 !// 6.3Tg for mean of 2000 and 2016
    else
       scaleEPOA = 1._r8  !// Not known yet; assume 1.
    end if

    !// --------------------------------------------------------------------
  end subroutine emissions_ocean_organiccarbon_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emissions_ocean_getChlA(LNEWM)
    !// --------------------------------------------------------------------
    !// Fetch chlorophyll A for this month.
    !// Uses MODIS chlorophyll A, which has been binned into 0.5x0.5 degree
    !// grid. Uses meteorological year for years 2003-2012, and otherwise
    !// a climatology of those years.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LBCOC, IPAR, JPAR, MPBLK
    use cmn_ctm, only: JMON, XDEDG, YDEDG, AREAXY, MPBLKIB, MPBLKIE, &
         MPBLKJB, MPBLKJE
    use cmn_met, only: MYEAR
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    logical, intent(in) :: LNEWM

    !// Locals
    real(r8) :: R8CTM(IPAR,JPAR), XYBOX
    integer           :: nLon, nLat, I,J, II,JJ, MP, getY
    character(len=7)  :: datestamp
    character(len=2)  :: CMON
    character(len=120) :: infile
    real(r8), allocatable, dimension(:) :: inLon, inLat, inLonE, inLatE
    real(r8), allocatable, dimension(:,:) :: R8XY
    !// --------------------------------------------------------------------

    !// Only update when new month
    if (.not. LNEWM) return

    !// Only read data if oceanic emissions are included
    if (.not. LOCEAN_OCARB) return

    !// Filename to read
    if (YEAR_CHLA .eq. 9999) then
       if (MYEAR .lt. 2003 .or. MYEAR .gt. 2014) then
          !// Use a climatology
          write(CMON(1:2),'(i2.2)') JMON
          infile = trim(chlaPATH)//'chlorophyllA/modis_chlorophyll_a_'// &
               'climatology_2003_2014_'//CMON//'.nc'
       else
          !// Use meteorological year
          write(datestamp(1:7),'(i4.4,a1,i2.2)') MYEAR,'_',JMON
          infile = trim(chlaPATH)//'chlorophyllA/modis_chlorophyll_a_'// &
               datestamp//'.nc'
       end if
    else if (YEAR_CHLA .gt. 0) then
       !// Use meteorological year
       write(datestamp(1:7),'(i4.4,a1,i2.2)') YEAR_CHLA,'_',JMON
       infile = trim(chlaPATH)//'chlorophyllA/modis_chlorophyll_a_'// &
            datestamp//'.nc'
    else
       !// Use the climatology
       write(CMON(1:2),'(i2.2)') JMON
       infile = trim(chlaPATH)//'chlorophyllA/modis_chlorophyll_a_'// &
            'climatology_2003_2014_'//CMON//'.nc'
    end if

    write(6,'(a)') '* Reading '//trim(infile)

    !// Check resolution (latitude/longitude)
    !// This routine allocates inLon/inLat
    call get_netcdf_var_1d( infile, 'lon',  inLon  )
    call get_netcdf_var_1d( infile, 'lat',  inLat  )
    call get_netcdf_var_1d( infile, 'lonE',  inLonE  )
    call get_netcdf_var_1d( infile, 'latE',  inLatE  )

    nLon  = SIZE( inLon  )
    nLat  = SIZE( inLat  )

    allocate( R8XY(nLon,nLat) )

    !// Get data, should be resolution 720x360
    call get_netcdf_var_2d(infile, 'chlor_a', R8XY, nlon, nlat)


    !// Regrid: must multiply by area
    do J = 1, nLat
       !// Area of emission grid box
       XYBOX =  A0*A0 * CPI180*(inLonE(2) - inLonE(1)) &
            * (sin(CPI180*inLatE(J+1)) - sin(CPI180*inLatE(J)))
       !// Convert to mg/m
       R8XY(:,J) = R8XY(:,J) * XYBOX
    end do
    call E_GRID(R8XY(:,:), inLonE, inLatE, nLon, nLat, &
         R8CTM, XDEDG, YDEDG, IPAR, JPAR, 1)
    !// Back to mg/m3
    R8CTM(:,:) = R8CTM(:,:) / AREAXY(:,:)

    !// Deallocate all local variables
    deallocate( inLon, inLat, inLonE, inLatE, R8XY )


    !// Convert to MP-block style
    do MP = 1, MPBLK
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             chlorophyllA(II,JJ,MP) = R8CTM(I,J)
          end do
       end do
    end do


    write(6,'(a,2es12.3)') '* Updated MODIS chlorophyll A ', &
         minval(R8CTM),maxval(R8CTM)

    !// --------------------------------------------------------------------
  end subroutine emissions_ocean_getChlA
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emissions_ocean_organiccarbon(NMET, NOPS, NSUB, CCYC, MP)
    !// --------------------------------------------------------------------
    !// Calculate emissions of oceanic organic carbon aerosols (POA).
    !// If sea salt module is included, use sea salt production rate,
    !// otherwise calculate it separately (slightly different method
    !// for now).
    !//
    !// Method after Gantt et al (2015, GMD, doi:10.5194/gmd-8-619-2015).
    !//
    !// This routine only depend on meteorological variables, and is called
    !// from update_emis_ii (emissions_oslo.f90), which again is called
    !// from pmain.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR_SALT, LSALT, LBCOC
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, PLAND, AREAXY
    use cmn_met, only: UMS, VMS, SFU, SFV, SFT, CI
    use cmn_sfc, only: LSMASK
    use cmn_parameters, only: CPI, TK_0C, LDEBUG
    use seasalt, only: seasalt_getflux
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
         fsaltwater, &
         wind10, &
         OMssa, chla, &
         flux, &
         sfcT, &
         rdry, Vpartdry, saltpartweight, AXY, dbin, &
         rbinStart, rbinEnd, rsubbin, dr_subbin, &
         A, B, A1, A2, r01, r02,f1,f2, &
         Epoa, rtmp
    real(r8) :: ssProd(N_SS_BINS_STANDALONE)
    integer :: I, J, II, JJ, N, NN, UPDATEEMIS
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='emissions_ocean_organiccarbon'
    !// ------------------------------------------------------------------

    !// To be emitted only at surface layer

    !// Only update if meteorology is updated.
    UPDATEEMIS = NOPS * NSUB * CCYC
    if (UPDATEEMIS .ne. 1) return

    !// Only calculate data if oceanic emissions are included
    if (.not.(LBCOC .and. LOCEAN_OCARB)) then
       !// No POA emissions
       EPOAseaspray(:,:,MP) = 0._r8
       return
    end if

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

          !// Skip calculations if no sea water
          if (fsaltwater .le. 0._r8) then
             EPOAseaspray(II,JJ,MP) = 0._r8
             !// Go to next (I,J)
             cycle
          end if


          !// 10m wind [m/s]
          wind10 = sqrt(SFU(I,J)**2 + SFV(I,J)**2)


          !// Chlorophyll A [mg/m3]
          chla = chlorophyllA(II,JJ,MP)

          !// Area of grid box
          AXY = AREAXY(I,J)

          !// Sea salt flux [kg/s]
          if (LSALT) then

             !// Get seasalt_flux(1:N_SS_BINS,II,JJ,MP)
             call seasalt_getflux(ssProd,N_SS_BINS,II,JJ,MP)

          else


             !// Calculate these based on N_SS_BINS, as for sea salt code.
             !// Use r80 as in sea salt code to generate rsalt(:N_SS_BINS)

             !// Calculate ssprod(1:N_SS_BINS) as in Gantt et al (2015)
             ssProd(:) = 0._r8

             do N = 1, N_SS_BINS
                !// Set rbin = r80 (r at RH 80%) [um]
                !// This is the mean radius of bin.
                rbinStart = rsalt80(N)
                rbinEnd   = rsalt80(N+1)

                if (SeaSaltScheme .eq. 1) then
                   !// Standard Oslo CTM3 method
                   call seasalt_production(wind10,fsaltwater,AXY, &
                        rbinStart,rbinEnd, rhosaltdry, I,J, flux)
                else if (SeaSaltScheme .eq. 2) then
                   !// Mårtensson et al 2003
                   call seasalt_production_martensson03(wind10, fsaltwater, &
                        AXY, rbinStart, rbinEnd, rhosaltdry, SFT(I,J), I, J, &
                        flux)
                else if (SeaSaltScheme .eq. 3) then
                   !// Gantt et al 2015.
                   call seasalt_production_gantt15(wind10,fsaltwater,AXY, &
                        rbinStart,rbinEnd, rhosaltdry, SFT(I,J), I,J, flux)
                else if (SeaSaltScheme .eq. 4) then
                   !// Witek et al 2016.
                   call seasalt_production_witek16(wind10,fsaltwater,AXY, &
                        rbinStart,rbinEnd, rhosaltdry, SFT(I,J), I,J, flux)
                else
                   !// Not defined
                   write(6,'(a,i3)') f90file//':'//subr// &
                        ': production scheme not defined', SeaSaltScheme
                   stop
                end if

                ssProd(N) = ssProd(N) + flux
             end do !// do N = 1, N_SS_BINS



          end if


          !// OMssa term not dependent on particle size
          rtmp = 1._r8 + exp(3._r8 * (-2.63_r8 * chla + 0.18_r8 * wind10))

          !// Loop through the bins to calculate emission (Epoa) of primary
          !// organic aerosols (POA):
          Epoa = 0._r8
          do N = 1, N_SS_BINS
             !// Calculate OMssa for each bin.
             !// Skip further division of bins into sub-bins (which is done
             !// in sea salt routine).

             !// Set rbin = r80 (r at RH 80%) [um]
             !// The mean radius of bin: 0.5_r8 * (rsalt80(N) + rsalt80(N+1))
             !// Mean diameter is then:
             dbin = (rsalt80(N) + rsalt80(N+1))

             !// Organic mass fraction of sea spray aerosols [0:1]
             OMssa = (1._r8 / rtmp) / (1._r8 + 0.03_r8 * exp(6.81_r8 * dbin)) &
                     + 0.03_r8 / rtmp

             !// Calculate organic carbon mass production [kg/s]
             !// -----------------------------------------------------------
             !// Fraction of open salt water has already been taken into
             !// account in ssProd.
             !//
             !// Should also take into account that density of organic
             !// matter is about 1g/cm3 while sea salt is 2.2g/cm3.
             !// A simpler method is to calculate the total mass of the
             !// particle, assumed it consists of sea salt and organic matter.
             !// Since ssProd is flux [kg/s], we have
             !//   ftot = ssProd(N) + fOM = ssProd(N) + ftot * OMssa
             !//   ftot = ssProd(N) / (1 - OMssa)
             !//
             !// Gantt et al (2015) multiplies by 6 to match observed values,
             !// but this produces too high POA emissions. They claim
             !// annual total emitted should be 6.3 Tg.
             !// The scaling depends on the choice of production, and is also
             !// somewhat dependent on resolution (3-5% change when halving
             !// or doubling resolution.
             Epoa = Epoa + scaleEPOA * OMssa * ssProd(N) / (1._r8 - OMssa)


             !write(6,'(3i4,3es12.4)') i,j,n,&
             !     6._r8 * OMssa * ssProd(N), OMssa,ssProd(N)


          end do !// do N = 1, N_SS_BINS

          !// Save emission array of POA [kg/s]
          EPOAseaspray(II,JJ,MP) = Epoa

          if (LDEBUG) then
             if (epoa .lt. 0._r8) then
                write(6,'(a,i3)') f90file//':'//subr//': neg Epoa'
                print*,Epoa,OMssa,ssProd,chla,wind10
                stop
             end if
          end if

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
 


    !// --------------------------------------------------------------------
  end subroutine emissions_ocean_organiccarbon
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine emissions_ocean_total()
    !// --------------------------------------------------------------------
    !// Print out current total POA flux as [Tg/yr].
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LBCOC
    use cmn_parameters, only: secYear
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    if (.not.LBCOC) return
    write(6,'(a,es12.3)') 'Instantaneous POA from ocean [Tg/yr]: ',&
         sum(EPOAseaspray) * secYear * 1.d-9
    !// --------------------------------------------------------------------
  end subroutine emissions_ocean_total
  !// ----------------------------------------------------------------------





  !// ----------------------------------------------------------------------
  subroutine add_oceanOCemis(BX,DT,MP)
    !// --------------------------------------------------------------------
    !// Routine to add oceanic carbon emissions.
    !// Routine is called from either emis4chem_oslo or SOURCE.
    !// IMPORTANT:
    !//   Note that the arguments are different things in the
    !//   two routines!:
    !//     emis4chem_oslo: BX is BEMIS [kg/s] and needs DT=1.d0
    !//     SOURCE:         BX is BTT [kg] and needs DT to be the time step.
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only:  MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TNAME
    use cmn_met, only: ZOFLE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: DT
    !// In/Out: Emissions (emis4chem_oslo) or tracer array (SOURCE)
    real(r8), intent(inout) :: BX(LPAR,NPAR,IDBLK,JDBLK)

    !// Local
    integer :: NTR, N, I,J, II,JJ
    real(r8) :: frac, DZ, ZH, topoH, plumeH, topZ3bot
    !// --------------------------------------------------------------------

    !// Skip volcanic emissions?
    if (.not. LOCEAN_OCARB) return

    do N = 1, NPAR_ocnPOA
       !// Species transport number
       NTR = ocnPOA_trsp_idx(N)

       !// Loop over latitude (J is global, JJ is block)
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ = J - MPBLKJB(MP) + 1

          !// Loop over longitude (I is global, II is block)
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II = I - MPBLKIB(MP) + 1

             BX(1,NTR,II,JJ) = BX(1,NTR,II,JJ) &
                  +  EPOAseaspray(II,JJ,MP) * ocnPOA_emisFraction(N) * DT

          end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    end do !// do N = 1, NPAR_ocnPOA

    !// --------------------------------------------------------------------
  end subroutine add_oceanOCemis
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module emissions_ocean
!//=========================================================================

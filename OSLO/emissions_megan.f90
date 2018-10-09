!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud,  October 2017
!//=========================================================================
!// MEGANv2.1 Biogenic emissions.
!//=========================================================================
module emissions_megan
  !// ----------------------------------------------------------------------
  !// MODULE: emissions_megan
  !// DESCRIPTION: MEGANv2.1 routines, fitted to CTM3.
  !// REFERENCE: Gunther et al, 2012 (G2012, doi:10.5194/gmd-5-1471-2012)
  !//
  !// Contains:
  !//   - subroutine megan_report
  !//   - subroutine add_meganBiogenic
  !//   - subroutine megan_input
  !//   - subroutine megan_get_co2
  !//   - subroutine megan_update_metdata
  !//   - subroutine megan_emis
  !//   - subroutine getWP
  !//   - subroutine megan_emis2file
  !//   MEGAN routines
  !//   - subroutine GAMMA_AGE
  !//   - subroutine GAMMA_CANOPY2 (replaces GAMMA_CANOPY)
  !//   - subroutine GAMMA_CANOPY
  !//   - function DIstomata
  !//   - function Ea1t99
  !//   - function Ealti99
  !//   - function Ea1p99
  !//   - function WaterVapPres
  !//   - function Stability
  !//   - subroutine GaussianIntegration
  !//   - subroutine SolarFractions
  !//   - subroutine WeightSLW
  !//   - subroutine CanopyRad
  !//   - subroutine CalcExtCoeff
  !//   - subroutine CalcRadComponents
  !//   - subroutine CanopyEB
  !//   - subroutine LeafEB
  !//   - function ConvertHumidityPa2kgm3
  !//   - function ResSC
  !//   - function LeafIROut
  !//   - function LHV
  !//   - function LeafLE
  !//   - function LeafBLC
  !//   - function LeafH
  !//   - function SvdTk
  !//   - function CalcEccentricity
  !//   - function UnexposedLeafIRin
  !//   - function ExposedLeafIRin
  !//   - subroutine SOILNOX
  !//   - function FERTLZ_ADJ
  !//   - subroutine GROWSEASON
  !//   - function prec_adj
  !//   - function getPulseType
  !//   - function precipfact
  !//
  !// Amund Sovde Haslerud, October/November 2017
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, NPAR, NRMETD
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Variables (2D) needed for MEGAN emissions.
  !// Only define what is not found in metdata (cmn_met.f90) or in
  !// other CTM3 fields.
  !// ----------------------------------------------------------------------
  !// Number of MEGAN species to be calculated
  integer, parameter :: N_MGN_SPC = 20
  !// Numper of Canopy types and characteristics
  !// Have included the 17th from MEGAN3
  integer, parameter :: N_MGN_PFT = 16, NrCha = 17
  !// Numper of speciated species in MEGAN. These will be lumped together
  !// in the Oslo CTM3. Have added NO2.
  integer, parameter :: N_SMAP_SPC = 151

  !// Canopy names
  character(len=16), dimension(N_MGN_PFT) :: MEGAN_CANOPY_NAME
  character(len=40), dimension(N_MGN_PFT) :: MEGAN_CANOPY_LONGNAME

  !// Canopy characteristics
  real(r8), dimension(NrCha, N_MGN_PFT) :: CanopyChar

  !// Emission factors for 20 species on N_MGN_PFT canopies
  real(r8), dimension(N_MGN_PFT, N_MGN_SPC) :: MEGAN_emisFactor
  !// Flag whether EFmaps will be used (not included)
  integer, dimension(N_MGN_SPC) :: use_EFmaps

  !// Molecular weight for the 20 species:
  real(r8), dimension(N_MGN_SPC) :: MGN_SPC_Mw

!//TO BE UPDATED
  logical, parameter :: SOIL_AVAIL = .false.

  !// MEGAN species names
  character(len=10), dimension(N_MGN_SPC) :: MEGAN_NAME

  !// Species data (Table 4)
  real(r8), dimension(N_MGN_SPC) :: &
       TDF_PRM, & ! Empirical constant (beta)
       LDF, &     ! light dependent fraction (LDF)
       Ctm1, &    ! Empirical constant (Ct1)
       CLeo, &    ! Emission class dependent empirical constant (Ceo)
       Anew, &    ! Empirical coefficients (new leaves)
       Agro, &    ! Empirical coefficients (growing leaves)
       Amat, &    ! Empirical coefficients (mature leaves)
       Aold, &    ! Empirical coefficients (old leaves)
       Rho        ! Density of some sort

  !// Speciated species
  real(r8), dimension(N_SMAP_SPC) :: MGN_SMAP_Mw
  integer, dimension(N_SMAP_SPC) :: MGN_SPC_MAP
  real(r8), dimension(N_MGN_PFT,N_SMAP_SPC) :: EFFS_all

  !// Map from speciated to CTM species
  integer, dimension(N_SMAP_SPC) :: MGN_SPC2CTM_MAP

  !// Flag for calculating MEGAN (initialise as false)
  logical :: LMEGAN = .false.

  !// Should consider writing T24,P24,T240,P240 to restart file:
  !// Average leaf temperature and PAR over 24 hours
  real(r8), dimension(NRMETD,IPAR,JPAR) :: TEMP24H, PAR24H
  integer :: init24 = 0
  integer :: INO = 0
  !// Average leaf temperature and PAR over 240 hours (10 days)
  !real(r8), dimension(IPAR,JPAR,10) :: TEMP10D, PAR10D

  !// 24-hour rain rate
  real(r8), dimension(NRMETD,IPAR,JPAR) :: PREC24H
  !// Precipitation pulse type
  integer, dimension(IPAR,JPAR) :: precPulseType
  !// Pulse start time [days]
  real(r8), dimension(IPAR,JPAR) :: precPulseStartTime
 

  !// Wilting point [m3/m3]
  real(r8), dimension(IPAR,JPAR) :: WPmap

  !// Monthly CO2 (vmr) for isoprene
  integer, parameter :: CO2nyears = 2018-1750
  real(r8), dimension(12,CO2nyears) :: CO2monthly
  real(r8) :: CO2thisMonth
  character(len=80) :: file_globalCO2

  !// Diagnostics (averages)
  real(r8), dimension(IPAR,JPAR,NPAR) :: AVG_MGN_CTM
  real(r8), dimension(IPAR,JPAR,N_MGN_SPC) :: AVG_MGN20
  real(r8) :: DTOPS

  !// Global scaling - set to give 572Tg isoprene for year 2000, using
  !// LAI climatology. Should be 0.466
  real(r8), parameter :: &
       Cce = 0.466_r8 !wrong first test 0.483_r8

  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'emissions_megan.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public megan_report, add_meganBiogenic, megan_input, megan_update_metdata, &
       megan_emis2file
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine megan_test()
    !// --------------------------------------------------------------------
    !// To be called from pmain.f90.
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    !use cmn_size, only: IDBLK, JDBLK, MPBLK
    use cmn_parameters, only: secYear
    use cmn_oslo, only: chem_idx, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !integer, intent(in) :: NDAY, NMET
    !// --------------------------------------------------------------------
    real(r8) :: MEGEM_tot(NPAR), MEGEM_CTM(NPAR)
    real(r8) :: MEGEM_tot2(N_MGN_SPC), MEGEM_ij(N_MGN_SPC)
    integer :: I, J
    !// --------------------------------------------------------------------

    MEGEM_tot(:) = 0._r8
    MEGEM_tot2(:) = 0._r8
    do J = 1, JPAR
       do I = 1, IPAR

          call megan_emis(I,J,MEGEM_CTM, MEGEM_ij)
          MEGEM_tot(:) = MEGEM_tot(:) + MEGEM_CTM(:)
          MEGEM_tot2(:) = MEGEM_tot2(:) + MEGEM_ij(:)
       end do
    end do

    write(6,'(a12,2es12.4)') 'MGN_ISOPRENE', &
         MEGEM_tot(trsp_idx(20)) * secYear * 1.e-9_r8, &
         MEGEM_tot2(1) * secYear * 1.e-9_r8

    !// --------------------------------------------------------------------
  end subroutine megan_test
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine megan_report(NDAY, NDAYI)
    !// --------------------------------------------------------------------
    !// To be called from daily_diag_output in diagnostics_general.f90
    !// (at the end of the day).
    !// Temporary routine to report global emissions from MEGAN
    !// for CTM-species.
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: secYear, secDay
    use cmn_chem, only: TNAME
    use cmn_diag, only: NDAY0
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NDAY, NDAYI
    !// --------------------------------------------------------------------
    real(r8) :: totEmis, dt, dtfac
    integer :: N
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_report'
    !// --------------------------------------------------------------------

    if (.not. LMEGAN) return

    !// Routine is called at the end of the day, so we add 1 to NDAY
    !// to get correct dt:
    dt = real(NDAY + 1 - NDAY0, r8) * secDay
    dtfac = secYear  / dt

    do N = 1, NPAR
       !// convert from kg to g
       totEmis = sum(AVG_MGN_CTM(:,:,N)) * 1.e3_r8
       if (totEmis .gt. 0._r8) then
          write(6,'(a,i3, 3es12.4)') f90file//':'//subr// &
               ':MEGAN_'//trim(TNAME(N))//' [g/yr]: ', &
               NDAY, totEmis*dtfac, totEmis, dt
       end if
    end do
    !// --------------------------------------------------------------------
  end subroutine megan_report
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine add_meganBiogenic(BEMIS, DT_IN, MP)
    !// --------------------------------------------------------------------
    !// To be called from emis4chem (emisdep4chem.f90) or source
    !// (source.f90).
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_size, only: LPAR, IDBLK, JDBLK
    use cmn_parameters, only: secYear
    use cmn_oslo, only: chem_idx, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    !//   BEMIS: unit is kg/s when called from emis4chem (emission rate),
    !//          but kg when called from source_uci.f90 (BTT is sent in)
    !//   DT_IN: Equals 1 if output is kg/s, and time step duration (s)
    !//          if output is kg.
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR, NPAR, IDBLK, JDBLK), intent(inout) :: BEMIS
    real(r8), intent(in) :: DT_IN
    integer, intent(in) :: MP
    !// --------------------------------------------------------------------
    real(r8) :: MEGEM_CTM(NPAR), MEGEM_SPC(N_MGN_SPC), DT4AVG
    integer :: I, J, II, JJ, N
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'add_meganBiogenic'
    !// --------------------------------------------------------------------

    !// Skip if MEGAN is not included
    if (.not. LMEGAN) return

    !// DT used for accumulating diagnose
    if (DT_IN .eq. 1._r8) then
       DT4AVG = DTOPS
    else
       DT4AVG = DT_IN
    end if

    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ = J - MPBLKJB(MP) + 1
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II = I - MPBLKIB(MP) + 1

          !// Find MEGAN emissions (kg/s)
          call megan_emis(I,J, MEGEM_CTM, MEGEM_SPC)

          !// Put emissions (kg/s) into surface level
          do N = 1, NPAR

             if (MEGEM_CTM(N) .gt. 0._r8) then
                !// If output is kg/s, DT_IN is 1, and if output
                !// is kg, DT_IN is time step duration.
                BEMIS(1,N,II,JJ) = BEMIS(1,N,II,JJ) &
                     + MEGEM_CTM(N) * DT_IN

                !// Diagnostic CTM species (kg)
                AVG_MGN_CTM(I,J,N) = AVG_MGN_CTM(I,J,N) &
                     + MEGEM_CTM(N) * DT4AVG
             end if

             !// Diagnostic MGN20 species (kg)
             AVG_MGN20(I,J,:) = AVG_MGN20(I,J,:) &
                  + MEGEM_SPC(:) * DT4AVG

          end do !// do N = 1, NPAR
       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine add_meganBiogenic
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine megan_input(LMEG, PATH)
    !// --------------------------------------------------------------------
    !// Initialise MEGAN - sets flag. Reads MEGAN input tables.
    !// PATH holds MEGAN specific files (e.g. PW).
    !// Called from emisutils based on STV entry.
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: secDay
    use cmn_size, only: LOSLOCTROP
    use cmn_ctm, only: NROPSM, IYEAR, JMON
    use cmn_sfc, only: StomRes
    use utilities, only: get_free_fileid
    use cmn_oslo, only: METHANEMIS, trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    logical, intent(in) :: LMEG
    character(len=*), intent(in) :: PATH

    character(len=160) :: filename
    character(len=160) :: LINE, CTMP
    character(len=130) :: BAR
    integer :: ifnr, I, N, TMP_IDX

    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_input'
    !// --------------------------------------------------------------------

    !// Make sure to turn off when chemistry is not included
    LMEGAN = LMEG .and. LOSLOCTROP

    !// Initialise diagnostics
    AVG_MGN_CTM(:,:,:) = 0._r8
    AVG_MGN20(:,:,:) = 0._r8
    DTOPS = secDay / real(NRMETD*NROPSM, r8)

    !// Initialise stomatal resistance
    StomRes(:,:,:) = 0._r8

    !// Skip this if MEGAN is not used
    if (.not. LMEGAN) return

    !// Read megan_tables.dat
    filename = 'tables/megan_tables.dat'
    ifnr = get_free_fileid()

    !//---------------------------------------------------------------------

    !// Make a bar of '>>>>'
    do I = 0, 12
       BAR((i*10+1):(i+1)*10) = '----------'
    end do
    write(6,'(a)') BAR

    open(ifnr,file=filename,form='formatted')

    !// Megan species
    read(ifnr,'(a80)') LINE !// Header
    write(6,'(a)') trim(LINE)
    write(6,'(a)') BAR
    do I = 1, 4
       read(ifnr,'(a80)') LINE !// Header
    end do
    !// Nr MEGAN-NAME beta  LDF   Ct1   Ceo   Anew  Agro  Amat  Aold  Rho
    do I = 1, N_MGN_SPC
       read(ifnr,*) N, MEGAN_NAME(I), TDF_PRM(I), LDF(I), Ctm1(I), CLeo(I), &
            Anew(I), Agro(I), Amat(I), Aold(I), Rho(I), MGN_SPC_Mw(I)
       write(6,'(i2,a15,10f9.2)') N, trim(MEGAN_NAME(I)), TDF_PRM(I), LDF(I), &
            Ctm1(I), CLeo(I), Anew(I), Agro(I), Amat(I), Aold(I), Rho(I), &
            MGN_SPC_Mw(I)
       !// Set index for NO
       if (trim(MEGAN_NAME(I)) .eq. 'NO') INO = I
    end do
    read(ifnr,*) LINE !// End
    if (LINE(1:5) .ne. '-END-') then
       write(6,'(a)') f90file//':'//subr//': Wrong end of table'
       stop
    end if
    write(6,'(a)') BAR

    !// Canopy types/names
    read(ifnr,'(a)') LINE !// Header
    write(6,'(a)') trim(LINE)
    write(6,'(a)') BAR
    do I = 1, 3
       read(ifnr,*) LINE !// Header
    end do
    do I = 1, N_MGN_PFT
       read(ifnr,*) N, MEGAN_CANOPY_NAME(I), MEGAN_CANOPY_LONGNAME(I)
       write(6,'(i2,1x,a16,a40)') N, MEGAN_CANOPY_NAME(I), MEGAN_CANOPY_LONGNAME(I)
    end do
    read(ifnr,*) LINE !// End
    if (LINE(1:5) .ne. '-END-') then
       write(6,'(a)') f90file//':'//subr//': Wrong end of table'
       stop
    end if
    write(6,'(a)') BAR

    !// Canopy characteristics
    read(ifnr,'(a)') LINE !// Header
    write(6,'(a)') trim(LINE)
    write(6,'(a)') BAR
    do I = 1, 2
       read(ifnr,*) LINE !// Header
    end do
    do I = 1, NrCha
       read(ifnr,*) N, canopyChar(I,:)
       write(6,'(i2,16f8.3)') N, canopyChar(I,:)
    end do
    read(ifnr,*) LINE !// End
    if (LINE(1:5) .ne. '-END-') then
       write(6,'(a)') f90file//':'//subr//': Wrong end of table'
       stop
    end if
    write(6,'(a)') BAR

    !// Emission factors
    read(ifnr,'(a)') LINE !// Header
    write(6,'(a)') trim(LINE)
    write(6,'(a)') BAR
    do I = 1, 5
       read(ifnr,*) LINE !// Header
    end do
    do I = 1, N_MGN_SPC
       read(ifnr,*) N, MEGAN_emisFactor(:,I), use_EFmaps(I)
       write(6,'(i2,16f9.2,i3)') N, MEGAN_emisFactor(:,I), use_EFmaps(I)
    end do
    read(ifnr,*) LINE !// End
    if (LINE(1:5) .ne. '-END-') then
       write(6,'(a)') f90file//':'//subr//': Wrong end of table'
       stop
    end if
    write(6,'(a)') BAR

    !// MEGAN speciated names, molecular weights and weighting factors
    !// Speciation based on species and canopy (EFFS_all) is listed after this loop,
    !// since the table is very wide for std.out.
    read(ifnr,'(a)') LINE !// Header
    write(6,'(a)') trim(LINE)
    read(ifnr,'(a)') LINE !// Header
    write(6,'(a)') trim(LINE)
    write(6,'(a)') BAR
    do I = 1, 3
       read(ifnr,*) LINE !// Header
    end do

    write(6,'(25x,a6,a4,2x,a9,2x,a12)') 'Mw    ','MG20','CTM #id  ','CTM #id used'

    do I = 1, N_SMAP_SPC
       read(ifnr,*) N, LINE, MGN_SMAP_Mw(I), MGN_SPC_MAP(I), CTMP, &
            TMP_IDX, EFFS_all(:,I)

       if (TMP_IDX .eq. 0) then
          !// No mapping tracer
          MGN_SPC2CTM_MAP(I) = 0
       else if (trsp_idx(TMP_IDX) .le. 0) then
          !// Not transported
          if (TMP_IDX .eq. 192 .and. &
               trsp_idx(12) .gt. 0) then
             !// Treat 192 as 12 when SOA not included
             MGN_SPC2CTM_MAP(I) = 12
          else
             !// Tracer not included in simulation - skip mapping
             MGN_SPC2CTM_MAP(I) = 0
          end if
       else if ( TMP_IDX .eq. 46 .and. (.not. METHANEMIS) ) then
          !// Special treatment for CH4
          MGN_SPC2CTM_MAP(I) = 0
       else
          !// Transported
          MGN_SPC2CTM_MAP(I) = TMP_IDX
       end if

       write(6,'(i3,1x,a20,1x,f6.2,2i4,a7,i4)') N, trim(LINE), MGN_SMAP_Mw(I), MGN_SPC_MAP(I), &
            TMP_IDX, '   --->', MGN_SPC2CTM_MAP(I)
    end do
    read(ifnr,*) LINE !// End
    if (LINE(1:5) .ne. '-END-') then
       write(6,'(a)') f90file//':'//subr//': Wrong end of table'
       stop
    end if
    write(6,'(a)') BAR

    !// MEGAN speciation species/canopy
    write(6,'(a)') 'Speciation based on species/canopy'
    do I = 1, N_SMAP_SPC
       write(6,'(i3,16f8.5)') I, EFFS_all(:,I)
    end do
    write(6,'(a)') BAR


    read(ifnr,*) LINE !// End
    read(ifnr,*) file_globalCO2
    write(6,'(a)') f90file//':'//subr//': global CO2 from '//trim(file_globalCO2)

    close(ifnr)


    !// Read PFTF - separate into MEGAN categories.

    !// Test MODIS first.

    !// Wilting Point in top soil layer (m2/mm3)
    !// Convert by 1/1000 to get m2/m3
    !// ECMWF 4-layer soil has top layer thickness of 0.07m, following by
    !// layer thicknesses of 0.21m, 0.27m and 1.89m.
    !// We will multiply top by 0.07 to get wilting point.

    call getWP(trim(PATH)//'wptop_150sec.nc')
    !// Convert to m3/m3
    WPmap(:,:) = 0.07e-3_r8 * WPmap(:,:)

    !// Initialise monthly CO2 (for isoprene)
    !// Should be read from file.
    CO2monthly(:,:) = 390._r8
    CO2thisMonth = 390._r8
    call megan_get_co2()

    CO2thisMonth = CO2monthly(JMON,IYEAR-1749)
    write(6,'(a,f9.3)') f90file//':'//subr//': CO2 now:',CO2thismonth
    !// --------------------------------------------------------------------
  end subroutine megan_input
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine megan_get_co2()
    !// --------------------------------------------------------------------
    !// Sets observed/estimated CO2.
    !//
    !// To be called from megan_input.
    !//
    !// Amund Sovde Haslerud, December 2017
    !// --------------------------------------------------------------------
    use ncutils, only: get_netcdf_dim, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer :: nTimeYears, nTimeMonths
    real(r8) :: oldCO2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_get_co2'
    !// --------------------------------------------------------------------

    CO2monthly(:,:) = 0.

    !// Use meteorological year or specific year (getYear)
    call get_netcdf_dim( file_globalCO2, 'years', nTimeYears)
    call get_netcdf_dim( file_globalCO2, 'months', nTimeMonths)

    call get_netcdf_var_2d(file_globalCO2, 'CO2', CO2monthly, nTimeMonths, nTimeYears)

    write(6,'(a)') f90file//':'//subr//': getting CO2 from file '// &
         trim(file_globalCO2)
    print*,CO2monthly(:,1)

    !// --------------------------------------------------------------------
  end subroutine megan_get_co2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine megan_update_metdata(NDAY, NMET, NDAYI)
    !// --------------------------------------------------------------------
    !// To be called from metdata_ecmwf.f90?
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: secDay
    use cmn_met, only: SFT, PhotActRad, PRECCNV, PRECLS, MYEAR
    use cmn_ctm, only: LFIXMET, JMON, JYEAR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NDAY, NMET, NDAYI
    !// --------------------------------------------------------------------
    integer :: N
    real(r8) :: oldCO2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_update_metdata'
    !// --------------------------------------------------------------------

    !// Skip if MEGAN is not included
    if (.not. LMEGAN) return

    if (NDAY .eq. NDAYI) then
       !// Possibly initialise TEMP24H, PAR24H
       if (init24 .eq. 0) then
          if (NMET .eq. 1) then
             !// Initialise all meteorological steps of TEMP24/PAR24
             do N = 1, NRMETD
                TEMP24H(N,:,:) = SFT(:,:)
                PAR24H(N,:,:) = PhotActRad(:,:)
             end do
             !// Zero accumulated rain - will set for step NMET below
             PREC24H(:,:,:) = 0._r8
             precPulseType(:,:) = 0
             precPulseStartTime(:,:) = 0._r8
          else
             !// Possibly keep initialising TEMP10D
             init24 = 1
          end if
       end if
    end if

    !// Set new value for step NMET
    TEMP24H(NMET,:,:) = SFT(:,:)
    PAR24H(NMET,:,:) = PhotActRad(:,:)
    PREC24H(NMET,:,:) = secDay/real(NRMETD,r8) * (PRECCNV(:,:,1) + PRECLS(:,:,1)) !// kg accumulated

    !// Set current CO2. Get year index starting from 1750:
    if (LFIXMET) then
       N = MYEAR - 1749
    else
       N = JYEAR - 1749
    end if
    if (N .lt. 1) then
       write(6,'(a)') f90file//':'//subr// &
            'Year index N is < 1 - needed to find CO2'
       stop
    else if (N .gt. CO2nyears) then
       write(6,'(a)') f90file//':'//subr// &
            'Year index N is too large - needed to find CO2'
       stop
    else
       oldCO2 = CO2thisMonth
       CO2thisMonth = CO2monthly(JMON,N)
       if (CO2thisMonth < 0._r8) then
          !// No CO2 was found for this month
          !// Keep old value
          CO2thisMonth = oldCO2
       end if
    end if

    !// --------------------------------------------------------------------
  end subroutine megan_update_metdata
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine megan_emis(I,J, EMIS_CTM, MGN20_EMIS)
    !// --------------------------------------------------------------------
    !// Driver for MEGAN emissions.
    !// The routine is a box model, first calculating MEGAN species, then
    !// doing the speciation into CTM3 tracers.
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code (emproc.F, ...)
    !// --------------------------------------------------------------------
    use cmn_parameters, only: CPI180
    use cmn_ctm, only: SOLDEC, SOLDIS, XDGRD, YDGRD, GMTAU, JDAY, JMON, &
         AREAXY, PLAND
    use cmn_chem, only: TMASS
    use cmn_met, only: PhotActRad, SFT, SFQ, SWVL1, STL1, P, SFU, SFV
    use cmn_sfc, only: LAI, ZOI, landSurfTypeFrac, LANDUSE_IDX, StomRes
    use utilities, only: LOCSZA
    use cmn_oslo, only: DINM, trsp_idx
    use utilities_oslo, only: GROWSEASON
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: I, J
    real(r8), dimension(NPAR), intent(out) :: EMIS_CTM
    real(r8), dimension(N_MGN_SPC), intent(out) :: MGN20_EMIS

    !// Locals
    real(r8), dimension(N_MGN_SPC) :: &
         EMIS_MEGAN, & ! Emission buffer
         GAM_I         ! Total correction factor
    real(r8), dimension(N_SMAP_SPC) :: &
         EMIS_MEGAN_SPC  ! Emission buffer
    real(r8), dimension(N_MGN_PFT) :: &
         canopyfraction ! canopy fractions
    real(r8) :: &
         LAT, LONG, &
         LAIp, &       ! Previous monthly LAI
         LAIc, &       ! Current monthly LAI
         dtLAIp, &     ! Duration of previous LAI dataset [days]
         TEMP, &       ! Temperature (K)
         PPFD, &       ! Calculated PAR [W/m2]
         GAM_LHT, &    ! LAI correction factor
         GAM_AGE, &    ! leaf age correction factor
         GAM_SM, &     ! Soil moisture correction factor
         GAM_CO2, &    ! CO2 correction factor
         GAM_TLD, &    ! Canopy correction factor (P,T) light dependent
         GAM_TLI, &    ! Canopy correction factor (P,T) light independent
         Ci, &         ! CO2 soil value
         CFNO, &    ! Correction factor for NO
         CFNOG, &   ! Correction factor for NO for grass
         WIND, & ! [m/s]
         PRES, & ! [Pa]
         QV, &   ! Humidity (mixing ratio, but assume specific humidity)
         DI, &   ! Drought-index (Palmer Drought Severity Index)
         PRECADJ, & ! precipitation adjustment factor
         COSSZA, &
         SZArad, & ! Solar zenith angle in radians
         SOLFX, &
         T24, T240, &
         P24, P240, &
         ADJUST_FACTOR_LI, ADJUST_FACTOR_LD, &
         GAMMA_TD, GAMMA_TI, &
         canopyFractionTotal, &
         landFrac, &
         soilMoisture, & !Soil Moisture (m3/m3)
         soilMoistureWiltingPoint, & !Soil Moisture Wilting Point [m2/m3]
         soilTemp, &  ! Soil Temperature (K)
         TMP1, TMP3, AREA, &
         accRain24, &
         ugm2hour_kgs

    !// Emission factors for 20 species on N_MGN_PFT canopies
    real(r8), dimension(N_MGN_PFT, N_MGN_SPC) :: emisFactor

    integer :: &
         soilType ! Soil Type (how many types?)

    integer :: S, I_PFT, N_SPEC, GDAY, GLEN
    integer :: N_CTM_CID, N_CTM

    !// Isoprene specific
    real(r8), parameter :: CO2conc = 390._r8 !// NEED to update to monthly since PI
    real(r8), parameter :: Ismax = 1.344_r8
    real(r8), parameter :: h = 1.4614_r8
    real(r8), parameter :: Cstar = 585._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_emis'
    !// --------------------------------------------------------------------

    !// If MEGAN is not included, the caller should test whether this
    !// routine should be called.

    EMIS_MEGAN(:) = 0._r8     !// Initialise emission array
    EMIS_MEGAN_SPC(:) = 0._r8 !// Initialise speciated emission array

    EMIS_CTM(:) = 0._r8       !// Initialise CTM species array (output)
    MGN20_EMIS(:) = 0._r8     !// Initialise species array (output)

    !// LAI current and previous
    LAIc = LAI(I,J,JMON)
    if (JMON .eq. 1) then
       LAIp = LAI(I,J,12)
       dtLAIp = 31._r8                 ! days in December
    else
       LAIp = LAI(I,J,JMON-1)
       dtLAIp = real(DINM(JMON-1), r8) ! days in previous month
    end if

    !// Initialise StomRes
    StomRes(:,I,J) = 0._r8

    !// Canopy fractions
    !01 'NT_EG_TEMP' 'Needleaf evergreen temperate tree  '
    !02 'NT_EG_BORL' 'Needleaf evergreen boreal tree     '
    !03 'NT_DC_BORL' 'Needleaf deciduous boreal tree     '
    !04 'BT_EG_TROP' 'Broadleaf evergreen tropical tree  '
    !05 'BT_EG_TEMP' 'Broadleaf evergreen temperate tree '
    !06 'BT_DC_TROP' 'Broadleaf deciduous tropical tree  '
    !07 'BT_DC_TEMP' 'Broadleaf deciduous temperate tree '
    !08 'BT_DC_BORL' 'Broadleaf deciduous boreal tree    '
    !09 'SB_EG_TEMP' 'Broadleaf evergreen temperate shrub'
    !10 'SB_DC_TEMP' 'Broadleaf deciduous temperate shrub' 
    !11 'SB_DC_BORL' 'Broadleaf deciduous boreal shrub   ' Could use poleward of 60
    !12 'GS_C3_COLD' 'Cold C3 grass                      ' 
    !13 'GS_C3_COOL' 'Cool C3 grass                      '
    !14 'GS_C3_WARM' 'Warm C3 grass                      '
    !15 'CORN      ' 'Corn                               ' Set to zero, emission factor and canopy specific are the same
    !16 'CROP      ' 'Other crops                        ' as in 16. Assume all crop, no corn.
    !canopyfraction(:) = PFTF(:,I,J)
    !// Start with MODIS (landSurfTypeFrac)
    !//  0=Water Bodies                        1=Evergreen Needleleaf Forests
    !//  2=Evergreen Broadleaf Forests         3=Deciduous Needleleaf Forests
    !//  4=Deciduous Broadleaf Forests         5=Mixed Forests
    !//  6=Closed Shrublands                   7=Open Shrublands
    !//  8=Woody Savannas                      9=Savannas
    !// 10=Grasslands                         11=Permanent Wetlands
    !// 12=Croplands                          13=Urban and Built-Up
    !// 14=Cropland/Natural Vegetation Mosaic 15=Permanent Snow and Ice
    !// 16=Barren or Sparsely Vegetated       17=Unclassified
    !01 -> 1 60S-60N  + 5 (25%)
    !02 -> 2 boreal (what about 60S-60N? add to 9 as broadleaf?)
    !03 -> 1 boreal  + 5 (25%)
    !04 -> 3 23S-23N  + 5 (25%)
    !05 -> 3 60S-23S, 23N-60N (check boreal?)  + 5 (25%)
    !06 -> 4 23S-23N + 5 (25%)
    !07 -> 4 60S-23S, 23N-60N + 5 (25%)
    !08 -> 4 90S-60S, 60N-90N + 5 (25%)
    !09 -> (6+7) 60S-60N assume 50% (possibly add 2 90S-60S, 60N-90N?) + 5 (25%MF)
    !10 -> (6+7) 60S-60N assume 50%
    !11 -> 3 90S-60S, 60N-90N (?) + (6+7) 90S-60S, 60N-90N
    !12 -> 10+8+9 90S-60S, 60N-90N
    !13 -> 10+8+9 60S-23S, 23N-60N
    !14 -> 10+8+9 23S-23N
    !15 -> zero
    !16 -> 12+14
    canopyfraction(:) = 0._r8
    if (LANDUSE_IDX .eq. 2) then
       if (abs(YDGRD(J)) .gt. 60._r8) then
          !// Set boreal stuff
          canopyfraction(2) = landSurfTypeFrac(1,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(3) = landSurfTypeFrac(3,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(8) = landSurfTypeFrac(4,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(11) = landSurfTypeFrac(6,I,J) &
                               + landSurfTypeFrac(7,I,J) &
                               + 0.25_r8 * landSurfTypeFrac(5,I,J) &
                               + landSurfTypeFrac(2,I,J)
          canopyfraction(12) = landSurfTypeFrac(8,I,J) &
                               + landSurfTypeFrac(9,I,J) &
                               + landSurfTypeFrac(10,I,J)
          canopyfraction(16) = landSurfTypeFrac(12,I,J) &
                               + landSurfTypeFrac(14,I,J)
       else if (abs(YDGRD(J)) .gt. 23._r8) then
          !// Set temperate stuff
          canopyfraction(1) = landSurfTypeFrac(1,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(5) = landSurfTypeFrac(2,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(7) = landSurfTypeFrac(4,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(9) = 0.5_r8 * (landSurfTypeFrac(6,I,J) + landSurfTypeFrac(7,I,J)) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(10) = 0.5_r8 * (landSurfTypeFrac(6,I,J) + landSurfTypeFrac(7,I,J))
          canopyfraction(13) = landSurfTypeFrac(8,I,J) &
                              + landSurfTypeFrac(9,I,J) &
                              + landSurfTypeFrac(10,I,J)
          canopyfraction(16) = landSurfTypeFrac(12,I,J) &
                              + landSurfTypeFrac(14,I,J)
       else
          !// Set tropical stuff
          canopyfraction(1) = landSurfTypeFrac(1,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(4) = landSurfTypeFrac(2,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(6) = landSurfTypeFrac(4,I,J) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(9) = 0.5_r8 * (landSurfTypeFrac(6,I,J) + landSurfTypeFrac(7,I,J)) &
                              + 0.25_r8 * landSurfTypeFrac(5,I,J)
          canopyfraction(10) = 0.5_r8 * (landSurfTypeFrac(6,I,J) + landSurfTypeFrac(7,I,J))
          canopyfraction(14) = landSurfTypeFrac(8,I,J) &
                               + landSurfTypeFrac(9,I,J) &
                               + landSurfTypeFrac(10,I,J)
          canopyfraction(16) = landSurfTypeFrac(12,I,J) + landSurfTypeFrac(14,I,J)
       end if
    else if (LANDUSE_IDX .eq. 3) then
       !// CLM4-PFT
       !// Note that the data from CLM4 has barren land as category 17.
       canopyfraction(1:N_MGN_PFT) = landSurfTypeFrac(1:N_MGN_PFT,I,J)
    else
       write(6,'(a,i5)') f90file//':'//subr// &
            ': No such LANDUSE_IDX available',LANDUSE_IDX
       stop
    end if


    CanopyFractionTotal = 0._r8
    do I_PFT = 1, N_MGN_PFT   !canopy types
       CanopyFractionTotal = CanopyFractionTotal + canopyfraction(I_PFT)
    end do !// do I_PFT = 1, N_MGN_PFT

    if (CanopyFractionTotal .gt. 0._r8 .and. LAIc .eq. 0._r8) then
       !write(6,'(a,2i4,2es12.3)') f90file//':'//subr// &
       !     ': Canopy/lai inconsistency',I,J, CanopyFractionTotal, LAIc
       !// Override so next test will return zero
       CanopyFractionTotal = 0._r8
    end if

    if (CanopyFractionTotal .eq. 0._r8) then
       return
    end if


    !// Emission factor
    emisFactor(:,:) = MEGAN_emisFactor(:,:)
    !// Possibly override if use_EFmaps is set

    !// Area
    AREA = AREAXY(I,J)

    !// Land fraction
    landFrac = PLAND(I,J)

    !// Surface (2m) temperature (K)
    TEMP = SFT(I,J)

    !// Surface pressure (Pa)
    PRES = P(I,J) * 100._r8

    !// Relative humidity [kg/kg, mixing ratio], but we assume specific humidity can
    !// be used (unsure about SFQ, use Q(I,J,1) directly?):
    QV = SFQ(I,J)

    !// Wind speed [m/s]
    WIND = sqrt(SFU(I,J)**2 + SFV(I,J)**2)

    !// Average temperature for past 24 hours
    T24 = sum(TEMP24H(:,I,J)) / real(NRMETD, r8)
    !// Average temperature for past 240 hours
    T240 = T24 ! sum(TEMP10D(:,I,J)) / 10._r8

    !// PAR [W/m2] - MEGAN converts to [umol/m2/s] later
    PPFD = PhotActRad(I,J)

    !// Average PAR for past 24 hours
    P24 = sum(PAR24H(:,I,J)) / real(NRMETD, r8)

    !// Average PAR for past 240 hours
    P240 = P24 ! sum(PAR10D(:,I,J)) / 10._r8

    !// Soil water volumetric level 1 [m3/m3]
    soilMoisture = SWVL1(I,J)
    !// Wilting point
    soilMoistureWiltingPoint = WPmap(I,J)

    !// Soil temperature level 1 [K]
    soilTemp = STL1(I,J)

    !// Precipitation factor for soil NOx
    !// When using soil temperature, moisture etc, it is possible to
    !// adjust soil NOx due to precipitation, depending on the pulse of
    !// rain.
    !// Get 24-hour accumulated rain [kg/day -> cm/day]
    !// Divide by water density to get [m3], then by area to get [m]. Density=1000, m2->cm2 = 100.
    accRain24 = sum(PREC24H(:,I,J)) * 1.e-1_r8 / AREA
    PRECADJ = prec_adj(I,J,accRain24)


    !// Drought index (PSDI)
    !// It seems the drought index was never used. It is left out of
    !// MEGAN3. I guess all the other gammas may cover some of the
    !// temperature effect, while soil moisture should also cover
    !// some of it.
    !// However, the soil moisture is simple and applies only to isoprene
    !// and soil NOx. Including DI will increase isoprene
    !// emisions due to the effect on leaf temperature, but will
    !// also lower the stomatal resistance (which reduces emissions).
    !// This was documented in the Final Report of Improving Modeled
    !// Biogenic Isoprene Emissions under Drought Conditions and
    !// Evaluating Their Impact on Ozone Formation, AQRP Project 14-030,
    !// November, 2015.
    !// Will set zero as assumed in MEGAN3.
    DI = 0._r8

    !// Solar zenith angle (rad) at UTC time GMTAU
    call LOCSZA(GMTAU,XDGRD(I),YDGRD(J),SOLDEC,SOLDIS, COSSZA,SOLFX)
    SZArad = acos(COSSZA)


    !// LAI correction factor - independent of species, calculate
    !// outside S-loop.
    !// Not mentioned in G2012, but in G2006 and used in MEGAN code.
    GAM_LHT = (0.49_r8 * LAIc) / ( (1._r8 + 0.2_r8 * (LAIc**2) )**0.5_r8 )

    !// Calculate activity factors for all MEGAN species
    do S = 1, N_MGN_SPC

       !// Light response
       !// Calculated through canopy.

       !// Temperature (light dependent, light independent, then total)
       !// Calculated through canopy.

       !// Leaf age activity
       !// For time step use duration of previous month
       CALL GAMMA_AGE(S, LAIp, LAIc, dtLAIp, T24, GAM_AGE)


       !// Soil moisture impact on activity - will get new value for isoprene
       GAM_SM = 1._r8
       !// CO2 impact on activity - will get new value for isoprene
       GAM_CO2 = 1._r8


       !// Isoprene specific adjustments
       !// Not in MEGAN code, only in G2012
       if (trim(MEGAN_NAME(S)) .eq. 'ISOP') then
          !// Soil moisture impact on activity
          if (soilMoisture .lt. soilMoistureWiltingPoint) then
             GAM_SM = 0._r8
          else if (soilMoisture .lt. (soilMoistureWiltingPoint + 0.04_r8)) then
             GAM_SM = (soilMoisture - soilMoistureWiltingPoint) / 0.04_r8
          else
             GAM_SM = 1._r8
          end if

          !// CO2 impact on activity
          !// Heald, C. L., Wilkinson, M. J., Monson, R. K., Alo, C. A., Wang,
          !// G., and Guenther, A.: Response of isoprene emission to ambient
          !// CO2 changes and implications for global budgets, Glob. Change
          !// Biol., 15, 1127-1140, 2009, doi:10.1111/j.1365-2486.2008.01802.x
          !// See also Wilkinson et al, 2008, doi:10.1111/j.1365-2486.2008.01803.x
          if (CO2thisMonth .ne. 400._r8) then
             Ci = 0.7_r8 * CO2thisMonth
             GAM_CO2 = Ismax - ( (Ismax * Ci**h) / (Cstar**h + Ci**h) )
          end if
       end if



       !// Canopy (sums up over PFTs)
       !// Done by separate routine (see below for the old code)
       call GAMMA_CANOPY2(JDAY, SZArad, TEMP, T24, T240, &
                     PPFD, P24, P240, &
                     WIND, QV, LAIc, Pres, DI, &
                     canopyFraction, canopyFractionTotal, &
                     S, MEGAN_NAME(S), I, J, &
                     GAM_TLD, GAM_TLI)

       !// Old canopy code
!       ADJUST_FACTOR_LD = 0._r8
!       ADJUST_FACTOR_LI = 0._r8
!       
!       if (CanopyFractionTotal .gt. 0._r8) then
!
!          do I_PFT = 1, N_MGN_PFT   !canopy types
!             if (canopyFraction(I_PFT) .ne. 0._r8) then
!                CALL GAMMA_CANOPY(JDAY, SZArad, TEMP, T24, T240, &
!                     PPFD, P24, P240, &
!                     WIND, QV, I_PFT, LAIc, PRES, DI, &
!                     S, MEGAN_NAME(S), &
!                     GAMMA_TD, GAMMA_TI)
!                ADJUST_FACTOR_LD = ADJUST_FACTOR_LD + &
!                     canopyFraction(I_PFT) * GAMMA_TD
!                ADJUST_FACTOR_LI = ADJUST_FACTOR_LI + &
!                     canopyFraction(I_PFT) * GAMMA_TI
!                if (PPFD .gt. 0._r8 .and. (GAMMA_TI+GAMMA_TD) .eq. 0._r8) then
!                   write(6,'(a,2i4,5es12.3)') f90file//':'//subr// &
!                        'Have light but zero TD',I,J,PPFD,P24,GAMMA_TD, GAMMA_TI, LAIc
!                   stop
!                end if
!             end if
!          end do !// do I_PFT = 1, N_MGN_PFT
!          GAM_TLD = ADJUST_FACTOR_LD / CanopyFractionTotal
!          GAM_TLI = ADJUST_FACTOR_LI / CanopyFractionTotal
!
!       else if (CanopyFractionTotal .eq. 0._r8) then
!          GAM_TLD = 1._r8
!          GAM_TLI = 1._r8
!       else if (CanopyFractionTotal .lt. 0._r8) then
!          write(6,'(a,2i5,es16.6)') f90file//':'//subr// &
!              'CanopyFractionTotal is less than 0. I,J: ',I,J,CanopyFractionTotal
!          stop
!       end if



       !// Total correction factor
       GAM_I(S) = GAM_AGE * GAM_SM * GAM_CO2 * RHO(S) &
            * ( max(0._r8, 1._r8 - LDF(S)) * GAM_TLI * GAM_LHT &
                + LDF(S) * GAM_TLD)


    end do !// do S = 1, N_MGN_SPC


    !// Estimate CFNO and CFNOG (Soil NOx adj)
    call SOILNOX( JDAY, TEMP, SOIL_AVAIL, soilType, soilMoisture, soilTemp, &
                 LAIc, YDGRD(J), PRECADJ, CFNO, CFNOG )


    !// --------------------------------------------------------------------
    !// Next step is to multiply with emission factors for each canopy type.
    !// This step is combined with the speciation into 150 species, because
    !// the speciation also depends on canopy types.
    !// Unit will be [ug/m2/h].
    !// --------------------------------------------------------------------

    !// Use EF from G2012 table 2. Could implement EFmaps, but they are very large.

    !// Converting from ug/m2/hour to kg/s
    ugm2hour_kgs = AREA / 3600._r8  *  1.e-9_r8

    !// --------------------------------------------------------------------
    !// Only the 20 species: This is not used for anything, and could be
    !// removed. However, it is nice to have the possibility of checking these
    !// totals later.
    do S = 1, N_MGN_SPC
       !// Using only EF from table 2
       do I_PFT = 1, N_MGN_PFT
          !// Multiply by GAM_I after loop
          EMIS_MEGAN(S) = EMIS_MEGAN(S) &
               + emisFactor(I_PFT,S) * canopyfraction(I_PFT)
       end do !// do T = 1, N_MGN_PFT
       EMIS_MEGAN(S) = EMIS_MEGAN(S) * GAM_I(S)
    end do !// do S = 1, N_MGN_SPC

    !// Convert from ug/m2/h to kg/s
    MGN20_EMIS(:) = EMIS_MEGAN(:) * ugm2hour_kgs
    !// --------------------------------------------------------------------




    !// Speciation
    !// --------------------------------------------------------------------
    do N_SPEC = 1, N_SMAP_SPC

       !// CTM chemical index number
       N_CTM_CID = MGN_SPC2CTM_MAP(N_SPEC)

       !// If species is not used in CTM, skip it
       if (N_CTM_CID .eq. 0) cycle

       !// CTM transport number
       N_CTM = trsp_idx(N_CTM_CID)

       !// Skip if species not transported
       if (N_CTM .le. 0) cycle

       !// index for the 20 species
       S = MGN_SPC_MAP(N_SPEC)

       !// Factor used in calculations
       if (use_EFmaps(S) .eq. 0) then
          TMP1 = 1._r8
       else
          !// For EFmaps, canopyFractionTotal is often not 1.
          !// In this case emisFactor and canopyfraction is used in
          !// a weighting factor, however, since EFmaps already has taken
          !// the canopy (land fraction) into account, it is not necessary
          !// to multiply by land fraction.
          TMP1 = canopyFractionTotal
       end if


       !// Treat all species but NO
       if (S .ne. INO) then
          !// Using EFmaps, emisFactor is independent of canopy (all values are the same),
          !// while for using G2012 table values, they depend on canopy.
          !// Here we have EFmaps values set in emisFactor, and all values are the
          !// same, so that it can be included in thet PFT-loop.
          !// The equation is thus the same for G2012 emission factors
          !// and EFmaps:
          TMP3 = 0._r8
          do I_PFT = 1, N_MGN_PFT
             TMP3 = TMP3 + emisFactor(I_PFT,S) &
                             * EFFS_ALL(I_PFT,N_SPEC) * canopyfraction(I_PFT)
          end do

          if (TMP1 .eq. 0._r8) then
             !// Only possible for EFmaps
             EMIS_MEGAN_SPC(N_SPEC) = 0._r8
          else
             !// For both emission factor (TMP1 = 1) and EFmaps (TMP1 = canopyFractionTotal/landFrac)
             EMIS_MEGAN_SPC(N_SPEC) = GAM_I(S) * TMP3 / TMP1
          end if

       else
          !// NO special treatment
          !// Note that emission factor for NO is actually NO, while original MEGAN code uses N.
          call GROWSEASON(JDAY, YDGRD(J), GDAY, GLEN)

          !// Check for growing season
          !// Always use CFNOG for 1-14. For growing season use CFNO for crop,
          !// but for non-growing season, use CFNOG for everywhere - override
          !// crop with grass warm = 14.

          !// Always
          TMP3 = 0._r8
          do I_PFT = 1, 14
             TMP3 = TMP3 + emisFactor(I_PFT,S) * EFFS_all(I_PFT,N_SPEC) &
                           * canopyFraction(I_PFT) * CFNOG
          end do

          if (GDAY .eq. 0) then

             !// Non-growing season - override crop with grass warm = 14
             do I_PFT = 15, N_MGN_PFT
                !// Override crop with grass warm
                TMP3 = TMP3 + emisFactor(14,S) * EFFS_all(I_PFT,N_SPEC) &
                              * canopyFraction(I_PFT) * CFNOG
             end do

          else if (GDAY .gt. 0 .and. GDAY .le. 366) then

             !// Growing season - use emisFactor for all I_PFT
             do I_PFT = 15, N_MGN_PFT
                TMP3 = TMP3 + emisFactor(I_PFT,S) * EFFS_all(I_PFT,N_SPEC) &
                                 * canopyFraction(I_PFT) * CFNO
             end do

          else
             write(6,'(a,i5)') f90file//':'//subr//': Bad GDAY: ',GDAY
             stop
          end if

          !// Set emissions [ug/m2/h]
          if (TMP1 .eq. 0._r8) then
             EMIS_MEGAN_SPC(N_SPEC) = 0._r8
          else
             EMIS_MEGAN_SPC(N_SPEC) = GAM_I(S) * TMP3 / TMP1
          end if

       end if !// if (S .ne. INO) then


       !// Add to CTM emission array as kg/s - convert using molecular mass.
       !// This may seem strange. MEGAN calculates a total mass, not
       !// number of molecules; the number of molecules are then found from
       !// the molecular mass of the 150 species.
       !// These are then converted to CTM species.
       EMIS_CTM(N_CTM) = EMIS_CTM(N_CTM) &
            + EMIS_MEGAN_SPC(N_SPEC) &
              * ugm2hour_kgs &
              / MGN_SMAP_Mw(N_SPEC) * TMASS(N_CTM)

       !// I did consider rather using / MGN_SPC_Mw(S) * TMASS(N_CTM) on the
       !// previous line, but decided to include conversions due to
       !// molecular weights in the emission factors.

    end do !// do N_SPEC = 1, N_SMAP_SPC
    !// --------------------------------------------------------------------



    !// --------------------------------------------------------------------
  end subroutine megan_emis
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getWP(wp_file)
    !// --------------------------------------------------------------------
    !// Reads wilting point file from MEGAN (http://cdp.ucar.edu) and
    !// interpolates to CTM resolution.
    !//
    !// Amund Sovde Haslerud, October 2017
    !// --------------------------------------------------------------------
    use cmn_parameters, only: CPI180, A0
    use cmn_ctm, only: AREAXY, XDEDG, YDEDG
    use regridding, only: e_grid
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_2d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), intent(in) :: wp_file
    !// --------------------------------------------------------------------
    integer            :: status, ncid, nDims, nVars, nAtts
    integer            :: nLon, nLat, I,J
    real(r8), allocatable, dimension(:) :: inLon, inLat, inLatE, inLonE, &
         XYBOX
    real(r8), allocatable, dimension(:,:) :: wp_in, wp_in2
    real(r8) :: dx, dy
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'getWP'
    !// --------------------------------------------------------------------

    !// Check resolution (latitude/longitude/time)
    !// This routine allocates inLon/inLat/inTime
    call get_netcdf_var_1d( wp_file, 'lon',  inLon  )
    call get_netcdf_var_1d( wp_file, 'lat',  inLat  )

    nLon  = size(inLon)
    nLat  = size(inLat)

    !// Deallocate all local variables
    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)

    allocate( inLonE(nLon+1), inLatE(nLat+1), XYBOX(nLat), &
         wp_in(nLon,nLat), wp_in2(nLon,nLat) )

    dx = 360._r8 / real(nLon, r8)
    dy = 180._r8 / real(nLat, r8)
    inLonE(1) = -180._r8
    do I = 2, nLon + 1
       inLonE(I) = inLonE(I-1) + dx
    end do

    inLatE(1) = -90._r8
    do J = 2, nLat + 1
       inLatE(J) = inLatE(J-1) + dy
    end do

    !// Emission variable (m2/mm3)
    call get_netcdf_var_2d( wp_file, &
         'Wilting_point_for_top_soil_layer_(m2_per_mm3)', &
         wp_in(:,:), nLon,nLat )

    do J = 1, nLat
       do I = 1, nLon
          if (wp_in(I,J) .lt. -1.e38_r8) wp_in(I,J) = 0._r8
       end do
    end do

    !// Area of input data
    do J = 1, nLat
       XYBOX(J) =  A0*A0 * CPI180*dx &
            * (sin(CPI180*inLatE(J+1)) - sin(CPI180*inLatE(J)))
    end do

    !// Upside down and multiply by area
    do J = 1, nLat
       wp_in2(:,J) = wp_in(:,nLat-J+1) * XYBOX(J)
    end do

    call E_GRID(wp_in2, inLonE, inLatE, nLon, nLat, &
         WPmap(:,:), XDEDG, YDEDG, IPAR, JPAR, 1)

    !// Divide by new area
    WPmap(:,:) = WPmap(:,:) / AREAXY(:,:)

    deallocate( inLonE, inLatE, XYBOX, wp_in, wp_in2 )

    !// --------------------------------------------------------------------
  end subroutine getWP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine megan_emis2file(NDAY)
    !// --------------------------------------------------------------------
    !// Use these files in addition to emis_accu*nc
    !// To be called from tnd_emis2file in diagnostics_general.f90.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_parameters, only: secDay
    use cmn_ctm, only: AREAXY, &
         JYEAR, JMON, JDATE, XDGRD, YDGRD, XDEDG, YDEDG
    use cmn_chem, only: TNAME
    use cmn_diag, only: NDAY0,JYEAR0,JMON0,JDATE0, &
         nc4deflate_global, nc4shuffle_global
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: NDAY
    !// --------------------------------------------------------------------
    integer :: N
    real(r8) :: DTDIAG
    integer :: &
         lat_dim_id, lon_dim_id, time_dim_id, &
         lat_id, lon_id, time_id, &
         ilat_dim_id, ilon_dim_id, &
         ilat_id, ilon_id, &
         areaxy_id, &
         date_size_dim_id, &
         start_time_id, end_time_id, & !IDs for start/end dates for average
         start_day_id, end_day_id, &   !IDs for start/end day (NDAY)
         status, ncid
    integer, dimension(NPAR) :: comps_id

    character(len=8) :: datestamp, datestamp0
    integer, dimension(6) :: start_time, end_time
    character(len=160) :: filename

    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'megan_emis2file'
    !// --------------------------------------------------------------------

    !// Skip if not included
    if (.not. LMEGAN) return

    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
    write(datestamp0(1:8),'(i4.4,2i2.2)') JYEAR0,JMON0,JDATE0
    start_time = (/JYEAR0,JMON0,JDATE0,0,0,0/)
    end_time   = (/JYEAR, JMON, JDATE, 0,0,0/)

    filename = 'emis_accumulated_3d_megan_'//datestamp0//'_'//datestamp//'.nc'
    write(6,'(a)') f90file//':'//subr//': creating file: '//trim(filename)
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')

    !//File headers
    status=nf90_put_att(ncid,nf90_global,'title', &
         'Accumulated MEGAN emissions in Oslo CTM3')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': file header')
    !//---------------------------------------------------------------------

    !// Define spatial dimensions (lat, lon, lev)
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')
    !// Define spatial dimensions (ilat, ilon, ilev)
    status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon dim')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,"date_size",size(start_time),date_size_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')

    !//---------------------------------------------------------------------

    !// Define the lon/lat/lev
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')

    !// Define ilon/ilat/ilev
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')

    !// Putting attributes to lon/lat/lev variables
    status = nf90_put_att(ncid,lon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    !// Putting attributes to ilon/ilat variables
    status = nf90_put_att(ncid,ilon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')

    !// Need info on accumulated time period
    status = nf90_def_var(ncid,"delta_time",nf90_double,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define time variable')
    status = nf90_put_att(ncid,time_id,'units', &
            'accumulation period given in seconds')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units time_delta')

    !// Start date for accumulating data - START_TIME
    status = nf90_def_var(ncid,"START_TIME", nf90_int, &
         date_size_dim_id, start_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define START_TIME variable')
    status = nf90_put_att(ncid,start_time_id,'description', &
         'Start date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description START_TIME')

    !// End date for accumulating data - END_TIME
    status = nf90_def_var(ncid,"END_TIME", nf90_int, &
         date_size_dim_id, end_time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define END_TIME variable')
    status = nf90_put_att(ncid,end_time_id,'description', &
         'End date [YYYY,MM,DD,hh,mm,ss] for accumulating data.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description END_TIME')

    !// Grid area (r8), deflate netcdf4
    status = nf90_def_var(ncid,"gridarea", nf90_double, &
         (/lon_dim_id, lat_dim_id/), areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')

    !// Define data component by component
    comps_id(:) = -1
    do N = 1, NPAR
       if (sum(AVG_MGN_CTM(:,:,N)) .gt. 0._r8) then
          status = nf90_def_var(ncid,trim(TNAME(N)), nf90_float, &
               (/lon_dim_id, lat_dim_id /), comps_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define '//trim(TNAME(N))//' variable')
          status = nf90_def_var_deflate(ncid, comps_id(N), &
               nc4shuffle, 1, nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate '//trim(TNAME(N))//' variable')
          status = nf90_put_att(ncid,comps_id(N),'units','accumulated kg of species')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute units '//trim(TNAME(N)))
       end if !// if (totemis(N) .gt. 0._r8) then
    end do

    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')

    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')

    !// delta-time in seconds
    DTDIAG = real(NDAY + 1 - NDAY0, r8) * secDay
    status = nf90_put_var(ncid,time_id,DTDIAG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NRAVG')
    !// Start time
    status = nf90_put_var(ncid,start_time_id,start_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting START_TIME')
    !// End time
    status = nf90_put_var(ncid,end_time_id,end_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting END_TIME')

    !// Grid area
    status = nf90_put_var(ncid, areaxy_id, AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting gridarea')

    !// Component by component
    do N = 1, NPAR
       if (comps_id(N) .gt. 0) then
          status = nf90_put_var(ncid, comps_id(N),AVG_MGN_CTM(:,:,N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting '//trim(TNAME(N)))
       end if
    end do !// do N = 1, NPAR

    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------

    !// Re-initialize
    AVG_MGN_CTM(:,:,:) = 0._r8
    AVG_MGN20(:,:,:) = 0._r8

    !// --------------------------------------------------------------------
  end subroutine megan_emis2file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  !// MEGAN subroutines, fitted to Oslo CTM3.
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine GAMMA_AGE( SPEC_NUM, &
                        LAIp, LAIc, TSTEP_LAI, T24, GAM_A )
    !// --------------------------------------------------------------------
    ! GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
    ! Gunther et al., 2006, Eq(16-19)
    !// --------------------------------------------------------------------
    !    where Fnew = new foliage fraction
    !          Fgro = growing foliage fraction
    !          Fmat = mature foliage fraction
    !          Fold = old foliage fraction
    !          Anew = relative emission activity for new foliage
    !          Agro = relative emission activity for growing foliage
    !          Amat = relative emission activity for mature foliage
    !          Aold = relative emission activity for old foliage
    !          ti = number of days between budbreak and the induction of
    !               emission
    !          tm = number of days between budbreak and the initiation of
    !               peak emissions rates
    !          Tt = average temperature (K) near top of the canopy during
    !               current time period (daily ave temp for this case)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, to fit Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: SPEC_NUM
    real(r8), intent(in) :: LAIp, LAIc, TSTEP_LAI, T24
    !// Output
    real(r8), intent(out) :: GAM_A

    !// Locals
    real(r8) :: Fnew, Fgro, Fmat, Fold
    real(r8) :: ti,tm        ! number of days between budbreak
                             ! and induction of emission, 
                             ! initiation of peak emissions rates

    !// Calculate foliage fraction
    if (LAIp .lt. LAIc) then

       if (T24 .le. 303._r8) then
          ti = 5._r8 + 0.7_r8 * (300._r8 - T24)
       else
          ti = 2.9_r8
       end if
       tm = 2.3_r8 * ti
 
       !// Calculate Fnew and Fmat, then Fgro and Fold
       !// Fnew
       if (ti .ge. TSTEP_LAI) then
          Fnew = max(0._r8, 1._r8 - LAIp / LAIc )
       else
          Fnew = ti / TSTEP_LAI * max(0._r8, 1._r8 - LAIp / LAIc )
       end if
 
       !// Fmat
       if (tm .GE. TSTEP_LAI) then
          Fmat = LAIp / LAIc
       else
          Fmat = LAIp / LAIc + (TSTEP_LAI - tm) / TSTEP_LAI &
                 * max(0._r8, 1._r8 - (LAIp / LAIc))
       end if

       Fgro = max(0._r8, 1._r8 - Fnew - Fmat)
       Fold = 0._r8

    else if (LAIp .eq. LAIc)  then

       Fnew = 0._r8
       Fgro = 0.1_r8
       Fmat = 0.8_r8
       Fold = 0.1_r8

    else if (LAIp .gt. LAIc) then
 
       Fnew = 0._r8
       Fgro = 0._r8
       Fold = ( LAIp - LAIc ) / LAIp
       Fmat = max(0._r8, 1._r8 - Fold)

    end if

    !// Calculate GAMMA_A
    GAM_A = Fnew * Anew(SPEC_NUM) &
            + Fgro * Agro(SPEC_NUM) &
            + Fmat * Amat(SPEC_NUM) &
            + Fold * Aold(SPEC_NUM)

    !// --------------------------------------------------------------------
  end subroutine GAMMA_AGE
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine GAMMA_CANOPY2(JDAY, SZArad, TairK0, T24, T240, &
                          PPFD0, PPFD24, PPFD240, Wind0, &
                          Humidity,LAI,Pres,DI, &
                          canopyFraction, canopyFractionTotal, &
                          SPEC_NUM, SPCNAME, ILOC, JLOC, &
                          GAM_TLD, GAM_TLI)
    !// --------------------------------------------------------------------
    !// Modified GAMME_CE from MEGANv2.10 code.
    !// Loops over all PFTs.
    !//
    !// Amund Sovde Haslerud, October-December 2017
    !//   Based on MEGANv2.10 code, to fit Oslo CTM3
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    use cmn_sfc, only: StomRes
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: JDAY, SPEC_NUM, ILOC, JLOC
    real(r8), intent(in) :: &
         SZArad, &   ! Solar zenith angle in radians
         TairK0, &   ! Surface temperature [K]
         Pres, &     ! Surface pressure [hPa]
         Wind0, &    ! Surface wind [m/s]
         Humidity, & ! Humidity, supposedly mixing ratio [kg(H2O)/kg(dry air)]
                     ! but using specific humidity [kg(H2O)/kg(moist air)].
         LAI, &      ! Leaf area index
         DI, &       ! Drought index (PSDI)
         PPFD0, &    ! PPFD current [W/m2]
         PPFD24, &   ! PPFD averaged over last 24 hours [W/m2]
         PPFD240, &  ! PPFD averaged over last 10 days [W/m2]
         T24, &      ! Surface temperature averaged over last 24 hours [K]
         T240, &     ! Surface temperature averaged over last 10 days [K]
         canopyFractionTotal
    real(r8), dimension(N_MGN_PFT), intent(in) :: canopyFraction
    character(len=*), intent(in) :: SPCNAME

    !// Output
    real(r8), intent(out) :: &
         GAM_TLD, & ! Light dependent
         GAM_TLI    ! Light independent

    real(r8) :: &
         Ea1Canopy, ADJUST_FACTOR_LD, & ! Light dependent local
         EatiCanopy, ADJUST_FACTOR_LI   ! Light independent local

    ! Local variables
    integer :: L, I_PFT

    real(r8) :: Sinbeta

    !// Number of canopy layers
    integer, parameter :: Layers = 5

    real(r8), parameter :: &
         !ConvertWm2toUmolm2s = 4.766_r8, &      ! Convert radiation from [W/m2] to [umol/m2/s1]
         SolarConstant = 1367._r8, &       ! Solar constant [W/m2]
         WaterAirRatio = 18.016_r8 / M_AIR, & ! Ratio between water and air molecules
         pstd_Sun = 200.0_r8, &
         pstd_Shade = 50.0_r8

    real(r8), dimension(Layers) ::  VPgausWt, VPgausDis2, &
         VPgausDis, VPslwWT, Sunfrac, QdAbsV, QsAbsV, QdAbsn, &
         QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD, &
         ShadePPFD, TairK, HumidairPa, Wind, Sunleaftk, SunleafSH, &
         SunleafLH,SunleafIR, Shadeleaftk, ShadeleafSH, &
         ShadeleafLH,ShadeleafIR,Ea1pLayer, Ea1tLayer, Ea1Layer, &
         EatiLayer

    real(r8) :: Solar, Maxsolar, PPFD, QbAbsn, &
         Trate, StomataDI, Qbeamv,Qdiffv, Qbeamn, Qdiffn, &
         QbAbsV,Ea1tCanopy, Ea1pCanopy, &
         HumidairPa0

    real(r8), dimension(Layers) :: sr_sun, sr_shade
    real(r8) :: sr_avg
                       
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'GAMMA_CANOPY2'
    !// --------------------------------------------------------------------


    !// Loop over canopy types
    if (CanopyFractionTotal .gt. 0._r8) then

       !// Need solar elevation angle (beta), which is related to solar zenith
       !// angle (SZA) as  sin(beta) = cos(SZA):
       sinBETA = cos(SZArad)

       !// CanopyRad will convert from W/m2 to umol/m2/s using different
       !// value for sunlit and shaded leaves. Thus, we have to assume
       !// PPFD should have units W/m2.
       PPFD       = PPFD0 
       StomataDI  = DIstomata(DI)
       !// MEGAN code did 
       !//    Solar      = PPFD / 2.25_r8
       !// and previously had commented out
       !//    !Solar      = PPFD /ConvertWm2toUmolm2s*2
       !//    !Solar      = PPFD/2.1
       !// which could be a way to convert umol/m2/s to W/m2 which can later
       !// be converted back to umol/m2/s in CanopyRad.
       Solar      = PPFD
       Maxsolar   = Sinbeta  * SolarConstant * CalcEccentricity(JDAY)

       call GaussianIntegration(VPgausWt, VPgausDis, VPgausDis2, Layers)

       call SolarFractions(Solar, Maxsolar, Qdiffv, Qbeamv, Qdiffn, Qbeamn)

       call WeightSLW(VPgausDis, VPgausWt, LAI, Layers, VPslwWT)

       HumidairPa0  =  WaterVapPres(Humidity, Pres, WaterAirRatio)

       ADJUST_FACTOR_LD = 0._r8
       ADJUST_FACTOR_LI = 0._r8

       do I_PFT = 1, N_MGN_PFT   !canopy types
          if (canopyFraction(I_PFT) .ne. 0._r8) then

             call CanopyRad(VPgausDis, Layers, LAI, Sinbeta, Qbeamv, &
                  Qdiffv, Qbeamn, Qdiffn, I_PFT, Sunfrac, &
                  QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv, &
                  ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD)

             Trate    =  Stability(I_PFT, Solar)

             !// Possible average of stomatal resistance [s/m]
             if (SPEC_NUM .eq. 1) then
                sr_avg = 0._r8
                do L = 1, layers    
                   sr_sun(L) = ResSC(SunPPFD(L), StomataDI)
                   sr_shade(L) = ResSC(ShadePPFD(L), StomataDI)
                   sr_avg = sr_avg &
                        + ( sr_sun(L) * Sunfrac(L) &
                              + sr_shade(L) * (1._r8 - Sunfrac(L)) ) &
                        * VPslwWT(L) * VPgausWt(L)
                end do
                !// Store in StomRes(NVGPAR,IPAR,JPAR), so we have one
                !// resistance per canopy type.

                !// Remember to initialise StomRes at the top of megan_emis
                !// Canopy is calculated for each SPEC_NUM, but shold not
                !// change from one to another.
                StomRes(I_PFT,ILOC,JLOC) = sr_avg
             end if

             call CanopyEB(Trate, Layers, VPgausDis, I_PFT, &
                  StomataDI, TairK, HumidairPa, Wind, SunPPFD, &
                  ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, &
                  Sunleaftk, SunleafSH, SunleafLH, SunleafIR, &
                  Shadeleaftk,ShadeleafSH,ShadeleafLH,ShadeleafIR, &
                  Wind0, TairK0, HumidairPa0)

             !// Initialize for this canopy
             Ea1tLayer(:) = 0._r8
             Ea1pLayer(:) = 0._r8
             Ea1Layer(:) = 0._r8
             EatiLayer(:) = 0._r8

             do L = 1, layers    

                Ea1tLayer(L) = &
                     Ea1t99(Sunleaftk(L), T24, T240,SPEC_NUM) * Sunfrac(L) &
                     + Ea1t99(Shadeleaftk(L),T24,T240,SPEC_NUM) * max(0._r8, 1._r8 - Sunfrac(L))

! pstd = 200 for sun leaves
! pstd = 50 for shade leaves
                Ea1pLayer(L) = &
                     Ea1p99(SunPPFD(L),PPFD24*0.5_r8,PPFD240*0.5_r8,pstd_Sun) &
                        * Sunfrac(L) &
                     + Ea1p99(ShadePPFD(L),PPFD24*0.16_r8,PPFD240*0.16_r8,pstd_Shade) &
                        * max(0._r8, 1._r8 - Sunfrac(L))

                Ea1Layer(L)  = &
                     Ea1t99(Sunleaftk(L), T24, T240,SPEC_NUM) &
                       * Ea1p99(SunPPFD(L),PPFD24*0.5_r8, PPFD240*0.5,pstd_Sun) &
                         * Sunfrac(L) &
                     + Ea1t99(Shadeleaftk(L), T24, T240,SPEC_NUM) & 
                       * Ea1p99(ShadePPFD(L), PPFD24*0.16_r8, PPFD240*0.16_r8, pstd_Shade) &
                         * max(0._r8, 1._r8 - Sunfrac(L))

                EatiLayer(L) = &
                     Ealti99(SPEC_NUM, Sunleaftk(L)) * Sunfrac(L)  &
                     + Ealti99(SPEC_NUM, Shadeleaftk(L)) * max(0._r8, 1._r8 - Sunfrac(L))

             end do !// do L = 1, layers

             Ea1pCanopy = SUM( Ea1pLayer(:) * VPslwWT(:) * VPgausWt(:) )
             Ea1tCanopy = SUM( Ea1tLayer(:) * VPslwWT(:) * VPgausWt(:) )
             Ea1Canopy  = SUM( Ea1Layer(:)  * VPslwWT(:) * VPgausWt(:) )
             EatiCanopy = SUM( EatiLayer(:) * VPslwWT(:) * VPgausWt(:) )
             Ea1Canopy = Ea1Canopy * Cce * LAI

             !// Sum up fractional contributions
             ADJUST_FACTOR_LD = ADJUST_FACTOR_LD + &
                  canopyFraction(I_PFT) * Ea1Canopy
             ADJUST_FACTOR_LI = ADJUST_FACTOR_LI + &
                  canopyFraction(I_PFT) * EatiCanopy
             if (PPFD0 .gt. 0._r8 .and. (Ea1Canopy+EatiCanopy) .eq. 0._r8) then
                write(6,'(a,2i4,5es12.3)') f90file//':'//subr// &
                     'Have light but zero TD',ILOC,JLOC,PPFD,PPFD24,Ea1Canopy,EatiCanopy, LAI
                stop
             end if
          end if
       end do !// do I_PFT = 1, N_MGN_PFT

       !// Set average GAMMA values
       GAM_TLD = ADJUST_FACTOR_LD / CanopyFractionTotal
       GAM_TLI = ADJUST_FACTOR_LI / CanopyFractionTotal

    else if (CanopyFractionTotal .eq. 0._r8) then
       !// If no canopy, GAMMAs could as well be zero
       GAM_TLD = 1._r8
       GAM_TLI = 1._r8
    else if (CanopyFractionTotal .lt. 0._r8) then
       write(6,'(a,2i5,es16.6)') f90file//':'//subr// &
            'CanopyFractionTotal is less than 0. I,J: ',ILOC,JLOC,CanopyFractionTotal
       stop
    end if



    !// --------------------------------------------------------------------
  end subroutine  GAMMA_CANOPY2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine GAMMA_CANOPY(JDAY, SZArad, Tc, T24, T240, &
                          PPFD0, PPFD24, PPFD240, Wind, &
                          Humidity,Cantype,LAI,Pres,DI, &
                          SPEC_NUM, SPCNAME, Ea1Canopy, EatiCanopy)
    !// --------------------------------------------------------------------
    !// Modified GAMME_CE from MEGANv2.10 code.
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, to fit Oslo CTM3
    !// --------------------------------------------------------------------
    use cmn_parameters, only: M_AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: JDAY, Cantype, SPEC_NUM
    real(r8), intent(in) :: &
         SZArad, &   ! Solar zenith angle in radians
         Tc, &       ! Surface temperature [K]
         Pres, &     ! Surface pressure [hPa]
         Wind, &     ! Surface wind [m/s]
         Humidity, & ! Humidity, supposedly mixing ratio [kg(H2O)/kg(dry air)]
                     ! but using specific humidity [kg(H2O)/kg(moist air)].
         LAI, &      ! Leaf area index
         DI, &       ! Drought index (PSDI)
         PPFD0, &    ! PPFD current [W/m2]
         PPFD24, &   ! PPFD averaged over last 24 hours [W/m2]
         PPFD240, &  ! PPFD averaged over last 10 days [W/m2]
         T24, &      ! Surface temperature averaged over last 24 hours [K]
         T240        ! Surface temperature averaged over last 10 days [K]
    character(len=*), intent(in) :: SPCNAME

    !// Output
    real(r8), intent(out) :: &
         Ea1Canopy, & ! Light dependent
         EatiCanopy   ! Light independent

    ! Local variables
    integer :: i, JJ, KK, K

    real(r8) :: Sinbeta, Beta

    !// Number of canopy layers
    integer, parameter :: Layers = 5

    real(r8), parameter :: &
         !ConvertWm2toUmolm2s = 4.766_r8, & ! Convert radiation from [W/m2] to [umol/m2/s1]
         SolarConstant = 1367._r8, &        ! Solar constant [W/m2]
         WaterAirRatio = 18.016_r8 / M_AIR, & ! Ratio between water and air molecules
         pstd_Sun = 200.0_r8, &
         pstd_Shade = 50.0_r8

    real(r8), dimension(Layers) ::  VPgausWt, VPgausDis2, &
         VPgausDis, VPslwWT, Sunfrac, QdAbsV, QsAbsV, QdAbsn, &
         QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD, &
         ShadePPFD, TairK, HumidairPa, Ws, Sunleaftk, SunleafSH, &
         SunleafLH,SunleafIR, Shadeleaftk, ShadeleafSH, &
         ShadeleafLH,ShadeleafIR,Ea1pLayer, Ea1tLayer, Ea1Layer, &
         EatiLayer

    real(r8) :: Solar, Maxsolar, PPFD, QbAbsn, &
         Trate, StomataDI, Qbeamv,Qdiffv, Qbeamn, Qdiffn, &
         QbAbsV,Ea1tCanopy, Ea1pCanopy, &
         TairK0, HumidairPa0,Ws0
                       
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'GAMMA_CANOPY'
    !// --------------------------------------------------------------------

    !// Need solar elevation angle (beta), which is related to solar zenith
    !// angle (SZA) as  sin(beta) = cos(SZA):
    sinBETA = cos(SZArad)

    TairK0     = Tc   !See the input unit(K)
    Ws0        = Wind 
    !// CanopyRad will convert from W/m2 to umol/m2/s using different
    !// value for sunlit and shaded leaves. Thus, we have to assume
    !// PPFD should have units W/m2.
    PPFD       = PPFD0 
    StomataDI  = DIstomata(DI)
    !// MEGAN code did 
    !//    Solar      = PPFD / 2.25_r8
    !// and previously had commented out
    !//    !Solar      = PPFD /ConvertWm2toUmolm2s*2
    !//    !Solar      = PPFD/2.1
    !// which could be a way to convert umol/m2/s to W/m2 which can later
    !// be converted back to umol/m2/s in CanopyRad.
    Solar      = PPFD
    Maxsolar   = Sinbeta  * SolarConstant * CalcEccentricity(JDAY)

    call GaussianIntegration(VPgausWt, VPgausDis, VPgausDis2, Layers)

    call SolarFractions(Solar, Maxsolar, Qdiffv, Qbeamv, Qdiffn, Qbeamn)

    call WeightSLW(VPgausDis, VPgausWt, LAI, Layers, VPslwWT)

    call CanopyRad(VPgausDis, Layers, LAI, Sinbeta, Qbeamv, &
         Qdiffv, Qbeamn, Qdiffn,Cantype, Sunfrac, &
         QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv, &
         ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD)


    HumidairPa0  =  WaterVapPres(Humidity , Pres , WaterAirRatio)

    Trate    =  Stability(Cantype, Solar)


    call CanopyEB(Trate, Layers, VPgausDis, Cantype, &
         StomataDI, TairK, HumidairPa, Ws, SunPPFD, &
         ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, &
         Sunleaftk, SunleafSH, SunleafLH, SunleafIR, &
         Shadeleaftk,ShadeleafSH,ShadeleafLH,ShadeleafIR, &
         Ws0, TairK0, HumidairPa0)


    do i = 1, layers    

       Ea1tLayer(i) = Ea1t99(Sunleaftk(i), T24, T240,SPEC_NUM) * Sunfrac(i) &
                    + Ea1t99(Shadeleaftk(i),T24,T240,SPEC_NUM) * max(0._r8, 1._r8 - Sunfrac(i))

! pstd = 200 for sun leaves
! pstd = 50 for shade leaves
       Ea1pLayer(i) = Ea1p99(SunPPFD(i),PPFD24*0.5_r8,PPFD240*0.5_r8,pstd_Sun) &
                            * Sunfrac(i) &
                    + Ea1p99(ShadePPFD(i),PPFD24*0.16_r8,PPFD240*0.16_r8,pstd_Shade) &
                        * max(0._r8, 1._r8 - Sunfrac(i))

       Ea1Layer(i)  = Ea1t99(Sunleaftk(i), T24, T240,SPEC_NUM) &
                      * Ea1p99(SunPPFD(i),PPFD24*0.5_r8, PPFD240*0.5,pstd_Sun) &
                        * Sunfrac(i) &
                    + Ea1t99(Shadeleaftk(i), T24, T240,SPEC_NUM) & 
                      * Ea1p99(ShadePPFD(i), PPFD24*0.16_r8, PPFD240*0.16_r8, pstd_Shade) &
                        * max(0._r8, 1._r8 - Sunfrac(i))

       EatiLayer(i) = Ealti99(SPEC_NUM, Sunleaftk(i)) * Sunfrac(i)  &
                    + Ealti99(SPEC_NUM, Shadeleaftk(i)) * max(0._r8, 1._r8 - Sunfrac(i))

    end do

    Ea1pCanopy = SUM( Ea1pLayer(:) * VPslwWT(:) * VPgausWt(:) )
    Ea1tCanopy = SUM( Ea1tLayer(:) * VPslwWT(:) * VPgausWt(:) )
    Ea1Canopy  = SUM( Ea1Layer(:)  * VPslwWT(:) * VPgausWt(:) )
    EatiCanopy = SUM( EatiLayer(:) * VPslwWT(:) * VPgausWt(:) )
    Ea1Canopy = Ea1Canopy * Cce* LAI

! this quantity is apparently not passed out of the subroutine
!      SH = SUM( (SunleafSH(:) * Sunfrac(:) +                 
!     &           ShadeleafSH(:) * (1 - Sunfrac(:))) *       
!     &           LAI  * VPgausWt(:)                   )

    !// --------------------------------------------------------------------
  end subroutine  GAMMA_CANOPY
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function DIstomata(DI)
    !// --------------------------------------------------------------------
    !// Calculate drought-induced effect on stomata (range 0-1).
    !// Input DI is the Palmer Severity Drought Index.
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: DI
    ! > -.5 incipient,  mild or no drought; < -4 extreme drought
    real(r8), parameter :: DIhigh = -0.5_r8, DIlow = -5._r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'DIstomata'
    !// --------------------------------------------------------------------

    if (DI .gt. DIhigh) then
       DIstomata = 1._r8  ! no drought
    else if (DI .gt. DIlow) then
       ! interpolate
       DIstomata = 1._r8 - (0.9_r8 * ((DI - DIhigh) / (DIlow - DIhigh))) 
    else
       DIstomata = 0._r8  ! Maximum drought, maximum stomatal reissstance
    end if

    !// --------------------------------------------------------------------
  end function DIstomata
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function Ea1t99(T1, T24, T240, SPEC_NUM)
    !// --------------------------------------------------------------------
    !// Temperature dependence activity factor for emission type 1 
    !//        (e.g. isoprene, MBO)
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), parameter :: Ctm2 = 230._r8, Ts  = 297._r8
    real(r8), intent(in) :: T1, T24, T240
    integer, intent(in) :: SPEC_NUM
    real(r8) :: Topt, X, Eopt
    !// --------------------------------------------------------------------

    if (T1 .lt. 260._r8) then
       Ea1t99 = 0._r8
    else
       ! Energy of activation and deactivation

       ! Temperature at which maximum emission occurs
       Topt = 312.5_r8 + 0.6_r8 * (T240 - Ts)
       X    = ((1._r8 / Topt) - (1._r8 / T1)) / 0.00831_r8

       ! Maximum emission (relative to emission at 30 C)
       Eopt = CLeo(SPEC_NUM) * exp(0.05_r8 * (T24 - Ts)) &
             * exp(0.05_r8 * (T240 - Ts))          

       Ea1t99 = Eopt * Ctm2 * exp(Ctm1(SPEC_NUM) * X) &
            / (Ctm2 - Ctm1(SPEC_NUM) * (1._r8 - exp(Ctm2 * X)))

    end if

    !// --------------------------------------------------------------------
  end function Ea1t99
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  real(r8) function Ealti99(SPEC_NUM, temp)
    !// --------------------------------------------------------------------
    ! calculate light indepent algorithms
    ! coded by Xuemei Wang 05 Nov. 2007
    !--   GAMMA_TLI =  exp[BETA*(T-Ts)]
    !           where BETA   = temperature dependent parameter
    !                 Ts     = standard temperature (normally 303K, 30C)
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: temp
    integer, intent(in) :: SPEC_NUM
    !// --------------------------------------------------------------------
    real(r8) ,parameter :: Ts = 303.15_r8
    !// --------------------------------------------------------------------
    Ealti99 = exp( TDF_PRM(SPEC_NUM) * (temp - Ts) )
    !// --------------------------------------------------------------------
  end function Ealti99
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function Ea1p99(PPFD1, PPFD24, PPFD240, PSTD)
    !// --------------------------------------------------------------------
    !// pstd = 200 for sun leaves and 50 for shade leaves 
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: PPFD1, PPFD24, PPFD240, PSTD
    real(r8) :: Alpha, C1, a2
    !// --------------------------------------------------------------------

    if (PPFD240 .lt. 0.01_r8) then
       Ea1p99 = 0._r8
    else
       Alpha  = 0.004_r8 - 0.0005_r8 * LOG(PPFD240)
       C1     = 0.0468_r8 * EXP(0.0005_r8 * (PPFD24 - PSTD)) &
                * (PPFD240 ** 0.6_r8)
       a2 = (Alpha * PPFD1)**2
       !Ea1p99 = (Alpha * C1 * PPFD1) &
       !     / ((1 + Alpha**2 * PPFD1**2.)**0.5_r8)
       Ea1p99 = A2 * C1 / ( (1._r8 + A2)**0.5_r8 )
    end if

    !// --------------------------------------------------------------------
  end function  Ea1p99
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function WaterVapPres(wmr, Pres, WaterAirRatio)
    !// --------------------------------------------------------------------
    !// Convert water mixing ratio (kg/kg) to water vapor pressure 
    !// Pa or Kpa depending on units of input )
    !// Mixing ratio wmr (kg/kg), pressure (Pa)
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: wmr, Pres, WaterAirRatio
    !// --------------------------------------------------------------------
    WaterVapPres = (wmr / (wmr + WaterAirRatio)) * Pres
    !// --------------------------------------------------------------------
  end function WaterVapPres
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  real(r8) function Stability(Cantype, Solar)
    !// --------------------------------------------------------------------
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: Cantype
    real(r8), intent(in) :: Solar
    real(r8), parameter :: Trateboundary = 500._r8
    !// --------------------------------------------------------------------

    if (Solar .gt. Trateboundary) then
       ! Daytime temperature lapse rate
       Stability = Canopychar(12, Cantype)
    else if (Solar .gt. 0._r8) then
       Stability = Canopychar(12, Cantype) &
            - ( (Trateboundary - Solar) / Trateboundary) &
                * (Canopychar(12, Cantype) - Canopychar(13, Cantype))
    else
       ! Nightime temperature lapse rate
       Stability = Canopychar(13, Cantype)
    end if

    !// --------------------------------------------------------------------
  end function Stability
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine GaussianIntegration(Weightgauss, Distgauss,Distgauss2, Layers)
    !// --------------------------------------------------------------------
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: Layers
    real(r8), dimension(Layers), intent(out) :: &
         Weightgauss, Distgauss, Distgauss2
    !/ local variables
    integer ::  i
    !// --------------------------------------------------------------------

    if (Layers .eq. 1) then
       Weightgauss(1) = 1._r8
       Distgauss(1)   = 0.5_r8
       Distgauss2(1)   = 1._r8
    else if (Layers .eq. 3) then
       Weightgauss(1) = 0.277778_r8
       Weightgauss(2) = 0.444444_r8
       Weightgauss(3) = 0.277778_r8
       Distgauss(1)   = 0.112702_r8
       Distgauss(2)   = 0.5_r8
       Distgauss(3)   = 0.887298_r8
       Distgauss2(1)   = 0.277778_r8
       Distgauss2(2)   = 0.722222_r8
       Distgauss2(3)   = 1._r8
    else if (Layers .eq. 5) then
       Weightgauss(1) = 0.1184635_r8
       Weightgauss(2) = 0.2393144_r8
       Weightgauss(3) = 0.284444444_r8
       Weightgauss(4) = 0.2393144_r8
       Weightgauss(5) = 0.1184635_r8
       Distgauss(1)   = 0.0469101_r8
       Distgauss(2)   = 0.2307534_r8
       Distgauss(3)   = 0.5_r8
       Distgauss(4)   = 0.7692465_r8
       Distgauss(5)   = 0.9530899_r8
       Distgauss2(1)   = 0.1184635_r8
       Distgauss2(2)   = 0.3577778_r8
       Distgauss2(3)   = 0.6422222_r8
       Distgauss2(4)   = 0.881536_r8
       Distgauss2(5)   = 1.0_r8
    else
       do i = 1, Layers
          Weightgauss(i) = 1._r8 / Layers
          Distgauss(i)   = (i - 0.5_r8) / Layers
          Distgauss2(i) = i/Layers
       end do
    end if

    !// --------------------------------------------------------------------
  end subroutine GaussianIntegration
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine SolarFractions( Solar, Maxsolar, Qdiffv,Qbeamv, Qdiffn, Qbeamn)
    !// --------------------------------------------------------------------
    !// Transmission, fraction of PPFD that is diffuse, 
    !// fraction of solar rad that is PPFD
    !// From MEGANv2.10 code, checked ok against MEGAN3.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Solar, Maxsolar
    real(r8), intent(out) :: Qdiffv,Qbeamv, Qdiffn, Qbeamn

    real(r8) :: FracDiff, PPFDfrac,PPFDdifFrac,Qv, Qn
    integer :: I,J
    real(r8) ::  Transmis 
    !// --------------------------------------------------------------------
         
    if (Maxsolar .le. 0._r8) then
       Transmis  = 0.5_r8
    else if (Maxsolar .lt. Solar) then
       Transmis  = 1.0_r8
    else
       Transmis  = Solar  / Maxsolar 
    end if

    !FracDiff is based on Lizaso 2005
    FracDiff = 0.156_r8 + 0.86_r8/(1._r8 + EXP(11.1_r8 * (Transmis - 0.53_r8)))

    !PPFDfrac is based on Goudrian and Laar 1994
    PPFDfrac  = 0.55_r8 -Transmis * 0.12_r8

    !PPFDdifFrac is based on data in Jacovides 2007
    PPFDdifFrac = FracDiff * (1.06_r8 + Transmis*0.4_r8)

    ! Calculte  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
    if (PPFDdifFrac .gt. 1._r8) then
       PPFDdifFrac  = 1._r8
    end if
 
    Qv  = PPFDfrac * Solar
    Qdiffv = Qv * PPFDdifFrac
    Qbeamv = Qv - Qdiffv
    Qn = Solar - Qv
    Qdiffn =  Qn * FracDiff
    Qbeamn =  Qn - Qdiffn

    !// --------------------------------------------------------------------
  end subroutine SolarFractions
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine WeightSLW(Distgauss, Weightgauss, LAI, Layers, SLW)
    !// --------------------------------------------------------------------
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: Layers
    real(r8), intent(in) :: LAI
    real(r8), dimension(Layers), intent(in) :: Distgauss, Weightgauss
    real(r8), dimension(Layers), intent(out) :: SLW

    !// local variables
    integer :: i
    !// --------------------------------------------------------------------

    do i = 1, Layers
       SLW(i) = 0.63_r8 + 0.37_r8 * EXP(-((LAI * Distgauss(i)) - 1._r8))
    end do

    !// --------------------------------------------------------------------
  end subroutine WeightSLW
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine CanopyRad(Distgauss, Layers, LAI, Sinbeta, &
       Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype, &
       Sunfrac, QbAbsV, QdAbsV, QsAbsV, &
       QbAbsn, QdAbsn, QsAbsn, SunQv, &
       ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD)
    !// --------------------------------------------------------------------
    !   Canopy light environment model
    !   Code developed by Alex Guenther, based on Spitters et al. (1986), 
    !   Goudrian and Laar (1994), Leuning (1997)
    !   Initial code 8-99, modified 7-2000 and 12-2001
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    ! input
    integer, intent(in) :: Layers, Cantype
    real(r8), intent(in) :: Qbeamv,Qdiffv,Sinbeta,LAI,Qbeamn,Qdiffn
    real(r8), dimension(Layers), intent(in) :: Distgauss

    ! output
    real(r8), intent(out) :: QbAbsV, QbAbsn
    real(r8), dimension(Layers),intent(out) :: ShadePPFD, SunPPFD, &
         QdAbsv, QsAbsv, QsAbsn, ShadeQv,  SunQn, &
         QdAbsn, SunQv, ShadeQn, Sunfrac

    ! internal variables
    integer :: i
    real(r8) :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN, &
         Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, Cluster, &
         QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL

    !REAL,PARAMETER :: ConvertPPFD = 4.766
    real(r8), parameter :: ConvertShadePPFD = 4.6_r8
    real(r8), parameter :: ConvertSunPPFD = 4.0_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CanopyRad'
    !// --------------------------------------------------------------------

    !// MEGAN3 checks LAIadj = LAI/(1 - canopy transparency) instead of LAI.
    !// MEGAN2.1 has only 16 canopy chars, but MEGAN3 includes also
    !// canopy transparency.
    !// adjust LAI for canopy transparency
    !if (LAI .eq. 0._r8) then
    !   LAI_adj = 0._r8
    !else
    !   !// In case you have set transparency to 1, avoid dividing by zero.
    !   LAIadj = LAI / max(1.e-5, 1._r8 - Canopychar(17,Cantype))
    !end if

    if (((Qbeamv  + Qdiffv ) .gt. 0.001_r8) .and. &
         (Sinbeta  .gt. 0.00002_r8) .AND. &
         (LAI .gt. 0.001_r8)) THEN       ! Daytime

       ! Scattering coefficients (scatV,scatN), diffuse and beam reflection 
       ! coefficients (ref..) for visible or near IR
       ScatV   = Canopychar(5,Cantype)
       ScatN   = Canopychar(6,Cantype)
       RefldV  = Canopychar(7,Cantype)
       RefldN  = Canopychar(8,Cantype)
       Cluster = Canopychar(9,Cantype)
       !        print*,'cluster',  Cluster
       ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
       Kb = Cluster * 0.5_r8 / Sinbeta 
       ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
       Kd = 0.8_r8 * Cluster              
       ! (0.8 assumes a spherical leaf angle distribution)

       Call CalcExtCoeff(Qbeamv,ScatV,Kb,Kd,ReflbV,KbpV,KdpV,QbAbsV)
       Call CalcExtCoeff(Qbeamn,ScatN,Kb,Kd,ReflbN,KbpN,KdpN,QbAbsn)

       do i = 1,layers
          ! LAI depth at this layer
          LAIdepth   = LAI  * Distgauss(i)  
          !fraction of leaves that are sunlit
          Sunfrac(i) = EXP(-Kb * LAIdepth)  

          Call CalcRadComponents(Qdiffv , Qbeamv , kdpV, &
                                kbpV, kb, scatV, refldV, &
                                reflbV, LAIdepth, QdAbsVL, QsAbsVL)

          Call CalcRadComponents(Qdiffn , Qbeamn , kdpN, &
                                kbpN, kb, scatN, refldN, &
                                reflbN, LAIdepth, QdAbsNL, QsAbsNL)
          !// Fail-safe
          if (scatV .eq. 1._r8) then
             write(6,'(a,i3,2es16.6)') f90file//':'//subr// &
                  ': scatV == 1. - stopping: ',i,LAIdepth, Sunfrac(i)
             stop
          end if
          ShadePPFD(i) = (QdAbsVL + QsAbsVL) * ConvertShadePPFD / (1._r8 - scatV)
          SunPPFD(i) = ShadePPFD(i) + (QbAbsV * ConvertSunPPFD / (1._r8 - scatV))
          QdAbsV(i) = QdAbsVL
          QsAbsV(i) = QsAbsVL
          QdAbsn(i) = QdAbsNL
          QsAbsn(i) = QsAbsNL
          ShadeQv(i) = QdAbsVL + QsAbsVL
          SunQv(i)   = ShadeQv(i) + QbAbsV 
          ShadeQn(i) = QdAbsNL + QsAbsNL
          SunQn(i)   = ShadeQn(i) + QbAbsn 
       end do

    else                           ! Night time

       QbAbsV  = 0._r8
       QbAbsn   = 0._r8

       Sunfrac(:)   = 0.2_r8
       SunQn(:)     = 0._r8
       ShadeQn(:)   = 0._r8
       SunQv(:)     = 0._r8
       ShadeQv(:)   = 0._r8
       SunPPFD(:)   = 0._r8
       ShadePPFD(:) = 0._r8
       QdAbsV(:)    = 0._r8
       QsAbsV(:)    = 0._r8
       QdAbsn(:)    = 0._r8
       QsAbsn(:)    = 0._r8

    end if

    !// --------------------------------------------------------------------
  end subroutine CanopyRad
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine CalcExtCoeff(Qbeam,scat,kb,kd,reflb,kbp,kdp,QbeamAbsorb)
    !// --------------------------------------------------------------------
    !// Calculate light extinction coefficients
    !// Code originally developed by Alex Guenther in 1990s
    !// From MEGANv2.10 code, checked ok against MEGAN3.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Qbeam, scat, Kb, Kd
    real(r8), intent(out) :: Reflb, Kbp, Kdp, QbeamAbsorb
    ! local variables
    real(r8) :: P
    !// --------------------------------------------------------------------
     
    P     = (1._r8 - scat)**0.5_r8
    Reflb = 1._r8 - exp((-2_r8 * ((1._r8 - P) / (1._r8 + P)) * kb) / (1._r8 + kb))

    ! Extinction coefficients
    Kbp   = Kb * P
    Kdp   = Kd * P
    ! Absorbed beam radiation
    QbeamAbsorb = kb * Qbeam * (1._r8 - scat) 

    !// --------------------------------------------------------------------
  end subroutine CalcExtCoeff
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, &
       scat, refld, reflb, LAIdepth, QdAbs, QsAbs)
    !// --------------------------------------------------------------------
    !// From MEGANv2.10 code, checked OK against MEGAN3.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Qdiff,Qbeam,kdp,kbp,kb,scat, &
                        refld,reflb,LAIdepth
    real(r8), intent(out) :: QdAbs, QsAbs
    !// --------------------------------------------------------------------

    QdAbs = Qdiff *    Kdp * (1._r8 - Refld) * exp(-Kdp * LAIdepth)
    QsAbs = Qbeam * ( (Kbp * (1._r8 - Reflb) * exp(-Kbp * LAIdepth)) &
                     - (Kb * (1._r8 - Scat ) * exp(-Kb * LAIdepth)) )

    !// --------------------------------------------------------------------
  end subroutine CalcRadComponents
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine CanopyEB(Trate, Layers, Distgauss, &
                  Cantype, StomataDI, TairK, HumidairPa, Ws, &
                  SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, &
                  Sunleaftk, SunleafSH, SunleafLH, &
                  SunleafIR, Shadeleaftk, ShadeleafSH, &
                  ShadeleafLH, ShadeleafIR, Ws0, &
                  TairK0, HumidairPa0)
    !// --------------------------------------------------------------------
    !// Canopy energy balance model for estimating leaf temperature
    !// Code developed by Alex Guenther, based on Goudrian and Laar (1994),
    !// and Leuning (1997).
    !// Initial code 8-99, modified 7-2000 and 12-2001
    !//
    !// Note: i denotes an array containing a vertical profile through the 
    !//       canopy with 0 
    !//       (above canopy conditions) plus 1 to number of canopy layers
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    ! inputs
    integer, intent(in) :: Layers, Cantype
    real(r8), intent(in) :: Trate, StomataDI, TairK0, HumidairPa0, Ws0
    real(r8), dimension(Layers), intent(in) :: Distgauss, SunQv,ShadeQv, &
         SunQn, ShadeQn, SunPPFD, ShadePPFD

    ! outputs
    real(r8), dimension(Layers), intent(out) :: HumidairPa, &
         Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, &
         Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

    ! local variables
    integer :: i
    real(r8) :: Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType, &
         Deltah, Rin,IRout, IRin
    real(r8), dimension(Layers) :: Ldepth, Wsh
    !// --------------------------------------------------------------------

         
    Cdepth        = Canopychar(1, Cantype)
    Lwidth        = Canopychar(2, Cantype)
    Llength       = Canopychar(3, Cantype)
    Cheight       = Canopychar(4, Cantype)
    Eps           = Canopychar(10,Cantype)
    TranspireType = Canopychar(11,Cantype)

    if (TairK0 .gt. 288._r8) then
       ! Pa m-1  (humidity profile for T < 288)
       Deltah =  Canopychar(14, Cantype) / Cheight
    else if (TairK0  .gt. 278._r8) then
       Deltah =( Canopychar(14,Cantype) - ((288._r8 - TairK0)/10._r8) &
             * (Canopychar(14,Cantype) - Canopychar(15,Cantype)) ) / Cheight
    else
       ! Pa m-1  (humidity profile for T <278)
       Deltah = Canopychar(15, Cantype) / Cheight
    end if

    Ldepth(:)     = Cdepth * Distgauss(:)
    TairK(:)      = TairK0  + (Trate  * Ldepth(:))      ! check this
    HumidairPa(:) = HumidairPa0  + (Deltah * Ldepth(:))

    Wsh(:) = (Cheight-Ldepth(:)) - (Canopychar(16,Cantype) * Cheight)
    Ws(:)  = (Ws0 * log(Wsh(:)) / log(Cheight - Canopychar(16,Cantype) * Cheight))
    where (Wsh(:) .lt. 0.001_r8) Ws(:) = 0.05_r8

    do i = 1, Layers
       IRin     = UnexposedLeafIRin(TairK(i), Eps)
       ShadeleafIR(i) = 2._r8 * IRin
       SunleafIR(i) = 0.5_r8 * ExposedLeafIRin(HumidairPa0,TairK0) + 1.5_r8 * IRin

       ! Sun
       call LeafEB(SunPPFD(i), SunQv(i) + SunQn(i), &
                    SunleafIR(i), Eps, TranspireType, Lwidth, Llength, &
                    TairK(i), HumidairPa(i), Ws(i), &
                    Sunleaftk(i), SunleafSH(i),SunleafLH(i), &
                    IRout,StomataDI )

       SunleafIR(i) = SunleafIR(i) - IRout

       ! Shade
       call LeafEB(ShadePPFD(i), ShadeQv(i)+ShadeQn(i), &
                    ShadeleafIR(i),Eps,TranspireType, Lwidth,Llength, &
                    TairK(i), HumidairPa(i), Ws(i), &
                    Shadeleaftk(i), ShadeleafSH(i),ShadeleafLH(i), &
                    IRout, StomataDI )

       ShadeleafIR(i) = ShadeleafIR(i) - IRout
    end do

    !// --------------------------------------------------------------------
  end subroutine CanopyEB
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine LeafEB(PPFD, Q, IRin, Eps, TranspireType, &
       Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf, &
       SH, LH, IRout, StomataDI)
    !// --------------------------------------------------------------------
    !// Leaf energy balance
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Eps, TranspireType, Lwidth, Llength,PPFD, Q, &
         IRin, TairK, HumidairPa, StomataDI, Ws
    real(r8), intent(out) :: IRout, Tleaf, SH, LH

    ! local variables
    integer :: i
    real(r8) :: HumidAirKgm3,GHforced,StomRes,IRoutairT,LatHv,Ws1, &
         LHairT,Tdelt,Balance, &
         GH1,SH1,LH1,E1,IRout1,GH
    !// --------------------------------------------------------------------

    if (Ws .le. 0._r8) then
       Ws1 = 0.001_r8
    else
       Ws1 = Ws
    end if

    ! Air vapor density kg m-3
    HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

    ! Heat convection coefficient (W m-2 K-1) for forced convection. 
    ! Nobel page 366
    GHforced = 0.0259_r8 / (0.004_r8 * ((Llength / Ws)**0.5_r8))

    ! Stomatal resistence s m-1
    StomRes  = ResSC(PPFD, stomataDI)

    IRoutairT = LeafIROut(tairK, eps)

    ! Latent heat of vaporization (J Kg-1)
    LatHv = LHV(TairK)

    ! Latent heat flux
    LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes,TranspireType)

    E1 = (Q + IRin - IRoutairT - LHairT)
    if (E1 .eq. 0._r8) then
       E1 = -1._r8
    end if

    Tdelt = 1._r8
    Balance = 10._r8
    do i = 1, 10
       if (abs(Balance) .gt. 2._r8) then
          ! Boundary layer conductance
          GH1 = LeafBLC(GHforced, Tdelt, Llength)
          ! Convective heat flux
          SH1 = LeafH(Tdelt, GH1)                
          ! Latent heat of vaporization (J Kg-1)
          LatHv = LHV(TairK + Tdelt)             
          LH = LeafLE(TairK + Tdelt, HumidAirKgm3, LatHv, GH1, StomRes, TranspireType)
          LH1 = LH - LHairT
          IRout  = LeafIROut(TairK + Tdelt, Eps)
          IRout1 = IRout - IRoutairT
          Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
          Balance = Q + IRin - IRout - SH1 - LH
       end if
    end do

    if (Tdelt .gt. 10)  Tdelt = 10._r8
    if (Tdelt .lt. -10) Tdelt = -10._r8

    Tleaf = TairK + Tdelt
    GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
    SH    = LeafH(Tleaf - TairK, GH)
    LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, GH, StomRes, TranspireType)
    IRout = LeafIROut(Tleaf, Eps)

    !// --------------------------------------------------------------------
  end subroutine LeafEB
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function ConvertHumidityPa2kgm3(Pa, Tk)
    !// --------------------------------------------------------------------
    !   Saturation vapor density  (kg/m3)
    !// From MEGANv2.10 code.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Pa, Tk
    !// --------------------------------------------------------------------
    ConvertHumidityPa2kgm3 = 0.002165_r8 * Pa / Tk
    !// --------------------------------------------------------------------
  end function ConvertHumidityPa2kgm3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function ResSC(Par, StomataDI)
    !// --------------------------------------------------------------------
    !// Leaf stomatal resistance [s m-1]
    !//
    !// From MEGANv2.10 code. MEGAN3 has removed StomataDI.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Par, StomataDI
    real(r8) :: SCadj
    !// --------------------------------------------------------------------

    !// StomataDI is multiplied by light dependence factor from
    !// Guenther et al. (1993, doi:10.1029/93JD00527)
    !// It is not explained why the G1993 equation can be used, but it seems
    !// a scaling factor important for sunrise/-set.
    !//
    !// PAR has units umol/m2/s, which is also often the case
    !// for stomatal conductance (also often converted to m/s).
    !// Resistance is the inverse of conductance.
    SCadj = StomataDI * ((0.0027_r8 * 1.066_r8 * Par) &
         / ((1._r8 + 0.0027_r8 * 0.0027_r8 * Par**2)**0.5_r8))

    if (SCadj .lt. 0.1_r8) then
       ResSC = 2000._r8
    else
       ResSC = 200._r8 / SCadj
    end if

    !// --------------------------------------------------------------------
  end function ResSC
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function LeafIROut(Tleaf, Eps)
    !// --------------------------------------------------------------------
    !// IR thermal radiation energy output by leaf
    !// Calculate IR transfer between leaf and air
    !// From MEGANv2.10 code. Called LeafIR in MEGAN3.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Tleaf, Eps
    ! Stefan-boltzman constant  W m-2 K-4
    real(r8), parameter :: Sb = 5.670367e-8_r8
    !// --------------------------------------------------------------------

    LeafIROut = Eps * Sb * 2_r8 * (Tleaf**4)

    !// --------------------------------------------------------------------
  end function LeafIROut
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function LHV(Tk)
    !// --------------------------------------------------------------------
    !// Latent Heat of vaporization(J Kg-1) from Stull p641
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    use cmn_parameters, only: Lv_0C, TK_0C
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Tk ! [K]
    !// --------------------------------------------------------------------

    LHV = Lv_0C - (2370._r8 * (Tk - TK_0C))

    !// --------------------------------------------------------------------
  end function LHV
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)
    !// --------------------------------------------------------------------
    !// Latent energy term in Energy balance
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: &
         Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType
    real(r8) :: LeafRes, Vapdeficit, LE
    !// --------------------------------------------------------------------

    LeafRes    = (1._r8 / (1.075_r8 * (GH / 1231._r8))) + StomRes
    Vapdeficit = (SvdTk(Tleaf) - Ambvap)

    ! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
    !                 leaf resistence (s m-1)
    LE = TranspireType * (1._r8 / LeafRes) * LatHv * Vapdeficit
    if (LE .lt. 0) then
       LeafLE = 0._r8
    else
       LeafLE = LE
    end if

    !// --------------------------------------------------------------------
  end function LeafLE
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function LeafBLC(GHforced, Tdelta, Llength)
    !// --------------------------------------------------------------------
    !// Boundary layer conductance
    !// This is based on Leuning 1995 p.1198 except using molecular 
    !// conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
    !// diffusivity so that you end up with a heat convection coefficient 
    !// (W m-2 K-1) instead of a conductance for free convection
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: GHforced, Tdelta, Llength
    real(r8) :: Ghfree
    !// --------------------------------------------------------------------

    if (Tdelta .ge. 0._r8) then
       GhFree = 0.5_r8 * 0.00253_r8 * ((160000000._r8 * Tdelta &
            / (Llength**3))**0.25_r8) / Llength
    else
       GhFree = 0._r8
    end if
    LeafBLC = GHforced + GhFree

    !// --------------------------------------------------------------------
  end function LeafBLC
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function LeafH(Tdelta, GH)
    !// --------------------------------------------------------------------
    !// Convective energy term in Energy balance (W m-2 heat flux from 
    !//    both sides of leaf)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Tdelta, GH
    !// --------------------------------------------------------------------
    ! 2 sides X conductance X Temperature gradient
    LeafH = 2._r8 * GH * Tdelta
    !// --------------------------------------------------------------------
  end function LeafH
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function SvdTk(Tk)
    !// --------------------------------------------------------------------
    !// Saturation vapor density  (kg/m3)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8) :: Tk, Svp
    !// --------------------------------------------------------------------

    ! Saturation vapor pressure (millibars)
    Svp = 10._r8**((-2937.4_r8 / Tk) - (4.9283_r8 * LOG10(Tk)) + 23.5518_r8)
    SvdTk = 0.2165_r8 * Svp / Tk
    !// --------------------------------------------------------------------
  end function  SvdTk
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function CalcEccentricity(Day)
    !// --------------------------------------------------------------------
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    use cmn_parameters, only: CPI
    use cmn_ctm, only: LYEAR
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: Day ! day of year
    !// --------------------------------------------------------------------

    if (.not. LYEAR) then
       CalcEccentricity = 1._r8 &
         + 0.033_r8 * cos(2._r8 * CPI * real(Day - 10, r8) / 365._r8)
    else
       CalcEccentricity = 1._r8 &
            + 0.033_r8 * cos(2._r8 * CPI * real(Day - 10, r8) / 366._r8)
    end if
    !// --------------------------------------------------------------------
  end function CalcEccentricity
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function UnexposedLeafIRin(Tk, Eps)
    !// --------------------------------------------------------------------
    !// Calculate IR into leaf that is not exposed to the sky
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Eps, Tk
    ! Stefan-boltzman constant  W m-2 K-4
    real(r8), parameter :: Sb = 5.670367e-8_r8
    !// --------------------------------------------------------------------
    UnexposedLeafIRin = Eps * Sb * (Tk**4)
    !// --------------------------------------------------------------------
  end function UnexposedLeafIRin
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function ExposedLeafIRin(HumidPa, Tk)
    !// --------------------------------------------------------------------
    !// Calculate IR into leaf that is exposed to the sky
    !// Apparent atmospheric emissivity for clear skies: 
    !// function of water vapor pressure (Pa) 
    !// and ambient Temperature (K) based on Brutsaert(1975) 
    !// referenced in Leuning (1997)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: Tk, HumidPa
    real(r8) :: EmissAtm
    ! Stefan-boltzman constant  W m-2 K-4
    real(r8), parameter :: Sb = 5.670367e-8_r8
    !// --------------------------------------------------------------------

    EmissAtm        = 0.642_r8 * (HumidPa / Tk)**(1._r8/7._r8)
    ExposedLeafIRin = EmissAtm * Sb * (Tk**4)

    !// --------------------------------------------------------------------
  end function ExposedLeafIRin
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  !// SOIL NOx
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine SOILNOX( JDAY, TA, LSOIL, SOILCAT, SOILM, SOILT, &
                      LAIc, LAT, &
                      PRECADJ, &
                      CFNO, CFNOG )
    !// --------------------------------------------------------------------
    !  DESCRIPTION:
    !  
    !     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
    !     to estimate NO emissions 
    !     Information needed to estimate NO emissions
    !     Julian Day          (integer)    JDATE
    !     Surface Temperature (MCIP field) TA    (K)
    !     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
    !          (ratio of volume of water per volume of soil)
    !     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
    !     Soil Type           (MCIP field) ISLTYP            (LSOIL)
    !
    !     saturation values for soil types (constants)       (LSOIL)
    !     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
    !                and  the precipitation adjustment factor accounts for saturated soils
    !     FOR the non-PX version, the basic algorithm remains with a temperature
    !            adjustment factor (dry soil) and no adjustment for saturated soils
    !
    !
    !     The following arrays are updated after each call to SOILNOX
    !     PULTYPE   type of NO emission pulse 
    !     PULSEDATE julian date for the beginning of an NO pulse 
    !     PULSETIME        time for the beginning of an NO pulse
    !  
    !     The calculation are based on the following paper by J.J. Yienger and H. Levy II
    !     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
    !
    !     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
    !       the following modification (PX version):
    !       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is 
    !       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
    !       factor is made using the relative amount of soil moisture in the top layer (1cm)
    !       as the interpolating factor.  The relative amount of soil moisture is determined by
    !       taking the MCIP soil moisture field and dividing by the saturation value defined for each
    !       soil type in the PX version of MCIP
    !       the soil temperature is used in PX version
    !
    !     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
    !       The rainrate is computed from the MCIP directly using a 24 hr daily total. 
    !       THe types of Pulses as described in YL95 were used to estimate the NO emission
    !       rate.  
    !
    !    Also see the following paper for more information:
    !    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
    !    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
    !    by Tom Pierce and Lucille Bender       
    !
    !    REFERENCES
    !
    !    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
    !    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
    !    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
    !       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
    !        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
    !
    !  PRECONDITIONS REQUIRED:
    !     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
    !     NO emission pulse type, soil moisture from previous time step, julian date
    !     of NO emission pulse start, time of NO emission pulse start,
    !     soil type, SOIL TYPES, Land use data
    !
    !  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
    !     FERTILIZER_ADJ computes fertlizer adjustment factor
    !     VEG_ADJ        computes vegatation adjustment factor
    !     GROWSEASON     computes day of growing season
    !     
    !  REVISION  HISTORY:
    !    10/01 : Prototype by GAP
    !    10/03 : modified transition to non growing season for jul-oct of the year
    !    08/04 : Converted to SMOKE code style by C. Seppanen
    !    07/21/11 : Imported form SMOKE-BEIS v3.14 for MEGAN v2.10
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// In/out
    integer, intent(in)  ::  JDAY        !  Day of year
    real(r8), intent(in)  ::  TA          !  air temperature (K)
    real(r8), intent(in)  ::  SOILM       !  soil moisture (m3/m3)
    real(r8), intent(in)  ::  SOILT       !  soil temperature (K)
    real(r8), intent(in)  ::  PRECADJ     !  precip adjustment
    real(r8), intent(in)  ::  LAIc        !  soil temperature (K)
    real(r8), intent(in)  ::  LAT         !  Latitude
    integer, intent(in)  ::  SOILCAT   !  soil type
    logical, intent(in) :: LSOIL              ! true: using PX version of MCIP

    real(r8), intent(inout)  ::  CFNO    !  NO correction factor
    real(r8), intent(inout)  ::  CFNOG   !  NO correction factor for grass

    !// Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
    !// PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
    !// See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
    integer, parameter :: MAXSTYPES = 11
    real(r8), dimension(MAXSTYPES) :: SATURATION = (/ &
         0.395_r8, 0.410_r8, 0.435_r8, 0.485_r8, 0.451_r8, 0.420_r8, &
         0.477_r8, 0.476_r8, 0.426_r8, 0.482_r8, 0.482_r8 /)

    real(r8) :: &
         CF, &           ! NO correction factor
         CFG, &          ! NO correction factor for grasslands
         TAIR, &         ! surface temperature
         TSOI, &         ! soil temperature
         CFNOWET, CFNODRY, RATIO, VEG_ADJ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'SOILNOX'
    !// --------------------------------------------------------------------

    TAIR = TA         ! unit in degree K

    !// Check max and min bounds for temperature
    if (TAIR .lt. 200._r8) then
       write(6,'(a,f12.3)') f90file//':'//subr//': TAIR<200: ',TAIR
       stop
    end if

    if (TAIR .gt. 315._r8) then
       TAIR = 315._r8
    end if

    !// CFNOG
    if ( TAIR .gt. 303._r8 ) TAIR = 303._r8

    if ( TAIR .gt. 268.8690_r8 ) then
       CFG = exp( 0.04686_r8 * TAIR - 14.30579_r8 ) ! grass (from BEIS2)
    else
       CFG = 0._r8
    end if

    CFNOG = CFG

    !// CFNO
    if ( .not. LSOIL ) then
       ! no soil field
       TSOI = 0.72_r8 * TAIR + 82.28_r8
    else
       ! soil
       TSOI = SOILT
    end if

    IF (TSOI .le. 273.16_r8) TSOI = 273.16_r8
    IF (TSOI .ge. 303.16_r8) TSOI = 303.16_r8

    ! see YL 1995 Equa 9a p. 11452
    CFNODRY = (1._r8/3._r8) * (1._r8/30._r8) * (TSOI - 273.16_r8)
    if (TSOI .le. 283.16_r8) then
       ! linear cold case
       CFNOWET = (TSOI - 273.16_r8)*exp(-0.103_r8*30._r8)*0.28_r8 ! see YL 1995 Equ 7b
    else
       ! exponential case
       CFNOWET = exp(0.103_r8 * (TSOI - 273.16_r8)) * exp(-0.103_r8 * 30._r8)
    end if

    
    if ( .not. LSOIL ) then
       !// no soil
       CF = 0.5_r8 * CFNOWET + 0.5_r8 * CFNODRY
    else
       !// soil
       if ( SOILCAT .gt. 0 .and. SOILCAT .le. MAXSTYPES ) then
          RATIO = SOILM / SATURATION( SOILCAT )
          CF = RATIO * CFNOWET + max(0._r8, 1._r8 - RATIO) * CFNODRY
       else
          CF = 0._r8
       end if
    end if

    !// This internal function computes a vegetation adjustment factor
    !// based on LAIv.  See Yienger and Levy 1995
    VEG_ADJ = ( exp(-0.24_r8 * LAIc) + exp(-0.0525_r8 * LAIc) ) * 0.5_r8

    CFNO = CF * FERTLZ_ADJ( JDAY, LAT ) &
              * VEG_ADJ * PRECADJ

    !// --------------------------------------------------------------------
  end subroutine SOILNOX
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function FERTLZ_ADJ( JDAY, LAT )
    !// --------------------------------------------------------------------
    !  DESCRIPTION:
    !    This internal function computes a fertilizer adjustment factor
    !    for the given day of year. If it is not growing 
    !    season, the adjustment factor is 0; otherwise, it ranges from
    !    0.0 to 1.0.
    !
    !  CALL:
    !    GROWSEASON
    !  HISTORY:
    !    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
    !//
    !// Amund Sovde Haslerud, October 2017
    !//   Based on MEGANv2.10 code, modified for Oslo CTM3
    !// --------------------------------------------------------------------
    use utilities_oslo, only: GROWSEASON
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: JDAY
    real(r8), intent(IN) :: LAT

    !// Local variables
    integer ::  GDAY, GLEN
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'FERTLZ_ADJ'
    !// --------------------------------------------------------------------

    call GROWSEASON( JDAY, LAT, GDAY, GLEN )

    if ( GDAY .eq. 0 ) then
       FERTLZ_ADJ = 0._r8
    else if ( GDAY .ge. 1 .AND. GDAY .le. 30 ) then
       ! first month of growing season
       FERTLZ_ADJ = 1._r8
    else if ( GDAY .ge. 30 .AND. GDAY .le. 366) then
       ! later month of growing season
       FERTLZ_ADJ = 1._r8 + 30._r8 / real(GLEN,r8) - real(GDAY,r8) / real(GLEN,r8)
    else
       write( 6,'(a)') f90file//':'//subr// &
            ': Invalid date specified; date: growing season day = ', GDAY
       stop
    end if


    !// --------------------------------------------------------------------
  end function FERTLZ_ADJ
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function prec_adj(I,J,accRain24)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    use cmn_ctm, only: modelTimeIntegrated
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: I,J
    real(r8), intent(in) :: accRain24
    !// --------------------------------------------------------------------
    real(r8) :: currentRain24
    integer :: prevPulseType, curPulseType
    real(r8) :: pulseTime ! [s]
    !// --------------------------------------------------------------------

    !// The routine megan_update_metdata will initialise
    !// accumulated rain and also pulseType, pulseTime

    !// Here we only consider the next time steps

    prevPulseType = precPulseType(i,j)
    pulseTime = precPulseStartTime(i,j)

    !// Find current pulse type
    curPulseType = getPulseType(accRain24)

    if (curPulseType .le. prevPulseType) then
       !// No new pulse
       !// The function will update prevPulseType if necessary - i.e.
       !// when there is no longer a pulse to consider.
       prec_adj = precipFact(prevPulseType, pulseTime, modelTimeIntegrated)
       precPulseType(i,j) = prevPulseType
       !// If prevPulseType was set to 0, we could update time - but this has no consequence
!       if (prec_adj .gt. 1._r8) write(6,'(a9,2i4,4es12.4,2i2)') 'precadjO:',i,j,&
!            prec_adj,accRain24,pulseTime,modelTimeIntegrated,prevPulseType,curPulseType

    else 
       !// heavier precipitation - consider a new pulse
       prec_adj = PRECIPFACT(curPulseType,modelTimeIntegrated, modelTimeIntegrated)
       !// Update time and type
       precPulseType(i,j) = curPulseType
       precPulseStartTime(i,j) = modelTimeIntegrated
!       if (prec_adj .gt. 1._r8) write(6,'(a9,2i4,4es12.4,2i2)') 'precadjN:',i,j,&
!            prec_adj,accRain24,modelTimeIntegrated,precPulseStartTime(i,j),prevPulseType,curPulseType
    end if



    !// --------------------------------------------------------------------
  end function prec_adj
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  integer function getPulseType(RRATE)
    !// --------------------------------------------------------------------
    !// Computes the pulse type from a rainfall rate (See YL 1995).
    !// Unit of RRATE is cm/day.
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in) :: RRATE
    !// --------------------------------------------------------------------
           
    if (RRATE .lt. 0.1_r8) then
       getPulseType = 0
    else if (RRATE .lt. 0.5_r8) then
       getPulseType = 1
    else if (RRATE .lt. 1.5_r8) then
       getPulseType = 2
    else
       getPulseType = 3
    end if

    !// --------------------------------------------------------------------
  end function getPulseType
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  real(r8) function precipfact( PULSETYPE,  pulseTime, CurrentTime)
    !// --------------------------------------------------------------------
    !// This internal function computes a precipitation adjustment
    !// factor from YL 1995 based on a rain rate. The pulse type is
    !// an integer ranging from 0 to 3 indicating the type of 
    !// rainfall rate.
    !//
    !// HISTORY:
    !// Amund Sovde Haslerud, November 2017
    !//   Modified for CTM3.
    !// 07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
    !// --------------------------------------------------------------------
    use cmn_parameters, only: secDay
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(inout) :: PULSETYPE
    real(r8), intent(IN) :: pulseTime, currentTime
    !// Locals
    real(r8) :: daysdiff
    !// --------------------------------------------------------------------

    daysDIFF = (currentTime - pulseTime) / secDay

    select case( PULSETYPE )
    case( 0 )
       precipFact = 1._r8
    case( 1 )
       if( daysDIFF .lt. 2._r8 ) then
          precipFact = 11.19_r8 * EXP(-0.805_r8*(daysDIFF + 1._r8))
       else
          PULSETYPE = 0
          precipFact = 1._r8
       end if
    case( 2 )
       IF( daysDIFF .lt. 6._r8 ) then
          precipFact = 14.68_r8 * EXP(-0.384_r8*(daysDIFF + 1._r8))
       else
          PULSETYPE = 0
          precipFact = 1._r8
       end if
    case DEFAULT
       if ( daysDIFF  < 13._r8 ) then
          precipFact = 18.46_r8 * EXP(-0.208_r8*(daysDIFF + 1._r8))
       else
          PULSETYPE = 0
          precipFact = 1._r8
       end if
    end select

    !// --------------------------------------------------------------------
  end function precipfact
  !// ----------------------------------------------------------------------











  !// ----------------------------------------------------------------------
end module emissions_megan
!//=========================================================================

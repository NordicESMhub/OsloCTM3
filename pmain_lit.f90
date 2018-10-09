!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Lightning generation.
!// Stripped version of Oslo CTM3 pmain.f90
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// p-main.f90 version of UCIrvine CTM p-main.f.
!//=========================================================================
!// Lightning generation: Stripped of transport and chemistry.
!//=========================================================================
program pmain
  !//-----------------------------------------------------------------------
  !// Modules to use
  !//-----------------------------------------------------------------------
  use cmn_precision, only: r8, r4, rTnd, rMom
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, LLINOZ, &
       LOSLOCHEM, LBCOC, LSOA, LSALT, LDUST
  use cmn_ctm, only: AIR, LYEAR, GMTAU, NRMETD, NROPSM, NRCHEM, JDATE, &
       JDAY, JMON, JYEAR, NTM, IDAY, IYEAR, LLPYR, LFIXMET, TMON, TMET, &
       SOLDEC, SOLDIS, LCONT, STT
  use cmn_chem, only: N_LZ, N_STE, LZLBO3, O3iso1,O3iso2
  use cmn_diag, only: NTDPAR, USTEP, VSTEP, WSTEP, JMON0, LFLXDG, LMXDG, &
       STTTND, NDAY0, JDO_A, JDO_C, JDO_X, JDO_T, &
       NTND_SOURCE, NTND_BNDLYR, NTND_DRYDEP, NTND_CHEM, &
       NTND_LSSCAV, NTND_CNSCAV, NTND_WADV, NTND_UVADV
  use cmn_fjx, only: JVN_
  use cmn_met, only: MYEAR,precls, CLDFR, CLDIWC, CLDLWC
  use cmn_parameters, only: LDEBUG
  !//-----------------------------------------------------------------------
  use averages, only: AVG_P1, AVG_WRT2, AVG_ADD2, AVG_CLR2
  use budgets, only: TBGT_G, TBGT_P0, TBGT_P1, TBGT_P2, TBGT_IJ, TBGT_L
  use cloudjx, only: cloud_init
  use convection, only: CONVW_OSLO
  use fastjx, only: SET_ATM
  use grid, only: AIRSET
  use initialize, only: input, setup_species, setup_unf_output
  use lightning, only: lightning_oas2015
  use metdata_ecmwf, only: update_metdata
  use omp, only: MPSPLIT, MPBIND
  use pbl_mixing, only: CNVDBL
  use scavenging_drydep_uci, only: DRYDEP
  use scavenging_largescale_uci, only: WASHO
  use source_uci, only: source
  use steflux, only: CHEMFLUX, CHEMFLUX_E90, dumpuvflux, DUMPTRMASS, &
       DUMPTRMASS_E90, SAVETRMASS, STEBGT_CLR, STEBGT_WRITE, &
       ctm3_pml, ctm3_o3scav
  use stt_save_load, only: oslo_con_sav
  use utilities, only: write_log, ctmExitC, LCM, CALENDR, CALENDL, &
       get_dinm, check_btt
  !//-----------------------------------------------------------------------
  use cmn_oslo, only: METHANEMIS, JVAL_IJ
  use emissions_aircraft, only: aircraft_h2o_zero
  use bcoc_oslo, only: bcsnow_master, bcsnow_save_restart
  use ch4routines, only: set_ch4_stt, setch4sfc
  use diagnostics_general, only: &
       init_daily_diag, daily_diag_output, nops_diag, &
       mp_diag, TBGT_2FILE, reports_chemistry, &
       tnd_emis2file, chembud_output
  use diagnostics_scavenging, only: &
       scav_diag_ls, scav_diag_cn, scav_diag_brd, scav_diag_2fileA
  use drydeposition_oslo, only: setdrydep
  use dust_oslo, only: dustbdg2d
  use emissions_ocean, only: emissions_ocean_total
  use emissions_oslo, only: update_emis, update_emis_ij
  use gmdump3hrs, only: dump3hrs
  use input_oslo, only: init_oslo
  use main_oslo, only: master_oslo, update_chemistry
  use physics_oslo, only: update_physics
  use seasalt, only: seasaltbdg2d, emissions_seasalt_total
  use utilities_oslo, only: get_chmcycles, &
       source_e90, decay_e90, &
       tpauseb_e90, tpauseb_o3
  !//-----------------------------------------------------------------------
  use stratchem_oslo, only: update_strat_boundaries
  use strat_h2o, only: set_h2_eurohydros
  use soa_oslo, only: soa_diag2file, soa_nopsdiag
  use fallingaerosols, only: aerosolsettling
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------

  !//---day number, starting and ending day
  integer :: NDAY, NDAYI, NDAYE
  !//---logical switches
  logical :: LNEWM, LEND
  logical :: LDO_ADV, LDO_CHM, LTBGT, LAVGS, LSERI, L2

  !//---time steps -TAU- (hr) and DT-- (s)
  real(r8) :: DTAU, DTAUMT, DTAUCH, DTAULCM, UTTAU, DAY
  real(r8) :: DTOPS, DTCHM, DTMET, DTADV, DTLCM
  !//---CTM3: time steps variables
  integer :: CHMCYCLES, CCYC
  real(r8) :: dtchm2

  !//---private arrays B--(reverse order) for IJ-block OMP calculation
  real(r8), dimension(LPAR,IDBLK,JDBLK)      :: AIRB, GAMAB, GAMACB, BTEM
  real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT, BTTBCK
  real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
       BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
  real(r8), dimension(MPBLK,NPAR,NTDPAR)   :: BTTTN0
  real(rTnd), dimension(LPAR,NPAR,IDBLK,JDBLK,NTDPAR)::  BTTTND
  logical, dimension(LPAR,IDBLK,JDBLK)     :: LPAUZB
  logical, dimension(LPAR,IPAR,JPAR,2)     :: LO3PAUZ

  !//---private array AIRUV for UV-advection
  real(r8), dimension(IPAR,JPAR)   ::  AIRUV
  real(r8), dimension(IPAR+1,JPAR) ::  ALFA2D
  real(r8), dimension(IPAR,JPAR+1) ::  BETA2D
  real(r8), dimension(LPAR,NPAR,NTDPAR) ::  STTTNL

  !//---private arrays for STE flux calculation
  real(r8), dimension(IPAR+1,JPAR,NPAR,2) :: QFU
  real(r8), dimension(IPAR,JPAR+1,NPAR,2) :: QFV
  real(r8), dimension(IPAR+1,JPAR+1,2)    :: CHMFLX
  real(r8), dimension(IPAR+1,JPAR,LPAR)   :: UFLX1, UFLX2, UFLX3, UFLX4
  real(r8), dimension(IPAR,JPAR+1,LPAR)   :: VFLX1, VFLX2, VFLX3, VFLX4

  !//---loop counter limits
  integer :: NMET,NOPS,NADV,NLCM,NSUB, KSERIES
  integer :: I, J, M, LUV,N, L, MSTEPU, MSTEPV,MSTEPW

  !//---random number for clouds
  integer :: NCLDRAN

  !//---Boolean table for Linoz
  !real(r4),  dimension(IPAR,JPAR,LPAR) ::  BLNAVG
  !integer, dimension(IPAR,JPAR,LPAR) ::  NRBLN

  !//---for timing the run
  character(len=10) :: BIG_BEN(3)
  integer :: DATE_TIME(8), START_TIME(8)
  real(r8) :: tot_dt, nops_dt, nday_dt
  !//-----------------------------------------------------------------------

  !// Start time
  call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
  start_time(:) = date_time(:)
  tot_dt = - (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
       + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)

  !// Write start info
  call write_log(0, start_time, start_time)

  !// Input (std in) and initialization, also calls calendar
  call INPUT(NDAYI, NDAYE)


!//UPDATE Will move to general init routine?
  !// Initialize clouds (Set random number RAN4, after INPUT)
  call cloud_init()


  !// Swiches for e.g. budgets
  LTBGT = .false.
  LSERI = .false.
  LAVGS = .false.
  L2    = .false.
  LNEWM = .true.

  !// Begin CTM calculation at 0000 UT on NDAYI relative to IYEAR
  !//-----------------------------------------------------------------------
  NDAY   = NDAYI
  IDAY   = NDAYI
  DAY    = real(NDAYI, r8)
  GMTAU  = 0._r8
  JMON0  = 0
  KSERIES = 0

  !// Time steps
  DTAUMT = 24._r8 / real(NRMETD, r8)  ! in hours
  DTAU   = DTAUMT / real(NROPSM, r8) ! in hours
  DTAUCH = DTAU / real(NRCHEM, r8)   ! in hours
  DTOPS  = 3600._r8 *DTAU    ! in seconds
  DTCHM  = 3600._r8 *DTAUCH  ! in seconds
  DTMET  = 3600._r8 *DTAUMT  ! in seconds

  !// Internal Oslo chemistry loop for CHEM/BLMIX (# cycles per NOPS)
  call get_chmcycles(CHMCYCLES, NRCHEM, NROPSM)

  !// Initialize metdata and air
  !//-----------------------------------------------------------------------
  call update_metdata(0, DTMET, 1)
  call AIRSET(DTMET)


  !// TRACER SPECIFIC INPUT
  !//-----------------------------------------------------------------------

  !// Separate read-in for Oslo chemistry
  call init_oslo(NDAY,DTCHM,CHMCYCLES)

  !// Set up species (initialize.f90)
  !call SETUP_SPECIES(NDAY,NDAYI)
  !call SETUP_UNF_OUTPUT(NDAY,NDAYI)

  !// RESET SELECTED TRACERS
  !// Set CH4 to HYMN-results
  !call set_ch4_stt()
  !// Set (non-transported) H2 to Eurohydros-results
  !call set_h2_eurohydros()
  !// Set aircraft H2O to zero
  !call aircraft_h2o_zero()

  !// Clear averages (AVG_CLR2 in oc_diagnoses.f)
  !call AVG_CLR2
  !// Initialize budgets
  !call TBGT_G(STTTNL,BTTTN0,0)

  !// STEFLUX: Initialize
  !if (LFLXDG) call STEBGT_CLR(NDAY,0) !// Initialize STE budgets

  !//-----------------------------------------------------------------------
  !// MAIN LOOP
  !//-----------------------------------------------------------------------
  !// NDAY(NDAYI:NDAYE-1)
  !//    Day number, fundamental timing unit, spans years (can be >365)
  !//    Note that run is from NDAYI.00 to NDAYE.00 = NDAYE-NDAYI days
  !//    NB: Leap Years are NOT correct although met field are avvailable
  !//  NMET(1:NRMETD)
  !//     Number of met field for the DAY (eg, 1:8 for 3-hr fields)
  !//     For UCI CTM we have required 8 to resolve dirnal cycle
  !//   NOPS(1:NROPSM)
  !//      Number of operator-split steps per met field (eg 3 1-hr steps)
  !//      This is arbitrary, but controls the finest time for
  !//      diagnostics.
  !//    NSUB(1:NLCM)
  !//       Sub-steps of op-split time to meet both CHEM or ADV criteria.
  !//       This can be asynchronous:  1/2-hr for ADV, 1/3-hr for CHEM.
  !//
  !// DTAU = 24 hr / (NRMETD * NROPSM) = master time step
  !//       This controls the frequency at which diagnostics can be made
  !//       shorter, computational time scales can be forced globally by 
  !//         (1) the CFL advection criteria (divegence from any grid box),
  !//             this is calculated dynamically for each met field.
  !//      or (2) the CHEMISTRY needs (at present set at runtime)
  !//      These shorter intervals (1/NADV, 1/NRCHEM) are run asynnchronously
  !//
  !// GMTAU = UT (Greenwich) time in hours (0.00 to 24.00)
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  do NDAY = NDAYI, NDAYE-1
    !//---------------------------------------------------------------------

    !// Date and time stuff
    call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
    nday_dt = - (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
         + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)
    if (NDAY.eq.NDAYI) tot_dt = tot_dt - nday_dt !// Time spent on init

    IDAY = NDAY !// Set IDAY

    !// Initialize LINOZ each day
    !if (LLINOZ)  call LNZ_SET(JMON,JYEAR)

    !// Daily initialization of diagnostics (oc_diagnostics.f)
    !call init_daily_diag(JYEAR,JMON,JDATE,NDAY,NDAYI)

    !// Set up number of days in months
    call get_dinm(LYEAR)


    !//---------------------------------------------------------------------
    do NMET = 1, NRMETD
      !//-------------------------------------------------------------------

      call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      nops_dt = - (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
           + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)

      !// Read met fields & set up surface
      call update_metdata(1, DTMET, NMET)

      if (LOSLOCHEM) then
        !// Set climatological mass, O3 & T above CTM
        call SET_ATM
        !// Find lightning
        call LIGHTNING_OAS2015(NDAY,NDAYI,DTMET,LNEWM)
        NCLDRAN   = 0
      end if
      call DYN0 (DTMET)             ! process wind/advection fields
      call CFLADV (DTOPS,NADV)      ! global cfl criteria
      call LCM (NRCHEM,NADV,NLCM)   ! find least common multiple

      !//---DTxxx = time step in seconds 
      DTLCM   = DTOPS / real(NLCM, r8)
      DTADV   = DTOPS / real(NADV, r8)
      DTCHM   = DTOPS / real(NRCHEM, r8) ! same as defined above
      DTAULCM = DTLCM / 3600._r8
      DTCHM2 = DTCHM / real(CHMCYCLES, r8) ! for Oslo chemistry internal loop

      !// Dump aerosol data each 3 hour? (set LDUMP3HRS in oc_gmdump3hrs.f90)
      !call dump3hrs(NDAYI,NDAY,NMET)

      !//-------------------------------------------------------------------
      do NOPS = 1, NROPSM
        !//-----------------------------------------------------------------
        !//---This loop is the core of the CTM, the choice of order is
        !//---critical here, but can be tested and should converge as
        !//---NROPSM ==> infinity.
        !//
        !// Parallelized over IJ-blocks
        !//---SOURCE:  2-D(surface) and 3-D emission sources
        !//---CNVDBL:  Boundary layer mixing
        !//---DRYDEP:  Dry deposition at the surface
        !//---CHEM:    Full chemistry package
        !//---WASHO:   Washout of tracers by large-scale precip
        !//---CONVW:   Convective transport (vertical, NOT subsidence)
        !//            and plume scavenging
        !//---DYN2W:   Vertical advection = large-scale (U&V) plus
        !//            convective subsidence
        !//
        !// Parallelized over N-tracer or L-levels
        !//---DYN2V:   N-S advection
        !//---DYN2U:   E-W advection
        !//
        !//---Note that time steps may be different:
        !//   1/NRCHEM = DTAU fraction for: SOURCE+CNVDBL+DRYDEP+CHEM
        !//   1/NADV   = DTAU fraction for: WASHO+CONVW(et)+DYN2W+DYN2V+DYN2U
        !//
        !//---Diagnoses (names set in LxxTyy-file)
        !//--- 0=INITIA, 1=SOURCE, 2=BndryL, 3=DRYDEP, 4=UV_ADV, 5=W_ADV_,
        !//--- 6=LSSCAV, 7=CHEMIS, 8=C_SCAV, 9=
        !call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
        !nops_dt = - (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
        !     + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)

        !// Update global variables needed for Oslo chemistry.
        !// Also initialises E90-tracer and tropopause.
        !call update_physics(NDAYI,NDAY,NMET,NOPS,JMON,LNEWM)

        !// Update parameters for chemistry and aerosols.
        !// Also initialises stratospheric H2O.
        !call update_chemistry(NDAYI,NDAY,NMET,NOPS,JMON,LNEWM)

        !// Update short term variations in emissions.
        !// Note that some parts may be updated within MP-block, using
        !// subroutine update_emis_ij.
        !call update_emis(NDAY,NMET,NOPS,NDAYI,LNEWM,DTOPS)


        !// Do NOPS-based diagnosics.
        !// Remember that for asynchronous chemistry/transport
        !// steps, the previous process MAY NOT BE transport.
        !call nops_diag(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI,LNEWM)
        


        !//-----------------------------------------------------------------
        do NSUB = 1, NLCM
          !//---------------------------------------------------------------
          !//---This loop provides an asynchronous sub-stepping between the
          !//---minimum times set by the global advection limitation
          !//---(dynamically calculated for each met field) and the
          !//---chemistry group which can demand a sub-split of NOPS.
          !//---GMTAU = UTTAU = Universal Time (in hr) at beginning of
          !//---this op-split step
          UTTAU  = GMTAU
          LDO_ADV = mod((NSUB-1)*NADV,NLCM).eq.0   !// Do advection?
          LDO_CHM = mod((NSUB-1)*NRCHEM,NLCM).eq.0 !// Do chemistry?
          if (LDO_CHM) NCLDRAN = NCLDRAN+1


          !// Accumulate global budgets
          !call TBGT_G (STTTNL,BTTTN0, 1)

          GMTAU = GMTAU + DTAULCM

          !//---------------------------------------------------------------
        end do !// NSUB                                                 NSUB
        !//-----------------------------------------------------------------

        !---have completed the basic op-split time step = DTAU
        DAY = real(NDAY, r8) + GMTAU / 24._r8

        !call DUMP3D (25, JYEAR,DAY)

        !// Add to averages (in averages.f90)
        !call AVG_ADD2

        !call TSER_1 (KSERIES)

        !// Report chemistry & physics
        !call REPORTS_CHEMISTRY(NDAY,NMET,NOPS,GMTAU, NDAYI,LNEWM)

        !// NOPS diag for SOA
        !if (LSOA) call soa_nopsdiag(NDAY,NMET,NOPS,JDAY,DTOPS)

        !// Done everything for the beginning of new month
        if (LNEWM) LNEWM = .false.

      !// Additional printout
      !call emissions_seasalt_total()
      !call emissions_ocean_total()

        !call DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
        !nops_dt = nops_dt + (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
        !     + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)
        !write(6,'(a,i4,2i3,f8.2)') 'WC NOPS [s]: ',NDAY,NMET,NOPS,nops_dt
        !//-----------------------------------------------------------------
      end do !// NOPS                                                   NOPS
      !//-------------------------------------------------------------------
      !//---have completed the time step of the met fields



      !if (LFLXDG .and. N_STE.ne.0)  then
      !  call DUMPTRMASS(LO3PAUZ,1)
      !  call DUMPTRMASS(LO3PAUZ,2)
      !  call DUMPTRMASS_E90(3,N_STE)
      !  if (LLINOZ) call DUMPTRMASS_E90(4,N_LZ)
      !end if

      call DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
      nops_dt = nops_dt + (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
           + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)
      write(6,'(a,i4,i3,f8.2)') 'WC NMET [s]: ',NDAY,NMET,nops_dt

      !//-------------------------------------------------------------------
    end do !// NMET                                                     NMET
    !//---------------------------------------------------------------------

    !// STEFLUX: Save daily tropospheric O3 mass
    !if (LFLXDG .and. N_STE .gt. 0) then
    !  !// Create the very initial trop column mass
    !  if (NDAY .eq. NDAYI)  call SAVETRMASS(0)
    !  call SAVETRMASS(1)
    !end if

    !//---have completed 24-hour (day) step
    GMTAU   = 0._r8
    LEND = NDAY.eq.NDAYE-1

    !// Daily/monthly diagnostics & saves
    !//---------------------------------------------------------------------

    !// Daily output of diagnostics; before calendar update (oc_diagnostics.f)
    !call daily_diag_output(JYEAR,JMON,JDATE,JDAY,NDAY,NDAYI)


    !//---update the calendar to the new, upcoming day
    IDAY = NDAY + 1
    call CALENDR (IYEAR,IDAY, LLPYR,LFIXMET,MYEAR, &
         JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
    LNEWM = JDATE.eq.1

    !// STEFLUX: flux dump share calendar JDO_X
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_X, LMXDG,L2)
    !if (LMXDG .and. LFLXDG)  then
    !  call STEBGT_WRITE
    !  call STEBGT_CLR(IDAY,1)
    !end if

    !//---TENDENCIES:  write/reset budget tendencies & time series
    !//---use calendar to trigger diagnostics/saves on/off
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_T, LTBGT,L2)

    !// Initialize 2D budgets files
    !if (NDAY .eq. NDAYI) call TBGT_2FILE(0)

    !if (LTBGT .or. LEND) then
      !// Write diagnostics: emissions
      !call tnd_emis2file(JYEAR,JMON,JDATE,NDAY)

      !// Print out tendencies 0D, 1D and 2D
      !call TBGT_P0 (IDAY)
      !call TBGT_P1 (IDAY)
      !call TBGT_P2 (IDAY) !// Unformatted output to 24 currently disabled

      !// Print out SALT budget (NDAY, not IDAY)
      !if (LSALT) call seasaltbdg2d(NDAY,NDAYI, STTTND, NDAY0)

      !// Print out DUST budget
      !if (LDUST) call dustbdg2d(NDAY,NDAYI, STTTND, NDAY0)

      !// Print 2D budgets to file
      !call TBGT_2FILE(1)

      !if (LSOA) call soa_diag2file(NDAY,NDAYI)

      !// Chemistry budgets
      !call chembud_output(JYEAR,JMON,JDATE,NDAY)

      !// Write diagnostics: scavenging
      !call scav_diag_2fileA(JYEAR,JMON,JDAY,NDAY)

      !// Initialize core tendencies for next average
      !call TBGT_G(STTTNL, BTTTN0, 0)

    !end if

    !//---SERIES:  write TRACER time series (unf unit=23)
    !//---use calendar to trigger diagnostics/saves on/off
    !//---reset (KSERIES=0) must be done every day since series is always accum.
    !call CALENDL (JYEAR,JDAY,LYEAR,  JDO_S, LSERI,L2)
    !if (LSERI) then
    !  call TSER_2 (23,KSERIES)
    !end if
    KSERIES = 0


    !//---write/reset average tracer
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_A, LAVGS,L2)
!    if (LAVGS .or. LEND) then
!      call AVG_WRT2    !// Write averages
!
!      if (L2) then
!        call AVG_P1
!      end if
!
!      call AVG_CLR2    !// Clear averages
!    end if

    !//---restart file save (at least at end of run)
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_C, LCONT,L2)
!    if (LCONT .or. LEND)  then
!      call OSLO_CON_SAV(NDAY+1)
!      if (LBCOC) call bcsnow_save_restart(NDAY+1)
!    end if

    call DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),DATE_TIME)
    nday_dt = nday_dt + (date_time(3)*86400._r8 + date_time(5)*3600._r8 &
         + date_time(6)*60._r8 + date_time(7) + date_time(8)*1.e-3_r8)
    tot_dt = tot_dt + nday_dt !// Summing up total time
    write(6,'(a,i4,f8.1,1x,f10.2,1x,f7.2)') &
         '** WC time NDAY [s] / total [hours/days]:', &
         NDAY,nday_dt, tot_dt/3600._r8, tot_dt/86400._r8
    !//---------------------------------------------------------------------
  end do !// NDAY                                                       NDAY
  !//-----------------------------------------------------------------------

  !//---Final: closeout the model run
  close(21)
  call write_log(1, start_time, date_time)
  call ctmExitC(' std exit')
  !//-----------------------------------------------------------------------
end program pmain
!//=========================================================================

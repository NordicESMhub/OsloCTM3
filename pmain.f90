!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Amund Sovde Haslerud, May 2017
!//=========================================================================
!// pmain.f90 version of UCIrvine CTM p-main.f.
!//=========================================================================
!// Driver for CHEMICAL TRANSPORT MODEL
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
       JDAY_NEXT, JDATE_NEXT, JMON_NEXT, JYEAR_NEXT, &
       SOLDEC, SOLDIS, LCONT, modelTimeIntegrated, &
       LSTOM1HRS, LDLYSCAV
  use cmn_chem, only: N_LZ, N_STE, LZLBO3, O3iso1,O3iso2
  use cmn_diag, only: NTDPAR, USTEP, VSTEP, WSTEP, JMON0, LFLXDG, LMXDG, &
       STTTND, NDAY0, JDO_A, JDO_C, JDO_X, JDO_T, &
       NTND_SOURCE, NTND_BNDLYR, NTND_DRYDEP, NTND_CHEM, &
       NTND_LSSCAV, NTND_CNSCAV, NTND_WADV, NTND_UVADV
  use cmn_fjx, only: JVN_
  use cmn_met, only: MYEAR
  use cmn_parameters, only: LDEBUG
  use cmn_sfc, only: LDDEPmOSaic
  !//-----------------------------------------------------------------------
  use averages, only: AVG_P1, AVG_WRT_NC4, AVG_ADD2, AVG_CLR2
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
  use stt_save_load, only: save_restart_file
  use utilities, only: write_log, LCM, calendar, get_soldecdis, &
       CALENDL, get_dinm, check_btt
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
       scav_diag_ls, scav_diag_cn, scav_diag_brd, scav_diag_2fileA, &
       scav_diag_nmet_output_nc
  use drydeposition_oslo, only: setdrydep
  use dust_oslo, only: dustbdg2file, dustInstBdg
  use emissions_ocean, only: emissions_ocean_total
  use emissions_oslo, only: update_emis, update_emis_ij
  use gmdump3hrs, only: dump3hrs
  use input_oslo, only: init_oslo
  use main_oslo, only: master_oslo, update_chemistry
  use physics_oslo, only: update_physics, set_blh_ij
  use seasalt, only: saltbdg2file, emissions_seasalt_total
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
  real(r8) :: dtchm2, nmetTimeIntegrated(MPBLK)

  !//---private arrays B--(reverse order) for IJ-block OMP calculation
  real(r8), dimension(LPAR,IDBLK,JDBLK)      :: AIRB, GAMAB, GAMACB, BTEM
  real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT, BTTBCK
  real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
       BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
  real(r8), dimension(MPBLK,NPAR,NTDPAR)   :: BTTTN0
  real(rTnd), dimension(LPAR,NPAR,IDBLK,JDBLK,NTDPAR)::  BTTTND
  logical, dimension(LPAR,IDBLK,JDBLK)     :: LSTRATAIR_E90B
  logical, dimension(LPAR,IPAR,JPAR,2)     :: LSTRATAIR_O3

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
  integer :: END_TIME(8), START_TIME(8)
  integer(kind=8) :: stime
  real(r8) :: tot_dt, nops_dt, nday_dt, nmp_dt, nuv_dt, rtime
  !//-----------------------------------------------------------------------

  !// Start time and system clock for timings
  call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),START_TIME)
  call SYSTEM_CLOCK(count=stime, count_rate=rtime)
  tot_dt = - real(stime, r8) / rtime

  !// Write start info
  call write_log(0, start_time, start_time, 0._r8)

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

  !// Separate read-in for Oslo chemistry (must be after AIRSET)
  call init_oslo(NDAY,DTCHM,CHMCYCLES)

  !// Set up species (initialize.f90)
  call SETUP_SPECIES(NDAY,NDAYI)
  call SETUP_UNF_OUTPUT(NDAY,NDAYI)

  !// RESET SELECTED TRACERS
  !// Set CH4 to HYMN-results
  !call set_ch4_stt()
  !// Set (non-transported) H2 to Eurohydros-results
  !call set_h2_eurohydros()
  !// Set aircraft H2O to zero
  !call aircraft_h2o_zero()

  !// Clear averages
  call AVG_CLR2
  !// Initialize budgets
  call TBGT_G(STTTNL,BTTTN0,0)

  !// STEFLUX: Initialize
  if (LFLXDG) call STEBGT_CLR(NDAY,0) !// Initialize STE budgets

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

    !// Keep track of seconds calculated per day
    call SYSTEM_CLOCK(count=stime)
    nday_dt = - real(stime, r8) / rtime
    if (NDAY.eq.NDAYI) tot_dt = tot_dt - nday_dt !// Time spent on init

    IDAY = NDAY !// Set IDAY

    !// Initialize LINOZ each day
    if (LLINOZ)  call LNZ_SET(JMON,JYEAR)

    !// Daily initialization of diagnostics (diagnostics_general.f90)
    call init_daily_diag(JYEAR,JMON,JDATE,NDAY,NDAYI)

    !// Set up number of days in months
    call get_dinm(LYEAR)


    !//---------------------------------------------------------------------
    do NMET = 1, NRMETD
      !//-------------------------------------------------------------------

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
      write(6,'(a,2i3)') 'NADV/NLCM: ',NADV,NLCM

      !//---DTxxx = time step in seconds 
      DTLCM   = DTOPS / real(NLCM, r8)
      DTADV   = DTOPS / real(NADV, r8)
      DTCHM   = DTOPS / real(NRCHEM, r8) ! same as defined above
      DTAULCM = DTLCM / 3600._r8
      DTCHM2 = DTCHM / real(CHMCYCLES, r8) ! for Oslo chemistry internal loop
      nmetTimeIntegrated(:) = 0._r8   ! sums up seconds calculated in CCYC-loop

      !// Dump aerosol data each 3 hour? (gmdump3hrs.f90)
      call dump3hrs(NDAYI,NDAY,NMET)

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
        call SYSTEM_CLOCK(count=stime)
        nops_dt = - real(stime, r8) / rtime

        !// Update global variables needed for Oslo chemistry.
        !// Also initialises E90-tracer and tropopause.
        call update_physics(NDAYI,NDAY,NMET,NOPS,LNEWM)

        !// Update parameters for chemistry and aerosols.
        !// Also initialises stratospheric H2O.
        call update_chemistry(NDAYI,NDAY,NMET,NOPS,JMON,LNEWM)

        !// Update short term variations in emissions.
        !// Note that some parts may be updated within MP-block, using
        !// subroutine update_emis_ij.
        call update_emis(NDAY,NMET,NOPS,NDAYI,LNEWM,DTOPS)


        !// Do NOPS-based diagnosics.
        !// Remember that for asynchronous chemistry/transport
        !// steps, the previous process MAY NOT BE transport.
        call nops_diag(JYEAR,JMON,JDATE,NDAY,NMET,NOPS,NDAYI,LNEWM)
        
        nmp_dt = 0._r8
        nuv_dt = 0._r8

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

!//-------------------------------------------------------------------------
          if (LDO_ADV .or. LDO_CHM) then
            call SYSTEM_CLOCK(count=stime)
            nmp_dt = nmp_dt - real(stime, r8) / rtime
!//-------------------------------------------------------------------------
!$omp parallel private(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,AIRB, &
!$omp                  BTEM,GAMAB,GAMACB,BTTBCK,BTTTND,LSTRATAIR_E90B, &
!$omp                  MSTEPW,M,CCYC) &
!$omp          shared (UTTAU,DTCHM,DTADV,DTCHM2,NTM,NCLDRAN,LSTRATAIR_O3, &
!$omp                  NDAY,NMET,NOPS,NSUB,CHMCYCLES,LDO_CHM,LDO_ADV, &
!$omp                  WSTEP,N_LZ,N_STE,LZLBO3,BTTTN0, &
!$omp                  NTND_SOURCE, NTND_BNDLYR, NTND_DRYDEP, NTND_CHEM, &
!$omp                  NTND_LSSCAV, NTND_CNSCAV, NTND_WADV, &
!$omp                  nmetTimeIntegrated, DTMET) &
!$omp          default(NONE)
!$omp do schedule(dynamic)
!//-------------------------------------------------------------------------
            do M = 1, MPBLK !// begin IJ-block

              !// Split into IJ-blocks / MP-blocks
              call MPSPLIT(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                            BTTTND,BTTBCK,AIRB,BTEM,GAMAB,LSTRATAIR_E90B,M)
 
              if (LDEBUG) call check_btt(BTT,M,'after mp-splitting')

              !// Diagnose total burden of species (for wet scav. diagnose)
              call scav_diag_brd(NSUB, NOPS, NMET, BTT, M)

!//-------------------------------------------------------------------------
              if (LDO_CHM) then
!//-------------------------------------------------------------------------

                !// Internal UiO loop to treat mixing and chemistry
                do CCYC = 1, CHMCYCLES

                  !// Calculate BLH for this time step
                  !call set_blh_ij(NMET, NOPS, NSUB, CCYC, DTCHM2, &
                  !     dtmet, nmetTimeIntegrated(M), M)

                  !// Update short term emissions in MP-block (emissions_oslo.f90)
                  call update_emis_ij(NMET,NOPS,NSUB,CCYC,M)

                  !// CTM3 diagnoses inside IJ-block (diagnostics_general.f90)
                  call mp_diag(BTT,NDAY,NMET,NOPS,NSUB,CCYC,M)
                  !// Reset CH4 at surface
                  if (.not.METHANEMIS) call setch4sfc(BTT,BTTBCK,AIRB,M)

                  !// Emissions
                  call SOURCE(BTT,BXT,BXX,BYT,BYY,BXY,BZT,DTCHM2,M)
                  call source_e90(BTT,DTCHM2,M)
                  call TBGT_IJ(BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_SOURCE)

                  !// Boundary layer mixing
                  call CNVDBL(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,&
                              AIRB,DTCHM2,M)
                  if (.not.METHANEMIS) call setch4sfc(BTT,BTTBCK,AIRB,M)
                  !// add to PML
                  if (N_STE .gt. 0) &
                       call ctm3_pml(BTT,BTTBCK,AIRB,LSTRATAIR_E90B,LSTRATAIR_O3,M)
                  call TBGT_IJ(BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_BNDLYR)


                  !// Update tropopause
                  call tpauseb_e90(BTT,AIRB,LSTRATAIR_E90B,M)
                  call tpauseb_o3(BTT,AIRB,LSTRATAIR_O3,M)

                  !// Linoz: calculate PML
                  if (LLINOZ) then
                    !// Also calculates PML for STE. When calculating
                    !// STE using CTM3 O3, this is done as tendency due to
                    !// chemistry (see master_oslo in main_oslo.f90)
                    call LNZ_PML(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                                 AIRB,BTEM,LSTRATAIR_E90B,UTTAU,DTCHM2,M,LZLBO3, N_LZ)
                  end if

                  !// Set dry deposition velocities (drydeposition_oslo.f90)
                  if (NSUB .eq. 1 .and. CCYC .eq. 1) &
                       call setdrydep(UTTAU, BTT, AIRB, BTEM, M)

                  !// Chemistry
                  call master_oslo(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                       BTEM,AIRB,BTTBCK,UTTAU,DTCHM,DTCHM2, &
                       NDAY,NMET,NOPS,NSUB,NCLDRAN,CCYC,M)

                  !// CTM3 O3: Find PML from surface up to tropopause
                  if (N_STE .gt. 0) &
                       call ctm3_pml(BTT,BTTBCK,AIRB,LSTRATAIR_E90B,LSTRATAIR_O3,M)

                  !// Decay the e90-tracer
                  call decay_e90(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                                 DTCHM2,M)
                  call TBGT_IJ(BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_CHEM)


                  !// Dry deposition
                  call DRYDEP(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                              UTTAU,DTCHM2,M)
                  call TBGT_IJ(BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_DRYDEP)

                  nmetTimeIntegrated(M) = nmetTimeIntegrated(M) + DTCHM2


                end do !// do CCYC = 1, CHMCYCLES

!//-------------------------------------------------------------------------
              end if !// end LDO_CHM
!//-------------------------------------------------------------------------


!//-------------------------------------------------------------------------
              if (LDO_ADV) then
!//-------------------------------------------------------------------------

                if (.not.METHANEMIS) call setch4sfc(BTT,BTTBCK,AIRB,M)

                !// Update boundary cond. for stratosphere before transport
                call update_strat_boundaries(BTT,BTTBCK,AIRB,BTEM,DTADV,M)

                !// Large scale washout
                call WASHO(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                           AIRB,BTEM,DTADV,M)
                !// Diagnose daily scavenging
                call scav_diag_ls(NMET,BTT,BTTBCK,M)

                !// Add washout of O3 to O3WSCAV
                if (N_STE .gt. 0) &
                     call ctm3_o3scav(BTT,BTTBCK,AIRB,LSTRATAIR_E90B,LSTRATAIR_O3,M)
                call TBGT_IJ (BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_LSSCAV)
 

                !// Convective transport and wet removal
                call CONVW_OSLO(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                              AIRB,BTEM,GAMACB,LSTRATAIR_E90B,DTADV,NOPS,NSUB,M)
                call scav_diag_cn(NMET,BTT,BTTBCK,M)
                call TBGT_IJ (BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_CNSCAV)


                !// BCsnow - master (correct place????)
                if (LBCOC) call bcsnow_master(DTADV,NDAY,NMET,NOPS,NSUB,M)

                !// Gravitational settling of aerosols
                !call aerosolsettling(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY, &
                !     BXZ,BYZ,AIRB,BTEM,DTADV,NDAY,NMET,NOPS,NSUB,M)

                !// Vert. advection & subsidence due to convection
                call DYN2W_OC(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                              AIRB,GAMACB,GAMAB,DTADV,M,MSTEPW)
                WSTEP(M) = WSTEP(M) + MSTEPW

                if (.not.METHANEMIS) call setch4sfc(BTT,BTTBCK,AIRB,M)
                call TBGT_IJ (BTT,BTTBCK,BTTTND,BTTTN0,M,1,NTM,NTND_WADV)

!//-------------------------------------------------------------------------
              end if !// end LDO_ADV
!//-------------------------------------------------------------------------

              !// Put IJ-blocks back into global arrays
              call MPBIND(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                          BTTTND,AIRB,LSTRATAIR_E90B,M)

            end do !// end of IJ-block
!//-------------------------------------------------------------------------
!$omp end do
!$omp end parallel
!//-------------------------------------------------------------------------
            call SYSTEM_CLOCK(count=stime)
            nmp_dt = nmp_dt + real(stime, r8) / rtime
            nuv_dt = nuv_dt - real(stime, r8) / rtime

            if (LDO_ADV) then
!//-------------------------------------------------------------------------
!$omp parallel private(I,J,LUV,AIRUV,ALFA2D,BETA2D,QFV,QFU, &
!$omp                  CHMFLX,MSTEPU,MSTEPV) &
!$omp          shared (DTADV,NTM,USTEP,VSTEP,AIR,STTTNL,NTND_UVADV, &
!$omp                  LFLXDG,N_LZ,N_STE,O3iso1,O3iso2,UFLX1,UFLX2, &
!$omp                  UFLX3,UFLX4,VFLX1,VFLX2,VFLX3,VFLX4) &
!$omp          default(NONE)
!$omp do schedule(dynamic)
!//-------------------------------------------------------------------------
              do LUV = 1,LPAR
                do J = 1,JPAR
                  do I = 1,IPAR
                    AIRUV(I,J) = AIR(I,J,LUV)
                  end do
                end do

                !// Horizontal advection

                !// Transport combines polar-pie with next lower latitude box,
                !// version 4 of qcode 61a.
                call POLES1(AIRUV,ALFA2D,BETA2D,LUV)
                call DYN2UL(DTADV,AIRUV,ALFA2D,LUV,QFU,MSTEPU)
                call DYN2VL(DTADV,AIRUV,BETA2D,LUV,QFV,MSTEPV)
                call POLES2(AIRUV,QFU,QFV,DTADV,LUV)
                !// Alternative transport - more accurate, retains polar box
                !// values. Increases time spent on transport. Instead of the
                !// 4 calls above, compile with p-dyn2-v2.f and p-dyn0-v2.f,
                !// and use these two calls:
                !//  call DYN2UL(DTADV,AIRUV,LUV,QFU,MSTEPU)
                !//  call DYN2VL(DTADV,AIRUV,LUV,QFV,MSTEPV)

                !// STEFLUX: Save horizontal fluxes
                if (LFLXDG .and. N_STE .gt. 0) then
                  !// Isopleth 1
                  call CHEMFLUX(DTADV,LUV,QFU,QFV,N_STE,o3iso1,CHMFLX)
                  UFLX1(:,:,LUV) = CHMFLX(:,1:JPAR,1)
                  VFLX1(:,:,LUV) = CHMFLX(1:IPAR,:,2)

                  !// Isopleth 2
                  call CHEMFLUX(DTADV,LUV,QFU,QFV,N_STE,o3iso2,CHMFLX)
                  UFLX2(:,:,LUV) = CHMFLX(:,1:JPAR,1)
                  VFLX2(:,:,LUV) = CHMFLX(1:IPAR,:,2)

                  !// From E90 tracer
                  call CHEMFLUX_E90(DTADV,LUV,QFU,QFV,CHMFLX,N_STE)
                  UFLX3(:,:,LUV) = CHMFLX(:,1:JPAR,1)
                  VFLX3(:,:,LUV) = CHMFLX(1:IPAR,:,2)

                  !// Using LINOZ and E90 tracer
                  if (LLINOZ) then
                    call CHEMFLUX_E90(DTADV,LUV,QFU,QFV,CHMFLX, N_LZ)
                    UFLX4(:,:,LUV) = CHMFLX(:,1:JPAR,1)
                    VFLX4(:,:,LUV) = CHMFLX(1:IPAR,:,2)
                  end if
                end if !// if (LFLXDG .and. N_STE.gt.0) then

                call TBGT_L(STTTNL,1,NTM,1,IPAR,1,JPAR,LUV,LUV,NTND_UVADV)

                USTEP(LUV) = USTEP(LUV) + MSTEPU
                VSTEP(LUV) = VSTEP(LUV) + MSTEPV

                do J = 1,JPAR
                  do I = 1,IPAR
                    AIR(I,J,LUV) = AIRUV(I,J)
                  end do
                end do

              end do !// do LUV = 1,LPAR
!//-------------------------------------------------------------------------
!$omp end do
!$omp end parallel
!//-------------------------------------------------------------------------

              !// STEFLUX: Sum up horizontal flux
              if (LFLXDG .and. N_STE .gt. 0) then
                call DUMPUVFLUX(UFLX1,VFLX1,1)
                call DUMPUVFLUX(UFLX2,VFLX2,2)
                call DUMPUVFLUX(UFLX3,VFLX3,3)
                if (LLINOZ) call DUMPUVFLUX(UFLX4,VFLX4,4)
              end if

!//-------------------------------------------------------------------------
            end if !// end of LDO_ADV block
!//-------------------------------------------------------------------------
            call SYSTEM_CLOCK(count=stime)
            nuv_dt = nuv_dt + real(stime, r8) / rtime

          end if !// end of LDO_ADV .or. LDO_CHM block
!//-------------------------------------------------------------------------

          !// Accumulate global budgets
          call TBGT_G (STTTNL,BTTTN0, 1)

          GMTAU = GMTAU + DTAULCM

          !// Model integrated time (seconds)
          modelTimeIntegrated = modelTimeIntegrated + DTLCM

          !//---------------------------------------------------------------
        end do !// NSUB                                                 NSUB
        !//-----------------------------------------------------------------

        !---have completed the basic op-split time step = DTAU
        DAY = real(NDAY, r8) + GMTAU / 24._r8

        !// Add to averages (in averages.f90)
        call AVG_ADD2

        call TSER_1 (KSERIES)

        !// Report chemistry & physics
        call REPORTS_CHEMISTRY(NDAY,NMET,NOPS,GMTAU, NDAYI,LNEWM)

        !// NOPS diag for SOA
        if (LSOA) call soa_nopsdiag(NDAY,NMET,NOPS,JDAY,DTOPS)

        if (LDUST) call dustInstBdg(NDAY,NDAY0,NMET,NOPS,DTOPS)

        !// Diag for stomatal conductance and flux
        if(LDDEPmOSaic .and. LSTOM1HRS) then
           call scav_diag_nmet_output_nc(JYEAR,JMON,JDATE,NDAY,NMET,NOPS)
        end if

        !// Done everything for the beginning of new month
        if (LNEWM) LNEWM = .false.

        !// Additional printouts
        !call emissions_seasalt_total()
        !call emissions_ocean_total()

        call SYSTEM_CLOCK(count=stime)
        nops_dt = nops_dt + real(stime, r8) / rtime
        write(6,'(a,i4,2i3,3f8.2)') 'WC NMP/NUV/NOPS [s]: ',NDAY,NMET,NOPS, &
             nmp_dt, nuv_dt, nops_dt
        !//-----------------------------------------------------------------
      end do !// NOPS                                                   NOPS
      !//-------------------------------------------------------------------
      !//---have completed the time step of the met fields



      if (LFLXDG .and. N_STE.ne.0)  then
        call DUMPTRMASS(LSTRATAIR_O3,1)
        call DUMPTRMASS(LSTRATAIR_O3,2)
        call DUMPTRMASS_E90(3,N_STE)
        if (LLINOZ) call DUMPTRMASS_E90(4,N_LZ)
      end if


      do M = 1, MPBLK
         if (nmetTimeIntegrated(M) .ne. DTMET) then
            print*,'pmain: Seriously wrong',M,nmetTimeIntegrated(M),DTMET
            stop
         end if
      end do

      !//-------------------------------------------------------------------
    end do !// NMET                                                     NMET
    !//---------------------------------------------------------------------

    !// STEFLUX: Save daily tropospheric O3 mass
    if (LFLXDG .and. N_STE .gt. 0) then
      !// Create the very initial trop column mass
      if (NDAY .eq. NDAYI)  call SAVETRMASS(0)
      call SAVETRMASS(1)
    end if

    !//---have completed 24-hour (day) step
    GMTAU   = 0._r8
    LEND = NDAY.eq.NDAYE-1

    !// Daily/monthly diagnostics & saves
    !//---------------------------------------------------------------------

    !// Daily output of diagnostics; before calendar update
    !// (diagnostics_general.f90)
    call daily_diag_output(JYEAR,JMON,JDATE,JDAY,NDAY,NDAYI)


    !//---update the calendar to the new, upcoming day
    IDAY = NDAY + 1
    !call CALENDR (IYEAR,IDAY, LLPYR,LFIXMET,MYEAR, &
    !     JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET,SOLDEC,SOLDIS)
    !//---update the calendar to the new, upcoming day plus more
    call calendar(IYEAR,IDAY, LLPYR,LFIXMET,MYEAR, &
         JYEAR,JDAY,JMON,TMON,JDATE,LYEAR,TMET, &
         JYEAR_NEXT,JDAY_NEXT,JMON_NEXT,JDATE_NEXT)
    call get_soldecdis(JDAY,SOLDEC,SOLDIS)
    LNEWM = JDATE.eq.1


    !// STEFLUX: flux dump share calendar JDO_X
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_X, LMXDG,L2)
    if (LMXDG .and. LFLXDG)  then
      call STEBGT_WRITE
      call STEBGT_CLR(IDAY,1)
    end if

    !//---TENDENCIES:  write/reset budget tendencies & time series
    !//---use calendar to trigger diagnostics/saves on/off
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_T, LTBGT,L2)

    !// Initialize 2D budgets files
    if (NDAY .eq. NDAYI) call TBGT_2FILE(0)

    if (LTBGT .or. LEND) then
      !// Write diagnostics: emissions
      call tnd_emis2file(NDAY)

      !// Print out tendencies 0D, 1D and 2D
      call TBGT_P0 (IDAY)
      call TBGT_P1 (IDAY)
      !call TBGT_P2 (IDAY) !// Unformatted output to 24 currently disabled

      !// Print out SALT budget (NDAY, not IDAY)
      if (LSALT) call saltbdg2file(NDAY,NDAYI,NDAY0)

      !// Print out DUST budget
      if (LDUST) call dustbdg2file(NDAY,NDAYI,NDAY0)

      !// Print 2D budgets to file
      call TBGT_2FILE(1)

      if (LSOA) call soa_diag2file(NDAY,NDAYI)

      !// Chemistry budgets
      !call chembud_output(JYEAR,JMON,JDATE,NDAY)

      !// Write diagnostics: scavenging daily totals
      if (LDLYSCAV(1)) then
         call scav_diag_2fileA(JYEAR,JMON,JDAY,NDAY)
      end if
      !// Initialize core tendencies for next average
      call TBGT_G(STTTNL, BTTTN0, 0)

    end if

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
    if (LAVGS .or. LEND) then
      call AVG_WRT_NC4    !// Write averages

      if (L2) then
        call AVG_P1
      end if

      call AVG_CLR2    !// Clear averages
    end if

    !//---restart file save (at least at end of run)
    call CALENDL (JYEAR,JDAY,LYEAR,  JDO_C, LCONT,L2)
    if (LCONT .or. LEND)  then
      call save_restart_file(NDAY+1,NDAYI)
      if (LBCOC) call bcsnow_save_restart(NDAY+1)
    end if

    call SYSTEM_CLOCK(count=stime)
    nday_dt = nday_dt + real(stime, r8) / rtime
    tot_dt = tot_dt + nday_dt !// Summing up total time
    write(6,'(a,i4,f8.1,1x,f10.2,1x,f7.2)') &
         '** WC time NDAY [s] / total [hours/days]:', &
         NDAY,nday_dt, tot_dt/3600._r8, tot_dt/86400._r8
    !//---------------------------------------------------------------------
  end do !// NDAY                                                       NDAY
  !//-----------------------------------------------------------------------

  !//---Final: closeout the model run
  close(21)
  call DATE_AND_TIME (BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),END_TIME)
  call write_log(1, start_time, end_time, tot_dt)
  !//-----------------------------------------------------------------------
end program pmain
!//=========================================================================

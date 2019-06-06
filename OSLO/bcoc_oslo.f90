!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routines for the Oslo Black and Organic Carbon (BCOC).
!//=========================================================================
module bcoc_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: bcoc_oslo
  !// DESCRIPTION: Routines for the Oslo Black and Organic Carbon (BCOC).
  !// ----------------------------------------------------------------------
  !// Oslo CTM3 - Module for BCOC package.
  !// Contains BCOC variables and routines.
  !//   - subroutine bcoc_init
  !//   - subroutine bcoc_master
  !//   - subroutine bcoc_setdrydep
  !//   - subroutine bcoc_vdep2
  !//   - subroutine bcoc_chetinit
  !// BCsnow
  !//   - subroutine bcsnow_init
  !//   - subroutine bcsnow_diagwetrm
  !//   - subroutine bcsnow_save_restart
  !//   - subroutine bcsnow_status
  !//   - subroutine bcsnow_check_snow
  !//   - subroutine bcsnow_getspringsummer
  !//   - subroutine bcsnow_nmet_output
  !//   - subroutine bcsnow_nmet_output_nc
  !//   - subroutine bcsnow_collect_ij
  !//   - subroutine bcsnow_meltevap_ij
  !//   - subroutine bcsnow_seaice_ij
  !//   - subroutine bcsnow_adjustment_ij
  !//   - subroutine bcsnow_master
  !//
  !// Amund Sovde, September 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, &
       TRACER_ID_MAX
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Parameters
  integer, parameter :: max_BC_comps = 6
  integer, parameter :: max_OM_comps = 8

  !// Variables
  !// List of all transport numbers used for BC and OM components
  integer   :: bc_trsp_idx(max_bc_comps)
  integer   :: om_trsp_idx(max_om_comps)

  integer   :: bc_idx_tot, om_idx_tot !// Total numbers of components in use

  !// BC names and component numbers in Oslo CTM3
  !// ----------------------------------------------------------------------
  character(len=8), dimension(max_bc_comps), parameter :: BC_NAMES=(/ &
       'bcBB1fob', 'bcBB1fil', 'bcFF1fob', 'bcFF1fil', &
       'bcBF1fob', 'bcBF1fil'/)
  integer, dimension(max_bc_comps), parameter :: BC_IDS = (/ &
          240,        241,        242,        243, &
          244,        245 /)
  !// Aerosol type information for J-values. Each bin corresponds to a type
  !// entry in Fast-JX tables.
  !// BC mainly absorbs, while OC mainly scatters.
  integer, dimension(max_bc_comps), parameter :: BC_JVAER_TYPE = (/ &
           33,         33,         31,         32, &
           31,         32 /)

  !// OM names and component numbers in Oslo CTM3
  !// ----------------------------------------------------------------------
  character(len=8), dimension(max_om_comps), parameter :: OM_NAMES=(/ &
       'omBB1fob', 'omBB1fil', 'omFF1fob', 'omFF1fil', &
       'omBF1fob', 'omBF1fil', 'omOCNfob', 'omOCNfil'/)
  integer, dimension(max_om_comps), parameter :: OM_IDS = (/ &
          230,        231,        232,        233, &
          234,        235,        236,        237 /)
  !// Aerosol type information for J-values. Each bin corresponds to a type
  !// entry in Fast-JX tables.
  !// BC mainly absorbs, while OC mainly scatters.
  integer, dimension(max_om_comps), parameter :: OM_JVAER_TYPE = (/ &
           33,         33,        -34,        -34, &
          -34,        -34,        -34,        -34/) ! NEED TO CHECK BF!!!!!!!!!!!!!!!!!!!!!!11


  !// Aerosol effective radius (used in oc_fallingaerosols)
  !// Aerosol density (used in oc_fallingaerosols)

  !// Actual BC/OM components in use
  integer, dimension(max_bc_comps) :: BC_IN_USE
  integer, dimension(max_om_comps) :: OM_IN_USE


  !// Dry deposition values (ocean/land/ice)
  real(r8), dimension(max_bc_comps,3) :: VBC
  real(r8), dimension(max_om_comps,3) :: VOM

  !// Transmission values hydrophob -> hydrophil; 12 months
  real(r8), dimension(JPAR,12) :: CHET


  !// BCsnow
  !// ----------------------------------------------------------------------
  !// Switch for turning on/off BCsnow (could be moved to Makefile?)
  logical, parameter :: LBCsnow = .true.
  !// Number of max snow layers
  integer, parameter :: ILMM = 10
  !// Snow layer and BC on snow
  real(r8), dimension(ILMM,IDBLK,JDBLK,MPBLK) :: &
       BSNOWL, BBCFFC, BBCBFC, BBCBIO
  !// Variables for collecting BC in drydep and precipitation
  real(r8), dimension(IDBLK,JDBLK,MPBLK) :: &
       bcsnow_dd_bio, bcsnow_dd_ffc, bcsnow_dd_bfc, &
       bcsnow_prec_bio, bcsnow_prec_ffc, bcsnow_prec_bfc
  !// Number of snow layers
  integer, dimension(IDBLK,JDBLK,MPBLK) :: LSNW_IJ
  !// Hours since last snowfall
  real(r8), dimension(IDBLK,JDBLK,MPBLK) :: tSinceLastSF
  !// Snow in spring
  real(r8), dimension(IDBLK,JDBLK,MPBLK) :: springmelt
  !// Snow layer thickness for a thin layer
  real(r8), parameter  :: thinsnow = 1.e-3_r8
  !// Snow layer thickness for a thin layer
  real(r8), parameter  :: thinsnowthreshold = 1.1_r8*thinsnow
  !// Minimum accepted snow fall, and also minimum accepted layer thickness
  real(r8), parameter :: SFLIM = 1.e-8_r8
  !// Limit of PLANT for sea grid box
  real(r8), parameter :: seafraclim = 0.25_r8
  !// Max snow depth to be used from metdata
  real(r8), parameter :: metdataSDmax = 0.2_r8
  !// DEBUG parameter BCsnow
  logical, parameter :: LDEBUG_BCsnow = .true.
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'bcoc_oslo.f90'
  !// ----------------------------------------------------------------------

  save !// All variables are to be saved.
  private
  public bc_trsp_idx, BC_idx_tot, BC_IDS, OM_IDS, &
       BC_JVAER_TYPE, OM_JVAER_TYPE, BC_IN_USE, OM_IN_USE, &
       bcoc_init, bcoc_master, &
       bcoc_setdrydep, bcoc_vdep2, bcoc_chetinit, &
       LBCsnow,  bcsnow_init, bcsnow_master, bcsnow_diagwetrm, &
       bcsnow_save_restart, bcsnow_status, bcsnow_nmet_output, &
       bcsnow_nmet_output_nc, bcsnow_diag_dd

  !// ----------------------------------------------------------------------
contains
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine bcoc_init()
    !// --------------------------------------------------------------------
    !// Initialize BCOC simulations. Only figures out the transport numbers
    !// and indices for BCOC components.
    !//
    !// Also sets dry deposition values, as they are not changed during the
    !// run, only modified by stability (in oc_drydeposition.f).
    !//
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR_BC, NPAR_OM
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, PLAND
    use cmn_chem, only: TNAME
    use cmn_oslo, only: chem_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Local variables
    integer :: I, II, J, JJ, N, K, MP, idx
    character(len=10) :: tracername !// Name of tracer

    !// Dry deposition; only hardcoded values, for carbonaceous aerosol
    !// on land, sea and ice (ref. Cooke et al., JGR, 22,137, 1999).
    !// They are modified due to stability later. 
    real(r8) :: &
         !// Black Carbon hydrophobic aerosol [cm/s]
         VBCDL = 0.025_r8, VBCDS = 0.025_r8, VBCDI = 0.025_r8, &
         !// Black Carbon hydrophilic aerosol [cm/s]
         VBCWL = 0.025_r8, VBCWS = 0.200_r8, VBCWI = 0.025_r8, &
         !// Organic Carbon hydrophobic aerosol [cm/s]
         VOCDL = 0.025_r8, VOCDS = 0.025_r8, VOCDI = 0.025_r8, &
         !// Organic Carbon hydrophilic aerosol [cm/s]
         VOCWL = 0.025_r8, VOCWS = 0.200_r8, VOCWI = 0.025_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcoc_init'
    !//---------------------------------------------------------------------

    if (NPAR_BC .gt. max_BC_comps) then
       write(6,'(a)') f90file//':'//subr//': NPAR_BC > max_BC_comps'
       stop
    end if
    if (NPAR_OM .gt. max_OM_comps) then
       write(6,'(a)') f90file//':'//subr//': NPAR_OM > max_OM_comps'
       stop
    end if

    !// Find transport numbers to corresponding names
    bc_idx_tot = 0 !// Initialize number of BC components
    BC_IN_USE(:) = 0

    om_idx_tot = 0 !// Initialize number of OM components
    OM_IN_USE(:) = 0

    !// Black carbon names, ids and transport numbers
    !// Assume that BC and OM can be listed arbitrarily, not increasing from
    !// 230-250. The order of local species (K) are listed by BC_IDS,
    !// and BC_IN_USE tells us if the ID is actually used.
    do K = 1, max_bc_comps
       do N = 1, NPAR
          if (TNAME(N) .eq. BC_NAMES(K)) then
             if (chem_idx(N) .ne. BC_IDS(K)) then
                write(6,'(a)') f90file//':'//subr//': Name does not match ID: '//trim(TNAME(N))
                write(6,'(a,2i7)') '  chem_idx/BC_IDS: ', chem_idx(N), BC_IDS(K)
                stop
             end if
             !// Save transport number for this BC component
             bc_trsp_idx(K) = N
             !// Count number of components
             bc_idx_tot = bc_idx_tot + 1
             !// Set in-use flag
             BC_IN_USE(K) = 1
          end if
       end do
    end do

    !// Repeat for organic matter
    do K = 1, max_om_comps
       do N = 1, NPAR
          if (TNAME(N) .eq. OM_NAMES(K)) then
             if (chem_idx(N) .ne. OM_IDS(K)) then
                write(6,'(a)') f90file//':'//subr//': Name does not match ID: '//trim(TNAME(N))
                write(6,'(a,2i7)') '  chem_idx/OM_IDS: ', chem_idx(N), OM_IDS(K)
                stop
             end if
             !// Save transport number for this BC component
             om_trsp_idx(K) = N
             !// Count number of components
             om_idx_tot = om_idx_tot + 1
             !// Set in-use flag
             OM_IN_USE(K) = 1
          end if
       end do
    end do



    !// Check if you read any BC/OC tracers at all
    if (bc_idx_tot .eq. 0) then
       write(*,'(a)') '* bcoc_oslo.f90: No BC components included'
       stop
    else if (bc_idx_tot .gt. max_bc_comps) then
       write(*,'(a,i3)') f90file//':'//subr//': Too many BC tracers found ',bc_idx_tot
       write(*,'(a,i3)') '  Max BC tracers ',max_bc_comps
       stop
    else
       write(*,'(a,i3)') f90file//':'//subr//': Number of BC tracers found ',bc_idx_tot
    end if
    if (om_idx_tot .eq. 0) then
       write(*,'(a)') f90file//':'//subr//': No OM components included'
       stop
    else if (om_idx_tot .gt. max_om_comps) then
       write(*,'(a,i3)') f90file//':'//subr//': Too many OM tracers found ',om_idx_tot
       write(*,'(a,i3)') '  Max OM tracers in bcoc_oslo ',max_om_comps
       stop
    else
       write(*,'(a,i3)') f90file//':'//subr//': Number of OM tracers found ',om_idx_tot
    end if


    !// Set dry deposition values - BC
    do K = 1, max_bc_comps
       if (bc_in_use(K)) then
          if (trim(bc_names(K)) .eq. 'bcBB1fob' .or. &
              trim(bc_names(K)) .eq. 'bcFF1fob' .or. &
              trim(bc_names(K)) .eq. 'bcBF1fob') then
             !// Hydrophobic BC
             VBC(K,1) = VBCDS !// sea
             VBC(K,2) = VBCDL !// land
             VBC(K,3) = VBCDI !// ice
          else if (trim(bc_names(K)) .eq. 'bcBB1fil' .or. &
                   trim(bc_names(K)) .eq. 'bcFF1fil' .or. &
                   trim(bc_names(K)) .eq. 'bcBF1fil') then
             !// Hydrophilic BC
             VBC(K,1) = VBCWS !// sea
             VBC(K,2) = VBCWL !// land
             VBC(K,3) = VBCWI !// ice
          end if
       end if
    end do

    !// Set dry deposition values - OM
    do K = 1, max_om_comps
       if (om_in_use(K)) then
          if (trim(om_names(K)) .eq. 'omBB1fob' .or. &
              trim(om_names(K)) .eq. 'omFF1fob' .or. &
              trim(om_names(K)) .eq. 'omBF1fob' .or. &
              trim(om_names(K)) .eq. 'omOCNfob') then
             !// Hydrophobic OM
             VOM(K,1) = VOCDS !// sea
             VOM(K,2) = VOCDL !// land
             VOM(K,3) = VOCDI !// ice
          else if (trim(om_names(K)) .eq. 'omBB1fil' .or. &
                   trim(om_names(K)) .eq. 'omFF1fil' .or. &
                   trim(om_names(K)) .eq. 'omBF1fil' .or. &
                   trim(om_names(K)) .eq. 'omOCNfil') then
             !// Hydrophilic OM
             VOM(K,1) = VOCWS !// sea
             VOM(K,2) = VOCWL !// land
             VOM(K,3) = VOCWI !// ice
          end if
       end if
    end do

    write(*,'(a)') f90file//':'//subr//': BC / OM dry deposition initialized'

    call bcoc_chetinit()


    !// --------------------------------------------------------------------
  end subroutine bcoc_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcoc_master(BTT,BEMIS,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Integrate Oslo Chemistry package for BC/OC aerosols using the QSSA
    !// integrator.
    !//
    !// Based on OSLO_CHEM_CA in Oslo CTM2, updated for Oslo CTM3.
    !//
    !// Due to the short BCOC code, I have put the whole B-array treatment
    !// in this routine. If the aerosol code grows, it may be wise to
    !// separate out the column treatment, as in oc_tropchem vs pchemc_ij.
    !//
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: LOSLOCSTRAT, LEMISDEP_INCHEM, LOSLOCTROP
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, JDAY, JMON
    use cmn_met, only: ZOFLE
    use cmn_parameters, only: AVOGNR
    use cmn_sfc, only: VDEP
    use cmn_oslo, only: LMTROP, trsp_idx, LVS2ADD2TC, DV_IJ, &
         SCAV_MAP_DRY, SCAV_DD, Xtrsp_idx, XSTT
    use qssa_integrator, only: qssa
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// In/Out parameters
    integer, intent(in) :: MP
    real(r8), intent(in) :: DTCHM
    real(r8), intent(in) :: BEMIS(LPAR,NPAR,IDBLK,JDBLK)
    real(r8), intent(inout) :: BTT(LPAR,NPAR,IDBLK,JDBLK)

    !// Locals
    !// Aerosol mass arrays
    real(r8) :: ZC_BC(max_bc_comps,LPAR)
    real(r8) :: ZC_OM(max_om_comps,LPAR)

    real(r8) :: bcBB1fob_old, bcFF1fob_old, bcBF1fob_old, &
         omBB1fob_old, omFF1fob_old, omBF1fob_old, omOCNfob_old

    real(r8) :: &
         PROD, LOSS, &                         !// production and loss rates
         bcBBhet, bcFFhet, omBBhet, omFFhet, & !// Conversion phob->phil
         DZ, fall_in, fall_out, life_grav, &   !// For gravitational settling
         frac1, rchet
    !real(r8) :: agingOH(LPAR), M_OH
    !integer :: NTR_OH, XNTR_OH

    !// Indices
    integer :: I,J,II,JJ,L,K,N,idx
    !// Indices loop up to
    integer :: NST, NCHEM_ITER, LMTP

    !// Emissions
    real(r8), dimension(max_bc_comps,LPAR) :: EMISX_BC
    real(r8), dimension(max_om_comps,LPAR) :: EMISX_OM

    !// Dry deposition for Oslo method
    real(r8), dimension(max_bc_comps) :: VDEP_BC
    real(r8), dimension(max_om_comps) :: VDEP_OM

    !// QSSA variables
    real(r8) :: ST, QLIN, DTCH, QTEST
    !// Max chemical time step in stratosphere is 15min = 900sec.
    real(r8), parameter :: DTCHM_MAX = 900._r8

    !// Convert kg OH to molec/cm3
    real(r8), parameter :: kg2molecOH = AVOGNR / 17._r8*1.e-3_r8
    !// ------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcoc_master'
    !//---------------------------------------------------------------------

    !// Want maximum time step given above (Oslo CTM2: 15min)

    !// Need a number of loops in chemistry to match the global chemical
    !// time step
    if (DTCHM .lt. DTCHM_MAX) then
       !// Global time step is shorter than DTCHM_MAX
       DTCH = DTCHM
       NCHEM_ITER = 1
    else
       NCHEM_ITER = Int(DTCHM/DTCHM_MAX + 0.5_r8)
       if (mod(DTCHM,DTCHM_MAX) .eq. 0._r8) then
          !// DTCH1 is given by max time step.
          DTCH = DTCHM_MAX
       else
          !// Must change DTCH (i.e. shorter than DTCH_MAX)
          DTCH = DTCHM / real(NCHEM_ITER, r8)
       end if
    end if


    !// Set QSSA parameters for this time step
    ST    = 10._r8 / DTCH
    QLIN  = 0.1_r8 / DTCH
    QTEST = 1.0_r8 / DTCH


    !// Conversion factors from hydrophobic to hydrophilic aerosols
    !C240HET  = 2.7778e-6_r8 ! BC 24% per day; Maria et al., Science 2004
    !C242HET  = 2.4306e-6_r8 ! OC 21% per day; Maria et al., Science 2004
    !C244HET  = 2.7778e-6_r8 ! BC 24% per day; Maria et al., Science 2004
    !C246HET  = 2.4306e-6_r8 ! OC 21% per day; Maria et al., Science 2004


    !// Initialize deposition
    VDEP_BC(:) = 0._r8
    VDEP_OM(:) = 0._r8

    !// Initialize emissions
    EMISX_BC(:,:) = 0._r8
    EMISX_OM(:,:) = 0._r8

    !// Get NTR/XNTR for OH
    !NTR_OH  = trsp_idx(40)
    !XNTR_OH = Xtrsp_idx(40)

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Set transmission values based on latitude and month
       bcBBhet  = CHET(J,JMON)
       bcFFhet  = CHET(J,JMON)
       omBBhet  = CHET(J,JMON)
       omFFhet  = CHET(J,JMON)

       !// Interpolate in time between seasons (not default)
       !// Could use call get_season(JDAY,s1,s2,frac1)
       !rchet = CHET(J,s1) * frac1 + CHET(J,s2) * (1._r8 - frac1)
       !if (rchet .lt. min(CHET(J,s1),CHET(J,s2)) .or. &
       !    rchet .gt. max(CHET(J,s1),CHET(J,s2))) then
       !   print*,'Wrong RCHET',rchet,CHET(J,s1), CHET(J,s2), frac1,season,s1,s2
       !   stop
       !end if
       !bcBBhet  = rchet
       !bcFFhet  = rchet
       !omBBhet  = rchet
       !omFFhet  = rchet

       !// For tropospheric chemistry only, integrate up to top of
       !// troposphere (in this latitude band) as in Oslo CTM2
       !// * Will be overwritten below if stratospheric chemistry is included.
       !// When stratchem is off, do tropchem LVS2ADD2TC levels above "normal":
       LMTP = maxval(LMTROP(:,J)) + LVS2ADD2TC

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1


          !// Tropopause height when using stratospheric chemistry
          if (LOSLOCSTRAT) LMTP = LMTROP(I,J)

          !// Fetch emissions (keep kg/s)
          if (LEMISDEP_INCHEM) then
             !// If source is a separate processes, it is handeled by SOURCE

             !// Emissions
             !// If source is a separate process, this is handeled by SOURCE
             do K = 1, BC_idx_tot 
                if (BC_IN_USE(K) .gt. 0) then
                   !// Find transport number of the BCOC component
                   N = bc_trsp_idx(K)
                   !// Emissions (all levels)
                   do L = 1, LPAR
                      EMISX_BC(K,L) = BEMIS(L,N,II,JJ)
                   end do
                end if
             end do
             do K = 1, OM_idx_tot 
                if (OM_IN_USE(K) .gt. 0) then
                   !// Find transport number of the BCOC component
                   N = om_trsp_idx(K)
                   !// Emissions (all levels)
                   do L = 1, LPAR
                      EMISX_OM(K,L) = BEMIS(L,N,II,JJ)
                   end do
                end if
             end do
          end if

          !// Get OH concentration [molec/cm3] and calculate aging rate
          !// according to Liu et al (2011, JGR, doi:10.1029/2010JD015145)
          !if (LOSLOCTROP) then
          !   if (XNTR_OH .gt. 0) then
          !      do L = 1, LPAR
          !         M_OH = XSTT(L,XNTR_OH,I,J) / DV_IJ(L,II,JJ,MP) * kg2molecOH
          !         agingOH(L) = 4.6e-12_r8 * M_OH + 5.8e-7_r8
          !      end do
          !   else if (NTR_OH .gt. 0) then
          !      do L = 1, LPAR
          !         M_OH = BTT(L,NTR_OH,II,JJ) / DV_IJ(L,II,JJ,MP) * kg2molecOH
          !         agingOH(L) = 4.6e-12_r8 * M_OH + 5.8e-7_r8
          !      end do
          !   end if
          !end if !// if (LOSLOCTROP) then

          !// Get component masses - BC
          do K = 1, BC_idx_tot
             if (BC_IN_USE(K) .gt. 0) then
                !// Find transport number of the BCOC component
                N = bc_trsp_idx(K)
                do L = 1, LPAR
                   !// Do the whole column
                   ZC_BC(K,L) = BTT(L,N,II,JJ)
                end do
             else
                ZC_BC(K,:) = 0._r8
             end if
          end do

          !// Get component masses - OM
          do K = 1, OM_idx_tot
             if (OM_IN_USE(K) .gt. 0) then
                !// Find transport number of the BCOC component
                N = om_trsp_idx(K)
                do L = 1, LPAR
                   !// Do the whole column
                   ZC_OM(K,L) = BTT(L,N,II,JJ)
                end do
             else
                ZC_OM(K,:) = 0._r8
             end if
          end do

          !// check consistency of BL and chemistry iteration steps has been
          !// removed.
          !// --------------------------------------------------------------
          do L = 1, LMTP !// Integrate up to LMTP
             !// -----------------------------------------------------------


             !// Get deposition rates (MUST convert from m/s to 1/s!)
             !// If drydep is a separate processes, it is handeled by DRYDEP
             if (LEMISDEP_INCHEM .and. L .eq. 1) then
                !// Get surface layer thickness
                DZ = ZOFLE(2,I,J) - ZOFLE(1,I,J)
                do K = 1, BC_idx_tot 
                   if (BC_IN_USE(K) .gt. 0) then
                      !// Find transport number of the BC component
                      N = bc_trsp_idx(K)
                      VDEP_BC(K) = VDEP(N,I,J) / DZ
                   end if
                end do
                do K = 1, OM_idx_tot 
                   if (OM_IN_USE(K) .gt. 0) then
                      !// Find transport number of the OM component
                      N = om_trsp_idx(K)
                      VDEP_OM(K) = VDEP(N,I,J) / DZ
                   end if
                end do
             else
                !// Deposition treated by DRYDEP
                VDEP_BC(:) = 0._r8
                VDEP_OM(:) = 0._r8
             end if

             !// Aging rate for this layer
             !if (LOSLOCTROP) then
             !   bcBBhet  = agingOH(L)
             !   bcFFhet  = agingOH(L)
             !   omBBhet  = agingOH(L)
             !   omFFhet  = agingOH(L)
             !end if

             !// Integrate BC aerosols
             !// -----------------------------------------------------------
             do NST = 1, NCHEM_ITER

                !// BCsnow: Store BC dry deposited on snow
                !// Only need one if-test for this, before ZC_BC is changed
                if (L .eq. 1) then
                   if (.not.LEMISDEP_INCHEM) then
                      print*,f90file//':'//subr//': DRYDEP should '// &
                           'not be fetched from bcoc_master: '// &
                           'LEMISDEP_INCHEM = .false.'
                      stop
                   end if

                   !// Snow diagnostic
                   if (LBCsnow) then
                      !// Collect drydep values
                      bcsnow_dd_bio(II,JJ,MP) = bcsnow_dd_bio(II,JJ,MP) &
                           !// bcBB1fob
                           + VDEP_BC(1) * ZC_BC(1,L) * DTCH &
                           !// bcBB1fil
                           + VDEP_BC(2) * ZC_BC(2,L) * DTCH
                      bcsnow_dd_ffc(II,JJ,MP) = bcsnow_dd_ffc(II,JJ,MP) &
                           !// bcFF1fob
                           + VDEP_BC(3) * ZC_BC(3,L) * DTCH &
                           !// bcFF1fil
                           + VDEP_BC(4) * ZC_BC(4,L) * DTCH
                      bcsnow_dd_bfc(II,JJ,MP) = bcsnow_dd_bfc(II,JJ,MP) &
                           !// bcBF1fob
                           + VDEP_BC(5) * ZC_BC(5,L) * DTCH &
                           !// bcBF1fil
                           + VDEP_BC(6) * ZC_BC(6,L) * DTCH
                   end if !// if (LBCsnow) then

                   !// Dry dep for all species
                   do K = 1, BC_idx_tot 
                      if (BC_IN_USE(K) .gt. 0) then
                         !// Find transport number of the BC component
                         N = bc_trsp_idx(K)
                         SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) &
                              + VDEP_BC(K) * DTCH * ZC_BC(K,L)
                         !// Save totals
                         SCAV_DD(N,MP) = SCAV_DD(N,MP) &
                              + VDEP_BC(K) * DTCH * ZC_BC(K,L)
                      end if
                   end do
                   do K = 1, OM_idx_tot 
                      if (OM_IN_USE(K) .gt. 0) then
                         !// Find transport number of the BC component
                         N = om_trsp_idx(K)
                         SCAV_MAP_DRY(N,II,JJ,MP) = SCAV_MAP_DRY(N,II,JJ,MP) &
                              + VDEP_OM(K) * DTCH * ZC_OM(K,L)
                         !// Save totals
                         SCAV_DD(N,MP) = SCAV_DD(N,MP) &
                              + VDEP_OM(K) * DTCH * ZC_OM(K,L)
                      end if
                   end do
                end if !// if (L .eq. 1) then


                !// Save the old values of hydrophobic aerosols
                bcBB1fob_old = ZC_BC(1,L)
                bcFF1fob_old = ZC_BC(3,L)
                bcBF1fob_old = ZC_BC(5,L)

                omBB1fob_old = ZC_OM(1,L)
                omFF1fob_old = ZC_OM(3,L)
                omBF1fob_old = ZC_OM(5,L)
                omOCNfob_old = ZC_OM(7,L)

                !// Black carbon

                !// bcBB1fob (240) BB = Biomass burning
                K = 1
                PROD = EMISX_BC(K,L)
                LOSS = bcBBhet + VDEP_BC(K)
                IF (loss .lt. 0._r8) then
                   print*,f90file//':'//subr//': 240; Negative loss:',I,J,L,LOSS
                end if
                CALL QSSA (240,'bcBB1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))

                !// bcBB1fil (241) BB = Biomass burning
                K = 2
                PROD = EMISX_BC(K,L) + bcBBhet * bcBB1fob_old
                LOSS = VDEP_BC(K)
                CALL QSSA (241,'bcBB1fil',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))

                !// bcFF1fob (242) FF = Fossil fuel
                K = 3
                PROD = EMISX_BC(K,L)
                LOSS = bcFFhet + VDEP_BC(K)
                if (loss .lt. 0._r8) then
                   print*,f90file//':'//subr//': 242; Negative loss:',I,J,L,LOSS
                end if
                CALL QSSA (242,'bcFF1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))

                !// bcFF1fil (243) FF = Fossil fuel
                K = 4
                PROD = EMISX_BC(K,L) + bcFFhet * bcFF1fob_old
                LOSS = VDEP_BC(K)
                CALL QSSA (243,'bcFF1fil',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))

                !// bcBF1fob (244) BF = Bio fuel
                K = 5
                PROD = EMISX_BC(K,L)
                LOSS = bcFFhet + VDEP_BC(K)
                if (loss .lt. 0._r8) then
                   print*,f90file//':'//subr//': 244; Negative loss:',I,J,L,LOSS
                end if
                CALL QSSA (244,'bcFF1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))
                if (bcBF1fob_old .gt. 1.) then
                   if (abs(ZC_BC(K,L) - bcBF1fob_old)/bcBF1fob_old .gt. 1000.) then
                      print*,'bcoc_master: large change',PROD,bcffhet,vdep_bc(K),ZC_BC(K,L), bcBF1fob_old
                      stop
                   end if
                end if

                !// bcBF1fil (245) BF = Bio fuel
                K = 6
                PROD = EMISX_BC(K,L) + bcFFhet * bcBF1fob_old
                LOSS = VDEP_BC(K)
                CALL QSSA (245,'bcFF1fil',DTCH,QLIN,ST,PROD,LOSS,ZC_BC(K,L))


                !// Organic matter

                !// omBB1fob (230)
                K = 1
                PROD = EMISX_OM(K,L)
                LOSS = omBBhet + VDEP_OM(K)
                CALL QSSA (230,'ocBB1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omBB1fil (231)
                K = 2
                PROD = EMISX_OM(K,L) + omBBhet * omBB1fob_old
                LOSS = VDEP_OM(K)
                CALL QSSA (231,'omBB1fil',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omFF1fob (232)
                K = 3
                PROD = EMISX_OM(K,L)
                LOSS = omFFhet + VDEP_OM(K)
                CALL QSSA (232,'omFF1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omFF1fil (233)
                K = 4
                PROD = EMISX_OM(K,L) + omFFhet * omFF1fob_old
                LOSS = VDEP_OM(K)
                CALL QSSA (233,'OCFFPHIL',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omBF1fob (234)
                K = 5
                PROD = EMISX_OM(K,L)
                LOSS = omFFhet + VDEP_OM(K)
                CALL QSSA (232,'omFF1fob',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omBF1fil (235)
                K = 6
                PROD = EMISX_OM(K,L) + omFFhet * omBF1fob_old
                LOSS = VDEP_OM(K)
                CALL QSSA (233,'OCFFPHIL',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omOCNfob (236) Oceanic organic carbon
                K = 7
                PROD = EMISX_OM(K,L)
                LOSS = omFFhet + VDEP_OM(K)
                CALL QSSA (234,'omOCNfob',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

                !// omOCNfil (237) Oceanic organic carbon
                K = 8
                PROD = EMISX_OM(K,L) + omFFhet * omOCNfob_old
                LOSS = VDEP_OM(K)
                CALL QSSA (235,'omOCNfil',DTCH,QLIN,ST,PROD,LOSS,ZC_OM(K,L))

             end do !// Do NST = 1, NCHEM_ITER
             !// -----------------------------------------------------------

             !// Save new values - BC
             do K = 1, BC_idx_tot
                if (BC_IN_USE(K) .gt. 0) then
                   !// Transport number of the BCOC component
                   N = bc_trsp_idx(K)
                   BTT(:,N,II,JJ) = ZC_BC(K,:)
                end if
             end do

             !// Save new values - OM
             do K = 1, OM_idx_tot
                if (OM_IN_USE(K) .gt. 0) then
                   !// Transport number of the BCOC component
                   N = om_trsp_idx(K)
                   BTT(:,N,II,JJ) = ZC_OM(K,:)
                end if
             end do

             !// -----------------------------------------------------------
          end do !// Do L=1,LMTP !// Integrate up to LMTP
          !// --------------------------------------------------------------


          !// Stratosphere: Gravitational settling; loop from top
          !// --------------------------------------------------------------
          !if (.false.) then
          do K = 1, BC_idx_tot
             if (BC_IN_USE(K) .le. 0) cycle
             !// Tracer K is in use; get transport number
             N = bc_trsp_idx(K)

             fall_in = 0._r8
             do L = LPAR, LMTP, -1

                !// Assume fall speed 2.5cm/hr = 6.944d-6m/s, as for SO4.
                !// 1/6.944d-6m/s = 144000s/m, so we get a life time of
                !// DZ / 6.944d-6 = DZ * 144000 (in seconds)
                life_grav = (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) * 144000._r8
                !// Fractional loss of BTT (will always be positive)
                LOSS = max(DTCH / life_grav,0._r8)

                !// Falling down
                fall_out = BTT(L,N,II,JJ) * LOSS

                BTT(L,N,II,JJ) = BTT(L,N,II,JJ) * (1._r8 - LOSS) + fall_in
                !// Coming in to the layer below
                fall_in = fall_out

             end do !// do L = LPAR, LMTP, -1
          end do !// do K = 1, BC_idx_tot

          do K = 1, OM_idx_tot
             if (OM_IN_USE(K) .le. 0) cycle
             !// Tracer K is in use; get transport number
             N = om_trsp_idx(K)

             fall_in = 0._r8
             do L = LPAR, LMTP, -1

                !// Assume fall speed 2.5cm/hr = 6.944d-6m/s, as for SO4.
                !// 1/6.944d-6m/s = 144000s/m, so we get a life time of
                !// DZ / 6.944d-6 = DZ * 144000 (in seconds)
                life_grav = (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) * 144000._r8
                !// Fractional loss of BTT (will always be positive)
                LOSS = max(DTCH / life_grav,0._r8)

                !// Falling down
                fall_out = BTT(L,N,II,JJ) * LOSS

                BTT(L,N,II,JJ) = BTT(L,N,II,JJ) * (1._r8 - LOSS) + fall_in
                !// Coming in to the layer below
                fall_in = fall_out

             end do !// do L = LPAR, LMTP, -1
          end do !// do K = 1, OM_idx_tot
          !// --------------------------------------------------------------
          !end if

       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)

    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine bcoc_master
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcoc_setdrydep(VDEP,RFR,MP)
    !// --------------------------------------------------------------------
    !// Set dry deposition for BCOC, in a given IJ-block.
    !// Called from subroutine oc_setdrydep, which also modifies deposition
    !// rate by stability.
    !//
    !// Amund Sovde, October 2013
    !//   Now takes land fractions (RFR) into account, and thereby also
    !//   treating ocean ice as ice. Old method (as in CTM2) used ocean
    !//   uptake over ice-covered ocean.
    !// Amund Sovde, September 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TNAME
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: RFR(5,IDBLK,JDBLK) !// Land fractions
    !// Output
    real(r8), intent(inout) :: VDEP(NPAR,IPAR,JPAR)

    !// Locals
    integer :: I,II,J,JJ,N, K
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Construct dry deposition based on land type fractions RFR
          !// Ocean: RFR(1,II,JJ)
          !// Land: RFR(2:4,II,JJ)
          !// Ice covered land/ocean: RFR(5,II,JJ)
          !// VDEP(N,I,J) = ( RFR(1,II,JJ)*VBCOC(K,1)
          !//                 + sum(RFR(2:4,II,JJ))*VBCOC(K,2)
          !//                 + RFR(5,II,JJ)*VBCOC(K,3) ) * 1.d-2
          do K = 1, BC_idx_tot
            if (BC_IN_USE(K) .gt. 0) then
               !// Find transport number of the BCOC component
               N = bc_trsp_idx(K)

               !// Convert from cm/s to m/s
               !// The deposition VBC is already sorted in the pre-defined
               !// BC order in this module (see the top if you are uncertain).
               VDEP(N,I,J) = ( RFR(1,II,JJ) * VBC(K,1) &
                               + sum(RFR(2:4,II,JJ)) * VBC(K,2) &
                               + RFR(5,II,JJ) * VBC(K,3) ) * 1.e-2_r8
            end if
          end do !// do K = 1, BC_idx_tot 

          do K = 1, OM_idx_tot
            if (OM_IN_USE(K) .gt. 0) then
               !// Find transport number of the BCOC component
               N = om_trsp_idx(K)

               !// Convert from cm/s to m/s
               !// The deposition VOM is already sorted in the pre-defined
               !// OC order in this module (see the top if you are uncertain).
               VDEP(N,I,J) = ( RFR(1,II,JJ) * VOM(K,1) &
                               + sum(RFR(2:4,II,JJ)) * VOM(K,2) &
                               + RFR(5,II,JJ) * VOM(K,3) ) * 1.e-2_r8
            end if
          end do !// do K = 1, OM_idx_tot 


       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine bcoc_setdrydep
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcoc_vdep2(VDEP,SCALESTABILITY,MP)
    !// --------------------------------------------------------------------
    !// Set dry deposition for BCOC, in a given IJ-block.
    !// Called from subroutine oc_setdrydep.
    !// This method is based on the treatment in EMEP model, Simpson etal
    !// (2012), ACP, doi:10.1054/acp-12-7825-2012, refered to as EMEP2012
    !// in this routine.
    !//
    !// Amund Sovde, February 2014
    !// --------------------------------------------------------------------
    use cmn_ctm, only: JMON, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, &
         YDGRD, PLAND
    use cmn_chem, only: TNAME
    use cmn_met, only: PRECLS, PRECCNV, MO_LENGTH, USTR, SFT, CI, SD, CLDFR
    use cmn_sfc, only: LANDUSE_IDX, landSurfTypeFrac, NVGPAR, LAI, NLCAT, &
         DDEP_PAR
    use utilities_oslo, only: landfrac2mosaic, set_vegetation_height
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    !// Output
    real(r8), intent(inout) :: VDEP(NPAR,IPAR,JPAR)
    integer, intent(inout) :: SCALESTABILITY(NPAR)

    !// Locals
    integer :: I,II,J,JJ,N, K, NN
    real(r8) :: SAI1, MOL, USR, RAIN, T2M, &
         a1L, a1W, a1I, a1Lfor, amol, Vtot, WETFRAC, fice, snowD, &
         focean
    real(r8),dimension(NLCAT) :: FL, VDLCAT, VEGH, fsnowC, tmpVEGH

    !// Ustar mean for year 2006 0.293m = 29.3cm
    real(r8), parameter :: ZmeanUSR = 1._r8/29.3_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcoc_vdep2'
    !// --------------------------------------------------------------------

    !// EMEP categories
    !//  1. Forrests, Mediterranean scrub
    !//  2. Crops
    !//  3. Moorland (savanna++)
    !//  4. Grassland
    !//  5. Wetlands
    !//  6. Tundra
    !//  7. Desert
    !//  8. Water
    !//  9. Urban
    !// 10. Ice+snow (strictly, EMEP has this as category 9 and urban as 10)
    !// mOSaic categories (reduced from Simpson et al., 2012)
    !//  1. Needleleaftree (temperated/boreal)
    !//  2. Deciduoustree (temperated/boral)
    !//  3. Needleleaftree (mediterranean)
    !//  4. Broadleaftree (mediterranean)
    !//  5. Crops <a,b,c>
    !//  6. Moorland (savanna++)
    !//  7. Grassland
    !//  8. Scrubs (med.)
    !//  9. Wetlands
    !//  10. Tundra
    !//  11. Desert
    !//  12. Water
    !//  13. Urban
    !//  14. Ice/Snow - is treated seperately

    !// Vegetation height
    VEGH    = DDEP_PAR(20,:) !// Will be modified latitude dependently
    tmpVEGH = DDEP_PAR(20,:) !// Temporare vegetation height
    
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      
      !// Set latitude dependent vegetation height
      !// The function is based on the latitude dependent description in Simpson et al. (2012) 
      !// and modified towards tropics (see utilities_oslo.f90).
      do NN=1, 4 ![5]
         call set_vegetation_height(tmpVEGH(NN),YDGRD(J),VEGH(NN))
      end do
      !// [Grass and] scrubs
      ! do NN=7, 8
      call set_vegetation_height(tmpVEGH(8),YDGRD(J),VEGH(8))
      ! end do
      
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Basis of method:
        !//   Vd/USTAR = a1
        !//   a1=0.002 (Wesely, 1985)
        !//     for Obukhov length L<0 also multiply with (1 + (-300/L)^(2/3))
        !//     with limit 1/L>-0.04; L<-25m
        !//   Forests: a1=0.008 SAI/10
        !// Instead we find a1 from our deposition velocities and assume
        !// they apply to the mean USTAR, which is 0.293m.

        !// Method needs surface area index for forests, i.e. LAI+1:
        SAI1 = max(0._r8, LAI(I,J,JMON)) + 1._r8

        !// Meteorological parameters
        RAIN = PRECLS(I,J,1) + PRECCNV(I,J,1) !// Rain at surface
        T2M = SFT(I,J)                        !// Surface (2m) temperature
        MOL    = MO_LENGTH(II,JJ,MP)          !// Obukhov lenght
        USR    = USTR(I,J)                    !// Friction velocity
        if (USR .le. 0._r8) USR = 5.e-3_r8    !// USR should not be zero

        !// Adjustment for negative Obukhov length
        if (MOL .lt. 0._r8) then
           !// Restriction: 1/L> -0.04, i.e. L<-25.
           if (MOL .le. -25._r8) then
              amol = 1._r8 + (-300._r8/MOL)**(2._r8/3._r8)
           else
              !// -300/-25=12
              amol = 1._r8 + (12._r8)**(2._r8/3._r8)
           end if
        else
           amol = 1._r8
        end if


        !// Set land fractions
        call landfrac2mosaic(FL,NLCAT,landSurfTypeFrac(:,I,J), &
             NVGPAR, YDGRD(J), LANDUSE_IDX)


        !// Snow cover for each vegetation type, calculate snowD as meter snow.
        !// See oc_drydeposition.f90 for further explanations.
        if (PLAND(I,J) .eq. 1._r8) then
           snowD = SD(I,J) * 10._r8
        else if (PLAND(I,J) .gt. 0._r8) then
           if (CI(I,J) .eq. 0._r8) then
              !// Assume snow depth covers only land
              snowD = SD(I,J) / PLAND(I,J) * 10._r8
           else
              if (CI(I,J) .ge. FL(12)) then
                 !// More sea ice than water; should not occur, but we accept
                 !// it for now and assume ocean to be fully covered by ice.
                 snowD = SD(I,J) * 10._r8
              else
                 !// Part of water is ice covered. The part of the gridbox
                 !// which is NOT covered by snow is (focean - CI(I,J)),
                 !// so we adjust snow depth for the covered part.
                 snowD = SD(I,J) / (1._r8 - (FL(12) - CI(I,J))) * 10._r8
              end if
           end if
        else
           !// All is water
           if (CI(I,J) .gt. 0._r8) then
              snowD = SD(I,J) / CI(I,J) * 10._r8 !// Snow depth on ice
           else
              snowD = 0._r8
           end if
        end if

        !// Find snow depth for all vegetation types
        do NN = 1, NLCAT-1
           if (snowD .eq. 0._r8) then
              fsnowC(NN) = 0._r8 !// No snow present
           else
              if (VEGH(NN) .eq. 0._r8) then
                 !// Snow and zero vegetation height is assumed to be snow
                 !// covered.
                 fsnowC(NN) = 1._r8
              else
                 !// Calculate snow cover for vegetation type
                 fsnowC(NN) = min(1._r8, snowD / (0.1_r8 * VEGH(NN)))
              end if
           end if
        end do
        !// Ocean snow cover is done below where fice is calculated.
        !// Snow land is of course snow covered
        fsnowC(14)= 1._r8


        !// Sea ice - check if we need to reduce fraction of ocean
        if (CI(I,J) .eq. 0._r8) then
           fice = 0._r8
        else
           !// Cannot have more sea ice than sea, so the fraction of sea
           !// covered by ice is
           if (FL(12) .gt. 0._r8) then
              fice = min(CI(I,J)/FL(12), 1._r8)
           else
              fice = 0._r8
           end if
        end if
        fsnowC(12) = fice

        !// Taking snow cover and sea ice into account.
        !// Above 0C, snow should be treated as a wet surface.
        !// Because ice is probably not wet until a few degrees above 0C, we
        !// should perhaps use 1-2C, but this must be tested.

        !// Change snow covered land fractions to snow/ice when T<0C
        if (T2M .le. 273.15_r8) then
           !// Check for snow cover at cold temperatures
           do NN = 1, NLCAT-1
              if (fsnowC(NN) .gt. 0._r8) then
                 !// Move the snow covered fraction of each land type to FL(14)
                 FL(14) = FL(14) + FL(NN)*fsnowC(NN)
                 FL(NN) = FL(NN) * (1._r8 - fsnowC(NN))
              end if
           end do
        else
           !// Snow cover for T>0C should be treated as wet surface.
           !// To simplify, we put all wet surfaces into FL(12).
           do NN = 1, 11
              if (fsnowC(NN) .gt. 0._r8) then
                 !// Move the snow covered fraction of each land type to FL(12)
                 FL(12) = FL(12) + FL(NN)*fsnowC(NN)
                 FL(NN) = FL(NN) * (1._r8 - fsnowC(NN))
              end if
           end do
           do NN = NLCAT-1, NLCAT
              if (fsnowC(NN) .gt. 0._r8) then
                 FL(12) = FL(12) + FL(NN)*fsnowC(NN)
                 FL(NN) = FL(NN) * (1._r8 - fsnowC(NN))
              end if
           end do
        end if


        !// Final check: wetland is wet when T>0C, but assume ice when T<0C
        if (T2M .lt. 273.15_r8 .and. FL(9) .gt. 0._r8) then
           FL(14) = FL(14) + FL(9)
           FL(9) = 0._r8
        end if



        !// Rainy surface, i.e. wet surface
        !// This fraction should be assumed to be distributed equally
        !// on all land types, in lack of other information
        if (RAIN .gt. 0._r8 .and. T2M .gt. 273.15_r8) then
           if (PLAND(I,J) .gt. 0._r8) then
              WETFRAC = maxval(CLDFR(I,J,:))
           else
              WETFRAC = 0._r8
           end if
        else
           WETFRAC = 0._r8
        end if



        !// Construct dry deposition based on land types - BC
        do K = 1, BC_idx_tot
          if (BC_IN_USE(K) .gt. 0) then
            !// Find transport number of the BCOC component
            N = bc_trsp_idx(K)

            !// Estimate a1 over ocean/land/ice using mean Ustar [cm/s], giving
            !// unitless a1. Multiplying by USR will give m/s. 
            a1W = VBC(K,1) * ZmeanUSR
            a1L = VBC(K,2) * ZmeanUSR
            a1I = VBC(K,3) * ZmeanUSR

            !// Forest; limit forest to our land value
            a1Lfor = max(0.008_r8 * SAI1*0.1_r8, a1L)

            !// If rain, perhaps wet forest/shrubs should increase?
            !// Try increasing 1.5 times.
            VDLCAT(1:4) = (a1Lfor * (1._r8 - WETFRAC) &
                          + 1.5_r8 * a1Lfor * WETFRAC) * USR * amol
            VDLCAT(8)   = (a1Lfor * (1._r8 - WETFRAC) &
                          + 1.5_r8 * a1Lfor * WETFRAC) * USR * amol


            !// Assume wetland is wet. Have already checked it for temperature,
            !// assuming T<0C is to treated as snow/ice.
            !// Also multiply with amol
            VDLCAT(9)  = a1W * USR * amol
            !// Ocean is wet
            VDLCAT(12) = a1W * USR * amol

            !// Non-wet surfaces - weight with WETFRAC
            !// Also multiply with amol
            VDLCAT(5:7)   = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            VDLCAT(10:11) = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            VDLCAT(13)    = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            
            !// Ice - weight with WETFRAC
            VDLCAT(14)    = (a1I * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol


            !// Make average velocity
            Vtot = 0._r8
            do NN = 1, NLCAT
               Vtot = Vtot + VDLCAT(NN) * FL(NN)
            end do
            VDEP(N,I,J) = Vtot

          end if !// if (BC_IN_USE(K) .gt. 0) then
        end do !// do K = 1, BC_idx_tot 


        !// Construct dry deposition based on land types - OM
        do K = 1, OM_idx_tot
          if (OM_IN_USE(K) .gt. 0) then
            !// Find transport number of the BCOC component
            N = om_trsp_idx(K)

            !// Estimate a1 over ocean/land/ice using mean Ustar [cm/s], giving
            !// unitless a1. Multiplying by USR will give m/s. 
            a1W = VOM(K,1) * ZmeanUSR
            a1L = VOM(K,2) * ZmeanUSR
            a1I = VOM(K,3) * ZmeanUSR

            !// Forest; limit forest to our land value
            a1Lfor = max(0.008_r8 * SAI1*0.1_r8, a1L)

            !// If rain, perhaps wet forest/shrubs should increase?
            !// Try increasing 1.5 times.
            VDLCAT(1:4) = (a1Lfor * (1._r8 - WETFRAC) &
                          + 1.5_r8 * a1Lfor * WETFRAC) * USR * amol
            VDLCAT(8)   = (a1Lfor * (1._r8 - WETFRAC) &
                          + 1.5_r8 * a1Lfor * WETFRAC) * USR * amol


            !// Assume wetland is wet. Have already checked it for temperature,
            !// assuming T<0C is to treated as snow/ice.
            !// Also multiply with amol
            VDLCAT(9)  = a1W * USR * amol
            !// Ocean is wet
            VDLCAT(12) = a1W * USR * amol

            !// Non-wet surfaces - weight with WETFRAC
            !// Also multiply with amol
            VDLCAT(5:7)   = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            VDLCAT(10:11) = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            VDLCAT(13)    = (a1L * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol
            
            !// Ice - weight with WETFRAC
            VDLCAT(14)    = (a1I * (1._r8 - WETFRAC) + a1W * WETFRAC) * USR * amol


            !// Make average velocity
            Vtot = 0._r8
            do NN = 1, NLCAT
               Vtot = Vtot + VDLCAT(NN) * FL(NN)
            end do
            VDEP(N,I,J) = Vtot

          end if !// if (OM_IN_USE(K) .gt. 0) then
        end do !// do K = 1, OM_idx_tot 

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)

    do K = 1, BC_idx_tot
      if (BC_IN_USE(K) .gt. 0) then
        !// Find transport number of the BC component
        N = bc_trsp_idx(K)
        SCALESTABILITY(N) = 0
      end if
    end do
    do K = 1, BC_idx_tot
      if (BC_IN_USE(K) .gt. 0) then
        !// Find transport number of the OM component
        N = om_trsp_idx(K)
        SCALESTABILITY(N) = 0
      end if
    end do

    !// --------------------------------------------------------------------
  end subroutine bcoc_vdep2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcoc_chetinit()
    !// --------------------------------------------------------------------
    !// Initialize CHET. Now monthly data.
    !//
    !// Amund Sovde Haslerud, September 2017,
    !// Amund Sovde, 2009
    !// --------------------------------------------------------------------
    use cmn_ctm, only: YDEDG
    use cmn_parameters, only: ZPI180, A0, CPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID_Y
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input

    !// Locals
    integer :: J, M
    !// For file reading and interpolation
    integer,parameter :: JG=64
    real(r8) :: CHET_T42(JG,12), CTMP(JG)
    character(len=80) :: FNAME
    logical :: fnr_ok, ex
    integer :: ifnr, io_err
    integer :: JGB2
    real(r8)  :: WGAULAT(JG),WGAUWT(JG), WYEDG1P1(JG+1),YBEDG(JG+1)
    real(r8)  :: awT42(JG), awCUR(JPAR)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcoc_chetinit'
    !//---------------------------------------------------------------------

    !// Find non-used file number for input file
    fnr_ok = .true.
    ifnr = 8
    do while (fnr_ok)
       ifnr = ifnr + 1
       inquire(ifnr,opened=fnr_ok)
    end do

    fname = 'Indata_CTM3/bcoc_trans_t42_monthly_v2016.dat'
    !// Read original T42 CHET and interpolate
    inquire (file=fname,exist=ex)
    if (ex) then
       open(ifnr,File=fname,Form='FORMATTED',Status='OLD')
    else
       write(6,'(a)') f90file//':'//subr//': no such file '//trim(fname)
       stop
    end if
    do J = 1, JG
       read(ifnr,'(12F13.9)') CHET_T42(J,:)
    end do
    close(ifnr)


    !// Set YBEDG
    call GAUSST2(JG,WGAULAT,WGAUWT)
    JGB2 = JG/2
    WYEDG1P1(1) = -1._r8
    do J = 2,JGB2
       WYEDG1P1(J) = WYEDG1P1(J-1) + WGAUWT(J-1)
    end do
    WYEDG1P1(JGB2+1) = 0._r8
    do J = JGB2+2,JG+1
       WYEDG1P1(J) = -WYEDG1P1(JG+2-J)
    end do
    do J = 1,JG+1
       YBEDG(J) = ZPI180 * asin(WYEDG1P1(J))
    end do

    !// Set CHET values, unit [1/s]
    if (jpar .eq. jg) then
       !// Current resolution is T42
       do J = 1, JPAR
          CHET(J,:) = CHET_T42(J,:)
       end do
    else
       !// Data are per seconds, but must be multiplied with area to be
       !// interpolated correctly, otherwise a T42 grid box having value
       !// 8 will be divided into the smaller grid boxes.
       !// 
       !// Grid box areas for T42 (dXBEDG = 2.8125) are defined by:
       !// AREA = A0*A0 * CPI180*dXBEDG &
       !//      * (sin(CPI180*YBEDG(J+1)) - sin(CPI180*YBEDG(J)))
       !// but in the case of interpolating only in meridional direction, the
       !// horizontal extent must be the whole latitude band (2PI):
       !// AREAband = A0*A0 * 2 * PI
       !//            * (sin(CPI180*YBEDG(J+1)) - sin(CPI180*YBEDG(J)))
       !// After interpolation, we divide by the current area of each band:
       !// CURRENT = A0*A0 * 2 * PI
       !//            * (sin(CPI180*YDEDG(J+1)) - sin(CPI180*YDEDG(J)))
       !// Because A0*A0*2*PI cancel, and we can use weighting factors
       !//   awT42(J) = (sin(CPI180*YBEDG(J+1)) - sin(CPI180*YBEDG(J)))
       !//   awCUR(J) = (sin(CPI180*YDEDG(J+1)) - sin(CPI180*YDEDG(J)))
       do J = 1, JG
          !// Weighting factor for T42 latitude band
          awT42(J) = (sin(CPI180*YBEDG(J+1)) - sin(CPI180*YBEDG(J)))
       end do
       do J = 1, JPAR
          !// Weighting factor for current resolution latitude band
          awCUR(J) = (sin(CPI180*YDEDG(J+1)) - sin(CPI180*YDEDG(J)))
       end do

       !// Interpolate to current resolution
       do M = 1, 12
          !// Must multiply with area weighting fraction of latitude band
          CTMP(:) = CHET_T42(:,M) * awT42(:)
          call E_GRID_Y(CTMP(:),YBEDG,JG, CHET(:,M),YDEDG,JPAR)
          !// And divide by new area weighting fraction of latitude band
          CHET(:,M) = CHET(:,M) / awCUR(:)
       end do
    end if
    write(*,'(a)') f90file//':'//subr//': BCOC CHET initialized as'
    !// Print out values
    write(6,'(6x,12(5x,a1,5x))') 'J','F','M','A','M','J','J','A','S','O','N','D'
    do J = 1, JPAR
       write(6,'(2x,i3,1x,12es11.4)') J, CHET(J,:)
    end do

    !// --------------------------------------------------------------------
  end subroutine bcoc_chetinit
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_season(JDAY,s1,s2,frac1)
    !// --------------------------------------------------------------------
    !// Find seasons s1 and s2, the seasons closest to current day.
    !// Temporal fraction from s1 is given by frac1.
    !//
    !// February 2014, Amund Sovde
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in)  :: JDAY
    integer, intent(out) :: s1, s2
    real(r8), intent(out)  :: frac1

    !// Mid season DJF: JDAY=15
    !// Mid season MAM: JDAY=105
    !// Mid season JJA: JDAY=197
    !// Mid season SON: JDAY=289
    integer, parameter :: midDJF=15, midMAM=105, midJJA=198, midSON=289
    integer, parameter :: sMAM=60, sJJA=152, sSON=244, sDJF=335
    !// --------------------------------------------------------------------

    if(JDAY .ge. sMAM .and. JDAY .lt. sJJA) then
       !// MAM
       if (JDAY .le. midMAM) then
          s1 = 1 !// DJF
          s2 = 2 !// MAM
          frac1 = 1._r8 - real(JDAY - midDJF, r8)/real(midMAM - midDJF, r8)
       else
          s1 = 2 !// MAM
          s2 = 3 !// JJA
          frac1 = 1._r8 - real(JDAY - midMAM, r8)/real(midJJA - midMAM, r8)
       end if
    else if(JDAY .ge. sJJA .and. JDAY .lt. sSON) then
       !// JJA
       if (JDAY .le. midJJA) then
          s1 = 2 !// MAM
          s2 = 3 !// JJA
          frac1 = 1._r8 - real(JDAY - midMAM, r8)/real(midJJA - midMAM, r8)
       else
          s1 = 3 !// JJA
          s2 = 4 !// SON
          frac1 = 1._r8 - real(JDAY - midJJA, r8)/real(midSON - midJJA, r8)
       end if
    else if(JDAY .ge. sSON .and. JDAY .lt. sDJF) then
       !// SON
       if (JDAY .le. midSON) then
          s1 = 3 !// JJA
          s2 = 4 !// SON
          frac1 = 1._r8 - real(JDAY - midJJA, r8)/real(midSON - midJJA, r8)
       else
          s1 = 4 !// SON
          s2 = 1 !// DJF
          frac1 = 1._r8 - real(JDAY - midSON, r8)/real(midDJF+365 - midSON, r8)
       end if
    else
       !// DJF (midDJF is 15)
       if (JDAY .lt. sMAM .and. JDAY .gt. midDJF) then
          s1 = 1 !// DJF
          s2 = 2 !// MAM
          frac1 = 1._r8 - real(JDAY - midDJF, r8)/real(midMAM - midDJF, r8)
       else
          s1 = 4 !// SON
          s2 = 1 !// DJF
          if (JDAY .lt. sMAM) then
             !// Day 1-59
             frac1 = 1._r8 - real(JDAY+365 - midSON, r8) &
                             / real(midDJF+365 - midSON, r8)
          else
             !// Day 335-365 (if leap 366 assume 1, i.e. same frac1 as 1)
             frac1 = 1._r8 - real(JDAY - midSON, r8) &
                             / real(midDJF+365 - midSON, r8)
          end if
       end if
    end if

    !// --------------------------------------------------------------------
  end subroutine get_season
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  !//                      BCsnow BCsnow BCsnow
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  !// Pt. 1 - initialisation
  !//   1. Read BSNOWL, BBCFFC, BBCBFC and BBCBIO from file
  !//   2. Initialise other stuff
  !// Pt. 1 - each NMET, build/evaporate snow layers
  !//   1. Snowfall
  !//   2. Evaporation
  !//   3. Diagnostic of 1-2
  !// Pt. 2 - each chem step
  !//   1. Dry deposition
  !//   2. Wet deposition
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine bcsnow_init(startday)
    !// --------------------------------------------------------------------
    !// Initialize BCsnow.
    !// Checking for tracers and reading restart-file specific for
    !// BCsnow.
    !//
    !// Amund Sovde, April 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: XDGRD, YDGRD, LCONT, &
         MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: SD
    use cmn_oslo, only: trsp_idx
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: startday
    !// Locals
    logical :: LTEST
    integer :: fnr, dday, I, J, II, JJ, MP, L
    real(r8) :: INFIELD3(ILMM,IPAR,JPAR), INFIELD2(IPAR,JPAR)
    integer :: INSNWIJ(IPAR,JPAR)
    character(len=80) :: filename
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_init'
    !//---------------------------------------------------------------------

    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return

    !// Test for all components necessary: 240, 241, 245, 246
    LTEST = trsp_idx(240) .gt. 0 .and. trsp_idx(241) .gt. 0 .and. &
            trsp_idx(242) .gt. 0 .and. trsp_idx(243) .gt. 0 .and. &
            trsp_idx(244) .gt. 0 .and. trsp_idx(244) .gt. 0
    if (.not.LTEST) then
       write(*,'(a)') f90file//':'//subr//': Not enough tracers to run BCsnow!'
       write(*,'(a)') '    trsp_idx(240)',trsp_idx(240)
       write(*,'(a)') '    trsp_idx(241)',trsp_idx(241)
       write(*,'(a)') '    trsp_idx(242)',trsp_idx(242)
       write(*,'(a)') '    trsp_idx(243)',trsp_idx(243)
       write(*,'(a)') '    trsp_idx(244)',trsp_idx(244)
       write(*,'(a)') '    trsp_idx(245)',trsp_idx(245)
       stop
    end if


    !// Continue from previous run?
    fnr = get_free_fileid()

    !// Does the file exist?
    filename = 'restart_bcsnow.sav'
    inquire(file=filename,exist=LTEST)
    !// Only if LCONT is true
    if (LCONT .and. LTEST) then
       !// Initialise from file
       !// Data are stored as (IPAR,JPAR), must be converted to
       !// (IDBLK,JDLK,MPBLK)
       open(fnr,file=filename,form='unformatted')
       read(fnr) dday
       read(fnr) INFIELD3
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = 1, ILMM
                   BSNOWL(L,II,JJ,MP) = INFIELD3(L,I,J)
                end do
             end do
          end do
       end do

       read(fnr) INFIELD3
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = 1, ILMM
                   BBCFFC(L,II,JJ,MP) = INFIELD3(L,I,J)
                end do
             end do
          end do
       end do

       read(fnr) INFIELD3
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = 1, ILMM
                   BBCBFC(L,II,JJ,MP) = INFIELD3(L,I,J)
                end do
             end do
          end do
       end do

       read(fnr) INFIELD3
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = 1, ILMM
                   BBCBIO(L,II,JJ,MP) = INFIELD3(L,I,J)
                end do
             end do
          end do
       end do

       read(fnr) INFIELD2
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                tSinceLastSF(II,JJ,MP) = INFIELD2(I,J)
             end do
          end do
       end do

       read(fnr) INSNWIJ
       do MP = 1, MPBLK
          !// Loop over latitude (J is global, JJ is block)
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             !// Loop over longitude (I is global, II is block)
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                LSNW_IJ(II,JJ,MP) = INSNWIJ(I,J)
             end do
          end do
       end do

       !// Close file
       close(fnr)
       write(*,'(a)') f90file//':'//subr// &
            ': Initialized BCsnow from restart file '//trim(filename)
       write(*,'(a,i4)') '  Day on file: ',dday
       write(*,'(a,i4)') '  Model day  : ',startday

    else

       if (LCONT .and. .not.LTEST) then
          !// Allow continuation for tracers, but zero for BCsnow.
          do L=1,3
             write(*,'(a)') '*** WARNING BCsnow WARNING BCsnow WARNING ***'
          end do
          write(*,'(a)') '* Starting from restart file, but not for BCsnow!'
          do L=1,3
             write(*,'(a)') '*** WARNING BCsnow WARNING BCsnow WARNING ***'
          end do
       end if

       BSNOWL(:,:,:,:) = 0._r8
       BBCFFC(:,:,:,:) = 0._r8
       BBCBFC(:,:,:,:) = 0._r8
       BBCBIO(:,:,:,:) = 0._r8
       LSNW_IJ(:,:,:) = 0
       !// Seconds since last snowfall, must be > 86400.
       tSinceLastSF(:,:,:) = 1.e5_r8
       !// Set snow depth from metdata
       do MP = 1, MPBLK
          do J = MPBLKJB(MP),MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP),MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                !// Do not allow more than 0.2m water equivalents; this
                !// is to avoid glacier grid boxes with thick ice/snow.
                if (SD(I,J) .gt. 0._r8) then
                   BSNOWL(1,II,JJ,MP) = min(SD(I,J), metdataSDmax)
                   LSNW_IJ(II,JJ,MP) = 1
                   !// Set up possible thinlayer for drydep of BC
                   if (SD(I,J) .gt. thinsnowthreshold) then
                      BSNOWL(1,II,JJ,MP) = min(SD(I,J), metdataSDmax) - thinsnow
                      BSNOWL(2,II,JJ,MP) = thinsnow
                      LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                   end if
                end if
             end do
          end do
       end do
    end if

    !// Initialize variables for collecting BC in drydep and precipitation
    bcsnow_dd_bio(:,:,:) = 0._r8
    bcsnow_dd_ffc(:,:,:) = 0._r8
    bcsnow_dd_bfc(:,:,:) = 0._r8
    bcsnow_prec_bio(:,:,:) = 0._r8
    bcsnow_prec_ffc(:,:,:) = 0._r8
    bcsnow_prec_bfc(:,:,:) = 0._r8


    write(*,'(a,3es12.3,i4)') f90file//':'//subr// &
         ': Init BCsnowFFC/BFC/BIO maxBSNOWL: ',&
         sum(BBCFFC),sum(BBCBFC),sum(BBCBIO),maxval(LSNW_IJ)

    !// --------------------------------------------------------------------
  end subroutine bcsnow_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_diagwetrm(BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// Diagnose wet removal of BC on snow. Uses BTT and BTTBCK, and is
    !// therefore coupled to the CTM diagnostics.
    !// Must therefore be called after each wet removal process in p-main.
    !//
    !// Amund Sovde, March 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LBCOC
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BTT, BTTBCK
    integer, intent(in) :: MP
    !// Locals
    integer :: I,J,II,JJ,SNL,TRNR
    real(r8)  :: total_lost, SLfact
    !// --------------------------------------------------------------------

    !// Check for BCOC and BCsnow
    if (.not. (LBCOC .and. LBCsnow)) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1
       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1


          SNL = LSNW_IJ(II,JJ,MP)
          !// If no snow layers, go to next I/J
          if (SNL .eq. 0) cycle


          !// Biomass BC hydrophobic
          TRNR = bc_trsp_idx(1)
          !// Total (positive) amount lost due to this process
          total_lost = sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Biomass BC hydrophilic
          TRNR = bc_trsp_idx(2)
          !// Total (positive) amount lost due to this process
          total_lost = total_lost &
               + sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Limit very small values
          if (abs(total_lost) .lt. 1.e-10_r8) then
             total_lost = 0._r8
          end if

          !// Collect BC from precipitation
          bcsnow_prec_bio(II,JJ,MP) = bcsnow_prec_bio(II,JJ,MP) + total_lost


          !// Fossil fuel BC hydrophobic
          TRNR = bc_trsp_idx(3)
          !// Total (positive) amount lost due to this process
          total_lost = sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Fossil fuel BC hydrophilic
          TRNR = bc_trsp_idx(4)
          !// Total (positive) amount lost due to this process
          total_lost = total_lost &
               + sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Limit very small values
          if (abs(total_lost) .lt. 1.e-10_r8) then
             total_lost = 0._r8
          end if

          !// Collect BC FF from precipitation
          bcsnow_prec_ffc(II,JJ,MP) = bcsnow_prec_ffc(II,JJ,MP) + total_lost


          !// Bio fuel BC hydrophobic
          TRNR = bc_trsp_idx(5)
          !// Total (positive) amount lost due to this process
          total_lost = sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Bio fuel BC hydrophilic
          TRNR = bc_trsp_idx(6)
          !// Total (positive) amount lost due to this process
          total_lost = total_lost &
               + sum(BTTBCK(:,TRNR,II,JJ) - BTT(:,TRNR,II,JJ))

          !// Limit very small values
          if (abs(total_lost) .lt. 1.e-10_r8) then
             total_lost = 0._r8
          end if

          !// Collect BC BF from precipitation
          bcsnow_prec_bfc(II,JJ,MP) = bcsnow_prec_bfc(II,JJ,MP) + total_lost

       end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine bcsnow_diagwetrm
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_save_restart(kday)
    !// --------------------------------------------------------------------
    !// Save restart-file for BCsnow.
    !//
    !// Amund Sovde, April 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_met, only: MYEAR
    use cmn_oslo, only: trsp_idx
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: kday
    !// Locals
    logical :: LTEST
    integer :: fnr, I, J, II, JJ, MP, L
    real(r8) :: OUTFIELD3(ILMM,IPAR,JPAR), OUTFIELD2(IPAR,JPAR)
    !logical :: OUTSFL8(8,IPAR,JPAR)
    integer :: OUTSNWIJ(IPAR,JPAR)
    character(len=80) :: filename
    character(len=5) :: cday
    character(len=4) :: cyear
    !// --------------------------------------------------------------------

    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return


    !// Continue from previous run?
    fnr = get_free_fileid()

    !// Filename
    write(cday(1:5),'(i5.5)') KDAY
    write(cyear(1:4),'(i4.4)') MYEAR
    filename = 'restart_bcsnow_'//cday//'_'//cyear//'.sav'


    !// Initialise file
    !// Data are stored as (IPAR,JPAR), must be converted from
    !// (IDBLK,JDLK,MPBLK)
    open(fnr,file=filename,form='unformatted')
    write(fnr) kday
    !// BSNOWL
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, ILMM
                OUTFIELD3(L,I,J) = BSNOWL(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCFFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, ILMM
                OUTFIELD3(L,I,J) = BBCFFC(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCBFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, ILMM
                OUTFIELD3(L,I,J) = BBCBFC(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCBIO
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, ILMM
                OUTFIELD3(L,I,J) = BBCBIO(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// tSinceLastSF
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTFIELD2(I,J) = tSinceLastSF(II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTFIELD2

    !// LSNW_IJ
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTSNWIJ(I,J) = LSNW_IJ(II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTSNWIJ

    close(fnr)
    !// --------------------------------------------------------------------
  end subroutine bcsnow_save_restart
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_status(NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// Print status of BCsnow.
    !//
    !// Amund Sovde, April 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: STT, MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_oslo, only: trsp_idx
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    !// Locals
    integer :: I, J, II, JJ, MP, L, fnr
    real(r8) :: totBCFFil, totBCBFil, totBCBBil
    real(r8) :: totBCFFob, totBCBFob, totBCBBob
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_status'
    !//---------------------------------------------------------------------

    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return

    !// In snow
    totBCBBil = sum(BBCBIO)
    totBCFFil = sum(BBCFFC)
    totBCBFil = sum(BBCBFC)
    write(6,'(a,3es11.3)') f90file//':'//subr// &
         ': Tot BCsnowFFC/BFC/BIO: ',totBCFFil,totBCBFil,totBCBBil

    !// In atmosphere
    totBCBBob = sum(STT(:,:,:,bc_trsp_idx(1)))
    totBCBBil = sum(STT(:,:,:,bc_trsp_idx(2)))
    totBCFFob = sum(STT(:,:,:,bc_trsp_idx(3)))
    totBCFFil = sum(STT(:,:,:,bc_trsp_idx(4)))
    totBCBFob = sum(STT(:,:,:,bc_trsp_idx(5)))
    totBCBFil = sum(STT(:,:,:,bc_trsp_idx(6)))
    write(6,'(a,3es10.2)') '  Tot BC FF/BF/BB (fil): ', &
         totBCFFil,totBCBFil,totBCBBil
    write(6,'(a,3es10.2)') '  Tot BC FF/BF/BB (fob): ', &
         totBCFFob,totBCBFob,totBCBBob

    call bcsnow_check_snow(NDAY,NMET,NOPS)

    !// --------------------------------------------------------------------
  end subroutine bcsnow_status
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_check_snow(NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// Check snow layers.
    !//
    !// Amund Sovde, April 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    !// Locals
    integer :: I, J, II, JJ, MP, L
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_check_snow'
    !//---------------------------------------------------------------------

    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return

    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             if (sum(BSNOWL(:,II,JJ,MP)) .gt. 0._r8 .and. &
                  LSNW_IJ(II,JJ,MP) .eq. 0) then
                write(6,'(a)') f90file//':'//subr// &
                     ': No snow layer, but snow exist!'
                write(6,'(a,5i7)') '  i,j,NDAY,NMET,NOPS:', &
                     i,j,NDAY,NMET,NOPS
                do L = ILMM,1,-1
                   write(6,'(i3,es20.12)') L, BSNOWL(L,II,JJ,MP)
                end do
                stop
             end if
             if (sum(BSNOWL(:,II,JJ,MP)) .gt. 0._r8 .and. &
                  LSNW_IJ(II,JJ,MP) .gt. 0) then
                do L = 1, LSNW_IJ(II,JJ,MP)
                   if (BSNOWL(L,II,JJ,MP) .eq. 0._r8) then
                      write(6,'(a,i3,es20.12)') f90file//':'//subr// &
                           ': No snow in layer ', L, BSNOWL(L,II,JJ,MP)
                      write(6,'(a,5i7)') '  i,j,NDAY,NMET,NOPS:', &
                           i,j,NDAY,NMET,NOPS
                      stop
                   end if
                end do
             end if
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine bcsnow_check_snow
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine bcsnow_getspringsummer(YDGRD, JDAY, &
       springday,springend,summermid,summerend,KDAY)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), intent(in)   :: ydgrd
    integer, intent(in)  :: JDAY
    integer, intent(out) :: KDAY,springday,springend,summermid,summerend
    !// --------------------------------------------------------------------

    !// Melting start days NH/SH
    !// NH: 136 16. May, SH: 320 16. November
    if (YDGRD .ge. 0._r8) then
       !// NH
       springday = 136
       springend = 172
       summermid = 219
       summerend = 244
       KDAY = JDAY
    else
       !// SH
       springday = 320
       springend = 354
       !// Adjust summermid and KDAY so days are monotonic
       summermid =  36 + 365
       summerend =  60 + 365
       if (JDAY .le. 60) then
          KDAY = JDAY + 365
       else
          KDAY = JDAY
       end if
    end if

    !// --------------------------------------------------------------------
  end subroutine bcsnow_getspringsummer
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_nmet_output(NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// Put BCsnow to file each NMET.
    !//
    !// Amund Sovde, April 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: STT, XDGRD, YDGRD, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_met, only: MYEAR
    use cmn_oslo, only: trsp_idx
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    !// Locals
    integer :: I, J, II, JJ, MP, fnr
    character(len=80) :: filename
    character(len=5) :: cday
    character(len=4) :: cyear
    real(r8) :: OUTFIELD3(ILMM,IPAR,JPAR)
    integer :: OUTSNWIJ(IPAR,JPAR)
    !// --------------------------------------------------------------------
    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return

    !// Only put out data each NMET
    if (NOPS.ne.1) return

    !// Filename
    write(cday(1:5),'(i5.5)') NDAY
    write(cyear(1:4),'(i4.4)') MYEAR
    filename = 'bcsnow_nmet_'//cday//'_'//cyear//'.dta'
    fnr = get_free_fileid()


    !// New file?
    if (NMET .eq. 1) then
       open(fnr,file=filename,form='unformatted')
       write(fnr) IPAR,JPAR,ILMM
       write(fnr) XDGRD,YDGRD
       close(fnr)
    end if

    !// IJ-bloc 2 IPARxJPAR
    open(fnr,file=filename,form='unformatted',position='append')
    write(fnr) NMET
    !// BSNOWL
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTFIELD3(:,I,J) = BSNOWL(:,II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCFFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTFIELD3(:,I,J) = BBCFFC(:,II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCBFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTFIELD3(:,I,J) = BBCBFC(:,II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// BBCBIO
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTFIELD3(:,I,J) = BBCBIO(:,II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTFIELD3

    !// LSNW_IJ
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTSNWIJ(I,J) = LSNW_IJ(II,JJ,MP)
          end do
       end do
    end do
    write(fnr) OUTSNWIJ
    !// Close file
    close(fnr)


    !// Print BCsnow status every 3 hours
    if (NMET .eq. 1) call bcsnow_status(NDAY,NMET,NOPS)

    !// --------------------------------------------------------------------
  end subroutine bcsnow_nmet_output
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_nmet_output_nc(NDAY,NMET,NOPS)
    !// --------------------------------------------------------------------
    !// Put BCsnow to file each NMET, but write to netcdf file.
    !//
    !// Amund Sovde, October 2014
    !// --------------------------------------------------------------------
    use netcdf
    use cmn_size, only: MPBLK
    use cmn_ctm, only: STT, XDGRD, YDGRD, AREAXY, &
         MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_diag, only: nc4deflate_global, nc4shuffle_global
    use cmn_met, only: MYEAR
    use cmn_oslo, only: trsp_idx
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS
    !// Locals
    integer :: I, J, II, JJ, MP, L
    character(len=80) :: filename
    character(len=3) :: cday
    character(len=4) :: cyear

    !// Dimensions ID and counters
    integer :: lat_dim_id              !Dimension ID for latitude
    integer :: lon_dim_id              !Dimension ID for longitude
    integer :: lev_dim_id              !Dimension ID for niv
    integer :: time_dim_id             !Dimension ID for time
    integer :: time_id                 !Variable ID for time
    integer :: lon_id                  !Variable ID for longitude
    integer :: lat_id                  !Variable ID for latitude
    integer :: lev_id                  !Variable ID for latitude
    integer :: niv(ILMM)
    integer :: dim_lon_lat_lev_time(4) !Dimension ID for processes
    integer :: srt_lon_lat_lev_time(4) !Start array for lon/lat/time
    integer :: cnt_lon_lat_lev_time(4) !Counting array for lon/lat/time
    integer :: dim_lon_lat_time(3)     !Dimension ID for processes
    integer :: srt_lon_lat_time(3)     !Start array for lon/lat/time
    integer :: cnt_lon_lat_time(3)     !Counting array for lon/lat/time
    integer :: srt_time(1)             !starting point for time array
    integer :: f1_label, f2_label, f2b_label, f3_label, f4_label, f5_label


    !// Other locals
    integer :: ncid                   !File ID for nc file
    integer :: status                 !Error status for nc file
    integer :: nlons                  !Number of longitudes found in file
    integer :: nlats                  !Number of latitudes found in file
    integer :: nlevs                  !Number of levels from nc-file
    integer :: nsteps                 !Number of steps found in file 
    character(len=80) :: time_label   !Label for variable "time"
    real(r8)  :: time                 !Time in this timestep
  
    real(r8) :: OUTFIELD3(IPAR,JPAR,ILMM)
    integer :: OUTSNWIJ(IPAR,JPAR)
    real(r8) :: xdgrd_loc(IPAR)
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = nc4deflate_global
    integer, parameter :: nc4shuffle = nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_nmet_output_nc'
    !//---------------------------------------------------------------------
    !// No need to initialize if BCsnow is not included
    if (.not.LBCsnow) return

    !// Note that CTM3 starts 3-hour output on 00UTC, not 03UTC as was done
    !// in CTM2. This is consistent with the of_gmdump3hrs.f90.
    if (NOPS.ne.1) return

    !// Filename
    write(cday(1:3),'(i3.3)') NDAY
    write(cyear(1:4),'(i4.4)') MYEAR
    filename = 'bcsnow_'//cyear//'_'//cday//'.nc'


    !write(6,*)'will put BCsnow to nc-file '//trim(filename)

    !// longitudes 0:360
    do i = 1, ipar
       if (XDGRD(i) .lt. 0._r8) then
          xdgrd_loc(i) = XDGRD(i) + 360._r8
       else
          xdgrd_loc(i) = XDGRD(i)
       end if
    end do

    !// No errors so far
    status=nf90_noerr 
       
    if (NMET .eq. 1) then
       !// New file; create file
       !// Clobber means you can overwrite existing data:
       status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
       if (status/=nf90_noerr) call handle_error(status,'in creating file')

       !//File headers
       status=nf90_put_att(ncid,nf90_global,'title','3h BC on snow output fields from Oslo CTM3')
       if (status/=nf90_noerr) call handle_error(status,'file header')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo1','Oslo CTM3 is a 3D Chemical Transport Model')
       if (status/=nf90_noerr) call handle_error(status,'modelifo1')
       status=nf90_put_att(ncid,nf90_global,'Modelinfo2','Oslo CTM3 is driven by ECMWF meteorological data')
       if (status/=nf90_noerr) call handle_error(status,'modelinfo2')
       status=nf90_put_att(ncid,nf90_global,'contactinfo','For errors, contact asovde@cicero.oslo.no')
       if (status/=nf90_noerr) call handle_error(status,'contactinfo')

       !// Define dimensions (JM, IM, LM, time)
       status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define lat dim')
       status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define lon dim')
       status = nf90_def_dim(ncid,"niv",ILMM,lev_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define niv dim')
       !// Defining dimension time of length unlimited
       status = nf90_def_dim(ncid,"time",nf90_unlimited,time_dim_id)
       if (status/=nf90_noerr) call handle_error(status,'define time dim')
     
       !// Defining the combined id for a field (lon / lev /lat /time)
       dim_lon_lat_lev_time(1)=lon_dim_id
       dim_lon_lat_lev_time(2)=lat_dim_id
       dim_lon_lat_lev_time(3)=lev_dim_id
       dim_lon_lat_lev_time(4)=time_dim_id
       dim_lon_lat_time(1)=lon_dim_id
       dim_lon_lat_time(2)=lat_dim_id
       dim_lon_lat_time(3)=time_dim_id
     
       !// Defining the lon/lat/lev-variable
       status = nf90_def_var(ncid,"lon",nf90_float,lon_dim_id,lon_id)
       if (status/=nf90_noerr) call handle_error(status,'define lon variable')
       status = nf90_def_var(ncid,"lat",nf90_float,lat_dim_id,lat_id)
       if (status/=nf90_noerr) call handle_error(status,'define lat variable')
       status = nf90_def_var(ncid,"niv",nf90_float,lev_dim_id,lev_id)
       if (status/=nf90_noerr) call handle_error(status,'define lev variable')
       status = nf90_def_var(ncid,"time",nf90_float,time_dim_id,time_id)
       if (status/=nf90_noerr) call handle_error(status,'define time variable')

       !// Putting attributes to /lon/lat/lev variables
       status = nf90_put_att(ncid,lon_id,'units','degree_east')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lon')
       status = nf90_put_att(ncid,lat_id,'units','degree_north')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lat')
       status = nf90_put_att(ncid,lev_id,'units','levels')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lev')

       time_label='hours since 00:00:00'

       status = nf90_put_att(ncid,time_id,'units',time_label)
       if (status/=nf90_noerr) call handle_error(status,'attribute units time')
       
       !// Define the tracer field variables

       !// Snow layer thickness
       status = nf90_def_var(ncid,'snowlayr',nf90_float, dim_lon_lat_lev_time, f1_label)
       if (status/=nf90_noerr) call handle_error(status,'define snowlayr')
       status = nf90_def_var_deflate(ncid,f1_label,nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate snowlayr variable')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,f1_label,'units','m water equiv.')
       if (status/=nf90_noerr) call handle_error(status,'attribute units snowlayr')

       !// bcffcsnw
       status = nf90_def_var(ncid,'bcffcsnw',nf90_float, dim_lon_lat_lev_time, f2_label)
       if (status/=nf90_noerr) call handle_error(status,'define bcffcsnw')
       status = nf90_def_var_deflate(ncid,f2_label,nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate bcffcsnw variable')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,f2_label,'units','kg/m2')
       if (status/=nf90_noerr) call handle_error(status,'attribute units bcffcsnw')

       !// bcbfcsnw
       status = nf90_def_var(ncid,'bcbfcsnw',nf90_float, dim_lon_lat_lev_time, f2b_label)
       if (status/=nf90_noerr) call handle_error(status,'define bcbfcsnw')
       status = nf90_def_var_deflate(ncid,f2b_label,nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate bcbfcsnw variable')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,f2b_label,'units','kg/m2')
       if (status/=nf90_noerr) call handle_error(status,'attribute units bcbfcsnw')

       !// bcbiosnw
       status = nf90_def_var(ncid,'bcbiosnw',nf90_float, dim_lon_lat_lev_time, f3_label)
       if (status/=nf90_noerr) call handle_error(status,'define bcbiosnw')
       status = nf90_def_var_deflate(ncid,f3_label,nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate bcbiosnw variable')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,f3_label,'units','kg/m2')
       if (status/=nf90_noerr) call handle_error(status,'attribute units bcbiosnw')

       !// lsnw
       status = nf90_def_var(ncid,'lsnw',nf90_int, dim_lon_lat_time, f4_label)
       if (status/=nf90_noerr) call handle_error(status,'define lsnw')
       status = nf90_def_var_deflate(ncid,f4_label,nc4shuffle,1,nc4deflate)
       if (status/=nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate lsnw variable')
       !// Add text descriptions and units
       status = nf90_put_att(ncid,f4_label,'units','levels')
       if (status/=nf90_noerr) call handle_error(status,'attribute units lsnw')


       !// End definition mode
       status = nf90_enddef(ncid)
       if (status/=nf90_noerr) call handle_error(status,'end defmode')
     
       !// Putting the lon/lat/lev variables
       !write(6,*)'Putting variables for lon/lev/lat. They do not change'
       status = nf90_put_var(ncid,lon_id,xdgrd_loc)
       if (status/=nf90_noerr) call handle_error(status,'putting lon')
       status = nf90_put_var(ncid,lat_id,ydgrd)
       if (status/=nf90_noerr) call handle_error(status,'putting lat')
       !// Levels array
       do L = 1, ILMM
          niv(L) = L
       end do
       status = nf90_put_var(ncid,lev_id,niv)
       if (status/=nf90_noerr) call handle_error(status,'putting niv')

    else
       !Open the existing file
       status=nf90_open(filename, nf90_write, ncid)
       if(status/=nf90_noerr)call handle_error(status,'opening file '//trim(filename))

       !Inquire dimension ids
       status = nf90_inq_dimid(ncid,"lat",lat_dim_id)
       if(status/=nf90_noerr)call handle_error(status,'retrieveing id lat')
       status = nf90_inq_dimid(ncid,"lon",lon_dim_id)
       if(status/=nf90_noerr)call handle_error(status,'retrieveing id lon')
       status = nf90_inq_dimid(ncid,"niv",lev_dim_id)
       if(status/=nf90_noerr)call handle_error(status,'retrieveing id lev')
       status = nf90_inq_dimid(ncid,"time",time_dim_id)
       if(status/=nf90_noerr)call handle_error(status,'retrieveing id time')
 
       !Inquire dimensions
       status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
       if(status/=nf90_noerr)call handle_error(status,'retrieving lats')
       if(nlats/=JPAR)then
          write(6,*)'File reports JPAR = ',nlats
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
       if(status/=nf90_noerr)call handle_error(status,'retrieving lons')
       if(nlons/=IPAR)then
          write(6,*)'File reports IPAR = ',nlons
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,lev_dim_id,len=nlevs)
       if(status/=nf90_noerr)call handle_error(status,'retrieving levs')
       if(nlevs/=ILMM)then
          write(6,*)'File reports ILMM = ',nlevs
          stop
       end if
       status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
       if(status/=nf90_noerr)call handle_error(status,'retrieving nsteps')
       if(nsteps .lt. 0)then
          write(6,*)'File reports already added nsteps = ',nsteps
          stop
       end if

       !Check variable ids
       status=nf90_inq_varid(ncid,'snowlayr',f1_label)
       if(status/=nf90_noerr)call handle_error(status,'retrieving f1_label')
       status=nf90_inq_varid(ncid,'bcffcsnw',f2_label)
       if(status/=nf90_noerr)call handle_error(status,'retrieving f2_label')
       status=nf90_inq_varid(ncid,'bcbfcsnw',f2b_label)
       if(status/=nf90_noerr)call handle_error(status,'retrieving f2b_label')
       status=nf90_inq_varid(ncid,'bcbiosnw',f3_label)
       if(status/=nf90_noerr)call handle_error(status,'retrieving f3_label')
       status=nf90_inq_varid(ncid,'lsnw',f4_label)
       if(status/=nf90_noerr)call handle_error(status,'retrieving f4_label')

       !Get variable id for time
       status=nf90_inq_varid(ncid,"time",time_id) 
       if(status/=nf90_noerr)call handle_error(status,'retrieving time id')

    end if


    !// Now write to file
    cnt_lon_lat_lev_time = (/IPAR, JPAR, ILMM, 1/)
    srt_lon_lat_lev_time = (/1, 1, 1, NMET/)

    cnt_lon_lat_time = (/IPAR, JPAR, 1/)
    srt_lon_lat_time = (/1, 1, NMET/)

    srt_time(1)=NMET             !Start value for new time step
    time=real(NMET-1, r8)*3._r8       !Time in real format



    status=nf90_put_var(ncid,time_id,time,start=srt_time)
    if(status/=nf90_noerr)call handle_error(status,'putting time')

    !// Initialise
    OUTFIELD3(:,:,:) = 0._r8

    !// BSNOWL
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, LSNW_IJ(II,JJ,MP)
                OUTFIELD3(I,J,LSNW_IJ(II,JJ,MP)+1-L) = BSNOWL(L,II,JJ,MP)
             end do
          end do
       end do
    end do
    status=nf90_put_var(ncid, f1_label, outfield3, &
          start=srt_lon_lat_lev_time,      &  !starting point for writing
          count=cnt_lon_lat_lev_time )        !Counts how many bytes written
    if(status/=nf90_noerr)call handle_error(status,'putting BSNOWL')

    !// BBCFFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, LSNW_IJ(II,JJ,MP)
                OUTFIELD3(I,J,LSNW_IJ(II,JJ,MP)+1-L) = &
                     BBCFFC(L,II,JJ,MP) / AREAXY(I,J)
             end do
          end do
       end do
    end do
    status=nf90_put_var(ncid, f2_label, outfield3, &
          start=srt_lon_lat_lev_time,      &  !starting point for writing
          count=cnt_lon_lat_lev_time )        !Counts how many bytes written
    if(status/=nf90_noerr)call handle_error(status,'putting BBCFFC')

    !// BBCBFC
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, LSNW_IJ(II,JJ,MP)
                OUTFIELD3(I,J,LSNW_IJ(II,JJ,MP)+1-L) = &
                     BBCBFC(L,II,JJ,MP) / AREAXY(I,J)
             end do
          end do
       end do
    end do
    status=nf90_put_var(ncid, f2b_label, outfield3, &
          start=srt_lon_lat_lev_time,      &  !starting point for writing
          count=cnt_lon_lat_lev_time )        !Counts how many bytes written
    if(status/=nf90_noerr)call handle_error(status,'putting BBCBFC')

    !// BBCBIO
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             do L = 1, LSNW_IJ(II,JJ,MP)
                OUTFIELD3(I,J,LSNW_IJ(II,JJ,MP)+1-L) = &
                     BBCBIO(L,II,JJ,MP) / AREAXY(I,J)
             end do
          end do
       end do
    end do
    status=nf90_put_var(ncid, f3_label, outfield3, &
          start=srt_lon_lat_lev_time,      &  !starting point for writing
          count=cnt_lon_lat_lev_time )        !Counts how many bytes written
    if(status/=nf90_noerr)call handle_error(status,'putting BBCBIO')

    !// LSNW_IJ
    do MP = 1, MPBLK
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             OUTSNWIJ(I,J) = LSNW_IJ(II,JJ,MP)
          end do
       end do
    end do
    status=nf90_put_var(ncid, f4_label, outsnwij, &
          start=srt_lon_lat_time,      &  !starting point for writing
          count=cnt_lon_lat_time )        !Counts how many bytes written
    if(status/=nf90_noerr)call handle_error(status,'putting LSNW_IJ')

    !close netcdf file
    status=nf90_close(ncid)
    if(status/=nf90_noerr) &
         call handle_error(status,'closing file'//trim(filename))


    !// Print BCsnow status every 3 hours
    if (NMET .eq. 1) call bcsnow_status(NDAY,NMET,NOPS)

    !// --------------------------------------------------------------------
  end subroutine bcsnow_nmet_output_nc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_collect_ij(DT,MP)
    !// --------------------------------------------------------------------
    !// Here we collect the deposited and precipitated BC and build
    !// snow layers from meteorological data.
    !//
    !// Originally written by Ragnhild B. Skeie for CTM2, re-structured to
    !// fit CTM3 in March-April 2013.
    !// Routine is called each operator split time step, whereas in
    !// CTM2 it was done every meteorological time step.
    !//
    !// Renamings
    !//   SL     => BSNOWL (snow level depth [m water equivalents])
    !//   SBCFFC => BBCFFC (BCFFC in each snow level [kg])
    !//   SBCBFC => BBCBFC (BCBFC in each snow level [kg])
    !//   SBCBIO => BBCBIO (BCBIO in each snow level [kg])
    !//
    !// Amund Sovde, March-April 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: NRMETD, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: SNFL, PRECCNV, PRECLS
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: DT
    integer, intent(in) :: MP

    !// Locals
    integer :: I, J, II, JJ, K, L, SNL
    real(r8)  :: snowfall, thinfrac
    logical :: oldL
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_collect_ij'
    !//---------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Snow fall (SF) at surface, should have units m/s
        if ( (PRECCNV(I,J,1) + PRECLS(I,J,1)) .eq. 0._r8 ) then
           !// No snow precipitation
           snowfall = 0._r8
        else if ((PRECCNV(I,J,1) + PRECLS(I,J,1)) .lt. 0._r8) then
           write(6,'(a,2i7,2es16.6)') f90file//':'//subr// &
                ': Negative precipitation!!',i,j,PRECCNV(I,J,1),PRECLS(I,J,1)
           stop
        else
           !// There is precipitation; get snowfall
           snowfall = SNFL(II,JJ,MP) * DT
        end if


        !// If there is no snowfall and no existing snowlayer, there is
        !// nothing to be done.
        if (snowfall .le. SFLIM) then
           !// Count seconds since last snowfall (i.e. after this DT is done).
           !// In case of snowfall, tSinceLastSF is set at the bottom of
           !// the loop.
           tSinceLastSF(II,JJ,MP) = tSinceLastSF(II,JJ,MP) + DT
           !// There is no snowfall, but it may still rain and we can
           !// have wet deposition on snow.
           !//
           !// If there is no snowfall, and no snowlayer exist, we
           !// move to next (II,JJ,MP).
           if (LSNW_IJ(II,JJ,MP) .eq. 0) then
              !// There is no snow layer, so BC cannot dry deposit on snow.
              !// Reset collected variables for next round:
              bcsnow_dd_bio(II,JJ,MP) = 0._r8
              bcsnow_dd_ffc(II,JJ,MP) = 0._r8
              bcsnow_dd_bfc(II,JJ,MP) = 0._r8
              !// And BC zero wet deposition.
              !// Reset collected variables for next round:
              bcsnow_prec_bio(II,JJ,MP) = 0._r8
              bcsnow_prec_ffc(II,JJ,MP) = 0._r8
              bcsnow_prec_bfc(II,JJ,MP) = 0._r8
              !// If there is no snow, go to next (II,JJ,MP)
              cycle
           end if
        end if !// if (snowfall .le. SFLIM) then



        !// Being here, we have 3 options:
        !// 1. no snowfall, but snow layer exist
        !//    - account for drydep on snow only
        !//    - put it in the top thinlayer
        !// 2. snowfall, and snow layer exist
        !//    - if precip more recently than 24 hours, collapse drydep
        !//      thinlayer with snow layer below. Otherwise, keep drydep
        !//      thinlayer and start new snow layer on top.
        !//    - account for drydep+precip in snow
        !// 3. snowfall, but no snow layer exist
        !//    - if precip ... same as for (2.)
        !//    - account for drydep+precip in snow


        !// Dry deposition is collected in the top thin snow layer.
        !// Snowfall adds BC to the top layer and at the end a new thin
        !// snow layer is generated at the top.

        if (snowfall .le. SFLIM) then
           !// Only DRYDEP is to be accounted for!
           !// -------------------------------------------------------------
           !// Add drydepped BC to thinlayer
           SNL = LSNW_IJ(II,JJ,MP) !// Top layer

           !// DEBUG
           if (LDEBUG_BCsnow) then
              if (SNL .eq. 0) then
                 write(6,'(a)') f90file//':'//subr// &
                      ': How did we get here? No snowlayer!'
                 stop
              end if
           end if

           !// Collect BC from drydep.
           !// Also collect precipitiation, since it might be
           !// raining on a snowy ground.
           BBCFFC(SNL,II,JJ,MP) = BBCFFC(SNL,II,JJ,MP) &
                                  + bcsnow_dd_ffc(II,JJ,MP) &
                                  + bcsnow_prec_ffc(II,JJ,MP)
           BBCBFC(SNL,II,JJ,MP) = BBCBFC(SNL,II,JJ,MP) &
                                  + bcsnow_dd_bfc(II,JJ,MP) &
                                  + bcsnow_prec_bfc(II,JJ,MP)
           BBCBIO(SNL,II,JJ,MP) = BBCBIO(SNL,II,JJ,MP) &
                                  + bcsnow_dd_bio(II,JJ,MP) &
                                  + bcsnow_prec_bio(II,JJ,MP)

           !// Reset collected variables for next round
           bcsnow_dd_bio(II,JJ,MP) = 0._r8
           bcsnow_dd_ffc(II,JJ,MP) = 0._r8
           bcsnow_dd_bfc(II,JJ,MP) = 0._r8
           !// Skip precipitation, since there is no snow
           bcsnow_prec_bio(II,JJ,MP) = 0._r8
           bcsnow_prec_ffc(II,JJ,MP) = 0._r8
           bcsnow_prec_bfc(II,JJ,MP) = 0._r8
           !// Go to next (II,JJ,MP)
           cycle

        else
           !// DRYDEP and PRECIPITATION is to be accounted for!
           !// -------------------------------------------------------------

           !// Check for existing snow layers
           !// 1. There is no snow: start generating layers.
           !// 2. Thers is snow: Add to it or start new layer.
           if (LSNW_IJ(II,JJ,MP) .eq. 0) then
              !// There is no snow yet: Start a new layer
              !// ----------------------------------------------------------

              !// Have already tested for snowfall > SFLIM,
              !// i.e., there is significant snowfall for a first layer
              LSNW_IJ(II,JJ,MP) = 1
              SNL = LSNW_IJ(II,JJ,MP)

              !// Snow layer thickness
              BSNOWL(SNL,II,JJ,MP) = snowfall

              !// Collect BC from both drydep and precipitation
              BBCFFC(SNL,II,JJ,MP) = bcsnow_dd_ffc(II,JJ,MP) &
                                     + bcsnow_prec_ffc(II,JJ,MP)
              BBCBFC(SNL,II,JJ,MP) = bcsnow_dd_bfc(II,JJ,MP) &
                                     + bcsnow_prec_bfc(II,JJ,MP)
              BBCBIO(SNL,II,JJ,MP) = bcsnow_dd_bio(II,JJ,MP) &
                                     + bcsnow_prec_bio(II,JJ,MP)

              !// Reset collected variables for next round
              bcsnow_dd_bio(II,JJ,MP) = 0._r8
              bcsnow_dd_ffc(II,JJ,MP) = 0._r8
              bcsnow_dd_bfc(II,JJ,MP) = 0._r8
              bcsnow_prec_bio(II,JJ,MP) = 0._r8
              bcsnow_prec_ffc(II,JJ,MP) = 0._r8
              bcsnow_prec_bfc(II,JJ,MP) = 0._r8


              !// Now, add possible thinlayer to account for drydep:
              !// After this snowfall has created a snow layer, BC may deposit
              !// on snow. The BC is only deposited on top of snow, so we
              !// have to generate a thin snow layer on top.
              !// The thin layer thickness is defined by the variable thinsnow,
              !// and if the newly formed layer is thicker than this limit,
              !// we generate an extra thin layer.
              !// Important:
              !// To avoid tiny tiny layers below, we do not do this until
              !// the layer thickness reaches thinsnowthreshold (typically
              !// 10% higher than thinsnow):
              if (snowfall .gt. thinsnowthreshold) then
                 !// Remove thinsnow from layer to make layer 2 (for deposition)
                 BSNOWL(SNL,II,JJ,MP) = BSNOWL(SNL,II,JJ,MP) - thinsnow
                 !// Make thin layer nr 2
                 BSNOWL(SNL+1,II,JJ,MP) = thinsnow
                 !// Increase snow layer
                 LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                 SNL = LSNW_IJ(II,JJ,MP)
              end if !// if (snowfall .gt. thinsnowthreshold) then

           else
              !// One or more snow layers already exist
              !// ---------------------------------------------------------

              !// Get current snow level
              SNL = LSNW_IJ(II,JJ,MP)

              !// If no snowfall during the last 24 hours, we start
              !// a new snow layer, otherwise we keep the old layer
              oldL = tSinceLastSF(II,JJ,MP) .le. 86400._r8
              !// oldL is .true. if there was snowfall during the last 24 hours


              !// Add to existing layer or start a new one
              if (oldL) then
                 !// Add to existing layer
                 !// -------------------------------------------------------

                 !// When adding to the existing layer, we first incorporate
                 !// the drydep layer to snow layer below it. The drydep
                 !// layer is only assumed to be separate if there is
                 !// no snow for 24 hours.
                 if (SNL .gt. 1) then
                    BSNOWL(SNL-1,II,JJ,MP) = BSNOWL(SNL-1,II,JJ,MP) &
                                             + BSNOWL(SNL,II,JJ,MP)
                    BBCFFC(SNL-1,II,JJ,MP) = BBCFFC(SNL-1,II,JJ,MP) &
                                             + BBCFFC(SNL,II,JJ,MP)
                    BBCBFC(SNL-1,II,JJ,MP) = BBCBFC(SNL-1,II,JJ,MP) &
                                             + BBCBFC(SNL,II,JJ,MP)
                    BBCBIO(SNL-1,II,JJ,MP) = BBCBIO(SNL-1,II,JJ,MP) &
                                             + BBCBIO(SNL,II,JJ,MP)
                    !// Zero current layer and remove it
                    BSNOWL(SNL,II,JJ,MP) = 0._r8
                    BBCFFC(SNL,II,JJ,MP) = 0._r8
                    BBCBFC(SNL,II,JJ,MP) = 0._r8
                    BBCBIO(SNL,II,JJ,MP) = 0._r8

                    !// Remove the zero-layer SNL; it will be re-built below
                    !// if necessary.
                    LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1
                    SNL = LSNW_IJ(II,JJ,MP)

                 end if

                 !// Update snow depth of current layer
                 BSNOWL(SNL,II,JJ,MP) = BSNOWL(SNL,II,JJ,MP) &
                                        + snowfall

                 !// Collect BC from both drydep and precipitation
                 BBCFFC(SNL,II,JJ,MP) = BBCFFC(SNL,II,JJ,MP) &
                                        + bcsnow_dd_ffc(II,JJ,MP) &
                                        + bcsnow_prec_ffc(II,JJ,MP)
                 BBCBFC(SNL,II,JJ,MP) = BBCBFC(SNL,II,JJ,MP) &
                                        + bcsnow_dd_bfc(II,JJ,MP) &
                                        + bcsnow_prec_bfc(II,JJ,MP)
                 BBCBIO(SNL,II,JJ,MP) = BBCBIO(SNL,II,JJ,MP) &
                                        + bcsnow_dd_bio(II,JJ,MP) &
                                        + bcsnow_prec_bio(II,JJ,MP)

                 !// Reset collected variables for next round
                 bcsnow_dd_bio(II,JJ,MP) = 0._r8
                 bcsnow_dd_ffc(II,JJ,MP) = 0._r8
                 bcsnow_dd_bfc(II,JJ,MP) = 0._r8
                 bcsnow_prec_bio(II,JJ,MP) = 0._r8
                 bcsnow_prec_ffc(II,JJ,MP) = 0._r8
                 bcsnow_prec_bfc(II,JJ,MP) = 0._r8

                 !// Extract a thinlayer for drydeposition next round.
                 if (BSNOWL(SNL,II,JJ,MP) .gt. thinsnowthreshold) then

                    !// Make thinlayer for deposition  by removing part of
                    !// existing snow layer.
                    LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                    SNL = LSNW_IJ(II,JJ,MP)

                    !// Find fraction of thinsnow layer to be created:
                    thinfrac = thinsnow/BSNOWL(SNL-1,II,JJ,MP)
                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (thinfrac .gt. 1._r8) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': thinfrac wrong 1',thinfrac
                          stop
                       end if
                    end if

                    !// Split snow layer
                    BSNOWL(SNL,II,JJ,MP) = thinsnow
                    BSNOWL(SNL-1,II,JJ,MP) = BSNOWL(SNL-1,II,JJ,MP) - thinsnow
                    !// And BC according to thinfrac
                    BBCFFC(SNL,II,JJ,MP) = thinfrac * BBCFFC(SNL-1,II,JJ,MP)
                    BBCBFC(SNL,II,JJ,MP) = thinfrac * BBCBFC(SNL-1,II,JJ,MP)
                    BBCBIO(SNL,II,JJ,MP) = thinfrac * BBCBIO(SNL-1,II,JJ,MP)
                    !// Subtract thin layer values from the layer below
                    BBCFFC(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCFFC(SNL-1,II,JJ,MP)
                    BBCBFC(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCBFC(SNL-1,II,JJ,MP)
                    BBCBIO(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCBIO(SNL-1,II,JJ,MP)

                 end if !// if (BSNOWL(SNL,II,JJ,MP) .gt. thinsnowthreshold) then

              else
                 !// Not snowfall during last 24 hours, start new layer
                 !// -------------------------------------------------------

                 !// Outline
                 !// 1. If we have reached max number of levels,
                 !//    collapse 2 bottom layers to 1.
                 !// 2. Create new layer on top.

                 !// Limit number of layers to ILMM. If more layers,
                 !// merge the two bottom layers
                 if (LSNW_IJ(II,JJ,MP) .eq. ILMM) then
                   BSNOWL(1,II,JJ,MP) = BSNOWL(1,II,JJ,MP) + BSNOWL(2,II,JJ,MP)
                   BBCFFC(1,II,JJ,MP) = BBCFFC(1,II,JJ,MP) + BBCFFC(2,II,JJ,MP)
                   BBCBFC(1,II,JJ,MP) = BBCBFC(1,II,JJ,MP) + BBCBFC(2,II,JJ,MP)
                   BBCBIO(1,II,JJ,MP) = BBCBIO(1,II,JJ,MP) + BBCBIO(2,II,JJ,MP)
                   !// Shift levels downwards
                   do L = 2, ILMM-1
                      BSNOWL(L,II,JJ,MP) = BSNOWL(L+1,II,JJ,MP)
                      BBCFFC(L,II,JJ,MP) = BBCFFC(L+1,II,JJ,MP)
                      BBCBFC(L,II,JJ,MP) = BBCBFC(L+1,II,JJ,MP)
                      BBCBIO(L,II,JJ,MP) = BBCBIO(L+1,II,JJ,MP)
                   end do
                   !// Zero BCsnow for next treatment
                   BSNOWL(ILMM,II,JJ,MP) = 0._r8
                   BBCFFC(ILMM,II,JJ,MP) = 0._r8
                   BBCBFC(ILMM,II,JJ,MP) = 0._r8
                   BBCBIO(ILMM,II,JJ,MP) = 0._r8
                   !// No need to change LSNW_IJ due to collapsing
                 else
                   !// Start NEW snow layer
                   LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                 end if

                 !// Snow layer
                 SNL = LSNW_IJ(II,JJ,MP)

                 !// Put snow into the new layer
                 BSNOWL(SNL,II,JJ,MP) = snowfall

                 !// Collect BC from both drydep and precipitation
                 BBCFFC(SNL,II,JJ,MP) = bcsnow_dd_ffc(II,JJ,MP) &
                                        + bcsnow_prec_ffc(II,JJ,MP)
                 BBCBFC(SNL,II,JJ,MP) = bcsnow_dd_bfc(II,JJ,MP) &
                                        + bcsnow_prec_bfc(II,JJ,MP)
                 BBCBIO(SNL,II,JJ,MP) = bcsnow_dd_bio(II,JJ,MP) &
                                        + bcsnow_prec_bio(II,JJ,MP)

                 !// Reset collected variables for next round
                 bcsnow_dd_bio(II,JJ,MP) = 0._r8
                 bcsnow_dd_ffc(II,JJ,MP) = 0._r8
                 bcsnow_dd_bfc(II,JJ,MP) = 0._r8
                 bcsnow_prec_bio(II,JJ,MP) = 0._r8
                 bcsnow_prec_ffc(II,JJ,MP) = 0._r8
                 bcsnow_prec_bfc(II,JJ,MP) = 0._r8

                 !// If snowfall created a layer thinner than thinsnow,
                 !// we are done and leave the layer as is for deposition.
                 !// If, however, the layer is thicker, we add a thin layer
                 !// on top (for deposition), possibly by repeating collapsing
                 !// of bottom layers.
                 if (BSNOWL(SNL,II,JJ,MP) .gt. thinsnowthreshold) then
                    !// Snowfall creates a layer thicker than thinsnow,
                    !// so we add a thin layer (for deposition) on top of it.
                    !// Unless, of course, the remaining layer does not get
                    !// too small.
                    LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                    SNL = LSNW_IJ(II,JJ,MP)
                    !// Re-do collapsing of bottom layers if necessary
                    if (SNL .gt. ILMM) then
                       BSNOWL(1,II,JJ,MP) = BSNOWL(1,II,JJ,MP)+BSNOWL(2,II,JJ,MP)
                       BBCFFC(1,II,JJ,MP) = BBCFFC(1,II,JJ,MP)+BBCFFC(2,II,JJ,MP)
                       BBCBFC(1,II,JJ,MP) = BBCBFC(1,II,JJ,MP)+BBCBFC(2,II,JJ,MP)
                       BBCBIO(1,II,JJ,MP) = BBCBIO(1,II,JJ,MP)+BBCBIO(2,II,JJ,MP)
                       !// Shift levels downwards
                       do L = 2, ILMM-1
                          BSNOWL(L,II,JJ,MP) = BSNOWL(L+1,II,JJ,MP)
                          BBCFFC(L,II,JJ,MP) = BBCFFC(L+1,II,JJ,MP)
                          BBCBFC(L,II,JJ,MP) = BBCBFC(L+1,II,JJ,MP)
                          BBCBIO(L,II,JJ,MP) = BBCBIO(L+1,II,JJ,MP)
                       end do
                       !// Zero BCsnow for next treatment
                       BBCFFC(ILMM,II,JJ,MP) = 0._r8
                       BBCBFC(ILMM,II,JJ,MP) = 0._r8
                       BBCBIO(ILMM,II,JJ,MP) = 0._r8
                       LSNW_IJ(II,JJ,MP) = ILMM
                       SNL = LSNW_IJ(II,JJ,MP)
                    end if
                    !// Remove thinsnow from SNL-1 and put it in SNL (on top)
                    !// Find fraction of thinsnow layer created from 
                    !// layer from below
                    thinfrac = thinsnow/BSNOWL(SNL-1,II,JJ,MP)
                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (thinfrac .gt. 1._r8) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': thinfrac wrong 2',thinfrac
                          stop
                       end if
                    end if

                    !// Split snow layer
                    BSNOWL(SNL,II,JJ,MP) = thinsnow
                    BSNOWL(SNL-1,II,JJ,MP) = BSNOWL(SNL-1,II,JJ,MP) - thinsnow
                    !// And BC according to thinfrac
                    BBCFFC(SNL,II,JJ,MP) = thinfrac * BBCFFC(SNL-1,II,JJ,MP)
                    BBCBFC(SNL,II,JJ,MP) = thinfrac * BBCBFC(SNL-1,II,JJ,MP)
                    BBCBIO(SNL,II,JJ,MP) = thinfrac * BBCBIO(SNL-1,II,JJ,MP)
                    !// Subtract thin layer values from the layer below
                    BBCFFC(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCFFC(SNL-1,II,JJ,MP)
                    BBCBFC(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCBFC(SNL-1,II,JJ,MP)
                    BBCBIO(SNL-1,II,JJ,MP) = &
                         (1._r8 - thinfrac) * BBCBIO(SNL-1,II,JJ,MP)


                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (BSNOWL(1,II,JJ,MP) .eq. 0._r8 .and. &
                            LSNW_IJ(II,JJ,MP) .gt. 0) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': Wrong 3: BSNOWL vs LSNW_IJ'
                          do L = ILMM,1,-1
                             write(6,'(a,i5,es16.6)') '  BSNOWL:',L, BSNOWL(L,II,JJ,MP)
                          end do
                          write(6,'(a,i5)') '  LSNW_IJ:',LSNW_IJ(II,JJ,MP)
                          stop
                       end if
                       if (minval(BSNOWL(:,II,JJ,MP)) .lt. 0._r8) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': Wrong 4: snow layers < 0'
                          do L = ILMM,1,-1
                             write(6,'(a,i5,es16.6)') '  BSNOWL:',L, BSNOWL(L,II,JJ,MP)
                          end do
                          stop
                       end if
                    end if

                 end if !// if (BSNOWL(SNL,II,JJ,MP) .gt. thinsnowthreshold) then

              end if !// if (oldL) then

              !// Reset tSinceLastSF
              tSinceLastSF(II,JJ,MP) = 0._r8

           end if !// if (LSNW_IJ(II,JJ,MP) .eq. 0) then

        end if !// if (snowfall .le. SFLIM) then


      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine bcsnow_collect_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_meltevap_ij(DT,MP)
    !// --------------------------------------------------------------------
    !// This routine evaporates snow layers.
    !//
    !// Originally written by Ragnhild B. Skeie for CTM2, modified to
    !// fit CTM3 in March-April 2013.
    !//
    !//
    !// Amund Sovde, March-April 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_met, only: SMLT, ES, SNFL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: DT
    integer, intent(in) :: MP

    !// Locals
    integer :: I, J, L, II, JJ, ILMX, SNL

    real(r8) :: evapmelt, evapfrac, current_evapmelt
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_meltevap_ij'
    !//---------------------------------------------------------------------


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// There will be no evaporation if there is no snow
        if (LSNW_IJ(II,JJ,MP) .eq. 0) cycle

        !// We do not distinguish between melting and evaporation
        !// Units on SMLT and ES are [m/s water equivalents]
        !// Note that ES may be negative, i.e. depositing snow on surface.
        evapmelt = (SMLT(II,JJ,MP) + ES(II,JJ,MP)) * DT

        !// Nothing to evaporate: go to next (II,JJ,MP)
        if (evapmelt .eq. 0._r8) cycle


        !// Basics of evaporation of snow layers
        !// -----------------------------------------------------------
        !// The uppermost layer is a thin layer for snow deposition.
        !// It does not evaporate until all snow is gone; what evaporates
        !// is the layer beneeth it.
        !//
        !// When the layer beneeth is fully evaporated, its BC is merged
        !// with the top layer, and in this way the top layer concentration
        !// of BC increases.
        !// Only when all layers beneeth are evaporated, the remaing thin
        !// layer may be allowed to evaporate.
        !//
        !// To solve this numerically, we loop downwards from the top layer,
        !// checking if the layer betneeth (if there is one) evaporates.

        !// Number of layers
        ILMX = LSNW_IJ(II,JJ,MP)
        !// Loop downwards
        do L = ILMX, 1, -1

           !// Note that in this loop, L is always the top thin layer, and
           !// melting/evaporation occurs in layer L-1.
           !//
           !// The melted/evaporated BC from layer L-1 is added to layer L,
           !// i.e. the top-most thin layer.

           !// If we only have 1 layer, we cannot evaporate layers beneeth
           !// it; hence, we first need to check for 1 layer and whether
           !// that should be evaporated.
           if (L .eq. 1) then
              !// There is no layer beneeth to evaporate;
              !// do evaporation of layer L only.

              !// For consistency, we double check that we actually have
              !// one layer left
              if (LSNW_IJ(II,JJ,MP) .ne. 1) then
                 write(6,'(a)') f90file//':'//subr//': Inconsistent LSNW_IJ != 1'
                 write(6,'(a,6i5,es16.6)') '  I,J,II,JJ,MP,LSNW_IJ,evapmelt: ', &
                      I,J,II,JJ,MP,LSNW_IJ(II,JJ,MP), evapmelt
                 stop
              end if


              !// If not everything melts/evaporates, we keep a thin layer;
              !// otherwise, all snow and BC is gone.
              !// (Note that this part works on layer L, which is 1.)
              if (evapmelt .ge. BSNOWL(L,II,JJ,MP)) then
                 !// All gone, no snow; BC is lost on ground
                 BSNOWL(L,II,JJ,MP) = 0._r8
                 LSNW_IJ(II,JJ,MP) = 0
                 BBCFFC(L,II,JJ,MP) = 0._r8
                 BBCBFC(L,II,JJ,MP) = 0._r8
                 BBCBIO(L,II,JJ,MP) = 0._r8
                 tSinceLastSF(II,JJ,MP) = 1.e5_r8
                 evapmelt = 0._r8
              else if (evapmelt .lt. 0._r8) then
                 !// Could possibly add to snow layer; skip for L=1 (as in CTM2)
                 evapmelt = 0._r8
              else
                 !// Melt part of layer
                 BSNOWL(L,II,JJ,MP) = BSNOWL(L,II,JJ,MP) - evapmelt
                 LSNW_IJ(II,JJ,MP) = 1
                 evapmelt = 0._r8
                 !// BC values are not changed before the whole layer is gone
                 !// but we will only accept layers thicker than SFLIM.
                 if (BSNOWL(L,II,JJ,MP) .le. SFLIM) then
                    BSNOWL(L,II,JJ,MP) = 0._r8
                    LSNW_IJ(II,JJ,MP) = 0
                    BBCFFC(L,II,JJ,MP) = 0._r8
                    BBCBFC(L,II,JJ,MP) = 0._r8
                    BBCBIO(L,II,JJ,MP) = 0._r8
                    tSinceLastSF(II,JJ,MP) = 1.e5_r8
                 end if
              end if


              !// Nothing left to do on this grid box, all that could be
              !// evaporated has evaporated. Go to next (II,JJ,MP), by using
              !// cycle; this tries first next L, which is 0, and then we
              !// get to next (II,JJ,MP).
              cycle
           end if !// if (L .eq. 1) then


           !// We only get here if L > 1
           !// ------------------------------------------------------------
           !// When L>1, we have one or more layers beneeth the topmost thin
           !// layer, and these layers are subject to evaporation.
           !// The layer beneeth (L-1) has two options
           !//   Case 1: Melt/evap > snow amount is in the layer:
           !//           Whole layer will melt/evap.
           !//   Case 2: Melt/evap < snow amount is in the layer:
           !//           Melt/evap a fraction of the layer.

           if (evapmelt .ge. BSNOWL(L-1,II,JJ,MP)) then
              !// The whole layer L-1 will melt/evaporate.

              !// Physically, layer L-1 melts and its BC merges with
              !// layer L. But since layer L-1 is consequently removed,
              !// layer L is moved downwards, becoming the new L-1.
              BBCFFC(L-1,II,JJ,MP) = BBCFFC(L-1,II,JJ,MP) &
                                     + BBCFFC(L,II,JJ,MP)
              BBCFFC(L,II,JJ,MP) = 0._r8

              BBCBFC(L-1,II,JJ,MP) = BBCBFC(L-1,II,JJ,MP) &
                                     + BBCBFC(L,II,JJ,MP)
              BBCBFC(L,II,JJ,MP) = 0._r8

              BBCBIO(L-1,II,JJ,MP) = BBCBIO(L-1,II,JJ,MP) &
                                     + BBCBIO(L,II,JJ,MP)
              BBCBIO(L,II,JJ,MP) = 0._r8

              !// The whole layer L-1 evaporates; store the amount of snow
              !// evaporated and move L downwards
              current_evapmelt = BSNOWL(L-1,II,JJ,MP)
              BSNOWL(L-1,II,JJ,MP) = BSNOWL(L,II,JJ,MP)
              BSNOWL(L,II,JJ,MP) = 0._r8

              !// We have one layer less
              LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1

              !// DEBUG
              if (LDEBUG_BCsnow) then
                 if (BSNOWL(L-1,II,JJ,MP) .lt. 1.e-10_r8) then
                    write(6,'(a,i5,2es16.6)') f90file//':'//subr// &
                         ': thin layer 1: ',L,BSNOWL(L-1,II,JJ,MP),evapmelt
                    stop
                 end if

                 !// Safety check
                 !// It should not be possible to get LSNW_IJ = 0 here.
                 if (LSNW_IJ(II,JJ,MP) .eq. 0) then
                    write(6,'(a,i5,2es16.6)') f90file//':'//subr// &
                         ': LSNW_IJ zero when evap/melt'
                    write(6,'(a,2i5,2es16.6)') '  I,ILMX,evapmelt,BSNOWL: ', &
                         L,ILMX,evapmelt,BSNOWL(L,II,JJ,MP)
                    stop
                 end if

                 if (LSNW_IJ(II,JJ,MP) .eq. 1 .and. .not.(L .eq. 2)) then
                    write(6,'(a)') f90file//':'//subr// &
                         ': Inconsistent: LSNW_IJ=1 and L!=2'
                    write(6,'(a,7i5)') '  I,J,II,JJ,MP,L,LSNW_IJ: ', &
                         I,J,II,JJ,MP,L,LSNW_IJ(II,JJ,MP)
                    write(6,'(a)') '  *** Should have had L=2 when LSNW_IJ=1 ***'
                    stop
                 end if
              end if


              !// Account for melt/evap, i.e. what is left to melt/evap
              !// Note that this can become slightly negative; it will then
              !// be treated as a source of snow in the layer below (see
              !// further down the code).
              evapmelt = evapmelt - current_evapmelt

           else

              !// Melt/evaporate partial layer
              if (evapmelt .gt. 0._r8) then

                 !// Important:
                 !// When here, we always have L > 1 

                 !// When the snow melts, BC remains at the top of the
                 !// snow. This BC is added to the thinsnow top layer.
                 !// (Have tested that evapmelt < BSNOWL)
                 evapfrac = evapmelt / BSNOWL(L-1,II,JJ,MP)
                 if (LDEBUG_BCsnow) then
                    if (evapfrac .gt. 1._r8) then
                       write(6,'(a,es16.6)') f90file//':'//subr// &
                            ': Wrong evapfrac 1',evapfrac
                       write(6,'(a,2es16.6)') '  Evapmelt/bsnowl', &
                            evapmelt,BSNOWL(L-1,II,JJ,MP)
                       stop
                    end if
                 end if
                 BBCFFC(L,II,JJ,MP) = BBCFFC(L,II,JJ,MP) &
                                      + evapfrac * BBCFFC(L-1,II,JJ,MP)
                 BBCFFC(L-1,II,JJ,MP) = &
                      BBCFFC(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                 BBCBFC(L,II,JJ,MP) = BBCBFC(L,II,JJ,MP) &
                                      + evapfrac * BBCBFC(L-1,II,JJ,MP)
                 BBCBFC(L-1,II,JJ,MP) = &
                      BBCBFC(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                 BBCBIO(L,II,JJ,MP) = BBCBIO(L,II,JJ,MP) &
                                      + evapfrac * BBCBIO(L-1,II,JJ,MP)
                 BBCBIO(L-1,II,JJ,MP) = &
                      BBCBIO(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                 !// Reduce snow layer depth of L-1
                 BSNOWL(L-1,II,JJ,MP) = BSNOWL(L-1,II,JJ,MP) - evapmelt

                 !// Must account for evapmelt, which is now zero
                 evapmelt = 0._r8

                 !// Need to check thickness of layer L-1; if it is smaller
                 !// than SFLIM we should remove it.
                 if (BSNOWL(L-1,II,JJ,MP) .le. SFLIM) then
                    !// Move BC mass from L to L-1 and reduce layers
                    BBCFFC(L-1,II,JJ,MP) = BBCFFC(L-1,II,JJ,MP) &
                                           + BBCFFC(L,II,JJ,MP)
                    BBCBFC(L-1,II,JJ,MP) = BBCBFC(L-1,II,JJ,MP) &
                                           + BBCBFC(L,II,JJ,MP)
                    BBCBIO(L-1,II,JJ,MP) = BBCBIO(L-1,II,JJ,MP) &
                                           + BBCBIO(L,II,JJ,MP)
                    !// Shift layers down
                    SNL = LSNW_IJ(II,JJ,MP)
                    BBCFFC(L:(SNL-1),II,JJ,MP) = BBCFFC((L+1):SNL,II,JJ,MP)
                    BBCBFC(L:(SNL-1),II,JJ,MP) = BBCBFC((L+1):SNL,II,JJ,MP)
                    BBCBIO(L:(SNL-1),II,JJ,MP) = BBCBIO((L+1):SNL,II,JJ,MP)
                    BSNOWL((L-1):(SNL-1),II,JJ,MP) = BSNOWL(L:SNL,II,JJ,MP)
                    BBCFFC(SNL,II,JJ,MP) = 0._r8
                    BBCBFC(SNL,II,JJ,MP) = 0._r8
                    BBCBIO(SNL,II,JJ,MP) = 0._r8
                    BSNOWL(SNL,II,JJ,MP) = 0._r8
                    LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1
                 end if

              else

                 !// Negative evaporation; add as simple snowfall in the
                 !// uppermost layer (which is L)
                 BSNOWL(L,II,JJ,MP) = BSNOWL(L,II,JJ,MP) - evapmelt
                 !// Must account for evapmelt, which is now zero
                 evapmelt = 0._r8

              end if !// if (evapmelt .gt. 0._r8) then

           end if !// if (evapmelt .gt. BSNOWL(L-1,II,JJ,MP)) then

           !// Check if evapmelt is zero; then exit loop
           if (evapmelt .eq. 0._r8) exit

        end do !// do L = ILMX, 1, -1


        !// DEBUG
        if (LDEBUG_BCsnow) then
           if (BSNOWL(1,II,JJ,MP) .eq. 0._r8 .and. &
                LSNW_IJ(II,JJ,MP) .gt. 0) then
              write(6,'(a,es16.6)') f90file//':'//subr// &
                   ': Wrong in evapmelt 1'
              do L = ILMM,1,-1
                 write(6,'(a,i5,es16.6)') '  BSNOWL:',L, BSNOWL(L,II,JJ,MP)
              end do
              write(6,'(a,i5)') '  LSNW_IJ:',LSNW_IJ(II,JJ,MP)
              stop
           end if
           if (minval(BSNOWL(:,II,JJ,MP)) .lt. 0._r8) then
              write(6,'(a,es16.6)') f90file//':'//subr// &
                   ': Wrong in evapmelt 2 (negative BSNOWL)'
              do L = ILMM,1,-1
                 write(6,'(a,i5,es16.6)') '  BSNOWL:',L, BSNOWL(L,II,JJ,MP)
              end do
              stop
           end if
        end if

        !// If all is evaporated, there is no snowlayer anymore, and we
        !// have to build a new one next time it snows.
        if (LSNW_IJ(II,JJ,MP) .eq. 0) tSinceLastSF(II,JJ,MP) = 1.e5_r8


      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine bcsnow_meltevap_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_seaice_ij(NMET,NOPS,NSUB,DT,MP)
    !// --------------------------------------------------------------------
    !// Snowdepth, snowmelt and evaporation on sea ice is not included in
    !// the ECMWF model. We parameterise this in the following way:
    !// 
    !// 1. Linear reduction of the snowlayer on ice, 15. april until
    !//    20. june.
    !// 2. Then there is no snow on the ice.
    !//
    !// 3. If not sea ice cover 30%, then remove the snow and black carbon.
    !//
    !// 4. We allow snowdepth data in an ocean grid box, but check metdata
    !//    for consistency.
    !//
    !// Originally written by Ragnhild B. Skeie for CTM2, modified to
    !// fit CTM3 in March-April 2013.
    !//
    !//
    !// Amund Sovde, March-April 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, YDGRD, JDAY, &
         NRMETD, NROPSM, PLAND
    use cmn_met, only: CI, SD
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NMET,NOPS,NSUB,MP
    real(r8), intent(in) :: DT

    !// Locals
    integer :: I, J, L, II, JJ, ILMX, SNL

    integer :: springday, springend, summermid, summerend, KDAY
    real(r8)  :: OLDBSNOWL(ILMM)

    real(r8)  :: evapmelt, evapfrac, tempSNOW, totBSNOWL, current_evapmelt
    !// Land-Sea Mask, calculated from PLAND
    logical :: LSEA_BOX
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_seaice_ij'
    !//---------------------------------------------------------------------


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Separate treatment for NH and SH; get day and
      !// info about spring and summer.
      call bcsnow_getspringsummer(YDGRD(J), JDAY, &
           springday,springend,summermid,summerend,KDAY)

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// No adjustment if no snow
        if (LDEBUG_BCsnow) then
           if (LSNW_IJ(II,JJ,MP) .eq. 0) then
              if (sum(BSNOWL(:,II,JJ,MP)) .ne. 0._r8) then
                 write(6,'(a,es16.6)') f90file//':'//subr// &
                      ': No snow layer, but snow???'
                 do L = ILMM,1,-1
                    write(6,'(a,i5,es16.6)') '  BSNOWL:',L, BSNOWL(L,II,JJ,MP)
                 end do
                 write(6,'(a,es16.6)') '  sum(BSNOWL):', sum(BSNOWL(:,II,JJ,MP))
                 stop
              end if
              cycle
           end if
        end if

        !// ECMWF Land-Sea Mask: 0 if land covers more than 50% of grid box,
        !//                      otherwise 1 (i.e. sea).
        !// Instead of using land-sea mask, we use land fraction (PLAND) and
        !// assume sea if PLAND < seafraclim
        LSEA_BOX = PLAND(I,J) .lt. seafraclim

        !// Cycle if not sea ice (i.e. sea plus snow depth)
        if (.not. (LSEA_BOX .and. (SD(I,J) .eq. 0._r8))) cycle

        !// For sea ice < 30%, remove BCsnow
        if (CI(I,J) .lt. 0.3_r8) then
           BSNOWL(:,II,JJ,MP) = 0._r8
           LSNW_IJ(II,JJ,MP) = 0
           BBCFFC(:,II,JJ,MP) = 0._r8
           BBCBFC(:,II,JJ,MP) = 0._r8
           BBCBIO(:,II,JJ,MP) = 0._r8
           !// Reset tSinceLastSF
           tSinceLastSF(II,JJ,MP) = 1.e5_r8
           !// Go to next II/JJ
           cycle
        end if !// if (CI(I,J) .lt. 0.3_r8) then

        !// Keep track of snow before correction
        OLDBSNOWL(:) = BSNOWL(:,II,JJ,MP)

        if (KDAY .ge. springday .and. KDAY .lt. springend) then
           if (KDAY .eq. springday .and. &
                NMET .eq. 1 .and. NOPS .eq. 1 .and. NSUB .eq. 1) then
              !// The variable springmelt is the amount of snow that
              !// needs to be evaporated/melted each second, to ensure
              !// a linear decrease reaching zero at the end of spring.
              !// Initialization is done at the first day of spring, with
              !// units snow/(timesteps to springend).
              springmelt(II,JJ,MP) = sum(BSNOWL(:,II,JJ,MP)) &
                                   / (real(springend - KDAY, r8)*86400._r8)
           end if

           !// Redefine/update springmelt if there has been snowfall
           if (tSinceLastSF(II,JJ,MP) .eq. 0._r8) then
              springmelt(II,JJ,MP) = sum(BSNOWL(:,II,JJ,MP)) &
                                       / (real(springend - KDAY, r8)*86400._r8)
           end if

           !// Melting this time step
           evapmelt = springmelt(II,JJ,MP) * DT

           !// Do melting as for regular melting

           !// Number of layers
           ILMX = LSNW_IJ(II,JJ,MP)
           !// Loop downwards
           do L = ILMX, 1, -1

              if (L .eq. 1) then
                 !// There is no layer beneeth to evaporate;
                 !// do evaporation of layer L only.

                 !// If not everything melts/evaporates, we keep a thin layer;
                 !// otherwise, all snow and BC is gone.
                 !// (Note that this part works on layer L, which is 1.)
                 if (evapmelt .ge. BSNOWL(L,II,JJ,MP)) then
                    !// All gone, no snow; BC is lost on ground
                    BSNOWL(L,II,JJ,MP) = 0._r8
                    LSNW_IJ(II,JJ,MP) = 0
                    BBCFFC(L,II,JJ,MP) = 0._r8
                    BBCBFC(L,II,JJ,MP) = 0._r8
                    BBCBIO(L,II,JJ,MP) = 0._r8
                    tSinceLastSF(II,JJ,MP) = 1.e5_r8
                    evapmelt = 0._r8
                 else if (evapmelt .lt. 0._r8) then
                    !// Could possibly add to snow layer; skip for L=1
                    !// (as in CTM2)
                    evapmelt = 0._r8
                 else
                    !// Melt part of layer
                    BSNOWL(L,II,JJ,MP) = BSNOWL(L,II,JJ,MP) - evapmelt
                    LSNW_IJ(II,JJ,MP) = 1
                    evapmelt = 0._r8
                    !// BC values are not changed before the whole layer is gone
                    !// but we will only accept layers thicker than SFLIM.
                    if (BSNOWL(L,II,JJ,MP) .le. SFLIM) then
                       write(6,'(a,2i7)') f90file//':'//subr//': SEAICE EVAPORATION: '// &
                            'ASSUME ZERO SNOW',I,J
                       BSNOWL(L,II,JJ,MP) = 0._r8
                       LSNW_IJ(II,JJ,MP) = 0
                       BBCFFC(L,II,JJ,MP) = 0._r8
                       BBCBFC(L,II,JJ,MP) = 0._r8
                       BBCBIO(L,II,JJ,MP) = 0._r8
                       tSinceLastSF(II,JJ,MP) = 1.e5_r8
                    end if
                 end if

                 !// Nothing left to do on this grid box.
                 cycle
              end if !// if (L .eq. 1) then


              !// We only get here if L > 1
              if (evapmelt .ge. BSNOWL(L-1,II,JJ,MP)) then
                 !// The whole layer L-1 will melt/evaporate.
                 BBCFFC(L-1,II,JJ,MP) = BBCFFC(L-1,II,JJ,MP) &
                                        + BBCFFC(L,II,JJ,MP)
                 BBCFFC(L,II,JJ,MP) = 0._r8

                 BBCBFC(L-1,II,JJ,MP) = BBCBFC(L-1,II,JJ,MP) &
                                        + BBCBFC(L,II,JJ,MP)
                 BBCBFC(L,II,JJ,MP) = 0._r8

                 BBCBIO(L-1,II,JJ,MP) = BBCBIO(L-1,II,JJ,MP) &
                                        + BBCBIO(L,II,JJ,MP)
                 BBCBIO(L,II,JJ,MP) = 0._r8

                 !// The whole layer L-1 evaporates; store the amount of snow
                 !// evaporated and move L downwards
                 current_evapmelt = BSNOWL(L-1,II,JJ,MP)
                 BSNOWL(L-1,II,JJ,MP) = BSNOWL(L,II,JJ,MP)
                 BSNOWL(L,II,JJ,MP) = 0._r8

                 !// We have one layer less
                 LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1

                 if (LDEBUG_BCsnow) then
                    if (LSNW_IJ(II,JJ,MP) .eq. 0) then
                       write(6,'(a)') f90file//':'//subr//': LSNW_IJ seaice 1'
                       write(6,'(2i5,2es16.6)') L,ILMX,evapmelt,BSNOWL(L,II,JJ,MP)
                       stop
                    end if
                 end if

                 !// Account for melt/evap, i.e. what is left to melt/evap
                 evapmelt = evapmelt - current_evapmelt

              else

                 !// Melt/evaporate partial layer
                 if (evapmelt .gt. 0._r8) then
                    !// When here, we always have L > 1
                    evapfrac = evapmelt / BSNOWL(L-1,II,JJ,MP)
                    if (LDEBUG_BCsnow) then
                       if (evapfrac .gt. 1._r8) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': Wrong evapfrac seaice 1: ',evapfrac
                          write(6,'(2es16.6)') evapmelt, BSNOWL(L-1,II,JJ,MP)
                          stop
                       end if
                    end if

                    BBCFFC(L,II,JJ,MP) = BBCFFC(L,II,JJ,MP) &
                                         + evapfrac * BBCFFC(L-1,II,JJ,MP)
                    BBCFFC(L-1,II,JJ,MP) = &
                         BBCFFC(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                    BBCBFC(L,II,JJ,MP) = BBCBFC(L,II,JJ,MP) &
                                         + evapfrac * BBCBFC(L-1,II,JJ,MP)
                    BBCBFC(L-1,II,JJ,MP) = &
                         BBCBFC(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                    BBCBIO(L,II,JJ,MP) = BBCBIO(L,II,JJ,MP) &
                                         + evapfrac * BBCBIO(L-1,II,JJ,MP)
                    BBCBIO(L-1,II,JJ,MP) = &
                         BBCBIO(L-1,II,JJ,MP) * (1._r8 - evapfrac)

                    !// Reduce snow layer depth of L-1
                    BSNOWL(L-1,II,JJ,MP) = BSNOWL(L-1,II,JJ,MP) - evapmelt

                    !// Must account for evapmelt, which is now zero
                    evapmelt = 0._r8

                    !// Need to check thickness of layer L-1; if it is smaller
                    !// than SFLIM we should remove it.
                    if (BSNOWL(L-1,II,JJ,MP) .le. SFLIM) then
                       !// Move BC mass from L-1 to L and reduce layers
                       BBCFFC(L-1,II,JJ,MP) = BBCFFC(L-1,II,JJ,MP) &
                                              + BBCFFC(L,II,JJ,MP)
                       BBCBFC(L-1,II,JJ,MP) = BBCBFC(L-1,II,JJ,MP) &
                                              + BBCBFC(L,II,JJ,MP)
                       BBCBIO(L-1,II,JJ,MP) = BBCBIO(L-1,II,JJ,MP) &
                                              + BBCBIO(L,II,JJ,MP)

                       !// Shift layers down
                       SNL = LSNW_IJ(II,JJ,MP)
                       BBCFFC(L:(SNL-1),II,JJ,MP) = BBCFFC((L+1):SNL,II,JJ,MP)
                       BBCBFC(L:(SNL-1),II,JJ,MP) = BBCBFC((L+1):SNL,II,JJ,MP)
                       BBCBIO(L:(SNL-1),II,JJ,MP) = BBCBIO((L+1):SNL,II,JJ,MP)
                       BSNOWL((L-1):(SNL-1),II,JJ,MP) = BSNOWL(L:SNL,II,JJ,MP)
                       BBCFFC(SNL,II,JJ,MP) = 0._r8
                       BBCBFC(SNL,II,JJ,MP) = 0._r8
                       BBCBIO(SNL,II,JJ,MP) = 0._r8
                       BSNOWL(SNL,II,JJ,MP) = 0._r8
                       LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1
                       write(6,'(a,3i7)') f90file//':'//subr// &
                            ': SEAICE REMOVING LAYER',L-1,I,J
                    end if

                 else

                    !// Negative evaporation; skip this here.
                    evapmelt = 0._r8

                 end if !// if (evapmelt .gt. 0._r8) then

              end if !// if (evapmelt .gt. BSNOWL(L-1,II,JJ,MP)) then

              !// Check if evapmelt is zero; then exit loop
              if (evapmelt .eq. 0._r8) exit

           end do !// do L = ILMX, 1, -1


           !// DEBUG
           if (LDEBUG_BCsnow) then
              if (BSNOWL(1,II,JJ,MP) .eq. 0._r8 .and. &
                   LSNW_IJ(II,JJ,MP) .gt. 0) then
                 write(6,'(a,3i7)') f90file//':'//subr// &
                      ': Wrong in seaice 1'
                 do L = ILMM, 1, -1
                    write(6,'(a,i5,es16.6)') 'BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
                 end do
                 write(6,'(a,es16.6)') 'LSNW_IJ: ',LSNW_IJ(II,JJ,MP)
                 stop
              end if
              if (minval(BSNOWL(:,II,JJ,MP)) .lt. 0._r8) then
                 write(6,'(a,3i7)') f90file//':'//subr// &
                      ': Wrong in seaice 2'
                 do L = ILMM, 1, -1
                    write(6,'(a,i5,es16.6)') 'BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
                 end do
                 stop
              end if
           end if

        else if (KDAY .ge. springend .and. KDAY .le. summermid) then
           !// Is there snow? Reduce it to very small depth.
           if (sum(BSNOWL(:,II,JJ,MP)) .gt. 0._r8) then
              BSNOWL(:,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = 1
              BSNOWL(1,II,JJ,MP) = thinsnow
              tSinceLastSF(II,JJ,MP) = 1.e5_r8
              BBCFFC(1,II,JJ,MP) = sum(BBCFFC(:,II,JJ,MP))
              BBCFFC(2:ILMM,II,JJ,MP) = 0._r8
              BBCBFC(1,II,JJ,MP) = sum(BBCBFC(:,II,JJ,MP))
              BBCBFC(2:ILMM,II,JJ,MP) = 0._r8
              BBCBIO(1,II,JJ,MP) = sum(BBCBIO(:,II,JJ,MP))
              BBCBIO(2:ILMM,II,JJ,MP) = 0._r8
           end if
           !// Zero arrays after summer
           if (KDAY .eq. summermid .and. NMET .eq. NRMETD .and. &
                NOPS .eq. NROPSM) then
              !// Now we remove the BC from last year, to make it possible
              !// to run several years
              BBCFFC(:,II,JJ,MP) = 0._r8
              BBCBFC(:,II,JJ,MP) = 0._r8
              BBCBIO(:,II,JJ,MP) = 0._r8
           end if

           !// DEBUG
           if (LDEBUG_BCsnow) then
              if (BSNOWL(1,II,JJ,MP) .eq. 0._r8 .and. &
                   LSNW_IJ(II,JJ,MP) .gt. 0) then
                 write(6,'(a,3i7)') f90file//':'//subr// &
                      ': Wrong in seaice 3'
                 do L = ILMM, 1, -1
                    write(6,'(a,i5,es16.6)') 'New BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
                 end do
                 write(6,'(a,es16.6)') 'LSNW_IJ: ',LSNW_IJ(II,JJ,MP)
                 do L = ILMM, 1, -1
                    write(6,'(a,i5,es16.6)') 'Old BSNOWL: ',L, OLDBSNOWL(L)
                 end do
                 stop
              end if
              if (minval(BSNOWL(:,II,JJ,MP)) .lt. 0._r8) then
                 write(6,'(a,3i7)') f90file//':'//subr// &
                      ': Wrong in seaice 4'
                 do L = ILMM, 1, -1
                    write(6,'(a,i5,es16.6)') 'BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
                 end do
                 stop
              end if
           end if
        end if !// if (KDAY .ge. springday .and. KDAY .lt. springend) then

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine bcsnow_seaice_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_adjustment_ij(NMET,NOPS,NSUB,MP)
    !// --------------------------------------------------------------------
    !// Adjustment of snowdepth from BSNOWL to metdata snowdepth SD.
    !//
    !// Amund Sovde, March-April 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, YDGRD, JDAY, PLAND
    use cmn_met, only: CI, SD
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NMET,NOPS,NSUB,MP
    !// Locals
    integer :: I, J, L, II, JJ, ILMX, SNL, LLL
    integer :: springday, springend, summermid, summerend, KDAY
    real(r8)  :: tempSNOW, thinfrac, frac, totBSNOWL
    real(r8)  :: OLDBSNOWL(ILMM)
    logical :: LSEA_BOX
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_adjustment_ij'
    !//---------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Get KDAY and summerend
      call bcsnow_getspringsummer(YDGRD(J), JDAY, &
           springday,springend,summermid,summerend,KDAY)

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Skip adjustment if there is no snow in neither model nor metdata
        if (LSNW_IJ(II,JJ,MP) .eq. 0 .and. SD(I,J) .eq. 0._r8) then
           !// DEBUG
           if (LDEBUG_BCsnow) then
              if (sum(BBCFFC(:,II,JJ,MP)) .gt. 0._r8) then
                 write(6,'(a)') f90file//':'//subr// &
                      ': NO SNOW, but BBCFFC:'
                 do L = ILMM, 1, -1
                    write(6,'(a)') '  BBCFFC: ',L,BBCFFC(L,II,JJ,MP)
                 end do
                 stop
              end if
           end if
           cycle
        end if

        !// Save old value before adjustment
        OLDBSNOWL(:) = BSNOWL(:,II,JJ,MP)

        !// To be able to run continius runs we will remove to BC in
        !// gridcels where there is whole-year snow 1st sept in NH and
        !// 1st March in SH.
        if (KDAY .eq. summerend .and. NMET .eq. 1 .and. NOPS .eq. 1 &
             .and. NSUB .eq. 1) then
           BBCFFC(:,II,JJ,MP) = 0._r8
           BBCBFC(:,II,JJ,MP) = 0._r8
           BBCBIO(:,II,JJ,MP) = 0._r8
           tempSNOW = sum(BSNOWL(:,II,JJ,MP))
           if (tempSNOW .gt. 0._r8) then
              BSNOWL(1,II,JJ,MP) = tempSNOW
              BSNOWL(2:ILMM,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = 1
              if (tempSNOW .gt. thinsnowthreshold) then
                 BSNOWL(1,II,JJ,MP) = BSNOWL(1,II,JJ,MP) - thinsnow
                 BSNOWL(2,II,JJ,MP) = thinsnow
                 LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
              end if
           else
              !// Fail-safe
              LSNW_IJ(II,JJ,MP) = 0
           end if
        end if !// if (KDAY .eq. summerend .and. NMET==1) then


        !// Adjust depth of snow layers to match metdata snow depth.

        !// Sea box?
        LSEA_BOX = PLAND(I,J) .lt. seafraclim


        !// Only scale levels below uppermost layer, for snow depth
        !// less than 0.2 m water equivalents.
        if (LSNW_IJ(II,JJ,MP) .gt. 0) then
           if (SD(I,J) .gt. 0._r8) then
              !// Adjust snow
              if (LSNW_IJ(II,JJ,MP) .gt. 1) then
                 !// Option 1: metdata > uppermost thin layer
                 !//           => scale layers below
                 !// Option 2: metdata <= uppermost thin layer
                 !//           => scale all layers (will always scale down)
                 SNL = LSNW_IJ(II,JJ,MP)
                 tempSNOW =  min(SD(I,J), metdataSDmax)
                 if ( (tempSNOW - BSNOWL(SNL,II,JJ,MP)) .gt. 0._r8) then
                    !// Scale levels below thinlayer (but not thinlayer)
                    frac = (tempSNOW - BSNOWL(SNL,II,JJ,MP)) / &
                             sum(BSNOWL(1:SNL-1,II,JJ,MP))
                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (frac .le. 0._r8) then
                          write(6,'(a,es16.6,5i7)') f90file//':'//subr// &
                               ': frac 1: ',frac,i,j,ii,jj,mp
                          write(6,'(a,3es16.6)') &
                               '  tempSNOW, BSNOWL(SNL), sum(BSNOWL(1:SNL)', &
                               tempSNOW, BSNOWL(SNL,II,JJ,MP), &
                               sum(BSNOWL(1:SNL-1,II,JJ,MP))
                          stop
                       end if
                    end if
                    BSNOWL(1:SNL-1,II,JJ,MP) = BSNOWL(1:SNL-1,II,JJ,MP) * frac
                 else
                    !// Scale all levels (will always scale down)
                    frac = tempSNOW / sum(BSNOWL(1:SNL,II,JJ,MP))
                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (frac .gt. 1._r8) then
                          write(6,'(a,es16.6,5i7)') f90file//':'//subr// &
                               ': frac 2: ',frac,i,j,ii,jj,mp
                         stop
                       end if
                    end if
                    BSNOWL(1:SNL,II,JJ,MP) = BSNOWL(1:SNL,II,JJ,MP) * frac
                 end if !// if ( (tempSNOW - BSNOWL(SNL,II,JJ,MP)) .gt. 0._r8) then

              else
                 !// One-layer only; update to match SD
                 BSNOWL(1,II,JJ,MP) = min(SD(I,J), metdataSDmax)
                 !// Add new layer for drydep
                 if (BSNOWL(1,II,JJ,MP) .gt. thinsnowthreshold) then
                    thinfrac = thinsnow / BSNOWL(1,II,JJ,MP)
                    !// DEBUG
                    if (LDEBUG_BCsnow) then
                       if (thinfrac .gt. 1._r8) then
                          write(6,'(a,es16.6)') f90file//':'//subr// &
                               ': thinfrac wrong 1', thinfrac
                          stop
                       end if
                    end if
                    BSNOWL(2,II,JJ,MP) = thinsnow
                    BSNOWL(1,II,JJ,MP) = BSNOWL(1,II,JJ,MP) - thinsnow
                    BBCFFC(2,II,JJ,MP) = thinfrac * BBCFFC(1,II,JJ,MP)
                    BBCBFC(2,II,JJ,MP) = thinfrac * BBCBFC(1,II,JJ,MP)
                    BBCBIO(2,II,JJ,MP) = thinfrac * BBCBIO(1,II,JJ,MP)
                    BBCFFC(1,II,JJ,MP) = (1._r8 - thinfrac) * BBCFFC(1,II,JJ,MP)
                    BBCBFC(1,II,JJ,MP) = (1._r8 - thinfrac) * BBCBFC(1,II,JJ,MP)
                    BBCBIO(1,II,JJ,MP) = (1._r8 - thinfrac) * BBCBIO(1,II,JJ,MP)
                    !// Increase by 1 levels (i.e. from 1 to 2)
                    LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
                 end if !// if (BSNOWL(1,II,JJ,MP) .gt. thinsnowthreshold) then
              end if !// if (LSNW_IJ(II,JJ,MP) .gt. 1) then

           else if (SD(I,J) .eq. 0._r8 .and. .not.LSEA_BOX) then
              !// Land only: Snow in the model, but not in metdata?
              !//            => Set model BCsnow to zero.
              !// Sea ice snow is treated in separate routine.
              BSNOWL(:,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = 0
              tSinceLastSF(II,JJ,MP) = 1.e5_r8
              BBCFFC(:,II,JJ,MP) = 0._r8
              BBCBFC(:,II,JJ,MP) = 0._r8
              BBCBIO(:,II,JJ,MP) = 0._r8
           end if

        end if !// if (LSNW_IJ(II,JJ,MP) .gt. 0) then

        !// Snow in metdata, but not in model: Put snow into BSNOWL.
        if (LSNW_IJ(II,JJ,MP) .eq. 0 .and. SD(I,J) .gt. 0._r8) then
           LSNW_IJ(II,JJ,MP) = 1
           SNL = LSNW_IJ(II,JJ,MP)

           BSNOWL(SNL,II,JJ,MP) = min(SD(I,J), metdataSDmax)
           !// Make a thin snow layer for drydeposition.
           if (BSNOWL(SNL,II,JJ,MP) .gt. thinsnowthreshold) then
              LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) + 1
              SNL = LSNW_IJ(II,JJ,MP)
              BSNOWL(SNL,II,JJ,MP) = thinsnow
              BSNOWL(SNL-1,II,JJ,MP) = BSNOWL(SNL-1,II,JJ,MP) - thinsnow
           end if
        end if

        !// Check snow depth (no need to check SD>metdataSDmax)
        if (SD(I,J) .gt. 0._r8 .and. SD(I,J) .lt. metdataSDmax) then
           if ( (SD(I,J) - sum(BSNOWL(:,II,JJ,MP))) .gt. 1.e-6_r8) then
              if (.not.LSEA_BOX) then
                 write(6,'(a,4es16.6)') f90file//':'//subr// &
                      ': snowdepth wrong ', &
                      SD(I,J), sum(BSNOWL(:,II,JJ,MP)), CI(I,J), PLAND(I,J)
              end if
           end if
        end if

        !// Final check in CTM2 (on top-most layer) is unnecessary in CTM3.

        !// Need to check thickness of layers; if it is smaller
        !// than SFLIM we should remove it.
        LLL = LSNW_IJ(II,JJ,MP) !// Loop may reduce LSNW_IJ
        do L = LLL-1, 2, -1
           if (BSNOWL(L,II,JJ,MP) .le. SFLIM) then
              SNL = LSNW_IJ(II,JJ,MP)
              !// Move BC mass from L to L-1 and remove layer L
              BBCFFC(L-1,II,JJ,MP) = BBCFFC(L-1,II,JJ,MP) + BBCFFC(L,II,JJ,MP)
              BBCBFC(L-1,II,JJ,MP) = BBCBFC(L-1,II,JJ,MP) + BBCBFC(L,II,JJ,MP)
              BBCBIO(L-1,II,JJ,MP) = BBCBIO(L-1,II,JJ,MP) + BBCBIO(L,II,JJ,MP)
              !// Shift layers above
              BBCFFC(L:(SNL-1),II,JJ,MP) = BBCFFC((L+1):SNL,II,JJ,MP)
              BBCBFC(L:(SNL-1),II,JJ,MP) = BBCBFC((L+1):SNL,II,JJ,MP)
              BBCBIO(L:(SNL-1),II,JJ,MP) = BBCBIO((L+1):SNL,II,JJ,MP)
              BSNOWL((L-1):(SNL-1),II,JJ,MP) = BSNOWL(L:SNL,II,JJ,MP)
              BBCFFC(SNL,II,JJ,MP) = 0._r8
              BBCBFC(SNL,II,JJ,MP) = 0._r8
              BBCBIO(SNL,II,JJ,MP) = 0._r8
              BSNOWL(SNL,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1
              !// Exit loop if LSNW_IJ(II,JJ,MP)==1
              if (LSNW_IJ(II,JJ,MP) .eq. 1) exit
           end if
        end do
        !// And finally check if layer 1 should be removed
        if (LSNW_IJ(II,JJ,MP) .gt. 0 .and. BSNOWL(1,II,JJ,MP) .le. SFLIM) then
           !// If only one layer, remove all
           if (LSNW_IJ(II,JJ,MP) .eq. 1) then
              BSNOWL(:,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = 0
              tSinceLastSF(II,JJ,MP) = 1.e5_r8
              BBCFFC(:,II,JJ,MP) = 0._r8
              BBCBFC(:,II,JJ,MP) = 0._r8
              BBCBIO(:,II,JJ,MP) = 0._r8
           else
              !// At least two layers; move BC mass from L to L+1 and reduce
              !// layers
              SNL = LSNW_IJ(II,JJ,MP)
              BBCFFC(1,II,JJ,MP) = BBCFFC(1,II,JJ,MP) + BBCFFC(2,II,JJ,MP)
              BBCBFC(1,II,JJ,MP) = BBCBFC(1,II,JJ,MP) + BBCBFC(2,II,JJ,MP)
              BBCBIO(1,II,JJ,MP) = BBCBIO(1,II,JJ,MP) + BBCBIO(2,II,JJ,MP)
              !// Shift layers
              BBCFFC(1:(SNL-1),II,JJ,MP) = BBCFFC(2:SNL,II,JJ,MP)
              BBCBFC(1:(SNL-1),II,JJ,MP) = BBCBFC(2:SNL,II,JJ,MP)
              BBCBIO(1:(SNL-1),II,JJ,MP) = BBCBIO(2:SNL,II,JJ,MP)
              BSNOWL(1:(SNL-1),II,JJ,MP) = BSNOWL(2:SNL,II,JJ,MP)
              BBCFFC(SNL,II,JJ,MP) = 0._r8
              BBCBFC(SNL,II,JJ,MP) = 0._r8
              BBCBIO(SNL,II,JJ,MP) = 0._r8
              BSNOWL(SNL,II,JJ,MP) = 0._r8
              LSNW_IJ(II,JJ,MP) = LSNW_IJ(II,JJ,MP) - 1
           end if
        end if


        !// DEBUG
        if (LDEBUG_BCsnow) then
           if (BSNOWL(1,II,JJ,MP) .eq. 0._r8 .and. &
                LSNW_IJ(II,JJ,MP) .gt. 0) then
              write(6,'(a,4es16.6)') f90file//':'//subr// &
                   ': Wrong in adjustment 3'
              do L = ILMM, 1, -1
                 write(6,'(a,i5,es16.6)') 'New BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
              end do
              write(6,'(a,es16.6)') 'LSNW_IJ: ',LSNW_IJ(II,JJ,MP)
              do L = ILMM, 1, -1
                 write(6,'(a,i5,es16.6)') 'Old BSNOWL: ',L, OLDBSNOWL(L)
              end do
              stop
           end if
           if (minval(BSNOWL(:,II,JJ,MP)) .lt. 0._r8) then
              write(6,'(a,4es16.6)') f90file//':'//subr// &
                   ': Wrong in adjustment 4'
              do L = ILMM, 1, -1
                 write(6,'(a,i5,es16.6)') 'BSNOWL: ',L, BSNOWL(L,II,JJ,MP)
              end do
              stop
           end if
        end if

      end do !// do I = MPBLKIB(MP),MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP),MPBLKJE(MP)
    !// --------------------------------------------------------------------
  end subroutine bcsnow_adjustment_ij
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine bcsnow_master(DT, NDAY, NMET, NOPS, NSUB, MP)
    !// --------------------------------------------------------------------
    !// This routine controls the main routines for BCsnow.
    !// It collects the deposited BC and distributes it in snow layers,
    !// then checks for evaporation or melting.
    !//
    !// Amund Sovde, March-April 2-013
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DT
    integer, intent(in) :: NDAY, NMET, NOPS, NSUB, MP
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_master'
    !//---------------------------------------------------------------------

    if (.not.LBCsnow) return

    call bcsnow_collect_ij(DT,MP)

    call bcsnow_meltevap_ij(DT,MP)

    call bcsnow_seaice_ij(NMET,NOPS,NSUB,DT,MP)

    call bcsnow_adjustment_ij(NMET,NOPS,NSUB,MP)

    !// --------------------------------------------------------------------
  end subroutine bcsnow_master
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine bcsnow_diag_dd(DDEP,COMPID,MP)
    !// --------------------------------------------------------------------
    !// Diagnose some amount of dry dep (DDEP) for component COMPID.
    !//
    !// Amund Sovde, September 2013
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DDEP(IDBLK,JDBLK)
    integer, intent(in) :: COMPID, MP
    !// Locals
    integer :: II,JJ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'bcsnow_diag_dd'
    !//---------------------------------------------------------------------

    if (.not. LBCsnow) return

    !// BB components
    if (COMPID .eq. BC_IDS(1) .or. COMPID .le. BC_IDS(2)) then
       !// Loop over latitude (J is global, JJ is block)
       do JJ = 1, JDBLK
          do II = 1, IDBLK
             !// Collect drydep values
             bcsnow_dd_bio(II,JJ,MP) = bcsnow_dd_bio(II,JJ,MP) + DDEP(II,JJ)
          end do
       end do
    end if

    !// FF components
    if (COMPID .eq. BC_IDS(3) .or. COMPID .eq. BC_IDS(4)) then
       !// Loop over latitude (J is global, JJ is block)
       do JJ = 1, JDBLK
          do II = 1, IDBLK
             !// Collect drydep values
             bcsnow_dd_ffc(II,JJ,MP) = bcsnow_dd_ffc(II,JJ,MP) + DDEP(II,JJ)
          end do
       end do
    end if

    !// BF components
    if (COMPID .eq. BC_IDS(5) .or. COMPID .eq. BC_IDS(6)) then
       !// Loop over latitude (J is global, JJ is block)
       do JJ = 1, JDBLK
          do II = 1, IDBLK
             !// Collect drydep values
             bcsnow_dd_bfc(II,JJ,MP) = bcsnow_dd_bfc(II,JJ,MP) + DDEP(II,JJ)
          end do
       end do
    end if

    !// --------------------------------------------------------------------
  end subroutine bcsnow_diag_dd
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
end module bcoc_oslo
!//=========================================================================

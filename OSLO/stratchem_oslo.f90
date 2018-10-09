!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Oslo stratospheric chemistry.
!//=========================================================================
module stratchem_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: stratchem_oslo
  !// DECRIPTION: Routine for driving the Oslo stratospheric chemistry.
  !//
  !// Contains:
  !//   subroutine oslochem_strat
  !//   subroutine read_oslo2d2
  !//   subroutine update_strat_boundaries
  !//   subroutine set_fam_in_trop
  !//
  !// Ole Amund Sovde, November 2014 (from .f to .f90)
  !//                  Updated March 2009 - October 2009
  !//                  September 2008
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_strat(BTT,BJV,BAIR_MOLEC,BVOL,BTEM,BEMIS,BTTBCK,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Stratospheric chemistry. Prepares the integration and calls
    !// integration routine for each column.
    !//
    !// This means that for each I,J (or II,JJ) the integration is done
    !// from the LMTROP(I,J)+1 and up to LPAR-1.
    !//
    !// Ole Amund Sovde, October 2008
    !//                  Updated March 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IDBLK, JDBLK, MPBLK, IPAR, JPAR, &
         LPAR, LPARW, LWEPAR, NPAR, NOTRPAR, TRACER_ID_MAX, &
         LEMISDEP_INCHEM, LSOA
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, ETAA, ETAB
    use cmn_chem, only: TMASS
    use cmn_fjx, only: JVN_
    use cmn_oslo, only: LMTROP, trsp_idx, chem_idx, &
         Xtrsp_idx, Xchem_idx, XTMASS, XSTT
    use cmn_met, only: P, Q
    use cmn_parameters, only: M_AIR, AVOGNR, R_AIR
    use chem_oslo_rates, only: TCRATE_TP_IJ_STR, TCRATE_HET_IJ, &
         TEMPRANGE, MINTEMP
    use pchemc_str_ij
    use utilities_oslo
    !use oc_pscvariables, only: SPS_PARTAREA, PSC1, PSC2
    use psc_microphysics, only: get_psc12_sad
    use emisdep4chem_oslo, only: getEMISX
    use strat_h2o, only: zc_strh2o, LOLD_H2OTREATMENT
    use soa_oslo, only: SOA_v9_separate, soa_diag_separate
    use diagnostics_general, only: nchemdiag, save_chemPL, save_chemOxPL
    !// ------------------------------------------------------------------
    implicit none
    !// ------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DTCHM
    integer, intent(in) :: MP
    real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK,MPBLK),intent(in) :: BJV
    real(r8), dimension(LPAR,IDBLK,JDBLK,MPBLK), intent(in) :: &
         BAIR_MOLEC, BVOL
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BEMIS
    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTTBCK

    !// Locals
    !// QSSA variables
    real(r8) :: ST_STR,QLIN_STR,QTEST_STR
    !// Pressure variables
    real(r8) :: P1, P2
    !// Column arrays -----------------------------------------------------
    !// Pressure & layer thickness
    real(r8),dimension(LPAR) :: PMIDL, PUEDG

    !// Air density [molecules / cm3]
    real(r8),dimension(LPAR) :: AIR_MOLEC

    !// Temperature [K]
    real(r8),dimension(LPAR) :: TEMP

    !// Volume [m3]
    real(r8),dimension(LPAR) :: DV

    !// PSC1, PSC2 and aerosol SAD [cm2/cm3]
    real(r8),dimension(LPAR) :: PSC1L, PSC2L, SADL

    !// Tracer concentration [molecules / cm3]
    real(r8),dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL, ZC0

    !// Photolytical reaction rates
    real(r8),dimension(JVN_,LPAR) :: JVALUES

    !// Emissions [molec/cm3/s]
    !// Includes ALL emissions, fetched from BEMIS
    real(r8), dimension(TRACER_ID_MAX,LPAR) :: EMISX

    !// Correct O3 BTTBCK
    real(r8),dimension(LPAR) :: deltaO3

    !// Chemistry diagnoses
    real(r8), dimension(nchemdiag,TRACER_ID_MAX,LPAR) :: CHEMLOSS
    real(r8), dimension(nchemdiag,TRACER_ID_MAX,LPAR) :: CHEMPROD
    real(r8), dimension(LPAR) :: OxCHEMLOSS
    real(r8), dimension(LPAR) :: OxCHEMPROD
    !// --------------------------------------------------------------------

    !// Reaction rates dependent on pressure and/or temperature
    real(r8), dimension(LPAR) :: &
         r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
         r_oh_hno3,   r_oh_co_a,   r_oh_co_b,     r_ho2_ho2_tot, &
         r_cl_o2_m,   r_cloo_heat, r_clo_clo_m,   r_cl2o2_m, &
         r_oh_no2_m,  r_op_no2_m,  r_oh_oh_m,     r_oh_no_m, &
         r_no2_clo_m, r_op_o2_m,   r_o2_h_m,      r_no2_bro_m, &
         r_no_ho2_b, &
         spsADPAR,spsADPARBK, &
         spsACPAR,spsACPARBK, &
         spsBDPAR,spsBDPARBK, &
         spsBCPAR,spsBCPARBK, &
         spsECPAR,spsECPARBK, &
         spsFDPAR,spsFDPARBK, &
         spsFCPAR,spsFCPARBK, &
         spsGCPAR,spsGCPARBK, &
         spsBHPAR

    !// --------------------------------------------------------------------

    !// Max chemical time step in stratosphere is 5min = 300sec.
    real(r8),parameter :: DTCHM_MAX = 300._r8

    !// Time step for stratospheric chemistry
    real(r8) :: DTCH_STR

    !// Number of iterations in stratosphere
    integer :: NCHEM_ITER

    !// For looping
    integer :: I, J, L, II, JJ

    !// Max height of tropopause & the actual tropopause level
    integer :: LMTP, LMTPP1

    logical, parameter :: verbose=.false.
    !// --------------------------------------------------------------------

    !// Want maximum 5min time step
    !DTCH_STR = DTCHM_MAX

    !// Need a number of loops in chemistry to match the global chemical
    !// time step
    if (DTCHM .lt. DTCHM_MAX) then
      !// Global time step is shorter than DTCHM_MAX
      DTCH_STR = DTCHM
      NCHEM_ITER = 1
    else
      NCHEM_ITER = int(DTCHM/DTCHM_MAX + 0.5_r8)
      if (mod(DTCHM,DTCHM_MAX) .eq. 0._r8) then
        !// DTCH_STR is given by max time step.
        DTCH_STR = DTCHM_MAX
      else
        !// Must change DTCH_STR to match chemical iterations.
        DTCH_STR = DTCHM / real(NCHEM_ITER, r8)
      end if
    end if

    !// Set QSSA parameters for this time step
    ST_STR    = 10._r8 / DTCH_STR
    QLIN_STR  = 0.1_r8 / DTCH_STR
    QTEST_STR = 1.0_r8 / DTCH_STR

    if (DTCH_STR .gt. DTCHM_MAX) then
      write(6,'(a)') '*** ERROR Something went wrong when determining'// &
            ' the chemistry time step, The program is stopped !'
      stop
    end if
    if (verbose) write(6,'(a,f12.6,a3,i3,a3,i3)') &
         '* Strat: DTCH_STR/NCHEM_ITER/MP',DTCH_STR,' / ',NCHEM_ITER,' / ',MP


    !// Zero emissions for chemistry if separate SOURCE process is treated
    if (.not. LEMISDEP_INCHEM) EMISX(:,:) = 0._r8

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1


        !// Tropopause height
        LMTP = LMTROP(I,J)
        LMTPP1 = LMTROP(I,J) + 1


        !// Generating vertical arrays
        do L = LMTP, LPAR

          !// Box thickness (for IJ build from ETAA/B)
          P2 = ETAA(L+1) + ETAB(L+1) * P(I,J)
          P1 = ETAA(L) + ETAB(L) * P(I,J)

          !// Box center pressure
          PMIDL(L) = (P1 + P2) * 0.5_r8
          !// Box top pressure
          PUEDG(L) = P2

          !// Box volume [m3]
          DV(L) = BVOL(L,II,JJ,MP)

          !// Calculate air density [molec/cm^3]
          AIR_MOLEC(L) = BAIR_MOLEC(L,II,JJ,MP)

          !// Temperature
          TEMP(L) = BTEM(L,II,JJ)

          !// J-values
          JVALUES(:,L) = BJV(:,L,II,JJ,MP)

        end do !// do L = LMTP,LPAR

        !// Get PSC1, PSC2 and aerosol SAD
        call get_psc12_sad(i,j,PSC1L, PSC2L, SADL)

        !// Treat emissions as production term in chemistry?
        if (LEMISDEP_INCHEM) then
          !// Get emissions (includes all emissions; surface, lightning, etc.)
          !// EMISX is initialzed inside.
          call getEMISX(BEMIS(:,:,II,JJ),DV,I,J,II,JJ,LMTPP1,LPAR, &
               LPAR, MP, NPAR, TRACER_ID_MAX, chem_idx, TMASS, EMISX)
        end if

        !// Set local arrays (ZC_LOCAL) [molec/cm3]
        call gotoZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,LMTPP1,LPAR, &
             IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
             NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)


        !// ZC from mass to concentration
        call ZC_MASS2CONC(ZC_LOCAL,LMTPP1,LPAR,LPAR,NPAR,NOTRPAR, &
             TRACER_ID_MAX, trsp_idx, Xtrsp_idx,TMASS,XTMASS,DV,AVOGNR)


        !// H2O treatment
        call zc_strh2o(ZC_LOCAL,II,JJ,MP,LMTPP1,LPAR, &
             LPAR,TRACER_ID_MAX,trsp_idx,DV)


        !// Calculate rates dependent on T & P
        call TCRATE_TP_IJ_STR(LPAR,  TRACER_ID_MAX, LMTP, TEMP, &
             AIR_MOLEC,PMIDL, &
             r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
             r_oh_hno3, r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
             r_cl_o2_m, r_cloo_heat, r_clo_clo_m, r_cl2o2_m, r_oh_no2_m, &
             r_op_no2_m, r_oh_oh_m, r_oh_no_m, r_no2_clo_m,  &
             r_op_o2_m, r_o2_h_m, r_no2_bro_m, &
             r_no_ho2_b )

        !// Rates for heterogeneous chemistry
        call TCRATE_HET_IJ(LPAR,  TRACER_ID_MAX, LMTP, TEMP, &
                           ZC_LOCAL, &
                           SADL, PSC1L, PSC2L, &
                           spsADPARBK, spsACPARBK, spsBDPARBK, &
                           spsBCPARBK, spsECPARBK, &
                           spsFDPARBK, spsFCPARBK, spsGCPARBK, &
                           spsADPAR, spsACPAR, spsBDPAR, &
                           spsBCPAR, spsECPAR, &
                           spsFDPAR, spsFCPAR, spsGCPAR, spsBHPAR )

          
        call OSLO_CHEM_STR_IJ ( &
             !//rate constants dependent on more than temperature
             r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
             r_oh_hno3, r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
             r_cl_o2_m, r_cloo_heat, r_clo_clo_m, r_cl2o2_m, r_oh_no2_m, &
             r_op_no2_m, r_oh_oh_m, r_oh_no_m, r_no2_clo_m,  &
             r_op_o2_m, r_o2_h_m, r_no2_bro_m, &
             r_no_ho2_b, &
             !// other parameters
             LMTP, ZC_LOCAL, JVALUES, &
             TEMP, AIR_MOLEC, &
             LPAR, TEMPRANGE, MINTEMP, trsp_idx, Xtrsp_idx, TRACER_ID_MAX, &
             JVN_, I, J, &
             DTCH_STR, ST_STR, QLIN_STR, QTEST_STR, NCHEM_ITER, &
             LOLD_H2OTREATMENT, &
             !// emissions
             EMISX, &
             !// deltaO3
             deltaO3, &
             !// heterogeneous chemistry rate constants
             spsADPARBK, spsADPAR, &
             spsACPARBK, spsACPAR, &
             spsBDPARBK, spsBDPAR, &
             spsBCPARBK, spsBCPAR, &
             spsECPARBK, spsECPAR, &
             spsFDPARBK, spsFDPAR, &
             spsFCPARBK, spsFCPAR, &
             spsGCPARBK, spsGCPAR, &
                         spsBHPAR, &
             nchemdiag,CHEMLOSS,CHEMPROD, &
             OxCHEMPROD, OxCHEMLOSS )


        !// Diagnose chemistry losses and prods
        call save_chemPL(CHEMLOSS,CHEMPROD,DV,PUEDG,I,J,II,JJ,MP, &
                         LMTPP1,LPAR,LMTP)
        call save_chemOxPL(OxCHEMLOSS,OxCHEMPROD,DV,I,J,LMTPP1,LPAR)

        !// ZC from concentration to mass
        call ZC_CONC2MASS(ZC_LOCAL,LMTPP1,LPAR,LPAR,NPAR,NOTRPAR, &
             TRACER_ID_MAX, trsp_idx, Xtrsp_idx,TMASS,XTMASS,DV,AVOGNR)

        if (LSOA) then
          !// Do SOA separation in stratosphere
          ZC0(:,:) = ZC_LOCAL(:,:)
          call SOA_v9_separate(ZC_LOCAL,TEMP,DV,I,J,LMTPP1,LPAR)
          !// Accumulate SOA production (from level LMTP+1 to LPAR)
          ZC0(:,:) = ZC_LOCAL(:,:) - ZC0(:,:)
          call soa_diag_separate(ZC0,I,J,LMTPP1,LPAR)
        end if !// if (LSOA) then

        !// Move local arrays back to BTT / XSTT
        call backfromZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,LMTPP1, &
             LPAR,IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
             NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)


        !// Adjust diagnose due to strat O3 treatment; molec/cm3 to kg/gridbox
        !// Only use stratospheric part where DV is set!
        if (trsp_idx(1) .gt. 0) then 
          deltaO3(LMTPP1:LPAR) = deltaO3(LMTPP1:LPAR) &
                * 1.e3_r8 * TMASS(trsp_idx(1)) / AVOGNR * DV(LMTPP1:LPAR)
          BTTBCK(LMTPP1:LPAR,trsp_idx(1),II,JJ) = &
               BTTBCK(LMTPP1:LPAR,trsp_idx(1),II,JJ) + deltaO3(LMTPP1:LPAR)
        end if

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)

    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine oslochem_strat
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine read_oslo2d2(LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Read boundary data from the Oslo 2D model.
    !// Updated routine reading original Oslo 2D data (sr-files) and
    !// interpolates to model resolution on the fly.
    !//
    !// Ole Amund Sovde, December 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPARW, JPAR
    use cmn_ctm, only: JMON, YDGRD
    use cmn_met, only: MYEAR
    use utilities, only: get_free_fileid
    use cmn_oslo, only: STT_2D_LB, STT_2D_LT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    logical, intent(in) :: LNEW_MONTH

    !// File name
    character(len=80) :: filename

    !// String year
    character(len=4) :: YEAR
    !// String month
    character(len=2) :: MON
    !// Counters:
    integer :: J, L, N, M, ifnr, ierr, NYR,IDAT,MND, Counter
    integer :: LEV2D
    real(r8) :: A
    !// Original data from Oslo 2D
    integer, parameter :: &
         J2D = 19, &        !// Number of latitudes
         J2D_edg = J2D+1, & !// Number of edges
         L2D = 26, &        !// Number of levels
         N2DL= 32, &        !// Longlived components
         N2DS= 31, &        !// Shortlived components
         N2D = N2DL+N2DS    !// Total components on 2D file
    real(r8) :: &
         XLONG(J2D,L2D,N2DL), &     !// Longlived components read from file
         XSHORT(J2D,L2D-1,N2DS), &  !// Shortlived components read from file
         AM(J2D,L2D,12)    !// Densities from OSLO 2-D model (one field per
    !// All species plus one (NOy)
    real(r8) :: STT2D(J2D,L2D,N2D+1)

    !// For interpolation
    real(r8) :: r8in(J2D), r8out(JPAR)
    !// For scaling (L60)
    real(r8) :: zfact, rscale(N2D),scalez
    !// Grid centers for Oslo 2D
    real(r8), dimension(J2D),parameter :: TWODGRID = &
          (/ -90._r8, -80._r8, -70._r8, -60._r8, -50._r8, -40._r8, -30._r8, &
             -20._r8, -10._r8,   0._r8,  10._r8,  20._r8,  30._r8,  40._r8, &
              50._r8,  60._r8,  70._r8,  80._r8,  90._r8 /)

    !// There are 64 components.
    !// 1) Components that are treated in the stratosphere AND NOT
    !//   in the troposphere (according to SCTM-1/OSLO-2D numbering):
    !// CTM3 number is STRAT_ONLY(N)+100, i.e. 101, 102, ..., 147
    integer, parameter :: NSTRAT = 47, NTOT=64
    integer, dimension(NSTRAT),parameter :: STRAT_ONLY= &
          (/  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
             14, 15, 16, 17, 18, 19, 20, 21, 22, 23, &
             24, 25, 26, 27, 28, 29, 30, 31, 32, 36, &
             44, 46, 47, 48, 49, 53, 54, 55, 56, 57, &
             58, 59, 60, 61, 62, 63, 64 /)

    !//  2) Components that are treated in the stratosphere (according
    !//     to SCTM-1/OSLO-2D numbering). The unusual order derives from
    !//     the order of components in CTM-2.
    !// This gives the order of STT_2D_LT
    integer, dimension(NTOT),parameter :: SCTM_NO= &       ! N on file
         (/ 35, 11, 12, 50, 39, 52, 45, 38, 51, 33, &      ! 01-10
            34, 37, 42, 43, 40, 41, 13, &                  ! 11-17
             1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &      ! 18-27
            14, 15, 16, 17, 18, 19, 20, 21, 22, 23, &      ! 28-37
            24, 25, 26, 27, 28, 29, 30, 31, 32, 36, &      ! 38-47
            44, 46, 47, 48, 49, 53, 54, 55, 56, 57, &      ! 48-57
            58, 59, 60, 61, 62, 63, 64 /)                  ! 58-64
    !// Their names (in file order) in CTM3 are 
    integer, dimension(NTOT),parameter :: CTM3NO = &       ! N on file
               (/  1,          4,      6,                & ! 01-03
                          13,     15, 16, 17,            & ! 04-07
                  21, 22,                                & ! 08-09
                                              38, 39, 40,& ! 10-12
                  41, 42, 43, 44,     46,                & ! 13-17
                  !// Components treated in stratospheric chemistry only:
                     101,102,103,104,105,106,107,108,109,& ! 18-26
                 110,111,112,113,114,115,116,117,118,119,& ! 27-36
                 120,121,122,123,124,125,126,127,128,129,& ! 37-46
                 130,131,132,133,134,135,136,137,138,139,& ! 47-56
                 140,141,142,143,144,145,146,147 /)        ! 57-64
    !//real(r8) :: STT_2D_LB2(JPAR,10,47), STT_2D_LT2(JPAR,64)
    !// In other words, the 2D data are
    !// SCTM  Name    CTM3  Place in STT_2D_LT  Place in STT_2D_BT
    !//   1   MCF      101    18                  CTM3-100 = 1
    !//   2   HCFC22   102    19                  CTM3-100 = 2
    !//   3  CFC11     103    20                    :
    !//   4  CFC12     104    21                    :
    !//   5  CCl4      105    22
    !//   6  CH3Cl     106    23
    !//   7  N2O       107    24
    !//   8  Clx       108    25
    !//   9  NOx_str   109    26
    !//  10  SO        110    27
    !//  11  HNO3        4    2
    !//  12  CO          6    3
    !//  13  CH4        46    17
    !//  14  HCl       111    28
    !//  15  Cly       112    29
    !//  16  H2        113    30
    !//  17  H2O       114    31
    !//  18  SH        115    32
    !//  19  CH3Br     116    33
    !//  20  H1211     117    34
    !//  21  H1301     118    35
    !//  22  Bry       119    36
    !//  23  H2402     120    37
    !//  24  F113      121    38
    !//  25  F114      122    39
    !//  26  F115      123    40
    !//  27  HNO3s     124    41
    !//  28  H2Os      125    42
    !//  29  HCls      126    43
    !//  30  HCFC123   127    44
    !//  31  HCFC141   128    45
    !//  32  HCFC142   129    46
    !//  33  OP         38    10
    !//  34  OD         39    11
    !//  35  O3          1     1
    !//  36  H         130    47
    !//  37  OH         40    12
    !//  38  HO2        21     8
    !//  39  H2O2       15     5
    !//  40  NO         43    15
    !//  41  NO2        44    16
    !//  42  NO3        41    13
    !//  43  N2O5       42    14
    !//  44  HNO2      131    48
    !//  45  HO2NO2     17     7
    !//  46  Cl        132    49
    !//  47  ClO       133    50
    !//  48  OHCl      134    51
    !//  49  ClONO2    135    52                  CTM3-100 = 35
    !//  50  HCHO       13     4
    !//  51  CH3O2      22     9
    !//  52  CH3O2H     16     6
    !//  53  Cl2       136    53                  CTM3-100 = 36
    !//  54  OClO      137    54
    !//  55  Br        138    55
    !//  56  BrO       139    56
    !//  57  HBr       140    57
    !//  58  BrONO2    141    58
    !//  59  OHBr      142    59
    !//  60  Br2       143    60
    !//  61  ClOO      144    61
    !//  62  Cl2O2     145    62
    !//  63  BrCl      146    63
    !// ------------------------------------------------------------------
    !// What is done in this routine:
    !// 1. Read 2D data
    !// 2. Interpolate
    !// 3. Put into STT_2d_* arrays
    !// 4. Scale and re-make sums if necessary
    !// ------------------------------------------------------------------

    !// Only do this every month
    if (.not. LNEW_MONTH) return

    !// Year and month as strings
    write(YEAR,'(I4)') MYEAR
    write(MON,'(I2.2)') JMON

    write(6,'(a)') '* Updating boundary conditions'

    !// Find file number to use
    ifnr = get_free_fileid()


    !// Read air density from Oslo 2D (AM = air density field, molec/cm3)
    !// Needed to convert short lived species to mixing ratio
    !// ------------------------------------------------------------------
    filename='Indata_CTM3/2d_data/am.dat'
    open(ifnr,file=filename,form='formatted',status='OLD',iostat=ierr)
    if (ierr.ne.0) then
       print*,'* No such file: '//trim(filename)
       stop
    end if
    write(6,'(a)') '  Reading '//trim(filename)
    do N = 12, 21, 3              ! season (in am.dat DJF comes first) ...
      MND = N
      if (MND.gt.12) MND = MND - 12
      do J = 1, J2D
        !// Stored from 80N to 80S; do not change this
        read(ifnr,'(8E10.4)') (AM(J,L,MND),L= 1, 8)
        read(ifnr,'(8E10.4)') (AM(J,L,MND),L= 9,16)
        read(ifnr,'(8E10.4)') (AM(J,L,MND),L=17,24)
        read(ifnr,'(8E10.4)') (AM(J,L,MND),L=25,26)
        do L = 1, 26  ! AM(J,L,1) for January, AM(J,L,12) for December
          if (MND.eq.12) then
            AM(J,L,1) = AM(J,L,MND)
            AM(J,L,2) = AM(J,L,MND)
          else
            AM(J,L,MND+1) = AM(J,L,MND)
            AM(J,L,MND+2) = AM(J,L,MND)
          end if
        end do
      end do
    end do
    close(ifnr)


    !// Read 2D data for current year
    !// ------------------------------------------------------------------
    !// File
    filename='Indata_CTM3/2d_data/sr'//YEAR
    open(ifnr,FILE=filename,form='formatted',STATUS='OLD',iostat=ierr)
    if (ierr.ne.0) then
       print*,'*** No 2D file in stratchem_oslo.f90:'//trim(filename)
       stop
    end if
    write(6,'(a)') '  Reading '//trim(filename)

    !// Read to get this month
    do M = 1, JMON
      read(ifnr,'((1X,6X,6X,I4,6X,i4))') NYR,IDAT
      do N = 1, N2DL
        do J = 1, J2D
          read(ifnr,'(8E10.4)') (XLONG(J,L,N),L= 1, 8)
          read(ifnr,'(8E10.4)') (XLONG(J,L,N),L= 9,16)
          read(ifnr,'(8E10.4)') (XLONG(J,L,N),L=17,24)
          read(ifnr,'(8E10.4)') (XLONG(J,L,N),L=25,26)
        end do
      end do
      do J = 2, J2D-1
        do L = 1, L2D-1
          read(ifnr,'(8E10.4)') (XSHORT(J,L,N),N= 1, 8)
          read(ifnr,'(8E10.4)') (XSHORT(J,L,N),N= 9,16)
          read(ifnr,'(8E10.4)') (XSHORT(J,L,N),N=17,24)
          read(ifnr,'(8E10.4)') (XSHORT(J,L,N),N=25,31)
        end do
      end do
    end do !// do M = 1, JMON
    close(ifnr)


    !// Check data
    do N = 1, N2DL
      do L = 1, L2D-1
        do J = 2, J2D-1
          if (XLONG(J,L,N) .LT. 0._r8) then
            A = XLONG(J,L-1,N) * XLONG(J,L+1,N)
            if (A .ge. 0._r8) XLONG(J,L,N) = sqrt(A)
            if (A .lt. 0._r8) XLONG(J,L,N) = &
                 (XLONG(J-1,L,N) + XLONG(J+1,L,N)) * 0.5_r8
          end if
          if (XLONG(J,L,N) .lt. 0._r8) XLONG(J,L,N) = &
               (XLONG(J,L-1,N) + XLONG(J,L+1,N)) * 0.5_r8
        end do
      end do
    end do

    !// sometimes layer 26 gets ridiculously small numbers. Fix this by
    !// checking 1 to 25 only and setting layer 26 to zero. in this layer
    !// Oslo 2d should give only zeros anyway.
    do N = 1, N2DL
      do L = 1, L2D-1
        do J = 2, J2D-1
          if (XLONG(J,L,N) .lt. 1.e-35_r8) XLONG(J,L,N) = 0._r8
        end do
      end do
      do L = L2D, L2D
        do J = 2, J2D-1
          XLONG(J,L,N) = 0._r8
        end do
      end do
    end do

    !// Put longlived into STT2D, switch from N-S to S-N
    !// ------------------------------------------------------------------
    do N = 1, N2DL
      do L = 1, L2D-1
        do J = 2, J2D-1 !// Poles treated afterwards
           STT2D(J2D-J+1,L,N) = XLONG(J,L,N)
        end do
      end do
    end do

    !// Check short lived components
    do N = 1, N2DS
      do L = 1, L2D-1
        do J = 2, J2D-1
         if (XSHORT(J,L,N) .lt. 1.e-35_r8) XSHORT(J,L,N) = 0._r8
        end do
      end do
    end do

    !// Put shortlived into STT2D, switch from N-S to S-N
    !// Also convert from concentration to mixing ratio (short lived
    !// species are given in concentration)
    !// ------------------------------------------------------------------
    do N = 1, N2DS
      do L = 1, L2D-1
        do J = 2, J2D-1 !// Poles treated afterwards
           STT2D(J2D-J+1,L,N2DL+N) = XSHORT(J,L,N) / AM(J,L,JMON)
        end do
      end do
    end do

    !// Set poles equal to neighbor latitude
    STT2D(1,:,:)   = STT2D(2,:,:)
    STT2D(J2D,:,:) = STT2D(J2D-1,:,:)


    !// Need some extra species, i.e. sums
    do L = 1, L2D
      do J = 1, J2D
        STT2D(J,L,15) = STT2D(J,L,8) &
                        +STT2D(J,L,14) ! Cly=Clxx+HCl
        STT2D(J,L,18) = STT2D(J,L,36) &
                        +STT2D(J,L,37) &
                        +STT2D(J,L,38) &
                        +STT2D(J,L,39) ! SH=H+OH+HO2+H2O2
        STT2D(J,L,53) = 1.e-18_r8  ! Set Cl=0.
        STT2D(J,L,64) = STT2D(J,L,9) &
                        +STT2D(J,L,11) ! NOy=NOx+HNO3
      end do
    end do


    !// Interpolation; keep old for simplicity (gives smoother curves)
    !// ------------------------------------------------------------------
    !// Lower boundary; i.e. STRAT_ONLY
    !// ------------------------------------------------------------------
    !// Interpolate for strat_only (10 layers starting from surface)
    do N = 1, 47
      do L = 1, 10
        !// Interpolate in N-S; Keep old interpolation
        r8in(:) = STT2D(:,L,STRAT_ONLY(N))
        do J = 1, JPAR
          Counter = 1
          do while (TwoDGrid(Counter) .le. YDGRD(J))
            Counter = Counter + 1
            if (counter .ge. J2D) exit
          end do
          r8out(J) = r8in(Counter-1) &
              * (TwoDGrid(Counter) - YDGRD(J)) &
              / (TwoDGrid(Counter) - TwoDGrid(Counter-1)) &
             + r8in(Counter) &
              * (YDGRD(J) - TwoDGrid(Counter-1)) &
              / (TwoDGrid(Counter) - TwoDGrid(Counter-1))
          !// Save interpolated value (set minimum to very small, so that
          !// mass is nonzero)
          STT_2D_LB(J,L,N) = max(1.e-18_r8,r8out(J))
        end do
      end do !// do L = 1, 10
      !// print*,'LB: STRAT ONLY',strat_only(N),'into',N+100
    end do !// do N = 1, 47
    !// Lower boundary ok.
    !// ------------------------------------------------------------------



    !// Upper boundary layer
    !// ------------------------------------------------------------------
    if (LPARW.eq.40) then
      !// For L40, pick 2D layer 16
      LEV2D = 16
    else if(LPARW.eq.60) then
      !// Use next to uppermost level (i.e. 25)
      LEV2D = L2D - 1
    else
      print*,'Wrong resolution in oc_read_oslo2d2',LPARW
      stop
    end if

    do N = 1, NTOT !// Include NOy
      !// Pick component from SCTM-list
      !// 1: O3; SCTM:35, CTM:1
      !// 2: HNO3: SCTM:11, CTM:4
      !// etc.
      r8in(:) = STT2D(:,LEV2D,SCTM_NO(N))
      do J = 1, JPAR
        Counter = 1
        do while (TwoDGrid(Counter) .le. YDGRD(J))
          Counter = Counter + 1
          if (counter .ge. J2D) exit
        end do
        r8out(J) = r8in(Counter-1) &
              * (TwoDGrid(Counter) - YDGRD(J)) &
              / (TwoDGrid(Counter) - TwoDGrid(Counter-1)) &
             + r8in(Counter) &
              * (YDGRD(J) - TwoDGrid(Counter-1)) &
              / (TwoDGrid(Counter) - TwoDGrid(Counter-1))
         !// Save interpolated value; scale later
        r8out(J) = max(1.e-18_r8,r8out(J))
        STT_2D_LT(J,N) = r8out(J)
      end do
    end do !// do N = 1, NTOT

    !// Upper boundary ok.
    !// ------------------------------------------------------------------
      
    !// Done with interpolations
    !// For Lyy>L40 we need to scale according to scale height
    !// Also need to remake sums (Cly, SH, Cl2, NOy, and HCl)
    !// ------------------------------------------------------------------
    if (LPARW .eq. 60) then
 
      !// Save Cly before scaling; assume constant in the scaled layer
      r8out(:) = STT_2D_LT(:,18) * 3._r8  & !MCF
                 + STT_2D_LT(:,19) * 1._r8 & !HCFC-22
                 + STT_2D_LT(:,20) * 3._r8 & !CFC-11
                 + STT_2D_LT(:,21) * 2._r8 & !CFC-12 
                 + STT_2D_LT(:,22) * 4._r8 & !CCl4
                 + STT_2D_LT(:,23) * 1._r8 & !CH3Cl
                 + STT_2D_LT(:,28) * 1._r8 & !HCl
                 + STT_2D_LT(:,38) * 3._r8 & !CFC-113
                 + STT_2D_LT(:,39) * 2._r8 & !CFC-114
                 + STT_2D_LT(:,40) * 1._r8 & !CFC-115
                 + STT_2D_LT(:,44) * 2._r8 & !HCFC-123
                 + STT_2D_LT(:,45) * 2._r8 & !HCFC-141
                 + STT_2D_LT(:,46) * 1._r8 & !HCFC-142
                 + STT_2D_LT(:,25) * 1._r8   !Clx

      !// Now do scaling
      do N = 1, N2D !// Not NOy (64)
        !// L60 is above 2D-data, so we need to scale according to a scale
        !// height. 2D layer 26 contains no data, use scale heights:
        !// Use Markku's scale height approach
        !//
        !// N is index of CTM3NO, i.e. 1, 4, 6, ...
        !// Use SCTM_NO(N) to locate the same components (35, 11, 12, ...)
        !// Select uppermost layer and scaling; pick component SCTM_NO(N)
        !// Scale from gradient at the top
        J = 10 ! i.e. tropical scale height at equator is computed
        if (STT2D(J,LEV2D-1,SCTM_NO(N)).le.0._r8 .or. &
               STT2D(J,LEV2D,SCTM_NO(N)).le.0._r8) then
          !// Species has zero or negative amount at the upper layers;
          !// set scale to zero
          rscale(N) = 0._r8
        else if (STT2D(J,LEV2D-1,SCTM_NO(N)) .eq. &
                 STT2D(J,LEV2D,SCTM_NO(N))) then
          !// Species has constant and non-zero mixing ratio at the
          !// upper layers
          rscale(N) = 1._r8
        else
          !// X is already converted into mixing ratio at this stage
          !// zfact which takes into account the decrease of air density
          !// between 47 and 49 km is not needed!
          zfact = 1._r8
          scalez = -2._r8/log((STT2D(J,LEV2D,SCTM_NO(N)) / &
                              STT2D(J,LEV2D-1,SCTM_NO(N)))*zfact)
          !// Use scale height:
          !// Assume mixing ratio to decrease/increase exponentially with
          !// altitude, characterized by a scale height. Let's use 61 km
          !// as 'characteristic' for layer 60. Can argue that 65 km is
          !// more appropriate ...
          rscale(N) = exp(-(63._r8 - 49._r8) / scalez)
        end if


        !// Scale STT_2D_LT2 with rscale
        !// NOy is not scaled, but set below
        do J = 1, JPAR
           STT_2D_LT(J,N) = STT_2D_LT(J,N) * rscale(N)
        end do
      end do


      do J = 1, JPAR

        !// Cly comes in 29.place among all CTM species in CTM-2, use
        !// SCTM_NO to count, answering the question where is 15 (Cly)
        STT_2D_LT(J,32) = STT_2D_LT(J,47) &
                          + STT_2D_LT(J,12) &
                          + STT_2D_LT(J,8) &
                          + 2._r8 * STT_2D_LT(J,5) ! SH=H+OH+HO2+H2O2

        STT_2D_LT(J,53)=1.e-18_r8   ! Set Cl2=0.

        !// new component (NOy) that is not included in m2ddata*
        !// set NOx=9ppb between 60 and 70 km
        STT_2D_LT(J,64) = 9.e-9_r8 & ! =NOx between 60 and 70 km
                         + STT_2D_LT(J,2) ! NOy=NOx+HNO3

        !// calculate HCL as sum(source gases+CLx+HCl) at uppermost 
        !// chemistry layer minus sum(source gases+CLx) at current layer
        STT_2D_LT(J,28) = r8out(J) & !// Cly before scaling
                     - STT_2D_LT(J,18) * 3._r8 - STT_2D_LT(J,19) * 1._r8 &
                     - STT_2D_LT(J,20) * 3._r8 - STT_2D_LT(J,21) * 2._r8 &
                     - STT_2D_LT(J,22) * 4._r8 - STT_2D_LT(J,23) * 1._r8 &
                     - STT_2D_LT(J,38) * 3._r8 &
                     - STT_2D_LT(J,39) * 2._r8 - STT_2D_LT(J,40) * 1._r8 &
                     - STT_2D_LT(J,44) * 2._r8 - STT_2D_LT(J,45) * 2._r8 &
                     - STT_2D_LT(J,46) * 1._r8 - STT_2D_LT(J,25) * 1._r8
        if (STT_2D_LT(J,28).lt.0._r8) then
           print*,'oc_read_oslo2d2: ups',STT_2D_LT(J,28)
           stop
        end if

        !// Here for Cly
        STT_2D_LT(J,29) = STT_2D_LT(J,25) &
                          + STT_2D_LT(J,28) ! Cly=Clxx+HCl

      end do !// do J = 1, JPAR
    end if !// if (LPARW .eq. 60) then

    write(6,'(a,i2)') '  Updated STT_2D_* for month: ',JMON
    !// --------------------------------------------------------------------
  end subroutine read_oslo2d2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_strat_boundaries(BTT,BTTBCK,AIRB,BTEM,DTADV,MP)
    !// --------------------------------------------------------------------
    !// Update BTT and XSTT upper and lower boundary using the Oslo 2D data.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    !// Update STT now:
    !// At the top level an upper boundary condition is applied to all
    !// components used in stratospheric chemistry (note that some of these
    !// components are treated in tropospheric chemistry, too). Mixing
    !// ratios are retrieved from OSLO 2D model (31km level !) and stored in
    !// STT_2D_LT(JCPAR,64) array (JCPAR longitudes, 64 components).
    !//
    !// At the bottom level a lower boundary condition is applied to all
    !// components used in stratospheric chemistry AND NOT in tropospheric
    !// chemistry. Mixing ratios are retrieved from OSLO 2D model and stored
    !// in STT_2D_LB(JCPAR,10,47) array (JCPAR longitudes, 10 layers, 47
    !// components). The 10 layers of archived Oslo 2D data are:
    !// 0-2, 2-4, 4-6, 6-8, 8-10, 10-12, 12-14, 14-16, 16-18, and 18-20 km.
    !//
    !// Exceptions: 114 (H2O), 124 (HNO3s), 125 (H2Os), and 126 (HCls) are
    !//           not treated here. H2O is represented by H2O_molec, while
    !//           the solid compounds are calculated in the heterogeneous
    !//           chemistry modules.
    !//
    !// Last modified: 08Mar2003, mga
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IDBLK, JDBLK, LPAR, LPARW, NPAR
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TMASS, TMASSMIX2MOLMIX, TMOLMIX2MASSMIX
    use cmn_oslo, only: LMTROP, trsp_idx, Xtrsp_idx, &
         XTMASSMIX2MOLMIX, XTMOLMIX2MASSMIX, XTMASS, &
         STT_2D_LB, STT_2D_LT, XSTT, METHANEMIS
    use strat_h2o, only: set_trop_h2o_b4trsp_clim, &
         strat_h2o_ubc2
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTTBCK
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in)  :: AIRB, BTEM
    real(r8), intent(in)  :: DTADV
    integer, intent(in) :: MP

    !// Locals
    !// components treated in strat. chemistry according to CTM2 numbering
    integer, dimension(64),parameter :: CTM2NO = &
         !// The following list contains all components that are included in
         !// stratospheric chemistry, with respect to CTM2 numbering
         !// (CTM2NO= CTM2-number). Components treated in both tropospheric
         !// and stratospheric chemistry: O3, HNO3, CO, CH2O, H2O2, etc.
         (/  1,          4,      6,                 & ! 01-03
                    13,     15, 16, 17,             & ! 04-07
            21, 22,                                 & ! 08-09
                                        38, 39, 40, & ! 10-12
            41, 42, 43, 44,     46,                 & ! 13-17
         !// Components treated in stratospheric chemistry only:
         !// MCF, HCFC22, CXFC11, CFC12, CCl4, etc.
                101,102,103,104,105,106,107,108,109, & ! 18-26
            110,111,112,113,114,115,116,117,118,119, & ! 27-36
            120,121,122,123,124,125,126,127,128,129, & ! 37-46
            130,131,132,133,134,135,136,137,138,139, & ! 47-56
            140,141,142,143,144,145,146,147 /)         ! 57-64

    integer :: I,J,L,N,MTC_COMP,XMTC_COMP,II,JJ,MTC_LOC,XMTC_LOC
    !// Surface conditions set up to layer LSFC
    integer,parameter :: LSFC=1

    !// List of components for families
    integer,parameter :: NR_FAM=7
    integer :: FAM_ID(NR_FAM)
    integer,parameter :: maxmembers=9
    integer :: FAM_MEM(maxmembers,NR_FAM),FAM_LOC
    real(r8)  :: MEM_FAC(maxmembers,NR_FAM)

    real(r8) :: RFRAC,OZONE
    !// --------------------------------------------------------------------

    !// Upper and lower boundary conditions for H2O is set at the
    !// bottom of this routine.
 
    !// If stratospheric chemistry is done:
    !// --------------------------------------------------------------------

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1


        !// First group: components that are used in stratospheric chemistry
        !// --------------------------------------------------------------
        !// Apply upper boundary condition:
        !// --------------------------------------------------------------
        !// (STT_2D_LT contains this group of components in correct order)

        do N = 1, 64

          !// Skip for H2O
          if (CTM2NO(N).eq.114) cycle
          !// Skip for H2 if transported
          if (trsp_idx(113).gt.0 .and. CTM2NO(N).eq.113) cycle
          !// Skip for CH4 if emissions are included
          if (METHANEMIS .and. CTM2NO(N).eq.46) cycle

          !// Transported species first
          MTC_COMP = trsp_idx(CTM2NO(N))
          if (MTC_COMP.gt.0) then
            BTT(LPAR,MTC_COMP,II,JJ) = STT_2D_LT(J,N) &
                 !* AIRB(LPAR,II,JJ) / TMASSMIX2MOLMIX(MTC_COMP)
                 * AIRB(LPAR,II,JJ) * TMOLMIX2MASSMIX(MTC_COMP)
           !// Reset for diagnostics
            BTTBCK(LPAR,MTC_COMP,II,JJ) = BTT(LPAR,MTC_COMP,II,JJ)
          end if
          !// Non-transported species
          XMTC_COMP = Xtrsp_idx(CTM2NO(N))
          if (XMTC_COMP.gt.0) then
            XSTT(LPAR,XMTC_COMP,I,J) = STT_2D_LT(J,N) &
                  !* AIRB(LPAR,II,JJ) / XTMASSMIX2MOLMIX(XMTC_COMP)
                  * AIRB(LPAR,II,JJ) * XTMOLMIX2MASSMIX(XMTC_COMP)
          end if

        end do

        !// Second group: Components included in stratospheric chemistry
        !// --------------------------------------------------------------
        !// Apply lower boundary condition:
        !// --------------------------------------------------------------
        !// AND NOT in tropospheric chemistry) need a lower boundary condition.
        !// Their numbers with respect to CTM2 numbering are simply 101, 102,
        !// ..., 147 (MCF, HCFC22, CXFC11, CFC12, CCl4, etc.
        !// - i.e. 47 components).
        !// Note that all these components are initialized throughout the
        !// troposphere at the start of the model run.
        do N = 101, 147
          !// No need to set H2O in lower boundary conditions; set from q.
          if (N.eq.114) cycle

          !// Get transport number
          MTC_COMP = trsp_idx(N)
          if (MTC_COMP.gt.0) then
            !// No need to set surface higher than layer 1.
            !// Setting it too high can cause problems in transport if the
            !// values are zero! It can then becom inconsistent with the
            !// moments and give negative values after transport.
            do L=1,LSFC !15       !// 0-2 km
              BTT(L,MTC_COMP,II,JJ) = STT_2D_LB(J,1,N-100) &
                   !* AIRB(L,II,JJ) / TMASSMIX2MOLMIX(MTC_COMP)
                   * AIRB(L,II,JJ) * TMOLMIX2MASSMIX(MTC_COMP)
              !// Reset for diagnostics
              BTTBCK(L,MTC_COMP,II,JJ) = BTT(L,MTC_COMP,II,JJ)
            end do
          end if
          !// Non-transported species (not needed, they are not transported)
          XMTC_COMP = Xtrsp_idx(N)
          if (XMTC_COMP.gt.0) then
            do L=1,LSFC       !// 0 - 2 km
              XSTT(L,XMTC_COMP,I,J) = STT_2D_LB(J,1,N-100) &
                   !*AIRB(L,II,JJ) / XTMASSMIX2MOLMIX(XMTC_COMP)
                   *AIRB(L,II,JJ) * XTMOLMIX2MASSMIX(XMTC_COMP)
            end do
          end if
        end do

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// Also set H2O-related lower and upper boundary conditions
    !// --------------------------------------------------------------------
    !// Update H2O before transport (strat_h2o.f90)
    call set_trop_h2o_b4trsp_clim(BTT,BTTBCK,AIRB,BTEM,MP)

    !// Upper boundary conditions for H2O AND H2 (strat_h2o.f90)
    call strat_h2o_ubc2(BTT,BTTBCK,AIRB,MP)

    !// --------------------------------------------------------------------
  end subroutine update_strat_boundaries
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_fam_in_trop(BTT,BTTBCK,MP)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IDBLK, JDBLK, LPAR, NPAR
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE
    use cmn_chem, only: TMASS
    use cmn_oslo, only: LMTROP, trsp_idx, Xtrsp_idx, XTMASS, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTTBCK
    integer, intent(in) :: MP

    integer :: I,J,L,II,JJ,F,M,MTC_LOC,XMTC_LOC

    !// List of components for families
    integer, parameter :: NR_FAM = 7
    integer :: FAM_ID(NR_FAM)
    integer, parameter :: maxmembers = 9
    integer :: FAM_MEM(maxmembers,NR_FAM),FAM_LOC
    real(r8)  :: MEM_FAC(maxmembers,NR_FAM)
    !// --------------------------------------------------------------------
    !// Families are now calculated from their members throughout the
    !// troposphere:
    !// --------------------------------------------------------------------
    !//
    !// Clx =  Cl + ClO + OHCl + ClONO2 + 2*Cl2 + OClO + BrCl + ClOO + 2*Cl2O2
    !// 108 = 132 + 133 +  134 +    135 + 2*136 +  137 +  146 +  144 +   2*145
    !//
    !// NOx_str = NO + NO2 + NO3 + 2*N2O5 + HO2NO2 + ClONO2 + BrONO2
    !//     109 = 43 +  44 +  41 +   2*42 +     17 +    135 +    141
    !//
    !// SO  = O3 + O(3P) + O(1D) - NO -  Cl -  Br
    !// 110 =  1 +    38 +    39 - 43 - 132 - 138
    !//
    !// Cly =  Cl + ClO + OHCl + ClONO2 + 2*Cl2 + OClO + BrCl + ClOO + 2*Cl2O2
    !//       + HCl + HCls = Clx + HCl + HCls
    !// 112 = 132 + 133 +  134 +    135 + 2*136 +  137 +  146 +  144 +   2*145
    !//       + 111 +  126 = 108 + 111 +  126
    !//
    !// SH  =   H + OH + HO2 + 2*H2O2
    !// 115 = 130 + 40 +  21 +   2*15
    !//
    !// Bry =  Br + BrO + BrONO2 + OHBr + HBr + 2*Br2 + BrCl
    !// 119 = 138 + 139 +    141 +  142 + 140 + 2*143 +  146
    !//
    !// NOy_str = NO + NO2 + NO3 + 2*N2O5 + HO2NO2 + ClONO2 + BrONO2
    !//          + HNO3 + HNO3s = NOx_str + HNO3 + HNO3s 
    !//     147 = 43 +  44 +  41 +   2*42 +     17 +    135 +    141
    !//          +    4 +   124 =     109 +    4 +   124
    !// --------------------------------------------------------------------

    !// FAM_MEM lists the members
    !// MEM_FAC lists the factors to account for the correct amount of
    !//         species. E.g. in Cly the factor for Cl2 is 2, but the
    !//         factor for HCl is 1.

    FAM_ID = (/108, 109, 110, 112, 115, 119, 147/)
    !// Clx (108)
    FAM_MEM(1: 6,1) = (/ 132,  133,  134,  135,  136,  137 /)
    FAM_MEM(7: 9,1) = (/ 146,  144,  145 /)
    MEM_FAC(1: 6,1) = (/1._r8, 1._r8, 1._r8, 1._r8, 2._r8, 1._r8 /)
    MEM_FAC(7: 9,1) = (/1._r8, 1._r8, 2._r8 /)
    !// NOx_str (109)
    FAM_MEM(1: 6,2) = (/  43,   44,   41,   42,   17,  135 /)
    FAM_MEM(7: 9,2) = (/ 141,    0,    0 /)
    MEM_FAC(1: 6,2) = (/1._r8, 1._r8, 1._r8, 2._r8, 1._r8, 1._r8 /)
    MEM_FAC(7: 9,2) = (/1._r8, 0._r8, 0._r8 /)
    !// SO (110)
    FAM_MEM(1: 6,3) = (/   1,   38,   39,   43,  132,  138 /)
    FAM_MEM(7: 9,3) = (/   0,    0,    0 /)
    MEM_FAC(1: 6,3) = (/1._r8, 1._r8, 1._r8,-1._r8,-1._r8,-1._r8 /)
    MEM_FAC(7: 9,3) = (/0._r8, 0._r8, 0._r8 /)
    !// Cly (112)
    FAM_MEM(1: 6,4) = (/ 108,  111,    0,    0,    0,    0 /)
    FAM_MEM(7: 9,4) = (/   0,    0,    0 /)
    MEM_FAC(1: 6,4) = (/1._r8, 1._r8, 0._r8, 0._r8, 0._r8, 0._r8 /)
    MEM_FAC(7: 9,4) = (/0._r8, 0._r8, 0._r8 /)
    !// SH (115)
    FAM_MEM(1: 6,5) = (/ 130,   40,   21,   15,    0,    0 /)
    FAM_MEM(7: 9,5) = (/   0,    0,    0 /)
    MEM_FAC(1: 6,5) = (/1._r8, 1._r8, 1._r8, 2._r8, 0._r8, 0._r8 /)
    MEM_FAC(7: 9,5) = (/0._r8, 0._r8, 0._r8 /)
    !// Bry (119)
    FAM_MEM(1: 6,6) = (/ 138,  139,  141,  142,  140,  143 /)
    FAM_MEM(7: 9,6) = (/ 146,    0,    0 /)
    MEM_FAC(1: 6,6) = (/1._r8, 1._r8, 1._r8, 1._r8, 1._r8, 2._r8 /)
    MEM_FAC(7: 9,6) = (/1._r8, 0._r8, 0._r8 /)
    !// NOy_str (147)
    FAM_MEM(1: 6,7) = (/ 109,    4,    0,    0,    0,    0 /)
    FAM_MEM(7: 9,7) = (/   0,    0,    0 /)
    MEM_FAC(1: 6,7) = (/1._r8, 1._r8, 0._r8, 0._r8, 0._r8, 0._r8 /)
    MEM_FAC(7: 9,7) = (/0._r8, 0._r8, 0._r8 /)


    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1
      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Loop over families
        do F=1,7

          !// Family must be transported
          FAM_LOC = trsp_idx(FAM_ID(F))
          if (FAM_LOC .le. 0) then
            print*,'Stratospheric families must be transported'
            stop
          end if

          do L=1, LMTROP(I,J)
            BTT(L,FAM_LOC,II,JJ) = 0._r8 !// Initialize tropospheric column
          end do

          !// Add up the members
          do M = 1, maxmembers
            if (FAM_MEM(M,F) .le. 0) exit !// Done with this family

            MTC_LOC  = trsp_idx(FAM_MEM(M,F))
            XMTC_LOC = Xtrsp_idx(FAM_MEM(M,F))
            if (MTC_LOC .gt. 0) then
              do L=1,LMTROP(I,J)
                !// Sum up members
                BTT(L,FAM_LOC,II,JJ) = BTT(L,FAM_LOC,II,JJ) &
                     + BTT(L,MTC_LOC,II,JJ) / TMASS(MTC_LOC) &
                       * MEM_FAC(M,F)
              end do
            else if (XMTC_LOC .gt. 0) then
              do L=1,LMTROP(I,J)
                !// Sum up members
                BTT(L,FAM_LOC,II,JJ) = BTT(L,FAM_LOC,II,JJ) &
                       + XSTT(L,XMTC_LOC,I,J) / XTMASS(XMTC_LOC) &
                       * MEM_FAC(M,F)
              end do
            end if
          end do !// do M = 1, maxmembers
          !// Scale it
          do L = 1, LMTROP(I,J)
            BTT(L,FAM_LOC,II,JJ) = BTT(L,FAM_LOC,II,JJ) * TMASS(FAM_LOC)
            !// Reset for diagnostics
            BTTBCK(L,FAM_LOC,II,JJ) = BTT(L,FAM_LOC,II,JJ)
          end do
        end do
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine set_fam_in_trop
  !// ----------------------------------------------------------------------

  !// ------------------------------------------------------------------
end module stratchem_oslo
!//=========================================================================

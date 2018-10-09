!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Oslo tropospheric chemistry.
!//=========================================================================
module tropchem_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: tropchem_oslo
  !// DECRIPTION: Routine for driving the Oslo tropospheric chemistry.
  !//
  !// Contains:
  !//   subroutine oslochem_trop
  !//
  !// Amund Sovde, November 2014 (from .f to .f90),
  !//              March-October 2009,
  !//              September-October 2008
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Global parameters
  !// Number of OH-iterations per chemical time step
  integer, parameter, private :: OH_HO2_ITER = 3
  !// ----------------------------------------------------------------------
  public oslochem_trop
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine oslochem_trop(BTT,BJV,BAIR_MOLEC,BVOL,BTEM,AIRB,BEMIS,DTCHM,MP)
    !// --------------------------------------------------------------------
    !// Tropospheric chemistry. Prepares the integration and calls integration
    !// routine for each column.
    !//
    !// This means that for each I,J (or II,JJ) the integration is done from
    !// the bottom of the column and up to LMTROP(I,J).
    !//
    !// Amund Sovde, October 2008
    !//              Updated March 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IDBLK, JDBLK, IPAR, JPAR, LPAR, &
         LWEPAR, NPAR, NOTRPAR, TRACER_ID_MAX, &
         LOSLOCSTRAT, LSULPHUR, LNITRATE, LSOA, LEMISDEP_INCHEM
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY, &
         ETAA, ETAB, YDGRD, PLAND
    use cmn_chem, only: TMASS
    use cmn_fjx, only: JVN_
    use cmn_met, only: P, Q, CLDFR, CLDLWC, CLDIWC, PRECLS
    use cmn_sfc, only: VDEP
    use cmn_oslo, only: LMTROP, PR42HET, trsp_idx, chem_idx, &
         Xtrsp_idx, Xchem_idx, XTMASS, XSTT, &
         LVS2ADD2TC, TROPCHEMnegO3
    use cmn_parameters, only: M_AIR, AVOGNR, R_AIR
    use chem_oslo_rates, only: TEMPRANGE, MINTEMP, TCRATE_TP_IJ_TRP, &
         TCRATE_onAER
    use sulphur_oslo, only: R3877, R4076, TCRATE_TP_S_IJ
    use utilities_oslo
    use pchemc_ij
    use emisdep4chem_oslo, only: getEMISX, getVDEP_oslo
    use soa_oslo, only: SOA_v9_separate, soa_diag_drydep, &
         soa_diag_separate, soa_diag_soagas
    use diagnostics_general, only: nchemdiag, save_chemPL, save_chemOxPL
    use diagnostics_scavenging, only: scav_diag_put_ddep
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)  :: DTCHM
    integer, intent(in) :: MP
    real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK),intent(in) :: BJV
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BAIR_MOLEC,BVOL
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: BTEM,AIRB
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(in) :: BEMIS
    !// In/out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT

    !// Locals
    !// QSSA variables
    real(r8) :: ST1,QLIN1,QTEST1,DTCH1, ST2,QLIN2,QTEST2,DTCH2

    !// Max chemical time step in troposphere is 15min = 900sec.
    real(r8),parameter :: DTCHM_MAX = 300._r8

    !// Number of iterations in troposphere
    integer :: NCHEM_ITER

    !// Pressure variables
    real(r8) :: P1, P2

    !// Column arrays ------------------------------------------------------
    !// Pressure & layer thickness
    real(r8),dimension(LPAR) :: PMIDL, PUEDG

    !// Air density [molecules / cm3]
    real(r8),dimension(LPAR) :: AIR_MOLEC

    !// H2O number density [molecules / cm3]
    real(r8),dimension(LPAR) :: H2O_MOLEC

    !// Temperature [K], box volume [m3] and mid-level height
    real(r8),dimension(LPAR) :: TEMP, DV

    !// Tracer concentration [molecules / cm3]
    real(r8),dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL, ZC0

    !// Photolytical reaction rates
    real(r8),dimension(JVN_,LPAR) :: JVALUES

    !// Dry deposition for Oslo method
    real(r8),dimension(TRACER_ID_MAX) :: VDEP_OC

    !// Reaction rates dependent on pressure and/or temperature
    real(r8), dimension(LPAR) :: &
          r_op_o2_m, &     !OP + O2 + M --> O3 + M
          r_ho2_ho2_tot, & !HO2 + HO2 (+ M) --> H2O2 + O2 (+ M)
          r_oh_oh_m, &     !OH + OH + M --> H2O2 + M
          r_no2_oh_m, &    !NO2 + OH + M --> HNO3 + M
          r_no2_no3_m, &   !NO2 + NO3 + M --> N2O5 + M
          r_n2o5_m, &      !N2O5 + M --> NO2 + NO3 + M
          r_ho2_no2_m, &   !NO2 + HO2 --> HO2NO2
          r_ho2no2_m, &    !HO2NO2 + M --> HO2 + NO2 + M
          r_oh_hno3, &     !OH + HNO3 --> NO3 (+ H2O)
          r_oh_co_a, &     !OH + CO --> HOCO --O2--> HO2 + CO2
          r_oh_co_b, &     !OH + CO --> H + CO2 --O2--> HO2 + CO2
          r_oh_c2h4_m, &   !OH + C2H4 --> (HOCH2CH2) --> CH3 + HCHO
          r_oh_c3h6_m, &   !OH + C3H6 + (O2) --> CH3XX, & CH3XX = CH3CH(O2)CH2OH
          r_ch3_o2_m, &    !CH3 + O2 +M --> CH3O2 + M
          r_oh_hcohco_m_a, & !OH + HCOHCO --> (HCO+CO)
          r_oh_hcohco_m_b, & !OH + HCOHCO --> CH3CO3 and 2CO + HO2
          r_no2_ch3x_m, &  !NO2 + CH3X (i.e. CH3C(O)O2) + M --> PAN + M
          r_pan_m, &       !PAN (i.e. CH3C(O)O2NO2) + M --> CH3X + NO2 + M
          r_no_ho2_b, &    !NO + HO2 --> HNO3
          r_op_no_m, &     !OP + NO + M --> NO2
          r_op_no2_m       !OP + NO2 + M --> NO3

    !// Sulfur T,p reaction rates
    real(r8), dimension(LPAR) :: &
         R4071b, &
         RTOT4072, &
         RAQ0172, RAQ1572,  RAQ1772, &
         RCATSO2

    !// Aerosol reaction rates
    real(r8), dimension(LPAR) :: RR_N2O5_H2O_AER, RR_HO2_AER, RR_QAER, &
         RR_NO2_SOOT

    !// Emissions [molec/cm3/s]
    !// Includes ALL emissions, fetched from BEMIS
    real(r8), dimension(TRACER_ID_MAX,LPAR) :: EMISX


    !// Other
    real(r8), dimension(nchemdiag,TRACER_ID_MAX,LPAR) :: CHEMLOSS
    real(r8), dimension(nchemdiag,TRACER_ID_MAX,LPAR) :: CHEMPROD
    real(r8), dimension(LPAR) :: OxCHEMLOSS
    real(r8), dimension(LPAR) :: OxCHEMPROD
    real(r8), dimension(LPAR) :: DTCHEMLOSS
    real(r8), dimension(TRACER_ID_MAX) :: drydepDIAG
    real(r8), dimension(LPAR) :: CIWC, CLWC, CFR !// cloud parameters
    real(r8), dimension(LPAR) :: COUNTnegO3

    logical, parameter :: LM7 = .false.


    !// --------------------------------------------------------------------

    !// Some helping constants
    real(r8), parameter :: Q2H2O  = M_AIR/18._r8  !// TMMVV for H2O

    !// For looping
    integer :: I,J,L,II,JJ,LL

    !// Height level of the tropopause and levels of cloud parameters
    integer :: LMTP, LCLD
    !// --------------------------------------------------------------------

    !// Find number of loops in chemistry to match the global
    !// chemical time step, i.e. DTCHM_MAX or shorter.
    if (DTCHM.lt.DTCHM_MAX) then
      !// Global time step is shorter than DTCHM_MAX
      DTCH1 = DTCHM
      NCHEM_ITER = 1
    else
      NCHEM_ITER = int(DTCHM/DTCHM_MAX + 0.5_r8)
      if (mod(DTCHM,DTCHM_MAX).eq.0._r8) then
        !// DTCH1 is given by max time step.
        DTCH1 = DTCHM_MAX
      else
        !// Must change DTCH (i.e. shorter than DTCH_MAX)
        DTCH1 = DTCHM / real(NCHEM_ITER, r8)
      end if
    end if


    !// Set QSSA parameters for this time step
    ST1    = 10._r8 / DTCH1
    QLIN1  = 0.1_r8 / DTCH1
    QTEST1 = 1.0_r8 / DTCH1

    !// For OH cycling (fast) we define other QSSA constants
    DTCH2 = DTCH1 /real(OH_HO2_ITER,r8)
    ST2   = ST1   *real(OH_HO2_ITER,r8) !// i.e. 1.d+1/(DTCH1/OH_HO2_ITER)
    QLIN2 = QLIN1 *real(OH_HO2_ITER,r8)
    QTEST2= QTEST1*real(OH_HO2_ITER,r8)

    !// Zero emissions for chemistry (mostly important if separate
    !// SOURCE process is treated)
    EMISX(:,:) = 0._r8
    VDEP_OC(:) = 0._r8

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// For tropospheric chemistry only, integrate up to max top of
      !// troposphere at this latitude (as in Oslo CTM2). If stratospheric
      !// chemistry is included, this is overwritten in the I-loop below.
      !// When stratchem is off, do tropchem LVS2ADD2TC levels above "normal":
      LMTP = maxval(LMTROP(:,J)) + LVS2ADD2TC

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Tropopause height when using stratospheric chemistry
        if (LOSLOCSTRAT) LMTP = LMTROP(I,J)


        !// Generating vertical arrays
        do L = 1, LMTP
          !// Box thickness (for IJ build from ETAA/B)
           P2 = ETAA(L+1) + ETAB(L+1) * P(I,J)
           P1 = ETAA(L) + ETAB(L) * P(I,J)

           !// Box center pressure
           PMIDL(L) = (P1 + P2) * 0.5_r8
           !// Box top pressure
           PUEDG(L) = P2
         end do !// do L = 1, LMTP

         !// Box volume [m3]
         DV(1:LMTP) = BVOL(1:LMTP,II,JJ)

         !// Calculate air density [molec/cm^3]
         AIR_MOLEC(1:LMTP) = BAIR_MOLEC(1:LMTP,II,JJ)

         !// H2O vapor density [molec/cm^3]
         H2O_MOLEC(1:LMTP) = Q(I,J,1:LMTP)*AIR_MOLEC(1:LMTP)*Q2H2O

         !// Temperature
         TEMP(1:LMTP) = BTEM(1:LMTP,II,JJ)

         !// J-values
         JVALUES(:,1:LMTP) = BJV(:,1:LMTP,II,JJ)

         !// Cloud parameters
         LCLD = min(LMTP,LWEPAR)
         do L = 1, LCLD
           CFR(L)  = CLDFR(I,J,L)
           CLWC(L) = CLDLWC(I,J,L)
           CIWC(L) = CLDIWC(I,J,L)
         end do
         do L = LCLD+1,LPAR
           CFR(L) = 0._r8
           CLWC(L)= 0._r8
           CIWC(L)= 0._r8
         end do

         !// Treat emissions/deposition as production/loss terms in chemistry?
         if (LEMISDEP_INCHEM) then
           !// Get emissions (incl. all emis.; surface, lightning, etc.)
           !// Routine sets EMISX, conversion from kg/s to molec/cm3/s
           !// is carried out.
           call getEMISX(BEMIS(:,:,II,JJ),DV,I,J,II,JJ,1,LMTP, &
                LPAR, MP, NPAR, TRACER_ID_MAX, chem_idx, TMASS, EMISX)
           !// Get dry deposition from VDEP, convert from m/s to 1/s
           call getVDEP_oslo(I,J,VDEP_OC,VDEP(:,I,J),DV(1)/AREAXY(I,J), &
                TRACER_ID_MAX,chem_idx,NPAR)
        end if

        !// Set local arrays (ZC_LOCAL); given in [kg/gridbox]
        call gotoZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,1,LMTP, &
             IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
             NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)


        !// Aerosol scavenging
        call TCRATE_onAER(I,J, TEMP, CIWC, CLWC, CFR, DV, PMIDL, &
             AIRB(:,II,JJ), AREAXY(I,J), LMTP, PLAND(I,J), ZC_LOCAL, &
             RR_N2O5_H2O_AER, RR_HO2_AER, RR_QAER, RR_NO2_SOOT)

        !// Save aerosol uptake rates
        !call diagnose_rr(I,J, RR_N2O5_H2O_AER, RR_HO2_AER, &
        !     RR_QAER, RR_NO2_SOOT)

        !// Save ZC
        ZC0(:,:) = ZC_LOCAL(:,:)

        !// ZC from mass to concentration [molec/cm3]
        call ZC_MASS2CONC(ZC_LOCAL,1,LMTP,LPAR,NPAR,NOTRPAR, &
             TRACER_ID_MAX, trsp_idx, Xtrsp_idx,TMASS,XTMASS,DV,AVOGNR)



        !// Calculate rates dependent on T & P
        Call TCRATE_TP_IJ_TRP(I, J, LPAR, LMTP, &
             TEMP, AIR_MOLEC, H2O_MOLEC, PMIDL, &
             r_op_o2_m, r_ho2_ho2_tot, r_oh_oh_m, r_no2_oh_m, r_no2_no3_m, &
             r_n2o5_m, r_ho2_no2_m, r_ho2no2_m, r_oh_hno3, &
             r_oh_co_a, r_oh_co_b, r_oh_c2h4_m, r_oh_c3h6_m, r_ch3_o2_m, &
             r_oh_hcohco_m_a, r_oh_hcohco_m_b, r_no2_ch3x_m, r_pan_m, &
             r_no_ho2_b, r_op_no_m, r_op_no2_m)

        !// Initialize/set sulphur reaction rates dependent on T & p
        !// and on dissolved fractions.
        call TCRATE_TP_S_IJ(LPAR, LMTP, TEMP, AIR_MOLEC, &
                          H2O_MOLEC, PMIDL, LSULPHUR, ydgrd(J), &
                          AREAXY(I,J), &
                          LWEPAR, CFR(:), CLWC(:), &
                          CIWC(:),PRECLS(I,J,:), &
                          DV, AIRB(:,II,JJ), &
                          R4071b, RTOT4072,RAQ0172, RAQ1572, RAQ1772, &
                          RCATSO2  )

        !// Integrate vertically
        call OSLO_CHEM ( &
             LPAR, TEMPRANGE, MINTEMP, JVN_, TRACER_ID_MAX, I, J, &
             NCHEM_ITER, OH_HO2_ITER, LMTP, LOSLOCSTRAT, &
             ZC_LOCAL, EMISX,  JVALUES,  TEMP, &
             AIR_MOLEC, H2O_MOLEC, &
             RR_QAER, RR_N2O5_H2O_AER, RR_HO2_AER, RR_NO2_SOOT, &
             !// T&p dependent
             r_op_o2_m, r_ho2_ho2_tot, r_oh_oh_m, r_no2_oh_m, r_no2_no3_m, &
             r_n2o5_m, r_ho2_no2_m, r_ho2no2_m, r_oh_hno3, &
             r_oh_co_a, r_oh_co_b, r_oh_c2h4_m, r_oh_c3h6_m, r_ch3_o2_m, &
             r_oh_hcohco_m_a, r_oh_hcohco_m_b, r_no2_ch3x_m, r_pan_m, &
             r_no_ho2_b, r_op_no_m, r_op_no2_m, &
             RAQ0172, RAQ1572,  RAQ1772, &
             R4071b,  RTOT4072, RCATSO2, &
             !// Constants etc
             NPAR, trsp_idx, TMASS, &
             DTCH1, ST1, QLIN1, QTEST1, &
             DTCH2, ST2, QLIN2, QTEST2, &
             VDEP_OC, &
             LSULPHUR, LNITRATE, LSOA, LM7, &
             drydepDIAG, nchemdiag, CHEMLOSS, CHEMPROD, &
             OxCHEMLOSS, OxCHEMPROD, &
             COUNTnegO3 )

        !// Save negative O3 for this MP. This is reported after
        !// chemistry in diagnostic subroutine OC_REPORTS.
        TROPCHEMnegO3(:,MP) = TROPCHEMnegO3(:,MP) + COUNTnegO3(:)

        !// Save dry deposited data
        call scav_diag_put_ddep(drydepDIAG,DV(1),II,JJ,MP)

        !// Diagnose chemistry losses and prods
        call save_chemPL(CHEMLOSS,CHEMPROD,DV,PUEDG,I,J,II,JJ,MP, &
                         1,LMTP,LMTP)
        call save_chemOxPL(OxCHEMLOSS,OxCHEMPROD,DV,I,J,1,LMTP)

        !// ZC from concentration to mass
        call ZC_CONC2MASS(ZC_LOCAL,1,LMTP,LPAR,NPAR,NOTRPAR, &
             TRACER_ID_MAX, trsp_idx, Xtrsp_idx,TMASS,XTMASS,DV,AVOGNR)

        if (LSOA) then
          !// Accumulate SOA loss
          call soa_diag_drydep(drydepDIAG,DV(1),I,J)

          !// Save change in SOAGS
          ZC0(:,:) = ZC_LOCAL(:,:) - ZC0(:,:)
          call soa_diag_soagas(ZC0,I,J,1,LMTP)

          !// SOA separation treated in mass space (will be repeated
          !// in stratchem)
          ZC0(:,:) = ZC_LOCAL(:,:)
          call SOA_v9_separate(ZC_LOCAL,TEMP,DV,I,J,1,LMTP)
          ZC0(:,:) = ZC_LOCAL(:,:) - ZC0(:,:)
          !// Accumulate SOA production (from level 1 to LMTP)
          call soa_diag_separate(ZC0,I,J,1,LMTP)
        end if !// if (LSOA) then

        !// Move local arrays back to BTT / XSTT
        call backfromZC_IJ(ZC_LOCAL,BTT,XSTT,I,J,II,JJ,1,LMTP, &
             IPAR,JPAR,LPAR,IDBLK,JDBLK,MPBLKIB,MPBLKIE,MP, &
             NPAR, NOTRPAR, TRACER_ID_MAX, chem_idx, Xchem_idx)

      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)

    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine oslochem_trop
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
end module tropchem_oslo
!//=========================================================================

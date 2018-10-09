!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, June 2016
!//=========================================================================
!// Column driver for Oslo tropospheric chemistry.
!//=========================================================================
module pchemc_ij
  !// ----------------------------------------------------------------------
  !// MODULE: pchemc_ij
  !// DECRIPTION: Module for OSLO_CHEM, the column driver for Oslo
  !//             tropospheric chemistry.
  !//
  !// Contains:
  !//   subroutine oslo_chem
  !//
  !// Amund Sovde Haslerud, June 2016
  !// Ole Amund Sovde, December 2014 (from .f to .f90)
  !//                  October 2009
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'pchemc_ij.f90'
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine OSLO_CHEM ( &
       LM, ITM, MINTEMP,NPHM, &
       TRACER_ID_MAX, ICOL, JCOL, NCHEM_ITER, OH_HO2_ITER, &
       LMTROP, LOSLOCSTRAT, &
       ZC,   EMISX,  JV,  TEMP, &
       AIR_MOLEC,H2O_MOLEC, &
       QAER, RR_N2O5_H2O_AER, RR_HO2_AER, RR_NO2_SOOT, &
       !// T & p dependent reactions
       r_op_o2_m, r_ho2_ho2_tot, r_oh_oh_m, r_no2_oh_m, r_no2_no3_m, &
       r_n2o5_m, r_ho2_no2_m, r_ho2no2_m, r_oh_hno3, &
       r_oh_co_a, r_oh_co_b, r_oh_c2h4_m, r_oh_c3h6_m, r_ch3_o2_m, &
       r_oh_hcohco_m_a, r_oh_hcohco_m_b, r_no2_ch3x_m, r_pan_m, &
       r_no_ho2_b, r_op_no_m, r_op_no2_m, &
       RAQ0172, RAQ1572,  RAQ1772, &
       R4071b,  RTOT4072, RCATSO2, &
       !// Transport index (trsp_idx) and tracer mass (TMASS)
       NPAR, trsp_idx, TMASS, &
       !// QSSA stuff
       DTCH1, ST1, QLIN1, QTEST1, &
       DTCH2, ST2, QLIN2, QTEST2, &
       VDEP_OC, &
       LSULPHUR, LNITRATE, LSOA, LM7, &
       DDDIAG,nchemdiag,CHEMLOSS,CHEMPROD, &
       OxCHEMLOSS, OxCHEMPROD, &
       COUNTnegO3)
    !// --------------------------------------------------------------------
    !//
    !// Column driver for integrating Oslo Chemistry in the troposphere using
    !// the QSSA integrator.
    !//
    !// Some important variable names:
    !// LM: The total number of layers.
    !// L: Layer loop index.
    !// LMTROP: Tropopause, i.e. the LM of TROPosphere.
    !// LOSLOCSTRAT: Logical to tell whether stratospheric chemistry is done.
    !// NCHEM_ITER: number of internal iterations to match chemistry time step.
    !// QSSA: Moved into separate module to be accessible in the stratosphere.
    !// TRACER_ID_MAX: Max number of tracers.
    !// ICOL, JCOL: the I and J for this column.
    !// Vxx: Dry depositions are available for all components in VDEP_OC.
    !// EMISX: Emission for each species (also lightning and aircraft).
    !//
    !// Ole Amund Sovde, October 2009; October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: LDEBUG
    use qssa_integrator, only: qssa
    use strat_h2o, only: LOLD_H2OTREATMENT
    use chem_oslo_rates, only: &
         !// Constant reaction rates
         r_od_ch4_a, r_od_ch4_b, r_od_ch4_c, r_od_h2, &
         r_oh_c6h14, r_oh_c6hxr, &
         r_cho_o2, &
         r_ch3co_o2, &
         r_ch3o2_c2h5o2, r_ch3o2_c3h7o2, r_ch3o2_c4h9o2, &
         r_ch3o2_c6h13o2, r_ch3o2_ch3cob, r_ch3o2_ch3cod, r_ch3o2_ch3xx, &
         r_ch3o2_isor1, r_ch3o2_isor2, r_oh_ar2, &
         r_no_c3h7o2, &
         !// Temperature dependent rates
         r_od_m, r_od_h2o, r_op_no2, &
         r_o3_no, r_o3_no2, r_o3_oh, r_o3_ho2, r_o3_c2h4, r_o3_c3h6, &
         r_oh_h2o2, r_oh_ho2, r_oh_h2, r_oh_ho2no2, r_oh_ch4, &
         r_oh_ch3o2h_a, r_oh_ch3o2h_b, r_oh_pan, r_oh_c2h6, r_oh_c3h8, &
         r_oh_c4h10, r_oh_isoprene, r_oh_ch2o, r_oh_ch3cho, &
         r_oh_rcohco, r_oh_isok, r_oh_aceton, r_oh_ch3cox, &
         r_oh_dms_a, r_oh_h2s, r_oh_ch3oh, &
         r_no_ho2, r_no_ch3o2, r_no_c2h5o2, r_no_c4h9o2, r_no_c6h13o2, &
         r_no_ar1, r_no_ar3, r_no_isor1, r_no_isor2, &
         r_no_ch3cob, r_no_ch3x, r_no_ch3cod, r_no_ch3xx, r_no_no3, &
         r_no3_dms, r_no3_ch3cho, r_no2_no3_b, &
         r_ho2_ch3o2, r_ho2_ch3x, &
         r_ch3o2_ch3o2, r_ch3o2_ch3x_a, r_ch3o2_ch3x_b, r_ch3x_ch3x, &
         r_ch3o_o2, r_ho2_radical, &
         r_o3_soaC1, r_oh_soaC1, r_no3_soaC1, &
         r_o3_soaC2, r_oh_soaC2, r_no3_soaC2, &
         r_o3_soaC3, r_oh_soaC3, r_no3_soaC3, &
         r_o3_soaC4, r_oh_soaC4, r_no3_soaC4, &
         r_o3_soaC5, r_oh_soaC5, r_no3_soaC5, &
         r_o3_soaC7, r_oh_soaC7, r_no3_soaC7, &
         r_o3_soaC8, r_oh_soaC8, r_no3_soaC8, &
         r_oh_benzene, &
         !// Branching ratios
         fa_no_ch3o2, fa_no_c2h5o2, fa_no_c3h7o2, fa_no_c4h9o2, &
         fa_no_c6h13o2, fa_no_ch3cob, fa_no_ch3cod, fa_no_isor1, &
         fb_no_c2h5o2, fb_no_c3h7o2, fb_no_c4h9o2

    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LM, ITM, TRACER_ID_MAX, NPHM  !// Array sizes
    Integer, intent(in) :: MINTEMP,ICOL,JCOL
    integer, intent(in) :: OH_HO2_ITER, NCHEM_ITER  !// Number of iterations
    integer, intent(in) :: LMTROP                   !// Tropopause level
    logical, intent(in) :: LOSLOCSTRAT !// Do stratospheric chemistry also?
    real(r8), dimension(TRACER_ID_MAX,LM), intent(in) :: EMISX !// emissions
    real(r8), dimension(NPHM,LM), intent(in)          :: JV    !// J-values
    real(r8), dimension(LM), intent(in)  :: &
         AIR_MOLEC,       & !// air conc.
         H2O_MOLEC,       & !// water conc.
         QAER,            & !// Aerosol uptake (old simple param.)
         RR_N2O5_H2O_AER, & !// Aerosol uptake N2O5+H2O
         RR_HO2_AER,      & !// Aerosol uptake HO2
         RR_NO2_SOOT        !// Soot uptake NO2
    !// Temperature and pressure dependent rates
    real(r8), dimension(LM), intent(in) :: &
         r_op_o2_m, r_ho2_ho2_tot, r_oh_oh_m, r_no2_oh_m, r_no2_no3_m, &
         r_n2o5_m, r_ho2_no2_m, r_ho2no2_m, r_oh_hno3, &
         r_oh_co_a, r_oh_co_b, r_oh_c2h4_m, r_oh_c3h6_m, r_ch3_o2_m, &
         r_oh_hcohco_m_a, r_oh_hcohco_m_b, r_no2_ch3x_m, r_pan_m, &
         r_no_ho2_b, r_op_no_m, r_op_no2_m, &
         RAQ0172, RAQ1572,  RAQ1772, &
         R4071b,  RTOT4072, RCATSO2
    !// Tracer info
    integer, intent(in) :: NPAR, trsp_idx(TRACER_ID_MAX)
    real(r8), intent(in)  :: TMASS(NPAR)
    !// QSSA parameters
    real(r8), intent(in) :: &
         DTCH1, ST1, QLIN1, QTEST1, &!// Regular chemistry
         DTCH2, ST2, QLIN2, QTEST2  !// OH/HO2 chemistry

    !// Meteorological variables
    real(r8), dimension(LM), intent(in) :: TEMP  !// Temperature

    !// Dry deposition ala CTM2
    real(r8), dimension(TRACER_ID_MAX), intent(in) :: VDEP_OC


    !// Logical switches
    logical, intent(in) :: &
         LSULPHUR,  &! Sulphur chemistry
         LNITRATE,  &! Nitrate chemistry
         LSOA,      &! Secondary organic aerosols
         LM7         ! True if using M7 module

    !// Input/Output
    real(r8), dimension(TRACER_ID_MAX,LM), intent(inout):: ZC  !// tracer conc.

    !// Diagnostics for production and loss [molec/cm3/s]
    integer, intent(in) :: nchemdiag
    real(r8),dimension(nchemdiag,TRACER_ID_MAX,LM), intent(out) :: CHEMLOSS
    real(r8),dimension(nchemdiag,TRACER_ID_MAX,LM), intent(out) :: CHEMPROD
    real(r8),dimension(LM), intent(out) :: OxCHEMLOSS
    real(r8),dimension(LM), intent(out) :: OxCHEMPROD
    real(r8),dimension(TRACER_ID_MAX), intent(out) :: DDDIAG

    !// Count number of corrected negative ozone
    real(r8), intent(out) :: COUNTnegO3(LM)
    !// --------------------------------------------------------------------
    !// Local parameters
    !// Chemical components names
    real(r8) :: &
         M_NO,  M_NO2,   M_HNO3, M_NO3,  M_N2O5, M_HO2NO2, M_PAN, &
         M_O3,  M_PANX,  M_CO,   M_H2O2, M_HO2,  &
         M_O3P, M_O1D,   M_OH,   M_CH4, &
         M_C2H4,  M_C2H6,   M_C3H6,  M_C3H8, M_C4H10, M_ARAD, &
         M_CH2O,  M_CH3CHO, M_HCOHCO, M_RCOHCO, &
         M_ISOPREN, M_ISOR1, M_ISOK,   M_ISOR2, &
         M_C6H14,  M_C6HXR, &
         M_CH3CO, &
         M_CH3O2H, M_CH3COY, M_CH3COX, &
         M_CH3O2, M_C2H5O2, M_C4H9O2, M_C6H13O2, M_CH2O2OH, M_CH3COB, &
         M_C3H7O2, M_ACETON,M_CH3XX,  M_AR1, M_AR2, M_AR3, &
         M_CH3X, M_CH3COD,  &
         M_O2,    M_H2O, &
         M_DMS, M_SO2, M_H2S, M_MSA, M_CS2, M_OCS, &
         M_SO4, SO4gas, SO4aq, &
         NH3gas, NH4fine, NH4coarse, NO3fine, NO3coarse, &
         M_Apine, M_Bpine, M_Limon, M_Myrcene, M_Sabine, &
         M_D3carene, M_Ocimene, &
         M_Trpolene, M_Trpinene, M_TrpAlc, M_Sestrp, M_Trp_Ket, M_Tolmatic, &
         M_Benzene, M_C6HXR_SOA, &
         M_H2, &
         AIRMOLEC, &
         !// Short-lived, steady state, not transported
         M_CH3, M_CH3O, M_CHO, M_O3NO, &
         XNOX,   XHO2NO2, &
         !// Family species calculated only in chemistry
         M_NOX, M_NOZ

    !// Local reaction rates
    real(r8) :: &
         !// Constant rates
         k_od_ch4_a, k_od_ch4_b, k_od_ch4_c, k_od_h2, &
         k_oh_c6h14, k_oh_c6hxr, &
         k_cho_o2, &
         k_ch3co_o2, &
         k_ch3o2_c2h5o2, k_ch3o2_c3h7o2, k_ch3o2_c4h9o2, k_ch3o2_c6h13o2, &
         k_ch3o2_ch3cob, k_ch3o2_ch3cod, k_ch3o2_isor1,  k_ch3o2_isor2, &
         k_ch3o2_ch3xx, &
         k_oh_ar2, &
         k_no_c3h7o2

    real(r8) :: &
         !// Temperature dependent rates
         k_od_m, k_od_h2o, k_op_no2, &
         k_o3_ho2, k_o3_oh, k_o3_no, k_o3_no2, k_o3_c2h4, k_o3_c3h6, &
         k_oh_h2o2, k_oh_ho2, k_oh_h2, k_oh_ho2no2, k_oh_ch4, &
         k_oh_ch3o2h_a, k_oh_ch3o2h_b, k_oh_pan, k_oh_c2h6, k_oh_c3h8, &
         k_oh_c4h10, k_oh_isoprene, k_oh_ch2o, k_oh_ch3cho, &
         k_oh_rcohco, k_oh_isok, k_oh_aceton, k_oh_ch3cox, &
         k_oh_dms_a, k_oh_h2s, k_oh_ch3oh, &
         k_no_ho2, k_no_ch3o2, k_no_c2h5o2, k_no_c4h9o2, k_no_c6h13o2, &
         k_no_ar1, k_no_ar3, k_no_isor1, k_no_isor2, &
         k_no_ch3cob, k_no_ch3x, k_no_ch3cod, k_no_ch3xx, k_no_no3, &
         k_no3_dms, k_no3_ch3cho, k_no2_no3_b, &
         k_ho2_ch3o2, k_ho2_ch3x, k_ho2_radical, &
         k_ch3o2_ch3o2, k_ch3o2_ch3x_a, k_ch3o2_ch3x_b, k_ch3x_ch3x, &
         k_ch3o_o2

    real(r8) :: &
         !// SOA Rates for hydrocarbon oxidation
         k_o3_soaC1, k_oh_soaC1, k_no3_soaC1, &
         k_o3_soaC2, k_oh_soaC2, k_no3_soaC2, &
         k_o3_soaC3, k_oh_soaC3, k_no3_soaC3, &
         k_o3_soaC4, k_oh_soaC4, k_no3_soaC4, &
         k_o3_soaC5, k_oh_soaC5, k_no3_soaC5, &
         !// SOA Aromatic oxidation rates
         k_o3_soaC7, k_oh_soaC7, k_no3_soaC7, &
         k_o3_soaC8, k_oh_soaC8, k_no3_soaC8, &
         k_oh_benzene

    real(r8) :: &
         !// Temperature and pressure dependent rates
         k_op_o2_m, k_ho2_ho2_tot, k_oh_oh_m, k_no2_oh_m, k_no2_no3_m, &
         k_n2o5_m, k_ho2_no2_m, k_ho2no2_m, k_oh_hno3, &
         k_oh_co_a, k_oh_co_b, k_oh_c2h4_m, k_oh_c3h6_m, k_ch3_o2_m, &
         k_oh_hcohco_m_a, k_oh_hcohco_m_b, k_oh_hcohco_m_c, &
         k_no2_ch3x_m, k_pan_m, &
         k_no_ho2_b, k_op_no_m, k_op_no2_m 

    real(r8) :: &
         !// Aerosol uptake/conversion rates
         k_n2o5_h2o_aer

    real(r8) :: &
         CAQ0172, &
         CAQ1572, &
         CAQ1772, &
         C4071b, CTOT4072, &
         CCATSO2, &
         !// J-values
         DO2, &
         DAO3,   DBO3,   DO3,    DNO2,   DH2O2, DHNO3,DACH2O, &
         DBCH2O, DCH3CHO,DCH3COX,DHO2NO2,DCH2O, DCH3O2H, &
         DCH3COY,DRCOHCO,DHCOHCO,DNO3,   DN2O5, &
         DACETON_A,      DACETON_B,      DPAN,   DOCS, &
         !// production and loss rates
         PROD,PROD_2,LOSS,LOSS_2,LOSS_3, &
         PROD_NO,LOSS_NO, PROD_OH,LOSS_OH, PROD_HO2,LOSS_HO2, &
         LOSS_O3,LOSS_NO2,TMP_PL,TMP_PL2, &
         LOSS_C1,LOSS_C2,LOSS_C3,LOSS_C4,LOSS_C5, &
         !// balancing Nitrogen components
         XJNO,XJNO2, RELATN,EKSTRSN, O3TEST, &
         !// help variables
         RCH3XP, &
         OH_OLD, OH_NEW, HO2_OLD, HO2_NEW, &
         RLIM1,RLIM2,TMP_ZC,ZTMP, &
         DTCH,ST,QLIN,QTEST
    !// Emissions
    real(r8) :: POLLX(TRACER_ID_MAX)

    !// Dry deposition array used in calculations (=VDEP_OC in layer 1, else 0)
    real(r8), dimension(TRACER_ID_MAX) :: VDEP_L

    !// loop counters
    integer :: L, N, IIJ, IJL, JTEMP, NST, SOA_N

    !// OH iteration
    !// Max number of iterations for OH/HO2 within OH_HO2_ITER loop
    integer, parameter :: ITER_OH_MAX = 99
    !// When to stop iterating: (OH_NEW-OH_OLD)/OH_OLD < limit_oh_ho2
    real(r8), parameter :: limit_oh_ho2 = 1.e-2_r8

    !// SOA Secondary organic aerosols: Mass based stoichiometric
    !// coefficients. To determine the amount of each SOA GAS produced
    !// from the oxidation of precursor hydrocarbons
    !//
    !// Note that the SOAGAS products do NOT account for all of these
    !// reactions: The mass stoichiometric coefficients makes only a part
    !// of the products as SOAGAS, while the rest of the products are
    !// not treated by the model. E.g., for Apine+OH:
    !//   Apine + OH -> a121*SOAGAS11 + a122*SOAGAS12 + other_products
    !// where a121 and a122 are the mass stoichiometric coefficients.
    !// We take other_products into account for OH, but not for SOAGAS
    !// species. This is similar for other SOAGAS species.
    real(r8), parameter :: &
         SOA_STOC111 = 0.067_r8, SOA_STOC211 = 0.239_r8, &
         SOA_STOC311 = 0.069_r8, SOA_STOC411 = 0.067_r8, &
         SOA_STOC511 = 1._r8, &
         SOA_STOC112 = 0.345_r8, SOA_STOC212 = 0.363_r8, &
         SOA_STOC312 = 0.201_r8, SOA_STOC412 = 0.135_r8, &
         SOA_STOC512 = 0._r8, &
         SOA_STOC121 = 0.067_r8, SOA_STOC221 = 0.239_r8, &
         SOA_STOC321 = 0.069_r8, SOA_STOC421 = 0.067_r8, &
         SOA_STOC521 = 1._r8, &
         SOA_STOC122 = 0.345_r8, SOA_STOC222 = 0.363_r8, &
         SOA_STOC322 = 0.201_r8, SOA_STOC422 = 0.135_r8, &
         SOA_STOC522 = 0._r8, &
         SOA_STOC131 = 1._r8,    SOA_STOC231 = 1._r8, &
         SOA_STOC331 = 1._r8,    SOA_STOC431 = 1._r8, &
         SOA_STOC531 = 1._r8, &
         SOA_STOC132 = 0._r8,    SOA_STOC232 = 0._r8, &
         SOA_STOC332 = 0._r8,    SOA_STOC432 = 0._r8, &
         SOA_STOC532 = 0._r8
    !// Fractions of ISOPREN+OH giving SOAGAS61 and SOAGAS62
    real(r8), parameter :: &
         f_oh_isoprene_sg61 = 0.232_r8, &
         f_oh_isoprene_sg62 = 0.0288_r8
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CHEM'
    !// --------------------------------------------------------------------

    !// Initialize
    DDDIAG(:) = 0._r8
    CHEMLOSS(:,:,:) = 0._r8
    CHEMPROD(:,:,:) = 0._r8
    OxCHEMLOSS(:) = 0._r8
    OxCHEMPROD(:) = 0._r8
    COUNTnegO3(:)=0._r8

    !// --------------------------------------------------------------------
    do L = 1, LMTROP !// Integrate up to LMTROP
      !// ------------------------------------------------------------------

        
      !// Assign chemical reaction rates
      !// 1. Constant reaction rates: given as arguments for the subroutine
      !// 2. Dependent on T only
      !// 3. Dependent on T and p

      !// Assign chemical reaction rates, depedent on T only
      JTEMP = Nint(TEMP(L)) - MINTEMP !// Temperature index

      k_od_m   = r_od_m(JTEMP)
      k_od_h2o = r_od_h2o(JTEMP)
      k_op_no2 = r_op_no2(JTEMP)
      k_o3_ho2 = r_o3_ho2(JTEMP)
      k_o3_oh  = r_o3_oh(JTEMP)
      k_o3_no  = r_o3_no(JTEMP)
      k_o3_no2  = r_o3_no2(JTEMP)
      k_o3_c2h4 = r_o3_c2h4(JTEMP)
      k_o3_c3h6 = r_o3_c3h6(JTEMP)
      k_oh_h2o2 = r_oh_h2o2(JTEMP)
      k_oh_ho2  = r_oh_ho2(JTEMP)
      k_oh_h2 = r_oh_h2(JTEMP)
      k_oh_ho2no2 = r_oh_ho2no2(JTEMP)
      k_oh_ch4  = r_oh_ch4(JTEMP)
      k_oh_ch3o2h_a = r_oh_ch3o2h_a(JTEMP)
      k_oh_ch3o2h_b = r_oh_ch3o2h_b(JTEMP)
      k_oh_pan = r_oh_pan(JTEMP)
      k_oh_c2h6 = r_oh_c2h6(JTEMP)
      k_oh_c3h8 = r_oh_c3h8(JTEMP)
      k_oh_c4h10 = r_oh_c4h10(JTEMP)
      k_oh_isoprene = r_oh_isoprene(JTEMP)
      k_oh_ch2o = r_oh_ch2o(JTEMP)
      k_oh_ch3cho = r_oh_ch3cho(JTEMP)
      k_oh_rcohco = r_oh_rcohco(JTEMP)
      k_oh_isok = r_oh_isok(JTEMP)
      k_oh_aceton = r_oh_aceton(JTEMP)
      k_oh_ch3cox = r_oh_ch3cox(JTEMP)
      k_oh_dms_a = r_oh_dms_a(JTEMP)
      k_oh_h2s = r_oh_h2s(JTEMP)
      k_oh_ch3oh = r_oh_ch3oh(JTEMP)
      k_no_ho2 = r_no_ho2(JTEMP)
      k_no_ch3o2 = r_no_ch3o2(JTEMP)
      k_no_c2h5o2 = r_no_c2h5o2(JTEMP)
      k_no_c4h9o2 = r_no_c4h9o2(JTEMP)
      k_no_c6h13o2 = r_no_c6h13o2(JTEMP)
      k_no_ar1 = r_no_ar1(JTEMP)
      k_no_ar3 = r_no_ar3(JTEMP)
      k_no_isor1 = r_no_isor1(JTEMP)
      k_no_isor2 = r_no_isor2(JTEMP)
      k_no_ch3cob = r_no_ch3cob(JTEMP)
      k_no_ch3x = r_no_ch3x(JTEMP)
      k_no_ch3cod = r_no_ch3cod(JTEMP)
      k_no_ch3xx = r_no_ch3xx(JTEMP)
      k_no_no3 = r_no_no3(JTEMP)
      k_no3_dms = r_no3_dms(JTEMP)
      k_no3_ch3cho = r_no3_ch3cho(JTEMP)
      k_no2_no3_b = r_no2_no3_b(JTEMP)
      k_ho2_ch3o2 = r_ho2_ch3o2(JTEMP)
      k_ho2_ch3x = r_ho2_ch3x(JTEMP)
      k_ho2_radical = r_ho2_radical(JTEMP)
      k_ch3o2_ch3o2 = r_ch3o2_ch3o2(JTEMP)
      k_ch3o2_ch3x_a = r_ch3o2_ch3x_a(JTEMP)
      k_ch3o2_ch3x_b = r_ch3o2_ch3x_b(JTEMP)
      k_ch3x_ch3x = r_ch3x_ch3x(JTEMP)
      k_ch3o_o2 = r_ch3o_o2(JTEMP)

      !// SOA JTEMP reactions
      k_o3_soaC1 = r_o3_soaC1(JTEMP)
      k_oh_soaC1 = r_oh_soaC1(JTEMP)
      k_no3_soaC1= r_no3_soaC1(JTEMP)
      k_o3_soaC2 = r_o3_soaC2(JTEMP)
      k_oh_soaC2 = r_oh_soaC2(JTEMP)
      k_no3_soaC2= r_no3_soaC2(JTEMP)
      k_o3_soaC3 = r_o3_soaC3(JTEMP)
      k_oh_soaC3 = r_oh_soaC3(JTEMP)
      k_no3_soaC3= r_no3_soaC3(JTEMP)
      k_o3_soaC4 = r_o3_soaC4(JTEMP)
      k_oh_soaC4 = r_oh_soaC4(JTEMP)
      k_no3_soaC4= r_no3_soaC4(JTEMP)
      k_o3_soaC5 = r_o3_soaC5(JTEMP)
      k_oh_soaC5 = r_oh_soaC5(JTEMP)
      k_no3_soaC5= r_no3_soaC5(JTEMP)
      !// SOA JTEMP aromatic oxidation reactions
      k_o3_soaC7 = r_o3_soaC7(JTEMP)
      k_oh_soaC7 = r_oh_soaC7(JTEMP)
      k_no3_soaC7= r_no3_soaC7(JTEMP) !// CRHdont really use this, its low..
      k_o3_soaC8 = r_o3_soaC8(JTEMP)
      k_oh_soaC8 = r_oh_soaC8(JTEMP)
      k_no3_soaC8= r_no3_soaC8(JTEMP) !// CRH dont really use this, its low..
      k_oh_benzene = r_oh_benzene(JTEMP)

      !// Assign chemical reaction rates, dependent on T and P
      k_op_o2_m   = r_op_o2_m(L)
      k_ho2_ho2_tot = r_ho2_ho2_tot(L)
      k_oh_oh_m = r_oh_oh_m(L)
      k_no2_oh_m  = r_no2_oh_m(L)
      k_no2_no3_m = r_no2_no3_m(L)
      k_n2o5_m   = r_n2o5_m(L)
      k_ho2_no2_m  = r_ho2_no2_m(L)
      k_ho2no2_m   = r_ho2no2_m(L)
      k_oh_hno3  = r_oh_hno3(L)
      k_oh_co_a  = r_oh_co_a(L) !Revise when H is included (gives HO2)
      k_oh_co_b =  r_oh_co_b(L) !Revise when H is included (gives HO2)
      k_oh_c2h4_m  = r_oh_c2h4_m(L)
      k_oh_c3h6_m  = r_oh_c3h6_m(L)
      k_ch3_o2_m  = r_ch3_o2_m(L)
      k_oh_hcohco_m_a = r_oh_hcohco_m_a(L)
      k_oh_hcohco_m_b = r_oh_hcohco_m_b(L)
      k_oh_hcohco_m_c = r_oh_hcohco_m_b(L) ! same rate as B
      k_no2_ch3x_m  = r_no2_ch3x_m(L)
      k_pan_m = r_pan_m(L)
      k_no_ho2_b = r_no_ho2_b(L)
      k_op_no_m = r_op_no_m(L)
      k_op_no2_m = r_op_no2_m(L)

      C4071b  = R4071b(L)
      CTOT4072= RTOT4072(L)
      CAQ0172= RAQ0172(L)
      CAQ1572 = RAQ1572(L)
      CAQ1772 = RAQ1772(L)
      CCATSO2 = RCATSO2(L)


      k_n2o5_h2o_aer = RR_N2O5_H2O_AER(L)

        
      !// Assign constant rates (could use constants directly, but with this
      !// setup it is easy to change the rate to some dependency).
      k_od_ch4_a = r_od_ch4_a
      k_od_ch4_b = r_od_ch4_b
      k_od_ch4_c = r_od_ch4_c
      k_od_h2    = r_od_h2
      k_oh_c6h14 = r_oh_c6h14
      k_oh_c6hxr = r_oh_c6hxr
      k_ch3co_o2 = r_ch3co_o2
      k_cho_o2   = r_cho_o2
      k_ch3o2_c2h5o2  = r_ch3o2_c2h5o2
      k_ch3o2_c3h7o2  = r_ch3o2_c3h7o2
      k_ch3o2_c4h9o2  = r_ch3o2_c4h9o2
      k_ch3o2_c6h13o2 = r_ch3o2_c6h13o2
      k_ch3o2_ch3cob  = r_ch3o2_ch3cob
      k_ch3o2_ch3cod  = r_ch3o2_ch3cod
      k_ch3o2_ch3xx   = r_ch3o2_ch3xx
      k_ch3o2_isor1   = r_ch3o2_isor1
      k_ch3o2_isor2   = r_ch3o2_isor2
      k_oh_ar2    = r_oh_ar2
      k_no_c3h7o2 = r_no_c3h7o2

      !// Emissions
      do N = 1,TRACER_ID_MAX
         POLLX(N) = EMISX(N,L)
      end do

      !// Wet deposition
      !// Wet deposition is treated in p-scav_oc.f,
      !// and should not be included here. The new scavenging scheme
      !// also takes into account subcloud scavenging.

      !// Dry deposition
      if (L .gt. 1) then
         !// Array of all dry deps: Zero above ground
         VDEP_L(:) = 0._r8
      else
         !// Array of all dry deps: Set from input values
         VDEP_L(:) = VDEP_OC(:)
      end if

      !// Assign photolysis rates
      !// These do not change during internal DTCH time step,
      !// so it is set before the NCHEM_ITER loop.
      DO3    = JV( 1,L)          !// O3 + hv  -> O3P + O2
      DBO3   = JV( 2,L)          !// O3 + hv  -> O1D + O2
      DNO2   = JV( 3,L)          !// NO2 + hv -> NO + O3P
      DH2O2  = JV( 4,L)          !// H2O2 + hv -> 2OH
      DHNO3  = JV( 5,L)          !// HNO3 + hv -> OH + NO2
      DACH2O = JV( 6,L)          !// HCHO + hv -> HCO + H -2O2-> CO + HO2 + HO2
      DBCH2O = JV( 7,L)          !// HCHO + hv -> CO + H2
      DCH2O  = DACH2O + DBCH2O   !// HCHO + hv total
      DCH3CHO= JV( 8,L)          !// CH3CHO + hv -O2-> CH3O2 + HO2 + CO
      DHCOHCO= JV( 9,L)
      DRCOHCO= JV(10,L)
      DNO3   = JV(11,L)         &!// NO3 + hv -> NO + O2
               + JV(12,L)        !// NO3 + hv -> NO2 + O3P
      DN2O5  = JV(13,L)          !// N2O5 + hv -> NO2 + NO3
      DCH3O2H= JV(14,L)         &!// CH3O2H + hv -O2-> HCHO + OH + HO2
               + JV(15,L)        !// CH3O2H + hv -O2-> HCHO + O(3P) + HO2
      DHO2NO2= JV(16,L)         &!// HO2NO2 + hv -> HO2 + NO2
               + JV(17,L)
      DPAN   = JV(18,L)
      DCH3COX= JV(19,L)          !// CH3COC2H5 + hv -> not incl: CH3 + C2H5CO
                                 !//                -> CH3CO + C2H5
                                 !//               -O2-> CH3CO + C2H5O2
      DCH3COY= JV(20,L)          !// CH3COCOCH3 + hv -> 2 CH3CO
                                 !// Uses J-value from CH3COCHO
      if (LOSLOCSTRAT) then
         DO2    = JV(21,L)    !// O2 + hv  -> O3P + O3P
      else
         DO2    = 0._r8
      end if

      !// Old treatment, lacking J-values for acetone
      DACETON_A = 0.7_r8*DCH3COX         !// CH3COCH3 + hv -> CH3CO + CH3
      DACETON_B = 0._r8                  !// not included: -> 2CH3 + CO
      !// New treatment: get numbers from fast-JX
      !// DACETON_A = JV(50,L)           !// CH3COCH3 + hv -> CH3CO + CH3
      !// DACETON_B = JV(51,L)           !// CH3COCH3 + hv -> 2CH3 + CO
 
      !// ------------------------------------------------------------------
      !// Sub time step for chemistry ...
      !// ------------------------------------------------------------------
      do NST = 1, NCHEM_ITER
        !// ----------------------------------------------------------------
        !// Assign chemical components
        M_O2     = 0.2095_r8 * AIR_MOLEC(L)
        M_H2O    = H2O_MOLEC(L)

        M_O3     = ZC(1,L)
        M_HNO3   = ZC(4,L)
        M_PANX   = ZC(5,L) !// PANX = CH3COO2 + PAN
        M_CO     = ZC(6,L)
        M_C2H4   = ZC(7,L)
        M_C2H6   = ZC(8,L)
        M_C3H6   = ZC(9,L)
        M_C4H10  = ZC(10,L)
        M_C6H14  = ZC(11,L)
        !// C6HXR is m-Xylene
        M_C6HXR  = ZC(12,L) !// If (LSOA) will zero this below
        M_CH2O   = ZC(13,L)
        M_CH3CHO = ZC(14,L)
        M_H2O2   = ZC(15,L)
        M_HO2NO2 = ZC(17,L)
        M_CH3COY = ZC(18,L) !// CH3COY = CH3COCOCH3
        M_CH3COX = ZC(19,L) !// CH3COX = CH3COC2H5
        M_ISOPREN = ZC(20,L)
        !M_ISOPREN_TOT = ZC(20,L)           !// Need full isoprene for SOA and OH
        !M_ISOPREN_GAS = ZC(20,L) * 0.44_r8 !// For chemistry use 44%
        M_HO2    = ZC(21,L)
        M_CH3O2  = ZC(22,L)
        M_C2H5O2 = ZC(23,L)
        M_C4H9O2 = ZC(24,L)
        M_C6H13O2= ZC(25,L)
        M_CH3COB = ZC(27,L) !// CH3COB = CH3COCH(O2)CH3
        M_CH3XX  = ZC(28,L) !// CH3XX =  CH3CH(O2)CH2OH
        !// Aromatic products
        M_AR1    = ZC(29,L)
        M_AR2    = ZC(30,L)
        M_AR3    = ZC(31,L)
        M_ISOR1  = ZC(32,L)
        M_ISOK   = ZC(33,L)
        M_ISOR2  = ZC(34,L)
        M_HCOHCO = ZC(35,L)
        M_RCOHCO = ZC(36,L)
        M_CH3X   = ZC(37,L) !// CH3X = CH3COO2
        M_O3P    = ZC(38,L)
        M_O1D    = ZC(39,L)
        M_OH     = ZC(40,L)
        M_NO3    = ZC(41,L)
        M_N2O5   = ZC(42,L)
        M_NO     = ZC(43,L)
        M_NO2    = ZC(44,L)
        M_CH4    = ZC(46,L)
        M_C3H8   = ZC(48,L)
        M_C3H7O2 = ZC(49,L)
        M_ACETON = ZC(50,L)
        M_CH3COD = ZC(51,L) !// CH3COD = CH3COCH2(O2)  


        !// Sulphur chemistry
        M_DMS    = ZC(71,L)
        M_SO2    = ZC(72,L)
        !// Split SO4 into gas/aquous for M7
        if (.not. LM7) then
           M_SO4  = ZC(73,L)
        else
           SO4AQ  = ZC(73,L)
           SO4Gas = ZC(80,L)
        end if
        M_H2S    = ZC(74,L)
        M_MSA    = ZC(75,L)
        !// Not included yet.
        !M_CS2  = ZC(76,L) !// Not transported; should be initialized
        !M_OCS  = ZC(77,L) !// Not transported; should be initialized


        !// Nitrate chemistry
        NH3gas   = ZC(61,L)
        NH4fine  = ZC(62,L)
        NH4coarse= ZC(63,L)
        NO3fine  = ZC(64,L)
        NO3coarse= ZC(65,L)


        !// Stratospheric components included in troposphere
        M_H2   = ZC(113,L)

        !// SOA chemistry: set concentrations
        if (LSOA) then
           !// The fast conversions to SOA components may produce very
           !// small values (~1.d-100), which may cause problems in transport
           !// due to the moments. Solution is to treat the SOA precursors
           !// as zero when very small.
           M_Apine     = ZC(280,L)
           if (M_Apine .lt. 1.e-20_r8) M_Apine = 0._r8
           M_Bpine     = ZC(281,L)
           if (M_Bpine .lt. 1.e-20_r8) M_Bpine = 0._r8
           M_Limon     = ZC(282,L)
           if (M_Limon .lt. 1.e-20_r8) M_Limon = 0._r8
           M_Myrcene   = ZC(283,L)
           if (M_Myrcene .lt. 1.e-20_r8) M_Myrcene = 0._r8
           M_Sabine    = ZC(284,L)
           if (M_Sabine .lt. 1.e-20_r8) M_Sabine = 0._r8
           M_D3carene  = ZC(285,L)
           if (M_D3carene .lt. 1.e-20_r8) M_D3carene = 0._r8
           M_Ocimene   = ZC(286,L)
           if (M_Ocimene .lt. 1.e-20_r8) M_Ocimene = 0._r8
           M_Trpolene  = ZC(287,L)
           if (M_Trpolene .lt. 1.e-20_r8) M_Trpolene = 0._r8
           M_Trpinene  = ZC(288,L)
           if (M_Trpinene .lt. 1.e-20_r8) M_Trpinene = 0._r8
           M_TrpAlc    = ZC(289,L)
           if (M_TrpAlc .lt. 1.e-20_r8) M_TrpAlc = 0._r8
           M_Sestrp    = ZC(290,L)
           if (M_Sestrp .lt. 1.e-20_r8) M_Sestrp = 0._r8
           M_Trp_Ket   = ZC(291,L)
           if (M_Trp_Ket .lt. 1.e-20_r8) M_Trp_Ket = 0._r8
           !// Set C6HXR_SOA and zero the non-used C6HXR
           M_C6HXR_SOA  = ZC(12,L)
           M_C6HXR  = 0._r8 !// MUST be zero when doing SOA
           M_Tolmatic = ZC(192,L)
           M_Benzene  = ZC(193,L)
        end if !// if (LSOA) then


        M_CH3O2H = ZC(16,L) - 0.5_r8*(M_OH + M_CH3O2) 
        if (M_CH3O2H .lt. 0._r8) then
           M_CH3O2H = 1.e5_r8
        end if

        M_O3NO     = M_O3 - M_NO


        !// PAN = CH3COO2NO2
        if (M_PANX .gt. 0._r8) then
           if (M_CH3X .lt. M_PANX .and. M_CH3X .gt. 0._r8) then
              M_PAN = M_PANX - M_CH3X
           else
              M_PAN = 100._r8
           end if
        else
           if (M_CH3X .gt. 0._r8) then
              M_PANX = min(1.e5_r8, M_CH3X)
              M_PAN  = M_PANX - M_CH3X
           else
              M_PANX = 1.e5_r8
              M_CH3X = 1.e4_r8
              M_PAN  = M_PANX - M_CH3X
           end if
        end if

        !// SUM of NOx
        M_NOX   = M_NO + M_NO2 + M_NO3 + 2._r8*M_N2O5 + M_HO2NO2 + M_PAN
  
        !// Some components are SUMs
        XNOX    = M_NO + M_NO2
        XHO2NO2 = M_NO2 + M_HO2NO2                            

        !// set NOZ
        !// August 2009: NO3 and N2O5 are transported; set NOZ from
        !// these and adjust NO3 and N2O5 based on NOZ at the end.
        M_NOZ = M_NO3 + M_N2O5
        !// NOZ should therefore never be zero. Check just in case,
        !// using the old adjustment (also done at the end)
        if (M_NOZ .eq. 0._r8) then
           M_NOZ = 100._r8
           !// set N2O5 to 99% and NO3 to be 1% of NOZ
           M_N2O5 = 0.99_r8*M_NOZ
           M_NO3  = M_NOZ - M_N2O5
        end if

        !// Initialize chemical integration parameters
        DTCH = DTCH1
        ST   = ST1
        QLIN = QLIN1
        QTEST= QTEST1

        !// ----------------------------------------------------------------
        !// -------- CHEMISTRY CALCULATIONS STARTS -------------------------
        !// ----------------------------------------------------------------
        !// The chemical equations are integrated in time with use of the
        !// QSSA method is described in:
        !// Hesstvedt, E., Hov, O., Isaksen, I.S.A. (1978) Quasi steady
        !// state approximation in air pollution modelling: Comparison
        !// of two numerical schemes for oxidant prediction.
        !// Int. J. chem. Kinet., 10, 971-994.
        !// ----------------------------------------------------------------


        !//..NOZ (NO3+N2O5)-------------------------------------------------
        PROD = &
             k_oh_hno3 * M_HNO3 * M_OH   &!HNO3 + OH -> NO3 + H2O
             + k_o3_no2 * M_O3 * M_NO2   &!NO2 + O3 -> NO3 + O2
             + POLLX(41)                 &!EMISX of NO3
             + POLLX(42)                 &!EMISX of N2O5
             + k_op_no2_m * M_O3P * M_NO2 !O3P + NO2 + M -> NO3
        !// Multiply loss rates with NO3 or N2O5 and sum up; then
        !// divide by (NO3+N2O5)
        LOSS = &
             ( k_no_no3 * M_NO &
               + k_no2_no3_b * M_NO2 &
               + DNO3 &
               + k_no3_ch3cho * M_CH3CHO &
               + k_no3_dms * M_DMS &
             ) * M_NO3 &
             + k_n2o5_h2o_aer * M_N2O5
        !// SOA Secondary organic aerosols
        if (LSOA) LOSS = LOSS &
             + ( k_no3_soaC1 * (M_Apine + M_Bpine + M_Sabine + &
                           M_D3Carene + M_Trp_Ket) &
                 + k_no3_soaC2 * M_Limon &
                 + k_no3_soaC3 * (M_Trpolene + M_Trpinene) &
                 + k_no3_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
                 + k_no3_soaC5 * M_Sestrp &
                 + k_no3_soaC7 * M_C6HXR_SOA &
                 + k_no3_soaC8 * M_tolmatic &
               ) * M_NO3
        !// Divide by NOZ to get correct units of LOSS
        LOSS = LOSS / (M_NO3 + M_N2O5)

        call QSSA(1,'NOZ',DTCH,QLIN,ST,PROD,LOSS,M_NOZ)

        !//..NO3, N2O5--v
        if (k_n2o5_m .lt. (k_no2_no3_m*M_NO2)) then
           PROD = &
                k_oh_hno3 * M_HNO3 * M_OH       &!HNO3 + OH -> NO3 + H2O
                + k_o3_no2 * M_O3 * M_NO2      &!NO2 + O3 -> NO3 + O2
                + k_n2o5_m * M_N2O5             &!N2O5 + -heat-> NO3 + NO2
                + POLLX(41)                 &!EMISX of NO3
                + k_op_no2_m * M_O3P * M_NO2  !// O3P + NO2 + M -> NO3
           LOSS = &
                k_no2_no3_m * M_NO2   &!NO3 + NO2 -M-> N2O5
                + k_no_no3 * M_NO   &!NO3 + NO -> 2NO2
                + k_no2_no3_b * M_NO2 &!NO3 + NO2 -> NO + NO2 + O2
                + DNO3           &!NO3 + hv
                + k_no3_ch3cho * M_CH3CHO &!CH3CHO + NO
                + k_no3_dms * M_DMS   !DMS + NO3
            !// SOA Secondary organic aerosols
           if (LSOA) LOSS = LOSS &
                + k_no3_soaC1 * (M_Apine + M_Bpine + M_Sabine &
                            + M_D3Carene + M_Trp_Ket) &
                + k_no3_soaC2 * M_Limon &
                + k_no3_soaC3 * (M_Trpolene + M_Trpinene) &
                + k_no3_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
                + k_no3_soaC5 * M_Sestrp &
                + k_no3_soaC7 * M_C6HXR_SOA &
                + k_no3_soaC8 * M_tolmatic
 
           call QSSA(2,'NO3 N2O5',DTCH,QLIN,ST,PROD,LOSS,M_NO3)
            
           M_N2O5 = M_NOZ - M_NO3
           if (M_N2O5 .lt. 0._r8) then
              RLIM1 = M_NOZ*0.7_r8
              RLIM2 = 100._r8
              M_N2O5  = Min(RLIM2,RLIM1)
              M_NO3   = M_NOZ - M_N2O5
           end if


        else
           PROD = &
                k_no2_no3_m * M_NO2 * M_NO3 &!NO3+NO2 -M-> N2O5
                + POLLX(42)             !EMISX of N2O5
           LOSS = &
                k_n2o5_m &
                + k_n2o5_h2o_aer

           call QSSA(3,'NO3 N2O5',DTCH,QLIN,ST,PROD,LOSS,M_N2O5)

           M_NO3 = M_NOZ - M_N2O5
           if (M_NO3 .lt. 0._r8) then
              RLIM1 = M_NOZ*0.7_r8
              RLIM2 = 100._r8
              M_NO3   = min(RLIM2,RLIM1)
              M_N2O5  = M_NOZ - M_NO3
           end if
        end if



        !//..NO (if NO < O3)------------------------------------------------
        if (M_O3NO .gt. 0._r8) then
           !// Loss of NO must be included in O3NO after HOx-loop!
           LOSS_NO = &
                k_o3_no * M_O3        &!O3 + NO -> NO2 + O2
                + k_no_ch3o2 * M_CH3O2   &!CH3O2 + NO -> prod
                + k_no_ho2 * M_HO2     &!NO + HO2 -> NO2 + OH
                + k_no_no3 * M_NO3     &!NO + NO3 -> 2NO2
                + k_no_ch3cob * M_CH3COB  &
                + k_no_ch3cod * M_CH3COD  &
                + k_no_ch3xx * M_CH3XX   &
                + k_no_ch3x * M_CH3X    &
                + k_no_c2h5o2 * M_C2H5O2  &
                + k_no_c3h7o2 * M_C3H7O2  &
                + k_no_c4h9o2 * M_C4H9O2  &
                + k_no_c6h13o2 * M_C6H13O2 &
                + k_no_ar1 * M_AR1     &
                + k_no_ar3 * M_AR3     &
                + k_no_isor1 * M_ISOR1   &
                + k_no_isor2 * M_ISOR2 * 1.5_r8 &
                + k_no_ho2_b * M_HO2    &!NO + HO2 -> HNO3
                + VDEP_L(43)        &!drydep
                + k_op_no_m * M_O3P !O3P + NO + M -> NO2

           !//..NO2---------------------------------------------------------
           LOSS_NO2 = &
                DNO2                &!NO2 + hv -> NO + O3P
                + k_op_no2 * M_O3P     &
                + k_no2_oh_m * M_OH      &
                + k_no2_no3_m * M_NO3    &
                + k_no2_ch3x_m * M_CH3X    &
                + k_ho2_no2_m * M_HO2     &
                + k_o3_no2 * M_O3      &
                + VDEP_L(44)        &!drydep
                + RR_NO2_SOOT(L)    &!NO2 uptake on soot
                + k_op_no2_m * M_O3P !// O3P + NO2 + M -> NO3

           !// For stability, production of NO2 from NO is treated as
           !// production rate times (NO+NO2) and also added to the
           !// NO2 loss.
           TMP_PL = &
                k_o3_no * M_O3 &
                + k_no_ho2 * M_HO2 &
                + 2._r8 * k_no_no3 * M_NO3 &
                + k_no_ch3xx * M_CH3XX &
                + k_no_ch3x * M_CH3X &
                + k_no_ar1 * M_AR1 &
                + k_no_ar3 * M_AR3 &
                + k_no_isor1 * M_ISOR1 * fa_no_isor1 &
                + k_no_isor2 * M_ISOR2 * 1.5_r8  &
                + k_no_ch3o2 * M_CH3O2 * fa_no_ch3o2 &
                + k_no_ch3cob * M_CH3COB * fa_no_ch3cob &
                + k_no_ch3cod * M_CH3COD * fa_no_ch3cod &
                + k_no_c2h5o2 * M_C2H5O2 * fa_no_c2h5o2 &
                + k_no_c3h7o2 * M_C3H7O2 * fa_no_c3h7o2 &
                + k_no_c4h9o2 * M_C4H9O2 * fa_no_c4h9o2 &
                + k_no_c6h13o2 * M_C6H13O2 * fa_no_c6h13o2


           !// Similarly, there are two other scalings
           TMP_PL2 = &
                k_pan_m               &!To scale PAN+NO2 in PROD
                + k_oh_ho2no2 * M_OH  &!To scale XHO2NO2 in PROD
                + k_ho2no2_m + DHO2NO2 !To scale XHO2NO2 in PROD

           PROD = &
                POLLX(44)              &!Emis NO2
                + k_n2o5_m * M_N2O5        &
                + DN2O5 * M_N2O5       &
                + DNO3 * M_NO3         &!NO3 + hv -> NO2 + O3P
                + k_oh_pan * M_OH * M_PAN &
                + DHNO3 * M_HNO3       &
                + k_op_no_m * M_O3P * M_NO &!O3P + NO + M -> NO2
                !// Add some terms for stability scaling
                + TMP_PL * XNOX            &!Add as loss of NO2
                + k_pan_m * (M_PAN + M_NO2)    &!Add as loss of NO2
                + (k_ho2no2_m + DHO2NO2 + k_oh_ho2no2*M_OH) * XHO2NO2 !Add as loss of NO2

           LOSS = &
                LOSS_NO2  &
                + TMP_PL  &!To scale XNOX in PROD
                + TMP_PL2  !To scale PAN+NO2 and XHO2NO2 in PROD


           if (L .eq. 1) DDDIAG(44) = DDDIAG(44) + VDEP_L(44) * M_NO2 * DTCH

           call QSSA(4,'NO2',DTCH,QLIN,ST,PROD,LOSS,ZC(44,L))

           !// Use the new value in following integration if lifetime is short
           if (LOSS_NO2 .gt. QTEST) M_NO2 = ZC(44,L)



           !//..NO----------------------------------------------------------
           !// FLY_NOX and LIGHTNING are included in POLLX
           !// Treat production from NO2 as rates times (NO+NO2) and
           !// add as loss of NO:
           PROD = &
                POLLX(43) &
                + ( DNO2 &
                    + k_op_no2 * M_O3P &
                    + k_no2_no3_b * M_NO3 &
                  ) * XNOX

           LOSS = &
                LOSS_NO &
                + DNO2 &
                + k_op_no2 * M_O3P &
                + k_no2_no3_b * M_NO3

           if (L .eq. 1) DDDIAG(43) = DDDIAG(43) + VDEP_L(43) * M_NO * DTCH

           call QSSA(5,'NO',DTCH,QLIN,ST,PROD,LOSS,ZC(43,L))


           !// Use the new value in following integration if lifetime is short
           if (LOSS_NO .gt. QTEST) M_NO = ZC(43,L)

           !// Do not allow negative NO or NO2<RLIM1
           RLIM1 = 1.e5_r8
           if (ZC(43,L) .lt. 0._r8)  ZC(43,L) = RLIM1
           if (ZC(44,L) .lt. RLIM1) ZC(44,L) = RLIM1

           !//..O3 (when NO < O3) -------------------------------------------
           ZC(1,L) = M_O3NO + ZC(43,L)

        else !// i.e. (NO .gt. O3)

           !// WARNING -- unsual condition
           !// In this case we do not calculate NO2 explicitly, but
           !// calculate it based on NOX at the end of chemistry.
           LOSS_NO  = 0._r8
           LOSS_NO2 = 0._r8

           !//..O3 (when NO > O3) ------------------------------------------
           PROD = DNO2 * M_NO2
           LOSS = &
                k_o3_no * M_NO &
                + k_o3_no2 * M_NO2 &
                + k_o3_ho2 * M_HO2 &
                + k_o3_c2h4 * M_C2H4 &
                + ( k_od_h2o * M_H2O * M_O1D &
                    + k_op_no2 * M_NO2 * M_O3P &
                  ) / M_O3 &
                + k_o3_c3h6 * M_C3H6 &
                + k_o3_oh * M_OH &
                + VDEP_L(1)
           if (LSULPHUR) LOSS = LOSS + CAQ0172 * M_SO2
           !// SOA Secondary organic aerosols
           if (LSOA) LOSS = LOSS &
                + k_o3_soaC1 * (M_Apine + M_Bpine + M_Sabine &
                           + M_D3Carene + M_Trp_Ket) &
                + k_o3_soaC2 * M_Limon &
                + k_o3_soaC3 * (M_Trpolene + M_Trpinene) &
                + k_o3_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
                + k_o3_soaC5 * M_Sestrp &
                + k_o3_soaC7 * M_C6HXR_SOA &
                + k_o3_soaC8 * M_tolmatic

           call QSSA(6,'O3',DTCH,QLIN,ST,PROD,LOSS,ZC(1,L))

           !//..NO----------------------------------------------------------
           ZC(43,L) = ZC(1,L) - M_O3NO

        end if !// if (M_O3NO > 0._r8) then

        !//..O'D------------------------------------------------------------
        PROD = DBO3 * M_O3
        LOSS = k_od_m * AIR_MOLEC(L) &
               + k_od_h2o * M_H2O                      

        M_O1D   = PROD / LOSS
        ZC(39,L) = M_O1D


        !// ----------------------------------------------------------------
        !// OH, HO2, AND RO2 CALCULATED BY ITERATIONS
        !// ----------------------------------------------------------------
        DTCH = DTCH2
        ST   = ST2
        QLIN = QLIN2
        QTEST= QTEST2

        !// Variables not initialized: ARAD, CH3O and CHO
        !// ----------------------------------------------------------------

        !//..Radicals
        M_ARAD = M_CH3XX  + M_C6H13O2 + M_C4H9O2 + M_CH3COB + M_C2H5O2 &
                 + M_C3H7O2 + M_CH3COD  + M_ISOR1  + M_ISOR2


        !//..CH3O-----------------------------------------------------------
        PROD = &
             ( 0.8_r8 * k_ch3o2_ch3o2 * M_CH3O2        &!Produced from CH3O2
               + k_no_ch3o2 * M_NO * fa_no_ch3o2    &!within this parenthesis
               + k_ch3o2_ch3xx * M_CH3XX              &
               + 0.5_r8 * (k_ch3o2_c4h9o2 * M_C4H9O2    &
                          + k_ch3o2_c2h5o2 * M_C2H5O2  &
                          + k_ch3o2_c3h7o2 * M_C3H7O2) &
               + k_ch3o2_ch3cob * M_CH3COB  &
               + k_ch3o2_c6h13o2 * M_C6H13O2 &
               + k_ch3o2_ch3x_a * M_CH3X   &
               + k_ch3o2_isor1 * M_ISOR1   &
               + k_ch3o2_isor2 * M_ISOR2   &
             ) * M_CH3O2 &
             + 0.03_r8 * k_o3_c3h6 * M_C3H6 * M_O3 &
             + DCH3O2H * M_CH3O2H
        if (.not. LOLD_H2OTREATMENT) PROD = PROD &
             + k_od_ch4_b * M_O1D * M_CH4 !O(1D) + CH4 -> (CH3O/CH2OH + H)

        LOSS = k_ch3o_o2 * M_O2

        M_CH3O = PROD / LOSS !// Steady-state for short-lived species


        !//..CHO------------------------------------------------------------
        PROD = &
             DACH2O * M_CH2O &
             + k_oh_ch2o * M_OH * M_CH2O &
             + DCH3CHO * M_CH3CHO &
             + k_oh_hcohco_m_a * M_OH * M_HCOHCO

        LOSS = k_cho_o2 * M_O2

        M_CHO = PROD / LOSS !// Steady-state for short-lived species



        !// ----------------------------------------------------------------
        !// HOx-LOOP (Time step DTCH2)
        !// ----------------------------------------------------------------
        do IIJ = 1, OH_HO2_ITER

          PROD_OH = &
               2._r8 * (DH2O2 * M_H2O2               &!H2O2 + hv -> 2OH
                        + k_od_h2o * M_O1D * M_H2O)  &!H2O + O1D -> 2OH
               + 0.15_r8 * k_o3_c3h6 * M_O3 * M_C3H6 &!C3H6 + O3 -> 0.15*OH
               + DCH3O2H * M_CH3O2H                !CH3O2H + hv -> OH + CH3O
          if (.not. LOLD_H2OTREATMENT) PROD_OH = PROD_OH &
               + k_od_ch4_a * M_O1D * M_CH4        &!O(1D) + CH4 -> OH + CH3
               + k_od_h2 * M_O1D * M_H2             !O(1D) + H2  -> OH + H

          LOSS_OH = &
               k_no2_oh_m * M_NO2       &!OH + NO2 -M-> HNO3
               + k_oh_c2h4_m * M_C2H4    &!OH + C2H4 -> (HOCH2CH2) --> CH3 + HCHO
               + k_oh_ch2o * M_CH2O    &!OH + HCHO -> 
               + (k_oh_co_a + k_oh_co_b) * M_CO      &!OH + CO   -> HO2 (CO2)
               + k_oh_ch4 * M_CH4     &!OH + CH4  -> CH3 + H2O
               + k_oh_h2o2 * M_H2O2    &!OH + H2O2 -> HO2 + H2O
               + k_o3_oh * M_O3      &!OH + O3   -> HO2 + O2
               + k_oh_hno3 * M_HNO3    &!OH + HNO3 -> NO3 + H2O
               + k_oh_ch3cox * M_CH3COX  &!OH + CH3COX -> CH3COB
               + k_oh_c3h6_m * M_C3H6    &!OH + C3H6 + O2 -> CH3XX
               + k_oh_ch3cho * M_CH3CHO  &!OH + CH3CHO -> CH3CO + H2O
               + k_oh_c4h10 * M_C4H10   &!OH + C4H10 + O2 -> C4H9O2 + H2O
               + k_oh_c6h14 * M_C6H14 &!OH + C6H14 + O2 -> C6H13O2 + H2O
               + k_oh_c6hxr * M_C6HXR &!OH + C6HXR -> AR1 (=0 if SOA included)
               + k_oh_c2h6 * M_C2H6    &!OH + C2H6 + O2 -> C2H5O2 + H2O
               + k_oh_ar2 * M_AR2     &!OH + AR2 -> AR3
               + k_oh_isok * M_ISOK    &!OH + ISOK -> ISOR2
               + k_oh_c3h8 * M_C3H8    &!OH + C3H8 + O2 -> C3H7O2 + H2O
               + k_oh_isoprene * M_ISOPREN &!OH + ISOPREN -> ISOR1
               !+ k_oh_isoprene * M_ISOPREN_TOT &!OH + ISOPREN -> products (44% ISOR1)
               + (k_oh_hcohco_m_a &
                  + k_oh_hcohco_m_b &
                  + k_oh_hcohco_m_c) * M_HCOHCO &!OH + HCOHCO
               + k_oh_rcohco * M_RCOHCO  &!OH + 
               + k_oh_ch3o2h_a * M_CH3O2H &!OH + 
               + k_oh_ho2no2 * M_HO2NO2  &!OH + HO2NO2 -> NO2 + H2O + O2
               + k_oh_pan * M_PAN     &!OH + PAN ->
               + k_oh_h2 * M_H2      &!OH + H2 -> H2O + H    c121205
               + 2._r8 * k_oh_oh_m * M_OH &!OH + OH + M -> H2O2 + M
               !// Extra terms for stability: Must be added to production below
               + k_no_ho2 * M_NO       !HO2 + NO -> OH + NO2 FOR STABILITY!

          !// Sulphur reactions
          if (LSULPHUR) LOSS_OH = LOSS_OH &
               + k_oh_dms_a * M_DMS  &!OH + DMS -> H2O + CH3SCH2
               + C4071b * M_DMS  &!OH + DMS -> 0.75 SO2 + 0.25 MSA
               + k_oh_h2s * M_H2S   &!OH + H2S + O2 -> SO2 + H2O
               + CTOT4072 * M_SO2 !OH + SO2 -> H2SO4

          !// SOA Secondary organic aerosols
          if (LSOA) LOSS_OH = LOSS_OH &
               + k_oh_soaC1 * (M_Apine + M_Bpine + M_Sabine &
                          + M_D3Carene + M_Trp_Ket) &
               + k_oh_soaC2 * M_Limon &
               + k_oh_soaC3 * (M_Trpolene + M_Trpinene) &
               + k_oh_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
               + k_oh_soaC5 * M_Sestrp &
               + k_oh_soaC7 * M_C6HXR_SOA &
               + k_oh_soaC8 * M_tolmatic &
               + k_oh_benzene * M_Benzene

          PROD_HO2 = &
               DACH2O * M_CH2O                &!HCHO + hv -> H + HCO -> HO2
               + k_ch3o_o2 * M_CH3O * M_O2       &!CH3O + O2 -> HO2 + HCHO
               + k_cho_o2 * M_CHO * M_O2         &!CHO + O2 -> HO2 + CO
               + 0.25_r8 * k_o3_c3h6 * M_C3H6 * M_O3   &!C3H6 + O3 -> 4HO2 ...
               + k_no_c4h9o2 * M_C4H9O2 * M_NO &
                       * fa_no_c4h9o2 * fb_no_c4h9o2          &!C4H9O2 + NO -> HO2 + ...
               + k_no_ch3cob * M_CH3COB * M_NO * fa_no_ch3cob &!CH3COB + NO -> HO2 + ...
               + k_no_ch3cod * M_CH3COD * M_NO * fa_no_ch3cod &!CH3COD + NO -> HO2 + ...
               + k_no_c2h5o2 * M_C2H5O2 * M_NO &
                       * fa_no_c2h5o2 * fb_no_c2h5o2          &!C2H5O2 + NO -> HO2 + ...
               + k_no_c3h7o2 * M_C3H7O2 * M_NO &
                       * fa_no_c3h7o2 * fb_no_c3h7o2 &!C3H7O2 + NO -> HO2 + ...
               + 0.12_r8 * k_o3_c2h4 * M_O3 * M_C2H4   &!O3+C2H4 -> HO2+HCHO+CH3O2
               + k_ch3o2_c4h9o2 * M_C4H9O2 * M_CH3O2 * 0.5_r8 * fb_no_c4h9o2 &
               + k_no_ch3xx * M_NO * M_CH3XX &
               + 0.5_r8 * fb_no_c2h5o2 * k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
               + k_ch3o2_ch3cob * M_CH3O2 * M_CH3COB &
               + k_ho2no2_m * M_HO2NO2            &!HO2NO2 -M-> HO2 + NO2
               + k_ch3o2_ch3xx * M_CH3O2 * M_CH3XX &
               + k_no_ar1 * M_NO * M_AR1 &
               + k_no_ar3 * M_NO * M_AR3 &
               + k_no_isor1 * M_NO * M_ISOR1 * fa_no_isor1 &
               + k_no_isor2 * M_NO * M_ISOR2 * 0.67_r8 &
               + (k_ch3o2_isor1 * M_ISOR1 &
                  + k_ch3o2_isor2 * M_ISOR2) * M_CH3O2 &
               + DHO2NO2 * M_HO2NO2         !HO2NO2 + hv -> HO2 + NO2

          LOSS_HO2 = &
               k_no_ho2 * M_NO &
               + k_ho2_no2_m * M_NO2 &
               + k_o3_ho2 * M_O3     &!O3 + HO2 -> OH + 2*O2
               + k_ho2_ch3o2 * M_CH3O2  &!HO2 + CH3O2 -> CH3O2H + (O2)
               + k_ho2_ch3x * M_CH3X   &!HO2 + CH3X --> CH3O2H
               + k_ho2_radical * M_ARAD   &!HO2 + RADICAL --> CH3O2H
               + RR_HO2_AER(L)    &!QAER(L) (uptake on aerosols)
               + k_no_ho2_b * M_NO    &!NO + HO2 -> HNO3
               !// Extra terms for stability: Must be added to production below
               + k_ho2no2_m             &!HO2NO2 -M-> HO2 + NO2
               + DHO2NO2           !HO2NO2 + hv -> HO2 + NO2



          !// OH/HO2 interation
          !// --------------------------------------------------------------

          !// Help constants for OH HO2 iteration

          !// OH produced by HO2
          RLIM1 = k_o3_ho2 * M_O3   &!O3 + HO2 -> OH + 2*O2
                  + k_no_ho2 * M_NO  !NO + HO2 -> OH + NO2

          !// HO2 produced by OH
          RLIM2 = &
               (k_oh_co_a + k_oh_co_b) * M_CO        &!OH + CO   -> HO2 (CO2)
               + k_oh_h2o2 * M_H2O2    &!OH + H2O2 -> HO2 + H2O
               + k_o3_oh * M_O3      &!OH + O3   -> HO2 + O2
               + k_oh_hcohco_m_c * M_HCOHCO &!OH + HCOHCO -O2-> 2CO + HO2
               + k_oh_ho2no2 * M_HO2NO2  &!HO2NO2 + OH -> NO2 + H2O + O2
               + k_oh_h2 * M_H2       !OH + H2 -> H2O + H    c121205

          !// Set separate OH and HO2 for iteration
          OH_NEW = M_OH
          HO2_NEW = M_HO2
          do IJL = 1, ITER_OH_MAX

            OH_OLD = OH_NEW          !// OH at last iteration
            HO2_OLD= HO2_NEW         !// HO2 at last iteration

            !//..OH---------------------------------------------------------
            PROD = &
                 PROD_OH                &!General production terms
                 + RLIM1 * HO2_OLD      &!Production due to HO2
                 + k_no_ho2 * M_NO * OH_OLD !Added to LOSS_OH above (for stability)
            LOSS = &
                 LOSS_OH           &!General loss terms
                 + k_oh_ho2 * HO2_OLD  !OH + HO2 -> H2O + O2
            !// Set new OH
            OH_NEW = PROD / LOSS



            !//..HO2--------------------------------------------------------
            PROD = &
                 PROD_HO2             &!General production terms
                 + RLIM2*OH_NEW       &!Production due to OH
                 + (k_ho2no2_m              &!Added to LOSS_HO2 above (stability)
                 + DHO2NO2) * HO2_OLD  !Added to LOSS_HO2 above (stability)
            LOSS = &
                 LOSS_HO2                 &!General loss terms
                 + k_oh_ho2 * OH_NEW      &!OH + HO2 -> H2O + O2
                 + 2._r8 * k_ho2_ho2_tot * HO2_OLD !HO2 + HO2 -> H2O2 + O2 (+/- M)
            !// Set new HO2
            HO2_NEW = PROD / LOSS


            !// Check for stable OH
            if (abs((OH_NEW - OH_OLD)/OH_OLD) .lt. limit_oh_ho2) exit

          end do !// do IJL = 1, ITER_OH_MAX

          !// Update OH and HO2---------------------------------------------
          M_OH = OH_NEW
          M_HO2 = HO2_NEW

          !// Did we get stable OH?
          if (IJL .gt. ITER_OH_MAX) write(6,'(a,3i4,2es12.5,es16.8)') &
               f90file//':'//subr//'OH not stable',ICOL,JCOL,L,M_OH,M_HO2, &
               abs((OH_NEW - OH_OLD)/OH_OLD)
          if (M_OH .ne. M_OH) then
             write(6,'(a,3i4,2es12.5,es16.8)') &
               f90file//':'//subr//': OH is NaN: ',ICOL,JCOL,L,M_OH,M_HO2, &
               abs((OH_NEW - OH_OLD)/OH_OLD)
             stop
          end if

          !//..O(3P)--------------------------------------------------------
          PROD = &
               DO3 * M_O3 &
               + DNO2 * M_NO2 &
               + DNO3 * M_NO3 &
               + 2._r8 * DO2 * M_O2
          LOSS = &
               k_op_o2_m * M_O2 &
               + k_op_no2 * M_NO2 &
               + k_op_no_m * M_NO   &!O3P + NO + M -> NO2
               + k_op_no2_m * M_NO2   !O3P + NO2 + M -> NO3

          M_O3P = PROD / LOSS


          !//..CH3CH(O2)CH2OH-----------------------------------------------
          PROD = k_oh_c3h6_m * M_OH * M_C3H6
          LOSS = &
               k_no_ch3xx * M_NO &
               + QAER(L) &
               + k_ch3o2_ch3xx * M_CH3O2 &
               + k_ho2_radical * M_HO2    !HO2 + RADICAL --> CH3O2H

          call QSSA(7,'CH3CH..',DTCH,QLIN,ST,PROD,LOSS,ZC(28,L))

          if (LOSS .gt. QTEST) M_CH3XX = ZC(28,L)


          !//..CH3CO---------------------------------------------------------
          PROD = &
               k_oh_ch3cho * M_OH * M_CH3CHO &
               + k_no_isor2 * M_NO * M_ISOR2 * 0.33_r8 &
               + DCH3COX * M_CH3COX &
               + 2._r8 * DCH3COY * M_CH3COY &
               + k_oh_hcohco_m_b * M_OH * M_HCOHCO &
               + k_oh_rcohco * M_OH * M_RCOHCO &
               + k_no3_ch3cho * M_NO3 * M_CH3CHO &
               + DACETON_A * M_ACETON

          LOSS = k_ch3co_o2 * M_O2

          M_CH3CO = PROD / LOSS


          !//..CH3COO2-------------------------------------------------------
          PROD = &
               k_ch3co_o2 * M_O2 * M_CH3CO &
               + k_pan_m * M_PANX
          LOSS = &
               k_no_ch3x * M_NO &
               + k_ch3o2_ch3x_a * M_CH3O2 &
               + k_ch3o2_ch3x_b * M_CH3O2 &
               + k_no2_ch3x_m * M_NO2 &
               + k_pan_m &
               + k_ho2_ch3x * M_HO2 &
               + 2._r8 * k_ch3x_ch3x * M_CH3X &
               + QAER(L)

          call QSSA(8,'CH3COO2',DTCH,QLIN,ST,PROD,LOSS,ZC(37,L))

          RCH3XP = M_CH3X

          if (LOSS .gt. QTEST) M_CH3X = ZC(37,L)


          !//..PANX----------------------------------------------------------
          PROD = &
               k_ch3co_o2 * M_O2 * M_CH3CO &
               + POLLX(5)              !EMISX of PANX
          LOSS = &
               ( ( VDEP_L(5) &
                   + k_oh_pan * M_OH &
                 ) * M_PAN &
                 + ( k_no_ch3x * M_NO &
                     + k_ch3o2_ch3x_a * M_CH3O2 &
                     + k_ch3o2_ch3x_b * M_CH3O2 &
                     + k_ho2_ch3x * M_HO2 &
                     + 2._r8 * k_ch3x_ch3x * M_CH3X &
                     + QAER(L) &
                   ) * RCH3XP &
               ) / M_PANX

          if (L .eq. 1) DDDIAG(5) = DDDIAG(5) + VDEP_L(5) * M_PAN * DTCH

          call QSSA(9,'PAN',DTCH,QLIN,ST,PROD,LOSS,ZC(5,L))

          M_PANX = ZC(5,L)
          M_PAN  = max(ZC(5,L) - ZC(37,L), 0._r8)



          !//..SEC-C6H13O2---------------------------------------------------
          PROD = k_oh_c6h14 * M_OH * M_C6H14
          LOSS = &
               k_no_c6h13o2 * M_NO &
               + QAER(L) &
               + k_ch3o2_c6h13o2 * M_CH3O2 &
               + k_ho2_radical * M_HO2    !HO2 + RADICAL --> CH3O2H

          call QSSA(10,'SEC-1..',DTCH,QLIN,ST,PROD,LOSS,ZC(25,L))

          if (LOSS .gt. QTEST) M_C6H13O2 = ZC(25,L)


          !//..SEC-C4H9O2----------------------------------------------------
          PROD = &
               k_oh_c4h10 * M_OH * M_C4H10 &
               + k_no_c6h13o2 * M_C6H13O2 * M_NO * fa_no_c6h13o2 &
               + k_ch3o2_c6h13o2 * M_CH3O2 * M_C6H13O2

          LOSS = &
               k_no_c4h9o2 * M_NO &
               + k_ch3o2_c4h9o2 * M_CH3O2 &
               + k_ho2_radical * M_HO2   &!HO2 + RADICAL --> CH3O2H
               + QAER(L)
            
          call QSSA(11,'SEC-2..',DTCH,QLIN,ST,PROD,LOSS,ZC(24,L))
            
          if (LOSS .gt. QTEST) M_C4H9O2 = ZC(24,L)


          !//..CH3COCH(O2)CH3------------------------------------------------
          PROD = k_oh_ch3cox * M_OH * M_CH3COX
          LOSS = &
               k_no_ch3cob * M_NO &
               + k_ch3o2_ch3cob * M_CH3O2 &
               + k_ho2_radical * M_HO2   &!HO2 + RADICAL --> CH3O2H
               + QAER(L)

          call QSSA(12,'CH3CO...',DTCH,QLIN,ST,PROD,LOSS,ZC(27,L))

          if (LOSS .gt. QTEST) M_CH3COB = ZC(27,L)


          !//..C2H5O2--------------------------------------------------------
          PROD = &
               DCH3COX * M_CH3COX &
               + k_no_c4h9o2 * M_C4H9O2 * M_NO &
                             * fa_no_c4h9o2 * (1._r8 - fb_no_c4h9o2) &
               + k_ch3o2_c4h9o2 * M_C4H9O2 * M_CH3O2 &
                                * (0.5_r8 * (1._r8 - fb_no_c4h9o2)) &
               + k_oh_c2h6 * M_OH * M_C2H6

          LOSS = &
               QAER(L) &
               + k_ch3o2_c2h5o2 * M_CH3O2 &
               + k_no_c2h5o2 * M_NO    &
               + k_ho2_radical * M_HO2    !HO2 + RADICAL --> CH3O2H

          call QSSA(13,'C2H5O2',DTCH,QLIN,ST,PROD,LOSS,ZC(23,L))

          if (LOSS .gt. QTEST) M_C2H5O2 = ZC(23,L)



          !//..C3H7O2--------------------------------------------------------
          PROD = k_oh_c3h8 * M_OH * M_C3H8
          LOSS = &
               QAER(L) &
               + k_ch3o2_c3h7o2 * M_CH3O2 &
               + k_no_c3h7o2 * M_NO    &
               + k_ho2_radical * M_HO2    !HO2 + RADICAL --> CH3O2H
            
          call QSSA(14,'C2H7O2',DTCH,QLIN,ST,PROD,LOSS,ZC(49,L))

          if (LOSS .gt. QTEST) M_C3H7O2 = ZC(49,L)


          !//..CH3COD (CH3COCH2(O2))-----------------------------------------
          PROD = k_oh_aceton * M_OH * M_ACETON
          LOSS = &
               QAER(L) &
               + k_no_ch3cod * M_NO &
               + k_ch3o2_ch3cod * M_CH3O2 &
               + k_ho2_radical * M_HO2    !HO2 + RADICAL --> CH3O2H

          call QSSA(15,'CH3COD',DTCH,QLIN,ST,PROD,LOSS,ZC(51,L))

          if (LOSS .gt. QTEST) M_CH3COD = ZC(51,L)


          !//..aromatic products 1-------------------------------------------
          if (LSOA) then
             !// SOA Secondary organic aerosols (use C6HXR_SOA)
             PROD = &
                  ( k_oh_benzene * M_Benzene &
                    + k_oh_soaC8 * M_tolmatic &
                    + k_oh_soaC7 * M_C6HXR_SOA &
                  ) * M_OH &
                  + ( k_o3_soaC8 * M_tolmatic &
                     + k_o3_soaC7 * M_C6HXR_SOA &
                    ) * M_O3 &
                  + ( k_no3_soaC8 * M_tolmatic &
                      + k_no3_soaC7 * M_C6HXR_SOA &
                    ) * M_NO3
          else
             !// Regular chemistry C6HXR
             PROD = 0.6_r8 * k_oh_c6hxr * M_OH * M_C6HXR
          end if
          LOSS = &
               k_no_ar1 * M_NO &
               + QAER(L)

          call QSSA(15,'AROMAT1',DTCH,QLIN,ST,PROD,LOSS,ZC(29,L))

          if (LOSS .gt. QTEST) M_AR1 = ZC(29,L)


          !//..aromatic products 2-------------------------------------------
          PROD = k_no_ar1 * M_NO * M_AR1
          LOSS = k_oh_ar2 * M_OH

          call QSSA(16,'AROMAT2',DTCH,QLIN,ST,PROD,LOSS,ZC(30,L))

          if (LOSS .gt. QTEST) M_AR2 = ZC(30,L)


          !//..aromatic products 3-------------------------------------------
          PROD = k_oh_ar2 * M_OH * M_AR2
          LOSS = &
               k_no_ar3 * M_NO &
               + QAER(L)

          call QSSA(16,'AROMAT3',DTCH,QLIN,ST,PROD,LOSS,ZC(31,L))

          if (LOSS .gt. QTEST) M_AR3 = ZC(31,L)


          !//..Peroxy radicals (RO2) from ISOPREN + OH-----------------------
          if (LSOA) then
             !// When SOA included, some of the products go elsewhere
             !PROD = k_oh_isoprene * M_ISOPREN_GAS * M_OH &
             PROD = k_oh_isoprene * M_ISOPREN * M_OH &
                         * (1._r8 - f_oh_isoprene_sg61 - f_oh_isoprene_sg62)
          else
             !PROD = k_oh_isoprene * M_ISOPREN_GAS * M_OH
             PROD = k_oh_isoprene * M_ISOPREN * M_OH
          end if
          LOSS = &
               k_no_isor1 * M_NO &
               + k_ch3o2_isor1 * M_CH3O2 &
               + k_ho2_radical * M_HO2   &!HO2 + RADICAL --> CH3O2H
               + QAER(L)

          call QSSA(17,'RO2-ISOP',DTCH,QLIN,ST,PROD,LOSS,ZC(32,L))

          if (LOSS .gt. QTEST) M_ISOR1 = ZC(32,L)


          !//..Peroxy radicals from ISOK + OH--------------------------------
          PROD = k_oh_isok * M_OH * M_ISOK
          LOSS = &
               k_no_isor2 * M_NO &
               + k_ch3o2_isor2 * M_CH3O2 &
               + k_ho2_radical * M_HO2   &!HO2 + RADICAL --> CH3O2H
               + QAER(L)

          call QSSA(18,'RO2-ISOK',DTCH,QLIN,ST,PROD,LOSS,ZC(34,L))

          if (LOSS .gt. QTEST) M_ISOR2 = ZC(34,L)


          !//..ISOK----------------------------------------------------------
          !// Sum of methylvinylketone and methacrolein (products of isoprene)
          PROD = &
               k_no_isor1 * M_ISOR1 * M_NO * fa_no_isor1 &
               + k_ch3o2_isor1 * M_CH3O2 * M_ISOR1
          LOSS = k_oh_isok * M_OH

          call QSSA(19,'ISOK',DTCH,QLIN,ST,PROD,LOSS,ZC(33,L))

          if (LOSS .gt. QTEST) M_ISOK = ZC(33,L)


          !//..Radicals revisited--------------------------------------------
          M_ARAD = M_CH3XX + M_C6H13O2 + M_C4H9O2 + M_CH3COB + M_C2H5O2 &
                 + M_C3H7O2 + M_CH3COD + M_ISOR1 + M_ISOR2

          !//..CH3-----------------------------------------------------------
          PROD = &
               k_oh_c2h4_m * M_OH * M_C2H4 &
               + k_oh_ch4 * M_OH * M_CH4 &
               + DCH3CHO * M_CH3CHO &
               + k_no_ch3x * M_NO * M_CH3X &
               + 0.31_r8 * k_o3_c3h6 * M_O3 * M_C3H6 &
               + k_no_c2h5o2 * M_C2H5O2 * M_NO &
                             * fa_no_c2h5o2 * (1._r8 - fb_no_c2h5o2) &
               + k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
                                * (1._r8 - fb_no_c2h5o2) * 0.5_r8 &
               + k_ch3o2_c3h7o2 * M_CH3O2 * M_C3H7O2 &
                                * (1._r8 - fb_no_c3h7o2) * 0.5_r8 &
               + k_no_c3h7o2 * M_C3H7O2 * M_NO &
                             * fa_no_c3h7o2 * (1._r8 - fb_no_c3h7o2) &
               + k_ch3o2_ch3x_a * M_CH3O2 * M_CH3X &
               + (DACETON_A + 2._r8 * DACETON_B) * M_ACETON
          if (.not. LOLD_H2OTREATMENT) PROD = PROD &
               + k_od_ch4_a * M_O1D * M_CH4 !O(1D) + CH4 -> OH + CH3

          LOSS = k_ch3_o2_m * M_O2

          M_CH3 = PROD / LOSS


          !//..CH3O2---------------------------------------------------------
          PROD = &
               k_ch3_o2_m * M_O2 * M_CH3 &
               + k_oh_ch3o2h_a * M_OH * M_CH3O2H 
          LOSS = &
               k_no_ch3o2 * M_NO &
               + 2._r8 * k_ch3o2_ch3o2 * M_CH3O2 &
               + k_ho2_ch3o2 * M_HO2 &
               + k_ch3o2_ch3xx * M_CH3XX &
               + k_ch3o2_c4h9o2 * M_C4H9O2 &
               + k_ch3o2_c2h5o2 * M_C2H5O2 &
               + k_ch3o2_c3h7o2 * M_C3H7O2 &
               + k_ch3o2_ch3cob * M_CH3COB &
               + k_ch3o2_c6h13o2 * M_C6H13O2 &
               + k_ch3o2_ch3x_a * M_CH3X &
               + k_ch3o2_ch3x_b * M_CH3X &
               + k_ch3o2_isor1 * M_ISOR1 &
               + k_ch3o2_isor2 * M_ISOR2 &
               + QAER(L)

          call QSSA(20,'CH3O2',DTCH,QLIN,ST,PROD,LOSS,ZC(22,L))

          if (LOSS .gt. QTEST) M_CH3O2 = ZC(22,L)


          !//..CH3O recalculated --------------------------------------------
          PROD = &
               0.8_r8 * k_ch3o2_ch3o2 * M_CH3O2 * M_CH3O2 &
               + k_no_ch3o2 * M_CH3O2 * M_NO * fa_no_ch3o2 &
               + k_ch3o2_ch3xx * M_CH3O2 * M_CH3XX &
               + ( 0.5_r8 * (k_ch3o2_c4h9o2 * M_C4H9O2 &
                            + k_ch3o2_c2h5o2 * M_C2H5O2 &
                            + k_ch3o2_c3h7o2 * M_C3H7O2) &
                   + k_ch3o2_ch3cob * M_CH3COB &
                   + k_ch3o2_c6h13o2 * M_C6H13O2 &
                 ) * M_CH3O2 &
               + 0.03_r8 * k_o3_c3h6 * M_C3H6 * M_O3 &
               + k_ch3o2_ch3x_a * M_CH3O2 * M_CH3X &
               + DCH3O2H * M_CH3O2H &
               + (k_ch3o2_isor1 * M_ISOR1 &
                  + k_ch3o2_isor2 * M_ISOR2) * M_CH3O2
          if (.not. LOLD_H2OTREATMENT) PROD = PROD &
               + k_od_ch4_b * M_O1D * M_CH4 !O(1D) + CH4 -> (CH3O/CH2OH + H)

          LOSS = k_ch3o_o2 * M_O2

          M_CH3O = PROD / LOSS


          !//..CHO recalculated ---------------------------------------------
          PROD = &
               DACH2O * M_CH2O &
               + k_oh_ch2o * M_OH * M_CH2O &
               + DCH3CHO * M_CH3CHO &
               + k_oh_hcohco_m_a * M_OH * M_HCOHCO

          LOSS = k_cho_o2 * M_O2

          M_CHO = PROD / LOSS



          !//..HO2 recalculated ---------------------------------------------
          PROD = &
               (k_oh_co_a + k_oh_co_b) * M_CO * M_OH &
               + DACH2O * M_CH2O &
               + k_ch3o_o2 * M_CH3O * M_O2 &
               + k_cho_o2 * M_CHO * M_O2 &
               + k_oh_h2o2 * M_OH * M_H2O2 &
               + k_o3_oh * M_OH * M_O3 &
               + 0.25_r8 * k_o3_c3h6 * M_C3H6 * M_O3 &
               + k_no_c4h9o2 * M_C4H9O2 * M_NO * fa_no_c4h9o2 * fb_no_c4h9o2 &
               + k_no_ch3cob * M_CH3COB * M_NO * fa_no_ch3cob &
               + k_no_ch3cod * M_CH3COD * M_NO * fa_no_ch3cod &
               + k_no_c2h5o2 * M_C2H5O2 * M_NO * fa_no_c2h5o2 * fb_no_c2h5o2 &
               + k_no_c3h7o2 * M_C3H7O2 * M_NO &
                             * fa_no_c3h7o2 * fb_no_c3h7o2 &
               + 0.12_r8 * k_o3_c2h4 * M_O3 * M_C2H4 &
               + k_ch3o2_c4h9o2 * M_C4H9O2 * M_CH3O2 * 0.5_r8 * fb_no_c4h9o2 &
               + k_ch3o2_c3h7o2 * M_C3H7O2 * M_CH3O2 * 0.5_r8 * fb_no_c3h7o2 &
               + k_no_ch3xx * M_NO * M_CH3XX &
               + 0.5_r8 * fb_no_c2h5o2 * k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
               + k_ch3o2_ch3cob * M_CH3O2 * M_CH3COB &
               + k_ho2no2_m * (M_HO2 + M_HO2NO2) &
               + k_ch3o2_ch3xx * M_CH3O2 * M_CH3XX &
               + k_no_ar1 * M_NO * M_AR1 &
               + k_no_ar3 * M_NO * M_AR3 &
               + k_no_isor1 * M_NO * M_ISOR1 * fa_no_isor1 &
               + k_no_isor2 * M_NO * M_ISOR2 * 0.67_r8 &
               + (k_ch3o2_isor1 * M_ISOR1 &
                  + k_ch3o2_isor2 * M_ISOR2) * M_CH3O2 &
               + DHO2NO2 * (M_HO2 + M_HO2NO2) &
               + k_oh_hcohco_m_c * M_OH * M_HCOHCO &
               + k_oh_ho2no2 * M_OH * M_HO2NO2 &
               + k_oh_h2 * M_H2 * M_OH !OH+H2

          LOSS = &
               k_no_ho2 * M_NO &
               + k_ho2_no2_m * M_NO2 &
               + 2._r8 * k_ho2_ho2_tot * M_HO2 &
               + k_o3_ho2 * M_O3 &
               + k_oh_ho2 * M_OH &
               + k_ho2_ch3o2 * M_CH3O2 &
               + k_ho2_radical * M_ARAD &!HO2 + RADICAL --> CH3O2H
               + k_ho2no2_m &
               + DHO2NO2 &
               + RR_HO2_AER(L) &!QAER(L)
               + k_no_ho2_b * M_NO  !NO + HO2 -> HNO3

          call QSSA(21,'HO2',DTCH,QLIN,ST,PROD,LOSS,ZC(21,L))

          if (LOSS .gt. QTEST) M_HO2 = ZC(21,L)               


          !//..HO2NO2--------------------------------------------------------
          PROD = &
               k_ho2_no2_m * M_HO2 * (M_NO2 + M_HO2NO2) &
               + POLLX(17)   !EMISX of HO2NO2
          LOSS = &
               k_ho2no2_m &
               + DHO2NO2 &
               + k_ho2_no2_m * M_HO2 &
               + k_oh_ho2no2 * M_OH
          !// Sulphur reactions
          if (LSULPHUR) LOSS = LOSS &
               + CAQ1772 * M_SO2

          call QSSA(22,'HO2NO2',DTCH,QLIN,ST,PROD,LOSS,ZC(17,L))

          if (LOSS .gt. QTEST) M_HO2NO2 = ZC(17,L)        


          !// ---------------------------------------------------------------
          !// SOA Secondary organic aerosols: SOA chemistry with shorter
          !// time step here. Still within HOx-loop!
          !// ---------------------------------------------------------------
          if (LSOA) then
             !// Class 5 loop -----------------------------------------------
             !// Time step: original (DTCH2) shortened by 15
             DTCH = DTCH2 / 15._r8
             ZTMP = 1._r8 / DTCH
             do SOA_N = 1, 15

                PROD = &
                     M_Sestrp * M_O3 * k_o3_soaC5 * SOA_STOC511 &
                              * (TMASS(trsp_idx(290))/TMASS(trsp_idx(154))) &
                     + M_Sestrp * M_OH * k_oh_soaC5 * SOA_STOC521 &
                              * (TMASS(trsp_idx(290))/TMASS(trsp_idx(154)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS51',DTCH,QLIN,ST,PROD,LOSS,ZC(154,L))

                PROD = &
                     M_Sestrp * M_O3 * k_o3_soaC5 * SOA_STOC512 &
                              * (TMASS(trsp_idx(290))/TMASS(trsp_idx(159))) &
                     + M_Sestrp * M_OH * k_oh_soaC5 * SOA_STOC522 &
                              * (TMASS(trsp_idx(290))/TMASS(trsp_idx(159)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS52',DTCH,QLIN,ST,PROD,LOSS,ZC(159,L))


                PROD = &
                     M_Sestrp * M_NO3 * k_no3_soaC5 * SOA_STOC531 &
                              * (TMASS(trsp_idx(290))/TMASS(trsp_idx(164)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS53',DTCH,QLIN,ST,PROD,LOSS,ZC(164,L))

                !// sesquiterpenes, class5
                PROD = POLLX(290)
                LOSS = M_O3 * k_o3_soaC5 &
                     + M_OH * k_oh_soaC5 &
                     + M_NO3 * k_no3_soaC5
                LOSS_C5 = LOSS
                TMP_ZC = M_Sestrp
                call QSSA(1,'Sestrp',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(290,L) = max(TMP_ZC, 1.e-21_r8) !Do not allow too small values
                if (LOSS .gt. ZTMP) M_Sestrp = ZC(290,L)

             end do !// do SOA_N = 1, 15 (class 1 loop)

             !// Reset time step back to HOx-loop
             DTCH = DTCH2
             ZTMP = 1._r8 / DTCH

             !// Class 1 ----------------------------------------------------
             !// Time step: original (DTCH2)
             PROD = &
                  M_Apine * M_O3 * k_o3_soaC1 * SOA_STOC111 &
                          * (TMASS(trsp_idx(280))/TMASS(trsp_idx(150))) &
                  + M_Apine * M_OH * k_oh_soaC1 * SOA_STOC121 &
                          * (TMASS(trsp_idx(280))/TMASS(trsp_idx(150))) &
                  + M_Bpine * M_O3 * k_o3_soaC1 * SOA_STOC111 &
                          * (TMASS(trsp_idx(281))/TMASS(trsp_idx(150))) &
                  + M_Bpine * M_OH * k_oh_soaC1 * SOA_STOC121 &
                          * (TMASS(trsp_idx(281))/TMASS(trsp_idx(150))) &
                  + M_Sabine * M_O3 * k_o3_soaC1 * SOA_STOC111 &
                          * (TMASS(trsp_idx(284))/TMASS(trsp_idx(150))) &
                  + M_Sabine * M_OH * k_oh_soaC1 * SOA_STOC121 &
                          * (TMASS(trsp_idx(284))/TMASS(trsp_idx(150))) &
                  + M_D3carene * M_O3 * k_o3_soaC1 * SOA_STOC111 &
                          * (TMASS(trsp_idx(285))/TMASS(trsp_idx(150))) &
                  + M_D3carene * M_OH * k_oh_soaC1 * SOA_STOC121 &
                          * (TMASS(trsp_idx(285))/TMASS(trsp_idx(150))) &
                  + M_Trp_Ket * M_O3 * k_o3_soaC1 * SOA_STOC111 &
                          * (TMASS(trsp_idx(291))/TMASS(trsp_idx(150))) &
                  + M_Trp_Ket * M_OH * k_oh_soaC1 * SOA_STOC121 &
                          * (TMASS(trsp_idx(291))/TMASS(trsp_idx(150)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS11',DTCH,QLIN,ST,PROD,LOSS,ZC(150,L))

             PROD = &
                  M_Apine * M_O3 * k_o3_soaC1 * SOA_STOC112 &
                          * (TMASS(trsp_idx(280))/TMASS(trsp_idx(155))) &
                  + M_Apine * M_OH * k_oh_soaC1 * SOA_STOC122 &
                          * (TMASS(trsp_idx(280))/TMASS(trsp_idx(155))) &
                  + M_Bpine * M_O3 * k_o3_soaC1 * SOA_STOC112 &
                          * (TMASS(trsp_idx(281))/TMASS(trsp_idx(155))) &
                  + M_Bpine * M_OH * k_oh_soaC1 * SOA_STOC122 &
                          * (TMASS(trsp_idx(281))/TMASS(trsp_idx(155))) &
                  + M_Sabine * M_O3 * k_o3_soaC1 * SOA_STOC112 &
                          * (TMASS(trsp_idx(284))/TMASS(trsp_idx(155))) &
                  + M_Sabine * M_OH * k_oh_soaC1 * SOA_STOC122 &
                          * (TMASS(trsp_idx(284))/TMASS(trsp_idx(155))) &
                  + M_D3carene * M_O3 * k_o3_soaC1 * SOA_STOC112 &
                          * (TMASS(trsp_idx(285))/TMASS(trsp_idx(155))) &
                  + M_D3carene * M_OH * k_oh_soaC1 * SOA_STOC122 &
                          * (TMASS(trsp_idx(285))/TMASS(trsp_idx(155))) &
                  + M_Trp_Ket * M_O3 * k_o3_soaC1 * SOA_STOC112 &
                          * (TMASS(trsp_idx(291))/TMASS(trsp_idx(155))) &
                  + M_Trp_Ket * M_OH * k_oh_soaC1 * SOA_STOC122 &
                          * (TMASS(trsp_idx(291))/TMASS(trsp_idx(155)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS12',DTCH,QLIN,ST,PROD,LOSS,ZC(155,L))


             !// NO3 oxidation products:
             PROD = &
                  M_Apine * M_NO3 * k_no3_soaC1 * 0.16_r8           &!Yield is ~0.16
                          * (TMASS(trsp_idx(280))/TMASS(trsp_idx(160))) &!(Spittler etal 2006, Atm.Env.)
                  + M_Bpine * M_NO3 * k_no3_soaC1 * SOA_STOC131 &
                          * (TMASS(trsp_idx(281))/TMASS(trsp_idx(160))) &
                  + M_Sabine * M_NO3 * k_no3_soaC1 * SOA_STOC131 &
                          * (TMASS(trsp_idx(284))/TMASS(trsp_idx(160) )) &
                  + M_D3carene * M_NO3 * k_no3_soaC1 * SOA_STOC131 &
                          * (TMASS(trsp_idx(285))/TMASS(trsp_idx(160))) &
                  + M_Trp_Ket * M_NO3 * k_no3_soaC1 * SOA_STOC131 &
                          * (TMASS(trsp_idx(291))/TMASS(trsp_idx(160)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS13',DTCH,QLIN,ST,PROD,LOSS,ZC(160,L))

             !// Loss class1
             LOSS_C1 = &
                  M_O3 * k_o3_soaC1 &
                  + M_OH * k_oh_soaC1 &
                  + M_NO3 * k_no3_soaC1

             !//..alpha pinene, class1 --------------------------------------
             PROD = POLLX(280)
             LOSS = LOSS_C1
             TMP_ZC = M_Apine
             call QSSA(1,'Apine',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
             ZC(280,L) = max(TMP_ZC, 1.e-21_r8)
             if (LOSS .gt. ZTMP) M_Apine = ZC(280,L)

             !//..beta pinene, class1 ---------------------------------------
             PROD = POLLX(281)
             LOSS = LOSS_C1
             TMP_ZC = M_Bpine
             call QSSA(1,'Bpine',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
             ZC(281,L) = max(TMP_ZC, 1.e-21_r8)
             if (LOSS .gt. ZTMP) M_Bpine = ZC(281,L)

             !//..Sabinene, class1 ------------------------------------------
             PROD = POLLX(284)
             LOSS = LOSS_C1
             TMP_ZC = M_Sabine
             call QSSA(1,'Sabine',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
             ZC(284,L) = max(TMP_ZC, 1.e-21_r8)
             if (LOSS .gt. ZTMP) M_Sabine = ZC(284,L)

             !//..d-3-carene, class1 ----------------------------------------
             PROD = POLLX(285)
             LOSS = LOSS_C1
             TMP_ZC = M_D3carene
             call QSSA(1,'D3Carene',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
             ZC(285,L) = max(TMP_ZC, 1.e-21_r8)
             if (LOSS .gt. ZTMP) M_D3carene = ZC(285,L)

             !//..terpinoid ketones, class1 ---------------------------------
             PROD = POLLX(291)
             LOSS = LOSS_C1
             TMP_ZC = M_Trp_Ket
             call QSSA(1,'Trp_Ket',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
             ZC(291,L) = max(TMP_ZC, 1.e-21_r8)
             if (LOSS .gt. ZTMP) M_Trp_Ket = ZC(291,L)


             !// Class 2 loop -----------------------------------------------
             !// Time step: original (DTCH2) shortened by 5
             DTCH = DTCH2 / 5._r8
             ZTMP = 1._r8 / DTCH
             do SOA_N = 1, 5

                PROD = &
                     M_Limon * M_O3 * k_o3_soaC2 * SOA_STOC211 &
                             * (TMASS(trsp_idx(282))/TMASS(trsp_idx(151))) &
                     + M_Limon * M_OH * k_oh_soaC2 * SOA_STOC221 &
                             * (TMASS(trsp_idx(282))/TMASS(trsp_idx(151)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS21',DTCH,QLIN,ST,PROD,LOSS,ZC(151,L))

                PROD = &
                     M_Limon * M_O3 * k_o3_soaC2 * SOA_STOC212 &
                             * (TMASS(trsp_idx(282))/TMASS(trsp_idx(156))) &
                     + M_Limon * M_OH * k_oh_soaC2 * SOA_STOC222 &
                             * (TMASS(trsp_idx(282))/TMASS(trsp_idx(156)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS22',DTCH,QLIN,ST,PROD,LOSS,ZC(156,L))

                PROD = &
                     M_Limon * M_NO3 * k_no3_soaC2 * 0.44_r8      &!Spittler et al.,
                             * (TMASS(trsp_idx(282))/TMASS(trsp_idx(161))) !Atm.Env.2006
                LOSS = 0._r8
                call QSSA(1,'SOAGAS23',DTCH,QLIN,ST,PROD,LOSS,ZC(161,L))

                !// Loss class2
                LOSS_C2 = &
                     M_O3 * k_o3_soaC2 &
                     + M_OH * k_oh_soaC2 &
                     + M_NO3 * k_no3_soaC2

                !//..Limonene, class2 ---------------------------------------
                PROD = POLLX(282)
                LOSS = LOSS_C2
                TMP_ZC = M_Limon
                call QSSA(1,'Limon',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(282,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_Limon = ZC(282,L)

             end do !// do SOA_N = 1, 5 (class 2 loop)



             !// Class 3 loop -----------------------------------------------
             !// Time step: original (DTCH2) shortened by 20
             DTCH = DTCH2 / 20._r8
             ZTMP = 1._r8 / DTCH
             do SOA_N = 1, 20

                PROD = &
                     M_Trpolene * M_O3 * k_o3_soaC3 * SOA_STOC311 &
                            * (TMASS(trsp_idx(287))/TMASS(trsp_idx(152))) &
                     + M_Trpolene * M_OH * k_oh_soaC3 * SOA_STOC321 &
                            * (TMASS(trsp_idx(287))/TMASS(trsp_idx(152))) &
                     + M_Trpinene * M_O3 * k_o3_soaC3 * SOA_STOC311 &
                            * (TMASS(trsp_idx(288))/TMASS(trsp_idx(152))) &
                     + M_Trpinene * M_OH * k_oh_soaC3 * SOA_STOC321 &
                            * (TMASS(trsp_idx(288))/TMASS(trsp_idx(152)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS31',DTCH,QLIN,ST,PROD,LOSS,ZC(152,L))
 
                PROD = &
                     M_Trpolene * M_O3 * k_o3_soaC3 * SOA_STOC312 &
                            * (TMASS(trsp_idx(287))/TMASS(trsp_idx(157))) &
                     + M_Trpolene * M_OH * k_oh_soaC3 * SOA_STOC322 &
                            * (TMASS(trsp_idx(287))/TMASS(trsp_idx(157))) &
                     + M_Trpinene * M_O3 * k_o3_soaC3 * SOA_STOC312 &
                            * (TMASS(trsp_idx(288))/TMASS(trsp_idx(157))) &
                     + M_Trpinene * M_OH * k_oh_soaC3 * SOA_STOC322 &
                            * (TMASS(trsp_idx(288))/TMASS(trsp_idx(157)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS32',DTCH,QLIN,ST,PROD,LOSS,ZC(157,L))

                PROD = &
                     M_Trpolene * M_NO3 * k_no3_soaC3 * SOA_STOC331 &
                            * (TMASS(trsp_idx(287))/TMASS(trsp_idx(162))) &
                     + M_Trpinene * M_NO3 * k_no3_soaC3 * SOA_STOC331 &
                            * (TMASS(trsp_idx(288))/TMASS(trsp_idx(162)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS33',DTCH,QLIN,ST,PROD,LOSS,ZC(162,L))

                !// Loss class3
                LOSS_C3 = &
                     M_O3 * k_o3_soaC3 &
                     + M_OH * k_oh_soaC3 &
                     + M_NO3 * k_no3_soaC3

                !//..terpinolene, class3 ------------------------------------
                PROD = POLLX(287)
                LOSS = LOSS_C3
                TMP_ZC = M_Trpolene
                call QSSA(1,'Trpolene',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(287,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_Trpolene = ZC(287,L)

                !//..alpha and gamma terpinene, class3 ----------------------
                PROD = POLLX(288)
                LOSS = LOSS_C3
                TMP_ZC = M_Trpinene
                call QSSA(1,'Trpinene',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(288,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_Trpinene = ZC(288,L)

             end do !// do SOA_N = 1, 20 (class 3 loop)


             !// Class 4 loop:
             !// Time step: original (DTCH2) shortened by 5
             DTCH = DTCH2 / 5._r8
             ZTMP = 1._r8 / DTCH
             do SOA_N = 1, 5

                PROD = &
                     M_Myrcene * M_O3 * k_o3_soaC4 * SOA_STOC411 &
                             * (TMASS(trsp_idx(283))/TMASS(trsp_idx(153))) &
                     + M_Myrcene * M_OH * k_oh_soaC4 * SOA_STOC421 &
                             * (TMASS(trsp_idx(283))/TMASS(trsp_idx(153))) &
                     + M_Ocimene * M_O3 * k_o3_soaC4 * SOA_STOC411 &
                             * (TMASS(trsp_idx(286))/TMASS(trsp_idx(153))) &
                     + M_Ocimene * M_OH * k_oh_soaC4 * SOA_STOC421 &
                             * (TMASS(trsp_idx(286))/TMASS(trsp_idx(153))) &
                     + M_TrpAlc * M_O3 * k_o3_soaC4 * SOA_STOC411 &
                             * (TMASS(trsp_idx(289))/TMASS(trsp_idx(153))) &
                     + M_TrpAlc * M_OH * k_oh_soaC4 * SOA_STOC421 &
                             * (TMASS(trsp_idx(289))/TMASS(trsp_idx(153)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS41',DTCH,QLIN,ST,PROD,LOSS,ZC(153,L))

                PROD = &
                     M_Myrcene * M_O3 * k_o3_soaC4 * SOA_STOC412 &
                             * (TMASS(trsp_idx(283))/TMASS(trsp_idx(158))) &
                     + M_Myrcene * M_OH * k_oh_soaC4 * SOA_STOC422 &
                             * (TMASS(trsp_idx(283))/TMASS(trsp_idx(158))) &
                     + M_Ocimene * M_O3 * k_o3_soaC4 * SOA_STOC412 &
                             * (TMASS(trsp_idx(286))/TMASS(trsp_idx(158))) &
                     + M_Ocimene * M_OH * k_oh_soaC4 * SOA_STOC422 &
                             * (TMASS(trsp_idx(286))/TMASS(trsp_idx(158))) &
                     + M_TrpAlc * M_O3 * k_o3_soaC4 * SOA_STOC412 &
                             * (TMASS(trsp_idx(289))/TMASS(trsp_idx(158))) &
                     + M_TrpAlc * M_OH * k_oh_soaC4 * SOA_STOC422 &
                             * (TMASS(trsp_idx(289))/TMASS(trsp_idx(158)))
                LOSS=0._r8
                call QSSA(1,'SOAGAS42',DTCH,QLIN,ST,PROD,LOSS,ZC(158,L))

                PROD = &
                     M_Myrcene * M_NO3 * k_no3_soaC4 * SOA_STOC431 &
                             * (TMASS(trsp_idx(283))/TMASS(trsp_idx(163))) &
                     + M_Ocimene * M_NO3 * k_no3_soaC4 * SOA_STOC431 &
                             * (TMASS(trsp_idx(286))/TMASS(trsp_idx(163))) &
                     + M_TrpAlc * M_NO3 * k_no3_soaC4 * SOA_STOC431 &
                             * (TMASS(trsp_idx(289))/TMASS(trsp_idx(163)))
                LOSS = 0._r8
                call QSSA(1,'SOAGAS43',DTCH,QLIN,ST,PROD,LOSS,ZC(163,L))

                !// Loss class4
                LOSS_C4 = &
                     M_O3 * k_o3_soaC4 &
                     + M_OH * k_oh_soaC4 &
                     + M_NO3 * k_no3_soaC4

                !//..Myrcene, class4 ----------------------------------------
                PROD = POLLX(283)
                LOSS = LOSS_C4
                TMP_ZC = M_Myrcene
                CALL QSSA(1,'Myrcene',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(283,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_Myrcene = ZC(283,L)

                !//..ocimene, class4 ----------------------------------------
                PROD = POLLX(286)
                LOSS = LOSS_C4
                TMP_ZC = M_Ocimene
                call QSSA(1,'Ocimene',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(286,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_Ocimene = ZC(286,L)

                !//..terpinoid alcohols, class4 -----------------------------
                PROD = POLLX(289)
                LOSS = LOSS_C4
                TMP_ZC = M_TrpAlc
                call QSSA(1,'TrpAlc',DTCH,QLIN,ST,PROD,LOSS,TMP_ZC)
                ZC(289,L) = max(TMP_ZC, 1.e-21_r8)
                if (LOSS .gt. ZTMP) M_TrpAlc = ZC(289,L)

             end do !// do SOA_N = 1, 5 (class 4 loop)


             !// Reset time step to original (DTCH2)
             DTCH = DTCH2
             ZTMP = 1._r8 / DTCH

             !//C6HXR_SOA products ------------------------------------------
             PROD = &
                  ( M_OH * k_oh_soaC7 &
                    + M_O3 * k_o3_soaC7 &
                  ) * M_C6HXR_SOA * 0.038_r8 &
                                  * (TMASS(trsp_idx(12))/TMASS(trsp_idx(184)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS71',DTCH,QLIN,ST,PROD,LOSS,ZC(184,L))

             PROD = &
                  ( M_OH * k_oh_soaC7 &
                    + M_O3 * k_o3_soaC7 &
                  ) * M_C6HXR_SOA * 0.167_r8 &
                                  * (TMASS(trsp_idx(12))/TMASS(trsp_idx(185)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS72',DTCH,QLIN,ST,PROD,LOSS,ZC(185,L))

             !//..C6HXR_SOA--------------------------------------------------
             !// aromats represented by m-xylene are xylene, and
             !// trimethyl benzene
             PROD = POLLX(12)
             LOSS = &
                  k_oh_soaC7 * M_OH &
                  + k_o3_soaC7 * M_O3 &
                  + k_no3_soaC7 * M_NO3
             call QSSA(41,'C6HXR_SOA',DTCH,QLIN,ST,PROD,LOSS,ZC(12,L))
             !// Originally, C6HXR was not re-set, but with the new definition
             !// its lifetime may become short enough so we should update
             !// the concentration accordingly:
             if (LOSS .gt. ZTMP) M_C6HXR_SOA = ZC(12,L)


             !//ISOPRENE PRODUCTS--------------------------------------------
             PROD = &
                  !M_ISOPREN_TOT * k_oh_isoprene * M_OH &
                  M_ISOPREN * k_oh_isoprene * M_OH &
                            * f_oh_isoprene_sg61 & ! constant is 0.232
                            * (TMASS(trsp_idx(20))/TMASS(trsp_idx(180)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS61',DTCH,QLIN,ST,PROD,LOSS,ZC(180,L))

             PROD = &
                  !M_ISOPREN_TOT * k_oh_isoprene * M_OH &
                  M_ISOPREN * k_oh_isoprene * M_OH &
                            * f_oh_isoprene_sg62 & ! constant is 0.0288
                            * (TMASS(trsp_idx(20))/TMASS(trsp_idx(181)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS62',DTCH,QLIN,ST,PROD,LOSS,ZC(181,L))

             !//..isopren when SOA is calculated-----------------------------
             !// When LSOA=.false., isoprene is treated below OH-loop
             PROD = POLLX(20)
             LOSS = &
                  k_oh_isoprene * M_OH &
                  + 1.e-5_r8
             call QSSA(43,'ISOPREN',DTCH,QLIN,ST,PROD,LOSS,ZC(20,L))

             if (ZC(20,L) .lt. 0._r8) then
                print*, 'OSLO_CHEM: negative ISOPREN in chemistry', &
                     ICOL,JCOL,L,ZC(20,L)
                print*, DTCH,QLIN,ST
                !print*, PROD,LOSS,M_ISOPREN_GAS,M_ISOPREN_TOT
                print*, PROD,LOSS,M_ISOPREN
             end if
             !// Update M_ISOPREN if loss is large
             if (LOSS .gt. ZTMP) M_ISOPREN = ZC(20,L)


             PROD = &
                  (M_OH * k_oh_soaC8 + M_O3 * k_o3_soaC8) * M_Tolmatic * 0.071_r8 &
                        * (TMASS(trsp_idx(192))/TMASS(trsp_idx(186))) &
                  + k_oh_benzene * M_OH * M_Benzene * 0.072_r8 &
                        * (TMASS(trsp_idx(193))/TMASS(trsp_idx(186)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS81',DTCH,QLIN,ST,PROD,LOSS,ZC(186,L))

             PROD = &
                  (M_OH * k_oh_soaC8 + M_O3 * k_o3_soaC8) * M_Tolmatic * 0.138_r8 &
                        * (TMASS(trsp_idx(192))/TMASS(trsp_idx(187))) &
                  + k_oh_benzene * M_OH * M_Benzene * 0.888_r8 &
                        * (TMASS(trsp_idx(193))/TMASS(trsp_idx(187)))
             LOSS = 0._r8
             call QSSA(1,'SOAGAS82',DTCH,QLIN,ST,PROD,LOSS,ZC(187,L))

             !//..tolmatic---------------------------------------------------
             PROD = POLLX(192)
             LOSS = &
                  M_NO3 * k_no3_soaC8 &
                  + M_O3 * k_o3_soaC8 &
                  + M_OH * k_oh_soaC8
             call QSSA(1,'tolmatic',DTCH,QLIN,ST,PROD,LOSS,ZC(192,L))
             if (LOSS .gt. ZTMP) M_TOLMATIC = ZC(192,L)


             !//..benzene----------------------------------------------------
             PROD = POLLX(193)
             LOSS = k_oh_benzene * M_OH
             call QSSA(1,'benzene',DTCH,QLIN,ST,PROD,LOSS,ZC(193,L))
             if (LOSS .gt. ZTMP) M_BENZENE = ZC(193,L)

          else
             !// Not doing SOA. Could treat isoprene and C6HXR here, but
             !// that is done further below.

          end if !// if (LSOA) then
          !// ---------------------------------------------------------------
          !// SOA Secondary organic aerosols: End of SOA chemistry with
          !// shorter time step.
          !// ---------------------------------------------------------------

        end do !// do IIJ = 1, OH_HO2_ITER   (500      CONTINUE)
        !// -----------------------------------------------------------------
        !// END doHOx-LOOP (Time step DTCH2)
        !// -----------------------------------------------------------------



        !// Done with OH-loop: Reset time step and QSSA parameters
        DTCH = DTCH1
        ST   = ST1
        QLIN = QLIN1
        QTEST= QTEST1
        ZTMP = 1._r8 / DTCH


        !// update O(3)P ----------------------------------------------------
        ZC(38,L) = M_O3P


        !//..O3-NO-----------------------------------------------------------
        PROD_NO = &
             POLLX(43) &
             + (k_op_no2 * M_O3P + k_no2_no3_b * M_NO3) * M_NO2

        LOSS_O3 = &
             k_o3_no2 * M_NO2 &
             + k_o3_ho2 * M_HO2 &
             + k_o3_c2h4 * M_C2H4 &
             + ( k_od_h2o * M_H2O * M_O1D &
                 + k_op_no2 * M_NO2 * M_O3P &
               ) / M_O3 &
             + k_o3_c3h6 * M_C3H6 &
             + k_o3_oh * M_OH &
             + VDEP_L(1)
        !// Sulphur reactions
        if (LSULPHUR) LOSS_O3 = LOSS_O3 &
             + CAQ0172 * M_SO2

        if (L .eq. 1) DDDIAG(1) = DDDIAG(1) + VDEP_L(1) * M_O3 * DTCH

        !// Loss NO
        LOSS_2 = &
             k_no_ch3o2 * M_CH3O2 &
             + k_no_ho2 * M_HO2    &!NO + HO2 -> OH + NO2
             + k_no_ho2_b * M_HO2  &!NO + HO2 -> HNO3
             + k_no_no3 * M_NO3     &!NO + NO3 -> 2*NO2
             + k_no_ch3cob * M_CH3COB &
             + k_no_ch3cod * M_CH3COD &
             + k_no_ch3xx * M_CH3XX &
             + k_no_ch3x * M_CH3X &
             + k_no_c2h5o2 * M_C2H5O2 &
             + k_no_c3h7o2 * M_C3H7O2 &
             + k_no_c4h9o2 * M_C4H9O2 &
             + k_no_c6h13o2 * M_C6H13O2 &
             + k_no_ar1 * M_AR1 &
             + k_no_ar3 * M_AR3 &
             + k_no_isor1 * M_ISOR1 &
             + k_no_isor2 * M_ISOR2 * 1.5_r8 &
             + VDEP_L(43) &
             + k_op_no_m * M_O3P !// O3P + NO + M -> NO2

        !// Prod/Loss for O3-NO
        PROD = LOSS_2 * M_NO
        LOSS = &
             LOSS_O3 * M_O3 &
             + PROD_NO

        !// Sum up prod and loss (both are per second here)
        OxCHEMPROD(L) = OxCHEMPROD(L) + PROD*DTCH
        OxCHEMLOSS(L) = OxCHEMLOSS(L) + LOSS*DTCH

        M_O3NO = M_O3NO + (PROD - LOSS) * DTCH


        !//..H2O2------------------------------------------------------------
        PROD = &
             k_ho2_ho2_tot * M_HO2 * M_HO2 &
             + k_oh_oh_m * M_OH * M_OH !OH + OH + M -> H2O2 + M

        LOSS = &
             DH2O2 &
             + k_oh_h2o2 * M_OH &
             + VDEP_L(15)
        !// Sulphur reactions
        if (LSULPHUR) LOSS = LOSS &
             + CAQ1572*M_SO2

        if (L .eq. 1) DDDIAG(15) = DDDIAG(15) + VDEP_L(15) * M_H2O2 * DTCH

        call QSSA(23,'H2O2',DTCH,QLIN,ST,PROD,LOSS,ZC(15,L))


        !//..CO--------------------------------------------------------------
        PROD = &
             POLLX(6) &
             + k_cho_o2 * M_O2 * M_CHO &
             + DBCH2O * M_CH2O &
             + 0.44_r8 * k_o3_c2h4 * M_O3 * M_C2H4 &
             + 0.4_r8 * k_o3_c3h6 * M_O3 * M_C3H6 &
             + DHCOHCO * M_HCOHCO &
             + DRCOHCO * M_RCOHCO &
             + (k_oh_hcohco_m_a + 2._r8*k_oh_hcohco_m_c) * M_OH * M_HCOHCO &
             + k_oh_rcohco * M_OH * M_RCOHCO &
             + DACETON_B * M_ACETON !ACETON + hv -> 2CH3 + CO

        LOSS = &
             (k_oh_co_a + k_oh_co_b) * M_OH &
             + VDEP_L(6)

        if (L .eq. 1) DDDIAG(6) = DDDIAG(6) + VDEP_L(6) * M_CO * DTCH

        call QSSA(24,'CO',DTCH,QLIN,ST,PROD,LOSS,ZC(6,L))



        !//..CH3O2H----------------------------------------------------------
        !// Tracer 16 is O.5*(OH + CH3O2 + CH3) + CH3O2H in the troposphere
        !// but not in stratosphere. Should figure out reason for this.
        !// for OH
        PROD_2 = &
             2._r8 * (DH2O2 * M_H2O2 &
                     + k_od_h2o * M_O1D * M_H2O) &
             + k_o3_ho2 * M_O3 * M_HO2 &
             + k_no_ho2 * M_HO2 * M_NO &
             + 0.15_r8 * k_o3_c3h6 * M_O3 * M_C3H6
        !// for CH3O2
        PROD_2 = PROD_2 &
             + DCH3CHO * M_CH3CHO &
             + k_no_ch3x * M_NO * M_CH3X &
             + k_no_c2h5o2 * M_C2H5O2 * M_NO &
                           * fa_no_c2h5o2 * (1._r8 - fb_no_c2h5o2) &
             + k_no_c3h7o2 * M_C3H7O2 * M_NO &
                           * fa_no_c3h7o2 * (1._r8 - fb_no_c3h7o2)
        !// for CH3
        PROD_2 = PROD_2 &
             + 0.31_r8 * k_o3_c3h6 * M_O3 * M_C3H6 &
             + k_oh_ch3o2h_a * M_OH * M_CH3O2H &
             + (DACETON_A + 2._r8 * DACETON_B) * M_ACETON

        !// for OH
        LOSS_2 = &
             k_no2_oh_m * M_NO2 &
             + k_oh_ch2o * M_CH2O &
             + (k_oh_co_a + k_oh_co_b) * M_CO &
             + k_oh_ho2 * M_HO2 &
             + k_oh_h2o2 * M_H2O2 &
             + k_o3_oh * M_O3 &
             + k_oh_hno3 * M_HNO3 &
             + k_oh_ch3cox * M_CH3COX &
             + k_oh_c3h6_m * M_C3H6 &
             + k_oh_ch3cho * M_CH3CHO &
             + k_oh_c4h10 * M_C4H10 &
             + k_oh_c6h14 * M_C6H14 &
             + k_oh_c6hxr * M_C6HXR &
             + k_oh_c2h6 * M_C2H6 &
             + k_oh_ar2 * M_AR2 &
             + k_oh_isok * M_ISOK &
             !+ k_oh_isoprene * M_ISOPREN_TOT &
             + k_oh_isoprene * M_ISOPREN &
             + k_oh_ho2no2 * M_HO2NO2 &
             + (k_oh_hcohco_m_a + k_oh_hcohco_m_b + k_oh_hcohco_m_c) * M_HCOHCO &
             + k_oh_rcohco * M_RCOHCO &
             + 2._r8 * k_oh_oh_m * M_OH  !OH + OH + M -> H2O2 + M
        !// Sulphur reactions
        if (LSULPHUR) LOSS_2 = LOSS_2 &
             + k_oh_dms_a * M_DMS &
             + C4071b * M_DMS &
             + k_oh_h2s * M_H2S &
             + CTOT4072 * M_SO2
        !// SOA Secondary organic aerosols
        !// Note that for SOA, C6HXR=0 and the contribution from the
        !// old (non-SOA) reaction k_oh_c6hrx above is also zero.
        !// Instead C6HXR_SOA is used, with a different reaction rate.
        if (LSOA) LOSS_2 = LOSS_2 &
              + k_oh_soaC1 * (M_Apine + M_Bpine + M_Sabine &
                         + M_D3Carene + M_Trp_Ket) &
              + k_oh_soaC2 * M_Limon &
              + k_oh_soaC3 * (M_Trpolene + M_Trpinene) &
              + k_oh_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
              + k_oh_soaC5 * M_Sestrp &
              + k_oh_soaC7 * M_C6HXR_SOA &
              + k_oh_soaC8 * M_tolmatic &
              + k_oh_benzene * M_Benzene
        !// for CH3O2
        LOSS_3 = &
             k_no_ch3o2 * M_NO &
             + 2._r8 * k_ch3o2_ch3o2 * M_CH3O2 &
             + k_ch3o2_ch3xx * M_CH3XX &
             + k_ch3o2_c4h9o2 * M_C4H9O2 &
             + k_ch3o2_c2h5o2 * M_C2H5O2 &
             + k_ch3o2_ch3cob * M_CH3COB &
             + k_ch3o2_c6h13o2 * M_C6H13O2 &
             + k_ch3o2_isor1 * M_ISOR1 &
             + k_ch3o2_isor2 * M_ISOR2 &
             + k_ch3o2_ch3x_b * M_CH3X &
             + QAER(L)

        PROD = &
             k_ho2_radical * M_HO2 * M_ARAD &!HO2 + RADICAL --> CH3O2H
             + 0.5_r8 * k_ho2_ch3o2 * M_HO2 * M_CH3O2 &
             + (k_ho2_ch3x * M_HO2 &
                + k_ch3x_ch3x * M_CH3X) * M_CH3X
        LOSS = &
             (k_oh_ch3o2h_a + k_oh_ch3o2h_b) * M_OH &
             + 0.5_r8 * DCH3O2H

        LOSS = ( 0.5_r8 * (LOSS_2 * M_OH &
                          + LOSS_3 * M_CH3O2) &
                 + LOSS * M_CH3O2H ) / ZC(16,L)

        PROD = PROD + 0.5_r8 * PROD_2

        call QSSA(25,'CH3O2H',DTCH,QLIN,ST,PROD,LOSS,ZC(16,L))


        !//..CH4-------------------------------------------------------------
        PROD = POLLX(46)
        LOSS = &
             k_oh_ch4 * M_OH    &!OH + CH4  -> CH3 + H2O
             + VDEP_L(46)     !Soil uptake / dry deposition
        if (.not. LOLD_H2OTREATMENT) LOSS = LOSS &
             + ( k_od_ch4_a   &!O(1D) + CH4 -> OH + CH3
                 + k_od_ch4_b &!O(1D)+CH4 -> (CH3O/CH2OH + H)  -O2-> CH2O+HO2+H
                 + k_od_ch4_c &!O(1D) + CH4 -> H2 + CH2O
               ) * M_O1D
        if (L .eq. 1) DDDIAG(46) = DDDIAG(46) + VDEP_L(46) * M_CH4 * DTCH

        !// Diagnose loss processes [molecules/cm3 in this time step]
        CHEMLOSS(1,46,L) = CHEMLOSS(1,46,L) + LOSS * M_CH4 * DTCH
        CHEMLOSS(2,46,L) = CHEMLOSS(2,46,L) + VDEP_L(46) * M_CH4 * DTCH
        CHEMLOSS(3,46,L) = CHEMLOSS(3,46,L) + k_oh_ch4 * M_OH * M_CH4 * DTCH
        if (.not. LOLD_H2OTREATMENT) then
           CHEMLOSS(4,46,L) = CHEMLOSS(4,46,L) + k_od_ch4_a*M_O1D*M_CH4*DTCH
           CHEMLOSS(5,46,L) = CHEMLOSS(5,46,L) + k_od_ch4_b*M_O1D*M_CH4*DTCH
           CHEMLOSS(6,46,L) = CHEMLOSS(6,46,L) + k_od_ch4_c*M_O1D*M_CH4*DTCH
        end if

        call QSSA(26,'CH4',DTCH,QLIN,ST,PROD,LOSS,ZC(46,L))



        !// -----------------------------------------------------------------
        !// SOA Secondary organic aerosols: SOA deposition
        !// -----------------------------------------------------------------
        if (LSOA .and. L .eq. 1) then
           PROD = 0._r8
           !// Do loop for SOA components
           do SOA_N = 165, 179
              LOSS = VDEP_L(SOA_N)
              RLIM1 = ZC(SOA_N,L) !// ZC before qssa
              call QSSA(SOA_N,'DEPSOA',DTCH,QLIN,ST,PROD,LOSS,ZC(SOA_N,L))
              !// Diagnose SOA loss due to drydeposition
              DDDIAG(SOA_N) = DDDIAG(SOA_N) + (RLIM1 - ZC(SOA_N,L))
           end do
           !// SOA the isoprene oxidation SOA:
           do SOA_N = 182, 183
              LOSS = VDEP_L(SOA_N)
              RLIM1 = ZC(SOA_N,L) !// ZC before qssa
              call QSSA(SOA_N,'DEPSOA',DTCH,QLIN,ST,PROD,LOSS,ZC(SOA_N,L))
              !// Diagnose SOA loss due to drydeposition
              DDDIAG(SOA_N) = DDDIAG(SOA_N) + (RLIM1 - ZC(SOA_N,L))
           end do

           !// SOA the aromatic oxidation product SOA
           do SOA_N = 188, 191
              LOSS = VDEP_L(SOA_N)
              RLIM1 = ZC(SOA_N,L) !// ZC before qssa
              call QSSA(SOA_N,'DEPSOA',DTCH,QLIN,ST,PROD,LOSS,ZC(SOA_N,L))
              !// Diagnose SOA loss due to drydeposition
              DDDIAG(SOA_N) = DDDIAG(SOA_N) + (RLIM1 - ZC(SOA_N,L))
           end do
        end if !// if (LSOA) then
        !// -----------------------------------------------------------------
        !// SOA Secondary organic aerosols: End SOA deposition
        !// -----------------------------------------------------------------
        if (.not. LOLD_H2OTREATMENT) then
           !//..H2-----------------------------------------------------------
           PROD = &
                DBCH2O * M_CH2O       &! CH2O + hv   -> H2 + CO
                + k_od_ch4_c * M_O1D *M_CH4 ! O(1D) + CH4 -> H2 + CH2O
           !// Not included (strat. reaction): ! H + HO2     -> H2 + O2

           LOSS = &
                k_od_h2 * M_O1D          &! O(1D) + H2  -> OH + H
                + k_oh_h2 * M_OH         ! OH + H2     -> H2O + H
           !// Not included (strat. reaction): ! H2 + Cl     -> HCl + H

           CHEMLOSS(1,113,L) = CHEMLOSS(1,113,L) + LOSS*M_H2*DTCH
           CHEMLOSS(2,113,L) = 0._r8 !// No drydep
           CHEMLOSS(3,113,L) = CHEMLOSS(3,113,L) + k_od_h2*M_O1D*M_H2*DTCH
           CHEMLOSS(4,113,L) = CHEMLOSS(4,113,L) + k_oh_h2*M_OH*M_H2*DTCH

           CHEMPROD(1,113,L) = CHEMPROD(1,113,L) + PROD*DTCH
           CHEMPROD(2,113,L) = CHEMPROD(2,113,L) + DBCH2O*M_CH2O*DTCH
           !CHEMPROD(2,113,L) not in strat
           CHEMPROD(4,113,L) = CHEMPROD(4,113,L) + k_od_ch4_c*M_O1D*M_CH4*DTCH

           call QSSA(27,'H2trop',DTCH,QLIN,ST,PROD,LOSS,ZC(113,L))

        end if !// if (.not. LOLD_H2OTREATMENT) then


        !//..HCHO------------------------------------------------------------
        PROD = &
             k_ch3o_o2 * M_O2 * M_CH3O &
             + M_C2H4 * (k_oh_c2h4_m * M_OH &
                         + k_o3_c2h4 * M_O3) &
             + 0.6_r8 * k_ch3o2_ch3o2 * M_CH3O2 * M_CH3O2 &
             + k_ch3o2_ch3x_b * M_CH3O2 * M_CH3X &
             + k_no_c2h5o2 * M_C2H5O2 * M_NO &
                           * fa_no_c2h5o2 * (1._r8 - fb_no_c2h5o2) &
             + 0.5_r8 * M_C3H6 * k_o3_c3h6 * M_O3 &
             + k_no_ch3xx * M_NO * M_CH3XX &
             + k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
                     * (0.25_r8 + 0.5_r8*(1._r8-fb_no_c2h5o2)) &
             + 0.25_r8 * (k_ch3o2_c4h9o2 * M_CH3O2 * M_C4H9O2 &
                         + k_ch3o2_c3h7o2 * M_CH3O2 * M_C3H7O2) &
             + DHCOHCO * M_HCOHCO &
             + k_no_isor1 * M_NO * M_ISOR1 * fa_no_isor1 &
             + k_no_isor2 * M_NO * M_ISOR2 * 0.43_r8 &
             + k_ch3o2_isor1 * M_ISOR1 * M_CH3O2 &
             + k_oh_pan * M_OH * M_PAN &
             + k_oh_ch3o2h_b * M_OH * M_CH3O2H &
             + k_ch3o2_ch3xx * M_CH3O2 * M_CH3XX &
             + POLLX(13)
        if (.not. LOLD_H2OTREATMENT) PROD = PROD &
             + k_od_ch4_c * M_O1D * M_CH4 !O(1D) + CH4 -> H2 + CH2O

        LOSS = &
             DCH2O &
             + k_oh_ch2o * M_OH &
             + VDEP_L(13)
        if (L .eq. 1) DDDIAG(13) = DDDIAG(13) + VDEP_L(13) * M_CH2O * DTCH

        CHEMLOSS(1,13,L) = CHEMLOSS(1,13,L) + LOSS*M_CH2O*DTCH
        CHEMLOSS(2,13,L) = CHEMLOSS(2,13,L) + VDEP_L(13)*DTCH
        CHEMLOSS(3,13,L) = CHEMLOSS(3,13,L) + DACH2O*M_CH2O*DTCH
        CHEMLOSS(4,13,L) = CHEMLOSS(4,13,L) + DBCH2O*M_CH2O*DTCH
        CHEMLOSS(5,13,L) = CHEMLOSS(5,13,L) + k_oh_ch2o*M_OH*M_CH2O*DTCH

        CHEMPROD(1,13,L) = CHEMPROD(1,13,L) + PROD*DTCH
        !CHEMPROD(2,13,L) = CHEMPROD(2,13,L) !Not in troposphere
        !CHEMPROD(3,13,L) = CHEMPROD(3,13,L) !Not in troposphere
        !CHEMPROD(4,13,L) = CHEMPROD(4,13,L) !Not in troposphere
        CHEMPROD(5,13,L) = CHEMPROD(5,13,L) + k_od_ch4_c*M_O1D*M_CH4*DTCH
        CHEMPROD(6,13,L) = CHEMPROD(6,13,L) + k_ch3o_o2*M_O2*M_CH3O*DTCH !not strat
        CHEMPROD(7,13,L) = CHEMPROD(7,13,L) &
             + (0.6_r8 * k_ch3o2_ch3o2 * M_CH3O2 * M_CH3O2 & 
                + k_ch3o2_ch3x_b * M_CH3O2 * M_CH3X &
                + k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
                                 * (0.25_r8 + 0.5_r8 * (1._r8 - fb_no_c2h5o2)) &
                + 0.25_r8 * (k_ch3o2_c4h9o2 * M_CH3O2 * M_C4H9O2 &
                            + k_ch3o2_c3h7o2 * M_CH3O2 * M_C3H7O2) &
             + k_ch3o2_isor1 * M_ISOR1 * M_CH3O2 &
             + k_ch3o2_ch3xx * M_CH3O2 * M_CH3XX ) * DTCH

        call QSSA(27,'HCHO',DTCH,QLIN,ST,PROD,LOSS,ZC(13,L))


        !//..CH3CHO----------------------------------------------------------
        PROD = &
             ( k_no_ch3xx * M_NO &
               + k_ch3o2_ch3xx * M_CH3O2 ) * M_CH3XX &
             + 0.5_r8 * k_o3_c3h6 * M_O3 * M_C3H6 &
             + k_no_c2h5o2 * M_C2H5O2 * M_NO * fa_no_c2h5o2 * fb_no_c2h5o2 &
             + k_no_c3h7o2 * M_C3H7O2 * M_NO &
                           * fa_no_c3h7o2 * (1._r8 - fb_no_c3h7o2) &
             + k_no_c4h9o2 * M_C4H9O2 * M_NO * fa_no_c4h9o2 &
                                      * (1._r8 - fb_no_c4h9o2) &
             + k_no_c6h13o2 * M_C6H13O2 * M_NO * fa_no_c6h13o2 &
             + k_ch3o2_c4h9o2 * M_C4H9O2 * M_CH3O2 &
                              * 0.5_r8 * (1._r8 - fb_no_c4h9o2) &
             + k_ch3o2_c3h7o2 * M_C3H7O2 * M_CH3O2 &
                              * 0.5_r8 * (1._r8 - fb_no_c3h7o2) &
             + k_ch3o2_c2h5o2 * M_CH3O2 * M_C2H5O2 &
                              * (0.25_r8 + 0.5_r8 * fb_no_c2h5o2) &
             + k_ch3o2_c6h13o2 * M_CH3O2 * M_C6H13O2 &
             + DRCOHCO * M_RCOHCO &
             + k_no_isor2 * M_ISOR2 * M_NO * 0.333333_r8 &
             + POLLX(14)

        LOSS = &
             k_oh_ch3cho * M_OH &
             + k_no3_ch3cho * M_NO3 &
             + DCH3CHO &
             + VDEP_L(14)

        if (L .eq. 1) DDDIAG(14) = DDDIAG(14) + VDEP_L(14) * M_CH3CHO * DTCH

        call QSSA(28,'CH3CHO',DTCH,QLIN,ST,PROD,LOSS,ZC(14,L))


        !//..ACETON (CH3COCH3)-----------------------------------------------
        PROD = &
             k_no_c3h7o2 * M_C3H7O2 * M_NO * fa_no_c3h7o2 * fb_no_c3h7o2 &
             + k_ch3o2_c3h7o2 * M_CH3O2 * M_C3H7O2 &
                              * (0.5_r8 * fb_no_c3h7o2 + 0.25_r8) &
             + POLLX(50)

        LOSS = &
             k_oh_aceton * M_OH &
             + DACETON_A &
             + DACETON_B

        call QSSA(29,'ACETON',DTCH,QLIN,ST,PROD,LOSS,ZC(50,L))


        !//..HCOHCO (glyoxal)------------------------------------------------
        PROD = &
             k_no_ar3 * M_AR3 * M_NO &
             + POLLX(35)
        LOSS = &
             DHCOHCO &
             + (k_oh_hcohco_m_a + k_oh_hcohco_m_b + k_oh_hcohco_m_c) * M_OH

        call QSSA(30,'HCOHCO',DTCH,QLIN,ST,PROD,LOSS,ZC(35,L))


        !//..RCOHCO (methyl glyoxal)-----------------------------------------
        PROD = &
             k_no_ar1 * M_NO * M_AR1 &
             + k_no_ar3 * M_NO * M_AR3 &
             + k_no_isor2 * M_NO * M_ISOR2 * 0.67_r8 &
             + ( k_no_ch3cod * M_NO &
                 + k_ch3o2_ch3cod * M_CH3O2) * M_CH3COD &
             + k_ch3o2_isor2 * M_ISOR2 * M_CH3O2

        LOSS = &
             DRCOHCO &
             + k_oh_rcohco * M_OH

        call QSSA(31,'RCOHCO',DTCH,QLIN,ST,PROD,LOSS,ZC(36,L))


        !//..CH3COC2H5 (2-butanone)------------------------------------------
        PROD = &
             k_no_c4h9o2 * M_C4H9O2 * M_NO * fa_no_c4h9o2 * fb_no_c4h9o2 &
             + k_ch3o2_c4h9o2 * M_C4H9O2 * M_CH3O2 &
                              * (0.25_r8 + 0.5_r8 * fb_no_c4h9o2)

        LOSS = &
             DCH3COX &
             + k_oh_ch3cox * M_OH

        call QSSA(32,'CH3CO...',DTCH,QLIN,ST,PROD,LOSS,ZC(19,L))


        !//..CH3COCOCH3 (biacetyl)-------------------------------------------
        PROD = &
             k_no_ch3cob * M_CH3COB * M_NO * fa_no_ch3cob &
             + k_ch3o2_ch3cob * M_CH3O2 * M_CH3COB
        LOSS = &
             DCH3COY &
             + 1.e-12_r8

        call QSSA(33,'CH3COCO.',DTCH,QLIN,ST,PROD,LOSS,ZC(18,L))


        !//..HNO3------------------------------------------------------------
        PROD = &
             k_no2_oh_m * M_OH * ZC(44,L) &
             + k_no3_ch3cho * M_NO3 * M_CH3CHO &
             + k_no3_dms * M_DMS * M_NO3 &
             + POLLX(4) &
             + 2._r8 * k_n2o5_h2o_aer * M_N2O5 &
             + k_no_ho2_b * M_NO * M_HO2     !NO + HO2 -> HNO3
        !// Sulphur reactions
        if (LSULPHUR) PROD = PROD &
             + CAQ1772 * M_HO2NO2 * M_SO2 !HO2NO2 + HSO3 -> HNO3 + SO4

        LOSS = &
             k_oh_hno3 * M_OH &
             + DHNO3 &
             + VDEP_L(4)

        if (L .eq. 1) DDDIAG(4) = DDDIAG(4) + VDEP_L(4) * M_HNO3 * DTCH

        call QSSA(34,'HNO3',DTCH,QLIN,ST,PROD,LOSS,ZC(4,L))


        !//..C2H6------------------------------------------------------------
        PROD = POLLX(8)
        LOSS = k_oh_c2h6 * M_OH 

        call QSSA(35,'C2H6',DTCH,QLIN,ST,PROD,LOSS,ZC(8,L))


        !//..C2H4------------------------------------------------------------
        PROD = POLLX(7)
        LOSS = &
             k_o3_c2h4 * M_O3 &
             + k_oh_c2h4_m * M_OH

        call QSSA(36,'C2H4',DTCH,QLIN,ST,PROD,LOSS,ZC(7,L))


        !//..C3H8------------------------------------------------------------
        PROD = POLLX(48)
        LOSS = k_oh_c3h8 * M_OH 

        call QSSA(37,'C3H8',DTCH,QLIN,ST,PROD,LOSS,ZC(48,L))


        !//..C3H6------------------------------------------------------------
        PROD = POLLX(9)
        LOSS = &
             k_oh_c3h6_m * M_OH &
             + k_o3_c3h6 * M_O3

        call QSSA(38,'C3H6',DTCH,QLIN,ST,PROD,LOSS,ZC(9,L))


        !//..C4H10-----------------------------------------------------------
        PROD = POLLX(10)
        LOSS = k_oh_c4h10 * M_OH

        call QSSA(39,'C4H10',DTCH,QLIN,ST,PROD,LOSS,ZC(10,L))


        !//..C6H14-----------------------------------------------------------
        PROD = POLLX(11)
        LOSS = k_oh_c6h14 * M_OH

        call QSSA(40,'C6H14',DTCH,QLIN,ST,PROD,LOSS,ZC(11,L))


        !//SOA Secondary organic aerosols------------------------------------
        if (.not. LSOA) then
           !// When using SOA, these components are calculated
           !// inside the OH-loop.

           !//..C6HXR--------------------------------------------------------
           !// aromats represented by m-xylene
           PROD = POLLX(12)
           LOSS = k_oh_c6hxr * M_OH
           CALL QSSA (41,'C6HXR',DTCH,QLIN,ST,PROD,LOSS,ZC(12,L))


           !//..isopren------------------------------------------------------
           PROD = POLLX(20)
           LOSS = &
                k_oh_isoprene * M_OH &
                + 1.e-5_r8
           call QSSA(43,'ISOPREN',DTCH,QLIN,ST,PROD,LOSS,ZC(20,L))        

           if (ZC(20,L) .lt. 0._r8) then
              print*, 'OSLO_CHEM: negative ISOPREN',ICOL,JCOL,L,ZC(20,L)
              print*, DTCH,QLIN,ST
              !print*, PROD,LOSS,M_ISOPREN_GAS,M_ISOPREN_TOT
              print*, PROD,LOSS,M_ISOPREN
           end if
           !// Update M_ISOPREN if loss is large
           if (LOSS .gt. ZTMP) M_ISOPREN = ZC(20,L)

        end if !// if (.not. LSOA) then

        !//..NOX-------------------------------------------------------------
        !// NO + NO2 + NO3 + 2*N2O5 + HO2NO2 + PAN
        PROD = &
             POLLX(5) &
             + POLLX(17) &
             + POLLX(41) &
             + 2._r8 * POLLX(42) &
             + POLLX(43) &
             + POLLX(44) &
             + ( k_oh_hno3 * M_OH &
                 + DHNO3) * M_HNO3

        !// Sum up LOSS and divide by NOX (after sums)
        LOSS = &
             + k_no2_oh_m * (ZC(40,L) + M_OH) * 0.5_r8 * ZC(44,L) &
             + ( k_no3_ch3cho * M_CH3CHO &
                 + k_no3_dms * M_DMS &
               ) * M_NO3 &
             + 2._r8 * k_n2o5_h2o_aer * M_N2O5 &
             + ( (1._r8 - fa_no_c2h5o2) * M_C2H5O2 * k_no_c2h5o2 &
                 + (1._r8 - fa_no_ch3o2) * M_CH3O2 * k_no_ch3o2 &
                 + (1._r8 - fa_no_c3h7o2) * M_C3H7O2 * k_no_c3h7o2 &
                 + (1._r8 - fa_no_c4h9o2) * M_C4H9O2 * k_no_c4h9o2 &
                 + (1._r8 - fa_no_c6h13o2) * M_C6H13O2 * k_no_c6h13o2 &
                 + (1._r8 - fa_no_ch3cob) * M_CH3COB * k_no_ch3cob &
                 + (1._r8 - fa_no_ch3cod) * M_CH3COD * k_no_ch3cod &
                 + (1._r8 - fa_no_isor1) * M_ISOR1 * k_no_isor1 &
               ) * ZC(43,L) &
             + VDEP_L(44) * M_NO2 &
             + VDEP_L(5) * M_PAN &
             + VDEP_L(43) * M_NO &
             + k_no_ho2_b * M_NO * M_HO2  !NO + HO2 -> HNO3
        !// Sulphur reactions
        if (LSULPHUR) LOSS = LOSS &
             + CAQ1772 * M_SO2 * M_HO2NO2
        !// SOA Secondary organic aerosols
        if (LSOA) LOSS = LOSS &
             + ( k_no3_soaC1 * (M_Apine + M_Bpine + M_Sabine &
                           + M_D3Carene + M_Trp_Ket) &
                 + k_no3_soaC2 * M_Limon &
                 + k_no3_soaC3 * (M_Trpolene + M_Trpinene) &
                 + k_no3_soaC4 * (M_Myrcene + M_Ocimene + M_TrpAlc) &
                 + k_no3_soaC5 * M_Sestrp &
                 + k_no3_soaC7 * M_C6HXR_SOA &
                 + k_no3_soaC8 * M_tolmatic &
               ) * M_NO3

        !// Divide by NOX
        LOSS = LOSS / M_NOX

        call QSSA(44,'NOX',DTCH,QLIN,ST,PROD,LOSS,M_NOX)


        !//..H2O from aircraft-----------------------------------------------
        !// (Note that H2O below 400hPa will be removed in separate routine).
        PROD = POLLX(148)
        LOSS = 0._r8
        call QSSA(148,'AC_H2O',DTCH,QLIN,ST,PROD,LOSS,ZC(148,L))


        if (LSULPHUR) then
           !//..DMS----------------------------------------------------------
           PROD = POLLX(71)
           LOSS = &
                k_no3_dms * M_NO3 &
                + k_oh_dms_a * M_OH &
                + C4071b * M_OH

           call QSSA(45,'DMS',DTCH,QLIN,ST,PROD,LOSS,ZC(71,L))


           !//..SO2----------------------------------------------------------
           PROD = &
                k_oh_dms_a * M_OH * M_DMS &
                + 0.75_r8 * C4071b * M_OH * M_DMS &
                + k_no3_dms * M_NO3 * M_DMS &
                + k_oh_h2s * M_OH * M_H2S &
                + POLLX(72)
           LOSS = &
                CTOT4072 * M_OH    &! gas phase
                + CAQ1572 * M_H2O2 &! aq. phase
                + CAQ1772 * M_HO2NO2 &
                + CAQ0172 * M_O3   &
                + VDEP_L(72)       &! dry dep.
                + CCATSO2           ! catalytic

           if (L .eq. 1) DDDIAG(72) = DDDIAG(72) + VDEP_L(72) * M_SO2 * DTCH

           call QSSA(46,'SO2',DTCH,QLIN,ST,PROD,LOSS,ZC(72,L))


           if (.not.LM7) then
              !//..SO4-------------------------------------------------------
              PROD = &
                   ( CTOT4072 * M_OH &
                     + CAQ1572 * M_H2O2 &
                     + CAQ1772 * M_HO2NO2 &
                     + CAQ0172 * M_O3 &
                     + CCATSO2 &
                   ) * M_SO2 &
                   + POLLX(73)
              LOSS = VDEP_L(73)

              call QSSA(47,'SO4',DTCH,QLIN,ST,PROD,LOSS,ZC(73,L))

              if (L .eq. 1) DDDIAG(73) = DDDIAG(73) + VDEP_L(73) * M_SO4 * DTCH

           else  !Using M7
              !// M7 splits gaseous and aquous SO4
              !// Gaseous SO4------------------------------------------------
              !// POLLX is done in prodlossm7
              PROD = CTOT4072 * M_OH * M_SO2  
              !CHEMPROD(L,80) = PROD
              LOSS = 0._r8

              call QSSA(47,'SO4Gas',DTCH,QLIN,ST,PROD,LOSS,ZC(80,L))

              !// Aquous SO4-------------------------------------------------
              !// POLLX and DEPSO4 is done in prodlossm7
              PROD = &
                    ( CAQ1572 * M_H2O2 &
                      + CAQ1772 * M_HO2NO2 &
                      + CAQ0172 * M_O3 &
                      + CCATSO2 &
                    ) * M_SO2
              LOSS = 0._r8

              call QSSA(47,'SO4AQ',DTCH,QLIN,ST,PROD,LOSS,ZC(73,L))

           end if !// if (.not.LM7) then


           !//..H2S----------------------------------------------------------
           PROD = POLLX(74)
           LOSS = k_oh_h2s * M_OH
           call QSSA(48,'H2S',DTCH,QLIN,ST,PROD,LOSS,ZC(74,L))


           !//..MSA----------------------------------------------------------
           PROD = 0.25_r8 * C4071b * M_OH * M_DMS
           LOSS = VDEP_L(75)
           if (L .eq. 1) DDDIAG(75) = DDDIAG(75) + VDEP_L(75) * M_MSA * DTCH
           call QSSA(49,'MSA',DTCH,QLIN,ST,PROD,LOSS,ZC(75,L))

        end if !// if (LSULPHUR) then


        if (LNITRATE) then
           !//..NH3gas-------------------------------------------------------
           PROD = POLLX(61)     !NH3 Produced by emissions
           LOSS = VDEP_L(61)     !NH3 lost by drydep
           if (L .eq. 1) DDDIAG(61) = DDDIAG(61) + VDEP_L(61) * NH3gas * DTCH
           call QSSA(61,'NH3',DTCH,QLIN,ST,PROD,LOSS,ZC(61,L))

           !//..NH4fine------------------------------------------------------
           PROD = 0._r8
           LOSS = VDEP_L(62)
           if (L .eq. 1) DDDIAG(62) = DDDIAG(62) + VDEP_L(62) * NH4fine * DTCH
           call QSSA(62,'NH4fine',DTCH,QLIN,ST,PROD,LOSS,ZC(62,L))

           !//..NH4coarse----------------------------------------------------
           PROD = 0._r8
           LOSS = VDEP_L(63)
           if (L .eq. 1) DDDIAG(63) = DDDIAG(63) + VDEP_L(63) * NH4coarse * DTCH
           call QSSA(63,'NH4coarse',DTCH,QLIN,ST,PROD,LOSS,ZC(63,L))

           !//..NO3fine------------------------------------------------------
           PROD = 0._r8
           LOSS = VDEP_L(64)
           if (L .eq. 1) DDDIAG(64) = DDDIAG(64) + VDEP_L(64) * NO3fine * DTCH
           call QSSA(64,'NO3fine',DTCH,QLIN,ST,PROD,LOSS,ZC(64,L))

           !//..NO3coarse--v
           PROD = 0._r8
           LOSS = VDEP_L(65)
           if (L .eq. 1) DDDIAG(65) = DDDIAG(65) + VDEP_L(65) * NO3coarse * DTCH
           call QSSA(65,'NO3coarse',DTCH,QLIN,ST,PROD,LOSS,ZC(65,L))
        end if !// if (LNITRATE) then

        !//..OH update OH----------------------------------------------------
        ZC(40,L) = M_OH


        !// Balancing the nitrogen components
        !// -----------------------------------------------------------------

        RLIM1 = M_NOX * 0.01_r8

        if (LOSS_NO .lt. LOSS_NO2 .or. M_O3NO .lt. 0._r8) then
           !// NO-loss > NO2 loss, but NO > O3

           !// Get NO from NOX and the NOx components
           XJNO  = M_NOX - ZC(44,L) - M_NO3 - 2._r8*M_N2O5 &
                   - ZC(17,L) - M_PAN

           !// Test value for O3
           O3TEST = M_O3NO + XJNO

           if (XJNO .lt. RLIM1 .or. O3TEST .lt. 0._r8) then
              !// Unusual condition: O3TEST < 0 or NO is small (< 1% of NOX)

              !// NO rather small or test-O3 negative; will check sum of NOx
              !// against NOX. Usually the sum of individually calculated
              !// NOx species will be larger, scaling up NO; then
              !// O3=O3NO+NO may be non-negative.
              !// HOWEVER, O3NO+NO may still be negative. This must be
              !// resolved, and is taken care of by the hack further down.

              !// Sum up present NOx
              EKSTRSN = ZC(43,L) + ZC(44,L) + 2._r8 * M_N2O5 + M_NO3 &
                        + M_PAN + ZC(17,L)
              !// Scale the produced NOx-tracer
              RELATN = M_NOX / EKSTRSN

              !// Scale NOx components
              ZC(43,L) = ZC(43,L) * RELATN !// NO
              ZC(44,L) = ZC(44,L) * RELATN !// NO2
              M_N2O5   = M_N2O5 * RELATN   !// N2O5
              ZC(17,L) = ZC(17,L) * RELATN !// HO2NO2
              M_NO3    = M_NO3 * RELATN    !// NO3
              M_NOZ    = M_NOZ * RELATN    !// NOZ
              ZC(5,L)  = ZC(5,L) * RELATN  !// PANX

              !// New O3
              O3TEST = M_O3NO + ZC(43,L)
           else
              !// Usual condition: O3TEST > 0

              !// No need for scaling NOx, so we use NO from the sum
              !// In this case we O3=O3TEST
              ZC(43,L) = XJNO
           end if
           !// Set O3
           ZC(1,L) = O3TEST

        else
           !// NO larger: (NO < O3 or NO-loss > NO2 loss)

           !// Get NO2 from NOX and the NOx components
           XJNO2 = M_NOX - ZC(43,L) - M_NO3 - 2._r8*M_N2O5 &
                - ZC(17,L) - M_PAN

           !// Is NO2 smaller than 1% of NOx?
           if (XJNO2 .lt. RLIM1) then
              !// Sum up present NOx
              EKSTRSN = ZC(44,L) + ZC(43,L) + 2._r8*M_N2O5 + M_NO3 &
                        + M_PAN + ZC(17,L)

              !// Scale the produced NOx-tracer
              RELATN  = M_NOX / EKSTRSN

              !// Scale NOx components
              ZC(44,L) = ZC(44,L) * RELATN !// NO2
              ZC(43,L) = ZC(43,L) * RELATN !// NO
              M_N2O5   = M_N2O5 * RELATN   !// N2O5
              ZC(17,L) = ZC(17,L) * RELATN !// HO2NO2
              M_NO3    = M_NO3 * RELATN    !// NO3
              M_NOZ    = M_NOZ * RELATN    !// NOZ
              ZC(5,L)  = ZC(5,L) * RELATN  !// PANX

           else
              !// NO2 > 1% NOX
              ZC(44,L) = XJNO2
           end if

           !// Set O3
           ZC(1,L) = M_O3NO + ZC(43,L)
        end if !// if (LOSS_NO .lt. LOSS_NO2 .or. O3NO .lt. 0._r8) then

        !// Hack to avoid negative ozone
        if ( ZC(1,L) .lt. 0._r8) then
           !// This scheme does from time to time produce some negative
           !// ozone. The reason for this is the O3-NO treatment, and that
           !// the time step is too long so that dry deposition removes
           !// too much O3. Possibly, this could also occur if NO production
           !// is very large.
           !// When negative O3 happens, it is reset to 1.1[NO]. There is no
           !// explanation of this value, but CTM2 skipped the warning.
           !// For CTM3 we count the number it happens, and also skip
           !// the warnings.
           COUNTnegO3(L) = COUNTnegO3(L) + 1._r8
           !// If you need to look into the specifics of the negative
           !// O3, you can include this print-outs:
           !write(*,'(a,2i4,i3,2es12.4)')'NEGATIVE O3',ICOL,JCOL,l,&
           !                             ZC(1,L),ZC(43,L)
           !// In CTM1 or older, O3 was set to 1.d+9: ZC(1,L) = 1.d+9
           !// but CTM2 always used 1.1[NO]
           ZC(1,L) = 1.1_r8 * ZC(43,L)
           if (LDEBUG) then
              if (zc(1,L).ne.zc(1,L)) then
                 print*,'OSLO_CHEM: O3 is NaN',ICOL,JCOL,L
                 print*,'check prod and loss for o3no, e.g. vdep - STOP'
                 stop
              end if
           end if
        end if

        !// Adjust NO3 and N2O5 according to NOZ
        if (M_NOZ .gt. 0._r8) then
           if ( (M_NOZ - M_N2O5) .gt. 0._r8) then
              M_NO3 = M_NOZ - M_N2O5
           else
              !// set N2O5 to 95% and NOx to be 5% of NOZ
              M_N2O5   = 0.99_r8 * M_NOZ
              M_NO3  = M_NOZ - M_N2O5
           end if
        else if (M_NOZ .lt. 0._r8) Then
           print*, 'OSLO_CHEM: NOZ is negative --- STOP ----'
           print*, 'NOZ,N2O5,NO3',M_NOZ,M_N2O5,M_NO3
           stop
        else
           M_NOZ = 1.d+2
           !// set N2O5 to 99% and NO3 to be 1% of NOZ
           M_N2O5 = 0.99_r8 * M_NOZ
           M_NO3  = M_NOZ - M_N2O5
        end if

        ZC(41,L) = M_NO3
        ZC(42,L) = M_N2O5


        !// ----------------------------------------------------------------
        !// Chemistry calculations finished --------------------------------
        !// ----------------------------------------------------------------

        !// ----------------------------------------------------------------
      end do !// do NST = 1, NCHEM_ITER
      !// ------------------------------------------------------------------
      !// ------------------------------------------------------------------
    end do !// do L = 1, LMTROP
    !// --------------------------------------------------------------------

    !// --------------------------------------------------------------------
  end subroutine OSLO_CHEM
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module pchemc_ij
!//=========================================================================

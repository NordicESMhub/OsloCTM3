!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, April 2017
!//=========================================================================
!// Chemistry reaction rates.
!//=========================================================================
module chem_oslo_rates
  !// ----------------------------------------------------------------------
  !// MODULE: chem_oslo_rates
  !// DESCRIPTION: Reaction rates for both troposphere and stratospherie.
  !//
  !// Contains:
  !//   subroutine TCRATE_CONST
  !//   subroutine set_pr42het
  !//   subroutine TCRATE_TP_IJ_TRP
  !//   subroutine TCRATE_TP_IJ_STR
  !//   subroutine TCRATE_onAER
  !//   subroutine TCRATE_HET_IJ
  !//   subroutine inSAD
  !//
  !// Amund Sovde Haslerud, April 2017
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_parameters, only: MINTEMP, MAXTEMP, TEMPRANGE, AVOGNR
  use cmn_size, only: IPAR, JPAR, LPAR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Constant rates (independent of T and p)
  !// ----------------------------------------------------------------------
  real(r8) :: &
       r_od_ch4_a, r_od_ch4_b, r_od_ch4_c, &
       r_od_h2, &
       r_oh_c6h14, &
       r_oh_c6hxr, &
       r_oh_ar2, &
       r_cho_o2, &
       r_ch3co_o2, &
       r_ch3o2_c2h5o2, &
       r_ch3o2_c4h9o2, &
       r_ch3o2_c6h13o2, &
       r_ch3o2_ch3cob, &
       r_ch3o2_ch3xx, &
       r_ch3o2_c3h7o2, &
       r_ch3o2_ch3cod, &
       r_ch3o2_isor1, &
       r_ch3o2_isor2, &
       r_no_c3h7o2, &
       r_h_ho2_a, r_h_ho2_b, r_h_ho2_c, &
       r_op_hno3, &
       r_od_cfc11_a, r_od_cfc11_b, &
       r_od_hcfc123, r_od_hcfc141, r_od_hcfc142, r_od_hcl, &
       r_h2o_clono2, r_hcl_clono2, &
       r_oh_bro_a



  !// Rates dependent on temperature only
  !// ----------------------------------------------------------------------
  real(r8),dimension(TEMPRANGE) :: &
       r_od_m, & 
       r_od_h2o, &
       r_op_no2, &
       r_o3_no, &
       r_o3_no2, &
       r_o3_oh, &
       r_o3_ho2, &
       r_o3_c2h4, &
       r_o3_c3h6, &
       r_oh_h2o2, &
       r_oh_ho2, &
       r_oh_ho2no2, &
       r_oh_h2, &
       r_oh_ch4, &
       r_oh_ch3o2h_a, r_oh_ch3o2h_b, &
       r_oh_pan, &
       r_oh_c2h6, &
       r_oh_c3h8, &
       r_oh_c4h10, &
       r_oh_isoprene, &
       r_oh_ch2o, &
       r_oh_ch3oh, &
       r_oh_ch3cho, &
       r_oh_rcohco, &
       r_oh_isok, &
       r_oh_aceton, &
       r_oh_ch3cox, &
       r_oh_dms_a, &
       r_oh_h2s, &
       r_oh_ch3ccl3, &
       !r_oh_nh3, &
       r_no_ho2, &
       r_no_ch3o2, &
       r_no_c2h5o2, &
       r_no_c4h9o2, &
       r_no_c6h13o2, &
       r_no_ar1, &
       r_no_ar3, &
       r_no_isor1, &
       r_no_isor2, &
       r_no_ch3cob, &
       r_no_ch3x, &
       r_no_ch3cod, &
       r_no_ch3xx, &
       r_no_no3, &
       r_no3_dms, &
       r_no3_ch3cho, &
       r_no2_no3_b, &
       r_ho2_ch3o2, &
       r_ho2_ch3x, &
       r_ho2_radical, &
       r_ch3o2_ch3o2, &
       r_ch3o2_ch3x_a, &
       r_ch3o2_ch3x_b, &
       r_ch3x_ch3x, &
       r_ch3o_o2, &
       !// SOA Rates for hydrocarbon oxidation
       r_o3_soaC1, r_oh_soaC1, r_no3_soaC1, &
       r_o3_soaC2, r_oh_soaC2, r_no3_soaC2, &
       r_o3_soaC3, r_oh_soaC3, r_no3_soaC3, &
       r_o3_soaC4, r_oh_soaC4, r_no3_soaC4, &
       r_o3_soaC5, r_oh_soaC5, r_no3_soaC5, &
       !// SOA Aromatic oxidation rates
       r_o3_soaC7, r_oh_soaC7, r_no3_soaC7, &
       r_o3_soaC8, r_oh_soaC8, r_no3_soaC8, &
       r_oh_benzene, &
       !// Stratospheric reactions
       r_n_no, &
       r_n_o2, &
       r_op_o3, r_op_clo, r_op_hcl, r_op_clono2, r_op_oh, r_op_ho2, &
       r_op_ch2o, r_op_oclo, r_op_hbr, r_op_bro, &
       r_od_n2, r_od_o2, r_od_n2o_a, r_od_n2o_b, &
       r_o3_h, r_o3_cl, r_o3_bro, &
       r_no_clo, r_no_bro, r_no_oclo, &
       r_oh_oh, r_oh_clo_a, r_oh_clo_b, r_oh_hocl, r_oh_hcl, r_oh_ch3cl, &
       r_oh_chclf2, r_oh_clono2_a, r_oh_clono2_b, &
       r_oh_bro_b, r_oh_hbr, r_oh_br2, r_oh_ch3br, &
       r_oh_hcfc123, r_oh_hcfc141, r_oh_hcfc142, &
       r_ho2_cl_a, r_ho2_cl_b, r_ho2_clo, &
       r_cl_h2, r_cl_h2o2, r_cl_ch4, r_cl_ch3oh, r_cl_ch2o, r_cl_clono2, &
       r_br_o3, r_br_h2o2, r_br_ch2o, r_br_ho2, &
       r_clo_co, r_bro_ho2,  &
       r_bro_clo_a, r_bro_clo_b, r_bro_clo_c, &
       r_bro_bro_a, r_bro_bro_b





  !// Nitrogen fractions for PAN reactions
  !// fa: fraction of NO going to NO2 in gas phase, and not staying as PAN
  !// fb: another fraction
  real(r8), parameter :: &
       fa_no_ch3o2 = 1._r8, &
       fa_no_c2h5o2 = 1._r8, &
       fa_no_c4h9o2 = 1._r8, &
       fa_no_c6h13o2 = 1._r8, &
       fa_no_ch3cob = 1._r8, &
       fa_no_isor1 = 1._r8, &
       fa_no_ch3cod = 1._r8, &
       fa_no_c3h7o2 = 1._r8, &
       fb_no_c2h5o2 = 0._r8, &
       fb_no_c3h7o2 = 0._r8, &
       fb_no_c4h9o2 = 0._r8

  !// Stratospheric branching ratios
  real(r8), parameter :: &
       fb_hv_ohcl = 0._r8

  !// Surface area density for tropospheric aerosols
  real(r8), dimension(LPAR,IPAR,JPAR,12) :: trop_SAD

  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'chem_oslo_rates.f90'
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
  save
  public
  private set_pr42het, getRQAER
  !// ----------------------------------------------------------------------

contains



  !// ----------------------------------------------------------------------
  subroutine TCRATE_CONST2()
    !// --------------------------------------------------------------------
    !// Find constant and temperature dependent reaction rates for
    !// the troposphere and stratosphere (i.e. not pressure dependent).
    !//
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    use cmn_oslo, only: PR42HET
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// local parameters
    real(r8) ::  ZTEM, TEM, Z298, FC
    integer :: I
    !// --------------------------------------------------------------------

    !//-------------------------------------------------------
    !// Reaction rates independent of temperature and pressure
    !//-------------------------------------------------------

    !// OD + CH4 -> OH + CH3
    !// JPL number: A11 (JPL06, 20080616)
    ! JPL06: 1.5E-10 and branching ratio of 75 +/- 15%.
    ! CH3 reacts with O2 instantaneously to form CH3O2
    r_od_ch4_a = 1.125e-10_r8 ! this is 75%

    !// OD + CH4 -> CH3O/CH2OH + H
    !// JPL number: A11 (JPL06, 20080616)
    ! Not included in the Oslo CTM2, but assumed to react instantaneously
    ! as CH3O + O2 -> HO2 + CH2O ! JPL06: Branching ratio of 20 +/- 7%
    r_od_ch4_b = 3.0e-11_r8 ! this is 20%

    !// OD + CH4 -> H2 + CH2O
    !// JPL number: A11 (JPL06, 20080616)
    ! JPL06: Channel to CH2O + H2 has the branching ratio of 5 +/- 5%.
    ! For CH2O formation the reaction rate is then C256B+C256C
    r_od_ch4_c = 7.5e-12_r8 ! this is 5%

    !// OD + H2 -> OH + H
    !// JPL number: A5 (JPL06, 20080617)
    r_od_h2 = 1.1e-10_r8 !JPL06, 20080617

    !// OP + NO3 -> NO2 + O2
    !// R3841 =  1.e-11_r8

    !// OH + C6H14 --> C6H13O2 (+ H2O)
    !// No reference so far
    r_oh_c6h14 = 5.61e-12_r8

    !// OH + C6HXR --> AR1
    !// No reference so far
    r_oh_c6hxr = 1.04e-11_r8

    !// CH3CO + O2 --> CH3X
    !// IUPAC number: R_Oxygen11 (IUPAC06, 20080612)
    r_ch3co_o2 = 5.10e-12_r8

    !// CHO + O2 --> HO2 + CO
    !// JPL number: D49 (JPL10, 20140612) (same as IUPAC, R_Oxygen10)
    r_cho_o2 = 5.20e-12_r8

    !// OH + AR2 --> AR3
    !// No reference so far
    r_oh_ar2 = 2.00e-11_r8

    !// OH + ISOK --> ISOR2 updated to T-dependent (20140612)
    !// OH + RCOHCO --> CH3CO+CO updated to T-dependent (20140612)
    !// NO + C2H5O2 --> fa4323*(NO2 + (1-fb4323)*(CH2O + CH3) + 
    !//                     fb4323*(HO_2 + CH3CHO))
    !//                updated to T-dependent (20140612)
    !// NO + CH3X --> NO2 + CH3 (+ O2) now T-dependent (JPL06)

    !// NO + C3H7O2 --> fa4349*(NO2 + (1-fb4349)*(CH3CHO + CH3) + 
    !//                 fb4349*(HO_2 + ACETON))
    !// IUPAC number: ROO_4 or ROO_5
    r_no_c3h7o2 = 4.90e-12_r8

    !// CH3O2 + C2H5O2 --> 0.5( CH3O + (1-fb4323)*(CH3 + HCHO)
    !//                         + fb4323(HO2 + CH3CHO) )
    !//                    + 0.25 CH2O + 0.25 CH3CHO
    !// IUPAC number: ROO_41 (they suggest 2.e-13_r8)
    r_ch3o2_c2h5o2 = 1.00e-13_r8

    !// CH3O2 + C4H9O2 --> 0.5( CH3O + (1-fb4324)*(C2H5O2 + CH3CHO)
    !//                         + fb(HO2 + CH3COX) )
    !//                    + 0.25 CH2O + 0.25 CH3COX
    !// IUPAC number: use ROO_41
    r_ch3o2_c4h9o2 = 1.00e-13_r8

    !// CH3O2 + C6H13O2 --> C4H9O2 + CH3O + CH3CHO
    !// IUPAC number: use ROO_41
    r_ch3o2_c6h13o2 = 1.00e-13_r8

    !// CH3O2 + CH3COB --> HO2 + CH3O + CH3COY
    !//   CH3COB = CH3COCH(O2)CH3
    !//   CH3COY = CH3COCOCH3
    !// IUPAC number: use ROO_41
    r_ch3o2_ch3cob = 1.00e-13_r8

    !// CH3O2 + CH3XX --> CH3O + CH3CHO + CH2O + HO2
    !//   CH3XX = CH3CH(O2)CH2OH
    !// IUPAC number: use ROO_41
    r_ch3o2_ch3xx = 1.00e-13_r8

    !// CH3O2 + C3H7O2 --> Products
    !// IUPAC number: use ROO_41
    r_ch3o2_c3h7o2 = 1.00e-13_r8

    !// CH3O2 + CH3COD --> HO2 + CH3O + RCOHCO
    !//   CH3COD = CH3COCH2(O2)
    !//   RCOCHO = CH3COCHO
    !// IUPAC number: use ROO_24 (IUPAC06, 20080612), which
    !// use suggests channels:
    !//   a: CH3OH+RCOCHO,
    !//   b: HCHO+CH3C(O)CH2OH (i.e. not CH3XX=CH3CH(O2)CH2OH),
    !//   c:CH3O+CH3COCH2O
    !// Should clarify our choice in products.
    r_ch3o2_ch3cod = 3.80e-12_r8


    !// CH3O2 + ISOR1,ISOR2 --> HO2 + CH3O + ALDEHYD/KETONE
    !// IUPAC number: use ROO_41
    r_ch3o2_isor1 = 1.00e-13_r8
    r_ch3o2_isor2 = 1.00e-13_r8

    !//---------------------------------------------
    !// Reaction rates dependent on temperature only
    !//---------------------------------------------
    Z298 = 1._r8 / 298._r8
    do I = 1, TEMPRANGE

       TEM  = real(I + MINTEMP, r8)
       ZTEM = 1._r8 / TEM

       !// OD + M --> OP + M  (M is N2 and O2, where N2 = 0.7809M O2 = 0.2095M)
       !// JPL06, 20080612
       r_od_m(I) = &
            !// O1D + O2 -> O3P + O2
            0.2095_r8 * 3.3e-11_r8 * exp(55._r8 * ZTEM) &
            !// O1D + N2 -> O3P + N2
            + 0.7809_r8 * 2.15e-11_r8 * exp(110._r8 * ZTEM)

       !// OD + H2O -> OH + OH
       !// JPL number: A6 (JPL06, 20080617)
       r_od_h2o(I) = 1.63e-10_r8 * exp(60._r8 * ZTEM)

       !// OP + NO2 --> NO (+ O2)
       r_op_no2(I) = 5.1e-12_r8 * exp(210._r8 * ZTEM)   !JPL06, 20080612

       !// O3 + NO --> NO2 + O2
       r_o3_no(I) = 3.0e-12_r8 * exp(-1500._r8 * ZTEM)   !JPL06, 20080612


       !// O3 + NO2 --> NO3 + O2
       r_o3_no2(I) = 1.2e-13_r8 * exp(-2450._r8 * ZTEM) !JPL06, 20080612

       !// O3 + OH --> HO2 + O2
       r_o3_oh(I) = 1.7e-12_r8 * exp(-940._r8 * ZTEM)   !JPL06, 20080612

       !// O3 + HO2 --> OH + 2*O2
       r_o3_ho2(I) = 1.0e-14_r8 * exp(-490._r8 * ZTEM)   !JPL06, 20080612

       !// O3 + C2H4 --> HO2 + HCHO + CH3O2
       !//    O3 + C2H4 --> HCHO + CH2COO*                      !// IUPAC
       !//      CH2COO* -> 37% CH3COO                          !// IUPAC
       !//              -> (CO2 + H2)
       !//              -> CO + H2O
       !//              -> 12-18% OH + HCO -O2-> OH + HO2 + CO !// IUPAC
       !//    O3 + C2H4 --> HCHO + CH2COO*                   !// Atkinson 1992
       !//      CH2COO* -> ??% CH2COO
       !//              -> (??% CO2 + H2)
       !//              -> 44% CO + H2O                     !// ???
       !//              -> 12% OH + HCO -O2-> OH + HO2 + CO !// Atkinson 1992
       r_o3_c2h4(I) = 1.2e-14_r8 * exp(-2630._r8 * ZTEM) !JPL06, 20080612

       !// O3 + C3H6 --> 1.5 HO2 + O.5*(CH3CHO + HCHO + CH3)
       r_o3_c3h6(I) = 6.5e-15_r8 * exp(-1900._r8 * ZTEM) !JPL06, 20080612


       !// OH + H2O2 --> HO2 (+ H2O)
       r_oh_h2o2(I) = 1.8e-12_r8 !JPL06, 20080612

       !// OH + HO2 --> H2O + O2
       r_oh_ho2(I) = 4.8e-11_r8 * exp(250._r8 * ZTEM) !JPL06, 20080612

       !// OH + H2 -->  H2O + H (cbr121205 OH+H2 in tropospheric chemistry)
       r_oh_h2(I)  = 2.8e-12_r8 * exp(-1800._r8 * ZTEM)    !JPL06, 20080612

       !// OH + HO2NO2 --> NO2 + H2O + O2
       r_oh_ho2no2(I) = 1.3e-12_r8 * exp(380._r8 * ZTEM) !JPL06, 20080612

       !// OH + CH4 --> CH3  + H2O
       r_oh_ch4(I) = 2.45e-12_r8 * exp(-1775._r8 * ZTEM) !JPL10, 20140610
       !// Non-standard treatment: Will change CH4 by ~0.5%
       !r_oh_ch4(I) = 2.8e-14_r8 * (TEM**0.667_r8) * exp(-1575._r8 * ZTEM) !JPL10, 20140610

       !// OH + CH3O2H --> CH3O2 (+ H2O)     A   !JPL06: 70%
       !// OH + CH3O2H --> CH2O + OH (+ H2O) B   !JPL06: 30%
       !// JPL number: D16 (JPL06, 20080612)
       r_oh_ch3o2h_a(I) = 0.7_r8 * 3.8e-12_r8 * exp(200._r8 * ZTEM)
       r_oh_ch3o2h_b(I) = 0.3_r8 * 3.8e-12_r8 * exp(200._r8 * ZTEM)

       !// OH + PAN --> HCHO + NO2 + (H2O + CO2)
       !//   PAN = CH3C(O)O2NO2
       !// IUPAC number: HOx_VOC44 (IUPAC06, 20080612)
       r_oh_pan(I)  = 3.0e-14_r8


       !// OH + C2H6 (+O2) --> C2H5O2 (+ H2O)
       r_oh_c2h6(I)  = 8.7e-12_r8 * exp(-1070._r8 * ZTEM) !JPL06, 20080612

       !// OH + C3H8 --> C3H7O2 + H2O
       r_oh_c3h8(I)  = 8.7e-12_r8 * exp(-615._r8 * ZTEM) !JPL06, 20080612


       !// OH + C4H10 (+ O2) --> C4H9O2 (+ H2O)
       r_oh_c4h10(I)  = 9.1e-12_r8 * exp(-405._r8 * ZTEM) !IUPAC06, 20080612

       !// OH + ISOPREN --> ISOR1
       !//   ISOPREN = CH2C(CH3)CHCH2
       !//   Reaction path: addition of OH then addition of O2, so that
       !//   ISOR1 = C5H9(O)O2
       r_oh_isoprene(I)  = 2.7e-11_r8 * exp(390._r8 * ZTEM) !IUPAC06, 20080612

       !// OH + CH3OH --> CH3O + H2O  -O2-> CH2O + HO2
       !//            --> CH2OH + H2O -O2-> CH2O + HO2
       r_oh_ch3oh(I)  = 2.9e-12_r8 * exp(-345._r8 * ZTEM) !JPL17 (10-6)


       !// OH + HCHO --> CHO + (H2O)
       r_oh_ch2o(I)  = 5.5e-12_r8 * exp(125._r8 * ZTEM) !JPL06, 20080612

       !// OH + CH3CHO --> CH3CO (+ H2O)
       !// IUPAC number: HOx_VOC12 (IUPAC06, 20080612)
       r_oh_ch3cho(I)  = 4.4e-12_r8 * exp(365._r8 * ZTEM)


       !// OH + RCOHCO --> CH3CO+CO
       !// Use IUPAC where R=CH3, i.e. OH + CH3C(O)CHO
       !// Taken from IUPAC HOx_VOC16 (20140612)
       r_oh_rcohco(I) = 1.9e-12_r8 * exp(575._r8 * ZTEM)

       !// OH + ISOK --> ISOR2
       !//   ISOK = C4H6O (methylvinylketone + Methacrolein       )
       !//                (but-3-en-2-one    + 2-Methyl-2-propenal)
       !//                (ka                + kb                 )
       !// IUPAC06 HOx_VOC15: ka = 8.0 x 10-12 exp(380/T)
       !// IUPAC06 HOx_VOC21: kb = 2.6 x 10-12 exp(610/T)
       !// Will use mean of these (20140612)
       r_oh_isok(I) = 0.5_r8 * (  8.0e-12_r8 * exp(380._r8 * ZTEM) &
                                + 2.6e-12_r8 * exp(610._r8 * ZTEM) )


       !// OH + ACETON  (ACETON = CH3COCH3)
       !// IUPAC number: HOx_VOC19 (IUPAC06, 20080612)
       r_oh_aceton(I)  = 8.8e-12_r8 * exp(-1320._r8 * ZTEM) &
                       + 1.7e-14_r8 * exp( 423._r8 * ZTEM) !IUPAC06, 20080612

       !// OH + CH3COX --> CH3COB
       !//   CH3COX = CH3C(O)CH2CH3
       !//   CH3COB = CH3COCH(O2)CH3
       r_oh_ch3cox(I)  = 1.3e-12_r8 * exp(-25._r8 * ZTEM) !IUPAC06, 20080612

       !// OH + DMS --> H2O + CH3SCH2 (abstraction)    JPL (1997), I15
       r_oh_dms_a(I) = 1.2e-11_r8 * exp(-260._r8 * ZTEM)
       !!!!  = 1.1e-11_r8 * exp(-240._r8 * ZTEM) !JPL06, 20080612

       !// OH + H2S --> --> SO2
       r_oh_h2s(I) = 6.e-12_r8 * exp(-75._r8 * ZTEM) !JPL (1997)

       !// OH + NH3 --> NH2 + H2O
       !r_oh_nh3(I) = 1.7e-12_r8 * exp(-710._r8 * ZTEM) !JPL (2011)

       !// OH + COS --> --> SO2
       !r_oh_cos(I) = 1.1e-13_r8 * exp(-1200._r8 * ZTEM) !JPL (1997)

       !// O3P + COS --> CO + SO --> SO2
       !r_op_cos(I) = 2.1e-11_r8 * exp(-2200._r8 * ZTEM) !JPL (1997)

       !// OH + CS2 --> --> SO2 + COS                see note I13, JPL (1997)
       !r_oh_cs2(I) = 0._r8 !JPL (1997)


       !// NO + HO2 --> NO2 + OH
       r_no_ho2(I)  = 3.5e-12_r8 * exp(250._r8 * ZTEM) !JPL06, 20080612

       !// NO + CH3O2 --> NO2 + CH3O
       !// JPL number: D51 (JPL06, 20080612)
       !// (IUPAC number: ROO_1)
       r_no_ch3o2(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + C2H5O2 --> fa4323*( NO2 + (1-fb4323)*(CH2O + CH3)
       !//                         + fb4323*(HO_2 + CH3CHO) )
       !// IUPAC number: ROO_2 (could also use JPL2010 number: D59/D8)
       !// Updated 20140612
       r_no_c2h5o2(I) = 2.55e-12_r8 * exp(380._r8 * ZTEM)

       !// NO + C4H9O2 --> fa4324*( NO2 + (1-fb4324(CH3CHO + C2H5O2)
       !//                          + fb4324*(HO2 + CH3COX) )
       !// IUPAC number: use nC3H7O2, ROO_4 (IUPAC06, 20080612)
       r_no_c4h9o2(I)  = 2.9e-12_r8 * exp(350._r8 * ZTEM)

       !// NO + C6H13O2 --> NO2 + C4H9O2 + CH3CHO
       !// IUPAC number: use nC3H7O2, ROO_4 (IUPAC06, 20080612)
       r_no_c6h13o2(I)  = 2.9e-12_r8 * exp(350._r8 * ZTEM)

       !// NO + AR1 --> NO2 + AR2 + HO2 + RCOHCO  (equal to R4322)
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_ar1(I)  =  2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + AR3 --> NO2 + HO2 + HCOHCO + RCOHCO (equal to R4322)
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_ar3(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + ISOR1 --> NO2 + ISOK + HO2 + HCHO   (equal to R4322)
       !//   ISOR1=C5H9(O)O2, so an additional O2 must be consumed after
       !//   this reaction.
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_isor1(I) = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + ISOR2 --> NO2 + HO2 + HCHO + RCOHCO  ( equal to R4322)
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_isor2(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)


       !// NO + CH3COB --> NO2 + HO2 + CH3COY   (equal to R4322)
       !//   CH3COB = CH3COCH(O2)CH3
       !//   CH3COY = CH3COCOCH3
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_ch3cob(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + CH3X --> NO2 + CH3 (+ O2)
       !//   CH3X = CH3C(O)O2 
       !// JPL number: D59 (JPL06, 20080612)
       !// (IUPAC number: ROO_7)
       r_no_ch3x(I) = 8.1e-12_r8 * exp(270._r8 * ZTEM)

       !// NO + CH3COD --> NO2 + HO2 + CH3COHCO   (equal to R4322)
       !//   CH3COD = CH3COCH2(O2)
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_ch3cod(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + CH3XX --> NO2 + HCHO + HO2 + CH3CHO   (equal to R4322)
       !//   CH3XX = CH3CH(O2)CH2OH
       !// JPL number: use D51, equal to R4322 (JPL06, 20080612)
       r_no_ch3xx(I)  = 2.8e-12_r8 * exp(300._r8 * ZTEM)

       !// NO + NO3 --> 2 NO2
       r_no_no3(I)  = 1.5e-11_r8 * exp(170._r8 * ZTEM) !JPL06, 20080612

       !// NO + NH2 --> N2 + H2O
       !r_no_nh2(I) = 4.0e-12_r8 * exp(450._r8 * ZTEM) !JPL (2011)

       !// NO3 + DMS --> HNO3 + .......
       r_no3_dms(I)  = 1.9e-13_r8 * exp(520._r8 * ZTEM) !IUPAC06, 20080612

       !// NO3 + CH3CHO --> CH3CO  + HNO3
       !// JPL number: D41 (JPL06, 20080612)
       r_no3_ch3cho(I)  = 1.4e-12_r8 * exp(-1900._r8 * ZTEM)


       !// NO2 + NO3 --> NO + NO2 + O2
       r_no2_no3_b(I) = 4.5e-14_r8 * exp(-1260._r8 * ZTEM) !JPL06, 20080612

       !// NO2 + NH2 --> N2O + H2O
       !r_no2_nh2(I) = 2.1e-12_r8 * exp(650._r8 * ZTEM) !JPL (2011)

       !// HO2 + CH3O2 --> CH3O2H + (O2)
       !// JPL number: D51 (JPL06, 20080617)
       r_ho2_ch3o2(I) = 4.1e-13_r8 * exp(750._r8 * ZTEM) !JPL06, 20080612


       !// HO2 + CH3X --> CH3O2H
       !//   CH3X = CH3C(O)O2
       !// IUPAC number: HOx_VOC54 (IUPAC06, 20080612)
       r_ho2_ch3x(I) = 4.3e-13_r8 * exp(1040._r8 * ZTEM)

       !// CH3O2 + CH3O2 --> HCHO  + (CH3OH + O2)     60 %
       !//               --> 2CH3O + (O2)             40 %
       r_ch3o2_ch3o2(I)  = 9.5e-14_r8 * exp(390._r8 * ZTEM) !JPL06, 20080612

       !// CH3O2 + CH3X --> CH3O + CH3 + (CO2 + O2)  R2237A
       !//              --> HCHO + (CH3COOH + O2)    R2237B
       !//   CH3X = CH3C(O)O2
       !// JPL number: D52 (JPL06, 20080612)
       !// IUPAC number: ROO_23
       FC = 2.0e-12_r8 * exp(500._r8 * ZTEM)
       r_ch3o2_ch3x_a(I) = 0.9_r8 * FC
       r_ch3o2_ch3x_b(I) = 0.1_r8 * FC


       !// CH3X + CH3X --> CH3O2H
       !//   CH3X = CH3C(O)O2
       r_ch3x_ch3x(I) = 2.9e-12_r8 * exp(500._r8 * ZTEM) !JPL06, 20080612

       !// HO2 + RADICAL --> CH3O2H
       !//   RADICAL = CH3XX+C6H13O2+C4H9O2+CH3COB+C2H5O2+ISOR1+ISOR2+CH3COD
       !//   CH3XX  = CH3CH(O2)CH2OH
       !//   CH3COB = CH3COCH(O2)CH3
       !//   CH3COD = CH3COCH2(O2)
       !// JPL number: use D36 (HO2 + C2H5O2) (JPL06, 20080612)
       !// (IUPAC number: HOx_VOC53)
       r_ho2_radical(I)  = 7.5e-13_r8 * exp(700._r8 * ZTEM)

       !// CH3O + O2 --> HO2 + HCHO
       r_ch3o_o2(I) = 3.9e-14_r8 * exp(-900._r8 * ZTEM) !JPL06, 20080612



       !// SOA
       !// Reaction rates for precursor hydrocarbon oxidation:
       !// The rates and expressions for temperature dependence come from
       !// Chung and Seinfeld 2002 (JGR, doi:10.1029/2001J_R801397)

       !// rate for class one hydrocarbon reaction with ozone:
       r_o3_soaC1(I) = 56.e-18_r8 * exp(-732._r8 * (ZTEM - Z298))
       !// rate for class one hydrocarbon reaction with OH:
       r_oh_soaC1(I) = 84.e-12_r8 * exp(400._r8 * (ZTEM - Z298))
       !// rate for class one hydrocarbon reaction with NO3:
       r_no3_soaC1(I) = 7.e-12_r8 * exp(490._r8 * (ZTEM - Z298))

       !// rate for class two hydrocarbon reaction with ozone:
       r_o3_soaC2(I) = 200.e-18_r8 * exp(-732._r8 * (ZTEM - Z298))
       !// rate for class two hydrocarbon reaction with OH:
       r_oh_soaC2(I) = 171.e-12_r8 * exp(400._r8 * (ZTEM - Z298))
       !// rate for class two hydrocarbon reaction with NO3:
       r_no3_soaC2(I) = 12.e-12_r8 * exp(490._r8 * (ZTEM - Z298))

       !// rate for class three hydrocarbon reaction with ozone:
       r_o3_soaC3(I) = 7700.e-18_r8 * exp(-732._r8 * (ZTEM - Z298))
       !// rate for class three hydrocarbon reaction with OH:
       r_oh_soaC3(I) = 255.e-12_r8 * exp(400._r8 * (ZTEM - Z298))
       !// rate for class three hydrocarbon reaction with NO3:
       r_no3_soaC3(I) = 89.e-12_r8 * exp(490._r8 * (ZTEM - Z298))

       !// rate for class four hydrocarbon reaction with ozone:
       r_o3_soaC4(I) = 423.e-18_r8 * exp(-732._r8 * (ZTEM - Z298))
       !// rate for class four hydrocarbon reaction with OH:
       r_oh_soaC4(I) = 199.e-12_r8 * exp(400._r8 * (ZTEM - Z298))
       !// rate for class four hydrocarbon reaction with NO3:
       r_no3_soaC4(I) = 15.e-12_r8 * exp(490._r8 * (ZTEM - Z298))

       !// rate for class five hydrocarbon reaction with ozone:
       r_o3_soaC5(I) = 11650.e-18_r8 * exp(-732._r8 * (ZTEM - Z298))
       !// rate for class five hydrocarbon reaction with OH:
       r_oh_soaC5(I) = 245.e-12_r8 * exp(400._r8 * (ZTEM - Z298))
       !// rate for class five hydrocarbon reaction with NO3:
       r_no3_soaC5(I) = 27.e-12_r8* exp(490._r8 * (ZTEM - Z298))

       !// rate for xylene and O3, Tsigaridis and Kanakidou (2003),
       !// average of rates for ortho, meta and para-isomers
       r_o3_soaC7(I) = ( 2.40e-13_r8*exp(-5586._r8 * ZTEM)  &
                       + 5.37e-13_r8*exp(-6039._r8 * ZTEM) & 
                       + 1.91e-13_r8*exp(-5586._r8 * ZTEM) ) / 3._r8

       !// rate for xylene and OH, Tsigaridis and Kanakidou (2003)
       r_oh_soaC7(I) = 1.72e-11_r8

       !//rate for xylene and NO3, Tsigaridis and Kanakidou (2003)
       r_no3_soaC7(I) = 3.54e-16_r8

       !// rate for toluene and O3, Tsigaridis and kanakidou 2003
       r_o3_soaC8(I) = 2.34e-12_r8 * exp(-6694._r8*ZTEM)

       !// rate for toluene and OH, Tsigaridis and Kanakidou (2003)
       r_oh_soaC8(I) = 5.96e-12_r8

       !// rate for toluene and NO3, Tsigaridis and Kanakidou (2003)
       r_no3_soaC8(I) = 6.8e-17_r8

       !// The rate of benzene reaction with OH (Semadeni et al. 1995, 
       !// Int. J. Chem. Kinet. vol.27 pg 287-304), measured over 274-363K.
       r_oh_benzene(I) = 2.57e-12_r8 * exp(-1.92_r8 * ZTEM)


    end do !// do I = 1, TEMPRANGE



    !// --------------------------------------------------------------------
    !// Stratospheric reactions
    !// --------------------------------------------------------------------

    !// H + HO2 -> 2OH      (a)
    !//         -> O + H2O  (b)
    !//         -> H2 + O2  (c)
    !// JPL number: B5 (JPL06, 20110222)
    r_h_ho2_a = 7.2e-11_r8
    r_h_ho2_b = 1.6e-12_r8
    r_h_ho2_c = 6.9e-12_r8

    !// OP + N2O5 -> prod.     < 3.0E-16 (omitted)
    !// JPL number: C3 (JPL06, 20080617)
    !r_op_n2o5 = 0._r8

    !// OP + HNO3 -> OH + NO3  < 3.E-17
    !// JPL number: C4 (JPL06, 20080617)
    r_op_hno3 = 0._r8

    !// OD + CFC11 -> prod.
    !// CFC11 = CFCl3
    !// JPL number: A19 (JPL06, 20080616)
    r_od_cfc11_a = 0.75_r8 * 2.3e-10_r8   !JPL06 quenching 75%

    !// OD + CFC12(CCl2X) -> prod.
    !// CFC12 = CF2Cl2
    !// JPL number: A20 (JPL06, 20080616)
    r_od_cfc11_b = 0.85_r8 * 1.4e-10_r8   !JPL06 quenching 85%


    !// Reaction constants for HCFCs and OD are from JPL, 1992 /MtR, 8/9/93.
    !// JPL number: A42 (JPL06, 20080617)
    r_od_hcfc123 = 2.0e-10_r8
    !// JPL number: A36 (JPL06, 20080617)
    r_od_hcfc141 = 2.6e-10_r8
    !// JPL number: A37 (JPL06, 20080617)
    r_od_hcfc142 = 2.2e-10_r8


    !// OD + HCl -> prod.
    !// JPL number: A12 (JPL06, 20080617)
    r_od_hcl = 1.5e-10_r8


    !// H2O + ClONO2 -> prod.  <2.0E-21 (JPL06, 20080617)
    !// JPL number: F47 (JPL06, 20080617)
    r_h2o_clono2 = 0._r8
    !// HCl + ClONO2 -> prod.  <1.0E-20 (JPL06, 20080617)
    !// JPL number: F51 (JPL06, 20080617)
    r_hcl_clono2 = 0._r8
    !// BrO + OH -> HBr + O2, see r_oh_bro_b for HO2 + Br branch
    !// JPL number: G6 (JPL06, 20080617)
    !JPL06: HBr yield less than 3%, using zero, 20080617
    r_oh_bro_a = 0._r8

    do I = 1, TEMPRANGE

       TEM  = real(I + MINTEMP, r8)
       ZTEM = 1._r8 / TEM

       !// N + NO = N2 + O
       !// JPL number: C18 (JPL06, 20080617)
       r_n_no(I) = 2.1e-11_r8 * exp(100._r8 * ZTEM)

       !// N + O2 = NO + O
       !// JPL number: C16 (JPL06, 20080617)
       r_n_o2(I) = 1.5e-11_r8 * exp(-3600._r8 * ZTEM)

       !// OP + O3 -> O2 + O2
       !// JPL number: A1 (JPL06, 20080617)
       r_op_o3(I) = 8.0e-12_r8 * exp(-2060._r8 * ZTEM)


       !// OP + ClO -> Cl + O2
       !// JPL number: F1 (JPL06, 20080617)
       r_op_clo(I) = 2.8e-11_r8 * exp(85._r8 * ZTEM)

       !// OP + HCl -> OH + Cl
       !// JPL number: F4 (JPL06, 20080617)
       r_op_hcl(I) = 1.0e-11_r8 * exp(-3300._r8 * ZTEM)

       !// OP + ClONO2 -> prod. (ClO + NO3)
       !// JPL number: F6 (JPL06, 20080617)
       r_op_clono2(I) = 2.9e-12_r8 * exp(-800._r8 * ZTEM)

       !// OP + OH -> O2 + H  
       !// JPL number: B1 (JPL06, 20080617)
       r_op_oh(I) = 2.2e-11_r8 * exp(120._r8 * ZTEM)

       !// OP + HO2 -> OH + O2
       !// JPL number: B2 (JPL06, 20080617)
       r_op_ho2(I) = 3.0e-11_r8 * exp(200._r8 * ZTEM)

       !// OP + CH2O -> prod.
       !// JPL number: D4 (JPL06, 20080617)
       r_op_ch2o(I) = 3.4e-11_r8 * exp(-1600._r8 * ZTEM)



       !// OD + N2 -> OP + N2
       !// JPL number: A7 (JPL06, 20080617)
       r_od_n2(I) = 2.15e-11_r8 * exp(110._r8 * ZTEM)

       !// OD + O2 -> OP + O2
       !// JPL number: A3 (JPL06, 20080617)
       r_od_o2(I) = 3.3e-11_r8 * exp(55._r8 * ZTEM)

       !// OD + N2O -> N2 + O2
       !// JPL number: A8 (JPL06, 20080617)
       r_od_n2o_a(I) = 4.7e-11_r8 * exp(20.0_r8 * ZTEM)

       !// OD + N2O -> 2 NO 
       !// JPL number: A7 (JPL06, 20080617)
       r_od_n2o_b(I) = 6.7e-11_r8 * exp(20.0_r8 * ZTEM)


       !// O3 + H -> OH + O2
       !// JPL number: B4 (JPL06, 20080617)
       r_o3_h(I) = 1.4e-10_r8 * exp(-470._r8 * ZTEM)

       !// O3 + Cl -> ClO + O2  
       !// JPL number: F52 (JPL06, 20080617)
       r_o3_cl(I) = 2.3e-11_r8 * exp(-200._r8 * ZTEM)

       !// NO + ClO -> NO2 + Cl
       !// JPL number: F111 (JPL06, 20080617)
       r_no_clo(I) = 6.4e-12_r8 * exp(290._r8 * ZTEM)


       !// OH + OH -> H2O + OP
       !// JPL number: B9 (JPL06, 20080617)
       r_oh_oh(I) = 1.8e-12_r8 * exp(0._r8 * ZTEM)

       !// OH + HNO2 -> H2O + NO2 
       !// r_oh_hno2(I) = 0._r8 !// HNO2 chemistry not included


       !// OH + ClO -> Cl + HO2
       !// JPL number: F10 (JPL06, 20080617)
       r_oh_clo_a(I) = 7.4e-12_r8 * exp(270._r8 * ZTEM)

       !// OH + ClO -> HCl + O2
       !// JPL number: F10 (JPL06, 20080617)
       r_oh_clo_b(I) = 6.0e-13_r8 * exp(230._r8 * ZTEM)

       !// OH + HOCl -> H2O + ClO
       !// JPL number: F13 (JPL06, 20080617)
       r_oh_hocl(I) = 3.0e-12_r8 * exp(-500._r8 * ZTEM)

       !// OH + HCl -> H2O + Cl
       !// JPL number: F12 (JPL06, 20080617)
       r_oh_hcl(I) = 2.6e-12_r8 * exp(-350._r8 * ZTEM)

       !// OH + CH3Cl -> CH2Cl + H2O
       !// JPL number: F16 (JPL06, 20080617)
       r_oh_ch3cl(I) = 2.4e-12_r8 * exp(-1250._r8 * ZTEM)

       !// OH + MCF -> CH2CCl3 + H2O   (MCF=CH3CCl3)
       !// (Also used for lifetime calculations in troposphere.)
       !// JPL number: F26 (JPL06, 20080617)
       r_oh_ch3ccl3(I) = 1.64e-12_r8 * exp(-1520._r8 * ZTEM)

       !// OH + CHClF2 -> CF2Cl + H2O (CHClF2 = HCFC22)
       !// JPL number: F22 (JPL06, 20080617)
       r_oh_chclf2(I) = 1.05e-12_r8 * exp(-1600._r8 * ZTEM)

       !// OH + ClONO2 -> prod.
       !//             -> OHCl + NO3  (A)
       !//             -> HO2 + ClONO (B)
       !//             -> HONO2 + ClO (C)
       !// There is no studies on the actual products: Assume all is A,
       !// and call B+C for B.
       !// JPL number: F15 (JPL06, 20080617)
       r_oh_clono2_a(I) = 1.2e-12_r8 * exp(-330._r8 * ZTEM)
       r_oh_clono2_b(I) = 0._r8



       !// HO2 + Cl -> HCl + O2 (A)
       !// HO2 + Cl -> OH + ClO (B)
       !// JPL number: F45 (JPL06, 20080617)
       r_ho2_cl_a(I) = 1.8e-11_r8 * exp(170._r8 * ZTEM)
       r_ho2_cl_b(I) = 4.1e-11_r8 * exp(-450._r8 * ZTEM)

       !// HO2 + ClO -> HOCl + O2  
       !// JPL number: F46 (JPL06, 20080617)
       r_ho2_clo(I) = 2.7e-12_r8 * exp(220._r8 * ZTEM)


       !// H2 + Cl -> HCl + H
       !// JPL number: F53 (JPL06, 20080617)
       r_cl_h2(I) = 3.05e-11_r8 * exp(-2270._r8 * ZTEM)

       !// H2O2 + Cl -> HCl + HO2
       !// JPL number: F54 (JPL06, 20080617)
       r_cl_h2o2(I) = 1.1e-11_r8 * EXP(-980._r8 * ZTEM)

       !// Cl + CH4 -> HCl + CH3
       !// JPL number: F59 (JPL06, 20080617)
       r_cl_ch4(I) = 7.3e-12_r8 * exp(-1280._r8 * ZTEM)

       !// Cl + CH3OH -> CH2OH + HCl -O2-> CH2O + HO2 + HCl
       r_cl_ch3oh(I) = 5.5e-11_r8 ! no temp dependence so far

       !// Cl + CH2O -> HCl + HCO  
       !// JPL number: F61 (JPL06, 20080617)
       r_cl_ch2o(I) = 8.1e-11_r8 * exp(-30._r8 * ZTEM)

       !// Cl + ClONO2 -> prod.
       !// JPL number: F85 (JPL06, 20080617)
       r_cl_clono2(I) = 6.5e-12_r8 * exp(135._r8 * ZTEM)


       !// ClO + CO -> prod.
       !// JPL number: F114 (JPL06, 20080617)
       r_clo_co(I) = 1.0e-12_r8 * exp(-3700._r8 * ZTEM)


       !// Br + O3 -> BrO + O2
       !// JPL number: G31 (JPL06, 20080617)
       r_br_o3(I) = 1.7e-11_r8 * exp(-800._r8 * ZTEM)

       !// Br + H2O2 -> HBr + HO2
       !// JPL number: G32 (JPL06, 20080617)
       r_br_h2o2(I) = 1.0e-11_r8 * exp(-3000._r8 * ZTEM)

       !// Br + CH2O -> HBr + HCO
       !// JPL number: G34 (JPL06, 20080617)
       r_br_ch2o(I) = 1.7e-11_r8 * exp(-800._r8 * ZTEM)

       !// Br + HO2 -> HBr + O2
       !// JPL number: G24 (JPL06, 20080617)
       r_br_ho2(I) = 4.8e-12_r8 * exp(-310._r8 * ZTEM)


       !// BrO + HO2 -> prod.
       !// JPL number: G25 (JPL06, 20080617)
       r_bro_ho2(I) = 4.5e-12_r8 * exp(460._r8 * ZTEM)

       !// BrO + OP -> Br + O2
       !// JPL number: G1 (JPL06, 20080617)
       r_op_bro(I) = 1.9e-11_r8 * exp(230._r8 * ZTEM)

       !// BrO + ClO = Br + OClO
       !// JPL number: G41 (JPL06, 20080617)
       r_bro_clo_a(I) = 9.5e-13_r8 * exp(550._r8 * ZTEM)

       !// BrO + ClO = BrCl + O2
       !// JPL number: G41 (JPL06, 20080617)
       r_bro_clo_b(I) = 4.1e-13_r8 * exp(290._r8 * ZTEM)

       !// BrO + ClO = Br + ClOO
       !// JPL number: G41 (JPL06, 20080617)
       r_bro_clo_c(I) = 2.3e-12_r8 * exp(260._r8 * ZTEM)

       !// BrO + NO -> NO2 + Br
       !// JPL number: G39 (JPL06, 20080617)
       r_no_bro(I) = 8.8e-12_r8 * exp(260._r8 * ZTEM)

       !// BrO + BrO -> 2Br + O2
       !// JPL number: G42 (JPL06, 20080617)
       r_bro_bro_a(I) = 2.4e-12_r8 * exp(40._r8 * ZTEM)

       !// BrO + BrO -> Br2 + O2
       !// JPL number: G42 (JPL06, 20080617)
       r_bro_bro_b(I) = 2.8e-14_r8 * exp(860._r8 * ZTEM)

       !// BrO + O3 -> Br + 2O2
       !// JPL number: G38 (JPL06, 20080617)
       r_o3_bro(I) = 1.0e-12_r8 * exp(-3200._r8 * ZTEM)

       !// BrO + OH -> prod.= Br + HO2
       !// JPL number: G6, see also C51J (JPL06, 20080617)
       r_oh_bro_b(I) = 1.7e-11_r8 * exp(250._r8 * ZTEM)

       !// OH + HBr -> H2O + Br
       !// JPL number: G7 (JPL06, 20080617)
       r_oh_hbr(I) = 5.5e-12_r8 * exp(200._r8 * ZTEM)

       !// Br2 + OH -> HOBr + Br
       !// JPL number: G5 (JPL06, 20080617)
       r_oh_br2(I) = 2.1e-11_r8 * exp(240._r8 * ZTEM)

       !// O + HBr -> OH + Br
       !// JPL number: G2 (JPL06, 20080617)
       r_op_hbr(I) = 5.8e-12_r8 * exp(-1500._r8 * ZTEM)

       !// OClO + NO -> NO2 + ClO
       !// JPL number: F48 (JPL06, 20080617)
       r_no_oclo(I) = 2.5e-12_r8 * exp(-600._r8 * ZTEM)

       !// OClO + OP -> ClO + O2
       !// JPL number: F2 (JPL06, 20080617)
       r_op_oclo(I) = 2.4e-12_r8 * exp(-960._r8 * ZTEM)

       !// OH + CH3Br -> CH2Br + H2O
       !// JPL number: G8 (JPL06, 20080617)
       r_oh_ch3br(I) = 2.35e-12_r8 * exp(-1300._r8 * ZTEM)




       !// HCFCs+OH = n*Clx
       !// JPL number: F33 (JPL06, 20080617)
       r_oh_hcfc123(I) = 6.3e-13_r8 * exp(-850._r8 * ZTEM)
       !// JPL number: F27 (JPL06, 20080617)
       r_oh_hcfc141(I) = 1.25e-12_r8 * exp(-1600._r8 * ZTEM)
       !// JPL number: F28 (JPL06, 20080617)
       r_oh_hcfc142(I) = 1.3e-12_r8 * exp(-1770._r8 * ZTEM)



    end do !// do I = 1, TEMPRANGE




    !// --------------------------------------------------------------
    !// The following is data taken from Dentener and Crutzen zonally
    !// and  yearly average for the 14 first model layers (from TKB).
    !// --------------------------------------------------------------
    !// N2O5 + aq(aerosol) --> NO3-   Dentener/Crutzen, 1993 (Fig. 1)
    call set_pr42het(PR42HET)

    !// --------------------------------------------------------------------

    write(*,'(a)') '** Constant tropospheric reaction rates are set.'

    !// --------------------------------------------------------------------
  end subroutine TCRATE_CONST2
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine set_pr42het(PR42HET) !,JPAR,LPAR,LPARW)
    !// --------------------------------------------------------------------
    !// Read R42HET and put it into the correct size (LPAR,JPAR).
    !// Modified from Oslo CTM2.
    !//
    !// Reference for dataset:
    !//   Dentener/Crutzen, JGR 1993, doi:10.1029/92JD02979, Fig. 1.
    !// --------------------------------------------------------------------
    use cmn_size, only: JPAR, LPAR, LPARW, NDGRD
    use cmn_ctm, only: LMMAP, XLMMAP, YDEDG
    use cmn_met, only: HnativeRES
    use cmn_chem, only: INFILE_RES
    use cmn_parameters, only: ZPI180, A0, CPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID_Y
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// input
    real(r8),intent(out) :: PR42HET(LPAR,JPAR)
    !// Locals
    integer,parameter :: J159 = 160,L42=14
    real(r8) :: R42HET_T159(J159,L42), R42_IN(J159)
    real(r8) :: R42HET_TMP(JPAR,L42), R42_OUT(JPAR)
    real(r8) :: P42H(LPARW)
    real(r8) :: WG(J159),WGT(J159), WY(J159+1),YBEDG(J159+1)
    real(r8) :: aw(J159), delta
    integer :: J,J2,L,LL,ierr,itest
    character(len=80) :: filename
    !// --------------------------------------------------------------------

    PR42HET(:,:) = 0._r8

    if (NDGRD .eq. 1) then
       !// Simulation is done in native resolution
       !//filename = './tables/r42het_'//trim(HnativeRES)//'.dat'
       filename = INFILE_RES
       !// Test if file profided is in right resolution
       itest = SCAN(INFILE_RES, trim(HnativeRES))
       if (itest .eq. 0) then
          !// wrong reolution or option!
          call ctmExitC('*** Wrong resolution or resolution option '//trim(filename))
       end if

       open(1,File=filename,Form='FORMATTED',Status='OLD',iostat=ierr)

       if (ierr .ne. 0) then
          !// file does not exist!
          write(*,'(a)') '*** Cannot find file '//trim(filename)
       end if

       !// Have file, read it
       read(1,'(15E10.2)') R42HET_TMP
       close(1)

    else
       !// Not native resolution; must interpolate
       filename = INFILE_RES
       open(1,File=filename,Form='FORMATTED',Status='OLD',iostat=ierr)
       if (ierr .ne. 0) then
          !// file does not exist!
          write(*,'(a)') '*** Cannot find file '//trim(filename)
       end if

       !// Have file, read it
       read(1,'(15E10.2)') R42HET_T159
       close(1)

       !// Set YBEDG
       call GAUSST2(J159,WG,WGT)
       J2 = J159 / 2
       WY(1) = -1._r8
       do J = 2, J2
          WY(J) = WY(J-1) + WGT(J-1)
       end do
       WY(J2+1) = 0._r8
       do J = J2+2, J159+1
          WY(J) = -WY(J159+2-J)
       end do
       do J = 1, J159+1
          YBEDG(J) = ZPI180 * asin(WY(J))
       end do

       !// Interpolate to current resolution
       do J = 1, J159
          !// Weighting factor for T42 latitude band
          aw(J) = (sin(CPI180 * YBEDG(J+1)) - sin(CPI180 * YBEDG(J)))
       end do
       !// Need T159 grid box size
       do L = 1, L42
          !// Multiply by T159 area of boxes
          R42_IN(:) = R42HET_T159(:,L) * aw(:)

          call E_GRID_Y(R42_IN,YBEDG,J159, R42_OUT,YDEDG,JPAR)

          !// Divide by new area
          do J=1,JPAR
             delta =  (sin(CPI180 * YDEDG(J+1)) - sin(CPI180 * YDEDG(J)))
             R42HET_TMP(J,L) = R42_OUT(J) / delta
          end do
       end do

    end if

    !cmga20feb03--v   ... only a temporary fix!!! mga 29AUG2003
    do J = 1, JPAR
       do L = 1, LPAR

          Select Case(L)

          Case(01:02)
             P42H(L) = R42HET_TMP(J,01)
          Case(03:05)
             P42H(L) = R42HET_TMP(J,02)
          Case(06:08)
             P42H(L) = R42HET_TMP(J,03)
          Case(09:10)
             P42H(L) = R42HET_TMP(J,04)
          Case(11:13)
             P42H(L) = R42HET_TMP(J,05)
          Case(14:15)
             P42H(L) = R42HET_TMP(J,06)
          Case(16:18)
             P42H(L) = R42HET_TMP(J,07)
          Case(19:20)
             P42H(L) = R42HET_TMP(J,08)
          Case(21:23)
             P42H(L) = R42HET_TMP(J,09)
          Case(24:25)
             P42H(L) = R42HET_TMP(J,10)
          Case(26:28)
             P42H(L) = R42HET_TMP(J,11)
          Case(29:30)
             P42H(L) = R42HET_TMP(J,12)
          Case(31:32)
             P42H(L) = R42HET_TMP(J,13)
          Case(33:34)
             P42H(L) = R42HET_TMP(J,14)
          Case DEFAULT
             P42H(L) = 0._r8
          End Select

       end do
       do L = 1, LPARW !// Not degraded
          LL = LMMAP(L) !// Degraded 
          PR42HET(LL,J) = PR42HET(LL,J) + P42H(L) * XLMMAP(L)
       end do
    end do


    !// --------------------------------------------------------------------
  end subroutine set_pr42het
  !// ----------------------------------------------------------------------


  !// ---------------------------------------------------------------------
  subroutine getRQAER(LM,ZMID,PLAND,PQAER)
    !// -------------------------------------------------------------------
    !// Find removal rate by aerosol, only dependent upon Land/Sea and
    !// with vertical variation. Equations are given in appendix A.4 in
    !// Berntsen and Isaksen, JGR, 1997, doi: 10.1029/97JD01140.
    !//
    !// Updated for CTM3 (from CTM2).
    !//
    !// THIS PARAMETERIZATION SHOULD BE PHASED OUT!
    !// -------------------------------------------------------------------
    implicit none
    !// -------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LM !// vertical resolution
    real(r8), intent(in) :: &
         ZMID(LM), &          !// model layer heigt [m]
         PLAND                !// land fraction
    !// Output
    real(r8), intent(out) :: PQAER(LM) !// Removal due to aerosols

    !// Locals
    integer :: L
    real(r8)  :: RHLP, CFREQL, CFREQS
    !// Parameter: sticking coefficient
    real(r8), parameter :: stcoeff = 0.01_r8
    !// -------------------------------------------------------------------
    do L = 1, LM
       !// Calculate below 16km (ZGRD is in meter, need km in calculations)
       if (ZMID(L) .lt. 16.e3_r8) then
          RHLP   = exp(-ZMID(L) * 0.5e-3_r8)
          CFREQL = 1.e-2_r8 * RHLP
          CFREQS = 1.e-3_r8 * RHLP
          PQAER(L) = STCOEFF * (CFREQS * (1._r8 - PLAND) + CFREQL * PLAND)
       else
          PQAER(L) = 0._r8
       end if
    end do
    !// -------------------------------------------------------------------
  end subroutine getRQAER
  !// ---------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TCRATE_TP_IJ_TRP(ICOL,JCOL, LM,LMTROP,&
       TEMP, AIR_MOLEC, H2O_MOLEC, PMIDL, &
       r_op_o2_m, r_ho2_ho2_tot, r_oh_oh_m, r_no2_oh_m, r_no2_no3_m, &
       r_n2o5_m, r_ho2_no2_m, r_ho2no2_m, r_oh_hno3, r_oh_co_a, r_oh_co_b, &
       r_oh_c2h4_m, r_oh_c3h6_m, r_ch3_o2_m, &
       r_oh_hcohco_m_a, r_oh_hcohco_m_b, r_no2_ch3x_m, r_pan_m, &
       r_no_ho2_b, r_op_no_m, r_op_no2_m)
    !// --------------------------------------------------------------------
    !// Find reaction rates, dependent upon pressure and temperature.
    !// Slightly modified TCRATE_TP, where J,L has been removed and
    !// where a test for LMTROP is added.
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use utilities_oslo, only: RATE3B
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ICOL, JCOL, LM, LMTROP
    real(r8), intent(in) :: AIR_MOLEC(LM), &
                            H2O_MOLEC(LM), &
                            PMIDL(LM), &
                            TEMP(LM)
    !// Output
    real(r8), intent(out) :: &
         r_op_o2_m(LM), r_ho2_ho2_tot(LM), r_oh_oh_m(LM), r_no2_oh_m(LM), &
         r_no2_no3_m(LM), r_ho2_no2_m(LM), r_ho2no2_m(LM), r_oh_hno3(LM), &
         r_oh_co_a(LM), r_oh_co_b(LM), r_oh_c2h4_m(LM), r_ch3_o2_m(LM), &
         r_oh_hcohco_m_a(LM), r_oh_hcohco_m_b(LM), &
         r_no2_ch3x_m(LM), r_pan_m(LM), r_n2o5_m(LM), r_oh_c3h6_m(LM), &
         r_no_ho2_b(LM), r_op_no_m(LM), r_op_no2_m(LM)

    !// Locals
    real(r8) :: &
         KZERO, KINF, FC, &
         TZ300, ZTEM, TEM, &
         R0, R2, R3M, O2CONC
    integer :: L
    !// --------------------------------------------------------------------

    !// Initialize
    if (.true.) then
       r_op_o2_m(:)   = 0._r8
       r_oh_oh_m(:)   = 0._r8
       r_no2_oh_m(:)  = 0._r8
       r_no2_no3_m(:) = 0._r8
       r_n2o5_m(:)    = 0._r8
       r_ho2_no2_m(:) = 0._r8
       r_ho2no2_m(:)  = 0._r8
       r_oh_hno3(:)   = 0._r8
       r_oh_co_a(:)   = 0._r8
       r_oh_co_b(:)   = 0._r8
       r_oh_c2h4_m(:) = 0._r8
       r_oh_c3h6_m(:) = 0._r8
       r_ch3_o2_m(:)  = 0._r8
       r_no2_ch3x_m(:)= 0._r8
       r_pan_m(:)     = 0._r8
       r_op_no_m(:)   = 0._r8
       r_op_no2_m(:)  = 0._r8

       r_ho2_ho2_tot(:)   = 0._r8
       r_oh_hcohco_m_a(:) = 0._r8
       r_oh_hcohco_m_b(:) = 0._r8
       r_no_ho2_b(:)      = 0._r8
    end if

    !// 3-body (T,p-dependent) calculations
    !// --------------------------------------------------------------------
    !// Only calculate up to the tropopause
    do L = 1, LMTROP

       TEM   = TEMP(L)
       ZTEM  = 1._r8 / TEM
       TZ300 = TEM / 300._r8

       !// OP + O2 + M --> O3 + M
       !// JPL number: A1 (JPL06, 20080618)
       r_op_o2_m(L) = 6.0e-34_r8*(TZ300**(-2.4_r8))*AIR_MOLEC(L)

       !// OP + NO + M --> NO2
       !// JPL number: C1
       !KZERO = 6.0e-32_r8 * (TZ300**(-1.5_r8)) * AIR_MOLEC(L) !JPL10, 20140610
       !KINF  = 3.0e-11_r8 !// * (TZ300**(-0._r8))            !JPL10, 20140610
       !R3843M(L) = TROE (KZERO,KINF,0.6_r8)                   !JPL10, 20140610
       r_op_no_m(L) = rate3B(101, TZ300, AIR_MOLEC(L), &
            6.0e-32_r8, 1.5_r8, 3.0e-12_r8, 0._r8, 0.6_r8, 0)

       !// OP + NO2 + M --> NO3
       !// JPL number: C2
       !KZERO    = 2.5e-31_r8 * (TZ300**(-1.8_r8)) * AIR_MOLEC(L) !JPL10, 20140610
       !KINF     = 2.2e-11_r8 * (TZ300**(-0.7_r8))               !JPL10, 20140610
       !R3844M(L) = TROE (KZERO,KINF,0.6_r8)                   !JPL10, 20140610
       r_op_no2_m(L) = rate3B(102, TZ300, AIR_MOLEC(L), &
            2.5e-31_r8, 1.8_r8, 2.2e-11_r8, 0.7_r8, 0.6_r8, 0)

       !// OH + NO + M --> HONO
       !// JPL number: C3, can be included later
       !KZERO    = 7.0e-31_r8 * (TZ300**(-2.6_r8)) * AIR_MOLEC(L) !JPL10, 20140610
       !KINF     = 3.6e-11_r8 * (TZ300**(-0.1_r8))               !JPL10, 20140610
       !R4043M(L) = TROE (KZERO,KINF,0.6_r8)                   !JPL10, 20140610
       !R4043M(L) = rate3B(103, TZ300, AIR_MOLEC(L), &
       !     7.0e-31_r8, 2.6_r8, 3.6e-11_r8, 0.1_r8, 0.6_r8, 0)

       !// OH + OH + M --> H2O2 + M
       !// JPL number: B2 (JPL06, 20080617)
       ! Should look more closely into this reaction.
       !KZERO    = 6.9e-31_r8 * (TZ300**(-1._r8)) * AIR_MOLEC(L) !JPL06, 20080619
       !KINF     = 2.6e-11_r8                                   !JPL06, 20080619
       !R4040M(L) = TROE (KZERO,KINF,0.6_r8)                   !JPL06, 20080619
       r_oh_oh_m(L) = rate3B(104, TZ300, AIR_MOLEC(L), &
            6.9e-31_r8, 1.0_r8, 2.6e-11_r8, 0._r8, 0.6_r8, 0)


       !// NO2 + OH + M --> HNO3 + M
       !// JPL number: C4 (JPL06, 20080617)
       !KZERO    = 1.8e-30_r8*(TZ300**(-3.0_r8))*AIR_MOLEC(L)
       !KINF     = 2.8e-11_r8*(TZ300**(0.0_r8))
       !R4440(L) = TROE (KZERO,KINF,0.6_r8)
       r_no2_oh_m(L) = rate3B(105, TZ300, AIR_MOLEC(L), &
            1.8e-30_r8, 3.0_r8, 2.8e-11_r8, 0._r8, 0.6_r8, 0)

       !// NO2 + NO3 + M --> N2O5 + M
       !// JPL number: C6 (JPL06, 20080617)
       !KZERO     = 2.0e-30_r8*(TZ300**(-4.4_r8))*AIR_MOLEC(L)
       !KINF      = 1.4e-12_r8*(TZ300**(-0.7_r8))
       !R4441A(L) = TROE (KZERO,KINF,0.6_r8)
       r_no2_no3_m(L) = rate3B(106, TZ300, AIR_MOLEC(L), &
            2.0e-30_r8, 4.4_r8, 1.4e-12_r8, 0.7_r8, 0.6_r8, 0)


       !// N2O5 + M --> NO2 + NO3 + M
       !// IUPAC number: NOx32 (IUPAC02, 20080617)
       !KZERO   = 1.3e-3_r8*(TZ300**(-3.5_r8))*Exp(-11000._r8*ZTEM) * AIR_MOLEC(L)
       !KINF    = 9.7e14_r8*(TZ300**(0.1_r8))*Exp(-11080._r8*ZTEM)
       !FC      = 0.9_r8*Exp(-TEM/430._r8) + 2.5_r8*Exp(-1950._r8*ZTEM)
       !R42T(L) = TROE (KZERO,KINF,FC)
       KZERO   = 1.3e-3_r8 * exp(-11000._r8 * ZTEM)
       KINF    = 9.7e14_r8 * exp(-11080._r8 * ZTEM)
       FC      = 0.9_r8 * exp(-TEM / 430._r8) &
                 + 2.5_r8 * exp(-1950._r8 * ZTEM)
       r_n2o5_m(L) = rate3B(107, TZ300, AIR_MOLEC(L), &
            KZERO, 3.5_r8, KINF, -0.1_r8, FC, 0)

       !// NO2 + HO2 --> HO2NO2
       !// JPL number: C5 (JPL06, 20080617)
       !KZERO    = 2.0e-31_r8*(TZ300**(-3.4_r8))*AIR_MOLEC(L)
       !KINF     = 2.9e-12_r8*(TZ300**(-1.1_r8))
       !R4421(L) = TROE (KZERO,KINF,0.6_r8)
       r_ho2_no2_m(L) = rate3B(108, TZ300, AIR_MOLEC(L), &
            2.0e-31_r8, 3.4_r8, 2.9e-12_r8, 1.1_r8, 0.6_r8, 0)


       !// HO2NO2 + M --> HO2 + NO2 + M
       !// IUPAC number: NOx17 (IUPAC02, 20080617)
       !KZERO   = 4.1e-5_r8 * Exp(-10650._r8*ZTEM)*AIR_MOLEC(L)
       !KINF    = 4.8e15_r8 * Exp(-11170._r8*ZTEM)
       !R17T(L) = TROE (KZERO,KINF,0.5_r8)
       KZERO   = 4.1e-5_r8 * exp(-10650._r8 * ZTEM)
       KINF    = 4.8e15_r8 * exp(-11170._r8 * ZTEM)
       r_ho2no2_m(L) = rate3B(109, TZ300, AIR_MOLEC(L), &
            KZERO, 0._r8, KINF, 0._r8, 0.5_r8, 0)


       !// OH + CO --> HO2 (CO2)
       !//         --> HOCO --O2--> HO2 + CO2    (A)
       !//         --> H + CO2 --O2--> HO2 + CO2 (B)
       !// JPL2010 number D1: Updated 20140610.
       !// Provides two channels (A/B); updating to these
       !// makes the reactions more consistent with other reactions.
       !// A-channel: Standard 3-body calculation:
       !KZERO    = 5.9e-33_r8*(TZ300**(-1.4_r8))*AIR_MOLEC(L)
       !KINF     = 1.1e-12_r8*(TZ300**(1.3_r8))
       !R4006A(L)= TROE (KZERO,KINF,0.6_r8)
       r_oh_co_a(L)= rate3B(110, TZ300, AIR_MOLEC(L), &
            5.9e-33_r8, 1.4_r8, 1.1e-12_r8, -1.3_r8, 0.6_r8, 0)
       !// B-channel: Here KINF is divided by AIR_MOLEC, i.e. it
       !// is calculated as an activation channel (k_f^ca):
       !KZERO    = 1.5e-13_r8*(TZ300**(0.6_r8))
       !KINF     = 2.1d9*(TZ300**(6.1_r8))/AIR_MOLEC(L)
       !R4006B(L)= TROE (KZERO,KINF,0.6_r8)
       r_oh_co_b(L)= rate3B(111, TZ300, AIR_MOLEC(L), &
            1.5e-13_r8, -0.6_r8, 2.9e9_r8, -6.1_r8, 0.6_r8, 1)
       !// Old value for historical reasons:
       !// IUPAC number: HOx_VOC10 (IUPAC06, 20080617)
       !// IUPAC assume pure [N2], hence we use AIR_MOLEC.
       !// R4006A(L) = 1.44e-13_r8*(1._r8 + AIR_MOLEC(L)/4.0d19)
       !// R4006B(L) = 0._r8
       !// (Actually, IUPAC recommends division by 4.2d19.)

       !// OH + C2H4 --> (HOCH2CH2) -M-> CH3 + HCHO
       !// JPL number: D5 (JPL06, 20080617)
       !KZERO    = 1.0e-28_r8*(TZ300**(-4.5_r8))*AIR_MOLEC(L)
       !KINF     = 8.8e-12_r8*(TZ300**(-0.85_r8))
       !R4007(L) = TROE (KZERO,KINF,0.6_r8)
       r_oh_c2h4_m(L) = rate3B(112, TZ300, AIR_MOLEC(L), &
            1.0e-28_r8, 4.5_r8, 8.8e-12_r8, -0.85_r8, 0.6_r8, 0)


       !// OH + C3H6 + (O2) --> CH3XX, CH3XX = CH3CH(O2)CH2OH
       !// IUPAC number: HOx_VOC15 (IUPAC06, 20080617)
       !KZERO    = 8.0e-27_r8*(TZ300**(-3.5_r8))*AIR_MOLEC(L)
       !KINF     = 3.0e-11_r8*(TZ300**(-1.0_r8))
       !R4009(L) = TROE (KZERO,KINF,0.5_r8)
       r_oh_c3h6_m(L) = rate3B(113, TZ300, AIR_MOLEC(L), &
            8.0e-27_r8, 3.5_r8, 3.0e-11_r8, 1._r8, 0.5_r8, 0)

       !// CH3 + O2 +M --> CH3O2 + M
       !// JPL number: D3 (JPL06, 20080617)
       !// (IUPAC number: R_Oxygen1)
       !KZERO    = 4.0e-31_r8*(TZ300**(-3.6_r8))*AIR_MOLEC(L)
       !KINF     = 1.2e-12_r8*(TZ300**(1.1_r8))
       !RCH3M(L) = TROE (KZERO,KINF,0.6_r8)
       r_ch3_o2_m(L) = rate3B(114, TZ300, AIR_MOLEC(L), &
            4.0e-31_r8, 3.6_r8, 1.2e-12_r8, -1.1_r8, 0.6_r8, 0)

       !// NO2 + CH3X + M --> PAN + M
       !//   CH3X = CH3C(O)O2
       !// JPL number: D12 (JPL06, 20080617)
       !KZERO    = 9.7e-29_r8*(TZ300**(-5.6_r8))*AIR_MOLEC(L)
       !KINF     = 9.3e-12_r8*(TZ300**(-1.5_r8))
       !R4437(L) = TROE (KZERO,KINF,0.6_r8)
       r_no2_ch3x_m(L) = rate3B(115, TZ300, AIR_MOLEC(L), &
            9.7e-29_r8, 5.6_r8, 9.3e-12_r8, 1.5_r8, 0.6_r8, 0)


       !// PAN + M --> CH3X + NO2 + M
       !//   PAN = CH3C(O)O2NO2
       !// IUPAC number: ROO_15 (IUPAC03, 20080617)
       !KZERO  = 4.9e-3_r8*Exp(-12100._r8*ZTEM)*AIR_MOLEC(L)
       !KINF   = 5.4e16_r8*Exp(-13830._r8*ZTEM)
       !R5T(L) = TROE (KZERO,KINF,0.3_r8)
       KZERO  = 4.9e-3_r8 * exp(-12100._r8 * ZTEM)
       KINF   = 5.4e16_r8 * exp(-13830._r8 * ZTEM)
       r_pan_m(L) = rate3B(116, TZ300, AIR_MOLEC(L), &
            KZERO, 0._r8, KINF, 0._r8, 0.3_r8, 0)


    end do !// do L = 1, LMTROP

    !// Other T,p-dependent calculations
    !// --------------------------------------------------------------------
    do L = 1, LMTROP

       TEM   = TEMP(L)
       ZTEM  = 1._r8 / TEM
       TZ300 = TEM / 300._r8
       O2CONC    = 0.21_r8 * AIR_MOLEC(L)


       !// HO2 + HO2 --> H2O2 (+ O2)
       !// JPL number: bi-molecular B13 (JPL06, 20080617)
       R0 = 3.5e-13_r8 * exp(430._r8 * ZTEM)

       !// HO2 + HO2 + M --> H2O2 (+ O2 + M)
       !// JPL number: bi-molecular B13 (JPL06, 20080617)
       R2 = 1.7e-33_r8 * exp(1000._r8 * ZTEM) * AIR_MOLEC(L)

       ! Uten H2O effekt ::  R = R0 + R2
       ! Med H2O effekt (Derwent-reaksjonen)
       r_ho2_ho2_tot(L) = (R0 + R2) &
            * (1._r8 + 1.4_r8 * 1.e-21_r8 * H2O_MOLEC(L) &
                       * exp(2200._r8 * ZTEM))


       !// OH + HNO3 --> NO3 (+ H2O)
       !// JPL number: bi-molecular C9 (JPL06, 20080617)
       R0       = 2.4e-14_r8 * exp(460._r8 * ZTEM)
       R2       = 2.7e-17_r8 * exp(2199._r8 * ZTEM)
       R3M      = 6.5e-34_r8 * exp(1335._r8 * ZTEM) * AIR_MOLEC(L)
       r_oh_hno3(L) = R0 + R3M / (1._r8 + (R3M/R2))


       !// OH + HCOHCO --> HC(O)CO
       !//           ----> (HCO+CO)         (a)
       !//           -O2-> CH3CO3           (b) 
       !//           -O2-> (2CO + HO2)      (c)
       !// IUPAC number: HOx_VOC16 (IUPAC06, 20080617)
       !//   k = 3.1e-12_r8*exp(340/T)
       !//   ka/kb = 3.5d18
       !//   kb/kc = 1.0, i.e. kb = kc
       !// The reaction is assumed to proceed to the three channels
       !// immediately, so that:
       !//      k = ka + kb[O2] + kc[O2] = ka + 2kb[O2]
       !//        = ka + 2ka/(ka/kb)[O2] = ka(1 + 2[O2]/(ka/kb))
       !//   ka/k = 1/(1 + 2[O2]/(ka/kb))
       !//   kb/k = ka/k*kb/ka = 1/(ka/kb) * 1/(1 + 2[O2]/(ka/kb))
       !//        = 1/(ka/kb + 2[O2])
       !// which gives
       !//    ka = k * (ka/kb) / (ka/kb + 2[O2])
       !//    kb = k * 1 / (ka/kb + 2[O2])
       !//    kc = kb
       FC        = 1._r8 / (3.5e18_r8 + 2._r8 * O2CONC)
       !// Updated 20140610 to use basis reaction rate from IUPAC,
       !// instead of JPL.
       r_oh_hcohco_m_a(L) = 3.1e-12_r8 * exp(340._r8 * ZTEM) * 3.5e18_r8 * FC
       !// Also multiply with O2, to skip this in chemistry
       r_oh_hcohco_m_b(L) = 3.1e-12_r8 * exp(340._r8 * ZTEM) * FC * O2CONC
       !// Old JPL (number: D18, JPL06, 20080612)
       !// R4035A(L) = 1.15e-11_r8 * 3.5d18 * FC
       !// R4035B(L) = 1.15e-11_r8 * O2 * FC


       !// NO + HO2 -> HNO3 (LeBras)
       !// 760Torr=1013.25mb, so we should multiply with 760./1013.25
       !// Their value is originally in percent!
       !// JPL2011 (nr 17) does not recommend this until other studies
       !// confirm the reaction, so we leave it out of the standard model runs.
       !// R0 = 1.e-2_r8*( 530._r8/TEMP(L) &
       !//                + 6.4e-4_r8*760._r8/1013.25_r8*PMIDL(L) - 1.73_r8 )
       !// if (R0 .lt. 0._r8) then
       !//    R0 = 0._r8
       !// else if (R0 .gt. 1._r8) then
       !//    R0 = 0.999_r8
       !// end if
       !// r_no_ho2_b(L) = r_no_ho2(NINT(TEMP(L)) - MINTEMP) * R0
       r_no_ho2_b(L) = 0._r8

    end do !// do L = 1, LMTROP



    !// Other reaction rates not included
    !//----------------------------------
    !// OP + C2H4 --> CH3 + CHO,
    !//   R3807(I,J,L) = 7.6e-13_r8
    !// NO3 + C2H6 --> C2H5O2  + HNO3  Ref.4     Wayne at al. (1990)
    !//   R0841(I,J,L) = 5.7e-12_r8*Exp(-4426._r8*ZTEM)
    !// NO + CH2O2OH --> NO2 + HO2 + ?   NASA/JPL (1990) (assumed equal to R4322)
    !//   R4326(I,J,L) = 4.2e-12_r8*Exp(180._r8*ZTEM)

    !// --------------------------------------------------------------------
  end subroutine TCRATE_TP_IJ_TRP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TCRATE_TP_IJ_STR ( &
          LM,   TRACER_ID_MAX, LMTROP, TEMP, &
          AIR_MOLEC,PMIDL, &
          r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
          r_oh_hno3, r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
          r_cl_o2_m, r_cloo_heat, r_clo_clo_m, r_cl2o2_m, r_oh_no2_m, &
          r_op_no2_m, r_oh_oh_m, r_oh_no_m, r_no2_clo_m,  &
          r_op_o2_m, r_o2_h_m, r_no2_bro_m, &
          r_no_ho2_b )
    !// --------------------------------------------------------------------
    !// Subroutine assigns values to reaction rates (3-body and
    !// reactions dependent on temperature and/or pressure).
    !// Temperature and density fields are sendt from calling routine.
    !//
    !// 3-body reactions follows JPL equations, either as
    !// standard calculation
    !//   kf([M],T) =
    !/        k0(T)*[M]
    !//    ------------------- * 0.6^{1+[log10(k0(T)*[M]/ki(T))]^2}^(-1)
    !//    1 + k0(T)*[M]/ki(T)
    !// or as activation channel
    !//   kf^ca([M],T) =
    !//       k0(T)
    !//    -------------------- * 0.6^{1+[log10(k0(T)/(ki(T)/[M]))]^2}^(-1)
    !//    1 + k0(T)/(ki(T)/[M])
    !// but it is easily seen that
    !//   kf([M],T) = kf^ca([M],T) * [M]
    !//
    !// Ole Amund Sovde, December 2014:
    !//              Updated with new routine calculating rates,
    !//              and also 3-body rates are now multiplied
    !//              by AIR_MOLEC here instead of in chemistry.
    !//              Rates updated to JPL10-6.
    !// Ole Amund Sovde, October 2008:
    !//              Updated for Oslo CTM3; IJ-block structure.
    !//              Heterogeneous reactions moved to routine
    !//              TCRATE_HET_IJ.
    !// March 2007: Updated for JPL06
    !// --------------------------------------------------------------------
    use utilities_oslo, only: RATE3B
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// Input
    integer, intent(in) :: LM, TRACER_ID_MAX, LMTROP
    real(r8), intent(in)  :: &
         TEMP(LM), &           !// Temperature
         AIR_MOLEC(LM), &      !// AIR molecular density (molec/cm3):
         PMIDL(LM)             !// Pressure, grid box center [hPa]

    !// Output
    real(r8), dimension(LM), intent(out) :: &
         !// Gas phase rate constants:
         r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
         r_oh_hno3, r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
         r_cl_o2_m, r_cloo_heat, r_clo_clo_m, r_cl2o2_m, r_oh_no2_m, &
         r_op_no2_m, r_oh_oh_m, r_oh_no_m, r_no2_clo_m,  &
         r_op_o2_m, r_o2_h_m, r_no2_bro_m, &
         r_no_ho2_b

    !// Local variables
    real(r8)  ::  THE, TZ300, C0,C2,C3M, KZERO, KINF, R0, R2, ZTEM
    integer :: LMAX !// Uppermost chemistry level
    integer :: L !// Loop index
    !// --------------------------------------------------------------------

    !// Calculate up to layer LM or LM-1, must correspond to chemistry
    LMAX = LM

    !// Initialise
    r_ho2_no2_m(:) = 0._r8
    r_no2_no3_m(:) = 0._r8
    r_ho2no2_heat(:) = 0._r8
    r_n2o5_heat(:) = 0._r8
    r_oh_hno3(:)   = 0._r8
    r_oh_co_a(:)   = 0._r8
    r_oh_co_b(:)   = 0._r8
    r_ho2_ho2_tot(:) = 0._r8
    r_cl_o2_m(:)   = 0._r8
    r_cloo_heat(:) = 0._r8
    r_clo_clo_m(:) = 0._r8
    r_cl2o2_m(:)   = 0._r8
    r_oh_no2_m(:)  = 0._r8
    r_op_no2_m(:)  = 0._r8
    r_oh_oh_m(:)   = 0._r8
    r_oh_no_m(:)   = 0._r8
    r_no2_clo_m(:) = 0._r8
    r_op_o2_m(:)   = 0._r8
    r_o2_h_m(:)    = 0._r8
    r_no2_bro_m(:) = 0._r8
    r_no_ho2_b(:)  = 0._r8


    !// 3-body (T,p-dependent) calculations
    !// ---------------------------------------------------------------
    !// Loop from start of stratosphere up to uppermost chemistry level
    do L = LMTROP+1, LMAX

       TZ300 = TEMP(L) / 300._r8


       !// NO2 + HO2 -M-> HO2NO2
       !// JPL number: C5 (JPL06, 20080617)
       r_ho2_no2_m(L) = rate3B(301, TZ300, AIR_MOLEC(L), &
            2.0e-31_r8, 3.4_r8, 2.9e-12_r8, 1.1_r8, 0.6_r8, 0)

       !// NO2 + NO3 -M-> N2O5
       !// JPL number: C6 (JPL06, 20080617)
       r_no2_no3_m(L) = rate3B(302, TZ300, AIR_MOLEC(L), &
            2.0e-30_r8, 4.4_r8, 1.4e-12_r8, 0.7_r8, 0.6_r8, 0)

       !// Cl + O2 -M-> ClOO
       !// JPL number: F1 (JPL06, 20080617)
       r_cl_o2_m(L) = rate3B(303, TZ300, AIR_MOLEC(L), &
            2.2e-33_r8, 3.1_r8, 1.8e-10_r8, 0._r8, 0.6_r8, 0)

       !// ClO + ClO -M-> Cl2O2
       !// JPL number: F10 (JPL10-6, 20110610)
       r_clo_clo_m(L) = rate3B(304, TZ300, AIR_MOLEC(L), &
            1.6e-32_r8, 4.5_r8, 3.0e-12_r8, 2._r8, 0.6_r8, 0)

       !// OH + CO -> CO2 + H
       !// There are two channels:
       !//         --M--> HOCO    --O2--> HO2 + CO2 (A)
       !//         --M--> H + CO2 --O2--> HO2 + CO2 (B)
       !//       As in old code, we assume (A) to also give H, not HO2.
       !// IUPAC number: HOx_VOC10 (IUPAC06)
       !// IUPAC uses N2, but it seems the experiments use pure N2, not
       !// air, so AIR_MOLEC should be applied when using this.
       !// There is only one channel for IUPAC
       !//   C1945A(L) = 0._r8
       !//   C1945B(L) = 1.44d-13 * (  1._r8 + AIR_MOLEC(L) / 4.0d19  )
       !// JPL2010 number D1: Updated 20140610.
       !// Provides two channels (A/B); updating to these
       !// makes the reactions more consistent with other reactions.
       !// A-channel: Standard 3-body calculation:
       r_oh_co_a(L) = rate3B(305, TZ300, AIR_MOLEC(L), &
            5.9e-33_r8, 1.4_r8, 1.1e-12_r8, -1.3_r8, 0.6_r8, 0)

       !// B-channel: Calculated as an activation channel (k_f^ca):
       r_oh_co_b(L) = rate3B(306, TZ300, AIR_MOLEC(L), &
            1.5e-13_r8,-0.6_r8, 2.1e9_r8, -6.1_r8, 0.6_r8, 1)


       !// NO2 + OH + M -> HNO3 + M
       !// JPL number: C4 (JPL06, 20080617)
       r_oh_no2_m(L) = rate3B(307, TZ300, AIR_MOLEC(L), &
            1.8e-30_r8, 3._r8, 2.8e-11_r8, 0._r8, 0.6_r8, 0)

       !// OP + NO2 + M-> NO3 + M
       !// JPL number: C2 (JPL06, 20080617)
       r_op_no2_m(L) = rate3B(308, TZ300, AIR_MOLEC(L), &
            2.5e-31_r8, 1.8_r8, 2.2e-11_r8, 0.7_r8, 0.6_r8, 0)

       !// OH + OH + M -> H2O2 + M
       !// JPL number: B2 (JPL06, 20080617)
       r_oh_oh_m(L) = rate3B(309, TZ300, AIR_MOLEC(L), &
            6.9e-31_r8, 1._r8, 2.6e-11_r8, 0._r8, 0.6_r8, 0)
       !C1919M(L) = C1919M(L) / AIR_MOLEC(L)

       !// NO + OH + M -> HONO + M
       !// JPL number: C3 (JPL06, 20080617)
       !C919M(L) = rate3B(310, TZ300, AIR_MOLEC(L), &
       !     7.0d-31, 2.6_r8, 1.5d-11, 0.5_r8, 0.6_r8, 0)
       r_oh_no_m(L) = 0._r8


       !// NO2 + ClO + M -> ClONO2 + M
       !// JPL number: F8 (JPL06, 20080617)
       r_no2_clo_m(L) = rate3B(311, TZ300, AIR_MOLEC(L), &
            1.8e-31_r8, 3.4_r8, 1.5e-11_r8, 1.9_r8, 0.6_r8, 0)

       !// OP + O2 + M -> O3 + M
       !// JPL number: A1 (JPL06, 20080617)
       r_op_o2_m(L) = rate3B(312, TZ300, AIR_MOLEC(L), &
            6.0e-34_r8, 2.4_r8, 0._r8, 0._r8, 0.6_r8, 0)

       !// O2 + H + M -> HO2 + M
       !// JPL number: B1 (JPL10-6, 20110610)
       r_o2_h_m(L) = rate3B(313, TZ300, AIR_MOLEC(L), &
            4.4e-32_r8, 1.3_r8, 7.5e-11_r8, -0.2_r8, 0.6_r8, 0)

       !// BrO + NO2 + M -> BrONO2 + M
       !// JPL number: G2 (JPL06, 20080617)
       r_no2_bro_m(L) = rate3B(314, TZ300, AIR_MOLEC(L), &
            5.2e-31_r8, 3.2_r8, 6.9e-12_r8, 2.9_r8, 0.6_r8, 0)

    end do !// do L = LMTROP+1, LMAX

    !// Other T,p-dependent calculations
    !// For a reaction A + B <--> C, there is typically a reaction rate
    !// (ka) given for the A + B -> C reaction, along with an equilibrium
    !// constant K. The decomposition of C to A+B is e.g. thermal, and is
    !// found from K = ka/kb
    !// ---------------------------------------------------------------
    !// Loop from start of stratosphere up to uppermost chemistry level
    do L = LMTROP+1, LMAX
 
       !// Temperature
       THE = TEMP(L)
       ZTEM = 1._r8 / THE

       !// HO2 + HO2 -> H2O2 + O2 (not standard 3-body)
       !// HO2 + HO2 --> H2O2 (+ O2)
       !// JPL number: bi-molecular B13 (JPL06, 20080617)
       R0 = 3.5e-13_r8 * exp(430._r8 * ZTEM)
       !// HO2 + HO2 + M --> H2O2 (+ O2 + M)
       !// JPL number: bi-molecular B13 (JPL06, 20080617)
       R2 = 1.7e-33_r8 * exp(1000._r8 * ZTEM) * AIR_MOLEC(L)

       ! Uten H2O effekt (Derwent-reaksjonen) ::  R = R0 + R2
       r_ho2_ho2_tot(L) = (R0 + R2)


       !// OH + HNO3 -> H2O + NO3 (not 3-body)
       !// JPL number: C9 (JPL06, 20080617)
       C0 = 2.4e-14_r8 * exp(  460._r8 * ZTEM )
       C2 = 2.7e-17_r8 * exp( 2199._r8 * ZTEM )
       C3M = 6.5e-34_r8 * exp( 1335._r8 * ZTEM ) * AIR_MOLEC(L)
       r_oh_hno3(L) = C0 + C3M / (  1._r8 + C3M / C2  ) !JPL06, 20080617


       !// HO2NO2 -> HO2 + NO2 thermal decomp (kb).
       !//   in equlibirum with NO2 + HO2 <--> HO2NO2 (ka)
       !// JPL number: C5 & Equlibrium constant no 2 (JPL06, 20080617)
       !// C1020 is already multiplied by AIR_MOLEC.
       r_ho2no2_heat(L) = r_ho2_no2_m(L) &
                        / ( 2.1e-27_r8 * exp(10900._r8 * ZTEM) )

       !// N2O5 -> NO2 + NO3 thermal decomp. (kb)
       !//   in equlibirum with NO2 + NO3 <--> N2O5 (ka)
       !// JPL number: C6 & Equlibrium constant no 5 (JPL06, 20080617)
       !// C1011X is already multiplied by AIR_MOLEC.
       r_n2o5_heat(L) = r_no2_no3_m(L) &
                      / ( 2.7e-27_r8 * exp(11000._r8 * ZTEM) )

       !// ClOO -> Cl + O2 thermal decomp. (ka)
       !//   in equlibrium with Cl + O2 <--> ClOO (kb)
       !// JPL number: F1 & Equlibrium constant no 11 (JPL06, 20080617)
       !// r_cl_o2_m is already multiplied by AIR_MOLEC.
       r_cloo_heat(L) = r_cl_o2_m(L) &
                      / ( 6.6e-25_r8 * exp(2502._r8 * ZTEM) )

       !// Cl2O2 + M  -> ClO + ClO + M (ka)
       !//   in equlibrium with ClO + ClO <--> Cl2O2 (kb)
       !// JPL number: F10 & Equlibrium constant no 14 (JPL06, 20080617)
       !// C36M is already multiplied by AIR_MOLEC.
       r_cl2o2_m(L) = r_clo_clo_m(L) &
                    / ( 9.3e-28_r8 * exp(8835._r8 * ZTEM) )


       !// NO + HO2 -> HNO3 (LeBras)
       !// 760Torr=1013.25mb, so we should multiply with 760./1013.25
       !// Their value is originally in percent!
       !// JPL2011 (nr 17) does not recommend this until other studies
       !// confirm the reaction, so we leave it out of the standard model runs.
       !//R0 = 1.e-2_r8 * (530._r8/TEMP(L) &
       !//     + 6.4e-4_r8*760._r8/1013.25_r8*PMIDL(L) &
       !//     - 1.73_r8)
       !//if (R0 .lt. 0._r8) then
       !//   R0 = 0._r8
       !//else if (R0 .gt. 1._r8) then
       !//   R0 = 0.999_r8
       !//end if
       !//r_no_ho2_b(L) = r_no_ho2(NINT(TEMP(L))-MINTEMP) * R0
       r_no_ho2_b(L) = 0._r8

    end do !// do L = LMTROP+1, LMAX


    !// --------------------------------------------------------------------
  end Subroutine TCRATE_TP_IJ_STR
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine TCRATE_onAER(ILOC,JLOC, TEM, CIWC, CLWC, CFR, DV, &
       PRES, AIR, AREAXY, LMT, PLAND, ZC_LOCAL, &
       RR_N2O5_H2O_AER, RR_HO2_AER, RR_QAER, RR_NO2_SOOT)
    !// --------------------------------------------------------------------
    !// Some species react on aerosol surfaces. E.g. N2O5 reacts with
    !// H2O on aerosols to produce HNO3. Also HO2 can react on aerosol
    !// surfaces.
    !// This routine provides calculations of the reaction rate
    !// for different aerosols; dust, ice, cloud droplets, etc.
    !//
    !// Reactions inside droplets is *NOT* taken into account.
    !//
    !// Ole Amund Sovde, May 2014
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, TRACER_ID_MAX, LDUST, LSALT, LBCOC, &
         LSULPHUR, LSOA, NPAR_DUST, NPAR_SALT
    use cmn_ctm, only: JMON
    use cmn_met, only: ZOFLE
    use cmn_parameters, only: CPI, R_UNIV
    use cmn_oslo, only: dustbinsradii, PR42HET
    use dust_oslo, only: dust_trsp_idx
    use seasalt, only: saltbinsradii, salt_trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ILOC, JLOC, LMT
    real(r8),dimension(LPAR), intent(in) :: &
         TEM, CLWC, CIWC, CFR, DV, PRES, AIR
    real(r8),dimension(TRACER_ID_MAX,LPAR), intent(in) :: ZC_LOCAL
    real(r8), intent(in) :: AREAXY, PLAND
    !// Output
    real(r8),dimension(LPAR), intent(out) :: &
         RR_N2O5_H2O_AER, &
         RR_HO2_AER, &
         RR_QAER, &
         RR_NO2_SOOT

    !// Locals
    real(r8) :: temp, dg, radius, iwc, lwc, liq_mass, &
         nr_cloud, nrconc_cloud, &
         sad_ice, sad_liq, reff_cloud, &
         cloud_drop_radius, cloud_drop_volume, cloud_drop_mass, r_ice, &
         kt1, kt2, &
         kt_n2o5_liq, kt_n2o5_ice, &
         kt_ho2_liq, kt_ho2_ice

    !// N2O5 + H2O + aerosol -> 2HNO3
    real(r8) :: kt_n2o5_dust, kt_n2o5_salt, &
              kt_n2o5_cloud_ice, kt_n2o5_cloud_liq, &
              kt_n2o5_soot, kt_n2o5_so4

    !// HO2 + aerosol -> scavenged or incorporated into aerosol
    real(r8) :: kt_ho2_dust, kt_ho2_salt, &
              kt_ho2_cloud_ice, kt_ho2_cloud_liq, &
              kt_ho2_soot, kt_ho2_oc, kt_ho2_so4, kt_ho2_soa

    !// NO2 + soot -> scavenged or incorporated into aerosol
    real(r8) :: kt_no2_soot

    !// particle calculations
    real(r8) :: aer_mass, aer_vol, aer_area, nr_aer, &
              particle_vol, &
              kt_qaer, so4wt

    real(r8) :: ZMID(LPAR)

    integer :: L, N, CID

    !// Uptake coefficients that can vary (fixed gammas are given below
    !// as parameters)
    real(r8) :: gamma_n2o5_salt, &
              gamma_n2o5_soot, &
              gamma_ho2_salt, &
              gamma_n2o5_so4, &
              gamma_ho2_so4

    !// Parameters
    !// Uptake coefficients on water ice and liquid droplets
    !// (Evans & Jacob, GRL, 32, L09813, doi:10.1029/2005GL022469, 2005)
    real(r8), parameter :: &
         gamma_n2o5_i = 0.02_r8, &
         gamma_n2o5_l = 0.02_r8
    !// Uptake on mineral dust
    !// Laboratory investigations have shown that N2O5 proceeds moderately
    !// efficiently on dust particles (N2O5) = 0.013 +/- 0.002
    !// (Wagner et al., Atmos. Chem. Phys., 8, 91-109, 2008)
    !// (RH=0%: 0.013, RH=28%:0.008, RH=59%:0.005)
    real(r8), parameter :: &
         gamma_n2o5_dust = 0.013_r8
    !// Uptake coefficients from JPL 2010:
    !// HO2 + h2o_liq -> gamma_ho2_liq = 0.020_r8 (mass accomodation)
    !// HO2 + h2o_ice -> gamma_ho2_ice = 0.025_r8 (gamma)
    !// HO2 + soot: gamma<1.d-2 on dry soot
    !// HO2 + h2so4 -> see gamma_ho2_so4
    real(r8), parameter :: &
         gamma_ho2_i = 0.025_r8, &
         gamma_ho2_l = 0.02_r8, &
         gamma_ho2_soot = 0.01_r8
    !// Matthews et al., ACPD 2014
    real(r8), parameter :: &
         gamma_ho2_dust = 0.018_r8
    !// Uptake on SOA and OC
    real(r8), parameter :: &
         gamma_ho2_soa = 1.e-4_r8

    !// NO2 + soot: gamma <1.d-4 on aged soot.
    !// Use 1.d-5 and assume also for dry.
    real(r8), parameter :: &
         gamma_no2_soot = 1.e-5_r8


    !// Uptake ala QAER at continental surface (for RR_QAER)
    !// JPL 2010 for CH3COO2 1.d-2, but that is maybe too large.
    !// JPL 2010 for CH3O2 ~4.d-3 (mass accomodation)
    real(r8), parameter :: &
         gamma_qaer = 4.e-3_r8

    !// Mean molecular speed [cm/s]
    real(r8), parameter :: vmms = 2.4e4_r8
    !// --------------------------------------------------------------------

    !// Initialise
    RR_N2O5_H2O_AER(:) = 0._r8
    RR_HO2_AER(:)      = 0._r8
    RR_QAER(:)         = 0._r8
    RR_NO2_SOOT(:)     = 0._r8

    !// If DUST is not included, these uptake rates will be too low.
    !// In that case we should either use the old method, or we should
    !// apply a climatology created with the full aerosol packages.
    !//
    !// Will apply old values for now; we need to generate a climatology
    !// to use when all aerosols are not included.
    !if (.not. LDUST) then
    if (.true.) then
       RR_N2O5_H2O_AER(:) = PR42HET(:,JLOC)
       !// PRQAER (simple removal by aerosol, as in CTM2)
       ZMID(:) = 0.5_r8 * (ZOFLE(1:LPAR,ILOC,JLOC) &
                          + ZOFLE(2:LPAR+1,ILOC,JLOC))
       ZMID(:) = ZMID(:) - ZMID(1) !// Height above ground
       call getRQAER(LPAR,ZMID,PLAND,RR_QAER)
       RR_HO2_AER(:) = RR_QAER(:)
       return
    end if

    !// UCI SAD
    if (.false.) then
       !// Need reference for this equation
       !RR_QAER(:) = trop_SAD(:,ILOC,JLOC,JMON) * 675._r8 * sqrt(TEM(:))
       !RR_N2O5_H2O_AER(:) = 0.05_r8 * RR_QAER(:)
       !RR_HO2_AER(:) = 0.01_r8 * RR_QAER(:)
       !RR_QAER(:) = 0.01_r8 * RR_QAER(:)

       !// Uptake rate is gamma * cg [cm/s] * A [cm2/cm3] * 1/4
       !//   cg = mean molecular speed
       !//   cg = 100 [cm/m] * sqrt(3RT/M) [m/s]  (RMS speed)
       !//  (cg = 100 [cm/m] * sqrt(2RT/M) [m/s]  (most probable speed))
       !// where M is molecular mass [kg/mol], so that rate unit is 1/s.
       !//   The factor 1/4 is the fraction of sweep-out area to the
       !//   surface area of a spherical particle.
       !// Split cg in: sqrt(T) * sqrt(3R) / sqrt(M)
       RR_QAER(:) = trop_SAD(:,ILOC,JLOC,JMON) &
                    * sqrt(TEM(:)) * sqrt(3._r8 * R_UNIV) * 25._r8
       !// N2O5 (0.0108kg/mol)
       RR_N2O5_H2O_AER(:) = 0.025_r8 * RR_QAER(:) / sqrt(0.0108_r8)
       !// HO2 (0.0033kg/mol)
       RR_HO2_AER(:) = 0.02_r8 * RR_QAER(:) / sqrt(0.0033_r8)
       !// Use CH3O2 for all (0.0047kg/mol)
       RR_QAER(:) = 0.004_r8 * RR_QAER(:) / sqrt(0.0047_r8)
       return
    end if
    


    !// About the method
    !// --------------------------------------------------------------------
    !// Reaction rate largely follows Dentener & Crutzen, JGR, 98(4),
    !// 7149-7163, doi:10.1029/92JD02979, 1993:
    !//   kt = 1 / (r/Dg + 4/(v*gamma)) * A(cm2/cm3)
    !// so that the loss of N2O5 is
    !//   d(N2O5)/dt = kt * [N2O5]
    !// Note that [H2O] is implicitly included in kt.
    !//
    !// In this routine we find kt = sum(kt_i) for all aerosol types i.


    !// General assumption
    !// --------------------------------------------------------------------
    !// It is assumed that the adsorbed gas stays with the aerosol until
    !// the aerosol is removed. In real life, there is an equlibrium
    !// between uptake and evaporation, and also there may be products
    !// that will evaporate. This is assumed for N2O5+H2O, where HNO3
    !// is produced will evaporate instantly.
    !// For HO2, there may be chemical conversion to H2O2, which could
    !// possibly evaporate (at least partly). However, this is not taken
    !// into account.


    !// Other gases which are removed in Oslo chemistry, using the
    !// old parameterisation QAER:
    !// CH3CH(O2)CH2OH (28)
    !// CH3COO2 (37) ~1.e-3_r8! JPL: prod may be HO2 and CH3COOH
    !// PANX (5)
    !// SEC-C6H13O2 (25)
    !// SEC-C4H9O2 (24)
    !// CH3COCH(O2)CH3 (27)
    !// C2H5O2 (23)
    !// C3H7O2 (49)
    !// CH3COD (CH3COCH2(O2)) (51)
    !// aromatic products 1 (29)
    !// aromatic products 3 (31)
    !// Peroxy radicals (RO2) (32)
    !// Peroxy radicals from ISOK (34)
    !// CH3O2 (22)

    !// Other possibilities
    !// H2O2 + h2o_liq: gamma_h2o2_liq = 0.18
    !// OH + h2o_liq: gamma = 0.1 or greater
    !// OH + h2o_ice: gamma = 0.1 or greater




    !// Loop through troposphere
    do L = 1, LMT

       !// May eventually use temperature
       TEMP = TEM(L)

       !// Simple approximation for diffusion coefficient [cm2/s]
       dg = 0.1_r8 * 1.e3_r8 / PRES(L)

       !// Initialise QAER uptake rate
       kt_qaer = 0._r8

       !// MINERAL DUST
       !// =================================================================
       kt_n2o5_dust = 0._r8
       kt_ho2_dust  = 0._r8
       !// Find SAD for dust
       !// Loop through size bins and calculate total area of dust, then
       !// divide by air volume.
       !sad_dust = 0._r8
       if (LDUST) then
          !// dustbinsradii(NPAR_DUST) from globalvariables [m]
          do N = 1, NPAR_DUST
             !// Get dust mass [kg]
             CID = dust_trsp_idx(N)
             aer_mass = ZC_LOCAL(CID,L)
             !// Total volume of dust particles, density 2600kg/m3
             aer_vol = aer_mass / 2.6e3_r8    !// [m3]
             !// Volume of one dust particle [m]
             particle_vol = 4._r8/3._r8 * CPI  * dustbinsradii(N)**3
             !// Number of dust particles of this size
             nr_aer = aer_vol / particle_vol
             !// Total dust area for this size bin [m2]
             aer_area = 3._r8 * CPI * dustbinsradii(N)**2 * nr_aer
             !// Divide by air volume and convert to [cm2/cm3]
             aer_area = aer_area / DV(L) * 1.e-2_r8
             !// Add up area [cm2/cm3]
             !sad_dust = sad_dust + aer_area
             !// Reaction constant k = A / (r/Dg + 4/(v*gamma))
             !// N2O5 + H2O + dust -> 2HNO3
             kt1 = aer_area / &
                  ((dustbinsradii(N)/dg) + (4._r8/(vmms * gamma_n2o5_dust)))
             !// Add up kt_n2o5_dust
             kt_n2o5_dust = kt_n2o5_dust + kt1

             !// HO2 + dust -> dust
             kt1 = aer_area / &
                  ((dustbinsradii(N)/dg) + (4._r8/(vmms * gamma_ho2_dust)))
             kt_ho2_dust = kt_ho2_dust + kt1

             !// qaer + dust ->
             kt1 = aer_area / &
                  ((dustbinsradii(N)/dg) + (4._r8/(vmms * gamma_qaer)))
             !kt_qaer = kt_qaer + kt1 * (PLAND + 0.1_r8*(1._r8 - PLAND))
             kt_qaer = kt_qaer + kt1

          end do !// do N = 1, NPAR_DUST
       end if !// if (LDUST) then


       !// SEA SALT
       !// =================================================================
       !// N2O5 + sea salt: Steward & Cox found 0.025, Hoffman et al
       !// found 0.038 which may be an upper limit, according to JPL 2010.
       gamma_n2o5_salt = 0.025_r8
       kt_n2o5_salt = 0._r8
       !// HO2 + sea salt: JPL 2010
       !// gamma_ho2_salt = 1.2d-2 (may use 5.7d-5exp(1560/T))
       gamma_ho2_salt = 0.012_r8
       kt_ho2_salt = 0._r8
       if (LSALT) then
          do N = 1, NPAR_SALT
             !// Get salt mass
             CID = salt_trsp_idx(N)
             aer_mass = ZC_LOCAL(CID,L)
             !// Total volume of salt particles, density 2200kg/m3 (dry)
             aer_vol = aer_mass / 2.2e3_r8    !// [m3]
             !// Volume of one salt particle [m]
             particle_vol = 4._r8/3._r8 * CPI  * saltbinsradii(N)**3
             !// Number of salt particles of this size
             nr_aer = aer_vol / particle_vol
             !// Total salt area for this size bin [m2]
             aer_area = 3._r8 * CPI * saltbinsradii(N)**2 * nr_aer
             !// Divide by air volume and convert to [cm2/cm3]
             aer_area = aer_area / DV(L) * 1.e-2_r8
             !// Reaction constant k = A / (r/Dg + 4/(v*gamma))

             !// N2O5 + H2O + salt -> 2HNO3
             kt1 = aer_area / &
                  ((saltbinsradii(N)/dg) + (4._r8/(vmms * gamma_n2o5_salt)))
             !// Add up kt_n2o5_salt
             kt_n2o5_salt = kt_n2o5_salt + kt1

             !// HO2 + sea salt
             kt2 = aer_area / &
                  ((saltbinsradii(N)/dg) + (4._r8/(vmms * gamma_ho2_salt)))
             kt_ho2_salt = kt_ho2_salt + kt2

             !// qaer + sea salt
             kt2 = aer_area / &
                  ((saltbinsradii(N)/dg) + (4._r8/(vmms * gamma_qaer)))
             kt_qaer = kt_qaer + kt2
          end do !// do N = 1, NPAR_SALT
       end if !// if (LSALT) then



       !// SOOT (BLACK CARBON)
       !// =================================================================
       !// N2O5 + H2O + soot: gamma < 4.d-5 on dry soot, 2.d-4 at RH50%
       !// N2O5 + soot -> NO2 + NO3: gamma < 4.d-6
       gamma_n2o5_soot = 1.e-5_r8
       kt_n2o5_soot = 0._r8
       kt_no2_soot = 0._r8
       kt_ho2_soot = 0._r8
       !// Assume hydrophobic to be dry, hydrophilic to be non-dry?
       if (LBCOC) then
          !// bcBB1fob:bcBB1fil 240-241, r=0.08d-6, rho=1.d3
          radius = 0.08e-6_r8 !// [m]

          !// Dry/new soot
          aer_mass = ZC_LOCAL(240,L) !// dry soot
          aer_vol = aer_mass / 1.e3_r8                   !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3  !// [m]
          nr_aer = aer_vol / particle_vol             !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer  !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8         !// [cm2/cm3]

          !// N2O5 + H2O + soot -> 2HNO3 : Not for dry soot
          !// HO2 + soot : Not for dry soot
          !// qaer + soot : Not for dry soot
          !// NO2 + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_no2_soot)))
          kt_no2_soot = kt_no2_soot + kt1


          !// Aged soot
          aer_mass = ZC_LOCAL(241,L)
          aer_vol = aer_mass / 1.e3_r8                   !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3  !// [m]
          nr_aer = aer_vol / particle_vol             !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer  !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8         !// [cm2/cm3]

          !// N2O5 + H2O + soot -> 2HNO3
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_n2o5_soot)))
          kt_n2o5_soot = kt_n2o5_soot + kt1
          !// HO2 + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_ho2_soot)))
          kt_ho2_soot = kt_ho2_soot + kt1
          !// qaer + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt1

          !// NO2 + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_no2_soot)))
          kt_no2_soot = kt_no2_soot + kt1


          !// bcFF1fob:bcFF1fil 242:243, r=0.039d-6, rho=1.d3
          !// bcBF1fob:bcBF1fil 244:245, r=0.039d-6, rho=1.d3
          radius = 0.039e-6_r8 !// [m]

          !// Dry/new soot
          aer_mass = ZC_LOCAL(242,L) + ZC_LOCAL(244,L) !// Dry soot
          aer_vol = aer_mass / 1.e3_r8                    !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3   !// [m]
          nr_aer = aer_vol / particle_vol            !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8        !// [cm2/cm3]

          !// N2O5 + H2O + soot -> 2HNO3 : Not for dry soot
          !// HO2 + soot : Not for dry soot
          !// qaer + soot : Not for dry soot
          !// NO2 + soot
          !// NO2 + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_no2_soot)))
          kt_no2_soot = kt_no2_soot + kt1

          !// Aged soot
          aer_mass = ZC_LOCAL(243,L) + ZC_LOCAL(245,L)
          aer_vol = aer_mass / 1.e3_r8                    !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3   !// [m]
          nr_aer = aer_vol / particle_vol            !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8        !// [cm2/cm3]

          !// N2O5 + H2O + soot -> 2HNO3
          kt2 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_n2o5_soot)))
          kt_n2o5_soot = kt_n2o5_soot + kt2
          !// HO2 + soot
          kt2 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_ho2_soot)))
          kt_ho2_soot = kt_ho2_soot + kt2
          !// qaer + soot
          kt2 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt2
          !// NO2 + soot
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_no2_soot)))
          kt_no2_soot = kt_no2_soot + kt1

       end if !// if (LBCOC) then



       !// ORGANIC CARBON
       !// =================================================================
       kt_ho2_oc = 0._r8
       if (LBCOC) then
          !// Only uptake on wet organic carbon

          !// ocBB1fil 231, r=0.08d-6, rho=1.35d3
          !// ocOCNfil 237, r=0.08d-6, rho=1.35d3
          radius = 0.08e-6_r8 !// [m]

          aer_mass = ZC_LOCAL(231,L) + ZC_LOCAL(237,L)
          aer_vol = aer_mass / 1.35e3_r8                 !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3  !// [m]
          nr_aer = aer_vol / particle_vol             !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer  !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8         !// [cm2/cm3]

          !// HO2 + OC (use gamma_ho2_soa)
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_ho2_soa)))
          kt_ho2_oc = kt_ho2_oc + kt1
          !// qaer + OC
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt1


          !// ocFF1fil: 233, r=0.039d-6, rho=1.d3
          !// ocBB1fil: 235, r=0.039d-6, rho=1.d3
          radius = 0.039e-6_r8 !// [m]

          !// Aged OC
          aer_mass = ZC_LOCAL(233,L) + ZC_LOCAL(235,L)
          aer_vol = aer_mass / 1.35e3_r8                 !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3  !// [m]
          nr_aer = aer_vol / particle_vol             !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer  !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8         !// [cm2/cm3]

          !// HO2 + OC (use gamma_ho2_soa)
          kt2 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_ho2_soa)))
          kt_ho2_oc = kt_ho2_oc + kt2
          !// qaer + soot
          kt2 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt2

       end if !// if (LBCOC) then


       !// SOA
       !// =================================================================
       kt_ho2_soa = 0._r8
       if (LSOA) then
          !// Assume uniform radius
          radius = 0.08e-6_r8 !// [m]
          !// Sum up mass of all SOA
          aer_mass = ZC_LOCAL(165,L)
          do N = 166,179
             aer_mass = aer_mass + ZC_LOCAL(N,L)
          end do
          aer_mass = aer_mass + ZC_LOCAL(182,L) + ZC_LOCAL(183,L) &
                     + ZC_LOCAL(188,L) + ZC_LOCAL(189,L) &
                     + ZC_LOCAL(190,L) + ZC_LOCAL(191,L)
          aer_vol = aer_mass / 1.35e3_r8                  !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3   !// [m]
          nr_aer = aer_vol / particle_vol              !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer   !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8          !// [cm2/cm3]

          !// HO2 + SOA
          kt1 = aer_area / &
              ((radius/dg) + (4._r8/(vmms * gamma_ho2_soa)))
          kt_ho2_soa = kt_ho2_soa + kt1
          !// qaer + SOA
          kt2 = aer_area / &
              ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt2

       end if !// if (LSOA) then



       !// SULPHATE
       !// =================================================================
       !// Assume 30% weight percent (not very sensitive to this)
       kt_n2o5_so4 = 0._r8
       kt_ho2_so4 = 0._r8
       if (LSULPHUR) then

          so4wt = 0.3_r8
          gamma_n2o5_so4 = exp( (-25.5265_r8 - 0.133188_r8*so4wt &
               + 0.00930846_r8*so4wt**2 - 0.90194e-5_r8*so4wt**3) &
               + (9283.76_r8 + 115.345_r8*so4wt - 5.19258_r8*so4wt**2 &
                  + 0.0483464_r8*so4wt**3)/TEMP &
               + (-851801._r8-22191.2_r8*so4wt + 766.916_r8*so4wt**2 &
                  - 6.85427_r8*so4wt**3)/(TEMP*TEMP) )
          gamma_ho2_so4 = 0.05_r8 !// 55%wt, T223, gamma>0.05

          radius = 0.08e-6_r8 !// [m]
          aer_mass = ZC_LOCAL(73,L)
          aer_vol = aer_mass / 1.e3_r8                    !// [m3]
          particle_vol = 4._r8/3._r8 * CPI * radius**3   !// [m]
          nr_aer = aer_vol / particle_vol              !// [#]
          aer_area = 3._r8 * CPI * radius**2 * nr_aer   !// [m2]
          aer_area = aer_area / DV(L) * 1.e-2_r8          !// [cm2/cm3]

          !// N2O5 + H2O + H2SO4xnH2O(l) -> 2HNO3
          !// Usually, HNO3 sticks to the sulphate aerosol, unless high
          !// sulphuric content. Let the nitrate package take care of
          !// this for now.
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_n2o5_so4)))
          kt_n2o5_so4 = kt_n2o5_so4 + kt1
          !// HO2 + H2SO2xH2O(l)
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_ho2_so4)))
          kt_ho2_so4 = kt_ho2_so4 + kt1
          !// qaer + H2SO2xH2O(l)
          kt1 = aer_area / &
               ((radius/dg) + (4._r8/(vmms * gamma_qaer)))
          kt_qaer = kt_qaer + kt1

       end if !// if (LSULPHUR) then



       !// OTHER AEROSOLS?
       !// =================================================================
       !// Nitrate? Probably covered by sea salt.
       !// SOA???



       !// CLOUD ICE AND WATER
       !// =================================================================

       !// SAD for cloud droplets (i.e. liquid)
       if (CLWC(L) .gt. 1.e-10_r8 .and. CFR(L) .gt. 0.01_r8) then
          !// First find droplet size and mass, then calculate SAD

          !// Mass of liquid [kg]
          liq_mass = CLWC(L) * AIR(L)

          !// Find LWC [g/m3]. Take cloud fraction into account.
          lwc = liq_mass * 1.e3_r8 / (CFR(L) * DV(L))   !// [g/m3]

          !// Effective cloud particle radius according to Fouquart el al.,
          !// Rev.Geophys., 28, 145-166, doi:10.1029/RG028i002p00145, 1990.
          !//   reff = 11*w + 4
          !// where w is liquid water content [g/m3]
          !// Another option is Martin et al., 1994.
          reff_cloud = 4._r8 + 11._r8 * lwc   !// [um]

          !// Typically the optical reff is somewhat larger than the physical
          !// r by 1-2um. Therefore downsize reff by 2.0um for droplets 9-13um
          !// and 1.0um for those between 6-9um.
          !// http://www-das.uwyo.edu/~geerts/cwx/notes/chap08/moist_cloud.html
          if (reff_cloud .ge. 9._r8) then
             cloud_drop_radius = (reff_cloud - 2._r8) * 1.e-4_r8   !// [cm]
          else if (reff_cloud .ge. 6._r8 .and. reff_cloud .lt. 9._r8) then
             cloud_drop_radius = (reff_cloud - 1._r8) * 1.e-4_r8   !// [cm]
          else if (reff_cloud .ge. 0.25_r8) then
             cloud_drop_radius = (reff_cloud - 0.25_r8) * 1.e-4_r8 !// [cm]
          else
             cloud_drop_radius = reff_cloud * 1.e-4_r8            !// [cm]
          end if

          !// Constrain droplet radius?
          !// LWC may also contain rain water, and we have no possibility to
          !// separate that from cloud water. A large radius will
          !// reduce the sad_liq, but sad_liq will still be large enough
          !// for our purpose. More suitable will probably be to constrain
          !// cloud droplet number (nr_cloud), and set a suitable droplet
          !// radius. See calculation of nr_cloud below.

          !// Droplet volume [cm3] and mass [g]
          cloud_drop_volume = (4._r8/3._r8) * CPI * cloud_drop_radius**3 
          cloud_drop_mass = cloud_drop_volume   !// Assuming rho_l = 1g/cm3

          if (cloud_drop_mass .gt. 0._r8) then
             !// Number of cloud droplets [#] (liq_mass is kg)
             nr_cloud = liq_mass * 1.e3_r8 / cloud_drop_mass
             !// Number concentration (#/cm3)
             nrconc_cloud = nr_cloud / (CFR(L)*DV(L)*1.e6_r8)

             !// Possible constraint on nr_cloud and cloud_drop_radius:
             !//   Marine air: < 100#/cm3
             !//   Continental air: < 900#/cm3
             !// If used, then also constrain cloud_drop_radius
             !//   Marine air: 15.d-4 cm
             !//   Continental air: 8.d-4 cm
             !// However, such a constraint will usually increase the total
             !// sad_liq, often by an order of magnitude. This should probably
             !// not be necessary for our purposes.

             !// SAD of cloud droplets [cm2/cm3]
             sad_liq = 4._r8 * CPI * cloud_drop_radius**2 &
                               * nr_cloud / (DV(L) * 1.e6_r8)
             !// Need to multiply by cloud fraction, taking into account only
             !// the grid box part in the cloud.
             sad_liq = sad_liq * CFR(L)
          else
             !// No droplets, no area
             sad_liq     = 0._r8
             nr_cloud    = 0._r8
          end if !// if (cloud_drop_mass .gt. 0._r8) then

       else
          !// No droplets, no area
          reff_cloud  = 0._r8
          cloud_drop_radius = 0._r8
          cloud_drop_volume = 0._r8
          cloud_drop_mass   = 0._r8
          sad_liq     = 0._r8
          nr_cloud    = 0._r8
       end if


       !// SAD for cloud ice
       if (CIWC(L) .gt. 1.e-10_r8) then
          !// Find IWC [g/m3], need conversion from kg/kg
          !// Take cloud fraction into account
          iwc = CIWC(L) * 1.e3_r8 * AIR(L) / (CFR(L) * DV(L))
          !// SAD according to Heymsfield & McFarquhar, J. Atmos. Sci.,
          !// vol 53(17), 2424-2451,
          !// doi:10.1175/1520-0469(1996)053<2424:HAOCIT>2.0.CO;2, 1996.
          sad_ice = 1.e-4_r8 * iwc**0.9_r8     !// [cm2/cm3]
          !// Calculate the r_eff using the relationship of Fu, J.Clim.,
          !// vol 9(9), 2058-2057,
          !// doi: 10.1175/1520-0442(1996)009<2058:AAPOTS>2.0.CO;2, 1996,
          !// i.e. r_eff = sqrt(3)/(3*rho_i) * IWC/Ac
          !// Unit of rho_i is g/m3, but we find r_eff in cm by multiplying
          !// iwc by 1.d-6m3/cm3:
          r_ice = (1.73205_r8 / (3._r8 * 0.917_r8)) * ((iwc*1.e-6_r8)/sad_ice)
          !// The value adopted in von Kuhlmann and Lawrence is too low,
          !// according to Schmitt and Heymsfield, J.Appl.Meteo.,
          !// vol. 44, 467-474, doi:10.1175/JAM2209.1, 2005.
          sad_ice = 10._r8 * sad_ice
          !// Need to multiply by cloud fraction, taking into account only
          !// the grid box part in the cloud.
          sad_ice = sad_ice * CFR(L)
       else
          !// No ice particles
          sad_ice = 0._r8
          r_ice = 0._r8
       end if



       !// Calculate reaction coefficient
       !// SAD for ice particles in now linked to the IWC by the
       !// parameterization of Heymsfield and McFarquar (1996)
       !// and the effective radii from Fu (1996)
       !//
       !// Return to the original formulation in dentener and crutzen
       !// of : kt = 1/ (r/Dg + 4/(v*gamma)) * A(cm2/cm3)
       if (cloud_drop_radius .gt. 0._r8) then
          !// N2O5 + H2O + liquid droplets -> 2HNO3
          kt_n2o5_liq = 1._r8 / &
               ((cloud_drop_radius/dg) + (4._r8/(vmms * gamma_n2o5_l)))
          !// HO2 + liquid droplets
          kt_ho2_liq  = 1._r8 / &
               ((cloud_drop_radius/dg) + (4._r8/(vmms * gamma_ho2_l)))
       else
          kt_n2o5_liq = 0._r8
          kt_ho2_liq  = 0._r8
       end if

!// REMOVE
       if (kt_n2o5_liq*sad_liq .gt. 1._r8) then
          print*,'tropchem_oslo_rates: large kt_n2o5',L, CLWC(L), &
               kt_n2o5_liq, sad_liq,cloud_drop_volume,lwc, &
               reff_cloud, cloud_drop_radius,nr_cloud
          stop
       end if

       if (r_ice .gt. 0._r8) then
          !// N2O5 + H2O + ice -> 2HNO3
          kt_n2o5_ice = 1._r8 / &
               ((r_ice/dg) + (4._r8/(vmms * gamma_n2o5_i)))
          !// HO2 + ice
          kt_ho2_ice  = 1._r8 / &
               ((r_ice/dg) + (4._r8/(vmms * gamma_ho2_i)))
       else
          kt_n2o5_ice = 0._r8
          kt_ho2_ice  = 0._r8
       end if

       !// Reaction rate for ice and cloud droplets
       if (sad_ice .gt. 0._r8) then
          !// N2O5 + H2O + ice
          kt_n2o5_cloud_ice = kt_n2o5_ice * sad_ice
          !// HO2 + ice
          kt_ho2_cloud_ice = kt_ho2_ice * sad_ice
          !// qaer + ice
          kt_qaer = kt_qaer + sad_ice / &
               ((r_ice/dg) + (4._r8/(vmms * gamma_qaer)))
       else
          !// No ice
          kt_n2o5_cloud_ice = 0._r8
          kt_ho2_cloud_ice = 0._r8
       end if

       if (sad_liq .gt. 0._r8) then
          !// N2O5 + H2O + liquid droplets
          kt_n2o5_cloud_liq = kt_n2o5_liq * sad_liq
          !// HO2 + liquid droplets
          kt_ho2_cloud_liq = kt_ho2_liq * sad_liq
          !// qaer + liquid droplets
          kt_qaer = kt_qaer + sad_liq / &
               ((cloud_drop_radius/dg) + (4._r8/(vmms * gamma_qaer)))
       else
          !// No liquid droplets
          kt_n2o5_cloud_liq = 0._r8
          kt_ho2_cloud_liq = 0._r8
       end if


       !// Calculate total loss rate
       RR_N2O5_H2O_AER(L) = kt_n2o5_dust + kt_n2o5_salt &
                          + kt_n2o5_soot &
                          + kt_n2o5_so4 &
                          + kt_n2o5_cloud_ice + kt_n2o5_cloud_liq

       RR_HO2_AER(L) = kt_ho2_dust + kt_ho2_salt &
                     + kt_ho2_soot + kt_ho2_oc &
                     + kt_ho2_so4 &
                     + kt_ho2_soa &
                     + kt_ho2_cloud_ice + kt_ho2_cloud_liq

       RR_QAER(L) = kt_qaer

       !// Will wait with this reaction: Not
       !RR_NO2_SOOT(L) = kt_no2_soot
    end do !// do L = 1, LMT


    !// May include eventually
    ! kn2o5aq and nh3so4 can be done implicitly, 
    ! it has occurred that these rates have
    ! become negative over antarctica (aug 1993), 
    ! therefore put minimum value of 0. (AJ jul1999)
    !
    !cmk    rr(i,j,kn2o5aq)=max(0.,het_n2o5(i,j)/
    !       1e-9*(y(jl,iso4)+y(jl,imsa))/aird
    !cmk multiplication moved to EBI
    !rr_n2o5aq = max(0._r8, het_n2o5(i,j))/1e-9/aird
    !
    ! knh3so4 is uptake coefficient on H2SO4. 
    ! 1 uptake of NH3 consumes 1 acid molecule.  
    ! 
    !rr(i,j,knh3so4)=max(0.,het_nh3(i,j))/1e-9/aird

    !// -------------------------------------------------------------------
  end subroutine TCRATE_onAER
  !// ---------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TCRATE_HET_IJ ( &
       LM,   TRACER_ID_MAX, LMTROP, TEMP, &
       ZC_LOCAL, &
       PARTAREA, PSC1, PSC2, &
       !// heterogeneous reaction rates
       spsADPARBK, spsACPARBK, spsBDPARBK, &
       spsBCPARBK, spsECPARBK, &
       spsFDPARBK, spsFCPARBK, spsGCPARBK, &
       spsADPAR, spsACPAR, spsBDPAR, &
       spsBCPAR, spsECPAR, &
       spsFDPAR, spsFCPAR, spsGCPAR, spsBHPAR )
    !// --------------------------------------------------------------------
    !// Subroutine assigns values to reaction rates (dependent on
    !// temperature and more).
    !// Temperature and density fields are sendt from calling routine.
    !//
    !// Including 3-body reactions (gas-phase) and heterogeneous reactions)
    !//
    !// October 2008: Updated for Oslo CTM3; IJ-block structure.
    !// March 2007: Updated for JPL06
    !// ..... 2006: New calculations of heterogeneous chemistry rates
    !//
    !// Name convention:
    !//      A: N2O5
    !//      B: ClONO2
    !//      C: HCL
    !//      D: H2O
    !//      E: HOCL
    !//      F: BRONO2
    !//      G: HOBr
    !//      H: HBr
    !// --------------------------------------------------------------------
    use psc_microphysics, only: LPSC, LAEROSOL
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input parameters
    integer, intent(in) :: LM, TRACER_ID_MAX, LMTROP
    real(r8), intent(in)  :: &
         TEMP(LM), &           !// Temperature
         PARTAREA(LM), &       !// SA surface density (background aerosol)
         PSC1(LM), &           !// PSC1 surface area density
         PSC2(LM), &           !// PSC2 surface area density
         ZC_LOCAL(TRACER_ID_MAX,LM) !// Tracer densities
    !// Output
    real(r8), intent(out) :: &
         !// Heterogenous rate constants (aerosols):
         spsADPARBK(LM), spsACPARBK(LM), spsBDPARBK(LM), &
         spsBCPARBK(LM), spsECPARBK(LM), &
         spsFDPARBK(LM), spsFCPARBK(LM), spsGCPARBK(LM), &
         !// Heterogenous rate constants (PSCs):
         spsADPAR(LM), spsACPAR(LM), spsBDPAR(LM), &
         spsBCPAR(LM), spsECPAR(LM), &
         spsFDPAR(LM), spsFCPAR(LM), spsGCPAR(LM), spsBHPAR(LM)

    !// Local variables
    real(r8) :: &
         GAMAD, GAMAC, GAMBD, GAMBC, GAMEC, GAMFD, GAMGC, &
         WN2O5, WCLONO2, WHOCL, WBRONO2, WHOBR, &
         PHCL,  PH2O,    aH2O,  g0,      HHCl, &
         Gs,    PPP,     Gclc,  adivl,   f, &
         gef,   z,       w,     HClstr,  x, &
         hd,    con_ndpp,  volHOCl,      r4, &
         M_HCl, &    ! HCl concentration
         THE, &      ! temperature
         H2O_MOLEC  !// H2O molecular density (molec/cm3):

    integer :: LMAX !// Uppermost chemistry level
    integer :: L !// Loop index

    !// Boltzmann's constant mbcm3/(K molecule)
    real(r8), parameter :: kB = 1.3806d-19  !// mbcm^3/(K molecule)
    !// --------------------------------------------------------------------

    !// For het. chem. calculations:
    !// A: N2O5, B: ClONO2, C: HCl, D: H2O, E:HOCl, F:BrONO2, G:OHBr

    !// Calculate up to layer LM or LM-1, must correspond to chemistry
    LMAX = LM

    !// Initialize PSC rates
    spsADPAR(:) = 0._r8
    spsACPAR(:) = 0._r8
    spsBDPAR(:) = 0._r8
    spsBCPAR(:) = 0._r8
    spsECPAR(:) = 0._r8
    spsFDPAR(:) = 0._r8
    spsFCPAR(:) = 0._r8
    spsGCPAR(:) = 0._r8
    spsBHPAR(:) = 0._r8

    if (LPSC) then

       !// Loop array length (could be longitude or height)
       do L = LMTROP+1, LMAX

          !// Only calculate when PSCs are present
          if ((PSC1(L).gt.0._r8).or.(PSC2(L).gt.0._r8)) then
             THE = TEMP(L)
             WN2O5    = 3636.73_r8 * sqrt(THE / 108._r8) ! N2O5
             WCLONO2  = 3636.73_r8 * sqrt(THE / 97.45_r8) ! CLONO2
             WHOCL    = 3636.73_r8 * sqrt(THE / 52.45_r8) ! HOCL
             WBRONO2  = 3636.73_r8 * sqrt(THE / 142._r8) ! BRONO2
             WHOBR    = 3636.73_r8 * sqrt(THE / 97._r8) ! HOBR

             !// THE RATES
             !// ------------------------------------------------

             !// N2O5+H2O -> HNO3+HNO3 
             !// For NAT depends on thickness - 0.0004 is average
             !// PSC1: JPL T5-2.26
             !// PSC2: almost JPL T5-2.24
             spsADPAR(L) = (0.0004_r8 * PSC1(L) &
                           + 0.025_r8 * PSC2(L)) * WN2O5

             !// N2O5+HCL -> CLONO+HNO3
             !// PSC1: JPL T5-2.31
             !// PSC2: JPL T5-2.30
             spsACPAR(L) = (0.003_r8 * PSC1(L) &
                           + 0.03_r8 * PSC2(L)) * WN2O5
  
             !// CLONO2+H2O -> HOCL+HNO3
             !// PSC1: JPL T5-2.72
             !// PSC2: JPL T5-2.70
             !// For NAT 0.004 is averaged value assuming no
             !// temperature dependence
             spsBDPAR(L) = (0.004_r8 * PSC1(L) &
                            + 0.3_r8 * PSC2(L)) * WCLONO2

             ! CLONO2+HCL -> CL2+HNO3
             !// PSC1: JPL T5-2.76
             !// PSC2: JPL T5-2.75
             spsBCPAR(L) = (0.23_r8 * PSC1(L) &
                            + 0.26_r8 * PSC2(L)) * WCLONO2

             ! HOCL+HCL -> CL2+H2O
             !// PSC1: JPL T5-2.61
             !// PSC2: JPL T5-2.61
             spsECPAR(L) = (0.14_r8 * PSC1(L) &
                            + 0.26_r8 * PSC2(L)) * WHOCL

             ! BRONO2+H2O -> HOBR+HNO3
             !// PSC1: no reference
             !// PSC2: JPL T5-2.97
             spsFDPAR(L) = (0.006_r8 * PSC1(L) &
                            + 0.26_r8 * PSC2(L)) * WBRONO2

             ! BRONO2+HCL -> BRCL+HNO3
             !// PSC1: Possibly ~JPL T5-2.99
             !// PSC2: Possibly ~JPL T5-2.99
             spsFCPAR(L) = (0.3_r8 * PSC1(L) &
                            + 0.5_r8 * PSC2(L)) * WBRONO2

             ! HOBR+HCL -> BRCL+H2O
             !// PSC1: JPL T5-2.89
             !// PSC2: JPL T5-2.89
             !// Should these be reversed?
             spsGCPAR(L) = (0.25_r8 * PSC1(L) &
                            + 0.1_r8 * PSC2(L)) * WHOBR

             ! CLONO2+HBR -> BRCL+HNO3   ! *** Added by SPS /12/03/07/
             !// PSC1: JPL T5-2.83
             !// PSC2: JPL T5-2.83
             spsBHPAR(L)= (0.3_r8 * PSC1(L) &
                           + 0.5_r8 * PSC2(L)) * WCLONO2

          end if !if ((PSC1(L).gt.0._r8).or.(PSC2(L).gt.0._r8)) then
       end do
    end if !if (LPSC) then


    !// Initialize aerosol rates
    spsADPARBK(:) = 0._r8
    spsACPARBK(:) = 0._r8
    spsBDPARBK(:) = 0._r8
    spsBCPARBK(:) = 0._r8
    spsECPARBK(:) = 0._r8
    spsFDPARBK(:) = 0._r8
    spsFCPARBK(:) = 0._r8
    spsGCPARBK(:) = 0._r8
    !// the heterogenous reactions
    !//     BRONO2+H2O   : FDPARBK
    !//     HOBr + HCl   : GCPARBK
    !//     ClONO2 + HCl : BCPARBK
    !//     HOCL + HCL   : ECPARBK
    !//     N2O5 + H2O   : ADPARBK
    !//     N2O5 + HCl   : ACPARBK
    !//     ClONO2 + H2O : BDPARBK


    if (LAEROSOL) then


       !// Loop from start of stratosphere up to uppermost chemistry level
       do L = LMTROP+1, LMAX 

          !// Only calculate when aerosols are present
          if (PARTAREA(L).gt.0._r8) then

             M_HCl = ZC_LOCAL(111,L) ! HCl [molec/cm3]
             H2O_MOLEC = ZC_LOCAL(114,L) ! H2O [molec/cm3]

             THE = TEMP(L)

             !// SPS velocities
             WN2O5  = 3636.73_r8 * SQRT(THE / 108._r8)
             WCLONO2= 3636.73_r8 * SQRT(THE / 97.45_r8)
             WHOCL  = 3636.73_r8 * SQRT(THE / 52.45_r8)
             WBRONO2= 3636.73_r8 * SQRT(THE / 142._r8)
             WHOBR  = 3636.73_r8 * SQRT(THE / 97._r8)


             !// N2O5 + H2O:
             !// Hanson, Ravishankara and Solomon, JGR, 99:D2, 3615-3629, 1994.
             GAMAD=0.1_r8

             !// N2O5 + HCl:
             GAMAC=0.0_r8

             !// BRONO2+H2O:
             !// From Hanson and Ravishankara, GRL 22:4, 385-388, 1995.
             GAMFD=0.4_r8

             if (PARTAREA(L) .gt. 1e-10_r8) then

                !// ClONO2 + HCl:
                !// From Hanson and Ravishankara, J. Phys. Chem 98:22,
                !// 5728-5735, 1994.
                !// With radius a=0.1e-6m=10^{-5}cm
                PHCL = M_HCl * kB * THE / 1013.25_r8 ! partial pressure (atm)
                PH2O = H2O_MOLEC * kB * THE          ! partial pressure (mb)
                ! mobility
                aH2O = pH2O / 10._r8**(9.217_r8 - 2190._r8 / (THE - 12.7_r8))
                ! g0
                g0   = 1.18e-4_r8 + 9.1e-3_r8 * aH2O + 0.5_r8 * aH2O * aH2O
                ! H*hcl (1/(Matm)
                HHCl = exp(6250._r8 / THE - 10.414_r8) * aH2O**3.49_r8
                Gs   = 576._r8 * aH2O * HHCl * pHCl ! Gamma_s
                PPP  = 2.e3_r8 * HHCl * pHCl / aH2O ! P
                Gclc = g0 * sqrt(1._r8 + PPP) ! Gamma_calc
                adivl= 0.1_r8 / (1.4e-2_r8 * sqrt(1._r8 / aH2O)) ! a/l
                f = 1._r8 / tanh(adivl) - 1._r8 / adivl          ! f
                ! gamma_e
                gef  = 1._r8 / (1._r8 / (Gs + f * Gclc) + 1._r8 / 0.3_r8)

                GAMBC = gef * (Gs + f * Gclc * PPP / (1._r8 + PPP)) &
                            / (Gs + f * Gclc) ! ClONO2+HCl
                GAMBD = gef - GAMBC ! ClONO2+H2O

                !// HOCl + HCl:
                !// Hanson, Ravishankara and Solomon, JGR. 99:D2,
                !// 3615-3629, 1994.
                z = log(pH2O)
                w = ((-14.458_r8 + 0.62456_r8 * z) * THE + 3565._r8) / &
                     (44.777_r8 + 1.3204_r8 * z - 0.19988_r8 * THE) ! w% H2SO4
                if (w .lt. 30._r8) w = 30._r8                ! restriction w
                HHCl   = 10._r8**(15.514_r8 - 0.1791_r8 * w) ! H*hcl (Table 3)
                HCLstr = HHCl * pHCl                         ! [HCl]*
                x = 60._r8 - w

                if ( x .ge. 0._r8) then
                   HD = 15._r8 + 3._r8 * x      ! H_HOCl sqrt(D_l) (Table 3)

                   !// HOCl+HCl(a) is the bulk reaction, since l=1.3mkm>raer,
                   !// so its rate is proportional to the volume of aerosol,
                   !// not the surface!

                   ! Using Eqn 2b and 4b from Hanson, Ravishankara 
                   !   and Solomon, JGR. 99:D2, 3615-3629, 1994.
                   ! Will find R4=g_e omega/4=HRTsqrt(k1Dl)f,
                   ! where f=coth(a/l)+l/a and express R (Latm/(K mol)) with kB:
                   ! R = kB[mbcm3/(K molecule)]*N_A[molecule/mol]
                   !     -----------------------------------------
                   !           1.d3 cm3/L 1013.25 mb/atm
                   !
                   !   = kB N_A/1000 Latm/(K mol)
                   !
                   ! Since l>a, we have f = coth(q)-1/q = a/(3sqrt(Dl/k1))
                   ! (according to Shi et al., JGR 106, D20, 24259-24274)
                   ! R4 = HRTsqrt(k1Dl) a/3 sqrt(k1)/sqrt(Dl)
                   !    = Hsqrt(D1) k1 a/3 1.d-3 kB N_A
                   !      -----------------------------
                   !            1013.25 sqrt(Dl)
                   ! k1 = k2[HCl]* = 1.d5 [HCl]*, Dl=1.d-8, a=1.d-5cm
                   ! Hsqrt(D1) = HD in this calculation (15+3(60-w))
                   ! R4 = HD 1.d5 [HCl]* 1.d-5/3 1.d-3 kB N_A
                   !      -----------------------------------
                   !                1013.25 1.d-4
                   !     = 10 HD/3 [HCl]* kB N_A /1013.25

                   ! convert. coef from number to partial pressure (atm)
                   con_ndpp = kB * THE / 1013.25_r8

                   volHOCl  = 10._r8 * con_ndpp * HClstr / 3._r8 * AVOGNR
                   R4       = HD * volHOCl
                   Gef      = R4 / WHOCL
                   GAMEC    = min(Gef, 1._r8)

                else

                   !// Do not allow gamma3=0. if wH2SO4>60% and fit it to
                   !// Fig.5 in D.R. Hanson, A.R. Ravishankara & S. Solomon,
                   !// JGR 99,3615,1994, but it seems unimportant for the
                   !// model results
                   GAMEC = 10._r8**(-0.408_r8 * w + 21.735_r8)

                end if

             else

                !// Ole Amund Sovde 20051111: Not sure where this comes from yet...
                GAMBC = 0.008_r8
                GAMBD = 0.010_r8
                GAMEC = 0.001_r8

             end if

             spsADPARBK(L) = GAMAD * PARTAREA(L) * WN2O5
             spsACPARBK(L) = GAMAC * PARTAREA(L) * WN2O5
             spsBDPARBK(L) = GAMBD * PARTAREA(L) * WCLONO2
             spsBCPARBK(L) = GAMBC * PARTAREA(L) * WCLONO2
             if (spsBCPARBK(L) .ne. spsBCPARBK(L)) then
                print*, 'THIS IS WRONG: ',spsBCPARBK(L),spsBDPARBK(L)
                print*, 'THE',THE
                print*, 'PHCL',PHCL
                print*, 'H2O',PH2O,aH2O
                print*, 'g0',g0
                print*, 'HHCL',HHCl
                print*, 'Gs',Gs
                print*, 'PPP',PPP
                print*, 'Gclc',Gclc
                print*, 'adivl',adivl
                print*, 'f',f
                print*, 'gef',gef
                print*, 'GAMBC',GAMBC
                print*, 'GAMBD',GAMBD
                print*, 'H2O_MOLEC',H2O_MOLEC,L
                print*, 'WCLONO2',WCLONO2
                print*, 'PARTAREA',PARTAREA(L)
             end if
             spsECPARBK(L) = GAMEC * PARTAREA(L) * WHOCL
             spsFDPARBK(L) = GAMFD * PARTAREA(L) * WBRONO2

             !// BrONO2+HCl:
             !// Hanson, Ravishankara and Lovejoy, JGR 101, D4, (9063-9070),
             !// 1996
             spsFCPARBK(L) = 0.8_r8 * PARTAREA(L) * WBRONO2

             !// HOBr + HCl:
             !// Hanson and Ravishankara, J. Phys. Chem 98:22, 5728-5735, 1994.
             GAMGC = 0.2_r8
             spsGCPARBK(L) = GAMGC * PARTAREA(L) * WHOBR

 
          end if !if (PARTAERA(I).gt.0._r8) then
       end do !// do L=LMTROP+1,LMAX

    end if ! If (LAEROSOL) Then

    !// End calculation of heterogeneous chemistry kinetics
 
    !// --------------------------------------------------------------------
  end Subroutine TCRATE_HET_IJ
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine inSAD(filename)
    !// --------------------------------------------------------------------
    !// Read Surface Area Density file SAD_12m_1x1.nc and map SAD to ctm
    !// grid. Based on UCI inSAD.
    !//
    !// Ole Amund Sovde, April 2016
    !// --------------------------------------------------------------------
    use cmn_parameters, only: A0, CPI, CPI180
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: PIJL, AREAXY, XDEDG, YDEDG
    use regridding, only: Regrid_Column_Weights, e_grid
    use netcdf
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_4d, &
         get_netcdf_att_char
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: filename

    !// Local variables - allocatable
    real(r8), allocatable, dimension(:) :: &
         inLon,  inLat,  inLevP, &
         inTime, inLonEdge,  inLatEdge, &
         inLevEdgeP, inPE, &
         in1D, inXYBOX
    real(r8), allocatable, dimension(:,:,:,:) :: &
         inProperty
    real(r8), allocatable, dimension(:,:,:) :: &
         outPropertyInLev
    real(r8), allocatable, dimension(:,:) :: R8XY, Weights

    !// Local variables - fixed resolution
    real(r8), dimension(IPAR,JPAR) :: EDXY

    integer :: &
         status, ncid, nLon, nLat, nLev, nTime, &
         I, J, L, LL, M

    real(r8) :: outLevEdgeP(LPAR+1), P1, P2, inP1, inP2, dp_in, dp_out
    real(r8) :: out1D(LPAR)

    character(len=25)  :: varname
    character(len=25)  :: units

    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'inSAD'
    !// --------------------------------------------------------------------

    if (maxval(PIJL) .eq. 0._r8) then
       write(6,'(a)') f90file//':'//subr//': PIJL is not set!'
       stop
    end if

    !// Initialise SAD for model run
    trop_SAD(:,:,:,:) = 0._r8

    !// Latitude, longitude, pressure (routine allocates)
    call get_netcdf_var_1d( filename, 'lon',  inLon  )
    call get_netcdf_var_1d( filename, 'lat',  inLat  )
    call get_netcdf_var_1d( filename, 'lev',  inLevP )
    call get_netcdf_var_1d( filename, 'nTime', inTime )
    call get_netcdf_var_1d( filename, 'ilon',  inLonEdge  )
    call get_netcdf_var_1d( filename, 'ilat',  inLatEdge  )
    call get_netcdf_var_1d( filename, 'ilev',  inLevEdgeP )
         
    !// Dimensions of inputs
    nLon  = size( inLon  )
    nLat  = size( inLat  )
    nLev  = size( inLevP )
    nTime = size( inTime )

    if (nTime .ne. 12) then
       write(6,'(a)') f90file//':'//subr// &
            ': File name is not a one year monthly file: '//trim(filename)
       stop ' EXECUTION STOPPED'
    end if


    !// Allocate array to be read from file
    varname = 'SAD_TROP'
    allocate( inProperty( nLon, nLat, nLev, nTime ) )
    inProperty(:,:,:,:) = 0._r8

    !// Get data for all months
    call get_netcdf_var_4d(filename,varname,inProperty, nlon,nlat,nlev,ntime)

    ! Check the units for varname
    call get_netcdf_att_char( filename, varname, 'units', units )

    write(6,'(4(A,X))') 'Units: ', &
                  trim(filename), trim(varname), trim(units)


    
    !// Grid box areas (all zonal boxes are of same size)
    allocate( inXYBOX(nLat) )
    inXYBOX(:) = 0._r8
    do J = 1, nLat
       inXYBOX(J) =  A0*A0 * CPI180 * (inLonEdge(2) - inLonEdge(1)) &
            * (sin(CPI180*inLatEdge(J+1)) - sin(CPI180*inLatEdge(J)))
    end do


    !// Horizontal regridding
    !// --------------------------------------------------------------------
    !// Array with horizontal resolution of output grid
    !// and vertical resolution of input grid
    allocate( outPropertyInLev(IPAR,JPAR,nLev ), R8XY(nLon,nLat), &
              in1D(nLev), Weights(nLev,LPAR), inPE(nLev+1) )

    ! Loop over time, altitudes
    do M = 1, nTime

       outPropertyInLev(:,:,:) = 0._r8

       do L = 1, nLev
          !// Multiply by area before interpolation
          do I = 1, nLon
             R8XY(I,:) = inProperty(I,:,L,M) * inXYBOX(:)
          end do
          !// Regrid emissions to CTM grid
          call e_grid(R8XY, inLonEdge, inLatEdge, nLon, nLat, &
               EDXY, XDEDG, YDEDG, IPAR, JPAR, 1)

          outPropertyInLev(:,:,L) = EDXY(:,:) / AREAXY(:,:)
       end do

       !// Vertical regridding for this month
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR+1
                outLevEdgeP(L) = PIJL(I,J,L)
             end do

             !// Note that for outP > inLevEdgeP(1), no data will be
             !// included. Therefore, I adjust (stretch) input pressure
             !// down to model surface:
             inPE(:) = inLevEdgeP(:)
             if (outLevEdgeP(1) .gt. inLevEdgeP(1)) then
                !// Stretch bottom layer down to model surface
                inPE(1) = outLevEdgeP(1)
             end if

             !// Find the weights
             Weights(:,:) = 0._r8
             call Regrid_Column_Weights(inPE, outLevEdgeP, Weights)

             !// Fetch column data, multiplied by pressure thickness
             do LL = 1, nLev
                in1D(LL) = outPropertyInLev(I,J,LL) &
                     * (inPE(LL) - inPE(LL+1))
             end do

             !// Regrid column by column
             out1D(:) = 0._r8 !// to be calculated
             do L = 1, LPAR
                do LL = 1, nLev
                   out1D(L) = out1D(L) + in1D(LL) * Weights(LL,L)
                end do !// do L = 1, nLev
                !// Divide by model pressure thickness
                out1D(L) = out1D(L) / ( outLevEdgeP(L) - outLevEdgeP(L+1) )
             end do !// do LL = 1, LPAR

             trop_SAD(:,I,J,M) = out1D(:)

          end do
       end do

    end do


    ! Deallocate all local variables
    if ( allocated(inProperty) ) deallocate(inProperty)
    if ( allocated(inXYBOX) ) deallocate(inXYBOX)
    if ( allocated(outPropertyInLev) ) deallocate(outPropertyInLev)
    if ( allocated(Weights) ) deallocate(Weights)
    if ( allocated(in1D) ) deallocate(in1D)
    if ( allocated(inPE) ) deallocate(inPE)

    if ( allocated(inLon) ) deallocate(inLon)
    if ( allocated(inLat) ) deallocate(inLat)
    if ( allocated(inLevP) ) deallocate(inLevP)
    if ( allocated(inTime) ) deallocate(inTime)

    if ( allocated(inLonEdge) )  deallocate(inLonEdge)
    if ( allocated(inLatEdge) )  deallocate(inLatEdge)
    if ( allocated(inLevEdgeP) ) deallocate(inLevEdgeP)

    !// --------------------------------------------------------------------
  end subroutine inSAD
  !// ----------------------------------------------------------------------





  !// ---------------------------------------------------------------------
end module chem_oslo_rates
!//=========================================================================

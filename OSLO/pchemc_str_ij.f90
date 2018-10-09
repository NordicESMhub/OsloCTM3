!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Column driver for Oslo stratospheric chemistry.
!//=========================================================================
module pchemc_str_ij
  !// ----------------------------------------------------------------------
  !// MODULE: pchemc_str_ij
  !// DECRIPTION: Module for OSLO_CHEM_STR_IJ, the column driver for Oslo
  !//             stratospheric chemistry.
  !//
  !// Module containg the columnwise integrator of stratospheric chemistry.
  !// Integrates from LMTROP+1 and upwards. Traditionally it stops at
  !// LPAR-1, but if upper boundary condition is skipped, it should be
  !// calculated to LPAR.
  !// 
  !// Ole Amund Sovde, December 2014 (from .f to .f90),
  !//                  October 2008
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  public
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine OSLO_CHEM_STR_IJ ( &
       !// rate constants dependent on more than temperature
       r_ho2_no2_m, r_no2_no3_m, r_ho2no2_heat, r_n2o5_heat, &
       r_oh_hno3, r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
       r_cl_o2_m, r_cloo_heat, r_clo_clo_m, r_cl2o2_m, r_oh_no2_m, &
       r_op_no2_m, r_oh_oh_m, r_oh_no_m, r_no2_clo_m,  &
       r_op_o2_m, r_o2_h_m, r_no2_bro_m, &
       r_no_ho2_b, &
       !// other parameters
       LMTROP, ZC_LOCAL, JV, &
       T, AIR_MOLEC, &
       LM, ITM, MINTEMP, trsp_idx, Xtrsp_idx, TRACER_ID_MAX, NPHM, ICOL, JCOL, &
       DTS, STEADYST, EULER, QTEST_STR, NCHEM_ITER, &
       LOLD_H2OTREATMENT, &
       !// lightning and aircraft source
       EMISX, &
       !// For adjusting change in O3 in diagnose BTTBCK
       deltaO3, &
       !// PSC rate constants
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
       OxCHEMLOSS, OxCHEMPROD )
    !// ---------------------------------------------------------------------
    !// Original author of the chemistry code in 2-D: Stordal et al., 1985.
    !//
    !// Modifications to 3-d started in June 21, 1993 by MtR...
    !// ... and to CTM21 ... August 9, 1995 -> SCTM-1 in fall/1996.
    !//
    !// Slightly modified to work in the vertical:
    !// The major difference is the IM-loop: It has canged to loop in vertical.
    !// IM: Is therefore the total number of layers, changed to LM.
    !// I: Level when integrating upwards, changed to L.
    !// LMSTRT: Now called LMTROP, since it is the LM of TROPosphere, and it is
    !// therefore not an array but a single integer.
    !// QSSA: Moved into separate module to be accessible in the stratosphere.
    !// IREPMX: Changed name to TRACER_ID_MAX as in parameter file.
    !// NTRACE: Not used; removed.
    !// J,L: Changed to the indices I,J for this column.
    !// NST, J_loop: The looping has changed. Only one loop is carried out,
    !//              NST=1,NCHEM_ITER
    !//              If we want to change the constant J-values, a J-loop may be
    !//              inserted later.
    !// EMISX: New array to treat emissions.
    !// FLY_NOX: Removed; included in EMISX.
    !// LIGHTN: Removed; included in EMISX.
    !//
    !// Ole Amund Sovde, October 2009; October 2008
    !// ---------------------------------------------------------------------
    !// TESTNO2, TESTNO = Used in partitioning NOX
    !// ---------------------------------------------------------------------
    use cmn_precision, only: r8
    use qssa_integrator, only: qssa, qssastr
    use chem_oslo_rates, only: &
         !// Constant reaction rates
         r_od_ch4_a, r_od_ch4_b, r_od_ch4_c, &
         r_od_h2, &
         r_h_ho2_a, r_h_ho2_b, r_h_ho2_c, &
         r_op_hno3,  &
         r_od_cfc11_a, r_od_cfc11_b, &
         r_od_hcfc123, r_od_hcfc141, r_od_hcfc142, &
         r_od_hcl, r_h2o_clono2, r_hcl_clono2, &
         r_oh_bro_a, &
         !// Temperature dependent rates
         r_oh_oh, r_oh_ho2, r_oh_ch4, r_oh_h2o2, r_oh_h2, &
         r_oh_ch2o, r_oh_ho2no2, r_oh_ch3oh, &
         r_oh_clo_a, r_oh_clo_b, r_oh_clono2_a, r_oh_clono2_b, &
         r_oh_hocl, r_oh_hcl, r_oh_ch3cl, &
         r_oh_ch3ccl3, r_oh_chclf2, &
         r_oh_bro_b, r_oh_hbr, r_oh_br2, r_oh_ch3br,  &
         r_oh_hcfc123, r_oh_hcfc141, r_oh_hcfc142, &
         r_op_o3, r_n_no, r_n_o2, &
         r_op_no2, r_op_clo, r_op_hcl, r_op_clono2, r_op_oh, r_op_ho2, &
         r_op_ch2o, r_op_oclo, r_op_hbr, r_op_bro, &
         r_od_n2, r_od_o2, r_od_h2o, r_od_n2o_a, r_od_n2o_b, &
         r_o3_no, r_o3_no2, r_o3_h, r_o3_oh, r_o3_ho2, r_o3_cl, &
         r_o3_bro, &
         r_no_no3, r_no_ho2, r_no_ch3o2, r_no_clo, r_no_bro, r_no_oclo, &
         r_no2_no3_b, &
         r_oh_ch3o2h_a, r_oh_ch3o2h_b, &
         r_ho2_ch3o2, &
         r_ho2_cl_a, r_ho2_cl_b, r_ho2_clo, &
         r_cl_h2, r_cl_h2o2, r_cl_ch4, r_cl_ch2o, r_cl_clono2, r_cl_ch3oh, &
         r_br_o3, r_br_h2o2, r_br_ch2o, r_br_ho2, &
         r_clo_co, r_bro_ho2, &
         r_bro_clo_a, r_bro_clo_b, r_bro_clo_c, &
         r_bro_bro_a, r_bro_bro_b, &
         !// Branching ratios
         fb_hv_ohcl
    !// ---------------------------------------------------------------------
    implicit none
    !// ---------------------------------------------------------------------

    !// In/Out parameters
    !// ---------------------------------------------------------------------
    integer, intent(in)   :: &
         ITM, LM, TRACER_ID_MAX, NPHM, & !// Array sizes
         MINTEMP, LMTROP, NCHEM_ITER, &
         ICOL, JCOL, &                   !// Global box indices
         trsp_idx(TRACER_ID_MAX), &           !// Transport nr
         Xtrsp_idx(TRACER_ID_MAX)             !// Non-transport nr

    real(r8), intent(inout) :: ZC_LOCAL(TRACER_ID_MAX,LM)     !// The tracer concentrations
    real(r8), intent(out) :: deltaO3(LM) !// For correcting BTTBCK
    real(r8), intent(in)    :: &
         JV(NPHM,LM), &             !// J-values
         AIR_MOLEC(LM), &           !// air molecules
         EMISX(TRACER_ID_MAX,LM), & !// emissions
         T(LM)                      !//Temperature

    logical, intent(in) :: LOLD_H2OTREATMENT

    !// Parameters used in the QSSA-method
    real(r8), intent(in) :: DTS,STEADYST,EULER,QTEST_STR

    !// Rate constants dep. on more than T (dimension LM)
    real(r8), dimension(LM), intent(in) :: &
         r_ho2_no2_m, &
         r_no2_no3_m, r_ho2no2_heat, &
         r_n2o5_heat, &
         r_oh_hno3, &
         r_oh_co_a, r_oh_co_b, r_ho2_ho2_tot, &
         r_cl_o2_m, &
         r_cloo_heat, &
         r_clo_clo_m, r_cl2o2_m, &
         r_oh_no2_m, r_oh_oh_m, r_oh_no_m, &
         r_op_no2_m, &
         r_no2_clo_m, r_op_o2_m, r_o2_h_m, &
         r_no2_bro_m, r_no_ho2_b
    !// PSC rate constants (dimension LM)
    real(r8), dimension(LM), intent(in) :: &
         spsADPARBK, spsADPAR, &
         spsACPARBK, spsACPAR, &
         spsBDPARBK, spsBDPAR, &
         spsBCPARBK, spsBCPAR, &
         spsECPARBK, spsECPAR, &
         spsFDPARBK, spsFDPAR, &
         spsFCPARBK, spsFCPAR, &
         spsGCPARBK, spsGCPAR, &
                     spsBHPAR
    !// Chemistry diagnoses
    integer, intent(in) :: nchemdiag
    real(r8),dimension(nchemdiag,TRACER_ID_MAX,LM), intent(out) :: CHEMLOSS
    real(r8),dimension(nchemdiag,TRACER_ID_MAX,LM), intent(out) :: CHEMPROD
    real(r8),dimension(LM), intent(out) :: OxCHEMLOSS
    real(r8),dimension(LM), intent(out) :: OxCHEMPROD

    !// Local variables
    !// ---------------------------------------------------------------------
    real(r8) :: &
         !// Local rate constants dep. on T only (extracted from 1D ITM arrays)
         k_n_no, k_n_o2, &
         k_od_cfc11_a, k_od_cfc11_b, &
         k_od_hcfc123, k_od_hcfc141, k_od_hcfc142, &
         k_od_ch4_a, k_od_ch4_b, k_od_ch4_c, &
         k_od_h2, k_od_hcl, k_od_n2, k_od_o2, k_od_h2o, &
         k_od_n2o_a, k_od_n2o_b, &
         k_op_o3, k_op_hno3, &
         k_op_oh, k_op_ho2, k_op_ch2o, &
         k_op_no2, k_op_clo, k_op_clono2, k_op_oclo, k_op_hcl, &
         k_op_bro, k_op_hbr, &
         k_h_ho2_a, k_h_ho2_b, k_h_ho2_c, &
         k_oh_h2o2, k_oh_h2, k_oh_ch2o, k_oh_ho2no2, k_oh_ch4, &
         k_oh_ch3o2h_a, k_oh_ch3o2h_b, k_oh_ch3oh, &
         k_oh_clo_a, k_oh_clo_b, k_oh_clono2_a, k_oh_clono2_b, &
         k_oh_hocl, k_oh_hcl, k_oh_ch3cl, &
         k_oh_hcfc123, k_oh_hcfc141, k_oh_hcfc142, &
         k_h2o_clono2, k_hcl_clono2, &
         k_oh_bro_a, k_oh_bro_b, k_oh_hbr, k_oh_br2, k_oh_ch3br, &
         k_o3_no, k_o3_no2, k_o3_h, k_o3_oh, k_o3_ho2, k_o3_cl, &
         k_no_no3, k_no_ho2, k_no_ch3o2, k_no_clo, k_no_bro, k_no_oclo, &
         k_no2_no3_b, k_oh_oh, k_oh_ho2, &
         k_oh_ch3ccl3, k_oh_chclf2, &
         k_ho2_ch3o2, &
         k_ho2_cl_a, k_ho2_cl_b, k_ho2_clo, &
         k_cl_h2, k_cl_h2o2, k_cl_ch4, k_cl_ch2o, k_cl_clono2, &
         k_br_o3, k_br_h2o2, k_br_ch2o, k_br_ho2, &
         k_clo_co, k_bro_ho2, &
         k_bro_clo_a, k_bro_clo_b, k_bro_clo_c, &
         k_bro_bro_a, k_bro_bro_b, k_o3_bro, &
         !// Local rate constants dep. on more than T
         k_ho2_no2_m, k_no2_no3_m, k_ho2no2_heat, k_n2o5_heat, &
         k_oh_hno3, k_oh_co_a, k_oh_co_b, k_ho2_ho2_tot, &
         k_cl_o2_m, k_cloo_heat, k_clo_clo_m, k_cl2o2_m, k_cl_ch3oh, &
         k_oh_no2_m, k_oh_oh_m, k_oh_no_m, &
         k_op_no2_m, k_op_o2_m, &
         k_no2_clo_m, k_o2_h_m, &
         k_no2_bro_m, k_no_ho2_b, &


         LC_ADPAR, &
         LC_BCPAR,   LC_ACPAR,   LC_BDPAR,  LC_ADPARBK, LC_BCPARBK, &
         LC_ACPARBK, LC_BDPARBK, LC_CEPAR,  LC_CEPARBK, LC_FDPARBK, &
         LC_GCPARBK, &
         !// PSC rate constants
         LC_spsADPARBK, LC_spsADPAR, &
         LC_spsACPARBK, LC_spsACPAR, &
         LC_spsBDPARBK, LC_spsBDPAR, &
         LC_spsBCPARBK, LC_spsBCPAR, &
         LC_spsECPARBK, LC_spsECPAR, &
         LC_spsFDPARBK, LC_spsFDPAR, &
         LC_spsFCPARBK, LC_spsFCPAR, &
         LC_spsGCPARBK, LC_spsGCPAR, &
                        LC_spsBHPAR, &
         !// J-values (extracted from JV array)
         J_O2,      J_O3_a,    J_O3_b,     J_NO2,      J_H2O2, &
         J_HNO3,    J_CH3O2H,  J_NO3_a,    J_NO3_b,    J_CH2O_a, &
         J_CH2O_b,  J_N2O5,    J_HO2NO2_a, J_HO2NO2_b, J_NO, &
         J_HNO2,    J_ClONO2,  J_Cl2,      J_HOCl,     J_OClO, &
         J_Cl2O2,   J_ClO,     J_BrO,      J_BrONO2,   J_HOBr, &
         J_N2O,     J_CFC11,   J_CFC12,    J_CFC113,   J_CFC114, &
         J_CFC115,  J_CCl4,    J_CH3Cl,    J_MCF,      J_CH3Br, &
         J_H1211,   J_H1301,   J_H2402,    J_HCFC22,   J_HCFC123, &
         J_HCFC141, J_HCFC142, J_CHBr3,    J_HBr,      J_Br2, &
         J_HCl,     J_H2O,     J_BrCl, &
         !// Components (extracted from ZC_LOCAL array)
         M_O3,        M_HNO3,     M_CO,      M_CH2O,      M_H2O2, &
         M_CH3O2H,    M_HO2NO2,   M_HO2,     M_CH3O2,     M_O3P, &
         M_O1D,       M_OH,       M_NO3,     M_N2O5,      M_NO, &
         M_NO2,       M_CH4,      M_MCF,     M_HCFC22,    M_CFC11, &
         M_CFC12,     M_CCl4,     M_CH3Cl,   M_N2O,       M_Clx, &
         M_NOx_str,   M_SO,       M_HCl,     M_Cly,M_H2,  M_H2O, &
         M_SH,        M_CH3Br,    M_H1211,   M_H1301,     M_Bry, &
         M_H2402,     M_CFC113,   M_CFC114,  M_CFC115,    M_HNO3s, &
         M_HCFC123,   M_HCFC141,  M_HCFC142, &
         M_H,M_HNO2,  M_Cl,M_ClO, M_HOCl,    M_ClONO2,    M_Cl2, &
         M_OClO,M_Br, M_BrO,      M_HBr,     M_BrONO2,    M_HOBr, &
         M_Br2,       M_ClOO,     M_Cl2O2,   M_BrCl,      M_NOy_str, &
         !// Other components
         M_O2, M_N2,  M_H2Os, M_H2O_ac

    !// Indices
    integer :: N, L, JTEMP, NST

    !// Integration variables
    real(r8) :: &
         PROD, LOSS, PROD_Cly, PROD_Clx, LOSS_HCl_het, PROD_Bry, &
         pHOx,qHOX,RHOX, pNO2,qNO2, &
         pBr2,qBr2, pOHBr,qOHBr, pBrNO3,qBrNO3, &
         pHBr,qHBr, pBrz,qBrz, pBr,qBr, pBrO,qBrO,  &
         pHNO3,qHNO3, pNOx,qNOx

    real(r8) :: &
         dpSH,dqSH,drSH, DO3,drkO3, y3,testfrac, &
         aNO3, hoxtest

    !// Used in scalings...
    real(r8) :: XBrX,YBrX, XVNOX, VNOy, FACN, XClx,YClx, SOhelp, xCly,yCly, xNOy,yNOy
    real(r8) :: HOx,xOH,xHO2,xH2O2,aH,aOH,aHO2,yOH,yHO2
    real(r8) :: xHO2NO2,zHO2NO2,HOxe,hside,HOxu,HOxyy,HOy,yH
    real(r8) :: xNO,xNO2,xNO3,xN2O5,xNOz,yNO3,yN2O5,zNO,zNO2,zN2O5,zNO3, &
         vNOz,yNOxnew,yNOxold,NOy
    real(r8) :: xCl,xClO,xOHCl,xClONO2,yClO, Clz,delClx,ClONO2g,yyClx
    real(r8) :: BrNO3x,Brzx,Brtot,HBrx,OHBrx,BrOx,Brz,Brx,Br2x,Brtst
    real(r8) :: fac,sc,scx,sq,x3,test,testNO,testNO2,scal,&
         tmp_O3,tmp_O3P,tmp_Y,tmp_YY


    !// Used in NOx-chemistry:
    real(r8) :: TR1,TR2
    real(r8) :: TOPNO2,AFAC,QNODIS

    real(r8) :: BrClx, BrClold

    !// Scaling index and moderator
    integer :: ISCALE, scalemod

    !// Parameterized washout rates: No washout in stratosphere.

    !// PSC calculations
    real(r8) :: pnoy, qnoy, phno3s, hno3tonox, rnoxtohno3, &
         dpkclx_het, old_v1, ktot,kpart1,kpart2,kpart3,kpart4,kpart5
    !// Debug stratospheric chemistry
    logical, parameter :: LDEBUG_STRCHEM = .false.
    !// ---------------------------------------------------------------------


    !// Define how often to do scaling
    !// Scaling is done for first NST step and then at the end
    !// according to scalemod. Should be similar to CTM2.
    if (NCHEM_ITER .eq. 12) then
       !// 1-hour time step, scale every 20 min, as in old code
       scalemod = 3
    else if (NCHEM_ITER .eq. 6) then
       !// 30-min time step, scale every 15 min
       scalemod = 2
    else if (NCHEM_ITER .le. 4) then
       !// 15-min time step or less, scale once
       scalemod = 1
    else
       print*, 'Weird scaling in stratchem'
       print*, 'DT',DTS
       print*, 'NCHEM_ITER',NCHEM_ITER
       stop
    end if

    !// Check for negative concentrations & NaNs:
    if (LDEBUG_STRCHEM) then
       do L = LMTROP+1, LM-1
          do N = 1,TRACER_ID_MAX
             if (trsp_idx(N) .gt. 0 .or. Xtrsp_idx(N) .gt. 0) then
                if (ZC_LOCAL(N,L) .lt. 0._r8) then
                   write(*,'(a,5i4,f20.6)') &
                        'OSLO_CHEM_STR_IJ: 1: NEG.CONC. FOR N,I,J,L: ', &
                        N,ICOL,JCOL,L,lmtrop,ZC_LOCAL(N,L)
                   stop
                end if
                if (1._r8*ZC_LOCAL(N,L) .ne. ZC_LOCAL(N,L)) then
                   print*,'NaNQ FOR N,I,J,L : ',N,ICOL,JCOL,L
                   stop
                end if
             end if
          end do
       end do
    end if

    !// Initialise
    TESTNO2 = 10._r8 !// Used in integrations
    TESTNO = 10._r8 !// Used in integrations
    deltaO3(:) = 0._r8
    CHEMLOSS(:,:,:) = 0._r8
    CHEMPROD(:,:,:) = 0._r8
    OxCHEMLOSS(:) = 0._r8
    OxCHEMPROD(:) = 0._r8

    !// Assign constant rates (could use constants directly, but with this
    !// setup it is easy to change the rate to some dependency).
    k_op_hno3 = r_op_hno3
    k_od_ch4_a = r_od_ch4_a
    k_od_ch4_b = r_od_ch4_b
    k_od_ch4_c = r_od_ch4_c
    k_od_h2    = r_od_h2
    k_od_cfc11_a = r_od_cfc11_a
    k_od_cfc11_b = r_od_cfc11_b
    k_od_hcfc123 = r_od_hcfc123
    k_od_hcfc141 = r_od_hcfc141
    k_od_hcfc142 = r_od_hcfc142
    k_od_hcl = r_od_hcl
    k_oh_bro_a = r_oh_bro_a
    k_h_ho2_a = r_h_ho2_a
    k_h_ho2_b = r_h_ho2_b
    k_h_ho2_c = r_h_ho2_c
    k_h2o_clono2 = r_h2o_clono2
    k_hcl_clono2 = r_hcl_clono2

    !// ---------------------------------------------------------------------
    do L = LMTROP+1, LM-1
    !// ---------------------------------------------------------------------

      !// the constant reaction rates are given as arguments for the subroutine
      JTEMP = nint(T(L)) - MINTEMP

      !// Assign chemical reaction rates, dependent on T only
      !// (prefix 'LC' for local constant)
      k_n_no   = r_n_no(JTEMP)
      k_n_o2   = r_n_o2(JTEMP)
      k_op_o3   = r_op_o3(JTEMP)
      k_op_oh = r_op_oh(JTEMP)
      k_op_ho2  = r_op_ho2(JTEMP)
      k_op_ch2o = r_op_ch2o(JTEMP)
      k_op_no2 = r_op_no2(JTEMP)
      k_op_clo = r_op_clo(JTEMP)
      k_op_hcl = r_op_hcl(JTEMP)
      k_op_clono2 = r_op_clono2(JTEMP)
      k_op_hbr = r_op_hbr(JTEMP)
      k_od_n2   = r_od_n2(JTEMP)
      k_od_o2   = r_od_o2(JTEMP)
      k_od_n2o_a = r_od_n2o_a(JTEMP)
      k_od_n2o_b = r_od_n2o_b(JTEMP)
      k_od_h2o = r_od_h2o(JTEMP)
      k_o3_no  = r_o3_no(JTEMP)
      k_o3_no2 = r_o3_no2(JTEMP)
      k_o3_oh = r_o3_oh(JTEMP)
      k_o3_ho2 = r_o3_ho2(JTEMP)
      k_o3_cl  = r_o3_cl(JTEMP)
      k_o3_h = r_o3_h(JTEMP)
      k_no_ho2  = r_no_ho2(JTEMP)
      k_no_ch3o2 = r_no_ch3o2(JTEMP)
      k_no_no3 = r_no_no3(JTEMP)
      k_no_clo = r_no_clo(JTEMP)
      k_no2_no3_b = r_no2_no3_b(JTEMP)
      k_oh_oh = r_oh_oh(JTEMP)
      k_oh_ho2 = r_oh_ho2(JTEMP)
      k_oh_h2 = r_oh_h2(JTEMP)
      k_oh_h2o2 = r_oh_h2o2(JTEMP)
      k_oh_ch2o = r_oh_ch2o(JTEMP)
      k_oh_ch3oh = r_oh_ch3oh(JTEMP)
      k_oh_ch4 = r_oh_ch4(JTEMP)
      k_oh_ho2no2 = r_oh_ho2no2(JTEMP)
      k_oh_ch3o2h_a = r_oh_ch3o2h_a(JTEMP)
      k_oh_ch3o2h_b = r_oh_ch3o2h_b(JTEMP)
      k_oh_hcl = r_oh_hcl(JTEMP)
      k_oh_hocl = r_oh_hocl(JTEMP)
      k_oh_ch3cl = r_oh_ch3cl(JTEMP)
      k_oh_clo_a = r_oh_clo_a(JTEMP)
      k_oh_clo_b = r_oh_clo_b(JTEMP)
      k_oh_clono2_a = r_oh_clono2_a(JTEMP)
      k_oh_clono2_b = r_oh_clono2_b(JTEMP)
      k_oh_ch3ccl3 = r_oh_ch3ccl3(JTEMP)
      k_oh_chclf2 = r_oh_chclf2(JTEMP)
      k_ho2_clo = r_ho2_clo(JTEMP)
      k_ho2_ch3o2 = r_ho2_ch3o2(JTEMP)
      k_ho2_cl_a = r_ho2_cl_a(JTEMP)
      k_ho2_cl_b = r_ho2_cl_b(JTEMP)

      k_cl_h2 = r_cl_h2(JTEMP)
      k_cl_h2o2 = r_cl_h2o2(JTEMP)
      k_cl_ch2o = r_cl_ch2o(JTEMP)
      k_cl_ch4 = r_cl_ch4(JTEMP)
      k_cl_clono2 = r_cl_clono2(JTEMP)
      k_cl_ch3oh = r_cl_ch3oh(JTEMP)

      k_br_o3  = r_br_o3(JTEMP)
      k_br_ho2 = r_br_ho2(JTEMP)
      k_br_h2o2  = r_br_h2o2(JTEMP)
      k_br_ch2o  = r_br_ch2o(JTEMP)

      k_clo_co = r_clo_co(JTEMP)

      k_op_bro = r_op_bro(JTEMP)
      k_bro_ho2 = r_bro_ho2(JTEMP)
      k_bro_clo_a = r_bro_clo_a(JTEMP)
      k_bro_clo_b = r_bro_clo_b(JTEMP)
      k_bro_clo_c = r_bro_clo_c(JTEMP)
      k_no_bro = r_no_bro(JTEMP)
      k_bro_bro_a = r_bro_bro_a(JTEMP)
      k_bro_bro_b = r_bro_bro_b(JTEMP)
      k_o3_bro = r_o3_bro(JTEMP)
      k_oh_bro_b = r_oh_bro_b(JTEMP)
      k_oh_hbr = r_oh_hbr(JTEMP)
      k_oh_br2 = r_oh_br2(JTEMP)

      k_no_oclo = r_no_oclo(JTEMP)
      k_op_oclo = r_op_oclo(JTEMP)

      k_oh_ch3Br = r_oh_ch3Br(JTEMP)

      k_oh_hcfc123 = r_oh_hcfc123(JTEMP)
      k_oh_hcfc141 = r_oh_hcfc141(JTEMP)
      k_oh_hcfc142 = r_oh_hcfc142(JTEMP)


      !// Assign chemical reaction rates, dependent on T and P
      k_ho2_no2_m = r_ho2_no2_m(L)
      k_no2_no3_m = r_no2_no3_m(L)
      k_ho2no2_heat = r_ho2no2_heat(L)
      k_n2o5_heat = r_n2o5_heat(L)
      k_oh_hno3 = r_oh_hno3(L)
      k_oh_co_a = r_oh_co_a(L)
      k_oh_co_b = r_oh_co_b(L)
      k_ho2_ho2_tot = r_ho2_ho2_tot(L)
      k_cl_o2_m = r_cl_o2_m(L)
      k_cloo_heat = r_cloo_heat(L)
      k_clo_clo_m = r_clo_clo_m(L)
      k_cl2o2_m = r_cl2o2_m(L)
      k_oh_no2_m = r_oh_no2_m(L)
      k_op_no2_m = r_op_no2_m(L)
      k_oh_oh_m = r_oh_oh_m(L)
      k_oh_no_m = r_oh_no_m(L)
      k_no2_clo_m = r_no2_clo_m(L)
      k_op_o2_m = r_op_o2_m(L)
      k_o2_h_m = r_o2_h_m(L)
      k_no2_bro_m = r_no2_bro_m(L)
      k_no_ho2_b = r_no_ho2_b(L)


      !// PSC rate constants
      LC_spsADPARBK = spsADPARBK(L)
      LC_spsACPARBK = spsACPARBK(L)
      LC_spsBDPARBK = spsBDPARBK(L)
      LC_spsBCPARBK = spsBCPARBK(L)
      LC_spsECPARBK = spsECPARBK(L)
      LC_spsFDPARBK = spsFDPARBK(L)
      LC_spsFCPARBK = spsFCPARBK(L)
      LC_spsGCPARBK = spsGCPARBK(L)
      LC_spsADPAR = spsADPAR(L)
      LC_spsACPAR = spsACPAR(L)
      LC_spsBDPAR = spsBDPAR(L)
      LC_spsBCPAR = spsBCPAR(L)
      LC_spsECPAR = spsECPAR(L)
      LC_spsFDPAR = spsFDPAR(L)
      LC_spsFCPAR = spsFCPAR(L)
      LC_spsGCPAR = spsGCPAR(L)
      LC_spsBHPAR = spsBHPAR(L)
      !// Sometimes there is very high conversion on PSCs, but this is not
      !// due to large rates; rather it is due to large surface area
      !// density calculated by microphysics.

      !// ------------------------------------------------------------------
      !// In the old days, the looping was divided into a NST loop and
      !// a loop for constant J-values. As the integration is carried
      !// out in vertical, the J-values are constant, and we loop
      !// over NST = 1, NCHEM_ITER.
      !// ------------------------------------------------------------------

      !// Generally, QSSA updates the M_* variables, with some exceptions.
      !// This is somewhat different from the troposphere, where they
      !// are only updated for short-lived species.

      !// Loop through the number of iterations to match the global DTCHM.
      do NST = 1, NCHEM_ITER

        !// Assign chemical components (prefix 'M' for molecular density)
        M_O3      = ZC_LOCAL( 1,L)
        M_HNO3    = ZC_LOCAL( 4,L)
        M_CO      = ZC_LOCAL( 6,L)
        M_CH2O    = ZC_LOCAL(13,L)
        M_H2O2    = ZC_LOCAL(15,L)
        M_CH3O2H  = ZC_LOCAL(16,L)
        M_HO2NO2  = ZC_LOCAL(17,L)
        M_HO2     = ZC_LOCAL(21,L)
        M_CH3O2   = ZC_LOCAL(22,L)
        M_O3P     = ZC_LOCAL(38,L)
        M_O1D     = ZC_LOCAL(39,L)
        M_OH      = ZC_LOCAL(40,L)
        M_NO3     = ZC_LOCAL(41,L)
        M_N2O5    = ZC_LOCAL(42,L)
        M_NO      = ZC_LOCAL(43,L)
        M_NO2     = ZC_LOCAL(44,L)
        M_CH4     = ZC_LOCAL(46,L)
        M_MCF     = ZC_LOCAL(101,L)
        M_HCFC22  = ZC_LOCAL(102,L)
        M_CFC11   = ZC_LOCAL(103,L)
        M_CFC12   = ZC_LOCAL(104,L)
        M_CCl4    = ZC_LOCAL(105,L)
        M_CH3Cl   = ZC_LOCAL(106,L)
        M_N2O     = ZC_LOCAL(107,L)
        M_Clx     = ZC_LOCAL(108,L)
        M_NOx_str = ZC_LOCAL(109,L)
        M_SO      = ZC_LOCAL(110,L)
        M_HCl     = ZC_LOCAL(111,L)
        M_Cly     = ZC_LOCAL(112,L)
        M_H2      = ZC_LOCAL(113,L)
        !// Make/Update the water vapour distribution:
        !// In the stratosphere it is recalcutated based on 'H2O+2CH4+H2=const.'
        !// before stratospheric chemistry is called.
        M_H2O     = ZC_LOCAL(114,L)
        M_SH      = ZC_LOCAL(115,L)
        M_CH3Br   = ZC_LOCAL(116,L)
        M_H1211   = ZC_LOCAL(117,L)
        M_H1301   = ZC_LOCAL(118,L)
        M_Bry     = ZC_LOCAL(119,L)
        M_H2402   = ZC_LOCAL(120,L)
        M_CFC113  = ZC_LOCAL(121,L)
        M_CFC114  = ZC_LOCAL(122,L)
        M_CFC115  = ZC_LOCAL(123,L)
        M_HNO3s   = ZC_LOCAL(124,L)
        M_H2Os    = ZC_LOCAL(125,L)
        M_HCFC123 = ZC_LOCAL(127,L)
        M_HCFC141 = ZC_LOCAL(128,L)
        M_HCFC142 = ZC_LOCAL(129,L)
        M_H       = ZC_LOCAL(130,L)
        M_HNO2    = ZC_LOCAL(131,L)
        M_Cl      = ZC_LOCAL(132,L)
        M_ClO     = ZC_LOCAL(133,L)
        M_HOCl    = ZC_LOCAL(134,L)
        M_ClONO2  = ZC_LOCAL(135,L)
        M_Cl2     = ZC_LOCAL(136,L)
        M_OClO    = ZC_LOCAL(137,L)
        M_Br      = ZC_LOCAL(138,L)
        M_BrO     = ZC_LOCAL(139,L)
        M_HBr     = ZC_LOCAL(140,L)
        M_BrONO2  = ZC_LOCAL(141,L)
        M_HOBr    = ZC_LOCAL(142,L)
        M_Br2     = ZC_LOCAL(143,L)
        M_ClOO    = ZC_LOCAL(144,L)
        M_Cl2O2   = ZC_LOCAL(145,L)
        M_BrCl    = ZC_LOCAL(146,L)
        M_NOy_str = ZC_LOCAL(147,L)
        M_H2O_ac  = ZC_LOCAL(148,L)
        !// Molecular oxygen concentration from dry air density...
        M_O2 = 0.2095_r8 * AIR_MOLEC(L)
        M_N2 = 0.7809_r8 * AIR_MOLEC(L)

        !// The calculations are done columnwise, and even though the sun moves
        !// 5 degrees in 20 minutes, we use the same J-value for the loops.
        !// If you want new J-values, they should be calculated first.
        !// Check the user manual for info on this.

        !// Assign photolysis rates (prefix 'J' for J-value)
        J_O2       = JV(21,L) ! O2 + hv     -> O + O
        J_O3_a     = JV( 1,L) ! O3 + hv     -> O2 + O(3P)
        J_O3_b     = JV( 2,L) ! O3 + hv     -> O2 + O(1D)
        J_NO2      = JV( 3,L) ! NO2 + hv    -> NO + O(3P)
        J_H2O2     = JV( 4,L) ! H2O2 + hv   -> 2OH
        J_HNO3     = JV( 5,L) ! HNO3 + hv   -> NO2 + OH
        J_CH3O2H   = JV(14,L) ! CH3O2H + hv -> CH2O + OH + H (->HO2)
        J_NO3_a    = JV(11,L) ! NO3 + hv    -> NO + O2
        J_NO3_b    = JV(12,L) ! NO3 + hv    -> NO2 + O(3P)
        J_CH2O_a   = JV( 6,L) ! CH2O + hv   -> H + CHO
        J_CH2O_b   = JV( 7,L) ! CH2O + hv   -> H2 + CO
        J_N2O5     = JV(13,L) ! N2O5 + hv   -> NO2 + NO3
        J_HO2NO2_a = JV(16,L) ! HO2NO2 + hv -> OH + NO3
        J_HO2NO2_b = JV(17,L) ! HO2NO2 + hv -> HO2 + NO2
        J_NO       = JV(22,L) ! NO + hv     -> N + O
        J_HNO2     = JV(23,L) ! HNO2 + hv   -> OH + NO (or H + NO2)
        J_ClONO2   = JV(24,L) ! ClONO2 + hv -> NO3 + Cl
        J_Cl2      = JV(25,L) ! Cl2 + hv    -> 2Cl
        J_HOCl     = JV(26,L) ! HOCl + hv   -> OH + Cl
                              !             -> HCl + O(3P)
        J_OClO     = JV(27,L) ! OClO + hv   -> ClO + O(3P)
        J_Cl2O2    = JV(28,L) ! Cl2O2 + hv  -> ClOO + Cl
        J_ClO      = JV(29,L) ! ClO + hv    -> Cl + O
        J_BrO      = JV(30,L) ! BrO + hv    -> Br + O
        J_BrONO2   = JV(31,L) ! BrONO2 + hv -> NO2 + BrO (0.71)
                              !             -> NO3 + Br  (0.29)
        J_HOBr     = JV(32,L) ! HOBr + hv   -> OH + Br
        J_N2O      = JV(33,L) ! N2O + hv    -> N2 + O(1D)
        J_CFC11    = JV(34,L) ! CFC11 + hv  -> prod.
        J_CFC12    = JV(35,L) ! CFC12 + hv  -> prod.
        J_CFC113   = JV(36,L) ! CFC113 + hv -> prod.
        J_CFC114   = JV(37,L) ! CFC114 + hv -> prod.
        J_CFC115   = JV(38,L) ! CFC115 + hv -> prod.
        J_CCl4     = JV(39,L) ! CCl4 + hv   -> prod.
        J_CH3Cl    = JV(40,L) ! CH3Cl + hv  -> prod.
        J_MCF      = JV(41,L) ! CH3CCl3 + hv -> prod.
        J_CH3Br    = JV(42,L) ! CH3Br + hv   -> prod. (CH3 + Br)
        J_H1211    = JV(43,L) ! CF2ClBr + hv -> prod.
        J_H1301    = JV(44,L) ! CF3Br + hv   -> prod.
        J_H2402    = JV(45,L) ! CF2BrCF2Br + hv -> prod.
        J_HCFC22   = JV(46,L) ! CHF2Cl + hv     -> prod.
        J_HCFC123  = JV(47,L) ! CF3CHCl2 + hv   -> prod.
        J_HCFC141  = JV(48,L) ! CF3CFCl2 + hv   -> prod.
        J_HCFC142  = 0._r8 ! ??? ! CF3CF2Cl + hv   -> prod.
        J_CHBr3    = JV(49,L) ! CHBr3 + hv      -> prod.
        J_HBr      = 0._r8 ! ??? ! HBr + hv -> H + Br
        J_Br2      = 0._r8 ! ??? ! Br2 + hv -> 2Br
        J_HCl      = 0._r8 ! ??? ! HCl + hv -> H + Cl
        if (LOLD_H2OTREATMENT) then
           J_H2O   = 0._r8     ! H2O + hv -> H + OH
        else
           J_H2O   = JV(52,L) ! H2O + hv -> H + OH
        end if
        !// When daylight --> J_BrCl=.1, else: J_BrCl=0
        if (J_O3_a .gt. 0._r8) then
           J_BrCl = 0.1_r8
           J_Br2  = 0.1_r8
        else
           J_BrCl = 0._r8
           J_Br2  = 0._r8
        end if
          

        !// ---------------------------------------------------------------
        !// CHEMISTRY CALCULATIONS START
        !// ---------------------------------------------------------------

        !// ---------------------------------------------------------------
        !// Scaling over families when starting chemistry right
        !// after transport. 950929, MtR, UiO/FMI
        !//
        !// The idea behind this scaling is to make sure that the sum of
        !// the concentrations of the components in a 'chemical family'
        !// does not exceed the concentration of the 'family' itself,
        !// calculated independently.
        !// As the family is more stable than the members of it, separate
        !// integration is performed for it and this result is taken as
        !// the new value. The members are then scaled with the ratio
        !// of this (individually integrated) family and the sum of
        !// the (individually integrated) members of the family.
        !//
        !// Bry = Br + BrO + BrONO2 + OHBr + HBr + 2*Br2 + BrCl
        !// SO  = O3 + OD + OP - NO - Cl - Br
        !// NOx_str = NO + NO2 + NO3 + 2*N2O5 + HNO2 + HO2NO2 + ClONO2 + BrONO2
        !// Clx = Cl + ClO + OHCl + ClONO2 + 2*Cl2 + OClO + BrCl + ClOO + 2*Cl2O2
        !// Cly = Clx + HCl
        !// NOy_str = NOx_str + HNO3
        !// ---------------------------------------------------------------

        !// Scale if this is the first loop with the latest J-values
        if (NST .eq. 1) then

          !// Do an iterative scaling, i.e. go through it three times
          do ISCALE = 1, 3

            !// Brx:
            xBrx = M_Bry 
            yBrx = M_Br + M_BrO + M_BrONO2 + M_HOBr + M_HBr &
                   + 2._r8*M_Br2 + M_BrCl

            FACN = xBrx/yBrx

            M_Br   = M_Br * FACN
            M_BrO  = M_BrO * FACN
            M_HBr  = M_HBr * FACN
            M_BrONO2 = M_BrONO2 * FACN
            M_HOBr = M_HOBr * FACN
            M_Br2  = M_Br2 * FACN
            M_BrCl = M_BrCl * FACN

            !// NOx:
            VNOy = M_NOx_str
            xVNOx = M_NO + M_NO2 + M_NO3 + 2._r8*M_N2O5 + M_HO2NO2 &
                    + M_ClONO2 + M_BrONO2

            FACN = VNOy/XVNOx

            M_NO   = M_NO * FACN
            M_NO2  = M_NO2 * FACN
            M_NO3  = M_NO3 * FACN
            M_N2O5 = M_N2O5 * FACN
            M_HO2NO2 = M_HO2NO2 * FACN
            M_ClONO2 = M_ClONO2 * FACN
            M_BrONO2 = M_BrONO2 * FACN

            !// Clx:
            xClx = M_Clx
            yClx = M_Cl + M_ClO + M_HOCl + M_ClONO2 + 2._r8*M_Cl2 &
                   + M_OClO + M_BrCl + M_ClOO + 2._r8*M_Cl2O2

            FACN = xClx/yClx

            M_Cl    = M_Cl * FACN
            M_ClO   = M_ClO * FACN
            M_HOCl  = M_HOCl * FACN
            M_ClONO2 = M_ClONO2 * FACN
            M_Cl2   = M_Cl2 * FACN
            M_OClO  = M_OClO * FACN
            M_Cl2O2 = M_Cl2O2 * FACN
            M_ClOO  = M_ClOO * FACN
            M_BrCl  = M_BrCl * FACN

          end do !// do ISCALE = 1, 3

          !// Families which do not need iteration:
          !// Cly:
          xCly = M_Cly
          yCly = M_Clx + M_HCl !// Do not calculate solid HCl

          FACN = xCly/yCly

          M_Clx = M_Clx * FACN
          M_HCl = M_HCl * FACN

          !// NOy:
          xNOy = M_NOy_str
          yNOy = M_NOx_str + M_HNO3 !// Solid HNO3 no longer in NOy

          FACN = xNOy/yNOy

          M_NOx_str = M_NOx_str * FACN
          M_HNO3 = M_HNO3 * FACN


          !// -----------------------------------------------------------
          !// Setting steady state O1D, O3P and calculating O3 from SO.
          !// -----------------------------------------------------------
          !// Ole Amund Sovde, May 2015
          !// -----------------------------------------------------------
          !// O1D: Steady state.
          !// O1D is very small, so we may set this using the old O3.
          !// It will not be critical that O3 is updated somewhat below.
          M_O1D = J_O3_b * M_O3 &              ! O3 + hv -> O1D
                  / ( k_od_n2 * M_N2   &       ! O1D + N2 -> O3P + N2
                      + k_od_o2 * M_O2 &       ! O1D + O2 -> O3P + O2
                      + k_od_h2o * M_H2O &     ! O1D + H2O -> 2OH
                      + k_od_n2o_a * M_N2O &   ! O1D + N2O -> N2 + O2
                      + k_od_n2o_b * M_N2O &   ! O1D + N2O -> 2NO
                    )

          !// O3P is larger, and at high altitudes it may be as large or
          !// larger than O3.
          !// In cases where O3P is (too) large, I tried assuming O3+O3P to
          !// be in steady state, but that produced too much O3 for some
          !// reason.
          !// A simpler method is to scale down the calculated O3P to match
          !// keep SO unchanged in the calculation. Note that SO is calculated
          !// further down, so it is allowed to change. If it is reduced too
          !// much, O3 may still become negative, but we assume SO will not
          !// change that much.

          !// O3P: Steady state.
          !// Note that we have to use some old non-transported species
          !// (Br, Cl, OH, ...) and some transported species (HO2, O3, ...).
          !// This is as correct as we can get it; the important thing
          !// is that J-values are updated so that O3P gets as correct
          !// as possible.
          if (J_O3_a .gt. 0._r8) then
             M_O3P = ( J_O3_a * M_O3 &                 ! O3 + hv -> O3P + O2
                       + J_NO2 * M_NO2 &               ! NO2 + hv -> O3P + NO
                       + J_NO3_b * M_NO3 &             ! NO3 + hv -> O3P + NO2
                       + k_od_n2 * M_O1D * M_N2   &    ! O1D + N2 -> O3P + N2
                       + k_od_o2 * M_O1D * M_O2 &      ! O1D + O2 -> O3P + O2
                       + 2._r8 * J_O2 * M_O2 &         ! O2 + hv -> 2O3P
                       + k_oh_oh * M_OH * M_OH &       ! OH + OH -> OP + H2O
                       + k_h_ho2_b * M_H * M_HO2 &     ! H + HO2   -> O + H2O
                     ) &
                     / ( k_op_o2_m * M_O2 &         ! O3P + O2 -M-> O3
                         + k_op_ch2o * M_CH2O &     ! O3P + CH2O -> prod.
                         + k_op_ho2 * M_HO2 &       ! O3P + HO2 -> OH + O2
                         + k_op_oh * M_OH &         ! O3P + OH -> O2 + H
                         + k_op_hcl * M_HCl &       ! O3P + HCl -> OH + Cl
                         + k_op_o3 * M_O3 &         ! O3P + O3 -> 2O2
                         + k_op_clo * M_ClO &       ! O3P + ClO    -> Cl + O2
                         + k_op_clono2 * M_ClONO2 & ! O3P + ClONO2 -> ClO + NO3
                         + k_op_bro * M_BrO &       ! BrO + O(3P)  -> Br + O2
                         + k_op_hbr * M_HBr &       ! O3P + HBr -> OH + Br
                         + k_op_oclo * M_OClO &     ! OClO + O(3P) -> ClO + O2 
                       )
          else
             !// No sunlight, set O3P and O1D to zero
             M_O3P = 0._r8
             M_O1D = 0._r8
          end if



          !// Now we need to update O3 also.
          tmp_O3 = M_O3  !// First; save transported value
          M_O3 = M_SO + M_NO + M_Cl + M_Br - M_O1D - M_O3P

          !// Check if O3 is smaller than O3P or even negative.
          !// If so, an adjustment of O3P must be carried out.
          if (M_O3 .lt. 1.e-2_r8 * M_O3P) then
             !// O3 will likely be negative or at least too small
             !// O3/O3P is about 1.d-6 at 100km altitude.
             !//
             !// This will occur sometimes at L=59, but more often at
             !// L=60, which is at about 65-75km.
             !// Assume O3 to be 0.1*M_O3P and scale down M_O3P to
             !// keep SO unchanged.
             tmp_Y = M_SO + M_NO + M_Cl + M_Br - M_O1D ! O3+O3P old
             tmp_YY = M_O3P + 0.1_r8 * M_O3P            ! Assume O3=0.1*O3P
             !// Scale down O3P
             tmp_O3P = M_O3P
             M_O3P = M_O3P * tmp_Y / tmp_YY

             !// Reset O3
             M_O3 = M_SO + M_NO + M_Cl + M_Br - M_O1D - M_O3P

             if (M_O3 .lt. 1.e-6_r8) then
                write(6,'(a,3i4,4es9.2)') 'pchemc_str_ij.f90: Neg/Tiny O3 A',&
                     L,ICOL,JCOL, M_O3P,M_O3,M_SO,M_O1D
                stop
             end if
          end if


          !// If O3 is updated, we do not want this change to be seen in
          !// BTTBCK, since it will appear as a sink or source.
          deltaO3(L) = M_O3 - tmp_O3

        end if !// if (NST .eq. 1) then
        !// ---------------------------------------------------------------


        !// Used in the loss of SO
        TOPNO2 = 2._r8 * k_op_no2 * M_O3P * M_NO2
        AFAC = J_NO2 / ( J_NO2 + k_o3_no * M_O3 )


        !// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !// Long lived species
        !// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        !// SO: SOURCE terms: (SO, not O3!!!  SO = O3 + OD + OP - NO - Cl - Br)
        PROD = k_no_ch3o2 * M_NO * M_CH3O2      & ! NO + CH3O2 -> CH3O + NO2
               + k_no_no3 * M_NO * M_NO3        & ! NO + NO3 -> 2NO2
               !//HNO2: + k_oh_no_m * M_NO * M_OH  ! NO + OH -> HNO2
               + k_no_ho2 * M_NO * M_HO2        & ! NO + HO2 -> NO2 + OH
               + k_no_ho2_b * M_NO * M_HO2      & ! NO + HO2 -> HNO3
               + k_cl_ch4 * M_Cl * M_CH4        & ! Cl + CH4 -> HCl + CH3
               - k_oh_hcl * M_OH * M_HCl        & ! OH + HCl -> H2O + Cl
               - k_oh_hbr * M_OH * M_HBr        & ! OH + HBr -> H2O + Br
               - 2._r8*k_hcl_clono2 * M_HCl * M_ClONO2& ! HCl + ClONO2 -> 2Cl+HNO3
               - 2._r8*J_Cl2 * M_Cl2            & ! Cl2 + hv -> 2Cl
               + 2._r8*J_O2 * M_O2              & ! O2 + hv  -> 2O
               + J_N2O * M_N2O                  & ! N2O + hv -> N2 + O(1D)
               !//HNO2: - J_HNO2 * M_HNO2 ! HNO2 + hv -> 
               + ( J_NO3_b                      & ! O3 + hv  -> O2 + O(1D)
                   - J_NO3_a ) * M_NO3          & ! NO3 + hv -> NO + O2
               - J_ClONO2 * M_ClONO2            & ! ClONO2 + hv -> NO3 + Cl
               !//CASoldPSC + (1._r8 - fb_hv_ohcl) * J_HOCl * M_HOCl
               - k_oh_clo_a * M_OH * M_ClO      & ! OH + ClO -> Cl + HO2
               - k_oh_bro_b * M_OH * M_BrO      & ! BrO + OH -> prod.= Br + HO2
               - k_cl_clono2 * M_Cl * M_ClONO2  & ! Cl + ClONO2 -> Cl2 + NO3
               + ( k_ho2_cl_a                   & ! HO2 + Cl -> HCl + O2
                   + k_ho2_cl_b ) * M_Cl * M_HO2 & ! HO2 + Cl -> OH + ClO
               + ( k_cl_h2 * M_H2               & ! H2 + Cl  -> HCl + H
                  + k_cl_h2o2 * M_H2O2          & ! H2O2 + Cl -> HCl + HO2
                  + k_cl_ch2o * M_CH2O ) * M_Cl & ! Cl + CH2O -> HCl + HCO
               - k_no2_no3_b * M_NO3 * M_NO2    & ! NO2 + NO3 -> NO + NO2 + O2
               + ( k_br_h2o2 * M_H2O2           & ! Br + H2O2 -> HBr + HO2
                   + k_br_ch2o * M_CH2O         & ! Br + CH2O -> HBr + HCO
                   + k_br_ho2 * M_HO2  ) * M_Br & ! Br + HO2  -> HBr + O2
               + J_OClO * M_OClO                & ! OClO + hv -> ClO + O(3P)
               + k_no_oclo * M_NO * M_OClO      & ! OClO + NO -> NO2 + ClO
               - 2._r8*J_Br2 * M_Br2            & ! Br2 + hv  -> 2Br
               - J_HBr * M_HBr                  & ! HBr + hv  -> H + Br
               - k_bro_clo_a * M_ClO * M_BrO    & ! BrO + ClO -> Br + OClO
               !// Assuming instant thermal decomp. of ClOO, omitting
               !//    - k_cloo_heat * M_ClOO   ! ClOO + M   -> Cl + O2 + M
               !// so that the expression becomes 2 * J_Cl2O2 * M_Cl2O2,
               !// resulting in 2Cl:
               - 2._r8*J_Cl2O2 * M_Cl2O2        & ! Cl2O2 + hv -> ClOO + Cl
               !// Assuming instant thermal decomp. of ClOO, omitting
               !//   - k_cloo_heat * M_ClOO
               !// above, so that the expression becomes - 2 * k_bro_clo_c * M_ClO * M_BrO,
               !// resulting in Br + Cl:
               - 2._r8*k_bro_clo_c * M_ClO * M_BrO & ! BrO + ClO  -> Br + ClOO
               - 2._r8*J_BrCl * M_BrCl          & ! BrCl + hv  -> Br + Cl
               + k_oh_oh * M_OH * M_OH          & ! OH + OH -> OP + H2O
               + k_h_ho2_b * M_H * M_HO2          ! H + HO2   -> O + H2O


        !// SO: LOSS terms:
        LOSS = TOPNO2 * (1._r8 - AFAC) &
             + k_op_no2_m * M_O3P * M_NO2      & ! OP + NO2 -M-> NO3
             + k_op_oh * M_O3P * M_OH          & ! O(3P) + OH  -> O2 + H
             + k_op_ho2 * M_O3P * M_HO2        & ! O(3P) + HO2 -> OH + O2
             + k_od_h2o * M_O1D * M_H2O        & ! O(1D) + H2O -> OH + OH
             + ( k_od_ch4_a                    & ! O(1D) + CH4 -> OH + CH3
                 + k_od_ch4_b  & ! O(1D) + CH4 -> (CH3O/CH2OH + H)  -O2-> CH2O + HO2 + H
                 + k_od_ch4_c                  & ! O(1D) + CH4 -> H2 + CH2O
               ) * M_O1D * M_CH4               &
             + k_o3_h * M_O3 * M_H             & ! O3 + H   -> OH + O2
             + k_o3_no2 * M_O3 * M_NO2         & ! O3 + NO2 -> NO3 + O2
             + k_o3_oh * M_O3 * M_OH           & ! O3 + OH  -> HO2 + O2
             + k_o3_ho2 * M_O3 * M_HO2         & ! O3 + HO2 -> OH + 2O2
             + 2._r8 * k_op_clo * M_O3P * M_ClO  & ! O(3P) + ClO    -> Cl + O2
             + k_op_clono2 * M_O3P * M_ClONO2  & ! O(3P) + ClONO2 -> ClO + NO3
             + 2._r8 * k_op_hcl * M_O3P * M_HCl & ! O(3P) + HCl  -> OH + Cl
             + 2._r8 * k_op_hbr * M_O3P * M_HBr & ! O(3P) + HBr  -> OH + Br
             + 2._r8 * k_od_hcl * M_O1D * M_HCl & ! O(1D) + HCl -> OH + Cl
             + 2._r8 * k_op_bro * M_O3P * M_BrO & ! BrO + O(3P)  -> Br + O2
             + 2._r8 * k_o3_bro * M_O3 * M_BrO  & ! BrO + O3     -> Br + 2O2
             + 2._r8 * k_bro_bro_a * M_BrO * M_BrO & ! BrO + BrO    -> 2Br + O2
             + k_op_oclo * M_O3P * M_OClO      & ! OClO + O(3P) -> ClO + O2 
             + J_HOBR * M_HOBr                 & ! HOBr + hv -> Br + OH
             + J_HOCL * M_HOCl                   ! HOCl + hv -> OH + Cl

        DRKO3 = 2._r8 * k_op_o3 * M_O3 * M_O3P & ! O(3P) + O3 -> O2 + O2
                + TOPNO2 * AFAC

        LOSS = LOSS + DRKO3


        OxCHEMPROD(L) = OxCHEMPROD(L) + PROD*DTS
        OxCHEMLOSS(L) = OxCHEMLOSS(L) + LOSS*DTS

        !// To get correct units, we divide SO loss by the concentration!
        LOSS = LOSS / M_SO

        !// QSSA does not allow for negative concentrations, which may
        !// occur for the SO family. Use the QSSASTR in stead:
        call QSSASTR(201,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_SO,ICOL,JCOL,L)



        !// NOy
        PNOy = 2._r8 * k_od_n2o_b * M_O1D * M_N2O & ! N2O + O(1D) -> 2NO
               !// Emissions sources (in general lightning, aircraft)
               !// HNO3, HO2NO2, NO3, 2*N2O5, NO, NO2, ClONO2, not PANX
               + EMISX(4,L) + EMISX(17,L) + EMISX(41,L) + 2._r8*EMISX(42,L) &
               + EMISX(43,L) + EMISX(44,L) + EMISX(135,L)
        TR1 = k_n_no * M_NO                  ! N + NO -> N2 + O
        TR2 = k_n_o2 * M_O2                  ! N + O2 -> NO + O
        QNODIS = 2.0_r8 * J_NO * ( TR1 / (TR1 + TR2) )

        !// HNO3 to NOx
        HNO3toNOX = &
             J_HNO3                         & ! HNO3 + hv -> NO2 + OH
             + k_oh_hno3 * M_OH              & ! HNO3 + OH -> NO3 + H2O
             + k_op_hno3 * M_O3P            ! OP + HNO3 -> OH + NO3 (very slow)

        !// NOx to HNO3
        RNOXtoHNO3 = &
             k_oh_no2_m * M_OH * M_NO2       & ! OH + NO2 + M -> HNO3 + M
             + k_no_ho2_b * M_HO2 * M_NO      & ! NO + HO2     -> HNO3
             + k_h2o_clono2 * M_H2O * M_ClONO2 & ! H2O + ClONO2 -> HNO3 + HOCl
             + k_hcl_clono2 * M_HCl * M_ClONO2 & ! HCl + ClONO2 -> 2Cl + HNO3
             !// Heterogeneous reactions on SAD
             + 2._r8 * LC_spsADPARBK * M_N2O5 &! N2O5 + H2O(sad)   -> 2HNO3
             + LC_spsACPARBK * M_N2O5       & ! N2O5 + HCl(sad)   -> ClONO + HNO3
             + LC_spsBDPARBK * M_ClONO2     & ! ClONO2 + H2O(sad) -> HOCl + HNO3
             + LC_spsBCPARBK * M_ClONO2     & ! ClONO2 + HCl(sad) -> Cl2 + HNO3
             + LC_spsFDPARBK * M_BrONO2     & ! BrONO2 + H2O(sad) -> HOCl + HNO3
             + LC_spsFCPARBK * M_BrONO2       ! BrONO2 + HCl(sad) -> BrCl + HNO3

        !// SOLID HNO3 PRODUCTION  (BOTH PSC1 AND PSC2)  ***    SPS /02/04/20
        PHNO3S = &
             2._r8 * LC_spsADPAR * M_N2O5  & ! N2O5 + H2O   -> 2HNO3s
             + LC_spsACPAR * M_N2O5       & ! N2O5 + HCl   -> ClONO + HNO3s
             + LC_spsBDPAR * M_ClONO2     & ! ClONO2 + H2O -> HOCl + HNO3s
             + LC_spsBCPAR * M_ClONO2     & ! ClONO2 + HCl -> Cl2 + HNO3s
             + LC_spsFDPAR * M_BrONO2     & ! BrONO2 + H2O -> HOBr + HNO3s
             + LC_spsFCPAR * M_BrONO2     & ! BrONO2 + HCl -> BrCl + HNO3s
             + LC_spsBHPAR * M_ClONO2       ! ClONO2 + HBr -> BrCl + HNO3s

        !// NOx PRODUCTION IS A SUM OF NOy PRODUCTION (above) AND HNO3toNOX
        PNOx = PNOy + HNO3toNOX * M_HNO3

        !// NOx DESTRUCTION IS A SUM OF NOy DESTRUCTION, RNOXtoHNO3 AND PHNO3S
        if (M_NOx_str .gt. 1.e-21_r8) then
           QNOx = ( QNODIS*M_NO + RNOXtoHNO3 + PHNO3S ) / M_NOx_str
        else
           !// Previously, an unknown bug, most likely in the psps_psc.f,
           !// sometimes led to NOx=NAN, so we check for this.
           if (M_NOx_str .ne. M_NOx_str) then
              write(*,'(A25,E10.3)') 'WARNING: NOx is not NOx: ',M_NOx_str
              write(*,'(A23)') 'Trying to set to 1.d-21'
           end if
           !// NOx is smaller than 1.d-21 (or NaN), so we set it to this value
           !// and the loss to 0.
           M_NOx_str = 1.e-21_r8
           QNOx      = 0._r8
        end if

        !// GASEOUS HNO3 PRODUCTION AND DESTRUCTION
        PHNO3 = RNOXtoHNO3 &
                + EMISX(4,L) !// EMIXS of HNO3
        QHNO3 = HNO3toNOX


        !// HNO3s is now considered to be outside the NOy family,
        !// so we should add PHNO3S to NOy loss
        if (M_NOy_str .gt. 1.e-21_r8) then
           qnoy = ( QNODIS * M_NO + PHNO3S ) / M_NOy_str
        else
           !// In case of very strong conversion of NOx into solid HNO3,
           !// NOy may eventually become small. It may be the calculation of
           !// VOLAA in the microphysics that causes unrealisticly low NOy.
           !// If this is encountered, tryto rescue with a lower limit of NOy
           write(*,'(A49)') 'WARNING: M_NOy_str < 1.d-21, repairing...'
           write(*,'(3I3,3E12.5,2I2)')ICOL,JCOL,L,M_NOy_str,qnoy,pnoy,NST
           M_NOy_str = 1.e-21_r8
           qnoy      = 0._r8
        end if

        if (LDEBUG_STRCHEM) then
           if (qnoy .ne. qnoy) then
              print*,'before',NST,M_NOy_str,pnoy,qnoy,QNODIS,M_NO,PHNO3S
              stop
           end if
        end if

        !// Save NOy
        old_v1 = M_NOy_str

        !// Integrate NOy
        call QSSA(202,'strat',DTS,EULER,STEADYST,pnoy,qnoy,M_NOy_str)

        !// Integrate NOx:
        call QSSA(203,'strat',DTS,EULER,STEADYST,pnox,qnox,M_NOx_str)

        !// Integrate GASEOUS HNO3
        call QSSA(204,'strat',DTS,EULER,STEADYST,phno3,qhno3,M_HNO3)

        !// Integrate SOLID HNO3
        !//call QSSA(205,'strat',DTS,EULER,STEADYST,phno3s,0.d0,M_HNO3s)
        !// HNO3s should not be set from production term only, but from
        !// the loss of NOy. This is because PHNO3s may be too large, producing
        !// more HNO3s than there are NOx available.
        !// NOy lost NO to NODIS and NOx to HNO3 and NOx to HNO3s and HNO3 to NOx
        !// NOx lost NO to NODIS and NOx to HNO3 and NOx to HNO3s
        !// HNO3 lost to NOx
        !// Calculate NOy using NO-dissociation as loss, i.e. no loss to HNO3s.
        !// Old NOy set before NOy integration
        if (old_v1 .gt. 1.e-21_r8) then
           qnoy = ( QNODIS * M_NO ) / M_NOy_str
        else
           old_v1=1.e-21_r8
           qnoy = 0._r8
        end if
        call QSSA(205,'strat',DTS,EULER,STEADYST,pnoy,qnoy,old_v1)
        !// The difference old_v1-M_NOy_str = produced HNO3s.
        M_HNO3s = M_HNO3s + max(0._r8,old_v1 - M_NOy_str)



        !// Integrate N2O
        PROD = 0._r8
        LOSS = k_od_n2o_a * M_O1D       & ! O(1D) + N2O -> N2 + O2
               + k_od_n2o_b * M_O1D     & ! O(1D) + N2O -> 2NO
               + J_N2O                  ! N2O + hv    -> N2 + O(1D)
        !// Diagnose loss processes [molecules/cm3 in this time step]
        CHEMLOSS(1,107,L) = CHEMLOSS(1,107,L) + LOSS*M_N2O*DTS
        CHEMLOSS(2,107,L) = 0._r8 !// No drydep
        CHEMLOSS(3,107,L) = CHEMLOSS(3,107,L) + k_od_n2o_a*M_O1D*M_N2O*DTS
        CHEMLOSS(4,107,L) = CHEMLOSS(4,107,L) + k_od_n2o_b*M_O1D*M_N2O*DTS
        CHEMLOSS(5,107,L) = CHEMLOSS(5,107,L) + J_N2O*M_N2O*DTS

        call QSSA(206,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_N2O)

        !// Integrate CH4
        PROD = 0._r8
        LOSS = ( k_od_ch4_a             & ! O(1D) + CH4 -> OH + CH3
                 + k_od_ch4_b   & ! O(1D) + CH4 -> (CH3O/CH2OH + H)  -O2-> CH2O + HO2 + H
                 + k_od_ch4_c ) * M_O1D & ! O(1D) + CH4 -> H2 + CH2O
               + k_oh_ch4 * M_OH      & ! OH + CH4    -> CH3 + H2O
               + k_cl_ch4 * M_Cl        ! Cl + CH4    -> HCl + CH3

        !// Diagnose loss processes [molecules/cm3 in this time step]
        CHEMLOSS(1,46,L) = CHEMLOSS(1,46,L) + LOSS * M_CH4*DTS
        CHEMLOSS(2,46,L) = 0._r8 !// No drydep
        CHEMLOSS(3,46,L) = CHEMLOSS(3,46,L) + k_oh_ch4*M_OH * M_CH4*DTS
        CHEMLOSS(4,46,L) = CHEMLOSS(4,46,L) + k_od_ch4_a*M_O1D * M_CH4*DTS
        CHEMLOSS(5,46,L) = CHEMLOSS(5,46,L) + k_od_ch4_b*M_O1D * M_CH4*DTS
        CHEMLOSS(6,46,L) = CHEMLOSS(6,46,L) + k_od_ch4_c*M_O1D * M_CH4*DTS
        CHEMLOSS(7,46,L) = CHEMLOSS(7,46,L) + k_cl_ch4*M_Cl * M_CH4*DTS

        call QSSA(207,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH4)

        !// Integrate CO
        PROD = ( k_oh_ch2o * M_OH         & ! OH + CH2O    -> H2O + HCO
                 + J_CH2O_a               & ! CH2O + hv    -> H + CHO
                 + J_CH2O_b               & ! CH2O + hv    -> H2 + CO
                 + k_cl_ch2o * M_Cl       & ! Cl + CH2O    -> HCl + HCO
                 + k_op_ch2o * M_O3P      & ! O(3P) + CH2O -> prod.
               ) * M_CH2O
        LOSS = k_clo_co * M_ClO           & ! ClO + CO  -> Cl + CO2
               + k_oh_co_b * M_OH         & ! OH + CO   -> CO2 + H
               + k_oh_co_a * M_OH           ! OH + CO   -> HOCO -O2-> HO2 + CO2

        call QSSA(208,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CO)


        !// F11, F12, CCl4, CH3Cl, F22, MCF, CFC113, CFC114, CFC115
        !// As the concentrations (and loss rates) of source gases for
        !// Clx are both needed, form first the loss rates for them,
        !// do then Clx and perform finally the integration of the
        !// source gas concentrations.

        !// F11, F12, CCl4, CH3Cl, F22, MCF, CFC113, CFC114, CFC115
        !// Sum up Cly production from the loss terms of CFCs.
        PROD_Cly = 0._r8

        !// Integrate CFC11
        PROD = 0._r8
        LOSS = k_od_cfc11_a * M_O1D         & ! O(1D) + CFC11 -> prod.
               + J_CFC11                      ! CFC11 + hv    -> prod.
        !// 3 Cl coming from CFC11
        PROD_Cly = PROD_Cly + 3._r8 * LOSS * M_CFC11

        call QSSA(213,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CFC11)


        !// Integrate CFC12
        PROD = 0._r8
        LOSS = k_od_cfc11_b * M_O1D         & ! O(1D) + CFC12(CCl2X) -> prod.
               + J_CFC12                      ! CFC12 + hv           -> prod.
        !// 2 Cl coming from CFC12
        PROD_Cly = PROD_Cly + 2._r8 * LOSS * M_CFC12

        call QSSA(214,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CFC12)


        !// Integrate CCl4
        PROD = 0._r8
        LOSS = J_CCl4                        ! CC4 + hv   -> prod.
        !// 4 Cl coming from CCl4
        PROD_Cly = PROD_Cly + 4._r8 * LOSS * M_CCl4

        call QSSA(215,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CCl4)


        !// Integrate CH3Cl
        PROD = 0._r8
        LOSS = k_oh_ch3cl * M_OH           & ! OH + CH3Cl -> CH2Cl + H2O
               + J_CH3Cl                    ! CH3Cl + hv -> prod.
        !// 1 Cl coming from CH3Cl
        PROD_Cly = PROD_Cly + LOSS * M_CH3Cl

        call QSSA(216,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH3Cl)


        !// Integrate HCFC22
        PROD = 0._r8
        LOSS = J_HCFC22                     & ! CHF2Cl + hv -> prod.
               + k_oh_chclf2 * M_OH           ! OH + CHF2Cl -> CF2Cl + H2O
        !// 1 Cl coming from HCFC22
        PROD_Cly = PROD_Cly + LOSS * M_HCFC22

        call QSSA(217,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_HCFC22)


        !// Integrate MCF
        LOSS = 0._r8
        LOSS = J_MCF                      & ! MCF + hv    -> prod.
               + k_oh_ch3ccl3 * M_OH        ! OH + CH3CCl3 -> CH2CCl3 + H2O
        !// 3 Cl coming from MCF
        PROD_Cly = PROD_Cly + 3._r8 * LOSS * M_MCF

        call QSSA(218,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_MCF)


        !// Integrate CFC113
        PROD = 0._r8
        LOSS = J_CFC113                       ! CF2ClCFCl2 + hv -> prod.
        !// 3 Cl coming from CFC113
        PROD_Cly = PROD_Cly + 3._r8 * LOSS * M_CFC113

        call QSSA(219,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CFC113)


        !// Integrate CFC114
        PROD = 0._r8
        LOSS = J_CFC114                       ! CF2ClCF2Cl + hv -> prod.
        !// 2 Cl coming from CFC114
        PROD_Cly = PROD_Cly + 2._r8 * LOSS * M_CFC114

        call QSSA(220,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CFC114)


        !// Integrate CFC115
        PROD = 0._r8
        LOSS = J_CFC115                       ! CF3CF2Cl + hv   -> prod.
        !// 1 Cl coming from CFC115
        PROD_Cly = PROD_Cly + LOSS * M_CFC115

        call QSSA(221,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CFC115)


        !// Integrate HCFC123
        PROD = 0._r8     
        LOSS = J_HCFC123                & ! CF3CHCl2 + hv    -> prod.
               + k_oh_hcfc123 * M_OH   & ! CF3CHCl2 + OH    -> prod.
               + k_od_hcfc123 * M_O1D    ! CF3CHCl2 + O(1D) -> prod.
        !// 2 Cl coming from HCFC123
        PROD_Cly = PROD_Cly + 2._r8 * LOSS * M_HCFC123

        !// Integrate the source gases for Clx: 
        call QSSA(222,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_HCFC123)


        !// Integrate HCFC141
        PROD = 0._r8
        LOSS = J_HCFC141                & ! CH3CFCl2 + hv    -> prod.
               + k_oh_hcfc141 * M_OH   & ! CH3CFCl2 + OH    -> prod.
               + k_od_hcfc141 * M_O1D    ! CH3CFCl2 + O(1D) -> prod.
        !// 2 Cl coming from HCFC141
        PROD_Cly = PROD_Cly + 2._r8 * LOSS * M_HCFC141

        call QSSA(223,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_HCFC141)


        !// Integrate HCFC142
        PROD = 0._r8
        LOSS = J_HCFC142                & ! CH3CF2Cl + hv    -> prod.
               + k_oh_hcfc142 * M_OH   & ! CH3CF2Cl + OH    -> prod.
               + k_od_hcfc142 * M_O1D    ! CH3CF2Cl + O(1D) -> prod.
        !// 1 Cl coming from HCFC142
        PROD_Cly = PROD_Cly + LOSS * M_HCFC141

        call QSSA(224,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_HCFC142)


        !// Integrate Cly
        PROD = PROD_Cly
        LOSS = 0._r8

        call QSSA(209,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_Cly)


        !// HCL AND Clx
        !// Integrate HCl
        PROD = ( k_cl_ch4 * M_CH4          &! Cl + CH4  -> HCl + CH3
                   + k_ho2_cl_a * M_HO2        &! Cl + HO2  -> HCl + O2
                   + k_cl_h2 * M_H2         &! Cl + H2   -> HCl + H
                   + k_cl_h2o2 * M_H2O2       &! Cl + H2O2 -> HCl + HO2
                   + k_cl_ch2o * M_CH2O )     &! Cl + H2CO -> HCl + HCO
                   * M_Cl                    &! Cl to HCl TRANSFORMATION
                 + k_oh_clo_b * M_OH * M_ClO  &! OH + ClO  -> HCl + O2
                 + fb_hv_ohcl * J_HOCL * M_HOCl  ! HOCl + hv -> HCl + O

        LOSS = k_oh_hcl * M_OH           &! HCl + OH     -> Cl + H2O
                 + k_op_hcl * M_O3P         &! HCl + O(3P)  -> OH + Cl
                 + k_hcl_clono2 * M_ClONO2 &! HCl + ClONO2 -> 2Cl + HONO2 (HNO3)
                 + J_HCL                    ! HCl + hv     -> H + Cl

        !// Heterogeneous HCl transformation
        LOSS_HCl_het = &
             + LC_spsBCPAR * M_ClONO2   &! ClONO2 + HCl(psc) -> Cl2 + HNO3s
             + LC_spsACPAR * M_N2O5     &! N2O5 + HCl(psc)   -> ClONO + HNO3s
             + LC_spsECPAR * M_HOCl     &! HOCL + HCl(psc)   -> Cl2 + H2Os
             + LC_spsFCPAR * M_BrONO2   &! BrONO2 + HCl(psc) -> BrCl + HNO3s
             + LC_spsGCPAR * M_HOBr     &! HOBr + HCl(psc)   -> BrCl + H2Os
             + LC_spsBCPARBK * M_ClONO2 &! ClONO2 + HCl(sad) -> Cl2 + HNO3
             + LC_spsACPARBK * M_N2O5   &! N2O5 + HCl(sad)   -> ClONO + HNO3
             + LC_spsECPARBK * M_HOCl   &! HOCl + HCl(sad)   -> Cl2 + H2O
             + LC_spsFCPARBK * M_BrONO2 &! BrONO2 + HCl(sad) -> BrCl + HNO3
             + LC_spsGCPARBK * M_HOBr    ! HOBr + HCl(sad)   -> BrCl + H2O

        !// In the original code, dpkclx_het was divided by HCl and added to
        !// DQKHCL, before multiplied with HCl, giving problems for small HCl.
        PROD_Clx = PROD_Cly + LOSS * M_HCl + LOSS_HCl_het

        !// Then remove LOSS_HCl_het from HCl
        if (LOSS_HCl_het .gt. 0._r8) then
           if (M_HCl .gt. 1.e-4_r8) then
              !// Skip this part for very small HCl.
              LOSS = LOSS + LOSS_HCl_het / M_HCl
           end if
        end if

        call QSSA(211,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_HCl)


        !// Integrate Clx:
        !// Loss of Clx is what is produced of HCl.
        !// Testing value of M_Clx
        if (M_Clx .gt. 1.e-21_r8) then
           LOSS = PROD / M_Clx
        else
           LOSS = 0._r8
        end if

        call QSSA(210,'strat',DTS,EULER,STEADYST,PROD_Clx,LOSS,M_Clx)


        !// Scaling of Clx and HCl with Cly:
        xCly = M_Clx + M_HCl
        M_Clx = M_Clx * M_Cly / xCly
        M_HCl = M_HCl * M_Cly / xCly


        !// Bry, MBr, H1211, H1301, H2402
        !// Define first the prod/loss rates and do the integrations
        !// after Bry has been handled too. C.f. Clx a bit earlier
        !// in the code.
        PROD_Bry = 0._r8

        !// Integrate CH3Br
        PROD = 0._r8
        LOSS = k_oh_ch3br * M_OH &! OH + CH3Br -> CH2Br + H2O
                 + J_CH3Br        ! CH3Br + hv -> CH3 + Br
        !// 1 Br in MBr 
        PROD_Bry = PROD_Bry + LOSS * M_CH3Br

        call QSSA(226,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH3Br)


        !// Integrate H1211
        PROD = 0._r8
        LOSS = J_H1211          ! CF2ClBr + hv -> prod.
        !// 1 Br in H1211
        PROD_Bry = PROD_Bry + LOSS * M_H1211

        call QSSA(227,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H1211)


        !// Integrate H1301
        PROD = 0._r8
        LOSS = J_H1301          ! CF3Br + hv   -> prod.
        !// 1 Br in H1301
        PROD_Bry = PROD_Bry + LOSS * M_H1301

        call QSSA(228,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H1301)


        !// Integrate H2402
        PROD = 0._r8
        LOSS = J_H2402         ! CF2BrCF2Br + hv -> prod.
        !// 2 Br in H2402
        PROD_Bry = PROD_Bry + 2._r8 * LOSS * M_H2402

        call QSSA(229,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H2402)


        !// Integrate Bry
        PROD = PROD_Bry
        !// The only loss process would be heterogeneous (in the
        !// troposphere), perhaps washout or aq.phase reactions?
        !// In any case, set zero
        LOSS = 0._r8

        !// Integrate Bry, MBr, H1211, H1311, H2402:
        call QSSA(225,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_Bry)

        !// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !// Calculation of the short lived species, begin:
        !// SH = H + OH + HO2 + 2*H2O2
        !// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        DPSH = ( 2._r8*k_op_ch2o * M_CH2O         &! O(3P) + CH2O -> prod.
                 + k_op_hcl * M_HCl             &! O(3P) + HCl -> OH + Cl
                 + k_op_hbr * M_HBr ) * M_O3P   &! O(3P) + HBr -> OH + Br
               + ( 2._r8*k_od_h2 * M_H2         &! O(1D) + H2 -> OH + H
                   + 2._r8*k_od_h2o * M_H2O      &! O(1D) + H2O -> OH + OH
                   + k_od_ch4_a * M_CH4          &! O(1D) + CH4 -> OH + CH3
                   + k_od_hcl * M_HCl ) * M_O1D &! O(1D) + HCl -> OH + Cl
               + k_ho2no2_heat * M_HO2NO2          &! HO2NO2 + heat -> HO2 + NO2
               + k_cl_h2 * M_H2 * M_Cl        &! H2 + Cl -> HCl + H
               + 2._r8*J_H2O * M_H2O            &! H2O + hv -> OH + H
               + J_HNO3 * M_HNO3               &! HNO3 + hv   -> NO2 + OH
               + 2._r8*J_CH2O_a * M_CH2O        &! CH2O + hv   -> H + CHO
               !//HNO2: + J_HNO2 * M_HNO2 ! HNO2 + hv ->
               + J_HCl * M_HCl                 &! HCl + hv    -> H + Cl
               + J_HOBr * M_HOBr               &! HOBr + hv   -> OH + Br
               + J_HBr * M_HBr                 &! HBr + hv    -> H + Br
               + (1._r8 - fb_hv_ohcl) * J_HOCl * M_HOCl & !HOCl + hv -> OH + Cl
               + ( J_HO2NO2_a                  &! HO2NO2 + hv -> OH + NO3
                   + J_HO2NO2_b ) * M_HO2NO2   &! HO2NO2 + hv -> HO2 + NO2
               + 2._r8*J_CH3O2H * M_CH3O2H  &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + k_no_ch3o2 * M_CH3O2 * M_NO &! NO + CH3O2 -> CH3O + NO2
               + ( k_cl_ch2o * M_Cl              &! Cl + CH2O -> HCl + HCO
                   + k_br_ch2o * M_Br ) * M_CH2O  &! Br + CH2O -> HBr + HCO
               + ( 1._r8 - fb_hv_ohcl ) * J_HOCl * ( M_OH + M_HO2 ) ! ???

        DQSH = ( k_cl_h2o2 * M_Cl               &! H2O2 + Cl -> HCl + HO2
                 + k_br_h2o2 * M_Br ) * M_H2O2      &! Br + H2O2 -> HBr + HO2
               + ( k_oh_hno3 * M_HNO3              &! OH + HNO3 -> H2O + NO3
               !//HNO2: + C1925 * M_HNO2
               !//HNO2: + C919M * M_NO
                   + k_oh_no2_m * M_NO2            &! NO2 + OH + M -> HNO3 + M
                   + k_oh_ch4 * M_CH4             &! OH + CH4    -> CH3 + H2O
                   + (k_oh_clono2_a + k_oh_clono2_b) &
                       * M_ClONO2             &! OH + ClONO2 -> Prod.
                   + k_oh_hbr*M_HBr                &! OH + HBr    -> H2O + Br
                   + k_oh_hcl * M_HCl             &! OH + HCl    -> H2O + Cl
                   + k_oh_ho2no2 * M_HO2NO2 ) * M_OH &! OH + HO2NO2 -> H2O + O2 + NO2
               + ( k_ho2_cl_a * M_Cl                  &! HO2 + Cl    -> HCl + O2
                   + k_br_ho2 * M_Br                 &! Br + HO2    -> HBr + O2
                   + k_ho2_clo * M_ClO               &! HO2 + ClO   -> HOCl + O2
                   + k_bro_ho2 * M_BrO                &! BrO + HO2   -> prod.
                   + k_ho2_ch3o2 * M_CH3O2       &! HO2 + CH3O2 -> CH3O2H + O2
                   + k_ho2_no2_m * M_NO2 ) * M_HO2     &! NO2 + HO2 + M -> HO2NO2 + M
               + ( k_oh_ch3o2h_a * M_CH3O2H &! OH + CH3OOH -> CH3O2 + H2O (70%)
                   + k_oh_hocl * M_HOCl ) * M_OH &! OH + HOCl   -> H2O + ClO
               + ( 1._r8 - fb_hv_ohcl ) * J_HOCl * ( M_OH + M_HO2 ) &! ???
               + k_oh_clo_b * M_ClO * M_OH       &! OH + ClO    -> HCl + O2
               + k_oh_bro_a * M_BrO * M_OH         &! BrO + OH    -> HBr + O2
               + k_oh_br2 * M_Br2 * M_OH          ! Br2 + OH    -> HOBr + Br

        DRSH = 2._r8*k_oh_ho2 * M_OH * M_HO2    &! OH + HO2    -> H2O + O2
               + 2._r8*k_oh_oh * M_OH * M_OH   &! OH + OH     -> H2O + O(3P)
               + 2._r8*k_oh_h2o2 * M_OH * M_H2O2  ! OH + H2O2   -> H2O + HO2

        PROD = DPSH
        LOSS = ( DQSH + DRSH ) / M_SH

        call QSSA(230,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_SH)

        !// OXYGEN
        !// OD & OP - Steady state assumption already here,  MtR, 20/9/93.

        !// Save Ox in case of possible adjustment
        tmp_YY = M_O3 + M_O3P + M_O1D

        !// Ole Amund Sovde, February 2015
        !// There are more sources and losses for O1D and O3P than was
        !// accounted for in CTM2 and first CTM3 versions; this will
        !// reduce O3P and O3 somewhat at 0.5-3hPa.
        !// M_O1D = J_O3_b * M_O3 &
        !//        / ( LC_R2M * AIR_MOLEC(L) + k_od_h2o * M_H2O )
        M_O1D = J_O3_b * M_O3 &               ! O3 + hv -> O1D
                / ( k_od_n2 * M_N2   &        ! O1D + N2 -> O3P + N2
                    + k_od_o2 * M_O2 &        ! O1D + O2 -> O3P + O2
                    + k_od_h2o * M_H2O &      ! O1D + H2O -> 2OH
                    + k_od_n2o_a * M_N2O &    ! O1D + N2O -> N2 + O2
                    + k_od_n2o_b * M_N2O &    ! O1D + N2O -> 2NO
                  )

        !// M_O3P = ( DO3 * M_O3 + J_NO2 * M_NO2 + J_NO3_b * M_NO3 ) &
        !//        / ( k_op_o2_m * M_O2 )
        M_O3P = ( J_O3_a * M_O3 &                 ! O3 + hv -> O3P + O2
                  + J_NO2 * M_NO2 &               ! NO2 + hv -> O3P + NO
                  + J_NO3_b * M_NO3 &             ! NO3 + hv -> O3P + NO2
                  + k_od_n2 * M_O1D * M_N2   &    ! O1D + N2 -> O3P + N2
                  + k_od_o2 * M_O1D * M_O2 &      ! O1D + O2 -> O3P + O2
                  + 2._r8 * J_O2 * M_O2 &         ! O2 + hv -> 2O3P
                  + k_oh_oh * M_OH * M_OH &       ! OH + OH -> OP + H2O
                  + k_h_ho2_b * M_H * M_HO2 &     ! H + HO2   -> O + H2O
                ) &
                / ( k_op_o2_m * M_O2 &         ! O3P + O2 -M-> O3
                    + k_op_ch2o * M_CH2O &     ! O3P + CH2O -> prod.
                    + k_op_ho2 * M_HO2 &       ! O3P + HO2 -> OH + O2
                    + k_op_oh * M_OH &         ! O3P + OH -> O2 + H
                    + k_op_hcl * M_HCl &       ! O3P + HCl -> OH + Cl
                    + k_op_o3 * M_O3 &         ! O3P + O3 -> 2O2
                    + k_op_clo * M_ClO &       ! O3P + ClO    -> Cl + O2
                    + k_op_clono2 * M_ClONO2 & ! O3P + ClONO2 -> ClO + NO3
                    + k_op_bro * M_BrO &       ! BrO + O(3P)  -> Br + O2
                    + k_op_hbr * M_HBr &       ! O3P + HBr -> OH + Br
                    + k_op_oclo * M_OClO &     ! OClO + O(3P) -> ClO + O2 
                  )


        if (LDEBUG_STRCHEM) then
           if (M_O3P .ne. M_O3P) then
              print*,'pchemc_str_ij.f90: O3P NaN:',M_O3, DO3
              print*,M_NO2,M_NO3,k_op_o2_m,AIR_MOLEC(L),M_O2
              print*,M_O1D,k_od_n2,k_od_o2,k_od_h2o
              stop
           end if
        end if

        !// O3P RESCUE
        !// Further down we set O3 from the sum, and if O3P is too
        !// large, O3 will be negative. We need to adjust O3P so
        !// that will not happen. I do not think we need to also set
        !// O3 when O3P is corrected; it will be set at the end
        !// of chemistry.
        !// The sum of oxygen is
        !//     [SO] = [O3] + [O3P] + [M_O1D] - [NO] - [Cl] - [Br]
        !//
        !// Assuming O3+O3P is in steady state seems to give too high O3.
        !// If O3P is found to be high above, it is better to
        !// scale down O3P. At 70km O3P is 0.1-0.2 times O3 (Brasseur
        !// and Solomon, 2005, Aeronomy of the middle atmosphere, table 5.2)

        !// Save Ox before possible adjustment
        tmp_YY = M_O3 + M_O3P + M_O1D

        !// New O3
        tmp_Y = M_SO + M_NO + M_Cl + M_Br - M_O1D - M_O3P
        tmp_O3P = M_O3P
        if (tmp_Y .lt. 1.e-2_r8 * M_O3P) then
           !// O3 will likely be negative or at least too small
           !// O3/O3P is about 1.d-6 at 100km altitude.
           !//
           !// Assume O3 to be 0.1*M_O3P and scale down M_O3P to
           !// keep SO unchanged.
           tmp_Y = M_SO + M_NO + M_Cl + M_Br - M_O1D ! old O3 + O3P
           tmp_YY = M_O3P + 0.1_r8 * M_O3P            ! Assume O3 = 0.1*O3P
           !// Scale down O3P
           tmp_O3P = M_O3P
           M_O3P = M_O3P * tmp_Y / tmp_YY

           tmp_O3 = M_O3
           M_O3 = M_SO + M_NO + M_Cl + M_Br - M_O1D - M_O3P
           !// This is the same as M_O3 = 0.1_r8 * M_O3P

           if (M_O3 .lt. 1.e-6_r8) then
              write(6,'(a,3i4,4es9.2)') 'pchemc_str_ij.f90: Neg/Tiny O3 B',&
                   L,ICOL,JCOL, M_O3P,M_O3,M_SO,M_O1D
              stop
           end if
        end if



        !// HYDROCARBONS
        !// Integrate CH3O2
        PROD = ( k_oh_ch4 * M_OH              &! OH + CH4    -> CH3 + H2O
                 + k_od_ch4_a * M_O1D           &! O(1D) + CH4 -> OH + CH3
                 + k_cl_ch4 * M_Cl ) * M_CH4  &! Cl + CH4    -> HCl + CH3
               + k_oh_ch3o2h_a * M_CH3O2H * M_OH ! OH + CH3OOH -> CH3O2 + H2O (70%)
        LOSS = k_no_ch3o2 * M_NO                &! NO + CH3O2  -> CH3O + NO2
               + k_ho2_ch3o2 * M_HO2             ! HO2 + CH3O2 -> CH3O2H + O2

        call QSSA(231,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH3O2)

        !// Integrate CH3OOH
        PROD = k_ho2_ch3o2 * M_CH3O2 * M_HO2 ! HO2 + CH3O2 -> CH3O2H + O2
        LOSS = J_CH3O2H                  &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + (k_oh_ch3o2h_a &
                  + k_oh_ch3o2h_b) * M_OH   ! OH + CH3OOH -> CH3O2/CH2OOH + H2O

        call QSSA(232,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH3O2H)

        !// Integrate CH2O
        PROD = k_no_ch3o2 * M_CH3O2 * M_NO      &! NO + CH3O2  -> CH3O + NO2
               + J_CH3O2H * M_CH3O2H          &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + ( k_od_ch4_b  &! O(1D) + CH4 -> (CH3O/CH2OH + H)  -O2-> CH2O + HO2 + H
                 + k_od_ch4_c ) * M_O1D * M_CH4 &! O(1D) + CH4 -> H2 + CH2O
               + k_oh_ch3o2h_b * M_CH3O2H * M_OH ! OH + CH3OOH -> CH2O + OH + H2O (30%)

        LOSS = J_CH2O_a               &! CH2O + hv   -> H + CHO
               + J_CH2O_b             &! CH2O + hv   -> H2 + CO
               + k_oh_ch2o * M_OH      &! OH + CH2O   -> H2O + HCO
               + k_op_ch2o * M_O3P      &! O(3P) + CH2O -> prod.
               + k_cl_ch2o * M_Cl      &! Cl + CH2O    -> HCl + HCO
               + k_br_ch2o * M_Br        ! Br + CH2O    -> HBr + HCO

        !// Diagnose loss processes [molecules/cm3 in this time step]
        CHEMLOSS(1,13,L) = CHEMLOSS(1,13,L) + LOSS * M_CH2O*DTS
        CHEMLOSS(2,13,L) = 0._r8 !// No drydep
        CHEMLOSS(3,13,L) = CHEMLOSS(3,13,L) + J_CH2O_a*M_CH2O*DTS
        CHEMLOSS(4,13,L) = CHEMLOSS(4,13,L) + J_CH2O_b*M_CH2O*DTS
        CHEMLOSS(5,13,L) = CHEMLOSS(5,13,L) + k_oh_ch2o*M_OH*M_CH2O*DTS
        CHEMLOSS(6,13,L) = CHEMLOSS(6,13,L) + k_op_ch2o*M_O3P*M_CH2O*DTS
        CHEMLOSS(7,13,L) = CHEMLOSS(7,13,L) + k_cl_ch2o*M_Cl*M_CH2O*DTS
        CHEMLOSS(8,13,L) = CHEMLOSS(8,13,L) + k_br_ch2o*M_Br*M_CH2O*DTS

        CHEMPROD(1,13,L) = CHEMPROD(1,13,L) + PROD*DTS
        CHEMPROD(2,13,L) = CHEMPROD(2,13,L) + k_no_ch3o2*M_CH3O2*M_NO*DTS
        CHEMPROD(3,13,L) = CHEMPROD(3,13,L) + J_CH3O2H*M_CH3O2H*DTS
        CHEMPROD(4,13,L) = CHEMPROD(4,13,L) + k_od_ch4_b*M_O1D*M_CH4*DTS
        CHEMPROD(5,13,L) = CHEMPROD(5,13,L) + k_od_ch4_c*M_O1D*M_CH4*DTS
        CHEMPROD(6,13,L) = CHEMPROD(6,13,L) + k_oh_ch3o2h_b*M_CH3O2H*M_OH*DTS

        call QSSA(233,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_CH2O)

        !// HYDROGEN
        !// Integrate the sum of H, OH and HO2:
        HOX = M_H + M_OH + M_HO2

        if (HOX .lt. 1.e-21_r8) then
           print *,'WARNING : HOX less than 1.e-20! at i,nst=',ICOL,JCOL,l,nst
           write(*,'(a9,3(1x,e15.9))') 'h,oh,ho2:', M_H, M_OH, M_HO2
        end if

        AH = M_H / HOX
        AOH = M_OH / HOX
        AHO2 = M_HO2 / HOX

        PHOX = k_no_ch3o2 * M_NO * M_CH3O2   &! NO + CH3O2  -> CH3O + NO2
               + 2._r8*J_CH3O2H * M_CH3O2H  &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + 2._r8*J_CH2O_a * M_CH2O    &! CH2O + hv   -> H + CHO
               + ( J_HO2NO2_a                 &! HO2NO2 + hv -> OH + NO3
                   + J_HO2NO2_b ) * M_HO2NO2  &! HO2NO2 + hv -> HO2 + NO2
               + 2._r8*J_H2O2 * M_H2O2         &! H2O2 + hv   -> 2OH
               + J_HNO3 * M_HNO3              &! HNO3 + hv   -> NO2 + OH
               !//HNO2: + J_HNO2 * M_HNO2
               + 2._r8*k_od_h2o * M_O1D * M_H2O       &! O(1D) + H2O   -> OH + OH
               + 2._r8*k_op_ch2o * M_CH2O * M_O3P      &! O(3P) + CH2O  -> prod.
               + k_ho2no2_heat * M_HO2NO2               &! HO2NO2 + heat -> HO2 + NO2
               + (1._r8 - fb_hv_ohcl) * J_HOCl * M_HOCl &!HOCl + hv -> OH + Cl
               + ( k_cl_ch2o * M_Cl &
                   + k_br_ch2o * M_Br ) * M_CH2O &
               + 2._r8*J_H2O * M_H2O          &! H2O + hv    -> H + OH
               + J_HCl * M_HCl               &! HCl + hv    -> H + Cl
               + k_cl_h2o2 * M_H2O2 * M_Cl    &! H2O2 + Cl -> HCl + HO2
               + J_HOBr * M_HOBr             &! HOBr + hv   -> OH + Br
               + J_HBr * M_HBr               &! HBr + hv    -> H + Br
               + k_br_h2o2 * M_H2O2 * M_Br     &! Br + H2O2   -> HBr + HO2
               + k_op_hcl * M_O3P * M_HCl     &! O(3P) + HCl -> OH + Cl
               + k_cl_h2 * M_Cl * M_H2      &! H2 + Cl     -> HCl + H
               + k_od_ch4_a * M_CH4 * M_O1D      &! O(1D) + CH4 -> OH + CH3
               + 2._r8*k_od_ch4_b * M_CH4 * M_O1D &! O(1D) + CH4 ->
                                                ! (CH3O/CH2OH + H) -O2-> CH2O + HO2 + H
               + k_op_hbr * M_O3P * M_HBr     &! O(3P) + HBr -> OH + Br
               + 2._r8*k_od_h2 * M_H2 * M_O1D &! O(1D) + H2  -> OH + H
               + k_od_hcl * M_O1D * M_HCl     &! O(1D) + HCl -> OH + Cl
               + (1._r8 - fb_hv_ohcl) * J_HOCl * ( M_OH + M_HO2 ) ! ???

        QHOX = ( k_ho2_ch3o2 * M_CH3O2         &! HO2 + CH3O2 -> CH3O2H + O2
                 + k_no_ho2_b * M_NO             &! NO + HO2 -> HNO3
                 + k_ho2_no2_m * M_NO2 ) * AHO2   &! NO2 + HO2 + M -> HO2NO2 + M
               + ( k_ho2_cl_a * M_Cl             &! HO2 + Cl -> HCl + O2
                   + k_br_ho2 * M_Br            &! Br + HO2 -> HBr + O2
                   + k_ho2_clo * M_ClO          &! HO2 + ClO -> HOCl + O2
                   + k_bro_ho2 * M_BrO ) * AHO2  &! BrO + HO2 -> prod.
               + ( k_oh_ch4 * M_CH4            &! OH + CH4 -> CH3 + H2O
                   + k_oh_hno3 * M_HNO3         &! OH + HNO3 -> H2O + NO3
                   + k_oh_hcl * M_HCl          &! OH + HCl -> H2O + Cl
                   + k_oh_hbr * M_HBr           &! OH + HBr -> H2O + Br
                   + k_oh_ho2no2 * M_HO2NO2     &! OH + HO2NO2 -> H2O + O2 + NO2
                   !//HNO2: + k_oh_hno2 * M_HNO2
                   !//HNO2: + k_oh_no_m * M_NO
                   + k_oh_no2_m * M_NO2 ) * AOH &! NO2 + OH + M -> HNO3 + M
               + ( k_oh_ch3o2h_a * M_CH3O2H &! OH + CH3OOH -> CH3O2 + H2O (70%)
                   + k_oh_hocl * M_HOCl         &! OH + HOCl    -> H2O + ClO
                   + (k_oh_clono2_a + k_oh_clono2_b) &
                       * M_ClONO2   &! OH + ClONO2 -> Prod.
                 ) * AOH &
               + ( 1._r8 - fb_hv_ohcl ) * J_HOCl * ( AOH + AHO2 ) &
               + ( k_oh_clo_b * M_ClO              &! OH + ClO    -> HCl + O2
                   + k_oh_bro_a * M_BrO              &! BrO + OH    -> HBr + O2
                   + k_oh_br2 * M_Br2 ) * AOH       ! Br2 + OH    -> HOBr + Br
        !// Add more when H2O chemistry is included
        if (.not. LOLD_H2OTREATMENT) then
           QHOX = QHOX &
                  + ( k_oh_hcfc123 * M_HCFC123   &! CF3CHCl2 + OH  -> H2O + ...
                      + k_oh_hcfc141 * M_HCFC141 &! CH3CFCl2 + OH  -> H2O + ...
                      + k_oh_hcfc142 * M_HCFC142 &! CH3CF2Cl + OH  -> H2O + ...
                      + k_oh_chclf2 * M_HCFC22      &! OH + CHF2Cl  -> CF2Cl + H2O
                      + k_oh_ch3ccl3 * M_MCF    &! OH + CH3CCl3 -> CH2CCl3 + H2O
                      + k_oh_ch3cl * M_CH3Cl    &! OH + CH3Cl   -> CH2Cl + H2O
                      + k_oh_ch3br * M_CH3Br ) * AOH! OH + CH3Br   -> CH2Br + H2O
        end if

        RHOX = 2._r8*k_ho2_ho2_tot * AHO2 * AHO2       &! HO2 + HO2   -> H2O2 + O2
               + 2._r8*k_oh_ho2 * AHO2 * AOH      &! OH + HO2    -> H2O + O2
               + 2._r8*k_oh_oh * AOH * AOH       &! OH +OH      -> H2O + O(3P)
               + 2._r8*k_oh_oh_m * AOH * AOH       ! OH + OH + M -> H2O2 + M

        SQ = DSQRT( QHOX*QHOX + 4._r8*RHOX*PHOX )
        HOXE= ( -QHOX + SQ ) / ( 2._r8*RHOX + 1.e-30_r8 )
        HOXU = ( -QHOX - SQ ) / ( 2._r8*RHOX + 1.e-30_r8 )
        HSIDE = ( HOX - HOXE ) / ( HOX - HOXU ) * DEXP(-SQ*DTS)
        HOXYY = ( HOXE - HOXU*HSIDE ) / ( 1._r8 - HSIDE )

        !// This is where Oslo 2-D code stops with HOX (set to HOXYY).
        !// The integration below follows the original 0-d code. MtR, 960507.
        HOXTEST = 0._r8
        if (HOXYY .lt. 0._r8) HOXTEST = 2._r8
        if (abs(1._r8 - HSIDE) .lt. 1.e-10_r8) HOXTEST = 2._r8
        if (abs(HOX - HOXU) .lt. 1.e-10_r8) HOXTEST = 2._r8
        if (HOXTEST .lt. 1._r8) HOX = HOXYY

        if (HOXTEST .gt. 1._r8) then
           call QSSA(234,'strat',DTS,EULER,STEADYST,PHOX,QHOX+RHOX*HOX,HOX)
        end if

        !// HOY is used only straight-forward integration of OH and HO2.
        !// It seems to be a balancing technique (follow e.g. C719 or
        !// R120).
        !// It is NOT used in scaling.
        HOY = M_OH + M_HO2

        !// Integrate OH
        PROD = J_CH3O2H * M_CH3O2H         &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + k_o3_ho2 * HOY * M_O3      &! O3 + HO2    -> OH + 2O2
               + k_no_ho2 * M_NO * HOY      &! NO + HO2    -> NO2 + OH
               + 2._r8*J_H2O2 * M_H2O2      &! H2O2 + hv   -> 2OH
               + J_HNO3 * M_HNO3           &! HNO3 + hv   -> NO2 + OH
               !//HNO2: + J_HNO2 * M_HNO2
               + k_op_ho2 * M_O3P * HOY           &! O(3P) + HO2  -> OH + O2
               + k_o3_h * M_O3 * M_H            &! O3 + H       -> OH + O2
               + 2._r8*k_od_h2o * M_O1D * M_H2O    &! O(1D) + H2O  -> OH + OH
               + k_op_ch2o * M_O3P * M_CH2O        &! O(3P) + CH2O -> prod.
               + (1._r8 - fb_hv_ohcl) * J_HOCl * M_HOCl &! HOCl + hv    -> Cl + OH
               + J_HOBr * M_HOBr                 &! HOBr + hv    -> Br + OH
               + k_ho2_cl_b * M_Cl * M_HO2        &! HO2 + Cl     -> OH + ClO
               + ( k_od_ch4_a * M_CH4              &! O(1D) + CH4  -> OH + CH3
                   + k_od_h2 * M_H2 ) * M_O1D    &! O(1D) + H2   -> OH + H
               + ( k_op_hcl * M_HCl               &! O(3P) + HCl  -> OH + Cl  
                   + k_op_hbr * M_HBr ) * M_O3P   &! O + HBr      -> OH + Br
               + J_HO2NO2_a * M_HO2NO2           &! HO2NO2 + hv  -> OH + NO3
               + 2._r8 * k_h_ho2_a * M_H * HOY     ! H + HO2 -> 2OH

        LOSS = k_oh_ch4 * M_CH4            &! OH + CH4  -> CH3 + H2O
               + k_oh_ch2o * M_CH2O         &! OH + CH2O -> H2O + HCO
               + k_oh_hno3 * M_HNO3         &! OH + HNO3 -> H2O + NO3
               !//HNO2: + k_oh_hno2 * M_HNO2
               + k_oh_h2o2 * M_H2O2         &! OH + H2O2 -> H2O + HO2
               + k_oh_ho2 * M_HO2          &! OH + HO2  -> H2O + O2
               + k_oh_co_b * M_CO          &! OH + CO   -> CO2 + H
               + k_o3_oh * M_O3            &! O3 + OH   -> HO2 + O2
               !//HNO2: + k_oh_no_m * M_NO
               + k_oh_no2_m * M_NO2         &! NO2 + OH + M -> HNO3 + M
               + k_oh_h2 * M_H2           &! OH + H2     -> H2O + H
               + k_op_oh * M_O3P           &! O(3P) + OH  -> O2 + H
               + k_oh_hcl * M_HCl          &! OH + HCl    -> H2O + Cl
               + k_oh_hbr * M_HBr           &! OH + HBr    -> H2O + Br
               + k_oh_ho2no2 * M_HO2NO2    &! OH + HO2NO2 -> H2O + O2 + NO2
               + (k_oh_clono2_a + k_oh_clono2_b) * M_ClONO2  &! OH + ClONO2 -> prod.
               + 2._r8*k_oh_oh * M_OH      &! OH + OH     -> H2O + O(3P)
               + k_oh_ch3o2h_a * M_CH3O2H &! OH + CH3OOH -> CH3O2 + H2O (70%)
               + k_oh_hocl * M_HOCl         &! OH + HOCl   -> H2O + ClO
               + ( k_oh_clo_a               &! OH + ClO    -> Cl + HO2
                   + k_oh_clo_b ) * M_ClO   &! OH + ClO    -> HCl + O2
               + k_oh_br2 * M_Br2           &! Br2 + OH    -> HOBr + Br
               + ( k_oh_bro_a                 &! BrO + OH    -> HBr + O2
                   + k_oh_bro_b ) * M_BrO     &! BrO + OH    -> prod.= Br + HO2
               + k_o3_ho2 * M_O3            &! O3 + HO2    -> OH + 2O2
               + k_no_ho2 * M_NO            &! NO + HO2    -> NO2 + OH
               + k_op_ho2 * M_O3P           &! O(3P) + HO2 -> OH + O2
               + 2._r8 * k_oh_oh_m * M_OH   &! OH + OH + M -> H2O2 + M
               + 2._r8 * k_h_ho2_a * M_H     ! H + HO2 -> 2OH
        !// Add more when H2O chemistry is included
        if (.not. LOLD_H2OTREATMENT) then
           LOSS = LOSS &
                  + k_oh_hcfc123 * M_HCFC123 &! CF3CHCl2 + OH  -> H2O + ...
                  + k_oh_hcfc141 * M_HCFC141 &! CH3CFCl2 + OH  -> H2O + ...
                  + k_oh_hcfc142 * M_HCFC142 &! CH3CF2Cl + OH  -> H2O + ...
                  + k_oh_chclf2 * M_HCFC22      &! OH + CHF2Cl  -> CF2Cl + H2O
                  + k_oh_ch3ccl3 * M_MCF      &! OH + CH3CCl3 -> CH2CCl3 + H2O
                  + k_oh_ch3cl * M_CH3Cl      &! OH + CH3Cl   -> CH2Cl + H2O
                  + k_oh_ch3br * M_CH3Br        ! OH + CH3Br   -> CH2Br + H2O
        end if
        !// In Frode's code, integration is called with OH, but the
        !// results are saved to XOH. This is implemented also for
        !// a number of other species later on in the code.
        XOH = M_OH
        call QSSA(235,'strat',DTS,EULER,STEADYST,PROD,LOSS,XOH)

        !// Integrate HO2
        PROD = k_no_ch3o2 * M_CH3O2 * M_NO   &! NO + CH3O2  -> CH3O + NO2
               + J_CH3O2H * M_CH3O2H       &! CH3O2H + hv -> CH2O + OH + H (->HO2)
               + J_CH2O_a * M_CH2O         &! CH2O + hv   -> H + CHO
               + k_oh_h2o2 * M_OH * M_H2O2  &! OH + H2O2   -> H2O +HO2
               + k_o2_h_m * M_O2 * M_H     &! O2 + H    -M-> HO2
               + k_o3_oh * M_O3 * HOY      &! O3 + OH     -> HO2 + O2
               + ( J_HO2NO2_b              &! HO2NO2 + hv   -> HO2 + NO2
                   + k_ho2no2_heat ) * M_HO2NO2 &! HO2NO2 + heat -> HO2 + NO2
               + k_oh_ch2o * M_CH2O * M_OH  &! OH + CH2O    -> H2O + HCO
               + ( k_op_ch2o * M_O3P              &! O(3P) + CH2O -> prod.
                   + k_cl_ch2o * M_Cl            &! Cl + CH2O    -> HCl + HCO
                   + k_br_ch2o * M_Br ) * M_CH2O  &! Br + CH2O -> HBr + HCO
               + ( k_oh_clo_a * M_ClO            &! OH + ClO  -> Cl + HO2
                   + k_oh_bro_b * M_BrO ) * M_OH   &! BrO + OH  -> prod.= Br + HO2
               + ( k_cl_h2o2 * M_Cl              &! H2O2 + Cl -> HCl + HO2
                   + k_br_h2o2 * M_Br ) * M_H2O2  &! Br + H2O2 -> HBr + HO2
               + k_od_ch4_b * M_CH4 * M_O1D       &! O(1D) + CH4 -> (CH3O/CH2OH + H)
                                                 !          -O2-> CH2O + HO2 + H
               + k_oh_co_a * M_OH * M_CO         ! OH + CO -> HOCO -O2-> HO2 + CO2

        LOSS = k_ho2_ch3o2 * M_CH3O2       &! HO2 + CH3O2   -> CH3O2H + O2
               + k_o3_ho2 * M_O3           &! O3 + HO2      -> OH + 2O2
               + 2._r8*k_ho2_ho2_tot * M_HO2    &! HO2 + HO2     -> H2O2 + O2
               + k_oh_ho2 * M_OH          &! OH + HO2      -> H2O + O2
               + k_no_ho2 * M_NO           &! NO + HO2      -> NO2 + OH
               + k_no_ho2_b * M_NO          &! NO + HO2      -> HNO3
               + k_ho2_no2_m * M_NO2         &! NO2 + HO2 + M -> HO2NO2 + M
               + k_op_ho2 * M_O3P          &! O(3P) + HO2   -> OH + O2
               + ( k_ho2_cl_a               &! HO2 + Cl      -> HCl + O2
                   + k_ho2_cl_b ) * M_Cl   &! HO2 + Cl      -> OH + ClO
               + k_br_ho2 * M_Br           &! Br + HO2      -> HBr + O2
               + k_ho2_clo * M_ClO         &! HO2 + ClO     -> HOCl + O2
               + k_bro_ho2 * M_BrO          &! BrO + HO2     -> prod.
               + k_o3_oh * M_O3           &! O3 + OH       -> HO2 + O2
               + ( k_h_ho2_a + k_h_ho2_b  &
                   + k_h_ho2_c ) * M_H     ! H + HO2       -> prods

        XHO2 = M_HO2
        call QSSA(236,'strat',DTS,EULER,STEADYST,PROD,LOSS,XHO2)


        !// Integrate H2O2
        PROD = k_ho2_ho2_tot * M_HO2 * M_HO2  &! HO2 + HO2 -> H2O2 + O2
               + k_oh_oh_m * M_OH * M_OH  ! OH + OH + M -> H2O2 + M
        LOSS = J_H2O2                    &! H2O2 + hv -> 2OH
               + k_cl_h2o2 * M_Cl         &! H2O2 + Cl -> HCl + HO2
               + k_br_h2o2 * M_Br          &! Br + H2O2 -> HBr + HO2
               + k_oh_h2o2 * M_OH          ! OH + H2O2 -> H2O +HO2

        xh2o2 = M_H2O2
        call QSSA(237,'strat',DTS,EULER,STEADYST,PROD,LOSS,xh2o2)


        if (.not. LOLD_H2OTREATMENT) then
           !// Integrate H2
           PROD = J_CH2O_b * M_CH2O         &! CH2O + hv   -> H2 + CO
                  + k_h_ho2_c * M_H * M_HO2 &! H + HO2     -> H2 + O2
                  + k_od_ch4_c * M_O1D *M_CH4  ! O(1D) + CH4 -> H2 + CH2O

           LOSS = k_od_h2 * M_O1D    &! O(1D) + H2  -> OH + H
                  + k_oh_h2 * M_OH  &! OH + H2     -> H2O + H
                  + k_cl_h2 * M_Cl   ! H2 + Cl     -> HCl + H
           !// Diagnose loss processes [molecules/cm3 in this time step]
           CHEMLOSS(1,113,L) = CHEMLOSS(1,113,L) + LOSS * M_H2*DTS
           CHEMLOSS(2,113,L) = 0._r8 !// No drydep
           CHEMLOSS(3,113,L) = CHEMLOSS(3,113,L) + k_od_h2*M_O1D*M_H2*DTS
           CHEMLOSS(4,113,L) = CHEMLOSS(4,113,L) + k_oh_h2*M_OH*M_H2*DTS
           CHEMLOSS(5,113,L) = CHEMLOSS(5,113,L) + k_cl_h2*M_Cl*M_H2*DTS
           !// Diagnose production [molecules/cm3 in this time step]
           CHEMPROD(1,113,L) = CHEMPROD(1,113,L) + PROD*DTS
           CHEMPROD(2,113,L) = CHEMPROD(2,113,L) + J_CH2O_b *M_CH2O*DTS
           CHEMPROD(3,113,L) = CHEMPROD(3,113,L) + k_h_ho2_c*M_H*M_HO2*DTS
           CHEMPROD(4,113,L) = CHEMPROD(4,113,L) + k_od_ch4_c *M_O1D*M_CH4*DTS

           call QSSA(270,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H2)


           !// Integrate H2O
           !// Included H + HO2 -> 2OH
           !//                  -> O + H2O
           !//                  -> H2 + O2
           PROD = EMISX(114,L)                 &!// Emissions, i.e. aircraft
                  + k_oh_oh * M_OH * M_OH     &! OH + OH      -> H2O + O(3P)
                  + k_oh_ho2 * M_OH * M_HO2    &! OH + HO2     -> H2O + O2
                  + k_oh_h2 * M_OH * M_H2     &! OH + H2      -> H2O + H
                  + k_oh_h2o2 * M_OH * M_H2O2   &! OH + H2O2    -> H2O + HO2
                  + k_oh_hno3 * M_HNO3 * M_OH   &! OH + HNO3    -> H2O + NO3
                  + k_oh_hocl * M_HOCl * M_OH   &! OH + HOCl    -> H2O + ClO
                  + k_oh_hcl * M_HCl * M_OH    &! OH + HCl     -> H2O + Cl
                  + k_oh_ch2o * M_CH2O * M_OH   &! OH + CH2O    -> H2O + HCO
                  + k_oh_ch4 * M_CH4 * M_OH    &! OH + CH4     -> CH3 + H2O
                  + (k_oh_ch3o2h_a + k_oh_ch3o2h_b) &
                       * M_CH3O2H * M_OH &! OH + CH3OOH  -> CH3O2/CH2OOH + H2O
                  + k_oh_hbr  * M_HBr * M_OH    &! OH + HBr     -> H2O + Br
                  + k_oh_ho2no2 * M_HO2NO2 * M_OH  &! OH + HO2NO2  -> H2O + O2 + NO2
                  + k_oh_hcfc123 * M_HCFC123 * M_OH &! CF3CHCl2 + OH  -> H2O + ...
                  + k_oh_hcfc141 * M_HCFC141 * M_OH &! CH3CFCl2 + OH  -> H2O + ...
                  + k_oh_hcfc142 * M_HCFC142 * M_OH &! CH3CF2Cl + OH  -> H2O + ...
                  + k_oh_chclf2 * M_OH * M_HCFC22 &! OH + CHF2Cl -> CF2Cl + H2O
                  + k_oh_ch3ccl3 * M_OH * M_MCF &! OH + CH3CCl3 -> CH2CCl3 + H2O
                  + k_oh_ch3cl * M_OH * M_CH3Cl   &! OH + CH3Cl   -> CH2Cl + H2O
                  + k_oh_ch3br * M_OH * M_CH3Br   &! OH + CH3Br   -> CH2Br + H2O
                  + LC_spsECPARBK * M_HOCl       &! HOCl + HCl(sad) -> Cl2 + H2O
                  + LC_spsGCPARBK * M_HOBr       &! HOBr + HCl(sad) -> BrCl + H2O
                  + k_h_ho2_b * M_H * M_HO2       ! H + HO2      -> O + H2O

           LOSS = k_od_h2o * M_O1D             &! O(1D) + H2O  -> OH + OH
                  + J_H2O                     &! H2O + hv     -> H + OH
                  + k_h2o_clono2 * M_ClONO2   &! H2O + ClONO2 -> HNO3 + HOCl
                  !//Gaseous reactions on PSCs
                  + LC_spsADPAR * M_N2O5/M_H2O     &! N2O5 + H2O   -> 2HNO3s
                  + LC_spsBDPAR * M_ClONO2/M_H2O   &! ClONO2 + H2O -> HOCl + HNO3s
                  + LC_spsFDPAR * M_BrONO2/M_H2O    ! BrONO2 + H2O -> HOBr + HNO3s

           !// Diagnose loss processes [molecules/cm3 in this time step]
           CHEMLOSS(1,114,L) = CHEMLOSS(1,114,L) + LOSS * M_H2O*DTS
           CHEMLOSS(2,114,L) = 0._r8 !// No drydep
           CHEMLOSS(3,114,L) = CHEMLOSS(3,114,L) + k_od_h2o * M_O1D*M_H2O*DTS
           CHEMLOSS(4,114,L) = CHEMLOSS(4,114,L) + J_H2O*M_H2O*DTS
           CHEMLOSS(5,114,L) = CHEMLOSS(5,114,L) + k_h2o_clono2*M_ClONO2*M_H2O*DTS
           CHEMLOSS(6,114,L) = CHEMLOSS(6,114,L) + (LC_spsADPAR * M_N2O5 &
                + LC_spsBDPAR * M_ClONO2 + LC_spsFDPAR * M_BrONO2) * DTS

           !// Diagnose production [molecules/cm3 in this time step]
           CHEMPROD(1,114,L) = CHEMPROD(1,114,L) + PROD*DTS
           CHEMPROD(2,114,L) = CHEMPROD(2,114,L) + k_oh_ho2*M_OH*M_HO2*DTS
           CHEMPROD(3,114,L) = CHEMPROD(3,114,L) + k_oh_ch4*M_CH4*M_OH*DTS
           CHEMPROD(4,114,L) = CHEMPROD(4,114,L) + k_oh_hcl*M_HCl*M_OH*DTS
           CHEMPROD(5,114,L) = CHEMPROD(5,114,L) + k_oh_hno3*M_HNO3*M_OH*DTS
           CHEMPROD(6,114,L) = CHEMPROD(6,114,L) + k_oh_h2*M_OH*M_H2*DTS
           CHEMPROD(7,114,L) = CHEMPROD(7,114,L) + k_oh_ch2o*M_CH2O*M_OH*DTS
           CHEMPROD(8,114,L) = CHEMPROD(8,114,L) + k_oh_ho2no2*M_HO2NO2*M_OH*DTS
           CHEMPROD(9,114,L) = CHEMPROD(9,114,L) + (k_oh_oh*M_OH*M_OH &
                                                   + k_oh_h2o2*M_OH*M_H2O2)*DTS
           CHEMPROD(10,114,L) = CHEMPROD(10,114,L) &
                + (k_oh_ch3o2h_a + k_oh_ch3o2h_b)*M_CH3O2H*M_OH*DTS

           call QSSA(271,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H2O)

           !// Integrate H2Os
           PROD = LC_spsECPAR * M_HOCl   &! HOCL + HCl(psc)   -> Cl2 + H2Os
                  + LC_spsGCPAR * M_HOBr  ! HOBr + HCl(psc)   -> BrCl + H2Os
           LOSS = &
                !//Reactions on aerosols (HCl/H2O are included in the uptake coeff.)
                + LC_spsADPARBK * M_N2O5/M_H2O   &! N2O5 + H2O(sad)   -> 2HNO3
                + LC_spsBDPARBK * M_ClONO2/M_H2O &! ClONO2 + H2O(sad) -> HOCl + HNO3
                + LC_spsFDPARBK * M_BrONO2/M_H2O  ! BrONO2 + H2O(sad) -> HOCl + HNO3

           call QSSA(272,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H2Os)

        end if !// if (.not. LOLD_H2OTREATMENT) then

        !// Integrate H2O from aircraft
        PROD = EMISX(148,L)
        LOSS = 0._r8
        call QSSA(148,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_H2O_ac)


        !// NITROGEN
        PROD = k_o3_no2 * M_O3 * M_NO2      &! O3 + NO2        -> NO3 + O2
               + k_n2o5_heat * M_N2O5          &! N2O5 + heat     -> NO2 + NO3 
               + k_op_no2_m * M_O3P * M_NO2  &! O(3P) + NO2 -M-> NO3
               + J_N2O5 * M_N2O5           &! N2O5 + hv       -> NO2 + NO3
               + ( k_op_clono2 * M_O3P        &! O(3P) + ClONO2  -> ClO + NO3
                   + J_CLONO2              &! ClONO2 + hv     -> NO3 + Cl
                   + k_cl_clono2 * M_Cl      &! Cl + ClONO2 -> Cl2 + NO3
                 ) * M_ClONO2              &
               + k_oh_clono2_a * M_OH * M_ClONO2 &! OH + ClONO2 -> HOCl + NO3
               + k_oh_hno3 * M_HNO3 * M_OH     &! HNO3 + OH -> NO3 + H2O
               + J_HO2NO2_a * M_HO2NO2        &! HO2NO2    -> OH + NO3
               + J_BRONO2 * M_BrONO2 * 0.29_r8  ! BrONO2    -> Br + NO3

        LOSS = k_no2_no3_m * M_NO2        &! NO2 + NO3 + M -> N2O5 + M
               + k_no_no3 * M_NO         &! NO + NO3      -> 2NO2
               + J_NO3_b                &! NO3 + hv      -> NO2 + O(3P)
               + J_NO3_a                &! NO3 + hv      -> NO + O2
               + k_no2_no3_b * M_NO2        ! NO2 + NO3     -> NO + NO2 + O2

        YNO3 = PROD / LOSS

        !// Apply PCE (P/Q) if loss is greater than 1.e-2,
        !// otherwise leave NO3 as it was.
        if (LOSS .ge. 1.e-2_r8) M_NO3 = YNO3
        if (LDEBUG_STRCHEM) then
           if (M_NO3 .ne. M_NO3) then
              print*,'OSLO_CHEM_STR_IJ: stop2',M_NO, M_O3P,M_NO3,M_HO2,M_O3,&
                   M_ClO,M_CH3O2, M_OCLO, M_BrO,QNODIS
              stop
           end if
        end if
        NOY = M_NO + M_NO2

        !// Integrate NO
        PROD = &
             EMISX(43,L) &!// Emissions NO e.g. aircraft and lightning
             !//HNO2: J_HNO2 * M_HNO2
             + ( k_op_no2 * M_O3P         &! O(3P) + NO2 -> NO + O2
                 + J_NO2                  &! NO2 + hv    -> NO + O(3P)
                 + k_no2_no3_b * M_NO3    &! NO2 + NO3   -> NO + NO2 + O2
               ) * NOY                    &
             + J_NO3_a * M_NO3             ! NO3 + hv    -> NO + O2

        LOSS = &
             !//HNO2: k_oh_no_m * M_OH +
             k_no_no3 * M_NO3         &! NO + NO3    -> 2NO2
             + k_no_ho2 * M_HO2       &! NO + HO2    -> NO2 + OH
             + k_no_ho2_b * M_NO      &! NO + HO2    -> HNO3
             + k_o3_no * M_O3         &! O3 + NO     -> NO2 + O2
             + k_no_ch3o2 * M_CH3O2   &! NO + CH3O2  -> CH3O + NO2
             + k_no_clo * M_ClO       &! NO + ClO    -> NO2 + Cl
             + k_no_bro * M_BrO       &! BrO + NO    -> NO2 + Br
             + k_no_oclo * M_OCLO     &! OClO + NO   -> NO2 + ClO
             + J_NO2                  &! NO2 + hv    -> NO + O(3P)
             + k_op_no2 * M_O3P       &! O(3P) + NO2 -> NO + O2
             + k_no2_no3_b * M_NO3    &! NO2 + NO3   -> NO + NO2 + O2
             + QNODIS                  ! NO dissociation

        XNO = M_NO
        call QSSA(238,'strat',DTS,EULER,STEADYST,PROD,LOSS,XNO)

        !// Integrate NO2
        PNO2 = &
             ( k_no_ch3o2 * M_CH3O2        &! NO + CH3O2  -> CH3O + NO2
               + 2._r8*k_no_no3 * M_NO3    &! NO + NO3    -> 2NO2
               + k_no_ho2 * M_HO2         &! NO + HO2    -> NO2 + OH
               + k_o3_no * M_O3          &! O3 + NO     -> NO2 + O2
               + k_no_clo * M_ClO         &! NO + ClO    -> NO2 + Cl
               + k_no_bro * M_BrO         &! BrO + NO    -> NO2 + Br
               + k_no_oclo * M_OCLO ) * M_NO  &! OClO + NO   -> NO2 + ClO
             !//HNO2: + k_oh_hno2 * M_HNO2 * M_OH
             + k_ho2no2_heat * M_HO2NO2         &! HO2NO2 +heat -> HO2 + NO2
             + k_n2o5_heat * M_N2O5             &! N2O5 + heat  -> NO2 + NO3
             + J_HNO3 * M_HNO3              &! HNO3 + hv    -> NO2 + OH
             + J_NO3_b * M_NO3              &! NO3 + hv     -> NO2 + O(3P)
             + J_HO2NO2_b * M_HO2NO2        &! HO2NO2 + hv  -> HO2 + NO2
             + J_N2O5 * M_N2O5              &! N2O5 + hv    -> NO2 + NO3
             + k_oh_ho2no2 * M_HO2NO2 * M_OH &! OH + HO2NO2  -> H2O + O2 + NO2
             + k_oh_clono2_b * M_ClONO2 * M_OH &! OH + ClONO2 -> Prod.
             + J_BrONO2 * M_BrONO2 * 0.71_r8  ! BrONO2 + hv  -> BrO + NO2

        !// Reaction AC : N2O5 + HCl -> 0.5 * Cl2 + NO2 + HNO3
        !// and NO2 is really NO2 + ClNO2
        PNO2 = PNO2 &
               + LC_spsACPAR * M_N2O5    &! N2O5 + HCl(psc) -> ClONO + HNO3s
               + LC_spsACPARBK * M_N2O5   ! N2O5 + HCl(sad) -> ClONO + HNO3

        QNO2 = k_o3_no2 * M_O3            &! NO2 + O3        -> NO3 + O2
               + k_no2_no3_m * M_NO3       &! NO2 + NO3 + M   -> N2O5 + M
               + k_oh_no2_m * M_OH        &! NO2 + OH + M    -> HNO3 + M
               + k_ho2_no2_m * M_HO2        &! NO2 + HO2 + M   -> HO2NO2 + M
               + k_no2_clo_m * M_ClO       &! NO2 + ClO     -M-> ClONO2
               + k_no2_bro_m * M_BrO         &! BrO + NO2 + M   -> BrONO2 + M
               + ( k_op_no2               &! O(3P) + NO2     -> NO + O2
                   + k_op_no2_m ) * M_O3P  &! O(3P) + NO2   -M-> NO3
               + J_NO2                    ! NO2 + hv        -> NO + O(3P)

        XNO2 = M_NO2
        call QSSA(239,'strat',DTS,EULER,STEADYST,pno2,qno2,XNO2)

        !// Integrate HO2NO2
        PROD = k_ho2_no2_m * M_HO2 * M_NO2 &! NO2 + HO2 + M -> HO2NO2 + M
               + EMISX(17,L)             ! EMISX of HO2NO2
               !//HNO2: + EMISX(131,L) ! EMISX of HNO2

        LOSS = k_ho2no2_heat             &! HO2NO2 + heat -> HO2 + NO2
               + J_HO2NO2_a              &! HO2NO2 + hv  -> OH + NO3
               + J_HO2NO2_b              &! HO2NO2 + hv  -> HO2 + NO2
               + k_oh_ho2no2 * M_OH       ! OH + HO2NO2   -> H2O + O2 + NO2

        XHO2NO2 = M_HO2NO2
        call QSSA(240,'strat',DTS,EULER,STEADYST,PROD,LOSS,XHO2NO2)

        !// CHLORINE
        !// The species of ClONO2, Cl, ClO and OHCl are first integrated
        !// into temporary variables "Xspecies".
        !// Integrate ClONO2
        PROD = k_no2_clo_m * M_ClO * M_NO2 ! NO2 + ClO -M-> ClONO2

        LOSS = &
             J_CLONO2              &! ClONO2 + hv -> NO3 + Cl
             + k_op_clono2 * M_O3P    &! O(3P) + ClONO2 -> ClO + NO3
             + (k_oh_clono2_a + k_oh_clono2_b) * M_OH &!OH+ClONO2 -> HOCl+NO3
             + k_cl_clono2 * M_Cl    &! Cl + ClONO2 -> Cl2 + NO3
             + k_h2o_clono2 * M_H2O   &! H2O + ClONO2 -> prod.  <2.0E-21
             + k_hcl_clono2 * M_HCl   &! HCl + ClONO2 -> 2Cl + HNO3
             + LC_spsBDPAR         &! ClONO2 + H2O(psc) -> HOCl + HNO3s
             + LC_spsBDPARBK       &! ClONO2 + H2O(SAD) -> HOCl + HNO3
             + LC_spsBCPAR         &! ClONO2 + HCl(psc) -> Cl2 + HNO3s
             + LC_spsBCPARBK       &! ClONO2 + HCl(SAD) -> Cl2 + HNO3
             + LC_spsBHPAR          ! ClONO2 + HBr(psc) -> BrCl + HNO3s

        XCLONO2 = M_ClONO2
        call QSSA(241,'strat',DTS,EULER,STEADYST,PROD,LOSS,XCLONO2)

        !// Integrate Cl
        PROD = &
             k_clo_co * M_ClO * M_CO     &! ClO + CO -> Cl + CO2
             + k_oh_hcl * M_OH * M_HCl   &! OH + HCl -> H2O + Cl
             + k_no_clo * M_NO * M_ClO    &! NO + ClO -> NO2 + Cl
             + k_op_hcl * M_O3P * M_HCl   &! O(3P) + HCl -> OH + Cl
             + k_op_clo * M_O3P * M_ClO   &! O(3P) + ClO -> Cl + O2
             + (1._r8 - fb_hv_ohcl) * J_HOCL * M_HOCl &! HOCl + hv -> OH + Cl
             + k_od_hcl * M_O1D * M_HCl   &! O(1D) + HCl -> OH + Cl
             + J_HCl * M_HCl             &! HCl +hv     -> H + Cl
             + J_CLO * M_ClO             &! ClO + hv    -> Cl + O
             + J_CLONO2 * M_ClONO2       &! ClONO2 + hv -> Cl + NO3
             + 2._r8*k_hcl_clono2 * M_HCl    &! HCl + ClONO2 -> prod. <1.0E-20
                     * M_ClONO2 &
             + k_oh_clo_a * M_ClO * M_OH  &! OH + ClO    -> Cl + HO2
             + 2._r8 * J_CL2 * M_Cl2      &! Cl2 + hv    -> 2Cl
             !// Assuming instant thermal decomp. of ClOO, omitting
             !//   + k_cloo_heat * M_ClOO
             !// so that the expression becomes 2._r8 * J_Cl2O2 * M_Cl2O2:
             + 2._r8 * J_CL2O2 * M_Cl2O2  &! Cl2O2 + hv  -> ClOO + Cl
             + J_BrCl * M_BrCl           &! BrCl + hv   -> Br + Cl
             !// This is also done for the reaction BrO + ClO -> ClOO + Br:
             + k_bro_clo_c * M_BrO * M_ClO    ! due to ClOO -> Cl + O2

        LOSS = &
             k_ho2_cl_a * M_HO2         &! HO2 + Cl    -> HCl + O2
             + k_o3_cl * M_O3         &! O3 + Cl     -> ClO + O2
             + k_cl_ch4 * M_CH4       &! Cl + CH4    -> HCl + CH3
             + k_cl_h2 * M_H2        &! H2 + Cl     -> HCl + H
             + k_cl_h2o2 * M_H2O2      &! H2O2 + Cl   -> HCl + HO2
             + k_cl_ch2o * M_CH2O      &! Cl + CH2O   -> HCl + HCO
             + k_ho2_cl_b * M_HO2      &! HO2 + Cl    -> OH + ClO
             + k_cl_clono2 * M_ClONO2    ! ClONO2 + Cl -> Cl2 + NO3
        ! And the thermal decomposition is why reaction 3504 is not
        ! included as loss here. But this is actually inconsistent
        ! with the integration of ClOO


        !// Assumption of steady state
        XCl = PROD / LOSS

        !// Integrate ClO
        PROD = &
             k_o3_cl * M_O3 * M_Cl         &! O3 + Cl        -> ClO + O2
             + k_op_clono2 * M_O3P * M_ClONO2 &! O(3P) + ClONO2 -> ClO + NO3
             + k_oh_hocl * M_OH * M_HOCl    &! OH + HOCl      -> H2O + ClO
             + k_ho2_cl_b * M_HO2 * M_Cl    &! HO2 + Cl       -> HCl + O2
             + ( k_no_oclo * M_NO            &! OClO + NO      -> NO2 + ClO
                 + k_op_oclo * M_O3P         &! OClO + O(3P)   -> ClO + O2
                 + J_OCLO ) * M_OCLO       &! OClO + hv      -> ClO + O(3P)
             + 2._r8*k_cl2o2_m * M_Cl2O2       ! Cl2O2 + M  -> ClO + ClO + M

        LOSS = &
             k_no_clo * M_NO                &! NO + ClO -> NO2 + Cl
             + k_op_clo * M_O3P             &! O(3P) + ClO -> Cl + O2
             + k_ho2_clo * M_HO2            &! HO2 + ClO -> HOCl + O2
             + k_no2_clo_m * M_NO2           &! NO2 + ClO -M-> ClONO2
             + k_clo_co * M_CO             &! ClO + CO -> Cl + CO2
             + J_ClO                       &! ClO + hv -> Cl + O
             + ( k_oh_clo_a                 &! OH + ClO -> Cl + HO2
                 + k_oh_clo_b ) * M_OH      &! OH + ClO -> HCl + O2
             + ( k_bro_clo_a                   &! BrO + ClO -> Br + OClO
                 + k_bro_clo_c                 &! BrO + ClO -> Br + ClOO
                 + k_bro_clo_b ) * M_BrO       &! BrO + ClO -> BrCl + O2
             + 2._r8*k_clo_clo_m * M_ClO         ! ClO + ClO -M-> Cl2O2

        XCLO = M_ClO
        call QSSA(242,'strat',DTS,EULER,STEADYST,PROD,LOSS,XCLO)


        !// Integrate OHCl
        PROD = &
             k_ho2_clo * M_HO2 * M_ClO        &! HO2 + ClO -> HOCl + O2
             + k_h2o_clono2 * M_H2O * M_ClONO2 &! H2O + ClONO2 -> prod.
             + k_oh_clono2_a * M_OH * M_ClONO2 &! OH + ClONO2 -> HOCl + NO3
             + LC_spsBDPAR * M_ClONO2        &! ClONO2 + H2O(psc) -> HOCl + HNO3s
             + LC_spsBDPARBK * M_ClONO2       ! ClONO2 + H2O(SAD) -> HOCl + HNO3

        LOSS = &
             J_HOCL                        &! HOCl + hv       -> OH + Cl
             + k_oh_hocl * M_OH             &! OH + HOCl       -> H2O + ClO
             + LC_spsECPAR                 &! HOCl + HCl(psc) -> Cl2 + H2Os
             + LC_spsECPARBK                ! HOCl + HCl(SAD) -> Cl2 + H2O

        XOHCl = M_HOCl
        call QSSA(243,'strat',DTS,EULER,STEADYST,PROD,LOSS,XOHCl)


        !// The species of CLZ (sum of Cl and ClO), Cl2, OClO, ClOO, Cl2O2,
        !// BrCl and Cl are integrated at once (there are some checking and
        !// scaling steps later):
        !// Integrate Clz
        ClZ = M_Cl + M_ClO

        PROD = &
             k_oh_hcl * M_OH * M_HCl              &! OH + HCl       -> H2O + Cl
             + k_op_hcl * M_O3P * M_HCl           &! O(3P) + HCl    -> OH + Cl
             + k_op_clono2 * M_O3P * M_ClONO2    &! O(3P) + ClONO2 -> ClO + NO3
             + (1._r8 - fb_hv_ohcl) * J_HOCL * M_HOCl &! HOCl + hv    -> Cl + OH
             + J_CLONO2 * M_ClONO2                 &! ClONO2 + hv  -> Cl + NO3
             + 2._r8*k_hcl_clono2 * M_HCl * M_ClONO2 &! HCl + ClONO2 -> 2Cl + HNO3
             + J_HCl * M_HCl                      &! HCl + hv     -> Cl + H
             + k_oh_hocl * M_OH * M_HOCl           &! OH + HOCl    -> H2O + ClO
             + k_od_hcl * M_O1D * M_HCl            &! O(1D) + HCl -> OH + Cl
             + 2._r8 * J_CL2 * M_Cl2                &! Cl2 + hv     -> 2Cl
             + ( k_no_oclo * M_NO                   &! NO + OClO    -> NO2 + ClO
                 + k_op_oclo * M_O3P                &! O(3P) + OClO -> ClO + O2
                 + J_OCLO ) * M_OCLO              &! OClO + hv    -> ClO + O(3P)
             !// Assuming instant thermal decomp. of ClOO, omitting
             !//   + k_cloo_heat * M_ClOO
             !// so that the expression becomes 2 * J_Cl2O2 * M_Cl2O2:
             + 2._r8 * J_CL2O2 * M_Cl2O2          &! Cl2O2 + hv  -> ClOO + Cl
             + J_BrCl * M_BrCl                    &! BrCl + hv   -> Br + Cl
             + 2._r8 * k_cl2o2_m * M_Cl2O2           ! Cl2O2 + M   -> ClO + ClO + M

        LOSS = &
             ( k_ho2_cl_a * M_HO2             &! HO2 + Cl -> ClO + OH
               + k_cl_ch4 * M_CH4           &! Cl + CH4 -> HCl + CH3
               + k_cl_h2 * M_H2            &! H2 + Cl  -> HCl + H
               + k_cl_h2o2 * M_H2O2          &! H2O2 + Cl -> HCl + HO2
               + k_cl_clono2 * M_ClONO2       &! Cl + ClONO2 -> Cl2 + NO3
               + k_cl_ch2o * M_CH2O          &! Cl + CH2O   -> HCl + HCO
             ) * M_Cl / CLZ &
             + ( k_ho2_clo * M_HO2           &! HO2 + ClO     -> HOCl + O2
                 + k_no2_clo_m * M_NO2        &! NO2 + ClO -M-> ClONO2
                 + k_oh_clo_b * M_OH         &! OH + ClO -> HCl + O2
                 + 2._r8*k_clo_clo_m * M_ClO     &! ClO + ClO -M-> Cl2O2
                 !// Assume thermal decomp. of ClOO, so we skip reaction
                 !// k_bro_clo_c where ClO+BrO->ClOO+Br; ClOO+heat->Cl+O2
                 + ( k_bro_clo_a                &! BrO + ClO -> Br + OClO
                     + k_bro_clo_b ) * M_BrO    &! BrO + ClO -> BrCl + O2
               ) * M_ClO / CLZ

        call QSSA(244,'strat',DTS,EULER,STEADYST,PROD,LOSS,Clz)

        !// Reaction AC : N2O5 + HCl -> 0.5 * Cl2 + NO2 + HNO3
        !// CL2 is really Cl2 + 2 * ClNO2


        !// Integrate Cl2
        PROD = &
             k_cl_clono2 * M_Cl * M_ClONO2   &! Cl + ClONO2 -> Cl2 + NO3
             + 0.5_r8*LC_spsACPAR * M_N2O5   &! N2O5 + HCl(psc)   -> ClONO + HNO3s
             + LC_spsBCPAR * M_ClONO2        &! CLONO2 + HCl(psc) -> Cl2 + HNO3s
             + LC_spsECPAR * M_HOCl          &! HOCL + HCl(psc)   -> Cl2 + H2Os
             + 0.5_r8*LC_spsACPARBK * M_N2O5 &! N2O5 + HCl(SAD)   -> ClONO + HNO3
             + LC_spsBCPARBK * M_ClONO2      &! CLONO2 + HCl(SAD) -> Cl2 + HNO3
             + LC_spsECPARBK * M_HOCl         ! HOCl + HCl(SAD)   -> Cl2 + H2O

        LOSS = J_CL2
        call QSSA(245,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_Cl2)


        !// Integrate OClO
        PROD = k_bro_clo_a * M_BrO * M_ClO  ! ClO + BrO    -> Br + OClO

        LOSS = J_OCLO                      &! OClO + hv    -> ClO + O(3P)
               + k_no_oclo * M_NO          &! OClO + NO    -> NO2 + ClO
               + k_op_oclo * M_O3P          ! OClO + O(3P) -> ClO + O2

        call QSSA(246,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_OClO)


        !// Integrate ClOO
        PROD = k_bro_clo_c * M_BrO * M_ClO       &! ClO + BrO   -> Br + ClOO
               + J_CL2O2 * M_Cl2O2           &! Cl2O2       -> ClOO + Cl
               + k_cl_o2_m * M_Cl * M_O2       ! Cl + O2   -M-> ClOO

        LOSS = k_cloo_heat                     ! ClOO + heat -> Cl + O2

        call QSSA(247,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_ClOO)


        !// Integrate Cl2O2
        PROD = k_clo_clo_m * M_ClO * M_ClO       ! ClO + ClO -M-> Cl2O2

        LOSS = J_CL2O2                      &! Cl2O2 + hv -> ClOO + Cl
               + k_cl2o2_m                     ! Cl2O2 + M  -> ClO + ClO + M

        call QSSA(248,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_Cl2O2)


        !// Integrate BrCl
        BrClold = M_BrCL
        PROD = &
             k_bro_clo_b * M_BrO * M_ClO    &! BrO + ClO         -> BrCl + O2
             + LC_spsFCPARBK * M_BrONO2 &! BrONO2 + HCl(SAD) -> BrCl + HNO3
             + LC_spsFCPAR * M_BrONO2   &! BrONO2 + HCl(psc) -> BrCl + HNO3s
             + LC_spsGCPARBK * M_HOBr   &! HOBr + HCl(SAD)   -> BrCl + H2O
             + LC_spsGCPAR * M_HOBr     &! HOBr + HCl(psc)   -> BrCl + H2Os
             + LC_spsBHPAR * M_ClONO2    ! ClONO2 + HBr(psc) -> BrCl + HNO3s

        LOSS = J_BrCl

        call QSSA(249,'strat',DTS,EULER,STEADYST,PROD,LOSS,M_BrCl)
        BrClx = M_BrCl

        !// Update Cl with the P and L calculated earlier (the old value
        !// for Cl (M_Cl) was used in the expressions above). Steady state Cl
        !// was found as XCl = PCL / QCL above.
        M_Cl = XCl

        !// Integrate NOZ
        !// NB: In the Oslo 2-d model, the section over NOz, N2O5 is
        !//     handled differently at layers up to about 18-20 km
        !//     and above!
        XNOZ = M_NO3 + M_N2O5
        ANO3 = M_NO3 / ( M_NO3 + M_N2O5 )

        PROD = &
             k_o3_no2 * M_O3 * M_NO2           &! O3 + NO2        -> NO3 + O2
             + k_op_no2_m * M_O3P * M_NO2       &! O(3P) + NO2   -M-> NO3
             + ( k_op_clono2 * M_O3P          &! O(3P) + ClONO2 -> ClO + NO3
                 + J_CLONO2                   &! ClONO2 + hv    -> NO3 + Cl
                 + k_cl_clono2 * M_Cl  &! Cl + ClONO2 -> Cl2 + NO3
               ) * M_ClONO2 &
             + k_oh_clono2_a * M_OH * M_ClONO2 &! OH + ClONO2 -> HOCl + NO3
             + k_oh_hno3 * M_OH * M_HNO3             &! OH + HNO3   -> H2O + NO3
             + J_HO2NO2_a * M_HO2NO2          &! HO2NO2 + hv -> OH + NO3
             + J_BrONO2 * M_BrONO2 * 0.29_r8   &! BrONO2 -> Br + NO3
             + EMISX(41,L)                    &! EMISX of NO3
             + EMISX(42,L)                     ! EMISX of N2O5

        LOSS = ( k_no_no3 * M_NO               &! NO + NO3  -> 2NO2
                 + k_no2_no3_b * M_NO2           &! NO2 + NO3 -> NO + NO2 + O2
                 + J_NO3_b                    &! NO3 + hv  -> NO2 + O(3P)
                 + J_NO3_a                    &! NO3 + hv  -> NO + O2
               ) * ANO3 &
               + LC_spsADPAR * ( 1._r8-ANO3 )   &! N2O5 + H2O(psc) -> 2HNO3s
               + LC_spsADPARBK * ( 1._r8-ANO3 ) &! N2O5 + H2O(SAD) -> 2HNO3
               + LC_spsACPAR * ( 1._r8-ANO3 )   &! N2O5 + HCl(psc) -> ClONO + HNO3s
               + LC_spsACPARBK * ( 1._r8-ANO3 )  ! N2O5 + HCl(SAD) -> ClONO + HNO3

        call QSSA(250,'strat',DTS,EULER,STEADYST,PROD,LOSS,xNOz)

        !// Integrate N2O5
        PROD = k_no2_no3_m * M_NO3 * M_NO2  &! NO2 + NO3 + M -> N2O5 + M
               + EMISX(42,L)                 ! EMISX of N2O5

        LOSS = k_n2o5_heat                  &! N2O5 + heat -> NO2 + NO3
               + J_N2O5                     &! N2O5 + hv   -> NO2 + NO3
               + LC_spsADPAR                &! N2O5 + H2O(psc) -> 2HNO3s
               + LC_spsADPARBK              &! N2O5 + H2O(SAD) -> 2HNO3
               + LC_spsACPAR                &! N2O5 + HCl(psc) -> ClONO + HNO3s
               + LC_spsACPARBK               ! N2O5 + HCl(SAD) -> ClONO + HNO3

        XN2O5 = M_N2O5
        call QSSA(251,'strat',DTS,EULER,STEADYST,PROD,LOSS,XN2O5)


        PROD = &
             k_o3_no2 * M_O3 * M_NO2      &! O3 + NO2        -> NO3 + O2
             + k_n2o5_heat * M_N2O5          &! N2O5 + heat     -> NO2 + NO3
             + k_op_no2_m * M_O3P * M_NO2  &! O(3P) + NO2 -M-> NO3
             + J_N2O5 * M_N2O5           &! N2O5 + hv       -> NO2 + NO3
             + ( k_op_clono2 * M_O3P     &! O(3P) + ClONO2  -> ClO + NO3
                 + J_CLONO2              &! ClONO2 + hv     -> NO3 + Cl
                 + k_cl_clono2 * M_Cl    &! Cl + ClONO2     -> Cl2 + NO3
               ) * M_ClONO2 &
             + k_oh_clono2_a * M_OH * M_ClONO2 &! OH + ClONO2 -> HOCl + NO3
             + k_oh_hno3 * M_HNO3 * M_OH             &! HNO3 + OH -> NO3 + H2O 
             + J_HO2NO2_a * M_HO2NO2                &! HO2NO2    -> OH + NO3
             + J_BrONO2 * M_BrONO2 * 0.29_r8         &! BrONO2    -> Br + NO3
             + EMISX(41,L)                           ! EMISX of NO3

        LOSS = &
             k_no2_no3_m * M_NO2         &! NO2 + NO3 + M -> N2O5 + M
             + k_no_no3 * M_NO          &! NO + NO3      -> 2NO2
             + J_NO3_b                 &! NO3 + hv      -> NO2 + O(3P)
             + J_NO3_a                 &! NO3 + hv      -> NO + O2
             + k_no2_no3_b * M_NO2         ! NO2 + NO3     -> NO + NO2 + O2

        !// Temporary PCE for NO3:
        XNO3 = PROD / LOSS

        !// Set smaller of NO3 and N2O5 into the individually integrated
        !// value and the other to the difference between the integrated
        !// sum of them and the individually integrated value.
        if (xNO3 .ge. xN2O5) then
           M_N2O5 = xN2O5
           M_NO3 = xNOz - xN2O5
        else
           M_N2O5 = xNOz - xNO3
           M_NO3 = xNO3
        end if

        !// If either N2O5 or NO3 acquired a negative value, undo the
        !// scaling above by scaling the individually integrated
        !// species with the integrated sum.
        SC = XNOZ / ( XN2O5 + xNO3 )
        zN2O5 = xN2O5 * SC
        zNO3 = xNO3 * SC
        X3 = min( M_N2O5,M_NO3 )
        if (X3 .lt. 0._r8) then
           M_N2O5 = zN2O5
           M_NO3  = zNO3
        end if

        !// Scaling of odd hydrogen 1-3:
        !// 1) Set the smaller of OH and HO2 to the individually
        !//    calculated value and the other of the two as the
        !//    difference from HOX-integr (HOX=OH+HO2+H).
        if (XHO2 .ge. XOH) then
           M_OH  = XOH
           M_HO2 = HOX - XOH - M_H
        else
           M_OH  = HOX - XHO2 - M_H
           M_HO2 = XHO2
        end if
        !// 2) Scale OH and HO2 if they have acquired negative values
        !//    or if the one which is greater is not large enough (5x).
        !//    Scaling is done with individually integrated values.
        if ( M_OH .lt. 0._r8 .or. M_HO2 .lt. 0._r8 .or. &
             (XHO2 .gt. XOH .and. XHO2 .lt. 5._r8*XOH) .or. &
             (XOH .gt. XHO2 .and. XOH .lt. 5._r8*XHO2)) then
           !// Test .and. should go before .or.
           if (.not. ( M_OH .lt. 0._r8 .or. M_HO2 .lt. 0._r8 .or. &
                XHO2 .gt. XOH .and. XHO2 .lt. 5._r8*XOH .or. &
                XOH .gt. XHO2 .and. XOH .lt. 5._r8*XHO2 ) ) then
              print*,'OSLO_CHEM_STR_IJ: ups for SC?'
           end if
           !// Scale HO2 and OH instead of using "old" values
           SC = HOX / ( XOH + XHO2 )
           M_HO2 = XHO2*SC
           M_OH = XOH*SC
        end if
        !// 3) Scale the individually integrated H2O2 and the others
        !//    with the integrated family amount.
        YH = 2._r8 * XH2O2 + M_HO2 + M_OH + M_H
        SC = M_SH / YH
        M_OH = M_OH * SC
        M_HO2 = M_HO2 * SC
        M_H2O2 = XH2O2 * SC
        !// H  Why only now? To use the present OH?
        PROD = &
             k_oh_co_b * M_OH * M_CO       &! OH + CO    -> CO2 + H
             + J_CH2O_a * M_CH2O           &! CH2O + hv  -> H + CHO
             + k_op_oh * M_O3P * M_OH      &! O(3P) + OH -> O2 + H
             + k_oh_h2 * M_H2 * M_OH      &! OH + H2    -> H2O + H
             + k_od_h2 * M_H2 * M_O1D      &! O(1D) + H2 -> OH + H
             + J_H2O * M_H2O               &! H2O + hv   -> OH + H
             + J_HCl * M_HCl               &! HCl + hv   -> Cl + H
             + J_HBr * M_HBr               &! HBr + hv   -> Br + H
             + k_cl_h2 * M_H2 * M_Cl      &! H2 + Cl    -> HCl + H
             + k_od_ch4_b * M_CH4 * M_O1D   ! O(1D) + CH4 -> (CH3O/CH2OH + H)
                                            !          -O2-> CH2O + HO2 + H

        LOSS = &
             k_o2_h_m * M_O2                &! O2 + H   -M-> HO2
             + k_o3_h * M_O3               &! O3 + H     -> OH + O2
             + ( k_h_ho2_a             &! H + HO2 -> prods
                 + k_h_ho2_b           &! H + HO2 -> prods
                 + k_h_ho2_c ) * M_HO2  ! H + HO2 -> prods

        !// Steady state H:
        M_H = PROD / LOSS

        !// In Oslo 2-d model, here is some "HOXMIN" stuff... i.e. if
        !// either OH or HO2 becomes less than 1.e+1, they are set to
        !// 10 and H2O2 is given the value of half of SH. Reason ?
        !// In Oslo 2-d model, there is also an additional treatment
        !// for tropospheric levels for odd hydrogen.

        !// Scaling of odd nitrogen... Note that BrONO2 is not included!
        !// Set HNO2 and ClONO2 to the individually integrated values:
        !//HNO2: M_HNO2 = XHNO2
        M_ClONO2 = XCLONO2
        !// Make new temporary N2O5, NO3 and NOZ
        XN2O5 = M_N2O5
        XNO3 = M_NO3
        VNOZ = XNO + XNO2 + XNO3 + 2.*XN2O5 + XHO2NO2 + XCLONO2 + M_BrONO2

        !// Set here the largest of HO2NO2, NO3, N2O5, NO2, and NO with
        !// the integrated family amount and the other members to their
        !// individually integrated values. NO2 and NO are multiplied by
        !// 10. (TESTNO=TESTNO2=10.) in the comparison. WHY? MtR, 960507.
        !// Note that only one of the species can become reset with the
        !// total family.
        !// Note also that NO3 and N2O5 have already been scaled with
        !// their new sum and that BrONO2 is not taken into account
        !// herein. WHY? Is it so small as to be ignored?
        if (XHO2NO2 .gt. TESTNO2*XNO2 .and. XHO2NO2 .gt. TESTNO*XNO .and. &
             XHO2NO2 .gt. 2._r8*XN2O5 .and. XHO2NO2 .gt. XNO3) then
           M_HO2NO2 = XHO2NO2 + M_NOx_str - VNOZ
        else
           M_HO2NO2 = XHO2NO2
        end if

        if (XNO3 .gt. TESTNO2*XNO2 .and. XNO3 .gt. TESTNO*XNO .and. &
         XNO3 .gt. 2._r8*XN2O5 .and. XNO3 .gt. XHO2NO2) then
           M_NO3 = XNO3 + M_NOx_str - VNOZ
        else
           M_NO3 = XNO3
        end if

        if (2._r8*XN2O5 .gt. TESTNO2*XNO2 .and. 2._r8*XN2O5 .gt. TESTNO*XNO .and. &
             2._r8*XN2O5 .gt. XNO3 .and. 2._r8*XN2O5 .gt. XHO2NO2) then
           M_N2O5 = XN2O5 + 0.5_r8*(M_NOx_str - VNOZ)
        else
           M_N2O5 = XN2O5
        end if

        if (TESTNO2*XNO2 .gt. XHO2NO2 .and. TESTNO2*XNO2 .gt. TESTNO*XNO .and. &
             TESTNO2*XNO2 .gt. 2._r8*XN2O5 .and. TESTNO2*XNO2 .gt. XNO3) then
           M_NO2 = XNO2 + M_NOx_str - VNOZ
        else
           M_NO2 = XNO2
        end if

        if (TESTNO*XNO .gt. XHO2NO2 .and. TESTNO*XNO .gt. TESTNO2*XNO2 .and. &
             TESTNO*XNO .gt. 2._r8*XN2O5 .and. TESTNO*XNO .gt. XNO3) then
           M_NO = XNO + M_NOx_str - VNOZ
        else
           M_NO = XNO
        end if

        !// If any of the NOx-members just treated acquired a negative
        !// value, "undo" the scaling and use the (integrated family)
        !// scaled individually integrated value.
        TEST = min(M_NO2,M_NO3,M_N2O5,M_NO,M_HO2NO2)

        if (TEST .lt. 0._r8) then
           FAC = M_NOx_str/VNOZ       
           M_NO = XNO*FAC
           M_NO2 = XNO2*FAC
           M_NO3 = XNO3*FAC
           M_N2O5 = XN2O5*FAC
           M_HO2NO2 = XHO2NO2*FAC
           !//HNO2: M_HNO2 = XHNO2*FAC
           M_ClONO2 = XCLONO2*FAC

           !// and it seems as Cl and ClO are mixed in here, reason?
           M_Cl = M_Cl + XClONO2*(1._r8 - FAC)*M_Cl/ClZ
           ClZ = ClZ + XClONO2*(1._r8 - FAC)
        end if

        !// But wouldn't this be better than the stuff with YClO below?
        M_HOCl = XOHCl
        ClONO2G = M_ClONO2

        !// Find ClO (can get negative, will check it in Y3 below)
        YClO = M_ClX - XClONO2 - M_HOCl - M_Cl - M_OClO &
               - 2._r8*M_Cl2 - M_ClOO - 2._r8*M_Cl2O2 - M_BrCl

        M_ClO = YClO

        !// If either ClONO2 or ClO has acquired a negative value:
        !// reset all members of Clx to (integrated sum) scaled,
        !// individually integrated values. Otherwise keep the values
        !// set so far.
        Y3 = min(M_ClONO2,M_ClO)

        !// No need to use CVMGP to test all species.
        if (Y3 .lt. 0._r8) then
           SC = M_ClX &
                / ( XCLONO2 + XOHCl + XCLO + XCL + M_BrCL + M_OCLO &
                    + 2._r8*M_Cl2 + M_ClOO + 2._r8*M_Cl2O2 )

           M_ClONO2 = XClONO2 * SC
           M_HOCl = XOHCl * SC
           M_ClO = XClO * SC
           M_Cl = XCl * SC
           M_Cl2 = M_Cl2 * SC
           M_OClO = M_OClO * SC
           M_Cl2O2 = M_Cl2O2 * SC
           M_ClOO = M_ClOO * SC
           M_BrCl = M_BrCl * SC    ! Incl: MtR, 950505
        end if

        !//****************************************************************
        !// Special case: ClONO2 exceeds available NOx
        !// In Oslo 2-d "TESTFRAC" is 0.7 
        !//****************************************************************
        TESTFRAC=0.97_r8
        if (M_ClONO2 .gt. TESTFRAC*M_NOx_str) then

           M_ClONO2 = TESTFRAC * M_NOx_str

           YYClX = M_ClX - M_ClONO2

           !// In Oslo 2-D this is done with "Z-species", both in the
           !// divisor to form SC and in applying SC below.
           SC = YYCLX &
                / ( M_HOCl + M_ClO + M_Cl + M_BrCL + 2._r8*M_Cl2 &
                    + M_OCLO + M_ClOO + 2._r8*M_Cl2O2 )

           M_HOCl = M_HOCl * SC
           M_ClO = M_ClO * SC
           M_Cl = M_Cl * SC
           M_Cl2 = M_Cl2 * SC
           M_OClO = M_OClO * SC
           M_Cl2O2 = M_Cl2O2 * SC
           M_ClOO = M_ClOO * SC
           M_BrCl = M_BrCl * SC
        end if

        !// DOES NOT TAKE INTO ACCOUNT BrONO2, WHY ?
        YNOXNEW = M_NOx_str - M_ClONO2
        YNOXOLD = M_NOx_str - ClONO2G
        SCX = YNOXNEW / YNOXOLD
        ZNO = M_NO * SCX
        ZNO2 = M_NO2 * SCX
        ZNO3 = M_NO3 * SCX
        ZN2O5 = M_N2O5 * SCX
        !//HNO2: ZHNO2 = M_HNO2 * SCX
        ZHO2NO2 = M_HO2NO2 * SCX

        !// Correct Nitrogen
        X3 = Y3 ! Y3 was set earlier
        if (X3 .lt. 0._r8) then
           if (LDEBUG_STRCHEM) then
              if (ZNO .ne. ZNO) then
                 print*,'pchemc_str_ij.f: NO NaN',ZNO,M_NO,X3
                 stop
              end if
           end if
           M_NO = ZNO
           M_NO2 = ZNO2
           M_NO3 = ZNO3
           M_N2O5 = ZN2O5
           !//HNO2: M_HNO2 = ZHNO2
           M_HO2NO2 = ZHO2NO2
        end if

        !// BROMINE (EXCEPT BrCl which was already done with Chlorine).
        !// Br2, OHBr, BrNO3, HBr, Brz, Br, BrO

        PBr2 = k_bro_bro_b * M_BrO * M_BrO       ! BrO + BrO -> Br2 + O2

        QBr2 = J_Br2                        &! Br2 + hv -> 2Br
               + k_oh_br2 * M_OH              ! Br2 + OH -> HOBr + Br

        POHBr = k_bro_ho2 * M_BrO * M_HO2    &! BrO + HO2 -> prod.
                + k_oh_br2 * M_Br2 * M_OH   &! Br2 + OH -> HOBr + Br
                + LC_spsFDPARBK * M_BrONO2 &! BrONO2 + H2O(SAD) -> HOBr + HNO3
                + LC_spsFDPAR * M_BrONO2    ! BrONO2 + H2O(psc) -> HOBr + HNO3s


        QOHBr = J_HOBr                     &! HOBr + hv   -> OH + Br
                + LC_spsGCPARBK            &! HOBr + HCl(SAD) -> BrCl + H2O
                + LC_spsGCPAR               ! HOBr + HCl(psc) -> BrCl + H2Os

        PBrNO3 = k_no2_bro_m * M_NO2           &! BrO + NO2 + M -> BrONO2 + M
              * ( M_BrO + M_BrONO2 )

        QBrNO3 = &
             J_BrONO2                   &! BrONO2 + hv   -> prod.
             + k_no2_bro_m * M_NO2          &! BrO + NO2 + M -> BrONO2 + M
             + LC_spsFDPARBK            &! BrONO2 + H2O(SAD) -> HOBr + HNO3
             + LC_spsFCPARBK            &! BrONO2 + HCl(SAD) -> BrCl + HNO3
             + LC_spsFDPAR              &! BrONO2 + H2O(psc) -> HOBr + HNO3s
             + LC_spsFCPAR               ! BrONO2 + HCl(psc) -> BrCl + HNO3s


        PHBr = &
             ( k_br_h2o2 * M_H2O2         &! Br + H2O2 -> HBr + HO2
               + k_br_ch2o * M_CH2O       &! Br + CH2O -> HBr + HCO
               + k_br_ho2 * M_HO2        &! Br + HO2  -> HBr + O2
             ) * M_Br &
             + k_oh_bro_a * M_BrO * M_OH    ! BrO + OH  -> HBr + O2

        QHBr = &
             J_HBr                      &! HBr + hv    -> H + Br
             + k_oh_hbr * M_OH           &! OH + HBr    -> H2O + Br
             + k_op_hbr * M_O3P           ! O(3P) + HBr -> OH + Br
        if (M_HBr .gt. 1.e-8_r8) then
           QHBr = QHBr &
                  + LC_spsBHPAR*M_ClONO2/M_HBR ! ClONO2 + HBr(psc) -> BrCl + HNO3s
        end if

        !// Model the sum of Br and BrO:
        BrZ = M_Br + M_BrO

        PBrZ = &
             J_BrONO2 * M_BrONO2         &! BrONO2 + hv -> prod.
             + J_HOBr * M_HOBr           &! HOBr + hv   -> OH + Br
             + 2._r8*J_Br2 * M_Br2        &! Br2 + hv    -> 2Br
             + k_oh_br2 * M_Br2 * M_OH    &! Br2 + OH    -> HOBr + Br
             + ( J_HBr                   &! HBr + hv    -> H + Br
                 + k_oh_hbr * M_OH        &! OH + HBr    -> H2O + Br
                 + k_op_hbr * M_O3P       &! O(3P) + HBr -> OH + Br
             ) * M_HBr &
             + J_BrCl * BrCLold            ! BrCl + hv   -> Br + Cl

        QBrZ = &
             ( k_br_h2o2 * M_H2O2          &! Br + H2O2     -> HBr + HO2
               + k_br_ch2o * M_CH2O        &! Br + CH2O     -> HBr + HCO
               + k_br_ho2 * M_HO2         &! Br + HO2      -> HBr + O2
             ) * M_Br / BrZ &
             + ( k_no2_bro_m * M_NO2         &! BrO + NO2 + M -> BrONO2 + M
                 + k_bro_ho2 * M_HO2       &! BrO + HO2     -> prod.
                 + 2._r8*k_bro_bro_b * M_BrO  &! BrO + BrO     -> Br2 + O2
                 + k_oh_bro_a * M_OH        &! BrO + OH      -> HBr + O2
                 + k_bro_clo_b * M_ClO       &! BrO + ClO     -> BrCl + O2
               ) * M_BrO / BrZ

        PBr = &
             ( J_BrO                     &! BrO + hv    -> Br + O
               + k_op_bro * M_O3P         &! BrO + O(3P) -> Br + O2
               + ( k_bro_clo_a               &! BrO + ClO   -> Br + OClO
                   + k_bro_clo_c ) * M_ClO   &! BrO + ClO   -> Br + ClOO
               + k_no_bro * M_NO          &! BrO + NO    -> NO2 + Br
               + 2._r8*k_bro_bro_a * M_BrO    &! BrO + BrO   -> 2Br + O2
               + k_oh_bro_b * M_OH          &! BrO + OH    -> prod. (Br + HO2)
               + k_o3_bro * M_O3          &! BrO + O3    -> Br + 2O2
             ) * M_BrO &
             + J_HOBr * M_HOBr           &! HOBr + hv   -> OH + Br
             + 2._r8*J_Br2 * M_Br2        &! Br2 + hv    -> 2Br
             + ( J_HBr                   &! HBr + hv    -> H + Br
                 + k_oh_hbr * M_OH        &! OH + HBr    -> H2O + Br
                 + k_op_hbr * M_O3P       &! O(3P) + HBr -> OH + Br
               ) * M_HBr &
             + k_oh_br2 * M_Br2 * M_OH    &! Br2 + OH    -> HOBr + Br
             + J_BrCl * BrClold          &! BrCl + hv   -> Br + Cl
             + J_BrONO2 * M_BrONO2 * 0.29_r8 ! BrONO2 + hv -> NO3 + Br

        QBr = &
             k_br_o3 * M_O3              &! Br + O3   -> BrO + O2
             + k_br_h2o2 * M_H2O2          &! Br + H2O2 -> HBr + HO2
             + k_br_ch2o * M_CH2O          &! Br + CH2O -> HBr + HCO
             + k_br_ho2 * M_HO2            ! Br + HO2  -> HBr + O2


        PBrO = &
             k_br_o3 * M_Br * M_O3         &! Br + O3     -> BrO + O2
             + J_BrONO2 * M_BrONO2 * 0.71_r8 ! BrONO2 + hv -> NO2 + BrO

        QBrO = &
             J_BrO                         &! BrO + hv      -> Br + O
             + k_op_bro * M_O3P            &! BrO + O(3P)   -> Br + O2
             + ( k_bro_clo_a               &! BrO + ClO     -> Br + OClO
                 + k_bro_clo_c                 &! BrO + ClO     -> Br + ClOO
                 + k_bro_clo_b ) * M_ClO       &! BrO + ClO     -> BrCl + O2
             + k_no_bro * M_NO              &! BrO + NO      -> NO2 + Br
             + k_no2_bro_m * M_NO2             &! BrO + NO2 + M -> BrONO2 + M
             + 2._r8 * ( k_bro_bro_a            &! BrO + BrO     -> 2Br + O2
                        + k_bro_bro_b ) * M_BrO &! BrO + BrO     -> Br2 + O2
             + k_o3_bro * M_O3              &! BrO + O3      -> Br + 2O2
             + k_bro_ho2 * M_HO2             &! BrO + HO2     -> prod.
             + ( k_oh_bro_a                   &! BrO + OH      -> HBr + O2
                 + k_oh_bro_b ) * M_OH         ! BrO + OH      -> prod. (Br + HO2)

        !// Integrations... Integrate OHBr, BrONO2, HBr, Br, BrO and
        !// BrZ into preliminary "X-species". Br2 can be set directly
        !//(as there is no code for it later).

        call QSSA(252,'strat',DTS,EULER,STEADYST,pBr2,qBr2,M_Br2)
        Br2x = M_Br2
         
        OHBrx = M_HOBr
        call QSSA(253,'strat',DTS,EULER,STEADYST,pOHBr,qOHBr,OHBrx)

        BrNO3x = M_BrONO2
        call QSSA(254,'strat',DTS,EULER,STEADYST,pBrNO3,qBrNO3,BrNO3x)

        HBrx = M_HBr
        call QSSA(255,'strat',DTS,EULER,STEADYST,pHBr,qHBr,HBrx)

        Brzx = Brz
        call QSSA(256,'strat',DTS,EULER,STEADYST,pBrz,qBrz,Brzx)

        Brx = M_Br
        call QSSA(257,'strat',DTS,EULER,STEADYST,pBr,qBr,Brx)

        BrOx = M_BrO
        call QSSA(258,'strat',DTS,EULER,STEADYST,pBrO,qBrO,BrOx)

        !// Set the largest of BrZ (Br+BrO), BrONO2, HBr, OHBr with the
        !// intgrated family amount and the others to their individually
        !// integrated values. As BrCl has a comparable conc. (?) to
        !// HBr and OHBr (at least under some conditions, it is thus
        !// included herein. Br2 seems to be excluded. MtR, 960507.
        !// Note that only one of the species underneath can become reset with
        !// Bry and BrTOT.
        BrTOT = BrZX + BrNO3X + 2._r8*Br2X + HBrX + OHBrX + M_BrCl

        if (BrZX .gt. M_BrCl .and. BrZX .gt. BrNO3X .and. &
             BrZX .gt. HBrX .and. BrZX .gt. OHBrX) then
           BrZ = BrZX + M_Bry - BrTOT
        else
           BrZ = BrZX
        end if
      
        if (BrNO3X .gt. M_BrCl .and. BrNO3X .gt. BrZX .and. &
             BrNO3X .gt. HBrX .and. BrNO3X .gt. OHBrX) then
           M_BrONO2 = BrNO3X + M_Bry - BrTOT
        else
           M_BrONO2 = BrNO3X
        end if

        if (HBrX .gt. M_BrCl .and. HBrX .gt. BrZX .and. &
             HBrX .gt. BrNO3X .and. HBrX .gt. OHBrX) then
           M_HBr = HBrX + M_Bry - BrTOT
        else
           M_HBr = HBrX
        end if

        if (OHBrX .gt. M_BrCl .and. OHBrX .gt. BrZX .and. &
             OHBrX .gt. BrNO3X .and. OHBrX .gt. HBrX) then
           M_HOBr = OHBrX + M_Bry - BrTOT
        else
           M_HOBr = OHBrX
        end if

        if (M_BrCl .gt. BrZX .and. M_BrCl .gt. BrNO3X .and. &
             M_BrCl .gt. HBrX .and. M_BrCl .gt. OHBrX) then
           M_BrCl = M_BrCl + M_Bry - BrTOT
           !// No need for else-statement
        end if


        !// If any of the Bry-members set in the routine above became
        !// negative, scale the individually integrated values with the
        !// integrated family. Br2 is again excluded from this.

        SCAL = ( M_Bry - 2._r8 * Br2x ) &
               / ( BrZX + BrNO3X + HBrX + OHBrx + BrClx )
        if (LDEBUG_STRCHEM) then
           if (SCAL .lt. 0._r8) then
              print*,'OSLO_CHEM_STR_IJ: Br SCAL',ICOL,JCOL,l,SCAL
              print*,M_Bry, Br2x, M_Br2, BrTOT
              print*,pbr2,qbr2,pbr2/qbr2
              stop
           end if
        end if

        if ( M_BrONO2 .lt. 0._r8 .or. M_HBr .lt. 0._r8 .or. &
             BrZ .lt. 0._r8 .or. &
             M_HOBr .lt. 0._r8 .or. M_BrCl .lt. 0._r8) then
           M_BrONO2 = BrNO3X*SCAL
           M_HBr = HBrX*SCAL
           BrZ = BrZX*SCAL
           M_HOBr = OHBrX*SCAL
           M_BrCl = BrClx*SCAL
           !// No need for else
        end if

        !// Largest Brz species is set with the
        !// integrated and checked sum of Br and BrO.
        !// The other of the two is set to its individually integrated value.
        if (BrX .gt. BrOX) then
           M_BrO = BrOX
           M_Br  = BrZ - BrOX
        else
           M_BrO = BrZ - BrX
           M_Br  = BrX
        end if

        !// Check that setting the largest BrZ-species did not
        !// produce negative values.
        BrTST = min(M_Br,M_BrO)
        if (BrTST .le. 0._r8) then
           M_BrO = BrZ * BrOX / (BrOX + BrX)
           M_Br  = BrZ * BrX / (BrOX + BrX)
           !// No need for else
        end if

        !// Extract O3 from SO and the other members of SO:
        M_O3 = M_SO - M_O3P - M_O1D + M_NO + M_Cl + M_Br
        if (M_O3 .lt. 1.e-6_r8) then
           write(6,'(a,3i4,4es9.2)') 'pchemc_str_ij.f90: Neg/Tiny O3 C',&
                L,ICOL,JCOL, M_O3P,M_O3,M_SO,M_O1D
           stop
        end if

        !// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !// Calculation of changes in the concentrations
        !// of short lived species, end
        !// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        !// ***************************************************************
        !//Scale the members in the chemical families 
        !//
        !// Scale again over the families, if this is the last round of
        !// chemistry for now and we go back to transport.
        !//
        !// Bry=Br+BrO+BrONO2+OHBr+HBr+2.*Br2+BrCl
        !// SO=O3+OD+OP-NO-Cl-Br
        !// NOx_str=NO+NO2+NO3+2.*N2O5+HNO2+HO2NO2+ClONO2+BrONO2
        !// Clx=Cl+ClO+OHCl+ClONO2+2.*Cl2+OClO+BrCl+ClOO+2.*Cl2O2
        !// ****************************************************************

        !// Iterative scaling (three rounds):
        if (mod(NST,scalemod) .eq. 0) then
           !// For NCHEM_ITER=12, scales at NST=3,6,9,12
           do ISCALE = 1, 3

              !// Brx:
              XBrx = M_Bry 
              YBrx = M_Br + M_BrO + M_BrONO2 + M_HOBr &
                     + M_HBr + 2._r8*M_Br2 + M_BrCL

              FACN = XBrx / YBrx

              M_Br   = M_Br * FACN
              M_BrO  = M_BrO * FACN
              M_HBr  = M_HBr * FACN
              M_BrONO2 = M_BrONO2 * FACN
              M_HOBr = M_HOBr * FACN
              M_Br2  = M_Br2 * FACN
              M_BrCl = M_BrCl * FACN

              !// NOx:
              XVNOX = M_NO + M_NO2 + M_NO3 + 2._r8*M_N2O5 &
                      + M_HO2NO2 + M_ClONO2 + M_BrONO2
              VNOY = M_NOx_str

              FACN = VNOy / XVNOx

              M_NO   = M_NO * FACN
              M_NO2  = M_NO2 * FACN
              M_NO3  = M_NO3 * FACN
              M_N2O5 = M_N2O5 * FACN
              !//HNO2: M_HNO2 = M_HNO2 *FACN
              M_HO2NO2 = M_HO2NO2 * FACN
              M_ClONO2 = M_ClONO2 * FACN
              M_BrONO2 = M_BrONO2 * FACN

              !// Clx:
              XClx = M_Clx
              YClx = M_Cl + M_ClO + M_HOCl + M_ClONO2 &
                     + 2._r8*M_Cl2 + M_OCLO + M_BrCL + M_ClOO + 2._r8*M_Cl2O2

              FACN = XClx / YClx

              M_Cl   = M_Cl * FACN
              M_ClO  = M_ClO * FACN
              M_HOCl = M_HOCl * FACN
              M_ClONO2 = M_ClONO2 * FACN
              M_Cl2  = M_Cl2 * FACN
              M_OClO = M_OClO * FACN
              M_ClOO = M_ClOO * FACN
              M_Cl2O2= M_Cl2O2 * FACN
              M_BrCl = M_BrCl * FACN

           end do !// do ISCALE = 1, 3

           xCly = M_Cly
           yCly = M_ClX + M_HCl

           FACN = xCly / yCly

           M_Clx = M_Clx * FACN
           M_HCl = M_HCl * FACN


           xNOy = M_NOy_str
           yNOy = M_NOx_str + M_HNO3

           FACN = xNOy / yNOy

           M_NOx_str = M_NOx_str * FACN
           M_HNO3 = M_HNO3 * FACN


           !// --------------------------------------------------------
           !// As the total families of NOx, Clx and Bry were calculated
           !// per se during the integration and those results are thought
           !// to be the 'final ones', they stay the same and are not
           !// affected by the scaling iteration. The members within
           !// the families are scaled, instead.
           !// --------------------------------------------------------
           !// Assumedly, one should NOT mess with the calculated SO!
           !// But, does the so family need to be included in the
           !// iterative scaling? Probably not in normal stratospheric
           !// chemistry, but the situation could be different under
           !// mesospheric conditions.
           !// Scaling of SO:
           !// Assumedly, the most correct value for SO is the
           !// integrated family concentration.
           !// At this point: this value is stored to SO.
           !//   O3, OP, OD are as in the end of the time integration.
           !//   Cl, Br and NO have been modified by scaling above.
           !// --------------------------------------------------------
           !// Extract O3 from SO and the other members of SO:
           M_O3 = M_SO - M_O3P - M_O1D + M_NO + M_Cl + M_Br
           if (M_O3 .lt. 1.e-6_r8) then
              write(6,'(a,3i4,4es9.2)') 'pchemc_str_ij.f90: Neg/Tiny O3 D',&
                   L,ICOL,JCOL, M_O3P,M_O3,M_SO,M_O1D
              stop
           end if

           !// There is no need to scale the SH family as the integrated
           !// SH is taken as the correct result and the members of the
           !// family have been scaled to the total family value in the
           !// course of the chemical calculations.
           !// --------------------------------------------------------
           !// End ifscaling over families...
        end if !// if (mod(NST,scalemod) .eq. 0) then

        !// Assign chemical components (prefix 'M' for molecular density)
        ZC_LOCAL( 1,L) = M_O3
        ZC_LOCAL( 4,L) = M_HNO3
        ZC_LOCAL( 6,L) = M_CO
        ZC_LOCAL(13,L) = M_CH2O
        ZC_LOCAL(15,L) = M_H2O2
        ZC_LOCAL(16,L) = M_CH3O2H
        ZC_LOCAL(17,L) = M_HO2NO2
        ZC_LOCAL(21,L) = M_HO2
        ZC_LOCAL(22,L) = M_CH3O2
        ZC_LOCAL(38,L) = M_O3P
        ZC_LOCAL(39,L) = M_O1D
        ZC_LOCAL(40,L) = M_OH
        ZC_LOCAL(41,L) = M_NO3
        ZC_LOCAL(42,L) = M_N2O5
        ZC_LOCAL(43,L) = M_NO
        ZC_LOCAL(44,L) = M_NO2
        ZC_LOCAL(46,L) = M_CH4
        ZC_LOCAL(101,L) = M_MCF
        ZC_LOCAL(102,L) = M_HCFC22
        ZC_LOCAL(103,L) = M_CFC11
        ZC_LOCAL(104,L) = M_CFC12
        ZC_LOCAL(105,L) = M_CCl4
        ZC_LOCAL(106,L) = M_CH3Cl
        ZC_LOCAL(107,L) = M_N2O
        ZC_LOCAL(108,L) = M_Clx
        ZC_LOCAL(109,L) = M_NOx_str
        ZC_LOCAL(110,L) = M_SO
        ZC_LOCAL(111,L) = M_HCl
        ZC_LOCAL(112,L) = M_Cly
        ZC_LOCAL(113,L) = M_H2
        if (.not. LOLD_H2OTREATMENT) ZC_LOCAL(114,L) = M_H2O
        ZC_LOCAL(115,L) = M_SH
        ZC_LOCAL(116,L) = M_CH3Br
        ZC_LOCAL(117,L) = M_H1211
        ZC_LOCAL(118,L) = M_H1301
        ZC_LOCAL(119,L) = M_Bry
        ZC_LOCAL(120,L) = M_H2402
        ZC_LOCAL(121,L) = M_CFC113
        ZC_LOCAL(122,L) = M_CFC114
        ZC_LOCAL(123,L) = M_CFC115
        ZC_LOCAL(124,L) = M_HNO3s
        ZC_LOCAL(125,L) = M_H2Os
        ZC_LOCAL(127,L) = M_HCFC123
        ZC_LOCAL(128,L) = M_HCFC141
        ZC_LOCAL(129,L) = M_HCFC142
        ZC_LOCAL(130,L) = M_H
        ZC_LOCAL(131,L) = M_HNO2
        ZC_LOCAL(132,L) = M_Cl
        ZC_LOCAL(133,L) = M_ClO
        ZC_LOCAL(134,L) = M_HOCl
        ZC_LOCAL(135,L) = M_ClONO2
        ZC_LOCAL(136,L) = M_Cl2
        ZC_LOCAL(137,L) = M_OClO
        ZC_LOCAL(138,L) = M_Br
        ZC_LOCAL(139,L) = M_BrO
        ZC_LOCAL(140,L) = M_HBr
        ZC_LOCAL(141,L) = M_BrONO2
        ZC_LOCAL(142,L) = M_HOBr
        ZC_LOCAL(143,L) = M_Br2
        ZC_LOCAL(144,L) = M_ClOO
        ZC_LOCAL(145,L) = M_Cl2O2
        ZC_LOCAL(146,L) = M_BrCl
        ZC_LOCAL(147,L) = M_NOy_str
        ZC_LOCAL(148,L) = M_H2O_ac

      end do ! NST=1,NCHEM_ITER

      if (LDEBUG_STRCHEM) then
         do N = 1,TRACER_ID_MAX
            if (trsp_idx(N) .gt. 0 .or. Xtrsp_idx(N) .gt. 0) Then
               !// Check for negative concentrations:
               if (ZC_LOCAL(N,L) .lt. 0._r8) Then
                  write(*,'(a,5i4,f20.6)') &
                       'OSLO_CHEM_STR_IJ: 2: NEG.CONC. FOR N,I,J,L: ', &
                       N,ICOL,JCOL,L,lmtrop,ZC_LOCAL(N,L)
                  stop
               end if
               !// Check for NANQs:
               if (1._r8*ZC_LOCAL(N,L) .ne. ZC_LOCAL(N,L)) then
                  print*,'OSLO_CHEM_STR_IJ: NaNQ FOR N,I,J,L : ',N,ICOL,JCOL,L
                  stop
               end if
            end if
         end do
      end if

      !// ------------------------------------------------------------------
    end do !// do L = LMTROP+1,LM-1
    !// --------------------------------------------------------------------

    !// --------------------------------------------------------------------
  end subroutine OSLO_CHEM_STR_IJ
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module pchemc_str_ij
!//=========================================================================

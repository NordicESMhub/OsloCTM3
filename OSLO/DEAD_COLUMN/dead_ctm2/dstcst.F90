! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstcst.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Fundamental physical constants needed by dust routines

! Usage: 
! use dstcst ! [mdl] Physical constants for dust routines

module dstcst ! [mdl] Physical constants for dust routines
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  save ! [stt] Changes to common variables are sticky
  
  ! These variables are initialized in dstmss:dst_cst_cmn_ini()
  ! dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun() 
  ! CCM:physics/inti() passes values from CCM:physics/comcon.h which were set by CCM:eul/initcom()
  real(r8) cst_von_krm          ! [frc] Von Karman's constant
  real(r8) cst_von_krm_rcp      ! [frc] Reciprocal of Von Karman's constant
  real(r8) eps_H2O              ! (0.622) [frc] mmw(H2O)/mmw(dry air)
  real(r8) eps_H2O_rcp_m1       ! [frc] (0.60777) Constant for virtual temperature
  real(r8) gas_cst_dry_air      ! [J kg-1 K-1] (287.05) Gas constant of dry air
  real(r8) gas_cst_unv          ! [J mol-1 K-1] Universal gas constant
  real(r8) grv_sfc              ! [m s-2] (9.80616) Gravity (mean surface)
  real(r8) grv_sfc_rcp          ! [s2 m-1] (0.101977) Reciprocal of gravity
  real(r8) kappa_dry_air        ! (0.286 = 2/7) [frc] Constant in potential temperature
  real(r8) mmw_H2O              ! [kg mol-1] Mean molecular weight of water
  real(r8) mmw_dry_air          ! [kg mol-1] Mean molecular weight of dry air
  real(r8) rds_Earth            ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
  
contains
  
  subroutine dst_cst_cmn_ini(gas_cst_dry_airx,grv_sfcx,rds_Earthx)
    ! Initialize time-invariant physical constants common block dstcstcmn
    ! dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
    ! CCM:physics/inti() passes values from CCM:physics/comcon.h which were set by CCM:eul/initcom()
    use precision ! [mdl] Precision r8, i8, ...
    implicit none
    ! Parameters for source and sink processes
    real(r8),parameter::gas_cst_unvx=8.31441 ! [J mol-1 K-1] Universal gas constant
    real(r8),parameter::mmw_H2Ox=1.8015259e-02 ! [kg mol-1] Mean molecular weight of water HITRAN
    real(r8),parameter::mmw_dry_airx=28.9644e-3 ! [kg mol-1] Mean molecular weight of dry air CCM2/3:physics/radcsw()
    real(r8),parameter::cst_von_krmx=0.4 ! [fraction] Von Karman constant
    ! Input
    real(r8),intent(in)::grv_sfcx             ! [m s-2] (9.80616) Mean gravitational acceleration at Earth's surface
    real(r8),intent(in)::gas_cst_dry_airx     ! [J kg-1 K-1] (287.05) Gas constant of dry air
    real(r8),intent(in)::rds_Earthx           ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
    ! Local
    real(r8) gas_cst_H2O          ! (461.65) [J kg-1 K-1] Gas constant of H2O
    ! Local
    real(r8) spc_heat_dry_air ! (1005.0) [J kg-1 K-1]
    ! Main code
    
    ! Parameters for source and sink processes
    grv_sfc=grv_sfcx          ! [m s-2] (9.80616) Mean gravitational acceleration at Earth's surface
    grv_sfc_rcp=1.0_r8/grv_sfcx  ! [s2 m-1] (0.101977) Reciprocal of gravity
    gas_cst_dry_air=gas_cst_dry_airx ! [J kg-1 K-1] (287.05) Gas constant of dry air
    rds_Earth=rds_Earthx      ! [m] (6.37122e+6) Radius of sphere of same volume as Earth
    mmw_H2O=mmw_H2Ox          ! [kg mol-1] Mean molecular weight of dry air
    mmw_dry_air=mmw_dry_airx  ! [kg mol-1] Mean molecular weight of dry air
    gas_cst_unv=gas_cst_unvx  ! [J mol-1 K-1] Universal gas constant
    cst_von_krm=cst_von_krmx  ! [frc] Von Karman's constant
    cst_von_krm_rcp=1.0_r8/cst_von_krmx ! [frc] Reciprocal of Von Karman's constant
    eps_H2O=mmw_H2Ox/mmw_dry_airx ! [frc] (0.622) mmw(H2O)/mmw(dry air)
    eps_H2O_rcp_m1=-1.0_r8+mmw_dry_airx/mmw_H2Ox ! [frc] (0.60777) Constant for virtual temperature

    ! Derived 
    spc_heat_dry_air=7.0_r8*gas_cst_dry_airx/2.0_r8 ! (1005.0) [J kg-1 K-1]
    kappa_dry_air=gas_cst_dry_airx/spc_heat_dry_air ! (0.286 = 2/7) [frc] Constant in potential temperature (IrG81 p. 25, Tre922 p. 72) 
    
    ! Sanity checks
    if (abs(0.28571_r8-kappa_dry_air)/0.28571_r8 > 1.0e-4_r8) stop 'dst_cst_cmn_ini(): kappa_dry_air error'
    if (abs(0.60777_r8-eps_H2O_rcp_m1)/0.60777_r8 > 1.0e-4_r8) stop 'dst_cst_cmn_ini(): eps_H2O_rcp_m1 error'
    
    return
  end subroutine dst_cst_cmn_ini                       ! end dst_cst_cmn_ini()
  
end module dstcst

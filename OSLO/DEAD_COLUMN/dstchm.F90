! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstchm.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Module dstchm contains chemical properties and subroutines for dust model

! dstchm MUST have access to dstgrd.F90 and to pmgrid.F90 to work

! Usage: 
! use dstchm ! [mdl] Chemical properties of dust

module dstchm ! [mdl] Chemical properties of dust
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none

  public::dst_chm_cmn_ini ! [sbr] Initialize chemistry common blocks
!  public::dst_chm_rxr ! [sbr] Heterogeneous chemical processes on dust
!  public::dst_chm_slv ! [sbr] Chemistry solver for heterogeneous processes
!  private::chm2nc ! [sbr] Write aerosol chemical properties to netCDF file
  
  integer,public,parameter::idx_H2O2_gas=1 ! [idx] Index of gaseous H2O2
  integer,public,parameter::idx_HNO3_gas=2 ! [idx] Index of gaseous HNO3
  integer,public,parameter::idx_HO2_gas=3 ! [idx] Index of gaseous HO2
  integer,public,parameter::idx_N2O5_gas=4 ! [idx] Index of gaseous N2O5
  integer,public,parameter::idx_NO3_gas=5 ! [idx] Index of gaseous NO3
  integer,public,parameter::idx_O3_gas=6 ! [idx] Index of gaseous O3
  integer,public,parameter::idx_OH_gas=7 ! [idx] Index of gaseous OH
  integer,public,parameter::idx_SO2_gas=8 ! [idx] Index of gaseous SO2
  integer,public,parameter::idx_SO4_aer=9 ! [idx] Index of particulate SO4
  integer,public,parameter::idx_NO3_aer=10 ! [idx] Index of particulate NO3
  integer,public,parameter::idx_chm_end=10 ! [idx] Last chemical index
  
  ! Chemical properties of dust
  ! These variables are initialized in dstchm.F:dst_chm_cmn_ini()
  ! dst_chm_cmn_ini() is called by dst_msc_cmn_ini()
  ! fxm: Eventually, uptake coefficients should be moved to dst_chm_rxr()
  real(r8),private::upt_cff_H2O2_dst     ! [frc] Uptake coefficient for H2O2 to dust
  real(r8),private::upt_cff_HNO3_dst     ! [frc] Uptake coefficient for HNO3 to dust
  real(r8),private::upt_cff_HO2_dst      ! [frc] Uptake coefficient for HO2 to dust
  real(r8),private::upt_cff_N2O5_dst     ! [frc] Uptake coefficient for N2O5 to dust
  real(r8),private::upt_cff_NO2_dst      ! [frc] Uptake coefficient for NO2 to dust
  real(r8),private::upt_cff_NO3_dst      ! [frc] Uptake coefficient for NO3 to dust
  real(r8),private::upt_cff_O3_dst       ! [frc] Uptake coefficient for O3 to dust
  real(r8),private::upt_cff_OH_dst       ! [frc] Uptake coefficient for OH to dust
  real(r8),private::upt_cff_SO2_dst      ! [frc] Uptake coefficient for SO2 to dust
  real(r8),private::vmr_HNO3_gas_mxm     ! [mlc mlc-1] Maximum HNO3 volume mixing ratio
  
contains
  
  subroutine dst_chm_cmn_ini()
    ! Purpose: Initialize chemistry common blocks
    ! dst_chm_cmn_ini() is called by dst_msc_cmn_ini()
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    ! Output
    ! Input/Output
    ! Local
    ! Main code
    
    ! Sanity check
    if (idx_chm_end /= chm_nbr) stop 'dst_chm_cmn_ini() reports idx_chm_end /= chm_nbr'
    
    ! H2O2 --> dust 
    ! ZhC99 p. 357 Tbl. 2 use gamma(H2O2->dust)=1.0e-4
    upt_cff_H2O2_dst=1.0e-4_r8   ! [frc] Uptake coefficient for H2O2 to dust
    
    ! HNO3 --> dust 
    ! Vicki Grassian (personal communication, 1999) suggests gamma(HNO3->dust) ~ 0.005
    ! DCZ96 p. 22872 (9) use gamma(HNO3->dust)=0.1 when soil balance is alkaline
    ! ZhC99 p. 357 Tbl. 2 use gamma(HNO3->dust)=0.01
    ! USP02 p.     Tbl. 2 use gamma(HNO3->China loess)=5.2e-5
    ! USP02 p.     Tbl. 2 use gamma(HNO3->Saharan sand)=2.0e-5
    ! USP02 p.     Tbl. 4 use gamma(HNO3->China loess)=1.1e-3
    upt_cff_HNO3_dst=1.1e-3 ! [frc] Uptake coefficient for HNO3 to dust
    
    ! HO2 --> dust 
    ! DCZ96 p. 22871 use gamma(HO2->dust)=0.1
    ! ZhC99 p. 357 Tbl. 2 use gamma(HO2->dust)=0.1
    upt_cff_HO2_dst=0.1_r8       ! [frc] Uptake coefficient for HO2 to dust
    
    ! N2O5 --> dust 
    ! DCZ96 p. 22871 use gamma(N2O5->dust)=0.1 but suggest ~0.001 for low RH
    ! ZhC99 p. 357 Tbl. 2 use gamma(N2O5->dust)=0.1
    upt_cff_N2O5_dst=0.001_r8    ! [frc] Uptake coefficient for N2O5 to dust
    
    ! NO2 --> dust 
    ! USP02 p.     Tbl. 2 use gamma(NO2->China loess)=2.1e-6
    ! USP02 p.     Tbl. 2 use gamma(NO2->Saharan sand)=1.2e-6
    ! USP02 p.     Tbl. 4 use gamma(NO2->China loess)=4.4e-5
    upt_cff_NO2_dst=4.4e-5_r8     ! [frc] Uptake coefficient for NO2 to dust

    ! NO3 --> dust 
    ! ZhC99 p. 357 Tbl. 2 use gamma(NO3->dust)=0.1
    upt_cff_NO3_dst=0.1_r8       ! [frc] Uptake coefficient for NO3 to dust
    
    ! O3 --> dust 
    ! DCZ96 p. 22874 use gamma(O3->dust)=5.0e-5 based on inversion of deposition measurements
    ! ZhC99 p. 357 Tbl. 2 use gamma(O3->dust)=1.0e-4
    ! MUG02 p.     use gamma(O3->China loess)=2.9e-5
    ! MUG02 p.     use gamma(O3->Saharan sand)=6.0e-5
    upt_cff_O3_dst=5.0e-5_r8     ! [frc] Uptake coefficient for O3 to dust
    
    ! OH --> dust 
    ! ZhC99 p. 357 Tbl. 2 use gamma(OH->dust)=0.1
    upt_cff_OH_dst=0.1_r8        ! [frc] Uptake coefficient for OH to dust
    
    ! SO2 --> dust 
    ! DCZ96 p. 22872 use gamma(SO2->dust)=3.0e-4 when soil balance is alkaline
    ! ZhC99 p. 357 Tbl. 2 use gamma(SO2->dust)=1.0e-4
    upt_cff_SO2_dst=3.0e-4_r8    ! [frc] Uptake coefficient for SO2 to dust
    
    vmr_HNO3_gas_mxm=1.0e-3_r8   ! [mlc mlc-1] Maximum gaseous HNO3 volume mixing ratio
    return
  end subroutine dst_chm_cmn_ini
  
  
end module dstchm

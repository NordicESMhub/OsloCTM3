! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstchm.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Module dstchm contains chemical properties and subroutines for dust model

! dstchm MUST have access to dstgrd.F90 and to pmgrid.F90 to work

! Usage: 
! use dstchm ! [mdl] Chemical properties of dust

module dstchm ! [mdl] Chemical properties of dust
  use precision ! [mdl] Precision r8, i8, ...
  implicit none

  public::dst_chm_cmn_ini ! [sbr] Initialize chemistry common blocks
  public::dst_chm_rxr ! [sbr] Heterogeneous chemical processes on dust
  public::dst_chm_slv ! [sbr] Chemistry solver for heterogeneous processes
  private::chm2nc ! [sbr] Write aerosol chemical properties to netCDF file
  
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
    use precision ! [mdl] Precision r8, i8, ...
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
  
  subroutine dst_chm_rxr(lchnk,ncol,obuf, &
       dns_mdp,             & ! I [kg m-3] Midlayer density
       lat_idx,             & ! I [idx] Model latitude index
       prs_dlt,             & ! I [Pa] Pressure thickness
       prs_mdp,             & ! I [Pa] Pressure
       q_dst,               & ! I [kg kg-1] Dust mixing ratio
       rxrc_HNO3_dst,       & ! O [s-1] Pseudo first order rate coefficient for HNO3
       rxrc_HO2_dst,        & ! O [s-1] Pseudo first order rate coefficient for HO2
       rxrc_N2O5_dst,       & ! O [s-1] Pseudo first order rate coefficient for N2O5
       rxrc_O3_dst,         & ! O [s-1] Pseudo first order rate coefficient for O3
       rxrc_SO2_dst,        & ! O [s-1] Pseudo first order rate coefficient for SO2
       tpt_mdp,             & ! I [K] Temperature
       tpt_vrt,             & ! I [K] Virtual temperature
       vmr_SO4_aer)         ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
    ! Heterogeneous chemical processes on dust
    ! dst_chm_rxr() will be called by MATCH:src/physlic()
    use dstaer ! [mdl] Aerosol microphysical properties
    use dstcst ! [mdl] Physical constants for dust routines
    use dstgrd ! [mdl] Dust grid sizes
    use dstmssutl,only:dst_add_lon_lev,dst_add_lon ! [mdl] Mass budget utilities
    use pmgrid ! [mdl] Spatial resolution parameters
    use precision ! [mdl] Precision r8, i8, ...
    implicit none
    ! Parameters
    real(r8),parameter::cst_Avagadro=6.022045e+23 ! [mlc mol-1] Avagadro's number
    real(r8),parameter::mmw_Ca=40.078e-03 ! [kg mol-1] Mean molecular weight of Ca
    real(r8),parameter::mmw_CaCO3=100.087e-03 ! [kg mol-1] Mean molecular weight of CaCO3
    real(r8),parameter::mmw_HNO3=62.995644e-03 ! [kg mol-1] Mean molecular weight of HNO3
    real(r8),parameter::mmw_HO2=33.0067e-03 ! [kg mol-1] Mean molecular weight of HO2
    real(r8),parameter::mmw_N2O5=108.01e-03 ! [kg mol-1] Mean molecular weight of N2O5
    real(r8),parameter::mmw_NO3=62.0049e-03 ! [kg mol-1] Mean molecular weight of NO3
    real(r8),parameter::mmw_O3=47.997832e-03 ! [kg mol-1] Mean molecular weight of O3
    real(r8),parameter::mmw_SO2=64.046674e-03 ! [kg mol-1] Mean molecular weight of SO2
    real(r8),parameter::mmw_SO4=96.0636e-03 ! [kg mol-1] Mean molecular weight of SO4
    ! Global continental mean CaCO3 content from IGBP-DIS is 2.5%
    ! ncwa -O -m lnd_msk -M 1.0 -o eq -a lat,lon -w gw -C -v mss_frc_CaCO3 dst_T42.nc foo.nc ; ncks -H foo.nc
    ! Approximate mean over dust source regions is 5%
    ! NB: DCZ96 quote sources saying Ca (not CaCO3) content is 3.6%, which clearly disagrees with IGBP-DIS
    real(r8),parameter::mss_frc_CaCO3=0.05 ! [frc] Soil calcium carbonate fraction 
    ! Input
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::dns_mdp(plond,plev) ! I [kg m-3] Midlayer density
    real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
    real(r8),intent(in)::prs_dlt(plond,plev) ! I [Pa] Pressure thickness
    real(r8),intent(in)::prs_mdp(plond,plev) ! I [Pa] Midlayer pressure
    real(r8),intent(in)::tpt_mdp(plond,plev) ! I [K] Temperature
    real(r8),intent(in)::tpt_vrt(plond,plev) ! I [K] Virtual temperature
    real(r8),intent(in)::q_dst(plond,plev,dst_nbr) ! I [kg kg-1] Dust mixing ratio
    real(r8),intent(in)::vmr_SO4_aer(plond,plev) ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
    ! Input/Output
    ! Output
    real(r8),intent(out)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
    real(r8),intent(out)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
    real(r8),intent(out)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
    real(r8),intent(out)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
    real(r8),intent(out)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
    ! Local Output
    ! Local
    character(80) fl_out       ! [sng] Name of netCDF output file
    real(r8) pi       ! [frc] 3
    integer i                 ! [idx] Counting index over lon
    integer k                 ! [idx] Counting index over lev
    integer m                 ! [idx] Counting index over cst
    real(r8) mpl_air(plond,plev)  ! [kg m-2] Air mass path in layer
    real(r8) sfc_ttl(plond,plev)  ! [m2 kg-1] Total surface area
    real(r8) tm_dlt               ! [s] Chemistry timestep
    real(r8) q_dst_ttl(plond,plev) ! [kg kg-1] Total dust mixing ratio
    real(r8) rxrc_fct(plond,plev) ! [frc] Rate coefficient factor = (1/4)*area*sqrt(T)
    real(r8) vlc_MWB_HNO3_fct     ! [m s-1] Thermal speed of HNO3 factor
    real(r8) vlc_MWB_HO2_fct      ! [m s-1] Thermal speed of HO2 factor
    real(r8) vlc_MWB_N2O5_fct     ! [m s-1] Thermal speed of N2O5 factor
    real(r8) vlc_MWB_O3_fct       ! [m s-1] Thermal speed of O3 factor
    real(r8) vlc_MWB_SO2_fct      ! [m s-1] Thermal speed of SO2 factor
    real(r8) vmr_CaCO3_aer(plond,plev) ! [mlc mlc-1] Particulate CaCO3 volume mixing ratio
    real(r8) vmr_CaCO3_fct        ! [kg kg-1] Factor needed for CaCO3 
    real(r8) vmr_NO3_aer(plond,plev) ! [mlc mlc-1] Particulate NO3 volume mixing ratio
    ! Local output
    ! GCM diagnostics
    
    ! Main Code
    
    ! Initialize fluxes and tendencies
    pi=4.0_r8*atan(1.0_r8)        ! [frc] 3
    sfc_ttl(:,:)=0.0_r8 ! [m2 kg-1] Total surface area
    rxrc_HNO3_dst(:,:)=0.0_r8 ! [s-1] Pseudo first order rate coefficient for HNO3
    rxrc_HO2_dst(:,:)=0.0_r8 ! [s-1] Pseudo first order rate coefficient for HO2
    rxrc_N2O5_dst(:,:)=0.0_r8 ! [s-1] Pseudo first order rate coefficient for N2O5
    rxrc_O3_dst(:,:)=0.0_r8 ! [s-1] Pseudo first order rate coefficient for O3
    rxrc_SO2_dst(:,:)=0.0_r8 ! [s-1] Pseudo first order rate coefficient for SO2
    vmr_NO3_aer(:,:)=0.0_r8 ! [mlc mlc-1] Particulate NO3 volume mixing ratio
    
    ! Time invariant scalar factor for thermal molecular speeds faciliates rate coefficient computation
    vlc_MWB_HNO3_fct=sqrt(8.0_r8*gas_cst_unv/(pi*mmw_HNO3)) ! [m s-1] Thermal speed of HNO3 factor
    vlc_MWB_HO2_fct=sqrt(8.0_r8*gas_cst_unv/(pi*mmw_HO2)) ! [m s-1] Thermal speed of HO2 factor
    vlc_MWB_N2O5_fct=sqrt(8.0_r8*gas_cst_unv/(pi*mmw_N2O5)) ! [m s-1] Thermal speed of N2O5 factor
    vlc_MWB_O3_fct=sqrt(8.0_r8*gas_cst_unv/(pi*mmw_O3)) ! [m s-1] Thermal speed of O3 factor
    vlc_MWB_SO2_fct=sqrt(8.0_r8*gas_cst_unv/(pi*mmw_SO2)) ! [m s-1] Thermal speed of SO2 factor
    
    ! Compute necessary derived fields
    ! Mass of air currently in gridbox
    mpl_air=prs_dlt*grv_sfc_rcp ! [kg m-2]
    
    ! Surface area and CaCO3 mixing ratio
    call dst_add_lon_lev(q_dst,q_dst_ttl) ! [kg kg-1] Total dust mixing ratio
    do m=1,dst_nbr
       do k=1,plev
          do i=1,plon
             sfc_ttl(i,k)=sfc_ttl(i,k) & ! [m2 kg-1] Total surface area
                  +q_dst(i,k,m)*mpl_air(i,k)*sfc_spc_rsl(m)
          end do              ! end loop over lon
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    
    ! Volume mixing ratios for electroneutrality condition
    vmr_CaCO3_fct=mmw_dry_air/mmw_CaCO3 ! [kg kg-1] Factor needed for CaCO3 volume mixing ratio
    vmr_CaCO3_aer=            & ! [mlc mlc-1] Particulate CaCO3 volume mixing ratio
         q_dst_ttl*mss_frc_CaCO3 & ! Converts dust to CaCO3 mass mixing ratio
         *vmr_CaCO3_fct       ! Converts CaCO3 mmr to CaCO3 vmr
    
    ! Pseudo first order reaction rate coefficients asumme Maxwell-Boltzmann statistics and kinetic regime (surface area limited) uptake
    ! k = ( Aerosol surface area * gamma * molecular velocity ) / 4.0
    rxrc_fct=0.25_r8*sfc_ttl*sqrt(tpt_mdp) ! [frc] Rate coefficient factor = (1/4)*area*sqrt(T)
    rxrc_HO2_dst=rxrc_fct*upt_cff_HO2_dst*vlc_MWB_HO2_fct ! [s-1] Pseudo first order rate coefficient for HO2
    rxrc_N2O5_dst=rxrc_fct*upt_cff_N2O5_dst*vlc_MWB_N2O5_fct ! [s-1] Pseudo first order rate coefficient for N2O5
    rxrc_O3_dst=rxrc_fct*upt_cff_O3_dst*vlc_MWB_O3_fct ! [s-1] Pseudo first order rate coefficient for O3
    
    ! Electroneutrality condition based on DCZ96 p. 22872:
    ! Uptake of SO2 and of HNO3 proceeds only when aerosol contains enough basic (alkaline) material to neutralize the acidity of the added species
    ! Stoichiometric relationship is based on mineral dust samples from Asia
    ! Reaction only proceeds when threshold is exceeded, so process causes some spikiness in output
    where (vmr_NO3_aer+2.0_r8*vmr_SO4_aer < 2.0_r8*vmr_CaCO3_aer) 
       rxrc_HNO3_dst=rxrc_fct*upt_cff_HNO3_dst*vlc_MWB_HNO3_fct ! [s-1] Pseudo first order rate coefficient for HNO3
       rxrc_SO2_dst=rxrc_fct*upt_cff_SO2_dst*vlc_MWB_SO2_fct ! [s-1] Pseudo first order rate coefficient for SO2
    end where                 ! end where electroneutrality condition satisfied
    
    return
  end subroutine dst_chm_rxr
  
  subroutine dst_chm_slv(lchnk,ncol,obuf, &
       lat_idx,             & ! I [idx] Model latitude index
       rxrc_HNO3_dst,       & ! I [s-1] Pseudo first order rate coefficient for HNO3
       rxrc_HO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for HO2
       rxrc_N2O5_dst,       & ! I [s-1] Pseudo first order rate coefficient for N2O5
       rxrc_O3_dst,         & ! I [s-1] Pseudo first order rate coefficient for O3
       rxrc_SO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for SO2
       tm_adj,              & ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
       vmr_HNO3_gas,        & ! I/O [mlc mlc-1] Gaseous HNO3 volume mixing ratio
       vmr_HO2_gas,         & ! I/O [mlc mlc-1] Gaseous HO2 volume mixing ratio
       vmr_N2O5_gas,        & ! I/O [mlc mlc-1] Gaseous N2O5 volume mixing ratio
       vmr_NO3_aer,         & ! I/O [mlc mlc-1] Particulate NO3 volume mixing ratio
       vmr_O3_gas,          & ! I/O [mlc mlc-1] Gaseous O3 volume mixing ratio
       vmr_SO2_gas,         & ! I/O [mlc mlc-1] Gaseous SO2 volume mixing ratio
       vmr_SO4_aer)         ! I/O [mlc mlc-1] Particulate SO4 volume mixing ratio
    ! Purpose: Chemistry solver for heterogeneous chemical processes on dust
    ! dst_chm_slv() will be called by MATCH:src/physlic()
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use precision ! [mdl] Precision r8, i8, ...
    use sng_mdl,only:ftn_strlen,ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Parameters
    ! Input
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    integer,intent(in)::lchnk ! I [id] Chunk identifier
    integer,intent(in)::ncol ! I [nbr] Number of atmospheric columns
    real(r8),intent(in)::obuf(*)  ! I [ptr] Output buffer
    real(r8),intent(in)::tm_adj   ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
    real(r8),intent(in)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
    real(r8),intent(in)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
    real(r8),intent(in)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
    real(r8),intent(in)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
    real(r8),intent(in)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
    ! Input/Output
    real(r8),intent(inout)::vmr_HNO3_gas(plond,plev) ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
    real(r8),intent(inout)::vmr_HO2_gas(plond,plev) ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
    real(r8),intent(inout)::vmr_N2O5_gas(plond,plev) ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
    real(r8),intent(inout)::vmr_O3_gas(plond,plev) ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
    real(r8),intent(inout)::vmr_SO2_gas(plond,plev) ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
    real(r8),intent(inout)::vmr_NO3_aer(plond,plev) ! I/O [mlc mlc-1] Particulate NO3 volume mixing ratio
    real(r8),intent(inout)::vmr_SO4_aer(plond,plev) ! I/O [mlc mlc-1] Particulate SO4 volume mixing ratio
    ! Output
    ! Local Output
    real(r8) rxr_HNO3_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
    real(r8) rxr_HO2_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of HO2 removal by dust
    real(r8) rxr_N2O5_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of N2O5 removal by dust
    real(r8) rxr_SO2_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of SO2 removal by dust
    real(r8) rxr_O3_gas_dst_vmr(plond,plev) ! O [mlc mlc-1 s-1] Mean rate of O3 removal by dust
    ! Local
    character(80) fl_out      ! [sng] Name of netCDF output file
    integer i                 ! [idx] Counting index over lon
    integer k                 ! [idx] Counting index over lev
    integer m                 ! [idx] Counting index over cst
    real(r8) tm_dlt               ! [s] Chemistry timestep
    real(r8) vmr_HNO3_gas_dlt(plond,plev) ! [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio
    real(r8) vmr_SO2_gas_dlt(plond,plev) ! [mlc mlc-1] Change in gaseous SO2 volume mixing ratio
    ! Local output
    ! GCM diagnostics
    
    ! Main Code
    ! Timesplit if desired
    tm_dlt=tm_adj             ! [s] Chemistry timestep
    
    ! Initialize fluxes and tendencies
    rxr_HNO3_gas_dst_vmr=0.0_r8  ! [mlc mlc-1 s-1] Rate of HNO3 removal by dust
    rxr_HO2_gas_dst_vmr=0.0_r8   ! [mlc mlc-1 s-1] Rate of HO2 removal by dust
    rxr_N2O5_gas_dst_vmr=0.0_r8  ! [mlc mlc-1 s-1] Rate of N2O5 removal by dust
    rxr_O3_gas_dst_vmr=0.0_r8    ! [mlc mlc-1 s-1] Rate of O3 removal by dust
    rxr_SO2_gas_dst_vmr=0.0_r8   ! [mlc mlc-1 s-1] Rate of SO2 removal by dust
    
    ! Compute heterogeneous removal rates and adjust mixing ratios
    
    ! HNO3
    ! 20030114: Possible floating point underflow here in single precision on ALPHA
    vmr_HNO3_gas_dlt=vmr_HNO3_gas*(1.0_r8-exp(-rxrc_HNO3_dst*tm_dlt)) ! [mlc mlc-1] Change in gaseous HNO3 volume mixing ratio
    rxr_HNO3_gas_dst_vmr=vmr_HNO3_gas_dlt/tm_dlt ! [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
    vmr_HNO3_gas=vmr_HNO3_gas-vmr_HNO3_gas_dlt ! [mlc mlc-1] Gaseous HNO3 volume mixing ratio
    vmr_NO3_aer=vmr_NO3_aer+vmr_HNO3_gas_dlt ! [mlc mlc-1] Particulate NO3 volume mixing ratio
    
    ! SO2
    vmr_SO2_gas_dlt=vmr_SO2_gas*(1.0_r8-exp(-rxrc_SO2_dst*tm_dlt)) ! [mlc mlc-1] Change in gaseous SO2 volume mixing ratio
    vmr_SO2_gas=vmr_SO2_gas-vmr_SO2_gas_dlt ! [mlc mlc-1] Gaseous SO2 volume mixing ratio
    vmr_SO4_aer=vmr_SO4_aer+vmr_SO2_gas_dlt ! [mlc mlc-1] Particulate SO4 volume mixing ratio
    
#ifdef DST_DBG
    do k=1,plev
       do i=1,plon
          if (vmr_HNO3_gas(i,k) > vmr_HNO3_gas_mxm) write(6,'(a,i2,a,i3,a,i2,a,es8.1,a)')  &
               'dst_chm_slv: lat = ',lat_idx,' vmr_HNO3_gas(',i,',',m,') = ',vmr_HNO3_gas(i,k),' mlc mlc-1'
       end do                 ! end loop over lon
    end do                    ! end loop over lev
#endif /* not DST_DBG */
    
#ifdef BXM
    ! Update netCDF file
    fl_out='aer.nc'           ! [sng] Name of netCDF output file
    call ftn_strnul(fl_out)
    call chm2nc(             &
         fl_out,              & ! I [sng] Name of netCDF output file
         lat_idx,             & ! I [idx] Model latitude index
         rxr_HNO3_gas_dst_vmr, & ! I [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
         rxrc_HNO3_dst,       & ! I [s-1] Pseudo first order rate coefficient for HNO3
         rxrc_HO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for HO2
         rxrc_N2O5_dst,       & ! I [s-1] Pseudo first order rate coefficient for N2O5
         rxrc_O3_dst,         & ! I [s-1] Pseudo first order rate coefficient for O3
         rxrc_SO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for SO2
         vmr_HNO3_gas,        & ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
         vmr_HO2_gas,         & ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
         vmr_N2O5_gas,        & ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
         vmr_NO3_aer,         & ! I [mlc mlc-1] Particulate NO3 volume mixing ratio
         vmr_O3_gas,          & ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
         vmr_SO2_gas,         & ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
         vmr_SO4_aer)         ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
#endif /* not BXM */
    
    return
  end subroutine dst_chm_slv
  
  subroutine chm2nc(             &
       fl_out,              & ! I [sng] Name of netCDF output file
       lat_idx,             & ! I [idx] Model latitude index
       rxr_HNO3_gas_dst_vmr, & ! I [mlc mlc-1 s-1] Mean rate of HNO3 removal by dust
       rxrc_HNO3_dst,       & ! I [s-1] Pseudo first order rate coefficient for HNO3
       rxrc_HO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for HO2
       rxrc_N2O5_dst,       & ! I [s-1] Pseudo first order rate coefficient for N2O5
       rxrc_O3_dst,         & ! I [s-1] Pseudo first order rate coefficient for O3
       rxrc_SO2_dst,        & ! I [s-1] Pseudo first order rate coefficient for SO2
       vmr_HNO3_gas,        & ! I [mlc mlc-1] Gaseous HNO3 volume mixing ratio
       vmr_HO2_gas,         & ! I [mlc mlc-1] Gaseous HO2 volume mixing ratio
       vmr_N2O5_gas,        & ! I [mlc mlc-1] Gaseous N2O5 volume mixing ratio
       vmr_NO3_aer,         & ! I [mlc mlc-1] Particulate NO3 volume mixing ratio
       vmr_O3_gas,          & ! I [mlc mlc-1] Gaseous O3 volume mixing ratio
       vmr_SO2_gas,         & ! I [mlc mlc-1] Gaseous SO2 volume mixing ratio
       vmr_SO4_aer          & ! I [mlc mlc-1] Particulate SO4 volume mixing ratio
       )                    
    ! Purpose: Output aerosol chemistry properties to netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use precision ! [mdl] Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstctl ! [mdl] Control variables, routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='chm2nc' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
    integer,intent(in)::lat_idx ! I [idx] Model latitude index
    real(r8),intent(in)::vmr_HNO3_gas(plond,plev) ! [s-1] Gaseous HNO3 volume mixing ratio
    real(r8),intent(in)::vmr_HO2_gas(plond,plev) ! [s-1] Gaseous HO2 volume mixing ratio
    real(r8),intent(in)::vmr_N2O5_gas(plond,plev) ! [s-1] Gaseous N2O5 volume mixing ratio
    real(r8),intent(in)::vmr_O3_gas(plond,plev) ! [s-1] Gaseous O3 volume mixing ratio
    real(r8),intent(in)::vmr_SO2_gas(plond,plev) ! [s-1] Gaseous SO2 volume mixing ratio
    real(r8),intent(in)::vmr_NO3_aer(plond,plev) ! [s-1] Particulate NO3 volume mixing ratio
    real(r8),intent(in)::vmr_SO4_aer(plond,plev) ! [s-1] Particulate SO4 volume mixing ratio
    real(r8),intent(in)::rxr_HNO3_gas_dst_vmr(plond,plev) ! [s-1] Mean rate of HNO3 removal by dust
    real(r8),intent(in)::rxrc_HNO3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HNO3
    real(r8),intent(in)::rxrc_HO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for HO2
    real(r8),intent(in)::rxrc_N2O5_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for N2O5
    real(r8),intent(in)::rxrc_O3_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for O3
    real(r8),intent(in)::rxrc_SO2_dst(plond,plev) ! [s-1] Pseudo first order rate coefficient for SO2
    ! Output
    ! Local
    ! File metadata and dimension IDs
    integer cnt_lon_lev_sz_time(4) ! Count array
    integer cnt_lon_lev_time(3) ! Count array
    integer dim_lon_lev_sz_time(4) ! [enm] Dimension IDs
    integer dim_lon_lev_time(3) ! [enm] Dimension IDs
    integer fll_mode_old      ! Old fill mode
    integer lev_dim_id        ! [enm] Dimension ID for lev
    integer lon_dim_id        ! [enm] Dimension ID for lon
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer srt_lon_lev_sz_time(4) ! Starting index array
    integer srt_lon_lev_time(3) ! Starting index array
    integer sz_dim_id         ! [enm] Dimension ID for sz
    integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
    integer time_dim_id       ! [enm] Dimension ID for time
    ! Variable IDs
    integer vmr_HNO3_gas_id   ! [enm] Variable ID
    integer vmr_HO2_gas_id    ! [enm] Variable ID
    integer vmr_N2O5_gas_id   ! [enm] Variable ID
    integer vmr_O3_gas_id     ! [enm] Variable ID
    integer vmr_SO2_gas_id    ! [enm] Variable ID
    integer vmr_NO3_aer_id    ! [enm] Variable ID
    integer vmr_SO4_aer_id    ! [enm] Variable ID
    integer rxr_HNO3_gas_dst_vmr_id ! [enm] Variable ID
    integer rxrc_HNO3_dst_id  ! [enm] Variable ID
    integer rxrc_HO2_dst_id   ! [enm] Variable ID
    integer rxrc_N2O5_dst_id  ! [enm] Variable ID
    integer rxrc_O3_dst_id    ! [enm] Variable ID
    integer rxrc_SO2_dst_id   ! [enm] Variable ID
    integer upt_cff_H2O2_dst_id ! [enm] Variable ID
    integer upt_cff_HNO3_dst_id ! [enm] Variable ID
    integer upt_cff_HO2_dst_id ! [enm] Variable ID
    integer upt_cff_N2O5_dst_id ! [enm] Variable ID
    integer upt_cff_NO3_dst_id ! [enm] Variable ID
    integer upt_cff_O3_dst_id ! [enm] Variable ID
    integer upt_cff_OH_dst_id ! [enm] Variable ID
    integer upt_cff_SO2_dst_id ! [enm] Variable ID
    ! Variable data
    ! Begin netCDF output routines
    rcd=nf90_noerr              ! NF90_STRERROR == 0
    rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
    rcd=rcd+nf90_redef(nc_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'sz',sz_dim_id)
    ! Add global attributes
    ! Define dimension IDs
    ! Assemble ID and count vectors for each multidimensional combination of dimensions
    dim_lon_lev_time=(/lon_dim_id,lev_dim_id,time_dim_id/)
    cnt_lon_lev_time=(/plon,plev,1/)
    srt_lon_lev_time=(/1,1,nstep/)
    
    dim_lon_lev_sz_time=(/lon_dim_id,lev_dim_id,sz_dim_id,time_dim_id/)
    cnt_lon_lev_sz_time=(/plon,plev,dst_nbr,1/)
    srt_lon_lev_sz_time=(/1,1,1,nstep/)
    
    if (nstep == 1) then
       ! Variable definitions
       rcd=rcd+nf90_def_var(nc_id,'vmr_HNO3_gas',nf90_float,dim_lon_lev_time,vmr_HNO3_gas_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_HO2_gas',nf90_float,dim_lon_lev_time,vmr_HO2_gas_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_N2O5_gas',nf90_float,dim_lon_lev_time,vmr_N2O5_gas_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_O3_gas',nf90_float,dim_lon_lev_time,vmr_O3_gas_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_SO2_gas',nf90_float,dim_lon_lev_time,vmr_SO2_gas_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_NO3_aer',nf90_float,dim_lon_lev_time,vmr_NO3_aer_id)
       rcd=rcd+nf90_def_var(nc_id,'vmr_SO4_aer',nf90_float,dim_lon_lev_time,vmr_SO4_aer_id)
       rcd=rcd+nf90_def_var(nc_id,'rxr_HNO3_gas_dst_vmr',nf90_float,dim_lon_lev_time,rxr_HNO3_gas_dst_vmr_id)
       rcd=rcd+nf90_def_var(nc_id,'rxrc_HNO3_dst',nf90_float,dim_lon_lev_time,rxrc_HNO3_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'rxrc_HO2_dst',nf90_float,dim_lon_lev_time,rxrc_HO2_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'rxrc_N2O5_dst',nf90_float,dim_lon_lev_time,rxrc_N2O5_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'rxrc_O3_dst',nf90_float,dim_lon_lev_time,rxrc_O3_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'rxrc_SO2_dst',nf90_float,dim_lon_lev_time,rxrc_SO2_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_H2O2_dst',nf90_float,upt_cff_H2O2_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_HNO3_dst',nf90_float,upt_cff_HNO3_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_HO2_dst',nf90_float,upt_cff_HO2_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_N2O5_dst',nf90_float,upt_cff_N2O5_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_NO3_dst',nf90_float,upt_cff_NO3_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_O3_dst',nf90_float,upt_cff_O3_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_OH_dst',nf90_float,upt_cff_OH_dst_id)
       rcd=rcd+nf90_def_var(nc_id,'upt_cff_SO2_dst',nf90_float,upt_cff_SO2_dst_id)
       ! Add english text descriptions
       rcd=rcd+nf90_put_att(nc_id,vmr_HNO3_gas_id,'long_name','Gaseous HNO3 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_HO2_gas_id,'long_name','Gaseous HO2 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_N2O5_gas_id,'long_name','Gaseous N2O5 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_O3_gas_id,'long_name','Gaseous O3 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_SO2_gas_id,'long_name','Gaseous SO2 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_NO3_aer_id,'long_name','Particulate NO3 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,vmr_SO4_aer_id,'long_name','Particulate SO4 volume mixing ratio')
       rcd=rcd+nf90_put_att(nc_id,rxr_HNO3_gas_dst_vmr_id,'long_name','Mean rate of HNO3 removal by dust')
       rcd=rcd+nf90_put_att(nc_id,rxrc_HNO3_dst_id,'long_name','Pseudo first order rate coefficient for HNO3')
       rcd=rcd+nf90_put_att(nc_id,rxrc_HO2_dst_id,'long_name','Pseudo first order rate coefficient for HO2')
       rcd=rcd+nf90_put_att(nc_id,rxrc_N2O5_dst_id,'long_name','Pseudo first order rate coefficient for N2O5')
       rcd=rcd+nf90_put_att(nc_id,rxrc_O3_dst_id,'long_name','Pseudo first order rate coefficient for O3')
       rcd=rcd+nf90_put_att(nc_id,rxrc_SO2_dst_id,'long_name','Pseudo first order rate coefficient for SO2')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_H2O2_dst_id,'long_name','Uptake coefficient for H2O2 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_HNO3_dst_id,'long_name','Uptake coefficient for HNO3 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_HO2_dst_id,'long_name','Uptake coefficient for HO2 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_N2O5_dst_id,'long_name','Uptake coefficient for N2O5 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_NO3_dst_id,'long_name','Uptake coefficient for NO3 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_O3_dst_id,'long_name','Uptake coefficient for O3 to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_OH_dst_id,'long_name','Uptake coefficient for OH to dust')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_SO2_dst_id,'long_name','Uptake coefficient for SO2 to dust')
       ! Add units
       rcd=rcd+nf90_put_att(nc_id,vmr_HNO3_gas_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_HO2_gas_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_N2O5_gas_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_O3_gas_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_SO2_gas_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_NO3_aer_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,vmr_SO4_aer_id,'units','molecule molecule-1')
       rcd=rcd+nf90_put_att(nc_id,rxr_HNO3_gas_dst_vmr_id,'units','molecule molecule-1 second-1')
       rcd=rcd+nf90_put_att(nc_id,rxrc_HNO3_dst_id,'units','second-1')
       rcd=rcd+nf90_put_att(nc_id,rxrc_HO2_dst_id,'units','second-1')
       rcd=rcd+nf90_put_att(nc_id,rxrc_N2O5_dst_id,'units','second-1')
       rcd=rcd+nf90_put_att(nc_id,rxrc_O3_dst_id,'units','second-1')
       rcd=rcd+nf90_put_att(nc_id,rxrc_SO2_dst_id,'units','second-1')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_H2O2_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_HNO3_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_HO2_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_N2O5_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_NO3_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_O3_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_OH_dst_id,'units','fraction')
       rcd=rcd+nf90_put_att(nc_id,upt_cff_SO2_dst_id,'units','fraction')
    else                      ! endif nstep == 1
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_HNO3_gas',vmr_HNO3_gas_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_HO2_gas',vmr_HO2_gas_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_N2O5_gas',vmr_N2O5_gas_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_O3_gas',vmr_O3_gas_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_SO2_gas',vmr_SO2_gas_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_NO3_aer',vmr_NO3_aer_id)
       rcd=nf90_wrp_inq_varid(nc_id,'vmr_SO4_aer',vmr_SO4_aer_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxr_HNO3_gas_dst_vmr',rxr_HNO3_gas_dst_vmr_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxrc_HNO3_dst',rxrc_HNO3_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxrc_HO2_dst',rxrc_HO2_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxrc_N2O5_dst',rxrc_N2O5_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxrc_O3_dst',rxrc_O3_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'rxrc_SO2_dst',rxrc_SO2_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_H2O2_dst',upt_cff_H2O2_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_HNO3_dst',upt_cff_HNO3_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_HO2_dst',upt_cff_HO2_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_N2O5_dst',upt_cff_N2O5_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_NO3_dst',upt_cff_NO3_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_O3_dst',upt_cff_O3_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_OH_dst',upt_cff_OH_dst_id)
       rcd=nf90_wrp_inq_varid(nc_id,'upt_cff_SO2_dst',upt_cff_SO2_dst_id)
    endif                     ! endif nstep /= 1
    ! End define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    if (nstep == 1) then
       rcd=rcd+nf90_put_var(nc_id,upt_cff_H2O2_dst_id,upt_cff_H2O2_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_HNO3_dst_id,upt_cff_HNO3_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_HO2_dst_id,upt_cff_HO2_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_N2O5_dst_id,upt_cff_N2O5_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_NO3_dst_id,upt_cff_NO3_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_O3_dst_id,upt_cff_O3_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_OH_dst_id,upt_cff_OH_dst)
       rcd=rcd+nf90_put_var(nc_id,upt_cff_SO2_dst_id,upt_cff_SO2_dst)
    endif                     ! endif nstep == 1
    rcd=rcd+nf90_put_var(nc_id,vmr_HNO3_gas_id,vmr_HNO3_gas,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_HO2_gas_id,vmr_HO2_gas,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_N2O5_gas_id,vmr_N2O5_gas,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_O3_gas_id,vmr_O3_gas,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_SO2_gas_id,vmr_SO2_gas,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_NO3_aer_id,vmr_NO3_aer,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,vmr_SO4_aer_id,vmr_SO4_aer,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxr_HNO3_gas_dst_vmr_id,rxr_HNO3_gas_dst_vmr,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxrc_HNO3_dst_id,rxrc_HNO3_dst,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxrc_HO2_dst_id,rxrc_HO2_dst,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxrc_N2O5_dst_id,rxrc_N2O5_dst,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxrc_O3_dst_id,rxrc_O3_dst,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    rcd=rcd+nf90_put_var(nc_id,rxrc_SO2_dst_id,rxrc_SO2_dst,start=srt_lon_lev_time,count=cnt_lon_lev_time)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
    if (nstep == 1) then
       write (6,'(a,a40,a)') prg_nm(1:ftn_strlen(prg_nm)),': Initialized chemistry data archive in ',fl_out(1:ftn_strlen(fl_out))
    endif                     ! endif
  end subroutine chm2nc
  
end module dstchm

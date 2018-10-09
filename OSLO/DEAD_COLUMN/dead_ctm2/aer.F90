! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/aer.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Offline driver for aerosol microphysics

! Copyright (C) 1998--2003 Charlie Zender

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! The author of this software, Charlie Zender, would like to receive
! your suggestions, improvements, bug-reports, and patches.
! Charlie Zender, zender@uci.edu
! Department of Earth System Science
! University of California at Irvine
! Irvine, CA 92697-3100

! Usage:
! aer --dbg=1 ! Turns on debugging output
! aer --time_nbr=10 ! Runs for 10 timesteps
#if ( !defined AIX ) && ( !defined CRAY )
! etags *.F90 *.h ~/c++/*.cc ~/c++/*.hh ~/mie/*.cc ~/mie/*.hh ~/map/*.F90 ~/f/*.F90
! */
#endif /* not CRAY */
! Compilation:
! Debugging             : cd ~/aer;make OPTS=D aer;cd - 
! Double precision reals: cd ~/aer;make OPTS=D PRC=D aer;cd -
! Single precision reals: cd ~/aer;make OPTS=D PRC=S aer;cd -
! Turn-on sandblasting  : cd ~/aer;make USR_TKN=-DAlG01 aer;cd -

#ifdef BXM
#include <params.h>
#include <dst.h> /* Dust preprocessor tokens */

program aer
  ! Purpose: Offline driver for aerosol microphysics
  use aernvr ! [mdl] Aerosol environmental properties
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use dstaer ! [mdl] Aerosol microphysical properties
  use dstbdg ! [mdl] Mass budget diagnostics
  use dstcmnini,only:dst_msc_cmn_ini ! [mdl] Module initialization
  use dstcst,only:dst_cst_cmn_ini ! [mdl] Physical constants for dust routines
  use dstctl ! [mdl] Control variables, routines
  use dstnm,only:dst_nm_cmn_ini ! [mdl] Nomenclature for outfld()
  use dstpsd,only:dst_psd_ini,dst_psd_src_ini ! [mdl] Dust particle size distributions
  use dstsfc,only:dst_sfc_set ! [mdl] 2-D surface fields on PLON x PLAT grid
  use phyzlic,only:phys_drv ! [mdl] Physics driver
  use pmgrid ! [mdl] Spatial resolution parameters
  use precision ! [mdl] Precision r8, i8, ...
  use sng_mdl ! [mdl] String manipulation
  use utl_mdl,only:date_time_get ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use xtr_mdl ! [mdl] Extrapolation constants
  implicit none
  ! Parameters
  character(len=*),parameter::CVS_Date='$Date: 2003/04/15 14:41:35 $' ! [sng] Date string
  character(len=*),parameter::CVS_Header='$Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/aer.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $' ! [sng] Full CVS Header
  character(len=*),parameter::CVS_Id='$Id: aer.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $' ! [sng] CVS Identification
  character(len=*),parameter::CVS_Name='$Name:  $' ! [sng] File name string
  character(len=*),parameter::CVS_Revision='$Revision: 1.1 $' ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  
  ! Physics parameters needed for call to dst_cst_cmn_ini()
  real(r8),parameter::gas_cst_unv=8.31441 ! [J mol-1 K-1] Universal gas constant
  real(r8),parameter::mmw_dry_air=28.9644e-3 ! [kg mol-1] (Source: radcsw.F in CCM2/3)
  real(r8),parameter::gas_cst_dry_air=gas_cst_unv/mmw_dry_air ! (287.05) [J kg-1 K-1] IrG81 p. 25, p. 245
  real(r8),parameter::grv_mean_sfc=9.80665 ! [m s-2] Mean gravitational acceleration at Earth's surface
  real(r8),parameter::rds_earth=6.370e+06 ! [m] Radius of sphere of same volume as Earth
  ! Commons
  ! Input
  ! Input/Output
  ! Output
  ! Local
  real(r8),dimension(:,:,:),allocatable::ext_dat_hst ! External data history (time_nbr,plond,ext_dat_nbr)
  
  ! Set defaults for command line options 
  character(80)::fl_in='in.nc'//nlc ! [sng] Name of netCDF input file
  character(80)::fl_out='aer.nc'//nlc ! [sng] Name of netCDF output file
  character(80)::drc_in='' ! [sng] Input directory
  character(80)::drc_out='' ! [sng] Output directory
  character(80)::fl_ext_dat='' ! [sng] Name of netCDF file with external forcing data (e.g. wind)
  integer::time_nbr=1 ! [nbr] Number of timesteps to simulate
  
  ! Derived fields
  
  ! Locals with simple initialization and no command line override
  integer::exit_status=0 ! [enm] Program exit status
  
  character(80)::arg_val ! [sng] Command line argument value
  character(80)::opt_sng ! [sng] Option string
  character(500)::cmd_ln ! [sng] Command line
  character(2)::dsh_key ! [sng] Command line dash and switch
  character(26)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(200)::prg_ID ! [sng] Program ID
  
  integer arg_nbr ! [nbr] Number of command line arguments
  integer arg_idx ! [idx] Counting index
  integer opt_lng ! [nbr] Length of option
  integer lat_idx ! [idx] Counting index
  integer::ext_dat_nbr=15 ! [nbr] Number of external forcing data variables
  real(r8) lat_wgt(plat) ! [frc] Latitude weights (currently must sum to 2.0)
  ! Main code
  ! Set defaults for command line arguments
  dbg_lvl=dbg_off ! dbg_lvl allocated in dbg.com, Option D
  ! tm_adj must be set before bdg_cmn_ini() is called
  tm_adj=2400.0_r8 ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
  ! [fnc] Initialize aer_nvr_cmn command-line fields
  call aer_nvr_cmn_cmd_ln_dfl() ! Not needed in global models
  
  ! Derived fields
  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command line into single string
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count() ! [nbr] Number of command line arguments
  arg_idx=1 ! [idx] Counting index
  do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
     dsh_key=arg_val(1:2) ! [sng] First two characters of option
     if (dsh_key == '--') then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop 'Long option has no name'
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (opt_sng == 'dbg' .or. opt_sng == 'dbg_lvl' ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == 'asp_rat_lps') then
           call ftn_arg_get(arg_idx,arg_val,asp_rat_lps) ! [frc] Ellipsoidal aspect ratio
        else if (opt_sng == 'drc_in') then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == 'drc_out') then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == 'fl_ext_dat') then
           call ftn_arg_get(arg_idx,arg_val,fl_ext_dat) ! [sng] Name of netCDF file with external forcing data
        else if (opt_sng == 'hgt_mdp') then
           call ftn_arg_get(arg_idx,arg_val,hgt_mdp(1,plev)) ! [m] Midlayer height above surface
        else if (opt_sng == 'oro') then
           call ftn_arg_get(arg_idx,arg_val,oro(1)) ! [frc] Orography
        else if (opt_sng == 'prs_mdp') then
           call ftn_arg_get(arg_idx,arg_val,prs_mdp(1,plev)) ! [Pa] Midlayer pressure
        else if (opt_sng == 'prs_ntf') then
           call ftn_arg_get(arg_idx,arg_val,prs_ntf(1,plevp)) ! [Pa] Interface pressure
        else if (opt_sng == 'q_H2O_vpr') then
           call ftn_arg_get(arg_idx,arg_val,q_H2O_vpr(1,plev)) ! [kg kg-1] Water vapor mixing ratio
        else if (opt_sng == 'sfc_typ') then
           call ftn_arg_get(arg_idx,arg_val,sfc_typ(1)) ! [idx] LSM surface type (0..28)
        else if (opt_sng == 'time_nbr') then
           call ftn_arg_get(arg_idx,arg_val,time_nbr) ! [nbr] Number of timesteps to simulate
        else if (opt_sng == 'tpt_gnd') then
           call ftn_arg_get(arg_idx,arg_val,tpt_gnd(1)) ! [K] Ground temperature
        else if (opt_sng == 'tpt_ice') then
           call ftn_arg_get(arg_idx,arg_val,tpt_ice(1)) ! [K] Ice temperature
        else if (opt_sng == 'tpt_mdp') then
           call ftn_arg_get(arg_idx,arg_val,tpt_mdp(1,plev)) ! [K] Temperature
        else if (opt_sng == 'tpt_soi') then
           call ftn_arg_get(arg_idx,arg_val,tpt_soi(1)) ! [K] Soil temperature
        else if (opt_sng == 'tpt_sst') then
           call ftn_arg_get(arg_idx,arg_val,tpt_sst(1)) ! [K] Sea surface temperature
        else if (opt_sng == 'vai_dst') then
           call ftn_arg_get(arg_idx,arg_val,vai_dst(1)) ! [m2 m-2] Vegetation area index, one-sided
        else if (opt_sng == 'wnd_mrd_mdp') then
           call ftn_arg_get(arg_idx,arg_val,wnd_mrd_mdp(1)) ! [m s-1] Meridional wind component
        else if (opt_sng == 'wnd_znl_mdp') then
           call ftn_arg_get(arg_idx,arg_val,wnd_znl_mdp(1)) ! [m s-1] Zonal wind component
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif option is recognized
        ! Jump to top of while loop
        cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
     endif ! endif long option
     ! Handle short options
     if (dsh_key == '-D') then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
     else if (dsh_key == '-d') then
        call ftn_arg_get(arg_idx,arg_val,tm_adj) ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
     else if (dsh_key == '-i') then
        call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
     else if (dsh_key == '-n') then
        call ftn_arg_get(arg_idx,arg_val,time_nbr) ! [nbr] Number of timesteps to simulate
     else if (dsh_key == '-o') then
        call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
     else if (dsh_key == '-v') then
        goto 1000 ! Goto exit with error status
     else ! Option not recognized
        arg_idx=arg_idx-1 ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif ! endif arg_val
  end do ! end while (arg_idx <= arg_nbr)
  
  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file
  
  ! Construct dummy Gaussian weights
  do lat_idx=1,plat
     lat_wgt(lat_idx)=2.0/plat ! [frc] Latitude weights (currently must sum to 2.0)
  end do ! end loop over lat
  ! Initialize nstep before using in dst_psd_ini()->aersz2nc()
  nstep=0 ! [nbr] 
  
  ! Initialize time-invariant physical constants
  ! dst_cst_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  call dst_cst_cmn_ini(gas_cst_dry_air,grv_mean_sfc,rds_Earth)
  ! Initialize dust names
  ! dst_nm_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  call dst_nm_cmn_ini()
  ! Initialize aerosol environmental properties common block aer_nvr_cmn
  call aer_nvr_cmn_ini() ! Not needed in global models
  ! Initialize size distributions at source
  ! dst_psd_src_ini is called by CCM:physics/inti(), MATCH:src/inirun
  call dst_psd_src_ini()
  ! Initialize size grid
  ! dst_psd_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  call dst_psd_ini()
  ! Initialize miscellaneous common blocks
  ! dst_msc_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  call dst_msc_cmn_ini()
#ifdef DST_MSS_BDG
  ! Initialize mass budget common block
  ! bdg_cmn_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
  call bdg_cmn_ini(grv_mean_sfc,lat_wgt,tm_adj)
#endif /* not DST_MSS_BDG */
  
  ! Box model method to set surface tvbds common block
  call dst_sfc_set(plond,   & ! Not needed in global models
       flx_LW_dwn_sfc,      & ! I [W m-2] Longwave downwelling flux at surface
       flx_SW_abs_sfc,      & ! I [W m-2] Solar flux absorbed by ground
       lnd_frc_dry,         & ! I [frc] Dry land fraction
       mbl_bsn_fct,         & ! I [frc] Erodibility factor
       mss_frc_CaCO3,       & ! I [frc] Mass fraction of CaCO3
       mss_frc_cly,         & ! I [frc] Mass fraction of clay
       mss_frc_snd,         & ! I [frc] Mass fraction of sand
       sfc_frc_bln,         & ! I [frc] Fraction of 4 available soil types
       sfc_typ,             & ! I [idx] LSM surface type (0..28)
       tpt_gnd,             & ! I [K] Ground temperature
       tpt_soi,             & ! I [K] Soil temperature
       vai_dst,             & ! I [m2 m-2] Vegetation area index, one-sided
       vwc_sfc)             ! I [m3 m-3] Volumetric water content

  ! Read external forcing data if specified
  if (ftn_strlen(fl_ext_dat) > 0) then
     allocate(ext_dat_hst(time_nbr,plond,ext_dat_nbr))
     call ext_dat_get(fl_ext_dat,time_nbr,plond,ext_dat_nbr,ext_dat_hst)
  endif
  
  ! Run physics on each timestep on each latitude slice
  do nstep=1,time_nbr
     if (dbg_lvl > dbg_off) write (6,'(a1)',advance="no") '.'

     ! Set time-varying external data
     if (ftn_strlen(fl_ext_dat) > 0) then
        !call ext_dat_set(time_nbr, plond, nstep, ext_dat_nbr, ext_dat_hst, 1, asp_rat_lps)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,2,hgt_mdp)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,3,oro)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,4,prs_mdp)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,5,prs_ntf)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,6,q_H2O_vpr)
        call ext_dat_set_i(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,7,sfc_typ)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,8,tpt_gnd)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,9,tpt_ice)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,10,tpt_mdp)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,11,tpt_soi)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,12,tpt_sst)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,13,vai_dst)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,14,wnd_mrd_mdp)
        call ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,15,wnd_znl_mdp)
     endif

     ! Simplest possible OpenMP parallelization
     !$omp parallel do
     do lat_idx=1,plat
        call phys_drv(lat_idx)
     end do                 ! end loop over lat
     !$omp end parallel do
#ifdef DST_MSS_BDG
     ! Output global average mass budget for this timestep
     ! bdg_update() is called by CCM:dynamics/eul/stepon(), MATCH:main()
     call bdg_update( nstep, mcdate, mcsec )
#endif /* not DST_MSS_BDG */
     
     ! Output time-varying variables
     call tm2nc(            & ! Not needed in global models
          fl_out,           & ! I [sng] Name of netCDF output file
          tm_adj)           ! I [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
  end do                    ! end loop over time
  
#ifdef DST_MSS_BDG
  ! bdg_close() is called by CCM:dynamics/eul/stepon(), MATCH:main()
  call bdg_close()
#endif /* not DST_MSS_BDG */
  
  ! Add environmental variables to file
  call nvr2nc(              & ! Not needed in global models
       CVS_Id,              & ! [sng] CVS Identification
       cmd_ln,              & ! [sng] Command line
       fl_out,              & ! [sng] Name of netCDF output file
       lcl_date_time,       & ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
       prg_ID)              ! [sng] Program ID
  
  if (dbg_lvl == dbg_fl) then
     write (6,'(a)') 'ncks -C -F -u -H aer.nc | m'
     write (6,'(a)') 'ncks -C -F -u -H -v nstep,mcdate,mcsec,nbdate,nbsec,ndcur,nscur,time,tm_adj aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v dst_slt_flx_rat_ttl,flx_mss_hrz_slt_ttl,'// &
          'flx_mss_vrt_dst_ttl,lnd_frc_mbl,lnd_frc_dry,mbl_bsn_fct,sfc_typ,lat_dgr,doy aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v mno_lng_mbl,mno_lng_dps,'// &
          'rgh_mmn_mbl,rgh_mmn_dps,hgt_zpd_mbl,hgt_zpd_dps,'// &
          'wnd_frc_mbl,wnd_frc_slt,wnd_frc_dps,'// &
          'wnd_frc_thr_slt,wnd_rfr_mbl,wnd_rfr_dps aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v mpc_dst,mpc_dst_ttl,odxc_dst,odxc_dst_ttl,q_dst,q_dst_ttl,mss_cnc_dst aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v flx_mss_dry_sfc,flx_mss_dry_sfc_ttl,'// &
          'flx_mss_grv_sfc,flx_mss_grv_sfc_ttl,flx_mss_trb_sfc,flx_mss_trb_sfc_ttl,'// &
          'q_dst_tnd_dry,q_dst_tnd_dry_ttl aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v flx_mss_pcp_sfc,flx_mss_pcp_sfc_ttl,'// &
          'spc_xsx_ncl_scv,spc_xsx_cll_scv,frc_trc_trc_cnv_ptn,pcp_flx_sfc,'// &
          'q_dst_tnd_ncl,q_dst_tnd_wet,q_dst_tnd_evp,q_dst_tnd_evp_ttl,q_dst_tnd_wet,q_dst_tnd_pcp_ttl aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v cld_frc,cld_frc_cnv,'// &
          'cld_vlm,q_H2O_cnd,q_H2O_cnd_cnv,q_H2O_cnd2pcp_tnd,q_H2O_cnd_pcp,'// &
          'q_H2O_cnd_tnd,q_H2O_pcp2vpr_tnd,q_H2O_vpr2pcp_cnv_tnd aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -d sz,2 -v mss_frc_cly,flx_mss_hrz_slt_ttl,'// &
          'flx_mss_vrt_dst_ttl,dst_slt_flx_rat_ttl,ovr_src_snk_mss,'// &
          'ovr_src_snk_frc,ovr_src_snk_mss_ttl,flx_mss_vrt_dst aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -v dns_mdp,gwc_sfc,mss_frc_cly,mss_frc_snd,'// &
          'wnd_frc_mbl,wnd_frc_slt,wnd_frc_thr_slt,wnd_rfr_mbl,wnd_mdp,'// &
          'wnd_rfr_thr_slt,frc_thr_ncr_drg,frc_thr_ncr_wtr,'// &
          'rgh_mmn_mbl,snw_frc,snw_hgt_lqd,vai_dst,vwc_sfc,hgt_zpd_mbl,'// &
          'hgt_mdp aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -v dns_mdp,mno_lng_mbl,prs_mdp,prs_ntf,q_H2O_vpr,'// &
          'sfc_typ,tpt_mdp,tpt_vrt,wnd_frc_mbl,wnd_frc_slt,wnd_mdp,hgt_mdp,'// &
          'rgh_mmn_mbl,hgt_zpd_mbl,wnd_rfr_mbl aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -v rss_lmn,shm_nbr,stk_nbr,vlc_grv,vlc_trb aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -d sz,2 -v nbr_spc_rsl,sfc_spc_rsl,vlm_spc_rsl aer.nc'
     write (6,'(a)')  &
          'ncks -C -H -F -u -d sz,2 -v ovr_src_snk_mss aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -d sz,2 -v dmt_nwr,dmt_swr,dmt_vwr,'// &
          'dmt_nmr,dmt_smr,dmt_vmr aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v '// &
          'rxrc_HNO3_dst,upt_cff_HNO3_dst,vmr_HNO3_gas,'// &
          'rxrc_N2O5_dst,upt_cff_N2O5_dst,vmr_N2O5_gas,'// &
          'rxrc_SO2_dst,upt_cff_SO2_dst,vmr_SO2_gas aer.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v tau_dps,tau_mbl,tau_pcp,tau_dry,'// &
          'tm_ttl,mpc_dst_mdl_avg,'// &
          'tau_grv,tau_trb dst_mss_bdg.nc'
     write (6,'(a)') 'ncks -C -H -F -u -v mpc_err,mpc_err_frc,'// &
          'mpc_err_tnd_frc dst_mss_bdg.nc'
  endif                     ! endif dbg
1000 continue
  call exit(exit_status)
end program aer                       ! end aer()

subroutine nvr2nc(             &
     CVS_Id,              & ! [sng] CVS Identification
     cmd_ln,              & ! [sng] Command line
     fl_out,              & ! [sng] Name of netCDF output file
     lcl_date_time,       & ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
     prg_ID)              ! [sng] Program ID
  ! Purpose: Output environmental properties to netCDF file
  use aernvr ! [mdl] Aerosol environmental properties
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use dstctl ! [mdl] Control variables, routines
  use dstgrd ! [mdl] Dust grid sizes
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use pmgrid ! [mdl] Spatial resolution parameters
  use precision ! [mdl] Precision r8, i8, ...
  use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
  implicit none
  ! Parameters
  character(len=*),parameter::sbr_nm="nvr2nc" ! [sng] Subroutine name
  ! Input
  character(len=*),intent(in)::CVS_Id ! [sng] CVS Identification
  character(len=*),intent(in)::cmd_ln ! [sng] Command line
  character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
  character(len=*),intent(in)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(len=*),intent(in)::prg_ID ! [sng] Program ID
  ! Output
  ! Local
  ! File metadata and dimension IDs
  integer dim_lon_lev_sz(3) ! [enm] Dimension IDs
  integer dim_lon_lev(2)    ! [enm] Dimension IDs
  integer dim_lon_levp(2)   ! [enm] Dimension IDs
  integer dim_lon_sz(2)     ! [enm] Dimension IDs
  integer dim_sz_src_sz(2)  ! [enm] Dimension IDs
  integer lev_dim_id        ! [enm] Dimension ID for lev
  integer levp_dim_id       ! [enm] Dimension ID for lev
  integer lon_dim_id        ! [enm] Dimension ID for lon
  integer lat_dim_id        ! [enm] Dimension ID for lat
  integer nc_id             ! File handle
  integer rcd               ! [rcd] Return success code
  integer sz_dim_id         ! [enm] Dimension ID for sz
  integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
  integer sz_src_dim_id     ! [enm] Dimension ID for sz_src
  ! Variable IDs
  integer dns_mdp_id        ! [enm] Variable ID
  integer doy_id            ! [enm] Variable ID
  integer flx_LW_dwn_sfc_id ! [enm] Variable ID
  integer flx_SW_abs_sfc_id ! [enm] Variable ID
  integer hgt_mdp_id        ! [enm] Variable ID
  integer lat_dgr_id        ! [enm] Variable ID
  integer lat_id            ! Coordinate ID
  integer lat_rdn_id        ! [enm] Variable ID
  integer lev_id            ! Coordinate ID
  integer levp_id           ! Coordinate ID
  integer lnd_frc_dry_id    ! [enm] Variable ID
  integer mbl_bsn_fct_id    ! [enm] Variable ID
  integer lon_id            ! Coordinate ID
  integer mpl_air_id        ! [enm] Variable ID
  integer mss_cnc_dst_id    ! [enm] Variable ID
  integer mss_frc_CaCO3_id  ! [enm] Variable ID
  integer mss_frc_cly_id    ! [enm] Variable ID
  integer mss_frc_slt_id    ! [enm] Variable ID
  integer mss_frc_snd_id    ! [enm] Variable ID
  integer oro_id            ! [enm] Variable ID
  integer prs_dlt_id        ! [enm] Variable ID
  integer prs_mdp_id        ! [enm] Variable ID
  integer prs_ntf_id        ! [enm] Variable ID
  integer q_H2O_cnd_pcp_id  ! [enm] Variable ID
  integer q_H2O_cnd_tnd_id  ! [enm] Variable ID
  integer q_H2O_pcp_lqd_id  ! [enm] Variable ID
  integer q_H2O_vpr_id      ! [enm] Variable ID
  integer sfc_typ_id        ! [enm] Variable ID
  integer snw_hgt_lqd_id    ! [enm] Variable ID
  integer tm_adj_id         ! [enm] Variable ID
  integer tpt_gnd_id        ! [enm] Variable ID
  integer tpt_mdp_id        ! [enm] Variable ID
  integer tpt_ptn_id        ! [enm] Variable ID
  integer tpt_sfc_id        ! [enm] Variable ID
  integer tpt_soi_id        ! [enm] Variable ID
  integer tpt_vrt_id        ! [enm] Variable ID
  integer vai_dst_id        ! [enm] Variable ID
  integer vwc_sfc_id        ! [enm] Variable ID
  integer wnd_mrd_mdp_id    ! [enm] Variable ID
  integer wnd_znl_mdp_id    ! [enm] Variable ID
  ! Variable data
  
  ! Begin netCDF output routines 
  rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
  rcd=nf90_wrp(nf90_redef(nc_id),sbr_nm//': nf90_redef')
  rcd=nf90_wrp_inq_dimid(nc_id,'levp',levp_dim_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dim_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lon',lon_dim_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'lat',lat_dim_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'sz',sz_dim_id)
  rcd=nf90_wrp_inq_dimid(nc_id,'sz_src',sz_src_dim_id)
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': put_att CVS_Id')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm//': put_att creation_date')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm//': put_att prg_ID')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm//': put_att cmd_ln')
  ! Define dimension IDs
  ! Assemble ID and count vectors for each multidimensional combination of dimensions
  dim_lon_sz=(/lon_dim_id,sz_dim_id/)
  dim_lon_lev=(/lon_dim_id,lev_dim_id/)
  dim_lon_levp=(/lon_dim_id,levp_dim_id/)
  dim_lon_lev_sz=(/lon_dim_id,lev_dim_id,sz_dim_id/)
  dim_sz_src_sz=(/sz_src_dim_id,sz_dim_id/)
  ! Variable definitions
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_mdp',nf90_float,dim_lon_lev,dns_mdp_id),sbr_nm//': def_var dns_mdp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'doy',nf90_float,doy_id),sbr_nm//': def_var doy')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_LW_dwn_sfc',nf90_float,lon_dim_id,flx_LW_dwn_sfc_id),sbr_nm//': def_var flx_LW_dwn_sfc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_SW_abs_sfc',nf90_float,lon_dim_id,flx_SW_abs_sfc_id),sbr_nm//': def_var flx_SW_abs_sfc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'hgt_mdp',nf90_float,dim_lon_lev,hgt_mdp_id),sbr_nm//': def_var hgt_mdp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat',nf90_float,lat_dim_id,lat_id),sbr_nm//': def_var lat')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat_dgr',nf90_float,lat_dim_id,lat_dgr_id),sbr_nm//': def_var lat_dgr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat_rdn',nf90_float,lat_dim_id,lat_rdn_id),sbr_nm//': def_var lat_rdn')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev',nf90_float,lev_dim_id,lev_id),sbr_nm//': def_var lev')
  rcd=nf90_wrp(nf90_def_var(nc_id,'levp',nf90_float,levp_dim_id,levp_id),sbr_nm//': def_var levp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lnd_frc_dry',nf90_float,lon_dim_id,lnd_frc_dry_id),sbr_nm//': def_var lnd_frc_dry')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mbl_bsn_fct',nf90_float,lon_dim_id,mbl_bsn_fct_id),sbr_nm//': def_var mbl_bsn_fct')
  rcd=nf90_wrp(nf90_def_var(nc_id,'lon',nf90_float,lon_dim_id,lon_id),sbr_nm//': def_var lon')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_air',nf90_float,dim_lon_lev,mpl_air_id),sbr_nm//': def_var mpl_air')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mss_cnc_dst',nf90_float,dim_lon_lev,mss_cnc_dst_id),sbr_nm//': def_var mss_cnc_dst')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mss_frc_CaCO3',nf90_float,lon_dim_id,mss_frc_CaCO3_id),sbr_nm//': def_var mss_frc_CaCO3')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mss_frc_cly',nf90_float,lon_dim_id,mss_frc_cly_id),sbr_nm//': def_var mss_frc_cly')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mss_frc_slt',nf90_float,lon_dim_id,mss_frc_slt_id),sbr_nm//': def_var mss_frc_slt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mss_frc_snd',nf90_float,lon_dim_id,mss_frc_snd_id),sbr_nm//': def_var mss_frc_snd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'oro',nf90_float,lon_dim_id,oro_id),sbr_nm//': def_var oro')
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_dlt',nf90_float,dim_lon_lev,prs_dlt_id),sbr_nm//': def_var prs_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_mdp',nf90_float,dim_lon_lev,prs_mdp_id),sbr_nm//': def_var prs_mdp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_ntf',nf90_float,dim_lon_levp,prs_ntf_id),sbr_nm//': def_var prs_ntf')
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2O_cnd_pcp',nf90_float,dim_lon_lev,q_H2O_cnd_pcp_id),sbr_nm//': def_var q_H2O_cnd_pcp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2O_cnd_tnd',nf90_float,dim_lon_lev,q_H2O_cnd_tnd_id),sbr_nm//': def_var q_H2O_cnd_tnd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2O_pcp_lqd',nf90_float,dim_lon_lev,q_H2O_pcp_lqd_id),sbr_nm//': def_var q_H2O_pcp_lqd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2O_vpr',nf90_float,dim_lon_lev,q_H2O_vpr_id),sbr_nm//': def_var q_H2O_vpr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'sfc_typ',nf90_int,lon_dim_id,sfc_typ_id),sbr_nm//': def_var sfc_typ')
  rcd=nf90_wrp(nf90_def_var(nc_id,'snw_hgt_lqd',nf90_float,lon_dim_id,snw_hgt_lqd_id),sbr_nm//': def_var snw_hgt_lqd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tm_adj',nf90_float,tm_adj_id),sbr_nm//': def_var tm_adj')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_gnd',nf90_float,lon_dim_id,tpt_gnd_id),sbr_nm//': def_var tpt_gnd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_mdp',nf90_float,dim_lon_lev,tpt_mdp_id),sbr_nm//': def_var tpt_mdp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ptn',nf90_float,dim_lon_lev,tpt_ptn_id),sbr_nm//': def_var tpt_ptn')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_sfc',nf90_float,lon_dim_id,tpt_sfc_id),sbr_nm//': def_var tpt_sfc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_soi',nf90_float,lon_dim_id,tpt_soi_id),sbr_nm//': def_var tpt_soi')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_vrt',nf90_float,dim_lon_lev,tpt_vrt_id),sbr_nm//': def_var tpt_vrt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'vwc_sfc',nf90_float,lon_dim_id,vwc_sfc_id),sbr_nm//': def_var vwc_sfc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wnd_mrd_mdp',nf90_float,lon_dim_id,wnd_mrd_mdp_id),sbr_nm//': def_var wnd_mrd_mdp')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wnd_znl_mdp',nf90_float,lon_dim_id,wnd_znl_mdp_id),sbr_nm//': def_var wnd_znl_mdp')
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_mdp_id,'long_name','Midlayer density'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,doy_id,'long_name','Day of year [1.0..366.0)'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_LW_dwn_sfc_id,'long_name','Longwave downwelling flux at surface'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_SW_abs_sfc_id,'long_name','Solar flux absorbed by ground'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,hgt_mdp_id,'long_name','Midlayer height above surface'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'long_name','Latitude'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'long_name','Latitude'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_rdn_id,'long_name','Latitude'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'long_name','Midlayer pressure'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'long_name','Interface pressure'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lnd_frc_dry_id,'long_name','Dry land fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mbl_bsn_fct_id,'long_name','Erodibility factor'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'long_name','Longitude'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_air_id,'long_name','Air mass path in layer'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_cnc_dst_id,'long_name','Mass concentration of dust'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_CaCO3_id,'long_name','Mass fraction CaCO3'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_cly_id,'long_name','Mass fraction clay'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_slt_id,'long_name','Mass fraction silt'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_snd_id,'long_name','Mass fraction sant'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,oro_id,'long_name','Orography: ocean=0.0, land=1.0, sea ice=2.0'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_dlt_id,'long_name','Pressure thickness'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_mdp_id,'long_name','Midlayer pressure'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_ntf_id,'long_name','Interface pressure'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_cnd_pcp_id,'long_name','H2O precipitation mixing ratio'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_cnd_tnd_id,'long_name','Net H2O condensate formation tendency'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_pcp_lqd_id,'long_name','Rain water mixing ratio'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_vpr_id,'long_name','Water vapor mixing ratio'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,sfc_typ_id,'long_name','LSM surface type (0..28)'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,snw_hgt_lqd_id,'long_name','Equivalent liquid water snow depth'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tm_adj_id,'long_name','Adjustment timestep (CCM: 2*dt, MATCH: dt'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_gnd_id,'long_name','Ground temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_mdp_id,'long_name','Midlayer temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_id,'long_name','Potential temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_sfc_id,'long_name','Surface temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_soi_id,'long_name','Soil temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_vrt_id,'long_name','Virtual temperature'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,vwc_sfc_id,'long_name','Volumetric water content'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,wnd_mrd_mdp_id,'long_name','Meridional wind component'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,wnd_znl_mdp_id,'long_name','Zonal wind component'),sbr_nm//': put_att')
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_mdp_id,'units','kilogram meter-3'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,doy_id,'units','day'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_LW_dwn_sfc_id,'units','watt meter-2'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_SW_abs_sfc_id,'units','watt meter-2'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,hgt_mdp_id,'units','meter'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'units','degrees north'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'units','radians north'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_rdn_id,'units','radians north'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'units','Pascal'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'units','Pascal'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lnd_frc_dry_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mbl_bsn_fct_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'units','degrees east'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_air_id,'units','kilogram meter-2'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_cnc_dst_id,'units','kilogram meter-3'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_CaCO3_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_cly_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_slt_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,mss_frc_snd_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,oro_id,'units','fraction'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_dlt_id,'units','pascal'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_mdp_id,'units','pascal'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_ntf_id,'units','pascal'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_cnd_pcp_id,'units','kilogram kilogram-1'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_cnd_tnd_id,'units','kilogram kilogram-1 second-1'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_pcp_lqd_id,'units','kilogram kilogram-1'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_vpr_id,'units','kilogram kilogram-1'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,sfc_typ_id,'units','index'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,snw_hgt_lqd_id,'units','meter'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tm_adj_id,'units','second'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_gnd_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_mdp_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_sfc_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_soi_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_vrt_id,'units','kelvin'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,vwc_sfc_id,'units','meter3 meter-3'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,wnd_mrd_mdp_id,'units','meter second-1'),sbr_nm//': put_att')
  rcd=nf90_wrp(nf90_put_att(nc_id,wnd_znl_mdp_id,'units','meter second-1'),sbr_nm//': put_att')
  ! End define mode
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//'enddef')
  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_mdp_id,dns_mdp),sbr_nm//': put_var dns_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,doy_id,doy),sbr_nm//': put_var doy')
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_LW_dwn_sfc_id,flx_LW_dwn_sfc),sbr_nm//': put_var flx_LW_dwn_sfc')
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_SW_abs_sfc_id,flx_SW_abs_sfc),sbr_nm//': put_var flx_SW_abs_sfc')
  rcd=nf90_wrp(nf90_put_var(nc_id,hgt_mdp_id,hgt_mdp),sbr_nm//': put_var hgt_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_dgr_id,lat_dgr),sbr_nm//': put_var lat_dgr')
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_id,lat),sbr_nm//': put_var lat')
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_rdn_id,lat_rdn),sbr_nm//': put_var lat_rdn')
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_id,lev),sbr_nm//': put_var lev')
  rcd=nf90_wrp(nf90_put_var(nc_id,levp_id,levp),sbr_nm//': put_var levp')
  rcd=nf90_wrp(nf90_put_var(nc_id,lnd_frc_dry_id,lnd_frc_dry),sbr_nm//': put_var lnd_frc_dry')
  rcd=nf90_wrp(nf90_put_var(nc_id,mbl_bsn_fct_id,mbl_bsn_fct),sbr_nm//': put_var mbl_bsn_fct')
  rcd=nf90_wrp(nf90_put_var(nc_id,lon_id,lon),sbr_nm//': put_var lon')
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_air_id,mpl_air),sbr_nm//': put_var mpl_air')
  rcd=nf90_wrp(nf90_put_var(nc_id,mss_cnc_dst_id,mss_cnc_dst),sbr_nm//': put_var mss_cnc_dst')
  rcd=nf90_wrp(nf90_put_var(nc_id,mss_frc_CaCO3_id,mss_frc_CaCO3),sbr_nm//': put_var mss_frc_CaCO3')
  rcd=nf90_wrp(nf90_put_var(nc_id,mss_frc_cly_id,mss_frc_cly),sbr_nm//': put_var mss_frc_cly')
  rcd=nf90_wrp(nf90_put_var(nc_id,mss_frc_slt_id,mss_frc_slt),sbr_nm//': put_var mss_frc_slt')
  rcd=nf90_wrp(nf90_put_var(nc_id,mss_frc_snd_id,mss_frc_snd),sbr_nm//': put_var mss_frc_snd')
  rcd=nf90_wrp(nf90_put_var(nc_id,oro_id,oro),sbr_nm//': put_var oro')
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_dlt_id,prs_dlt),sbr_nm//': put_var prs_dlt')
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_mdp_id,prs_mdp),sbr_nm//': put_var prs_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_ntf_id,prs_ntf),sbr_nm//': put_var prs_ntf')
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2O_cnd_pcp_id,q_H2O_cnd_pcp),sbr_nm//': put_var q_H2O_cnd_pcp')
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2O_cnd_tnd_id,q_H2O_cnd_tnd),sbr_nm//': put_var q_H2O_cnd_tnd')
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2O_pcp_lqd_id,q_H2O_pcp_lqd),sbr_nm//': put_var q_H2O_pcp_lqd')
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2O_vpr_id,q_H2O_vpr),sbr_nm//': put_var q_H2O_vpr')
  rcd=nf90_wrp(nf90_put_var(nc_id,snw_hgt_lqd_id,snw_hgt_lqd),sbr_nm//': put_var snw_hgt_lqd')
  rcd=nf90_wrp(nf90_put_var(nc_id,tm_adj_id,tm_adj),sbr_nm//': put_var tm_adj')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_gnd_id,tpt_gnd),sbr_nm//': put_var tpt_gnd')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_mdp_id,tpt_mdp),sbr_nm//': put_var tpt_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ptn_id,tpt_ptn),sbr_nm//': put_var tpt_ptn')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_sfc_id,tpt_sfc),sbr_nm//': put_var tpt_sfc')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_soi_id,tpt_soi),sbr_nm//': put_var tpt_soi')
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_vrt_id,tpt_vrt),sbr_nm//': put_var tpt_vrt')
  rcd=nf90_wrp(nf90_put_var(nc_id,vwc_sfc_id,vwc_sfc),sbr_nm//': put_var vwc_sfc')
  rcd=nf90_wrp(nf90_put_var(nc_id,wnd_mrd_mdp_id,wnd_mrd_mdp),sbr_nm//': put_var wnd_mrd_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,wnd_znl_mdp_id,wnd_znl_mdp),sbr_nm//': put_var wnd_znl_mdp')
  rcd=nf90_wrp(nf90_put_var(nc_id,sfc_typ_id,sfc_typ),sbr_nm//': put_var sfc_typ')
  ! Close output file
  rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
  return 
end subroutine nvr2nc ! end nvr2nc()

subroutine tm2nc(             &
     fl_out,              & ! [sng] Name of netCDF output file
     tm_adj)              ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
  ! Purpose: Output a single time step to netCDF file
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use precision ! [mdl] Precision r8, i8, ...
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use dstctl ! [mdl] Control variables, routines
  use sng_mdl,only:ftn_strlen,ftn_strnul,ftn_date2sng,ftn_sec2sng ! [mdl] String manipulation
  implicit none
  ! Parameters
  character(len=*),parameter::sbr_nm='tm2nc' ! [sng] Subroutine name
  ! Input
  character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
  real(r8),intent(in)::tm_adj   ! [s] Adjustment timestep (CCM: 2*dt, MATCH: dt)
  ! Output
  ! Local
  character(31) time_unit_sng ! [sng] Units attribute for time coordinate
  real(DBLKIND) time     ! [day] Days since simulation start
  ! File metadata and dimension IDs
  integer time_dim_id       ! [enm] Dimension ID for time
  integer nc_id             ! File handle
  integer rcd               ! [rcd] Return success code
  integer srt(1)            ! [idx] Starting offsets
  ! Variable IDs
  integer nstep_id          ! [enm] Variable ID
  integer mcdate_id         ! [enm] Variable ID
  integer mcsec_id          ! [enm] Variable ID
  integer nbdate_id         ! [enm] Variable ID
  integer nbsec_id          ! [enm] Variable ID
  integer ndcur_id          ! [enm] Variable ID
  integer nscur_id          ! [enm] Variable ID
  integer time_id           ! [enm] Variable ID
  ! Variable data
 
  srt(1)=nstep              ! [idx] Starting offsets
  
  ! Must be careful that variable "time" is defined in routine without access to F77 intrinsic time() function
  if (nstep == 1) then
     nbdate=640312          ! [day] Simulation start date in YYMMDD format
     nbsec=0                ! [s] Simulation start second relative to nbdate
     mcdate=nbdate          ! Current date in YYMMDD format
     mcsec=nbsec            ! [s] Seconds past current date at 0Z
     ndcur=0                ! [day] Current day number of simulation[day] 
     nscur=0                ! [s] Seconds relative to ndcur
     time_unit_sng(:)=char(0) ! [sng] Units attribute for time coordinate
     ! ftn_date2sng ! [sng] Convert YYYYMMDD integer to YYYY-MM-DD string
     ! ftn_sec2sng  ! [sng] Convert integer to 'HH:MM:SS' string
     time_unit_sng='days since '//ftn_date2sng(mcdate)//' '//ftn_sec2sng(mcsec)
  endif                     ! endif nstep == 1
  
  ! Increment time coordinates
  mcsec=mcsec+tm_adj        ! [s] Seconds past current date at 0Z
  nscur=nscur+tm_adj        ! [s] Seconds relative to ndcur
  
  ! Take care of odometer issues
  if (nscur >= 86400.0_r8) then
     ndcur=ndcur+nscur/86400 ! [day] Current day number of simulation
     nscur=mod(nscur,86400) ! [s] Seconds relative to ndcur
  end if                    ! endif nscur >= 86400.0
  if (mcsec >= 86400.0_r8) then
     mcdate=mcdate+mcsec/86400 ! [day] Current date in YYMMDD format
     mcsec=mod(mcsec,86400) ! [s] Seconds past current date at 0Z
  end if                    ! endif nscur >= 86400.0
  ! Construct time from ndcur and nscur
  time=ndcur+nscur/86400.0_r8  ! [day] Days since simulation start
  
  ! Begin netCDF output routines
  rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
  rcd=nf90_wrp(nf90_redef(nc_id),sbr_nm//': redef')
  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
  if (nstep == 1) then
     ! Variable definitions
     rcd=nf90_wrp(nf90_def_var(nc_id,'nstep',nf90_int,time_dim_id,nstep_id),sbr_nm//': def_var nstep')
     rcd=nf90_wrp(nf90_def_var(nc_id,'mcdate',nf90_int,time_dim_id,mcdate_id),sbr_nm//': def_var mcdate')
     rcd=nf90_wrp(nf90_def_var(nc_id,'mcsec',nf90_int,time_dim_id,mcsec_id),sbr_nm//': def_var mcsec')
     rcd=nf90_wrp(nf90_def_var(nc_id,'nbdate',nf90_int,nbdate_id),sbr_nm//': def_var nbdate')
     rcd=nf90_wrp(nf90_def_var(nc_id,'nbsec',nf90_int,nbsec_id),sbr_nm//': def_var nbsec')
     rcd=nf90_wrp(nf90_def_var(nc_id,'ndcur',nf90_int,time_dim_id,ndcur_id),sbr_nm//': def_var ndcur')
     rcd=nf90_wrp(nf90_def_var(nc_id,'nscur',nf90_int,time_dim_id,nscur_id),sbr_nm//': def_var nscur')
     rcd=nf90_wrp(nf90_def_var(nc_id,'time',nf90_float,time_dim_id,time_id),sbr_nm//': def_var time')
     ! Add english text descriptions
     rcd=nf90_wrp(nf90_put_att(nc_id,nstep_id,'long_name','Timestep'),sbr_nm//': put_att nstep')
     rcd=nf90_wrp(nf90_put_att(nc_id,mcdate_id,'long_name','Current date in YYMMDD format'),sbr_nm//': put_att mcdate')
     rcd=nf90_wrp(nf90_put_att(nc_id,mcsec_id,'long_name','Seconds past current date at 0Z'),sbr_nm//': put_att mcsec')
     rcd=nf90_wrp(nf90_put_att(nc_id,nbdate_id,'long_name','Simulation start date in YYMMDD format'),sbr_nm//': put_att nbdate')
     rcd=nf90_wrp(nf90_put_att(nc_id,nbsec_id,'long_name','Simulation start second relative to nbdate'),sbr_nm//': put_att nbsec')
     rcd=nf90_wrp(nf90_put_att(nc_id,ndcur_id,'long_name','Current day number of simulation'),sbr_nm//': put_att ndcur')
     rcd=nf90_wrp(nf90_put_att(nc_id,nscur_id,'long_name','Seconds relative to ndcur'),sbr_nm//': put_att nscur')
     rcd=nf90_wrp(nf90_put_att(nc_id,time_id,'long_name','Days since simulation start'),sbr_nm//': put_att time')
     ! Add units
     rcd=nf90_wrp(nf90_put_att(nc_id,nstep_id,'units','index'),sbr_nm//': put_att nstep')
     rcd=nf90_wrp(nf90_put_att(nc_id,mcdate_id,'units','day'),sbr_nm//': put_att mcdate')
     rcd=nf90_wrp(nf90_put_att(nc_id,mcsec_id,'units','second'),sbr_nm//': put_att mcsec')
     rcd=nf90_wrp(nf90_put_att(nc_id,nbdate_id,'units','day'),sbr_nm//': put_att nbdate')
     rcd=nf90_wrp(nf90_put_att(nc_id,nbsec_id,'units','second'),sbr_nm//': put_att nbsec')
     rcd=nf90_wrp(nf90_put_att(nc_id,ndcur_id,'units','day'),sbr_nm//': put_att ndcur')
     rcd=nf90_wrp(nf90_put_att(nc_id,nscur_id,'units','second'),sbr_nm//': put_att nscur')
     rcd=nf90_wrp(nf90_put_att(nc_id,time_id,'units',time_unit_sng),sbr_nm//': put_att time')
  else                      ! endif nstep == 1
     ! Get variable IDs
     rcd=nf90_wrp_inq_varid(nc_id,'nstep',nstep_id)
     rcd=nf90_wrp_inq_varid(nc_id,'mcdate',mcdate_id)
     rcd=nf90_wrp_inq_varid(nc_id,'mcsec',mcsec_id)
     rcd=nf90_wrp_inq_varid(nc_id,'nbdate',nbdate_id)
     rcd=nf90_wrp_inq_varid(nc_id,'nbsec',nbsec_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ndcur',ndcur_id)
     rcd=nf90_wrp_inq_varid(nc_id,'nscur',nscur_id)
     rcd=nf90_wrp_inq_varid(nc_id,'time',time_id)
  endif                     ! endif nstep /= 1
  ! End define mode
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef') 
  ! Write data
  if (nstep == 1) then
     rcd=nf90_wrp(nf90_put_var(nc_id,nbdate_id,nbdate),sbr_nm//': put_var nbdate')
     rcd=nf90_wrp(nf90_put_var(nc_id,nbsec_id,nbsec),sbr_nm//': put_var nbsec')
  endif                     ! endif nstep == 1
  rcd=nf90_wrp(nf90_put_var(nc_id,nstep_id,nstep,start=srt),sbr_nm//': put_var nstep')
  rcd=nf90_wrp(nf90_put_var(nc_id,mcdate_id,mcdate,start=srt),sbr_nm//': put_var mcdate')
  rcd=nf90_wrp(nf90_put_var(nc_id,mcsec_id,mcsec,start=srt),sbr_nm//': put_var mcsec')
  rcd=nf90_wrp(nf90_put_var(nc_id,ndcur_id,ndcur,start=srt),sbr_nm//': put_var ndcur')
  rcd=nf90_wrp(nf90_put_var(nc_id,nscur_id,nscur,start=srt),sbr_nm//': put_var nscur')
  rcd=nf90_wrp(nf90_put_var(nc_id,time_id,time,start=srt),sbr_nm//': put_var time')
  ! Close output file
  rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
  return 
end subroutine tm2nc            ! end tm2nc()

! Read external forcing data (currently just zonal and meridional winds)
subroutine ext_dat_get(fl_ext_dat,time_nbr,plond,ext_dat_nbr,ext_dat_hst)
  ! Purpose: set time-varying zonal or meridional wind speed
  use precision ! [mdl] Precision r8, i8, ...
  use netcdf    ! [mdl] netCDF interface
  use nf90_utl  ! [mdl] netCDF utilities
  use utl_mdl   ! [mdl] Utility functions (date_time_get,mnt_chk...)
  implicit none
  character(len=*),intent(in)::fl_ext_dat ! [sng] File name
  integer,intent(in)::time_nbr
  integer,intent(in)::plond
  integer,intent(in)::ext_dat_nbr
  real(r8),intent(out)::ext_dat_hst(time_nbr,plond,ext_dat_nbr) ! O External forcing data
  
  integer nc_id             ! [id] File ID
  integer rcd               ! [rcd] Return success code
  integer time_dmn_id       ! [enm] Dimension ID for time
  integer time_id           ! [enm] Variable ID
  integer time_nbr_in       ! [nbr] Dimension size
  integer ext_dat_hst_id(15)
  real(r8) time(time_nbr)   ! [day] Time coordinate (day of year)
  integer idx

  ! open netCDF file
  rcd=0
  rcd=rcd+nf90_wrp_open(fl_ext_dat,nf90_nowrite,nc_id,"ext_dat_get")

  ! Get dimension IDs
  rcd=rcd+nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
  rcd=rcd+nf90_wrp(nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr_in),"inquire_dim time")
  if (time_nbr > time_nbr_in) stop "time_nbr > time_nbr_in"
  if (ext_dat_nbr /= 15) stop "ERROR: ext_dat_nbr /= 15"

  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,"time",      time_id)
  rcd=nf90_wrp_inq_varid(nc_id,"asp_rat_lps",ext_dat_hst_id(1))
  rcd=nf90_wrp_inq_varid(nc_id,"hgt_mdp",   ext_dat_hst_id(2))
  rcd=nf90_wrp_inq_varid(nc_id,"oro",       ext_dat_hst_id(3))
  rcd=nf90_wrp_inq_varid(nc_id,"prs_mdp",   ext_dat_hst_id(4))
  rcd=nf90_wrp_inq_varid(nc_id,"prs_ntf",   ext_dat_hst_id(5))
  rcd=nf90_wrp_inq_varid(nc_id,"q_H2O_vpr", ext_dat_hst_id(6))
  rcd=nf90_wrp_inq_varid(nc_id,"sfc_typ",   ext_dat_hst_id(7))
  rcd=nf90_wrp_inq_varid(nc_id,"tpt_gnd",   ext_dat_hst_id(8))
  rcd=nf90_wrp_inq_varid(nc_id,"tpt_ice",   ext_dat_hst_id(9))
  rcd=nf90_wrp_inq_varid(nc_id,"tpt_mdp",   ext_dat_hst_id(10))
  rcd=nf90_wrp_inq_varid(nc_id,"tpt_soi",   ext_dat_hst_id(11))
  rcd=nf90_wrp_inq_varid(nc_id,"tpt_sst",   ext_dat_hst_id(12))
  rcd=nf90_wrp_inq_varid(nc_id,"vai_dst",   ext_dat_hst_id(13))
  rcd=nf90_wrp_inq_varid(nc_id,"wnd_mrd_mdp",ext_dat_hst_id(14))
  rcd=nf90_wrp_inq_varid(nc_id,"wnd_znl_mdp",ext_dat_hst_id(15))

  ! Get time and wind data
  rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time),"get_var time")
  
  do idx=1,ext_dat_nbr
     rcd=nf90_wrp(nf90_get_var(nc_id,ext_dat_hst_id(idx),ext_dat_hst(:,1,idx)),"get_var idx")
  end do
  
  ! Verify time coordinate is monotonic increasing
  if(.not.mnt_ncr_chk(time,time_nbr)) stop "time coordinate not monotonic increasing"
  if (rcd > 0) write(6,'(a)') "WARNING: rcd > 0 in ext_dat_get"

end subroutine ext_dat_get

! Read external forcing data (currently just zonal and meridional winds)
subroutine ext_dat_get_eureka(fl_ext_dat,time_nbr,plond,ext_dat_nbr,ext_dat_hst)
  ! Purpose: set time-varying zonal or meridional wind speed
  use precision ! [mdl] Precision r8, i8, ...
  use netcdf    ! [mdl] netCDF interface
  use nf90_utl  ! [mdl] netCDF utilities
  use utl_mdl   ! [mdl] Utility functions (date_time_get,mnt_chk...)
  implicit none
  character(len=*),intent(in)::fl_ext_dat ! [sng] File name
  integer,intent(in)::time_nbr
  integer,intent(in)::plond
  integer,intent(in)::ext_dat_nbr
  real(r8),intent(out)::ext_dat_hst(time_nbr,plond,ext_dat_nbr) ! O External forcing data
  
  integer nc_id             ! [id] File ID
  integer rcd               ! [rcd] Return success code
  integer time_dmn_id       ! [enm] Dimension ID for time
  integer time_id           ! [enm] Variable ID
  integer time_nbr_in       ! [nbr] Dimension size
  integer wnd_mrd_mdp_hst_id
  integer wnd_znl_mdp_hst_id
  real(r8) time(time_nbr)   ! [day] Time coordinate (day of year)

  ! Open netCDF file
  rcd=nf90_wrp_open(fl_ext_dat,nf90_nowrite,nc_id,"ext_dat_get")

  ! Get dimension IDs
  rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
  rcd=nf90_wrp(nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr_in),"inquire_dim time")
  if (time_nbr > time_nbr_in) stop "time_nbr > time_nbr_in"
  
  ! Get variable IDs
  rcd=nf90_wrp_inq_varid(nc_id,"time",time_id)
  rcd=nf90_wrp_inq_varid(nc_id,"wnd_mrd_mdp_hst",wnd_mrd_mdp_hst_id)
  rcd=nf90_wrp_inq_varid(nc_id,"wnd_znl_mdp_hst",wnd_znl_mdp_hst_id)

  ! Get time and wind data
  rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time),"get_var time")
  rcd=nf90_wrp(nf90_get_var(nc_id,wnd_mrd_mdp_hst_id,ext_dat_hst(:,1,1)),"get_var wnd_mrd_mdp_hst")
  rcd=nf90_wrp(nf90_get_var(nc_id,wnd_znl_mdp_hst_id,ext_dat_hst(:,1,2)),"get_var wnd_znl_mdp_hst")
  
  ! Verify time coordinate is monotonic increasing
  if(.not.mnt_ncr_chk(time,time_nbr)) stop "time coordinate not monotonic increasing"
  print *,"Ingested ",fl_ext_dat
  print *,"time_nbr     = ",time_nbr
  print *,"max mrd wind = ",maxval(ext_dat_hst(:,:,1))
  print *,"min mrd wind = ",minval(ext_dat_hst(:,:,1))
  print *,"max znl wind = ",maxval(ext_dat_hst(:,:,2))
  print *,"min znl wind = ",minval(ext_dat_hst(:,:,2))
end subroutine ext_dat_get_eureka

! Set time-varying forcing data
subroutine ext_dat_set(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,ext_dat_idx,ext_dat_now)
  ! Purpose: set time-varying zonal or meridional wind speed
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  integer,intent(in)::time_nbr
  integer,intent(in)::plond
  integer,intent(in)::nstep
  integer,intent(in)::ext_dat_nbr
  real(r8),intent(in)::ext_dat_hst(time_nbr,plond,ext_dat_nbr) ! I External forcing data
  integer,intent(in)::ext_dat_idx
  real(r8),intent(out)::ext_dat_now(plond) ! Instantaneous forcing data
  ext_dat_now = ext_dat_hst(nstep,:,ext_dat_idx)
end subroutine ext_dat_set

! Set time-varying forcing data of integer type
subroutine ext_dat_set_i(time_nbr,plond,nstep,ext_dat_nbr,ext_dat_hst,ext_dat_idx,ext_dat_now)
  ! Purpose: set time-varying zonal or meridional wind speed
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  integer,intent(in)::time_nbr
  integer,intent(in)::plond
  integer,intent(in)::nstep
  integer,intent(in)::ext_dat_nbr
  real(r8),intent(in)::ext_dat_hst(time_nbr,plond,ext_dat_nbr) ! I External forcing data
  integer,intent(in)::ext_dat_idx
  integer,intent(out)::ext_dat_now(plond) ! Instantaneous forcing data
  ext_dat_now = nint(ext_dat_hst(nstep,:,ext_dat_idx))
end subroutine ext_dat_set_i

#ifdef CCM
subroutine outfld(name,data,pcols,lchnk)
  ! Purpose: Dummy routine for history tape write (CCM)
  use precision ! [mdl] Precision r8, i8, ...
  character(len=*),intent(in)::name ! [sng] Variable name
  integer,intent(in)::pcols ! [nbr] Number of columns
  integer,intent(in)::lchnk ! [idx] Chunk index
  real(r8),intent(in)::data(pcols,*) ! [frc] Data
  return 
end subroutine outfld            ! end outfld()
#else /* not CCM */
subroutine outfld(name,data,plond,lat_idx,obuf)
  ! Purpose: Dummy routine for history tape write (MATCH)
  use precision ! [mdl] Precision r8, i8, ...
  character(len=*),intent(in)::name ! [sng] Variable name
  integer,intent(in)::plond ! [nbr] Number of longitudes
  integer,intent(in)::lat_idx ! [idx] Latitude index
  real(r8),intent(in)::data(plond,*) ! [frc] Data
  real(r8),intent(inout):obuf(*) ! [bfr] Output buffer
  return 
end subroutine outfld            ! end outfld()
#endif /* not CCM */
#else /* not BXM */
! The following is just to allow the compiler to compile something
! when it attempts to compile this routine in nested mode.
! This routine is not needed or used when Dust is run stand-alone
! But, if the user tells the compiler to compile this code without
! anything here -- the compiler complains.  A more sophisticated
! build mechanism could decide which files should and should not
! be built, but the simple solution is to do as follows...
subroutine aer_stb
  stop 'aer_stb() error: BXM not set, this routine should not be called'
end subroutine aer_stb            ! end aer_stb()
#endif  /* not BXM */

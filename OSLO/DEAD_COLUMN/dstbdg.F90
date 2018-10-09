! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstbdg.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Mass budget diagnostic subroutines for the dust model
! Most routines are based on Brian Eaton's massbgt.F

! Usage:
! use dstbdg ! [mdl] Mass budget diagnostics

! params.h required for DST_MSS_BDG
!#include <params_dust.h>
!#include <dst.h> /* Dust preprocessor tokens */

module dstbdg ! [mdl] Mass budget diagnostics
!  use cmn_dust !// CTM3 params_dust is now module 20150213
!  use oc_dust, only: PLAT, !// CTM3 params_dust is now module 20150213
  use dead_precision ! [mdl] Precision r8, i8, ...
  use pmgrid,only:plat ! [mdl] Spatial resolution parameters
  implicit none
  save ! [stt] Changes to common variables are sticky
  private ! [stt] Symbols are private unless individually qualified as public
  public::bdg_cmn_ini ! [sbr] Initialize budget module
!  public::bdg_gam_dry ! [sbr] Add latitude of dry mass mixing ratio to 3D tracer budget
!  public::bdg_gam_wet ! [sbr] Add latitude of moist mass mixing ratio to 3D tracer budget
!  public::bdg_gam_wet_2d ! [sbr] Add latitude of moist mass mixing ratio to 2D tracer budget
!  public::bdg_aal ! [sbr] Add latitude to area average of 2D field
  public::bdg_update ! [sbr] Write budget terms to netCDF file
  public::bdg_close ! [sbr] Close mass budget netCDF file
!#ifndef BXM
!  public::qneg3_dgn ! [sbr] Custom replacement for qneg3()
!#endif /* BXM */
  
  ! Fields needed by mass diagnostic routines for dust parameterization
  ! These variables are initialized with routine dst_bdg_cmn_ini() 
  ! dst_bdg_cmn_ini() is called by BXM:aer(), CCM:physics/inti(), MATCH:src/inirun()
  integer,parameter::var_nbr=19 ! [nbr] Number of variables in mass budget file
  
  real(r8) grv_sfc_rcp          ! [s2 m-1] Reciprocal of gravity
  real(r8) grv_sfc_plon_rcp     ! [s2 m-1] 1.0/(grv_sfc*plon)
  real(r8) lat_wgt(plat)        ! [frc] Latitude weights (currently must sum to 2.0)
  real(r8) mss_lat(plat,var_nbr) ! [kg m-2] Mass path of tracer (modulo pre-factors)
  
  integer time_srt(1)  ! [idx] Record index for netCDF calls in subroutine bdg_update
  integer time_idx          ! [idx] Current record index
  integer nc_id             ! [id] File handle
  integer var_id(3+var_nbr) ! [enm] Variable IDs
  
  character(64) var_lst(var_nbr) ! [sng] Variable names
  character(80) fl_out       ! [sng] Output file
  
  ! tm_dlt must be carried in a common block since it is not passed each timestep
  real(r8) tm_dlt               ! [s] Timestep duration
  logical uninitialized     ! [flg] True until first bdg_update() call of run
  
  ! Variable IDs for derived (non-enumerated) variables are carried separately from var_id array
  integer dst_dry_dlt_id    ! [enm] Variable ID
  integer dst_mbl_dlt_id    ! [enm] Variable ID
  integer dst_pcp_dlt_id    ! [enm] Variable ID
  integer dst_sf_dps_id     ! [enm] Variable ID
  integer dst_sf_net_id     ! [enm] Variable ID
  integer flx_mss_dps_id    ! [enm] Variable ID
  integer flx_mss_dry_id    ! [enm] Variable ID
  integer flx_mss_grv_id    ! [enm] Variable ID
  integer flx_mss_mbl_id    ! [enm] Variable ID
  integer flx_mss_pcp_id    ! [enm] Variable ID
  integer flx_mss_trb_id    ! [enm] Variable ID
  integer mpc_dst_dgn_id    ! [enm] Variable ID
  integer mpc_dst_dlt_id    ! [enm] Variable ID
  integer mpc_dst_mdl_avg_id ! [enm] Variable ID
  integer mpc_err_id        ! [enm] Variable ID
  integer mpc_err_frc_id    ! [enm] Variable ID
  integer mpc_err_tnd_frc_id ! [enm] Variable ID
  integer tau_dps_id        ! [enm] Variable ID
  integer tau_dry_id        ! [enm] Variable ID
  integer tau_grv_id        ! [enm] Variable ID
  integer tau_mbl_id        ! [enm] Variable ID
  integer tau_pcp_id        ! [enm] Variable ID
  integer tau_trb_id        ! [enm] Variable ID
  integer tm_ttl_id         ! [enm] Variable ID
  
contains
  
  subroutine bdg_cmn_ini(grv_sfc,lat_wgt_tmp,tm_dlt_tmp)
    ! Initialize parameter-dependent fields in dust common block
    ! bdg_cmn_ini() is called by CCM:physics/inti(), MATCH:inirun()
    ! NB: MATCH reads delt as integer, CCM reads dtime as real
    ! Be sure calling routine passes a real to bdg_cmn_ini()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use pmgrid,only:plat,plon ! [mdl] Spatial resolution parameters
    use sng_mdl,only:ftn_strlen,ftn_strnul ! [mdl] String manipulation
    use vec_mdl,only:vec_set ! [mdl] Vector manipulation, interpolation, rebinning
    implicit none
    ! Input
    real(r8),intent(in)::grv_sfc  ! [m s-2] Gravity at surface of Earth
    real(r8),intent(in)::lat_wgt_tmp(plat) ! [frc] Latitude weights (currently must sum to 2.0)
    real(r8),intent(in)::tm_dlt_tmp ! [s] Timestep duration
    ! Local
    integer var_idx           ! [idx] Counting index
    integer lat_idx           ! [idx] Counting index
    integer rcd               ! [enm] Return success code
    integer time_dim_id       ! [enm] Dimension ID for time
    integer::nf90_r8 ! [enm] External netCDF type for r8 kind
    real(r8) lat_wgt_ttl ! [frc] Total of latitude weights (currently must sum to 2.0)
    ! Main code
    
    ! Initialize defaults
    uninitialized=.true.      ! [flg] True until first bdg_update() call of run
    
    ! Initializing the time index to 1 seems to leave record 0 blank in dst_mss_bdg.nc
    time_idx=1                ! [idx] Current record index
    
    rcd=nf90_noerr              ! [enm] nf90_noerr == 0
    fl_out='dst_mss_bdg.nc'   ! [sng] Output file
    
    ! Set common block variables from inputs
    grv_sfc_rcp=1.0_r8/grv_sfc   ! [s2 m-1] Reciprocal of gravity 
    tm_dlt=tm_dlt_tmp         ! [s] Timestep duration
    grv_sfc_plon_rcp=1.0_r8/(plon*grv_sfc) ! [s2 m-1] 1.0/(grv_sfc*plon)
    
    do lat_idx=1,plat
       ! Latitude weights should be (at least) double precision in production runs
       lat_wgt(lat_idx)=lat_wgt_tmp(lat_idx) ! [frc] Latitude weights (currently must sum to 2.0_r8)
       
    end do                    ! end loop over lat
    
    call vec_set(mss_lat,plat*var_nbr,0.0_r8) ! [kg m-2]
    
    ! Set up identifiers for mass calculations
    var_lst( 1)='dst_sf_dry'
    var_lst( 2)='dst_sf_mbl'
    var_lst( 3)='dst_sf_grv'
    var_lst( 4)='dst_sf_trb'
    var_lst( 5)='dst_sf_pcp'
    var_lst( 6)='dst_ss_dry'
    var_lst( 7)='dst_ss_pcp'
    var_lst( 8)='dst_ss_evp'
    var_lst( 9)='dst_ss_mbl'
    
    ! Which column mass path diagnostics are output is model-specific:
    ! MATCH currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_trn,mpc_dst_a_trn_qng,mpc_dst_a_dry,mpc_dst_a_mbl,mpc_dst_a_wet
    ! CCM currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_flt,mpc_dst_a_tphys
    ! BXM currently outputs: mpc_dst_mdl,mpc_dst_a_fxr,mpc_dst_a_flt,mpc_dst_a_tphys,mpc_dst_a_dry,mpc_dst_a_mbl,mpc_dst_a_wet
    var_lst(10)='mpc_dst_mdl' ! [kg m-2] Column mass path of dust at end of timestep 
    var_lst(11)='mpc_dst_a_fxr' ! [kg m-2] Column mass path of dust after mass fixer
    var_lst(12)='mpc_dst_a_flt' ! [kg m-2] Column mass path of dust after time filter
    var_lst(13)='mpc_dst_b_tphys' ! [kg m-2] Column mass path of dust after tphys()
    var_lst(14)='mpc_dst_a_trn' ! [kg m-2] Column mass path of dust after transport
    var_lst(15)='mpc_dst_a_trn_qng' ! [kg m-2] Column mass path of dust after qneg check for transport
    var_lst(16)='mpc_dst_a_dry' ! [kg m-2] Column mass path of dust after dry deposition
    var_lst(17)='mpc_dst_a_mbl' ! [kg m-2] Column mass path of dust after mobilization
    var_lst(18)='mpc_dst_a_wet' ! [kg m-2] Column mass path of dust after wet deposition
    var_lst(19)='mpc_H2O'     ! [kg m-2] Column mass path of H2O vapor
    
    ! Sanity checks
    if (var_nbr /= 19) stop 'dst: var_nbr /= 19 in bdg_cmn_ini()'
    if (tm_dlt <= 0.0_r8 .or. tm_dlt > 3600.0_r8) then
       write(6,'(a,f9.3,a2)') 'dst: bdg_cmn_ini() reports unreasonable timestep duration tm_dlt = ',tm_dlt,' s'
       stop
    endif                     ! endif err
    lat_wgt_ttl=0.0_r8           ! [frc] Total of latitude weights (currently must sum to 2.0)
    do lat_idx=1,plat
       lat_wgt_ttl=lat_wgt_ttl+lat_wgt(lat_idx) ! [frc] Total of latitude weights (currently must sum to 2.0)
    end do                    ! end loop over lat
    if (abs(2.0_r8-lat_wgt_ttl)/2.0_r8 > 1.0e-6) then
       write (6,'(a,f9.3,a)') 'dst: bdg_cmn_ini() reports lat_wgt_ttl = ',lat_wgt_ttl,' != 2.0'
       stop
    endif                     ! endif err
    
    ! Create netCDF output file
    call ftn_strnul(fl_out)
    rcd=rcd+nf90_create(fl_out,nf90_clobber,nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'Unable to open file '//fl_out)
    write(6,'(a,a28,1x,a)') prg_nm(1:ftn_strlen(prg_nm)),': Created netCDF output file',fl_out(1:ftn_strlen(fl_out))
    
    ! Define dimensions
    rcd=rcd+nf90_def_dim(nc_id,'time',nf90_unlimited,time_dim_id)
    ! Define variables
    nf90_r8=nf90_xtype_r8_get() ! [enm] External netCDF type for r8 kind
    rcd=rcd+nf90_def_var(nc_id,'nstep',nf90_int,time_dim_id,var_id(1))
    rcd=rcd+nf90_def_var(nc_id,'date',nf90_int,time_dim_id,var_id(2))
    rcd=rcd+nf90_def_var(nc_id,'datesec',nf90_int,time_dim_id,var_id(3))
    do var_idx=1,var_nbr
       rcd=rcd+nf90_def_var(nc_id,var_lst(var_idx),nf90_r8,time_dim_id,var_id(3+var_idx))
    end do                    ! end loop over var
    ! Derived budget fields defined in bdg_update()
    rcd=rcd+nf90_def_var(nc_id,'dst_dry_dlt',nf90_r8,time_dim_id,dst_dry_dlt_id)
    rcd=rcd+nf90_def_var(nc_id,'dst_mbl_dlt',nf90_r8,time_dim_id,dst_mbl_dlt_id)
    rcd=rcd+nf90_def_var(nc_id,'dst_pcp_dlt',nf90_r8,time_dim_id,dst_pcp_dlt_id)
    rcd=rcd+nf90_def_var(nc_id,'dst_sf_dps',nf90_r8,time_dim_id,dst_sf_dps_id)
    rcd=rcd+nf90_def_var(nc_id,'dst_sf_net',nf90_r8,time_dim_id,dst_sf_net_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_dps',nf90_r8,time_dim_id,flx_mss_dps_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_dry',nf90_r8,time_dim_id,flx_mss_dry_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_grv',nf90_r8,time_dim_id,flx_mss_grv_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_mbl',nf90_r8,time_dim_id,flx_mss_mbl_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_pcp',nf90_r8,time_dim_id,flx_mss_pcp_id)
    rcd=rcd+nf90_def_var(nc_id,'flx_mss_trb',nf90_r8,time_dim_id,flx_mss_trb_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_dst_dgn',nf90_r8,time_dim_id,mpc_dst_dgn_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_dst_dlt',nf90_r8,time_dim_id,mpc_dst_dlt_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_dst_mdl_avg',nf90_r8,time_dim_id,mpc_dst_mdl_avg_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_err',nf90_r8,time_dim_id,mpc_err_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_err_frc',nf90_r8,time_dim_id,mpc_err_frc_id)
    rcd=rcd+nf90_def_var(nc_id,'mpc_err_tnd_frc',nf90_r8,time_dim_id,mpc_err_tnd_frc_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_dps',nf90_r8,time_dim_id,tau_dps_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_dry',nf90_r8,time_dim_id,tau_dry_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_grv',nf90_r8,time_dim_id,tau_grv_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_mbl',nf90_r8,time_dim_id,tau_mbl_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_pcp',nf90_r8,time_dim_id,tau_pcp_id)
    rcd=rcd+nf90_def_var(nc_id,'tau_trb',nf90_r8,time_dim_id,tau_trb_id)
    rcd=rcd+nf90_def_var(nc_id,'tm_ttl',nf90_r8,time_dim_id,tm_ttl_id)
    
    ! Add descriptions
    ! rcd=rcd+nf90_put_att(nc_id,var_id(3+var_idx),'units','kilogram meter-2')
    ! End definitions
    rcd=rcd+nf90_enddef(nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'bdg_cmn_ini(): Unable to end define mode')
    return
  end subroutine bdg_cmn_ini                       ! end bdg_cmn_ini()
  
  integer function bdg_idx_get(var_nm)
    ! Purpose: Return index of var_nm string in var_lst
    implicit none
    ! Input
    character,intent(in)::var_nm*(*) ! [sng] Variable name
    ! Local
    integer var_idx           ! [idx] Counting index
    ! Main code
    bdg_idx_get=0
    do var_idx=1,var_nbr
       if (var_nm == var_lst(var_idx)) then
          bdg_idx_get=var_idx
          return
       end if                 ! endif
    end do                    ! end loop over var
    write(6,'(2a)') 'ERROR bdg_idx_get(): Variable name not found:',var_nm
    stop
  end function bdg_idx_get                       ! end bdg_idx_get()
  
  
  real(r8) function bdg_gam(var_nm)
    ! Purpose: Return global average mass [kg m-2] of a tracer field
    ! bdg_gam() is called by bdg_update() once for each 3D field each timestep
    ! Each latitude's contribution must be in memory prior to this call
    ! Routine usually used to convert 3D MMR field to global average MPC field
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    character,intent(in)::var_nm*(*) ! [sng] Variable name
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer var_idx           ! [idx] Counting index
    ! Main Code
    var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
    bdg_gam=0.0_r8
    do lat_idx=1,plat
       bdg_gam=bdg_gam+mss_lat(lat_idx,var_idx) ! [kg m-2]
    end do                    ! end loop over lat
    ! NB: 0.5 factor assumes lat_wgt sums to 2.0
    bdg_gam=bdg_gam*grv_sfc_plon_rcp*0.5_r8 ! [kg m-2]
    return
  end function bdg_gam                       ! end bdg_gam()


  real(r8) function bdg_aa(var_nm)
    ! Purpose: Return area average of field
    ! Convert 2D MPC to global average MPC
    ! bdg_aa() is called by bdg_update() once for each field each timestep
    ! Each latitude's contribution must be in memory prior to this call
    use pmgrid ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    character,intent(in)::var_nm*(*) ! [sng] Variable name
    ! Local
    integer var_idx           ! [idx] Variable index
    integer lat_idx
    ! Main code
    var_idx=bdg_idx_get(var_nm) ! [idx] Counting index
    bdg_aa=0.0_r8                ! [kg m-2]
    do lat_idx=1,plat
       bdg_aa=bdg_aa+mss_lat(lat_idx,var_idx) ! [kg m-2]
    end do                    ! end loop over lat
    ! 0.5 factor assumes lat_wgt sums to 2.0
    bdg_aa=bdg_aa*0.5_r8/plon    ! [kg m-2]
    return
  end function bdg_aa                       ! end bdg_aa()
  
  subroutine bdg_update( &
       nstep,               & ! I [idx] Timestep of host model
       date,                & ! I [YYMMDD] Current date in YYMMDD format
       datesec)             ! I [s] Seconds past current date at 0Z
    ! Purpose: Write global average budget terms to netCDF file
    ! bdg_update() is called by CCM:dynamics/eul/stepon(), MATCH:main()
    ! bdg_update() should be called each timestep
    ! A few diagnostic fields in bdg_update() have hysterisis and are incorrect after restarts
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    implicit none
    ! Parameters
    real(r8),parameter::mss_val_dbl=1.0e36 ! [frc] Missing value in double precision
    ! Input variables
    integer,intent(in)::nstep ! I [idx] Timestep of host model
    integer,intent(in)::date  ! I [YYMMDD] Current date in yymmdd format
    integer,intent(in)::datesec ! I [s] Seconds past current date at 0Z
    ! Local
    integer hdr_nbr                       ! [nbr] Number of header fields
    integer rcd                           ! [rcd] Return success code
    real(r8) wgt_pst  ! [frc] Weight of past in running average (fragile)
    real(r8) wgt_crr  ! [frc] Weight of current timestep in running average (fragile)
    ! Fields stored in mass budget arrays
    real(r8) mpc_dst_mdl ! [kg m-2] bdg_gam('DSTQ')
    real(r8) dst_sf_dry ! [kg m-2 s-1] bdg_aa('DSTSFDRY')
    real(r8) dst_sf_grv ! [kg m-2 s-1] bdg_aa('DSTSFGRV')
    real(r8) dst_sf_mbl ! [kg m-2 s-1] bdg_aa('DSTSFMBL')
    real(r8) dst_sf_pcp ! [kg m-2 s-1] bdg_aa('DSTSFPCP')
    real(r8) dst_sf_trb ! [kg m-2 s-1] bdg_aa('DSTSFTRB')
    real(r8) dst_ss_evp ! [kg m-2] bdg_gam('DSTSSEVP')
    real(r8) dst_ss_dry ! [kg m-2] bdg_gam('DSTSSDRY')
    real(r8) dst_ss_mbl ! [kg m-2] bdg_gam('DSTSSMBL')
    real(r8) dst_ss_pcp ! [kg m-2] bdg_gam('DSTSSPCP')
    ! Derived mass budget fields
    real(r8) dst_dry_dlt ! [kg m-2] dst_ss_dry-dst_sf_dry
    real(r8) dst_mbl_dlt ! [kg m-2] dst_ss_mbl-dst_sf_mbl
    real(r8) dst_pcp_dlt ! [kg m-2] dst_ss_pcp-dst_ss_evp-dst_sf_pcp
    real(r8) dst_sf_dps ! [kg m-2 s-1] dst_sf_dry+dst_sf_pcp (Positive definite) (fragile)
    real(r8) mpc_dst_dlt ! [kg m-2] Approximate error defined as mpc_dst_dgn-mpc_dst_mdl (fragile)

    real(r8) mpc_err  ! [kg m-2] Mass imbalance
    real(r8) mpc_err_frc ! [frc] Mass imbalance as fraction of mean column burden
    real(r8) mpc_err_tnd_frc ! [frc] Mean imbalance tendency as fraction of mean source tendency
    
    real(r8) tau_dps  ! [s] Timescale for deposition (fragile)
    real(r8) tau_dry  ! [s] Timescale for dry deposition (fragile)
    real(r8) tau_grv  ! [s] Timescale for gravitational deposition (fragile)
    real(r8) tau_mbl  ! [s] Timescale for mobilization (fragile)
    real(r8) tau_pcp  ! [s] Timescale for wet deposition (fragile)
    real(r8) tau_trb  ! [s] Timescale for turbulent deposition (fragile)

    ! Initialize fragile variables to give them save attribute
    real(r8),save::dst_sf_net=0.0 ! [kg m-2 s-1] dst_sf_mbl-dst_sf_dps (Positive into atmosphere) (fragile)
    real(r8),save::flx_mss_dps=0.0 ! [kg m-2 s-1] Running average deposition (fragile)
    real(r8),save::flx_mss_dry=0.0 ! [kg m-2 s-1] Running average dry deposition (fragile)
    real(r8),save::flx_mss_grv=0.0 ! [kg m-2 s-1] Running average gravitational deposition (fragile)
    real(r8),save::flx_mss_mbl=0.0 ! [kg m-2 s-1] Running average mobilization (fragile)
    real(r8),save::flx_mss_pcp=0.0 ! [kg m-2 s-1] Running average wet deposition (fragile)
    real(r8),save::flx_mss_trb=0.0 ! [kg m-2 s-1] Running average turbulent deposition (fragile)    
    real(r8),save::mpc_dst_dgn=0.0 ! [kg m-2] Time integral of dst_sf_net (fragile)
    real(r8),save::mpc_dst_mdl_avg=0.0 ! [kg m-2] Running average column mass path (fragile)
    real(r8),save::tm_ttl=0.0   ! [s] Elapsed time of simulation (fragile)

    ! Main Code
    
    ! Initialize defaults
    time_srt(1)=time_idx  ! [idx] Record index for netCDF calls in subroutine bdg_update
    hdr_nbr=3                 ! [nbr] Number of header fields
    rcd=nf90_noerr              ! nf90_noerr == 0
    
    ! Convert surface fluxes to global average surface fluxes
    dst_sf_mbl=bdg_aa('dst_sf_mbl') ! [kg m-2 s-1]
    dst_sf_grv=bdg_aa('dst_sf_grv') ! [kg m-2 s-1]
    dst_sf_dry=bdg_aa('dst_sf_dry') ! [kg m-2 s-1]
    dst_sf_trb=bdg_aa('dst_sf_trb') ! [kg m-2 s-1]
    dst_sf_pcp=bdg_aa('dst_sf_pcp') ! [kg m-2 s-1]
    mpc_dst_mdl=bdg_gam('mpc_dst_mdl') ! [kg m-2]
    
    ! Convert source/sink tendencies to global average surface fluxes
    dst_ss_mbl=bdg_gam('dst_ss_mbl') ! [kg m-2 s-1]
    dst_ss_dry=bdg_gam('dst_ss_dry') ! [kg m-2 s-1]
    dst_ss_pcp=bdg_gam('dst_ss_pcp') ! [kg m-2 s-1]
    dst_ss_evp=bdg_gam('dst_ss_evp') ! [kg m-2 s-1]
    
    ! Difference between diagnosed and modeled surface flux tendencies
    dst_mbl_dlt=dst_ss_mbl-dst_sf_mbl ! [kg m-2 s-1]
    dst_dry_dlt=dst_ss_dry-dst_sf_dry ! [kg m-2 s-1]
    dst_pcp_dlt=dst_ss_pcp-dst_ss_evp-dst_sf_pcp ! [kg m-2 s-1]
    
    ! Because of timesplitting, position of next two lines is somewhat model-dependent
    ! There seems to be no elegant way to code this to work satisfactorily on both MATCH and CCM
    ! This is largely due to hairy CCM timesplitting
    ! For CCM these lines should appear after mpc_dst_dlt computation so that this timestep's dst_sf_dps and dst_sf_net contribute to next timestep's mpc_dst_dgn and mpc_dst_dlt
    ! For MATCH (and box model) these lines should appear before mpc_dst_dlt computation so that this timestep's dst_sf_dps and dst_sf_net contribute to this timestep's mpc_dst_dgn and mpc_dst_dlt
    dst_sf_dps=dst_sf_dry+dst_sf_pcp ! [kg m-2 s-1] dst_sf_dry+dst_sf_pcp (Positive definite)
    dst_sf_net=dst_sf_mbl-dst_sf_dps ! [kg m-2 s-1] dst_sf_mbl-dst_sf_dps (Positive into atmosphere)
    
    ! In CCM, mpc_dst_mdl lags dst_mbl() and dst_sf_net by 1 timestep
    ! In MATCH, mpc_dst_mdl is in sync with dst_mbl() and dst_sf_net
    ! Using static local storage for mpc_dst_dgn is OK only for initial runs
    ! Variables with unpredicatble behavior on restart are marked "fragile"
    ! On restart runs, the diagnostics which rely on integration will break
    ! because they are not stored in the restart file and so are re-initialized
    ! to 0.0 in bdg_cmn_ini().
    ! Uninitialized flag is .true. only on first call to bdg_update()
    if (uninitialized) then
       ! Initialize running sums
       ! Uncomment following line for CCM (and move dst_sf_net definition after this block)
       ! mpc_dst_dgn=0.0_r8        ! [kg m-2] Time integral of dst_sf_net
       mpc_dst_dgn=tm_dlt*dst_sf_net ! [kg m-2] Time integral of dst_sf_net (fragile)
       mpc_dst_mdl_avg=0.0_r8    ! [kg m-2] Running average column mass path (fragile)
       tm_ttl=0.0_r8             ! [s] Elapsed time of simulation (fragile)
       uninitialized=.false.  ! [flg] True until first bdg_update() call of run
    else
       mpc_dst_dgn=mpc_dst_dgn+tm_dlt*dst_sf_net ! [kg m-2] Time integral of dst_sf_net (fragile)
    endif                     ! endif
    ! Difference between diagnosed and modeled atmospheric burden
    ! Positive when integral of dst_sf_net exceeds mpc_dst_mdl
    mpc_dst_dlt=mpc_dst_dgn-mpc_dst_mdl ! [kg m-2] Approximate error defined as mpc_dst_dgn-mpc_dst_mdl (fragile)
    
    ! Running averages
    ! Integration method employed here may accumulate significant roundoff error 
    tm_ttl=tm_ttl+tm_dlt      ! [s] Elapsed time of simulation (fragile)
    wgt_pst=(time_idx-1.0_r8)/time_idx ! [frc] Weight of past in running average (fragile)
    wgt_crr=1.0_r8/time_idx      ! [frc] Weight of current timestep in running average (fragile)
    mpc_dst_mdl_avg=wgt_pst*mpc_dst_mdl_avg+wgt_crr*mpc_dst_mdl ! [kg m-2] Running average column mass path (fragile)
    flx_mss_dps=wgt_pst*flx_mss_dps+wgt_crr*dst_sf_dps ! [kg m-2 s-1] Running average deposition (fragile)
    flx_mss_dry=wgt_pst*flx_mss_dry+wgt_crr*dst_sf_dry ! [kg m-2 s-1] Running average dry deposition (fragile)
    flx_mss_grv=wgt_pst*flx_mss_grv+wgt_crr*dst_sf_grv ! [kg m-2 s-1] Running average gravitational deposition (fragile)
    flx_mss_mbl=wgt_pst*flx_mss_mbl+wgt_crr*dst_sf_mbl ! [kg m-2 s-1] Running average mobilization (fragile)
    flx_mss_pcp=wgt_pst*flx_mss_pcp+wgt_crr*dst_sf_pcp ! [kg m-2 s-1] Running average wet deposition (fragile)
    flx_mss_trb=wgt_pst*flx_mss_trb+wgt_crr*dst_sf_trb ! [kg m-2 s-1] Running average turbulent deposition (fragile)
    
    ! Mass conservation
    ! Mass imbalance is difference between time integrated mean net surface flux and current atmospheric burden
    ! Mass imbalance is deviation from perfect mass conservation
    mpc_err=tm_ttl*(flx_mss_mbl-flx_mss_dps)-mpc_dst_mdl ! [kg m-2] Mass imbalance
    if (mpc_dst_mdl_avg /= 0.0_r8) then
       mpc_err_frc=mpc_err/mpc_dst_mdl_avg ! [frc] Mass imbalance as fraction of mean column burden
    else 
       mpc_err_frc=0.0_r8        ! [frc] Mass imbalance as fraction of mean column burden
    endif ! endif
    if (flx_mss_mbl /= 0.0_r8) then
       mpc_err_tnd_frc=mpc_err/(flx_mss_mbl*tm_ttl) ! [frc] Mean imbalance tendency as fraction of mean source tendency
    else 
       mpc_err_tnd_frc=0.0_r8 ! [frc] Mean imbalance tendency as fraction of mean source tendency
    endif ! endif
    
    ! Diagnostic timescales
    tau_dps=mss_val_dbl       ! [s] Timescale for deposition (fragile)
    tau_dry=mss_val_dbl       ! [s] Timescale for dry deposition (fragile)
    tau_grv=mss_val_dbl       ! [s] Timescale for gravitational deposition (fragile)
    tau_mbl=mss_val_dbl       ! [s] Timescale for mobilization (fragile)
    tau_pcp=mss_val_dbl       ! [s] Timescale for wet deposition (fragile)
    tau_trb=mss_val_dbl       ! [s] Timescale for turbulent deposition (fragile)
    if (flx_mss_dps /= 0.0_r8) tau_dps=mpc_dst_mdl_avg/flx_mss_dps ! [s] Timescale for deposition (fragile)
    if (flx_mss_dry /= 0.0_r8) tau_dry=mpc_dst_mdl_avg/flx_mss_dry ! [s] Timescale for dry deposition (fragile)
    if (flx_mss_grv /= 0.0_r8) tau_grv=mpc_dst_mdl_avg/flx_mss_grv ! [s] Timescale for gravitational deposition (fragile)
    if (flx_mss_mbl /= 0.0_r8) tau_mbl=mpc_dst_mdl_avg/flx_mss_mbl ! [s] Timescale for mobilization (fragile)
    if (flx_mss_pcp /= 0.0_r8) tau_pcp=mpc_dst_mdl_avg/flx_mss_pcp ! [s] Timescale for wet deposition (fragile)
    if (flx_mss_trb /= 0.0_r8) tau_trb=mpc_dst_mdl_avg/flx_mss_trb ! [s] Timescale for turbulent deposition (fragile)
    
    rcd=rcd+nf90_put_var(nc_id,var_id(1),nstep,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(2),date,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(3),datesec,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_mbl')),dst_sf_mbl,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_dry')),dst_sf_dry,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_grv')),dst_sf_grv,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_trb')),dst_sf_trb,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_sf_pcp')),dst_sf_pcp,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_mbl')),dst_ss_mbl,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_dry')),dst_ss_dry,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_pcp')),dst_ss_pcp,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('dst_ss_evp')),dst_ss_evp,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_mdl')),mpc_dst_mdl,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_fxr')),bdg_gam('mpc_dst_a_fxr'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_flt')),bdg_gam('mpc_dst_a_flt'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_b_tphys')),bdg_gam('mpc_dst_b_tphys'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_trn')),bdg_gam('mpc_dst_a_trn'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_trn_qng')),bdg_gam('mpc_dst_a_trn_qng'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_dry')),bdg_gam('mpc_dst_a_dry'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_mbl')),bdg_gam('mpc_dst_a_mbl'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_dst_a_wet')),bdg_gam('mpc_dst_a_wet'),start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,var_id(hdr_nbr+bdg_idx_get('mpc_H2O')),bdg_gam('mpc_H2O'),start=time_srt)
    
    ! Output derived diagnostic constants
    rcd=rcd+nf90_put_var(nc_id,dst_dry_dlt_id,dst_dry_dlt,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,dst_mbl_dlt_id,dst_mbl_dlt,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,dst_pcp_dlt_id,dst_pcp_dlt,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,dst_sf_dps_id,dst_sf_dps,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,dst_sf_net_id,dst_sf_net,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_dps_id,flx_mss_dps,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_dry_id,flx_mss_dry,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_grv_id,flx_mss_grv,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_mbl_id,flx_mss_mbl,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_pcp_id,flx_mss_pcp,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,flx_mss_trb_id,flx_mss_trb,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_dst_dgn_id,mpc_dst_dgn,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_dst_dlt_id,mpc_dst_dlt,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_dst_mdl_avg_id,mpc_dst_mdl_avg,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_err_id,mpc_err,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_err_frc_id,mpc_err_frc,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,mpc_err_tnd_frc_id,mpc_err_tnd_frc,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_dps_id,tau_dps,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_dry_id,tau_dry,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_grv_id,tau_grv,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_mbl_id,tau_mbl,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_pcp_id,tau_pcp,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tau_trb_id,tau_trb,start=time_srt)
    rcd=rcd+nf90_put_var(nc_id,tm_ttl_id,tm_ttl,start=time_srt)
    
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'bdg_update()')
    
    ! Update record index
    time_idx=time_idx+1       ! [idx] Current record index
    
    return
  end subroutine bdg_update                       ! end bdg_update()
  
  subroutine bdg_close()
    ! Purpose: Close mass budget netCDF file
    ! bdg_close() is called by CCM:stepon(), MATCH:main()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen,ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Local
    integer rcd               ! [rcd] Return success code
    ! Main code
    ! Initialize defaults
    rcd=nf90_noerr              ! [enm] nf90_noerr == 0
    rcd=rcd+nf90_close(nc_id)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'Unable to close file '//fl_out)
    call ftn_strnul(fl_out)
    write (6,'(a,a20,1x,a)') prg_nm(1:ftn_strlen(prg_nm)),': Closed netCDF file',fl_out(1:ftn_strlen(fl_out))
    return
  end subroutine bdg_close                       ! end bdg_close()
  
  
end module dstbdg

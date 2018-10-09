! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstgrd.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! dstgrd.F90 is similar in function and usage to pmgrid.F90 
! dstgrd.F90 should be declared after dst.F90 but ahead of other dust headers

! Usage:
! use dstgrd ! [mdl] Dust grid sizes

! params.h needed for DST_NBR and DST_IDX_SRT
!#include <params_dust.h>

module dstgrd             ! [mdl] Dust grid sizes
  use pmgrid, only: dst_idx_srt, dst_nbr
  use dead_precision             ! [mdl] Precision r8, i8, ...
  implicit none
  save                      ! [stt] Changes to common variables are sticky
  
  ! Size grid info
!//CTM3: These are defined in cmn_dust instead of params_dust.h
!  integer,parameter::dst_idx_srt=DST_IDX_SRT ! Constituent index of smallest dust size
!  integer,parameter::dst_nbr=DST_NBR ! [nbr] Number of dust constituents
  integer,parameter::dst_idx_end=dst_idx_srt+dst_nbr-1 ! Constituent index of largest dust size
  
  ! Fixed size grid info
  integer,parameter::dst_src_nbr=3 ! [nbr] Number of size distributions in source soil
  integer,parameter::sgs_nbr=200 ! [nbr] Number of sub-gridscale bins in large bin

  !alf ++ 20020405
  ! Dimensions of look up tables for dust emissions
  integer,parameter::wnd_frc_nbr=100 ![nbr] number of u* for which we know dust size distr
  integer,parameter::bln_nbr=4 ![nbr] number of natural soils size distributions available
  !alf -- 20020405

  ! Chemistry grid
  integer,parameter::chm_nbr=10 ! [nbr] Number of chemical species
  
  ! Shortwave spectral grid
  integer,parameter::bnd_nbr_SW=19 ! [nbr] Number of CCM SW spectral bands (must be same as radcsw.F:nspint)
  ! NB: Use of dst_odxc_chn_idx is currently deprecated in favor of narrow band optical depths simulating, e.g., AVHRR channels
  integer,parameter::dst_odxc_chn_idx=8 ! CCM SW band index for optical depth output
  
  ! GSFC Longwave spectral grid
  integer,parameter::bnd_nbr_GSFC_LW=10 ! [nbr] Number of GSFC LW spectral bands
  
  ! CCM Longwave spectral grid
  integer,parameter::bnd_nbr_LW=6 ! [nbr] Number of CCM LW spectral bands
  integer,parameter::idx_LW_0000_0800=1 ! [idx] Index of dust cnt. abs. 0000--0800 cm-1
  integer,parameter::idx_LW_0500_0650=2 ! [idx] Index of dust cnt. abs. 0500--0650 cm-1
  integer,parameter::idx_LW_0650_0800=3 ! [idx] Index of dust cnt. abs. 0650--0800 cm-1
  integer,parameter::idx_LW_0800_1000=4 ! [idx] Index of dust cnt. abs. 0800--1000 cm-1
  integer,parameter::idx_LW_1000_1200=5 ! [idx] Index of dust cnt. abs. 1000--1200 cm-1
  integer,parameter::idx_LW_1200_2000=6 ! [idx] Index of dust cnt. abs. 1200--2000 cm-1
  ! Combination of sub-bands
  integer,parameter::idx_LW_0500_0800=3 ! [idx] Index of dust cnt. abs. 0500--0800 cm-1
  integer,parameter::idx_LW_0800_1200=5 ! [idx] Index of dust cnt. abs. 0800--1200 cm-1
  
  ! LSM grid size
  ! NB: These numbers must agree with parameters in dstlsm.F90
  integer,parameter::mpy_nbr=12 ! [nbr] Number of months per year
  integer,parameter::pln_nbr_LSM=3 ! [nbr] Number of plant types comprising a surface blend in LSM
  integer,parameter::sfc_typ_LSM_max=28 ! Maximum value of LSM surface type
  integer,parameter::sfc_typ_LSM_nbr=29 ! [nbr] Number of surface blends in LSM
  integer,parameter::vgi_nbr_LSM=14 ! [nbr] Number of vegetation types in LSM (mvt)
  
end module dstgrd

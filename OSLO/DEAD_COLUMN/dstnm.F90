! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstnm.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Common block dstnm stores variables names used in outfld()
! dst_nm_cmn is initialized in dst_nm_cmn_ini() which is called by dst_msc_cmn_ini()

! Usage: 
! use dstnm ! [mdl] Nomenclature for outfld()

module dstnm              ! [mdl] Nomenclature for outfld()
  use dstgrd,only:dst_nbr ! [mdl] Dust grid sizes
  implicit none
  save                      ! [stt] Changes to common variables are sticky
  
  ! Nomenclature info for outfld() initialized in dst_nm_cmn_ini()
  character(8) flx_mss_dry_sfc_nm(dst_nbr) ! [sng] Name of surface dry deposition flux
  character(8) flx_mss_grv_sfc_nm(dst_nbr) ! [sng] Name of surface gravitational deposition flux
  character(8) flx_mss_mbl_sfc_nm(dst_nbr) ! [sng] Name of surface mobilization flux
  character(8) flx_mss_pcp_sfc_nm(dst_nbr) ! [sng] Name of surface precipitation flux
  character(8) flx_mss_trb_sfc_nm(dst_nbr) ! [sng] Name of surface turbulent deposition flux
  character(8) mpc_dst_nm(dst_nbr) ! [sng] Name of column mass path of dust
  character(8) odxc_dst_nm(dst_nbr) ! [sng] Name of column extinction optical depth of dust
  character(8) flat_nm      ! [sng] Flux Longwave Absorbed Atmosphere
  character(8) flatc_nm     ! [sng] Flux Longwave Absorbed Atmosphere, clear sky
  character(8) flds_nm      ! [sng] Flux Longwave Downwelling Surface
  character(8) fsat_nm      ! [sng] Flux Shortwave Absorbed Atmosphere
  character(8) fsatc_nm     ! [sng] Flux Shortwave Absorbed Atmosphere, clear sky
  character(8) ftat_nm      ! [sng] Flux Total Absorbed Atmosphere
  character(8) ftns_nm      ! [sng] Flux Total Net Surface
  character(8) ftnt_nm      ! [sng] Flux Total Net TOA
  character(8) qrt_nm       ! [sng] Total radiative heating rate
  
  ! Dust forcing nomenclature is initialized in dst_nm_cmn_ini()
  character(8) flat_dst_nm  ! [sng]
  character(8) flat_frc_nm  ! [sng]
  character(8) flatc_dst_nm ! [sng]
  character(8) flatc_frc_nm ! [sng]
  character(8) flds_dst_nm  ! [sng]
  character(8) flds_frc_nm  ! [sng]
  character(8) flns_dst_nm  ! [sng]
  character(8) flns_frc_nm  ! [sng]
  character(8) flnsc_dst_nm ! [sng]
  character(8) flnsc_frc_nm ! [sng]
  character(8) flnt_dst_nm  ! [sng]
  character(8) flnt_frc_nm  ! [sng]
  character(8) flntc_dst_nm ! [sng]
  character(8) flntc_frc_nm ! [sng]
  character(8) fsat_dst_nm  ! [sng]
  character(8) fsat_frc_nm  ! [sng]
  character(8) fsatc_dst_nm ! [sng]
  character(8) fsatc_frc_nm ! [sng]
  character(8) fsds_dst_nm  ! [sng]
  character(8) fsds_frc_nm  ! [sng]
  character(8) fsns_dst_nm  ! [sng]
  character(8) fsns_frc_nm  ! [sng]
  character(8) fsnsc_dst_nm ! [sng]
  character(8) fsnsc_frc_nm ! [sng]
  character(8) fsnt_dst_nm  ! [sng]
  character(8) fsnt_frc_nm  ! [sng]
  character(8) fsntc_dst_nm ! [sng]
  character(8) fsntc_frc_nm ! [sng]
  character(8) ftat_dst_nm  ! [sng]
  character(8) ftat_frc_nm  ! [sng]
  character(8) ftns_dst_nm  ! [sng]
  character(8) ftns_frc_nm  ! [sng]
  character(8) ftnt_dst_nm  ! [sng]
  character(8) ftnt_frc_nm  ! [sng]
  character(8) qrl_dst_nm   ! [sng]
  character(8) qrl_frc_nm   ! [sng]
  character(8) qrs_dst_nm   ! [sng]
  character(8) qrs_frc_nm   ! [sng]
  character(8) qrt_dst_nm   ! [sng]
  character(8) qrt_frc_nm   ! [sng]
  
contains
  
  subroutine dst_nm_cmn_ini()
    ! Purpose: Initialize names of dust fields NOT stored in comtracnm
    ! dst_nm_cmn_ini() is called by CCM:physics/inti(), MATCH:inirun()
    ! NB: The call to dst_nm_cmn_ini() could be lumped into dst_msc_cmn_ini()
    ! except that MATCH requires dst_nm_cmn_ini() be called "earlier" than CCM
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstgrd ! [mdl] Dust grid sizes
    implicit none
    ! Local
    character(3) trc_nbr      ! Advected species number (Character)
    integer idx               ! cst index 
    
    ! Constituent-dependent arrays in dst_nm_cmn
    do idx=1,dst_nbr
       write (unit=trc_nbr,fmt='(i3)') idx+100
       flx_mss_mbl_sfc_nm(idx)='DSTSFM'//trc_nbr(2:3)
       flx_mss_pcp_sfc_nm(idx)='DSTSFP'//trc_nbr(2:3)
       flx_mss_trb_sfc_nm(idx)='DSTSFT'//trc_nbr(2:3)
       flx_mss_grv_sfc_nm(idx)='DSTSFG'//trc_nbr(2:3)
       flx_mss_dry_sfc_nm(idx)='DSTSFD'//trc_nbr(2:3)
       mpc_dst_nm(idx)='DSTMPC'//trc_nbr(2:3)
       odxc_dst_nm(idx)='DSTODX'//trc_nbr(2:3)
    end do                    ! end loop over cst
    
    ! Nomenclature for standard radiative fields (currently also in dst_nm_cmn)
    ! NB: These names are used by CCM, but not MATCH or MOZART
    flat_nm='FLAT    '
    flatc_nm='FLATC   '
    flds_nm='FLDS    '
    fsat_nm='FSAT    '
    fsatc_nm='FSATC   '
    ftat_nm='FTAT    '
    ftns_nm='FTNS    '
    ftnt_nm='FTNT    '
    qrt_nm='QRT     '
    
    ! Nomenclature for radiative forcing fields in dst_frc_nm_cmn
    ! NB: These names are used by CCM, but not MATCH or MOZART
    flat_dst_nm='FLATDST '
    flat_frc_nm='FLATFRC '
    flatc_dst_nm='FLATCDST'
    flatc_frc_nm='FLATCFRC'
    flds_dst_nm='FLDSDST '
    flds_frc_nm='FLDSFRC '
    flns_dst_nm='FLNSDST '
    flns_frc_nm='FLNSFRC '
    flnsc_dst_nm='FLNSCDST'
    flnsc_frc_nm='FLNSCFRC'
    flnt_dst_nm='FLNTDST '
    flnt_frc_nm='FLNTFRC '
    flntc_dst_nm='FLNTCDST'
    flntc_frc_nm='FLNTCFRC'
    fsat_dst_nm='FSATDST '
    fsat_frc_nm='FSATFRC '
    fsatc_dst_nm='FSATCDST'
    fsatc_frc_nm='FSATCFRC'
    fsds_dst_nm='FSDSDST '
    fsds_frc_nm='FSDSFRC '
    fsns_dst_nm='FSNSDST '
    fsns_frc_nm='FSNSFRC '
    fsnsc_dst_nm='FSNSCDST'
    fsnsc_frc_nm='FSNSCFRC'
    fsnt_dst_nm='FSNTDST '
    fsnt_frc_nm='FSNTFRC '
    fsntc_dst_nm='FSNTCDST'
    fsntc_frc_nm='FSNTCFRC'
    ftat_dst_nm='FTATDST '
    ftat_frc_nm='FTATFRC '
    ftns_dst_nm='FTNSDST '
    ftns_frc_nm='FTNSFRC '
    ftnt_dst_nm='FTNTDST '
    ftnt_frc_nm='FTNTFRC '
    qrl_dst_nm='QRLDST  '
    qrl_frc_nm='QRLFRC  '
    qrs_dst_nm='QRSDST  '
    qrs_frc_nm='QRSFRC  '
    qrt_dst_nm='QRTDST  '
    qrt_frc_nm='QRTFRC  '
    
    return
  end subroutine dst_nm_cmn_ini                       ! end dst_nm_cmn_ini()
  
end module dstnm          ! [mdl] Nomenclature for outfld()

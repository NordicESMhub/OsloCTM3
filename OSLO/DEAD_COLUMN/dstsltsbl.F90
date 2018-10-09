! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstsltsbl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Information needed to compute saltation sandblasting
! These variables are initialized in dst_slt_sbl_cmn_ini() which is called by dst_msc_cmn_ini() 

! Usage:
! use dstsltsbl ! [mdl] Saltation sandblasting physics

module dstsltsbl ! [mdl] Dust saltation sandblasting physics
  use dstgrd,only:dst_nbr,dst_src_nbr,bln_nbr,wnd_frc_nbr ! [mdl] Dust grid sizes,soil types available, wnd_frict
  ! ions speeds for which soil properties are given 
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  save ! [stt] Changes to common variables are sticky
  
  ! Saltation sandblasting info initialized in dst_slt_sbl_cmn_ini
  real(r8) mss_frc_src_lut(wnd_frc_nbr,bln_nbr,dst_src_nbr) ! [frc] lut for mass fraction of source
  real(r8) dst_slt_flx_rat_ttl_lut(wnd_frc_nbr,bln_nbr) ! [m-1] lut for ratio of vertical to horizontal dust flux
  real(r8) flx_mss_hrz_slt_ttl_lut(wnd_frc_nbr,bln_nbr) ! [kg m-1] lut for horizontal saltation flux (complex)
  
contains
  
  subroutine dst_slt_sbl_cmn_ini()
    !Purpose: Initialize look up tables for saltation and sandblasting
    !dst_slt_sbl_cmn_ini() is called from dst_msc_cmn_ini()
    !Each wind speed will give a distribution between the source modes
    !Each wind speed will give a alpha (vertical to horizontal flux)
    !Theory is given in Alfaro/Gomes, 2001(JGR), and Grini et al 2002 (GRL).
    !Dependencies: None, exept constants
    !Source code for inputs accessible by cvs
    !cvs -d dust.ess.uci.edu:/home/zender/cvs co -kk sltsbl
    !Yes, the numbers are generated from a program. I did not actually write them in by hand...

    !Alf Grini: alf.grini@geofysikk.uio.no

    use dstgrd, only:dst_src_nbr,bln_nbr,wnd_frc_nbr !(see above)
    implicit none
    !LOCAL VARIABLE:
    real(r8) :: mss_frc_src_lutx(wnd_frc_nbr,bln_nbr,dst_src_nbr)  !look up table for mss_frc_src
    real(r8) :: dst_slt_flx_rat_ttl_lutx(wnd_frc_nbr,bln_nbr)          !look up table for "alpha"
    real(r8) :: flx_mss_hrz_slt_ttl_lutx(wnd_frc_nbr,bln_nbr)       !Loop up table for horizontal mass flux
    integer src_idx                                            !Counter for source modes


        ! soil type    1  mmd=   125.0  um  sigma=     1.8
            ! first index is wind friction speed in cm/s
                            ! second index is soil type 
                  ! third index is number of source mode
      data (mss_frc_src_lutx(  1,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  2,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  3,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  4,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  5,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  6,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  7,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  8,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  9,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 10,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 11,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 12,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 13,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 14,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 15,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 16,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 17,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 18,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 19,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 20,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 21,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 22,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 23,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 24,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 25,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 26,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 27,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 28,1,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 29,1,src_idx),src_idx=1,dst_src_nbr) /   0.961074E-03,   0.798130E-01,   0.919226E+00/
      data (mss_frc_src_lutx( 30,1,src_idx),src_idx=1,dst_src_nbr) /   0.241670E-02,   0.104497E+00,   0.893086E+00/
      data (mss_frc_src_lutx( 31,1,src_idx),src_idx=1,dst_src_nbr) /   0.361964E-02,   0.115535E+00,   0.880845E+00/
      data (mss_frc_src_lutx( 32,1,src_idx),src_idx=1,dst_src_nbr) /   0.462186E-02,   0.122360E+00,   0.873018E+00/
      data (mss_frc_src_lutx( 33,1,src_idx),src_idx=1,dst_src_nbr) /   0.537328E-02,   0.125058E+00,   0.869569E+00/
      data (mss_frc_src_lutx( 34,1,src_idx),src_idx=1,dst_src_nbr) /   0.603354E-02,   0.127700E+00,   0.866263E+00/
      data (mss_frc_src_lutx( 35,1,src_idx),src_idx=1,dst_src_nbr) /   0.666543E-02,   0.131124E+00,   0.862207E+00/
      data (mss_frc_src_lutx( 36,1,src_idx),src_idx=1,dst_src_nbr) /   0.716101E-02,   0.132948E+00,   0.859888E+00/
      data (mss_frc_src_lutx( 37,1,src_idx),src_idx=1,dst_src_nbr) /   0.753756E-02,   0.133526E+00,   0.858934E+00/
      data (mss_frc_src_lutx( 38,1,src_idx),src_idx=1,dst_src_nbr) /   0.799124E-02,   0.136091E+00,   0.855914E+00/
      data (mss_frc_src_lutx( 39,1,src_idx),src_idx=1,dst_src_nbr) /   0.838104E-02,   0.137991E+00,   0.853624E+00/
      data (mss_frc_src_lutx( 40,1,src_idx),src_idx=1,dst_src_nbr) /   0.872147E-02,   0.139422E+00,   0.851853E+00/
      data (mss_frc_src_lutx( 41,1,src_idx),src_idx=1,dst_src_nbr) /   0.882945E-02,   0.137490E+00,   0.853680E+00/
      data (mss_frc_src_lutx( 42,1,src_idx),src_idx=1,dst_src_nbr) /   0.910555E-02,   0.138474E+00,   0.852421E+00/
      data (mss_frc_src_lutx( 43,1,src_idx),src_idx=1,dst_src_nbr) /   0.937286E-02,   0.139495E+00,   0.851128E+00/
      data (mss_frc_src_lutx( 44,1,src_idx),src_idx=1,dst_src_nbr) /   0.964377E-02,   0.140702E+00,   0.849649E+00/
      data (mss_frc_src_lutx( 45,1,src_idx),src_idx=1,dst_src_nbr) /   0.993009E-02,   0.142227E+00,   0.847839E+00/
      data (mss_frc_src_lutx( 46,1,src_idx),src_idx=1,dst_src_nbr) /   0.100268E-01,   0.141180E+00,   0.848789E+00/
      data (mss_frc_src_lutx( 47,1,src_idx),src_idx=1,dst_src_nbr) /   0.103753E-01,   0.143734E+00,   0.845887E+00/
      data (mss_frc_src_lutx( 48,1,src_idx),src_idx=1,dst_src_nbr) /   0.105460E-01,   0.143875E+00,   0.845576E+00/
      data (mss_frc_src_lutx( 49,1,src_idx),src_idx=1,dst_src_nbr) /   0.107717E-01,   0.144843E+00,   0.844382E+00/
      data (mss_frc_src_lutx( 50,1,src_idx),src_idx=1,dst_src_nbr) /   0.108256E-01,   0.143576E+00,   0.845597E+00/
      data (mss_frc_src_lutx( 51,1,src_idx),src_idx=1,dst_src_nbr) /   0.111850E-01,   0.146386E+00,   0.842427E+00/
      data (mss_frc_src_lutx( 52,1,src_idx),src_idx=1,dst_src_nbr) /   0.111443E-01,   0.144010E+00,   0.844845E+00/
      data (mss_frc_src_lutx( 53,1,src_idx),src_idx=1,dst_src_nbr) /   0.114278E-01,   0.145862E+00,   0.842710E+00/
      data (mss_frc_src_lutx( 54,1,src_idx),src_idx=1,dst_src_nbr) /   0.115629E-01,   0.145853E+00,   0.842583E+00/
      data (mss_frc_src_lutx( 55,1,src_idx),src_idx=1,dst_src_nbr) /   0.117967E-01,   0.147089E+00,   0.841114E+00/
      data (mss_frc_src_lutx( 56,1,src_idx),src_idx=1,dst_src_nbr) /   0.118886E-01,   0.146606E+00,   0.841502E+00/
      data (mss_frc_src_lutx( 57,1,src_idx),src_idx=1,dst_src_nbr) /   0.120920E-01,   0.147500E+00,   0.840406E+00/
      data (mss_frc_src_lutx( 58,1,src_idx),src_idx=1,dst_src_nbr) /   0.121594E-01,   0.146776E+00,   0.841064E+00/
      data (mss_frc_src_lutx( 59,1,src_idx),src_idx=1,dst_src_nbr) /   0.123514E-01,   0.147563E+00,   0.840086E+00/
      data (mss_frc_src_lutx( 60,1,src_idx),src_idx=1,dst_src_nbr) /   0.126800E-01,   0.149973E+00,   0.837347E+00/
      data (mss_frc_src_lutx( 61,1,src_idx),src_idx=1,dst_src_nbr) /   0.126119E-01,   0.147710E+00,   0.839678E+00/
      data (mss_frc_src_lutx( 62,1,src_idx),src_idx=1,dst_src_nbr) /   0.129578E-01,   0.150298E+00,   0.836750E+00/
      data (mss_frc_src_lutx( 63,1,src_idx),src_idx=1,dst_src_nbr) /   0.129102E-01,   0.148346E+00,   0.838751E+00/
      data (mss_frc_src_lutx( 64,1,src_idx),src_idx=1,dst_src_nbr) /   0.132940E-01,   0.151342E+00,   0.835369E+00/
      data (mss_frc_src_lutx( 65,1,src_idx),src_idx=1,dst_src_nbr) /   0.132820E-01,   0.149836E+00,   0.836888E+00/
      data (mss_frc_src_lutx( 66,1,src_idx),src_idx=1,dst_src_nbr) /   0.134362E-01,   0.150225E+00,   0.836346E+00/
      data (mss_frc_src_lutx( 67,1,src_idx),src_idx=1,dst_src_nbr) /   0.137664E-01,   0.152569E+00,   0.833672E+00/
      data (mss_frc_src_lutx( 68,1,src_idx),src_idx=1,dst_src_nbr) /   0.136962E-01,   0.150483E+00,   0.835826E+00/
      data (mss_frc_src_lutx( 69,1,src_idx),src_idx=1,dst_src_nbr) /   0.138078E-01,   0.150421E+00,   0.835778E+00/
      data (mss_frc_src_lutx( 70,1,src_idx),src_idx=1,dst_src_nbr) /   0.141049E-01,   0.152354E+00,   0.833548E+00/
      data (mss_frc_src_lutx( 71,1,src_idx),src_idx=1,dst_src_nbr) /   0.143047E-01,   0.153241E+00,   0.832461E+00/
      data (mss_frc_src_lutx( 72,1,src_idx),src_idx=1,dst_src_nbr) /   0.143997E-01,   0.152992E+00,   0.832614E+00/
      data (mss_frc_src_lutx( 73,1,src_idx),src_idx=1,dst_src_nbr) /   0.143974E-01,   0.151755E+00,   0.833852E+00/
      data (mss_frc_src_lutx( 74,1,src_idx),src_idx=1,dst_src_nbr) /   0.145978E-01,   0.152638E+00,   0.832769E+00/
      data (mss_frc_src_lutx( 75,1,src_idx),src_idx=1,dst_src_nbr) /   0.147035E-01,   0.152539E+00,   0.832757E+00/
      data (mss_frc_src_lutx( 76,1,src_idx),src_idx=1,dst_src_nbr) /   0.150246E-01,   0.154649E+00,   0.830327E+00/
      data (mss_frc_src_lutx( 77,1,src_idx),src_idx=1,dst_src_nbr) /   0.149431E-01,   0.152638E+00,   0.832424E+00/
      data (mss_frc_src_lutx( 78,1,src_idx),src_idx=1,dst_src_nbr) /   0.150822E-01,   0.152888E+00,   0.832032E+00/
      data (mss_frc_src_lutx( 79,1,src_idx),src_idx=1,dst_src_nbr) /   0.154496E-01,   0.155424E+00,   0.829127E+00/
      data (mss_frc_src_lutx( 80,1,src_idx),src_idx=1,dst_src_nbr) /   0.154125E-01,   0.153908E+00,   0.830682E+00/
      data (mss_frc_src_lutx( 81,1,src_idx),src_idx=1,dst_src_nbr) /   0.156094E-01,   0.154721E+00,   0.829671E+00/
      data (mss_frc_src_lutx( 82,1,src_idx),src_idx=1,dst_src_nbr) /   0.157228E-01,   0.154711E+00,   0.829569E+00/
      data (mss_frc_src_lutx( 83,1,src_idx),src_idx=1,dst_src_nbr) /   0.160832E-01,   0.157110E+00,   0.826810E+00/
      data (mss_frc_src_lutx( 84,1,src_idx),src_idx=1,dst_src_nbr) /   0.160288E-01,   0.155458E+00,   0.828515E+00/
      data (mss_frc_src_lutx( 85,1,src_idx),src_idx=1,dst_src_nbr) /   0.162264E-01,   0.156250E+00,   0.827526E+00/
      data (mss_frc_src_lutx( 86,1,src_idx),src_idx=1,dst_src_nbr) /   0.163462E-01,   0.156299E+00,   0.827358E+00/
      data (mss_frc_src_lutx( 87,1,src_idx),src_idx=1,dst_src_nbr) /   0.163878E-01,   0.155611E+00,   0.828005E+00/
      data (mss_frc_src_lutx( 88,1,src_idx),src_idx=1,dst_src_nbr) /   0.166909E-01,   0.157375E+00,   0.825935E+00/
      data (mss_frc_src_lutx( 89,1,src_idx),src_idx=1,dst_src_nbr) /   0.169249E-01,   0.158495E+00,   0.824581E+00/
      data (mss_frc_src_lutx( 90,1,src_idx),src_idx=1,dst_src_nbr) /   0.170833E-01,   0.158890E+00,   0.824030E+00/
      data (mss_frc_src_lutx( 91,1,src_idx),src_idx=1,dst_src_nbr) /   0.171654E-01,   0.158566E+00,   0.824270E+00/
      data (mss_frc_src_lutx( 92,1,src_idx),src_idx=1,dst_src_nbr) /   0.171760E-01,   0.157608E+00,   0.825217E+00/
      data (mss_frc_src_lutx( 93,1,src_idx),src_idx=1,dst_src_nbr) /   0.174696E-01,   0.159229E+00,   0.823300E+00/
      data (mss_frc_src_lutx( 94,1,src_idx),src_idx=1,dst_src_nbr) /   0.173364E-01,   0.156978E+00,   0.825687E+00/
      data (mss_frc_src_lutx( 95,1,src_idx),src_idx=1,dst_src_nbr) /   0.174899E-01,   0.157329E+00,   0.825183E+00/
      data (mss_frc_src_lutx( 96,1,src_idx),src_idx=1,dst_src_nbr) /   0.179400E-01,   0.160322E+00,   0.821739E+00/
      data (mss_frc_src_lutx( 97,1,src_idx),src_idx=1,dst_src_nbr) /   0.179566E-01,   0.159425E+00,   0.822621E+00/
      data (mss_frc_src_lutx( 98,1,src_idx),src_idx=1,dst_src_nbr) /   0.179094E-01,   0.157994E+00,   0.824097E+00/
      data (mss_frc_src_lutx( 99,1,src_idx),src_idx=1,dst_src_nbr) /   0.181624E-01,   0.159200E+00,   0.822639E+00/
      data (mss_frc_src_lutx(100,1,src_idx),src_idx=1,dst_src_nbr) /   0.183526E-01,   0.159841E+00,   0.821808E+00/
        ! soil type    2  mmd=   210.0  um  sigma=     1.6
            ! first index is wind friction speed in cm/s
                            ! second index is soil type 
                  ! third index is number of source mode
      data (mss_frc_src_lutx(  1,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  2,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  3,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  4,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  5,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  6,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  7,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  8,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  9,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 10,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 11,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 12,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 13,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 14,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 15,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 16,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 17,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 18,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 19,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 20,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 21,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 22,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 23,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 24,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 25,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 26,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 27,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 28,2,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 29,2,src_idx),src_idx=1,dst_src_nbr) /   0.100510E-02,   0.815546E-01,   0.917439E+00/
      data (mss_frc_src_lutx( 30,2,src_idx),src_idx=1,dst_src_nbr) /   0.263153E-02,   0.108297E+00,   0.889072E+00/
      data (mss_frc_src_lutx( 31,2,src_idx),src_idx=1,dst_src_nbr) /   0.408181E-02,   0.120952E+00,   0.874968E+00/
      data (mss_frc_src_lutx( 32,2,src_idx),src_idx=1,dst_src_nbr) /   0.537310E-02,   0.129073E+00,   0.865554E+00/
      data (mss_frc_src_lutx( 33,2,src_idx),src_idx=1,dst_src_nbr) /   0.641896E-02,   0.132776E+00,   0.860803E+00/
      data (mss_frc_src_lutx( 34,2,src_idx),src_idx=1,dst_src_nbr) /   0.738243E-02,   0.136290E+00,   0.856325E+00/
      data (mss_frc_src_lutx( 35,2,src_idx),src_idx=1,dst_src_nbr) /   0.832974E-02,   0.140525E+00,   0.851143E+00/
      data (mss_frc_src_lutx( 36,2,src_idx),src_idx=1,dst_src_nbr) /   0.912325E-02,   0.143039E+00,   0.847835E+00/
      data (mss_frc_src_lutx( 37,2,src_idx),src_idx=1,dst_src_nbr) /   0.977472E-02,   0.144209E+00,   0.846013E+00/
      data (mss_frc_src_lutx( 38,2,src_idx),src_idx=1,dst_src_nbr) /   0.105285E-01,   0.147416E+00,   0.842051E+00/
      data (mss_frc_src_lutx( 39,2,src_idx),src_idx=1,dst_src_nbr) /   0.112092E-01,   0.149906E+00,   0.838880E+00/
      data (mss_frc_src_lutx( 40,2,src_idx),src_idx=1,dst_src_nbr) /   0.118333E-01,   0.151889E+00,   0.836273E+00/
      data (mss_frc_src_lutx( 41,2,src_idx),src_idx=1,dst_src_nbr) /   0.121581E-01,   0.150339E+00,   0.837501E+00/
      data (mss_frc_src_lutx( 42,2,src_idx),src_idx=1,dst_src_nbr) /   0.127073E-01,   0.151825E+00,   0.835463E+00/
      data (mss_frc_src_lutx( 43,2,src_idx),src_idx=1,dst_src_nbr) /   0.132507E-01,   0.153342E+00,   0.833403E+00/
      data (mss_frc_src_lutx( 44,2,src_idx),src_idx=1,dst_src_nbr) /   0.138055E-01,   0.155047E+00,   0.831144E+00/
      data (mss_frc_src_lutx( 45,2,src_idx),src_idx=1,dst_src_nbr) /   0.143882E-01,   0.157079E+00,   0.828527E+00/
      data (mss_frc_src_lutx( 46,2,src_idx),src_idx=1,dst_src_nbr) /   0.147143E-01,   0.156400E+00,   0.828883E+00/
      data (mss_frc_src_lutx( 47,2,src_idx),src_idx=1,dst_src_nbr) /   0.153973E-01,   0.159508E+00,   0.825090E+00/
      data (mss_frc_src_lutx( 48,2,src_idx),src_idx=1,dst_src_nbr) /   0.158372E-01,   0.160071E+00,   0.824087E+00/
      data (mss_frc_src_lutx( 49,2,src_idx),src_idx=1,dst_src_nbr) /   0.163611E-01,   0.161501E+00,   0.822134E+00/
      data (mss_frc_src_lutx( 50,2,src_idx),src_idx=1,dst_src_nbr) /   0.166422E-01,   0.160565E+00,   0.822789E+00/
      data (mss_frc_src_lutx( 51,2,src_idx),src_idx=1,dst_src_nbr) /   0.173737E-01,   0.163938E+00,   0.818684E+00/
      data (mss_frc_src_lutx( 52,2,src_idx),src_idx=1,dst_src_nbr) /   0.175231E-01,   0.161819E+00,   0.820655E+00/
      data (mss_frc_src_lutx( 53,2,src_idx),src_idx=1,dst_src_nbr) /   0.181575E-01,   0.164177E+00,   0.817662E+00/
      data (mss_frc_src_lutx( 54,2,src_idx),src_idx=1,dst_src_nbr) /   0.185765E-01,   0.164557E+00,   0.816862E+00/
      data (mss_frc_src_lutx( 55,2,src_idx),src_idx=1,dst_src_nbr) /   0.191510E-01,   0.166260E+00,   0.814587E+00/
      data (mss_frc_src_lutx( 56,2,src_idx),src_idx=1,dst_src_nbr) /   0.195151E-01,   0.166132E+00,   0.814352E+00/
      data (mss_frc_src_lutx( 57,2,src_idx),src_idx=1,dst_src_nbr) /   0.200562E-01,   0.167466E+00,   0.812476E+00/
      data (mss_frc_src_lutx( 58,2,src_idx),src_idx=1,dst_src_nbr) /   0.203913E-01,   0.167076E+00,   0.812531E+00/
      data (mss_frc_src_lutx( 59,2,src_idx),src_idx=1,dst_src_nbr) /   0.209269E-01,   0.168294E+00,   0.810778E+00/
      data (mss_frc_src_lutx( 60,2,src_idx),src_idx=1,dst_src_nbr) /   0.216879E-01,   0.171239E+00,   0.807071E+00/
      data (mss_frc_src_lutx( 61,2,src_idx),src_idx=1,dst_src_nbr) /   0.218207E-01,   0.169202E+00,   0.808977E+00/
      data (mss_frc_src_lutx( 62,2,src_idx),src_idx=1,dst_src_nbr) /   0.226264E-01,   0.172338E+00,   0.805033E+00/
      data (mss_frc_src_lutx( 63,2,src_idx),src_idx=1,dst_src_nbr) /   0.227983E-01,   0.170624E+00,   0.806576E+00/
      data (mss_frc_src_lutx( 64,2,src_idx),src_idx=1,dst_src_nbr) /   0.236848E-01,   0.174198E+00,   0.802115E+00/
      data (mss_frc_src_lutx( 65,2,src_idx),src_idx=1,dst_src_nbr) /   0.239235E-01,   0.172958E+00,   0.803115E+00/
      data (mss_frc_src_lutx( 66,2,src_idx),src_idx=1,dst_src_nbr) /   0.244439E-01,   0.173741E+00,   0.801813E+00/
      data (mss_frc_src_lutx( 67,2,src_idx),src_idx=1,dst_src_nbr) /   0.252695E-01,   0.176617E+00,   0.798109E+00/
      data (mss_frc_src_lutx( 68,2,src_idx),src_idx=1,dst_src_nbr) /   0.254213E-01,   0.174748E+00,   0.799829E+00/
      data (mss_frc_src_lutx( 69,2,src_idx),src_idx=1,dst_src_nbr) /   0.258870E-01,   0.175043E+00,   0.799068E+00/
      data (mss_frc_src_lutx( 70,2,src_idx),src_idx=1,dst_src_nbr) /   0.266818E-01,   0.177481E+00,   0.795835E+00/
      data (mss_frc_src_lutx( 71,2,src_idx),src_idx=1,dst_src_nbr) /   0.273163E-01,   0.178793E+00,   0.793886E+00/
      data (mss_frc_src_lutx( 72,2,src_idx),src_idx=1,dst_src_nbr) /   0.277736E-01,   0.178888E+00,   0.793337E+00/
      data (mss_frc_src_lutx( 73,2,src_idx),src_idx=1,dst_src_nbr) /   0.280623E-01,   0.177915E+00,   0.794021E+00/
      data (mss_frc_src_lutx( 74,2,src_idx),src_idx=1,dst_src_nbr) /   0.287182E-01,   0.179221E+00,   0.792057E+00/
      data (mss_frc_src_lutx( 75,2,src_idx),src_idx=1,dst_src_nbr) /   0.292108E-01,   0.179469E+00,   0.791317E+00/
      data (mss_frc_src_lutx( 76,2,src_idx),src_idx=1,dst_src_nbr) /   0.301033E-01,   0.182096E+00,   0.787796E+00/
      data (mss_frc_src_lutx( 77,2,src_idx),src_idx=1,dst_src_nbr) /   0.302648E-01,   0.180282E+00,   0.789451E+00/
      data (mss_frc_src_lutx( 78,2,src_idx),src_idx=1,dst_src_nbr) /   0.308372E-01,   0.180903E+00,   0.788257E+00/
      data (mss_frc_src_lutx( 79,2,src_idx),src_idx=1,dst_src_nbr) /   0.318447E-01,   0.183988E+00,   0.784164E+00/
      data (mss_frc_src_lutx( 80,2,src_idx),src_idx=1,dst_src_nbr) /   0.321010E-01,   0.182703E+00,   0.785193E+00/
      data (mss_frc_src_lutx( 81,2,src_idx),src_idx=1,dst_src_nbr) /   0.328050E-01,   0.183929E+00,   0.783262E+00/
      data (mss_frc_src_lutx( 82,2,src_idx),src_idx=1,dst_src_nbr) /   0.333573E-01,   0.184264E+00,   0.782374E+00/
      data (mss_frc_src_lutx( 83,2,src_idx),src_idx=1,dst_src_nbr) /   0.343942E-01,   0.187198E+00,   0.778406E+00/
      data (mss_frc_src_lutx( 84,2,src_idx),src_idx=1,dst_src_nbr) /   0.346354E-01,   0.185761E+00,   0.779600E+00/
      data (mss_frc_src_lutx( 85,2,src_idx),src_idx=1,dst_src_nbr) /   0.353736E-01,   0.186960E+00,   0.777663E+00/
      data (mss_frc_src_lutx( 86,2,src_idx),src_idx=1,dst_src_nbr) /   0.359668E-01,   0.187356E+00,   0.776676E+00/
      data (mss_frc_src_lutx( 87,2,src_idx),src_idx=1,dst_src_nbr) /   0.364101E-01,   0.186952E+00,   0.776632E+00/
      data (mss_frc_src_lutx( 88,2,src_idx),src_idx=1,dst_src_nbr) /   0.373842E-01,   0.189201E+00,   0.773409E+00/
      data (mss_frc_src_lutx( 89,2,src_idx),src_idx=1,dst_src_nbr) /   0.382306E-01,   0.190749E+00,   0.771016E+00/
      data (mss_frc_src_lutx( 90,2,src_idx),src_idx=1,dst_src_nbr) /   0.389328E-01,   0.191514E+00,   0.769551E+00/
      data (mss_frc_src_lutx( 91,2,src_idx),src_idx=1,dst_src_nbr) /   0.394860E-01,   0.191502E+00,   0.769008E+00/
      data (mss_frc_src_lutx( 92,2,src_idx),src_idx=1,dst_src_nbr) /   0.398971E-01,   0.190801E+00,   0.769297E+00/
      data (mss_frc_src_lutx( 93,2,src_idx),src_idx=1,dst_src_nbr) /   0.409027E-01,   0.192888E+00,   0.766204E+00/
      data (mss_frc_src_lutx( 94,2,src_idx),src_idx=1,dst_src_nbr) /   0.410224E-01,   0.190786E+00,   0.768188E+00/
      data (mss_frc_src_lutx( 95,2,src_idx),src_idx=1,dst_src_nbr) /   0.417494E-01,   0.191497E+00,   0.766748E+00/
      data (mss_frc_src_lutx( 96,2,src_idx),src_idx=1,dst_src_nbr) /   0.431178E-01,   0.195063E+00,   0.761812E+00/
      data (mss_frc_src_lutx( 97,2,src_idx),src_idx=1,dst_src_nbr) /   0.435716E-01,   0.194426E+00,   0.762000E+00/
      data (mss_frc_src_lutx( 98,2,src_idx),src_idx=1,dst_src_nbr) /   0.438907E-01,   0.193204E+00,   0.762902E+00/
      data (mss_frc_src_lutx( 99,2,src_idx),src_idx=1,dst_src_nbr) /   0.448672E-01,   0.194835E+00,   0.760291E+00/
      data (mss_frc_src_lutx(100,2,src_idx),src_idx=1,dst_src_nbr) /   0.457169E-01,   0.195855E+00,   0.758423E+00/
        ! soil type    3  mmd=   520.0  um  sigma=     1.6
            ! first index is wind friction speed in cm/s
                            ! second index is soil type 
                  ! third index is number of source mode
      data (mss_frc_src_lutx(  1,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  2,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  3,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  4,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  5,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  6,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  7,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  8,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  9,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 10,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 11,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 12,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 13,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 14,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 15,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 16,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 17,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 18,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 19,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 20,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 21,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 22,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 23,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 24,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 25,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 26,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 27,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 28,3,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 29,3,src_idx),src_idx=1,dst_src_nbr) /   0.117581E-02,   0.879637E-01,   0.910860E+00/
      data (mss_frc_src_lutx( 30,3,src_idx),src_idx=1,dst_src_nbr) /   0.359946E-02,   0.123497E+00,   0.872903E+00/
      data (mss_frc_src_lutx( 31,3,src_idx),src_idx=1,dst_src_nbr) /   0.648268E-02,   0.144134E+00,   0.849384E+00/
      data (mss_frc_src_lutx( 32,3,src_idx),src_idx=1,dst_src_nbr) /   0.983114E-02,   0.159374E+00,   0.830795E+00/
      data (mss_frc_src_lutx( 33,3,src_idx),src_idx=1,dst_src_nbr) /   0.134465E-01,   0.169089E+00,   0.817463E+00/
      data (mss_frc_src_lutx( 34,3,src_idx),src_idx=1,dst_src_nbr) /   0.175471E-01,   0.177972E+00,   0.804481E+00/
      data (mss_frc_src_lutx( 35,3,src_idx),src_idx=1,dst_src_nbr) /   0.222499E-01,   0.187134E+00,   0.790620E+00/
      data (mss_frc_src_lutx( 36,3,src_idx),src_idx=1,dst_src_nbr) /   0.271135E-01,   0.193783E+00,   0.779104E+00/
      data (mss_frc_src_lutx( 37,3,src_idx),src_idx=1,dst_src_nbr) /   0.319297E-01,   0.198385E+00,   0.769686E+00/
      data (mss_frc_src_lutx( 38,3,src_idx),src_idx=1,dst_src_nbr) /   0.373649E-01,   0.204979E+00,   0.757657E+00/
      data (mss_frc_src_lutx( 39,3,src_idx),src_idx=1,dst_src_nbr) /   0.429426E-01,   0.210411E+00,   0.746646E+00/
      data (mss_frc_src_lutx( 40,3,src_idx),src_idx=1,dst_src_nbr) /   0.486826E-01,   0.214970E+00,   0.736346E+00/
      data (mss_frc_src_lutx( 41,3,src_idx),src_idx=1,dst_src_nbr) /   0.536869E-01,   0.215193E+00,   0.731117E+00/
      data (mss_frc_src_lutx( 42,3,src_idx),src_idx=1,dst_src_nbr) /   0.597522E-01,   0.218724E+00,   0.721522E+00/
      data (mss_frc_src_lutx( 43,3,src_idx),src_idx=1,dst_src_nbr) /   0.660944E-01,   0.222091E+00,   0.711812E+00/
      data (mss_frc_src_lutx( 44,3,src_idx),src_idx=1,dst_src_nbr) /   0.727747E-01,   0.225484E+00,   0.701738E+00/
      data (mss_frc_src_lutx( 45,3,src_idx),src_idx=1,dst_src_nbr) /   0.798623E-01,   0.229068E+00,   0.691066E+00/
      data (mss_frc_src_lutx( 46,3,src_idx),src_idx=1,dst_src_nbr) /   0.860921E-01,   0.229429E+00,   0.684476E+00/
      data (mss_frc_src_lutx( 47,3,src_idx),src_idx=1,dst_src_nbr) /   0.941669E-01,   0.233880E+00,   0.671948E+00/
      data (mss_frc_src_lutx( 48,3,src_idx),src_idx=1,dst_src_nbr) /   0.101379E+00,   0.235348E+00,   0.663269E+00/
      data (mss_frc_src_lutx( 49,3,src_idx),src_idx=1,dst_src_nbr) /   0.109217E+00,   0.237620E+00,   0.653160E+00/
      data (mss_frc_src_lutx( 50,3,src_idx),src_idx=1,dst_src_nbr) /   0.116055E+00,   0.237215E+00,   0.646729E+00/
      data (mss_frc_src_lutx( 51,3,src_idx),src_idx=1,dst_src_nbr) /   0.125375E+00,   0.241282E+00,   0.633337E+00/
      data (mss_frc_src_lutx( 52,3,src_idx),src_idx=1,dst_src_nbr) /   0.131870E+00,   0.239436E+00,   0.628693E+00/
      data (mss_frc_src_lutx( 53,3,src_idx),src_idx=1,dst_src_nbr) /   0.141094E+00,   0.242156E+00,   0.616745E+00/
      data (mss_frc_src_lutx( 54,3,src_idx),src_idx=1,dst_src_nbr) /   0.149337E+00,   0.242707E+00,   0.607952E+00/
      data (mss_frc_src_lutx( 55,3,src_idx),src_idx=1,dst_src_nbr) /   0.158597E+00,   0.244476E+00,   0.596918E+00/
      data (mss_frc_src_lutx( 56,3,src_idx),src_idx=1,dst_src_nbr) /   0.166829E+00,   0.244307E+00,   0.588862E+00/
      data (mss_frc_src_lutx( 57,3,src_idx),src_idx=1,dst_src_nbr) /   0.176166E+00,   0.245479E+00,   0.578353E+00/
      data (mss_frc_src_lutx( 58,3,src_idx),src_idx=1,dst_src_nbr) /   0.184393E+00,   0.244903E+00,   0.570703E+00/
      data (mss_frc_src_lutx( 59,3,src_idx),src_idx=1,dst_src_nbr) /   0.193866E+00,   0.245769E+00,   0.560361E+00/
      data (mss_frc_src_lutx( 60,3,src_idx),src_idx=1,dst_src_nbr) /   0.204752E+00,   0.248089E+00,   0.547155E+00/
      data (mss_frc_src_lutx( 61,3,src_idx),src_idx=1,dst_src_nbr) /   0.211991E+00,   0.245804E+00,   0.542199E+00/
      data (mss_frc_src_lutx( 62,3,src_idx),src_idx=1,dst_src_nbr) /   0.223300E+00,   0.248043E+00,   0.528653E+00/
      data (mss_frc_src_lutx( 63,3,src_idx),src_idx=1,dst_src_nbr) /   0.230878E+00,   0.245956E+00,   0.523162E+00/
      data (mss_frc_src_lutx( 64,3,src_idx),src_idx=1,dst_src_nbr) /   0.242795E+00,   0.248287E+00,   0.508914E+00/
      data (mss_frc_src_lutx( 65,3,src_idx),src_idx=1,dst_src_nbr) /   0.250874E+00,   0.246495E+00,   0.502629E+00/
      data (mss_frc_src_lutx( 66,3,src_idx),src_idx=1,dst_src_nbr) /   0.260679E+00,   0.246297E+00,   0.493017E+00/
      data (mss_frc_src_lutx( 67,3,src_idx),src_idx=1,dst_src_nbr) /   0.272341E+00,   0.247637E+00,   0.480017E+00/
      data (mss_frc_src_lutx( 68,3,src_idx),src_idx=1,dst_src_nbr) /   0.280003E+00,   0.245209E+00,   0.474783E+00/
      data (mss_frc_src_lutx( 69,3,src_idx),src_idx=1,dst_src_nbr) /   0.289560E+00,   0.244393E+00,   0.466043E+00/
      data (mss_frc_src_lutx( 70,3,src_idx),src_idx=1,dst_src_nbr) /   0.301075E+00,   0.245065E+00,   0.453853E+00/
      data (mss_frc_src_lutx( 71,3,src_idx),src_idx=1,dst_src_nbr) /   0.311629E+00,   0.244790E+00,   0.443577E+00/
      data (mss_frc_src_lutx( 72,3,src_idx),src_idx=1,dst_src_nbr) /   0.321136E+00,   0.243584E+00,   0.435275E+00/
      data (mss_frc_src_lutx( 73,3,src_idx),src_idx=1,dst_src_nbr) /   0.329662E+00,   0.241605E+00,   0.428732E+00/
      data (mss_frc_src_lutx( 74,3,src_idx),src_idx=1,dst_src_nbr) /   0.340295E+00,   0.241096E+00,   0.418609E+00/
      data (mss_frc_src_lutx( 75,3,src_idx),src_idx=1,dst_src_nbr) /   0.349956E+00,   0.239819E+00,   0.410224E+00/
      data (mss_frc_src_lutx( 76,3,src_idx),src_idx=1,dst_src_nbr) /   0.361820E+00,   0.239944E+00,   0.398234E+00/
      data (mss_frc_src_lutx( 77,3,src_idx),src_idx=1,dst_src_nbr) /   0.369564E+00,   0.237293E+00,   0.393147E+00/
      data (mss_frc_src_lutx( 78,3,src_idx),src_idx=1,dst_src_nbr) /   0.379547E+00,   0.236067E+00,   0.384389E+00/
      data (mss_frc_src_lutx( 79,3,src_idx),src_idx=1,dst_src_nbr) /   0.391803E+00,   0.236159E+00,   0.372042E+00/
      data (mss_frc_src_lutx( 80,3,src_idx),src_idx=1,dst_src_nbr) /   0.399957E+00,   0.233736E+00,   0.366307E+00/
      data (mss_frc_src_lutx( 81,3,src_idx),src_idx=1,dst_src_nbr) /   0.410417E+00,   0.232642E+00,   0.356943E+00/
      data (mss_frc_src_lutx( 82,3,src_idx),src_idx=1,dst_src_nbr) /   0.419999E+00,   0.231018E+00,   0.348986E+00/
      data (mss_frc_src_lutx( 83,3,src_idx),src_idx=1,dst_src_nbr) /   0.431918E+00,   0.230624E+00,   0.337455E+00/
      data (mss_frc_src_lutx( 84,3,src_idx),src_idx=1,dst_src_nbr) /   0.439772E+00,   0.228040E+00,   0.332190E+00/
      data (mss_frc_src_lutx( 85,3,src_idx),src_idx=1,dst_src_nbr) /   0.449979E+00,   0.226680E+00,   0.323343E+00/
      data (mss_frc_src_lutx( 86,3,src_idx),src_idx=1,dst_src_nbr) /   0.459374E+00,   0.224901E+00,   0.315725E+00/
      data (mss_frc_src_lutx( 87,3,src_idx),src_idx=1,dst_src_nbr) /   0.467968E+00,   0.222744E+00,   0.309288E+00/
      data (mss_frc_src_lutx( 88,3,src_idx),src_idx=1,dst_src_nbr) /   0.478886E+00,   0.221678E+00,   0.299433E+00/
      data (mss_frc_src_lutx( 89,3,src_idx),src_idx=1,dst_src_nbr) /   0.489048E+00,   0.220245E+00,   0.290708E+00/
      data (mss_frc_src_lutx( 90,3,src_idx),src_idx=1,dst_src_nbr) /   0.498423E+00,   0.218452E+00,   0.283121E+00/
      data (mss_frc_src_lutx( 91,3,src_idx),src_idx=1,dst_src_nbr) /   0.507041E+00,   0.216341E+00,   0.276619E+00/
      data (mss_frc_src_lutx( 92,3,src_idx),src_idx=1,dst_src_nbr) /   0.514960E+00,   0.213971E+00,   0.271071E+00/
      data (mss_frc_src_lutx( 93,3,src_idx),src_idx=1,dst_src_nbr) /   0.525174E+00,   0.212568E+00,   0.262261E+00/
      data (mss_frc_src_lutx( 94,3,src_idx),src_idx=1,dst_src_nbr) /   0.531718E+00,   0.209714E+00,   0.258571E+00/
      data (mss_frc_src_lutx( 95,3,src_idx),src_idx=1,dst_src_nbr) /   0.540553E+00,   0.207806E+00,   0.251641E+00/
      data (mss_frc_src_lutx( 96,3,src_idx),src_idx=1,dst_src_nbr) /   0.551618E+00,   0.206754E+00,   0.241630E+00/
      data (mss_frc_src_lutx( 97,3,src_idx),src_idx=1,dst_src_nbr) /   0.559104E+00,   0.204375E+00,   0.236522E+00/
      data (mss_frc_src_lutx( 98,3,src_idx),src_idx=1,dst_src_nbr) /   0.566011E+00,   0.201841E+00,   0.232153E+00/
      data (mss_frc_src_lutx( 99,3,src_idx),src_idx=1,dst_src_nbr) /   0.575101E+00,   0.200119E+00,   0.224779E+00/
      data (mss_frc_src_lutx(100,3,src_idx),src_idx=1,dst_src_nbr) /   0.583567E+00,   0.198202E+00,   0.218232E+00/
        ! soil type    4  mmd=   690.0  um  sigma=     1.6
            ! first index is wind friction speed in cm/s
                            ! second index is soil type 
                  ! third index is number of source mode
      data (mss_frc_src_lutx(  1,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  2,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  3,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  4,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  5,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  6,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  7,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  8,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx(  9,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 10,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 11,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 12,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 13,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 14,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 15,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 16,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 17,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 18,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 19,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 20,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 21,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 22,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 23,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 24,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 25,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 26,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 27,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 28,4,src_idx),src_idx=1,dst_src_nbr) /   0.000000E+00,   0.000000E+00,   0.000000E+00/
      data (mss_frc_src_lutx( 29,4,src_idx),src_idx=1,dst_src_nbr) /   0.123570E-02,   0.901227E-01,   0.908642E+00/
      data (mss_frc_src_lutx( 30,4,src_idx),src_idx=1,dst_src_nbr) /   0.399021E-02,   0.129109E+00,   0.866901E+00/
      data (mss_frc_src_lutx( 31,4,src_idx),src_idx=1,dst_src_nbr) /   0.759188E-02,   0.153428E+00,   0.838981E+00/
      data (mss_frc_src_lutx( 32,4,src_idx),src_idx=1,dst_src_nbr) /   0.121697E-01,   0.172432E+00,   0.815398E+00/
      data (mss_frc_src_lutx( 33,4,src_idx),src_idx=1,dst_src_nbr) /   0.176021E-01,   0.185754E+00,   0.796645E+00/
      data (mss_frc_src_lutx( 34,4,src_idx),src_idx=1,dst_src_nbr) /   0.242598E-01,   0.198130E+00,   0.777613E+00/
      data (mss_frc_src_lutx( 35,4,src_idx),src_idx=1,dst_src_nbr) /   0.324155E-01,   0.210617E+00,   0.756968E+00/
      data (mss_frc_src_lutx( 36,4,src_idx),src_idx=1,dst_src_nbr) /   0.414584E-01,   0.220157E+00,   0.738383E+00/
      data (mss_frc_src_lutx( 37,4,src_idx),src_idx=1,dst_src_nbr) /   0.509224E-01,   0.227178E+00,   0.721899E+00/
      data (mss_frc_src_lutx( 38,4,src_idx),src_idx=1,dst_src_nbr) /   0.617996E-01,   0.235918E+00,   0.702282E+00/
      data (mss_frc_src_lutx( 39,4,src_idx),src_idx=1,dst_src_nbr) /   0.734209E-01,   0.243068E+00,   0.683510E+00/
      data (mss_frc_src_lutx( 40,4,src_idx),src_idx=1,dst_src_nbr) /   0.858085E-01,   0.248932E+00,   0.665257E+00/
      data (mss_frc_src_lutx( 41,4,src_idx),src_idx=1,dst_src_nbr) /   0.975354E-01,   0.250041E+00,   0.652423E+00/
      data (mss_frc_src_lutx( 42,4,src_idx),src_idx=1,dst_src_nbr) /   0.111387E+00,   0.254137E+00,   0.634476E+00/
      data (mss_frc_src_lutx( 43,4,src_idx),src_idx=1,dst_src_nbr) /   0.126119E+00,   0.257678E+00,   0.616202E+00/
      data (mss_frc_src_lutx( 44,4,src_idx),src_idx=1,dst_src_nbr) /   0.141793E+00,   0.260840E+00,   0.597365E+00/
      data (mss_frc_src_lutx( 45,4,src_idx),src_idx=1,dst_src_nbr) /   0.158473E+00,   0.263765E+00,   0.577760E+00/
      data (mss_frc_src_lutx( 46,4,src_idx),src_idx=1,dst_src_nbr) /   0.174066E+00,   0.263315E+00,   0.562617E+00/
      data (mss_frc_src_lutx( 47,4,src_idx),src_idx=1,dst_src_nbr) /   0.192861E+00,   0.266222E+00,   0.540919E+00/
      data (mss_frc_src_lutx( 48,4,src_idx),src_idx=1,dst_src_nbr) /   0.210439E+00,   0.266071E+00,   0.523486E+00/
      data (mss_frc_src_lutx( 49,4,src_idx),src_idx=1,dst_src_nbr) /   0.229116E+00,   0.266270E+00,   0.504613E+00/
      data (mss_frc_src_lutx( 50,4,src_idx),src_idx=1,dst_src_nbr) /   0.246277E+00,   0.263933E+00,   0.489785E+00/
      data (mss_frc_src_lutx( 51,4,src_idx),src_idx=1,dst_src_nbr) /   0.267258E+00,   0.264915E+00,   0.467824E+00/
      data (mss_frc_src_lutx( 52,4,src_idx),src_idx=1,dst_src_nbr) /   0.283840E+00,   0.260973E+00,   0.455184E+00/
      data (mss_frc_src_lutx( 53,4,src_idx),src_idx=1,dst_src_nbr) /   0.304402E+00,   0.260294E+00,   0.435303E+00/
      data (mss_frc_src_lutx( 54,4,src_idx),src_idx=1,dst_src_nbr) /   0.323289E+00,   0.257754E+00,   0.418952E+00/
      data (mss_frc_src_lutx( 55,4,src_idx),src_idx=1,dst_src_nbr) /   0.343373E+00,   0.255856E+00,   0.400766E+00/
      data (mss_frc_src_lutx( 56,4,src_idx),src_idx=1,dst_src_nbr) /   0.361698E+00,   0.252459E+00,   0.385842E+00/
      data (mss_frc_src_lutx( 57,4,src_idx),src_idx=1,dst_src_nbr) /   0.381066E+00,   0.249850E+00,   0.369083E+00/
      data (mss_frc_src_lutx( 58,4,src_idx),src_idx=1,dst_src_nbr) /   0.398427E+00,   0.246098E+00,   0.355476E+00/
      data (mss_frc_src_lutx( 59,4,src_idx),src_idx=1,dst_src_nbr) /   0.416878E+00,   0.243187E+00,   0.339937E+00/
      data (mss_frc_src_lutx( 60,4,src_idx),src_idx=1,dst_src_nbr) /   0.436483E+00,   0.241018E+00,   0.322496E+00/
      data (mss_frc_src_lutx( 61,4,src_idx),src_idx=1,dst_src_nbr) /   0.451137E+00,   0.236282E+00,   0.312580E+00/
      data (mss_frc_src_lutx( 62,4,src_idx),src_idx=1,dst_src_nbr) /   0.470061E+00,   0.233942E+00,   0.295992E+00/
      data (mss_frc_src_lutx( 63,4,src_idx),src_idx=1,dst_src_nbr) /   0.484215E+00,   0.229383E+00,   0.286401E+00/
      data (mss_frc_src_lutx( 64,4,src_idx),src_idx=1,dst_src_nbr) /   0.502614E+00,   0.226982E+00,   0.270403E+00/
      data (mss_frc_src_lutx( 65,4,src_idx),src_idx=1,dst_src_nbr) /   0.516397E+00,   0.222636E+00,   0.260967E+00/
      data (mss_frc_src_lutx( 66,4,src_idx),src_idx=1,dst_src_nbr) /   0.531509E+00,   0.219051E+00,   0.249440E+00/
      data (mss_frc_src_lutx( 67,4,src_idx),src_idx=1,dst_src_nbr) /   0.547899E+00,   0.216119E+00,   0.235984E+00/
      data (mss_frc_src_lutx( 68,4,src_idx),src_idx=1,dst_src_nbr) /   0.560003E+00,   0.211658E+00,   0.228343E+00/
      data (mss_frc_src_lutx( 69,4,src_idx),src_idx=1,dst_src_nbr) /   0.573467E+00,   0.207908E+00,   0.218627E+00/
      data (mss_frc_src_lutx( 70,4,src_idx),src_idx=1,dst_src_nbr) /   0.588184E+00,   0.204748E+00,   0.207069E+00/
      data (mss_frc_src_lutx( 71,4,src_idx),src_idx=1,dst_src_nbr) /   0.601560E+00,   0.201257E+00,   0.197182E+00/
      data (mss_frc_src_lutx( 72,4,src_idx),src_idx=1,dst_src_nbr) /   0.613645E+00,   0.197488E+00,   0.188868E+00/
      data (mss_frc_src_lutx( 73,4,src_idx),src_idx=1,dst_src_nbr) /   0.624582E+00,   0.193528E+00,   0.181888E+00/
      data (mss_frc_src_lutx( 74,4,src_idx),src_idx=1,dst_src_nbr) /   0.636765E+00,   0.190108E+00,   0.173124E+00/
      data (mss_frc_src_lutx( 75,4,src_idx),src_idx=1,dst_src_nbr) /   0.647832E+00,   0.186505E+00,   0.165664E+00/
      data (mss_frc_src_lutx( 76,4,src_idx),src_idx=1,dst_src_nbr) /   0.660055E+00,   0.183369E+00,   0.156571E+00/
      data (mss_frc_src_lutx( 77,4,src_idx),src_idx=1,dst_src_nbr) /   0.669091E+00,   0.179497E+00,   0.151411E+00/
      data (mss_frc_src_lutx( 78,4,src_idx),src_idx=1,dst_src_nbr) /   0.679319E+00,   0.176099E+00,   0.144577E+00/
      data (mss_frc_src_lutx( 79,4,src_idx),src_idx=1,dst_src_nbr) /   0.690623E+00,   0.173105E+00,   0.136273E+00/
      data (mss_frc_src_lutx( 80,4,src_idx),src_idx=1,dst_src_nbr) /   0.699048E+00,   0.169527E+00,   0.131426E+00/
      data (mss_frc_src_lutx( 81,4,src_idx),src_idx=1,dst_src_nbr) /   0.708572E+00,   0.166351E+00,   0.125081E+00/
      data (mss_frc_src_lutx( 82,4,src_idx),src_idx=1,dst_src_nbr) /   0.717282E+00,   0.163113E+00,   0.119606E+00/
      data (mss_frc_src_lutx( 83,4,src_idx),src_idx=1,dst_src_nbr) /   0.726980E+00,   0.160220E+00,   0.112803E+00/
      data (mss_frc_src_lutx( 84,4,src_idx),src_idx=1,dst_src_nbr) /   0.734208E+00,   0.156905E+00,   0.108890E+00/
      data (mss_frc_src_lutx( 85,4,src_idx),src_idx=1,dst_src_nbr) /   0.742429E+00,   0.153926E+00,   0.103649E+00/
      data (mss_frc_src_lutx( 86,4,src_idx),src_idx=1,dst_src_nbr) /   0.749983E+00,   0.150926E+00,   0.990928E-01/
      data (mss_frc_src_lutx( 87,4,src_idx),src_idx=1,dst_src_nbr) /   0.756928E+00,   0.147922E+00,   0.951542E-01/
      data (mss_frc_src_lutx( 88,4,src_idx),src_idx=1,dst_src_nbr) /   0.764737E+00,   0.145193E+00,   0.900722E-01/
      data (mss_frc_src_lutx( 89,4,src_idx),src_idx=1,dst_src_nbr) /   0.771947E+00,   0.142454E+00,   0.855983E-01/
      data (mss_frc_src_lutx( 90,4,src_idx),src_idx=1,dst_src_nbr) /   0.778593E+00,   0.139713E+00,   0.816947E-01/
      data (mss_frc_src_lutx( 91,4,src_idx),src_idx=1,dst_src_nbr) /   0.784720E+00,   0.136980E+00,   0.783003E-01/
      data (mss_frc_src_lutx( 92,4,src_idx),src_idx=1,dst_src_nbr) /   0.790389E+00,   0.134270E+00,   0.753415E-01/
      data (mss_frc_src_lutx( 93,4,src_idx),src_idx=1,dst_src_nbr) /   0.796808E+00,   0.131782E+00,   0.714112E-01/
      data (mss_frc_src_lutx( 94,4,src_idx),src_idx=1,dst_src_nbr) /   0.801623E+00,   0.129122E+00,   0.692559E-01/
      data (mss_frc_src_lutx( 95,4,src_idx),src_idx=1,dst_src_nbr) /   0.807175E+00,   0.126674E+00,   0.661509E-01/
      data (mss_frc_src_lutx( 96,4,src_idx),src_idx=1,dst_src_nbr) /   0.813374E+00,   0.124410E+00,   0.622161E-01/
      data (mss_frc_src_lutx( 97,4,src_idx),src_idx=1,dst_src_nbr) /   0.818118E+00,   0.122005E+00,   0.598770E-01/
      data (mss_frc_src_lutx( 98,4,src_idx),src_idx=1,dst_src_nbr) /   0.822543E+00,   0.119637E+00,   0.578235E-01/
      data (mss_frc_src_lutx( 99,4,src_idx),src_idx=1,dst_src_nbr) /   0.827582E+00,   0.117438E+00,   0.549798E-01/
      data (mss_frc_src_lutx(100,4,src_idx),src_idx=1,dst_src_nbr) /   0.832278E+00,   0.115264E+00,   0.524592E-01/
                      !alpha values, vert/hor flux [m-1]
                             ! first index is u* in cm/s
                    ! second index soil type (see above)
      data dst_slt_flx_rat_ttl_lutx (  1,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  2,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  3,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  4,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  5,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  6,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  7,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  8,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  9,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 10,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 11,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 12,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 13,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 14,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 15,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 16,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 17,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 18,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 19,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 20,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 21,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 22,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 23,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 24,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 25,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 26,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 27,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 28,1) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 29,1) /   0.420805E-06/
      data dst_slt_flx_rat_ttl_lutx ( 30,1) /   0.922056E-06/
      data dst_slt_flx_rat_ttl_lutx ( 31,1) /   0.135488E-05/
      data dst_slt_flx_rat_ttl_lutx ( 32,1) /   0.172794E-05/
      data dst_slt_flx_rat_ttl_lutx ( 33,1) /   0.209171E-05/
      data dst_slt_flx_rat_ttl_lutx ( 34,1) /   0.241046E-05/
      data dst_slt_flx_rat_ttl_lutx ( 35,1) /   0.267653E-05/
      data dst_slt_flx_rat_ttl_lutx ( 36,1) /   0.294595E-05/
      data dst_slt_flx_rat_ttl_lutx ( 37,1) /   0.322349E-05/
      data dst_slt_flx_rat_ttl_lutx ( 38,1) /   0.343558E-05/
      data dst_slt_flx_rat_ttl_lutx ( 39,1) /   0.364810E-05/
      data dst_slt_flx_rat_ttl_lutx ( 40,1) /   0.385997E-05/
      data dst_slt_flx_rat_ttl_lutx ( 41,1) /   0.415984E-05/
      data dst_slt_flx_rat_ttl_lutx ( 42,1) /   0.436797E-05/
      data dst_slt_flx_rat_ttl_lutx ( 43,1) /   0.456631E-05/
      data dst_slt_flx_rat_ttl_lutx ( 44,1) /   0.475042E-05/
      data dst_slt_flx_rat_ttl_lutx ( 45,1) /   0.491564E-05/
      data dst_slt_flx_rat_ttl_lutx ( 46,1) /   0.516647E-05/
      data dst_slt_flx_rat_ttl_lutx ( 47,1) /   0.528040E-05/
      data dst_slt_flx_rat_ttl_lutx ( 48,1) /   0.547704E-05/
      data dst_slt_flx_rat_ttl_lutx ( 49,1) /   0.563793E-05/
      data dst_slt_flx_rat_ttl_lutx ( 50,1) /   0.588353E-05/
      data dst_slt_flx_rat_ttl_lutx ( 51,1) /   0.595894E-05/
      data dst_slt_flx_rat_ttl_lutx ( 52,1) /   0.624569E-05/
      data dst_slt_flx_rat_ttl_lutx ( 53,1) /   0.634874E-05/
      data dst_slt_flx_rat_ttl_lutx ( 54,1) /   0.652916E-05/
      data dst_slt_flx_rat_ttl_lutx ( 55,1) /   0.664894E-05/
      data dst_slt_flx_rat_ttl_lutx ( 56,1) /   0.684435E-05/
      data dst_slt_flx_rat_ttl_lutx ( 57,1) /   0.697150E-05/
      data dst_slt_flx_rat_ttl_lutx ( 58,1) /   0.717333E-05/
      data dst_slt_flx_rat_ttl_lutx ( 59,1) /   0.729813E-05/
      data dst_slt_flx_rat_ttl_lutx ( 60,1) /   0.733869E-05/
      data dst_slt_flx_rat_ttl_lutx ( 61,1) /   0.760881E-05/
      data dst_slt_flx_rat_ttl_lutx ( 62,1) /   0.762952E-05/
      data dst_slt_flx_rat_ttl_lutx ( 63,1) /   0.788181E-05/
      data dst_slt_flx_rat_ttl_lutx ( 64,1) /   0.787143E-05/
      data dst_slt_flx_rat_ttl_lutx ( 65,1) /   0.809533E-05/
      data dst_slt_flx_rat_ttl_lutx ( 66,1) /   0.821617E-05/
      data dst_slt_flx_rat_ttl_lutx ( 67,1) /   0.822717E-05/
      data dst_slt_flx_rat_ttl_lutx ( 68,1) /   0.847784E-05/
      data dst_slt_flx_rat_ttl_lutx ( 69,1) /   0.861551E-05/
      data dst_slt_flx_rat_ttl_lutx ( 70,1) /   0.863526E-05/
      data dst_slt_flx_rat_ttl_lutx ( 71,1) /   0.871245E-05/
      data dst_slt_flx_rat_ttl_lutx ( 72,1) /   0.885072E-05/
      data dst_slt_flx_rat_ttl_lutx ( 73,1) /   0.904722E-05/
      data dst_slt_flx_rat_ttl_lutx ( 74,1) /   0.911484E-05/
      data dst_slt_flx_rat_ttl_lutx ( 75,1) /   0.923894E-05/
      data dst_slt_flx_rat_ttl_lutx ( 76,1) /   0.922638E-05/
      data dst_slt_flx_rat_ttl_lutx ( 77,1) /   0.946187E-05/
      data dst_slt_flx_rat_ttl_lutx ( 78,1) /   0.955736E-05/
      data dst_slt_flx_rat_ttl_lutx ( 79,1) /   0.950769E-05/
      data dst_slt_flx_rat_ttl_lutx ( 80,1) /   0.970787E-05/
      data dst_slt_flx_rat_ttl_lutx ( 81,1) /   0.975978E-05/
      data dst_slt_flx_rat_ttl_lutx ( 82,1) /   0.986164E-05/
      data dst_slt_flx_rat_ttl_lutx ( 83,1) /   0.980835E-05/
      data dst_slt_flx_rat_ttl_lutx ( 84,1) /   0.100091E-04/
      data dst_slt_flx_rat_ttl_lutx ( 85,1) /   0.100518E-04/
      data dst_slt_flx_rat_ttl_lutx ( 86,1) /   0.101408E-04/
      data dst_slt_flx_rat_ttl_lutx ( 87,1) /   0.102764E-04/
      data dst_slt_flx_rat_ttl_lutx ( 88,1) /   0.102476E-04/
      data dst_slt_flx_rat_ttl_lutx ( 89,1) /   0.102606E-04/
      data dst_slt_flx_rat_ttl_lutx ( 90,1) /   0.103180E-04/
      data dst_slt_flx_rat_ttl_lutx ( 91,1) /   0.104197E-04/
      data dst_slt_flx_rat_ttl_lutx ( 92,1) /   0.105634E-04/
      data dst_slt_flx_rat_ttl_lutx ( 93,1) /   0.105327E-04/
      data dst_slt_flx_rat_ttl_lutx ( 94,1) /   0.107607E-04/
      data dst_slt_flx_rat_ttl_lutx ( 95,1) /   0.108114E-04/
      data dst_slt_flx_rat_ttl_lutx ( 96,1) /   0.106808E-04/
      data dst_slt_flx_rat_ttl_lutx ( 97,1) /   0.108107E-04/
      data dst_slt_flx_rat_ttl_lutx ( 98,1) /   0.109785E-04/
      data dst_slt_flx_rat_ttl_lutx ( 99,1) /   0.109621E-04/
      data dst_slt_flx_rat_ttl_lutx (100,1) /   0.109830E-04/
                             ! first index is u* in cm/s
                    ! second index soil type (see above)
      data dst_slt_flx_rat_ttl_lutx (  1,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  2,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  3,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  4,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  5,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  6,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  7,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  8,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  9,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 10,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 11,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 12,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 13,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 14,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 15,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 16,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 17,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 18,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 19,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 20,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 21,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 22,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 23,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 24,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 25,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 26,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 27,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 28,2) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 29,2) /   0.243914E-05/
      data dst_slt_flx_rat_ttl_lutx ( 30,2) /   0.498042E-05/
      data dst_slt_flx_rat_ttl_lutx ( 31,2) /   0.686192E-05/
      data dst_slt_flx_rat_ttl_lutx ( 32,2) /   0.825433E-05/
      data dst_slt_flx_rat_ttl_lutx ( 33,2) /   0.946807E-05/
      data dst_slt_flx_rat_ttl_lutx ( 34,2) /   0.103858E-04/
      data dst_slt_flx_rat_ttl_lutx ( 35,2) /   0.110211E-04/
      data dst_slt_flx_rat_ttl_lutx ( 36,2) /   0.116261E-04/
      data dst_slt_flx_rat_ttl_lutx ( 37,2) /   0.122221E-04/
      data dst_slt_flx_rat_ttl_lutx ( 38,2) /   0.125491E-04/
      data dst_slt_flx_rat_ttl_lutx ( 39,2) /   0.128585E-04/
      data dst_slt_flx_rat_ttl_lutx ( 40,2) /   0.131472E-04/
      data dst_slt_flx_rat_ttl_lutx ( 41,2) /   0.136957E-04/
      data dst_slt_flx_rat_ttl_lutx ( 42,2) /   0.139294E-04/
      data dst_slt_flx_rat_ttl_lutx ( 43,2) /   0.141190E-04/
      data dst_slt_flx_rat_ttl_lutx ( 44,2) /   0.142553E-04/
      data dst_slt_flx_rat_ttl_lutx ( 45,2) /   0.143295E-04/
      data dst_slt_flx_rat_ttl_lutx ( 46,2) /   0.146277E-04/
      data dst_slt_flx_rat_ttl_lutx ( 47,2) /   0.145484E-04/
      data dst_slt_flx_rat_ttl_lutx ( 48,2) /   0.146810E-04/
      data dst_slt_flx_rat_ttl_lutx ( 49,2) /   0.147147E-04/
      data dst_slt_flx_rat_ttl_lutx ( 50,2) /   0.149470E-04/
      data dst_slt_flx_rat_ttl_lutx ( 51,2) /   0.147652E-04/
      data dst_slt_flx_rat_ttl_lutx ( 52,2) /   0.150706E-04/
      data dst_slt_flx_rat_ttl_lutx ( 53,2) /   0.149490E-04/
      data dst_slt_flx_rat_ttl_lutx ( 54,2) /   0.149968E-04/
      data dst_slt_flx_rat_ttl_lutx ( 55,2) /   0.149106E-04/
      data dst_slt_flx_rat_ttl_lutx ( 56,2) /   0.149799E-04/
      data dst_slt_flx_rat_ttl_lutx ( 57,2) /   0.149052E-04/
      data dst_slt_flx_rat_ttl_lutx ( 58,2) /   0.149759E-04/
      data dst_slt_flx_rat_ttl_lutx ( 59,2) /   0.148923E-04/
      data dst_slt_flx_rat_ttl_lutx ( 60,2) /   0.146518E-04/
      data dst_slt_flx_rat_ttl_lutx ( 61,2) /   0.148359E-04/
      data dst_slt_flx_rat_ttl_lutx ( 62,2) /   0.145648E-04/
      data dst_slt_flx_rat_ttl_lutx ( 63,2) /   0.147037E-04/
      data dst_slt_flx_rat_ttl_lutx ( 64,2) /   0.143870E-04/
      data dst_slt_flx_rat_ttl_lutx ( 65,2) /   0.144691E-04/
      data dst_slt_flx_rat_ttl_lutx ( 66,2) /   0.143766E-04/
      data dst_slt_flx_rat_ttl_lutx ( 67,2) /   0.141103E-04/
      data dst_slt_flx_rat_ttl_lutx ( 68,2) /   0.142235E-04/
      data dst_slt_flx_rat_ttl_lutx ( 69,2) /   0.141567E-04/
      data dst_slt_flx_rat_ttl_lutx ( 70,2) /   0.139142E-04/
      data dst_slt_flx_rat_ttl_lutx ( 71,2) /   0.137619E-04/
      data dst_slt_flx_rat_ttl_lutx ( 72,2) /   0.136993E-04/
      data dst_slt_flx_rat_ttl_lutx ( 73,2) /   0.137168E-04/
      data dst_slt_flx_rat_ttl_lutx ( 74,2) /   0.135548E-04/
      data dst_slt_flx_rat_ttl_lutx ( 75,2) /   0.134714E-04/
      data dst_slt_flx_rat_ttl_lutx ( 76,2) /   0.132095E-04/
      data dst_slt_flx_rat_ttl_lutx ( 77,2) /   0.132727E-04/
      data dst_slt_flx_rat_ttl_lutx ( 78,2) /   0.131544E-04/
      data dst_slt_flx_rat_ttl_lutx ( 79,2) /   0.128593E-04/
      data dst_slt_flx_rat_ttl_lutx ( 80,2) /   0.128741E-04/
      data dst_slt_flx_rat_ttl_lutx ( 81,2) /   0.127101E-04/
      data dst_slt_flx_rat_ttl_lutx ( 82,2) /   0.126074E-04/
      data dst_slt_flx_rat_ttl_lutx ( 83,2) /   0.123295E-04/
      data dst_slt_flx_rat_ttl_lutx ( 84,2) /   0.123427E-04/
      data dst_slt_flx_rat_ttl_lutx ( 85,2) /   0.121800E-04/
      data dst_slt_flx_rat_ttl_lutx ( 86,2) /   0.120702E-04/
      data dst_slt_flx_rat_ttl_lutx ( 87,2) /   0.120111E-04/
      data dst_slt_flx_rat_ttl_lutx ( 88,2) /   0.117819E-04/
      data dst_slt_flx_rat_ttl_lutx ( 89,2) /   0.116010E-04/
      data dst_slt_flx_rat_ttl_lutx ( 90,2) /   0.114686E-04/
      data dst_slt_flx_rat_ttl_lutx ( 91,2) /   0.113818E-04/
      data dst_slt_flx_rat_ttl_lutx ( 92,2) /   0.113361E-04/
      data dst_slt_flx_rat_ttl_lutx ( 93,2) /   0.111257E-04/
      data dst_slt_flx_rat_ttl_lutx ( 94,2) /   0.111597E-04/
      data dst_slt_flx_rat_ttl_lutx ( 95,2) /   0.110292E-04/
      data dst_slt_flx_rat_ttl_lutx ( 96,2) /   0.107397E-04/
      data dst_slt_flx_rat_ttl_lutx ( 97,2) /   0.106863E-04/
      data dst_slt_flx_rat_ttl_lutx ( 98,2) /   0.106654E-04/
      data dst_slt_flx_rat_ttl_lutx ( 99,2) /   0.104876E-04/
      data dst_slt_flx_rat_ttl_lutx (100,2) /   0.103448E-04/
                             ! first index is u* in cm/s
                    ! second index soil type (see above)
      data dst_slt_flx_rat_ttl_lutx (  1,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  2,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  3,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  4,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  5,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  6,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  7,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  8,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  9,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 10,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 11,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 12,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 13,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 14,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 15,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 16,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 17,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 18,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 19,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 20,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 21,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 22,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 23,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 24,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 25,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 26,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 27,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 28,3) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 29,3) /   0.163314E-04/
      data dst_slt_flx_rat_ttl_lutx ( 30,3) /   0.250721E-04/
      data dst_slt_flx_rat_ttl_lutx ( 31,3) /   0.264813E-04/
      data dst_slt_flx_rat_ttl_lutx ( 32,3) /   0.249222E-04/
      data dst_slt_flx_rat_ttl_lutx ( 33,3) /   0.227558E-04/
      data dst_slt_flx_rat_ttl_lutx ( 34,3) /   0.202400E-04/
      data dst_slt_flx_rat_ttl_lutx ( 35,3) /   0.177260E-04/
      data dst_slt_flx_rat_ttl_lutx ( 36,3) /   0.156852E-04/
      data dst_slt_flx_rat_ttl_lutx ( 37,3) /   0.140726E-04/
      data dst_slt_flx_rat_ttl_lutx ( 38,3) /   0.125407E-04/
      data dst_slt_flx_rat_ttl_lutx ( 39,3) /   0.112776E-04/
      data dst_slt_flx_rat_ttl_lutx ( 40,3) /   0.102156E-04/
      data dst_slt_flx_rat_ttl_lutx ( 41,3) /   0.946706E-05/
      data dst_slt_flx_rat_ttl_lutx ( 42,3) /   0.866136E-05/
      data dst_slt_flx_rat_ttl_lutx ( 43,3) /   0.795038E-05/
      data dst_slt_flx_rat_ttl_lutx ( 44,3) /   0.731466E-05/
      data dst_slt_flx_rat_ttl_lutx ( 45,3) /   0.673992E-05/
      data dst_slt_flx_rat_ttl_lutx ( 46,3) /   0.631249E-05/
      data dst_slt_flx_rat_ttl_lutx ( 47,3) /   0.581958E-05/
      data dst_slt_flx_rat_ttl_lutx ( 48,3) /   0.544525E-05/
      data dst_slt_flx_rat_ttl_lutx ( 49,3) /   0.508712E-05/
      data dst_slt_flx_rat_ttl_lutx ( 50,3) /   0.481472E-05/
      data dst_slt_flx_rat_ttl_lutx ( 51,3) /   0.447938E-05/
      data dst_slt_flx_rat_ttl_lutx ( 52,3) /   0.427808E-05/
      data dst_slt_flx_rat_ttl_lutx ( 53,3) /   0.401460E-05/
      data dst_slt_flx_rat_ttl_lutx ( 54,3) /   0.380684E-05/
      data dst_slt_flx_rat_ttl_lutx ( 55,3) /   0.359634E-05/
      data dst_slt_flx_rat_ttl_lutx ( 56,3) /   0.342906E-05/
      data dst_slt_flx_rat_ttl_lutx ( 57,3) /   0.325606E-05/
      data dst_slt_flx_rat_ttl_lutx ( 58,3) /   0.311836E-05/
      data dst_slt_flx_rat_ttl_lutx ( 59,3) /   0.297257E-05/
      data dst_slt_flx_rat_ttl_lutx ( 60,3) /   0.282024E-05/
      data dst_slt_flx_rat_ttl_lutx ( 61,3) /   0.272899E-05/
      data dst_slt_flx_rat_ttl_lutx ( 62,3) /   0.259519E-05/
      data dst_slt_flx_rat_ttl_lutx ( 63,3) /   0.251396E-05/
      data dst_slt_flx_rat_ttl_lutx ( 64,3) /   0.239403E-05/
      data dst_slt_flx_rat_ttl_lutx ( 65,3) /   0.232004E-05/
      data dst_slt_flx_rat_ttl_lutx ( 66,3) /   0.223554E-05/
      data dst_slt_flx_rat_ttl_lutx ( 67,3) /   0.214228E-05/
      data dst_slt_flx_rat_ttl_lutx ( 68,3) /   0.208589E-05/
      data dst_slt_flx_rat_ttl_lutx ( 69,3) /   0.201906E-05/
      data dst_slt_flx_rat_ttl_lutx ( 70,3) /   0.194364E-05/
      data dst_slt_flx_rat_ttl_lutx ( 71,3) /   0.187945E-05/
      data dst_slt_flx_rat_ttl_lutx ( 72,3) /   0.182528E-05/
      data dst_slt_flx_rat_ttl_lutx ( 73,3) /   0.177943E-05/
      data dst_slt_flx_rat_ttl_lutx ( 74,3) /   0.172506E-05/
      data dst_slt_flx_rat_ttl_lutx ( 75,3) /   0.167856E-05/
      data dst_slt_flx_rat_ttl_lutx ( 76,3) /   0.162454E-05/
      data dst_slt_flx_rat_ttl_lutx ( 77,3) /   0.159144E-05/
      data dst_slt_flx_rat_ttl_lutx ( 78,3) /   0.155044E-05/
      data dst_slt_flx_rat_ttl_lutx ( 79,3) /   0.150272E-05/
      data dst_slt_flx_rat_ttl_lutx ( 80,3) /   0.147282E-05/
      data dst_slt_flx_rat_ttl_lutx ( 81,3) /   0.143595E-05/
      data dst_slt_flx_rat_ttl_lutx ( 82,3) /   0.140381E-05/
      data dst_slt_flx_rat_ttl_lutx ( 83,3) /   0.136564E-05/
      data dst_slt_flx_rat_ttl_lutx ( 84,3) /   0.134177E-05/
      data dst_slt_flx_rat_ttl_lutx ( 85,3) /   0.131182E-05/
      data dst_slt_flx_rat_ttl_lutx ( 86,3) /   0.128545E-05/
      data dst_slt_flx_rat_ttl_lutx ( 87,3) /   0.126226E-05/
      data dst_slt_flx_rat_ttl_lutx ( 88,3) /   0.123388E-05/
      data dst_slt_flx_rat_ttl_lutx ( 89,3) /   0.120860E-05/
      data dst_slt_flx_rat_ttl_lutx ( 90,3) /   0.118621E-05/
      data dst_slt_flx_rat_ttl_lutx ( 91,3) /   0.116636E-05/
      data dst_slt_flx_rat_ttl_lutx ( 92,3) /   0.114872E-05/
      data dst_slt_flx_rat_ttl_lutx ( 93,3) /   0.112666E-05/
      data dst_slt_flx_rat_ttl_lutx ( 94,3) /   0.111305E-05/
      data dst_slt_flx_rat_ttl_lutx ( 95,3) /   0.109511E-05/
      data dst_slt_flx_rat_ttl_lutx ( 96,3) /   0.107337E-05/
      data dst_slt_flx_rat_ttl_lutx ( 97,3) /   0.105920E-05/
      data dst_slt_flx_rat_ttl_lutx ( 98,3) /   0.104648E-05/
      data dst_slt_flx_rat_ttl_lutx ( 99,3) /   0.103014E-05/
      data dst_slt_flx_rat_ttl_lutx (100,3) /   0.101537E-05/
                             ! first index is u* in cm/s
                    ! second index soil type (see above)
      data dst_slt_flx_rat_ttl_lutx (  1,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  2,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  3,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  4,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  5,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  6,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  7,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  8,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx (  9,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 10,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 11,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 12,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 13,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 14,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 15,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 16,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 17,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 18,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 19,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 20,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 21,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 22,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 23,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 24,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 25,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 26,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 27,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 28,4) /   0.000000E+00/
      data dst_slt_flx_rat_ttl_lutx ( 29,4) /   0.232842E-04/
      data dst_slt_flx_rat_ttl_lutx ( 30,4) /   0.320103E-04/
      data dst_slt_flx_rat_ttl_lutx ( 31,4) /   0.304793E-04/
      data dst_slt_flx_rat_ttl_lutx ( 32,4) /   0.260441E-04/
      data dst_slt_flx_rat_ttl_lutx ( 33,4) /   0.217250E-04/
      data dst_slt_flx_rat_ttl_lutx ( 34,4) /   0.177746E-04/
      data dst_slt_flx_rat_ttl_lutx ( 35,4) /   0.144189E-04/
      data dst_slt_flx_rat_ttl_lutx ( 36,4) /   0.119114E-04/
      data dst_slt_flx_rat_ttl_lutx ( 37,4) /   0.100731E-04/
      data dst_slt_flx_rat_ttl_lutx ( 38,4) /   0.853273E-05/
      data dst_slt_flx_rat_ttl_lutx ( 39,4) /   0.733357E-05/
      data dst_slt_flx_rat_ttl_lutx ( 40,4) /   0.637719E-05/
      data dst_slt_flx_rat_ttl_lutx ( 41,4) /   0.568271E-05/
      data dst_slt_flx_rat_ttl_lutx ( 42,4) /   0.502753E-05/
      data dst_slt_flx_rat_ttl_lutx ( 43,4) /   0.447766E-05/
      data dst_slt_flx_rat_ttl_lutx ( 44,4) /   0.401032E-05/
      data dst_slt_flx_rat_ttl_lutx ( 45,4) /   0.360892E-05/
      data dst_slt_flx_rat_ttl_lutx ( 46,4) /   0.330152E-05/
      data dst_slt_flx_rat_ttl_lutx ( 47,4) /   0.299197E-05/
      data dst_slt_flx_rat_ttl_lutx ( 48,4) /   0.275160E-05/
      data dst_slt_flx_rat_ttl_lutx ( 49,4) /   0.253483E-05/
      data dst_slt_flx_rat_ttl_lutx ( 50,4) /   0.236427E-05/
      data dst_slt_flx_rat_ttl_lutx ( 51,4) /   0.218351E-05/
      data dst_slt_flx_rat_ttl_lutx ( 52,4) /   0.205994E-05/
      data dst_slt_flx_rat_ttl_lutx ( 53,4) /   0.192405E-05/
      data dst_slt_flx_rat_ttl_lutx ( 54,4) /   0.181434E-05/
      data dst_slt_flx_rat_ttl_lutx ( 55,4) /   0.171046E-05/
      data dst_slt_flx_rat_ttl_lutx ( 56,4) /   0.162569E-05/
      data dst_slt_flx_rat_ttl_lutx ( 57,4) /   0.154464E-05/
      data dst_slt_flx_rat_ttl_lutx ( 58,4) /   0.147868E-05/
      data dst_slt_flx_rat_ttl_lutx ( 59,4) /   0.141436E-05/
      data dst_slt_flx_rat_ttl_lutx ( 60,4) /   0.135180E-05/
      data dst_slt_flx_rat_ttl_lutx ( 61,4) /   0.130874E-05/
      data dst_slt_flx_rat_ttl_lutx ( 62,4) /   0.125678E-05/
      data dst_slt_flx_rat_ttl_lutx ( 63,4) /   0.122068E-05/
      data dst_slt_flx_rat_ttl_lutx ( 64,4) /   0.117656E-05/
      data dst_slt_flx_rat_ttl_lutx ( 65,4) /   0.114565E-05/
      data dst_slt_flx_rat_ttl_lutx ( 66,4) /   0.111352E-05/
      data dst_slt_flx_rat_ttl_lutx ( 67,4) /   0.108060E-05/
      data dst_slt_flx_rat_ttl_lutx ( 68,4) /   0.105759E-05/
      data dst_slt_flx_rat_ttl_lutx ( 69,4) /   0.103308E-05/
      data dst_slt_flx_rat_ttl_lutx ( 70,4) /   0.100752E-05/
      data dst_slt_flx_rat_ttl_lutx ( 71,4) /   0.985373E-06/
      data dst_slt_flx_rat_ttl_lutx ( 72,4) /   0.966200E-06/
      data dst_slt_flx_rat_ttl_lutx ( 73,4) /   0.949491E-06/
      data dst_slt_flx_rat_ttl_lutx ( 74,4) /   0.931520E-06/
      data dst_slt_flx_rat_ttl_lutx ( 75,4) /   0.915784E-06/
      data dst_slt_flx_rat_ttl_lutx ( 76,4) /   0.898985E-06/
      data dst_slt_flx_rat_ttl_lutx ( 77,4) /   0.886994E-06/
      data dst_slt_flx_rat_ttl_lutx ( 78,4) /   0.873781E-06/
      data dst_slt_flx_rat_ttl_lutx ( 79,4) /   0.859605E-06/
      data dst_slt_flx_rat_ttl_lutx ( 80,4) /   0.849362E-06/
      data dst_slt_flx_rat_ttl_lutx ( 81,4) /   0.838055E-06/
      data dst_slt_flx_rat_ttl_lutx ( 82,4) /   0.827979E-06/
      data dst_slt_flx_rat_ttl_lutx ( 83,4) /   0.817027E-06/
      data dst_slt_flx_rat_ttl_lutx ( 84,4) /   0.809070E-06/
      data dst_slt_flx_rat_ttl_lutx ( 85,4) /   0.800193E-06/
      data dst_slt_flx_rat_ttl_lutx ( 86,4) /   0.792209E-06/
      data dst_slt_flx_rat_ttl_lutx ( 87,4) /   0.785012E-06/
      data dst_slt_flx_rat_ttl_lutx ( 88,4) /   0.777061E-06/
      data dst_slt_flx_rat_ttl_lutx ( 89,4) /   0.769864E-06/
      data dst_slt_flx_rat_ttl_lutx ( 90,4) /   0.763351E-06/
      data dst_slt_flx_rat_ttl_lutx ( 91,4) /   0.757447E-06/
      data dst_slt_flx_rat_ttl_lutx ( 92,4) /   0.752065E-06/
      data dst_slt_flx_rat_ttl_lutx ( 93,4) /   0.746056E-06/
      data dst_slt_flx_rat_ttl_lutx ( 94,4) /   0.741620E-06/
      data dst_slt_flx_rat_ttl_lutx ( 95,4) /   0.736563E-06/
      data dst_slt_flx_rat_ttl_lutx ( 96,4) /   0.730991E-06/
      data dst_slt_flx_rat_ttl_lutx ( 97,4) /   0.726790E-06/
      data dst_slt_flx_rat_ttl_lutx ( 98,4) /   0.722917E-06/
      data dst_slt_flx_rat_ttl_lutx ( 99,4) /   0.718549E-06/
      data dst_slt_flx_rat_ttl_lutx (100,4) /   0.714528E-06/
        
 !Horizontal flux for different soils
    
      data  flx_mss_hrz_slt_ttl_lutx(  1,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  2,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  3,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  4,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  5,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  6,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  7,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  8,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  9,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 10,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 11,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 12,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 13,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 14,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 15,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 16,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 17,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 18,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 19,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 20,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 21,1) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 22,1) /   0.160285E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 23,1) /   0.514955E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 24,1) /   0.998590E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 25,1) /   0.159247E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 26,1) /   0.228625E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 27,1) /   0.307329E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 28,1) /   0.394910E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 29,1) /   0.491064E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 30,1) /   0.595604E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 31,1) /   0.708420E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 32,1) /   0.829475E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 33,1) /   0.958783E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 34,1) /   0.109639E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 35,1) /   0.124240E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 36,1) /   0.139688E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 37,1) /   0.155991E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 38,1) /   0.173151E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 39,1) /   0.191188E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 40,1) /   0.210119E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 41,1) /   0.229963E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 42,1) /   0.250737E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 43,1) /   0.272462E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 44,1) /   0.295154E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 45,1) /   0.318833E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 46,1) /   0.343518E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 47,1) /   0.369227E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 48,1) /   0.395980E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 49,1) /   0.423795E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 50,1) /   0.452694E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 51,1) /   0.482689E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 52,1) /   0.513807E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 53,1) /   0.546060E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 54,1) /   0.579472E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 55,1) /   0.614060E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 56,1) /   0.649841E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 57,1) /   0.686838E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 58,1) /   0.725066E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 59,1) /   0.764551E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 60,1) /   0.805305E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 61,1) /   0.847349E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 62,1) /   0.890703E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 63,1) /   0.935384E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 64,1) /   0.981415E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 65,1) /   0.102881E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 66,1) /   0.107759E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 67,1) /   0.112778E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 68,1) /   0.117939E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 69,1) /   0.123244E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 70,1) /   0.128695E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 71,1) /   0.134294E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 72,1) /   0.140044E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 73,1) /   0.145945E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 74,1) /   0.152000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 75,1) /   0.158211E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 76,1) /   0.164579E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 77,1) /   0.171108E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 78,1) /   0.177796E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 79,1) /   0.184649E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 80,1) /   0.191667E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 81,1) /   0.198851E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 82,1) /   0.206205E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 83,1) /   0.213729E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 84,1) /   0.221427E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 85,1) /   0.229299E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 86,1) /   0.237348E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 87,1) /   0.245575E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 88,1) /   0.253984E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 89,1) /   0.262574E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 90,1) /   0.271348E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 91,1) /   0.280310E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 92,1) /   0.289459E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 93,1) /   0.298797E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 94,1) /   0.308328E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 95,1) /   0.318052E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 96,1) /   0.327972E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 97,1) /   0.338089E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 98,1) /   0.348406E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 99,1) /   0.358925E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(100,1) /   0.369646E+00  /
    
      data  flx_mss_hrz_slt_ttl_lutx(  1,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  2,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  3,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  4,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  5,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  6,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  7,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  8,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  9,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 10,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 11,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 12,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 13,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 14,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 15,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 16,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 17,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 18,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 19,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 20,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 21,2) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 22,2) /   0.558435E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 23,2) /   0.211781E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 24,2) /   0.471580E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 25,2) /   0.844190E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 26,2) /   0.133506E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 27,2) /   0.194595E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 28,2) /   0.267567E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 29,2) /   0.352108E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 30,2) /   0.447793E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 31,2) /   0.554170E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 32,2) /   0.670813E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 33,2) /   0.797359E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 34,2) /   0.933515E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 35,2) /   0.107908E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 36,2) /   0.123381E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 37,2) /   0.139753E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 38,2) /   0.157033E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 39,2) /   0.175235E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 40,2) /   0.194371E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 41,2) /   0.214458E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 42,2) /   0.235508E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 43,2) /   0.257540E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 44,2) /   0.280569E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 45,2) /   0.304611E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 46,2) /   0.329687E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 47,2) /   0.355811E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 48,2) /   0.383002E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 49,2) /   0.411280E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 50,2) /   0.440661E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 51,2) /   0.471166E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 52,2) /   0.502810E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 53,2) /   0.535616E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 54,2) /   0.569600E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 55,2) /   0.604783E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 56,2) /   0.641180E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 57,2) /   0.678813E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 58,2) /   0.717702E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 59,2) /   0.757864E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 60,2) /   0.799317E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 61,2) /   0.842080E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 62,2) /   0.886172E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 63,2) /   0.931618E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 64,2) /   0.978432E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 65,2) /   0.102663E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 66,2) /   0.107624E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 67,2) /   0.112728E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 68,2) /   0.117976E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 69,2) /   0.123370E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 70,2) /   0.128913E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 71,2) /   0.134606E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 72,2) /   0.140452E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 73,2) /   0.146452E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 74,2) /   0.152607E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 75,2) /   0.158920E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 76,2) /   0.165394E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 77,2) /   0.172029E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 78,2) /   0.178828E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 79,2) /   0.185793E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 80,2) /   0.192924E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 81,2) /   0.200226E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 82,2) /   0.207699E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 83,2) /   0.215345E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 84,2) /   0.223166E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 85,2) /   0.231164E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 86,2) /   0.239342E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 87,2) /   0.247701E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 88,2) /   0.256242E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 89,2) /   0.264967E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 90,2) /   0.273880E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 91,2) /   0.282982E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 92,2) /   0.292274E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 93,2) /   0.301758E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 94,2) /   0.311436E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 95,2) /   0.321312E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 96,2) /   0.331384E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 97,2) /   0.341657E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 98,2) /   0.352132E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 99,2) /   0.362811E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(100,2) /   0.373696E+00  /
    
      data  flx_mss_hrz_slt_ttl_lutx(  1,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  2,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  3,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  4,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  5,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  6,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  7,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  8,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  9,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 10,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 11,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 12,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 13,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 14,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 15,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 16,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 17,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 18,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 19,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 20,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 21,3) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 22,3) /   0.551503E-06  /
      data  flx_mss_hrz_slt_ttl_lutx( 23,3) /   0.352267E-05  /
      data  flx_mss_hrz_slt_ttl_lutx( 24,3) /   0.121150E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 25,3) /   0.315729E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 26,3) /   0.695055E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 27,3) /   0.136003E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 28,3) /   0.243469E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 29,3) /   0.406169E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 30,3) /   0.639560E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 31,3) /   0.959472E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 32,3) /   0.138128E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 33,3) /   0.191914E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 34,3) /   0.258539E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 35,3) /   0.339016E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 36,3) /   0.432986E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 37,3) /   0.538779E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 38,3) /   0.656756E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 39,3) /   0.787267E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 40,3) /   0.930619E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 41,3) /   0.108706E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 42,3) /   0.125681E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 43,3) /   0.144005E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 44,3) /   0.163691E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 45,3) /   0.184752E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 46,3) /   0.207196E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 47,3) /   0.231034E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 48,3) /   0.256271E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 49,3) /   0.282913E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 50,3) /   0.310964E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 51,3) /   0.340432E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 52,3) /   0.371316E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 53,3) /   0.403625E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 54,3) /   0.437363E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 55,3) /   0.472534E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 56,3) /   0.509141E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 57,3) /   0.547067E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 58,3) /   0.586267E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 59,3) /   0.626762E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 60,3) /   0.668570E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 61,3) /   0.711706E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 62,3) /   0.756196E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 63,3) /   0.802048E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 64,3) /   0.849291E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 65,3) /   0.897937E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 66,3) /   0.948008E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 67,3) /   0.999521E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 68,3) /   0.105249E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 69,3) /   0.110695E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 70,3) /   0.116290E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 71,3) /   0.122036E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 72,3) /   0.127937E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 73,3) /   0.133992E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 74,3) /   0.140205E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 75,3) /   0.146576E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 76,3) /   0.153109E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 77,3) /   0.159805E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 78,3) /   0.166665E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 79,3) /   0.173692E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 80,3) /   0.180886E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 81,3) /   0.188252E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 82,3) /   0.195789E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 83,3) /   0.203500E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 84,3) /   0.211387E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 85,3) /   0.219451E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 86,3) /   0.227695E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 87,3) /   0.236120E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 88,3) /   0.244728E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 89,3) /   0.253522E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 90,3) /   0.262502E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 91,3) /   0.271671E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 92,3) /   0.281031E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 93,3) /   0.290583E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 94,3) /   0.300330E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 95,3) /   0.310272E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 96,3) /   0.320413E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 97,3) /   0.330755E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 98,3) /   0.341296E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 99,3) /   0.352041E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(100,3) /   0.362992E+00  /
    
      data  flx_mss_hrz_slt_ttl_lutx(  1,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  2,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  3,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  4,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  5,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  6,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  7,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  8,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(  9,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 10,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 11,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 12,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 13,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 14,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 15,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 16,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 17,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 18,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 19,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 20,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 21,4) /   0.000000E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 22,4) /   0.647393E-07  /
      data  flx_mss_hrz_slt_ttl_lutx( 23,4) /   0.504983E-06  /
      data  flx_mss_hrz_slt_ttl_lutx( 24,4) /   0.204424E-05  /
      data  flx_mss_hrz_slt_ttl_lutx( 25,4) /   0.613537E-05  /
      data  flx_mss_hrz_slt_ttl_lutx( 26,4) /   0.153201E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 27,4) /   0.336091E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 28,4) /   0.668183E-04  /
      data  flx_mss_hrz_slt_ttl_lutx( 29,4) /   0.122800E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 30,4) /   0.211509E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 31,4) /   0.344885E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 32,4) /   0.536526E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 33,4) /   0.801224E-03  /
      data  flx_mss_hrz_slt_ttl_lutx( 34,4) /   0.115435E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 35,4) /   0.161125E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 36,4) /   0.217554E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 37,4) /   0.283636E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 38,4) /   0.360070E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 39,4) /   0.447575E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 40,4) /   0.546823E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 41,4) /   0.658439E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 42,4) /   0.783002E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 43,4) /   0.921037E-02  /
      data  flx_mss_hrz_slt_ttl_lutx( 44,4) /   0.107303E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 45,4) /   0.123940E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 46,4) /   0.142053E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 47,4) /   0.161678E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 48,4) /   0.182843E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 49,4) /   0.205574E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 50,4) /   0.229894E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 51,4) /   0.255821E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 52,4) /   0.283371E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 53,4) /   0.312558E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 54,4) /   0.343391E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 55,4) /   0.375882E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 56,4) /   0.410025E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 57,4) /   0.445514E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 58,4) /   0.482202E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 59,4) /   0.520111E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 60,4) /   0.559255E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 61,4) /   0.599650E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 62,4) /   0.641317E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 63,4) /   0.684272E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 64,4) /   0.728530E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 65,4) /   0.774111E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 66,4) /   0.821030E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 67,4) /   0.869303E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 68,4) /   0.918953E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 69,4) /   0.969991E-01  /
      data  flx_mss_hrz_slt_ttl_lutx( 70,4) /   0.102244E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 71,4) /   0.107631E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 72,4) /   0.113162E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 73,4) /   0.118839E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 74,4) /   0.124664E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 75,4) /   0.130638E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 76,4) /   0.136763E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 77,4) /   0.143041E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 78,4) /   0.149472E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 79,4) /   0.156061E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 80,4) /   0.162808E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 81,4) /   0.169714E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 82,4) /   0.176781E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 83,4) /   0.184011E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 84,4) /   0.191407E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 85,4) /   0.198968E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 86,4) /   0.206698E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 87,4) /   0.214598E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 88,4) /   0.222670E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 89,4) /   0.230914E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 90,4) /   0.239334E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 91,4) /   0.247931E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 92,4) /   0.256706E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 93,4) /   0.265662E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 94,4) /   0.274799E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 95,4) /   0.284119E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 96,4) /   0.293626E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 97,4) /   0.303319E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 98,4) /   0.313202E+00  /
      data  flx_mss_hrz_slt_ttl_lutx( 99,4) /   0.323275E+00  /
      data  flx_mss_hrz_slt_ttl_lutx(100,4) /   0.333539E+00  /


      !BEGIN CODE
      !PUTTING VALUES FROM ALL ARRAYS DECLARED HERE INTO REAL ARRAYS USED
      !BY ZENDER DUST MODULE
      
      !Mass fraction of source (look up table)
      mss_frc_src_lut(:,:,:)=mss_frc_src_lutx(:,:,:)
      
      !Alpha (mass sandblasting efficiency) (look up table)
      dst_slt_flx_rat_ttl_lut(:,:)=dst_slt_flx_rat_ttl_lutx(:,:)
      
      !Horizontal flux look up table:
      flx_mss_hrz_slt_ttl_lut(:,:)=flx_mss_hrz_slt_ttl_lutx(:,:)
      
      return
  end subroutine dst_slt_sbl_cmn_ini           ! end dst_slt_sbl_cmn_ini()
  
end module dstsltsbl

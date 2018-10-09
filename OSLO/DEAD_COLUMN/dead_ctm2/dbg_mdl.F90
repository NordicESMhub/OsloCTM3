! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dbg_mdl.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $ -*-f90-*- 

! Purpose: Systematic nomenclature for debugging/verbosity levels

! Copyright (C) 1997--2003 Charlie Zender
! License: This file is in the Public Domain

! Usage:
!use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl

module dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  implicit none
 
  integer dbg_lvl ! [enm] Debugging level, initialized in main()
  character(80) prg_nm ! [sng] Program name, initialized in ftn_cmd_ln_sng()

  ! Enumerate bebugging levels
  integer,parameter::dbg_nbr=9 ! [nbr] Number of different debugging levels

  integer,parameter::dbg_off=0 ! [enm] Production mode. Debugging is turned off.
  integer,parameter::dbg_fl=1 ! [enm] Filenames
  integer,parameter::dbg_scl=2 ! [enm] Scalars
  integer,parameter::dbg_crr=3 ! [enm] Current task
  integer,parameter::dbg_sbr=4 ! [enm] Subroutine names on entry and exit
  integer,parameter::dbg_io=5 ! [enm] Subroutine I/O
  integer,parameter::dbg_vec=6 ! [enm] Entire vectors
  integer,parameter::dbg_vrb=7 ! [enm] Everything
  integer,parameter::dbg_old=8 ! [enm] Old debugging blocks not used anymore

end module dbg_mdl

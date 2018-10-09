! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/precision.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $ -*-f90-*-

! Purpose: Portable precision definitions for floating point and integer variables

! Usage:
! use precision ! [mdl] Precision r8, i8, ...

! Compilation
! cd ~/f;f90 -c -I$HOME/include -o $MY_OBJ_DIR/precision.o precision.F90
! cd ~/f;pgf90 -c -I$HOME/include -o $MY_OBJ_DIR/precision.o precision.F90
! cd ~/f;pgf90 -c -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign -Ktrap=fp -fast -DLINUX -I. -I$HOME/include -I/usr/local/include -o $MY_OBJ_DIR/precision.o precision.F90

! Standard precisions may be selected with predefined kind numbers
! Real kinds:
! 4 = Single precision values with approximately  7 significant digits (4 bytes)
! 8 = Double precision values with approximately 14 significant digits (8 bytes)
! IEEE and Cray storage both return 8 byte reals for precision=12 (but not p=15)
! Integer kinds:
! 1 =  8-bit integers [-2^7,2^7-1]
! 2 = 16-bit integers [-2^15,2^15-1]
! 3 = 32-bit integers [-2^31,2^31-1]
! 4 = 64-bit integers [-2^63,2^63-1]

module precision

  ! Definitions are adopted from CCM and MATCH
  integer,parameter::r4=selected_real_kind(p=6) ! r4: 4B (C float) default, 8B (C double) possible
#ifdef PRC_DBL
  integer,parameter::r8=selected_real_kind(p=12) ! r8: 8B (C double) default, 4B (C float) possible
#endif /* !PRC_DBL */
#ifdef PRC_FLT
  integer,parameter::r8=selected_real_kind(p=6) ! r8: 8B (C double) default, 4B (C float) possible
#endif /* !PRC_FLT */

  integer,parameter::rntv=kind(1.0) ! rntv: Native real kind
  integer,parameter::i8=selected_int_kind(13) ! i8: 8B (C long long) default, 4B (C int) possible
  integer,parameter::i4=selected_int_kind(6) ! i4: 4B (C int) default, 8B (C long long) possible
  integer,parameter::intv=selected_int_kind(1) ! intv: Native integer kind
  integer,parameter::DBLKIND=selected_real_kind(p=12) ! DBLKIND: Fixed 8B (double precision) always

end module precision


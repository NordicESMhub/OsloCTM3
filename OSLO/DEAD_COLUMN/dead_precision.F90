

module dead_precision

  !// DEAD precision for floats are taken from CTM3.
  use cmn_precision, only: r8, r4
  implicit none

  ! r4: 4B (C float) default, 8B (C double) possible
  !integer,parameter::r4=selected_real_kind(p=6)

  ! r8: 8B (C double) default, 4B (C float) possible
  !integer,parameter::r8=selected_real_kind(p=12)


  ! rntv: Native real kind
  integer, parameter::rntv=kind(1.0)
  ! i8: 8B (C long long) default, 4B (C int) possible
  integer, parameter::i8=selected_int_kind(13)
  ! i4: 4B (C int) default, 8B (C long long) possible
  integer,parameter::i4=selected_int_kind(6)
  ! intv: Native integer kind
  integer,parameter::intv=selected_int_kind(1)
  ! DBLKIND: Fixed 8B (double precision) always
  integer,parameter::DBLKIND=selected_real_kind(p=12)

end module dead_precision

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


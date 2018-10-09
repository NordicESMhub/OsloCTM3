! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/utl_mdl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Library of general purpose utility routines

! Usage: 
! use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)

module utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  
contains 
  
  logical function mnt_chk(grd,grd_nbr)
    ! Purpose: Return .true. if array is monotonic, .false. if array is not
    ! Prototype:
    ! logical mnt_chk ! [flg] Array is monotonic (libcsz_f90)
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Input
    integer,intent(in)::grd_nbr
    real(r8),intent(in)::grd(grd_nbr)
    ! Local
    integer idx
    logical dcr
    logical ncr
    
    ! Main Code
    ncr=.true.
    if (grd_nbr > 1) then 
       if (grd(2)-grd(1) < 0.0_r8) ncr=.false.
    endif                     ! endif
    dcr=.not.ncr
    
    if (ncr) then
       do idx=2,grd_nbr
          if (grd(idx)-grd(idx-1) < 0.0_r8) goto 100
       end do                 ! end loop over characters
    endif                     ! endif
    if (dcr) then
       do idx=2,grd_nbr
          if (grd(idx)-grd(idx-1) > 0.0_r8) goto 100
       end do                 ! end loop over characters
    endif                     ! endif
    
100 continue
    if (idx == grd_nbr+1) then
       mnt_chk=.true.
    else
       mnt_chk=.false.
       if (dbg_lvl > dbg_io) then
          ! Print non-monotonic points
          write (6,'(a)') 'mnt_chk() reports non-monotonic array:'
          if (ncr) write (6,'(a)') 'Array should be increasing'
          if (dcr) write (6,'(a)') 'Array should be decreasing'
          write (6,'(2(a,i4,a,es10.3))')  &
               'grd(',idx-1,') = ',grd(idx-1), &
               ', grd(',idx,') = ',grd(idx)
       endif                  ! endif dbg
    endif                     ! endif
    
    return
  end function mnt_chk
  
  logical function mnt_ncr_chk(grd,grd_nbr)
    ! Purpose: Return .true. if array monotonically increases, .false. if not
    ! Prototype:
    ! logical mnt_ncr_chk ! [flg] Array is monotonic increasing (libcsz_f90)
    implicit none
    ! Input
    integer,intent(in)::grd_nbr
    real(r8),intent(in)::grd(grd_nbr)
    ! Local
    logical mnt               ! Monotonicity flag
    ! Main Code
    if (grd_nbr <= 1) then
       mnt_ncr_chk=.true.
       return
    endif ! endif
    mnt=mnt_chk(grd,grd_nbr)
    mnt_ncr_chk=.false.
    if (mnt) then 
       if (grd(1) < grd(2)) mnt_ncr_chk=.true.
    endif ! endif
    return
  end function mnt_ncr_chk
  
  logical function mnt_dcr_chk(grd,grd_nbr)
    ! Purpose: Return .true. if array monotonically decreases, .false. if not
    ! Prototype:
    ! logical mnt_ncr_chk ! [flg] Array is monotonic decreasing (libcsz_f90)
    implicit none
    ! Input
    integer,intent(in)::grd_nbr
    real(r8),intent(in)::grd(grd_nbr)
    ! Local
    logical mnt               ! Monotonicity flag
    ! Main Code
    if (grd_nbr <= 1) then
       mnt_dcr_chk=.true.
       return
    endif ! endif
    mnt=mnt_chk(grd,grd_nbr)
    mnt_dcr_chk=.false.
    if (mnt) then 
       if (grd(1) > grd(2)) mnt_dcr_chk=.true.
    endif                     ! endif
    return
  end function mnt_dcr_chk
  
  subroutine date_time_get( & ! [sbr] Format time as Day Mth DD HH:MM:SS TZ YYYY
       lcl_date_time) ! O [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
    ! Purpose: Return time formatted as Day Mth DD HH:MM:SS TZ YYYY
    ! Usage call lcl_date_time_get(lcl_date_time)
    implicit none
    ! Parameters
    ! Input
    ! Output 
    character(26),intent(out)::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
    ! Local workspace
    character(8) date
    character(10) time
    character(5) zone
    integer time_val(8)
    ! Main code
    call date_and_time(date,time,zone,time_val) ! Fortran90 DATE_AND_TIME intrinsic
    write (lcl_date_time,'(i4.4,5(a1,i2.2),a1,i3.3)') &
         time_val(1),'/',time_val(2),'/',time_val(3),' ', &
         time_val(5),':',time_val(6),':',time_val(7),'.',time_val(8)
    return
  end subroutine date_time_get
  
end module utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)

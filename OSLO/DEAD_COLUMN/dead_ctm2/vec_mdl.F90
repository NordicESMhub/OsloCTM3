! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/vec_mdl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Vector manipulation routines

! fxm: vec routines may segfault if compiled with bounds checking/electric fence
! Reason appears to be that ntp_vec() and rbn_vec() _appear_ to access illegal memory
! _appear_ because I checked these routines so many times, bugs are hard to imagine
! Three such potential segfaults are marked in the code

! ntp_vec() and rbn_vec() contain some while() conditions that depend on previous conditions being true
! Compiler optimization level determines whether subsequent conditions are tested when one is false
! If previous condition is false, optimized code never checks subsequent conditions, so everything works
! If previous condition is false, unoptimized code checks subsequent conditions which may be illegal

! Usage: 
! use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning

module vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  private::trp_area_vec
  
contains
  
  subroutine ntp_vec(in_nbr,crd_in,dat_in,out_nbr,crd_out,dat_out,xtr_typ_LHS,xtr_typ_RHS)
    ! Purpose: Project a vector onto another vector
    ! Input and output coordinate vectors must be monotonic
    
    ! crd_in defines the coordinate grid associated with the values dat_in
    ! crd_out defines a new (user) coordinate grid onto which dat_in is linearly interpolated
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use xtr_mdl ! [mdl] Extrapolation constants
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Commons
    ! Input
    integer,intent(in)::in_nbr
    integer,intent(in)::out_nbr
    integer,intent(in)::xtr_typ_LHS
    integer,intent(in)::xtr_typ_RHS
    real(r8),intent(in)::crd_in(in_nbr)
    real(r8),intent(in)::crd_out(out_nbr)
    real(r8),intent(in)::dat_in(in_nbr)
    ! Output
    real(r8),intent(out)::dat_out(out_nbr) ! 
    ! Local
    integer brk_lft_idx
    integer brk_rgt_idx
    integer out_idx
    logical in_dcr
    logical in_ncr
    logical out_dcr
    logical out_ncr
    real(r8) crd_in_mnt(in_nbr) ! [frc] Input coordinate monotonically increasing
    real(r8) crd_out_mnt(out_nbr) ! [frc] Output coordinate monotonically increasing
    real(r8) dat_in_mnt(in_nbr) ! [frc] Input data on monotonically increasing coordinate
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering ntp_vec()'
    
    ! Set extrapolation flags and strings
    call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
    
    ! Vet input
    if (in_nbr <= 1) stop 'nbr_dat_in <= 1 in ntp_vec()'
    if (out_nbr < 1) stop 'nbr_dat_out < 1 in ntp_vec()'
    
    ! Find direction of monotonicity
    if (crd_in(2)-crd_in(1) > 0.0_r8) then 
       in_ncr=.true. 
    else 
       in_ncr=.false.
    endif                     ! endif
    in_dcr=.not.in_ncr
    out_ncr=.true. 
    if (out_nbr > 1) then 
       if (crd_out(2)-crd_out(1) < 0.0_r8) out_ncr=.false.
    endif                     ! endif
    out_dcr=.not.out_ncr
    
    ! Convert input and output arrays to monotonic increasing arrays
    crd_in_mnt(:)=crd_in(:) ! [frc] Input coordinate monotonically increasing
    crd_out_mnt(:)=crd_out(:) ! [frc] Output coordinate monotonically increasing
    dat_in_mnt(:)=dat_in(:) ! [frc] Input data on monotonically increasing coordinate
    if (in_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing input grid in ntp_vec()'
       call rvr_vec(crd_in_mnt,in_nbr)
       call rvr_vec(dat_in_mnt,in_nbr)
    endif                     ! endif in_dcr
    if (out_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing output grid in ntp_vec()'
       call rvr_vec(crd_out_mnt,out_nbr)
    endif                     ! endif in_dcr
    
    ! Initialize bracketing index
    brk_lft_idx=1
    ! Loop over desired output coordinates
    do out_idx=1,out_nbr
       ! Order of conditions is important since second condition is illegal if brk_lft_idx > in_nbr
       do while((brk_lft_idx <= in_nbr).and.(crd_in_mnt(brk_lft_idx) < crd_out_mnt(out_idx)))
          brk_lft_idx=brk_lft_idx+1
       end do                 ! end while
       brk_lft_idx=brk_lft_idx-1
       ! Handle identity interpolation separately to preserve symmetry in extrapolation code 
       if (brk_lft_idx /= in_nbr) then
          if (crd_in_mnt(brk_lft_idx+1) == crd_out_mnt(out_idx)) then
             dat_out(out_idx)=dat_in_mnt(brk_lft_idx+1)
             if (brk_lft_idx == 0) brk_lft_idx=1 ! Reset brk_lft_idx to 1 so next while loop works
             ! NB: In C or F90, use "continue" or "cycle", respectively, instead of "goto" 
             goto 100         ! Jump to next iteration
          endif               ! endif
       endif                  ! endif
       if (brk_lft_idx == 0) then
          ! LHS Extrapolation required
          ! Degenerate case: crd_out_mnt(out_idx) < crd_in_mnt(1)
          brk_lft_idx=1       ! Reset brk_lft_idx to 1 so next while loop works
          if (xtr_flg_LHS(xtr_vrb_msk)) then
             write (6,'(a,i5,a,i5,a,es10.3,a,i5,a,es10.3,a,i5,a,es10.3)')  &
                  'WARNING: ntp_vec() output value dat_out(', &
                  out_idx,') at coordinate crd_out_mnt(', &
                  out_idx,') = ', &
                  crd_out_mnt(out_idx),' requires LHS extrapolation beyond leftmost valid coordinate at crd_in_mnt(', &
                  brk_lft_idx,') = ', &
                  crd_in_mnt(brk_lft_idx),'. Nearest valid datum is dat_in_mnt(', &
                  brk_lft_idx,') = ', &
                  dat_in_mnt(brk_lft_idx)
          endif               ! endif
          ! Extrapolation options are presented in decreasing order of preference
          if (.not.xtr_flg_LHS(xtr_fll_msk)) then
             stop 'ERROR: Full LHS extrapolation required but not permitted in ntp_vec()'
          else if (xtr_flg_LHS(xtr_fll_nil_msk)) then
             dat_out(out_idx)=0.0_r8
          else if (xtr_flg_LHS(xtr_fll_ngh_msk)) then
             dat_out(out_idx)=dat_in_mnt(1)
          else if (xtr_flg_LHS(xtr_fll_lnr_msk)) then
             dat_out(out_idx)=dat_in_mnt(1)- &
                  (crd_in_mnt(1)-crd_out_mnt(out_idx))* &
                  (dat_in_mnt(2)-dat_in_mnt(1))/(crd_in_mnt(2)-crd_in_mnt(1))
          else
             stop 'Unknown xtr_typ_LHS in ntp_vec()'
          endif               ! endif xtr_typ_LHS
          if (xtr_flg_LHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,es10.3)')  &
               xtr_sng_LHS(1:ftn_strlen(xtr_sng_LHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
       else if (brk_lft_idx < in_nbr) then
          ! Normal case: crd_out_mnt is interpolable
          brk_rgt_idx=brk_lft_idx+1
          ! Linearly interpolate
          dat_out(out_idx)= &
               dat_in_mnt(brk_lft_idx)+ &
               (crd_out_mnt(out_idx)-crd_in_mnt(brk_lft_idx))* &
               (dat_in_mnt(brk_rgt_idx)-dat_in_mnt(brk_lft_idx))/ &
               (crd_in_mnt(brk_rgt_idx)-crd_in_mnt(brk_lft_idx))
       else if (brk_lft_idx == in_nbr) then
          ! RHS Extrapolation required
          ! Degenerate case: brk_lft_idx is last element of crd_in_mnt 
          brk_rgt_idx=brk_lft_idx 
          if (xtr_flg_RHS(xtr_vrb_msk)) then
             write (6,'(a,i5,a,i5,a,es10.3,a,i5,a,es10.3,a,i5,a,es10.3)')  &
                  'WARNING: ntp_vec() output value dat_out(', &
                  out_idx,') at coordinate crd_out_mnt(', &
                  out_idx,') = ', &
                  crd_out_mnt(out_idx),' requires RHS extrapolation beyond rightmost valid coordinate at crd_in_mnt(', &
                  brk_rgt_idx,') = ', &
                  crd_in_mnt(brk_rgt_idx),'. Nearest valid datum is dat_in_mnt(', &
                  brk_rgt_idx,') = ', &
                  dat_in_mnt(brk_rgt_idx)
          endif               ! endif
          ! Extrapolation options are presented in decreasing order of preference
          if (.not.xtr_flg_RHS(xtr_fll_msk)) then
             stop 'ERROR: Full RHS extrapolation required but not permitted in ntp_vec()'
          else if (xtr_flg_RHS(xtr_fll_nil_msk)) then
             dat_out(out_idx)=0.0_r8
          else if (xtr_flg_RHS(xtr_fll_ngh_msk)) then
             dat_out(out_idx)=dat_in_mnt(in_nbr)
          else if (xtr_flg_RHS(xtr_fll_lnr_msk)) then
             dat_out(out_idx)=dat_in_mnt(in_nbr)+ &
                  (crd_out_mnt(out_idx)-crd_in_mnt(in_nbr))* &
                  (dat_in_mnt(in_nbr)-dat_in_mnt(in_nbr-1))/ &
                  (crd_in_mnt(in_nbr)-crd_in_mnt(in_nbr-1))
          else
             stop 'Unknown xtr_typ_RHS in ntp_vec()'
          endif               ! endif xtr_typ_RHS
          if (xtr_flg_RHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,es10.3)')  &
               xtr_sng_RHS(1:ftn_strlen(xtr_sng_RHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
       else
          stop 'ERROR: Unforeseen value of brk_lft_idx in ntp_vec()'
       endif
100    continue               ! Branch from identity interpolation
    end do                     ! end loop over output coordinates
    
    ! Convert input and output arrays to original direction of monotonicity
    ! if (in_dcr) then
    !       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing input grid in ntp_vec()'
    !   call rvr_vec(crd_in_mnt,in_nbr)
    !   call rvr_vec(dat_in_mnt,in_nbr)
    ! endif                     ! endif in_dcr
    if (out_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing output grid in ntp_vec()'
       !       call rvr_vec(crd_out_mnt,out_nbr)
       call rvr_vec(dat_out,out_nbr)
    endif                     ! endif in_dcr
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting ntp_vec()'
    return
  end subroutine ntp_vec
  
  subroutine rbn_vec(in_nbr,grd_in,dat_in,  &
       out_nbr,grd_out,dat_out, &
       xtr_typ_LHS,xtr_typ_RHS)
    ! Purpose: Rebin a vector from one grid onto another grid
    ! Input and output coordinate vectors must both be monotonic 
    ! Since they are each monotonic, they are each specified with a single array containing the interface coordinate values
    
    ! crd_min_in and crd_max_in define the coordinate grid associated with the values dat_in
    ! crd_min_out and crd_max_out define a new (user) coordinate grid onto which dat_in is interpolated
    
    ! The output (interpolated) values are set so that the output data has the same integral properties as the input data
    ! In plain English this means the interpolated values are the average values in some neighborhood of the specified grid point
    ! Use this option to interpolate noisy input data, e.g., solar specta or NO2 absorption cross sections, to more regular grids
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl ! [mdl] Utility functions
    use xtr_mdl ! [mdl] Extrapolation constants
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::in_nbr            ! [nbr] max(in_nbr)=49934 for Kur95 1 cm-1 spectrum
    integer,intent(in)::out_nbr           ! [nbr] max(out_nbr)=8192 for TGC98 H2OH2O data
    integer,intent(in)::xtr_typ_LHS
    integer,intent(in)::xtr_typ_RHS
    real(r8),intent(in)::grd_in(in_nbr+1)
    real(r8),intent(in)::grd_out(out_nbr+1)
    real(r8),intent(in)::dat_in(in_nbr)
    ! Output
    real(r8),intent(out)::dat_out(out_nbr) ! Input data rebinned to output grid
    ! Local
    integer brk_lft_idx(out_nbr+1)
    integer brk_rgt_idx(out_nbr+1)
    integer grd_in_nbr
    integer grd_out_idx
    integer grd_out_nbr
    integer idx
    integer in_idx
    integer out_idx
    integer out_idxp1
    integer trp_idx
    integer trp_srt_idx
    integer trp_nbr
    logical mnt
    logical in_dcr
    logical in_ncr
    logical out_dcr
    logical out_ncr
    real(r8) crd_in(in_nbr)       ! Input coordinates (midpoints of grid)
    real(r8) crd_out(out_nbr)     ! Output coordinates (midpoints of grid)
    real(r8) dat_ntp_in(in_nbr+1) ! Linearly interpolated data on grd_in grid
    real(r8) dat_ntp_out(out_nbr+1) ! Linearly interpolated data on grd_out grid
    real(r8) trp_area(in_nbr)     ! Area of each trapezoid defined by input data
    real(r8) grd_in_mnt(in_nbr+1) ! [frc] Input grid monotonically increasing
    real(r8) grd_out_mnt(out_nbr+1) ! [frc] Output grid monotonically increasing
    real(r8) dat_in_mnt(in_nbr) ! [frc] Input data on monotonically increasing grid
    real(r8) crd_dlt_lft
    real(r8) crd_dlt_rgt
    real(r8) crd_dlt_ttl
    real(r8) ovl_avg              ! Overlap average
    real(r8) ovl_frc              ! Overlap fraction
    real(r8) trp_area_lft
    real(r8) trp_area_rgt
    real(r8) trp_area_ttl
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering rbn_vec()'
    
    ! Vet input
    if (in_nbr < 1) stop 'in_nbr < 1 in rbn_vec()'
    if (out_nbr < 1) stop 'out_nbr < 1 in rbn_vec()'
    
    ! Initialize default values
    grd_in_nbr=in_nbr+1
    grd_out_nbr=out_nbr+1
    
    ! Set extrapolation flags and strings
    call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
    
    ! Check for monotonicity
    mnt=mnt_chk(grd_in,in_nbr+1)
    if (.not.mnt) stop 'grd_in not monotonic in rbn_vec()'
    mnt=mnt_chk(grd_out,out_nbr+1)
    if (.not.mnt) stop 'grd_out not monotonic in rbn_vec()'
    
    ! Find direction of monotonicity
    if (grd_in(2)-grd_in(1) > 0.0_r8) then 
       in_ncr=.true. 
    else 
       in_ncr=.false.
    endif
    in_dcr=.not.in_ncr
    if (grd_out(2)-grd_out(1) > 0.0_r8) then 
       out_ncr=.true. 
    else 
       out_ncr=.false.
    endif
    out_dcr=.not.out_ncr
    
    ! Convert both input and output arrays to monotonic increasing arrays
    grd_in_mnt(:)=grd_in(:) ! [frc] Input grid monotonically increasing
    grd_out_mnt(:)=grd_out(:) ! [frc] Output grid monotonically increasing
    dat_in_mnt(:)=dat_in(:) ! [frc] Input data on monotonically increasing grid
    if (in_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing input grid in rbn_vec()'
       call rvr_vec(grd_in_mnt,grd_in_nbr)
       call rvr_vec(dat_in_mnt,in_nbr)
    endif                     ! endif in_dcr
    if (out_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Reversing output grid in rbn_vec()'
       call rvr_vec(grd_out_mnt,grd_out_nbr)
    endif                     ! endif in_dcr
    
    ! Input and output coordinates both increase monotonically
    if (xtr_typ_LHS == xtr_err.and.xtr_typ_RHS == xtr_err) then
       ! Vet input more stringently
       if (grd_out_mnt(grd_out_nbr) < grd_in_mnt(1)) then
          write (6,'(a,i4,a,es10.3,a,es10.3)')  &
               'ERROR: rbn_vec() has grd_out_mnt(grd_out_nbr=',grd_out_nbr, &
               ') = ',grd_out_mnt(grd_out_nbr), &
               ' < grd_in_mnt(1) = ',grd_in_mnt(1)
          stop
       endif
       if (grd_out_mnt(1) > grd_in_mnt(grd_in_nbr)) then
          write (6,'(a,es10.3,a,i4,a,es10.3)')  &
               'ERROR: rbn_vec() has grd_out_mnt(1) = ',grd_out_mnt(1), &
               ' > grd_in_mnt(grd_in_nbr = ',grd_in_nbr, &
               ') = ',grd_in_mnt(grd_in_nbr)
          stop
       endif
    endif                     ! endif xtr_err
    
    ! Initialize coordinate arrays
    do in_idx=1,in_nbr
       crd_in(in_idx)=0.5*(grd_in_mnt(in_idx)+grd_in_mnt(in_idx+1))
    end do                     ! end loop over input coordinates
    do out_idx=1,out_nbr
       crd_out(out_idx)=0.5*(grd_out_mnt(out_idx)+grd_out_mnt(out_idx+1))
    end do                     ! end loop over output coordinates
    
    ! Assemble LHS bracketing indices for each output gridpoint
    do grd_out_idx=1,grd_out_nbr
       brk_lft_idx(grd_out_idx)=1
       ! Initialize bracketing index with previous value if possible
       if (grd_out_idx > 1) then
          if (brk_lft_idx(grd_out_idx-1) >= 1.and. &
               brk_lft_idx(grd_out_idx-1) <= grd_in_nbr)  &
               brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx-1)
       endif                  ! endif
       ! Order of conditions is important since second condition is illegal if brk_lft_idx > grd_in_nbr
       do while((brk_lft_idx(grd_out_idx) <= grd_in_nbr).and. &
            (grd_in_mnt(brk_lft_idx(grd_out_idx)) <= grd_out_mnt(grd_out_idx)))
          brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx)+1
       end do                 ! end while
       ! Decrement bracketing index since last loop iteration overshoots target
       brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx)-1
       
       ! ! fxm: 20000819 potential segfault here when brk_lft_idx[grd_out_idx] == grd_in_nbr+1
       ! if (grd_in_mnt(brk_lft_idx(grd_out_idx)) /= grd_out_mnt(grd_out_idx)) then
       ! ! Decrement bracketing index unless while loop was broken by an equality
       ! brk_lft_idx(grd_out_idx)=brk_lft_idx(grd_out_idx)-1
       ! endif                  ! endif
    end do                     ! end loop over output coordinates
    
    ! Assemble RHS bracketing indices for each output gridpoint
    do grd_out_idx=grd_out_nbr,1,-1
       brk_rgt_idx(grd_out_idx)=grd_in_nbr
       ! Initialize bracketing index with previous value if possible
       if (grd_out_idx < grd_out_nbr) then
          if (brk_rgt_idx(grd_out_idx+1) <= grd_in_nbr.and. &
               brk_rgt_idx(grd_out_idx+1) >= 1)  &
               brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx+1)
       endif                  ! endif
       ! Order of conditions is important since second condition is illegal if brk_rgt_idx < 1
       do while((brk_rgt_idx(grd_out_idx) >= 1).and. &
            (grd_in_mnt(brk_rgt_idx(grd_out_idx)) >= grd_out_mnt(grd_out_idx)))
          brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx)-1
       end do                 ! end while
       ! Increment bracketing index since last loop iteration overshoots target
       brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx)+1
       
       ! ! fxm: 20000819 potential segfault here when brk_rgt_idx[grd_out_idx] == 0
       ! if (grd_in_mnt(brk_rgt_idx(grd_out_idx)) /= grd_out_mnt(grd_out_idx)) then
       ! ! Increment bracketing index unless while loop was broken by an equality
       ! brk_rgt_idx(grd_out_idx)=brk_rgt_idx(grd_out_idx)+1
       ! endif                  ! endif
       
    end do                     ! end loop over output coordinates
    
    ! Vet bracketing indices
    do grd_out_idx=1,grd_out_nbr
       if (brk_lft_idx(grd_out_idx) > brk_rgt_idx(grd_out_idx)) then
          stop 'brk_lft_idx(grd_out_idx) > brk_rgt_idx(grd_out_idx) in rbn_vec()'
       endif                  ! endif
    end do                     ! end loop over output coordinates
    
    ! Precompute all linearly interpolated values for use below
    ! Input values interpolated to input grid (Recall values are input on coordinates, not grid)
    call ntp_vec(in_nbr,crd_in,dat_in_mnt,grd_in_nbr,grd_in_mnt,dat_ntp_in,xtr_fll_lnr,xtr_fll_lnr)
    ! Input values interpolated to output grid
    call ntp_vec(in_nbr,crd_in,dat_in_mnt,grd_out_nbr,grd_out_mnt,dat_ntp_out,xtr_fll_lnr,xtr_fll_lnr)
    ! Precompute area of all trapezoids defined by input grid for use below
    call trp_area_vec(in_nbr,grd_in_mnt,dat_ntp_in,trp_area)
    
    ! Reset extrapolation flags and strings
    call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
    
    ! Assign an output value to each output region
    do out_idx=1,out_nbr
       out_idxp1=out_idx+1
       ! Traverse truth table by examining all possible values of brk_lft_idx(out_idx)
       if (brk_lft_idx(out_idx) < 0) then
          stop 'ERROR: brk_lft_idx(out_idx) < 0 in rbn_vec()'
       else if (brk_lft_idx(out_idx) == 0) then
          ! LHS extrapolation required
          if (xtr_flg_LHS(xtr_vrb_msk)) then
             write (6,'(a,i5,a,i5,a,es10.3,a,i5,a,es10.3,a,es10.3,a,es10.3,a,es10.3)')  &
                  'WARNING: rbn_vec() output bin value centered at crd_out(', &
                  out_idx,') between interface coordinates grd_out_mnt(', &
                  out_idx,') = ', &
                  grd_out_mnt(out_idx),' and grd_out_mnt(', &
                  out_idxp1,') = ', &
                  grd_out_mnt(out_idxp1), &
                  ' requires LHS extrapolation beyond leftmost valid interface coordinates at grd_in_mnt(1) = ', &
                  grd_in_mnt(1),' and grd_in_mnt(2) = ', &
                  grd_in_mnt(2),'. Nearest valid datum is dat_in_mnt(1) = ', &
                  dat_in_mnt(1)
          endif               ! endif xtr_flg_vrb
          if (brk_lft_idx(out_idxp1) == 0) then
             ! Fully degenerate case: current output bin has no overlap with input bins
             if (xtr_flg_LHS(xtr_vrb_msk)) then
                write (6,'(2(a,es10.3))') &
                     'Full LHS extrpolation required. Distances of bin endpoints from valid data are ', &
                     grd_in_mnt(1)-grd_out_mnt(out_idx),' and ',grd_in_mnt(1)-grd_out_mnt(out_idxp1)
             endif            ! endif xtr_flg_vrb
             ! Extrapolation options are presented in decreasing order of preference
             if (.not.xtr_flg_LHS(xtr_fll_msk)) then
                stop 'ERROR: Full LHS extrapolation required but not permitted in rbn_vec()'
             else if (xtr_flg_LHS(xtr_fll_nil_msk)) then
                dat_out(out_idx)=0.0_r8
             else if (xtr_flg_LHS(xtr_fll_ngh_msk)) then
                dat_out(out_idx)=dat_in_mnt(1)
             else if (xtr_flg_LHS(xtr_fll_lnr_msk)) then
                stop 'xtr_fll_lnr not implemented yet in rbn_vec()'
             else
                stop 'Unknown xtr_typ_LHS in fll branch of rbn_vec()'
             endif            ! endif xtr_typ_LHS
          else if (brk_lft_idx(out_idxp1) > 0) then ! brk_lft_idx(out_idxp1) == 0
             ! Half degenerate case: current output bin has partial overlap with input bins
             ovl_frc=(grd_out_mnt(out_idxp1)-grd_in_mnt(1))/(grd_out_mnt(out_idxp1)-grd_out_mnt(out_idx))
             if (xtr_flg_LHS(xtr_vrb_msk)) then
                write (6,'(a,f9.3,a)') 'Partial LHS extrapolation required. Overlap fraction is ',ovl_frc
             endif            ! endif xtr_flg_vrb
             ! Extrapolation options are presented in decreasing order of preference
             if (.not.xtr_flg_LHS(xtr_prt_msk)) then
                stop 'ERROR: Partial LHS extrapolation required but not permitted in rbn_vec()'
             else if (xtr_flg_LHS(xtr_prt_frc_msk)) then
                ! Compute average of overlap input bins
                crd_dlt_rgt=grd_out_mnt(out_idxp1)-grd_in_mnt(brk_lft_idx(out_idxp1))
                trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out(out_idxp1)+dat_ntp_in(brk_lft_idx(out_idxp1)))
                crd_dlt_ttl=crd_dlt_rgt
                trp_area_ttl=trp_area_rgt
                trp_nbr=brk_lft_idx(out_idxp1)-1
                do trp_idx=1,trp_nbr
                   in_idx=trp_idx
                   trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                   crd_dlt_ttl=crd_dlt_ttl+grd_in_mnt(in_idx+1)-grd_in_mnt(in_idx)
                end do        ! end loop over inner trapezoids
                if (crd_dlt_ttl <= 0.0_r8) then 
                   stop 'ERROR: crd_dlt_ttl <= 0.0 in xtr_LHS branch of rbn_vec()'
                endif         ! endif crd_dlt_ttl <= 0.0
                ovl_avg=trp_area_ttl/crd_dlt_ttl 
                if (xtr_flg_LHS(xtr_prt_wgt_msk)) then
                   dat_out(out_idx)=ovl_frc*ovl_avg
                else
                   dat_out(out_idx)=ovl_avg
                endif         ! endif wgt
             else if (xtr_flg_LHS(xtr_prt_lnr_msk)) then
                stop 'xtr_prt_lnr not implemented yet in rbn_vec()'
             else if (xtr_flg_LHS(xtr_prt_ngh_msk)) then
                dat_out(out_idx)=dat_in_mnt(1)
             else if (xtr_flg_LHS(xtr_prt_nil_msk)) then
                dat_out(out_idx)=0.0_r8
             else
                stop 'Unknown xtr_typ_LHS in prt branch of rbn_vec()'
             endif            ! endif xtr_typ_LHS
          else                ! brk_lft_idx(out_idxp1) == 0
             write (6,'(a)') 'Unforeseen brk_lft_idx in xtr_typ_LHS in rbn_vec()'
          endif               ! brk_lft_idx(out_idxp1) == 0
          if (xtr_flg_LHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,es10.3)')  &
               xtr_sng_LHS(1:ftn_strlen(xtr_sng_LHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
       else if ((brk_lft_idx(out_idx) <= grd_in_nbr).and. &
            (brk_rgt_idx(out_idxp1) <= grd_in_nbr)) then ! brk_lft_idx(out_idx) == 0
          ! Normal case: 99% of execution time for large arrays should be in this block
          ! Current and next output gridpoints are both interpolable
          ! Find out how many input crds lay between current and next output gridpoints
          trp_nbr=brk_lft_idx(out_idxp1)-brk_rgt_idx(out_idx)
          if (trp_nbr < 0) then
             ! Case A: Current and next bracketing indices are the same
             ! This means the output crd resolution is finer than the input crd resolution
             ! Use trapezoidal are current point (this is the same as using just a LHS trapezoid)
             dat_out(out_idx)=0.5*(dat_ntp_out(out_idx)+dat_ntp_out(out_idxp1))
          else if (trp_nbr <= in_nbr) then
             ! Case B: Current and previous bracketing indices differ by one
             ! This means the output crd resolution roughly equals the input crd resolution
             ! This case can be handled by the following code since the number of inner trapezoids
             ! is 0 and the loop will not be executed. 
             ! The input/output area will just be the sum of the LHS and RHS trapezoids
             
             ! Case C: Current and previous bracketing indices differ by many.
             ! This means the output crd resolution is coarser than the input crd resolution.
             ! First, figure out area under the RHS and LHS bookend trapezoids:
             crd_dlt_lft=grd_in_mnt(brk_rgt_idx(out_idx))-grd_out_mnt(out_idx)
             trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out(out_idx)+dat_ntp_in(brk_rgt_idx(out_idx)))
             crd_dlt_rgt=grd_out_mnt(out_idxp1)-grd_in_mnt(brk_lft_idx(out_idxp1))
             trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out(out_idxp1)+dat_ntp_in(brk_lft_idx(out_idxp1)))
             ! NB: We use trapezoids to integrate. But we are unable to use the classical 
             ! Trapezoidal rule because the abscissas are, in general, unevenly spaced.
             ! Therefore we use the trapezoidal rule generalized to unevenly spaced abscissas.
             ! If the input grid equals the output grid, there is one inner trapezoid.
             
             ! Trapezoids are indexed with the Fortran (1-based) convention where
             ! trapezoid 1 is bounded by the abscissae crd_in(1) and crd_in(2) and
             ! trapezoid N-1 is bounded by the abscissae crd_in(N-1) and crd_in(N).
             crd_dlt_ttl=crd_dlt_lft+crd_dlt_rgt
             trp_area_ttl=trp_area_lft+trp_area_rgt
             if (brk_lft_idx(out_idx) /= brk_rgt_idx(out_idx)) then
                trp_srt_idx=brk_lft_idx(out_idx)
             else
                trp_srt_idx=brk_lft_idx(out_idx)-1
             endif            ! endif
             do trp_idx=1,trp_nbr
                in_idx=trp_srt_idx+trp_idx
                trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                crd_dlt_ttl=crd_dlt_ttl+grd_in_mnt(in_idx+1)-grd_in_mnt(in_idx)
             end do           ! end loop over inner trapezoids
             ! This relationship results from requiring output rectangle to be equal in area to 
             ! total of input trapezoids (between current and next gridpoints).
             if (crd_dlt_ttl <= 0.0) then 
                write (6,'(a)') 'ERROR: crd_dlt_ttl <= 0.0 in main branch of rbn_vec()'
             else
                dat_out(out_idx)=trp_area_ttl/crd_dlt_ttl 
             endif            ! endif crd_dlt_ttl <= 0.0
          else if (trp_nbr > in_nbr) then
             write (6,'(2(a,i5,a,i5))')  &
                  'brk_lft_idx(',out_idx, &
                  ') = ',brk_lft_idx(out_idx), &
                  ', brk_rgt_idx(',out_idx, &
                  ') = ',brk_rgt_idx(out_idx)
             write (6,'(2(a,i5,a,i5))')  &
                  'brk_lft_idx(',out_idxp1, &
                  ') = ',brk_lft_idx(out_idxp1), &
                  ', brk_rgt_idx(',out_idxp1, &
                  ') = ',brk_rgt_idx(out_idxp1)
             write (6,'(a,i5)') 'trp_nbr = ',trp_nbr
             stop 'Unforeseen trp_nbr in rbn_vec()' 
          endif               ! endif trp_nbr == 0
       else if (brk_rgt_idx(out_idxp1) > grd_in_nbr) then ! brk_lft_idx(out_idx) == 0
          ! RHS extrapolation required
          if (xtr_flg_RHS(xtr_vrb_msk)) then
             write (6,'(a,i5,a,i5,a,es10.3,a,i5,a,es10.3,a,i5,a,es10.3,a,i5,a,es10.3,a,i5,a,es10.3)')  &
                  'WARNING: rbn_vec() output bin value centered at crd_out(', &
                  out_idx,') between interface coordinates grd_out_mnt(', &
                  out_idx,') = ', &
                  grd_out_mnt(out_idx),' and grd_out_mnt(', &
                  out_idxp1,') = ', &
                  grd_out_mnt(out_idxp1), &
                  ' requires RHS extrapolation beyond rightmost valid interface coordinates at grd_in_mnt(', &
                  grd_in_nbr-1,') = ', &
                  grd_in_mnt(grd_in_nbr-1),' and grd_in_mnt(', &
                  grd_in_nbr,') = ', &
                  grd_in_mnt(grd_in_nbr),'. Nearest valid datum is dat_in_mnt(', &
                  in_nbr,') = ', &
                  dat_in_mnt(in_nbr)
          endif               ! endif xtr_flg_vrb
          if (brk_rgt_idx(out_idx) == grd_in_nbr+1) then
             ! Fully degenerate case: current output bin has no overlap with input bins
             if (xtr_flg_RHS(xtr_vrb_msk)) then
                write (6,'(2(a,es10.3))') &
                     'Full RHS extrpolation required. Distances of bin endpoints from valid data are ', &
                     grd_out_mnt(out_idx)-grd_in_mnt(grd_in_nbr),' and ',grd_out_mnt(out_idxp1)-grd_in_mnt(grd_in_nbr)
             endif            ! endif xtr_flg_vrb
             ! Extrapolation options are presented in decreasing order of preference
             if (.not.xtr_flg_RHS(xtr_fll_msk)) then
                stop 'ERROR: Full RHS extrapolation required but not permitted in rbn_vec()'
             else if (xtr_flg_RHS(xtr_fll_nil_msk)) then
                dat_out(out_idx)=0.0_r8
             else if (xtr_flg_RHS(xtr_fll_ngh_msk)) then
                dat_out(out_idx)=dat_in_mnt(in_nbr)
             else if (xtr_flg_RHS(xtr_fll_lnr_msk)) then
                stop 'xtr_fll_lnr not implemented yet in rbn_vec()'
             else
                stop 'Unknown xtr_typ_RHS in rbn_vec()'
             endif            ! endif xtr_typ_RHS
          else                ! brk_rgt_idx(out_idx) == grd_in_nbr+1
             ! Half degenerate case: current output bin has partial overlap with input bins
             ovl_frc=(grd_in_mnt(grd_in_nbr)-grd_out_mnt(out_idx))/(grd_out_mnt(out_idxp1)-grd_out_mnt(out_idx))
             if (xtr_flg_RHS(xtr_vrb_msk)) then
                write (6,'(a,f9.3,a)') 'Partial RHS extrapolation required. Overlap fraction is ',ovl_frc
             endif            ! endif xtr_flg_vrb
             ! Extrapolation options are presented in decreasing order of preference
             if (.not.xtr_flg_RHS(xtr_prt_msk)) then
                stop 'ERROR: Partial RHS extrapolation required but not permitted in rbn_vec()'
             else if (xtr_flg_RHS(xtr_prt_frc_msk)) then
                ! Compute average of overlap input bins
                crd_dlt_lft=grd_in_mnt(brk_rgt_idx(out_idx))-grd_out_mnt(out_idx)
                trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out(out_idx)+dat_ntp_in(brk_rgt_idx(out_idx)))
                crd_dlt_ttl=crd_dlt_lft
                trp_area_ttl=trp_area_lft
                trp_nbr=grd_in_nbr-brk_rgt_idx(out_idx)
                do trp_idx=1,trp_nbr
                   in_idx=brk_lft_idx(out_idx)+trp_idx
                   trp_area_ttl=trp_area_ttl+trp_area(in_idx)
                   crd_dlt_ttl=crd_dlt_ttl+grd_in_mnt(in_idx+1)-grd_in_mnt(in_idx)
                end do        ! end loop over inner trapezoids
                if (crd_dlt_ttl <= 0.0_r8) then 
                   stop 'ERROR: crd_dlt_ttl <= 0.0 in xtr_RHS branch of rbn_vec()'
                endif         ! endif crd_dlt_ttl <= 0.0
                ovl_avg=trp_area_ttl/crd_dlt_ttl 
                if (xtr_flg_RHS(xtr_prt_wgt_msk)) then
                   dat_out(out_idx)=ovl_frc*ovl_avg
                else
                   dat_out(out_idx)=ovl_avg
                endif         ! endif wgt
             else if (xtr_flg_RHS(xtr_prt_lnr_msk)) then
                stop 'xtr_prt_lnr not implemented yet in rbn_vec()'
             else if (xtr_flg_RHS(xtr_prt_ngh_msk)) then
                dat_out(out_idx)=dat_in_mnt(in_nbr)
             else if (xtr_flg_RHS(xtr_prt_nil_msk)) then
                dat_out(out_idx)=0.0_r8
             else
                stop 'Unknown xtr_typ_RHS in prt branch of rbn_vec()'
             endif            ! endif xtr_typ_RHS
          endif               ! brk_rgt_idx(out_idx) == grd_in_nbr+1
          if (xtr_flg_RHS(xtr_vrb_msk)) write (6,'(a,a,i4,a,es10.3)')  &
               xtr_sng_RHS(1:ftn_strlen(xtr_sng_RHS)),' yields dat_out(',out_idx,') = ',dat_out(out_idx)
       else
          write (6,'(a)') 'ERROR: Incomplete truth table rbn_vec():'
          write (6,'(2(a,i5,a,i5))')  &
               'brk_lft_idx(',out_idx, &
               ') = ',brk_lft_idx(out_idx), &
               ', brk_rgt_idx(',out_idx, &
               ') = ',brk_rgt_idx(out_idx)
          stop 'Unforeseen bracketing values in rbn_vec()' 
       endif                  ! endelse RHS extrapolate
    end do                     ! end loop over output coordinates
    
    if (dbg_lvl >= dbg_io) then
       ! Print input crd
       write (6,'(a)') 'Input crd_in to rbn_vec():'
       write (6,'(a,i4,a,a)') 'xtr_typ_LHS = ',xtr_typ_LHS,' = ',xtr_sng_LHS(1:ftn_strlen(xtr_sng_LHS))
       write (6,'(a,i4,a,a)') 'xtr_typ_RHS = ',xtr_typ_RHS,' = ',xtr_sng_RHS(1:ftn_strlen(xtr_sng_RHS))
       write (6,'(5(a,1x))') 'idx','crd_in','crd_in_min','crd_in_max','dat_in_mnt'
       do idx=1,in_nbr
          write (6,'(i4,1x,4(es10.3,1x))')  &
               idx,crd_in(idx),grd_in_mnt(idx),grd_in_mnt(idx+1),dat_in_mnt(idx)
       end do                 ! end loop over crd_out
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_io) then
       ! Print input grd
       write (6,'(a)') 'Input grd_in_mnt rbn_vec():'
       write (6,'(3(a,1x))') 'idx','grd_in_mnt','dat_ntp_in'
       do idx=1,grd_in_nbr
          write (6,'(i4,1x,2(es10.3,1x))')  &
               idx,grd_in_mnt(idx),dat_ntp_in(idx)
       end do                 ! end loop over crd_out
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_io) then
       ! Print bracketing indices
       write (6,'(/,a)') 'Grid output from rbn_vec():'
       write (6,'(5(a,1x))') 'idx','grd_out_mnt','brk_lft','brk_rgt','dat_ntp_out'
       do idx=1,grd_out_nbr
          write (6,'(i4,1x,1(es10.3,1x,i4,1x,i4,1x),es10.3)')  &
               idx, &
               grd_out_mnt(idx),brk_lft_idx(idx),brk_rgt_idx(idx), &
               dat_ntp_out(idx)
       end do                 ! end loop over crd_out
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_io) then
       ! Print results of rebinning
       write (6,'(/,a)') 'Coordinate output from rbn_vec():'
       write (6,'(5(a,1x))') 'idx','crd_out','crd_out_min','crd_out_max','dat_out'
       do idx=1,out_nbr
          write (6,'(i4,1x,4(es10.3,1x))')  &
               idx,crd_out(idx),grd_out_mnt(idx),grd_out_mnt(idx+1),dat_out(idx)
       end do                 ! end loop over crd_out
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_io) then
       ! Print trapezoid info
       write (6,'(/,a)') 'Trapezoid output from rbn_vec():'
       write (6,'(5(a,1x))') 'idx','grd_out_mnt','grd_out_mntp1','trp_nbr','trp_idx'
       do out_idx=1,out_nbr
          trp_nbr=brk_lft_idx(out_idx+1)-brk_rgt_idx(out_idx)
          write (6,'(i4,1x,2(es10.3,1x),i2,a2)',advance="no")  &
               out_idx,grd_out_mnt(out_idx),grd_out_mnt(out_idx+1),trp_nbr,': '
          if (brk_lft_idx(out_idx) /= brk_rgt_idx(out_idx)) then
             trp_srt_idx=brk_lft_idx(out_idx)
          else
             trp_srt_idx=brk_lft_idx(out_idx)-1
          endif               ! endif
          do trp_idx=1,trp_nbr-1
             in_idx=trp_srt_idx+trp_idx
             write (6,'(i2,a1)',advance="no") in_idx,','
          end do              ! end loop over trp
          in_idx=trp_srt_idx+trp_nbr
          if (trp_nbr > 0) then 
             write (6,'(i2)') trp_srt_idx+trp_nbr 
          else 
             write (6,'(a4)') 'none'
          endif               ! endif
       end do                 ! end loop over crd_out
    endif                     ! endif dbg
    
    ! Convert input and output arrays to original direction of monotonicity
    if (in_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing input grid in rbn_vec()'
       call rvr_vec(grd_in_mnt,grd_in_nbr)
       call rvr_vec(dat_in_mnt,in_nbr)
    endif                     ! endif in_dcr
    if (out_dcr) then
       if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Un-reversing output grid in rbn_vec()'
       call rvr_vec(grd_out_mnt,grd_out_nbr)
       call rvr_vec(dat_out,out_nbr)
    endif                     ! endif in_dcr
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting rbn_vec()'
    return
  end subroutine rbn_vec
  
  subroutine rvr_vec(dat_out,out_nbr)
    ! Purpose: Reverse array in place
    ! Order is reversed, but values are unchanged
    implicit none
    ! Parameters
    ! Input
    integer,intent(in)::out_nbr ! [nbr] max(out_nbr)=205001 for NO2.F max(bnd_nbr_JPL,bnd_nbr_NCAR,bnd_nbr_NOAA)
    ! Input/Output
    real(r8),intent(inout)::dat_out(out_nbr) ! Input grid
    ! Local
    integer idx
    real(r8) dat_swp(out_nbr)     ! Swap array
    ! Main Code
    ! Copy input array into local swap array
    do idx=1,out_nbr
       dat_swp(idx)=dat_out(idx)
    end do                    ! end loop over out_idx
    ! Reverse local swap array into output array
    do idx=1,out_nbr
       dat_out(idx)=dat_swp(out_nbr-idx+1)
    end do                    ! end loop over out_idx
    return
  end subroutine rvr_vec
  
  subroutine trp_area_vec(in_nbr,grd_in,dat_in,trp_area)
    ! Purpose: Compute area of each trapezoid in input array of trapezoids
    ! N trapezoids are centered within N+1 gridpoints
    ! Trapezoids are indexed with the Fortran (1-based) convention where
    ! Trapezoid 1 is bounded by the abscissae grd_in(1) and grd_in(2)
    ! Trapezoid 2 is bounded by the abscissae grd_in(2) and grd_in(3)
    ! Trapezoid k is bounded by the abscissae grd_in(k) and grd_in(k+1)
    ! Trapezoid N-1 is bounded by the abscissae grd_in(N-1) and grd_in(N)
    ! Trapezoid N is bounded by the abscissae grd_in(N) and grd_in(N+1)
    implicit none
    ! Input
    integer in_nbr
    real(r8),intent(in)::grd_in(in_nbr+1)     ! Input grid
    real(r8),intent(in)::dat_in(in_nbr+1)     ! Input data
    ! Output
    real(r8),intent(out)::trp_area(in_nbr)     ! Area of each trapezoid defined by input data
    ! Local
    integer in_idx
    ! Main Code
    do in_idx=1,in_nbr
       trp_area(in_idx)=0.5* &
            (grd_in(in_idx+1)-grd_in(in_idx))* &
            (dat_in(in_idx+1)+dat_in(in_idx))
    end do                    ! end loop over trapezoids
    return
  end subroutine trp_area_vec                       ! end trp_area_vec()
  
  real(r8) function ntp_vec_one(crd_nbr,crd_in,dat_in,crd_out_scl)
    ! Purpose: Set up and perform a one point interpolation call to ntp_vec()
    ! This function is fairly useful, e.g., at inverting functions
    ! Function takes care of messy overhead needed to call ntp_vec()
    ! Interpolation technique is set to linear
    ! Prototype:
    ! real(r8) ntp_vec_one       ! Scalar interpolation (libcsz_f90)
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use xtr_mdl ! [mdl] Extrapolation constants
    implicit none
    ! Commons
    ! Input
    integer,intent(in)::crd_nbr
    real(r8),intent(in)::crd_out_scl
    real(r8),intent(in)::crd_in(crd_nbr)
    real(r8),intent(in)::dat_in(crd_nbr)
    ! Output
    ! Local
    real(r8) crd_out(1)
    integer out_nbr
    integer xtr_typ_LHS
    integer xtr_typ_RHS
    real(r8)::dat_out(1)
    ! Main Code
    ! Initialize
    out_nbr=1
    crd_out(1)=crd_out_scl
    ! No extrapolation allowed under any circumstances
    xtr_typ_LHS=xtr_err
    xtr_typ_RHS=xtr_err
    call ntp_vec(crd_nbr,crd_in,dat_in,out_nbr,crd_out,dat_out,xtr_typ_LHS,xtr_typ_RHS)
    ntp_vec_one=dat_out(1)
    return
  end function ntp_vec_one                       ! end ntp_vec_one()
  
  integer function vec_val2idx(crd,crd_nbr,val)
    ! Purpose: Locate index of array member closest to specified value
    ! Prototype:
    ! integer vec_val2idx       ! Locate index of array member closest to specified value (libcsz_f90)
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Commons
    ! Input
    integer,intent(in)::crd_nbr           ! Size of array
    real(r8),intent(in)::crd(crd_nbr)         ! Array
    real(r8),intent(in)::val                  ! Value to search for
    ! Local
    integer idx               ! [idx] Counting index
    integer idx_val           ! Index ov value
    real(r8) dst_new              ! Distance from current point to value
    real(r8) dst_old              ! Closest distance yet to value
    ! Main code
    ! Vet input:
    if (crd_nbr < 1) write (6,'(2a,i4)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR vec_val2idx() reports crd_nbr = ',crd_nbr
    idx_val=1
    dst_old=abs(crd(1)-val)
    do idx=1,crd_nbr
       dst_new=abs(crd(idx)-val)
       if (dst_new < dst_old) then
          idx_val=idx
          dst_old=dst_new
       endif                  ! endif
    end do                     ! end loop over crd
    vec_val2idx=idx_val
    return
  end function vec_val2idx
  
  subroutine vec_set(vec,vec_nbr,val)
    ! Purpose: Set a vector vec of vec_nbr real values to the value val
    ! NB: Routine has same interface as CCM:control/resetr()
    implicit none
    ! Input
    integer,intent(in)::vec_nbr           ! Size of input vector
    real(r8),intent(in)::val                  ! Value to fill vector with
    ! Input/Output
    real(r8),intent(out)::vec(vec_nbr)         ! Input vector
    ! Local
    integer idx               ! [idx] Counting index
    do idx=1,vec_nbr
       vec(idx)=val
    end do                    ! end loop over vec
    return
  end subroutine vec_set
  
#if (defined BXM) && (defined SGI)
  ! Routines which are in default Cray and/or CCM libraries 
  subroutine whenflt(chk_nbr,crd,idx_srd,val,vld_idx,vld_nbr)
    ! Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) < val
    ! Only every idx_srdth value in input array will be checked
    ! Use idx_srd=1 to check every element
    ! chk_nbr is the number of elements which are compared
    ! chk_nbr does not necessarily equal the size of crd or vld_idx
    ! When idx_srd != 1, only a hyperslab of crd is checked
    ! Function is rewrite of CCM:srchutil/whenflt()
    ! This routine is very flexible, but somewhat confusing
    implicit none
    ! Input
    integer,intent(in)::idx_srd           ! Stride to employ
    integer,intent(in)::chk_nbr           ! [nbr] Number of elements to check
    real(r8),dimension(:),intent(in)::crd ! Values of candidates
    real(r8),intent(in)::val              ! Target value to use in comparison
    ! Output
    integer,intent(out)::vld_nbr           ! [nbr] Number of valid points
    integer,dimension(:),intent(out)::vld_idx ! Array of valid indices
    ! Local
    integer idx               ! [idx] Counting index
    integer idx_crr           ! Current index into input array
    ! Main Code
    ! Initialize counter and initial index
    vld_nbr=0
    idx_crr=1
    if (idx_srd < 0) idx_crr=(-idx_srd)*(chk_nbr-1)+1
    do idx=1,chk_nbr
       if (crd(idx_crr) < val) then
          vld_nbr=vld_nbr+1
          vld_idx(vld_nbr)=idx
       end if                 ! endif criteria satisfied
       idx_crr=idx_crr+idx_srd
    end do                    ! end loop over idx
    return
  end subroutine whenflt
  
  subroutine whenfgt(chk_nbr,crd,idx_srd,val,vld_idx,vld_nbr)
    ! Purpose: Return array of indices vld_idx into crd which satisfy crd(vld_idx) > val
    ! Only every idx_srdth value in the input array will be checked
    ! Use idx_srd=1 to check every element
    ! chk_nbr is the number of elements which are compared
    ! chk_nbr does not necessarily equal the size of crd or vld_idx
    ! When idx_srd != 1, only a hyperslab of crd is checked
    ! Function is rewrite of CCM:srchutil/whenfgt()
    ! This routine is very flexible, but somewhat confusing
    implicit none
    ! Input
    integer,intent(in)::idx_srd           ! Stride to employ
    integer,intent(in)::chk_nbr           ! [nbr] Number of elements to check
    real(r8),dimension(:),intent(in)::crd ! Values of candidates
    real(r8),intent(in)::val              ! Target value to use in comparison
    ! Output
    integer,intent(out)::vld_nbr           ! [nbr] Number of valid points
    integer,dimension(:),intent(out)::vld_idx ! Array of valid indices
    ! Local
    integer idx               ! [idx] Counting index
    integer idx_crr           ! Current index into input array
    ! Main Code
    ! Initialize counter and initial index
    vld_nbr=0
    idx_crr=1
    if (idx_srd < 0) idx_crr=(-idx_srd)*(chk_nbr-1)+1
    do idx=1,chk_nbr
       if (crd(idx_crr) > val) then
          vld_nbr=vld_nbr+1
          vld_idx(vld_nbr)=idx
       end if                 ! endif criteria satisfied
       idx_crr=idx_crr+idx_srd
    end do                    ! end loop over idx
    return
  end subroutine whenfgt
#endif /* CRAY or CCM */
  
end module vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning

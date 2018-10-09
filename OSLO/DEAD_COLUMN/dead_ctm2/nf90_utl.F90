! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/nf90_utl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Utilities for netCDF Fortran9x interface

! Usage:
! use nf90_utl ! [mdl] netCDF utilities

module nf90_utl ! [mdl] netCDF utilities
  implicit none
  
contains
  
  integer function nf90_wrp(rcd,msg,rcd_opt)
    ! Purpose: Wrap netCDF Fortran90 interface calls
    ! If no error is indicated, return silently
    ! If error is indicated, print corresponding nf90 message, user message, and exit
    ! If rcd_opt option is present, then 
    ! -61 < rcd_opt < 0: Specified error code is not fatal, return rcd to calling program
    ! rcd_opt == 0: Behaves identically to no rcd_opt
    ! rcd_opt > 0: No error is fatal, return rcd to calling program
    use netcdf ! [mdl] netCDF interface
    use dbg_mdl,only:prg_nm ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Input
    character(len=*),intent(in)::msg ! [sng] Diagnostic message
    integer,intent(in)::rcd ! [rcd] Return code
    integer,optional,intent(in)::rcd_opt ! [rcd] Non-fatal return code
    ! Local
    integer rcd_xcp ! [rcd] Non-fatal return code
    ! Main Code
    ! If no error or error is non-fatal, return value of rcd
    nf90_wrp=rcd ! [enm] Return code
    if(present(rcd_opt)) then
       ! Set exception rcd to user-defined rcd
       rcd_xcp=rcd_opt ! [rcd] Non-fatal return code
    else 
       ! Setting exception rcd to nf90_noerr is same as no exception
       rcd_xcp=nf90_noerr ! [rcd] Non-fatal return code
    endif ! endif
    ! Return if no error or if error matches exception rcd
    if (rcd == nf90_noerr .or. rcd == rcd_xcp) return
    ! Houston, we have a problem...
    write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),nf90_strerror(rcd)
    write (6,'(3a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR ',msg(1:ftn_strlen(msg))
    ! Die, monster, die!
    stop
  end function nf90_wrp ! end nf90_wrp()
  
  subroutine nf90_err_exit(rcd,msg)  
    ! Purpose: Handle netCDF return codes
    ! If no error is indicated, return silently
    ! If error is indicated, print corresponding nf90 message, user message, and exit
    use netcdf ! [mdl] netCDF interface
    use dbg_mdl,only:prg_nm ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Input
    integer,intent(in)::rcd ! [rcd] Return code
    character(len=*),intent(in)::msg ! [sng] Diagnostic message
    ! Main Code
    if (rcd == nf90_noerr) return
    write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),nf90_strerror(rcd)
    write (6,'(3a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR ',msg(1:ftn_strlen(msg))
    stop
  end subroutine nf90_err_exit ! end nf90_err_exit()
  
  integer function nf90_wrp_create(fl_nm,mode,nc_id,sbr_nm)  
    ! Purpose: Create netCDF file
    ! Prototype:
    ! integer nf90_wrp_create ! [fnc] Create new netCDF file
    use netcdf ! [mdl] netCDF interface
    use dbg_mdl,only:prg_nm ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none  
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_nm ! I [sng] Filename
    character(len=*),optional,intent(in)::sbr_nm ! [sng] Subroutine name
    integer,intent(in)::mode ! I [enm] netCDF create mode
    integer,intent(out)::nc_id ! I [ptr] netCDF file handle 
    ! Local
    integer::rcd ! [rcd] Return success code 
    ! Main code
    rcd=nf90_noerr ! [rcd] nf90_noerr == 0
    rcd=rcd+nf90_create(fl_nm,mode,nc_id)  
    nf90_wrp_create=rcd  
    if (rcd /= nf90_noerr) then  
       if(present(sbr_nm)) then 
          write (6, '(6a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR creating ', &
               fl_nm(1:ftn_strlen(fl_nm)),' in ',sbr_nm,'()'
       else
          write (6, '(3a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR creating ',fl_nm(1:ftn_strlen(fl_nm))
       endif ! endif
       write (6,'(2a,i3,2a)') prg_nm(1:ftn_strlen(prg_nm)),': rcd = ',rcd,', ',nf90_strerror(rcd)
       stop  
    endif ! endif err  
    return  
  end function nf90_wrp_create ! end nf90_wrp_create()
  
  integer function nf90_wrp_open(fl_nm,mode,nc_id,sbr_nm)  
    ! Purpose: Open netCDF file
    ! Prototype:
    ! integer nf90_wrp_open ! [fnc] Open existing netCDF file
    use netcdf ! [mdl] netCDF interface
    use dbg_mdl,only:prg_nm ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none  
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_nm ! I [sng] Filename
    character(len=*),optional,intent(in)::sbr_nm ! [sng] Subroutine name
    integer,intent(in)::mode ! I [enm] netCDF read/write mode
    integer,intent(out)::nc_id ! I [ptr] netCDF file handle 
    ! Local
    integer::rcd ! [rcd] Return success code 
    ! Main code
    rcd=nf90_noerr ! [rcd] nf90_noerr == 0
    rcd=rcd+nf90_open(fl_nm,mode,nc_id)  
    nf90_wrp_open=rcd  
    if (rcd /= nf90_noerr) then  
       if(present(sbr_nm)) then 
          write (6, '(6a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR opening ', &
               fl_nm(1:ftn_strlen(fl_nm)),' in ',sbr_nm,'()'
       else
          write (6, '(3a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR opening ',fl_nm(1:ftn_strlen(fl_nm))
       endif ! endif
       write (6,'(2a,i3,2a)') prg_nm(1:ftn_strlen(prg_nm)),': rcd = ',rcd,', ',nf90_strerror(rcd)
       stop  
    endif ! endif err  
    return  
  end function nf90_wrp_open ! end nf90_wrp_open()
  
  integer function nf90_wrp_close(nc_id,fl_nm,vrb_sng,sbr_nm)
    ! Purpose: Close netCDF file
    ! Prototype:
    ! integer nf90_wrp_close ! [fnc] Close existing netCDF file
    use netcdf ! [mdl] netCDF interface
    use dbg_mdl,only:prg_nm ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none  
    ! Commons
    ! Input
    character(len=*),optional,intent(in)::fl_nm ! I [sng] Filename
    character(len=*),optional,intent(in)::vrb_sng ! I [sng] Verb string
    character(len=*),optional,intent(in)::sbr_nm ! [sng] Subroutine name
    integer,intent(in)::nc_id ! I [ptr] netCDF file handle 
    ! Local
    integer::rcd ! [rcd] Return success code 
    ! Main code
    rcd=nf90_noerr ! [rcd] nf90_noerr == 0
    rcd=rcd+nf90_close(nc_id) 
    nf90_wrp_close=rcd  
    if (rcd /= nf90_noerr) then
       if(present(fl_nm)) write (6, '(3a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': ERROR closing ',fl_nm(1:ftn_strlen(fl_nm))
       write (6,'(2a,i3,2a)') prg_nm(1:ftn_strlen(prg_nm)),': rcd = ',rcd,', ',nf90_strerror(rcd)
       stop  
    endif ! endif err  
    if (present(vrb_sng) .and. vrb_sng /= '') &
         write (6,'(a,a2,a,1x,a)') prg_nm(1:ftn_strlen(prg_nm)), &
         ': ',vrb_sng,fl_nm(1:ftn_strlen(fl_nm))
    return  
  end function nf90_wrp_close ! end nf90_wrp_close()
  
  integer function nf90_xtype_r8_get()
    ! Purpose: Return netCDF external type corresponding to computational precision
    ! Prototype:
    ! integer::nf90_xtype_r8_get ! [enm] External netCDF type for r8 kind
    use netcdf,only:nf90_double,nf90_float ! [mdl] netCDF interface
    use precision,only:r8,DBLKIND ! [mdl] Precision r8, i8, ...
    implicit none
    if (kind(r8) == kind(DBLKIND)) then
       nf90_xtype_r8_get=nf90_double 
    else
       nf90_xtype_r8_get=nf90_float
    endif ! endif
    return
  end function nf90_xtype_r8_get
  
  ! Begin wrappers
  integer function nf90_wrp_inq_dimid(nc_id,dmn_nm,dmn_id,rcd_opt)
    ! Purpose: Wrapper for nf90_inq_dimid()
    use netcdf ! [mdl] netCDF interface
    implicit none
    ! Input
    character(len=*),intent(in)::dmn_nm ! [sng] Dimension name
    integer,intent(in)::nc_id ! [id] netCDF file ID
    integer,intent(out)::dmn_id ! [id] Dimension ID
    integer,optional,intent(in)::rcd_opt ! [rcd] Non-fatal return code
    ! Local
    integer rcd ! [rcd] Return code
    integer rcd_xcp ! [rcd] Non-fatal return code
    ! Main Code
    rcd=nf90_inq_dimid(nc_id,dmn_nm,dmn_id) ! [enm] Return code
    nf90_wrp_inq_dimid=rcd
    if(present(rcd_opt)) then
       ! Set exception rcd to user-defined rcd
       rcd_xcp=rcd_opt ! [rcd] Non-fatal return code
    else 
       ! Setting exception rcd to nf90_noerr is same as no exception
       rcd_xcp=nf90_noerr ! [rcd] Non-fatal return code
    endif ! endif
    ! If no error or error is non-fatal, return value of rcd
    if (rcd == nf90_noerr .or. rcd == rcd_xcp) return
    write (6,'(a)') nf90_strerror(rcd) ! Houston, we have a problem...
    write (6,'(2a)') 'ERROR nf90_inq_dimid() failed for dimension',dmn_nm
    stop ! Die, monster, die!
  end function nf90_wrp_inq_dimid

  integer function nf90_wrp_inq_varid(nc_id,var_nm,var_id,rcd_opt)
    ! Purpose: Wrapper for nf90_inq_varid()
    use netcdf ! [mdl] netCDF interface
    implicit none
    ! Input
    character(len=*),intent(in)::var_nm ! [sng] Variable name
    integer,intent(in)::nc_id ! [id] netCDF file ID
    integer,intent(out)::var_id ! [id] Variable ID
    integer,optional,intent(in)::rcd_opt ! [rcd] Non-fatal return code
    ! Local
    integer rcd ! [rcd] Return code
    integer rcd_xcp ! [rcd] Non-fatal return code
    ! Main Code
    rcd=nf90_inq_varid(nc_id,var_nm,var_id) ! [enm] Return code
    nf90_wrp_inq_varid=rcd
    if(present(rcd_opt)) then
       ! Set exception rcd to user-defined rcd
       rcd_xcp=rcd_opt ! [rcd] Non-fatal return code
    else 
       ! Setting exception rcd to nf90_noerr is same as no exception
       rcd_xcp=nf90_noerr ! [rcd] Non-fatal return code
    endif ! endif
    ! If no error or error is non-fatal, return value of rcd
    if (rcd == nf90_noerr .or. rcd == rcd_xcp) return
    write (6,'(a)') nf90_strerror(rcd) ! Houston, we have a problem...
    write (6,'(2a)') 'ERROR nf90_inq_varid() failed for variable',var_nm
    stop ! Die, monster, die!
  end function nf90_wrp_inq_varid
  
end module nf90_utl

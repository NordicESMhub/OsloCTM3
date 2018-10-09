!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// netCDF utilities.
!//=========================================================================
module ncutils
  !// ----------------------------------------------------------------------
  !// MODULE: ncutils
  !// DESCRIPTION: Routines for extracting different fields from
  !//              netcdf-files.
  !//
  !// Contains
  !//   subroutine readnc_3d_from4d
  !//   subroutine readnc_2d_from3d
  !//   subroutine readnc_1d
  !//   subroutine readnc_1d_from2d
  !//   subroutine get_netcdf_var_1d
  !//   subroutine get_netcdf_var_2d
  !//   subroutine get_netcdf_var_3d
  !//   subroutine get_netcdf_var_4d
  !//   subroutine get_netcdf_att_char
  !//   subroutine get_netcdf_var_dims
  !//   subroutine get_netcdf_dim
  !//   subroutine get_netcdf_r4var_1d_from_2d
  !//   subroutine get_netcdf_r4var_2d_from_3d
  !//   subroutine get_netcdf_r4var_3d_from_4d
  !//   subroutine handle_err
  !//   subroutine handle_error
  !//
  !// Ole Amund Sovde, March 2011
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'ncutils.f90'
  !//-----------------------------------------------------------------------
  public
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine readnc_3d_from4d(INFILE, DIMNAME1, DIM1, DIMNAME2, DIM2, &
       DIMNAME3, DIM3, DIMNAME4, GET4ENTRY, GETFIELDNAME, R8XYZ)
    !// --------------------------------------------------------------------
    !// Read netCDF file containing 4D fields (DIM1,DIM2,DIM3,DIM4), and
    !// returns 3D-field for DIM4-entry TSTEP. Typically, we have
    !//   DIM1: longitude
    !//   DIM2: latitude
    !//   DIM3: height
    !//   DIM4: time
    !//
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in)   :: infile     !// Name of netCDFfile
    character(len=*), intent(in)   :: dimname1   !// Name of dimension1
    integer,intent(in)          :: DIM1       !// Dimension of dimname1
    character(len=*), intent(in)   :: dimname2   !// Name of dimension2
    integer,intent(in)          :: DIM2       !// Dimension of dimname2
    character(len=*), intent(in)   :: dimname3   !// Name of dimension3
    integer,intent(in)          :: DIM3       !// Dimension of dimname3
    character(len=*), intent(in)   :: dimname4   !// Name of dimension4
                                              !// No need for DIM4
    integer,intent(in)          :: GET4ENTRY  !// Entry to get from dimension4
    character(len=*), intent(in)   :: getfieldname  !// Name of field to get

    !// Output
    real(r8), intent(out)        :: R8XYZ(DIM1,DIM2,DIM3)

    !// Locals
    integer  :: field_id      !// Id for 1d-field
    integer  :: dim1_id       !// Id for dimension1
    integer  :: dim2_id       !// Id for dimension2
    integer  :: dim3_id       !// Id for dimension3
    integer  :: dim4_id       !// Id for dimension4
    integer  :: nsize1        !// Size of dim1
    integer  :: nsize2        !// Size of dim2
    integer  :: nsize3        !// Size of dim3
    integer  :: nsize4        !// Size of dim4

    integer  :: srt_entry(4)  !// Start point 
    integer  :: cnt_entry(4)  !// Count indexes
    integer  :: status        !// status of process (0=OK)
    integer  :: ncid          !// file id 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='readnc_3d_from4d'
    !// --------------------------------------------------------------------

    status=nf90_noerr  !// Status is 0 and should be kept that way !!

    !// Array which tells you where to start picking your 3D field
    srt_entry = (/ 1, 1, 1, GET4ENTRY /)    !Start array
    !// Array which tells you how far to count when picking it
    cnt_entry = (/ DIM1, DIM2, DIM3, 1 /)         !Count array


    !// Open the existing file
    status = nf90_open(infile, nf90_nowrite, ncid)
    if (status /= NF90_NOERR ) call handle_error(status, &
         subr//':'//trim(infile)//':'//TRIM(getfieldname)//':open')

    !// Inquire dimension ids
    !// DIM1
    status = nf90_inq_dimid(ncid,dimname1,dim1_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 1')
    !// DIM2
    status = nf90_inq_dimid(ncid,dimname2,dim2_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 2')
    !// DIM3
    status = nf90_inq_dimid(ncid,dimname3,dim3_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 3')
    !// DIM4
    status = nf90_inq_dimid(ncid,dimname4,dim4_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 4')


    !// Inquire dimensions
    status = nf90_Inquire_Dimension(ncid,dim1_id,len=nsize1)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 1')
    if (nsize1/=DIM1) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM1 = ',nsize1
       write(6,'(a,i7)') 'your array has dimension ',DIM1
       stop
    endif

    status = nf90_Inquire_Dimension(ncid,dim2_id,len=nsize2)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 2')
    if (nsize2/=DIM2) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM2 = ',nsize2
       write(6,'(a,i7)') 'your array has dimension',DIM2
       stop
    endif
 
    status = nf90_Inquire_Dimension(ncid,dim3_id,len=nsize3)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 3')
    if (nsize3/=DIM3) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM3 = ',nsize3
       write(6,'(a,i7)') '    your array has dimension',DIM3
       stop
    endif
 
    status = nf90_Inquire_Dimension(ncid,dim4_id,len=nsize4)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 4')
    !// Check GET4ENTRY against nsize4
    if (GET4ENTRY .gt. nsize4 .or. nsize4 .le.0) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM4 = ',nsize4
       write(6,'(a,i7)') 'you try to get entry ',GET4ENTRY
       stop
    endif


    !// Number of variables
    !status = nf90_Inquire(ncid,nDims,nVars)

    !// Get variable ID from name
    status=nf90_inq_varid(ncid,getfieldname,field_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_varid')

    !// Get the variable you want and put it in the threeDfield array
    status=nf90_get_var(ncid,field_id,R8XYZ, &
         start=srt_entry, count=cnt_entry )
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':get_var')

    write(6,'(a)') f90file//':'//subr//': Got variable '//trim(getfieldname)

    !// Close file
    status=nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':close')

    !// --------------------------------------------------------------------
  end subroutine readnc_3d_from4d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine readnc_2d_from3d(INFILE, DIMNAME1, DIM1, DIMNAME2, DIM2, &
       DIMNAME3, GET3ENTRY, GETFIELDNAME, R8XY)
    !// --------------------------------------------------------------------
    !// Read netCDF file containing 4D fields (DIM1,DIM2,DIM3), and
    !// returns 2D-field for DIM3-entry TSTEP. Typically, we have
    !//   DIM1: longitude
    !//   DIM2: latitude
    !//   DIM3: time
    !//
    !// Amund Sovde, November 2014
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: infile     !// Name of netCDFfile
    character(len=*), intent(in) :: dimname1   !// Name of dimension1
    integer, intent(in)          :: DIM1       !// Dimension of dimname1
    character(len=*), intent(in) :: dimname2   !// Name of dimension2
    integer, intent(in)          :: DIM2       !// Dimension of dimname2
    character(len=*), intent(in) :: dimname3   !// Name of dimension3

    integer, intent(in)          :: GET3ENTRY  !// Entry to get from dimension3
    character(len=*), intent(in) :: getfieldname  !// Name of field to get

    !// Output
    real(r8), intent(out)        :: R8XY(DIM1,DIM2)

    !// Locals
    integer  :: field_id      !// Id for 1d-field
    integer  :: dim1_id       !// Id for dimension1
    integer  :: dim2_id       !// Id for dimension2
    integer  :: dim3_id       !// Id for dimension3
    integer  :: nsize1        !// Size of dim1
    integer  :: nsize2        !// Size of dim2
    integer  :: nsize3        !// Size of dim3

    integer  :: srt_entry(3)  !// Start point 
    integer  :: cnt_entry(3)  !// Count indexes
    integer  :: status        !// status of process (0=OK)
    integer  :: ncid          !// file id 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='readnc_2d_from3d'
    !// --------------------------------------------------------------------

    status=nf90_noerr  !// Status is 0 and should be kept that way !!

    !// Array which tells you where to start picking your 3D field
    srt_entry = (/ 1, 1, GET3ENTRY /)    !Start array
    !// Array which tells you how far to count when picking it
    cnt_entry = (/ DIM1, DIM2, 1 /)      !Count array


    !// Open the existing file
    status=nf90_open(infile, nf90_nowrite, ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':open')

    !// Inquire dimension ids
    !// DIM1
    status = nf90_inq_dimid(ncid,dimname1,dim1_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 1')
    !// DIM2
    status = nf90_inq_dimid(ncid,dimname2,dim2_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 2')
    !// DIM3
    status = nf90_inq_dimid(ncid,dimname3,dim3_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 3')


    !// Inquire dimensions
    status = nf90_Inquire_Dimension(ncid,dim1_id,len=nsize1)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 1')
    if (nsize1/=DIM1) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM1 = ',nsize1
       write(6,'(a,i7)') 'your array has dimension ',DIM1
       stop
    endif

    status = nf90_Inquire_Dimension(ncid,dim2_id,len=nsize2)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 2')
    if (nsize2/=DIM2) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM2 = ',nsize2
       write(6,'(a,i7)') 'your array has dimension',DIM2
       stop
    endif
 
    status = nf90_Inquire_Dimension(ncid,dim3_id,len=nsize3)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 3')
    !// Check GET3ENTRY against nsize3
    if (GET3ENTRY .gt. nsize3 .or. nsize3 .le.0) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM3 = ',nsize3
       write(6,'(a,i7)') 'you try to get entry ',GET3ENTRY
       stop
    endif


    !// Number of variables
    !status = nf90_Inquire(ncid,nDims,nVars)

    !// Get variable ID from name
    status=nf90_inq_varid(ncid,getfieldname,field_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_varid')

    !// Get the variable you want and put it in the threeDfield array
    status=nf90_get_var(ncid,field_id,R8XY, &
         start=srt_entry, count=cnt_entry )
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':get_var')

    write(6,'(a)') f90file//':'//subr//': Got variable '//trim(getfieldname)

    !// Close file
    status=nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':close')

    !// --------------------------------------------------------------------
  end subroutine readnc_2d_from3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine readnc_sub3d_from3d(INFILE, DIMNAME1, DIM1, DIMNAME2, DIM2, &
       DIMNAME3, DIM3, GETFIELDNAME, srt_entry, cnt_entry, R8XY)
    !// --------------------------------------------------------------------
    !// Read netCDF file containing 4D fields (DIM1,DIM2,DIM3), and
    !// returns 2D-field for DIM3-entry TSTEP. Typically, we have
    !//   DIM1: longitude
    !//   DIM2: latitude
    !//   DIM3: time
    !//
    !// Amund Sovde Haslerud, January 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: infile     !// Name of netCDFfile
    character(len=*), intent(in) :: dimname1   !// Name of dimension1
    integer, intent(in)          :: DIM1       !// Dimension of dimname1
    character(len=*), intent(in) :: dimname2   !// Name of dimension2
    integer, intent(in)          :: DIM2       !// Dimension of dimname2
    character(len=*), intent(in) :: dimname3   !// Name of dimension3
    integer, intent(in)          :: DIM3       !// Dimension of dimname3

    character(len=*), intent(in) :: getfieldname  !// Name of field to get
    integer, intent(in)  :: srt_entry(3)  !// Start point 
    integer, intent(in)  :: cnt_entry(3)  !// Count indexes

    !// Output (note the third dimension)
    real(r8), intent(out)        :: R8XY(DIM1,DIM2,CNT_ENTRY(3))

    !// Locals
    integer  :: field_id      !// Id for 1d-field
    integer  :: dim1_id       !// Id for dimension1
    integer  :: dim2_id       !// Id for dimension2
    integer  :: dim3_id       !// Id for dimension3
    integer  :: nsize1        !// Size of dim1
    integer  :: nsize2        !// Size of dim2
    integer  :: nsize3        !// Size of dim3

    integer  :: status        !// status of process (0=OK)
    integer  :: ncid          !// file id 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='readnc_sub3d_from3d'
    !// --------------------------------------------------------------------

    status=nf90_noerr  !// Status is 0 and should be kept that way !!


    !// Open the existing file
    status=nf90_open(infile, nf90_nowrite, ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':open')

    !// Inquire dimension ids
    !// DIM1
    status = nf90_inq_dimid(ncid,dimname1,dim1_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 1')
    !// DIM2
    status = nf90_inq_dimid(ncid,dimname2,dim2_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 2')
    !// DIM3
    status = nf90_inq_dimid(ncid,dimname3,dim3_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_dimid 3')


    !// Inquire dimensions
    status = nf90_Inquire_Dimension(ncid,dim1_id,len=nsize1)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 1')
    if (nsize1/=DIM1) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM1 = ',nsize1
       write(6,'(a,i7)') 'your array has dimension ',DIM1
       stop
    endif

    status = nf90_Inquire_Dimension(ncid,dim2_id,len=nsize2)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 2')
    if (nsize2/=DIM2) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM2 = ',nsize2
       write(6,'(a,i7)') 'your array has dimension',DIM2
       stop
    endif
 
    status = nf90_Inquire_Dimension(ncid,dim3_id,len=nsize3)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(infile)//':'//trim(getfieldname)//':inquire_dimension 3')
    !// Check DIM3 against size of srt_entry and cnt_entry
    !// Remember that cnt_entry(3) is 12 when reading 12 datasets,
    !// so that when checking against DIM3, we have to use
    !// srt_entry(3)+cnt_entry(3)-1:
    if (srt_entry(3) .gt. DIM3 .or. &
         (srt_entry(3)+cnt_entry(3)-1) .gt. DIM3) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(infile)
       write(6,'(a,i7)') '    reports DIM3 = ',DIM3
       write(6,'(a,i7,a3,i7)') 'you try to get entries ',&
            srt_entry(3), ' : ', (srt_entry(3)+cnt_entry(3)-1)
       stop
    endif


    !// Number of variables
    !status = nf90_Inquire(ncid,nDims,nVars)

    !// Get variable ID from name
    status=nf90_inq_varid(ncid,getfieldname,field_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':inq_varid')

    !// Get the variable you want and put it in the threeDfield array
    status=nf90_get_var(ncid,field_id,R8XY, &
         start=srt_entry, count=cnt_entry )
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':get_var')

    write(6,'(a)') f90file//':'//subr//': Got variable '//trim(getfieldname)

    !// Close file
    status=nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(infile)//':'//trim(getfieldname)//':close')

    !// --------------------------------------------------------------------
  end subroutine readnc_sub3d_from3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine readnc_1d(oneDfield, getfieldname, dimname1, DIM1, filename)
    !// --------------------------------------------------------------------
    !// Open a netCDF file and get 1D field from the file.
    !// Reads _only_ 1D fields, i.e., not by extracting 1D from 2D fields
    !// (see readnc_1d_from2d for that).
    !//
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in)       :: getfieldname !// Name of field to get
    character(len=*), intent(in)       :: filename     !// Name of netCDFfile
    character(len=*), intent(in)       :: dimname1     !// Name of dimension1
    integer,intent(in)              :: DIM1         !// Dimension of dimname1
                                                    !//    and getfieldname

    !// Output
    real(r8), intent(out)   :: oneDfield(DIM1)        !// 1D field

    !// Locals
    integer  :: dim_id       !// Id for dimension
    integer  :: field_id     !// Id for 1d-field
    integer  :: nsize        !// Size of 1d-field, read from file

    integer  :: srt_entry(1) !// Start point 
    integer  :: cnt_entry(1) !// Count indexes
    integer  :: status       !// status of process (0=OK)
    integer  :: ncid         !// file id 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='readnc_1d'
    !// --------------------------------------------------------------------

    !// Array which tells you where to start picking your field
    srt_entry = (/ 1 /)    !Start array
    !// Array which tells you how far to count when picking it
    cnt_entry = (/ DIM1 /)   !Count array
  
    status=nf90_noerr  !Status is 0 and should be kept that way !!

    !// Open the existing file
    status=nf90_open(filename, nf90_nowrite, ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':open')

    !// Inquire dimension ID
    status = nf90_inq_dimid(ncid,dimname1,dim_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inq_dimid')

    !// Inquire dimension (i.e. size)
    status = nf90_Inquire_Dimension(ncid,dim_id,len=nsize)
    if (status /= NF90_NOERR) call handle_error(status, &
       subr//':'//trim(filename)//':'//trim(getfieldname)//':inquire_dimension')
    if (nsize /= DIM1) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(filename)
       write(6,'(a,i7)') '    reports DIM1 = ',nsize
       write(6,'(a,i7)') 'your array has dimension ',DIM1
       stop
    endif

    !// Get variable ID for 1D field
    status=nf90_inq_varid(ncid,getfieldname,field_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inq_varid')

    !// Get the variable you want and put it in the oneDfield array
    status=nf90_get_var(ncid,field_id,oneDfield, &
         start=srt_entry, count=cnt_entry )
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':get_var')

    write(6,'(a)') f90file//':'//subr//': Got variable ',trim(getfieldname)

    !Closing file
    status=nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':close')

    !// --------------------------------------------------------------------
  end subroutine readnc_1d
  !// ----------------------------------------------------------------------

  

  !// ----------------------------------------------------------------------
  subroutine readnc_1d_from2d(oneDfield, getfieldname, dimname1, DIM1, &
       dimname2, GETENTRY, filename )
    !// --------------------------------------------------------------------
    !// Open a netCDF file and get 1D field from a 2D field on the file.
    !// Reads a 1D field for specified entry of dimension 2 (DIM2).
    !// DIM2 does not have to be sent in to routine.
    !// In other words, we extract the whole DIM1 for a specified entry
    !// of DIM2. This entry is specified by GETENTRY.
    !//
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in)  :: filename      !// Name of netCDFfile
    character(len=*), intent(in)  :: getfieldname  !// Name of field to get
    character(len=*), intent(in)  :: dimname1      !// Name of dimension1
    integer,intent(in)         :: DIM1          !// Dimension of dimname1
    character(len=*), intent(in)  :: dimname2      !// Name of dimension2
                                                !// No need for DIM2
    integer,intent(in)         :: GETENTRY      !// Entry to get from dimension2

    !// Output
    real(r8), intent(out)  :: oneDfield(DIM1)  !// 1D field

    !// Locals
    integer  :: dim1_id       !// Id for dimension1
    integer  :: dim2_id       !// Id for dimension2
    integer  :: field_id      !// Id for 1d-field
    integer  :: nsize1        !// Size of dim1
    integer  :: nsize2        !// Size of dim2

    integer  :: srt_entry(2)  !// Start point 
    integer  :: cnt_entry(2)  !// Count indexes
    integer  :: status        !// status of process (0=OK)
    integer  :: ncid          !// file id 
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='readnc_1d_from2d'
    !// --------------------------------------------------------------------

    !// Array which tells you where to start picking your field
    srt_entry = (/ 1, GETENTRY /)    !Start array
    !// Array which tells you how far to count when picking it
    cnt_entry = (/ DIM1, 1 /)   !Count array
  
    status=nf90_noerr  !Status is 0 and should be kept that way !!

    !// Open the existing file
    status=nf90_open(filename, nf90_nowrite, ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':open')

    !// Inquire dimension1 ID
    status = nf90_inq_dimid(ncid,dimname1,dim1_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inq_dimid 1')

    !// Inquire dimension1 (i.e. size)
    status = nf90_Inquire_Dimension(ncid,dim1_id,len=nsize1)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inquire_dimension 1')
    if (nsize1 /= DIM1) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(filename)
       write(6,'(a,i7)') '    reports DIM1 = ',nsize1
       write(6,'(a,i7)') 'your array has dimension ',DIM1
       stop
    endif


    !// Inquire dimension2 ID
    status = nf90_inq_dimid(ncid,dimname2,dim2_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inq_dimid 2')

    !// Inquire dimension1 (i.e. size)
    status = nf90_Inquire_Dimension(ncid,dim2_id,len=nsize2)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inquire_dimension 2')
    !// Test against GETENTRY
    if (GETENTRY .gt. nsize2 .or. nsize2 .le.0) then
       write(6,'(a)') f90file//':'//subr//': File '//trim(filename)
       write(6,'(a,i7)')'    reports DIM2 = ',nsize2
       write(6,'(a,i7)')'you try to get entry ',GETENTRY
       stop
    endif


    !// Get variable ID for 2D field
    status=nf90_inq_varid(ncid,getfieldname,field_id)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':inq_varid')

    !// Get the variable you want and put it in the oneDfield array
    status=nf90_get_var(ncid,field_id,oneDfield, &
         start=srt_entry, count=cnt_entry )
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':get_var')

    write(6,'(a)') f90file//':'//subr//': Got variable ',trim(getfieldname)

    !Closing file
    status=nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//trim(getfieldname)//':close')

    !// --------------------------------------------------------------------
  end subroutine readnc_1d_from2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_0d( filename, varname, data )
    !// --------------------------------------------------------------------
    !// Reads 0D variable (varname) from netCDF file.
    !//
    !// Amund Sovde Haslerud, January 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    !// Output
    real(r8), intent(out):: data

    !// Local variables
    integer :: status, ncID, varID, dimLen, xtype, ndims
    integer :: dimIDs(1)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_0d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':inq_varid')

    !// Variable information
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs,xtype=xtype,ndims=ndims )

!    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
!         subr//':'//trim(filename)//':'//TRIM(varname)//':variable')

!    !// Dimension of the variable
!    status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(1), len=dimLen )
!    !IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
!    IF (status /= NF90_EBADDIM) CALL HANDLE_ERROR(status, &
!         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension')

!    if (dimLen .ne. 0) then
    if (ndims .ne. 0) then
       write(6,'(a,i3)') f90file//':'//subr//': Not zero dimension length: ',ndims
       stop
    end if

    !// Read the variable
    status = NF90_GET_VAR( ncid, varID, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':get_var')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_0d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_1d( filename, varname, data )
    !// --------------------------------------------------------------------
    !// Reads 1D variable (varname) from netCDF file.
    !// From Chris D. Holmes, UCI
    !//
    !// Ole Amund Sovde, March 2012
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    !// Output
    real(r8),allocatable,intent(out):: data(:)

    !// Local variables
    integer :: status, ncID, varID, dimLen
    integer :: dimIDs(1)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_1d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':inq_varid')

    !// Variable information
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':variable')

    !// Dimension of the variable
    status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(1), len=dimLen )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 1')

    !// Allocate the variable
    ALLOCATE( data(dimLen), STAT=status )
    IF (status < 0 ) then
       write(6,'(a)') f90file//':'//subr//': data allocation failed'
       stop
    end if

    !// Read the variable
    status = NF90_GET_VAR( ncid, varID, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':get_var')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_1d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_2d( filename, varname, data, nlon,nlat )
    !// --------------------------------------------------------------------
    !// Reads 2D variable (varname) from netCDF file.
    !// From Chris D. Holmes, UCI
    !//
    !// Ole Amund Sovde, September 2012
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    integer, intent(in)           :: nlon,nlat
    !// Output
    real(r8),intent(out)            :: data(nlon,nlat)
          
    !// Local variables
    integer :: status, ncID, varID
    integer :: dimIDs(3), dimLen(3)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_2d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':inq_varid')

    !// Variable information
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':variable')
          
    !// 1st Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(1),len=dimLen(1))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 1')

    !// 2nd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(2),len=dimLen(2))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 2')

    !// Read the variable
    status = NF90_GET_VAR( ncid, varID, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':get_var')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_3d( filename, varname, data, nlon,nlat,ntime )
    !// --------------------------------------------------------------------
    !// Reads 3D variable (varname) from netCDF file.
    !// From Chris D. Holmes, UCI
    !//
    !// Ole Amund Sovde, March 2012
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in)  :: filename, varname
    integer, intent(in)           :: nlon,nlat,ntime
    !// Output
    real(r8), intent(out)         :: data(nlon,nlat,ntime)

    !// Local variables
    integer :: status, ncID, varID
    integer :: dimIDs(3), dimLen(3)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_3d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':inq_varid')

    !// Variable information
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':variable')
          
    !// 1st Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(1),len=dimLen(1))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 1')

    !// 2nd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(2),len=dimLen(2))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 2')

    !// 3rd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(3),len=dimLen(3))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 3')

    !// Read the variable
    status = NF90_GET_VAR( ncid, varID, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':get_var')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine inq_netcdf_var(filename, varname, error)
    !// --------------------------------------------------------------------
    !// Inquires if a variable exist on file. Return error code.
    !//
    !// Amund Sovde Haslerud, August 2017
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: filename, varname
    !// Output
    integer, intent(out)         :: error

    !// Local variables
    integer :: status, ncID, varID
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='inq_nc_var'
    !// --------------------------------------------------------------------

    !// Open file
    status = nf90_open(trim(filename), NF90_NoWrite, ncid)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//': open')

    !// Inquire variable id
    error = nf90_inq_varid(ncid, trim(varname), varID)
    if (.not. (status .eq. NF90_NOERR .or. status .eq. NF90_ENOTVAR)) &
         call handle_error(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//': inq_varid')

    !// Close
    status = nf90_close(ncid)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//': close')

    !// --------------------------------------------------------------------
  end subroutine inq_netcdf_var
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_4d(filename,varname,data, nlon,nlat,nlev,ntime)
    !// --------------------------------------------------------------------
    !// Reads 4D variable (varname) from netCDF file.
    !//
    !// Ole Amund Sovde, April 2016
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    integer, intent(in)           :: nlon,nlat,nlev,ntime
    !// Output
    real(r8),intent(out)            :: data(nlon,nlat,nlev,ntime)

    !// Local variables
    integer :: status, ncID, varID
    integer :: dimIDs(4), dimLen(4)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_4d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':inq_varid')

    !// Variable information
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':variable')
          
    !// 1st Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(1),len=dimLen(1))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 1')

    !// 2nd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(2),len=dimLen(2))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 2')

    !// 3rd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(3),len=dimLen(3))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 3')

    !// 4rd Dimension of the variable
    status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(4),len=dimLen(4))
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':dimension 4')

    !// Read the variable
    status = NF90_GET_VAR( ncid, varID, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':get_var')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, &
         subr//':'//trim(filename)//':'//TRIM(varname)//':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_4d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_att_char( filename, varname, attname, data )
    !// --------------------------------------------------------------------
    !// Reads attribute (attname) data for variable (varname) from netCDF file.
    !// From Chris D. Holmes, UCI
    !//
    !// Ole Amund Sovde, March 2012
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname, attname
    !// Output
    character(len=*),intent(out)  :: data
          
    !// Local variables
    integer :: status, ncID, varID, dimLen
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_att_char'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status,subr // ':open')

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status,subr // ':inq_varid')

    !// Read the attribute
    status = NF90_GET_ATT( ncid, varID, attname, data )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status,subr // ':get_att')

    !// Trim character string
    data = TRIM( data )

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, subr // ':close')

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_att_char
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_var_dims( filename, varname, dims )
    !// --------------------------------------------------------------------
    !// Reads dimensions of variable (varname) from netCDF file.
    !//
    !// Ole Amund Sovde, July 2015
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    !// Output
    integer,allocatable,intent(out):: dims(:)

    !// Local variables
    integer :: status, ncID, varID, I, ND
    integer, parameter :: maxdims = 20 ! should be more than enough
    integer, dimension(maxdims) :: dimIDs, dimlen
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_dims'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr // ':open '//TRIM(filename))

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    if (status /= NF90_NOERR) call handle_error(status,subr // ':inq_varid')

    !// Variable information
    dimIDs(:) = 0
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    if (status /= NF90_NOERR) call handle_error(status,subr // ':variable')

    !// Dimension of the variable
    ND = 0 ! number of dimensions
    dimLen(:) = 0
    do I = 1, maxdims
       if (dimIDs(I) .gt. 0) then
          status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(I), len=dimLen(I) )
          IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status,subr//':dimension')
          ND = ND + 1
       end if
    end do

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, subr // ':close')

    !// Allocate the variable
    ALLOCATE( dims(ND), STAT=status )
    IF (status < 0 ) then
       write(6,'(a)') f90file//':'//subr//': data allocation failed'
       stop
    end if

    !// Return the dimensions
    dims(1:ND) = dimLen(1:ND)

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_var_dims
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_dim( filename, dimname, dimLen)
    !// --------------------------------------------------------------------
    !// Finds the size of a given dimension on netCDF file.
    !//
    !// Amund Sovde Haslerud, November 2017
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, dimname
    !// Output
    integer, intent(out):: dimLen

    !// Local variables
    integer :: status, ncID
    integer  :: dimID
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='get_netcdf_var_dims'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr // ':open '//TRIM(filename))

    !// Locate dimension
    status = NF90_INQ_DIMID( ncid, TRIM(dimname), dimID )
    if (status /= NF90_NOERR) call handle_error(status,subr // ':inq_dimid')

    !// Dimension length
    status = NF90_INQUIRE_DIMENSION( ncid, dimID, len=dimLen )
    if (status /= NF90_NOERR) call handle_error(status,subr // ':inq_dimension')

    !// Close
    status = NF90_CLOSE( ncid )
    IF (status /= NF90_NOERR) CALL HANDLE_ERROR(status, subr // ':close')


    !// --------------------------------------------------------------------
  end subroutine get_netcdf_dim
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_r4var_1d_from_2d(filename, varname, data, &
       n1, n2, nstep)
    !// --------------------------------------------------------------------
    !// Reads netCDF4 file.
    !// Retrieves 1D variable varname(n1) from a 2D field of size (n1,n2).
    !// The n2 step defined by nstep is returned.
    !//
    !// Ole Amund Sovde, July 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    integer, intent(in)           :: n1, n2, nstep
    !// Output
    real(r4),intent(out)          :: data(n1)

    !// Local variables
    integer :: status, ncID, varID, I, ND
    integer, parameter :: maxdims = 20 ! should be more than enough
    integer, dimension(maxdims) :: dimIDs, dimLen
    !// --------------------------------------------------------------------
    character(len=*),parameter :: subr='get_netcdf_var_1d_from_2d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':open: '//trim(filename))

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':inq_varid: '//trim(varname))

    !// Variable information
    dimIDs(:) = 0
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':variable: '//trim(varname))

    !// Find the dimensions n1 and n2
    ND = 0 ! number of dimensions
    dimLen(:) = 0
    do I = 1, maxdims
       if (dimIDs(I) .gt. 0) then
          status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(I), len=dimLen(I) )
          if (status /= NF90_NOERR) &
               call handle_error(status,subr//':dimension: '//trim(varname))
          ND = ND + 1
       end if
    end do

    if (ND .ne. 2) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong dimension on field: '//trim(varname)
       write(6,'(a,i6)') ' Expected ND=2, found: ',ND
       stop
    end if

    !// Check dimensions of field against expected (n1, n2)
    if (dimLen(1) .ne. n1 .or. (dimLen(2) .ne. n2)) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong dimension on field: '//trim(varname)
       write(6,'(a,2i6)') ' Expected/found n1: ',n1,dimLen(1)
       write(6,'(a,2i6)') ' Expected/found n2: ',n2,dimLen(2)
       stop
    end if

    !// Read the variable for nstep
    if (nstep .gt. n2) then
       write(6,'(a,i7)') f90file//':'//subr//': nstep > n2',nstep
       stop
    end if
    status = NF90_GET_VAR( ncid, varID, data, &
         start = (/1, nstep/), count=(/n1, 1 /) )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':get_var: '//trim(varname))

    !// Close
    status = NF90_CLOSE( ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':close: '//trim(filename))

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_r4var_1d_from_2d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_r4var_2d_from_3d(filename, varname, data, &
       n1, n2, n3, nstep)
    !// --------------------------------------------------------------------
    !// Reads netCDF4 file.
    !// Retrieves 2D variable varname(n1,n2) from a 3D field of size
    !// (n1,n2,n3).
    !// The n3 step defined by nstep is returned.
    !//
    !// Ole Amund Sovde, July 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    integer, intent(in)           :: n1, n2, n3, nstep
    !// Output
    real(r4),intent(out)            :: data(n1,n2)

    !// Local variables
    integer :: status, ncID, varID, I, ND
    integer, parameter :: maxdims = 20 ! should be more than enough
    integer, dimension(maxdims) :: dimIDs, dimLen
    !// --------------------------------------------------------------------
    character(len=*),parameter :: subr='get_netcdf_var_2d_from_3d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':open: '//trim(filename))

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':inq_varid: '//trim(varname))

    !// Variable information
    dimIDs(:) = 0
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':variable: '//trim(varname))
          
    !// Find the dimensions
    ND = 0 ! number of dimensions
    dimLen(:) = 0
    do I = 1, maxdims
       if (dimIDs(I) .gt. 0) then
          status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(I), len=dimLen(I) )
          if (status /= NF90_NOERR) &
               call handle_error(status,subr//':dimension: '//trim(varname))
          ND = ND + 1
       end if
    end do

    if (ND .ne. 3) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong dimension on field: '//trim(varname)
       write(6,'(a,i6)') ' Expected ND=3, found: ',ND
       stop
    end if

    !// Check dimensions of field against expected (n1, n2, n3)
    if (dimLen(1) .ne. n1 .or. (dimLen(2) .ne. n2) &
         .or. (dimLen(3) .ne. n3)) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong size on field: '//trim(varname)
       write(6,'(a,2i6)') ' Expected/found n1: ',n1,dimLen(1)
       write(6,'(a,2i6)') ' Expected/found n2: ',n2,dimLen(2)
       write(6,'(a,2i6)') ' Expected/found n3: ',n3,dimLen(3)
       stop
    end if

    !// Read the variable for nstep
    if (nstep .gt. n3) then
       write(6,'(a,i7)') f90file//':'//subr//': nstep > n3',nstep
       stop
    end if
    status = NF90_GET_VAR( ncid, varID, data, &
         start = (/1, 1, nstep/), count=(/n1, n2, 1 /) )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':get_var: '//trim(varname))

    !// Close
    status = NF90_CLOSE( ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':close: '//trim(filename))

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_r4var_2d_from_3d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_netcdf_r4var_3d_from_4d( filename, varname, data, &
       n1, n2, n3, n4, nstep )
    !// --------------------------------------------------------------------
    !// Reads netCDF4 file.
    !// Retrieves 3D variable varname(n1,n2,n3) from a 4D field of size
    !// (n1,n2,n3,n4).
    !// The n4 step defined by nstep is returned.
    !//
    !// Ole Amund Sovde, July 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*),intent(in)   :: filename, varname
    integer, intent(in)           :: n1, n2, n3, n4, nstep
    !// Output
    real(r4),intent(out)            :: data(n1,n2,n3)

    !// Local variables
    integer :: status, ncID, varID, I, ND
    integer, parameter :: maxdims = 20 ! should be more than enough
    integer, dimension(maxdims) :: dimIDs, dimLen
    !// --------------------------------------------------------------------
    character(len=*),parameter :: subr='get_netcdf_var_3d_from_4d'
    !// --------------------------------------------------------------------

    !// Open file
    status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':open: '//trim(filename))

    !// Locate variable
    status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':inq_varid: '//trim(varname))

    !// Variable information
    dimIDs(:) = 0
    status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
    if (status /= NF90_NOERR) &
         call handle_error(status,subr//':variable: '//trim(varname))
          
    !// Find the dimensions
    ND = 0 ! number of dimensions
    dimLen(:) = 0
    do I = 1, maxdims
       if (dimIDs(I) .gt. 0) then
          status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(I), len=dimLen(I) )
          if (status /= NF90_NOERR) &
               call handle_error(status,subr//':dimension: '//trim(varname))
          ND = ND + 1
       end if
    end do

    if (ND .ne. 4) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong dimension on field: '//trim(varname)
       write(6,'(a,i6)') ' Expected ND=4, found: ',ND
       stop
    end if

    !// Check dimensions of field against expected (n1, n2, n3)
    if (dimLen(1) .ne. n1 .or. (dimLen(2) .ne. n2) &
         .or. (dimLen(3) .ne. n3) .or. (dimLen(4) .ne. n4)) then
       write(6,'(a)') f90file//':'//subr// &
            ': Wrong size on field: '//trim(varname)
       write(6,'(a,2i6)') ' Expected/found n1: ',n1,dimLen(1)
       write(6,'(a,2i6)') ' Expected/found n2: ',n2,dimLen(2)
       write(6,'(a,2i6)') ' Expected/found n3: ',n3,dimLen(3)
       write(6,'(a,2i6)') ' Expected/found n4: ',n4,dimLen(4)
       stop
    end if
          
    !// Read the variable for nstep
    if (nstep .gt. n4) then
       write(6,'(a,i7)') f90file//':'//subr//': nstep > n4',nstep
       stop
    end if
    status = NF90_GET_VAR( ncid, varID, data, &
         start = (/1, 1, 1, nstep/), count=(/n1, n2, n3, 1 /) )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':get_var: '//trim(varname))

    !// Close
    status = NF90_CLOSE( ncid )
    if (status /= NF90_NOERR) &
         call handle_error(status, subr//':close: '//trim(filename))

    !// --------------------------------------------------------------------
  end subroutine get_netcdf_r4var_3d_from_4d
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine handle_err(status)
    !// --------------------------------------------------------------------
    !// Handle netCDF errors. This routine does not print out much
    !// helpful information.
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    integer, intent(in)   :: status  !Error status
    !// --------------------------------------------------------------------
    character(len=*),parameter :: subr='handle_err'
    !// --------------------------------------------------------------------
    if (status/=nf90_noerr) then
       write(6,'(a,i7)') f90file//':'//subr// &
            ': ERROR: '//trim(nf90_strerror(status))
       stop "stopped"
    end if
    !// --------------------------------------------------------------------
  end subroutine handle_err
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine handle_error(ECODE,EMSG)
    !// --------------------------------------------------------------------
    !// Handle netCDF errors slightly better than only handle_err.
    !// Ole Amund Sovde, March 2011
    !// --------------------------------------------------------------------
    use netcdf
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    integer, intent(in) :: ECODE
    character(len=*) :: EMSG
    !// --------------------------------------------------------------------
    character(len=*),parameter :: subr='handle_error'
    !// --------------------------------------------------------------------
    if (ECODE /= nf90_noerr) then
       write(6,'(a)') f90file//':'//subr//': netCDF error!'
       write(*,'(a,i5)') '  Error code: ',ECODE
       write(*,'(a)')    '  Error message: '//EMSG
       write(*,'(a)')    '  nf90_strerror: '//trim(nf90_strerror(ECODE))
       stop "STOPPED"
    end if
    !// --------------------------------------------------------------------
  end subroutine handle_error
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module ncutils
!//=========================================================================

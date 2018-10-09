!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, August 2017
!//=========================================================================
!// Initialize STT and model fields.
!//=========================================================================
module stt_save_load
  !// ----------------------------------------------------------------------
  !// MODULE: stt_save_load
  !// DESCRIPTION: Routines for setting STT and relevant fields from
  !//              restart files or monthly averages.
  !//              Also routine for saving restart file.
  !//
  !// Contains
  !//   subroutine load_restart_file
  !//   subroutine getField_and_interpolate
  !//   subroutine save_restart_file
  !// Old versions
  !//   subroutine OSLO_CON_SAV
  !//   subroutine OSLO_CON_RUN
  !//   subroutine restart_from_CTM3avg
  !//   subroutine restart_from_CTM3avg_T42
  !//   subroutine OSLO_CON_RUN42
  !//   subroutine OSLO_CON_RUNxx
  !//   subroutine OSLO_RESTARTFILE_INFO
  !//
  !// Amund Sovde Haslerud, August 2017
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'stt_save_load.f90'
  !// ----------------------------------------------------------------------
  public
  private getField_and_interpolate
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine load_restart_file(FILENAME, MODE)
    !// --------------------------------------------------------------------
    !// Read CTM3 restart file (netCDF4).
    !// UCI CON_RUN modified to include non-transported Oslo components.
    !// Reads component by component, and places the data in the current
    !// array. Works fine when current components are less than in file.
    !//
    !// Can read ADDITIONAL restart files to fill in non-initialized
    !// tracers. This is controlled by MODE:
    !// MODE:  0: First call; initializes zeroinit.
    !//      !=0: Will add non-initialized tracers that are
    !//           available on the specified file. Skips AIR.
    !//
    !// Amund Sovde Haslerud, 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_ctm, only: IYEAR, IDAY, JYEAR, JDAY, JMON, JDATE, NTM, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW
    use cmn_chem, only: TNAME, TMASS
    use utilities, only: get_free_fileid
    use cmn_oslo, only: ZEROINIT, XZEROINIT, trsp_idx, chem_idx, &
          XTNAME, XTMASS, Xtrsp_idx, Xchem_idx, XSTT
    use netcdf
    use ncutils, only: inq_netcdf_var, handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in):: FILENAME
    integer, intent(in) :: MODE

    !// Locals
    real(r8) :: SUMT
    real(r8), dimension(IPAR,JPAR,LPAR) :: R8XYZ

    integer :: &
         version, IN_IPAR, IN_JPAR, IN_LPAR, IN_NDAYS, IN_IYEAR, &
         I, J, L, N, M, ND, XND
    logical :: LREGRID
    character(len=80) :: RTITLE2

    integer :: &
         var_id, &
         lat_dim_id, lon_dim_id, lev_dim_id, &
         lat_id, lon_id, lev_id, &
         ilat_id, ilon_id, ilev_id, &
         ihya_id, ihyb_id, &
         version_id, &
         time_id, ndays_id, iyear_id, &
         status, ncid
    integer, dimension(6) :: in_time
    integer, dimension(3) :: dimIDs3d, dimLen3d
    integer, allocatable, dimension(:) :: Xcomp_id
    character(len=20) :: varname

    !// The moments
    character(len=3), dimension(9) :: &
         moments = (/'SUT','SVT','SWT','SUU','SVV','SWW','SUV','SUW','SVW'/)

    !// Fields read from file
    real(r8), dimension(:), allocatable :: &
         IN_XDEDG, IN_YDEDG, IN_ETAA, IN_ETAB

    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'load_restart_file'
    !// --------------------------------------------------------------------

    !// Initialize arrays; assume all are zero
    if (MODE .eq. 0) then
       ZEROINIT(:)  = 1
       if (NOTRPAR .gt. 0) XZEROINIT(:) = 1
    end if

    !//read in continuation restart file:
    if (MODE .eq. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': reading restart file: '//trim(FILENAME)
    else
       write(6,'(a)') f90file//':'//subr// &
            ': reading ADDITIONAL restart file: '//trim(FILENAME)
    end if

    !//---------------------------------------------------------------------
    status = nf90_open(filename, nf90_NoWrite, ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in opening file '//trim(filename))
    !//---------------------------------------------------------------------
    write(6,'(a)') 'Restart file info:'

    !// Fetch model info
    status = nf90_get_att(ncid, NF90_GLOBAL, 'runtitle', RTITLE2)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get runtitle')
    write(6,'(a)') '  Runtitle: '//trim(RTITLE2)

    !// Model date
    status = nf90_inq_varid(ncid, 'output_time', time_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get output_time varID')
    status = nf90_get_var(ncid, time_id, in_time)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get output_time')
    write(6,'(a,6i6)') '  file time:  ',in_time
    write(6,'(a,6i6)') '  model time: ',JYEAR,JMON,JDATE,0,0,0

    !// Find NDAYS (NDAY+1 - NDAYI) = # integrated days
    status = nf90_inq_varid(ncid, 'NDAYS', ndays_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get NDAYS ID')
    status = nf90_get_var(ncid, ndays_id, IN_NDAYS)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get NDAYS')
    write(6,'(a,i5)') '  Written after # days: ',IN_NDAYS

    !// Find IYEAR
    status = nf90_inq_varid(ncid, 'IYEAR', iyear_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get IYEAR ID')
    status = nf90_get_var(ncid, iyear_id, IN_IYEAR)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get IYEAR')
    write(6,'(a,i5)') '  IYEAR: ',IN_IYEAR


    !// Fetch VERSION
    status = nf90_inq_varid(ncid, 'VERSION', version_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get VERSION ID')
    status = nf90_get_var(ncid, version_id, version)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get VERSION')
    write(6,'(a,i5)') '  VERSION: ',VERSION
    !//---------------------------------------------------------------------

    !// Fetch resolution
    status = nf90_inq_dimid(ncid, 'lon', lon_dim_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get lon_dim_id')
    status = nf90_inquire_dimension(ncid, lon_dim_id, len=IN_IPAR)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get dimension IN_IPAR')

    status = nf90_inq_dimid(ncid, 'lat', lat_dim_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get lat_dim_id')
    status = nf90_inquire_dimension(ncid, lat_dim_id, len=IN_JPAR)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get dimension IN_JPAR')

    status = nf90_inq_dimid(ncid, 'lev', lev_dim_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get lev_dim_id')
    status = nf90_inquire_dimension(ncid, lev_dim_id, len=IN_LPAR)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get dimension IN_LPAR')

    !// Resolution - print-out
    write(6,'(a,3i5)') '  Resolution on file (IPAR,JPAR,LPAR): ', &
         IN_IPAR, IN_JPAR, IN_LPAR
    write(6,'(a,3i5)') '  Model resolution (IPAR,JPAR,LPAR):   ', &
         IPAR, JPAR, LPAR


    !// Check if input data must be regridded
    !//---------------------------------------------------------------------
    if (IN_IPAR .eq. IPAR .and. IN_JPAR .eq. JPAR) then
       LREGRID = .false.
    else
       !// Will only use vertical moments.
       LREGRID = .true.
       !// Print warning about it
       write(6,'(a71)') '--------------------------------------------'// &
            '---------------------------'
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
       write(6,'(a)') '   '//f90file//':'//subr//':'
       write(6,'(a)') '   Restarting from a different horizontal resolution!'
       write(6,'(a)') '   Only STT and vertical moments (SWT and SWW) '// &
            'are interpolated!'
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
       write(6,'(a71)') '--------------------------------------------'// &
            '---------------------------'
    end if

    !//---------------------------------------------------------------------
    !// Allocate read-in array
    allocate( IN_XDEDG(IN_IPAR+1), IN_YDEDG(IN_JPAR+1), &
              IN_ETAA(IN_LPAR+1), IN_ETAB(LPAR+1) )
    !//---------------------------------------------------------------------

    !// IN_XDEDG
    status = nf90_inq_varid(ncid, 'ilon', ilon_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ilon ID')
    status = nf90_get_var(ncid, ilon_id, IN_XDEDG)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get XDEDG')
    !// IN_YDEDG
    status = nf90_inq_varid(ncid, 'ilat', ilat_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ilat ID')
    status = nf90_get_var(ncid, ilat_id, IN_YDEDG)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get YDEDG')
    !// IN_ETAA
    status = nf90_inq_varid(ncid, 'ihya', ihya_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ihya ID')
    status = nf90_get_var(ncid, ihya_id, IN_ETAA)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ETAA')
    !// IN_ETAB
    status = nf90_inq_varid(ncid, 'ihyb', ihyb_id)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ihyb ID')
    status = nf90_get_var(ncid, ihyb_id, IN_ETAB)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)//': get ETAB')


    !// Get variables
    !//---------------------------------------------------------------------

    !// AIR
    if (MODE .eq. 0) then
       !// Check variable
       status = nf90_inq_varid(ncid, 'AIR', var_id)
       if (status .ne. NF90_NOERR) call handle_error(status, &
            f90file//':'//subr//':'//trim(filename)//': get AIR id')
       !// Get field and interpolate
       call getField_and_interpolate(ncid, var_id, filename, &
            'AIR', f90file//':'//subr, &
            IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
            IN_ETAA, IN_ETAB, LREGRID, R8XYZ)

          AIR(:,:,:) = R8XYZ(:,:,:)

       write(6,'(a)') f90file//':'//subr//': got AIR'

    end if !// if (MODE .eq. 0) then


    !// STT
    !//---------------------------------------------------------------------
    do N = 1, NPAR

       if (ZEROINIT(N) .eq. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': already initialised: '//trim(TNAME(N))
          cycle
       end if

       !// Look for STT variable in file STT
       varname = 'STT_'//trim(TNAME(N))

       status = nf90_inq_varid(ncid, TRIM(varname), var_id)

       if (status .eq. NF90_NOERR) then
          !// Variable exists - get it
          call getField_and_interpolate(ncid, var_id, filename, &
               varname, f90file//':'//subr, &
               IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
               IN_ETAA, IN_ETAB, LREGRID, R8XYZ)

          !// Set STT and ZEROINIT
          STT(:,:,:,N) = R8XYZ(:,:,:)
          ZEROINIT(N) = 0

          !// Get moments if they exist
          do m = 1, 9

             !// If regridding, only do vertical moments
             if (LREGRID .and. .not. (moments(m) .eq. 'SWT' .or. &
                                      moments(m) .eq. 'SWW')) cycle

             !// Variable name
             varname = moments(m)//'_'//trim(TNAME(N))

             !// Check variable
             status = nf90_inq_varid(ncid, TRIM(varname), var_id)
             if (status .eq. NF90_NOERR) then

                if (moments(m) .eq. 'SUT') then
                   SUT(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SVT') then
                   SVT(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SWT') then
                   SWT(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SUU') then
                   SUU(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SVV') then
                   SVV(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SWW') then
                   SWW(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SUV') then
                   SUV(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SUW') then
                   SUW(:,:,:,N) = R8XYZ(:,:,:)
                else if (moments(m) .eq. 'SVW') then
                   SVW(:,:,:,N) = R8XYZ(:,:,:)
                end if

             else
                !// This should never happen. If a tracer was in STT,
                !// it should have moments.
                write(6, '(a)') f90file//':'//subr// &
                     'Eroor: Found no moments for species '//&
                     trim(TNAME(N))//' - stopping'
                stop
             end if !// if (status .eq. NF90_NOERR) then

          end do !// do m = 1, 9

          write(6,'(a)') f90file//':'//subr// &
               ': got variable & moments: '//trim(TNAME(N))

       else if (status .eq. NF90_ENOTVAR) then
          !// Variable not found - look in XSTT
          varname = 'XSTT_'//trim(TNAME(N))

          status = nf90_inq_varid(ncid, TRIM(varname), var_id)
          if (status .eq. NF90_NOERR) then
             !// Variable exists - get it
             call getField_and_interpolate(ncid, var_id, filename, &
                  varname, f90file//':'//subr, &
                  IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
                  IN_ETAA, IN_ETAB, LREGRID, R8XYZ)

             !// Set STT and ZEROINIT
             STT(:,:,:,N) = R8XYZ(:,:,:)
             ZEROINIT(N) = 0

             write(6,'(a)') f90file//':'//subr// &
                  ': got variable (XSTT->STT) : '//trim(TNAME(N))

          end if !// if (status .eq. NF90_NOERR) then

       else
          !// Variable does not exist on file - print info
          write(6,'(a)') f90file//':'//subr// &
               ': no such variable: '//trim(TNAME(N))
       end if

       !// If not in file STT or XSTT, component is not initialised.

    end do !// do N = 1, NPAR


    !// XSTT
    do N = 1, NOTRPAR

       if (XZEROINIT(N) .eq. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': already initialised: '//trim(XTNAME(N))
          cycle
       end if

       !// Look for XSTT variable in file XSTT
       varname = 'XSTT_'//trim(XTNAME(N))

       status = nf90_inq_varid(ncid, TRIM(varname), var_id)
       if (status .eq. NF90_NOERR) then
          !// Variable exists - get it
          call getField_and_interpolate(ncid, var_id, filename, &
               varname, f90file//':'//subr, &
               IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
               IN_ETAA, IN_ETAB, LREGRID, R8XYZ)

          !// Set XSTT and XZEROINIT
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XSTT(L,N,I,J) = R8XYZ(I,J,L)
                end do
             end do
          end do
          XZEROINIT(N) = 0

          write(6,'(a)') f90file//':'//subr// &
               ': got variable (XSTT) : '//trim(XTNAME(N))

       else if (status .eq. NF90_ENOTVAR) then
          !// Variable not found - look in STT
          varname = 'STT_'//trim(XTNAME(N))

          status = nf90_inq_varid(ncid, TRIM(varname), var_id)
          if (status .eq. NF90_NOERR) then
             !// Variable exists - get it
             call getField_and_interpolate(ncid, var_id, filename, &
                  varname, f90file//':'//subr, &
                  IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
                  IN_ETAA, IN_ETAB, LREGRID, R8XYZ)

             !// Set XSTT and XZEROINIT
             do J = 1, JPAR
                do I = 1, IPAR
                   do L = 1, LPAR
                      XSTT(L,N,I,J) = R8XYZ(I,J,L)
                   end do
                end do
             end do
             XZEROINIT(N) = 0

             write(6,'(a)') f90file//':'//subr// &
                  ': got variable (STT->ZSTT) : '//trim(XTNAME(N))

          end if !// if (status .eq. NF90_NOERR) then
       else
          !// Variable does not exist on file - print info
          write(6,'(a)') f90file//':'//subr// &
               ': no such variable: '//trim(XTNAME(N))
       end if

    end do !// do N = 1, NOTRPAR


    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------
    deallocate( IN_XDEDG, IN_YDEDG, IN_ETAA, IN_ETAB )
    !//---------------------------------------------------------------------


    !// Print out info about the initialization field
    !// --------------------------------------------------------------------
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// Transported species
    write(6,'(a)') 'Transported species:'
    do N = 1, NPAR
       !// Print out only if initialized
       if (ZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                SUMT = SUMT + STT(I,J,L,N)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,TNAME(N),TMASS(N),SUMT
    end do

    !// Non-transported species
    write(6,'(a)') 'Non-transported species:'
    do N = 1, NOTRPAR
       !// Print out only if initialized
       if (XZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                SUMT = SUMT + XSTT(L,N,I,J)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,XTNAME(N),XTMASS(N),SUMT
    end do

    !// Which components are initialized as zero?
    do N = 1, NPAR
       if (ZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** Component trsp_idx/chem_idx ', &
               N,chem_idx(N),' is initialized as ZERO!'
       end if
    end do
    do N = 1, NOTRPAR
       if (XZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** X-Component Xtrsp_idx/Xchem_idx ', &
               N,Xchem_idx(N),' is initialized as ZERO!'
       end if
    end do


    !// Info if your restart file conains more/less species
    ND = sum(ZEROINIT)
    if (NOTRPAR .gt. 0) then
       XND = sum(XZEROINIT)
    else
       XND = 0
    endif
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    if ((ND + XND) .gt. 0) then
       write(6,'(a,i3)') '  Transported species not initialized: ', ND
       write(6,'(a,i3)') '  Non-Transported species not initialized: ', XND
       write(6,'(a)') 'Program _may_ stop if critical components'// &
            ' are not initialized!'
    else
       write(6,'(a)') '  All species are initialized!'
    end if
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'


    !// --------------------------------------------------------------------
  end subroutine load_restart_file
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine getField_and_interpolate(ncid, var_id, filename, &
               varname, caller, &
               IN_IPAR, IN_JPAR, IN_LPAR, IN_XDEDG, IN_YDEDG, &
               IN_ETAA, IN_ETAB, LREGRID, R8XYZ)
    !// --------------------------------------------------------------------
    !// This routine is bound to subroutine load_restart_file, and reads
    !// a 3d field from a netCDF file, with ncid as file id and
    !// var_id as variable id.
    !// If resolution on file differs from model resolution, the field is
    !// interpolated to model resoluion.
    !// Vertical interpolation is not possible yet.
    !//
    !// Amund Sovde Haslerud, August 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR
    use cmn_ctm, only: XDEDG, YDEDG
    use regridding, only: E_GRID
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ncid, var_id, IN_IPAR, IN_JPAR, IN_LPAR
    character(len=*), intent(in) :: filename, varname, caller
    real(r8), intent(in) :: IN_XDEDG(IN_IPAR+1), IN_YDEDG(IN_JPAR+1)
    real(r8), intent(in) :: IN_ETAA(IN_LPAR+1), IN_ETAB(IN_LPAR+1)
    logical, intent(in) :: LREGRID
    !// Output
    real(r8), intent(out) :: R8XYZ(IPAR,JPAR,LPAR)

    !// Local variables
    integer :: status, I, J, L
    integer, dimension(3) :: dimIDs3d, dimLen3d
    real(r8) :: IN_FIELD(IN_IPAR, IN_JPAR, IN_LPAR)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'getField_and_interpolate'
    !// --------------------------------------------------------------------

    !// Initialise
    R8XYZ(:,:,:) = 0._r8

    !// Variable information
    status = nf90_inquire_variable(ncid, var_id, dimIDs=dimIDs3d)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)// &
         ': inquire variable '//trim(varname)//' called by '//caller)

    !// 1. Dimension of the variable
    status = nf90_inquire_dimension(ncid,dimIDs3d(1),len=dimLen3d(1))
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)// &
         ': inquire dimesions 1 '//trim(varname)//' called by '//caller)
    !// 2. Dimension of the variable
    status = nf90_inquire_dimension(ncid,dimIDs3d(2),len=dimLen3d(2))
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)// &
         ': inquire dimesions 2 '//TRIM(varname)//' called by '//caller)
    !// 3. Dimension of the variable
    status = nf90_inquire_dimension(ncid,dimIDs3d(3),len=dimLen3d(3))
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)// &
         ': inquire dimesions 3 '//TRIM(varname)//' called by '//caller)

    !// Check the dimensions
    if ( (dimLen3d(1) .ne. IN_IPAR) .or. (dimLen3d(2) .ne. IN_JPAR) .or. &
         (dimLen3d(3) .ne. IN_LPAR) ) then
       write(6,'(a)') f90file//':'//subr//': resolution mismatch - stop'
       write(6,'(a,3i5)') 'On file',dimLen3d
       write(6,'(a,3i5)') 'On file',IN_IPAR, IN_JPAR, IN_LPAR
       write(6,'(a)') 'Routine called by '//caller
       stop
    end if

    !// Get variable
    status = nf90_get_var(ncid, var_id, IN_FIELD)
    if (status .ne. NF90_NOERR) call handle_error(status, &
         f90file//':'//subr//':'//trim(filename)// &
         ': get '//trim(varname)//' called by '//caller)


    !// Will so far only handle files where vertical resolution is the same
    !// as in model run.
    if (IN_LPAR .ne. LPAR) then
       write(6,'(a,2(1x,i3))') f90file//':'//subr// &
            ': Vertical resolution of restart file must match! - stopping'
       write(6,'(a,2(1x,i3))') 'restart file/IN_LPAR/LPAR: '//trim(filename),&
            IN_LPAR, LPAR
       write(6,'(a)') 'Routine called by '//caller
       stop
    end if



    !// Regrid horizontally?
    !// IMPORTANT: L-loop must be revised if vertical interpolation is
    !//            included.
    if (.not. LREGRID) then
       !// No regridding
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                R8XYZ(I,J,L) = IN_FIELD(I,J,L)
             end do
          end do
       end do
    else
       !// Must regrid horizontally
!$omp parallel private (L) &
!$omp          shared (IN_XDEDG,IN_YDEDG,XDEDG,YDEDG,IN_FIELD,R8XYZ,&
!$omp                  IN_IPAR, IN_JPAR, IN_LPAR) &
!$omp          default(NONE)
!$omp do
       do L = 1, IN_LPAR
          !// Interpolate
          call E_GRID(IN_FIELD(:,:,L),IN_XDEDG,IN_YDEDG,IN_IPAR,IN_JPAR, &
               R8XYZ(:,:,L),XDEDG,YDEDG,IPAR,JPAR,1)
       end do
!$omp end do
!$omp end parallel

       !// Possible future vertical regridding


    end if

    !// --------------------------------------------------------------------
  end subroutine getField_and_interpolate
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine save_restart_file(NDAY, NDAYI)
    !// --------------------------------------------------------------------
    !// Save restart file for exact continuation run.
    !// Written as netCDF4.
    !//
    !// Amund Sovde Haslerud, August 2017
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPARW, JPARW, LPARW, IPAR, JPAR, LPAR, NPAR, NOTRPAR, &
         TNAMELEN
    use cmn_met, only: METTYPE, metCYCLE, metREVNR, MET_ROOT, MPATH1, MPATH2, &
         MYEAR
    use cmn_ctm, only: IYEAR, JYEAR, JDAY, JMON, JDATE, IDAY, NTM, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         XDGRD, YDGRD, XDEDG, YDEDG, ZGRD, ZEDG, ETAA, ETAB, AREAXY
    use cmn_chem, only: TNAME, TMASS
    use cmn_diag, only: RUNTITLE, metTypeInfo, resolutionInfo, &
         nc4deflate_global, nc4shuffle_global
    use cmn_oslo, only: chem_idx, Xchem_idx, XSTT, XTNAME, XTMASS
    use netcdf
    use ncutils, only: handle_error
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NDAYI

    !// Locals
    character(len=8) :: datestamp
    character(len=10) :: CDATE
    character(len=80) :: filename
    integer :: N, I, J, L
    real(r8),dimension(IPAR,JPAR,LPAR)  :: XFIELD
    integer, dimension(6) :: output_time

    integer :: &
         lat_dim_id, lon_dim_id, lev_dim_id, npar_dim_id, notrpar_dim_id, &
         lat_id, lon_id, lev_id, &
         ilat_dim_id, ilon_dim_id, ilev_dim_id, &
         ilat_id, ilon_id, ilev_id, &
         ihya_dim_id, ihyb_dim_id, &
         ihya_id, ihyb_id, &
         native_lon_id, native_lat_id, native_lev_id, &
         areaxy_id, &
         version_id, &
         date_size_dim_id, &
         time_id, nday_id, iyear_id, &
         tracer_name_len_dim_id, tracer_name_id, xtracer_name_id, &
         chemid_id, xchemid_id, molw_id, xmolw_id, &
         air_id, &
         status, ncid
    integer, dimension(NPAR) :: comp_id, &
         compSUT_id, compSVT_id, compSWT_id, &
         compSUU_id, compSVV_id, compSWW_id, &
         compSUV_id, compSUW_id, compSVW_id
    integer, dimension(NOTRPAR) :: Xcomp_id
    character(len=20) :: varname

    !// Version number
    integer, parameter :: version = 4
    !// --------------------------------------------------------------------
    integer, parameter :: nc4deflate = 9 !// Could use nc4deflate_global
    integer, parameter :: nc4shuffle = 1 !// Could use nc4shuffle_global
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CON_SAV_NC4'
    !// --------------------------------------------------------------------

    !// Date stamp
    write(datestamp(1:8),'(i4.4,2i2.2)') JYEAR,JMON,JDATE
    !// File to save
    filename = 'ctm3_restart_'//datestamp//'.nc'

    output_time = (/JYEAR,JMON,JDATE,0,0,0/)


    !//---------------------------------------------------------------------
    write(6,'(a)') f90file//':'//subr//': creating file: '//trim(filename)
    status=nf90_create(path=filename,cmode=nf90_netcdf4,ncid=ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': in creating file')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !// Headers
    !//---------------------------------------------------------------------
    status=nf90_put_att(ncid,nf90_global,'model_info', &
         'Oslo CTM3 is a 3D Chemical Transport Model')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': model_info')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology',metTypeInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology')
    status=nf90_put_att(ncid,nf90_global,'driving_meteorology_path',MPATH1)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': driving_meteorology_path')
    status=nf90_put_att(ncid,nf90_global,'resolution_info',resolutionInfo)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': resolution_info')
    status=nf90_put_att(ncid,nf90_global,'runtitle',trim(runtitle))
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': runtitle')
    status=nf90_put_att(ncid,nf90_global,'contact_info','CICERO')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': contact_info')
    !//---------------------------------------------------------------------

    !//---------------------------------------------------------------------
    !// Define dimensions
    !//---------------------------------------------------------------------
    !// Define spatial dimensions (lat, lon, lev)
    status = nf90_def_dim(ncid,"lat",JPAR,lat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat dim')
    status = nf90_def_dim(ncid,"lon",IPAR,lon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon dim')
    status = nf90_def_dim(ncid,"lev",LPAR,lev_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev dim')
    !// Define spatial dimensions (ilat, ilon, ilev)
    status = nf90_def_dim(ncid,"ilat",JPAR+1,ilat_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat dim')
    status = nf90_def_dim(ncid,"ilon",IPAR+1,ilon_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon dim')
    status = nf90_def_dim(ncid,"ilev",LPAR+1,ilev_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev dim')

    !// NPAR
    status = nf90_def_dim(ncid,"NPAR",NPAR,npar_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NPAR dim')
    !// NOTRPAR
    status = nf90_def_dim(ncid,"NOTRPAR",NOTRPAR,notrpar_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NOTRPAR dim')

    !// Define size of date stamps
    status = nf90_def_dim(ncid,"date_size",size(output_time),date_size_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define date_size dim')

    !// Define length of tracer name string
    status = nf90_def_dim(ncid,"tracer_name_len", &
         TNAMELEN, tracer_name_len_dim_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name_len dim')


    !//---------------------------------------------------------------------
    !// Define variables
    !//---------------------------------------------------------------------
    !// Define the lon/lat/lev
    status = nf90_def_var(ncid,"lon",nf90_double,lon_dim_id,lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lon variable')
    status = nf90_def_var(ncid,"lat",nf90_double,lat_dim_id,lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lat variable')
    status = nf90_def_var(ncid,"lev",nf90_double,lev_dim_id,lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define lev variable')

    !// Define ilon/ilat/ilev
    status = nf90_def_var(ncid,"ilon",nf90_double,ilon_dim_id,ilon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilon variable')
    status = nf90_def_var(ncid,"ilat",nf90_double,ilat_dim_id,ilat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilat variable')
    status = nf90_def_var(ncid,"ilev",nf90_double,ilev_dim_id,ilev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ilev variable')

    !// Putting attributes to lon/lat/lev variables
    status = nf90_put_att(ncid,lon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lon')
    status = nf90_put_att(ncid,lat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lat')
    status = nf90_put_att(ncid,lev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units lev')
    !// Putting attributes to ilon/ilat variables
    status = nf90_put_att(ncid,ilon_id,'units','degree_east')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilon')
    status = nf90_put_att(ncid,ilat_id,'units','degree_north')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilat')
    status = nf90_put_att(ncid,ilev_id,'units','pressure [hPa]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ilev')

    !// Defining hybrid sigma coordinate A
    status = nf90_def_var(ncid,"ihya",nf90_double,ilev_dim_id,ihya_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihya variable')
    !// Attributes
    status = nf90_put_att(ncid,ihya_id,'units','hPa')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihya')
    status = nf90_put_att(ncid,ihya_id,'description', &
         'Sigma hybrid coordinate A.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihya')
    status = nf90_put_att(ncid,ihya_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihya')

    !// Defining hybrid sigma coordinate B
    status = nf90_def_var(ncid,"ihyb",nf90_double,ilev_dim_id,ihyb_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define ihyb variable')
    !// Attributes
    status = nf90_put_att(ncid,ihyb_id,'units','1')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'description', &
         'Sigma hybrid coordinate B.')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description ihyb')
    status = nf90_put_att(ncid,ihyb_id, 'usage', &
         'p_box_bottom(L) = ihya(L) + ihyb(L)*p_surface(I,J)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute usage ihyb')


    !// Define file version number
    status = nf90_def_var(ncid,"VERSION",nf90_int,version_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define VERSION variable')
    status = nf90_put_att(ncid,version_id,'description', &
         'Restart file version number')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description VERSION')

    !// Define output time
    status = nf90_def_var(ncid,"output_time",nf90_int,date_size_dim_id,time_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define output_time variable')
    status = nf90_put_att(ncid,time_id,'description', &
         'Output time [YYYY,MM,DD,hh,mm,ss]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description output_time')

    !// Put out NDAYS and IYEAR to indicate what kind of spinup
    status = nf90_def_var(ncid,"NDAYS",nf90_int,nday_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define NDAYS variable')
    status = nf90_put_att(ncid,nday_id,'description', &
         'Numper of days since model start (NDAY-NDAYI+1)')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description NDAYS')

    status = nf90_def_var(ncid,"IYEAR",nf90_int,iyear_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IYEAR variable')
    status = nf90_put_att(ncid,iyear_id,'description', &
         'Initial year of model simulation')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IYEAR')

    !// Info about native size IPARW
    status = nf90_def_var(ncid,"IPARW",nf90_int,native_lon_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define IPARW variable')
    status = nf90_put_att(ncid,native_lon_id,'description', &
         'Meteorological data native longitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description IPARW')
    !// Info about native size JPARW
    status = nf90_def_var(ncid,"JPARW",nf90_int,native_lat_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define JPARW variable')
    status = nf90_put_att(ncid,native_lat_id,'description', &
         'Meteorological data native latitudinal resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description JPARW')
    !// Info about native size LPARW
    status = nf90_def_var(ncid,"LPARW",nf90_int,native_lev_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define LPARW variable')
    status = nf90_put_att(ncid,native_lev_id,'description', &
         'Meteorological data vertical resolution')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description LPARW')

    !// Tracer names
    status = nf90_def_var(ncid, "tracer_name", nf90_char, &
         (/tracer_name_len_dim_id, npar_dim_id/), tracer_name_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_name variable')
    status = nf90_put_att(ncid,tracer_name_id,'description', &
         'Tracer names - transported species')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_name')

    !// Component IDs
    status = nf90_def_var(ncid, "chem_idx", nf90_int, &
         (/npar_dim_id/), chemid_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define chem_idx variable')
    status = nf90_put_att(ncid,chemid_id,'description', &
         'Component ID in Oslo CTM3 - transported species')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description chem_idx')

    !// Tracer mass (molecular weight)
    status = nf90_def_var(ncid, "tracer_molweight", nf90_double, &
         (/npar_dim_id/), molw_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define tracer_molweight variable')
    status = nf90_put_att(ncid,molw_id,'description', &
         'Tracer mlecular weights - transported species')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description tracer_molweight')
    status = nf90_put_att(ncid,molw_id,'units','g/mol')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units tracer_molweight')

    if (NOTRPAR .gt. 0) then
       !// Tracer names - non-transported species
       status = nf90_def_var(ncid, "Xtracer_name", nf90_char, &
         (/tracer_name_len_dim_id, notrpar_dim_id/), xtracer_name_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define Xtracer_name variable')
       status = nf90_put_att(ncid,xtracer_name_id,'description', &
            'Tracer names - non-transported species')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description Xtracer_name')

       !// Component IDs - non-transported species
       status = nf90_def_var(ncid, "Xchem_idx", nf90_int, &
            (/notrpar_dim_id/), xchemid_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define Xchem_idx variable')
       status = nf90_put_att(ncid,xchemid_id,'description', &
            'Component ID in Oslo CTM3 - non-transported species')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description Xchem_idx')

       !// Tracer mass (molecular weight) - non-transported species
       status = nf90_def_var(ncid, "Xtracer_molweight", nf90_double, &
            (/notrpar_dim_id/), xmolw_id)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define Xtracer_molweight variable')
       status = nf90_put_att(ncid,xmolw_id,'description', &
            'Tracer mlecular weights - non-transported species')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute description Xtracer_molweight')
       status = nf90_put_att(ncid,xmolw_id,'units','g/mol')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units Xtracer_molweight')

    end if !// if (NOTRPAR .gt. 0) then

    !// Grid area
    status = nf90_def_var(ncid, "gridarea", nf90_double, &
         (/lon_dim_id, lat_dim_id /), areaxy_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define gridarea variable')
    status = nf90_def_var_deflate(ncid,areaxy_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate gridarea variable')
    status = nf90_put_att(ncid,areaxy_id,'unit','m2')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute unit gridarea')

    !//---------------------------------------------------------------------
    !// Define AIR and tracers
    !//---------------------------------------------------------------------

    !// Define AIR
    status = nf90_def_var(ncid, "AIR", nf90_double, &
         (/ lon_dim_id, lat_dim_id, lev_dim_id /), air_id)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define AIR variable')
    status = nf90_def_var_deflate(ncid,air_id,nc4shuffle,1,nc4deflate)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': define deflate AIR variable')
    status = nf90_put_att(ncid, air_id, 'units', 'kg')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute units AIR')
    status = nf90_put_att(ncid, air_id, 'description', &
         'Air mass [kg of dry air]')
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': attribute description AIR')

    !// Transported species - STT
    do N = 1, NPAR
       varname = 'STT_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_double, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), comp_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,comp_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, comp_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, comp_id(N), 'long_name', &
            'Tracer mass [kg] of '//trim(TNAME(N))//' - STT array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SUT
    do N = 1, NPAR
       varname = 'SUT_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSUT_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSUT_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSUT_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSUT_id(N), 'long_name', &
            'Tracer moment SUT [kg] of '//trim(TNAME(N))//' - SUT array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SVT
    do N = 1, NPAR
       varname = 'SVT_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSVT_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSVT_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSVT_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSVT_id(N), 'long_name', &
            'Tracer moment SVT [kg] of '//trim(TNAME(N))//' - SVT array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SWT
    do N = 1, NPAR
       varname = 'SWT_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSWT_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSWT_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSWT_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSWT_id(N), 'long_name', &
            'Tracer moment SWT [kg] of '//trim(TNAME(N))//' - SWT array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SUU
    do N = 1, NPAR
       varname = 'SUU_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSUU_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSUU_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSUU_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSUU_id(N), 'long_name', &
            'Tracer moment SUU [kg] of '//trim(TNAME(N))//' - SUU array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SVV
    do N = 1, NPAR
       varname = 'SVV_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSVV_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSVV_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSVV_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSVV_id(N), 'long_name', &
            'Tracer moment SVV [kg] of '//trim(TNAME(N))//' - SVV array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SWW
    do N = 1, NPAR
       varname = 'SWW_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSWW_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSWW_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSWW_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSWW_id(N), 'long_name', &
            'Tracer moment SWW [kg] of '//trim(TNAME(N))//' - SWW array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SUV
    do N = 1, NPAR
       varname = 'SUV_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSUV_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSUV_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSUV_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSUV_id(N), 'long_name', &
            'Tracer moment SUV [kg] of '//trim(TNAME(N))//' - SUV array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SUW
    do N = 1, NPAR
       varname = 'SUW_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSUW_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSUW_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSUW_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSUW_id(N), 'long_name', &
            'Tracer moment SUW [kg] of '//trim(TNAME(N))//' - SUW array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR

    !// Transported species - SVW
    do N = 1, NPAR
       varname = 'SVW_'//trim(TNAME(N))
       status = nf90_def_var(ncid, trim(varname), nf90_float, &
            (/ lon_dim_id, lat_dim_id, lev_dim_id /), compSVW_id(N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define variable '//trim(varname))
       status = nf90_def_var_deflate(ncid,compSVW_id(N),nc4shuffle,1,nc4deflate)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': define deflate variable '//trim(varname))
       status = nf90_put_att(ncid, compSVW_id(N), 'units', 'kg')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute units '//trim(varname))
       status = nf90_put_att(ncid, compSVW_id(N), 'long_name', &
            'Tracer moment SVW [kg] of '//trim(TNAME(N))//' - SVW array')
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': attribute long_name '//trim(varname))
    end do !// do N = 1, NPAR



    !// Non-transported species
    if (NOTRPAR .gt. 0) then
       !// Only put out if NOTRPAR > 0
       do N = 1, NOTRPAR
          varname = 'XSTT_'//trim(XTNAME(N))
          status = nf90_def_var(ncid, trim(varname), nf90_double, &
               (/ lon_dim_id, lat_dim_id, lev_dim_id /), xcomp_id(N))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define variable '//trim(varname))
          status = nf90_def_var_deflate(ncid,xcomp_id(N),nc4shuffle,1,nc4deflate)
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': define deflate variable '//trim(varname))
          status = nf90_put_att(ncid, xcomp_id(N), 'units', 'kg')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute units '//trim(varname))
          status = nf90_put_att(ncid, xcomp_id(N), 'long_name', &
               'Tracer mass [kg] of '//trim(XTNAME(N))//' - XSTT array')
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': attribute long_name '//trim(varname))
       end do !// do N = 1, NOTRPAR
    end if !// if (NOTRPAR .gt. 0) then



    !//---------------------------------------------------------------------
    !// End definition mode
    status = nf90_enddef(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': end defmode')
    !//---------------------------------------------------------------------

    !// Putting the lon/lat/lev variables
    status = nf90_put_var(ncid,lon_id,XDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lon')
    status = nf90_put_var(ncid,lat_id,YDGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lat')
    status = nf90_put_var(ncid,lev_id,ZGRD)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting lev')

    !// Putting the ilon/ilat/ilev variables
    status = nf90_put_var(ncid,ilon_id,XDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilon')
    status = nf90_put_var(ncid,ilat_id,YDEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilat')
    status = nf90_put_var(ncid,ilev_id,ZEDG)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ilev')

    !// Putting the ihya/ihyb variables
    status = nf90_put_var(ncid,ihya_id,ETAA)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihya')
    status = nf90_put_var(ncid,ihyb_id,ETAB)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting ihyb')

    !// Putting file version number
    status = nf90_put_var(ncid, version_id, version)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting VERSION')

    !// Output time
    status = nf90_put_var(ncid,time_id,output_time)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting output_time')

    !// NDAYS
    status = nf90_put_var(ncid,nday_id,NDAY+1-NDAYI)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting NDAYS')
    !// IYEAR
    status = nf90_put_var(ncid,iyear_id,IYEAR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IYEAR')

    !// Info about native size
    status = nf90_put_var(ncid,native_lon_id,IPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting IPARW')
    status = nf90_put_var(ncid,native_lat_id,JPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting JPARW')
    status = nf90_put_var(ncid,native_lev_id,LPARW)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting LPARW')

    !// Tracer names
    status = nf90_put_var(ncid,tracer_name_id,TNAME)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_name')

    !// Component ID
    status = nf90_put_var(ncid,chemid_id,chem_idx)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting chem_idx')

    !// Tracer molecular weights
    status = nf90_put_var(ncid,molw_id,TMASS)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting tracer_molweight')

    if (NOTRPAR .gt. 0) then
       !// Tracer names - non-transported species
       status = nf90_put_var(ncid,xtracer_name_id,XTNAME)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting Xtracer_name')

       !// Component ID
       status = nf90_put_var(ncid,xchemid_id,Xchem_idx)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting Xchem_idx')

       !// Tracer molecular weights
       status = nf90_put_var(ncid,xmolw_id,XTMASS)
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting Xtracer_molweight')

    end if !// if (NOTRPAR .gt. 0) then

    !// Put Grid area
    status = nf90_put_var(ncid,areaxy_id,AREAXY)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AREAXY')

    !// Put AIR
    status = nf90_put_var(ncid,air_id,AIR)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': putting AIR')

    !// Put STT
    do N = 1, NPAR
       status = nf90_put_var(ncid,comp_id(N),STT(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting STT_'//trim(TNAME(N)))
    end do

    !// Put SUT moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSUT_id(N),SUT(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SUT_'//trim(TNAME(N)))
    end do

    !// Put SVT moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSVT_id(N),SVT(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SVT_'//trim(TNAME(N)))
    end do

    !// Put SWT moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSWT_id(N),SWT(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SWT_'//trim(TNAME(N)))
    end do

    !// Put SUU moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSUU_id(N),SUU(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SUU_'//trim(TNAME(N)))
    end do

    !// Put SVV moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSVV_id(N),SVV(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SVV_'//trim(TNAME(N)))
    end do

    !// Put SWW moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSWW_id(N),SWW(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SWW_'//trim(TNAME(N)))
    end do

    !// Put SUV moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSUV_id(N),SUV(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SUV_'//trim(TNAME(N)))
    end do

    !// Put SUW moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSUW_id(N),SUW(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SUW_'//trim(TNAME(N)))
    end do

    !// Put SVW moment
    do N = 1, NPAR
       status = nf90_put_var(ncid,compSVW_id(N),SVW(:,:,:,N))
       if (status .ne. nf90_noerr) call handle_error(status, &
            f90file//':'//subr//': putting SVW_'//trim(TNAME(N)))
    end do


    !// Put XSTT - non-transported species
    if (NOTRPAR .gt. 0) then
       do N = 1, NOTRPAR
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XFIELD(I,J,L) = XSTT(L,N,I,J)
                end do
             end do
          end do
          status = nf90_put_var(ncid,xcomp_id(N),XFIELD(:,:,:))
          if (status .ne. nf90_noerr) call handle_error(status, &
               f90file//':'//subr//': putting XSTT_'//trim(XTNAME(N)))
       end do
    end if !// if (NOTRPAR .gt. 0) then

    


    !//---------------------------------------------------------------------
    !// close netcdf file
    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_error(status, &
         f90file//':'//subr//': close file: '//trim(filename))
    !//---------------------------------------------------------------------

    write(6,'(a)') f90file//':'//subr//': wrote file: '//trim(filename)

    !// --------------------------------------------------------------------
  end subroutine save_restart_file
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine OSLO_CON_SAV(NDAY)
    !// --------------------------------------------------------------------
    !// Save restart file for exact continuation run.
    !// UCI CON_SAV modified to include non-transported Oslo components.
    !// This version writes species by species, for a more flexible
    !// treatment.
    !//
    !// Ole Amund Sovde, October 2015, January 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPARW, JPARW, LPARW, IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_met, only: MYEAR
    use cmn_ctm, only: IYEAR, JYEAR, JDAY, JMON, JDATE, IDAY, NTM, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         XDEDG, YDEDG, ETAA, ETAB
    use cmn_diag, only: RUNTITLE
    use utilities, only: get_free_fileid
    use cmn_oslo, only: chem_idx, Xchem_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY

    !// Locals
    character(len=10) :: CDATE
    character(len=80) :: FN_CON, TITLE
    integer :: ioerr, N, file_id
    real(r8),dimension(IPAR,JPAR,LPAR)  :: XFIELD

    !// Version number
    integer, parameter :: version = 3
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CON_SAV'
    !// --------------------------------------------------------------------

    !// Define file name for restart/avgs from NDAYW (allow for 100-yr run !)
    write (CDATE(1:10),'(i5.5,a1,i4.4)')  NDAY,'_',MYEAR
    !// File to save
    FN_CON = 'ctm3_restart_'//CDATE//'.sav'
        
    !// Find a file number
    file_id = get_free_fileid()

    open(file_id,file=FN_CON,form='UNFORMATTED',status='new',iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr//&
            ': Error opening file: '//trim(fn_con)
       stop
    end if

    write(file_id) RUNTITLE

    !// General information and the transported data
    write(file_id) IYEAR,NDAY, JYEAR,JDAY,JMON,JDATE,version
    write(file_id) IPARW,JPARW,LPARW, IPAR,JPAR,LPAR,NPAR,NOTRPAR
    write(file_id) XDEDG, YDEDG, ETAA, ETAB
    write(file_id) AIR
    !// Write species by species. Write one array at a time (avoids possible
    !// too-large records)
    do N = 1, NPAR
       write(file_id) chem_idx(N)
       write(file_id) STT(:,:,:,N)
       write(file_id) SUT(:,:,:,N)
       write(file_id) SVT(:,:,:,N)
       write(file_id) SWT(:,:,:,N)
       write(file_id) SUU(:,:,:,N)
       write(file_id) SVV(:,:,:,N)
       write(file_id) SWW(:,:,:,N)
       write(file_id) SUV(:,:,:,N)
       write(file_id) SUW(:,:,:,N)
       write(file_id) SVW(:,:,:,N)
    end do

    !// Non-transported data
    do N = 1,NOTRPAR
       write(file_id) Xchem_idx(N)
       write(file_id) XSTT(:,N,:,:)
    end do

    close(file_id)
    write(6,'(a,i8)') 'Wrote restart file: '//trim(FN_CON),IDAY

    !// --------------------------------------------------------------------
  end subroutine OSLO_CON_SAV
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine OSLO_CON_RUN(FN_CON, MODE)
    !// --------------------------------------------------------------------
    !// Read restart file and set up correct day information.
    !// UCI CON_RUN modified to include non-transported Oslo components.
    !// Reads component by component, and places the data in the current
    !// array. Works fine when current components are less than in file.
    !//
    !// Can read ADDITIONAL restart files to fill in non-initialized
    !// tracers. This is controlled by MODE:
    !// MODE:  0: First call; initializes zeroinit.
    !//      !=0: Will add non-initialized tracers that are
    !//           available on the specified file. Skips AIR.
    !//
    !// Ole Amund Sovde, October 2015, June 2010, January 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_ctm, only: IYEAR, IDAY, JYEAR, JDAY, JMON, JDATE, NTM, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW
    use cmn_chem, only: TNAME, TMASS
    use utilities, only: get_free_fileid
    use cmn_oslo, only: ZEROINIT, XZEROINIT, trsp_idx, chem_idx, &
          XTNAME, XTMASS, Xtrsp_idx, Xchem_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in):: FN_CON
    integer, intent(in) :: MODE

    !// Locals
    real(r8) :: SUMT
    real(r8), dimension(IPAR,JPAR,LPAR) :: INFIELD
    real(r8), dimension(LPAR,IPAR,JPAR) :: XINFIELD
    real(rMom), dimension(IPAR,JPAR,LPAR) :: &
         TSUT, TSVT, TSWT, TSUU, TSVV, TSWW, TSUV, TSUW, TSVW
    real(r8), dimension(IPAR+1) :: INXDEDG
    real(r8), dimension(JPAR+1) :: INYDEDG
    real(r8), dimension(LPAR+1) :: INETAA, INETAB

    integer :: &
         IIYEAR, IIDAY, JJYEAR, JJDAY, JJMON, JJDATE, version, &
         I,J,L,N, IDW,JDW,LDW, ID,JD,LD, ND, XND, &
         ioerr, file_id, CID, TRID, XTRID
    logical :: LERR1, LERR2, LERR3, file_io
    character(len=80) :: RTITLE2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CON_RUN'
    !// --------------------------------------------------------------------

    !// Initialize arrays; assume all are zero
    if (MODE .eq. 0) then
       ZEROINIT(:)  = 1
       if (NOTRPAR .gt. 0) XZEROINIT(:) = 1
    end if

    !//read in continuation restart file:
    if (MODE .eq. 0) then
       write(6,'(a)') '* Read restart file: '//trim(FN_CON)
    else
       write(6,'(a)') '* Read ADDITIONAL restart file: '//trim(FN_CON)
    end if

    !// Find a file number
    file_id = get_free_fileid()

    open(file_id,file=FN_CON,status='OLD',form='UNFORMATTED', &
         iostat=ioerr,action='read')
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Error opening continuation file: '//trim(FN_CON)
       stop
    end if

    read(file_id) RTITLE2
    write(6,'(2a)') 'CTM run: ',trim(RTITLE2)

    !// Read dimensions and transported tracers
    read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE,version
    if (ioerr .ne. 0) then
       !// Version is not specified, meaning version 2
       version = 2
       backspace(file_id)
       read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    end if
    write(6,'(a,6i6)') 'Model date:',IYEAR,IDAY,JYEAR,JDAY,JMON,JDATE
    write(6,'(a,6i6)') 'File date :',IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    write(6,'(a,i6)')  'File version :',version
    !// It is possible to override the date information here,
    !// but I don't see why you would.

    if (version .le. 2) then
       read(file_id,iostat=ioerr) IDW,JDW,LDW,I,ID,JD,LD,ND,XND
    else
       !// Read native resolution and model resolution
       read(file_id,iostat=ioerr) IDW,JDW,LDW, ID,JD,LD, ND,XND
    end if

    if (ioerr .ne. 0) then
       write(6,'(a,8(1x,i3))') f90file//':'//subr// &
            ': Problems reading IDW,JDW,LDW,ID,JD,LD,ND,XND', &
            IDW,JDW,LDW, ID,JD,LD, ND,XND
       stop
    end if
    !// Check dimensions
    if (ID .ne. IPAR .or. JD .ne. JPAR .or. LD .ne. LPAR) then
       write(6,'(a,6(1x,i3))') f90file//':'//subr// &
            ': Wrong resolution of restart file: '//trim(FN_CON), &
            ID,JD,LD,IPAR,JPAR,LPAR
       stop
    end if

    write(6,'(a,i4)') 'Transported species on file :    ',ND
    write(6,'(a,i4)') 'Non-transported species on file: ',XND

    !// Read grid box info for version > 2
    if (version .gt. 2) then
       read(file_id,iostat=ioerr) INXDEDG, INYDEDG, INETAA, INETAB
       if (ioerr .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Problems reading XDEDG/YDEDG/ETAA/ETAB'
          stop
       end if
    end if !// if (version .gt. 2) then

    !// Read air
    read(file_id,iostat=ioerr) INFIELD
    if (ioerr .ne. 0) then
       write(6,'(a)') 'utilities_oslo.f90: OSLO_CON_RUN: Problems reading AIR'
       stop
    end if
    if (MODE .eq. 0) then
       AIR(:,:,:) = INFIELD(:,:,:)
       write(6,'(a)') '* AIR read and set!'
    else
       write(6,'(a)') '* AIR read and skipped!'
    end if


    !// Read each of the transported components
    do N = 1, ND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading CID',CID
          stop
       end if
       !// Fetch transport id and non-transport id
       TRID  = trsp_idx(CID)
       XTRID = Xtrsp_idx(CID)

       !// STT
       read(file_id) INFIELD
       read(file_id) TSUT(:,:,:)
       read(file_id) TSVT(:,:,:)
       read(file_id) TSWT(:,:,:)
       read(file_id) TSUU(:,:,:)
       read(file_id) TSVV(:,:,:)
       read(file_id) TSWW(:,:,:)
       read(file_id) TSUV(:,:,:)
       read(file_id) TSUW(:,:,:)
       read(file_id) TSVW(:,:,:)


       !// Place component if included in the run
       if (TRID .gt. 0) then
          !// Check if already initialized
          if (ZEROINIT(TRID) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          STT(:,:,:,TRID) = INFIELD(:,:,:)
          SUT(:,:,:,TRID) = TSUT(:,:,:)
          SVT(:,:,:,TRID) = TSVT(:,:,:)
          SWT(:,:,:,TRID) = TSWT(:,:,:)
          SUU(:,:,:,TRID) = TSUU(:,:,:)
          SVV(:,:,:,TRID) = TSVV(:,:,:)
          SWW(:,:,:,TRID) = TSWW(:,:,:)
          SUV(:,:,:,TRID) = TSUV(:,:,:)
          SUW(:,:,:,TRID) = TSUW(:,:,:)
          SVW(:,:,:,TRID) = TSVW(:,:,:)
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (transported)'
          !// Keep track of initialized components
          ZEROINIT(TRID) = 0
       else if (XTRID .gt. 0) then
          !// Check if already initialized
          if (XZEROINIT(XTRID) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          !// Put transported component into non-transported array
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XSTT(L,XTRID,I,J) = INFIELD(I,J,L)
                end do
             end do
          end do
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (not transported)'
          !// Keep track of initialized components
          XZEROINIT(XTRID) = 0
       else
          write(6,'(a,i3,a)') '  Skipping component ',CID,' (not included)'
       end if

    end do !//  do N = 1,ND

    !// Read each of the non-transported components
    do N = 1, XND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading X CID',CID
          stop
       end if
       !// Fetch transport id and non-transport id
       TRID  = trsp_idx(CID)
       XTRID = Xtrsp_idx(CID)

       !// XSTT
       read(file_id) XINFIELD

       !// Place component if included in the run
       if (XTRID .gt. 0) then
          !// Check if already initialized
          if (XZEROINIT(XTRID) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ',CID,' (not transported)'
          XSTT(:,XTRID,:,:) = XINFIELD(:,:,:)
          !// Keep track of initialized components
          XZEROINIT(XTRID) = 0
       else if (TRID .gt. 0) then
          !// Check if already initialized
          if (ZEROINIT(TRID) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ',CID,' (transported, zeroed moments)'
          do L = 1, LPAR
             do J = 1, JPAR
                do I = 1, IPAR
                   STT(I,J,L,TRID) = XINFIELD(L,I,J)
                end do
             end do
          end do
          !// Moments are set to zero
          !// Keep track of initialized components
          ZEROINIT(TRID) = 0
       else
          !// Ah, tracer is not included in this run
          write(6,'(a,i3,a)') '  Skipping X component ',CID,' (not included)'
       end if

    end do !//  do N = 1,XND

    close(file_id)
    write(6,'(a)') '* Finished reading restart file: '//trim(FN_CON)


    !// Print out info about the initialization field
    !// ------------------------------------------------------------------
    write(6,'(A,7I6)') 'JYEAR/JDATE/IDAY/IYEAR/JDAY/JMON', &
         JYEAR,JDATE,IDAY,IYEAR,JDAY,JMON
    write(6,'(a71)') '--------------------------------------------' &
         //'---------------------------'

    !// Transported species
    write(6,'(a)') 'Transported species:'
    do N = 1,NTM
       !// Print out only if initialized
       if (ZEROINIT(N) .eq. 1) cycle

       SUMT = 0._r8
       do L = 1,LPAR
          do J = 1,JPAR
             do I = 1,IPAR
                SUMT = SUMT + STT(I,J,L,N)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,TNAME(N),TMASS(N),SUMT
    end do

    !// Non-transported species
    write(6,'(a)') 'Non-transported species:'
    do N = 1, NOTRPAR
       !// Print out only if initialized
       if (XZEROINIT(N) .eq. 1) cycle

       SUMT = 0._r8
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                SUMT = SUMT + XSTT(L,N,I,J)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,XTNAME(N),XTMASS(N),SUMT
    end do

    !// Which components are initialized as zero?
    do N = 1, NTM
       if (ZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** Component trsp_idx/chem_idx ', &
               N,chem_idx(N),' is initialized as ZERO!'
       end if
    end do
    do N = 1, NOTRPAR
       if (XZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** X-Component Xtrsp_idx/Xchem_idx ', &
               N,Xchem_idx(N),' is initialized as ZERO!'
       end if
    end do


    !// Info if your restart file conains more/less species
    ND = sum(ZEROINIT)
    if (NOTRPAR .gt. 0) then
       XND = sum(XZEROINIT)
    else
       XND = 0
    end if
    write(6,'(a71)') '--------------------------------------------' &
         //'---------------------------'
    if ((ND+XND) .gt. 0) then
       write(6,'(a,i3)') '  Transported species not initialized: ', ND
       write(6,'(a,i3)') '  Non-Transported species not initialized: ', XND
       write(6,'(a)') 'Program _may_ stop if critical components'// &
            ' are not initialized!'
    else
       write(6,'(a)') '  All species are initialized!'
    end if
    write(6,'(a71)') '--------------------------------------------' &
         //'---------------------------'

    !// --------------------------------------------------------------------
  end subroutine OSLO_CON_RUN
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine restart_from_CTM3avg(zeroinit)
    !// --------------------------------------------------------------------
    !// Restart from Oslo CTM3 avgsav file, version 7.
    !// Reads a CTM3 file (L37/40/57/60) and collapses it simply if
    !// necessary (when possible): only possible for L37.
    !//
    !// Ole Amund Sovde, November 2013
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR, LPARW
    use cmn_ctm, only: STT, LCONT
    use cmn_chem, only: TNAME
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx, XTNAME, Xtrsp_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Output
    integer, intent(out) :: zeroinit(NPAR+NOTRPAR)

    !// Locals
    real(r4) :: RINDATA(IPAR,JPAR,LPAR)
    integer :: &
         IN_RES(4), IN_TIME(3), &
         version, IN_NOTRPAR, I, J, L, N, &
         TRSPNR, ioerr, II, JJ, LL, TRACER_ID

    character(len=80) :: CINPUT_RUN

    real(r8) :: RFAC !// for unit conversion

    !// for file handling
    logical :: ex
    integer :: fnr, err
    character(len=80) :: FN_AVG

    real(r4) :: intmp1
    logical :: inlog

    !// Read-in data, possibly different resolution
    integer, allocatable :: IN_IDMTC(:), IN_XIDMTC(:)
    real(r4), allocatable :: RXSTT(:,:,:,:)
    real(r8), allocatable :: IN_TMASS(:), IN_XTMASS(:)
    character(len=10), allocatable :: IN_TNAME(:), IN_XTNAME(:)
    !// --------------------------------------------------------------------

    !// Will read from this file
    FN_AVG = 'avgsav_input'

    if (LCONT) then
       write(6,'(a)') 'Check LCONT!'
       write(6,'(a)') 'Trying to read '//trim(FN_AVG)//' when LCONT=.true.'
       stop
    end if

    inquire(File=trim(FN_AVG),exist=ex)
    if (.not. ex) then
       write(6,'(a)') '** The file '//trim(FN_AVG)//' does not exist **'
       stop
    end if

    fnr = get_free_fileid()

    write(6,'(/a)') 'Reading initial tracer data from: '//trim(FN_AVG)
    !// Open the file
    open(fnr, File=trim(FN_AVG), Form='UNFORMATTED', Status='OLD', iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') '*** error opening avg-file '//trim(fn_avg)
       stop
    end if

    !// Header
    read(fnr) CINPUT_RUN
    !// Resolution (IPAR,JPAR,LPAR,NPAR)
    read(fnr) IN_RES
    !// Average info and version
    read(fnr) IN_TIME, version

    if (version .lt. 7) then
       write(6,'(a)') '* Cannot start from avgsav_input older than version 7'
       stop
    end if

    !// Check resolution
    if (IN_RES(1).ne.IPAR) then
       write(6,'(a,2i5)') '*** IPAR differ from input file',IN_RES(1), IPAR
       stop
    end if
    if (IN_RES(2).ne.JPAR) then
       write(6,'(a,2i5)') '*** JPAR differ from input file',IN_RES(2), JPAR
       stop
    end if
    if (IN_RES(3).ne.LPAR) then
       write(6,'(a,2i5)') '*** LPAR differ from input file',IN_RES(3), LPAR
       if (IN_RES(3).ne.LPARW) then
          write(6,'(a,3i5)') '*** LPAR not correct', IN_RES(3), LPAR, LPARW
       end if
       stop
    end if

    !// Allocate arrays to be read
    allocate(IN_IDMTC(IN_RES(4)))
    allocate(IN_TMASS(IN_RES(4)))
    allocate(IN_TNAME(IN_RES(4)))

    !// All the transported components
    read(fnr) IN_IDMTC
    !// Molecular masses of transported components (r8)
    read(fnr) IN_TMASS
    !// Tracer names
    read(fnr) IN_TNAME

    !// Non-transported
    read(fnr) IN_NOTRPAR
    if (IN_NOTRPAR .gt. 0) then
       allocate(IN_XIDMTC(IN_NOTRPAR))
       allocate(IN_XTMASS(IN_NOTRPAR))
       allocate(IN_XTNAME(IN_NOTRPAR))
       read(fnr) IN_XIDMTC
       read(fnr) IN_XTMASS
       read(fnr) IN_XTNAME
    end if

    !// DUMMY READING
    !// Grid area (r8)
    read(fnr) intmp1
    !// Longitude centers and edges (r8)
    read(fnr) intmp1
    !// Latitude centers and edges (r8)
    read(fnr) intmp1
    !// Hybrid sigma/pressure coordinates (r8)
    read(fnr) intmp1
    !// Average surface pressure (r4)
    read(fnr) intmp1
    !// Air (r4)
    read(fnr) intmp1
    !// Volume (r4)
    read(fnr) intmp1
    !// Air density (r4)
    read(fnr) intmp1
    !// Tempterature (r4)
    read(fnr) intmp1
    !// Height (r4)
    read(fnr) intmp1
    !// LMTROP
    read(fnr) intmp1
    !// H2O (r4) and Q (r4) if wanted
    read(fnr) inlog
    if (inlog) then
       read(fnr) intmp1 !// H2OAVG
       read(fnr) intmp1 !// QAVG
    end if

    write(6,'(2x,A2,1x,A6,4x,2x,A2,1x,A8,1x,A9,1x,A15)') &
         'NR','Tracer','ID','trsp_idx','Xtrsp_idx','TOTAL MASS (kg)'

    !// Component by component
    do N = 1, IN_RES(4)
       read(fnr) TRACER_ID
       read(fnr) RINDATA(:,:,:)
       !// Transported tracers are stored in STT, while
       !// non-transported are stored in XSTT.
       if (trsp_idx(TRACER_ID) .gt. 0) then
          TRSPNR = trsp_idx(TRACER_ID)
          STT(:,:,:,TRSPNR) = RINDATA(:,:,:)
          write(6,'(1x,i3,1x,a10,1x,i3,1x,i8,6x,es12.6)') &
               N,TNAME(TRSPNR),TRACER_ID,TRSPNR, sum(STT(:,:,:,TRSPNR))
          !// Non-zero initialization
          zeroinit(TRSPNR) = 0
       else if (Xtrsp_idx(TRACER_ID) .gt. 0) then
          !// Tracer is not transported in this simulation;
          !// put tracer in XSTT
          TRSPNR = Xtrsp_idx(TRACER_ID)
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XSTT(:,TRSPNR,I,J) = RINDATA(I,J,L)
                end do
             end do
          end do
          write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
               N,XTNAME(TRSPNR),TRACER_ID,TRSPNR, sum(XSTT(:,TRSPNR,:,:))
          !// Non-zero initialization
          zeroinit(NPAR+TRSPNR) = 0
       else
          !// Tracer is not defined
          write(6,'(a,i3,a)') 'Tracer number ',TRACER_ID, &
               ' is not defined; skipping!'
       end if
    end do
    deallocate(IN_IDMTC,IN_TMASS,IN_TNAME)


    !// Read non-transported data
    !// For reading no-transported data into STT, make a new read-in
    !// routine.
    if (IN_NOTRPAR .gt. 0) THEN
       !// Not checking for NOTRPAR>0 generates floating point exception
       !// during compiling; then crash.

       allocate(RXSTT(IN_RES(3),IN_NOTRPAR,IN_RES(1),IN_RES(2)))
       !// All components (LPAR,NOTRPAR,IPAR,JPAR)
       read(fnr) RXSTT

       do N = 1, IN_NOTRPAR

          TRACER_ID = IN_XIDMTC(N)
          !// Tracer is not transported, put it into XSTT
          TRSPNR = Xtrsp_idx(TRACER_ID)
          if (TRSPNR .gt. 0) then
             do J = 1, JPAR
                do I = 1, IPAR
                   XSTT(:,TRSPNR,I,J) = RXSTT(:,N,I,J)
                end do
             end do
             write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
                  N,XTNAME(TRSPNR),TRACER_ID,TRSPNR, sum(XSTT(:,TRSPNR,:,:))
             !// Non-zero initialization
             zeroinit(NPAR+TRSPNR) = 0
          else if (trsp_idx(TRACER_ID) .gt. 0) then
             TRSPNR = trsp_idx(TRACER_ID)
             do L = 1, LPAR
                do J = 1, JPAR
                   do I = 1, IPAR
                      STT(I,J,L,TRSPNR) = RXSTT(L,N,I,J)
                   end do
                end do
             end do
             write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
                  N,TNAME(TRSPNR),TRACER_ID,TRSPNR, sum(STT(:,:,:,TRSPNR))
             !// Non-zero initialization
             zeroinit(TRSPNR) = 0
          else
             !// Tracer is not defined
             write(6,'(a,i3,a)') 'Tracer number ',TRACER_ID, &
                  ' is not defined; skipping!'
          end if

       end do
       close(fnr)
       deallocate(RXSTT, IN_XIDMTC,IN_XTMASS,IN_XTNAME)

       !// NO3 to something very small
       if (Xtrsp_idx(41).gt.0) RXSTT(:,Xtrsp_idx(41),:,:) = 1.e-4_r8

    end if !// if (NOTRPAR .gt. 0)

    !// NO3 to something very small
    if (trsp_idx(41).gt.0) STT(:,:,:,trsp_idx(41)) = 1.e-4_r8
    !// Check for negative
    do N = 1, NPAR
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (STT(I,J,L,N).le.0._r8) STT(I,J,L,N) = 1.e-20_r8
             end do
          end do
       end do
    end do

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    !// --------------------------------------------------------------------
  end subroutine restart_from_CTM3avg
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  !// Next follows routines to restart from another horizontal resolution,
  !// currently set up to be T42.
  !// ----------------------------------------------------------------------
  subroutine restart_from_CTM3avg_T42(zeroinit)
    !// --------------------------------------------------------------------
    !// Restart from Oslo CTM3 avgsav file, version 7 and newer.
    !// Reads a CTM3 file (L37/40/57/60) and collapses it simply if
    !// necessary (when possible): only possible for L37.
    !//
    !// Ole Amund Sovde, November 2013
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR, LPARW
    use cmn_ctm, only: LCONT, STT, XDEDG, YDEDG
    use cmn_chem, only: TNAME
    use cmn_parameters, only: ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    use cmn_oslo, only: trsp_idx, XTNAME, Xtrsp_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Output
    integer, intent(out) :: zeroinit(NPAR+NOTRPAR)

    !// Locals
    integer :: &
         IN_RES(4), IN_TIME(3), &
         version, IN_NOTRPAR, I, J, L, N, &
         TRSPNR, ioerr, II, JJ, LL, TRACER_ID

    character(len=80) :: CINPUT_RUN

    real(r8) :: RFAC !// for unit conversion

    !// for file handling
    logical :: ex
    integer :: fnr, err
    character(len=80) :: FN_AVG

    real(r4) :: intmp1
    logical :: inlog

    !// Read-in data, possibly different resolution
    integer, allocatable :: IN_IDMTC(:), IN_XIDMTC(:)
    real(r4), allocatable :: RXSTT(:,:,:,:)
    real(r8), allocatable :: IN_TMASS(:), IN_XTMASS(:)
    character(len=10), allocatable :: IN_TNAME(:), IN_XTNAME(:)

    integer,parameter :: I42=128, J42=64
    real(r4) :: RINDATA(I42,J42,LPAR)
    real(r8) :: R8IN(I42,J42), R8OUT(IPAR,JPAR)

    real(r8) :: WGL42(J42),WGW42(J42), WYEDG1P1(J42+1)
    real(r8), dimension(J42+1) :: YDEDG42
    real(r8), dimension(I42+1) :: XDEDG42
    real(r8) :: DEL_ID
    real(r8),parameter :: GM0000 = 1.5_r8
    !// --------------------------------------------------------------------

    !// Will read from this file
    FN_AVG = 'avgsav_input'

    if (LCONT) then
       write(6,'(a)') 'Check LCONT!'
       write(6,'(a)') 'Trying to read '//trim(FN_AVG)//' when LCONT=.true.'
       stop
    end if

    inquire(File=trim(FN_AVG),exist=ex)
    if (.not. ex) then
       write(6,'(a)') '** The file '//trim(FN_AVG)//' does not exist **'
       stop
    end if

    fnr = get_free_fileid()

    write(6,'(/a)') 'Reading initial tracer data from: '//trim(FN_AVG)
    open(fnr,File=FN_AVG,Form='UNFORMATTED',Status='OLD',iostat=ioerr)
    if (ioerr .ne. 0) then
       write(6,'(a)') '*** error opening avg-file '//trim(fn_avg)
       stop
    end if

    !// Header
    read(fnr) CINPUT_RUN
    !// Resolution (IPAR,JPAR,LPAR,NPAR)
    read(fnr) IN_RES
    !// Average info and version
    read(fnr) IN_TIME, version

    if (version .lt. 7) then
       write(6,'(a)') '* Cannot start from avgsav_input older than version 7'
       stop
    end if

    !// Check resolution
    if (IN_RES(1).ne.I42) then
       write(6,'(a,2i5)') '*** IN_RES(1) /= 128 ',IN_RES(1), I42
       stop
    end if
    if (IN_RES(2).ne.J42) then
       write(6,'(a,2i5)') '*** IN_RES(1) /= 64 ',IN_RES(1), J42
       stop
    end if
    if (IN_RES(3).ne.LPAR) then
       write(6,'(a,2i5)') '*** IN_RES(3) /= LPAR ',IN_RES(3), LPAR
       stop
    end if

    !// Need XBEDGE and YBEDGE for T42
    DEL_ID = 360._r8 / real(IN_RES(1), r8)
    do I = 1, IN_RES(1)+1
       XDEDG42(I) = (real(I, r8) - GM0000) * DEL_ID
    end do
      
    call GAUSST2(IN_RES(2),WGL42,WGW42)
    WYEDG1P1(1) = -1._r8
    do J = 2, IN_RES(2)/2
       WYEDG1P1(J) = WYEDG1P1(J-1) + WGW42(J-1)
    end do
    WYEDG1P1(IN_RES(2)/2+1) = 0._r8
    do J = IN_RES(2)/2+2, IN_RES(2)+1
       WYEDG1P1(J) = -WYEDG1P1(IN_RES(2)+2-J)
    end do
    do J = 1, IN_RES(2)+1
       YDEDG42(J) = ZPI180 * asin(WYEDG1P1(J))
    end do


    !// Allocate arrays to be read
    allocate(IN_IDMTC(IN_RES(4)))
    allocate(IN_TMASS(IN_RES(4)))
    allocate(IN_TNAME(IN_RES(4)))

    !// All the transported components
    read(fnr) IN_IDMTC
    !// Molecular masses of transported components (r8)
    read(fnr) IN_TMASS
    !// Tracer names
    read(fnr) IN_TNAME

    !// Non-transported
    read(fnr) IN_NOTRPAR
    if (IN_NOTRPAR .gt. 0) then
       allocate(IN_XIDMTC(IN_NOTRPAR))
       allocate(IN_XTMASS(IN_NOTRPAR))
       allocate(IN_XTNAME(IN_NOTRPAR))
       read(fnr) IN_XIDMTC
       read(fnr) IN_XTMASS
       read(fnr) IN_XTNAME
    end if

    !// DUMMY READING
    !// Grid area (r8)
    read(fnr) intmp1
    !// Longitude centers and edges (r8)
    read(fnr) intmp1
    !// Latitude centers and edges (r8)
    read(fnr) intmp1
    !// Hybrid sigma/pressure coordinates (r8)
    read(fnr) intmp1
    !// Average surface pressure (r4)
    read(fnr) intmp1
    !// Air (r4)
    read(fnr) intmp1
    !// Volume (r4)
    read(fnr) intmp1
    !// Air density (r4)
    read(fnr) intmp1
    !// Tempterature (r4)
    read(fnr) intmp1
    !// Height (r4)
    read(fnr) intmp1
    !// LMTROP
    read(fnr) intmp1
    !// H2O (r4) and Q (r4) if wanted
    read(fnr) inlog
    if (inlog) then
       read(fnr) intmp1 !// H2OAVG
       read(fnr) intmp1 !// QAVG
    end if

    write(6,'(2x,A2,1x,A6,4x,2x,A2,1x,A8,1x,A9,1x,A15)') &
         'NR','Tracer','ID','trsp_idx','Xtrsp_idx','TOTAL MASS (kg)'

    !// Component by component
    do N = 1, IN_RES(4)
       read(fnr) TRACER_ID
       read(fnr) RINDATA(:,:,:)
       !// Transported tracers are stored in STT, while
       !// non-transported are stored in XSTT.
       if (trsp_idx(TRACER_ID) .gt. 0) then
          TRSPNR = trsp_idx(TRACER_ID)
!$omp parallel private (L,R8IN,R8OUT) &
!$omp          shared (STT,XDEDG42,YDEDG42,XDEDG,YDEDG,RINDATA, &
!$omp                  TRSPNR,IN_RES) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, IN_RES(2)
                do I = 1, IN_RES(1)
                   R8IN(I,J) = RINDATA(I,J,L)
                end do
             end do
             !// "Interpolate"
             call E_GRID(R8IN,XDEDG42,YDEDG42,I42,J42, &
                  R8OUT,XDEDG,YDEDG,IPAR,JPAR,1)
             do J = 1, JPAR
                do I = 1, IPAR
                   STT(I,J,L,TRSPNR) = R8OUT(I,J)
                end do
             end do
          end do
!$omp end do
!$omp end parallel
          write(6,'(1x,i3,1x,a10,1x,i3,1x,i8,6x,es12.6)') &
               N,TNAME(TRSPNR),TRACER_ID,TRSPNR, &
               sum(STT(:,:,:,TRSPNR))
          !// Non-zero initialization
          zeroinit(TRSPNR) = 0
       else if (Xtrsp_idx(TRACER_ID) .gt. 0) then
          !// Tracer is not transported in this simulation;
          !// put tracer in XSTT
          TRSPNR = Xtrsp_idx(TRACER_ID)
!$omp parallel private (L,R8IN,R8OUT) &
!$omp          shared (XSTT,XDEDG42,YDEDG42,XDEDG,YDEDG,RINDATA, &
!$omp                  TRSPNR,IN_RES) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, IN_RES(2)
                do I = 1, IN_RES(1)
                   R8IN(I,J) = RINDATA(I,J,L)
                end do
             end do
             !// "Interpolate"
             call E_GRID(R8IN,XDEDG42,YDEDG42,I42,J42, &
                  R8OUT,XDEDG,YDEDG,IPAR,JPAR,1)
             do J = 1, JPAR
                do I = 1, IPAR
                   XSTT(L,TRSPNR,I,J) = R8OUT(I,J)
                end do
             end do
          end do
!$omp end do
!$omp end parallel
          write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
               N,XTNAME(TRSPNR),TRACER_ID,TRSPNR, sum(XSTT(:,TRSPNR,:,:))
          !// Non-zero initialization
          zeroinit(NPAR+TRSPNR) = 0
       else
          !// Tracer is not defined
          write(6,'(a,i3,a)') 'Tracer number ',TRACER_ID, &
               ' is not defined; skipping!'
       end if
    end do
    deallocate(IN_IDMTC,IN_TMASS,IN_TNAME)


    !// Read non-transported data
    !// For reading no-transported data into STT, make a new read-in
    !// routine.
    if (IN_NOTRPAR .gt. 0) THEN
       !// Not checking for NOTRPAR>0 generates floating point exception
       !// during compiling; then crash.

       allocate(RXSTT(IN_RES(3),IN_NOTRPAR,IN_RES(1),IN_RES(2)))
       !// All components (LPAR,NOTRPAR,IPAR,JPAR)
       read(fnr) RXSTT

       do N = 1, IN_NOTRPAR

          TRACER_ID = IN_XIDMTC(N)
          !// Tracer is not transported, put it into XSTT
          TRSPNR = Xtrsp_idx(TRACER_ID)
          if (TRSPNR .gt. 0) then
!$omp parallel private (L,R8IN,R8OUT) &
!$omp          shared (XSTT,XDEDG42,YDEDG42,XDEDG,YDEDG,RXSTT, &
!$omp                  TRSPNR,IN_RES,N) &
!$omp          default(NONE)
!$omp do
             do L = 1, LPAR
                do J = 1, IN_RES(2)
                   do I = 1, IN_RES(1)
                      R8IN(I,J) = RXSTT(L,N,I,J)
                   end do
                end do
                !// "Interpolate"
                call E_GRID(R8IN,XDEDG42,YDEDG42,I42,J42, &
                     R8OUT,XDEDG,YDEDG,IPAR,JPAR,1)
                do J = 1, JPAR
                   do I = 1, IPAR
                      XSTT(L,TRSPNR,I,J) = R8OUT(I,J)
                   end do
                end do
             end do
!$omp end do
!$omp end parallel
             write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
                  N,XTNAME(TRSPNR),TRACER_ID,TRSPNR, sum(XSTT(:,TRSPNR,:,:))
             !// Non-zero initialization
             zeroinit(NPAR+TRSPNR) = 0
          else if (trsp_idx(TRACER_ID) .gt. 0) then
             TRSPNR = trsp_idx(TRACER_ID)
!$omp parallel private (L,R8IN,R8OUT) &
!$omp          shared (STT,XDEDG42,YDEDG42,XDEDG,YDEDG,RXSTT, &
!$omp                  TRSPNR,IN_RES,N) &
!$omp          default(NONE)
!$omp do
             do L = 1, LPAR
                do J = 1, IN_RES(2)
                   do I = 1, IN_RES(1)
                      R8IN(I,J) = RXSTT(L,N,I,J)
                   end do
                end do
                !// "Interpolate"
                call E_GRID(R8IN,XDEDG42,YDEDG42,I42,J42, &
                     R8OUT,XDEDG,YDEDG,IPAR,JPAR,1)
                do J = 1, JPAR
                   do I = 1, IPAR
                      STT(I,J,L,TRSPNR) = R8OUT(I,J)
                   end do
                end do
             end do
!$omp end do
!$omp end parallel
             write(6,'(1x,i3,1x,a10,1x,i3,6x,i8,1x,es12.6)') &
                  N,TNAME(TRSPNR),TRACER_ID,TRSPNR, sum(STT(:,:,:,TRSPNR))
             !// Non-zero initialization
             zeroinit(TRSPNR) = 0
          else
             !// Tracer is not defined
             write(6,'(a,i3,a)') 'Tracer number ',TRACER_ID, &
                  ' is not defined; skipping!'
          end if

       end do
       close(fnr)
       deallocate(RXSTT, IN_XIDMTC,IN_XTMASS,IN_XTNAME)

       !// NO3 to something very small
       if (Xtrsp_idx(41).gt.0) RXSTT(:,Xtrsp_idx(41),:,:) = 1.e-4_r8

    end if !// if (NOTRPAR .gt. 0)

    !// NO3 to something very small
    if (trsp_idx(41).gt.0) STT(:,:,:,trsp_idx(41)) = 1.e-4_r8
    !// Check for negative
    do N = 1, NPAR
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (STT(I,J,L,N).le.0._r8) STT(I,J,L,N) = 1.e-20_r8
             end do
          end do
       end do
    end do

    write(6,'(a71)') '--------------------------------------------'// &
          '---------------------------'
    !// ------------------------------------------------------------------
  end subroutine restart_from_CTM3avg_T42
  !// --------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine OSLO_CON_RUN42(FN_CON, MODE)
    !// --------------------------------------------------------------------
    !// Read restart file and set up correct day information.
    !// Reads T42 and interpolates. Interpolates tracer mass and vertical
    !// moments (TW and WW), while zeroing the other moments.
    !//
    !// Can read ADDITIONAL restart files to fill in non-initialized
    !// tracers. This is controlled by MODE:
    !// MODE:  0: First call; initializes zeroinit.
    !//      !=0: Will add non-initialized tracers that are
    !//           available on the specified file. Skips AIR.
    !//
    !// Ole Amund Sovde, December 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_ctm, only: IYEAR, IDAY, JYEAR, JDAY, JMON, JDATE, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         XDEDG, YDEDG, NTM
    use cmn_chem, only: TNAME, TMASS
    use cmn_parameters, only: ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    use cmn_oslo, only: ZEROINIT, XZEROINIT, trsp_idx, chem_idx, &
         XTNAME, XTMASS, Xtrsp_idx, Xchem_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in):: FN_CON
    integer, intent(in) :: MODE

    !// Locals
    integer,parameter :: I42=128,J42=64
    real(r4) :: R4XYZ(I42,J42,LPAR)
    real(r8) :: R8XY(I42,J42),R8XYZ(I42,J42,LPAR)
    real(r8) :: C3XY(IPAR,JPAR), C3XYZ(IPAR,JPAR,LPAR)

    real(r8)  :: SUMT
    real(r8),dimension(I42,J42,LPAR)  :: INFIELD,INAIR
    real(r8),dimension(LPAR,I42,J42) :: XINFIELD
    real(rMom),dimension(I42,J42,LPAR) :: &
         TSUT, TSVT, TSWT, &
         TSUU, TSVV, TSWW, &
         TSUV, TSUW, TSVW
    real(r8), dimension(I42+1) :: INXDEDG
    real(r8), dimension(J42+1) :: INYDEDG
    real(r8), dimension(LPAR+1) :: INETAA, INETAB

    real(r8)  WGL42(J42),WGW42(J42), WYEDG1P1(J42+1)
    real(r8) :: DEL_ID
    real(r8),parameter :: GM0000=1.5_r8

    integer :: IIYEAR,IIDAY, JJYEAR,JJDAY,JJMON,JJDATE, version
    integer :: I,J,L,N,IDW,JDW,LDW, ID,JD,LD
    integer :: ioerr, file_id, ND, XND, CID
    logical :: LERR1, LERR2, LERR3
    character(len=80) :: RTITLE2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CON_RUN42'
    !// --------------------------------------------------------------------

    !// Initialize arrays; assume all are zero
    if (MODE .eq. 0) then
       ZEROINIT(:)  = 1
       if (NOTRPAR .gt. 0) XZEROINIT(:) = 1
    end if

    !// WARNING
    write(6,*)
    do I = 1, 3
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
    end do
    write(6,'(a)') '   Using OSLO_CON_RUN42!'
    write(6,'(a)') '   Restarting from a T42 restart.sav!'
    write(6,'(a)') '   Only STT and vertical moments are interpolated!'
    do I = 1, 3
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
    end do
    write(6,*)

    !//read in continuation restart file:
    if (MODE .eq. 0) then
       write(6,'(a)') '* Read restart file: '//trim(FN_CON)
    else
       write(6,'(a)') '* Read ADDITIONAL restart file: '//trim(FN_CON)
    end if

    !// Find a file number
    file_id = get_free_fileid()

    open(file_id,file=FN_CON,status='OLD',form='UNFORMATTED', &
         iostat=ioerr,action='read')
    if (ioerr .ne. 0) then
       write(6,'(a)') ' >>>>error opening continuation file: ',trim(FN_CON)
       stop
    end if

    read(file_id) RTITLE2
    write(6,'(2a)') 'CTM run: ',trim(RTITLE2)

    !// Read dimensions and transported tracers
    read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE,version
    if (ioerr .ne. 0) then
       !// Version is not specified, meaning version 2
       version = 2
       backspace(file_id)
       read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    end if
    write(6,'(a,6i6)') 'Model date:',IYEAR,IDAY,JYEAR,JDAY,JMON,JDATE
    write(6,'(a,6i6)') 'File date :',IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    write(6,'(a,i6)')  'File version :',version

    !// It is possible to override the date information here
    !IYEAR = IIYEAR
    !IDAY  = IIDAY
    !JYEAR = JJYEAR
    !JDAY  = JJDAY
    !JMON  = JJMON
    !JDATE = JJDATE

    if (version .le. 2) then
       read(file_id,iostat=ioerr) IDW,JDW,LDW,I,ID,JD,LD,ND,XND
    else
       !// Read native resolution and model resolution
       read(file_id,iostat=ioerr) IDW,JDW,LDW, ID,JD,LD, ND,XND
    end if

    if (ioerr .ne. 0) then
       write(6,'(a,8(1x,i3))') f90file//':'//subr// &
            ': Problems reading IDW,JDW,LDW,ID,JD,LD,ND,XND', &
            IDW,JDW,LDW, ID,JD,LD, ND,XND
       stop
    end if

    !// Check dimensions: vertical resolution (LD vs LPAR)
    if (LD .ne. LPAR) then
       write(6,'(a,2(1x,i3))') f90file//':'//subr// &
            ': Wrong vertical resolution of restart file '//trim(FN_CON),&
            LD,LPAR
       stop
    end if

    !// Do not use this routine if ID.ne.I42: that will end up bad
    if (ID .ne. I42) then
       write(6,'(a,4(1x,i3))') f90file//':'//subr// &
            ': File resolution is not T42: '//trim(FN_CON),ID,I42,JD,J42
       stop
    end if

    write(6,'(a,i4)') 'Transported species on file :    ',ND
    write(6,'(a,i4)') 'Non-transported species on file: ',XND

    !// Read grid box info for version > 2
    if (version .gt. 2) then
       read(file_id,iostat=ioerr) INXDEDG, INYDEDG, INETAA, INETAB
       if (ioerr .ne. 0) then
          write(6,'(a)') f90file//':'//subr// &
               ': Problems reading XDEDG/YDEDG/ETAA/ETAB'
          stop
       end if
    else
       !// Calculate XBEDGE and YBEDGE for T42
       DEL_ID = 360._r8 / real(I42, r8)
       do I = 1, I42+1
          INXDEDG(I) = (real(I, r8) - GM0000) * DEL_ID
       end do
      
       call GAUSST2(J42,WGL42,WGW42)
       WYEDG1P1(1) = -1._r8
       do J = 2, J42/2
          WYEDG1P1(J) = WYEDG1P1(J-1) + WGW42(J-1)
       end do
       WYEDG1P1(J42/2+1) = 0._r8
       do J = J42/2+2, J42+1
          WYEDG1P1(J) = -WYEDG1P1(J42+2-J)
       end do
       do J = 1, J42+1
          INYDEDG(J) = ZPI180 * asin(WYEDG1P1(J))
       end do
    end if !// if (version .gt. 2) then

    !// Read air
    read(file_id,iostat=ioerr) INAIR
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': Problems reading AIR'
       stop
    end if

    if (MODE .eq. 0) then
       !// Interpolate horizontally
!$omp parallel private (L) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INAIR,AIR) &
!$omp          default(NONE)
!$omp do
       do L = 1, LPAR
          !// Interpolate
          call E_GRID(INAIR(:,:,L),INXDEDG,INYDEDG,I42,J42, &
               AIR(:,:,L),XDEDG,YDEDG,IPAR,JPAR,1)
       end do
!$omp end do
!$omp end parallel
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (AIR(I,J,L).ne.AIR(I,J,L)) then
                   write(6,'(a,4i5)') f90file//':'//subr// &
                        ': AIR NAN',I,J,L
                   stop
                end if
             end do
          end do
       end do
       print*
       write(6,'(a,2es20.12)') '* AIR read and set!', sum(air),sum(inair)
    else
       write(6,'(a)') '* AIR read and skipped!'
    end if


    !// Read each of the transported components
    do N = 1, ND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading CID',CID
          stop
       end if
       !// STT
       read(file_id) INFIELD
       read(file_id) TSUT(:,:,:)
       read(file_id) TSVT(:,:,:)
       read(file_id) TSWT(:,:,:)
       read(file_id) TSUU(:,:,:)
       read(file_id) TSVV(:,:,:)
       read(file_id) TSWW(:,:,:)
       read(file_id) TSUV(:,:,:)
       read(file_id) TSUW(:,:,:)
       read(file_id) TSVW(:,:,:)


       !// Place component if included in the run
       if (trsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (ZEROINIT(trsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          !// Interpolate horizontally
          !//STT(:,:,:,trsp_idx(CID)) = INFIELD(:,:,:)
!$omp parallel private (L,C3XY,R8XY) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INFIELD,TSWT,TSWW, &
!$omp                  STT,SWT,SWW,trsp_idx,CID) &
!$omp          default(NONE)
!$omp do
          do L=1,LPAR
             !// Interpolate
             call E_GRID(INFIELD(:,:,L),INXDEDG,INYDEDG,I42,J42, &
                  STT(:,:,L,trsp_idx(CID)),XDEDG,YDEDG,IPAR,JPAR,1)

             !// Interpolate vertical moments only (moments are real*4)
             R8XY(:,:) = real(TSWT(:,:,L), r8)
             call E_GRID(R8XY,INXDEDG,INYDEDG,I42,J42, &
                  C3XY,XDEDG,YDEDG,IPAR,JPAR,1)
             SWT(:,:,L,trsp_idx(CID)) = real(C3XY(:,:), rMom)

             R8XY(:,:) = real(TSWW(:,:,L), r8)
             call E_GRID(R8XY,INXDEDG,INYDEDG,I42,J42, &
                  C3XY,XDEDG,YDEDG,IPAR,JPAR,1)
             SWW(:,:,L,trsp_idx(CID)) = real(C3XY(:,:), rMom)
          end do
!$omp end do
!$omp end parallel
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (transported)'
          !// Keep track of initialized components
          ZEROINIT(trsp_idx(CID)) = 0
       else if (Xtrsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (XZEROINIT(Xtrsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          !// Put transported component into non-transported array
!$omp parallel private (L) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INFIELD,C3XYZ) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             !// Interpolate
             call E_GRID(INFIELD(:,:,L),INXDEDG,INYDEDG,I42,J42, &
                  C3XYZ(:,:,L),XDEDG,YDEDG,IPAR,JPAR,1)
          end do
!$omp end do
!$omp end parallel
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XSTT(L,trsp_idx(CID),I,J) = C3XYZ(I,J,L)
                end do
             end do
          end do
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (not transported)'
          !// Keep track of initialized components
          XZEROINIT(Xtrsp_idx(CID)) = 0
       else
          write(6,'(a,i3,a)') '  Skipping component ',CID,' (not included)'
       end if

    end do !//  do N = 1, ND

    do N = 1, NPAR
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (STT(I,J,L,N).ne.STT(I,J,L,N)) then
                   write(6,'(a,4i5)') f90file//':'//subr// &
                        ': STT NAN',I,J,L,N
                   stop
                end if
             end do
          end do
       end do
    end do



    !// Read each of the non-transported components
    do N = 1, XND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading X CID',CID
          stop
       end if
       !// XSTT
       read(file_id) XINFIELD

       !// Place component if included in the run
       if (Xtrsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (XZEROINIT(Xtrsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ',CID,' (not transported)'
          !// XSTT(:,Xtrsp_idx(CID),:,:) = XINFIELD(:,:,:)
!$omp parallel private (L,I,J,R8XY,C3XY) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG) &
!$omp          shared (XINFIELD,XSTT,Xtrsp_idx,CID) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, J42
                do I = 1, I42
                   R8XY(I,J) = XINFIELD(L,I,J)
                end do
             end do
             !// Interpolate
             call E_GRID(R8XY(:,:),INXDEDG,INYDEDG,I42,J42, &
                  C3XY(:,:),XDEDG,YDEDG,IPAR,JPAR,1)

             do J = 1, JPAR
                do I = 1, IPAR
                   XSTT(L,Xtrsp_idx(CID),I,J) = C3XY(I,J)
                end do
             end do
          end do
!$omp end do
!$omp end parallel

          !// Keep track of initialized components
          XZEROINIT(Xtrsp_idx(CID)) = 0
       else if (trsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (ZEROINIT(trsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ', &
               CID,' (transported, zeroed moments)'
!$omp parallel private (L,I,J,R8XY) &
!$omp          shared (XINFIELD,INXDEDG,INYDEDG,XDEDG,YDEDG,STT,trsp_idx,CID) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, J42
                do I = 1, I42
                   R8XY(I,J) = XINFIELD(L,I,J)
                end do
             end do
             !// Interpolate
             call E_GRID(R8XY(:,:),INXDEDG,INYDEDG,I42,J42, &
                  STT(:,:,L,trsp_idx(CID)),XDEDG,YDEDG,IPAR,JPAR,1)
          end do
!$omp end do
!$omp end parallel
          !// Moments are set to zero
          !// Keep track of initialized components
          ZEROINIT(trsp_idx(CID)) = 0
       else
          !// Ah, tracer is not included in this run
          write(6,'(a,i3,a)') '  Skipping X component ',CID,' (not included)'
       end if

    end do !//  do N = 1, XND

    close(file_id)
    write(6,'(a)') '* Finished reading restart file: '//trim(FN_CON)


    !// Print out info about the initialization field
    !// --------------------------------------------------------------------
    write(6,'(A,7I6)') 'JYEAR/JDATE/IDAY/IYEAR/JDAY/JMON', &
         JYEAR,JDATE,IDAY,IYEAR,JDAY,JMON
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// Transported species
    write(6,'(a)') 'Transported species:'
    do N = 1, NTM
       !// Print out only if initialized
       if (ZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                SUMT = SUMT + STT(I,J,L,N)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,TNAME(N),TMASS(N),SUMT
    end do

    !// Non-transported species
    write(6,'(a)') 'Non-transported species:'
    do N = 1, NOTRPAR
       !// Print out only if initialized
       if (XZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                SUMT = SUMT + XSTT(L,N,I,J)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,XTNAME(N),XTMASS(N),SUMT
    end do

    !// Which components are initialized as zero?
    do N = 1, NTM
       if (ZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** Component trsp_idx/chem_idx ', &
               N,chem_idx(N),' is initialized as ZERO!'
       end if
    end do
    do N = 1, NOTRPAR
       if (XZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** X-Component Xtrsp_idx/Xchem_idx ', &
               N,Xchem_idx(N),' is initialized as ZERO!'
       end if
    end do


    !// Info if your restart file conains more/less species
    ND = sum(ZEROINIT)
    if (NOTRPAR .gt. 0) then
       XND = sum(XZEROINIT)
    else
       XND = 0
    endif
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    if ((ND+XND) .gt. 0) then
       write(6,'(a,i3)') '  Transported species not initialized: ', ND
       write(6,'(a,i3)') '  Non-Transported species not initialized: ', XND
       write(6,'(a)') 'Program _may_ stop if critical components'// &
            ' are not initialized!'
    else
       write(6,'(a)') '  All species are initialized!'
    end if
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine OSLO_CON_RUN42
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine OSLO_CON_RUNxx(FN_CON, MODE)
    !// --------------------------------------------------------------------
    !// Read restart file and set up correct day information.
    !// Reads any horizontal resolution sav-file, as long as it is
    !// version 3 or higher.
    !// Interpolates tracer mass and vertical
    !// moments (TW and WW), while zeroing the other moments.
    !//
    !// Can read ADDITIONAL restart files to fill in non-initialized
    !// tracers. This is controlled by MODE:
    !// MODE:  0: First call; initializes zeroinit.
    !//      !=0: Will add non-initialized tracers that are
    !//           available on the specified file. Skips AIR.
    !//
    !// Ole Amund Sovde, November 2015
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, r4, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_ctm, only: IYEAR, IDAY, JYEAR, JDAY, JMON, JDATE, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         XDEDG, YDEDG, NTM
    use cmn_chem, only: TNAME, TMASS
    use cmn_parameters, only: ZPI180
    use grid, only: GAUSST2
    use regridding, only: E_GRID
    use utilities, only: get_free_fileid
    use cmn_oslo, only: ZEROINIT, XZEROINIT, trsp_idx, chem_idx, &
         XTNAME, XTMASS, Xtrsp_idx, Xchem_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in):: FN_CON
    integer, intent(in) :: MODE

    !// Locals
    integer :: INIPAR, INJPAR
    real(r8) :: C3XY(IPAR,JPAR), C3XYZ(IPAR,JPAR,LPAR)
    real(r8) :: SUMT
    real(r8), dimension(LPAR+1) :: INETAA, INETAB

    integer :: IIYEAR,IIDAY, JJYEAR,JJDAY,JJMON,JJDATE, version
    integer :: I,J,L,N,IDW,JDW,LDW, ID,JD,LD
    integer :: ioerr, file_id, ND, XND, CID
    logical :: LERR1, LERR2, LERR3
    character(len=80) :: RTITLE2

    !// Allocatable arrays to be read from file
    real(r4), dimension(:,:,:), allocatable :: R4XYZ
    real(r8), dimension(:,:), allocatable :: R8XY
    real(r8), dimension(:,:,:), allocatable :: R8XYZ
    real(r8), dimension(:,:,:), allocatable :: INFIELD,INAIR
    real(r8), dimension(:,:,:), allocatable :: XINFIELD
    real(rMom), dimension(:,:,:), allocatable :: &
         TSUT, TSVT, TSWT, &
         TSUU, TSVV, TSWW, &
         TSUV, TSUW, TSVW
    real(r8), dimension(:), allocatable :: INXDEDG
    real(r8), dimension(:), allocatable :: INYDEDG
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_CON_RUNxx'
    !// --------------------------------------------------------------------

    !// Initialize arrays; assume all are zero
    if (MODE .eq. 0) then
       ZEROINIT(:)  = 1
       if (NOTRPAR .gt. 0) XZEROINIT(:) = 1
    end if

    !// WARNING
    write(6,*)
    do I = 1, 3
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
    end do
    write(6,'(a)') '   Using '//subr//'!'
    write(6,'(a)') '   Restarting from a different horizontal resolution!'
    write(6,'(a)') '   Only STT and vertical moments are interpolated!'
    do I = 1, 3
       write(6,'(a)') 'WARNING WARNING WARNING WARNING WARNING WARNING'
    end do
    write(6,*)

    !//read in continuation restart file:
    if (MODE .eq. 0) then
       write(6,'(a)') '* Read restart file: '//trim(FN_CON)
    else
       write(6,'(a)') '* Read ADDITIONAL restart file: '//trim(FN_CON)
    end if

    !// Find a file number
    file_id = get_free_fileid()

    open(file_id,file=FN_CON,status='OLD',form='UNFORMATTED', &
         iostat=ioerr,action='read')
    if (ioerr .ne. 0) then
       write(6,'(a)') ' >>>>error opening continuation file: '//trim(FN_CON)
       stop
    end if

    read(file_id) RTITLE2
    write(6,'(a)') 'CTM run: '//trim(RTITLE2)

    !// Read dimensions and transported tracers
    read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE,version
    if (ioerr .ne. 0) then
       !// Cannot use routine for version < 3
       write(6,'(a,8(1x,i3))') f90file//':'//subr// &
            ': Cannot use routine for version < 3.'
       stop
    end if
    write(6,'(a,6i6)') 'Model date:',IYEAR,IDAY,JYEAR,JDAY,JMON,JDATE
    write(6,'(a,6i6)') 'File date :',IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    write(6,'(a,i6)')  'File version :',version

    !// It is possible to override the date information here
    !IYEAR = IIYEAR
    !IDAY  = IIDAY
    !JYEAR = JJYEAR
    !JDAY  = JJDAY
    !JMON  = JJMON
    !JDATE = JJDATE

    !// Read native resolution and model resolution
    read(file_id,iostat=ioerr) IDW,JDW,LDW, ID,JD,LD, ND,XND
    if (ioerr .ne. 0) then
       write(6,'(a,8(1x,i3))') f90file//':'//subr// &
            ': Problems reading IDW,JDW,LDW,ID,JD,LD,ND,XND', &
            IDW,JDW,LDW, ID,JD,LD, ND,XND
       stop
    end if

    !// Check dimensions: vertical resolution (LD vs LPAR)
    if (LD .ne. LPAR) then
       write(6,'(a,2(1x,i3))') f90file//':'//subr// &
            ': Wrong vertical resolution of restart file '//trim(FN_CON),&
            LD,LPAR
       stop
    end if

    !// Do not use this routine if ID.eq.IPAR
    if (ID .eq. IPAR) then
       write(6,'(a)') f90file//':'//subr// &
            ': File resolution is the same as model resolution!'
       write(6,'(a)') ' You should use OSLO_CON_RUN instead!'
       write(6,'(a)') ' Tried to read file: '//trim(FN_CON)
       stop
    end if

    write(6,'(a,i4)') 'Transported species on file :    ',ND
    write(6,'(a,i4)') 'Non-transported species on file: ',XND

    !// Resolution is ID/JD/LD
    !// Allocate arrays
    allocate( R4XYZ(ID,JD,LPAR), R8XY(ID,JD), R8XYZ(ID,JD,LPAR), &
              INFIELD(ID,JD,LPAR), INAIR(ID,JD,LPAR), XINFIELD(LPAR,ID,JD), &
              TSUT(ID,JD,LPAR), TSVT(ID,JD,LPAR), TSWT(ID,JD,LPAR), &
              TSUU(ID,JD,LPAR), TSVV(ID,JD,LPAR), TSWW(ID,JD,LPAR), &
              TSUV(ID,JD,LPAR), TSUW(ID,JD,LPAR), TSVW(ID,JD,LPAR), &
              INXDEDG(ID+1), INYDEDG(JD+1) )


    read(file_id,iostat=ioerr) INXDEDG, INYDEDG, INETAA, INETAB
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Problems reading XDEDG/YDEDG/ETAA/ETAB'
       stop
    end if

    !// Read air
    read(file_id,iostat=ioerr) INAIR
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr//': Problems reading AIR'
       stop
    end if

    if (MODE .eq. 0) then
       !// Interpolate horizontally
!$omp parallel private (L) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INAIR,AIR,ID,JD) &
!$omp          default(NONE)
!$omp do
       do L = 1, LPAR
          !// Interpolate
          call E_GRID(INAIR(:,:,L),INXDEDG,INYDEDG,ID,JD, &
               AIR(:,:,L),XDEDG,YDEDG,IPAR,JPAR,1)
       end do
!$omp end do
!$omp end parallel
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (AIR(I,J,L).ne.AIR(I,J,L)) then
                   write(6,'(a,4i5)') f90file//':'//subr// &
                        ': AIR NAN',I,J,L
                   stop
                end if
             end do
          end do
       end do
       print*
       write(6,'(a,2es20.12)') '* AIR read and set!', sum(air),sum(inair)
    else
       write(6,'(a)') '* AIR read and skipped!'
    end if


    !// Read each of the transported components
    do N = 1, ND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading CID',CID
          stop
       end if
       !// STT
       read(file_id) INFIELD
       read(file_id) TSUT(:,:,:)
       read(file_id) TSVT(:,:,:)
       read(file_id) TSWT(:,:,:)
       read(file_id) TSUU(:,:,:)
       read(file_id) TSVV(:,:,:)
       read(file_id) TSWW(:,:,:)
       read(file_id) TSUV(:,:,:)
       read(file_id) TSUW(:,:,:)
       read(file_id) TSVW(:,:,:)


       !// Place component if included in the run
       if (trsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (ZEROINIT(trsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          !// Interpolate horizontally
          !//STT(:,:,:,trsp_idx(CID)) = INFIELD(:,:,:)
!$omp parallel private (L,C3XY,R8XY) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INFIELD,TSWT,TSWW, &
!$omp                  STT,SWT,SWW,trsp_idx,CID,ID,JD) &
!$omp          default(NONE)
!$omp do
          do L=1,LPAR
             !// Interpolate
             call E_GRID(INFIELD(:,:,L),INXDEDG,INYDEDG,ID,JD, &
                  STT(:,:,L,trsp_idx(CID)),XDEDG,YDEDG,IPAR,JPAR,1)

             !// Interpolate vertical moments only (moments are real*4)
             R8XY(:,:) = real(TSWT(:,:,L), r8)
             call E_GRID(R8XY,INXDEDG,INYDEDG,ID,JD, &
                  C3XY,XDEDG,YDEDG,IPAR,JPAR,1)
             SWT(:,:,L,trsp_idx(CID)) = real(C3XY(:,:), rMom)

             R8XY(:,:) = real(TSWW(:,:,L), r8)
             call E_GRID(R8XY,INXDEDG,INYDEDG,ID,JD, &
                  C3XY,XDEDG,YDEDG,IPAR,JPAR,1)
             SWW(:,:,L,trsp_idx(CID)) = real(C3XY(:,:), rMom)
          end do
!$omp end do
!$omp end parallel
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (transported)'
          !// Keep track of initialized components
          ZEROINIT(trsp_idx(CID)) = 0
       else if (Xtrsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (XZEROINIT(Xtrsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          !// Put transported component into non-transported array
!$omp parallel private (L) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG,INFIELD,C3XYZ,ID,JD) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             !// Interpolate
             call E_GRID(INFIELD(:,:,L),INXDEDG,INYDEDG,ID,JD, &
                  C3XYZ(:,:,L),XDEDG,YDEDG,IPAR,JPAR,1)
          end do
!$omp end do
!$omp end parallel
          do J = 1, JPAR
             do I = 1, IPAR
                do L = 1, LPAR
                   XSTT(L,Xtrsp_idx(CID),I,J) = C3XYZ(I,J,L)
                end do
             end do
          end do
          write(6,'(a,i3,a)') '  Placing component  ',CID,' (not transported)'
          !// Keep track of initialized components
          XZEROINIT(Xtrsp_idx(CID)) = 0
       else
          write(6,'(a,i3,a)') '  Skipping component ',CID,' (not included)'
       end if

    end do !//  do N = 1, ND

    do N = 1, NPAR
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                if (STT(I,J,L,N).ne.STT(I,J,L,N)) then
                   write(6,'(a,4i5)') f90file//':'//subr// &
                        ': STT NAN',I,J,L,N
                   stop
                end if
             end do
          end do
       end do
    end do



    !// Read each of the non-transported components
    do N = 1, XND
       !// Component number
       read(file_id,iostat=ioerr) CID
       if (ioerr .ne. 0) then
          write(6,'(a,i5)') f90file//':'//subr// &
               ': Problems reading X CID',CID
          stop
       end if
       !// XSTT
       read(file_id) XINFIELD

       !// Place component if included in the run
       if (Xtrsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (XZEROINIT(Xtrsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ',CID,' (not transported)'
          !// XSTT(:,Xtrsp_idx(CID),:,:) = XINFIELD(:,:,:)
!$omp parallel private (L,I,J,R8XY,C3XY) &
!$omp          shared (INXDEDG,INYDEDG,XDEDG,YDEDG) &
!$omp          shared (XINFIELD,XSTT,Xtrsp_idx,CID,ID,JD) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, JD
                do I = 1, ID
                   R8XY(I,J) = XINFIELD(L,I,J)
                end do
             end do
             !// Interpolate
             call E_GRID(R8XY(:,:),INXDEDG,INYDEDG,ID,JD, &
                  C3XY(:,:),XDEDG,YDEDG,IPAR,JPAR,1)

             do J = 1, JPAR
                do I = 1, IPAR
                   XSTT(L,Xtrsp_idx(CID),I,J) = C3XY(I,J)
                end do
             end do
          end do
!$omp end do
!$omp end parallel

          !// Keep track of initialized components
          XZEROINIT(Xtrsp_idx(CID)) = 0
       else if (trsp_idx(CID).gt.0) then
          !// Check if already initialized
          if (ZEROINIT(trsp_idx(CID)) .ne. 1) then
             write(6,'(a,i3,a)') '  Component ',CID,' is already initialized'
             cycle
          end if
          write(6,'(a,i3,a)') '  Placing X component  ', &
               CID,' (transported, zeroed moments)'
!$omp parallel private (L,I,J,R8XY) &
!$omp          shared (XINFIELD,INXDEDG,INYDEDG,XDEDG,YDEDG,STT) &
!$omp          shared (trsp_idx,CID,ID,JD) &
!$omp          default(NONE)
!$omp do
          do L = 1, LPAR
             do J = 1, JPAR
                do I = 1, IPAR
                   R8XY(I,J) = XINFIELD(L,I,J)
                end do
             end do
             !// Interpolate
             call E_GRID(R8XY(:,:),INXDEDG,INYDEDG,ID,JD, &
                  STT(:,:,L,trsp_idx(CID)),XDEDG,YDEDG,IPAR,JPAR,1)
          end do
!$omp end do
!$omp end parallel
          !// Moments are set to zero
          !// Keep track of initialized components
          ZEROINIT(trsp_idx(CID)) = 0
       else
          !// Ah, tracer is not included in this run
          write(6,'(a,i3,a)') '  Skipping X component ',CID,' (not included)'
       end if

    end do !//  do N = 1, XND

    close(file_id)
    write(6,'(a)') '* Finished reading restart file: '//trim(FN_CON)

    !// Deallocate arrays
    deallocate( R4XYZ, R8XY, R8XYZ, INFIELD,INAIR,XINFIELD, &
                TSUT, TSVT, TSWT, TSUU, TSVV, TSWW, &
                TSUV, TSUW, TSVW, INXDEDG, INYDEDG )

    !// Print out info about the initialization field
    !// --------------------------------------------------------------------
    write(6,'(A,7I6)') 'JYEAR/JDATE/IDAY/IYEAR/JDAY/JMON', &
         JYEAR,JDATE,IDAY,IYEAR,JDAY,JMON
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// Transported species
    write(6,'(a)') 'Transported species:'
    do N = 1, NTM
       !// Print out only if initialized
       if (ZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do L = 1, LPAR
          do J = 1, JPAR
             do I = 1, IPAR
                SUMT = SUMT + STT(I,J,L,N)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,TNAME(N),TMASS(N),SUMT
    end do

    !// Non-transported species
    write(6,'(a)') 'Non-transported species:'
    do N = 1, NOTRPAR
       !// Print out only if initialized
       if (XZEROINIT(N) .eq. 1) cycle

       SUMT   = 0._r8
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, LPAR
                SUMT = SUMT + XSTT(L,N,I,J)
             end do
          end do
       end do
       write(6,'(a,i4,1x,a10,0P,f8.2,1P,es14.6)') &
            'tracer',N,XTNAME(N),XTMASS(N),SUMT
    end do

    !// Which components are initialized as zero?
    do N = 1, NTM
       if (ZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** Component trsp_idx/chem_idx ', &
               N,chem_idx(N),' is initialized as ZERO!'
       end if
    end do
    do N = 1, NOTRPAR
       if (XZEROINIT(N) .eq. 1) then
          write(6,'(a,i3,1x,i3,a)') '*** X-Component Xtrsp_idx/Xchem_idx ', &
               N,Xchem_idx(N),' is initialized as ZERO!'
       end if
    end do


    !// Info if your restart file conains more/less species
    ND = sum(ZEROINIT)
    if (NOTRPAR .gt. 0) then
       XND = sum(XZEROINIT)
    else
       XND = 0
    endif
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    if ((ND+XND) .gt. 0) then
       write(6,'(a,i3)') '  Transported species not initialized: ', ND
       write(6,'(a,i3)') '  Non-Transported species not initialized: ', XND
       write(6,'(a)') 'Program _may_ stop if critical components'// &
            ' are not initialized!'
    else
       write(6,'(a)') '  All species are initialized!'
    end if
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine OSLO_CON_RUNxx
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine OSLO_RESTARTFILE_INFO(FN_CON, VERSION, &
       IDW,JDW,LDW,ID,JD,LD,ND,XND,TITLE, readInType)
    !// --------------------------------------------------------------------
    !// Read resolution and version of restart file.
    !//
    !// Ole Amund Sovde, March 2016
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: IPAR, JPAR, LPAR, NPAR, NOTRPAR
    use cmn_ctm, only: IYEAR, IDAY, JYEAR, JDAY, JMON, JDATE, NTM, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW
    use cmn_chem, only: TNAME, TMASS
    use utilities, only: get_free_fileid
    use cmn_oslo, only: ZEROINIT, XZEROINIT, trsp_idx, chem_idx, &
          XTNAME, XTMASS, Xtrsp_idx, Xchem_idx, XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in):: FN_CON
    integer, intent(out) :: VERSION,IDW,JDW,LDW,ID,JD,LD,ND,XND,readInType
    character(len=*), intent(out):: TITLE

    !// Locals
    real(r8) :: SUMT
    real(r8), dimension(IPAR,JPAR,LPAR) :: INFIELD
    real(r8), dimension(LPAR,IPAR,JPAR) :: XINFIELD
    real(rMom), dimension(IPAR,JPAR,LPAR) :: &
         TSUT, TSVT, TSWT, TSUU, TSVV, TSWW, TSUV, TSUW, TSVW
    real(r8), dimension(IPAR+1) :: INXDEDG
    real(r8), dimension(JPAR+1) :: INYDEDG
    real(r8), dimension(LPAR+1) :: INETAA, INETAB

    integer :: &
         IIYEAR, IIDAY, JJYEAR, JJDAY, JJMON, JJDATE, I, &
         ioerr, file_id
    character(len=80) :: RTITLE2
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'OSLO_RESTARTFILE_INFO'
    !// --------------------------------------------------------------------

    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'
    write(6,'(a)') f90file//':'//subr// &
         ': Inspecting restart file: '//trim(FN_CON)

    !// Find a file number
    file_id = get_free_fileid()

    open(file_id,file=FN_CON,status='OLD',form='UNFORMATTED', &
         iostat=ioerr,action='read')
    if (ioerr .ne. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Error opening continuation file: '//trim(FN_CON)
       stop
    end if

    read(file_id) RTITLE2
    !// Return value for title
    TITLE = trim(RTITLE2)

    !// Read dimensions and transported tracers
    read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE,version
    if (ioerr .ne. 0) then
       !// Version is not specified, meaning version 2
       version = 2
       backspace(file_id)
       read(file_id,iostat=ioerr) IIYEAR,IIDAY,JJYEAR,JJDAY,JJMON,JJDATE
    end if

    if (version .le. 2) then
       read(file_id,iostat=ioerr) IDW,JDW,LDW,I,ID,JD,LD,ND,XND
    else
       !// Read native resolution and model resolution
       read(file_id,iostat=ioerr) IDW,JDW,LDW, ID,JD,LD, ND,XND
    end if

    close(file_id)

    !// Return values
    write(6,'(a,3i7)') ' Native resolution IDW/JDW/LDW:      ',IDW,JDW,LDW
    write(6,'(a,3i7)') ' Simulation resolution ID/JD/LD:     ',ID,JD,LD
    write(6,'(a,2i7)') ' Transported/non-transported ND/XND: ',ND,XND

    !// Check against current resolution to find read-in type
    !// 1: STD, 2: Txx, 3: T42
    readInType = -1
    if (VERSION .le. 2) then
       if (IDW.eq.128 .and. JDW.eq.64 .and. LDW.eq.60 .and. &
            ID.eq.128 .and. JD.eq.64 .and. LD.eq.60) then
          if (ID.eq.IPAR .and. JD.eq.JPAR .and. LD.eq.LPAR) then
             !// T42 resolution is same as current run
             readInType = 1
          else
             !// Interpolate from T42 run
             readInType = 3
          end if
       else
          !// Not set up to interpolate from other than T42
          write(6,'(a)') f90file//':'//subr// &
               ': Restart file is old and there is no read-in for '// &
               ' this resolution:'
          write(6,'(a)') 'File (window) resolution:  ',ID,JD,LD
          write(6,'(a)') 'Model (window) resolution: ',IPAR,JPAR,LPAR
          stop
       end if
    else
       !// Version > 2
       !// Not T42 - check against current run
       if (ID.eq.IPAR .and. JD.eq.JPAR .and. LD.eq.LPAR) then
          !// Resolution is same as current run
          readInType = 1
       else
          !// Interpolate from XX
          readInType = 2
       end if
    end if
    write(6,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// --------------------------------------------------------------------
  end subroutine OSLO_RESTARTFILE_INFO
  !// ----------------------------------------------------------------------








  !// ----------------------------------------------------------------------
end module stt_save_load
!//=========================================================================

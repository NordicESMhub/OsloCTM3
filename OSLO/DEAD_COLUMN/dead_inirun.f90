!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, March 2016
!//=========================================================================
!// Mineral dust initialisation.
!//=========================================================================
module dead_inirun
  !// ----------------------------------------------------------------------
  !// MODULE: dead_inirun
  !// DESCRIPTION: Initialises DEAD using settings from CTM3
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'dead_inirun.f90'
  !// ----------------------------------------------------------------------
  private
  public inirun
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine inirun( &
       delt, &          ! I [s] Time step
       latwts, &        ! I [-] latitude weighting
       calday, &        ! I [1-366]
       mbl_name, &       ! I Name of mobilisation map to read
       mbl_fudgefactor)  ! I Fudge factor matching mobilisation map
    !// --------------------------------------------------------------------
    !// Purpose: Initialize dust stuff
    !//
    !// History: Taken from match_dst/src/inirun routine (in main.F). 
    !//
    !// Rewritten for Oslo CTM2: Alf Grini (2002/2003)
    !// Checked for Oslo CTM3, Amund Sovde, October 2009
    !// Converted to f90 March 2016
    !// --------------------------------------------------------------------
    ! Choice of real8 /real4
    use dead_precision
    ! Grid definitions
    use pmgrid
    ! [mdl] Control variables, routines
    use dstctl
    ! [mdl] Nomenclature for outfld()
    use dstnm, only: dst_nm_cmn_ini
    ! [mdl] Physical constants for dust routines
    use dstcst, only:dst_cst_cmn_ini
    ! [mdl] Module initialization
    use dstcmnini, only:dst_msc_cmn_ini
    ! [mdl] Dust particle size distributions
    use dstpsd, only:dst_psd_ini,dst_psd_src_ini
    ! [mdl] Time-invariant boundary data sets
    use dsttibds, only:dst_tibds_ini
    ! [mdl] Time-varying boundary data sets
    use dsttvbds, only:dst_tvbds_ini
    ! Initialization of budget stuff
    use dstbdg, only:bdg_cmn_ini
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    !// INPUT
    real(r8), intent(in) :: &
         delt, &              !Timestep
         latwts(plat), &      !Latitude weights (0-1)
         calday               !Days for reading time varying input (real)
    character(len=*), intent(in) :: mbl_name ! Variable name for mobilisation
    real(r8), intent(in) :: mbl_fudgefactor  ! Fudge factor for mobilisation

    character(len=50) :: dstfile

    !// LOCAL VARIABLES 
    integer :: J
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'inirun'
    !// --------------------------------------------------------------------
      
    write(6,'(a)') f90file//':'//subr//': DEAD INIRUN'
    write(6,'(a,f12.3)') 'timestep:   ',delt
    write(6,'(a,f12.3)') 'calday:     ',calday
    !// Check if mbl_name is reasonable
    if (trim(mbl_name) .eq. 'NOT_SET') then
       write(6,'(a)') f90file//':'//subr//': mbl_name is not defined'
       write(6,'(a)') '  You need to define correct entry in the STV section'
       write(6,'(a)') '  of emission list, which sets the mbl_name.'
       stop
    end if
    if (mbl_fudgefactor .eq. 0._r8) then
       write(6,'(a)') f90file//':'//subr//': mbl_fudgefactor is not defined'
       write(6,'(a)') '  You need to define EFACT entry in the STV section'
       write(6,'(a)') '  of emission list.'
       stop
    end if
    write(6,'(a,es12.6)') 'mbl_fudgefactor: ',mbl_fudgefactor
    write(6,'(a)')        'mbl_name:        '//trim(mbl_name)

    do J = 1, plat
       write(6,'(a,i5,es24.15)') ' latwts',J,latwts(J)
    end do
    !// Initialize dust names (checked ok for CTM3)
    !// NB: dst_nm_cmn_ini() must be called before MATCH:src/deffld()
    call dst_nm_cmn_ini()

    !// Initialize time-invariant physical constants
    !// Variables are input to subroutine where variables are set for later use
    !// in dust subroutines. (checked ok for CTM3)
    call dst_cst_cmn_ini(mbl_fudgefactor)


    !//Initialize dust size distributions and miscellaneaous
    !// (checked ok for CTM3)
    call dst_psd_src_ini()

    !// Initialize size grid (checked ok for CTM3)
    call dst_psd_ini()

    !// Initialize miscellaneous common blocks (checked ok for CTM3)
    call dst_msc_cmn_ini()

    !// Initialize time-invariant boundary data from netCDF file
    !// (checked ok for CTM3)
    if (IPARW .eq. 128) then
       dstfile = 'dst_T42.nc'
    else if (IPARW .eq. 320) then
       dstfile = 'dst_T159.nc'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': Not defined dst_Txxx.nc file for this resolution! - Stopping'
       stop
    end if

    !// Adjusted to read mobilisation name from input file
    call dst_tibds_ini(trim(dstfile),mbl_name)

    !// NB: dst_tibds_ini() opens and closes the file, while
    !// dst_tvbds_ini() opens the file and leaves it open
    !// Initialize time-varying seasonal cycle data from netCDF file
    !// (checked ok for CTM3)
    call dst_tvbds_ini(trim(dstfile),calday)
    !//Initialize fields derived from external datasets
    !// call dst_lsm_ini()

    !// --------------------------------------------------------------------
  end subroutine inirun
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module dead_inirun
!//=========================================================================

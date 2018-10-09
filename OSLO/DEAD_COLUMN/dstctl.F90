! Purpose: Control variables, routines for dust model
! Currently nstep is only variable used outside of BXM:aer()
! Currently global models use only dst_bnr() routine
! Eventually this module is intended to contain process flags, e.g., dps_wet_flg=.true.

! Usage:
! use dstctl ! [mdl] Control variables, routines

! Requires dst.h for DST_MSS_BDG and many, many control tokens
!#include <dst.h> /* Dust preprocessor tokens */

module dstctl ! [mdl] Control routines for dust model
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::dst_bnr ! [fnc] Print banner with dust parameterization information
  save ! [stt] Changes to common variables are sticky
  
  ! Dust model control variables
  ! These variables are initialized in BXM:aer()
  ! CCM/MATCH control variables
  integer,public::mcdate ! [day] Current date in YYMMDD format
  integer,public::mcsec ! [s] Seconds past mcdate at 0Z
  integer,public::nbdate ! [day] Simulation start date in YYMMDD format
  integer,public::nbsec ! [s] Simulation start second relative to nbdate
  integer,public::ndcur ! [day] Current day number of simulation
  integer,public::nscur ! [s] Seconds relative to ndcur
  integer,public::nstep ! [idx] Current timestep number
  
contains
  
  subroutine dst_bnr() ! [fnc] Print banner with dust parameterization information
    ! Purpose: Print banner with dust parameterization information
    ! dst_bnr() is called by CCM:dst/dst_psd_ini(), MATCH:main()
    use dstgrd ! [mdl] Dust grid sizes
    use sng_mdl,only:ftn_strlen,ftn_cmd_ln_sng,ftn_prg_ID_mk ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::CVS_Date='$Date: 2003/04/15 14:41:35 $' ! [sng] Date string
    character(len=*),parameter::CVS_Id='$Id: dstctl.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $' ! [sng] CVS Identification
    character(len=*),parameter::CVS_Name='$Name:  $' ! [sng] File name string
    character(len=*),parameter::CVS_Revision='$Revision: 1.1 $' ! [sng] File revision string
    ! Local
    character(500)::cmd_ln ! [sng] Command line
    character(200)::prg_ID ! [sng] Program ID
    ! Initialize defaults
    ! dbg_lvl=dbg_off           ! dbg_lvl allocated in dbg.com, Option D
    call ftn_cmd_ln_sng(cmd_ln)
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
    ! Describe dust parameterization being used
    write (6,'(a)') '------------------------------------------------------------'
    write (6,'(a)') 'Mineral dust parameterization is active:'
    write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
    write (6,'(a)') '------------------------------------------------------------'
    ! write (6,'(a)') 'Mineral Dust Version ' // CVS_Name(8:30) // 'dated ' // CVS_Date(8:26) // ' GMT'
    write (6,'(i2,a,i2,a,i2)') &
         dst_nbr,' mineral dust bins are advected ' // &
         'species ',dst_idx_srt,' through ',dst_idx_end
    write (6,'(a)') 'Mobilization parameterization:'
    write (6,'(2x,a)') &
#ifdef AlG01
         'Saltation by Whi79 sandblasts dynamic size distribution of AlG01'
#else /* not AlG01 */
    'Saltation by Whi79 generates prescribed dust size distribution of Dal87'
#endif /* not AlG01 */

!// These are not used for Oslo CTM3. Skip printout.
!    write (6,'(a)') 'Radiative mode:'
!    write (6,'(2x,a)') &
!#ifdef DST_RAD
!         'Radiative effects of dust are computed'
!#else /* not DST_RAD */
!    'Radiative effects of dust are not computed'
!#endif /* not DST_RAD */
    write (6,'(2x,a)') &
!#ifdef DST_FDB
!         'Dust radiative effects are interactive, not diagnostic'
!#else /* not DST_FDB */
!    'Dust radiative effects are diagnostic, not interactive'
!#endif /* not DST_FDB */
!    write (6,'(2x,a)') 'Diagnostic optical depths will be archived'
!    write (6,'(a)') 'Land surface mode:'
!    write (6,'(2x,a)') &
!#ifdef DST_LSM
!         'Interactive LSM surface soil moisture is active'
!#else /* not DST_LSM */
!    'Interactive LSM surface soil moisture is turned off'
!#endif /* not DST_LSM */
!    write (6,'(a)') 'Chemistry mode:'
!    write (6,'(2x,a)') &
!#ifdef DST_CHM
!         'Heterogeneous chemistry is active (but incomplete)'
!#else /* not DST_CHM */
!    'Heterogeneous chemistry is turned off'
!#endif /* not DST_CHM */
!    write (6,'(a)') 'Debugging modes:'
!    write (6,'(2x,a)') &
#ifdef DST_DBG
         'Unvectorized run-time debugging blocks are active'
#else /* not DST_DBG */
    'Unvectorized run-time debugging blocks are turned off'
#endif /* not DST_DBG */
!    write (6,'(2x,a)') &
!#ifdef DST_MSS_BDG
!         'Diagnostic global mass budget is active'
!#else /* not DST_MSS_BDG */
!    'Diagnostic global mass budget is turned off'
!#endif /* not DST_MSS_BDG */
    return
  end subroutine dst_bnr                       ! end dst_bnr()
  
end module dstctl

      subroutine inirun(
     &     delt           ! I [s] Time step
     &     ,gravit        ! I [m/s2] gravity
     &     ,rair          ! I [J/K/kg] gas constant 
     &     ,latwts        ! I [-] latitude weighting
     &     ,rearth        ! I [m] Radius of earth
     &     ,calday        ! I [1-366]
     &     )
C
C     Purpose: Initialize dust stuff
C
C     History: Taken from match_dst/src/inirun routine (in main.F). 
C     In match, this is called from main.F. Here it is an own file
C     This means that a lot of the "use" statements in main.F in match
C     have to be repeated here
C
C     Rewritten for Oslo CTM2: Alf Grini (2002/2003)
C
      use precision  !Choice of real8 /real4
      use pmgrid     !Grid definitions (set in params.h in dust module)

      use dstctl                ! [mdl] Control variables, routines
      use dstnm,only:dst_nm_cmn_ini ! [mdl] Nomenclature for outfld()
      use dstcst,only:dst_cst_cmn_ini ! [mdl] Physical constants for dust routines
      use dstcmnini,only:dst_msc_cmn_ini ! [mdl] Module initialization
      use dstpsd,only:dst_psd_ini,dst_psd_src_ini ! [mdl] Dust particle size distributions
      use dsttibds,only:dst_tibds_ini ! [mdl] Time-invariant boundary data sets
      use dsttvbds,only:dst_tvbds_ini ! [mdl] Time-varying boundary data sets
C     previously only in match_dst/src/main.F :
      use dstbdg,only:bdg_cmn_ini     !Initialization
      
      implicit none

C     previously only in match_dst/src/main.F
#include <dst.h>

C     INPUT
      real(r8), intent(in)   :: delt              !Timestep
      real(r8), intent(in)   :: gravit            !Gravity (9.81 m s^(-2))
      real(r8), intent(in)   :: rair              !Gas constant air (J/kg/K)
      real(r8), intent(in)   :: latwts(plat)      !Latitude weights (0-1)
      real(r8), intent(in)   :: rearth            !Radius of earth (m)
      real(r8), intent(in)   :: calday            !Days for reading time varying input (real)

C     OUTPUT
C     NO PUTPUT
C     LOCAL VARIABLES 
c++added by Dave Newman 
c++alfgr: Probably not needed since we have delt in real(r8) already
C      real(r8) mydelt
c++csz
      
      write(6,*)'INIRUN',delt,rair,gravit,rearth,calday,latwts
c     Initialize dust names
c     NB: dst_nm_cmn_ini() must be called before MATCH:src/deffld()
      write(6,*)'timestep ',delt
      call dst_nm_cmn_ini()
c     Initialize time-invariant physical constants
C     Variables are input to subroutine where variables are set for later use
C     in dust subroutines.
      call dst_cst_cmn_ini(rair,gravit,rearth)
c++alfgr added 24/01/03
C     Initialize dust size distributions and miscellaneaous common blocks
      call dst_psd_src_ini()
c     Initialize size grid
      call dst_psd_ini()
c     Initialize miscellaneous common blocks
      call dst_msc_cmn_ini()
#ifdef DST_MSS_BDG
c     NB: MATCH reads delt as integer, CCM reads dtime as real(r8)
c     Initialize mass budget common block
      call bdg_cmn_ini(gravit,latwts,delt)
#endif /* not DST_MSS_BDG */
c     Initialize time-invariant boundary data from netCDF file
      call dst_tibds_ini('dst_bnd.nc')
c     NB: dst_tibds_ini() opens and closes the file, while
c     dst_tvbds_ini() opens the file and leaves it open
c     Initialize time-varying seasonal cycle data from netCDF file
      call dst_tvbds_ini('dst_bnd.nc',calday)
c     Initialize fields derived from external datasets
c     call dst_lsm_ini()
c--csz

      end subroutine inirun

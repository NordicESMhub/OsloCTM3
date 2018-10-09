c     $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/pmgrid.F,v 1.1 2003/04/15 14:41:36 alfgr Exp $ -*-fortran-*- 
      
c     Purpose: Supply spatial resolution parameters as Fortran variables
c     Pre-processor defined tokens are stored in Fortran variables here
      
c     Usage: 
c     use pmgrid ! [mdl] Spatial resolution parameters
      
      module pmgrid ! [mdl] Spatial resolution parameters
!//CTM3: Include parameters which used to be in params.h.
!//CTM3: Parameters in pmgrid.h were not used.
      use cmn_size, only: IPAR, JPAR, LPAR, NPAR_DUST, IDGRD, JDGRD,
     &     IPARW, JPARW !// For possible interpolation of input data
      use cmn_ctm, only: AREAXY, XDGRD, YDGRD, XDEDG, YDEDG,
     &     IMAP, JMAP
      implicit none
      save ! [stt] Changes to common variables are sticky

      !// Parameters from params_dust.h
      integer, parameter ::
     &     PLON = IPAR,
     &     PLAT = JPAR,
     &     PLEV = LPAR,
     &     PNATS= 0,
     &     DST_IDX_SRT = 1

      integer, parameter :: DST_NBR = NPAR_DUST
      integer, parameter :: PCNST = NPAR_DUST

!//CTM3: Not needed  
!#include <pmgrid.h> /* Spatial resolution parameters */
      integer, parameter ::
c     $     plon = PLON          ! number of longitudes
c     $     , plat = PLAT        ! number of latitudes
c     $     , plev = PLEV        ! number of vertical levels
c     $     , pcnst = PCNST      ! number of advected constituents
c     $     , pnats = PNATS      ! number of non-advected trace species
     $      plevp = plev + 1
     $     , plevd = 2*plev
     $     , nxpt = 1           ! no. of pts outside active domain of interpolant
     $     , jintmx = 1         ! number of extra latitudes in polar region
     $     , plond = plon !++alfgr + 1 + 2*nxpt
     $     , platd = plat !++alfgr + 2*nxpt + 2*jintmx
     $     , i1 = nxpt + 1      ! model starting longitude (3-d)
     $     , j1 = jintmx + nxpt + 1 ! model starting latitude (3-d)
     $     , j1m = j1 - 1       ! model starting offset (3-d)
     $     , padv = 2           ! number of miscellaneous advected fields
     $     , mxdynflds = 42     ! maximum number of dynamics input fields
c++csz
     $     , mxoutflds = 9*pcnst + pnats + 140 ! maximum number of history output fields
      
      logical, parameter ::
     $     masterproc = .true.  ! for CCM compatibility.  true => shared memory code


      end module pmgrid

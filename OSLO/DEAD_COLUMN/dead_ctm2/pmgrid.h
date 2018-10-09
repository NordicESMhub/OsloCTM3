c     $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/pmgrid.h,v 1.1 2003/04/15 14:41:36 alfgr Exp $ -*-fortran-*- 

c     Purpose: pmgrid.h sets all spatial resolution parameters
c     pmgrid.h is being deprecated in favor of pmgrid.F90
c     MATCH 3.X and 4.X contain fixed format code so MATCH portion of 
c     pmgrid.h must be syntactically valid fixed format Fortran
c     in order to work in files that #include it.
c     pmgrid.F is a wrapper that places pmgrid.h in a module so
c     it can be use'd by any pure, free-format code as well

c     Usage: 
c     #include <pmgrid.h> /* Spatial resolution parameters */ 

#include <params.h> /* Preprocessor tokens */

#ifdef BXM
c     Customized for offline aerosol driver
      integer,parameter::pcnst= PCNST ! [nbr] number of advected constituents
      integer,parameter::plon= PLON ! [nbr] number of longitudes
      integer,parameter::plat= PLAT ! [nbr] number of latitudes
      integer,parameter::plev= PLEV ! [nbr] number of vertical levels
      integer,parameter::plevp= plev + 1 ! [nbr] plev + 1
      integer,parameter::plond= plon ! [nbr] slt extended domain longitude
#endif /* endif BXM */

#if (defined CCM) && (!defined BXM)
c     
c     $Id: pmgrid.h,v 1.1 2003/04/15 14:41:36 alfgr Exp $
c     $Author: alfgr $
c     
c     
c     Grid point resolution parameters
      integer plon              ! number of longitudes
      integer plev              ! number of vertical levels
      integer plat              ! number of latitudes
      integer pcnst             ! number of constituents (including water vapor)
      integer pnats             ! number of non-advected trace species
      integer plevmx            ! number of subsurface levels

      integer plevp             ! plev + 1
      integer nxpt              ! no.of pts outside active domain of interpolant
      integer jintmx            ! number of extra latitudes in polar region
      integer plond             ! slt extended domain longitude
      integer platd             ! slt extended domain lat.
      integer p3d               ! dimensioning construct: num. of 3-d flds in /com3d/
      integer plevd             ! fold plev,pcnst indices into one
      integer i1                ! model starting longitude index
      integer j1                ! model starting latitude index
      integer numbnd            ! no.of latitudes passed N and S of forecast lat

      integer beglat            ! beg. index for latitudes owned by a given proc
      integer endlat            ! end. index for latitudes owned by a given proc
      integer beglatex          ! extended grid beglat
      integer endlatex          ! extended grid endlat
      integer numlats           ! number of latitudes owned by a given proc

      logical masterproc        ! Flag for (iam eq 0)

      parameter (plon   = PLON)
      parameter (plev   = PLEV)
      parameter (plat   = PLAT)
      parameter (pcnst  = PCNST)
      parameter (pnats  = PNATS)
      parameter (plevmx = 4)
      parameter (plevp  = plev + 1)
      parameter (nxpt   = 1)
      parameter (jintmx = 1)
      parameter (plond  = plon + 1 + 2*nxpt)
      parameter (platd  = plat + 2*nxpt + 2*jintmx)
      parameter (p3d    = 3 + pcnst + pnats)
      parameter (plevd  = plev*p3d)
      parameter (i1     = 1 + nxpt)
      parameter (j1     = 1 + nxpt + jintmx)
      parameter (numbnd = nxpt + jintmx)
                                ! 
#if ( defined SPMD )
      common/spmdlats/beglat  ,endlat  ,numlats ,beglatex,endlatex 
      common/spmdlats/masterproc
#else 
      parameter (beglat   = 1)
      parameter (endlat   = plat)
      parameter (numlats  = plat)
      parameter (beglatex = 1)
      parameter (endlatex = platd)
      parameter (masterproc = .true.)
#endif
                                ! 
#endif /* endif CCM */
      
#if (!defined CCM) && (!defined BXM)

#ifndef PARAMS_H
#include <params.h>
#endif
      
c     Set kind of real variables to have at least 12 digits of precision.
c     integer, parameter :: REALKIND = selected_real_kind( p = 12 )
c     integer, parameter :: REALKIND = selected_real_kind( p = RPREC )
      
c     Basic grid point resolution parameters
      
      integer, parameter ::
     $     plon = PLON          ! number of longitudes
     $     , plat = PLAT        ! number of latitudes
     $     , plev = PLEV        ! number of vertical levels
     $     , pcnst = PCNST      ! number of advected constituents
     $     , pnats = PNATS      ! number of non-advected trace species
     $     , plevp = plev + 1
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
#ifdef DST
     $     , mxoutflds = 9*pcnst + pnats + 140 ! maximum number of history output fields
#else 
     $     , mxoutflds = 5*pcnst + pnats + 70 ! maximum number of history output fields
#endif /* not DST */
c--csz
      
      logical, parameter ::
     $     masterproc = .true.  ! for CCM compatibility.  true => shared memory code
      
#endif /* endif MATCH */

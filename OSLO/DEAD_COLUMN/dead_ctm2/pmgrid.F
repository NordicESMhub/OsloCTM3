c     $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/pmgrid.F,v 1.1 2003/04/15 14:41:36 alfgr Exp $ -*-fortran-*- 
      
c     Purpose: Supply spatial resolution parameters as Fortran variables
c     Pre-processor defined tokens are stored in Fortran variables here
      
c     Usage: 
c     use pmgrid ! [mdl] Spatial resolution parameters
      
      module pmgrid ! [mdl] Spatial resolution parameters
      implicit none
      save ! [stt] Changes to common variables are sticky
  
#include <pmgrid.h> /* Spatial resolution parameters */

      end module pmgrid

      subroutine qneg3(subnam  ,lat     ,q       )

      use pmgrid
      use dstgrd
      use precision

C     This subroutine is really here just to have something to compile
C     It does not do anything, but DEAD calls it all the time

C     Alf Grini, april 2003

      implicit none

      character(len=80)     :: subnam
      integer     :: lat
      real(r8)    :: q(plon,plev,dst_nbr)

      end subroutine

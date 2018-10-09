
c-----------------------------------------------------------------------
      Program PMEANC
c-----------------------------------------------------------------------
      include 'cmn.f'

      real*8  ZNN,ZIJ, SOLDEC,SOLDIS, TL(LPAR),TSUM,ZOFL(LPAR+1),DELZ
c---day number, starting and ending day
      integer NDAY, NDAYI, NDAYE, NMET, NN,I,J,L
      character*80 INFILE1
      logical  LYEAR

C------- input and initialization, includes Chem input -----------------
      INFILE1 = 'P_meanT42L60.dat'

      open (10, file=INFILE1, form='formatted',status='old')
      read (10,*)
      read (10,'(10f8.2)') ((PMEAN(I,J),I=1,IPAR),J=1,JPAR)
      read (10,*)
      read (10,'(10f8.4)') (AREAXY(1,J),J=1,JPAR)
      read (10,*)
      read (10,'(10f8.2)') (XDGRD(I),I=1,IPAR)
      read (10,*)
      read (10,'(10f8.2)') (YDGRD(J),J=1,JPAR)
      read (10,*)
      read (10,'(10f8.3)') (ZOFL(L),L=1,LPAR+1)
      read (10,*)
      read (10,'(10f8.2)') (TL(L),L=1,LPAR)
      read (10,*)
      read (6,'(10f8.3)') (ETAA(L),L=1,LPAR+1)
      read (10,*)
      read (6,'(10f8.5)') (ETAB(L),L=1,LPAR+1)
      close(10)

      stop
      end

!------------------------------------------------------------------------------
!          UCI CTM core p-7.1 (1/2013)                                        !
!------------------------------------------------------------------------------
!
! !MODULE: CMN_FLAG
!
! !DESCRIPTION: CMN_FLAG contains run time settings in CTM
!
!
! !INTERFACE:
!
      MODULE CMN_FLAG_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : NPAR

      IMPLICIT NONE
      PUBLIC

!------------------------------------------------------------------------------
      ! input from Lxxx.inp
!-----------------------------------------------------------------------

      ! optin for using met fields
      !LOGICAL  LLPYR      ! allow leap year (true: looks for day=901 for Feb 29)
      !LOGICAL  LFIXMET    ! annually recycle met field (do not allow leap year)

      ! optin for writing out flux diagnostics
      !LOGICAL  LFLXDG

!      INTEGER  CLDFLAG, NRANDO  ! # of ICA's used for LCDFLAG(5)
!      1 : Clear sky J's
!      2 : Averaged cloud cover
!      3 : cloud-fract**3/2, then average cloud cover
!      4 : Average direct solar beam over all ICAs, invert to get clouds
!      5 : Random select NRANDO ICA's from all(Independent Column Atmos.)
!      6 : Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!      7 : Use all (up to 4) QCAs (average clouds within each Q-bin)
!      8 : Calcluate J's for ALL ICAs (up to 20,000 per cell!)

!-----------------------------------------------------------------------
      ! input from Ltracer.inp
!-----------------------------------------------------------------------
      !LOGICAL  LCONT    ! continuation run, read restart file
      !LOGICAL  LINIT    ! if LCONT=T && LINIT=T, use restart file as
                        ! initialization, date info from restart is not used
      !LOGICAL  LSTRAT   ! use Linoz for stratospheric chemistry
      LOGICAL  LINOZ3   ! True: use Linoz version 3 (including O3S,n2o,noy,ch4)
                        ! False: Linoz version 2 (O3strat only)
      LOGICAL  LN2OCH4E ! for Linoz version 3
                        ! True:  use N2O and CH4 emission from Ltracer.inp
                        ! False: keep constant N2O & CH4 vmr at sfc L=1,LZLBO3
      !/ /CTM3: disabled ASAD
      !LOGICAL  LASAD    ! if true: full ASAD tropospheric chemistry
      LOGICAL  LZONE(NPAR)  ! True: scale tracer at beginning of the run

      END MODULE CMN_FLAG_MOD

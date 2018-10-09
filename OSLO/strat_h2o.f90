!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Stratpspheric H2O.
!//=========================================================================
module strat_h2o
  !// ----------------------------------------------------------------------
  !// MODULE: strat_h2o
  !// DESCRIPTION: The purpose of this module is to combine a prognostic
  !// H2O in the stratospehere with a static H2O (meteorological data) in
  !// the troposphere.
  !//
  !// From the EU project EUROHYDROS, it was recognized that defining a
  !// flux from the troposphere, based on metdata specific humidity fields,
  !// was troublesome. It caused too large input to the stratosphere during
  !// convection, and trying to constrain it became rather unphysical.
  !//
  !// A new method is now implemented, where the H2O tracer 114 is a purely
  !// stratospheric tracer, set to a fixed mixing ratio at the tropopause
  !// level.
  !//
  !//
  !// A switch (LOLD_H2OTREATMENT) controls whether to use old or new
  !// treatment of stratospheric H2O. Either way, H2O_MOLEC is now
  !// not in use in the stratosphere. We use component
  !// 114 and 125, both for old and new treatment.
  !//
  !// Contains:
  !//   subroutine set_trop_h2o_b4trsp_clim
  !//   subroutine reset_trop_h2o
  !//   subroutine set_strat_h2o_b4chem
  !//   subroutine set_d_h2o
  !//   subroutine strat_h2o_init
  !//   subroutine strat_h2o_ubc
  !//   subroutine strat_h2o_ubc2
  !//   subroutine strat_h2o_max
  !//   subroutine strh2o
  !//
  !// Amund Sovde, September 2013
  !//   New treatment of transported H2O. Tracer 114 is now only stratospheric
  !//   H2O, and is set to a fixed mixing ratio at tropopause level, and
  !//   zeroed in the lower troposphere.
  !// Amund Sovde, December 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, NPAR, IDBLK, JDBLK, MPBLK, LWEPAR
  use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, YDGRD
  use cmn_met, only: CWETE, Q, T, CLDFR
  use cmn_parameters, only: AVOGNR, M_AIR
  use cmn_oslo, only: trsp_idx, Xtrsp_idx, DV_IJ, AIRMOLEC_IJ, LMTROP
  use utilities_oslo, only: H2O_SAT
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

  !// Switch to use old stratospheric H2O treatment.
  logical, parameter :: LOLD_H2OTREATMENT = .true.

  !// Sum of H2+2CH4+H2O:e sumH2 = 7.72ppm (Zöger et al, JGR 1999, vol 104,
  !// D1, pp 1817-1825) OLD VALUE: 6.97d-6
  real(r8), parameter :: sumH2 = 7.72e-6_r8

  !// For converting from mass to concentration
  real(r8), parameter :: &
       Z_CH4  = AVOGNR * 1.e-3_r8 / 16._r8, &
       Z_H2   = AVOGNR * 1.e-3_r8 /  2._r8, &
       Z_H2O  = AVOGNR * 1.e-3_r8 / 18._r8
  !// For converting from concentration to mass
  real(r8), parameter :: X_H2O  = 18._r8 / (AVOGNR * 1.e-3_r8)

  !// For old treatment
  real(r8) :: str_h2o(LPAR,IDBLK,JDBLK,MPBLK)
  real(r8) :: d_H2O(LPAR,IDBLK,JDBLK,MPBLK)

  !// Climatological H2O
  real(r8) :: H2O_TP_CLIM(IPAR,JPAR,12)
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'strat_h2o.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public strat_h2o_init, set_strat_h2o_b4chem, &
       reset_trop_h2o, strat_h2o_ubc2, strat_h2o_max, &
       set_d_h2o, set_h2_eurohydros, &
       LOLD_H2OTREATMENT, str_h2o, zc_strh2o, sumH2, &
       strat_h2o_init_clim, set_trop_h2o_b4trsp_clim
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine set_trop_h2o_b4trsp_clim(BTT,BTTBCK,AIRB,BTEM,MP)
    !// --------------------------------------------------------------------
    !// Set H2O in troposphere before transport.
    !// In the 5 uppermost tropospheric layers, H2O is set from climatological
    !// values. These values are currently set in routine
    !// strat_h2o_init_clim, and can be set from either model data or
    !// as fixed mixing ratios.
    !//
    !// Carried out in IJ-blocks.
    !//
    !// Amund Sovde, May 2013
    !// --------------------------------------------------------------------
    use cmn_ctm, only: JMON
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8),dimension(LPAR,IDBLK,JDBLK),intent(in)  :: AIRB, BTEM
    !// Input/Output
    real(r8),dimension(LPAR,NPAR,IDBLK,JDBLK),intent(inout) :: BTT, BTTBCK

    !// Locals
    integer :: I,J,L,II,JJ,L1,L2,LTPMIN,LTPMAX
    real(r8) :: c_h2o_max, m_h2o_max, m_h2o
    !// --------------------------------------------------------------------

    !// Skip this if old H2O treatment
    if (LOLD_H2OTREATMENT) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       LTPMIN = minval(LMTROP(:,J))
       LTPMAX = maxval(LMTROP(:,J))

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Set some levels to be used, currently 5 layers
          L2 = max(1,LMTROP(I,J)-4)
          L1 = LMTROP(I,J)

          !// Zero (almost) troposphere to avoid large spikes into
          !// the stratosphere
          do L = 1, L2-1
             BTT(L,trsp_idx(114),II,JJ) = 1.e-10_r8
          end do

          !// Set tropopause values
          do L = L2, L1
             !// Set H2O mass from vmr
             m_h2o = H2O_TP_CLIM(I,J,JMON) * AIRB(L,II,JJ) * 18._r8 / M_AIR
             !// But limit to saturation value
             c_h2o_max = H2O_SAT(BTEM(L,II,JJ))
             !// Convert to mass
             m_h2o_max = c_h2o_max * DV_IJ(L,II,JJ,MP) * X_H2O

             BTT(L,trsp_idx(114),II,JJ) = min(m_h2o, m_h2o_max)
          end do
       
       end do
    end do

    !// Update BTTBCK
    BTTBCK(:,trsp_idx(114),:,:) = BTT(:,trsp_idx(114),:,:)

    !// --------------------------------------------------------------------
  end subroutine set_trop_h2o_b4trsp_clim
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine reset_trop_h2o(LEV,AIRUV)
    !// --------------------------------------------------------------------
    !// Reset H2O in troposphere after transport, i.e. to Q*AIR plus
    !// aircraft H2O tracer.
    !//
    !// Amund Sovde, January 2011
    !// --------------------------------------------------------------------
    use cmn_precision, only: rMom
    use cmn_ctm, only: STT, SUT,SVT,SWT, SUU,SVV,SWW, SUV,SUW,SVW
    use cmn_met, only: T
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LEV
    real(r8), intent(in)  :: AIRUV(IPAR,JPAR)
    !// Locals
    integer :: I,J,MP,II,JJ, TNR
    real(r8) :: mfrac,m_h2o_max,c_h2o_max
    !// --------------------------------------------------------------------

    !// This routine will produce way too much H2O in stratosphere,
    !// and should not be used.
    stop f90file//': reset_trop_h2o: OLD ROUTINE - SHOULD NOT BE USED!'

    !// Skip this if old H2O treatment
    if (LOLD_H2OTREATMENT) return

    !// Set 114 if included
    if (trsp_idx(114) .gt. 0) then
       do J = 1, JPAR
          do I = 1, IPAR
             if (LEV .le. LMTROP(I,J)) then
                STT(I,J,LEV,trsp_idx(114)) = Q(I,J,LEV) * AIRUV(I,J)
             end if
          end do
       end do

       !// Aircraft emissions included?
       if (trsp_idx(148) .gt. 0) then
          do J = 1, JPAR
             do I = 1, IPAR
                if (LEV .le. LMTROP(I,J)) then
                   STT(I,J,LEV,trsp_idx(114)) = STT(I,J,LEV,trsp_idx(114)) &
                        + STT(I,J,LEV,trsp_idx(148))
                end if
             end do
          end do
       end if
    end if


    !// Limit stratospheric H2O to 2% supersat, after transport
    if (trsp_idx(114) .gt. 0) then
      TNR = trsp_idx(114)
      do MP = 1, MPBLK
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1

            if (LEV .gt. LMTROP(I,J) .and. LEV.lt.LMTROP(I,J)+5) then
              !// Find saturation value for the given temperature
              c_h2o_max = H2O_SAT(T(I,J,LEV))

              !// to mass
              m_h2o_max = 1.02_r8 * c_h2o_max * DV_IJ(LEV,II,JJ,MP) * X_H2O

              mfrac = min(1._r8, m_h2o_max/STT(I,J,LEV,TNR))

              STT(I,J,LEV,TNR) = mfrac * STT(I,J,LEV,TNR)
              SUT(I,J,LEV,TNR) = real(mfrac,rMom) * SUT(I,J,LEV,TNR)
              SVT(I,J,LEV,TNR) = real(mfrac,rMom) * SVT(I,J,LEV,TNR)
              SWT(I,J,LEV,TNR) = real(mfrac,rMom) * SWT(I,J,LEV,TNR)
              SUU(I,J,LEV,TNR) = real(mfrac,rMom) * SUU(I,J,LEV,TNR)
              SVV(I,J,LEV,TNR) = real(mfrac,rMom) * SVV(I,J,LEV,TNR)
              SWW(I,J,LEV,TNR) = real(mfrac,rMom) * SWW(I,J,LEV,TNR)
              SUV(I,J,LEV,TNR) = real(mfrac,rMom) * SUV(I,J,LEV,TNR)
              SUW(I,J,LEV,TNR) = real(mfrac,rMom) * SUW(I,J,LEV,TNR)
              SVW(I,J,LEV,TNR) = real(mfrac,rMom) * SVW(I,J,LEV,TNR)

            end if
          end do
        end do
      end do
    end if

    !// --------------------------------------------------------------------
  end subroutine reset_trop_h2o
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_strat_h2o_b4chem(AIRB,BTT,MP)
    !// --------------------------------------------------------------------
    !// Non-prognostic stratospheric H2O treatment; set before chemistry,
    !// before PSC physics.
    !// Based on sum of H2O+2CH4+H2. Updated from the old Oslo CTM2 sum
    !// (6.97ppmv) to the new defined at the top of the module.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LOSLOCSTRAT
    use cmn_oslo, only: XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MP
    real(r8), intent(in) :: AIRB(LPAR,IDBLK,JDBLK)
    !// Input/Output
    real(r8), intent(inout) :: BTT(LPAR,NPAR,IDBLK,JDBLK)

    !// Locals
    real(r8) :: c_h2o_max, m_h2o_max, m_h2o_met, maxfrac, TEMP
    real(r8)  :: m_h2o, m_h2, m_ch4, m_h2os
    integer :: I,J,L,II,JJ,L1,L2,XID
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'set_strat_h2o_b4chem'
    !// --------------------------------------------------------------------

    !// If no stratospheric chemistry, skip H2O
    if (.not. LOSLOCSTRAT) return

    !// Skip this if new H2O treatment
    if (.not. LOLD_H2OTREATMENT) return

    !// Initialize
    if (trsp_idx(114) .gt. 0) then
       write(6,'(a)') f90file//':'//subr// &
            ': H2O should not be transported when using old strat. H2O'
       stop
    end if

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          do L = 1, LPAR
             !// Save H2O as mass [kg]
             STR_H2O(L,II,JJ,MP) = Q(I,J,L) * AIRB(L,II,JJ)
          end do
       end do
    end do


    !// Stratosphere: as for old treatment in CTM2
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          do L = LMTROP(I,J)+1, LPAR-1
             !// CH4
             if (trsp_idx(46) .gt. 0) then
                M_CH4 = BTT(L,trsp_idx(46),II,JJ)
             else if (Xtrsp_idx(46).gt.0) then
                M_CH4 = XSTT(L,Xtrsp_idx(46),I,J)
             else
                write(6,'(a)') f90file//':'//subr//': CH4 is not simulated!'
                stop
             end if
             !// H2
             if (trsp_idx(113).gt.0) then
                M_H2 = BTT(L,trsp_idx(113),II,JJ)
             else if (Xtrsp_idx(113).gt.0) then
                M_H2 = XSTT(L,Xtrsp_idx(113),I,J)
             else
                write(6,'(a)') f90file//':'//subr//': H2 is not simulated!'
                stop
             end if

             !// From mass to concentration
             M_CH4 = M_CH4 * Z_CH4 / DV_IJ(L,II,JJ,MP)
             M_H2  = M_H2 * Z_H2 / DV_IJ(L,II,JJ,MP)

             !// Set H2O from sum of H2.
             !// Do not subtract H2Os, since it is a static field, 
             !// calculated from H2O and T.
             m_h2o = sumH2 * AIRMOLEC_IJ(L,II,JJ,MP) &
                      - 2._r8 * M_CH4 & !// - 2*CH4
                      - M_H2            !// - H2

             if (m_h2o .le. 0._r8) then
                write(6,'(a,3i5)') f90file//':'//subr// &
                     ': H2O <= 0 in stratosphere A',i,j,l
                write(*,'(a8,2es20.10)') '  CH4:  ', M_CH4, &
                     BTT(L,trsp_idx(46),II,JJ) / AIRB(L,II,JJ) * M_AIR / 16._r8
                write(*,'(a8,es20.10)') '  H2:   ', M_H2
                write(*,'(a8,es20.10)') '  sumH2:', &
                     sumH2 * AIRMOLEC_IJ(L,II,JJ,MP)
                stop
             end if
             !// Update mass [kg]
             STR_H2O(L,II,JJ,MP) = m_h2o * DV_IJ(L,II,JJ,MP) * X_H2O
          end do
       end do
    end do

    !// Adjust for frozen H2O
    d_H2O(LPAR,:,:,MP) = 0._r8 !// Just to make sure (in LPAR)
    !// d_H2O is calculated in psps_psc and passed to subroutine
    !// set_d_h2o here.

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Loop through stratosphere, omit uppermost layer
          do L = LMTROP(I,J)+1, LPAR-1

             m_h2o = STR_H2O(L,II,JJ,MP)
             !// Update H2O mass [kg]
             STR_H2O(L,II,JJ,MP) = &
                  max(0._r8, STR_H2O(L,II,JJ,MP) + d_H2O(L,II,JJ,MP))


             if (STR_H2O(L,II,JJ,MP) .eq. 0._r8) then
                write(6,'(a,i3,i4,3i3,3es13.5)') f90file//':'//subr// &
                     ': ZERO H2O: ',&
                     MP,I,J,L,LMTROP(I,J),m_h2o,d_H2O(L,II,JJ,MP)
                stop
             end if

          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine set_strat_h2o_b4chem
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_d_h2o(ii,jj,mp,LMAX,dh2o)
    !// --------------------------------------------------------------------
    !// Set d_H2O [molec/cm3]. Called from PSC physics.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ii,jj,mp, LMAX
    real(r8), dimension(LPAR), intent(in) :: dh2o
    !// Locals
    integer :: L
    !// --------------------------------------------------------------------
    !// Concentration of sedimented H2O [molec/cm3]
    d_H2O(1:LMAX,II,JJ,MP) = dh2o(1:LMAX) * DV_IJ(1:LMAX,II,JJ,MP) * X_H2O
    !// --------------------------------------------------------------------
  end subroutine set_d_h2o
  !// ----------------------------------------------------------------------



  !// --------------------------------------------------------------------
  subroutine strat_h2o_init()
    !// --------------------------------------------------------------------
    !// Initialize H2 and H2O globally.
    !// H2 based on annual mean of the last Eurohydros simulation.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_size, only: LOSLOCSTRAT
    use cmn_ctm, only: AIR, STT, AREAXY
    use cmn_met, only: ZOFLE
    use cmn_oslo, only: XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8)  :: m_h2o, m_h2, m_ch4, old_h2, metdata_h2o, new_h2o, &
         max_h2o, c_h2o_max, maxfrac, rsum, am_box, dv_box,sh2
    integer :: I,J,L,II,JJ,MP

    !// H2 from file
    character(len=80) :: FNAME
    logical :: fnr_ok
    integer :: ifnr, io_err, IR, JR, LR
    real(r4)  :: STTH2(IPAR,JPAR,LPAR)
    real(r8), parameter :: AIR2AM = AVOGNR * 1.e-3_r8 / M_AIR
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'strat_h2o_init'
    !// --------------------------------------------------------------------

    !// If no stratospheric chemistry, skip H2O
    if (.not. LOSLOCSTRAT) return

    !// Test if H2O, CH4 and H2 are included and how they are treated.
    if (.not.LOLD_H2OTREATMENT .and. &
         (trsp_idx(113).le.0 .or. trsp_idx(114).le.0 .or. trsp_idx(46).le.0 &
          .or. trsp_idx(125).le.0) ) then
       write(6,'(a)') f90file//':'//subr//': Initialization problem!'
       print*,'  New H2O treatment needs CH4, H2, H2O and H2Os to be transported!'
       print*,'  Xtrsp_idx( 46)/trsp_idx( 46):',Xtrsp_idx(46),trsp_idx(46)
       print*,'  Xtrsp_idx(113)/trsp_idx(113):',Xtrsp_idx(113),trsp_idx(113)
       print*,'  Xtrsp_idx(114)/trsp_idx(114):',Xtrsp_idx(114),trsp_idx(114)
       print*,'  Xtrsp_idx(125)/trsp_idx(125):',Xtrsp_idx(125),trsp_idx(125)
       stop
    else if (LOLD_H2OTREATMENT .and. &
         (trsp_idx(46).le.0 .or. &
          Xtrsp_idx(113).le.0 .or. Xtrsp_idx(125).le.0) ) then
       write(6,'(a)') f90file//':'//subr//': Initialization problem!'
       print*,'  Old H2O treatment needs CH4 to be transported, but'
       print*,'  H2, and H2Os must not be transported!'
       print*,'  Xtrsp_idx( 46)/trsp_idx( 46):',Xtrsp_idx(46),trsp_idx(46)
       print*,'  Xtrsp_idx(113)/trsp_idx(113):',Xtrsp_idx(113),trsp_idx(113)
       print*,'  Xtrsp_idx(125)/trsp_idx(125):',Xtrsp_idx(125),trsp_idx(125)
       stop
    end if


    !// H2 at surface is continuously set by 2d-data, to mimick emissions.
    !// Emissions can be included later.
    write(6,'(a)') f90file//':'//subr//': Initializing H2O'

    !// Initialize diagnose of sedimented H2O
    d_H2O(:,:,:,:) = 0._r8

    !// If STT is already set, skip this part
    if (trsp_idx(114) .gt. 0) then
       if (maxval(STT(:,:,:,trsp_idx(114))).gt.1._r8) then
          write(*,'(a)') '* H2O has already been initialized'
          return
       end if
    end if


    !// Initialize H2O as mass
    if (trsp_idx(114) .gt. 0) then
       !// Tracer 114 is purely stratospheric H2O, with zero values
       !// in troposphere.
       STT(:,:,:,trsp_idx(114)) = 1.e-10_r8
    else 
       !// H2O not transported, set str_h2o.
       do MP=1,MPBLK
         do J = MPBLKJB(MP), MPBLKJE(MP)
           JJ = J - MPBLKJB(MP) + 1
           do I = MPBLKIB(MP), MPBLKIE(MP)
             II = I - MPBLKIB(MP) + 1
             do L = 1, LPAR
               str_h2o(L,II,JJ,MP) = max(1.e-15_r8, Q(I,J,L) * AIR(I,J,L))
             end do
           end do
         end do
       end do
    end if

    !// For old treatment, set solid H2O to sero
    if (Xtrsp_idx(125) .gt. 0) XSTT(:,Xtrsp_idx(125),:,:) = 0._r8


    !// Stratosphere as for old treatment in CTM2
    !// What is the change in stratospheric H2O burden?
    metdata_h2o= 0._r8
    new_h2o  = 0._r8
    max_h2o = 0._r8
    do MP = 1, MPBLK

      do J = MPBLKJB(MP), MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = LMTROP(I,J)+1, LPAR-1

            !// CH4
            if (trsp_idx(46).gt.0) then
               M_CH4 = STT(I,J,L,trsp_idx(46))
            else
               M_CH4 = XSTT(L,Xtrsp_idx(46),I,J)
            end if
            !// H2
            if (trsp_idx(113).gt.0) then
               M_H2 = STT(I,J,L,trsp_idx(113))
            else
               M_H2 = XSTT(L,Xtrsp_idx(113),I,J)
            end if
            !// No need for H2Os

            !// Box volume and air density are not set yet, they are
            !// set in IJ-block loop.
            !// Box volume [m3]
            DV_BOX = AREAXY(I,J) * (ZOFLE(L+1,I,J) - ZOFLE(L,I,J))
            !// Calculate air density [molec/cm^3]
            AM_BOX = AIR(I,J,L) / DV_BOX * AIR2AM


            !// From mass to concentration
            M_CH4 = M_CH4 * Z_CH4 / DV_BOX
            M_H2  = M_H2 * Z_H2 / DV_BOX

            !// Set H2O from sum of H2.
            m_h2o = sumH2 * AM_BOX &
                      - 2._r8 * M_CH4 & !// - 2*CH4
                      - M_H2            !// - H2

            if (m_h2o .le. 0._r8) then
               write(6,'(a,3i5)') f90file//':'//subr// &
                    ': H2O <= 0 in stratosphere B',i,j,l
               write(6,'(a8,2es20.10)') '  CH4:  ',M_CH4, M_CH4 / AM_BOX
               write(6,'(a8,2es20.10)') '  H2:   ',M_H2, M_H2 / AM_BOX
               write(6,'(a8,2es20.10)') '  sumH2:',sumH2 * AM_BOX, sumH2
               stop
            end if


            !// To mass
            m_h2o = m_h2o * DV_BOX * X_H2O
            if (trsp_idx(114).gt.0) then
               STT(I,J,L,trsp_idx(114)) = m_h2o
            else
               str_h2o(L,II,JJ,MP) = m_h2o
            end if
            metdata_h2o = metdata_h2o + Q(I,J,L)*AIR(I,J,L)
            new_h2o = new_h2o + m_h2o
            if (m_h2o/air(i,j,l) .gt. max_h2o) then
               max_h2o = m_h2o / air(i,j,l)
            end if
          end do
        end do
      end do
    end do
    if (trsp_idx(114).gt.0) then
       m_h2o = minval(STT(:,:,:,trsp_idx(114)))
    else
       m_h2o = minval(str_h2o)
    end if
    if (m_h2o .le. 0._r8) then
       write(6,'(a,es16.6)') f90file//':'//subr// &
            ': H2O should be positive!',m_h2o
       stop
    end if


    !// Additional print-outs
    write(*,'(a,es14.7,a,es14.7)') &
         '* METDATA strat.H2O: ',metdata_h2o,' New strat.H2O: ',new_h2o
    write(*,'(a,es14.7)') &
         '* Max stratospheric H2O [ppmv]: ',max_h2o * M_AIR / 18._r8 * 1.e6_r8


    !// --------------------------------------------------------------------
  end subroutine strat_h2o_init
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine strat_h2o_ubc(BTT,BTTBCK,AIRB,MP)
    !// --------------------------------------------------------------------
    !// Set upper boundary condition of H2O.
    !// As for old CTM2 treatment, but updated sumH2.
    !// Requires transported CH4 and H2.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/Output
    integer, intent(in) :: MP
    real(r8), intent(in) :: AIRB(LPAR,IDBLK,JDBLK)
    real(r8),dimension(LPAR,NPAR,IDBLK,JDBLK),intent(inout) :: BTT, BTTBCK
    !// Locals
    real(r8) :: m_h2o, m_h2, m_ch4, msed
    integer :: I,J,II,JJ
    !// --------------------------------------------------------------------

    !// Skip new UBC for old H2O treatment
    if (LOLD_H2OTREATMENT) return

    !// In case stratospheric chemistry is off.
    if (trsp_idx(114) .lt. 0) return

    if (trsp_idx(46).le.0 .or. trsp_idx(113).le.0) then
       print*,'* strat_h2o.f90: Upper stratospheric boundary condition of'
       print*,'  H2O needs CH4 and H2 to be transported!'
       print*,'  trsp_idx(46):',trsp_idx(46)
       print*,'  trsp_idx(113):',trsp_idx(113)
       stop
    end if

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// CH4
          M_CH4 = BTT(LPAR,trsp_idx( 46),II,JJ)
          !// H2
          M_H2  = BTT(LPAR,trsp_idx(113),II,JJ)

          !// From mass to concentration
          M_CH4 = M_CH4 * Z_CH4 / DV_IJ(LPAR,II,JJ,MP)
          M_H2  = M_H2 * Z_H2 / DV_IJ(LPAR,II,JJ,MP)

          !// Set H2O from sum of H2 [molec/cm3].
          m_h2o = sumH2 * AIRMOLEC_IJ(LPAR,II,JJ,MP) &
                - 2._r8 * M_CH4 & !// - 2*CH4
                - M_H2            !// - H2
          if (m_h2o .le. 0._r8) then
             print*,'* H2O <= 0 in stratosphere C',i,j,lpar
             print*,'  CH4:  ',M_CH4
             print*,'  H2:   ',M_H2
             print*,'  sumH2:',sumH2 * AIRMOLEC_IJ(Lpar,II,JJ,MP)
             stop
          end if
          !// To mass
          BTT(LPAR,trsp_idx(114),II,JJ) = m_h2o * DV_IJ(LPAR,II,JJ,MP) * X_H2O
       end do
    end do

    !// Update BTTBCK
    BTTBCK(LPAR,trsp_idx(114),:,:) = BTT(LPAR,trsp_idx(114),:,:)

    !// --------------------------------------------------------------------
  end subroutine strat_h2o_ubc
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine strat_h2o_ubc2(BTT,BTTBCK,AIRB,MP)
    !// --------------------------------------------------------------------
    !// Set upper boundary condition of H2O/H2/CH4.
    !// 
    !// Assumption is to use same mixing ratio as level below for H2.
    !// First approach was to use 2D CH4, and assume (H2O+H2+2xCH4) to
    !// be constant from LPAR-1 to LPAR, and calculate H2O.
    !// However, when CH4 increase below (due to emissions), the sum
    !// gets larger. Because CH4 is unchanged from 2D data, H2O then
    !// increased in LPAR. This was not correct.
    !//
    !// So the second approach is to set all three with same mixing ratio at
    !// LPAR as in LPAR-1. I.e. it is not necessary to introduce the sum.
    !//
    !// Amund Sovde, May-August 2013
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input/Output
    integer, intent(in) :: MP
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: AIRB
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT, BTTBCK
    !// Locals
    real(r8)  :: m_h2o, sumLP1
    integer :: I,J,II,JJ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'strat_h2o_ubc2'
    !// --------------------------------------------------------------------

    !// Skip new UBC for old H2O treatment
    if (LOLD_H2OTREATMENT) return

    !// In case H2O is not transported, return
    if (trsp_idx(114) .lt. 0) return

    if (trsp_idx(46).le.0 .or. trsp_idx(113).le.0) then
       write(6,'(a)') f90file//':'//subr// &
            ': Upper stratospheric boundary condition of'
       write(6,'(a)') '  H2O needs CH4 and H2 to be transported!'
       write(6,'(a,i5)') '  trsp_idx(46):',trsp_idx(46)
       write(6,'(a,i5)') '  trsp_idx(113):',trsp_idx(113)
       stop
    end if

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// Set H2 at model top = H2 at LPAR-1
          BTT(LPAR,trsp_idx(113),II,JJ) = BTT(LPAR-1,trsp_idx(113),II,JJ) &
                                    / AIRB(LPAR-1,II,JJ) * AIRB(LPAR,II,JJ)

          !// Set CH4 at model top = CH4 at LPAR-1
          BTT(LPAR,trsp_idx(46),II,JJ) = BTT(LPAR-1,trsp_idx(46),II,JJ) &
                                    / AIRB(LPAR-1,II,JJ) * AIRB(LPAR,II,JJ)

          !// Set H2O at model top = H2O at LPAR-1
          BTT(LPAR,trsp_idx(114),II,JJ) = BTT(LPAR-1,trsp_idx(114),II,JJ) &
                                    / AIRB(LPAR-1,II,JJ) * AIRB(LPAR,II,JJ)

       end do
    end do

    !// Update BTTBCK
    BTTBCK(LPAR,trsp_idx(114),:,:) = BTT(LPAR,trsp_idx(114),:,:)
    BTTBCK(LPAR,trsp_idx(113),:,:) = BTT(LPAR,trsp_idx(113),:,:)
    BTTBCK(LPAR,trsp_idx( 46),:,:) = BTT(LPAR,trsp_idx( 46),:,:)

    !// --------------------------------------------------------------------
  end subroutine strat_h2o_ubc2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine strat_h2o_max(message)
    !// --------------------------------------------------------------------
    !// Print max H2O in stratosphere.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    use cmn_ctm, only: STT, AIR
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: message
    !// Locals
    real(r8) :: m_h2o,max_h2o,max_h2o_lpar, st_h2o, st_h2o_L36
    integer :: I,J,L,mi,mj,ml
    !// --------------------------------------------------------------------
    if (LOLD_H2OTREATMENT) return

    max_h2o = 0._r8
    st_h2o = 0._r8
    do J = 1, JPAR
       do I = 1, IPAR
          do L = LMTROP(I,J)+1, LPAR-1

             st_h2o = st_h2o + STT(I,J,L,trsp_idx(114))

             m_h2o = STT(I,J,L,trsp_idx(114))/AIR(I,J,L)
             if (m_h2o .gt. max_h2o) then
                max_h2o = m_h2o
                mi=i
                mj=j
                ml=l
                if (max_h2o .gt. 1) then
                   print*,'weird h2o',m_h2o,i,j,l
                   print*,STT(I,J,L,trsp_idx(114)), AIR(I,J,L)
                   stop
                end if
             end if
          end do
       end do
    end do
    max_h2o_lpar = maxval(STT(:,:,LPAR,trsp_idx(114))/AIR(:,:,LPAR))
    st_h2o_L36 = sum(STT(:,:,36:LPAR,trsp_idx(114)))
    write(6,'(a,2es11.3,4i4)') &
          '* Max H2O strat/LPAR [ppmv]: '//trim(message)//': ', &
          max_h2o * M_AIR / 18._r8 * 1.e6_r8,&
          max_h2o_lpar * M_AIR / 18._r8 * 1.e6_r8, &
          mi,mj,ml,lmtrop(mi,mj)
    write(6,'(a,2es14.5)') '* Total H2O strat/L36[kg]: ', &
         st_h2o, st_h2o_L36
    !// --------------------------------------------------------------------
  end subroutine strat_h2o_max
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine zc_strh2o(ZC_LOCAL,II,JJ,MP,L_START,L_END,LPAR,TRACER_ID_MAX, &
       trsp_idx,DV)
    !// --------------------------------------------------------------------
    !// SPECIAL: H2O is put into 114 also when it is not transported, i.e.
    !//          when H2O is set from the sum of H2 in stratosphere.
    !//
    !// The conversion uses molecular weight (MW [g/mol]), volume (DV [m3]),
    !// and Avogadros number (Na):
    !// ZC(molec/cm3) = ZC(kg) * 1/(MW * kg/(1d3g) * Na / (DV * 1.d6cm3/m3)
    !//               = ZC(kg) * 1d3/MW * Na / DV / 1d6
    !//               = ZC(kg) * 1d-3/MW * Na / DV
    !//
    !// Amund Sovde, February 2011
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: LPAR,TRACER_ID_MAX
    integer, intent(in) :: II,JJ, MP, L_START,L_END
    integer, intent(in) :: trsp_idx(TRACER_ID_MAX)
    real(r8), intent(in) :: DV(LPAR)
    !// Input/output
    real(r8), intent(inout), dimension(TRACER_ID_MAX,LPAR) :: ZC_LOCAL

    !// Locals
    integer :: L
    real(r8)  :: RDUM
    !// --------------------------------------------------------------------

    !// If H2O is not transported, H2O treatment is from sum of H2,
    !// and we put the STR_H2O into ZC_LOCAL.
    if (trsp_idx(114) .gt. 0) return

    !// Convert from mass to concentration.
    RDUM = 1.e-3_r8 / 18._r8 * AVOGNR
    do L = L_START, L_END
       !// kg * 1d-3molec/kg / cm3 = molec/cm3
       ZC_LOCAL(114,L) = STR_H2O(L,II,JJ,MP) * RDUM / DV(L)
    end do


    !// --------------------------------------------------------------------
  end subroutine zc_strh2o
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine set_h2_eurohydros()
    !// --------------------------------------------------------------------
    !// Initialize H2 and H2O globally.
    !// H2 based on annual mean of the last Eurohydros simulation.
    !//
    !// Amund Sovde, December 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use cmn_size, only: LOSLOCSTRAT
    use cmn_ctm, only: AIR, STT, AREAXY
    use cmn_met, only: ZOFLE
    use utilities, only: get_free_fileid
    use cmn_oslo, only: XSTT
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8)  :: old_h2, rsum
    integer :: I,J,L,II,JJ,MP
    !// H2 from file
    character(len=80) :: FNAME
    integer :: ifnr, io_err, IR, JR, LR
    real(r4)  :: STTH2(IPAR,JPAR,LPAR)
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'set_h2_eurohydros'
    !// --------------------------------------------------------------------

    !// If no stratospheric chemistry, skip H2O
    if (.not. LOSLOCSTRAT) return

    !// Test if H2O, CH4 and H2 are included and how they are treated.
    if (.not.LOLD_H2OTREATMENT .and. trsp_idx(113).le.0) then
       write(6,'(a)') f90file//':'//subr// &
            ': H2-initialization problem!'
       write(6,'(a)') '  New H2O treatment needs H2 to be transported!'
       write(6,'(a,2i5)') '  Xtrsp_idx(113)/trsp_idx(113):',Xtrsp_idx(113),trsp_idx(113)
       stop
    else if (LOLD_H2OTREATMENT .and. Xtrsp_idx(113).le.0 ) then
       write(6,'(a)') f90file//':'//subr// &
            ': H2-initialization problem!'
       write(6,'(a)') '  Old H2O treatment needs H2 not to be transported'
       write(6,'(a,2i5)') '  Xtrsp_idx(113)/trsp_idx(113):',Xtrsp_idx(113),trsp_idx(113)
       stop
    end if


    !// Initialize H2 from EUROHYDROS annual mean
    write(*,'(a)') '* Initializing H2'
    !// Find non-used file number for input file
    ifnr = get_free_fileid()

    FNAME = 'Indata_CTM3/annualmean_h2.dta'
    open(ifnr,file=FNAME,status='old',form='unformatted',iostat=io_err)
    if (io_err .ne. 0) then
       print*,'* Could not open '//trim(fname)
       stop
    end if
    read(ifnr) IR,JR,LR
    if (IR.ne.IPAR.or.JR.ne.JPAR.or.LR.ne.LPAR) then
       print*,'* Resolution problems for '//trim(FNAME)
       print*,'  file:',IR,JR,LR
       print*,'  run: ',IPAR,JPAR,LPAR
       stop
    end if
    read(ifnr) STTH2
    close(ifnr)

    if (trsp_idx(113) .gt. 0) then
       old_h2 = sum(STT(:,:,:,trsp_idx(113)))
    else 
       old_h2 = sum(XSTT(:,Xtrsp_idx(113),:,:))
    end if

    !// Data is volume mixing ratio; convert to mass
    if (trsp_idx(113) .gt. 0) then
       STT(:,:,:,trsp_idx(113)) = real(STTH2(:,:,:), r8) &
                                  * AIR(:,:,:) * 2._r8 / M_AIR
       rsum = sum(STT(:,:,:,trsp_idx(113)))
    else if (Xtrsp_idx(113) .gt. 0) then
       do MP=1,MPBLK
          do J = MPBLKJB(MP), MPBLKJE(MP)
             JJ    = J - MPBLKJB(MP) + 1
             do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                do L = 1, LPAR
                   XSTT(L,Xtrsp_idx(113),I,J) = real(STTH2(I,J,L), r8) &
                                *AIR(I,J,L) * 2._r8 / M_AIR
                end do
             end do
          end do
       end do
       rsum = sum(XSTT(:,Xtrsp_idx(113),:,:))
    else
       write(6,'(a)') f90file//':'//subr//': How dit we get here???'
       stop
    end if

    write(*,'(a,es14.7,a,es14.7)') '* Old H2: ',old_h2,' New H2: ',rsum

    !// --------------------------------------------------------------------
  end subroutine set_h2_eurohydros
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine strat_h2o_init_clim()
    !// --------------------------------------------------------------------
    !// Initialize tropopause mixing ratio of H2O used as climatology.
    !// Default is to set a fixed mixing ratio, but could be set from
    !// model climatologi or the like.
    !//
    !// Amund Sovde, May 2013
    !// --------------------------------------------------------------------
    use cmn_size, only: LOSLOCSTRAT
    use ncutils, only: get_netcdf_var_1d, get_netcdf_var_3d
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=80)  :: infile   !// Name of netCDFfile
    integer            :: status, ncid, nDims, nVars, nAtts
    integer            :: nLon, nLat, nLev, nTime
    real(r8), allocatable, dimension(:) :: inLon, inLat, inTime
    !// --------------------------------------------------------------------

    !// If no stratospheric chemistry, skip H2O
    if (.not. LOSLOCSTRAT) return

    !// Skip if old H2O treatment
    if (LOLD_H2OTREATMENT) return

    if (.false.) then
       !// Read data from netcdf file
       infile = 'ctm3_h2o_tropopause_1998-2005.nc'
       write(*,'(a)') '* Reading H2O@TP: '//trim(infile)

       !// Check resolution (latitude/longitude/time)
       !// This routine allocates inLon/inLat/inTime
       call get_netcdf_var_1d( infile, 'lon',  inLon  )
       call get_netcdf_var_1d( infile, 'lat',  inLat  )
       call get_netcdf_var_1d( infile, 'time', inTime )

       nLon  = SIZE( inLon  )
       nLat  = SIZE( inLat  )
       nTime = SIZE( inTime )

       !// Deallocate all local variables
       IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
       IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
       IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)

       if (nLon.ne.IPAR .or. nLat.ne.JPAR) then
          write(*,'(a)') '* Data resolution != model resolution'
          write(*,'(a,2(i5))') '  Specified:',IPAR,JPAR
          write(*,'(a,2(i5))') '  On file:  ',nLon,nLat
          stop
       end if

       !// Get H2O [vmr]
       call get_netcdf_var_3d( infile, 'H2OtpCLIM', H2O_TP_CLIM, &
                               nlon,nlat,ntime )

       write(*,'(a,2es11.3)') '* Got climatology of H2O@TP',&
            minval(H2O_TP_CLIM),maxval(H2O_TP_CLIM)
    else
       !// Set fixed H2O in LS
       H2O_TP_CLIM(:,:,:) = 3.7e-6_r8
       write(*,'(a,2es11.3)') '* Setting fixed value as climatology of H2O@TP',&
            minval(H2O_TP_CLIM), maxval(H2O_TP_CLIM)
    end if


    !// --------------------------------------------------------------------
  end subroutine strat_h2o_init_clim
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module strat_h2o
!//=========================================================================

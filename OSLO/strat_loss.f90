!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, May 2015
!//=========================================================================
!// Stratospheric loss of species.
!//=========================================================================
module strat_loss
  !// ----------------------------------------------------------------------
  !// MODULE: strat_aerosols
  !// DESCRIPTION: Calculates loss of tracers in the stratosphere, i.e.
  !//              tracers that are calculated in the troposphere but not
  !//              in the stratosphere.
  !//
  !// Amund Sovde, June 2010
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  private
  public stratloss_oslo
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine stratloss_oslo(BTT,DTOPS,MP)
    !// --------------------------------------------------------------------
    !// To make sure you have control on the tracers in the stratosphere
    !// that are not handled by the stratospheric chemistry module.
    !//
    !// Introduces simple loss mechanisms of tropospheric components in the
    !// stratosphere. Stratospheric life times are assumed and used for
    !// e-folding losses.
    !//
    !// Note that this treatment is unnecessary for non-transported
    !// species, and will not be done.
    !//
    !// Handles:
    !//   Tropospheric chemistry
    !//   Stratospheric chemistry
    !//   Sulphur
    !//
    !// Amund Sovde, October 2008 - October 2009
    !// ------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LSULPHUR, LOSLOCSTRAT, LSOA
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_chem, only: TMASS
    use cmn_met, only: ZOFLE
    use cmn_oslo, only: LMTROP, trsp_idx, LVS2ADD2TC
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out parameters
    real(r8), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK):: BTT
    integer, intent(in) :: MP
    real(r8), intent(in)  :: DTOPS
    !// local parameters
    real(r8), parameter :: &
          ZTAUDAY   = 1.157407e-5_r8, &
          ZTAUWEEK  = 1.653439e-6_r8, &
          ZTAUMONTH = 3.858024e-7_r8, &
          ZTAUYEAR  = 3.170979e-8_r8, &
          ZTAUCH4   = 3.170979e-9_r8
    real(r8) :: LOSS, fall_in, fall_out, TAUPARTICLE, DZ, KEEP_FRAC
    integer :: I, J, L, N, II, JJ

    !// Tropospheric components lost in stratosphere
    integer, parameter :: CNR = 25
    integer, dimension(CNR), parameter :: COMP_LOST = &
         (/ 7,  8,  9, 10, 11, 12,     14,             18, 19, &
                                   26, 27, 28, 29, 30, 31, 32, &
           33, 34, 35, 36, 37, &
                       48, 49, 50, 51/)
    !// Tracer 13,15,16,17 are also lost when stratospheric chemistry
    !// is NOT calculated.
    integer, parameter :: CNR_ADD = 4
    integer, dimension(CNR_ADD), parameter :: &
         ADD_LOST = (/13, 15, 16, 17/)
    real(r8), dimension(LPAR,IDBLK,JDBLK) :: StoSO4
    integer :: LVSADD, TRNR
    !// --------------------------------------------------------------------

    !// Most of the tropospheric tracers have a lifetime of about one week
    !// in the stratosphere
    LOSS = -DTOPS * ZTAUWEEK

    if (LOSLOCSTRAT) then
       LVSADD = 0          !// Do not add levels when STRATCHEM is used
    else
       LVSADD = LVS2ADD2TC !// Add levels when STRATCHEM is not used
    end if

    !// These tracers are lost also when stratospheric chemistry is done
    KEEP_FRAC = exp(LOSS)
    do N = 1, CNR
      if (trsp_idx(COMP_LOST(N)).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              BTT(L,trsp_idx(COMP_LOST(N)),II,JJ) = &
                   BTT(L,trsp_idx(COMP_LOST(N)),II,JJ) * KEEP_FRAC
              !// Moments are scaled in oc_main at the bottom
            end do
          end do
        end do
      end if
    end do

    !// If only tropospheric chemistry, some additional tracers are lost
    if (.not. LOSLOCSTRAT) then
      do N = 1, CNR_ADD
        if (trsp_idx(ADD_LOST(N)).gt.0) then
          !// Losing transported species
          do J = MPBLKJB(MP), MPBLKJE(MP)
            JJ    = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP), MPBLKIE(MP)
              II    = I - MPBLKIB(MP) + 1
              do L = LMTROP(I,J)+1+LVSADD, LPAR
                BTT(L,trsp_idx(ADD_LOST(N)),II,JJ) = &
                     BTT(L,trsp_idx(ADD_LOST(N)),II,JJ) * KEEP_FRAC
              end do
            end do
          end do
        end if
      end do
    end if !// if (.not. LOSLOCSTRAT) then



    !// Section with explicit Loss
    !// ------------------------------------------------------------------
    !// PANX (tracer ID = 5)
    if (trsp_idx(5).gt.0) then
      !// Losing transported species
      do J = MPBLKJB(MP), MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = LMTROP(I,J)+1+LVSADD, LPAR
            BTT(L,trsp_idx(5),II,JJ) = BTT(L,trsp_idx(5),II,JJ) * KEEP_FRAC
          end do
        end do
      end do
    end if
      
    !// Isoprene (tracer ID = 20)
    KEEP_FRAC = exp(-1._r8)
    if (trsp_idx(20).gt.0) then
      !// Losing transported species
      do J = MPBLKJB(MP), MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          do L = LMTROP(I,J)+1+LVSADD, LPAR
            BTT(L,trsp_idx(20),II,JJ) = BTT(L,trsp_idx(20),II,JJ) * KEEP_FRAC
          end do
        end do
      end do
    end if


    !// When running only tropospheric chemistry
    !// ------------------------------------------------------------------
    if (.not. LOSLOCSTRAT) then
      !// CO (tracer ID = 6)
      KEEP_FRAC = exp(-DTOPS * ZTAUMONTH)
      if (trsp_idx(6).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              BTT(L,trsp_idx(6),II,JJ) = BTT(L,trsp_idx(6),II,JJ) * KEEP_FRAC
            end do
          end do
        end do
      end if
      !// Methane (tracer ID = 46)
      KEEP_FRAC = exp(-DTOPS * ZTAUCH4)
      if (trsp_idx(46).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              BTT(L,trsp_idx(46),II,JJ) = BTT(L,trsp_idx(46),II,JJ) * KEEP_FRAC
            end do
          end do
        end do
      end if

    end if !// if (.not. LOSLOCSTRAT) then


    !// Sulphur stratloss
    !// ------------------------------------------------------------------
    if (LSULPHUR) then

      StoSO4(:,:,:) = 0._r8
      !// DMS
      KEEP_FRAC = exp(-DTOPS*ZTAUWEEK)
      if (trsp_idx(71).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              StoSO4(L,II,JJ) = StoSO4(L,II,JJ) &
                   + BTT(L,trsp_idx(71),II,JJ) * (1._r8 - KEEP_FRAC) &
                     / TMASS(trsp_idx(71)) * TMASS(trsp_idx(73)) !// to mass SO4
              BTT(L,trsp_idx(71),II,JJ) = BTT(L,trsp_idx(71),II,JJ) * KEEP_FRAC
            end do
          end do
        end do
      end if
      !// SO2
      KEEP_FRAC = exp(-DTOPS*ZTAUWEEK)
      if (trsp_idx(72).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              StoSO4(L,II,JJ) = StoSO4(L,II,JJ) &
                   + BTT(L,trsp_idx(72),II,JJ) * (1._r8 - KEEP_FRAC) &
                     / TMASS(trsp_idx(72)) * TMASS(trsp_idx(73)) !// to mass SO4
              BTT(L,trsp_idx(72),II,JJ) = BTT(L,trsp_idx(72),II,JJ) * KEEP_FRAC
            end do
          end do
        end do
      end if
      !// H2S
      KEEP_FRAC = exp(-DTOPS * ZTAUWEEK)
      if (trsp_idx(74).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
              StoSO4(L,II,JJ) = StoSO4(L,II,JJ) &
                   + BTT(L,trsp_idx(74),II,JJ) * (1._r8 - KEEP_FRAC) &
                     / TMASS(trsp_idx(74)) * TMASS(trsp_idx(73)) !// to mass SO4
              BTT(L,trsp_idx(74),II,JJ) = BTT(L,trsp_idx(74),II,JJ) * KEEP_FRAC
            end do
          end do
        end do
      end if

      !// SO4 & MSA, assume fall velocity of 2,5 cm^-hr
      if (trsp_idx(73).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            !// Loop from top
            fall_in = 0._r8
            do L = LPAR, LMTROP(I,J)+1+LVSADD, -1
              DZ = ZOFLE(L+1,I,J) - ZOFLE(L,I,J)
              TAUPARTICLE = 3600._r8 * DZ / 0.025_r8 ! lifetime seconds
              LOSS = DTOPS / TAUPARTICLE

              !// Falling down
              fall_out = BTT(L,trsp_idx(73),II,JJ) * LOSS

              BTT(L,trsp_idx(73),II,JJ) = &
                   BTT(L,trsp_idx(73),II,JJ) * (1._r8 - LOSS) + fall_in
              !// Coming in to the layer below
              fall_in = fall_out
            end do
          end do
        end do
        !// Add DMS, SO2 and H2S to H2SO4 [kg]
        do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            do L = LMTROP(I,J)+1+LVSADD, LPAR
                 BTT(L,trsp_idx(73),II,JJ) = BTT(L,trsp_idx(73),II,JJ) &
                                             + StoSO4(L,II,JJ)
            end do
          end do
        end do
      end if

      if (trsp_idx(75).gt.0) then
        !// Losing transported species
        do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP),MPBLKIE(MP)
            II    = I - MPBLKIB(MP) + 1
            !// Loop from top
            fall_in = 0._r8
            do L = LPAR,LMTROP(I,J)+1+LVSADD, -1
              DZ = ZOFLE(L+1,I,J) - ZOFLE(L,I,J)
              TAUPARTICLE = 3600._r8 * DZ / 0.025_r8 ! lifetime seconds
              LOSS = DTOPS / TAUPARTICLE

              !// Falling down
              fall_out = BTT(L,trsp_idx(75),II,JJ) * LOSS

              BTT(L,trsp_idx(75),II,JJ) = BTT(L,trsp_idx(75),II,JJ) &
                   * max(0._r8, (1._r8 - LOSS)) + fall_in
              !// Coming in to the layer below
              fall_in = fall_out
            end do
          end do
        end do
      end if

    end if !// if (LSULPHUR) then
    !// ------------------------------------------------------------------


    !// NITRATE stratloss
    !// ------------------------------------------------------------------
    !// Nitrate stratloss is treated in nitrate routines.


    !// BCOC stratloss
    !// ------------------------------------------------------------------
    !// June 2012: BCOC stratloss treated in BCOC routines.


    !// SOA stratloss
    !// ------------------------------------------------------------------
    !// SOA should only be lost when products are not treated in the CTM,
    !// e.g. CO2. Sedimentation could also be a loss process, but not as
    !// stratloss.
    !// June 2012:
    !// CTM2 used 1 week e-folding lifetime, which is very short. Up it to
    !// one month.
    !// When running stratospheric chemistry, we do SOA separation there as
    !// well.
    !//
    !// Stratloss applies for SOAGASes, but when not running stratospheric
    !// chemistry, we do stratloss for aerosols also.
    !// Aerosols are subject to gravitational settling, treated here.
    if (LSOA) then
      do N = 150, 193
        TRNR = trsp_idx(N)
        !// Go to next N if tracer not included
        if (TRNR .le. 0) cycle

        !// Losing transported species
        select case(N)

        case(150:164,180:181,184:187,192:193)
          !// SOA gases
          KEEP_FRAC = exp(-DTOPS * ZTAUMONTH)
          do J = MPBLKJB(MP), MPBLKJE(MP)
            JJ    = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP), MPBLKIE(MP)
              II    = I - MPBLKIB(MP) + 1
              do L = LMTROP(I,J)+1+LVSADD, LPAR
                BTT(L,TRNR,II,JJ) = BTT(L,TRNR,II,JJ) * KEEP_FRAC
              end do
            end do
          end do

        case default
          !// SOA aerosols
          !// First do loss, then gravitational settling
          KEEP_FRAC = Exp(-DTOPS*ZTAUMONTH)
          do J = MPBLKJB(MP),MPBLKJE(MP)
            JJ    = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP),MPBLKIE(MP)
              II    = I - MPBLKIB(MP) + 1
              do L = LMTROP(I,J)+1+LVSADD, LPAR
                BTT(L,TRNR,II,JJ) = BTT(L,TRNR,II,JJ) * KEEP_FRAC
              end do
            end do
          end do

          !// Only tropospheric chemistry: Also lose aerosols
          !// When stratospheric chemistry, include gravitational
          !// settling of aerosols.
          !if (LOSLOCSTRAT) then
          if (.false.) then
            do J = MPBLKJB(MP), MPBLKJE(MP)
              JJ    = J - MPBLKJB(MP) + 1
              do I = MPBLKIB(MP), MPBLKIE(MP)
                II    = I - MPBLKIB(MP) + 1
                !// Loop downwards
                fall_in = 0._r8
                do L = LPAR, LMTROP(I,J)+1+LVSADD, -1
                  !// Add from above
                  BTT(L,TRNR,II,JJ) = BTT(L,TRNR,II,JJ) + fall_in
                  !// Gravitational settling
                  !// Assume fall speed 2.5cm/hr = 6.944d-6m/s, as for SO4.
                  !// 1/6.944d-6m/s = 144000s/m, so we get a life time of
                  !// DZ / 6.944d-6 = DZ * 144000 (in seconds)
                  LOSS = 1._r8 / ((ZOFLE(L+1,I,J) - ZOFLE(L,I,J))*144000._r8)
                  !// Falling down during this time step
                  fall_out = BTT(L,TRNR,II,JJ) * LOSS * DTOPS
                  BTT(L,TRNR,II,JJ) = BTT(L,TRNR,II,JJ) &
                       * max(0._r8, (1._r8 - LOSS * DTOPS))
                  !// Input to layer below
                  fall_in = fall_out
                end do
              end do
            end do
          end if

        end select
      end do


      !// Lose all SOA primary gases
      KEEP_FRAC = Exp(-DTOPS*ZTAUMONTH)
      do N = 280, 291
        TRNR = trsp_idx(N)
        if (TRNR.gt.0) then
          !// Losing transported species
          do J = MPBLKJB(MP), MPBLKJE(MP)
            JJ    = J - MPBLKJB(MP) + 1
            do I = MPBLKIB(MP), MPBLKIE(MP)
              II    = I - MPBLKIB(MP) + 1
              do L = LMTROP(I,J)+1+LVSADD, LPAR
                BTT(L,TRNR,II,JJ) = BTT(L,TRNR,II,JJ) * KEEP_FRAC
              end do
            end do
          end do
        end if
      end do
    end if !// if (LSOA) then

    !// --------------------------------------------------------------------
  end subroutine stratloss_oslo
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module strat_loss
!//=========================================================================

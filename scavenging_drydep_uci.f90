!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, June 2015
!//=========================================================================
!// UCI dry deposition routines.
!//=========================================================================
module scavenging_drydep_uci
  !// ----------------------------------------------------------------------
  !// MODULE: scavenging_drydep_uci
  !// DESCRIPTION: Treat dry deposition as a separate process.
  !//              Converted from UCI p-scav.f and modified for CTM3.
  !//
  !// ---(p-scav.f)    p-code 5.6 (Neu, Tang, Prather Aug 2008)
  !// Version: qcode_56d; 20090318
  !//     CTM3: modified VDEP
  !//           DRYDEP: Added LEMISDEP_INCHEM switch
  !//
  !// Contains
  !//   subroutine DRYSET_CTM3
  !//   subroutine DRYDEP
  !//   subroutine DEP_V2
  !//
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  private
  public DRYSET_CTM3, DRYDEP
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine DRYSET_CTM3(INFILE)
    !// --------------------------------------------------------------------
    !// Sets UCI dry dep velocities.
    !// These deposition rates will not be used in CTM3; CTM3 uses a more
    !// sophisticated parameterisation to set VDEP.
    !//
    !// Modified for CTM3.
    !// Ole Amund Sovde, June 2012
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, NPAR
    use cmn_ctm, only: NTM, YDGRD, PLAND, NSCX
    use cmn_chem, only: TNAME
    use cmn_sfc, only: VDEPIN, VDEP
    use utilities, only: get_free_fileid
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) ::  INFILE
    !// Locals
    character(len=80) :: TITLE, SCAV_FMT
    character(len=10) :: XTRACR
    logical :: LXX_AER
    real(r8) :: XXDEP1, XXDEP2, XXDEP3, XXLIFE(9)
    integer :: N_TRACERS, I, J, N, NN, LSURF, FUN
    !// --------------------------------------------------------------------

    VDEPIN(:,:) = 0._r8

    write(*,'(a71)') '--------------------------------------------'// &
         '---------------------------'

    !// Open dry scavenging data file
      write(6,*) ' open file: '//trim(INFILE)
    fun = get_free_fileid()
    open(fun,file=INFILE,status='old',form='formatted')
    read(fun,'(a60)') TITLE
      write(6,*) trim(TITLE)
    read(fun,'(a)') SCAV_FMT
    read(fun,SCAV_FMT) N_TRACERS
    do N = 1, N_TRACERS
      read(fun,SCAV_FMT) I,XTRACR,XXDEP1,XXDEP2,XXDEP3
      !// Search for current tracers:
      do NN = 1, NTM
        if (TNAME(NN) .eq. XTRACR) then
          write(6,*) ' dry scavenging match: ', XTRACR,I,NN
          VDEPIN(NN,1) = XXDEP1 * 1.e-2_r8      ! v-dep in m/s
          VDEPIN(NN,2) = XXDEP2 * 1.e-2_r8
          VDEPIN(NN,3) = XXDEP3 * 1.e-2_r8
        end if
      end do
    end do
    close(fun)

    !// Map deposition velocities based on land/ocean/cryo surface
    do J = 1, JPAR
      do I = 1, IPAR
        if (abs(YDGRD(J)) .lt. 60._r8) then
          if (PLAND(I,J) .gt. 0.5_r8) then
            LSURF=1
          else
            LSURF=2
          end if
        else
          LSURF=3
        end if
        do N = 1, NTM
          !// CTM3: Changed order of indices
          VDEP(N,I,J) = VDEPIN(N,LSURF)
        end do
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine DRYSET_CTM3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DRYDEP(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
       UTTAU,DTDRYD,MP)
    !// --------------------------------------------------------------------
    !// all of this needs to be updated 
    !// this method calculated the depostion velocity VDEP for each species
    !// and then just removes tracer & moments in the lowest part of LAYER=1
    !// using VDEP*DTDRYD to determine how much ***only layer L=1 affected!
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LEMISDEP_INCHEM
    use cmn_ctm, only: NTM, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, JMON, NDPX
    use cmn_met, only: ZOFLE
    use cmn_sfc, only: VDEP
    use utilities, only: ctmExitC
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) ::  UTTAU, DTDRYD
    integer, intent(in) ::  MP
    !// Input/Output
    real(r8), intent(inout) ::  BTT(LPAR,NPAR,IDBLK,JDBLK)
    real(rMom), intent(inout), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    !// Locals
    real(r8) :: ETT, EZT, EZZ, EXZ, EYZ, EXT, EXX, EYT, EYY, EXY
    real(r8) :: XL, XDEP, DELZ, LSTAU
    integer :: I, J, II, JJ, M, N
    !// --------------------------------------------------------------------

    !//---collect column variables to calculate height of layers
    !//---R = 287.*(1-Q) + 461.5*Q -- assume 0.5% bl w.v.==> R = 288.
    !//---delta-z = dln(P) * R * T / g   where R/g = 288/9.81 = 29.36
    !// --------------------------------------------------------------------
    !do L = 1,LM+1
    !   POFL(L) = ETAA(L) + ETAB(L)*P(I,J)
    !end do
    !ZOFL(1) = 0._r8
    !do L = 1,LM
    !   DELZ = log(POFL(L)/POFL(L+1))*T(I,J,L)*29.36_r8
    !   ZOFL(L+1) = ZOFL(L) + DELZ
    !end do
    !DELZ = ZOFL(2)
    !// --------------------------------------------------------------------


    !// Swich for Oslo chemistry treatment of emissions. If they are to 
    !// be treated as production terms in chemistry, they will be fetched
    !// in oc_main.f, so we just return from here.
    if (LEMISDEP_INCHEM) return


    M = JMON        ! may need to know month
    LSTAU = UTTAU   ! may need Local Solar TAU for each I
    !// begin super loop over OpenMP IJ-block
    do J = MPBLKJB(MP),MPBLKJE(MP)
      JJ   = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP),MPBLKIE(MP)
        II   = I - MPBLKIB(MP) + 1

        !// NB dep velocity cannot reach beyond L=1, height = DELZ
        !// at 1 cm/s (fast) ==> 36 m in 1 hr, should be OK
        DELZ = ZOFLE(2,I,J) - ZOFLE(1,I,J)

        do N = 1, NTM
          if (NDPX .eq.0) then
            !// Fraction to be removed
            XDEP = VDEP(N,I,J) * DTDRYD / DELZ
            XDEP = min(0.99_r8, XDEP)
          else
            !// updated dry dep code should be function of veg, LAI, month...
            call ctmExitC('>>>>cannot run NDPX>0 dry dep with this code')
          end if
          if (XDEP .gt. 1.e-8_r8) then
            ETT   = BTT(1,N,II,JJ)
            EZT   = BZT(1,N,II,JJ)
            EZZ   = BZZ(1,N,II,JJ)
            EXZ   = BXZ(1,N,II,JJ)
            EYZ   = BYZ(1,N,II,JJ)
            EXT   = BXT(1,N,II,JJ)
            EXX   = BXX(1,N,II,JJ)
            EYT   = BYT(1,N,II,JJ)
            EYY   = BYY(1,N,II,JJ)
            EXY   = BXY(1,N,II,JJ)

            call DEP_V2(ETT,EZT,EZZ,EXZ,EYZ,EXT,EXX,EYT,EYY,EXY,XDEP,XL)

            BTT(1,N,II,JJ) = ETT
            BZT(1,N,II,JJ) = EZT
            BZZ(1,N,II,JJ) = EZZ
            BXZ(1,N,II,JJ) = EXZ
            BYZ(1,N,II,JJ) = EYZ
            BXT(1,N,II,JJ) = EXT
            BXX(1,N,II,JJ) = EXX
            BYT(1,N,II,JJ) = EYT
            BYY(1,N,II,JJ) = EYY
            BXY(1,N,II,JJ) = EXY
          end if
        end do

      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine DRYDEP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine DEP_V2(ETT,EZT,EZZ,EXZ,EYZ,EXT,EXX,EYT,EYY,EXY,XDEP,XL)
    !// --------------------------------------------------------------------
    !// Remove tracer from lowest XDEP (fraction) of layer; total loss (XL)
    !// correct moments to reflect deposition.    Oliver(27/7/99)+MJP(7/02)
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in)    :: XDEP
    !// Input/Output
    real(r8), intent(inout) :: ETT,EZT,EZZ,EXZ,EYZ,EXT,EXX,EYT,EYY,EXY
    !// Output
    real(r8), intent(out)   :: XL
    !// Locals
    real(r8) :: EM, EPS, EPS1, EPS1Q, ETT0
    !// --------------------------------------------------------------------

    !// calling qlimit to restrain Z moments (mass EM not needed for Lim=2)
    EM = 0._r8
    call QLIMIT2 (ETT,EZT,EZZ,EXZ,EYZ)

    !// compute moments of tracer in fraction of layer not deposited
    ETT0  = ETT
    EPS   = min(1._r8, max(0._r8, XDEP))
    EPS1  = 1._r8 - EPS
    EPS1Q = EPS1 * EPS1

    ETT = EPS1 * (ETT + EPS * (EZT - (EPS1 - EPS) * EZZ))
    EZT = EPS1Q * (EZT + 3._r8 * EPS * EZZ)
    EZZ = EPS1 * EPS1Q * EZZ
    EXT = EPS1 * (EXT + EPS * EXZ)
    EYT = EPS1 * (EYT + EPS * EYZ)
    EXX = EPS1 * EXX
    EYY = EPS1 * EYY
    EXY = EPS1 * EXY
    EXZ = EPS1Q * EXZ
    EYZ = EPS1Q * EYZ

    !// put Z-moments back into full layer with bottom of layer zeroed
    XL   = ETT - ETT0
    EZZ  = EPS1Q * EZZ+ 5._r8 * EPS * (EPS1 * EZT - (EPS1 - EPS) * ETT)
    EZT  = EPS1 * EZT + 3._r8 * EPS * ETT
    EXZ  = EPS1 * EXZ + 3._r8 * EPS * EXT
    EYZ  = EPS1 * EYZ + 3._r8 * EPS * EYT

    !// --------------------------------------------------------------------
  end subroutine DEP_V2
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module scavenging_drydep_uci
!//=========================================================================

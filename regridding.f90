!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routines for regridding datasets,
!//=========================================================================
module regridding
  !//-----------------------------------------------------------------------
  !// MODULE: regridding
  !// DESCRIPTION: Routines for regridding datasets.
  !//
  !// Contains
  !//   subroutine E_GRID
  !//   subroutine E_GRID_Y
  !//   subroutine TRUNG8
  !//   subroutine TRUNG4
  !//   subroutine Regrid_Column_Weights
  !//
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'regridding.f90'
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains



  !//-----------------------------------------------------------------------
  subroutine E_GRID(EBOX,XBEDG,YBEDG,IG,JG, EDXY,XDEDG,YDEDG,IM,JM,NEDXY)
    !//---------------------------------------------------------------------
    !// generic regridder of 
    !//        emissions EBOX(1:IG,1:JG) 
    !//        on regular grid XBEDG(1:IG+1),YBEDG(1:JG+1) 
    !//
    !// onto CTM grid   EDXY(1:IM,1:JM)
    !//        EBOX = emissions or absolute quantity (kg or kg/yr, not kg/m2)
    !// IF original data is %-coverage or kg/m2, THEN must be pre-mult by area
    !//
    !//    X's Y's are the grid edges in degrees
    !//    note that the Longitude (X) edges can jump by 360 deg
    !//    XEDG(1,IM+1) is created = CTM longitude-grid edges
    !//                              (radians, monotonic)
    !//
    !//    EDXY(1:IM, 1:JM, 1:6) = emissions partitioned into box (I,J)
    !//                            + moments
    !//           EDXY(:,:,1) = total emissions
    !//           EDXY(:,:,2) = Sx  moment of emissions
    !//           EDXY(:,:,3) = Sxx moment of emissions
    !//           EDXY(:,:,4) = Sy  moment of emissions
    !//           EDXY(:,:,5) = Syy moment of emissions
    !//           EDXY(:,:,6) = Sxy moment of emissions
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: CPI, C2PI, CPI180
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IG,JG, IM,JM, NEDXY
    real(r8), intent(in) :: EBOX(IG,JG),XBEDG(IG+1),YBEDG(JG+1)
    real(r8), intent(in) :: XDEDG(IM+1),YDEDG(JM+1)
    !// Output
    real(r8), intent(out) :: EDXY(IM,JM,NEDXY)
    !// Locals
    real(r8) :: Y1,Y2,SY1,SY2,DYBOX,YG,  YB1,YB2,SYB1,SYB2
    real(r8) :: B_GYA(JG,JM), B_GYY(JG,JM)
    real(r8) :: B_GXA(IG,IM), B_GXX(IG,IM)
    real(r8) :: KX0,KXX,KY0,KYY,KXY, EMKIKJ
    real(r8) :: SYG(JM+1), XEDG(IM+1)
    real(r8) :: X1,X2,DXBOX,XG,  XB1,XB2
    integer :: J,J1,J2,I,I1,I2, KI,KJ, I2M
    !//---------------------------------------------------------------------

    XEDG(1) = CPI180 * XDEDG(1)
    do I = 2, IM+1
       if (XDEDG(I) .gt. 0._r8) then
          XEDG(I) = CPI180 * XDEDG(I)
       else
          XEDG(I) = CPI180 * (XDEDG(I) + 360._r8)
       end if
    end do

    !// Initialize
    B_GXA(:,:) = 0._r8
    B_GXX(:,:) = 0._r8
    B_GYA(:,:) = 0._r8
    B_GYY(:,:) = 0._r8

    !// for relative area use sin(LATITUDE) x LONGITUDE (degrees)
    do J = 1, JM+1
       SYG(J) = sin(YDEDG(J) * CPI180)
    end do

    !// setup the LATITUDE mapping and weights for the emission boxes
    do KJ = 1, JG
       YB1 = YBEDG(KJ)
       YB2 = YBEDG(KJ+1)
       SYB1 = sin(CPI180 * YB1)
       SYB2 = sin(CPI180 * YB2)
       J1 = JM !// Initialize with J1=JM
       do J = 2, JM
          if (YDEDG(J) .gt. YB1) then
             J1 = J - 1
             exit !// exit J-loop (goto 12)
          end if
       end do
       !//J1 = JM
       !//12   continue

       J2 = JM !// Initialize with J2=JM
       do J = J1+1, JM+1
          if (YDEDG(J) .ge. YB2) then
             J2 = J - 1
             exit !// exit J-loop (goto 14)
          end if
       end do
       !//J2 = JM
       !//14   continue

       !// The emission box KJ=1:JG is contained in CTM grid boxes J1 to J2
       !// save the indices and fractions and weightings
       do J = J1, J2
          Y1 = max(YB1, YDEDG(J))
          Y2 = min(YB2, YDEDG(J+1))
          SY1 = sin(CPI180 * Y1)
          SY2 = sin(CPI180 * Y2)
          DYBOX = (SY2 - SY1) / (SYB2 - SYB1) ! fractional area of emission box
          !// mean value of Y in grid box
          YG = (0.5_r8 * (SY1 + SY2) - SYG(J)) / (SYG(J+1) - SYG(J))

          B_GYA(KJ,J) = DYBOX
          B_GYY(KJ,J) = YG

       end do
    end do

    !// Now setup the LONGITUDE mapping and weights for the emission boxes.
    !// Require XEDG(I=1:IM) to be monotonic, 
    !// but XBEDG(KI=1:IG) may have shift at dateline
    do KI = 1, IG
       XB1 = XBEDG(KI) * CPI180
       XB2 = XBEDG(KI+1) * CPI180
       if (XB1 .gt. XB2) XB2 = XB2 + C2PI
       if (XB1 .gt. XB2) then
          !// ERROR
          write(6,'(a)') 'regrid.f90:E_GRID: err#1 emis. X-grid'
          stop
       end if
       if (XB1 .lt. XEDG(1)) then
          XB1 = XB1 + C2PI
          XB2 = XB2 + C2PI
          if (XB1 .lt. XEDG(1)) then
             write(6,'(I4,6F9.3)') KI,XBEDG(KI)*CPI180,XBEDG(KI+1)*CPI180, &
                  XB1,XB2,XEDG(1)
             !// ERROR
             write(6,'(a)') 'regrid.f90:E_GRID: err#2 emis. X-grid'
             stop
          end if
       end if
       if (XB1 .gt. XEDG(IM+1)) then
          XB1 = XB1 - C2PI
          XB2 = XB2 - C2PI
          if (XB1 .gt. XEDG(IM+1)) then
             !// ERROR
             write(6,'(a)') 'regrid.f90:E_GRID: err#3 emis. X-grid'
             stop
          end if
       end if

       I1 = IM !// Initialise with I1=IM
       do I = 2, IM
          if (XEDG(I) .gt. XB1) then
             I1 = I - 1
             exit !// exit loop (goto 16)
          end if
       end do
       !//I1 = IM
       !//16   continue

       I2 = 0 !// Initialize with no value
       do I = I1+1, IM+1
          if (XEDG(I) .ge. XB2) then
             I2 = I - 1
             exit !// Exit this loop (goto 18)
          end if
       end do
       !// If I2 not found, allow for XB2 to extend (wrap) beyond XEDG(IM+1)
       if (I2 .eq. 0) then
          do I = 2, IM
             if ((XEDG(I) + C2PI) .ge. XB2) then
                I2 = IM + I-1
                exit !// Exit this loop (also goto 18)
             end if
          end do
       end if
       !// If I2 still not found, it is 2*IM
       if (I2 .eq. 0) I2 = IM + IM
       !//I2 = IM + IM
       !//18   continue

       !// The emission box K=1:IG is contained in CTM grid boxes I1 to I2
       !// but the I2 may be > IM to wrap around the Longitude circle.
       !// save the indices and fractions and weightings
       I2M = min(I2, IM)
       do I = I1, I2M
          X1 = max(XB1, XEDG(I))
          X2 = min(XB2, XEDG(I+1))
          DXBOX = (X2 - X1) / (XB2 - XB1) ! fractional area of emission box
          XG = (0.5_r8 * (X1 + X2) - XEDG(I)) &
                       / (XEDG(I+1) - XEDG(I)) ! mean X value 

          B_GXA(KI,I) = DXBOX
          B_GXX(KI,I) = XG
       end do

       !// while XB1 is constrained ot lie within 1:IM, XB2 may wrap around:
       IF (I2 .gt. IM) then
          I2M = I2 - IM
          do I = 1, I2M
             X1 = max(XB1, XEDG(I) + C2PI)
             X2 = min(XB2, XEDG(I+1) + C2PI)
             DXBOX = (X2 - X1) / (XB2 - XB1)
             XG = (0.5_r8 * (X1 + X2) - XEDG(I) - C2PI) &
                          / (XEDG(I+1) - XEDG(I))

             B_GXA(KI,I) = DXBOX
             B_GXX(KI,I) = XG
          end do
       end if
    end do !// do KI = 1, IG

    !// sum contributions (& moments) to CTM grid (I,J) from emission
    !// boxes (KI,KJ)
    !// SOM XY moments:
    !//     k00 = 1.d0
    !//     kX0 = XG - 0.5d0
    !//     kXX = XG*(XG - 1.d0) + C16TH
    !//     kY0 = YG - 0.5d0
    !//     kYY = YG*(YG - 1.d0) + C16TH
    !//     kXY = (XG - 0.5d0)*(YG - 0.5d0)
    !// for Kronecker-delat approx:  f(x,y) = DXBOX*DYBOX*EMI-BOX * k__(XG,YG)
    !//     So = INTEG[ f(x,y) ] x=0:1, y=0:1
    !//     Sx = 6 * INTEG[ f(x,y) * kX(x) ]
    !//     Sxy = 36 * INTEG[ f(x,y) * kXY(x,y) ]
    !//     Sxx = 30 * INTEG[ f(x,y) * kXX(x) ]
    !// adjust K_s to get So,Sx,Sxx,Sy,Syy,Sxy  (1:6)
    !//     K00 = 1.d0
    !//     KX0 = 6.d0*XG - 3.0d0
    !//     KXX = 30.d0*XG*(XG - 1.d0) + 5.d0
    !//     KY0 = 6.d0*YG - 3.0d0
    !//     KYY = 30.d0*YG*(YG - 1.d0) + 5.d0
    !//     KXY = KX0 * KY0

    if (NEDXY .eq.  6) then
      do J = 1, JM
        do I = 1, IM
          EDXY(I,J,:) = 0._r8
          do KJ = 1, JG
            if (B_GYA(KJ,J) .gt. 1.e-5_r8) then
              do KI = 1, IG
                if (B_GXA(KI,I) .gt. 1.e-5_r8) then
                  EMKIKJ = EBOX(KI,KJ) * B_GXA(KI,I) * B_GYA(KJ,J)
                  XG = B_GXX(KI,I)
                  YG = B_GYY(KJ,J)
                  KX0 = 6._r8 * XG - 3._r8
                  KXX = 30._r8 * XG * (XG - 1._r8) + 5._r8
                  KY0 = 6._r8 * YG - 3._r8
                  KYY = 30._r8 * YG * (YG - 1._r8) + 5._r8
                  KXY = KX0 * KY0
                  EDXY(I,J,1) = EDXY(I,J,1) + EMKIKJ
                  EDXY(I,J,2) = EDXY(I,J,2) + EMKIKJ * KX0
                  EDXY(I,J,3) = EDXY(I,J,3) + EMKIKJ * KXX
                  EDXY(I,J,4) = EDXY(I,J,4) + EMKIKJ * KY0
                  EDXY(I,J,5) = EDXY(I,J,5) + EMKIKJ * KYY
                  EDXY(I,J,6) = EDXY(I,J,6) + EMKIKJ * KXY
                end if
              end do
            end if
          end do
        end do
      end do !// do J = 1, JM

      !// 1st variable - EDXY(I,J,1) has real units of EBOX (kg or m2)
      !// 2nd-6th variables are the XY-moments and dimensionless.
      do J = 1, JM
        do I = 1, IM
          if (EDXY(I,J,1) .gt. 1.e-30_r8) then
            EDXY(I,J,2) = EDXY(I,J,2) / EDXY(I,J,1)
            EDXY(I,J,3) = EDXY(I,J,3) / EDXY(I,J,1)
            EDXY(I,J,4) = EDXY(I,J,4) / EDXY(I,J,1)
            EDXY(I,J,5) = EDXY(I,J,5) / EDXY(I,J,1)
            EDXY(I,J,6) = EDXY(I,J,6) / EDXY(I,J,1)
          end if
        end do
      end do
    else
      !// Do not consider moments (NEDXY == 1)
      do J = 1, JM
        do I = 1, IM
          EDXY(I,J,1) = 0._r8
          do KJ = 1,JG
            if (B_GYA(KJ,J) .gt. 1.e-5_r8) then
              do KI = 1,IG
                if (B_GXA(KI,I) .gt. 1.e-5_r8) then
                  EMKIKJ = EBOX(KI,KJ) * B_GXA(KI,I) * B_GYA(KJ,J)
                  EDXY(I,J,1) = EDXY(I,J,1) + EMKIKJ
                end if
              end do
            end if
          end do
        end do
      end do !// do J = 1, JM
    end if !// if (NEDXY .eq.  6) then

    !// Sum check for complete overlaps (done in calling program)
    !SUMIJ = 0._r8
    !do J = 1, JM
    !   do I = 1, IM
    !      SUMIJ = SUMIJ + EDXY(I,J,1)
    !   end do
    !end do
    !SUMK  = 0._r8
    !do KJ = 1, JG
    !   do KI = 1, IG
    !      SUMK = SUMK + EBOX(KI,KJ)
    !   end do
    !end do
    !if ( abs(SUMIJ - SUMK) .gt. 1.d-6 * SUMK) then
    !   write (6,'(a,1p,e12.5,a,e12.5)') &
    !        ' emissions mis-mapping: ',SUMK,' to ',SUMIJ
    !end if

    !//---------------------------------------------------------------------
  end subroutine E_GRID
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine E_GRID_Y(EBOX,YBEDG,JG, EDY,YDEDG,JM)
    !// --------------------------------------------------------------------
    !// generic regridder of 
    !//    emissions EBOX(1:JG) 
    !// on regular grid YBEDG(1:JG+1) 
    !// onto CTM grid   EDY(1:JM)
    !//
    !// EBOX = emissions or absolute quantity (kg or kg/yr, not kg/m2)
    !// IF original data is %-coverage or kg/m2, THEN must be pre-mult by area
    !//
    !// Y's are the grid edges in degrees
    !//     note that the Longitude (X) edges can jump by 360 deg
    !//
    !// EDY(1:JM) = emissions partitioned into box (J)
    !//
    !// Based on UCI E_GRID, but only with Y-component
    !// Ole Amund Sovde, December 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: CPI, C2PI, CPI180
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: JG, JM
    real(r8), intent(in) :: EBOX(JG),YBEDG(JG+1)
    real(r8), intent(in) :: YDEDG(JM+1)
    !// Output
    real(r8), intent(out) :: EDY(JM)
    !// Locals
    real(r8) :: Y1,Y2,SY1,SY2,DYBOX,YG,  YB1,YB2,SYB1,SYB2
    real(r8) :: B_GYA(JG,JM)
    real(r8) :: EMKIKJ
    real(r8) :: SYG(JM+1)
    integer :: J,J1,J2,KJ
    !// --------------------------------------------------------------------

    !// Fraction of JG inside each JM
    B_GYA(:,:) = 0._r8

    !// For relative area use sin(LATITUDE) x LONGITUDE (degrees)
    do J = 1,JM+1
       SYG(J) = sin(YDEDG(J) * CPI180)
    end do

    !// Setup the LATITUDE mapping and weights for the emission boxes
    do KJ = 1, JG
       YB1 = YBEDG(KJ)
       YB2 = YBEDG(KJ+1)
       SYB1 = sin(CPI180 * YB1)
       SYB2 = sin(CPI180 * YB2)

       !// Have removed original goto statements. If there was no exit
       !// of loop, J1 was set to JM, so let us initialize it with that.
       J1 = JM !// In case no exit of loop.
       do J = 2, JM
          if (YDEDG(J) .gt. YB1) then
             J1 = J - 1
             exit !// exit J-loop (goto 12)
          end if
       end do
       !J1 = JM
       !12 continue

       !// Same for J2
       J2 = JM
       do J = J1+1, JM+1
          if (YDEDG(J) .ge. YB2) then
             J2 = J - 1
             exit !// exit J-loop (goto 14)
          end if
       end do
       !J2 = JM
       !14 continue

       !// The emission box KJ=1:JG is contained in CTM grid boxes J1 to J2
       !// save the indices and fractions and weightings
       do J = J1, J2
          Y1 = max(YB1, YDEDG(J))
          Y2 = min(YB2, YDEDG(J+1))
          SY1 = sin(CPI180*Y1)
          SY2 = sin(CPI180*Y2)
          DYBOX = (SY2-SY1)/(SYB2-SYB1)     ! fractional area of emission box
          !// mean value of Y in grid box
          YG = (0.5_r8 * (SY1 + SY2) - SYG(J))/(SYG(J+1) - SYG(J))

          B_GYA(KJ,J) = DYBOX
       end do

    end do !// do KJ = 1,JG


    !// Sum contributions to CTM grid (J) from emission boxes (KJ)
    do J = 1, JM
       EDY(J) = 0._r8
       do KJ = 1, JG
          if (B_GYA(KJ,J) .gt. 1.e-5_r8) then
             EMKIKJ = EBOX(KJ) * B_GYA(KJ,J)
             EDY(J) = EDY(J) + EMKIKJ
          end if
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine E_GRID_Y
  !// ----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine TRUNG8(SPHI, SPLO, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                    JDGRD, IPARW, JPARW, IM, JM, LM, LV)
    !//---------------------------------------------------------------------
    !// Double precision!
    !// Average gridpoint data to go with truncation of spectral fields.
    !// Or, in other words, map high resolution SPHI onto low resolution
    !// SPLO.
    !//
    !// Will now also work if IPARW==IM and JPARW=JM
    !//
    !// Ole Amund Sovde, July 2015
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPARW, JPARW, IM, JM, LM, LV
    integer, intent(in) :: IDGRD, JDGRD
    integer, intent(in) :: IMAP(IDGRD,IM), JMAP(JDGRD,JM)
    real(r8), intent(in)  :: SPHI(IPARW,JPARW,LM,LV)
    real(r8), intent(in)  :: ZDEGI(IDGRD,IM),ZDEGJ(JDGRD,IM)
    !// Output
    real(r8), intent(out)  :: SPLO(IM,JM,LM,LV)
    !// Locals
    integer :: I,J,K,L,ix,jx
    !//---------------------------------------------------------------------

    !// If IPARW == IM and JPARW==JM, the resolution is the same!
    if (IPARW .eq. IM .and. JPARW .eq. JM) then
      do K = 1, LV
        do L = 1, LM
          do J = 1, JM
            do I = 1, IM
              SPLO(I,J,L,K) = SPHI(I,J,L,K)
            end do
          end do
        end do
      end do
    else
      !// Zero output array
      SPLO(:,:,:,:) = 0._r8
      !// Map SPHI onto SPLO
      do K = 1, LV
        do L = 1, LM
          do J = 1, JM
            do I = 1, IM
              do jx = 1, jdgrd
                do ix = 1, idgrd
                  SPLO(I,J,L,K) = SPLO(I,J,L,K) &
                                + SPHI(imap(ix,i), jmap(jx,j), L, K) &
                                  * ZDEGI(ix,i) * ZDEGJ(jx,j)
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    !//---------------------------------------------------------------------
  end subroutine TRUNG8
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine TRUNG4(SPHI, SPLO, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                    JDGRD, IPARW, JPARW, IM, JM, LM, LV)
    !//---------------------------------------------------------------------
    !// Single precision!
    !// Average gridpoint data to go with truncation of spectral fields.
    !// Or, in other words, map high resolution SPHI onto low resolution
    !// SPLO.
    !//
    !// Will now also work if IPARW==IM and JPARW=JM
    !//
    !// Ole Amund Sovde, July 2015
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, r4
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: IPARW, JPARW, IM, JM, LM, LV
    integer, intent(in) :: IDGRD, JDGRD
    integer, intent(in) :: IMAP(IDGRD,IM), JMAP(JDGRD,JM)
    real(r4), intent(in) :: SPHI(IPARW,JPARW,LM,LV)
    real(r8), intent(in) :: ZDEGI(IDGRD,IM),ZDEGJ(JDGRD,IM)
    !// Output
    real(r4), intent(out) :: SPLO(IM,JM,LM,LV)
    !// Locals
    integer :: I,J,K,L,ix,jx
    !//---------------------------------------------------------------------

    !// If IPARW == IM and JPARW==JM, the resolution is the same!
    if (IPARW .eq. IM .and. JPARW .eq. JM) then
      do K = 1, LV
        do L = 1, LM
          do J = 1, JM
            do I = 1, IM
              SPLO(I,J,L,K) = SPHI(I,J,L,K)
            end do
          end do
        end do
      end do
    else
      !// Zero output array
      SPLO(:,:,:,:) = 0._r4
      !// Map SPHI onto SPLO
      do K = 1, LV
        do L = 1, LM
          do J = 1, JM
            do I = 1, IM
              do jx = 1, jdgrd
                do ix = 1, idgrd
                  SPLO(I,J,L,K) = SPLO(I,J,L,K) &
                                + SPHI(imap(ix,i), jmap(jx,j), L, K) &
                                  * ZDEGI(ix,i) * ZDEGJ(jx,j)
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    !//---------------------------------------------------------------------
  end subroutine TRUNG4
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine Regrid_Column_Weights( PEdge1, PEdge2, Weights )
    !//---------------------------------------------------------------------
    !// Based on UCI qcode 72b.
    !//
    !// Ole Amund Sovde, April 2016
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    ! Input/Output
    real(r8), dimension(:), intent(in) :: PEDGE1, PEDGE2
    real(r8), dimension(:,:), intent(out):: Weights

    !// Local variables
    integer :: LM1, LM2, L, K
    real(r8) :: rsum
    logical, parameter :: CHECK=.true., DEBUG=.false.
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'Regrid_Column_Weights'
    !// --------------------------------------------------------------------

    !// Initialise
    LM1      = size( PEdge1 ) - 1
    LM2      = size( PEdge2 ) - 1
    Weights(:,:) = 0._r8


    !// Determine fraction of each INPUT box
    !// which contributes to each OUTPUT box

    !// Loop over INPUT layers
    do L = 1, LM1

       !// If top-of-box INPUT pressure is located below model
       !// surface, skip the data. I.e. do not put below-surface
       !// data into surface layer of model.
       if ( ( PEdge1(L+1) .ge. Pedge2(1) ) ) cycle

       !// Loop over OUTPUT layers
       do K = 1, LM2

          !==============================================================
          ! No contribution if:
          ! -------------------
          ! Bottom of OUTPUT layer above Top    of INPUT layer  OR
          ! Top    of OUTPUT layer below Bottom of INPUT layer
          !==============================================================
          if ( ( PEdge2(K)   .lt. PEdge1(L+1)) .or. &
               ( PEdge2(K+1) .gt. PEdge1(L)  ) ) cycle

          !==============================================================
          ! Contribution if: 
          ! ----------------
          ! Entire INPUT layer in OUTPUT layer
          !==============================================================
          if ( (PEdge2(K)   .ge. PEdge1(L))   .and. &
               (PEdge2(K+1) .le. PEdge1(L+1)) ) then 

             Weights(L,K) = 1._r8

             !// Go to next K iteration
             cycle
          end if

          !==============================================================
          ! Contribution if: 
          ! ----------------
          ! Top of OUTPUT layer in INPUT layer
          !==============================================================
          if ( (PEdge2(K+1) .le. PEdge1(L)) .and.  &
               (PEdge2(K)   .ge. PEdge1(L)) ) then 

             Weights(L,K) = ( PEdge1(L) - PEdge2(K+1) ) /  &
                            ( PEdge1(L) - PEdge1(L+1) ) 

             !// Go to next K iteration
             cycle
          end if

          !==============================================================
          ! Contribution if: 
          ! ----------------
          ! Entire OUTPUT layer in INPUT layer
          !==============================================================
          if ( (PEdge2(K)   .le. PEdge1(L))   .and.  &
               (PEdge2(K+1) .ge. PEdge1(L+1)) ) then

             Weights(L,K) = ( PEdge2(K) - PEdge2(K+1) ) /  &
                            ( PEdge1(L) - PEdge1(L+1) )

             !// Go to next K iteration
             cycle
          end if

          !==============================================================
          ! Contribution if: 
          ! ----------------
          ! Bottom of OUTPUT layer in INPUT layer
          !==============================================================
          if ( (PEdge2(K)   .ge. PEdge1(L+1)) .and. &
               (PEdge2(K+1) .le. PEdge1(L+1)) ) then

             Weights(L,K) = ( PEdge2(K) - PEdge1(L+1) ) /   &
                            ( PEdge1(L) - PEdge1(L+1) )

             !// Go to next K iteration
             cycle
          end if

       end do !// do K = 1, LM2


       !// Consistency Check:
       !// If SUM( WEIGHTS(L,:) ) does not = 1, there may be a problem.
       !//
       !// This is actually only true when Pedge1(1) .le. Pedge2(1).
       !// Otherwise, a only a partial of input is put in output,
       !// the rest going "below surface".
       rsum = sum( Weights(L,:) )
       if ( rsum .gt. 0._r8 .AND. Check ) then 
          if ( Pedge1(1) .le. Pedge2(1) .and. &
               abs( 1._r8 - rsum ) .ge. 1.e-4_r8 ) then
             write(6,'(a,i4,1x,f13.7,es20.16)') f90file//':'//subr// &
                  ': Weights does not add to 1:',L,rsum, &
                  Abs( 1._r8 - rsum )
             do K = 1, LM2
                print*,K,Pedge2(K),weights(L,K)
             end do
             print*,Pedge1
             stop
          end if
       end if

    end do !// do L = 1, LM1

    !//---------------------------------------------------------------------
  end subroutine Regrid_Column_Weights
  !//-----------------------------------------------------------------------




  !//-----------------------------------------------------------------------
end module regridding
!//=========================================================================

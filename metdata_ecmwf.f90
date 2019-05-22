!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, November 2015
!//=========================================================================
!// Routine for reading meteorological data from ECMWF, from nedCDF4 files.
!//=========================================================================
module metdata_ecmwf
  !//-----------------------------------------------------------------------
  !// MODULE: metdata_ecmwf
  !// DESCRIPTION: Routine for reading meteorological data from ECMWF,
  !//              stored on nedCDF4 format.
  !//
  !// Contains
  !//   subroutine update_metdata
  !//   subroutine fluxfilter2
  !//   subroutine data2mpblocks
  !//   subroutine gotData
  !//   subroutine skipData
  !//-----------------------------------------------------------------------
  use cmn_size, only: LPAR
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  integer :: LMAP(LPAR+1)
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'metdata_ecmwf.f90'
  !//-----------------------------------------------------------------------
  save
  private
  public update_metdata
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine update_metdata(ILOOP, DTMET, NMET)
    !//---------------------------------------------------------------------
    !// Read input data. This version reads netCDF4 files, where all data
    !// for one time step (NMET) is located in a single file.
    !//
    !// This routine assumes meteorology is to be updated each NMET.
    !//
    !// Ole Amund Sovde, April 2015
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, r4
    use cmn_size, only: IPARW, JPARW, LPARW, LWEPARW, &
         IPAR, JPAR, IDGRD, JDGRD, &
         LPAR, LWEPAR, LWDPAR, &
         LOSLOCHEM, LOSLOCTROP, NRMETD, LDUST
    use cmn_ctm, only: JYEAR, JMON, JDAY, JDATE, GMTAU, &
         JDATE_NEXT, JMON_NEXT, JYEAR_NEXT, &
         LMMAP, XLMMAP, TMON, TMET, &
         ALP, ALPV, TRIG, IFAX, LDEG, ZDEGI, ZDEGJ, IMAP, JMAP, &
         XYZA, XYZB, IMEPZ, DD, SS, ETAAW, ETABW, ETAA, ETAB, &
         WGLYE, WGLYG, YDGRD, AREAXYW, AREAXY, PLAND, LFIXMET
    use cmn_met, only: P, U, V, T, Q, CWETN, CWETE, CWETD, CENTU, CENTD, &
         PRECCNV, PRECLS, ZOFLE, &
         CLDFR, CLDLWC, CLDIWC, &
         SLH, SHF, SMF, SFT, SFQ, SFU, SFV, BLH, BLH_CUR, BLH_NEXT, LBLH, &
         SA, PVU, UMS, VMS, CI, SD, SMLT, ES, SNFL, PhotActRad, MSLP, USTR, &
         SWVL1, SWVL3, STL1, &
         PMEAN, PMEANW, MYEAR, &
         metTYPE, MPATH1, MPATH2, MFILE3
    use cmn_parameters, only: M_AIR, R_UNIV, R_AIR, R_H2O, A0, CPI
    use cloudjx, only: SWSTORE, CLDSTORE, TYPSTORE, RAN4, IRAN0, &
                       LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
    use regridding, only: TRUNG4, TRUNG8
    use utilities, only: ctmExitC, CFRMIN,CIWMIN
    use ncutils, only: get_netcdf_var_2d, get_netcdf_var_3d
    use dust_oslo, only: dust_set_ssrd, dust_set_strd, dust_set_SWVL1
    use physics_oslo, only: get_pvu, ijlw2lij
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: ILOOP, & ! <=0, init, > 0 ordinary run
                           NMET     ! the meteorological timestep
    real(r8), intent(in) :: DTMET     ! meteorological time step [s]
    !//---------------------------------------------------------------------
    logical, parameter :: VERBOSE = .false.
    real(r8), parameter :: EPS = 0.01_r8 ! minimum cloud fraction
    !// Limits for convective fluxes
    real(r8), parameter :: mincwetelim =  3.8e-4_r8
    real(r8), parameter :: maxcwetdlim = -2.0e-5_r8
    real(r8), parameter :: mincdetulim =  1.e-4_r8
    real(r8), parameter :: mincdetdlim =  1.e-6_r8

    !// Local parameters
    real(r8) :: &
         VEDGE, BAND, DETA, DETB, DELP, &
         ZCOS, ZDT, SFTD, LV, ESAT, SDEN, &
         DELZ, QMIN, ZTOP, ZBOT, ZMID, &
         DMASS, PSRF, &
         POFLE(LPAR+1)

    !// Indices
    integer :: I,J,II,JJ,L,LL

    !// To be read from file
    logical :: fex

    !// Filename for metdata
    character(len=160) :: FNAME, FNAME_NEXT
    character(len=2) :: CMON, CDATE, CUTC
    integer :: NMET_NEXT

    !//---------------------------------------------------------------------

    !// Allocatable arrays - double precision
    real(r8), dimension(:), allocatable :: &
         UTMP(:)
    real(r8), dimension(:,:), allocatable :: &
         VTMP, PW, EWSS, NSSS, W2D, R8XY
    !real(r8), dimension(:,:), allocatable :: &
    !     LSPREC, CNVPREC
    real(r8), dimension(:,:,:), allocatable :: &
         TW, QW, CLDFRW, CLDIWCW, CLDLWCW, CDETU, CDETD, &
         ZOFLEW, W3Da, W3Db, R8XYZ
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'update_metdata'
    !//---------------------------------------------------------------------
         
    !// Allocate 3D arrays - native resolution
    allocate( VTMP(IPARW,JPARW), UTMP(IPARW), &
              W3Da(IPARW,JPARW,LPARW), W3Db(IPARW,JPARW,LPARW), &
              CLDFRW(IPARW,JPARW,LWEPAR), &
              CLDIWCW(IPARW,JPARW,LWEPAR), CLDLWCW(IPARW,JPARW,LWEPAR) )

    !// Allocate 3D arrays - native horizontal resolution
    allocate( ZOFLEW(LPAR+1,IPARW,JPARW), QW(IPARW,JPARW,LPAR), &
              TW(IPARW,JPARW,LPAR) )

    !// Allocate 3D arrays - window resolution (IPAR/JPAR)
    allocate( CDETU(IPAR,JPAR,LWEPAR), CDETD(IPAR,JPAR,LWDPAR), &
              R8XYZ(IPAR,JPAR,LPARW) )

    !// Allocate 2D arrays - native resolution
    allocate( PW(IPARW,JPARW), W2D(IPARW,JPARW) )

    !// Allocate 2D arrays - window resolution (IPAR/JPAR)
    allocate( EWSS(IPAR,JPAR), NSSS(IPAR,JPAR), R8XY(IPAR,JPAR) )
    !allocate( LSPREC(IPAR,JPAR),CNVPREC(IPAR,JPAR) )
    !//---------------------------------------------------------------------

    !// Initialize
    !// Time step for meteorological data is DTMET
    ZDT    = 1._r8 / DTMET

    !// Define minimum humidity QMIN in kg/kg
    QMIN = 3.e-6_r8 * 18._r8 / M_AIR

!//UPDATE
    !// locate the position of random number sequence based on year/day/hour
    IRAN0 = 1 + 3*(JYEAR - 1900) + 7*JDAY + 11*nint(GMTAU)


    !// File name for this DAY/NMET
    write(CUTC(1:2),'(i2.2)') (NMET - 1) * 24/NRMETD  ! Time UTC
    write(CDATE(1:2),'(i2.2)') JDATE                  ! Date
    write(CMON(1:2),'(i2.2)') JMON                    ! Month
    write(MPATH2(1:4),'(i4.4)') MYEAR                 ! Year
    if (trim(metTYPE) .eq. 'ECMWF_oIFSnc4') then
       !// ECMWF Open IFS netcdf4 file
       !// Date yYYYYmMMdDDhHH
       !// Filename: ECopenIFSc38r1_yYYYYmMMdDDhHH_T159N80L60.nc
       write(MFILE3(16:29),'(a1,i4.4,a1,i2.2,a1,a2,a1,a2)') &
            'y',MYEAR,'m',JMON,'d',CDATE,'h',CUTC

       FNAME  = trim(MPATH1)//trim(MPATH2)//'/'//CMON//'/'// &
            trim(MFILE3)//'.nc'

       !// File for next time step
       if (NMET .eq. 8) then
          NMET_NEXT = 1
          write(CUTC(1:2),'(i2.2)') (NMET_NEXT - 1) * 24/NRMETD    ! Time UTC
          write(CDATE(1:2),'(i2.2)') JDATE_NEXT                    ! Date
          write(CMON(1:2),'(i2.2)') JMON_NEXT                      ! Month
          if (.not.LFIXMET) then
             i = JYEAR_NEXT
          else
             i = MYEAR
          end if
          write(MPATH2(1:4),'(i4.4)') i ! Year
       else
          NMET_NEXT = NMET + 1
          if (.not.LFIXMET) then
             i = JYEAR
          else
             i = MYEAR
          end if
          write(CUTC(1:2),'(i2.2)') (NMET_NEXT - 1) * 24/NRMETD    ! Time UTC
       end if
       write(MFILE3(16:29),'(a1,i4.4,a1,i2.2,a1,a2,a1,a2)') &
            'y',i,'m',JMON_NEXT,'d',CDATE,'h',CUTC
       FNAME_NEXT  = trim(MPATH1)//trim(MPATH2)//'/'//CMON//'/'// &
            trim(MFILE3)//'.nc'

    else
       write(6,'(a)') f90file//':'//subr// &
            ': Not set up for metTYPE: '//trim(metTYPE)
       if (trim(metTYPE) .eq. 'ECMWF_oIFS') then
          write(6,'(a)') '* If you do not want netcdf4 files, you '// &
               'should use metdata_ecmwf_uioformat.f90 instead'
       end if
       stop
    end if

    !// Check if files exist
    inquire(FILE=trim(FNAME), exist=fex)
    if (.not. fex) then
       write(6,'(a)') f90file//':'//subr// &
            ': No such file: '//trim(FNAME)
       stop
    end if



    !//---------------------------------------------------------------------
    !// Initial step - setup P, T and Q
    !//---------------------------------------------------------------------
    if (ILOOP .le. 0) then

       !// Set up level weightings if vertical resolution degraded
       do LL = LPARW+1, 1, -1
          L   = LMMAP(LL)
          LMAP(L) = LL
       end do

       write(6,'(a)') 'Initializing meteorological data'
       write(6,'(2x,a)') trim(FNAME)


       !// Pressure field --------------------------------------------------
       !// -----------------------------------------------------------------
       call get_netcdf_var_2d(FNAME, 'pres_sfc',PW, IPARW, JPARW)
       !// Put PW into P, degrading or not
       call TRUNG8(PW, P, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)

       if (verbose) call gotData('SFC','Surface pressure (P)')


       !// Temperature -----------------------------------------------------
       !// -----------------------------------------------------------------
       call get_netcdf_var_3d(FNAME, 'temperature', W3Da, &
            IPARW, JPARW, LPARW)
       !// Collapse layers
       TW(:,:,:) = 0._r8
       do L = 1, LPAR
         do LL = LMAP(L), LMAP(L+1) - 1
           do J = 1, JPARW
             do I = 1, IPARW
               TW(I,J,L) = TW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
             end do
           end do
         end do
       end do
       !// Put TW into T, degrading or not
       call TRUNG8(TW, T, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)

       !// Test to check temperature field
       if (minval(W3Da(:,:,1)) .le. 0._r8) then
          write(6,'(a)')  f90file//':'//subr// &
               ': Surface temperature is <= 0?'
          write(6,'(a,f9.3)') '  MIN sfc T: ',minval(W3Da(:,:,1))
          stop
       end if

       if (verbose) call gotData('3di','Temperature (T)')


       !// Water vapour ----------------------------------------------------
       !// -----------------------------------------------------------------
       call get_netcdf_var_3d(FNAME, 'qhum', W3Da, &
            IPARW, JPARW, LPARW)
       !// Collapse layers
       QW(:,:,:) = 0._r8
       do L = 1,LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do J = 1,JPARW
             do I = 1,IPARW
               QW(I,J,L) = QW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
             end do
           end do
         end do
       end do
       !// There may be some entries of small negative numbers
       if (minval(QW) .lt. 0._r8) then
          print*,'update_metdata: min QW:',minval(QW),', setting to',QMIN
          do L = 1, LPAR
            do J = 1, JPARW
              do I = 1, IPARW
                if (QW(I,J,L) .lt. 0._r8) QW(I,J,L) = QMIN
              end do
            end do
          end do
       end if

       !// Put QW into Q, degrading or not
       call TRUNG8(QW, Q, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR,1)

       if (verbose) call gotData('3di','Specific humidity (Q)')

       !// Polar cap filtering
       call EPZ_TQ(T, Q, XYZA, XYZB, P, IMEPZ,IPAR,JPAR,LPAR,IPAR,JPAR,LPAR)
       call EPZ_P(P, IMEPZ, IPAR,JPAR, IPAR, JPAR)

       !// Initialization is done - these steps are for all time steps
       return
    end if !// if (ILOOP .le. 0) then


    !//---------------------------------------------------------------------
    !// All time steps - update all meteorological fields
    !//---------------------------------------------------------------------
    if (verbose) then
       write(6,'(a,i5)') f90file//':'//subr// &
            ': Reading new metdata JDAY: '//TMET//', NMET:',NMET
       write(6,'(2x,a)') trim(FNAME)
    end if

    !// Clear arrays
    U(:,:,:) = 0._r8
    V(:,:,:) = 0._r8
    T(:,:,:) = 0._r8
    Q(:,:,:) = 0._r8
    CWETN(:,:,:) = 0._r8
    CWETE(:,:,:) = 0._r8
    CWETD(:,:,:) = 0._r8
    CENTU(:,:,:) = 0._r8
    CENTD(:,:,:) = 0._r8
    CDETU(:,:,:) = 0._r8
    CDETD(:,:,:) = 0._r8
    PRECCNV(:,:,:) = 0._r8
    PRECLS(:,:,:) = 0._r8
    CLDFR(:,:,:) = 0._r8
    CLDLWC(:,:,:) = 0._r8
    CLDIWC(:,:,:) = 0._r8
    !// No need to clear 2D arrays; they are fully updated below
    SLH(:,:) = 0._r8
    SHF(:,:) = 0._r8
    SMF(:,:) = 0._r8
    SFT(:,:) = 0._r8
    SFQ(:,:) = 0._r8
    SFU(:,:) = 0._r8
    SFV(:,:) = 0._r8
    BLH(:,:) = 0._r8
    MSLP(:,:) = 0._r8
    SA(:,:) = 0._r8



    !//---------------------------------------------------------------------
    !// 3D GRIDDED DATA - INSTANTANEOUS
    !//---------------------------------------------------------------------

    !// Pressure field -----------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_2d(FNAME, 'pres_sfc',PW, IPARW, JPARW)
    !// Put PW into P, degrading or not
    call TRUNG8(PW, P, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)

    if (verbose) call gotData('SFC','Surface pressure (P)')


    !// Temperature --------------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_3d(FNAME, 'temperature', W3Da, &
         IPARW, JPARW, LPARW)
    !// Collapse layers
    TW(:,:,:) = 0._r8
    do L = 1, LPAR
      do LL = LMAP(L), LMAP(L+1) - 1
        do J = 1, JPARW
          do I = 1, IPARW
            TW(I,J,L) = TW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    !// Put TW into T, degrading or not
    call TRUNG8(TW, T, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)

    if (verbose) call gotData('3di','Temperature (T)')


    !// Zonal wind (U) -----------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_3d(FNAME, 'U', W3Da, IPARW, JPARW, LPARW)

    !// scale from U=u*cos() to u and to flux [kg/s]
    do L = 1, LPARW
       DETA = ETAAW(L) - ETAAW(L+1)
       DETB = ETABW(L) - ETABW(L+1)
       do J = 1, JPARW
          BAND = A0 * (WGLYE(J+1) - WGLYE(J))
          ZCOS = 1._r8 / COS(WGLYG(J))
          do I = 1, IPARW
             UTMP(I) = W3Da(I,J,L) * ZCOS * BAND
             !// Save center values for m/s
             W3Db(I,J,L) = W3Da(I,J,L) * ZCOS !// Save U(m/s)
          end do
          !// Get edge values (unit conversion needs pressure on edge)
          do I = 2, IPARW
             DELP = DETA + DETB * (PW(I-1,J) + PW(I,J)) * 0.5_r8
             W3Da(I,J,L) = (UTMP(I-1) + UTMP(I)) * 0.5_r8 * DELP
          end do
          DELP = DETA + DETB * (PW(1,J) + PW(IPARW,J)) * 0.5_r8
          W3Da(1,J,L) = (UTMP(IPARW) + UTMP(1)) * 0.5_r8 * DELP
       end do
    end do

    if (ldeg) then
       do L = 1, LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do I = 1, IPAR
             II = IMAP(1,I)
             do J = 1, JPAR
               do JJ = 1, JDGRD
                 U(I,J,L) = U(I,J,L) + W3Da(II,JMAP(JJ,J),LL)
               end do
             end do
           end do
         end do
       end do
    else
       do L = 1, LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do J = 1, JPAR
             do I = 1, IPAR
               U(I,J,L) = U(I,J,L) + W3Da(I,J,LL)
             end do
           end do
         end do
       end do
    end if

    if (verbose) call gotData('3di','Zonal wind (U)')

    !// Save UMS -----------------------------------------------------------
    !//---------------------------------------------------------------------
    !// Degrade vertically and horizontally, put into (LPAR,IPAR,JPAR)
    call ijlw2lij(W3Db, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
         JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, UMS)


    !// Meridional wind (V) ------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_3d(FNAME, 'V', W3Db, IPARW, JPARW, LPARW)

    !// scale from V=v*cos() to v (remember V is on edge already)
    do L = 1, LPARW
       DETA = ETAAW(L) - ETAAW(L+1)
       DETB = ETABW(L) - ETABW(L+1)
       VTMP(:,:) = 0._r8 !// Needed for getting VMSW
       !// W3Db does not start at SP, but at the J=2. It means
       !// W3Db is flux out of box J, not into J.
       !// Given all edge points, the number of boxes should be
       !// JPARW+1, covering 1:JPARW+1. But both the poles should
       !// have zero wind, so we only need JPARW-1 boxes.
       !//
       !// But in addition, the values at JPARW/2 and JPARW/2+1
       !// are duplicated, both representing wind across Equator.
       !// So W3Db has the size JPARW instead of JPARW-1.
       !//
       !// A bit confusing, the CTM uses V as flux into J-box instead of
       !// out of it. This will be taken care of at the end.
       !// SH flux out of J-box:
       do J = 1, JPARW/2
          ZCOS = 1._r8 / COS(WGLYE(J+1))
          do I = 1, IPARW
             VTMP(I,J) = W3Db(I,J,L) * ZCOS
          end do
       end do
       !// Now VTMP covers the V up to Equator, starting at the
       !// models index J=2, but stored in VTMP at J=1.
       !// NH flux out of J-box (JPARW/2 and JPARW/2+1 are duplicated):
       do J = JPARW/2 + 2, JPARW
          ZCOS = 1._r8 / COS(WGLYE(J))
          do I = 1, IPARW
             VTMP(I,J-1) = W3Db(I,J,L) * ZCOS
          end do
       end do
       !// Thus, at J=JPARW, VTMP (with its indicing) represents NP
       !// and is therefore zero.

       !//change unit, m/s ==> Kg/s (100./G0 is done in pdyn0.f)
       !// Save center values for m/s J=1 (assume V=0 at SP)
       W3Da(:,1,L) = VTMP(:,1) * 0.5_r8
       !// Map the VTMP back to correct model edge grid (from J to J+1),
       !// i.e. converting from flux out of J-1 to flux into J.
       do J = 2, JPARW
          VEDGE  = 2._r8 * CPI * A0 * COS(WGLYE(J)) / real(IPARW, r8)
          do I = 1, IPARW
             DELP = DETA + DETB * (PW(I,J-1) + PW(I,J)) * 0.5_r8
             !// Save center values for m/s J>1
             W3Da(I,J,L) = (VTMP(I,J) + VTMP(I,J-1)) * 0.5_r8
             W3Db(I,J,L) = VTMP(I,J-1) * VEDGE * DELP
          end do
       end do
    end do

    if (ldeg) then
       do L = 1, LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do J = 1, JPAR
             JJ = JMAP(1,J)
             do I = 1, IPAR
               do II = 1, IDGRD
                 V(I,J,L) = V(I,J,L) + W3Db(IMAP(II,I),JJ,LL)
               end do
             end do
           end do
         end do
       end do
    else
       do L = 1, LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do J = 1, JPAR
             do I = 1, IPAR
               V(I,J,L) = V(I,J,L) + W3Db(I,J,LL)
             end do
           end do
         end do
       end do
    end if

    if (verbose) call gotData('3di','Meridional wind (V)')


    !// Save VMS -----------------------------------------------------------
    !//---------------------------------------------------------------------
    !// Degrade vertically and horizontally, put into (LPAR,IPAR,JPAR)
    call ijlw2lij(W3Da, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
         JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, VMS)




    !//---------------------------------------------------------------------
    !// 3-d GRID POINT DATA
    !//---------------------------------------------------------------------


    !// Water vapour -------------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_3d(FNAME, 'qhum', W3Da, IPARW, JPARW, LPARW)
    !// Collapse layers
    QW(:,:,:) = 0._r8
    do L = 1, LPAR
      do LL = LMAP(L), LMAP(L+1)-1
        do J = 1, JPARW
          do I = 1, IPARW
            QW(I,J,L) = QW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    if (minval(QW) .lt. 0._r8) then
       print*,'update_metdata: min QW:',minval(QW),', setting to',QMIN
       do L = 1, LPAR
         do J = 1, JPARW
           do I = 1, IPARW
             if (QW(I,J,L) .lt. 0._r8) QW(I,J,L) = QMIN
           end do
         end do
       end do
    end if

    !// Put QW into Q, degrading or not
    call TRUNG8(QW, Q, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)

    if (verbose) call gotData('3di','Specific humidity (Q)')


    !// Altitudes - calculated ---------------------------------------------
    !// --------------------------------------------------------------------
    !// R = 287.*(1-Q) + 461.5*Q -- assume 0.5% bl w.v.==> R = 288.
    !// delta-z (m) = dln(P) * R * T / g   where R/g = 288/9.81 = 29.36
    !// CTM3: We use Q directly instead of assuming 0.5%;
    !// Using Rd=287 and Tv=T*(1 + 0.6*q); 287/9.80665=29.26586
    do J = 1, JPAR
      do I = 1, IPAR
        PSRF  = P(I,J)
        !// Surface pressure; must be set for DELZ to be calculated
        POFLE(1)  = PSRF
        !// Height of box bottom, above sea level (acts as topography)
        ZOFLE(1,I,J)  = 16.e3_r8 * log10(1013.25_r8 / PMEAN(I,J))
        if (ZOFLE(1,I,J) .ne. ZOFLE(1,I,J)) then
           write(6,'(a,i3,2i5,2es12.2)')  f90file//':'//subr// &
                ': ZOFLE a',1,i,j,PMEAN(I,J),PMEANW(I,J)
           write(6,*) pmean(:,j)
           stop
        end if
        do L = 2, LPAR + 1
          !// Pressure of box bottom
          POFLE(L)  = ETAA(L) + ETAB(L) * PSRF
          !DELZ  = -29.36_r8 * T(I,J,L-1) * log(POFLE(L)/POFLE(L-1))
          !// Thickness of layer (L-1) (remember ZOFLE starts at topography)
          DELZ = -29.26586_r8 * T(I,J,L-1) * (1._r8 + 0.6_r8 * Q(I,J,L-1)) &
                             * log(POFLE(L)/POFLE(L-1))
          !// Add DELZ of (L-1) to get box bottom height of L
          ZOFLE(L,I,J) = ZOFLE(L-1,I,J) + DELZ
          if (ZOFLE(L,I,J) .ne. ZOFLE(L,I,J)) then
             write(6,'(a,i3,2i5,4es12.2)')  f90file//':'//subr// &
                  ': ZOFLE b',l,i,j,delz,q(i,j,l-1),&
                  T(I,J,L-1),psrf
             stop
          end if
        end do
      end do
    end do

    !// Also set up ZOFLEW for metdata grid
    do J = 1, JPARW
      do I = 1, IPARW
        PSRF  = PW(I,J)
        POFLE(1)  = PSRF
        ZOFLEW(1,I,J)  = 16e3_r8 * log10(1013.25_r8 / PMEANW(I,J))
        do L = 2, LPAR + 1
          POFLE(L)  = ETAA(L) + ETAB(L)*PSRF
          !DELZ  = -29.36_r8 * TW(I,J,L-1) * log(POFLE(L)/POFLE(L-1))
          !// Thickness of layer below
          DELZ = -29.26586_r8 * TW(I,J,L-1)*(1._r8 + 0.6_r8*QW(I,J,L-1)) &
                             * log(POFLE(L)/POFLE(L-1))
          ZOFLEW(L,I,J) = ZOFLEW(L-1,I,J) + DELZ
        end do
      end do
    end do


    !// Potantial vorticity ------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_3d(FNAME, 'PV', W3Da, IPARW, JPARW, LPARW)
    !// Convert to PVU
    W3Da(:,:,:) = W3Da(:,:,:) * 1.e6_r8
    if (maxval(W3Da) .eq. 0._r8) then
       !// PV was not on available.
       !// Generate PV on model resolution and convert to PVU
       call get_pvu()
    else
       !// Transform into PVU-array
       call ijlw2lij(W3Da, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
            JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, PVU)
    end if


    !// Cloud Liquid Water Content [Kg/Kg] - midpoint value ----------------
    !// --------------------------------------------------------------------
    !// Cloud routine needs CLDLWCW in any case, so we calculate
    !// that first. Could repeat with R8XYZ, but rather do
    !// degradation of CLDLWCW.
    call get_netcdf_var_3d(FNAME, 'lwc', W3Da, IPARW, JPARW, LPARW)

    CLDLWCW(:,:,:)  = 0._r8
    do L = 1, LWEPAR
      do LL = LMAP(L), LMAP(L+1)-1
        do J = 1, JPARW
          do I = 1, IPARW
            CLDLWCW(I,J,L) = CLDLWCW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    !// Put CLDLWCW into LCDLWC, degrading or not
    call TRUNG8(CLDLWCW, CLDLWC, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)


    !// Cloud Ice Water Content [Kg/Kg] - midpoint value -------------------
    !// --------------------------------------------------------------------
    !// Cloud routine needs CLDIWCW in any case, so we calculate
    !// that first.
    call get_netcdf_var_3d(FNAME, 'iwc', W3Da, IPARW, JPARW, LPARW)

    CLDIWCW(:,:,:)  = 0._r8
    do L = 1, LWEPAR
      do LL = LMAP(L), LMAP(L+1)-1
        do J = 1, JPARW
          do I = 1, IPARW
            CLDIWCW(I,J,L) = CLDIWCW(I,J,L) + W3Da(I,J,LL) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    !// Put CLDIWCW into LCDIWC, degrading or not
    call TRUNG8(CLDIWCW, CLDIWC, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)


    !// Cloud Fraction [0,1] - midpoint value ------------------------------
    !// --------------------------------------------------------------------
    !// Cloud routine needs CLDFRW in any case, so we calculate
    !// that first.
    call get_netcdf_var_3d(FNAME, 'cfr', W3Da, IPARW, JPARW, LPARW)

    CLDFRW(:,:,:) = 0._r8
    do L = 1, LWEPAR
      do LL = LMAP(L), LMAP(L+1)-1
        do J = 1, JPARW
          do I = 1, IPARW
            CLDFRW(I,J,L) = CLDFRW(I,J,L) &
                    + max(min(W3Da(I,J,LL),1._r8),0._r8) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    !// Put CLDFRW into LCDFR, degrading or not
    call TRUNG8(CLDFRW, CLDFR, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)


    !//---------------------------------------------------------------------
    !// 3D GRIDDED DATA - ACCUMULATED
    !//---------------------------------------------------------------------

    !// Mass Flux updrafts [accumulated kg/(m^2*s)] ------------------------
    !// --------------------------------------------------------------------
    !// - Mass flux is through BOTTOM EDGE of grid box (i.e. model
    !//   half layers).
    !// - EC-data are stored from half layer 2, since the flux into
    !//   layer 1 is always zero (no flux in through surface).
    !// - Scale with area and dt  --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'cflxu', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    CWETE(:,:,:) = 0._r8 !// No flux up from surface
    !// Start retrieval at layer 2
    do L = 2, LWEPAR
       !// LMAP is full-level, but for since mass fluxes are not stored
       !// for surface, and starts at layer 2, we need to subtract 1.
       !// What comes into L=2, comes from LMAP(2)-1.
       !// For L40, L=2 and LMAP(2)-1 = 1, while
       !// for L37, L=2 and LMAP(2)-1 = 3
       LL = LMAP(L) - 1
       do J = 1, JPAR
          do I = 1, IPAR
             !// Filter values. File data has minimum values less than zero.
             !// Hence we treat values less that that as as zero
             if (R8XYZ(I,J,LL) .gt. mincwetelim) then
                CWETE(I,J,L) = R8XYZ(I,J,LL) * AREAXY(I,J) * ZDT
             else
                CWETE(I,J,L) = 0._r8
             end if
          end do
       end do
    end do


    !// Mass Flux downdrafts [accumulated kg/(m^2*s)] ----------------------
    !// --------------------------------------------------------------------
    !// - Mass flux is through BOTTOM EDGE of grid box (i.e. model
    !//   half layers).
    !//   This means that the flux for layer 1 is zero, and that
    !//   downward flux is NEGATIVE, going out at the bottom of the box.
    !// - EC-data are stored from half layer 2, since the flux into
    !//   layer 1 is always zero (no flux in through surface).
    !// - For a grid box in layer L, -FD(L) goes out at bottom
    !//   and -FD(L+1) comes in from above (FD is negative). The
    !//   balance with entrainment E and detrainment D is:
    !//   FD(L) - FD(L+1) = E - D
    !// - Scale with area and dt  --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'cflxd', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    CWETD(:,:,:) = 0._r8 !// No flux down to surface
    !// Start retrieval at layer 2
    do L = 2, LWDPAR
       !// See Case(212) for comment on LMAP
       LL = LMAP(L) - 1
       do J = 1, JPAR
          do I = 1, IPAR
             if (R8XYZ(I,J,LL) .lt. maxcwetdlim) then
                CWETD(I,J,L) = R8XYZ(I,J,LL) * AREAXY(I,J) * ZDT
             else
                CWETD(I,J,L) = 0._r8
             end if
          end do
       end do
    end do



    !// Updrafts detrainment rate [accumulated kg/(m3*s)] ------------------
    !// [i.e.  kg/(m2*s) per gridbox height]              ------------------
    !// --------------------------------------------------------------------
    !// - Entrainment has to be built from detrainment.
    !// - The rate is per height, so we have to multiply with dZ and area.
    !// - Detrainment/entrainment rates are given at GRID CENTER
    !//   (i.e. model full layers).
    !// - Detrainment must be summed up when collapsing layers!
    !// - Scale with box height, area and dt --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'cdetu', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    CDETU(:,:,:) = 0._r8
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LWEPAR
             do LL = LMAP(L), LMAP(L+1) - 1
                !// Detrainment must be summed up when collapsing layers
                if (R8XYZ(I,J,LL) .gt. mincdetulim) then
                   CDETU(I,J,L) = CDETU(I,J,L) &
                        + R8XYZ(I,J,LL) * XLMMAP(LL) &
                          * (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) &
                          * AREAXY(I,J) * ZDT
                end if
             end do
          end do
       end do
    end Do


    !// Downdrafts detrainment rate [acc. kg/(m3*s)] -----------------------
    !// [i.e. kg/(m2*s) per gridbox height]          -----------------------
    !// --------------------------------------------------------------------
    !// - Entrainment has to be built from detrainment.
    !// - The rate is per height, so we have to multiply with dZ and
    !//   area.
    !// - Detrainment/entrainment rates are given at GRID CENTER
    !//   (i.e. model full layers).
    !// - Detrainment must be summed up when collapsing layers!
    !// - Scale with box height and dt --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'cdetd', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    CDETD(:,:,:) = 0._r8
    do J = 1, JPAR
       do I = 1, IPAR
          do L = 1, LWDPAR
             do LL = LMAP(L),LMAP(L+1)-1
                if (R8XYZ(I,J,LL) .gt. mincdetdlim) then
                   CDETD(I,J,L) = CDETD(I,J,L) &
                        + R8XYZ(I,J,LL) * XLMMAP(LL) &
                          * (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) &
                          * AREAXY(I,J) * ZDT
                end if
             end do
          end do
       end do
    end do



    !// CONVECTIVE RAINFALL [kg/(m^2)] (accumulated kg/(m^2*s)) ------------
    !// --------------------------------------------------------------------
    !// EDGE VALUE
    !// Scale with area and dt or re-initialize --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'convrain', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    PRECCNV(:,:,:) = 0._r8
    do L = 1, LWEPAR
       LL = LMAP(L)
       do J = 1, JPAR
          do I = 1, IPAR
             PRECCNV(I,J,L) = &
                  max(0._r8, R8XYZ(I,J,LL) * AREAXY(I,J) * ZDT)
          end do
       end do
    end do



    !// LARGE SCALE RAINFALL [kg/(m^2)] (accumulated kg/(m^2*s)) -----------
    !// --------------------------------------------------------------------
    !// EDGE VALUE
    !// Scale with area and dt or re-initialize --> [kg/s]
    call get_netcdf_var_3d(FNAME, 'lsrain', W3Da, IPARW, JPARW, LPARW)
    W3Da(:,:,LWEPARW+1:LPARW) = 0._r8
    !// Possibly degrade without collapsing
    call TRUNG8(W3Da, R8XYZ, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)

    PRECLS(:,:,:) = 0._r8
    do L = 1, LWEPAR
       LL = LMAP(L)
       do J = 1, JPAR
          do I = 1, IPAR
             PRECLS(I,J,L) = &
                  max(0._r8, R8XYZ(I,J,LL) * AREAXY(I,J) * ZDT)
          end do
       end do
    end do




    !//---------------------------------------------------------------------
    !// Cloud cover routines - Only needed for J-values, i.e. chemistry.
!    if (LOSLOCHEM) &
    if (LOSLOCTROP) &
         !// New cloud treatment (qcode_60a)
         call CLOUD(CLDFRW,CLDIWCW,CLDLWCW,PW,TW,ETAA,ETAB,AREAXYW, &
                    ZOFLEW,ZDEGI,ZDEGJ,IMAP,JMAP, &
                    SWSTORE,CLDSTORE,TYPSTORE,RAN4,IRAN0, &
                    LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)
    !//---------------------------------------------------------------------





    !//---------------------------------------------------------------------
    !// 2-d SURFACE GRID POINT DATA
    !//---------------------------------------------------------------------

    !// SEA ICE COVER (CI) -------------------------------------------------
    !// --------------------------------------------------------------------
    call get_netcdf_var_2d(FNAME, 'CI',W2D, IPARW, JPARW)
    call TRUNG8(W2D, CI, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','Sea ice (CI)')


    !// Snow Evaporation (ES) (accumulated) --------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m/s water equivalent, accumulated (m w.eq.)
    !// Snow evaporation may actually also be negative, for some reason
    call get_netcdf_var_2d(FNAME, 'ES',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    R8XY(:,:) = R8XY(:,:) * ZDT
    call data2mpblocks(R8XY, ES)
    if (verbose) call gotData('2da','Snow Evaporation (ES)')


    !// Snow Melt (SMLT) (accumulated) -------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m/s water equivalent, accumulated (m w.eq.)
    call get_netcdf_var_2d(FNAME, 'SMLT',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    R8XY(:,:) = max(0._r8, R8XY(:,:) * ZDT)  !// Limit to positive values
    call data2mpblocks(R8XY, SMLT)
    if (verbose) call gotData('2da','Snow Melt (SMLT)')


    !// Photosynthetically active radiation @ sfc (PhotActRad) (accumulated)
    !// --------------------------------------------------------------------
    !// Unit: (W/m2)*s, accumulated. Divide by ZDT to get W/m2.
    call get_netcdf_var_2d(FNAME, 'PAR',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    PhotActRad(:,:) = max(0._r8, R8XY(:,:) * ZDT)  !// Limit to positive values
    if (verbose) call gotData('2da','Photosyn. act. rad. sfc (PhotActRad)')



    !// Snow depth (SD) ----------------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m water equivalent
    call get_netcdf_var_2d(FNAME, 'SD',W2D, IPARW, JPARW)
    call TRUNG8(W2D, SD, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','Snow depth (SD)')


    !// Large Scale Precipitation (stratiform) (accumulated) ---------------
    !// --------------------------------------------------------------------
    !call get_netcdf_var_2d(FNAME, 'LSPREC',W2D, IPARW, JPARW)
    !call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
    !     JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    !LSPREC(:,:) = R8XY(:,:) * ZDT
    !if (verbose) call gotData('2da','Large Scale Precip. (LSPREC)')


    !// Convective Precipitation (accumulated) -----------------------------
    !// --------------------------------------------------------------------
    !call get_netcdf_var_2d(FNAME, 'CONVPREC',W2D, IPARW, JPARW)
    !call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
    !     JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    !CNVPREC(:,:) = R8XY(:,:) * ZDT
    !if (verbose) call gotData('2da','Conv. Precip. (CONVPREC)')


    !// SnowFall SF (accumulated) ------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m/s, accumulated (m water equivalent)
    call get_netcdf_var_2d(FNAME, 'SF',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    R8XY(:,:) = R8XY(:,:) * ZDT
    call data2mpblocks(R8XY, SNFL)
    if (verbose) call gotData('2da','Snow fall (SNFL)')


    !// Surface Sensible Heat Flux SSHF (accumulated) ----------------------
    !// --------------------------------------------------------------------
    !// Unit: W/m2, accumulated (W/m2*s)
    call get_netcdf_var_2d(FNAME, 'SSHF',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    SHF(:,:) = -R8XY(:,:) * ZDT
    !// Fluxes of zero cause grief with deposition, PBL.
    !// SHF numbers from -400 to 100, set a small number to avoid
    !// overflow, not worry about small negative going to small positive
    do J = 1, JPAR
       do I = 1, IPAR
          if (abs(SHF(I,J)) .lt. 1.e-6_r8) SHF(I,J) = 1.e-6_r8
       end do
    end do
    if (verbose) call gotData('2da','Sfc. sens. heat flux (SHF)')


    !// Surface Latent Heat Flux SLHF (accumulated) ------------------------
    !// --------------------------------------------------------------------
    !// Unit: W/m2, accumulated (W/m2*s)
    call get_netcdf_var_2d(FNAME, 'SLHF',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    SLH(:,:) = -R8XY(:,:) * ZDT
    !// Fluxes of zero cause grief with deposition, PBL.
    !// SLH: set small value if zero.
    do J = 1, JPAR
       do I = 1, IPAR
          if (SLH(I,J) .eq. 0._r8) SLH(I,J) = 1.e-30_r8
       end do
    end do
    if (verbose) call gotData('2da','Sfc. lat. heat flux (SLH)')


    !// Mean Sealevel Pressure (Diagnostic only - not in old 19-layer) -----
    !// --------------------------------------------------------------------
    !// Unit: Pa -> hPa
    call get_netcdf_var_2d(FNAME, 'MSLP',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    MSLP(:,:) = R8XY(:,:) * 1.e-2_r8
    if (verbose) call gotData('2di','MSLP')


    !// Boundary Layer Height BLH ------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m
    !// Here we first set BLH, and then BLH_CUR=BLH.
    !// Second we rerad BLH_NEXT from the next time step, if it exists. If
    !// it does not exist, BLH_NEXT is also set to BLH.
    call get_netcdf_var_2d(FNAME, 'BLH',W2D, IPARW, JPARW)
    call TRUNG8(W2D, BLH, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    !// Average BLH at poles
    do J = 1, JPAR, JPAR-1
       DMASS = 0._r8
       do I = 1, IPAR
          DMASS = DMASS + BLH(I,J)
       end do
       DMASS = DMASS / real(IPAR, r8)
       do I = 1, IPAR
          BLH(I,J) = DMASS
       end do
    end do
    BLH_CUR(:,:) = BLH(:,:)
    if (verbose) call gotData('2di','Boundary Layer Height (BLH)')

    !// Check if file exist for next time step
    inquire(FILE=trim(FNAME_NEXT), exist=fex)
    if (.not. fex ) then
       write(6,'(a)') f90file//':'//subr// &
            ': No such file: '//trim(FNAME_NEXT)// &
            ', using BLH_NEXT=BLH_CUR'
       BLH_NEXT(:,:) = BLH_CUR(:,:)
    else
       call get_netcdf_var_2d(FNAME_NEXT, 'BLH',W2D, IPARW, JPARW)
       call TRUNG8(W2D, BLH_NEXT, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
       !// Average BLH at poles
       do J = 1, JPAR, JPAR-1
          DMASS = 0._r8
          do I = 1, IPAR
             DMASS = DMASS + BLH_NEXT(I,J)
          end do
          DMASS = DMASS / real(IPAR, r8)
          do I = 1, IPAR
             BLH_NEXT(I,J) = DMASS
          end do
       end do
       if (verbose) call gotData('2di','Boundary Layer Height (BLH_NEXT)')
    end if

    !// Find current level completely within BLH
    do J = 1, JPAR
       do I = 1, IPAR
          !// Find model level for this
          !// It is at least level 1
          LL = 1
          ZTOP = ZOFLE(2,I,J) - ZOFLE(1,I,J)
          ZBOT = 0._r8
          !// Check if it is higher up
          do L = 2, LPAR
             ZBOT = ZTOP !ZOFLE(L,I,J) - ZOFLE(1,I,J)
             ZTOP = ZOFLE(L+1,I,J) - ZOFLE(1,I,J)
             ZMID = 0.5_r8 * (ZTOP + ZBOT)
             if (ZTOP .gt. BLH(I,J)) then
                LL = L - 1
                exit
             end if
          end do
          LBLH(I,J) = LL
          if (LL .eq. 0) then
             write(6,'(a)') f90file//':'//subr// &
                  'LBLH is 1!!!!',BLH(I,J),ZBOT,ZMID,ZTOP
             stop
          end if
       end do
    end do


    !// 10m U wind component U10M ------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m/s
    call get_netcdf_var_2d(FNAME, 'U10M',W2D, IPARW, JPARW)
    call TRUNG8(W2D, SFU, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','10m U-wind (SFU)')


    !// 10m V wind component V10M ------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m/s
    call get_netcdf_var_2d(FNAME, 'V10M',W2D, IPARW, JPARW)
    call TRUNG8(W2D, SFV, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','10m V-wind (SFV)')


    !// 2m Temperature T2M -------------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: K
    call get_netcdf_var_2d(FNAME, 'T2M',W2D, IPARW, JPARW)
    call TRUNG8(W2D, SFT, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (minval(SFT) .eq. 0._r8) then
       write(6,'(a)')  f90file//':'//subr//': SFT zero: VERY WRONG!'
       stop
    end if
    if (verbose) call gotData('SFC','2m Temperature (SFT)')


    !// 2m Dewpoint Temperature TDEW2M -------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: K
    !// Read in Td, convert to mixing ratio using Clausius-Clapeyron
    !// SFQ is surface specific humidity.
    call get_netcdf_var_2d(FNAME, 'TDEW2M',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    do J = 1, JPAR
       do I = 1, IPAR
          SFTD = R8XY(I,J)
          !// Stull (1988):
          LV = (2.501_r8 - (2.37e-3_r8 * (SFT(I,J) - 273.16_r8))) * 1.e6_r8
          ESAT = 610.78_r8 * EXP( (-LV/R_H2O) &
                                  * ((1._r8/SFTD) - (1._r8/273.16_r8)) )
          !// The equation above is not valid below 0C, and should
          !// rather be calculated from another equation.
          !// The vapor pressure at T is equal to saturation pressure
          !// at Td.
          SFQ(I,J) = (R_AIR/R_H2O) * (ESAT/(P(I,J) * 100._r8))
       end do
    end do
    if (verbose) call gotData('2di','2m dew point -> SFQ')


    !// E/W Surface Stress EWSS (accumulated) ------------------------------
    !// --------------------------------------------------------------------
    !// Unit: N/m2, accumulated (N/m2*s)
    call get_netcdf_var_2d(FNAME, 'EWSS',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    EWSS(:,:) = R8XY(:,:) * ZDT
    if (verbose) call gotData('2da','E/W sfc stress (EWSS)')


    !// N/S Surface Stress NSSS (accumulated) ------------------------------
    !// --------------------------------------------------------------------
    !// Unit: N/m2, accumulated (N/m2*s)
    call get_netcdf_var_2d(FNAME, 'NSSS',W2D, IPARW, JPARW)
    call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    NSSS(:,:) = R8XY(:,:) * ZDT
    if (verbose) call gotData('2da','N/S sfc stress (NSSS)')


    !// Forecast Albedo ----------------------------------------------------
    !// --------------------------------------------------------------------
    !// Unit: fraction
    call get_netcdf_var_2d(FNAME, 'SA',W2D, IPARW, JPARW)
    call TRUNG8(W2D, SA, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','Forecast Albedo (SA)')


    !// --------------------------------------------------------------------
    !// Read some additional fields for DUST module
    !// --------------------------------------------------------------------

    !// Surface Solar Radiation Downwards (SSRD) ---------------------------
    !// --------------------------------------------------------------------
    !// Unit: W/m2, accumulated (W/m2*s)
    if (LDUST) then
       call get_netcdf_var_2d(FNAME, 'SSRD',W2D, IPARW, JPARW)
       !// Limit to positive values just in case
       W2D(:,:) = max(W2D(:,:), 0._r8)
       call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
       R8XY(:,:) = R8XY(:,:) * ZDT
       call dust_set_ssrd(SWVL1)
       if (verbose) call gotData('2da','Surface Solar Radiation Downwards (SSRD)')
    else
       !// Surface Solar Radiation Downwards (SSRD) (accumulated)
       if (verbose) call skipData('SFC','Surface Solar Radiation Downwards (SSRD)')
    end if

    !// Surface Thermal Radiation Downwards (STRD) ---------------------------
    !// --------------------------------------------------------------------
    !// Unit: W/m2, accumulated (W/m2*s)
    if (LDUST) then
       call get_netcdf_var_2d(FNAME, 'STRD',W2D, IPARW, JPARW)
       !// Limit to positive values just in case
       W2D(:,:) = max(W2D(:,:), 0._r8)
       call TRUNG8(W2D, R8XY, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
            JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
       R8XY(:,:) = R8XY(:,:) * ZDT
       call dust_set_strd(R8XY)
       if (verbose) call gotData('2da','Surface Thermal Radiation Downwards (SSRD)')
    else
       !// Surface Thermal Radiation Downwards (STRD) (accumulated)
       if (verbose) call skipData('SFC','Surface Thermal Radiation Downwards (STRD)')
    end if

    !// Volumetric soil water layers          ------------------------------
    !// --------------------------------------------------------------------
    !// SWVL1 - 0-7 cm
    !// SWVL2 - 7-28 cm
    !// SWVL3 - 28-100 cm 
    !// SWVL4 - 100-255 cm

    !// Volumetric soil water layer 1 (SWVL1) ------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m3/m3
    call get_netcdf_var_2d(FNAME, 'SWVL1',W2D, IPARW, JPARW)
    !// Limit to positive values just in case
    W2D(:,:) = max(W2D(:,:), 0._r8)
    call TRUNG8(W2D, SWVL1, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (LDUST) then
       !// DUST code has option to calculate its own SWVL1; set it here
       call dust_set_SWVL1(R8XY)
    end if
    if (verbose) call gotData('2di','Volumetric soil water layer 1 (SWVL1)')

    !// Volumetric soil water layer  (SWVL3) ------------------------------
    !// --------------------------------------------------------------------
    !// Unit: m3/m3
    call get_netcdf_var_2d(FNAME, 'SWVL3',W2D, IPARW, JPARW)
    !// Limit to positive values just in case
    W2D(:,:) = max(W2D(:,:), 0._r8)
    call TRUNG8(W2D, SWVL3, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','Volumetric soil water layer 3 (SWVL3)')


    !// Soil Temperature Level 1 (STL1 - ups, it is called SLT1 on file)
    !// Unit: K
    call get_netcdf_var_2d(FNAME, 'SLT1',W2D, IPARW, JPARW)
    call TRUNG8(W2D, STL1, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
         JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    if (verbose) call gotData('2di','Soil temp. lev. 1 (SLT1)')


    !// Other fields which may or may not be available ---------------------
    !// --------------------------------------------------------------------
    !// Snow Albedo (ASN)
    if (verbose) call skipData('SFC','Snow Albedo (ASN)')
    !// Snow Density (RSN)
    if (verbose) call skipData('SFC','Snow Density (RSN)')
    !// Sea Surface Temperature-Kelvin (SSTK)
    if (verbose) call skipData('SFC','SSTK')
    !// Ice surface temperature Layer 1 (ISTL1)
    if (verbose) call skipData('SFC','ISTL1')
    !// Ice surface temperature Layer 2 (ISTL2) 
    if (verbose) call skipData('SFC','ISTL2')
    !// Ice surface temperature Layer 3 (ISTL3) 
    if (verbose) call skipData('SFC','ISTL3')
    !// Ice surface temperature Layer 4 (ISTL4) 
    if (verbose) call skipData('SFC','ISTL4')
    !// Volumetric soil water layer 2 (SWVL2)
    if (verbose) call skipData('SFC','SWVL2')
    !// Volumetric soil water layer 3 (SWVL3)
    !// if (verbose) call skipData('SFC','SWVL3')
    !// Volumetric soil water layer 4 (SWVL4)
    if (verbose) call skipData('SFC','SWVL4')
    !// 10m wind gust
    if (verbose) call skipData('SFC','10m wind gust')
    !// Large-scale precipitation fraction (accumulated)
    if (verbose) call skipData('SFC','Large-scale prec. frac.')
    !// Downward UV radiation at the surface
    if (verbose) call skipData('SFC','Downw. UV rad. sfc')
    !// Convective available potential energy
    if (verbose) call skipData('SFC','CAPE')
    !// Total column liquid water
    if (verbose) call skipData('SFC','Tot. col. liq. water')
    !// Total column ice water
    if (verbose) call skipData('SFC','Tot. col. ice water')
    !// Surface Geopotential (Z)
    if (verbose) call skipData('SFC','Sfc. geopot. (Z)')
    !// Total Column Water (TCW)
    if (verbose) call skipData('SFC','Tot. col. water')
    !// Total Column Water Vapor (TCWV)
    if (verbose) call skipData('SFC','Tot. col. water vap.')
    !// Surface stress (eller Charnock)
    if (verbose) call skipData('SFC','Sfc stress/Charnock')
    !// Boundary Layer Dissipation (BLD) (accumulated)
    if (verbose) call skipData('SFC','Boundary Layer Dissipation')
    !// Total Cloud Cover
    if (verbose) call skipData('SFC','Tot. cloud cover')
    !// Soil Temperature Level 2 (STL2)
    if (verbose) call skipData('SFC','Soil temp. lev. 2')
    !// Land/Sea mask (LSM), NOTE: {0,1}
    !// Not so useful; it is 1 if land fraction > 0.5.
    if (verbose) call skipData('SFC','Land/sea mask')
    !// Surface Thermal Radiation Downwards (STRD) (accumulated)
    if (verbose) call skipData('SFC','Sfc. thermal rad. downwards')
    !// Surface Solar Radiation (SSR) (accumulated)
    if (verbose) call skipData('SFC','Sfc. solar rad.')
    !// Surface Thermal Radiation (STR) (accumulated)
    if (verbose) call skipData('SFC','Sfc. thermal rad.')
    !// Top Solar Radiation (TSR) (accumulated)
    if (verbose) call skipData('SFC','Top Solar Radiation (TSR)')
    !// Top Thermal Radiation (TTR) (accumulated)
    if (verbose) call skipData('SFC','Top Thermal Radiation (TTR)')
    !// Evaporation (E) (accumulated)
    if (verbose) call skipData('SFC','Evaporation')
    !// Soil Temperature Level 3 (STL3)
    if (verbose) call skipData('SFC','Soil temp. lev. 3')
    !// Soil Wetness Level 3 (SWL3)
    if (verbose) call skipData('SFC','Soil wetness level 3')
    !// Convective Cloud Cover (CCC)
    if (verbose) call skipData('SFC','Convective cloud cover')
    !// Low Cloud Cover (LCC)
    if (verbose) call skipData('SFC','Low cloud cover')
    !// Medium Cloud Cover (MCC)
    if (verbose) call skipData('SFC','Medium cloud cover')
    !// High Cloud Cover (HCC)
    if (verbose) call skipData('SFC','High cloud cover')
    !// Sun Shine Duration (SUND) (accumulated)
    if (verbose) call skipData('SFC','Sunshine duration')
    !// Latitudinal Gravity Wave Stress (LGWS) (accumulated)
    if (verbose) call skipData('SFC','LGWS')
    !// Meridional Gravity Wave Stress (MGWS)
    if (verbose) call skipData('SFC','MGWS')
    !// Gravity Wave Dissipation (GWD)
    if (verbose) call skipData('SFC','GWD')
    !// Skin Reservoir Content (SRC)
    if (verbose) call skipData('SFC','Skin Reservoir Content (SRC)')
    !// Maximum Temperature at 2m (MX2T)
    if (verbose) call skipData('SFC','MX2T')
    !// Minimum Temperature at 2m (MN2T)
    if (verbose) call skipData('SFC','MN2T')
    !// RunOff (RO)
    if (verbose) call skipData('SFC','RunOff')
    !// Top Net Solar Radiation Clear Sky (TSRC) (accumulated)
    if (verbose) call skipData('SFC','TSRC')
    !// Top Net Thermal Radiation Clear Sky (TTRC) (accumulated)
    if (verbose) call skipData('SFC','TTRC')
    !// Surface net solar radiation, clear sky (accumulated)
    if (verbose) call skipData('SFC','SNSRCS')
    !// Surface net thermal radiation, clear sky (accumulated)
    if (verbose) call skipData('SFC','SNTRCS')
    !// Log Surface roughness for Heat
    if (verbose) call skipData('SFC','Log(sfc. roughness for heat)')
    !// Skin Temperature (SKT)
    if (verbose) call skipData('SFC','Skin Temperature (SKT)')
    !// Soil Temperature Level 4 (STL4)
    if (verbose) call skipData('SFC','Soil temp. lev. 4')
    !// Temperature of Snow Layer (TSN)
    if (verbose) call skipData('SFC','Temperature of snow layer')
    !// Forecast Surface Roughness (FSR) 
    if (verbose) call skipData('SFC','Forecast sfc. roughness')
    !// Forecast Log(Surface Roughness Heat)
    if (verbose) call skipData('SFC','Forecast log(Sfc. roughness heat)')





    !// Additional calculations
    !//---------------------------------------------------------------------

    !// USTAR from Surface Stress NSSS and EWSS
    do J = 1, JPAR
      do I = 1, IPAR
        SMF(I,J) = sqrt(NSSS(I,J)*NSSS(I,J) + EWSS(I,J)*EWSS(I,J))
        !// ideal gas density [kg/m3]
        SDEN      = M_AIR * 100._r8*P(I,J) / (R_UNIV * 1.e3_r8 * SFT(I,J))
        USTR(I,J) = sqrt(SMF(I,J)/SDEN)
        if (USTR(I,J) .eq. 0._r8) then
          !// When zero, override with minimum value
          USTR(I,J) = max(5.e-3_r8, USTR(I,J))
          !// Also override SMF
          SMF(I,J) = (USTR(I,J)**2) * SDEN
        end if
      end do
    end do

    !// Filter fluxes for wrong signs etc.
    call fluxfilter2(cdetu, cdetd)

    !// Put low level limits to CLDFR, CLDLWC, and CLDIWC, and set junk
    !// values to zero.
    do L = 1, LWEPAR
      do J = 1, JPAR
        do I = 1, IPAR
          if (CLDFR(I,J,L) .gt. EPS) then
            CLDLWC(I,J,L) = CLDLWC(I,J,L) / CLDFR(I,J,L)
            CLDIWC(I,J,L) = CLDIWC(I,J,L) / CLDFR(I,J,L)

            call CIWMIN(T(I,J,L), CLDLWC(I,J,L), CLDIWC(I,J,L))

            CLDLWC(I,J,L) = CLDLWC(I,J,L) * CLDFR(I,J,L)
            CLDIWC(I,J,L) = CLDIWC(I,J,L) * CLDFR(I,J,L)
          else
            CLDFR(I,J,L) = 0._r8
            call CFRMIN(T(I,J,L), CLDFR(I,J,L), CLDLWC(I,J,L), &
                        CLDIWC(I,J,L), EPS)
          end if
        end do
      end do
    end do

    !// Done reading all meteorological data
    deallocate( VTMP, UTMP, W3Da, W3Db, ZOFLEW, QW, TW, R8XYZ, &
         CLDFRW, CLDIWCW, CLDLWCW, &
         CDETU, CDETD, &
         PW, EWSS, NSSS, W2D, R8XY )
    !deallocate( LSPREC, CNVPREC )

    !//---------------------------------------------------------------------
  end subroutine update_metdata
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine fluxfilter2(cdetu,cdetd)
    !//---------------------------------------------------------------------
    !// Filter convective mass fluxes.
    !// Based on filtering in p-wind_ec.f.
    !//
    !// Ole Amund Sovde, March 2015, March 2010
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LWEPAR, LWDPAR
    use cmn_met, only: CENTU, CENTD, CWETE, CWETD
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(inout) :: cdetu(ipar,jpar,lwepar)
    real(r8), intent(inout) :: cdetd(ipar,jpar,lwdpar)

    !// Local variables
    real(r8) :: sden
    integer :: KCTOP(ipar,jpar), KDBOT(ipar,jpar), i,j,l,LL
    !//---------------------------------------------------------------------

    !// CWETE and CWETD are already filtered for wrong signs.
    !// There may still be some places where a grid box have
    !// flux, but there is none above or below. This is a bit strange,
    !// but fairly ok; such a small convective flux will then be
    !// detrained above. No need to filter these away.

    !// CDETU and CDETD are also already filtered for wrong signs.


    !// More CDETU and CDETD filters

    !// Find top of convection
    KCTOP(:,:) = 1
    do J = 1, JPAR
       do I = 1, IPAR
          do L = LWEPAR-1,1,-1
             !// Loop downwards
             if (CWETE(i,j,l) .gt. 0._r8) then
                KCTOP(i,j) = L
                exit
             else
                !// No flux; no detrainment
                cdetu(i,j,l) = 0._r8
             end if
          end do
       end do
    end do

    !// Find bottom of downdrafts
    KDBOT(:,:) = LWDPAR-1
    do J = 1,JPAR
       do I = 1, IPAR
          do L = 1, LWDPAR-1
             !// Loop downwards
             if (cwetd(i,j,l+1) .lt. 0._r8) then
                !// Mass comes in at top
                KDBOT(i,j) = L
                exit
             else
                !// No flux in at top; no detrainment
                cdetd(i,j,l) = 0._r8
             end if
          end do
       end do
    end do

    !// Detrainment cannot be larger than flux in at bottom
    do J = 1,JPAR
       do I = 1, IPAR
          do l = 1, KCTOP(I,J)
             !// Loop downwards; cdetu and cwete are non-negative
             !// due to filtering above
             if (cdetu(i,j,l) .gt. cwete(i,j,l)) then
                !// CDETU cannot be larger than the flux in!
                if (cwete(i,j,l) .eq. 0._r8) then
                   CDETU(I,J,L) = 0._r8
                else
                   CDETU(i,j,l) = CWETE(i,j,l)
                end if
             end if
          end do
       end do
    end do

    !// Detrainment cannot be larger than flux in at top
    do J = 1,JPAR
       do I = 1, IPAR
          !// Make sure top is zero
          CDETD(i,j,LWDPAR) = 0._r8
          do L = LWDPAR-1, KDBOT(I,J), -1
             !// Loop downwards; cdetd is non-negative and cwetd is
             !// non-positive due to filtering above
             if (CDETD(i,j,l) .gt. -CWETD(i,j,l+1)) then
                if (CWETD(i,j,l+1) .lt. 0._r8) then
                   CDETD(i,j,l) = -CWETD(i,j,l+1)
                else
                   CDETD(i,j,l) = 0._r8
                end if
             end if
          end do
       end do
    end do

    !// FINALLY
    !// Build entrainment - updrafts (detrainment is positive)
    do J = 1,JPAR
       do I = 1, IPAR
          do L = 1, KCTOP(I,J)
             !// F(L) + E = F(L+1) + D
             sden = cwete(i,j,l+1) - cwete(i,j,l) + cdetu(i,j,l)
             !// Only allow positive entrainment
             centu(i,j,l) = max(0._r8, sden)
          end do
       end do
    end do
    !// Build entrainment - downdrafts (flux is negative)
    do L = 1, LWDPAR-1
       do J = 1, JPAR
          do I = 1, IPAR
             !// -F(L+1) + E = -F(L) + D
             sden = cwetd(i,j,l+1) - cwetd(i,j,l) + cdetd(i,j,l)
             !// Only allow positive entrainment
             centd(i,j,l) = max(0._r8, sden)
          end do
       end do
    end do

    !//---------------------------------------------------------------------
  end subroutine fluxfilter2
  !//-----------------------------------------------------------------------




  !//-----------------------------------------------------------------------
  subroutine data2mpblocks(r8data, mpdata)
    !//---------------------------------------------------------------------
    !// Puts real*8 data of size (ipar,jpar) into real*8 MP-block structure
    !// mpdata(idblk,jdblk,mpblk).
    !//
    !// Ole Amund Sovde, April 2013
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: r8data(ipar,jpar)
    !// Output
    real(r8), intent(out) :: mpdata(idblk,jdblk,mpblk)
    !// Locals
    integer :: I, J, II, JJ, MP
    !//---------------------------------------------------------------------
    do MP = 1, MPBLK
      !// Loop over latitude (J is global, JJ is block)
      do J = MPBLKJB(MP),MPBLKJE(MP)
        JJ    = J - MPBLKJB(MP) + 1
        !// Loop over longitude (I is global, II is block)
        do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Change structure
          mpdata(II,JJ,MP) = r8data(I,J)
        end do
      end do
    end do
    !//---------------------------------------------------------------------
  end subroutine data2mpblocks
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine gotData(TYP,LABEL)
    !//---------------------------------------------------------------------
    !// Print info about 2D field read.
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=3), intent(in) :: TYP
    character(len=*), intent(in) :: LABEL
    !//---------------------------------------------------------------------
    write(6,'(a,1x,a)') ' update_metdata: Read '//TYP//':   ', &
         trim(LABEL)
    !//---------------------------------------------------------------------
  end subroutine gotData
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  subroutine skipData(TYP,LABEL)
    !//---------------------------------------------------------------------
    !// Print info about 2D field skipped.
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    character(len=3), intent(in) :: TYP
    character(len=*), intent(in) :: LABEL
    !//---------------------------------------------------------------------
    write(6,'(a,i5,1x,a)') ' update_metdata: Skipped '//TYP//':', &
         trim(LABEL)
    !//---------------------------------------------------------------------
  end subroutine skipData
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
end module metdata_ecmwf
!//=========================================================================

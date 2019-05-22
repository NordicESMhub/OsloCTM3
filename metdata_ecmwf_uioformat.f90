!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routine for reading meteorological data from ECMWF, in UiO EC-files.
!//=========================================================================
module metdata_ecmwf
  !//-----------------------------------------------------------------------
  !// MODULE: metdata_ecmwf
  !// DESCRIPTION: Routine for reading meteorological data from ECMWF,
  !//              stored on UiO EC-file format.
  !//
  !// Contains
  !//   subroutine update_metdata
  !//   subroutine CFRMIN
  !//   subroutine CIWMIN
  !//   subroutine fluxfilter2
  !//   subroutine r4data2mpblocks
  !//-----------------------------------------------------------------------
  use cmn_size, only: LPAR
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  integer :: LMAP(LPAR+1)
  !//-----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'metdata_ecmwf.f90'
  !//-----------------------------------------------------------------------
  private
  public update_metdata
  save
  !//-----------------------------------------------------------------------

contains


  !//-----------------------------------------------------------------------
  subroutine update_metdata(ILOOP, DTMET, NMET)
    !//---------------------------------------------------------------------
    !// Read input data. ECMWF data has three files, one for spectral 
    !// data (.b01), one for 2-d gridpoint data (.b03) and one for 
    !// 3-d gridpoint data (.b02).
    !// Spectral transformation is done in this routine, and sumation of
    !// physical fields from Txx ---> T42, T21 or other resolution.
    !// Allow several different formats for 2-d gridpoint files by labelling
    !// the record numbers we need (we need N15R records)
    !//
    !// This routine assumes meteorology is to be updated each NMET. E.g.
    !// ERA-40 had 6-hour fields, and was tweaked to be read every other
    !// NMET, but with the double time step. This has been removed for now.
    !//
    !// Ole Amund Sovde, April 2015
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, r4
    use cmn_size, only: IPARW, JPARW, LPARW, IPAR, JPAR, IDGRD, JDGRD, &
         LPAR, LWEPAR, LWDPAR, &
         FFTPTS, NRFFT, NRNM, NMMAX, NTRUNW, LOSLOCHEM, LDUST
    use cmn_ctm, only: JYEAR, JMON, JDAY, GMTAU, LMMAP, XLMMAP, TMET, &
         ALP, ALPV, TRIG, IFAX, LDEG, ZDEGI, ZDEGJ, IMAP, JMAP, &
         XYZA, XYZB, IMEPZ, DD, SS, ETAAW, ETABW, ETAA, ETAB, &
         WGLYE, WGLYG, YDGRD, AREAXYW, AREAXY, PLAND
    use cmn_met, only: P, U, V, T, Q, CWETN, CWETE, CWETD, CENTU, CENTD, &
         PRECCNV, PRECLS, ZOFLE, &
         CLDFR, CLDLWC, CLDIWC, &
         SLH, SHF, SMF, SFT, SFQ, SFU, SFV, BLH, SA, &
         PVU, UMS, VMS, CI, SD, SMLT, ES, SNFL, PhotActRad, MSLP, USTR, &
         PMEAN, PMEANW, MPATH1, MPATH2, MFILE3, MYEAR, &
         metTYPE, metCYCLE, metREVNR, HnativeRES, VRESW
    use cmn_parameters, only: M_AIR, R_UNIV, R_AIR, R_H2O, A0, CPI
    use cloudjx, only: SWSTORE, CLDSTORE, TYPSTORE, RAN4, IRAN0, &
                       LCLDAVG, LCLDQMD, LCLDQMN, LCLDRANA, LCLDRANQ
    use regridding, only: TRUNG4, TRUNG8
    use utilities, only: ctmExitC, CFRMIN, CIWMIN
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
         PT, PB, &
         DELZ, LWC, IWC, TFACT, QMIN, &
         DMASS, ZAREA, PSRF, &
         POFLE(LPAR+1)

    !// Indices
    integer :: I,J,II,JJ,L,LL
    integer :: NF, NrOf3Dfields, NrOf2Dfields, ioerr

    !// To be read from file
    integer :: FLD_ID, ISEC(16)
    logical :: fex

    logical :: LUV  !// Flag if U or V are calculated by SPE2GP
    logical :: LPVU !// Need to calculate PVU?

    character(len=120) :: FNAME1, FNAME2, FNAME3
    character(len=4) :: CYYx

    !// To check if RAIN is split into CONVRAIN/LSRAIN or not
    !// If 217 and 218 are found instead of 216 we have CONVRAIN/LSRAIN:
    logical :: found216                !// field 216 exist; only total rain
    real(r8) :: CNVFRAC, CONVTOT, RAINTOT


    !//---------------------------------------------------------------------
    !// Allocatable arrays - single precision
    real(r4), allocatable :: &
         GRDATA2(:,:), GRDATA2HI(:,:), SPWK5(:,:), SPWK6(:,:)
    real(r4), allocatable :: &
         GRDATA3(:,:,:), GRDATA3HI(:,:,:)
    real(r4), allocatable :: &
         PSPEHI(:)

    !// Allocatable arrays - double precision
    real(r8), allocatable :: &
         FFTAR(:,:), WORK(:), &
         VTMP(:,:), UTMP(:)
    real(r8), dimension(:,:), allocatable :: &
         PW, EWSS, NSSS, LSPREC, CNVPREC, R8XY
    real(r8), dimension(:,:,:), allocatable :: &
         TW, QW, CLDFRW, CLDIWCW, CLDLWCW, CDETU, CDETD, &
         ZOFLEW, WK3DT, WK3DU, WK3DV, UMSW, VMSW, &
         TMPCRAIN, TMPLRAIN
    !//---------------------------------------------------------------------
    character(len=*), parameter :: subr = 'update_metdata'
    !//---------------------------------------------------------------------

    !// Allocate spectral arrays
    allocate( FFTAR(FFTPTS,NRFFT), WORK(FFTPTS*NRFFT), &
              SPWK5(NRNM,LPARW), SPWK6(NRNM,LPARW), PSPEHI(NRNM) )

    !// Allocate 3D arrays - native resolution
    allocate( VTMP(IPARW,JPARW), UTMP(IPARW), &
              WK3DU(IPARW,JPARW,LPARW), WK3DV(IPARW,JPARW,LPARW), &
              WK3DT(IPARW,JPARW,LPARW), &
              ZOFLEW(LPAR+1,IPARW,JPARW), QW(IPARW,JPARW,LPAR), &
              TW(IPARW,JPARW,LPAR), CLDFRW(IPARW,JPARW,LWEPAR), &
              CLDIWCW(IPARW,JPARW,LWEPAR), CLDLWCW(IPARW,JPARW,LWEPAR), &
              GRDATA3HI(IPARW,JPARW,LPARW), &
              UMSW(IPARW,JPARW,LPARW), VMSW(IPARW,JPARW,LPARW) )

    !// Allocate 3D arrays - window resolution (IPAR/JPAR)
    allocate( TMPCRAIN(IPAR,JPAR,LPARW), TMPLRAIN(IPAR,JPAR,LPARW), &
              CDETU(IPAR,JPAR,LWEPAR), CDETD(IPAR,JPAR,LWDPAR), &
              GRDATA3(IPAR,JPAR,LPARW) )

    !// Allocate 2D arrays - native resolution
    allocate( PW(IPARW,JPARW), GRDATA2HI(IPARW,JPARW) )

    !// Allocate 2D arrays - window resolution (IPAR/JPAR)
    allocate( LSPREC(IPAR,JPAR),CNVPREC(IPAR,JPAR), &
              EWSS(IPAR,JPAR), NSSS(IPAR,JPAR), &
              GRDATA2(IPAR,JPAR), R8XY(IPAR,JPAR) )
    !//---------------------------------------------------------------------

    !// Initialize
    found216 = .false.
    LUV  = .false.

    !// Time step for meteorological data is DTMET
    ZDT    = 1._r8 / DTMET

    !//Define minimum humidity QMIN in kg/kg
    QMIN = 3.e-6_r8 * 18._r8 / M_AIR

!//UPDATE
    !// locate the position of random number sequence based on year/day/hour
    IRAN0 = 1 + 3*(JYEAR - 1900) + 7*JDAY + 11*nint(GMTAU)


    if (metTYPE(1:10) .eq. 'ECMWF_oIFS') then
       !// ECMWF openIFS EC-format
       write(MPATH2(1:4),'(i4.4)') MYEAR
       write(MFILE3(7:13),'(i4.4,A3)') MYEAR,TMET
       if (metREVNR .eq. 1) then
          write(CYYx(1:4),'(a1,i2.2,a1)') 'C',metCYCLE,'a'
       else
          write(6,'(a)')     '*** Unknown metdata revision number!'
          write(6,'(a,2i4)') '    metCYCLE / metREVNR: ',metCYCLE,metREVNR
          stop
       end if

       FNAME3 = trim(MFILE3)
       FNAME1 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b01'
       FNAME2 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b02'
       FNAME3 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b03'

    else if (trim(metTYPE) .eq. 'ECMWF_IFS') then

       !// Old style metdata names
       write(MPATH2(1:4),'(i4)') MYEAR
       write(MFILE3(3:5),'(A3)') TMET
       FNAME3 = trim(MFILE3)
       FNAME1 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b01'
       FNAME2 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b02'
       FNAME3 = trim(MPATH1)//trim(MPATH2)//trim(FNAME3)//'.b03'
    else
       write(6,'(a)') f90file//':'//subr// &
            ': Not set up for metTYPE: '//trim(metTYPE)
       stop
    end if

    !//---------------------------------------------------------------------
    !// Initial step - setup dry airmass (P() and Q())
    !//---------------------------------------------------------------------
    if (ILOOP .le. 0) then

       !// Set up level weightings if vertical resolution degraded
       do LL = LPARW+1, 1, -1
          L   = LMMAP(LL)
          LMAP(L) = LL
       end do

       !//print out info, open files and read in data, transform and sum up
       write(6,'(a)') 'Initializing meteorological data'
       write(6,'(2x,a)') trim(FNAME1)
       write(6,'(2x,a)') trim(FNAME2)

       inquire(FILE=trim(FNAME1), exist=fex)
       if (.not. fex) call ctmExitC('update_metdata: No such file: '//trim(FNAME1))
       close(13)
       open(13,FILE=trim(FNAME1),FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
       if (ioerr .ne. 0) call ctmExitC('update_metdata: ERROR on .b01 file')

       inquire(FILE=trim(FNAME2), exist=fex)
       if (.not. fex) call ctmExitC('update_metdata: No such file: '//trim(FNAME2))
       close(14)
       open(14,FILE=trim(FNAME2),FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
       if (ioerr .ne. 0) call ctmExitC('update_metdata: ERROR on .b02 file')

       !// Pressure field
       read(13) FLD_ID
       if (FLD_ID .ne. 152) &
            call ctmExitC('update_metdata: ERROR, input Data PSPEHI')
       read(13) PSPEHI

       !// Pressure; transform from spectral to grid point data
       call SPE2GP(PSPEHI, NRNM, LUV, PW, IPARW, JPARW, 1, ALP, NMMAX, &
                   FFTAR, WORK, TRIG, IFAX, (NTRUNW+1), FFTPTS, NRFFT)

       !// save high-res pressure Log(Ps) --> hPa
       do J = 1, JPARW
          do I = 1, IPARW
             PW(I,J) = EXP(PW(I,J)) * 1.e-2_r8
          end do
       end do
       if (ldeg) then
          !// Truncate to window resolution
          call TRUNG8(PW, P, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                      JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
       else
          do J = 1,JPAR
            do I = 1,IPAR
              P(I,J) = PW(I,J)
            end do
          end do
       end if


       !// Temperature (originally not in the UCI code)
       read(13) FLD_ID
       if (FLD_ID .ne. 130) &
            call ctmExitC('update_metdata: ERROR, input Data T')
       read(13) SPWK5
       !// Temperature; transform from spectral to grid point data
       call SPE2GP(SPWK5, NRNM, LUV, WK3DT, IPARW, JPARW, LPARW, ALP, &
                   NMMAX, FFTAR, WORK, TRIG, IFAX, (NTRUNW+1), FFTPTS, NRFFT)
       TW(:,:,:)  = 0._r8
       do L = 1,LPAR
          do LL = LMAP(L),LMAP(L+1)-1
            do J = 1,JPARW
              do I = 1,IPARW
                TW(I,J,L) = TW(I,J,L) + WK3DT(I,J,LL) * XLMMAP(LL)
              end do
            end do
          end do
       end do
       if (ldeg) then
          call TRUNG8(TW, T, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                      JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)
       else
          do L = 1,LPAR
            do J = 1,JPAR
              do I = 1,IPAR
                T(I,J,L) = TW(I,J,L)
              end do
            end do
          end do
       end if
       !// Test to check temperature field
       if (minval(WK3DT(:,:,1)) .le. 0.) then
          print*,'update_metdata: Temperature is wrong'
          print*,'MIN T',minval(WK3DT(:,:,1))
          print*,'Maybe ALP=0? min/max ALP',minval(alp), maxval(alp)
          stop
       end if

       !// Close file
       close(13)

       !// Water vapour
       read(14) FLD_ID
       if (FLD_ID .ne. 133)  &
            call ctmExitC('update_metdata: ERROR, input Data Q')

       read(14) GRDATA3HI
       !// Collapse layers
       QW(:,:,:)  = 0._r8
       do L = 1,LPAR
          do LL = LMAP(L),LMAP(L+1)-1
            do J = 1,JPARW
              do I = 1,IPARW
                QW(I,J,L) = QW(I,J,L) + GRDATA3HI(I,J,LL) * XLMMAP(LL)
              end do
            end do
          end do
       end do
       !// There may be some entries of small negative numbers
       if (minval(QW) .lt. 0._r8) then
          print*,'update_metdata: min QW:',minval(QW),', setting to',QMIN
          QW(:,:,:) = max(QW(:,:,:), QMIN)
       end if

       if (ldeg) then
          call TRUNG8(QW, Q, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                      JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR,1)
       else
          Q(:,:,:) = QW(:,:,:)
       end if

       !// Close file
       close(14)

       !// Polar cap filtering
       call EPZ_TQ(T, Q, XYZA, XYZB, P, IMEPZ,IPAR,JPAR,LPAR,IPAR,JPAR,LPAR)
       call EPZ_P(P, IMEPZ, IPAR,JPAR, IPAR, JPAR)

       !// Initialization is done - these steps are for all time steps
       return
    end if !// if (ILOOP .le. 0) then


    !//---------------------------------------------------------------------
    !// All time steps - update all meteorological fields
    !//---------------------------------------------------------------------
    if (verbose) write(6,'(a,i5)') f90file//':'//subr// &
         ': Reading new metdata JDAY: '//TMET//', NMET:',NMET
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

    !// New day? If so, open new files
    if (NMET .eq. 1) then

       inquire(FILE=trim(FNAME1), exist=fex)
       if (.not. fex) call ctmExitC('update_metdata: No such file: '//trim(FNAME1))
       close(13)
       open(13,FILE=trim(FNAME1),FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
       if (ioerr .ne. 0) call ctmExitC('update_metdata: ERROR on .b01 file')

       inquire(FILE=trim(FNAME2), exist=fex)
       if (.not. fex) call ctmExitC('update_metdata: No such file: '//trim(FNAME2))
       close(14)
       open(14,FILE=trim(FNAME2),FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
       if (ioerr .ne. 0) call ctmExitC('update_metdata: ERROR on .b02 file')

       inquire(FILE=trim(FNAME3), exist=fex)
       if (.not. fex) call ctmExitC('update_metdata: No such file: '//trim(FNAME3))
       close(15)
       open(15,FILE=trim(FNAME3),FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
       if (ioerr .ne. 0) call ctmExitC('update_metdata: ERROR on .b03 file')

    end if



    !//---------------------------------------------------------------------
    !// SPECTRAL DATA
    !//---------------------------------------------------------------------
    !// Read in spectral data - error check only on Pressure since it 
    !// is assumed that the ordering of fields is fine.

    !// Pressure field
    read(13) FLD_ID
    if (FLD_ID .ne. 152) &
         call ctmExitC('update_metdata: ERROR, input Data PSPEHI')
    read(13) PSPEHI

    !// Pressure; transform from spectral to grid point data
    call SPE2GP(PSPEHI, NRNM, LUV, PW, IPARW, JPARW, 1, ALP, &
                NMMAX, FFTAR, WORK, TRIG, IFAX, (NTRUNW+1), FFTPTS, NRFFT)
    !// save high-res pressure Log(Ps) --> hPa
    do J = 1,JPARW
       do I = 1,IPARW
          PW(I,J) = EXP(PW(I,J))*1.e-2_r8
       end do
    end do

    if (ldeg) then
       Call TRUNG8(PW, P, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                   JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
    else
       do J = 1, JPAR
         do I = 1, IPAR
           P(I,J) = PW(I,J)
         end do
       end do
    end if

    !// Temperature
    read(13) FLD_ID
    If (FLD_ID .ne. 130) &
         call ctmExitC('update_metdata: ERROR, input Data PSPEHI')

    read(13) SPWK5
    !// Temperature; transform from spectral to grid point data
    call SPE2GP(SPWK5, NRNM, LUV, WK3DT, IPARW, JPARW, LPARW, ALP, &
                NMMAX, FFTAR, WORK, TRIG, IFAX, (NTRUNW+1), FFTPTS, NRFFT)
    TW(:,:,:)  = 0._r8 !// Clear array
    do L = 1,LPAR
       do LL = LMAP(L),LMAP(L+1)-1
         do J = 1,JPARW
           do I = 1,IPARW
             TW(I,J,L) = TW(I,J,L) + WK3DT(I,J,LL) * XLMMAP(LL)
           end do
         end do
       end do
    end do
    if (ldeg) then
       call TRUNG8(TW, T, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                   JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)
    else
       do L = 1, LPAR
         do J = 1, JPAR
           do I = 1, IPAR
             T(I,J,L) = TW(I,J,L)
           end do
         end do
       end do
    end if


    !// Vorticity and Divergence - find U and V
    read(13) FLD_ID
    if (FLD_ID .ne. 138) &
         call ctmExitC('update_metdata: ERROR, input Data 138')
    read(13) SPWK5

    read(13) FLD_ID
    if (FLD_ID .ne. 155) &
         call ctmExitC('update_metdata: ERROR, input Data 155')
    read(13) SPWK6

    !// Vorticity and Divergence - find U and V,
    !// then do spectral to GP transform
    call ZD2UV(SPWK5, SPWK6, WK3DU, WK3DV, DD, SS, IPARW, JPARW, LPARW, &
               ALP, ALPV, NMMAX, NRNM, FFTAR, WORK, TRIG, IFAX, &
               NTRUNW, NRFFT, FFTPTS)
    !// scale from U=u*cos() to u
    do L = 1, LPARW
       DETA = ETAAW(L) - ETAAW(L+1)
       DETB = ETABW(L) - ETABW(L+1)
       do J = 1, JPARW
          BAND = A0 * (WGLYE(J+1) - WGLYE(J))
          ZCOS = 1._r8 / COS(WGLYG(J))
          do I = 1, IPARW
             UTMP(I) = WK3DU(I,J,L) * ZCOS * BAND
             !// Save center values for m/s
             UMSW(I,J,L) = WK3DU(I,J,L)*ZCOS !// Save U(m/s)
          end do
          !// Get edge values (unit conversion needs pressure on edge)
          do I = 2, IPARW
             DELP = DETA + DETB * (PW(I-1,J) + PW(I,J)) * 0.5_r8
             WK3DU(I,J,L) = (UTMP(I-1) + UTMP(I)) * 0.5_r8 * DELP
          end do
          DELP = DETA + DETB * (PW(1,J) + PW(IPARW,J)) * 0.5_r8
          WK3DU(1,J,L) = (UTMP(IPARW) + UTMP(1)) * 0.5_r8 * DELP
       end do
    end do

    !// scale from V=v*cos() to v (remember V is on edge already)
    do L = 1, LPARW
       DETA = ETAAW(L) - ETAAW(L+1)
       DETB = ETABW(L) - ETABW(L+1)
       VTMP(:,:) = 0._r8 !// Needed for getting VMSW
       !// WK3DV does not start at SP, but at the J=2. It means
       !// WK3DV is flux out of box J, not into J.
       !// Given all edge points, the number of boxes should be
       !// JPARW+1, covering 1:JPARW+1. But both the poles should
       !// have zero wind, so we only need JPARW-1 boxes.
       !//
       !// But in addition, the values at JPARW/2 and JPARW/2+1
       !// are duplicated, both representing wind across Equator.
       !// So WK3DV has the size JPARW instead of JPARW-1.
       !//
       !// A bit confusing, the CTM uses V as flux into J-box instead of
       !// out of it. This will be taken care of at the end.
       !// SH flux out of J-box:
       do J = 1, JPARW/2
          ZCOS = 1._r8 / COS(WGLYE(J+1))
          do I = 1, IPARW
             VTMP(I,J) = WK3DV(I,J,L) * ZCOS
          end do
       end do
       !// Now VTMP covers the V up to Equator, starting at the
       !// models index J=2, but stored in VTMP at J=1.
       !// NH flux out of J-box (JPARW/2 and JPARW/2+1 are duplicated):
       do J = JPARW/2 + 2, JPARW
          ZCOS = 1._r8 / COS(WGLYE(J))
          do I = 1, IPARW
             VTMP(I,J-1) = WK3DV(I,J,L) * ZCOS
          end do
       end do
       !// Thus, at J=JPARW, VTMP (with its indicing) represents NP
       !// and is therefore zero.

       !//change unit, m/s ==> Kg/s (100./G0 is done in pdyn0.f)
       !// Save center values for m/s J=1 (assume V=0 at SP)
       VMSW(:,1,L) = VTMP(:,1) * 0.5_r8
       !// Map the VTMP back to correct model edge grid (from J to J+1),
       !// i.e. converting from flux out of J-1 to flux into J.
       do J = 2, JPARW
          VEDGE  = 2._r8 * CPI * A0 * COS(WGLYE(J))/real(IPARW, r8)
          do I = 1, IPARW
             DELP = DETA + DETB * (PW(I,J-1) + PW(I,J)) * 0.5_r8
             !// Save center values for m/s J>1
             VMSW(I,J,L) = (VTMP(I,J) + VTMP(I,J-1)) * 0.5_r8
             WK3DV(I,J,L) = VTMP(I,J-1) * VEDGE * DELP
          end do
       end do
    end do

    if (ldeg) then
       do L = 1, LPAR
         do LL = LMAP(L),LMAP(L+1)-1
           do I = 1, IPAR
             II = IMAP(1,I)
             do J = 1, JPAR
               do JJ = 1, JDGRD
                 U(I,J,L) = U(I,J,L) + WK3DU(II,JMAP(JJ,J),LL)
               end do
             end do
           end do
           do J = 1, JPAR
             JJ = JMAP(1,J)
             do I = 1, IPAR
               do II = 1, IDGRD
                 V(I,J,L) = V(I,J,L) + WK3DV(IMAP(II,I),JJ,LL)
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
               U(I,J,L) = U(I,J,L) + WK3DU(I,J,LL)
               V(I,J,L) = V(I,J,L) + WK3DV(I,J,LL)
             end do
           end do
         end do
       end do
    end if

    !// VERTICAL VELOCITY [Pa/s] phony read in if present
    !// Some of the oldest metdata also contain W (FLD_ID 135).
    !// W is not used, so I have included a dummy read-in for that.
    !read(13, iostat=ioerr) FLD_ID
    !if (ioerr .eq. 0) then
    !   if (FLD_ID .ne. 135) then
    !      read(13) RSUM
    !   else
    !      backspace(13)
    !   end if
    !end if !// if (ioerr .eq. 0) then


    !//---------------------------------------------------------------------
    !// 3-d GRID POINT DATA
    !//---------------------------------------------------------------------
    !// Read in 3-d Gridpoint data, check on field type


    !// WATER VAPOUR
    read(14) FLD_ID
    if (FLD_ID .ne. 133)  &
         call ctmExitC('update_metdata: ERROR, input Data Q')

    read(14) GRDATA3HI
    !// Collapse layers
    QW(:,:,:) = 0._r8
    do L = 1, LPAR
      do LL = LMAP(L), LMAP(L+1)-1
        do J = 1, JPARW
          do I = 1, IPARW
            QW(I,J,L) = QW(I,J,L) + GRDATA3HI(I,J,LL) * XLMMAP(LL)
          end do
        end do
      end do
    end do
    if (minval(QW) .lt. 0._r8) then
       print*,'update_metdata: min QW:',minval(QW),', setting to',QMIN
       QW(:,:,:) = max(QW(:,:,:), QMIN)
    end if
    if (ldeg) then
       call TRUNG8(QW, Q, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                   JDGRD, IPARW, JPARW, IPAR, JPAR, LPAR, 1)
    else
       do L = 1, LPAR
         do J = 1, JPAR
           do I = 1, IPAR
             Q(I,J,L) = QW(I,J,L)
           end do
         end do
       end do
    end if


    !// Altitudes
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
           print*,'update_metdata: ZOFLE a',1,i,j,PMEAN(I,J),PMEANW(I,J)
           print*,pmean(:,j)
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
             print*,'update_metdata: ZOFLE b',l,i,j,delz,q(i,j,l-1),&
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


    if (LPARW .eq. 40) then
       if (MYEAR .lt. 2001) then
          NrOf3Dfields = 8
       else
          ! Splitting rain in large scale and convective
          NrOf3Dfields = 9
       end if
    else if (LPARW .eq. 60) then
       NrOf3Dfields = 9
       !// Note that ERA-40 has 8 3D fields
    end if


    do NF = 1, NrOf3Dfields

       read(14) FLD_ID, ISEC
       !// Read native resolution
       GRDATA3HI(:,:,:) = 0._r4
       read(14) (((GRDATA3HI(I,J,L),I=1,IPARW),J=1,JPARW),L=ISEC(13),ISEC(14))

       !// Degrade or not
       GRDATA3(:,:,:) = 0._r4
       if (ldeg) then
         call TRUNG4(GRDATA3HI, GRDATA3, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                     JDGRD,IPARW,JPARW,IPAR,JPAR,LPARW,1)
       else
         do L = 1, LPAR
            do J = 1, JPAR
               do I = 1, IPAR
                  GRDATA3(I,J,L) = GRDATA3HI(I,J,L)
               end do
            end do
         end do
       end if

       select case(FLD_ID)

       case (212)
          !// Mass Flux updrafts [accumulated kg/(m^2*s)]
          !// - Mass flux is through BOTTOM EDGE of grid box (i.e. model
          !//   half layers).
          !// - EC-data are stored from half layer 2, since the flux into
          !//   layer 1 is always zero (no flux in through surface).
          !// - Scale with area and dt  --> [kg/s]
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
                if (GRDATA3(I,J,LL) .gt. mincwetelim) then
                  CWETE(I,J,L) = GRDATA3(I,J,LL) * AREAXY(I,J) * ZDT
                else
                  CWETE(I,J,L) = 0._r8
                end if
              end do
            end do
          end do


       case (213)
          !// Mass Flux downdrafts [accumulated kg/(m^2*s)]
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
          CWETD(:,:,:) = 0._r8 !// No flux down to surface
          !// Start retrieval at layer 2
          do L = 2, LWDPAR
            !// See Case(212) for comment on LMAP
            LL = LMAP(L) - 1
            do J = 1, JPAR
              do I = 1, IPAR
                if (GRDATA3(I,J,LL) .lt. maxcwetdlim) then
                  CWETD(I,J,L) = GRDATA3(I,J,LL) * AREAXY(I,J) * ZDT
                else
                  CWETD(I,J,L) = 0._r8
                end if
              end do
            end do
          end do


       case (214)
          !// Updrafts detrainment rate [accumulated kg/(m3*s); kg/(m2*s)
          !// per gridbox height]
          !// - Entrainment has to be built from detrainment.
          !// - The rate is per height, so we have to multiply with dZ and area.
          !// - Detrainment/entrainment rates are given at GRID CENTER
          !//   (i.e. model full layers).
          !// - Detrainment must be summed up when collapsing layers!
          !// - Scale with box height, area and dt --> [kg/s]
          CDETU(:,:,:) = 0._r8
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, LWEPAR
                do LL = LMAP(L), LMAP(L+1) - 1
                  !// Detrainment must be summed up when collapsing layers
                  if (GRDATA3(I,J,LL) .gt. mincdetulim) then
                    CDETU(I,J,L) = CDETU(I,J,L) &
                         + GRDATA3(I,J,LL) * XLMMAP(LL) &
                           * (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) &
                           * AREAXY(I,J) * ZDT
                  end if
                end do
              end do
            end do
          end Do


       case (215)
          !// Downdrafts detrainment rate [acc. kg/(m3*s); kg/(m2*s) per
          !// gridbox height]
          !// - Entrainment has to be built from detrainment.
          !// - The rate is per height, so we have to multiply with dZ and
          !//   area.
          !// - Detrainment/entrainment rates are given at GRID CENTER
          !//   (i.e. model full layers).
          !// - Detrainment must be summed up when collapsing layers!
          !// - Scale with box height and dt --> [kg/s]
          CDETD(:,:,:) = 0._r8
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, LWDPAR
                do LL = LMAP(L),LMAP(L+1)-1
                  if (GRDATA3(I,J,LL) .gt. mincdetdlim) then
                    CDETD(I,J,L) = CDETD(I,J,L) &
                         + GRDATA3(I,J,LL) * XLMMAP(LL) &
                           * (ZOFLE(L+1,I,J) - ZOFLE(L,I,J)) &
                           * AREAXY(I,J) * ZDT
                  end if
                end do
              end do
            end do
          end do


       case (216)
          !/ Rainfall [kg/(m^2)] (accumulated kg/(m^2*s)) - edge value
          !// scale with area and dt --> [kg/s]
          do L = 1, LWEPAR
            LL = LMAP(L)
            do J = 1, JPAR
              do I = 1, IPAR
                PRECLS(I,J,L) = max(0._r8, GRDATA3(I,J,LL) * AREAXY(I,J) * ZDT)
              end do
            end do
          end do
          found216 = .TRUE.


       case (217)
          !// CONVECTIVE RAINFALL [kg/(m^2)] (accumulated kg/(m^2*s))
          !// EDGE VALUE
          !// Scale with area and dt or re-initialize --> [kg/s]
          TMPCRAIN(:,:,:) = 0._r8
          do L = 1, LWEPAR
            LL = LMAP(L)
            do J = 1, JPAR
              do I = 1, IPAR
                TMPCRAIN(I,J,L) = max(0._r8, GRDATA3(I,J,LL) * AREAXY(I,J) * ZDT)
              end do
            end do
          end do
          PRECCNV(:,:,:) = TMPCRAIN(:,:,:)

       case (218)
          !// LARGE SCALE RAINFALL 2 [kg/(m^2)] (accumulated kg/(m^2*s))
          !// EDGE VALUE
          !// Scale with area and dt or re-initialize --> [kg/s]
          TMPLRAIN(:,:,:) = 0._r8
          do L = 1, LWEPAR
            LL = LMAP(L)
            do J = 1, JPAR
              do I = 1, IPAR
                TMPLRAIN(I,J,L) = max(0._r8, GRDATA3(I,J,LL) * AREAXY(I,J) * ZDT)
              end do
            end do
          end do

       case (246)
          !// Cloud Liquid Water Content [Kg/Kg] - midpoint value
          !// Cloud routine needs CLDLWCW in any case, so we calculate
          !// that first. Could repeat with GRDATA3, but rather do
          !// degradation of CLDLWCW.
          CLDLWCW(:,:,:)  = 0._r8
          do L = 1, LWEPAR
            do LL = LMAP(L), LMAP(L+1)-1
              do J = 1, JPARW
                do I = 1, IPARW
                  CLDLWCW(I,J,L) = CLDLWCW(I,J,L) &
                                   + GRDATA3HI(I,J,LL) * XLMMAP(LL)
                end do
              end do
            end do
          end do
          if (ldeg) then
            call TRUNG8(CLDLWCW, CLDLWC, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                        JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)
          else
            CLDLWC(:,:,:) = 0._r8
            do L = 1, LWEPAR
              do J = 1,JPAR
                do I = 1,IPAR
                   CLDLWC(I,J,L) = CLDLWCW(I,J,L)
                end do
              end do
            end do
          end if

       case (247)
          !// Cloud Ice Water Content [Kg/Kg] - midpoint value
          !// Cloud routine needs CLDIWCW in any case, so we calculate
          !// that first.
          CLDIWCW(:,:,:)  = 0._r8
          do L = 1, LWEPAR
            do LL = LMAP(L), LMAP(L+1)-1
              do J = 1, JPARW
                do I = 1, IPARW
                  CLDIWCW(I,J,L) = CLDIWCW(I,J,L) &
                                   + GRDATA3HI(I,J,LL) * XLMMAP(LL)
                end do
              end do
            end do
          end do

          if (ldeg) then
            call TRUNG8(CLDIWCW, CLDIWC, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                        JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)
          else
            CLDIWC(:,:,:) = 0._r8
            do L = 1, LWEPAR
              do J = 1, JPAR
                do I = 1, IPAR
                  CLDIWC(I,J,L) = CLDIWCW(I,J,L)
                end do
              end do
            end do
          end if

       case (248)
          !// Cloud Fraction [0,1] - midpoint value
          !// Cloud routine needs CLDFRW in any case, so we calculate
          !// that first.
          CLDFRW(:,:,:) = 0._r8
          do L = 1, LWEPAR
            do LL = LMAP(L), LMAP(L+1)-1
              do J = 1, JPARW
                do I = 1, IPARW
                  CLDFRW(I,J,L) = CLDFRW(I,J,L) &
                       + XLMMAP(LL) * max(min(GRDATA3HI(I,J,LL),1._r4),0._r4)
                end do
              end do
            end do
          end do

          if (ldeg) then
            call TRUNG8(CLDFRW, CLDFR, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                        JDGRD, IPARW, JPARW, IPAR, JPAR, LWEPAR, 1)
          else
            do L = 1, LWEPAR
              do J = 1, JPAR
                do I = 1, IPAR
                  CLDFR(I,J,L) = CLDFRW(I,J,L)
                end do
              end do
            end do
          end if

          !// Cloud cover routines
          if (LOSLOCHEM) &
             !// New cloud treatment (qcode_60a)
             call CLOUD(CLDFRW,CLDIWCW,CLDLWCW,PW,TW,ETAA,ETAB,AREAXYW, &
                        ZOFLEW,ZDEGI,ZDEGJ,IMAP,JMAP, &
                        SWSTORE,CLDSTORE,TYPSTORE,RAN4,IRAN0, &
                        LCLDAVG,LCLDQMD,LCLDQMN,LCLDRANA,LCLDRANQ)

       case DEFAULT
          write(6,'(a,i7)') f90file//':'//subr//': UNKNOWN 3D FIELD',FLD_ID
          STOP
       end select

    end do !// Do NF = 1, NrOf3Dfields

    !// POTENTIAL VORTICITY
    !// PV was not available for some older metdata. If not available,
    !// it will be calculated.
    read(14,iostat=ioerr) FLD_ID
    GRDATA3(:,:,:)   = 0._r4
    GRDATA3HI(:,:,:) = 0._r4
    LPVU = .true. !// Assume PVU must be calculated
    if (ioerr .eq. 0) then
       if (FLD_ID .eq. 60) then
          read(14) GRDATA3HI
          LPVU = .false. !// No need to calculate PVU
       else
          backspace(14) !// No PV, take one step back.
       end if
    else
       !// Last time step and PVU was not on file
    end if

    !// Save UMS and VMS
    !// Degrade vertically and horizontally, put into (LPAR,IPAR,JPAR)
    call ijlw2lij(UMSW, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
         JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, UMS)
    call ijlw2lij(VMSW, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
         JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, VMS)


    if (LPVU) then
       !// Generate PV on model resolution and convert to PVU
       call get_pvu()
    else
       !// Degrade into PVU
       WK3DU(:,:,:) =real(GRDATA3HI, r8) * 1.e6_r8
       call ijlw2lij(WK3DU, LMAP, XLMMAP, IPAR, JPAR, LPAR, IPARW, &
               JPARW, LPARW, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, JDGRD, PVU)
    end if


    !//---------------------------------------------------------------------
    !// 2-d GRID POINT DATA
    !//---------------------------------------------------------------------
    NrOf2Dfields = 74
    if (LPARW .eq. 40) then
       if (MYEAR .lt. 2001) then
          NrOf2Dfields = 68
          ! Missing fields in these data
          !  57 Downward UV radiation at the surface
          !  58 Photosynthetically active radiation at the surface
          !  59 Convective available potential energy
          !  78 Total column liquid water
          !  79 Total column ice water
          ! 234 Log Surface roughness for Heat
       end if
    end if


    !// for ease check on FLD_ID
    do NF = 1, NrOf2Dfields

       read(15) FLD_ID
       if (ldeg) then
          read(15) GRDATA2HI
          call TRUNG4(GRDATA2HI, GRDATA2, ZDEGI, ZDEGJ, IMAP, JMAP, IDGRD, &
                      JDGRD, IPARW, JPARW, IPAR, JPAR, 1, 1)
       else
          read(15) GRDATA2
       end if

       !// use only fields we want
       select case(FLD_ID)

       case (31)
          !// SEA ICE COVER (CI)
          CI(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','Sea ice (CI)')

       case (32)
          !// Snow Albedo (ASN)
          if (verbose) call skipData(FLD_ID,'SFC','Snow Albedo (ASN)')

       case (33)
          !// Snow Density (RSN)
          if (verbose) call skipData(FLD_ID,'SFC','Snow Density (RSN)')

       case (34)
          !// Sea Surface Temperature-Kelvin (SSTK)
          if (verbose) call skipData(FLD_ID,'SFC','SSTK')

       case (35)
          !// Ice surface temperature Layer 1 (ISTL1)
          if (verbose) call skipData(FLD_ID,'SFC','ISTL1')

       case (36)
          !// Ice surface temperature Layer 2 (ISTL2) 
          if (verbose) call skipData(FLD_ID,'SFC','ISTL2')

       case (37)
          !// Ice surface temperature Layer 3 (ISTL3) 
          if (verbose) call skipData(FLD_ID,'SFC','ISTL3')

       case (38)
          !// Ice surface temperature Layer 4 (ISTL4) 
          if (verbose) call skipData(FLD_ID,'SFC','ISTL4')

       case (39)
          !// Volumetric soil water layer 1 (SWVL1)         
          if (LDUST) then
             !// Limit to positive values, just in case.
             R8XY(:,:) = max(real(GRDATA2(:,:), r8), 0._r8)
             call dust_set_SWVL1(R8XY)
             if (verbose) call gotData(FLD_ID,'SFC','Volumetric soil water layer 1 (SWVL1)')
          else
             if (verbose) call skipData(FLD_ID,'SFC','Volumetric soil water layer 1 (SWVL1)')
          end if
       case (40)
          !// Volumetric soil water layer 2 (SWVL2)
          if (verbose) call skipData(FLD_ID,'SFC','SWVL2')

       case (41)
          !// Volumetric soil water layer 3 (SWVL3)
          if (LmOSaic) then
             if (verbose) call gotData(FLD_ID,'SFC','Volumetric soil water layer 3 (SWVL3)')
          else
             if (verbose) call skipData(FLD_ID,'SFC','SWVL3')
          end if

       case (42)
          !// Volumetric soil water layer 4 (SWVL4)
          if (verbose) call skipData(FLD_ID,'SFC','SWVL4')

       case (44)
          !// Snow Evaporation (ES) (accumulated, [m water equivalent])
          !// Snow evaporation may actually also be negative, for some reason
          GRDATA2(:,:) = GRDATA2(:,:) * ZDT
          call r4data2mpblocks(GRDATA2, ES)
          if (verbose) call gotData(FLD_ID,'SFC','Snow Evaporation (ES)')

       case (45)
          !// Snow Melt (SMLT) (accumulated, [m water equivalent])
          !// Limit to positive values
          GRDATA2(:,:) = max(0._r8, GRDATA2(:,:)) * ZDT
          call r4data2mpblocks(GRDATA2, SMLT)
          if (verbose) call gotData(FLD_ID,'SFC','Snow Melt (SMLT)')

       case (49)
          !// 10m wind gust
          if (verbose) call skipData(FLD_ID,'SFC','10m wind gust')

       case (50)
          !// Large-scale precipitation fraction (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Large-scale prec. frac.')

       case (57)
          !// Downward UV radiation at the surface
          if (verbose) call skipData(FLD_ID,'SFC','Downw. UV rad. sfc')

       case (58)
          !// Photosynthetically active radiation (PhotActRad) at the surface
          !// (accumulated)
          !// Unit: (W/m2)*s, accumulated. Divide by ZDT to get W/m2.
          PhotActRad(:,:) = max(0._r8, GRDATA2(:,:)) * ZDT
          if (verbose) call gotData(FLD_ID,'SFC','Photosyn. act. rad. sfc (PAR)')

       case (59)
          !// Convective available potential energy
          if (verbose) call skipData(FLD_ID,'SFC','CAPE')

       case (78)
          !// Total column liquid water
          if (verbose) call skipData(FLD_ID,'SFC','Tot. col. liq. water')

       case (79)
          !// Total column ice water
          if (verbose) call skipData(FLD_ID,'SFC','Tot. col. ice water')

       case (129)
          !// Surface Geopotential (Z)
          if (verbose) call skipData(FLD_ID,'SFC','Sfc. geopot. (Z)')

       case (136)
          !// Total Column Water (TCW)
          if (verbose) call skipData(FLD_ID,'SFC','Tot. col. water')

       case (137)
          !// Total Column Water Vapor (TCWV)
          if (verbose) call skipData(FLD_ID,'SFC','Tot. col. water vap.')

       case (139)
          !// Surface Temperature (ST)
          if (verbose) call skipData(FLD_ID,'SFC','Sfc. temp.')

       case (141)
          !// SNOW DEPTH (SD)
          SD(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','Snow depth (SD)')

       case (142)
          !// Large Scale Precipitation (stratiform) (accumulated)
          LSPREC(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','Large Scale Precip. (LSPREC)')

       case (143)
          !// Convective Precipitation (accumulated)
          CNVPREC(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','Conv. Precip. (CONVPREC)')

       case (144)
          !// SnowFall (SF) (accumulated, [m water equivalent])
          GRDATA2(:,:) = GRDATA2(:,:) * ZDT
          call r4data2mpblocks(GRDATA2,SNFL)
          if (verbose) call gotData(FLD_ID,'SFC','Snow fall (SNFL)')

       case (145)
          !// Boundary Layer Dissipation (BLD) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Boundary Layer Dissipation')


       case (146)
          !// Surface Sensible Heat Flux SSHF
          SHF(:,:) = -GRDATA2(:,:) * ZDT
          !// Fluxes of zero cause grief with deposition, PBL.
          !// SHF numbers from -400 to 100, set a small number to avoid
          !// overflow, not worry about small negative going to small positive
          do J = 1, JPAR
            do I = 1, IPAR
              if (abs(SHF(I,J)) .lt. 1.e-6_r8) SHF(I,J) = 1.e-6_r8
            end do
          end do
          if (verbose) call gotData(FLD_ID,'SFC','Sfc. sens. heat flux (SHF)')

       case (147)
          !// Surface Latent Heat Flux SLHF
          SLH(:,:) = -GRDATA2(:,:) * ZDT
          !// Fluxes of zero cause grief with deposition, PBL.
          !// SLH: set small value if zero.
          do J = 1, JPAR
            do I = 1, IPAR
              if (SLH(I,J) .eq. 0._r8) SLH(I,J) = 1.e-30_r8
            end do
          end do
          if (verbose) call gotData(FLD_ID,'SFC','Sfc. lat. heat flux (SLH)')

       case (148)
          !// Surface stress (eller Charnock)
          if (verbose) call skipData(FLD_ID,'SFC','Sfc stress/Charnock')

       case (151)
          !// Mean Sealevel Pressure (Diagnostic only - not in old 19-layer)
          MSLP(:,:) = GRDATA2(:,:) * 1.e-2_r8
          if (verbose) call gotData(FLD_ID,'SFC','MSL pres.')

       case (159)
          !// Boundary Layer Height BLH
          BLH(:,:) = GRDATA2(:,:)
          !// Average BLH at poles
          do J = 1, JPAR, JPAR-1
            DMASS = 0._r8
            do I = 1, IPAR
              DMASS = DMASS + GRDATA2(I,J)
            end do
            DMASS = DMASS / real(IPAR, r8)
            do I = 1, IPAR
              BLH(I,J) = DMASS
            end do
          end do
          if (verbose) call gotData(FLD_ID,'SFC','Boundary Layer Height (BLH)')

       case (164)
          !// Total Cloud Cover
          if (verbose) call skipData(FLD_ID,'SFC','Tot. cloud cover')

       case (165)
          !// 10m U wind component U10M
          SFU(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','10m U-wind (SFU)')

       case (166)
          !// 10m V wind component V10M
          SFV(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','10m V-wind (SFV)')

       case (167)
          !// 2m Temperature T2M
          SFT(:,:) = GRDATA2(:,:)
          if (minval(SFT) .eq. 0._r8) then
             write(6,'(a)')  f90file//':'//subr//': SFT zero: VERY WRONG!'
             stop
          end if
          if (verbose) call gotData(FLD_ID,'SFC','2m Temperature (SFT)')

       case (168)
          !// 2m Dewpoint Temperature TDEW2M
          !// Read in Td, convert to mixing ratio using Clausius-Clapeyron
          !// SFQ is surface specific humidity.
          do J = 1, JPAR
            do I = 1, IPAR
              SFTD = GRDATA2(I,J)
              !// Stull (1988):
              LV = (2.501_r8 - (2.37e-3_r8 * (SFT(I,J) - 273.16_r8)))*1.e6_r8
              ESAT = 610.78_r8 * EXP( (-LV/R_H2O) &
                                     * ((1._r8/SFTD) - (1._r8/273.16_r8)) )
              !// The equation above is not valid below 0C, and should
              !// rather be calculated from another equation.
              !// The vapor pressure at T is equal to saturation pressure
              !// at Td.
              SFQ(I,J) = (R_AIR/R_H2O) * (ESAT/(P(I,J) * 100._r8))
            end do
          end do
          if (verbose) call gotData(FLD_ID,'SFC','2m dew point -> SFQ')

       case (169)
          !// Surface Solar Radiation Downwards (SSRD) (accumulated)
          !// Unit: W/m2, accumulated (W/m2*s)
          if (LDUST) then
             R8XY(:,:) = real(GRDATA2(:,:), r8) * ZDT
             R8XY(:,:) = max(R8XY(:,:), 0._r8)
             call dust_set_ssrd(R8XY)
             if (verbose) call gotData(FLD_ID,'SFC','Surface Solar Radiation Downwards (SSRD)')
          else
             !// Surface Solar Radiation Downwards (SSRD) (accumulated)
             if (verbose) call skipData(FLD_ID,'SFC','Surface Solar Radiation Downwards (SSRD)')
          end if


       case (170)
          !// Soil Temperature Level 2 (STL2)
          if (verbose) call skipData(FLD_ID,'SFC','Soil temp. lev. 2')

       case (172)
          !// Land/Sea mask (LSM), NOTE: {0,1}
          !// Not so useful; it is 1 if land fraction > 0.5.
          if (verbose) call skipData(FLD_ID,'SFC','Land/sea mask')

       case (175)
       !// Surface Thermal Radiation Downwards (STRD) (accumulated)
          !// Unit: W/m2, accumulated (W/m2*s)
          if (LDUST) then
             !// Limit to positive values just in case
             R8XY(:,:) = max(real(GRDATA2(:,:), r8), 0._r8) * ZDT
             call dust_set_strd(R8XY)
             if (verbose) call gotData(FLD_ID,'SFC','Surface Thermal Radiation Downwards (STRD)')
          else
             if (verbose) call skipData(FLD_ID,'SFC','Surface Thermal Radiation Downwards (STRD)')
          end if

       case (176)
          !// Surface Solar Radiation (SSR) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Sfc. solar rad.')

       case (177)
          !// Surface Thermal Radiation (STR) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Sfc. thermal rad.')

       case (178)
          !// Top Solar Radiation (TSR) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Top Solar Radiation (TSR)')

       case (179)
          !// Top Thermal Radiation (TTR) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Top Thermal Radiation (TTR)')

       case (180)
          !// E/W Surface Stress EWSS (read in as Ns/m2; divide by
          !// timestep to get mean)
          EWSS(:,:) = GRDATA2(:,:) * ZDT
          if (verbose) call gotData(FLD_ID,'SFC','E/W sfc stress (EWSS)')

       case (181)
          !// N/S Surface Stress NSSS (read in as Ns/m2; divide by
          !// timestep to get mean)
          NSSS(:,:) = GRDATA2(:,:) * ZDT
          if (verbose) call gotData(FLD_ID,'SFC','N/S sfc stress (NSSS)')

       case (182)
          !// Evaporation (E) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Evaporation')

       case (183)
          !// Soil Temperature Level 3 (STL3)
          if (verbose) call skipData(FLD_ID,'SFC','Soil temp. lev. 3')

       case (184)
          !// Soil Wetness Level 3 (SWL3)
          if (verbose) call skipData(FLD_ID,'SFC','Soil wetness level 3')

       case (185)
          !// Convective Cloud Cover (CCC)
          if (verbose) call skipData(FLD_ID,'SFC','Convective cloud cover')

       case (186)
          !// Low Cloud Cover (LCC)
          if (verbose) call skipData(FLD_ID,'SFC','Low cloud cover')

       case (187)
          !// Medium Cloud Cover (MCC)
          if (verbose) call skipData(FLD_ID,'SFC','Medium cloud cover')

       case (188)
          !// High Cloud Cover (HCC)
          if (verbose) call skipData(FLD_ID,'SFC','High cloud cover')

       case (189)
          !// Sun Shine Duration (SUND) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','Sunshine duration')

       case (195)
          !// Latitudinal Gravity Wave Stress (LGWS) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','LGWS')

       case (196)
          !// Meridional Gravity Wave Stress (MGWS)
          if (verbose) call skipData(FLD_ID,'SFC','MGWS')

       case (197)
          !// Gravity Wave Dissipation (GWD)
          if (verbose) call skipData(FLD_ID,'SFC','GWD')

       case (198)
          !// Skin Reservoir Content (SRC)
          if (verbose) call skipData(FLD_ID,'SFC','Skin Reservoir Content (SRC)')

       case (201)
          !// Maximum Temperature at 2m (MX2T)
          if (verbose) call skipData(FLD_ID,'SFC','MX2T')

       case (202)
          !// Minimum Temperature at 2m (MN2T)
          if (verbose) call skipData(FLD_ID,'SFC','MN2T')

       case (205)
          !// RunOff (RO)
          if (verbose) call skipData(FLD_ID,'SFC','RunOff')

       case (208)
          !// Top Net Solar Radiation Clear Sky (TSRC) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','TSRC')

       case (209)
          !// Top Net Thermal Radiation Clear Sky (TTRC) (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','TTRC')

       case (210)
          !// Surface net solar radiation, clear sky (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','SNSRCS')

       case (211)
          !// Surface net thermal radiation, clear sky (accumulated)
          if (verbose) call skipData(FLD_ID,'SFC','SNTRCS')

       case (234)
          !// Log Surface roughness for Heat
          if (verbose) call skipData(FLD_ID,'SFC','Log(sfc. roughness for heat)')

       case (235)
          !// Skin Temperature (SKT)
          if (verbose) call skipData(FLD_ID,'SFC','Skin Temperature (SKT)')

       case (236)
          !// Soil Temperature Level 4 (STL4)
          if (verbose) call skipData(FLD_ID,'SFC','Soil temp. lev. 4')

       case (238)
          !// Temperature of Snow Layer (TSN)
          if (verbose) call skipData(FLD_ID,'SFC','Temperature of snow layer')

       case (243)
          !// Forecast Albedo (Not in old 19-layer)
          SA(:,:) = GRDATA2(:,:)
          if (verbose) call gotData(FLD_ID,'SFC','Forecast Albedo (SA)')

       case (244)
          !// Forecast Surface Roughness (FSR) 
          if (verbose) call skipData(FLD_ID,'SFC','Forecast sfc. roughness')

       case (245)
          !// Forecast Log(Surface Roughness Heat)
          if (verbose) call skipData(FLD_ID,'SFC','Forecast log(Sfc. roughness heat)')

       case DEFAULT
          if (verbose) call skipData(FLD_ID,'SFC','Default')

       end select

    end do !// do NF = 1, NrOf2Dfields


    !// set correct values of PRECLS/PRECCNV
    if (found216) then
       !// Old RAIN = LS + CNV. Build CNV and LS from 2D LSPREC:
       do L = 1, LWEPAR
         do J = 1, JPAR
           do I = 1, IPAR
             if ((CNVPREC(I,J)+LSPREC(I,J)).gt.0._r8) then
               CNVFRAC = min(1._r8, CNVPREC(I,J) / (CNVPREC(I,J) + LSPREC(I,J)))
             else
               CNVFRAC = 0._r8
             end if
             RAINTOT = PRECLS(I,J,L)
             PRECCNV(I,J,L) = RAINTOT * CNVFRAC
             PRECLS(I,J,L)  = RAINTOT * (1._r8 - CNVFRAC)
           end do
         end do
       end do
    else
       !// PRECLS is large scale rain, not total rain
       do L = 1, LWEPAR
         do J = 1, JPAR
           do I = 1, IPAR
             PRECLS(I,J,L) = max(0._r8, TMPLRAIN(I,J,L))
           end do
         end do
       end do
    end if



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


    !// Derive for 19-Layer  (Assume Land 0.15, Ocean 0.05, Ice 0.65)
    if (LPARW .eq. 19) then
      do J = 1, JPAR
        do I = 1, IPAR
          SA(I,J) = 0.05_r8 + 0.10_r8 * PLAND(I,J)   ! Land and Ocean
          if(SFT(I,J) .lt. 263._r8) SA(I,J) = 0.65_r8   ! Snow/Ice Cover
          if(YDGRD(J) .lt. -65._r8 .or. YDGRD(J) .gt. 75._r8) &
               SA(I,J) = 0.65_r8                    ! Polar or Sea Ice
        end do
      end do
    end if

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
    deallocate( FFTAR, WORK, SPWK5, SPWK6, PSPEHI, &
         VTMP, UTMP, WK3DU, WK3DV, WK3DT, ZOFLEW, QW, TW, &
         CLDFRW, CLDIWCW, CLDLWCW, GRDATA3HI, UMSW, VMSW, &
         TMPCRAIN, TMPLRAIN, CDETU, CDETD, GRDATA3, &
         PW, GRDATA2HI, LSPREC, CNVPREC, EWSS, NSSS, GRDATA2, R8XY )

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
  subroutine r4data2mpblocks(rdata, mpdata)
    !//---------------------------------------------------------------------
    !// Puts real*4 data of size (ipar,jpar) into real*8 MP-block structure
    !// mpdata(idblk,jdblk,mpblk).
    !//
    !// Ole Amund Sovde, April 2013
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, r4
    use cmn_size, only: IPAR, JPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r4), intent(in) :: rdata(ipar,jpar)
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
          mpdata(II,JJ,MP) = rdata(I,J)
        end do
      end do
    end do
    !//---------------------------------------------------------------------
  end subroutine r4data2mpblocks
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine gotData(ID,TYP,LABEL)
    !//---------------------------------------------------------------------
    !// Print info about 2D field read.
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: ID
    character(len=3) :: TYP
    character(len=*) :: LABEL
    !//---------------------------------------------------------------------
    write(6,'(a,i5,1x,a)') ' update_metdata: Read '//TYP//':   ', &
         ID, trim(LABEL)
    !//---------------------------------------------------------------------
  end subroutine gotData
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
  subroutine skipData(ID,TYP,LABEL)
    !//---------------------------------------------------------------------
    !// Print info about 2D field skipped.
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) :: ID
    character(len=3) :: TYP
    character(len=*) :: LABEL
    !//---------------------------------------------------------------------
    write(6,'(a,i5,1x,a)') ' update_metdata: Skipped '//TYP//':', &
         ID,trim(LABEL)
    !//---------------------------------------------------------------------
  end subroutine skipData
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
end module metdata_ecmwf
!//=========================================================================

!//=========================================================================
!// Summary of EC fields
!//  *.b01  Spectral Data read on unit 13 - one record per met record
!//
!//          152     Surface Pressure        (nrnm)
!//          130     Temperature             (nrnm,lpar)
!//          138     Vorticity               (nrnm,lpar)
!//          155     Divergence              (nrnm,lpar)
!//          135     Vertical Velocity       (nrnm,lpar)   [not used]
!//
!//  *.b02  Gridpoint Data read on unit 14 - one record per day (recl varies!)
!//
!//          133      Specific humidity (Q, kg/kg)       (ipar,jpar,lpar)
!//          212      Mass Flux updrafts [Kg/(m^2*s)]    (ipar,jpar,*)
!//          213      Mass Flux downdrafts               (ipar,jpar,*)
!//          214      Mass Flux updrafts entrainment     (ipar,jpar,*)
!//          215      Mass Flux downdrafts entrainment   (ipar,jpar,*)
!//          216      Rainfall                           (ipar,jpar,*)
!//          217      Large scale rainfall               (ipar,jpar,*)
!//          218      Convective rainfall                (ipar,jpar,*)
!//          246      Cloud Liquid Water Content [Kg/Kg] (ipar,jpar,*)
!//          247      Cloud Ice Water Content [Kg/Kg]    (ipar,jpar,*)
!//          248      Cloud Fraction [0,1]               (ipar,jpar,*)
!//
!//  *.b03  Gridpoint Data read on unit 15 - one record per variable
!//
!//          TCW      Total Column Water             (ipar,jpar)   [not used]
!//          LSP      Large Scale Precip.            (ipar,jpar)   [not used]
!//          CP       Convective Precip.             (ipar,jpar)   [not used]
!//          BLD      Boundary Layer Dissipation     (ipar,jpar)
!//          SSHF     Surface Sensible Heat Flux     (ipar,jpar)
!//          SLHF     Surface Latent Heat FLux SLHF  (ipar,jpar)
!//          U10M     10m U wind component           (ipar,jpar)
!//          V10M     10m V wind component           (ipar,jpar)
!//          T2M      2m Temperature                 (ipar,jpar)
!//          TDEW2M   2m Dewpoint Temperature        (ipar,jpar)
!//          EWSS     E/W Surface Stress             (ipar,jpar)
!//          NSSS     N/S Surface Stress             (ipar,jpar)
!//          RUNOFF   Runoff                         (ipar,jpar)
!//          CLDTYP   Accumulated cloud type         (ipar,jpar)   [not used]
!//          BLH      Boundary Layer Height          (ipar,jpar)
!//=========================================================================


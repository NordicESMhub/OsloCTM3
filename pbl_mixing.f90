!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Amund Sovde Haslerud, November 2017
!//=========================================================================
!// Planetary boundary layer mixing.
!//=========================================================================
module pbl_mixing
  !// ----------------------------------------------------------------------
  !// MODULE: pbl_mixing
  !// DECRIPTION: Routine for PBL mixing.
  !//
  !// This is the UCI p-pbl.f rewritten to f90, and including Holtslag
  !// kprof from an earlier UCI version.
  !// It uses e.g. KC(2), PRANDTL(1) and MO_LEN
  !//
  !// Contains
  !//   subroutine CNVDBL
  !//   subroutine BULK
  !//   subroutine get_keddy_L1
  !//   subroutine KPROF2
  !//   subroutine KPROF
  !//   real(r8) function PHIM
  !//   real(r8) function PHIH
  !//   subroutine TRIDIAG
  !//   subroutine INTERP
  !//   subroutine QXZON
  !//   subroutine QZONX
  !// ----------------------------------------------------------------------
  use cmn_size, only: LPAR
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  integer, parameter :: NSUB = 3, BLPAR = LPAR*NSUB
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'pbl_mixing.f90'
  !// ----------------------------------------------------------------------
  private
  public CNVDBL
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine CNVDBL(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,AIRB,DTPBL,MP)
    !// --------------------------------------------------------------------
    !// new subs for planetary boundary layer mixing of tracer:
    !// step = DTPBL(sec)
    !// UCI method currently only uses a simple e-fold to well mixed profile.
    !// CTM3 uses Holtslag.
    !//
    !// Amund Sovde Haslerud, November 2017
    !//   Updated to make NBLX1 possible with Oslo chemistry.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: NPAR, IDBLK, JDBLK
    use cmn_ctm, only: NBLX, MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         ETAA, ETAB, NTM, DISTX, DISTY, JMON
    use cmn_met, only: P, SHF, SFT, SFQ, SFU, SFV, SLH, SMF, USTR, &
         BLH, PBL_KEDDY, MO_LENGTH, PRANDTLL1, ZOFLE, U, V, T, Q
    use cmn_parameters, only: R_AIR
    use cmn_sfc, only: ZOI
    use utilities, only: moninobukhov_length
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer,intent(in) :: MP
    real(r8), intent(in) :: DTPBL
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) :: AIRB
    !// Input/Output
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ

    !// Parameters
    real(r8), parameter ::  DT3HR = 10800._r8     ! 3-hr e-fold

    !// Locals
    integer :: NBL, I,II,J,JJ,L,N
    real(r8), dimension(LPAR+1) :: POFL, DZL
    real(r8), dimension(LPAR) :: &
         QM, QTT, QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY, QM2
    real(r8) :: ZBL,ZCL,ZLL,ZF,ZPP,ZMF,XBL, UMIX, FMIX
    real(r8) :: SFCP, SFCS, SFCT, SFCQ, SFCU, SFCV, SFCL, SFCM, SFCD
    real(r8) :: SFMU, SFNU, USTAR, MO_LEN, PRTL1
    !real(r8) :: ZO, ZOW

    !// Additions for KPROF
    real(r8), parameter ::  VK = 0.4_r8
    real(r8), dimension(LPAR) :: DP, UVEL, VVEL, TOFL,QVAP,QTMPT,QTMPM
    real(r8) :: SFCTH, CNST, LV, TSV, MF, DZM
    integer :: NBLP1, NGMAX, NG, MID, M, K
    real(r8), dimension(NSUB) :: &              ! used for QZONX/QXZON
         PM, PTT, PXT, PXX, PYT, PYY, PXY, PZT, PZZ, PXZ, PYZ
    real(r8), dimension(BLPAR) :: &
         GM, GTT, GXT, GXX, GYT, GYY, GXY, &
         ZMD, ZD, DUDZ, DVDZ, DQDZ, DTHDZ, &
         KC, KCZ, &
         GOUT, &     ! new tracer
         BU, BD, BL  ! matrix coefficients
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'CNVDBL'
    !// --------------------------------------------------------------------

    !// Initialise prandtL1
    PrandtlL1(:,:,MP) = 0.85_r8
    PBL_KEDDY(:,:,MP) = 1.e-19_r8

    if (NBLX .eq. 1) then
       !// NBLX=1 - Prather_s bulk scheme (on mean grid)
       !//          calculate e-fold to bulk-mixed in 3 hour
       !FMIX = 1._r8 - exp(-DTPBL / DT3HR)
       UMIX = exp(-DTPBL / DT3HR)
       FMIX = 1._r8 - UMIX

       do J = MPBLKJB(MP), MPBLKJE(MP)
         JJ   = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II   = I - MPBLKIB(MP) + 1


           !// Surface properties, all in SI units from met fields.
           SFCP  = 100._r8 * P(I,J)   ! P(mbar) to SFCP (pascals)
           SFCS  = SHF(I,J)           ! sensible heat flux (W/m2)
           SFCT  = SFT(I,J)           ! surface temperature (K)
           SFCQ  = SFQ(I,J)           ! surface (2-m?) Q, u, v 
           !SFCU  = SFU(I,J)
           !SFCV  = SFV(I,J)
           SFCL  = SLH(I,J)           ! latent heat flux (W/m2)
           !SFCM  = SMF(I,J)           ! momentum flux = (u*)^2 * density
           !// surface momentum flux (SMF) = density * ustar^2    ! m/s
           USTAR = USTR(I,J)
           if (USTAR.le.0._r8) USTAR = 5.e-3_r8 ! USTAR should not be zero
           SFCD  = SFCP/(SFCT * R_AIR) ! density                   ! kg/m^3
           SFMU  = 6.2e-8_r8*SFCT ! absol.visc. 6.2d-8*T (lin fit:-30C to +40C)
           SFNU  = SFMU / SFCD      ! kinematic visc(nu) = mu/density  (m*m/s)

           MO_LENGTH(II,JJ,MP) = moninobukhov_length(SFCD,SFCT,USTAR,SFCS) !  M-O length

           !ZO = ZOI(I,J,JMON)       ! surface roughness (m)
           !!// correct water surf roughness for wind/waves:
           !ZOW = min(0.135_r8*SFNU/USTAR + 1.83e-3_r8*USTAR**2, 2.e-3_r8)

           !// For keddy
           !// surface momentum flux (SMF) = density * ustar^2    ! m/s
           SFCTH = SFCT * ((1.e5_r8 / SFCP)**0.286_r8)
           LV    = (2.501_r8 &
                - (2.37e-3_r8 * (SFCT - 273.16_r8)) ) * 1.e6_r8
           TSV   = SFCTH * (1._r8 + 0.61_r8 * SFCQ)



           !// determine number (fraction) of CTM levels to mix (XBL)
           !// - do local heights
           ZBL    = BLH(I,J)
           if (ZBL .gt. 0._r8)  then
              do L = 1, LPAR+1
                 POFL(L)  = ETAA(L) + ETAB(L)*P(I,J)
              end do
              do L = 1, LPAR
                 DZL(L) = ZOFLE(L+1,I,J) - ZOFLE(L,I,J)
              end do
              ZCL  = min(ZOFLE(LPAR+1,I,J)-ZOFLE(1,I,J),ZBL)
              ZLL = 0._r8
              NBL = 1
              do L = 1,LPAR
                 if (ZCL .ge. 0._r8) then
                    ZLL = ZCL
                    NBL = L
                 end if
                 ZCL = ZCL - DZL(L)
              end do
              !// altitude fraction of layer = NBL within the mixed PBL,
              !// find mass fraction 
              ZF   = ZLL / DZL(NBL)
              ZPP = POFL(NBL+1) / POFL(NBL)
              ZMF   = (1._r8 - ZPP**ZF) / (1._r8 - ZPP)
              !// NBL = no. layers completely in PBL, XBL = NBL + fraction
              !//       above NBL
              NBL  = NBL - 1
              XBL  = real(NBL, r8) + ZMF

              !// Get PBL_KEDDY for model surface level, but use the height of 2/3
              !// as when calculating KPROF in Holtslag-scheme.
              call get_keddy_L1(PBL_KEDDY(II,JJ,MP), PrandtlL1(II,JJ,MP), USTAR, &
                          SFCS, (SFCL/LV), MO_LEN, ZBL, VK, 0.6666667_r8*DZL(1), TSV)


              QM2(:)  = AIRB(:,II,JJ)
              do N = 1, NTM
              !  do L = 1, LPAR
              !    QM(L)  = AIRB(L,II,JJ)
              !    QTT(L) = BTT(L,N,II,JJ)
              !    QZT(L) = BZT(L,N,II,JJ)
              !    QZZ(L) = BZZ(L,N,II,JJ)
              !    QXZ(L) = BXZ(L,N,II,JJ)
              !    QYZ(L) = BYZ(L,N,II,JJ)
              !    QXT(L) = BXT(L,N,II,JJ)
              !    QXX(L) = BXX(L,N,II,JJ)
              !    QYT(L) = BYT(L,N,II,JJ)
              !    QYY(L) = BYY(L,N,II,JJ)
              !    QXY(L) = BXY(L,N,II,JJ)
              !  end do
                !// ASH 2017: I think this should be faster
                QM(:)  = QM2(:)
                QTT(:) = BTT(:,N,II,JJ)
                QZT(:) = BZT(:,N,II,JJ)
                QZZ(:) = BZZ(:,N,II,JJ)
                QXZ(:) = BXZ(:,N,II,JJ)
                QYZ(:) = BYZ(:,N,II,JJ)
                QXT(:) = BXT(:,N,II,JJ)
                QXX(:) = BXX(:,N,II,JJ)
                QYT(:) = BYT(:,N,II,JJ)
                QYY(:) = BYY(:,N,II,JJ)
                QXY(:) = BXY(:,N,II,JJ)
            
                call BULK(QTT,QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY,QM,XBL,LPAR)

                !do L = 1, LPAR
                !  BTT(L,N,II,JJ) = QTT(L)*FMIX + BTT(L,N,II,JJ)*(1._r8-FMIX)
                !  BZT(L,N,II,JJ) = QZT(L)*FMIX + BZT(L,N,II,JJ)*(1._r8-FMIX)
                !  BZZ(L,N,II,JJ) = QZZ(L)*FMIX + BZZ(L,N,II,JJ)*(1._r8-FMIX)
                !  BXZ(L,N,II,JJ) = QXZ(L)*FMIX + BXZ(L,N,II,JJ)*(1._r8-FMIX)
                !  BYZ(L,N,II,JJ) = QYZ(L)*FMIX + BYZ(L,N,II,JJ)*(1._r8-FMIX)
                !  BXT(L,N,II,JJ) = QXT(L)*FMIX + BXT(L,N,II,JJ)*(1._r8-FMIX)
                !  BXX(L,N,II,JJ) = QXX(L)*FMIX + BXX(L,N,II,JJ)*(1._r8-FMIX)
                !  BYT(L,N,II,JJ) = QYT(L)*FMIX + BYT(L,N,II,JJ)*(1._r8-FMIX)
                !  BYY(L,N,II,JJ) = QYY(L)*FMIX + BYY(L,N,II,JJ)*(1._r8-FMIX)
                !  BXY(L,N,II,JJ) = QXY(L)*FMIX + BXY(L,N,II,JJ)*(1._r8-FMIX)
                !end do
                !// ASH 2017: I think this should be faster
                BTT(:,N,II,JJ) = QTT(:)*FMIX + BTT(:,N,II,JJ)*UMIX
                BZT(:,N,II,JJ) = QZT(:)*FMIX + BZT(:,N,II,JJ)*UMIX
                BZZ(:,N,II,JJ) = QZZ(:)*FMIX + BZZ(:,N,II,JJ)*UMIX
                BXZ(:,N,II,JJ) = QXZ(:)*FMIX + BXZ(:,N,II,JJ)*UMIX
                BYZ(:,N,II,JJ) = QYZ(:)*FMIX + BYZ(:,N,II,JJ)*UMIX
                BXT(:,N,II,JJ) = QXT(:)*FMIX + BXT(:,N,II,JJ)*UMIX
                BXX(:,N,II,JJ) = QXX(:)*FMIX + BXX(:,N,II,JJ)*UMIX
                BYT(:,N,II,JJ) = QYT(:)*FMIX + BYT(:,N,II,JJ)*UMIX
                BYY(:,N,II,JJ) = QYY(:)*FMIX + BYY(:,N,II,JJ)*UMIX
                BXY(:,N,II,JJ) = QXY(:)*FMIX + BXY(:,N,II,JJ)*UMIX

              end do
              !Should we set AIRB(:,II,JJ) = QM(:)
           end if     ! end of ZBL>0

         end do      ! J
       end do       ! I

    else if (NBLX .eq. 5) then
       !// Holtslag kprof
       MID  = int(NSUB/2) + 1

       do J = MPBLKJB(MP), MPBLKJE(MP)
         JJ   = J - MPBLKJB(MP) + 1
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II   = I - MPBLKIB(MP) + 1


           !// Surface properties (all in SI units from EC met)
           SFCP  = 100._r8 * P(I,J)    ! P(mbar) to SFCP (pascals)
           SFCS  = SHF(I,J)           ! sensible heat flux (W/m2)
           SFCT  = SFT(I,J)           ! surface temperature (K)
           SFCQ  = SFQ(I,J)           ! surface (2-m?) Q, u, v 
           SFCU  = SFU(I,J)
           SFCV  = SFV(I,J)
           SFCL  = SLH(I,J)           ! latent heat flux (W/m2)
           SFCM  = SMF(I,J)           ! momentum flux = (u*)^2 * density
           !// surface momentum flux (SMF) = density * ustar^2    ! m/s
           USTAR = USTR(I,J)
           if (USTAR.le.0._r8) USTAR = 5.e-3_r8 ! USTAR should not be zero
           SFCD  = SFCP / (SFCT * R_AIR) ! density (kg/m^3)
           SFMU  = 6.2e-8_r8*SFCT ! absol.visc. 6.2d-8*T (lin fit:-30C to +40C)
           SFNU  = SFMU / SFCD    ! kinematic visc(nu) = mu/density  (m*m/s)
           MO_LEN = -256._r8 * SFCD * SFCT * (USTAR**3) / SFCS !  M-O length
           if (MO_LEN .eq. 0._r8) then
              write(6,'(a,2i5,4es16.6)') f90file//':'//subr//': MO_LEN is 0! -> 0.001: '// &
                   'NBLX,I,J,SFCD,SFCT,USTAR,SFCS',NBLX,I,J,SFCD,SFCT,USTAR,SFCS
              MO_LEN = 1.e-3_r8
           end if
           if (MO_LEN .ne. MO_LEN) then
              write(6,'(a,2i5,4es16.6)') f90file//':'//subr//': MO_LEN is 0! -> 1000: '// &
                   'NBLX,I,J,SFCD,SFCT,USTAR,SFCS',NBLX,I,J,SFCD,SFCT,USTAR,SFCS
              MO_LEN = 1.e3_r8
           end if
           MO_LENGTH(II,JJ,MP) = MO_LEN !// Set MO-length
           !ZO = ZOI(I,J,JMON)       ! surface roughness (m)
           !// correct water surf roughness for wind/waves:
           !ZOW = 0.135_r8 * SFNU / USTAR + 1.83e-3_r8 * USTAR**2
           !ZOW = min(ZOW, 2.e-3_r8)

           !frac_n = DTCHM2 * (real(NSTEP, r8) - 0.5_r8) / real(DTMET, r8)
           !ZBL   = frac_c * BLH_CUR(I,J) + (1._r8 - frac_c) * BLH_NEXT(I,J)
           !// check for non-zero PBL height (< 1 m)
           ZBL   = BLH(I,J)
           if (ZBL .ge. 1._r8) then
              !// Double the BLH; originally =min(BLH(I,J)*2._r8,9000._r8)
              !// There may be problems with ZBL>8500, where NG may become too
              !// large. 
              ZBL = min(ZBL * 2._r8, 8000._r8)

              !// determine number of CTM levels to mix (XBL)
              !// Collect column variables to calculate height of layers
              !// POFL(L) is in mbar.
              do L = 1, LPAR+1
                 POFL(L)  = ETAA(L) + ETAB(L) * P(I,J)
              end do
              ZLL  = 0._r8
              do L = 1, LPAR
                 DP(L)  = POFL(L) - POFL(L+1)
                 !// R = 287.*(1-Q) + 461.5*Q; assume 0.5% bl w.v.==> R = 288.
                 !// delta-z = dln(P) * R * T / g; where R/g = 288/9.81 = 29.36
                 DZL(L) = ZOFLE(L+1,I,J) - ZOFLE(L,I,J)
                 ZLL    = ZOFLE(L+1,I,J) !// Top of layer
              end do

              !// determine number of CTM layers (1:NBL) to reach ht XBL
              !// then subgrid levels (1:NG, based on subdividing CTM layers)
              !// ZLL is top of layer L (in meters)
              ZBL  = min(ZLL, ZBL)
              CNST = ZBL
              NBL = 1
              do L = 1, LPAR
                 if (CNST .ge. 0._r8) then
                    ZLL = CNST
                    NBL = L
                 end if
                 CNST = CNST - DZL(L)
              end do
              !// altitude fraction of boundary layer in layer NBL
              !// containing PBL top
              ZF   = ZLL / DZL(NBL)
              !// convert altitude fraction ZF to mass fraction MF
              CNST = POFL(NBL+1) / POFL(NBL)
              MF   = (1._r8 - CNST**ZF) / (1._r8 - CNST)
              !// NBL is no. of layers completely in PBL, XBL includes the
              !// fraction above NBL
              NBL  = NBL - 1
              XBL  = real(NBL, r8) + MF

              !// if boundary layer scheme = 2,3, or 4, then need to
              !// subdivide and get data.

              !// collect column variables for use in creating PBL profiles
              do L = 1, LPAR
                 UVEL(L) = U(I,J,L) / (DP(L) * DISTY(J))
                 VVEL(L) = V(I,J,L) / (DP(L) * DISTX(J))
                 TOFL(L) = T(I,J,L)
                 QVAP(L) = Q(I,J,L)
              end do
              !// surface momentum flux (SMF) = density * ustar^2    ! m/s
              SFCTH = SFCT * ((1.e5_r8 / SFCP)**0.286_r8)
              LV    = (2.501_r8 &
                   - (2.37e-3_r8 * (SFT(I,J) - 273.16_r8)) ) * 1.e6_r8
              TSV   = SFCTH * (1._r8 + 0.61_r8 * SFCQ)

              !// assume all CTM layers divided into NSUB mass-equal levels
              !// include partial layers in top box (NBL+1) only if mixed
              !// through, i.e., keep only the NSUB layers fully included
              !// in PBL.
              NG  = (NSUB*NBL) + int( (XBL - real(NBL,r8)) * real(NSUB,r8) )
              if (NG .gt. 1) then

                 if  (ng.gt.100) then
                    write(6,'(a,i4,es16.6,i4,3es16.6)') f90file//':'//subr// &
                         ': HOLTSLAG: NG>100',ng,xbl,nbl,mf,zbl,zf
                    CNST = ZBL
                    NBL = 1
                    do L = 1, LPAR
                       if (CNST .ge. 0._r8) then
                          ZLL = CNST
                          NBL = L
                       end if
                       CNST = CNST - DZL(L)
                       if (NBL+1 .ge. L) &
                            write(6,'(a,3i4,es16.6,i4,2es16.6)') f90file//':'//subr// &
                            ': NBL+1 > L:',I,J,L,DZL(L),NBL,CNST,dp(l)
                    end do
                    write(6,'(a,3i4,es16.6,i4,2es16.6)') 'Consider changing 100 closer' &
                         //' to BLPAR or lowering max BL height'
                    if (NG .gt. BLPAR) then
                       write(6,'(a,3i4,es16.6,i4,2es16.6)') f90file//':'//subr// &
                            ': NG > BLPAR; indicating something is wrong.'
                       stop
                    end if
                 end if
                 !// note NMAX changed to NGMAX >= NG, which includes all
                 !// levels*NSUB
                 NBLP1 = min(LPAR, NBL+1)
                 NGMAX = NSUB * NBLP1

                 !// determine PBL profiles from similarity; Ekman solution
                 call INTERP(POFL,TOFL,QVAP,UVEL,VVEL,SFCQ,SFCU,SFCV,SFCTH, &
                      ZMD,ZD,DUDZ,DVDZ,DTHDZ,DQDZ,NSUB,NBLP1)

                 !// loop over tracers
                 QM2(:)  = AIRB(:,II,JJ)
                 do N = 1, NTM
                   !do L = 1, LPAR
                   !  QM(L)  = AIRB(L,II,JJ)
                   !  !// if there is a production P (kg/s): QTT=BTT+P*DTURB
                   !  QTT(L) = BTT(L,N,II,JJ)
                   !  QZT(L) = BZT(L,N,II,JJ)
                   !  QZZ(L) = BZZ(L,N,II,JJ)
                   !  QXZ(L) = BXZ(L,N,II,JJ)
                   !  QYZ(L) = BYZ(L,N,II,JJ)
                   !  QXT(L) = BXT(L,N,II,JJ)
                   !  QXX(L) = BXX(L,N,II,JJ)
                   !  QYT(L) = BYT(L,N,II,JJ)
                   !  QYY(L) = BYY(L,N,II,JJ)
                   !  QXY(L) = BXY(L,N,II,JJ)
                   !end do
                   !// ASH 2017: I think this should be faster
                   QM(:)  = QM2(:)
                   QTT(:) = BTT(:,N,II,JJ)
                   QZT(:) = BZT(:,N,II,JJ)
                   QZZ(:) = BZZ(:,N,II,JJ)
                   QXZ(:) = BXZ(:,N,II,JJ)
                   QYZ(:) = BYZ(:,N,II,JJ)
                   QXT(:) = BXT(:,N,II,JJ)
                   QXX(:) = BXX(:,N,II,JJ)
                   QYT(:) = BYT(:,N,II,JJ)
                   QYY(:) = BYY(:,N,II,JJ)
                   QXY(:) = BXY(:,N,II,JJ)

                   !// #2, 3, or 4 subdivide CTM layers
                   !// call QZONX to slice up tracer/moments into NSUB
                   !// sub-grid layers; ignore the Z-moments in the sub-grid
                   !// layers
                   M  = 0
                   do L = 1, NBLP1
                     call QZONX(QTT(L),QZT(L),QZZ(L),QXZ(L),QYZ(L),QXT(L), &
                          QXX(L),QYT(L),QYY(L),QXY(L),QM(L),PTT,PZT,PZZ, &
                          PXZ,PYZ,PXT,PXX,PYT,PYY,PXY,PM,NSUB)
                     do K = 1, NSUB
                       M      = M+1
                       GM(M)  = PM(K)
                       GTT(M) = PTT(K)
                       GXT(M) = PXT(K)
                       GXX(M) = PXX(K)
                       GXY(M) = PXY(K)
                       GYT(M) = PYT(K)
                       GYY(M) = PYY(K)
                     end do
                     ! put PM(1:3) into GM.
                     !// Set M = 0 before L-loop
                     !K = M + NSUB !// End of GM-index for this layer (3,6,9,...)
                     !M = M + 1    !// Beginning (1,4,7,...)
                     !GM(M:K)  = PM(:)
                     !GTT(M:K) = PTT(:)
                     !GXT(M:K) = PXT(:)
                     !GXX(M:K) = PXX(:)
                     !GXY(M:K) = PXY(:)
                     !GYT(M:K) = PYT(:)
                     !GYY(M:K) = PYY(:)
                     
                   end do

                   !// CGC, CGQ and CGH are removed, since they are all set
                   !// to zero.
                   !// CGH: COUNTER-GRADIENT TERM FOR HEAT
                   !// CGQ: COUNTER-GRADIENT TERM FOR MOISTURE
                   !// CGC: COUNTER-GRADIENT TERM FOR SCALAR C
                   !// call closure model to obtain Kc, countergradient terms
                   !// Kc is evaluated at where DUDZ and DVDZ defined, at
                   !// the edge of a box (ZD).
                   !// ZMD is at the middle of a box, where U, V and
                   !// THETA are defined.
                   call KPROF2 (KC, PRTL1, NGMAX, USTAR, &
                          SFCS, (SFCL/LV), 0._r8, MO_LEN, ZBL, VK, ZD, TSV)

                   !// Save PBL_KEDDY for deposition scalings
                   PBL_KEDDY(II,JJ,MP) = KC(2)

                   !// Adjustment to max from 1:MID-1
                   do L = 1, MID-1
                     KC(L) = max(KC(L), KC(MID))
                   end do

                   !// Save PrandtlL1
                   PrandtlL1(II,JJ,MP) = PRTL1

                   !// KC(M=1:NG-1) = diffusion coefficient (m*m/s)
                   !// between layer M and M+1
                   !do M = 1, NG-1
                   !  DZM(M) = ZMD(M+1) - ZMD(M)
                   !end do
                   do M=1, NG-1
                      DZM = ZMD(M+1) - ZMD(M)
                     KCZ(M) = 0.5_r8 * DTPBL * KC(M) * (GM(M)+GM(M+1)) &
                          / (DZM * DZM * GM(M) * GM(M+1))
                   end do
                   !// set up tridiagonal coeffs as dimensionless
                   !// (not 1/sec), so cannot do s-s
                   !// boundary conditions are zero-flux across layer edges
                   BL(1)  = 0._r8
                   do M = 2, NG-1
                     BL(M) =        -KCZ(M-1) * GM(M)
                   end do
                   BL(NG) =       -KCZ(NG-1) * GM(NG)
                   BD(1)  = 1._r8 + KCZ(1) * GM(2)
                   do M = 2, NG-1
                     BD(M) = 1._r8 + (KCZ(M-1) * GM(M-1) + KCZ(M) * GM(M+1))
                   end do
                   BD(NG) = 1._r8 + KCZ(NG-1) * GM(NG-1)
                   BU(1)  =       -KCZ(1) * GM(1)
                   do M = 2, NG-1
                     BU(M) =        -KCZ(M) * GM(M)
                   end do
                   BU(NG) = 0._r8

                   !// the right-hand side is just the initial values of
                   !// GTT (in kg)

                   !// call tridiagonal solver to vertically diffuse mass
                   !// and 5 XY moments
                   call TRIDIAG (BL,BD,BU,GTT,GOUT,NG)
                   do M = 1, NG
                      GTT(M) = GOUT(M)
                   end do
                   call TRIDIAG (BL,BD,BU,GXT,GOUT,NG)
                   do M = 1, NG
                      GXT(M) = GOUT(M)
                   end do
                   call TRIDIAG (BL,BD,BU,GXX,GOUT,NG)
                   do M = 1, NG
                      GXX(M) = GOUT(M)
                   end do
                   call TRIDIAG (BL,BD,BU,GYT,GOUT,NG)
                   do M = 1, NG
                      GYT(M) = GOUT(M)
                   end do
                   call TRIDIAG (BL,BD,BU,GYY,GOUT,NG)
                   do M = 1, NG
                      GYY(M) = GOUT(M)
                   end do
                   call TRIDIAG (BL,BD,BU,GXY,GOUT,NG)
                   do M = 1, NG
                      GXY(M) = GOUT(M)
                   end do
 
                   !// call QXZON to recombine turbulent grid to large-scale
                   !// grid w/moments
                   QTMPT(:) = QTT(:)
                   QTMPM(:) = QM(:)
                   M    = 0
                   do L = 1, NBLP1
                     do K = 1, NSUB
                       M = M+1
                       PM(K)  = GM(M)
                       PTT(K) = GTT(M)
                       PXT(K) = GXT(M)
                       PXX(K) = GXX(M) 
                       PXY(K) = GXY(M)
                       PYT(K) = GYT(M)
                       PYY(K) = GYY(M)
                       PZT(K) = 0._r8
                       PZZ(K) = 0._r8
                       PXZ(K) = 0._r8
                       PYZ(K) = 0._r8
                     end do
                     call QXZON(PTT,PZT,PZZ,PXZ,PYZ,PXT,PXX,PYT,PYY,PXY,PM, &
                          NSUB,QTT(L),QZT(L),QZZ(L),QXZ(L),QYZ(L), &
                          QXT(L),QXX(L),QYT(L),QYY(L),QXY(L),QM(L), &
                          NBLP1,L,N)
                   end do

                   !// pass resultant concentrations/moments back to CTM grid
                   !do L = 1, LPAR
                   !  BTT(L,N,II,JJ) = QTT(L)
                   !  BZT(L,N,II,JJ) = QZT(L)
                   !  BZZ(L,N,II,JJ) = QZZ(L)
                   !  BXZ(L,N,II,JJ) = QXZ(L)
                   !  BYZ(L,N,II,JJ) = QYZ(L)
                   !  BXT(L,N,II,JJ) = QXT(L)
                   !  BXX(L,N,II,JJ) = QXX(L)
                   !  BYT(L,N,II,JJ) = QYT(L)
                   !  BYY(L,N,II,JJ) = QYY(L)
                   !  BXY(L,N,II,JJ) = QXY(L)
                   !end do            ! L 
                   !// ASH 2017: I think this should be faster
                   BTT(:,N,II,JJ) = QTT(:)
                   BZT(:,N,II,JJ) = QZT(:)
                   BZZ(:,N,II,JJ) = QZZ(:)
                   BXZ(:,N,II,JJ) = QXZ(:)
                   BYZ(:,N,II,JJ) = QYZ(:)
                   BXT(:,N,II,JJ) = QXT(:)
                   BXX(:,N,II,JJ) = QXX(:)
                   BYT(:,N,II,JJ) = QYT(:)
                   BYY(:,N,II,JJ) = QYY(:)
                   BXY(:,N,II,JJ) = QXY(:)


                 end do !// N

              end if !// if (NG .gt. 1) then
           end if !// if (ZBL .ge. 1._r8) then

         end do !// J
       end do !// I
    else
       write(6,'(a,i4,a)') f90file//':'//subr//': NBLX=',NBLX, &
            ' not defined - stopping!'
       stop
    end if


    !// --------------------------------------------------------------------
  end subroutine CNVDBL
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine BULK(QTT,QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY,QM,XBL,NQ)
    !// --------------------------------------------------------------------
    !// Computes tracer profile (w/moments) for a fully mixed PBL up to
    !// layer = XBL
    !// e.g., XBL = 1.25 --> mixed all of L=1 and lowest 25% of L=2 by mass,
    !// erase all vertical moments (QZT,QZZ,QXZ,QYZ) in mixed regions
    !// keep all horizontal information (QXT,QYT,QXX,QYY,QXY)
    !// dates to MJP 7/97
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NQ
    real(r8), intent(in) :: XBL
    !// Input/Output
    real(r8), dimension(NQ), intent(inout) :: &
         QTT,QZT,QZZ,QXZ,QYZ,QXT,QXX,QYT,QYY,QXY,QM

    !// Locals
    integer :: L, NBL
    real(r8) :: FTT,FXT,FXX,FYT,FYY,FXY,FMM,FBL,FMIN
    real(r8) :: FTTL,FXTL,FXXL,FYTL,FYYL,FXYL,FMML
    real(r8) :: EPS,EPS1,EPS1Q,EPS0
    !// --------------------------------------------------------------------

    if (XBL .gt. 0.01_r8) then

       NBL   = int(XBL)
       FBL   = XBL - real(NBL,r8)
       FTT   = 0._r8
       FXT   = 0._r8
       FXX   = 0._r8
       FYT   = 0._r8
       FYY   = 0._r8
       FXY   = 0._r8
       FMM   = 0._r8
       FMML  = 0._r8

       !// Add up air mass and tracer moments in boxes 1:NBL
       do L = 1, NBL
          FTT = FTT + QTT(L)
          FXT = FXT + QXT(L)
          FXX = FXX + QXX(L)
          FYT = FYT + QYT(L)
          FYY = FYY + QYY(L)
          FXY = FXY + QXY(L)
          FMM = FMM + QM(L)
       end do

       !// do fractional layer >NBL only if >1% of layer affected
       if (FBL .gt. 0.01_r8)  then

          !// restrain moments in fractional box (LIMIT=2)
          L = NBL+1
          call QLIMIT2 (QTT(L),QZT(L),QZZ(L),QXZ(L),QYZ(L))
        
          !// compute air mass, tracer moments in lower parcel of L
          !// (into mixed layer)
          FMML  = FBL*QM(L)          ! same as advective step QU = FBL*QM(L)
          EPS   = FBL
          EPS1  = 1._r8 - EPS
          EPS1Q = EPS1*EPS1
          FTTL  = EPS*(QTT(L) - EPS1*(QZT(L) - (EPS1-EPS)*QZZ(L)))
          FXTL  = EPS*(QXT(L) - EPS1*QXZ(L))
          FYTL  = EPS*(QYT(L) - EPS1*QYZ(L))
          FXXL  = EPS*QXX(L)
          FYYL  = EPS*QYY(L)
          FXYL  = EPS*QXY(L)

          !// readjust air mass, tracer moments remaining in top of layer L
          QM(L)  = QM(L) - FMML
          QTT(L) = QTT(L) - FTTL
          QZT(L) = EPS1Q*(QZT(L) + 3._r8*EPS*QZZ(L))
          QZZ(L) = EPS1*EPS1Q*QZZ(L)
          QXT(L) = QXT(L) - FXTL
          QXX(L) = QXX(L) - FXXL
          QYT(L) = QYT(L) - FYTL
          QYY(L) = QYY(L) - FYYL
          QXZ(L) = EPS1Q*QXZ(L)
          QYZ(L) = EPS1Q*QYZ(L)
          QXY(L) = QXY(L) - FXYL

          !// add partial layer to sum of lower layers (1:NBL)
          FTT = FTT + FTTL
          FXT = FXT + FXTL
          FXX = FXX + FXXL
          FYT = FYT + FYTL
          FYY = FYY + FYYL
          FXY = FXY + FXYL
          FMM = FMM + FMML

       end if

       !// calculate mean mass and moments in mixed PBL & distribute
       FMIN  = 1._r8 / FMM
       FTT   = FTT * FMIN
       FXT   = FXT * FMIN
       FXX   = FXX * FMIN
       FYT   = FYT * FMIN
       FYY   = FYY * FMIN
       FXY   = FXY * FMIN
       do L = 1, NBL
          QTT(L) = FTT * QM(L)
          QXT(L) = FXT * QM(L)
          QXX(L) = FXX * QM(L)
          QYT(L) = FYT * QM(L)
          QYY(L) = FYY * QM(L)
          QXY(L) = FXY * QM(L)
          QZT(L) = 0._r8
          QZZ(L) = 0._r8
          QXZ(L) = 0._r8
          QYZ(L) = 0._r8
       end do

       if (FBL .gt. 0.01_r8) then
          L = NBL+1
          FTTL = FTT * FMML
          FXTL = FXT * FMML
          FXXL = FXX * FMML
          FYTL = FYT * FMML
          FYYL = FYY * FMML
          FXYL = FXY * FMML
          QM(L)  = QM(L) + FMML
          EPS    = FMML / QM(L)
          EPS1   = 1._r8 - EPS
          EPS0   = EPS*QTT(L) - EPS1*FTTL
          QZZ(L) = EPS1*EPS1*QZZ(L) &
                   + 5._r8*(EPS*EPS1*QZT(L) - (EPS1-EPS)*EPS0)
          QZT(L) = EPS1*QZT(L) + 3._r8*EPS0
          QTT(L) = QTT(L) + FTTL
          QXZ(L) = EPS1*QXZ(L) + 3._r8*(-EPS1*FXTL + EPS*QXT(L))
          QYZ(L) = EPS1*QYZ(L) + 3._r8*(-EPS1*FYTL + EPS*QYT(L))
          QXT(L) = QXT(L) + FXTL
          QXX(L) = QXX(L) + FXXL
          QYT(L) = QYT(L) + FYTL
          QYY(L) = QYY(L) + FYYL
          QXY(L) = QXY(L) + FXYL
       end if

    end if

    !// --------------------------------------------------------------------
  end subroutine BULK
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_keddy_L1(keddy, PRTL1, USTAR, &
       HEATV, EVAP, MOL, PBLHGT, VK, ZT, TSV )
    !// --------------------------------------------------------------------
    !// Description: 
    !//  Calculates diffusivities and countergradient terms used in 
    !//  PBL closure.   
    !//
    !// Amund Sovde Haslerud, November 2017
    !//  The original KPROF differs from this routine: In KPROF2, onlt
    !//  KC is calculated, not KH, KM, KQ. This is to make the routine
    !//  faster.
    !//
    !// Method:
    !//   Holtslag and Boville (Mon. Wea. Rvw., 1990)  
    !//
    !// History: 
    !// Version   Date     Comment 
    !// -------   ----     ------- 
    !//  1.0    06/07/99   Original code. Bryan Hannegan 
    !//
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: G0
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: &
         ZT, PBLHGT, & ! layer top height, pbl height
         USTAR, HEATV, EVAP, MOL, VK, TSV
    !// Output
    real(r8), intent(out) :: &
         keddy, &  !// KH for surface layer
         PRTL1     !// Return PRANDTL for surface layer

    !// Local parameters
    real(r8) :: &
         PRANDTL, &
         PBLK, WSTAR, WTURB, &
         ZH, ZL, ZZH, ZLSURF, ZPHIM, &
         ZMOL, ZBL
    integer :: K
    logical :: LUNSTABLE

    real(r8), parameter :: &
         CONST = 7.2_r8, &
         BGK = 1.e-19_r8, &
         GRAV = G0
    !// --------------------------------------------------------------------

    !// calculate inverse of Monin-Obukov Length
    ZMOL = 1._r8 / MOL
    !// calculate inverser of BL Height
    ZBL  = 1._r8 / PBLHGT
    !// Unstable or not?
    LUNSTABLE = HEATV .gt. 0._r8
    !// in case of unstable BL find W^*
    if (LUNSTABLE) &
         WSTAR = (GRAV * HEATV * PBLHGT / TSV)**(1._r8 / 3._r8)

    !// Initialise PRANDTL for K=1 (PRTL1)
    PRTL1 = 0.85_r8


    !// Estimate K in lowest model layer. Always within PBL height.
    ZH = ZT * ZBL
    if (ZH .gt. 1._r8) then
       ZZH = 0._r8
    else
       ZZH = VK * ZT * (1._r8 - ZH) * (1._r8 - ZH)
    end if

    !// evaluate based on stability
    !// unstable PBL
    if (LUNSTABLE) then
       if (ZH .lt. 0.1_r8) then  
          !// surface layer, compute corrections at actual height
          ZL      = ZT * ZMOL
          ZPHIM   = 1._r8 / PHIM(ZL)

          PBLK    = ZZH * USTAR * ZPHIM
          PRANDTL = PHIH(ZL) * ZPHIM

          !// Set PRTL1
          PRTL1 = PRANDTL

       else
          !// outer layer, compute stability corrections at top of
          !// sfc layer (1/10 PBL)  
          ZLSURF  = 0.1_r8 * PBLHGT * ZMOL
          ZPHIM   = 1._r8 / PHIM(ZLSURF)
          WTURB   = USTAR * ZPHIM
          PRANDTL = PHIH(ZLSURF) * ZPHIM &
               + VK * 0.1_r8 * (CONST * WSTAR) / WTURB
          PBLK    = WTURB * ZZH

          !// Set PRTL1
          PRTL1 = PRANDTL
       end if

       keddy = max( PBLK / PRANDTL, BGK )
       !// stable and neutral PBL
    else
       ZL   = ZT * ZMOL
       PBLK = ZZH * USTAR / PHIM(ZL)
       keddy = max( PBLK, BGK )

       !// Set PRTL1
       !// For stable/neutral we approximate KM=KH, i.e. PRANDTL=1
       PRTL1 = 1._r8

    end if

    !// --------------------------------------------------------------------
  end subroutine get_keddy_L1
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine KPROF2(KC, PRTL1, NTOP, USTAR, &
       HEATV, EVAP, WCSFC, MOL, PBLHGT, VK, ZT, TSV )
    !// --------------------------------------------------------------------
    !// Description: 
    !//  Calculates diffusivities and countergradient terms used in 
    !//  PBL closure.   
    !//
    !// Method:
    !//   Holtslag and Boville (Mon. Wea. Rvw., 1990)  
    !//
    !// Amund Sovde Haslerud, November 2017
    !//  The original KPROF differs from this routine: In KPROF2, onlt
    !//  KC is calculated, not KH, KM, KQ. This is to make the routine
    !//  faster.
    !//
    !// History: 
    !// Version   Date     Comment 
    !// -------   ----     ------- 
    !//  1.0    06/07/99   Original code. Bryan Hannegan 
    !//
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: G0
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  NTOP
    real(r8), intent(in) :: &
         ZT(NTOP), PBLHGT, &
         USTAR, HEATV, EVAP, WCSFC, MOL, VK, TSV
    !// Output
    real(r8), intent(out) :: &
         KC(NTOP)
    real(r8), intent(out) :: PRTL1 !// Return PRANDTL for surface layer

    !// Local parameters
    real(r8) :: &
         PRANDTL, &
         PBLK, WSTAR, WTURB, &
         ZH, ZL, ZZH, ZLSURF, ZPHIM, &
         ZMOL, ZBL
    integer :: K
    logical :: LUNSTABLE

    real(r8), parameter :: &
         CONST = 7.2_r8, &
         BGK = 1.e-19_r8, &
         GRAV = G0
    !// --------------------------------------------------------------------

    !// Initialize
    !// set diffusivities to background
    !KC(:)  = BGK
    
    !// calculate inverse of Monin-Obukov Length
    ZMOL = 1._r8 / MOL
    !// calculate inverser of BL Height
    ZBL  = 1._r8 / PBLHGT
    !// Unstable or not?
    LUNSTABLE = HEATV .gt. 0._r8
    !// in case of unstable BL find W^*
    if (LUNSTABLE) &
         WSTAR = (GRAV * HEATV * PBLHGT / TSV)**(1._r8 / 3._r8)

    !// Initialise PRANDTL for K=1 (PRTL1)
    PRTL1 = 0.85_r8


    !// set up K-profile function
    KC(NTOP) = BGK
    do K = 1, NTOP - 1

       !// Only calculate as long as within PBL
       if (ZT(K) .lt. PBLHGT) then
          ZH = ZT(K) * ZBL
          if (ZH .gt. 1._r8) then
             ZZH = 0._r8
          else
             ZZH = VK * ZT(K) * (1._r8 - ZH) * (1._r8 - ZH)
          end if

          !// evaluate based on stability
          !// unstable PBL
          if (LUNSTABLE) then
             if (ZH .lt. 0.1_r8) then  
                !// surface layer, compute corrections at actual height
                ZL      = ZT(K) * ZMOL
                ZPHIM   = 1._r8 / PHIM(ZL)

                PBLK    = ZZH * USTAR * ZPHIM
                PRANDTL = PHIH(ZL) * ZPHIM

                !// Set PRTL1
                if (K .eq. 1) PRTL1 = PRANDTL

             else
                !// outer layer, compute stability corrections at top of
                !// sfc layer (1/10 PBL)  
                ZLSURF  = 0.1_r8 * PBLHGT * ZMOL
                ZPHIM   = 1._r8 / PHIM(ZLSURF)
                WTURB   = USTAR * ZPHIM
                PRANDTL = PHIH(ZLSURF) * ZPHIM &
                     + VK * 0.1_r8 * (CONST * WSTAR) / WTURB
                PBLK    = WTURB * ZZH

                !// Set PRTL1
                if (K .eq. 1) PRTL1 = PRANDTL
             end if

             KC(K) = max( PBLK / PRANDTL, BGK )
             !// stable and neutral PBL
          else
             ZL   = ZT(K) * ZMOL
             PBLK = ZZH * USTAR / PHIM(ZL)
             KC(K) = max( PBLK, BGK )

             !// Set PRTL1
             !// For stable/neutral we approximate KM=KH, i.e. PRANDTL=1
             if (K .eq. 1) PRTL1 = 1._r8

          end if

       else
          !// Set background value above PBL.
          KC(K) = BGK
       end if
    end do
    !// --------------------------------------------------------------------
  end subroutine KPROF2
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine KPROF(KC, PRTL1, NTOP, USTAR, &
       HEATV, EVAP, WCSFC, MOL, PBLHGT, VK, ZT, TSV )
    !// --------------------------------------------------------------------
    !// Description: 
    !//  Calculates diffusivities and countergradient terms used in 
    !//  PBL closure.   
    !//
    !// Method:
    !//   Holtslag and Boville (Mon. Wea. Rvw., 1990)  
    !//
    !// History: 
    !// Version   Date     Comment 
    !// -------   ----     ------- 
    !//  1.0    06/07/99   Original code. Bryan Hannegan 
    !//
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_parameters, only: G0
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  NTOP
    real(r8), intent(in) :: &
         ZT(NTOP), PBLHGT, &
         USTAR, HEATV, EVAP, WCSFC, MOL, VK, TSV
    !// Output
    real(r8), intent(out) :: &
         KC(NTOP)
    real(r8), intent(out) :: PRTL1 !// Return PRANDTL for surface layer

    !// Local parameters
    real(r8) :: &
         KM(BLPAR), KH(BLPAR), KQ(BLPAR), &
         CGHA, CGWA, CGCA, &
         PRANDTL, &
         PBLK, WSTAR, WTURB, &
         ZH, ZL, ZZH, ZLSURF, ZPHIM, &
         ZMOL, ZBL
    integer :: K

    real(r8), parameter :: &
         CONST = 7.2_r8, &
         BGK = 1.e-19_r8, &
         GRAV = G0
    !// --------------------------------------------------------------------

    !// Initialize
    !// set diffusivities to background; counter-gradient terms to zero
    KM(:)  = BGK
    KH(:)  = BGK
    KQ(:)  = BGK
    KC(:)  = BGK
    !// CGH(:) = 0._r8
    !// CGQ(:) = 0._r8
    !// CGC(:) = 0._r8
    
    !// calculate inverse of Monin-Obukov Length
    ZMOL = 1._r8 / MOL
    !// calculate inverser of BL Height
    ZBL  = 1._r8 / PBLHGT
    !// in case of unstable BL find W^*
    if (HEATV .gt. 0._r8) &
         WSTAR = (GRAV * HEATV * PBLHGT / TSV)**(1._r8 / 3._r8)

    !// Initialise PRANDTL for K=1 (PRTL1)
    PRTL1 = 0.85_r8

    !// set up K-profile function
    do K = 1, NTOP-1

       if (ZT(K) .lt. PBLHGT) then
          ZH = ZT(K) * ZBL
          if (ZH .gt. 1._r8) then
             ZZH = 0._r8
          else
             ZZH = VK * ZT(K) * (1._r8 - ZH) * (1._r8 - ZH)
          end if

          !// evaluate based on stability
          !// unstable PBL
          if (HEATV .gt. 0._r8) then
             if (ZH .lt. 0.1_r8) then  
                !// surface layer, compute corrections at actual height
                ZL      = ZT(K) * ZMOL
                ZPHIM   = 1._r8 / PHIM(ZL)

                PBLK    = ZZH * USTAR * ZPHIM
                PRANDTL = PHIH(ZL) * ZPHIM
                !// CGHA    = 0._r8
                !// CGWA    = 0._r8
                !// CGCA    = 0._r8        

                !// Set PRTL1
                if (K .eq. 1) PRTL1 = PRANDTL

             else
                !// outer layer, compute stability corrections at top of
                !// sfc layer (1/10 PBL)  
                ZLSURF  = 0.1_r8 * PBLHGT * ZMOL
                ZPHIM   = 1._r8 / PHIM(ZLSURF)
                WTURB   = USTAR * ZPHIM
                PRANDTL = PHIH(ZLSURF) * ZPHIM &
                     + VK * 0.1_r8 * (CONST * WSTAR) / WTURB
                PBLK    = WTURB * ZZH
                !// CGHA    = HEATV*(CONST*WSTAR)/(WTURB*WTURB*PBLHGT)
                !// CGWA    = EVAP *(CONST*WSTAR)/(WTURB*WTURB*PBLHGT)
                !// CGCA    = WCSFC*(CONST*WSTAR)/(WTURB*WTURB*PBLHGT)

                !// Set PRTL1
                if (K .eq. 1) PRTL1 = PRANDTL
             end if

             if (PBLK .gt. BGK) KM(K)  = PBLK
             if ((PBLK / PRANDTL) .gt. BGK) KH(K)  = (PBLK/PRANDTL)
             KQ(K)  = KH(K)
             KC(K)  = KH(K)
             !// CGH(K) = CGHA
             !// CGQ(K) = CGWA
             !// CGC(K) = CGCA

             !// stable and neutral PBL
          else
             ZL   = ZT(K) * ZMOL
             PBLK = ZZH * USTAR / PHIM(ZL)
             if (PBLK .gt. BGK) then
                KM(K)  = PBLK
                KH(K)  = PBLK
                KQ(K)  = PBLK
                KC(K)  = PBLK
             end if
             !// CGH(K) = 0._r8
             !// CGQ(K) = 0._r8
             !// CGC(K) = 0._r8

             !// Set PRTL1
             !// For stable/neutral we approximate KM=KH, i.e. PRANDTL=1
             if (K .eq. 1) PRTL1 = 1._r8

          end if
          
       end if
    end do
    !// --------------------------------------------------------------------
  end subroutine KPROF
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  real(r8) function PHIM(ZETA)
    !// --------------------------------------------------------------------
    !// Description: 
    !//   Calculates similarity theory stability correction for momentum  
    !//
    !// Method:
    !//   From Monin-Obukhov length, use Businger-Dyer relationship 
    !//
    !// Owner:  Bryan Hannegan (bjhanneg@uci.edu)
    !//
    !// History: 
    !//     Version   Date     Comment 
    !//     -------   ----     ------- 
    !//      1.0    05/24/99   Original code. Bryan Hannegan 
    !//
    !// Code Description: 
    !//     Language:           Fortran 77. 
    !//     Software Standards: Adapted from "European Standards for 
    !//       Writing and Documenting Exchangeable Fortran 90 Code". 
    !//         (though not entirely!) 
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// variables passed in function call
    real(r8) :: ZETA         ! M-O length
    !// calculate stability correction (c.f. Stull, p. 383-384)
    !// scaling parameters BETAM and GAMM specified in params.h
    real(r8) , parameter :: &
         GAMM  = 15.0_r8, &
         BETAM =  4.7_r8 
    !// --------------------------------------------------------------------
    if (ZETA .lt. 0._r8) then
       PHIM = 1._r8 / ((1._r8 - (GAMM * ZETA))**0.25_r8)
    else
       PHIM = (1._r8 + (BETAM * ZETA))
    end if
    !// --------------------------------------------------------------------
  end function PHIM
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  real(r8) function PHIH(ZETA)
    !// --------------------------------------------------------------------
    !// Description: 
    !//   Calculates similarity theory stability correction for heat 
    !//
    !// Method:
    !//   From Monin-Obukhov length, use Businger-Dyer relationship 
    !//
    !// Owner:  Bryan Hannegan (bjhanneg@uci.edu)
    !//
    !// History: 
    !//     Version   Date     Comment 
    !//     -------   ----     ------- 
    !//      1.0    05/24/99   Original code. Bryan Hannegan 
    !//
    !// Code Description: 
    !//     Language:           Fortran 77. 
    !//     Software Standards: Adapted from "European Standards for 
    !//       Writing and Documenting Exchangeable Fortran 90 Code". 
    !//         (though not entirely!) 
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// variables passed in function call
    real(r8) :: ZETA         ! M-O length
    !// calculate stability correction (cf. Stull, p. 383-384)
    !// scaling parameters BETAH and GAMH specified in params.h
    !// Prandtl number PRT specified in params.h
    real(r8), parameter :: &
         PRT = 0.74_r8, &
         GAMH = 9._r8, &
         BETAH = 4.7_r8
    !// --------------------------------------------------------------------
    if (ZETA .lt. 0._r8) then
       PHIH = PRT / sqrt(1._r8 - (GAMH * ZETA))
    else
       PHIH = PRT * (1._r8 + (BETAH * ZETA))
    end if
    !// --------------------------------------------------------------------
  end function PHIH
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine TRIDIAG(A,B,C,R,U,N)
    !// --------------------------------------------------------------------
    !// Solves tridiagonal system according to Numerical Recipes (1992)    
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: N
    real(r8), intent(in) :: A(*),B(*),C(*),R(*)
    !// Output
    real(r8), intent(out) :: U(*)
    !// Locals
    real(r8) :: GAM(N), BET
    integer :: J
    !// --------------------------------------------------------------------

    BET  = 1._r8 / B(1)
    U(1) = R(1) * BET
    do J = 2, N
       GAM(J) = C(J-1) * BET
       BET    = 1._r8 / (B(J) - A(J) * GAM(J))
       U(J) = (R(J) - A(J) * U(J-1)) * BET
    end do
    do J = N-1, 1, -1
       U(J) = U(J) - (GAM(J+1) * U(J+1))
    end do
    !// --------------------------------------------------------------------
  end subroutine TRIDIAG
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine INTERP(POFL,TOFL,QVAP,UVEL,VVEL,SFCQ,SFCU,SFCV,SFCTH, &
       ZMD,ZD,DUDZ,DVDZ,DTHDZ,DQDZ,NSUB,NBLP1)
    !// --------------------------------------------------------------------
    !// mass weighted interpalation routine
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit  none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NSUB, NBLP1
    real(r8), dimension(LPAR+1), intent(in) :: POFL(LPAR+1)
    real(r8), dimension(LPAR), intent(in) :: TOFL, QVAP, UVEL, VVEL
    real(r8), intent(in) :: SFCQ, SFCU, SFCV, SFCTH
    !// Output
    real(r8), dimension(BLPAR), intent(out) :: &
         DUDZ, DVDZ, DTHDZ, DQDZ, ZMD, ZD

    !// Locals
    integer, parameter ::  LD = BLPAR, LMAX = LPAR
    real(r8), dimension(LD) ::  THTAT, QT, UT, VT, TVT, DZ, DZM
    real(r8), dimension(LMAX) ::  MWTL, MWTR, DP, THTA, TV
    real(r8) :: YQ, YU, YV, YTV, ZSUB, ZDZ1, P1, P2, PMID
    integer :: L, M, K, ML, MLOW, MUP, MID
    !// --------------------------------------------------------------------

    do L = 1, NSUB * NBLP1
       DUDZ(L)   = 0._r8
       DVDZ(L)   = 0._r8
       DTHDZ(L)  = 0._r8
       DQDZ(L)   = 0._r8
       ZMD(L)    = 0._r8
       ZD(L)     = 0._r8
    end do
    ZSUB  = 1._r8 / real(NSUB, r8)
    MID   = int(NSUB/2) + 1
    K     = 0
    do L = 1,NBLP1
       DP(L)   = POFL(L) - POFL(L+1)
       PMID    = 0.5_r8 * (POFL(L) + POFL(L+1))
       THTA(L) = TOFL(L) * ((1.e3_r8 / PMID)**0.286_r8)
       TV(L)   = THTA(L) * (1._r8 + 0.61_r8 * QVAP(L))
       do M = 1, NSUB
          K     = K + 1
          P1    = POFL(L) - real(M, r8) * ZSUB * DP(L)
          P2    = POFL(L) - real(M-1, r8) * ZSUB * DP(L)
          !// R = 287.*(1-Q) + 461.5*Q -- assume 0.5% bl w.v.==> R = 288.
          !// delta-z = dln(P) * R * T / g   where R/g = 288/9.81 = 29.36
          DZ(K) = log(P2/P1) * TOFL(L) * 29.36_r8
       end do
    end do
    DZM(1)    = 0.5_r8 * DZ(1)
    ZMD(1)    = 0.5_r8 * DZ(1)
    ZD(1)     = DZ(1)
    do K = 2, NBLP1 * NSUB
       DZM(K)  = 0.5_r8 * (DZ(K-1) + DZ(K))
       ZMD(K)  = ZMD(K-1) + DZM(K)
       ZD(K)   = ZD(K-1) + DZ(K)
    end do
    do L = 1, NBLP1-1
       MWTR(L) = DP(L) / (DP(L) + DP(L+1))
       MWTL(L) = DP(L+1) / (DP(L) + DP(L+1))
    end do

    do L = 1, NBLP1-1
       ML    = L * NSUB
       !// calculate boundary values at edge of two layers --------------
       YQ    = QVAP(L)*MWTL(L) + QVAP(L+1)*MWTR(L)
       YU    = UVEL(L)*MWTL(L) + UVEL(L+1)*MWTR(L)
       YV    = VVEL(L)*MWTL(L) + VVEL(L+1)*MWTR(L)
       YTV   = TV(L)*MWTL(L) + TV(L+1)*MWTR(L)
       !// MID value in layer L (NSUB=odd assumption is used here) -------
       K        = ML + 1 - MID
       QT(K)    = QVAP(L)
       UT(K)    = UVEL(L)
       VT(K)    = VVEL(L)
       TVT(K)   = TV(L)
       K        = 0
       do M = 1,NSUB-1,2
          K      = K + 1
          !// upper part sub-layers in layer L -----------------------------
          MUP       = ML + 1 - K
          QT(MUP)   = YQ + real(M,r8)*ZSUB*(QVAP(L)-YQ)
          UT(MUP)   = YU + real(M,r8)*ZSUB*(UVEL(L)-YU)
          VT(MUP)   = YV + real(M,r8)*ZSUB*(VVEL(L)-YV)
          TVT(MUP)  = YTV + real(M,r8)*ZSUB*(TV(L)-YTV)
          !// lower part sub-layers in layer L+1 ---------------------------
          MLOW      = ML + K
          QT(MLOW)  = YQ + real(M,r8)*ZSUB*(QVAP(L+1)-YQ)
          UT(MLOW)  = YU + real(M,r8)*ZSUB*(UVEL(L+1)-YU)
          VT(MLOW)  = YV + real(M,r8)*ZSUB*(VVEL(L+1)-YV)
          TVT(MLOW) = YTV + real(M,r8)*ZSUB*(TV(L+1)-YTV)
       end do
    end do
    !// upper boundary sub-layers ------------------------------------
    K         = -1
    ML        = NBLP1 * NSUB
    do M = 1, MID
       K         = K + 1
       QT(ML-K)  = QVAP(NBLP1)
       UT(ML-K)  = UVEL(NBLP1)
       VT(ML-K)  = VVEL(NBLP1)
       TVT(ML-K) = TV(NBLP1)
    end do
    do M = MID,NBLP1*NSUB
       THTAT(M) = TVT(M) / (1._r8 + 0.61_r8 * QT(M))
    end do
    !// use surface value to interpalate the bottom sub-layer mid point value
    do M = 1, MID-1
       P2       = ZMD(M) / ZMD(MID)
       P1       = (ZMD(MID) - ZMD(M)) / ZMD(MID)
       QT(M)    = P1*SFCQ + P2*QT(MID)
       UT(M)    = P1*SFCU + P2*UT(MID)
       VT(M)    = P1*SFCV + P2*VT(MID)
       THTAT(M) = P1*SFCTH + P2*THTAT(MID)
    end do
    !// derivatives are defined at the top edge of sub-layers
    do M = 1, NBLP1 * NSUB-1
       ZDZ1     = 1._r8 / DZM(M+1)
       DQDZ(M)  = (QT(M+1)-QT(M)) * ZDZ1
       DUDZ(M)  = (UT(M+1)-UT(M)) * ZDZ1
       DVDZ(M)  = (VT(M+1)-VT(M)) * ZDZ1
       DTHDZ(M) = (THTAT(M+1)-THTAT(M)) * ZDZ1
    end do

    !// --------------------------------------------------------------------
  end subroutine INTERP
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine QXZON(QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ,QM, NQ, &
                   ETT,EXT,EXX,EXY,EXZ,EYT,EYY,EZT,EZZ,EYZ,EM, &
                   NBLP1,LL,NN)
    !// --------------------------------------------------------------------
    !// Combines NQ equal mass boxes (QTT,QXT,QXX,...,QYX) aligned in
    !// into 1 Extended zone (ETT,EXT,EXX,...,EYZ)
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: &
          NQ, &           ! # of boxes in initial distribution to be combined
          LL,NN,NBLP1
    real(r8), dimension(NQ), intent(in) :: &
         QTT, &             ! tracer mass in box [I]
         QM, &              ! dry air mass in box [I]
         QXT,QXX, &      ! 1st & 2nd moments of tracer in X-direct
         QXY,QXZ, &      ! coupled 1st moments including X component
         QYZ,QYT,QYY, &
         QZT,QZZ       ! other SOMs without X
    !// Output
    real(r8), intent(out) :: &
         ETT, &                ! tracer mass in extended zone
         EXT,EXX,EXY,EXZ,EYZ,EYT,EYY,EZT,EZZ, & ! moments in extended zone
         EM                   ! dry air mass in extended zone
    !// Locals
    real(r8) :: &
         QTTL(BLPAR),SUMT,SUMF, ZNQ,FLTIQ
    integer :: N1N2, IQ
    !// --------------------------------------------------------------------

    !// Initialize
    N1N2 = (NQ+1) * (NQ+2)
    ZNQ  = 1._r8 / real(NQ, r8)
    ETT  = 0._r8
    EXT  = 0._r8
    EXX  = 0._r8
    EXY  = 0._r8
    EXZ  = 0._r8
    EYZ  = 0._r8
    EYT  = 0._r8
    EYY  = 0._r8
    EZT  = 0._r8
    EZZ  = 0._r8
    EM   = 0._r8
    !// since assume equal-mass boxes, rescale STT() so that no spurious
    !// moments from variations in AIRD() with constant STT/AIR
    SUMT = 0._r8
    SUMF = 0._r8
    do IQ = 1, NQ
       SUMT     = SUMT + QTT(IQ)
       QTTL(IQ) = QTT(IQ)
       if (QM(IQ) .gt. 0._r8) then
          QTTL(IQ) = QTT(IQ) / QM(IQ)
       else
          write(6,*) '*** QM(IQ).LT.0. in QXZON ***',IQ,NQ
          stop
       end if
       SUMF    = SUMF + QTTL(IQ)
    end do
    if (SUMF .gt. 0._r8) then
       SUMT = SUMT / SUMF
    else if (sumf .lt. 0._r8) then
       !// Allow tracer to be negative
       if (SUMT .le. 0) then
          SUMT = SUMT / SUMF
       else
          write(6,*) '*** SUMF.LE.0. in QXZON ***', SUMF,NQ,LL,NN,NBLP1
       end if
    else
       !// Do nothing for sumf=0, then sumt=0 also
    end if
    do IQ = 1, NQ
       QTTL(IQ) = QTTL(IQ) * SUMT
    end do

    !// sum up for all boxes in EPZ
    do IQ = 1, NQ
       FLTIQ = real(2*IQ - NQ - 1, r8)
       ETT   = ETT + QTTL(IQ)
       EXT   = EXT + QXT(IQ) + (3._r8*FLTIQ)*QTTL(IQ)
       EYT   = EYT + QYT(IQ)
       EZT   = EZT + QZT(IQ)
       EXX   = EXX + QXX(IQ) + 5._r8*QXT(IQ)*FLTIQ &
               + QTTL(IQ) * real(5*(6*IQ*(IQ-NQ-1)+N1N2), r8)
       EYY   = EYY + QYY(IQ)
       EZZ   = EZZ + QZZ(IQ)
       EXY   = EXY + QXY(IQ) + (3._r8*FLTIQ)*QYT(IQ)
       EXZ   = EXZ + QXZ(IQ) + (3._r8*FLTIQ)*QZT(IQ)
       EYZ   = EYZ + QYZ(IQ)
       EM    = EM + QM(IQ)
    end do
    !// average
    EXT = EXT*ZNQ
    EXX = EXX*ZNQ*ZNQ
    EXY = EXY*ZNQ
    EXZ = EXZ*ZNQ

    !// --------------------------------------------------------------------
  end subroutine QXZON
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine QZONX (ETT,EXT,EXX,EXY,EXZ,EYT,EYY,EZT,EZZ,EYZ,EM, &
                    QTT,QXT,QXX,QXY,QXZ,QYT,QYY,QZT,QZZ,QYZ,QM, &
                    NQ )
    !// --------------------------------------------------------------------
    !// Divide one Extended zone (ETT,EXT,EXX,...,EYZ) into equal-mass
    !// boxes (QTT,QXT,QXX,...,QYX) aligned in X-direction
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: &
         NQ                  !# of boxes in init. distr. to be combined & lim
    real(r8), intent(in) :: &
         ETT, &                ! tracer mass in extended zone
         EXT,EXX,EXY,EXZ,EYZ,EYT,EYY,EZT,EZZ, & ! moments in extended zone
         EM                   ! dry air mass in extended zone
    !// Output
    real(r8), dimension(NQ), intent(out) :: &
         QTT, &            ! tracer mass in box [I]
         QM, &             ! dry air mass in box [I]
         QXT,QXX, &     ! 1st & 2nd moments of tracer in X-direct
         QXY,QXZ, &     ! coupled 1st moments including X component
         QYZ,QYT,QYY, &
         QZT,QZZ       ! other SOMs without X

    !// Locals
    real(r8) :: &
         ETT0,EXT0,EXX0,EXY0,EXZ0,EYZ0,EYT0,EYY0,EZT0,EZZ0,EM0, &
         FLTIQ,FLTIQ2,ZNQ,ZNQ2
    integer :: N1N2,IQ
    !// --------------------------------------------------------------------

    !// Initialize
    ZNQ  = 1._r8 / real(NQ, r8)
    ZNQ2 = ZNQ * ZNQ
    N1N2 = (NQ+1) * (NQ+2)
    if (ETT .lt. 1.e-30_r8)  then
       EM0  = EM * ZNQ
       ETT0 = ETT * ZNQ
       !// remap to individual boxes in EPZ
       do IQ = 1, NQ
          QTT(IQ) = ETT0
          QXT(IQ) = 0._r8
          QXX(IQ) = 0._r8
          QYT(IQ) = 0._r8
          QZT(IQ) = 0._r8
          QXY(IQ) = 0._r8
          QXZ(IQ) = 0._r8
          QYY(IQ) = 0._r8
          QZZ(IQ) = 0._r8
          QYZ(IQ) = 0._r8
          QM(IQ)  = EM0
       end do
       return
    end if
    !// N.B. have not scaled Y,Z moments (EYT,...) to ETT, because not
    !// reversible
    !// must call positive limiter to ensure So > 0 in sub-boxes
    !Call QLIMIT (ETT,EXT,EXX,EXY,EXZ,EM,FMIN,FMAX,LIM)
    call QLIMIT2 (ETT,EXT,EXX,EXY,EXZ)

    !// average
    ETT0 = ETT*ZNQ
    EXT0 = EXT*ZNQ2
    EXX0 = EXX*ZNQ2*ZNQ
    EXY0 = EXY*ZNQ2
    EXZ0 = EXZ*ZNQ2
    EYT0 = EYT*ZNQ
    EZT0 = EZT*ZNQ
    EYY0 = EYY*ZNQ
    EZZ0 = EZZ*ZNQ
    EYZ0 = EYZ*ZNQ
    EM0  = EM*ZNQ
    !// remap to individual boxes in EPZ
    do IQ = 1, NQ
       FLTIQ   = real(2*IQ-NQ-1, r8)
       FLTIQ2  = real(6*IQ*(IQ-NQ-1)+N1N2, r8)
       QTT(IQ) = ETT0 + FLTIQ*EXT0 + FLTIQ2*EXX0
       QXT(IQ) = EXT0 + 3._r8*FLTIQ*EXX0
       QXX(IQ) = EXX0
       QYT(IQ) = EYT0 + FLTIQ*EXY0
       QZT(IQ) = EZT0 + FLTIQ*EXZ0
       QXY(IQ) = EXY0
       QXZ(IQ) = EXZ0
       QYY(IQ) = EYY0
       QZZ(IQ) = EZZ0
       QYZ(IQ) = EYZ0
       QM(IQ)  = EM0
    end do

    !// --------------------------------------------------------------------
  end subroutine QZONX
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module pbl_mixing
!//=========================================================================

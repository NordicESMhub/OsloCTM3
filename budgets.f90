!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// Budgets: Tendency budget diagnostics subroutines.
!//=========================================================================
module budgets
  !//-----------------------------------------------------------------------
  !// MODULE: budgets
  !// DESCRIPTION: Tendency budget diagnostics subroutines.
  !//
  !// Contains
  !//   subroutine TBGT_G (STTTNL,BTTTN0, NTDIAG)
  !//   subroutine TBGT_L (STTTNL,N1,N2,I1,I2,J1,J2,L1,L2,NTDIAG)
  !//   subroutine TBGT_IJ (BTT,BTTBCK,BTTTND,BTTTN0,MP,N1,N2, NTDIAG)
  !//   subroutine TBGT_P0(NDAY)
  !//   subroutine TBGT_P1(NDAY)
  !//   subroutine TBGT_P2(NDAY)
  !//
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------
  public
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine TBGT_G (STTTNL,BTTTN0, NTDIAG)
    !//---------------------------------------------------------------------
    !// Calculate/accumulate global(G) 3-D tendency budgets for
    !// process = NTDIAG. Carried out each NSUB.
    !//
    !// Resets all tendencies globally with NTDIAG=0.
    !//
    !// BTTTN0 and STTTNL are always reset! This is because they need to be
    !// zero for next NSUB-loop.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rTnd
    use cmn_size, only: LPAR, NPAR, MPBLK
    use cmn_ctm, only: NTM, GMTAU, JDAY, JYEAR, JMON, JDATE, TMON, IDAY, &
         STT
    use cmn_diag, only: NTND, NTDPAR, STTBCK, STTTND, STTTN0, &
         NDAY0, TAU0, JDAY0, JMON0, JDATE0, TMON0, JYEAR0, &
         USTEP, VSTEP, WSTEP, MNADV
    use cmn_oslo, only: CONVWASHOUT, DIAGEMIS_IJ, emisTotalsOld
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NTDIAG  ! index to initialize/reset the tendency
    real(r8), intent(inout) :: STTTNL(LPAR,NPAR,NTDPAR)
    real(r8), intent(inout) :: BTTTN0(MPBLK,NPAR,NTDPAR)
    !// Locals
    integer :: I, J, K, L, M, N
    !//---------------------------------------------------------------------

    if (NTDIAG .lt. 1)  then
       !// Initialize/ re-zero the budget tendencies and the ref STT=STTBCK
       STTTND(:,:,:,:,:) = 0._rTnd
       STTTN0(:,:) = 0._r8
       TAU0  = GMTAU
       NDAY0 = IDAY
       JDAY0 = JDAY
       JYEAR0= JYEAR
       JMON0 = JMON
       JDATE0 = JDATE
       TMON0 = TMON
       USTEP(:) = 0
       VSTEP(:) = 0
       WSTEP(:) = 0
       MNADV    = 0

       STTBCK(:,:,:,:) = STT(:,:,:,:)

       !// Initialize CTM3 convective washout diagnostic
       CONVWASHOUT(:,:,:,:,:) = 0._r8

       !// Initialize CTM3 emission diagnostic
       DIAGEMIS_IJ(:,:,:,:,:) = 0._r8
       emisTotalsOld(:) = 0._r8

       !// Initialise BTTTN0 and STTTNL
       BTTTN0(:,:,:)  = 0._r8
       STTTNL(:,:,:)  = 0._r8

    else
       !// ---accumulate 3D changes in STTTN0 under NTDIAG>0
       do N = 1, NTM
         do K = 1, NTND
           do M = 1, MPBLK
             STTTN0(N,K) = STTTN0(N,K) + BTTTN0(M,N,K)
           end do
           do L = 1, LPAR
             STTTN0(N,K) = STTTN0(N,K) + STTTNL(L,N,K)
           end do
         end do
       end do

       !// Initialise BTTTN0 and STTTNL each NSUB
       BTTTN0(:,:,:)  = 0._r8
       STTTNL(:,:,:)  = 0._r8

    end if

    !//---------------------------------------------------------------------
  end subroutine TBGT_G
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine TBGT_L(STTTNL,N1,N2,I1,I2,J1,J2,L1,L2, NTDIAG)
    !//---------------------------------------------------------------------
    !// Calculate/accumulate global(G) 3-D tendency budgets for
    !// process = NTDIAG.
    !// Usually called for global (I1=1,I2=IM; J1=1,J2=JM, ...) but can do
    !// IJ-domain.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, NPAR
    use cmn_ctm, only: STT
    use cmn_diag, only: NTDPAR, STTBCK, STTTND
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) ::  N1, N2, I1, I2, J1, J2, L1, L2
    integer, intent(in) ::  NTDIAG  ! index to store tendencies
    !// Input/output
    real(r8), intent(inout) :: STTTNL(LPAR,NPAR,NTDPAR)
    !// Locals
    integer :: I, J, L, N
    real(r8)  :: SUM8
    !//---------------------------------------------------------------------

    if (NTDIAG .gt. 0)  then

       !// Accumulate 3D changes in STT (tendencies) in CTM operation
       !// under NTDIAG>0
       do N = N1,N2
         do L = L1,L2
           SUM8  = 0._r8
           do J = J1,J2
             do I = I1,I2
               STTTND(I,J,L,N,NTDIAG) = STTTND(I,J,L,N,NTDIAG) &
                    + real(STT(I,J,L,N) - STTBCK(I,J,L,N), r8)
               SUM8  = SUM8 + (STT(I,J,L,N) - STTBCK(I,J,L,N))
             end do
           end do
           STTTNL(L,N,NTDIAG) = STTTNL(L,N,NTDIAG) + SUM8
         end do
       end do
       do N = N1,N2
         do L = L1,L2
           do J = J1,J2
             do I = I1,I2
               STTBCK(I,J,L,N) = STT(I,J,L,N)
             end do
           end do
         end do
       end do

    end if

    !//---------------------------------------------------------------------
  end subroutine TBGT_L
  !//-----------------------------------------------------------------------

      
  !//-----------------------------------------------------------------------
  subroutine TBGT_IJ (BTT,BTTBCK,BTTTND,BTTTN0,MP,N1,N2, NTDIAG)
    !//---------------------------------------------------------------------
    !// Calculate/accumulate global(G) 3-D tendency budgets for
    !// process = NTDIAG. Works within OpenMP IJ-block for tracers
    !// N=N1:N2, does not do AIR (N=0)
    !// use TBGT_G to initialize/reset the tendency diagnostics.
    !// Must call MPBIND to move IJ-block BTTTND to global STTTND,
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rTnd
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, MPBLK
    use cmn_diag, only: NTDPAR, NTDIAG
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: BTT(LPAR,NPAR,IDBLK,JDBLK)
    integer, intent(in)  :: MP      ! OpenMP block index
    integer, intent(in)  :: N1, N2  ! tracer index
    integer, intent(in)  :: NTDIAG  ! index to store tendencies
    !// Input/output
    real(r8), intent(inout) ::  BTTBCK(LPAR,NPAR,IDBLK,JDBLK)
    real(r8), intent(inout) ::  BTTTN0(MPBLK,NPAR,NTDPAR)
    real(rTnd), intent(inout) ::  BTTTND(LPAR,NPAR,IDBLK,JDBLK,NTDPAR)

    !// Locals
    real(r8) :: BTND_T(NPAR)
    integer ::  I, J, L, N
    !//---------------------------------------------------------------------

    if (NTDIAG > 0) then
       BTND_T(:)  = 0._r8
       do J = 1,JDBLK
         do I = 1,IDBLK
           do N = N1,N2
             do L = 1,LPAR
               BTTTND(L,N,I,J,NTDIAG) = BTTTND(L,N,I,J,NTDIAG) &
                    + real(BTT(L,N,I,J) - BTTBCK(L,N,I,J),r8)
               BTND_T(N) = BTND_T(N) + (BTT(L,N,I,J) - BTTBCK(L,N,I,J))
             end do
           end do
         end do
       end do

       do N = N1,N2
         BTTTN0(MP,N,NTDIAG) = BTTTN0(MP,N,NTDIAG) + BTND_T(N)
       end do

       do J = 1,JDBLK
         do I = 1,IDBLK
           do N = N1,N2
             do L = 1,LPAR
               BTTBCK(L,N,I,J) = BTT(L,N,I,J)
             end do
           end do
         end do
       end do
    end if

    !//---------------------------------------------------------------------
  end subroutine TBGT_IJ
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine TBGT_P0(NDAY)
    !//---------------------------------------------------------------------
    !// Processes/prints 0-D tendency budgets for globe, stddout only
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: IPAR, JPAR, LPAR, IMDIV, NTDPAR
    use cmn_ctm, only: AIR, STT, NTM, GMTAU, JDATE, TMON, JYEAR, &
         NRMETD, NROPSM, ANAME
    use cmn_chem, only: TNAME
    use cmn_diag, only: NDAY0, TAU0, JDATE0, TMON0, JYEAR0, &
         NTND, LBGA1, NBOXD, USTEP, VSTEP, WSTEP, MNADV, &
         TLDIAG, STTTN0
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY
    !// Locals
    real(r8) :: DT_HR, D_DAY, ANADV,AUSTEP,AVSTEP,AWSTEP
    integer :: I, J, L, M, M1, M2, N
    character(len=80) :: BTITLE
    !// Will also write out sum of tendecies ('>SUMS<') and
    !// instant mass ('INST-M')
    real(r8) :: BGTBX(NTDPAR+2)
    character(len=6) :: TLOUTPUT(NTDPAR+2)
    !//---------------------------------------------------------------------

    D_DAY = real(NDAY - NDAY0, r8)
    DT_HR = 24._r8 * D_DAY
    write(BTITLE(1:25),'(a25)') 'Global Tendency Bdgts(kg)'
    write(BTITLE(26:45),'(f8.3,i3,1x,a3,i5)') TAU0,JDATE0,TMON0,JYEAR0
    write(BTITLE(46:50),'(a5)') ' -to-'
    write(BTITLE(51:70),'(f8.3,i3,1x,a3,i5)') GMTAU,JDATE,TMON,JYEAR
    write(BTITLE(71:80),'(f10.2)') DT_HR

    M1 = NTND + 1
    M2 = NTND + 2
    TLOUTPUT(1:NTDPAR) = TLDIAG(:)
    TLOUTPUT(M1) = '>SUMS<'
    TLOUTPUT(M2) = 'INST-M'

    !// Calculate AIR mass (no tendencies are recorded)
    BGTBX(:)    = 0._r8
    do L = 1,LPAR
      do J = 1,JPAR
        do I = 1,IPAR
          BGTBX(M2) = BGTBX(M2) + AIR(I,J,L)
        end do
      end do
    end do

    ANADV = real(MNADV, r8) / real(NRMETD*D_DAY, r8)
    AUSTEP = real(sum(USTEP), r8) &
             / real(D_DAY*NRMETD*NROPSM*LPAR*NTM*(JPAR-2), r8)
    AVSTEP = real(sum(VSTEP), r8) &
             / real(D_DAY*NRMETD*NROPSM*LPAR*NTM*IPAR/2, r8)
    AWSTEP = real(sum(WSTEP), r8) &
             / real(D_DAY*NRMETD*NROPSM*IPAR*JPAR*NTM/IMDIV, r8)
    N=0
    write(6,'(a80)') BTITLE
    write(6,'(A,F6.2,A,3F7.3)') ' Averaged(per operator-split) NADV=', &
         ANADV,'  U, V & W steps=',AUSTEP,AVSTEP,AWSTEP
    write(6,'(a12,15(3x,a9))') '  tendencies', (TLOUTPUT(M),M=1,M2)
    write(6,'(i3,1x,a8,1p,15e12.4)') N,ANAME   ,(BGTBX(M),M=1,M2)

    !// ---all tracers
    do N = 1, NTM
       BGTBX(:)    = 0._r8
       !// Sum of all tendencies
       do M = 1,NTND
          BGTBX(M1) = BGTBX(M1) + STTTN0(N,M)
       end do
       !// Total tracer mass
       do L = 1,LPAR
         do J = 1,JPAR
           do I = 1,IPAR
             BGTBX(M2) = BGTBX(M2) + STT(I,J,L,N)
           end do
         end do
       end do

        write(6,'(i3,1x,a8,1p,15e12.4)') N,TNAME(N), &
            (STTTN0(N,M),M=1,NTND),(BGTBX(M),M=M1,M2)

    end do
        
    !//---------------------------------------------------------------------
  end subroutine TBGT_P0
  !//-----------------------------------------------------------------------


  !//-----------------------------------------------------------------------
  subroutine TBGT_P1(NDAY)
    !//---------------------------------------------------------------------
    !// Processes/prints 0-D & 1-D tendency budgets for specified boxes
    !// unformatted dump to UNIT=24
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rTnd
    use cmn_size, only: IPAR, JPAR, LPAR, NTDPAR
    use cmn_ctm, only: AIR, STT, NTM, GMTAU,JDATE,TMON,JYEAR, ANAME, &
         TALT, TLNG, TLAT
    use cmn_chem, only: TNAME
    use cmn_diag, only: NDAY0, TAU0, JDATE0, TMON0, JYEAR0,&
         NTND, IBOXD, JBOXD, KBOXD, TBOXD, NBOXD, TLDIAG, &
         LBGT1, LBGA1, STTTN0, STTTND
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) ::  NDAY
    real(r8) :: DT_HR
    integer :: I,II,J,K,L,M,M1,M2,N, IB0,IB1,JB0,JB1,LB0,LB1
    character(len=80) :: BTITLE
    !// Will also write out sum of tendecies ('>SUMS<') and
    !// instant mass ('INST-M')
    real(rTnd) :: BGTBX(NTDPAR+2),BGTBXL(LPAR,NTDPAR+2)
    real(rTnd) :: BGTBXJ(JPAR,NTDPAR+2),BGTBXI(IPAR,NTDPAR+2)
    character(len=6) :: TLOUTPUT(NTDPAR+2)
    !//---------------------------------------------------------------------

    DT_HR = 24._r8 * real(NDAY - NDAY0, r8)
    write(BTITLE(1:25),'(a25)') '1-D Tendency Budgets (kg)'
    write(BTITLE(26:45),'(f8.3,i3,1x,a3,i5)') TAU0,JDATE0,TMON0,JYEAR0
    write(BTITLE(46:50),'(a5)') ' -to-'
    write(BTITLE(51:70),'(f8.3,i3,1x,a3,i5)') GMTAU,JDATE,TMON,JYEAR
    write(BTITLE(71:80),'(f10.2)') DT_HR

    M1 = NTND + 1
    M2 = NTND + 2
    TLOUTPUT(1:NTDPAR) = TLDIAG(:)
    TLOUTPUT(M1) = '>SUMS<'
    TLOUTPUT(M2) = 'INST-M'

    !// Calculate 0-D and 1-D budgets (if .true.) for all boxes for AIR (N=0)
    if (LBGA1) then
       N=0
       write(6,'(a10,a80)') ANAME, BTITLE
       do K = 1,NBOXD
         IB0 = IBOXD(1,K)
         IB1 = IBOXD(2,K)
         JB0 = JBOXD(1,K)
         JB1 = JBOXD(2,K)
         LB0 = KBOXD(1,K)
         LB1 = KBOXD(2,K)
         BGTBX(:)    = 0._rTnd
         BGTBXL(:,:) = 0._rTnd
         BGTBXJ(:,:) = 0._rTnd
         BGTBXI(:,:) = 0._rTnd

         do L = LB0,LB1
           do J = JB0,JB1
             do II = IB0,IB1
               I = mod(II-1,IPAR) + 1
               BGTBX(   M2) = BGTBX(   M2) + AIR(I,J,L)
               BGTBXL(L,M2) = BGTBXL(L,M2) + AIR(I,J,L)
               BGTBXJ(J,M2) = BGTBXJ(J,M2) + AIR(I,J,L)
               BGTBXI(I,M2) = BGTBXI(I,M2) + AIR(I,J,L)
            end do
          end do
        end do

        !// Print/dump the diagnostics for AIR (N=0) for each box K
        write(6,'(a4,i3,2x,a15)') ' box',K, TBOXD(K)
        write(6,'(2i3,10(3x,a9))') N,K,TLOUTPUT(M2)
        write(6,'(a6,1p,10e12.5)') '  0-D ',BGTBX(M2)
        do L = LPAR, 1, -1
          write(6,'(a2,a4,1p,15e12.4)') 'L ',TALT(L),BGTBXL(L,M2)
        end do
        do J = JPAR, 1, -1
          write(6,'(a2,a4,1p,15e12.4)') 'J ',TLAT(J),BGTBXJ(J,M2)
        end do
        do I = 1, IPAR
          write(6,'(a2,a4,1p,15e12.4)') 'I ',TLNG(I),BGTBXI(I,M2)
        end do

      end do
    end if

    !// Calculate 0-D and 1-D budgets (if .true.) for all boxes for tracer N
    do N = 1, NTM
       if (LBGT1(N)) then
         write(6,'(a10,a80)') TNAME(N), BTITLE
         do K = 1,NBOXD
           IB0 = IBOXD(1,K)
           IB1 = IBOXD(2,K)
           JB0 = JBOXD(1,K)
           JB1 = JBOXD(2,K)
           LB0 = KBOXD(1,K)
           LB1 = KBOXD(2,K)
           BGTBX(:)    = 0._rTnd
           BGTBXL(:,:) = 0._rTnd
           BGTBXJ(:,:) = 0._rTnd
           BGTBXI(:,:) = 0._rTnd

           !// Total tracer mass (0-D & 1-D) in box
           do L = LB0,LB1
              do J = JB0,JB1
                 do II = IB0,IB1
                    I = mod(II-1,IPAR) + 1
                    BGTBX(   M2) = BGTBX(   M2) + STT(I,J,L,N)
                    BGTBXL(L,M2) = BGTBXL(L,M2) + STT(I,J,L,N)
                    BGTBXJ(J,M2) = BGTBXJ(J,M2) + STT(I,J,L,N)
                    BGTBXI(I,M2) = BGTBXI(I,M2) + STT(I,J,L,N)
                 end do
              end do
           end do

           !// Load all tendencies 1:NTND (0-D & 1-D) for the box
           do M = 1,NTND
              BGTBX(   M) = real(STTTN0(N,M),rTnd)
              do L = LB0,LB1
                do J = JB0,JB1
                  do II = IB0,IB1
                    I = mod(II-1,IPAR) + 1
                    BGTBXL(L,M) = BGTBXL(L,M) + STTTND(I,J,L,N,M)
                    BGTBXJ(J,M) = BGTBXJ(J,M) + STTTND(I,J,L,N,M)
                    BGTBXI(I,M) = BGTBXI(I,M) + STTTND(I,J,L,N,M)
                  end do
                end do
              end do
           end do

           !// Sum of all tendencies
           do M = 1,NTND
              BGTBX(   M1) = BGTBX(   M1) + BGTBX(   M)
              do L = 1,LPAR
                 BGTBXL(L,M1) = BGTBXL(L,M1) + BGTBXL(L,M)
              end do
              do J = 1,JPAR                         
                 BGTBXJ(J,M1) = BGTBXJ(J,M1) + BGTBXJ(J,M)
              end do
              do I = 1,IPAR
                 BGTBXI(I,M1) = BGTBXI(I,M1) + BGTBXI(I,M)
              end do
           end do

           !// Print/dump the diagnostics for tracer N for each box K
           write(6,'(a4,i3,2x,a15)') ' box',K, TBOXD(K)
           write(6,'(2i3,15(3x,a9))') N,K, (TLOUTPUT(M),M=1,M2)
           write(6,'(a6, 1p,15e12.5)') '  0-D ',(BGTBX(M),M=1,M2)
           do L = LB1, LB0, -1
              write(6,'(a2,a4,1p,15e12.4)')'L ',TALT(L),(BGTBXL(L,M),M=1,M2)
           end do
           write(6,'(2i3,15(3x,a9))') N,K, (TLOUTPUT(M),M=1,M2)
           do J = JB1, JB0, -1
              write(6,'(a2,a4,1p,15e12.4)')'J ',TLAT(J),(BGTBXJ(J,M),M=1,M2)
           end do
           write(6,'(2i3,15(3x,a9))') N,K, (TLOUTPUT(M),M=1,M2)
           do II=IB0,IB1
              I = mod(II-1,IPAR) + 1
              write(6,'(a2,a4,1p,15e12.4)')'I ',TLNG(I),(BGTBXI(I,M),M=1,M2)
           end do

         end do
       end if
    end do
        
    !//---------------------------------------------------------------------
  end subroutine TBGT_P1
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine TBGT_P2(NDAY)
    !// Processes/dumps(unf only)  2-D tendency budgets for specified boxes
    !// unformatted dump to UNIT=24
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rTnd
    use cmn_size, only: IPAR, JPAR, LPAR, NTDPAR
    use cmn_ctm, only: AIR, STT, NTM, GMTAU, JDATE, TMON, JYEAR, ANAME
    use cmn_chem, only: TNAME
    use cmn_diag, only: NDAY0, TAU0, JDATE0, TMON0, JYEAR0,&
         NTND, IBOXD, JBOXD, KBOXD, TBOXD, NBOXD, TLDIAG, &
         LBGT2, LBGA2, STTTND
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    integer, intent(in) ::  NDAY
    real(r8) :: DT_HR
    integer :: I,II,J,K,L,M,M1,M2,N, IB0,IB1,JB0,JB1,LB0,LB1
    character(len=80) :: BTITLE
    !// Will also write out sum of tendecies ('>SUMS<') and
    !// instant mass ('INST-M')
    real(rTnd) :: BGTBXJL(JPAR,LPAR,NTDPAR+2), &
                  BGTBXIJ(IPAR,JPAR,NTDPAR+2), &
                  BGTBXIL(IPAR,LPAR,NTDPAR+2)
    character(len=6) :: TLOUTPUT(NTDPAR+2)
    !//---------------------------------------------------------------------

    DT_HR = 24._r8 * real(NDAY - NDAY0, r8)
    write(BTITLE(1:25),'(a25)') '2-D Tendency Budgets (kg)'
    write(BTITLE(26:45),'(f8.3,i3,1x,a3,i5)') TAU0,JDATE0,TMON0,JYEAR0
    write(BTITLE(46:50),'(a5)') ' -to-'
    write(BTITLE(51:70),'(f8.3,i3,1x,a3,i5)') GMTAU,JDATE,TMON,JYEAR
    write(BTITLE(71:80),'(f10.2)') DT_HR

    M1 = NTND+1
    M2 = NTND+2
    TLOUTPUT(1:NTDPAR) = TLDIAG(:)
    TLOUTPUT(M1) = '>SUMS<'
    TLOUTPUT(M2) = 'INST-M'

    !// AIR (N=0) calculate 2-D budget tendencies = ONLY total inst mass
    if (LBGA2) then
      N=0
      write(6,'(a10,a80)') ANAME, BTITLE
      do K = 1,NBOXD
        IB0 = IBOXD(1,K)
        IB1 = IBOXD(2,K)
        JB0 = JBOXD(1,K)
        JB1 = JBOXD(2,K)
        LB0 = KBOXD(1,K)
        LB1 = KBOXD(2,K)
        BGTBXJL(:,:,:) = 0._rTnd
        BGTBXIJ(:,:,:) = 0._rTnd
        BGTBXIL(:,:,:) = 0._rTnd

        do L = LB0,LB1
          do J = JB0,JB1
            do II = IB0,IB1
              I = mod(II-1,IPAR) + 1
              BGTBXJL(J,L,M2) = BGTBXJL(J,L,M2) + real(AIR(I,J,L),rTnd)
              BGTBXIJ(I,J,M2) = BGTBXIJ(I,J,M2) + real(AIR(I,J,L),rTnd)
              BGTBXIL(I,L,M2) = BGTBXIL(I,L,M2) + real(AIR(I,J,L),rTnd)
            end do
          end do
        end do

        !// note for AIR (N=0) the unf write includes all processes
        !//but =0 only M=NTND+2
        !//write(24) N,K,BTITLE,TBOXD(K) 
        !//write(24) BGTBXJL,BGTBXIJ,BGTBXIL

      end do
    end if

    !// Calculate 2-D tendency budgets (if .true.) for all boxes for tracer N
    do N = 1, NTM
      if (LBGT2(N)) then
        write(6,'(a10,a80)') TNAME(N), BTITLE
        do K = 1,NBOXD
          IB0 = IBOXD(1,K)
          IB1 = IBOXD(2,K)
          JB0 = JBOXD(1,K)
          JB1 = JBOXD(2,K)
          LB0 = KBOXD(1,K)
          LB1 = KBOXD(2,K)
          BGTBXJL(:,:,:) = 0._rTnd
          BGTBXIJ(:,:,:) = 0._rTnd
          BGTBXIL(:,:,:) = 0._rTnd

          !// Total instantaneous tracer mass in box (M2)
          do L = LB0,LB1
            do J = JB0,JB1
              do II = IB0,IB1
                I = mod(II-1,IPAR) + 1
                BGTBXJL(J,L,M2) = BGTBXJL(J,L,M2) + real(STT(I,J,L,N),rTnd)
                BGTBXIJ(I,J,M2) = BGTBXIJ(I,J,M2) + real(STT(I,J,L,N),rTnd)
                BGTBXIL(I,L,M2) = BGTBXIL(I,L,M2) + real(STT(I,J,L,N),rTnd)
              end do
            end do
          end do

          !// cccumulate all 2-D tendencies 1:NTND for the box
          do M = 1, NTND
            do L = LB0,LB1
              do J = JB0,JB1
                do II = IB0,IB1
                  I = mod(II-1,IPAR) + 1
                  BGTBXJL(J,L,M) = BGTBXJL(J,L,M) + STTTND(I,J,L,N,M)
                  BGTBXIJ(I,J,M) = BGTBXIJ(I,J,M) + STTTND(I,J,L,N,M)
                  BGTBXIL(I,L,M) = BGTBXIL(I,L,M) + STTTND(I,J,L,N,M)
                end do
              end do
            end do
          end do

          !// Sum of all tendencies (M1)
          do M = 1, NTND
            do L = 1,LPAR
              do J = 1,JPAR
                BGTBXJL(J,L,M1) = BGTBXJL(J,L,M1) + BGTBXJL(J,L,M)
              end do
            end do
            do J = 1,JPAR                           
              do I = 1,IPAR
                BGTBXIJ(I,J,M1) = BGTBXIJ(I,J,M1) + BGTBXIJ(I,J,M)
              end do
            end do
            do L = 1,LPAR
              do I = 1,IPAR
                BGTBXIL(I,L,M1) = BGTBXIL(I,L,M1) + BGTBXIL(I,L,M)
              end do
            end do
          end do

          !//unformatted dump diagnostics for tracer N for each box K
          !//write(24) N,K,BTITLE,TBOXD(K) 
          !//write(24) BGTBXJL,BGTBXIJ,BGTBXIL

        end do
      end if
    end do
    !//---------------------------------------------------------------------
  end subroutine TBGT_P2
  !//-----------------------------------------------------------------------

  !//-----------------------------------------------------------------------
end module budgets
!//-------------------------------------------------------------------------

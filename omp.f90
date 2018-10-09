!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, March 2015
!//=========================================================================
!// OpenMP related routines.
!//=========================================================================
module omp
  !//-----------------------------------------------------------------------
  !// MODULE: omp
  !// DESCRIPTION: Routnes for OpenMP.
  !//
  !// Contains:
  !//   - subroutine MPSPLIT
  !//   - subroutine MPBIND
  !//   - subroutine get_iijjmp
  !//   - subroutine get_all_mpind
  !//-----------------------------------------------------------------------
  implicit none
  !//-----------------------------------------------------------------------

contains

  !//-----------------------------------------------------------------------
  subroutine MPSPLIT(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                     BTTTND,BTTBCK,AIRB,BTEM,GAMAB,LSTRATAIR_E90B,MP)
    !//---------------------------------------------------------------------
    !// Splits global arrays into private B-arrays.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rMom, rTnd
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW, &
         GAMA
    use cmn_chem, only: LSTRATAIR_E90
    use cmn_diag, only: NTDPAR
    use cmn_met, only: T
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    real(r8), intent(out), dimension(LPAR,IDBLK,JDBLK) :: AIRB, BTEM, GAMAB
    real(r8), intent(out), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT, BTTBCK
    real(rTnd), intent(out), dimension(LPAR,NPAR,IDBLK,JDBLK,NTDPAR) :: &
         BTTTND
    real(rMom), intent(out), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    logical, intent(out), dimension(LPAR,IDBLK,JDBLK) :: LSTRATAIR_E90B
    integer, intent(in) :: MP
    integer :: N, L, I, J, II, JJ, I1, I2, J1, J2
    !//---------------------------------------------------------------------

    J1 = MPBLKJB(MP)
    J2 = MPBLKJE(MP)
    I1 = MPBLKIB(MP)
    I2 = MPBLKIE(MP)

    do L = 1, LPAR
      do J = J1, J2
        JJ   = J - J1 + 1
        do I = I1, I2
          II = I - I1 + 1
          GAMAB(L,II,JJ) = GAMA(I,J,L)
          AIRB(L,II,JJ) = AIR(I,J,L)
          BTEM(L,II,JJ) = T(I,J,L)
        end do
      end do
    end do

    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ   = J - MPBLKJB(MP) + 1
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II = I - MPBLKIB(MP) + 1
        do L = 1, LPAR
          LSTRATAIR_E90B(L,II,JJ) = LSTRATAIR_E90(L,I,J)
        end do
      end do
    end do

    do N = 1, NPAR
      do L = 1, LPAR
        do J = J1, J2
          JJ   = J - J1 + 1
          do I = I1, I2
            II = I - I1 + 1
            BTT(L,N,II,JJ) = STT(I,J,L,N)
            BTTBCK(L,N,II,JJ) = STT(I,J,L,N)
            BXT(L,N,II,JJ) = SUT(I,J,L,N)
            BXX(L,N,II,JJ) = SUU(I,J,L,N)
            BYT(L,N,II,JJ) = SVT(I,J,L,N)
            BYY(L,N,II,JJ) = SVV(I,J,L,N)
            BZT(L,N,II,JJ) = SWT(I,J,L,N)
            BZZ(L,N,II,JJ) = SWW(I,J,L,N)
            BXY(L,N,II,JJ) = SUV(I,J,L,N)
            BXZ(L,N,II,JJ) = SUW(I,J,L,N)
            BYZ(L,N,II,JJ) = SVW(I,J,L,N)
          end do
        end do
      end do
    end do

    BTTTND(:,:,:,:,:) = 0._rTnd

    !//---------------------------------------------------------------------
  end subroutine MPSPLIT
  !//-----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
  subroutine MPBIND(BTT,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ, &
                    BTTTND,AIRB,LSTRATAIR_E90B,MP)
    !//---------------------------------------------------------------------
    !// Puts private B-arrays back into global arrays.
    !//---------------------------------------------------------------------
    use cmn_precision, only: r8, rMom, rTnd
    use cmn_size, only: LPAR, NPAR, NTDPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE, &
         AIR, STT, SUT, SVT, SWT, SUU, SVV, SWW, SUV, SUW, SVW
    use cmn_diag, only: NTND, STTBCK, STTTND
    use cmn_chem, only: LSTRATAIR_E90
    !//---------------------------------------------------------------------
    implicit none
    !//---------------------------------------------------------------------
    !// Input
    real(r8), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK) :: BTT
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: AIRB
    real(rTnd), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK,NTDPAR) :: &
         BTTTND
    real(rMom), intent(in), dimension(LPAR,NPAR,IDBLK,JDBLK) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    logical, intent(in), dimension(LPAR,IDBLK,JDBLK) :: LSTRATAIR_E90B
    integer, intent(in) :: MP
    integer :: N, L, I, J, II, JJ, M
    integer :: K, NM, JS, J1, J2, I1, I2, II1, II2
    !//---------------------------------------------------------------------

    J1 = MPBLKJB(MP)
    J2 = MPBLKJE(MP)
    I1 = MPBLKIB(MP)
    I2 = MPBLKIE(MP)
    do N = 1,NPAR
      do L = 1,LPAR
        do J = J1, J2
          JJ   = J - J1 + 1
          do I = I1, I2
            II = I - I1 + 1
            !// Assumes tendency cleared and logged to some process
            STTBCK(I,J,L,N) = BTT(L,N,II,JJ)
            STT(I,J,L,N) = BTT(L,N,II,JJ)
            SUT(I,J,L,N) = BXT(L,N,II,JJ)
            SUU(I,J,L,N) = BXX(L,N,II,JJ)
            SVT(I,J,L,N) = BYT(L,N,II,JJ)
            SVV(I,J,L,N) = BYY(L,N,II,JJ)
            SWT(I,J,L,N) = BZT(L,N,II,JJ)
            SWW(I,J,L,N) = BZZ(L,N,II,JJ)
            SUV(I,J,L,N) = BXY(L,N,II,JJ)
            SUW(I,J,L,N) = BXZ(L,N,II,JJ)
            SVW(I,J,L,N) = BYZ(L,N,II,JJ)
          end do
        end do
      end do
    end do

    !// M is outer index
    do M = 1, NTND
      do N = 1, NPAR
        do L = 1, LPAR
          do J = J1, J2
            JJ   = J - J1 + 1
            do I = I1, I2
              II = I - I1 + 1
              STTTND(I,J,L,N,M) = STTTND(I,J,L,N,M) + BTTTND(L,N,II,JJ,M) 
            end do
          end do
        end do
      end do
    end do


    do L = 1, LPAR
      do J = J1, J2
        JJ   = J - J1 + 1
        do I = I1, I2
          II = I - I1 + 1
          AIR(I,J,L) = AIRB(L,II,JJ)
        end do
      end do
    end do
    do J = J1, J2
      JJ   = J - J1 + 1
      do I = I1, I2
        II = I - I1 + 1
        do L = 1, LPAR
          LSTRATAIR_E90(L,I,J) = LSTRATAIR_E90B(L,II,JJ)
        end do
      end do
    end do

    !//---------------------------------------------------------------------
  end subroutine MPBIND
  !//-----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_iijjmp(iget,jget,ii,jj,mp)
    !// --------------------------------------------------------------------
    !// Send in global indices iget,jget, and get the MP-block indices.
    !//
    !// Ole Amund Sovde, February 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: MPBLKJB,MPBLKJE, MPBLKIB,MPBLKIE
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: iget, jget
    !// Output
    integer, intent(out) :: ii,jj,mp
    !// Locals
    integer :: I,J
    logical :: foundII, foundJJ
    !// ------------------------------------------------------------------

    !// Only necessary to initialize foundII once
    foundII = .false.

    do MP = 1, MPBLK

       !// Should initialize foundJJ each MP.
       foundJJ = .false.

       !// Check first for J. Then check if also I is in the current MP.
       !// When both are found, we are at the correct MP.
       do J = MPBLKJB(MP),MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          !// Found J
          if (J == jget) then
             foundJJ = .true.
             exit
          end if
       end do
       do I = MPBLKIB(MP),MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1
          !// Found I
          if (foundJJ .and. I == iget) then
             foundII = .true.
             exit
          end if
       end do

       !// When we have I and J, we have MP
       if (foundJJ .and. foundII) exit
    end do

    if (.not.foundJJ .or. .not.foundII) then
       !// This should never happen
       print*,'utilities_oslo.f90: get_iijjmp: Cannot find II/JJ',iget,jget
       print*,'II/JJ/MP:',ii,jj,mp
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine get_iijjmp
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_all_mpind()
    !// --------------------------------------------------------------------
    !// Find all MP-block indices as function of all i,j.
    !//
    !// Ole Amund Sovde, July 2011
    !// --------------------------------------------------------------------
    use cmn_size, only: MPBLK
    use cmn_ctm, only: MPBLKJB,MPBLKJE, MPBLKIB,MPBLKIE, &
         all_mp_indices
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Locals
    integer :: I,J,MP,II,JJ
    !// --------------------------------------------------------------------

    do MP = 1, MPBLK
       do J = MPBLKJB(MP), MPBLKJE(MP)
          JJ    = J - MPBLKJB(MP) + 1
          do I = MPBLKIB(MP), MPBLKIE(MP)
             II    = I - MPBLKIB(MP) + 1
             all_mp_indices(1,i,j) = II
             all_mp_indices(2,i,j) = JJ
             all_mp_indices(3,i,j) = MP
          end do
       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine get_all_mpind
  !// ----------------------------------------------------------------------



  !//-----------------------------------------------------------------------
end module omp
!//=========================================================================

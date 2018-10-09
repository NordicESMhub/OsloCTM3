!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// QSSA integrator for Oslo chemistry.
!//=========================================================================
module qssa_integrator
  !// ----------------------------------------------------------------------
  !// MODULE: qssa_integrator
  !// DESCRIPTION: Subroutines for Quasi Steady State Approximation (QSSA).
  !//
  !// Contains:
  !//   subroutine qssa
  !//   subroutine qssastr (allows for negative PROD)
  !//
  !// Ole Amund Sovde, September 2008
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine qssa(N,NAME,TIMESTEP,EULER,STEADYST,PROD,LOSS,CONC)
    !// --------------------------------------------------------------------
    !// This routine integrates the chemical term. The integration is
    !// based on the lifetime, or the Q-term, of the component.
    !//
    !// A large Q (compared with the inverse of the time step) implies
    !// a steady state (C = PROD/Q).
    !//
    !// A small Q allows a linear approximation of the exponential equation.
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)          :: N
    character(len=*), intent(in) :: NAME
    real(r8),intent(in)          :: TIMESTEP, EULER,STEADYST, PROD,LOSS
    !// Input/Output
    real(r8),intent(inout)       :: CONC

    !// Locals
    real(r8) :: RFAST, RELTAU
    logical  :: LONGLIVED, LSHTLIVED

    !// Parameters
    real(r8), parameter :: EPS_LOSS = 1.e-30_r8, EPS_PROD = 1.e-30_r8
    !// --------------------------------------------------------------------

    if (LOSS .gt. EPS_LOSS) then       !// If (LOSS.GT.0.d0) Then
       LONGLIVED = LOSS .lt. EULER
       LSHTLIVED = LOSS .gt. STEADYST
       !// RELTAU    = Min(1.d0,LOSS*TIMESTEP)
       RELTAU    = min(0.99999999_r8, LOSS*TIMESTEP)

       if (PROD .gt. EPS_PROD) then     !// If (PROD.GT.0.0) Then
          RFAST = PROD / LOSS
          if (LONGLIVED) then
            CONC = CONC + (RFAST - CONC) * RELTAU
          else if (LSHTLIVED) then
            CONC = RFAST
          else 
            CONC = RFAST + (CONC - RFAST) * exp(-RELTAU)
          endif
       else if (prod.lt.-EPS_PROD) then
          write(*,'(a,es12.5,i4,1x,a)') 'QSSA: PROD negative',PROD,N,NAME
       else
          !// PROD = 0. therefore expressions are somewhat different
          if (LONGLIVED) then
            CONC = CONC * (1._r8 - RELTAU)
          else
            CONC = CONC * exp(-RELTAU)
          end if 
       end if
    else if (LOSS .LT. -EPS_LOSS) then  !// ElseIf (LOSS.LT.0.) Then
       !// Negative loss should not happen
       write(6,*) 'WARNING LOSS .LT. 0.: Treated as zero'
       write(6,'(3E12.4)') LOSS,PROD,CONC
       write(6,'(A,I4,X,A)') 'Reaction number and name are ::',N,NAME
       CONC = CONC + PROD * TIMESTEP
    else
       CONC = CONC + PROD * TIMESTEP
    end if

    !// --------------------------------------------------------------------
  end subroutine qssa
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine qssastr (N,NAME,TIMESTEP,EULER,STEADYST,PROD,LOSS,CONC,I,J,L)
    !// --------------------------------------------------------------------
    !// This routine integrates the chemical term. The integration is
    !// based on the lifetime, or the Q-term, of the component.
    !//
    !// A large Q (compared with the inverse of the time step) implies
    !// a steady state (C = PROD/Q).
    !//
    !// A small Q allows a linear approximation of the exponential equation.
    !//
    !// Allows for negative PROD, which may happen in rare cases for families
    !// like sum of odd oxygen (SO).
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)          :: N, I, J, L
    character(len=*), intent(in) :: NAME
    real(r8),intent(in)          :: TIMESTEP, EULER,STEADYST, PROD,LOSS
    !// Input/Output
    real(r8),intent(inout)       :: CONC

    !// Locals
    real(r8) :: RFAST, RELTAU, OLDCONC
    logical :: LONGLIVED, LSHTLIVED

    !// Parameters
    real(r8), parameter :: EPS_LOSS = 1.e-30_r8
    !// --------------------------------------------------------------------

    OLDCONC = CONC

    if (LOSS .gt. EPS_LOSS) then
       LONGLIVED = LOSS .lt. EULER
       LSHTLIVED = LOSS .gt. STEADYST
       RELTAU    = Min(0.99999999_r8,LOSS*TIMESTEP)

       if (PROD .eq. 0._r8) then
          !// PROD = 0. therefore expressions are somewhat different
          if (LONGLIVED) then
            CONC = CONC * (1._r8 - RELTAU)
          else
            CONC = CONC * exp(-RELTAU)
          end if 

       else
          !// PROD .ne. 0.d0 In stratosphere, families may have negative PROD:
          RFAST = PROD / LOSS
          if (LONGLIVED) then
             !// LOSS.LT.EULER
             CONC = CONC + (RFAST - CONC) * RELTAU
          else if (LSHTLIVED) then
             !// LOSS.GT.STEADYST
             CONC = RFAST
          else
             !// LOSS.GE.EULER
             CONC = RFAST + (CONC - RFAST) * exp(-RELTAU)
          end if
       end if
    else if (LOSS .lt. -EPS_LOSS) then
       !// Negative LOSS should never occur, but for SO family it may happen
       !// when CONC is negative.
       if (CONC .gt. 0._r8) then
          print*,'QSSASTR: CONC > 0 and LOSS<0'
          print*,CONC,LOSS
          print*,'IJL: ',I,J,L
          stop
       end if
       write(6,'(A20,3E12.4)') 'WARNING LOSS .LT. 0:',LOSS,PROD,CONC
       write(6,'(A,I4,X,A)') 'Reaction number and name are ::',N,NAME
       write(6,'(a,3i5)') 'IJL: ',I,J,L
       !// Will try to use negative loss and concentration
       !// The integration for LOSS>0 can still be used, since it is LOSS*CONC
       !// that will be positive
       RFAST = PROD/LOSS ! LOSS negative, so this term is negative
       RELTAU    = Max(-0.99999999_r8, LOSS*TIMESTEP) !// This will be negative
       CONC = RFAST + (CONC - RFAST) * exp(-RELTAU) !// Exp-term >1.
       write(6,'(A29,E12.4)') 'WARNING LOSS .LT. 0: New CONC',CONC
    else
       !// Loss = zero
       CONC = CONC + PROD * TIMESTEP
    end If

    !// Integration done, check for negative
    if (CONC .lt. 0._r8) then
       print*,'QSSASTR: NEGATIVE CONC',N,NAME, CONC
       print*,'QSSASTR: OLD CONC',OLDCONC
       print*,'QSSASTR: TIMESTEP',TIMESTEP
       print*,'QSSASTR: PROD*TIMESTEP',PROD*TIMESTEP
       print*,'QSSASTR: LOSS',LOSS,I,J,L
    end if

    !// --------------------------------------------------------------------
  end subroutine qssastr
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module qssa_integrator
!//=========================================================================

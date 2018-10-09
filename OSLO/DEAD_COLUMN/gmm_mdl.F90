! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/gmm_mdl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Portable implementation of gamma functions
! Routines are external subroutines so that host program may
! link to machine's native erf functions if available
! fxm: Add incomplete gamma function gammq() to gmm_mdl

! Usage:
! use gmm_mdl ! [mdl] Gamma function gamma()

! Source: libspecfun.a specfun/src.dp/gamma.f
! http://scicomp.ewha.ac.kr/netlib/specfun/
! Modifications by Alf Grini and Charlie Zender
! 20030121: AG Turn into Fortran 90 free-format module
! 20030121: AG Use generic precision
! 20030219: CZ Clean up
! 20030219: CZ Inline use of real()


module gmm_mdl ! [mdl] Gamma function gamma()
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::gamma ! [fnc] Compute gamma function gamma(x)
! public::gammq ! [fnc] Compute incomplete gamma function gammq(x)
  
contains
  
  real(r8) function gamma(x) ! [fnc] Compute gamma function gamma(x)
    ! This routine calculates the GAMMA function for a real argument X.
    ! Computation is based on an algorithm outlined in reference 1.
    ! The program uses rational functions that approximate the GAMMA
    ! function to at least 20 significant decimal digits.  Coefficients
    ! for the approximation over the interval (1,2) are unpublished.
    ! Those for the approximation for X .GE. 12 are from reference 2.
    ! The accuracy achieved depends on the arithmetic system, the
    ! compiler, the intrinsic functions, and proper selection of the
    ! machine-dependent constants.
    
    ! Explanation of machine-dependent constants
    ! beta   - radix for the floating-point representation
    ! maxexp - the smallest positive power of beta that overflows
    ! XBIG   - the largest argument for which GAMMA(X) is representable
    !          in the machine, i.e., the solution to the equation
    !                  GAMMA(XBIG) = beta**maxexp
    ! XINF   - the largest machine representable floating-point number;
    !          approximately beta**maxexp
    ! EPS    - the smallest positive floating-point number such that
    !          1.0+EPS .GT. 1.0
    ! XMININ - the smallest positive floating-point number such that
    !          1/XMININ is machine representable
    
    !     Approximate values for some important machines are:
    !                            beta       maxexp        XBIG
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! Cyber 180/855
    !   under NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-Format   (D.P.)        2          127        34.844
    ! VAX G-Format   (D.P.)        2         1023        171.489
    !                            XINF         EPS        XMININ
    ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
    ! Cyber 180/855
    !   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
    ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
    ! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
    ! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
    
    ! Error returns
    ! The program returns the value XINF for singularities or
    ! when overflow would occur.  The computation is believed
    ! to be free of underflow and overflow.
    
    ! Intrinsic functions required are:
    ! INT, DBLE, EXP, LOG, REAL, SIN
    
    ! References: "An Overview of Software Development for Special
    ! Functions", W. J. Cody, Lecture Notes in Mathematics,
    ! 506, Numerical Analysis Dundee, 1975, G. A. Watson
    ! (ed.), Springer Verlag, Berlin, 1976.
    ! Computer Approximations, Hart, Et. Al., Wiley and
    ! sons, New York, 1968.
    ! Latest modification: October 12, 1989
    ! Authors: W. J. Cody and L. Stoltz
    ! Applied Mathematics Division
    ! Argonne National Laboratory
    ! Argonne, IL 60439
    integer I,N
    logical PARITY
    real(r8) &
         ! C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         C,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
    
    dimension C(7),P(8),Q(8)
    ! Mathematical constants
    data ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
         SQRTPI/0.9189385332046727417803297E0_r8/, &
         PI/3.1415926535897932384626434E0_r8/
    ! Machine dependent parameters
    data XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/,  &
         XINF/3.4E38_r8/
    ! Numerator and denominator coefficients for rational minimax
    ! approximation over (1,2).
    data P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,  &
         -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,    &
         8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,    &
         -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/  
    data Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8, & 
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
         2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8, &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
    ! Coefficients for minimax approximation over (12, INF).
    data C/-1.910444077728E-03_r8,8.4171387781295E-04_r8,                  &
         -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,              &
         -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8, &
         5.7083835261E-03_r8/
    ! csz++ Get rid of statement function by placing directly inline
    ! Statement functions for conversion between integer and float
    ! CONV(I) = REAL(I,kind=r8)
    ! CD    CONV(I) = DBLE(I)
    ! csz-- 
    PARITY = .FALSE.
    FACT = ONE
    N = 0
    Y = X
    IF (Y .LE. ZERO) THEN
       !  Argument is negative
       Y = -X
       Y1 = AINT(Y,kind=r8)
       RES = Y - Y1
       IF (RES .NE. ZERO) THEN
          IF (Y1 .NE. AINT(Y1*HALF,kind=r8)*TWO) PARITY = .TRUE.
          FACT = -PI / SIN(PI*RES)
          Y = Y + ONE
       ELSE
          RES = XINF
          GO TO 900
       END IF
    END IF
    !  Argument is positive
    IF (Y .LT. EPS) THEN
       !  Argument .LT. EPS
       IF (Y .GE. XMININ) THEN
          RES = ONE / Y
       ELSE
          RES = XINF
          GO TO 900
       END IF
    ELSE IF (Y .LT. TWELVE) THEN
       Y1 = Y
       IF (Y .LT. ONE) THEN
          !  0.0 .LT. argument .LT. 1.0
          Z = Y
          Y = Y + ONE
       ELSE
          !  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
          N = INT(Y) - 1
          Y = Y - REAL(N,kind=r8)
          Z = Y - ONE
       END IF
       !  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
       XNUM = ZERO
       XDEN = ONE
       DO 260 I = 1, 8
          XNUM = (XNUM + P(I)) * Z
          XDEN = XDEN * Z + Q(I)
260       CONTINUE
          RES = XNUM / XDEN + ONE
          IF (Y1 .LT. Y) THEN
             !  Adjust result for case  0.0 .LT. argument .LT. 1.0
             RES = RES / Y1
          ELSE IF (Y1 .GT. Y) THEN
             !  Adjust result for case  2.0 .LT. argument .LT. 12.0
             DO 290 I = 1, N
                RES = RES * Y
                Y = Y + ONE
290             CONTINUE
             END IF
          ELSE
             !  Evaluate for argument .GE. 12.0,
             IF (Y .LE. XBIG) THEN
                YSQ = Y * Y
                SUM = C(7)
                DO 350 I = 1, 6
                   SUM = SUM / YSQ + C(I)
350                CONTINUE
                   SUM = SUM/Y - Y + SQRTPI
                   SUM = SUM + (Y-HALF)*LOG(Y)
                   RES = EXP(SUM)
                ELSE
                   RES = XINF
                   GO TO 900
                END IF
             END IF
             !  Final adjustments and return
             IF (PARITY) RES = -RES
             IF (FACT .NE. ONE) RES = FACT / RES
900          GAMMA = RES
             !CD900 DGAMMA = RES
             RETURN
             
           end function gamma
           
         end module gmm_mdl ! [mdl] Gamma function gamma()

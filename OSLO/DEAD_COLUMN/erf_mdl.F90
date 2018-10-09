! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/erf_mdl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Portable error function family erf(), erfc(), erfcx()
! SGI and AIX machines both have problematic implementations of erf
! in some versions of their mathematical libraries. 

! Usage:
! use erf_mdl ! [mdl] Error functions erf(), erfc(), erfcx()

! Source: libspecfun.a specfun/src.dp/erf.f
! Modifications by Charlie Zender
! 1998xxxx: Add preprocessor directives for global models on NCAR supercomputers
! 20010927: support Fortran 90 precision statements, implicit none
! 20010928: Turn into Fortran 90 free-format module

module erf_mdl ! [mdl] Error functions erf(), erfc(), erfcx()
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  private::calerf ! [sbr] Internal routine for error function evaluation
  public::erf ! [fnc] Compute erf(x)
  public::erfc ! [fnc] Compute erfc(x)

contains
  
#ifndef CRAY /* Cray has a good erf() intrinsic, others do not */
  
  subroutine calerf(arg,result,jint) ! [sbr] Internal routine for error function evaluation
    ! This packet evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
    ! for a real argument  x.  It contains three function type
    ! subprograms: erf, erfc, and erfcx (or derf, derfc, and derfcx),
    ! and one subroutine type subprogram, calerf.  The calling
    ! statements for the primary entries are:
    
    ! y=erf(x)     (or   y=derf(x)),
    ! y=erfc(x)    (or   y=derfc(x)),
    ! and
    ! y=erfcx(x)   (or   y=derfcx(x)).
    
    ! The routine  calerf  is intended for internal packet use only,
    ! all computations within the packet being concentrated in this
    ! routine.  The function subprograms invoke  calerf  with the
    ! statement
    ! call calerf(arg,result,jint)
    ! where the parameter usage is as follows
    
    ! Function                     Parameters for calerf
    ! Call              Arg                  Result          Jint
    ! 
    ! erf(arg)      any real argument         erf(arg)          0
    ! erfc(arg)     abs(arg)  <  xbig        erfc(arg)          1
    ! erfcx(arg)    xneg  <  arg  <  xmax   erfcx(arg)          2
    
    ! The main computation evaluates near-minimax approximations:
    ! from "Rational Chebyshev Approximations for the Error Function"
    ! by W. J. Cody, Math. Comp., 1969, pp. 631-638.  This
    ! transportable program uses rational functions that theoretically
    ! approximate  erf(x)  and  erfc(x)  to at least 18 significant
    ! decimal digits.  The accuracy achieved depends on the arithmetic
    ! system, the compiler, the intrinsic functions, and proper
    ! selection of the machine-dependent constants.
    
    ! Explanation of machine-dependent constants:
    ! xmin   = The smallest positive floating-point number.
    ! xinf   = The largest positive finite floating-point number.
    ! xneg   = The largest negative argument acceptable to erfcx;
    ! the negative of the solution to the equation
    ! 2*exp(x*x) = xinf.
    ! xsmall = Argument below which erf(x) may be represented by
    ! 2*x/sqrt(pi)  and above which  x*x  will not underflow.
    ! A conservative value is the largest machine number x
    ! such that   1.0 + x = 1.0   to machine precision.
    ! xbig   = Largest argument acceptable to erfc;  solution to
    ! the equation:  w(x)* (1-0.5/x**2) = xmin,  where
    ! w(x) = exp(-x*x)/[x*sqrt(pi)].
    ! xhuge  = Argument above which  1.0 - 1/(2*x*x) = 1.0  to
    ! machine precision.  a conservative value is
    ! 1/[2*sqrt(xsmall)]
    ! xmax   = Largest acceptable argument to erfcx; the minimum
    ! of xinf and 1/[sqrt(pi)*xmin].
    
    ! Approximate values for some important machines are:
    !                       xmin       xinf        xneg     xsmall
    ! CDC 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
    ! Cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
    ! IEEE (IBM/XT,
    ! Sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
    ! IEEE (IBM/XT,
    ! Sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
    ! IBM 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
    ! Univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
    ! Vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
    ! Vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16
    
    !                        xbig       xhuge       xmax
    ! CDC 7600      (s.p.)  25.922      8.39e+6     1.80x+293
    ! Cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
    ! IEEE (IBM/XT,
    ! Sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
    ! IEEE (IBM/XT,
    ! Sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
    ! IBM 195       (d.p.)  13.306      1.90d+8     7.23e+75
    ! Univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
    ! Vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
    ! Vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307
    
    ! Error returns:
    ! The program returns  erfc = 0      for  arg  >=  xbig;
    ! erfcx = xinf  for  arg  <  xneg;
    ! and
    ! erfcx = 0     for  arg  >=  xmax.
    
    ! Intrinsic functions required are:
    ! abs, aint, exp
    
    ! Author: W. J. Cody
    ! Mathematics And Computer Science Division
    ! Argonne National Laboratory
    ! Argonne, IL 60439
    ! Latest modification: March 19, 1990
    implicit none
    integer i,jint
    real(r8) &
         a,arg,b,c,d,del,four,half,p,one,q,result,sixten,sqrpi, &
         two,thresh,x,xbig,xden,xhuge,xinf,xmax,xneg,xnum,xsmall, &
         y,ysq,zero
    dimension a(5),b(4),c(9),d(8),p(6),q(5)
    
    ! Mathematical constants
    data four,one,half,two,zero/4.0e0_r8,1.0e0_r8,0.5e0_r8,2.0e0_r8,0.0e0_r8/, &
         sqrpi/5.6418958354775628695e-1_r8/,thresh/0.46875e0_r8/, &
         sixten/16.0e0_r8/
    
    ! Machine-dependent constants
    data xinf,xneg,xsmall/3.40e+38_r8,-9.382e0_r8,5.96e-8_r8/, &
         xbig,xhuge,xmax/9.194e0_r8,2.90e3_r8,4.79e37_r8/
    
    ! Coefficients for approximation to  erf  in first interval
    data a/3.16112374387056560e00,1.13864154151050156e02, &
         3.77485237685302021e02,3.20937758913846947e03, &
         1.85777706184603153e-1/
    data b/2.36012909523441209e01,2.44024637934444173e02, &
         1.28261652607737228e03,2.84423683343917062e03/
    
    ! Coefficients for approximation to  erfc  in second interval
    data c/5.64188496988670089e-1,8.88314979438837594e0, &
         6.61191906371416295e01,2.98635138197400131e02, &
         8.81952221241769090e02,1.71204761263407058e03, &
         2.05107837782607147e03,1.23033935479799725e03, &
         2.15311535474403846e-8/
    data d/1.57449261107098347e01,1.17693950891312499e02, &
         5.37181101862009858e02,1.62138957456669019e03, &
         3.29079923573345963e03,4.36261909014324716e03, &
         3.43936767414372164e03,1.23033935480374942e03/
    
    ! Coefficients for approximation to  erfc  in third interval
    data p/3.05326634961232344e-1,3.60344899949804439e-1, &
         1.25781726111229246e-1,1.60837851487422766e-2, &
         6.58749161529837803e-4,1.63153871373020978e-2/
    data q/2.56852019228982242e00,1.87295284992346047e00, &
         5.27905102951428412e-1,6.05183413124413191e-2, &
         2.33520497626869185e-3/
    
    ! Main Code
    x=arg
    y=abs(x)
    if (y <= thresh) then
       ! Evaluate  erf  for  |x| <= 0.46875
       ysq=zero
       if (y > xsmall) ysq=y*y
       xnum=a(5)*ysq
       xden=ysq
       do i=1,3
          xnum=(xnum+a(i))*ysq
          xden=(xden+b(i))*ysq
       end do
       result=x*(xnum+a(4))/(xden+b(4))
       if (jint /= 0) result=one-result
       if (jint == 2) result=exp(ysq)*result
       go to 800
       ! Evaluate  erfc  for 0.46875 <= |x| <= 4.0
    else if (y <= four) then
       xnum=c(9)*y
       xden=y
       do i=1,7
          xnum=(xnum+c(i))*y
          xden=(xden+d(i))*y
       end do
       result=(xnum+c(8))/(xden+d(8))
       if (jint /= 2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=exp(-ysq*ysq)*exp(-del)*result
       end if
       ! Evaluate  erfc  for |x| > 4.0
    else
       result=zero
       if (y >= xbig) then
          if ((jint /= 2).or.(y >= xmax)) go to 300
          if (y >= xhuge) then
             result=sqrpi/y
             go to 300
          end if
       end if
       ysq=one/(y*y)
       xnum=p(6)*ysq
       xden=ysq
       do i=1,4
          xnum=(xnum+p(i))*ysq
          xden=(xden+q(i))*ysq
       end do
       result=ysq*(xnum+p(5))/(xden+q(5))
       result=(sqrpi-result)/y
       if (jint /= 2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=exp(-ysq*ysq)*exp(-del)*result
       end if
    end if
    ! Fix up for negative argument, erf, etc.
300 if (jint == 0) then
       result=(half-result)+half
       if (x < zero) result=-result
    else if (jint == 1) then
       if (x < zero) result=two-result
    else
       if (x < zero) then
          if (x < xneg) then
             result=xinf
          else
             ysq=aint(x*sixten)/sixten
             del=(x-ysq)*(x+ysq)
             y=exp(ysq*ysq)*exp(del)
             result=(y+y)-result
          end if
       end if
    end if
800 return
  end subroutine calerf
  
  function erf(x) ! [fnc] Compute erf(x)
    implicit none
    ! This subprogram computes approximate values for erf(x)
    ! (see comments heading calerf).
    ! Author/Date: W. J. Cody, January 8, 1985
    integer jint
    real(r8) x,result,erf
    jint=0
    call calerf(x,result,jint)
    erf=result
    return
  end function erf
  
  function erfc(x) ! [fnc] Compute approximate values for erfc(x)
    implicit none
    ! This subprogram computes approximate values for erfc(x).
    ! (see comments heading calerf).
    ! Author/Date: W. J. Cody, January 8, 1985
    integer jint
    real(r8) x,result,erfc
    jint=1
    call calerf(x,result,jint)
    erfc=result
    return
  end function erfc
  
  function erfcx(x)
    ! This subprogram computes approximate values for exp(x*x)*erfc(x)
    ! (see comments heading calerf)
    ! Author/Date: W. J. Cody, March 30, 1987
    integer jint
    real(r8) x,result,erfcx
    jint=2
    call calerf(x,result,jint)
    erfcx=result
    return
  end function erfcx
  
#else /* CRAY */
  
  ! Ensure file has something to compile
  ! Above routines may be superceded by any local erf implementation
  ! SGI and AIX machines have both been shown to have problematic 
  ! implementations of erf in some versions of their mathematical libraries.
  ! Other architectures may be fine and use intrinsic erf implementation.
  ! But, if user tells compiler to compile this code without
  ! anything here -- the compiler complains.  A more sophisticated
  ! build mechanism could decide which files should and should not
  ! be built, but the simple solution is to do as follows...
  subroutine erf_stb
    stop 'erf_stb() ERROR: This routine should not be called'
  end subroutine erf_stb
#endif  /* not CRAY */
  
end module erf_mdl

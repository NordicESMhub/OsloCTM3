! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/wbl_mdl.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Utilities for Weibull PDF of wind speeds
! References:
! K. G. Justus et al., J. Appl. Met., 17(3), 350-354 (JHM78)
! Author: Alf Grini, alf.grini@geofysikk.uio.no
! History:
! 20030120 A. Grini  Original Version
! 20030218 C. Zender Cleaned up
! Usage: 
! use wbl_mdl ! [mdl] Weibull wind speed distribution

module wbl_mdl ! [mdl] Weibull wind speed distribution
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  private::wbl_prm_rfr_get ! [sbr] Compute Weibull c & k parameters at reference height
  private::wbl_rfr2mdp ! [sbr] Move Weibull parameters c & k to midpoint height
  private::wnd_mdp_wgt_wbl_get ! [sbr] Determine weights of discrete wind bins
  private::wnd_min_max_wbl_get ! [sbr] Determine max, min discrete wind speed
  public::wbl_wnd ! [sbr] Weibull winds should be known to other modules
  public::wnd_mdp_bin_get ! [sbr] Weibull winds should be known to other modules
  
contains
  
  !// ------------------------------------------------------------------
  subroutine wbl_wnd( & ! [sbr] Generate Weibull wind speed PDF abscissae and weights
       hgt_mdp, & ! I [m] Height of layer midpoint
       hgt_rfr, & ! I [m] Reference height
       wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
       wnd_mdp, & ! I [m s-1] Wind speed at midpoint
       wnd_mdp_min, & ! O [m s-1] Min wind speed at midpoint
       wnd_mdp_max, & ! O [m s-1] Max wind speed at midpoint 
       wnd_mdp_nbr, & ! I [nbr] Number of wind speeds
       wnd_mdp_wgt, & ! O [frc] Wind speed bin weight
       wnd_rfr) ! I [m s-1] Wind speed at reference height
    !// ------------------------------------------------------------------
    ! Purpose: Generate Weibull wind speed PDF abscissae and weights
    ! wbl_wnd() is called by dst_mbl()
    ! Resolve wind PDF into wbl_wnd_nbr bins with "wbl_wnd_nbr+1" limits
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    real(r8),intent(in)::hgt_mdp ! [m] Midpoint height
    real(r8),intent(in)::hgt_rfr ! [m] Reference height
    real(r8),intent(in)::wnd_frc_rsl ! [frc] Fraction of wind PDF to resolve
    integer,intent(in)::wnd_mdp_nbr ! [nbr] Number of "bins" in which we split the wind distr.
    real(r8),intent(in)::wnd_rfr ! [m s-1] Wind speed at reference height
    real(r8),intent(in)::wnd_mdp ! [m s-1] Wind speed at midpoint
    
    ! Output
    real(r8),intent(out)::wnd_mdp_max ! [m s-1] Maximum wind speed in grid midpoint
    real(r8),intent(out)::wnd_mdp_min ! [m s-1] Minimum wind speed in grid midpoint
    real(r8),intent(out)::wnd_mdp_wgt(wnd_mdp_nbr) ! [frc] Wind speed bin weight
    
    ! Local
    real(r8)::c_wbl_mdp ! Weibull distribution scale factor (midpoint)
    real(r8)::c_wbl_rfr ! Weibull distribution scale factor (refence height)
    real(r8)::k_wbl_mdp ! Weibull distribution shape factor (midpoint)
    real(r8)::k_wbl_rfr ! Weibull distribution shape factor (reference height)
    !// ------------------------------------------------------------------
    
    ! Compute Weibull c & k parameters at reference height
    call wbl_prm_rfr_get( &
         c_wbl_rfr, & ! O [frc] Weibull scale factor at reference height
         k_wbl_rfr, & ! O [frc] Weibull shape factor at reference height
         wnd_rfr) ! I [m s-1] Wind speed at reference heightk_wbl_rfr 
    
    ! Move Weibull parameters from c & k from reference to midpoint height
    call wbl_rfr2mdp( &
         c_wbl_mdp, & ! O [frc] Weibull distribution scale factor (midpoint)
         c_wbl_rfr, & ! I [frc] Weibull distribution shape factor (reference height)
         hgt_mdp, & ! I [m] Midpoint height
         hgt_rfr, & ! I [m] Reference height
         k_wbl_mdp, & ! O [frc] Weibull distribution shape factor (midpoint)
         k_wbl_rfr, & ! I [frc] Weibull distribution shape factor (reference height)
         wnd_mdp) ! I [m s-1] wind at midpoint
    
    ! Determine max, min discrete wind speed at layer midpoint
    call wnd_min_max_wbl_get( &
         c_wbl_mdp, & ! I [frc] Weibull distribution scale factor (midpoint)
         k_wbl_mdp, & ! I [frc] Weibull distribution 
         wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
         wnd_mdp_max, & ! O [m s-1] Max wind speed to calculate
         wnd_mdp_min) ! O [m s-1] Min wind speed to calculate
    
    ! Determine weights of discrete wind bins
    call wnd_mdp_wgt_wbl_get( &
         c_wbl_mdp, & ! I [frc] Weibull distribution scale factor (midpoint)
         k_wbl_mdp, & ! I [frc] Weibull distribution shape factor (midpoint)
         wnd_mdp_max, & ! I [m s-1] max wind speed which we calculate
         wnd_mdp_min, & ! I [m s-1] min wind speed which we calculate
         wnd_mdp_nbr, & ! I [nbr] number for which we calculate wind speeds
         wnd_mdp_wgt) ! O [frc] weight given to each wind speed "bin"
    
    !// ------------------------------------------------------------------
  end subroutine wbl_wnd
  !// ------------------------------------------------------------------
  
  !// ------------------------------------------------------------------
  subroutine wbl_prm_rfr_get( & ! [sbr] Compute Weibull c & k parameters at reference height
       c_wbl_rfr, & ! O Weibull distribution scale factor
       k_wbl_rfr, & ! O Weibull distribution shape factor
       wnd_rfr) ! I Mean wind speed at reference height
    ! Purpose: Compute c & k of Weibull distribution winds at anemometer height (10 m)
    ! Author: Alf Grini (alf.grini@geofysikk.uio.no)
    !//
    !// Rewritten for Oslo CTM3, works on single values.
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    use gmm_mdl,only:gamma ! [mdl] Gamma function gamma()
    implicit none
    ! Input:
    real(r8),intent(in)::wnd_rfr ! I Reference height wind speed
    ! Output:
    real(r8),intent(out)::c_wbl_rfr ! O Weibull distribution scale factor
    real(r8),intent(out)::k_wbl_rfr ! O Weibull distribution shape factor 
    ! Local:
    real(r8),parameter::cst_k_wbl=0.94_r8 ! Constant in eq. 20 Justus for average variablility
    real(r8)::gmm_prm ! Parameter related to gamma distribution
    
    ! Get Weibull shape factor, assume average variablility in winds
    ! Justus, eq. 20
    k_wbl_rfr = cst_k_wbl*sqrt(wnd_rfr)
       
    ! write(6,*)'reference height wind',wnd_rfr(lon_idx)
    ! write(6,*)'k_wbl_rfr',k_wbl_rfr(lon_idx)
       
    ! Get gammafunction of (1+1/k) 
    gmm_prm = gamma((1.0_r8 + 1.0_r8/k_wbl_rfr))
    ! write(6,*)'gmm_prm',gmm_prm
       
    ! Get Weibull scale factor from mean wind speed
    ! Justus, eq. 16.
    c_wbl_rfr = wnd_rfr/gmm_prm
    ! write(6,*)'c_wbl_rfr',c_wbl_rfr
    
    ! write(6,*)'gamma 0.5',gamma(0.5_r8)
    ! write(6,*)'sqrt (pi)',sqrt(3.141592654_r8)
    
    ! do i=1,24
    !   xx=real(i,r8)/6.0_r8
    !   write(6,*)'gamma',xx,gamma(xx)
    ! enddo
    
  end subroutine wbl_prm_rfr_get
  !// ------------------------------------------------------------------
  
  !// ------------------------------------------------------------------
  subroutine wbl_rfr2mdp( &
       c_wbl_mdp, & ! O [frc] Weibull scale parameter (midpoint height)
       c_wbl_rfr, & ! I [frc] Weibull scale parameter (reference height)
       hgt_mdp, & ! I [m] midpoint height  
       hgt_rfr, & ! I [m] reference height
       k_wbl_mdp, & ! O [frc] Weibull shape parameter (midpoint height)
       k_wbl_rfr, & ! I [frc] Weibull shape parameter (reference height)
       wnd_mdp) ! I [m s-1] wind at layer midpoint
    ! Purpose: 
    ! Interpolate Weibull shape and scale factors from reference height (where they are calculated)
    ! to midpoint height (where they are used)
    ! Reference 
    ! K.G. Justus et. al. , J. appl. met, vol 17, no 3, pp 350-354
    ! Author: Alf Grini, alf.grini@geofysikk.uio.no
    !//
    !// Rewritten for Oslo CTM3, works on single values.
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    use gmm_mdl
    implicit none
    
    ! INPUT 
    real(r8),intent(in)::c_wbl_rfr ! I [frc] Weibull distribution scale factor (reference height)
    real(r8),intent(in)::hgt_mdp ! I [m] Height of layer midpoint
    real(r8),intent(in)::hgt_rfr ! I [m] Reference height
    real(r8),intent(in)::k_wbl_rfr ! I [frc] Weibull distribution shape factor (reference height)
    real(r8),intent(in)::wnd_mdp ! I [m s-1] wind speed at layer midpoint
    
    ! OUTPUT
    real(r8),intent(out)::c_wbl_mdp ! O [frc] Weibull distribution scale factor (midpoint)
    real(r8),intent(out)::k_wbl_mdp ! O [frc] Weibull distribution shale factor (midpoint)
    
    ! LOCAL
    real(r8)::n_wbl ! [-] Factor in Justus, eq. 9
    ! real(r8)::xx
    
    ! Equation 9 in the Justus paper
    ! n_wbl=(0.37_r8-0.088_r8*log(c_wbl_rfr(i))) &
    !     / (1-0.088*log(hgt_rfr/10.0_r8)) 
       
    ! Equation 7 in the Justus paper proposes
    ! c_wbl_mdp(i)=c_wbl_rfr(i)*(hgt_mdp(i)/hgt_rfr)**n_wbl
    ! where n_wbl is given by eq. 9 :
    ! n_wbl=(0.37_r8-0.088_r8*log(c_wbl_rfr(i))) &
    !     / (1-0.088*log(hgt_rfr/10.0_r8))
       
    ! The problem is over-determined since we also have the wind in the midpoint, and since we 
    ! assume the wind is Weibull distributed. Then we instead use the relation between mean and shape factor
       
    ! Equation 8 in the Justus paper adjusts the shape of the wind distribution to another height
    k_wbl_mdp = k_wbl_rfr*(1._r8 - 0.088_r8*log(hgt_rfr/10.0_r8)) &
         /(1._r8 - 0.088_r8*log(hgt_mdp/10.0_r8))
       
    ! Relation between mean (wnd_mdp) and shape factor (k) for Weibull distribution
    c_wbl_mdp = wnd_mdp/gamma(1.0_r8 + 1.0_r8/k_wbl_mdp)
    ! write(6,*)'c mdp /rfr ',c_wbl_mdp,c_wbl_rfr
    ! write(6,*)'k mdp /rfr ',k_wbl_mdp,k_wbl_rfr
    
  end subroutine wbl_rfr2mdp
  !// ------------------------------------------------------------------
  
  !// ------------------------------------------------------------------
  subroutine wnd_min_max_wbl_get( &
       c_wbl_mdp, & ! I [frc] Weibull distribution scale factor
       k_wbl_mdp, & ! I [frc] Weibull distribution shape factor
       wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
       wnd_mdp_max, & ! O [m s-1] maximum wind speed for which we do calc
       wnd_mdp_min) ! O [m s-1] minimum wind speed for which we do calc
    ! Purpose: 
    ! Find out what winds we want to calculate dust emissions for.
    ! If for example wnd_frc_rsl=0.95, the max wind we will use is the wind speed
    ! For which 95 % of the wind speeds is smaller and the lowest wind speed is 
    ! the one for which 5 % is lower. Wnd_frc_rsl thus has to be a number between
    ! zero and one.
    ! Author: Alf Grini, alf.grini@geofysikk.uio.no
    !//
    !// Rewritten for Oslo CTM3, works on single values.
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    real(r8),intent(in)::c_wbl_mdp ! Weibull distribution scale factor
    real(r8),intent(in)::k_wbl_mdp ! Weibull distribution shape factor
    real(r8),intent(in)::wnd_frc_rsl ! Fraction of wind PDF to resolve
    ! Output
    real(r8),intent(out)::wnd_mdp_max ! Highest wind speed for which we do calc
    real(r8),intent(out)::wnd_mdp_min ! Lowest wind speed for which we do calc
    ! Local
    real(r8)::Fmin ! [frc] Fraction of winds slower than wnd_mdp_min_wbl
    real(r8)::Fmax ! [frc] Fraction of winds slower than wnd_mdp_max_wbl
    real(r8)::xpn_fct_max ! The factor (x_max/c)**k in the cumulative distribution
    real(r8)::xpn_fct_min ! The factor (x_min/c)**k in the cumulative distribution

    ! Begin code
    Fmin = 1.0_r8-wnd_frc_rsl ! [frc] Fraction of winds slower than wnd_mdp_min_wbl 
    Fmax = wnd_frc_rsl ! [frc] Fraction of winds slower than wnd_mdp_max_wbl
    
    if(Fmin >= Fmax)then
       write(6,*)'ERROR in wnd_min_max_wbl_get(): '
       write(6,*)'Specify fraction larger than 0.5 to resolve'
       write(6,*)'Doing otherwise is no better than using mean wind speed'
       stop
    endif ! endif err

    ! Cumulative distribution with parameter c & k
    ! F(x)= 1-exp(-(x/c)**k) (where x is the wind speed we are interested in)
       
    ! Get the factor (x_min/c)**k
    xpn_fct_min = -1.0_r8*log(1.0_r8-Fmin) 
       
    ! Get the factor (x_max/c)**k
    xpn_fct_max = -1.0_r8*log(1.0_r8-Fmax)
       
    ! The minimum wind speed of our distribution
    wnd_mdp_min = &
         (xpn_fct_min)**(1.0_r8/k_wbl_mdp) &
         *c_wbl_mdp
       
    ! The maximum wind speed of our distribution
    wnd_mdp_max = &
         (xpn_fct_max)**(1.0_r8/k_wbl_mdp)*c_wbl_mdp
       
    ! write(6,*)'wnd_mdp_min ',wnd_mdp_min
    ! write(6,*)'wnd_mdp_max ',wnd_mdp_max

  end subroutine wnd_min_max_wbl_get
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine wnd_mdp_wgt_wbl_get( &
       c_wbl_mdp, & ! I [frc] Weibull distribution scale factor
       k_wbl_mdp, & ! I [frc] Weibull distribution shape factor
       wnd_mdp_max, & ! I [m s-1] max wind speed which we want to calculate
       wnd_mdp_min, & ! I [m s-1] min wind speed which we want to calculate
       wnd_mdp_nbr, & ! I [nbr] Number for which we calculate winds (number of wind "bins")
       wnd_mdp_wgt) ! O [frc] Wind speed bin weight
    ! Purpose: Given the Weibull distribution factors and the min and max wind speeds
    ! which we want to calculate, we get an array of weights for the wind speeds.
    ! Author: Alf Grini (alf.grini@geofysikk.uio.no)
    !//
    !// Rewritten for Oslo CTM3, works on single values.
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    real(r8),intent(in)::c_wbl_mdp ! I [frc] Weibull distribution scale factor (midpoint)
    real(r8),intent(in)::k_wbl_mdp ! I [frc] Weibull distribution shape factor (midpoint)
    integer,intent(in)::wnd_mdp_nbr ! I [nbr] Number of wind speed calculations to do
    real(r8),intent(in)::wnd_mdp_max ! I [m s-1] max wind speed for which we calculate
    real(r8),intent(in)::wnd_mdp_min ! I [m s-1] min wind speed for which we calculate
    ! Output
    real(r8),intent(out)::wnd_mdp_wgt(wnd_mdp_nbr) ! O [frc] Wind speed bin weight
    ! Local
    real(r8)::Fmin ! Cumulative distribution of lowest wind in bin
    real(r8)::Fmax ! Cumulative distribution of highest wind in bin
    integer::i ! [idx] Counting index for longitude
    real(r8)::wnd_inc ! [m s-1] Wind increment for a longitude
    real(r8)::wnd_max_bin ! [m s-1] Max wind speed for which we calculate
    integer::wnd_mdp_idx ! [idx] Counting index for wind midpoint
    real(r8)::wnd_min_bin ! [m s-1] Min wind in a bin
    real(r8)::wnd_srt ! [m s-1] Start wind speed for a bin
    real(r8)::sumfracs ! [frac] Sum of all fractions for a longitude
    ! Begin code
    
    ! Initialize
    wnd_srt = wnd_mdp_min ! [m s-1] Lowest wind speed for which we do calculations
    ! Increment at all  longitude
    wnd_inc = (wnd_mdp_max-wnd_mdp_min)/real(wnd_mdp_nbr, r8) 
    
    do wnd_mdp_idx=1,wnd_mdp_nbr
       wnd_min_bin = wnd_srt ! [m s-1] Lowest wind speed in bin
       wnd_max_bin = wnd_srt+wnd_inc ! [m s-1] Highest wind speed in bin
          
       ! Fraction of winds which are below lowest wind speed in bin
       Fmin = 1.0_r8-exp(-1.0_r8*(wnd_min_bin/c_wbl_mdp)**k_wbl_mdp)
          
       ! Fraction of winds which are below highest wind speeds in bin
       Fmax = 1.0_r8-exp(-1.0_r8*(wnd_max_bin/c_wbl_mdp)**k_wbl_mdp)
          
       ! Weight of current wind speed bin
       wnd_mdp_wgt(wnd_mdp_idx) = Fmax-Fmin ! [frc] Wind speed bin weight
          
       ! New start point for next wind speed bin
       wnd_srt = wnd_max_bin
       
    enddo ! end loop over winds
    
    ! Normalize all sums to one
    ! Sum of all wind weights at a longitude
    sumfracs = sum(wnd_mdp_wgt(:)) 
    ! Normalize weights
    wnd_mdp_wgt(:)=wnd_mdp_wgt(:)/sumfracs
    
  end subroutine wnd_mdp_wgt_wbl_get
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine wnd_mdp_bin_get( &
       wnd_mdp_bin, & ! O [m s-1] wind speed in bin
       wnd_mdp_idx, & ! I [nbr] The bin we are interested in
       wnd_mdp_max, & ! I [m s-1] largest wind speed which we calculate
       wnd_mdp_min, & ! I [m s-1] smallest wind speed which we calculate
       wnd_mdp_nbr) ! I [nbr] Weibull wind number
    ! Purpose: 
    ! Given the max and min speed we care about, the max number of bins 
    ! and the bin we want to calculate, we get, from this subroutine
    ! The wind speed in this bin
    !//
    !// Rewritten for Oslo CTM3, works on single values.
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    integer,intent(in)::wnd_mdp_idx ! [idx] The number of the wind bin we are interested in
    integer,intent(in)::wnd_mdp_nbr ! [nbr] Number of winds to calculate
    real(r8),intent(in)::wnd_mdp_max ! [m s-1] Max wind speed for which we calculate
    real(r8),intent(in)::wnd_mdp_min ! [m s-1] Min wind speed for which we calculate
    ! Output
    real(r8),intent(out)::wnd_mdp_bin ! [m s-1] Representative wind speed of bin
    ! Local 
    real(r8)::onebin ! [m s-1] wind range in one wind bin
    real(r8)::totrange ! [m s-1] total range of winds we care about
    real(r8)::wnd_max_bin ! [m s-1] Highest wind in bin
    real(r8)::wnd_min_bin ! [m s-1] Lowest wind in bin
    ! Begin code

    ! Total wind range
    totrange = wnd_mdp_max-wnd_mdp_min ! [m s-1] 
    ! write(6,*)'totrange',totrange
       
    ! Total wind range divided by wind bins
    onebin = 1.0_r8/real(wnd_mdp_nbr, r8)*totrange
    ! write(6,*)'onebin ',onebin
       
    ! Minimum wind speed in this wind bin
    wnd_min_bin = real(wnd_mdp_idx - 1._r8, r8)*onebin + wnd_mdp_min
       
    ! Maximum wind speed in this wind bin
    wnd_max_bin = wnd_min_bin+onebin
       
    ! Average wind speed in this wind bin
    wnd_mdp_bin = 0.5_r8*(wnd_min_bin+wnd_max_bin)
    
  end subroutine wnd_mdp_bin_get
  !// ------------------------------------------------------------------

end module wbl_mdl ! [mdl] Weibull wind speed distribution

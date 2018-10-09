! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dpsdryutl.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Routines to compute dry deposition tendencies

! Usage: 
! use dpsdryutl ! [mdl] Dry deposition utilities

module dpsdryutl ! [mdl] Dry deposition utilities
  use dead_precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  private::CEISK_fst_scl ! [fnc] Complete Elliptic Integral of the Second Kind
  private::cff_drg_Boo71_fst_scl ! [fnc] Drag coefficient for aspherical particles
  private::cff_drg_fst_scl ! [fnc] Drag coefficient
  private::dmt_eqv_lps_fst_scl ! [fnc] Diameter of sphere with same surface area as ellipsoid
  private::lpt_fct_fst_scl ! [fnc] Ellipticity factor of given ellipsoidal aspect ratio
  public::psi_lps_fst_scl ! [fnc] Equivalent diameter divided by semi-minor axis
  private::sph_fct_fst_scl ! [fnc] Sphericity factor of given ellipsoidal aspect ratio
  public::rss_aer_get ! [sbr] Compute aerodynamic resistance
  public::rss_lmn_get ! [sbr] Compute resistance of quasi-laminar layer 
  public::stk_crc_get ! [sbr] Determine bias in Stokes settling velocity
  public::vlc_grv_get ! [sbr] Compute terminal fall speed
  
contains


  !// ------------------------------------------------------------------
  real(r8) function CEISK_fst_scl( & ! [fnc] Complete Elliptic Integral of the Second Kind
       kkk) ! [frc] Modulus k (0 <= k < 1)
    ! Purpose: Complete Elliptic Integral of the Second Kind (CEISK)
    ! CEISK is often denoted K(k), and CEIFK is often denoted E(k)
    ! Source: Polynomial Approximation from AbS64 p. 590 (17.3.34)
    ! Approximation is accurate to <= 2.0e-8
    ! Alternates:
    ! Numerical Recipes First Edition: cel() function
    ! Numerical Recipes Second Edition: rd_[sv]()
    ! Computation of Special Functions: comelp(k,ck,ce)
    ! http://www.esrf.fr/computing/expg/libraries/smf/PROGRAMS/MCOMELP.FOR
    implicit none
    real(r8),intent(in)::kkk ! [frc] Modulus k (0 <= k < 1)
    real(r8) one_mns_kkk_sqr ! [frc] 1-k^2
    real(r8) ak ! [frc] 
    real(r8) bk ! [frc] 
    one_mns_kkk_sqr=1.0_r8-kkk*kkk
    if (kkk == 1.0) then
       CEISK_fst_scl=1.0_r8+300
    else
       ak=(((0.01451196212_r8*one_mns_kkk_sqr+0.03742563713_r8)*one_mns_kkk_sqr &
            +0.03590092383_r8)*one_mns_kkk_sqr+0.09666344259_r8)*one_mns_kkk_sqr+ &
            1.38629436112_r8
       bk=(((0.00441787012_r8*one_mns_kkk_sqr+0.03328355346_r8)*one_mns_kkk_sqr+ &
            0.06880248576_r8)*one_mns_kkk_sqr+0.12498593597_r8)*one_mns_kkk_sqr+0.5_r8
       CEISK_fst_scl=ak-bk*log(one_mns_kkk_sqr)
    endif ! kkk != 1.0
  end function CEISK_fst_scl
  !// ------------------------------------------------------------------
  
  !// ------------------------------------------------------------------
  real(r8) function cff_drg_Boo71_fst_scl( & ! [fnc] Drag coefficient for aspherical particles
       ryn_nbr, & ! I [frc] Reynolds number
       sph_fct) ! I [frc] Sphericity factor
    ! Purpose: "Fast scalar" drag coefficient for aspherical particles at given Reynolds number
    implicit none
    real(r8),intent(in)::ryn_nbr ! I [frc] Reynolds number
    real(r8),intent(in)::sph_fct ! I [frc] Sphericity factor
    cff_drg_Boo71_fst_scl=cff_drg_fst_scl(ryn_nbr)
    if (sph_fct /= 1.0_r8) cff_drg_Boo71_fst_scl=cff_drg_Boo71_fst_scl+ &
         (24.0_r8/ryn_nbr)*10.0_r8*(1.0_r8-sph_fct)*ryn_nbr**0.35_r8/sph_fct ! Gin03 p. 2 (4)
  end function cff_drg_Boo71_fst_scl
  !// ------------------------------------------------------------------
  
  !// ------------------------------------------------------------------
  real(r8) function cff_drg_fst_scl( & ! [fnc] Drag coefficient
       ryn_nbr) ! I [frc] Reynolds number
    ! Purpose: "Fast scalar" drag coefficient for given Reynolds number
    ! Taken from Seinfeld and Pandis (1997) as reported in SeP97 p. 463 (8.32)
    ! Range of validity is < Re < 2.0e5
    ! In order to maximize chances of compiler inlining this function,
    ! it avoids error checking and diagnostics

    ! SeP97 expressions are not smoothly matched at boundaries
    ! This causes finite C_D jumps between continuously varying particle sizes
    ! For dust, jump occurs near D=80 um because Re changes from < 2 to > 2
    ! This jump nearly doubles sedimentation speed and is very unrealistic
    ! One approach to solve this is to blend solutions over limited range

    ! Pruppacher and Klett offer alternative formulations for intermediate C_D
    ! PrK78 p. 294 (10-50)--(10.52) and PrK98 p. 373 (10-51)--(10-53) have ln(Re/2)
    ! but SeP97 p. 463 (8.32) has ln(Re*2).
    ! Thus it appears one or the other has a typo

    ! GinO3 implies that Re < 2 for D < 100 um but really Re < 2 for D < 80 um
    ! so his figures are misleading.
    implicit none
    real(r8),parameter::cst_M_LN2l=0.6931471805599453094172321214581766 ! [frc] log_e(2.0)
    real(r8),parameter::cst_M_EULERl=0.577215664901532860606512090082 ! [frc] Euler's constant
    real(r8),parameter::ryn_nbr_ntr_min=2.0 ! [frc] Minimum Reynolds number for intermediate regime
    real(r8),parameter::ryn_nbr_ntr_max=5.0 ! [frc] Maximum Reynolds number for intermediate regime
    real(r8),intent(in)::ryn_nbr ! I [frc] Reynolds number
    real(r8) wgt_mpr ! [frc] Weight of empirical parameterization regime
    if (ryn_nbr < 0.1_r8) then
       cff_drg_fst_scl=24.0_r8/ryn_nbr ! Stokes' law Sep97 p. 463 (8.32)
    else if (ryn_nbr < ryn_nbr_ntr_min) then
       cff_drg_fst_scl=(24.0/ryn_nbr)* &
            (1.0_r8+3.0_r8*ryn_nbr/16.0_r8+ &
            9.0_r8*ryn_nbr*ryn_nbr* &
            (log(0.5_r8*ryn_nbr)+cst_M_EULERl+5.0_r8*cst_M_LN2l/3.0_r8-323.0_r8/360.0_r8) &
            /160.0_r8+ &
            27.0_r8*ryn_nbr*ryn_nbr*ryn_nbr*log(0.5_r8*ryn_nbr)/640.0_r8)
    else if (ryn_nbr < ryn_nbr_ntr_max) then
       wgt_mpr=(ryn_nbr-ryn_nbr_ntr_min)/(ryn_nbr_ntr_max-ryn_nbr_ntr_min) ! [frc] Weight of empirical parameterization regime
       cff_drg_fst_scl=(1.0-wgt_mpr)*(24.0/ryn_nbr)* &
            (1.0_r8+3.0_r8*ryn_nbr/16.0_r8+ &
            9.0_r8*ryn_nbr*ryn_nbr* &
            (log(0.5_r8*ryn_nbr)+cst_M_EULERl+5.0_r8*cst_M_LN2l/3.0_r8-323.0_r8/360.0_r8) &
            /160.0_r8+ &
            27.0_r8*ryn_nbr*ryn_nbr*ryn_nbr*log(0.5_r8*ryn_nbr)/640.0_r8) &
            +wgt_mpr*(24.0_r8/ryn_nbr)*(1.0_r8+0.15_r8*ryn_nbr**0.687_r8)
    else if (ryn_nbr < 500.0_r8) then
       cff_drg_fst_scl=(24.0_r8/ryn_nbr)* &
            (1.0_r8+0.15_r8*ryn_nbr**0.687_r8) ! Sep97 p. 463 (8.32)
    else if (ryn_nbr < 2.0e5_r8) then
       cff_drg_fst_scl=0.44_r8 ! Sep97 p. 463 (8.32)
    else
       write (6,'(a,es9.2)') "ryn_nbr = ",ryn_nbr
       stop 'ERROR: Reynolds number too large in cff_drg_fst_scl()'
    endif               ! end else
  end function cff_drg_fst_scl ! [fnc] Drag coefficient
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function dmt_eqv_lps_fst_scl( & ! [fnc] Diameter of sphere with same surface area as ellipsoid
       asp_rat_lps, & ! I [frc] Ellipsoidal aspect ratio
       rds_lps_b) ! I [m] Semi-minor axis (b) of ellipsoid
    ! Purpose: "Fast scalar" equivalent diameter for given ellipsoidal aspect ratio
    ! Equivalent diameter is diameter of sphere with same surface area as given ellipsoid
    ! Routine computes dmt_eqv_lps from semi-minor axis b and aspect ratio asp_rat_lps = a/b
    ! Taken from Ginoux (2003) Gin03 p. 2 (10)
    ! In order to maximize chances of compiler inlining this function,
    ! it avoids error checking and diagnostics
    implicit none
    real(r8),intent(in)::asp_rat_lps ! I [frc] Ellipsoidal aspect ratio
    real(r8),intent(in)::rds_lps_b ! I [m] Semi-minor axis (b) of ellipsoid
    dmt_eqv_lps_fst_scl=rds_lps_b*psi_lps_fst_scl(asp_rat_lps) ! [frc] Diameter of sphere with same surface area as ellipsoid Gin03 p. 2 (10)
  end function dmt_eqv_lps_fst_scl
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function lpt_fct_fst_scl( & ! [fnc] Ellipticity factor of given ellipsoidal aspect ratio
       asp_rat_lps) ! I [frc] Ellipsoidal aspect ratio
    ! Purpose: "Fast scalar" ellipticity factor for a given ellipsoidal aspect ratio
    ! The ellipticity factor lpt_fct of a prolate spheroid is sqrt(a^2-b^2)/a
    ! This routine computes lpt_fct from the aspect ratio asp_rat_lps = a/b
    ! Taken from Ginoux (2003) Gin03 p. 2 (7), http:!mathworld.wolfram.com/ProlateSpheroid.html
    ! In order to maximize chances of compiler inlining this function,
    ! it avoids error checking and diagnostics
    implicit none
    real(r8),intent(in)::asp_rat_lps ! I [frc] Ellipsoidal aspect ratio
    lpt_fct_fst_scl=sqrt(asp_rat_lps*asp_rat_lps-1.0_r8)/asp_rat_lps ! [frc] Ellipticity factor Gin03 p. 2 (8)
  end function lpt_fct_fst_scl
  !// ------------------------------------------------------------------
  !// ------------------------------------------------------------------
  real(r8) function psi_lps_fst_scl( & ! [fnc] Equivalent diameter divided by semi-minor axis
       asp_rat_lps) ! I [frc] Ellipsoidal aspect ratio
    ! Purpose: "Fast scalar" equivalent diameter divided by semi-minor axis for given ellipsoidal aspect ratio
    ! Equivalent diameter is diameter of sphere with same surface area as given ellipsoid
    ! Psi is non-dimensional function defined such that Psi*b=dmt_eqv, so that
    ! Psi is equivalent diameter divided by semi-minor axis
    ! Routine computes psi_lps from aspect ratio asp_rat_lps = a/b
    ! Ginoux (2003) Gin03 p. 2 (10) has a typo/error in Psi definition
    ! Denominator under radical is "asp_rat_lps^2-1" not "1-asp_rat_lps^2"
    ! Former works fine for asp_rat >= 1 but latter is imaginary!
    implicit none
    real(r8),intent(in)::asp_rat_lps ! I [frc] Ellipsoidal aspect ratio
    real(r8) asp_rat_lps_sqr ! [frc] Ellipsoidal aspect ratio squared
    if(asp_rat_lps == 1.0_r8) then 
       ! Return limiting value for spherical case
       psi_lps_fst_scl=2.0_r8
       return
    endif ! asp_rat_lps != 1.0
    asp_rat_lps_sqr=asp_rat_lps*asp_rat_lps ! [frc] Ellipsoidal aspect ratio squared
    psi_lps_fst_scl=sqrt(2.0_r8+2.0_r8*asp_rat_lps_sqr/sqrt(asp_rat_lps_sqr-1.0_r8)*asin(sqrt(1.0_r8-1.0_r8/asp_rat_lps_sqr))) ! [frc] Equivalent diameter divided by semi-minor axis Gin03 p. 2 (10)
  end function psi_lps_fst_scl
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  real(r8) function sph_fct_fst_scl( & ! [fnc] Sphericity factor of given ellipsoidal aspect ratio
       asp_rat_lps) ! I [frc] Ellipsoidal aspect ratio
    ! Purpose: "Fast scalar" sphericity factor for a given ellipsoidal aspect ratio
    ! Taken from Ginoux (2003) Gin03 p. 2 (8)
    ! In order to maximize chances of compiler inlining this function,
    ! it avoids error checking and diagnostics */
    implicit none
    real(r8),intent(in)::asp_rat_lps ! I [frc] Ellipsoidal aspect ratio
    real(r8) asp_rat_lps_sqr ! [frc] Ellipsoidal aspect ratio squared
    if(asp_rat_lps == 1.0_r8) then
       sph_fct_fst_scl=1.0_r8
       return
    endif
    asp_rat_lps_sqr=asp_rat_lps*asp_rat_lps ! [frc] Ellipsoidal aspect ratio squared
    sph_fct_fst_scl=2.0_r8*asp_rat_lps**(2.0_r8/3.0_r8) ! [frc] Sphericity factor Gin03 p. 2 (8)
    sph_fct_fst_scl=sph_fct_fst_scl/(1.0_r8+(asp_rat_lps_sqr/ &
         sqrt(asp_rat_lps_sqr-1.0_r8)*asin(sqrt(1.0_r8-1.0_r8/asp_rat_lps_sqr)))) ! [frc] Sphericity factor Gin03 p. 2 (8)
  end function sph_fct_fst_scl
  !// ------------------------------------------------------------------

  !// ------------------------------------------------------------------
  subroutine rss_aer_get( & ! [sbr] Compute aerodynamic resistance
       hgt_mdp,             & ! I [m] Midlayer height above surface
       hgt_zpd,             & ! I [m] Zero plane displacement height
       mno_lng,             & ! I [m] Monin-Obukhov length
       rgh_mmn,             & ! I [m] Roughness length momentum
       rss_aer,             & ! O [s m-1] Aerodynamic resistance
       wnd_frc)             ! I [m s-1] Surface friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given the surface layer structure and dynamics properties,
    ! compute and return the aerodynamic resistance
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    real(r8),intent(in) :: mno_lng, &       ! [m] Monin-Obukhov length
                           rgh_mmn, &       ! [m] Roughness length momentum
                           wnd_frc, &       ! [m s-1] Surface friction velocity from pbldif()
                           hgt_mdp, &       ! [m] Midlayer height above surface
                           hgt_zpd          ! [m] Zero plane displacement
    ! Output
    real(r8),intent(out) :: rss_aer         ! [s m-1] Aerodynamic resistance

    ! Local
    real(r8) eta_sqr_rlm          ! Eta squared at roughness height
    real(r8) eta_sqr_gcm          ! Eta squared at GCM layer height
    real(r8) eta_rlm              ! Eta at roughness height
    real(r8) eta_gcm              ! Eta at GCM layer height
    real(r8) nmr_rlm              ! Numerator
    real(r8) dnm_gcm              ! Denominator
    real(r8) tmp4                 ! Correction to neutral atmosphere factor
    real(r8) tmp5                 ! Neutral atmosphere factor
    real(r8) mno_prm_rlm          ! [frc] Monin-Obukhov parameter at roughness height
    real(r8) mno_prm_gcm          ! [frc] Monin-Obukhov parameter at GCM layer height
    !// ------------------------------------------------------------------
    ! Main Code
    
    ! Sanity check
#if 0
    if (rgh_mmn <= 0.0.or.hgt_mdp <= hgt_zpd) then
       write(6,'(a,3(a,es8.1),a)') 'rss_aer_get: ','rgh_mmn = ',rgh_mmn, &
            ', hgt_mdp = ',hgt_mdp,' m, hgt_zpd = ',hgt_zpd,' m s-1'
    endif                  ! endif err
#endif /* !0 */
      
    ! Compute stability parameter
    ! Maximum value is 1 because stability correction function is valid only for zeta < 1, e.g., Bon96 p. 52, Bru82 p. 71, SeP97 p. 963
    mno_prm_rlm = min(rgh_mmn/mno_lng,1.0_r8) ! [frc]
    mno_prm_gcm = min((hgt_mdp-hgt_zpd)/mno_lng,1.0_r8) ! [frc]
    if (mno_lng < 0.0_r8) then
       ! Difference between unstable corrections
       eta_sqr_rlm = sqrt(1.0_r8-16.0_r8*mno_prm_rlm)
       eta_sqr_gcm = sqrt(1.0_r8-16.0_r8*mno_prm_gcm)
       eta_rlm = sqrt(eta_sqr_rlm)
       eta_gcm = sqrt(eta_sqr_gcm)
       nmr_rlm = (eta_sqr_rlm+1.0_r8)*(eta_rlm+1.0_r8)*(eta_rlm+1.0_r8)
       dnm_gcm = (eta_sqr_gcm+1.0_r8)*(eta_gcm+1.0_r8)*(eta_gcm+1.0_r8)
       tmp4 = log(nmr_rlm/dnm_gcm)+2.0_r8*(atan(eta_gcm)-atan(eta_rlm)) ! [frc]
    else                   ! not stable
       ! Difference between stable corrections
       tmp4 = 5.0_r8*(mno_prm_gcm-mno_prm_rlm) ! [frc]
    endif                  ! not stable
    tmp5 = log((hgt_mdp-hgt_zpd)/rgh_mmn)
    rss_aer = (tmp4+tmp5)/(cst_von_krm*wnd_frc) ! [s m-1] Bon96 p. 54
  
    !// ------------------------------------------------------------------
  end subroutine rss_aer_get                       ! end rss_aer_get()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine rss_lmn_get( & ! [sbr] Compute resistance of quasi-laminar layer 
       dmt_aer,             & ! I [m] Particle diameter
       dns_aer,             & ! I [kg m-3] Particle density
       dns_mdp,             & ! I [kg m-3] Midlayer density
       prs_mdp,             & ! I [Pa] Pressure
       rss_lmn,             & ! O [s m-1] Quasi-laminar layer resistance
       shm_nbr,             & ! O [frc] Schmidt number
       stk_crc,             & ! I [frc] Correction to Stokes settling velocity
       stk_nbr,             & ! O [frc] Stokes number
       tpt_mdp,             & ! I [K] Temperature
       wnd_frc)             ! I [m s-1] Friction velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given size and density of aerosol, and surface environmental properties, compute quasi-laminar layer resistance
    ! rss_lmn_get() describes physics occuring within a few meters (at most) of the zero plane displacement height
    ! Thus rss_lmn_get() requires the aerodynamic pressure and temperature which is closer to lowest interface than to lowest midlayer
    ! Currently we ignore this subtlety and use lowest midlayer quantities
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstgrd, only: dst_nbr ! [mdl] Dust grid sizes
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    ! Parameters
    ! Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
    real(r8),parameter :: &
         cst_Boltzmann = 1.38063d-23, &  ![J K-1] Boltzmann's constant
         shm_nbr_xpn_lnd = -2.d0/3.d0, & ![frc] Exponent for aerosol-diffusion
                                         !dependence on Schmidt number over land
         shm_nbr_xpn_ocn = -0.5d0        ![frc] Exponent for aerosol-diffusion
                                         !dependence on Schmidt number over ocean
    ! Input
    real(r8),intent(in) :: &
         dns_aer(dst_nbr), &  ! I [kg m-3] Particle density
         dmt_aer(dst_nbr), &  ! I [m] Particle diameter
         stk_crc(dst_nbr), &  ! I [frc] Correction to Stokes settling velocity
         dns_mdp, &           ! I [kg m-3] Midlayer density
         prs_mdp, &           ! I [Pa] Pressure 
         tpt_mdp, &           ! I [K] Temperature
         wnd_frc              ! I [m s-1] Friction velocity
    ! Output
    real(r8), intent(out) :: &
         rss_lmn(dst_nbr), &  ! O [s m-1] Quasi-laminar layer resistance
         stk_nbr(dst_nbr), &  ! O [frc] Stokes number
         shm_nbr(dst_nbr)     ! O [frc] Schmidt number
    ! Local
    real(r8) :: pi               ! [frc] 3
    integer :: k                 ! [idx] Counting index
    integer :: m                 ! [idx] Counting index
    real(r8) :: dff_aer(dst_nbr) ! [m2 s-1] Brownian diffusivity of particle
    real(r8) :: mfp_atm          ! [m] Mean free path of air
    real(r8) :: shm_nbr_xpn      ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt number
    real(r8) :: slp_crc(dst_nbr) ! [frc] Slip correction factor
    real(r8) :: tmp              ! [frc] Factor in rss_lmn computation
    real(r8) :: vlc_grv(dst_nbr) ! [m s-1] Settling velocity
    real(r8) :: vsc_dyn_atm      ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm      ! [m2 s-1] Kinematic viscosity of air
    !// ------------------------------------------------------------------
    ! Main Code
    ! Initialize some variables
    pi=4.0_DBLKIND*atan(1.0_DBLKIND)        ! [frc] 3

    ! Size-independent thermokinetic properties
    vsc_dyn_atm = 1.72e-5_r8*((tpt_mdp/273.0_r8)**1.5_r8)*393.0_r8/(tpt_mdp+120.0_r8) ! [kg m-1 s-1] RoY94  p. 102
    mfp_atm = 2.0_r8*vsc_dyn_atm/(prs_mdp*sqrt(8.0_r8*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp))) ! [m] Mean free path of air SeP97 p. 455
    vsc_knm_atm = vsc_dyn_atm/dns_mdp ! [m2 s-1] Kinematic viscosity of air

    do m=1,dst_nbr
       slp_crc(m) = 1.0_r8+2.0_r8*mfp_atm*(1.257_r8+0.4_r8*exp(-1.1_r8*dmt_aer(m)/(2.0_r8*mfp_atm)))/dmt_aer(m) ! [frc] Slip correction factor SeP97 p. 464
       vlc_grv(m) = (1.0_r8/18.0_r8)*dmt_aer(m)*dmt_aer(m)*dns_aer(m)*grv_sfc*slp_crc(m)/vsc_dyn_atm ! [m s-1] Stokes' settling velocity SeP97 p. 466
       vlc_grv(m) = vlc_grv(m)*stk_crc(m) ! [m s-1] Correction to Stokes settling velocity
    end do                     ! end loop over size
    
    do m=1,dst_nbr
       stk_nbr(m) = vlc_grv(m)*wnd_frc*wnd_frc/(grv_sfc*vsc_knm_atm) ! [frc] Stokes number SeP97 p. 965
       dff_aer(m) = cst_Boltzmann*tpt_mdp*slp_crc(m)/(3.0_r8*pi*vsc_dyn_atm*dmt_aer(m)) ! [m2 s-1] Brownian diffusivity of particle SeP97 p. 474
       shm_nbr(m) = vsc_knm_atm/dff_aer(m) ! [frc] Schmidt number SeP97 p. 972
       shm_nbr_xpn = shm_nbr_xpn_lnd ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt number          
       ! fxm: Turning this on dramatically reduces deposition velocity in low wind regimes
       ! Schmidt number exponent is -2/3 over solid surfaces and -1/2 over liquid surfaces SlS80 p. 1014
       ! if (oro(i) == 0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else shm_nbr_xpn=shm_nbr_xpn_lnd ! [frc] Surface-dependent exponent for aerosol-diffusion dependence on Schmidt number
       tmp = shm_nbr(m)**shm_nbr_xpn+10.0_r8**(-3.0_r8/stk_nbr(m))
       rss_lmn(m) = 1.0_r8/(tmp*wnd_frc) ! [s m-1] SeP97 p. 972, 965
    end do                     ! end loop over size

    !// ------------------------------------------------------------------
  end subroutine rss_lmn_get                      ! end rss_lmn_get()
  !// ------------------------------------------------------------------
  


  !// ------------------------------------------------------------------
  subroutine stk_crc_get( & ! [sbr] Determine bias in Stokes settling velocity
       asp_rat_lps,         & ! I [frc] Ellipsoidal aspect ratio
       dmt_prt,             & ! I [m] Particle diameter
       dns_prt,             & ! I [kg m-3] Particle density
       stk_crc,             & ! O [frc] Correction to Stokes settling velocity
       sz_nbr)              ! I [nbr] Number of particle sizes
    !// ------------------------------------------------------------------
    ! Purpose: Compute factor which corrects Stokes settling velocity to actual
    ! settling velocity when Reynolds' number is greater than 0.1
    ! stk_crc_get() is called by dst_psd_ini()
    ! This implementation computes a time-invariant, size-dependent ratio which
    ! corrects the solution of the (linear) Stokes settling velocity to the full
    ! (non-linear, iterative solution) for the settling velocity.
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstgrd ! [mdl] Dust grid sizes
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
    real(r8),parameter::cst_Boltzmann=1.38063e-23 ! [J K-1] Boltzmann's constant
    real(r8),parameter::prs_mdp=100000.0 ! [Pa] Pressure 
    real(r8),parameter::tpt_mdp=295.0 ! [K] Temperature
    real(r8),parameter::tpt_vrt=295.0 ! [K] Virtual temperature
    ! Input
    integer,intent(in)::sz_nbr ! I [nbr] Number of particle sizes
    real(r8),intent(in)::asp_rat_lps(sz_nbr) ! I [frc] Ellipsoidal aspect ratio
    real(r8),intent(in)::dmt_prt(sz_nbr) ! I [m] Particle diameter
    real(r8),intent(in)::dns_prt(sz_nbr) ! I [kg m-3] Particle density
    ! Output
    real(r8),intent(out)::stk_crc(sz_nbr) ! O [frc] Correction to Stokes settling velocity
    ! Local
    real(r8) pi       ! [frc] 3
    integer itr_idx           ! [idx] Counting index
    integer sz_idx            ! [idx] Counting index
    real(r8) cff_drg_grv(sz_nbr)  ! [frc] Drag coefficient at terminal velocity
    real(r8) dns_mdp              ! [kg m-3] Midlayer density
    real(r8) eps_crr              ! [frc] Current relative accuracy
    real(r8) eps_max              ! [frc] Relative accuracy for convergence
    real(r8) mfp_atm              ! [m] Mean free path of air
    real(r8) ryn_nbr_grv(sz_nbr)  ! [frc] Reynolds number at terminal velocity
    real(r8) slp_crc(sz_nbr)      ! [frc] Slip correction factor
    real(r8) vlc_grv(sz_nbr)      ! [m s-1] Settling velocity
    real(r8) vlc_grv_old          ! [m s-1] Previous gravitational settling velocity
    real(r8) vlc_stk(sz_nbr)      ! [m s-1] Stokes settling velocity
    real(r8) vsc_dyn_atm          ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) vsc_knm_atm          ! [m2 s-1] Kinematic viscosity of atmosphere 
    ! Variable required by aspherical considerations
    real(r8) CEISK ! [frc] Complete Elliptic Integral of the Second Kind
    real(r8) fct_asp ! [frc] Factor by which spherical vlc_grv^2 differs from spherical
    real(r8) lpt_fct ! [frc] Ellipticity factor
    real(r8) psi_lps ! [frc] Equivalent diameter divided by semi-minor axis
    real(r8) sph_fct ! [frc] Sphericity factor
    !// ------------------------------------------------------------------
    ! Initialize some variables
    pi=4.0_DBLKIND*atan(1.0_DBLKIND)        ! [frc] 3

    eps_max=1.0e-4_r8            ! [frc] Relative accuracy for convergence
    dns_mdp=prs_mdp/(tpt_vrt*gas_cst_dry_air) ! [kg m-3]
    ! Size-independent thermokinetic properties
    vsc_dyn_atm=1.72e-5_r8*((tpt_mdp/273.0_r8)**1.5_r8)*393.0_r8/(tpt_mdp+120.0_r8) ! [kg m-1 s-1] RoY94  p. 102
    mfp_atm=2.0_r8*vsc_dyn_atm/(prs_mdp*sqrt(8.0_r8*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp))) ! [m] Mean free path of air SeP97 p. 455
    vsc_knm_atm=vsc_dyn_atm/dns_mdp ! [m2 s-1] Kinematic viscosity of air
    
    do sz_idx=1,sz_nbr
       slp_crc(sz_idx)=1.0_r8+2.0_r8*mfp_atm* &
       (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_prt(sz_idx)/(2.0_r8*mfp_atm)))/ &
       dmt_prt(sz_idx) ! [frc] Slip correction factor SeP97 p. 464
       vlc_stk(sz_idx)=(1.0_r8/18.0_r8)*dmt_prt(sz_idx)*dmt_prt(sz_idx)* &
            dns_prt(sz_idx)*grv_sfc*slp_crc(sz_idx)/vsc_dyn_atm ! [m s-1] Stokes settling velocity SeP97 p. 466
    end do                     ! end loop over size
    
    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for vlc_grv SeP97 p. 466 (8.42)
    ! For larger Re, inertial effects become important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical drag coefficient causes 60% errors for D = 200 um SeP97 p. 468
    
    ! Iterative solution for drag coefficient, Reynolds number, and terminal velocity
    do sz_idx=1,sz_nbr
       ! Initialize accuracy and counter
       eps_crr=eps_max+1.0_r8    ! [frc] Current relative accuracy
       itr_idx=0              ! [idx] Counting index
       ! Initial guess for vlc_grv is exact for Re < 0.1
       vlc_grv(sz_idx)=vlc_stk(sz_idx)  ! [m s-1]
       lpt_fct=lpt_fct_fst_scl(asp_rat_lps(sz_idx)) ! [frc] Ellipticity factor Gin03 p. 2 (8)
       psi_lps=psi_lps_fst_scl(asp_rat_lps(sz_idx)) ! [frc] Equivalent diameter divided by semi-minor axis Gin03 p. 2 (10)
       sph_fct=sph_fct_fst_scl(asp_rat_lps(sz_idx)) ! [frc] Sphericity factor Gin03 p. 2 (8)
       do while(eps_crr > eps_max)
          ! Save terminal velocity for convergence test
          vlc_grv_old=vlc_grv(sz_idx) ! [m s-1] 
          ryn_nbr_grv(sz_idx)=vlc_grv(sz_idx)*dmt_prt(sz_idx)/vsc_knm_atm ! [frc] SeP97 p. 460
          ! Update drag coefficient based on new Reynolds number
          ! cff_drg_grv(sz_idx)=cff_drg_fst_scl(ryn_nbr_grv(sz_idx))
          cff_drg_grv(sz_idx)=cff_drg_Boo71_fst_scl(ryn_nbr_grv(sz_idx),sph_fct) ! [frc] Drag coefficient at terminal velocity
          CEISK=CEISK_fst_scl(lpt_fct) ! [frc] Complete Elliptic Integral of the Second Kind
          fct_asp=pi/(CEISK*psi_lps) ! [frc] Factor by which spherical vlc_grv^2 differs from spherical
          ! Update terminal velocity based on new Reynolds number and drag coefficient
          vlc_grv(sz_idx)=sqrt(4.0_r8*grv_sfc*dmt_prt(sz_idx)*slp_crc(sz_idx)*(dns_prt(sz_idx)-dns_mdp)*fct_asp/ &
               (3.0_r8*cff_drg_grv(sz_idx)*dns_mdp)) ! [m s-1] Terminal velocity SeP97 p. 467 (8.44)
          eps_crr=abs((vlc_grv(sz_idx)-vlc_grv_old)/vlc_grv(sz_idx)) ! Relative convergence
          if (itr_idx == 12) then
             ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0 due to discontinuities in derivative of drag coefficient
             vlc_grv(sz_idx)=0.5_r8*(vlc_grv(sz_idx)+vlc_grv_old) ! [m s-1]
          endif               ! endif
          if (itr_idx > 20) then
             write (6,'(a)') 'dst: Terminal velocity not converging in stk_crc_get(), breaking loop...'
             goto 100         ! Jump to next iteration
          endif               ! endif
          itr_idx=itr_idx+1
       end do                 ! end while
100    continue               ! Label to jump to when iteratio does not converge
    end do                    ! end loop over size
    
    ! Compute factors to convert Stokes' settling velocities to actual settling velocities
    do sz_idx=1,sz_nbr
       stk_crc(sz_idx)=vlc_grv(sz_idx)/vlc_stk(sz_idx) ! [frc] Correction to Stokes settling velocity
    end do                     ! end loop over size

  end subroutine stk_crc_get                       ! end stk_crc_get()
  !// ------------------------------------------------------------------



  !// ------------------------------------------------------------------
  subroutine vlc_grv_get( & ! [sbr] Compute terminal fall speed
       dmt_prt,             & ! I [m] Particle diameter
       dns_prt,             & ! I [kg m-3] Particle density
       prs_mdp,             & ! I [Pa] Pressure
       stk_crc,             & ! I [frc] Correction to Stokes settling velocity
       tpt_mdp,             & ! I [K] Temperature
       vlc_grv)             ! O [m s-1] Settling velocity
    !// ------------------------------------------------------------------
    ! Purpose: Given size and density of particle, and column thermodynamic profile,
    ! compute terminal fall speed
    !//
    !// Rewritten for Oslo CTM3, works on only one (lon_idx,lat_idx).
    !// Removed longitude loop.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstgrd, only: dst_nbr ! [mdl] Dust grid sizes
    use pmgrid, only: plev ! [mdl] Spatial resolution parameters
    use dstcst ! [mdl] Physical constants for dust routines
    implicit none
    !// ------------------------------------------------------------------
    ! Parameters
    ! Apparently Boltzmann's constant is not set anywhere by CCM or MATCH
    real(r8),parameter :: cst_Boltzmann=1.38063e-23 ! [J K-1] Boltzmann's constant
    ! Input
    real(r8),intent(in) :: dns_prt(dst_nbr), &     ! I [kg m-3] Particle density
                           dmt_prt(dst_nbr), &     ! I [m] Particle diameter
                           prs_mdp(plev), &  ! I [Pa] Pressure 
                           tpt_mdp(plev), &  ! I [K] Temperature
                           stk_crc(dst_nbr)     ! I [frc] Correction to Stokes settling velocity
    ! Output
    real(r8),intent(out) :: vlc_grv(plev,dst_nbr) ! O [m s-1] Settling velocity

    ! Local
    real(r8) :: pi                ! [frc] 3
    integer :: k                  ! [idx] Counting index
    integer :: m                  ! [idx] Counting index
    real(r8) :: mfp_atm(plev)     ! [m] Mean free path of air
    real(r8) :: slp_crc           ! [frc] Slip correction factor
    real(r8) :: vsc_dyn_atm(plev) ! [kg m-1 s-1] Dynamic viscosity of air
    !// ------------------------------------------------------------------
    ! Main Code
    ! Initialize some variables
    pi=4.0_DBLKIND*atan(1.0_DBLKIND)        ! [frc] 3
    
    ! Size-independent thermokinetic properties
    do k = 1,plev
       vsc_dyn_atm(k) = 1.72e-5_r8*((tpt_mdp(k)/273.0_r8)**1.5_r8)*393.0_r8/(tpt_mdp(k)+120.0_r8) ! [kg m-1 s-1] RoY94  p. 102
       mfp_atm(k) = 2.0_r8*vsc_dyn_atm(k)/(prs_mdp(k)*sqrt(8.0_r8*mmw_dry_air/(pi*gas_cst_unv*tpt_mdp(k)))) ! [m] Mean free path of air SeP97 p. 455
    end do                    ! end loop over lev
    
    do m=1,dst_nbr
       do k=1,plev
          slp_crc = 1.0_r8+2.0_r8*mfp_atm(k)*(1.257_r8+0.4_r8*exp(-1.1_r8*dmt_prt(m)/(2.0_r8*mfp_atm(k))))/dmt_prt(m) ! [frc] Slip correction factor SeP97 p. 464
          vlc_grv(k,m) = (1.0_r8/18.0_r8)*dmt_prt(m)*dmt_prt(m)*dns_prt(m)*grv_sfc*slp_crc/vsc_dyn_atm(k) ! [m s-1] Stokes settling velocity SeP97 p. 466
          vlc_grv(k,m) = vlc_grv(k,m)*stk_crc(m) ! [m s-1] Corrected Stokes settling velocity
       end do                 ! end loop over lev
    end do                     ! end loop over size

    !// ------------------------------------------------------------------
  end subroutine vlc_grv_get                       ! end vlc_grv_get()
  !// ------------------------------------------------------------------
 
end module dpsdryutl

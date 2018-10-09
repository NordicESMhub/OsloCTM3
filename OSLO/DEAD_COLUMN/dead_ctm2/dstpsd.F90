! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstpsd.F90,v 1.2 2003/07/11 14:01:26 alfgr Exp $

! Purpose: Characterize and initialize mineral dust particle size distributions
! Variables in dstpsd are not needed every timestep by source and sink routines
! Related microphysical quantities, needed every timestep, are stored in dstaer.F90
! Most common variables in dstpsd are initialized in dst_psd_ini(),
! dst_psd_ini() is called by CCM:physics/inti(), MATCH:src/inirun()
! The three modes of the source distribution are initialized (and saved) 
! in the first routine to use dstpsd.

! dstpsd MUST have access to modules dstgrd and pmgrid

! Usage:
! use dstpsd ! [mdl] Dust particle size distributions

module dstpsd ! [mdl] Dust particle size distributions
  use precision ! [mdl] Precision r8, i8, ...
  use dstgrd,only:dst_nbr,dst_src_nbr ! [mdl] Dust grid sizes
  use dpsdryutl,only:psi_lps_fst_scl ! [mdl] Dry deposition utilities
  use pmgrid,only:plond,plon ! [mdl] Spatial resolution parameters 
  implicit none
  save ! [stt] Changes to common variables are sticky
  private ! [stt] Symbols are private unless individually qualified as public
  public::asp_rat_lps_set ! [fnc] Set ellipsoidal aspect ratio
  public::dst_psd_ini ! [sbr] Initialize particle size distributions
  public::dst_psd_src_ini ![sbr] Initialize psd used in transport
  
  real(r8)::dmt_vma_src(dst_src_nbr) ! [m] Mass median diameter of source modes
  real(r8)::gsd_anl_src(dst_src_nbr) ! [m] Geometric standard deviation

  ! Global model uses this value, BXM may set this from command line
  real(r8)::asp_rat_lps_dfl=1.0 ! [frc] Ellipsoidal aspect ratio

  ! Computed at run time in dst_psd_ini()
  real(r8) asp_rat_lps(dst_nbr) ! [frc] Ellipsoidal aspect ratio
  real(r8) dmt_grd(dst_nbr+1)   ! [m] Particle diameter grid
  real(r8) dmt_max(dst_nbr)     ! [m] Maximum diameter in bin
  real(r8) dmt_min(dst_nbr)     ! [m] Minimum diameter in bin
  real(r8) dmt_mjr(dst_nbr)     ! [m] Major axis of spheroid
  real(r8) dmt_mnr(dst_nbr)     ! [m] Minor axis of spheroid
  real(r8) dmt_naa(dst_nbr)     ! [m] Number mean particle size
  real(r8) dmt_nma(dst_nbr)     ! [m] Number median particle diameter
  real(r8) dmt_nmr(dst_nbr)     ! [m] Number median diameter resolved
  real(r8) dmt_nwr(dst_nbr)     ! [m] Number mean diameter resolved
  real(r8) dmt_sma(dst_nbr)     ! [m] Surface median particle size
  real(r8) dmt_smr(dst_nbr)     ! [m] Surface area median diameter resolved
  real(r8) dmt_swa(dst_nbr)     ! [m] Surface area weighted mean particle size
  real(r8),public::dmt_swr(dst_nbr) ! [m] Surface area weighted diameter resolved
  real(r8) dmt_vma(dst_nbr)     ! [m] Mass median diameter
  real(r8) dmt_vmr(dst_nbr)     ! [m] Mass median diameter resolved
  real(r8) dmt_vwa(dst_nbr)     ! [m] Mass weighted mean particle size
  real(r8) gsd_anl(dst_nbr)     ! [frc] Geometric standard deviation
  real(r8) nbr_spc_rsl(dst_nbr) ! [# kg-1] Specific concentration resolved 
  real(r8) vlm_spc_rsl(dst_nbr) ! [m3 kg-1] Specific volume resolved 

  real(r8),public::ovr_src_snk_frc(dst_src_nbr,dst_nbr) ! [frc] Overlap of src with snk

contains

  subroutine dst_psd_src_ini()
    !Purpose: Set values of the size distribution at source
    !Perform choice of size distribution here
    use dstaer,only:mss_frc_src
    implicit none
    real(r8) dmt_vma_srcx(dst_src_nbr)       ! [m] Mass median diameter of source mode
    real(r8) gsd_anl_srcx(dst_src_nbr)       ! [frc] geometric standard deviation of source mode
    real(r8) mss_frc_srcx(plond,dst_src_nbr) ! [frc] Mass fraction of source mode
    integer::lon_idx ! [idx] Counting index for lon

  ! Properties of source soil distribution
  ! Time/space invariant characteristics of source size distribution:
  ! dmt_vma_src ! [m] Mass median diameter
  ! gsd_anl_src ! [frc] Geometric standard deviation
  ! mss_frc_src ! [frc] Mass fraction
  ! Pick one of the following sets of observed distributions or create your own
  ! Currently initialized in BXM:aer(), CCM:control/ccm3(), MATCH:src/main()

#if 0
  ! Shettle's "Background Desert Model" modes
  ! She84 p. 75 Table 1 with gsd = 2.0 rather than 3.2 modification from SBG98 p. 10581 Table 1
    do lon_idx=1,plon
       mss_frc_srcx(lon_idx,:)=(/2.6e-6_r8,0.781_r8,0.219_r8/) ! [frc] Mass fraction She84 p. 75 Table 1
    enddo
    dmt_vma_srcx(:) = (/ 0.0111e-6,  2.524e-6, 42.10e-6 /)  ! [m] Mass median diameter She84 p. 75 Table 1
    gsd_anl_srcx(:) = (/ 1.89     ,  2.0     ,  2.13    /)  ! [frc] Geometric standard deviation She84 p. 75 Table 1
#endif

#ifndef AlG01
! 20020316: Switch from She84 to Dal87 as default
  ! D'Almeida's (1987) "Background" modes
  ! These modes also summarized in BSM96 p. 73 Table 2 
    do lon_idx=1,plon
       mss_frc_srcx(lon_idx,:)=(/0.036_r8,0.957_r8,0.007_r8/) ! [frc] Mass fraction BSM96 p. 73 Table 2
    enddo
    dmt_vma_srcx(:)=(/ 0.832e-6 ,  4.82e-6 , 19.38e-6 /) ! [m] Mass median diameter BSM96 p. 73 Table 2
    gsd_anl_srcx(:)=(/ 2.10     ,  1.90    ,  1.60    /) ! [frc] Geometric standard deviation BSM96 p. 73 Table 2
#endif /* Dal87 */

#ifdef AlG01
  ! Alfaro and Gomes production modes AlG01 p. 18076 Table 1
  ! Mass fraction taken from BSM96 for lack of better alternative
    do lon_idx=1,plon
       mss_frc_srcx(lon_idx,:)=(/0.036_r8,0.957_r8,0.007_r8/) ! [frc] Mass fraction BSM96 p. 73 Table 2 (ad hoc)
    enddo
    dmt_vma_srcx(:) = (/ 1.5e-6   ,  6.7e-6  , 14.2e-6  /) ! [m] Mass median diameter AlG01 p. 18076 Table 1
    gsd_anl_srcx(:) = (/ 1.7      ,  1.6     ,  1.5     /) ! [frc] Geometric standard deviation AlG01 p. 18076 Table 1
#endif /* AlG01 */
  
#if 0
  ! Patterson and Gillette
  ! PaG77 p. 2080 Table 1 (converted from radius to diameter)
  ! Source soils assumed to have three modes dubbed, in increasing size order, modes C, A, and B
  ! Mode A is mineral dust transport mode, seen in source regions and downwind
  ! Mode B is seen in the source soil itself, and in the atmosphere during dust events
  ! Mode C is seen most everywhere, but does not usually correlate with local dust amount
  ! Mode C is usually a global, aged, background, anthropogenic aerosol, typically rich in sulfate and black carbon
  ! Sometimes, however, mode C has a mineral dust component
  ! Modes C and B are averages from PaG77 Table 1 p. 2080
  ! Mode A is based on summary recommendation: rds_mdn_sfc = 1.5 and gsd = 2.2
  ! mss_frc_src is unavailable in PaG77 because N_0 is not given
  ! Use mss_frc_src from BSM96 p. 73 Table 2 instead
  ! PaG77 nomenclature:  / Mode C   ,  Mode A  ,  Mode B  /
    do lon_idx=1,plon
       mss_frc_srcx(lon_idx,:)=(/0.036_r8,0.957_r8,0.007_r8/) ! [frc] Mass fraction BSM96 p. 73 Table 2 (ad hoc)
    enddo
    dmt_vma_srcx(:) = (/ 0.27e-6  ,  5.6e-6  ,  57.6e-6 /) ! [m] Mass median diameters PaG77 p. 2080 Table 1 
    gsd_anl_srcx(:) = (/ 1.88     ,  2.2     ,  1.62    /) ! [frc] Geometric standard deviation PaG77 p. 2080 Table 1

#endif /* PaG77 */

    mss_frc_src(:,:)=mss_frc_srcx(:,:) ! [frc] Mass fraction of source mode
    dmt_vma_src(:)=dmt_vma_srcx(:) ! [m] Mass median diameter of source mode
    gsd_anl_src(:)=gsd_anl_srcx(:) ! [frc] geometric standard deviation of source mode

    return
    end subroutine dst_psd_src_ini
  
  subroutine asp_rat_lps_set( & ! [fnc] Set ellipsoidal aspect ratio
       asp_rat_lps) ! [frc] Ellipsoidal aspect ratio
    real(r8),intent(in)::asp_rat_lps ! [frc] Ellipsoidal aspect ratio
    asp_rat_lps_dfl=asp_rat_lps ! [frc] Ellipsoidal aspect ratio
  end subroutine asp_rat_lps_set ! [fnc] Set ellipsoidal aspect ratio

  subroutine dst_psd_ini() ! [fnc] Initialize particle size distributions
    ! Purpose: Get dust size grid from data file and set other
    ! important statistics of size distributions
    ! dst_psd_ini() is called by CCM:physics/inti(), MATCH:inirun()
    use dpsdryutl,only:stk_crc_get ! [mdl] Dry deposition utilities
    use dstaer ! [mdl] Aerosol microphysical properties
    use dstctl,only:dst_bnr ! [mdl] Control variables, routines
    use dstgrd,only:dst_nbr ! [mdl] Dust grid sizes
    use psdlgn          ! [mdl] Lognormal particle size distributions
    implicit none
    ! Parameters
    integer,parameter::grd_typ_lgr=0 ! [enm] Logarithmic grid
    integer,parameter::grd_typ_lnr=1 ! [enm] Linear grid
    ! Size range to consider--computed bins in dmt_grd will fill this size range
    ! Range set here has no effect if preset dmt_grd configuration is used
    !real(r8),parameter::dmt_min_min=0.1e-6; ! [m] Minimum grid diameter
    !real(r8),parameter::dmt_max_max=100.0e-6; ! [m] Maximum grid diameter
    real(r8),parameter::dmt_min_min=0.06e-6; ! [m] Minimum grid diameter
    real(r8),parameter::dmt_max_max=50.0e-6; ! [m] Maximum grid diameter
    ! Preset boundaries for common bin configurations
    ! Values are used when a preset configurations exists for dst_nbr bins
    ! Bin boundaries at 1.0, 2.5, and 10.0 microns are nice for comparison to observations
    real(r8),dimension(5),parameter::dmt_grd_04=(/ 0.1e-6 ,  1.0e-6,  2.5e-6,  5.0e-6, 10.0e-6 /) ! [m] Particle diameter grid
    ! real(r8),dimension(5),parameter::dmt_grd_04=(/ 0.1e-6 ,  0.5e-6,  1.0e-6,  5.0e-6, 10.0e-6 /) ! [m] Particle diameter grid
    
    ! Geometric standard deviation (GSD)
    ! She84 p. 75 Table 1, BSM96 p. 73 Table 2, PaG77 p. 2080 Table 1, SBG98 p. 10581 Table 1
    ! She84, BSM96, PaG77, and SBG98 show GSD is always near 2.0
    ! SBG98 shows GSD >> 2.0 is unsupportable for long range transport
    ! SBG98 shows GSD =~ 2.0 is optimal for long range transport mode
    real(r8),parameter::gsd_anl_dfl=2.0 ! [frc] Geometric standard deviation

    ! real(r8),parameter::gsd_anl_dfl=2.0       ! [frc] Geometric standard deviation SBG98 p. 75 Table 1
    ! real(r8),parameter::gsd_anl_dfl=1.90       ! [frc] Geometric standard deviation BSM96 p. 73 Table 2
    ! real(r8),parameter::gsd_anl_dfl=2.2       ! [frc] Geometric standard deviation PaG77 p. 2080 Table 1
    
    ! Aerosol density 
    ! PaG77 p. 2076 say dns_aer = 2500 kg m-3, DKS91 p. 118 say dns_aer = 1600 kg m-3
    real(r8),parameter::dns_aer_dfl=2.5e+3 ! [kg m-3] Aerosol density
    
    ! Commons
    ! Input
    ! Output
    ! Local
    character(80) fl_out      ! [sng] Name of netCDF output file
    integer lon_idx           ! [idx] Counting index for lon
    integer src_idx           ! [idx] Counting index for src
    integer snk_idx           ! [idx] Counting index for snk
    integer idx               ! [idx] Counting index
    integer sz_grd_typ        ! [enm] Size grid type
    real(r8) dmt_ctr(dst_nbr)     ! [m] Diameter at bin center
    real(r8) dmt_dlt(dst_nbr)     ! [m] Width of size bin
    ! Main Code
    
    ! Set grid type
    sz_grd_typ=grd_typ_lgr    ! [enm] Size grid type
    asp_rat_lps(:)=asp_rat_lps_dfl ! [frc] Ellipsoidal aspect ratio
    write (6,'(a,f6.4)') 'asp_rat_lps_dfl = ',asp_rat_lps_dfl
    
    ! Create size grid
    call grd_mk(sz_grd_typ,   & ! I
         dmt_min_min,dmt_max_max,dst_nbr, & ! I
         dmt_ctr,dmt_dlt,dmt_min,dmt_max,dmt_grd) ! O
    
    ! Override automatic grid with preset grid if available
    if (dst_nbr == 4) then
       do idx=1,dst_nbr+1     ! NB: dst_nbr+1
          dmt_grd(idx)=dmt_grd_04(idx) ! [m] Particle diameter grid
       end do                  ! end loop over grd
       do idx=1,dst_nbr
          dmt_min(idx)=dmt_grd_04(idx) ! [m] Maximum diameter in bin
          dmt_max(idx)=dmt_grd_04(idx+1) ! [m] Minimum diameter in bin
          dmt_ctr(idx)=0.5_r8*(dmt_min(idx)+dmt_max(idx)) ! [m] Diameter at bin center
          dmt_dlt(idx)=dmt_max(idx)-dmt_min(idx) ! [m] Width of size bin
       end do                  ! end loop over grd
    endif                     ! endif dst_nbr == 4
    
    ! Bin physical properties
    do idx=1,dst_nbr
       gsd_anl(idx)=gsd_anl_dfl ! [frc] Geometric standard deviation
       dns_aer(idx)=dns_aer_dfl ! [kg m-3] Aerosol density
    end do                     ! end loop over grd
    
    ! Set a fundamental statistic for each bin
    ! For a large number of bins, statistic is irrelevent since distribution is resolved anyway
    ! Choosing "right" statistic is key to value of method for relatively small number of bins (less than about 10) 
    ! Setting bin center equal to number median diameter may be too simplistic
    ! Setting dmt_vma to observed size distribution improves quadrature points when modeled distribution approaches observations
    ! Fit to long-lived transport mode in order to capture long-range effects better
    do idx=1,dst_nbr            
       dmt_vma(idx)=2.524e-6  ! [m] Mass median diameter analytic She84 p. 75 Table 1
       ! dmt_vma(idx)=4.82e-6  ! [m] Mass median diameter analytic BSM96 p. 73 Table 2
       ! dmt_vma(idx)=5.6e-6  ! [m] Mass median diameter analytic PaG77 p. 2080 Table 1
       ! dmt_vma(idx)=1.5e-6  ! [m] Mass median diameter analytic AlG01 p. 18076 Table 1
       ! rds_nma = 0.5_r8*dmt_nma is dst_a input to mie() program
       ! dmt_nma(idx)=dmt_ctr(idx) ! [m] Number median particle size
    end do                     ! end loop over grd
    
    ! Compute analytic size statistics
    
    ! Convert mass median diameter to number median diameter
    call vma2nma(dst_nbr,dmt_vma,gsd_anl, & ! I
         dmt_nma)             ! O [m] Number median particle size
    
    ! ! Convert number median diameter to mass median diameter
    ! call nma2vma(dst_nbr,dmt_nma,gsd_anl, & ! I
    !      dmt_vma)             ! O [m] Mass median diameter
    
    ! Compute resolved size statistics for each size distribution
    do idx=1,dst_nbr
       call dst_sz_rsl( &
            asp_rat_lps(idx), & ! I
            dmt_min(idx),dmt_max(idx), & ! I
            dns_aer(idx),dmt_nma(idx),gsd_anl(idx), & ! I
            dmt_nwr(idx),dmt_swr(idx),dmt_vwr(idx), & ! O
            dmt_nmr(idx),dmt_smr(idx),dmt_vmr(idx), & ! O
            nbr_spc_rsl(idx),sfc_spc_rsl(idx),vlm_spc_rsl(idx) & ! O
            )
    end do                     ! end loop over sz
    
    ! Source/sink fractional overlap factors 
    call ovr_src_snk_frc_get(dst_src_nbr,dmt_vma_src,gsd_anl_src, & ! I
         dst_nbr,dmt_min,dmt_max, & ! I
         ovr_src_snk_frc)     ! O
    
    !APRIL 2002 THIS FOLLOWING CODE IS ALSO DONE IN DSTMBL FOR ALG01 PARAMETERIZATION
    !WE NEED TO COMPUTE mss_frc_src AT RUN TIME BECAUSE IT IS f(u*).
    !THE OVERLAP IS STILL COMPUTED HERE BECAUSE FOR THE NON DYNAMIC SIZE
    !DISTRIBUTIONS, THEY NEVER CHANGES AFTER THIS POINT 
    !AS OF APRIL 2002 THEY ARE FUNCTIONS OF LON TO BE CONSISTENT WITH AlG01 PARAMETERIZATION
    ! 
    ! Multiply ovr_src_snk_frc(src_idx,*) by mss_frc(src_idx) to obtain
    ! absolute mass fraction mapping from source dists. to sink bins
    ! \sum_{i=1}^{i=src_nbr} \sum_{j=1}^{j=snk_nbr} ovr_src_snk_mss(i,j) <= 1
    ! with equality holding iff sources and sinks are completely overlapped
    ! Sum of ovr_src_snk_mss(*,*) is fraction of vertical dust flux which is transported
    ovr_src_snk_mss_ttl(:)=0.0_r8   ! [frc]
    mss_frc_trn_dst_src(:,:)=0.0_r8 ! [frc] Fraction of transported dust mass at source
    do snk_idx=1,dst_nbr
       do src_idx=1,dst_src_nbr
          do lon_idx=1,plon
          ovr_src_snk_mss(lon_idx,src_idx,snk_idx)= & ! [frc]
               ovr_src_snk_frc(src_idx,snk_idx)*mss_frc_src(lon_idx,src_idx) ! [frc]
          mss_frc_trn_dst_src(lon_idx,snk_idx)= & ! [frc] Fraction of transported dust mass at source
               mss_frc_trn_dst_src(lon_idx,snk_idx)+ovr_src_snk_mss(lon_idx,src_idx,snk_idx)
          ovr_src_snk_mss_ttl(lon_idx)= & ! [frc]
               ovr_src_snk_mss_ttl(lon_idx)+ovr_src_snk_mss(lon_idx,src_idx,snk_idx)
          enddo                ! end loo  over lon
       end do                  ! end loop over src
    end do                     ! end loop over snk
   ! Convert fraction of mobilized mass to fraction of transported mass
    do lon_idx=1,plon
       mss_frc_trn_dst_src(lon_idx,:)= & ! [frc] Fraction of transported dust mass at source
            mss_frc_trn_dst_src(lon_idx,:)/ovr_src_snk_mss_ttl(lon_idx)
    enddo
    ! Corrections from Stokes' settling velocities to terminal velocities
    call stk_crc_get( &
         asp_rat_lps,         & ! I [frc] Ellipsoidal aspect ratio
         dmt_vwr,             & ! I [m] Mass-weighted mean diameter resolved
         dns_aer,             & ! I [kg m-3] Aerosol density
         stk_crc,             & ! O [frc] Correction to Stokes settling velocity
         dst_nbr)             ! I [nbr] Number of dust constituents
    
    ! Particles are ellipsoidal when asp_rat_lps != 1.0
    ! In this case, diameter does not suffice to characterize particle dimension
    ! First, dmt_ctr remains the canonical particle size and so must be assigned a physically meaningful definition
    ! Options are to set dmt_ctr to one of 
    ! dmt_stk: Diameter of sphere with same terminal settling velocity and density
    ! dmt_aer: Diameter of sphere with same terminal settling velocity but unit density
    ! dmt_eqv_sfc: Diameter of sphere with same surface area
    ! dmt_eqv_vlm: Diameter of sphere with same surface volume
    ! dmt_mjr: Major axis of ellipsoid is diameter along "a" axis, i.e., 2*a
    ! dmt_mnr: Minor axis of ellipsoid is diameter along "b" axis, i.e., 2*b
    ! We set dmt_ctr equal to diameter of sphere with same surface area as in Gin03 p. 2 (10)
    ! This is sometimes called the surface equivalent diameter
    do idx=1,dst_nbr
       if (asp_rat_lps(idx) == 1.0_r8) then
          dmt_mjr(idx)=dmt_ctr(idx) ! [m] Major axis of ellipsoid
          dmt_mnr(idx)=dmt_ctr(idx) ! [m] Minor axis of ellipsoid
       else
          ! Minor axis of ellipsoid is diameter along "b" axis, i.e., 2*b
          dmt_mnr(idx)=dmt_ctr(idx)/psi_lps_fst_scl(asp_rat_lps(idx)) ! [m] Minor axis of ellipsoid Gin03 p.2 (10)
          ! Major axis of ellipsoid is diameter along "a" axis, i.e., 2*a
          dmt_mjr(idx)=dmt_mnr(idx)*asp_rat_lps(idx) ! [m] Major axis of ellipsoid
       endif ! endif ellipsoid
       ! Stokes diameter is diameter of sphere of same density with same terminal settling velocity
       ! Stokes diameter equals particle diameter when particle is spherical
       ! fxm: Now that ellipsoids are allowed, compute Stokes diameter correctly!
       ! dmt_stk(idx)=dmt_ctr(idx) ! [m] Stokes diameter
    end do ! end loop over sz

    ! Size distribution diagnostics
    ! Convert number median diameter to number mean diameter
    call nma2naa(dst_nbr,dmt_nma,gsd_anl, & ! I
         dmt_naa)             ! O [m] Number mean diameter analytic
    ! Convert number median diameter to surface median diameter
    call nma2sma(dst_nbr,dmt_nma,gsd_anl, & ! I
         dmt_sma)             ! O [m] Surface median diameter analytic
    ! Convert number median diameter to surface area weighted mean diameter
    call nma2swa(dst_nbr,dmt_nma,gsd_anl, & ! I
         dmt_swa)             ! O [m] Surface weighted mean diameter analytic
    ! Convert number median diameter to mass weighted mean diameter
    call nma2vwa(dst_nbr,dmt_nma,gsd_anl, & ! I
         dmt_vwa)             ! O [m] Mass-weighted mean diameter analytic
    
    ! Display banner
    call dst_bnr()
    
    ! Update netCDF file
    fl_out='aer.nc' ! [sng] Name of netCDF output file
    call aersz2nc( &
         fl_out) ! [sng] Name of netCDF output file
    
    write (6,'(a,f4.3)') 'Total transported mass fraction of dust flux = ',ovr_src_snk_mss_ttl
    write (6,'(a)') 'Transported size classes:'
    write (6,'(a3,7(a10,1x))') 'idx','dmt_min','dmt_max','dmt_ctr','dmt_nmr','dmt_nwr','dmt_vmr','dmt_vwr'
    write (6,'(a3,7(a10,1x))') '','um','um','um','um','um','um','um'
    do idx=1,dst_nbr
       write (6,'(i3,7(es10.3,1x))') idx, &
            dmt_min(idx)*1.0e6_r8,dmt_max(idx)*1.0e6_r8,dmt_ctr(idx)*1.0e6_r8, &
            dmt_nmr(idx)*1.0e6_r8,dmt_nwr(idx)*1.0e6_r8, &
            dmt_vmr(idx)*1.0e6_r8,dmt_vwr(idx)*1.0e6_r8
    end do                     ! end loop over cst
    
    do lon_idx=1,plon
    write (6,'(a)') 'Size distribution information:'
    write (6,'(a3,5(a10,1x))') 'idx','dmt_nma','dmt_vma','gsd','dns','mss_frc'
    write (6,'(a3,5(a10,1x))') '','um','um','frc','kg m-3','frc'
    do idx=1,dst_nbr
       write (6,'(i3,5(es10.3,1x))') idx, &
            dmt_nma(idx)*1.0e6_r8,dmt_vma(idx)*1.0e6_r8, &
            gsd_anl(idx),dns_aer(idx),mss_frc_trn_dst_src(lon_idx,idx)
    end do                     ! end loop over cst
    end do ! end loop over longtitude
    
    write (6,'(a)') 'Aspherical particle information:'
    write (6,'(a3,4(a10,1x))') 'idx','asp_rat','dmt_mnr','dmt_mjr','dmt_ctr'
    write (6,'(a3,4(a10,1x))') '','','um','um','um'
    do idx=1,dst_nbr
       write (6,'(i3,4(es10.3,1x))') idx, &
            asp_rat_lps(idx), &
            dmt_mnr(idx)*1.0e6_r8,dmt_mjr(idx)*1.0e6_r8, &
            dmt_ctr(idx)*1.0e6_r8
    end do                     ! end loop over cst
    
    write (6,'(a)') 'Chemistry information:'
    write (6,'(a3,3(a10,1x))') 'idx','stk_crc','nbr_spc','sfc_spc'
    write (6,'(a3,3(a10,1x))') '','frc','# kg-1','m2 kg-1'
    do idx=1,dst_nbr
       write (6,'(i3,3(es10.3,1x))') idx, &
            stk_crc(idx),nbr_spc_rsl(idx),sfc_spc_rsl(idx)
    end do                     ! end loop over cst
    
    write (6,'(a)') '------------------------------------------------------------'
    
    ! Sanity checks
    do idx=1,dst_nbr
       if (dmt_min(idx) > dmt_nmr(idx).or.dmt_max(idx) < dmt_nmr(idx)) stop &
            'dst: dst_psd_ini() reports dmt_min(idx) > dmt_nmr(idx).or.dmt_max(idx) < dmt_nmr(idx)'
       if (dmt_min(idx) > dmt_vwr(idx).or.dmt_max(idx) < dmt_vwr(idx)) stop &
            'dst: dst_psd_ini() reports dmt_min(idx) > dmt_vwr(idx).or.dmt_max(idx) < dmt_vwr(idx)'
       if (gsd_anl(idx) < 1.0.or.gsd_anl(idx) > 3.0_r8) stop &
            'dst: dst_psd_ini() reports gsd_anl(idx) < 1.0.or.gsd_anl(idx) > 3.0'
       if (dns_aer(idx) < 500.0.or.dns_aer(idx) > 3000.0_r8) stop &
            'dst: dst_psd_ini() reports dns_aer(idx) < 500.0.or.dns_aer(idx) > 3000.0'
    end do                     ! end loop over cst
    
    return
  end subroutine dst_psd_ini                       ! end dst_psd_ini()
  
  subroutine dst_sz_rsl( &
       asp_rat_lps, & ! I
       dmt_min,dmt_max, & ! I
       dns_aer,dmt_nma,gsd, & ! I
       dmt_nwr,dmt_swr,dmt_vwr, & ! O
       dmt_nmr,dmt_smr,dmt_vmr, & ! O
       nbr_spc_rsl,sfc_spc_rsl,vlm_spc_rsl & ! O
       )
    ! Compute resolved size statistics for a single size distribution
    use dstcst ! [mdl] Physical constants for dust routines
    use dstgrd ! [mdl] Dust grid sizes
    use psdlgn ! [mdl] Lognormal particle size distributions
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    implicit none
    ! Input
    real(r8),intent(in)::asp_rat_lps ! [frc] Ellipsoidal aspect ratio
    real(r8),intent(in)::dns_aer ! [kg m-3] Particle density
    real(r8),intent(in)::dmt_nma ! [m] Number median particle size
    real(r8),intent(in)::dmt_max ! [m] Size grid maximum
    real(r8),intent(in)::dmt_min ! [m] Size grid minimum
    real(r8),intent(in)::gsd ! [frc] Geometric standard deviation
    ! Output
    real(r8),intent(out)::dmt_nmr ! [m] Number median diameter resolved
    real(r8),intent(out)::dmt_smr ! [m] Surface area median diameter resolved
    real(r8),intent(out)::dmt_vmr ! [m] Mass median diameter resolved
    real(r8),intent(out)::dmt_nwr ! [m] Number weighted diameter resolved
    real(r8),intent(out)::dmt_swr ! [m] Surface area weighted diameter resolved
    real(r8),intent(out)::dmt_vwr ! [m] Mass weighted diameter resolved
    real(r8),intent(out)::nbr_spc_rsl ! [# kg-1] Specific concentration resolved 
    real(r8),intent(out)::sfc_spc_rsl ! [m2 kg-1] Specific surface area resolved 
    real(r8),intent(out)::vlm_spc_rsl ! [m3 kg-1] Specific volume resolved 
    ! Local
    real(r8) pi       ! [frc] 3
    integer sz_grd_typ        ! Size grid type
    integer grd_typ_lgr       ! Logarithmic
    integer grd_typ_lnr       ! Linear
    integer idx               ! [idx] Counting index
    real(r8) cnc(sgs_nbr)         ! [# m-3] Lognormal concentration at sz_ctr
    real(r8) cnc_nbr_rsl          ! [# m-3] Number concentration resolved
    real(r8) dst(sgs_nbr)         ! [# m-3 m-1] Lognormal distribution at sz_ctr
    real(r8) sz_ctr(sgs_nbr)      ! [m] Size Bin centers
    real(r8) sz_dlt(sgs_nbr)      ! [m] Size Bin widths
    real(r8) sz_grd(sgs_nbr+1)    ! [m] Size Grid
    real(r8) sz_max(sgs_nbr)      ! [m] Size Bin maxima
    real(r8) sz_min(sgs_nbr)      ! [m] Size Bin minima
    real(r8) vlm(sgs_nbr)         ! [m3] Volume
    real(r8) vlm_rsl              ! [m3 m-3] Volume concentration resolved
    real(r8) sfc(sgs_nbr)         ! [m2] Surface area
    real(r8) sfc_rsl              ! [m2 m-3] Surface area concentration resolved
    real(r8) nbr_prt_rsl(sgs_nbr) ! [# m-3] Number concentration of smaller particles
    real(r8) nbr_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of number concentration from smaller particles
    real(r8) sfc_prt_rsl(sgs_nbr) ! [m2 m-3] Surface area concentration of smaller particles
    real(r8) sfc_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of surface area concentration from smaller particles
    real(r8) vlm_prt_rsl(sgs_nbr) ! [m3 m-3] Volume concentration of smaller particles
    real(r8) vlm_prt_rsl_frc(sgs_nbr) ! [frc] Fraction of volume concentration from smaller particles
    
    ! Main Code
    ! Initialize
    pi=4.0*atan(1.0_DBLKIND)        ! [frc] 3
    ! Set grid types
    grd_typ_lgr=0             ! Logarithmic
    grd_typ_lnr=1             ! Linear
    sz_grd_typ=grd_typ_lgr
    
    ! Create size grid
    call grd_mk(sz_grd_typ, & ! I
         dmt_min,dmt_max,sgs_nbr, & ! I
         sz_ctr,sz_dlt,sz_min,sz_max,sz_grd) ! O
    
    ! Evaluate lognormal distribution for these sizes
    call lgn_evl(sgs_nbr,sz_ctr,dmt_nma,gsd, & ! I
         dst) ! O
    
    ! Integrate moments of size distribution
    cnc_nbr_rsl=0.0_r8 ! [# m-3] Number concentration resolved
    sfc_rsl=0.0_r8 ! [m2 m-3] Surface area concentration resolved
    vlm_rsl=0.0_r8 ! [m3 m-3] Volume concentration resolved
    dmt_nwr=0.0_r8 ! [m] Number weighted diameter resolved
    dmt_swr=0.0_r8 ! [m] Surface area weighted diameter resolved
    dmt_vwr=0.0_r8 ! [m] Mass weighted diameter resolved
    do idx=1,sgs_nbr
       cnc(idx)=dst(idx)*sz_dlt(idx) ! [# m-3] Number concentration
       cnc_nbr_rsl=cnc_nbr_rsl+cnc(idx) ! [# m-3] Number concentration resolved
       nbr_prt_rsl(idx)=cnc_nbr_rsl ! [# m-3] Number concentration of smaller particles
       sfc(idx)=pi*sz_ctr(idx)*sz_ctr(idx) ! [m2] Surface area 
       sfc_rsl=sfc_rsl+sfc(idx)*cnc(idx) ! [m2 m-3] Surface area distribution resolved so far
       sfc_prt_rsl(idx)=sfc_rsl ! [m2 m-3] Surface area of smaller particles
       vlm(idx)=(1.0/6.0_r8)*pi*(sz_ctr(idx)**3.0_r8) ! [m3] Volume 
       vlm_rsl=vlm_rsl+vlm(idx)*cnc(idx) ! [m3 m-3] Volume distribution resolved so far
       vlm_prt_rsl(idx)=vlm_rsl ! [m3 m-3] Volume of smaller particles
       dmt_nwr=dmt_nwr+sz_ctr(idx)*cnc(idx) ! [m] Number weighted diameter resolved
       dmt_swr=dmt_swr+sz_ctr(idx)*sfc(idx)*cnc(idx) ! [m] Surface area weighted diameter resolved
       dmt_vwr=dmt_vwr+sz_ctr(idx)*vlm(idx)*cnc(idx) ! [m] Mass weighted diameter resolved
    end do                     ! end loop over size
    
    ! Normalize moment weighted sizes by moment
    dmt_nwr=dmt_nwr/cnc_nbr_rsl ! [m] Number weighted diameter resolved
    dmt_swr=dmt_swr/sfc_rsl   ! [m] Surface area weighted diameter resolved
    dmt_vwr=dmt_vwr/vlm_rsl   ! [m] Mass weighted diameter resolved
    
    ! Statistics per unit mass
    nbr_spc_rsl=cnc_nbr_rsl/(vlm_rsl*dns_aer) ! [# kg-1] Specific concentration resolved 
    sfc_spc_rsl=sfc_rsl/(vlm_rsl*dns_aer) ! [m2 kg-1] Specific surface area resolved 
    vlm_spc_rsl=vlm_rsl/(vlm_rsl*dns_aer) ! [m3 kg-1] Specific volume resolved 
    
    do idx=1,sgs_nbr
       nbr_prt_rsl_frc(idx)=nbr_prt_rsl(idx)/cnc_nbr_rsl ! [frc] Fraction of number concentration from smaller particles
       sfc_prt_rsl_frc(idx)=sfc_prt_rsl(idx)/sfc_rsl ! [frc] Fraction of surface area concentration from smaller particles
       vlm_prt_rsl_frc(idx)=vlm_prt_rsl(idx)/vlm_rsl ! [frc] Fraction of volume concentration from smaller particles
    end do                     ! end loop over size
    
    ! Find resolved median sizes by interpolating inverted size -> [number,surface,volume] relationships
    if (sgs_nbr > 1) then
       dmt_nmr=ntp_vec_one(sgs_nbr,nbr_prt_rsl_frc,sz_ctr,0.5_r8) ! [m] Number median diameter resolved
       dmt_smr=ntp_vec_one(sgs_nbr,sfc_prt_rsl_frc,sz_ctr,0.5_r8) ! [m] Surface area median diameter resolved
       dmt_vmr=ntp_vec_one(sgs_nbr,vlm_prt_rsl_frc,sz_ctr,0.5_r8) ! [m] Volume median diameter resolved
    else
       dmt_nmr=sz_ctr(1)      ! [m] Number median diameter resolved
       dmt_smr=sz_ctr(1)      ! [m] Surface area median diameter resolved
       dmt_vmr=sz_ctr(1)      ! [m] Volume median diameter resolved
    endif                     ! endif
    
    return
  end subroutine dst_sz_rsl
  
  subroutine grd_mk(sz_grd_typ, & ! I
       sz_grd_min,sz_grd_max,sz_nbr, & ! I
       sz_ctr,sz_dlt,sz_min,sz_max,sz_grd) ! O
    ! Purpose: Create a size grid
    ! Routine uses sz nomenclature but works for any (linear) grid
    implicit none
    ! Input
    integer,intent(in)::sz_grd_typ ! [enm] Type of grid
    integer,intent(in)::sz_nbr ! [nbr] Number of bins
    real(r8),intent(in)::sz_grd_min ! Grid minimum
    real(r8),intent(in)::sz_grd_max ! Grid maximum
    ! Output
    real(r8),intent(out)::sz_grd(sz_nbr+1) ! Grid
    real(r8),intent(out)::sz_ctr(sz_nbr) ! Bin centers
    real(r8),intent(out)::sz_dlt(sz_nbr) ! Bin widths
    real(r8),intent(out)::sz_min(sz_nbr) ! Bin minima
    real(r8),intent(out)::sz_max(sz_nbr) ! Bin maxima
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) sz_grd_rsn           ! Resolution for linear grid
    real(r8) max_min_ratio        ! Factor for logarithmic grid
    real(r8) series_ratio         ! Factor for logarithmic grid
    integer grd_typ_lgr       ! Logarithmic
    integer grd_typ_lnr       ! Linear
    ! Main Code
    ! Initialize
    grd_typ_lgr=0             ! Logarithmic
    grd_typ_lnr=1             ! Linear
    
    if (sz_grd_typ == grd_typ_lnr) then
       sz_grd_rsn=(sz_grd_max-sz_grd_min)/sz_nbr
       do idx=1,sz_nbr
          sz_min(idx)=sz_grd_min+(idx-1)*sz_grd_rsn ! [m]
       end do                  ! end loop over grd
    else if (sz_grd_typ == grd_typ_lgr) then
       max_min_ratio=sz_grd_max/sz_grd_min
       series_ratio=max_min_ratio**(1.0/sz_nbr)
       if (sz_grd_min == 0.0_r8) stop 'dst: grd_mk() reports sz_grd_min = 0.0'
       sz_min(1)=sz_grd_min
       do idx=2,sz_nbr        ! Loop starts at 2
          sz_min(idx)=sz_min(idx-1)*series_ratio
       end do                  ! end loop over grd
    else
       stop 'dst: grd_mk() reports unknown grid type'
    endif                     ! endif
    
    ! Derived grid values 
    do idx=1,sz_nbr-1         ! Loop ends at sz_nbr-1
       sz_max(idx)=sz_min(idx+1) ! [m]
    end do                     ! end loop over grd
    sz_max(sz_nbr)=sz_grd_max ! [m]
    do idx=1,sz_nbr
       sz_grd(idx)=sz_min(idx) ! [m]
    end do                     ! end loop over grd
    sz_grd(sz_nbr+1)=sz_max(sz_nbr)
    ! These may be parameters so do not change them
    ! sz_grd_max=sz_max(sz_nbr) ! [m]
    ! sz_grd_min=sz_min(1)      ! [m]
    
    ! Final derived grid values 
    do idx=1,sz_nbr
       sz_ctr(idx)=0.5_r8*(sz_min(idx)+sz_max(idx))
       sz_dlt(idx)=sz_max(idx)-sz_min(idx)
    end do                     ! end loop over grd
    
    return
  end subroutine grd_mk
  
  subroutine aersz2nc(             &
       fl_out               & ! [sng] Name of netCDF output file
       )                    
    ! Purpose: Output time-constant aerosol properties to netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstgrd ! [mdl] Dust grid sizes
    use dstaer ! [mdl] Aerosol microphysical properties
    use pmgrid ! [mdl] Spatial resolution parameters
    use dstctl,only:nstep ! [mdl] Control variables, routines
    use sng_mdl,only:ftn_strlen,ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='aersz2nc' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_out ! [sng] Name of netCDF output file
    ! Output
    ! Local
    integer idx               ! [idx] Counting index
    ! File metadata and dimension IDs
    integer dim_lon_lev_sz(3) ! [enm] Dimension IDs
    integer dim_lon_sz(2)     ! [enm] Dimension IDs
    integer dim_sz_src_sz(2)  ! [enm] Dimension IDs
    integer fll_mode_old      ! Old fill mode
    integer lev_dim_id        ! [enm] Dimension ID for lev
    integer levp_dim_id       ! [enm] Dimension ID for levp
    integer lon_dim_id        ! [enm] Dimension ID for lon
    integer lat_dim_id        ! [enm] Dimension ID for lat
    integer time_dim_id       ! [enm] Dimension ID for time
    integer nc_id             ! File handle
    integer rcd               ! [rcd] Return success code
    integer sz_dim_id         ! [enm] Dimension ID for sz
    integer sz_grd_dim_id     ! [enm] Dimension ID for sz grid
    integer sz_src_dim_id     ! [enm] Dimension ID for sz_src
    ! Variable IDs
    integer dmt_naa_id        ! [enm] Variable ID
    integer dmt_vma_id        ! [enm] Variable ID
    integer dmt_vma_src_id    ! [enm] Variable ID
    integer dmt_nma_id        ! [enm] Variable ID
    integer dmt_sma_id        ! [enm] Variable ID
    integer dns_aer_id        ! [enm] Variable ID
    integer dmt_vwa_id        ! [enm] Variable ID
    integer dmt_swa_id        ! [enm] Variable ID
    integer gsd_anl_id        ! [enm] Variable ID
    integer gsd_anl_src_id    ! [enm] Variable ID
    integer asp_rat_lps_id    ! [enm] Variable ID
    integer dmt_max_id        ! [enm] Variable ID
    integer dmt_min_id        ! [enm] Variable ID
    integer dmt_mjr_id        ! [enm] Variable ID
    integer dmt_mnr_id        ! [enm] Variable ID
    integer dmt_nwr_id        ! [enm] Variable ID
    integer dmt_swr_id        ! [enm] Variable ID
    integer dmt_vwr_id        ! [enm] Variable ID
    integer dmt_nmr_id        ! [enm] Variable ID
    integer dmt_smr_id        ! [enm] Variable ID
    integer dmt_vmr_id        ! [enm] Variable ID
    integer mss_frc_src_id    ! [enm] Variable ID
    integer mss_frc_trn_dst_src_id ! [enm] Variable ID
    integer ovr_src_snk_frc_id ! [enm] Variable ID
    integer ovr_src_snk_mss_id ! [enm] Variable ID
    integer ovr_src_snk_mss_ttl_id ! [enm] Variable ID
    integer nbr_spc_rsl_id    ! [enm] Variable ID
    integer sfc_spc_rsl_id    ! [enm] Variable ID
    integer vlm_spc_rsl_id    ! [enm] Variable ID
    integer sz_grd_id         ! Coordinate ID
    integer sz_id             ! Coordinate ID
    integer sz_src_id         ! Coordinate ID
    ! Variable data
    real(r8) sz(dst_nbr)          ! [m] Coordinate variable
    real(r8) sz_src(dst_src_nbr)  ! [m] Coordinate variable
    real(r8) sz_grd(dst_nbr+1)    ! [m] Coordinate variable
    ! Initialize coordinates
    do idx=1,dst_src_nbr
       sz_src(idx)=dmt_vma_src(idx)
    end do                     ! end loop over sz
    do idx=1,dst_nbr
       sz(idx)=0.5_r8*(dmt_min(idx)+dmt_max(idx))
       sz_grd(idx)=dmt_min(idx)
    end do                     ! end loop over sz
    sz_grd(dst_nbr+1)=dmt_max(dst_nbr)
    ! Begin netCDF output routines
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=nf90_wrp_create(fl_out,nf90_clobber,nc_id,sbr_nm=sbr_nm)
    ! Add global attributes
    ! Define dimension IDs
    rcd=rcd+nf90_def_dim(nc_id,'time',nf90_unlimited,time_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'sz_src',dst_src_nbr,sz_src_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'sz',dst_nbr,sz_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'lat',plat,lat_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'lon',plon,lon_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'lev',plev,lev_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'levp',plevp,levp_dim_id)
    rcd=rcd+nf90_def_dim(nc_id,'sz_grd',dst_nbr+1,sz_grd_dim_id)
    ! Assemble ID and count vectors for each multidimensional combination of dimensions
    dim_lon_sz(1)=lon_dim_id
    dim_lon_sz(2)=sz_dim_id
    dim_lon_lev_sz(1)=lon_dim_id
    dim_lon_lev_sz(2)=lev_dim_id
    dim_lon_lev_sz(3)=sz_dim_id
    dim_sz_src_sz(1)=sz_src_dim_id
    dim_sz_src_sz(2)=sz_dim_id
    ! Variable definitions
    rcd=rcd+nf90_def_var(nc_id,'asp_rat_lps',nf90_float,sz_dim_id,asp_rat_lps_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_max',nf90_float,sz_dim_id,dmt_max_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_min',nf90_float,sz_dim_id,dmt_min_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_mjr',nf90_float,sz_dim_id,dmt_mjr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_mnr',nf90_float,sz_dim_id,dmt_mnr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_naa',nf90_float,sz_dim_id,dmt_naa_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_nma',nf90_float,sz_dim_id,dmt_nma_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_nmr',nf90_float,sz_dim_id,dmt_nmr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_nwr',nf90_float,sz_dim_id,dmt_nwr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_sma',nf90_float,sz_dim_id,dmt_sma_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_smr',nf90_float,sz_dim_id,dmt_smr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_swa',nf90_float,sz_dim_id,dmt_swa_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_swr',nf90_float,sz_dim_id,dmt_swr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_vma',nf90_float,sz_dim_id,dmt_vma_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_vma_src',nf90_float,sz_src_dim_id,dmt_vma_src_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_vmr',nf90_float,sz_dim_id,dmt_vmr_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_vwa',nf90_float,sz_dim_id,dmt_vwa_id)
    rcd=rcd+nf90_def_var(nc_id,'dmt_vwr',nf90_float,sz_dim_id,dmt_vwr_id)
    rcd=rcd+nf90_def_var(nc_id,'dns_aer',nf90_float,sz_dim_id,dns_aer_id)
    rcd=rcd+nf90_def_var(nc_id,'gsd_anl',nf90_float,sz_dim_id,gsd_anl_id)
    rcd=rcd+nf90_def_var(nc_id,'gsd_anl_src',nf90_float,sz_src_dim_id,gsd_anl_src_id)
    rcd=rcd+nf90_def_var(nc_id,'nbr_spc_rsl',nf90_float,sz_dim_id,nbr_spc_rsl_id)
    rcd=rcd+nf90_def_var(nc_id,'ovr_src_snk_frc',nf90_float,dim_sz_src_sz,ovr_src_snk_frc_id)
    rcd=rcd+nf90_def_var(nc_id,'sfc_spc_rsl',nf90_float,sz_dim_id,sfc_spc_rsl_id)
    rcd=rcd+nf90_def_var(nc_id,'sz',nf90_float,sz_dim_id,sz_id)
    rcd=rcd+nf90_def_var(nc_id,'sz_grd',nf90_float,sz_grd_dim_id,sz_grd_id)
    rcd=rcd+nf90_def_var(nc_id,'sz_src',nf90_float,sz_src_dim_id,sz_src_id)
    rcd=rcd+nf90_def_var(nc_id,'vlm_spc_rsl',nf90_float,sz_dim_id,vlm_spc_rsl_id)
    ! Add english text descriptions
    rcd=rcd+nf90_put_att(nc_id,asp_rat_lps_id,'long_name','Ellipsoidal aspect ratio')
    rcd=rcd+nf90_put_att(nc_id,dmt_max_id,'long_name','Maximum diameter in bin')
    rcd=rcd+nf90_put_att(nc_id,dmt_min_id,'long_name','Minimum diameter in bin')
    rcd=rcd+nf90_put_att(nc_id,dmt_mjr_id,'long_name','Major axis of ellipsoid')
    rcd=rcd+nf90_put_att(nc_id,dmt_mnr_id,'long_name','Minor axis of ellipsoid')
    rcd=rcd+nf90_put_att(nc_id,dmt_naa_id,'long_name','Number mean diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_nma_id,'long_name','Number median diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_nmr_id,'long_name','Number median diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dmt_nwr_id,'long_name','Number weighted diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dmt_sma_id,'long_name','Surface median diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_smr_id,'long_name','Surface area median diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dmt_swa_id,'long_name','Surface area weighted mean diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_swr_id,'long_name','Surface area weighted diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dmt_vma_id,'long_name','Mass median diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_vma_src_id,'long_name','Mass median diameter of source distribution')
    rcd=rcd+nf90_put_att(nc_id,dmt_vmr_id,'long_name','Mass median diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dmt_vwa_id,'long_name','Mass weighted mean diameter analytic')
    rcd=rcd+nf90_put_att(nc_id,dmt_vwr_id,'long_name','Mass weighted diameter resolved')
    rcd=rcd+nf90_put_att(nc_id,dns_aer_id,'long_name','Particle density')
    rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'long_name','Geometric standard deviation')
    rcd=rcd+nf90_put_att(nc_id,gsd_anl_src_id,'long_name','Geometric standard deviation of source distribution')
    rcd=rcd+nf90_put_att(nc_id,nbr_spc_rsl_id,'long_name','Specific number concentration resolved')
    rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_frc_id,'long_name','Overlap of src dist. i with sink bin j')
    rcd=rcd+nf90_put_att(nc_id,sfc_spc_rsl_id,'long_name','Specific surface area resolved')
    rcd=rcd+nf90_put_att(nc_id,sz_grd_id,'long_name','Size grid interfaces')
    rcd=rcd+nf90_put_att(nc_id,sz_id,'long_name','Nominal size bin center')
    rcd=rcd+nf90_put_att(nc_id,sz_src_id,'long_name','Mass median diameter of source distribution')
    rcd=rcd+nf90_put_att(nc_id,vlm_spc_rsl_id,'long_name','Specific volume resolved')
    ! Add units
    rcd=rcd+nf90_put_att(nc_id,asp_rat_lps_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_max_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_min_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_mjr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_mnr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_naa_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_nma_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_nmr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_nwr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_sma_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_smr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_swa_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_swr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_vma_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_vma_src_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_vmr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_vwa_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dmt_vwr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,dns_aer_id,'units','kilogram meter-3')
    rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'units','fraction')
    rcd=rcd+nf90_put_att(nc_id,gsd_anl_src_id,'units','fraction')
    rcd=rcd+nf90_put_att(nc_id,nbr_spc_rsl_id,'units','number kilogram-1')
    rcd=rcd+nf90_put_att(nc_id,ovr_src_snk_frc_id,'units','fraction')
    rcd=rcd+nf90_put_att(nc_id,sfc_spc_rsl_id,'units','meter2 kilogram-1')
    rcd=rcd+nf90_put_att(nc_id,sz_grd_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,sz_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,sz_src_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,vlm_spc_rsl_id,'units','meter3 kilogram-1')
    ! End define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    rcd=rcd+nf90_put_var(nc_id,asp_rat_lps_id,asp_rat_lps)
    rcd=rcd+nf90_put_var(nc_id,dmt_max_id,dmt_max)
    rcd=rcd+nf90_put_var(nc_id,dmt_min_id,dmt_min)
    rcd=rcd+nf90_put_var(nc_id,dmt_mjr_id,dmt_mjr)
    rcd=rcd+nf90_put_var(nc_id,dmt_mnr_id,dmt_mnr)
    rcd=rcd+nf90_put_var(nc_id,dmt_naa_id,dmt_naa)
    rcd=rcd+nf90_put_var(nc_id,dmt_nma_id,dmt_nma)
    rcd=rcd+nf90_put_var(nc_id,dmt_nmr_id,dmt_nmr)
    rcd=rcd+nf90_put_var(nc_id,dmt_nwr_id,dmt_nwr)
    rcd=rcd+nf90_put_var(nc_id,dmt_sma_id,dmt_sma)
    rcd=rcd+nf90_put_var(nc_id,dmt_smr_id,dmt_smr)
    rcd=rcd+nf90_put_var(nc_id,dmt_swa_id,dmt_swa)
    rcd=rcd+nf90_put_var(nc_id,dmt_swr_id,dmt_swr)
    rcd=rcd+nf90_put_var(nc_id,dmt_vma_id,dmt_vma)
    rcd=rcd+nf90_put_var(nc_id,dmt_vma_src_id,dmt_vma_src)
    rcd=rcd+nf90_put_var(nc_id,dmt_vmr_id,dmt_vmr)
    rcd=rcd+nf90_put_var(nc_id,dmt_vwa_id,dmt_vwa)
    rcd=rcd+nf90_put_var(nc_id,dmt_vwr_id,dmt_vwr)
    rcd=rcd+nf90_put_var(nc_id,dns_aer_id,dns_aer)
    rcd=rcd+nf90_put_var(nc_id,gsd_anl_id,gsd_anl)
    rcd=rcd+nf90_put_var(nc_id,gsd_anl_src_id,gsd_anl_src)
    rcd=rcd+nf90_put_var(nc_id,nbr_spc_rsl_id,nbr_spc_rsl)
    rcd=rcd+nf90_put_var(nc_id,ovr_src_snk_frc_id,ovr_src_snk_frc)
    rcd=rcd+nf90_put_var(nc_id,sfc_spc_rsl_id,sfc_spc_rsl)
    rcd=rcd+nf90_put_var(nc_id,sz_grd_id,sz_grd)
    rcd=rcd+nf90_put_var(nc_id,sz_id,sz)
    rcd=rcd+nf90_put_var(nc_id,sz_src_id,sz_src)
    rcd=rcd+nf90_put_var(nc_id,vlm_spc_rsl_id,vlm_spc_rsl)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'',sbr_nm=sbr_nm)
    if (nstep == 0) then
       write (6,'(a,a47,1x,a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': Initialized size distribution data archive in',fl_out(1:ftn_strlen(fl_out))
    endif                     ! endif
    return 
  end subroutine aersz2nc                       ! end aersz2nc()
  
end module dstpsd

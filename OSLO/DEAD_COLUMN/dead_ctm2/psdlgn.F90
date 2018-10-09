! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/psdlgn.F90,v 1.1 2003/04/15 14:41:36 alfgr Exp $

! Purpose: Initialize, compute, and manipulate lognormal particle size distributions

! Usage:
! use psdlgn ! [mdl] Lognormal particle size distributions

! NB: Unless otherwise specified, all functions assume standard, untruncated, lognormal distributions 

! Size distribution nomenclature:
! Three letter acronym: 
! First letter [nsvx] = number (nbr), surface area (sfc), volume (vlm), cross-sectional area (xsa)
! Second letter [amw] = average (i.e., mean) (avg), median (mdn), weighted (wgt)
! Third letter [ar] = analytic (anl), resolved (rsl)
! If necessary, a fourth letter option ("o"?) could be used 
! in the second letter slot to denote true "mode" averages 
! Not all 24 permutations are unique (much less useful):
! 1. Number mean diameter equals number weighted diameter: naa=nwa, nar=nwr
! 2. Surface area "mean" diameters equal cross-sectional area "mean" diameters 
! for spheres since factor of 4 difference on top and bottom cancel eachother: 
! saa=xaa,sma=xma,sar=xar,smr=xmr

! naa = nbr avg anl 
! nar = nbr avg rsl
! nma = nbr mdn anl
! nmr = nbr mdn rsl
! nwa = nbr wgt anl
! nwr = nbr wgt rsl
! saa = sfc avg anl
! sar = sfc avg rsl
! sma = sfc mdn anl
! smr = sfc mdn rsl
! swa = sfc wgt anl
! swr = sfc wgt rsl
! vaa = vlm avg anl
! var = vlm avg rsl
! vma = vlm mdn anl
! vmr = vlm mdn rsl
! vwa = vlm wgt anl
! vwr = vlm wgt rsl
! xaa = xsa avg anl
! xar = xsa avg rsl
! xma = xsa mdn anl
! xmr = xsa mdn rsl
! xwa = xsa wgt anl
! xwr = xsa wgt rsl

module psdlgn ! [mdl] Lognormal particle size distributions
  use precision ! [mdl] Precision r8, i8, ...
contains
  
  subroutine lgn_evl(abc_nbr,abc_ctr,mdn,gsd, & ! I
       lgn_dst)             ! O
    ! Purpose: Given the median value and geometric standard deviation
    ! for a single lognormal size distribution, compute and return the value of 
    ! the distribution for each of the input values.
    ! This function works for any lognormal size distribution, thus when input
    ! values are diameters, this function returns n=dN/dD
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer abc_nbr           ! [nbr] Dimension size
    real(r8) mdn              ! Median value of distribution
    real(r8) abc_ctr(abc_nbr) ! Abscissae
    real(r8) gsd              ! [frc] Geometric standard deviation
    ! Output
    real(r8) lgn_dst(abc_nbr) ! Lognormal distribution at abc_ctr
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    real(DBLKIND) pi          ! [frc] 3
    real(DBLKIND) lngsdsqrttwopi_rcp ! Factor in lognormal distribution
    real(DBLKIND) tmp         ! Factor in lognormal distribution
    real(DBLKIND) xpn         ! Factor in lognormal distribution
    ! Main Code
    pi=4.0_DBLKIND*atan(1.0_DBLKIND)        ! [frc] 3
    ln_gsd=log(gsd)
    lngsdsqrttwopi_rcp=1.0_DBLKIND/(ln_gsd*sqrt(2.0_DBLKIND*pi))
    do idx=1,abc_nbr
       tmp=log(abc_ctr(idx)/mdn)/ln_gsd
       xpn=exp(-0.5_DBLKIND*tmp*tmp)
       lgn_dst(idx)=lngsdsqrttwopi_rcp*xpn/abc_ctr(idx)
    end do                    ! end loop over abc
    return
  end subroutine lgn_evl
  
  subroutine naa2nma(aer_nbr,dmt_avg,gsd, & ! I
       dmt_mdn)             ! O
    ! Purpose: Given the mean values and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the median
    ! values of each distribution 
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_avg(aer_nbr) ! Mean value
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_mdn(aer_nbr) ! Median value
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_mdn(idx)=dmt_avg(idx)*exp(-ln_gsd*ln_gsd/2.0_r8)
    end do                     ! end loop over aer
    return
  end subroutine naa2nma
  
  subroutine nma2naa(aer_nbr,dmt_mdn,gsd, & ! I
       dmt_avg)             ! O
    ! Purpose: Given the median values and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the mean
    ! value of each distribution 
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_mdn(aer_nbr) ! Median value
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_avg(aer_nbr) ! Mean value
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_avg(idx)=dmt_mdn(idx)*exp(ln_gsd*ln_gsd/2.0_r8)
    end do                    ! end loop over aer
    return
  end subroutine nma2naa
  
  subroutine nma2swa(aer_nbr,dmt_nma,gsd, & ! I
       dmt_swa)             ! O
    ! Purpose: Given the number median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the 
    ! area-weighted mean size of each distribution
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_swa(aer_nbr) ! [m] Effective (area weighted) mean particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_swa(idx)=dmt_nma(idx)*exp(5.0_r8*ln_gsd*ln_gsd/2.0_r8) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine nma2swa
  
  subroutine sma2nma(aer_nbr,dmt_sma,gsd, & ! I
       dmt_nma)             ! O
    ! Purpose: Given the surface median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the number median size of each distribution 
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_sma(aer_nbr) ! [m] Surface median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_nma(idx)=dmt_sma(idx)*exp(-2.0_r8*ln_gsd*ln_gsd) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine sma2nma
  
  subroutine vma2nma(aer_nbr,dmt_vma,gsd, & ! I
       dmt_nma)             ! O
    ! Purpose: Given the mass median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the number median size of each distribution 
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_vma(aer_nbr) ! [m] Mass median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_nma(idx)=dmt_vma(idx)*exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine vma2nma
  
  subroutine nma2sma(aer_nbr,dmt_nma,gsd, & ! I
       dmt_sma)             ! O
    ! Purpose: Given the number median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the surface median size of each distribution
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_sma(aer_nbr) ! [m] Surface median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_sma(idx)=dmt_nma(idx)*exp(2.0_r8*ln_gsd*ln_gsd) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine nma2sma
  
  subroutine saa2nma(aer_nbr,dmt_das,gsd, & ! I
       dmt_nma)             ! O
    ! Purpose: Given the surface average diameter and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the number median size of each distribution
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_das(aer_nbr) ! [m] Surface average particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_nma(idx)=dmt_das(idx)*exp(-ln_gsd*ln_gsd) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine saa2nma
  
  subroutine vaa2nma(aer_nbr,dmt_dam,gsd, & ! I
       dmt_nma)             ! O
    ! Purpose: Given the mass average diameter and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the number median size of each distribution
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_dam(aer_nbr) ! [m] Mass average particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_nma(idx)=dmt_dam(idx)*exp(-3.0_r8*ln_gsd*ln_gsd/2.0_r8) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine vaa2nma
  
  subroutine nma2vma(aer_nbr,dmt_nma,gsd, & ! I
       dmt_vma)             ! O
    ! Purpose: Given the number median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the mass median (same as volume-median) size of each distribution 
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_vma(aer_nbr) ! [m] Mass median particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_vma(idx)=dmt_nma(idx)*exp(3.0_r8*ln_gsd*ln_gsd) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine nma2vma
  
  subroutine nma2vwa(aer_nbr,dmt_nma,gsd, & ! I
       dmt_vwa)             ! O
    ! Purpose: Given the number median particle sizes and geometric standard deviations
    ! for a number of lognormal size distributions, compute and return the mass-weighted (same as volume-weighted) mean size of each distribution 
    ! Definition of gsd conforms to SeP97
    use precision             ! [mdl] Precision r8, i8, ...
    implicit none
    ! Input
    integer aer_nbr           ! [nbr] Dimension size
    real(r8) dmt_nma(aer_nbr) ! [m] Number median particle size
    real(r8) gsd(aer_nbr)     ! [frc] Geometric standard deviation
    ! Output
    real(r8) dmt_vwa(aer_nbr) ! [m] Mass-weighted mean particle size
    ! Local
    integer idx               ! [idx] Counting index
    real(r8) ln_gsd           ! [frc] ln(gsd)
    ! Main Code
    do idx=1,aer_nbr
       ln_gsd=log(gsd(idx))
       dmt_vwa(idx)=dmt_nma(idx)*exp(7.0_r8*ln_gsd*ln_gsd/2.0_r8) ! [m]
    end do                     ! end loop over aer
    return
  end subroutine nma2vwa
  
  subroutine ovr_src_snk_frc_get(src_nbr,mdn_src,gsd_src, & ! I
       snk_nbr,dmt_min_snk,dmt_max_snk, & ! I
       ovr_src_snk_frc)     ! O
    ! Purpose: Given one set (the "source") of lognormal distributions, 
    ! and one set of bin boundaries (the "sink"), compute and return 
    ! the overlap factors between the source distributions and the sink bins.
    ! The first group, src, contains src_nbr distributions
    ! Each src distribution consists of a specified median diameter and geometric standard deviation
    ! The second group, snk, contains snk_nbr bins 
    ! Each sink bin consists of a specified minimum and maximum diameter
    ! The output is a matrix ovr_src_snk_frc(src_nbr,snk_nbr)
    ! Element ovr_src_snk_frc(i,j) is the fraction of size distribution i in group src 
    ! that overlaps sink bin j
    ! If the i'th source is completely bracketed by sink bins, then
    ! \sum_{j=1}^{j=snk_nbr} ovr_src_snk_frc(i,j) = 1
    ! If the j'th sink is completely bounded by the sources, then
    ! \sum_{i=1}^{i=src_nbr} ovr_src_snk_frc(i,j) = 1
    ! Note that input are generic median diameters of lognormal distributions
    ! Thus if routine is called with number median diameters it will compute
    ! the overlap factors for number concentration, and if it is called with mass
    ! median diameters then it will compute the overlap factors for mass.
    ! When ovr_src_snk_frc_get() is called with mass distributions parameters as arguments,
    ! then ovr_src_snk_frc(src_idx,*) should be multiplied by mss_frc(src_idx) to obtain
    ! the absolute mass fraction mapping from source bins to sink bins.
    use precision ! [mdl] Precision r8, i8, ...
    use erf_mdl,only:erf ! [mdl] Error functions erf(), erfc(), erfcx()
    implicit none
    ! Input
    integer snk_nbr           ! [nbr] Dimension size
    integer src_nbr           ! [nbr] Dimension size
    real(r8) dmt_max_snk(snk_nbr) ! [m] Maximum diameter in bin
    real(r8) dmt_min_snk(snk_nbr) ! [m] Minimum diameter in bin
    real(r8) gsd_src(src_nbr) ! [frc] Geometric standard deviation
    real(r8) mdn_src(src_nbr) ! [m] Mass median particle size
    ! Output
    real(r8) ovr_src_snk_frc(src_nbr,snk_nbr) ! [frc] Fractional overlap of src with snk
    ! Local
    integer src_idx           ! [idx] Counting index for src
    integer snk_idx           ! [idx] Counting index for snk
    real(r8) ln_gsd           ! [frc] ln(gsd)
    real(r8) sqrt2lngsdi      ! [frc] Factor in erf() argument
    real(r8) lndmaxjovrdmdni  ! [frc] Factor in erf() argument
    real(r8) lndminjovrdmdni  ! [frc] Factor in erf() argument
    ! Main Code
    ! Sanity check
    ! 19990913: erf() in SGI /usr/lib64/mips4/libftn.so is bogus
    if (abs(0.8427_r8-erf(1.0_r8))/0.8427_r8 > 0.001_r8) then
       write (6,'(a,f12.10)') 'erf(1.0_r8) = ',erf(1.0_r8)
       stop 'dst: ovr_src_snk_frc_get() reports Error function error'
    endif                     ! endif
    if (erf(0.0_r8) /= 0.0_r8) then
       write (6,'(a,f12.10)') 'erf(0.0_r8) = ',erf(0.0_r8)
       stop 'dst: ovr_src_snk_frc_get() reports Error function error'
    end if                     ! endif
    do src_idx=1,src_nbr
       sqrt2lngsdi=sqrt(2.0_r8)*log(gsd_src(src_idx)) ! [frc]
       do snk_idx=1,snk_nbr
          lndmaxjovrdmdni=log(dmt_max_snk(snk_idx)/mdn_src(src_idx)) ! [frc]
          lndminjovrdmdni=log(dmt_min_snk(snk_idx)/mdn_src(src_idx)) ! [frc]
          ovr_src_snk_frc(src_idx,snk_idx)= & ! [frc]
               0.5_r8*(erf(lndmaxjovrdmdni/sqrt2lngsdi)- &
               erf(lndminjovrdmdni/sqrt2lngsdi))
       end do                  ! end loop over snk
    end do                     ! end loop over src
        
    return
  end subroutine ovr_src_snk_frc_get
      
end module psdlgn

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstmssutl.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $

! Purpose: Mass budget utilities for the dust model

! Usage: 
! use dstmssutl ! [mdl] Mass budget utilities

module dstmssutl ! [mdl] Mass budget utilities
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  
contains
  
  
  !// ------------------------------------------------------------------
  subroutine dst_add_lev(q,q_ttl)
    !// ------------------------------------------------------------------
    ! Purpose: Given a 2-D array of an additive property (e.g., mixing ratio, flux)
    ! which is dimensioned (plev,dst_nbr), dst_add_lev() computes 
    ! and returns the total property (e.g., mixing ratio, flux), obtained 
    ! by simply adding along the (dust) constituent dimension.
    !//
    !// Added for Oslo CTM3.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstgrd ! [mdl] Dust grid sizes
    use pmgrid ! [mdl] Spatial resolution parameters
    use dead_precision ! [mdl] Precision r8, i8, ...
    use vec_mdl,only:vec_set ! [mdl] Vector manipulation, interpolation, rebinning
    implicit none
    ! Input
    real(r8),intent(in) :: q(plev,dst_nbr) ! Total property
    ! Output
    real(r8),intent(out) :: q_ttl(plev)    ! Property for each size class
    ! Local
    integer k                 ! lev index
    integer m                 ! cst index 
    !// ------------------------------------------------------------------
    ! Initialize total property to zero then integrate
    call vec_set(q_ttl,plev,0.0_r8)
    do m=1,dst_nbr
       do k=1,plev
          q_ttl(k) = q_ttl(k)+q(k,m)
       end do                 ! end loop over lev
    end do                    ! end loop over cst
    !// ------------------------------------------------------------------
  end subroutine dst_add_lev                       ! end dst_add_lon_lev()
  !// ------------------------------------------------------------------
 
  !// ------------------------------------------------------------------
  subroutine dst_add_nbr(q,q_ttl)
    !// ------------------------------------------------------------------
    ! Purpose: Given a 1-D array of an additive property (e.g., mixing ratio, flux)
    ! which is dimensioned (plond,dst_nbr), dst_add_lon() computes 
    ! and returns the total property (e.g., mixing ratio, flux), obtained 
    ! by simply adding along the (dust) constituent dimension.
    !//
    !// Added for Oslo CTM3.
    !//
    !// Amund Sovde, October 2009
    !// ------------------------------------------------------------------
    use dstgrd ! [mdl] Dust grid sizes
    !use pmgrid ! [mdl] Spatial resolution parameters
    use dead_precision ! [mdl] Precision r8, i8, ...
    implicit none
    !// ------------------------------------------------------------------
    ! Input
    real(r8),intent(in)::q(dst_nbr)     ! Total property
    ! Output
    real(r8),intent(out)::q_ttl         ! Property for each size class
    ! Local
    integer m                 ! cst index 
    !// ------------------------------------------------------------------
    ! Initialize total property to zero then integrate
    q_ttl = 0.0_r8
    do m=1,dst_nbr
          q_ttl = q_ttl + q(m)
    end do                    ! end loop over cst
    !// ------------------------------------------------------------------
  end subroutine dst_add_nbr                       ! end dst_add_lon()
  !// ------------------------------------------------------------------
  

end module dstmssutl ! [mdl] Mass budget utilities

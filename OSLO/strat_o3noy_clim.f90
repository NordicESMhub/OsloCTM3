!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Stratospheric O3 and NOy when running only tropchem.
!//=========================================================================
module strat_o3noy_clim
  !// ----------------------------------------------------------------------
  !// MODULE: strat_o3noy_clim
  !// DECRIPTION: Routines for setting stratosphere O3 and NOy when
  !//             stratchem is not included.
  !//
  !// This module has 2 purposes:
  !// 1. CTM3 climatology for stratospheric O3 (volume mixing ratio)
  !//    and for HNO3, PAN, HO2NO2, NO3, N2O5, NO and NO2, when
  !//    stratospheric chemistry is not calculated.
  !// 2. Update DO3 used in calculation of J-values, which is done also
  !//    when stratospheric chemistry is used.
  !//
  !// Contains
  !//   subroutine read_o3clim
  !//   subroutine get_strato3noy_clim
  !//   subroutine stratO3_interp
  !//   subroutine update_stratO3
  !//   subroutine stratNOY_interp
  !//   subroutine update_stratNOX
  !//   subroutine update_stratNOX2
  !//   subroutine update_stratNOY
  !//   subroutine update_strato3ctm2
  !//
  !//
  !// Amund Sovde, May - June 2010
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: LPAR, IPAR, JPAR, LOSLOCSTRAT, LOSLOCTROP
  use cmn_ctm, only: JDATE
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// O3 climatology (vmr) interpolated
  real(r8), dimension(LPAR,IPAR,JPAR) :: STRATO3

  !// O3 climatology (vmr) for two months
  real(r8), dimension(2,LPAR,IPAR,JPAR) :: O3CLIM

  !// NOy climatology (vmr)
  integer, parameter :: NXNY=7
  integer,dimension(NXNY), parameter :: NYID=(/4,5,17,41,42,43,44/)
  real(r8), dimension(NXNY,LPAR,IPAR,JPAR) :: STRATNOY
  real(r8), dimension(2,NXNY,LPAR,IPAR,JPAR) :: NOYCLIM


  !// ----------------------------------------------------------------------
  !// All variables are to be saved.
  save
  !// All is private
  private
  !// except
  public get_strato3noy_clim, update_stratO3, update_stratNOX, &
       update_stratNOX2, update_stratNOY
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine read_o3clim(MONTH)
    !// --------------------------------------------------------------------
    !// Read O3 (vmr) climatology from file. Input is T42L60.
    !// Handles collapsed layers, but reads data from non-collapsed CTM levels.
    !//
    !// Amund Sovde, May 2010
    !// --------------------------------------------------------------------
    use cmn_precision, only: r4
    use netcdf
    use cmn_size, only: LPARW
    use cmn_ctm, only: XLMMAP, LMMAP, XDEDG, YDEDG, AREAXY
    use cmn_parameters, only: A0, CPI180
    use regridding, only: E_GRID
    use ncutils, only: handle_err
    use cmn_oslo, only: trsp_idx
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MONTH

    !// File treatment
    character(len=80) :: filename
    real(r4), dimension(:,:,:,:),allocatable :: o3_in
    real(r4), dimension(:,:,:,:,:),allocatable :: noy_in
    real(r8), dimension(:,:),allocatable :: tmp2d
    real(r8), dimension(:),allocatable :: o3_xdedg, o3_ydedg, o3_xybox

    !// netcdf integers
    integer                  :: lon_dim_id   !// Id for longitude dimension
    integer                  :: lon_id       !// Id for variable longitude
    integer                  :: lat_dim_id   !// Id for latitude dimension
    integer                  :: lat_id       !// Id for latitude
    integer                  :: lev_dim_id   !// Id for level dimension
    integer                  :: lev_id       !// Id for level
    integer                  :: time_dim_id  !// Id for time dimension
    integer                  :: time_id      !// Id for time
    integer                  :: field_id     !// Variable id for field
    integer                  :: nlons        !// Longitudes in file
    integer                  :: nlats        !// Latitudes in file
    integer                  :: nlevs        !// Levels in file
    integer                  :: nsteps       !// Timesteps avaiable in file
    integer                  :: status       !// status of process (0=OK)
    integer                  :: ncid         !// file id 
    integer                  :: nnoy !// Number of NOy components on file
    integer                  :: nnoy_dim_id !// Id for time dimension

    !// For interpolation
    logical :: LINTERP
    real(r8), dimension(IPAR, JPAR) :: R8XY
    integer :: mon1, mon2, APAR, I, J, L, LL, NY, LMAP(LPAR+1)
    !// ------------------------------------------------------------------

    if (.not.(trsp_idx(1).gt.0 .and. trsp_idx(4).gt.0 .and. trsp_idx(43).gt.0 .and. &
              trsp_idx(5).gt.0 .and. trsp_idx(17).gt.0 .and.trsp_idx(41).gt.0 .and. &
              trsp_idx(44).gt.0 )) then
       print*,'*** oc_strato3clim.f90: Cannot run tropospheric chemistry,'// &
            'since one component is missing/not transported:'
       print*,'    O3/HNO3/NO',trsp_idx(1).gt.0,trsp_idx(4).gt.0,trsp_idx(43).gt.0
       print*,'    PANX/HO2NO2/NO3',trsp_idx(5).gt.0,trsp_idx(17).gt.0,trsp_idx(41).gt.0
       print*,'    N2O5/NO2',trsp_idx(42).gt.0,trsp_idx(44).gt.0
       stop
    end if
    if (LPARW .ne. 60) then
       print*,'*** oc_strato3clim.f90: Not set up for other vertical '//&
            'resolutions than L60.'
       print*,'    You use:',LPAR
       stop
    end if
      

    !// File to be read is netcdf
    filename = 'Indata_CTM3/ctm2_o3noy_climatology.nc'
    status = nf90_noerr  !Status is 0 and should be kept that way !!
    write(*,'(a)') '* Reading CTM2 produced O3 & NOy monthly climatology'
    write(*,'(a)') '  File: '//trim(filename)

    !// Open the existing file
    status = nf90_open(filename, nf90_nowrite, ncid)
    if (status /= nf90_noerr) call handle_err(status)

    !// Inquire dimension ids
    status = nf90_inq_dimid(ncid,'Latitudes',lat_dim_id)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid(ncid,'Longitudes',lon_dim_id)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid(ncid,'Levels',lev_dim_id)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_dimid(ncid,'Months',time_dim_id)
    if (status /= nf90_noerr) call handle_err(status)

    !// NOy
    status = nf90_inq_dimid(ncid,'NOYcomps',nnoy_dim_id)
    if (status /= nf90_noerr) call handle_err(status)



    !// Inquire dimensions
    status = nf90_Inquire_Dimension(ncid,lat_dim_id,len=nlats)
    if (status /= nf90_noerr) call handle_err(status)
    if (nlats == JPAR) then
       LINTERP = .false.
    else
       LINTERP = .true.
    end if

    status = nf90_Inquire_Dimension(ncid,lon_dim_id,len=nlons)
    if (status /= nf90_noerr) call handle_err(status)
    if (nlons /= IPAR) then
       LINTERP = .true.
    endif

    status = nf90_Inquire_Dimension(ncid,lev_dim_id,len=nlevs)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_Inquire_Dimension(ncid,time_dim_id,len=nsteps)
    if (status /= nf90_noerr) call handle_err(status)
    if (nsteps /= 12) then
       write(*,'(a)') '* Not 12 months of data in '//trim(filename)
       stop
    end if

    !// NOy components
    status = nf90_Inquire_Dimension(ncid,nnoy_dim_id,len=nnoy)
    if (status /= nf90_noerr) call handle_err(status)


    !// Allocate
    allocate( o3_in(nlons,nlats,nlevs,2), noy_in(nlons,nlats,nlevs,nnoy,2), &
         tmp2d(nlons,nlats), o3_xdedg(nlons+1), o3_ydedg(nlats+1), o3_xybox(nlats) )


    !// Get grid size
    status = nf90_inq_varid(ncid,'XDEDG',field_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var( ncid,field_id,o3_xdedg ) !, &

    status = nf90_inq_varid(ncid,'YDEDG',field_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_get_var( ncid,field_id,o3_ydedg ) !, &


    !// Get data for this month and the next
    mon1 = MONTH
    if (MONTH .lt. 12) then
       mon2 = MONTH + 1
    else
       mon2 = 1
    end if


    !// May handle collapsed layers
    !// Number of levels in CTM; less than nlevs when collapsing layers:
    APAR = maxval(LMMAP(1:nlevs))
    !// Each full (non-collapsed) level starting points are given in LMAP
    !// LMAP(1)=1, and if layers 1:3 are collapsed, LMAP(2) = 4.
    do LL = LPAR+1,1,-1
       L   = LMMAP(LL)
       LMAP(L) = LL    !// The start levels (non-collapsed) of each CTM level
    end do


    !// GET AREA of grid boxes =====================================
    do J=1,nlats
       o3_xybox(J) =  A0*A0 * CPI180*(o3_xdedg(2)-o3_xdedg(1)) &
            * (sin(CPI180*o3_ydedg(J+1)) - sin(CPI180*o3_ydedg(J)))
    end do


    !// GET O3 =====================================================

    !// Get variable ID
    status = nf90_inq_varid(ncid,'O3CLIM',field_id)
    if (status /= nf90_noerr) call handle_err(status)
       
    !// Get the variable month 1
    status = nf90_get_var( ncid,field_id,o3_in(:,:,:,1), &
         start=(/1, 1, 1, mon1/), count=(/nlons, nlats, nlevs, 1/) )
    if (status /= nf90_noerr) call handle_err(status)

    !// Get the variable month 2
    status = nf90_get_var( ncid,field_id,o3_in(:,:,:,2), &
         start=(/1, 1, 1, mon2/), count=(/nlons, nlats, nlevs, 1/) )
    if (status /= nf90_noerr) call handle_err(status)

    write(*,'(a)') '  Got variable: O3CLIM'



    !// GET NOy ====================================================

    !// Get variable ID
    status = nf90_inq_varid(ncid,'NOYCLIM',field_id)
    if (status /= nf90_noerr) call handle_err(status)

    !// Get the variable month 1
    status = nf90_get_var( ncid,field_id,noy_in(:,:,:,:,1), &
         start=(/1, 1, 1, 1, mon1/), count=(/nlons, nlats, nlevs, nnoy, 1/) )
    if (status /= nf90_noerr) call handle_err(status)

    !// Get the variable month 2
    status = nf90_get_var( ncid,field_id,noy_in(:,:,:,:,2), &
         start=(/1, 1, 1, 1, mon2/), count=(/nlons, nlats, nlevs, nnoy, 1/) )
    if (status /= nf90_noerr) call handle_err(status)

    write(*,'(a)') '  Got variable: NOYCLIM'

    !// Closing file
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)



    !// Store data =================================================
    !// Initialize
    O3CLIM(:,:,:,:) = 0._r8
    NOYCLIM(:,:,:,:,:) = 0._r8

    !// O3
    if (LINTERP) then
       write(*,'(a)') '  - will interpolate'

       !// CTM layers corresponding to NLEVS
       do L = 1, APAR

          !// Loop over all full levels for each (collapsed) CTM levels
          do LL = LMAP(L),LMAP(L+1)-1
             
             !// Interpolate month 1 (this month: o3_in(:,:,LL,1))
             do J = 1, nlats
               do I = 1, nlons
                 tmp2d(I,J) = o3_in(I,J,LL,1) * o3_xybox(J)
               end do
             end do
             call E_GRID(tmp2d, o3_xdedg, o3_ydedg, nlons, nlats, &
                         R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
             !// Put into CTM array. Make average if collapsed layers.
             do J = 1, JPAR
               do I = 1, IPAR
                 O3CLIM(1,L,I,J) = O3CLIM(1,L,I,J) &
                                   + R8XY(I,J)/AREAXY(I,J)*XLMMAP(LL)
               end do
             end do

             !// Interpolate month 2 (next month: o3_in(:,:,LL,2))
             do J = 1, nlats
               do I = 1, nlons
                 tmp2d(I,J) = o3_in(I,J,LL,2) * o3_xybox(J)
               end do
             end do
             call E_GRID(tmp2d, o3_xdedg, o3_ydedg, nlons, nlats, &
                         R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
             !// Put into CTM array
             do J = 1, JPAR
               do I = 1, IPAR
                 O3CLIM(2,L,I,J) = O3CLIM(2,L,I,J) &
                                   + R8XY(I,J)/AREAXY(I,J)*XLMMAP(LL)
               end do
             end do

          end do !// do LL = LMAP(L),LMAP(L+1)-1
       end do !// do L = 1, APAR
    else
       write(*,'(a)') '  - In current resolution'

       !// This month
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  O3CLIM(1,L,I,J) = O3CLIM(1,L,I,J) + &
                                    o3_in(I,J,LL,1)*XLMMAP(LL)
                end do
             end do
          end do
       end do
       !// Next month
       do J = 1, JPAR
          do I = 1, IPAR
             do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  O3CLIM(2,L,I,J) = O3CLIM(2,L,I,J) + &
                                    o3_in(I,J,L,2)*XLMMAP(LL)
                end do
             end do
          end do
       end do
    end if !// if (LINTERP) then


    !// NOy
    do NY = 1, NNOY

       if (LINTERP) then
          write(*,'(a,i2,a)') '  - NOy ',ny,' interpolated'

          !// CTM layers corresponding to NLEVS
          do L = 1, APAR

             !// Loop over all full levels for each (collapsed) CTM levels
             do LL = LMAP(L),LMAP(L+1)-1
             
                !// Need to multiply by area before E_GRID and divide by
                !// new area afterwards.

                !// Interpolate month 1 (this month: noy_in(:,:,LL,NY,1))
                do J = 1, nlats
                   do I = 1, nlons
                      tmp2d(I,J) = noy_in(I,J,LL,NY,1) * o3_xybox(J)
                   end do
                end do
                call E_GRID(tmp2d, o3_xdedg, o3_ydedg, nlons, nlats, &
                     R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
                !// Put into CTM array. Make average if collapsed layers.
                do J = 1, JPAR
                   do I = 1, IPAR
                      NOYCLIM(1,NY,L,I,J) = NOYCLIM(1,NY,L,I,J) &
                                            + R8XY(I,J)/AREAXY(I,J)*XLMMAP(LL)
                   end do
                end do

                !// Interpolate month 2 (next month: noy_in(:,:,LL,NY,2))
                do J = 1, nlats
                   do I = 1, nlons
                      tmp2d(I,J) = noy_in(I,J,LL,NY,2) * o3_xybox(J)
                   end do
                end do
                call E_GRID(tmp2d, o3_xdedg, o3_ydedg, nlons, nlats, &
                         R8XY, XDEDG, YDEDG, IPAR, JPAR, 1)
                !// Put into CTM array
                do J = 1, JPAR
                   do I = 1, IPAR
                      NOYCLIM(2,NY,L,I,J) = NOYCLIM(2,NY,L,I,J) &
                                            + R8XY(I,J)/AREAXY(I,J)*XLMMAP(LL)
                   end do
                end do
             end do !// do LL = LMAP(L),LMAP(L+1)-1
          end do !// do L = 1, APAR
       else

          write(*,'(a,i2,a)') '  - NOy ',ny,' in current resolution'

          !// This month
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  NOYCLIM(1,NY,L,I,J) = NOYCLIM(1,NY,L,I,J) + &
                                        noy_in(I,J,LL,NY,1)*XLMMAP(LL)
                end do
              end do
            end do
          end do
          !// Next month
          do J = 1, JPAR
            do I = 1, IPAR
              do L = 1, APAR
                !// Loop over all full levels for each (collapsed) CTM levels
                do LL = LMAP(L),LMAP(L+1)-1
                  NOYCLIM(2,NY,L,I,J) = NOYCLIM(2,NY,L,I,J) + &
                                         noy_in(I,J,LL,NY,2)*XLMMAP(LL)
                end do
              end do
            end do
          end do

       end if
    end do !// do NY = 1, NNOY
 
    !// Deallocate
    deallocate(o3_in, noy_in, tmp2d, o3_xdedg, o3_ydedg, o3_xybox)


    !// --------------------------------------------------------------------
  end subroutine read_o3clim
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine get_strato3noy_clim(NDAY,NMET,NOPS,MONTH,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Read a CTM2 climatology of O3, on monthly mean basis.
    !// To be used in stratosphere when running only tropospheric chemistry.
    !//
    !// Amund Sovde, May 2010
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY,NMET,NOPS,MONTH
    logical, intent(in) :: LNEW_MONTH
    !// --------------------------------------------------------------------

    !// Should not be here if stratospheric chemistry is included
    !// and not if tropospheric chemistry is not included.
    if (LOSLOCSTRAT .or. (.not.LOSLOCTROP)) return


    !// Read climatology from file
    if (LNEW_MONTH) call read_o3clim(MONTH)


    !// Interpolate each day at 00UTC
    if (NMET.eq.1 .and. NOPS.eq.1) then
       call stratO3_interp(MONTH)
       call stratNOY_interp(MONTH)
    end if

    !// --------------------------------------------------------------------
  end subroutine get_strato3noy_clim
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
  subroutine stratO3_interp(MONTH)
    !// --------------------------------------------------------------------
    !// Interpolate climatology linearly in time between two monthly means.
    !// Done one time each day, outside parallel region.
    !//
    !// Amund Sovde, May 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MONTH

    !// Locals
    real(r8) :: dtfrac
    integer :: I, J, L
    !// --------------------------------------------------------------------

    !// Interpolate linearly in time between two monthly means,
    !// updated every day. Fraction of first (this) month:
    dtfrac = 1._r8 - real(jdate-1, r8) / DINM(MONTH)

    !// Loop over latitude
    do J = 1, JPAR
      !// Loop over longitude
      do I = 1, IPAR
        do L = 1, LPAR
          !// Set O3 from climatology, linearly interpolated between two
          !// monthly averages.
          STRATO3(L,I,J) = O3CLIM(1,L,I,J) * dtfrac + &
                           O3CLIM(2,L,I,J) * (1._r8 - dtfrac)
        end do
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine stratO3_interp
  !// ----------------------------------------------------------------------





  !// ----------------------------------------------------------------------
  subroutine update_stratO3(BTT,BTTBCK,AIRB,DTADV,MP)
    !// --------------------------------------------------------------------
    !// Set stratospheric O3 when running with only tropospheric chemistry,
    !// and DO3 for J-values also when stratospheric chemistry is included.
    !//
    !// This routine works in parallel region. To affect photolytical
    !// calculations, we also need to adjust DO3 (in cmn_jv.f) accordingly
    !// (it is set each NMET in p-setc_oc.f, routine set_atm, which operates
    !// outside parallel region).
    !//
    !// L40 uses old method for now.
    !//
    !// Amund Sovde, May - June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPARW, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY
    use cmn_chem, only: TMASS
    use cmn_fjx, only: DO3
    use cmn_parameters, only: LDEBUG, M_AIR, AVOGNR
    use cmn_oslo, only: LMTROP, trsp_idx, DINM, LVS2ADD2TC
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(out) :: BTTBCK
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK), DTADV
    integer, intent(in) :: MP

    !// Locals
    integer :: I,J,L,N,II,JJ, LMIN
    real(r8) :: R1
    !// --------------------------------------------------------------------

    !// Should not be here if tropospheric chemistry is not included.
    if (.not.LOSLOCTROP) return


    !// If tropospheric chemistry, but no stratospheric chemistry, then go on
    !// to update stratospheric O3 from climatology.
    !// After this if-block, DO3 for J-values will be updated.
    if (.not.LOSLOCSTRAT) then

       if (LPARW .eq. 40) then
          !// For now, use old routine for L40:
          call update_strato3ctm2(BTT,BTTBCK,AIRB,DTADV,MP)
          return
       end if

       !// Conversion factors for vmr to mmr, for O3
       R1  = TMASS(trsp_idx(1)) / M_AIR

    
       !// Loop over latitude (J is global, JJ is block)
       do J = MPBLKJB(MP), MPBLKJE(MP)
         JJ    = J - MPBLKJB(MP) + 1

         !// Update O3 from slightly above tropopause
         LMIN = maxval(LMTROP(:,J)) + 1 + LVS2ADD2TC

         !// Loop over longitude (I is global, II is block)
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II    = I - MPBLKIB(MP) + 1

           do L = LMIN, LPAR

             !// Set O3 from climatology, linearly interpolated between two
             !// monthly averages. From vmr to mass;
             BTT(L,trsp_idx(1),II,JJ) = STRATO3(L,I,J) * AIRB(L,II,JJ) * R1

             !// Update diagnostic array
             BTTBCK(L,trsp_idx(1),II,JJ) = BTT(L,trsp_idx(1),II,JJ)


           end do
 
         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    end if !// if (.not.LOSLOCSTRAT) then


    !// Update DO3. This is also done when stratospheric chemistry is included!
    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        do L = 1, LPAR
          !// Update DO3 
          DO3(L,I,J) = 1000._r8 * AVOGNR * BTT(L,trsp_idx(1),II,JJ) &
                       / 48._r8 * 1.e-4_r8 / AREAXY(I,J)
        end do !// do L = 1, LPAR
      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine update_stratO3
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine stratNOY_interp(MONTH)
    !// --------------------------------------------------------------------
    !// Interpolate climatology linearly in time between two monthly means.
    !// Done one time each day, outside parallel region.
    !//
    !// Amund Sovde, June 2010
    !// --------------------------------------------------------------------
    use cmn_oslo, only: DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: MONTH

    !// Locals
    real(r8) :: dtfrac
    integer :: I, J, L, NY
    !// --------------------------------------------------------------------


    !// Interpolate linearly in time between two monthly means,
    !// updated every day. Fraction of first (this) month:
    dtfrac = 1._r8 - real(jdate-1, r8) / DINM(MONTH)

    !// Loop over latitude (J is global, JJ is block)
    do J = 1, JPAR
      !// Loop over longitude (I is global, II is block)
      do I = 1, IPAR
        do L = 1, LPAR
          !// Set NOy from climatology, linearly interpolated between two
          !// monthly averages.
          do NY = 1,NXNY
            STRATNOY(NY,L,I,J) = NOYCLIM(1,NY,L,I,J) * dtfrac + &
                                 NOYCLIM(2,NY,L,I,J) * (1._r8 - dtfrac)
          end do
        end do
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine stratNOY_interp
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_stratNOX(BTT,BTTBCK,AIRB,DTADV,MP)
    !// --------------------------------------------------------------------
    !// Set stratospheric HNO3 and NOX when running with only tropospheric
    !// chemistry.
    !//
    !// Should NOT be used when NOX (trsp_idx(2)) is not used. Setting NO=NOX as
    !// in this routine was ok when NOX = NO in stratosphere (STRAT_SETNOX
    !// in oc_stratchem.f).
    !//
    !// Better alternative: See oc_update_stratNOX2
    !//
    !// This routine works in parallel region.
    !//
    !// Amund Sovde, May - June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPARW, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY
    use cmn_chem, only: TMASS
    use cmn_parameters, only: LDEBUG, M_AIR
    use cmn_oslo, only: LMTROP, trsp_idx, DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(out) :: BTTBCK
    real(r8), intent(in) :: AIRB(LPAR,IDBLK,JDBLK), DTADV
    integer, intent(in) :: MP

    !// Locals
    integer :: I,J,L,N,II,JJ, LMIN
    real(r8) :: R1, R4, R43
    !// --------------------------------------------------------------------

    !// Should not be here if tropospheric chemistry is not included.
    if (.not.LOSLOCTROP) return


    !// If tropospheric chemistry, but no stratospheric chemistry, then go on
    !// to update stratospheric O3 from climatology.
    !// After this if-block, DO3 for J-values will be updated.
    if (.not.LOSLOCSTRAT) then



       !// As in old CTM2 code, we set stratospheric NO and HNO3.
       !// This is done from a correlation of NOy with O3 in the lower
       !// stratosphere (i.e. volume mixing ratio).
       !// Old comment states that:
       !//   Try to set only NOX+HNO3 = NOY, and let chemistry partition.
       !//   NOY is very good corr. with O3 in lower stratosphere.
       !// Since STT is in mass, we need to convert O3 to vmr and then NO and
       !// HNO3 back to mass.
       R4  = TMASS(trsp_idx( 4)) / TMASS(trsp_idx(1))
       R43 = TMASS(trsp_idx(43)) / TMASS(trsp_idx(1))

    
       !// Loop over latitude (J is global, JJ is block)
       do J = MPBLKJB(MP), MPBLKJE(MP)
         JJ    = J - MPBLKJB(MP) + 1

         !// Update O3 from slightly above tropopause
         LMIN = maxval(LMTROP(:,J))

         !// Loop over longitude (I is global, II is block)
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II    = I - MPBLKIB(MP) + 1

 
           !// Update HNO3 and NO
           do L = LMTROP(I,J)+1, LPAR
             !// BTT is mass
             !// Convert HNO3 and NO mass from O3 mass, via vmr equivalents
             BTT(L,trsp_idx( 4),II,JJ) = 0.6_r8 * 4.e-3_r8 &
                                    * BTT(L,trsp_idx(1),II,JJ) * R4
             BTT(L,trsp_idx(43),II,JJ) = 0.4_r8 * 4.e-3_r8 &
                                    * BTT(L,trsp_idx(1),II,JJ) * R43
             !// Update diagnostic arrays
             BTTBCK(L,trsp_idx( 4),II,JJ) = BTT(L,trsp_idx( 4),II,JJ)
             BTTBCK(L,trsp_idx(43),II,JJ) = BTT(L,trsp_idx(43),II,JJ)
           end do

         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    end if !// if (.not.LOSLOCSTRAT) then


    !// --------------------------------------------------------------------
  end subroutine update_stratNOX
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_stratNOX2(BTT,BTTBCK,AIRB,MP)
    !// --------------------------------------------------------------------
    !// Set stratospheric HNO3 and NOX when running with only tropospheric
    !// chemistry, and when NOX (trsp_idx(2)) is removed.
    !//
    !// Alternative to this routine is to let NOX species be inert in
    !// stratosphere.
    !//
    !// This routine works in parallel region.
    !//
    !// Amund Sovde, May - June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPARW, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY
    use cmn_chem, only: TMASS
    use cmn_parameters, only: LDEBUG, M_AIR
    use cmn_oslo, only: LMTROP, trsp_idx, DINM
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(out) :: BTTBCK
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK)
    integer, intent(in) :: MP

    !// Locals
    integer :: I,J,L,N,II,JJ, LMIN
    real(r8) :: R1, R4, MNOX, CNOX, FRAC, MNOX2, XPAN
    !// --------------------------------------------------------------------

    !// Should not be here if tropospheric chemistry is not included.
    if (.not.LOSLOCTROP) return


    !// If tropospheric chemistry, but no stratospheric chemistry, then go on
    !// to update stratospheric O3 from climatology.
    !// After this if-block, DO3 for J-values will be updated.
    if (.not.LOSLOCSTRAT) then



       !// In old CTM2 code, we set stratospheric NO and HNO3. However, the
       !// NO value was set to some NOX value based on a correlation of NOy
       !// with O3 in the lower stratosphere (i.e. volume mixing ratio).
       !// Old comment states that:
       !//   Try to set only NOX+HNO3 = NOY, and let chemistry partition.
       !//   NOY is very good corr. with O3 in lower stratosphere.
       !// Setting NO=NOX worked, because in the stratosphere CTM2 also set NOX=NO
       !// when using only tropospheric chemistry.
       !// This is not good, and instead we should either let NOy be inert in the
       !// stratosphere or scale the existing NOX & HNO3 to the values from O3.
       !// Which is carried out below:

       !// Since STT is in mass, we need to convert O3 to vmr and then NO and
       !// HNO3 back to mass.
       R4  = TMASS(trsp_idx( 4)) / TMASS(trsp_idx(1))

    
       !// Loop over latitude (J is global, JJ is block)
       do J = MPBLKJB(MP), MPBLKJE(MP)
         JJ    = J - MPBLKJB(MP) + 1

         !// Update O3 from slightly above tropopause
         LMIN = maxval(LMTROP(:,J))

         !// Loop over longitude (I is global, II is block)
         do I = MPBLKIB(MP), MPBLKIE(MP)
           II    = I - MPBLKIB(MP) + 1

 
           !// Update HNO3 and NO
           do L = LMTROP(I,J)+1,LPAR
             !// BTT is mass
             !// Convert HNO3 mass from O3 mass, via vmr equivalents
             BTT(L,trsp_idx( 4),II,JJ) = 0.6_r8 * 4.e-3_r8 &
                                    * BTT(L,trsp_idx(1),II,JJ) * R4


             !// Molecular equivalents of strat NOX:
             CNOX = 0.4_r8 * 4.e-3_r8 * BTT(L,trsp_idx(1),II,JJ) &
                                        / TMASS(trsp_idx(1))
             !// Modeled molecular equivalents of NOX
             MNOX = BTT(L,trsp_idx(17),II,JJ) / TMASS(trsp_idx(17)) + &
                    BTT(L,trsp_idx(41),II,JJ) / TMASS(trsp_idx(41)) + &
                    2._r8 * BTT(L,trsp_idx(42),II,JJ) / TMASS(trsp_idx(42)) + &
                    BTT(L,trsp_idx(43),II,JJ) / TMASS(trsp_idx(43)) + &
                    BTT(L,trsp_idx(44),II,JJ) / TMASS(trsp_idx(44)) + &
                    !// PANX = PAN + CH3X
                    BTT(L,trsp_idx( 5),II,JJ) / TMASS(trsp_idx( 5)) &
                       - BTT(L,trsp_idx(37),II,JJ) / TMASS(trsp_idx(37))
             if (MNOX .gt. 1.e-10_r8) then
                FRAC = CNOX/MNOX
             else
                FRAC = 0._r8
             end if
             !// Remember CH3X
             XPAN= BTT(L,trsp_idx(37),II,JJ)
             !// Fraction to change components
             BTT(L,trsp_idx(17),II,JJ) = BTT(L,trsp_idx(17),II,JJ) * FRAC
             BTT(L,trsp_idx(41),II,JJ) = BTT(L,trsp_idx(41),II,JJ) * FRAC
             BTT(L,trsp_idx(42),II,JJ) = BTT(L,trsp_idx(42),II,JJ) * FRAC
             BTT(L,trsp_idx(43),II,JJ) = BTT(L,trsp_idx(43),II,JJ) * FRAC
             BTT(L,trsp_idx(44),II,JJ) = BTT(L,trsp_idx(44),II,JJ) * FRAC
             BTT(L,trsp_idx( 5),II,JJ) = (BTT(L,trsp_idx( 5),II,JJ) - XPAN) &
                                         * FRAC + XPAN

             MNOX2= BTT(L,trsp_idx(17),II,JJ) / TMASS(trsp_idx(17)) + &
                    BTT(L,trsp_idx(41),II,JJ) / TMASS(trsp_idx(41)) + &
                    2._r8 * BTT(L,trsp_idx(42),II,JJ) / TMASS(trsp_idx(42)) + &
                    BTT(L,trsp_idx(43),II,JJ) / TMASS(trsp_idx(43)) + &
                    BTT(L,trsp_idx(44),II,JJ) / TMASS(trsp_idx(44)) + &
                    BTT(L,trsp_idx( 5),II,JJ) / TMASS(trsp_idx( 5)) &
                       - BTT(L,trsp_idx(37),II,JJ) / TMASS(trsp_idx(37))
             FRAC = ((MNOX2-CNOX)/CNOX)*100._r8
             if (CNOX .gt. 1.e-10_r8 .and. abs(FRAC).gt.1.e-1_r8) then
                print*,'NOX not scaled correctly',FRAC
                print*,CNOX
                print*,MNOX
                print*,MNOX2
                stop
             end if

             !// Update diagnostic arrays
             BTTBCK(L,trsp_idx( 4),II,JJ) = BTT(L,trsp_idx( 4),II,JJ)
             BTTBCK(L,trsp_idx( 5),II,JJ) = BTT(L,trsp_idx( 5),II,JJ)
             BTTBCK(L,trsp_idx(17),II,JJ) = BTT(L,trsp_idx(17),II,JJ)
             BTTBCK(L,trsp_idx(41),II,JJ) = BTT(L,trsp_idx(41),II,JJ)
             BTTBCK(L,trsp_idx(42),II,JJ) = BTT(L,trsp_idx(42),II,JJ)
             BTTBCK(L,trsp_idx(43),II,JJ) = BTT(L,trsp_idx(43),II,JJ)
             BTTBCK(L,trsp_idx(44),II,JJ) = BTT(L,trsp_idx(44),II,JJ)
           end do

         end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
       end do !// do J = MPBLKJB(MP), MPBLKJE(MP)
    end if !// if (.not.LOSLOCSTRAT) then


    !// --------------------------------------------------------------------
  end subroutine update_stratNOX2
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_stratNOY(BTT,BTTBCK,AIRB,MP)
    !// --------------------------------------------------------------------
    !// Updates tropospheric NOy in the stratosphere, when using only
    !// tropopsheric chemistry, from climatology of CTM2 full chemistry.
    !//
    !//
    !// This routine works in parallel region.
    !//
    !// Amund Sovde, June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPARW, NPAR, IDBLK, JDBLK
    use cmn_ctm, only: MPBLKJB, MPBLKJE, MPBLKIB, MPBLKIE, AREAXY
    use cmn_chem, only: TMASS
    use cmn_parameters, only: LDEBUG, M_AIR
    use cmn_oslo, only: LMTROP, trsp_idx, DINM, LVS2ADD2TC
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// In/Out
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(out) :: BTTBCK
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK)
    integer, intent(in) :: MP

    !// Locals
    integer :: I,J,L,N,II,JJ, LMIN
    real(r8) :: R4, R5, R17, R41, R42,R43,R44, XPAN
    !// --------------------------------------------------------------------


    !// Should not be here if tropospheric chemistry is not included.
    if (.not.LOSLOCTROP) return

    !// If tropospheric chemistry, but no stratospheric chemistry, then
    !// set NOy.
    if (LOSLOCSTRAT) return

    !// For L40, use the old CTM2 treatment for now, until we can
    !// interpolate from L60-climatology.
    if (LPARW.eq.40) then
       call update_stratNOX2(BTT,BTTBCK,AIRB,MP)
       return
    end if

    !// Conversion factor for vmr to mmr, for all NOy species
    R4  = TMASS(trsp_idx( 4)) / M_AIR
    R5  = TMASS(trsp_idx( 5)) / M_AIR
    R17 = TMASS(trsp_idx(17)) / M_AIR
    R41 = TMASS(trsp_idx(41)) / M_AIR
    R42 = TMASS(trsp_idx(42)) / M_AIR
    R43 = TMASS(trsp_idx(43)) / M_AIR
    R44 = TMASS(trsp_idx(44)) / M_AIR

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
      JJ    = J - MPBLKJB(MP) + 1

      !// Update O3 from slightly above tropopause. First level
      !// to set climatology will be where the lowest stratosphere
      !// is highest in this latitude band:
      LMIN = maxval(LMTROP(:,J)) + 1 + LVS2ADD2TC

      !// Loop over longitude (I is global, II is block)
      do I = MPBLKIB(MP), MPBLKIE(MP)
        II    = I - MPBLKIB(MP) + 1

        !// Update NOy in stratosphere
        do L = LMIN, LPAR     !// Could use LMTROP(I,J)+1,LPAR
          !// BTT is mass
          !// HNO3
          BTT(L,trsp_idx( 4),II,JJ) = STRATNOY(1,L,I,J) * AIRB(L,II,JJ) * R4
          !// PANX = PAN + CH3X
          !// Find PAN mass equivalent of CH3X
          XPAN = BTT(L,trsp_idx(37),II,JJ) &
                 / TMASS(trsp_idx(37)) * TMASS(trsp_idx(5))
          BTT(L,trsp_idx( 5),II,JJ) = &
               STRATNOY(2,L,I,J) * AIRB(L,II,JJ) * R5 + XPAN


          !// HO2NO2
          BTT(L,trsp_idx(17),II,JJ) = STRATNOY(3,L,I,J) * AIRB(L,II,JJ) * R17
          !// NO3
          BTT(L,trsp_idx(41),II,JJ) = STRATNOY(4,L,I,J) * AIRB(L,II,JJ) * R41
          !// N2O5
          BTT(L,trsp_idx(42),II,JJ) = STRATNOY(5,L,I,J) * AIRB(L,II,JJ) * R42
          !// NO
          BTT(L,trsp_idx(43),II,JJ) = STRATNOY(6,L,I,J) * AIRB(L,II,JJ) * R43
          !// NO2
          BTT(L,trsp_idx(44),II,JJ) = STRATNOY(7,L,I,J) * AIRB(L,II,JJ) * R44


          !// Update diagnostic arrays
          BTTBCK(L,trsp_idx( 4),II,JJ) = BTT(L,trsp_idx( 4),II,JJ)
          BTTBCK(L,trsp_idx( 5),II,JJ) = BTT(L,trsp_idx( 5),II,JJ)
          BTTBCK(L,trsp_idx(17),II,JJ) = BTT(L,trsp_idx(17),II,JJ)
          BTTBCK(L,trsp_idx(41),II,JJ) = BTT(L,trsp_idx(41),II,JJ)
          BTTBCK(L,trsp_idx(42),II,JJ) = BTT(L,trsp_idx(42),II,JJ)
          BTTBCK(L,trsp_idx(43),II,JJ) = BTT(L,trsp_idx(43),II,JJ)
          BTTBCK(L,trsp_idx(44),II,JJ) = BTT(L,trsp_idx(44),II,JJ)
        end do !// do L = LMTROP(I,J)+1,LPAR
      end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)


    !// --------------------------------------------------------------------
  end subroutine update_stratNOY
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_strato3ctm2(BTT,BTTBCK,AIRB,DT,MP)
    !// --------------------------------------------------------------------
    !// Update O3 in stratosphere when stratosphere is not
    !// calculated.
    !//
    !// Old routine for L40! Should be abandoned!
    !//
    !// Amund Sovde, November 2009 - June 2010
    !// --------------------------------------------------------------------
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LOSLOCSTRAT, LOSLOCTROP
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_oslo, only: LMTROP, trsp_idx, O3_UPBND
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: BTT
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(out) :: BTTBCK
    real(r8), intent(in)  :: AIRB(LPAR,IDBLK,JDBLK), DT
    integer, intent(in) :: MP

    !// Locals
    Integer :: I,J,II,JJ
    real(r8) :: RFRAC
    !// --------------------------------------------------------------------

    !// Should not be here if stratospheric chemistry is included
    !// and not if tropospheric chemistry is not included.
    if (LOSLOCSTRAT .or. (.not.LOSLOCTROP)) return

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       RFRAC = O3_UPBND(J) * DT

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          BTT(LPAR,trsp_idx(1),II,JJ) = BTT(LPAR,trsp_idx(1),II,JJ) + RFRAC
          !// Update diagnostic array
          BTTBCK(LPAR,trsp_idx(1),II,JJ) = BTT(LPAR,trsp_idx(1),II,JJ)

       end do
    end do

    !// --------------------------------------------------------------------
  end subroutine update_strato3ctm2
  !// ----------------------------------------------------------------------

  !// ----------------------------------------------------------------------
end module strat_o3noy_clim
!//=========================================================================

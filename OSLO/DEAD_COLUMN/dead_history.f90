!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, March 2016
!//=========================================================================
!// DEAD history fields.
!//=========================================================================
module dead_history
  !//=======================================================================
  !// Purpose: 
  !// Keep track of history fields, such as budget terms.
  !// These are instant values, to be available after all dust calculations
  !// have been carried out.
  !//
  !// Except dustprod and dustdrydep, none of the fields were used by
  !// the Oslo CTM3. I have therefore commented out most of the fields.
  !// This is also done in subroutine outfld_1 and in dstmbl.F90.
  !// Ole Amund Sovde, March 2016
  !// ----------------------------------------------------------------------
  !// Physical grid
  use pmgrid, only: plon, plat
  !// Precision
  use dead_precision, only: r8
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  !// Could in principle have these local
  real(r8), dimension(plon,plat) :: &
       dustdrydep, &         !Dry dep in kg/m2/s
       dustprod !, &           !PROD in kg/m2/s
       !dns_mdp, &            !rho (kg/m3)
       !flx_SW_abs, &         !Shortwave abs
       !flx_LW_dwn, &         !Longwave down
       !hgt_mdp, &            !Midpoint height
       !hgt_zpd, &            !Zero plane displ
       !lnd_frc_dry, &        !Dry land fraction
       !mbl_bsn_fct, &        !Mobilizaion basin factor
       !mno_lng, &            !Monin obukhov length
       !oro, &                !oreography
       !prs_dlt, &            !Pressure thickness
       !prs_mdp, &            !Midpoint pressure (Pa)
       !prs_sfc, &            !Surface pressure (Pa)
       !q_h2o_vpr, &          !Q
       !rgh_mmn_dep, &        !z0 in dry deposition
       !snw_frc, &            !snow fraction
       !tpt_gnd, &            !Ground temp
       !tpt_mdp, &            !T
       !tpt_ptn, &            !Potential temp
       !wnd_frc, &            !friction wind speed
       !wnd_rfr, &            !wind @ 10m height
       !wnd_znl_mdp, &        !zonal wind speed
       !wnd_mrd_mdp           !meridional wind speed
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='dead_history.f90'
  !// ----------------------------------------------------------------------
  save
  private
  public outfld_1, gethistfld_ij, gethistfld_2d
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine outfld_1(string,field,lon_idx,lat_idx,obuf)
    !// --------------------------------------------------------------------
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    character(len=*), intent(in) :: string !Name of field we want to set
    integer, intent(in)  :: lon_idx        !longitudes is number of longitudes,
    integer, intent(in)  :: lat_idx        !latitude is current latitude
    real(r8), intent(in) :: field          !Field to set
    real(r8), intent(in) :: obuf(1)        !Variable which is not used

    !// Local
    integer :: lenstring  !Length of string
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='outfld_1.f90'
    !// --------------------------------------------------------------------
    !// Find length of string
    lenstring = len(string)

    if (string(1:lenstring).eq.'DSTSFMBL') then
       !// setting dust total production :
       dustprod(lon_idx,lat_idx) = field
    else if (string(1:lenstring).eq.'DSTSFDRY') then !Dry dep flux
       !// Setting dust dry deposition
       dustdrydep(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'DNS_MDP') then
!       dns_mdp(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'FLX_SWA') then
!       flx_SW_abs(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'FLX_LWD') then
!       flx_LW_dwn(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'HGT_MDP') then
!       hgt_mdp(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'HGT_ZPD') then
!       hgt_zpd(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'LND_MBL') then
!       !// Land fraction for mobilization
!       lnd_frc_dry(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'BSN_FCT') then
!       mbl_bsn_fct(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'MNO_LNG') then
!       mno_lng(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'ORO') then
!       oro(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'PRS_DLT') then
!       prs_dlt(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'PRS_MDP') then
!       prs_mdp(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'PRS_SFC') then
!       prs_sfc(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'Q_H2O') then
!       q_h2O_vpr(lon_idx,lat_idx) = field 
!    else if (string(1:lenstring).eq.'RGH_MMN_DEP') then
!       rgh_mmn_dep(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'SNW_FRC') then
!       snw_frc(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'TPT_MDP') then
!       tpt_mdp(lon_idx,lat_idx) = field   
!    else if (string(1:lenstring).eq.'TPT_GND') then
!       tpt_gnd(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'TPT_PTN') then
!       tpt_ptn(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'WND_FRC_MBL') then
!       !// Wind friction speed
!       wnd_frc(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'WND_RFR_MBL') then
!       !// Wind at reference height
!       wnd_rfr(lon_idx,lat_idx) = field 
!    else if (string(1:lenstring).eq.'WNDZNLMDP') then
!       !// Zonal wind at layer midpoint
!       wnd_znl_mdp(lon_idx,lat_idx) = field
!    else if (string(1:lenstring).eq.'WNDMRDMDP') then
!       !// Meridional wind at layer midpoint
!       wnd_mrd_mdp(lon_idx,lat_idx) = field
    else
       !// Error; not implemented
       write(6,'(a)') f90file//':'//subr// &
            ': Not implemented for field name '//trim(string)
       stop
    end if

    !// --------------------------------------------------------------------
  end subroutine outfld_1
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine gethistfld_ij(string,field,lon_idx,lat_idx)
    !// --------------------------------------------------------------------
    !// Retrieve a history field for a given grid box.
    !// Not fully set up, but included so it is easy to modify.
    !// See also gethistfield_2d to retrieve a 2D field.
    !//
    !// Ole Amund Sovde, March 2016
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), intent(in) :: string !Name of field to retrieve
    integer, intent(in) :: lon_idx         !longitudes is number of longitudes,
    integer, intent(in) :: lat_idx         !latitude is current latitude
    real(r8), intent(out) :: field         !Retrieved value

    !// Local
    integer :: lenstring !Length of string
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='gethistfld_ij.f90'
    !// --------------------------------------------------------------------

    !// Find length of string
    lenstring = len(string)

    if (string(1:lenstring).eq.'DSTSFMBL') then
       !// dust production
       field = dustprod(lon_idx,lat_idx)
    else if (string(1:lenstring).eq.'DSTSFDRY') then
       !// dust dry deposition
       field = dustdrydep(lon_idx,lat_idx)
    else
       !// Error; not implemented
       write(6,'(a)') f90file//':'//subr// &
            ': Not implemented for field name '//trim(string)
       stop
    end if
    !// ----------------------------------------------------------------------
  end subroutine gethistfld_ij
  !// ------------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine gethistfld_2d(string,field)
    !// --------------------------------------------------------------------
    !// Retrieve a history 2D field.
    !// Not fully set up, but included so it is easy to modify.
    !// See also gethistfield_ij to retrieve a single grid box field.
    !//
    !// Ole Amund Sovde, March 2016
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    character(len=*), intent(in) :: string    !Name of field to retrieve
    real(r8), intent(out) :: field(plon,plat) !Retrieved value

    !// Local
    integer :: lenstring !Length of string
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr='gethistfld_2d.f90'
    !// --------------------------------------------------------------------

    !// Find length of string
    lenstring = len(string)

    if (string(1:lenstring).eq.'DSTSFMBL') then
       !// dust production
       field(:,:) = dustprod(:,:)
    else if (string(1:lenstring).eq.'DSTSFDRY') then
       !// dust dry deposition
       field(:,:) = dustdrydep(:,:)
    else
       !// Error; not implemented
       write(6,'(a)') f90file//':'//subr// &
            ': Not implemented for field name '//trim(string)
       stop
    end if
    !// ----------------------------------------------------------------------
  end subroutine gethistfld_2d
  !// ------------------------------------------------------------------------


  !//=======================================================================
end module dead_history
!//=========================================================================

      module histout
C
C     usage: use histout, only...(selected fields or subroutines)
C
C     Purpose: 
C     Set field from dust module
C     Main program should be able to use this module pick the fields it wants for averaging
C     NOTE: There should be a way to know how many times fields have been set and average  !
C
      use pmgrid        !Physical grid
      use precision     !real4/real8 definition

      implicit none

C     PRIVATE !All variables/subroutines are privat unless decleared otherwise
      PUBLIC                    !:: outfld, dustfield2main
C     These fields are local in this module
      real(r8)        :: dustdrydep(plon,plat)    !Dry dep in kg/m2/s
      real(r8)        :: dustprod(plon,plat)      !PROD in kg/m2/s
      real(r8)        :: lnd_frc_dry(plon,plat)   !Dry land fraction
      real(r8)        :: wnd_rfr(plon,plat)       !wind @ 10m height
      real(r8)        :: mbl_bsn_fct(plon,plat)   !Mobilizaion basin factor
      real(r8)        :: wnd_znl_mdp(plon,plat)   !zonal wind speed
      real(r8)        :: wnd_mrd_mdp(plon,plat)   !meridional wind speed
      real(r8)        :: wnd_frc(plon,plat)       !friction wind speed
      real(r8)        :: hgt_mdp(plon,plat)       !Midpoint height
      real(r8)        :: hgt_zpd(plon,plat)       !Zero plane displ
      real(r8)        :: mno_lng(plon,plat)       !Monin obukhov length
      real(r8)        :: prs_dlt(plon,plat)       !Pressure thickness
      real(r8)        :: oro(plon,plat)           !oreography
      real(r8)        :: q_h2o_vpr(plon,plat)     !Q
      real(r8)        :: tpt_mdp(plon,plat)       !T
      real(r8)        :: dns_mdp(plon,plat)       !rho (kg/m3)
      real(r8)        :: flx_SW_abs(plon,plat)    !Shortwave abs
      real(r8)        :: flx_LW_dwn(plon,plat)    !Longwave down
      real(r8)        :: snw_frc(plon,plat)       !snow fraction
      real(r8)        :: prs_mdp(plon,plat)       !Midpoint pressure (Pa)
      real(r8)        :: prs_sfc(plon,plat)       !Surface pressure (Pa)
      real(r8)        :: tpt_gnd(plon,plat)       !Ground temp
      real(r8)        :: tpt_ptn(plon,plat)       !Potential temp
      real(r8)        :: mbl_frc_log(plon,plat)   !fraction with logical true
      real(r8)        :: vdep1(plon,plat)         !Dry deposition 1
      real(r8)        :: vdep2(plon,plat)         !Dry deposition 2
      real(r8)        :: vdep3(plon,plat)         !Dry deposition 3
      real(r8)        :: vdep4(plon,plat)         !Dry deposition 4
      real(r8)        :: rgh_mmn_dep(plon,plat)   !z0 in dry deposition
C
      contains

c*******************************************************************
      subroutine outfld(string,field,longitudes,latitude,obuf)
c
      implicit none
c
      character,intent(in)    ::string*(*)              !Name of field we want to set
      integer,intent(in)      :: longitudes             !longitudes is number of longitudes,
      integer, intent(in)     :: latitude               !latitude is current latitude
      real(r8),intent(in)     ::field(longitudes)       !Field is the fiels which we want to set
      real(r8),intent(in)     ::obuf(1)                 !Variable which is not used
C
C     LOCAL
      integer                 :: lenstring              !Length of string
C
C      write(6,*)'inside outfld'
C      write(6,*)string,longitudes,latitude,obuf
C      write(6,*)'lbound',LBOUND(field)
C      write(6,*)'ubound',UBOUND(field)
C      write(6,*)'string(1:7)',string
C      write(6,*)'what is plon ?',plon
C      write(6,*)'what is plat ?',plat
C      write(6,*)'Length of string 1 ',len_trim(string)
C      write(6,*)'length of string 2',len(string)

C     Find length of string
      lenstring=len(string)

c     Begin code
C     setting dust production :
      if(string(1:lenstring).eq.'DSTSFMBL')then
         dustprod(:,latitude)=field(:)
         !write(6,*)'setting production at latitude',latitude
C         write(6,*)'maxval ',maxval(dustprod)
C     Setting dust dry deposition
      elseif(string(1:lenstring).eq.'DSTSFDRY')then !Dry dep flux
         dustdrydep(:,latitude)=field(:)
         !write(6,*)'setting dustdrydep at latitude',latitude
      elseif(string(1:lenstring).eq.'LND_MBL')then !Land fraction for mobilization
         !write(6,*)'setting lndfrcdry at latitude',latitude
         lnd_frc_dry(:,latitude)=field(:)
         !write(6,*)'setting lnd_mbl at latitude',latitude
      elseif(string(1:lenstring).eq.'WND_RFR_MBL')then      !Wind at reference height
         !write(6,*)'setting wnd_frc for lat = ',latitude
         wnd_rfr(:,latitude)=field(:) 
      elseif(string(1:lenstring).eq.'WND_FRC_MBL')then      !Wind friction speed
         wnd_frc(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'WNDZNLMDP')then !Zonal wind at layer midpoint
         wnd_znl_mdp(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'WNDMRDMDP')then !Meridional wind at layer midpoint
         wnd_mrd_mdp(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'PRS_DLT')then
         prs_dlt(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'HGT_MDP')then
         hgt_mdp(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'ORO')then
         oro(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'Q_H2O')then
         q_h2O_vpr(:,latitude)=field(:) 
      elseif(string(1:lenstring).eq.'TPT_MDP')then
         tpt_mdp(:,latitude)=field(:)   
      elseif(string(1:lenstring).eq.'DNS_MDP')then
         dns_mdp(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'MNO_LNG')then
         mno_lng(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'HGT_ZPD')then
         hgt_zpd(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'FLX_SWA')then
         flx_SW_abs(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'FLX_LWD')then
         flx_LW_dwn(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'SNW_FRC')then
         snw_frc(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'PRS_MDP')then
         prs_mdp(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'PRS_SFC')then
         prs_sfc(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'TPT_GND')then
         tpt_gnd(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'TPT_PTN')then
         tpt_ptn(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'MBLFRCL')then
         mbl_frc_log(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'VDEP1')then
         vdep1(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'VDEP2')then
         vdep2(:,latitude)=field(:)
      elseif(string(1:lenstring).eq.'VDEP3')then
         vdep3(:,latitude)=field(:)
C         write(6,*)'setting vdep at latitude',latitude
      elseif(string(1:lenstring).eq.'VDEP4')then
         vdep4(:,latitude)=field(:)
C         write(6,*)'setting vdep4 at latitude',latitude
      elseif(string(1:lenstring).eq.'RGH_MMN_DEP')then
         rgh_mmn_dep(:,latitude)=field(:)
C         write(6,*)'setting rgh_mmn at latitude',latitude
C
      else
C         write(6,*)'doing nothing'
c     DO NOTHING
      endif

C      write(6,*)'DONE outfld'
c      
      end subroutine outfld
c****************************************************************
      subroutine dustfield2main(string, field, longitudes,J)
C
C     Purpose: should be called from main program
C     A random field, "field" in main program will get the value of one of
C     the local fields in this routine. The fields to access are fields in the DEAD
C     dust model. In dead, a call to outfld, will set one of the fields in this module

      implicit none
C
C     INPUTS
      character,intent(in)    :: string*(*)  !Name of field 
      integer, intent(in)     :: longitudes  !Number of longitudes
      integer, intent(in)     :: J           !Latitude in question

C     OUTPUTS
      real(r8), intent(out)   :: field(longitudes)  !The field which will go to main
C      
C     LOCAL
      integer                :: lenstring  !Length of string

C     Check that longitudes is OK
      if(plon.ne.longitudes)then
         write(6,*)'not the same longitudes in CTM and dust model'
         stop
      endif

C     Get length of string with variable name
      lenstring=len_trim(string)
C
C     Start getting the Output field
      if(string(1:lenstring).eq.'DSTSFMBL')then
         field(:)=dustprod(:,J)
      elseif(string(1:lenstring).eq.'DSTSFDRY')then
         field(:)=dustdrydep(:,J)
      elseif(string(1:lenstring).eq.'LND_MBL')then
         field(:)=lnd_frc_dry(:,J)
      elseif(string(1:lenstring).eq.'WND_RFR_MBL')then
         field(:)=wnd_rfr(:,J)
      elseif(string(1:lenstring).eq.'WND_FRC_MBL')then      !Wind friction speed
         !write(6,*)'getting wnd_frc for lat/lon =',longitudes,latitudes
         field(:)=wnd_frc(:,J)
      elseif(string(1:lenstring).eq.'WNDZNLMDP')then     !Zonal wind at layer midpoint
         field(:)=wnd_znl_mdp(:,J)
      elseif(string(1:lenstring).eq.'WNDMRDMDP')then     !Meridional wind at layer midpoint
         field(:)=wnd_mrd_mdp(:,J)   
      elseif(string(1:lenstring).eq.'PRS_DLT')then
         field(:)=prs_dlt(:,J)
      elseif(string(1:lenstring).eq.'HGT_MDP')then
         field(:)=hgt_mdp(:,J)
      elseif(string(1:lenstring).eq.'ORO')then
         field(:)=oro(:,J)
      elseif(string(1:lenstring).eq.'Q_H2O')then
         field(:)=q_h2O_vpr(:,J) 
      elseif(string(1:lenstring).eq.'TPT_MDP')then
         field(:)=tpt_mdp(:,J)
      elseif(string(1:lenstring).eq.'DNS_MDP')then
         field(:)=dns_mdp(:,J)
      elseif(string(1:lenstring).eq.'MNO_LNG')then
         field(:)=mno_lng(:,J)
      elseif(string(1:lenstring).eq.'HGT_ZPD')then
         field(:)=hgt_zpd(:,J)
      elseif(string(1:lenstring).eq.'FLX_SWA')then
         field(:)=flx_SW_abs(:,J)
      elseif(string(1:lenstring).eq.'FLX_LWD')then
         field(:)=flx_LW_dwn(:,J)
      elseif(string(1:lenstring).eq.'SNW_FRC')then
         field(:)=snw_frc(:,J)
      elseif(string(1:lenstring).eq.'PRS_MDP')then
         field(:)=prs_mdp(:,J)
      elseif(string(1:lenstring).eq.'PRS_SFC')then
         field(:)=prs_sfc(:,J)
      elseif(string(1:lenstring).eq.'TPT_GND')then
         field(:)=tpt_gnd(:,J)
      elseif(string(1:lenstring).eq.'TPT_PTN')then
         field(:)=tpt_ptn(:,J)
      elseif(string(1:lenstring).eq.'MBLFRCL')then
         field(:)=mbl_frc_log(:,J)
      elseif(string(1:lenstring).eq.'VDEP1')then
         field(:)=vdep1(:,J)   
      elseif(string(1:lenstring).eq.'VDEP2')then
         field(:)=vdep2(:,J)
      elseif(string(1:lenstring).eq.'VDEP3')then
         field(:)=vdep3(:,J)   
      elseif(string(1:lenstring).eq.'VDEP4')then
         field(:)=vdep4(:,J)
      elseif(string(1:lenstring).eq.'RGH_MMN_DEP')then
         field(:)=rgh_mmn_dep(:,J)

      else
         write(6,*)'ERROR IN dustfield2main: '
         write(6,*)'You are trying to access a field which is not set'
         write(6,*)'This field name is wrong ',string
         STOP
      endif
C
      end subroutine dustfield2main
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end module histout

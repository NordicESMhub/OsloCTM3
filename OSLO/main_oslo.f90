!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Main driver for Oslo chemistry
!//=========================================================================
module main_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: main_oslo
  !// DESCRIPTION: Routines for controlling Oslo chemistry/physics.
  !//
  !// Contains:
  !//   subroutine master_oslo
  !//   subroutine jvalues_oslo
  !//   subroutine jv_column
  !//   subroutine update_chemistry_oslo
  !//
  !// Ole Amund Sovde, September/October 2008
  !//                  Updated March 2009
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file='main_oslo.f90'
  !// ----------------------------------------------------------------------
  public
  private jvalues_oslo, jv_column
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine master_oslo(BTT, BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ, &
       BTEM, AIRB, BTTBCK, UTTAU, ORG_DTCHM, DTCHM2, &
       NDAY, NMET, NOPS, NSUB, NCLDRAN, CCYC, MP)
    !// --------------------------------------------------------------------
    !// Master routine to control Oslo chemistry.
    !//
    !// Oslo chemistry has an internal time step of DTCHM2, whereas the
    !// overall chemical time step ORG_DTCHM is set by NRCHEM. The internal
    !// loop is carried out in p-main and loops CHMCYCLES times, so we have
    !//   DTCHM2 = ORG_DTCHM/CHMCYCLES
    !// Shorter time steps may be taken in the following calls.
    !// Remember that
    !//   ORG_DTCHM = DTOPS/NRCHEM
    !//
    !// Processes for each IJ-block is called separately, but in principle
    !// the processes could be carried out for a column (inside a IJ-loop).
    !// That could possibly save some CPU time, but probably not very much.
    !//
    !// Ole Amund Sovde, September/October 2008
    !//                  Updated March 2009
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8, rMom
    use cmn_size, only: LPAR, NPAR, IDBLK, JDBLK, LEMISDEP_INCHEM, &
         LOSLOCTROP, LOSLOCSTRAT, LSALT, LBCOC, LDUST, LNITRATE
    use cmn_ctm, only: JMON
    use cmn_met, only: Q
    use cmn_sfc, only: ZOI
    use averages, only: AVG_ADD2_H2O
    use cloudjx, only:  LCLDRANA, LCLDRANQ
    use utilities, only: check_btt
    !// --------------------------------------------------------------------
    use cmn_oslo, only: JVAL_IJ, AIRMOLEC_IJ, DV_IJ, EMIS_IJ, LJCCYC
    use aerosols2fastjx, only: set_aer4fjx, set_aer4fjx_ctm2, LJV_AEROSOL
    use bcoc_oslo, only: bcoc_master
    use diagnostics_general, only: sumup_burden_and_lifetimes, ch4n2o_burden
    use dust_oslo, only: dust_master
    use emisdep4chem_oslo, only: emis4chem
    use emisutils_oslo, only: emis_diag
    use psc_microphysics, only: oslochem_psc
    use nitrate, only: nitrate_master
    use physics_oslo, only: metdata_ij
    use seasalt, only: seasalt_master
    use stratchem_oslo, only: oslochem_strat, set_fam_in_trop
    use strat_loss, only: stratloss_oslo
    use strat_h2o, only: set_strat_h2o_b4chem
    use strat_o3noy_clim, only: update_stratO3, update_stratNOY
    use tropchem_oslo, only: oslochem_trop
    use utilities, only: adjust_moments
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAY, NMET, NOPS, NSUB, NCLDRAN, CCYC, MP
    real(r8), dimension(LPAR,IDBLK,JDBLK), intent(in) ::  AIRB, BTEM
    real(r8), intent(in)  ::  UTTAU, ORG_DTCHM, DTCHM2
    !// Inout/Output
    real(r8), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BTT, BTTBCK
    real(rMom), dimension(LPAR,NPAR,IDBLK,JDBLK), intent(inout) :: &
         BXT, BXX, BYT, BYY, BZT, BZZ, BXY, BXZ, BYZ
    !// --------------------------------------------------------------------
    !// Locals
    real(r8) :: dtjval !// Time span for which to calculate J-values
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'master_oslo'
    !// --------------------------------------------------------------------


    !// Store IJ-block arrays which changes only each meteorological time step
    !// --------------------------------------------------------------------
    if (NOPS.eq.1 .and. NSUB.eq.1 .and. CCYC.eq.1) then
       !// AIRMOLEC_IJ: air number density (molec/cm3)
       !// DV_IJ: Volume of grid boxes
       !// Note the new structure of AIRMOLEC_IJ and DV_IJ,
       !// with indexing (LPAR,NPAR,IDBLK,JDBLK,MPBLK). These are therefore
       !// not private B-arrays, but on B-array dimensions for all MPBLK.
       !//
       !// It can be argued that since AIR is changed by transport, we should
       !// calculate new pressure and then new DV and AIR_MOLEC_IJ after each
       !// transport. That may be done later if necessary (should be done if
       !// routine is called for each NOPS).
       call metdata_ij(AIRB, BTEM, MP)
    end if
    !// Update arrays that can change more frequently than each NOPS
    !// --------------------------------------------------------------------
    if (LJCCYC .or. (NSUB.eq.1.and.CCYC.eq.1)) then
       !// Update stratospheric O3 **before** J-values are calculated,
       !// and **before** NOY is updated (in case of L40). When L40 NOY
       !// is interpolated from climatology, this routine should be moved
       !// down to J-values section.
       call update_stratO3(BTT,BTTBCK,AIRB,DTCHM2,MP)
    end if

    !// Update arrays changing only each NOPS
    !// --------------------------------------------------------------------
    if (NSUB.eq.1 .and. CCYC.eq.1) then
       !// Generate Oslo style emissions? (Oslo style deposition is treated
       !// in chemistry as source terms). Get emissions in [kg/s], convert
       !// later.
       if (LEMISDEP_INCHEM) then
          !// Only send in this IJ-block (emisdep4chem.f90)
          call emis4chem(EMIS_IJ(:,:,:,:,MP), BTEM, UTTAU, MP)

          !// Accumulate the emissions for diagnose:
          !// NOTE that in this routine ORG_DTCHM is multiplied
          !// with NRCHEM to account for the whole NOPS step!
          call emis_diag(EMIS_IJ(:,:,:,:,MP), ORG_DTCHM, MP)
       end if
       !// Update stratospheric NOy when only tropospheric chemistry
       if (LOSLOCTROP .and. .not.LOSLOCSTRAT) &
            call update_stratNOY(BTT,BTTBCK,AIRB,MP)
    end if


    !// If stratospheric H2O is not transported, but is calculated from the
    !// sum of H2, it must be set before chemistry.
    call set_strat_h2o_b4chem(AIRB,BTT,MP)


    !// Calculate J-values.
    !// --------------------------------------------------------------------
    !// There are two options, controlled by LJCCYC in the input file.
    !// LJCCYC = .false.:
    !//   Calculated once every NRCHEM; only done once during each internal
    !//   loop, i.e. when CCYC=1 and NSUB=1.
    !// LJCCYC = .true.:
    !//   Calculated for ALL CCYC.
    if (LOSLOCTROP) then
       if (LJCCYC .or. (NSUB.eq.1.and.CCYC.eq.1)) then

          !// Calculate at UTTAU + DTJVAL/3600.
          if (.not. LJCCYC) then
             !// Same JVALs for all internal chemistry loops during one NRCHEM.
             !// Calculate JVALs in the middle of two NRCHEMs
             dtjval = ORG_DTCHM*0.5_r8
             dtjval = 0._r8
          else
             !// JVALs every chemical step (May try DTCHM2*0.5_r8)
             dtjval = 0._r8
          end if
          if (LJV_AEROSOL) then
             !call set4fjx_aer(BTT,MP)
             call set_aer4fjx_ctm2(MP) !// Simple profile as in CTM2
          end if
          call jvalues_oslo(JVAL_IJ,UTTAU,DTjval,NCLDRAN,BTEM,MP)
       end if
    end if

    !// Loss of tropospheric components in stratosphere
    !// All packages (e.g. sulphur) are located here.
    call stratloss_oslo(BTT,DTCHM2,MP)


    !// PSC microphysics for stratosphere
    call oslochem_psc(BTT,BTEM,DTCHM2,MP)

    !// Sum up H2O in separate aerage evey NOPS (as in the other
    !// averages)
    if (NSUB.eq.1.and.CCYC.eq.1) call AVG_ADD2_H2O(BTT,AIRB,MP)


    !// SALT - Do SEA SALT before chemistry (seasalt.f90)
    if (LSALT) call seasalt_master(AIRB, BTEM, DV_IJ(:,:,:,MP), BTT, &
         Q, ZOI, JMON, DTCHM2, MP)


    !// BC/OC (bcoc_oslo.f90)
    if (LBCOC) call bcoc_master(BTT, EMIS_IJ(:,:,:,:,MP), DTCHM2, MP)


    !// DUST (dust_oslo.f90)
    if (LDUST) call dust_master(BTT, AIRB, BTEM, DTCHM2, NOPS, MP)

 
    !// Find CH4 and N2O burden for lifetime diagnoses
    call ch4n2o_burden(BTT,NSUB.eq.1.and.CCYC.eq.1,DTCHM2,MP)

    !// Do tropospheric chemistry
    if (LOSLOCTROP) call oslochem_trop(BTT, JVAL_IJ(:,:,:,:,MP), &
         AIRMOLEC_IJ(:,:,:,MP), &
         DV_IJ(:,:,:,MP), BTEM, AIRB, EMIS_IJ(:,:,:,:,MP), DTCHM2, MP)

    !// NITRATE aerosols (chemistry in oslo_chem/pchemc_ij.f, called from
    !// oslochem_trop)
    if (LNITRATE) &
         call nitrate_master(BTT, BTEM, Q, DV_IJ(:,:,:,MP), AIRB, MP)


    !// Stratospheric chemistry
    if (LOSLOCSTRAT) call oslochem_strat(BTT, JVAL_IJ, &
         AIRMOLEC_IJ, DV_IJ, BTEM, &
         EMIS_IJ(:,:,:,:,MP), BTTBCK, DTCHM2, MP)


    !// Calculate stratospheric families in troposphere, so that transport
    !// has correct families
    if (LOSLOCSTRAT) call set_fam_in_trop(BTT,BTTBCK,MP)

    !// Calculate some burdens and lifetimes for the troposphere
    call sumup_burden_and_lifetimes(BTT,AIRB,DV_IJ,BTEM,MP)


    call check_btt(BTT,MP,'end of '//f90file//':'//subr)

    !// Reduce value of moments if tracer has been reduced
    call adjust_moments(BTT,BTTBCK,BXT,BXX,BYT,BYY,BZT,BZZ,BXY,BXZ,BYZ,MP)

    !// --------------------------------------------------------------------
  end subroutine master_oslo
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine jvalues_oslo(BJV,UTTAU,DT_ADD,NCLDRAN,BTEM,MP)
    !// --------------------------------------------------------------------
    !// Collects J-values in a B-array, since tropospheric and stratospheric
    !// chemistry is called separately.
    !//
    !// Ole Amund Sovde, September 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR, IDBLK, JDBLK, MPBLK
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_fjx, only: NQD_, JVN_
    use cloudjx, only: SWSTORE, CLDSTORE, TYPSTORE, LCLDRANA, LCLDRANQ
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NCLDRAN, MP
    real(r8), intent(in)  :: UTTAU, DT_ADD
    real(r8), intent(in), dimension(LPAR,IDBLK,JDBLK) :: BTEM
    !// Output
    real(r8), dimension(JVN_,LPAR,IDBLK,JDBLK,MPBLK),intent(out):: BJV
    !// Locals
    real(r8), dimension(JVN_,LPAR) :: JV_COL
    integer :: NCLDX(LPAR+1), NODT, NRAN, I, II, J, JJ, L, N
    real(r8) :: ODCLD(LPAR+1,NQD_), SWIJ(NQD_), PHOTAU
    logical :: LCLDRAN
    !// --------------------------------------------------------------------

    !// Random cloud cover?
    LCLDRAN = LCLDRANA .or. LCLDRANQ
    !// Random number for clouds
    NRAN  = mod(NCLDRAN-1,3) + 1

    !// Loop over latitude (J is global, JJ is block)
    do J = MPBLKJB(MP), MPBLKJE(MP)
       JJ    = J - MPBLKJB(MP) + 1

       !// Loop over longitude (I is global, II is block)
       do I = MPBLKIB(MP), MPBLKIE(MP)
          II    = I - MPBLKIB(MP) + 1

          !// FastJX; cloud cover
          do NODT = 1, NQD_
             SWIJ(NODT)  = SWSTORE(NODT,I,J)
             do L = 1, LPAR+1
                ODCLD(L,NODT) = CLDSTORE(L,NODT,I,J)
             end do
          end do
          do L = 1,LPAR+1
             NCLDX(L) = TYPSTORE(L,I,J)
          end do

          !// Time for J-values
          PHOTAU = UTTAU + DT_ADD / 3600._r8

          !// Get column values (replaces PHOTOL in UCI code)
          call jv_column(JV_COL, PHOTAU, SWIJ, ODCLD, BTEM(1,II,JJ), &
               NCLDX, LCLDRAN, NRAN, I, J, LPAR+1)

          !// Put values into BJV
          do L = 1, LPAR
             do N = 1, JVN_
                BJV(N,L,II,JJ,MP) = JV_COL(N,L)
             end do
          end do

       end do !// do I = MPBLKIB(MP), MPBLKIE(MP)
    end do !// do J = MPBLKJB(MP), MPBLKJE(MP)

    !// --------------------------------------------------------------------
  end subroutine jvalues_oslo
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine jv_column(JVALUES, PHOTAU, SWIJ, ODCLDQ, TEML, &
       NCLDX, LCLDRAN, NRAN, NSLON, NSLAT, LD)
    !// --------------------------------------------------------------------
    !//
    !// Oslo CTM3
    !// Sets photolysis rates in the column. Based on UCI subroutine PHOTOL.
    !//
    !// JVALUES: The J-values in the column.
    !// photau:  time at the end of time interval (hr)
    !//
    !// Ole Amund Sovde, October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LPAR
    use cmn_ctm, only: MPBLKIB, MPBLKIE, MPBLKJB, MPBLKJE
    use cmn_fjx, only: NQD_, JVN_, JIND
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// input
    integer, intent(in) :: NSLON, NSLAT, NRAN, LD
    integer, intent(in) :: NCLDX(LD)
    real(r8), intent(in) :: SWIJ(NQD_), ODCLDQ(LD,NQD_), PHOTAU, TEML(LPAR)
    logical, intent(in) :: LCLDRAN
    !// Output
    real(r8), intent(out) :: JVALUES(JVN_,LPAR)
    !// Locals
    integer :: jp, jl, n, fj_o1d, fj_o3
    real(r8) :: zpj(lpar,JVN_), zpjquad(lpar,JVN_), &
         zpjo(lpar), zpjquado(lpar), odcld(lpar+1), &
         frefl, u0, &            ! fraction of flux(energy-wtd) reflected
         sza
    !// --------------------------------------------------------------------

    !// Initialize
    JVALUES(:,:) = 0._r8

    !// Calculate photolysis rates for 4 quadrature atmospheres
    zpjquad(:,:) = 0._r8
    zpjquado(:)  = 0._r8

    if (LCLDRAN)  then
      do JL = 1, LD
        ODCLD(JL) = ODCLDQ(JL,NRAN)
      end do
      call PHOTOJ(PHOTAU,NSLON,NSLAT,SZA,U0,FREFL,ODCLD,TEML,NCLDX,ZPJ,ZPJO)

      do JP = 1, JVN_
        do JL = 1, LPAR
          ZPJQUAD(JL,JP) = ZPJQUAD(JL,JP) + ZPJ(JL,JP)
        end do
      end do

      !//Do not need this in CTM3
      !//do JL = 1,LPAR
      !//   ZPJQUADO(JL) = ZPJQUADO(JL) + ZPJO(JL)
      !// end do
    else

      do N = 1, NQD_
        if (SWIJ(N) .gt. 0._r8)  then

          do JL = 1,LD
            ODCLD(JL) = ODCLDQ(JL,N)
          end do
          !// ZPJO is not used and can be deleted, but then also in
          !// p-phot.f.
          !// Fast-JX, rate 2 (for O3 + hv), is for total O(3P) and O(1D)
          !// channels. That will be separated below.
          call PHOTOJ(PHOTAU,NSLON,NSLAT,SZA,U0,FREFL,ODCLD,TEML,NCLDX,ZPJ,ZPJO)

          do JP = 1,JVN_
            do JL = 1,LPAR
              ZPJQUAD(JL,JP) = ZPJQUAD(JL,JP) + SWIJ(N)*ZPJ(JL,JP)
            end do
          end do
          !// Do not need this in CTM3
          !// do JL = 1,LPAR
          !//   ZPJQUADO(JL) = ZPJQUADO(JL) + SWIJ(N)*ZPJO(JL)
          !// end do
        end if
      end do

    end if


    !// PHOTOJ has already mapped the FastJX indices to the listing indexing
    !// in ratj_oc.d, so we loop directly through the arrays:
    !// Fast-JX rate 2, for O3 + hv, is for total O(3P) and O(1D) channels.
    !// The O(1D) channel is rate number 3.
    do fj_o1d = 1, JVN_
      !// Find the rate index for the O(1D) channel.
      if (jind(fj_o1d) .eq. 3) exit
    end do

    do jp = 1, JVN_
      do jl = 1, LPAR
        !//rate_c(jl,nprkx(jp)) = zpjquad(jl,jp)
        !// NPRKX maps 1 to 1, 2 to 2, etc. NPRKX is set in phot_in (p-setc.f).

        !// Note the change in indices
        if (jind(jp) .eq. 2) then
          !// This value is to be subtracted from the Fast-JX rate 2,
          !// to get the O(3P) channel instead of the total.
          !print*,'subtracted J',jp,n,jind(n)
          JVALUES(jp,jl) = zpjquad(jl,jp) - zpjquad(jl,fj_o1d)
        else
          JVALUES(jp,jl) = zpjquad(jl,jp)
        end if
      end do
    end do

    !// --------------------------------------------------------------------
  end subroutine jv_column
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine update_chemistry(NDAYI,NDAY,NMET,NOPS,MONTH,LNEW_MONTH)
    !// --------------------------------------------------------------------
    !// Update chemistry, e.g. boundary conditions.
    !// Also initialises stratospheric H2O. This initialisation needs
    !// the LMTROP to be defined, and has to be inside NOPS-loop.
    !// Done outside parallell loop.
    !//
    !// Ole Amund Sovde, November 2014, October 2009, October 2008
    !// --------------------------------------------------------------------
    use cmn_precision, only: r8
    use cmn_size, only: LDUST, LOSLOCSTRAT
    use ch4routines, only: update_ch4surface
    use drydeposition_oslo, only: update_drydepvariables
    use dust_oslo, only: dust_globalupdate
    use psc_microphysics, only: LPSC, LAEROSOL
    use stratchem_oslo, only: read_oslo2d2
    use strat_aerosols, only: update_strat_backaer
    use aerosols2fastjx, only: update_tropaerosols4fjx
    use strat_h2o, only: strat_h2o_init, strat_h2o_init_clim
    use strat_o3noy_clim, only: get_strato3noy_clim
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in) :: NDAYI, NDAY, NMET, NOPS, MONTH
    logical, intent(in) :: LNEW_MONTH
    !// --------------------------------------------------------------------

    !// Update boundary conditions for stratospheric chemistry
    if (LNEW_MONTH) call read_oslo2d2(LNEW_MONTH)

    !// Update global DUST stuff
    if (LDUST .and. NOPS.eq.1 .and. NMET.eq.1) then
       call dust_globalupdate(NDAY)
    end if


    !// Set background aerosol surface area density for stratosphere
    if (LOSLOCSTRAT .and. LNEW_MONTH .and. NMET.eq.1. .and. NOPS.eq.1) &
         call update_strat_backaer(NDAY,LPSC.or.LAEROSOL,LNEW_MONTH)


    !// Get CTM climatology of tropospheric aerosol paths
    if (NMET.eq.1 .and. NOPS.eq.1) &
         call update_tropaerosols4fjx(MONTH, LNEW_MONTH)

    !// When running without stratospheric chemistry, read O3 climatology
    !// for use in the stratosphere. The fields read are used in oc_master.
    if (.not. LOSLOCSTRAT) &
         call get_strato3noy_clim(NDAY,NMET,NOPS,MONTH,LNEW_MONTH)


    !// Update surface CH4FIELD to be used in oc_setch4.
    if (LNEW_MONTH) call update_ch4surface()

    !// Update dry deposition for special cases
    call update_drydepvariables(LNEW_MONTH,NDAYI,NDAY,NMET,NOPS)


    !// Initialize species that needs LMTROP
    if (NDAY.eq.NDAYI .and.NMET.eq.1.and.NOPS.eq.1) then
       !// Initialize H2O: ECMWF in troposphere, from sumH2 in stratosphere
       !// Skipping if trsp_idx(114) is read from restart file!
       call strat_h2o_init()

       !// Initialize climatological H2O at tropopause
       call strat_h2o_init_clim()
    end if

    !// --------------------------------------------------------------------
  end subroutine update_chemistry
  !// ----------------------------------------------------------------------



  !// ------------------------------------------------------------------
end module main_oslo
!//=========================================================================

!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Ole Amund Sovde, February 2016
!//=========================================================================
!// Sea salt production module.
!//=========================================================================
module seasaltprod
  !// ----------------------------------------------------------------------
  !// MODULE: seasaltprod
  !// DESCRIPTION: Module for SEA SALT production.
  !//              Separate routine because it is used by both SALT
  !//              and BCOC packages.
  !//
  !// Contains SALT variables and the following routines:
  !//   subroutine seasalt_production
  !//   subroutine seasalt_production_maartenson03
  !//   subroutine seasalt_production_gantt15
  !//
  !// Ole Amund Sovde, February 2016
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: IPAR, JPAR, LPAR, IDBLK, JDBLK, MPBLK, NPAR_SALT
  use cmn_parameters, only: CPI, TK_0C
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter, private :: f90file = 'seasaltprod.f90'
  !// ----------------------------------------------------------------------
  save    !// All variables are to be saved.
  public  !// All is public
  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine seasalt_production(wind10,focean,AREAXY,rbinStart,rbinEnd, &
       rhosaltdry, I,J, PRODSEASALT)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Calculate production of seasalt particles at 80%RH at sea surface,
    !// and transfer this to kg of dry seasalt for a given bin with
    !// radius ranging from rbinStart to rbinEnd.
    !//
    !// Theory:
    !//   Gong et al, JGR vol 102, 1997, pp 3805
    !//   Smith et al, QJR met soc, vol 119,1993, pp 809-824
    !//   Fitzgerald, J. appl. met., vol 14, 1975, pp 1044-1049
    !//
    !// Original author: Alf Grini
    !//
    !// Adjusted to CTM3:
    !// Ole Amund Sovde, June 2015, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: wind10    !// 10m wind speed (m/s)
    real(r8), intent(in) :: focean    !// Fraction of salt water (0-1)
    real(r8), intent(in) :: AREAXY    !// Area of a grid cell (m2)
    real(r8), intent(in) :: rbinStart !// Radius of bin low limit (um)
    real(r8), intent(in) :: rbinEnd   !// Radius of bin high limit (um)
    real(r8), intent(in) :: rhosaltdry!// Density of dry salt particle [kg/m3]
    integer, intent(in) :: I,J        !// Not currently used

    !// Output
    real(r8), intent(out) :: prodseasalt !// Prod. of seasalt (kg/s) for bin

    !// Local
    integer  :: NN              ! Couting index for substep loop in the bin
    real(r8) :: saltpartweight  ! Weight of one particle (kg)
    real(r8) :: deltar          ! Length of bin (m)
    real(r8) :: flux1           ! Flux of small particles (#/m/s) 
    real(r8) :: flux2           ! Flux of large particles (spume) (#/m/s)
    real(r8) :: rsaltbin        ! radius of bin for which calc. are made (um)
    real(r8) :: Bsalt           ! Parameter in Monahan production
    real(r8) :: A1              ! Parameter in Smith production
    real(r8) :: A2              ! Parameter in Smith production
    real(r8) :: rbindry         ! radius of corresponding dry sea salt
    real(r8) :: Vpartdry        ! Volume of dry particle

    !// Local parameters
    integer, parameter :: Nsubsteps = 10 !Max iterations inside each bin
    real(r8), parameter :: &
         r01 = 2.1_r8, &    !Parameter in Smith formulation (um)
         r02 = 9.2_r8, &    !Parameter in Smith formulation (um)
         f1 = 3.1_r8, &     !Parameter in Smith formulation
         f2 = 3.3_r8        !Parameter in Smith formulation
    logical, parameter :: LDEBUG_SLT_PROD=.FALSE. !Debug /print statements
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'seasalt_production'
    !// --------------------------------------------------------------------
    !// Initializing for this bin
    PRODSEASALT = 0._r8

    !// Only sea salt over (open) ocean
    if (focean .le. 0._r8) return

    !// Parameters given by Smith et al 1993
    A1 = 10._r8**(0.0676_r8 * wind10 + 2.43_r8)
    A2 = 10._r8**(0.959_r8 * wind10**0.5_r8 - 1.476_r8)

    deltar = (rbinEnd - rbinStart) / real(Nsubsteps, r8)

    do NN = 1, Nsubsteps

       !// The radius of the "active" bin
       if (NN .eq. 1) then
          rsaltbin = rbinStart + 0.5_r8 * deltar
       else
          rsaltbin = rsaltbin + deltar
       end if

       !// Comment by Alf:
       !// In the following we ignore the dV(mix) which occurs when salt
       !// and water is mixed. Vsalt + Vwater + dV(mix) = newVolume.
       !// We set dV(mix)=0, but it certainly may be negative.
       !// I don't think this really comes out as an error in the
       !// calculations since we use Empirical functions for aerosol
       !// growth....

       !//Fitzgerald 1975 (Dry radius is half the radius at 80%RH) [um]
       rbindry = 0.5_r8 * rsaltbin 

       !// Volume of dry particle (ignoring change in total volume due
       !// to mixing) [m3]
       Vpartdry = rbindry**3 * 4._r8 / 3._r8 * CPI * 1.e-18_r8
     
       !// Weight of dry salt in particle [kg]
       saltpartweight = Vpartdry * rhosaltdry

       if (LDEBUG_SLT_PROD) then
          write(6,*)'rsaltbin, rbindry',rsaltbin,rbindry
          write(6,*)'Vpartdry,        ',Vpartdry
          write(6,*)'weight of dry part',saltpartweight
       end if

       !// CALCULATING DF/DR (USING Monahan for r<10 and Smith for r>10)    
       !// --------------------------------------------------------------
       if (rsaltbin .le. 7._r8) then
          !// MONAHAN, but following Gong et al 1997 suggestions
          !// -----------------------------------------------------------
          Bsalt = (0.380_r8 - log10(Rsaltbin)) / 0.650_r8

          !// Small particles
          flux1 = 1.37_r8 * wind10**3.41_r8 * Rsaltbin**(-3.0_r8) &
                  * (1.0_r8 + 0.057_r8 * Rsaltbin**1.05_r8) &
                  * 10.0_r8**(1.19_r8 * EXP(-1.0_r8 * BSALT * BSALT))

          !// Large (spume) particles, Monahan overestimates them
          flux2 = 0._r8

       else !// if (rsaltbin.gt.7.) then
          !// SMITH
          !// -----------------------------------------------------------
          flux1 = A1 * exp(-f1 * (log(rsaltbin / r01))**2)

          flux2 = A2 * exp(-f2 * (log(rsaltbin / r02))**2)

       end if  !Check on radius > 7um

       !// Converting to kg(dry salt)/sec from #/m2/sec/um
       flux1 = flux1 * AREAXY * saltpartweight * deltar
                  
       !// Converting to kg(dry salt)/sec from #/m2/sec/um
       flux2 = flux2 * AREAXY * saltpartweight * deltar

       !// Correction for percent of land in grid-cell
       Prodseasalt = Prodseasalt + &
            (flux1 + flux2) &
            * focean !// Correct for ocean/salt water fraction
                  
    end do !// do NN = 1, Nsubsteps

    !// --------------------------------------------------------------------
  end subroutine seasalt_production
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_production_martensson03(wind10,focean,AREAXY, &
       rbinStart,rbinEnd, rhosaltdry, SFT, I,J, PRODSEASALT)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Calculate production of seasalt particles at 80%RH and transfer
    !// this to kg of dry seasalt, for a given size bin.
    !//
    !// Theory:
    !//   Mårtensson et al, JGR vol 108, no D9, doi:10.1029/2002JD002263
    !//   Gong et al, JGR vol 102, 1997, pp 3805
    !//   Smith et al, QJR met soc, vol 119,1993, pp 809-824
    !//   Fitzgerald, J. appl. met., vol 14, 1975, pp 1044-1049
    !//
    !// Note that martensson03 use dry diameter (m), but both monahan86 and
    !// smith93 use radius (um) at 80 % RH in their parameterizations.
    !// If you mix this up, you are in big trouble.
    !//
    !// Original author: Alf Grini
    !//
    !// Ole Amund Sovde, June 2015, September 2009
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: wind10    !// 10m wind speed (m/s)
    real(r8), intent(in) :: focean    !// Fraction of salt water (0-1)
    real(r8), intent(in) :: AREAXY    !// Area of a grid cell (m2)
    real(r8), intent(in) :: rbinStart !// Radius of bin low limit (um)
    real(r8), intent(in) :: rbinEnd   !// Radius of bin high limit (um)
    real(r8), intent(in) :: rhosaltdry!// Density of dry salt particle [kg/m3]
    real(r8), intent(in) :: SFT       !// Surface temperature (K)
    integer, intent(in) :: I,J        !// Not currently used

    !// Output
    real(r8), intent(out) :: prodseasalt !// Prod. of seasalt (kg/s) for bin

    !// LOCAL
    integer  :: NN  !// Counting index for substeps
    real(r8) :: &
         saltpartweight, &  !// Weight of one particle (kg)
         deltar, &          !// Length of bin (m)
         flux1, &           !// Flux due to Martensson formulation (1/m/s) 
         flux2, &           !// Flux due to Smith formulation   (1/m/s)
         rsaltbin, &        !// radius of bin for which calc. are made (um)
         Bsalt, &           !// Parameter in Monahan production
         A1, &              !// Parameter in Smith production
         A2, &              !// Parameter in Smith production
         rbindry, &         !// radius of corresponding dry sea salt
         Vpartdry           !// Volume of dry particle

    !// Variables used in martensson03
    real(r8) :: &
         Ak,Bk, &        !// Parameters in martensson parameterizaion    
         W, &            !// Fraction of ocean covered by whitecaps
         phi, &          !// Production of aerosols in whitecaps 
         Dbindry, &      !// size (um) of dry seasaltaerosol
         dlog10D, &      !// logarithmic size increment of bin
         sft_bnd         !// Limited value of surface temp (K)
    integer :: k         !// variable deciding size interval in martensson03

    !// LOCAL PARAMETERS
    integer, parameter :: Nsubsteps = 10  !// Max iterations inside each bin
    real(r8), parameter :: &
         r01 = 2.1_r8, &    !Parameter in Smith formulation (um)
         r02 = 9.2_r8, &    !Parameter in Smith formulation (um)
         f1 = 3.1_r8, &     !Parameter in Smith formulation
         f2 = 3.3_r8        !Parameter in Smith formulation
    logical, parameter :: LDEBUG_SLT_PROD2 = .FALSE. !// Debug /print statements

    real(r8)             :: c(5,3)         !// Parameters for martensson03
    real(r8)             :: d(5,3)         !// Parameters for martensson03
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'seasalt_production_martensson03'
    !// --------------------------------------------------------------------

    !// Data for c(k): 3 size intervals and 5 parameters martensson03, table 1
    c(:,1) = (/-2.881e6_r8,-3.003e13_r8,-2.867e21_r8, 5.932e28_r8,-2.576e35_r8/)
    c(:,2) = (/-6.743e6_r8, 1.183e14_r8,-8.148e20_r8, 2.404e27_r8,-2.452e33_r8/)
    c(:,3) = (/ 2.181e6_r8,-4.165e12_r8, 3.132e18_r8,-9.841e23_r8, 1.085e29_r8/)

    !// Data for d(k): 3 size intervals and 5 parameters martensson03, table 2
    d(:,1) = (/ 7.609e8_r8, 1.829e16_r8, 6.791e23_r8,-1.616e31_r8, 7.188e37_r8/)
    d(:,2) = (/ 2.279e9_r8,-3.787e16_r8, 2.528e23_r8,-7.310e29_r8, 7.368e35_r8/)
    d(:,3) = (/-5.800e8_r8, 1.105e15_r8,-8.297e20_r8, 2.601e26_r8,-2.859e31_r8/)

    !// ------------------------------------------------------------------

    !// Initialize
    PRODSEASALT = 0._r8

    !// Only sea salt over (open) ocean
    if (focean .le. 0._r8) return

    !// Limit surface temperature to 298K for which martensson is useful
    sft_bnd = min(SFT, 298._r8)

    !// Parameters given by Smith et al 1993
    A1 = 10._r8**(0.0676_r8 * wind10 + 2.43_r8)
    A2 = 10._r8**(0.959_r8 * wind10**0.5_r8 - 1.476_r8)

    !// Fraction of ocean covered with whitecaps
    !// Note that Martenson et al (2003) (eq.2) give this in PERCENT
    W = 3.84e-4_r8 * wind10**3.41_r8 * 1.e-2_r8

    deltar = (rbinEnd - rbinStart) / real(Nsubsteps, r8)

    do NN = 1, Nsubsteps

       !// The radius of the "active" bin
       !// Here d(log10(r)) is calculated at the centers of the bins,
       !// which means that the first bin needs to be treated somewhat
       !// specially.
       !// For other bins, d(log10(r)) is calculated at the current bin
       !// using radius at previous bin, i.e. it is calculated before
       !// rsaltbin is updated.
       if (NN .eq. 1) then
          rsaltbin = rbinStart + 0.5_r8 * deltar
          !dlog10D = log10((rbinStart + deltar) / rsaltbin)
          !// Assume same increment at center (rsaltbin) as from
          !// rbinStart to rbinEnd
          dlog10D = log10((rbinStart + deltar) / rbinStart)
       else
          dlog10D = log10((rsaltbin + deltar) / rsaltbin)
          rsaltbin = rsaltbin + deltar
       end if


       if (LDEBUG_SLT_PROD2) then
          write(6,*)f90file//':'//subr//': substep',NN,rbinStart,rbinEnd
       end if

     
       !// Comment by Alf:
       !// In the following we ignore the dV(mix) which occurs when salt
       !// and water is mixed: Vsalt + Vwater + dV(mix) = newVolume.
       !// We set dV(mix)=0, but it certainly be negative.
       !// I don't think this really comes out as an error in the
       !// calculations since we use Empirical functions for aerosol
       !// growth....

       !// Fitzgerald 1975 (Half the radius at 80%RH)
       rbindry = 0.5_r8 * rsaltbin 

       !// The dry diameter is 2*dry radius (um)
       Dbindry = 2._r8 * rbindry

       !// Volume of dry particle (ignoring change in total volume due
       !// to mixing)
       Vpartdry = rbindry**3 * 4._r8 / 3._r8 * CPI * 1.e-18_r8

       !// Weight of dry salt in particle
       saltpartweight = Vpartdry * rhosaltdry

       if (LDEBUG_SLT_PROD2) then
          write(6,*)'rsaltbin, rbindry',rsaltbin,rbindry
          write(6,*)'Vpartdry,        ',Vpartdry
          write(6,*)'weight of dry part',saltpartweight
       end if

       !// CALCULATING DF/DR (USING martensson for rdry<1.4 and
       !// Smith for r>1.4)    
       !// -----------------------------------------------------------
       !// martensson03 for Dsaltdry smaller than...
       !// -----------------------------------------------------------
       IF (rbindry .le. 1.4_r8) THEN
          !// Martensson03 is valid up to Ddry=2.8 um

          !// Find out what regime you are in: (what is k in equation 5)
          if (Dbindry .lt. 0.145_r8) then
             k = 1
          else if (Dbindry .lt. 0.419_r8) then
             k = 2
          else if (Dbindry .lt. 2.8001_r8) then
             k = 3
          else
             write(6,'(a)') f90file//':'//subr// &
                  ': error in production'
             stop
          end if

          !// Get Ak in equation 5 (correct for um-->m)
          Ak = c(5,k) * Dbindry**4 * 1.e-24_r8 &
               + c(4,k) * Dbindry**3 * 1.e-18_r8 &
               + c(3,k) * Dbindry**2 * 1.e-12_r8 &
               + c(2,k) * Dbindry * 1.e-6_r8 &
               + c(1,k)

          !// Get Bk in equation 5 (correct for um-->m)
          Bk = d(5,k) * Dbindry**4 * 1.e-24_r8 &
               + d(4,k) * Dbindry**3 * 1.e-18_r8 &
               + d(3,k) * Dbindry**2 * 1.e-12_r8 &
               + d(2,k) * Dbindry * 1.e-6_r8 &
               + d(1,k)

          !// Get phi (=dF/dlog10D(dry)) for all whitecaps
          !// 
          phi = Ak * sft_bnd + Bk

          !// Multiply flux with fraction of whitecaps
          flux1 = W * phi * AREAXY * saltpartweight * dlog10D

          !// Maartensson03 does not have a flux2
          flux2 = 0._r8

       else if (rsaltbin .le. 7._r8) then
          !// monahan is valid approx. up to r=7um
          !// --------------------------------------------------------
          !// MONAHAN
          !// --------------------------------------------------------
          !// Get Bsalt in Monahan formulation
          Bsalt = (0.380_r8 - log10(Rsaltbin))/0.650_r8   
     
          !// Small particles (kg dry salt)
          flux1 = 1.37_r8 * wind10**3.41_r8 * Rsaltbin**(-3.0_r8) &
                  * (1.0_r8 + 0.057_r8 * Rsaltbin**1.05_r8) &
                  * 10.0_r8**(1.19_r8 * exp(-1.0_r8 * BSALT * BSALT)) &
                  * AREAXY * saltpartweight * deltar

          !// Large (spume) particles, Monahan overestimates them
          flux2 = 0._r8

       else
          !// 80 % radius is larger than 7 um
          !// --------------------------------------------------------
          !// SMITH
          !// --------------------------------------------------------

          !// Get flux1 in kg(dry salt)
          Flux1 = A1 * exp(-f1 * (log(rsaltbin / r01))**2) &
                  * AREAXY * saltpartweight * deltar

          Flux2 = A2 * exp(-f2 * (log(rsaltbin / r02))**2) &
                  * AREAXY * saltpartweight * deltar

       end if !// Check on radius sizes 

       !// No production if ocean is frozen (i.e. sft<273K)
       if (sft_bnd .lt. TK_0C) then
          flux1 = 0._r8
          flux2 = 0._r8
       end if

       !// Sum up the production for this bin
       Prodseasalt = Prodseasalt &
            + (Flux1 + flux2) * focean ! Correct for ocean/salt water fraction

    end do !// do NN = 1, Nsubsteps
         
    !// --------------------------------------------------------------------
  end subroutine seasalt_production_martensson03
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_production_gantt15(wind10,fsaltwater,AXY,rbinStart, &
       rbinEnd, rhosaltdry, surfaceT, I,J, PRODSEASALT)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Calculate production of seasalt particles at 80%RH at sea surface,
    !// and transfer this to kg of dry seasalt for a given bin with
    !// radius ranging from rbinStart to rbinEnd.
    !//
    !// According to Gantt et al (GMD, 2015, doi:10.5194/gmd-8-619-2015).
    !// Follows Gong et al (Global Biogeochem. Cy., 2003,
    !//         doi:10.1029/2003GB002079)
    !// with SST dependence as in Jaegle et al (ACP, 2011,
    !//         doi:10.5194/acp-11-3137-2011).
    !//
    !// Ole Amund Sovde, June 2015
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: wind10    !// 10m wind speed (m/s)
    real(r8), intent(in) :: fsaltwater!// Fraction of salt water (0-1)
    real(r8), intent(in) :: AXY       !// Area of a grid cell (m2)
    real(r8), intent(in) :: rbinStart !// Radius of bin low limit (um)
    real(r8), intent(in) :: rbinEnd   !// Radius of bin high limit (um)
    real(r8), intent(in) :: rhosaltdry!// Density of dry salt particle [kg/m3]
    real(r8), intent(in) :: surfaceT  !// Surface temperature [K]
    integer, intent(in) :: I,J        !// Not currently used

    !// Output
    real(r8), intent(out) :: prodseasalt !// Prod. of seasalt (kg/s) for bin

    !// Local
    integer  :: NN              ! Couting index for substep loop in the bin
    real(r8) :: saltpartweight  ! Weight of one particle [kg]
    real(r8) :: dr_subbin       ! Length of sub-bin [um]
    real(r8) :: rsubbin         ! Radius of actual sub-bins [um]
    real(r8) :: flux            ! Flux [particles/(m2*s*um)]
    real(r8) :: A,B             ! Parameters in Gantt production
    real(r8) :: sfcTC           ! Surface temperature in Celcius
    real(r8) :: jsst            ! Jaegle sea surface temperature factor [0-1]
    real(r8) :: rdry            ! Radius of corresponding dry sea salt [um]
    real(r8) :: Vpartdry        ! Volume of dry particle [m3]

    !// Local parameters
    integer, parameter :: Nsubsteps = 10 !Max iterations inside each bin
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'seasalt_production_gantt15'
    !// --------------------------------------------------------------------
    !// Initializing for this bin
    PRODSEASALT = 0._r8

    !// Only sea salt over (open) ocean
    if (fsaltwater .le. 0._r8) return


    !// Emissions of sea spray aerosols [kg/s]
    !// According to Gong et al (2003) with SST dependence as in
    !// Jaegle et al (2011).
    sfcTC = surfaceT - TK_0C !// Temperature in Celsius
    jsst = 0.3_r8 + 0.1_r8 * sfcTC &
           - 0.0076_r8 * sfcTC**2 &
           + 0.00021_r8 * sfcTC**3
    !// Do not allow negative jsst (for about sfcTC<-3.5);
    !// that means no production, i.e. return from routine
    if (jsst .lt. 0._r8) return


    !// dr for each substeps [um]
    dr_subbin = (rbinEnd - rbinStart) / real(Nsubsteps, r8)

    do NN = 1, Nsubsteps

       !// The radius sub-bin [um]
       if (NN .eq. 1) then
          rsubbin = rbinStart + 0.5_r8 * dr_subbin
       else
          rsubbin = rsubbin + dr_subbin
       end if

       !// Gong et al (2003) coefficients
       A = 4.7_r8 * (1._r8 &
            + 30._r8 * rsubbin)**(-0.017_r8 * rsubbin**(-1.44_r8))
       B = (0.433_r8 - log10(rsubbin)) / 0.433_r8

       !// Gong et al (2003) dF/dr80 has units particles/(m2*s*um)
       flux = 1.373_r8 * wind10**3.41_r8 * rsubbin**(-A) &
            * (1._r8 + 0.057_r8 * rsubbin**3.45_r8) &
            * 10._r8**(1.607_r8 * exp(-(B**2)))

       !// flux has units [#/m2/sec/um], which we convert
       !// to [kg/s]. Sea salt code calculates dry mass, which
       !// may present some inconsistency, but I think this
       !// may work. Gantt uses a volume and apparent density,
       !// but in the end we need the mass of sea salt.

       !//Fitzgerald 1975 (dry radius half the 80%RH radius) [um]
       rdry = 0.5_r8 * rsubbin 

       !// Volume of dry particle (ignoring change in total volume due
       !// to mixing) [m3] (converted from um)
       Vpartdry = rdry**3 * 4._r8 / 3._r8 * CPI * 1.e-18_r8
     
       !// Weight of dry salt in particle [kg]
       saltpartweight = Vpartdry * rhosaltdry

       !// From particles/(m2*s*um) to kg/s
       !// First step: sum up as kg/(m2*s), multiply by constants
       !// at the end of Nsubsteps loop.
       prodseasalt = prodseasalt + &
            + flux * saltpartweight * dr_subbin
                  
    end do !// do NN = 1, Nsubsteps

    !// Take JSST, AXY and fsaltwater into account
    !// (i.e. conversion from kg/(m2*s) to kg/s)
    prodseasalt = prodseasalt * jsst * AXY * fsaltwater


    !// --------------------------------------------------------------------
  end subroutine seasalt_production_gantt15
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine seasalt_production_witek16(wind10,fsaltwater,AXY,rbinStart, &
       rbinEnd, rhosaltdry, surfaceT, I,J, PRODSEASALT)
    !// --------------------------------------------------------------------
    !// PURPOSE:
    !// Calculate production of seasalt particles at 80%RH at sea surface,
    !// and transfer this to kg of dry seasalt for a given bin with
    !// radius ranging from rbinStart to rbinEnd.
    !//
    !// According to Witek et al (JGR, 2016, doi:10.1002/2015JD023726).
    !// Follows Sofiev et al (JGR, 2011, doi:10.1029/2010JD014713)
    !// omitting salinity effects and with SST dependence from
    !// Jaegle et al (ACP, 2011, doi:10.5194/acp-11-3137-2011).
    !//
    !// Ole Amund Sovde, February 2016
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    real(r8), intent(in) :: wind10    !// 10m wind speed (m/s)
    real(r8), intent(in) :: fsaltwater!// Fraction of salt water (0-1)
    real(r8), intent(in) :: AXY       !// Area of a grid cell (m2)
    real(r8), intent(in) :: rbinStart !// Radius of bin low limit (um)
    real(r8), intent(in) :: rbinEnd   !// Radius of bin high limit (um)
    real(r8), intent(in) :: rhosaltdry!// Density of dry salt particle [kg/m3]
    real(r8), intent(in) :: surfaceT  !// Surface temperature [K]
    integer, intent(in) :: I,J        !// Not currently used

    !// Output
    real(r8), intent(out) :: prodseasalt !// Prod. of seasalt (kg/s) for bin

    !// Local
    integer  :: NN              ! Couting index for substep loop in the bin
    real(r8) :: &
         saltpartweight, &  ! Weight of one particle [kg]
         dr_subbin, &       ! Length of sub-bin [um]
         rsubbin, &         ! Radius of actual sub-bins [um]
         flux25, &          ! Flux at 25C
         B, &               ! Parameters in Gantt production
         WU10, &            ! Fraction of ocean covered by whitecaps
         sfcTC, &           ! Surface temperature in Celcius
         jsst, &            ! Jaegle sea surface temperature factor [0-1]
         rdry, &            ! Radius of corresponding dry sea salt [um]
         Vpartdry           ! Volume of dry particle [m3]

    !// Local parameters
    integer, parameter :: Nsubsteps = 10 !Max iterations inside each bin
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'seasalt_production_witek16'
    !// --------------------------------------------------------------------
    !// Initializing for this bin
    PRODSEASALT = 0._r8

    !// Only sea salt over (open) ocean
    if (fsaltwater .le. 0._r8) return


    !// Emissions of sea spray aerosols [kg/s]
    !// According to Gong et al (2003) with SST dependence as in
    !// Jaegle et al (2011).
    sfcTC = surfaceT - TK_0C !// Temperature in Celsius
    jsst = 0.3_r8 + 0.1_r8 * sfcTC &
           - 0.0076_r8 * sfcTC**2 &
           + 0.00021_r8 * sfcTC**3
    !// Do not allow negative jsst (for about sfcTC<-3.5);
    !// that means no production, i.e. return from routine
    if (jsst .lt. 0._r8) return


    !// Fraction of ocean covered with whitecaps
    !// Note that Martenson et al (2003) (eq.2) give this in PERCENT
    WU10 = 3.84e-6_r8 * wind10**3.41_r8


    !// dr for each substeps [um]
    dr_subbin = (rbinEnd - rbinStart) / real(Nsubsteps, r8)

    do NN = 1, Nsubsteps

       !// The radius sub-bin [um]
       if (NN .eq. 1) then
          rsubbin = rbinStart + 0.5_r8 * dr_subbin
       else
          rsubbin = rsubbin + dr_subbin
       end if

       !// 

       !// Sofiev et al (2011)
       B = (0.27_r8 - log10(rsubbin)) / 1.1_r8
       !// Flux (dF0/dr80)_(SST=25C) has units particles/(m2*um) during
       !// whitecap decay time (3.53s)
       flux25 = 3.53e6_r8 * exp(-0.09_r8/(rsubbin + 0.003_r8)) &
                         / (2._r8 + exp(-5._r8/rsubbin)) &
               * (1._r8 + 0.05_r8*rsubbin**1.05_r8) / rsubbin**3 &
               * 10._r8**(1.05_r8*exp(-B**2))

       !// Fhe final flux is
       !// flux = WU10 / 3.53_r8 * flux25 * jsst
       !//        [0-1]   [s]      #/(m2*um)  = #/(m*s*um)
       !// Will do multiplication with WU10*jsst/3.53 at the end since
       !// they are constant for substep

       !// flux has units [#/m2/sec/um], which we convert
       !// to [kg/s]. Sea salt code calculates dry mass, which
       !// may present some inconsistency, but I think this
       !// may work. Gantt uses a volume and apparent density,
       !// but in the end we need the mass of sea salt.

       !//Fitzgerald 1975 (dry radius half the 80%RH radius) [um]
       rdry = 0.5_r8 * rsubbin 

       !// Volume of dry particle (ignoring change in total volume due
       !// to mixing) [m3] (converted from um)
       Vpartdry = rdry**3 * 4._r8 / 3._r8 * CPI * 1.e-18_r8
     
       !// Weight of dry salt in particle [kg]
       saltpartweight = Vpartdry * rhosaltdry

       !// From particles/(m2*s*um) to kg/s
       !// First step: sum up as kg/(m2*s), multiply by constants
       !// at the end of Nsubsteps loop.
       prodseasalt = prodseasalt + &
            + flux25 * saltpartweight * dr_subbin
                  
    end do !// do NN = 1, Nsubsteps

    !// Take JSST, WU10, decay time, AXY and fsaltwater into account
    !// (i.e. conversion from kg/(m2*s) to kg/s)
    prodseasalt = prodseasalt * WU10 * jsst/3.53_r8 * AXY * fsaltwater


    !// --------------------------------------------------------------------
  end subroutine seasalt_production_witek16
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
end module seasaltprod
!//=========================================================================

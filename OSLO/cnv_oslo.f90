!//=========================================================================
!// Oslo CTM3
!//=========================================================================
!// Based on UCI CTM core p-7.1 (1/2013).
!//
!// Ole Amund Sovde, April 2015
!//=========================================================================
!// Routines for oslo convective plume calculations.
!//=========================================================================
module cnv_oslo
  !// ----------------------------------------------------------------------
  !// MODULE: cnv_oslo
  !// DESCRIPTION: Routines for oslo convective plume calculations.
  !// ----------------------------------------------------------------------
  !// Oslo wet removal due to convective rain is treated similarly to CTM2,
  !// although there are some differences.
  !// The method consists of two parts;
  !// 1. Calculating the fraction of convective rain to cloudwater in the
  !//    elevator (QFRAC)
  !// 2 . Calculating the fraction of tracer in cloudwater (DISSOLVEDFRAC)
  !//
  !// For a tracer, the fraction removed is then
  !//       WETL(L) =           !// Fraction of tracer removed
  !//             DISSOLVEDFRAC !// Fraction of tracer in cloudwater
  !//             * QFRAC       !// Fraction of cloudwater removed by
  !//                           !// convective rain
  !//
  !// - QFRAC is global, depending on meteorological data.
  !// - DISSOLVEDFRAC is calculated based on the liquid water volume
  !//   concentration, i.e. the ratio of elevator liquid water volume to
  !//   elevator volume.
  !//
  !// Changes from CTM2 code:
  !// - FAQ was 4D, but we only need column values in the new treatment, so
  !//   FAQ has been removed.
  !// - Fluxes are treated slightly differently.
  !// - In the input (p-wind_ec.f) they are filtered for small values.
  !// - The lifting/condensation/rainout/entrainment/detrainment is treated
  !//   in two steps, similarly to the convective transport; first entrainment
  !//   (ENT_U, or CENTU) and then detrainment is calculated to balance the
  !//    mass flux (FLUX_E[L+1] - (FLUX_E[L] + ENT_U[L])). Both entrainment and
  !//    detrainment are assumed to possibly be positive or negative, although
  !//    the first should be non-negative due to read-in, while the second may
  !//    be both.
  !// - Also the liquid water mass in the elevator is allowed to be detrained.
  !//
  !// Contains
  !//   - subroutine elevator_fractions
  !//   - subroutine liquid_fractions
  !//   - subroutine wf_henry
  !//
  !// Ole Amund Sovde, April 2016, September 2014, April-November 2009
  !// ----------------------------------------------------------------------
  use cmn_precision, only: r8
  use cmn_size, only: LPAR, LWEPAR, IDBLK, JDBLK, MPBLK
  use cmn_ctm, only: ETAA, ETAB
  use cmn_met, only: Q, P
  use cmn_oslo, only: QFRAC, LW_VOLCONC, LELEVTEMP
  !// ----------------------------------------------------------------------
  implicit none
  !// ----------------------------------------------------------------------
  character(len=*), parameter :: f90file = 'cnv_oslo.f90'
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------

  private
  public elevator_fractions, liquid_fractions

  !// ----------------------------------------------------------------------

contains

  !// ----------------------------------------------------------------------
  subroutine elevator_fractions(I,J,II,JJ,MP, LE, &
          BTEM,CNVRAIN, ENTU_DT, FLUXE_DT, DTCONV)
    !// --------------------------------------------------------------------
    !// Compute
    !// 1. volume fraction of liquid water in the convective plume (elevator),
    !//    i.e. volume of liquid water to total elevator volume: LW_VOLCONC
    !// 2. fraction of precipitated mass to liquid mass in elevator: QFRAC
    !//
    !// Important: No new precipitation is formed in this scheme! This should
    !//            perhaps be revised. Fluxes and rain are accumulated fields,
    !//            and should be consistent. If new rain is to be included,
    !//            then rain and cloud water should be treated separately.
    !//            However, convective rain and mass fluxes should be fairly
    !//            consistent, so produced rain due to condensation should
    !//            probably already be accounted for.
    !// IMPORTANT: UNITS: Units must match!
    !//                   CNVRAIN: [kg/s]
    !//                   ENT_U:   [kg/s]
    !//                   FLUX_E:  [kg/s]
    !//
    !// The volume of liquid water in the elevator is
    !//   V_LW_elevator = liquid_water_mass / rho_LW
    !// where rho_LW is the density of water (which is 10^3kg/m3).
    !//
    !// The volume concentration is then
    !//   LW_VOLCONC = V_LW_elevator / V_elevator
    !//
    !// QFRAC and LW_VOLCONC are global arrays, of size (LPAR,IDBLK,JDBLK,MPBLK)
    !// to make them quickly accessible in IJ-blocks.
    !//
    !// For the convective lifitng, we follow the air from the surface and take
    !// entrainment/detrainment and mass fluxes into account.
    !//
    !// Mass flux in + entrainment = mass flux out + detrainment.
    !//   F_e(L) + E_u(L) = F_e(L+1) + D_u(L)
    !// E_u is non-negative, and D_u is calculated to balance the mass flux
    !// F_e in and out of each grid box. Here we calculate
    !//   -D_u(L) =  F_e(L+1) - (F_e(L) + E_u(L))
    !// and name it an additional entrainment.
    !//
    !// Outline:
    !//   1. Find mass and elevator properties coming up
    !//   2. Find entrained mass and its properties (also extra entrainent to
    !//      balance mass flux)
    !//   3. Mix updraft and entrainment, condense water in the mix due to
    !//      change in temperature
    !//   4. Possibly treat evaporation
    !//   5. Find amount of liquid water condensed
    !//      Add to the liquid water flux from below
    !//   6. Treat detrainment
    !//   7. Remove if net rain out; skip possible evaporation
    !//   8. Move rest to next level
    !//
    !// Ole Amund Sovde, November 2009, Updated November 2011
    !// --------------------------------------------------------------------
    use cmn_parameters, only: R_AIR, CP_AIR, Lv_0C
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer , intent(in) :: I, J, &  !// Global indices of the column
                            II, JJ, &!// IJ-block indices of the column
                            MP,&     !// IJ-block index
                            LE       !// Max height level of convection

    real(r8) , intent(in)  :: &
          BTEM(LPAR,IDBLK,JDBLK), &  !// Temperature [K]
          CNVRAIN(LWEPAR), &         !// Convective precipitiation [kg/s]
          FLUXE_DT(LPAR+1), &        !// Mass fluxes during time step [kg]
          ENTU_DT(LPAR+1), &         !// Entrainment during time step [kg]
          DTCONV                     !// Time step [s]

    !// Output
    !// No ouput to be returned; routine will update the global fields
    !// available in oc_globalvariables.f:
    !//   QFRAC: Fraction conv. rain/Cloudwater mass in elevator
    !//   LW_VOLCONC: elevator liquid water volume (m3) / elevator volume (m3)


    !// Local paramters
    real(r8) , parameter :: &
          TLIM = 0.01_r8, &      !// Minimum dT allowed due to condensation
          Rd = R_AIR, &          !// Gas constant for dry air.
          EPSILON = 0.622_r8, &  !// Epsilon = Rd/Rv
          CP = CP_AIR, &         !// Specific heat capacity for air at const. p
          Lw = Lv_0C, &          !// Latent heat for water.
          alpha_LW = 1.e-3_r8, & !// Spec. volume (1/density) liq.water[m3/kg]
          firstguess = 0.2_r8, & !// Assume 20% condensation as a first guess.
          secondguess = 0.2_r8   !// Assume 20% also in each iteration step.

    !// Local variables
    real(r8), dimension(LPAR+1) :: FLUX_E, ENT_U !// Mass fluxes [kg/s]

    real(r8), dimension(LWEPAR) :: &
          pres, &               !// Center pressure [hPa]
          ambient_t, &          !// Environment temperature [K]
          ambient_q, &          !// Environment specific humididy [kg/kg]
          DET_U, &              !// Detrained air mass [kg/s]
          cprecip, &            !// Convective rain formed at each level [kg/s]
          elev_mass_air, &      !// Mass of air in elevator [kg/s]
          elev_mass_vapor, &    !// Mass of vapor in elevator [kg/s]
          elev_theta, &         !// Potential temperature in elevator [K]
          elev_t, &             !// Temperature in elevator [K]
          elev_mass_lw          !// Liquid water mass in elevator [kg/s]


    real(r8) :: &
          elev_vol, &                      !// Elevator volume
          elev_q, &                        !// Elevator specific humidity
          ma_vapor, ma_air, &              !// Intermediate mass numbers
          cond_water, &                    !// condensed water [kg]
          entrainment, &  !// Entrained air
          detrainment, &  !// Detrined air
          air_up, &       !// Air coming up from below (updraft)
          vap_up, &       !// Vapor coming up from below (updraft)
          evap_rain, &    !// Evaproated rain
          elev_mmr_lw, &  !// Elevator mass mixing ratio of liquid water
          liq_in, &       !// Liquid coming in from above (rain)
          rainout, &      !// Rainout (net rain out of bottom)
          rin_top, &      !// Rain in at gridbox top
          tmpv            !// Temporary variable

    integer :: &
          L, LL, &                    !// Loop indices
          LCONVMAX,LCONVMAX1,LCONVMAX2, & !// Highest levels of convection
          LCONVMIN,LCONVMIN1,LCONVMIN2    !// Lowest levels of convection
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'elevator_fractions'
    !//---------------------------------------------------------------------


    !// Initialize Output for this column
    QFRAC(:,II,JJ,MP)      = 0._r8
    LW_VOLCONC(:,II,JJ,MP) = 0._r8

    LELEVTEMP(:,II,JJ,MP) = 0

    !// If no entrainment and no flux, there is no convection and we return.
    !// Both ENTU_DT and FLUXE_DT are positive by definition.
    if (maxval(ENTU_DT) .eq. 0._r8 .and. maxval(FLUXE_DT) .eq. 0._r8) return

    !// Convert from [kg] to [kg/s]
    ENT_U(:) = ENTU_DT(:) / DTCONV
    FLUX_E(:) = FLUXE_DT(:) / DTCONV

    !// Initialize local variables; all have units [kg/s]
    cprecip(:)          = 0._r8 !// Convective precipitation
    elev_mass_air(:)    = 0._r8 !// Air mass in elevator
    elev_mass_vapor(:)  = 0._r8 !// Water vapor mass in elevator
    elev_theta(:)       = 0._r8 !// Potential temperature
    elev_t(:)           = 0._r8 !// Temperatur in elevator
    elev_mass_lw(:)     = 0._r8 !// Liquid water content (in elevator)


    !// Calculate detrainment to balance mass fluxes
    !// F_E(L) + E(L) = F_E(L+1) + D
    DET_U(:) = 0._r8
    do L = 1, LWEPAR-1
       DET_U(L) = FLUX_E(L) + ENT_U(L) - FLUX_E(L+1)
    end do
    !// Remember that flux/entrainment/detrainment balance is given by
    !//   ENT_U(L) + FLUX_E(L) = FLUX_E(L+1) - DET_U(L)
    !//   DET_U(L) = FLUX_E(L+1) - FLUX_E(L) - ENT_U(L)
    !// If the balance requires degative detrainment, we have ADDITIONAL
    !// entrainment: extra_entrainment = -DET_U(L)


    !// Set vertical arrays from global arrays
    do L = 1, LWEPAR-1
       ! Specific humidity
       ambient_q(L)  = Q(I,J,L)
       !// Convective precipitaiton (not the same as PCNV in CTM2)
       cprecip(L) = CNVRAIN(L)
       !// Set vertical arrays from IJ-block arrays
       ambient_t(L)  = BTEM(L,II,JJ)
       !// Set center pressure [hPa]
       PRES(L)  = (ETAA(L) + ETAA(L+1) + P(I,J)*(ETAB(L) + ETAB(L+1)))*0.5_r8
    end do

    !// Check that we have rain in the column, if not leave
    if (maxval(abs(cprecip)).eq.0._r8) return
    !// The net sum should not be negative; we cannot evaporate more than
    !// is generated.
    if (sum(cprecip) .lt. 0._r8) return


    !// At what level does convection start?
    !// Find lowest level of convection for entrainment (can be layer 1)
    LCONVMIN1 = LE
    do L = 1, LE
       if (ENT_U(L) .gt. 0._r8) then
          LCONVMIN1 = L
          exit
       end if
    end do
    !// Find lowest level of the total mass flux in updrafts
    !// (level 1 is zero, so we start at 2)
    LCONVMIN2 = LE
    do L = 2, LE
       if (FLUX_E(L) .gt. 0._r8) then
          LCONVMIN2 = L
          exit
       end if
    end do


    !// How high does convection start?
    !// Find highest level of convection for entrainment
    LCONVMAX1 = 1
    do L = LE, 1, -1
       if (ENT_U(L) .ne. 0._r8) then
          LCONVMAX1 = L
          exit
       end if
    end do
    !// Find highest level of the total mass flux in updrafts
    LCONVMAX2 = 1
    do L = LE, 1, -1
       if (FLUX_E(L) .ne. 0._r8) then
          !// flux into this box, so we need to calculate for this level
          LCONVMAX2 = L
          exit
       end if
    end do

    !// Find min/max level for calculations
    LCONVMIN = min(LCONVMIN1,LCONVMIN2)
    LCONVMAX = max(LCONVMAX1,LCONVMAX2)


    !// Check for inconsistencies in input data; is the elevator starting
    !// with FLUX_E(L)>0  but where ENT_U(L-1)=0.? It should not happen, as
    !// this should be taken care of during metdata read-in.
    !// If is happens we set the starting mass of the elevator in layer L-1
    !// from FLUX_E(L).
    if (FLUX_E(LCONVMIN) .gt. 0._r8) then
       !// Treat the lowest level of convection
       L = LCONVMIN
       !// mass in from below is FLUX_E(L)
       if (L.eq.1) then
          write(6,'(a,es16.4)') f90file//':'//subr// &
               ': Flux inconsistencies at surface:',FLUX_E(L)
          stop
       end if
       !// Elevator air mass in the box below, to be lifted
       elev_mass_air(L-1) = FLUX_E(L)
       !// Vapor mass in the box below, to be lifted. Comes from entrainment
       !// below, but where flux (FLUX_E(L-1)) is zero. Entrained air will
       !// have specific humidity as ambient air in that layer.
       elev_mass_vapor(L-1) = FLUX_E(L)*ambient_q(L-1)
    end if
    !// Set initial temperature and potential temperature
    do L = 1, LWEPAR
       elev_t(L)     = ambient_t(L)
       elev_theta(L) = ambient_t(L)*(1000._r8/pres(L))**(Rd/cp)
    end do



    do L = LCONVMIN, LCONVMAX

       !// Outline;
       !//   1. Step 1:
       !//     a. Air rises and cools
       !//     b. Entrained air mixes with elevator air coming from below
       !//   2. Step 2:
       !//     a. Assume simultaneous mixing of 1.a and 1.b.
       !//        - New potential temperature, mass weighted
       !//        - Calculate condensation.
       !//   3. Treat detrainment of air, vapor and liquid water
       !//   4. Rain out

       !// Lifted air will decrease the temperature, whereas entrained air
       !// should increase the elevator temperature (if ambient air is warmer
       !// than the lifted air, i.e. more stable).
       !// After lifting and simultaneous entrainment, the condensed water
       !// can be calculated.

       !// At each level the dry air mass coming from below is given by the
       !// mass fluxes. For clarity we use fluxes rather than entrainment and
       !// detrainment to balance the fluxes.

       !// First locate the air masses relevant for the elevator
       !// If DET_U(L)<0, we have to add entrainment, so we include that:
       entrainment = ENT_U(L) - min(DET_U(L),0._r8)
       !// Detrainment
       detrainment = max(DET_U(L),0._r8)
       !// Air from below (equals elev_mass_air(L-1))
       air_up = FLUX_E(L)


       !// Total elevator air before mixing/condensation: lifted + entrained
       elev_mass_air(L) = air_up + entrainment

       !// The elevator mass may be zero; skip and go to next level
       if (elev_mass_air(L) .le. 0._r8) then
          if (elev_mass_air(L) .lt. 0._r8) then
             write(6,'(a,i5,es16.4)') f90file//':'//subr// &
                  ': NEGATIVE AIR IN ELEVATOR!',L,elev_mass_air(L)
             stop
          end if

          !// zero vapor and liquid water as well
          elev_mass_vapor(L) = 0._r8
          elev_mass_lw(L) = 0._r8

          !// Assume rain falls right through, and that net cprecip is zero
          !// (consistent with no elevator), even though this may not be the
          !// case (rain in may differ slightly from rain out, but differences
          !// should be small).

       else

          !// There is air in the plume

          !// Total vapor mass before condensation;
          !//   a. possibly coming from below
          !//   b. and from entrainment of ambient air
          !//   c. perhaps also evaporation of rain from above, but that would
          !//      change temperature:
          !//      evap_rain = max(0._r8, cprecip(L+1) - cprecip(L))
          if (L .gt. 1) then

             !// is there elevator air below at all?
             if (elev_mass_air(L-1) .gt. 0._r8) then
                !// Could use mass mixing ratio of vapor and find water
                !//   elev_q = elev_mass_vapor(L-1)/elev_mass_air(L-1)
                !//   vap_up = air_up * elev_q
                !// but this vapor should be equal to elev_mass_vapor(L-1)
                !// since air_up=FLUX_E(L)=elev_mass_air(L-1)
                vap_up = elev_mass_vapor(L-1)
             else
                !// No air or vapor coming from below
                vap_up = 0._r8
             end if

          else
             !// For layer 1 there is no air/vapor coming in.
             vap_up = 0._r8
          end if

          !// Vapor from below + entrained
          elev_mass_vapor(L) = vap_up + entrainment * ambient_q(L)


          !// Simple evaporation of rain
          !// Skip for now
          !liq_in = max(0._r8, cprecip(L+1) - cprecip(L))
          !evap_rain = 0._r8
          !if (liq_in .gt. 0._r8) then
          !   !// If so, evaporate, change temperature, liq_in, and find
          !   !// evap_rain
          !   call evaporate_rain(liq_in, elev_mass_air(L), &
          !     elev_mass_vapor(L), pres(L), L, elev_t(L), ambient_t(L+1), &
          !     evap_rain)
          !end if
       


          !// Adjust elevator temperature according to dry adiabatic lifting
          !// followed by entrainment. Apply mass weighting.
          !// Remember: elev_theta(L) has not been changed and is for ambient
          !// air.
          if (L .gt. 1) then
             !// Mass weighted theta for lifted air in elevator
             tmpv = elev_theta(L-1)*air_up
          else
             !// No contribution from flux in layer 1 (FLUX_E(1) is zero)
             tmpv = 0._r8
          end if
          !// Dividing by air_up is ok since tmpv is consistent with it.
          !// We have checked that the flux plus entrainment is positive
          elev_theta(L) = (tmpv + elev_theta(L)*entrainment) &
                           / ( air_up + entrainment )
          !// Calculate temperature from potential temperature
          elev_t(L) = elev_theta(L)*(pres(L)*1.e-3_r8)**(Rd/cp)

          !// These temperatures are the new elevator temperatures unless there
          !// is condensation, so then we check for condensation.
          !// Condense the lifted water and adjust elevator temperatures
          call condense_water(elev_mass_air(L), elev_mass_vapor(L), pres(L), &
               L, elev_t(L), elev_theta(L), cond_water)


          !// Liquid water in elevator
          !// 1. Coming from below; use mass mixing ratio of liquid water
          !//    below together with flux: elev_q_lw(L-1)*FLUX_E(L)
          !// 2. Condensed
          !// 3. Detrained (later)
          !// 4. Rained out (later)

          !// 1. Liquid water in updraft (coming from below)
          !//    Since FLUX_E(L) = elev_mass_air(L-1), the liquid water from
          !//    below is just elev_mass_lw(L-1). No need to find
          !//    "mixing ratio" below and multiply with air_up.
          if (air_up .gt. 0._r8) then
             elev_mass_lw(L) = elev_mass_lw(L-1)
          else
             elev_mass_lw(L) = 0._r8
          end if

          !// 2. Condensed
          elev_mass_lw(L) = elev_mass_lw(L) + cond_water

          !// 3. Detrainment
          !//    Elevator mass mixing ratio for vapor (specific humidity)
          !//    and liquid water, now given at current level
          elev_q  = elev_mass_vapor(L) / elev_mass_air(L)
          elev_mmr_lw = elev_mass_lw(L) / elev_mass_air(L)

          !//    Remove elevator dry air mass
          elev_mass_air(L) = max(elev_mass_air(L) - detrainment, 0._r8)

          !//    Remove elevator vapor; need elevator specific humidity
          elev_mass_vapor(L) = &
               max(elev_mass_vapor(L) - detrainment * elev_q, 0._r8)


          !//    Remove elevator liquid mass
          !//    This may end up to be unphysical when rain is added;
          !//    Rain will most likely not be detrained
          elev_mass_lw(L) = elev_mass_lw(L) - detrainment * elev_mmr_lw

          !//    The liquid mass will be used at the next lavel.

          !// 4. Rainout; this may be a little dubious.
          !//    If net rain is positive, more rains out at the bottom than
          !//    comes in at top.
          !//    If, however, the net is negative, more comes in at the top
          !//    than goes out at bottom. This is because of
          !//      a. There is production of rain, but possibly not matched
          !//         by b:
          !//      b. Evaproation of rain coming from above.
          !//      c. Rain is caught in the updraft, staying liquid.
          !//    When the net rain is positive, we assume all is due to
          !//    production, and remove it from the liquid water created.
          if (L .lt. LWEPAR) then
             rin_top = cprecip(L+1)
          else
             rin_top = 0._r8
          end if
          rainout = max(0._r8, cprecip(L) - rin_top)

          !// Set QFRAC and LW_VOLCONC.

          !// QFRAC is the ratio rainout / elev_mass_lw.
          !// For X kg of a water soluable component in the elevator, e.g.
          !// with 20% solubility, the amount rained out is
          !// X*0.2*rainout/elev_mass_lw(L) kg.
          if (elev_mass_lw(L) .gt. 1.e-4_r8 .and. rainout .gt. 0._r8) then
             !// Fraction of convective rain [kg/s] to elevator liquid
             !// water [kg/s]. Skip the smallest liquid water masses, to
             !// avoid Q=1 when as good as no water is present.
             !// One can argue that Q could be 1 even though there is almost
             !// no liquid water, but since the dissolved fraction will be
             !// small it does not matter much. Allow max 1.
             QFRAC(L,II,JJ,MP) = min( rainout / elev_mass_lw(L), 1._r8 )
          end if

          !// Do liquid water removal due to rain
          elev_mass_lw(L) = max(0._r8, elev_mass_lw(L) - rainout)
          !//    It is not easy to estimate the amount of evaporation; you
          !//    would have to build a more sophisticated cloud model to do
          !//    that. You could assume all would evaporate into the mix of
          !//    air coming up and entrained air, but that would lower the
          !//    temperature. Anyway, that part should be considered along
          !//    with entrainment at the top, or as a separate process after
          !//    lifting/entrainment.



          !// Use equation of state to calculate volume in the elevator;
          !// V = Rd T Ma / p
          elev_vol = Rd * elev_T(L) * elev_mass_air(L) / (100._r8 * pres(L))
          !// Only treat volumes greater than 0.1 m3 (arbitrary value to
          !// filter the smallest elevator volumes)
          if (elev_vol .gt. 1.e-1_r8) then
             !// The second key variable we need is volume concentration of
             !// liquid water in elevator. The volume of liquid water is
             !// elev_mass_lw/rho, where rho is the liquid water density
             !// (rho = 1.d3 kg/m3). 1/rho is better known as specific volume,
             !// and is given as the parameter alpha_LW.
             LW_VOLCONC(L,II,JJ,MP) = elev_mass_lw(L) * alpha_LW / elev_vol
          end if

          !// Skip LSEXTRA.
          !// If (rainout-elev_mass_lw(L)) > 0, then one could assume that
          !// some of the convective rainout is really large scale by keeping
          !// QFRAC=1 and adding (rainout-elev_mass_lw(L)) to large scale rain.
          !// HOWEVER; It is not very physical, since elev_mass_lw(L) already
          !// should be larger than rainout, since it also covers cloud
          !// droplets.
          !// Excess convective rain to be added to large scale rain
          !LSEXTRA(L,II,JJ,MP) = max(0._r8, rainout - elev_mass_lw(L))


       end if !// if (elev_mass_air(L) .le. 0._r8) then

    end do !// do L = LCONVMIN, LCONVMAX

    !// Save elevator temperature?
    !// Min temperature in elevator < 258K
    if (minval(elev_t(LCONVMIN:LCONVMAX)) .le. 258._r8) then
       LELEVTEMP(1,II,JJ,MP) = 1
    end if
    !// Max temperature in elevator < 273.15K
    if (maxval(elev_t(LCONVMIN:LCONVMAX)) .le. 273.15_r8) then
       LELEVTEMP(2,II,JJ,MP) = 1
    end if

    !// --------------------------------------------------------------------
  end subroutine elevator_fractions
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine condense_water(emass_d, emass_v, pres, L, &
       elev_t, elev_the, cond_water)
    !// --------------------------------------------------------------------
    !// Calculate liquid water condensed after lifted and entrained
    !// air is mixed.
    !// Modifies elevator temperature and potential temperature.
    !//
    !// Taken from iterative loop in old CTM2 code.
    !//
    !// Ole Amund Sovde, November 2009
    !// --------------------------------------------------------------------
    use cmn_parameters, only: R_AIR, CP_AIR, Lv_0C
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    integer, intent(in)   :: L
    real(r8), intent(in)    :: emass_d, & !// Elevator dry mass 
                             emass_v, & !// Elevator liquid mass
                             pres       !// Pressure [hPa]
    real(r8), intent(inout) :: elev_the, elev_t
    real(r8), intent(out)   :: cond_water


    real(r8) , parameter :: &
          TLIM = 0.01_r8, &     !// Minimum dT allowed due to condensation
          Rd = R_AIR, &         !// Gas constant for dry air.
          EPSILON = 0.622_r8, & !// Epsilon = Rd/Rv
          CP = CP_AIR, &        !// Specific heat capacity for air at
                                !// constant pressure.
          Lw = Lv_0C, &         !// Latent heat for water.
          alpha_LW = 1.e-3_r8, &!// Specific volume (inverse density) of
                                !// liquid water [m3/kg]
          firstguess = 0.2_r8, &!// Assume 20% condensation as a first guess.
          secondguess = 0.2_r8  !// Assume 20% more condensation in each
                                !// iteration step.

    real(r8) :: &
          es, esn, e_h2o, e_gas, de, dex, &!// Vapor pressure variables
          TH, Te, tx, dT, dtt, T1          !// Temperature variables

    integer :: IT
    integer , parameter :: MAXITER = 100
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'condense_water'
    !//---------------------------------------------------------------------


    !// Partial pressure of water vapor in the elevator after lifting and
    !// entrainment, but before possible condensation
    e_h2o = (pres/epsilon) * emass_v / emass_d

    !// Temperature of elevator after dry adiabatic lifting and entrainment
    !// (without condensation) has already been calculated.
    Te = elev_t

    !// Saturation partial pressure of water after lifting.
    !// Rogers and Yau, A short course in Cloud Physics,
    !// Butterworth & Heinemann, 
    !// International series in natural philosophy, vol 113, 1989 (3rd edition),
    !// page 16.
    tx = Te - 273.15_r8 !// Celcius
    es = 6.112_r8 * exp(17.67_r8 * Tx / (Tx + 243.5_r8))

    If (e_h2o .gt. es) Then
       !// Supersaturation; assume condensation

       !// Will save original vapor pressure, use e_gas in iteration
       e_gas = e_h2o

       !// Max possible condensed water
       de = e_h2o - es

       !// Not everything should condensate; first try with the fraction
       !// firstguess
       dex = firstguess * de

       !// Increase in T of the elevator due to release of latent heat
       dT = Lw * dex * epsilon / (pres * cp)

       !// Iterative procedure to calculate the true temperature of the
       !// elevator, including the release of latent heat.
       !// Will iterate, and stop iteration when DTT < TLIM

       dtt = dt
       IT = 0

       do while (abs(DTT) .gt. TLIM)

          !// Count number of iterations
          IT = IT + 1

          !// Remaining water in gas phase after condensation
          e_gas = e_gas - dex

          !// New temperature after condensation (Celcius)
          T1 = Te + dT - 273.15_r8

          !// Saturation vapor pressure after condensation
          esn = 6.112_r8 * exp(17.67_r8 * T1 / (T1 + 243.5_r8))

          !// Check if still supersaturation
          de = e_gas - esn
          if (de .gt. 0._r8) then
             !// Not all should evaporate, will try with secondguess
             dex = secondguess * de
             !// Increase in T of the elevator due to release of latent heat
             dTT = Lw*dex*epsilon/(pres*cp)
             !// Add the additional increase to dT
             dT = dT + dTT
          else
             !// No longer supersaturation; exit
             !// Do not change dT anymore.
             !// Note that dT may be large enough to give de<0, which is not
             !// a problem.
             exit
          End If

          if (IT .gt. MAXITER) then
             write(6,'(a)') f90file//':'//subr// &
                  ': Too many iterations!'
             write(6,*) 'IT >  ',MAXITER
             write(6,*) 'esn  : ',esn
             write(6,*) 'e_gas: ',e_gas
             write(6,*) 'dTT  : ',dTT
             write(6,*) 'dT   : ',dT
             stop
          End If
       end do  !// do while (abs(DTT) .gt. TLIM)


       !// Have found temperature of elevator with mixture of vapor and
       !// droplets.
       !// Update temperature and potential temperature of the elevator:
       Te = Te + dT
       !// Adjust the elevator temperatures
       elev_the  = Te * (1000._r8/pres)**(Rd/cp)
       elev_t    = Te


       !// Calculate liquid water mass in the elevator [kg/s].
       !//    Get saturation pressure:
       Tx = Te - 273.15_r8
       es = 6.112_r8 * exp(17.67_r8 * Tx / (Tx + 243.5_r8))
       !//    Condensed water formed due to temperature drop (partial pressure):
       de = e_h2o - es
       cond_water = emass_v * de/e_h2o

       if (cond_water .lt. 0._r8) then
          write(6,'(a)') f90file//':'//subr//': Condensed water negative!'
          print*,Te,e_h2o,es, de
          stop
       end if

       if (abs(de/e_h2o) .gt. 1.) then
          write(6,'(a)') f90file//':'//subr//': Wrong vapor!'
          print*,'de',de
          print*,'e_h2o'
          stop
       end if

    else
       !// No condensation
       cond_water = 0._r8
       !// No changes in temperatures (input values were after dry adiabatic
       !// lifting and then entrainment.
    end if !// Check on condensatiion



    !// --------------------------------------------------------------------
  end subroutine condense_water
  !// ----------------------------------------------------------------------



  !// ----------------------------------------------------------------------
  subroutine evaporate_rain(liq_in, emass_d, emass_v, pres, L, &
       elev_t, rain_t, evap_rain)
    !// --------------------------------------------------------------------
    !// Calculate evaporation of rain from above.
    !// Modifies elevator temperature, potential temperature and incoming rain.
    !// --------------------------------------------------------------------
    use cmn_parameters, only: R_AIR, CP_AIR, Lv_0C
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------

    integer, intent(in) :: L
    real(r8), intent(in) :: emass_d, emass_v, pres, rain_t
    real(r8), intent(inout) :: elev_t, liq_in
    real(r8), intent(out) :: evap_rain


    real(r8) , parameter :: &
         min_dt = -5._r8, &    !// Max absolute change allowed in temperature
         max_rh = 1._r8, &     !// Max RH due to evap. and temperature change
         Rd = R_AIR, &         !// Gas constant for dry air.
         EPSILON = 0.622_r8, & !// Epsilon = Rd/Rv
         cp_d = CP_AIR, &      !// Specific heat capacity (cp) for air at
                               !// constant pressure.
         cp_v = 1859._r8, &    !// cp for vapor at constant pressure.
         cp_lw = 4186._r8, &   !// cp for liquid water at constant pressure.
         Lw = Lv_0C            !// Latent heat for water.

    real(r8) :: e_1, e_s, e_max, &
         Tx, T1, T2, Tn, dT, dT1, dT2

    logical, parameter :: debug=.true.
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'evaporate_rain'
    !//---------------------------------------------------------------------



    !// Properties of rain from above
    !//    mass of rain [kg/s]: mrw
    !//    temperature:         T2
    !// Dry mass in elevator
    !//    mass:        md
    !//    temperature: T1
    !// Vapor mass in elevator
    !//    mass:        mv
    !//    temperature: T1
    !// New temperature after evaporation: Tn

    !// Energy required to change temperature of rain from T2 to T1
    !// (T2<T1: removing heat from surroundings, T2>T1: adding heat):
    !//   dErw = mrw * cp_lw * (T2-T1)    (positive if temperature increases)
    !// Energy to evaporate rain (removing heat from surroundings)
    !//   dEvap = mrw * Lw                (positive; removing heat)
    !// Energy to change dry mass from T1 to Tn
    !//   dEd = md * cp_d * (Tn-T1)       (negative for temperature decrease)
    !// Energy to change vapor mass
    !//   dEv1 = mv * cp_v * (Tn-T1)      (negative for temperature decrease)
    !//   dEv2 = mrw * cp_v * (Tn-T1)     (negative for temperature decrease)
    !// The total change in energy is assumed to be unchanged
    !//   dErw + dEvap + dEd + dEv1 + dEv2 = 0

    !// Catches:
    !// 1. Allow max dT = -5K
    !// 2. Allow max X% saturation after evaporation

    !// Current vapor pressure
    e_1 = (pres/epsilon) * emass_v / emass_d

    !// Current saturation vapor pressure
    Tx = elev_t - 273.15_r8
    e_s = 6.112_r8 * exp(17.67_r8 * Tx / (Tx + 243.5_r8))

    !// If all rain evaporates, the partial pressure of vapor is
    e_max = (pres/epsilon) * (emass_v + liq_in)/emass_d

    !// Change air temperature due to increase/decrease of rain temperature
    !// i.e. total dE without dEvap
    !//   mrw*cp_lw*(Tn-T2) + md*cp_d*(Tn-T1) + mv*cp_v*(Tn-T1) = 0
    !//   Tn = mrw*cp_lw*T2 + md*cp_d*T1 + mv*cp_v*T1
    !//        --------------------------------------
    !//            mrw*cp_lw + md*cp_d + mv*cp_v
    T1 = elev_t
    T2 = rain_t

    Tn = (liq_in * cp_lw*T2 + emass_d*cp_d*T1 + emass_v*cp_v*T1) &
          / (liq_in * cp_lw + emass_d*cp_d + emass_v*cp_v)
    dT1 = Tn - T1
    if(dT1 .ne. dT1) then
       write(6,'(a,6es16.4)') f90file//':'//subr// &
            ': dT error', dt1, t1, tn, liq_in,emass_d,emass_v
       stop
    end if
    !// abs(dT) can not be unrealisticly large. Allow dT=max(min_dt,dT).
    !// Evaporation will lower temperature, so for now it is ok with dT>min_dt.
    if (dT1 .lt. min_dt) then
       !// Find Tn to match dT limit
       Tn = T1 + min_dt
       !// Skip evaporation; too cold for evaporation
       evap_rain = 0._r8
       !// Rain is kept as liquid water
    else
       !// Have new elevator temperature T1
       T1 = Tn
       !// Assume all evaporates
       !//   mrw*cp_v*(Tn-T1) + md*cp_d*(Tn-T1) + mv*cp_v*(Tn-T1) + mrw*Lw = 0
       dT2 = -liq_in * Lw / (emass_d*cp_d + (emass_v + liq_in)*cp_v)
       !// Test for dT; perhaps adjust evaporation
       !// Use total change in dT = dT1+dT2
       dT = dT1 + dT2
       if(dT .ne. dT) then
          write(6,'(a,3es16.4)') f90file//':'//subr// &
               ': dT error2', dt1, t1, tn
          stop
       end if

       if (dT .lt. min_dt) then
          !// Too much evaporated, limit by min_dt
          dT = min_dt
          !// m*cp_v*dT + m*Lw = -md*cp_d*dT - mv*cp_v*dT
          evap_rain = -(emass_d*cp_d + emass_v*cp_v)*dT / (cp_v*dT + Lw)
          !// do not allow negative evaporation
          evap_rain = max(evap_rain, 0._r8)

       else if (dT .gt. -min_dt) then
          !// The change in T is too large! T2 must be much larger than T1?
          !// I don't think this should happen, so we stop for now
          !// More mass of rain than air, perhaps?
          write(6,'(a,es16.4)') f90file//':'//subr// &
               ': dT from change of rain_t',dT1
          write(6,*) 'dT due to evaporation   ',dT2
          write(6,*) 'elevator temperature',elev_t
          write(6,*) 'rain temperature    ',rain_t
          write(6,*) 'elevator mass',emass_d
          write(6,*) 'rain mass    ',liq_in
          stop
       else
          !// All is evaporated within specified dT criteria
          evap_rain = liq_in
          !// And no additional liquid water to the elevator
          liq_in = 0._r8
       end if

       !// Check for relative humidity
       !// New temperature in elevator after changing rain temperature and evap.
       Tn = elev_t + dT
       Tx = Tn - 273.15_r8
       e_s = 6.112_r8 * exp(17.67_r8 * Tx / (Tx + 243.5_r8))
       !// Current vapor pressure
       e_1 = (pres/epsilon) * (emass_v + evap_rain) / emass_d
       !// Check if condensation occurs... Do not allow that yet
       if (e_1/e_s .gt. max_rh) then
          !// Evaproated rain to match max_rh
          evap_rain = max_rh * emass_d * epsilon / pres - emass_v
          !// Did that change the temperature? Assume first order increase
          !// and leave.
          !// m*cp_v*dT + m*Lw = -md*cp_d*dT - mv*cp_v*dT
          !// dT = -m*Lw / (m*cp_v + md*cp_d + mv*cp_v)
          dT2 = -evap_rain * Lw / (emass_d*cp_d + (emass_v + evap_rain)*cp_v)
          !// Re-adjust dT and assume evaporation/condensation is ok
          dT = dT1 + dT2
          if (debug) then
            print*,'Should have condensation after evap of rain'
            print*,elev_t,dT,e_s,e_1,e_1/e_s
          end if
       end if

       !// Adjust incoming rain
       liq_in = liq_in - evap_rain

    end if

    !// Adjust elevator temperature
    elev_t = Tn

    !// --------------------------------------------------------------------
  end subroutine evaporate_rain
  !// ----------------------------------------------------------------------


  !// ----------------------------------------------------------------------
  subroutine liquid_fractions(I,J,II,JJ,MP,BTEM, CNV_WETL)
    !// --------------------------------------------------------------------
    !// Calculate mass fraction of tracer removed by convective precipitation
    !// (i.e. the fraction of the tracer solved in the elevator plume).
    !// This is done in two steps for each tracer:
    !//
    !// 1. Calculate mass fraction of the tracer dissolved in cloud droplets,
    !//    in the elevator (not whole grid box) namely the DISSOLVEDFRAC.
    !// 2. Only a fraction of the droplets will be removed by precipitation,
    !//    namely QFRAC.
    !//
    !// Units
    !//   QFRAC:      Fraction conv. rain/Cloudwater in elevator
    !//   LW_VOLCONC: Liquid water (kg) / elevator volume (m3)
    !//
    !// The mass fraction of tracer which is washed out is therefore
    !//   CNV_WETL = DISSOLVEDFRAC * QFRAC
    !//
    !// Calculations done in the column, and are only stored for the column.
    !//
    !// Ole Amund Sovde, April-November 2009
    !// --------------------------------------------------------------------
    use cmn_size, only: NPAR
    use cmn_chem, only: TCHENA, TCHENB, TCKAQA, TCWETL, TNAME
    use cmn_parameters, only: R_UNIV
    use cmn_oslo, only: chem_idx, TCCNVHENRY
    !// --------------------------------------------------------------------
    Implicit NONE
    !// --------------------------------------------------------------------
    !// Input
    integer , intent(in) :: I, J, &   !// Global indices of the column
                            II, JJ, & !// IJ-block indices of the column
                            MP        !// IJ-block index

    !// B-array input
    real(r8) , intent(in) :: BTEM(LPAR,IDBLK,JDBLK)   !// Temperature [K]

    !// Output
    real(r8),intent(out) :: CNV_WETL(LWEPAR,NPAR) !// Fraction of tracer in cloud 

    !// Locals
    integer :: comp_idx, &          !// Counter for chemical components
               N, L                 !// Counter for transport number and layers
    real(r8) :: DISSOLVEDFRAC(LWEPAR)   !// Fraction dissolved in cloudwater
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'liquid_fractions'
    !//---------------------------------------------------------------------


    !// Initialize
    CNV_WETL(:,:) = 0._r8

    !// Go through all transported species
    !// Only transported species can be washed out!
    do N = 1, NPAR

       !// Component number
       comp_idx = chem_idx(N)

       !// Check for convective removal
       if (TCCNVHENRY(N) .eq. 0 .or. TCWETL(N) .eq. 0._r8) cycle

       !// TCWETL: Fraction of gridbox allowed to be dissolved
       !// TCCNVHENRY: Determines how to calculate the fraction removed.

       !// If TCWETL=0 or TCCNVHENRY=0, there will be no convective washout!

       if (TCCNVHENRY(N) .eq. 1) then
          !// Use Henry coefficients, possibly modified/hard coded
          call wf_henry(LPAR,LWEPAR,chem_idx(N),BTEM(:,II,JJ), &
               LW_VOLCONC(:,II,JJ,MP),TCHENA(N),TCHENB(N),TCKAQA(N), &
               R_UNIV,DISSOLVEDFRAC)

          !// Set convective wetloss
          CNV_WETL(:,N) = TCWETL(N) * DISSOLVEDFRAC(:) * QFRAC(:,II,JJ,MP)

       else if (TCCNVHENRY(N) .eq. 2) then
          write(6,'(a)') f90file//':'//subr// &
               ': TCCNVHENRY(N) .eq. 2 is not in use!'
          stop

       else if (TCCNVHENRY(N) .eq. 3) then
          !// Use the constant WETL parameter, but multiply with QFRAC
          CNV_WETL(:,N) = TCWETL(N) * QFRAC(:,II,JJ,MP)

       else if (TCCNVHENRY(N) .eq. 4) then
          !// Use the constant WETL parameter as in UCI
          !// It should not be multiplied with QFRAC unless it is defined
          !// otherwise.
          CNV_WETL(:,N) = TCWETL(N)

       else if (TCCNVHENRY(N) .eq. 5) then
          !// Same as 3, but no large scale scavenging!
          CNV_WETL(:,N) = TCWETL(N) * QFRAC(:,II,JJ,MP)

       else if (TCCNVHENRY(N) .eq. 6) then
          !// Same as 3, but also removes tracer due to large scale ice
          CNV_WETL(:,N) = TCWETL(N) * QFRAC(:,II,JJ,MP)

       else if (TCCNVHENRY(N) .eq. 7) then
          !// Convective removal when Tplume<258K (i.e. on ice)
          !// This option also removes tracer due to large scale ice
          if (LELEVTEMP(1,II,JJ,MP) .eq. 1) then
             CNV_WETL(:,N) = TCWETL(N) * QFRAC(:,II,JJ,MP)
          else
             CNV_WETL(:,N) = 0._r8
          end if

       else if (TCCNVHENRY(N) .eq. 8) then
          !// Convective removal when min(Tplume)<258K (i.e. assume ice)
          !// and max(Tplume)<273K (i.e. assume all is ice).
          !//   Assume 20% reduction in this case.
          !// This option also removes tracer due to large scale ice
          if (LELEVTEMP(1,II,JJ,MP) .eq. 1 .and. LELEVTEMP(2,II,JJ,MP) .eq. 1) then
             CNV_WETL(:,N) = TCWETL(N) * QFRAC(:,II,JJ,MP) * 0.2_r8
          else
             CNV_WETL(:,N) = 0._r8
          end if
       else
          write(6,'(a,es16.4)') f90file//':'//subr// &
               ': CANNOT TREAT TCCNVHENRY=',TCCNVHENRY(N)
          stop
       end if

    end do !// do N = 1, NPAR


    !// --------------------------------------------------------------------
  end subroutine liquid_fractions
  !// ----------------------------------------------------------------------





    !// --------------------------------------------------------------------
  subroutine wf_henry(LPAR,LWEPAR,chem_idx,TEM,LW_VOLCONC, &
                                 KHA,KHB,HFLAG,R_UNIV,DISFRAC)
    !// --------------------------------------------------------------------
    !// Calculate fraction of mass washed out, using regular Henry constants.
    !// HFLAG denotes whether hard-coded effective constants should be
    !// calculated.
    !//
    !// DISFRAC is the fraction of tracer mass dissolved in cloud water
    !// in the grid cell; thus the fraction removed if cloud water disappears
    !// through rainout.
    !//
    !// Henry's law coefficient for any gas is defined as
    !//    P(gas) = k_H * X                                             (1)
    !// where P(gas) is the partial pressure of the gas above the solution,
    !// and X is the molar fraction of the dissolved gas in the solution:
    !//    X = n_aq / (n_aq + n_solvent)
    !// Re-writing to concentration (dividing by volume of solution) we have
    !//    X = C_aq / (C_aq + C_solvent)
    !// and since we always have C_aq << C_solvent, this can be approximated
    !//    X = C_aq / C_solvent                                         (2)
    !//
    !// Since the concentration of the solvent (water) is approximately
    !// constant, we arrive at the other common form of Henry's law
    !//    P(gas) = k * C_aq                                            (3)
    !//
    !// Units: k [atm*L(solvent)/mol], P(gas) [atm]
    !//
    !// If k is high, it means the component prefers thermodynamically
    !// to be dissolved in the liquid phase.
    !//
    !// We are interested in calculating C_aq from P(gas), so we introduce
    !//   HSTAR = 1/k  [mol/(atm*L(solvent))]                           (4)
    !// so that
    !//   C_aq = HSTAR * P(gas)
    !//
    !// If we want to apply the calculations to molar concentration (mol/L),
    !// so we have to change some units:
    !//    P(gas) = C_g * R * T                                         (5)
    !// given correct units of R, i.e. [atm*L /(mol*K):
    !//    J = kgm^2/s^2 = Pa * m^3
    !//    R = 8.31451 [J/(mol*K) / 101325 Pa/atm * 1000 L/m3]
    !//      = 0.0820578 [atm*L / (mol*K)]
    !//
    !// Henry's law therefore implies that the concentration in the solution
    !// is proportional to the atmospheric concentration: 
    !//   C_aq = HSTAR * R * T * C_g = HH * C_g                         (6)
    !// where HH has the units [mol/L(solvent) / (mol/L(air))].
    !//
    !// We want the mass fraction of the dissolved gas, which equals the
    !// molar fraction.
    !//   fraction = n_aq / (n_aq + n_g)                                (7)
    !//
    !// Number of moles in air (C_g has units [mol/L(air)])
    !//   n_g  = C_g * V_elevator_air * 10^3L(air)/m3(air)              (8)
    !//
    !// Number of moles in solution (C_aq has units [mol/L(solvent)])
    !//   n_aq = C_aq * V_elevator_solvent * 10^3L(solvent)/m3(solvent) (9)
    !//
    !// If you are confused about "solvent": The solvent is liquid water.
    !// The volume of the solvent is given by the liquid water content, and
    !// is calculated in elevator_fractions as volume concentration, i.e.
    !// volume of liquid water to total volume. I repeat the calculation
    !// here; the volume of liquid water is
    !//   V_elevator_solvent = liquid water / rho
    !// where rho is the density of water (which is 10^3kg/m3). The volume
    !// concentration is then
    !//   LW_VOLCONC = V_elevator_solvent / V_elevator                 (10)
    !//
    !// The moles of gas dissolved in water is therefore
    !//   n_aq = C_aq * LW_VOLCONC * V_elevator * 10^3L(solvent)/m3(solvent)
    !//
    !// The mass fraction (this fraction equals mole fraction) solved in the
    !// droplets, which are subject to washout, are therefore
    !//   DISFRAC = n_aq / (n_aq + n_g)
    !//           = C_aq*LW_VOLCONC / (C_aq*LW_VOLCONC + C_g)          (11)
    !//
    !// and by equation (6) we get
    !//   DISFRAC = HH*LW_VOLCONC / (HH*LW_VOLCONC + 1)                (12)
    !//
    !// Remark: The change from L to m3 cancels, even though one is for air
    !// and the other for solvent.
    !//
    !// Ole Amund Sovde, April 2016, April 2009
    !// --------------------------------------------------------------------
    use scavenging_largescale_uci, only: getHstar
    !// --------------------------------------------------------------------
    implicit none
    !// --------------------------------------------------------------------
    !// Input
    integer, intent(in)  :: LPAR,LWEPAR    !// Size of column
    integer, intent(in)  :: chem_idx       !// Tracer ID
    real(r8), intent(in) :: TEM(LPAR), &
                          LW_VOLCONC(LWEPAR), &
                          KHA,KHB,HFLAG, &
                          R_UNIV

    !// Output
    real(r8), intent(out) :: DISFRAC(LWEPAR)

    !// Locals
    real(r8), parameter  :: INV298 = 1._r8/298._r8
    real(r8) :: HSTAR, LW_VC, HH, RGAS_MOD
    integer :: L
    integer, parameter :: HCALLER = 2 !// caller flag for getHstar
    !// --------------------------------------------------------------------
    character(len=*), parameter :: subr = 'wf_henry'
    !//---------------------------------------------------------------------

    !// Initialize
    DISFRAC(:) = 0._r8

    !// Modified universal gas constant, scaled from [J/(mole*K)]
    !// to [atm*L/(mol*K)]
    RGAS_MOD = R_UNIV / 1.01325e5_r8 * 1.e3_r8

    !// Can take ice temperature into account, as for large scale.
    do L = 1, LWEPAR
       LW_VC = LW_VOLCONC(L)
       !// The regular Henry expression [mol/(atm*L)]
       !HSTAR = KAQ * KHA * exp(KHB/HTEMP - KHB*INV298)
       call getHstar(chem_idx, KHA, KHB,HFLAG, TEM(L), HCALLER, HSTAR)
       !// Modify HSTAR to be used with concentration instead of pressure,
       !// i.e. C_aq = HSTAR*R*T * C_g. This modified expression HH is
       !// unit-less.
       HH = HSTAR * RGAS_MOD * TEM(L)

       !// Mass fraction dissolved in liquid water
       DISFRAC(L) = HH * LW_VC / (HH * LW_VC + 1._r8)
    end do

    !// --------------------------------------------------------------------
  end subroutine wf_henry
  !// ----------------------------------------------------------------------




  !// ----------------------------------------------------------------------
end module cnv_oslo
!//=========================================================================

! $Header: /mn/hox/ozon/ctm2_model/ctm2_src/dst_dst/dstlsm.F90,v 1.1 2003/04/15 14:41:35 alfgr Exp $ -*-f90-*-

! Purpose: LSM information needed by dust parameterization
! dstlsm is initialized the first time it is used
! dstlsm is first used in BXM:aer(), CCM:control/ccm3(), MATCH:src/main()
! dstlsm most emphatically does NOT need params.h or lsmpar.h to work

! First half of module is Fortran 90 translation of CCM:lsm/vegtyp.h
! Second half of module is Fortran 90 translation of CCM:lsm/vegtypi.F

! ------------------------ code history ---------------------------
! source file:       vegtypi.F
! purpose:           subgrid plant type and fractional area for surface types
! date last revised: March 1996 - lsm version 1
! author:            Gordon Bonan
! standardized:      J. Truesdale, Feb. 1996
! reviewed:          G. Bonan, Feb. 1996
! -----------------------------------------------------------------

! Usage:
! use dstlsm ! [mdl] LSM data

module dstlsm ! [mdl] LSM data
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  
  ! Initialize fractional plant cover corresponding to LSM surface types
  ! dstlsm is currently initialized in BXM:aer(), CCM:control/ccm3(), MATCH:src/main()
  
  ! Contents of lsmpar.h
  integer mband             !number of solar radiation bands: vis, nir
  integer msl               !number of soil layers
  integer mst               !number of "soil" types (soil, ice, 2 lakes, wetland)
  integer mvt               !number of plant types
  integer msc               !number of soil color types
  
  parameter (mband=2, msl=6, mst=5, mvt=14, msc=9)
  
  ! Contents of vegtyp.h with cover <-> pln_frc and plant <-> pln_typ
  ! and vegcon_i <-> dstvegcon_i and vegcon_r <-> dstvegcon_r
  ! dstvegtyp_i and dstvegtyp_r
  real(r8),public::pln_frc(0:28,3) ! [frc] Weight of corresponding plant type (sums to 1.0)
  integer,public::pln_typ(0:28,3) ! [idx] LSM plant type (1..14 = nbr_LSM_pln_typ)
  
  ! ------------------------ code history ---------------------------
  ! source file:       vegcon.h
  ! purpose:           vegetation type constants 
  ! date last revised: March 1996 - lsm version 1
  ! author:            Gordon Bonan
  ! standardized:      J. Truesdale, Feb 1996
  ! reviewed:          G. Bonan, Feb 1996
  ! -----------------------------------------------------------------
  
  ! dstvegcon_i
  integer nic               !value for irrigated crop 
  integer noveg             !value for not vegetated 
  
  ! dstvegcon_r
  real(r8) vw(mvt)          !btran exponent: [(h2osoi-watdry)/(watopt-watdry)]**vw 
  real(r8) rdp(mvt)         !defines root fraction decrease with depth
  real(r8) ch2op(mvt)       !maximum intercepted h2o per unit lai+sai (mm) 
  real(r8) dleaf(mvt)       !characteristic leaf dimension (m) 
  
  ! dstvegcon_r
  real(r8) c3psn(mvt)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8) kc25(mvt)        !co2 michaelis-menten constant at 25c (pa) 
  real(r8) akc(mvt)         !q10 for kc25 
  real(r8) ko25(mvt)        !o2 michaelis-menten constant at 25c (pa) 
  real(r8) ako(mvt)         !q10 for ko25 
  real(r8) vcmx25(mvt)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
  real(r8) avcmx(mvt)       !q10 for vcmx25 
  real(r8) bp(mvt)          !minimum leaf conductance (umol/m**2/s) 
  real(r8) mp(mvt)          !slope of conductance-to-photosynthesis relationship 
  real(r8) qe25(mvt)        !quantum efficiency at 25c (umol co2 / umol photon)
  real(r8) aqe(mvt)         !q10 for qe25 
  real(r8) rmf25(mvt)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
  real(r8) rms25(mvt)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
  real(r8) rmr25(mvt)       !root maintenance respiration at 25c (umol co2/kg bio/s)
  real(r8) arm(mvt)         !q10 for maintenance respiration
  real(r8) dmcf(mvt)        !co2-to-biomass conversion factor (ug biomass/umol co2) 
  real(r8) folnmx(mvt)      !foliage nitrogen concentration when f(n)=1 (%) 
  real(r8) tmin(mvt)        !minimum temperature for photosynthesis (kelvin) 
  
  ! dstvegcon_r
  real(r8) xl(mvt)          !leaf/stem orientation index
  real(r8) rhol(mvt,mband)  !leaf reflectance: 1=vis, 2=nir 
  real(r8) rhos(mvt,mband)  !stem reflectance: 1=vis, 2=nir 
  real(r8) taul(mvt,mband)  !leaf transmittance: 1=vis, 2=nir 
  real(r8) taus(mvt,mband)  !stem transmittance: 1=vis, 2=nir 
  
  ! dstvegcon_r
  real(r8) binvvt(mvt)      !1/vkc*ln(z0m/z0h)
  real(r8),public::z0mvt(mvt)       !momentum roughness length (m) 
  real(r8),public::zpdvt(mvt)       !displacement height (m) 
  real(r8) stembvt(mvt)     !stem biomass (kg /m**2)
  real(r8) rootbvt(mvt)     !root biomass (kg /m**2)
  real(r8) folnvt(mvt)      !foliage nitrogen concentration (%)
  real(r8) mrp(mvt)         !microbial respiration parameter (umol co2 /kg c/ s)
  real(r8) soilcvt(mvt)     !soil carbon (kg c /m**2)
  real(r8) cwpvt(mvt)       !empirical canopy wind parameter
  
  ! dstvegcon_r
  real(r8),public::tai(mvt,12)      !monthly leaf area index + stem area index, one-sided
  real(r8),public::gai(mvt,12)      !monthly leaf area index, one-sided
  real(r8) hvt(mvt)         !top of canopy (m)
  real(r8) hvb(mvt)         !bottom of canopy (m)
  
  ! ------------------------ end vegcon.h ---------------------------
  
  !++csz
  !#include <vartyp.h>
  !--csz
  
  ! ------------------------ code history ---------------------------
  ! source file:       vegtypi.F
  ! purpose:           subgrid plant type and fractional area for surface types
  ! date last revised: March 1996 - lsm version 1
  ! author:            Gordon Bonan
  ! standardized:      J. Truesdale, Feb. 1996
  ! reviewed:          G. Bonan, Feb. 1996
  ! -----------------------------------------------------------------
  
  !++csz
  !#include <vegtyp.h>
  !--csz
  
  integer i          ! loop index
  
  ! there are 29 land surface types: 0 = ocean, 1 to 28 = land. each
  ! land point has up to three vegetation types, ranging in value from
  ! 1 to 14. [plant] contains the vegetation type of the 3 subgrid points 
  ! for each surface type. [cover] contains the fractional area of the 3 
  ! subgrid points for each surface type.
  
  data (pln_typ(i,1),i=0,28) /   0, &
       &                            14,  14,   1,   2,   4,   1  , 1, &
       &                             4,   1,   3,   5,  13,   1,   2, &
       &                            11,  11,   6,  13,   9,   7,   8, &
       &                             8,  12,  11,  12,  11,   3,  14/
  data (pln_frc(i,1),i=0,28) /0.00, &
       &                          1.00,1.00,0.75,0.50,0.75,0.37,0.75, &
       &                          0.75,0.37,0.95,0.75,0.70,0.25,0.25, &
       &                          0.40,0.40,0.60,0.60,0.30,0.80,0.80, &
       &                          0.10,0.85,0.85,0.85,0.85,0.80,1.00/
  data (pln_typ(i,2),i=0,28) /   0, &
       &                            14,  14,  14,  14,  14,   4  ,14, &
       &                            14,   4,  14,  14,   5,  10,  10, &
       &                             4,   4,  13,   6,  10,  14,  14, &
       &                            14,  14,  14,  14,  14,  14,  14/
  data (pln_frc(i,2),i=0,28) /0.00, &
       &                          0.00,0.00,0.25,0.50,0.25,0.37,0.25, &
       &                          0.25,0.37,0.05,0.25,0.30,0.25,0.25, &
       &                          0.30,0.30,0.20,0.20,0.30,0.20,0.20, &
       &                          0.90,0.15,0.15,0.15,0.15,0.20,0.00/
  data (pln_typ(i,3),i=0,28) /   0, &
       &                            14,  14,  14,  14,  14,  14,  14, &
       &                            14,  14,  14,  14,  14,  14,  14, &
       &                             1,   1,  14,  14,  14,  14,  14, &
       &                            14,  14,  14,  14,  14,  14,  14/
  data (pln_frc(i,3),i=0,28) /0.00, &
       &                          0.00,0.00,0.00,0.00,0.00,0.26,0.00, &
       &                          0.00,0.26,0.00,0.00,0.00,0.50,0.50, &
       &                          0.30,0.30,0.20,0.20,0.40,0.00,0.00, &
       &                          0.00,0.00,0.00,0.00,0.00,0.00,0.00/
  
  ! ------------------------------------------------------------------
  ! description of the 29 surface types 
  ! ------------------------------------------------------------------
  
  ! no vegetation
  ! -------------
  !  0 ocean                                 
  !  1 land ice (glacier)                             
  !  2 desert                                
  
  ! forest vegetation
  ! -----------------
  !  3 cool needleleaf evergreen tree           
  !  4 cool needleleaf deciduous tree           
  !  5 cool broadleaf  deciduous tree           
  !  6 cool mixed needleleaf evergreen and broadleaf deciduous tree    
  !  7 warm needleleaf evergreen tree           
  !  8 warm broadleaf  deciduous tree            
  !  9 warm mixed needleleaf evergreen and broadleaf deciduous tree    
  ! 10 tropical broadleaf evergreen tree  
  ! 11 tropical seasonal deciduous tree         
  
  ! interrupted woods
  ! ----------------
  ! 12 savanna                               
  ! 13 evergreen forest tundra               
  ! 14 deciduous forest tundra               
  ! 15 cool forest crop                           
  ! 16 warm forest crop                           
  
  ! non-woods
  ! ---------
  ! 17 cool grassland                             
  ! 18 warm grassland                            
  ! 19 tundra                              
  ! 20 evergreen shrub                   
  ! 21 deciduous shrub                 
  ! 22 semi-desert                     
  ! 23 cool irrigated crop                
  ! 24 cool non-irrigated crop                
  ! 25 warm irrigated crop               
  ! 26 warm non-irrigated crop               
  
  ! wetlands
  ! --------
  ! 27 forest (mangrove)                    
  ! 28 non-forest                          
  ! ------------------------------------------------------------------
  
  ! ------------------------------------------------------------------
  ! description of the 14 plant types. see vegconi.F for parameters
  ! that depend on vegetation type
  ! ------------------------------------------------------------------
  
  !  1 = needleleaf evergreen tree
  !  2 = needleleaf deciduous tree
  !  3 = broadleaf evergreen tree
  !  4 = broadleaf deciduous tree
  !  5 = tropical seasonal tree
  !  6 = cool grass (c3)
  !  7 = evergreen shrub
  !  8 = deciduous shrub
  !  9 = arctic deciduous shrub
  ! 10 = arctic grass
  ! 11 = crop
  ! 12 = irrigated crop
  ! 13 = warm grass (c4)
  ! 14 = not vegetated
  
  ! Contents of vegconi.F
  
  !++csz
  !#include <preproc.h>
  !  block data vegconi
  !--csz
  
  !++csz
  !#include <vartyp.h>
  !#include <lsmpar.h>
  !--csz
  
  ! ------------------------ code history ---------------------------
  ! source file:       vegconi.F
  ! purpose:           initialize vegetation constants 
  ! date last revised: March 1996 - lsm version 1
  ! author:            Gordon Bonan
  ! standardized:      J. Truesdale, Feb. 1996
  ! reviewed:          G. Bonan, Feb. 1996
  ! -----------------------------------------------------------------
  
  !++csz
  !#include <vegcon.h>
  
  !  integer i         !loop index
  !--csz
  
  ! plant types are:
  !  1 = needleleaf evergreen tree
  !  2 = needleleaf deciduous tree
  !  3 = broadleaf evergreen tree
  !  4 = broadleaf deciduous tree
  !  5 = tropical seasonal tree
  !  6 = cool grass (c3)
  !  7 = evergreen shrub
  !  8 = deciduous shrub
  !  9 = arctic deciduous shrub
  ! 10 = arctic grass
  ! 11 = crop
  ! 12 = irrigated crop
  ! 13 = warm grass (c4)
  ! 14 = not vegetated
  
  ! note: there are two types of crops: non-irrigated and irrigated.
  ! this is strictly for the purposes of soil hydrology (need to know 
  ! when to irrigate). for convenience, this is keyed off vegetation 
  ! type. however, plant physiology does not vary between irrigated 
  ! and non-irrigated varieties. 
  
  ! value for irrigated crop 
  data nic /12/
  
  ! value for not vegetated 
  data noveg /14/
  
  ! exponent when calculating btran: [(h2osoi-watdry)/(watopt-watdry)]**vw 
  data vw /14*1.0/
  
  ! binvvt = 1/vkc*ln(z0m/z0h)
  data binvvt /14* 0.0/
  
  ! momentum roughness length (m) 
  data z0mvt /0.94,0.77,2.62,1.10,0.99,0.06,0.06, &
       &            0.06,0.06,0.06,0.06,0.06,0.06,0.00/
  
  ! displacement height (m) 
  data zpdvt /11.39,9.38,23.45,13.40,12.06,0.34,0.34, &
       &             0.34,0.34, 0.34, 0.34, 0.34,0.34,0.00/
  
  ! characteristic leaf dimension (m) 
  data dleaf /13*0.04, 0.00/
  
  ! photosynthetic pathway: c3 = 1, c4 = 0
  data c3psn /12*1., 0., 1./
  
  ! co2 michaelis-menten constant at 25c (pa) 
  data kc25 /14*30./
  
  ! q10 for kc25 
  data akc /14*2.1/
  
  ! o2 michaelis-menten constant at 25c (pa) 
  data ko25 /14*30000./
  
  ! q10 for ko25 
  data ako /14*1.2/
  
  ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
  data vcmx25 /33.,33.,50.,33.,50.,33.,17.,17., &
       &             33.,33.,50.,50.,33., 0./
  
  ! q10 for vcmx25 
  data avcmx /14*2.4/
  
  ! minimum leaf conductance (umol/m**2/s) 
  data bp /13*2000., 1.e15/
  
  ! slope for conductance-to-photosynthesis relationship 
  data mp /6.,6.,9.,9.,9.,9.,9.,9.,9.,9.,9.,9.,5.,9./
  
  ! quantum efficiency at 25c (umol co2 / umol photon)
  data qe25 /12*0.06, 0.04, 0.00/
  
  ! q10 for qe25 
  data aqe /14*1.0/
  
  ! foliage maintenance respiration rate at 25c (umol co2 /m**2 /s)
  data rmf25 /0.50,0.50,0.75,0.50,0.75,0.50,0.26, &
       &            0.26,0.50,0.50,0.75,0.75,0.82,0.00/
  
  ! stem maintenance respiration at 25c (umol co2/kg biomass/s)
  data rms25 /0.9396,0.1364,0.1622,0.0198,0.0227,0.0000,0.0000, &
       &            0.0000,1.0230,1.0230,0.0000,0.0000,0.0000,0.0000/
  
  ! stem biomass (kg /m**2)
  data stembvt / 3.6, 3.6, 9.0, 6.2, 4.5, 0.0, 0.0, &
       &               0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0/
  
  ! root maintenance respiration at 25c (umol co2/kg biomass/s)
  data rmr25 /0.3637,0.0530,0.0455,0.0088,0.2091,0.5911,0.0000, &
       &            0.0000,2.1142,2.1142,0.0000,0.0000,2.2733,0.0000/
  
  ! root biomass (kg /m**2)
  data rootbvt /7.2, 7.2,18.0,12.4, 9.0, 0.3, 0.0, &
       &              0.0, 0.4, 0.4, 0.0, 0.0, 0.3, 0.0/
  
  ! q10 for maintenance respiration
  data arm /14*2.0/
  
  ! co2-to-biomass conversion factor (ug biomass / umol co2) 
  data dmcf /13*28.5, 0./
  
  ! foliage nitrogen concentration (%)
  data folnvt /13*2.0, 0.0/
  
  ! foliage nitrogen concentration when f(n)=1 (%) 
  data folnmx /13*1.5, 0.0/
  
  ! minimum temperature for photosynthesis (kelvin) 
  data tmin /268.16,268.16,278.16,273.16,278.16,273.16,268.16, &
       &           273.16,273.16,273.16,273.16,273.16,273.16,  0.00/
  
  ! empirical parameter that defines root fraction decrease with depth 
  data rdp /0.94,0.94,0.94,0.94,0.97,0.97,0.97, &
       &          0.97,0.94,0.94,0.94,0.94,0.97,1.00/
  
  ! microbial respiration parameter (umol co2 /kg c /s)
  data mrp /0.3727, 0.3727, 0.2333, 0.4000, 0.1250, 0.1700, 0.1909, &
       &          0.1909, 0.0500, 0.0500, 0.2273, 0.2273, 0.1700, 0.0000/
  
  ! soil carbon (kg c /m**2)
  data soilcvt /11.0, 11.0, 15.0, 11.0,  8.0, 10.0, 11.0, &
       &              11.0, 18.0, 18.0, 11.0, 11.0, 10.0,  0.0/
  
  ! leaf reflectance: 1=vis, 2=nir 
  data (rhol(i,1),i=1,mvt) /0.07,0.07,0.10,0.10,0.10,0.11, &
       &                          0.07,0.10,0.10,0.11,0.11,0.11, &
       &                          0.11,0.00/
  data (rhol(i,2),i=1,mvt) /0.35,0.35,0.45,0.45,0.45,0.58, &
       &                          0.35,0.45,0.45,0.58,0.58,0.58, &
       &                          0.58,0.00/
  
  ! stem reflectance: 1=vis, 2=nir 
  data (rhos(i,1),i=1,mvt) /0.16,0.16,0.16,0.16,0.16,0.36, &
       &                          0.16,0.16,0.16,0.36,0.36,0.36, &
       &                          0.36,0.00/
  data (rhos(i,2),i=1,mvt) /0.39,0.39,0.39,0.39,0.39,0.58, &
       &                          0.39,0.39,0.39,0.58,0.58,0.58, &
       &                          0.58,0.00/
  
  ! leaf transmittance: 1=vis, 2=nir 
  data (taul(i,1),i=1,mvt) /0.05,0.05,0.05,0.05,0.05,0.07, &
       &                          0.05,0.05,0.05,0.07,0.07,0.07, &
       &                          0.07,0.00/
  data (taul(i,2),i=1,mvt) /0.10,0.10,0.25,0.25,0.25,0.25, &
       &                          0.10,0.25,0.25,0.25,0.25,0.25, &
       &                          0.25,0.00/
  
  ! stem transmittance: 1=vis, 2=nir 
  data (taus(i,1),i=1,mvt) /0.001,0.001,0.001,0.001,0.001, &
       &                          0.220,0.001,0.001,0.001,0.220, &
       &                          0.220,0.220,0.220,0.000/
  data (taus(i,2),i=1,mvt) /0.001,0.001,0.001,0.001,0.001, &
       &                          0.380,0.001,0.001,0.001,0.380, &
       &                          0.380,0.380,0.380,0.000/
  
  ! leaf/stem orientation index: valid range = -0.4 to 0.6 
  data xl / 0.01, 0.01, 0.10, 0.25, 0.01,-0.30, 0.01, &
       &          0.25, 0.25,-0.30,-0.30,-0.30,-0.30, 0.00/
  
  ! empirical canopy wind parameter
  data cwpvt /14*3.0/
  
  ! maximum intercepted h2o per unit lai+sai (mm) 
  data ch2op /14*0.1/
  
  ! monthly leaf area index + stem area index, one-sided
  data (tai(1,i),i=1,12) &
       /4.5,4.7,5.0,5.1,5.3,5.5,5.3,5.3,5.2,4.9,4.6,4.5/
  data (tai(2,i),i=1,12) &
       /0.3,0.3,0.3,1.0,1.6,2.4,4.3,2.9,2.0,1.3,0.8,0.5/
  data (tai(3,i),i=1,12) &
       /5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0/
  data (tai(4,i),i=1,12) &
       /0.4,0.4,0.7,1.6,3.5,5.1,5.4,4.8,3.8,1.7,0.6,0.4/
  data (tai(5,i),i=1,12) &
       /1.2,1.0,0.9,0.8,0.8,1.0,2.0,3.7,3.2,2.7,1.9,1.2/
  data (tai(6,i),i=1,12) &
       /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(7,i),i=1,12) &
       /1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3/
  data (tai(8,i),i=1,12) &
       /1.0,1.0,0.8,0.3,0.6,0.0,0.1,0.3,0.5,0.6,0.7,0.9/
  data (tai(9,i),i=1,12) &
       /0.1,0.1,0.1,0.1,0.1,0.3,1.5,1.7,1.4,0.1,0.1,0.1/
  data (tai(10,i),i=1,12) &
       /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(11,i),i=1,12) &
       /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (tai(12,i),i=1,12) &
       /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (tai(13,i),i=1,12) &
       /0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8/
  data (tai(14,i),i=1,12) &
       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
  
  ! monthly leaf area index, one-sided
  data (gai(1,i),i=1,12) &
       /4.1,4.2,4.6,4.8,4.9,5.0,4.8,4.7,4.6,4.2,4.0,4.0/
  data (gai(2,i),i=1,12) &
       /0.0,0.0,0.0,0.6,1.2,2.0,2.6,1.7,1.0,0.5,0.2,0.0/
  data (gai(3,i),i=1,12) &
       /4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5/
  data (gai(4,i),i=1,12) &
       /0.0,0.0,0.3,1.2,3.0,4.7,4.5,3.4,1.2,0.3,0.0,0.0/
  data (gai(5,i),i=1,12) &
       /0.8,0.7,0.4,0.5,0.5,0.7,1.7,3.0,2.5,1.6,1.0,1.0/
  data (gai(6,i),i=1,12) &
       /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(7,i),i=1,12) &
       /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
  data (gai(8,i),i=1,12) &
       /0.9,0.8,0.2,0.2,0.0,0.0,0.0,0.2,0.4,0.5,0.6,0.8/
  data (gai(9,i),i=1,12) &
       /0.0,0.0,0.0,0.0,0.0,0.2,1.4,1.2,0.0,0.0,0.0,0.0/
  data (gai(10,i),i=1,12) &
       /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(11,i),i=1,12) &
       /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (gai(12,i),i=1,12) &
       /0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0/
  data (gai(13,i),i=1,12) &
       /0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4/
  data (gai(14,i),i=1,12) &
       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
  
  ! top of canopy (m)
  data hvt /17.0,14.0,35.0,20.0,18.0, 0.5, 0.5, &
       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0/
  
  ! bottom of canopy (m)
  data hvb / 8.50, 7.00, 1.00,11.50,10.00, 0.01, 0.10, &
       0.10, 0.10, 0.01, 0.01, 0.01, 0.01, 0.00/
  
end module dstlsm ! [mdl] LSM data

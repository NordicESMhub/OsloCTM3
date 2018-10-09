# ASAD-style ratefile for defining STRATOSPHERIC loss parameters
# Read in in the INPTR subroutine (pchem.f)    Oliver (22/10/97)
# Method Flag:   1  Frozen, no strat. chemistry for this species
#                2  First-order loss to fixed mixing ratio (or zero)
#                3  Half frozen, half scale
#                4  Auto scale to O3Strat distribution (no chemistry)
#                5  Initial scaling and strat flux
# Species   Meth  Tendency    Strat.    (Strat.=strat VMR, or frac of O3Strat)
  O3         3    2.315E-06   0.d0       ! ???? 3
  HO2        1    4.630E-05   0.d0       ! 6 hour loss rate (2 steps) *
  NO         1    2.315E-06   0.d0
  NO3        1    4.630E-05   0.d0       ! 6 hour loss rate (2 steps)
  NO2        1    2.315E-06   0.d0
  N2O5       1    2.315E-06   0.d0       ! 10% of NOy
  HONO2      1    2.315E-06   0.d0
  HO2NO2     1    2.315E-06   0.d0       ! 10% of NOy
  H2O2       1    5.787E-06   0.d0       ! 2 day loss rate *
  CH4        2    3.964E-09   0.d0       ! 8 year loss rate - estimate!
  MeOO       2    3.858E-07   0.d0       ! 30 day loss rate
  MeOOH      2    3.858E-07   0.d0       ! 30 day loss rate
  HCHO       2    3.858E-07   0.d0       ! 30 day loss rate
  CO         2    3.858E-07   0.d0       ! 30 day loss rate
  C2H6       2    3.858E-07   0.d0       ! 30 day loss rate
  EtOO       2    3.858E-07   0.d0       ! 30 day loss rate
  EtOOH      2    3.858E-07   0.d0       ! 30 day loss rate
  MeCHO      2    3.858E-07   0.d0       ! 30 day loss rate
  MeCO3      2    3.858E-07   0.d0       ! 30 day loss rate
  PAN        2    3.858E-07   0.d0       ! 30 day loss rate
  Alkane     2    3.858E-07   0.d0       ! 30 day loss rate
  Alkene     2    3.858E-07   0.d0       ! 30 day loss rate
  Aromatic   2    3.858E-07   0.d0       ! 30 day loss rate
  Isoprene   2    3.858E-07   0.d0       ! 30 day loss rate
  MVKMACR    2    3.858E-07   0.d0       ! 30 day loss rate
  MeONO2     1    4.630E-05   0.d0
  EtONO2     1    4.630E-05   0.d0
  O3Strat    3    0.0         0.d0
  e30strat   1    0.0         0.d0

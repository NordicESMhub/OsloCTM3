# RATEFILE GENERATED USING *SELECT* AND ACMSU REACTION DATABASE                                     
# MASTER RATEFILE: trimol.d                                                                         
# REACTION NETWORK: OPEN                                                                            
# TRIMOLECULAR REACTIONS - MASTER RATEFILE - Paul Brown, Oliver Wild, David Rowley & Glenn Carver   
# Centre for Atmospheric Science, Cambridge, U.K.  Release date:  11th May 1994                     
# SCCS version information: @(#)trimol.d	1.6 5/12/94                                                
    1 HO2        HO2        H2O2       O2            0.00 1.70E-33  0.00  -1000.0 0.00E+00  0.00 0.00E+00  D1,JPL         
    2 HO2        NO2        HO2NO2     m             0.60 2.00E-31 -3.40      0.0 2.90E-12 -1.10 0.00E+00  JPL         
    3 HO2NO2     m          HO2        NO2           0.60 2.00E-31 -3.40  10900.0 2.90E-12 -1.10 2.10E-27  JPL
    4 N2O5       m          NO2        NO3           0.60 2.00E-30 -4.40  11000.0 1.40E-12 -0.70 2.70E-27  JPL
    5 NO2        NO3        N2O5       m             0.60 2.00E-30 -4.40      0.0 1.40E-12 -0.70 0.00E+00  JPL         
    6 OH         NO2        HONO2      m             0.60 1.80E-30 -3.00      0.0 2.80E-11  0.00 0.00E+00  JPL         
    7 OH         OH         H2O2       m             0.60 6.90E-31 -1.00      0.0 2.60E-11  0.00 0.00E+00  JPL         
    8 MeCO3      NO2        PAN        m             0.60 9.70E-29 -5.60      0.0 9.30E-12 -1.50 0.00E+00  JPL         
    9 PAN        m          MeCO3      NO2           0.60 9.70E-29 -5.60  14000.0 9.30E-12 -1.50 9.00E-29  JPL
   10 Alkene     OH         ROHOO      HO2           0.60 1.00E-28 -4.50      0.0 8.80E-12  0.85 0.00E+00  JPL         
 9999                                                0.00 0.00E+00  0.00      0.0 0.00E-00  0.00    
                                                                                 
                                                                                 
The following removed - non-modelled and rapidly decomposing products                                             
    1 HO2        HCHO       HOCH2OO                  0.00 0.00E+00  0.00      0.0 0.00E+00  0.00      0.0           
    5 MeOO       NO2        MeO2NO2    m             0.36 2.50E-30 -5.50      0.0 7.50E-12  0.00      0.0           
   13 OH         NO         HONO       m          1420.00 7.40E-31 -2.40      0.0 3.30E-11  0.00      0.0           
                                                                                 
                                                                                 
                                                                                 
                                                                                 
 NOTES:                                                                          
 -----                                                                           
  All reaction data taken from IUPAC supplement IV unless                        
  otherwise indicated.                                                           
                                                                                 
  JPL - data from JPL 06-2
                                                                                 
  ? - reaction products unknown                                                  
  * - user strongly advised to consult source material                           
  B - branching ratio assumed equal for all channels in the                      
       absence of more information                                               
  F - parameters given are not sufficient for a full                             
       calculation of Fc - see source material -                                 
       however the additional term is small in comparison                        
       (eg. NO2+NO3->N2O5 and reverse reaction                                   
         Full calculation: Fc=exp(-T/250)+exp(-1050/T)                           
          at 200 K this term=0.45 additional term=0.005                          
          at 300 K          =0.30                =0.03  )                        
  M - reaction does not actually involve use of third body                       
  PD - rate at 1 bar; pressure dependence unknown                                
  PB - rate at 1 bar; pressure dependence unknown; these                         
       reactions are written as bimolecular reactions in                         
       the IUPAC assessment                                                      
  U - upper limit for rate coefficient                                           
  A1 - two rates are given for [O2] & [N2}. Rate is calculated by                
       ( ko2x[0.21] + kn2x[0.78] ) / ( [0.21]+[0.78] )                           
  Alkene - C2H4
                                                                                 
  Dependent reactions:                                                           
                                                                                 
  D1 - Depends on the concentration of H2O. See JPL 1992 and the paper           
       it references: R.R.Lii et al, J.Phys.Chem 85, 1981, p2833. This           
       reaction (and the bimolecular branch) need to be multiplied by            
       the following factor:                                                     
       (1 + 1.4E-21[H2O]exp(2200/T).                                             
                                                                                 
 Changes since 22/11/93 release:                                                 
 (1) Corrected virtually all rates. Previously they had been multiplied          
     by [N2] when written as such in IUPAC, but this is incorrect as             
     the rate is assumed to apply equally in O2 unless otherwise specified.      
 (2) Added note that HO2+HO2 reaction is H2O dependent.                          
 (3) Used high pressure limit for O3(3P)+NO2 from IUPAC J. Phys. Ref. Data.       
     This is different from that reported in Atm. Env.                           
 Changes since 08/3/93 release:                                                  
 (1) O now written as O3(3P)                                                      
 (2) some comments altered                                                       


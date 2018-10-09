# RATEFILE GENERATED USING *SELECT* AND ACMSU REACTION DATABASE                                     
# MASTER RATEFILE: bimol.d                                                                          
# REACTION NETWORK: OPEN                                                                            
# BIMOLECULAR REACTIONS - MASTER RATEFILE - Paul Brown, Oliver Wild & David Rowley                  
# Centre for Atmospheric Science, Cambridge, U.K.  Release date:  10th May 1994                     
# SCCS version information: @(#)bimol.d	1.4 5/11/94
# Updated with JPL 06-2 data on July 24, 2008 by Qi Tang                                                 
    1 HO2        HO2        H2O2       O2                   3.50E-13  0.00   -430.0  D1,JPL       
    2 HO2        MeOO       O2         MeOOH                4.10E-13  0.00   -750.0  JPL         
    3 HO2        NO         OH         NO2                  3.50E-12  0.00   -250.0  JPL         
    4 HO2        O3         OH         O2         O2        1.00E-14  0.00    490.0  JPL         
    5 MeOO       MeOO       HCHO       HO2        O2        9.50E-14  0.00   -390.0  b3       
    $                       2.0        2.0        1.0                                               
    6 MeOO       NO         HCHO       NO2        HO2       2.80E-12  0.00   -300.0  JPL         
    7 N2O5       H2O        HONO2      HONO2                2.00E-21  0.00      0.0  U        
    8 NO         NO3        NO2        NO2                  1.50E-11  0.00   -170.0  JPL         
    9 NO         O3         NO2        O2                   3.00E-12  0.00   1500.0  JPL         
   10 NO2        O3         NO3        O2                   1.20E-13  0.00   2450.0  JPL         
   11 O3         H2O        OH         OH                   1.63E-10  0.00    -60.0  JPL         
   12 OH         CH4        H2O        MeOO                 2.45E-12  0.00   1775.0  JPL         
   13 OH         CO         HO2        CO2                  1.50E-13  0.00      0.0  P1       
   14 OH         H2O2       H2O        HO2                  1.80E-12  0.00      0.0  JPL         
   15 OH         HCHO       H2O        CO         HO2       5.50E-12  0.00   -125.0  JPL      
   16 OH         HO2        H2O        O2                   4.80E-11  0.00   -250.0           
#   17 OH         HO2NO2     ?                               1.50E-12  0.00   -360.0  ?        
   17 OH         HO2NO2     H2O        NO2                  1.30E-12  0.00   -380.0  JPL         
   18 OH         HONO2      H2O        NO3                  2.06E-13  0.00      0.0  P2,JPL       
   19 OH         MeOOH      H2O        HCHO       OH        1.14E-12  0.00   -200.0  JPL        
   20 OH         MeOOH      H2O        MeOO                 2.66E-12  0.00   -200.0  JPL        
   21 OH         NO3        HO2        NO2                  2.20E-11  0.00      0.0  JPL        
   22 OH         O3         HO2        O2                   1.70E-12  0.00    940.0  JPL         
   23 OH         OH         H2O        O3                   1.80E-12  0.00      0.0  JPL         
   24 C2H6       OH         EtOO                            8.70E-12  0.00   1070.0  JPL         
   25 EtOO       NO         MeCHO      HO2        NO2       2.60E-12  0.00   -365.0  JPL         
   26 EtOO       EtOO       MeCHO      HO2                  6.80E-14  0.00      0.0  JPL         
    $                       2.0        2.0                                                    
   27 EtOO       MeOO       MeCHO      HCHO       HO2       2.50E-14  0.00      0.0 ?         
    $                       1.0        1.0        2.0                                         
   28 MeCHO      OH         MeCO3      MeOO       CO        5.60E-12  0.00   -310.0           
    $                       0.50       0.5        0.5                               de-PAN!  
   29 MeCHO      NO3        MeCO3      HONO2                1.40E-12  0.00   1900.0  JPL         
   30 MeCO3      NO         MeOO       NO2                  8.10E-12  0.00   -270.0  JPL         
   31 MeOO       MeCO3      HCHO       MeOO       HO2       2.00E-12  0.00   -500.0  JPL        
   32 MeCO3      MeCO3      MeOO       MeOO                 2.90E-12  0.00   -500.0  JPL        
   33 EtOO       HO2        EtOOH                           7.50E-13  0.00   -700.0  JPL        
   34 EtOOH      OH         EtOO       MeCHO      OH        1.40E-11  0.00      0.0          
    $                       0.25       0.75       0.75                                       
#   35 Alkane     OH         EtOO       EtOO                 1.55E-11  0.00    540.0  Harwell 
   35 Alkane     OH         EtOO       MeOO       CO        8.70E-11  0.00    615.0  JPL
   36 Alkene     O3         HCHO       MeCHO      CO        1.20E-14  0.00   2630.0  JPL
    $                       1.0        0.5        1.0                                        
   37 ROHOO      NO         HCHO       MeCHO      NO2       4.20E-12  0.00   -180.0  H gen   
   38 ROHOO      HO2        HCHO       MeCHO      OH        6.50E-13  0.00   -650.0  H gen   
   39 ROHOO      MeOO       HCHO       MeCHO      HO2       1.00E-13  0.00      0.0  BI gen  
    $                       2.0        1.0        1.0                                        
   40 Aromatic   OH         ArOO       HO2                  1.40E-11  0.00      0.0  MCM     
   41 ArOO       NO         MeCO3      CO         NO2       4.20E-12  0.00   -180.0  H gen   
    $                       1.0        4.0        1.0                                de-PAN! 
   42 ArOO       HO2        MeCO3      CO         OH        6.50E-13  0.00   -650.0  H gen   
    $                       1.0        4.0        1.0                                de-PAN! 
   43 Isoprene   O3         MVKMACR    HCHO       OH        7.86E-15  0.00   1913.0  MCM     
   44 Isoprene   OH         IsopOO                          2.54E-11  0.00   -410.0  MCM     
   45 IsopOO     NO         MVKMACR    HCHO       NO2       4.20E-12  0.00   -180.0  H gen   
   46 IsopOO     HO2        EtOOH      CO         CO        6.50E-13  0.00   -650.0  H gen   
   47 IsopOO     MeOO       MVKMACR    HCHO       HO2       1.00E-13  0.00      0.0  BI gen  
    $                       1.0        2.0        2.0                                        
   48 MVKMACR    O3         MeCO3      HCHO       MeOO      7.51E-16  0.00   1521.0  MCM MVK 
    $                       0.5        1.0        1.0                                de-PAN! 
   49 MVKMACR    OH         MVKOO                           4.13E-12  0.00   -452.0  MCM MVK 
   50 MVKOO      NO         MeCO3      HCHO       NO2       4.20E-12  0.00   -180.0  H gen   
    $                       0.5        1.0        1.0                                de-PAN! 
   51 MVKOO      HO2        EtOOH      CO                   6.50E-13  0.00   -650.0  H gen   
   52 MVKOO      MeOO       MeCO3      HCHO       HO2       1.00E-13  0.00      0.0  BI gen  
    $                       1.0        2.0        1.0                                        
   53 MeONO2     OH         NO2        HCHO       HO2       4.00E-13  0.00    845.0  IUPAC
   54 EtONO2     OH         NO2        MeCHO      HO2       6.70e-13  0.00    395.0  IUPAC
 9999                                                       0.00E+00  0.00      0.0                 


The following removed - non-modelled and decomposing products
    1 HO2        HCHO       MeOO                            9.70E-15  0.00   -625.0           
    7 MeOO       MeOO       MeOH       HCHO       O2        6.60E-14  0.00   -365.0  b3       
    9 MeOO       MeOO       MeOOMe     O2                   1.10E-14  0.00   -365.0  b3       


Former NMHC scheme neglecting peroxy radicals
   56 Alkane     OH         EtOO       EtOO                 1.55E-11  0.00    540.0          
   57 Alkene     OH         HCHO       MeCHO      HO2       4.10E-12  0.00   -545.0          
   58 Alkene     O3         HCHO       MeCHO      CO        1.20E-14  0.00   2633.0          
    $                       1.0        0.5        1.0                                        
   59 Aromatic   OH         MeCO3      CO         HO2       2.10E-12  0.00   -322.0   1.10E-11  0.00      0.0 ??
    $                       2.0        4.0        1.0                                        
   60 Isoprene   OH         MeCO3      HCHO       CO        1.50E-11  0.00   -500.0          
    $                       1.0        2.0        1.0                                        
   61 Isoprene   O3         MeCO3      HCHO       CO        7.00E-15  0.00   1900.0          
    $                       1.0        2.0        1.0                                        
 9999                                                       0.00E+00  0.00      0.0                 
                                                                                 
                                                                                 
                                                                                 
 NOTES:                                                                          
 -----                                                                           
  All reaction data taken from IUPAC supplement IV unless                        
  otherwise indicated.                                                           
                                                                                 
 General comments                                                                
  JPL - data from JPL (latest assessment used:  2006)                            
  ? - reaction products unknown                                                  
  * - user strongly advised to consult source material                           
  B - branching ratio assumed equal for all channels                             
       in the absence of more information                                        
  C - temperature dependent branching ratio based on                             
       ratio st 298 K                                                            
  U - upper limit for rate coefficient                                           
  n - normal form assumed (rather than iso-)                                     
       in the absence of more information
  Alkane - C3H8
  Alkene - C2H4                                        
                                                                                 
 Specific comments                                                               
  r1 - OH + H1403t (methyl chloroform) - data from Talukdar et al. 1992          
  r2 - BrO + HO2 - rate data from Poulet et al. GRL 19, 2305. See also JPL 1992. 
        Branch giving HBr likely to be small, assumed zero.                      
  b1 - branching ratio for ClO+HO2->HOCL+O2/HCl+O3 - second branch               
       has an upper limit of between 0.3-1.5 %. Set to zero here but             
       retained in ratefile in case a non-zero value is determined later.        
  b2 - branching ratio for ClO+NO3->ClOO+NO2/OClO+NO2 - IUPACIV states           
       that the first branch will dominate at temperatures less than 300 K       
       Second branch set to zero here but retained in ratefile in case a         
       non-zero value is determined later.                                       
  b3 - brancing ratio for MeOO+MeOO->MeOH+HCHO+O2/MeO+MeO+O2/MeOOMe+O2           
       is 0.6:0.3:0.1 at 298 K according to JPL92. IUPAC states the ratio        
       is temperature dependent but the JPL ratios retained in the absence       
       of further information.                                                   
  b4 - branching ratio for OH+ClO->HO2+Cl/HCl+O2. IUPACIV and JPL92 give         
       upper limits for second branch of ~2-14 %, and state that zero            
       yield is possible within experimental errors. Set to zero here but        
       retained in ratefile in case a non-zero value is determined later.        
                                                                                 
  Pressure dependent reactions - require extra code in model                     
  to calculate rate coefficients - see source material                           
                                                                                 
  P1 - OH + CO -> H + CO2  k=k(1+0.6*P/bar)                                      
  P2 - OH + HONO2 -> H2O + NO3  (2.06E-13 = rate at 1 bar & 298 K)               
       Generally, k=k1(T)+k2(T,M) where k2(M,T)=k3[M]/(1+k3[M]/k4)               
       and k1=2.4e-14*exp(460/T)  cm3s-1                                         
           k3=6.5e-34*exp(1335/T) cm6s-1    T is temperature (K)                 
           k4=2.7e-17*exp(2199/T) cm3s-1   [M] is tot.no.density (cm-3)          
  P3 - Et + O2 -> C2H4 + HO2                                                     
       k=3.80e-15  at 1 bar                                                      
       k=1.90e-14  at 0.133 bar                                                  
                                                                                 
  Other dependencies:                                                            
                                                                                 
  D1 - Depends on the concentration of H2O. See JPL 1992 and the paper           
       it references: R.R.Lii et al, J.Phys.Chem 85, 1981, p2833. This           
       reaction (and the trimolecular branch) need to be multiplied by           
       the following factor:                                                     
       (1 + 1.4E-21[H2O]exp(2200/T).                                             
                                                                                 
                                                                                 
  Note: there are some other pressure dependent reactions which are              
  listed as though bimolecular reactions in the IUPAC assessment,                
  but which have been included in the trimolecular ratefile                      
  (marked 'PB'). Their pressure dependence is generally unknown.                 
                                                                                 
                                                                                 
 Changes since 22/11/93 release:                                                 
 (1) Added note that HO2+HO2 reaction depends on H2O.                            
 (2) Corrected BrO+HO2 rate.                                                     
 Changes since 11/10/93 release:                                                 
 (1) Revised rate/branching ratio for BrO+HO2                                    
 Changes since 24/8/93 release:                                                  
 (1) some branching ratios included                                              
 Changes since 08/3/93 release:                                                  
 (1) O now written as O3(3P)                                                     
 (2) some comments altered                                                       

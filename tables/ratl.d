# ASAD-style ratefile for defining extra products for NMHC reactions
# Read in from inrats                              Oliver (18/02/98)
#
# Define species and factors to use; reactions are defined by the
# reactants, not by the reaction numbers.
#
No. Species:  4             NO         NO2        OH         HO2       
   56 Alkane     OH          0.0        0.0       -1.0        0.0      
   61 Aromatic   OH          0.0        0.0       -3.0       +2.0      
   62 ArOO       NO         -1.0       +1.0        0.0        0.0      
   63 ArOO       HO2         0.0        0.0       +1.0       -1.0      
   66 IsopOO     NO          0.0        0.0        0.0       +1.0      
   70 MVKMACR    OH          0.0        0.0       -1.0       +1.0      
 9999                                                                  

Old scheme without peroxyradicals:

No. Species:  4             NO         NO2        OH         HO2       
   59 Alkene     OH         -1.0       +1.0        0.0        0.0      
   61 Aromatic   OH         -2.0       +2.0       -4.0       +2.0      
   62 Isoprene   OH         -2.0       +2.0       -2.0       +1.0      
   63 Isoprene   O3         -1.0       +1.0       -2.0       +1.0      
 9999                                                                  

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 20:25:59 2020

@author: wangh
"""


import numpy as np




    
    

# 17 variables are assigned for each canopy type 
# 1  = canopy depth
# 2  = leaf width
# 3  = leaf length
# 4  = canopy height
# 5  = scattering coefficient for PPFD
# 6  = scattering coefficient for near IR
# 7  = reflection coefficient for diffuse PPFD
# 8  = reflection coefficient for diffuse near IR
# 9  = clustering coefficient (accounts for leaf clumping influence on mean
#    projected leaf area in the direction of the suns beam)
# 10 = leaf IR emissivity
# 11 = leaf stomata and cuticle factor: 1=hypostomatous, 2=amphistomatous,
#     1.25=hypostomatous but with some transpiration through cuticle
# 12 = daytime temperature lapse rate (K m-1)
# 13 = nighttime temperature lapse rate (K m-1)
# 14 = warm (>283K) canopy total humidity change (Pa)
# 15 = cool (>= 283K) canopy total humidity change (Pa)
# 16 = normalized canopy depth where wind is negligible
# 17 = canopy transparency
#
# Six canopy types currently used in MEGCAN: 
# 1  = Needleleaf trees
# 2  = Tropical forest trees, 
# 3  = Temperate broadleaf trees
# 4  = shrubs
# 5  = herbaceous
# 6  = crops    
      
Canopychar = np.array(\
[[  16.   , 16.   , 16.   ,  1.   ,  0.5  ,  1.  ],\
[  0.005 ,  0.05 ,  0.05 ,  0.015,  0.01 ,  0.02 ],\
[   0.1  ,  0.1  ,  0.1  ,  0.1  ,  0.15 ,  0.15 ],\
[  24.   , 24.   , 24.   ,  2.   ,  0.5  ,  1.0  ],\
[   0.2  ,  0.2  ,  0.2  ,  0.2  ,  0.2  ,  0.2  ],\
[   0.8  ,  0.8  ,  0.8  ,  0.8  ,  0.8  ,  0.8  ],\
[   0.057,  0.057,  0.057,  0.057,  0.057,  0.057],\
[   0.389,  0.389,  0.389,  0.389,  0.389,  0.389],\
[   0.85 ,  1.1  ,  0.9  ,  0.85 ,  0.7  ,  0.65],\
[   0.95 ,  0.95 ,  0.95 ,  0.95 ,  0.95 ,  0.95 ],\
[   1.25 ,  1.25 ,  1.25 ,  1.   ,  1.25 ,  1.25 ],\
[   0.06 ,  0.06 ,  0.06 ,  0.06 ,  0.06 ,  0.06 ],\
[  -0.06 , -0.06 , -0.06 , -0.06 , -0.06 , -0.06 ],\
[ 700.   ,700.   ,700.   ,700.   ,700.   ,700.   ],\
[ 150.   ,150.   ,150.   ,150.   ,150.   ,150.   ],\
[   0.7  ,  0.7  ,  0.7  ,  0.7  ,  0.7  ,  0.7  ],\
[   0.2  ,  0.2  ,  0.2  ,  0.2  ,  0.2  ,  0.2  ]])
# Canopy characteristics for MEGCAN canopy types 
NrTyp = 6  	# Number of canopy types
NrCha = 17  	# Number of canopy characteristics


def GaussianDist(NLayers):
    #Number of layers
    
    Distgauss = np.zeros(NLayers)
    if NLayers == 1:
        Distgauss[0]   = 0.5
    elif NLayers == 3:
        Distgauss[0]   = 0.112702
        Distgauss[1]   = 0.5
        Distgauss[2]   = 0.887298
    elif NLayers == 5:
        Distgauss[0]   = 0.0469101
        Distgauss[1]   = 0.2307534
        Distgauss[2]   = 0.5
        Distgauss[3]   = 0.7692465
        Distgauss[4]   = 0.9530899
    else:
        for i in range(NLayers):
          Distgauss[i]   = (i + 1 - 0.5)/NLayers

    return Distgauss

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION Stability
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def Stability(Cantype, Solar):

#----------------------------------------------------------------
    Trateboundary = 500

    if Solar > Trateboundary:
        # Daytime temperature lapse rate
        Stability = Canopychar[11, Cantype]
    elif Solar > 0 and Solar <= Trateboundary:
        Stability = Canopychar[11, Cantype] - \
        ((Trateboundary - Solar) / Trateboundary)*(Canopychar[11, Cantype] - Canopychar[12, Cantype])
    else:
        # Nightime temperature lapse rate
        Stability = Canopychar[12, Cantype]

    return Stability 


#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION WaterVapPres
#
#   Convert water mixing ratio (kg/kg) to water vapor pressure 
#   (Pa or Kpa depending on units of input )
#   Mixing ratio (kg/kg), temp (C), pressure (KPa)
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def WaterVapPres(Dens, Pres):
    WaterAirRatio       = 18.016/28.97
#----------------------------------------------------------------
    WaterVapPres = (Dens / (Dens + WaterAirRatio)) * Pres
    return WaterVapPres

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION RHtoWaterVapPres
#
#   Convert water mixing ratio (kg/kg) to water vapor pressure 
#   (Pa or Kpa depending on units of input )
#   Mixing ratio (kg/kg), temp (C), pressure (KPa)
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



def RHtoWaterVapPres(RH,T):
    t = T- 273.15
    #Tetens equation
    es= 6.112*np.exp(17.67*t/(t+243.5))
#----------------------------------------------------------------
    WaterVapPres = es*RH*0.01*1000
    return WaterVapPres

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION RHtoWaterVapPres
#
#   Convert water mixing ratio (kg/kg) to water vapor pressure 
#   (Pa or Kpa depending on units of input )
#   Mixing ratio (kg/kg), temp (C), pressure (KPa)
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



def WaterVapPrestoRH(WVP,T):
    t = T- 273.15
    #Tetens equation
    es= 6.112*np.exp(17.67*t/(t+243.5))
#----------------------------------------------------------------
    RH = WVP/es
    return RH



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   Subroutine CalcExtCoeff
#   Calculate light extinction coefficients
#   Code originally developed by Alex Guenther in 1990s
#   Coded into FORTRAN by Xuemei Wang
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def CalcExtCoeff(Qbeam,scat,kb,kd):
      
    P     = (1 - scat)**0.5
    reflb = 1 - np.exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))

# Extinction coefficients
    kbp   = kb * P
    kdp   = kd * P
# Absorbed beam radiation
    QbeamAbsorb = kb * Qbeam * (1 - scat) 
    Out_List = [QbeamAbsorb,reflb,kbp,kdp]
    return Out_List

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   Subroutine CalcRadComponents
#   Code originally developed by Alex Guenther in 1990s
#   Coded into FORTRAN by Xuemei Wang
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def CalcRadComponents(Qdiff, Qbeam, Kdp, Kbp, Kb,Scat, Refld, Reflb, LAIdepth):
#-------------------------------------------------------------------

    QdAbs = Qdiff * Kdp * (1 - Refld) * np.exp(-Kdp * LAIdepth)
    QsAbs = Qbeam * ((Kbp * (1 - Reflb) * np.exp(-Kbp * LAIdepth)) - \
                  (Kb * (1 - Scat) * np.exp(-Kb * LAIdepth)))
    return QdAbs,QsAbs






   
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo#
#   Subroutine CanopyRad
#
#   Canopy light environment model
#   Code originally developed by Alex Guenther in 1990s
#   Coded into FORTRAN by Xuemei Wang
#   based on Spitters et al. (1986), 
#   Goudrian and van Laar (1994), Leuning (1997)
#   Initial code 8-99, 
#   modified 7-2000, 12-2001, 1-2017
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#def CanopyRad(Distgauss, Layers, LAI, Sinbeta,\
#                Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype,\    
#                 Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV,\ 
#                 QbAbsn, QdAbsn, QsAbsn, SunQv,\        
#                 ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, \
#                 NrCha, NrTyp)
def CanopyRad(Distgauss, layers, LAI, Sinbeta,Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype):   
                 

                 
#var
#Parameter
    #1.Canopychar 
    #2.NrCha
    #3.NrTyp
#    Input
    #Distgauss, Layers, LAI, Sinbeta, Qbeamv,Qdiffv, Qbeamn, Qdiffn, Cantype 
#    Output  
    #Sunfrac, QbAbsV, QdAbsV, QsAbsV,\ 
    #QbAbsn, QdAbsn, QsAbsn, SunQv,\        
    #ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD
# Stefan-boltzman constant  W m-2 K-4
    ConvertShadePPFD = 4.6
    ConvertSunPPFD = 4.0
#---------------------------------------------------------------------
    Sunfrac = np.zeros(layers)
    QbAbsV  = np.zeros(layers)
    QdAbsV  = np.zeros(layers)
    QsAbsV  = np.zeros(layers)
    QbAbsN  = np.zeros(layers)
    QdAbsN  = np.zeros(layers)
    QsAbsN  = np.zeros(layers)
    SunQv   = np.zeros(layers)     
    ShadeQv = np.zeros(layers)
    SunQn   = np.zeros(layers)
    ShadeQn = np.zeros(layers)
    SunPPFD = np.zeros(layers)
    ShadePPFD = np.zeros(layers)
        
# adjust LAI for canopy transparency
    CANTRAN = Canopychar[16,Cantype]
    LAIadj = LAI / ( 1 - CANTRAN )
    Cluster = 0
    if(((Qbeamv  + Qdiffv ) > 0.001) and (Sinbeta  > 0.002) and (LAIadj  > 0.001)):
        #Daytime = True

# Scattering coefficients (scatV,scatN), diffuse and beam reflection 
# coefficients (ref..) for visible or near IR
        ScatV   = Canopychar[4,Cantype]
        ScatN   = Canopychar[5,Cantype]
        RefldV  = Canopychar[6,Cantype]
        RefldN  = Canopychar[7,Cantype]
        Cluster = Canopychar[8,Cantype]

# Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
        Kb = Cluster * 0.5 / Sinbeta 
# (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
        Kd = 0.8 * Cluster              
# (0.8 assumes a spherical leaf angle distribution)
        #Out_List = [QbeamAbsorb,reflb,kbp,kdp]
        Ext_list = CalcExtCoeff(Qbeamv,ScatV,Kb,Kd)
        QbAbsV =Ext_list[0]
        ReflbV =Ext_list[1]
        KbpV   =Ext_list[2]
        KdpV   =Ext_list[3]
        Ext_list = CalcExtCoeff(Qbeamn,ScatN,Kb,Kd)
        QbAbsN =Ext_list[0]
        ReflbN =Ext_list[1]
        KbpN   =Ext_list[2]
        KdpN   =Ext_list[3]


        for i in range(layers):
            # LAI depth at this layer
            LAIdepth   = LAIadj  * Distgauss[i]
            #fraction of leaves that are sunlit
            Sunfrac[i] = np.exp(-Kb * LAIdepth)
            QdAbsVL,QsAbsVL = CalcRadComponents(Qdiffv , Qbeamv , KdpV,\
                                                KbpV, Kb, ScatV, RefldV,\
                                                    ReflbV, LAIdepth)           
            QdAbsNL,QsAbsNL = CalcRadComponents(Qdiffn , Qbeamn , KdpN,\
                                                KbpN, Kb, ScatN, RefldN,\
                                                    ReflbN, LAIdepth)
            ShadePPFD[i] = (QdAbsVL + QsAbsVL) * ConvertShadePPFD/(1 - ScatV)  
            SunPPFD[i] = ShadePPFD[i] + (QbAbsV* ConvertSunPPFD/(1 - ScatV))
            QdAbsV[i] = QdAbsVL
            QsAbsV[i] = QsAbsVL
            QdAbsN[i] = QdAbsNL
            QsAbsN[i] = QsAbsNL
            ShadeQv[i] = QdAbsVL + QsAbsVL
            SunQv[i]   = ShadeQv[i] + QbAbsV 
            ShadeQn[i] = QdAbsNL + QsAbsNL
            SunQn[i]   = ShadeQn[i] + QbAbsN 
    else: # Night time

            QbAbsV   = 0
            QbAbsn   = 0
    
            Sunfrac[:]   = 0.2
            SunQn[:]     = 0
            ShadeQn[:]   = 0
            SunQv[:]     = 0
            ShadeQv[:]   = 0
            SunPPFD[:]   = 0
            ShadePPFD[:] = 0
            QdAbsV[:]    = 0
            QsAbsV[:]    = 0
            QdAbsN[:]    = 0
            QsAbsN[:]    = 0
            #QbAbsV
    ListAll = [Sunfrac,QbAbsV,QbAbsN,SunQn,ShadeQn,SunQv,ShadeQv,\
               SunPPFD,ShadePPFD,QdAbsV,QsAbsV,QdAbsN,QsAbsN]
           #SunPPFD,ShadePPFD,QdAbsV,QsAbsV,QdAbsn,QsAbsn]
    return ListAll


#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION LeafLE
#
#   Latent energy term in Energy balance
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType):

    LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
    
    Svp = 10**((-2937.4 / Tleaf) - (4.9283 * np.log10(Tleaf)) + 23.5518)  
    SvdTk = 0.2165 * Svp / Tleaf
    Vapdeficit = (SvdTk - Ambvap)

# Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
#                 leaf resistence (s m-1)
    LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
    if LE < 0:
        LeafLE = 0
    else:
        LeafLE = LE


    return  LeafLE


#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION LeafBLC
#
#   Boundary layer conductance
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def LeafBLC(GHforced, Tdelta, Llength):

#----------------------------------------------------------------

# This is based on Leuning 1995 p.1198 except using molecular 
# conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
# diffusivity so that you end up with a heat convection coefficient 
# (W m-2 K-1) instead of a conductance for free convection

      if Tdelta >= 0:
          GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta/(Llength**3.))**0.25)/Llength
      else:
          GhFree = 0
      LeafBLC = GHforced + GhFree

      return LeafBLC

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION ResSC
#
#   Leaf stomatal cond. resistance s m-1
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def ResSC(PPFD):
    SCadj = ((0.0027*1.066*PPFD)/((1 + 0.0027 * 0.0027 * PPFD**2.)**0.5))

    if SCadj < 0.1:
        StomRes = 2000
    else:
        StomRes = 200 / SCadj
    return StomRes


#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION LeafH
#
#   Convective energy term in Energy balance (W m-2 heat flux 
#      from both sides of leaf)
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def LeafH(Tdelta, GH):
#----------------------------------------------------------------

# 2 sides X conductance X Temperature gradient
    LeafH = 2 * GH * Tdelta
    return LeafH
      
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION LHV
#
#   Latent Heat of vaporization(J Kg-1) from Stull p641
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def LHV(Tk):

# REVISE - Replace 273 with 273.15

    LHV = 2501000 - (2370 * (Tk - 273.15))

    return LHV

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION LeafIR
#
#   Calculate IR transfer between leaf and air
#   Added by Alex Guenther and Ling Huang to replace previous
#   MEGAN2.1 IR balance functions
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def LeafIR(Tk, Eps):
# Stefan-boltzman constant  W m-2 K-4
    Sb = 0.0000000567 
#----------------------------------------------------------------
    LeafIR = Eps * Sb * (2 * (Tk**4.)) 
    return LeafIR



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   Subroutine LeafEB
#
#   Leaf energy balance
#   Code originally developed by Alex Guenther in 1990s
#   Coded into FORTRAN by Xuemei Wang
#   Coded into Python by Hui Wang
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#      SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType,             
#     &         Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf,           
#     &         SH, LH, IRout)

def LeafEB(PPFD, Q, IRin, Eps, TranspireType,Lwidth, Llength, TairK,\
           HumidairPa, Ws):

    if Ws <= 0:
        Ws1 = 0.001
    else:
        Ws1 = Ws

#   Air vapor density kg m-3
#   HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)
    HumidAirKgm3 = 0.002165 * HumidairPa / TairK


#   Heat convection coefficient (W m-2 K-1) for forced convection. 
#   Nobel page 366
    GHforced = 0.0259 / (0.004 * ((Llength / Ws1)**0.5))

#   Stomatal resistence s m-1
    StomRes  = ResSC(PPFD)
 
    

#   REVISE - Replace LeafIRout with LeafIR
#   IRoutairT = LeafIROut(tairK, eps)
#XJ IRoutairT  = LeafIR(TairK + Tdelt, Eps)        
    IRoutairT = LeafIR(TairK, Eps)
 
#   Latent heat of vaporization (J Kg-1)
    LatHv = LHV(TairK)
    
#   Latent heat flux
    LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes,TranspireType)

    E1 = (Q + IRin - IRoutairT - LHairT)
    if E1 == 0.:
        E1 = -1.


    Tdelt = 1.
    Balance = 10.
    for i in np.arange(1,11,1):
        
        if abs(Balance) > 2:
            
          # Boundary layer conductance
            GH1 = LeafBLC(GHforced, Tdelt, Llength)
          # Convective heat flux
            SH1 = LeafH(Tdelt, GH1)                
          # Latent heat of vaporization (J Kg-1)
            LatHv = LHV(TairK + Tdelt)             
            LH = LeafLE(TairK + Tdelt, HumidAirKgm3,LatHv, GH1, StomRes, TranspireType)
            LH1 = LH - LHairT

# REVISE - Replace LeafIROut with LeafIR
#           IRout  = LeafIROut(TairK + Tdelt, Eps)
            IRout  = LeafIR(TairK + Tdelt, Eps)
            IRout1 = IRout - IRoutairT
            Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
            Balance = Q + IRin - IRout - SH1 - LH


    if Tdelt > 10:
        Tdelt = 10
    if Tdelt < -10:
        Tdelt = -10

    Tleaf = TairK + Tdelt
    GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
    SH    = LeafH(Tleaf - TairK, GH)
    LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, GH, StomRes, TranspireType)

# REVISE - Replace LeafIROut with LeafIR
#      IRout = LeafIROut(Tleaf, Eps)
    IRout = LeafIR(Tleaf, Eps)
    List_out = [IRout, Tleaf, SH, LH, StomRes]
    #List_out = [PPFD, Tleaf, SH, LH, StomRes]
    #REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH, StomRes


    return List_out
      











#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   Subroutine CanopyEB
#
#   Canopy energy balance model for estimating leaf temperature
#   Coded into FORTRAN by Xuemei Wang
#   Code developed by Alex Guenther in 1990s
#   based on Goudrian and Laar (1994) and Leuning (1997)
#   Initial code 8-99, modified 7-2000 and 12-2001
#   Modified in 1-2017 by Alex Guenther and Ling Huang
#   to correct IR balance and atmos. emissivity
#   Note: i denotes an array containing a vertical profile 
#         through the canopy with 0 (above canopy conditions) 
#         plus 1 to number of canopy layers
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#      SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar,    
#     &             Cantype, TairK, HumidairPa, Ws,       
#     &             SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,  
#     &             Sunleaftk, SunleafSH, SunleafLH,                 
#     &             SunleafIR, Shadeleaftk, ShadeleafSH,             
#     &             ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0,     
#     &             TairK0, HumidairPa0)



def CanopyEB(Trate, Layers, Distgauss, Cantype, TairK0, Ws0,\
             SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,HumidairPa0):
    
     #            HumidairPa,Sunleaftk, SunleafSH, SunleafLH, \                 
     #            SunleafIR, Shadeleaftk, ShadeleafSH, \             
     #            ShadeleafLH, ShadeleafIR, Ws0,  \   
     #            Ws, TairK0, )


# inputs
#      INTEGER,INTENT(IN) :: NrCha, NrTyp, Layers, Cantype
#     REAL,INTENT(IN) :: Trate, TairK0, HumidairPa0, Ws0
#      REAL,DIMENSION(Layers),INTENT(IN) ::  Distgauss, SunQv,ShadeQv, 
#     &            SunQn, ShadeQn, SunPPFD, ShadePPFD
#      REAL,DIMENSION(NrCha, NrTyp),INTENT(IN)  :: Canopychar

# outputs
#      REAL,DIMENSION(Layers),INTENT(OUT) :: HumidairPa,            
#     &       Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, 
#     &       Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

#output vars
    HumidairPa= np.zeros(Layers)
    Ws        = np.zeros(Layers)
    Sunleaftk = np.zeros(Layers)
    SunleafSH = np.zeros(Layers)
    SunleafLH = np.zeros(Layers)
    SunleafIR = np.zeros(Layers)
    TairK     = np.zeros(Layers) 
    Shadeleaftk= np.zeros(Layers)
    ShadeleafSH= np.zeros(Layers)
    ShadeleafLH= np.zeros(Layers)
    ShadeleafIR= np.zeros(Layers)
    SunStomRes = np.zeros(Layers)
    ShadeStomRes = np.zeros(Layers)

#temporal vars   
    Wsh= np.zeros(Layers)
    Ldepth    = np.zeros(Layers)   
    
    
    
    
         
    Cdepth        = Canopychar[0, Cantype]
    Lwidth        = Canopychar[1, Cantype]
    Llength       = Canopychar[2, Cantype]
    Cheight       = Canopychar[3, Cantype]
    Eps           = Canopychar[9, Cantype]
    TranspireType = Canopychar[10,Cantype]



# 14 = warm (>283K) canopy total humidity change (Pa)
    WarmHumiChange = Canopychar[13,Cantype]
# 15 = cool (>= 283K) canopy total humidity change (Pa)
    CoolHumiChange = Canopychar[14,Cantype]
# 16 = normalized canopy depth where wind is negligible
    NormCanWind = Canopychar[15,Cantype]

    
    if TairK0  > 288:
# Pa m-1  (humidity profile for T < 288)
        Deltah =  WarmHumiChange / Cheight
    elif TairK0  > 278:
        Deltah =(WarmHumiChange - ((288-TairK0)/10) *\
                 (WarmHumiChange - CoolHumiChange))/Cheight
    else:
# Pa m-1  (humidity profile for T <278)
        Deltah = CoolHumiChange / Cheight


    


    Ldepth[:]     = Cdepth * Distgauss[:]
    TairK[:]      = TairK0  + (Trate  * Ldepth[:])      # check this
    HumidairPa[:] = HumidairPa0  + (Deltah * Ldepth[:])

    Wsh[:] = (Cheight-Ldepth[:]) - (NormCanWind * Cheight)
    Wsindex1  = np.where(Wsh < 0.001)
    Wsindex2  = np.where(Wsh >= 0.001)
    Ws[Wsindex1] = 0.05
    Ws[Wsindex2] = (Ws0*np.log(Wsh[Wsindex2])/np.log(Cheight-NormCanWind*Cheight))
    #Ws[:]  = (Ws0*np.log(Wsh[:])/np.log(tmp))
    
    
    
#parameter for IR transfer    
    for i in range(Layers):        
# REVISE - Replace UnexposedLeafIR with LeafIR
#        IRin     = UnexposedLeafIRin(TairK(i), Eps)
#        ShadeleafIR(i) = 2 * IRin
#        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin
# Apparent atmospheric emissivity for clear skies: 
# function of water vapor pressure (Pa) 
# and ambient Temperature (K) based on Brutsaert(1975) 
# referenced in Leuning (1997)

        #EmissAtm        = 0.642 * (HumidairPa[i] / TairK[i])**(1./7.)   
        EmissAtm        = 0.642 * (HumidairPa0 / TairK[i])**(1./7.)               
        IRin            = LeafIR(TairK[i], EmissAtm)
        #IRin            = EmissAtm * Sb * (2 * (TairK[i]**4.)) 
        ShadeleafIR[i]  = IRin
        SunleafIR[i]    = IRin 
        
      # Sun

        List_LeafEB = LeafEB(SunPPFD[i],SunQv[i]+SunQn[i],\
                             SunleafIR[i], Eps, TranspireType, Lwidth, Llength,\
                                 TairK[i], HumidairPa[i], Ws[i]) 
                
        #REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH
        IRout        = List_LeafEB[0]
        Sunleaftk[i] = List_LeafEB[1]
        SunleafSH[i] = List_LeafEB[2]
        SunleafLH[i] = List_LeafEB[3]
        SunStomRes[i]= List_LeafEB[4]
        #Stomatal Conductance
        SunleafIR[i] = SunleafIR[i] - IRout

      # Shade
        List_LeafEB = LeafEB(ShadePPFD[i], ShadeQv[i]+ShadeQn[i],\
                             ShadeleafIR[i],Eps,TranspireType, Lwidth,Llength,\
                            TairK[i], HumidairPa[i], Ws[i])
            
        IRout         = List_LeafEB[0]
        Shadeleaftk[i]= List_LeafEB[1]
        ShadeleafSH[i]= List_LeafEB[2]
        ShadeleafLH[i]= List_LeafEB[3]
        #Stomatal Conductance
        ShadeStomRes[i] = List_LeafEB[4]

        ShadeleafIR[i] = ShadeleafIR[i] - IRout
        #test
        #SunleafIR[i] = IRout
    #end for   
    List_out = [HumidairPa,Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK,\
                Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR, SunStomRes, ShadeStomRes]
    return List_out  

def laiadj(laic,S):
    laiadj = laic / (1 - Canopychar[16,S])
    return laiadj

def CalcwaterVPpa(RH, tk): #input temperature in K
    CalcwaterVPpa = RH * 0.01 * (0.6112 * np.exp((17.67 * (tk - 273.16)) / (tk - 29.66)))
    return CalcwaterVPpa

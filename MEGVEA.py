#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 22:10:24 2020

@author: wangh
Based On Code by Alex Guenther
"""

#response functions:
# 1. gamma_a: EA leaf age response
# 2. gamma_cd: EA canopy depath response
# 3. gamma_laibidir: EA bidirectional exchange LAI response    
# 4. gamma_aq: EA Air Quality response    
# 5. gamma_ht: EA High Temperature response 
# 6. gamma_lt: EA Low Temparature response
# 7. gamma_hw: EA High Wind response
# 8. gamma_co2: CO2 concentration response for isoprene only
# 9. gamtld: EA Temperature response (light dependent emission)
#10. gamp: EA Light response
#11. gamtli: EA Temperature response (light independent emission)
#12. gamhcho: EA HCHO response algorithm by Josh and Josh


import numpy as np

# leaf age response
def gamma_a(LAIp,LAIc,d_temp):
    #LAIp: pass LAI
    #LAIc: current LAI
    #d_temp: 
    #t: time step of LAI (8 for MODIS)
#Parameters for leaf age algorithm for each emission activity classes
    Anew = np.array([0.05, 0.05,  2., 2.,  2., \
                      2.0,  2.0,  2., 0.4, 0.4, \
                      3.5,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0])
    Agro = np.array([ 0.6,  0.6, 1.8, 1.8, 1.8, \
                      1.8,  1.8, 1.8, 0.6, 0.6, \
                      3.0,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0])
        
    Amat = np.array([ 1.0,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0])

    Aold = np.array([ 0.9,  0.9,1.05,1.05,1.05, \
                     1.05, 1.05,1.05,0.95,0.95, \
                      1.2,  1.0, 1.0, 1.0, 1.0, \
                      1.0,  1.0, 1.0, 1.0, 1.0])
   
    
    Tt = d_temp
    t = 8
    
    #****************************** 
    if LAIp < LAIc :
        #Calculate ti and tm
        if Tt <= 303.0:
          ti = 5.0 + 0.7*(300-Tt)
          
        else:
          ti = 2.9

        tm = 2.3*ti
        #get Fnew
        if ti >= t:
            Fnew = 1.0 - (LAIp/LAIc)
        else:
            Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
        #get Fmat
        if tm >= t:
            Fmat = LAIp/LAIc
        else:
            Fmat = (LAIp/LAIc) + ( (t-tm)/t )*( 1-(LAIp/LAIc) )
        #get grow
        Fgro = 1.- Fnew - Fmat
        Fold = 0.
            
    elif LAIp == LAIc:
        Fnew = 0.0
        Fgro = 0.1
        Fmat = 0.8
        Fold = 0.1

    else: #for LAIp > LAIc
        Fnew = 0.0
        Fgro = 0.0
        Fold = ( LAIp-LAIc ) / LAIp
        Fmat = 1-Fold
    #return Fnew,Fgro,Fold,Fmat
    #****************************** 
    #        Calculate GAMLA        
    gamla = Fnew*Anew + Fgro*Agro+Fmat*Amat+ Fold*Aold
    return gamla



#       Emission response to canopy depath
def gamma_cd(Layers,LAI):
    #Layers: number of layers
    #CDEA: output facotr
# canopy depth emission response

    CCD1 = -0.2
    CCD2 = 1.3
    
    
    Cdepth = np.zeros(Layers)
    CDEA = np.zeros(Layers)
    if Layers == 5:
        Cdepth[0]   = 0.0469101
        Cdepth[1]   = 0.2307534
        Cdepth[2]   = 0.5
        Cdepth[3]   = 0.7692465
        Cdepth[4]   = 0.9530899
    else:
        for k in range(Layers):
            kd = k + 1#do not work
            Cdepth[k] = (kd-0.5)/Layers
#********calculate LAI fraction
    for k in range(Layers):

        LAIdepth = LAI*Cdepth[k]
        if LAIdepth > 3:
            LAIdepth = 3.0
        CDEA[k] = CCD1 * LAIdepth + CCD2
    return CDEA
    


# EA bidirectional exchange LAI response    
def gamma_laibidir(LAI):
    if LAI <2:
        gambd = 0.5 * LAI
    elif LAI <= 6 and LAI >= 2:
        gambd = 1 - 0.0625 * (LAI - 2)
    else:
        gambd = 0.75
    return gambd



# EA Air Quality response    
def gamma_aq(S,AQI):
    #S: species Index
    #CAQ:coefficient for poor Air Quality stress
    CAQ = np.array([1., 1.,  1., 5.,  1., \
                    1., 1.,  1., 5.,  5., \
                    1., 1.,  1., 1.,  1.,  \
                    1., 5.,  1., 1.,  1. ])
    #TAQ:threshold for poor Air Quality stress (ppm-hours)
    TAQ = np.array([20., 20.,  20., 20.,  20., \
                    20., 20.,  20., 20.,  20.,\
                    20., 20.,  20., 20.,  20.,\
                    20., 20.,  20., 20.,  20.]) 
    #DTAQ:delta threshold for poor Air Quality stress (ppm-hours)    
    DTAQ = np.array([30., 30.,  30., 30.,  30., \
                    30., 30.,  30., 30.,  30.,\
                    30., 30.,  30., 30.,  30.,\
                    30., 30.,  30., 30.,  30.])

       
    t1 = TAQ[S] + DTAQ[S]
    if AQI <= TAQ[S]:
        gamaq = 1.0
    elif AQI > TAQ[S] and AQI < t1:
        gamaq = 1.0 + (CAQ[S] - 1.0)*(AQI - TAQ[S])/DTAQ[S]
    else:
        gamaq = CAQ[S]
        
    return gamaq



# EA response to high temperature
def gamma_ht(S,MaxT):
    #S: species Index
    #MaxT
    #CHT: coefficient for high temperature stress
    CHT = np.array([1., 1.,  1., 5.,  1., \
                    1., 1.,  1., 5.,  5., \
                    1., 1.,  1., 1.,  1., \
                    1., 5.,  1., 1.,  1. ])
    #THT: threshold for high temperature stress (Celsius degree)
    THT = np.array([40., 40.,  40., 40.,  40.,\
                    40., 40.,  40., 40.,  40.,\
                    40., 40.,  40., 40.,  40.,\
                    40., 40.,  40., 40.,  40.]) 
    #DTHT: delta threshold for high temperature stress (Celsius degree)  
    DTHT = np.array([8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8. ])

    
    THTK = 273.15+THT[S]
    t1 = THTK + DTHT[S]

    if MaxT <= THTK:
        gamht = 1.0
    elif MaxT > THTK and MaxT < t1 :
        gamht = 1.0 + (CHT[S] - 1.0)*(MaxT- THTK)/DTHT[S]
        
    else:
        gamht = CHT[S]
        
    return gamht       
    

    
# EA response to low temperature
def gamma_lt(S,MinT):
    #S: species Index
    #MinT: Minimum Tem
    #CLT: coefficient for low temperature stress
    CLT = np.array([1., 1.,  1., 5.,  1., \
                    1., 1.,  1., 5.,  5., \
                    1., 1.,  1., 1.,  1., \
                    1., 5.,  1., 1.,  1. ])
    #TLT: threshold for low temperature stress (Celsius degree)
    TLT = np.array([10., 10.,  10., 10.,  10.,\
                    10., 10.,  10., 10.,  10.,\
                    10., 10.,  10., 10.,  10.,\
                    10., 10.,  10., 10.,  10.]) 
    #DTLT: delta threshold for low temperature stress (Celsius degree)  
    DTLT = np.array([8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8. ])

    
    TLTK = 273.15+TLT[S]
    t1 = TLTK - DTLT[S]

    if MinT >= TLTK:
        gamlt = 1.0
    elif MinT < TLTK and MinT > t1 :
        gamlt = 1.0 + (CLT[S] - 1.0)*(TLTK-MinT)/DTLT[S]
    else:
        gamlt = CLT[S]
        
    return gamlt       
    
    
    
# EA response to high wind
def gamma_hw(S,MaxWS):
    #S: species Index
    #MaxWS
    #CHW: coefficient for high wind speed stress
    CHW = np.array([1., 1.,  5., 5.,  5., \
                    5., 5.,  5., 5.,  5., \
                    1., 1.,  1., 1.,  1., \
                    1., 5.,  1., 1.,  1. ])
    #THW: threshold for high wind speed stress (m/s)
    THW = np.array([12., 12.,  12., 12.,  12.,\
                    12., 12.,  12., 12.,  12.,\
                    12., 12.,  12., 12.,  12.,\
                    12., 12.,  12., 12.,  12.]) 
    #DTHW: delta threshold for high wind speed stress (m/s)  
    DTHW = np.array([8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8., \
                     8., 8.,  8., 8.,  8. ])

    

    t1 = THW[S] + DTHW[S]

    if MaxWS <= THW[S]:
        gamhw = 1.0
    elif MaxWS > THW[S] and MaxWS < t1 :
        gamhw = 1.0 + (CHW[S] - 1.0)*(MaxWS-THW[S])/DTHW[S]
    else:
        gamhw = CHW[S]
        
    return gamhw      
    


# Subroutine GAMMA_CO2
#CO2 stress on isoprene only
#-----------------------------------------------------------------------
#From Alex Guenther 2017-03-11
def gamma_co2(CO2):
    #CO2: co2 concentration (ppm)
    ISmax = 1.344
    CO2h = 1.4614
    Cstar = 585
    
    CO2temp = CO2
    Ci = 0.7* CO2

    if CO2temp == 400.0:
        gamco2 = 1.0
    else:
        gamco2 = ISmax- ((ISmax*Ci**CO2h) /(Cstar**CO2h+Ci**CO2h))
    return gamco2

     
#gamtld: EA Temperature response (light dependent emission)
def gamtld(T1,T24,T240,S):
    
    #T240=T24
    #Cleo: temperature coefficient (emission type 1: light dependent)
    Cleo = np.array([  2.,  2., 1.83, 1.83, 1.83,\
                     1.83,1.83, 1.83, 2.37, 2.37,\
                     1.60,1.83,   2.,   2., 1.83,\
                     1.83,1.83, 1.83, 1.60, 1.83 ])
    #CT1: temperature coefficient (emission type 1: light dependent)
    Ct1 = np.array([ 95., 95., 80., 80., 80.,\
                     80., 80., 80.,130.,130.,\
                     60., 80., 95., 95., 80.,\
                     80., 80., 80., 60., 80.])        
    Ct2 = 230

    if T1 < 260:
        gamtld = 0
    else:
    #Temperature at which maximum emission occurs
        Topt = 312.5 + 0.6 * (T240 - 297)
        X    = ((1/Topt) - (1/T1))/0.00831
    #Maximum emission (relative to emission at 30 C)
        Eopt = Cleo[S]*np.exp(0.05 * (T24 - 297))*np.exp(0.05*(T240-297))  
        gamtld= Eopt*Ct2*np.exp(Ct1[S]*X)/(Ct2 - Ct1[S]*(1 - np.exp(Ct2 * X)))
    return gamtld
      

#gamp: EA light response
def gamp(PPFD1,PPFD24):
    
      
    C1 = 1.03   
    Alpha  = 0.004
    if PPFD24 < 0.01:
        gamp= 0.
    else:
        
#       C1 = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD)) \
#        * (PPFD24 ** 0.6)
        gamp= (Alpha * C1 * PPFD1)/((1 + Alpha**2. * PPFD1**2.)**0.5)  
    return gamp

#gamp: EA temperature response (light independent emission)
def gamtli(temp,S):
    Ts = 303.15
    #beta: temperature coefficient (emission type 2: light independent)
    beta = np.array([0.13,0.13, 0.10, 0.10, 0.10,\
                     0.10,0.10, 0.10, 0.17, 0.17,\
                     0.08,0.10, 0.13, 0.13, 0.10,\
                     0.10,0.10, 0.10, 0.08, 0.10 ])    
    
    gamtli = np.exp(beta[S]*(temp-Ts) )       
    return gamtli

#gamhcho: EA HCHO response
#*************************************************    
def ResRC(PPFD):
    
    SCadj = ((0.0027*1.066*PPFD)/((1 + 0.0027 * 0.0027 * PPFD**2.)**0.5))
 
    return SCadj

#*************************************************  
#HCHO Resistance (s cm-1)
def HCHOLeafRes(StoCond):
    
    Res= 0.32+50.8*np.exp(-(StoCond-7.5)/18.884)
    return Res
#
def gamhcho(Tem,PPFD,HCHO):
    
    scadj = ResRC(PPFD)
    scadj = scadj*65
    # s cm-1 
    Res = HCHOLeafRes(scadj)
   
    
    #Compensation point of HCHO
    CP = 0.552+0.0368*(Tem -30-273.15)+\
        0.00363*(Tem- 30-273.15)**2+\
        0.000189*(Tem - 30-273.15)**3
    direc = 1
    if HCHO > CP:
        direc = -1
    
    return Res,direc,CP
    
        
      
#*************for drought*********************8
def gamma_sm(SM,WT):
    #SM: Soil moisture at 10 cm deep
    #WT: Wilting Point
    theta = 0.04
    T1 = WT + theta
    if SM > T1:
        gamma_sm = 1.
    elif SM <= T1 and SM > WT:
        gamma_sm = (SM - WT)/theta
    else:
        gamma_sm = 0.
    
    return gamma_sm
def gamma_sm_kc(Kc, Kc_max, Kc_min):
    #Kc: 7 day running averaged Kc
    k1 = -7.45
    b1 = 3.26
    k2 = -28.76
    b2 = 2.35e6
    #gmax = 1.4
    frac = 1.4

    if Kc > Kc_max:
        Kc = Kc_max
    Kc_in = (Kc - Kc_min)/(Kc_max- Kc_min)
    #gsub = 1/(1+b1*np.exp(a1*(Kc-0.2)))
    #glt  = 1/gmax + (1-1/gmax)/(1+b2*np.exp(a2*(1.3-Kc)))
    #gamma_sm = gmax*gsub*glt
    gamma_sm = frac*(1/(1+b1*np.exp(k1*(Kc_in-0.2))))*((1-1/frac)/(1+b2*np.exp(k2*(1.3-Kc_in)))+1/frac)
    return gamma_sm
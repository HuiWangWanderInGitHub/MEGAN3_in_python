#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 19:48:45 2020

@author: wangh
"""


import numpy as np
import MEGVEA as VEA
import TIMEFUNC as TF
import MEGCAN as CAN
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

#import random
#parallel
import numba
from numba import jit,prange


print("MEGANv3.1 in Python")

#Number of species
NEMIS = 19
nmol2mg  = 1e-9*68.12*1000*3600
#Meterology input
Inputfile = "./1.Met_HourlyData_2012_moflux_Kc.csv"
#Inputfile = "./1.Met_HourlyData_2012_moflux.csv"
#output name
Output_name_all = 'Moflux_2012_simulaiton.csv'
Output_name_isop = 'Moflux_2012_simulation_isop.csv'

data1 = np.array(pd.read_csv(Inputfile,header=0))

ParaInputFile = "./2.Parameter.csv"

data2 = np.array(pd.read_csv(ParaInputFile,header=0))

EmisInputFile = "./3.EF_LDF.csv"

data3 = np.array(pd.read_csv(EmisInputFile,header=0))

PFTInputFile = "./4.PFT_Fraction.csv"

data4 = np.array(pd.read_csv(PFTInputFile,header=0))



#=================================================================
#**********MEGAN Control Parameters***********
Para_name= data2[:,0]
Para_value= data2[:,1]

ind= np.where(Para_name == "GAMBD_YN")
para = Para_value[ind]
GAMBD_YN   = para[0] #EA bidirectional exchange LAI response

ind= np.where(Para_name == "GAMAQ_YN")
para = Para_value[ind]
GAMAQ_YN   = para[0] #EA emission activity response to air quality

ind= np.where(Para_name == "GAMHT_YN")
para = Para_value[ind]
GAMHT_YN   = para[0] #EA response to high temperature

ind= np.where(Para_name == "GAMLT_YN")
para = Para_value[ind]
GAMLT_YN   = para[0] #EA response to low temperature

ind= np.where(Para_name == "GAMHW_YN")
para = Para_value[ind]
GAMHW_YN   = para[0] #EA response to high wind speed

ind= np.where(Para_name == "GAMSM_YN")
para = Para_value[ind]
GAMSM_YN   = para[0] #EA response to drought (soil moisture)

GAMHCHO_YN = 0#Para_value[ind] #EA response to HCHO concentration by Josh and Josh

#**********Canopy model inputs**********************
ind= np.where(Para_name == "Latitude")
lat     = Para_value[ind] # 40°30'

#******************RH or QV*************************
ind= np.where(Para_name == "RH_QV")
para = Para_value[ind]
RH_QV   = para # 1 for Relative Humidity; 0: water vapor mixing ratio

#******************Wilting point for the drought algorithm *************************
ind= np.where(Para_name == "WT")
para = Para_value[ind]
WT   = para
#***************************************************
#=================================================================



#=================================================================
NLayers = 5 #Number of Layers
DOY  = data1[:,0]
Hour = data1[:,1]          #hour
TEMP = data1[:,2]+273.15   #Temperature (K)
PPFD = data1[:,4]    #Incoming PPFD (umol/m2/s1)
LAI_array =  data1[:,5]    #LAI
PRES = data1[:,6]             #Pa
WIND = data1[:,7]  #Wind speed  (m s-1)
ISOP_OBS = data1[:,8]
SWC_10CM = data1[:,9]
Kc = data1[:,11]
Kc_max = 0.82
Kc_min = 0.


#*******************************************************
if RH_QV == 1:
    RH   = data1[:,3]        #Relative humidity (%)
else:
    QV   = np.zeros(TEMP.size)
    QV[:]= data1[:,3] # Kg/Kg    
#=================================================================


# Six Canopy Type
# 1  = Needleleaf trees
# 2  = Tropical forest trees, 
# 3  = Temperate broadleaf trees
# 4  = shrubs
# 5  = herbaceous
# 6  = crops 
#CTF  = np.array([0,0,100,0,0,0])
CTF = data4[:,1]

#Stress Indexes
AQI     = 40   #Air Quality Index
CO2     = 403  #CO2 concentration

#=================Read EF and LDF =========================
#emission factor
EF = data3[:,1]
#LDF: light dependent fraction 
LDF = data3[:,2]
#=================================================================



SPC= ["isoprene","MBO","pinenes","ocimenes",\
      "carene","limonene","cymene","camphor",\
      "b-caryophyllene","longifolene","methanol",\
      "acetone","acetaldehyde and ethanol",\
      "formic acid; acetic acid; pyruvic acid",\
      "ethene; ethane","methacrolein","linalool","other VOC","CO"]
       
#***********************************************************
#*******************Calculation Start***********************
#***********************************************************
#output
EMIS = np.zeros([TEMP.size,21])
EMIS[:,:] = np.nan
EMIS[:,0] = data1[:,0]
EMIS[:,1] = data1[:,1]
NRTYP = CTF.size

#Average
NT = TEMP.size
D_PPFD_ARY = np.zeros([TEMP.size])
D_TEMP_1D_ARY = np.zeros([TEMP.size])
D_TEMP_10D_ARY = np.zeros([TEMP.size])
SWC_10D_ARY = np.zeros([TEMP.size])
MaxT_ARY = np.zeros([TEMP.size])
MinT_ARY = np.zeros([TEMP.size])
MaxWS_ARY    = np.zeros([TEMP.size])

# Solar constant [W/m2]
SolarConstant       = 1361.5
# PPB to nmole/m-3; 10^-9*1000/22.4 mole m^-3*10^9 nmole/mole
ppb2nmol = 1000/22.4
SRAD2PPFD = 2.1


TotalCT = np.sum(CTF)*0.01

#start Date
SDAY = DOY[0]


#Time steps
print("========reading inputs successfully=============")
for T in range(TEMP.size) :
    if not np.isnan(PPFD[T]) and not np.isnan(TEMP[T]) and not np.isnan(RH[T]):
        #Average
        Day     = DOY[T]  # Julia Day
        DOY_index = np.where(DOY == Day)
        #D_TEMP = np.average(TEMP[DOY_index])
      
        D_PPFD_ARY[T] = np.nanmean(PPFD[DOY_index])
        D_TEMP_1D_ARY[T] = np.nanmean(TEMP[DOY_index])
        #D_TEMP_10D = D_TEMP_1D
        D_TEMP_10D_ARY[T] = TF.T10(SDAY,Day,DOY,TEMP)  
        SWC_10D_ARY[T] = TF.SWC30D(SDAY,DOY[T],DOY,SWC_10CM)
        MaxT_ARY[T] =   np .max(TEMP[DOY_index])
        MinT_ARY[T] =   np.min(TEMP[DOY_index])
        MaxWS_ARY[T] =  np.max(WIND[DOY_index])      
#calculate the dependant vars

print("========      start running        =============")


#for T in range(5608):

for T in prange(NT):
    
    if not np.isnan(PPFD[T]) and not np.isnan(TEMP[T]) and not np.isnan(RH[T]):
        #print(str(DOY[T])+":"+str(Hour[T]))
#Read input parameters
        if T == 1:
                LAIc    = LAI_array[T]  #current  LAI 
                LAIp    = LAI_array[T] #previous LAI
        else:
                LAIc    = LAI_array[T-1]  #current  LAI 
                LAIp    = LAI_array[T] #previous LAI


        Day     = DOY[T]  # Julia Day

        #Average
        #DOY_index = np.where(DOY == Day)
        #D_TEMP = np.average(TEMP[DOY_index])
      
        #D_PPFD = np.mean(PPFD[DOY_index])
        #D_TEMP_1D = np.mean(TEMP[DOY_index])

        #D_TEMP_10D = TF.T10(SDAY,Day,DOY,TEMP)  
        #SWC_10D = TF.SWC30D(SDAY,DOY[T],DOY,SWC_10CM)
        #MaxT =   np.max(TEMP[DOY_index])
        #MinT =   np.min(TEMP[DOY_index])
        #MaxWS =  np.max(WIND[DOY_index])      


        D_PPFD = D_PPFD_ARY[T]
        D_TEMP_1D = D_TEMP_1D_ARY[T]
        #D_TEMP_10D = D_TEMP_1D
        D_TEMP_10D = D_TEMP_10D_ARY[T]
        SWC_10D = SWC_10D_ARY[T]
        MaxT =   MaxT_ARY[T]
        MinT =   MinT_ARY[T]
        MaxWS =  MaxWS_ARY[T]         
    
        SunleafTK   = np.zeros(NLayers)   
        ShadeleafTK = np.zeros(NLayers) 
        SunPPFD     = np.zeros(NLayers)     
        ShadePPFD   = np.zeros(NLayers)   
        SunFrac     = np.zeros(NLayers)
        #Stomatal Conductance s cm-1
        sunStomRec  = np.zeros(NLayers)
        shadeStomRec= np.zeros(NLayers)  
 
        sun_ppfd_total     = np.zeros(NLayers)
        shade_ppfd_total   = np.zeros(NLayers)
        sun_tk_total       = np.zeros(NLayers)
        shade_tk_total     = np.zeros(NLayers)
        sun_frac_total     = np.zeros(NLayers)
        #Stomatal Conductance s cm-1
        sun_StomRec_total  = np.zeros(NLayers)
        shade_StomRec_total  = np.zeros(NLayers)


        sun_ppfd     = np.zeros(NLayers)
        shade_ppfd   = np.zeros(NLayers)
        sun_tk       = np.zeros(NLayers)
        shade_tk     = np.zeros(NLayers)
        sun_frac     = np.zeros(NLayers)
    
    
#****************Canopy Model***************************************
        #SunleafTK   = TEMP[T]
        #ShadeleafTK = TEMP[T]
        #SunPPFD     = PPFD[T]
        #ShadePPFD   = PPFD[T]
        #SunFrac     = 1.0
    
       
        #Solar angle
        Beta   = TF.Calcbeta(Day , lat[0] , Hour[T] )
        Sinbeta    = np.sin(Beta/ 57.29578)
        TairK0 = TEMP[T]
        Ws0    = WIND[T]
        Solar  = PPFD[T]/SRAD2PPFD
        SM     = SWC_10CM[T]
            
        Maxsolar = Sinbeta  * SolarConstant * TF.CalcEccentricity(Day)
    
        #Gaussian Distribution
        VPgausDis = CAN.GaussianDist(NLayers)
        List_SolarFrac = TF.SolarFractions(Solar, Maxsolar)
        Qdiffv = List_SolarFrac[0]
        Qbeamv = List_SolarFrac[1]
        Qdiffn = List_SolarFrac[2]
        Qbeamn = List_SolarFrac[3]
         
    
        sun_ppfd_total     = 0.
        shade_ppfd_total   = 0.
        sun_tk_total       = 0.
        shade_tk_total     = 0.
        sun_frac_total     = 0.
    
        if TotalCT > 0.0 and LAIc > 0.0:
    
            for i in range(NRTYP):
                #for i in range(NRTYP):
                sun_ppfd     = 0.
                shade_ppfd   = 0.
                sun_tk       = 0.
                shade_tk     = 0.
                sun_frac     = 0.
        
            
            
                          #def CanopyRad(Distgauss, layers, LAI, Sinbeta,Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype):   
                List_CanRad = CAN.CanopyRad(VPgausDis, NLayers, LAIc, Sinbeta,Qbeamv, Qdiffv, Qbeamn, Qdiffn, i)
                sun_frac   = List_CanRad[0]
                QbAbsV     = List_CanRad[1]
                QbAbsn     = List_CanRad[2]
                SunQn      = List_CanRad[3]
                ShadeQn    = List_CanRad[4]
                SunQv      = List_CanRad[5]
                ShadeQv    = List_CanRad[6]
                sun_ppfd   = List_CanRad[7]
                shade_ppfd = List_CanRad[8]
                QdAbsV     = List_CanRad[9]
                QsAbsV     = List_CanRad[10]
                QdAbsn     = List_CanRad[11]
                QsAbsn     = List_CanRad[12]

#**************Relative Humidity or Water Vapor Mixing Ratio****************
                if RH_QV == 1:
                    #HumidairPa0= CAN.CalcwaterVPpa(RH[T], TEMP[T])
                    HumidairPa0= CAN.RHtoWaterVapPres(RH[T], TEMP[T])
                    RH_test = CAN.WaterVapPrestoRH(HumidairPa0/1000, TEMP[T])
                else:
                    HumidairPa0= CAN.WaterVapPres(QV[T], PRES[T])
                    RH_test = CAN.WaterVapPrestoRH(HumidairPa0/1000, TEMP[T])
#***************************************************************************
            #if i == 2:
               # print(HumidairPa0,RH_test)
            
            
                Trate      = CAN.Stability(i, Solar)

            
                List_CanEB =  CAN.CanopyEB(Trate, NLayers, VPgausDis,  i,    TairK0, Ws0, \
                          sun_ppfd,shade_ppfd, SunQv, ShadeQv, SunQn, ShadeQn,HumidairPa0)
   
             
                HumidairPa = List_CanEB[0]
                Ws         = List_CanEB[1]
                sun_tk     = List_CanEB[2]
                SunleafSH  = List_CanEB[3]
                SunleafLH  = List_CanEB[4]
                SunleafIR  = List_CanEB[5]
                TairK      = List_CanEB[6]
                shade_tk   = List_CanEB[7]
                ShadeleafSH= List_CanEB[8]
                ShadeleafLH= List_CanEB[9]
                ShadeleafIR= List_CanEB[10]
                Sun_StomRes = List_CanEB[11]
                Shade_StomRes = List_CanEB[12]
            

                sun_ppfd_total   = sun_ppfd_total + 0.01*CTF[i]*sun_ppfd
                shade_ppfd_total = shade_ppfd_total + 0.01*CTF[i]*shade_ppfd
                sun_tk_total     = sun_tk_total + 0.01*CTF[i]*sun_tk
                shade_tk_total   = shade_tk_total + 0.01*CTF[i]*shade_tk
                sun_frac_total   = sun_frac_total + 0.01*CTF[i]*sun_frac
                sun_StomRec_total = sun_StomRec_total + (0.01*CTF[i]/TotalCT)*Sun_StomRes
                shade_StomRec_total = shade_StomRec_total + 0.01*CTF[i]/TotalCT*Shade_StomRes
            
            
            
            SunleafTK   = sun_tk_total/TotalCT
            ShadeleafTK = shade_tk_total/TotalCT
            SunPPFD     = sun_ppfd_total/TotalCT
            ShadePPFD   = shade_ppfd_total/TotalCT
            SunFrac     = sun_frac_total/TotalCT

        # Leaf age activity factor
        GAMLA  = VEA.gamma_a(LAIp, LAIc, D_TEMP_1D)
        
#****************MEGAN Model***************************************
        #ER  = np.zeros(NEMIS)         
        for S in range(NEMIS): 
        
        #LDF
            LDFMAP = LDF[S]
        
        #EA resonse to soil moilsture
            if S == 0 and GAMSM_YN == 1: 
                #GAMSM  = VEA.gamma_sm(SM,WT)
                GAMSM  = VEA.gamma_sm_kc(Kc[T], Kc_max, Kc_min)
                
            else:
                GAMSM  = 1
        
        # Emission response to canopy depth
            CDEA   = VEA.gamma_cd(NLayers,LAIc)
        
        #EA bidirectional exchange LAI response
            if GAMBD_YN == 1:
                GAMBD = VEA.gamma_laibidir(LAIc)
            else:
                GAMBD = 1.
        
        #emission activity response to air quality
            if GAMAQ_YN == 1:
                GAMAQ = VEA.gamma_aq(S,AQI)
            else:
                GAMAQ = 1.
        
        #EA response to high temperature
            if GAMHT_YN == 1:
                GAMHT = VEA.gamma_ht(S,MaxT)
            else:
                GAMHT = 1.
            
        #EA response to low temperature
            if GAMLT_YN == 1:
                GAMLT = VEA.gamma_lt(S,MinT)
            else:
                GAMLT = 1.
    
        #EA response to high wind speed
            if GAMHW_YN == 1:
                GAMHW = VEA.gamma_hw(S,MaxWS)
            else:
                GAMHW = 1.
            
        #EA response to CO2 only for isoprene
            GAMCO2 = 1.
            if S == 0 and GAMHW_YN == 1:
                GAMCO2 = VEA.gamma_co2(CO2)
                

 
# EA response to canopy temperature/light
            VPGWT = np.zeros(NLayers)
            if NLayers == 5:
                VPGWT[0] = 0.1184635
                VPGWT[1] = 0.2393144
                VPGWT[2] = 0.284444444
                VPGWT[3] = 0.2393144
                VPGWT[4] = 0.1184635
            else:
                VPGWT = 1/NLayers

        
            Ea1L = np.zeros(NLayers)
            Ea2L = np.zeros(NLayers)
            
            
            
            for K in range(NLayers):
                
                    
                Ea1L[K]=CDEA[K]*\
                    (VEA.gamtld(SunleafTK[K],D_TEMP_1D,D_TEMP_10D,S)*\
                     VEA.gamp(SunPPFD[K],D_PPFD)*SunFrac[K]+\
                    VEA.gamtld(ShadeleafTK[K],D_TEMP_1D,D_TEMP_10D,S)*\
                    VEA.gamp(ShadePPFD[K],D_PPFD)*\
                    (1-SunFrac[K]))
           
                
                Ea2L[K]=VEA.gamtli(SunleafTK[K],S)*SunFrac[K]+\
                        VEA.gamtli(ShadeleafTK[K],S)*(1-SunFrac[K])
                
            GAMTP = np.sum((Ea1L*LDFMAP+Ea2L*(1-LDFMAP))*VPGWT)
            #****************TEST*******************
 
                               
            if S == 0:
                ER = LAIc * GAMTP * GAMCO2 * GAMLA[S] * \
                    GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM
                EMIS[T,S+2] = ER*EF[S]
            
#       GAMBD only applied to ethanol and acetaldehyde
            elif S== 12:
                ER = LAIc * GAMTP * GAMBD * GAMLA[S] * \
                    GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM
                EMIS[T,S+2] = ER*EF[S]

            else:
                ER = LAIc * GAMTP * GAMLA[S] *\
                    GAMHW * GAMAQ * GAMHT * GAMLT * GAMSM
                EMIS[T,S+2] = ER*EF[S]




#output Results
item = ["day","hour"]

first_line = item + SPC
head = TF.join_l(first_line,",")

DOY  = data1[:,0]
Hour = data1[:,1]          #hour
TEMP = data1[:,2]+273.15   #Temperature (K)
PPFD = data1[:,4]    #Incoming PPFD (umol/m2/s1)
LAI_array =  data1[:,5]    #LAI
PRES = data1[:,6]             #Pa
WIND = data1[:,7]  #Wind speed  (m s-1)
ISOP_OBS = data1[:,8]
SWC_10CM = data1[:,9]


Isop_head = ["Isop_obs","Isop_mod","Tem(Celcius)","PPFD(umol/m2/s)","LAI","Swc_10cm"]
first_line_isop = item + Isop_head
head_isop = TF.join_l(first_line_isop,",")
fmt_all = '%7.3f'
Output_Isop = np.zeros([TEMP.size,8])
Output_Isop[:,0] = data1[:,0]
Output_Isop[:,1] = data1[:,1]
Output_Isop[:,2] = ISOP_OBS
Output_Isop[:,3] = EMIS[:,2]*nmol2mg
Output_Isop[:,4] = data1[:,2]
Output_Isop[:,5] = data1[:,4]#PPFD
Output_Isop[:,6] = data1[:,7]#LAI
Output_Isop[:,7] = data1[:,9]#SWC10cm

#print(head)
encd =0
print("========        writing outputs    =============")
#np.savetxt("output.csv", EMIS, delimiter=',', header=head, encoding=encd)
np.savetxt(Output_name_all, EMIS, fmt=fmt_all, delimiter=',', header=head) 
np.savetxt(Output_name_isop, Output_Isop, fmt=fmt_all, delimiter=',', header=head_isop) 
#print('encd',encd)

print("========         finished          =============")

isop_obs_ra = ISOP_OBS
isop_mod_ra = EMIS[:,2]*nmol2mg

Time_ins = np.where((Hour >= 9) & (Hour <=17))
Time_all = DOY+Hour
Plot_time = Time_all[Time_ins]
isop_obs_plot = isop_obs_ra[Time_ins]
isop_mod_plot = isop_mod_ra[Time_ins]

isop_obs_p1 = isop_obs_plot[~np.isnan(isop_obs_plot)]
isop_mod_p1 = isop_mod_plot[~np.isnan(isop_obs_plot)]

isop_obs = isop_obs_p1[~np.isnan(isop_mod_p1)]
isop_mod = isop_mod_p1[~np.isnan(isop_mod_p1)]

dif = isop_obs - isop_mod

slope, intercept, r_value, p_value, std_err = stats.linregress(isop_obs,isop_mod)
print(r_value,slope,intercept)
slope_fmt= "%.2f" %slope 
intercept_fmt= "%.2f" %intercept 
r2 = "%0.2f" %(r_value*r_value)
Text1 = "y = "+slope_fmt+"x + "+intercept_fmt
Text2 = "R$^2$ = "+r2
#axs[0].plot(data1[:,0],EMIS[:,2]*nmol2mg,'go',label = "Isoprene_no_drought")
x_d = np.arange(100)

fig, axs = plt.subplots(1,2,figsize=(12, 6))
axs[0].plot(Plot_time,isop_obs_plot,'ro',markersize=3,label = "Observation")
axs[0].plot(Plot_time,isop_mod_plot,'bo',markersize=3,label = "Model")
axs[0].legend()
axs[0].set_xlabel("Time")
axs[0].set_ylim(0,50)
axs[0].set_ylabel("Isoprene Flux(mg m-2 h-1)")
axs[1].plot(isop_obs_plot,isop_mod_plot,'ro',label = "Isoprene")
axs[1].plot(x_d,x_d,'-',)
#axs[1].plot(isop_obs_nonan,isop_mod_nonan,'ro',label = "Isoprene")

axs[1].set_xlabel("Isoprene Obs(mg m-2 h-1)")
axs[1].set_ylabel("Isoprene Sim(mg m-2 h-1)")
axs[1].text(35,47,Text1)
axs[1].text(35,43,Text2)

axs[1].set_xlim(0,50)
axs[1].set_ylim(0,50)


plt.show()











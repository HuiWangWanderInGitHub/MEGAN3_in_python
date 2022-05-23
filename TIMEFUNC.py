#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 19:45:43 2020

@author: wangh
time related functions
"""

import numpy as np
import numba
from numba import jit,prange



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION Calcbeta
#   Calculates the solar zenith angle
#   Code originally developed by Alex Guenther in 1990s
#   Coded into FORTRAN by Xuemei Wang
#   Coded into Python by Hui Wang
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
@jit(nopython=True)
def Calcbeta(Day, Lat, Hour):
      Pi = 3.14159
      Rpi180 = 57.29578
#--------------------------------------------------------------------
      SinDelta = -np.sin(0.40907) * np.cos(6.28 * (Day + 10) / (365))
      CosDelta = (1 - SinDelta**2.)**0.5

      A = np.sin(Lat / Rpi180) * SinDelta
      B = np.cos(Lat / Rpi180) * CosDelta
      Sinbeta = A + B * np.cos(2 * Pi * (Hour - 12) / 24)
      Calcbeta = np.arcsin(Sinbeta) * 57.29578
      return Calcbeta

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   FUNCTION CalcEccentricity
#
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
@jit(nopython=True)
def CalcEccentricity(Day):

#----------------------------------------------------------------
      CalcEccentricity = 1 + 0.033 * np.cos(2*3.14*(Day-10)/365)
    
      return CalcEccentricity
  
    
    
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#   SUBROUTINE SolarFractions
#   Based on actual and potential max solar radiation:
#   Determine the fraction of solar radiation that is 
#   diffuse PPFD, direct PPFD, diffuse near IR, direct near IR 
#
#   Originally developed by Alex Guenther in 1990s
#   Modified in 2010
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
@jit(nopython=True)
def SolarFractions( Solar, Maxsolar):
#-----------------------------------------------------
    if Maxsolar  <= 0:
        Transmis  = 0.5
    elif Maxsolar < Solar:
        Transmis  = 1.
    else:
        Transmis  = Solar/Maxsolar

#FracDiff is based on Lizaso 2005
    FracDiff = 0.156 + 0.86/(1 + np.exp(11.1*(Transmis -0.53)))

#PPFDfrac is based on Goudrian and Laar 1994
    PPFDfrac  = 0.55 -Transmis*0.12

#PPFDdifFrac is based on data in Jacovides 2007
    PPFDdifFrac = FracDiff *(1.06 + Transmis*0.4)  

# Calculate  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
    if PPFDdifFrac > 1.0:
        PPFDdifFrac  = 1.0

    Qv  = PPFDfrac * Solar
    Qdiffv = Qv * PPFDdifFrac
    Qbeamv = Qv - Qdiffv
    Qn = Solar - Qv
    Qdiffn =  Qn * FracDiff
    Qbeamn =  Qn - Qdiffn
    List_out = [Qdiffv,Qbeamv,Qdiffn,Qbeamn]
    return List_out
 

#*********file function***********************
def join_l(l, sep):
    li = iter(l)
    string = str(next(li))
    for i in li:
        string += str(sep) + str(i)
    return string

#*************get 10-day temperature**********
@jit(nopython=True)
def T10(STIME,CTIME,TIME,TEM):
    #STIME: start date
    #CTIME: current date
    #TIME: julian day array
    #TEM: Temperature array 
    DAY_RANGE = CTIME - STIME
    #print(STIME,CTIME)
    if DAY_RANGE >= 10.:
        STIME = CTIME - 10
        TIME_ARR = np.arange(STIME,CTIME)
        TEM_DAILY = np.zeros(TIME_ARR.size)
        for i in range(TIME_ARR.size):
            indexes = np.where(TIME == TIME_ARR[i])
            
            TEM_TMP = TEM[indexes]
            TEM_DAILY[i] = np.nanmean(TEM_TMP)
        
        TEM_OUT = np.nanmean(TEM_DAILY)
        #print("10 days,",TEM_OUT)
    elif DAY_RANGE == 0.:

        indexes = np.where(TIME == CTIME)
        TEM_TMP = TEM[indexes]
        
        TEM_OUT = np.nanmean(TEM_TMP)
        #print("First day,",TEM_OUT)
    else:
        
        TIME_ARR = np.arange(STIME,CTIME)
        TEM_DAILY = np.zeros(TIME_ARR.size)
        for i in range(TIME_ARR.size):
            indexes = np.where(TIME == TIME_ARR[i])
            TEM_TMP = TEM[indexes]
            TEM_DAILY[i] = np.nanmean(TEM_TMP)
            #print(indexes)
        #print(TEM_DAILY)
        TEM_OUT = np.nanmean(TEM_DAILY)
        #print(">10 days,",TEM_OUT)
    return TEM_OUT

def SWC30D(STIME,CTIME,TIME,SWC):
    #STIME: start date
    #CTIME: current date
    #TIME: julian day array
    #SWC: SWC at 10cm array 
    
    DAY_RANGE = CTIME - STIME
    #print(STIME,CTIME)
    if DAY_RANGE >= 30.:
        STIME = CTIME - 30
        TIME_ARR = np.arange(STIME,CTIME)
        SWC_DAILY = np.zeros(TIME_ARR.size)
        for i in range(TIME_ARR.size):
            indexes = np.where(TIME == TIME_ARR[i])
            
            SWC_TMP = SWC[indexes]
            SWC_DAILY[i] = np.nanmean(SWC_TMP)
        
        SWC_OUT = np.nanmean(SWC_DAILY)
        
    elif DAY_RANGE == 0.:

        indexes = np.where(TIME == CTIME)
        SWC_TMP = SWC[indexes]
        
        SWC_OUT = np.nanmean(SWC_TMP)
        #print("First day,",TEM_OUT)
    else:
        
        TIME_ARR = np.arange(STIME,CTIME)
        SWC_DAILY = np.zeros(TIME_ARR.size)
        for i in range(TIME_ARR.size):
            indexes = np.where(TIME == TIME_ARR[i])
            SWC_TMP = SWC[indexes]
            SWC_DAILY[i] = np.nanmean(SWC_TMP)
            #print(indexes)
        #print(TEM_DAILY)
        SWC_OUT = np.nanmean(SWC_DAILY)
        #print(">10 days,",TEM_OUT)
    return SWC_OUT   
    
    
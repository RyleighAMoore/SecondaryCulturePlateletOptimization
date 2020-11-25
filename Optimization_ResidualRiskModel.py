"""
@author: Ryleigh Moore
rmoore@math.utah.edu
"""

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings("ignore")  # Suppress overflow in exponential. We can get exp(-inf) = 0 which is okay for our case.


ts1= 36 # First sample time in hours after blood draw
tE = 120 # time of expiration
LE = 1 # Exposure threshold
Np0 = 1 # Inoculum size
Vp = 300 # Volume of platelet component
Vs = 8 # Volume of sample for culture
h1 = 12 # hold time after first smple (hours)
h2 = 12 # hold time after second smple (hours)
Nmax = Vp*10**9 # Upper limit of bacterial growth
Tmax = 10*10**9 # Time that the upper limit of bacterial growth is reached


# Discretization of ts2, TD, and TL
ts2 =  np.arange(0,169,1)
deltaTD = 1
deltaTL = 1
TDList = np.arange(0.001,40,deltaTD)
TDList = TDList.tolist()
TLList = np.arange(0,170,deltaTL)

    
def harm(t, TL, TD, Vp, LE, Np0):
    Npt = Np(t, TL, TD, Np0)
    if Npt/Vp < LE:
        return 0
    else:
        return 1
    
    
def Np(t, TL, TD, Np0):
    if t < TL:
        return Np0
    if TL <= t and t <= Tmax:
        return Np0*np.exp((0.69/TD)*(t-TL))
    if t > Tmax:
        return Nmax


ts2Constrained = []
TotalRisks = []
TotalRisksNormalized = []
x=ts1 # Take x as ts1 and y as ts2 values
A = tE - x - h1 - h2
for y in ts2:
    countNumRisky = 0
    if (x + h1 < y) and (y + h2 < tE) and (x < y) and (x >= 0) and (y >= 0): # Constraint
        ts2Constrained.append(y) # Track which ts2 values pass constraint
        TotalRisk = 0    
        RiskTDNonzero = True # There are still risks involved
        TDListCopy = TDList.copy()
        while RiskTDNonzero and len(TDListCopy) > 0:
            td = TDListCopy.pop(0)
            totalTDRisk = 0 # Risk for a TD over all TL
            for tl in TLList:
                p2 = (tE-y-h2)/A
                p1 = 1-p2
                pFgivenSts1 = np.exp(-(Vs/Vp)*Np(x,tl,td, Np0))
                assert pFgivenSts1 <=1 and pFgivenSts1>=0
                pFgivenSts2 = np.exp(-(Vs/Vp)*Np(y,tl,td, Np0))
                assert pFgivenSts2 <=1 and pFgivenSts2>=0
                
                harmx = harm(y, tl, td, Vp, LE, Np0)
                harmy = harm(tE, tl, td, Vp, LE, Np0)
                
                Risk = (1-p2)*pFgivenSts1*harmx + p2*pFgivenSts2*harmy
                Risk = Risk*deltaTD*deltaTL
                
                TotalRisk = TotalRisk + Risk
                
                if Risk > 0:
                    countNumRisky = countNumRisky+1

                if tl == np.max(TLList):
                    assert Risk ==0
                    
                totalTDRisk = totalTDRisk + Risk
            
            if totalTDRisk == 0: # No need to look at larger TD values for risk, risk will be 0
                RiskTDNonzero = False
            
            if len(TDListCopy) == 0: 
                if totalTDRisk!=0: # Need to add more TD values, the largest value still has risk
                    TDListCopy.append(td+deltaTD)
                    
        TotalRisks.append(TotalRisk) # Save the total risk of associated with a ts2 value
        TotalRisksNormalized.append(TotalRisk*(1/countNumRisky))

# Find optimal value of ts2, if there are multiple values, the discretization of ts2 or TL and TD could be made finer
TotalRisks = np.asarray(TotalRisks)
TotalRisksNormalized = np.asarray(TotalRisksNormalized)

minIndex = np.asarray(np.where(TotalRisks == TotalRisks.min()))[0]
minIndexNormalized = np.asarray(np.where(TotalRisksNormalized == TotalRisksNormalized.min()))[0]
assert minIndexNormalized == minIndex

for i in range(len(minIndex)): # In case there are multiple minimums due to the discretization being too coarse.
    idx = minIndex[i]
    print('The optimal second testing time is', ts2Constrained[idx], 'hours.')
    print('The total risk is', np.round(TotalRisks[idx],2))
    print('The total risk divided by the number of risky scenarios is', np.round(TotalRisksNormalized[idx],4))
     


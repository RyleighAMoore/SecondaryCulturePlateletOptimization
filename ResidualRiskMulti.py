import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D


Vp = 300
Vs = 8
# LD = 0.99
h1 = 12
h2 = 12
Nmax = Vp*10**9
# tE = 168
Tmax = 10*10**9

# ts1 = np.arange(0,72+1,4)
ts1 = [24, 36]
ts2 = np.arange(0,169,1)
# ts2 = np.append(ts2,np.max(ts2))
deltaTD = 1
deltaTL = 1

TDList = np.arange(0.001,40,deltaTD)
TDList = TDList.tolist()
TLList = np.arange(0,170,deltaTL)

# TLList = np.arange(0,200,1)
# TDList = np.arange(0,20,0.1)
    
LEs = []
Np0s = []
tE2s = []
ts1Vals = []
ts2Stars = []
ts2Vals = []
RiskStars = []
AniX = []
AniY =[]
AnoColor = []
AllAreas = []
AllAreasScaled = []
NumRiskyVals = []

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
        # print(Np0*np.exp((0.69/TD)*(t-TL)))
        return Np0*np.exp((0.69/TD)*(t-TL))
    if t > Tmax:
        return Nmax

for LE in tqdm([1,10**3, 10**5]): #tqdm([1, 10, 10**2, 10**3, 10**4, 10**5])
    for Np0 in [1, 50, 100]:
        for tE in [120, 144, 168]:#range(120,168+1,8): 
            gridX = []
            gridY = []
            AllRisks =[]
            for x in tqdm(ts1):
                xxs = []
                yys = []
                TotalRisks = []
                for y in ts2:
                    if (x + h1 < y) and (y + h2 < tE) and (x < y) and (x >= 0) and (y >= 0):
                        xxs.append(x)
                        yys.append(y)
                        gridX.append(x)
                        gridY.append(y)
                        TotalRisk = 0
                        RiskyTLsx = []
                        RiskyTDsx = []
                        RiskyTLsy = []
                        RiskyTDsy = []
                        NoRiskyTLs = []
                        NoRiskyTDs = []
                        RiskVals= []
                        Alltd = []
                        Alltl = []
                        countNumRisky = 0
                        RowNonzero = True # There are still risks involved
                        TDListCopy = TDList.copy()
                        while RowNonzero and len(TDListCopy)>0:
                            td = TDListCopy.pop(0)
                            rowRisk = 0
                            for tl in TLList:
                                Alltd.append(td)
                                Alltl.append(tl)
                                A = tE - x - h1 - h2
                                p2 = (tE-y-h2)/A
                                p1 = 1-p2
                                pFgivenSts1 = np.exp(-(Vs/Vp)*Np(x,tl,td, Np0))
                                assert pFgivenSts1 <=1 and pFgivenSts1>=0
                                pFgivenSts2 = np.exp(-(Vs/Vp)*Np(y,tl,td, Np0))
                                assert pFgivenSts2 <=1 and pFgivenSts2>=0
                                
                                harmx = harm(y, tl, td, Vp, LE, Np0)
                                harmy = harm(tE, tl, td, Vp, LE, Np0)
                                
                                if harmx == 1:
                                    RiskyTLsx.append(tl)
                                    RiskyTDsx.append(td)
                                if harmy == 1:
                                    RiskyTLsy.append(tl)
                                    RiskyTDsy.append(td)
                                else:
                                    NoRiskyTLs.append(tl)
                                    NoRiskyTDs.append(td)
                                
                                Risk = (1-p2)*pFgivenSts1*harmx + p2*pFgivenSts2*harmy
                                Risk = Risk*deltaTD*deltaTL
                                
                                if Risk > 0:
                                    countNumRisky = countNumRisky+1
                                
                                TotalRisk = TotalRisk + Risk

                                
                                if tl == np.max(TLList):
                                    assert Risk == 0
                                    
                                RiskVals.append(Risk)
                                rowRisk = rowRisk + Risk #Risk over a row in TL vs TD space
                            
                            # print(rowRisk)
                            if rowRisk == 0:
                                RowNonzero = False
                            
                            if len(TDListCopy) == 0:
                                if rowRisk!=0:
                                    TDListCopy.append(td+deltaTD)
                                    
                        assert(TotalRisk>0)
                        
                        LEs.append(LE)
                        Np0s.append(Np0)
                        tE2s.append(tE)
                        ts1Vals.append(x)
                        ts2Vals.append(y)
                        AllAreas.append(TotalRisk)
                        AllAreasScaled.append(TotalRisk*1/countNumRisky)
                        NumRiskyVals.append(countNumRisky)
                        
                        

 

import pandas as pd
df = pd.DataFrame()
df['LE'] = LEs
df['Np0'] = Np0s
df['tE'] = tE2s
df['ts1'] = ts1Vals
df['ts2'] = ts2Vals
df['Risk'] = AllAreas
df['Risk Scaled'] = AllAreasScaled
df['Number Risky Values'] = NumRiskyVals
df.isnull().values.any()
df.to_excel("dataForBobScaledMoreData.xlsx") 
                   
                    
        
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.scatter(np.asarray(gridX), np.asarray(gridY), np.asarray(AllRisks),  c='k', marker='o')
# plt.xlabel("ts1")
# plt.ylabel("ts2")
# plt.show()

# fig = plt.figure()
# plt.scatter(ts1Vals, ts2Stars)
# plt.xlabel("ts1")
# plt.ylabel("ts2*")
# plt.show()

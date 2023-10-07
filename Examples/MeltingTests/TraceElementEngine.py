import numpy as np
import pandas as pd
import pyMelt
from pyMelt import chemistry as c

D = {}
D['ol'] = c.olv_D
D['cpx'] = c.cpx_D
D['opx'] = c.opx_D
D['g'] = c.grt_D
D['spl'] = c.spn_D
# D['Other'] = D['ol'].copy()
# for e in D['Other'].keys():
#     D['Other'][e] = 0.0
Ds = pd.DataFrame(columns = list(D['ol'].keys()), index = list(D.keys()), data = np.zeros((len(D.keys()), len(D['ol'].keys()))))
for d in D:
    Ds.loc[d] = D[d]

Co = c.workman05_dmm

def accum_frac_melt(Mass, source = None):
    if source is None:
        source = pd.DataFrame(data = Co, index = [0])

    F_inst = Mass['liq'].copy()
    for i in range(len(F_inst)):
        if i != 0:
            F_inst.loc[i] = Mass['liq'].loc[i] - Mass['liq'].loc[i-1]
    
    Mass['liq'] = F_inst
    Mass_sum = Mass.copy()
    Mass_sum = Mass_sum.drop(columns = 'liq')
    Tot = Mass_sum.sum(axis = 1)
    for g in Mass_sum.keys():
        Mass_sum[g] = Mass_sum[g]/Tot

    # calculate bulk D
    Kd = pd.DataFrame(columns = list(D['ol'].keys()), data = np.zeros((len(F_inst),len(list(D['ol'].keys())))))
    C_residual = Kd.copy()
    C_melt = Kd.copy()
    for i in range(len(Kd['La'])):
        for p in Mass_sum:
            if p != "liq":
                if p in D.keys():
                    Kd.loc[i] = Kd.loc[i] + Ds.loc[p]*Mass_sum.loc[i,p]

    for i in range(len(C_melt['La'])):
        if i == 0:
            if F_inst[i] > 0.0:
                C_melt.loc[i, list(C_melt.keys())] = (source.loc[0, list(C_melt.keys())]*Mass.sum(axis = 1)[i])/(Kd.loc[i, list(C_melt.keys())]*(Mass.sum(axis = 1)[i] - F_inst[i])+F_inst[i])
                C_residual.loc[i, list(C_melt.keys())] = C_melt.loc[i, list(C_melt.keys())]*Kd.loc[i, list(C_melt.keys())]
            else:
                C_residual.loc[i, list(C_melt.keys())] = source.loc[0, list(C_melt.keys())]
        else:
            if F_inst[i] > 0.0:
                C_melt.loc[i, list(C_melt.keys())] = (C_residual.loc[i-1, list(C_melt.keys())]*Mass.sum(axis = 1)[i])/(Kd.loc[i, list(C_melt.keys())]*(Mass.sum(axis = 1)[i] - F_inst[i])+F_inst[i])
                C_residual.loc[i, list(C_melt.keys())] = C_melt.loc[i, list(C_melt.keys())]*Kd.loc[i, list(C_melt.keys())]
            else:
                C_residual.loc[i, list(C_melt.keys())] = C_residual.loc[i-1, list(C_melt.keys())]

    C_accumulated_melt = C_melt.copy()
    for i in range(len(C_melt['La'])):
        A = (Mass['liq'].loc[0:i]/np.nansum(F_inst[0:i]))
        A = A.values.reshape(i+1,1)
        C_accumulated_melt.loc[i, list(C_melt.keys())] = (C_melt.loc[0:i]*A).sum()
        
    return C_melt
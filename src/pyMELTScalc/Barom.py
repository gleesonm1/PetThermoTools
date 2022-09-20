import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Thermobar as pt
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.MELTS import *
try:
    from pyMELTScalc.Holland import *
except:
    pass

def SatPress(P, Model, bulk, T_initial = None, Phases = None, Fe3 = None, H2O = None, fO2 = None, dt = None, T_step = None, T_cut = None, Plot = None, findRange = None, findMin = None, cores = None):
    '''
    Find the location in P (+ H2O + fO2) where there is the minimum misfit between the saturation curves of the listed phases.
    '''

    if T_initial is None:
        T_initial = 1200

    if cores is None:
        cores = 4

    if Phases is None:
        Phases = ['quartz1', 'plagioclase1', 'k-feldspar1']

    if dt is None:
        dt = 25

    if T_step is None:
        T_step= 1

    if T_cut is None:
        T_cut = 10.0001
    else:
        T_cut = T_cut + 0.0001

    if Fe3 is None and H2O is None:
        Results = {}

        if Plot is not None:
            f = plt.figure(figsize = plt.figaspect(0.5))
            a0 = f.add_subplot(1,2,1)
            a1 = f.add_subplot(1,2,2)

        if len(Phases) == 3:
            a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)
        else:
            a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)

        if Plot is not None:
            a0.plot(P[np.where(a_Sat>0)], a_Sat[np.where(a_Sat>0)], '-k', label = 'a_sat')
            a0.plot(P[np.where(b_Sat>0)], b_Sat[np.where(b_Sat>0)], '-r', label = 'b_sat')
            if len(Phases) == 3:
                a0.plot(P[np.where(c_Sat>0)], c_Sat[np.where(c_Sat>0)], '-b', label = 'c_sat')

            a0.legend()

        # calc residuals
        if len(Phases) == 3:
            Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
        else:
            Res = findResiduals(a_Sat, b_Sat)

        if Plot is not None:
            if len(Phases) == 3:
                a1.plot(P, Res['abc'], '-k', label = 'abc')
                a1.plot(P, Res['ab'], '-r', label = 'ab')
                a1.plot(P, Res['ac'], '-b', label = 'ac')
                a1.plot(P, Res['bc'], '-y', label = 'bc')
            else:
                a1.plot(P, Res['ab'], '-r', label = 'ab')

            a1.legend()

        Results = Res.copy()
        Results['H2O_sat'] = H2O_Melt
        Results['T_Liq'] = T_Liq

        # return residual results in range or minimum not specified
        if findRange is None and findMin is None:
            if Plot is None:
                return Results
            else:
                return Results, f

        if findRange is not None:
            Results_new = detRange(P, Model, bulk, Results, Phases, T_initial = T_initial, T_cut = T_cut, cores = cores)
            Results['Range'] = Results_new

            if Plot is None:
                return Results
            else:
                return Results, f

        if findMin is not None:
            PolyMin = findMinimum(P, Results, Phases, T_cut = T_cut)
            Results['PolyMin'] = PolyMin

            if Plot is None:
                return Results
            else:
                return Results, f

    if H2O is None and Fe3 is not None:
        if Plot is not None:
            f = plt.figure(figsize = plt.figaspect(0.5))
            a0 = f.add_subplot(1,2,1, projection = '3d')
            a1 = f.add_subplot(1,2,2, projection = '3d')

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(Fe3), len(P)))
            Res_ac_mat = np.zeros((len(Fe3), len(P)))
            Res_bc_mat = np.zeros((len(Fe3), len(P)))
            c_Sat_mat = np.zeros((len(Fe3), len(P)))

        Res_ab_mat = np.zeros((len(Fe3), len(P)))
        H2O_mat = np.zeros((len(Fe3), len(P)))
        T_Liq_mat = np.zeros((len(Fe3), len(P)))
        a_Sat_mat = np.zeros((len(Fe3), len(P)))
        b_Sat_mat = np.zeros((len(Fe3), len(P)))

        for i in range(len(Fe3)):
            Bulk_new = bulk.copy()
            Bulk_new[3] = Fe3[i]*((159.69/2)/71.844)*bulk[5]
            Bulk_new[5] = (1-Fe3[i])*bulk[5]

            if len(Phases) == 3:
                a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)
            else:
                a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)

            # calc residuals
            if len(Phases) == 3:
                Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                Res_abc_mat[i, :] = Res['abc']
                Res_ac_mat[i, :] = Res['ac']
                Res_bc_mat[i, :] = Res['bc']
                c_Sat_mat[i, :] = c_Sat
            else:
                Res = findResiduals(a_Sat, b_Sat)

            T_Liq_mat[i, :] = T_Liq
            H2O_mat[i, :] = H2O_Melt
            a_Sat_mat[i, :] = a_Sat
            b_Sat_mat[i, :] = b_Sat
            Res_ab_mat[i, :] = Res['ab']

        if len(Phases) == 3:
            Results = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}
        else:
            Results = {'ab': Res_ab_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}

        if Plot is not None:
            X, Y = np.meshgrid(P, Fe3)
            a0.plot_surface(X[np.where(a_Sat_mat > 0)], Y[np.where(a_Sat_mat > 0)], a_Sat[np.where(a_Sat_mat > 0)], color = 'k', label = 'a_Sat')
            a0.plot_surface(X[np.where(b_Sat_mat > 0)], Y[np.where(b_Sat_mat > 0)], b_Sat[np.where(b_Sat_mat > 0)], color = 'r', label = 'b_Sat')
            if len(Phases) == 3:
                a0.plot_surface(X[np.where(c_Sat_mat > 0)], Y[np.where(c_Sat_mat > 0)], c_Sat[np.where(c_Sat_mat > 0)], color = 'r', label = 'b_Sat')

            a0.legend()

            if len(Phases) == 3:
                a1.plot_surface(X[np.where(~np.isnan(Res['abc']))],Y[np.where(~np.isnan(Res['abc']))], Res['abc'][np.where(~np.isnan(Res['abc']))], color = 'k', label = 'abc')
                a1.plot_surface(X[np.where(~np.isnan(Res['ab']))],Y[np.where(~np.isnan(Res['ab']))], Res['ab'][np.where(~np.isnan(Res['ab']))], color = 'r', label = 'ab')
                a1.plot_surface(X[np.where(~np.isnan(Res['ac']))],Y[np.where(~np.isnan(Res['ac']))], Res['ac'][np.where(~np.isnan(Res['ac']))], color = 'b', label = 'ac')
                a1.plot_surface(X[np.where(~np.isnan(Res['bc']))],Y[np.where(~np.isnan(Res['bc']))], Res['bc'][np.where(~np.isnan(Res['bc']))], color = 'y', label = 'bc')
            else:
                a1.plot_surface(X[np.where(~np.isnan(Res['ab']))],Y[np.where(~np.isnan(Res['ab']))], Res['ab'][np.where(~np.isnan(Res['ab']))], color = 'r', label = 'ab')

            a1.legend()

        # return residual results in range or minimum not specified
        if findRange is None and findMin is None:
            if Plot is None:
                return Results
            else:
                return Results, f

        if findRange is not None:
            Results_new = detRange(P, Model, bulk, Results, Phases, T_initial = T_initial, Fe3 = Fe3, T_cut = T_cut, cores = cores)
            Results['Range'] = Results_new

            if Plot is None:
                return Results
            else:
                return Results, f

    if H2O is not None and Fe3 is None:
        if Plot is not None:
            f = plt.figure(figsize = plt.figaspect(0.5))
            a0 = f.add_subplot(1,2,1, projection = '3d')
            a1 = f.add_subplot(1,2,2, projection = '3d')

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(H2O), len(P)))
            Res_ac_mat = np.zeros((len(H2O), len(P)))
            Res_bc_mat = np.zeros((len(H2O), len(P)))
            c_Sat_mat = np.zeros((len(H2O), len(P)))

        Res_ab_mat = np.zeros((len(H2O), len(P)))
        H2O_mat = np.zeros((len(H2O), len(P)))
        T_Liq_mat = np.zeros((len(H2O), len(P)))
        a_Sat_mat = np.zeros((len(H2O), len(P)))
        b_Sat_mat = np.zeros((len(H2O), len(P)))


        for i in range(len(H2O)):
            Bulk_new = bulk.copy()
            Bulk_new[14] = H2O[i]

            if len(Phases) == 3:
                a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)
            else:
                a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)

            # calc residuals
            if len(Phases) == 3:
                Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                Res_abc_mat[i, :] = Res['abc']
                Res_ac_mat[i, :] = Res['ac']
                Res_bc_mat[i, :] = Res['bc']
                c_Sat_mat[i, :] = c_Sat
            else:
                Res = findResiduals(a_Sat, b_Sat)

            Res_ab_mat[i, :] = Res['ab']
            T_Liq_mat[i, :] = T_Liq
            H2O_mat[i, :] = H2O_Melt
            a_Sat_mat[i, :] = a_Sat
            b_Sat_mat[i, :] = b_Sat

        if len(Phases) == 3:
            Results = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}
        else:
            Results = {'ab': Res_ab_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}

        if Plot is not None:
            X, Y = np.meshgrid(P, H2O)
            a0.plot_surface(X[np.where(a_Sat_mat > 0)], Y[np.where(a_Sat_mat > 0)], a_Sat[np.where(a_Sat_mat > 0)], color = 'k', label = 'a_Sat')
            a0.plot_surface(X[np.where(b_Sat_mat > 0)], Y[np.where(b_Sat_mat > 0)], b_Sat[np.where(b_Sat_mat > 0)], color = 'r', label = 'b_Sat')
            if len(Phases) == 3:
                a0.plot_surface(X[np.where(c_Sat_mat > 0)], Y[np.where(c_Sat_mat > 0)], c_Sat[np.where(c_Sat_mat > 0)], color = 'r', label = 'b_Sat')

            a0.legend()

            if len(Phases) == 3:
                a1.plot_surface(X[np.where(~np.isnan(Res['abc']))],Y[np.where(~np.isnan(Res['abc']))], Res['abc'][np.where(~np.isnan(Res['abc']))], color = 'k', label = 'abc')
                a1.plot_surface(X[np.where(~np.isnan(Res['ab']))],Y[np.where(~np.isnan(Res['ab']))], Res['ab'][np.where(~np.isnan(Res['ab']))], color = 'r', label = 'ab')
                a1.plot_surface(X[np.where(~np.isnan(Res['ac']))],Y[np.where(~np.isnan(Res['ac']))], Res['ac'][np.where(~np.isnan(Res['ac']))], color = 'b', label = 'ac')
                a1.plot_surface(X[np.where(~np.isnan(Res['bc']))],Y[np.where(~np.isnan(Res['bc']))], Res['bc'][np.where(~np.isnan(Res['bc']))], color = 'y', label = 'bc')
            else:
                a1.plot_surface(X[np.where(~np.isnan(Res['ab']))],Y[np.where(~np.isnan(Res['ab']))], Res['ab'][np.where(~np.isnan(Res['ab']))], color = 'r', label = 'ab')

            a1.legend()

        # return residual results in range or minimum not specified
        if findRange is None and findMin is None:
            if Plot is None:
                return Results
            else:
                return Results, f

        if findRange is not None:
            Results_new = detRange(P, Model, bulk, Results, Phases, T_initial = T_initial, H2O = H2O, T_cut = T_cut, cores = cores)
            Results['Range'] = Results_new

            if Plot is None:
                return Results
            else:
                return Results, f

    if H2O is not None and Fe3 is not None:
        if Plot is not None:
            if len(Phases) == 2:
                f = plt.figure(figsize = plt.figaspect(0.33))
                a0 = f.add_subplot(1,3,1, projection = '3d')
                a1 = f.add_subplot(1,3,2, projection = '3d')
                a2 = f.add_subplot(1,3,3, projection = '3d')
            else:
                f = plt.figure(figsize = plt.figaspect(0.33))
                a0 = f.add_subplot(1,3,1, projection = '3d')
                a1 = f.add_subplot(1,3,2, projection = '3d')
                a2 = f.add_subplot(1,3,3, projection = '3d')
                a3 = f.add_subplot(2,3,4, projection = '3d')
                a4 = f.add_subplot(2,3,5, projection = '3d')
                a5 = f.add_subplot(2,3,6, projection = '3d')

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(H2O), len(Fe3), len(P)))
            Res_ac_mat = np.zeros((len(H2O), len(Fe3), len(P)))
            Res_bc_mat = np.zeros((len(H2O), len(Fe3), len(P)))

        Res_ab_mat = np.zeros((len(H2O), len(Fe3), len(P)))
        T_Liq_3Dmat = np.zeros((len(H2O), len(Fe3), len(P)))
        H2O_Liq_3Dmat = np.zeros((len(H2O), len(Fe3), len(P)))

        for i in range(len(H2O)):
            Bulk_new = bulk.copy()
            Bulk_new[14] = H2O[i]

            if len(Phases) == 3:
                Res_abc_mat_2D = np.zeros((len(Fe3), len(P)))
                Res_ac_mat_2D = np.zeros((len(Fe3), len(P)))
                Res_bc_mat_2D = np.zeros((len(Fe3), len(P)))

            Res_ab_mat_2D = np.zeros((len(Fe3), len(P)))
            T_Liq_mat_2D = np.zeros((len(Fe3), len(P)))
            H2O_Liq_mat_2D = np.zeros((len(Fe3), len(P)))

            for j in range(len(Fe3)):

                Bulk_new[3] = Fe3[j]*((159.69/2)/71.844)*bulk[5]
                Bulk_new[5] = (1-Fe3[j])*bulk[5]

                if len(Phases) == 3:
                    a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)
                else:
                    a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = dt, T_step = T_step, cores = cores)

                # calc residuals
                if len(Phases) == 3:
                    Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                    Res_abc_mat_2D[j, :] = Res['abc']
                    Res_ac_mat_2D[j, :] = Res['ac']
                    Res_bc_mat_2D[j, :] = Res['bc']
                else:
                    Res = findResiduals(a_Sat, b_Sat)

                Res_ab_mat_2D[j, :] = Res['ab']
                T_Liq_mat_2D[j, :] = T_Liq
                H2O_Liq_mat_2D[j, :] = H2O_Melt

            if len(Phases) == 3:
                Res_abc_mat[i,:,:] = Res_abc_mat_2D
                Res_ac_mat[i,:,:] = Res_ac_mat_2D
                Res_bc_mat[i,:,:] = Res_bc_mat_2D

            Res_ab_mat[i,:,:] = Res_ab_mat_2D
            T_Liq_3Dmat[i,:,:] = T_Liq_mat_2D
            H2O_Liq_3Dmat[i,:,:] = H2O_Liq_mat_2D

        if len(Phases) == 3:
            Results = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_Liq_3Dmat, 'T_Liq': T_Liq_3Dmat}
        else:
            Results = {'ab': Res_ab_mat, 'H2O_Sat': H2O_Liq_3Dmat, 'T_Liq': T_Liq_3Dmat}

        if findRange is None and findMin is None:
            if Plot is None:
                return Results
            else:
                return Results, f

        if findRange is not None:
            Results_new = detRange(P, Model, bulk, Results, Phases, T_initial = T_initial, Fe3 = Fe3, H2O = H2O, T_cut = T_cut, findRange = findRange, cores = cores)
            Results['Range'] = Results_new

            try:
                print('Max P = ' + str(Results['Range']['abc_max_P']))
                print('Min P = ' + str(Results['Range']['abc_min_P']))
            except:
                print('No Convergence')

            if Plot is None:
                return Results
            else:
                return Results, f

def findMinimum(P, Results, Phases, Fe3 = None, H2O = None, T_cut = None):
    '''
    Take the results of SatPress and search for the minimum point using a parabolic fit.
    '''
    if T_cut is None:
        T_cut = 5

    if H2O is None and Fe3 is None:
        Res_key = list(Results.keys())[:-2]
        PolyMin = {}
        for Res in Res_key:
            if len(Results[Res][~np.isnan(Results[Res])]) > 2:
                p, poly = polymin(P, Results[Res])
                if poly.size == 0:
                    poly = np.array([np.nan])

                elif p[0]*poly[0]**2 + p[1]*poly[0] + p[2] > T_cut:
                    poly = np.array([np.nan])
            else:
                poly = np.array([np.nan])

            PolyMin[Res] = poly

        return PolyMin

def detRange(P, Model, bulk, Results, Phases, T_initial = None, Fe3 = None, H2O = None, T_cut = None, findRange = None, cores = None):
    '''
    re-run the SatPress calculations using a higher resolution to determine the range of P (+H2O +fO2) where the maximum difference between the target saturation curve is below some reference value.
    '''
    if T_cut is None:
        T_cut = 5

    if T_initial is None:
        T_initial = 1200

    Res_key = list(Results.keys())[:-2]

    if Fe3 is None and H2O is None:
        P_min = np.zeros(len(Res_key))
        P_max = np.zeros(len(Res_key))
        i = 0
        for Res in Res_key:
            if len(Results[Res][np.where(Results[Res]<T_cut*2)]) > 0:
                P_min[i] = np.nanmin(P[np.where(Results[Res]<T_cut*2)])
                P_max[i] = np.nanmax(P[np.where(Results[Res]<T_cut*2)])
            else:
                P_min[i] = np.nan
                P_max[i] = np.nan

            i = i + 1

        P_start = np.nanmin(P_min)
        P_end = np.nanmax(P_max)

        P_space = 1 + 3*(P_end - P_start)/((np.max(P) - np.min(P))/(len(P)-1))

        P = np.linspace(P_start, P_end, int(P_space))

        with open('bulk.obj', 'wb') as f:
            pickle.dump(bulk, f)

        if len(Phases) == 3:
            a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)
        else:
            a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)

        # calc residuals
        if len(Phases) == 3:
            Results_new = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
        else:
            Results_new = findResiduals(a_Sat, b_Sat)

        Results_new['H2O_sat'] = H2O_Melt
        Results_new['T_Liq'] = T_Liq

        for Res in Res_key:
            if len(Results_new[Res][np.where(~np.isnan(Results_new[Res]))]) > 0:
                Results_new[Res + '_min_P'] = np.nanmin(P[~np.isnan(Results_new[Res])])
                Results_new[Res + '_max_P'] = np.nanmin(P[~np.isnan(Results_new[Res])])
            else:
                Results_new[Res + '_min_P'] = np.nan
                Results_new[Res + '_max_P'] = np.nan

        return Results_new

    if Fe3 is not None and H2O is None:
        P_min = np.zeros(len(Res_key))
        P_max = np.zeros(len(Res_key))

        Fe3_min = np.zeros(len(Res_key))
        Fe3_max = np.zeros(len(Res_key))

        P_mat, Fe3_mat = np.meshgrid(P, Fe3)

        i = 0
        for Res in Res_key:
            if len(Results[Res][np.where(Results[Res]<T_cut*2)]) > 0:
                P_min[i] = np.nanmin(P_mat[np.where(Results[Res]<T_cut*2)])
                P_max[i] = np.nanmax(P_mat[np.where(Results[Res]<T_cut*2)])

                Fe3_min[i] = np.nanmin(Fe3_mat[np.where(Results[Res]<T_cut*2)])
                Fe3_max[i] = np.nanmax(Fe3_mat[np.where(Results[Res]<T_cut*2)])
            else:

                P_min[i] = np.nan
                P_max[i] = np.nan

                Fe3_min[i] = np.nan
                Fe3_max[i] = np.nan

            i = i + 1

        P_start = np.nanmin(P_min)
        P_end = np.nanmax(P_max)

        Fe3_start = np.nanmin(Fe3_min)
        Fe3_end = np.nanmax(Fe3_max)

        P = np.linspace(P_start, P_end, int(1 + 3*(P_end - P_start)/((np.max(P) - np.min(P))/(len(P)-1))))
        Fe3 = np.linspace(Fe3_start, Fe3_end, int(1 + 3*(Fe3_end - Fe3_start)/((np.max(Fe3) - np.min(Fe3))/(len(Fe3)-1))))

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(Fe3), len(P)))
            Res_ac_mat = np.zeros((len(Fe3), len(P)))
            Res_bc_mat = np.zeros((len(Fe3), len(P)))

        Res_ab_mat = np.zeros((len(Fe3), len(P)))
        H2O_mat = np.zeros((len(Fe3), len(P)))
        T_Liq_mat = np.zeros((len(Fe3), len(P)))

        for i in range(len(Fe3)):
            Bulk_new = bulk.copy()
            Bulk_new[3] = Fe3[i]*((159.69/2)/71.844)*bulk[5]
            Bulk_new[5] = (1-Fe3[i])*bulk[5]

            with open('bulk.obj', 'wb') as f:
                pickle.dump(Bulk_new, f)

            if len(Phases) == 3:
                a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)
            else:
                a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)

            # calc residuals
            if len(Phases) == 3:
                Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                Res_abc_mat[i, :] = Res['abc']
                Res_ac_mat[i, :] = Res['ac']
                Res_bc_mat[i, :] = Res['bc']
            else:
                Res = findResiduals(a_Sat, b_Sat)

            T_Liq_mat[i, :] = T_Liq
            H2O_mat[i, :] = H2O_Melt
            Res_ab_mat[i, :] = Res['ab']

        if len(Phases) == 3:
            Results_new = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}
        else:
            Results_new = {'ab': Res_ab_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}

        P_mat, Fe3_mat = np.meshgrid(P, Fe3)
        for Res in Res_key:
            if len(Results_new[Res][np.where(~np.isnan(Results_new[Res]))]) > 0:
                Results_new[Res + '_min_P'] = np.nanmin(P_mat[np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_P'] = np.nanmin(P_mat[np.where(~np.isnan(Results_new[Res]))])

                Results_new[Res + '_min_Fe3'] = np.nanmin(Fe3_mat[np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_Fe3'] = np.nanmin(Fe3_mat[np.where(~np.isnan(Results_new[Res]))])
            else:
                Results_new[Res + '_min_P'] = np.nan
                Results_new[Res + '_max_P'] = np.nan

                Results_new[Res + '_min_Fe3'] = np.nan
                Results_new[Res + '_max_Fe3'] = np.nan

        return Results_new

    if Fe3 is None and H2O is not None:
        P_min = np.zeros(len(Res_key))
        P_max = np.zeros(len(Res_key))

        H2O_min = np.zeros(len(Res_key))
        H2O_max = np.zeros(len(Res_key))

        P_mat, H2O_mat = np.meshgrid(P, H2O)

        i = 0
        for Res in Res_key:
            if len(Results[Res][np.where(Results[Res]<T_cut*2)]) > 0:
                P_min[i] = np.nanmin(P_mat[np.where(Results[Res]<T_cut*2)])
                P_max[i] = np.nanmax(P_mat[np.where(Results[Res]<T_cut*2)])

                H2O_min[i] = np.nanmin(H2O_mat[np.where(Results[Res]<T_cut*2)])
                H2O_max[i] = np.nanmax(H2O_mat[np.where(Results[Res]<T_cut*2)])
            else:

                P_min[i] = np.nan
                P_max[i] = np.nan

                H2O_min[i] = np.nan
                H2O_max[i] = np.nan

            i = i + 1

        P_start = np.nanmin(P_min)
        P_end = np.nanmax(P_max)

        H2O_start = np.nanmin(H2O_min)
        H2O_end = np.nanmax(H2O_max)

        P = np.linspace(P_start, P_end, int(1 + 3*(P_end - P_start)/((np.max(P) - np.min(P))/(len(P)-1))))
        H2O = np.linspace(H2O_start, H2O_end, int(1 + 3*(H2O_end - H2O_start)/((np.max(H2O) - np.min(H2O))/(len(H2O)-1))))

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(H2O), len(P)))
            Res_ac_mat = np.zeros((len(H2O), len(P)))
            Res_bc_mat = np.zeros((len(H2O), len(P)))

        Res_ab_mat = np.zeros((len(H2O), len(P)))
        H2O_mat = np.zeros((len(H2O), len(P)))
        T_Liq_mat = np.zeros((len(H2O), len(P)))

        for i in range(len(H2O)):
            Bulk_new = bulk.copy()
            Bulk_new[14] = H2O[i]

            with open('bulk.obj', 'wb') as f:
                pickle.dump(Bulk_new, f)

            if len(Phases) == 3:
                a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)
            else:
                a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, bulk, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)

            # calc residuals
            if len(Phases) == 3:
                Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                Res_abc_mat[i, :] = Res['abc']
                Res_ac_mat[i, :] = Res['ac']
                Res_bc_mat[i, :] = Res['bc']
            else:
                Res = findResiduals(a_Sat, b_Sat)

            T_Liq_mat[i, :] = T_Liq
            H2O_mat[i, :] = H2O_Melt
            Res_ab_mat[i, :] = Res['ab']

        if len(Phases) == 3:
            Results_new = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}
        else:
            Results_new = {'ab': Res_ab_mat, 'H2O_Sat': H2O_mat, 'T_Liq': T_Liq_mat}

        P_mat, H2O_mat = np.meshgrid(P, H2O)
        for Res in Res_key:
            if len(Results_new[Res][np.where(~np.isnan(Results_new[Res]))]) > 0:
                Results_new[Res + '_min_P']= np.nanmin(P_mat[np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_P'] = np.nanmin(P_mat[np.where(~np.isnan(Results_new[Res]))])

                Results_new[Res + '_min_H2O'] = np.nanmin(Results_new['H2O_Sat'][np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_H2O'] = np.nanmin(Results_new['H2O_Sat'][np.where(~np.isnan(Results_new[Res]))])
            else:
                Results_new[Res + '_min_P']= np.nan
                Results_new[Res + '_max_P'] = np.nan

                Results_new[Res + '_min_H2O'] = np.nan
                Results_new[Res + '_max_H2O'] = np.nan

        return Results_new

    if H2O is not None and Fe3 is not None:
        Results_new = {}
        P_min = np.zeros(len(Res_key))
        P_max = np.zeros(len(Res_key))

        Fe3_min = np.zeros(len(Res_key))
        Fe3_max = np.zeros(len(Res_key))

        H2O_min = np.zeros(len(Res_key))
        H2O_max = np.zeros(len(Res_key))

        Fe3_mat, H2O_mat, P_mat = np.meshgrid(Fe3, H2O, P)

        i = 0
        for Res in Res_key:
            if len(Results[Res][np.where(Results[Res]<T_cut*2)]) > 0:
                P_min[i] = np.nanmin(P_mat[np.where(Results[Res]<T_cut*2)])
                P_max[i] = np.nanmax(P_mat[np.where(Results[Res]<T_cut*2)])

                H2O_min[i] = np.nanmin(H2O_mat[np.where(Results[Res]<T_cut*2)])
                H2O_max[i] = np.nanmax(H2O_mat[np.where(Results[Res]<T_cut*2)])

                Fe3_min[i] = np.nanmin(Fe3_mat[np.where(Results[Res]<T_cut*2)])
                Fe3_max[i] = np.nanmax(Fe3_mat[np.where(Results[Res]<T_cut*2)])
            else:

                P_min[i] = np.nan
                P_max[i] = np.nan

                H2O_min[i] = np.nan
                H2O_max[i] = np.nan

                Fe3_min[i] = np.nan
                Fe3_max[i] = np.nan

            i = i + 1

        if findRange == 'abc':
            P_start = P_min[0] - (P[1] - P[0])/3
            P_end = P_max[0] + (P[1] - P[0])/3

            H2O_start = H2O_min[0] - (H2O[1] - H2O[0])/3
            H2O_end = H2O_max[0] + (H2O[1] - H2O[0])/3

            Fe3_start = Fe3_min[0] - (Fe3[1] - Fe3[0])/3
            Fe3_end = Fe3_max[0] + (Fe3[1] - Fe3[0])/3
        else:
            P_start = np.nanmin(P_min)
            P_end = np.nanmax(P_max)

            H2O_start = np.nanmin(H2O_min)
            H2O_end = np.nanmax(H2O_max)

            Fe3_start = np.nanmin(Fe3_min)
            Fe3_end = np.nanmax(Fe3_max)

        if P_end - P_start > 0:
            P = np.linspace(P_start, P_end, int(1 + 3*(P_end - P_start)/((np.max(P) - np.min(P))/(len(P)))))
        else:
            return Results_new

        if H2O_end - H2O_start > 0:
            H2O = np.linspace(H2O_start, H2O_end, int(1 + 3*(H2O_end - H2O_start)/((np.max(H2O) - np.min(H2O))/(len(H2O)))))
        else:
            return Results_new

        if Fe3_end - Fe3_start > 0:
            Fe3 = np.linspace(Fe3_start, Fe3_end, int(1 + 3*(Fe3_end - Fe3_start)/((np.max(Fe3) - np.min(Fe3))/(len(Fe3)))))
        else:
            return Results_new

        if len(Phases) == 3:
            Res_abc_mat = np.zeros((len(H2O), len(Fe3), len(P)))
            Res_ac_mat = np.zeros((len(H2O), len(Fe3), len(P)))
            Res_bc_mat = np.zeros((len(H2O), len(Fe3), len(P)))

        Res_ab_mat = np.zeros((len(H2O), len(Fe3), len(P)))
        T_Liq_3Dmat = np.zeros((len(H2O), len(Fe3), len(P)))
        H2O_Liq_3Dmat = np.zeros((len(H2O), len(Fe3), len(P)))

        for i in range(len(H2O)):
            Bulk_new = bulk.copy()
            Bulk_new[14] = H2O[i]

            if len(Phases) == 3:
                Res_abc_mat_2D = np.zeros((len(Fe3), len(P)))
                Res_ac_mat_2D = np.zeros((len(Fe3), len(P)))
                Res_bc_mat_2D = np.zeros((len(Fe3), len(P)))

            Res_ab_mat_2D = np.zeros((len(Fe3), len(P)))
            T_Liq_mat_2D = np.zeros((len(Fe3), len(P)))
            H2O_Liq_mat_2D = np.zeros((len(Fe3), len(P)))

            for j in range(len(Fe3)):

                Bulk_new[3] = Fe3[j]*((159.69/2)/71.844)*bulk[5]
                Bulk_new[5] = (1-Fe3[j])*bulk[5]

                with open('bulk.obj', 'wb') as f:
                    pickle.dump(Bulk_new, f)

                if len(Phases) == 3:
                    a_Sat, b_Sat, c_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)
                else:
                    a_Sat, b_Sat, T_Liq, H2O_Melt = findSat(P, Model, Phases, Bulk_new, T_initial = T_initial, dt = T_cut, T_step = T_cut, cores = cores)

                # calc residuals
                if len(Phases) == 3:
                    Res = findResiduals(a_Sat, b_Sat, c_Sat = c_Sat)
                    Res_abc_mat_2D[j, :] = Res['abc']
                    Res_ac_mat_2D[j, :] = Res['ac']
                    Res_bc_mat_2D[j, :] = Res['bc']
                else:
                    Res = findResiduals(a_Sat, b_Sat)

                Res_ab_mat_2D[j, :] = Res['ab']
                T_Liq_mat_2D[j, :] = T_Liq
                H2O_Liq_mat_2D[j, :] = H2O_Melt

            if len(Phases) == 3:
                Res_abc_mat[i,:,:] = Res_abc_mat_2D
                Res_ac_mat[i,:,:] = Res_ac_mat_2D
                Res_bc_mat[i,:,:] = Res_bc_mat_2D

            Res_ab_mat[i,:,:] = Res_ab_mat_2D
            T_Liq_3Dmat[i,:,:] = T_Liq_mat_2D
            H2O_Liq_3Dmat[i,:,:] = H2O_Liq_mat_2D

        if len(Phases) == 3:
            Results_new = {'abc': Res_abc_mat, 'ab': Res_ab_mat, 'ac': Res_ac_mat, 'bc': Res_bc_mat, 'H2O_Sat': H2O_Liq_3Dmat, 'T_Liq': T_Liq_3Dmat}
        else:
            Results_new = {'ab': Res_ab_mat, 'H2O_Sat': H2O_Liq_3Dmat, 'T_Liq': T_Liq_3Dmat}


        Fe3_mat, H2O_mat, P_mat = np.meshgrid(Fe3, H2O, P)

        for Res in Res_key:
            if len(Results_new[Res][np.where(~np.isnan(Results_new[Res]))]) > 0:
                Results_new[Res + '_min_P'] = np.nanmin(P_mat[np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_P'] = np.nanmax(P_mat[np.where(~np.isnan(Results_new[Res]))])

                Results_new[Res + '_min_H2O'] = np.nanmin(Results_new['H2O_Sat'][np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_H2O'] = np.nanmax(Results_new['H2O_Sat'][np.where(~np.isnan(Results_new[Res]))])

                Results_new[Res + '_min_Fe3'] = np.nanmin(Fe3_mat[np.where(~np.isnan(Results_new[Res]))])
                Results_new[Res + '_max_Fe3'] = np.nanmax(Fe3_mat[np.where(~np.isnan(Results_new[Res]))])
            else:
                Results_new[Res + '_min_P'] = np.nan
                Results_new[Res + '_max_P'] = np.nan

                Results_new[Res + '_min_H2O'] = np.nan
                Results_new[Res + '_max_H2O'] = np.nan

                Results_new[Res + '_min_Fe3'] = np.nan
                Results_new[Res + '_max_Fe3'] = np.nan

        return Results_new


def findResiduals(a_Sat, b_Sat, c_Sat = None):
    '''
    Determine the temperature offse between the different saturation curves.
    '''
    Res = {}

    if c_Sat is not None:
        Res_abc = np.nanmax(np.array([abs(a_Sat-b_Sat), abs(a_Sat - c_Sat), abs(b_Sat - c_Sat)]), axis = 0)
        Res_abc[np.where(a_Sat == 0)] = np.nan
        Res_abc[np.where(b_Sat == 0)] = np.nan
        Res_abc[np.where(c_Sat == 0)] = np.nan

        Res_ac = abs(a_Sat - c_Sat)
        Res_ac[np.where(a_Sat == 0)] = np.nan
        Res_ac[np.where(c_Sat == 0)] = np.nan

        Res_bc = abs(b_Sat - c_Sat)
        Res_bc[np.where(b_Sat == 0)] = np.nan
        Res_bc[np.where(c_Sat == 0)] = np.nan

        Res['abc'] = Res_abc
        Res['ac'] = Res_ac
        Res['bc'] = Res_bc

    Res_ab = abs(a_Sat - b_Sat)
    Res_ab[np.where(a_Sat == 0)] = np.nan
    Res_ab[np.where(b_Sat == 0)] = np.nan

    Res['ab'] = Res_ab

    return Res

def polymin(P, Res):
    '''
    Finds the minimum residual temperature using a 2nd degree polynomial.
    '''
    arr = np.sort(Res)
    Ind = np.where(Res == arr[0])[0][0]

    if P[Ind] == np.nanmax(P):
        p = np.array([0,0,0])
        p_min = np.array([P[Ind]])
    elif P[Ind] == np.nanmin(P):
        p = np.array([0,0,0])
        p_min = np.array([P[Ind]])
    else:
        p = np.polyfit(P[np.array([Ind-1, Ind, Ind+1])],Res[np.array([Ind-1, Ind, Ind+1])],2)

        x = np.linspace(np.nanmin(P),np.nanmax(P),501)
        y = p[0]*x**2 + p[1]*x + p[2]

        p_min = x[np.where(y == np.nanmin(y))]

    return p, p_min

def findSat(P, Model, Phases, bulk, T_initial = None, dt = None, T_step = None, cores = None):
    '''
    Calculation to find saturation temperatures. Can be run in parallel.
    '''
    if dt is None:
        dt = 25

    if T_step is None:
        T_step = 1

    if T_initial is None:
        T_initial = 1200

    if cores is None:
        cores = 20

    a_Sat = np.zeros(len(P))
    b_Sat = np.zeros(len(P))
    T_Liq = np.zeros(len(P))
    H2O_Melt = np.zeros(len(P))
    P_out = np.zeros(len(P))

    if len(Phases) == 3:
        c_Sat = np.zeros(len(P))

    A = len(P)//cores
    B = len(P) % cores

    Group = np.zeros(A) + cores
    Group = np.append(Group, B)

    qs = []
    q = Queue()

    for j in range(len(Group)):
        ps = []
        if Model == "MELTS":
            for i in range(int(cores*j), int(cores*j + Group[j])):
                p = Process(target = satTemperature_MELTS, args = (q, Model, Phases, P[i], T_initial, bulk, dt, T_step))
                ps.append(p)
                p.start()

            for p in ps:
                try:
                    ret = q.get(timeout = 30)
                except:
                    ret = []

                qs.append(ret)

            TIMEOUT = 120
            start = time.time()
            for p in ps:
                if p.is_alive():
                    while time.time() - start <= TIMEOUT:
                        if not p.is_alive():
                            p.join()
                            p.terminate()
                            break
                        time.sleep(.1)
                    else:
                        p.terminate()
                        p.join(5)
                else:
                    p.join()
                    p.terminate()

            time.sleep(.1)

    for i in range(len(qs)):
        if len(qs[i])>0:
            if len(Phases) == 3:
                a_Sat[i], b_Sat[i], c_Sat[i], T_Liq[i], H2O_Melt[i], P_out[i] = qs[i]
            else:
                a_Sat[i], b_Sat[i], T_Liq[i], H2O_Melt[i], P_out[i] = qs[i]

    a_Sat_new = np.zeros(len(a_Sat))
    b_Sat_new = np.zeros(len(a_Sat))
    T_Liq_new = np.zeros(len(a_Sat))
    H2O_Melt_new = np.zeros(len(a_Sat))
    if len(Phases) == 3:
        c_Sat_new = np.zeros(len(a_Sat))

    for i in range(len(P)):
        if len(np.where(P_out == P[i])[0]) > 0:
            a_Sat_new[i] = a_Sat[np.where(P_out == P[i])]
            b_Sat_new[i] = b_Sat[np.where(P_out == P[i])]
            T_Liq_new[i] = T_Liq[np.where(P_out == P[i])]
            H2O_Melt_new[i] = H2O_Melt[np.where(P_out == P[i])]
            if len(Phases) == 3:
                c_Sat_new[i] = c_Sat[np.where(P_out == P[i])]

    if len(Phases) == 3:
        return a_Sat_new, b_Sat_new, c_Sat_new, T_Liq_new, H2O_Melt_new
    else:
        return a_Sat_new, b_Sat_new, T_Liq_new, H2O_Melt_new
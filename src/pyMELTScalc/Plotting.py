import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint, Point, Polygon
import matplotlib.pyplot as plt
from pyMELTScalc.GenFuncs import *

def harker(Results = None, x_axis = None, y_axis = None, phase = None, line_style = None, line_color = None, data = None, d_color = None, d_marker = None):
    '''
    Construct harker plots.

    Parameters:
    -----------
    Results: dict
        Contains DataFrames with the results of the MELTS or MAGEMin calculations.

    x_axis: str
        Oxide to be placed on the x-axis, Default = "MgO".

    y_axis: list
        Oxides to be displayed on the y-axes.

    phase: str
        Phase compositions to be plotted. Plots the liquid component by default.

    line_style: str
        Line style to use for the MELTS/MAGEMin results. Default = '-'.

    line_color: str or tuple
        Color of the line used to display the MELTS/MAGEMin results. Black as default.

    data: DataFrame
        Optional. Include natural or experimental data to plot against the calculation results.

    d_color: str or tuple
        Color of the symbols used to display the natural/experimental data.

    d_marker: str
        Marker style for the natural/experimental data.
    '''
    if Results is None:
        raise Exception("We need some results to work with!")

    if x_axis is None:
        x_axis = "MgO"

    if y_axis is None:
        y_axis = ["SiO2", "TiO2", "Al2O3", "FeOt", "CaO", "Na2O"]

    if phase is None:
        phase = "liquid1"

    if line_style is None:
        line_style = '-'

    if line_color is None:
        line_color = 'k'

    if data is not None:
        if d_color is None:
            d_color = 'red'

        if d_marker is None:
            d_marker = 'o'

    if type(y_axis) == list:
        y_axis = np.array(y_axis)

    if len(y_axis) % 3 == 0.0:
        y_axis = y_axis.reshape(len(y_axis)//3, 3)
    else:
        y_axis = np.append(y_axis, np.array(["None"] * (3 - (len(y_axis) % 3))))
        y_axis = y_axis.reshape(len(y_axis)//3, 3)

    f, a = plt.subplots(np.shape(y_axis)[0], np.shape(y_axis)[1], figsize = (3.5 * np.shape(y_axis)[1], 3.5 * np.shape(y_axis)[0]))
    if type(phase) == str:
        for i in range(np.shape(y_axis)[0]):
            for j in range(np.shape(y_axis)[1]):
                if y_axis[i,j] != "None":
                    if data is not None:
                        a[i][j].plot(data.loc[:,data.columns.str.contains(x_axis)], data.loc[:,data.columns.str.contains(y_axis[i,j])], d_marker, markerfacecolor = d_color, markeredgecolor = 'k', markersize = 4)

                    a[i][j].plot(Results['All'][x_axis + Names[phase]], Results['All'][y_axis[i,j] + Names[phase]], line_style, linewidth = 2, color = line_color)
                    a[i][j].set_ylabel(y_axis[i][j] + " wt%")
                    if i != np.shape(y_axis)[0] - 1:
                        if i == np.shape(y_axis)[0] - 2:
                            if y_axis[i+1,j] == "None":
                                a[i][j].set_xlabel(x_axis + " wt%")
                    else:
                        a[i][j].set_xlabel(x_axis + " wt%")

                else:
                    a[i][j].axis('off')

    f.tight_layout()

def residualT_plot(Results = None, P_bar = None, phases = None, H2O_Liq = None, Fe3Fet_Liq = None, T_cut_C = None):
    if T_cut_C is None:
        T_cut_C = 12

    if H2O_Liq is None and Fe3Fet_Liq is None:
        if len(Results['CurveMin']) == 4:
            f, a = plt.subplots(2,2, figsize = (10,8), sharex = True, sharey = True)
            f.tight_layout()
            a[1][0].set_xlabel('P (bars)')
            a[1][1].set_xlabel('P (bars)')
            a[0][0].set_ylabel('Residual T ($\degree$C)')
            a[1][0].set_ylabel('Residual T ($\degree$C)')

            m = np.array([['Res_abc', 'Res_ab'], ['Res_ac', 'Res_bc']])
            Name = np.array([['Three phase saturation', phases[0] + ' - ' + phases[1]],
                [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])

            for i in range(2):
                for j in range(2):
                    a[i][j].set_title(Name[i,j])
                    if ~np.isnan(Results['CurveMin'][m[i,j]]['P_min']):
                        a[i][j].plot(P_bar, Results[m[i,j]][0,0,:], 'ok', markerfacecolor="b", label="original", markersize = 8)
                        a[i][j].plot(Results['CurveMin'][m[i,j]]['P_new'], Results['CurveMin'][m[i,j]]['y_new'],
                                    '-', c="r", label="spline fit")
                        a[i][j].plot([np.nanmin(Results['CurveMin'][m[i,j]]['P_new']), np.nanmax(Results['CurveMin'][m[i,j]]['P_new'])],
                                    [Results['CurveMin'][m[i,j]]['Res_min'], Results['CurveMin'][m[i,j]]['Res_min']], ':k')
                        a[i][j].plot([Results['CurveMin'][m[i,j]]['P_min'], Results['CurveMin'][m[i,j]]['P_min']],
                                    [np.nanmin(Results['CurveMin'][m[i,j]]['y_new']) - 5,
                                    np.nanmax(Results['CurveMin'][m[i,j]]['y_new']) + 5], ':k')
        else:
            f, a = plt.subplots(1,1, figsize = (5,4))
            a.set_xlabel('P (bars)')
            a.set_ylabel('Residual T ($\degree$C)')
            a.set_title(phases[0] + ' - ' + phases[1])
            if ~np.isnan(Results['CurveMin']['Res_ab']['P_min']):
                a.plot(P_bar, Results['Res_ab'][0,0,:], 'ok', markerfacecolor="b", label="original", markersize = 8)
                a.plot(Results['CurveMin']['Res_ab']['P_new'], Results['CurveMin']['Res_ab']['y_new'],
                            '-', c="r", label="spline fit")
                a.plot([np.nanmin(Results['CurveMin']['Res_ab']['P_new']), np.nanmax(Results['CurveMin']['Res_ab']['P_new'])],
                            [Results['CurveMin']['Res_ab']['Res_min'], Results['CurveMin']['Res_ab']['Res_min']], ':k')
                a.plot([Results['CurveMin']['Res_ab']['P_min'], Results['CurveMin']['Res_ab']['P_min']],
                            [np.nanmin(Results['CurveMin']['Res_ab']['y_new']) - 5,
                            np.nanmax(Results['CurveMin']['Res_ab']['y_new']) + 5], ':k')

    if H2O_Liq is not None and Fe3Fet_Liq is None:
        if len(Results['CurveMin']) == 4:
            f = plt.figure(figsize = (10,8))
            a1 = f.add_subplot(2,2,1, projection = '3d')
            a2 = f.add_subplot(2,2,2, projection = '3d')
            a3 = f.add_subplot(2,2,3, projection = '3d')
            a4 = f.add_subplot(2,2,4, projection = '3d')

            a = np.array([[a1,a2], [a3, a4]])

            m = np.array([['Res_abc', 'Res_ab'], ['Res_ac', 'Res_bc']])
            Name = np.array([['Three phase saturation', phases[0] + ' - ' + phases[1]],
                [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])

            X, Y = np.meshgrid(P_bar, H2O_Liq)
            Y = Results['H2O_melt'][:,0,:].copy()

            for i in range(2):
                for j in range(2):
                    a[i][j].set_xlabel('P (bars)')
                    a[i][j].set_ylabel('H$_{2}$O (wt%)')
                    a[i][j].set_zlabel('Residual T ($\degree$C)')
                    a[i][j].set_title(Name[i,j])

            for i in range(2):
                for j in range(2):
                    if ~np.isnan(Results['CurveMin'][m[i,j]]['P_min']):
                        A = Results[m[i][j]][:,0,:].copy()
                        X, Y = np.meshgrid(P_bar, H2O_Liq)
                        Y = Results['H2O_melt'][:,0,:].copy()

                        Y_save = Y.copy()

                        a[i][j].scatter(X, Y, A+1, marker = 'o', facecolor = 'red')
                        a[i][j].scatter(Results['CurveMin'][m[i,j]]['P_min'], Results['CurveMin'][m[i,j]]['H2O_min'], Results['CurveMin'][m[i,j]]['Res_min'], marker = '^', facecolor = 'yellow')

                        A[np.where(A > T_cut_C*2)] = np.nan

                        H2O_new = np.linspace(np.nanmin(Y[np.where(~np.isnan(A))]), np.nanmax(Y[np.where(~np.isnan(A))]), 200)
                        P_new = np.linspace(np.nanmin(X[np.where(~np.isnan(A))]), np.nanmax(X[np.where(~np.isnan(A))]), 200)

                        X_new, Y_new = np.meshgrid(P_new, H2O_new)
                        x = X[~np.isnan(A)].flatten()
                        y = Y[~np.isnan(A)].flatten()

                        MyPoly = MultiPoint(list(zip(x, y))).convex_hull

                        points = list(zip(X_new.flatten(), Y_new.flatten()))
                        Include = np.zeros(len(X_new.flatten()))
                        for k in range(len(points)):
                            p = Point(points[k])
                            Include[k] = p.within(MyPoly)

                        YayNay = Include.reshape(X_new.shape)
                        x_new = X_new[np.where(YayNay == True)].flatten()
                        y_new = Y_new[np.where(YayNay == True)].flatten()

                        z_plot = X_new.copy() * 0
                        z_plot = z_plot.flatten()
                        z_plot[np.where(Include == True)] = Results['CurveMin'][m[i][j]]['z_new'](x_new, y_new, grid = False)
                        z_plot = z_plot.reshape(X_new.shape)
                        X_new[np.where(YayNay == False)] = np.nan
                        Y_new[np.where(YayNay == False)] = np.nan

                        a[i][j].plot_surface(X_new, Y_new, z_plot, cmap = 'viridis')
                        a[i][j].set_zlim([0,50])

    if H2O_Liq is not None and Fe3Fet_Liq is not None:
        if len(Results['CurveMin']) == 4:
            f = plt.figure(figsize = (10,8))
            a1 = f.add_subplot(2,2,1, projection = '3d')
            a2 = f.add_subplot(2,2,2, projection = '3d')
            a3 = f.add_subplot(2,2,3, projection = '3d')
            a4 = f.add_subplot(2,2,4, projection = '3d')

            a = np.array([[a1,a2], [a3, a4]])

            m = np.array([['Res_abc', 'Res_ab'], ['Res_ac', 'Res_bc']])
            Name = np.array([['Three phase saturation', phases[0] + ' - ' + phases[1]],
                [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])

            X, Y = np.meshgrid(P_bar, H2O_Liq)
            Y = Results['H2O_melt'][:,0,:].copy()

            for i in range(2):
                for j in range(2):
                    a[i][j].set_xlabel('P (bars)')
                    a[i][j].set_ylabel('H$_{2}$O (wt%)')
                    a[i][j].set_zlabel('Residual T ($\degree$C)')
                    a[i][j].set_title(Name[i,j])

            for i in range(2):
                for j in range(2):
                    if ~np.isnan(Results['CurveMin'][m[i,j]]['P_min']):
                        loc = np.where(Fe3Fet_Liq == Results['CurveMin'][m[i,j]]['Fe3Fet_Liq'])[0][0]
                        A = Results[m[i][j]][:,loc,:].copy()
                        A[np.where(A > T_cut_C)] = np.nan
                        X, Y = np.meshgrid(P_bar, H2O_Liq)
                        Y = Results['H2O_melt'][:,loc,:].copy()

                        a[i][j].scatter(X, Y, A+1, marker = 'o', facecolor = 'red')

                        H2O_new = np.linspace(np.nanmin(Y[np.where(~np.isnan(A))]), np.nanmax(Y[np.where(~np.isnan(A))]), 200)
                        P_new = np.linspace(np.nanmin(X[np.where(~np.isnan(A))]), np.nanmax(X[np.where(~np.isnan(A))]), 200)

                        X_new, Y_new = np.meshgrid(P_new, H2O_new)
                        x = X[~np.isnan(A)].flatten()
                        y = Y[~np.isnan(A)].flatten()

                        MyPoly = MultiPoint(list(zip(x, y))).convex_hull

                        points = list(zip(X_new.flatten(), Y_new.flatten()))
                        Include = np.zeros(len(X_new.flatten()))
                        for k in range(len(points)):
                            p = Point(points[k])
                            Include[k] = p.within(MyPoly)

                        YayNay = Include.reshape(X_new.shape)
                        x_new = X_new[np.where(YayNay == True)].flatten()
                        y_new = Y_new[np.where(YayNay == True)].flatten()

                        z_plot = X_new.copy() * 0
                        z_plot = z_plot.flatten()
                        z_plot[np.where(Include == True)] = Results['CurveMin'][m[i][j]]['z_new'](x_new, y_new, grid = False)
                        z_plot = z_plot.reshape(X_new.shape)
                        X_new[np.where(YayNay == False)] = np.nan
                        Y_new[np.where(YayNay == False)] = np.nan

                        a[i][j].plot_surface(X_new, Y_new, z_plot, cmap = 'viridis')
                        a[i][j].set_zlim([0,50])







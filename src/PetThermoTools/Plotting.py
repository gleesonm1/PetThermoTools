import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint, Point, Polygon
import matplotlib.pyplot as plt
from PetThermoTools.GenFuncs import *

def harker(Results = None, x_axis = None, y_axis = None, phase = None, line_style = None, line_color = None, data = None, d_color = None, d_marker = None, label = None):
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
    if 'All' in list(Results.keys()):
        if type(phase) == str:
            for i in range(np.shape(y_axis)[0]):
                for j in range(np.shape(y_axis)[1]):
                    if y_axis[i,j] != "None":
                        if data is not None:
                            a[i][j].plot(data.loc[:,data.columns.str.contains(x_axis)].values, data.loc[:,data.columns.str.contains(y_axis[i,j])].values, d_marker, markerfacecolor = d_color, markeredgecolor = 'k', markersize = 4)

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

    else:
        if type(phase) == str:
            for i in range(np.shape(y_axis)[0]):
                for j in range(np.shape(y_axis)[1]):
                    if y_axis[i,j] != "None":
                        if data is not None:
                            a[i][j].plot(data.loc[:,data.columns.str.contains(x_axis)], data.loc[:,data.columns.str.contains(y_axis[i,j])], d_marker, markerfacecolor = d_color, markeredgecolor = 'k', markersize = 4, label = "Data")

                        a[i][j].set_ylabel(y_axis[i][j] + " wt%")
                        if i != np.shape(y_axis)[0] - 1:
                            if i == np.shape(y_axis)[0] - 2:
                                if y_axis[i+1,j] == "None":
                                    a[i][j].set_xlabel(x_axis + " wt%")
                        else:
                            a[i][j].set_xlabel(x_axis + " wt%")

                    else:
                        a[i][j].axis('off')

        for r in Results:
            Res = Results[r].copy()
            if type(phase) == str:
                for i in range(np.shape(y_axis)[0]):
                    for j in range(np.shape(y_axis)[1]):
                        if y_axis[i,j] != "None":
                            a[i][j].plot(Res['All'][x_axis + Names[phase]], Res['All'][y_axis[i,j] + Names[phase]], line_style, linewidth = 2, label = r)

        if label is not None:
            a[0][0].legend()

    f.tight_layout()

def plot_surfaces(Results = None, P_bar = None, phases = None, H2O_Liq = None):
    if H2O_Liq is None:
        f, a = plt.subplots(1,1, figsize = (5,4))
        a.set_xlabel('P (bars)')
        a.set_ylabel('T ($\degree$C)')
        for i in range(len(phases)):
            if i == 0:
                a.plot(P_bar, Results[phases[0]][0,0,:], '-r', linewidth = 2, label = phases[i])
            if i == 1:
                a.plot(P_bar, Results[phases[1]][0,0,:], '-b', linewidth = 2, label = phases[i])
            if i == 2:
                a.plot(P_bar, Results[phases[2]][0,0,:], '-k', linewidth = 2, label = phases[i])

        a.legend()

    if H2O_Liq is not None:
        X, Y = np.meshgrid(P_bar, H2O_Liq)
        Y = Results['H2O_melt'][:,0,:].copy()

        f = plt.figure(figsize = (5,4))
        a = f.add_subplot(1,1,1, projection = '3d')
        for i in range(len(phases)):
            if i == 0:
                a.plot_surface(X, Y, Results[phases[0]][:,0,:], color = 'r', label = phases[i])
            if i == 1:
                a.plot_surface(X, Y, Results[phases[1]][:,0,:], color = 'b', label = phases[i])
            if i == 2:
                a.plot_surface(X, Y, Results[phases[2]][:,0,:], color = 'k', label = phases[i])

    return f, a

def residualT_plot(Results = None, P_bar = None, phases = None, H2O_Liq = None, Fe3Fet_Liq = None, T_cut_C = None, interpolate = True, xlim=None, ylim=None):

    if T_cut_C is None:
        T_cut_C = 12

    if phases is None:
        phases = ['quartz1', 'plagioclase1', 'alkali-feldspar1']

    if H2O_Liq is None and Fe3Fet_Liq is None:
        if len(phases) == 3:
            f, a = plt.subplots(2,2, figsize = (10,8), sharex = True, sharey = True)
            f.tight_layout()
            a[1][0].set_xlabel('P (bars)')
            a[1][1].set_xlabel('P (bars)')
            a[0][0].set_ylabel('Residual T ($\degree$C)')
            a[1][0].set_ylabel('Residual T ($\degree$C)')
            if xlim is not None:
                a[0][0].set_xlim(xlim)
            if ylim is not None:
                a[0][0].set_ylim(ylim)

            m = np.array([['3 Phase Saturation', phases[0] + ' - ' + phases[1]], [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])
            Name = np.array([['Three phase saturation', phases[0] + ' - ' + phases[1]],
                [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])

            for i in range(2):
                for j in range(2):
                    a[i][j].set_title(Name[i,j])
                    a[i][j].plot(P_bar, Results[m[i,j]][0,0,:], 'ok', markerfacecolor="b", label="original", markersize = 8)
                    if interpolate is True:
                        if ~np.isnan(Results['CurveMin'][m[i,j]]['P_min']):
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
            if xlim is not None:
                a.set_xlim(xlim)
            if ylim is not None:
                a.set_ylim(ylim)

            if ~np.isnan(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_min']):
                a.plot(P_bar, Results[phases[0] + ' - ' + phases[1]][0,0,:], 'ok', markerfacecolor="b", label="original", markersize = 8)
                if interpolate is True:
                    a.plot(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_new'], Results['CurveMin'][phases[0] + ' - ' + phases[1]]['y_new'],
                                '-', c="r", label="spline fit")
                    a.plot([np.nanmin(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_new']), np.nanmax(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_new'])],
                                [Results['CurveMin'][phases[0] + ' - ' + phases[1]]['Res_min'], Results['CurveMin'][phases[0] + ' - ' + phases[1]]['Res_min']], ':k')
                    a.plot([Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_min'], Results['CurveMin'][phases[0] + ' - ' + phases[1]]['P_min']],
                                [np.nanmin(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['y_new']) - 5,
                                np.nanmax(Results['CurveMin'][phases[0] + ' - ' + phases[1]]['y_new']) + 5], ':k')

    if H2O_Liq is not None and Fe3Fet_Liq is None:
        if len(Results['CurveMin']) == 4:
            f = plt.figure(figsize = (10,8))
            a1 = f.add_subplot(2,2,1, projection = '3d')
            a2 = f.add_subplot(2,2,2, projection = '3d')
            a3 = f.add_subplot(2,2,3, projection = '3d')
            a4 = f.add_subplot(2,2,4, projection = '3d')

            a = np.array([[a1,a2], [a3, a4]])

            m = np.array([['3 Phase Saturation', phases[0] + ' - ' + phases[1]], [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])
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

                        #A[np.where(A > T_cut_C*2)] = np.nan

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

            m = np.array([['3 Phase Saturation', phases[0] + ' - ' + phases[1]], [phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]])
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

def plot_phaseDiagram(Model = "Holland", Combined = None, P_units = "bar", T_units = "C", 
                      lines = None, T_C = None, P_bar = None, label = True, colormap = None):
    '''
    This function plots the phase diagrams based on the results obtained from thermodynamic models.
    The data should be organized in a pandas dataframe that contains two
    columns correspond to temperature (in Celsius) and pressure (in bars), and the following
    columns correspond to the mass fractions of the phases present.
    The function outputs a phase diagram with different colors for
    different phases, as well as a legend with the names of the phases and a black solid line
    for the phase boundaries (if requested).

    Input Parameters:
        Model: string specifying which thermodynamic model was used to obtain the results. Default
               is "Holland". Other option is "MELTS".
        Combined: pandas dataframe with the results. The default value is None.
        P_units: string specifying the units used for pressure. The default is "bar". Other options
                 are "MPa", "kbar", and "GPa".
        T_units: string specifying the units used for temperature. The default is "C". The other
                 option is "K".
        lines: True/False. Plot phase boundaries.
        T_C: array with the temperature values in Celsius to be plotted. The default is None, and
             the function will use all the temperature values provided in the Combined dataframe.
        P_bar: array with the pressure values in bars to be plotted. The default is None, and the
               function will use all the pressure values provided in the Combined dataframe.

    Output:
        A plot of the phase diagram.
    '''
    Results = Combined.copy()
    Phases = list(Results.columns[Results.columns.str.contains('mass')])

    if T_C is None:
        T_C = np.unique(Results['T_C'])
    if P_bar is None:
        P_bar = np.unique(Results['P_bar'])

    A = [None]*len(Results['T_C'])
    for i in range(len(Results['T_C'])):
        for p in Phases:
            if (~np.isnan(Results[p].loc[i])) & (Results[p].loc[i] > 0.0):
                if A[i] is None:
                    A[i] = p[5:]
                else:
                    A[i] = A[i] + ' ' + p[5:]

        if A[i] is None:
            A[i] = "None"

    Results['Phase'] = A

    Results = Results.sort_values(['T_C', 'P_bar'], ascending=[True, True])

    B = np.zeros(len(Results['Phase']))
    CC = np.sort(np.unique(A))
    C = CC[np.where(CC != 'None')]
    for i in range(len(C)):
        B[np.where(Results['Phase'] == C[i])] = i

    Results['PhaseNo.'] = B
    Results.loc[Results['Phase'] == "None", 'PhaseNo.'] = np.nan

    PT_Results = {}
    PT_Results['T_C'] = Results['T_C'].values.reshape((len(T_C),len(P_bar)))
    PT_Results['P_bar'] = Results['P_bar'].values.reshape((len(T_C),len(P_bar)))
    PT_Results['Phase'] = Results['Phase'].values.reshape((len(T_C),len(P_bar)))
    PT_Results['PhaseNo.'] = Results['PhaseNo.'].values.reshape((len(T_C),len(P_bar)))

    if T_units == "K":
        PT_Results['T_C'] = PT_Results['T_C'] + 273.15

    if P_units == "MPa":
        PT_Results['P_bar'] = PT_Results['P_bar']/10
    elif P_units == "kbar":
        PT_Results['P_bar'] = PT_Results['P_bar']/1000
    elif P_units == "GPa":
        PT_Results['P_bar'] = PT_Results['P_bar']/10000

    f, a = plt.subplots(1,2, figsize = (10,6), gridspec_kw = {'width_ratios': [8,2]})
    a[1].axis("off")
    # im = a[0].pcolormesh(PT_Results['T_C'],
    #                 PT_Results['P_bar'],
    #                 PT_Results['PhaseNo.'], cmap = "Reds", zorder = 2, shading = 'auto')
    im = a[0].imshow(PT_Results['PhaseNo.'].T,
                 extent=[np.nanmin(PT_Results['T_C']), np.nanmax(PT_Results['T_C']),
                         np.nanmin(PT_Results['P_bar']), np.nanmax(PT_Results['P_bar'])],
                 cmap="Reds", zorder=2, aspect='auto', origin='lower')

    #f.colorbar(im, ax = a[1])

    if T_units == "K":
        a[0].set_xlabel('Temperature ($\degree$K)')
    else:
        a[0].set_xlabel('Temperature ($\degree$C)')

    if P_units == "MPa":
        a[0].set_ylabel('Pressure (MPa)')
    elif P_units == "kbar":
        a[0].set_ylabel('Pressure (kbar)')
    elif P_units == "GPa":
        a[0].set_ylabel('Pressure (GPa)')
    else:
        a[0].set_ylabel('Pressure (bar)')

    #j = len(np.unique(np.unique(Results['Phase'])))
    j = len(C)
    if colormap is None:
        cmap = plt.get_cmap('Reds')
    else:
        cmap = plt.get_cmap(colormap)

    #cmap = cmap(np.linspace(1,0,len(np.unique(np.unique(Results['Phase'])))+1))
    cmap = cmap(np.linspace(1,0,len(C)+1))
    for i in range(len(C)): #np.unique(np.unique(Results['Phase'])):
        # if len(Results['PhaseNo.'].values[np.where(Results['Phase'].values == i)]) > 1200:
        T_print = np.nanmedian(Results['T_C'][Results['Phase'] == C[i]])
        P_print = np.nanmedian(Results['P_bar'][Results['Phase'] == C[i]])
        #
        #     p = np.polyfit(Results['T_C'][Results['Phase'] == i], Results['P_bar'][Results['Phase'] == i], 1)
        #
        if label is True:
            a[0].text(T_print, P_print, str(round(i)),
                    horizontalalignment='center',
                    verticalalignment='center',
                fontsize = 8)

        a[1].plot(0,j, 'sk', markerfacecolor = 'w', mec = 'none', ms = 10)
        a[1].text(1, j-0.2, str(round(i)) + '. ' + C[i], fontsize = 10)
        # a[1].set_ylim([1.1,-0.1])
        j = j-1

    a[1].set_xlim([-0.2, 4])

    if lines is not None:
        A = np.zeros((len(T_C)+len(T_C)-1,len(P_bar)+len(P_bar)-1))
        T_mid = np.zeros((len(T_C)+len(T_C)-1,len(P_bar)+len(P_bar)-1))
        P_mid = np.zeros((len(T_C)+len(T_C)-1,len(P_bar)+len(P_bar)-1))
        tol = 1e-9
        for i in range(np.shape(A[0])[0]):
            for j in range(np.shape(A[1])[0]):
                if i % 2 == 0:
                    if j % 2 == 0:
                        T_mid[i,j] = Results['T_C'][int(i/2),int(j/2)]
                        P_mid[i,j] = Results['P_bar'][int(i/2), int(j/2)]
                    elif j % 2 != 0:
                        T_mid[i,j] = np.nanmean(np.array([PT_Results['T_C'][int(i/2),int(j/2-0.5)],
                                                        PT_Results['T_C'][int(i/2),int(j/2+0.5)]]))
                        P_mid[i,j] = np.nanmean(np.array([PT_Results['P_bar'][int(i/2),int(j/2-0.5)],
                                                        PT_Results['P_bar'][int(i/2),int(j/2+0.5)]]))
                        identical = np.allclose(PT_Results['PhaseNo.'][int(i/2),int(j/2-0.5)],
                                                PT_Results['PhaseNo.'][int(i/2),int(j/2+0.5)], rtol=tol, atol=tol)

                        if identical == False:
                            A[i,j] = 1

                elif i % 2 != 0:
                    if j % 2 == 0:
                        T_mid[i,j] = np.nanmean(np.array([PT_Results['T_C'][int(i/2-0.5),int(j/2)],
                                                        PT_Results['T_C'][int(i/2+0.5),int(j/2)]]))
                        P_mid[i,j] = np.nanmean(np.array([PT_Results['P_bar'][int(i/2-0.5),int(j/2)],
                                                        PT_Results['P_bar'][int(i/2+0.5),int(j/2)]]))
                        identical = np.allclose(PT_Results['PhaseNo.'][int(i/2-0.5),int(j/2)],
                                                PT_Results['PhaseNo.'][int(i/2+0.5),int(j/2)], rtol=tol, atol=tol)

                        if identical == False:
                            A[i,j] = 1

                    elif j % 2 != 0:
                        T_mid[i,j] = np.nanmean(np.array([PT_Results['T_C'][int(i/2-0.5),int(j/2-0.5)],
                                                        PT_Results['T_C'][int(i/2-0.5),int(j/2+0.5)],
                                                        PT_Results['T_C'][int(i/2+0.5),int(j/2-0.5)],
                                                        PT_Results['T_C'][int(i/2+0.5),int(j/2+0.5)]]))
                        P_mid[i,j] = np.nanmean(np.array([PT_Results['P_bar'][int(i/2-0.5),int(j/2-0.5)],
                                                        PT_Results['P_bar'][int(i/2-0.5),int(j/2+0.5)],
                                                        PT_Results['P_bar'][int(i/2+0.5),int(j/2-0.5)],
                                                        PT_Results['P_bar'][int(i/2+0.5),int(j/2+0.5)]]))
                        identical = np.allclose(PT_Results['PhaseNo.'][int(i/2-0.5),int(j/2-0.5)], PT_Results['PhaseNo.'][int(i/2+0.5),int(j/2-0.5)], rtol=tol, atol=tol) \
                                    and np.allclose(PT_Results['PhaseNo.'][int(i/2-0.5),int(j/2-0.5)], PT_Results['PhaseNo.'][int(i/2-0.5),int(j/2+0.5)], rtol=tol, atol=tol) \
                                    and np.allclose(PT_Results['PhaseNo.'][int(i/2-0.5),int(j/2-0.5)], PT_Results['PhaseNo.'][int(i/2+0.5),int(j/2+0.5)], rtol=tol, atol=tol)

                        if identical == False:
                            A[i,j] = 1

        a[0].contour(T_mid, P_mid/10, A, colors = 'k', linewidths = 0.5)

    plt.show()

    return f, a

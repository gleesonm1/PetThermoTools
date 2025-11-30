import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint, Point, Polygon
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mc
from petthermotools.GenFuncs import *
from itertools import zip_longest

labelling = {'SiO2': 'SiO$_2$',
             'TiO2': 'TiO$_2$',
             'Al2O3': 'Al$_2$O$_3$',
             'FeOt': 'FeO$_t$',
             'Cr2O3': 'Cr$_2$O$_3$',
             'MnO': 'MnO',
             'MgO': 'MgO',
             'CaO': 'CaO',
             'Na2O': 'Na$_2$O',
             'K2O': 'K$_2$O',
             'H2O': 'H$_2$O',
             'CO2': 'CO$_2$'}

def harker(Results=None, x_axis="MgO", y_axis=("SiO2", "TiO2", "Al2O3", "Cr2O3", "FeOt", "CaO", "Na2O", "K2O"),
    phase="liquid1", line_color=None, data=None, d_color=None, d_marker=None,
    legend=True, legend_loc=None,
    xlim=None, ylim=None):
    """
    Generates a Harker diagram (oxide vs. MgO or another oxide) from model results.

    Parameters
    ----------
    Results : dict or None, optional
        MELTS or MAGEMin results dictionary. If multiple results, they are overlaid.
    x_axis : str, default="MgO"
        Oxide to use for the x-axis.
    y_axis : tuple of str, default=("SiO2","TiO2","Al2O3","Cr2O3","FeOt","CaO","Na2O","K2O")
        Oxides to plot on y-axes. Plotted in groups of three per row.
    phase : str, default="liquid1"
        Phase to plot (must be in petthermotools.GenFuncs.Names).
    line_color : str or None, optional
        Line color. If None, uses matplotlib's default cycle.
    data : pd.DataFrame, dict of DataFrames, or str, optional
        External data to plot for comparison. If str, treated as CSV path.
        If dict, plots each DataFrame with different marker/color.
    d_colors : list of str or None, optional
        Colors for external datasets. Cycles if fewer than number of datasets.
    d_markers : list of str or None, optional
        Markers for external datasets. Cycles if fewer than number of datasets.
    legend : bool, default=True
        Whether to display a legend.
    legend_loc : tuple(int, int) or None, optional
        Location of legend in subplot grid, as (row, col). If None, placed on last axis.
    xlim, ylim : tuple or None, optional
        Limits for x and y axes.

    Returns
    -------
    f : matplotlib.figure.Figure
        Figure object.
    a : np.ndarray of Axes
        Array of Axes objects.
    """
    if Results is None:
        raise ValueError("Results cannot be None. Provide MELTS or MAGEMin results dictionary.")

    if isinstance(y_axis, str):
        y_axis = (y_axis,)

    y_axis = list(y_axis)

    # Data loading
    if data is not None:
        if isinstance(data, str):
            if data.contains('.csv'):
                data = pd.read_csv(data)
            else:
                data = pd.read_excel(data)
        elif isinstance(data, dict):
            # validate
            for k, v in data.items():
                if not isinstance(v, pd.DataFrame):
                    raise ValueError(f"data['{k}'] must be a DataFrame, got {type(v)}")
                
    # -------------------- filter y_axis to only include variables present --------------------
    valid_y = []
    for y in y_axis:
        found = False

        # Check Results
        for r in (Results.keys() if isinstance(Results, dict) else [Results]):
            Res = Results[r] if isinstance(Results, dict) else Results
            if "All" in Res:
                cols = Res["All"].columns
                if (y + Names[phase]) in cols or y in cols:
                    if np.nanmax(Res["All"][y + Names[phase]]) > 0.0:
                        found = True
                    break

        # Check data (dict or DataFrame)
        if not found and data is not None:
            if isinstance(data, dict):
                for df in data.values():
                    cols = df.columns
                    if (y + Names[phase]) in cols or y in cols:
                        found = True
                        break
            else:
                cols = data.columns
                if (y + Names[phase]) in cols or y in cols:
                    found = True

        if found:
            valid_y.append(y)

    # Replace y_axis with only valid variables
    y_axis = valid_y

    # -------------------- subplot grid handling --------------------
    n_vars = len(y_axis)

    if n_vars == 4:  # special case: 2x2 grid
        ncols, nrows = 2, 2
        rows = list(zip_longest(*[iter(y_axis)] * ncols, fillvalue=None))
    else:
        ncols = 3
        nrows = int(np.ceil(n_vars / ncols))
        rows = list(zip_longest(*[iter(y_axis)] * ncols, fillvalue=None))

    # scale figure size dynamically
    fig_width = 3.2 * ncols
    fig_height = 3.0 * nrows
    f, a = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
    a = np.atleast_2d(a)

    # Color handling for models
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    if line_color is None:
        def pick_color(i): return color_cycle[i % len(color_cycle)]
    else:
        def pick_color(i): return line_color

    # Color/marker cycles for user data
    d_color_cycle = d_color if d_color else color_cycle
    d_marker_cycle = d_marker if d_marker else ["o", "s", "^", "D", "v", "P", "*"]

    # -------------------- helper functions --------------------
    def plot_panel(ax, x_var, y_var, idx=None, Res = None):
        """Plot one panel (model) on given Axes."""
        suffix = Names[phase] if phase in Names else ""
        xcol, ycol = x_var + suffix, y_var + suffix

        # plot model data
        if Res is not None:
            if xcol in Results[Res]["All"] and ycol in Results[Res]["All"]:
                x = Results[Res]["All"][xcol].values
                y = Results[Res]["All"][ycol].values
                ax.plot(x, y, color=pick_color(idx), label=Res)
        else:
            if xcol in Results["All"] and ycol in Results["All"]:
                x = Results["All"][xcol].values
                y = Results["All"][ycol].values
                ax.plot(x, y, color='k')

    def plot_data(ax, x_var, y_var):
        # plot external data
        suffix = Names[phase] if phase in Names else ""
        xcol, ycol = x_var + suffix, y_var + suffix

        if data is not None:
            if isinstance(data, pd.DataFrame):
                datasets = {"data": data}
            else:
                datasets = data

            for k, df in datasets.items():
                c = d_color_cycle[list(datasets.keys()).index(k) % len(d_color_cycle)]
                m = d_marker_cycle[list(datasets.keys()).index(k) % len(d_marker_cycle)]

                # try with suffix first
                if xcol in df.columns:
                    dx = df[xcol].values
                elif x_var in df.columns:
                    dx = df[x_var].values
                else:
                    dx = None
                    print(f"x axis variable  {x_var} not found in data")

                if ycol in df.columns:
                    dy = df[ycol].values
                elif y_var in df.columns:
                    dy = df[y_var].values
                else:
                    dy = None
                    print(f"y axis variable  {y_var} not found in data")

                if dx is not None and dy is not None:
                    ax.plot(dx, dy, m,
                            markerfacecolor=c,
                            markeredgecolor="k",
                            markersize=4,
                            linestyle="None",
                            label=k)

        ax.set(xlabel=labelling[x_var], ylabel=labelling[y_var])

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # -------------------- main plotting --------------------
    for i, row in enumerate(rows):
        for j, y in enumerate(row):
            if y is None:
                a[i, j].axis("off")
                continue
            if 'All' in Results.keys():
                if data is not None:
                    plot_data(a[i,j], x_axis, y)
                plot_panel(a[i,j], x_axis, y)
            else:
                if data is not None:
                    plot_data(a[i,j], x_axis, y)
                for idx, Res in enumerate(Results):
                    plot_panel(a[i, j], x_axis, y, idx=idx, Res = Res)

    # -------------------- legend --------------------
    empty_axes = []
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if i * a.shape[1] + j >= len(y_axis):  # beyond valid y variables
                empty_axes.append(a[i][j])

    handles, labels = a[0][0].get_legend_handles_labels()
    if empty_axes and handles:
        empty_axes[0].legend(handles, labels, loc = "center")
        empty_axes[0].axis("off")
    else:
        if legend:
            if legend_loc is None:
                loc_i, loc_j = len(rows) - 1, 0
            else:
                loc_i, loc_j = legend_loc
            a[loc_i, loc_j].legend()

    f.tight_layout()

    return f, a

def plot_surfaces(Results = None, P_bar = None, phases = None, H2O_Liq = None):
    if H2O_Liq is None:
        f, a = plt.subplots(1,1, figsize = (5,4))
        a.set_xlabel('P (bars)')
        a.set_ylabel('T ($\degree$C)')
        for i in range(len(phases)):
            try:
                if i == 0:
                    a.plot(P_bar, Results['Output'][phases[0]], '-r', linewidth = 2, label = phases[i])
                if i == 1:
                    a.plot(P_bar, Results['Output'][phases[1]], '-b', linewidth = 2, label = phases[i])
                if i == 2:
                    a.plot(P_bar, Results['Output'][phases[2]], '-k', linewidth = 2, label = phases[i])
            except:
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
                    try:
                        a[i][j].plot(P_bar, Results['Output'][m[i,j]], 'ok', markerfacecolor="b", label="original", markersize = 8)
                    except:
                        print("You are using the old find_mineral_saturation function.\n This will soon be removed, please transition to the mineral_cosaturation function.")
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


def phase_plot(Results, x_axis = None, y_axis = None, cmap = "Reds"):
    """
    Create stacked phase mass-fraction plots from thermodynamic model results.

    This function generates diagrams of phase proportions from the results of 
    crystallization or melting simulations. Mass fractions of crystalline and 
    liquid phases are stacked either along the x-axis or y-axis, depending on 
    user specification, to visualize how phase proportions evolve with pressure, 
    temperature, or another variable.

    Parameters
    ----------
    Results : dict
        Dictionary containing model outputs. It should have one of the following structures:
        - **Single-run results**:
          ```
          Results = {
              "Mass": pandas.DataFrame,   # columns = phase names or phase_cumsum
              "All":  pandas.DataFrame    # contains the axis variable (e.g., T_C, P_bar)
          }
          ```
        - **Multi-run results**:
          ```
          Results = {
              run_label: {
                  "Mass": pandas.DataFrame,
                  "All":  pandas.DataFrame
              },
              ...
          }
          ```

        The `"Mass"` DataFrame must include mass fractions for each phase
        (either raw values or cumulative values with `"_cumsum"` suffix).
        The `"All"` DataFrame must contain the column specified by `x_axis` or `y_axis`.

    x_axis : str, optional
        Column name in `Results["All"]` (or equivalent) to plot on the x-axis. 
        If provided, `y_axis` must be `None`. Typically something like `"T_C"`.
    
    y_axis : str, optional
        Column name in `Results["All"]` (or equivalent) to plot on the y-axis. 
        If provided, `x_axis` must be `None`. Typically `"P_bar"`.
    
    cmap : str, default = "Reds"
        Matplotlib colormap used to assign colors to phases.

    Returns
    -------
    fig : matplotlib.figure.Figure or list of Figures
        - If `Results` contains a single run: a single `Figure`.
        - If `Results` contains multiple runs: a list of `Figure` objects, one per run.

    axes : matplotlib.axes.Axes or list of Axes
        - If `Results` contains a single run: a single `Axes` object.
        - If `Results` contains multiple runs: a list of `Axes` objects, one per run.

    Notes
    -----
    - If both `x_axis` and `y_axis` are provided, the function raises a `ValueError`.
    - Phases are automatically ordered:
        1. By the index where they first appear (highest pressure or temperature).
        2. By the order specified in the petthermotools dictionaries `Names` and `Names_MM`, 
           if available.
    - The liquid phase (`"liq1"` or `"liquid1"`) is always plotted last.
    - A legend with readable phase labels is added outside 
      the plot area.

    Examples
    --------
    >>> fig, ax = phase_plot(Results, y_axis="P_bar")
    # Plots stacked phase proportions vs. pressure

    >>> fig, axes = phase_plot(MultiRunResults, x_axis="T_C")
    # Creates one stacked phase plot per run, vs. temperature
    """
    if x_axis is not None and y_axis is not None:
        raise ValueError("Please provide either a x-axis or y-axis parameter to plot the mass fractions against")
    
    def makeplot(Mass, Res, title = None):
        # --- Identify whether _cumsum columns exist ---
        use_cumsum = any(col.endswith("_cumsum") for col in Mass.columns)
        suffix = "_cumsum" if use_cumsum else ""

        # --- Identify phases ---
        exclude_cols = {'T_C', None}
        phases = [
            col.replace(suffix, "")
            for col in Mass.columns
            if col.endswith(suffix) or (not use_cumsum and not col.endswith("_cumsum") and col not in exclude_cols)
        ]

        # --- Handle liquid phase ---
        liquid_name = None
        for liq_name in ["liq1", "liquid1"]:
            if liq_name in Mass.columns:
                liquid_name = liq_name
                if liq_name in phases:
                    phases.remove(liq_name)
                break  # use the first match

        # --- Create dictionary priority (order in Names/Names_MM) ---
        phase_priority = {}
        for order_dict in [Names, Names_MM]:
            for i, key in enumerate(order_dict.keys()):
                phase_priority[key] = i  # smaller i = higher priority in tie

        # --- Sort crystalline phases ---
        phase_first_index = {}
        for p in phases:
            col = p + suffix if (use_cumsum and p + suffix in Mass.columns) else p
            vals = Mass[col].values
            nz = np.flatnonzero(vals > 0)
            phase_first_index[p] = nz[0] if len(nz) > 0 else np.inf

        def sort_key(p):
            return (phase_first_index[p], phase_priority.get(p, 1e6))

        phases = sorted(phases, key=sort_key)

        # --- Append liquid last ---
        if liquid_name is not None:
            phases.append(liquid_name)

        # --- Assign colors ---
        c = cm.get_cmap(cmap, len(phases))
        PhaseColors = {p: c(i) for i, p in enumerate(phases)}

        # --- Setup figure ---
        if y_axis == "P_bar":
            f, a = plt.subplots(figsize=(3.5, 5))
        else:
            f, a = plt.subplots(figsize=(4, 3.5))

        # --- Determine orientation ---
        coord = Res[y_axis] if y_axis else Res[x_axis]
        horizontal = y_axis is not None

        Stop = np.zeros(len(Mass))
        for p in phases:
            col = p + suffix if (use_cumsum and p + suffix in Mass.columns) else p
            vals = Mass[col].values

            # --- Legend label mapping ---
            label = Names.get(p, Names_MM.get(p, p))
            if '_' in label:
                label = label[1:]

            if horizontal:
                a.fill_betweenx(coord, Stop, Stop + vals, color=PhaseColors[p], alpha=0.75, lw=0, label=label)
            else:
                a.fill_between(coord, Stop, Stop + vals, color=PhaseColors[p], alpha=0.75, lw=0, label=label)
            Stop += vals

        # --- Labels ---
        if horizontal:
            if y_axis == "P_bar":
                a.set_ylabel("Pressure (bar)")
                a.invert_yaxis()
            else:
                a.set_ylabel(y_axis)
            a.set_xlabel("Mass fraction")
        else:
            a.set_xlabel(x_axis)
            a.set_ylabel("Mass fraction")

        if title is not None:
            a.set_title(title)

        # --- Remove whitespace ---
        a.margins(0)
        f.tight_layout()

        a.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))

        return f, a
    fig = []
    axes = []

    if 'All' in Results.keys():
        fig, axes = makeplot(Results['mass_g'], Results['All'])
    else:
        for val in Results.keys():
            f, a = makeplot(Results[val]['mass_g'],
                            Results[val]['All'],
                            title = val)
            fig.append(f)
            axes.append(a)
    return fig, axes


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
    Phases = list(Results.columns[(Results.columns.str.contains('mass_g') & (~Results.columns.str.contains('%')) & (~Results.columns.str.contains('mole')))])

    if T_C is None:
        T_C = np.unique(Results['T_C'])
    if P_bar is None:
        P_bar = np.unique(Results['P_bar'])

    A = [None]*len(Results['T_C'])
    for i in range(len(Results['T_C'])):
        for p in Phases:
            if (~np.isnan(Results[p].loc[i])) & (Results[p].loc[i] > 0.0):
                if A[i] is None:
                    A[i] = p[7:]
                else:
                    A[i] = A[i] + ' ' + p[7:]

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

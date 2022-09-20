import numpy as np
import pandas as pd
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




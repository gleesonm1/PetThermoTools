import multiprocessing

# ---- Guard required for ProcessPoolExecutor on Windows ------------------- #
if __name__ == "__main__":
    multiprocessing.freeze_support()

    import numpy as np
    import pandas as pd
    import petthermotools as ptt
    from time import time

    Kil = pd.read_excel('Kilauea_olivines.xlsx')
    LL8 = Kil.loc[Kil['Label'].str.contains("LL8")].reset_index(drop = True)

    times = [time()]
    nGibbs_equil = ptt.equilibrate_multi(Model = "nMELTSv1.0.2NoProp",
                                bulk = LL8,
                                P_bar = np.linspace(50,900,LL8.shape[0]),
                                T_C = LL8['Temp'])

    print(nGibbs_equil['mass_g'])
    times.append(time())
    print(f"nGibbs equilibration time for {LL8.shape[0]} samples: {times[-1] - times[0]:.2f} seconds")

    Kil_starting = Kil.loc[Kil['MgO'] == Kil['MgO'].max()].squeeze()

    isobaric = ptt.isobaric_crystallisation(Model = "nMELTSv1.0.2NoProp",
                                            bulk = Kil_starting,
                                            P_bar = [250,500,1000,2000],
                                            find_liquidus = True,
                                            T_end_C = 1050,
                                            dt_C = 2,
                                            label = "P_bar")

    print(f"Keys in isobaric: {list(isobaric['Pressure=2000.00'].keys())}")
    print(f"Liquid Table of isobaric: {isobaric['Pressure=250.00']['mass_g']}")
    times.append(time())
    print(f"nGibbs isobaric crystallisation time for 4 paths: {times[-1] - times[-2]:.2f} seconds")

    oxSearch = ptt.isobaric_crystallisation(Model = "nMELTSv1.0.2",
                                            bulk = Kil_starting,
                                            P_bar = np.linspace(200,2000,10),
                                            T_start_C = 1500,
                                            T_end_C = 1050,
                                            dt_C = 2,
                                            fO2_offset = [-2,-1,0,1,2])

    print(f"Keys in oxSearch: {list(oxSearch['logfO2-QFM=0.00'].keys())}")
    print(f"Liquid Table of oxSearch: {oxSearch['logfO2-QFM=0.00']['liquid1_prop']}")
    times.append(time())
    print(f"nGibbs properties for {5*10*((1500-1050)/2)} samples: {times[-1] - times[-2]:.2f} seconds")

    watergrid = ptt.isobaric_crystallisation(Model = "nMELTSv1.0.2NoProp",
                                            bulk = Kil_starting,
                                            P_bar = np.linspace(200,2000,10),
                                            T_start_C = 1500,
                                            T_end_C = 1050,
                                            dt_C = 2,
                                            H2O_init = [0,0.1,0.5,1,2])

    print(f"Keys in watergrid: {list(watergrid['H2O=0.00'].keys())}")
    print(f"Liquid Table of watergrid: {watergrid['H2O=0.00']['liquid1']}")
    times.append(time())
    print(f"nGibbs water grid time for 10 paths: {times[-1] - times[-2]:.2f} seconds")

    waterpath = ptt.multi_path(Model = "nMELTSv1.0.2",
                                            bulk = Kil_starting,
                                            P_path_bar = np.linspace(2000,200,500),
                                            T_path_C = np.linspace(1500,950,500),
                                            H2O_init = 0.1)

    print(f"Keys in waterpath: {list(waterpath.keys())}")
    print(f"waterpath liquid Fe3FeT ratio: {waterpath['liquid1']['Fe3Fet_Liq']}")
    times.append(time())
    print(f"nGibbs water path time for 500 steps: {times[-1] - times[-2]:.2f} seconds")
    print(Kil_starting)
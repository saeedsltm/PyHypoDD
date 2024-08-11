import os

import proplot as plt
import seaborn as sns
from numpy import linspace
from pandas import Series, read_csv
from pyproj import Proj

from core.Extra import loadxyzm


def plotSeismicityMap(config):
    print("+++ Plotting seismicity map ...")
    clat = config["Region"]["CentralLat"]
    clon = config["Region"]["CentralLon"]
    EventsMaxDepth = config["Figures"]["EventsMaxDepth"]
    outName = f"{config['Region']['RegionName']}"
    proj = Proj(f"+proj=sterea\
            +lon_0={clon}\
            +lat_0={clat}\
            +units=km")
    catalog_ini_path = os.path.join(
        "results", f"xyzm_{outName}_initial.dat")
    catalog_hdd_path = os.path.join(
        "results", f"xyzm_{outName}_hypodd.dat")
    report_ini, report_hdd = loadxyzm(catalog_ini_path,
                                      catalog_hdd_path)
    conds = (report_hdd.ORT.notna()) & (
        report_hdd.Lon.notna()) & (report_hdd.Lat.notna())
    report_ini = report_ini[conds]
    for db in [report_ini, report_hdd]:
        db[["x", "y"]] = db.apply(
            lambda x: Series(
                proj(longitude=x.Lon, latitude=x.Lat)), axis=1)
        db["z"] = db["Dep"]
    stationPath = os.path.join("stations", "stations.csv")
    stations_df = read_csv(stationPath)

    axShape = [
        [1, 2],
        [3, 4]
    ]
    plt.rc.update(
        {'fontsize': 7, 'legend.fontsize': 6, 'label.weight': 'bold'})
    fig, axs = plt.subplots(axShape, share=True)
    axs.format(
        xlabel="Easting (km)",
        ylabel="Northing (km)",
        xlocator=("maxn", 5),
        ylocator=("maxn", 5))
    [ax.grid(True, ls=":") for ax in axs]

    sc = axs[0].scatter(report_ini.x, report_ini.y, s=10*10**report_ini.Mag,
                        cmap="plasma_r", c=report_ini.Dep,
                        mec="k", mew=0.2, vmax=EventsMaxDepth)
    sc = axs[1].scatter(report_hdd.x, report_hdd.y, s=10*10**report_hdd.Mag,
                        cmap="plasma_r", c=report_hdd.Dep,
                        mec="k", mew=0.2, vmax=EventsMaxDepth)
    axs[1].colorbar(
        sc, loc="r", label="Depth (km)", ticks=linspace(0, EventsMaxDepth, 10),
        format="%.0f")

    sns.kdeplot(x=report_ini.x, y=report_ini.y, shade=True,
                bw_adjust=0.5, ax=axs[2], cmap="plasma_r")
    sns.kdeplot(x=report_hdd.x, y=report_hdd.y, shade=True, bw_adjust=0.5,
                ax=axs[3], cmap="plasma_r", cbar=True,
                cbar_kws={"title": "Density", "format": "%.1e"})

    for ax in axs:
        ax.plot(stations_df.x, stations_df.y, marker="^", ms=3,
                ls="", c="white", mec="k", mew=1.0, autoreverse=False)

    fig.save(os.path.join("results", "seismicity.png"))  # type: ignore

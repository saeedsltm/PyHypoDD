import os
from pathlib import Path

from obspy import read_events
from obspy.geodetics.base import gps2dist_azimuth as gps
from pandas import DataFrame, Series
from pyproj import Proj
from yaml import dump, safe_load

from core.GetStationInfo import download_IRSSI


def GetStationListFromCatalog(config):
    print("+++ Generating list of used stations from input catalog ...")
    Path("stations").mkdir(parents=True, exist_ok=True)
    catalogPath = config["Files"]["InputCatalogFileName"]
    catalog = read_events(catalogPath)
    stationsList = []
    for event in catalog:
        picks = event.picks
        codes = [pick.waveform_id.station_code for pick in picks]
        for code in codes:
            if code not in stationsList:
                stationsList.append(code)
    stationsList = sorted(stationsList, key=lambda x: (len(x), x))
    with open(os.path.join("stations", "usedStations.yml"), "w") as outfile:
        dump({"usedStations": stationsList},
             outfile,
             default_flow_style=False,
             sort_keys=False)


def CreatInputStationFile(config):
    print("+++ Creating HypoDD station file ...")
    clat = config["Region"]["CentralLat"]
    clon = config["Region"]["CentralLon"]
    radius = config["Region"]["Radius"]
    proj = Proj(f"+proj=sterea\
            +lon_0={clon}\
            +lat_0={clat}\
            +units=km")
    data = []
    missedStations = []
    with open(os.path.join("stations", "usedStations.yml")) as infile:
        usedStations = safe_load(infile)
    statioFileNames = download_IRSSI("stations")
    stationsInfo = {}
    for name in statioFileNames:
        with open(os.path.join("stations", name)) as infile:
            info = safe_load(infile)
            stationsInfo.update(info)
    for station in usedStations["usedStations"]:
        if station in stationsInfo:
            info = {"code": station,
                    "lat": stationsInfo[station][-1]["latitude"],
                    "lon": stationsInfo[station][-1]["longitude"],
                    "elv": stationsInfo[station][-1]["elevation"]}
            data.append(info)
        else:
            missedStations.append(station)
    stations_df = DataFrame(data)
    stations_df[["x", "y"]] = stations_df.apply(
        lambda x: Series(
            proj(longitude=x.lon, latitude=x.lat)), axis=1)
    stations_df[["r"]] = stations_df.apply(
        lambda x: Series(
            gps(clat, clon, x.lat, x.lon)[0]*1e-3), axis=1)
    stations_df["z"] = stations_df["elv"]
    stations_df.sort_values(by=["r"], inplace=True)
    unusedStations_df = stations_df[stations_df.r > radius]
    stations_df = stations_df[stations_df.r <= radius]
    stations_df.to_csv(os.path.join("stations", "stations.csv"),
                       index=False, float_format="%8.3f")
    unusedStations_df.to_csv(os.path.join("stations", "unusedStations.csv"),
                             index=False, float_format="%8.3f")
    with open(os.path.join("stations", "missedStations.yml"), "w") as outfile:
        dump({"missedStations": missedStations},
             outfile,
             default_flow_style=False,
             sort_keys=False)

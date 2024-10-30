import os
from pathlib import Path
import requests
from bs4 import BeautifulSoup

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
    print("+++ Reading catalog using Obspy ...")
    catalog = read_events(catalogPath)
    stationsList = []
    for event in catalog:
        picks = event.picks
        codes = [pick.waveform_id.station_code for pick in picks]
        for code in codes:
            if code not in stationsList:
                stationsList.append(code)
    stationsList = sorted(stationsList, key=lambda x: (len(x), x))
    with open(os.path.join("stations", "stationsInCatlog.yml"), "w") as outfile:
        dump({"catalogStations": stationsList},
             outfile,
             default_flow_style=False,
             sort_keys=False)


def downloadMissedStationFromISC(missedStations):
    sta_list = "%2C".join([f"{code}" for code in missedStations])
    data = []
    foundedStations = []
    # URL for the station search
    url = f"https://www.isc.ac.uk/cgi-bin/stations?stnsearch=STN&sta_list={sta_list}&stn_ctr_lat=&stn_ctr_lon=&stn_radius=&max_stn_dist_units=deg&stn_bot_lat=&stn_top_lat=&stn_left_lon=&stn_right_lon=&stn_srn=&stn_grn="
    # Send a GET request to the URL
    response = requests.get(url)
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the HTML content
        soup = BeautifulSoup(response.content, 'html.parser')
        text = soup.text.splitlines()
        for line in text:
            code = line[:5].strip()
            if code in missedStations:
                lat = float(line[59:67])
                lon = float(line[69:77])
                elv = float(line[79:88])
                info = {
                    "code": code,
                    "lat": lat,
                    "lon": lon,
                    "elv": elv
                }
                foundedStations.append(code)
                data.append(info)
    missedStations = list(set(missedStations) - set(foundedStations))
    return data, missedStations


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
    with open(os.path.join("stations", "stationsInCatlog.yml")) as infile:
        usedStations = safe_load(infile)
    statioFileNames = download_IRSSI("stations")
    stationsInfo = {}
    for name in statioFileNames:
        with open(os.path.join("stations", name)) as infile:
            info = safe_load(infile)
            stationsInfo.update(info)
    for station in usedStations["catalogStations"]:
        if station in stationsInfo:
            info = {"code": station,
                    "lat": stationsInfo[station][-1]["latitude"],
                    "lon": stationsInfo[station][-1]["longitude"],
                    "elv": stationsInfo[station][-1]["elevation"]}
            data.append(info)
        else:
            missedStations.append(station)
    newData, missedStations = downloadMissedStationFromISC(missedStations)
    data.extend(newData)
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
    stations_df.to_csv(os.path.join("stations", "usedStations.csv"),
                       index=False, float_format="%8.3f")
    unusedStations_df.to_csv(os.path.join("stations", "unusedStations.csv"),
                             index=False, float_format="%8.3f")
    with open(os.path.join("stations", "missedStations.yml"), "w") as outfile:
        dump({"missedStations": missedStations},
             outfile,
             default_flow_style=False,
             sort_keys=False)

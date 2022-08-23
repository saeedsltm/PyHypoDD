import os
import warnings
from math import sqrt

from LatLon import lat_lon as ll
from numpy import loadtxt, mean, nan, round_
from obspy import read_events
from obspy.geodetics.base import degrees2kilometers as d2k
from pandas import DataFrame

from utils.extra import handleNone

warnings.filterwarnings("ignore")


def getHer(event):
    """Get horizontal error of event
    Args:
        event (obspy.event): an obspy event
    Returns:
        float: event horizontal error
    """
    if event.origins[0].latitude_errors.uncertainty:
        x = event.origins[0].latitude_errors.uncertainty
        y = event.origins[0].longitude_errors.uncertainty
        return round(d2k(sqrt(x**2 + y**2)), 1)
    else:
        return None


def getZer(event):
    """Get depth error of event
    Args:
        event (obspy.event): an obspy event
    Returns:
        float: event depth error
    """
    if event.origins[0].depth_errors.uncertainty:
        return event.origins[0].depth_errors.uncertainty*0.001
    else:
        return None


def readVelocityFile(stationFile):
    """reading velocity model from "STATION0.HYP" file 
    Args:
        stationFile (str): station file in "NORDIC" format
    Returns:
        (dict): a dictionary contains velocity model
    """
    emptyLines = 0
    velocityModelDict = {"Vp": [], "Z": [], "VpVs": 1.73, "Interfaces": []}
    with open(stationFile) as f:
        for l in f:
            if not l.strip():
                emptyLines += 1
            if emptyLines == 2 and l.strip():
                Vp, Z = [float(x) for x in l.split()[:2]]
                velocityModelDict["Vp"].append(Vp)
                velocityModelDict["Z"].append(Z)
            if emptyLines == 2 and len(l) > 20:
                # _, Z = [float(x) for x in l.split()[:2]]
                velocityModelDict["Interfaces"].append(l[21])
            if emptyLines == 3 and l.strip():
                VpVs = float(l[16:20])
                velocityModelDict["VpVs"] = VpVs
                break
    return velocityModelDict


def readStationFile(stationFile):
    """read station information from "STATION0.HYP" file
    Args:
        stationFile (str): station file in "NORDIC" format
    Returns:
        dict: a dictionary contains stations information
    """
    emptyLines = 0
    stationsDict = {}
    with open(stationFile) as f:
        for l in f:
            if not l.strip():
                emptyLines += 1
            if emptyLines == 1 and l.strip():
                code = l[:6].strip()
                lat = ll.Latitude(degree=int(
                    l[6:8]), minute=float(l[8:13])).decimal_degree  # type: ignore
                lon = ll.Longitude(degree=int(
                    l[15:17]), minute=float(l[17:22])).decimal_degree  # type: ignore
                elv = float(l[23:27])
                stationsDict[code] = {"Lat": lat, "Lon": lon, "Elv": elv}
    return stationsDict


def writeStationFile(stationFile, stationCodes, stationLats, stationLons, stationElvs, Vp, Z, Interfaces, trialDepth, xNear, xFar, VpVs):
    """Write STATION0.HYP file
    Args:
        stationFile (str): output file name,
        stationCodes (list): a list contains stations code,
        stationLats (list): a list contains stations latitude,
        stationLons (list): a list contains stations longitude,
        stationElvs (list): a list contains stations elevation,
        Vp (list): a list contains layers velocity,
        Z (list): a list contains layers depth,
        Interfaces (list): a list contains layers interface,
        trialDepth (float): starting trial depth,
        xNear (float): distance where the weight is maximum within,
        xFar (float): distance where the weight is 0 beyond,
        VpVs (float): Vp to Vs ratio,
    Returns:
        str: saved station file name.
    """
    with open(os.path.join("files", "resets.dat")) as f, open(stationFile, "w") as g:
        for l in f:
            g.write(l)
        g.write("\n\n")
        for code, lat, lon, elv in zip(stationCodes, stationLats, stationLons, stationElvs):
            lat = ll.Latitude(lat)
            lon = ll.Longitude(lon)
            g.write("  {code:4s}{latDeg:2.0f}{latMin:05.2f}N {lonDeg:2.0f}{lonMin:05.2f}E{elv:4.0f}\n".format(
                code=code,
                latDeg=lat.degree, latMin=lat.decimal_minute,
                lonDeg=lon.degree, lonMin=lon.decimal_minute,
                elv=elv
            ))
        g.write("\n")
        for v, z, i in zip(Vp, Z, Interfaces):
            g.write(" {v:5.2f}  {z:6.3f}       {i:1s}     \n".format(
                v=v, z=z, i=i
            ))
        g.write("\n")
        g.write("{trialDepth:4.0f}.{xNear:4.0f}.{xFar:4.0f}. {VpVs:4.2f}".format(
            trialDepth=trialDepth, xNear=xNear, xFar=xFar, VpVs=VpVs
        ))
        g.write("\nNew")
    return stationFile


def catalog2xyzm(hypInp, catalogFileName):
    """Convert catalog to xyzm file format
    Args:
        hypInp (str): file name of NORDIC file
        catalogFileName (str): file name of xyzm.dat file
    """
    cat = read_events(hypInp)
    magnitudes = loadtxt("magnitudes.dat")
    outputFile = "xyzm_{catalogFileName:s}.dat".format(
        catalogFileName=catalogFileName.split(".")[0])
    catDict = {}
    for i, event in enumerate(cat):
        arrivals = event.origins[0].arrivals
        ort = event.origins[0].time
        lat = event.origins[0].latitude
        lon = event.origins[0].longitude
        try:
            dep = event.origins[0].depth*0.001
        except TypeError:
            dep = nan
        try:
            nus = event.origins[0].quality.used_station_count
        except AttributeError:
            nus = nan
        nuP = len(
            [arrival.phase for arrival in arrivals if "P" in arrival.phase.upper()])
        nuS = len(
            [arrival.phase for arrival in arrivals if "S" in arrival.phase.upper()])
        mds = handleNone(
            min([handleNone(arrival.distance) for arrival in event.origins[0].arrivals]), degree=True)
        ads = round_(handleNone(
            mean([handleNone(arrival.distance) for arrival in event.origins[0].arrivals]), degree=True), 2)
        try:
            gap = handleNone(event.origins[0].quality.azimuthal_gap)
        except AttributeError:
            gap = nan
        try:
            rms = handleNone(event.origins[0].quality.standard_error)
        except AttributeError:
            rms = nan
        erh = getHer(event)
        erz = getZer(event)
        catDict[i] = {
            "ORT": ort,
            "Lon": lon,
            "Lat": lat,
            "Dep": dep,
            # "mag":mag,
            "Nus": nus,
            "NuP": nuP,
            "NuS": nuS,
            "ADS": ads,
            "MDS": mds,
            "GAP": gap,
            "RMS": rms,
            "ERH": erh,
            "ERZ": erz,
        }
    df = DataFrame(catDict).T
    df["MAG"] = magnitudes
    df = df.replace({"None": nan})
    with open(outputFile, "w") as f:
        df.to_string(f, index=False)
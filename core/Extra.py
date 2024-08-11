import os
import sys
import time
import warnings
from pathlib import Path

from numpy import array, diff, max, mean, nan, round_, sqrt
from obspy import UTCDateTime as utc
from obspy import read_events
from obspy.core.event import Catalog
from obspy.geodetics.base import degrees2kilometers as d2k
from obspy.geodetics.base import gps2dist_azimuth as gps
from obspy.geodetics.base import kilometers2degrees as k2d
from pandas import DataFrame, Series, read_csv, to_datetime
from yaml import SafeLoader, load

warnings.filterwarnings("ignore")


def logger(message, mode="a"):
    """
    Function for generating logs.

    Parameters
    ----------
    message : str
        message to be loged.
    mode : str, optional
        loging mode. The default is "a".

    Returns
    -------
    None.

    """
    Path("results").mkdir(parents=True, exist_ok=True)
    logPath = os.path.join("results", "running.log")
    message = time.strftime("%d %b %Y %H:%M:%S - ") + message + "\n"
    with open(logPath, mode) as f:
        f.write(message)


def readConfiguration():
    """
    Read configuration file

    Returns
    -------
    config : dict
        configuration parameters.

    """
    if not os.path.exists("config.yml"):
        msg = "+++ Could not find configuration file! Aborting ..."
        print(msg)
        logger(msg, mode="w")
        sys.exit()
    with open("config.yml") as f:
        config = load(f, Loader=SafeLoader)
    msg = "+++ Configuration file was loaded successfully ..."
    print(msg)
    logger(msg, mode="w")
    return config


def handleNone(value, degree=False, dtype="float"):
    """Handle missing values

    Args:
        value (float): a float value
        degree (bool, optional): whether convert to degree or not.
        Defaults to False.

    Returns:
        float: handled value
    """
    if value is None:
        return nan
    else:
        if degree:
            return d2k(value)
        return int(value) if dtype == "int" else value


def getRMS(arrivals):
    time_residuals = array([
        arrival.time_residual for arrival in arrivals if isinstance(
            arrival.time_residual, float)
    ])
    time_weights = array([
        arrival.time_weight if isinstance(arrival.time_weight, float) else nan for arrival in arrivals
    ])
    if time_residuals.size:
        weighted_rms = sum(time_weights * time_residuals **
                           2) / sum(time_weights)
        weighted_rms = sqrt(weighted_rms)
    else:
        weighted_rms = nan
    return weighted_rms


def getHer(event):
    """Get horizontal error of event

    Args:
        event (obspy.event): an obspy event

    Returns:
        float: event horizontal error
    """
    try:
        x = event.origins[0].latitude_errors.uncertainty
        y = event.origins[0].longitude_errors.uncertainty
        return round(d2k(sqrt(x**2 + y**2)), 1)
    except TypeError:
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


def readHypoddConfig():
    hypoddConfigPath = os.path.join("files", "hypodd.yml")
    if not os.path.exists(hypoddConfigPath):
        msg = "+++ Could not find hypoDD configuration file! Aborting ..."
        print(msg)
        sys.exit()
    with open(hypoddConfigPath) as f:
        config = load(f, Loader=SafeLoader)
    msg = "+++ HypoDD Configuration file was loaded successfully ..."
    print(msg)
    return config


def loadHypoDDRelocFile():
    names = ["ID",  "LAT",  "LON",  "DEPTH",
             "X",  "Y",  "Z",
             "EX",  "EY",  "EZ",
             "YR",  "MO",  "DY",  "HR",  "MI",  "SC",
             "MAG",
             "NCCP",  "NCCS",
             "NCTP",  "NCTS",
             "RCC",  "RCT",
             "CID "]
    hypodd_df = read_csv("hypoDD.reloc", delim_whitespace=True, names=names)
    hypodd_df.sort_values(by=["ID"], inplace=True)
    hypodd_df.set_index(["ID"], inplace=True, drop=False)
    return hypodd_df


def loadVelocityFile(config):
    Pvel = config["VelocityModel"]["Pvel"]
    Deps = config["VelocityModel"]["Deps"]
    VpVs = config["VelocityModel"]["VpVs"]
    VpVs = [VpVs] * len(Pvel)
    velocity_df = DataFrame({"vp": Pvel, "depth": Deps, "vpvs": VpVs})
    return velocity_df


def writexyzm(outName, nEvents):
    hypodd_df = loadHypoDDRelocFile()
    outputFile = f"xyzm_{outName}.dat"
    hypodd_df["year"] = hypodd_df.YR
    hypodd_df["month"] = hypodd_df.MO.replace(0, 1)
    hypodd_df["day"] = hypodd_df.DY.replace(0, 1)
    hypodd_df["hour"] = hypodd_df.HR
    hypodd_df["minute"] = hypodd_df.MI
    hypodd_df["second"] = hypodd_df.SC
    hypodd_df["ORT"] = to_datetime(hypodd_df[["year",
                                              "month",
                                              "day",
                                              "hour",
                                              "minute",
                                              "second"]])
    hypodd_df["ORT"] = hypodd_df["ORT"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    hypodd_df["Lon"] = hypodd_df.LON
    hypodd_df["Lat"] = hypodd_df.LAT
    hypodd_df["Dep"] = hypodd_df.DEPTH
    hypodd_df["Mag"] = hypodd_df.MAG
    hypodd_df["Nus"] = hypodd_df.NCTP
    hypodd_df["NuP"] = hypodd_df.NCTP
    hypodd_df["NuS"] = hypodd_df.NCTS
    hypodd_df["ADS"] = nan
    hypodd_df["MDS"] = nan
    hypodd_df["GAP"] = nan
    hypodd_df["RMS"] = hypodd_df.RCT
    hypodd_df["ERH"] = sqrt(hypodd_df.EX**2 + hypodd_df.EY**2)*1e-3
    hypodd_df["ERZ"] = hypodd_df.EZ*1e-3
    columns = ["ORT", "Lon", "Lat", "Dep", "Mag",
               "Nus", "NuP", "NuS", "ADS", "MDS", "GAP", "RMS", "ERH", "ERZ"]
    UnLocatedEventsID = set(range(1, nEvents+1)) - set(hypodd_df.index)
    for i in UnLocatedEventsID:
        hypodd_df.loc[i] = nan
    hypodd_df.sort_index(inplace=True)
    with open(outputFile, "w") as f:
        hypodd_df.to_string(f, columns=columns, index=False, formatters={
            "ORT": "{:}".format,
            "Lon": "{:7.3f}".format,
            "Lat": "{:7.3f}".format,
            "Dep": "{:7.3f}".format,
            "Mag": "{:4.1f}".format,
            "Nus": "{:3.0f}".format,
            "NuP": "{:3.0f}".format,
            "NuS": "{:3.0f}".format,
            "ADS": "{:5.1f}".format,
            "MDS": "{:5.1f}".format,
            "GAP": "{:3.0f}".format,
            "RMS": "{:5.2f}".format,
            "ERH": "{:7.3f}".format,
            "ERZ": "{:7.3f}".format,
        })
    return len(hypodd_df)


def computeExtraInfo(evLon, evLat, arrivals, stationFile):
    station_df = read_csv(stationFile)
    codes = list(set([arv.pick.waveform_id.station_code for arv in arrivals]))
    station_df = station_df[station_df.code.isin(codes)]
    station_df[["Dist"]] = station_df.apply(
        lambda x: Series(gps(evLat, evLon, x.lat, x.lon)[0]*1e-3), axis=1)
    station_df[["Azim"]] = station_df.apply(
        lambda x: Series(gps(evLat, evLon, x.lat, x.lon)[1]), axis=1)
    azimuths = station_df["Azim"].sort_values()
    Dist = station_df["Dist"]
    GAP = int(max(diff(azimuths)))
    ADS = Dist.mean()
    MDS = Dist.min()
    Nus = len(codes)
    NuP = sum([1 for arrival in arrivals if arrival.phase.upper().startswith("P")])
    NuS = sum([1 for arrival in arrivals if arrival.phase.upper().startswith("S")])
    return Nus, NuP, NuS, ADS, MDS, GAP


def hypoDD2nordic(outName, stationFile):
    print(f"+++ Reading & Updating catalog for {outName} ...")
    catalog = read_events(f"{outName}.out")
    hypodd_df = read_csv(f"xyzm_{outName}.dat", delim_whitespace=True)
    hypodd_df_out = hypodd_df.copy()
    hypodd_df.ERH = k2d(hypodd_df.ERH)
    hypodd_df.ERZ = hypodd_df.ERZ*1e3
    hypodd_df.replace(nan, None, inplace=True)
    outCatalog = Catalog()
    for r, row in hypodd_df.iterrows():
        event = catalog[r]
        preferred_origin = event.preferred_origin()
        arrivals = preferred_origin.arrivals
        picks = {pick.resource_id: pick for pick in event.picks}
        for arrival in arrivals:
            arrival.update({"pick": picks[arrival.pick_id]})
        eOrt = utc(row.ORT) if row.ORT else preferred_origin.time
        eLat = row.Lat
        erLat = row.ERH
        eLon = row.Lon
        erLon = row.ERH
        eDep = row.Dep
        erDep = row.ERZ
        preferred_origin.time = eOrt
        preferred_origin.latitude = eLat
        preferred_origin.longitude = eLon
        preferred_origin.depth = eDep*1e3 if eDep else None
        preferred_origin.latitude_errors.uncertainty = erLat
        preferred_origin.longitude_errors.uncertainty = erLon
        preferred_origin.depth_errors.uncertainty = erDep
        preferred_origin.quality.azimuthal_gap = row.GAP
        if None not in [eLat, eLon]:
            Nus, NuP, NuS, ADS, MDS, GAP = computeExtraInfo(
                eLon, eLat, arrivals, stationFile)
            hypodd_df_out.loc[r, "Nus"] = Nus
            hypodd_df_out.loc[r, "NuP"] = NuP
            hypodd_df_out.loc[r, "NuS"] = NuS
            hypodd_df_out.loc[r, "ADS"] = ADS
            hypodd_df_out.loc[r, "MDS"] = MDS
            hypodd_df_out.loc[r, "GAP"] = GAP
            preferred_origin.quality.azimuthal_gap = GAP
        outCatalog.append(event)
        hypodd_df_out.loc[r, "Mag"] = row.Mag
    outCatalog.write(f"{outName}_hypodd.out",
                     format="nordic", high_accuracy=False)
    columns = ["ORT", "Lon", "Lat", "Dep", "Mag",
               "Nus", "NuP", "NuS", "ADS", "MDS", "GAP", "RMS", "ERH", "ERZ"]
    with open(f"xyzm_{outName}_hypodd.dat", "w") as f:
        hypodd_df_out.to_string(
            f, columns=columns, index=False, formatters={
                "ORT": "{:}".format,
                "Lon": "{:7.3f}".format,
                "Lat": "{:7.3f}".format,
                "Dep": "{:7.3f}".format,
                "Mag": "{:4.1f}".format,
                "Nus": "{:3.0f}".format,
                "NuP": "{:3.0f}".format,
                "NuS": "{:3.0f}".format,
                "ADS": "{:5.1f}".format,
                "MDS": "{:5.1f}".format,
                "GAP": "{:3.0f}".format,
                "RMS": "{:5.2f}".format,
                "ERH": "{:7.3f}".format,
                "ERZ": "{:7.3f}".format,
            })


def catalog2xyzm(hypInp, outName):
    """Convert catalog to xyzm file format

    Args:
        hypInp (str): file name of NORDIC file
        catalogFileName (str): file name of xyzm.dat file
    """
    cat = read_events(hypInp)
    outputFile = f"xyzm_{outName:s}_initial.dat"
    catDict = {}
    for i, event in enumerate(cat):
        preferred_origin = event.preferred_origin()
        preferred_magnitude = event.preferred_magnitude()
        arrivals = preferred_origin.arrivals
        ort = preferred_origin.time
        lat = preferred_origin.latitude
        lon = preferred_origin.longitude
        mag = preferred_magnitude.mag if preferred_magnitude else nan
        try:
            dep = preferred_origin.depth*0.001
        except TypeError:
            dep = nan
        try:
            nus = handleNone(
                preferred_origin.quality.used_station_count, dtype="int")
        except AttributeError:
            nus = nan
        nuP = len(
            [arrival.phase for arrival in arrivals if "P" in arrival.phase.upper()])
        nuS = len(
            [arrival.phase for arrival in arrivals if "S" in arrival.phase.upper()])
        mds = handleNone(
            min([handleNone(arrival.distance) for arrival in preferred_origin.arrivals]), degree=True)
        ads = round_(handleNone(
            mean([handleNone(arrival.distance) for arrival in preferred_origin.arrivals]), degree=True), 2)
        try:
            gap = handleNone(
                preferred_origin.quality.azimuthal_gap, dtype="int")
        except AttributeError:
            gap = nan
        rms = getRMS(preferred_origin.arrivals)
        erh = getHer(event)
        erz = getZer(event)
        catDict[i] = {
            "ORT": ort,
            "Lon": lon,
            "Lat": lat,
            "Dep": dep,
            "Mag": mag,
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
    df = df.replace({"None": nan})
    with open(outputFile, "w") as f:
        df.to_string(f, index=False, formatters={
            "ORT": "{:}".format,
            "Lon": "{:7.3f}".format,
            "Lat": "{:7.3f}".format,
            "Dep": "{:7.3f}".format,
            "Mag": "{:4.1f}".format,
            "Nus": "{:3.0f}".format,
            "NuP": "{:3.0f}".format,
            "NuS": "{:3.0f}".format,
            "ADS": "{:5.1f}".format,
            "MDS": "{:5.1f}".format,
            "GAP": "{:3.0f}".format,
            "RMS": "{:5.2f}".format,
            "ERH": "{:5.1f}".format,
            "ERZ": "{:5.1f}".format,
        })


def loadxyzm(*xyzmPaths):
    reports = []
    for xyzmPath in xyzmPaths:
        report = read_csv(xyzmPath, delim_whitespace=True)
        reports.append(report)
    return reports

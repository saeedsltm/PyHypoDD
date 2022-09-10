import os
import random
from math import sqrt

from numpy import nan, sqrt
from obspy.geodetics.base import degrees2kilometers as d2k
from obspy.geodetics.base import locations2degrees as l2d
from obspy import read_events
from pandas import read_json, read_csv
from LatLon import lat_lon as ll


def parseStationFile(stationFileName):
    """Parse station information

    Args:
        stationFileName (str): path to station file

    Returns:
        DataFrame: a panda dataFrame contains stations information
    """
    stationsDF = read_csv(stationFileName, delim_whitespace=True)
    return stationsDF


def parseVelocityFile(velocityFileName):
    """Parse velocity information

    Args:
        velocityFileName (str): path to velocity file

    Returns:
        DataFrame: a panda dataFrame contains velocity information
    """
    velocityDF = read_csv(velocityFileName, delim_whitespace=True)
    return velocityDF


def distanceDiff(xa, xb, ya, yb):
    """Compute distance between two lists of points

    Args:
        xa (array): X array of first points
        xb (array): X array of second points
        ya (array): Y array of first points
        yb (array): Y array of second points

    Returns:
        array: an array contains distances
    """
    return sqrt((xa-xb)**2+(ya-yb)**2)


def handleNone(value, degree=False):
    """Handle missing values

    Args:
        value (float): a float value
        degree (bool, optional): whether convert to degree or not. Defaults to False.

    Returns:
        float: handled value
    """
    if value == None:
        return nan
    else:
        if degree:
            return d2k(value)
        return value


def setTTError(tt, m=0.0, std=1.0, addWeight=True):
    """Setting error with travel times.

    Args:
        tt (float): travel time of a desired phase (s),
        m (float, optional): mean of Gaussian distribution for generating error. Defaults to 0.0.
        std (float, optional): standard deviation of Gaussian distribution for generating error. Defaults to 1.0.
        addWeight (bool, optional): whether or not considering weights.. Defaults to True.

    Returns:
        tuple: a tuple contains of travel time and associated weight.
    """
    err = random.gauss(m, std)
    tt = tt + err
    w = mapError2Weight(err, m, std, addWeight)
    return tt, w


def mapError2Weight(err, m, std, addWeight=True):
    """Map error magnitude to weight.

    Args:
        err (float): error in seconds,
        m (float, optional): mean of Gaussian distribution for generating error. Defaults to 0.0.
        std (float, optional): standard deviation of Gaussian distribution for generating error. Defaults to 1.0.
        addWeight (bool, optional): whether or not considering weights.. Defaults to True.

    Returns:
        _type_: _description_
    """
    if not addWeight:
        return 0
    absError = abs(abs(std) - abs(m))
    if 0.0 <= abs(err) <= 0.25*absError:
        return 0
    elif 0.25*absError <= abs(err) <= 0.5*absError:
        return 1
    elif 0.5*absError <= abs(err) <= 0.75*absError:
        return 2
    elif 0.75*absError <= abs(err) <= 1.0*absError:
        return 3
    elif 1.0*absError <= abs(err):
        return 4


def ReadExtra(pick, reverseWeight=False):
    """reads pick weight 

    Args:
        pick (event.pick): an obspy event pick
        reverseWeight (bool, optional): if True reverse weight is on [4>0, 3>0.75... 0>1.0]. Defaults to False.


    Returns:
        str: pick weight
    """
    mapWeight = {"4": 0, "3": 0.25, "2": 0.5, "1": 0.75, "0": 1.0}
    try:
        weight = pick.extra.get("nordic_pick_weight", "0")
        if isinstance(weight, type({})):
            weight = weight["value"]
    except AttributeError:
        weight = "0"
    return float(mapWeight[weight])


def hypoDD2xyzm():
    """Convert HypoDD "reloc" file to "xyzm"
    """
    originEvents = read_json("events.json")
    columns = ["ID", "LAT", "LON", "DEPTH", "X", "Y", "Z", "EX", "EY", "EZ", "YYYY", "MM",
               "DD", "HH", "MN", "SEC", "MAG", "NCCP", "NCCS", "NPHASEP", "NPHASES", "RCC", "RMS", "CID"]
    relocatedEvents = read_csv(
        "hypoDD.reloc", delim_whitespace=True, names=columns)
    relocatedEvents["SEH"] = 0.5*(sqrt(relocatedEvents["EX"]**2 + relocatedEvents["EY"]**2)*1e-3)  # type: ignore
    relocatedEvents["SEZ"] = relocatedEvents["EZ"]*1e-3
    for i, relocatedEvent in relocatedEvents.iterrows():
        mag = originEvents[originEvents.ID == relocatedEvent.ID].Mag
        relocatedEvents.at[i, "MAG"] = mag
    for NoneValue in ["NSTUSED", "MIND", "GAP"]:
        relocatedEvents[NoneValue] = nan
    formatters = {
        "LON": '{:.3f}'.format,
        "LAT": '{:.3f}'.format,
        "DEPTH": '{:.1f}'.format,
        "MAG": '{:.1f}'.format,
        "NSTUSED": '{:.0f}'.format,
        "NPHASEP": '{:.0f}'.format,
        "NPHASES": '{:.0f}'.format,
        "MIND": '{:.0f}'.format,
        "GAP": '{:.0f}'.format,
        "RMS": '{:.2f}'.format,
        "SEH": '{:.1f}'.format,
        "SEZ": '{:.1f}'.format,
        "YYYY": '{:.0f}'.format,
        "MM": '{:.0f}'.format,
        "DD": '{:.0f}'.format,
        "HH": '{:.0f}'.format,
        "MN": '{:.0f}'.format,
        "SEC": '{:.2f}'.format
    }
    columns = list(formatters.keys())
    relocatedEvents.to_string(
        "xyzm.dat", columns=columns, formatters=formatters, index=False)
    cmd = "rm *reloc*"
    os.system(cmd)

def writeStationZeroFile(stationFile, stationsDF, velocityDF, trialDepth, xNear, xFar):
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
    stationCodes = stationsDF["CODE"].values
    stationLats = stationsDF["LAT"].values
    stationLons = stationsDF["LON"].values
    stationElvs = stationsDF["ELV"].values
    Vp = velocityDF["Vp"].values
    Z = velocityDF["Z"].values
    VpVs = velocityDF["VpVs"].values[0]
    Interfaces = velocityDF["Interface"].values
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
                v=v, z=z, i=i.replace("nul", " ")
            ))
        g.write("\n")
        g.write("{trialDepth:4.0f}.{xNear:4.0f}.{xFar:4.0f}. {VpVs:4.2f}".format(
            trialDepth=trialDepth, xNear=xNear, xFar=xFar, VpVs=VpVs
        ))
        g.write("\nNew")
    return stationFile    

def filterStations(configs, stationsDF):
    """Check if station is within the radius defined by user

    Args:
        configs (dict): a dictionary contains configurations
        stationsDF (DataFrame): a pandas data-frame contains stations information

    Returns:
        DataFrame: a pandas data-frame contains filtered stations
    """
    criticalDistance = d2k(l2d(configs["Region"]["CentralLat"], configs["Region"]["CentralLon"],
                               stationsDF.LAT.values, stationsDF.LON.values))
    stationsDF = stationsDF[criticalDistance <= configs["PH2DT"]["MAXDIST"]]
    return stationsDF

def OutputCatalog(configs, resultDir, CatalogFile, stationsDF, velocityDF):
    catalog = read_events(CatalogFile)
    outputCatalogName = os.path.join(resultDir, "selectedCatalog.dat")
    catalog.write(outputCatalogName, format=configs["OutPut"]["CatalogFormat"])
    if configs["OutPut"]["CatalogFormat"] == "NORDIC":
        outputStationName = os.path.join(resultDir, "STATION0.HYP")
        writeStationZeroFile(outputStationName, stationsDF, velocityDF, trialDepth=10, xNear=50, xFar=500)


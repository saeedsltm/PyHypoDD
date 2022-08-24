import os
import random
import time
from datetime import datetime as dt
from datetime import timedelta as td
from math import sqrt

from numpy import nan, sqrt
from obspy.geodetics.base import degrees2kilometers as d2k
from obspy.geodetics.base import locations2degrees as l2d
from pandas import read_json, read_csv, to_datetime

def parseStationFile(stationFileName):
    stationsDF = read_csv(stationFileName, delim_whitespace=True)
    return stationsDF

def parseVelocityFile(velocityFileName):
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
    events = read_json("events.json")
    ots = to_datetime(events["OT"])
    with open("hypoDD.reloc") as f, open("xyzm.dat", "w") as g:
        header = "     LON     LAT   DEPTH     MAG    PHUSD   NO_ST   MIND     GAP     RMS     SEH     SEZ  YYYY MM DD HH MN SEC"
        g.write(header+"\n")
        for l in f:
            if "NaN" in l.split():
                continue
            _, lat, lon, dep, _, _, _, ex, ey, ez, yr, mo, dy, hr, mn, sc, mag, _, _, _, _, _, rms, _ = [
                float(_) for _ in l.split()]
            if sc == 60.0:
                sc = 59.99
            yr, mo, dy, hr, mn = [int(_) for _ in [yr, mo, dy, hr, mn]]
            msc = int((sc - int(sc))*1e6)
            sc = int(sc)
            ort = dt(yr, mo, dy, hr, mn, sc, msc)
            mag = events[abs(ots-ort) < td(seconds=5)  # type: ignore
                         ].Mag.values[0]
            if not mag:
                mag = 0.0
            ort = ort.strftime("  %Y %m %d %H %M %S.%f")[:24]
            nop = 99
            nst = 99
            mds = 99
            gap = 360
            ex, ey, ez = ex/1000., ey/1000., ez/1000.
            seh = sqrt((ex)**2 + (ey)**2)
            sez = ez
            fmt = '%8.3f%8.3f%8.1f%8.1f%8d%8d%8.1f%8d%8.3f%8.1f%8.1f'
            fmt = fmt + ort+'\n'
            # lon lat dep mag phusd no_st mds gap rms seh sez
            g.write(fmt % (lon, lat, dep, mag, nop,
                    nst, mds, gap, rms, seh, sez))
    cmd = "rm *reloc*"
    os.system(cmd)


def checkStationDist(configs, stationsDF):
    """Check if station is within the radius defined by user

    Args:
        configs (dict): a dictionary contains configurations
        staLat (float): latitude of station
        staLon (float): longitude of station

    Returns:
        bool: True if the station lies within the pre-defined radius
    """

    d = d2k(l2d(configs["Region"]["CentralLat"], configs["Region"]["CentralLon"],
            stationsDF.LAT.values, stationsDF.LON.values))
    stationsDF = stationsDF[d <= configs["PH2DT"]["MAXDIST"]]
    return stationsDF
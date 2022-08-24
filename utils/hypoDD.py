import os
from json import dumps
from shutil import copy

from utils.extra import ReadExtra
from utils.nordic import readVelocityFile


def writeHeader(preferred_origin, eventNumber, fo):
    """Write header in HypoDD phase file

    Args:
        preferred_origin (event.origin): preferred origin of the event
        eventNumber (int): number of the event
        fo (file): an open file ready to be updated
    """
    header = "# {year:4d} {month:02d} {day:02d} {hour:002d} {minute:002d} {seconds:05.2f}  {lat:6.3f}  {lon:6.3f}  {dep:4.1f} 0.00 0.00 0.00 {rms:4.1f} {eventID:9d}\n".format(
        year=preferred_origin.time.year,
        month=preferred_origin.time.month,
        day=preferred_origin.time.day,
        hour=preferred_origin.time.hour,
        minute=preferred_origin.time.minute,
        seconds=preferred_origin.time.second + preferred_origin.time.microsecond*1e-6,
        lat=preferred_origin.latitude,
        lon=preferred_origin.longitude,
        dep=preferred_origin.depth*1e-3,
        rms=preferred_origin.quality.standard_error,
        eventID=int(eventNumber+1e5))
    fo.write(header)

def writeEventsInfo(eventInfoDict, resultDir):
    """Write event information

    Args:
        eventInfoDict (dict): a dictionary contains event information
        resultDir (str): path to the HypoDD input files
    """
    json_eventsID = dumps(eventInfoDict, indent=4)
    with open(os.path.join(resultDir, "events.json"), "w") as fo:
        fo.write(json_eventsID)


def writePhase(preferred_origin, picks, arrivals, usedStation, fo):
    """Write phase readings in HypoDD format

    Args:
        preferred_origin (event.origin): preferred origin of the event
        picks (event.picks): an obspy event picks object
        arrivals (event.origin.arrivals): an obspy event origin arrivals object
        usedStation (list): a list contains of used stations
        fo (file): an open file ready to be updated
    """
    picksSta = {pick.resource_id: pick.waveform_id.station_code for pick in picks}
    picksTime = {pick.resource_id: pick.time-preferred_origin.time for pick in picks}
    picksWeight = {pick.resource_id: ReadExtra(pick) for pick in picks}
    for arrival in arrivals:
        if picksSta[arrival.pick_id] in usedStation:
            fo.write("{station:4s} {time:7.3f}  {wet:4.1f} {phase:1s}\n".format(
                station=picksSta[arrival.pick_id],
                time=picksTime[arrival.pick_id],
                wet=picksWeight[arrival.pick_id],
                phase=arrival.phase[0]))

def writeStation(stationCode, lat, lon, elv, fo):
    """Write stations info in HypoDD format

    Args:
        stationCode (str): station code
        lat (float): station latitude
        lon (float): station longitude
        elv (float): station elevation
        fo (file): an open file ready to be updated
    """
    fo.write("{stationCode:4s}     {lat:06.3f}  {lon:06.3f} {elv:00004.0f}\n".format(
        stationCode=stationCode,
        lat=lat,
        lon=lon,
        elv=elv))

def prepareHypoDDConfigFiles(configs, velocityDF, resultDir):
    """Prepare HypoDD configuration files

    Args:
        configs (dict): a dictionary contains configurations
        resultDir (str): path to HypoDD input files
    """
    copy(os.path.join("files", "ph2dt.inp"), resultDir)
    Vp = velocityDF["Vp"]
    Z = velocityDF["Z"]
    VpVs = velocityDF["VpVs"][0]
    with open(os.path.join("files", "ph2dt.inp")) as fo, open(os.path.join(resultDir, "ph2dt.inp"), "w") as go:
        for l in fo:
            if "MINWGHT MAXDIST MAXSEP MAXNGH MINLNKS MINOBS MAXOBS" in l:
                go.write(l)
                next(fo)
                l = "  ".join(map(str, configs["PH2DT"].values()))
            go.write(l)
    with open(os.path.join("files", "hypoDD.inp")) as fo, open(os.path.join(resultDir, "hypoDD.inp"), "w") as go:
        for l in fo:
            if "* IDAT   IPHA   DIST" in l and ":" not in l:
                go.write(l)
                next(fo)
                l = "    2     3    {MAXDIST}\n".format(MAXDIST=configs["HypoDD"]["DIST"])
            elif "* OBSCC  OBSCT" in l and ":" not in l:
                go.write(l)
                next(fo)                
                l = "     0    {OBSCT}\n".format(OBSCT=configs["HypoDD"]["OBSCT"])                
            elif "* NLAY  RATIO" in l and ":" not in l:
                go.write(l)
                next(fo)
                l = "   {NLAY}      {VpVs}\n".format(NLAY=len(Z), VpVs=VpVs)
            elif "* TOP" in l and ":" not in l:
                go.write(l)
                next(fo)
                l = " ".join(map(str, Z))+"\n"
            elif "* VEL" in l and ":" not in l:
                go.write(l)
                next(fo)
                l = " ".join(map(str, Vp))+"\n"
            go.write(l)




    
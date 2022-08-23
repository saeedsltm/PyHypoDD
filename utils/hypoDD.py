
from json import dumps
from utils.extra import ReadExtra
import os
from shutil import copy
from utils.nordic import readVelocityFile


def writeHeader(preferred_origin, eventNumber, fo):
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

def writeEventsID(eventIDDict, hypoDDInputs):
    json_eventsID = dumps(eventIDDict, indent=4)
    with open(os.path.join(hypoDDInputs, "events.json"), "w") as fo:
        fo.write(json_eventsID)


def writePhase(preferred_origin, picks, arrivals, usedStation, fo):
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
    fo.write("{stationCode:4s}     {lat:06.3f}  {lon:06.3f} {elv:00004.0f}\n".format(
        stationCode=stationCode,
        lat=lat,
        lon=lon,
        elv=elv))

def prepareHypoDDConfigFiles(configs, hypoDDInputs):
    copy(os.path.join("files", "ph2dt.inp"), hypoDDInputs)
    velocityFileDict = readVelocityFile(configs["InputStationFileName"])
    Vp = velocityFileDict["Vp"]
    Z = velocityFileDict["Z"]
    VpVs = velocityFileDict["VpVs"]
    MINWGHT = configs["MINWGHT"]
    MAXDIST = configs["MAXDIST"]
    MAXSEP = configs["MAXSEP"]
    MAXNGH = configs["MAXNGH"]
    MINLNKS = configs["MINLNKS"]
    MINOBS = configs["MINOBS"]
    MAXOBS = configs["MAXOBS"]
    OBSCT = configs["OBSCT"]
    with open(os.path.join("files", "ph2dt.inp")) as fo, open(os.path.join(hypoDDInputs, "ph2dt.inp"), "w") as go:
        for l in fo:
            if "MINWGHT MAXDIST MAXSEP MAXNGH MINLNKS MINOBS MAXOBS" in l:
                go.write(l)
                next(fo)
                l = "  ".join(map(str, [MINWGHT, MAXDIST, MAXSEP, MAXNGH, MINLNKS, MINOBS, MAXOBS]))
            go.write(l)
    with open(os.path.join("files", "hypoDD.inp")) as fo, open(os.path.join(hypoDDInputs, "hypoDD.inp"), "w") as go:
        for l in fo:
            if "* IDAT   IPHA   DIST" in l and ":" not in l:
                go.write(l)
                next(fo)
                l = "    2     3    {MAXDIST}\n".format(MAXDIST=MAXDIST)
            elif "* OBSCC  OBSCT" in l and ":" not in l:
                go.write(l)
                next(fo)                
                l = "     0    {OBSCT}\n".format(OBSCT=OBSCT)                
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




    
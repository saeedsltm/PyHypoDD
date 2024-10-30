import os
from glob import glob
from pathlib import Path
from shutil import copy
from time import time

from core.Extra import (catalog2xyzm, hypoDD2nordic, loadVelocityFile, logger,
                        readHypoddConfig, hypoddReloc2xyzm, mergeDFs)
from core.Input import prepareHypoddInputs
from obspy import read_events

def locateHypoDD(config):
    hypoddConfig = readHypoddConfig()
    outName = f"{config['Region']['RegionName']}"
    stationPath = os.path.join("stations", "usedStations.csv")
    stationFile = os.path.abspath(stationPath)
    locationPath = os.path.abspath("results")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    velocity_df = loadVelocityFile(config)
    catalogFile = config["Files"]["InputCatalogFileName"]
    copy(catalogFile, os.path.join(locationPath, f"{outName}.out"))
    root = os.getcwd()
    os.chdir(locationPath)
    catalog = read_events(f"{outName}.out")
    maxAllowdedEventsPerChunk = 6e3
    nEvents = len(catalog)
    nChunks = int(nEvents//maxAllowdedEventsPerChunk)
    for nChunk in range(nChunks+1):
        print(f"+++ Relocating chunk {nChunk+1} ...")
        s = int(nChunk*maxAllowdedEventsPerChunk)
        if nChunk != nChunks:
            e = int((nChunk+1)*maxAllowdedEventsPerChunk)
            selectedCatalog = catalog[s:e]
        else:
            selectedCatalog = catalog[s:]
        nEvents = len(selectedCatalog)
        chunkPath = os.path.join(f"chunk_{nChunk+1}")
        Path(chunkPath).mkdir(parents=True, exist_ok=True)
        os.chdir(chunkPath)
        prepareHypoddInputs(config,
                            hypoddConfig,
                            selectedCatalog,
                            stationFile,
                            velocity_df,
                            locationPath)
        cmd = "ph2dt ph2dt.inp >/dev/null 2>/dev/null"
        os.system(cmd)
        cmd = "hypoDD hypoDD.inp >/dev/null 2>/dev/null"
        st = time()
        os.system(cmd)
        et = time()
        print("+++ Making summary files ...")
        nEvents = hypoddReloc2xyzm(nEvents, outName)
        hypoDD2nordic(selectedCatalog, stationFile, outName)
        for f in glob("hypoDD.reloc*"):
            os.remove(f)
        catalog2xyzm(selectedCatalog, outName)
        os.chdir(locationPath)
    mergeDFs(nChunks, outName)
    os.chdir(root)
    logger(f"Processing time for relocating {nEvents} events using HypoDD is: \
{et-st:.3f} s")
import os
from glob import glob
from pathlib import Path
from shutil import copy
from time import time

from core.Extra import (catalog2xyzm, hypoDD2nordic, loadVelocityFile, logger,
                        readHypoddConfig, writexyzm)
from core.Input import prepareHypoddInputs


def locateHypoDD(config):
    hypoddConfig = readHypoddConfig()
    outName = f"{config['Region']['RegionName']}"
    stationPath = os.path.join("stations", "stations.csv")
    stationFile = os.path.abspath(stationPath)
    locationPath = os.path.join("results")
    Path(locationPath).mkdir(parents=True, exist_ok=True)
    velocity_df = loadVelocityFile(config)
    catalogFile = config["Files"]["InputCatalogFileName"]
    copy(catalogFile, os.path.join(locationPath, f"{outName}.out"))
    root = os.getcwd()
    os.chdir(locationPath)
    nEvents = prepareHypoddInputs(config,
                                  hypoddConfig,
                                  outName,
                                  stationFile,
                                  velocity_df,
                                  locationPath)
    print("+++ Relocating ...")
    cmd = "ph2dt ph2dt.inp >/dev/null 2>/dev/null"
    os.system(cmd)
    cmd = "hypoDD hypoDD.inp >/dev/null 2>/dev/null"
    st = time()
    os.system(cmd)
    et = time()
    print("+++ Making summary files ...")
    nEvents = writexyzm(outName, nEvents)
    hypoDD2nordic(outName, stationFile)
    for f in glob("hypoDD.reloc*"):
        os.remove(f)
    catalog2xyzm(f"{outName}.out", outName)
    os.chdir(root)
    logger(f"Processing time for relocating {nEvents} events using HypoDD is: \
{et-st:.3f} s")

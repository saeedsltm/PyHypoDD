import os
import warnings
from datetime import datetime as dt
from pathlib import Path
from shutil import copy

from obspy import read_events
from tqdm import tqdm
from yaml import SafeLoader, load

from utils.extra import (filterStations, hypoDD2xyzm, parseStationFile,
                         parseVelocityFile)
from utils.hypoDD import (prepareHypoDDConfigFiles, writeEventsInfo,
                          writeHeader, writePhase)

warnings.filterwarnings("ignore")


class HypoDDRunner():
    def __init__(self):
        if not os.path.exists("config.yml"):
            copy(os.path.join("files", "defaultConfig.yml"), "config.yml")
        with open("config.yml") as f:
            self.configs = load(f, Loader=SafeLoader)
        self.resultDir = os.path.join(
            "results", dt.now().strftime("%y-%m-%dT%H%M%S"))
        self.createResultFolder()

    def createResultFolder(self):
        """Create folder for outputs
        """
        Path(self.resultDir).mkdir(parents=True, exist_ok=True)

    def prepareHypoDDInputFiles(self, catalog, stationsDF, velocityDF):
        """Prepare input files for running HypoDD
        """
        usedStationsDF = filterStations(self.configs, stationsDF)
        usedStationsDF.to_csv(os.path.join(
            self.resultDir, "station.dat"), sep=' ', index=False, header=False)
        eventsInfo = []
        with open(os.path.join(self.resultDir, "phase.dat"), "w") as f:
            for eventNumber, event in enumerate(tqdm(catalog)):
                picks = event.picks
                preferred_origin = event.preferred_origin()
                preferred_magnitude = event.preferred_magnitude()
                eventsInfo.append({
                    "OT": preferred_origin.time.strftime("%Y-%m-%d %H:%M:%S.%f"),
                    "Lat": preferred_origin.latitude,
                    "Lon": preferred_origin.longitude,
                    "Dep": preferred_origin.depth*1e-3,
                    "Mag": preferred_magnitude.mag,
                    "ErH": preferred_origin.longitude_errors.uncertainty,
                    "ErZ": preferred_origin.depth_errors.uncertainty,
                    "RMS": preferred_origin.quality.standard_error,
                    "ID": int(eventNumber+1e5),
                    "SMI": event.resource_id.id})
                arrivals = preferred_origin.arrivals
                writeHeader(preferred_origin, eventNumber, f)
                writePhase(preferred_origin, picks, arrivals,
                           usedStationsDF.CODE.values, f)
        writeEventsInfo(eventsInfo, self.resultDir)
        prepareHypoDDConfigFiles(self.configs, velocityDF, self.resultDir)

    def runHypoDD(self):
        """Run HypoDD programs
        """
        cmd01 = "ph2dt ph2dt.inp > /dev/null"
        cmd02 = "hypoDD hypoDD.inp > /dev/null"
        root = os.getcwd()
        os.chdir(self.resultDir)
        os.system(cmd01)
        os.system(cmd02)
        hypoDD2xyzm()
        os.chdir(root)


if __name__ == "__main__":
    app = HypoDDRunner()
    catalog = read_events(app.configs["IO"]["InputCatalogFileName"])
    stationsDF = parseStationFile(app.configs["IO"]["InputStationFileName"])
    velocityDF = parseVelocityFile(app.configs["IO"]["InputVelocityFileName"])
    app.prepareHypoDDInputFiles(catalog, stationsDF, velocityDF)
    app.runHypoDD()

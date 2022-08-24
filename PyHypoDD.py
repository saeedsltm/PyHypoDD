import os
from pathlib import Path
from datetime import datetime as dt

from obspy import read_events
from tqdm import tqdm
from yaml import SafeLoader, load

from utils.extra import parseStationFile, parseVelocityFile, filterStations, hypoDD2xyzm
from utils.hypoDD import (prepareHypoDDConfigFiles, writeEventsInfo, writeHeader,
                          writePhase)


class HypoDD():
    def __init__(self):
        with open("config.yml") as f:
            self.configs = load(f, Loader=SafeLoader)
        self.catalog = read_events(self.configs["IO"]["InputCatalogFileName"])
        self.stationsDF = parseStationFile(self.configs["IO"]["InputStationFileName"])
        self.velocityDF = parseVelocityFile(self.configs["IO"]["InputVelocityFileName"])
        self.resultDir = os.path.join("results", dt.now().strftime("%y-%m-%dT%H%M%S"))
        self.createResultFolder()

    def createResultFolder(self):
        """Create folder for outputs
        """
        Path(self.resultDir).mkdir(parents=True, exist_ok=True)

    def prepareHypoDDInputFiles(self):
        """Prepare input files for running HypoDD
        """
        usedStationsDF = filterStations(self.configs, self.stationsDF)
        usedStationsDF.to_csv(os.path.join(self.resultDir, "station.dat"), sep=' ', index=False, header=False)
        eventsInfo = []
        with open(os.path.join(self.resultDir, "phase.dat"), "w") as f:
            for eventNumber, event in enumerate(tqdm(self.catalog)):
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
                writePhase(preferred_origin, picks, arrivals, usedStationsDF.CODE.values, f)
        writeEventsInfo(eventsInfo, self.resultDir)
        prepareHypoDDConfigFiles(self.configs, self.velocityDF, self.resultDir)

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
    app = HypoDD()
    app.prepareHypoDDInputFiles()
    app.runHypoDD()

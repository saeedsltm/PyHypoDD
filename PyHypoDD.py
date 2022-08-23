from obspy import read_events
from tqdm import tqdm
from json import load
from utils.nordic import readStationFile
from utils.hypoDD import writeHeader, writePhase, writeEventsID, writeStation, prepareHypoDDConfigFiles
from utils.extra import hypoDD2xyzm, checkStationDist
from pathlib import Path
import os

class HypoDD():
    def __init__(self):
        with open("config.json") as f:
            self.configs = load(f)
        self.catalog = read_events(self.configs["InputCatalogFileName"])
        self.stationsDict = readStationFile(self.configs["InputStationFileName"])
        self.hypoDDInputs = "hypoDDInputs"
        self.createResultFolder()

    def createResultFolder(self):
        hypoDDInputs = os.path.join(self.hypoDDInputs)
        Path(hypoDDInputs).mkdir(parents=True, exist_ok=True)        

    def prepareHypoDDInputFiles(self):
        usedStation = []
        with open(os.path.join(self.hypoDDInputs, "station.dat"), "w") as fo:
            for station in self.stationsDict:
                stationCode=station
                staLat = self.stationsDict[station]["Lat"]
                staLon = self.stationsDict[station]["Lon"]
                staElv = self.stationsDict[station]["Elv"]
                if checkStationDist(self.configs, staLat, staLon):
                    usedStation.append(stationCode)
                    writeStation(stationCode, staLat, staLon, staElv, fo)        
        eventsInfo = []    
        with open(os.path.join(self.hypoDDInputs, "phase.dat"), "w") as f:
            for eventNumber,event in enumerate(tqdm(self.catalog)):
                picks = event.picks
                preferred_origin = event.preferred_origin()
                preferred_magnitude = event.preferred_magnitude()
                eventsInfo.append({
                    "OT":preferred_origin.time.strftime("%Y-%m-%d %H:%M:%S.%f"),
                    "Lat":preferred_origin.latitude,
                    "Lon":preferred_origin.longitude,
                    "Dep":preferred_origin.depth*1e-3,
                    "Mag":preferred_magnitude.mag,
                    "ErH":preferred_origin.longitude_errors.uncertainty,
                    "ErZ":preferred_origin.depth_errors.uncertainty,
                    "RMS":preferred_origin.quality.standard_error,
                    "ID":int(eventNumber+1e5),
                    "SMI":event.resource_id.id})
                arrivals = preferred_origin.arrivals
                writeHeader(preferred_origin, eventNumber, f)
                writePhase(preferred_origin, picks, arrivals, usedStation, f)
        writeEventsID(eventsInfo, self.hypoDDInputs)

        prepareHypoDDConfigFiles(self.configs, self.hypoDDInputs)
    
    def runHypoDD(self):
        cmd01 = "ph2dt ph2dt.inp > /dev/null"
        cmd02 = "hypoDD hypoDD.inp > /dev/null"
        root = os.getcwd()
        os.chdir(self.hypoDDInputs)
        os.system(cmd01)
        os.system(cmd02)
        hypoDD2xyzm()
        os.chdir(root)


if __name__ == "__main__":
    app = HypoDD()
    app.prepareHypoDDInputFiles()
    app.runHypoDD()
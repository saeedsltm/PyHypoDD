from core.Extra import readConfiguration
from core.Locate import locateHypoDD
from core.PrepareInputs import CreatInputStationFile, GetStationListFromCatalog
from core.Visulize import plotSeismicityMap


class Main():
    def __init__(self):
        self.config = readConfiguration()

    def prepareStations(self):
        GetStationListFromCatalog(self.config)
        CreatInputStationFile(self.config)

    def locate(self):
        locateHypoDD(self.config)

    def visulize(self):
        plotSeismicityMap(self.config)


if "__main__" == __name__:
    app = Main()
    app.prepareStations()
    app.locate()
    app.visulize()

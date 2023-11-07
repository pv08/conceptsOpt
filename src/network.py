from src.utils.functions import mkdir_if_not_exists
from src.utils.bus import BBus

class Network:
    def __init__(self):
        mkdir_if_not_exists('results/')
        self.h = 0.25
        self.delta = 0.01
        self.ref_bus_i = 0
        self.cost_deficit = 10000
        self.DGER = [{'BARRA': 0, 'CUSTO': 10, 'MAX': 30, 'MIN': 5, 'RAMPA': 10, 'CO2': 90},
                      {'BARRA': 1, 'CUSTO': 30, 'MAX': 40, 'MIN': 15, 'RAMPA': 5, 'CO2': 10},
                      {'BARRA': 2, 'CUSTO': 100, 'MAX': 40, 'MIN': 0, 'RAMPA': 3, 'CO2': 70}]



        self.DLIN =  [{'DE': 0, 'PARA': 1, 'SUSCEPTANCIA': 33, 'CONDUTANCIA': 25, 'LIMITES': 20},
                      {'DE': 0, 'PARA': 2, 'SUSCEPTANCIA': 50, 'CONDUTANCIA': 20, 'LIMITES': 25},
                      {'DE': 0, 'PARA': 2, 'SUSCEPTANCIA': 50, 'CONDUTANCIA': 20, 'LIMITES': 25},
                      {'DE': 1, 'PARA': 2, 'SUSCEPTANCIA': 50, 'CONDUTANCIA': 20, 'LIMITES': 30}]

        self.DEMANDA = [{'HORA': 0, 'BARRA0': 0, 'BARRA1': 40, 'BARRA2': 30},
                                 {'HORA': 1, 'BARRA0': 0, 'BARRA1': 43, 'BARRA2': 25},
                                 {'HORA': 2, 'BARRA0': 0, 'BARRA1': 25, 'BARRA2': 25}]

        self.bus = BBus(DBAR=self.DGER, DLIN=self.DLIN, ref_bar=self.ref_bus_i)

        self.tempDim = [value for value in range(len(self.DEMANDA))]
        self.gerDim = [value for value in range(len(self.DGER))]


    def minimize(self):
        raise NotImplementedError

    def obj_fn(self, model):
        raise NotImplementedError

    def saveModelDetails(self, name):
        raise NotImplementedError


    def ger_limits(self, model, i, t):
        return (self.DGER[i]['MIN'], self.DGER[i]['MAX'])

    def deficit_lim(self, model, i, t):
        return (0, self.DEMANDA[t][f'BARRA{i}'])



import pandas as pd
from src.utils.bus import BBus

class Network:
    def __init__(self):
        self.h = 0.25
        self.delta = 0.01
        self.base = 100
        self.ref_bus_i = 0
        self.cost_deficit = 10000
        self.DGER = [{'BARRA': 0, 'CUSTO': 10, 'MAX': 30, 'MIN': 5, 'RAMPA': 10, 'CO2': 90},
                      {'BARRA': 1, 'CUSTO': 30, 'MAX': 40, 'MIN': 15, 'RAMPA': 5, 'CO2': 10},
                      {'BARRA': 2, 'CUSTO': 100, 'MAX': 40, 'MIN': 0, 'RAMPA': 3, 'CO2': 70}]



        self.DLIN = [{'DE': 0, 'PARA': 1, 'SUSCEPTANCIA': 33/self.base, 'CONDUTANCIA': 25/self.base, 'LIMITES': 20/self.base},
                      {'DE': 0, 'PARA': 2, 'SUSCEPTANCIA': 50/self.base, 'CONDUTANCIA': 20/self.base, 'LIMITES': 25/self.base},
                      {'DE': 0, 'PARA': 2, 'SUSCEPTANCIA': 50/self.base, 'CONDUTANCIA': 20/self.base, 'LIMITES': 25/self.base},
                      {'DE': 1, 'PARA': 2, 'SUSCEPTANCIA': 50/self.base, 'CONDUTANCIA': 20/self.base, 'LIMITES': 30/self.base}]

        self.DEMANDA = [{'HORA': 0, 'BARRA0': 0, 'BARRA1': 40, 'BARRA2': 30},
                                 {'HORA': 1, 'BARRA0': 0, 'BARRA1': 43, 'BARRA2': 25},
                                 {'HORA': 2, 'BARRA0': 0, 'BARRA1': 25, 'BARRA2': 25}]

        self.bus = BBus(DBAR=self.DGER, DLIN=self.DLIN, ref_bar=self.ref_bus_i)

        self.tempDim = [value for value in range(len(self.DEMANDA))]
        self.gerDim = [value for value in range(len(self.DGER))]

        self.fromLin = [value['DE'] for value in self.DLIN]
        self.toLin = [value['PARA'] for value in self.DLIN]

    def minimize(self):
        raise NotImplementedError

    def obj_fn(self, model):
        raise NotImplementedError


    def ger_limits(self, model, i, t):
        return (self.DGER[i]['MIN'], self.DGER[i]['MAX'])

    def deficit_lim(self, model, i, t):
        return (0, self.DEMANDA[t][f'BARRA{i}'])



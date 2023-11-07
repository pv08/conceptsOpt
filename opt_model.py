
from pyomo.environ import *
from src.network import Network
class OptModel(Network):
    def __init__(self):
        super(OptModel, self).__init__()
        self.network = Network()
        self.model = ConcreteModel("CustoTermo")


        self.model.pger = Var(self.gerDim, self.tempDim, within=NonNegativeReals, bounds=self.ger_limits)
        self.model.p_emissao = Var(self.gerDim, self.tempDim, within=NonNegativeReals, bounds=self.ger_limits)
        self.model.deficit = Var(self.tempDim, within=NonNegativeReals)

        self.model.thetas = Var(self.gerDim, self.tempDim, within=Reals)


        #Problem Constraints
        self.model.consGerEqualDem = ConstraintList()
        self.model.consTheta = ConstraintList()
        self.model.consPowerflow = ConstraintList()
        self.model.consTechnicalLineLimit = ConstraintList()
        self.model.consUpRamp = ConstraintList()
        self.model.consDownRamp = ConstraintList()

        self.solver = SolverFactory('glpk')




    def saveModelDetails(self, name):
        with open(f'results/{name}.txt', 'w') as file:
            self.model.pprint(ostream=file)


    def minimize(self):
        self.model.obj = Objective(rule=self.obj_fn, sense=minimize)
        self.createConstraints()
        self.results = self.solver.solve(self.model, tee=False)
        print(self.results)
        self.saveModelDetails(name='CustoTermo')


    def obj_fn(self, model):
        C_P = 0
        E_P = 0
        for t, temp in enumerate(self.tempDim):
            for k, ger in enumerate(self.gerDim):
                #Custo
                C_P += (self.network.DGER[k]['CUSTO'] * self.model.pger[ger, temp])
                #Emissão
                E_P += (self.network.DGER[k]['CO2'] * self.model.p_emissao[ger, temp])
            C_P += self.cost_deficit * self.model.deficit[temp]
            E_P += (1 - self.delta) * self.h * E_P
        fob = (self.delta * C_P) + ((1 - self.delta) * self.h * E_P)
        return fob

    def createConstraints(self):
        #Colocar os angulos da barra slack como o
        for temp in self.tempDim:
            self.model.consTheta.add(self.model.thetas[self.ref_bus_i, temp] == 0.0)

        #Constraint que define que a geração tem que ser igual a demanda
        for t, temp in enumerate(self.tempDim):
            _ger = 0
            for k, ger in enumerate(self.gerDim):
                _ger += self.model.pger[ger, t]
            _ger += self.model.deficit[t]
            self.model.consGerEqualDem.add(_ger == sum(self.network.DEMANDA[temp][f"BARRA{ger}"] for ger in self.gerDim))

        #Define os limites das linhas e o fluxo de potência
        for time, load_data in zip(self.tempDim, self.DEMANDA):
            for line in self.DLIN:
                from_bus = line['DE']
                to_bus = line['PARA']
                diff_thetas = (self.model.thetas[from_bus, time] - self.model.thetas[to_bus, time])
                thermal_sum = 0
                for thermal_idx, thermal_data in enumerate(self.DGER):
                    if thermal_data['BARRA'] == from_bus:
                        thermal_sum += self.model.pger[thermal_idx, time]

                load_sum = load_data[f'BARRA{from_bus}']


                loss = .5*(line['CONDUTANCIA']/2) * diff_thetas
                power_flow = diff_thetas / (1/line['SUSCEPTANCIA'])

                self.model.consPowerflow.add(thermal_sum - load_sum - loss == power_flow)
                self.model.consTechnicalLineLimit.add(power_flow <= line['LIMITES'])

        # Define as constraints de rampa
        for temp in range(1, len(self.tempDim)):
            for ger in self.gerDim:
                self.model.consUpRamp.add(self.model.pger[ger, temp] - self.model.pger[ger, temp - 1] <= self.DGER[ger]['RAMPA'])
                self.model.consDownRamp.add(self.model.pger[ger, temp - 1] - self.model.pger[ger, temp] <= self.DGER[ger]['RAMPA'])

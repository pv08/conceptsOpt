from pyomo.core.kernel.constraint import constraint
from pyomo.environ import *
from src.network import Network
class OptModel(Network):
    def __init__(self):
        super(OptModel, self).__init__()
        self.network = Network()
        self.model = ConcreteModel("CustoTermo")

        self.model.EqCons = ConstraintList()
        self.model.FlowLineLimitCons = ConstraintList()
        self.model.KirchhoffCons = ConstraintList()
        self.solver = SolverFactory('glpk')



        self.model.pger = Var(self.gerDim, self.tempDim, within=NonNegativeReals, bounds=self.ger_limits)
        self.model.thetas = Var(self.gerDim, self.tempDim, within=NonNegativeReals)
        self.model.powerflow = Var(list(zip(self.fromLin, self.toLin)), self.tempDim)

        self.model.p_emissao = Var(self.gerDim, self.tempDim, within=NonNegativeReals, bounds=self.ger_limits)
        self.model.deficit = Var(self.tempDim, within=NonNegativeReals)



    def minimize(self):
        self.model.obj = Objective(rule=self.obj_fn, sense=minimize)
        self.createConstraints()


    def obj_fn(self, model):
        C_P = 0
        E_P = 0
        for t, temp in enumerate(self.tempDim):
            for k, ger in enumerate(self.gerDim):
                #Custo
                C_P += (self.network.DGER[k]['CUSTO'] * self.model.pger[ger, temp])
                #Emissão
                E_P += (self.network.DGER[k]['CO2'] * self.model.p_emissao[ger, temp])
            # TODO: {deficits fica aqui?}
            C_P += self.cost_deficit * self.model.deficit[temp]
            E_P += (1 - self.delta) * self.h * E_P
        fob = (self.delta * C_P) + ((1 - self.delta) * self.h * E_P)
        return fob

    def createConstraints(self):
        #Constraint que define que a geração tem que ser igual a demanda
        for t, temp in enumerate(self.tempDim):
            _ger = 0
            for k, ger in enumerate(self.gerDim):
                _ger += self.model.pger[ger, t]
            _ger += self.model.deficit[t]
            self.model.EqCons.add(_ger == sum(self.network.DEMANDA[temp][f"BARRA{ger}"] for ger in self.gerDim))

        #Define os thetas e acha os fluexo de potência
        for temp in self.tempDim:
            thetas = []
            powerflow = []
            losses = []
            for l, line in enumerate(self.bus.b_):
                summation = 0
                for cel, ger in zip(line, self.gerDim):
                    summation += cel * (self.model.pger[ger, temp] - self.DEMANDA[temp][f'BARRA{ger}'])
                thetas.append(summation)
            #Define o fluxo de potência de cada linha
            for conn_i, conn in enumerate(zip(self.fromLin, self.toLin)):
                k, m = conn
                flow = (thetas[k] - thetas[m]) / (1/self.DLIN[conn_i]['SUSCEPTANCIA'])
                powerflow.append({conn: flow})
                self.model.FlowLineLimitCons.add(flow <= self.DLIN[conn_i]['LIMITES'])


            #Constraint para os fluxos que entram e todos o que saem
            for bar in self.gerDim:
                sum_flow = 0
                #Todos mundo que ta indo
                for flow in powerflow:
                    conn = flow.keys()
                    k, _ = conn
                    if bar == k:
                        sum_flow += flow
                #Todos mundo que ta voltando
                for flow in powerflow:
                    conn = flow.keys()
                    _, m = conn
                    if bar == m:
                        sum_flow -= powerflow
                        #TODO {Colocar a demanda}
                self.model.KirchhoffCons.add(sum_flow == self.model.pger[bar, temp] - self.DEMANDA[temp][f'BARRA{bar}'])
